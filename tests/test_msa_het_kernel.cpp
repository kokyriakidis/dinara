// Unit tests for the pure-matrix MSA het kernels (src/MsaHetKernel.hpp).
//
// These exercise the correctness-critical coordinate logic of the abPOA-MSA
// het-detection / 2-base-allele-anchor plan with synthetic MSA matrices, with
// no abPOA or Assembler dependency. The matrices below use:
//   base encoding 0=A 1=C 2=G 3=T, gap = 5 (matching abPOA abpt->m default).

#include "../external/catch2/catch.hpp"
#include "../src/MsaHetKernel.hpp"

#include <vector>

using namespace dinara::msahet;
using std::vector;

namespace {

constexpr int GAP = 5;

// Adapter: vector<vector<uint8_t>> -> const uint8_t* const* for the kernels.
struct Mat {
    vector<vector<uint8_t>> rows;
    vector<const uint8_t*> ptrs;
    const uint8_t* const* data() {
        ptrs.clear();
        for (auto& r : rows) ptrs.push_back(r.data());
        return ptrs.data();
    }
    int nSeq() const { return int(rows.size()); }
    int msaLen() const { return rows.empty() ? 0 : int(rows[0].size()); }
};

} // namespace

// ---------------------------------------------------------------------------
// Pass 1: column <-> backbone maps.
// ---------------------------------------------------------------------------
TEST_CASE("buildColumnMaps: no insertions, identity mapping", "[msahet]") {
    //      col: 0 1 2 3
    // backbone: A C G T
    vector<uint8_t> row0 = {0, 1, 2, 3};
    const ColumnMaps m = buildColumnMaps(row0.data(), int(row0.size()), GAP);

    REQUIRE(m.backboneLen == 4);
    for (int c = 0; c < 4; ++c) {
        REQUIRE(m.colToBb[c] == c);
        REQUIRE(m.bbToCol[c] == c);
    }
}

TEST_CASE("buildColumnMaps: member-insertion gap columns in row 0", "[msahet]") {
    //      col: 0 1 2 3 4
    // backbone: A C - G T   (gap at col 2 = a member insertion)
    vector<uint8_t> row0 = {0, 1, GAP, 2, 3};
    const ColumnMaps m = buildColumnMaps(row0.data(), int(row0.size()), GAP);

    REQUIRE(m.backboneLen == 4);
    REQUIRE(m.colToBb[0] == 0);
    REQUIRE(m.colToBb[1] == 1);
    REQUIRE(m.colToBb[2] == -1);   // insertion column: not a backbone base
    REQUIRE(m.colToBb[3] == 2);
    REQUIRE(m.colToBb[4] == 3);

    // bb -> col must skip the inserted column.
    REQUIRE(m.bbToCol[0] == 0);
    REQUIRE(m.bbToCol[1] == 1);
    REQUIRE(m.bbToCol[2] == 3);    // backbone base 2 lives at column 3
    REQUIRE(m.bbToCol[3] == 4);
}

// ---------------------------------------------------------------------------
// Pass 2: column tally + alt base.
// ---------------------------------------------------------------------------
TEST_CASE("tallyColumn: het column with clear alt and strand split", "[msahet]") {
    // backbone base = A (0) at col 0. Members: 3 carry A (ref), 3 carry G (alt).
    Mat mat;
    mat.rows = {
        {0}, // row 0 backbone = A
        {0}, {0}, {0}, // ref members
        {2}, {2}, {2}, // alt members = G
    };
    // strands for the 6 members: first three fwd, last three: 2 fwd 1 rev.
    vector<bool> strandIsRev = {false, false, false, false, true, false};

    const ColumnTally t = tallyColumn(mat.data(), mat.nSeq(), 0, GAP, strandIsRev);
    REQUIRE(t.count[0] == 3); // A
    REQUIRE(t.count[2] == 3); // G
    REQUIRE(t.covered == 6);
    REQUIRE(t.altBase(/*refBase=*/0) == 2); // alt = G

    // strand split for the alt allele (G): 2 fwd, 1 rev.
    REQUIRE(t.fwd[2] == 2);
    REQUIRE(t.rev[2] == 1);
}

TEST_CASE("tallyColumn: gaps and N are non-informative", "[msahet]") {
    Mat mat;
    mat.rows = {
        {0},        // backbone
        {0}, {0},   // ref
        {2},        // alt
        {GAP},      // deletion: gap
        {4},        // N
    };
    vector<bool> strandIsRev = {false, false, false, false, false};
    const ColumnTally t = tallyColumn(mat.data(), mat.nSeq(), 0, GAP, strandIsRev);
    REQUIRE(t.count[0] == 2);
    REQUIRE(t.count[2] == 1);
    REQUIRE(t.covered == 3); // gap + N excluded
}

// ---------------------------------------------------------------------------
// Pass 3: read-position recovery along an MSA row.
// ---------------------------------------------------------------------------
TEST_CASE("recoverReadPositions: non-gaps advance, gaps do not", "[msahet]") {
    // Row: A C - G T   (gap at col 2). rowStartReadPos = 100.
    // read positions:  100 101 _ 102 103  at cols 0,1,3,4.
    vector<uint8_t> row = {0, 1, GAP, 2, 3};
    vector<int> cols = {0, 1, 2, 3, 4}; // ask for all columns, sorted
    std::unordered_map<uint64_t, uint32_t> out;
    recoverReadPositions(row.data(), int(row.size()), GAP,
                         /*rowStartReadPos=*/100, /*row=*/7, cols, out);

    auto key = [](int r, int c) { return (uint64_t(uint32_t(r)) << 32) | uint32_t(c); };
    REQUIRE(out[key(7, 0)] == 100);
    REQUIRE(out[key(7, 1)] == 101);
    REQUIRE(out.find(key(7, 2)) == out.end()); // gap column: no position recorded
    REQUIRE(out[key(7, 3)] == 102);
    REQUIRE(out[key(7, 4)] == 103);
}

TEST_CASE("recoverReadPositions: leading gaps before first base", "[msahet]") {
    // Row: - - A C   start read pos 50 -> A at col2 = 50, C at col3 = 51.
    vector<uint8_t> row = {GAP, GAP, 0, 1};
    vector<int> cols = {0, 1, 2, 3};
    std::unordered_map<uint64_t, uint32_t> out;
    recoverReadPositions(row.data(), int(row.size()), GAP, 50, 3, cols, out);

    auto key = [](int r, int c) { return (uint64_t(uint32_t(r)) << 32) | uint32_t(c); };
    REQUIRE(out.find(key(3, 0)) == out.end());
    REQUIRE(out.find(key(3, 1)) == out.end());
    REQUIRE(out[key(3, 2)] == 50);
    REQUIRE(out[key(3, 3)] == 51);
}

// ---------------------------------------------------------------------------
// Agreed-neighbor selection for the 2-base allele anchor.
// ---------------------------------------------------------------------------
TEST_CASE("chooseNeighbor: right neighbor unanimous, picks (bb,bb+1)", "[msahet]") {
    // backbone: A C G T  (cols 0..3). Het at bb=1 (backbone C). Alt allele = A.
    // Alt-group members all read A at col1 and agree on G at col2 (right nbr).
    Mat mat;
    mat.rows = {
        {0, 1, 2, 3}, // backbone
        {0, 0, 2, 3}, // alt member: col1=A, col2=G
        {0, 0, 2, 3},
        {0, 0, 2, 3},
    };
    vector<uint8_t> row0 = {0, 1, 2, 3};
    const ColumnMaps maps = buildColumnMaps(row0.data(), 4, GAP);
    vector<int> group = {1, 2, 3};

    const NeighborChoice nc = chooseNeighbor(mat.data(), GAP, maps,
                                             /*bb=*/1, /*alleleBase=*/0, group);
    REQUIRE(nc.valid);
    REQUIRE(nc.leftBb == 1);       // (bb, bb+1) => left base is the het base
    REQUIRE(nc.leftCol == 1);
    REQUIRE(nc.base0 == 0);        // A (allele) at bb
    REQUIRE(nc.base1 == 2);        // G (neighbor consensus) at bb+1
    REQUIRE(nc.matchCount == 3);
}

TEST_CASE("chooseNeighbor: right neighbor noisy, falls back to left", "[msahet]") {
    // Het at bb=2 (col2). Alt allele = A at col2.
    // Right neighbor (col3 = bb+1) disagrees across the group (A/C/G);
    // left neighbor (col1 = bb-1) is unanimous G -> left pair wins on matches.
    Mat mat;
    mat.rows = {
        {3, 2, 1, 0, 3}, // backbone: T G C A T (cols 0..4)
        {3, 2, 0, 0, 3}, // alt: col2=A, left col1=G (agree), right col3=A
        {3, 2, 0, 1, 3}, // alt: col2=A, left col1=G (agree), right col3=C
        {3, 2, 0, 2, 3}, // alt: col2=A, left col1=G (agree), right col3=G
    };
    vector<uint8_t> row0 = {3, 2, 1, 0, 3};
    const ColumnMaps maps = buildColumnMaps(row0.data(), 5, GAP);
    vector<int> group = {1, 2, 3};

    const NeighborChoice nc = chooseNeighbor(mat.data(), GAP, maps,
                                             /*bb=*/2, /*alleleBase=*/0, group);
    REQUIRE(nc.valid);
    REQUIRE(nc.leftBb == 1);       // (bb-1, bb) chosen
    REQUIRE(nc.leftCol == 1);
    REQUIRE(nc.base0 == 2);        // G (left neighbor consensus) at bb-1
    REQUIRE(nc.base1 == 0);        // A (allele) at bb
    REQUIRE(nc.matchCount == 3);   // all 3 agree on (G,A)
    // Right pair would have only 1 match (col3 majority is A, only row1 has it).
}

TEST_CASE("chooseNeighbor: dissenting read dropped from match count", "[msahet]") {
    // Het at bb=0 (no left neighbor, so the right pair is forced). Alt=A at col0.
    // Right neighbor col1 majority = G, but one read has T and one has a gap ->
    // match count includes only the reads agreeing on the (A,G) 2-mer.
    Mat mat;
    mat.rows = {
        {1, 2},   // backbone: C G
        {0, 2},   // alt A, neighbor G  (match)
        {0, 2},   // alt A, neighbor G  (match)
        {0, 3},   // alt A, neighbor T  (dissent)
        {0, GAP}, // alt A, neighbor gap (dropped)
    };
    vector<uint8_t> row0 = {1, 2};
    const ColumnMaps maps = buildColumnMaps(row0.data(), 2, GAP);
    vector<int> group = {1, 2, 3, 4};

    const NeighborChoice nc = chooseNeighbor(mat.data(), GAP, maps,
                                             /*bb=*/0, /*alleleBase=*/0, group);
    REQUIRE(nc.valid);
    REQUIRE(nc.leftBb == 0);
    REQUIRE(nc.base0 == 0);        // A (allele) at bb
    REQUIRE(nc.base1 == 2);        // majority neighbor = G
    REQUIRE(nc.matchCount == 2);   // only the two (A,G) reads agree
}

TEST_CASE("chooseNeighbor: het at backbone end has only one neighbor", "[msahet]") {
    // Het at bb=0 (no left neighbor). Must use right pair (0,1).
    Mat mat;
    mat.rows = {
        {1, 2}, // backbone C G
        {0, 2}, // alt A, right neighbor G
        {0, 2},
    };
    vector<uint8_t> row0 = {1, 2};
    const ColumnMaps maps = buildColumnMaps(row0.data(), 2, GAP);
    vector<int> group = {1, 2};

    const NeighborChoice nc = chooseNeighbor(mat.data(), GAP, maps,
                                             /*bb=*/0, /*alleleBase=*/0, group);
    REQUIRE(nc.valid);
    REQUIRE(nc.leftBb == 0);
    REQUIRE(nc.base0 == 0); // A
    REQUIRE(nc.base1 == 2); // G
    REQUIRE(nc.matchCount == 2);
}

// ---------------------------------------------------------------------------
// End-to-end coordinate round-trip: the readPos recovered for an allele anchor
// must point at the allele base in the read.
// ---------------------------------------------------------------------------
TEST_CASE("anchor coordinate round-trip across an insertion column", "[msahet]") {
    // backbone: A C G T  but member row has an insertion between C and G, so the
    // MSA has 5 columns and a gap appears in row 0 at col 2.
    //   col:        0 1 2 3 4
    //   row0 (bb):  A C - G T
    //   member:     A A T G T   (insertion 'T' at col2; het A at bb=1/col1)
    Mat mat;
    mat.rows = {
        {0, 1, GAP, 2, 3},
        {0, 0, 3,   2, 3},
    };
    const ColumnMaps maps = buildColumnMaps(mat.rows[0].data(), 5, GAP);
    REQUIRE(maps.backboneLen == 4);
    REQUIRE(maps.bbToCol[1] == 1); // het backbone base 1 at col 1
    REQUIRE(maps.bbToCol[2] == 3); // backbone base 2 at col 3 (after insertion)

    // Choose neighbor for the member's alt allele A at bb=1.
    vector<int> group = {1};
    const NeighborChoice nc = chooseNeighbor(mat.data(), GAP, maps, 1, 0, group);
    REQUIRE(nc.valid);
    REQUIRE(nc.leftBb == 1);
    REQUIRE(nc.leftCol == 1);

    // Recover the member's read position at the anchor's left column.
    // Member row non-gaps: cols 0,1,2,3,4 are all non-gap -> read pos = col+start.
    std::unordered_map<uint64_t, uint32_t> out;
    vector<int> cols = {nc.leftCol};
    recoverReadPositions(mat.rows[1].data(), 5, GAP, /*start=*/200, /*row=*/1, cols, out);
    const uint64_t key = (uint64_t(1u) << 32) | uint32_t(nc.leftCol);
    REQUIRE(out.count(key) == 1);
    const uint32_t readPos = out[key];
    REQUIRE(readPos == 201); // start 200 at col0, col1 -> 201

    // Round-trip: the member's base at that read position must equal base0 (A).
    // In this synthetic row the read sequence (non-gap) is A A T G T, so index 1
    // (=201-200) is 'A' (0) == nc.base0.
    vector<uint8_t> readSeq;
    for (uint8_t v : mat.rows[1]) if (int(v) != GAP) readSeq.push_back(v);
    REQUIRE(readSeq[readPos - 200] == nc.base0);
}
