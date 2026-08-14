/**
 * @file test_phasing_coordinates.cpp
 * @brief Tests that CIGAR coordinates and base lookups are correct for
 *        both same-strand and reverse-complement overlaps.
 *
 * Verifies the coordinate contract between:
 *   - AlignmentData (ad.qs/qe/ts/te: marker-based forward-strand coordinates)
 *   - OverlapCigarStore (yk in oriented/RC coordinates for strand-1 reads)
 *   - RC conversion: oriented_start = targetLen - ad.te
 *   - getBaseAtPosition() base extraction
 */

#include "../external/catch2/catch.hpp"
#include "../src/Assembler.hpp"
#include "../src/OverlapCigarStore.hpp"
#include "../src/Alignment.hpp"
#include "../src/Base.hpp"
#include "../src/Reads.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using namespace dinara;

// ============================================================================
// Helpers
// ============================================================================

static std::string reverseComplement(const std::string& seq) {
    std::string rc = seq;
    std::reverse(rc.begin(), rc.end());
    for (auto& b : rc) {
        if (b == 'A') b = 'T';
        else if (b == 'C') b = 'G';
        else if (b == 'G') b = 'C';
        else if (b == 'T') b = 'A';
    }
    return rc;
}

static std::string randomSequence(size_t length, uint32_t seed = 42) {
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dis(0, 3);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string seq(length, 'N');
    for (size_t i = 0; i < length; i++) seq[i] = bases[dis(gen)];
    return seq;
}

/// Map ASCII base to Base::value (A=0, C=1, G=2, T=3).
static uint8_t charToBaseValue(char c) {
    switch (c) {
        case 'A': return 0; case 'C': return 1;
        case 'G': return 2; case 'T': return 3;
    }
    return 255;
}

/// Complement of a base value.
static uint8_t complementValue(uint8_t v) { return 3 - v; }

/// Silence cout/cerr during a block.
class ScopedSilence {
    std::ofstream null{"/dev/null"};
    std::streambuf *ob, *eb;
public:
    ScopedSilence() : ob(std::cout.rdbuf(null.rdbuf())),
                      eb(std::cerr.rdbuf(null.rdbuf())) {}
    ~ScopedSilence() { std::cout.rdbuf(ob); std::cerr.rdbuf(eb); }
};

/// Create a temp directory that is cleaned up on destruction.
struct TempDir {
    fs::path path;
    TempDir(const std::string& prefix) {
        char tmpl[256];
        snprintf(tmpl, sizeof(tmpl), "/tmp/%s_XXXXXX", prefix.c_str());
        path = mkdtemp(tmpl);
    }
    ~TempDir() { fs::remove_all(path); }
};

/// Write sequences as FASTQ.
static void writeFastq(const fs::path& p, const std::vector<std::string>& seqs) {
    std::ofstream out(p);
    for (size_t i = 0; i < seqs.size(); i++) {
        out << "@read_" << i << "\n"
            << seqs[i] << "\n+\n"
            << std::string(seqs[i].size(), '~') << "\n";
    }
}
// ============================================================================
// Unit test: manually constructed CIGAR, verify walkRange coordinates
// ============================================================================

TEST_CASE("CIGAR walk: same-strand overlap coordinates are forward",
          "[phasing][coordinates][same-strand]")
{
    // Simulate: read0 (len=1000) overlaps read1 (len=1000) on same strand.
    // Overlap region: read0 [200, 500), read1 [300, 600).
    // CIGAR: 280 match, 5 mismatch, 15 match = 300 bases on each read.
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(280);
    store.pushMismatch(5);
    store.pushMatch(15);
    uint32_t cnt = store.tokensSince(off);

    // For same-strand: ad.ts == oriented start == 300.
    const uint64_t read0Start = 200;
    const uint64_t read1Start = 300; // forward == oriented for same-strand

    // Walk full range and collect all (op, len, xk, yk).
    struct Event { uint8_t op; uint32_t len; uint64_t xk; uint64_t yk; };
    std::vector<Event> events;
    store.walkRange(off, cnt, read0Start, read1Start, 0, UINT32_MAX,
        [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
            events.push_back({op, len, xk, yk});
        });

    REQUIRE(events.size() == 3);

    // Match: starts at (200, 300), len=280
    CHECK(events[0].op == 0);
    CHECK(events[0].len == 280);
    CHECK(events[0].xk == 200);
    CHECK(events[0].yk == 300);

    // Mismatch: starts at (480, 580), len=5
    CHECK(events[1].op == 1);
    CHECK(events[1].len == 5);
    CHECK(events[1].xk == 480);
    CHECK(events[1].yk == 580);

    // Match: starts at (485, 585), len=15
    CHECK(events[2].op == 0);
    CHECK(events[2].len == 15);
    CHECK(events[2].xk == 485);
    CHECK(events[2].yk == 585);

    // Final positions: xk=500, yk=600 — matches overlap end.
    uint64_t finalXk = events.back().xk + events.back().len;
    uint64_t finalYk = events.back().yk + events.back().len;
    CHECK(finalXk == 500);
    CHECK(finalYk == 600);
}

TEST_CASE("CIGAR walk: RC overlap needs oriented read1Start",
          "[phasing][coordinates][rc]")
{
    // Simulate: read0 (len=1000) overlaps read1 (len=800) on opposite strands.
    // Forward coords: read0 [100, 400), read1 [200, 500) (forward).
    // RC coords of read1: rcStart = 800 - 500 = 300, rcEnd = 800 - 200 = 600.
    // The CIGAR was built with read1 in RC orientation, so yk starts at 300.
    //
    // ad.ts = 200 (forward), ad.te = 500 (forward).
    // cigarRead1Start should return 800 - 500 = 300 (RC start).

    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(200);
    store.pushMismatch(10);
    store.pushMatch(90);
    uint32_t cnt = store.tokensSince(off);

    const uint64_t read0Start = 100;
    const uint32_t targetLen = 800;
    const uint32_t adTs = 200; // forward start
    const uint32_t adTe = 500; // forward end

    // cigarRead1Start conversion: targetLen - adTe = 800 - 500 = 300
    const uint64_t read1StartOriented = targetLen - adTe;
    CHECK(read1StartOriented == 300);

    // Walk with oriented start.
    struct Event { uint8_t op; uint32_t len; uint64_t xk; uint64_t yk; };
    std::vector<Event> events;
    store.walkRange(off, cnt, read0Start, read1StartOriented, 0, UINT32_MAX,
        [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
            events.push_back({op, len, xk, yk});
        });

    REQUIRE(events.size() == 3);

    // yk starts at 300 (RC coordinate), advances to 600.
    CHECK(events[0].yk == 300);
    CHECK(events[1].yk == 500); // 300 + 200
    CHECK(events[2].yk == 510); // 300 + 200 + 10

    uint64_t finalYk = events.back().yk + events.back().len;
    CHECK(finalYk == 600); // RC end = targetLen - adTs = 800 - 200 = 600

    // WRONG: if we used ad.ts (200) as read1Start, yk would start at 200.
    std::vector<Event> wrongEvents;
    store.walkRange(off, cnt, read0Start, adTs, 0, UINT32_MAX,
        [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
            wrongEvents.push_back({op, len, xk, yk});
        });
    // This would give yk starting at 200 — wrong for RC overlaps.
    CHECK(wrongEvents[0].yk == 200); // demonstrates the bug if we don't convert
    CHECK(wrongEvents[0].yk != events[0].yk); // different from correct value
}

TEST_CASE("CIGAR walk: base lookup with oriented yk for RC overlap",
          "[phasing][coordinates][base-lookup]")
{
    // Create two known sequences.
    // read0 (forward): ACGTACGTAC (len=10)
    // read1 (forward): TTTTACGTAA (len=10), RC = TTACGTAAAA
    //
    // Overlap: read0[2..6) = "GTAC" aligns to read1-RC[0..4) = "TTAC"
    //   Position 0: G vs T = mismatch
    //   Position 1: T vs T = match
    //   Position 2: A vs A = match
    //   Position 3: C vs C = match
    //
    // CIGAR: mismatch(1), match(3)
    // read0Start = 2, read1Start (RC) = 0
    //
    // ad.ts (forward) = 10 - 4 = 6, ad.te (forward) = 10 - 0 = 10
    // cigarRead1Start = 10 - 10 = 0 ✓

    const std::string read0_fwd = "ACGTACGTAC";
    const std::string read1_fwd = "TTTTACGTAA";
    const std::string read1_rc  = reverseComplement(read1_fwd);
    // read1_rc = "TTACGTAAAA"

    REQUIRE(read1_rc == "TTACGTAAAA");

    // Verify the overlap: read0[2..6) vs read1_rc[0..4)
    CHECK(read0_fwd[2] == 'G'); CHECK(read1_rc[0] == 'T'); // mismatch
    CHECK(read0_fwd[3] == 'T'); CHECK(read1_rc[1] == 'T'); // match
    CHECK(read0_fwd[4] == 'A'); CHECK(read1_rc[2] == 'A'); // match
    CHECK(read0_fwd[5] == 'C'); CHECK(read1_rc[3] == 'C'); // match

    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMismatch(1);
    store.pushMatch(3);
    uint32_t cnt = store.tokensSince(off);

    const uint64_t read0Start = 2;
    const uint64_t read1StartOriented = 0; // RC coordinate

    // Walk and verify each base.
    struct Event { uint8_t op; uint32_t len; uint64_t xk; uint64_t yk; };
    std::vector<Event> events;
    store.walkRange(off, cnt, read0Start, read1StartOriented, 0, UINT32_MAX,
        [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
            events.push_back({op, len, xk, yk});
        });

    REQUIRE(events.size() == 2);

    // For each position, verify base lookup.
    // read0 base: sequence[xk] (always forward)
    // read1 base (RC): sequence[baseCount - 1 - yk].complement()
    const uint32_t read1Len = 10;

    for (const auto& ev : events) {
        for (uint32_t b = 0; b < ev.len; b++) {
            uint32_t xpos = uint32_t(ev.xk) + b;
            uint32_t ypos = uint32_t(ev.yk) + b;

            uint8_t read0Base = charToBaseValue(read0_fwd[xpos]);
            // RC lookup: sequence[baseCount-1-ypos].complement()
            uint8_t read1Base = complementValue(
                charToBaseValue(read1_fwd[read1Len - 1 - ypos]));

            if (ev.op == 0) { // match
                CHECK(read0Base == read1Base);
            } else { // mismatch
                CHECK(read0Base != read1Base);
            }
        }
    }
}

TEST_CASE("CIGAR walk: indels advance only one coordinate",
          "[phasing][coordinates][indels]")
{
    // CIGAR: match(10), insertion(3), match(10), deletion(5), match(10)
    // SAM/PAF convention:
    //   Insertion: xk (query) advances by 3, yk stays.
    //   Deletion:  yk (target) advances by 5, xk stays.
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(10);
    store.pushInsertion(3);
    store.pushMatch(10);
    store.pushDeletion(5);
    store.pushMatch(10);
    uint32_t cnt = store.tokensSince(off);

    struct Event { uint8_t op; uint32_t len; uint64_t xk; uint64_t yk; };
    std::vector<Event> events;
    store.walkRange(off, cnt, 100, 200, 0, UINT32_MAX,
        [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
            events.push_back({op, len, xk, yk});
        });

    REQUIRE(events.size() == 5);

    // Event positions are the op's start coordinates.
    // match(10) at xk=100,yk=200 → ends xk=110, yk=210
    CHECK(events[1].xk == 110); CHECK(events[1].yk == 210);
    // insertion(3) at xk=110,yk=210 → xk advances to 113, yk unchanged
    CHECK(events[2].xk == 113); CHECK(events[2].yk == 210);
    // match(10) at xk=113,yk=210 → ends xk=123, yk=220
    CHECK(events[3].xk == 123); CHECK(events[3].yk == 220);
    // deletion(5) at xk=123,yk=220 → yk advances to 225, xk unchanged
    CHECK(events[4].xk == 123); CHECK(events[4].yk == 225);
    // Final: match(10) at xk=123,yk=225 → xk=133, yk=235
    uint64_t finalXk = events[4].xk + events[4].len;
    uint64_t finalYk = events[4].yk + events[4].len;
    CHECK(finalXk == 133);
    CHECK(finalYk == 235);
}
// ============================================================================
// Integration test: full pipeline, verify CIGAR bases match read sequences
// ============================================================================

TEST_CASE("Integration: CIGAR walk bases match read sequences for all overlaps",
          "[phasing][coordinates][integration]")
{
    // Create reads with known overlapping regions.
    // Use a long enough shared region to produce markers and candidates.
    const std::string shared = randomSequence(2000, 55555);

    // read0: prefix + shared region
    const std::string prefix0 = randomSequence(300, 10001);
    const std::string suffix0 = randomSequence(300, 10002);
    const std::string r0 = prefix0 + shared + suffix0;

    // read1: same-strand overlap with read0 (shares the same region)
    const std::string prefix1 = randomSequence(400, 10003);
    const std::string suffix1 = randomSequence(400, 10004);
    const std::string r1 = prefix1 + shared + suffix1;

    // read2: RC overlap with read0 (RC of the shared region)
    const std::string prefix2 = randomSequence(350, 10005);
    const std::string suffix2 = randomSequence(350, 10006);
    const std::string r2 = prefix2 + reverseComplement(shared) + suffix2;

    // read3: introduce 2 SNPs in the shared region for mismatch testing
    std::string sharedMut = shared;
    sharedMut[500] = (sharedMut[500] == 'A') ? 'C' : 'A';
    sharedMut[1500] = (sharedMut[1500] == 'G') ? 'T' : 'G';
    const std::string prefix3 = randomSequence(280, 10007);
    const std::string suffix3 = randomSequence(280, 10008);
    const std::string r3 = prefix3 + sharedMut + suffix3;

    const std::vector<std::string> seqs = {r0, r1, r2, r3};

    TempDir tmp("dinara_phasing_coord_test");
    const auto fastqPath = tmp.path / "reads.fastq";
    writeFastq(fastqPath, seqs);

    // Run the pipeline.
    std::unique_ptr<Assembler> assembler;
    {
        ScopedSilence silence;

        auto prevDir = fs::current_path();
        fs::current_path(tmp.path);

        assembler = std::make_unique<Assembler>(
            tmp.path.string() + "/", true, 0, 4096);
        assembler->addReads(fastqPath.string(), 0, true, 1);
        assembler->computeReadIdsSortedByName();
        assembler->findMarkersSimdMinimizers(1, 16, 4, /*useHifiasm*/ false);
        assembler->countKmersFromMarkerKmerIds(1);
        assembler->applyKmerCountFilter(1, 1000, 1);

        OverlapCandidatesOptions candOpts;
        candOpts.method = "InvertedIndex";
        candOpts.driftRateTolerance = 0.1;
        candOpts.minChainMarkerCount = 4;
        assembler->findAlignmentCandidatesInvertedIndex(
            0.1, 100, candOpts, 1);

        AlignOptions alignOpts;
        alignOpts.alignMethod = 6;
        alignOpts.maxSkip = 100;
        alignOpts.maxDrift = 100;
        alignOpts.maxTrim = 10000;
        alignOpts.minAlignedMarkerCount = 2;
        alignOpts.minAlignedFraction = 0.0;
        alignOpts.maxMarkerFrequency = 1000;
        alignOpts.matchScore = 3;
        alignOpts.mismatchScore = -1;
        alignOpts.gapScore = -1;
        alignOpts.downsamplingFactor = 0.1;
        alignOpts.bandExtend = 10;
        alignOpts.maxBand = 1000;
        alignOpts.sameChannelReadAlignmentSuppressDeltaThreshold = 0;
        alignOpts.suppressContainments = false;
        alignOpts.align4DeltaX = 200;
        alignOpts.align4DeltaY = 10;
        alignOpts.align4MinEntryCountPerCell = 10;
        alignOpts.align4MaxDistanceFromBoundary = 100;
        alignOpts.align5DriftRateTolerance = 0.02;
        alignOpts.align5MinBandExtend = 10;
        alignOpts.maxErrorRate = 0.3;
        alignOpts.overlapDpMatchScore = 2;
        alignOpts.overlapDpMismatchScore = -4;
        alignOpts.overlapDpGapOpen1 = 4;
        alignOpts.overlapDpGapExtend1 = 2;
        alignOpts.overlapDpGapOpen2 = 24;
        alignOpts.overlapDpGapExtend2 = 1;

        assembler->computeBaseAlignmentsAndStore(alignOpts, 1);

        fs::current_path(prevDir);
    }

    // Now verify: for every alignment with a CIGAR, walk the CIGAR and
    // check that match positions have identical bases and mismatch positions
    // have different bases.
    const auto& cigarStore = assembler->getOverlapCigarStore();
    const auto& reads = assembler->getReads();
    size_t alignmentsChecked = 0;
    size_t basesVerified = 0;
    size_t rcAlignmentsChecked = 0;

    for (size_t ai = 0; ai < assembler->alignmentData.size(); ai++) {
        const auto& ad = assembler->alignmentData[ai];
        const auto& info = ad.info;

        if (info.cigarOffset == uint32_t(-1)) continue;
        if (info.cigarTokenCount == 0) continue;

        const ReadId read0Id = ad.readIds[0];
        const ReadId read1Id = ad.readIds[1];
        const auto read0Seq = reads.getRead(read0Id);
        const auto read1Seq = reads.getRead(read1Id);
        const uint32_t read1Len = uint32_t(read1Seq.baseCount);

        // ad.qs/ts now store marker-based (non-extended) forward coordinates.
        // For the CIGAR cursor, read0Start = ad.qs (always forward).
        // For read1: same-strand → ad.ts, RC → targetLen - ad.te.
        const uint64_t read0Start = ad.qs;
        uint64_t read1Start;
        if (ad.isSameStrand) {
            read1Start = ad.ts;
        } else {
            read1Start = read1Len - ad.te;
            rcAlignmentsChecked++;
        }

        // Walk the CIGAR.
        cigarStore.walkRange(
            info.cigarOffset, info.cigarTokenCount,
            read0Start, read1Start,
            0, UINT32_MAX,
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (op == 2 || op == 3) return; // skip indels

                for (uint32_t b = 0; b < len; b++) {
                    uint32_t xpos = uint32_t(xk) + b;
                    uint32_t ypos = uint32_t(yk) + b;

                    REQUIRE(xpos < read0Seq.baseCount);

                    // read0 is always strand 0: direct lookup.
                    uint8_t base0 = read0Seq[xpos].value;

                    // read1: strand 0 → direct, strand 1 → RC lookup.
                    uint8_t base1;
                    if (ad.isSameStrand) {
                        REQUIRE(ypos < read1Seq.baseCount);
                        base1 = read1Seq[ypos].value;
                    } else {
                        REQUIRE(ypos < read1Seq.baseCount);
                        base1 = read1Seq[read1Len - 1 - ypos].complement().value;
                    }

                    if (op == 0) { // match
                        CHECK(base0 == base1);
                    } else { // mismatch
                        CHECK(base0 != base1);
                    }
                    basesVerified++;
                }
            });

        alignmentsChecked++;
    }

    // Verify we actually tested something meaningful.
    INFO("Alignments checked: " << alignmentsChecked);
    INFO("Bases verified: " << basesVerified);
    INFO("RC alignments checked: " << rcAlignmentsChecked);
    REQUIRE(alignmentsChecked > 0);
    REQUIRE(basesVerified > 100);
    // We should have at least one RC alignment (read2 is RC of shared).
    CHECK(rcAlignmentsChecked > 0);
}

TEST_CASE("Diagnostic: determine correct read1Start for RC overlaps",
          "[phasing][coordinates][diagnostic]")
{
    const std::string shared = randomSequence(2000, 88888);
    const std::string r0 = randomSequence(300, 30001) + shared + randomSequence(300, 30002);
    const std::string r1 = randomSequence(350, 30003) + reverseComplement(shared) + randomSequence(350, 30004);

    const std::vector<std::string> seqs = {r0, r1};
    TempDir tmp("dinara_diag_test");
    writeFastq(tmp.path / "reads.fastq", seqs);

    std::unique_ptr<Assembler> assembler;
    {
        ScopedSilence silence;
        auto prevDir = fs::current_path();
        fs::current_path(tmp.path);
        assembler = std::make_unique<Assembler>(tmp.path.string() + "/", true, 0, 4096);
        assembler->addReads((tmp.path / "reads.fastq").string(), 0, true, 1);
        assembler->computeReadIdsSortedByName();
        assembler->findMarkersSimdMinimizers(1, 16, 4, /*useHifiasm*/ false);
        assembler->countKmersFromMarkerKmerIds(1);
        assembler->applyKmerCountFilter(1, 1000, 1);
        OverlapCandidatesOptions candOpts;
        candOpts.method = "InvertedIndex";
        candOpts.driftRateTolerance = 0.1;
        candOpts.minChainMarkerCount = 4;
        assembler->findAlignmentCandidatesInvertedIndex(0.1, 100, candOpts, 1);
        AlignOptions alignOpts;
        alignOpts.alignMethod = 6;
        alignOpts.maxSkip = 100; alignOpts.maxDrift = 100; alignOpts.maxTrim = 10000;
        alignOpts.minAlignedMarkerCount = 2; alignOpts.minAlignedFraction = 0.0;
        alignOpts.maxMarkerFrequency = 1000;
        alignOpts.matchScore = 3; alignOpts.mismatchScore = -1; alignOpts.gapScore = -1;
        alignOpts.downsamplingFactor = 0.1; alignOpts.bandExtend = 10; alignOpts.maxBand = 1000;
        alignOpts.sameChannelReadAlignmentSuppressDeltaThreshold = 0;
        alignOpts.suppressContainments = false;
        alignOpts.align4DeltaX = 200; alignOpts.align4DeltaY = 10;
        alignOpts.align4MinEntryCountPerCell = 10; alignOpts.align4MaxDistanceFromBoundary = 100;
        alignOpts.align5DriftRateTolerance = 0.02; alignOpts.align5MinBandExtend = 10;
        alignOpts.maxErrorRate = 0.3;
        alignOpts.overlapDpMatchScore = 2; alignOpts.overlapDpMismatchScore = -4;
        alignOpts.overlapDpGapOpen1 = 4; alignOpts.overlapDpGapExtend1 = 2;
        alignOpts.overlapDpGapOpen2 = 24; alignOpts.overlapDpGapExtend2 = 1;
        assembler->computeBaseAlignmentsAndStore(alignOpts, 1);
        fs::current_path(prevDir);
    }

    const auto& cigarStore = assembler->getOverlapCigarStore();
    const auto& reads = assembler->getReads();

    for (size_t ai = 0; ai < assembler->alignmentData.size(); ai++) {
        const auto& ad = assembler->alignmentData[ai];
        if (ad.isSameStrand) continue;
        if (ad.info.cigarOffset == uint32_t(-1)) continue;
        if (ad.info.cigarTokenCount == 0) continue;

        const auto read0Seq = reads.getRead(ad.readIds[0]);
        const auto read1Seq = reads.getRead(ad.readIds[1]);
        const uint32_t read1Len = uint32_t(read1Seq.baseCount);

        // Compute oriented read1Start from ad.ts/te.
        uint64_t orientedR1Start = read1Len - ad.te; // RC conversion

        // Try oriented start with forward lookup.
        size_t matchesFwd = 0, totalFwd = 0;
        cigarStore.walkRange(ad.info.cigarOffset, ad.info.cigarTokenCount,
            ad.qs, orientedR1Start, 0, UINT32_MAX,
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (op != 0 && op != 1) return;
                for (uint32_t b = 0; b < len; b++) {
                    uint32_t xpos = uint32_t(xk) + b;
                    uint32_t ypos = uint32_t(yk) + b;
                    if (xpos >= read0Seq.baseCount || ypos >= read1Seq.baseCount) continue;
                    uint8_t b0 = read0Seq[xpos].value;
                    // Try forward lookup
                    uint8_t b1 = read1Seq[ypos].value;
                    if (op == 0 && b0 == b1) matchesFwd++;
                    totalFwd++;
                }
            });

        // Try oriented start with RC lookup.
        size_t matchesFwdRc = 0, totalFwdRc = 0;
        cigarStore.walkRange(ad.info.cigarOffset, ad.info.cigarTokenCount,
            ad.qs, orientedR1Start, 0, UINT32_MAX,
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (op != 0 && op != 1) return;
                for (uint32_t b = 0; b < len; b++) {
                    uint32_t xpos = uint32_t(xk) + b;
                    uint32_t ypos = uint32_t(yk) + b;
                    if (xpos >= read0Seq.baseCount || ypos >= read1Seq.baseCount) continue;
                    uint8_t b0 = read0Seq[xpos].value;
                    uint8_t b1 = read1Seq[read1Len - 1 - ypos].complement().value;
                    if (op == 0 && b0 == b1) matchesFwdRc++;
                    totalFwdRc++;
                }
            });

        // Try ad.ts (forward) directly with RC lookup.
        size_t matchesRcRc = 0, totalRcRc = 0;
        cigarStore.walkRange(ad.info.cigarOffset, ad.info.cigarTokenCount,
            ad.qs, ad.ts, 0, UINT32_MAX,
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (op != 0 && op != 1) return;
                for (uint32_t b = 0; b < len; b++) {
                    uint32_t xpos = uint32_t(xk) + b;
                    uint32_t ypos = uint32_t(yk) + b;
                    if (xpos >= read0Seq.baseCount || ypos >= read1Seq.baseCount) continue;
                    uint8_t b0 = read0Seq[xpos].value;
                    uint8_t b1 = read1Seq[read1Len - 1 - ypos].complement().value;
                    if (op == 0 && b0 == b1) matchesRcRc++;
                    totalRcRc++;
                }
            });

        // Try ad.ts (forward) with forward lookup.
        size_t matchesRcFwd = 0, totalRcFwd = 0;
        cigarStore.walkRange(ad.info.cigarOffset, ad.info.cigarTokenCount,
            ad.qs, ad.ts, 0, UINT32_MAX,
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (op != 0 && op != 1) return;
                for (uint32_t b = 0; b < len; b++) {
                    uint32_t xpos = uint32_t(xk) + b;
                    uint32_t ypos = uint32_t(yk) + b;
                    if (xpos >= read0Seq.baseCount || ypos >= read1Seq.baseCount) continue;
                    uint8_t b0 = read0Seq[xpos].value;
                    uint8_t b1 = read1Seq[ypos].value;
                    if (op == 0 && b0 == b1) matchesRcFwd++;
                    totalRcFwd++;
                }
            });

        INFO("ad.ts=" << ad.ts << " ad.te=" << ad.te
             << " read1Len=" << read1Len
             << " orientedR1Start=" << orientedR1Start);
        INFO("oriented+fwd: " << matchesFwd << "/" << totalFwd);
        INFO("oriented+rc:  " << matchesFwdRc << "/" << totalFwdRc);
        INFO("adTs+rc:      " << matchesRcRc << "/" << totalRcRc);
        INFO("adTs+fwd:     " << matchesRcFwd << "/" << totalRcFwd);

        // Exactly one combination should have near-100% match rate.
        // The correct one will have matchCount == totalCount for match ops.
        size_t bestMatches = std::max({matchesFwd, matchesFwdRc, matchesRcRc, matchesRcFwd});
        size_t bestTotal = 0;
        std::string bestLabel;
        if (bestMatches == matchesFwdRc) { bestTotal = totalFwdRc; bestLabel = "oriented+rc"; }
        else if (bestMatches == matchesFwd) { bestTotal = totalFwd; bestLabel = "oriented+fwd"; }
        else if (bestMatches == matchesRcRc) { bestTotal = totalRcRc; bestLabel = "adTs+rc"; }
        else { bestTotal = totalRcFwd; bestLabel = "adTs+fwd"; }

        INFO("Best: " << bestLabel << " (" << bestMatches << "/" << bestTotal << ")");
        // The correct combination should match nearly all bases.
        REQUIRE(bestTotal > 0);
        CHECK(bestMatches == bestTotal);
    }
}

TEST_CASE("Integration: ad.ts stores marker-based forward coords, RC conversion works",
          "[phasing][coordinates][integration][cigarRead1Start]")
{
    // Same setup as above but focused on verifying the coordinate conversion.
    const std::string shared = randomSequence(2000, 77777);
    const std::string r0 = randomSequence(300, 20001) + shared + randomSequence(300, 20002);
    const std::string r1 = randomSequence(350, 20003) + reverseComplement(shared) + randomSequence(350, 20004);

    const std::vector<std::string> seqs = {r0, r1};

    TempDir tmp("dinara_cigarstart_test");
    writeFastq(tmp.path / "reads.fastq", seqs);

    std::unique_ptr<Assembler> assembler;
    {
        ScopedSilence silence;
        auto prevDir = fs::current_path();
        fs::current_path(tmp.path);

        assembler = std::make_unique<Assembler>(
            tmp.path.string() + "/", true, 0, 4096);
        assembler->addReads((tmp.path / "reads.fastq").string(), 0, true, 1);
        assembler->computeReadIdsSortedByName();
        assembler->findMarkersSimdMinimizers(1, 16, 4, /*useHifiasm*/ false);
        assembler->countKmersFromMarkerKmerIds(1);
        assembler->applyKmerCountFilter(1, 1000, 1);

        OverlapCandidatesOptions candOpts;
        candOpts.method = "InvertedIndex";
        candOpts.driftRateTolerance = 0.1;
        candOpts.minChainMarkerCount = 4;
        assembler->findAlignmentCandidatesInvertedIndex(0.1, 100, candOpts, 1);

        AlignOptions alignOpts;
        alignOpts.alignMethod = 6;
        alignOpts.maxSkip = 100;
        alignOpts.maxDrift = 100;
        alignOpts.maxTrim = 10000;
        alignOpts.minAlignedMarkerCount = 2;
        alignOpts.minAlignedFraction = 0.0;
        alignOpts.maxMarkerFrequency = 1000;
        alignOpts.matchScore = 3;
        alignOpts.mismatchScore = -1;
        alignOpts.gapScore = -1;
        alignOpts.downsamplingFactor = 0.1;
        alignOpts.bandExtend = 10;
        alignOpts.maxBand = 1000;
        alignOpts.sameChannelReadAlignmentSuppressDeltaThreshold = 0;
        alignOpts.suppressContainments = false;
        alignOpts.align4DeltaX = 200;
        alignOpts.align4DeltaY = 10;
        alignOpts.align4MinEntryCountPerCell = 10;
        alignOpts.align4MaxDistanceFromBoundary = 100;
        alignOpts.align5DriftRateTolerance = 0.02;
        alignOpts.align5MinBandExtend = 10;
        alignOpts.maxErrorRate = 0.3;
        alignOpts.overlapDpMatchScore = 2;
        alignOpts.overlapDpMismatchScore = -4;
        alignOpts.overlapDpGapOpen1 = 4;
        alignOpts.overlapDpGapExtend1 = 2;
        alignOpts.overlapDpGapOpen2 = 24;
        alignOpts.overlapDpGapExtend2 = 1;

        assembler->computeBaseAlignmentsAndStore(alignOpts, 1);
        fs::current_path(prevDir);
    }

    const auto& reads = assembler->getReads();
    bool foundRcAlignment = false;

    for (size_t ai = 0; ai < assembler->alignmentData.size(); ai++) {
        const auto& ad = assembler->alignmentData[ai];
        if (ad.isSameStrand) continue;
        if (ad.info.cigarOffset == uint32_t(-1)) continue;
        if (ad.info.cigarTokenCount == 0) continue;

        foundRcAlignment = true;

        const uint32_t read1Len = uint32_t(
            reads.getRead(ad.readIds[1]).baseCount);

        // ad.ts and ad.te are marker-based forward-strand coordinates.
        CHECK(ad.ts < ad.te);
        CHECK(ad.te <= read1Len);

        // For RC overlaps, convert to oriented (RC) coordinates.
        uint64_t orientedStart = read1Len - ad.te;
        uint64_t orientedEnd = read1Len - ad.ts;

        // The oriented interval length should match the forward interval.
        CHECK((orientedEnd - orientedStart) == (ad.te - ad.ts));

        // Verify: walking with the oriented start produces valid yk values
        // that stay within [orientedStart, orientedEnd).
        const auto& cigarStore = assembler->getOverlapCigarStore();
        uint64_t minYk = UINT64_MAX, maxYk = 0;
        cigarStore.walkRange(
            ad.info.cigarOffset, ad.info.cigarTokenCount,
            ad.qs, orientedStart,
            0, UINT32_MAX,
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (yk < minYk) minYk = yk;
                uint64_t ykEnd = yk;
                // Target (yk) advances on match/mismatch and deletion (op3).
                if (opConsumesTarget(op)) ykEnd += len;
                if (ykEnd > maxYk) maxYk = ykEnd;
            });

        // yk values should be within the oriented interval.
        CHECK(minYk >= orientedStart);
        CHECK(maxYk <= orientedEnd);

        // WRONG: if we used ad.ts directly, yk would be in forward coords.
        uint64_t wrongMinYk = UINT64_MAX, wrongMaxYk = 0;
        cigarStore.walkRange(
            ad.info.cigarOffset, ad.info.cigarTokenCount,
            ad.qs, ad.ts,
            0, UINT32_MAX,
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (yk < wrongMinYk) wrongMinYk = yk;
                uint64_t ykEnd = yk;
                if (opConsumesTarget(op)) ykEnd += len;
                if (ykEnd > wrongMaxYk) wrongMaxYk = ykEnd;
            });

        // The wrong start produces yk in [ad.ts, ad.te) — forward coords.
        // These are different from the oriented coords.
        if (orientedStart != ad.ts) {
            CHECK(wrongMinYk != minYk);
        }
    }

    REQUIRE(foundRcAlignment);
}
