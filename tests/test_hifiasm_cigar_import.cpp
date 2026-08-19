/**
 * @file test_hifiasm_cigar_import.cpp
 * @brief Unit tests for normalizeHifiasmCigar in src/HifiasmCigarImport.hpp.
 *
 * hifiasm emits one CIGAR per overlap in the ALIGNMENT frame (query forward,
 * target in alignment orientation; op I consumes query, op D consumes target).
 * dinara canonicalizes every pair to read0 = min(ReadId), read1 = max(ReadId),
 * read0 forward, read1 in alignment orientation. normalizeHifiasmCigar reframes
 * the hifiasm CIGAR into that canonical frame.
 *
 * These tests build a synthetic overlap with a known op script, materialize the
 * two reads (with flanks so anchors are non-zero), normalize for all four
 * cases, and verify that every match/mismatch column pairs the correct bases
 * and that indels advance the correct side:
 *   (query<target | query>target) x (same-strand | reverse-strand).
 *
 * Pure logic over the header-only OverlapCigarStore, so this target links only
 * DINARA_ASSERT.cpp (no astarpa/poasta/shasta2).
 */

#include "catch.hpp"
#include "HifiasmCigarImport.hpp"
#include "HifiasmImportedCigarStore.hpp"
#include "PafImport.hpp"

#include <string>
#include <vector>

using namespace dinara;
using std::string;
using std::vector;

namespace {

char complementBase(char b) {
    switch(b) { case 'A': return 'T'; case 'C': return 'G';
                case 'G': return 'C'; case 'T': return 'A'; }
    return 'N';
}
string reverseComplement(const string& s) {
    string o; o.reserve(s.size());
    for(size_t i = s.size(); i-- > 0; ) o.push_back(complementBase(s[i]));
    return o;
}

struct Col { uint8_t op; uint32_t len; };

vector<CigarToken> toTokens(const vector<Col>& script) {
    vector<CigarToken> t;
    for(const auto& c : script) t.emplace_back(c.op, uint16_t(c.len));
    return t;
}

// Walk canonical tokens from (r0s,r1s) and check base correspondence against
// read0 forward and read1 in alignment orientation. Asserts expected op counts.
void verifyWalk(const vector<CigarToken>& toks,
                uint32_t r0s, uint32_t r1s,
                const string& read0fwd, const string& read1aln,
                int expM, int expMM, int expIns, int expDel)
{
    uint64_t x = r0s, y = r1s;
    int m = 0, mm = 0, ins = 0, del = 0;
    for(auto tk : toks) {
        const uint8_t op = tk.op();
        const uint32_t len = tk.len();
        for(uint32_t k = 0; k < len; ++k) {
            if(op == CigarOpMatch)      { REQUIRE(read0fwd[x] == read1aln[y]); ++m;  ++x; ++y; }
            else if(op == CigarOpMismatch) { REQUIRE(read0fwd[x] != read1aln[y]); ++mm; ++x; ++y; }
            else if(op == CigarOpIns)   { ++ins; ++x; }   // consumes read0
            else if(op == CigarOpDel)   { ++del; ++y; }   // consumes read1
        }
    }
    REQUIRE(m == expM);
    REQUIRE(mm == expMM);
    REQUIRE(ins == expIns);
    REQUIRE(del == expDel);
}

// Shared synthetic overlap. Script columns in the hifiasm (query,target) frame:
//   5 M, 1 X, 2 M, 1 I(query-only), 3 M, 1 D(target-only), 2 M
// query-consuming len = 5+1+2+1+3+2 = 14 ; target-consuming = 5+1+2+3+1+2 = 14
struct Overlap {
    vector<Col> script;
    string qOv;        // query overlap bases
    string tAln;       // target bases in alignment orientation
    int EM, EMM, EINS, EDEL;
    Overlap() {
        script = {
            {CigarOpMatch,5},{CigarOpMismatch,1},{CigarOpMatch,2},
            {CigarOpIns,1},{CigarOpMatch,3},{CigarOpDel,1},{CigarOpMatch,2}
        };
        auto pick = [](int seed){ const char* B = "ACGT"; return B[seed & 3]; };
        int s = 0;
        for(const auto& c : script) {
            for(uint32_t k = 0; k < c.len; ++k) {
                const char b = pick(s++);
                if(c.op == CigarOpMatch)         { qOv.push_back(b);   tAln.push_back(b); }
                else if(c.op == CigarOpMismatch) { qOv.push_back('A'); tAln.push_back('C'); }
                else if(c.op == CigarOpIns)      { qOv.push_back(b); }
                else if(c.op == CigarOpDel)      { tAln.push_back(b); }
            }
        }
        EM = 5 + 2 + 3 + 2; EMM = 1; EINS = 1; EDEL = 1;
    }
};

} // namespace

TEST_CASE("normalizeHifiasmCigar: query<target, same strand (no reframe)",
          "[hifiasm][cigar][import]") {
    Overlap ov;
    const string qL = "GGGG", qR = "TTT", tL = "CC", tR = "AAAAA";
    const uint32_t qStart = qL.size(), qEnd = qStart + ov.qOv.size();
    const uint32_t tStart = tL.size(), tEnd = tStart + ov.tAln.size();
    const string qFull = qL + ov.qOv + qR;
    const string tFull = tL + ov.tAln + tR;

    auto tokens = toTokens(ov.script);
    auto n = normalizeHifiasmCigar(
        span<const CigarToken>(tokens.data(), tokens.size()),
        /*readIdQ*/1, /*readIdT*/2, qStart, qEnd, tStart, tEnd,
        uint32_t(qFull.size()), uint32_t(tFull.size()), /*isSameStrand*/true);

    REQUIRE(n.read0Start == qStart); REQUIRE(n.read0End == qEnd);
    REQUIRE(n.read1Start == tStart); REQUIRE(n.read1End == tEnd);
    verifyWalk(n.tokens, n.read0Start, n.read1Start, qFull, tFull,
               ov.EM, ov.EMM, ov.EINS, ov.EDEL);
}

TEST_CASE("normalizeHifiasmCigar: query<target, reverse strand (read1 RC frame)",
          "[hifiasm][cigar][import]") {
    Overlap ov;
    const string qL = "GGGG", qR = "TTT", tL = "CC", tR = "AAAAA";
    const uint32_t qStart = qL.size(), qEnd = qStart + ov.qOv.size();
    const uint32_t tStart = tL.size(), tEnd = tStart + ov.tAln.size();
    // Target forward has rc(tAln) at [tStart,tEnd]; alignment frame = rc(tForward).
    const string tForward = tL + reverseComplement(ov.tAln) + tR;
    const string qFull = qL + ov.qOv + qR;
    const string read1aln = reverseComplement(tForward);

    auto tokens = toTokens(ov.script);
    auto n = normalizeHifiasmCigar(
        span<const CigarToken>(tokens.data(), tokens.size()),
        1, 2, qStart, qEnd, tStart, tEnd,
        uint32_t(qFull.size()), uint32_t(tForward.size()), /*isSameStrand*/false);

    REQUIRE(n.read0Start == qStart); REQUIRE(n.read0End == qEnd);
    REQUIRE(n.read1Start == uint32_t(tForward.size()) - tEnd);
    REQUIRE(n.read1End   == uint32_t(tForward.size()) - tStart);
    verifyWalk(n.tokens, n.read0Start, n.read1Start, qFull, read1aln,
               ov.EM, ov.EMM, ov.EINS, ov.EDEL);
}

TEST_CASE("normalizeHifiasmCigar: query>target, same strand (transpose)",
          "[hifiasm][cigar][import]") {
    Overlap ov;
    const string qL = "GGGG", qR = "TTT", tL = "CC", tR = "AAAAA";
    const uint32_t qStart = qL.size(), qEnd = qStart + ov.qOv.size();
    const uint32_t tStart = tL.size(), tEnd = tStart + ov.tAln.size();
    const string qFull = qL + ov.qOv + qR;
    const string tFull = tL + ov.tAln + tR;

    auto tokens = toTokens(ov.script);
    // readIdQ > readIdT: dinara read0 = target, read1 = query.
    auto n = normalizeHifiasmCigar(
        span<const CigarToken>(tokens.data(), tokens.size()),
        /*readIdQ*/9, /*readIdT*/2, qStart, qEnd, tStart, tEnd,
        uint32_t(qFull.size()), uint32_t(tFull.size()), /*isSameStrand*/true);

    REQUIRE(n.read0Start == tStart); REQUIRE(n.read0End == tEnd);
    REQUIRE(n.read1Start == qStart); REQUIRE(n.read1End == qEnd);
    // Transpose swaps I<->D: query insertions become target deletions.
    verifyWalk(n.tokens, n.read0Start, n.read1Start, tFull, qFull,
               ov.EM, ov.EMM, /*ins*/ov.EDEL, /*del*/ov.EINS);
}

TEST_CASE("normalizeHifiasmCigar: query>target, reverse strand (reverse+transpose)",
          "[hifiasm][cigar][import]") {
    Overlap ov;
    const string qL = "GGGG", qR = "TTT", tL = "CC", tR = "AAAAA";
    const uint32_t qStart = qL.size(), qEnd = qStart + ov.qOv.size();
    const uint32_t tStart = tL.size(), tEnd = tStart + ov.tAln.size();
    const string tForward = tL + reverseComplement(ov.tAln) + tR;
    const string qFull = qL + ov.qOv + qR;
    const string read0fwd = tForward;              // read0 = target forward
    const string read1aln = reverseComplement(qFull); // read1 = query RC

    auto tokens = toTokens(ov.script);
    auto n = normalizeHifiasmCigar(
        span<const CigarToken>(tokens.data(), tokens.size()),
        9, 2, qStart, qEnd, tStart, tEnd,
        uint32_t(qFull.size()), uint32_t(tForward.size()), /*isSameStrand*/false);

    REQUIRE(n.read0Start == tStart); REQUIRE(n.read0End == tEnd);
    REQUIRE(n.read1Start == uint32_t(qFull.size()) - qEnd);
    REQUIRE(n.read1End   == uint32_t(qFull.size()) - qStart);
    verifyWalk(n.tokens, n.read0Start, n.read1Start, read0fwd, read1aln,
               ov.EM, ov.EMM, /*ins*/ov.EDEL, /*del*/ov.EINS);
}

// ---------------------------------------------------------------------------
// Association: the CIGAR stored for a (key,strand) must belong to the overlap
// whose interval survives dedup, and the key the consumer reconstructs from the
// candidate must equal the key used at store time.
//
// importAlignmentCandidatesFromMemory tags each PafEntry with sourceIndex (the
// overlap it came from), runs dedupPafEntriesKeepBestScore, then builds the
// CIGAR store from the survivors. computeBaseAlignmentsAndStore rebuilds the
// lookup key as (candidate.readIds[0]<<32)|readIds[1] where candidate comes from
// OrientedReadPair(key>>32, key&mask, isSameStrand). These tests pin both the
// "best score wins, sourceIndex rides along" and the "key round-trips"
// invariants.
// ---------------------------------------------------------------------------

TEST_CASE("import association: dedup keeps the best-scoring overlap's sourceIndex",
          "[hifiasm][cigar][import][association]") {
    using namespace dinara;
    // Three overlaps on the same pair+strand with increasing chain score; the
    // highest-scoring (sourceIndex 2) must be the survivor whose CIGAR we store.
    std::vector<PafEntry> entries;
    for(uint32_t i = 0; i < 3; ++i) {
        PafEntry e = makePafEntry(/*readId0*/5, /*readId1*/9,
                                  /*q*/0, /*qEnd*/10 + i, /*t*/0, /*tEnd*/10 + i,
                                  /*blockLen*/10 + i, /*sharedSeedScore*/10 + i,
                                  /*isSameStrand*/true);
        e.sourceIndex = i;
        entries.push_back(e);
    }
    // A different-strand overlap on the same pair survives independently.
    {
        PafEntry e = makePafEntry(5, 9, 0, 7, 0, 7, 7, 7, /*isSameStrand*/false);
        e.sourceIndex = 99;
        entries.push_back(e);
    }

    dedupPafEntriesKeepBestScore(entries);

    // One survivor per (key,strand).
    const uint64_t key = (uint64_t(5) << 32) | uint64_t(9);
    const PafEntry* same = nullptr;
    const PafEntry* diff = nullptr;
    for(const auto& e : entries) {
        REQUIRE(e.key == key);
        if(e.iv.isSameStrand) { REQUIRE(same == nullptr); same = &e; }
        else                  { REQUIRE(diff == nullptr); diff = &e; }
    }
    REQUIRE(same != nullptr);
    REQUIRE(diff != nullptr);
    REQUIRE(same->sourceIndex == 2);   // highest-scoring same-strand overlap
    REQUIRE(same->iv.sharedSeedScore == 12);
    REQUIRE(same->iv.blockLen == 12);
    REQUIRE(diff->sourceIndex == 99);  // the lone reverse overlap
}

TEST_CASE("import association: makePafEntry key survives the candidate round-trip",
          "[hifiasm][cigar][import][association]") {
    using namespace dinara;
    // Regardless of which read hifiasm called query vs target, the canonical
    // key and the interval orientation are stable, and the consumer's
    // (readIds[0]<<32)|readIds[1] over the canonical pair equals that key.
    struct Case { uint32_t a, b; bool sameStrand; };
    const Case cases[] = { {5, 9, true}, {9, 5, true}, {5, 9, false}, {9, 5, false} };

    for(const auto& c : cases) {
        PafEntry e = makePafEntry(c.a, c.b, /*q*/1, 11, /*t*/2, 12, 10, 10, c.sameStrand);

        const uint32_t lo = std::min(c.a, c.b);
        const uint32_t hi = std::max(c.a, c.b);
        const uint64_t expectedKey = (uint64_t(lo) << 32) | uint64_t(hi);
        REQUIRE(e.key == expectedKey);

        // Consumer side: OrientedReadPair(key>>32, key&mask) yields readIds
        // {lo,hi}; the rebuilt lookup key must match the stored key.
        const uint32_t readId0 = uint32_t(e.key >> 32);
        const uint32_t readId1 = uint32_t(e.key & 0xffffffffULL);
        REQUIRE(readId0 == lo);
        REQUIRE(readId1 == hi);
        const uint64_t rebuiltKey = (uint64_t(readId0) << 32) | uint64_t(readId1);
        REQUIRE(rebuiltKey == e.key);

        // Strand is carried verbatim (same map bucket at store and lookup).
        REQUIRE(e.iv.isSameStrand == c.sameStrand);
    }
}

// HifiasmImportedCigarStore::add() must transpose op2<->op3 at ingest.
//
// hifiasm's exported CIGAR uses op2 = target-consuming and op3 = query-consuming
// (the transpose of dinara's OverlapCigarStore convention, where op2/CigarOpIns
// consumes the query and op3/CigarOpDel consumes the target). The store fixes
// this at the single ingest boundary so the stored tokens follow dinara's
// convention and stay consistent with the recorded qStart/qEnd/tStart/tEnd.
//
// This test would fail if the transpose in add() were removed: the stored
// query/target consumption would no longer match the recorded overlap spans.
TEST_CASE("HifiasmImportedCigarStore transposes op2/op3 at ingest") {
    // A raw hifiasm CIGAR (its frame): 5 match, then op2 x3, then op3 x2, 4 match.
    // Raw op2 consumes TARGET, raw op3 consumes QUERY. So on the raw tokens:
    //   query-consuming = op0(5) + op3(2) + op0(4) = 11  -> qSpan
    //   target-consuming = op0(5) + op2(3) + op0(4) = 12 -> tSpan
    const vector<uint16_t> rawTokens = {
        CigarToken(CigarOpMatch, 5).data,
        CigarToken(CigarOpIns,   3).data,   // raw op2 = target-consuming
        CigarToken(CigarOpDel,   2).data,   // raw op3 = query-consuming
        CigarToken(CigarOpMatch, 4).data,
    };
    const uint32_t qStart = 0, qEnd = 11;   // query span (raw op3 + matches)
    const uint32_t tStart = 0, tEnd = 12;   // target span (raw op2 + matches)

    HifiasmImportedCigarStore store;
    const uint64_t pairKey = (uint64_t(0) << 32) | uint64_t(1);
    store.add(pairKey, /*isSameStrand*/ true,
              span<const uint16_t>(rawTokens.data(), rawTokens.size()),
              /*readIdQ*/ 0, /*readIdT*/ 1, qStart, qEnd, tStart, tEnd);

    const auto* rec = store.find(pairKey, /*isSameStrand*/ true);
    REQUIRE(rec != nullptr);

    // Sum query/target consumption of the STORED tokens using dinara's helpers.
    uint64_t storedQ = 0, storedT = 0;
    OverlapCigarStore::forEachOp(store.tokensOf(*rec),
        [&](uint8_t op, uint32_t len){
            if(opConsumesQuery(op))  storedQ += len;
            if(opConsumesTarget(op)) storedT += len;
        });

    // After the ingest transpose the stored tokens must be consistent with the
    // recorded spans under dinara's convention.
    REQUIRE(storedQ == qEnd - qStart);   // 11
    REQUIRE(storedT == tEnd - tStart);   // 12

    // Concretely, the two indel ops must have been swapped: op2 (raw target)
    // became CigarOpDel and op3 (raw query) became CigarOpIns.
    vector<CigarToken> stored(store.tokensOf(*rec).begin(), store.tokensOf(*rec).end());
    REQUIRE(stored.size() == 4);
    REQUIRE(stored[0].op() == CigarOpMatch); REQUIRE(stored[0].len() == 5);
    REQUIRE(stored[1].op() == CigarOpDel);   REQUIRE(stored[1].len() == 3);
    REQUIRE(stored[2].op() == CigarOpIns);   REQUIRE(stored[2].len() == 2);
    REQUIRE(stored[3].op() == CigarOpMatch); REQUIRE(stored[3].len() == 4);
}
