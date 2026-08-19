/**
 * @file test_paf_import.cpp
 * @brief Unit tests for the overlap import helpers in src/PafImport.hpp.
 *
 * These tests exercise the pure, allocation-free building blocks used by
 * Assembler::importAlignmentCandidatesFromMemory:
 *   - makePafEntry / pafEntryLess / dedupPafEntriesKeepBestScore: canonical
 *     keying and highest-chain-score-wins deduplication (matching hifiasm's
 *     oreg_ss_lt selection), independent of input order.
 *
 * They also pin the semantics of the in-memory hifiasm overlap path (canonical
 * keying, both-orientation retention, block-length floor, self-overlap filter).
 */

#include "catch.hpp"
#include "PafImport.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

using namespace dinara;
using std::vector;

namespace {

// Mirror of one hifiasm in-memory overlap (see hifiasm_overlap_t in
// external/hifiasmCandidates/hifiasm_overlaps.h). Duplicated here so the pure
// PafImport unit tests do not need to link the submodule; the field semantics
// must stay in sync with that header.
struct MemOverlap {
    uint32_t q_id, t_id;
    uint32_t q_start, q_end;
    uint32_t t_start, t_end;
    uint32_t n_match, block_len;
    uint32_t shared_seed;
    bool     is_same_strand;
};

// Build a deduplicated PafEntry set from in-memory overlaps, exactly as
// Assembler::importAlignmentCandidatesFromMemory does (same filters, same
// makePafEntry call, same dedup). Read ids are the overlap's own q_id/t_id here
// (the real importer resolves them by name first, which is orthogonal to the
// dedup/keying property under test).
vector<PafEntry> importMem(
    const vector<MemOverlap>& ovs, uint64_t minOverlapLength = 1000)
{
    vector<PafEntry> entries;
    for (const MemOverlap& o : ovs) {
        // Single overlap-length gate: block_len (span between the first and last
        // shared minimizer) must reach minOverlapLength. Mirrors production.
        if (o.block_len < minOverlapLength) continue;
        if (o.q_id == o.t_id) continue;
        entries.push_back(makePafEntry(
            ReadId(o.q_id), ReadId(o.t_id),
            o.q_start, o.q_end, o.t_start, o.t_end,
            o.block_len, o.shared_seed, o.is_same_strand));
    }
    dedupPafEntriesKeepBestScore(entries);
    return entries;
}

} // namespace


TEST_CASE("makePafEntry: canonicalizes read pair and swaps intervals", "[paf]") {
    // readId0 > readId1 -> must swap so key uses (min<<32|max) and q* = min id.
    PafEntry e = makePafEntry(
        /*readId0=*/5, /*readId1=*/2,
        /*qStart=*/100, /*qEnd=*/900,   // these are read 5's coords
        /*tStart=*/50,  /*tEnd=*/850,   // these are read 2's coords
        /*blockLen=*/700, /*sharedSeedScore=*/700, /*isSameStrand=*/true);
    REQUIRE(e.key == ((uint64_t(2) << 32) | 5));
    // After swap: q* refers to read 2 (=old target), t* to read 5 (=old query).
    REQUIRE(e.iv.qStart == 50);
    REQUIRE(e.iv.qEnd == 850);
    REQUIRE(e.iv.tStart == 100);
    REQUIRE(e.iv.tEnd == 900);
    REQUIRE(e.iv.blockLen == 700);
}

TEST_CASE("makePafEntry: reciprocal A->B and B->A produce identical entry", "[paf]") {
    // A=5, B=2. Record A->B: q=A coords, t=B coords.
    PafEntry ab = makePafEntry(5, 2, 100, 900, 50, 850, 700, 700, true);
    // Record B->A: q=B coords, t=A coords.
    PafEntry ba = makePafEntry(2, 5, 50, 850, 100, 900, 700, 700, true);
    REQUIRE(ab.key == ba.key);
    REQUIRE(ab.iv.qStart == ba.iv.qStart);
    REQUIRE(ab.iv.qEnd == ba.iv.qEnd);
    REQUIRE(ab.iv.tStart == ba.iv.tStart);
    REQUIRE(ab.iv.tEnd == ba.iv.tEnd);
}

TEST_CASE("dedupPafEntriesKeepBestScore: keeps the highest chain score per pair", "[paf]") {
    // makePafEntry args: readId0, readId1, qs, qe, ts, te, blockLen,
    // sharedSeedScore, isSameStrand. Selection is by sharedSeedScore, NOT span.
    vector<PafEntry> entries = {
        makePafEntry(1, 2, 0, 300, 0, 300, 300, 40, true),
        makePafEntry(1, 2, 0, 700, 0, 700, 700, 90, true),   // best score for (1,2)
        makePafEntry(1, 2, 0, 150, 0, 150, 150, 10, true),
        makePafEntry(3, 4, 0, 500, 0, 500, 500, 55, false),
    };
    dedupPafEntriesKeepBestScore(entries);
    REQUIRE(entries.size() == 2);
    // Sorted by key ascending: (1,2) then (3,4).
    REQUIRE(entries[0].key == ((uint64_t(1) << 32) | 2));
    REQUIRE(entries[0].iv.sharedSeedScore == 90);
    REQUIRE(entries[0].iv.blockLen == 700);
    REQUIRE(entries[0].iv.qEnd == 700);
    REQUIRE(entries[1].key == ((uint64_t(3) << 32) | 4));
    REQUIRE(entries[1].iv.sharedSeedScore == 55);
}

TEST_CASE("dedupPafEntriesKeepBestScore: highest score wins over longer span", "[paf]") {
    // The longer overlap (span 700) has the LOWER chain score; hifiasm selects
    // by shared_seed, so the shorter, higher-scoring overlap must survive.
    vector<PafEntry> entries = {
        makePafEntry(1, 2, 0, 700, 0, 700, 700, 30, true),   // long span, low score
        makePafEntry(1, 2, 0, 300, 0, 300, 300, 95, true),   // short span, high score
    };
    dedupPafEntriesKeepBestScore(entries);
    REQUIRE(entries.size() == 1);
    REQUIRE(entries[0].iv.sharedSeedScore == 95);
    REQUIRE(entries[0].iv.blockLen == 300);
}

TEST_CASE("dedupPafEntriesKeepBestScore: ties on score fall back to shorter span", "[paf]") {
    // hifiasm's per-target collapse (anchor.cpp:756) keeps the SMALLER overlapLen
    // on a shared_seed tie, so the shorter span must survive.
    vector<PafEntry> entries = {
        makePafEntry(1, 2, 0, 700, 0, 700, 700, 60, true),   // same score, longer span
        makePafEntry(1, 2, 0, 300, 0, 300, 300, 60, true),   // same score, shorter span (kept)
    };
    dedupPafEntriesKeepBestScore(entries);
    REQUIRE(entries.size() == 1);
    REQUIRE(entries[0].iv.sharedSeedScore == 60);
    REQUIRE(entries[0].iv.blockLen == 300);
}

TEST_CASE("dedupPafEntriesKeepBestScore: result is order-independent", "[paf]") {
    auto build = [](bool reversed) {
        vector<PafEntry> v = {
            makePafEntry(1, 2, 0, 300, 0, 300, 300, 40, true),
            makePafEntry(1, 2, 0, 700, 0, 700, 700, 90, true),
            makePafEntry(2, 5, 0, 400, 0, 400, 400, 55, false),
        };
        if (reversed) std::reverse(v.begin(), v.end());
        dedupPafEntriesKeepBestScore(v);
        return v;
    };
    const auto a = build(false);
    const auto b = build(true);
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        REQUIRE(a[i].key == b[i].key);
        REQUIRE(a[i].iv.sharedSeedScore == b[i].iv.sharedSeedScore);
        REQUIRE(a[i].iv.blockLen == b[i].iv.blockLen);
        REQUIRE(a[i].iv.qStart == b[i].iv.qStart);
        REQUIRE(a[i].iv.qEnd == b[i].iv.qEnd);
    }
}


// --------------------------------------------------------------------------
// In-memory overlap import semantics (the production path).
//
// The in-memory bridge hands back each overlap in its own (q=x, t=y)
// orientation. makePafEntry canonicalizes it to q*=min(readId) coordinates so
// downstream candidates are keyed identically regardless of which read hifiasm
// labeled as query.
// --------------------------------------------------------------------------

TEST_CASE("mem import: canonicalizes to min-id query coordinates", "[paf][mem]") {
    // In-memory route (record's own orientation): q_id=0 [100,900),
    // t_id=1 [50,850), same strand.
    vector<MemOverlap> mem = {
        {0, 1, 100, 900, 50, 850, 700, 800, 800, true},
    };
    // Length floor disabled: this test is about coordinate canonicalization.
    vector<PafEntry> fromMem = importMem(mem, /*minOverlapLength=*/0);

    REQUIRE(fromMem.size() == 1);
    // The canonical interval is q=min(id) coords.
    REQUIRE(fromMem[0].key == ((uint64_t(0) << 32) | 1));
    REQUIRE(fromMem[0].iv.qStart == 100);  // read 0 (min id) coords
    REQUIRE(fromMem[0].iv.qEnd   == 900);
    REQUIRE(fromMem[0].iv.tStart == 50);   // read 1 coords
    REQUIRE(fromMem[0].iv.tEnd   == 850);
    REQUIRE(fromMem[0].iv.isSameStrand == true);

    // Same overlap with q/t swapped must canonicalize identically.
    vector<MemOverlap> swapped = {
        {1, 0, 50, 850, 100, 900, 700, 800, 800, true},
    };
    vector<PafEntry> fromSwapped = importMem(swapped, /*minOverlapLength=*/0);
    REQUIRE(fromSwapped.size() == 1);
    REQUIRE(fromSwapped[0].key == fromMem[0].key);
    REQUIRE(fromSwapped[0].iv.qStart == fromMem[0].iv.qStart);
    REQUIRE(fromSwapped[0].iv.qEnd   == fromMem[0].iv.qEnd);
    REQUIRE(fromSwapped[0].iv.tStart == fromMem[0].iv.tStart);
    REQUIRE(fromSwapped[0].iv.tEnd   == fromMem[0].iv.tEnd);
}

TEST_CASE("mem import: both orientations of a pair are kept", "[paf][mem]") {
    // Reads 0 and 1 overlap both same-strand and reverse-complement (e.g. an
    // inverted repeat). Both must survive as distinct entries.
    vector<MemOverlap> mem = {
        {0, 1, 100, 900, 50, 850, 700, 800, 800, true},   // +
        {0, 1, 100, 900, 50, 850, 700, 800, 800, false},  // -
    };
    // Length floor disabled: this test is about keeping both orientations.
    vector<PafEntry> got = importMem(mem, /*minOverlapLength=*/0);
    REQUIRE(got.size() == 2);
    // Deterministic order: same-strand before reverse within a key.
    REQUIRE(got[0].key == got[1].key);
    REQUIRE(got[0].iv.isSameStrand == true);
    REQUIRE(got[1].iv.isSameStrand == false);
}

TEST_CASE("mem import: highest chain score wins per (pair, strand)", "[paf][mem]") {
    // MemOverlap fields: q_id,t_id, q_start,q_end, t_start,t_end, n_match,
    // block_len, shared_seed, is_same_strand. Selection is by shared_seed.
    // Length floor disabled so only the score-selection behavior is exercised.
    vector<MemOverlap> mem = {
        {0, 1, 100, 500, 50, 450, 380, 400, 420, true},   // lower score +
        {0, 1, 100, 900, 50, 850, 700, 800, 900, true},   // higher score + (kept)
        {2, 3, 0, 400, 0, 400, 380, 150, 150, true},      // distinct pair, kept
    };
    vector<PafEntry> got = importMem(mem, /*minOverlapLength=*/0);
    REQUIRE(got.size() == 2);
    // Entries are sorted by key ascending: (0,1) first, then (2,3).
    REQUIRE(got[0].iv.sharedSeedScore == 900);
    REQUIRE(got[0].iv.blockLen == 800);
    REQUIRE(got[0].iv.qStart == 100);
    REQUIRE(got[0].iv.qEnd   == 900);
    REQUIRE(got[1].key == ((uint64_t(2) << 32) | 3));
}

TEST_CASE("mem import: minOverlapLength floor and self-overlaps filtered", "[paf][mem]") {
    // The single overlap-length gate is block_len >= minOverlapLength. Using the
    // production default (1000), an overlap whose block_len is one below is
    // dropped and one exactly at the threshold is kept.
    vector<MemOverlap> mem = {
        {0, 0, 0, 2000, 0, 2000, 380, 1500, 1500, true},  // self overlap: dropped
        {0, 1, 0, 2000, 0, 2000, 380, 999,  999,  true},  // below floor: dropped
        {0, 1, 0, 2000, 0, 2000, 380, 1000, 1000, true},  // exactly floor: kept
    };
    vector<PafEntry> got = importMem(mem);  // default minOverlapLength = 1000
    REQUIRE(got.size() == 1);
    REQUIRE(got[0].iv.blockLen == 1000);
}

TEST_CASE("mem import: minOverlapLength=0 disables the length floor", "[paf][mem]") {
    // With the gate disabled, even a tiny overlap survives (self-overlaps are
    // still dropped independently of the length gate).
    vector<MemOverlap> mem = {
        {0, 0, 0, 400, 0, 400, 380, 50, 50, true},   // self overlap: dropped
        {0, 1, 0, 400, 0, 400, 380, 50, 50, true},   // tiny, but kept (floor off)
    };
    vector<PafEntry> got = importMem(mem, /*minOverlapLength=*/0);
    REQUIRE(got.size() == 1);
    REQUIRE(got[0].iv.blockLen == 50);
}
