/**
 * @file test_overlap_cigar_store.cpp
 * @brief Unit tests for OverlapCigarStore.
 *
 * Tests the hifiasm-style uint16_t packed CIGAR store, covering:
 * - Token encoding/decoding
 * - All 6 segment cases from constructQuickRawSparse
 * - Multi-segment stitching
 * - Long run splitting (>16383)
 * - Thread-local merge with offset adjustment
 * - forEachOp coalescing across segment boundaries
 * - Coordinate consumption verification
 */

#include "../external/catch2/catch.hpp"
#include "../src/OverlapCigarStore.hpp"

using namespace dinara;

// Helper: walk a CIGAR and return (consumed_read0, consumed_read1).
static std::pair<uint64_t, uint64_t> consumedBases(
    const OverlapCigarStore& store, uint32_t offset, uint32_t count)
{
    uint64_t c0 = 0, c1 = 0;
    store.forEachOp(offset, count, [&](uint8_t op, uint32_t len) {
        switch(op) {
            case 0: // match
            case 1: // mismatch
                c0 += len; c1 += len; break;
            case 2: // insertion (extra in read1)
                c1 += len; break;
            case 3: // deletion (extra in read0)
                c0 += len; break;
        }
    });
    return {c0, c1};
}

// Helper: collect all (op, len) pairs from forEachOp.
static std::vector<std::pair<uint8_t, uint32_t>> collectOps(
    const OverlapCigarStore& store, uint32_t offset, uint32_t count)
{
    std::vector<std::pair<uint8_t, uint32_t>> ops;
    store.forEachOp(offset, count, [&](uint8_t op, uint32_t len) {
        ops.emplace_back(op, len);
    });
    return ops;
}

// ============================================================================
// Token encoding
// ============================================================================

TEST_CASE("CigarToken encoding and decoding", "[CigarStore]") {
    CigarToken t0(0, 100);
    REQUIRE(t0.op() == 0);
    REQUIRE(t0.len() == 100);

    CigarToken t1(1, 1);
    REQUIRE(t1.op() == 1);
    REQUIRE(t1.len() == 1);

    CigarToken t2(2, 500);
    REQUIRE(t2.op() == 2);
    REQUIRE(t2.len() == 500);

    CigarToken t3(3, CigarToken::MAX_LEN);
    REQUIRE(t3.op() == 3);
    REQUIRE(t3.len() == CigarToken::MAX_LEN);
}

// ============================================================================
// Case 1: Both lengths zero (empty segment)
// ============================================================================

TEST_CASE("Segment case 1: len0=0, len1=0 — no tokens emitted", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    // Simulate: both sides empty, nothing emitted (just like the code does).
    // No pushOp calls.

    auto tokens = store.getTokens(off, store.tokensSince(off));
    REQUIRE(tokens.size() == 0);

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 0);
    REQUIRE(c1 == 0);
}

// ============================================================================
// Case 2: len0 > 0, len1 = 0 — deletion
// ============================================================================

TEST_CASE("Segment case 2: len0>0, len1=0 — deletion", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    uint32_t len0 = 42;
    store.pushDeletion(len0);

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0].first == 3);  // deletion
    REQUIRE(ops[0].second == 42);

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 42);
    REQUIRE(c1 == 0);
}

// ============================================================================
// Case 3: len0 = 0, len1 > 0 — insertion
// ============================================================================

TEST_CASE("Segment case 3: len0=0, len1>0 — insertion", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    uint32_t len1 = 37;
    store.pushInsertion(len1);

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0].first == 2);  // insertion
    REQUIRE(ops[0].second == 37);

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 0);
    REQUIRE(c1 == 37);
}

// ============================================================================
// Case 4: len0 > 0, len1 > 0, identical sequences — match run
// ============================================================================

TEST_CASE("Segment case 4: identical sequences — match run", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    uint32_t len = 1000;
    store.pushMatch(len);

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0].first == 0);  // match
    REQUIRE(ops[0].second == 1000);

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 1000);
    REQUIRE(c1 == 1000);
}

// ============================================================================
// Case 5: len0 == len1, sequences differ — match/mismatch split
// ============================================================================

TEST_CASE("Segment case 5: same length, with mismatches", "[CigarStore]") {
    // Simulate: ACGTACGT vs ACGAACGT (mismatch at position 3)
    // The code would emit: match(3), mismatch(1), match(4)
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    store.pushMatch(3);
    store.pushMismatch(1);
    store.pushMatch(4);

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 3);
    REQUIRE(ops[0] == std::make_pair(uint8_t(0), uint32_t(3)));
    REQUIRE(ops[1] == std::make_pair(uint8_t(1), uint32_t(1)));
    REQUIRE(ops[2] == std::make_pair(uint8_t(0), uint32_t(4)));

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 8);
    REQUIRE(c1 == 8);
}

// ============================================================================
// Case 6: len0 != len1, sequences differ — match/mismatch/indel
// ============================================================================

TEST_CASE("Segment case 6: different lengths, with indels", "[CigarStore]") {
    // Simulate A*PA2 CIGAR: 50M 2I 3D 45M (len0=98, len1=97)
    // M block with all matches for simplicity.
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    store.pushMatch(50);
    store.pushInsertion(2);
    store.pushDeletion(3);
    store.pushMatch(45);

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 4);
    REQUIRE(ops[0] == std::make_pair(uint8_t(0), uint32_t(50)));
    REQUIRE(ops[1] == std::make_pair(uint8_t(2), uint32_t(2)));
    REQUIRE(ops[2] == std::make_pair(uint8_t(3), uint32_t(3)));
    REQUIRE(ops[3] == std::make_pair(uint8_t(0), uint32_t(45)));

    // read0 consumes: 50 (match) + 3 (del) + 45 (match) = 98
    // read1 consumes: 50 (match) + 2 (ins) + 45 (match) = 97
    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 98);
    REQUIRE(c1 == 97);
}

// ============================================================================
// Multi-segment stitching
// ============================================================================

TEST_CASE("Multi-segment stitching: identical + aligned + identical", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    // Segment 0: identical, 500 bases
    store.pushMatch(500);

    // Segment 1: aligned with mismatches and an insertion
    // A*PA2 produced: 100M (with 2 mismatches at pos 30 and 70) then 3I then 97M
    store.pushMatch(30);
    store.pushMismatch(1);
    store.pushMatch(39);
    store.pushMismatch(1);
    store.pushMatch(29);
    store.pushInsertion(3);
    store.pushMatch(97);

    // Segment 2: identical, 300 bases
    store.pushMatch(300);

    // forEachOp should coalesce the match at end of seg0 with match at start of seg1,
    // and the match at end of seg1 with match at start of seg2.
    auto ops = collectOps(store, off, store.tokensSince(off));

    // Expected after coalescing:
    // match(500+30=530), mismatch(1), match(39), mismatch(1), match(29),
    // ins(3), match(97+300=397)
    REQUIRE(ops.size() == 7);
    REQUIRE(ops[0] == std::make_pair(uint8_t(0), uint32_t(530)));
    REQUIRE(ops[1] == std::make_pair(uint8_t(1), uint32_t(1)));
    REQUIRE(ops[2] == std::make_pair(uint8_t(0), uint32_t(39)));
    REQUIRE(ops[3] == std::make_pair(uint8_t(1), uint32_t(1)));
    REQUIRE(ops[4] == std::make_pair(uint8_t(0), uint32_t(29)));
    REQUIRE(ops[5] == std::make_pair(uint8_t(2), uint32_t(3)));
    REQUIRE(ops[6] == std::make_pair(uint8_t(0), uint32_t(397)));

    // Total: read0 = 530+1+39+1+29+397 = 997, read1 = 997+3 = 1000
    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 997);
    REQUIRE(c1 == 1000);
}

TEST_CASE("Multi-segment: deletion segment between two match segments", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    // Segment 0: identical, 200 bases
    store.pushMatch(200);

    // Segment 1: empty on read1 side (len0=5, len1=0) → deletion
    store.pushDeletion(5);

    // Segment 2: identical, 200 bases
    store.pushMatch(200);

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 3);
    REQUIRE(ops[0] == std::make_pair(uint8_t(0), uint32_t(200)));
    REQUIRE(ops[1] == std::make_pair(uint8_t(3), uint32_t(5)));
    REQUIRE(ops[2] == std::make_pair(uint8_t(0), uint32_t(200)));

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 405);
    REQUIRE(c1 == 400);
}

// ============================================================================
// Long run splitting
// ============================================================================

TEST_CASE("Long match run exceeding MAX_LEN is split and coalesced", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    uint32_t longLen = 20000; // > 16383
    store.pushMatch(longLen);

    // Raw tokens should be split.
    auto tokens = store.getTokens(off, store.tokensSince(off));
    REQUIRE(tokens.size() == 2);
    REQUIRE(tokens[0].op() == 0);
    REQUIRE(tokens[0].len() == CigarToken::MAX_LEN);
    REQUIRE(tokens[1].op() == 0);
    REQUIRE(tokens[1].len() == longLen - CigarToken::MAX_LEN);

    // forEachOp should coalesce back.
    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0] == std::make_pair(uint8_t(0), uint32_t(20000)));

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 20000);
    REQUIRE(c1 == 20000);
}

TEST_CASE("Very long run: 50000 bases", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    store.pushMatch(50000);

    // Should need ceil(50000/16383) = 4 tokens.
    auto tokens = store.getTokens(off, store.tokensSince(off));
    REQUIRE(tokens.size() == 4);

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0].second == 50000);
}

// ============================================================================
// Thread-local merge
// ============================================================================

TEST_CASE("Merge two thread-local stores", "[CigarStore]") {
    OverlapCigarStore store0, store1;

    // Thread 0: one alignment with match(100)
    uint32_t off0 = store0.beginAlignment();
    store0.pushMatch(100);

    // Thread 1: one alignment with match(50), ins(3), match(47)
    uint32_t off1 = store1.beginAlignment();
    store1.pushMatch(50);
    store1.pushInsertion(3);
    store1.pushMatch(47);

    // Merge into global store.
    OverlapCigarStore global;
    uint32_t base0 = global.merge(store0);
    uint32_t base1 = global.merge(store1);

    REQUIRE(base0 == 0);
    REQUIRE(base1 == 1);
    
    // Verify alignment 0 (from thread 0).
    auto ops0 = collectOps(global, base0 + off0, store0.tokensSince(off0));
    REQUIRE(ops0.size() == 1);
    REQUIRE(ops0[0] == std::make_pair(uint8_t(0), uint32_t(100)));

    // Verify alignment 1 (from thread 1).
    auto ops1 = collectOps(global, base1 + off1, store1.tokensSince(off1));
    REQUIRE(ops1.size() == 3);
    REQUIRE(ops1[0] == std::make_pair(uint8_t(0), uint32_t(50)));
    REQUIRE(ops1[1] == std::make_pair(uint8_t(2), uint32_t(3)));
    REQUIRE(ops1[2] == std::make_pair(uint8_t(0), uint32_t(47)));
}

TEST_CASE("Merge multiple threads with multiple alignments each", "[CigarStore]") {
    OverlapCigarStore t0, t1;

    // Thread 0: 2 alignments
    uint32_t t0_a0 = t0.beginAlignment();
    t0.pushMatch(10);
    uint32_t t0_a1 = t0.beginAlignment();
    t0.pushMismatch(5);

    // Thread 1: 3 alignments
    uint32_t t1_a0 = t1.beginAlignment();
    t1.pushDeletion(7);
    uint32_t t1_a1 = t1.beginAlignment();
    t1.pushInsertion(3);
    uint32_t t1_a2 = t1.beginAlignment();
    t1.pushMatch(20);

    OverlapCigarStore global;
    uint32_t base0 = global.merge(t0);
    uint32_t base1 = global.merge(t1);

        REQUIRE(base0 == 0);
    REQUIRE(base1 == 2);

    // Spot-check: thread1's third alignment (global id = 2 + 2 = 4)
    auto ops = collectOps(global, base1 + t1_a2, t1.tokensSince(t1_a2));
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0] == std::make_pair(uint8_t(0), uint32_t(20)));
}

// ============================================================================
// Coalescing across segment boundaries
// ============================================================================

TEST_CASE("Coalescing: consecutive match segments merge", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    // 5 consecutive identical segments of 100 bases each.
    for(int i = 0; i < 5; i++) {
        store.pushMatch(100);
    }

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0] == std::make_pair(uint8_t(0), uint32_t(500)));
}

TEST_CASE("Coalescing: consecutive deletions from adjacent empty segments", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    // Two adjacent segments where read1 is empty.
    store.pushDeletion(10);
    store.pushDeletion(15);

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 1);
    REQUIRE(ops[0] == std::make_pair(uint8_t(3), uint32_t(25)));
}

TEST_CASE("No coalescing across different ops", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    store.pushMatch(100);
    store.pushMismatch(1);
    store.pushMatch(100);

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 3);
}

// ============================================================================
// Edge cases
// ============================================================================

TEST_CASE("Single-base mismatch", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMismatch(1);

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 1);
    REQUIRE(c1 == 1);
}

TEST_CASE("Alternating match/mismatch (worst case for token count)", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    // 100 bases, alternating match/mismatch.
    for(int i = 0; i < 100; i++) {
        if(i % 2 == 0) store.pushMatch(1);
        else store.pushMismatch(1);
    }

    // No coalescing possible — all ops alternate.
    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 100);

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 100);
    REQUIRE(c1 == 100);
}

TEST_CASE("Empty alignment (no segments)", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    auto tokens = store.getTokens(off, store.tokensSince(off));
    REQUIRE(tokens.size() == 0);

    auto ops = collectOps(store, off, store.tokensSince(off));
    REQUIRE(ops.size() == 0);

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == 0);
    REQUIRE(c1 == 0);
}

TEST_CASE("getTokens with invalid offset returns empty", "[CigarStore]") {
    OverlapCigarStore store;
    auto tokens = store.getTokens(uint32_t(-1), uint32_t(-1));
    REQUIRE(tokens.size() == 0);
}

// ============================================================================
// forEachOpWithPositions
// ============================================================================

TEST_CASE("forEachOpWithPositions: simple match", "[CigarStore][Positions]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(100);

    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> ops;
    store.forEachOpWithPositions(off, store.tokensSince(off), 1000, 2000, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        ops.emplace_back(op, len, r0, r1);
    });

    REQUIRE(ops.size() == 1);
    REQUIRE(std::get<0>(ops[0]) == 0);   // match
    REQUIRE(std::get<1>(ops[0]) == 100); // len
    REQUIRE(std::get<2>(ops[0]) == 1000); // read0 start
    REQUIRE(std::get<3>(ops[0]) == 2000); // read1 start
}

TEST_CASE("forEachOpWithPositions: match + ins + match", "[CigarStore][Positions]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushInsertion(3);
    store.pushMatch(50);

    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> ops;
    store.forEachOpWithPositions(off, store.tokensSince(off), 100, 200, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        ops.emplace_back(op, len, r0, r1);
    });

    REQUIRE(ops.size() == 3);
    // match(50) at r0=100, r1=200
    REQUIRE(std::get<0>(ops[0]) == 0);
    REQUIRE(std::get<1>(ops[0]) == 50);
    REQUIRE(std::get<2>(ops[0]) == 100);
    REQUIRE(std::get<3>(ops[0]) == 200);

    // ins(3) at r0=150, r1=250 (read0 didn't advance)
    REQUIRE(std::get<0>(ops[1]) == 2);
    REQUIRE(std::get<1>(ops[1]) == 3);
    REQUIRE(std::get<2>(ops[1]) == 150);
    REQUIRE(std::get<3>(ops[1]) == 250);

    // match(50) at r0=150, r1=253
    REQUIRE(std::get<0>(ops[2]) == 0);
    REQUIRE(std::get<1>(ops[2]) == 50);
    REQUIRE(std::get<2>(ops[2]) == 150);
    REQUIRE(std::get<3>(ops[2]) == 253);
}

TEST_CASE("forEachOpWithPositions: match + del + match", "[CigarStore][Positions]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(40);
    store.pushDeletion(5);
    store.pushMatch(40);

    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> ops;
    store.forEachOpWithPositions(off, store.tokensSince(off), 0, 0, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        ops.emplace_back(op, len, r0, r1);
    });

    REQUIRE(ops.size() == 3);
    // match(40) at r0=0, r1=0
    REQUIRE(std::get<2>(ops[0]) == 0);
    REQUIRE(std::get<3>(ops[0]) == 0);

    // del(5) at r0=40, r1=40 (read1 doesn't advance)
    REQUIRE(std::get<0>(ops[1]) == 3);
    REQUIRE(std::get<1>(ops[1]) == 5);
    REQUIRE(std::get<2>(ops[1]) == 40);
    REQUIRE(std::get<3>(ops[1]) == 40);

    // match(40) at r0=45, r1=40
    REQUIRE(std::get<0>(ops[2]) == 0);
    REQUIRE(std::get<1>(ops[2]) == 40);
    REQUIRE(std::get<2>(ops[2]) == 45);
    REQUIRE(std::get<3>(ops[2]) == 40);
}

// ============================================================================
// queryToTarget / targetToQuery
// ============================================================================

TEST_CASE("queryToTarget: positions in match region", "[CigarStore][CoordMap]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(100);

    // read0 starts at 1000, read1 starts at 5000.
    REQUIRE(store.queryToTarget(off, store.tokensSince(off), 1000, 5000, 1000) == 5000);
    REQUIRE(store.queryToTarget(off, store.tokensSince(off), 1000, 5000, 1050) == 5050);
    REQUIRE(store.queryToTarget(off, store.tokensSince(off), 1000, 5000, 1099) == 5099);
}

TEST_CASE("queryToTarget: position in deletion returns -1", "[CigarStore][CoordMap]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushDeletion(10);
    store.pushMatch(50);

    // Position 55 is inside the deletion (read0 positions 50..59).
    REQUIRE(store.queryToTarget(off, store.tokensSince(off), 0, 0, 55) == uint64_t(-1));
    // Position 60 is in the second match block.
    REQUIRE(store.queryToTarget(off, store.tokensSince(off), 0, 0, 60) == 50);
    // Position 70 maps to 60.
    REQUIRE(store.queryToTarget(off, store.tokensSince(off), 0, 0, 70) == 60);
}

TEST_CASE("queryToTarget: insertion shifts target positions", "[CigarStore][CoordMap]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushInsertion(5);
    store.pushMatch(50);

    // Before insertion: 1:1 mapping.
    REQUIRE(store.queryToTarget(off, store.tokensSince(off), 0, 0, 25) == 25);
    // After insertion: target is shifted by 5.
    REQUIRE(store.queryToTarget(off, store.tokensSince(off), 0, 0, 50) == 55);
    REQUIRE(store.queryToTarget(off, store.tokensSince(off), 0, 0, 75) == 80);
}

TEST_CASE("targetToQuery: positions in match region", "[CigarStore][CoordMap]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(100);

    REQUIRE(store.targetToQuery(off, store.tokensSince(off), 0, 0, 0) == 0);
    REQUIRE(store.targetToQuery(off, store.tokensSince(off), 0, 0, 50) == 50);
    REQUIRE(store.targetToQuery(off, store.tokensSince(off), 0, 0, 99) == 99);
}

TEST_CASE("targetToQuery: position in insertion returns -1", "[CigarStore][CoordMap]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushInsertion(10);
    store.pushMatch(50);

    // Position 55 is inside the insertion (read1 positions 50..59).
    REQUIRE(store.targetToQuery(off, store.tokensSince(off), 0, 0, 55) == uint64_t(-1));
    // Position 60 is in the second match block.
    REQUIRE(store.targetToQuery(off, store.tokensSince(off), 0, 0, 60) == 50);
}

TEST_CASE("targetToQuery: deletion shifts query positions", "[CigarStore][CoordMap]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushDeletion(5);
    store.pushMatch(50);

    REQUIRE(store.targetToQuery(off, store.tokensSince(off), 0, 0, 25) == 25);
    // After deletion: query is shifted by 5.
    REQUIRE(store.targetToQuery(off, store.tokensSince(off), 0, 0, 50) == 55);
    REQUIRE(store.targetToQuery(off, store.tokensSince(off), 0, 0, 75) == 80);
}

TEST_CASE("queryToTarget and targetToQuery are inverses", "[CigarStore][CoordMap]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(100);
    store.pushMismatch(5);
    store.pushMatch(200);
    store.pushInsertion(3);
    store.pushMatch(100);
    store.pushDeletion(7);
    store.pushMatch(100);

    // For every query position in a match/mismatch region,
    // targetToQuery(queryToTarget(pos)) should return pos.
    uint64_t r0 = 500, r1 = 1000;
    for(uint64_t qpos = r0; qpos < r0 + 100 + 5 + 200; qpos++) {
        uint64_t tpos = store.queryToTarget(off, store.tokensSince(off), r0, r1, qpos);
        REQUIRE(tpos != uint64_t(-1));
        uint64_t back = store.targetToQuery(off, store.tokensSince(off), r0, r1, tpos);
        REQUIRE(back == qpos);
    }
}

// ============================================================================
// walkRange
// ============================================================================

TEST_CASE("walkRange: full range equals forEachOpWithPositions", "[CigarStore][WalkRange]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushMismatch(2);
    store.pushMatch(48);

    // Walk the full range [0, 100).
    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> rangeOps;
    store.walkRange(off, store.tokensSince(off), 0, 0, 0, 100, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        rangeOps.emplace_back(op, len, r0, r1);
    });

    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> fullOps;
    store.forEachOpWithPositions(off, store.tokensSince(off), 0, 0, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        fullOps.emplace_back(op, len, r0, r1);
    });

    REQUIRE(rangeOps == fullOps);
}

TEST_CASE("walkRange: sub-range clips ops", "[CigarStore][WalkRange]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(100);
    store.pushMismatch(10);
    store.pushMatch(100);

    // Walk [90, 115) — clips into the first match, all of the mismatch,
    // and the start of the second match.
    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> ops;
    store.walkRange(off, store.tokensSince(off), 0, 0, 90, 115, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        ops.emplace_back(op, len, r0, r1);
    });

    REQUIRE(ops.size() == 3);
    // Clipped match: [90, 100), len=10
    REQUIRE(std::get<0>(ops[0]) == 0);
    REQUIRE(std::get<1>(ops[0]) == 10);
    REQUIRE(std::get<2>(ops[0]) == 90);
    REQUIRE(std::get<3>(ops[0]) == 90);

    // Full mismatch: [100, 110), len=10
    REQUIRE(std::get<0>(ops[1]) == 1);
    REQUIRE(std::get<1>(ops[1]) == 10);
    REQUIRE(std::get<2>(ops[1]) == 100);
    REQUIRE(std::get<3>(ops[1]) == 100);

    // Clipped match: [110, 115), len=5
    REQUIRE(std::get<0>(ops[2]) == 0);
    REQUIRE(std::get<1>(ops[2]) == 5);
    REQUIRE(std::get<2>(ops[2]) == 110);
    REQUIRE(std::get<3>(ops[2]) == 110);
}

TEST_CASE("walkRange: insertion within range is reported", "[CigarStore][WalkRange]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushInsertion(5);
    store.pushMatch(50);

    // Walk [40, 60) — should see tail of first match, the insertion, and
    // start of second match.
    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> ops;
    store.walkRange(off, store.tokensSince(off), 0, 0, 40, 60, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        ops.emplace_back(op, len, r0, r1);
    });

    REQUIRE(ops.size() == 3);
    // match [40, 50), len=10
    REQUIRE(std::get<0>(ops[0]) == 0);
    REQUIRE(std::get<1>(ops[0]) == 10);
    REQUIRE(std::get<2>(ops[0]) == 40);

    // insertion at r0=50, r1=50, len=5
    REQUIRE(std::get<0>(ops[1]) == 2);
    REQUIRE(std::get<1>(ops[1]) == 5);
    REQUIRE(std::get<2>(ops[1]) == 50);
    REQUIRE(std::get<3>(ops[1]) == 50);

    // match [50, 60), len=10, r1=55 (shifted by insertion)
    REQUIRE(std::get<0>(ops[2]) == 0);
    REQUIRE(std::get<1>(ops[2]) == 10);
    REQUIRE(std::get<2>(ops[2]) == 50);
    REQUIRE(std::get<3>(ops[2]) == 55);
}

TEST_CASE("walkRange: deletion within range", "[CigarStore][WalkRange]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushDeletion(5);
    store.pushMatch(50);

    // Walk [45, 60) — tail of first match, deletion, start of second match.
    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> ops;
    store.walkRange(off, store.tokensSince(off), 0, 0, 45, 60, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        ops.emplace_back(op, len, r0, r1);
    });

    REQUIRE(ops.size() == 3);
    // match [45, 50), len=5
    REQUIRE(std::get<0>(ops[0]) == 0);
    REQUIRE(std::get<1>(ops[0]) == 5);
    REQUIRE(std::get<2>(ops[0]) == 45);
    REQUIRE(std::get<3>(ops[0]) == 45);

    // deletion [50, 55), len=5, r1 stays at 50
    REQUIRE(std::get<0>(ops[1]) == 3);
    REQUIRE(std::get<1>(ops[1]) == 5);
    REQUIRE(std::get<2>(ops[1]) == 50);
    REQUIRE(std::get<3>(ops[1]) == 50);

    // match [55, 60), len=5, r1=50
    REQUIRE(std::get<0>(ops[2]) == 0);
    REQUIRE(std::get<1>(ops[2]) == 5);
    REQUIRE(std::get<2>(ops[2]) == 55);
    REQUIRE(std::get<3>(ops[2]) == 50);
}

// ============================================================================
// OverlapCigarStore::Cursor — resumable walking
// ============================================================================

// Helper: collect ops from walkRangeWithCursor.
static std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> collectCursorOps(
    OverlapCigarStore& store, OverlapCigarStore::Cursor& cur, uint64_t qs, uint64_t qe)
{
    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> ops;
    store.walkRangeWithCursor(cur, qs, qe, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        ops.emplace_back(op, len, r0, r1);
    });
    return ops;
}

// Helper: collect ops from walkRange (non-cursor).
static std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> collectWalkRangeOps(
    const OverlapCigarStore& store, uint32_t off, uint32_t cnt, uint64_t r0s, uint64_t r1s, uint64_t qs, uint64_t qe)
{
    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> ops;
    store.walkRange(off, cnt, r0s, r1s, qs, qe, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        ops.emplace_back(op, len, r0, r1);
    });
    return ops;
}

TEST_CASE("OverlapCigarStore::Cursor: matches walkRange for full range", "[CigarStore][Cursor]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushMismatch(2);
    store.pushInsertion(3);
    store.pushMatch(48);
    store.pushDeletion(5);
    store.pushMatch(50);

    OverlapCigarStore::Cursor cur;
    cur.reset(off, store.tokensSince(off), 0, 0, store);

    auto cursorOps = collectCursorOps(store, cur, 0, 155);
    auto walkOps = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, 0, 155);
    REQUIRE(cursorOps == walkOps);
}

TEST_CASE("OverlapCigarStore::Cursor: sliding window forward", "[CigarStore][Cursor]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    // 500 bases: match(200), mismatch(5), match(100), ins(3), match(192)
    store.pushMatch(200);
    store.pushMismatch(5);
    store.pushMatch(100);
    store.pushInsertion(3);
    store.pushMatch(192);

    OverlapCigarStore::Cursor cur;
    cur.reset(off, store.tokensSince(off), 0, 0, store);

    // Slide a 100-base window across the CIGAR.
    for(uint64_t s = 0; s < 497; s += 50) {
        uint64_t e = s + 100;
        if(e > 497) e = 497;

        auto cursorOps = collectCursorOps(store, cur, s, e);
        auto walkOps = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, s, e);
        REQUIRE(cursorOps == walkOps);
    }
}

TEST_CASE("OverlapCigarStore::Cursor: backward seek for overlapping windows", "[CigarStore][Cursor]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(100);
    store.pushMismatch(10);
    store.pushMatch(100);

    OverlapCigarStore::Cursor cur;
    cur.reset(off, store.tokensSince(off), 0, 0, store);

    // Walk [50, 150) first.
    auto ops1 = collectCursorOps(store, cur, 50, 150);
    auto expected1 = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, 50, 150);
    REQUIRE(ops1 == expected1);

    // Now walk [30, 130) — overlaps and starts before the cursor.
    auto ops2 = collectCursorOps(store, cur, 30, 130);
    auto expected2 = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, 30, 130);
    REQUIRE(ops2 == expected2);
}

TEST_CASE("OverlapCigarStore::Cursor: with non-zero start positions", "[CigarStore][Cursor]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(100);
    store.pushDeletion(5);
    store.pushMatch(100);

    OverlapCigarStore::Cursor cur;
    cur.reset(off, store.tokensSince(off), 1000, 5000, store);

    // Walk [1050, 1150).
    auto cursorOps = collectCursorOps(store, cur, 1050, 1150);
    auto walkOps = collectWalkRangeOps(store, off, store.tokensSince(off), 1000, 5000, 1050, 1150);
    REQUIRE(cursorOps == walkOps);

    // Walk [1100, 1200).
    auto cursorOps2 = collectCursorOps(store, cur, 1100, 1200);
    auto walkOps2 = collectWalkRangeOps(store, off, store.tokensSince(off), 1000, 5000, 1100, 1200);
    REQUIRE(cursorOps2 == walkOps2);
}

TEST_CASE("OverlapCigarStore::Cursor: sliding window across insertion", "[CigarStore][Cursor]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushInsertion(10);
    store.pushMatch(50);

    OverlapCigarStore::Cursor cur;
    cur.reset(off, store.tokensSince(off), 0, 0, store);

    // Window [0, 40) — before insertion.
    auto ops1 = collectCursorOps(store, cur, 0, 40);
    auto exp1 = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, 0, 40);
    REQUIRE(ops1 == exp1);

    // Window [40, 80) — spans the insertion point.
    auto ops2 = collectCursorOps(store, cur, 40, 80);
    auto exp2 = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, 40, 80);
    REQUIRE(ops2 == exp2);

    // Window [60, 100) — after insertion.
    auto ops3 = collectCursorOps(store, cur, 60, 100);
    auto exp3 = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, 60, 100);
    REQUIRE(ops3 == exp3);
}

TEST_CASE("OverlapCigarStore::Cursor: sliding window across deletion", "[CigarStore][Cursor]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(50);
    store.pushDeletion(10);
    store.pushMatch(50);

    OverlapCigarStore::Cursor cur;
    cur.reset(off, store.tokensSince(off), 0, 0, store);

    // Window [0, 40).
    auto ops1 = collectCursorOps(store, cur, 0, 40);
    auto exp1 = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, 0, 40);
    REQUIRE(ops1 == exp1);

    // Window [40, 80) — spans the deletion.
    auto ops2 = collectCursorOps(store, cur, 40, 80);
    auto exp2 = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, 40, 80);
    REQUIRE(ops2 == exp2);

    // Window [70, 110) — after deletion.
    auto ops3 = collectCursorOps(store, cur, 70, 110);
    auto exp3 = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, 70, 110);
    REQUIRE(ops3 == exp3);
}

TEST_CASE("OverlapCigarStore::Cursor: many small windows (phasing-like pattern)", "[CigarStore][Cursor]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    // Realistic-ish CIGAR: ~5000 bases with scattered variants.
    store.pushMatch(800);
    store.pushMismatch(1);
    store.pushMatch(400);
    store.pushInsertion(2);
    store.pushMatch(600);
    store.pushDeletion(3);
    store.pushMatch(500);
    store.pushMismatch(1);
    store.pushMatch(700);
    store.pushInsertion(1);
    store.pushMatch(300);
    store.pushMismatch(1);
    store.pushMatch(690);
    // Total read0: 800+1+400+600+3+500+1+700+300+1+690 = 3996

    OverlapCigarStore::Cursor cur;
    cur.reset(off, store.tokensSince(off), 0, 0, store);

    // Slide 200-base windows with 100-base step.
    for(uint64_t s = 0; s < 3996; s += 100) {
        uint64_t e = s + 200;
        if(e > 3996) e = 3996;
        if(s >= e) break;

        auto cursorOps = collectCursorOps(store, cur, s, e);
        auto walkOps = collectWalkRangeOps(store, off, store.tokensSince(off), 0, 0, s, e);
        REQUIRE(cursorOps == walkOps);
    }
}

TEST_CASE("walkRange: range entirely within one op", "[CigarStore][WalkRange]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();
    store.pushMatch(1000);

    std::vector<std::tuple<uint8_t, uint32_t, uint64_t, uint64_t>> ops;
    store.walkRange(off, store.tokensSince(off), 0, 0, 400, 600, [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        ops.emplace_back(op, len, r0, r1);
    });

    REQUIRE(ops.size() == 1);
    REQUIRE(std::get<0>(ops[0]) == 0);
    REQUIRE(std::get<1>(ops[0]) == 200);
    REQUIRE(std::get<2>(ops[0]) == 400);
    REQUIRE(std::get<3>(ops[0]) == 400);
}

// ============================================================================
// Realistic multi-segment overlap
// ============================================================================

TEST_CASE("Realistic overlap: 10 segments with mixed cases", "[CigarStore]") {
    OverlapCigarStore store;
    uint32_t off = store.beginAlignment();

    // Simulate a realistic overlap with 10 segments between 11 aligned markers.
    // Track expected consumed bases manually.
    uint64_t expected0 = 0, expected1 = 0;

    // Seg 0: identical, 800 bases
    store.pushMatch(800);
    expected0 += 800; expected1 += 800;

    // Seg 1: 1 mismatch in 600 bases (same length)
    store.pushMatch(299);
    store.pushMismatch(1);
    store.pushMatch(300);
    expected0 += 600; expected1 += 600;

    // Seg 2: identical, 1200 bases
    store.pushMatch(1200);
    expected0 += 1200; expected1 += 1200;

    // Seg 3: 2bp insertion in 500 bases (len0=500, len1=502)
    store.pushMatch(250);
    store.pushInsertion(2);
    store.pushMatch(250);
    expected0 += 500; expected1 += 502;

    // Seg 4: identical, 900 bases
    store.pushMatch(900);
    expected0 += 900; expected1 += 900;

    // Seg 5: 3bp deletion (len0=403, len1=400)
    store.pushMatch(200);
    store.pushDeletion(3);
    store.pushMatch(200);
    expected0 += 403; expected1 += 400;

    // Seg 6: empty segment (len0=0, len1=0) — nothing emitted

    // Seg 7: identical, 700 bases
    store.pushMatch(700);
    expected0 += 700; expected1 += 700;

    // Seg 8: len0=3, len1=0 — pure deletion
    store.pushDeletion(3);
    expected0 += 3; expected1 += 0;

    // Seg 9: identical, 500 bases
    store.pushMatch(500);
    expected0 += 500; expected1 += 500;

    auto [c0, c1] = consumedBases(store, off, store.tokensSince(off));
    REQUIRE(c0 == expected0);
    REQUIRE(c1 == expected1);

    // Verify total makes sense.
    REQUIRE(expected0 == 5606);
    REQUIRE(expected1 == 5602);
}
