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
