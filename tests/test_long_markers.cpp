/**
 * @file test_long_markers.cpp
 * @brief Validation of the DINARA_LONG_MARKERS (capacity-128 Kmer, 256-bit
 *        KmerId) path.
 *
 * These tests only run when the translation unit is compiled with
 * -DDINARA_LONG_MARKERS (see the test_long_markers target in CMakeLists.txt).
 * They exercise the primitives that make long marker k-mers (k > 64) possible:
 *   - Uint256 shift/mask/compare/cast arithmetic
 *   - ShortBaseSequence<__uint128_t> set/get round-trip at k > 64
 *   - id(k) -> ShortBaseSequence(id, k) round-trip (the packed-KmerId path)
 *   - reverseComplement correctness at k > 64 (fast path vs base-by-base)
 *   - extractKmer from a LongBaseSequenceView at k > 64, including k-mers that
 *     straddle three 64-base storage blocks
 */

#include "../external/catch2/catch.hpp"

#ifdef DINARA_LONG_MARKERS

#include "../src/Uint256.hpp"
#include "../src/Kmer.hpp"
#include "../src/ShortBaseSequence.hpp"
#include "../src/extractKmer.hpp"
#include "../src/LongBaseSequence.hpp"
#include "../src/Base.hpp"

#include <random>
#include <vector>

using namespace dinara;
using std::vector;

TEST_CASE("Uint256 arithmetic", "[long_markers]") {
    const Uint256 one(uint64_t(1));

    SECTION("shift round-trip across 64-bit boundary") {
        REQUIRE(((one << 70) >> 70) == one);
        REQUIRE(((one << 63) << 1) == (one << 64));
    }
    SECTION("high-word placement") {
        const Uint256 c = one << 192;
        REQUIRE(c.w[3] == 1ULL);
        REQUIRE(c.w[0] == 0ULL);
        REQUIRE(c.w[1] == 0ULL);
        REQUIRE(c.w[2] == 0ULL);
    }
    SECTION("mask (1<<n)-1") {
        const Uint256 mask = (one << 130) - one;
        REQUIRE(mask.w[0] == ~0ULL);
        REQUIRE(mask.w[1] == ~0ULL);
        REQUIRE((mask.w[2] & 0x3ULL) == 0x3ULL);
        REQUIRE((mask.w[2] & ~0x3ULL) == 0ULL);
        REQUIRE(mask.w[3] == 0ULL);
    }
    SECTION("ordering across words") {
        REQUIRE((one << 65) > (one << 64));
        REQUIRE_FALSE((one << 64) > (one << 65));
    }
    SECTION("truncating cast to __uint128_t keeps low 128 bits") {
        const Uint256 e = (one << 100) | Uint256(uint64_t(5));
        const __uint128_t lo = __uint128_t(e);
        REQUIRE(uint64_t(lo) == 5ULL);
        REQUIRE(uint64_t(lo >> 64) == (1ULL << (100 - 64)));
    }
}

TEST_CASE("Kmer128 build shape", "[long_markers]") {
    REQUIRE(Kmer::capacity == 128);
    REQUIRE(sizeof(KmerId) == 32);
}

static Base randomBase(std::mt19937& rng) {
    return Base::fromInteger(uint8_t(rng() & 3));
}

TEST_CASE("ShortBaseSequence set/get and id round-trip at k=70", "[long_markers]") {
    std::mt19937 rng(12345);
    const uint64_t k = 70;
    for(int trial = 0; trial < 2000; ++trial) {
        Kmer kmer;
        vector<Base> bases(k);
        for(uint64_t i = 0; i < k; ++i) {
            const Base b = randomBase(rng);
            bases[i] = b;
            kmer.set(i, b);
        }
        for(uint64_t i = 0; i < k; ++i) {
            REQUIRE(kmer[i].value == bases[i].value);
        }
        const KmerId id = kmer.id(k);
        const Kmer reconstructed(id, k);
        for(uint64_t i = 0; i < k; ++i) {
            REQUIRE(reconstructed[i].value == bases[i].value);
        }
        REQUIRE(kmer.id(k) == id); // deterministic
    }
}

TEST_CASE("reverseComplement at k=70", "[long_markers]") {
    std::mt19937 rng(999);
    const uint64_t k = 70;
    for(int trial = 0; trial < 1000; ++trial) {
        Kmer kmer;
        vector<Base> bases(k);
        for(uint64_t i = 0; i < k; ++i) {
            const Base b = randomBase(rng);
            bases[i] = b;
            kmer.set(i, b);
        }
        const Kmer rc = kmer.reverseComplement(k);
        for(uint64_t i = 0; i < k; ++i) {
            REQUIRE(rc[i].value == bases[k - 1 - i].complement().value);
        }
        const Kmer rcrc = rc.reverseComplement(k);
        for(uint64_t i = 0; i < k; ++i) {
            REQUIRE(rcrc[i].value == bases[i].value);
        }
    }
}

TEST_CASE("extractKmer at k=70 across block boundaries", "[long_markers]") {
    std::mt19937 rng(2024);
    const uint64_t k = 70;
    const uint64_t readLen = 300;
    const uint64_t blocks = (readLen + 63) / 64;
    vector<uint64_t> storage(blocks * 2, 0);
    LongBaseSequenceView read(storage.data(), readLen);
    vector<Base> readBases(readLen);
    for(uint64_t i = 0; i < readLen; ++i) {
        const Base b = randomBase(rng);
        readBases[i] = b;
        read.set(i, b);
    }
    for(uint64_t i = 0; i < readLen; ++i) {
        REQUIRE(read[i].value == readBases[i].value);
    }
    // Step by 1 so we hit every intra-block offset, including the worst case
    // where a k=70 k-mer straddles three storage blocks.
    for(uint64_t pos = 0; pos + k <= readLen; ++pos) {
        Kmer kmer;
        extractKmer(read, pos, k, kmer);
        for(uint64_t i = 0; i < k; ++i) {
            REQUIRE(kmer[i].value == readBases[pos + i].value);
        }
    }
}

#endif // DINARA_LONG_MARKERS
