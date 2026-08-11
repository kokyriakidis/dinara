// Unit tests for the myloasm chaining DP port (src/MyloasmChainKernel.hpp),
// exercised directly on hand-built synthetic anchor sets, with no Assembler
// dependency. Values are chosen so the DP arithmetic can be checked by hand
// (see the accompanying implementation plan for the worked-out expected scores).

#include "../external/catch2/catch.hpp"
#include "../src/MyloasmChainKernel.hpp"

using namespace dinara;

namespace {

// Default params match myloasm's own defaults (matchScore=11, gapCost=1,
// maxGap=200, doubleGap=10000, maxSkip=10, maxIter=50, minChainLength=3).
MyloasmChainParams defaultParams()
{
    return MyloasmChainParams{};
}

} // namespace

TEST_CASE("myloasmChainAnchors: same-strand collinear chain", "[myloasm]") {
    // Perfect diagonal: dist1 == dist2 == 10 at every step, strand 0 throughout.
    // By hand: f = [11, 21, 31, 41, 51] (each step adds min(10,10,matchScore=11)=10).
    MyloasmAnchorInput anchors;
    anchors.selfOffset = {0, 10, 20, 30, 40};
    anchors.offset     = {100, 110, 120, 130, 140};
    anchors.strand     = {0, 0, 0, 0, 0};

    const auto chains = myloasmChainAnchors(anchors, defaultParams());

    REQUIRE(chains.size() == 1);
    CHECK(chains[0].score == 51);
    CHECK(chains[0].originalIndices == std::vector<int64_t>{0, 1, 2, 3, 4});
}

TEST_CASE("myloasmChainAnchors: reverse-strand chain uses dinara's offset field directly", "[myloasm]") {
    // Same relative dist1/dist2 shape as the forward-strand test (step 10/10),
    // but built the way Dinara actually encodes a genuine reverse-strand
    // overlap: as self_offset increases, the raw target start position
    // decreases, and HifiasmKmerHit::offset stores (readLenB - 1 - posB) for
    // strand=1 hits (AssemblerInvertedIndex.cpp "offDiff"). This is a real
    // reverse-strand hit layout, not a synthetic pos2 value -- if the "offset
    // field can be used directly, no ones'-complement needed" derivation in
    // MyloasmChainKernel.hpp were wrong, this chain would fail to form.
    constexpr uint32_t readLenB = 1000;
    const std::vector<uint32_t> posB = {200, 190, 180, 170}; // decreasing raw start position.

    MyloasmAnchorInput anchors;
    anchors.selfOffset = {0, 10, 20, 30};
    anchors.strand     = {1, 1, 1, 1};
    for(const uint32_t p : posB) {
        anchors.offset.push_back(readLenB - 1U - p);
    }
    REQUIRE(anchors.offset == std::vector<uint32_t>{799, 809, 819, 829});

    const auto chains = myloasmChainAnchors(anchors, defaultParams());

    REQUIRE(chains.size() == 1);
    CHECK(chains[0].score == 41);
    CHECK(chains[0].originalIndices == std::vector<int64_t>{0, 1, 2, 3});
}

TEST_CASE("myloasmChainAnchors: two far-apart groups yield two independent chains", "[myloasm]") {
    // Two diagonals far enough apart (in both self_offset and offset) that no
    // anchor in one group can ever chain to the other (gap imbalance exceeds
    // maxGap), so the greedy claim-based extraction must recover both as
    // separate chains rather than merging or dropping one.
    MyloasmAnchorInput anchors;
    anchors.selfOffset = {0, 10, 20,    5000, 5010, 5020};
    anchors.offset     = {100, 110, 120, 9100, 9110, 9120};
    anchors.strand     = {0, 0, 0, 0, 0, 0};

    const auto chains = myloasmChainAnchors(anchors, defaultParams());

    REQUIRE(chains.size() == 2);
    CHECK(chains[0].score == 31);
    CHECK(chains[0].originalIndices == std::vector<int64_t>{0, 1, 2});
    CHECK(chains[1].score == 31);
    CHECK(chains[1].originalIndices == std::vector<int64_t>{3, 4, 5});
}

TEST_CASE("myloasmChainAnchors: chain shorter than minChainLength is dropped", "[myloasm]") {
    // Only two collinear anchors, with the default minChainLength=3: the one
    // candidate endpoint that clears the score heuristic traces back to a
    // 2-anchor chain, which is below minChainLength and must be dropped.
    MyloasmAnchorInput anchors;
    anchors.selfOffset = {0, 10};
    anchors.offset     = {100, 110};
    anchors.strand     = {0, 0};

    const auto chains = myloasmChainAnchors(anchors, defaultParams());

    CHECK(chains.empty());
}
