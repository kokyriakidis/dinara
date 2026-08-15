/**
 * @file test_integration.cpp
 * @brief Full integration tests for marker generation and filtering using real Assembler.
 * 
 * These tests verify the complete dinara marker pipeline:
 * 1. Read loading from FASTQ
 * 2. Marker generation using SIMD closed syncmers
 * 3. K-mer counting from marker KmerIds
 * 4. Marker filtering by frequency
 */

#include "../external/catch2/catch.hpp"

// Dinara headers
#include "../src/Assembler.hpp"
#include "../src/Reads.hpp"
#include "../src/KmerCounter.hpp"
#include "../src/ProjectedAlignment.hpp"
#include "../src/DINARA_ASSERT.hpp"
#include "../src/markerAccessFunctions.hpp"
#include "../src/MarkerGraph.hpp"
#include "../src/mode3-DirectedAnchorGraph.hpp"

// Standard library
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <tuple>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

using namespace dinara;

namespace {
class ScopedSilenceIo {
public:
    ScopedSilenceIo()
        : nullStream("/dev/null")
        , oldCoutBuf(std::cout.rdbuf())
        , oldCerrBuf(std::cerr.rdbuf())
    {
        std::cout.rdbuf(nullStream.rdbuf());
        std::cerr.rdbuf(nullStream.rdbuf());
    }

    ~ScopedSilenceIo()
    {
        std::cout.rdbuf(oldCoutBuf);
        std::cerr.rdbuf(oldCerrBuf);
    }

    ScopedSilenceIo(const ScopedSilenceIo&) = delete;
    ScopedSilenceIo& operator=(const ScopedSilenceIo&) = delete;

private:
    std::ofstream nullStream;
    std::streambuf* oldCoutBuf;
    std::streambuf* oldCerrBuf;
};

template<class F>
decltype(auto) withSilencedIo(F&& f)
{
    ScopedSilenceIo silence;
    return std::forward<F>(f)();
}

class ScopedCurrentPath {
public:
    explicit ScopedCurrentPath(const fs::path& p) : old(fs::current_path())
    {
        fs::current_path(p);
    }
    ~ScopedCurrentPath() { fs::current_path(old); }

    ScopedCurrentPath(const ScopedCurrentPath&) = delete;
    ScopedCurrentPath& operator=(const ScopedCurrentPath&) = delete;

private:
    fs::path old;
};

template<class F>
decltype(auto) withSilencedIoInDir(const fs::path& dir, F&& f)
{
    ScopedCurrentPath scoped(dir);
    return withSilencedIo(std::forward<F>(f));
}

fs::path makeUniqueTempDir(const std::string& prefix)
{
    static std::atomic<uint64_t> counter{0};
    const auto base = fs::temp_directory_path();

    for (int attempt = 0; attempt < 100; ++attempt) {
        const auto suffix =
            std::to_string(::getpid()) + "_" +
            std::to_string(counter.fetch_add(1)) + "_" +
            std::to_string(uint64_t(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::path dir = base / (prefix + suffix);
        std::error_code ec;
        if (fs::create_directory(dir, ec)) {
            return dir;
        }
    }
    throw std::runtime_error("Failed to create unique temp directory under " + base.string());
}

std::string toString(LongBaseSequenceView view)
{
    std::string s;
    s.resize(view.baseCount);
    for (size_t i = 0; i < view.baseCount; ++i) {
        static const char map[] = {'A', 'C', 'G', 'T'};
        s[i] = map[view[i].value & 3U];
    }
    return s;
}

template<class AlignmentDataContainer>
uint32_t findAlignmentEvidenceId(
    const AlignmentDataContainer& alignmentData,
    ReadId a,
    ReadId b,
    std::optional<bool> expectedSameStrand = std::nullopt)
{
    for (const auto& ad : alignmentData) {
        if (ad.isDeleted()) continue;
        const bool matchesPair =
            (ad.readIds[0] == a && ad.readIds[1] == b) ||
            (ad.readIds[0] == b && ad.readIds[1] == a);
        if (!matchesPair) continue;
        if (expectedSameStrand && ad.isSameStrand != *expectedSameStrand) continue;
        return uint32_t(ad.info.alignmentId);
    }
    return invalid<uint32_t>;
}

template<class AlignmentDataContainer>
const AlignmentData* findAlignmentDataPtr(
    const AlignmentDataContainer& alignmentData,
    ReadId a,
    ReadId b,
    std::optional<bool> expectedSameStrand = std::nullopt)
{
    for (const auto& ad : alignmentData) {
        if (ad.isDeleted()) continue;
        const bool matchesPair =
            (ad.readIds[0] == a && ad.readIds[1] == b) ||
            (ad.readIds[0] == b && ad.readIds[1] == a);
        if (!matchesPair) continue;
        if (expectedSameStrand && ad.isSameStrand != *expectedSameStrand) continue;
        return &ad;
    }
    return nullptr;
}

AlignOptions makeDefaultAlignOptionsForIntegration()
{
    AlignOptions options;
    options.alignMethod = 6;
    options.maxSkip = 100;
    options.maxDrift = 100;
    options.maxTrim = 10000;
    options.minAlignedMarkerCount = 4;
    options.minAlignedFraction = 0.0;
    options.maxMarkerFrequency = 1000;
    options.matchScore = 3;
    options.mismatchScore = -1;
    options.gapScore = -1;
    options.downsamplingFactor = 0.1;
    options.bandExtend = 10;
    options.maxBand = 1000;
    options.sameChannelReadAlignmentSuppressDeltaThreshold = 0;
    options.suppressContainments = false;
    options.align4DeltaX = 200;
    options.align4DeltaY = 10;
    options.align4MinEntryCountPerCell = 10;
    options.align4MaxDistanceFromBoundary = 100;
    options.align5DriftRateTolerance = 0.02;
    options.align5MinBandExtend = 10;
    options.maxErrorRate = 0.3;
    // Base-level overlap DP scoring defaults (hifiasm/minimap2 parity).
    options.overlapDpMatchScore = 2;
    options.overlapDpMismatchScore = -4;
    options.overlapDpGapOpen1 = 4;
    options.overlapDpGapExtend1 = 2;
    options.overlapDpGapOpen2 = 24;
    options.overlapDpGapExtend2 = 1;
    return options;
}

} // namespace

TEST_CASE("Marker-vertex overlap splitting removes bridge reads and peels dense cores", "[anchors][vertexSplit]")
{
    using dinara::testing::splitVertexByOverlapSupportForTesting;

    // Build an undirected adjacency list with sorted unique neighbors.
    auto makeAdj = [](uint32_t n) {
        return std::vector<std::vector<uint32_t>>(n);
    };
    auto addEdge = [](std::vector<std::vector<uint32_t>>& adj, uint32_t u, uint32_t v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    };
    auto normalize = [](std::vector<std::vector<uint32_t>>& adj) {
        for(auto& nbr : adj) {
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
        }
    };

    const uint64_t minCov = 3;
    const uint64_t maxCov = 100;

    SECTION("Single bridge read (articulation) is dropped and vertex splits into two anchors")
    {
        // Clique1: 0..4, Clique2: 5..8, bridge: 9 connects to everyone.
        auto adj = makeAdj(10);
        for(uint32_t i=0; i<5; ++i) {
            for(uint32_t j=i+1; j<5; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=5; i<9; ++i) {
            for(uint32_t j=i+1; j<9; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=0; i<9; ++i) addEdge(adj, i, 9);
        normalize(adj);

        const auto groups = splitVertexByOverlapSupportForTesting(adj, {}, false, minCov, maxCov);
        REQUIRE(groups.size() == 2);
        std::vector<uint32_t> sizes;
        for(const auto& g : groups) sizes.push_back(uint32_t(g.size()));
        std::sort(sizes.begin(), sizes.end());
        CHECK(sizes[0] == 4);
        CHECK(sizes[1] == 5);
        // Bridge read 9 should not be present.
        for(const auto& g : groups) {
            CHECK(std::find(g.begin(), g.end(), 9) == g.end());
        }
    }

    SECTION("Two bridge reads (no articulation) are handled by greedy dense-core split + peeling")
    {
        // Clique1: 0..4, Clique2: 5..8, bridges: 9 and 10 connect to everyone and to each other.
        auto adj = makeAdj(11);
        for(uint32_t i=0; i<5; ++i) {
            for(uint32_t j=i+1; j<5; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=5; i<9; ++i) {
            for(uint32_t j=i+1; j<9; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=0; i<9; ++i) {
            addEdge(adj, i, 9);
            addEdge(adj, i, 10);
        }
        addEdge(adj, 9, 10);
        normalize(adj);

        const auto groups = splitVertexByOverlapSupportForTesting(adj, {}, false, minCov, maxCov);
        REQUIRE(groups.size() >= 2);

        // Expect at least the two cliques to be extracted as dense cores.
        // Bridge reads may be peeled away or attached to one group; ensure cliques are present.
        auto containsAll = [](const std::vector<uint32_t>& g, uint32_t lo, uint32_t hiExclusive) {
            for(uint32_t x=lo; x<hiExclusive; ++x) {
                if(std::find(g.begin(), g.end(), x) == g.end()) return false;
            }
            return true;
        };
        bool foundClique1 = false;
        bool foundClique2 = false;
        for(const auto& g : groups) {
            if(containsAll(g, 0, 5)) foundClique1 = true;
            if(containsAll(g, 5, 9)) foundClique2 = true;
        }
        CHECK(foundClique1);
        CHECK(foundClique2);
    }
}

TEST_CASE("Vertex splitting can use non-contained cores and drop contained bridges", "[anchors][vertexSplit][contained]")
{
    using dinara::testing::splitVertexByOverlapSupportWithCoreMaskForTesting;

    auto makeAdj = [](uint32_t n) {
        return std::vector<std::vector<uint32_t>>(n);
    };
    auto addEdge = [](std::vector<std::vector<uint32_t>>& adj, uint32_t u, uint32_t v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    };
    auto normalize = [](std::vector<std::vector<uint32_t>>& adj) {
        for(auto& nbr : adj) {
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
        }
    };

    const uint64_t minCov = 3;
    const uint64_t maxCov = 100;

    // Two equal cliques connected only via "contained" bridge nodes.
    // Core cliques: 0..4 and 5..9. Bridges (contained): 10 and 11 connect to everyone in both cliques.
    auto adj = makeAdj(12);
    for(uint32_t i=0; i<5; ++i) {
        for(uint32_t j=i+1; j<5; ++j) addEdge(adj, i, j);
    }
    for(uint32_t i=5; i<10; ++i) {
        for(uint32_t j=i+1; j<10; ++j) addEdge(adj, i, j);
    }
    for(uint32_t b=10; b<12; ++b) {
        for(uint32_t i=0; i<10; ++i) addEdge(adj, b, i);
    }
    normalize(adj);

    std::vector<uint8_t> isCore(12, 1);
    isCore[10] = 0;
    isCore[11] = 0;

    const auto groups = splitVertexByOverlapSupportWithCoreMaskForTesting(
        adj, {}, false, isCore,
        /*coreMinSize*/ 3,
        /*attachMinSupport*/ 1,
        minCov, maxCov);

    REQUIRE(groups.size() == 2);
    for(const auto& g : groups) {
        CHECK(std::find(g.begin(), g.end(), 10) == g.end());
        CHECK(std::find(g.begin(), g.end(), 11) == g.end());
    }
}

TEST_CASE("Vertex splitting with non-contained cores: additional edge cases", "[anchors][vertexSplit][contained][core]")
{
    using dinara::testing::splitVertexByOverlapSupportForTesting;
    using dinara::testing::splitVertexByOverlapSupportWithCoreMaskForTesting;

    auto makeAdj = [](uint32_t n) {
        return std::vector<std::vector<uint32_t>>(n);
    };
    auto addEdge = [](std::vector<std::vector<uint32_t>>& adj, uint32_t u, uint32_t v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    };
    auto normalize = [](std::vector<std::vector<uint32_t>>& adj) {
        for(auto& nbr : adj) {
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
        }
    };

    const uint64_t minCov = 3;
    const uint64_t maxCov = 100;

    SECTION("Small core (2 nodes) can still drive a split after attaching contained reads")
    {
        // Core singletons: 0 and 5 (no edge). Contained nodes 1..4 attach to 0; 6..9 attach to 5.
        auto adj = makeAdj(10);
        for(uint32_t u=1; u<=4; ++u) addEdge(adj, 0, u);
        for(uint32_t u=6; u<=9; ++u) addEdge(adj, 5, u);
        // Add internal support among contained nodes in each side (to avoid being isolated after attach).
        for(uint32_t i=1; i<=4; ++i) {
            for(uint32_t j=i+1; j<=4; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=6; i<=9; ++i) {
            for(uint32_t j=i+1; j<=9; ++j) addEdge(adj, i, j);
        }
        normalize(adj);

        std::vector<uint8_t> isCore(10, 0);
        isCore[0] = 1;
        isCore[5] = 1;

        const auto groups = splitVertexByOverlapSupportWithCoreMaskForTesting(
            adj, {}, false, isCore,
            /*coreMinSize*/ 2,
            /*attachMinSupport*/ 1,
            minCov, maxCov);

        REQUIRE(groups.size() == 2);
        std::vector<uint32_t> sizes;
        for(const auto& g : groups) sizes.push_back(uint32_t(g.size()));
        std::sort(sizes.begin(), sizes.end());
        CHECK(sizes[0] >= 5);
        CHECK(sizes[1] >= 5);
        // Each core should remain present.
        bool has0 = false;
        bool has5 = false;
        for(const auto& g : groups) {
            if(std::find(g.begin(), g.end(), 0) != g.end()) has0 = true;
            if(std::find(g.begin(), g.end(), 5) != g.end()) has5 = true;
        }
        CHECK(has0);
        CHECK(has5);
    }

    SECTION("Ambiguous contained bridge (ties) is dropped instead of gluing clusters")
    {
        // Core cliques: 0..2 and 3..5. Contained bridge: 6 connects equally to both cliques.
        auto adj = makeAdj(7);
        for(uint32_t i=0; i<3; ++i) {
            for(uint32_t j=i+1; j<3; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=3; i<6; ++i) {
            for(uint32_t j=i+1; j<6; ++j) addEdge(adj, i, j);
        }
        // Bridge edges: 6 connects to one node in each clique (tie support 1 vs 1).
        addEdge(adj, 6, 0);
        addEdge(adj, 6, 3);
        normalize(adj);

        std::vector<uint8_t> isCore(7, 1);
        isCore[6] = 0;

        const auto groups = splitVertexByOverlapSupportWithCoreMaskForTesting(
            adj, {}, false, isCore,
            /*coreMinSize*/ 2,
            /*attachMinSupport*/ 1,
            /*minAnchorCoverage*/ 2,
            maxCov);

        REQUIRE(groups.size() == 2);
        for(const auto& g : groups) {
            CHECK(std::find(g.begin(), g.end(), 6) == g.end());
        }
    }

    SECTION("When all reads are core (no contained), core-mode falls back to normal splitting")
    {
        // Two cliques connected by a single bridge node (articulation).
        auto adj = makeAdj(10);
        for(uint32_t i=0; i<5; ++i) {
            for(uint32_t j=i+1; j<5; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=5; i<9; ++i) {
            for(uint32_t j=i+1; j<9; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=0; i<9; ++i) addEdge(adj, i, 9);
        normalize(adj);

        // All core.
        std::vector<uint8_t> isCore(10, 1);

        const auto groupsCoreMode = splitVertexByOverlapSupportWithCoreMaskForTesting(
            adj, {}, false, isCore,
            /*coreMinSize*/ 2,
            /*attachMinSupport*/ 1,
            /*minAnchorCoverage*/ 3,
            maxCov);

        const auto groupsBaseline = splitVertexByOverlapSupportForTesting(
            adj, {}, false, /*minAnchorCoverage*/ 3, maxCov);

        // Both paths should agree on producing a split (2 groups) and dropping the bridge.
        REQUIRE(groupsCoreMode.size() == 2);
        REQUIRE(groupsBaseline.size() == 2);
        for(const auto& g : groupsCoreMode) {
            CHECK(std::find(g.begin(), g.end(), 9) == g.end());
        }
    }
}

TEST_CASE("Vertex splitting: dynamic attachment threshold blocks high-degree weak assignments", "[anchors][vertexSplit][contained][core][dynamic]")
{
    using dinara::testing::splitVertexByOverlapSupportWithCoreMaskForTesting;

    auto makeAdj = [](uint32_t n) {
        return std::vector<std::vector<uint32_t>>(n);
    };
    auto addEdge = [](std::vector<std::vector<uint32_t>>& adj, uint32_t u, uint32_t v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    };
    auto normalize = [](std::vector<std::vector<uint32_t>>& adj) {
        for(auto& nbr : adj) {
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
        }
    };

    // Core: nodes 0..11, with a small core cluster {0,1} and 10 isolated core nodes.
    // Non-core (contained): node 12 connects to both 0 and 1 (bestCount=2) and to each isolated core (secondCount=1).
    // degToCore=12 triggers dynamicMinSupport=3, so node 12 should NOT attach anywhere.
    auto adj = makeAdj(13);
    addEdge(adj, 0, 1);
    for(uint32_t u=2; u<=11; ++u) {
        addEdge(adj, 12, u);
    }
    addEdge(adj, 12, 0);
    addEdge(adj, 12, 1);
    normalize(adj);

    std::vector<uint8_t> isCore(13, 1);
    isCore[12] = 0;

    const auto groups = splitVertexByOverlapSupportWithCoreMaskForTesting(
        adj, {}, false, isCore,
        /*coreMinSize*/ 2,
        /*attachMinSupport*/ 1,
        /*minAnchorCoverage*/ 1,
        /*maxAnchorCoverage*/ 1000);

    // All core groups are kept at minCov=1; ensure node 12 is not present.
    REQUIRE(!groups.empty());
    for(const auto& g : groups) {
        CHECK(std::find(g.begin(), g.end(), 12) == g.end());
    }
}

TEST_CASE("Vertex splitting: core split fallback preserves baseline splits when core cannot split", "[anchors][vertexSplit][contained][core][fallback]")
{
    using dinara::testing::splitVertexByOverlapSupportForTesting;
    using dinara::testing::splitVertexByOverlapSupportWithCoreMaskForTesting;

    auto makeAdj = [](uint32_t n) {
        return std::vector<std::vector<uint32_t>>(n);
    };
    auto addEdge = [](std::vector<std::vector<uint32_t>>& adj, uint32_t u, uint32_t v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    };
    auto normalize = [](std::vector<std::vector<uint32_t>>& adj) {
        for(auto& nbr : adj) {
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
        }
    };

    // Core nodes 0..3 form one clique (so core split cannot produce 2+ groups).
    // Non-core nodes 4..5 attach to core; non-core nodes 6..8 form a separate clique disconnected from core.
    // Baseline splitter should return 2 components; core-mode should fall back and do the same.
    auto adj = makeAdj(9);
    for(uint32_t i=0; i<4; ++i) {
        for(uint32_t j=i+1; j<4; ++j) addEdge(adj, i, j);
    }
    addEdge(adj, 0, 4);
    addEdge(adj, 1, 5);
    addEdge(adj, 4, 5);
    for(uint32_t i=6; i<9; ++i) {
        for(uint32_t j=i+1; j<9; ++j) addEdge(adj, i, j);
    }
    normalize(adj);

    std::vector<uint8_t> isCore(9, 0);
    for(uint32_t i=0; i<4; ++i) isCore[i] = 1;

    const auto baseline = splitVertexByOverlapSupportForTesting(adj, {}, false, /*minCov*/ 2, /*max*/ 1000);
    const auto coreMode = splitVertexByOverlapSupportWithCoreMaskForTesting(
        adj, {}, false, isCore,
        /*coreMinSize*/ 2,
        /*attachMinSupport*/ 1,
        /*minCov*/ 2,
        /*max*/ 1000);

    REQUIRE(baseline.size() == 2);
    REQUIRE(coreMode.size() == 2);

    // Ensure the disconnected clique {6,7,8} is preserved in coreMode.
    auto containsAll = [](const std::vector<uint32_t>& g, std::initializer_list<uint32_t> nodes) {
        for(uint32_t x : nodes) {
            if(std::find(g.begin(), g.end(), x) == g.end()) return false;
        }
        return true;
    };
    bool found = false;
    for(const auto& g : coreMode) {
        if(containsAll(g, {6,7,8})) found = true;
    }
    CHECK(found);
}

TEST_CASE("Clique-cover vertex splitter: splits disjoint core cliques and drops ambiguous bridges", "[anchors][vertexSplit][cliqueCover]")
{
    using dinara::testing::splitVertexByCliqueCoverForTesting;

    auto makeAdj = [](uint32_t n) {
        return std::vector<std::vector<uint32_t>>(n);
    };
    auto addEdge = [](std::vector<std::vector<uint32_t>>& adj, uint32_t u, uint32_t v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    };
    auto normalize = [](std::vector<std::vector<uint32_t>>& adj) {
        for(auto& nbr : adj) {
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
        }
    };

    const uint64_t minCov = 3;
    const uint64_t maxCov = 100;

    SECTION("Contained bridge excluded from core -> two core cliques split")
    {
        // Core clique A: 0..4, core clique B: 5..8, contained bridge: 9 connects to all.
        auto adj = makeAdj(10);
        for(uint32_t i=0; i<5; ++i) {
            for(uint32_t j=i+1; j<5; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=5; i<9; ++i) {
            for(uint32_t j=i+1; j<9; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=0; i<9; ++i) addEdge(adj, i, 9);
        normalize(adj);

        std::vector<uint8_t> isCore(10, 1);
        isCore[9] = 0;

        const auto groups = splitVertexByCliqueCoverForTesting(adj, {}, false, isCore, /*attachMinSupport*/ 1, minCov, maxCov);
        REQUIRE(groups.size() == 2);
        for(const auto& g : groups) {
            CHECK(std::find(g.begin(), g.end(), 9) == g.end());
        }
    }

    SECTION("If a bridge is treated as core, clique-cover returns empty (forces fallback in pipeline)")
    {
        auto adj = makeAdj(10);
        for(uint32_t i=0; i<5; ++i) {
            for(uint32_t j=i+1; j<5; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=5; i<9; ++i) {
            for(uint32_t j=i+1; j<9; ++j) addEdge(adj, i, j);
        }
        for(uint32_t i=0; i<9; ++i) addEdge(adj, i, 9);
        normalize(adj);

        std::vector<uint8_t> isCore(10, 1); // bridge is core too

        const auto groups = splitVertexByCliqueCoverForTesting(adj, {}, false, isCore, /*attachMinSupport*/ 1, minCov, maxCov);
        CHECK(groups.empty());
    }
}

TEST_CASE("Vertex splitting: MCL branch is attempted when enabled and suspicious", "[anchors][vertexSplit][mcl][auto]")
{
    using dinara::testing::autoSplitVertexWithMclTriedFlagForTesting;

    auto addClique = [](std::vector<std::vector<uint32_t>>& adj, uint32_t lo, uint32_t hiExclusive) {
        for(uint32_t i=lo; i<hiExclusive; ++i) {
            for(uint32_t j=i+1; j<hiExclusive; ++j) {
                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }
    };
    auto normalize = [](std::vector<std::vector<uint32_t>>& adj) {
        for(auto& nbr : adj) {
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
        }
    };

    // Use a clique so the baseline splitter returns a single group.
    std::vector<std::vector<uint32_t>> adj(6);
    addClique(adj, 0, 6);
    normalize(adj);

    std::vector<uint8_t> isCore(6, 1);
    const uint64_t minCov = 2;
    const uint64_t maxCov = 1000;

    const auto [noMclGroups, noMclTried] = autoSplitVertexWithMclTriedFlagForTesting(
        adj, {}, false, isCore,
        /*useNonContainedCores*/ false,
        /*coreMinSize*/ 2,
        /*attachMinSupport*/ 1,
        /*useMclSecondary*/ false,
        /*mclMinVertexSize*/ 0,
        /*mclInflation*/ 2.2,
        /*mclMaxIterations*/ 80,
        /*suspiciousMaxDensity*/ 1.0,
        /*suspiciousMaxAverageClustering*/ 1.0,
        /*minCov*/ minCov,
        /*max*/ maxCov);
    REQUIRE(noMclGroups.size() == 1);
    REQUIRE(noMclTried == false);

    const auto [withMclGroups, withMclTried] = autoSplitVertexWithMclTriedFlagForTesting(
        adj, {}, false, isCore,
        /*useNonContainedCores*/ false,
        /*coreMinSize*/ 2,
        /*attachMinSupport*/ 1,
        /*useMclSecondary*/ true,
        /*mclMinVertexSize*/ 0,
        /*mclInflation*/ 2.2,
        /*mclMaxIterations*/ 80,
        /*suspiciousMaxDensity*/ 1.0,
        /*suspiciousMaxAverageClustering*/ 1.0,
        /*minCov*/ minCov,
        /*max*/ maxCov);
    REQUIRE(withMclGroups.size() == 1);
    REQUIRE(withMclTried == true);
}

TEST_CASE("MCL secondary clustering can separate weakly-connected communities", "[anchors][mcl]")
{
    using dinara::testing::mclClusterForTesting;

    auto makeAdj = [](uint32_t n) {
        return std::vector<std::vector<uint32_t>>(n);
    };
    auto addEdge = [](std::vector<std::vector<uint32_t>>& adj, uint32_t u, uint32_t v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    };
    auto normalize = [](std::vector<std::vector<uint32_t>>& adj) {
        for(auto& nbr : adj) {
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
        }
    };

    // Two cliques connected by two weak cross edges.
    // Clique1: 0..4, Clique2: 5..9, cross edges: (0,5) and (1,6).
    auto adj = makeAdj(10);
    for(uint32_t i=0; i<5; ++i) {
        for(uint32_t j=i+1; j<5; ++j) addEdge(adj, i, j);
    }
    for(uint32_t i=5; i<10; ++i) {
        for(uint32_t j=i+1; j<10; ++j) addEdge(adj, i, j);
    }
    addEdge(adj, 0, 5);
    addEdge(adj, 1, 6);
    normalize(adj);

    const auto groups = mclClusterForTesting(adj, /*inflation*/ 1.8, /*maxIterations*/ 50);
    REQUIRE(groups.size() >= 2);

    auto containsAll = [](const std::vector<uint32_t>& g, uint32_t lo, uint32_t hiExclusive) {
        for(uint32_t x=lo; x<hiExclusive; ++x) {
            if(std::find(g.begin(), g.end(), x) == g.end()) return false;
        }
        return true;
    };
    bool foundClique1 = false;
    bool foundClique2 = false;
    for(const auto& g : groups) {
        if(containsAll(g, 0, 5)) foundClique1 = true;
        if(containsAll(g, 5, 10)) foundClique2 = true;
    }
    CHECK(foundClique1);
    CHECK(foundClique2);
}

TEST_CASE("MarkerGraph vertex coverage histogram counts canonical vertices", "[markerGraph][histogram]")
{
    MarkerGraph markerGraph;
    markerGraph.constructVertices();
    markerGraph.vertices().createNew("", 4096);

    // Create 4 vertices with coverages: [3,3,1,2].
    markerGraph.vertices().appendVector(3);
    markerGraph.vertices().appendVector(3);
    markerGraph.vertices().appendVector(1);
    markerGraph.vertices().appendVector(2);

    // Reverse-complement pairs: (0,1) and (2,3).
    markerGraph.reverseComplementVertex.createNew("", 4096);
    markerGraph.reverseComplementVertex.resize(4);
    markerGraph.reverseComplementVertex[0] = 1;
    markerGraph.reverseComplementVertex[1] = 0;
    markerGraph.reverseComplementVertex[2] = 3;
    markerGraph.reverseComplementVertex[3] = 2;

    const auto all = markerGraph.computeVertexCoverageHistogram(false);
    REQUIRE(all.size() >= 4);
    CHECK(all[1] == 1);
    CHECK(all[2] == 1);
    CHECK(all[3] == 2);

    const auto canonical = markerGraph.computeVertexCoverageHistogram(true);
    REQUIRE(canonical.size() >= 4);
    CHECK(canonical[1] == 1); // vertex 2
    CHECK(canonical[2] == 0);
    CHECK(canonical[3] == 1); // vertex 0
}

// =============================================================================
// HELPER: Generate random DNA sequence
// =============================================================================
std::string randomSequence(size_t length, uint32_t seed = 42) {
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dis(0, 3);
    std::string seq(length, 'N');
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (size_t i = 0; i < length; i++) {
        seq[i] = bases[dis(gen)];
    }
    return seq;
}

// =============================================================================
// FIXTURE: Creates temp directory with test FASTQ and initializes Assembler
// =============================================================================
class AssemblerIntegrationFixture {
public:
    fs::path testDir;
    fs::path fastqPath;
    std::unique_ptr<Assembler> assembler;
    
    AssemblerIntegrationFixture() {
        testDir = makeUniqueTempDir("dinara_integration_");
        fastqPath = testDir / "test_reads.fastq";
    }
    
    ~AssemblerIntegrationFixture() {
        assembler.reset();
        try {
            fs::remove_all(testDir);
        } catch (...) {}
    }
    
    void createFastq(const std::vector<std::string>& sequences) {
        std::ofstream out(fastqPath);
        if (!out) {
            throw std::runtime_error("Failed to open " + fastqPath.string());
        }
        for (size_t i = 0; i < sequences.size(); i++) {
            out << "@read_" << i << "\n";
            out << sequences[i] << "\n";
            out << "+\n";
            out << std::string(sequences[i].size(), '~') << "\n";
        }
    }
    
    void initAssembler() {
        withSilencedIoInDir(testDir, [&] {
            assembler = std::make_unique<Assembler>(
                testDir.string() + "/",
                true,   // createNew
                0,      // readRepresentation = raw
                4096    // pageSize
            );
        });
    }
    
    void loadReads(uint64_t minReadLength = 0) {
        withSilencedIoInDir(testDir, [&] {
            assembler->addReads(fastqPath.string(), minReadLength, true, 1);
            assembler->computeReadIdsSortedByName();
        });
    }
    
    void generateMarkers(int k = 16, int s = 4) {
        // The simd-minimizers backend requires the window length l = k + w - 1
        // to be odd (k and w same parity). Nudge w up by one when needed so the
        // various (k, s) combinations used by these tests stay valid.
        const int w = ((k + s) % 2 == 0) ? s : s + 1;
        withSilencedIoInDir(testDir, [&] {
            assembler->findMarkersSimdMinimizers(1, k, w, /*useHifiasm*/ false);
        });
    }
    
    void countKmers() {
        withSilencedIoInDir(testDir, [&] { assembler->countKmersFromMarkerKmerIds(1); });
    }
    
    void applyFilter(uint64_t minFreq, uint64_t maxFreq) {
        withSilencedIoInDir(testDir, [&] { assembler->applyKmerCountFilter(minFreq, maxFreq, 1); });
    }
};

TEST_CASE("ReadGraph.filterSecondaryRequireNonRedundantOnBothReads removes symmetric-only overlaps",
    "[integration][readgraph][secondary][symmetry]")
{
    AssemblerIntegrationFixture fixture;

    const std::string r0 = randomSequence(3000, 9301);
    const std::string r1 = randomSequence(3000, 9302);
    fixture.createFastq({r0, r1});
    fixture.initAssembler();
    fixture.loadReads();

    fixture.assembler->alignmentData.createNew("", 4096);

    AlignmentInfo infoA;
    infoA.dpScore = 1000;
    infoA.mismatchCount = 0;
    infoA.gapCount = 0;
    infoA.gapEventCount = 0;
    AlignmentData adA(array<ReadId, 2>{ReadId(0), ReadId(1)}, true, infoA);
    adA.qs = 0;
    adA.qe = 1000;
    adA.ts = 500;
    adA.te = 1500;
    fixture.assembler->alignmentData.push_back(adA);

    AlignmentInfo infoB = infoA;
    infoB.dpScore = 900; // within 80% of best
    AlignmentData adB(array<ReadId, 2>{ReadId(0), ReadId(1)}, true, infoB);
    adB.qs = 1600;
    adB.qe = 2600;
    // Highly redundant on the partner read interval, but non-overlapping on r0.
    adB.ts = 520;
    adB.te = 1520;
    fixture.assembler->alignmentData.push_back(adB);

    fixture.assembler->computeAlignmentTableForTesting();

    // Legacy behavior (r0-only redundancy): both are kept.
    fixture.assembler->filterSecondaryAlignmentsPerReadPair(1, false);
    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted0());
    CHECK_FALSE(fixture.assembler->alignmentData[1].isDeleted0());

    // Symmetric behavior: second overlap is removed due to partner-read redundancy only.
    fixture.assembler->filterSecondaryAlignmentsPerReadPair(1, true);
    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted0());
    CHECK(fixture.assembler->alignmentData[1].isDeleted0());
    CHECK(fixture.assembler->getRemovedSecondaryAlignmentBySymmetryOnlyCountForTesting() == 1);
}

// =============================================================================
// HELPER: Reverse complement sequence
// =============================================================================
std::string reverseComplement(const std::string& seq) {
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


// =============================================================================
// HELPER: Toggle base (for SNPs)
// =============================================================================
char otherBase(char b) {
    if (b == 'A') return 'C';
    return 'A';
}

// =============================================================================
// TEST: Basic Read Loading
// =============================================================================
TEST_CASE("Integration: Assembler read loading", "[integration][reads]") {
    AssemblerIntegrationFixture fixture;
    
    std::string seq1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";  // 64 bases
    std::string seq2 = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";  // 64 bases
    
    fixture.createFastq({seq1, seq2});
    fixture.initAssembler();
    fixture.loadReads();
    
    const auto& reads = fixture.assembler->getReads();
    REQUIRE(reads.readCount() == 2);
    CHECK(toString(reads.getRead(0)) == seq1);
    CHECK(toString(reads.getRead(1)) == seq2);
    CHECK(std::string(reads.getReadName(0).begin(), reads.getReadName(0).end()) == "read_0");
    CHECK(std::string(reads.getReadName(1).begin(), reads.getReadName(1).end()) == "read_1");
}

// =============================================================================
// TEST: Marker Generation
// =============================================================================
TEST_CASE("Integration: Marker generation from syncmers", "[integration][markers]") {
    AssemblerIntegrationFixture fixture;
    
    // 100 base sequence - should produce several markers
    std::string seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
                      "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";
    
    fixture.createFastq({seq});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(15, 5);
    
    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    // Read 0 strand 0
    const auto read0Strand0 = OrientedReadId(0, 0).getValue();
    const auto read0Strand1 = OrientedReadId(0, 1).getValue();
    const size_t count0 = markers->size(read0Strand0);
    const size_t count1 = markers->size(read0Strand1);
    REQUIRE(count0 > 0);
    CHECK(count0 == count1);

    // Markers are sorted and within bounds.
    const auto readMarkers = (*markers)[read0Strand0];
    for (size_t i = 0; i < readMarkers.size(); ++i) {
        CHECK(uint32_t(readMarkers[i].position) + 15 <= seq.size());
        if (i > 0) {
            CHECK(uint32_t(readMarkers[i - 1].position) <= uint32_t(readMarkers[i].position));
        }
    }
}

// =============================================================================
// TEST: K-mer Counting
// =============================================================================
TEST_CASE("Integration: K-mer counting", "[integration][counting]") {
    AssemblerIntegrationFixture fixture;
    
    // Create reads with varied k-mers to satisfy getHistogramHistogram logic
    // We need more repeated k-mers than unique ones to create a valley in the histogram.
    std::string seqSH = randomSequence(200, 123);
    std::string seqLG = randomSequence(1000, 456); 
    
    fixture.createFastq({seqLG, seqLG, seqSH});  // 2x Long (Freq 2), 1x Short (Freq 1)
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(15, 5);
    fixture.countKmers();
    
    auto& kmerCounter = *(fixture.assembler->kmerCounter);
    REQUIRE(kmerCounter.kmerIdFrequencies.size() > 0);

    bool foundFreq1 = false;
    bool foundFreq2Plus = false;
    uint64_t maxFreq = 0;
    for (auto& kv : kmerCounter.kmerIdFrequencies) {
        const uint64_t freq = kv.second;
        maxFreq = std::max(maxFreq, freq);
        if (freq == 1) foundFreq1 = true;
        if (freq >= 2) foundFreq2Plus = true;
    }
    CHECK(foundFreq1);      // seqSH contributes rare markers
    CHECK(foundFreq2Plus);  // seqLG duplicated
    CHECK(maxFreq >= 2);
}

// =============================================================================
// TEST: Marker Filtering
// =============================================================================
TEST_CASE("Integration: Marker filtering by frequency", "[integration][filtering]") {
    AssemblerIntegrationFixture fixture;
    
    // Create reads with varied k-mers to satisfy getHistogramInfo valley search
    // seqLG (repeated twice) creates Freq(2). seqSH (repeated once) creates Freq(1).
    // Ensure Freq(2) > Freq(1) by making seqLG much longer.
    std::string seqSH = randomSequence(200, 789);
    std::string seqLG = randomSequence(1000, 101112);
    
    fixture.createFastq({seqLG, seqLG, seqSH});
    fixture.initAssembler();
    fixture.loadReads();
    const uint64_t k = 15;
    fixture.generateMarkers(int(k), 5);
    fixture.countKmers();
    
    auto* markersBeforePtr = fixture.assembler->markers.get();
    REQUIRE(markersBeforePtr != nullptr);

    const uint64_t markersBefore = markersBeforePtr->totalSize();
    const auto read0s0 = markersBeforePtr->size(OrientedReadId(0, 0).getValue());
    const auto read2s0 = markersBeforePtr->size(OrientedReadId(2, 0).getValue());
    
    // Filter: keep only k-mers with frequency >= 2
    fixture.applyFilter(2, 1000);
    
    auto* markersAfterPtr = fixture.assembler->markers.get();
    REQUIRE(markersAfterPtr != nullptr);
    const uint64_t markersAfter = markersAfterPtr->totalSize();
    const auto read0s0After = markersAfterPtr->size(OrientedReadId(0, 0).getValue());
    const auto read2s0After = markersAfterPtr->size(OrientedReadId(2, 0).getValue());
    
    CHECK(markersAfter <= markersBefore);
    REQUIRE(markersAfter > 0);
    REQUIRE(read0s0After > 0);
    CHECK(read0s0After <= read0s0);

    // The unique read (read_2) should lose essentially all its markers under minFreq=2.
    CHECK(read2s0After == 0);

    // Every remaining marker must have k-mer frequency within the requested bounds.
    const auto& reads = fixture.assembler->getReads();
    auto& kmerCounter = *(fixture.assembler->kmerCounter);
    kmerCounter.buildFrequencyLUT();

    for (uint32_t orientedReadIdValue = 0; orientedReadIdValue < markersAfterPtr->size(); ++orientedReadIdValue) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
        const auto orientedMarkers = (*markersAfterPtr)[orientedReadIdValue];
        for (uint32_t ordinal = 0; ordinal < orientedMarkers.size(); ++ordinal) {
            const Kmer kmer = dinara::getOrientedReadMarkerKmer(
                orientedReadId, ordinal, k, reads, *markersAfterPtr);
            const KmerId id = KmerId(kmer.id(k));
            const KmerId rc = KmerId(kmer.reverseComplement(k).id(k));
            const KmerId canonical = std::min(id, rc);
            const uint64_t freq = kmerCounter.getFrequencyFast(canonical);
            CHECK(freq >= 2);
            CHECK(freq <= 1000);
        }
    }
}

// =============================================================================
// TEST: Forward/Reverse Complement Symmetry
// =============================================================================
TEST_CASE("Integration: F/RC marker symmetry", "[integration][symmetry]") {
    AssemblerIntegrationFixture fixture;
    
    std::string seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    
    fixture.createFastq({seq});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(15, 5);
    
    auto* markers = fixture.assembler->markers.get();
    
    size_t strand0 = markers->size(OrientedReadId(0, 0).getValue());
    size_t strand1 = markers->size(OrientedReadId(0, 1).getValue());
    
    CHECK(strand0 == strand1);
}
// =============================================================================
// TEST: Projected Alignment and Evidence Storage
// =============================================================================


TEST_CASE("Integration: Overlap-event anchors select vertices uniquely", "[integration][anchors][events]") {
    AssemblerIntegrationFixture fixture;

    // Two reads with enough markers. We will manually construct:
    // - a single read-graph overlap (and its RC twin)
    // - a minimal marker graph with two canonical vertices (and their RC vertices)
    // such that the overlap start/end events map to those vertices.
    fixture.createFastq({randomSequence(4000, 501), randomSequence(4000, 502)});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);

    withSilencedIoInDir(fixture.testDir, [&] {
        // -----------------------
        // Alignment + ReadGraph
        // -----------------------
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        AlignmentInfo info;
        info.alignmentId = 0;
        info.isInReadGraph = 1;
        info.markerCount = lastOrdinal + 1 - firstOrdinal;
        info.data[0] = AlignmentInfo::Data(m0, firstOrdinal, lastOrdinal);
        info.data[1] = AlignmentInfo::Data(m1, firstOrdinal, lastOrdinal);
        info.minOrdinalOffset = 0;
        info.maxOrdinalOffset = 0;
        info.averageOrdinalOffset = 0;
        info.maxSkip = 0;
        info.maxDrift = 0;

        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);

        fixture.assembler->readGraph.edges.resize(2);
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0p, r1p};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[0] = e;
        }
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0m, r1m};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[1] = e;
        }

        // connectivity is indexed by orientedReadIdValue = readId*2 + strand.
        // Index 0: r0+, 1: r0-, 2: r1+, 3: r1-.
        auto& conn = fixture.assembler->readGraph.connectivity;
        conn.appendVector(); // 0
        conn.append(0);
        conn.appendVector(); // 1
        conn.append(1);
        conn.appendVector(); // 2
        conn.append(0);
        conn.appendVector(); // 3
        conn.append(1);

        // -----------------------
        // Minimal MarkerGraph
        // -----------------------
        auto& mg = fixture.assembler->markerGraph;
        mg.constructVertices();
        mg.vertices().createNew("", 4096);
        mg.vertexTable.createNew("", 4096);
        mg.vertexTable.resize(markers->totalSize());
        std::fill(mg.vertexTable.begin(), mg.vertexTable.end(), MarkerGraph::invalidCompressedVertexId);

        // Two canonical vertices and their reverse complements:
        //  vertex0 <-> vertex2, vertex1 <-> vertex3.
        const MarkerGraphVertexId v0 = 0;
        const MarkerGraphVertexId v1 = 1;
        const MarkerGraphVertexId v2 = 2;
        const MarkerGraphVertexId v3 = 3;

        mg.reverseComplementVertex.createNew("", 4096);
        mg.reverseComplementVertex.resize(4);
        mg.reverseComplementVertex[v0] = v2;
        mg.reverseComplementVertex[v2] = v0;
        mg.reverseComplementVertex[v1] = v3;
        mg.reverseComplementVertex[v3] = v1;

        const MarkerId r0p_first = fixture.assembler->getMarkerId(r0p, firstOrdinal);
        const MarkerId r1p_first = fixture.assembler->getMarkerId(r1p, firstOrdinal);
        const MarkerId r0p_after = fixture.assembler->getMarkerId(r0p, lastOrdinal + 1);
        const MarkerId r1p_after = fixture.assembler->getMarkerId(r1p, lastOrdinal + 1);

        const uint32_t rcFirstOrdinal0 = uint32_t(m0) - 1 - firstOrdinal;
        const uint32_t rcFirstOrdinal1 = uint32_t(m1) - 1 - firstOrdinal;
        const uint32_t rcAfterOrdinal0 = uint32_t(m0) - 1 - (lastOrdinal + 1);
        const uint32_t rcAfterOrdinal1 = uint32_t(m1) - 1 - (lastOrdinal + 1);

        const MarkerId r0m_first = fixture.assembler->getMarkerId(r0m, rcFirstOrdinal0);
        const MarkerId r1m_first = fixture.assembler->getMarkerId(r1m, rcFirstOrdinal1);
        const MarkerId r0m_after = fixture.assembler->getMarkerId(r0m, rcAfterOrdinal0);
        const MarkerId r1m_after = fixture.assembler->getMarkerId(r1m, rcAfterOrdinal1);

        // Vertex 0: overlap "start" event markers.
        {
            vector<MarkerId> ids = {r0p_first, r1p_first};
            std::sort(ids.begin(), ids.end());
            mg.vertices().appendVector();
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(v0));
            }
        }
        // Vertex 1: overlap "end" event markers (last+1).
        {
            vector<MarkerId> ids = {r0p_after, r1p_after};
            std::sort(ids.begin(), ids.end());
            mg.vertices().appendVector();
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(v1));
            }
        }
        // Vertex 2: RC of vertex 0.
        {
            vector<MarkerId> ids = {r0m_first, r1m_first};
            std::sort(ids.begin(), ids.end());
            mg.vertices().appendVector();
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(v2));
            }
        }
        // Vertex 3: RC of vertex 1.
        {
            vector<MarkerId> ids = {r0m_after, r1m_after};
            std::sort(ids.begin(), ids.end());
            mg.vertices().appendVector();
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(v3));
            }
        }
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromMarkerGraphVerticesAtOverlapEvents(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() == 4); // 2 canonical vertices * (forward+RC)

    // Verify that each (forward, rc) pair is consecutive and consistent.
    for(mode3::AnchorId a = 0; a < anchors->size(); a += 2) {
        const auto forward = (*anchors)[a];
        const auto reverse = (*anchors)[a + 1];
        REQUIRE(forward.size() == reverse.size());
        // The reverse-complement anchor is stored sorted by OrientedReadId (like all anchors).
        std::vector<dinara::mode3::AnchorMarkerInterval> expectedReverse;
        expectedReverse.reserve(forward.size());
        for(const auto& f : forward) {
            auto expected = f;
            expected.orientedReadId.flipStrand();
            const uint64_t markerCount = fixture.assembler->markers->size(expected.orientedReadId.getValue());
            expected.ordinal0 = uint32_t(markerCount) - 1 - f.ordinal0;
            expectedReverse.push_back(expected);
        }
        std::sort(expectedReverse.begin(), expectedReverse.end(),
            [](const dinara::mode3::AnchorMarkerInterval& x, const dinara::mode3::AnchorMarkerInterval& y) {
                return x.orientedReadId < y.orientedReadId;
            });

        for(size_t i = 0; i < reverse.size(); ++i) {
            CHECK(reverse[i].orientedReadId == expectedReverse[i].orientedReadId);
            CHECK(reverse[i].ordinal0 == expectedReverse[i].ordinal0);
        }
    }
}

TEST_CASE("Integration: Best-per-interval anchors pick the max-coverage vertex inside the interval", "[integration][anchors][events][best]") {
    AssemblerIntegrationFixture fixture;

    // Four reads. We manually construct readGraph overlaps so that read 0 has one active overlap interval,
    // and construct a marker graph where:
    // - The highest-coverage interior vertex contains both strands for read 0 (duplicate ReadId) and must be skipped.
    // - The next-best interior vertex is clean and should be selected.
    fixture.createFastq({
        randomSequence(4000, 601),
        randomSequence(4000, 602),
        randomSequence(4000, 603),
        randomSequence(4000, 604),
    });
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const OrientedReadId r2p(ReadId(2), 0);
    const OrientedReadId r2m(ReadId(2), 1);
    const OrientedReadId r3p(ReadId(3), 0);
    const OrientedReadId r3m(ReadId(3), 1);

    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    const uint32_t m2 = uint32_t(markers->size(r2p.getValue()));
    const uint32_t m3 = uint32_t(markers->size(r3p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);
    REQUIRE(m2 > 16);
    REQUIRE(m3 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);
    REQUIRE(lastOrdinal + 1 < m2);
    REQUIRE(lastOrdinal + 1 < m3);

    const uint32_t interiorBadOrdinal = 2;
    const uint32_t interiorGoodOrdinal = 3;
    REQUIRE(interiorBadOrdinal >= firstOrdinal);
    REQUIRE(interiorBadOrdinal < lastOrdinal + 1);
    REQUIRE(interiorGoodOrdinal >= firstOrdinal);
    REQUIRE(interiorGoodOrdinal < lastOrdinal + 1);

    withSilencedIoInDir(fixture.testDir, [&] {
        // -----------------------
        // AlignmentData: (0,1), (0,2), (0,3) all same-strand, same ordinal interval.
        // -----------------------
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        auto fillAlignment = [&](size_t idx, ReadId otherReadId) {
            AlignmentInfo info;
            info.alignmentId = idx;
            info.isInReadGraph = 1;
            info.markerCount = lastOrdinal + 1 - firstOrdinal;

            info.data[0] = AlignmentInfo::Data(m0, firstOrdinal, lastOrdinal);
            uint32_t otherMarkerCount = 0;
            if(otherReadId == ReadId(1)) otherMarkerCount = m1;
            if(otherReadId == ReadId(2)) otherMarkerCount = m2;
            if(otherReadId == ReadId(3)) otherMarkerCount = m3;
            info.data[1] = AlignmentInfo::Data(otherMarkerCount, firstOrdinal, lastOrdinal);

            info.minOrdinalOffset = 0;
            info.maxOrdinalOffset = 0;
            info.averageOrdinalOffset = 0;
            info.maxSkip = 0;
            info.maxDrift = 0;

            AlignmentData ad(OrientedReadPair(ReadId(0), otherReadId, true), info);
            fixture.assembler->alignmentData[idx] = ad;
        };
        fillAlignment(0, ReadId(1));
        fillAlignment(1, ReadId(2));
        fillAlignment(2, ReadId(3));

        // -----------------------
        // ReadGraph edges: each alignment yields an edge for + strand and an edge for - strand.
        // We only need connectivity for read 0 (both strands) but keep it consistent for endpoints too.
        // -----------------------
        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);
        fixture.assembler->readGraph.edges.resize(6);

        auto setEdge = [&](size_t edgeIndex, OrientedReadId a, OrientedReadId b, uint64_t alignmentId) {
            ReadGraphEdge e;
            e.orientedReadIds = {a, b};
            e.alignmentId = alignmentId;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[edgeIndex] = e;
        };
        setEdge(0, r0p, r1p, 0);
        setEdge(1, r0m, r1m, 0);
        setEdge(2, r0p, r2p, 1);
        setEdge(3, r0m, r2m, 1);
        setEdge(4, r0p, r3p, 2);
        setEdge(5, r0m, r3m, 2);

        auto& conn = fixture.assembler->readGraph.connectivity;
        // orientedReadCount = 8 (4 reads * 2 strands).
        // Build vectors in order 0..7.
        conn.appendVector(); // 0 r0+
        conn.append(0); conn.append(2); conn.append(4);
        conn.appendVector(); // 1 r0-
        conn.append(1); conn.append(3); conn.append(5);

        conn.appendVector(); // 2 r1+
        conn.append(0);
        conn.appendVector(); // 3 r1-
        conn.append(1);

        conn.appendVector(); // 4 r2+
        conn.append(2);
        conn.appendVector(); // 5 r2-
        conn.append(3);

        conn.appendVector(); // 6 r3+
        conn.append(4);
        conn.appendVector(); // 7 r3-
        conn.append(5);

        // -----------------------
        // Minimal MarkerGraph:
        // - boundary at firstOrdinal maps to low-coverage vertex v0 (2 markers: r0,r1)
        // - boundary at afterLast maps to low-coverage vertex v2 (2 markers: r0,r1)
        // - interior ordinal maps to:
        //     v4: highest coverage but *invalid* (contains r0+ and r0- -> duplicate ReadId)
        //     v6: next-best coverage and valid (unique ReadIds)
        // RC pairs are (0<->1), (2<->3), (4<->5), (6<->7).
        // -----------------------
        auto& mg = fixture.assembler->markerGraph;
        mg.constructVertices();
        mg.vertices().createNew("", 4096);
        mg.vertexTable.createNew("", 4096);
        mg.vertexTable.resize(markers->totalSize());
        std::fill(mg.vertexTable.begin(), mg.vertexTable.end(), MarkerGraph::invalidCompressedVertexId);

        mg.reverseComplementVertex.createNew("", 4096);
        mg.reverseComplementVertex.resize(8);
        mg.reverseComplementVertex[0] = 1;
        mg.reverseComplementVertex[1] = 0;
        mg.reverseComplementVertex[2] = 3;
        mg.reverseComplementVertex[3] = 2;
        mg.reverseComplementVertex[4] = 5;
        mg.reverseComplementVertex[5] = 4;
        mg.reverseComplementVertex[6] = 7;
        mg.reverseComplementVertex[7] = 6;

        const MarkerId r0_first = fixture.assembler->getMarkerId(r0p, firstOrdinal);
        const MarkerId r1_first = fixture.assembler->getMarkerId(r1p, firstOrdinal);
        const MarkerId r0_after = fixture.assembler->getMarkerId(r0p, lastOrdinal + 1);
        const MarkerId r1_after = fixture.assembler->getMarkerId(r1p, lastOrdinal + 1);

        const MarkerId r0_bad = fixture.assembler->getMarkerId(r0p, interiorBadOrdinal);
        const MarkerId r1_bad = fixture.assembler->getMarkerId(r1p, interiorBadOrdinal);
        const MarkerId r2_bad = fixture.assembler->getMarkerId(r2p, interiorBadOrdinal);
        const MarkerId r3_bad = fixture.assembler->getMarkerId(r3p, interiorBadOrdinal);

        const uint32_t rcBadOrdinal0 = uint32_t(m0) - 1 - interiorBadOrdinal;
        const MarkerId r0m_bad = fixture.assembler->getMarkerId(r0m, rcBadOrdinal0);

        const MarkerId r0_good = fixture.assembler->getMarkerId(r0p, interiorGoodOrdinal);
        const MarkerId r1_good = fixture.assembler->getMarkerId(r1p, interiorGoodOrdinal);
        const MarkerId r2_good = fixture.assembler->getMarkerId(r2p, interiorGoodOrdinal);
        const MarkerId r3_good = fixture.assembler->getMarkerId(r3p, interiorGoodOrdinal);

        // Build 8 vertices in order.
        mg.vertices().clear();

        // Vertex 0: start boundary (low coverage, 2 markers).
        mg.vertices().appendVector();
        {
            std::vector<MarkerId> ids = {r0_first, r1_first};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(0));
            }
        }

        // Vertex 1: rc partner (unused content for this test).
        mg.vertices().appendVector();

        // Vertex 2: end boundary (low coverage, 2 markers).
        mg.vertices().appendVector();
        {
            std::vector<MarkerId> ids = {r0_after, r1_after};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(2));
            }
        }

        // Vertex 3: rc partner (unused content).
        mg.vertices().appendVector();

        // Vertex 4: interior (highest coverage but invalid due to duplicate ReadId: r0+ and r0-).
        mg.vertices().appendVector();
        {
            std::vector<MarkerId> ids = {r0_bad, r0m_bad, r1_bad, r2_bad, r3_bad};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(4));
            }
        }

        // Vertex 5: rc partner (unused content).
        mg.vertices().appendVector();

        // Vertex 6: interior (next-best, valid, unique ReadIds).
        mg.vertices().appendVector();
        {
            std::vector<MarkerId> ids = {r0_good, r1_good, r2_good, r3_good};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(6));
            }
        }

        // Vertex 7: rc partner for vertex 6 (valid).
        mg.vertices().appendVector();
        {
            const uint32_t rcGoodOrdinal0 = uint32_t(m0) - 1 - interiorGoodOrdinal;
            const uint32_t rcGoodOrdinal1 = uint32_t(m1) - 1 - interiorGoodOrdinal;
            const uint32_t rcGoodOrdinal2 = uint32_t(m2) - 1 - interiorGoodOrdinal;
            const uint32_t rcGoodOrdinal3 = uint32_t(m3) - 1 - interiorGoodOrdinal;

            const MarkerId r0m_good = fixture.assembler->getMarkerId(r0m, rcGoodOrdinal0);
            const MarkerId r1m_good = fixture.assembler->getMarkerId(r1m, rcGoodOrdinal1);
            const MarkerId r2m_good = fixture.assembler->getMarkerId(r2m, rcGoodOrdinal2);
            const MarkerId r3m_good = fixture.assembler->getMarkerId(r3m, rcGoodOrdinal3);

            std::vector<MarkerId> ids = {r0m_good, r1m_good, r2m_good, r3m_good};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(7));
            }
        }
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() == 2); // one selected vertex => (forward + rc)

    const auto a0 = (*anchors)[0];
    REQUIRE(a0.coverage() == 4);
    // All plus-strand read ids appear at the interior ordinal.
    {
        std::vector<OrientedReadId> ids;
        for(const auto& mi : a0) {
            ids.push_back(mi.orientedReadId);
        }
        std::sort(ids.begin(), ids.end(), [](const OrientedReadId& a, const OrientedReadId& b) {
            return a.getValue() < b.getValue();
        });
        REQUIRE(ids.size() == 4);
        CHECK(ids[0] == OrientedReadId(ReadId(0), 0));
        CHECK(ids[1] == OrientedReadId(ReadId(1), 0));
        CHECK(ids[2] == OrientedReadId(ReadId(2), 0));
        CHECK(ids[3] == OrientedReadId(ReadId(3), 0));
    }
}

TEST_CASE("Integration: Best-per-interval anchors fall back to duplicate-ReadId vertex when no clean vertex exists", "[integration][anchors][events][best][fallback]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(4000, 701), randomSequence(4000, 702)});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);

    const uint32_t interiorOrdinal = 2;
    const uint32_t rcInterior0 = uint32_t(m0) - 1 - interiorOrdinal;

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        AlignmentInfo info;
        info.alignmentId = 0;
        info.isInReadGraph = 1;
        info.markerCount = lastOrdinal + 1 - firstOrdinal;
        info.data[0] = AlignmentInfo::Data(m0, firstOrdinal, lastOrdinal);
        info.data[1] = AlignmentInfo::Data(m1, firstOrdinal, lastOrdinal);
        info.minOrdinalOffset = 0;
        info.maxOrdinalOffset = 0;
        info.averageOrdinalOffset = 0;
        info.maxSkip = 0;
        info.maxDrift = 0;

        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);
        fixture.assembler->readGraph.edges.resize(2);
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0p, r1p};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[0] = e;
        }
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0m, r1m};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[1] = e;
        }
        auto& conn = fixture.assembler->readGraph.connectivity;
        conn.appendVector(); // r0+
        conn.append(0);
        conn.appendVector(); // r0-
        conn.append(1);
        conn.appendVector(); // r1+
        conn.append(0);
        conn.appendVector(); // r1-
        conn.append(1);

        // Marker graph with only one interior vertex in range, but it has duplicate ReadId (r0+ and r0-).
        auto& mg = fixture.assembler->markerGraph;
        mg.constructVertices();
        mg.vertices().createNew("", 4096);
        mg.vertexTable.createNew("", 4096);
        mg.vertexTable.resize(markers->totalSize());
        std::fill(mg.vertexTable.begin(), mg.vertexTable.end(), MarkerGraph::invalidCompressedVertexId);

        mg.reverseComplementVertex.createNew("", 4096);
        mg.reverseComplementVertex.resize(2);
        mg.reverseComplementVertex[0] = 1;
        mg.reverseComplementVertex[1] = 0;

        mg.vertices().clear();
        mg.vertices().appendVector();
        {
            const MarkerId r0p_mid = fixture.assembler->getMarkerId(r0p, interiorOrdinal);
            const MarkerId r0m_mid = fixture.assembler->getMarkerId(r0m, rcInterior0);
            const MarkerId r1p_mid = fixture.assembler->getMarkerId(r1p, interiorOrdinal);
            std::vector<MarkerId> ids = {r0p_mid, r0m_mid, r1p_mid};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(0));
            }
        }
        mg.vertices().appendVector(); // rc partner unused
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() == 2); // selected vertex => (forward + rc)

    const auto a0 = (*anchors)[0];
    // Must contain both strands for read 0 (duplicate ReadId), because no clean alternative exists.
    bool hasR0p = false;
    bool hasR0m = false;
    for(const auto& mi : a0) {
        if(mi.orientedReadId == r0p) hasR0p = true;
        if(mi.orientedReadId == r0m) hasR0m = true;
    }
    CHECK(hasR0p);
    CHECK(hasR0m);
}

TEST_CASE("Integration: Best-per-interval anchor decomposition splits a vertex bridged by a single read", "[integration][anchors][events][best][decompose]") {
    AssemblerIntegrationFixture fixture;

    // Five reads. We construct a single marker graph vertex containing markers from all five,
    // but only one read (read 4) connects the two groups via overlaps. The decomposed anchor
    // method must split this into two anchors (plus their reverse complements).
    fixture.createFastq({
        randomSequence(4000, 801),
        randomSequence(4000, 802),
        randomSequence(4000, 803),
        randomSequence(4000, 804),
        randomSequence(4000, 805),
    });
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const OrientedReadId r2p(ReadId(2), 0);
    const OrientedReadId r2m(ReadId(2), 1);
    const OrientedReadId r3p(ReadId(3), 0);
    const OrientedReadId r3m(ReadId(3), 1);
    const OrientedReadId r4p(ReadId(4), 0);
    const OrientedReadId r4m(ReadId(4), 1);

    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    const uint32_t m2 = uint32_t(markers->size(r2p.getValue()));
    const uint32_t m3 = uint32_t(markers->size(r3p.getValue()));
    const uint32_t m4 = uint32_t(markers->size(r4p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);
    REQUIRE(m2 > 16);
    REQUIRE(m3 > 16);
    REQUIRE(m4 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    const uint32_t interiorOrdinal = 2;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);
    REQUIRE(lastOrdinal + 1 < m2);
    REQUIRE(lastOrdinal + 1 < m3);
    REQUIRE(lastOrdinal + 1 < m4);

    withSilencedIoInDir(fixture.testDir, [&] {
        // Alignments:
        // - Group A: (0,1)
        // - Group B: (2,3)
        // - Bridge read 4 overlaps to all in both groups.
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(6);

        auto fillAlignment = [&](size_t idx, ReadId a, ReadId b, uint32_t ma, uint32_t mb) {
            AlignmentInfo info;
            info.alignmentId = idx;
            info.isInReadGraph = 1;
            info.markerCount = lastOrdinal + 1 - firstOrdinal;
            info.data[0] = AlignmentInfo::Data(ma, firstOrdinal, lastOrdinal);
            info.data[1] = AlignmentInfo::Data(mb, firstOrdinal, lastOrdinal);
            info.minOrdinalOffset = 0;
            info.maxOrdinalOffset = 0;
            info.averageOrdinalOffset = 0;
            info.maxSkip = 0;
            info.maxDrift = 0;

            AlignmentData ad(OrientedReadPair(a, b, true), info);
            fixture.assembler->alignmentData[idx] = ad;
        };

        fillAlignment(0, ReadId(0), ReadId(1), m0, m1);
        fillAlignment(1, ReadId(2), ReadId(3), m2, m3);
        fillAlignment(2, ReadId(0), ReadId(4), m0, m4);
        fillAlignment(3, ReadId(1), ReadId(4), m1, m4);
        fillAlignment(4, ReadId(2), ReadId(4), m2, m4);
        fillAlignment(5, ReadId(3), ReadId(4), m3, m4);

        // ReadGraph edges: each alignment yields an edge for + and an edge for -.
        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);
        fixture.assembler->readGraph.edges.resize(12);

        auto setEdge = [&](size_t edgeIndex, OrientedReadId a, OrientedReadId b, uint64_t alignmentId) {
            ReadGraphEdge e;
            e.orientedReadIds = {a, b};
            e.alignmentId = alignmentId;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[edgeIndex] = e;
        };
        // Alignment 0: (0,1)
        setEdge(0, r0p, r1p, 0);
        setEdge(1, r0m, r1m, 0);
        // Alignment 1: (2,3)
        setEdge(2, r2p, r3p, 1);
        setEdge(3, r2m, r3m, 1);
        // Alignment 2: (0,4)
        setEdge(4, r0p, r4p, 2);
        setEdge(5, r0m, r4m, 2);
        // Alignment 3: (1,4)
        setEdge(6, r1p, r4p, 3);
        setEdge(7, r1m, r4m, 3);
        // Alignment 4: (2,4)
        setEdge(8, r2p, r4p, 4);
        setEdge(9, r2m, r4m, 4);
        // Alignment 5: (3,4)
        setEdge(10, r3p, r4p, 5);
        setEdge(11, r3m, r4m, 5);

        auto& conn = fixture.assembler->readGraph.connectivity;
        // orientedReadCount = 10 (5 reads * 2 strands). Build vectors in order 0..9.
        conn.appendVector(); // 0 r0+
        conn.append(0); conn.append(4);
        conn.appendVector(); // 1 r0-
        conn.append(1); conn.append(5);

        conn.appendVector(); // 2 r1+
        conn.append(0); conn.append(6);
        conn.appendVector(); // 3 r1-
        conn.append(1); conn.append(7);

        conn.appendVector(); // 4 r2+
        conn.append(2); conn.append(8);
        conn.appendVector(); // 5 r2-
        conn.append(3); conn.append(9);

        conn.appendVector(); // 6 r3+
        conn.append(2); conn.append(10);
        conn.appendVector(); // 7 r3-
        conn.append(3); conn.append(11);

        conn.appendVector(); // 8 r4+
        conn.append(4); conn.append(6); conn.append(8); conn.append(10);
        conn.appendVector(); // 9 r4-
        conn.append(5); conn.append(7); conn.append(9); conn.append(11);

        // MarkerGraph: one vertex that contains all five reads on the + strand at interiorOrdinal.
        auto& mg = fixture.assembler->markerGraph;
        mg.constructVertices();
        mg.vertices().createNew("", 4096);
        mg.vertexTable.createNew("", 4096);
        mg.vertexTable.resize(markers->totalSize());
        std::fill(mg.vertexTable.begin(), mg.vertexTable.end(), MarkerGraph::invalidCompressedVertexId);

        mg.reverseComplementVertex.createNew("", 4096);
        mg.reverseComplementVertex.resize(2);
        mg.reverseComplementVertex[0] = 1;
        mg.reverseComplementVertex[1] = 0;

        mg.vertices().clear();
        mg.vertices().appendVector();
        {
            const MarkerId r0_mid = fixture.assembler->getMarkerId(r0p, interiorOrdinal);
            const MarkerId r1_mid = fixture.assembler->getMarkerId(r1p, interiorOrdinal);
            const MarkerId r2_mid = fixture.assembler->getMarkerId(r2p, interiorOrdinal);
            const MarkerId r3_mid = fixture.assembler->getMarkerId(r3p, interiorOrdinal);
            const MarkerId r4_mid = fixture.assembler->getMarkerId(r4p, interiorOrdinal);
            std::vector<MarkerId> ids = {r0_mid, r1_mid, r2_mid, r3_mid, r4_mid};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(0));
            }
        }
        mg.vertices().appendVector(); // rc partner unused
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromMarkerGraphVerticesBestPerOverlapIntervalDecomposed(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() == 4); // two decomposed anchors => (forward + rc) each

    bool foundGroupA = false;
    bool foundGroupB = false;
    for(dinara::mode3::AnchorId anchorId=0; anchorId<anchors->size(); ++anchorId) {
        const auto a = (*anchors)[anchorId];
        bool hasR0 = false, hasR1 = false, hasR2 = false, hasR3 = false, hasR4 = false;
        for(const auto& mi : a) {
            if(mi.orientedReadId == r0p) hasR0 = true;
            if(mi.orientedReadId == r1p) hasR1 = true;
            if(mi.orientedReadId == r2p) hasR2 = true;
            if(mi.orientedReadId == r3p) hasR3 = true;
            if(mi.orientedReadId == r4p) hasR4 = true;
        }
        if(hasR0) {
            CHECK(hasR1);
            CHECK_FALSE(hasR2);
            CHECK_FALSE(hasR3);
            CHECK_FALSE(hasR4);
            foundGroupA = true;
        }
        if(hasR2) {
            CHECK(hasR3);
            CHECK_FALSE(hasR0);
            CHECK_FALSE(hasR1);
            CHECK_FALSE(hasR4);
            foundGroupB = true;
        }
    }
    CHECK(foundGroupA);
    CHECK(foundGroupB);
}

TEST_CASE("Integration: Best-per-interval anchor decomposition splits a vertex connected by multiple weak cross-edges", "[integration][anchors][events][best][decompose][multibridge]") {
    AssemblerIntegrationFixture fixture;

    // Six reads. Two tight clusters (0,1,2) and (3,4,5) with two weak cross-edges (0-3 and 1-4).
    // This graph has no articulation point, but cross-edges have divergent neighbor sets and should be dropped.
    fixture.createFastq({
        randomSequence(4000, 901),
        randomSequence(4000, 902),
        randomSequence(4000, 903),
        randomSequence(4000, 904),
        randomSequence(4000, 905),
        randomSequence(4000, 906),
    });
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const OrientedReadId r2p(ReadId(2), 0);
    const OrientedReadId r2m(ReadId(2), 1);
    const OrientedReadId r3p(ReadId(3), 0);
    const OrientedReadId r3m(ReadId(3), 1);
    const OrientedReadId r4p(ReadId(4), 0);
    const OrientedReadId r4m(ReadId(4), 1);
    const OrientedReadId r5p(ReadId(5), 0);
    const OrientedReadId r5m(ReadId(5), 1);

    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    const uint32_t m2 = uint32_t(markers->size(r2p.getValue()));
    const uint32_t m3 = uint32_t(markers->size(r3p.getValue()));
    const uint32_t m4 = uint32_t(markers->size(r4p.getValue()));
    const uint32_t m5 = uint32_t(markers->size(r5p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);
    REQUIRE(m2 > 16);
    REQUIRE(m3 > 16);
    REQUIRE(m4 > 16);
    REQUIRE(m5 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    const uint32_t interiorOrdinal = 2;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);
    REQUIRE(lastOrdinal + 1 < m2);
    REQUIRE(lastOrdinal + 1 < m3);
    REQUIRE(lastOrdinal + 1 < m4);
    REQUIRE(lastOrdinal + 1 < m5);

    withSilencedIoInDir(fixture.testDir, [&] {
        // Alignments for cluster A.
        // Alignments for cluster B.
        // Two cross edges: (0,3) and (1,4).
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(8);

        auto fillAlignment = [&](size_t idx, ReadId a, ReadId b, uint32_t ma, uint32_t mb) {
            AlignmentInfo info;
            info.alignmentId = idx;
            info.isInReadGraph = 1;
            info.markerCount = lastOrdinal + 1 - firstOrdinal;
            info.data[0] = AlignmentInfo::Data(ma, firstOrdinal, lastOrdinal);
            info.data[1] = AlignmentInfo::Data(mb, firstOrdinal, lastOrdinal);
            info.minOrdinalOffset = 0;
            info.maxOrdinalOffset = 0;
            info.averageOrdinalOffset = 0;
            info.maxSkip = 0;
            info.maxDrift = 0;

            AlignmentData ad(OrientedReadPair(a, b, true), info);
            fixture.assembler->alignmentData[idx] = ad;
        };

        // A clique: 0-1, 0-2, 1-2
        fillAlignment(0, ReadId(0), ReadId(1), m0, m1);
        fillAlignment(1, ReadId(0), ReadId(2), m0, m2);
        fillAlignment(2, ReadId(1), ReadId(2), m1, m2);
        // B clique: 3-4, 3-5, 4-5
        fillAlignment(3, ReadId(3), ReadId(4), m3, m4);
        fillAlignment(4, ReadId(3), ReadId(5), m3, m5);
        fillAlignment(5, ReadId(4), ReadId(5), m4, m5);
        // Cross edges
        fillAlignment(6, ReadId(0), ReadId(3), m0, m3);
        fillAlignment(7, ReadId(1), ReadId(4), m1, m4);

        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);
        fixture.assembler->readGraph.edges.resize(16);

        auto setEdge = [&](size_t edgeIndex, OrientedReadId a, OrientedReadId b, uint64_t alignmentId) {
            ReadGraphEdge e;
            e.orientedReadIds = {a, b};
            e.alignmentId = alignmentId;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[edgeIndex] = e;
        };

        // For each alignment, add + and -.
        setEdge(0, r0p, r1p, 0); setEdge(1, r0m, r1m, 0);
        setEdge(2, r0p, r2p, 1); setEdge(3, r0m, r2m, 1);
        setEdge(4, r1p, r2p, 2); setEdge(5, r1m, r2m, 2);
        setEdge(6, r3p, r4p, 3); setEdge(7, r3m, r4m, 3);
        setEdge(8, r3p, r5p, 4); setEdge(9, r3m, r5m, 4);
        setEdge(10, r4p, r5p, 5); setEdge(11, r4m, r5m, 5);
        setEdge(12, r0p, r3p, 6); setEdge(13, r0m, r3m, 6);
        setEdge(14, r1p, r4p, 7); setEdge(15, r1m, r4m, 7);

        auto& conn = fixture.assembler->readGraph.connectivity;
        // orientedReadCount = 12 (6 reads * 2 strands) in order 0..11.
        conn.appendVector(); // r0+
        conn.append(0); conn.append(2); conn.append(12);
        conn.appendVector(); // r0-
        conn.append(1); conn.append(3); conn.append(13);

        conn.appendVector(); // r1+
        conn.append(0); conn.append(4); conn.append(14);
        conn.appendVector(); // r1-
        conn.append(1); conn.append(5); conn.append(15);

        conn.appendVector(); // r2+
        conn.append(2); conn.append(4);
        conn.appendVector(); // r2-
        conn.append(3); conn.append(5);

        conn.appendVector(); // r3+
        conn.append(6); conn.append(8); conn.append(12);
        conn.appendVector(); // r3-
        conn.append(7); conn.append(9); conn.append(13);

        conn.appendVector(); // r4+
        conn.append(6); conn.append(10); conn.append(14);
        conn.appendVector(); // r4-
        conn.append(7); conn.append(11); conn.append(15);

        conn.appendVector(); // r5+
        conn.append(8); conn.append(10);
        conn.appendVector(); // r5-
        conn.append(9); conn.append(11);

        // MarkerGraph: one vertex contains all six reads on + strand at interiorOrdinal.
        auto& mg = fixture.assembler->markerGraph;
        mg.constructVertices();
        mg.vertices().createNew("", 4096);
        mg.vertexTable.createNew("", 4096);
        mg.vertexTable.resize(markers->totalSize());
        std::fill(mg.vertexTable.begin(), mg.vertexTable.end(), MarkerGraph::invalidCompressedVertexId);

        mg.reverseComplementVertex.createNew("", 4096);
        mg.reverseComplementVertex.resize(2);
        mg.reverseComplementVertex[0] = 1;
        mg.reverseComplementVertex[1] = 0;

        mg.vertices().clear();
        mg.vertices().appendVector();
        {
            const MarkerId r0_mid = fixture.assembler->getMarkerId(r0p, interiorOrdinal);
            const MarkerId r1_mid = fixture.assembler->getMarkerId(r1p, interiorOrdinal);
            const MarkerId r2_mid = fixture.assembler->getMarkerId(r2p, interiorOrdinal);
            const MarkerId r3_mid = fixture.assembler->getMarkerId(r3p, interiorOrdinal);
            const MarkerId r4_mid = fixture.assembler->getMarkerId(r4p, interiorOrdinal);
            const MarkerId r5_mid = fixture.assembler->getMarkerId(r5p, interiorOrdinal);
            std::vector<MarkerId> ids = {r0_mid, r1_mid, r2_mid, r3_mid, r4_mid, r5_mid};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(0));
            }
        }
        mg.vertices().appendVector(); // rc partner unused
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromMarkerGraphVerticesBestPerOverlapIntervalDecomposed(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() == 4); // two anchors (A and B) + RC each

    bool foundA = false;
    bool foundB = false;
    for(dinara::mode3::AnchorId anchorId=0; anchorId<anchors->size(); ++anchorId) {
        const auto a = (*anchors)[anchorId];
        bool has0=false, has1=false, has2=false, has3=false, has4=false, has5=false;
        for(const auto& mi : a) {
            if(mi.orientedReadId == r0p) has0 = true;
            if(mi.orientedReadId == r1p) has1 = true;
            if(mi.orientedReadId == r2p) has2 = true;
            if(mi.orientedReadId == r3p) has3 = true;
            if(mi.orientedReadId == r4p) has4 = true;
            if(mi.orientedReadId == r5p) has5 = true;
        }
        if(has0 || has1 || has2) {
            CHECK(has0);
            CHECK(has1);
            CHECK(has2);
            CHECK_FALSE(has3);
            CHECK_FALSE(has4);
            CHECK_FALSE(has5);
            foundA = true;
        }
        if(has3 || has4 || has5) {
            CHECK_FALSE(has0);
            CHECK_FALSE(has1);
            CHECK_FALSE(has2);
            CHECK(has3);
            CHECK(has4);
            CHECK(has5);
            foundB = true;
        }
    }
    CHECK(foundA);
    CHECK(foundB);
}

TEST_CASE("Integration: Overlap-only best-per-interval anchors work without marker graph vertices", "[integration][anchors][overlaps][best]") {
    AssemblerIntegrationFixture fixture;

    // Two identical reads => markers (k-mers) match at the same ordinals.
    const std::string seq = randomSequence(4000, 1001);
    fixture.createFastq({seq, seq});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);

    withSilencedIoInDir(fixture.testDir, [&] {
        // Alignment + ReadGraph (same as the overlap-event test, but we do not create a marker graph).
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        AlignmentInfo info;
        info.alignmentId = 0;
        info.isInReadGraph = 1;
        info.markerCount = lastOrdinal + 1 - firstOrdinal;
        info.data[0] = AlignmentInfo::Data(m0, firstOrdinal, lastOrdinal);
        info.data[1] = AlignmentInfo::Data(m1, firstOrdinal, lastOrdinal);
        info.minOrdinalOffset = 0;
        info.maxOrdinalOffset = 0;
        info.averageOrdinalOffset = 0;
        info.maxSkip = 0;
        info.maxDrift = 0;

        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);
        fixture.assembler->readGraph.edges.resize(2);
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0p, r1p};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[0] = e;
        }
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0m, r1m};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[1] = e;
        }

        auto& conn = fixture.assembler->readGraph.connectivity;
        conn.appendVector(); // r0+
        conn.append(0);
        conn.appendVector(); // r0-
        conn.append(1);
        conn.appendVector(); // r1+
        conn.append(0);
        conn.appendVector(); // r1-
        conn.append(1);
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromOverlapsBestPerOverlapInterval(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() >= 2); // at least one anchor + its RC anchor

    // Verify we can find at least one anchor containing both reads on strand 0 and its RC on strand 1.
    bool foundForward = false;
    bool foundReverse = false;
    for(uint64_t i=0; i<anchors->size(); ++i) {
        const auto a = (*anchors)[i];
        if(a.coverage() != 2) {
            continue;
        }
        if(a[0].orientedReadId == r0p && a[1].orientedReadId == r1p) {
            foundForward = true;
        }
        if(a[0].orientedReadId == r0m && a[1].orientedReadId == r1m) {
            foundReverse = true;
        }
        if(foundForward && foundReverse) {
            break;
        }
    }
    CHECK(foundForward);
    CHECK(foundReverse);
}

TEST_CASE("Integration: ma_hit_flt keeps contained overlaps", "[integration][hifiasm][filter]") {
    AssemblerIntegrationFixture fixture;

    // 3 reads: 0 is contained in 1, and 2 has an internal overlap to 1.
    fixture.createFastq({randomSequence(1000, 1), randomSequence(1500, 2), randomSequence(1000, 3)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(2);

        // Alignment 0: QCONT (query read 0 is fully contained in target read 1).
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 1000;
            ad.ts = 250; ad.te = 1250;
            fixture.assembler->alignmentData[0] = ad;
        }

        // Alignment 1: internal match (should be removed by ma_hit_flt).
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(2), ReadId(1), true), info);
            ad.qs = 200; ad.qe = 800;   // internal on query
            ad.ts = 100; ad.te = 700;   // internal on target
            fixture.assembler->alignmentData[1] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();

        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->applyCoverageCuts(50, 1);
        fixture.assembler->filterHangingOverlaps(1000, 0.8, 50, 1);
    });

    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted());
    CHECK(fixture.assembler->alignmentData[1].isDeleted());
    CHECK((fixture.assembler->alignmentData[1].deleteReasons0 & AlignmentData::DeleteReasonHanging) != 0);
    CHECK((fixture.assembler->alignmentData[1].deleteReasons1 & AlignmentData::DeleteReasonHanging) != 0);
}

TEST_CASE("Integration: ma_hit_cut skips deleted endpoints", "[integration][hifiasm][filter]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 11), randomSequence(1000, 12)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        AlignmentInfo info;
        info.alignmentId = 0;
        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
        ad.qs = 0; ad.qe = 1000;
        ad.ts = 0; ad.te = 1000;
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->computeAlignmentTableForTesting();

        // With min_dp=2, a single overlap yields depth 1, so both reads get deleted in ma_hit_sub.
        fixture.assembler->filterLocalSegments(2, 1);

        // Hifiasm parity: ma_hit_cut skips overlaps if either endpoint read is deleted; it does not force-delete the edge.
        fixture.assembler->applyCoverageCuts(50, 1);
    });

    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted());
    CHECK(fixture.assembler->alignmentData[0].deleteReasons0 == AlignmentData::DeleteReasonNone);
    CHECK(fixture.assembler->alignmentData[0].deleteReasons1 == AlignmentData::DeleteReasonNone);
}

TEST_CASE("Integration: ma_hit_sub selects the max-depth interval", "[integration][hifiasm][filter][sub]") {
    AssemblerIntegrationFixture fixture;

    // Read 0 has depth>=2 only on [400,600); reads 1 and 2 have only depth 1 everywhere.
    fixture.createFastq({randomSequence(1000, 201), randomSequence(1000, 202), randomSequence(1000, 203)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(2);

        // Overlap 0: [0,600) on read 0.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 600;
            ad.ts = 0; ad.te = 600;
            fixture.assembler->alignmentData[0] = ad;
        }

        // Overlap 1: [400,1000) on read 0.
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(2), true), info);
            ad.qs = 400; ad.qe = 1000;
            ad.ts = 0; ad.te = 600;
            fixture.assembler->alignmentData[1] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(2, 1);
    });

    const auto v0 = fixture.assembler->getValidReadIntervalForTesting(ReadId(0));
    CHECK_FALSE(v0.isDeleted);
    CHECK(v0.start == 400);
    CHECK(v0.end == 600);

    const auto v1 = fixture.assembler->getValidReadIntervalForTesting(ReadId(1));
    const auto v2 = fixture.assembler->getValidReadIntervalForTesting(ReadId(2));
    CHECK(v1.isDeleted);
    CHECK(v2.isDeleted);
}

TEST_CASE("Integration: ma_hit_cut trims and normalizes overlaps to valid segments", "[integration][hifiasm][filter][cut]") {
    AssemblerIntegrationFixture fixture;

    // Build a case where read 0's valid region is [100,900] and read 1's valid region is [200,800].
    // The overlap (0,1) is full-length in raw coordinates and must be trimmed + normalized.
    fixture.createFastq({
        randomSequence(1000, 301),
        randomSequence(1000, 302),
        randomSequence(1000, 303),
        randomSequence(1000, 304),
    });
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        // Full-length overlap between reads 0 and 1.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 1000;
            ad.ts = 0; ad.te = 1000;
            fixture.assembler->alignmentData[0] = ad;
        }

        // Second overlap to give read 0 depth>=2 only on [100,900].
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(2), true), info);
            ad.qs = 100; ad.qe = 900;
            ad.ts = 0; ad.te = 800;
            fixture.assembler->alignmentData[1] = ad;
        }

        // Second overlap to give read 1 depth>=2 only on [200,800].
        {
            AlignmentInfo info;
            info.alignmentId = 2;
            AlignmentData ad(OrientedReadPair(ReadId(1), ReadId(3), true), info);
            ad.qs = 200; ad.qe = 800;
            ad.ts = 0; ad.te = 600;
            fixture.assembler->alignmentData[2] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(2, 1);

        const auto v0 = fixture.assembler->getValidReadIntervalForTesting(ReadId(0));
        const auto v1 = fixture.assembler->getValidReadIntervalForTesting(ReadId(1));
        REQUIRE_FALSE(v0.isDeleted);
        REQUIRE_FALSE(v1.isDeleted);
        REQUIRE(v0.start == 100);
        REQUIRE(v0.end == 900);
        REQUIRE(v1.start == 200);
        REQUIRE(v1.end == 800);

        fixture.assembler->applyCoverageCuts(50, 1);
    });

    const AlignmentData& a01 = fixture.assembler->alignmentData[0];
    REQUIRE_FALSE(a01.isDeleted());
    CHECK(a01.deleteReasons0 == AlignmentData::DeleteReasonNone);
    CHECK(a01.deleteReasons1 == AlignmentData::DeleteReasonNone);

    // After trimming and normalization:
    // - read0 segment is [100,900], read1 segment is [200,800]
    // - intersection in raw coordinates is [200,800]
    // - normalized: on read0 => [100,700], on read1 => [0,600]
    CHECK(a01.qs == 100);
    CHECK(a01.qe == 700);
    CHECK(a01.ts == 0);
    CHECK(a01.te == 600);
}

TEST_CASE("Integration: ma_hit_contained_advance does not RC-map swapped intervals", "[integration][hifiasm][filter]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 41), randomSequence(1000, 42)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        // One reverse-strand overlap. When viewed from read 1's perspective, the query/target
        // intervals must be swapped without RC-mapping (hifiasm set_reverse_overlap parity).
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), false), info);
            ad.qs = 0; ad.qe = 900;
            ad.ts = 0; ad.te = 900;
            fixture.assembler->alignmentData[0] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();

        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->applyCoverageCuts(50, 1);
        fixture.assembler->filterHangingOverlaps(1000, 0.8, 50, 1);
        fixture.assembler->removeContainedReads(1000, 0.8, 50, 1);
    });

    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted());
    CHECK((fixture.assembler->alignmentData[0].deleteReasons0 & AlignmentData::DeleteReasonContained) == 0);
    CHECK((fixture.assembler->alignmentData[0].deleteReasons1 & AlignmentData::DeleteReasonContained) == 0);
}

TEST_CASE("Integration: try_rescue_overlaps clears directional phase deletions on consensus span", "[integration][hifiasm][filter][rescue]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({
        randomSequence(1000, 501),
        randomSequence(1000, 502),
        randomSequence(1000, 503),
        randomSequence(1000, 504),
        randomSequence(1000, 505),
    });
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(4);

        // Four conflict overlaps: read 0 deletes (phase) but other endpoints keep.
        // All intervals share the consensus intersection [300,650).
        const std::array<std::pair<uint32_t, uint32_t>, 4> qIntervals = {{
            {0, 700},
            {100, 900},
            {200, 800},
            {300, 650},
        }};
        for (uint32_t i = 0; i < 4; ++i) {
            AlignmentInfo info;
            info.alignmentId = i;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1 + i), true), info);
            ad.qs = qIntervals[i].first;
            ad.qe = qIntervals[i].second;
            ad.ts = 0;
            ad.te = qIntervals[i].second - qIntervals[i].first;
            ad.addDeleteReasonsFromReadPerspective(ReadId(0), AlignmentData::DeleteReasonPhase);
            fixture.assembler->alignmentData[i] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->rescuePhasedOverlaps(4, 1);
    });

    for (uint32_t i = 0; i < 4; ++i) {
        const auto& ad = fixture.assembler->alignmentData[i];
        CHECK_FALSE(ad.isDeleted0());
        CHECK_FALSE(ad.isDeleted1());
        CHECK(ad.deleteReasons0 == AlignmentData::DeleteReasonNone);
        CHECK(ad.deleteReasons1 == AlignmentData::DeleteReasonNone);
    }
}

TEST_CASE("Integration: ma_hit_contained_advance compresses containment chains", "[integration][hifiasm][filter][contained][chain]") {
    AssemblerIntegrationFixture fixture;

    // Chain: 0 contained in 1, 1 contained in 2. Keep (2,3) as a dovetail so container stays alive.
    fixture.createFastq({
        randomSequence(500, 601),
        randomSequence(800, 602),
        randomSequence(1200, 603),
        randomSequence(1200, 604),
    });
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        // 0 in 1.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 500;
            ad.ts = 100; ad.te = 600;
            fixture.assembler->alignmentData[0] = ad;
        }
        // 1 in 2.
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(1), ReadId(2), true), info);
            ad.qs = 0; ad.qe = 800;
            ad.ts = 200; ad.te = 1000;
            fixture.assembler->alignmentData[1] = ad;
        }
        // 2 dovetail 3 (kept, and NOT a containment).
        {
            AlignmentInfo info;
            info.alignmentId = 2;
            AlignmentData ad(OrientedReadPair(ReadId(2), ReadId(3), true), info);
            ad.qs = 0; ad.qe = 1100;
            ad.ts = 100; ad.te = 1200;
            fixture.assembler->alignmentData[2] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->applyCoverageCuts(50, 1);
        fixture.assembler->filterHangingOverlaps(1000, 0.8, 50, 1);
        fixture.assembler->removeContainedReads(1000, 0.8, 50, 1);
    });

    CHECK(fixture.assembler->getContainmentRootForTesting(ReadId(0)) == ReadId(2));
    CHECK(fixture.assembler->getContainmentRootForTesting(ReadId(1)) == ReadId(2));
    CHECK(fixture.assembler->getContainmentRootForTesting(ReadId(2)) == ReadId(invalidReadId));

    // Contained overlaps are deleted with the contained reason.
    CHECK((fixture.assembler->alignmentData[0].deleteReasons0 & AlignmentData::DeleteReasonContained) != 0);
    CHECK((fixture.assembler->alignmentData[0].deleteReasons1 & AlignmentData::DeleteReasonContained) != 0);
    CHECK((fixture.assembler->alignmentData[1].deleteReasons0 & AlignmentData::DeleteReasonContained) != 0);
    CHECK((fixture.assembler->alignmentData[1].deleteReasons1 & AlignmentData::DeleteReasonContained) != 0);

    // The dovetail edge (2,3) remains.
    CHECK_FALSE(fixture.assembler->alignmentData[2].isDeleted());
}

TEST_CASE("Integration: detect_chimeric_reads deletes edges for simple chimera", "[integration][hifiasm][filter]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 21), randomSequence(1000, 22), randomSequence(1000, 23)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        // For read 0 (length 1000), create a left anchor (qs==0) and right anchor (qe==rLen),
        // with a gap between them -> simple chimera.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 400;
            ad.ts = 0; ad.te = 400;
            fixture.assembler->alignmentData[0] = ad;
        }
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(2), true), info);
            ad.qs = 600; ad.qe = 1000;
            ad.ts = 0; ad.te = 400;
            fixture.assembler->alignmentData[1] = ad;
        }
        // Control overlap not involving read 0.
        {
            AlignmentInfo info;
            info.alignmentId = 2;
            AlignmentData ad(OrientedReadPair(ReadId(1), ReadId(2), true), info);
            ad.qs = 0; ad.qe = 1000;
            ad.ts = 0; ad.te = 1000;
            fixture.assembler->alignmentData[2] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->detectChimericReads(1);
    });

    CHECK(fixture.assembler->alignmentData[0].isDeleted());
    CHECK(fixture.assembler->alignmentData[1].isDeleted());
    CHECK_FALSE(fixture.assembler->alignmentData[2].isDeleted());
    CHECK((fixture.assembler->alignmentData[0].deleteReasons0 & AlignmentData::DeleteReasonChimeric) != 0);
    CHECK((fixture.assembler->alignmentData[0].deleteReasons1 & AlignmentData::DeleteReasonChimeric) != 0);
    CHECK((fixture.assembler->alignmentData[1].deleteReasons0 & AlignmentData::DeleteReasonChimeric) != 0);
    CHECK((fixture.assembler->alignmentData[1].deleteReasons1 & AlignmentData::DeleteReasonChimeric) != 0);
}

TEST_CASE("Integration: ONT chemical arc mask deletes overlaps for low-depth reads", "[integration][hifiasm][filter][ont][chemical]") {
    AssemblerIntegrationFixture fixture;

    // 4 reads. Read 0 has two full-length overlaps (min depth >=2) => not flagged.
    // Read 3 has two non-overlapping intervals => min depth == 0 => flagged.
    fixture.createFastq({
        randomSequence(2000, 31),
        randomSequence(2000, 32),
        randomSequence(2000, 33),
        randomSequence(2000, 34),
    });
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(4);

        // Read 0 overlaps (0,1) and (0,2) spanning the full read.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 2000;
            ad.ts = 0; ad.te = 2000;
            fixture.assembler->alignmentData[0] = ad;
        }
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(2), true), info);
            ad.qs = 0; ad.qe = 2000;
            ad.ts = 0; ad.te = 2000;
            fixture.assembler->alignmentData[1] = ad;
        }

        // Read 3 has two disjoint overlaps (gap in the middle).
        {
            AlignmentInfo info;
            info.alignmentId = 2;
            AlignmentData ad(OrientedReadPair(ReadId(1), ReadId(3), true), info);
            // Interval on read 3: [0,800)
            ad.qs = 0; ad.qe = 2000;
            ad.ts = 0; ad.te = 800;
            fixture.assembler->alignmentData[2] = ad;
        }
        {
            AlignmentInfo info;
            info.alignmentId = 3;
            AlignmentData ad(OrientedReadPair(ReadId(2), ReadId(3), true), info);
            // Interval on read 3: [1200,2000)
            ad.qs = 0; ad.qe = 2000;
            ad.ts = 1200; ad.te = 2000;
            fixture.assembler->alignmentData[3] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();

        // Use a small flank and cov=1 to make the test robust on short reads.
        fixture.assembler->applyOntChemicalArcMask(1, 0, 0.02, 1);
    });

    // Overlaps incident to read 3 are chemically deleted.
    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted());
    CHECK_FALSE(fixture.assembler->alignmentData[1].isDeleted());
    CHECK(fixture.assembler->alignmentData[2].isDeleted());
    CHECK(fixture.assembler->alignmentData[3].isDeleted());
    CHECK((fixture.assembler->alignmentData[2].deleteReasons0 & AlignmentData::DeleteReasonChemical) != 0);
    CHECK((fixture.assembler->alignmentData[2].deleteReasons1 & AlignmentData::DeleteReasonChemical) != 0);
    CHECK((fixture.assembler->alignmentData[3].deleteReasons0 & AlignmentData::DeleteReasonChemical) != 0);
    CHECK((fixture.assembler->alignmentData[3].deleteReasons1 & AlignmentData::DeleteReasonChemical) != 0);
}




















namespace {

uint64_t addDagSimpleNode(
    dinara::mode3::DirectedAnchorGraph& dag,
    uint64_t markerId,
    uint64_t lengthBp,
    double coverage)
{
    dinara::mode3::DagNodeInfo info;
    info.anchorChain.push_back(markerId);
    info.lengthBp = lengthBp;
    info.coverage = coverage;
    info.removed = false;
    return dag.addNode(info);
}

bool dagHasOutEdge(
    const dinara::mode3::DirectedAnchorGraph& dag,
    dinara::mode3::DagNodeId from,
    dinara::mode3::DagNodeId to)
{
    const auto& out = dag.getOutEdges(from);
    return std::find(out.begin(), out.end(), to) != out.end();
}

uint64_t findDagSegmentByAnchorChain(
    const dinara::mode3::DirectedAnchorGraph& dag,
    const std::vector<dinara::mode3::DagNodeId>& chain)
{
    for(uint64_t segId : dag.getActiveNodeIds()) {
        if(dag.getNode(segId).anchorChain == chain) {
            return segId;
        }
    }
    return std::numeric_limits<uint64_t>::max();
}

} // namespace

TEST_CASE("DirectedAnchorGraph unitigifyAll collapses MBG-style linear chains and rewrites paths",
    "[integration][dag][unitigify][mbg]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segR = addDagSimpleNode(dag, 90, 7, 5.0);
    const uint64_t segS = addDagSimpleNode(dag, 91, 7, 5.0);
    const uint64_t segP = addDagSimpleNode(dag, 100, 8, 10.0);
    const uint64_t segA = addDagSimpleNode(dag, 101, 10, 20.0);
    const uint64_t segB = addDagSimpleNode(dag, 102, 12, 30.0);
    const uint64_t segC = addDagSimpleNode(dag, 103, 14, 40.0);
    const uint64_t segQ = addDagSimpleNode(dag, 104, 9, 50.0);
    const uint64_t segT = addDagSimpleNode(dag, 110, 7, 5.0);
    const uint64_t segU = addDagSimpleNode(dag, 111, 7, 5.0);

    dag.addEdge(fwdNodeId(segR), fwdNodeId(segP), 1);
    dag.addEdge(fwdNodeId(segS), fwdNodeId(segP), 1);
    dag.addEdge(fwdNodeId(segP), fwdNodeId(segA), 2);
    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 3);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segC), 4);
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segQ), 1);
    dag.addEdge(fwdNodeId(segQ), fwdNodeId(segT), 1);
    dag.addEdge(fwdNodeId(segQ), fwdNodeId(segU), 1);

    const uint64_t p0 = dag.addPath({
        fwdNodeId(segP),
        fwdNodeId(segA),
        fwdNodeId(segB),
        fwdNodeId(segC),
        fwdNodeId(segQ)
    }, 1);
    const uint64_t p1 = dag.addPath({
        fwdNodeId(segB),
        fwdNodeId(segC),
        fwdNodeId(segQ)
    }, 1);
    const uint64_t p2 = dag.addPath({
        fwdNodeId(segP),
        fwdNodeId(segA),
        fwdNodeId(segB)
    }, 1);

    REQUIRE(dag.nodeCount() == 9);
    REQUIRE(dag.pathCount() == 3);

    dag.unitigifyAll();

    REQUIRE(dag.nodeExists(segP) == false);
    REQUIRE(dag.nodeExists(segA) == false);
    REQUIRE(dag.nodeExists(segB) == false);
    REQUIRE(dag.nodeExists(segC) == false);
    REQUIRE(dag.nodeExists(segQ) == false);
    REQUIRE(dag.nodeExists(segR));
    REQUIRE(dag.nodeExists(segS));
    REQUIRE(dag.nodeExists(segT));
    REQUIRE(dag.nodeExists(segU));
    REQUIRE(dag.nodeCount() == 5);
    REQUIRE(dag.pathCount() == 3);

    uint64_t mergedSeg = std::numeric_limits<uint64_t>::max();
    for(uint64_t segId : dag.getActiveNodeIds()) {
        const auto& node = dag.getNode(segId);
        if(node.anchorChain == std::vector<DagNodeId>{100, 101, 102, 103, 104}) {
            mergedSeg = segId;
            break;
        }
    }
    REQUIRE(mergedSeg != std::numeric_limits<uint64_t>::max());

    const auto& mergedInfo = dag.getNode(mergedSeg);
    REQUIRE(mergedInfo.lengthBp == (8 + (10 - 2) + (12 - 3) + (14 - 4) + (9 - 1)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(segR), fwdNodeId(mergedSeg)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(segS), fwdNodeId(mergedSeg)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedSeg), fwdNodeId(segT)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedSeg), fwdNodeId(segU)));

    REQUIRE(dag.getPath(p0) ==
        std::vector<DagNodeId>{fwdNodeId(mergedSeg)});
    REQUIRE(dag.getPath(p1) ==
        std::vector<DagNodeId>{fwdNodeId(mergedSeg)});
    REQUIRE(dag.getPath(p2) ==
        std::vector<DagNodeId>{fwdNodeId(mergedSeg)});
}

TEST_CASE("DirectedAnchorGraph unitigifyAll handles circular unitigs in MBG style",
    "[integration][dag][unitigify][mbg][circular]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segA = addDagSimpleNode(dag, 201, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 202, 11, 11.0);
    const uint64_t segC = addDagSimpleNode(dag, 203, 12, 12.0);

    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 2);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segC), 2);
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segA), 2);

    const uint64_t p0 = dag.addPath({
        fwdNodeId(segA),
        fwdNodeId(segB),
        fwdNodeId(segC),
        fwdNodeId(segA)
    }, 1);

    REQUIRE(dag.nodeCount() == 3);
    dag.unitigifyAll();
    REQUIRE(dag.nodeCount() == 1);
    REQUIRE(dag.pathCount() == 1);

    const auto active = dag.getActiveNodeIds();
    REQUIRE(active.size() == 1);
    const uint64_t mergedSeg = active.front();

    const auto& rewritten = dag.getPath(p0);
    REQUIRE(rewritten.empty() == false);
    for(DagNodeId n : rewritten) {
        REQUIRE(segmentOf(n) == mergedSeg);
    }
}

TEST_CASE("DirectedAnchorGraph unitigifyAll leaves non-unitigifiable branch nodes unchanged",
    "[integration][dag][unitigify][mbg][no-merge]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segA = addDagSimpleNode(dag, 401, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 402, 11, 11.0);
    const uint64_t segC = addDagSimpleNode(dag, 403, 12, 12.0);
    const uint64_t segD = addDagSimpleNode(dag, 404, 13, 13.0);

    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 2);
    dag.addEdge(fwdNodeId(segA), fwdNodeId(segC), 2);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segD), 2);
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segD), 2);

    const uint64_t p0 = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segB), fwdNodeId(segD)
    }, 1);
    const uint64_t p1 = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segC), fwdNodeId(segD)
    }, 1);

    dag.unitigifyAll();

    REQUIRE(dag.nodeCount() == 4);
    REQUIRE(dag.nodeExists(segA));
    REQUIRE(dag.nodeExists(segB));
    REQUIRE(dag.nodeExists(segC));
    REQUIRE(dag.nodeExists(segD));

    REQUIRE(dag.getPath(p0) == std::vector<DagNodeId>{
        fwdNodeId(segA), fwdNodeId(segB), fwdNodeId(segD)
    });
    REQUIRE(dag.getPath(p1) == std::vector<DagNodeId>{
        fwdNodeId(segA), fwdNodeId(segC), fwdNodeId(segD)
    });

    for(uint64_t segId : dag.getActiveNodeIds()) {
        REQUIRE(dag.getNode(segId).anchorChain.size() == 1);
    }
}

TEST_CASE("DirectedAnchorGraph unitigifyAll merges disjoint chains and is idempotent",
    "[integration][dag][unitigify][mbg][idempotent]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segA = addDagSimpleNode(dag, 501, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 502, 12, 12.0);
    const uint64_t segC = addDagSimpleNode(dag, 503, 14, 14.0);
    const uint64_t segD = addDagSimpleNode(dag, 504, 16, 16.0);
    const uint64_t segE = addDagSimpleNode(dag, 505, 9, 9.0);

    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 3);
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segD), 4);

    const uint64_t p0 = dag.addPath({fwdNodeId(segA), fwdNodeId(segB)}, 1);
    const uint64_t p1 = dag.addPath({fwdNodeId(segC), fwdNodeId(segD)}, 1);

    dag.unitigifyAll();

    REQUIRE(dag.nodeCount() == 3);
    REQUIRE(dag.nodeExists(segE));
    REQUIRE_FALSE(dag.nodeExists(segA));
    REQUIRE_FALSE(dag.nodeExists(segB));
    REQUIRE_FALSE(dag.nodeExists(segC));
    REQUIRE_FALSE(dag.nodeExists(segD));

    const uint64_t mergedAB =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{501, 502});
    const uint64_t mergedCD =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{503, 504});
    REQUIRE(mergedAB != std::numeric_limits<uint64_t>::max());
    REQUIRE(mergedCD != std::numeric_limits<uint64_t>::max());

    REQUIRE(dag.getPath(p0) == std::vector<DagNodeId>{fwdNodeId(mergedAB)});
    REQUIRE(dag.getPath(p1) == std::vector<DagNodeId>{fwdNodeId(mergedCD)});

    const uint64_t nodesAfterFirst = dag.nodeCount();
    const uint64_t edgesAfterFirst = dag.edgeCount();
    const auto path0AfterFirst = dag.getPath(p0);
    const auto path1AfterFirst = dag.getPath(p1);

    dag.unitigifyAll();

    REQUIRE(dag.nodeCount() == nodesAfterFirst);
    REQUIRE(dag.edgeCount() == edgesAfterFirst);
    REQUIRE(dag.getPath(p0) == path0AfterFirst);
    REQUIRE(dag.getPath(p1) == path1AfterFirst);
}

TEST_CASE("DirectedAnchorGraph unitigifyAll preserves MBG endpoint rewiring corner cases",
    "[integration][dag][unitigify][mbg][rewire]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segA = addDagSimpleNode(dag, 601, 20, 8.0);
    const uint64_t segB = addDagSimpleNode(dag, 602, 21, 8.0);
    const uint64_t segC = addDagSimpleNode(dag, 603, 22, 8.0);
    const uint64_t segX = addDagSimpleNode(dag, 604, 10, 5.0);
    const uint64_t segY = addDagSimpleNode(dag, 605, 10, 5.0);

    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 3);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segC), 4);
    dag.addEdge(revNodeId(segA), revNodeId(segC), 5);  // rc(first)->rc(last)
    dag.addEdge(revNodeId(segA), fwdNodeId(segA), 2);  // rc(first)->first
    dag.addEdge(revNodeId(segA), fwdNodeId(segX), 1);  // rc(first)->outside
    dag.addEdge(fwdNodeId(segC), revNodeId(segC), 2);  // last->rc(last)
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segY), 1);  // last->outside

    const uint64_t pFwd = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segB), fwdNodeId(segC)
    }, 1);
    const uint64_t pRev = dag.addPath({
        revNodeId(segC), revNodeId(segB), revNodeId(segA)
    }, 1);

    dag.unitigifyAll();

    const uint64_t mergedABC =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{601, 602, 603});
    REQUIRE(mergedABC != std::numeric_limits<uint64_t>::max());

    REQUIRE(dag.getPath(pFwd) == std::vector<DagNodeId>{fwdNodeId(mergedABC)});
    REQUIRE(dag.getPath(pRev) == std::vector<DagNodeId>{revNodeId(mergedABC)});

    REQUIRE(dagHasOutEdge(dag, revNodeId(mergedABC), revNodeId(mergedABC)));
    REQUIRE(dagHasOutEdge(dag, revNodeId(mergedABC), fwdNodeId(mergedABC)));
    REQUIRE(dagHasOutEdge(dag, revNodeId(mergedABC), fwdNodeId(segX)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedABC), revNodeId(mergedABC)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedABC), fwdNodeId(segY)));
}

TEST_CASE("DirectedAnchorGraph unitigifyAll rewrites internal and reverse-entry path orientations",
    "[integration][dag][unitigify][mbg][path-rewrite]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segA = addDagSimpleNode(dag, 701, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 702, 11, 11.0);
    const uint64_t segC = addDagSimpleNode(dag, 703, 12, 12.0);
    const uint64_t segZ = addDagSimpleNode(dag, 704, 8, 6.0);
    const uint64_t segW = addDagSimpleNode(dag, 705, 8, 6.0);

    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 2);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segC), 2);
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segW), 1);
    dag.addEdge(fwdNodeId(segZ), revNodeId(segC), 1);

    const uint64_t p0 = dag.addPath({fwdNodeId(segB), fwdNodeId(segC)}, 1);
    const uint64_t p1 = dag.addPath({revNodeId(segB), revNodeId(segA)}, 1);
    const uint64_t p2 = dag.addPath({
        fwdNodeId(segZ), revNodeId(segC), revNodeId(segB)
    }, 1);
    const uint64_t p3 = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segB), fwdNodeId(segC)
    }, 1);

    dag.unitigifyAll();

    const uint64_t mergedABC =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{701, 702, 703});
    REQUIRE(mergedABC != std::numeric_limits<uint64_t>::max());

    REQUIRE(dag.getPath(p0) == std::vector<DagNodeId>{fwdNodeId(mergedABC)});
    REQUIRE(dag.getPath(p1) == std::vector<DagNodeId>{revNodeId(mergedABC)});
    REQUIRE(dag.getPath(p2) == std::vector<DagNodeId>{
        fwdNodeId(segZ), revNodeId(mergedABC)
    });
    REQUIRE(dag.getPath(p3) == std::vector<DagNodeId>{fwdNodeId(mergedABC)});

    const auto& crossing = dag.getPathsCrossingNode(mergedABC);
    auto hasPath = [&](uint64_t pathIdx) {
        for(const auto& occ : crossing) {
            if(occ.pathIdx == pathIdx) return true;
        }
        return false;
    };
    REQUIRE(hasPath(p0));
    REQUIRE(hasPath(p1));
    REQUIRE(hasPath(p2));
    REQUIRE(hasPath(p3));
}

// Regression test: batch unitigifyAll must preserve edges between
// adjacent chains separated by a hub. Two chains [A,B] and [C,D]
// connected through hub H (H has degree > 1 in both directions)
// should produce merged nodes M_AB and M_CD with edges M_AB→H and H→M_CD.
TEST_CASE("DirectedAnchorGraph unitigifyAll preserves cross-chain edges",
    "[integration][dag][unitigify][mbg][cross-chain]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    // Hub P has out-degree 2 → can't be part of any chain.
    const uint64_t segP = addDagSimpleNode(dag, 800, 10, 5.0);
    const uint64_t segX = addDagSimpleNode(dag, 806, 8, 4.0);
    // Chain 1: A→B
    const uint64_t segA = addDagSimpleNode(dag, 801, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 802, 12, 12.0);
    // Hub H between chains.
    // H needs in-degree >= 2 (so B can't extend into H)
    // and out-degree >= 2 (so C can't extend backward into H).
    const uint64_t segH = addDagSimpleNode(dag, 808, 10, 6.0);
    const uint64_t segZ = addDagSimpleNode(dag, 809, 8, 4.0);  // H→Z
    const uint64_t segW = addDagSimpleNode(dag, 810, 8, 4.0);  // W→H
    // Chain 2: C→D
    const uint64_t segC = addDagSimpleNode(dag, 803, 14, 14.0);
    const uint64_t segD = addDagSimpleNode(dag, 804, 11, 11.0);
    // Hub Q has in-degree 2.
    const uint64_t segQ = addDagSimpleNode(dag, 805, 10, 5.0);
    const uint64_t segY = addDagSimpleNode(dag, 807, 8, 4.0);

    dag.addEdge(fwdNodeId(segP), fwdNodeId(segA), 2);
    dag.addEdge(fwdNodeId(segP), fwdNodeId(segX), 1);  // P out-degree 2
    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 3);   // chain 1 internal
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segH), 2);   // chain 1 → hub
    dag.addEdge(fwdNodeId(segW), fwdNodeId(segH), 1);   // H in-degree 2
    dag.addEdge(fwdNodeId(segH), fwdNodeId(segC), 2);   // hub → chain 2
    dag.addEdge(fwdNodeId(segH), fwdNodeId(segZ), 1);   // H out-degree 2
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segD), 3);   // chain 2 internal
    dag.addEdge(fwdNodeId(segD), fwdNodeId(segQ), 2);
    dag.addEdge(fwdNodeId(segY), fwdNodeId(segQ), 1);   // Q in-degree 2

    // Path crossing both chains through hub.
    const uint64_t p0 = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segB),
        fwdNodeId(segH),
        fwdNodeId(segC), fwdNodeId(segD)
    }, 1);

    dag.unitigifyAll();

    // Find merged nodes by anchor chain.
    const uint64_t mergedAB =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{801, 802});
    const uint64_t mergedCD =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{803, 804});
    REQUIRE(mergedAB != std::numeric_limits<uint64_t>::max());
    REQUIRE(mergedCD != std::numeric_limits<uint64_t>::max());

    // Edges through hub must exist.
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedAB), fwdNodeId(segH)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(segH), fwdNodeId(mergedCD)));

    // External edges preserved.
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(segP), fwdNodeId(mergedAB)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedCD), fwdNodeId(segQ)));

    // Path must be rewritten to use both merged nodes.
    REQUIRE(dag.getPath(p0) == std::vector<DagNodeId>{
        fwdNodeId(mergedAB), fwdNodeId(segH), fwdNodeId(mergedCD)
    });
}

// Regression test: batch unitigifyAll must perform secondary merges.
// When chain [A,B] and chain [C,D] are adjacent (B→C) with both
// endpoints having degree 1 after merging, the merged nodes M1→M2
// should themselves be merged into a single node M3.
TEST_CASE("DirectedAnchorGraph unitigifyAll performs secondary merges",
    "[integration][dag][unitigify][mbg][secondary-merge]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    // Hub P with out-degree 2 prevents A from extending backward.
    const uint64_t segP = addDagSimpleNode(dag, 900, 10, 5.0);
    const uint64_t segX = addDagSimpleNode(dag, 906, 8, 4.0);
    // Chain 1: A→B (P has out-degree 2, so A can't extend backward past P)
    const uint64_t segA = addDagSimpleNode(dag, 901, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 902, 12, 12.0);
    // Chain 2: C→D (Q has in-degree 2, so D can't extend forward past Q)
    const uint64_t segC = addDagSimpleNode(dag, 903, 14, 14.0);
    const uint64_t segD = addDagSimpleNode(dag, 904, 11, 11.0);
    // Hub Q with in-degree 2 prevents D from extending forward.
    const uint64_t segQ = addDagSimpleNode(dag, 905, 10, 5.0);
    const uint64_t segY = addDagSimpleNode(dag, 907, 8, 4.0);

    dag.addEdge(fwdNodeId(segP), fwdNodeId(segA), 2);
    dag.addEdge(fwdNodeId(segP), fwdNodeId(segX), 1);  // P out-degree 2
    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 3);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segC), 2);   // cross-chain
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segD), 3);
    dag.addEdge(fwdNodeId(segD), fwdNodeId(segQ), 2);
    dag.addEdge(fwdNodeId(segY), fwdNodeId(segQ), 1);   // Q in-degree 2

    // Path spanning all four chain members.
    const uint64_t p0 = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segB),
        fwdNodeId(segC), fwdNodeId(segD)
    }, 1);

    dag.unitigifyAll();

    // After batch merge: M_AB and M_CD both have degree 1 connecting
    // to each other. The secondary pass should merge them into M_ABCD.
    const uint64_t mergedABCD =
        findDagSegmentByAnchorChain(dag,
            std::vector<DagNodeId>{901, 902, 903, 904});
    REQUIRE(mergedABCD != std::numeric_limits<uint64_t>::max());

    // The path should reference the fully-merged node.
    REQUIRE(dag.getPath(p0) == std::vector<DagNodeId>{
        fwdNodeId(mergedABCD)
    });

    // External edges must reach the final merged node.
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(segP), fwdNodeId(mergedABCD)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedABCD), fwdNodeId(segQ)));
}


// =============================================================================
// Test: .fastq.gz loading produces the same reads as .fastq
// =============================================================================
TEST_CASE("ReadLoader handles gzip-compressed fastq files", "[readloader][gz]") {
    // Create a small fastq in memory, write both plain and gzipped versions.
    auto dir = makeUniqueTempDir("dinara_gz_test_");
    auto plainPath = dir / "reads.fastq";
    auto gzPath = dir / "reads.fastq.gz";

    // Write plain fastq.
    {
        std::ofstream out(plainPath);
        for(int i = 0; i < 100; i++) {
            std::string seq(200, 'A');
            // Make each read slightly different.
            for(int j = 0; j < 200; j++) {
                seq[j] = "ACGT"[(i * 7 + j * 3) % 4];
            }
            out << "@read_" << i << "\n" << seq << "\n+\n" << std::string(200, '~') << "\n";
        }
    }

    // Create gzipped version using system gzip.
    {
        std::string cmd = "gzip -k " + plainPath.string();
        REQUIRE(std::system(cmd.c_str()) == 0);
        REQUIRE(fs::exists(gzPath));
    }

    // Load from plain fastq.
    uint64_t countPlain = 0;
    {
        auto subdir = dir / "plain";
        fs::create_directories(subdir);
        std::string prefix = subdir.string() + "/";
        withSilencedIoInDir(subdir, [&] {
            Assembler assembler(prefix, true, 0, 4096);
            assembler.addReads(plainPath.string(), 0, false, 1);
            countPlain = assembler.getReads().readCount();
        });
    }

    // Load from gzipped fastq.
    uint64_t countGz = 0;
    {
        auto subdir = dir / "gz";
        fs::create_directories(subdir);
        std::string prefix = subdir.string() + "/";
        withSilencedIoInDir(subdir, [&] {
            Assembler assembler(prefix, true, 0, 4096);
            assembler.addReads(gzPath.string(), 0, false, 1);
            countGz = assembler.getReads().readCount();
        });
    }

    REQUIRE(countPlain == 100);
    REQUIRE(countGz == 100);
    REQUIRE(countPlain == countGz);

    fs::remove_all(dir);
}


TEST_CASE("ReadLoader: GIAB fastq.gz matches plain fastq", "[readloader][gz][giab]") {
    // Use the real GIAB test file if available.
    fs::path plainPath = fs::path(__FILE__).parent_path() / "GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq";
    if(!fs::exists(plainPath)) {
        WARN("Skipping: GIAB test file not found");
        return;
    }

    auto dir = makeUniqueTempDir("dinara_gz_giab_");

    // Create gzipped version.
    fs::path gzPath = dir / "reads.fastq.gz";
    {
        std::string cmd = "gzip -c " + plainPath.string() + " > " + gzPath.string();
        REQUIRE(std::system(cmd.c_str()) == 0);
        REQUIRE(fs::exists(gzPath));
    }

    uint64_t countPlain = 0;
    {
        auto subdir = dir / "plain";
        fs::create_directories(subdir);
        std::string prefix = subdir.string() + "/";
        withSilencedIoInDir(subdir, [&] {
            Assembler assembler(prefix, true, 0, 4096);
            assembler.addReads(plainPath.string(), 0, false, 1);
            countPlain = assembler.getReads().readCount();
        });
    }

    uint64_t countGz = 0;
    {
        auto subdir = dir / "gz";
        fs::create_directories(subdir);
        std::string prefix = subdir.string() + "/";
        withSilencedIoInDir(subdir, [&] {
            Assembler assembler(prefix, true, 0, 4096);
            assembler.addReads(gzPath.string(), 0, false, 1);
            countGz = assembler.getReads().readCount();
        });
    }

    REQUIRE(countPlain > 0);
    REQUIRE(countPlain == countGz);
    INFO("GIAB: " << countPlain << " reads from both plain and gz");

    fs::remove_all(dir);
}
