#ifndef DINARA_MODE3_DIRECTED_ANCHOR_GRAPH_HPP
#define DINARA_MODE3_DIRECTED_ANCHOR_GRAPH_HPP

/******************************************************************************

DirectedAnchorGraph — Verkko-style directed graph resolution for dinara.

Node representation (matching Verkko exactly):
  Each segment (bidirected anchor / merged unitig) has TWO oriented node IDs:
    - Forward:  segmentId * 2     (even, displayed as ">segId")
    - Reverse:  segmentId * 2 + 1 (odd,  displayed as "<segId")
  RC operation is a single XOR:  rcNode(n) = n ^ 1

  Edges, paths, and all graph operations use oriented DagNodeIds (plain uint64_t).
  Segment metadata (anchor chain, coverage) is stored once per segment.

Key Verkko concepts → dinara types:
  ">X" / "<X"         →  fwdNodeId(X) / revNodeId(X)
  revnode(">X")="<X"  →  rcNode(n) = n ^ 1
  canon(e)             →  canonEdge(from, to)
  edges[">X"]          →  getOutEdges(fwdNodeId(X))
  paths_crossing[X]    →  pathsCrossing index (by segment ID)
  node_seqs[X]         →  DagNodeInfo::anchorChain (per segment, forward orientation)

******************************************************************************/

// Forward declaration — full definition in mode3-Anchor.hpp.
namespace dinara { namespace mode3 { class Anchors; class DirectedAnchors; } }

#include "Base.hpp"
#include "LongBaseSequence.hpp"
#include "MultithreadedObject.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <tuple>

namespace dinara {
namespace mode3 {


// ============================================================================
// DagNodeId — oriented node identifier.
//
// Even IDs = forward (">"), odd IDs = reverse complement ("<").
// Segment ID = DagNodeId >> 1.
// No struct, no custom hash — just uint64_t with XOR for RC.
// ============================================================================

using DagNodeId = uint64_t;


// ---- Oriented node helpers (Verkko-style) ----

// Reverse complement: flip the lowest bit.
inline DagNodeId rcNode(DagNodeId n) { return n ^ 1; }

// Extract the canonical segment ID.
inline uint64_t segmentOf(DagNodeId n) { return n >> 1; }

// Check forward orientation.
inline bool isForward(DagNodeId n) { return !(n & 1); }

// Create oriented node IDs from a segment ID.
inline DagNodeId fwdNodeId(uint64_t segId) { return segId << 1; }
inline DagNodeId revNodeId(uint64_t segId) { return (segId << 1) | 1; }


// ---- Tip sentinel (no predecessor or successor in a path) ----

constexpr DagNodeId dagTipSentinel = UINT64_MAX;


// ---- Canonical edge ----

inline std::pair<DagNodeId, DagNodeId>
canonEdge(DagNodeId from, DagNodeId to) {
    DagNodeId rTo = to ^ 1;
    DagNodeId rFrom = from ^ 1;
    if(rTo < from || (rTo == from && rFrom < to)) {
        return {rTo, rFrom};
    }
    return {from, to};
}


// ---- Reverse-complement a path ----

inline std::vector<DagNodeId>
rcPath(const std::vector<DagNodeId>& path) {
    std::vector<DagNodeId> result;
    result.reserve(path.size());
    for(auto it = path.rbegin(); it != path.rend(); ++it) {
        result.push_back(*it ^ 1);
    }
    return result;
}


// ============================================================================
// DagTriplet — (predecessor, through, successor) for resolution.
// ============================================================================

struct DagTriplet {
    DagNodeId from;     // dagTipSentinel if path starts here
    DagNodeId through;  // forward-oriented ID of the through-segment
    DagNodeId to;       // dagTipSentinel if path ends here
    uint64_t support = 0;

    bool isFromTip() const { return from == dagTipSentinel; }
    bool isToTip() const { return to == dagTipSentinel; }
};


// ============================================================================
// DagNodeInfo — per-segment metadata.
// Stored once per canonical segment (not per oriented node).
// ============================================================================

struct DagNodeInfo {
    // Chain of base anchor IDs (in forward orientation of this segment).
    // Each element is an AnchorId (= DagNodeId = uint64_t); even=forward, odd=reverse.
    // Initially single-element; grows after unitigification.
    std::vector<DagNodeId> anchorChain;

    // Estimated length in base pairs.
    uint64_t lengthBp = 0;

    // Average read coverage.
    double coverage = 0.0;

    // Whether this segment has been removed.
    bool removed = false;
};


// ============================================================================
// DagPathsCrossingIndex — maps segment ID → set of path indices.
// ============================================================================

class DagPathsCrossingIndex {
public:
    void clear() { index.clear(); }

    void addPath(uint64_t pathIdx, const std::vector<DagNodeId>& path);
    void removePath(uint64_t pathIdx, const std::vector<DagNodeId>& path);
    void rebuild(const std::vector<std::vector<DagNodeId>>& paths,
                 const std::vector<bool>& pathRemoved);

    const std::unordered_set<uint64_t>& getPathsCrossing(uint64_t segId) const;
    bool hasPathsCrossing(uint64_t segId) const;

private:
    std::unordered_map<uint64_t, std::unordered_set<uint64_t>> index;
    static const std::unordered_set<uint64_t> emptySet;
};


// ============================================================================
// DirectedAnchorGraphConfig
// ============================================================================

struct DirectedAnchorGraphConfig {
    std::vector<uint64_t> resolveSteps = {20, 10, 5};
    double minAllowedCoverage = 5.0;
    uint64_t maxResolveLength = 500000;
    uint64_t maxUnconditionalResolveLength = 0;
    uint64_t maxLocalResolve = 0;
    uint64_t minEdgeSupport = 3;
    std::string outputGfa = "DirectedAnchorGraph.gfa";
    std::string outputPaths = "DirectedAnchorGraph.paths.gaf";

    // MBG-style cleaning parameters
    bool doCleaning = true;
    bool doGuesswork = false;
    bool copycountFilterHeuristic = false;
    bool resolvePalindromesGlobal = false;
    double tipMinCoverage = 3.0;
    uint64_t tipMinLength = 10;
    uint64_t tipMaxLength = 10000;
    double tipSecondPassMinCoverage = 2.0;
    uint64_t tipSecondPassMinLength = 5;
    double crosslinkMinCoverage = 2.0;
    uint64_t crosslinkMinLength = 10;
    double crosslinkFirstPassMinCoverage = 1.0;
    uint64_t crosslinkFirstPassMinLength = 5;

    // MBG resolveUnitigs does a second round with minCoverage=1.
    // In this DAG pipeline, we emulate that with minEdgeSupport=1.
    bool doFinalLowSupportResolvePass = true;
};

struct DagPathReadSupport {
    uint64_t readId = 0;
    uint64_t readNameIndex = 0;
    uint64_t readInfoIndex = 0;
    uint64_t readPosZeroOffset = 0;
    uint64_t leftClip = 0;
    uint64_t rightClip = 0;
    uint64_t readPosStartIndex = 0;
    uint64_t readPosEndIndex = 0;
    // Per-path-node read positions (relative to readPosZeroOffset).
    // Size is path.size()+1 for the corresponding path: node starts + final end.
    std::vector<uint64_t> readPoses;
};

struct DagReadInterval {
    uint64_t start = 0;
    uint64_t end = 0;
    bool valid = false;
};


// ============================================================================
// DirectedAnchorGraph — Verkko-style directed resolution graph.
//
// Each segment has two oriented node IDs (even=fwd, odd=rev).
// Edges and paths use oriented DagNodeIds.
// Segment metadata is stored once per canonical segment.
// ============================================================================

class DirectedAnchorGraph : public MultithreadedObject<DirectedAnchorGraph> {
public:

    DirectedAnchorGraph() : MultithreadedObject(*this) {}

    // Build from mode3::Anchors (multithreaded).
    // Each anchor pair (2*m, 2*m+1) maps directly to a DAG segment m
    // with fwdNodeId = 2*m, revNodeId = 2*m+1.
    // Build from mode3::Anchors (multithreaded).
    // Each anchor pair (2*m, 2*m+1) maps directly to a DAG segment m.
    void buildFromAnchors(
        const Anchors& anchors,
        uint64_t threadCount = 1,
        const std::vector<DagReadInterval>* readIntervals = nullptr);

    // Optimized build from DirectedAnchors.
    // This overload leverages precomputed base positions for O(1) overlap
    // calculation and avoids the costly search for anchor coordinates.
    void buildFromAnchors(
        const DirectedAnchors& anchors,
        uint64_t threadCount = 1);


    // ---- Segment accessors ----

    uint64_t nodeCount() const;          // active segments
    uint64_t totalNodeCount() const;     // total including removed
    uint64_t edgeCount() const;          // directed edges (each canonical edge counted twice)
    uint64_t pathCount() const;          // active paths

    bool nodeExists(uint64_t segId) const;
    const DagNodeInfo& getNode(uint64_t segId) const { return nodes[segId]; }

    // Collect active segment IDs (sorted).
    std::vector<uint64_t> getActiveNodeIds() const;


    // ---- Oriented node operations ----

    const std::vector<DagNodeId>& getOutEdges(DagNodeId from) const;
    uint64_t outDegree(DagNodeId from) const;
    uint64_t inDegree(DagNodeId from) const;


    // ---- Edge manipulation (maintains RC invariant) ----
    // addEdge(from, to) also inserts rcNode(to) → rcNode(from).

    void addEdge(
        DagNodeId from,
        DagNodeId to,
        uint64_t overlapBp = 0);
    void removeEdge(DagNodeId from, DagNodeId to);
    void removeAllEdges(uint64_t segId);
    uint64_t getBpOverlap(DagNodeId from, DagNodeId to) const;
    void setBpOverlap(DagNodeId from, DagNodeId to, uint64_t overlapBp);
    uint64_t nodeBpLength(uint64_t segId) const;
    uint64_t pathBpLength(const std::vector<DagNodeId>& path) const;


    // ---- Segment manipulation ----
    // addNode returns the segment ID.
    // Forward oriented ID = segId << 1, reverse = (segId << 1) | 1.

    uint64_t addNode(const DagNodeInfo& info);
    void removeNode(uint64_t segId);


    // ---- Path manipulation ----

    uint64_t addPath(
        std::vector<DagNodeId> path,
        uint64_t weight = 1,
        std::vector<DagPathReadSupport> readSupport = {});
    void removePath(uint64_t pathIdx);
    const std::vector<DagNodeId>& getPath(uint64_t pathIdx) const;
    void setPath(uint64_t pathIdx, std::vector<DagNodeId>&& newPath);


    // ---- Graph cleaning (MBG-style) ----

    struct CleaningStats {
        uint64_t nodesRemoved = 0;
        uint64_t edgesRemoved = 0;
        std::unordered_set<uint64_t> maybeUnitigifiable;
    };

    struct ResolveStats {
        uint64_t nodesResolved = 0;
        uint64_t nodesAdded = 0;
        std::unordered_set<uint64_t> newSegs;
        std::unordered_set<uint64_t> maybeUnitigifiable;
        std::unordered_map<DagNodeId, uint64_t> maybeTrimmable;
    };

    CleaningStats removeLowCoverageTips(
        double maxRemovableCoverage,
        double minSafeCoverage,
        uint64_t maxRemovableLength,
        const std::unordered_set<uint64_t>& maybeUntippable = {});

    CleaningStats removeLowCoverageCrosslinks(
        double minCoverage,
        uint64_t minLength);

    CleaningStats cleanComponentsByCopynumber(
        double averageCoverage);

    CleaningStats cleanComponentsByCopynumber(
        double averageCoverage,
        uint64_t minLongLength,
        uint64_t minUnresolvableLength,
        uint64_t maxUnresolvableLength,
        const std::unordered_set<uint64_t>& checkThese,
        uint64_t minNew);

    CleaningStats removeLowCoverageSegments(
        double minCoverage);


    // ---- Resolution (Verkko's resolve_triplets_kmerify.py) ----

    void runResolution(const DirectedAnchorGraphConfig& config);

    std::vector<DagTriplet> getValidTriplets(
        uint64_t segId,
        uint64_t minEdgeSupport) const;

    ResolveStats resolveNodes(
        const std::vector<uint64_t>& candidates,
        uint64_t minEdgeSupport,
        bool unconditional = false,
        bool copycountFilterHeuristic = false,
        bool guesswork = false,
        double averageCoverage = 1.0);

    std::unordered_set<uint64_t> resolveHairpins(
        const std::vector<uint64_t>& candidates);

    void resolveRound(
        uint64_t minEdgeSupport,
        uint64_t maxResolveLength,
        bool doCleaning = false,
        bool doGuesswork = false,
        uint64_t maxUnconditionalResolveLength = 0,
        bool copycountFilterHeuristic = false,
        uint64_t maxLocalResolve = 0,
        bool resolvePalindromesGlobal = false);


    // ---- Unitigification ----

    uint64_t unitigifyOne(uint64_t segId);
    void unitigifyAll();
    std::vector<DagNodeId> extendForward(DagNodeId start) const;


    // ---- Path maintenance ----

    void splitPathsAtBreaks();

    void replacePathNodes(
        const std::unordered_map<uint64_t,
            std::unordered_map<DagNodeId,
                std::unordered_map<DagNodeId, uint64_t>>>& resolutionMap);

    // MBG-style path-derived coverage helpers.
    double getPathCoverage(uint64_t segId) const;
    uint64_t getEdgePathCoverage(
        DagNodeId from,
        DagNodeId to,
        uint64_t maxCount = UINT64_MAX) const;


    // ---- GFA I/O and Consistency ----

    std::vector<Base> getSegmentSequence(uint64_t segId) const;
    void verifyEdgeConsistency() const;

    void writeGfa(const std::string& filename, bool omitSequences = false) const;
    void writePaths(const std::string& filename) const;
    void writeSummary(std::ostream& out) const;


    // ---- HTTP server accessors ----

    const std::unordered_set<uint64_t>& getPathsCrossingNode(uint64_t segId) const {
        return pathsCrossing.getPathsCrossing(segId);
    }

    bool pathExists(uint64_t pathIdx) const {
        return pathIdx < paths.size() && !pathRemoved[pathIdx];
    }

    uint64_t totalPathCount() const { return paths.size(); }

private:
    struct DagEdgePairHash {
        size_t operator()(const std::pair<DagNodeId, DagNodeId>& p) const
        {
            return std::hash<uint64_t>()(p.first) ^
                (std::hash<uint64_t>()(p.second) << 1);
        }
    };

    // Segment storage — indexed by segment ID (= DagNodeId >> 1).
    // Verkko-style: flat vector for O(1) access, cache-friendly iteration.
    std::vector<DagNodeInfo> nodes;

    // Directed adjacency — indexed by oriented DagNodeId.
    // Verkko-style: flat vector of adjacency lists for O(1) access.
    // Size = 2 * nodes.size() (even IDs = forward, odd IDs = reverse).
    std::vector<std::vector<DagNodeId>> edges;
    // Canonical edge overlap in bp (same value for both orientations).
    std::unordered_map<
        std::pair<DagNodeId, DagNodeId>,
        uint64_t,
        DagEdgePairHash> edgeOverlaps;

    // Paths — sequences of oriented DagNodeIds.
    std::vector<std::vector<DagNodeId>> paths;
    std::vector<bool> pathRemoved;
    // Path-group size (MBG PathGroup-like read multiplicity).
    std::vector<uint64_t> pathWeights;
    // MBG PathGroup-like per-read metadata for each path group.
    std::vector<std::vector<DagPathReadSupport>> pathReadIds;

    // Path-level read support stores MBG-like clipping and positional metadata.

    // Paths-crossing index (by segment ID).
    DagPathsCrossingIndex pathsCrossing;

    const Anchors* anchorsPtr = nullptr;

    // MBG-style incremental tip candidate bookkeeping.
    std::vector<uint64_t> everTippable;
    uint64_t lastTippableChecked = 0;

    // Merge a linear chain into a single segment.
    uint64_t replaceUnitig(const std::vector<DagNodeId>& chain);

    // Empty vector for missing edge lookups.
    static const std::vector<DagNodeId> emptyEdgeVec;

    // ---- Multithreaded buildFromAnchors internals ----

    // Per-thread collected data during parallel journey processing.
    std::vector<std::vector<std::vector<DagNodeId>>> threadPaths;
    std::vector<std::vector<DagPathReadSupport>> threadPathReadSupport;
    std::vector<std::vector<std::pair<DagNodeId, DagNodeId>>> threadEdgePairs;
    std::vector<std::vector<std::tuple<DagNodeId, DagNodeId, uint64_t>>> threadEdgeOverlapSamples;
    const Anchors* buildAnchorsPtr = nullptr;
    const DirectedAnchors* buildDirectedAnchorsPtr = nullptr;
    const std::vector<DagReadInterval>* buildReadIntervalsPtr = nullptr;
    std::vector<uint64_t> strand0Indices;  // even oridValues with journey.size() >= 2

    void buildFromAnchorsThreadFunction(size_t threadId);
    void buildFromDirectedAnchorsThreadFunction(size_t threadId);

    uint64_t pathGroupWeight(uint64_t pathIdx) const;
    uint64_t getCrossingCount(uint64_t segId) const;
    std::vector<DagTriplet> getRawTriplets(
        uint64_t segId,
        uint64_t minEdgeSupport,
        bool partTriplets) const;
    std::vector<DagTriplet> getReadSupportedTriplets(
        const std::unordered_set<uint64_t>& resolvables,
        uint64_t segId,
        uint64_t minEdgeSupport,
        bool unconditional,
        bool guesswork) const;
    std::vector<DagTriplet> getGuessworkTriplets(
        const std::unordered_set<uint64_t>& resolvables,
        uint64_t segId,
        uint64_t minEdgeSupport,
        bool unconditional,
        double averageCoverage) const;
    void filterCopyCountTriplets(
        uint64_t segId,
        std::vector<DagTriplet>& triplets,
        double averageCoverage) const;
    std::vector<DagTriplet> getValidTripletsForResolve(
        const std::unordered_set<uint64_t>& resolvables,
        uint64_t segId,
        uint64_t minEdgeSupport,
        bool unconditional,
        bool guesswork,
        bool copycountFilterHeuristic,
        double averageCoverage) const;
    bool nodeIsPalindrome(uint64_t segId) const;
    bool isLocallyRepetitive(uint64_t segId) const;
    std::unordered_set<uint64_t> filterToOnlyLocallyRepetitives(
        const std::unordered_set<uint64_t>& unfilteredResolvables,
        uint64_t maxDist,
        bool resolvePalindromesGlobal) const;
    std::unordered_set<uint64_t> trimNodes(
        const std::unordered_map<DagNodeId, uint64_t>& maybeTrimmable);
    void addPathButFirstMaybeTrim(
        std::vector<DagNodeId>&& path,
        std::vector<DagPathReadSupport>&& reads);
    void checkValidity() const;
    uint64_t getTrimAmountToCheck(
        DagNodeId from,
        DagNodeId to) const;
    uint64_t getAnchorSize(
        DagNodeId from,
        DagNodeId to) const;
    uint64_t createFakeEdgeNode(
        DagNodeId from,
        DagNodeId to,
        const std::unordered_set<uint64_t>& resolvables,
        const std::unordered_set<uint64_t>& unresolvables);
    uint64_t createEdgeNode(
        DagNodeId from,
        DagNodeId to,
        const std::unordered_set<uint64_t>& resolvables,
        const std::unordered_set<uint64_t>& unresolvables);
    void replacePathsFromEdgeNodes(
        const std::unordered_set<uint64_t>& actuallyResolvables,
        const std::unordered_map<std::pair<DagNodeId, DagNodeId>, uint64_t, DagEdgePairHash>& newEdgeNodes);
    void compactGraph();

    std::vector<uint64_t> computePathNodeBoundaries(
        const std::vector<DagNodeId>& path) const;
    bool reconcileReadSupportToPath(
        DagPathReadSupport& read,
        uint64_t pathBpLength) const;
    bool hasUsableReadPoses(
        const DagPathReadSupport& read,
        const std::vector<uint64_t>& pathBoundaries) const;
    uint64_t mapPathPositionFromSupport(
        const DagPathReadSupport& read,
        const std::vector<uint64_t>& pathBoundaries,
        uint64_t pathPos) const;
    std::optional<DagPathReadSupport> projectReadSupportToBoundaries(
        const DagPathReadSupport& read,
        const std::vector<uint64_t>& oldPathBoundaries,
        const std::vector<uint64_t>& projectedBoundaries) const;
    bool isReadSupportStrictForPath(
        const DagPathReadSupport& read,
        const std::vector<uint64_t>& pathBoundaries,
        uint64_t pathBpLength) const;
    void refreshReadPosesForPath(uint64_t pathIdx);
};


} // namespace mode3
} // namespace dinara

#endif
