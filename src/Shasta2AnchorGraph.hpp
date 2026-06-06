#pragma once

// Shasta2AnchorGraph.hpp

#include "AnchorWindows.hpp"
#include "BidirectedAnchorGraph.hpp"
#include "DinaraDetangle.hpp"
#include "Reads.hpp"
#include "Shasta2AnchorPair.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/base_object.hpp>

// Standard library.
#include <map>
#include <set>
#include "utility.hpp"
#include "vector.hpp"

namespace dinara {

        class Shasta2AnchorGraph;
        class Shasta2AnchorGraphEdge;
        using Shasta2AnchorGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            boost::no_property,
            Shasta2AnchorGraphEdge>;
}

class dinara::Shasta2AnchorGraphEdge {
public:
    Shasta2AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    uint64_t id = invalid<uint64_t>;
    bool useForAssembly = false;

    // Inter-window edge attributes.
    // supportingSpanPrev/Next: base-level span of the selected read
    // in the source (prev) and destination (next) window.
    // 0 for intra-window edges.
    uint64_t supportingSpanPrev = 0;
    uint64_t supportingSpanNext = 0;
    // sharedReadCount: number of distinct reads that touch both windows
    // (not just the reads using this specific anchor pair).
    uint64_t sharedReadCount = 0;

    Shasta2AnchorGraphEdge(const Shasta2AnchorPair& anchorPair, uint64_t offset, uint64_t id) :
        anchorPair(anchorPair),
        offset(offset),
        id(id)
    {}

    Shasta2AnchorGraphEdge() {}

    uint64_t coverage() const {return anchorPair.size();}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorPair;
        ar & offset;
        ar & id;
        ar & useForAssembly;
        ar & supportingSpanPrev;
        ar & supportingSpanNext;
        ar & sharedReadCount;
    }
};



class dinara::Shasta2AnchorGraph :
    public Shasta2AnchorGraphBaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<Shasta2AnchorGraph> {
public:
    using AnchorPairKey = std::pair<Shasta2AnchorId, Shasta2AnchorId>;

    // Construct the AnchorGraph from the Journeys using the same edge
    // creation rule as mode3::AnchorGraph: for each anchor, call
    // Shasta2Anchors::findChildren and create one edge per child that
    // satisfies minEdgeCoverage. The threadCount parameter is accepted
    // for API compatibility but is not used.
    Shasta2AnchorGraph(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        uint64_t minEdgeCoverage,
        uint64_t threadCount);

    // Construct from anchor windows: each window becomes a chain of its
    // backbone anchors, and inter-window edges are discovered by walking
    // read journeys.
    // If bypassEdges is provided, additional edges are created to bypass
    // detangled windows.
    Shasta2AnchorGraph(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        const vector<AnchorWindow>& anchorWindows,
        uint64_t minInterWindowCoverage,
        uint64_t minInterWindowEdgeCoverage,
        uint64_t threadCount,
        const Reads* reads = nullptr,
        const vector<DetangleBypassEdge>* bypassEdges = nullptr,
        const std::set<std::pair<uint32_t, uint32_t>>* detourWindowPairs = nullptr);

    // Default constructor (empty graph).
    Shasta2AnchorGraph() : MultithreadedObject<Shasta2AnchorGraph>(*this) {}

    // Constructor from binary data.
    Shasta2AnchorGraph(const MappedMemoryOwner&, const string& name);

    uint64_t nextEdgeId = 0;

    void transitiveReduction(
        uint64_t transitiveReductionMaxEdgeCoverage,
        uint64_t maxDistance);
    uint64_t cutWeakStalksLeadingToBranch(
        const Shasta2Anchors& anchors,
        uint64_t maxTipReadCount);
    // Disable an edge and its RC mirror (dst^1 -> src^1).
    void disableEdge(edge_descriptor e);

    // Find detour window pairs: (W, X) where window X enters W at
    // backbone position i and exits at position j > i. These pairs
    // are used during journey walks to suppress W→X→W transitions.
    std::set<std::pair<uint32_t, uint32_t>> findDetourWindowPairs(
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys) const;

    // Trim dangling backbone ends beyond outermost inter-window edges.
    // Returns the number of trimmed vertices.
    uint64_t trimBackbones(
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys);

    // Remove inter-window edges that land internally on a window's backbone
    // (between the first and last inter-window connection points).
    // For reads traversing internal connections, create bypass edges that
    // skip the internal windows. Returns the number of edges removed.
    uint64_t removeInternalConnections(
        const Shasta2Anchors& anchors,
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys);

    // Remove edges between a window and its RC counterpart.
    uint64_t removeRcWindowConnections();

    // Window-level transitive reduction: if A→B→C exists and A→C also exists,
    // remove A→C (the direct edge is redundant).
    uint64_t windowTransitiveReduction();

    // Dump detailed connection statistics for the largest window.
    void writeWindowConnectionStats(
        const Shasta2Anchors& anchors,
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys) const;

    // Disable all edges of windows that have no active inter-window edges.
    // Returns the number of isolated windows removed.
    uint64_t removeIsolatedWindows(
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys);

    // Remove short tip chains from the window graph.
    // A tip is a dead-end linear chain of windows (connected on only one side).
    // Chains of length <= maxTipWindows are removed (all their edges disabled).
    // This mirrors hifiasm's asg_arc_cut_tips: short tips are artifacts,
    // long tips are legitimate subregion ends.
    // Returns the number of windows removed.
    uint64_t removeTipWindows(
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys,
        uint32_t maxTipWindows = 3);

    // Pop superbubbles: find superbubbles in the active-edge subgraph,
    // choose the best path through each (prefer the source window's
    // backbone path), and disable edges on all other paths.
    // maxSize limits the number of internal vertices per superbubble.
    // Returns the number of superbubbles popped.
    uint64_t popSuperbubbles(
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys,
        uint64_t maxSize = 100);

private:
    bool transitiveReductionCanRemove(edge_descriptor, uint64_t transitiveReductionMaxDistance) const;
public:

    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<Shasta2AnchorGraphBaseClass>(*this);
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    void save(const string& name) const;
    void load(const string& name);

    // Write the graph to GFA format.
    // If anchorWindows is provided, each link line gets a tp:Z: tag
    // Convert to a Verkko-style bidirected graph. Each anchor pair
    // (anchorId / 2) becomes one node. RC mirror edges are collapsed
    // into single bidirected links. Only useForAssembly edges are included.
    // Normalizes orientations so backbone chains appear as +/+ links.
    BidirectedAnchorGraph toBidirected(
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys) const;

    // classifying it as "intra", "endpoint", or "internal".
    void writeGfa(const string& fileName,
                  const vector<AnchorWindow>* anchorWindows = nullptr) const;
    void writeBubbleFinderGraph(const string& fileName, bool useForAssemblyOnly = true) const;

    // Write Bandage color CSV: each anchor colored by its window.
    void writeCsv(const string& fileName) const;

    // Save binary data (dinara's own format).
    void saveAnchorGraph(const string& name = "Shasta2AnchorGraph") const { save(name); }

    // Export in shasta2-compatible MemoryMapped::Vector<char> format.
    // The output file can be passed to shasta2 via --external-anchor-graph-name.
    void saveForShasta2(const string& fileName) const;

    // Per-anchor window assignment (populated by the anchor-window constructor).
    // Maps anchorId -> windowId. noWindow means unmapped.
    // Includes RC mirror windows (windowId >= windowCount).
    static constexpr uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    vector<uint32_t> anchorToWindow;
    uint32_t windowCount = 0;

    // Per-window: distinct oriented read IDs that touch the window.
    std::map<uint32_t, std::set<uint32_t>> windowReads;

    // Per-read: ordered window sequence (consecutive duplicates removed).
    // Gives the order in which the read visits windows.
    std::map<uint32_t, vector<uint32_t>> readWindows;

};
