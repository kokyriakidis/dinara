#pragma once

// Shasta2AnchorGraph.hpp

#include "AnchorWindows.hpp"
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
    bool isEndpointAnchorPrev = false;  // source anchor is an endpoint anchor
    bool isEndpointAnchorNext = false;  // target anchor is an endpoint anchor

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
        uint64_t threadCount,
        const Reads* reads = nullptr,
        const vector<DetangleBypassEdge>* bypassEdges = nullptr);

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

    // Normalized endpoint window pairs from the two-pass edge creation.
    // Stored as (min, max) for direction-independent lookup.
    // A pair (A, B) means some backbone read transitions between A and B.
    std::set<std::pair<uint32_t, uint32_t>> endpointWindowPairs;

    // Anchors used by endpoint edges (and their RC mirrors).
    // An anchor in this set is at an endpoint position of its window.
    std::set<uint64_t> endpointAnchors;

    // Per-anchor endpoint flag, indexed by anchor ID.
    // true = endNode (endpoint anchor), false = intraNode.
    // Populated after pass 1 and chain-end promotion.
    vector<bool> isEndpointAnchor;
};
