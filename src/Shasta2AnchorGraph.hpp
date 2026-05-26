#pragma once

// Shasta2AnchorGraph.hpp

#include "AnchorWindows.hpp"
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
    Shasta2AnchorGraph(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        const vector<AnchorWindow>& anchorWindows,
        uint64_t minInterWindowCoverage,
        uint64_t threadCount);

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
    void writeGfa(const string& fileName) const;
    void writeBubbleFinderGraph(const string& fileName, bool useForAssemblyOnly = true) const;

    // Save binary data (Shasta2 compatibility).
    void saveAnchorGraph(const string& name = "Shasta2AnchorGraph") const { save(name); }


};
