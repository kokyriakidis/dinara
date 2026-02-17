#pragma once

// A Shasta2RestrictedAnchorGraph is a small AnchorGraph constructed
// using only selected portions of the Shasta2Journeys of a given set of OrientedReadIds.

// Shasta
#include "Shasta2AnchorPair.hpp"
#include "shasta2/CycleAvoider.hpp"
#include "orderPairs.hpp"
#include "ReadId.hpp"
#include "shasta2/SimpleMap.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>

// Standard library.
#include "utility.hpp"
#include "vector.hpp"



namespace dinara {

    class Shasta2RestrictedAnchorGraph;
    class Shasta2RestrictedAnchorGraphVertex;
    class Shasta2RestrictedAnchorGraphEdge;

    using Shasta2RestrictedAnchorGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        Shasta2RestrictedAnchorGraphVertex,
        Shasta2RestrictedAnchorGraphEdge>;
    class Shasta2RestrictedAnchorGraph;

    class Shasta2JourneyPortion;

    class Shasta2TangleMatrix1;
}



class dinara::Shasta2JourneyPortion {
public:
    OrientedReadId orientedReadId;
    uint32_t begin;
    uint32_t end;
    Shasta2JourneyPortion(OrientedReadId orientedReadId, uint32_t begin, uint32_t end) :
        orientedReadId(orientedReadId), begin(begin), end(end) {}
};



class dinara::Shasta2RestrictedAnchorGraphVertex : public shasta2::CycleAvoiderVertex {
public:
    Shasta2AnchorId anchorId;
    vector<OrientedReadId> orientedReadIds;
    Shasta2RestrictedAnchorGraphVertex(Shasta2AnchorId anchorId = invalid<Shasta2AnchorId>) : anchorId(anchorId) {}
};



class dinara::Shasta2RestrictedAnchorGraphEdge {
public:
    Shasta2AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;

    bool isOptimalPathEdge = false;

    // Field used by approximateTopologicalSort.
    bool isDagEdge = false;

    // Field used by removeLowCoverageEdges.
    bool wasRemoved = false;
};



class dinara::Shasta2RestrictedAnchorGraph : public Shasta2RestrictedAnchorGraphBaseClass {
public:

    // Constructor using a Shasta2TangleMatrix1.
    Shasta2RestrictedAnchorGraph(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        const Shasta2TangleMatrix1&,
        uint64_t iEntrance,
        uint64_t iExit,
        ostream& html);
    void constructFromTangleMatrix1(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        const Shasta2TangleMatrix1&,
        uint64_t iEntrance,
        uint64_t iExit,
        ostream& html);

    // The journey portions that define this Shasta2RestrictedAnchorGraph.
    vector<Shasta2JourneyPortion> journeyPortions;

    // Fill the journey portions using a Shasta2TangleMatrix1.
    void fillJourneyPortions(
        const Shasta2Journeys&,
        const Shasta2TangleMatrix1&,
        uint64_t iEntrance,
        uint64_t iExit,
        ostream& html);

    // A SimpleMap containing information about all the N distinct AnchorIds
    // used in this Shasta2RestrictedAnchorGraph.
    // The anchorIndex of an Shasta2AnchorId is a number in [0, N) that uniquely
    // identifies each of these distinct AnchorIds.
    // In other words, it is a perfect hash function for these N AnchorIds.
    class AnchorInformation {
    public:
        uint64_t anchorIndex = invalid<uint64_t>;
        vertex_descriptor v = null_vertex();
    };
    shasta2::SimpleMap<Shasta2AnchorId, AnchorInformation> anchorMap;
    vector<Shasta2AnchorId> anchorIds;    // Indexed by anchorIndex.

    void writeAnchorMap() const;
    void writeAnchorIds() const;



    // Create a new vertex and add it to the anchorMap.
    vertex_descriptor addVertex(Shasta2AnchorId);

    // Return the vertex_descriptor corresponding to an Shasta2AnchorId.
    // This asserts if there is not such vertex.
    vertex_descriptor getExistingVertex(Shasta2AnchorId anchorId)
    {
        const auto p = anchorMap.getExisting(anchorId);
        DINARA_ASSERT(p);
        const vertex_descriptor v = p->second.v;
        DINARA_ASSERT(v != null_vertex());
        return v;
    }



    // Return the vertex_descriptor corresponding to an Shasta2AnchorId.
    // This returns null_vertex() there is not such vertex.
    vertex_descriptor getVertex(Shasta2AnchorId anchorId)
    {
        const auto p = anchorMap.getExisting(anchorId);
        DINARA_ASSERT(p);
        const AnchorInformation& anchorInformation = p->second;
        return anchorInformation.v;
    }



    // Find out if a vertex with the given Shasta2AnchorId exists.
    bool vertexExists(Shasta2AnchorId anchorId)
    {
        const auto p = anchorMap.getExisting(anchorId);
        DINARA_ASSERT(p);
        const AnchorInformation& anchorInformation = p->second;
        return anchorInformation.v != null_vertex();
    }



    // A pointer to CycleAvoider, if we are using one.
    // I was not able to get this to compile with a shared_ptr.
    shasta2::CycleAvoider<Shasta2RestrictedAnchorGraph>* cycleAvoider = 0;

    // Only keep vertices that are forward reachable from the
    // vertex at anchorId0 and backward reachable from the vertex at anchorId1.
    void keepBetween(Shasta2AnchorId anchorId0, Shasta2AnchorId anchorId1);

    // Approximate topological sort.
    void approximateTopologicalSort();

    // Find the optimal assembly path.
    // This also sets the isOptimalPathEdge on the edges of the optimal path path.
    void findOptimalPath(
        Shasta2AnchorId anchorId0,
        Shasta2AnchorId anchorId1,
        vector<edge_descriptor>&);

    // Write a table showing which OrientedReadIds are in each vertex.
    // Vertices are written out in rank order.
    void writeOrientedReadsInVertices(ostream& html) const;

    void writeGraphviz(const string& fileName, const vector<Shasta2AnchorId>& highlightVertices) const;
    void writeGraphviz(ostream&, const vector<Shasta2AnchorId>& highlightVertices) const;
    void writeHtml(ostream&, const vector<Shasta2AnchorId>& highlightVertices) const;



    // Gather all the distinct AnchorIds that appear in the JourneyPortions
    // and store them in the anchorMap.
    void gatherAllAnchorIds(const Shasta2Journeys&);

    uint64_t getAnchorIndex(Shasta2AnchorId anchorId)
    {
        const auto p = anchorMap.getExisting(anchorId);
        DINARA_ASSERT(p);
        return p->second.anchorIndex;
    }

    // The anchorIndexes for each Anchor of the JourneyPortions.
    vector< vector<uint64_t> > journeyPortionsAnchorIndexes;
    void fillJourneyPortionsAnchorIndexes(const Shasta2Journeys&);

    // Gather all transitions(anchorIndex0, anchorIndex1) for consecutive
    // anchors in the journey portions. The number of times each
    // transition appears in the journeys is its coverage.
    // Store the transitions in a vector indexed by coverage.
    class Transition {
    public:
        uint64_t anchorIndex0;
        uint64_t anchorIndex1;
    };
    vector< vector<Transition> > transitions;
    void gatherTransitions(ostream& html);

    // Failure modes.
    class NoTransitions : public std::exception {
        const char* what() const noexcept
        {
            return "Shasta2RestrictedAnchorGraph: NoTransitions";
        }
    };
    class Unreachable : public std::exception {
        const char* what() const noexcept
        {
            return "Shasta2RestrictedAnchorGraph: Unreachable";
        }
    };

};
