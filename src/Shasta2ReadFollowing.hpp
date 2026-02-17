#pragma once

// Read following in the Shasta2AssemblyGraph.

// Shasta.
#include "Shasta2AssemblyGraph.hpp"
#include "Shasta2SegmentStepSupport.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <array>
#include <map>
#include <random>
#include <set>



namespace dinara {

    namespace Shasta2ReadFollowing {

        class Graph;
        class Vertex;
        class Edge;
        using GraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            Vertex,
            Edge>;

        class PathGraph;
        class PathGraphVertex;
        class PathGraphEdge;
        using PathGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            PathGraphVertex,
            PathGraphEdge>;

        // A Segment is an edge of the Shasta2AssemblyGraph.
        using Segment = Shasta2AssemblyGraph::edge_descriptor;
    }
}



class dinara::Shasta2ReadFollowing::Vertex {
public:
    // A Segment is an edge of the Shasta2AssemblyGraph.
    Segment segment;

    // The sequence length or estimated offset of this Shasta2AssemblyGraph edge.
    uint64_t length = invalid<uint64_t>;

    double coverage = 0.;

    Vertex(const Shasta2AssemblyGraph&, Segment);

    vector<Shasta2SegmentStepSupport> initialSupport;
    vector<Shasta2SegmentStepSupport> finalSupport;

    // These are needed for efficient creation of random paths.
    // We cannot use vecS for the adjacency_list because we
    // need to remove vertices and edges (for pruning).
    // These are filled by fillConnectivity after
    // pruning is done and the graph will no longer change.
    vector<GraphBaseClass::edge_descriptor> outEdges;
    vector<GraphBaseClass::edge_descriptor> inEdges;

};



class dinara::Shasta2ReadFollowing::Edge {
public:
    Edge(const Shasta2AssemblyGraph&, Segment, Segment);
    Shasta2SegmentPairInformation segmentPairInformation;
};



class dinara::Shasta2ReadFollowing::Graph : public GraphBaseClass {
public:
    Graph(const Shasta2AssemblyGraph&);

private:
    const Shasta2AssemblyGraph& assemblyGraph;

    // Initial creation.
    std::map<Segment, vertex_descriptor> vertexMap;
    void createVertices();
    void createEdges();

    // Prune short leaves.
    void prune();
    bool pruneIteration();

    void fillConnectivity();

    // Output.
    void write(const string& name) const;
    void writeCsv(const string& name) const;
    void writeVerticesCsv(const string& name) const;
    void writeEdgesCsv(const string& name) const;
    void writeGraphviz(const string& name) const;

public:



    // Random paths.
    // These functions find a random path starting at the given vertex.
    // Direction is 0 for forward and 1 backward.
    // The path ends when one of the stop vertices is encountered,
    // but can end sooner.
    // Note these are paths in the Shasta2ReadFollowing::Graph
    // but not in the Shasta2AssemblyGraph.
    template<std::uniform_random_bit_generator RandomGenerator> void findRandomPath(
        vertex_descriptor, uint64_t direction,
        RandomGenerator&,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& stopVertices);
    template<std::uniform_random_bit_generator RandomGenerator> void findRandomForwardPath(
        vertex_descriptor,
        RandomGenerator&,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& stopVertices);
    template<std::uniform_random_bit_generator RandomGenerator> void findRandomBackwardPath(
        vertex_descriptor,
        RandomGenerator&,
        vector<vertex_descriptor>& path,
        const std::set<vertex_descriptor>& stopVertices);



    // Find assembly paths.
    void findPaths(vector< vector<Segment> >& assemblyPaths);

    // Python callable.
    void writeRandomPath(Segment, uint64_t direction);
    void writePaths();
};



// Each PathGraph vertex corresponds to a long segment.
class dinara::Shasta2ReadFollowing::PathGraphVertex {
public:
    Segment segment;
};



class dinara::Shasta2ReadFollowing::PathGraphEdge {
public:

    // Store information for each direction.
    class Info {
    public:

        // The number of paths between these two vertices found in each direction.
        uint64_t pathCount = 0;

        // The longest of the paths.
        // These are paths in the Shasta2ReadFollowing::Graph
        // but not in the Shasta2AssemblyGraph.
        vector<Graph::vertex_descriptor> path;
    };

    array<Info, 2> infos;

    uint64_t forwardPathCount() const
    {
        return infos[0].pathCount;
    }
    uint64_t backwardPathCount() const
    {
        return infos[1].pathCount;
    }

    const vector<Graph::vertex_descriptor>& longestPath() const
    {
        return
            (infos[0].path.size() >= infos[1].path.size()) ?
            infos[0].path :
            infos[1].path;

    }
};



// Class used to store paths between long segments.
class dinara::Shasta2ReadFollowing::PathGraph : public PathGraphBaseClass {
public:
    PathGraph(const Shasta2AssemblyGraph&);
    void removeNonBestEdges();

    // Graphviz output.
    void writeGraphviz() const;
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

private:
    const Shasta2AssemblyGraph& assemblyGraph;

    uint64_t forwardPathCount(vertex_descriptor) const;
    uint64_t backwardPathCount(vertex_descriptor) const;

    double forwardPathFraction(edge_descriptor) const;
    double backwardPathFraction(edge_descriptor) const;

    edge_descriptor bestInEdge(vertex_descriptor) const;
    edge_descriptor bestOutEdge(vertex_descriptor) const;

};

namespace dinara {
    class Shasta2Journeys;
    class Shasta2ReadFollowingGraph {
    public:
        explicit Shasta2ReadFollowingGraph(const Shasta2AssemblyGraph& assemblyGraph) :
            assemblyGraph(assemblyGraph)
        {}

        void findPaths(
            const Shasta2Journeys& journeys,
            uint64_t minPathLength,
            uint64_t maxPathCount,
            vector< vector<Shasta2AssemblyGraph::edge_descriptor> >& assemblyPaths) const;

    private:
        const Shasta2AssemblyGraph& assemblyGraph;
    };
}
