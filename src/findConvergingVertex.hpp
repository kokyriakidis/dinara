#pragma once

// Given a vertex vA of a directed graph, this finds the first vertex vB
// (first in topological ordering) such that:
// - Distance(vA, vB) <= maxDistance.
// - All paths that start at vA go through vB.
// If such a vertex is not found, returns Graph::null_vertex().

// The "General" version handles cycles by contracting
// strongly connected components (SCC condensation) before
// computing flow.

// Ported from shasta2.

// Dinara.
#include "CondensedGraph.hpp"
#include "invalid.hpp"
#include "LocalSubgraph.hpp"
#include "DINARA_ASSERT.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/rational.hpp>

// Standard library.
#include "cstdint.hpp"
#include <map>
#include <queue>
#include "utility.hpp"
#include "vector.hpp"



namespace dinara {

   // General version: handles cycles via SCC condensation + flow computation.
   template<class Graph> typename Graph::vertex_descriptor
       findConvergingVertexGeneral(const Graph&, typename Graph::vertex_descriptor, uint64_t maxDistance);
}



// Three graphs are involved here:
// - The input Graph.
// - The LocalGraph created by the BFS.
// - The CondensedGraph created by condensing strongly connected components of the LocalGraph.
template<class Graph> typename Graph::vertex_descriptor dinara::findConvergingVertexGeneral(
    const Graph& graph,
    typename Graph::vertex_descriptor vStart,
    uint64_t maxDistance)
{

    // Do a BFS to create the LocalGraph.
    using V = typename Graph::vertex_descriptor;
    class LocalGraphVertex {
    public:
        V v = Graph::null_vertex();
        uint64_t distance = invalid<uint64_t>;
    };
    using LocalGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, LocalGraphVertex>;
    using VL = typename LocalGraph::vertex_descriptor;
    const LocalGraph localGraph = createLocalSubgraph<Graph, LocalGraph>(graph, vStart, maxDistance);

    // Sanity check.
    const VL vLStart = 0;
    DINARA_ASSERT(localGraph[vLStart].v == vStart);



    // Create the CondensedGraph, in which each strongly connected component is collapsed
    // into a single vertex.
    using Rational = boost::rational<uint64_t>;
    class CondensedGraphVertex {
    public:
        vector<VL> vertices;
        Rational flow = 0;
    };
    using CondensedGraphType =
        boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, CondensedGraphVertex>;
    using VC = typename CondensedGraphType::vertex_descriptor;
    std::map<VL, VC> vertexMap;
    CondensedGraphType condensedGraph = createCondensedGraph<LocalGraph, CondensedGraphType>(localGraph, vertexMap);



    // If our starting vertex is in a non-trivial connected component, just return null_vertex().
    const VC vCStart = vertexMap[vLStart];
    const CondensedGraphVertex& condensedVertexStart = condensedGraph[vCStart];
    if(condensedVertexStart.vertices.size() > 1) {
        return Graph::null_vertex();
    }



    // Do a topological sort of the CondensedGraph.
    vector<VC> topologicalOrder;
    try {
        boost::topological_sort(condensedGraph, back_inserter(topologicalOrder));
    } catch(const boost::not_a_dag&) {
        // This cannot happen because the condensed local subgraph is guaranteed
        // to be acyclic.
        DINARA_ASSERT(0);
    }
    reverse(topologicalOrder.begin(), topologicalOrder.end());

    // Sanity check on the topological ordering.
    // We already checked that the start vertex is not in a non-trivial
    // strong component. So it must correspond to the first CondensedGraph vertex
    // in topological order.
    DINARA_ASSERT(topologicalOrder.size() == num_vertices(condensedGraph));
    DINARA_ASSERT(topologicalOrder.front() == vCStart);



    // Compute flow in topological order.
    // Source gets flow = 1. At each vertex, flow is distributed equally
    // among out-edges. A vertex with flow == 1 means all paths converge there.
    condensedGraph[topologicalOrder.front()].flow = 1;
    for(uint64_t i=1; i<topologicalOrder.size(); i++) {
        const VC vC1 = topologicalOrder[i];
        CondensedGraphVertex&  condensedGraphVertex1 = condensedGraph[vC1];
        condensedGraphVertex1.flow = 0;

        BGL_FORALL_INEDGES_T(vC1, e, condensedGraph, CondensedGraphType){
            const VC vC0 = source(e, condensedGraph);
            // Defensive guard: avoid division by zero in release builds.
            // In a consistent DAG this should never be zero for a parent that has an outgoing edge to vC1.
            const auto outDeg = out_degree(vC0, condensedGraph);
            if(outDeg == 0) {
                continue;
            }
            condensedGraphVertex1.flow += condensedGraph[vC0].flow / outDeg;
        }
    }



    // Find the first vertex, in topological order, that has flow equal to 1
    // and corresponds to a single vertex of the LocalSubgraph and the original Graph.
    for(uint64_t i=1; i<topologicalOrder.size(); i++) {
        const VC vC = topologicalOrder[i];
        const auto& condensedVertex = condensedGraph[vC];
        if(condensedVertex.flow == 1 and condensedVertex.vertices.size() == 1) {
            const VL vL = condensedVertex.vertices.front();
            return localGraph[vL].v;
        }
    }


    // Not found.
    return Graph::null_vertex();
}
