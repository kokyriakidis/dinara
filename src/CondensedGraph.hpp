#pragma once

// This computes the condensed graph of a directed graph.
// The condensed graph is obtained by contracting each
// strongly connected component of the original graph
// into a single vertex.
// The condensed graph is guaranteed to be acyclic.
// https://en.wikipedia.org/wiki/Strongly_connected_component

// This is intrusive: the CondensedGraph type is required
// to have a field that will contain the vertices of the
// strongly connected component of the original graph it
// corresponds to:
// vector<Graph::vertex_descriptor> vertices.

// Ported from shasta2.

#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>

#include <map>
#include "vector.hpp"


namespace dinara {

    template<class Graph, class CondensedGraph> CondensedGraph createCondensedGraph(
        const Graph&,
        std::map<typename Graph::vertex_descriptor, typename CondensedGraph::vertex_descriptor>& vertexMap);
}



template<class Graph, class CondensedGraph> CondensedGraph dinara::createCondensedGraph(
    const Graph& graph,
    std::map<typename Graph::vertex_descriptor, typename CondensedGraph::vertex_descriptor>& vertexMap)
{
    using VertexDescriptor = typename Graph::vertex_descriptor;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    std::map<VertexDescriptor, uint64_t> vertexIndexMap;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        vertexIndexMap.insert({v, vertexIndex++});
    }

    // Compute strong components.
    std::map<VertexDescriptor, uint64_t> componentMap;
    boost::strong_components(
        graph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexIndexMap)));

    // Gather the vertices in each strong component.
    std::map<uint64_t, vector<VertexDescriptor> > componentVertices;
    for(const auto& p: componentMap) {
        componentVertices[p.second].push_back(p.first);
    }

    // Each strong component generates a vertex of the condensed graph.
    CondensedGraph condensedGraph;
    vertexMap.clear();
    for(const auto& p: componentVertices) {
        const vector<VertexDescriptor>& vertices = p.second;
        const typename CondensedGraph::vertex_descriptor cv = add_vertex(condensedGraph);
        condensedGraph[cv].vertices = vertices;
        for(const VertexDescriptor v: vertices) {
            vertexMap.insert({v, cv});
        }
    }

    // Now create edges of the condensed graph.
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        const VertexDescriptor v0 = source(e, graph);
        const VertexDescriptor v1 = target(e, graph);
        const typename CondensedGraph::vertex_descriptor cv0 = vertexMap[v0];
        const typename CondensedGraph::vertex_descriptor cv1 = vertexMap[v1];
        if(cv1 != cv0) {
            add_edge(cv0, cv1, condensedGraph);
        }
    }

    return condensedGraph;
}
