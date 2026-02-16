#pragma once

// Given a directed graph, this does a forward BFS from a given start vertex
// up to a specified distance to construct a local subgraph.

// This is intrusive: the LocalSubgraph vertex is required
// to have the following two fields:
// - Graph::vertex_descriptor v
// - uint64_t distance

// Ported from shasta2.

#include "DINARA_ASSERT.hpp"

#include "cstdint.hpp"
#include <map>
#include <queue>
#include "vector.hpp"


namespace dinara {

    template<class Graph, class LocalSubgraph> LocalSubgraph createLocalSubgraph(
        const Graph&,
        typename Graph::vertex_descriptor startVertex,
        uint64_t maxDistance);
}



template<class Graph, class LocalSubgraph> LocalSubgraph dinara::createLocalSubgraph(
    const Graph& graph,
    typename Graph::vertex_descriptor vStart,
    uint64_t maxDistance)
{
    using VertexDescriptor = typename Graph::vertex_descriptor;
    using LocalVertexDescriptor = typename LocalSubgraph::vertex_descriptor;
    using OutEdgeIterator = typename Graph::out_edge_iterator;

    LocalSubgraph localSubgraph;

    // Create the local start vertex.
    const LocalVertexDescriptor lvStart = add_vertex(localSubgraph);
    auto& localVertex = localSubgraph[lvStart];
    localVertex.v = vStart;
    localVertex.distance = 0;

    // Map Graph vertices to LocalGraph vertices.
    std::map<VertexDescriptor, LocalVertexDescriptor> m;
    m.insert({vStart, lvStart});

    // Initialize the BFS.
    std::queue<VertexDescriptor> q;
    q.push(vStart);



    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a global vertex.
        const VertexDescriptor v0 = q.front();
        q.pop();
        auto it0 = m.find(v0);
        DINARA_ASSERT(it0 != m.end());

        // Get the corresponding local vertex.
        const LocalVertexDescriptor& lv0 = it0->second;
        const auto& localVertex0 = localSubgraph[lv0];
        const uint64_t distance0 = localVertex0.distance;
        const uint64_t distance1 = distance0 + 1;

        // Loop over its out-edges.
        OutEdgeIterator itOut, itOutEnd;
        tie(itOut, itOutEnd) = out_edges(v0, graph);
        for(; itOut!=itOutEnd; ++itOut) {
            const auto e = *itOut;
            const VertexDescriptor v1 = target(e, graph);

            // Check if we already encountered v1.
            auto it1 = m.find(v1);

            if(it1 == m.end()) {

                // We have not yet encountered v1.
                const LocalVertexDescriptor lv1 = add_vertex(localSubgraph);
                auto& localVertex1 = localSubgraph[lv1];
                localVertex1.v = v1;
                localVertex1.distance = distance1;
                m.insert({v1, lv1});
                add_edge(lv0, lv1, localSubgraph);
                if(distance1 < maxDistance) {
                   q.push(v1);
                }
            } else {

                // We already encountered v1. Just add this edge.
                LocalVertexDescriptor lv1 = it1->second;
                add_edge(lv0, lv1, localSubgraph);
            }

        }
    }

    return localSubgraph;
}
