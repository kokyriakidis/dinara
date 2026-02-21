#pragma once

// Given a directed graph and a start vertex on it:
// - findReachableVertices find all vertices that are reachable from the start vertex
//   moving forward (if direction is 0) or backward (if direction is 1).
// - isReachable returns true if the end vertex is reachable from the start vertex.

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <queue>
#include <set>

namespace dinara {
    template<class Graph> void findReachableVertices(
        const Graph&,
        typename Graph::vertex_descriptor,
        uint64_t direction,
        std::set<typename Graph::vertex_descriptor>&
        );

    template<class Graph> bool isReachable(
        const Graph&,
        typename Graph::vertex_descriptor vStart,
        typename Graph::vertex_descriptor vEnd,
        uint64_t direction);
}



template<class Graph> void dinara::findReachableVertices(
    const Graph& graph,
    typename Graph::vertex_descriptor vStart,
    uint64_t direction,
    std::set<typename Graph::vertex_descriptor>& reachableVertices)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;

    reachableVertices.clear();
    reachableVertices.insert(vStart);

    // Initialize the BFS queue.
    std::queue<vertex_descriptor> q;
    q.push(vStart);

    // Main BFS loop.
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        if(direction == 0) {
            BGL_FORALL_OUTEDGES_T(v0, e, graph, Graph) {
                const vertex_descriptor v1 = target(e, graph);
                if(not reachableVertices.contains(v1)) {
                    reachableVertices.insert(v1);
                    q.push(v1);
                }
            }
        } else {
            BGL_FORALL_INEDGES_T(v0, e, graph, Graph) {
                const vertex_descriptor v1 = source(e, graph);
                if(not reachableVertices.contains(v1)) {
                    reachableVertices.insert(v1);
                    q.push(v1);
                }
            }
        }
    }
}



template<class Graph> bool dinara::isReachable(
    const Graph& graph,
    typename Graph::vertex_descriptor vStart,
    typename Graph::vertex_descriptor vEnd,
    uint64_t direction)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;

    std::set<vertex_descriptor> reachableVertices;
    findReachableVertices(graph, vStart, direction, reachableVertices);

    return reachableVertices.contains(vEnd);

}

