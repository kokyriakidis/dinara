#pragma once

// Superbubble detection using the Onodera et al. 2013 algorithm.
// Reference: "Detecting Superbubbles in Assembly Graphs", Onodera et al. 2013.
//
// Given a directed graph and a start vertex s with out-degree >= 2,
// finds the exit vertex t such that all paths from s converge at t
// and the subgraph between s and t is acyclic.
//
// Unlike findConvergingVertexGeneral (which uses SCC condensation + flow),
// this algorithm directly checks the superbubble definition:
// a vertex becomes processable only when all its incoming edges come
// from already-visited vertices. Cycles cause the algorithm to stall
// and return null_vertex().

#include <boost/graph/iteration_macros.hpp>

#include "cstdint.hpp"
#include <unordered_map>
#include <unordered_set>
#include "vector.hpp"

namespace dinara {

    // Find the exit vertex of a superbubble starting at vStart.
    // Returns Graph::null_vertex() if no superbubble is found.
    // maxCount limits the number of vertices explored (0 = no limit).
    template<class Graph>
    typename Graph::vertex_descriptor findSuperbubbleOnodera(
        const Graph& graph,
        typename Graph::vertex_descriptor vStart,
        uint64_t maxCount = 0);
}



template<class Graph>
typename Graph::vertex_descriptor dinara::findSuperbubbleOnodera(
    const Graph& graph,
    typename Graph::vertex_descriptor vStart,
    uint64_t maxCount)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;

    // Start vertex must have out-degree >= 2 (excluding self-loops).
    {
        uint64_t nonSelfOutDegree = 0;
        BGL_FORALL_OUTEDGES_T(vStart, e, graph, Graph) {
            if(target(e, graph) != vStart) {
                ++nonSelfOutDegree;
            }
        }
        if(nonSelfOutDegree < 2) {
            return Graph::null_vertex();
        }
    }

    // Vertices that have been fully processed (all out-edges explored).
    std::unordered_set<vertex_descriptor> visited;

    // For each reached-but-not-yet-processable vertex,
    // track how many incoming edges remain unvisited.
    std::unordered_map<vertex_descriptor, uint64_t> remainingIncoming;

    // Vertices with all incoming edges from visited vertices (ready to process).
    vector<vertex_descriptor> ready;

    // Count of reached vertices that are not yet ready.
    uint64_t notReadyCount = 0;

    // Seed with the start vertex.
    ready.push_back(vStart);
    // Mark vStart as reached with 0 remaining incoming
    // (we don't count edges from outside the bubble).
    remainingIncoming[vStart] = 0;

    while(!ready.empty()) {

        // Check vertex count limit.
        if(maxCount > 0 && visited.size() + notReadyCount > maxCount) {
            return Graph::null_vertex();
        }

        const vertex_descriptor v = ready.back();
        ready.pop_back();
        visited.insert(v);

        // Dead-end: not a superbubble.
        if(out_degree(v, graph) == 0) {
            return Graph::null_vertex();
        }

        // Explore out-edges.
        BGL_FORALL_OUTEDGES_T(v, e, graph, Graph) {
            const vertex_descriptor w = target(e, graph);

            // Self-loop on an internal vertex: reject.
            if(w == v) {
                return Graph::null_vertex();
            }

            // Back-edge to start: cycle, reject.
            if(w == vStart) {
                return Graph::null_vertex();
            }

            if(remainingIncoming.find(w) == remainingIncoming.end()) {
                // First time reaching w.
                // Initialize its remaining incoming count to its full in-degree.
                notReadyCount++;
                remainingIncoming[w] = in_degree(w, graph);
            }

            // Account for this edge.
            auto& rem = remainingIncoming[w];
            rem--;

            // If all incoming edges are now from visited vertices, w is ready.
            if(rem == 0) {
                ready.push_back(w);
                notReadyCount--;
            }
        }

        // Check superbubble exit condition:
        // exactly one vertex ready, nothing pending.
        if(ready.size() == 1 && notReadyCount == 0) {
            const vertex_descriptor t = ready.back();

            // Check that exit doesn't have an edge back to start.
            BGL_FORALL_OUTEDGES_T(t, e, graph, Graph) {
                if(target(e, graph) == vStart) {
                    return Graph::null_vertex();
                }
            }

            return t;
        }
    }

    // Ran out of processable vertices (cycle or disconnected interior).
    return Graph::null_vertex();
}
