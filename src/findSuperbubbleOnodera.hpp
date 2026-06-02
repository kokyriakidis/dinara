#pragma once

// Superbubble detection using the Onodera et al. 2013 algorithm.
// Reference: "Detecting Superbubbles in Assembly Graphs", Onodera et al. 2013.
// Ported from verkko's pop_bubbles_coverage_based.py find_bubble().
//
// Given a directed graph and a start vertex s with out-degree >= 2,
// finds the exit vertex t such that all paths from s converge at t
// and the subgraph between s and t is acyclic.
//
// The algorithm maintains:
//   - visited: vertices that have been fully processed
//   - seen: vertices reached but not yet processable
//   - S (ready): vertices with all incoming edges from visited vertices
//
// A vertex becomes ready when all its in-graph predecessors are visited.
// The superbubble exit is found when exactly one vertex is ready and
// no other vertices are pending.

#include <boost/graph/iteration_macros.hpp>

#include "cstdint.hpp"
#include <unordered_set>
#include "vector.hpp"

namespace dinara {

    // Find the exit vertex of a superbubble starting at vStart.
    // Returns Graph::null_vertex() if no superbubble is found.
    // maxSize limits the number of internal vertices (0 = no limit).
    template<class Graph>
    typename Graph::vertex_descriptor findSuperbubbleOnodera(
        const Graph& graph,
        typename Graph::vertex_descriptor vStart,
        uint64_t maxSize = 0);
}



template<class Graph>
typename Graph::vertex_descriptor dinara::findSuperbubbleOnodera(
    const Graph& graph,
    typename Graph::vertex_descriptor vStart,
    uint64_t maxSize)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;

    // Start vertex must have out-degree >= 2.
    if(out_degree(vStart, graph) < 2) {
        return Graph::null_vertex();
    }

    // S: stack of vertices ready to process (all predecessors visited).
    vector<vertex_descriptor> S;
    S.push_back(vStart);

    // visited: fully processed vertices.
    std::unordered_set<vertex_descriptor> visited;

    // seen: reached but not yet ready vertices.
    std::unordered_set<vertex_descriptor> seen;
    seen.insert(vStart);

    while(!S.empty()) {
        const vertex_descriptor v = S.back();
        S.pop_back();

        // v must be in seen and not in visited.
        seen.erase(v);
        visited.insert(v);

        // Size limit on internal vertices.
        if(v != vStart && maxSize > 0 && visited.size() > maxSize) {
            return Graph::null_vertex();
        }

        // Dead end: not a superbubble.
        if(out_degree(v, graph) == 0) {
            return Graph::null_vertex();
        }

        // Explore out-edges.
        BGL_FORALL_OUTEDGES_T(v, e, graph, Graph) {
            const vertex_descriptor u = target(e, graph);

            // Self-loop: reject.
            if(u == v) {
                return Graph::null_vertex();
            }

            // Edge back to start: cycle, reject.
            if(u == vStart) {
                return Graph::null_vertex();
            }

            // Already visited: back-edge, reject.
            if(visited.count(u)) {
                return Graph::null_vertex();
            }

            seen.insert(u);

            // Check if all predecessors of u are visited.
            bool hasNonvisitedParent = false;
            BGL_FORALL_INEDGES_T(u, ie, graph, Graph) {
                const vertex_descriptor parent = source(ie, graph);
                if(!visited.count(parent)) {
                    hasNonvisitedParent = true;
                    break;
                }
            }

            // If all predecessors visited, u is ready.
            if(!hasNonvisitedParent) {
                S.push_back(u);
            }
        }

        // Superbubble exit condition:
        // exactly one vertex ready, no other vertices pending.
        if(S.size() == 1 && seen.size() == 1 && S.back() == *seen.begin()) {
            const vertex_descriptor t = S.back();

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
