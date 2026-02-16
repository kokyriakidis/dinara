#pragma once

// Superbubble detection for the AnchorGraph using flow-based
// convergence detection (ported from shasta2).
//
// A superbubble is a subgraph bounded by a source and target vertex
// where all forward paths from source converge at target, and
// all backward paths from target converge at source.

#include "mode3-AnchorGraph.hpp"

// Boost libraries.
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_traits.hpp>

// Standard library.
#include "vector.hpp"
#include <set>
#include <stack>
#include <unordered_set>
#include <utility>

namespace dinara {
    namespace mode3 {
        class AnchorGraphSuperbubble;
        class AnchorGraphSuperbubbleChain;

        // Edge predicate that keeps only non-flagged edges
        // (edges where isNonTransitiveReductionEdge == false).
        struct KeptEdgePredicate {
            const AnchorGraphBaseClass* graph = nullptr;
            KeptEdgePredicate() = default;
            explicit KeptEdgePredicate(const AnchorGraphBaseClass& g) : graph(&g) {}
            bool operator()(AnchorGraphBaseClass::edge_descriptor e) const {
                if(!graph) return true;
                return !(*graph)[e].isNonTransitiveReductionEdge;
            }
        };

        // Verkko-style Onodera bubble end detection
        // (ported from scripts/estimate_unique_local.py:find_bubble_end and
        // scripts/pop_bubbles_coverage_based.py:find_bubble).
        //
        // The graph wrapper is expected to provide:
        // - using vertex_descriptor
        // - out_neighbors(v)
        // - in_neighbors(v)
        //
        // Returns true and sets bubbleEnd when a bubble is found from start.
        template<class Graph>
        bool find_bubble_end(
            const Graph& g,
            typename Graph::vertex_descriptor start,
            typename Graph::vertex_descriptor& bubbleEnd)
        {
            using Vertex = typename Graph::vertex_descriptor;

            const auto startOut = g.out_neighbors(start);
            if(startOut.size() < 2) {
                return false;
            }

            std::vector<Vertex> stack;
            stack.push_back(start);
            std::unordered_set<Vertex> visited;
            std::unordered_set<Vertex> seen;
            seen.insert(start);

            while(!stack.empty()) {
                const Vertex v = stack.back();
                stack.pop_back();

                auto itSeen = seen.find(v);
                if(itSeen == seen.end()) {
                    return false;
                }
                seen.erase(itSeen);

                if(visited.contains(v)) {
                    return false;
                }
                visited.insert(v);

                const auto out = g.out_neighbors(v);
                if(out.empty()) {
                    return false;
                }
                for(const Vertex u: out) {
                    // Reject self-loops and back-edges to the start.
                    if(u == v || u == start) {
                        return false;
                    }
                    // Verkko rejects back-edges into already visited region.
                    if(visited.contains(u)) {
                        return false;
                    }
                    seen.insert(u);

                    const auto parents = g.in_neighbors(u);
                    if(parents.empty()) {
                        return false;
                    }
                    bool hasNonVisitedParent = false;
                    for(const Vertex parent: parents) {
                        if(!visited.contains(parent)) {
                            hasNonVisitedParent = true;
                            break;
                        }
                    }
                    if(!hasNonVisitedParent) {
                        stack.push_back(u);
                    }
                }

                if(stack.size() == 1 && seen.size() == 1) {
                    const Vertex candidate = stack.back();
                    if(*seen.begin() != candidate) {
                        continue;
                    }
                    // Reject cycle closure candidate->start.
                    const auto candidateOut = g.out_neighbors(candidate);
                    for(const Vertex x: candidateOut) {
                        if(x == start) {
                            return false;
                        }
                    }
                    bubbleEnd = candidate;
                    return true;
                }
            }

            return false;
        }

        // Enumerate all superbubbles using Verkko-style Onodera bubble-end search.
        // The graph wrapper is expected to provide:
        // - using vertex_descriptor
        // - vertices()
        // - out_neighbors(v)
        // - in_neighbors(v)
        template<class Graph>
        std::vector<std::pair<typename Graph::vertex_descriptor, typename Graph::vertex_descriptor> >
        find_superbubbles(const Graph& g)
        {
            using Vertex = typename Graph::vertex_descriptor;
            std::set<std::pair<Vertex, Vertex> > uniquePairs;
            for(const Vertex start: g.vertices()) {
                Vertex bubbleEnd{};
                if(find_bubble_end(g, start, bubbleEnd)) {
                    uniquePairs.insert(std::make_pair(start, bubbleEnd));
                }
            }

            std::vector<std::pair<Vertex, Vertex> > result;
            result.reserve(uniquePairs.size());
            for(const auto& p: uniquePairs) {
                result.push_back(p);
            }
            return result;
        }

        // Adapter helper to run Onodera detection directly on Boost graphs.
        template<class BoostGraph>
        std::vector<std::pair<
            typename boost::graph_traits<BoostGraph>::vertex_descriptor,
            typename boost::graph_traits<BoostGraph>::vertex_descriptor> >
        find_superbubbles_boost(const BoostGraph& boostGraph)
        {
            class Adapter {
            public:
                using vertex_descriptor =
                    typename boost::graph_traits<BoostGraph>::vertex_descriptor;

                explicit Adapter(const BoostGraph& g) : g(g) {}

                std::vector<vertex_descriptor> vertices() const
                {
                    std::vector<vertex_descriptor> v;
                    auto [it, end] = boost::vertices(g);
                    for(; it != end; ++it) {
                        v.push_back(*it);
                    }
                    return v;
                }

                std::vector<vertex_descriptor> out_neighbors(vertex_descriptor v) const
                {
                    std::vector<vertex_descriptor> n;
                    auto [it, end] = boost::out_edges(v, g);
                    for(; it != end; ++it) {
                        n.push_back(boost::target(*it, g));
                    }
                    return n;
                }

                std::vector<vertex_descriptor> in_neighbors(vertex_descriptor v) const
                {
                    std::vector<vertex_descriptor> n;
                    auto [it, end] = boost::in_edges(v, g);
                    for(; it != end; ++it) {
                        n.push_back(boost::source(*it, g));
                    }
                    return n;
                }

            private:
                const BoostGraph& g;
            };

            const Adapter adapter(boostGraph);
            return find_superbubbles(adapter);
        }
    }
}



class dinara::mode3::AnchorGraphSuperbubble {
public:
    using vertex_descriptor = AnchorGraphBaseClass::vertex_descriptor;
    using edge_descriptor = AnchorGraphBaseClass::edge_descriptor;

    vertex_descriptor sourceVertex;
    vertex_descriptor targetVertex;

    // The internal vertices, excluding source and target.
    vector<vertex_descriptor> internalVertices;

    // Source edges = out-edges of source vertex (only kept edges).
    // Target edges = in-edges of target vertex (only kept edges).
    // Internal edges = all kept edges within the superbubble.
    vector<edge_descriptor> sourceEdges;
    vector<edge_descriptor> targetEdges;
    vector<edge_descriptor> internalEdges;

    AnchorGraphSuperbubble(
        const AnchorGraph&,
        vertex_descriptor sourceVertex,
        vertex_descriptor targetVertex);

    bool isBubble() const {
        return internalVertices.empty();
    }

    uint64_t sourcePloidy() const {
        return sourceEdges.size();
    }

    uint64_t targetPloidy() const {
        return targetEdges.size();
    }

    uint64_t ploidy() const {
        DINARA_ASSERT(isBubble());
        DINARA_ASSERT(sourcePloidy() == targetPloidy());
        return sourcePloidy();
    }

    bool isTrivial() const {
        return internalVertices.empty() and internalEdges.size() == 1;
    }

private:
    void gatherInternalVertices(const AnchorGraph&);
    void gatherEdges(const AnchorGraph&);
};



class dinara::mode3::AnchorGraphSuperbubbleChain : public vector<AnchorGraphSuperbubble> {
public:
    using vertex_descriptor = AnchorGraphBaseClass::vertex_descriptor;
};
