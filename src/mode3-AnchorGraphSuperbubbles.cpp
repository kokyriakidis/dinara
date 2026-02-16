// Dinara.
#include "mode3-AnchorGraphSuperbubbles.hpp"
#include "findConvergingVertex.hpp"

using namespace dinara;
using namespace mode3;

// Boost libraries.
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/reverse_graph.hpp>

// Standard library.
#include <map>
#include <queue>
#include <set>
#include <unordered_set>

namespace {

template<class SetType, class Vertex>
bool containsVertex(const SetType& s, const Vertex v)
{
    return s.find(v) != s.end();
}

} // namespace


// Construct a superbubble given source and target vertices.
AnchorGraphSuperbubble::AnchorGraphSuperbubble(
    const AnchorGraph& anchorGraph,
    vertex_descriptor sourceVertex,
    vertex_descriptor targetVertex) :
    sourceVertex(sourceVertex),
    targetVertex(targetVertex)
{
    gatherInternalVertices(anchorGraph);
    gatherEdges(anchorGraph);
}


// BFS from source, stopping at target, collecting internal vertices.
// Only follows kept edges (isNonTransitiveReductionEdge == false).
void AnchorGraphSuperbubble::gatherInternalVertices(const AnchorGraph& anchorGraph)
{
    internalVertices.clear();
    std::set<vertex_descriptor> visited;

    std::queue<vertex_descriptor> q;
    visited.insert(sourceVertex);
    q.push(sourceVertex);

    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        BGL_FORALL_OUTEDGES(v0, e, anchorGraph, AnchorGraph) {
            if(anchorGraph[e].isNonTransitiveReductionEdge) {
                continue;
            }
            const vertex_descriptor v1 = target(e, anchorGraph);
            if(v1 == sourceVertex || v1 == targetVertex) {
                continue;
            }
            if(visited.insert(v1).second) {
                internalVertices.push_back(v1);
                q.push(v1);
            }
        }
    }
}


void AnchorGraphSuperbubble::gatherEdges(const AnchorGraph& anchorGraph)
{
    sourceEdges.clear();
    targetEdges.clear();
    internalEdges.clear();

    std::unordered_set<vertex_descriptor> internalSet;
    internalSet.reserve(internalVertices.size());
    for(const vertex_descriptor v: internalVertices) {
        internalSet.insert(v);
    }

    const auto inSubgraph = [this, &internalSet](vertex_descriptor v) {
        return
            v == sourceVertex ||
            v == targetVertex ||
            containsVertex(internalSet, v);
    };

    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        if(anchorGraph[e].isNonTransitiveReductionEdge) {
            continue;
        }
        const vertex_descriptor v0 = source(e, anchorGraph);
        const vertex_descriptor v1 = target(e, anchorGraph);
        const bool sourceIn = inSubgraph(v0);
        const bool targetIn = inSubgraph(v1);

        if(v0 == sourceVertex && targetIn) {
            sourceEdges.push_back(e);
        }
        if(v1 == targetVertex && sourceIn) {
            targetEdges.push_back(e);
        }
        if(sourceIn && targetIn) {
            internalEdges.push_back(e);
        }
    }
}


void AnchorGraph::findSuperbubbles(
    vector<AnchorGraphSuperbubble>& superbubbles,
    uint64_t maxDistance) const
{
    superbubbles.clear();
    const AnchorGraph& anchorGraph = *this;

    using FilteredGraph =
        boost::filtered_graph<const AnchorGraph, KeptEdgePredicate>;
    const FilteredGraph keptGraph(anchorGraph, KeptEdgePredicate(anchorGraph));
    const boost::reverse_graph<const FilteredGraph> reverseKeptGraph(keptGraph);

    BGL_FORALL_VERTICES(vSource, anchorGraph, AnchorGraph) {
        if(out_degree(vSource, keptGraph) < 2) {
            continue;
        }

        const vertex_descriptor vTarget =
            findConvergingVertexGeneral(keptGraph, vSource, maxDistance);
        if(vTarget == AnchorGraph::null_vertex() || vTarget == vSource) {
            continue;
        }

        const vertex_descriptor vBackward =
            findConvergingVertexGeneral(reverseKeptGraph, vTarget, maxDistance);
        if(vBackward != vSource) {
            continue;
        }

        superbubbles.emplace_back(anchorGraph, vSource, vTarget);
    }
}


void AnchorGraph::removeContainedSuperbubbles(
    vector<AnchorGraphSuperbubble>& superbubbles) const
{
    if(superbubbles.size() <= 1) {
        return;
    }

    // Remove duplicate source/target pairs first.
    {
        std::set<std::pair<vertex_descriptor, vertex_descriptor> > seen;
        vector<AnchorGraphSuperbubble> unique;
        unique.reserve(superbubbles.size());
        for(const auto& superbubble: superbubbles) {
            const auto key =
                std::make_pair(superbubble.sourceVertex, superbubble.targetVertex);
            if(seen.insert(key).second) {
                unique.push_back(superbubble);
            }
        }
        superbubbles.swap(unique);
    }

    if(superbubbles.size() <= 1) {
        return;
    }

    vector<std::unordered_set<vertex_descriptor> > vertexSets(superbubbles.size());
    for(uint64_t i = 0; i < superbubbles.size(); i++) {
        auto& s = vertexSets[i];
        s.reserve(superbubbles[i].internalVertices.size() + 2);
        s.insert(superbubbles[i].sourceVertex);
        s.insert(superbubbles[i].targetVertex);
        for(const vertex_descriptor v: superbubbles[i].internalVertices) {
            s.insert(v);
        }
    }

    vector<bool> isContained(superbubbles.size(), false);
    for(uint64_t i = 0; i < superbubbles.size(); i++) {
        for(uint64_t j = 0; j < superbubbles.size(); j++) {
            if(i == j) {
                continue;
            }
            if(vertexSets[i].size() > vertexSets[j].size()) {
                continue;
            }
            bool subset = true;
            for(const vertex_descriptor v: vertexSets[i]) {
                if(not containsVertex(vertexSets[j], v)) {
                    subset = false;
                    break;
                }
            }
            if(subset) {
                isContained[i] = true;
                break;
            }
        }
    }

    vector<AnchorGraphSuperbubble> kept;
    kept.reserve(superbubbles.size());
    for(uint64_t i = 0; i < superbubbles.size(); i++) {
        if(!isContained[i]) {
            kept.push_back(superbubbles[i]);
        }
    }
    superbubbles.swap(kept);
}


void AnchorGraph::findSuperbubbleChains(
    const vector<AnchorGraphSuperbubble>& superbubbles,
    vector<AnchorGraphSuperbubbleChain>& chains) const
{
    chains.clear();
    if(superbubbles.empty()) {
        return;
    }

    std::map<vertex_descriptor, vector<uint64_t> > sourceToBubble;
    std::map<vertex_descriptor, vector<uint64_t> > targetToBubble;
    for(uint64_t i = 0; i < superbubbles.size(); i++) {
        sourceToBubble[superbubbles[i].sourceVertex].push_back(i);
        targetToBubble[superbubbles[i].targetVertex].push_back(i);
    }

    vector<bool> used(superbubbles.size(), false);

    // Start chains from bubbles that do not have predecessors.
    for(uint64_t i = 0; i < superbubbles.size(); i++) {
        if(used[i]) {
            continue;
        }
        auto itPred = targetToBubble.find(superbubbles[i].sourceVertex);
        if(itPred != targetToBubble.end() && !itPred->second.empty()) {
            continue;
        }

        AnchorGraphSuperbubbleChain chain;
        uint64_t current = i;
        while(true) {
            if(used[current]) {
                break;
            }
            used[current] = true;
            chain.push_back(superbubbles[current]);

            const vertex_descriptor nextSource = superbubbles[current].targetVertex;
            auto itSucc = sourceToBubble.find(nextSource);
            if(itSucc == sourceToBubble.end() || itSucc->second.size() != 1) {
                break;
            }
            const uint64_t next = itSucc->second.front();
            if(used[next]) {
                break;
            }
            current = next;
        }

        if(!chain.empty()) {
            chains.push_back(chain);
        }
    }

    // Any remaining bubbles become singleton chains.
    for(uint64_t i = 0; i < superbubbles.size(); i++) {
        if(used[i]) {
            continue;
        }
        AnchorGraphSuperbubbleChain chain;
        chain.push_back(superbubbles[i]);
        chains.push_back(chain);
        used[i] = true;
    }
}
