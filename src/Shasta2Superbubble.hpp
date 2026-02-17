#pragma once

#include "Shasta2AssemblyGraph.hpp"
#include "DINARA_ASSERT.hpp"

namespace dinara {
    class Shasta2Superbubble;
}

class dinara::Shasta2Superbubble {
public:
    using vertex_descriptor = Shasta2AssemblyGraph::vertex_descriptor;
    using edge_descriptor = Shasta2AssemblyGraph::edge_descriptor;

    const Shasta2AssemblyGraph& assemblyGraph;

    vertex_descriptor sourceVertex;
    vertex_descriptor targetVertex;

    // Sorted by vertex id.
    vector<vertex_descriptor> internalVertices;

    // Sorted by edge id.
    vector<edge_descriptor> sourceEdges;
    vector<edge_descriptor> targetEdges;
    vector<edge_descriptor> internalEdges;

    Shasta2Superbubble(
        const Shasta2AssemblyGraph&,
        vertex_descriptor sourceVertex,
        vertex_descriptor targetVertex);

    void gatherInternalVertices();
    void gatherEdges();

    bool contains(vertex_descriptor v) const;
    bool isBubble() const
    {
        return internalVertices.empty();
    }
    uint64_t sourcePloidy() const
    {
        return sourceEdges.size();
    }
    uint64_t targetPloidy() const
    {
        return targetEdges.size();
    }
    uint64_t ploidy() const
    {
        DINARA_ASSERT(isBubble());
        DINARA_ASSERT(sourcePloidy() == targetPloidy());
        return sourcePloidy();
    }
    bool isTrivial() const
    {
        return internalVertices.empty() && (internalEdges.size() == 1);
    }
};
