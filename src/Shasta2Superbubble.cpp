#include "Shasta2Superbubble.hpp"

#include <boost/graph/iteration_macros.hpp>

#include <algorithm>
#include <queue>
#include <set>

using namespace dinara;
using namespace std;

Shasta2Superbubble::Shasta2Superbubble(
    const Shasta2AssemblyGraph& assemblyGraph,
    const vertex_descriptor sourceVertex,
    const vertex_descriptor targetVertex) :
    assemblyGraph(assemblyGraph),
    sourceVertex(sourceVertex),
    targetVertex(targetVertex)
{
    gatherInternalVertices();
    gatherEdges();
}

void Shasta2Superbubble::gatherInternalVertices()
{
    queue<vertex_descriptor> q;
    set<vertex_descriptor> internalVerticesSet;
    q.push(sourceVertex);
    while(!q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, Shasta2AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(v1 == targetVertex) {
                continue;
            }
            if(internalVerticesSet.insert(v1).second) {
                q.push(v1);
            }
        }
    }

    copy(
        internalVerticesSet.begin(),
        internalVerticesSet.end(),
        back_inserter(internalVertices));
    sort(
        internalVertices.begin(),
        internalVertices.end(),
        [this](const vertex_descriptor v0, const vertex_descriptor v1) {
            return assemblyGraph[v0].id < assemblyGraph[v1].id;
        });
}

void Shasta2Superbubble::gatherEdges()
{
    BGL_FORALL_OUTEDGES(sourceVertex, e, assemblyGraph, Shasta2AssemblyGraph) {
        sourceEdges.push_back(e);
    }
    sort(
        sourceEdges.begin(),
        sourceEdges.end(),
        [this](const edge_descriptor e0, const edge_descriptor e1) {
            return assemblyGraph[e0].id < assemblyGraph[e1].id;
        });

    BGL_FORALL_INEDGES(targetVertex, e, assemblyGraph, Shasta2AssemblyGraph) {
        targetEdges.push_back(e);
    }
    sort(
        targetEdges.begin(),
        targetEdges.end(),
        [this](const edge_descriptor e0, const edge_descriptor e1) {
            return assemblyGraph[e0].id < assemblyGraph[e1].id;
        });

    BGL_FORALL_OUTEDGES(sourceVertex, e, assemblyGraph, Shasta2AssemblyGraph) {
        internalEdges.push_back(e);
    }
    for(const vertex_descriptor v: internalVertices) {
        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
            internalEdges.push_back(e);
        }
    }
    sort(
        internalEdges.begin(),
        internalEdges.end(),
        [this](const edge_descriptor e0, const edge_descriptor e1) {
            return assemblyGraph[e0].id < assemblyGraph[e1].id;
        });
}

bool Shasta2Superbubble::contains(const vertex_descriptor v) const
{
    return binary_search(
        internalVertices.begin(),
        internalVertices.end(),
        v,
        [this](const vertex_descriptor x, const vertex_descriptor y) {
            return assemblyGraph[x].id < assemblyGraph[y].id;
        });
}
