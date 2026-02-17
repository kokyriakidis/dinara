#include "Shasta2Tangle1.hpp"

#include "Shasta2RestrictedAnchorGraph.hpp"
#include "Shasta2TangleMatrix1.hpp"
#include "DINARA_ASSERT.hpp"

#include <boost/graph/iteration_macros.hpp>

#include <algorithm>
#include <exception>

using namespace dinara;
using namespace std;

Shasta2Tangle1::Shasta2Tangle1(
    Shasta2AssemblyGraph& assemblyGraph,
    const vector<vertex_descriptor>& tangleVerticesArgument) :
    assemblyGraph(assemblyGraph),
    tangleVertices(tangleVerticesArgument)
{
    sort(
        tangleVertices.begin(),
        tangleVertices.end(),
        [&assemblyGraph](const vertex_descriptor v0, const vertex_descriptor v1) {
            return assemblyGraph[v0].id < assemblyGraph[v1].id;
        });

    findEntrances();
    findExits();
    findTangleEdges();

    ostream html(0);
    tangleMatrixPointer = make_shared<Shasta2TangleMatrix1>(assemblyGraph, entrances, exits, html);
}

Shasta2Tangle1::Shasta2Tangle1(
    Shasta2AssemblyGraph& assemblyGraph,
    const vertex_descriptor v) :
    Shasta2Tangle1(assemblyGraph, vector<vertex_descriptor>(1, v))
{}

Shasta2Tangle1::Shasta2Tangle1(
    Shasta2AssemblyGraph& assemblyGraph,
    const edge_descriptor e) :
    Shasta2Tangle1(
        assemblyGraph,
        vector<vertex_descriptor>({source(e, assemblyGraph), target(e, assemblyGraph)}))
{}

void Shasta2Tangle1::findEntrances()
{
    entrances.clear();
    for(const vertex_descriptor v1: tangleVertices) {
        BGL_FORALL_INEDGES(v1, e, assemblyGraph, Shasta2AssemblyGraph) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            if(!isTangleVertex(v0)) {
                entrances.push_back(e);
            }
        }
    }

    sort(
        entrances.begin(),
        entrances.end(),
        [this](const edge_descriptor e0, const edge_descriptor e1) {
            return assemblyGraph[e0].id < assemblyGraph[e1].id;
        });
}

void Shasta2Tangle1::findExits()
{
    exits.clear();
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, Shasta2AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(!isTangleVertex(v1)) {
                exits.push_back(e);
            }
        }
    }

    sort(
        exits.begin(),
        exits.end(),
        [this](const edge_descriptor e0, const edge_descriptor e1) {
            return assemblyGraph[e0].id < assemblyGraph[e1].id;
        });
}

void Shasta2Tangle1::findTangleEdges()
{
    tangleEdges.clear();
    for(const vertex_descriptor v0: tangleVertices) {
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, Shasta2AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            if(isTangleVertex(v1)) {
                tangleEdges.push_back(e);
            }
        }
    }
}

bool Shasta2Tangle1::isTangleVertex(const vertex_descriptor v) const
{
    return binary_search(
        tangleVertices.begin(),
        tangleVertices.end(),
        v,
        [this](const vertex_descriptor x, const vertex_descriptor y) {
            return assemblyGraph[x].id < assemblyGraph[y].id;
        });
}

bool Shasta2Tangle1::addConnectPair(const uint64_t entranceIndex, const uint64_t exitIndex)
{
    DINARA_ASSERT(entranceIndex < entrances.size());
    DINARA_ASSERT(exitIndex < exits.size());

    connectPairs.emplace_back(entranceIndex, exitIndex);
    Shasta2AssemblyGraphEdge& newEdge = connectPairs.back().newEdge;

    const Shasta2AssemblyGraph::vertex_descriptor v0 = target(entrances[entranceIndex], assemblyGraph);
    const Shasta2AssemblyGraph::vertex_descriptor v1 = source(exits[exitIndex], assemblyGraph);
    const Shasta2AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const Shasta2AnchorId anchorId1 = assemblyGraph[v1].anchorId;

    if(anchorId0 == anchorId1) {
        return true;
    }

    try {
        ostream html(0);
        Shasta2RestrictedAnchorGraph restrictedAnchorGraph(
            *assemblyGraph.getAnchorsPointer(),
            *assemblyGraph.getJourneysPointer(),
            tangleMatrix(),
            entranceIndex,
            exitIndex,
            html);

        vector<Shasta2RestrictedAnchorGraph::edge_descriptor> longestPath;
        restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

        for(const Shasta2RestrictedAnchorGraph::edge_descriptor re: longestPath) {
            const auto& rEdge = restrictedAnchorGraph[re];
            if(rEdge.anchorPair.size() == 0) {
                newEdge.clear();
                return false;
            }
            newEdge.push_back(Shasta2AssemblyGraphEdgeStep(rEdge.anchorPair, rEdge.offset));
        }

    } catch(const std::exception&) {
        return false;
    }

    return true;
}

void Shasta2Tangle1::detangle()
{
    vector<vertex_descriptor> newEntranceVertices;
    rerouteEntrances(newEntranceVertices);

    vector<vertex_descriptor> newExitVertices;
    rerouteExits(newExitVertices);

    for(ConnectPair& connectPair: connectPairs) {
        const uint64_t entranceIndex = connectPair.entranceIndex;
        const uint64_t exitIndex = connectPair.exitIndex;
        reconnect(connectPair, newEntranceVertices[entranceIndex], newExitVertices[exitIndex]);
    }

    removedVertices = tangleVertices;
    for(const vertex_descriptor v: tangleVertices) {
        clear_vertex(v, assemblyGraph);
        remove_vertex(v, assemblyGraph);
    }
}

void Shasta2Tangle1::rerouteEntrances(vector<vertex_descriptor>& newEntranceVertices) const
{
    newEntranceVertices.clear();

    for(const edge_descriptor& eOld: entrances) {
        const vertex_descriptor v0Old = source(eOld, assemblyGraph);
        Shasta2AssemblyGraphEdge& edgeOld = assemblyGraph[eOld];
        const Shasta2AnchorId lastAnchorId = edgeOld.back().anchorPair.anchorIdB;

        const vertex_descriptor v1 = add_vertex(
            Shasta2AssemblyGraphVertex(lastAnchorId, assemblyGraph.nextVertexId++),
            assemblyGraph);
        newEntranceVertices.push_back(v1);

        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0Old, v1, Shasta2AssemblyGraphEdge(edgeOld.id), assemblyGraph);
        Shasta2AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        copy(edgeOld.begin(), edgeOld.end(), back_inserter(edgeNew));
        edgeNew.wasAssembled = edgeOld.wasAssembled;
    }
}

void Shasta2Tangle1::rerouteExits(vector<vertex_descriptor>& newExitVertices) const
{
    newExitVertices.clear();

    for(const edge_descriptor& eOld: exits) {
        const vertex_descriptor v1Old = target(eOld, assemblyGraph);
        Shasta2AssemblyGraphEdge& edgeOld = assemblyGraph[eOld];
        const Shasta2AnchorId firstAnchorId = edgeOld.front().anchorPair.anchorIdA;

        const vertex_descriptor v0 = add_vertex(
            Shasta2AssemblyGraphVertex(firstAnchorId, assemblyGraph.nextVertexId++),
            assemblyGraph);
        newExitVertices.push_back(v0);

        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(v0, v1Old, Shasta2AssemblyGraphEdge(edgeOld.id), assemblyGraph);
        Shasta2AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        copy(edgeOld.begin(), edgeOld.end(), back_inserter(edgeNew));
        edgeNew.wasAssembled = edgeOld.wasAssembled;
    }
}

void Shasta2Tangle1::reconnect(
    ConnectPair& connectPair,
    const vertex_descriptor v0,
    const vertex_descriptor v1) const
{
    Shasta2AssemblyGraph::edge_descriptor e;
    tie(e, ignore) = add_edge(v0, v1, assemblyGraph);
    Shasta2AssemblyGraphEdge& newEdge = assemblyGraph[e];
    newEdge = connectPair.newEdge;
    newEdge.id = assemblyGraph.nextEdgeId++;
}
