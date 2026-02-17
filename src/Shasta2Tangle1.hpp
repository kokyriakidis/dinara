#pragma once

#include "Shasta2AssemblyGraph.hpp"

#include "memory.hpp"
#include "vector.hpp"

namespace dinara {
    class Shasta2Tangle1;
    class Shasta2TangleMatrix1;
}

class dinara::Shasta2Tangle1 {
public:
    using vertex_descriptor = Shasta2AssemblyGraph::vertex_descriptor;
    using edge_descriptor = Shasta2AssemblyGraph::edge_descriptor;

    Shasta2Tangle1(Shasta2AssemblyGraph&, const vector<vertex_descriptor>&);
    Shasta2Tangle1(Shasta2AssemblyGraph&, vertex_descriptor);
    Shasta2Tangle1(Shasta2AssemblyGraph&, edge_descriptor);

    Shasta2AssemblyGraph& assemblyGraph;

    vector<vertex_descriptor> tangleVertices;
    bool isTangleVertex(vertex_descriptor) const;

    vector<edge_descriptor> entrances;
    void findEntrances();

    vector<edge_descriptor> exits;
    void findExits();

    vector<edge_descriptor> tangleEdges;
    void findTangleEdges();

    shared_ptr<Shasta2TangleMatrix1> tangleMatrixPointer;
    const Shasta2TangleMatrix1& tangleMatrix() const
    {
        return *tangleMatrixPointer;
    }

    class ConnectPair {
    public:
        uint64_t entranceIndex;
        uint64_t exitIndex;

        ConnectPair(
            uint64_t entranceIndex,
            uint64_t exitIndex) :
            entranceIndex(entranceIndex),
            exitIndex(exitIndex)
        {}

        Shasta2AssemblyGraphEdge newEdge;
    };

    vector<ConnectPair> connectPairs;
    bool addConnectPair(uint64_t entranceIndex, uint64_t exitIndex);
    void detangle();

    vector<vertex_descriptor> removedVertices;

    void rerouteEntrances(vector<vertex_descriptor>& newEntranceVertices) const;
    void rerouteExits(vector<vertex_descriptor>& newExitVertices) const;
    void reconnect(
        ConnectPair& connectPair,
        vertex_descriptor v0,
        vertex_descriptor v1) const;
};
