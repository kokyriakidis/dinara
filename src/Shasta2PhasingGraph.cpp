#include "Shasta2PhasingGraph.hpp"

#include "DINARA_ASSERT.hpp"
#include "orderVectors.hpp"

#include <boost/graph/iteration_macros.hpp>

#include <algorithm>
#include <fstream>

using namespace dinara;
using namespace std;

void Shasta2PhasingGraph::addVertex(const uint64_t position)
{
    if(position >= vertexTable.size()) {
        vertexTable.resize(position + 1, null_vertex());
    }

    DINARA_ASSERT(vertexTable[position] == null_vertex());

    const vertex_descriptor v = boost::add_vertex(Shasta2PhasingGraphVertex(position), *this);
    vertexTable[position] = v;
}

void Shasta2PhasingGraph::addEdge(
    const uint64_t position0,
    const uint64_t position1,
    const Shasta2GTest::Hypothesis& bestHypothesis)
{
    DINARA_ASSERT(position0 < vertexTable.size());
    const vertex_descriptor v0 = vertexTable[position0];
    DINARA_ASSERT(v0 != null_vertex());

    DINARA_ASSERT(position1 < vertexTable.size());
    const vertex_descriptor v1 = vertexTable[position1];
    DINARA_ASSERT(v1 != null_vertex());

    boost::add_edge(v0, v1, Shasta2PhasingGraphEdge(bestHypothesis), *this);
}

void Shasta2PhasingGraph::writeGraphviz(const string& fileName) const
{
    const Shasta2PhasingGraph& phasingGraph = *this;

    ofstream dot(fileName);
    dot << "digraph Shasta2PhasingGraph {\n";

    BGL_FORALL_VERTICES(v, phasingGraph, Shasta2PhasingGraph) {
        const Shasta2PhasingGraphVertex& vertex = phasingGraph[v];
        dot << vertex.position << " [label=\"" <<
            vertex.position << "\\n" <<
            vertex.componentId << "\\n" <<
            vertex.pathLength << "\"]";
        dot << ";\n";
    }

    BGL_FORALL_EDGES(e, phasingGraph, Shasta2PhasingGraph) {
        const vertex_descriptor v0 = source(e, phasingGraph);
        const vertex_descriptor v1 = target(e, phasingGraph);

        dot << phasingGraph[v0].position << "->";
        dot << phasingGraph[v1].position;
        if(phasingGraph[e].isShortestPathEdge) {
            dot << " [color=DarkOrange]";
        }
        dot << ";\n";
    }

    dot << "}\n";
}

uint64_t Shasta2PhasingGraph::removeLowDegreeVertices(const uint64_t minDegree)
{
    Shasta2PhasingGraph& phasingGraph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, phasingGraph, Shasta2PhasingGraph) {
        if(in_degree(v, phasingGraph) + out_degree(v, phasingGraph) < minDegree) {
            verticesToBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        vertexTable[phasingGraph[v].position] = null_vertex();
        boost::clear_vertex(v, phasingGraph);
        boost::remove_vertex(v, phasingGraph);
    }

    return verticesToBeRemoved.size();
}

void Shasta2PhasingGraph::computeConnectedComponents()
{
    Shasta2PhasingGraph& phasingGraph = *this;
    const uint64_t n = vertexTable.size();

    vector<uint64_t> parent(n);
    for(uint64_t position=0; position<n; position++) {
        parent[position] = position;
    }
    const auto findRoot = [&parent](uint64_t x) {
        while(parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    const auto unite = [&findRoot, &parent](uint64_t x, uint64_t y) {
        const uint64_t rootX = findRoot(x);
        const uint64_t rootY = findRoot(y);
        if(rootX != rootY) {
            parent[rootY] = rootX;
        }
    };

    BGL_FORALL_EDGES(e, phasingGraph, Shasta2PhasingGraph) {
        const vertex_descriptor v0 = source(e, phasingGraph);
        const vertex_descriptor v1 = target(e, phasingGraph);

        const uint64_t position0 = phasingGraph[v0].position;
        const uint64_t position1 = phasingGraph[v1].position;

        DINARA_ASSERT(position0 < n);
        DINARA_ASSERT(position1 < n);

        unite(position0, position1);
    }

    vector< vector<uint64_t> > componentTable(n);
    for(uint64_t position=0; position<n; position++) {
        const uint64_t componentId = findRoot(position);
        componentTable[componentId].push_back(position);
    }

    components.clear();
    for(const vector<uint64_t>& component: componentTable) {
        if(component.size() > 1) {
            components.push_back(component);
        }
    }

    sort(components.begin(), components.end(), OrderVectorsByDecreasingSize<uint64_t>());

    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<uint64_t>& component = components[componentId];
        for(const uint64_t position: component) {
            const vertex_descriptor v = vertexTable[position];
            DINARA_ASSERT(v != null_vertex());
            phasingGraph[v].componentId = componentId;
        }
    }
}

void Shasta2PhasingGraph::findLongestPaths()
{
    Shasta2PhasingGraph& phasingGraph = *this;

    for(uint64_t position0=0; position0<vertexTable.size(); position0++) {
        const vertex_descriptor v0 = vertexTable[position0];
        if(v0 == null_vertex()) {
            continue;
        }

        Shasta2PhasingGraphVertex& vertex0 = phasingGraph[v0];
        if(in_degree(v0, phasingGraph) == 0) {
            vertex0.pathLength = 0;
        }

        const uint64_t pathLength0 = vertex0.pathLength;
        const uint64_t pathLength1 = pathLength0 + 1;

        BGL_FORALL_OUTEDGES(v0, e, phasingGraph, Shasta2PhasingGraph) {
            const vertex_descriptor v1 = target(e, phasingGraph);
            Shasta2PhasingGraphVertex& vertex1 = phasingGraph[v1];

            if(vertex1.pathLength == invalid<uint64_t>) {
                vertex1.pathLength = pathLength1;
            } else {
                vertex1.pathLength = max(vertex1.pathLength, pathLength1);
            }
        }
    }

    longestPaths.resize(components.size());
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<uint64_t>& component = components[componentId];
        vector<edge_descriptor>& path = longestPaths[componentId];

        uint64_t maxPathLength = 0;
        vertex_descriptor vLast = null_vertex();
        for(const uint64_t position: component) {
            const vertex_descriptor v = vertexTable[position];
            const Shasta2PhasingGraphVertex& vertex = phasingGraph[v];
            DINARA_ASSERT(v != null_vertex());

            if(vLast == null_vertex()) {
                vLast = v;
                maxPathLength = vertex.pathLength;
            } else if(vertex.pathLength > maxPathLength) {
                vLast = v;
                maxPathLength = vertex.pathLength;
            }
        }

        vertex_descriptor v1 = vLast;
        while(in_degree(v1, phasingGraph) > 0) {
            const Shasta2PhasingGraphVertex& vertex1 = phasingGraph[v1];
            const uint64_t pathLength1 = vertex1.pathLength;
            const uint64_t pathLength0 = pathLength1 - 1;
            bool found = false;
            BGL_FORALL_INEDGES(v1, e, phasingGraph, Shasta2PhasingGraph) {
                const vertex_descriptor v0 = source(e, phasingGraph);
                if(phasingGraph[v0].pathLength == pathLength0) {
                    phasingGraph[e].isShortestPathEdge = true;
                    path.push_back(e);
                    v1 = v0;
                    found = true;
                    break;
                }
            }
            DINARA_ASSERT(found);
        }
        reverse(path.begin(), path.end());
    }
}
