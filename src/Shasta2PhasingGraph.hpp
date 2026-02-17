#pragma once

#include "Shasta2GTest.hpp"
#include "invalid.hpp"

#include <boost/graph/adjacency_list.hpp>

#include "string.hpp"
#include "vector.hpp"

namespace dinara {
    class Shasta2PhasingGraph;
    class Shasta2PhasingGraphVertex;
    class Shasta2PhasingGraphEdge;

    using Shasta2PhasingGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Shasta2PhasingGraphVertex,
        Shasta2PhasingGraphEdge>;
}

class dinara::Shasta2PhasingGraphVertex {
public:
    uint64_t position;
    uint64_t componentId = invalid<uint64_t>;
    uint64_t pathLength = invalid<uint64_t>;

    Shasta2PhasingGraphVertex(uint64_t position) :
        position(position)
    {}
};

class dinara::Shasta2PhasingGraphEdge {
public:
    bool isShortestPathEdge = false;
    Shasta2GTest::Hypothesis bestHypothesis;

    Shasta2PhasingGraphEdge(const Shasta2GTest::Hypothesis& bestHypothesis) :
        bestHypothesis(bestHypothesis)
    {}
};

class dinara::Shasta2PhasingGraph : public Shasta2PhasingGraphBaseClass {
public:
    void addVertex(uint64_t position);

    void addEdge(
        uint64_t position0,
        uint64_t position1,
        const Shasta2GTest::Hypothesis& bestHypothesis);

    uint64_t removeLowDegreeVertices(uint64_t minDegree);

    vector< vector<uint64_t> > components;
    void computeConnectedComponents();

    vector< vector<edge_descriptor> > longestPaths;
    void findLongestPaths();

    void writeGraphviz(const string& fileName) const;

private:
    vector<vertex_descriptor> vertexTable;
};
