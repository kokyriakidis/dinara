#ifndef DINARA_LOCAL_UNITIG_GRAPH_HPP
#define DINARA_LOCAL_UNITIG_GRAPH_HPP

// A local subgraph of the global UnitigGraph, starting at one or more oriented unitigs
// and extending out to a specified distance (number of arcs).

// Dinara.
#include <computeLayout.hpp>
#include "OrientedUnitigId.hpp"
#include "iostream.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "array.hpp"
#include <map>
#include "string.hpp"

namespace dinara {
    class LocalUnitigGraphVertex;
    class LocalUnitigGraphEdge;
    class LocalUnitigGraph;

    using LocalUnitigGraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::listS,
        boost::directedS,
        LocalUnitigGraphVertex,
        LocalUnitigGraphEdge
        >;

    class Assembler;
}


class dinara::LocalUnitigGraphVertex {
public:
    OrientedUnitigId orientedUnitigId;
    uint32_t orientedValue;
    uint32_t readCount;
    bool isCircular;
    uint32_t distance;

    LocalUnitigGraphVertex(
        OrientedUnitigId orientedUnitigId,
        uint32_t readCount,
        bool isCircular,
        uint32_t distance) :
        orientedUnitigId(orientedUnitigId),
        orientedValue(orientedUnitigId.getValue()),
        readCount(readCount),
        isCircular(isCircular),
        distance(distance)
    {}

    array<double, 2> position;
};


class dinara::LocalUnitigGraphEdge {
public:
    uint32_t overlapLen = 0;
    uint32_t len = 0;
    uint64_t globalArcId = 0;
    string color;

    LocalUnitigGraphEdge(
        uint32_t overlapLen,
        uint32_t len,
        uint64_t globalArcId) :
        overlapLen(overlapLen),
        len(len),
        globalArcId(globalArcId)
    {}
};


class dinara::LocalUnitigGraph :
    public LocalUnitigGraphBaseClass {
public:
    void addVertex(
        OrientedUnitigId,
        uint32_t readCount,
        bool isCircular,
        uint32_t distance);

    void addEdge(
        OrientedUnitigId from,
        OrientedUnitigId to,
        uint32_t overlapLen,
        uint32_t len,
        uint64_t globalArcId);

    bool vertexExists(OrientedUnitigId) const;
    uint32_t getDistance(OrientedUnitigId) const;

    ComputeLayoutReturnCode computeLayout(
        const string& layoutMethod,
        double timeout);

    void writeSvg(
        const string& svgId,
        uint64_t width,
        uint64_t height,
        double vertexScalingFactor,
        double edgeThicknessScalingFactor,
        uint64_t maxDistance,
        const Assembler&,
        ostream& svg) const;

private:
    std::map<OrientedUnitigId, vertex_descriptor> vertexMap;
};

#endif
