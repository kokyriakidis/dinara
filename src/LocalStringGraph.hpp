#ifndef DINARA_LOCAL_STRING_GRAPH_HPP
#define DINARA_LOCAL_STRING_GRAPH_HPP

// A local subgraph of the global StringGraph, starting at one or more oriented reads
// and extending out to a specified distance (number of arcs).

// Dinara.
#include <computeLayout.hpp>
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "array.hpp"
#include <map>
#include "string.hpp"

namespace dinara {

    class LocalStringGraphVertex;
    class LocalStringGraphEdge;
    class LocalStringGraph;

    using LocalStringGraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::listS,
        boost::directedS,
        LocalStringGraphVertex,
        LocalStringGraphEdge
        >;

    class Assembler;
}


class dinara::LocalStringGraphVertex {
public:
    OrientedReadId orientedReadId;
    uint32_t orientedReadIdValue;
    uint32_t markerCount;
    bool isChimeric;
    uint32_t distance;

    LocalStringGraphVertex(
        OrientedReadId orientedReadId,
        uint32_t markerCount,
        bool isChimeric,
        uint32_t distance) :
        orientedReadId(orientedReadId),
        orientedReadIdValue(uint32_t(orientedReadId.getValue())),
        markerCount(markerCount),
        isChimeric(isChimeric),
        distance(distance)
    {}

    array<double, 2> position;
};


class dinara::LocalStringGraphEdge {
public:
    uint32_t overlapLen = 0;
    uint32_t len = 0;
    uint64_t globalArcId = 0;
    string color;

    LocalStringGraphEdge(
        uint32_t overlapLen,
        uint32_t len,
        uint64_t globalArcId) :
        overlapLen(overlapLen),
        len(len),
        globalArcId(globalArcId)
    {}
};


class dinara::LocalStringGraph :
    public LocalStringGraphBaseClass {
public:
    void addVertex(
        OrientedReadId,
        uint32_t markerCount,
        bool isChimeric,
        uint32_t distance);

    void addEdge(
        OrientedReadId from,
        OrientedReadId to,
        uint32_t overlapLen,
        uint32_t len,
        uint64_t globalArcId);

    bool vertexExists(OrientedReadId) const;
    uint32_t getDistance(OrientedReadId) const;

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
    std::map<OrientedReadId, vertex_descriptor> vertexMap;
};

#endif
