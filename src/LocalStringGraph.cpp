// Dinara.
#include "LocalStringGraph.hpp"
#include "Assembler.hpp"
#include "writeGraph.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <cmath>



void LocalStringGraph::addVertex(
    OrientedReadId orientedReadId,
    uint32_t markerCount,
    bool isChimeric,
    uint32_t distance)
{
    DINARA_ASSERT(vertexMap.find(orientedReadId) == vertexMap.end());
    const vertex_descriptor v = add_vertex(LocalStringGraphVertex(
        orientedReadId, markerCount, isChimeric, distance), *this);
    vertexMap.insert(make_pair(orientedReadId, v));
}



void LocalStringGraph::addEdge(
    OrientedReadId from,
    OrientedReadId to,
    uint32_t overlapLen,
    uint32_t len,
    uint64_t globalArcId)
{
    const auto itFrom = vertexMap.find(from);
    DINARA_ASSERT(itFrom != vertexMap.end());
    const vertex_descriptor vFrom = itFrom->second;
    const auto itTo = vertexMap.find(to);
    DINARA_ASSERT(itTo != vertexMap.end());
    const vertex_descriptor vTo = itTo->second;

    add_edge(vFrom, vTo,
        LocalStringGraphEdge(overlapLen, len, globalArcId),
        *this);
}



bool LocalStringGraph::vertexExists(OrientedReadId orientedReadId) const
{
    return vertexMap.find(orientedReadId) != vertexMap.end();
}



uint32_t LocalStringGraph::getDistance(OrientedReadId orientedReadId) const
{
    const auto it = vertexMap.find(orientedReadId);
    DINARA_ASSERT(it != vertexMap.end());
    const vertex_descriptor v = it->second;
    return (*this)[v].distance;
}



ComputeLayoutReturnCode LocalStringGraph::computeLayout(
    const string& layoutMethod,
    double timeout)
{
    LocalStringGraph& graph = *this;
    std::map<vertex_descriptor, array<double, 2> > positionMap;
    const ComputeLayoutReturnCode returnCode =
        dinara::computeLayoutGraphviz(graph, layoutMethod, timeout, positionMap);
    if (returnCode != ComputeLayoutReturnCode::Success) {
        return returnCode;
    }

    BGL_FORALL_VERTICES(v, graph, LocalStringGraph) {
        const auto it = positionMap.find(v);
        DINARA_ASSERT(it != positionMap.end());
        graph[v].position = it->second;
    }
    return ComputeLayoutReturnCode::Success;
}



void LocalStringGraph::writeSvg(
    const string& svgId,
    uint64_t width,
    uint64_t height,
    double vertexScalingFactor,
    double edgeThicknessScalingFactor,
    uint64_t maxDistance,
    const Assembler& assembler,
    ostream& svg) const
{
    using Graph = LocalStringGraph;
    using VertexAttributes = WriteGraph::VertexAttributes;
    using EdgeAttributes = WriteGraph::EdgeAttributes;
    const Graph& graph = *this;

    std::map<vertex_descriptor, VertexAttributes> vertexAttributes;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        const auto& vertex = graph[v];
        const OrientedReadId orientedReadId = vertex.orientedReadId;
        VertexAttributes attributes;
        attributes.radius = vertexScalingFactor * sqrt(3.e-7 * double(vertex.markerCount));
        attributes.id = "Vertex-" + orientedReadId.getString();

        if (vertex.distance == 0) {
            attributes.color = "lime";
        } else if (vertex.distance == maxDistance) {
            attributes.color = "cyan";
        } else if (vertex.isChimeric) {
            attributes.color = "red";
        }

        attributes.tooltip =
            "Read " + orientedReadId.getString() + ", " + to_string(vertex.markerCount) +
            " markers, distance " + to_string(vertex.distance);
        attributes.url = "exploreRead?readId=" + to_string(orientedReadId.getReadId()) +
            "&strand=" + to_string(orientedReadId.getStrand());

        vertexAttributes.insert(make_pair(v, attributes));
    }

    std::map<edge_descriptor, EdgeAttributes> edgeAttributes;
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        const auto& edge = graph[e];
        const auto vFrom = source(e, graph);
        const auto vTo = target(e, graph);
        const OrientedReadId from = graph[vFrom].orientedReadId;
        const OrientedReadId to = graph[vTo].orientedReadId;

        EdgeAttributes attributes;
        attributes.thickness = edgeThicknessScalingFactor * (1.e-4 * double(edge.overlapLen));
        attributes.tooltip =
            "Arc " + from.getString() + " -> " + to.getString() +
            ", ol=" + to_string(edge.overlapLen) +
            ", len=" + to_string(edge.len) +
            ", arcId=" + to_string(edge.globalArcId) +
            ", alignmentId=" + to_string(assembler.stringGraph.arcs[edge.globalArcId].alignmentId);
        if (!edge.color.empty()) {
            attributes.color = edge.color;
        }
        edgeAttributes.insert(make_pair(e, attributes));
    }

    WriteGraph::writeSvg(graph, svgId, width, height, vertexAttributes, edgeAttributes, svg);
}
