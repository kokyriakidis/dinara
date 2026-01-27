// Dinara.
#include "LocalUnitigGraph.hpp"
#include "Assembler.hpp"
#include "writeGraph.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include <cmath>



void LocalUnitigGraph::addVertex(
    OrientedUnitigId orientedUnitigId,
    uint32_t readCount,
    bool isCircular,
    uint32_t distance)
{
    DINARA_ASSERT(vertexMap.find(orientedUnitigId) == vertexMap.end());
    const vertex_descriptor v = add_vertex(LocalUnitigGraphVertex(
        orientedUnitigId, readCount, isCircular, distance), *this);
    vertexMap.insert(make_pair(orientedUnitigId, v));
}



void LocalUnitigGraph::addEdge(
    OrientedUnitigId from,
    OrientedUnitigId to,
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
        LocalUnitigGraphEdge(overlapLen, len, globalArcId),
        *this);
}



bool LocalUnitigGraph::vertexExists(OrientedUnitigId orientedUnitigId) const
{
    return vertexMap.find(orientedUnitigId) != vertexMap.end();
}



uint32_t LocalUnitigGraph::getDistance(OrientedUnitigId orientedUnitigId) const
{
    const auto it = vertexMap.find(orientedUnitigId);
    DINARA_ASSERT(it != vertexMap.end());
    const vertex_descriptor v = it->second;
    return (*this)[v].distance;
}



ComputeLayoutReturnCode LocalUnitigGraph::computeLayout(
    const string& layoutMethod,
    double timeout)
{
    LocalUnitigGraph& graph = *this;
    std::map<vertex_descriptor, array<double, 2> > positionMap;
    const ComputeLayoutReturnCode returnCode =
        dinara::computeLayoutGraphviz(graph, layoutMethod, timeout, positionMap);
    if (returnCode != ComputeLayoutReturnCode::Success) {
        return returnCode;
    }

    BGL_FORALL_VERTICES(v, graph, LocalUnitigGraph) {
        const auto it = positionMap.find(v);
        DINARA_ASSERT(it != positionMap.end());
        graph[v].position = it->second;
    }
    return ComputeLayoutReturnCode::Success;
}



void LocalUnitigGraph::writeSvg(
    const string& svgId,
    uint64_t width,
    uint64_t height,
    double vertexScalingFactor,
    double edgeThicknessScalingFactor,
    uint64_t maxDistance,
    const Assembler& assembler,
    ostream& svg) const
{
    using Graph = LocalUnitigGraph;
    using VertexAttributes = WriteGraph::VertexAttributes;
    using EdgeAttributes = WriteGraph::EdgeAttributes;
    const Graph& graph = *this;

    std::map<vertex_descriptor, VertexAttributes> vertexAttributes;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        const auto& vertex = graph[v];
        const OrientedUnitigId orientedUnitigId = vertex.orientedUnitigId;
        VertexAttributes attributes;
        attributes.radius = vertexScalingFactor * sqrt(0.05 * double(std::max<uint32_t>(vertex.readCount, 1U)));
        attributes.id = "Vertex-" + orientedUnitigId.getString();

        if (vertex.distance == 0) {
            attributes.color = "lime";
        } else if (vertex.distance == maxDistance) {
            attributes.color = "cyan";
        } else if (vertex.isCircular) {
            attributes.color = "orange";
        }

        attributes.tooltip =
            "Unitig " + orientedUnitigId.getString() +
            ", reads=" + to_string(vertex.readCount) +
            ", distance=" + to_string(vertex.distance) +
            (vertex.isCircular ? ", circular" : "");
        attributes.url = "exploreUnitigGraph?unitigId=" + orientedUnitigId.getString();

        vertexAttributes.insert(make_pair(v, attributes));
    }

    std::map<edge_descriptor, EdgeAttributes> edgeAttributes;
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        const auto& edge = graph[e];
        const auto vFrom = source(e, graph);
        const auto vTo = target(e, graph);
        const OrientedUnitigId from = graph[vFrom].orientedUnitigId;
        const OrientedUnitigId to = graph[vTo].orientedUnitigId;

        EdgeAttributes attributes;
        attributes.thickness = edgeThicknessScalingFactor * (1.e-4 * double(edge.overlapLen));
        attributes.tooltip =
            "Arc " + from.getString() + " -> " + to.getString() +
            ", ol=" + to_string(edge.overlapLen) +
            ", len=" + to_string(edge.len) +
            ", arcId=" + to_string(edge.globalArcId);
        if (!edge.color.empty()) {
            attributes.color = edge.color;
        }
        edgeAttributes.insert(make_pair(e, attributes));
    }

    WriteGraph::writeSvg(graph, svgId, width, height, vertexAttributes, edgeAttributes, svg);
}

