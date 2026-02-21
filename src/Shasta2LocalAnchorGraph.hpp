#pragma once

#include "Shasta2Anchors.hpp"
#include "Shasta2AnchorGraph.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <map>
#include "vector.hpp"

namespace dinara {
    class Shasta2LocalAnchorGraph;
    class Shasta2LocalAnchorGraphVertex;
    class Shasta2LocalAnchorGraphEdge;
    class Shasta2LocalAnchorGraphDisplayOptions;
    class Shasta2AssemblyGraphPostprocessor;

    using Shasta2LocalAnchorGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Shasta2LocalAnchorGraphVertex,
        Shasta2LocalAnchorGraphEdge>;
}



class dinara::Shasta2LocalAnchorGraphDisplayOptions {
public:
    uint64_t sizePixels;
    string layoutMethod;

    // Vertices.
    double vertexSize;
    bool vertexSizeByCoverage;
    bool vertexLabels;
    string vertexColoring;
    string similarityMeasure;
    string referenceAnchorIdString;
    string assemblyStage;

    // Edges.
    string edgeColoring;
    double edgeThickness;
    double minimumEdgeLength;
    double additionalEdgeLengthPerKb;
    double arrowSize;
    bool edgeLabels;

    // Construct from an html request.
    Shasta2LocalAnchorGraphDisplayOptions(const vector<string>& request);

    // Write the form.
    void writeForm(ostream& html) const;
};



class dinara::Shasta2LocalAnchorGraphVertex {
public:
    Shasta2AnchorId anchorId;
    uint64_t distance;
    Shasta2LocalAnchorGraphVertex(
        Shasta2AnchorId anchorId,
        uint64_t distance) :
        anchorId(anchorId),
        distance(distance)
    {}
};



class dinara::Shasta2LocalAnchorGraphEdge {
public:
    // The edge of the global AnchorGraph that corresponds to this Shasta2LocalAnchorGraphEdge.
    Shasta2AnchorGraph::edge_descriptor eG;
};



class dinara::Shasta2LocalAnchorGraph : public Shasta2LocalAnchorGraphBaseClass {
public:

    Shasta2LocalAnchorGraph(
        const Shasta2Anchors&,
        const Shasta2AnchorGraph&,
        const vector<Shasta2AnchorId>&,
        uint64_t maxDistance,
        uint64_t minCoverage,
        bool edgesMarkedForAssembly);

    const Shasta2Anchors& anchors;
    const Shasta2AnchorGraph* anchorGraphPointer = nullptr;
    uint64_t maxDistance = 0;
    std::map<Shasta2AnchorId, vertex_descriptor> vertexMap;

    void writeHtml(
        ostream& html,
        const Shasta2LocalAnchorGraphDisplayOptions&,
        const Shasta2AssemblyGraphPostprocessor*);
    void writeHtml1(
        ostream& html,
        const Shasta2LocalAnchorGraphDisplayOptions&,
        const Shasta2AssemblyGraphPostprocessor*) const;

    void writeGraphviz(
        const string& fileName,
        const Shasta2LocalAnchorGraphDisplayOptions&,
        const Shasta2AssemblyGraphPostprocessor*) const;
    void writeGraphviz(
        ostream&,
        const Shasta2LocalAnchorGraphDisplayOptions&,
        const Shasta2AssemblyGraphPostprocessor*) const;

private:
    // Html/svg output without using svg output created by Graphviz.
    void writeHtml2(
        ostream& html,
        const Shasta2LocalAnchorGraphDisplayOptions&,
        const Shasta2AssemblyGraphPostprocessor*);

    // The position of each vertex in the computed layout.
    std::map<vertex_descriptor, array<double, 2> > layout;
    void computeLayout(const Shasta2LocalAnchorGraphDisplayOptions&);

    // The bounding box of the computed layout.
    class Box {
    public:
        double xMin;
        double xMax;
        double yMin;
        double yMax;
        double xSize() {return xMax - xMin;}
        double ySize() {return yMax - yMin;}
        void makeSquare();
        void extend(double factor);
    };
    Box boundingBox;
    void computeLayoutBoundingBox();

    void writeVertices(
        ostream& html,
        const Shasta2LocalAnchorGraphDisplayOptions&,
        const Shasta2AssemblyGraphPostprocessor*) const;
    void writeEdges(
        ostream& html,
        const Shasta2LocalAnchorGraphDisplayOptions&) const;
    void writeSvgControls(
        ostream& html,
        const Shasta2LocalAnchorGraphDisplayOptions&) const;
};

