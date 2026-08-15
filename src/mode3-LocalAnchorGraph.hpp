#pragma once

#include "mode3-Anchor.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library;
#include <map>
#include <unordered_map>
#include "vector.hpp"

namespace dinara {
    namespace mode3 {
        class LocalAnchorGraph;
        class LocalAnchorGraphVertex;
        class LocalAnchorGraphEdge;

        using LocalAnchorGraphAssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            LocalAnchorGraphVertex,
            LocalAnchorGraphEdge>;

        class AnchorGraph;
        class LocalAnchorGraphDisplayOptions;
    }
}



class dinara::mode3::LocalAnchorGraphDisplayOptions {
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

    // Edges.
    string edgeColoring;
    double edgeThickness;
    double minimumEdgeLength;
    double additionalEdgeLengthPerKb;
    double arrowSize;
    bool edgeLabels;

    // Construct from an html request.
    LocalAnchorGraphDisplayOptions(const vector<string>& request);

    // Write the form.
    void writeForm(ostream& html) const;
};



class dinara::mode3::LocalAnchorGraphVertex {
public:
    AnchorId anchorId;
    uint64_t distance;
    LocalAnchorGraphVertex(
        AnchorId anchorId,
        uint64_t distance) :
        anchorId(anchorId),
        distance(distance)
    {}
};


class dinara::mode3::LocalAnchorGraphEdge {
public:
    AnchorPairInfo info;
    uint64_t coverage;

    double coverageLoss() const
    {
        return double(info.common - coverage) / double(info.common);
    }
};



class dinara::mode3::LocalAnchorGraph : public LocalAnchorGraphAssemblyGraphBaseClass {
public:
    // Construct by BFS on the Anchors directly (existing behavior).
    LocalAnchorGraph(
        const Anchors&,
        const vector<AnchorId>&,
        uint64_t maxDistance,
        bool filterEdgesByCoverageLoss,
        double maxCoverageLoss,
        uint64_t minCoverage);

    // Construct by BFS on a pre-computed global AnchorGraph.
    LocalAnchorGraph(
        const Anchors&,
        const AnchorGraph&,
        const vector<AnchorId>&,
        uint64_t maxDistance,
        uint64_t minCoverage);

    const Anchors* anchorsPtr = nullptr;

    uint64_t maxDistance;
    std::map<AnchorId, vertex_descriptor> vertexMap;

    void writeHtml(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&);
    void writeHtml1(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&) const;


    void writeGraphviz(
        const string& fileName,
        const LocalAnchorGraphDisplayOptions&
        ) const;
    void writeGraphviz(
        ostream&,
        const LocalAnchorGraphDisplayOptions&
        ) const;

private:

    // Helper methods that dispatch based on data source.
    uint64_t getVertexCoverage(AnchorId anchorId) const;
    string getAnchorIdString(AnchorId anchorId) const;
    string getAnchorUrl(AnchorId anchorId) const;
    string getEdgeUrl(AnchorId anchorId0, AnchorId anchorId1) const;

    // Get similarity info between referenceAnchorId and anchorId.
    void getAnchorPairInfo(AnchorId referenceAnchorId, AnchorId anchorId, AnchorPairInfo& info) const;

    // Html/svg output without using svg output created by Graphviz.
    void writeHtml2(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&);

    // The position of each vertex in the computed layout.
    std::map<vertex_descriptor, array<double, 2> > layout;
    void computeLayout(const LocalAnchorGraphDisplayOptions&);

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
        const LocalAnchorGraphDisplayOptions&) const;

    void writeEdges(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&) const;

    void writeSvgControls(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&) const;
};