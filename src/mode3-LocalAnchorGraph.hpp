#pragma once

// Dinara.
#include "mode3-Anchor.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library;
#include <map>
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

        class LocalAnchorGraphDisplayOptions;
    }
}



class dinara::mode3::LocalAnchorGraphDisplayOptions {
public:
    uint64_t sizePixels;
    string layoutMethod;
    double layoutTimeoutSeconds;

    // Optional graph processing for visualization/debugging.
    // raw: no additional processing.
    // cleaned: hide edges without valid offsets, normalize directions by offset sign.
    // reduced: cleaned + transitive reduction.
    // bubble: reduced + highlight detected bubble edges.
    string processingMode;
    uint64_t transitiveFuzzBases;
    double transitiveCoverageFactor;
    uint64_t bubbleMaxDepth;

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
    uint8_t displayFlags = 0;

    enum DisplayFlag : uint8_t {
        Hidden = 1,
        Transitive = 2,
        Bubble = 4,
        FlippedByOffset = 8
    };

    bool isHidden() const
    {
        return (displayFlags & Hidden) != 0;
    }

    double coverageLoss() const
    {
        // Robust against incomplete/invalid AnchorPairInfo.
        // In normal operation, coverage <= common and common > 0.
        if(info.common == 0) {
            return 1.;
        }
        if(coverage >= info.common) {
            return 0.;
        }
        return double(info.common - coverage) / double(info.common);
    }
};



class dinara::mode3::LocalAnchorGraph : public LocalAnchorGraphAssemblyGraphBaseClass {
public:
    LocalAnchorGraph(
        const Anchors&,
        const vector<AnchorId>&,
        uint64_t maxDistance,
        bool filterEdgesByCoverageLoss,
        double maxCoverageLoss,
        uint64_t minCoverage);

    const Anchors& anchors;
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
    void prepareForDisplay(const LocalAnchorGraphDisplayOptions&);
    void normalizeDirectionsByOffset();
    void reduceTransitiveEdges(uint64_t fuzzBases, double coverageFactor);
    void highlightBubbles(uint64_t maxDepth);

    uint64_t hiddenEdgeCount = 0;
    uint64_t transitiveEdgeCount = 0;
    uint64_t flippedEdgeCount = 0;
    uint64_t bubbleEdgeCount = 0;

    // Html/svg output without using svg output created by Graphviz.
    void writeHtml2(
        ostream& html,
        const LocalAnchorGraphDisplayOptions&);

    // The position of each vertex in the computed layout.
    std::map<vertex_descriptor, array<double, 2> > layout;
    void computeLayout(const LocalAnchorGraphDisplayOptions&);
    string layoutStatusMessage;

    // The bounding box of the computed layout.
    class Box {
    public:
        double xMin;
        double xMax;
        double yMin;
        double yMax;
        double xSize() const {return xMax - xMin;}
        double ySize() const {return yMax - yMin;}
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
