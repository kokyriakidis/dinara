#pragma once

// The read-anchor graph is an undirected bipartite graph
// in which each vertex represents an oriented read or an anchor.

// Shasta.
#include "Shasta2Anchors.hpp"
#include "invalid.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include <map>
#include "vector.hpp"
#include "string.hpp"



namespace dinara {
    class Shasta2LocalReadAnchorGraph;

    class Shasta2LocalReadAnchorGraphVertex;

    using Shasta2LocalReadAnchorGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::undirectedS,
        Shasta2LocalReadAnchorGraphVertex>;

    class Shasta2Journeys;
}



class dinara::Shasta2LocalReadAnchorGraphVertex
{
public:
    Shasta2AnchorId anchorId = invalid<Shasta2AnchorId>;
    OrientedReadId orientedReadId = OrientedReadId(invalid<ReadId>, 0);
    uint64_t distance = 0;

    // The position of this vertex in the computed layout.
    double x;
    double y;

    Shasta2LocalReadAnchorGraphVertex(Shasta2AnchorId anchorId, uint64_t distance) :
        anchorId(anchorId),
        distance(distance)
    {}

    Shasta2LocalReadAnchorGraphVertex(OrientedReadId orientedReadId, uint64_t distance) :
        orientedReadId(orientedReadId),
        distance(distance)
    {}

    bool isAnchor() const
    {
        return anchorId != invalid<Shasta2AnchorId>;
    }

    bool isOrientedRead() const
    {
        return not isAnchor();
    }


};



class dinara::Shasta2LocalReadAnchorGraph : public Shasta2LocalReadAnchorGraphBaseClass{
public:

    // The constructor parses the request, creates the Shasta2LocalReadAnchorGraph,
    // and displays it to html.
    Shasta2LocalReadAnchorGraph(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        const vector<string>& request,
        ostream& html);

private:
    const Shasta2Anchors& anchors;
    const Shasta2Journeys& journeys;
    ostream& html;

    // Figure out if command "customLayout" is available.
    bool customLayoutIsAvailable = false;
    void findOutIfCustomLayoutIsAvailable();

    // The request options.
    string startVerticesString;
    uint64_t maxDistance = 3;
    uint64_t sizePixels = 600;
    string layoutMethod = "sfdp";
    void getRequestOptions(const vector<string>& request);

    // Write the page header.
    void writeHeader();

    // Write the form to enter the options.
    void writeForm();

    // The AnchorIds and OrientedReadIds for the starting vertices.
    vector<Shasta2AnchorId> startAnchorIds;
    vector<OrientedReadId> startOrientedReadIds;
    void parseStartVertices();

    // Create the graph.
    void createVertices();
    void createEdges();
    std::map<Shasta2AnchorId, vertex_descriptor> anchorVertexMap;
    std::map<OrientedReadId, vertex_descriptor> orientedReadVertexMap;

    // Compute the layout of the graph.
    // This stores a layout position in each vertex.
    void computeLayout();

    // The bounding box of the layout.
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
    void computeBoundingBox();

    // Display the graph to html.
    void write() const;
    void writeVertices() const;
    void writeEdges() const;
    void writeSvgControls() const;

};
