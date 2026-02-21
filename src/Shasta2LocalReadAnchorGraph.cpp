// Shasta.
#include "Shasta2LocalReadAnchorGraph.hpp"
#include "computeLayout.hpp"
#include "deduplicate.hpp"
#include "html.hpp"
#include "HttpServer.hpp"
#include "Shasta2Journeys.hpp"
#include "Reads.hpp"
#include "DINARA_ASSERT.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/tokenizer.hpp>

// Standard library.
#include "iostream.hpp"
#include <queue>



Shasta2LocalReadAnchorGraph::Shasta2LocalReadAnchorGraph(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const vector<string>& request,
    ostream& html) :
    anchors(anchors),
    journeys(journeys),
    html(html)
{
    Shasta2LocalReadAnchorGraph& graph = *this;

    findOutIfCustomLayoutIsAvailable();

    writeHeader();

    getRequestOptions(request);
    boost::trim(startVerticesString);

    writeForm();
    parseStartVertices();

    createVertices();
    createEdges();
    html << "<br>The local read anchor graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;
    html << "<br>There are " << orientedReadVertexMap.size() <<
        " vertices representing oriented reads and " <<
        anchorVertexMap.size() << " vertices representing anchors.\n";

    computeLayout();
    computeBoundingBox();
    write();
}



void Shasta2LocalReadAnchorGraph::findOutIfCustomLayoutIsAvailable()
{
    const int commandStatus = std::system("which customLayout > /dev/null");
    DINARA_ASSERT(WIFEXITED(commandStatus));
    const int returnCode = WEXITSTATUS(commandStatus);

    customLayoutIsAvailable = (returnCode == 0);
}



void Shasta2LocalReadAnchorGraph::getRequestOptions(const vector<string>& request)
{
    HttpServer::getParameterValue(request, "startVerticesString", startVerticesString);
    HttpServer::getParameterValue(request, "maxDistance", maxDistance);
    HttpServer::getParameterValue(request, "sizePixels", sizePixels);

    // The default layout is custom, if available, instead of sfdp.
    if(customLayoutIsAvailable) {
        layoutMethod = "custom";
    }
    HttpServer::getParameterValue(request, "layoutMethod", layoutMethod);
}



void Shasta2LocalReadAnchorGraph::writeForm()
{
    html <<
        "<br><br><form><table>"
        "<tr>"
        "<th class=left>Starting oriented read ids and/or anchor ids"
        "<td class=centered><input type=text name=startVerticesString style='text-align:center' required";
    if(not startVerticesString.empty()) {
        html << " value='" << startVerticesString + "'";
    }
    html <<
        " size=20>";

    html << "<tr>"
        "<th class=left>Distance"
        "<td class=centered>"
        "<input type=text name=maxDistance style='text-align:center' required size=8 value=" <<
        maxDistance << ">";

    html <<
        "<tr>"
        "<th class=left title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "Graphics size in pixels"
        "<td class=centered><input type=text required name=sizePixels size=8 style='text-align:center'" <<
        " value='" << sizePixels <<
        "'>";

    html <<
        "<tr>"
        "<th>Layout method"
        "<td class=left>"
        "<input type=radio required name=layoutMethod value='sfdp'" <<
        (layoutMethod == "sfdp" ? " checked=on" : "") <<
        ">sfdp"
        "<br><input type=radio required name=layoutMethod value='fdp'" <<
        (layoutMethod == "fdp" ? " checked=on" : "") <<
        ">fdp"
        "<br><input type=radio required name=layoutMethod value='neato'" <<
        (layoutMethod == "neato" ? " checked=on" : "") <<
        ">neato"
        "<br><input type=radio required name=layoutMethod value='dot'" <<
        (layoutMethod == "dot" ? " checked=on" : "") <<
        ">dot";

    // If command "customLayout" is available, add an option for that.
    if(customLayoutIsAvailable) {
        html <<
            "<br><input type=radio required name=layoutMethod value='custom'" <<
            (layoutMethod == "custom" ? " checked=on" : "") <<
            ">custom";
    }

    html <<
        "</table>"
        "<br>"
        "<input type=submit value='Create the local read-anchor graph'>"
        "</form>";

}



void Shasta2LocalReadAnchorGraph::writeHeader()
{
    html <<
        "<h2>Local read-anchor graph</h2>"
        "<p>The read-anchor graph is an undirected bipartite graph "
        "in which each vertex represents an oriented read or an anchor.";
}



void Shasta2LocalReadAnchorGraph::parseStartVertices()
{
    // Loop over tokens in the startVerticesString.
    boost::trim(startVerticesString);
    boost::tokenizer< boost::char_separator<char> > tokenizer(startVerticesString, boost::char_separator<char>(", "));
    for(const string& startVertexString: tokenizer) {

        // Try to interpret it as an Shasta2AnchorId.
        const Shasta2AnchorId anchorId = shasta2AnchorIdFromString(startVertexString);
        if((anchorId == invalid<Shasta2AnchorId>)) {

            // Not a valid Shasta2AnchorId. Try to interpret it as an OrientedReadId.
            const OrientedReadId orientedReadId = OrientedReadId(startVertexString);

            // The token represents an OrientedReadId, but we have to check that it is valid.
            if(orientedReadId.getValue() >= journeys.size()) {
                throw runtime_error(startVertexString + " is not a valid oriented read id for this assembly.");
            }
            startOrientedReadIds.push_back(orientedReadId);

        } else {

            // The token represents an Shasta2AnchorId, but we have to check that it is valid.
            if(anchorId >= anchors.size()) {
                throw runtime_error(startVertexString + " is not a valid anchor for this assembly.");
            }
            startAnchorIds.push_back(anchorId);
        }

    }

    // Sort and remove duplicates.
    deduplicate(startAnchorIds);
    deduplicate(startOrientedReadIds);
}



void Shasta2LocalReadAnchorGraph::createVertices()
{
    Shasta2LocalReadAnchorGraph& graph = *this;

    // BFS queue.
    std::queue<vertex_descriptor> q;

    // Create the vertices corresponding to the initial AnchorIds and add them to the queue.
    for(const Shasta2AnchorId anchorId: startAnchorIds) {
        const vertex_descriptor v = add_vertex(Shasta2LocalReadAnchorGraphVertex(anchorId, 0), graph);
        anchorVertexMap.insert(make_pair(anchorId, v));
        q.push(v);
    }

    // Create the vertices corresponding to the initial OrientedReadIds and add them to the queue.
    for(const OrientedReadId orientedReadId: startOrientedReadIds) {
        const vertex_descriptor v = add_vertex(Shasta2LocalReadAnchorGraphVertex(orientedReadId, 0), graph);
        orientedReadVertexMap.insert(make_pair(orientedReadId, v));
        q.push(v);
    }



    // BFS loop.
    while(not q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const Shasta2LocalReadAnchorGraphVertex& vertex0 = graph[v0];
        const uint64_t distance0 = vertex0.distance;
        const uint64_t distance1 = distance0 + 1;

        if(vertex0.isAnchor()) {

            // This vertex corresponds to an Shasta2AnchorId.
            // Loop over the AnchorMarkerInfos for this Shasta2Anchor.
            const Shasta2AnchorId anchorId0 = vertex0.anchorId;
            const Shasta2Anchor anchor0 = anchors[anchorId0];
            for(const Shasta2AnchorMarkerInfo& anchorMarkerInfo1: anchor0) {
                const OrientedReadId orientedReadId1 = anchorMarkerInfo1.orientedReadId;
                if(not orientedReadVertexMap.contains(orientedReadId1)) {

                    // We don't have a vertex for this OrientedReadId. Add it.
                    const vertex_descriptor v1 =
                        add_vertex(Shasta2LocalReadAnchorGraphVertex(orientedReadId1, distance1), graph);
                    orientedReadVertexMap.insert(make_pair(orientedReadId1, v1));
                    if(distance1 < maxDistance) {
                        q.push(v1);
                    }

                }
            }

        } else {

            // This vertex corresponds to an OrientedReadId.
            // Loop over its journey.
            const OrientedReadId orientedReadId0 = vertex0.orientedReadId;
            const Shasta2Journey journey0 = journeys[orientedReadId0];
            for(const Shasta2AnchorId anchorId1: journey0) {
                if(not anchorVertexMap.contains(anchorId1)) {

                    // We don't have a vertex for this Shasta2AnchorId. Add it.
                    const vertex_descriptor v1 =
                        add_vertex(Shasta2LocalReadAnchorGraphVertex(anchorId1, distance1), graph);
                    anchorVertexMap.insert(make_pair(anchorId1, v1));
                    if(distance1 < maxDistance) {
                        q.push(v1);
                    }

                }
            }
        }

    }
}



void Shasta2LocalReadAnchorGraph::createEdges()
{
    Shasta2LocalReadAnchorGraph& graph = *this;

    for(const auto& p: orientedReadVertexMap) {
        const OrientedReadId orientedReadId0 = p.first;
        const vertex_descriptor v0 = p.second;

        const Shasta2Journey journey0 = journeys[orientedReadId0];
        for(const Shasta2AnchorId anchorId1: journey0) {
            const auto it1 = anchorVertexMap.find(anchorId1);
            if(it1 != anchorVertexMap.end()) {
                const vertex_descriptor v1 = it1->second;
                add_edge(v0, v1, graph);
            }
        }
    }
}



void Shasta2LocalReadAnchorGraph::computeLayout()
{
    Shasta2LocalReadAnchorGraph& graph = *this;

    // Create a map containing the desired length for each edge.
    std::map<edge_descriptor, double> edgeLengthMap;
    BGL_FORALL_EDGES(e, graph, Shasta2LocalReadAnchorGraph) {
        edgeLengthMap.insert(make_pair(e, 1));
    }



    // Compute the graph layout.
    std::map<vertex_descriptor, array<double, 2> > layout;
    const double timeout = 60.;
    if(layoutMethod == "custom") {
        const int quality = 2;
        computeLayoutCustom(
            graph,
            edgeLengthMap,
            layout,
            quality,
            timeout);
    } else {
        const string additionalOptions = "";
        computeLayoutGraphviz(
            graph,
            layoutMethod,
            timeout,
            layout,
            additionalOptions,
            &edgeLengthMap);
    }



    // Store the layout in the vertices.
    BGL_FORALL_VERTICES(v, graph, Shasta2LocalReadAnchorGraph) {
        Shasta2LocalReadAnchorGraphVertex& vertex = graph[v];
        const array<double, 2>& position = layout[v];
        vertex.x = position[0];
        vertex.y = position[1];
    }
}



void Shasta2LocalReadAnchorGraph::Box::extend(double factor)
{
    const double extend = factor * max(xSize(), ySize());
    xMin -= extend;
    xMax += extend;
    yMin -= extend;
    yMax += extend;
}



void Shasta2LocalReadAnchorGraph::Box::makeSquare()
{
    if(xSize() > ySize()) {
        const double delta = (xSize() - ySize()) / 2.;
        yMin -= delta;
        yMax += delta;
    } else {
        const double delta = (ySize() - xSize()) / 2.;
        xMin -= delta;
        xMax += delta;
    }
}



void Shasta2LocalReadAnchorGraph::computeBoundingBox()
{
    Shasta2LocalReadAnchorGraph& graph = *this;

    boundingBox.xMin = std::numeric_limits<double>::max();
    boundingBox.xMax = std::numeric_limits<double>::min();
    boundingBox.yMin = std::numeric_limits<double>::max();
    boundingBox.yMax = std::numeric_limits<double>::min();
    BGL_FORALL_VERTICES(v, graph, Shasta2LocalReadAnchorGraph) {
        Shasta2LocalReadAnchorGraphVertex& vertex = graph[v];
        const double x = vertex.x;
        const double y = vertex.y;
        boundingBox.xMin = min(boundingBox.xMin, x);
        boundingBox.xMax = max(boundingBox.xMax, x);
        boundingBox.yMin = min(boundingBox.yMin, y);
        boundingBox.yMax = max(boundingBox.yMax, y);
    }
}



void Shasta2LocalReadAnchorGraph::write() const
{
    // Use scientific notation because svg does not accept floating points
    // ending with a decimal point.
    html << std::scientific;

    // Define the viewport box.
    Box viewportBox = boundingBox;
    viewportBox.extend(0.05);
    viewportBox.makeSquare();

    // Begin the svg.
    const string svgId = "Shasta2LocalReadAnchorGraph";
    html <<
        "\n<br><br><div style='display:inline-block;vertical-align:top;'>"
        "<svg id='" << svgId <<
        "' width='" <<  sizePixels <<
        "' height='" << sizePixels <<
        "' viewbox='" << viewportBox.xMin << " " << viewportBox.yMin << " " <<
        viewportBox.xSize() << " " <<
        viewportBox.ySize() << "'"
        " style='background-color:#f0f0f0'"
        ">\n";

    // Write the edges first to avoid obscuring the vertices.
    writeEdges();
    writeVertices();

    // Finish the svg.
    html << "</svg></div>";
    writeSvgControls();
}



void Shasta2LocalReadAnchorGraph::writeSvgControls() const
{
    // Add drag and zoom.
    addSvgDragAndZoom(html);

    // Side panel.
    html << "<div style='display:inline-block;margin-left:20px'>"
        "<p><table>";

    // Buttons to change vertex size.
    html << R"stringDelimiter(
    <tr><th class=left>Shasta2Anchor vertex size<td>
    <button type='button' onClick='changeVertexSize("anchorVertices", 0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeVertexSize("anchorVertices", 0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeVertexSize("anchorVertices", 0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeVertexSize("anchorVertices", 1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeVertexSize("anchorVertices", 2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeVertexSize("anchorVertices", 10.)' style='width:3em'>+++</button>

    <tr><th class=left>Oriented read vertex size<td>
    <button type='button' onClick='changeVertexSize("orientedReadVertices", 0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeVertexSize("orientedReadVertices", 0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeVertexSize("orientedReadVertices", 0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeVertexSize("orientedReadVertices", 1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeVertexSize("orientedReadVertices", 2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeVertexSize("orientedReadVertices", 10.)' style='width:3em'>+++</button>

    <tr><th class=left>Edge thickness<td>
    <button type='button' onClick='changeThickness(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeThickness(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeThickness(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeThickness(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeThickness(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeThickness(10.)' style='width:3em'>+++</button>

    <tr><th class=left>Zoom<td>
    <button type='button' onClick='zoomSvg(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='zoomSvg(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='zoomSvg(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='zoomSvg(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='zoomSvg(2.)' style='width:3em'>++</button>
    <button type='button' onClick='zoomSvg(10.)' style='width:3em'>+++</button>

    <script>
        function changeVertexSize(group, factor)
        {
            var vertexGroup = document.getElementById(group);
            var vertices = vertexGroup.getElementsByTagName('circle');
            for(i=0; i<vertices.length; i++) {
                v = vertices[i];
                v.setAttribute('r', factor * v.getAttribute('r'));
            }
        }
        function changeThickness(factor)
        {
            var edgeGroup = document.getElementById('edges');
            numericWidth = parseFloat(edgeGroup.style.strokeWidth);
            edgeGroup.style.strokeWidth = numericWidth * factor + 'px';
        }
    </script>
    )stringDelimiter";


    // End of side panel.
    html << "</table>";
    html << "</div>";

    // Scroll down to the svg.
    const string svgId = "Shasta2LocalReadAnchorGraph";
    html <<
        "<script>"
        "document.getElementById('" << svgId << "').scrollIntoView({block:'center'});"
        "</script>";

    html <<
        "<br>Use Ctrl+Click to pan."
        "<br>Use Ctrl-Wheel or the above buttons to zoom.";
}



void Shasta2LocalReadAnchorGraph::writeVertices() const
{
    const Shasta2LocalReadAnchorGraph& graph = *this;


    // There are many more anchor vertices than read vertices.
    // Write the anchor vertices first, to avoid obscuring the
    // read vertices.

    html << "\n<g id='anchorVertices' style='stroke:none'>";
    const double anchorVertexRadius = 0.02;
    const string anchorVertexColor = "Red";
    BGL_FORALL_VERTICES(v, graph, Shasta2LocalReadAnchorGraph) {
        const Shasta2LocalReadAnchorGraphVertex& vertex = graph[v];
        if(vertex.isAnchor()) {
            html << "\n<circle cx='" << vertex.x << "' cy='" << vertex.y <<
                "' fill='" << anchorVertexColor <<
                "' r='" << anchorVertexRadius << "'>"
                "<title>" <<
                shasta2AnchorIdToString(vertex.anchorId) <<
                "</title>"
                "</circle>";
        }
    }
    html << "\n</g>";

    html << "\n<g id='orientedReadVertices' style='stroke:none'>";
    const double orientedReadVertexRadius = 0.05;
    const string orientedReadVertexColor = "Green";
    BGL_FORALL_VERTICES(v, graph, Shasta2LocalReadAnchorGraph) {
        const Shasta2LocalReadAnchorGraphVertex& vertex = graph[v];
        if(vertex.isOrientedRead()) {
            html << "\n<circle cx='" << vertex.x << "' cy='" << vertex.y <<
                "' fill='" << orientedReadVertexColor <<
                "' r='" << orientedReadVertexRadius << "' >"
                "<title>" <<
                vertex.orientedReadId <<
                "</title>"
                "</circle>";
        }
    }
    html << "\n</g>";

}



void Shasta2LocalReadAnchorGraph::writeEdges() const
{
    const Shasta2LocalReadAnchorGraph& graph = *this;

    html << "\n<g id='edges' style='stroke:Black;stroke-width:0.005'>";

    BGL_FORALL_EDGES(e, graph, Shasta2LocalReadAnchorGraph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Shasta2LocalReadAnchorGraphVertex& vertex0 = graph[v0];
        const Shasta2LocalReadAnchorGraphVertex& vertex1 = graph[v1];

        html <<
            "\n<line "
            "x1='" << vertex0.x << "' y1='" << vertex0.y <<
            "' x2='" << vertex1.x << "' y2='" << vertex1.y <<
            "' />";
    }

    html << "\n</g>";
}
