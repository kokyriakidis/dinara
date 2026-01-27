// Dinara.
#include "mode3-LocalAnchorGraph.hpp"
#include "computeLayout.hpp"
#include "html.hpp"
#include "HttpServer.hpp"
#include "MurmurHash2.hpp"
#include "platformDependent.hpp"
#include "runCommandWithTimeout.hpp"
using namespace dinara;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"
#include <cmath>
#include <algorithm>
#include <queue>
#include <set>



namespace {
    std::string computeLayoutReturnCodeName(dinara::ComputeLayoutReturnCode code)
    {
        using dinara::ComputeLayoutReturnCode;
        switch(code) {
        case ComputeLayoutReturnCode::Success:
            return "Success";
        case ComputeLayoutReturnCode::Error:
            return "Error";
        case ComputeLayoutReturnCode::Timeout:
            return "Timeout";
        case ComputeLayoutReturnCode::Signal:
            return "Signal";
        default:
            return "Unknown";
        }
    }
}


LocalAnchorGraph::LocalAnchorGraph(
    const Anchors& anchors,
    const vector<AnchorId>& anchorIds,
    uint64_t maxDistance,
    bool filterEdgesByCoverageLoss,
    double maxCoverageLoss,
    uint64_t minCoverage) :
    anchors(anchors),
    maxDistance(maxDistance)
{
    LocalAnchorGraph& graph = *this;

    // Initialize a BFS from these AnchorIds.
    std::queue<vertex_descriptor> q;
    for(const AnchorId anchorId: anchorIds) {
        DINARA_ASSERT(not vertexMap.contains(anchorId));
        const vertex_descriptor v = boost::add_vertex(LocalAnchorGraphVertex(anchorId, 0), graph);
        vertexMap.insert({anchorId, v});
        q.push(v);
    }

    // BFS to find the vertices. We will add the edges later.
    vector<AnchorId> neighbors;
    vector<uint64_t> coverage;
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        const LocalAnchorGraphVertex& vertex0 = graph[v0];
        const AnchorId anchorId0 = vertex0.anchorId;
        const uint64_t distance0 = vertex0.distance;
        const uint64_t distance1 = distance0 + 1;

        anchors.findChildren(anchorId0, neighbors, coverage, minCoverage);
        for(uint64_t i=0; i<neighbors.size(); i++) {
            const AnchorId anchorId1 = neighbors[i];
            auto it1 = vertexMap.find(anchorId1);
            if(it1 != vertexMap.end()) {
                continue;
            }

            // Filter by coverage loss, if requested.
            if(filterEdgesByCoverageLoss) {
                LocalAnchorGraphEdge edge;
                edge.coverage = coverage[i];
                anchors.analyzeAnchorPair(anchorId0, anchorId1, edge.info);
                if(edge.coverageLoss() > maxCoverageLoss) {
                    continue;
                }
            }

            const vertex_descriptor v1 = boost::add_vertex(LocalAnchorGraphVertex(anchorId1, distance1), graph);
            vertexMap.insert({anchorId1, v1});
            if(distance1 < maxDistance) {
                q.push(v1);
            }
        }

        anchors.findParents(anchorId0, neighbors, coverage, minCoverage);
        for(uint64_t i=0; i<neighbors.size(); i++) {
            const AnchorId anchorId1 = neighbors[i];
            auto it1 = vertexMap.find(anchorId1);
            if(it1 != vertexMap.end()) {
                continue;
            }

            // Filter by coverage loss, if requested.
            if(filterEdgesByCoverageLoss) {
                LocalAnchorGraphEdge edge;
                edge.coverage = coverage[i];
                anchors.analyzeAnchorPair(anchorId0, anchorId1, edge.info);
                if(edge.coverageLoss() > maxCoverageLoss) {
                    continue;
                }
            }

            const vertex_descriptor v1 = boost::add_vertex(LocalAnchorGraphVertex(anchorId1, distance1), graph);
            vertexMap.insert({anchorId1, v1});
            if(distance1 < maxDistance) {
                q.push(v1);
            }
        }
    }



    // Now add the edges.
    BGL_FORALL_VERTICES(v0, graph, LocalAnchorGraph) {
        const AnchorId anchorId0 = graph[v0].anchorId;
        anchors.findChildren(anchorId0, neighbors, coverage);
        for(uint64_t i=0; i<neighbors.size(); i++) {
            if(coverage[i] < minCoverage) {
                continue;
            }
            const AnchorId& anchorId1 = neighbors[i];
            auto it1 = vertexMap.find(anchorId1);
            if(it1 == vertexMap.end()) {
                continue;
            }
            const vertex_descriptor v1 = it1->second;

            // Create the tentative edge.
            LocalAnchorGraphEdge edge;
            edge.coverage = coverage[i];
            anchors.analyzeAnchorPair(anchorId0, anchorId1, edge.info);

            // Add it if requested.
            if((not filterEdgesByCoverageLoss) or
                (edge.coverageLoss() <= maxCoverageLoss)) {
                edge_descriptor e;
                tie(e, ignore) = add_edge(v0, v1, edge, graph);
            }
        }
    }
}



void LocalAnchorGraph::writeGraphviz(
    const string& fileName,
    const LocalAnchorGraphDisplayOptions& options) const
{
    ofstream file(fileName);
    writeGraphviz(file, options);
}



void LocalAnchorGraph::writeGraphviz(
    ostream& s,
    const LocalAnchorGraphDisplayOptions& options) const
{
    const LocalAnchorGraph& graph = *this;

    AnchorId referenceAnchorId = invalid<AnchorId>;
    if(options.vertexColoring == "byReadComposition") {
        referenceAnchorId = anchorIdFromString(options.referenceAnchorIdString);
        if((referenceAnchorId == invalid<AnchorId>) or (referenceAnchorId >= anchors.size())) {
            throw runtime_error("Invalid reference anchor id " + options.referenceAnchorIdString +
                ". Must be a number between 0 and " +
                to_string(anchors.size() / 2 - 1) + " followed by + or -.");
        }
    }
    const uint64_t referenceAnchorIdCoverage = anchors[referenceAnchorId].coverage();

    s << "digraph LocalAnchorGraph {\n";

    // Write the vertices.
    BGL_FORALL_VERTICES(v, graph, LocalAnchorGraph) {
        const LocalAnchorGraphVertex& vertex = graph[v];
        const AnchorId anchorId = vertex.anchorId;
        const string anchorIdString = anchorIdToString(anchorId);
        const uint64_t coverage = anchors[anchorId].coverage();

        // Vertex name.
        s << "\"" << anchorIdString << "\"";

        // Begin vertex attributes.
        s << "[";

        // URL
        s << "URL=\"exploreAnchor?anchorIdString=" << HttpServer::urlEncode(anchorIdString) << "\"";

        // Tooltip.
        s << " tooltip=\"" << anchorIdString << " " << coverage << "\"";

        // Label.
        if(options.vertexLabels) {
            s << " label=\"" << anchorIdString << "\\n" << coverage << "\"";
        }



        // Color.
        if(vertex.distance == 0) {
            s << " color=blue";
        } else if(vertex.distance == maxDistance) {
            s << " color=cyan";
        } else {

            // Color by similarity of read composition with the reference Anchor.
            if(options.vertexColoring == "byReadComposition") {
                AnchorPairInfo info;
                anchors.analyzeAnchorPair(referenceAnchorId, anchorId, info);

                double hue = 1.;    // 0=red, 1=green.
                if(options.similarityMeasure == "commonCount") {
                    // By common count.
                    hue = double(info.common) / double(referenceAnchorIdCoverage);

                } else if(options.similarityMeasure == "jaccard") {
                    // By Jaccard similarity.
                    hue = info.jaccard();
                } else {
                    // By corrected Jaccard similarity.
                    hue = info.correctedJaccard();
                 }

                const string colorString = "\"" + to_string(hue / 3.) + " 1. 1.\"";
                if(options.vertexLabels) {
                    s << " style=filled fillcolor=" << colorString;
                } else {
                    s << " color=" << colorString;
                    s << " fillcolor=" << colorString;
                }
             }
        }



        // Size.
        if(not options.vertexLabels) {
            const double displaySize =
                (options.vertexSizeByCoverage ?
                options.vertexSize * sqrt(0.1 * double(coverage)) :
                options.vertexSize
                ) / 72.;
            s << " width=" << displaySize ;
            s << " penwidth=" << 0.5 * displaySize;
        }

        // End vertex attributes.
        s << "]";

        // End the line for this vertex.
        s << ";\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, graph, LocalAnchorGraph) {
        if(graph[e].isHidden()) {
            continue;
        }
        const LocalAnchorGraphEdge& edge = graph[e];
        const double loss = edge.coverageLoss();

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        const LocalAnchorGraphVertex& vertex0 = graph[v0];
        const LocalAnchorGraphVertex& vertex1 = graph[v1];

        const AnchorId anchorId0 = vertex0.anchorId;
        const AnchorId anchorId1 = vertex1.anchorId;

        const string anchorId0String = anchorIdToString(anchorId0);
        const string anchorId1String = anchorIdToString(anchorId1);

        s << "\"" << anchorId0String << "\"->";
        s << "\"" << anchorId1String << "\"";

        // Begin edge attributes.
        s << " [";

        // URL
        s << "URL=\"exploreAnchorPair?"
            "anchorIdAString=" << HttpServer::urlEncode(anchorId0String) << "&"
            "anchorIdBString=" << HttpServer::urlEncode(anchorId1String) << "\"";

        // Tooltip.
        s << " tooltip="
            "\"" << anchorId0String << " to "
            << anchorId1String <<
            " " << edge.coverage << "/" << edge.info.common <<
            " loss " << std::fixed << std::setprecision(2) << loss <<
            " offset " << edge.info.offsetInBases << "\"";

        // Label.
        if(options.edgeLabels) {
            s << " label=\"" <<
                edge.coverage << "/" << edge.info.common <<
                "\\nLoss " << std::fixed << std::setprecision(2) << loss <<
                "\\nOffset " << edge.info.offsetInBases << "\"";
        }

        // Color.
        if(options.edgeColoring == "byCoverageLoss") {
            const double hue = (1. - loss) / 3.;
            s << " color=\"" << std::fixed << std::setprecision(2) << hue << " 1. 1.\"";
        } else if(options.edgeColoring == "random") {
            // To decide the color, hash the AnchorIds.
            // This way we always get the same color for the same edge.
            const auto p = make_pair(anchorId0, anchorId1);
            const uint32_t hashValue = MurmurHash2(&p, sizeof(p), 759);
            const double hue = double(hashValue % 360) / 360.;
            s << " color=\"" << std::fixed << std::setprecision(2) << hue << " 1. 1.\"";
        }

        // Thickness.
        s << " penwidth=" << 0.5 * options.edgeThickness * double(edge.coverage);

        // Arrow size.
        s << " arrowsize=" << 0.5 * options.arrowSize;

        // Length. Only use by fdp and neato layouts.
        int64_t offsetInBases = edge.info.offsetInBases;
        if(offsetInBases < 0) {
            offsetInBases = 10;
        }
        const double displayLength =
            (options.minimumEdgeLength +
                options.additionalEdgeLengthPerKb * 0.001 * double(offsetInBases)) / 72.;
        s << " len=" << displayLength;



        // End edge attributes.
        s << "]";

        // End the line for this edge.
        s << ";\n";
    }


    s << "}\n";
}



LocalAnchorGraphDisplayOptions::LocalAnchorGraphDisplayOptions(const vector<string>& request)
{
    // Figure out if command "customLayout" is available.
    const int commandStatus = std::system("which customLayout > /dev/null");
    DINARA_ASSERT(WIFEXITED(commandStatus));
    const int returnCode = WEXITSTATUS(commandStatus);
    const bool customLayoutIsAvailable = (returnCode == 0);


    sizePixels = 600;
    HttpServer::getParameterValue(request, "sizePixels", sizePixels);

    layoutMethod = (customLayoutIsAvailable ? "custom" : "sfdp");
    HttpServer::getParameterValue(request, "layoutMethod", layoutMethod);

    layoutTimeoutSeconds = 30.;
    HttpServer::getParameterValue(request, "layoutTimeoutSeconds", layoutTimeoutSeconds);

    processingMode = "raw";
    HttpServer::getParameterValue(request, "processingMode", processingMode);
    transitiveFuzzBases = 1000;
    HttpServer::getParameterValue(request, "transitiveFuzzBases", transitiveFuzzBases);
    transitiveCoverageFactor = 0.8;
    HttpServer::getParameterValue(request, "transitiveCoverageFactor", transitiveCoverageFactor);
    bubbleMaxDepth = 50;
    HttpServer::getParameterValue(request, "bubbleMaxDepth", bubbleMaxDepth);

    vertexColoring = "black";
    HttpServer::getParameterValue(request, "vertexColoring", vertexColoring);

    similarityMeasure = "commonCount";
    HttpServer::getParameterValue(request, "similarityMeasure", similarityMeasure);

    referenceAnchorIdString = "";
    HttpServer::getParameterValue(request, "referenceAnchorId", referenceAnchorIdString);

    edgeColoring = "random";
    HttpServer::getParameterValue(request, "edgeColoring", edgeColoring);

    vertexSize =  1.;
    HttpServer::getParameterValue(request, "vertexSize", vertexSize);

    string vertexSizeByCoverageString;
    vertexSizeByCoverage = HttpServer::getParameterValue(request,
        "vertexSizeByCoverage", vertexSizeByCoverageString);

    string vertexLabelsString;
    vertexLabels = HttpServer::getParameterValue(request,
        "vertexLabels", vertexLabelsString);

    minimumEdgeLength = 1.;
    HttpServer::getParameterValue(request, "minimumEdgeLength", minimumEdgeLength);

    additionalEdgeLengthPerKb = 1.;
    HttpServer::getParameterValue(request, "additionalEdgeLengthPerKb", additionalEdgeLengthPerKb);

    edgeThickness = 1.;
    HttpServer::getParameterValue(request, "edgeThickness", edgeThickness);

    arrowSize = 1.;
    HttpServer::getParameterValue(request, "arrowSize", arrowSize);

    string edgeLabelsString;
    edgeLabels = HttpServer::getParameterValue(request,
        "edgeLabels", edgeLabelsString);
}



void LocalAnchorGraphDisplayOptions::writeForm(ostream& html) const
{
    // Figure out if command "customLayout" is available.
    const int commandStatus = std::system("which customLayout > /dev/null");
    DINARA_ASSERT(WIFEXITED(commandStatus));
    const int returnCode = WEXITSTATUS(commandStatus);
    const bool customLayoutIsAvailable = (returnCode == 0);

    html <<
        "<tr>"
        "<th title='Graphics size in pixels. "
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
        "<tr>"
        "<th>Processing"
        "<td class=left>"
        "<input type=radio required name=processingMode value='raw'" <<
        (processingMode == "raw" ? " checked=on" : "") << ">raw"
        "<br><input type=radio required name=processingMode value='cleaned'" <<
        (processingMode == "cleaned" ? " checked=on" : "") << ">cleaned (hide edges without offsets; normalize directions)"
        "<br><input type=radio required name=processingMode value='reduced'" <<
        (processingMode == "reduced" ? " checked=on" : "") <<
        ">reduced (cleaned + transitive reduction)"
        "<br><input type=radio required name=processingMode value='bubble'" <<
        (processingMode == "bubble" ? " checked=on" : "") <<
        ">bubble (reduced + highlight bubbles)"
        "<hr>"
        "<input type=text name=transitiveFuzzBases style='text-align:center' required size=8 value=" <<
        transitiveFuzzBases << "> Transitive fuzz (bases)"
        "<br><input type=text name=transitiveCoverageFactor style='text-align:center' required size=6 value=" <<
        transitiveCoverageFactor << "> Transitive support factor (0..1)"
        "<br><input type=text name=bubbleMaxDepth style='text-align:center' required size=8 value=" <<
        bubbleMaxDepth << "> Bubble search depth";

    html <<
        "<tr>"
        "<th title='Maximum time allowed for graph layout. "
        "If layout times out or fails, a deterministic fallback layout is used.'>"
        "Layout timeout (seconds)"
        "<td class=centered><input type=text required name=layoutTimeoutSeconds size=8 style='text-align:center'" <<
        " value='" << layoutTimeoutSeconds <<
        "'>";

    html <<
        "<tr>"
        "<th>Vertices"
        "<td class=left>"
        "<input type=text name=vertexSize style='text-align:center' required size=6 value=" <<
        vertexSize << "> Vertex size (arbitrary units)"
        "<br><input type=checkbox name=vertexSizeByCoverage" <<
        (vertexSizeByCoverage ? " checked" : "") <<
        "> Size proportional to coverage"

        "<hr>"
        "<input type=checkbox name=vertexLabels" <<
        (vertexLabels ? " checked" : "") <<
        "> Labels (dot layout only)"

        "<hr>"
        "<b>Vertex coloring</b>"

        "<br><input type=radio required name=vertexColoring value='black'" <<
        (vertexColoring == "black" ? " checked=on" : "") << ">Black"
        "<br><input type=radio required name=vertexColoring value='byReadComposition'" <<
        (vertexColoring == "byReadComposition" ? " checked=on" : "") <<
        "> By similarity of read composition using similarity measure:"

        "<div style='padding-left:50px'>"
        "<input type=radio required name=similarityMeasure value='commonCount'" <<
        (similarityMeasure == "commonCount" ? " checked=on" : "") << ">Number of common oriented reads"
        "<br><input type=radio required name=similarityMeasure value='jaccard'" <<
        (similarityMeasure == "jaccard" ? " checked=on" : "") << ">Jaccard similarity"
        "<br><input type=radio required name=similarityMeasure value='correctedJaccard'" <<
        (similarityMeasure == "correctedJaccard" ? " checked=on" : "") << ">Corrected Jaccard similarity"
        "</div>"

        "<input type=text name=referenceAnchorId size=6 style='text-align:center'";
        if(not referenceAnchorIdString.empty()) {
            html << " value='" << referenceAnchorIdString + "'";
        }
        html << "> Reference anchor id";



    html <<
        "<tr>"
        "<th>Edges"
        "<td class=left>"

        "<b>Edge coloring</b>"
        "<br><input type=radio required name=edgeColoring value='black'" <<
        (edgeColoring == "black" ? " checked=on" : "") << ">Black"
        "<br><input type=radio required name=edgeColoring value='random'" <<
        (edgeColoring == "random" ? " checked=on" : "") << ">Random"
        "<br><input type=radio required name=edgeColoring value='byCoverageLoss'" <<
        (edgeColoring == "byCoverageLoss" ? " checked=on" : "") << ">By coverage loss"
        "<hr>"

        "<b>Edge graphics</b>"

        "<br><input type=text name=edgeThickness style='text-align:center' required size=6 value=" <<
        edgeThickness << "> Thickness (arbitrary units)"

        "<br><input type=text name=minimumEdgeLength style='text-align:center' required size=6 value=" <<
        minimumEdgeLength << "> Minimum edge length (arbitrary units)"

        "<br><input type=text name=additionalEdgeLengthPerKb style='text-align:center' required size=6 value=" <<
        additionalEdgeLengthPerKb << "> Additional edge length per Kb (arbitrary units)"

        "<br><input type=text name=arrowSize style='text-align:center' required size=6 value=" <<
        arrowSize << "> Arrow size (arbitrary units)"

        "<hr>"
        "<input type=checkbox name=edgeLabels" <<
        (edgeLabels ? " checked" : "") <<
        "> Labels (dot layout only)";

}



void LocalAnchorGraph::writeHtml(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& options)
{
    prepareForDisplay(options);

    if((options.layoutMethod == "dot") and (options.vertexLabels or options.edgeLabels)) {

        // Use svg output from graphviz.
        writeHtml1(html, options);

    } else {

        // Compute graph layout and use it to generate svg.
        writeHtml2(html, options);

    }
}



void LocalAnchorGraph::prepareForDisplay(const LocalAnchorGraphDisplayOptions& options)
{
    hiddenEdgeCount = 0;
    transitiveEdgeCount = 0;
    flippedEdgeCount = 0;
    bubbleEdgeCount = 0;

    BGL_FORALL_EDGES(e, *this, LocalAnchorGraph) {
        (*this)[e].displayFlags = 0;
    }

    if(options.processingMode == "raw") {
        return;
    }

    // Hide edges that don't have a valid offset (common==0).
    BGL_FORALL_EDGES(e, *this, LocalAnchorGraph) {
        auto& edge = (*this)[e];
        if(edge.info.common == 0) {
            edge.displayFlags |= LocalAnchorGraphEdge::Hidden;
        }
    }

    // Normalize directions: make all active edges have non-negative offset.
    normalizeDirectionsByOffset();

    if((options.processingMode == "reduced") or (options.processingMode == "bubble")) {
        reduceTransitiveEdges(options.transitiveFuzzBases, options.transitiveCoverageFactor);
    }
    if(options.processingMode == "bubble") {
        highlightBubbles(options.bubbleMaxDepth);
    }

    BGL_FORALL_EDGES(e, *this, LocalAnchorGraph) {
        const auto& edge = (*this)[e];
        if(edge.isHidden()) {
            ++hiddenEdgeCount;
        }
        if((edge.displayFlags & LocalAnchorGraphEdge::Transitive) != 0) {
            ++transitiveEdgeCount;
        }
        if((edge.displayFlags & LocalAnchorGraphEdge::FlippedByOffset) != 0) {
            ++flippedEdgeCount;
        }
        if((edge.displayFlags & LocalAnchorGraphEdge::Bubble) != 0) {
            ++bubbleEdgeCount;
        }
    }
}



void LocalAnchorGraph::normalizeDirectionsByOffset()
{
    // Collect edges to flip first, then apply changes, to avoid invalidating iterators.
    struct FlipItem {
        edge_descriptor e;
        vertex_descriptor v0;
        vertex_descriptor v1;
    };
    vector<FlipItem> toFlip;
    BGL_FORALL_EDGES(e, *this, LocalAnchorGraph) {
        auto& edge = (*this)[e];
        if(edge.isHidden()) {
            continue;
        }
        if(edge.info.common == 0) {
            continue;
        }
        if(edge.info.offsetInBases < 0) {
            toFlip.push_back({e, source(e, *this), target(e, *this)});
        }
    }

    for(const auto& item: toFlip) {
        auto& edge = (*this)[item.e];
        if(edge.isHidden()) {
            continue;
        }

        // Hide the negative-offset edge.
        edge.displayFlags |= LocalAnchorGraphEdge::Hidden;
        edge.displayFlags |= LocalAnchorGraphEdge::FlippedByOffset;

        // Create or update the reversed edge.
        LocalAnchorGraphEdge reversed = edge;
        reversed.displayFlags = LocalAnchorGraphEdge::FlippedByOffset;
        reversed.info.reverse();
        reversed.displayFlags &= ~LocalAnchorGraphEdge::Hidden;

        const auto existing = boost::edge(item.v1, item.v0, *this);
        if(existing.second) {
            auto& existingEdge = (*this)[existing.first];
            // Prefer keeping a valid-offset edge, and prefer higher coverage.
            if(existingEdge.info.common == 0 and reversed.info.common > 0) {
                existingEdge.info = reversed.info;
            }
            if(existingEdge.coverage < reversed.coverage) {
                existingEdge.coverage = reversed.coverage;
            }
            existingEdge.displayFlags &= ~LocalAnchorGraphEdge::Hidden;
            existingEdge.displayFlags |= LocalAnchorGraphEdge::FlippedByOffset;
        } else {
            add_edge(item.v1, item.v0, reversed, *this);
        }
    }
}



void LocalAnchorGraph::reduceTransitiveEdges(const uint64_t fuzzBases, const double coverageFactor)
{
    if(num_edges(*this) == 0) {
        return;
    }

    // Build active adjacency and a lookup for direct edges (u,v) -> edge_descriptor.
    std::map<std::pair<vertex_descriptor, vertex_descriptor>, edge_descriptor> edgeMap;
    std::map<vertex_descriptor, vector<edge_descriptor> > outEdges;
    std::map<vertex_descriptor, uint64_t> outCount;
    std::map<vertex_descriptor, uint64_t> inCount;

    BGL_FORALL_EDGES(e, *this, LocalAnchorGraph) {
        auto& edge = (*this)[e];
        if(edge.isHidden()) {
            continue;
        }
        const vertex_descriptor u = source(e, *this);
        const vertex_descriptor v = target(e, *this);
        const auto key = std::make_pair(u, v);
        const auto it = edgeMap.find(key);
        if(it == edgeMap.end()) {
            edgeMap.insert({key, e});
        } else {
            // Keep only one active edge per (u,v): hide the weaker duplicate.
            const edge_descriptor existing = it->second;
            if((*this)[existing].coverage >= edge.coverage) {
                edge.displayFlags |= LocalAnchorGraphEdge::Hidden;
                continue;
            } else {
                (*this)[existing].displayFlags |= LocalAnchorGraphEdge::Hidden;
                it->second = e;
            }
        }
    }

    // Build counts and outgoing lists after duplicate collapsing.
    BGL_FORALL_EDGES(e, *this, LocalAnchorGraph) {
        const auto& edge = (*this)[e];
        if(edge.isHidden()) {
            continue;
        }
        const vertex_descriptor u = source(e, *this);
        const vertex_descriptor v = target(e, *this);
        outEdges[u].push_back(e);
        ++outCount[u];
        ++inCount[v];
    }

    vector<edge_descriptor> candidates;
    candidates.reserve(edgeMap.size() / 4);

    for(const auto& p: outEdges) {
        const vertex_descriptor u = p.first;
        const auto& outU = p.second;
        for(const edge_descriptor e_uv: outU) {
            const auto& edge_uv = (*this)[e_uv];
            if(edge_uv.info.common == 0) {
                continue;
            }
            const int64_t off_uv = edge_uv.info.offsetInBases;
            if(off_uv < 0) {
                continue;
            }
            const vertex_descriptor v = target(e_uv, *this);
            const auto itOutV = outEdges.find(v);
            if(itOutV == outEdges.end()) {
                continue;
            }
            for(const edge_descriptor e_vw: itOutV->second) {
                const auto& edge_vw = (*this)[e_vw];
                if(edge_vw.info.common == 0) {
                    continue;
                }
                const int64_t off_vw = edge_vw.info.offsetInBases;
                if(off_vw < 0) {
                    continue;
                }
                const vertex_descriptor w = target(e_vw, *this);
                const auto itUw = edgeMap.find(std::make_pair(u, w));
                if(itUw == edgeMap.end()) {
                    continue;
                }
                const edge_descriptor e_uw = itUw->second;
                if(e_uw == e_uv or e_uw == e_vw) {
                    continue;
                }
                const auto& edge_uw = (*this)[e_uw];
                if(edge_uw.info.common == 0) {
                    continue;
                }
                const int64_t off_uw = edge_uw.info.offsetInBases;
                if(off_uw < 0) {
                    continue;
                }
                const int64_t implied = off_uv + off_vw;
                const int64_t diff = (off_uw > implied) ? (off_uw - implied) : (implied - off_uw);
                if(uint64_t(diff) > fuzzBases) {
                    continue;
                }

                // Support guard: don't delete a strong direct edge due to a weak 2-hop path.
                const double minPathCoverage = double(min(edge_uv.coverage, edge_vw.coverage));
                const double directCoverage = double(edge_uw.coverage);
                if(minPathCoverage + 1e-9 < coverageFactor * directCoverage) {
                    continue;
                }

                candidates.push_back(e_uw);
            }
        }
    }

    if(candidates.empty()) {
        return;
    }

    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

    // Prefer removing weaker edges first to preserve high-coverage connectivity.
    std::sort(candidates.begin(), candidates.end(),
        [&](const edge_descriptor a, const edge_descriptor b) {
            return (*this)[a].coverage < (*this)[b].coverage;
        });

    for(const edge_descriptor e: candidates) {
        auto& edge = (*this)[e];
        if(edge.isHidden()) {
            continue;
        }
        const vertex_descriptor u = source(e, *this);
        const vertex_descriptor v = target(e, *this);
        if(outCount[u] <= 1 or inCount[v] <= 1) {
            continue;
        }
        edge.displayFlags |= LocalAnchorGraphEdge::Hidden;
        edge.displayFlags |= LocalAnchorGraphEdge::Transitive;
        --outCount[u];
        --inCount[v];
    }
}



void LocalAnchorGraph::highlightBubbles(const uint64_t maxDepth)
{
    const LocalAnchorGraph& graph = *this;
    if(num_edges(graph) == 0) {
        return;
    }

    // Build active edge lookup and outgoing adjacency.
    std::map<std::pair<vertex_descriptor, vertex_descriptor>, edge_descriptor> edgeMap;
    std::map<vertex_descriptor, vector<vertex_descriptor> > out;
    std::map<vertex_descriptor, uint64_t> inCount;
    BGL_FORALL_EDGES(e, graph, LocalAnchorGraph) {
        if(graph[e].isHidden()) {
            continue;
        }
        const vertex_descriptor u = source(e, graph);
        const vertex_descriptor v = target(e, graph);
        edgeMap.insert({{u, v}, e});
        out[u].push_back(v);
        ++inCount[v];
    }

    // Simple bubble heuristic: for each branching vertex s, look for the closest
    // common descendant t reachable from all outgoing neighbors.
    BGL_FORALL_VERTICES(s, graph, LocalAnchorGraph) {
        const auto itOutS = out.find(s);
        if(itOutS == out.end()) {
            continue;
        }
        auto neighbors = itOutS->second;
        if(neighbors.size() < 2) {
            continue;
        }
        if(neighbors.size() > 4) {
            // Avoid overly-branching regions for this lightweight visualization heuristic.
            continue;
        }

        vector< std::map<vertex_descriptor, uint64_t> > dists;
        vector< std::map<vertex_descriptor, vertex_descriptor> > parents;
        vector< std::set<vertex_descriptor> > reach;

        dists.resize(neighbors.size());
        parents.resize(neighbors.size());
        reach.resize(neighbors.size());

        for(size_t i=0; i<neighbors.size(); i++) {
            const vertex_descriptor root = neighbors[i];
            std::queue<vertex_descriptor> q;
            dists[i].insert({root, 0});
            q.push(root);

            while(not q.empty()) {
                const vertex_descriptor x = q.front();
                q.pop();
                const uint64_t dx = dists[i][x];
                if(dx >= maxDepth) {
                    continue;
                }
                const auto itOutX = out.find(x);
                if(itOutX == out.end()) {
                    continue;
                }
                for(const vertex_descriptor y: itOutX->second) {
                    if(dists[i].contains(y)) {
                        continue;
                    }
                    dists[i].insert({y, dx + 1});
                    parents[i].insert({y, x});
                    q.push(y);
                }
            }
            for(const auto& p: dists[i]) {
                reach[i].insert(p.first);
            }
        }

        // Intersection of reachable sets.
        std::set<vertex_descriptor> intersection = reach[0];
        for(size_t i=1; i<reach.size(); i++) {
            std::set<vertex_descriptor> tmp;
            std::set_intersection(
                intersection.begin(), intersection.end(),
                reach[i].begin(), reach[i].end(),
                std::inserter(tmp, tmp.begin()));
            intersection.swap(tmp);
            if(intersection.empty()) {
                break;
            }
        }
        if(intersection.empty()) {
            continue;
        }

        vertex_descriptor bestT{};
        bool foundBestT = false;
        uint64_t bestScore = std::numeric_limits<uint64_t>::max();
        for(const vertex_descriptor t: intersection) {
            if(t == s) {
                continue;
            }
            if(inCount[t] < 2) {
                continue;
            }
            uint64_t maxD = 0;
            bool ok = true;
            for(size_t i=0; i<dists.size(); i++) {
                const auto it = dists[i].find(t);
                if(it == dists[i].end()) {
                    ok = false;
                    break;
                }
                maxD = max(maxD, it->second);
            }
            if(not ok) {
                continue;
            }
            if(maxD < bestScore) {
                bestScore = maxD;
                bestT = t;
                foundBestT = true;
            }
        }
        if(not foundBestT) {
            continue;
        }

        // Mark edges on the shortest paths from each neighbor to bestT, plus the branching edges.
        for(size_t i=0; i<neighbors.size(); i++) {
            const vertex_descriptor n = neighbors[i];
            const auto itSn = edgeMap.find({s, n});
            if(itSn != edgeMap.end()) {
                (*this)[itSn->second].displayFlags |= LocalAnchorGraphEdge::Bubble;
            }

            vertex_descriptor cur = bestT;
            while(cur != n) {
                const auto itP = parents[i].find(cur);
                if(itP == parents[i].end()) {
                    break;
                }
                const vertex_descriptor prev = itP->second;
                const auto itE = edgeMap.find({prev, cur});
                if(itE != edgeMap.end()) {
                    (*this)[itE->second].displayFlags |= LocalAnchorGraphEdge::Bubble;
                }
                cur = prev;
            }
        }
    }
}



// This is the code that uses svg output from graphviz.
void LocalAnchorGraph::writeHtml1(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& options) const
{


        // Write it out in graphviz format.
        const string uuid = to_string(boost::uuids::random_generator()());
        const string dotFileName = tmpDirectory() + uuid + ".dot";
        writeGraphviz(dotFileName, options);

        // Use graphviz to compute the layout.
        const string svgFileName = dotFileName + ".svg";
        const string shape = options.vertexLabels ? "rectangle" : "point";
        string command =
            options.layoutMethod +
            " -T svg " + dotFileName + " -o " + svgFileName +
            " -Nshape=" + shape +
            " -Gsize=" + to_string(options.sizePixels/72) + " -Gratio=expand ";
        if(options.vertexLabels) {
            command += " -Goverlap=false";
        }
        // cout << "Running command: " << command << endl;
        const int timeout = std::max(1, int(std::ceil(options.layoutTimeoutSeconds)));
        bool timeoutTriggered = false;
        bool signalOccurred = false;
        int returnCode = 0;
        runCommandWithTimeout(command, timeout, timeoutTriggered, signalOccurred, returnCode);
        if(signalOccurred) {
            html << "Error during graph layout. Command was<br>" << endl;
            html << command;
            return;
        }
        if(timeoutTriggered) {
            html << "Timeout during graph layout." << endl;
            return;
        }
        if(returnCode!=0 ) {
            html << "Error during graph layout. Command was<br>" << endl;
            html << command;
            return;
        }
        std::filesystem::remove(dotFileName);



        // Write the svg to html.
        html << "<p><div style='border:solid;display:inline-block'>";
        ifstream svgFile(svgFileName);
        html << svgFile.rdbuf();
        svgFile.close();
        html << "</div>";

        // Remove the .svg file.
        std::filesystem::remove(svgFileName);

        // Add drag and zoom.
        addSvgDragAndZoom(html);
    }



// This is the code that computes the graph layout,
// then creates the svg.
void LocalAnchorGraph::writeHtml2(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& options)
{
    // Use scientific notation because svg does not accept floating points
    // ending with a decimal point.
    html << std::scientific;

    computeLayout(options);
    if(not layoutStatusMessage.empty()) {
        html << "<p style='max-width:800px'><b>Layout:</b> " << layoutStatusMessage << "</p>";
    }
    if(options.processingMode != "raw") {
        html << "<p style='max-width:800px'><b>Processing:</b> mode=" << options.processingMode <<
            ", hiddenEdges=" << hiddenEdgeCount <<
            ", flippedEdges=" << flippedEdgeCount <<
            ", transitiveEdges=" << transitiveEdgeCount <<
            ", bubbleEdges=" << bubbleEdgeCount << "</p>";
    }
    if(layout.empty()) {
        html << "<p>Unable to display graph: no layout coordinates available.</p>";
        return;
    }
    computeLayoutBoundingBox();

    Box viewportBox = boundingBox;
    viewportBox.extend(0.05);
    viewportBox.makeSquare();

    // Begin the svg.
    const string svgId = "LocalAnchorGraph";
    html <<
        "\n<br><div style='display:inline-block;vertical-align:top;'>"
        "<svg id='" << svgId <<
        "' width='" <<  options.sizePixels <<
        "' height='" << options.sizePixels <<
        "' viewBox='" << viewportBox.xMin << " " << viewportBox.yMin << " " <<
        viewportBox.xSize() << " " <<
        viewportBox.ySize() << "'"
        " style='background-color:#f0f0f0'"
        ">\n";

    // Write the edges first so they don't obscure the vertices.
    writeEdges(html, options);

    // Write the vertices.
    writeVertices(html, options);

    // Finish the svg.
    html << "</svg></div>";

    // Side panel.
    html << "<div style='display:inline-block;margin-left:20px'>";
    writeSvgControls(html, options);
    html << "</div>";
}



void LocalAnchorGraph::computeLayout(const LocalAnchorGraphDisplayOptions& options)
{
    const LocalAnchorGraph& graph = *this;
    layoutStatusMessage.clear();


    // Compute layout on an induced subgraph containing only visible (non-hidden) edges.
    // This keeps hidden/no-offset edges from influencing the layout.
    using LayoutGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS>;
    using LayoutVertex = LayoutGraph::vertex_descriptor;
    using LayoutEdge = LayoutGraph::edge_descriptor;

    vector<vertex_descriptor> originalVertices;
    originalVertices.reserve(num_vertices(graph));
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexCount = 0;
    BGL_FORALL_VERTICES(v, graph, LocalAnchorGraph) {
        originalVertices.push_back(v);
        vertexIndexMap.insert({v, vertexCount++});
    }

    LayoutGraph layoutGraph(vertexCount);
    std::map<LayoutEdge, double> edgeLengthMap;
    constexpr int64_t maxOffsetInBasesForLayout = 2'000'000; // Keep layout coordinates sane.
    BGL_FORALL_EDGES(e, graph, LocalAnchorGraph) {
        if(graph[e].isHidden()) {
            continue;
        }
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const uint64_t i0 = vertexIndexMap[v0];
        const uint64_t i1 = vertexIndexMap[v1];

        int64_t offsetInBases = graph[e].info.offsetInBases;
        if((graph[e].info.common == 0) or (offsetInBases < 0)) {
            offsetInBases = 10;
        } else if(offsetInBases > maxOffsetInBasesForLayout) {
            offsetInBases = maxOffsetInBasesForLayout;
        }
        const double displayLength =
            options.minimumEdgeLength +
            options.additionalEdgeLengthPerKb * 0.001 * double(offsetInBases);

        LayoutEdge le{};
        bool inserted = false;
        tie(le, inserted) = add_edge(LayoutVertex(i0), LayoutVertex(i1), layoutGraph);
        if(inserted) {
            edgeLengthMap.insert({le, displayLength});
        }
    }

    // Compute the graph layout.
    layout.clear();
    const double timeout = std::max(0.1, options.layoutTimeoutSeconds);
    ComputeLayoutReturnCode rc = ComputeLayoutReturnCode::Error;
    std::map<LayoutVertex, array<double, 2> > layoutPositions;
    if(options.layoutMethod == "custom") {
        const int quality = 2;
        rc = computeLayoutCustom(
            layoutGraph,
            edgeLengthMap,
            layoutPositions,
            quality,
            timeout);
    } else {
        const string additionalOptions = "";
        rc = computeLayoutGraphviz(
            layoutGraph,
            options.layoutMethod,
            timeout,
            layoutPositions,
            additionalOptions,
            &edgeLengthMap);
    }

    // If the layout is dot, reverse the y coordinates so the arrows point down.
    if((rc == ComputeLayoutReturnCode::Success) and (options.layoutMethod == "dot")) {
        for(auto& p: layoutPositions) {
            auto& y = p.second[1];
            y = -y;
        }
    }

    // Map back to the original graph vertex descriptors.
    if(rc == ComputeLayoutReturnCode::Success) {
        for(const auto& p: layoutPositions) {
            const uint64_t i = uint64_t(p.first);
            if(i >= originalVertices.size()) {
                continue;
            }
            layout.insert({originalVertices[i], p.second});
        }
    }

    // Fallback layout if external layout fails or returns nothing.
    if((rc != ComputeLayoutReturnCode::Success) or layout.empty()) {
        layoutStatusMessage =
            "External layout (" + options.layoutMethod + ") failed (" +
            computeLayoutReturnCodeName(rc) + "); using fallback layout.";

        // Deterministic BFS-radial layout by stored vertex distance.
        layout.clear();
        std::map<uint64_t, vector<vertex_descriptor> > rings;
        BGL_FORALL_VERTICES(v, graph, LocalAnchorGraph) {
            rings[graph[v].distance].push_back(v);
        }
        for(auto& ring: rings) {
            auto& vertices = ring.second;
            std::sort(vertices.begin(), vertices.end(),
                [&graph](const vertex_descriptor a, const vertex_descriptor b) {
                    return graph[a].anchorId < graph[b].anchorId;
                });
        }

        const double radiusStep = 100.;
        const double twoPi = 2. * std::acos(-1.);
        for(const auto& ring: rings) {
            const uint64_t distance = ring.first;
            const auto& vertices = ring.second;
            const uint64_t k = vertices.size();
            if(k == 0) {
                continue;
            }
            double radius = radiusStep * double(distance);
            if(distance == 0 and k > 1) {
                radius = 0.2 * radiusStep;
            }
            for(uint64_t i=0; i<k; i++) {
                const double angle = (k == 1) ? 0. : twoPi * double(i) / double(k);
                const double x = radius * std::cos(angle);
                const double y = radius * std::sin(angle);
                layout.insert({vertices[i], array<double, 2>{x, y}});
            }
        }
    }

    // Normalize to a render-friendly coordinate scale (SVG renderers can struggle with ~1e16).
    // Keep relative geometry but translate/scale so spans are ~O(1e3).
    if(not layout.empty()) {
        double xMin = std::numeric_limits<double>::max();
        double xMax = -std::numeric_limits<double>::max();
        double yMin = xMin;
        double yMax = xMax;
        for(const auto& p: layout) {
            const double x = p.second[0];
            const double y = p.second[1];
            if(not (std::isfinite(x) and std::isfinite(y))) {
                continue;
            }
            xMin = min(xMin, x);
            xMax = max(xMax, x);
            yMin = min(yMin, y);
            yMax = max(yMax, y);
        }
        const double xSpan = xMax - xMin;
        const double ySpan = yMax - yMin;
        const double span = max(xSpan, ySpan);
        if(std::isfinite(span) and span > 1e6) {
            const double xCenter = 0.5 * (xMin + xMax);
            const double yCenter = 0.5 * (yMin + yMax);
            const double targetSpan = 1000.;
            const double scale = span / targetSpan;
            for(auto& p: layout) {
                p.second[0] = (p.second[0] - xCenter) / scale;
                p.second[1] = (p.second[1] - yCenter) / scale;
            }
            if(layoutStatusMessage.empty()) {
                layoutStatusMessage = "Layout coordinates normalized for SVG rendering.";
            } else {
                layoutStatusMessage += " Layout coordinates normalized for SVG rendering.";
            }
        }
    }
}



void LocalAnchorGraph::computeLayoutBoundingBox()
{

    if(layout.empty()) {
        boundingBox.xMin = -1.;
        boundingBox.xMax =  1.;
        boundingBox.yMin = -1.;
        boundingBox.yMax =  1.;
        return;
    }
    boundingBox.xMin = std::numeric_limits<double>::max();
    boundingBox.xMax = -std::numeric_limits<double>::max();
    boundingBox.yMin = boundingBox.xMin;
    boundingBox.yMax = boundingBox.xMax;

    bool any = false;
    for(const auto& p: layout) {
        const array<double, 2>& xy = p.second;
        const double x = xy[0];
        const double y = xy[1];
        if(not (std::isfinite(x) and std::isfinite(y))) {
            continue;
        }
        any = true;
        boundingBox.xMin = min(boundingBox.xMin, x);
        boundingBox.xMax = max(boundingBox.xMax, x);
        boundingBox.yMin = min(boundingBox.yMin, y);
        boundingBox.yMax = max(boundingBox.yMax, y);
    }

    // If all layout coordinates are non-finite or identical, ensure a valid viewBox.
    if(not any) {
        boundingBox.xMin = -1.;
        boundingBox.xMax =  1.;
        boundingBox.yMin = -1.;
        boundingBox.yMax =  1.;
        return;
    }
    if(boundingBox.xSize() == 0.) {
        boundingBox.xMin -= 1.;
        boundingBox.xMax += 1.;
    }
    if(boundingBox.ySize() == 0.) {
        boundingBox.yMin -= 1.;
        boundingBox.yMax += 1.;
    }

}



void LocalAnchorGraph::Box::makeSquare()
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



void LocalAnchorGraph::Box::extend(double factor)
{
    const double extend = factor * max(xSize(), ySize());
    xMin -= extend;
    xMax += extend;
    yMin -= extend;
    yMax += extend;
}



void LocalAnchorGraph::writeVertices(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& options) const
{
    const LocalAnchorGraph& graph = *this;

    // SVG uses the layout coordinate system (via viewBox). Convert desired pixel sizes
    // into layout units so vertices/edges remain visible even when the viewBox spans
    // a large coordinate range.
    const double layoutSpan = max(boundingBox.xSize(), boundingBox.ySize());
    const double unitsPerPixel = (layoutSpan > 0.) ?
        (layoutSpan / double(options.sizePixels)) : 1.;
    const double minRadius = 1.5 * unitsPerPixel; // ~1.5px minimum radius for visibility.
    const double scalingFactor = 0.002; // Legacy scaling for coverage-based sizing.

    // Get the reference anchor, if needed.
    AnchorId referenceAnchorId = invalid<AnchorId>;
    if(options.vertexColoring == "byReadComposition") {
        referenceAnchorId = anchorIdFromString(options.referenceAnchorIdString);
        if((referenceAnchorId == invalid<AnchorId>) or (referenceAnchorId >= anchors.size())) {
            throw runtime_error("Invalid reference anchor id " + options.referenceAnchorIdString +
                ". Must be a number between 0 and " +
                to_string(anchors.size() / 2 - 1) + " followed by + or -.");
        }
    }
    const uint64_t referenceAnchorIdCoverage = anchors[referenceAnchorId].coverage();

    html << "\n<g id='vertices' style='stroke:none'>";

    BGL_FORALL_VERTICES(v, graph, LocalAnchorGraph) {
        const LocalAnchorGraphVertex& vertex = graph[v];
        const AnchorId anchorId = vertex.anchorId;
        const string anchorIdString = anchorIdToString(anchorId);
        const uint64_t coverage = anchors[anchorId].coverage();

        // Get the position of this vertex in the computed layout.
        const auto it = layout.find(v);
        DINARA_ASSERT(it != layout.end());
        const auto& p = it->second;
        const double x = p[0];
        const double y = p[1];

        AnchorPairInfo info;
        if(options.vertexColoring == "byReadComposition") {
            anchors.analyzeAnchorPair(referenceAnchorId, anchorId, info);
        }

        // Choose the color for this vertex.
        string color;
        if(vertex.distance == maxDistance) {
            color = "Cyan";
        } else if(vertex.distance == 0) {
            color = "Blue";
        } else {

            // Color by similarity of read composition with the reference Anchor.
            if(options.vertexColoring == "byReadComposition") {

                double hue = 1.;    // 0=red, 1=green.
                if(options.similarityMeasure == "commonCount") {
                    // By common count.
                    hue = double(info.common) / double(referenceAnchorIdCoverage);

                } else if(options.similarityMeasure == "jaccard") {
                    // By Jaccard similarity.
                    hue = info.jaccard();
                } else {
                    // By corrected Jaccard similarity.
                    hue = info.correctedJaccard();
                }
                color = "hsl(" + to_string(uint32_t(std::round(hue * 120.))) +
                    ",100%,50%)";

            } else {

                color = "Black";
            }
        }

        // Hyperlink.
        html << "\n<a href='exploreAnchor?anchorIdString=" <<
            HttpServer::urlEncode(anchorIdString) << "'>";

        const double coverageFactor = options.vertexSizeByCoverage ? double(coverage) : 1.;
        const double radius = max(
            options.vertexSize * (scalingFactor * coverageFactor),
            minRadius);

        // Write the vertex.
        html << "<circle cx='" << x << "' cy='" << y <<
            "' fill='" << color <<
            "' r='" << radius <<
            "' id='" << anchorIdString << "'>"
            "<title>" << anchorIdString << ", coverage " << coverage;
        if(options.vertexColoring == "byReadComposition") {
            html << ", common " << info.common << ", J " <<
                std::fixed << std::setprecision(2) << info.jaccard() <<
                ", J' " << info.correctedJaccard();
            if(info.common > 0) {
                html << ", offset " << info.offsetInBases;
            }
        }
        html << "</title></circle>";

        // End the hyperlink.
        html << "</a>";
    }
    html << "\n</g>";
}




void LocalAnchorGraph::writeEdges(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& options) const
{
    const LocalAnchorGraph& graph = *this;

    const double layoutSpan = max(boundingBox.xSize(), boundingBox.ySize());
    const double unitsPerPixel = (layoutSpan > 0.) ?
        (layoutSpan / double(options.sizePixels)) : 1.;
    const double minStrokeWidth = 0.6 * unitsPerPixel; // ~0.6px minimum for visibility.
    const double scalingFactor = 0.001; // Legacy scaling for coverage-based sizing.

    html << "\n<g id=edges>";

    BGL_FORALL_EDGES(e, graph, LocalAnchorGraph) {
        const LocalAnchorGraphEdge& edge = graph[e];
        if(edge.isHidden()) {
            continue;
        }
        const uint64_t coverage = edge.coverage;
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        // Get the position of these vertices in the computed layout.
        const auto it0 = layout.find(v0);
        DINARA_ASSERT(it0 != layout.end());
        const auto& p0 = it0->second;
        const double x0 = p0[0];
        const double y0 = p0[1];
        const auto it1 = layout.find(v1);
        DINARA_ASSERT(it1 != layout.end());
        const auto& p1 = it1->second;
        const double x1 = p1[0];
        const double y1 = p1[1];

        const LocalAnchorGraphVertex& vertex0 = graph[v0];
        const AnchorId anchorId0 = vertex0.anchorId;
        const string anchorIdString0 = anchorIdToString(anchorId0);
        const LocalAnchorGraphVertex& vertex1 = graph[v1];
        const AnchorId anchorId1 = vertex1.anchorId;
        const string anchorIdString1 = anchorIdToString(anchorId1);

        string color = "Black";

        if(options.edgeColoring == "random") {
            // To decide the color, hash the AnchorIds.
            // This way we always get the same color for the same edge.
            const auto p = make_pair(anchorId0, anchorId1);
            const uint32_t hashValue = MurmurHash2(&p, sizeof(p), 759);
            const uint32_t hue = hashValue % 360;
            color = "hsl(" + to_string(hue) + ",50%,50%)";
        }
        if((options.processingMode == "bubble") and
            ((edge.displayFlags & LocalAnchorGraphEdge::Bubble) != 0)) {
            color = "Magenta";
        }

        // Hyperlink.
        html << "\n<a href='exploreAnchorPair?"
            "anchorIdAString=" << HttpServer::urlEncode(anchorIdString0) << "&"
            "anchorIdBString=" << HttpServer::urlEncode(anchorIdString1) << "'>";

        const double strokeWidth = max(
            scalingFactor * options.edgeThickness * double(coverage),
            minStrokeWidth);

        html <<
            "\n<line x1='" << x0 << "' y1='" << y0 <<
            "' x2='" << x1 << "' y2='" << y1 <<
            "' stroke='" << color <<
            "' stroke-width='" << strokeWidth <<
            "'>"
            "<title>" <<
            anchorIdString0 << " to " << anchorIdString1 <<
            ", coverage " << coverage << "/" << edge.info.common <<
            ", loss ";
        const auto oldPrecision = html.precision(2);
        const auto oldFlags = html.setf(std::ios_base::fixed, std::ios_base::floatfield);
        html << edge.coverageLoss();
        html.precision(oldPrecision);
        html.flags(oldFlags);
        html << "</title>""</line>";

        // End the hyperlink.
        html << "</a>";
    }
    html << "</g>";



    // Write the "arrows" to show edge directions.
    // They are just short lines near the target vertex of each edge.
    html << "\n<g id=arrows";
    if(options.edgeColoring == "black") {
        html << " stroke=white";
    } else {
        html << " stroke=black";
    }
    html << ">";
    BGL_FORALL_EDGES(e, graph, LocalAnchorGraph) {
        if(graph[e].isHidden()) {
            continue;
        }
        const uint64_t coverage = graph[e].coverage;
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        // Get the position of these vertices in the computed layout.
        const auto it0 = layout.find(v0);
        DINARA_ASSERT(it0 != layout.end());
        const auto& p0 = it0->second;
        const double x0 = p0[0];
        const double y0 = p0[1];
        const auto it1 = layout.find(v1);
        DINARA_ASSERT(it1 != layout.end());
        const auto& p1 = it1->second;
        const double x1 = p1[0];
        const double y1 = p1[1];

        const double relativeArrowLength = 0.3;
        const double x2 = (1. - relativeArrowLength) * x1 + relativeArrowLength * x0;
        const double y2 = (1. - relativeArrowLength) * y1 + relativeArrowLength * y0;

        const double arrowStrokeWidth = max(
            0.2 * scalingFactor * options.edgeThickness * double(coverage),
            0.3 * unitsPerPixel);
        html <<
            "\n<line x1='" << x1 << "' y1='" << y1 <<
            "' x2='" << x2 << "' y2='" << y2 <<
            "' stroke-width='" << arrowStrokeWidth <<
            "' />";

    }
    html << "</g>";
}



void LocalAnchorGraph::writeSvgControls(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& /* options */) const
{
    html <<
        "<p><table>";

    // Add drag and zoom.
    addSvgDragAndZoom(html);



    // Buttons to change vertex size.
    html << R"stringDelimiter(
    <tr><th class=left>Vertex size<td>
    <button type='button' onClick='changeVertexSize(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeVertexSize(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeVertexSize(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeVertexSize(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeVertexSize(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeVertexSize(10.)' style='width:3em'>+++</button>
        <script>
        function changeVertexSize(factor)
        {
            var vertexGroup = document.getElementById('vertices');
            var vertices = vertexGroup.getElementsByTagName('circle');
            for(i=0; i<vertices.length; i++) {
                v = vertices[i];
                v.setAttribute('r', factor * v.getAttribute('r'));
            }
        }
        </script>
        )stringDelimiter";



    // Buttons to change edge thickness.
    html << R"stringDelimiter(
    <tr><th class=left>Edge thickness<td>
    <button type='button' onClick='changeThickness(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeThickness(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeThickness(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeThickness(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeThickness(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeThickness(10.)' style='width:3em'>+++</button>
        <script>
        function changeThickness(factor)
        {
            var edgeGroup = document.getElementById('edges');
            var edges = edgeGroup.getElementsByTagName('line');
            for(i=0; i<edges.length; i++) {
                e = edges[i];
                e.setAttribute('stroke-width', factor * e.getAttribute('stroke-width'));
            }

            var arrowsGroup = document.getElementById('arrows');
            var arrows = arrowsGroup.getElementsByTagName('line');
            for(i=0; i<arrows.length; i++) {
                a = arrows[i];
                a.setAttribute('stroke-width', factor * a.getAttribute('stroke-width'));
            }
        }
        </script>
        )stringDelimiter";



    // Zoom buttons.
    html << R"stringDelimiter(
        <tr title='Or use the mouse wheel.'><th class=left>Zoom<td>
        <button type='button' onClick='zoomSvg(0.1)' style='width:3em'>---</button>
        <button type='button' onClick='zoomSvg(0.5)' style='width:3em'>--</button>
        <button type='button' onClick='zoomSvg(0.8)' style='width:3em'>-</button>
        <button type='button' onClick='zoomSvg(1.25)' style='width:3em'>+</button>
        <button type='button' onClick='zoomSvg(2.)' style='width:3em'>++</button>
        <button type='button' onClick='zoomSvg(10.)' style='width:3em'>+++</button>
    )stringDelimiter";



    // Buttons to highlight an anchor and zoom to an anchor.
    html << R"stringDelimiter(
        <tr><td colspan=2>
        <button onClick='highlightAnchor()'>Highlight</button>
        <button onClick='zoomToAnchor()'>Zoom to</button>anchor
        <input id=selectedAnchorId type=text size=10 style='text-align:center'>
    <script>
    function zoomToAnchor()
    {
        // Get the anchor id from the input field.
        var anchorId = document.getElementById("selectedAnchorId").value;
        zoomToGivenAnchor(anchorId);
    }
    function zoomToGivenAnchor(anchorId)
    {
        var element = document.getElementById(anchorId);
        // Find the bounding box and its center.
        var box = element.getBBox();
        var xCenter = box.x + 0.5 * box.width;
        var yCenter = box.y + 0.5 * box.height;

        // Change the viewbox of the svg to be a bit larger than a square
        // containing the bounding box.
        var enlargeFactor = 5.;
        var size = enlargeFactor * Math.max(box.width, box.height);
        var factor = size / width;
        width = size;
        height = size;
        x = xCenter - 0.5 * size;
        y = yCenter - 0.5 * size;
        var svg = document.querySelector('svg');
        svg.setAttribute('viewBox', `${x} ${y} ${size} ${size}`);
        ratio = size / svg.getBoundingClientRect().width;
        svg.setAttribute('font-size', svg.getAttribute('font-size') * factor);
    }
    function highlightAnchor()
    {
        // Get the anchor id  from the input field.
        var anchorId = document.getElementById("selectedAnchorId").value;
        var element = document.getElementById(anchorId);

        element.style.fill = "Magenta";
    }
    </script>
    )stringDelimiter";


    html << "</table>";

    // Scroll down to the svg.
    const string svgId = "LocalAnchorGraph";
    html <<
        "<script>"
        "document.getElementById('" << svgId << "').scrollIntoView({block:'center'});"
        "</script>";

    html <<
        "<p>Use Ctrl+Click to pan."
        "<p>Use Ctrl-Wheel or the above buttons to zoom.";
}
