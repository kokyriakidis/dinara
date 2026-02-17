
#include "Assembler.hpp"
#include "Shasta2AnchorGraph.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "HttpServer.hpp"
#include "invalid.hpp"

using namespace dinara;
using namespace std;

namespace {

void accessShasta2HttpData(Assembler& assembler)
{
    const MappedMemoryOwner shasta2Owner = assembler.shasta2MappedMemoryOwner();

    // Anchors are needed by all Shasta2 HTTP views.
    if(!assembler.shasta2Anchors) {
        assembler.shasta2Anchors = make_shared<Shasta2Anchors>(
            "",
            shasta2Owner,
            assembler.getReads(),
            assembler.assemblerInfo->k,
            *assembler.markers);
    }

    // Load the graph on demand if needed.
    if(!assembler.shasta2AnchorGraph) {
        try {
            assembler.shasta2AnchorGraph = make_shared<Shasta2AnchorGraph>(
                shasta2Owner,
                "Shasta2AnchorGraphAfterTransitiveReduction");
        } catch(const exception&) {
            assembler.shasta2AnchorGraph = make_shared<Shasta2AnchorGraph>(
                shasta2Owner,
                "Shasta2AnchorGraph");
        }
    }
}

}

// ============================================================================
// Shasta2 Anchor Graph Summary
// ============================================================================

void Assembler::exploreShasta2AnchorGraph(const vector<string>& request, ostream& html)
{
    (void)request;

    try {
        accessShasta2HttpData(*this);
    } catch(const exception&) {
    }

    if(!shasta2AnchorGraph) {
        html << "<p>The Shasta2 Anchor Graph is not available.";
        return;
    }

    const auto& graph = *shasta2AnchorGraph;

    html << "<h1>Shasta2 Anchor Graph Summary</h1>";

    html <<
        "<table>"
        "<tr><th class=left>Vertices<td class=centered>" << num_vertices(graph) <<
        "<tr><th class=left>Edges<td class=centered>" << num_edges(graph) <<
        "</table>";

    // Add degree distribution later if needed.
}



// ============================================================================
// Shasta2 Anchor Details
// ============================================================================

void Assembler::exploreShasta2Anchor(const vector<string>& request, ostream& html)
{
    try {
        accessShasta2HttpData(*this);
    } catch(const exception&) {
    }

    if(!shasta2Anchors || !shasta2AnchorGraph) {
        html << "<p>Shasta2 Anchors or Graph not available.";
        return;
    }

    uint64_t anchorId = invalid<uint64_t>;
    const bool anchorIdIsPresent =
        HttpServer::getParameterValue(request, "anchorId", anchorId);

    html <<
        "<h2>Shasta2 Anchor</h2>"
        "<form>"
        "<table>"
        "<tr><th class=left>Anchor Id"
        "<td class=centered>"
        "<input type=text name=anchorId required style='text-align:center'";
    if(anchorIdIsPresent) {
        html << " value='" << anchorId << "'";
    }
    html <<
        " size=8 title='Enter an anchor id.'>"
        "</table>"
        "<input type=submit value='Show anchor details'> "
        "</form>";

    if(!anchorIdIsPresent) {
        return;
    }

    if(anchorId >= shasta2Anchors->size()) {
        html << "<p>Invalid Anchor Id " << anchorId;
        return;
    }

    const auto anchor = (*shasta2Anchors)[anchorId];

    html << "<h1>Shasta2 Anchor " << anchorId << "</h1>";

    // Basic info.
    html <<
        "<table>"
        "<tr><th class=left>Coverage (markers)<td class=centered>" << anchor.coverage() <<
        "</table>";
    
    // Neighbors in the graph.
    // The Shasta2AnchorGraph vertices correspond 1:1 to Shasta2Anchors (if constructed that way).
    // The graph uses Boost Graph Library.
    // Vertex descriptor is uint64_t, usually same as AnchorId.
    
    // Check if anchorId exists in graph as a vertex.
    // Assuming num_vertices(graph) == num_anchors.
    
    using Graph = Shasta2AnchorGraphBaseClass;
    Graph& g = *shasta2AnchorGraph;
    
    if(anchorId >= num_vertices(g)) {
        html << "<p>Anchor Id " << anchorId << " is not in the graph (or graph is smaller).";
    } else {
        
        // Out edges.
        html << "<h2>Children (Out-Edges)</h2>";
        auto outEdges = boost::out_edges(anchorId, g);
        if(outEdges.first == outEdges.second) {
             html << "<p>None.";
        } else {
            html << "<table><tr><th>Target Anchor<th>Coverage<th>Analyze";
            for(auto it = outEdges.first; it != outEdges.second; ++it) {
                auto edge = *it;
                auto target = boost::target(edge, g);
                auto& edgeProp = g[edge];
                
                html << "<tr><td class=centered><a href='exploreShasta2Anchor?anchorId=" << target << "'>" << target << "</a>"
                     << "<td class=centered>" << edgeProp.coverage(); // from Shasta2AnchorPair
                html << "<td class=centered><a href='exploreShasta2AnchorPair?anchorId0=" << anchorId << "&anchorId1=" << target << "'>Analyze Pair</a>";

            }
            html << "</table>";
        }

        // In edges.
        html << "<h2>Parents (In-Edges)</h2>";
        auto inEdges = boost::in_edges(anchorId, g);
        if(inEdges.first == inEdges.second) {
             html << "<p>None.";
        } else {
            html << "<table><tr><th>Source Anchor<th>Coverage<th>Analyze";
            for(auto it = inEdges.first; it != inEdges.second; ++it) {
                auto edge = *it;
                auto source = boost::source(edge, g);
                auto& edgeProp = g[edge];
                
                html << "<tr><td class=centered><a href='exploreShasta2Anchor?anchorId=" << source << "'>" << source << "</a>"
                     << "<td class=centered>" << edgeProp.coverage();
                html << "<td class=centered><a href='exploreShasta2AnchorPair?anchorId0=" << source << "&anchorId1=" << anchorId << "'>Analyze Pair</a>";
            }
            html << "</table>";
        }
    }


    // Markers.
    html << "<h2>Markers</h2>";
    html << "<table><tr><th>Read Id<th>Strand<th>Ordinal<th>Position in Journey";
    const uint64_t maxShow = 100;
    uint64_t displayed = 0;
    for(const auto& markerInfo : anchor) {
        if(displayed++ >= maxShow) {
            html << "<tr><td colspan=4 class=centered>... (" << (anchor.size() - maxShow) << " more) ...";
            break;
        }
        html << "<tr><td class=centered><a href='exploreRead?readId=" << markerInfo.orientedReadId.getReadId() << "'>" << markerInfo.orientedReadId.getReadId() << "</a>"
             << "<td class=centered>" << markerInfo.orientedReadId.getStrand()
             << "<td class=centered>" << markerInfo.ordinal
             << "<td class=centered>" << markerInfo.positionInJourney;
    }
    html << "</table>";

}


// ============================================================================
// Shasta2 Anchor Pair Analysis
// ============================================================================

void Assembler::exploreShasta2AnchorPair(const vector<string>& request, ostream& html)
{
    try {
        accessShasta2HttpData(*this);
    } catch(const exception&) {
    }

    if(!shasta2Anchors) {
         html << "<p>Shasta2 Anchors not available.";
         return;
    }

    uint64_t anchorId0 = invalid<uint64_t>;
    uint64_t anchorId1 = invalid<uint64_t>;
    bool present0 = HttpServer::getParameterValue(request, "anchorId0", anchorId0);
    bool present1 = HttpServer::getParameterValue(request, "anchorId1", anchorId1);

    html << "<h2>Analyze Shasta2 Anchor Pair</h2>";
    html << "<form><table>";
    html << "<tr><th class=left>Anchor Id 0<td class=centered><input type=text name=anchorId0 size=8" << (present0 ? (" value='" + to_string(anchorId0) + "'") : "") << ">";
    html << "<tr><th class=left>Anchor Id 1<td class=centered><input type=text name=anchorId1 size=8" << (present1 ? (" value='" + to_string(anchorId1) + "'") : "") << ">";
    html << "</table><input type=submit value='Analyze'></form>";

    if(!present0 || !present1) return;

    if(anchorId0 >= shasta2Anchors->size() || anchorId1 >= shasta2Anchors->size()) {
        html << "<p>Invalid anchor IDs.";
        return;
    }

    html << "<h1>Pair " << anchorId0 << " -> " << anchorId1 << "</h1>";

    Shasta2AnchorPairInfo info;
    shasta2Anchors->analyzeAnchorPair(anchorId0, anchorId1, info);

    html << "<table>";
    html << "<tr><th class=left>Total Reads in Anchor " << anchorId0 << "<td class=centered>" << info.totalA;
    html << "<tr><th class=left>Total Reads in Anchor " << anchorId1 << "<td class=centered>" << info.totalB;
    html << "<tr><th class=left>Common Reads<td class=centered>" << info.common;
    html << "<tr><th class=left>Only in " << anchorId0 << "<td class=centered>" << (info.totalA - info.common);
    html << "<tr><th class=left>Only in " << anchorId1 << "<td class=centered>" << (info.totalB - info.common);
    html << "<tr><th class=left>Jaccard<td class=centered>" << info.jaccard();
    html << "</table>";

    // Could show list of common reads here if needed.
}
