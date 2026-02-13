// HTTP server handlers for the DirectedAnchorGraph (Verkko-style resolution).
//
// Two-node representation: each segment has oriented DagNodeIds
// (even=fwd, odd=rev).  Display uses segmentOf() and isForward().

#include "Assembler.hpp"
#include "mode3-DirectedAnchorGraph.hpp"
#include "HttpServer.hpp"
#include "invalid.hpp"

using namespace dinara;
using namespace mode3;
using namespace std;



// ============================================================================
// Graph summary page.
// ============================================================================

void Assembler::exploreDirectedAnchorGraph(const vector<string>& request, ostream& html)
{
    (void)request;

    if(!directedAnchorGraph) {
        html << "<p>The Directed Anchor Graph is not available. "
            "Run <code>runDirectedAnchorGraphResolution</code> first.";
        return;
    }

    const auto& dag = *directedAnchorGraph;

    html << "<h1>Directed Anchor Graph Summary</h1>";

    html <<
        "<table>"
        "<tr><th class=left>Active nodes<td class=centered>" << dag.nodeCount() <<
        "<tr><th class=left>Total nodes (incl. removed)<td class=centered>" << dag.totalNodeCount() <<
        "<tr><th class=left>Directed edges<td class=centered>" << dag.edgeCount() <<
        "<tr><th class=left>Active paths<td class=centered>" << dag.pathCount() <<
        "</table>";


    // Degree distribution.
    const auto activeIds = dag.getActiveNodeIds();

    map<uint64_t, uint64_t> outDegreeDist;
    map<uint64_t, uint64_t> inDegreeDist;
    for(const uint64_t segId : activeIds) {
        DagNodeId fwd = fwdNodeId(segId);
        DagNodeId rev = revNodeId(segId);
        outDegreeDist[dag.outDegree(fwd)]++;
        outDegreeDist[dag.outDegree(rev)]++;
        inDegreeDist[dag.inDegree(fwd)]++;
        inDegreeDist[dag.inDegree(rev)]++;
    }

    html << "<h2>Out-degree distribution</h2>"
        "<table><tr><th>Degree<th>Count";
    for(const auto& [deg, cnt] : outDegreeDist) {
        html << "<tr><td class=centered>" << deg << "<td class=centered>" << cnt;
    }
    html << "</table>";

    html << "<h2>In-degree distribution</h2>"
        "<table><tr><th>Degree<th>Count";
    for(const auto& [deg, cnt] : inDegreeDist) {
        html << "<tr><td class=centered>" << deg << "<td class=centered>" << cnt;
    }
    html << "</table>";


    // Chain length distribution.
    map<uint64_t, uint64_t> chainLengthDist;
    for(const uint64_t segId : activeIds) {
        chainLengthDist[dag.getNode(segId).anchorChain.size()]++;
    }

    html << "<h2>Anchor chain length distribution</h2>"
        "<p>Number of base anchors merged into each node."
        "<table><tr><th>Chain length<th>Count";
    for(const auto& [len, cnt] : chainLengthDist) {
        html << "<tr><td class=centered>" << len << "<td class=centered>" << cnt;
    }
    html << "</table>";


    // Coverage distribution.
    map<uint64_t, uint64_t> coverageDist;
    for(const uint64_t segId : activeIds) {
        uint64_t covBucket = uint64_t(dag.getNode(segId).coverage);
        coverageDist[covBucket]++;
    }

    html << "<h2>Coverage distribution</h2>"
        "<table><tr><th>Coverage<th>Count";
    for(const auto& [cov, cnt] : coverageDist) {
        html << "<tr><td class=centered>" << cov << "<td class=centered>" << cnt;
    }
    html << "</table>";


    // List first 200 active nodes.
    const uint64_t maxList = 200;
    const uint64_t showCount = min(uint64_t(activeIds.size()), maxList);

    html << "<h2>Active nodes";
    if(activeIds.size() > maxList) {
        html << " (first " << maxList << " of " << activeIds.size() << ")";
    }
    html << "</h2>"
        "<table>"
        "<tr><th>Node<th>Chain length<th>Coverage"
        "<th>Out-deg (&gt;)<th>In-deg (&gt;)"
        "<th>Out-deg (&lt;)<th>In-deg (&lt;)<th>Paths crossing";
    for(uint64_t i = 0; i < showCount; i++) {
        const uint64_t segId = activeIds[i];
        const auto& info = dag.getNode(segId);
        DagNodeId fwd = fwdNodeId(segId);
        DagNodeId rev = revNodeId(segId);
        const uint64_t crossingCount = dag.getPathsCrossingNode(segId).size();

        html <<
            "<tr>"
            "<td class=centered><a href='exploreDirectedAnchorGraphNode?nodeId=" << segId << "'>"
            << segId << "</a>"
            "<td class=centered>" << info.anchorChain.size() <<
            "<td class=centered>" << info.coverage <<
            "<td class=centered>" << dag.outDegree(fwd) <<
            "<td class=centered>" << dag.inDegree(fwd) <<
            "<td class=centered>" << dag.outDegree(rev) <<
            "<td class=centered>" << dag.inDegree(rev) <<
            "<td class=centered>" << crossingCount;
    }
    html << "</table>";
}



// ============================================================================
// Node explorer page.
// ============================================================================

void Assembler::exploreDirectedAnchorGraphNode(const vector<string>& request, ostream& html)
{
    if(!directedAnchorGraph) {
        html << "<p>The Directed Anchor Graph is not available.";
        return;
    }

    const auto& dag = *directedAnchorGraph;

    uint64_t nodeIdValue = invalid<uint64_t>;
    const bool nodeIdIsPresent =
        HttpServer::getParameterValue(request, "nodeId", nodeIdValue);

    html <<
        "<h2>Directed Anchor Graph Node</h2>"
        "<p>Each node (segment) has two oriented versions: "
        "&gt;segId (forward) and &lt;segId (reverse complement)."
        "<form>"
        "<table>"
        "<tr><th class=left>Segment id"
        "<td class=centered>"
        "<input type=text name=nodeId required style='text-align:center'";
    if(nodeIdIsPresent) {
        html << " value='" << nodeIdValue << "'";
    }
    html <<
        " size=8 title='Enter a segment id.'>"
        "</table>"
        "<input type=submit value='Show node details'> "
        "</form>";

    if(!nodeIdIsPresent) {
        return;
    }

    const uint64_t segId = nodeIdValue;
    if(!dag.nodeExists(segId)) {
        html << "<p>Segment " << segId << " does not exist or has been removed.";
        return;
    }

    const auto& info = dag.getNode(segId);
    DagNodeId fwd = fwdNodeId(segId);
    DagNodeId rev = revNodeId(segId);

    html << "<h1>Directed Anchor Graph Segment " << segId << "</h1>";

    html <<
        "<table>"
        "<tr><th class=left>Anchor chain length<td class=centered>" << info.anchorChain.size() <<
        "<tr><th class=left>Length (bp)<td class=centered>" << info.lengthBp <<
        "<tr><th class=left>Coverage<td class=centered>" << info.coverage <<
        "<tr><th class=left>Out-degree (&gt;" << segId << ")<td class=centered>" << dag.outDegree(fwd) <<
        "<tr><th class=left>In-degree (&gt;" << segId << ")<td class=centered>" << dag.inDegree(fwd) <<
        "<tr><th class=left>Out-degree (&lt;" << segId << ")<td class=centered>" << dag.outDegree(rev) <<
        "<tr><th class=left>In-degree (&lt;" << segId << ")<td class=centered>" << dag.inDegree(rev) <<
        "</table>";


    // Anchor chain.
    html << "<h2>Anchor chain</h2>"
        "<p>The base anchors composing this segment (forward orientation). "
        "Each anchor ID is an AnchorId; even=forward, odd=reverse."
        "<table>"
        "<tr><th>Position<th>Anchor pair<th>Orientation<th>AnchorId";
    for(uint64_t i = 0; i < info.anchorChain.size(); i++) {
        DagNodeId aid = info.anchorChain[i];
        html <<
            "<tr>"
            "<td class=centered>" << i <<
            "<td class=centered>" << segmentOf(aid) <<
            "<td class=centered>" << (isForward(aid) ? "&gt; (fwd)" : "&lt; (rev)") <<
            "<td class=centered>" << aid;
    }
    html << "</table>";


    // Edge tables — display neighbors as oriented segment references.
    auto writeEdgeTable = [&](const string& title, DagNodeId orientedNode, bool outgoing) {
        // For outgoing: getOutEdges(orientedNode).
        // For incoming: predecessors of orientedNode = { Y^1 : Y ∈ getOutEdges(orientedNode^1) }.
        html << "<h2>" << title << "</h2>";

        if(outgoing) {
            const auto& nbrs = dag.getOutEdges(orientedNode);
            if(nbrs.empty()) {
                html << "<p>None.";
                return;
            }
            html <<
                "<table>"
                "<tr><th>Segment<th>Orientation";
            for(DagNodeId nbr : nbrs) {
                html <<
                    "<tr>"
                    "<td class=centered><a href='exploreDirectedAnchorGraphNode?nodeId="
                    << segmentOf(nbr) << "'>" << segmentOf(nbr) << "</a>"
                    "<td class=centered>" << (isForward(nbr) ? "&gt; (fwd)" : "&lt; (rev)");
            }
            html << "</table>";
        } else {
            // Incoming: predecessors = { Y^1 : Y ∈ getOutEdges(orientedNode^1) }
            const auto& revNbrs = dag.getOutEdges(orientedNode ^ 1);
            if(revNbrs.empty()) {
                html << "<p>None.";
                return;
            }
            html <<
                "<table>"
                "<tr><th>Segment<th>Orientation";
            for(DagNodeId nbr : revNbrs) {
                DagNodeId pred = nbr ^ 1;
                html <<
                    "<tr>"
                    "<td class=centered><a href='exploreDirectedAnchorGraphNode?nodeId="
                    << segmentOf(pred) << "'>" << segmentOf(pred) << "</a>"
                    "<td class=centered>" << (isForward(pred) ? "&gt; (fwd)" : "&lt; (rev)");
            }
            html << "</table>";
        }
    };

    writeEdgeTable("Incoming edges (&gt;" + to_string(segId) + ")", fwd, false);
    writeEdgeTable("Outgoing edges (&gt;" + to_string(segId) + ")", fwd, true);
    writeEdgeTable("Incoming edges (&lt;" + to_string(segId) + ")", rev, false);
    writeEdgeTable("Outgoing edges (&lt;" + to_string(segId) + ")", rev, true);


    // Paths crossing this segment.
    const auto& crossingPathIndices = dag.getPathsCrossingNode(segId);
    html << "<h2>Paths crossing this segment (" << crossingPathIndices.size() << ")</h2>";

    if(crossingPathIndices.empty()) {
        html << "<p>No paths cross this segment.";
    } else {
        vector<uint64_t> sortedIndices(crossingPathIndices.begin(), crossingPathIndices.end());
        sort(sortedIndices.begin(), sortedIndices.end());

        const uint64_t maxShow = 100;
        const uint64_t showCount = min(uint64_t(sortedIndices.size()), maxShow);

        html <<
            "<table>"
            "<tr><th>Path<th>Length<th>Context (5 nodes around this segment)";

        for(uint64_t i = 0; i < showCount; i++) {
            const uint64_t pathIdx = sortedIndices[i];
            if(!dag.pathExists(pathIdx)) continue;
            const auto& path = dag.getPath(pathIdx);

            // Find this segment's position in the path.
            int64_t posInPath = -1;
            for(uint64_t j = 0; j < path.size(); j++) {
                if(segmentOf(path[j]) == segId) {
                    posInPath = int64_t(j);
                    break;
                }
            }

            html <<
                "<tr>"
                "<td class=centered><a href='exploreDirectedAnchorGraphPath?pathId=" << pathIdx << "'>"
                << pathIdx << "</a>"
                "<td class=centered>" << path.size() <<
                "<td class=centered>";

            const int64_t contextRadius = 2;
            const int64_t start = max(int64_t(0), posInPath - contextRadius);
            const int64_t end = min(int64_t(path.size()) - 1, posInPath + contextRadius);

            if(start > 0) html << "... ";
            for(int64_t j = start; j <= end; j++) {
                if(j > start) html << " &rarr; ";
                DagNodeId n = path[j];
                if(j == posInPath) html << "<b>";
                html << (isForward(n) ? "&gt;" : "&lt;") << segmentOf(n);
                if(j == posInPath) html << "</b>";
            }
            if(end < int64_t(path.size()) - 1) html << " ...";
        }
        html << "</table>";

        if(sortedIndices.size() > maxShow) {
            html << "<p>Showing " << maxShow << " of " << sortedIndices.size() << " paths.";
        }
    }
}



// ============================================================================
// Path explorer page.
// ============================================================================

void Assembler::exploreDirectedAnchorGraphPath(const vector<string>& request, ostream& html)
{
    if(!directedAnchorGraph) {
        html << "<p>The Directed Anchor Graph is not available.";
        return;
    }

    const auto& dag = *directedAnchorGraph;

    uint64_t pathIdValue = invalid<uint64_t>;
    const bool pathIdIsPresent =
        HttpServer::getParameterValue(request, "pathId", pathIdValue);

    html <<
        "<h2>Directed Anchor Graph Path</h2>"
        "<p>A path is a read's traversal through the directed anchor graph. "
        "Each step is an oriented segment (>&gt;segId or &lt;segId)."
        "<form>"
        "<table>"
        "<tr><th class=left>Path id"
        "<td class=centered>"
        "<input type=text name=pathId required style='text-align:center'";
    if(pathIdIsPresent) {
        html << " value='" << pathIdValue << "'";
    }
    html <<
        " size=8 title='Enter a path id.'>"
        "</table>"
        "<input type=submit value='Show path details'> "
        "</form>";

    if(!pathIdIsPresent) {
        return;
    }

    const uint64_t pathIdx = pathIdValue;
    if(!dag.pathExists(pathIdx)) {
        html << "<p>Path " << pathIdx << " does not exist or has been removed.";
        return;
    }

    const auto& path = dag.getPath(pathIdx);

    html << "<h1>Directed Anchor Graph Path " << pathIdx << "</h1>";

    html <<
        "<table>"
        "<tr><th class=left>Path length (nodes)<td class=centered>" << path.size() <<
        "</table>";

    if(path.empty()) {
        html << "<p>This path is empty.";
        return;
    }

    // Main path table.
    html <<
        "<table>"
        "<tr>"
        "<th>Position"
        "<th>Segment"
        "<th>Orientation"
        "<th>Chain length"
        "<th>Coverage"
        "<th>Out-degree"
        "<th>In-degree"
        "<th>Edge";

    for(uint64_t pos = 0; pos < path.size(); pos++) {
        DagNodeId n = path[pos];
        uint64_t segId = segmentOf(n);

        uint64_t chainLen = 0;
        double cov = 0.0;
        uint64_t outDeg = 0;
        uint64_t inDeg = 0;
        bool segOk = dag.nodeExists(segId);
        if(segOk) {
            const auto& info = dag.getNode(segId);
            chainLen = info.anchorChain.size();
            cov = info.coverage;
            outDeg = dag.outDegree(n);
            inDeg = dag.inDegree(n);
        }

        html <<
            "<tr>"
            "<td class=centered>" << pos <<
            "<td class=centered>";

        if(segOk) {
            html << "<a href='exploreDirectedAnchorGraphNode?nodeId=" << segId << "'>"
                << segId << "</a>";
        } else {
            html << segId << " (removed)";
        }

        html <<
            "<td class=centered>" << (isForward(n) ? "&gt; (fwd)" : "&lt; (rev)") <<
            "<td class=centered>" << chainLen <<
            "<td class=centered>" << cov <<
            "<td class=centered>" << outDeg <<
            "<td class=centered>" << inDeg;

        // Edge indicator to next step.
        html << "<td class=centered>";
        if(pos + 1 < path.size()) {
            DagNodeId next = path[pos + 1];
            bool sameStrand = (isForward(n) == isForward(next));
            if(sameStrand) {
                html << "same-strand &rarr;";
            } else {
                html << "diff-strand &harr;";
            }
        }
    }
    html << "</table>";


    // Compact path notation (Verkko GFA style).
    html << "<h2>Compact path</h2><p style='font-family:monospace; word-break:break-all'>";
    for(uint64_t pos = 0; pos < path.size(); pos++) {
        if(pos > 0) html << ",";
        DagNodeId n = path[pos];
        html << (isForward(n) ? "&gt;" : "&lt;") << segmentOf(n);
    }
    html << "</p>";
}
