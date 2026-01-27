// Dinara.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "LocalStringGraph.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace dinara;

// Standard library.
#include "chrono.hpp"
#include <unordered_set>



void Assembler::exploreStringGraph(
    const vector<string>& request,
    ostream& html)
{
    try {
        checkStringGraphIsOpen();
    } catch (...) {
        html << "The string graph is not available." << endl;
        return;
    }
    const bool incomingAdjacencyAvailable = stringGraph.incoming.isOpen();

    // Parameters.
    vector<OrientedReadId> readIds;
    string readIdsString;
    const bool readIdsArePresent = getParameterValue(request, "readId", readIdsString);
    const bool readStringsAreValid = parseCommaSeparatedReadIDs(readIdsString, readIds, html);

    string pickActiveString;
    const bool pickActive = getParameterValue(request, "pickActive", pickActiveString);

    string listActiveString;
    const bool listActive = getParameterValue(request, "listActive", listActiveString);
    uint32_t activeLimit = 30;
    getParameterValue(request, "activeLimit", activeLimit);

    uint32_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    string allowChimericReadsString;
    const bool allowChimericReads = getParameterValue(request, "allowChimericReads", allowChimericReadsString);

    const auto pickActiveStringGraphStartVertices = [&](uint32_t limit) -> vector<OrientedReadId> {
        vector<OrientedReadId> result;
        if (!stringGraph.arcs.isOpen) {
            return result;
        }
        result.reserve(limit);
        std::unordered_set<uint32_t> seen;
        seen.reserve(limit * 2);

        const uint64_t arcCount = stringGraph.arcs.size();
        for (uint64_t arcId = 0; arcId < arcCount && result.size() < limit; ++arcId) {
            const auto& arc = stringGraph.arcs[arcId];
            if (arc.del) continue;

            const auto considerVertex = [&](uint32_t v) {
                if (result.size() >= limit) return;
                if (!seen.insert(v).second) return;
                const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
                if (!allowChimericReads && reads->getFlags(orientedReadId.getReadId()).isChimeric) {
                    return;
                }
                if (stringGraph.readDeleted.isOpen &&
                    orientedReadId.getReadId() < stringGraph.readDeleted.size() &&
                    stringGraph.readDeleted[orientedReadId.getReadId()]) {
                    return;
                }
                result.push_back(orientedReadId);
            };

            considerVertex(arc.from);
            considerVertex(arc.to);
        }
        return result;
    };

    string followOutgoingString;
    const bool followOutgoing = getParameterValue(request, "followOutgoing", followOutgoingString);
    string followIncomingString;
    const bool followIncoming = incomingAdjacencyAvailable && getParameterValue(request, "followIncoming", followIncomingString);
    const bool followNoneSpecified = (!followOutgoing && !followIncoming);

    string layoutMethod = "sfdp";
    getParameterValue(request, "layoutMethod", layoutMethod);

    uint32_t sizePixels = 600;
    getParameterValue(request, "sizePixels", sizePixels);

    double vertexScalingFactor = 1.;
    getParameterValue(request, "vertexScalingFactor", vertexScalingFactor);

    double edgeThicknessScalingFactor = 1.;
    getParameterValue(request, "edgeThicknessScalingFactor", edgeThicknessScalingFactor);

    double timeout = 30;
    getParameterValue(request, "timeout", timeout);



    // Form.
    html <<
        "<h3>Display a local subgraph of the string graph</h3>"
        "<p>Tip: following both incoming and outgoing arcs can explode the neighborhood. "
        "For path-like exploration, use outgoing-only.</p>"
        "<p>"
        "<a href='exploreStringGraph?pickActive=on'>Pick an active start vertex</a>"
        "&nbsp;&nbsp;"
        "<a href='exploreStringGraph?listActive=on'>List some active vertices</a>"
        "</p>"
        "<form>"
        "<div style='clear:both; display:table;'>"
        "<div style='float:left;margin:10px;'>"
        "<table>"
        "<tr title='Read id between 0 and " << reads->readCount() - 1 << "'>"
        "<td style=\"white-space:pre-wrap; word-wrap:break-word\">"
        "Start vertex reads\n"
        "The oriented read should be in the form <code>readId-strand</code>\n"
        "where strand is 0 or 1. For example, <code>\"1345871-1</code>\".\n"
        "To add multiple start points, use a comma separator."
        "<td><input type=text required name=readId size=8 style='text-align:center'"
        << (readIdsArePresent ? ("value='" + readIdsString + "'") : "") <<
        ">"

        "<tr title='Maximum distance from start vertex (number of arcs)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance <<
        "'>"

        "<tr title='Allow reads marked as chimeric to be included in the local string graph.'>"
        "<td>Allow chimeric reads"
        "<td class=centered><input type=checkbox name=allowChimericReads" <<
        (allowChimericReads ? " checked" : "") <<
        ">"

        "<tr title='Follow outgoing arcs from each visited vertex.'>"
        "<td>Follow outgoing arcs"
        "<td class=centered><input type=checkbox name=followOutgoing" <<
        (followOutgoing || followNoneSpecified ? " checked" : "") <<
        ">"

        "<tr title='Follow incoming arcs into each visited vertex.'>"
        "<td>Follow incoming arcs"
        "<td class=centered><input type=checkbox name=followIncoming" <<
        (followIncoming ? " checked" : "") <<
        (incomingAdjacencyAvailable ? ">" : " disabled>") <<

        "<tr>"
        "<td>Layout method"
        "<td class=centered>"
        "<input type=radio required name=layoutMethod value='sfdp'" <<
        (layoutMethod == "sfdp" ? " checked=on" : "") <<
        ">sfdp"
        "<br><input type=radio required name=layoutMethod value='fdp'" <<
        (layoutMethod == "fdp" ? " checked=on" : "") <<
        ">fdp"
        "<br><input type=radio required name=layoutMethod value='neato'" <<
        (layoutMethod == "neato" ? " checked=on" : "") <<
        ">neato"

        "<tr title='Graphics size in pixels. Changing this works better than zooming.'>"
        "<td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'" <<
        " value='" << sizePixels <<
        "'>"

        "<tr>"
        "<td>Vertex scaling factor"
        "<td><input type=text required name=vertexScalingFactor size=8 style='text-align:center'" <<
        " value='" << vertexScalingFactor <<
        "'>"

        "<tr>"
        "<td>Edge thickness scaling factor"
        "<td><input type=text required name=edgeThicknessScalingFactor size=8 style='text-align:center'" <<
        " value='" << edgeThicknessScalingFactor <<
        "'>"

        "<tr title='Maximum time (in seconds) allowed for graph creation and layout'>"
        "<td>Timeout (seconds) for graph layout"
        "<td><input type=text required name=timeout size=8 style='text-align:center'" <<
        " value='" << timeout <<
        "'>"

        "</table>"
        "</div>"
        "</div>"
        "<br><input type=submit value='Display'>"
        "</form>";

    // If no start vertices are specified, default to listing some active vertices so the
    // user can click into an existing component of the cleaned graph.
    if (!readIdsArePresent && !pickActive && !listActive) {
        listActiveString = "on";
    }

    if (listActive || !listActiveString.empty()) {
        const vector<OrientedReadId> picked = pickActiveStringGraphStartVertices(activeLimit);
        html << "<h3>Active vertices</h3>";
        if (picked.empty()) {
            html << "<p>No active vertices found (string graph has no non-deleted arcs).</p>";
        } else {
            html << "<p>Showing up to " << picked.size() << " vertices found in the arc list.</p>";
            html << "<table><tr><th>Oriented read<th>Outgoing arcs";
            if (incomingAdjacencyAvailable) {
                html << "<th>Incoming arcs";
            }
            for (const OrientedReadId v : picked) {
                const uint32_t vv = uint32_t(v.getValue());
                const uint32_t outDegree = uint32_t(stringGraph.outgoing.size(vv));
                const uint32_t inDegree = incomingAdjacencyAvailable ? uint32_t(stringGraph.incoming.size(vv)) : 0;
                html << "<tr><td class=centered><a href='exploreStringGraph?readId=" << v.getString() <<
                    "&maxDistance=" << maxDistance <<
                    "&layoutMethod=" << layoutMethod <<
                    "&sizePixels=" << sizePixels <<
                    "&vertexScalingFactor=" << vertexScalingFactor <<
                    "&edgeThicknessScalingFactor=" << edgeThicknessScalingFactor <<
                    "&timeout=" << timeout <<
                    (allowChimericReads ? "&allowChimericReads=on" : "") <<
                    (followOutgoing ? "&followOutgoing=on" : "") <<
                    (followIncoming ? "&followIncoming=on" : "") <<
                    "'>" << v.getString() << "</a>"
                    "<td class=centered>" << outDegree;
                if (incomingAdjacencyAvailable) {
                    html << "<td class=centered>" << inDegree;
                }
            }
            html << "</table>";
        }
        return;
    }

    if (pickActive && !readIdsArePresent) {
        readIds = pickActiveStringGraphStartVertices(1);
        if (readIds.empty()) {
            html << "<p>No active vertices found (string graph has no non-deleted arcs).</p>";
            return;
        }
        readIdsString = readIds.front().getString();
    }

    if (!readIdsArePresent && !pickActive) {
        return;
    }
    if (!readStringsAreValid) {
        return;
    }
    for (const auto& readId : readIds) {
        if (readId.getReadId() > reads->readCount()) {
            html << "<p>Invalid read id " << readId;
            html << ". Must be between 0 and " << reads->readCount() - 1 << ".";
            return;
        }
    }

    // Default to outgoing-only if nothing is selected (more string-graph-like).
    const bool followOutgoingEffective = followOutgoing || followNoneSpecified;
    const bool followIncomingEffective = followIncoming;

    LocalStringGraph graph;
    if (!createLocalStringGraph(readIds, maxDistance, allowChimericReads,
        followOutgoingEffective, followIncomingEffective, timeout, graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    html << "<p>The local string graph has " << num_vertices(graph);
    html << " vertices and " << num_edges(graph) << " arcs.";

    if (num_edges(graph) == 0) {
        html << "<p>Start vertex has no arcs in the cleaned string graph. "
            "Try <a href='exploreStringGraph?pickActive=on'>Pick an active start vertex</a> or "
            "<a href='exploreStringGraph?listActive=on'>List some active vertices</a>.</p>";
    }


    html << "<h1 style='line-height:10px'>String graph near oriented read(s) " << readIdsString << "</h1>";

    html << R"stringDelimiter(
        <script>
        function highlight_vertex()
        {
            vertex = document.getElementById("highlight").value;
            document.getElementById("highlight").value = "";
            element = document.getElementById("Vertex-" + vertex);
            element.setAttribute("fill", "#ff00ff");
        }
        </script>
        <p>
        <input id=highlight type=text onchange="highlight_vertex()" size=10>
        Enter an oriented read to highlight, then press Enter. The oriented read should be
        in the form <code>readId-strand</code> where strand is 0 or 1 (for example, <code>"1345871-1</code>").
        To highlight multiple oriented reads, enter them one at a time in the same way.
        <p>
        )stringDelimiter";

    addScaleSvgButtons(html, sizePixels);

    const ComputeLayoutReturnCode returnCode = graph.computeLayout(layoutMethod, timeout);
    if (returnCode == ComputeLayoutReturnCode::Timeout) {
        html << "<p>Timeout exceeded for computing graph layout. Try longer timeout or different parameters.</p>";
        return;
    }
    if (returnCode != ComputeLayoutReturnCode::Success) {
        html << "<p>ERROR: graph layout failed </p>";
        return;
    }

    graph.writeSvg("svg", sizePixels, sizePixels, vertexScalingFactor, edgeThicknessScalingFactor,
        maxDistance, *this, html);
}
