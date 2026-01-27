// Dinara.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "LocalUnitigGraph.hpp"
#include "timestamp.hpp"
using namespace dinara;

// Standard library.
#include "chrono.hpp"
#include <unordered_set>


namespace {
    bool parseCommaSeparatedUnitigIDs(
        const std::string& commaSeparated,
        std::vector<OrientedUnitigId>& ids,
        std::ostream& html)
    {
        ids.clear();
        std::string token;
        for (char c : commaSeparated) {
            if (c == ',') {
                if (!token.empty()) {
                    try {
                        ids.emplace_back(token);
                    } catch (...) {
                        html << "<p>Invalid oriented unitig id: '" << token << "'</p>";
                        html << "<p>Specify one or more comma separated oriented unitig ids of the form "
                            "<code>unitigId-strand</code> where strand is 0 or 1.</p>";
                        return false;
                    }
                    token.clear();
                }
            } else {
                token += c;
            }
        }
        if (!token.empty()) {
            try {
                ids.emplace_back(token);
            } catch (...) {
                html << "<p>Invalid oriented unitig id: '" << token << "'</p>";
                html << "<p>Specify one or more comma separated oriented unitig ids of the form "
                    "<code>unitigId-strand</code> where strand is 0 or 1.</p>";
                return false;
            }
        }
        return true;
    }
}



void Assembler::exploreUnitigGraph(
    const std::vector<std::string>& request,
    std::ostream& html)
{
    try {
        checkUnitigGraphIsOpen();
    } catch (...) {
        html << "The unitig graph is not available." << std::endl;
        return;
    }

    // Parameters.
    std::vector<OrientedUnitigId> unitigIds;
    std::string unitigIdsString;
    const bool unitigIdsArePresent = getParameterValue(request, "unitigId", unitigIdsString);
    const bool unitigStringsAreValid = parseCommaSeparatedUnitigIDs(unitigIdsString, unitigIds, html);

    std::string pickActiveString;
    const bool pickActive = getParameterValue(request, "pickActive", pickActiveString);

    std::string listActiveString;
    const bool listActive = getParameterValue(request, "listActive", listActiveString);
    uint32_t activeLimit = 30;
    getParameterValue(request, "activeLimit", activeLimit);

    uint32_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    std::string followOutgoingString;
    const bool followOutgoing = getParameterValue(request, "followOutgoing", followOutgoingString);
    std::string followIncomingString;
    const bool followIncoming = getParameterValue(request, "followIncoming", followIncomingString);
    const bool followNoneSpecified = (!followOutgoing && !followIncoming);

    std::string layoutMethod = "sfdp";
    getParameterValue(request, "layoutMethod", layoutMethod);

    uint32_t sizePixels = 600;
    getParameterValue(request, "sizePixels", sizePixels);

    double vertexScalingFactor = 1.;
    getParameterValue(request, "vertexScalingFactor", vertexScalingFactor);

    double edgeThicknessScalingFactor = 1.;
    getParameterValue(request, "edgeThicknessScalingFactor", edgeThicknessScalingFactor);

    double timeout = 30;
    getParameterValue(request, "timeout", timeout);

    const auto pickActiveStartVertices = [&](uint32_t limit) -> std::vector<OrientedUnitigId> {
        std::vector<OrientedUnitigId> result;
        if (!unitigGraph.arcs.isOpen) return result;
        result.reserve(limit);
        std::unordered_set<uint32_t> seen;
        seen.reserve(limit * 2);

        const uint64_t arcCount = unitigGraph.arcs.size();
        for (uint64_t arcId = 0; arcId < arcCount && result.size() < limit; ++arcId) {
            const auto& arc = unitigGraph.arcs[arcId];
            if (arc.del) continue;

            const auto considerVertex = [&](uint32_t v) {
                if (result.size() >= limit) return;
                if (!seen.insert(v).second) return;
                const OrientedUnitigId u = OrientedUnitigId::fromValue(v);
                if (unitigGraph.unitigDeleted.isOpen &&
                    u.getUnitigId() < unitigGraph.unitigDeleted.size() &&
                    unitigGraph.unitigDeleted[u.getUnitigId()]) {
                    return;
                }
                result.push_back(u);
            };

            considerVertex(arc.from);
            considerVertex(arc.to);
        }
        return result;
    };

    // Form.
    html <<
        "<h3>Display a local subgraph of the unitig graph</h3>"
        "<p>"
        "<a href='exploreUnitigGraph?pickActive=on'>Pick an active start vertex</a>"
        "&nbsp;&nbsp;"
        "<a href='exploreUnitigGraph?listActive=on'>List some active vertices</a>"
        "</p>"
        "<form>"
        "<div style='clear:both; display:table;'>"
        "<div style='float:left;margin:10px;'>"
        "<table>"
        "<tr title='Unitig id between 0 and " << (unitigGraph.unitigs.size() ? unitigGraph.unitigs.size() - 1 : 0) << "'>"
        "<td style=\"white-space:pre-wrap; word-wrap:break-word\">"
        "Start unitig(s)\n"
        "The oriented unitig should be in the form <code>unitigId-strand</code>\n"
        "where strand is 0 or 1. For example, <code>\"42-1</code>\".\n"
        "To add multiple start points, use a comma separator."
        "<td><input type=text required name=unitigId size=10 style='text-align:center'"
        << (unitigIdsArePresent ? ("value='" + unitigIdsString + "'") : "") <<
        ">"

        "<tr title='Maximum distance from start vertex (number of arcs)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance <<
        "'>"

        "<tr title='Follow outgoing arcs from each visited vertex.'>"
        "<td>Follow outgoing arcs"
        "<td class=centered><input type=checkbox name=followOutgoing" <<
        (followOutgoing || followNoneSpecified ? " checked" : "") <<
        ">"

        "<tr title='Follow incoming arcs into each visited vertex.'>"
        "<td>Follow incoming arcs"
        "<td class=centered><input type=checkbox name=followIncoming" <<
        (followIncoming ? " checked" : "") <<
        ">"

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

    // If no start vertices are specified, default to listing active vertices.
    if (!unitigIdsArePresent && !pickActive && !listActive) {
        listActiveString = "on";
    }

    if (listActive || !listActiveString.empty()) {
        const std::vector<OrientedUnitigId> picked = pickActiveStartVertices(activeLimit);
        html << "<h3>Active vertices</h3>";
        if (picked.empty()) {
            html << "<p>No active vertices found (unitig graph has no non-deleted arcs).</p>";
        } else {
            html << "<p>Showing up to " << picked.size() << " vertices found in the arc list.</p>";
            html << "<table><tr><th>Oriented unitig<th>Outgoing arcs<th>Incoming arcs";
            for (const OrientedUnitigId u : picked) {
                const uint32_t vv = u.getValue();
                const uint32_t outDegree = uint32_t(unitigGraph.outgoing.size(vv));
                const uint32_t inDegree = uint32_t(unitigGraph.incoming.size(vv));
                html << "<tr><td class=centered><a href='exploreUnitigGraph?unitigId=" << u.getString() <<
                    "&maxDistance=" << maxDistance <<
                    "&layoutMethod=" << layoutMethod <<
                    "&sizePixels=" << sizePixels <<
                    "&vertexScalingFactor=" << vertexScalingFactor <<
                    "&edgeThicknessScalingFactor=" << edgeThicknessScalingFactor <<
                    "&timeout=" << timeout <<
                    (followOutgoing ? "&followOutgoing=on" : "") <<
                    (followIncoming ? "&followIncoming=on" : "") <<
                    "'>" << u.getString() << "</a>"
                    "<td class=centered>" << outDegree <<
                    "<td class=centered>" << inDegree;
            }
            html << "</table>";
        }
        return;
    }

    if (pickActive && !unitigIdsArePresent) {
        unitigIds = pickActiveStartVertices(1);
        if (unitigIds.empty()) {
            html << "<p>No active vertices found (unitig graph has no non-deleted arcs).</p>";
            return;
        }
        unitigIdsString = unitigIds.front().getString();
    }

    if (!unitigIdsArePresent && !pickActive) return;
    if (!unitigStringsAreValid) return;

    for (const auto& u : unitigIds) {
        if (u.getUnitigId() >= unitigGraph.unitigs.size()) {
            html << "<p>Invalid unitig id " << u.getUnitigId() << ".";
            return;
        }
    }

    const bool followOutgoingEffective = followOutgoing || followNoneSpecified;
    const bool followIncomingEffective = followIncoming;

    LocalUnitigGraph graph;
    if (!createLocalUnitigGraph(unitigIds, maxDistance, followOutgoingEffective, followIncomingEffective, timeout, graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.</p>";
        return;
    }
    html << "<p>The local unitig graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " arcs.</p>";

    html << "<h1 style='line-height:10px'>Unitig graph near " << unitigIdsString << "</h1>";

    addScaleSvgButtons(html, sizePixels);

    const ComputeLayoutReturnCode returnCode = graph.computeLayout(layoutMethod, timeout);
    if (returnCode == ComputeLayoutReturnCode::Timeout) {
        html << "<p>Timeout exceeded for computing graph layout.</p>";
        return;
    }
    if (returnCode != ComputeLayoutReturnCode::Success) {
        html << "<p>ERROR: graph layout failed.</p>";
        return;
    }

    graph.writeSvg("svg", sizePixels, sizePixels, vertexScalingFactor, edgeThicknessScalingFactor,
        maxDistance, *this, html);
}

