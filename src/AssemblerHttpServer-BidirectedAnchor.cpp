// HTTP server handlers for BRG-native anchor visualization.

#include "Assembler.hpp"
#include "mode3-BidirectedAnchor.hpp"
#include "mode3-LocalAnchorGraph.hpp"
#include "HttpServer.hpp"
#include "invalid.hpp"
#include "Marker.hpp"
#include "Reads.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include <cctype>
#include <cmath>
#include <iomanip>
#include <queue>
#include <set>
#include <unordered_map>

using namespace dinara;
using namespace mode3;
using namespace std;

using LocalBidirectedEdge = mode3::LocalBidirectedEdge;

namespace {

uint64_t orientedId(const OrientedBidirectedAnchor oa)
{
    return 2 * oa.anchorId + uint64_t(oa.strand);
}

OrientedBidirectedAnchor orientedAnchorFromId(const uint64_t nodeId)
{
    return OrientedBidirectedAnchor(nodeId / 2, Strand(nodeId & 1));
}

string orientedAnchorToString(const OrientedBidirectedAnchor oa)
{
    return string(1, oa.strand == 0 ? '>' : '<') + to_string(oa.anchorId);
}

bool parseOrientedAnchorToken(
    const string& token,
    OrientedBidirectedAnchor& oa)
{
    if(token.empty()) {
        return false;
    }

    string s = token;
    boost::trim(s);
    if(s.empty()) {
        return false;
    }

    Strand strand = 0;
    if(s.front() == '>' || s.front() == '<') {
        strand = (s.front() == '>') ? 0 : 1;
        s.erase(s.begin());
    } else if(s.back() == '+' || s.back() == '-') {
        strand = (s.back() == '+') ? 0 : 1;
        s.pop_back();
    }

    if(s.empty()) {
        return false;
    }

    for(const char c : s) {
        if(!std::isdigit(static_cast<unsigned char>(c))) {
            return false;
        }
    }

    try {
        oa = OrientedBidirectedAnchor(stoull(s), strand);
    } catch(...) {
        return false;
    }
    return true;
}

double bidirectedEdgeCoverageLoss(
    const BidirectedAnchors& graph,
    const BidirectedEdge& edge)
{
    const uint64_t common = graph.countCommon(edge.from.anchorId, edge.to.anchorId);
    if(common == 0) {
        return 1.;
    }
    const double loss = double(common) - double(edge.coverage);
    return loss <= 0. ? 0. : loss / double(common);
}

// BFS on canonical BidirectedAnchorIds.
// Each vertex is a single canonical anchor (not doubled by strand).
// At each anchor, out-edges from both strands are followed to discover neighbors.
// Edges are deduplicated per canonical anchor pair, using max coverage.
void computeLocalBidirectedSubgraph(
    const BidirectedAnchors& graph,
    const vector<BidirectedAnchorId>& seeds,
    const uint64_t maxDistance,
    const uint64_t minCoverage,
    const bool filterEdgesByCoverageLoss,
    const double maxCoverageLoss,
    unordered_map<uint64_t, uint64_t>& nodeDistance,
    vector<LocalBidirectedEdge>& edges)
{
    nodeDistance.clear();
    edges.clear();

    queue<uint64_t> q;
    for(const BidirectedAnchorId seed : seeds) {
        if(nodeDistance.contains(seed)) {
            continue;
        }
        nodeDistance.insert({seed, 0});
        q.push(seed);
    }

    // Safety limit to avoid accidental huge layouts from broad queries.
    const uint64_t maxNodeCount = 5000;

    auto keepEdge = [&](const BidirectedEdge& e) {
        if(e.coverage < minCoverage) {
            return false;
        }
        if(filterEdgesByCoverageLoss &&
            bidirectedEdgeCoverageLoss(graph, e) > maxCoverageLoss) {
            return false;
        }
        return true;
    };

    while(!q.empty()) {
        const BidirectedAnchorId anchorId0 = q.front();
        q.pop();

        const uint64_t distance0 = nodeDistance.at(anchorId0);
        if(distance0 >= maxDistance) {
            continue;
        }

        // Follow out-edges from both strands to discover canonical neighbors.
        // By the symmetry of the bidirected edge table, out-edges from both
        // strands cover all connectivity (forward and backward).
        for(Strand s = 0; s <= 1; s++) {
            const OrientedBidirectedAnchor oa(anchorId0, s);
            for(const BidirectedEdge& e : graph.getOutEdges(oa)) {
                if(!keepEdge(e)) {
                    continue;
                }
                const BidirectedAnchorId neighborId = e.to.anchorId;
                if(nodeDistance.contains(neighborId)) {
                    continue;
                }
                nodeDistance.insert({neighborId, distance0 + 1});
                if(nodeDistance.size() <= maxNodeCount) {
                    q.push(neighborId);
                }
            }
        }

        if(nodeDistance.size() > maxNodeCount) {
            break;
        }
    }

    // Collect edges between canonical anchors in the subgraph.
    // Use canonical pair (min, max) as key to deduplicate symmetric edges.
    // Use max coverage across all strand combinations.
    map<pair<uint64_t, uint64_t>, uint64_t> edgeCoverage;
    for(const auto& [anchorId0, distance] : nodeDistance) {
        (void)distance;
        for(Strand s = 0; s <= 1; s++) {
            const OrientedBidirectedAnchor oa(anchorId0, s);
            for(const BidirectedEdge& e : graph.getOutEdges(oa)) {
                if(e.coverage < minCoverage) {
                    continue;
                }
                if(filterEdgesByCoverageLoss &&
                    bidirectedEdgeCoverageLoss(graph, e) > maxCoverageLoss) {
                    continue;
                }
                const BidirectedAnchorId anchorId1 = e.to.anchorId;
                if(!nodeDistance.contains(anchorId1)) {
                    continue;
                }
                const auto key = make_pair(
                    min(anchorId0, anchorId1),
                    max(anchorId0, anchorId1));
                auto& cov = edgeCoverage[key];
                cov = max(cov, e.coverage);
            }
        }
    }

    edges.reserve(edgeCoverage.size());
    for(const auto& [p, cov] : edgeCoverage) {
        LocalBidirectedEdge edge;
        edge.from = p.first;
        edge.to = p.second;
        edge.coverage = cov;
        edges.push_back(edge);
    }
}

} // namespace



void Assembler::exploreBidirectedAnchor(const vector<string>& request, ostream& html)
{
    if(!bidirectedAnchors) {
        html << "<p>BRG-native anchors are not available.";
        return;
    }

    // Get the request parameters.
    uint64_t anchorIdValue = invalid<uint64_t>;
    const bool anchorIdIsPresent =
        HttpServer::getParameterValue(request, "anchorId", anchorIdValue);

    // Begin the form.
    html <<
        "<h2>BRG Anchor information</h2>"
        "<p>BRG anchors are bidirectional (no forward/reverse complement doubling). "
        "Each anchor is identified by a simple numeric id. "
        "An <i>oriented anchor</i> is (anchorId, strand), representing one end of "
        "the bidirected node."
        "<form>"
        "<table>"

        "<tr><th class=left>BRG Anchor id"
        "<td class=centered>"
        "<input type=text name=anchorId required style='text-align:center'";
    if(anchorIdIsPresent) {
        html << " value='" << anchorIdValue << "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        bidirectedAnchors->size() - 1 << ".'>"

        "</table>"
        "<input type=submit value='Show BRG anchor details'> "
        "</form>";

    // If the anchor id is missing or invalid, stop here.
    if(!anchorIdIsPresent) {
        return;
    }

    const BidirectedAnchorId anchorId = anchorIdValue;
    if(anchorId >= bidirectedAnchors->size()) {
        html << "<p>Invalid BRG anchor id. Must be between 0 and " <<
            bidirectedAnchors->size() - 1 << ".";
        return;
    }


    html << "<h1>BRG Anchor " << anchorId << "</h1>";

    const auto intervals = (*bidirectedAnchors)[anchorId];
    const uint64_t cov = intervals.size();

    // Use OrientedBidirectedAnchor for the two directions.
    const OrientedBidirectedAnchor oa0(anchorId, 0);
    const OrientedBidirectedAnchor oa1(anchorId, 1);

    // Forward/backward neighbors using edge table.
    vector<BidirectedAnchors::Neighbor> fwdNeighbors0, fwdNeighbors1;
    bidirectedAnchors->findForwardNeighbors(oa0, fwdNeighbors0);
    bidirectedAnchors->findForwardNeighbors(oa1, fwdNeighbors1);

    vector<BidirectedAnchors::Neighbor> bwdNeighbors0, bwdNeighbors1;
    bidirectedAnchors->findBackwardNeighbors(oa0, bwdNeighbors0);
    bidirectedAnchors->findBackwardNeighbors(oa1, bwdNeighbors1);


    // Summary table.
    html <<
        "<table>"
        "<tr><th class=left>Coverage<td class=centered>" << cov <<
        "<tr><th class=left>Forward neighbors (strand 0)<td class=centered>" << fwdNeighbors0.size() <<
        "<tr><th class=left>Forward neighbors (strand 1)<td class=centered>" << fwdNeighbors1.size() <<
        "<tr><th class=left>Backward neighbors (strand 0)<td class=centered>" << bwdNeighbors0.size() <<
        "<tr><th class=left>Backward neighbors (strand 1)<td class=centered>" << bwdNeighbors1.size() <<
        "</table>";


    // Neighbor tables.
    auto writeNeighborTable = [&](const string& title, const vector<BidirectedAnchors::Neighbor>& neighbors) {
        html << "<h2>" << title << "</h2>";
        if(neighbors.empty()) {
            html << "<p>None.";
            return;
        }
        html <<
            "<table>"
            "<tr><th>Anchor<th>Strand<th>Coverage";
        for(const auto& n : neighbors) {
            html <<
                "<tr>"
                "<td class=centered><a href='exploreBidirectedAnchor?anchorId=" << n.anchorId << "'>"
                << n.anchorId << "</a>"
                "<td class=centered>" << uint32_t(n.strand) <<
                "<td class=centered>" << n.count;
        }
        html << "</table>";
    };

    writeNeighborTable("Backward neighbors (strand 0 reads)", bwdNeighbors0);
    writeNeighborTable("Forward neighbors (strand 0 reads)", fwdNeighbors0);
    writeNeighborTable("Backward neighbors (strand 1 reads)", bwdNeighbors1);
    writeNeighborTable("Forward neighbors (strand 1 reads)", fwdNeighbors1);


    // Marker intervals table.
    html <<
        "<h2>Marker intervals</h2>"
        "<table>"
        "<tr>"
        "<th>Index"
        "<th>Read id"
        "<th>Strand"
        "<th>Ordinal<br>(strand 0)"
        "<th>Position"
        "<th>Previous<br>anchor<br>in journey"
        "<th>Next<br>anchor<br>in journey";

    for(uint64_t i = 0; i < cov; i++) {
        const BidirectedAnchorMarkerInterval& interval = intervals[i];
        const ReadId readId = interval.readId;
        const uint32_t ordinal = interval.ordinal;
        const Strand strand = interval.strand;

        // Get position from strand-0 markers.
        const OrientedReadId orientedReadId0(readId, 0);
        const auto readMarkers = (*markers)[orientedReadId0.getValue()];
        const uint32_t position = readMarkers[ordinal].position;

        // Find this anchor's position in the journey and get prev/next.
        const auto jrn = bidirectedAnchors->journey(readId);
        OrientedBidirectedAnchor prevOa;
        OrientedBidirectedAnchor nextOa;
        for(uint64_t j = 0; j < jrn.size(); j++) {
            if(jrn[j].anchorId == anchorId) {
                if(j > 0) {
                    prevOa = jrn[j - 1].oriented();
                }
                if(j + 1 < jrn.size()) {
                    nextOa = jrn[j + 1].oriented();
                }
                break;
            }
        }

        html <<
            "<tr>"
            "<td class=centered>" << i <<
            "<td class=centered><a href='exploreBidirectedJourney?readId=" << readId << "'>"
            << readId << "</a>"
            "<td class=centered>" << uint32_t(strand) <<
            "<td class=centered>" << ordinal <<
            "<td class=centered>" << position;

        // Previous anchor in journey.
        html << "<td class=centered>";
        if(prevOa.isValid()) {
            html << "<a href='exploreBidirectedAnchor?anchorId=" << prevOa.anchorId << "'>"
                << prevOa.anchorId << "</a>"
                " (s" << uint32_t(prevOa.strand) << ")";
        }

        // Next anchor in journey.
        html << "<td class=centered>";
        if(nextOa.isValid()) {
            html << "<a href='exploreBidirectedAnchor?anchorId=" << nextOa.anchorId << "'>"
                << nextOa.anchorId << "</a>"
                " (s" << uint32_t(nextOa.strand) << ")";
        }
    }
    html << "</table>";


    // Tangle matrix — now keyed by (OrientedBidirectedAnchor prev, OrientedBidirectedAnchor next).
    map<pair<OrientedBidirectedAnchor, OrientedBidirectedAnchor>, uint64_t> tangleMatrix;
    for(uint64_t i = 0; i < cov; i++) {
        const BidirectedAnchorMarkerInterval& interval = intervals[i];
        const auto jrn = bidirectedAnchors->journey(interval.readId);
        OrientedBidirectedAnchor prev;
        OrientedBidirectedAnchor next;
        for(uint64_t j = 0; j < jrn.size(); j++) {
            if(jrn[j].anchorId == anchorId) {
                if(j > 0) prev = jrn[j - 1].oriented();
                if(j + 1 < jrn.size()) next = jrn[j + 1].oriented();
                break;
            }
        }
        ++tangleMatrix[make_pair(prev, next)];
    }

    html << "<h2>Tangle matrix</h2>"
        "<table><tr><th>In<th>In strand<th>Out<th>Out strand<th>Coverage";
    for(const auto& [key, count] : tangleMatrix) {
        const auto& [prev, next] = key;
        html << "<tr><td class=centered>";
        if(prev.isValid()) {
            html << "<a href='exploreBidirectedAnchor?anchorId=" << prev.anchorId << "'>"
                << prev.anchorId << "</a>";
        }
        html << "<td class=centered>";
        if(prev.isValid()) {
            html << uint32_t(prev.strand);
        }
        html << "<td class=centered>";
        if(next.isValid()) {
            html << "<a href='exploreBidirectedAnchor?anchorId=" << next.anchorId << "'>"
                << next.anchorId << "</a>";
        }
        html << "<td class=centered>";
        if(next.isValid()) {
            html << uint32_t(next.strand);
        }
        html << "<td class=centered>" << count;
    }
    html << "</table>";
}


void Assembler::exploreBidirectedAnchorGraph(const vector<string>& request, ostream& html)
{
    if(!bidirectedAnchors) {
        html << "<p>BRG-native anchors are not available.";
        return;
    }

    const uint64_t anchorCount = bidirectedAnchors->size();
    const uint64_t orientedNodeCount = 2 * anchorCount;

    string anchorIdsString;
    HttpServer::getParameterValue(request, "anchorIdsString", anchorIdsString);
    boost::trim(anchorIdsString);

    uint64_t distance = 10;
    HttpServer::getParameterValue(request, "distance", distance);

    uint64_t minCoverage = 0;
    HttpServer::getParameterValue(request, "minCoverage", minCoverage);

    string filterEdgesByCoverageLossString;
    const bool filterEdgesByCoverageLoss = HttpServer::getParameterValue(
        request, "filterEdgesByCoverageLoss", filterEdgesByCoverageLossString);
    double maxCoverageLoss = 0.5;
    HttpServer::getParameterValue(request, "maxCoverageLoss", maxCoverageLoss);

    const LocalAnchorGraphDisplayOptions displayOptions(request);

    html <<
        "<h2>Local Bidirected Anchor Graph</h2>"
        "<p>Enter starting anchor ids to display a local subgraph. "
        "Each anchor is a canonical bidirected node (like a GFA segment). "
        "Enter numeric ids, comma or space separated."
        "<form><table>"

        "<tr><th class=left>Starting anchor ids"
        "<td class=centered><input type=text name=anchorIdsString style='text-align:center' required";
    if(!anchorIdsString.empty()) {
        html << " value='" << anchorIdsString << "'";
    }
    html <<
        " size=32 title='Enter anchor ids, for example: 120,455,900'>"

        "<tr><th class=left>Distance"
        "<td class=centered><input type=text name=distance style='text-align:center' required size=8 value='" <<
        distance << "'>"

        "<tr><th>Edge filtering"
        "<td>"
        "Minimum coverage "
        "<input type=text name=minCoverage style='text-align:center' required size=8 value='" <<
        minCoverage << "'>"
        "<br>"
        "<input type=checkbox name=filterEdgesByCoverageLoss" <<
        (filterEdgesByCoverageLoss ? " checked" : "") <<
        ">Filter edges by coverage loss"
        "<br><input type=text name=maxCoverageLoss style='text-align:center' required size=6 value=" <<
        maxCoverageLoss << "> Maximum coverage loss"
        "<hr>";

    displayOptions.writeForm(html);

    html <<
        "</table>"
        "<input type=submit value='Display local bidirected anchor graph'>"
        "</form>";

    if(!anchorIdsString.empty()) {
        vector<BidirectedAnchorId> seedAnchors;
        {
            boost::tokenizer<boost::char_separator<char>> tokenizer(
                anchorIdsString, boost::char_separator<char>(", "));
            for(const string& token : tokenizer) {
                string s = token;
                boost::trim(s);
                if(s.empty()) {
                    continue;
                }
                BidirectedAnchorId anchorId;
                try {
                    anchorId = stoull(s);
                } catch(...) {
                    html << "<p>Invalid anchor id: " << token << ".";
                    html << " Enter plain numeric anchor ids.";
                    return;
                }
                if(anchorId >= anchorCount) {
                    html << "<p>Invalid anchor id " << anchorId
                        << ". Must be in [0, " <<
                        (anchorCount == 0 ? 0 : anchorCount - 1) << "].";
                    return;
                }
                seedAnchors.push_back(anchorId);
            }
        }

        if(seedAnchors.empty()) {
            html << "<p>No starting anchors were provided.";
        } else {
            sort(seedAnchors.begin(), seedAnchors.end());
            seedAnchors.erase(
                unique(seedAnchors.begin(), seedAnchors.end()),
                seedAnchors.end());

            unordered_map<uint64_t, uint64_t> nodeDistance;
            vector<LocalBidirectedEdge> localEdges;
            computeLocalBidirectedSubgraph(
                *bidirectedAnchors,
                seedAnchors,
                distance,
                minCoverage,
                filterEdgesByCoverageLoss,
                maxCoverageLoss,
                nodeDistance,
                localEdges);

            // Create a LocalAnchorGraph from the BidirectedAnchors data
            // and use the same rendering pipeline as the AnchorGraph.
            LocalAnchorGraph graph(
                *bidirectedAnchors,
                nodeDistance,
                localEdges,
                distance);

            html << "<h1>Local bidirected anchor graph</h1>"
                "<p>The local bidirected anchor graph has "
                << num_vertices(graph) << " vertices and "
                << num_edges(graph) << " edges.";
            if(nodeDistance.size() >= 5000) {
                html << " (Reached safety limit of 5000 nodes.)";
            }

            graph.writeHtml(html, displayOptions);
        }
    }

    html << "<h1>Bidirected Anchor Graph Summary</h1>";

    uint64_t directedEdgeCount = 0;
    uint64_t activePathCount = 0;

    map<uint64_t, uint64_t> outDegreeDist;
    map<uint64_t, uint64_t> inDegreeDist;
    map<uint64_t, uint64_t> coverageDist;

    for(BidirectedAnchorId anchorId = 0; anchorId < anchorCount; ++anchorId) {
        const uint64_t cov = bidirectedAnchors->coverage(anchorId);
        coverageDist[cov]++;

        for(Strand strand = 0; strand <= 1; ++strand) {
            const OrientedBidirectedAnchor oa(anchorId, strand);
            const uint64_t outDeg = bidirectedAnchors->getOutEdges(oa).size();
            const uint64_t inDeg = bidirectedAnchors->getInEdges(oa).size();
            outDegreeDist[outDeg]++;
            inDegreeDist[inDeg]++;
            directedEdgeCount += outDeg;
        }
    }

    for(ReadId readId = 0; readId < reads->readCount(); ++readId) {
        if(!bidirectedAnchors->journey(readId).empty()) {
            ++activePathCount;
        }
    }

    html <<
        "<table>"
        "<tr><th class=left>Anchors<td class=centered>" << anchorCount <<
        "<tr><th class=left>Oriented nodes<td class=centered>" << orientedNodeCount <<
        "<tr><th class=left>Directed edges<td class=centered>" << directedEdgeCount <<
        "<tr><th class=left>Active paths (reads with non-empty journeys)<td class=centered>" << activePathCount <<
        "</table>";

    html << "<h2>Out-degree distribution (oriented nodes)</h2>"
        "<table><tr><th>Degree<th>Count";
    for(const auto& [deg, cnt] : outDegreeDist) {
        html << "<tr><td class=centered>" << deg << "<td class=centered>" << cnt;
    }
    html << "</table>";

    html << "<h2>In-degree distribution (oriented nodes)</h2>"
        "<table><tr><th>Degree<th>Count";
    for(const auto& [deg, cnt] : inDegreeDist) {
        html << "<tr><td class=centered>" << deg << "<td class=centered>" << cnt;
    }
    html << "</table>";

    html << "<h2>Coverage distribution (anchors)</h2>"
        "<table><tr><th>Coverage<th>Count";
    for(const auto& [cov, cnt] : coverageDist) {
        html << "<tr><td class=centered>" << cov << "<td class=centered>" << cnt;
    }
    html << "</table>";

    const uint64_t maxList = 200;
    const uint64_t showCount = min(anchorCount, maxList);

    html << "<h2>Anchors";
    if(anchorCount > maxList) {
        html << " (first " << maxList << " of " << anchorCount << ")";
    }
    html << "</h2>"
        "<table>"
        "<tr><th>Anchor<th>Coverage"
        "<th>Out-deg (s0)<th>In-deg (s0)"
        "<th>Out-deg (s1)<th>In-deg (s1)"
        "<th>Paths crossing";

    for(uint64_t anchorId = 0; anchorId < showCount; ++anchorId) {
        const uint64_t cov = bidirectedAnchors->coverage(anchorId);
        const uint64_t out0 = bidirectedAnchors->getOutEdges(OrientedBidirectedAnchor(anchorId, 0)).size();
        const uint64_t in0 = bidirectedAnchors->getInEdges(OrientedBidirectedAnchor(anchorId, 0)).size();
        const uint64_t out1 = bidirectedAnchors->getOutEdges(OrientedBidirectedAnchor(anchorId, 1)).size();
        const uint64_t in1 = bidirectedAnchors->getInEdges(OrientedBidirectedAnchor(anchorId, 1)).size();

        html <<
            "<tr>"
            "<td class=centered><a href='exploreBidirectedAnchorGraphNode?nodeId=" << anchorId << "'>"
            << anchorId << "</a>"
            "<td class=centered>" << cov <<
            "<td class=centered>" << out0 <<
            "<td class=centered>" << in0 <<
            "<td class=centered>" << out1 <<
            "<td class=centered>" << in1 <<
            "<td class=centered>" << cov;
    }
    html << "</table>";
}


void Assembler::exploreBidirectedAnchorGraphNode(const vector<string>& request, ostream& html)
{
    if(!bidirectedAnchors) {
        html << "<p>The Bidirected Anchor Graph is not available.";
        return;
    }

    uint64_t nodeIdValue = invalid<uint64_t>;
    const bool nodeIdIsPresent =
        HttpServer::getParameterValue(request, "nodeId", nodeIdValue);

    html <<
        "<h2>Bidirected Anchor Graph Node</h2>"
        "<p>Each node is one BRG anchor id, with two oriented versions: "
        "(nodeId, strand 0) and (nodeId, strand 1)."
        "<form>"
        "<table>"
        "<tr><th class=left>Node id"
        "<td class=centered>"
        "<input type=text name=nodeId required style='text-align:center'";
    if(nodeIdIsPresent) {
        html << " value='" << nodeIdValue << "'";
    }
    html <<
        " size=8 title='Enter an anchor id.'>"
        "</table>"
        "<input type=submit value='Show node details'> "
        "</form>";

    if(!nodeIdIsPresent) {
        return;
    }

    const BidirectedAnchorId anchorId = nodeIdValue;
    if(anchorId >= bidirectedAnchors->size()) {
        html << "<p>Node " << anchorId << " does not exist.";
        return;
    }

    const uint64_t cov = bidirectedAnchors->coverage(anchorId);
    const OrientedBidirectedAnchor oa0(anchorId, 0);
    const OrientedBidirectedAnchor oa1(anchorId, 1);

    html << "<h1>Bidirected Anchor Graph Node " << anchorId << "</h1>";
    html <<
        "<table>"
        "<tr><th class=left>Coverage<td class=centered>" << cov <<
        "<tr><th class=left>Out-degree (strand 0)<td class=centered>" << bidirectedAnchors->getOutEdges(oa0).size() <<
        "<tr><th class=left>In-degree (strand 0)<td class=centered>" << bidirectedAnchors->getInEdges(oa0).size() <<
        "<tr><th class=left>Out-degree (strand 1)<td class=centered>" << bidirectedAnchors->getOutEdges(oa1).size() <<
        "<tr><th class=left>In-degree (strand 1)<td class=centered>" << bidirectedAnchors->getInEdges(oa1).size() <<
        "</table>";

    auto writeEdgeTable = [&](const string& title, const span<const BidirectedEdge>& edges, bool outgoing) {
        html << "<h2>" << title << "</h2>";
        if(edges.empty()) {
            html << "<p>None.";
            return;
        }
        html << "<table><tr><th>Anchor<th>Strand<th>Coverage";
        for(const auto& e : edges) {
            const OrientedBidirectedAnchor other = outgoing ? e.to : e.from;
            html <<
                "<tr>"
                "<td class=centered><a href='exploreBidirectedAnchorGraphNode?nodeId=" << other.anchorId << "'>"
                << other.anchorId << "</a>"
                "<td class=centered>" << uint32_t(other.strand) <<
                "<td class=centered>" << e.coverage;
        }
        html << "</table>";
    };

    writeEdgeTable("Outgoing edges (strand 0)", bidirectedAnchors->getOutEdges(oa0), true);
    writeEdgeTable("Incoming edges (strand 0)", bidirectedAnchors->getInEdges(oa0), false);
    writeEdgeTable("Outgoing edges (strand 1)", bidirectedAnchors->getOutEdges(oa1), true);
    writeEdgeTable("Incoming edges (strand 1)", bidirectedAnchors->getInEdges(oa1), false);

    const auto intervals = (*bidirectedAnchors)[anchorId];
    vector<ReadId> readIds;
    readIds.reserve(intervals.size());
    for(const auto& interval : intervals) {
        readIds.push_back(interval.readId);
    }
    sort(readIds.begin(), readIds.end());
    readIds.erase(unique(readIds.begin(), readIds.end()), readIds.end());

    html << "<h2>Paths crossing this node (" << readIds.size() << ")</h2>";
    if(readIds.empty()) {
        html << "<p>No paths cross this node.";
        return;
    }

    const uint64_t maxShow = 100;
    const uint64_t showCount = min(uint64_t(readIds.size()), maxShow);
    html <<
        "<table>"
        "<tr><th>Path<th>Length<th>Context (5 nodes around this node)";

    for(uint64_t i = 0; i < showCount; ++i) {
        const ReadId readId = readIds[i];
        const auto path = bidirectedAnchors->journey(readId);
        if(path.empty()) {
            continue;
        }

        int64_t posInPath = -1;
        for(uint64_t j = 0; j < path.size(); ++j) {
            if(path[j].anchorId == anchorId) {
                posInPath = int64_t(j);
                break;
            }
        }

        html <<
            "<tr>"
            "<td class=centered><a href='exploreBidirectedAnchorGraphPath?pathId=" << readId << "'>"
            << readId << "</a>"
            "<td class=centered>" << path.size() <<
            "<td class=centered>";

        if(posInPath >= 0) {
            const int64_t radius = 2;
            const int64_t start = max(int64_t(0), posInPath - radius);
            const int64_t end = min(int64_t(path.size()) - 1, posInPath + radius);
            if(start > 0) {
                html << "... ";
            }
            for(int64_t j = start; j <= end; ++j) {
                if(j > start) {
                    html << " &rarr; ";
                }
                const auto& step = path[j];
                if(j == posInPath) {
                    html << "<b>";
                }
                html << (step.strand == 0 ? "&gt;" : "&lt;") << step.anchorId;
                if(j == posInPath) {
                    html << "</b>";
                }
            }
            if(end < int64_t(path.size()) - 1) {
                html << " ...";
            }
        }
    }
    html << "</table>";

    if(readIds.size() > maxShow) {
        html << "<p>Showing " << maxShow << " of " << readIds.size() << " paths.";
    }
}


void Assembler::exploreBidirectedAnchorGraphPath(const vector<string>& request, ostream& html)
{
    if(!bidirectedAnchors) {
        html << "<p>The Bidirected Anchor Graph is not available.";
        return;
    }

    uint64_t pathIdValue = invalid<uint64_t>;
    const bool pathIdIsPresent =
        HttpServer::getParameterValue(request, "pathId", pathIdValue);

    html <<
        "<h2>Bidirected Anchor Graph Path</h2>"
        "<p>A path is a physical read journey through BRG anchors. "
        "Each step is an oriented anchor (&gt;anchorId for strand 0, &lt;anchorId for strand 1)."
        "<form>"
        "<table>"
        "<tr><th class=left>Path id (Read id)"
        "<td class=centered>"
        "<input type=text name=pathId required style='text-align:center'";
    if(pathIdIsPresent) {
        html << " value='" << pathIdValue << "'";
    }
    html <<
        " size=8 title='Enter a read id.'>"
        "</table>"
        "<input type=submit value='Show path details'> "
        "</form>";

    if(!pathIdIsPresent) {
        return;
    }

    const ReadId readId = pathIdValue;
    if(readId >= reads->readCount()) {
        html << "<p>Path " << pathIdValue << " is invalid (read id out of range).";
        return;
    }

    const auto path = bidirectedAnchors->journey(readId);
    html << "<h1>Bidirected Anchor Graph Path " << readId << "</h1>";
    html << "<table>"
        "<tr><th class=left>Path length (nodes)<td class=centered>" << path.size() <<
        "</table>";

    if(path.empty()) {
        html << "<p>This path is empty.";
        return;
    }

    html <<
        "<table>"
        "<tr>"
        "<th>Position"
        "<th>Node"
        "<th>Orientation"
        "<th>Coverage"
        "<th>Out-degree"
        "<th>In-degree"
        "<th>Edge";

    for(uint64_t pos = 0; pos < path.size(); ++pos) {
        const auto& step = path[pos];
        const OrientedBidirectedAnchor oa(step.anchorId, step.strand);
        const uint64_t cov = bidirectedAnchors->coverage(step.anchorId);
        const uint64_t outDeg = bidirectedAnchors->getOutEdges(oa).size();
        const uint64_t inDeg = bidirectedAnchors->getInEdges(oa).size();

        html <<
            "<tr>"
            "<td class=centered>" << pos <<
            "<td class=centered><a href='exploreBidirectedAnchorGraphNode?nodeId=" << step.anchorId << "'>" <<
            step.anchorId << "</a>"
            "<td class=centered>" << (step.strand == 0 ? "&gt; (s0)" : "&lt; (s1)") <<
            "<td class=centered>" << cov <<
            "<td class=centered>" << outDeg <<
            "<td class=centered>" << inDeg;

        html << "<td class=centered>";
        if(pos + 1 < path.size()) {
            const Strand nextStrand = path[pos + 1].strand;
            if(step.strand == nextStrand) {
                html << "same-strand &rarr;";
            } else {
                html << "diff-strand &harr;";
            }
        }
    }
    html << "</table>";

    html << "<h2>Compact path</h2><p style='font-family:monospace; word-break:break-all'>";
    for(uint64_t pos = 0; pos < path.size(); ++pos) {
        if(pos > 0) {
            html << ",";
        }
        html << (path[pos].strand == 0 ? "&gt;" : "&lt;") << path[pos].anchorId;
    }
    html << "</p>";
}



void Assembler::exploreBidirectedJourney(const vector<string>& request, ostream& html)
{
    if(!bidirectedAnchors) {
        html << "<p>BRG-native anchors are not available.";
        return;
    }

    html <<
        "<h2>Physical read journey (BRG)</h2>"
        "<p>The BRG journey of a physical read is the sequence of BRG anchors it visits, "
        "with the strand (direction) at each anchor. "
        "Consecutive entries (A,sA)&rarr;(B,sB) define directed edges in the "
        "bidirected anchor graph.";

    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = HttpServer::getParameterValue(request, "readId", readId);

    // Write the form.
    html <<
        "<form>"
        "<table>"
        "<tr>"
        "<th class=left>Read id"
        "<td><input type=text name=readId" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " title='Enter a read id between 0 and " << reads->readCount() - 1 << "'>"
        "</table>"
        "<input type=submit value='Display'>"
        "</form>";

    if(!readIdIsPresent) {
        html << "Specify a numeric read id.";
        return;
    }

    if(readId >= reads->readCount()) {
        html << "<p>Invalid read id.";
        return;
    }

    const auto jrn = bidirectedAnchors->journey(readId);

    html << "<h2>BRG Journey of read " << readId << "</h2>"
        "<p>Journey length: " << jrn.size() << " anchors.";

    if(jrn.empty()) {
        html << "<p>This read has no BRG journey entries.";
        return;
    }

    // Main table.
    html <<
        "<table>"
        "<tr>"
        "<th>Position<br>in journey"
        "<th>BRG Anchor"
        "<th>Strand"
        "<th>Coverage"
        "<th>Ordinal<br>(strand 0)"
        "<th>Position"
        "<th>Edge type";

    // To show ordinals and positions, we need the intervals of each anchor
    // to find this read's entry.
    const OrientedReadId orientedReadId0(readId, 0);
    const auto readMarkers = (*markers)[orientedReadId0.getValue()];

    for(uint64_t pos = 0; pos < jrn.size(); pos++) {
        const BidirectedJourneyEntry& entry = jrn[pos];
        const BidirectedAnchorId anchorId = entry.anchorId;
        const Strand strand = entry.strand;
        const uint64_t anchorCoverage = bidirectedAnchors->coverage(anchorId);

        // Find this read's ordinal in the anchor.
        const auto intervals = (*bidirectedAnchors)[anchorId];
        uint32_t ordinal = invalid<uint32_t>;
        for(const auto& interval : intervals) {
            if(interval.readId == readId) {
                ordinal = interval.ordinal;
                break;
            }
        }

        uint32_t position = invalid<uint32_t>;
        if(ordinal != invalid<uint32_t> && ordinal < readMarkers.size()) {
            position = readMarkers[ordinal].position;
        }

        html <<
            "<tr>"
            "<td class=centered>" << pos <<
            "<td class=centered><a href='exploreBidirectedAnchor?anchorId=" << anchorId << "'>"
            << anchorId << "</a>"
            "<td class=centered>" << uint32_t(strand) <<
            "<td class=centered>" << anchorCoverage;

        if(ordinal != invalid<uint32_t>) {
            html << "<td class=centered>" << ordinal;
        } else {
            html << "<td class=centered>";
        }

        if(position != invalid<uint32_t>) {
            html << "<td class=centered>" << position;
        } else {
            html << "<td class=centered>";
        }

        // Edge type: compare strand of this entry with strand of next entry.
        html << "<td class=centered>";
        if(pos + 1 < jrn.size()) {
            const Strand nextStrand = jrn[pos + 1].strand;
            if(strand == nextStrand) {
                html << "same-strand &rarr;";
            } else {
                html << "diff-strand &harr;";
            }
        }
    }
    html << "</table>";
}
