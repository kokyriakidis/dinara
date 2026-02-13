// HTTP server handlers for BRG-native anchor visualization.

#include "Assembler.hpp"
#include "mode3-BidirectedAnchor.hpp"
#include "HttpServer.hpp"
#include "invalid.hpp"
#include "Marker.hpp"
#include "Reads.hpp"

using namespace dinara;
using namespace mode3;
using namespace std;



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
