// Shasta.
#include "Shasta2AnchorPair.hpp"
#include "Shasta2Anchors.hpp"
#include "HttpServer.hpp"
#include "Shasta2Journeys.hpp"
#include "Reads.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/iterator/function_output_iterator.hpp>

// Standard library.
#include <cmath>
#include <iomanip>



Shasta2AnchorPair::Shasta2AnchorPair(
    const Shasta2Anchors& anchors,
    Shasta2AnchorId anchorIdA,
    Shasta2AnchorId anchorIdB,
    bool adjacentInJourney) :
    anchorIdA(anchorIdA),
    anchorIdB(anchorIdB)
{
    const Shasta2Anchor anchorA = anchors[anchorIdA];
    const Shasta2Anchor anchorB = anchors[anchorIdB];

    // Loop over common oriented reads between these two anchors.
    // If adjacentInJourney is false, the journey offset is required to be positive.
    // If adjacentInJourney is true, the journey offset is required to be exactly 1.

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    while(itA != endA and itB != endB) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        const OrientedReadId orientedReadId = itA->orientedReadId;
        DINARA_ASSERT(orientedReadId == itB->orientedReadId);

        if(adjacentInJourney) {
            if(itB->positionInJourney == itA->positionInJourney + 1 &&
               itB->position >= itA->position) {
                orientedReadIds.push_back(orientedReadId);
            }
        } else {
            if(itB->positionInJourney > itA->positionInJourney &&
               itB->position >= itA->position) {
                orientedReadIds.push_back(orientedReadId);
            }
        }

        ++itA;
        ++itB;
    }
}



// Get positions in journey and base positions
// for each of the two reads and for each of the two anchors.
// The positions returned are the midpoint of the markers
// corresponding to anchorIdA and anchorIdB.
void Shasta2AnchorPair::get(
    const Shasta2Anchors& anchors,
    vector< pair<Positions, Positions> >& positions) const
{

    const uint32_t kHalf = uint32_t(anchors.k / 2);
    positions.clear();

    const Shasta2Anchor anchorA = anchors[anchorIdA];
    const Shasta2Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        DINARA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const uint32_t positionInJourneyA = itA->positionInJourney;
            const uint32_t positionInJourneyB = itB->positionInJourney;
            DINARA_ASSERT(positionInJourneyB >= positionInJourneyA);    // Allow degenerate Shasta2AnchorPair witn anchorIdA==anchorIdB

            positions.push_back(make_pair(
                Positions(positionInJourneyA, itA->position),
                Positions(positionInJourneyB, itB->position)
                ));
        }

        ++itA;
        ++itB;
    }
}



// Remove from the Shasta2AnchorPair OrientedReadIds that have negative offsets.
void Shasta2AnchorPair::removeNegativeOffsets(const Shasta2Anchors& anchors)
{
    vector< pair<uint32_t, uint32_t> > anchorPositions;
    getAnchorPositions(anchors, anchorPositions);
    DINARA_ASSERT(anchorPositions.size() == orientedReadIds.size());

    // Also get journey positions for diagnostic output.
    vector< pair<uint32_t, uint32_t> > journeyPositions;
    getPositionsInJourneys(anchors, journeyPositions);

    vector<OrientedReadId> newOrientedReadIds;
    uint64_t removedCount = 0;
    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const auto& p = anchorPositions[i];
        if(p.second >= p.first) {
            newOrientedReadIds.push_back(orientedReadIds[i]);
        } else {
            if(removedCount == 0) {
                cout << "Negative base offsets on edge "
                     << anchorIdA << " -> " << anchorIdB << ":" << endl;
            }
            if(removedCount < 5) {
                cout << "  " << orientedReadIds[i]
                     << " basePos A=" << p.first << " B=" << p.second
                     << " (offset " << int64_t(p.second) - int64_t(p.first) << ")"
                     << " journeyPos A=" << journeyPositions[i].first
                     << " B=" << journeyPositions[i].second
                     << " ordinal A=" << anchors[anchorIdA][0].ordinal
                     << " B=" << anchors[anchorIdB][0].ordinal
                     << endl;
            }
            ++removedCount;
        }
    }
    if(removedCount > 0) {
        cout << "  Removed " << removedCount << " / " << orientedReadIds.size()
             << " reads with negative base offset from edge "
             << anchorIdA << " -> " << anchorIdB << endl;
    }

    orientedReadIds.swap(newOrientedReadIds);
}




void Shasta2AnchorPair::assertNoNegativeOffsets(const Shasta2Anchors& anchors) const
{
    vector< pair<uint32_t, uint32_t> > anchorPositions;
    getAnchorPositions(anchors, anchorPositions);
    DINARA_ASSERT(anchorPositions.size() == orientedReadIds.size());

    for(uint64_t i = 0; i < orientedReadIds.size(); i++) {
        const auto& p = anchorPositions[i];
        if(p.second < p.first) {
            cout << "Negative base offset on edge "
                 << anchorIdA << " -> " << anchorIdB
                 << ": " << orientedReadIds[i]
                 << " basePos A=" << p.first << " B=" << p.second
                 << " (offset " << int64_t(p.second) - int64_t(p.first) << ")"
                 << endl;
            DINARA_ASSERT(false);
        }
    }
}

bool Shasta2AnchorPair::hasNegativeOffsets(const Shasta2Anchors& anchors) const
{
    vector< pair<uint32_t, uint32_t> > anchorPositions;
    getAnchorPositions(anchors, anchorPositions);
    for(const auto& p : anchorPositions) {
        if(p.second < p.first) return true;
    }
    return false;
}



// Same as the above, but also returns compute the sequences.
void Shasta2AnchorPair::get(
    const Shasta2Anchors& anchors,
    vector< pair<Positions, Positions> >& positions,
    vector< vector<Base> >& sequences) const
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);
    const Reads& reads = anchors.reads;

    get(anchors, positions);

    sequences.clear();
    sequences.resize(orientedReadIds.size());

    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const auto& positionsAB = positions[i];
        vector<Base>& sequence = sequences[i];

        const uint32_t positionA = positionsAB.first.basePosition;
        const uint32_t positionB = positionsAB.second.basePosition;

        for(uint32_t position=positionA; position!=positionB; position++) {
            sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
        }
    }
}



// Same as the above, but only compute the anchor positions (midpoints).
void Shasta2AnchorPair::getAnchorPositions(
    const Shasta2Anchors& anchors,
    vector< pair<uint32_t, uint32_t> >& anchorPositions) const
{
    anchorPositions.clear();

    const Shasta2Anchor anchorA = anchors[anchorIdA];
    const Shasta2Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        DINARA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            anchorPositions.push_back(make_pair(itA->position, itB->position));
        }

        ++itA;
        ++itB;
    }

    DINARA_ASSERT(it == orientedReadIds.end());
}



// Just return the journey positions.
void Shasta2AnchorPair::getPositionsInJourneys(
    const Shasta2Anchors& anchors,
    vector< pair<uint32_t, uint32_t> >& positionsInJourneys) const
{
    positionsInJourneys.clear();

    const Shasta2Anchor anchorA = anchors[anchorIdA];
    const Shasta2Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        DINARA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const uint32_t positionInJourneyA = itA->positionInJourney;
            const uint32_t positionInJourneyB = itB->positionInJourney;

            positionsInJourneys.push_back(make_pair(positionInJourneyA, positionInJourneyB));
        }

        ++itA;
        ++itB;
    }

    DINARA_ASSERT(it == orientedReadIds.end());

}



// This finds AnchorPairs as follows:
// - anchorIdA is as specified.
// - Coverage is at least minCoverage.
// - All oriented reads have a journey offset equal to 1.
void Shasta2AnchorPair::createChildren(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    Shasta2AnchorId anchorIdA,
    uint64_t minCoverage,
    vector<Shasta2AnchorPair>& anchorPairs
    )
{
    // Find possible choices for anchorIdB.
    vector<Shasta2AnchorId> anchorIdsB;
    vector<uint64_t> coverage;
    anchors.findChildren(journeys, anchorIdA, anchorIdsB, coverage, minCoverage);

    anchorPairs.clear();
    for(const Shasta2AnchorId anchorIdB: anchorIdsB) {
        anchorPairs.emplace_back(anchors, anchorIdA, anchorIdB, true);
    }
}



uint32_t Shasta2AnchorPair::getAverageOffset(const Shasta2Anchors& anchors) const
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    uint64_t sumBaseOffset = 0;

    const Shasta2Anchor anchorA = anchors[anchorIdA];
    const Shasta2Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        DINARA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            if(ordinalB < ordinalA) {       // Degenerate Shasta2AnchorPair with AnchorIdA==AnchorIdB is ok.
                throw runtime_error(
                    "Order violation at anchor pair " +
                    shasta2AnchorIdToString(anchorIdA) + " " +
                    shasta2AnchorIdToString(anchorIdB) + " " +
                    orientedReadId.getString() + " ordinals " +
                    to_string(ordinalA) + " " +
                    to_string(ordinalB));
            }
            const uint32_t positionA = orientedReadMarkers[ordinalA].position + kHalf;
            const uint32_t positionB = orientedReadMarkers[ordinalB].position + kHalf;
            DINARA_ASSERT(positionB >= positionA);      // Degenerate Shasta2AnchorPair with AnchorIdA==AnchorIdB is ok.

            const uint32_t offset = positionB - positionA;
            sumBaseOffset += offset;
        }

        ++itA;
        ++itB;
    }

    DINARA_ASSERT(it == orientedReadIds.end());

    return uint32_t(std::round(double(sumBaseOffset) / double(size())));
}



void Shasta2AnchorPair::getOffsets(
    const Shasta2Anchors& anchors,
    uint32_t& averageBaseOffset,
    uint32_t& minBaseOffset,
    uint32_t& maxBaseOffset) const
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    uint64_t sumBaseOffset = 0;
    minBaseOffset = std::numeric_limits<uint32_t>::max();
    maxBaseOffset = 0;

    const Shasta2Anchor anchorA = anchors[anchorIdA];
    const Shasta2Anchor anchorB = anchors[anchorIdB];

    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    auto itA = beginA;
    auto itB = beginB;
    auto it = orientedReadIds.begin();
    const auto itEnd = orientedReadIds.end();
    while(itA != endA and itB != endB and it != itEnd) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        const OrientedReadId orientedReadId = itA->orientedReadId;
        DINARA_ASSERT(orientedReadId == itB->orientedReadId);

        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            if(ordinalB < ordinalA) {          // Degenerate Shasta2AnchorPair with AnchorIdA==AnchorIdB is ok.
                throw runtime_error(
                    "Order violation at anchor pair " +
                    shasta2AnchorIdToString(anchorIdA) + " " +
                    shasta2AnchorIdToString(anchorIdB) + " " +
                    orientedReadId.getString() + " ordinals " +
                    to_string(ordinalA) + " " +
                    to_string(ordinalB));
            }
            const uint32_t positionA = orientedReadMarkers[ordinalA].position + kHalf;
            const uint32_t positionB = orientedReadMarkers[ordinalB].position + kHalf;
            DINARA_ASSERT(positionB > positionA);

            const uint32_t offset = positionB - positionA;
            sumBaseOffset += offset;
            minBaseOffset = min(minBaseOffset, offset);
            maxBaseOffset = max(maxBaseOffset, offset);
        }

        ++itA;
        ++itB;
    }

    DINARA_ASSERT(it == orientedReadIds.end());

    averageBaseOffset = uint32_t(std::round(double(sumBaseOffset) / double(size())));
}



// Count OrientedReadIds in common with another Shasta2AnchorPair.
uint64_t Shasta2AnchorPair::countCommon(const Shasta2AnchorPair& y) const
{
    const Shasta2AnchorPair& x = *this;

    uint64_t n = 0;
    auto counter = [&n](auto){++n;};

    std::set_intersection(
        x.orientedReadIds.begin(),
        x.orientedReadIds.end(),
        y.orientedReadIds.begin(),
        y.orientedReadIds.end(),
        boost::make_function_output_iterator(counter)
    );

    return n;

}



bool Shasta2AnchorPair::contains(OrientedReadId orientedReadId) const
{
    const auto it = std::lower_bound(orientedReadIds.begin(), orientedReadIds.end(), orientedReadId);
    return it != orientedReadIds.end() and (*it == orientedReadId);
}



// Return the url for the exploreShasta2AnchorPair2 page for this Shasta2AnchorPair.
string Shasta2AnchorPair::url() const
{
    string s =
        "exploreShasta2AnchorPair2?"
        "anchorIdAString=" + HttpServer::urlEncode(shasta2AnchorIdToString(anchorIdA)) +
        "&anchorIdBString=" + HttpServer::urlEncode(shasta2AnchorIdToString(anchorIdB)) +
        "&orientedReadIdsString=";

    for(const OrientedReadId orientedReadId: orientedReadIds) {
        s += orientedReadId.getString();
        s += ",";
    }

    return s;
}



void Shasta2AnchorPair::writeAllHtml(
    ostream& html,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys) const
{
    // Write the summary and oriented reads.
    writeSummaryHtml(html, anchors);
    writeOrientedReadIdsHtml(html, anchors);

    // Get the positions of anchorIdA and anchorIdB in the
    // journeys of all OrientedReadids.
    vector< pair<uint32_t, uint32_t> > positionsInJourneys;
    getPositionsInJourneys(anchors, positionsInJourneys);

    // Write the journey portions between anchorIdA and anchorIdB.
    writeJourneysHtml(html, journeys, positionsInJourneys);

}



void Shasta2AnchorPair::writeSummaryHtml(ostream& html, const Shasta2Anchors& anchors) const
{
    Shasta2AnchorPairInfo info;
    anchors.analyzeAnchorPair(anchorIdA, anchorIdB, info);

    html <<
        "<table>"
        "<tr><th>Shasta2Anchor A<td class=centered>" << shasta2AnchorIdToString(anchorIdA) <<
        "<tr><th>Shasta2Anchor B<td class=centered>" << shasta2AnchorIdToString(anchorIdB) <<
        "<tr><th>Coverage<td class=centered>" << size() <<
        "<tr><th>Total A<td class=centered>" << info.totalA <<
        "<tr><th>Total B<td class=centered>" << info.totalB <<
        "<tr><th>Common<td class=centered>" << info.common <<
        "<tr><th>Only A<td class=centered>" << info.onlyA <<
        "<tr><th>Only B<td class=centered>" << info.onlyB;

    if(info.common != 0) {
        html <<
            "<tr><th>Only A, short<td class=centered>" << info.onlyAShort <<
            "<tr><th>Only B, short<td class=centered>" << info.onlyBShort <<
            "<tr><th>Only A, missing<td class=centered>" << (info.onlyA - info.onlyAShort) <<
            "<tr><th>Only B, missing<td class=centered>" << (info.onlyB - info.onlyBShort) <<
            "<tr><th>Jaccard<td class=centered>" << std::fixed << std::setprecision(2) << info.jaccard() <<
            "<tr><th>Corrected Jaccard<td class=centered>" << std::fixed << std::setprecision(2) << info.correctedJaccard() <<
            "<tr><th>Offset in markers<td class=centered>" << info.offsetInMarkers <<
            "<tr><th>Offset in bases<td class=centered>" << info.offsetInBases;
    }

    html <<
        "</table>";
}



void Shasta2AnchorPair::writeOrientedReadIdsHtml(ostream& html, const Shasta2Anchors& anchors) const
{
    Shasta2AnchorPairInfo info;
    anchors.analyzeAnchorPair(anchorIdA, anchorIdB, info);

    const Shasta2Anchor anchorA = anchors[anchorIdA];
    const Shasta2Anchor anchorB = anchors[anchorIdB];
    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    html << "<h3>Oriented reads</h3>";
    html <<
        "<p>In the following table, positions in red are hypothetical, based on the estimated base offset."
        "<table>"
        "<tr>"
        "<th rowspan=2>Oriented<br>read id"
        "<th colspan=2>Length"
        "<th colspan=3>Anchor A"
        "<th colspan=3>Anchor B"
        "<th rowspan=2>Journey<br>offset"
        "<th rowspan=2>Ordinal<br>offset"
        "<th rowspan=2>Base<br>offset"
        "<th rowspan=2>Classification"
        "<tr>"
        "<th>Markers"
        "<th>Bases"
        "<th>Journey<br>position"
        "<th>Ordinal"
        "<th>Middle<br>position"
        "<th>Journey<br>position"
        "<th>Ordinal"
        "<th>Middle<br>position";

    auto itA = beginA;
    auto itB = beginB;
    while(true) {
        if(itA == endA and itB == endB) {
            break;
        }

        if(itB == endB or ((itA != endA) and (itA->orientedReadId < itB->orientedReadId))) {
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(anchors.reads.getReadRawSequenceLength(orientedReadId.getReadId()));
            const uint32_t positionInJourneyA = itA->positionInJourney;
            const uint32_t ordinalA = itA->ordinal;
            const int64_t positionA = int64_t(orientedReadMarkers[ordinalA].position) + int64_t(anchors.k / 2);

            const int64_t positionB = positionA + info.offsetInBases;
            const bool isShort = positionB < 0 or positionB >= lengthInBases;

            html <<
                "<tr><td class=centered><a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>" <<
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << positionInJourneyA <<
                "<td class=centered>" << ordinalA <<
                "<td class=centered>" << positionA <<
                "<td><td><td class=centered style='color:Red'>" << positionB <<
                "<td><td><td>" <<
                "<td class=centered>OnlyA, " << (isShort ? "short" : "missing");

            ++itA;
            continue;
        }

        if(itA == endA or ((itB != endB) and (itB->orientedReadId < itA->orientedReadId))) {
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(anchors.reads.getReadRawSequenceLength(orientedReadId.getReadId()));
            const uint32_t positionInJourneyB = itB->positionInJourney;
            const uint32_t ordinalB = itB->ordinal;
            const int64_t positionB = int64_t(orientedReadMarkers[ordinalB].position) + int64_t(anchors.k / 2);

            const int64_t positionA = positionB - info.offsetInBases;
            const bool isShort = positionA < 0 or positionA >= lengthInBases;

            html <<
                "<tr><td class=centered><a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>" <<
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << lengthInBases <<
                "<td><td><td class=centered style='color:Red'>" << positionA <<
                "<td class=centered>" << positionInJourneyB <<
                "<td class=centered>" << ordinalB <<
                "<td class=centered>" << positionB <<
                "<td><td><td>" <<
                "<td class=centered>OnlyB, " << (isShort ? "short" : "missing");

            ++itB;
            continue;
        }

        const OrientedReadId orientedReadId = itA->orientedReadId;
        DINARA_ASSERT(orientedReadId == itB->orientedReadId);
        const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
        const int64_t lengthInBases = int64_t(anchors.reads.getReadRawSequenceLength(orientedReadId.getReadId()));
        const uint32_t positionInJourneyA = itA->positionInJourney;
        const uint32_t positionInJourneyB = itB->positionInJourney;
        const uint32_t ordinalA = itA->ordinal;
        const uint32_t ordinalB = itB->ordinal;
        const int64_t positionA = int64_t(orientedReadMarkers[ordinalA].position) + int64_t(anchors.k / 2);
        const int64_t positionB = int64_t(orientedReadMarkers[ordinalB].position) + int64_t(anchors.k / 2);

        html <<
            "<tr><td class=centered><a href='exploreRead?readId=" << orientedReadId.getReadId() <<
            "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>" <<
            "<td class=centered>" << orientedReadMarkers.size() <<
            "<td class=centered>" << lengthInBases <<
            "<td class=centered>" << positionInJourneyA <<
            "<td class=centered>" << ordinalA <<
            "<td class=centered>" << positionA <<
            "<td class=centered>" << positionInJourneyB <<
            "<td class=centered>" << ordinalB <<
            "<td class=centered>" << positionB <<
            "<td class=centered>" << (positionInJourneyB - positionInJourneyA) <<
            "<td class=centered>" << (ordinalB - ordinalA) <<
            "<td class=centered>" << (positionB - positionA) <<
            "<td class=centered>Common";

        ++itA;
        ++itB;
    }
    html << "</table>";

}



void Shasta2AnchorPair::writeJourneysHtml(
    ostream& html,
    const Shasta2Journeys& journeys,
    const vector< pair<uint32_t, uint32_t> >& positionsInJourneys   // As computed by getPositionsInJourneys.
    ) const
{
    html << "<h3>Shasta2Journey portions within this anchor pair</h3>";
    html << "<table>";

    for(uint64_t i=0; i<size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const Shasta2Journey& journey = journeys[orientedReadId];

        const auto& positionInJourney = positionsInJourneys[i];
        const auto positionInJourneyA = positionInJourney.first;
        const auto positionInJourneyB = positionInJourney.second;

        if(html) {
            html << "<tr><th class=centered>" << orientedReadId;
        }
        for(auto position=positionInJourneyA+1; position<positionInJourneyB; position++) {
            const Shasta2AnchorId anchorId = journey[position];
            html << "<td class=centered>" << shasta2AnchorIdToString(anchorId);
        }
    }

    html << "</table>";

}
