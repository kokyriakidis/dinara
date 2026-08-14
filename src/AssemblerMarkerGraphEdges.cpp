#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace dinara;



// Analyze and compare the read compositions of two marker graph edges.
// This can only be done if the two edges have no duplicate OrientedReadIds
// in the markers. In that case, each OrientedReadId of an edge
// corresponds to one and only one markerInterval for each edge.
bool Assembler::analyzeMarkerGraphEdgePair(
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB,
    MarkerGraphEdgePairInfo& info
    ) const
{

    // Check for duplicate OrientedReadIds on the two edges.
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdA)) {
        return false;
    }
    if(markerGraph.edgeHasDuplicateOrientedReadIds(edgeIdB)) {
        return false;
    }

    // Prepare for the joint loop over OrientedReadIds of the two edges.
    const auto markerIntervalsA = markerGraph.edgeMarkerIntervals[edgeIdA];
    const auto markerIntervalsB = markerGraph.edgeMarkerIntervals[edgeIdB];
    const auto beginA = markerIntervalsA.begin();
    const auto beginB = markerIntervalsB.begin();
    const auto endA = markerIntervalsA.end();
    const auto endB = markerIntervalsB.end();

    // Store the total number of OrientedReadIds on the two edges.
    info.totalA = endA - beginA;
    info.totalB = endB - beginB;



    // Joint loop over the MarkerIntervals of the two edges,
    // to count the common  reads and compute average offsets.
    info.common = 0;
    int64_t sumMarkerOffsets = 0;
    int64_t sumTwiceBaseOffsets = 0;
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

        // We found a common OrientedReadId.
        ++info.common;
        const OrientedReadId orientedReadId = itA->orientedReadId;
        const auto orientedReadMarkers = (*markers)[orientedReadId.getValue()];

        // Compute the offset in markers.
        DINARA_ASSERT(itA->ordinals[1] == itA->ordinals[0] + 1);
        DINARA_ASSERT(itB->ordinals[1] == itB->ordinals[0] + 1);
        const uint32_t ordinalA = itA->ordinals[0];
        const uint32_t ordinalB = itB->ordinals[0];
        const int64_t markerOffset = int64_t(ordinalB) - int64_t(ordinalA);
        sumMarkerOffsets += markerOffset;

        // Compute the offset in bases.
        const int64_t positionA0 = int64_t(orientedReadMarkers[ordinalA].position);
        const int64_t positionA1 = int64_t(orientedReadMarkers[ordinalA+1].position);
        const int64_t positionB0 = int64_t(orientedReadMarkers[ordinalB].position);
        const int64_t positionB1 = int64_t(orientedReadMarkers[ordinalB+1].position);
        sumTwiceBaseOffsets -= positionA0;
        sumTwiceBaseOffsets -= positionA1;
        sumTwiceBaseOffsets += positionB0;
        sumTwiceBaseOffsets += positionB1;

        // Continue the joint loop.
        ++itA;
        ++itB;

    }
    info.onlyA = info.totalA - info.common;
    info.onlyB = info.totalB - info.common;

    // If there are no common reads, this is all we can do.
    if(info.common == 0) {
        info.offsetInMarkers = invalid<int64_t>;
        info.offsetInBases = invalid<int64_t>;
        info.onlyAShort = invalid<uint64_t>;
        info.onlyBShort = invalid<uint64_t>;
        return true;
    }

    // Compute the estimated offsets.
    info.offsetInMarkers = int64_t(std::round(double(sumMarkerOffsets) / double(info.common)));
    info.offsetInBases = int64_t(0.5 * std::round(double(sumTwiceBaseOffsets) / double(info.common)));

    // Now do the joint loop again, and count the onlyA and onlyB oriented reads
    // that are too short to appear in the other edge.
    itA = beginA;
    itB = beginB;
    uint64_t onlyACheck = 0;
    uint64_t onlyBCheck = 0;
    info.onlyAShort = 0;
    info.onlyBShort = 0;
    while(true) {
        if(itA == endA and itB == endB) {
            break;
        }

        else if(itB == endB or ((itA!=endA) and (itA->orientedReadId < itB->orientedReadId))) {
            // This oriented read only appears in edge A.
            ++onlyACheck;
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = (*markers)[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(getReads().getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge A in this oriented read.
            const uint32_t ordinalA0 = itA->ordinals[0];
            const uint32_t ordinalA1 = itA->ordinals[1];
            const int64_t positionA0 = int64_t(orientedReadMarkers[ordinalA0].position);
            const int64_t positionA1 = int64_t(orientedReadMarkers[ordinalA1].position);

            // Find the hypothetical positions of edge B, assuming the estimated base offset.
            const int64_t positionB0 = positionA0 + info.offsetInBases;
            const int64_t positionB1 = positionA1 + info.offsetInBases;

            // If this ends up outside the read, this counts as onlyAShort.
            if(positionB0 < 0 or positionB1 >= lengthInBases) {
                ++info.onlyAShort;
            }

            ++itA;
            continue;
        }

        else if(itA == endA or ((itB!=endB) and (itB->orientedReadId < itA->orientedReadId))) {
            // This oriented read only appears in edge B.
            ++onlyBCheck;
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const auto orientedReadMarkers = (*markers)[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(getReads().getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge B in this oriented read.
            const uint32_t ordinalB0 = itB->ordinals[0];
            const uint32_t ordinalB1 = itB->ordinals[1];
            const int64_t positionB0 = int64_t(orientedReadMarkers[ordinalB0].position);
            const int64_t positionB1 = int64_t(orientedReadMarkers[ordinalB1].position);

            // Find the hypothetical positions of edge A, assuming the estimated base offset.
            const int64_t positionA0 = positionB0 - info.offsetInBases;
            const int64_t positionA1 = positionB1 - info.offsetInBases;

            // If this ends up outside the read, this counts as onlyBShort.
            if(positionA0 < 0 or positionA1 >= lengthInBases) {
                ++info.onlyBShort;
            }

            ++itB;
            continue;
        }

        else {
            // This oriented read appears in both edges. In this loop, we
            // don't need to do anything.
            ++itA;
            ++itB;
        }
    }
    DINARA_ASSERT(onlyACheck == info.onlyA);
    DINARA_ASSERT(onlyBCheck == info.onlyB);


    return true;
}



#if 0
// More detailed analysis for a pair of marker graph edges,
// both of which must be primary.
void Assembler::analyzePrimaryMarkerGraphEdgePair(
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB) const
{
    cout << "analyzePrimaryMarkerGraphEdgePair begins for " << edgeIdA << " " << edgeIdB << endl;

    // Sanity checks.
    DINARA_ASSERT(markerGraph.edges[edgeIdA].isPrimary == 1);
    DINARA_ASSERT(markerGraph.edges[edgeIdB].isPrimary == 1);

    // The MarkerIntervals on these two edges.
    const auto markerIntervalsA = markerGraph.edgeMarkerIntervals[edgeIdA];
    const auto markerIntervalsB = markerGraph.edgeMarkerIntervals[edgeIdB];

    // Find the position of edgeA on the primary journey of each oriented read on edgeA.
    vector<uint64_t> positionInJourneyA(markerIntervalsA.size(), invalid<uint64_t>);
    for(uint64_t i=0; i<markerIntervalsA.size(); i++) {
        const OrientedReadId orientedReadId = markerIntervalsA[i].orientedReadId;
        const auto journey = markerGraph.primaryJourneys[orientedReadId.getValue()];
        for(uint64_t position=0; position<journey.size(); position++) {
            if(journey[position].edgeId == edgeIdA) {
                positionInJourneyA[i] = position;
                break;
            }
        }
        DINARA_ASSERT(positionInJourneyA[i] != invalid<uint64_t>);
    }

    // Find the position of edgeB on the primary journey of each oriented read on edgeB.
    vector<uint64_t> positionInJourneyB(markerIntervalsB.size(), invalid<uint64_t>);
    for(uint64_t i=0; i<markerIntervalsB.size(); i++) {
        const OrientedReadId orientedReadId = markerIntervalsB[i].orientedReadId;
        const auto journey = markerGraph.primaryJourneys[orientedReadId.getValue()];
        for(uint64_t position=0; position<journey.size(); position++) {
            if(journey[position].edgeId == edgeIdB) {
                positionInJourneyB[i] = position;
                break;
            }
        }
        DINARA_ASSERT(positionInJourneyB[i] != invalid<uint64_t>);
    }


    // The MarkerGraphEdgeIds that we encountered so far by moving forward from edgeA on
    // the primary journeys of oriented reads on edgeA.
    std::set<MarkerGraphEdgeId> edgeIdsForwardA;

    // The MarkerGraphEdgeIds that we encountered so far by moving backward from edgeB on
    // the primary journeys of oriented reads on edgeB.
    std::set<MarkerGraphEdgeId> edgeIdsBackwardB;

    // Iterate over offsets in the primary journeys.
    // For journeys of the oriented reads on edgeA, we use positive offsets.
    // For journeys of the oriented reads on edgeB, we use negative offsets.
    for(uint64_t offset=1; ; ++offset) {

        uint64_t activeCountA = 0;
        for(uint64_t i=0; i<markerIntervalsA.size(); i++) {
            const OrientedReadId orientedReadId = markerIntervalsA[i].orientedReadId;
            const auto journey = markerGraph.primaryJourneys[orientedReadId.getValue()];
            const uint64_t position = positionInJourneyA[i] + offset;
            if(position >= journey.size()) {
                continue;
            }
            ++activeCountA;
            const MarkerGraphEdgeId edgeId = journey[position].edgeId;

            if(not edgeIdsForwardA.contains(edgeId)) {
                edgeIdsForwardA.insert(edgeId);

                if(edgeIdsBackwardB.contains(edgeId)) {
                    MarkerGraphEdgePairInfo infoA;
                    analyzeMarkerGraphEdgePair(edgeIdA, edgeId, infoA);
                    MarkerGraphEdgePairInfo infoB;
                    analyzeMarkerGraphEdgePair(edgeId, edgeIdB, infoB);
                    cout << "At offset " << offset << " found " << edgeId <<
                        ", common " << infoA.common << " " << infoB.common << ", total offset " <<
                        infoA.offsetInBases+ infoB.offsetInBases << endl;
                }
            }
        }

        uint64_t activeCountB = 0;
        for(uint64_t i=0; i<markerIntervalsB.size(); i++) {
            const OrientedReadId orientedReadId = markerIntervalsB[i].orientedReadId;
            const auto journey = markerGraph.primaryJourneys[orientedReadId.getValue()];
            if(offset > positionInJourneyB[i]) {
                continue;
            }
            const uint64_t position = positionInJourneyB[i] - offset;
            ++activeCountB;
            const MarkerGraphEdgeId edgeId = journey[position].edgeId;

            if(not edgeIdsBackwardB.contains(edgeId)) {
                edgeIdsBackwardB.insert(edgeId);

                if(edgeIdsForwardA.contains(edgeId)) {
                    MarkerGraphEdgePairInfo infoA;
                    analyzeMarkerGraphEdgePair(edgeIdA, edgeId, infoA);
                    MarkerGraphEdgePairInfo infoB;
                    analyzeMarkerGraphEdgePair(edgeId, edgeIdB, infoB);
                    cout << "At offset " << offset << " found " << edgeId <<
                        ", common " << infoA.common << " " << infoB.common << endl;
                }
            }
        }

        if(activeCountA == 0 or activeCountB == 0) {
            break;
        }
    }
}
#endif



// Estimate the offset, in bases, between two marker graph edges.
// This assumes, WITHOUT CHECKING, that each of the two edges has no duplicate
// oriented reads. This assumption is satisfied for primary marker graph edges
// in Mode 3 assembly.
// If there are common oriented reads between the two edges, this uses
// analyzeMarkerGraphEdgePair.
// This can fail, in which case it returns invalid<uint64_t>.
uint64_t Assembler::estimateBaseOffsetUnsafe(
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB) const
{
    // If there are common oriented reads between the two edges, use
    // analyzeMarkerGraphEdgePair. This is the most common case.
    if(countCommonOrientedReadsUnsafe(edgeIdA, edgeIdB) > 0) {
        MarkerGraphEdgePairInfo info;
        DINARA_ASSERT(analyzeMarkerGraphEdgePair(edgeIdA, edgeIdB, info));
        if(info.offsetInBases >= 0) {
            return info.offsetInBases;
        } else {
            return invalid<uint64_t>;
        }
    } else {
        return invalid<uint64_t>;
    }

#if 0
    // There are no common oriented reads between the two edges.
    // Find a primary marker graph edge in-between that has common
    // oriented reads with both edgeIdA and edgeIdB.

    // Sanity checks.
    DINARA_ASSERT(markerGraph.edges[edgeIdA].isPrimary == 1);
    DINARA_ASSERT(markerGraph.edges[edgeIdB].isPrimary == 1);

    // The MarkerIntervals on these two edges.
    const auto markerIntervalsA = markerGraph.edgeMarkerIntervals[edgeIdA];
    const auto markerIntervalsB = markerGraph.edgeMarkerIntervals[edgeIdB];

    // Find the position of edgeA on the primary journey of each oriented read on edgeA.
    vector<uint64_t> positionInJourneyA(markerIntervalsA.size(), invalid<uint64_t>);
    for(uint64_t i=0; i<markerIntervalsA.size(); i++) {
        const OrientedReadId orientedReadId = markerIntervalsA[i].orientedReadId;
        const auto journey = markerGraph.primaryJourneys[orientedReadId.getValue()];
        for(uint64_t position=0; position<journey.size(); position++) {
            if(journey[position].edgeId == edgeIdA) {
                positionInJourneyA[i] = position;
                break;
            }
        }
        DINARA_ASSERT(positionInJourneyA[i] != invalid<uint64_t>);
    }

    // Find the position of edgeB on the primary journey of each oriented read on edgeB.
    vector<uint64_t> positionInJourneyB(markerIntervalsB.size(), invalid<uint64_t>);
    for(uint64_t i=0; i<markerIntervalsB.size(); i++) {
        const OrientedReadId orientedReadId = markerIntervalsB[i].orientedReadId;
        const auto journey = markerGraph.primaryJourneys[orientedReadId.getValue()];
        for(uint64_t position=0; position<journey.size(); position++) {
            if(journey[position].edgeId == edgeIdB) {
                positionInJourneyB[i] = position;
                break;
            }
        }
        DINARA_ASSERT(positionInJourneyB[i] != invalid<uint64_t>);
    }


    // The MarkerGraphEdgeIds that we encountered so far by moving forward from edgeA on
    // the primary journeys of oriented reads on edgeA.
    std::set<MarkerGraphEdgeId> edgeIdsForwardA;

    // The MarkerGraphEdgeIds that we encountered so far by moving backward from edgeB on
    // the primary journeys of oriented reads on edgeB.
    std::set<MarkerGraphEdgeId> edgeIdsBackwardB;

    // The best edgeId we found, and the lowest of its common oriented reads
    // with edgeIdA and edgeIdB.
    uint64_t edgeIdBest = invalid<uint64_t>;
    uint64_t commonBest = 0;

    // Iterate over offsets in the primary journeys.
    // For journeys of the oriented reads on edgeA, we use positive offsets.
    // For journeys of the oriented reads on edgeB, we use negative offsets.
    for(uint64_t offset=1; ; ++offset) {

        uint64_t activeCountA = 0;
        for(uint64_t i=0; i<markerIntervalsA.size(); i++) {
            const OrientedReadId orientedReadId = markerIntervalsA[i].orientedReadId;
            const auto journey = markerGraph.primaryJourneys[orientedReadId.getValue()];
            const uint64_t position = positionInJourneyA[i] + offset;
            if(position >= journey.size()) {
                continue;
            }
            ++activeCountA;
            const MarkerGraphEdgeId edgeId = journey[position].edgeId;

            if(not edgeIdsForwardA.contains(edgeId)) {
                edgeIdsForwardA.insert(edgeId);

                if(edgeIdsBackwardB.contains(edgeId)) {
                    const uint64_t commonCountA = countCommonOrientedReadsUnsafe(edgeIdA, edgeId);
                    const uint64_t commonCountB = countCommonOrientedReadsUnsafe(edgeId, edgeIdB);
                    const uint64_t commonCountMin = min(commonCountA, commonCountB);
                    if(commonCountMin > commonBest) {
                        edgeIdBest = edgeId;
                        commonBest = commonCountMin;
                    }
                }
            }
        }

        uint64_t activeCountB = 0;
        for(uint64_t i=0; i<markerIntervalsB.size(); i++) {
            const OrientedReadId orientedReadId = markerIntervalsB[i].orientedReadId;
            const auto journey = markerGraph.primaryJourneys[orientedReadId.getValue()];
            if(offset > positionInJourneyB[i]) {
                continue;
            }
            const uint64_t position = positionInJourneyB[i] - offset;
            ++activeCountB;
            const MarkerGraphEdgeId edgeId = journey[position].edgeId;

            if(not edgeIdsBackwardB.contains(edgeId)) {
                edgeIdsBackwardB.insert(edgeId);

                if(edgeIdsForwardA.contains(edgeId)) {
                    const uint64_t commonCountA = countCommonOrientedReadsUnsafe(edgeIdA, edgeId);
                    const uint64_t commonCountB = countCommonOrientedReadsUnsafe(edgeId, edgeIdB);
                    const uint64_t commonCountMin = min(commonCountA, commonCountB);
                    if(commonCountMin > commonBest) {
                        edgeIdBest = edgeId;
                        commonBest = commonCountMin;
                    }
                }
            }
        }

        if(activeCountA == 0 or activeCountB == 0) {
            break;
        }
    }

    if(commonBest == 0) {
        return invalid<uint64_t>;
    }

    // edgeIdBest has common oriented reads with both edgeIdA and edgeIdB.
    MarkerGraphEdgePairInfo infoA;
    MarkerGraphEdgePairInfo infoB;
    DINARA_ASSERT(analyzeMarkerGraphEdgePair(edgeIdA, edgeIdBest, infoA));
    DINARA_ASSERT(analyzeMarkerGraphEdgePair(edgeIdBest, edgeIdB, infoB));
    DINARA_ASSERT(infoA.common > 0);
    DINARA_ASSERT(infoB.common > 0);
    return infoA.offsetInBases + infoB.offsetInBases;
#endif
}



// Count the number of common oriented reads between two marker graph edges.
// This assumes, WITHOUT CHECKING, that each of the two edges has no duplicate
// oriented reads. This assumption is satisfied for primary marker graph edges
// in Mode 3 assembly.
uint64_t Assembler::countCommonOrientedReadsUnsafe(
    MarkerGraphEdgeId edgeIdA,
    MarkerGraphEdgeId edgeIdB) const
{
    // Prepare for the joint loop over OrientedReadIds of the two edges.
    const auto markerIntervalsA = markerGraph.edgeMarkerIntervals[edgeIdA];
    const auto markerIntervalsB = markerGraph.edgeMarkerIntervals[edgeIdB];
    const auto beginA = markerIntervalsA.begin();
    const auto beginB = markerIntervalsB.begin();
    const auto endA = markerIntervalsA.end();
    const auto endB = markerIntervalsB.end();


    // Joint loop over the MarkerIntervals of the two edges,
    // to count the common  reads and compute average offsets.
    // This assumes that there are no duplicate oriented reads
    // on the two edges.
    uint64_t n = 0;
    auto itA = beginA;
    auto itB = beginB;
    while(itA != endA and itB != endB) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
        } else if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        } else {
            // We found a common OrientedReadId.
            ++n;
            ++itA;
            ++itB;
        }

    }
    return n;
}
