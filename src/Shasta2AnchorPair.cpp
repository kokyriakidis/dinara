#include "Shasta2AnchorPair.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "Marker.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <iostream>

using namespace std;
using namespace dinara;

// Helper to get string representation of Shasta2AnchorId
string dinara::shasta2AnchorIdToString(Shasta2AnchorId anchorId)
{
    return std::to_string(anchorId);
}


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
        // SHASTA2_ASSERT(orientedReadId == itB->orientedReadId); // Assuming standard assert or just skip
        if (orientedReadId != itB->orientedReadId) {
             throw runtime_error("OrientedReadId mismatch in Shasta2AnchorPair constructor");
        }

        if(adjacentInJourney) {
            if(itB->positionInJourney == itA->positionInJourney + 1) {
                orientedReadIds.push_back(orientedReadId);
            }
        } else {
            if(itB->positionInJourney > itA->positionInJourney) {
                orientedReadIds.push_back(orientedReadId);
            }
        }

        ++itA;
        ++itB;
    }
}


void Shasta2AnchorPair::get(
    const Shasta2Anchors& anchors,
    vector< pair<Positions, Positions> >& positions) const
{

    const uint32_t kHalf = uint32_t(anchors.k / 2); // anchors.markers.k is not directly available, but Anchors stores k.
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
        
        // Only process is this is one of our OrientedReadIds;
        if(orientedReadId == *it) {
            ++it;

            // Access markers via anchors.
            // In Dinara, `anchors` has a reference to `markers`.
            // markers is MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>
            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t positionInJourneyA = itA->positionInJourney;
            const uint32_t positionInJourneyB = itB->positionInJourney;
            // SHASTA2_ASSERT(positionInJourneyB >= positionInJourneyA);    // Allow degenerate AnchorPair witn anchorIdA==anchorIdB
            
            const uint32_t ordinalA = itA->ordinal; 
            const uint32_t ordinalB = itB->ordinal;
            
            // Note: CompressedMarker has `position`.
            const uint32_t positionA = orientedReadMarkers[ordinalA].position + kHalf;
            const uint32_t positionB = orientedReadMarkers[ordinalB].position + kHalf;

            positions.push_back(make_pair(
                Positions(positionInJourneyA, ordinalA, positionA),
                Positions(positionInJourneyB, ordinalB, positionB)
                ));
        }

        ++itA;
        ++itB;
    }
}

void Shasta2AnchorPair::removeNegativeOffsets(const Shasta2Anchors& anchors)
{
    vector< pair<uint32_t, uint32_t> > ordinals;
    getOrdinals(anchors, ordinals);
    
    vector<OrientedReadId> newOrientedReadIds;
    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const auto& p = ordinals[i];
        if(p.second >= p.first) {
            newOrientedReadIds.push_back(orientedReadIds[i]);
        }
    }

    orientedReadIds.swap(newOrientedReadIds);
}

void Shasta2AnchorPair::get(
    const Shasta2Anchors& anchors,
    vector< pair<Positions, Positions> >& positions,
    vector< vector<Base> >& sequences) const
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);
    // Reads access is available via anchors.reads
    const Reads& reads = anchors.reads;

    get(anchors, positions);

    sequences.clear();
    sequences.resize(orientedReadIds.size());

    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const auto& positionsAB = positions[i];
        vector<Base>& sequence = sequences[i];

        const uint32_t positionA = positionsAB.first.basePosition + kHalf;
        const uint32_t positionB = positionsAB.second.basePosition + kHalf;

        // reads.getOrientedReadBase signature check needed. assuming it works.
        // In Dinara/Reads.hpp: Base getOrientedReadBase(OrientedReadId, uint32_t position) const
        for(uint32_t position=positionA; position!=positionB; position++) {
            sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
        }
    }
}

void Shasta2AnchorPair::getOrdinals(
    const Shasta2Anchors& anchors,
    vector< pair<uint32_t, uint32_t> >& ordinals) const
{
    ordinals.clear();

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

        const OrientedReadId orientedReadId = itA->orientedReadId;
        
        if(orientedReadId == *it) {
            ++it;

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;

            ordinals.push_back(make_pair(ordinalA, ordinalB));
        }

        ++itA;
        ++itB;
    }
}

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

        const OrientedReadId orientedReadId = itA->orientedReadId;
        
        if(orientedReadId == *it) {
            ++it;

            const uint32_t positionInJourneyA = itA->positionInJourney;
            const uint32_t positionInJourneyB = itB->positionInJourney;

            positionsInJourneys.push_back(make_pair(positionInJourneyA, positionInJourneyB));
        }

        ++itA;
        ++itB;
    }
}

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
    
    // In Dinara, findChildren signature:
    // void findChildren(Shasta2AnchorId, vector<Shasta2AnchorId>&, vector<uint64_t>& count, uint64_t minCoverage = 0) const;
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

        const OrientedReadId orientedReadId = itA->orientedReadId;
        
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            if(ordinalB < ordinalA) {       // Degenerate AnchorPair with AnchorIdA==AnchorIdB is ok.
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
            
            const uint32_t offset = positionB - positionA;
            sumBaseOffset += offset;
        }

        ++itA;
        ++itB;
    }

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

        const OrientedReadId orientedReadId = itA->orientedReadId;
        
        if(orientedReadId == *it) {
            ++it;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            const uint32_t ordinalA = itA->ordinal;
            const uint32_t ordinalB = itB->ordinal;
            if(ordinalB < ordinalA) {          // Degenerate AnchorPair with AnchorIdA==AnchorIdB is ok.
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
            
            const uint32_t offset = positionB - positionA;
            sumBaseOffset += offset;
            minBaseOffset = min(minBaseOffset, offset);
            maxBaseOffset = max(maxBaseOffset, offset);
        }

        ++itA;
        ++itB;
    }

    averageBaseOffset = uint32_t(std::round(double(sumBaseOffset) / double(size())));
}

uint64_t Shasta2AnchorPair::countCommon(const Shasta2AnchorPair& y) const
{
    const Shasta2AnchorPair& x = *this;

    uint64_t n = 0;
   
    auto xIt = x.orientedReadIds.begin();
    auto yIt = y.orientedReadIds.begin();
    const auto xEnd = x.orientedReadIds.end();
    const auto yEnd = y.orientedReadIds.end();

    while(xIt != xEnd && yIt != yEnd) {
        if(*xIt < *yIt) {
            ++xIt;
        } else if(*yIt < *xIt) {
            ++yIt;
        } else {
            ++n;
            ++xIt;
            ++yIt;
        }
    }
    
    return n;
}

bool Shasta2AnchorPair::contains(OrientedReadId orientedReadId) const
{
    const auto it = std::lower_bound(orientedReadIds.begin(), orientedReadIds.end(), orientedReadId);
    return it != orientedReadIds.end() and (*it == orientedReadId);
}
