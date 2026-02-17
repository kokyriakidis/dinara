#pragma once

#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "Reads.hpp"

#include <vector>
#include <utility>
#include <iostream>

namespace dinara {
        class Shasta2AnchorPair;
}

class dinara::Shasta2AnchorPair {
public:

    // The constructor creates an AnchorPair between anchorIdA and anchorIdB.
    // - If adjacentInJourney is false, it includes all OrientedReadIds
    //  that visit anchorIdA and then, later in their Journey, anchorIdB.
    // - If adjacentInJourney is true, only OrientedReadIds that visit anchorIdB
    //   immediately after visiting anchorIdA are included.
    //   That is, the journey offset between anchorIdA and anchorIdB
    //   for OrientedReadIds that are included must be exactly 1.
    Shasta2AnchorPair(
        const Shasta2Anchors&,
        Shasta2AnchorId anchorIdA,
        Shasta2AnchorId anchorIdB,
        bool adjacentInJourney);

    Shasta2AnchorPair() {}

    Shasta2AnchorPair(const Shasta2AnchorPair& that) :
        anchorIdA(that.anchorIdA),
        anchorIdB(that.anchorIdB),
        orientedReadIds(that.orientedReadIds)
    {}
    Shasta2AnchorPair& operator=(const Shasta2AnchorPair&) = default;

    // This finds AnchorPairs as follows:
    // - anchorIdA is as specified.
    // - Coverage is at least minCoverage.
    // - All oriented reads have a journey offset equal to 1.
    static void createChildren(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        Shasta2AnchorId anchorIdA,
        uint64_t minCoverage,
        vector<Shasta2AnchorPair>&
        );

    Shasta2AnchorId anchorIdA = invalid<Shasta2AnchorId>;
    Shasta2AnchorId anchorIdB = invalid<Shasta2AnchorId>;

    vector<OrientedReadId> orientedReadIds;
    uint64_t size() const
    {
        return orientedReadIds.size();
    }

    uint32_t getAverageOffset(const Shasta2Anchors&) const;
    void getOffsets(
        const Shasta2Anchors&,
        uint32_t& averageBaseOffset,
        uint32_t& minBaseOffset,
        uint32_t& maxBaseOffset) const;

    // Get positions in journey, ordinals, base positions
    // for each of the two reads and for each of the two anchors.
    // The positions returned are the midpoint of the markers
    // corresponding to anchorIdA and anchorIdB.
    class Positions {
    public:
        uint32_t positionInJourney;
        uint32_t ordinal;
        uint32_t basePosition;
        Positions(
            uint32_t positionInJourney,
            uint32_t ordinal,
            uint32_t basePosition) :
            positionInJourney(positionInJourney),
            ordinal(ordinal),
            basePosition(basePosition)
        {}
        Positions() {}

        template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
        {
            ar & positionInJourney;
            ar & ordinal;
            ar & basePosition;
        }
    };
    void get(
        const Shasta2Anchors&,
        vector< pair<Positions, Positions> >& positions) const;

    // Same as the above, but also returns compute the sequences.
    void get(
        const Shasta2Anchors&,
        vector< pair<Positions, Positions> >& positions,
        vector< vector<Base> >&) const;

    // Just return the ordinals.
    void getOrdinals(const Shasta2Anchors&, vector< pair<uint32_t, uint32_t> >&) const;

    // Just return the positions in journeys.
    void getPositionsInJourneys(const Shasta2Anchors&, vector< pair<uint32_t, uint32_t> >&) const;

    // Count OrientedReadIds in common with another AnchorPair.
    uint64_t countCommon(const Shasta2AnchorPair&) const;

    // Remove from the AnchorPair OrientedReadIds that have negative offsets.
    void removeNegativeOffsets(const Shasta2Anchors&);

    bool contains(OrientedReadId) const;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorIdA;
        ar & anchorIdB;
        ar & orientedReadIds;
    }
};
