#pragma once

// Shasta2Journeys.hpp

#include "Shasta2Anchors.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "Reads.hpp"
#include "span.hpp"

#include <vector>
#include <memory>

namespace dinara {
    class Shasta2Journeys;
    using Shasta2Journey = span<const Shasta2AnchorId>;
}


class dinara::Shasta2Journeys :
    public MultithreadedObject<Shasta2Journeys>,
    public MappedMemoryOwner {
public:

    // Initial creation.
    // This sets the positionInJourney for every AnchorMarkerInfo
    // stored in the Anchors, and for this reason the Anchors
    // are not passed in as const.
    Shasta2Journeys(
        uint64_t orientedReadCount,
        shared_ptr<Shasta2Anchors>,
        uint64_t threadCount,
        const MappedMemoryOwner&);

    // Access from binary data.
    Shasta2Journeys(const MappedMemoryOwner&);

    // Rewrite each oriented read's journey in place, keeping only the longest
    // chain of anchors where every consecutive pair in the chain has at least
    // minCommonForBackbone common reads (forward flow), looking back at most
    // maxSkipForBackbone positions. Each read is filtered independently -- the
    // decision for one read never depends on any other read's journey (it uses
    // only pairwise anchor countCommon, which is a global anchor-pair property).
    // No endpoints are forced: a read whose anchors are weakly linked collapses
    // to its longest well-supported subsequence. The stored journeys and every
    // anchor's positionInJourney are rebuilt to match; anchors dropped from a
    // read's chain get positionInJourney = invalid for that read.
    // Requires the Anchors pointer retained from initial creation.
    void filterByAnchorChaining(
        uint64_t minCommonForBackbone,
        uint64_t maxSkipForBackbone,
        uint64_t threadCount);

    // Rebuild the journeys and every anchor's positionInJourney from scratch,
    // keyed by each occurrence's `position` field instead of `ordinal`.
    // Unlike ordinal, position is defined for anchors appended after initial
    // construction (e.g. Shasta2Anchors::appendHetAnchorPair het anchors,
    // whose ordinal is always invalid) and is RC-symmetric by construction, so
    // this lets newly-appended anchors take their correct place in every
    // affected read's journey with no special-casing.
    //
    // Only two anchor ranges are included: primary anchors (id <
    // anchors.hetAnchorFirstId) -- i.e. exactly what was already in the
    // journeys -- and the new anchors in [newAnchorsBegin, anchors.size()).
    // Any het/hom anchor in between (id in [hetAnchorFirstId,
    // newAnchorsBegin)) is skipped: those would be anchors some OTHER caller
    // appended for a purpose that does not include internal journey/anchor-
    // graph participation (e.g. an export-only bolt-on not built with the
    // position-margin safety this rebuild's callers rely on -- see
    // Shasta2AnchorGraphHetOnGraph.cpp's file header). Pass newAnchorsBegin ==
    // the anchor count captured right before appending the anchors this
    // rebuild is meant to fold in.
    //
    // Call after appending new anchors and before rebuilding any
    // Shasta2AnchorGraph from this object. Requires the Anchors pointer
    // retained from initial creation.
    void rebuildAfterNewAnchors(Shasta2AnchorId newAnchorsBegin, uint64_t threadCount);

    // Which anchor survives when two independent anchors land on the SAME base
    // of the same read. Keeping both is impossible (shasta2 asserts strictly
    // increasing positions), so one occurrence must be dropped from that
    // read's journey. true = keep the het anchor, false = keep the primary.
    // Defaults to true: preferring primary makes het anchors lose every
    // collision they are in, measured at 11.1% of all het-anchor occurrences,
    // whereas the primary that loses instead carries ~2x the coverage and so
    // pays proportionally much less. See Shasta2Journeys.cpp for the
    // measurement.
    bool journeyTiePreferHet = true;

    // The (canonical anchor, ReadId) units this rebuild removed from journeys.
    //
    // The external-anchor export must omit exactly these members, and it used to
    // work that out for itself -- first by re-resolving the ties (a second
    // implementation of the policy, which disagreed), then by rebuilding the
    // whole journey membership as a ~1.1M-entry hash set and diffing every
    // anchor member against it. Both are roundabout: the rebuild HAS this list
    // at the moment it decides, and it is tiny (457 units on the 989-read
    // fixture against 576k exported members). Recording it here makes the export
    // O(drops) instead of O(journey entries + anchor members), and removes the
    // last place where the export could infer something different from what the
    // rebuild actually did.
    vector< pair<Shasta2AnchorId, ReadId> > journeyTieDrops;

    // No exported anchor may fall below this many members. Tie resolution
    // removes members, so without the bound it can thin an anchor past the
    // coverage floor the anchor was selected under -- and an anchor exported
    // with one read is a spurious branch, not a weak one.
    uint64_t minAnchorCoverage = 2;

    // Return the Journey for an oriented read.
    Shasta2Journey operator[](OrientedReadId orientedReadId) const
    {
        return journeys[orientedReadId.getValue()];
    }

    uint64_t size() const
    {
        return journeys.size();
    }

    bool isOpen() const
    {
        return journeys.isOpen();
    }

private:

    MemoryMapped::VectorOfVectors<Shasta2AnchorId, uint64_t> journeys;

    // This is only stored during initial creation.
    shared_ptr<Shasta2Anchors> anchorsPointer;

    void threadFunction1(uint64_t threadId);
    void threadFunction2(uint64_t threadId);
    void threadFunction12(uint64_t pass);
    void threadFunction3(uint64_t threadId);
    void threadFunction4(uint64_t threadId);

    // Temporary storage of journeys with ordinals.
    MemoryMapped::VectorOfVectors<pair<uint64_t, uint32_t>, uint64_t> journeysWithOrdinals;

    // Transient state for filterByAnchorChaining. filteredJourneys[oidValue]
    // holds the surviving anchor ids for that oriented read; each entry is
    // written by exactly one thread, so no locking is needed.
    std::vector<std::vector<Shasta2AnchorId>> filteredJourneys;
    uint64_t filterMinCommon = 0;
    uint64_t filterMaxSkip = 0;
    void filterThreadFunction(uint64_t threadId);

};
