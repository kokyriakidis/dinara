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
