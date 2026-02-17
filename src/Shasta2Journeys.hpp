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

};
