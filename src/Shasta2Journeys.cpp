#include "Shasta2Journeys.hpp"
#include "Shasta2Anchors.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace dinara;

// Explicit instantiation.
#include "MultithreadedObject.tpp"
namespace dinara {
    template class MultithreadedObject<Shasta2Journeys>;
}



// Initial creation.
// This sets the positionInJourney for every AnchorMarkerInfo
// stored in the Anchors, and for this reason the Anchors
// are not passed in as const.
Shasta2Journeys::Shasta2Journeys(
    uint64_t orientedReadCount,
    shared_ptr<Shasta2Anchors> anchorsPointer,
    uint64_t threadCount,
    const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Shasta2Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    anchorsPointer(anchorsPointer)

{
    performanceLog << timestamp << "Journeys creation begins." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t anchorCount = anchorsPointer->size();
    const uint64_t anchorBatchCount = 1000;
    const uint64_t orientedReadBatchCount = 1000;

    // Pass1: make space for the journeysWithOrdinals.
    journeysWithOrdinals.createNew(largeDataName("tmp-JourneysWithOrdinals"), largeDataPageSize);
    journeysWithOrdinals.beginPass1(orientedReadCount);
    setupLoadBalancing(anchorCount, anchorBatchCount);
    runThreads(&Shasta2Journeys::threadFunction1, threadCount);

    // Pass2: store the unsorted journeysWithOrdinals.
    journeysWithOrdinals.beginPass2();
    setupLoadBalancing(anchorCount, anchorBatchCount);
    runThreads(&Shasta2Journeys::threadFunction2, threadCount);
    journeysWithOrdinals.endPass2();

    // Pass 3:sort the journeysWithOrdinals and make space for the journeys
    journeys.createNew(largeDataName("Journeys"), largeDataPageSize);
    journeys.beginPass1(orientedReadCount);
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Shasta2Journeys::threadFunction3, threadCount);

    // Pass 4: copy the sorted journeysWithOrdinals to the journeys.
    journeys.beginPass2();
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Shasta2Journeys::threadFunction4, threadCount);
    journeys.endPass2(false, true);

    journeysWithOrdinals.remove();

    performanceLog << timestamp << "Journeys creation ends." << endl;
}



void Shasta2Journeys::threadFunction1(uint64_t /* threadId */)
{
    threadFunction12(1);
}



void Shasta2Journeys::threadFunction2(uint64_t /* threadId */)
{
    threadFunction12(2);
}



void Shasta2Journeys::threadFunction12(uint64_t pass)
{
    const Shasta2Anchors& anchors = *anchorsPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all AnchorIds in this batch.
        for(Shasta2AnchorId anchorId=begin; anchorId!=end; anchorId++) {
            Shasta2Anchor anchor = anchors[anchorId];

            // Loop over the marker intervals of this Anchor.
            for(const auto& anchorMarkerInterval: anchor) {
                const auto orientedReadIdValue = anchorMarkerInterval.orientedReadId.getValue();

                if(pass == 1) {
                    journeysWithOrdinals.incrementCountMultithreaded(orientedReadIdValue);
                } else {
                    journeysWithOrdinals.storeMultithreaded(
                        orientedReadIdValue, {anchorId, anchorMarkerInterval.ordinal});
                }
            }
        }
    }
}



void Shasta2Journeys::threadFunction3(uint64_t /* threadId */)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            auto v = journeysWithOrdinals[orientedReadValue];
            sort(v.begin(), v.end(), OrderPairsBySecondOnly<uint64_t, uint32_t>());
            journeys.incrementCountMultithreaded(orientedReadValue, v.size());
        }
    }
}



void Shasta2Journeys::threadFunction4(uint64_t /* threadId */)
{
    Shasta2Anchors& anchors = *anchorsPointer;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(orientedReadValue));

            // Copy the journeysWithOrdinals to the journeys.
            const auto v = journeysWithOrdinals[orientedReadValue];
            const auto journey = journeys[orientedReadValue];
            DINARA_ASSERT(journey.size() == v.size());
            for(uint64_t i=0; i<v.size(); i++) {
                journey[i] = v[i].first;
            }

            // Store journey information for this oriented read in the marker interval.
            for(uint64_t position=0; position<journey.size(); position++) {
                const Shasta2AnchorId anchorId = journey[position];
                span<Shasta2AnchorMarkerInfo> markerInfos = anchors.anchorMarkerInfos[anchorId];
                bool found = false;
                for(Shasta2AnchorMarkerInfo& markerInfo: markerInfos) {
                    if(markerInfo.orientedReadId == orientedReadId) {
                        markerInfo.positionInJourney = uint32_t(position);
                        found = true;
                        break;
                    }
                }
                DINARA_ASSERT(found);
            }
        }
    }
}



// Access from binary data.
Shasta2Journeys::Shasta2Journeys(const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Shasta2Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner)
{
    journeys.accessExistingReadOnly(largeDataName("Journeys"));
}
