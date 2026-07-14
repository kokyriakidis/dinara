#include "Shasta2Journeys.hpp"
#include "Shasta2Anchors.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace dinara;

#include <algorithm>
#include <thread>

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
    journeysWithOrdinals.createNew(largeDataName("tmp-Shasta2JourneysWithOrdinals"), largeDataPageSize);
    journeysWithOrdinals.beginPass1(orientedReadCount);
    setupLoadBalancing(anchorCount, anchorBatchCount);
    runThreads(&Shasta2Journeys::threadFunction1, threadCount);

    // Pass2: store the unsorted journeysWithOrdinals.
    journeysWithOrdinals.beginPass2();
    setupLoadBalancing(anchorCount, anchorBatchCount);
    runThreads(&Shasta2Journeys::threadFunction2, threadCount);
    journeysWithOrdinals.endPass2();

    // Pass 3:sort the journeysWithOrdinals and make space for the journeys
    journeys.createNew(largeDataName("Shasta2Journeys"), largeDataPageSize);
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



// Per-read journey filtering: keep the longest chain of anchors where every
// consecutive pair has >= minCommonForBackbone common reads (forward flow),
// with a bounded look-back of maxSkipForBackbone. Each read is handled
// independently; the surviving anchor ids are written to filteredJourneys.
void Shasta2Journeys::filterThreadFunction(uint64_t /* threadId */)
{
    const Shasta2Anchors& anchors = *anchorsPointer;
    const uint64_t minCommon = filterMinCommon;
    const uint64_t maxSkip = std::max<uint64_t>(filterMaxSkip, 1);

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t oidValue = begin; oidValue != end; oidValue++) {
            const auto journey = journeys[oidValue];
            const uint32_t n = uint32_t(journey.size());
            std::vector<Shasta2AnchorId>& out = filteredJourneys[oidValue];
            out.clear();
            if(n == 0) continue;
            if(n == 1) { out.push_back(journey[0]); continue; }

            // DP over the journey: dp[i] = length of the longest chain ending at
            // position i; prev[i] = predecessor position (-1 if the chain starts
            // at i). A pair (j, i) is chainable when countCommon(journey[j],
            // journey[i]) >= minCommon, considering only j in [i-maxSkip, i-1].
            std::vector<uint32_t> dp(n, 1);
            std::vector<int32_t> prev(n, -1);
            uint32_t bestEnd = 0;

            for(uint32_t i = 1; i < n; i++) {
                const Shasta2AnchorId anchorI = journey[i];
                const uint32_t lookBack = uint32_t(std::min<uint64_t>(i, maxSkip));
                for(uint32_t step = 1; step <= lookBack; step++) {
                    const uint32_t j = i - step;
                    if(dp[j] + 1 <= dp[i]) continue;  // cannot improve
                    const Shasta2AnchorId anchorJ = journey[j];
                    if(anchors.countCommon(anchorJ, anchorI) >= minCommon) {
                        dp[i] = dp[j] + 1;
                        prev[i] = int32_t(j);
                    }
                }
                if(dp[i] > dp[bestEnd]) bestEnd = i;
            }

            // Reconstruct the best chain (positions), then map to anchor ids.
            std::vector<uint32_t> chain;
            for(int32_t idx = int32_t(bestEnd); idx >= 0; idx = prev[uint32_t(idx)]) {
                chain.push_back(uint32_t(idx));
            }
            out.reserve(chain.size());
            for(auto it = chain.rbegin(); it != chain.rend(); ++it) {
                out.push_back(journey[*it]);
            }
        }
    }
}



void Shasta2Journeys::filterByAnchorChaining(
    uint64_t minCommonForBackbone,
    uint64_t maxSkipForBackbone,
    uint64_t threadCount)
{
    performanceLog << timestamp << "Journey filtering (anchor chaining) begins." << endl;

    // The Anchors pointer is retained only when this object was built via the
    // initial-creation constructor. Filtering rewrites both the journeys and
    // the per-anchor positionInJourney, so it needs mutable Anchors access.
    DINARA_ASSERT(anchorsPointer);
    DINARA_ASSERT(journeys.isOpen());

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t orientedReadCount = journeys.size();
    const uint64_t orientedReadBatchCount = 1000;

    // Pass A: compute the filtered chain for every oriented read in parallel.
    filterMinCommon = minCommonForBackbone;
    filterMaxSkip = maxSkipForBackbone;
    filteredJourneys.assign(orientedReadCount, {});
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Shasta2Journeys::filterThreadFunction, threadCount);

    // Pass B: rebuild the journeys VectorOfVectors in place from
    // filteredJourneys. The old content has already been captured into
    // filteredJourneys (pass A), so the fixed-size mmap storage is dropped and
    // recreated under the same name (two-pass count-then-store). The storage is
    // not movable, so it is rebuilt in place rather than swapped.
    journeys.remove();
    journeys.createNew(largeDataName("Shasta2Journeys"), largeDataPageSize);
    journeys.beginPass1(orientedReadCount);
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        journeys.incrementCount(oidValue, filteredJourneys[oidValue].size());
    }
    journeys.beginPass2();
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        for(const Shasta2AnchorId anchorId : filteredJourneys[oidValue]) {
            journeys.store(oidValue, anchorId);
        }
    }
    journeys.endPass2();
    filteredJourneys.clear();
    filteredJourneys.shrink_to_fit();

    // Pass C: reconcile positionInJourney for every anchor's marker infos. Reset
    // all of this run's reads to invalid first, then set the surviving positions
    // from the rebuilt journeys. Serial to keep the reset/set ordering simple;
    // this is a single linear sweep over anchor marker infos plus the journeys.
    Shasta2Anchors& anchors = *anchorsPointer;
    for(uint64_t anchorId = 0; anchorId < anchors.anchorMarkerInfos.size(); anchorId++) {
        for(Shasta2AnchorMarkerInfo& markerInfo : anchors.anchorMarkerInfos[anchorId]) {
            markerInfo.positionInJourney = invalid<uint32_t>;
        }
    }
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(oidValue));
        const auto journey = journeys[oidValue];
        for(uint64_t position = 0; position < journey.size(); position++) {
            const Shasta2AnchorId anchorId = journey[position];
            span<Shasta2AnchorMarkerInfo> markerInfos = anchors.anchorMarkerInfos[anchorId];
            for(Shasta2AnchorMarkerInfo& markerInfo : markerInfos) {
                if(markerInfo.orientedReadId == orientedReadId) {
                    markerInfo.positionInJourney = uint32_t(position);
                    break;
                }
            }
        }
    }

    performanceLog << timestamp << "Journey filtering (anchor chaining) ends." << endl;
}



// Access from binary data.
Shasta2Journeys::Shasta2Journeys(const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Shasta2Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner)
{
    journeys.accessExistingReadOnly(largeDataName("Shasta2Journeys"));
}
