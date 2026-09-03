#include "Shasta2Journeys.hpp"
#include "Shasta2Anchors.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace dinara;

#include <algorithm>
#include <iostream>
#include <thread>
using std::cout;
using std::endl;

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
            // Anchor members are stored sorted ascending by OrientedReadId (see
            // Shasta2Anchors), so binary search instead of scanning the whole anchor.
            for(uint64_t position=0; position<journey.size(); position++) {
                const Shasta2AnchorId anchorId = journey[position];
                span<Shasta2AnchorMarkerInfo> markerInfos = anchors.anchorMarkerInfos[anchorId];
                const auto it = std::lower_bound(markerInfos.begin(), markerInfos.end(), orientedReadId,
                    [](const Shasta2AnchorMarkerInfo& info, OrientedReadId oid) {
                        return info.orientedReadId < oid;
                    });
                DINARA_ASSERT(it != markerInfos.end() and it->orientedReadId == orientedReadId);
                it->positionInJourney = uint32_t(position);
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
            // Strand symmetry: the journey of oriented read (R,1) must be the
            // exact reverse complement of (R,0) -- journey(R,1)[k] ==
            // journey(R,0)[n-1-k] ^ 1. The DP tie-breaking is not mirror-
            // symmetric, so filtering the two strands independently would
            // diverge and break the coverage(a) == coverage(a^1) invariant that
            // downstream stages assert. Filter only strand 0 (even oidValue)
            // and mirror the result into strand 1 (oidValue ^ 1).
            if((oidValue & 1ULL) != 0ULL) continue;   // strand 1 done below.

            const auto journey = journeys[oidValue];
            const uint32_t n = uint32_t(journey.size());
            std::vector<Shasta2AnchorId>& out = filteredJourneys[oidValue];
            std::vector<Shasta2AnchorId>& outRc = filteredJourneys[oidValue ^ 1ULL];
            out.clear();
            outRc.clear();

            // Emit both strands: strand 0 = fwd, strand 1 = reversed + RC.
            auto emitMirror = [&]() {
                outRc.resize(out.size());
                for(size_t k = 0; k < out.size(); k++) {
                    outRc[k] = out[out.size() - 1 - k] ^ 1ULL;
                }
            };

            if(n == 0) { continue; }
            if(n == 1) { out.push_back(journey[0]); emitMirror(); continue; }

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
            emitMirror();
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

    // With minCommonForBackbone == 0, filterThreadFunction's chainability test
    // (countCommon(...) >= minCommon) is trivially true for any pair, for any
    // input -- not just on the data we happened to test. Its DP always ends up
    // picking the immediately preceding position (dp[i-1]+1 strictly beats any
    // dp[j]+1 for j further back, since dp is non-decreasing), so the "longest
    // chain" is always the full, untouched original journey. The whole call
    // (per-read DP, a full drop+recreate+two-pass rebuild of the journeys
    // VectorOfVectors, and an anchor-store-wide positionInJourney reconciliation
    // pass) is therefore provably a no-op here, not just empirically one -- skip
    // it rather than pay for a rewrite that changes nothing.
    if(minCommonForBackbone == 0) {
        cout << timestamp << "Journey filtering (anchor chaining): minCommonForBackbone=0, "
                "filter is a no-op by construction -- skipping." << endl;
        performanceLog << timestamp << "Journey filtering (anchor chaining) skipped "
                "(minCommonForBackbone=0)." << endl;
        return;
    }

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

    // Report how much the filter removed.
    {
        uint64_t anchorsBefore = 0, anchorsAfter = 0, collapsed = 0;
        for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
            anchorsBefore += journeys[oidValue].size();
            anchorsAfter += filteredJourneys[oidValue].size();
            if(journeys[oidValue].size() >= 2 && filteredJourneys[oidValue].size() < 2) collapsed++;
        }
        cout << timestamp << "Journey filtering: anchors " << anchorsBefore
            << " -> " << anchorsAfter << ", " << collapsed
            << " journeys collapsed below 2 anchors (minCommon=" << minCommonForBackbone
            << " maxSkip=" << maxSkipForBackbone << ")." << endl;
    }

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
    // Fill by direct index assignment, NOT store(): store() writes each vector
    // back-to-front (it decrements the per-index count), which would reverse
    // every journey and break the ordinal monotonicity window construction
    // relies on. beginPass2 has already allocated the exact space.
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        const auto journey = journeys[oidValue];
        const std::vector<Shasta2AnchorId>& filtered = filteredJourneys[oidValue];
        DINARA_ASSERT(journey.size() == filtered.size());
        for(uint64_t i = 0; i < filtered.size(); i++) {
            journey[i] = filtered[i];
        }
    }
    // count was consumed by incrementCount but store() never ran, so skip the
    // all-zero check (check=false); free the count vector (free=true).
    journeys.endPass2(false, true);
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
            // Anchor members are stored sorted ascending by OrientedReadId (see
            // Shasta2Anchors), so binary search instead of scanning the whole anchor.
            const auto it = std::lower_bound(markerInfos.begin(), markerInfos.end(), orientedReadId,
                [](const Shasta2AnchorMarkerInfo& info, OrientedReadId oid) {
                    return info.orientedReadId < oid;
                });
            if(it != markerInfos.end() and it->orientedReadId == orientedReadId) {
                it->positionInJourney = uint32_t(position);
            }
        }
    }

    performanceLog << timestamp << "Journey filtering (anchor chaining) ends." << endl;
}



// Rebuild the journeys and positionInJourney from scratch using the current
// anchor set, keyed by `position` instead of `ordinal` -- see the .hpp
// comment for why. Strand 1 is never sorted independently: it is derived by
// reversing strand 0's sorted list and flipping each anchor id, exactly like
// filterThreadFunction's emitMirror, so the journey(R,1) ==
// reverse(journey(R,0)) invariant holds regardless of how same-position ties
// on strand 0 happen to be broken (ascending anchorId here, arbitrary but
// deterministic). Runs serially: the anchor count here is small (tens of
// thousands), so a second full pass is cheap and avoids any thread-safety
// concerns in a rebuild that runs once, off the hot path.
void Shasta2Journeys::rebuildAfterNewAnchors(Shasta2AnchorId newAnchorsBegin, uint64_t threadCount)
{
    performanceLog << timestamp << "Journeys rebuild (new anchors) begins." << endl;
    DINARA_ASSERT(anchorsPointer);
    DINARA_ASSERT(journeys.isOpen());
    static_cast<void>(threadCount);

    Shasta2Anchors& anchors = *anchorsPointer;
    const uint64_t anchorCount = anchors.size();
    const uint64_t orientedReadCount = journeys.size();
    const uint64_t readCount = orientedReadCount / 2;
    const Shasta2AnchorId hetFirst = anchors.hetAnchorFirstId;

    // Pass A: collect, per read (strand 0 only), (position, anchorId) pairs by
    // scanning primary anchors and the new anchor range only -- see the .hpp
    // comment for why anchors in between (e.g. export-only SNPmer anchors)
    // are excluded.
    vector<vector<pair<uint32_t, Shasta2AnchorId>>> strand0(readCount);
    for(Shasta2AnchorId anchorId = 0; anchorId < anchorCount; anchorId++) {
        const bool isPrimary = (hetFirst == invalid<Shasta2AnchorId>) || (anchorId < hetFirst);
        const bool isNew = (anchorId >= newAnchorsBegin);
        if(!isPrimary && !isNew) continue;
        for(const Shasta2AnchorMarkerInfo& info : anchors.anchorMarkerInfos[anchorId]) {
            if((info.orientedReadId.getValue() & 1U) != 0U) continue;   // strand 1 done via mirror.
            strand0[info.orientedReadId.getReadId()].push_back({info.position, anchorId});
        }
    }

    // Pass B: sort each read's strand-0 list by position (ties broken by
    // ascending anchorId -- correctness does not depend on which tied anchor
    // sorts first, only that strand 1 is mirrored from strand 0, not sorted
    // independently), then derive strand 1 by reversing + flipping.
    filteredJourneys.assign(orientedReadCount, {});
    for(uint64_t readIdValue = 0; readIdValue < readCount; readIdValue++) {
        auto& v = strand0[readIdValue];
        std::sort(v.begin(), v.end());
        std::vector<Shasta2AnchorId>& out = filteredJourneys[2 * readIdValue];
        std::vector<Shasta2AnchorId>& outRc = filteredJourneys[2 * readIdValue + 1];
        out.reserve(v.size());
        for(const auto& [position, anchorId] : v) {
            static_cast<void>(position);
            out.push_back(anchorId);
        }
        outRc.resize(out.size());
        for(size_t k = 0; k < out.size(); k++) {
            outRc[k] = out[out.size() - 1 - k] ^ 1ULL;
        }
    }
    strand0.clear();
    strand0.shrink_to_fit();

    // Pass C: rebuild the journeys VectorOfVectors in place from
    // filteredJourneys -- identical shape to filterByAnchorChaining's Pass B.
    journeys.remove();
    journeys.createNew(largeDataName("Shasta2Journeys"), largeDataPageSize);
    journeys.beginPass1(orientedReadCount);
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        journeys.incrementCount(oidValue, filteredJourneys[oidValue].size());
    }
    journeys.beginPass2();
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        const auto journey = journeys[oidValue];
        const std::vector<Shasta2AnchorId>& filtered = filteredJourneys[oidValue];
        DINARA_ASSERT(journey.size() == filtered.size());
        for(uint64_t i = 0; i < filtered.size(); i++) {
            journey[i] = filtered[i];
        }
    }
    journeys.endPass2(false, true);
    filteredJourneys.clear();
    filteredJourneys.shrink_to_fit();

    // Pass D: reconcile positionInJourney for every anchor's marker infos
    // (reset all to invalid, then set from the rebuilt journeys), identical
    // shape to filterByAnchorChaining's Pass C.
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
            const auto it = std::lower_bound(markerInfos.begin(), markerInfos.end(), orientedReadId,
                [](const Shasta2AnchorMarkerInfo& info, OrientedReadId oid) {
                    return info.orientedReadId < oid;
                });
            if(it != markerInfos.end() and it->orientedReadId == orientedReadId) {
                it->positionInJourney = uint32_t(position);
            }
        }
    }

    performanceLog << timestamp << "Journeys rebuild (new anchors) ends." << endl;
}



// Access from binary data.
Shasta2Journeys::Shasta2Journeys(const MappedMemoryOwner& mappedMemoryOwner) :
    MultithreadedObject<Shasta2Journeys>(*this),
    MappedMemoryOwner(mappedMemoryOwner)
{
    journeys.accessExistingReadOnly(largeDataName("Shasta2Journeys"));
}
