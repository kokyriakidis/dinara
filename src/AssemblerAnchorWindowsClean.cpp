
// Dinara.
#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "Kmer.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"
using namespace dinara;

// Standard libraries.
#include "chrono.hpp"
#include <algorithm>
#include <iomanip>
#include <limits>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>



namespace {

// Priority queue candidate: a contiguous unclaimed interval on a read's journey.
struct CleanWindowCandidate {
    OrientedReadId backboneOrientedReadId;
    uint32_t begin;
    uint32_t end;
    uint64_t baseSpan;
    uint64_t readLength;
    uint32_t generation;
};

struct CleanWindowCandidateLess {
    bool operator()(const CleanWindowCandidate& a, const CleanWindowCandidate& b) const {
        if(a.baseSpan != b.baseSpan) return a.baseSpan < b.baseSpan;
        if(a.readLength != b.readLength) return a.readLength < b.readLength;
        return a.backboneOrientedReadId.getValue() > b.backboneOrientedReadId.getValue();
    }
};

// Whole-journey window claiming (one hop). When true, a window claims not just
// the backbone-overlapping span of each touching read, but the read's ENTIRE
// journey (clamped only by what is still unclaimed, since ownership is
// exclusive). This pulls each window's left/right dovetails fully into the
// window instead of leaving them as the inter-window "mess", so downstream
// path-finding can be done inside the window. The toucher's readInterval (used
// for per-window MSA/consensus) stays at the backbone-overlapping span; only
// the ANCHOR CLAIMING span (which drives window membership / topology) is
// extended. When false, the original convergent-span claiming is used.
//
// NOTE: exclusive ownership + whole-journey claiming disconnects neighbors
// (connecting reads get fully absorbed by whichever window claims first) and
// assigns dovetails by an arbitrary claim-race. The intended model is a SHARED
// dovetail halo (a separate multi-owner structure) layered on disjoint
// exclusive cores; this exclusive toggle is kept off pending that halo work.
constexpr bool claimWholeJourneyDovetails = false;

} // anonymous namespace



void Assembler::computeAnchorWindowsClean(
    shared_ptr<Shasta2Anchors> shasta2Anchors,
    shared_ptr<Shasta2Journeys> shasta2Journeys,
    const vector<ReadId>& readIdsSortedByLength,
    vector<AnchorWindow>& anchorWindows,
    uint64_t threadCount,
    uint64_t minCommonForBackbone,
    uint64_t maxSkipForBackbone,
    uint64_t minWindowBaseSpan,
    vector<uint32_t>* anchorDovetailWindow)
{
    cout << timestamp << "computeAnchorWindowsClean begins." << endl;
    const auto t0 = steady_clock::now();

    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    DINARA_ASSERT(shasta2Anchors);
    DINARA_ASSERT(shasta2Journeys);
    DINARA_ASSERT(shasta2Journeys->isOpen());
    DINARA_ASSERT(reads->readCount() > 0);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    const uint64_t anchorCount = shasta2Anchors->size();

    constexpr uint32_t minBackboneWindowAnchors = 2;
    const uint32_t anchorUnclaimed = std::numeric_limits<uint32_t>::max();

    vector<uint32_t> anchorOwner(anchorCount, anchorUnclaimed);
    anchorWindows.clear();

    // Forward-oriented dovetail window membership, persisted to the caller so
    // the anchor graph can map claimed dovetail anchors to their owning window
    // (anchorOwner alone is local scratch and orientation-collapsed). Sized
    // only when whole-journey claiming is enabled and the caller requested it;
    // otherwise left empty so the old behavior is byte-identical downstream.
    constexpr uint32_t dovetailNoWindow = std::numeric_limits<uint32_t>::max();
    const bool recordDovetails =
        (anchorDovetailWindow != nullptr) && claimWholeJourneyDovetails;
    if(anchorDovetailWindow != nullptr) {
        anchorDovetailWindow->clear();
        if(claimWholeJourneyDovetails) {
            anchorDovetailWindow->assign(anchorCount, dovetailNoWindow);
        }
    }

    // Per-oriented-read scratch for tracking which reads touch the current backbone.
    vector<uint32_t> touchedEpoch(orientedReadCount, 0);
    vector<uint32_t> touchedOrientedReads;
    uint32_t epoch = 0;

    // Priority queue and generation tracking for stale-candidate detection.
    std::priority_queue<
        CleanWindowCandidate,
        vector<CleanWindowCandidate>,
        CleanWindowCandidateLess> candidateHeap;
    vector<uint32_t> candidateGeneration(readCount, 0);

    // Counters.
    uint64_t claimedAnchors = 0;
    uint64_t totalReadIntervals = 0;


    // Compute the base span of a journey interval [begin, end).
    auto intervalBaseSpan = [&](OrientedReadId oid, const auto& journey, uint32_t begin, uint32_t end) {
        if(begin >= end || oid.getValue() >= markers->size()) {
            return uint64_t(end - begin);
        }
        const auto orientedReadMarkers = (*markers)[oid.getValue()];
        if(orientedReadMarkers.empty()) {
            return uint64_t(end - begin);
        }

        const Shasta2AnchorId leftAnchorId = journey[begin];
        const Shasta2AnchorId rightAnchorId = journey[end - 1];
        const uint32_t leftOrdinal = shasta2Anchors->getOrdinal(leftAnchorId, oid);
        const uint32_t rightOrdinal = shasta2Anchors->getOrdinal(rightAnchorId, oid);
        if(leftOrdinal == invalid<uint32_t> || rightOrdinal == invalid<uint32_t>) {
            return uint64_t(end - begin);
        }

        const uint64_t leftOrdClamped =
            min<uint64_t>(leftOrdinal, orientedReadMarkers.size() - 1);
        const uint64_t leftPosition = orientedReadMarkers[leftOrdClamped].position;
        uint64_t rightPosition = reads->getRead(oid.getReadId()).baseCount;
        if(rightOrdinal < orientedReadMarkers.size()) {
            rightPosition = orientedReadMarkers[rightOrdinal].position;
        }
        return rightPosition > leftPosition ? rightPosition - leftPosition : uint64_t(end - begin);
    };

    // Push a candidate interval onto the heap if it meets minimum anchor count
    // and minimum base span.
    auto pushCandidate = [&](OrientedReadId oid, const auto& journey, uint32_t begin, uint32_t end) {
        if(end - begin < minBackboneWindowAnchors) {
            return;
        }
        const uint64_t span = intervalBaseSpan(oid, journey, begin, end);
        if(minWindowBaseSpan > 0 && span < minWindowBaseSpan) {
            return;
        }
        candidateHeap.push(CleanWindowCandidate{
            oid,
            begin,
            end,
            span,
            reads->getReadRawSequenceLength(oid.getReadId()),
            candidateGeneration[uint64_t(oid.getReadId())]});
    };

    // Re-push all current unclaimed intervals for a read.
    auto pushCurrentUnclaimedIntervals = [&](OrientedReadId oid) {
        if(oid.getValue() >= shasta2Journeys->size()) {
            return;
        }
        const auto journey = (*shasta2Journeys)[oid];
        if(journey.empty()) {
            return;
        }

        uint32_t position = 0;
        while(position < journey.size()) {
            while(position < journey.size() &&
                  anchorOwner[uint64_t(journey[position])] != anchorUnclaimed) {
                ++position;
            }
            const uint32_t runBegin = position;
            while(position < journey.size() &&
                  anchorOwner[uint64_t(journey[position])] == anchorUnclaimed) {
                ++position;
            }
            const uint32_t runEnd = position;
            if(runBegin != runEnd) {
                pushCandidate(oid, journey, runBegin, runEnd);
            }
        }
    };

    // ========================================================================
    // LIS helper: given a sequence of backbone positions (uint32_t),
    // return the indices (into the input) of the longest increasing subsequence.
    // Uses patience sorting (O(n log n)).
    // ========================================================================

    // Claim an anchor and its reverse complement.
    auto claimAnchor = [&](uint64_t anchorId, uint32_t windowId) {
        if(anchorOwner[anchorId] == anchorUnclaimed) {
            anchorOwner[anchorId] = windowId;
            const uint64_t rcAnchorId = anchorId ^ 1ULL;
            if(rcAnchorId < anchorCount) {
                anchorOwner[rcAnchorId] = windowId;
            }
        }
    };

    // ========================================================================
    // Create a window from a backbone interval.
    // For each touched read, find first/last shared anchor with backbone
    // and claim all anchors in between.
    // ========================================================================
    // Find k-mers that appear more than once in a journey interval.
    auto findDuplicateKmers = [&](const auto& journey, uint32_t begin, uint32_t end) {
        std::set<Kmer> seen;
        std::set<Kmer> duplicateKmers;
        for(uint32_t pos = begin; pos < end; pos++) {
            const Kmer kmer = shasta2Anchors->anchorKmer(journey[pos]);
            if(!seen.insert(kmer).second) {
                duplicateKmers.insert(kmer);
            }
        }
        return duplicateKmers;
    };

    // Find anchor IDs in a journey interval whose k-mer appears more than once.
    // Anchors at the first and last positions of the interval are never marked
    // as duplicates, to avoid disconnecting chain endpoints.
    auto findDuplicateKmerAnchors = [&](const auto& journey, uint32_t begin, uint32_t end) {
        const auto duplicateKmers = findDuplicateKmers(journey, begin, end);
        std::unordered_set<uint64_t> duplicateAnchorIds;
        if(!duplicateKmers.empty()) {
            for(uint32_t pos = begin; pos < end; pos++) {
                if(pos == begin || pos == end - 1) continue;
                const Kmer kmer = shasta2Anchors->anchorKmer(journey[pos]);
                if(duplicateKmers.count(kmer)) {
                    duplicateAnchorIds.insert(uint64_t(journey[pos]));
                }
            }
        }
        return duplicateAnchorIds;
    };

    // Filter a backbone journey interval to keep the longest subsequence
    // where every consecutive pair has at least minCommonForBackbone common reads.
    // Uses dynamic programming on a DAG with limited forward look.
    auto filterBackboneJourney = [&](
        OrientedReadId backboneOid,
        uint32_t seedBegin,
        uint32_t seedEnd,
        vector<uint32_t>& filteredPositions)
    {
        filteredPositions.clear();
        const auto journey = (*shasta2Journeys)[backboneOid];
        const uint32_t n = seedEnd - seedBegin;
        if(n == 0) return;
        if(n <= 2) {
            for(uint32_t pos = seedBegin; pos < seedEnd; pos++) {
                filteredPositions.push_back(pos);
            }
            return;
        }

        // Run LIS on the interior positions [1..n-2], keeping only
        // consecutive pairs with >= minCommonForBackbone common reads.
        // First and last positions are always kept.
        const uint32_t interiorN = n - 2;  // number of interior positions

        // dp[i] = length of longest consistent subsequence ending at interior position i.
        // prev[i] = predecessor index in the subsequence (-1 if none).
        vector<uint32_t> dp(interiorN, 0);
        vector<int32_t> prev(interiorN, -1);

        // Seed: interior positions reachable from the first anchor.
        const Shasta2AnchorId firstAnchor = journey[seedBegin];
        for(uint32_t i = 0; i < interiorN && i < maxSkipForBackbone; i++) {
            const Shasta2AnchorId anchorI = journey[seedBegin + 1 + i];
            if(shasta2Anchors->countCommon(firstAnchor, anchorI) >= minCommonForBackbone) {
                dp[i] = 1;
            }
        }

        // Fill DP for interior positions.
        for(uint32_t i = 1; i < interiorN; i++) {
            const Shasta2AnchorId anchorI = journey[seedBegin + 1 + i];
            const uint32_t lookBack = min(uint64_t(i), maxSkipForBackbone);
            for(uint32_t step = 1; step <= lookBack; step++) {
                const uint32_t j = i - step;
                if(dp[j] == 0) continue;  // j not reachable from first anchor
                const Shasta2AnchorId anchorJ = journey[seedBegin + 1 + j];
                if(shasta2Anchors->countCommon(anchorJ, anchorI) >= minCommonForBackbone) {
                    if(dp[j] + 1 > dp[i]) {
                        dp[i] = dp[j] + 1;
                        prev[i] = int32_t(j);
                    }
                }
            }
        }

        // Find the best interior endpoint that also connects to the last anchor.
        const Shasta2AnchorId lastAnchor = journey[seedEnd - 1];
        int32_t bestEnd = -1;
        for(int32_t i = int32_t(interiorN) - 1;
            i >= 0 && i >= int32_t(interiorN) - int32_t(maxSkipForBackbone);
            i--)
        {
            if(dp[i] == 0) continue;
            const Shasta2AnchorId anchorI = journey[seedBegin + 1 + uint32_t(i)];
            if(shasta2Anchors->countCommon(anchorI, lastAnchor) >= minCommonForBackbone) {
                if(bestEnd < 0 || dp[i] > dp[bestEnd]) {
                    bestEnd = i;
                }
            }
        }

        // Build result: first anchor + interior LIS + last anchor.
        filteredPositions.push_back(seedBegin);
        if(bestEnd >= 0) {
            vector<uint32_t> reversePath;
            for(int32_t idx = bestEnd; idx >= 0; idx = prev[idx]) {
                reversePath.push_back(seedBegin + 1 + uint32_t(idx));
            }
            filteredPositions.insert(filteredPositions.end(),
                reversePath.rbegin(), reversePath.rend());
        }
        filteredPositions.push_back(seedEnd - 1);
    };

    auto createWindow = [&](OrientedReadId backboneOid, uint32_t seedBegin, uint32_t seedEnd) {
        // Use the full backbone range — no filtering.
        const uint32_t n = seedEnd - seedBegin;
        if(n < minBackboneWindowAnchors) {
            return; // Too few anchors.
        }

        const auto fullJourney = (*shasta2Journeys)[backboneOid];

        vector<uint32_t> filteredPositions;
        filteredPositions.reserve(n);
        for(uint32_t pos = seedBegin; pos < seedEnd; pos++) {
            filteredPositions.push_back(pos);
        }

        const uint32_t windowId = uint32_t(anchorWindows.size());
        AnchorWindow window;
        window.windowId = windowId;
        window.backboneOrientedReadId = backboneOid;
        window.backboneBegin = seedBegin;
        window.backboneEnd = seedEnd;
        window.filteredBackbonePositions = filteredPositions;

        // Get backbone journey and build a lookup: anchorId -> position in backbone.
        // Only include anchors that survived filtering and don't have duplicate k-mers.
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];
        std::unordered_map<uint64_t, uint32_t> backboneAnchorToPos;
        backboneAnchorToPos.reserve(filteredPositions.size());

        // Find backbone anchors with duplicate k-mers and exclude them.
        const auto backboneDuplicates = findDuplicateKmerAnchors(
            backboneJourney, window.backboneBegin, window.backboneEnd);
        // Build the lookup from filtered positions only.
        std::unordered_set<uint64_t> filteredAnchorSet;
        for(const uint32_t pos : filteredPositions) {
            filteredAnchorSet.insert(uint64_t(backboneJourney[pos]));
        }
        for(const uint32_t pos : filteredPositions) {
            const uint64_t aid = uint64_t(backboneJourney[pos]);
            if(backboneDuplicates.count(aid) == 0) {
                backboneAnchorToPos[aid] = pos;
            }
        }

        // Add backbone read interval (covers the full filtered span).
        window.readIntervals.push_back(AnchorWindowReadInterval{
            backboneOid,
            window.backboneBegin,
            window.backboneEnd,
            uint32_t(filteredPositions.size())});

        // Collect all anchors to claim at the end.
        // Each entry is (orientedReadId, journey span begin, journey span end).
        // We store the spans so we can do a single claiming pass at the end.
        struct ReadSpan {
            OrientedReadId oid;
            uint32_t begin;
            uint32_t end;
        };
        vector<ReadSpan> readSpans;

        // The backbone span.
        readSpans.push_back(ReadSpan{backboneOid, seedBegin, seedEnd});

        // Find all other oriented reads that share anchors with the backbone.
        ++epoch;
        touchedOrientedReads.clear();
        for(uint32_t pos = seedBegin; pos < seedEnd; pos++) {
            const Shasta2AnchorId anchorId = backboneJourney[pos];
            const auto anchor = (*shasta2Anchors)[anchorId];
            for(const Shasta2AnchorMarkerInfo& ami : anchor) {
                const OrientedReadId oid = ami.orientedReadId;
                if(oid == backboneOid) continue;
                const uint32_t oidValue = uint32_t(oid.getValue());
                if(touchedEpoch[oidValue] != epoch) {
                    touchedEpoch[oidValue] = epoch;
                    touchedOrientedReads.push_back(oidValue);
                }
            }
        }

        // Build the set of direct cis overlap oriented read IDs for the backbone
        // using the read graph (which contains only cis overlaps).
        std::unordered_set<uint32_t> directCisOverlapReads;
        {
            const auto backboneEdges = readGraph.connectivity[backboneOid.getValue()];
            for(const uint32_t edgeId : backboneEdges) {
                const auto& edge = readGraph.edges[edgeId];
                const OrientedReadId partner = edge.getOther(backboneOid);
                directCisOverlapReads.insert(partner.getValue());
            }
        }

        // For each touched read:
        // 1. Find shared anchors (present in both the read's journey and the backbone).
        // 2. Map them to backbone positions.
        // 3. Compute LIS of backbone positions to enforce backbone order.
        // 4. Record the span for later claiming.
        // 5. For non-direct overlaps, extract alternate paths between LIS pillars.
        for(const uint32_t oidValue : touchedOrientedReads) {
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            if(oid.getValue() >= shasta2Journeys->size()) continue;
            const auto journey = (*shasta2Journeys)[oid];
            if(journey.empty()) continue;

            const bool isDirectCis = (directCisOverlapReads.find(oidValue) != directCisOverlapReads.end());

            // Find anchors with duplicate k-mers in this read's journey.
            const auto readDuplicates = findDuplicateKmerAnchors(journey, 0, uint32_t(journey.size()));

            // Step 1-2: Find shared anchors and their backbone positions.
            vector<uint32_t> sharedReadPositions;
            vector<uint32_t> sharedBackbonePositions;

            for(uint32_t readPos = 0; readPos < uint32_t(journey.size()); readPos++) {
                const uint64_t anchorId = uint64_t(journey[readPos]);
                if(readDuplicates.count(anchorId)) continue;
                auto it = backboneAnchorToPos.find(anchorId);
                if(it != backboneAnchorToPos.end()) {
                    if(anchorOwner[anchorId] == anchorUnclaimed) {
                        sharedReadPositions.push_back(readPos);
                        sharedBackbonePositions.push_back(it->second);
                    }
                }
            }

            if(sharedReadPositions.empty()) continue;

            // Use first and last shared positions to define the interval.
            // Claim all anchors in between regardless of order.
            const uint32_t convergentBegin = sharedReadPositions.front();
            const uint32_t convergentEnd = sharedReadPositions.back() + 1;
            const uint32_t convergentCount = uint32_t(sharedReadPositions.size());

            // Claiming span: whole journey (one hop) when enabled, so the
            // window absorbs this toucher's left/right dovetails. The
            // readInterval below (for MSA/consensus) stays at the convergent
            // backbone-overlapping span regardless.
            const uint32_t claimBegin = claimWholeJourneyDovetails ?
                0u : convergentBegin;
            const uint32_t claimEnd = claimWholeJourneyDovetails ?
                uint32_t(journey.size()) : convergentEnd;
            readSpans.push_back(ReadSpan{oid, claimBegin, claimEnd});

            window.readIntervals.push_back(AnchorWindowReadInterval{
                oid,
                convergentBegin,
                convergentEnd,
                convergentCount});

        }

        // Claim all anchors across all spans (backbone + touched reads) and their RC.
        for(const ReadSpan& span : readSpans) {
            const auto journey = (*shasta2Journeys)[span.oid];
            // Forward orientation of this read relative to the strand-0 backbone:
            // a strand-1 toucher's journey holds RC anchors, so the forward
            // anchor is anchorId^1. Used to record window membership in the
            // backbone-forward frame (matching anchorToWindow's convention).
            const bool readIsStrand1 = (span.oid.getStrand() == 1);
            for(uint32_t readPos = span.begin; readPos < span.end; readPos++) {
                const uint64_t anchorId = uint64_t(journey[readPos]);
                if(anchorOwner[anchorId] == anchorUnclaimed) {
                    claimAnchor(anchorId, windowId);
                    ++window.claimedAnchorCount;
                    // Record forward-oriented dovetail membership (toggle only).
                    if(recordDovetails) {
                        const uint64_t forwardAid =
                            readIsStrand1 ? (anchorId ^ 1ULL) : anchorId;
                        if(forwardAid < anchorCount) {
                            (*anchorDovetailWindow)[forwardAid] = windowId;
                        }
                    }
                }
            }
        }

        // Now that claiming is done, bump generations and re-push unclaimed
        // intervals for all touched reads.
        for(const uint32_t oidValue : touchedOrientedReads) {
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            ++candidateGeneration[uint64_t(oid.getReadId())];
            pushCurrentUnclaimedIntervals(oid);
        }

        claimedAnchors += window.claimedAnchorCount;
        totalReadIntervals += window.readIntervals.size();
        anchorWindows.push_back(std::move(window));
    };

    // Initialize the heap: all strand-0 reads, ordered by length (longest first).
    for(const ReadId readId : readIdsSortedByLength) {
        const OrientedReadId backboneOid(readId, 0);
        if(backboneOid.getValue() >= shasta2Journeys->size()) {
            continue;
        }
        const auto journey = (*shasta2Journeys)[backboneOid];
        if(journey.empty()) {
            continue;
        }
        pushCandidate(backboneOid, journey, 0, uint32_t(journey.size()));
    }

    // Process the heap. For each candidate (longest base span first):
    // - If the full journey is unclaimed, create one window for the whole journey.
    // - Otherwise, find contiguous unclaimed intervals within the journey
    //   and create a separate window for each.
    uint64_t fullJourneyWindows = 0;
    uint64_t fragmentWindows = 0;

    while(!candidateHeap.empty()) {
        const CleanWindowCandidate candidate = candidateHeap.top();
        candidateHeap.pop();

        const ReadId readId = candidate.backboneOrientedReadId.getReadId();
        if(candidate.generation != candidateGeneration[uint64_t(readId)]) {
            continue;
        }
        const auto journey = (*shasta2Journeys)[candidate.backboneOrientedReadId];

        // Find contiguous unclaimed intervals within the full journey.
        // We always scan the full journey regardless of candidate.begin/end,
        // since claimed regions may have appeared since this candidate was pushed.
        vector<pair<uint32_t, uint32_t>> unclaimedIntervals;
        {
            uint32_t pos = 0;
            const uint32_t journeySize = uint32_t(journey.size());
            while(pos < journeySize) {
                while(pos < journeySize &&
                      anchorOwner[uint64_t(journey[pos])] != anchorUnclaimed) {
                    ++pos;
                }
                const uint32_t runBegin = pos;
                while(pos < journeySize &&
                      anchorOwner[uint64_t(journey[pos])] == anchorUnclaimed) {
                    ++pos;
                }
                if(pos > runBegin && (pos - runBegin) >= minBackboneWindowAnchors) {
                    unclaimedIntervals.push_back({runBegin, pos});
                }
            }
        }

        if(unclaimedIntervals.empty()) continue;

        // Create a window for each unclaimed interval that passes the
        // base span threshold.
        const bool isFullJourney = (unclaimedIntervals.size() == 1 &&
            unclaimedIntervals[0].first == 0 &&
            unclaimedIntervals[0].second == uint32_t(journey.size()));

        for(const auto& [intervalBegin, intervalEnd] : unclaimedIntervals) {
            if(minWindowBaseSpan > 0 && (intervalEnd - intervalBegin) >= 2) {
                const Shasta2AnchorId firstAnchor = journey[intervalBegin];
                const Shasta2AnchorId lastAnchor = journey[intervalEnd - 1];
                const uint32_t firstPos = shasta2Anchors->getPosition(
                    firstAnchor, candidate.backboneOrientedReadId);
                const uint32_t lastPos = shasta2Anchors->getPosition(
                    lastAnchor, candidate.backboneOrientedReadId);
                const uint64_t baseSpan = (lastPos >= firstPos) ? (lastPos - firstPos) : 0;
                if(baseSpan < minWindowBaseSpan) {
                    continue;
                }
            }

            createWindow(candidate.backboneOrientedReadId, intervalBegin, intervalEnd);

            if(isFullJourney) {
                ++fullJourneyWindows;
            } else {
                ++fragmentWindows;
            }
        }
    }

    // Count unclaimed anchors.
    uint64_t unclaimedAnchorCount = 0;
    for(const uint32_t owner : anchorOwner) {
        if(owner == anchorUnclaimed) {
            ++unclaimedAnchorCount;
        }
    }

    const auto t1 = steady_clock::now();
    const double elapsedSeconds = seconds(t1 - t0);

    cout << timestamp << "computeAnchorWindowsClean ends."
         << " windows=" << anchorWindows.size()
         << " (fullJourney=" << fullJourneyWindows
         << " fragment=" << fragmentWindows << ")"
         << " anchors=" << anchorCount
         << " claimedAnchors=" << claimedAnchors
         << " unclaimedAnchors=" << unclaimedAnchorCount
         << " readIntervals=" << totalReadIntervals
         << " seconds=" << std::fixed << std::setprecision(2) << elapsedSeconds
         << std::defaultfloat << endl;
}
