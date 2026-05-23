
// Dinara.
#include "Assembler.hpp"
#include "AnchorWindows.hpp"
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

} // anonymous namespace



void Assembler::computeAnchorWindowsClean(
    shared_ptr<Shasta2Anchors> shasta2Anchors,
    shared_ptr<Shasta2Journeys> shasta2Journeys,
    const vector<ReadId>& readIdsSortedByLength,
    vector<AnchorWindow>& anchorWindows,
    uint64_t threadCount)
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

    // Push a candidate interval onto the heap if it meets the minimum anchor count.
    auto pushCandidate = [&](OrientedReadId oid, const auto& journey, uint32_t begin, uint32_t end) {
        if(end - begin < minBackboneWindowAnchors) {
            return;
        }
        candidateHeap.push(CleanWindowCandidate{
            oid,
            begin,
            end,
            intervalBaseSpan(oid, journey, begin, end),
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
    auto longestIncreasingSubsequence = [](const vector<uint32_t>& seq) -> vector<uint32_t> {
        const uint32_t n = uint32_t(seq.size());
        if(n == 0) return {};

        // tails[i] = smallest ending value of an increasing subsequence of length i+1
        vector<uint32_t> tails;
        // tailIdx[i] = index in seq of the element at tails[i]
        vector<uint32_t> tailIdx;
        // prev[i] = index in seq of the predecessor of seq[i] in the LIS
        vector<int32_t> prev(n, -1);

        for(uint32_t i = 0; i < n; i++) {
            // Binary search for the position where seq[i] should go.
            auto it = std::lower_bound(tails.begin(), tails.end(), seq[i]);
            uint32_t pos = uint32_t(it - tails.begin());
            if(pos == tails.size()) {
                tails.push_back(seq[i]);
                tailIdx.push_back(i);
            } else {
                tails[pos] = seq[i];
                tailIdx[pos] = i;
            }
            if(pos > 0) {
                prev[i] = int32_t(tailIdx[pos - 1]);
            }
        }

        // Reconstruct the LIS.
        vector<uint32_t> result(tails.size());
        int32_t idx = int32_t(tailIdx.back());
        for(int32_t k = int32_t(result.size()) - 1; k >= 0; k--) {
            result[k] = uint32_t(idx);
            idx = prev[idx];
        }
        return result;
    };

    // ========================================================================
    // Create a window from a backbone interval.
    // For each touched read, keep only anchors shared with the backbone,
    // enforce backbone order via LIS, discard the rest.
    // ========================================================================
    auto createWindow = [&](OrientedReadId backboneOid, uint32_t seedBegin, uint32_t seedEnd) {
        const uint32_t windowId = uint32_t(anchorWindows.size());
        AnchorWindow window;
        window.windowId = windowId;
        window.backboneOrientedReadId = backboneOid;
        window.backboneBegin = seedBegin;
        window.backboneEnd = seedEnd;

        // Get backbone journey and build a lookup: anchorId -> position in backbone.
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];
        std::unordered_map<uint64_t, uint32_t> backboneAnchorToPos;
        backboneAnchorToPos.reserve(seedEnd - seedBegin);
        for(uint32_t pos = seedBegin; pos < seedEnd; pos++) {
            backboneAnchorToPos[uint64_t(backboneJourney[pos])] = pos;
        }

        // Claim backbone anchors.
        for(uint32_t pos = seedBegin; pos < seedEnd; pos++) {
            const Shasta2AnchorId anchorId = backboneJourney[pos];
            if(anchorOwner[uint64_t(anchorId)] == anchorUnclaimed) {
                anchorOwner[uint64_t(anchorId)] = windowId;
                ++window.claimedAnchorCount;
            }
        }

        // Add backbone read interval.
        window.readIntervals.push_back(AnchorWindowReadInterval{
            backboneOid,
            seedBegin,
            seedEnd,
            uint32_t(seedEnd - seedBegin)});

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

        // For each touched read:
        // 1. Find shared anchors (present in both the read's journey and the backbone).
        // 2. Map them to backbone positions.
        // 3. Compute LIS of backbone positions to enforce backbone order.
        // 4. Claim only the LIS anchors.
        for(const uint32_t oidValue : touchedOrientedReads) {
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            if(oid.getValue() >= shasta2Journeys->size()) continue;
            const auto journey = (*shasta2Journeys)[oid];
            if(journey.empty()) continue;

            // Step 1-2: Find shared anchors and their backbone positions.
            // sharedReadPositions[i] = position in the read's journey
            // sharedBackbonePositions[i] = corresponding position in the backbone
            vector<uint32_t> sharedReadPositions;
            vector<uint32_t> sharedBackbonePositions;

            for(uint32_t readPos = 0; readPos < uint32_t(journey.size()); readPos++) {
                const uint64_t anchorId = uint64_t(journey[readPos]);
                auto it = backboneAnchorToPos.find(anchorId);
                if(it != backboneAnchorToPos.end()) {
                    // This anchor is shared with the backbone and unclaimed.
                    if(anchorOwner[anchorId] == anchorUnclaimed ||
                       anchorOwner[anchorId] == windowId) {
                        sharedReadPositions.push_back(readPos);
                        sharedBackbonePositions.push_back(it->second);
                    }
                }
            }

            if(sharedReadPositions.empty()) continue;

            // Step 3: LIS of backbone positions to enforce backbone order.
            const auto lisIndices = longestIncreasingSubsequence(sharedBackbonePositions);

            if(lisIndices.empty()) continue;

            // Step 4: Claim the LIS anchors and record the interval.
            uint32_t convergentCount = 0;
            uint32_t convergentBegin = sharedReadPositions[lisIndices.front()];
            uint32_t convergentEnd = sharedReadPositions[lisIndices.back()] + 1;

            for(const uint32_t li : lisIndices) {
                const uint32_t readPos = sharedReadPositions[li];
                const uint64_t anchorId = uint64_t(journey[readPos]);
                if(anchorOwner[anchorId] == anchorUnclaimed) {
                    anchorOwner[anchorId] = windowId;
                    ++window.claimedAnchorCount;
                }
                ++convergentCount;
            }

            window.readIntervals.push_back(AnchorWindowReadInterval{
                oid,
                convergentBegin,
                convergentEnd,
                convergentCount});

            // Bump generation so stale candidates for this read are discarded.
            ++candidateGeneration[uint64_t(oid.getReadId())];
            // Re-push unclaimed intervals for this read.
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

    // Process the heap.
    while(!candidateHeap.empty()) {
        const CleanWindowCandidate candidate = candidateHeap.top();
        candidateHeap.pop();

        const ReadId readId = candidate.backboneOrientedReadId.getReadId();
        if(candidate.generation != candidateGeneration[uint64_t(readId)]) {
            continue;
        }
        if(candidate.backboneOrientedReadId.getValue() >= shasta2Journeys->size()) {
            continue;
        }
        const auto journey = (*shasta2Journeys)[candidate.backboneOrientedReadId];
        if(journey.empty() || candidate.end > journey.size()) {
            continue;
        }

        // Check if the interval is still fully unclaimed.
        bool isStillUnclaimed = true;
        for(uint32_t position = candidate.begin; position < candidate.end; position++) {
            if(anchorOwner[uint64_t(journey[position])] != anchorUnclaimed) {
                isStillUnclaimed = false;
                break;
            }
        }

        if(!isStillUnclaimed) {
            // Some anchors were claimed since this candidate was pushed.
            // Bump generation and re-push remaining unclaimed intervals.
            ++candidateGeneration[uint64_t(readId)];
            pushCurrentUnclaimedIntervals(candidate.backboneOrientedReadId);
            continue;
        }

        createWindow(candidate.backboneOrientedReadId, candidate.begin, candidate.end);
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
         << " anchors=" << anchorCount
         << " claimedAnchors=" << claimedAnchors
         << " unclaimedAnchors=" << unclaimedAnchorCount
         << " readIntervals=" << totalReadIntervals
         << " seconds=" << std::fixed << std::setprecision(2) << elapsedSeconds
         << std::defaultfloat << endl;
}
