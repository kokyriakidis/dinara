
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
#include <vector>



namespace {

// Priority queue candidate: a contiguous unclaimed interval on a backbone read's journey.
struct AnchorWindowCandidate {
    OrientedReadId backboneOrientedReadId;
    uint32_t begin;
    uint32_t end;
    uint64_t baseSpan;
    uint64_t readLength;
    uint32_t generation;
};

struct AnchorWindowCandidateLess {
    bool operator()(const AnchorWindowCandidate& a, const AnchorWindowCandidate& b) const {
        if(a.baseSpan != b.baseSpan) return a.baseSpan < b.baseSpan;
        if(a.readLength != b.readLength) return a.readLength < b.readLength;
        return a.backboneOrientedReadId.getValue() > b.backboneOrientedReadId.getValue();
    }
};

} // anonymous namespace



/*
Partition anchor journeys into disjoint windows.

Algorithm (from computeTheseusReadWindowMSAPrototype):
1. Seed a max-heap with full strand-0 journey intervals, prioritized by base span.
2. Pop the largest unclaimed interval. If any of its anchors have been claimed
   since it was pushed, re-push the remaining unclaimed sub-intervals and retry.
3. Create a window: claim the backbone anchors, then find all other reads
   touching those anchors and claim their unclaimed anchors too.
4. Repeat until the heap is empty.

Each window contains:
- The backbone read and its journey interval
- Read intervals for all other reads sharing anchors with the backbone
*/
void Assembler::computeAnchorWindows(
    shared_ptr<Shasta2Anchors> shasta2Anchors,
    shared_ptr<Shasta2Journeys> shasta2Journeys,
    const vector<ReadId>& readIdsSortedByLength,
    vector<AnchorWindow>& anchorWindows,
    uint64_t threadCount)
{
    cout << timestamp << "computeAnchorWindows begins." << endl;
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

    // Per-oriented-read scratch for tracking which reads touch the current backbone interval.
    vector<uint32_t> touchedEpoch(orientedReadCount, 0);
    vector<uint32_t> touchedMin(orientedReadCount, std::numeric_limits<uint32_t>::max());
    vector<uint32_t> touchedMax(orientedReadCount, 0);
    vector<uint32_t> touchedCount(orientedReadCount, 0);
    vector<uint32_t> touchedOrientedReads;
    uint32_t epoch = 0;

    // Priority queue and generation tracking for stale-candidate detection.
    std::priority_queue<
        AnchorWindowCandidate,
        vector<AnchorWindowCandidate>,
        AnchorWindowCandidateLess> candidateHeap;
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
        candidateHeap.push(AnchorWindowCandidate{
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

    // Create a window from a backbone interval and expand to touching reads.
    auto createWindow = [&](OrientedReadId backboneOid, uint32_t seedBegin, uint32_t seedEnd) {
        const uint32_t windowId = uint32_t(anchorWindows.size());
        AnchorWindow window;
        window.windowId = windowId;
        window.backboneOrientedReadId = backboneOid;
        window.backboneBegin = seedBegin;
        window.backboneEnd = seedEnd;
        window.readIntervals.push_back(AnchorWindowReadInterval{
            backboneOid,
            seedBegin,
            seedEnd,
            uint32_t(seedEnd - seedBegin)});

        // Claim backbone anchors.
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];
        for(uint32_t position = seedBegin; position < seedEnd; position++) {
            const Shasta2AnchorId anchorId = backboneJourney[position];
            if(anchorOwner[uint64_t(anchorId)] == anchorUnclaimed) {
                anchorOwner[uint64_t(anchorId)] = windowId;
                ++window.claimedAnchorCount;
            }
        }

        // Find all other reads touching the backbone anchors.
        ++epoch;
        touchedOrientedReads.clear();
        for(uint32_t position = seedBegin; position < seedEnd; position++) {
            const Shasta2AnchorId anchorId = backboneJourney[position];
            const Shasta2Anchor anchor = (*shasta2Anchors)[anchorId];
            for(const Shasta2AnchorMarkerInfo& ami : anchor) {
                const OrientedReadId oid = ami.orientedReadId;
                if(oid == backboneOid || ami.positionInJourney == invalid<uint32_t>) {
                    continue;
                }
                const uint32_t oidValue = uint32_t(oid.getValue());
                if(touchedEpoch[oidValue] != epoch) {
                    touchedEpoch[oidValue] = epoch;
                    touchedMin[oidValue] = ami.positionInJourney;
                    touchedMax[oidValue] = ami.positionInJourney;
                    touchedCount[oidValue] = 1;
                    touchedOrientedReads.push_back(oidValue);
                } else {
                    touchedMin[oidValue] = min(touchedMin[oidValue], ami.positionInJourney);
                    touchedMax[oidValue] = max(touchedMax[oidValue], ami.positionInJourney);
                    ++touchedCount[oidValue];
                }
            }
        }

        // For each touched read, claim unclaimed anchors in the touched range
        // and add read intervals.
        for(const uint32_t oidValue : touchedOrientedReads) {
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            if(oid.getValue() >= shasta2Journeys->size()) {
                continue;
            }
            const auto journey = (*shasta2Journeys)[oid];
            if(journey.empty()) {
                continue;
            }
            const uint32_t begin = touchedMin[oidValue];
            const uint32_t end = min<uint32_t>(touchedMax[oidValue] + 1, uint32_t(journey.size()));
            uint32_t position = begin;
            while(position < end) {
                while(position < end && anchorOwner[uint64_t(journey[position])] != anchorUnclaimed) {
                    ++position;
                }
                const uint32_t runBegin = position;
                while(position < end && anchorOwner[uint64_t(journey[position])] == anchorUnclaimed) {
                    anchorOwner[uint64_t(journey[position])] = windowId;
                    ++window.claimedAnchorCount;
                    ++position;
                }
                if(runBegin != position) {
                    window.readIntervals.push_back(AnchorWindowReadInterval{
                        oid,
                        runBegin,
                        position,
                        touchedCount[oidValue]});
                }
            }
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
        const AnchorWindowCandidate candidate = candidateHeap.top();
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

    cout << timestamp << "computeAnchorWindows ends."
         << " windows=" << anchorWindows.size()
         << " anchors=" << anchorCount
         << " claimedAnchors=" << claimedAnchors
         << " unclaimedAnchors=" << unclaimedAnchorCount
         << " readIntervals=" << totalReadIntervals
         << " seconds=" << std::fixed << std::setprecision(2) << elapsedSeconds
         << std::defaultfloat << endl;
}
