// AssemblerTheseusReadWindowMSA.cpp
//
// Diagnostic prototype: partition physical reads into disjoint one-hop overlap
// windows. The Theseus MSA execution block is currently disabled while we
// verify window creation.

#include "Assembler.hpp"
#include "Reads.hpp"
#include "mode3-Anchor.hpp"
#include "timestamp.hpp"

#if 0
#include <theseus/heuristics.h>
#include <theseus/penalties.h>
#include <theseus/theseus_msa_aligner.h>
#endif

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <sstream>
#include <thread>
#include <vector>

using namespace dinara;
using namespace dinara::mode3;
using namespace std;

namespace {

constexpr array<uint64_t, 14> histogramUpperBounds = {
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
    numeric_limits<uint64_t>::max()
};

constexpr bool runTheseusMsa = false;
constexpr uint32_t minSharedAnchorsForRescue = 4;

struct ReadWindowTask {
    uint32_t windowId = 0;
    ReadId backboneReadId;
    vector<OrientedReadId> orientedReads;
    vector<ReadId> claimedReads;
    vector<uint32_t> alignmentIds;
};

struct EvidenceOccurrence {
    uint32_t windowId = 0;
    uint32_t row = 0;
    uint32_t alignmentId = 0;
};

struct CandidateEvidence {
    OrientedReadId orientedReadId;
    uint32_t alignmentId = 0;
};

struct AnchorWindowReadInterval {
    OrientedReadId orientedReadId;
    uint32_t begin = 0; // Inclusive position in the oriented read journey.
    uint32_t end = 0;   // Exclusive position in the oriented read journey.
    uint32_t sharedBackboneAnchors = 0;
};

struct AnchorWindowTask {
    uint32_t windowId = 0;
    OrientedReadId backboneOrientedReadId;
    uint32_t backboneBegin = 0;
    uint32_t backboneEnd = 0;
    uint64_t claimedAnchorCount = 0;
    vector<AnchorWindowReadInterval> readIntervals;
};

struct AnchorWindowCandidate {
    OrientedReadId backboneOrientedReadId;
    uint32_t begin = 0;
    uint32_t end = 0;
    uint64_t baseSpan = 0;
    uint64_t readLength = 0;
    uint32_t generation = 0;
};

struct AnchorWindowCandidateLess {
    bool operator()(const AnchorWindowCandidate& a, const AnchorWindowCandidate& b) const
    {
        if(a.baseSpan != b.baseSpan) {
            return a.baseSpan < b.baseSpan;
        }
        const uint32_t anchorCountA = a.end - a.begin;
        const uint32_t anchorCountB = b.end - b.begin;
        if(anchorCountA != anchorCountB) {
            return anchorCountA < anchorCountB;
        }
        if(a.readLength != b.readLength) {
            return a.readLength < b.readLength;
        }
        return a.backboneOrientedReadId.getReadId() > b.backboneOrientedReadId.getReadId();
    }
};

struct ThreadCounters {
    uint64_t windows = 0;
    uint64_t skippedSmallWindows = 0;
    uint64_t rows = 0;
    uint64_t bases = 0;
    double msaSeconds = 0.;
};

template<class T> void printWrappedItems(
    ostream& s,
    const string& prefix,
    const vector<T>& items,
    uint64_t itemsPerLine)
{
    s << timestamp << prefix;
    for(uint64_t i=0; i<items.size(); i++) {
        if((i % itemsPerLine) == 0) {
            s << "\n" << timestamp << "  ";
        }
        s << items[i] << " ";
    }
    s << endl;
}

void addToHistogram(array<uint64_t, histogramUpperBounds.size()>& histogram, uint64_t value)
{
    for(uint64_t i=0; i<histogramUpperBounds.size(); i++) {
        if(value <= histogramUpperBounds[i]) {
            ++histogram[i];
            return;
        }
    }
}

string histogramToString(const array<uint64_t, histogramUpperBounds.size()>& histogram)
{
    ostringstream s;
    uint64_t previous = 0;
    for(uint64_t i=0; i<histogram.size(); i++) {
        if(i != 0) {
            s << ",";
        }
        const uint64_t upper = histogramUpperBounds[i];
        if(upper == numeric_limits<uint64_t>::max()) {
            s << ">" << previous;
        } else if(previous + 1 == upper) {
            s << upper;
        } else {
            s << (previous + 1) << "-" << upper;
        }
        s << ":" << histogram[i];
        previous = upper;
    }
    return s.str();
}

string extractWholeOrientedReadSequence(const Reads& reads, OrientedReadId oid)
{
    const uint32_t length = uint32_t(reads.getRead(oid.getReadId()).baseCount);
    string sequence;
    sequence.reserve(length);
    for(uint32_t pos=0; pos<length; pos++) {
        sequence.push_back(reads.getOrientedReadBase(oid, pos).character());
    }
    return sequence;
}

} // namespace



void Assembler::computeTheseusReadWindowMSAPrototype(uint64_t threadCount)
{
    cout << timestamp << "[TheseusReadWindowMSA] Prototype begins." << endl;
    const auto totalBegin = chrono::steady_clock::now();

    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();
    checkReadGraphIsOpen();
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    DINARA_ASSERT(assemblerInfo.isOpen);
    DINARA_ASSERT((assemblerInfo->k % 2) == 0);
    DINARA_ASSERT(reads->readCount() > 0);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    threadCount = max<uint64_t>(1, threadCount);

    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    const uint32_t invalidAlignmentId = numeric_limits<uint32_t>::max();

    auto anchors = make_shared<Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        markerGraph,
        2,
        numeric_limits<uint64_t>::max(),
        threadCount,
        true);
    anchors->computeJourneys(threadCount);

    {
    // Anchor-interval window prototype.
    // Claim anchor ids, not whole reads. Each accepted window is seeded by a
    // contiguous unclaimed interval on a strand-0 backbone journey, then expanded
    // to contiguous unclaimed intervals on reads touching the seed anchors.
    constexpr uint32_t minBackboneWindowAnchors = 2;
    const uint32_t anchorUnclaimed = numeric_limits<uint32_t>::max();
    const uint64_t anchorCount = anchors->size();
    vector<uint32_t> anchorOwner(anchorCount, anchorUnclaimed);
    vector<AnchorWindowTask> anchorWindows;
    vector<uint32_t> touchedEpoch(orientedReadCount, 0);
    vector<uint32_t> touchedMin(orientedReadCount, numeric_limits<uint32_t>::max());
    vector<uint32_t> touchedMax(orientedReadCount, 0);
    vector<uint32_t> touchedCount(orientedReadCount, 0);
    vector<uint32_t> touchedOrientedReads;
    uint32_t epoch = 0;
    priority_queue<
        AnchorWindowCandidate,
        vector<AnchorWindowCandidate>,
        AnchorWindowCandidateLess> candidateHeap;
    vector<uint32_t> candidateGeneration(readCount, 0);

    vector<ReadId> anchorReadsByLength;
    anchorReadsByLength.reserve(readCount);
    for(uint64_t readId=0; readId<readCount; readId++) {
        anchorReadsByLength.push_back(ReadId(readId));
    }
    sort(anchorReadsByLength.begin(), anchorReadsByLength.end(),
        [&](ReadId a, ReadId b) {
            const uint64_t lengthA = reads->getRead(a).baseCount;
            const uint64_t lengthB = reads->getRead(b).baseCount;
            if(lengthA != lengthB) {
                return lengthA > lengthB;
            }
            return a < b;
        });

    uint64_t anchorWindowClaimedAnchors = 0;
    uint64_t anchorWindowReadIntervals = 0;
    uint64_t anchorWindowSkippedNoJourney = 0;
    uint64_t anchorWindowSkippedShortRuns = 0;
    uint64_t anchorWindowBackboneIntervals = 0;
    uint64_t anchorWindowTouchedReads = 0;
    uint64_t anchorWindowSplitIntervals = 0;
    uint64_t candidateIntervalsPushed = 0;
    uint64_t candidateIntervalsPopped = 0;
    uint64_t staleCandidateIntervals = 0;
    uint64_t discardedOldGenerationCandidates = 0;
    uint64_t maxAnchorWindowReadIntervals = 0;
    uint64_t maxAnchorWindowClaimedAnchors = 0;
    const auto anchorWindowBegin = chrono::steady_clock::now();

    auto intervalBaseSpan = [&](OrientedReadId oid, const auto& journey, uint32_t begin, uint32_t end) {
        if(begin >= end || oid.getValue() >= markers->size()) {
            return uint64_t(end - begin);
        }
        const auto orientedReadMarkers = (*markers)[oid.getValue()];
        if(orientedReadMarkers.empty()) {
            return uint64_t(end - begin);
        }

        const AnchorId leftAnchor = journey[begin];
        const AnchorId rightAnchor = journey[end - 1];
        uint64_t leftOrdinal = anchors->getFirstOrdinal(leftAnchor, oid);
        uint64_t rightOrdinal = uint64_t(anchors->getFirstOrdinal(rightAnchor, oid)) +
            anchors->ordinalOffset(rightAnchor);

        leftOrdinal = min<uint64_t>(leftOrdinal, orientedReadMarkers.size() - 1);
        const uint64_t leftPosition = orientedReadMarkers[leftOrdinal].position;
        uint64_t rightPosition = reads->getRead(oid.getReadId()).baseCount;
        if(rightOrdinal < orientedReadMarkers.size()) {
            rightPosition = orientedReadMarkers[rightOrdinal].position;
        }
        return rightPosition > leftPosition ? rightPosition - leftPosition : uint64_t(end - begin);
    };

    auto pushCandidate = [&](OrientedReadId oid, const auto& journey, uint32_t begin, uint32_t end) {
        if(end - begin < minBackboneWindowAnchors) {
            ++anchorWindowSkippedShortRuns;
            return;
        }
        candidateHeap.push(AnchorWindowCandidate{
            oid,
            begin,
            end,
            intervalBaseSpan(oid, journey, begin, end),
            reads->getRead(oid.getReadId()).baseCount,
            candidateGeneration[uint64_t(oid.getReadId())]});
        ++candidateIntervalsPushed;
    };

    auto pushCurrentUnclaimedIntervals = [&](OrientedReadId oid) {
        if(oid.getValue() >= anchors->journeys.size()) {
            ++anchorWindowSkippedNoJourney;
            return;
        }
        const auto journey = anchors->journeys[oid.getValue()];
        if(journey.empty()) {
            ++anchorWindowSkippedNoJourney;
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

    auto createAnchorWindow = [&](OrientedReadId backboneOid, uint32_t seedBegin, uint32_t seedEnd) {
        const uint32_t windowId = uint32_t(anchorWindows.size());
        AnchorWindowTask task;
        task.windowId = windowId;
        task.backboneOrientedReadId = backboneOid;
        task.backboneBegin = seedBegin;
        task.backboneEnd = seedEnd;
        task.readIntervals.push_back(AnchorWindowReadInterval{
            backboneOid,
            seedBegin,
            seedEnd,
            uint32_t(seedEnd - seedBegin)});

        const auto backboneJourney = anchors->journeys[backboneOid.getValue()];
        for(uint32_t position=seedBegin; position<seedEnd; position++) {
            const AnchorId anchorId = backboneJourney[position];
            if(anchorOwner[uint64_t(anchorId)] == anchorUnclaimed) {
                anchorOwner[uint64_t(anchorId)] = windowId;
                ++task.claimedAnchorCount;
            }
        }

        ++epoch;
        touchedOrientedReads.clear();
        for(uint32_t position=seedBegin; position<seedEnd; position++) {
            const AnchorId anchorId = backboneJourney[position];
            const Anchor anchor = (*anchors)[anchorId];
            for(const AnchorMarkerInterval& ami: anchor) {
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
        anchorWindowTouchedReads += touchedOrientedReads.size();

        for(const uint32_t oidValue: touchedOrientedReads) {
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            if(oid.getValue() >= anchors->journeys.size()) {
                continue;
            }
            const auto journey = anchors->journeys[oid.getValue()];
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
                    ++task.claimedAnchorCount;
                    ++position;
                }
                if(runBegin != position) {
                    task.readIntervals.push_back(AnchorWindowReadInterval{
                        oid,
                        runBegin,
                        position,
                        touchedCount[oidValue]});
                    ++anchorWindowSplitIntervals;
                }
            }
        }

        anchorWindowClaimedAnchors += task.claimedAnchorCount;
        anchorWindowReadIntervals += task.readIntervals.size();
        maxAnchorWindowReadIntervals = max<uint64_t>(
            maxAnchorWindowReadIntervals,
            task.readIntervals.size());
        maxAnchorWindowClaimedAnchors = max<uint64_t>(
            maxAnchorWindowClaimedAnchors,
            task.claimedAnchorCount);
        anchorWindows.push_back(std::move(task));
    };

    // Initialize the heap with one complete strand-0 journey candidate per read.
    for(const ReadId readId: anchorReadsByLength) {
        const OrientedReadId backboneOid(readId, 0);
        if(backboneOid.getValue() >= anchors->journeys.size()) {
            ++anchorWindowSkippedNoJourney;
            continue;
        }
        const auto journey = anchors->journeys[backboneOid.getValue()];
        if(journey.empty()) {
            ++anchorWindowSkippedNoJourney;
            continue;
        }
        pushCandidate(backboneOid, journey, 0, uint32_t(journey.size()));
    }

    while(!candidateHeap.empty()) {
        const AnchorWindowCandidate candidate = candidateHeap.top();
        candidateHeap.pop();
        ++candidateIntervalsPopped;

        const ReadId readId = candidate.backboneOrientedReadId.getReadId();
        if(candidate.generation != candidateGeneration[uint64_t(readId)]) {
            ++discardedOldGenerationCandidates;
            continue;
        }
        if(candidate.backboneOrientedReadId.getValue() >= anchors->journeys.size()) {
            ++anchorWindowSkippedNoJourney;
            continue;
        }
        const auto journey = anchors->journeys[candidate.backboneOrientedReadId.getValue()];
        if(journey.empty() || candidate.end > journey.size()) {
            ++anchorWindowSkippedNoJourney;
            continue;
        }

        bool isStillUnclaimed = true;
        for(uint32_t position=candidate.begin; position<candidate.end; position++) {
            if(anchorOwner[uint64_t(journey[position])] != anchorUnclaimed) {
                isStillUnclaimed = false;
                break;
            }
        }

        if(!isStillUnclaimed) {
            ++staleCandidateIntervals;
            ++candidateGeneration[uint64_t(readId)];
            pushCurrentUnclaimedIntervals(candidate.backboneOrientedReadId);
            continue;
        }

        ++anchorWindowBackboneIntervals;
        createAnchorWindow(candidate.backboneOrientedReadId, candidate.begin, candidate.end);
    }

    uint64_t unclaimedAnchorCount = 0;
    for(const uint32_t owner: anchorOwner) {
        if(owner == anchorUnclaimed) {
            ++unclaimedAnchorCount;
        }
    }

    array<uint64_t, histogramUpperBounds.size()> readIntervalsPerWindowHistogram = {};
    array<uint64_t, histogramUpperBounds.size()> claimedAnchorsPerWindowHistogram = {};
    array<uint64_t, histogramUpperBounds.size()> backboneAnchorSpanHistogram = {};
    for(const AnchorWindowTask& task: anchorWindows) {
        addToHistogram(readIntervalsPerWindowHistogram, task.readIntervals.size());
        addToHistogram(claimedAnchorsPerWindowHistogram, task.claimedAnchorCount);
        addToHistogram(backboneAnchorSpanHistogram, task.backboneEnd - task.backboneBegin);
    }

    cout << timestamp << "[TheseusReadWindowMSA] Top 10 read journey anchor usage:" << endl;
    for(uint64_t i=0; i<min<uint64_t>(10, anchorReadsByLength.size()); i++) {
        const ReadId readId = anchorReadsByLength[i];
        const OrientedReadId oid(readId, 0);
        uint64_t journeyAnchorCount = 0;
        uint64_t claimedJourneyAnchors = 0;
        uint64_t claimedRuns = 0;
        uint64_t unclaimedRuns = 0;
        if(oid.getValue() < anchors->journeys.size()) {
            const auto journey = anchors->journeys[oid.getValue()];
            journeyAnchorCount = journey.size();
            bool previousClaimed = false;
            bool previousUnclaimed = false;
            for(const AnchorId anchorId: journey) {
                const bool isClaimed = anchorOwner[uint64_t(anchorId)] != anchorUnclaimed;
                if(isClaimed) {
                    ++claimedJourneyAnchors;
                    if(!previousClaimed) {
                        ++claimedRuns;
                    }
                } else if(!previousUnclaimed) {
                    ++unclaimedRuns;
                }
                previousClaimed = isClaimed;
                previousUnclaimed = !isClaimed;
            }
        }
        cout << timestamp << "  rank=" << i
             << " readId=" << readId
             << " length=" << reads->getRead(readId).baseCount
             << " journeyAnchors=" << journeyAnchorCount
             << " claimedJourneyAnchors=" << claimedJourneyAnchors
             << " unclaimedJourneyAnchors=" << (journeyAnchorCount - claimedJourneyAnchors)
             << " claimedFraction=" << (journeyAnchorCount == 0 ? 0. : double(claimedJourneyAnchors) / double(journeyAnchorCount))
             << " claimedRuns=" << claimedRuns
             << " unclaimedRuns=" << unclaimedRuns
             << endl;
    }

    const auto anchorWindowEnd = chrono::steady_clock::now();
    const double anchorWindowSeconds = chrono::duration<double>(anchorWindowEnd - anchorWindowBegin).count();
    const auto totalEndAnchorPrototype = chrono::steady_clock::now();
    const double totalAnchorPrototypeSeconds = chrono::duration<double>(totalEndAnchorPrototype - totalBegin).count();

    cout << timestamp << "[TheseusReadWindowMSA] Anchor-window prototype ends."
         << " reads=" << readCount
         << " anchors=" << anchorCount
         << " windows=" << anchorWindows.size()
         << " minBackboneWindowAnchors=" << minBackboneWindowAnchors
         << " claimedAnchors=" << anchorWindowClaimedAnchors
         << " unclaimedAnchors=" << unclaimedAnchorCount
         << " backboneIntervals=" << anchorWindowBackboneIntervals
         << " candidateIntervalsPushed=" << candidateIntervalsPushed
         << " candidateIntervalsPopped=" << candidateIntervalsPopped
         << " staleCandidateIntervals=" << staleCandidateIntervals
         << " discardedOldGenerationCandidates=" << discardedOldGenerationCandidates
         << " readIntervals=" << anchorWindowReadIntervals
         << " splitReadIntervals=" << anchorWindowSplitIntervals
         << " touchedReads=" << anchorWindowTouchedReads
         << " maxReadIntervalsPerWindow=" << maxAnchorWindowReadIntervals
         << " avgReadIntervalsPerWindow=" << (anchorWindows.empty() ? 0. : double(anchorWindowReadIntervals) / double(anchorWindows.size()))
         << " maxClaimedAnchorsPerWindow=" << maxAnchorWindowClaimedAnchors
         << " avgClaimedAnchorsPerWindow=" << (anchorWindows.empty() ? 0. : double(anchorWindowClaimedAnchors) / double(anchorWindows.size()))
         << " readIntervalsPerWindowHistogram=" << histogramToString(readIntervalsPerWindowHistogram)
         << " claimedAnchorsPerWindowHistogram=" << histogramToString(claimedAnchorsPerWindowHistogram)
         << " backboneAnchorSpanHistogram=" << histogramToString(backboneAnchorSpanHistogram)
         << " skippedNoJourney=" << anchorWindowSkippedNoJourney
         << " skippedShortRuns=" << anchorWindowSkippedShortRuns
         << " anchorWindowSeconds=" << fixed << setprecision(6) << anchorWindowSeconds
         << " totalSeconds=" << totalAnchorPrototypeSeconds
         << defaultfloat << endl;
    return;
    }

    vector<ReadId> readsByLength;
    readsByLength.reserve(readCount);
    for(uint64_t readId=0; readId<readCount; readId++) {
        readsByLength.push_back(ReadId(readId));
    }
    sort(readsByLength.begin(), readsByLength.end(),
        [&](ReadId a, ReadId b) {
            const uint64_t lengthA = reads->getRead(a).baseCount;
            const uint64_t lengthB = reads->getRead(b).baseCount;
            if(lengthA != lengthB) {
                return lengthA > lengthB;
            }
            return a < b;
        });

    const uint32_t unclaimed = numeric_limits<uint32_t>::max();
    vector<uint32_t> readOwner(readCount, unclaimed);
    vector<ReadWindowTask> windows;
    uint64_t crossWindowEdgeCount = 0;
    uint64_t claimedReadCount = 0;
    uint64_t scannedReadGraphEdges = 0;
    uint64_t skippedCrossStrandEdges = 0;
    uint64_t skippedInconsistentEdges = 0;
    uint64_t skippedSelfEdges = 0;
    uint64_t borrowedReadCount = 0;
    uint64_t rejectedBackboneCandidates = 0;
    uint64_t rejectedClaimedNeighborEdges = 0;

    const auto planBegin = chrono::steady_clock::now();
    for(const ReadId seedReadId: readsByLength) {
        if(readOwner[uint64_t(seedReadId)] != unclaimed) {
            continue;
        }

        vector<CandidateEvidence> candidateEvidence;
        bool touchesClaimedRead = false;
        const OrientedReadId seedOid(seedReadId, 0);
        for(const uint32_t edgeId: readGraph.connectivity[seedOid.getValue()]) {
            ++scannedReadGraphEdges;
            const ReadGraphEdge& edge = readGraph.edges[edgeId];
            if(edge.crossesStrands) {
                ++skippedCrossStrandEdges;
                continue;
            }
            if(edge.hasInconsistentAlignment) {
                ++skippedInconsistentEdges;
                continue;
            }
            const OrientedReadId other = edge.getOther(seedOid);
            const uint64_t otherReadId = uint64_t(other.getReadId());
            if(otherReadId == uint64_t(seedReadId)) {
                ++skippedSelfEdges;
                continue;
            }

            if(readOwner[otherReadId] != unclaimed) {
                touchesClaimedRead = true;
                ++rejectedClaimedNeighborEdges;
            }
            candidateEvidence.push_back(CandidateEvidence{
                other,
                uint32_t(edge.alignmentId)});
        }

        if(touchesClaimedRead) {
            ++rejectedBackboneCandidates;
            continue;
        }

        const uint32_t windowId = uint32_t(windows.size());
        ReadWindowTask task;
        task.windowId = windowId;
        task.backboneReadId = seedReadId;
        task.orientedReads.push_back(OrientedReadId(seedReadId, 0));
        task.claimedReads.push_back(seedReadId);
        readOwner[uint64_t(seedReadId)] = windowId;
        ++claimedReadCount;

        for(const CandidateEvidence& evidence: candidateEvidence) {
            task.orientedReads.push_back(evidence.orientedReadId);
            task.alignmentIds.push_back(evidence.alignmentId);

            const ReadId otherReadId = evidence.orientedReadId.getReadId();
            if(readOwner[uint64_t(otherReadId)] == unclaimed) {
                readOwner[uint64_t(otherReadId)] = windowId;
                ++claimedReadCount;
                task.claimedReads.push_back(otherReadId);
            }
        }

        windows.push_back(std::move(task));
    }
    const auto planEnd = chrono::steady_clock::now();
    const double planSeconds = chrono::duration<double>(planEnd - planBegin).count();

    const auto anchorRescueBegin = chrono::steady_clock::now();
    vector<uint32_t> sharedAnchorCount(orientedReadCount, 0);
    vector<uint32_t> touchedOrientedReads;
    uint64_t anchorRescueRows = 0;
    uint64_t anchorRescueWindows = 0;
    uint64_t anchorRescueSharedAnchorHits = 0;
    uint64_t anchorRescueSkippedNoJourney = 0;
    uint64_t anchorRescueSkippedClaimedReads = 0;
    uint64_t anchorRescueSkippedLowSharedAnchors = 0;

    for(ReadWindowTask& task: windows) {
        const OrientedReadId backboneOid(task.backboneReadId, 0);
        if(backboneOid.getValue() >= anchors->journeys.size()) {
            ++anchorRescueSkippedNoJourney;
            continue;
        }

        const auto journey = anchors->journeys[backboneOid.getValue()];
        if(journey.empty()) {
            ++anchorRescueSkippedNoJourney;
            continue;
        }

        touchedOrientedReads.clear();
        for(const AnchorId anchorId: journey) {
            const Anchor anchor = (*anchors)[anchorId];
            for(const AnchorMarkerInterval& ami: anchor) {
                const OrientedReadId oid = ami.orientedReadId;
                if(oid == backboneOid) {
                    continue;
                }
                const ReadId readId = oid.getReadId();
                if(readOwner[uint64_t(readId)] != unclaimed) {
                    continue;
                }
                const uint64_t oidValue = uint64_t(oid.getValue());
                if(sharedAnchorCount[oidValue] == 0) {
                    touchedOrientedReads.push_back(uint32_t(oidValue));
                }
                ++sharedAnchorCount[oidValue];
                ++anchorRescueSharedAnchorHits;
            }
        }

        uint64_t rescuedInWindow = 0;
        for(const uint32_t oidValue: touchedOrientedReads) {
            const uint32_t shared = sharedAnchorCount[oidValue];
            sharedAnchorCount[oidValue] = 0;
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            if(readOwner[uint64_t(oid.getReadId())] != unclaimed) {
                ++anchorRescueSkippedClaimedReads;
                continue;
            }
            if(shared < minSharedAnchorsForRescue) {
                ++anchorRescueSkippedLowSharedAnchors;
                continue;
            }

            task.orientedReads.push_back(oid);
            task.alignmentIds.push_back(invalidAlignmentId);
            ++anchorRescueRows;
            ++rescuedInWindow;
        }
        if(rescuedInWindow != 0) {
            ++anchorRescueWindows;
        }
    }
    const auto anchorRescueEnd = chrono::steady_clock::now();
    const double anchorRescueSeconds = chrono::duration<double>(anchorRescueEnd - anchorRescueBegin).count();

    uint64_t singletonWindowCount = 0;
    uint64_t maxClaimedReadCount = 0;
    uint64_t totalEvidenceReadCount = 0;
    uint64_t maxEvidenceReadCount = 0;
    uint64_t ownerMismatchCount = 0;
    vector<bool> isBackboneRead(readCount, false);
    for(const ReadWindowTask& task: windows) {
        maxClaimedReadCount = max<uint64_t>(maxClaimedReadCount, task.claimedReads.size());
        maxEvidenceReadCount = max<uint64_t>(maxEvidenceReadCount, task.orientedReads.size());
        totalEvidenceReadCount += task.orientedReads.size();
        if(task.claimedReads.size() == 1) {
            ++singletonWindowCount;
        }
        isBackboneRead[uint64_t(task.backboneReadId)] = true;
        for(const ReadId readId: task.claimedReads) {
            if(readOwner[uint64_t(readId)] != task.windowId) {
                ++ownerMismatchCount;
            }
        }
    }

    uint64_t backboneConflictEdgeCount = 0;
    for(const ReadWindowTask& task: windows) {
        const OrientedReadId seedOid(task.backboneReadId, 0);
        for(const uint32_t edgeId: readGraph.connectivity[seedOid.getValue()]) {
            const ReadGraphEdge& edge = readGraph.edges[edgeId];
            if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                continue;
            }
            const ReadId otherReadId = edge.getOther(seedOid).getReadId();
            if(isBackboneRead[uint64_t(otherReadId)] &&
                readOwner[uint64_t(otherReadId)] != task.windowId) {
                ++backboneConflictEdgeCount;
            }
        }
    }
    backboneConflictEdgeCount /= 2;

    vector<vector<EvidenceOccurrence> > orientedReadOccurrences(2 * readCount);
    for(const ReadWindowTask& task: windows) {
        for(uint64_t row=0; row<task.orientedReads.size(); row++) {
            const OrientedReadId oid = task.orientedReads[row];
            const uint32_t alignmentId = (row == 0) ?
                invalidAlignmentId :
                task.alignmentIds[row - 1];
            orientedReadOccurrences[uint64_t(oid.getValue())].push_back(EvidenceOccurrence{
                task.windowId,
                uint32_t(row),
                alignmentId});
        }
    }

    uint64_t indexedOrientedReadCount = 0;
    uint64_t sharedOrientedReadCount = 0;
    uint64_t maxOrientedReadOccurrenceCount = 0;
    for(const vector<EvidenceOccurrence>& occurrences: orientedReadOccurrences) {
        if(occurrences.empty()) {
            continue;
        }
        ++indexedOrientedReadCount;
        maxOrientedReadOccurrenceCount = max<uint64_t>(
            maxOrientedReadOccurrenceCount,
            occurrences.size());
        if(occurrences.size() > 1) {
            ++sharedOrientedReadCount;
        }
    }

    const ReadId debugReadId = ReadId(2109);
    const uint32_t debugOwner = readOwner[uint64_t(debugReadId)];
    cout << timestamp << "[TheseusReadWindowMSA] Debug physical read " << debugReadId << "\n"
         << timestamp << "  length=" << reads->getRead(debugReadId).baseCount << "\n"
         << timestamp << "  ownerWindow=" << debugOwner;
    if(debugOwner != unclaimed) {
        const ReadWindowTask& ownerTask = windows[debugOwner];
        cout << "\n"
             << timestamp << "  ownerBackbone=" << OrientedReadId(ownerTask.backboneReadId, 0) << "\n"
             << timestamp << "  ownerBackboneLength=" << reads->getRead(ownerTask.backboneReadId).baseCount << "\n"
             << timestamp << "  ownerClaimedReads=" << ownerTask.claimedReads.size() << "\n"
             << timestamp << "  ownerEvidenceRows=" << ownerTask.orientedReads.size();
    }
    cout << endl;

    vector<bool> printedDebugWindow(windows.size(), false);
    uint64_t debugOccurrenceCount = 0;
    for(Strand strand=0; strand<2; strand++) {
        const OrientedReadId debugOrientedReadId(debugReadId, strand);
        const vector<EvidenceOccurrence>& debugOccurrences =
            orientedReadOccurrences[uint64_t(debugOrientedReadId.getValue())];
        cout << timestamp << "[TheseusReadWindowMSA] Occurrences of "
             << debugOrientedReadId
             << " count=" << debugOccurrences.size();
        for(const EvidenceOccurrence& occurrence: debugOccurrences) {
            ++debugOccurrenceCount;
            cout << "\n"
                 << timestamp << "  windowId=" << occurrence.windowId
                 << " row=" << occurrence.row;
            if(occurrence.alignmentId == invalidAlignmentId) {
            cout << " noAlignmentId";
            } else {
                cout << " alignmentId=" << occurrence.alignmentId;
            }
        }
        cout << endl;
    }

    for(Strand strand=0; strand<2; strand++) {
        const OrientedReadId debugOrientedReadId(debugReadId, strand);
        const vector<EvidenceOccurrence>& debugOccurrences =
            orientedReadOccurrences[uint64_t(debugOrientedReadId.getValue())];
        for(const EvidenceOccurrence& occurrence: debugOccurrences) {
            if(printedDebugWindow[occurrence.windowId]) {
                continue;
            }
            printedDebugWindow[occurrence.windowId] = true;

            const ReadWindowTask& task = windows[occurrence.windowId];
            cout << timestamp << "[TheseusReadWindowMSA] Window containing physical read "
                 << debugReadId << "\n"
                 << timestamp << "  windowId=" << task.windowId << "\n"
                 << timestamp << "  backbone=" << OrientedReadId(task.backboneReadId, 0) << "\n"
                 << timestamp << "  backboneLength=" << reads->getRead(task.backboneReadId).baseCount << "\n"
                 << timestamp << "  claimedReads=" << task.claimedReads.size() << "\n"
                 << timestamp << "  evidenceRows=" << task.orientedReads.size() << "\n"
                 << timestamp << "  debugReadOwner=" << debugOwner << "\n"
                 << timestamp << "  isOwnerWindow=" << (task.windowId == debugOwner)
                 << endl;

            printWrappedItems(cout, "[TheseusReadWindowMSA] Claimed reads:", task.claimedReads, 24);

            cout << timestamp << "[TheseusReadWindowMSA] Evidence rows:";
            for(uint64_t row=0; row<task.orientedReads.size(); row++) {
                if((row % 8) == 0) {
                    cout << "\n" << timestamp << "  ";
                }
                cout << " " << row << ":" << task.orientedReads[row];
                if(row > 0) {
                    if(task.alignmentIds[row - 1] == invalidAlignmentId) {
                        cout << "(anchorRescue)";
                    } else {
                        cout << "(alignmentId=" << task.alignmentIds[row - 1] << ")";
                    }
                } else {
                    cout << "(backbone)";
                }
            }
            cout << endl;
        }
    }
    if(debugOccurrenceCount == 0) {
        cout << timestamp << "[TheseusReadWindowMSA] No window contains physical read "
             << debugReadId
             << " as an evidence row."
             << endl;
    }

    vector<ThreadCounters> threadCounters(threadCount);
    double msaWallSeconds = 0.;
#if 0
    atomic<uint64_t> nextWindow(0);
    const auto msaBegin = chrono::steady_clock::now();
    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        threads.emplace_back([&, threadId]() {
            ThreadCounters& counters = threadCounters[threadId];
            theseus::Penalties penalties(0, 2, 3, 1);
            theseus::Heuristics heuristics(false, false);
            while(true) {
                const uint64_t windowIndex = nextWindow.fetch_add(1);
                if(windowIndex >= windows.size()) {
                    break;
                }
                const ReadWindowTask& task = windows[windowIndex];
                if(task.orientedReads.size() < 2) {
                    ++counters.skippedSmallWindows;
                    continue;
                }

                vector<string> sequences;
                sequences.reserve(task.orientedReads.size());
                uint64_t baseCount = 0;
                for(const OrientedReadId oid: task.orientedReads) {
                    sequences.push_back(extractWholeOrientedReadSequence(getReads(), oid));
                    baseCount += sequences.back().size();
                }
                if(sequences.empty() || sequences.front().empty()) {
                    ++counters.skippedSmallWindows;
                    continue;
                }

                const auto begin = chrono::steady_clock::now();
                theseus::TheseusMSA aligner(
                    penalties,
                    heuristics,
                    sequences.front(),
                    1,
                    false);
                for(size_t i=1; i<sequences.size(); i++) {
                    if(!sequences[i].empty()) {
                        aligner.align(sequences[i], 1, false, true);
                    }
                }
                ostringstream discard;
                aligner.print_as_msa(discard);
                const auto end = chrono::steady_clock::now();

                ++counters.windows;
                counters.rows += sequences.size();
                counters.bases += baseCount;
                counters.msaSeconds += chrono::duration<double>(end - begin).count();
            }
        });
    }
    for(thread& t: threads) {
        t.join();
    }
    const auto msaEnd = chrono::steady_clock::now();
    msaWallSeconds = chrono::duration<double>(msaEnd - msaBegin).count();
#endif

    ThreadCounters totalCounters;
    for(const ThreadCounters& counters: threadCounters) {
        totalCounters.windows += counters.windows;
        totalCounters.skippedSmallWindows += counters.skippedSmallWindows;
        totalCounters.rows += counters.rows;
        totalCounters.bases += counters.bases;
        totalCounters.msaSeconds += counters.msaSeconds;
    }

    const auto totalEnd = chrono::steady_clock::now();
    const double totalSeconds = chrono::duration<double>(totalEnd - totalBegin).count();
    cout << timestamp << "[TheseusReadWindowMSA] Prototype ends."
         << " reads=" << readCount
         << " alignments=" << alignmentData.size()
         << " readGraphEdges=" << readGraph.edges.size()
         << " scannedReadGraphEdges=" << scannedReadGraphEdges
         << " skippedCrossStrandEdges=" << skippedCrossStrandEdges
         << " skippedInconsistentEdges=" << skippedInconsistentEdges
         << " skippedSelfEdges=" << skippedSelfEdges
         << " windows=" << windows.size()
         << " claimedReads=" << claimedReadCount
         << " unclaimedReads=" << (readCount - claimedReadCount)
         << " rejectedBackboneCandidates=" << rejectedBackboneCandidates
         << " rejectedClaimedNeighborEdges=" << rejectedClaimedNeighborEdges
         << " anchorRescueMinSharedAnchors=" << minSharedAnchorsForRescue
         << " anchorRescueRows=" << anchorRescueRows
         << " anchorRescueWindows=" << anchorRescueWindows
         << " anchorRescueSharedAnchorHits=" << anchorRescueSharedAnchorHits
         << " anchorRescueSkippedNoJourney=" << anchorRescueSkippedNoJourney
         << " anchorRescueSkippedClaimedReads=" << anchorRescueSkippedClaimedReads
         << " anchorRescueSkippedLowSharedAnchors=" << anchorRescueSkippedLowSharedAnchors
         << " singletonWindows=" << singletonWindowCount
         << " maxClaimedReadsPerWindow=" << maxClaimedReadCount
         << " avgClaimedReadsPerWindow=" << (windows.empty() ? 0. : double(claimedReadCount) / double(windows.size()))
         << " evidenceReads=" << totalEvidenceReadCount
         << " borrowedEvidenceReads=" << borrowedReadCount
         << " maxEvidenceReadsPerWindow=" << maxEvidenceReadCount
         << " avgEvidenceReadsPerWindow=" << (windows.empty() ? 0. : double(totalEvidenceReadCount) / double(windows.size()))
         << " indexedOrientedReads=" << indexedOrientedReadCount
         << " sharedOrientedReads=" << sharedOrientedReadCount
         << " maxOrientedReadOccurrences=" << maxOrientedReadOccurrenceCount
         << " crossWindowEdges=" << crossWindowEdgeCount
         << " backboneConflictEdges=" << backboneConflictEdgeCount
         << " ownerMismatches=" << ownerMismatchCount
         << " runMsa=" << runTheseusMsa
         << " processedWindows=" << totalCounters.windows
         << " skippedSmallWindows=" << totalCounters.skippedSmallWindows
         << " rows=" << totalCounters.rows
         << " bases=" << totalCounters.bases
         << " planSeconds=" << fixed << setprecision(6) << planSeconds
         << " anchorRescueSeconds=" << anchorRescueSeconds
         << " msaThreadSeconds=" << totalCounters.msaSeconds
         << " msaWallSeconds=" << msaWallSeconds
         << " totalSeconds=" << totalSeconds
         << " avgSecondsPerWindow=" << (windows.empty() ? 0. : totalSeconds / double(windows.size()))
         << " avgMsaThreadSecondsPerWindow=" << (totalCounters.windows ? totalCounters.msaSeconds / double(totalCounters.windows) : 0.)
         << " threadCount=" << threadCount
         << defaultfloat << endl;
}
