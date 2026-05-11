// AssemblerTheseusReadWindowMSA.cpp
//
// Diagnostic prototype: partition physical reads into disjoint one-hop overlap
// windows. The Theseus MSA execution block is currently disabled while we
// verify window creation.

#include "Assembler.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"

#if 0
#include <theseus/heuristics.h>
#include <theseus/penalties.h>
#include <theseus/theseus_msa_aligner.h>
#endif

#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <thread>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

constexpr bool runTheseusMsa = false;

struct ReadWindowTask {
    uint32_t windowId = 0;
    ReadId backboneReadId;
    vector<OrientedReadId> orientedReads;
    vector<ReadId> claimedReads;
    vector<uint32_t> alignmentIds;
};

struct ThreadCounters {
    uint64_t windows = 0;
    uint64_t skippedSmallWindows = 0;
    uint64_t rows = 0;
    uint64_t bases = 0;
    double msaSeconds = 0.;
};

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
    DINARA_ASSERT(reads->readCount() > 0);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    threadCount = max<uint64_t>(1, threadCount);

    const uint64_t readCount = reads->readCount();

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

    const auto planBegin = chrono::steady_clock::now();
    for(const ReadId seedReadId: readsByLength) {
        if(readOwner[uint64_t(seedReadId)] != unclaimed) {
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

            task.orientedReads.push_back(other);
            task.alignmentIds.push_back(uint32_t(edge.alignmentId));

            const uint32_t owner = readOwner[otherReadId];
            if(owner == unclaimed) {
                readOwner[otherReadId] = windowId;
                ++claimedReadCount;
                task.claimedReads.push_back(other.getReadId());
            } else if(owner != windowId) {
                ++crossWindowEdgeCount;
                ++borrowedReadCount;
            }
        }

        windows.push_back(std::move(task));
    }
    const auto planEnd = chrono::steady_clock::now();
    const double planSeconds = chrono::duration<double>(planEnd - planBegin).count();

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
         << " singletonWindows=" << singletonWindowCount
         << " maxClaimedReadsPerWindow=" << maxClaimedReadCount
         << " avgClaimedReadsPerWindow=" << (windows.empty() ? 0. : double(claimedReadCount) / double(windows.size()))
         << " evidenceReads=" << totalEvidenceReadCount
         << " borrowedEvidenceReads=" << borrowedReadCount
         << " maxEvidenceReadsPerWindow=" << maxEvidenceReadCount
         << " avgEvidenceReadsPerWindow=" << (windows.empty() ? 0. : double(totalEvidenceReadCount) / double(windows.size()))
         << " crossWindowEdges=" << crossWindowEdgeCount
         << " backboneConflictEdges=" << backboneConflictEdgeCount
         << " ownerMismatches=" << ownerMismatchCount
         << " runMsa=" << runTheseusMsa
         << " processedWindows=" << totalCounters.windows
         << " skippedSmallWindows=" << totalCounters.skippedSmallWindows
         << " rows=" << totalCounters.rows
         << " bases=" << totalCounters.bases
         << " planSeconds=" << fixed << setprecision(6) << planSeconds
         << " msaThreadSeconds=" << totalCounters.msaSeconds
         << " msaWallSeconds=" << msaWallSeconds
         << " totalSeconds=" << totalSeconds
         << " avgSecondsPerWindow=" << (windows.empty() ? 0. : totalSeconds / double(windows.size()))
         << " avgMsaThreadSecondsPerWindow=" << (totalCounters.windows ? totalCounters.msaSeconds / double(totalCounters.windows) : 0.)
         << " threadCount=" << threadCount
         << defaultfloat << endl;
}
