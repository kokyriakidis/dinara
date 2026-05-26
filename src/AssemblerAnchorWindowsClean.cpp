
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

} // anonymous namespace



void Assembler::computeAnchorWindowsClean(
    shared_ptr<Shasta2Anchors> shasta2Anchors,
    shared_ptr<Shasta2Journeys> shasta2Journeys,
    const vector<ReadId>& readIdsSortedByLength,
    vector<AnchorWindow>& anchorWindows,
    uint64_t threadCount,
    uint64_t minCommonForBackbone,
    uint64_t maxSkipForBackbone)
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
    uint64_t totalAltPathsBeforeFilter = 0;
    uint64_t totalAltPathsAfterFilter = 0;
    uint64_t totalIntermediatesConsidered = 0;
    uint64_t totalIntermediatesKept = 0;

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
    // For each touched read, keep only anchors shared with the backbone,
    // enforce backbone order via LIS, discard the rest.
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
        if(n == 1) {
            filteredPositions.push_back(seedBegin);
            return;
        }

        // dp[i] = length of longest consistent subsequence ending at position i.
        // prev[i] = predecessor index in the subsequence (-1 if none).
        vector<uint32_t> dp(n, 1);
        vector<int32_t> prev(n, -1);

        for(uint32_t i = 1; i < n; i++) {
            const Shasta2AnchorId anchorI = journey[seedBegin + i];
            // Look back up to maxSkipForBackbone positions.
            const uint32_t lookBack = min(uint64_t(i), maxSkipForBackbone);
            for(uint32_t step = 1; step <= lookBack; step++) {
                const uint32_t j = i - step;
                const Shasta2AnchorId anchorJ = journey[seedBegin + j];
                if(shasta2Anchors->countCommon(anchorJ, anchorI) >= minCommonForBackbone) {
                    if(dp[j] + 1 > dp[i]) {
                        dp[i] = dp[j] + 1;
                        prev[i] = int32_t(j);
                    }
                }
            }
        }

        // Find the position with the longest subsequence.
        uint32_t bestEnd = 0;
        for(uint32_t i = 1; i < n; i++) {
            if(dp[i] > dp[bestEnd]) {
                bestEnd = i;
            }
        }

        // Reconstruct the subsequence.
        vector<uint32_t> reversePath;
        for(int32_t idx = int32_t(bestEnd); idx >= 0; idx = prev[idx]) {
            reversePath.push_back(seedBegin + uint32_t(idx));
        }
        filteredPositions.assign(reversePath.rbegin(), reversePath.rend());
    };

    auto createWindow = [&](OrientedReadId backboneOid, uint32_t seedBegin, uint32_t seedEnd) {
        // Filter the backbone journey to keep only the longest subsequence
        // where consecutive pairs have sufficient common read support.
        vector<uint32_t> filteredPositions;
        filterBackboneJourney(backboneOid, seedBegin, seedEnd, filteredPositions);

        if(filteredPositions.size() < minBackboneWindowAnchors) {
            return; // Too few anchors after filtering.
        }

        const uint32_t windowId = uint32_t(anchorWindows.size());
        AnchorWindow window;
        window.windowId = windowId;
        window.backboneOrientedReadId = backboneOid;
        // Use the filtered positions: first and last define the span.
        window.backboneBegin = filteredPositions.front();
        window.backboneEnd = filteredPositions.back() + 1;
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

            // Step 3: LIS of backbone positions to enforce backbone order.
            const auto lisIndices = longestIncreasingSubsequence(sharedBackbonePositions);

            if(lisIndices.empty()) continue;

            // Step 4: Record the span and the read interval.
            uint32_t convergentCount = uint32_t(lisIndices.size());
            uint32_t convergentBegin = sharedReadPositions[lisIndices.front()];
            uint32_t convergentEnd = sharedReadPositions[lisIndices.back()] + 1;

            readSpans.push_back(ReadSpan{oid, convergentBegin, convergentEnd});

            window.readIntervals.push_back(AnchorWindowReadInterval{
                oid,
                convergentBegin,
                convergentEnd,
                convergentCount});

            // Step 5: For non-direct overlaps, extract alternate paths.
            // Between consecutive LIS pillars, collect intermediate anchors
            // from the read's journey, but only keep those that have at least
            // one oriented read that is not a direct cis overlap of the backbone.
            // This filters out anchors that merely echo the backbone haplotype.
            if(!isDirectCis && lisIndices.size() >= 2) {
                for(uint64_t li = 0; li + 1 < lisIndices.size(); li++) {
                    const uint32_t readPosA = sharedReadPositions[lisIndices[li]];
                    const uint32_t readPosB = sharedReadPositions[lisIndices[li + 1]];
                    if(readPosB > readPosA + 1) {
                        AnchorWindowAlternatePath altPath;
                        altPath.anchorIdA = journey[readPosA];
                        altPath.anchorIdB = journey[readPosB];
                        for(uint32_t rp = readPosA + 1; rp < readPosB; rp++) {
                            const Shasta2AnchorId midAnchorId = journey[rp];
                            // Skip anchors with duplicate k-mers in this read.
                            if(readDuplicates.count(uint64_t(midAnchorId))) continue;
                            ++totalIntermediatesConsidered;
                            // Keep only anchors that have at least one read
                            // that is not the backbone and not a direct cis overlap.
                            const auto anchor = (*shasta2Anchors)[midAnchorId];
                            bool hasNonCisRead = false;
                            for(const Shasta2AnchorMarkerInfo& ami : anchor) {
                                const uint32_t readOidValue = ami.orientedReadId.getValue();
                                if(readOidValue == backboneOid.getValue()) continue;
                                if(directCisOverlapReads.find(readOidValue) == directCisOverlapReads.end()) {
                                    hasNonCisRead = true;
                                    break;
                                }
                            }
                            if(hasNonCisRead) {
                                altPath.intermediateAnchorIds.push_back(midAnchorId);
                                ++totalIntermediatesKept;
                            }
                        }
                        ++totalAltPathsBeforeFilter;
                        if(!altPath.intermediateAnchorIds.empty()) {
                            window.alternatePaths.push_back(std::move(altPath));
                            ++totalAltPathsAfterFilter;
                        }
                    }
                }
            }
        }

        // TODO: revisit alternate path filtering:
        // 1. The direct-cis filter removes intermediates that only have
        //    direct-cis reads. This may be too aggressive or too lenient
        //    depending on the phasing context.
        // 2. The deduplication assigns each intermediate to one path
        //    (furthest pillar B). This is arbitrary and may discard
        //    valid alternate paths.
        // 3. Paths with no surviving intermediates are dropped entirely.
        // These heuristics should be revisited once phasing-aware
        // alternate paths are working.

        // Deduplicate alternate path intermediates: each intermediate anchor
        // must appear in exactly one alternate path. If the same intermediate
        // appears in multiple paths, assign it to the path with the furthest
        // pillar B (highest backbone position). This maximizes forward
        // connectivity for long-range phasing decisions.
        {
            auto pillarBPos = [&](const AnchorWindowAlternatePath& p) -> uint32_t {
                auto it = backboneAnchorToPos.find(uint64_t(p.anchorIdB));
                return (it != backboneAnchorToPos.end()) ? it->second : 0;
            };

            // Sort paths by pillar B position (furthest first).
            std::sort(window.alternatePaths.begin(), window.alternatePaths.end(),
                [&](const AnchorWindowAlternatePath& a, const AnchorWindowAlternatePath& b) {
                    return pillarBPos(a) > pillarBPos(b);
                });

            // Assign each intermediate to the first (tightest) path that contains it.
            std::unordered_set<uint64_t> claimedIntermediates;
            vector<AnchorWindowAlternatePath> filteredPaths;
            for(AnchorWindowAlternatePath& altPath : window.alternatePaths) {
                vector<Shasta2AnchorId> filtered;
                for(const Shasta2AnchorId mid : altPath.intermediateAnchorIds) {
                    if(claimedIntermediates.insert(uint64_t(mid)).second) {
                        filtered.push_back(mid);
                    }
                }
                if(!filtered.empty()) {
                    altPath.intermediateAnchorIds = std::move(filtered);
                    filteredPaths.push_back(std::move(altPath));
                }
            }
            window.alternatePaths = std::move(filteredPaths);
        }

        // Claim all anchors across all spans (backbone + touched reads) and their RC.
        for(const ReadSpan& span : readSpans) {
            const auto journey = (*shasta2Journeys)[span.oid];
            for(uint32_t readPos = span.begin; readPos < span.end; readPos++) {
                const uint64_t anchorId = uint64_t(journey[readPos]);
                if(anchorOwner[anchorId] == anchorUnclaimed) {
                    claimAnchor(anchorId, windowId);
                    ++window.claimedAnchorCount;
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

    // Process the heap.
    while(!candidateHeap.empty()) {
        const CleanWindowCandidate candidate = candidateHeap.top();
        candidateHeap.pop();

        const ReadId readId = candidate.backboneOrientedReadId.getReadId();
        if(candidate.generation != candidateGeneration[uint64_t(readId)]) {
            continue;
        }
        const auto journey = (*shasta2Journeys)[candidate.backboneOrientedReadId];

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
    cout << timestamp << "Alternate path filter:"
         << " pathsBefore=" << totalAltPathsBeforeFilter
         << " pathsAfter=" << totalAltPathsAfterFilter
         << " intermediatesConsidered=" << totalIntermediatesConsidered
         << " intermediatesKept=" << totalIntermediatesKept
         << endl;
}
