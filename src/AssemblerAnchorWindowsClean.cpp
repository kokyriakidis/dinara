
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
    uint64_t anchorCount;
    uint64_t readLength;
    uint32_t generation;
};

// Priority key for the seeding heap. When false (default), candidates are popped
// largest base span first (the original behavior). When true, they are popped
// largest anchor count first, so windows are seeded from the interval covering
// the most anchors rather than the most base pairs. Base span and anchor count
// swap primary/secondary roles so each remains the tiebreak of the other.
constexpr bool orderByAnchorCount = true;

struct CleanWindowCandidateLess {
    bool operator()(const CleanWindowCandidate& a, const CleanWindowCandidate& b) const {
        if(orderByAnchorCount) {
            if(a.anchorCount != b.anchorCount) return a.anchorCount < b.anchorCount;
            if(a.baseSpan != b.baseSpan) return a.baseSpan < b.baseSpan;
        } else {
            if(a.baseSpan != b.baseSpan) return a.baseSpan < b.baseSpan;
            if(a.anchorCount != b.anchorCount) return a.anchorCount < b.anchorCount;
        }
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

// Window seeding now runs in two priority-queue passes (see runPqPass in
// computeAnchorWindowsClean): a first pass that seeds pristine full-journey
// cores only (one window = one whole read journey), and an optional second
// pass (tileUnclaimedIntervals) that tiles the leftover unclaimed runs into
// fragment windows, longest base span first. The second pass is off by default,
// so the single-pass full-journey-only behavior is the default.

// Full-journey-only window seeding (legacy toggle). When true, a read seeds a
// window ONLY if its entire journey is still unclaimed (one window = one whole
// read journey); partially-claimed reads are skipped, leaving the dovetail/
// overlap seams between cores unclaimed. When false, the original behavior was
// used: each contiguous unclaimed run of a read's journey seeds its own window.
//
// Superseded by the two-pass scheme above (pass 1 == this flag true, pass 2 ==
// tiling the seams). Kept (currently unused) so the single-pass gate can be
// restored quickly while window construction is still being experimented with.
[[maybe_unused]] constexpr bool fullJourneyWindowsOnly = true;

} // anonymous namespace



void Assembler::computeAnchorWindowsClean(
    shared_ptr<Shasta2Anchors> shasta2Anchors,
    shared_ptr<Shasta2Journeys> shasta2Journeys,
    const vector<ReadId>& readIdsSortedByLength,
    vector<AnchorWindow>& anchorWindows,
    uint64_t /* threadCount */,
    uint64_t minCommonForBackbone,
    uint64_t maxSkipForBackbone,
    uint64_t minWindowBaseSpan,
    vector<uint32_t>* anchorDovetailWindow,
    bool tileUnclaimedIntervals)
{
    // threadCount is intentionally unused: window claiming is a single
    // priority-queue sweep over shared, mutually-exclusive anchor ownership
    // (anchorOwner), with lazy stale-candidate invalidation
    // (candidateGeneration) driving a strict largest-first claim order.
    // That's inherently serial, not just not-yet-parallelized -- two threads
    // racing to claim overlapping journey spans would break both the
    // disjoint-windows invariant and the deterministic tiebreak order the
    // priority comparator (CleanWindowCandidateLess) is built around. Kept
    // in the signature for parity with sibling mode3-pipeline calls.
    cout << timestamp << "computeAnchorWindowsClean begins." << endl;
    const auto t0 = steady_clock::now();

    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    DINARA_ASSERT(shasta2Anchors);
    DINARA_ASSERT(shasta2Journeys);
    DINARA_ASSERT(shasta2Journeys->isOpen());
    DINARA_ASSERT(reads->readCount() > 0);

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

    // Temporary phase-timing breakdown (env DINARA_WINDOW_CLAIM_PROFILE=1), to
    // find where computeAnchorWindowsClean's wall time actually goes. Overhead
    // is a few steady_clock::now() calls per createWindow invocation (negligible
    // next to the work being timed).
    const bool profileClaiming = (getenv("DINARA_WINDOW_CLAIM_PROFILE") != nullptr);
    double profTouchedDiscovery = 0.0, profDirectCis = 0.0, profBackboneDup = 0.0;
    double profSharedScan = 0.0, profLis = 0.0, profClaim = 0.0, profRepush = 0.0;
    double profInitSeed = 0.0, profUnclaimedScan = 0.0, profCreateWindowTotal = 0.0, profReseed = 0.0;
    double profFilteredPosBuild = 0.0, profBackboneMapBuild = 0.0, profDupCacheLookup = 0.0;
    double profDupCacheHitTime = 0.0, profDupCacheMissTime = 0.0;
    uint64_t profDupCacheHits = 0, profDupCacheMisses = 0;
    double profFdkaBuildKmers = 0.0, profFdkaSeenSet = 0.0, profFdkaMarkDup = 0.0;

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
            uint64_t(end - begin),
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
    // LIS helper: given a sequence of values (uint32_t), return the indices
    // (into the input) of a strictly-increasing longest subsequence.
    // Uses patience sorting (O(n log n)).
    // ========================================================================
    auto longestIncreasingSubsequence = [](const vector<uint32_t>& values) {
        vector<uint32_t> result;
        const uint32_t n = uint32_t(values.size());
        if(n == 0) return result;

        // tailsIdx[len-1] = index (into values) of the smallest tail value of
        // any strictly-increasing subsequence of length len found so far.
        // predecessor[i] = index of the element before values[i] in the LIS
        // that ends at i.
        vector<uint32_t> tailsIdx;
        vector<int32_t> predecessor(n, -1);

        for(uint32_t i = 0; i < n; i++) {
            // Binary search for the first tail whose value is >= values[i]
            // (strict increase: replace that tail, extend otherwise).
            uint32_t lo = 0;
            uint32_t hi = uint32_t(tailsIdx.size());
            while(lo < hi) {
                const uint32_t mid = (lo + hi) / 2;
                if(values[tailsIdx[mid]] < values[i]) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            if(lo > 0) {
                predecessor[i] = int32_t(tailsIdx[lo - 1]);
            }
            if(lo == uint32_t(tailsIdx.size())) {
                tailsIdx.push_back(i);
            } else {
                tailsIdx[lo] = i;
            }
        }

        // Reconstruct from the last tail backwards.
        for(int32_t k = int32_t(tailsIdx.back()); k >= 0; k = predecessor[uint32_t(k)]) {
            result.push_back(uint32_t(k));
        }
        std::reverse(result.begin(), result.end());
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
    // For each touched read, find first/last shared anchor with backbone
    // and claim all anchors in between.
    // ========================================================================
    // KmerId lookup for an anchor's representative marker (front() of its
    // sorted member list), reusing the marker's OWN precomputed kmerId
    // (populated pipeline-wide, once, by computeMarkerKmerIds -- called well
    // before window creation, for both the SIMD and non-SIMD marker paths)
    // instead of reconstructing the k-mer base by base from the read sequence
    // via Shasta2Anchors::anchorKmer(). Valid only for PRIMARY anchors, whose
    // Shasta2AnchorMarkerInfo::ordinal is always a valid index into markers /
    // markerKmerIds (het/hom anchors leave ordinal invalid -- see
    // Shasta2AnchorMarkerInfo -- but window creation runs before any het/hom
    // anchor is appended, so every anchor here is primary).
    //
    // markerKmerIds[oid][ordinal] is populated via the same extraction and
    // strand handling (map strand-1 to strand-0's raw sequence, extract,
    // reverse-complement) that Shasta2Anchors::getKmerAtPosition/anchorKmer
    // perform -- see getOrientedReadMarkerKmerStrand0/1 in
    // markerAccessFunctions.cpp -- so this is bit-for-bit the same value
    // anchorKmer() would return for the same (orientedReadId, ordinal), just
    // an O(1) lookup instead of an O(k) base-by-base reconstruction.
    // Profiling showed anchorKmer() was ~59% of window creation's remaining
    // wall time after the caching fixes above.
    auto anchorKmerId = [&](Shasta2AnchorId anchorId) -> KmerId {
        const Shasta2Anchor anchor = (*shasta2Anchors)[anchorId];
        if(anchor.empty()) {
            return KmerId(0);
        }
        const Shasta2AnchorMarkerInfo& mi = anchor.front();
        DINARA_ASSERT(mi.ordinal != invalid<uint32_t>);
        return (*markerKmerIds)[mi.orientedReadId.getValue()][mi.ordinal];
    };

    // Find anchor IDs in a journey interval whose k-mer appears more than once.
    // Anchors at the first and last positions of the interval are never marked
    // as duplicates, to avoid disconnecting chain endpoints.
    auto findDuplicateKmerAnchors = [&](const auto& journey, uint32_t begin, uint32_t end) {
        std::unordered_set<uint64_t> duplicateAnchorIds;
        if(end <= begin) {
            return duplicateAnchorIds;
        }
        const auto profFdkaT0 = steady_clock::now();
        vector<KmerId> kmers(end - begin);
        for(uint32_t pos = begin; pos < end; pos++) {
            kmers[pos - begin] = anchorKmerId(journey[pos]);
        }
        if(profileClaiming) profFdkaBuildKmers += seconds(steady_clock::now() - profFdkaT0);

        const auto profFdkaT1 = steady_clock::now();
        std::unordered_set<KmerId> seen;
        std::unordered_set<KmerId> duplicateKmers;
        seen.reserve(end - begin);
        for(uint32_t pos = begin; pos < end; pos++) {
            if(!seen.insert(kmers[pos - begin]).second) {
                duplicateKmers.insert(kmers[pos - begin]);
            }
        }
        if(profileClaiming) profFdkaSeenSet += seconds(steady_clock::now() - profFdkaT1);

        const auto profFdkaT2 = steady_clock::now();
        if(!duplicateKmers.empty()) {
            for(uint32_t pos = begin; pos < end; pos++) {
                if(pos == begin || pos == end - 1) continue;
                if(duplicateKmers.count(kmers[pos - begin])) {
                    duplicateAnchorIds.insert(uint64_t(journey[pos]));
                }
            }
        }
        if(profileClaiming) profFdkaMarkDup += seconds(steady_clock::now() - profFdkaT2);
        return duplicateAnchorIds;
    };

    // Cache of findDuplicateKmerAnchors(journey, 0, journey.size()) for a read's
    // FULL own journey, keyed by oidValue. This depends only on the read's own
    // (fixed, read-only during window creation) journey, not on which window is
    // asking, but the same oriented read can be "touched" by many different
    // windows over the course of this sweep (a read stays touchable until every
    // anchor it shares with some backbone is claimed) -- each touch previously
    // recomputed this from scratch by rescanning the read's whole journey with a
    // std::set. Sparse map: only reads actually touched pay for an entry, so
    // memory stays bounded by the number of distinct reads ever touched (no
    // worse than before), while repeat touches of the same read become O(1).
    std::unordered_map<uint32_t, std::unordered_set<uint64_t>> readDuplicatesCache;
    auto getFullJourneyDuplicateKmerAnchors = [&](OrientedReadId oid, const auto& journey)
        -> const std::unordered_set<uint64_t>&
    {
        const uint32_t oidValue = uint32_t(oid.getValue());
        if(profileClaiming) {
            const auto profT0hit = steady_clock::now();
            const auto existing = readDuplicatesCache.find(oidValue);
            if(existing != readDuplicatesCache.end()) {
                profDupCacheHitTime += seconds(steady_clock::now() - profT0hit);
                ++profDupCacheHits;
                return existing->second;
            }
            const auto& result = readDuplicatesCache.emplace(oidValue,
                findDuplicateKmerAnchors(journey, 0, uint32_t(journey.size()))).first->second;
            profDupCacheMissTime += seconds(steady_clock::now() - profT0hit);
            ++profDupCacheMisses;
            return result;
        }
        const auto existing = readDuplicatesCache.find(oidValue);
        if(existing != readDuplicatesCache.end()) {
            return existing->second;
        }
        return readDuplicatesCache.emplace(oidValue,
            findDuplicateKmerAnchors(journey, 0, uint32_t(journey.size()))).first->second;
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

        const auto profT0start = steady_clock::now();
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
        if(profileClaiming) profFilteredPosBuild += seconds(steady_clock::now() - profT0start);

        // Get backbone journey and build a lookup: anchorId -> position in backbone.
        // Only include anchors that survived filtering and don't have duplicate k-mers.
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];
        std::unordered_map<uint64_t, uint32_t> backboneAnchorToPos;
        backboneAnchorToPos.reserve(filteredPositions.size());

        // Find backbone anchors with duplicate k-mers and exclude them.
        const auto profT0a = steady_clock::now();
        const auto backboneDuplicates = findDuplicateKmerAnchors(
            backboneJourney, window.backboneBegin, window.backboneEnd);
        if(profileClaiming) profBackboneDup += seconds(steady_clock::now() - profT0a);
        // Build the lookup from filtered positions only.
        const auto profT0map = steady_clock::now();
        for(const uint32_t pos : filteredPositions) {
            const uint64_t aid = uint64_t(backboneJourney[pos]);
            if(backboneDuplicates.count(aid) == 0) {
                backboneAnchorToPos[aid] = pos;
            }
        }
        if(profileClaiming) profBackboneMapBuild += seconds(steady_clock::now() - profT0map);

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
        const auto profT0b = steady_clock::now();
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
        if(profileClaiming) profTouchedDiscovery += seconds(steady_clock::now() - profT0b);

        // Build the set of direct cis overlap oriented read IDs for the backbone
        // using the read graph (which contains only cis overlaps).
        const auto profT0c = steady_clock::now();
        std::unordered_set<uint32_t> directCisOverlapReads;
        {
            const auto backboneEdges = readGraph.connectivity[backboneOid.getValue()];
            for(const uint32_t edgeId : backboneEdges) {
                const auto& edge = readGraph.edges[edgeId];
                const OrientedReadId partner = edge.getOther(backboneOid);
                directCisOverlapReads.insert(partner.getValue());
            }
        }
        if(profileClaiming) profDirectCis += seconds(steady_clock::now() - profT0c);

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
            const auto profT0dup = steady_clock::now();
            const auto& readDuplicates = getFullJourneyDuplicateKmerAnchors(oid, journey);
            if(profileClaiming) profDupCacheLookup += seconds(steady_clock::now() - profT0dup);

            // Step 1-2: Find shared anchors and their backbone positions.
            vector<uint32_t> sharedReadPositions;
            vector<uint32_t> sharedBackbonePositions;

            const auto profT0d = steady_clock::now();
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
            if(profileClaiming) profSharedScan += seconds(steady_clock::now() - profT0d);

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

            AnchorWindowReadInterval readInterval{
                oid,
                convergentBegin,
                convergentEnd,
                convergentCount};

            // Persist the shared-anchor pins (read/backbone journey positions),
            // ordered by backbone position so consecutive pins delimit the
            // inter-anchor segments used by the per-segment abPOA MSA.
            //
            // sharedReadPositions is strictly increasing (collected in read
            // order). The shared anchors are NOT guaranteed colinear with the
            // backbone: a repeat can pair an anchor out of diagonal order, so
            // sharedBackbonePositions may dip. The per-segment MSA requires
            // pins where both read and backbone positions increase together
            // (a proper alignment's diagonal). Keep the longest such colinear
            // subset via LIS over the backbone positions (read order is the
            // implicit increasing axis), discarding off-diagonal pins.
            {
                const auto profT0e = steady_clock::now();
                const vector<uint32_t> keep =
                    longestIncreasingSubsequence(sharedBackbonePositions);
                vector<AnchorWindowSharedPin>& pins = readInterval.sharedPins;
                pins.reserve(keep.size());
                for(const uint32_t s : keep) {
                    pins.push_back(AnchorWindowSharedPin{
                        sharedReadPositions[s], sharedBackbonePositions[s]});
                }
                // LIS over backbone positions (read order increasing) yields
                // pins already sorted by backbone position with strictly
                // increasing read positions, matching the consumer's asserts.
                if(profileClaiming) profLis += seconds(steady_clock::now() - profT0e);
            }

            window.readIntervals.push_back(std::move(readInterval));

        }

        // Claim all anchors across all spans (backbone + touched reads) and their RC.
        const auto profT0f = steady_clock::now();
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

        if(profileClaiming) profClaim += seconds(steady_clock::now() - profT0f);

        // Now that claiming is done, bump generations and re-push unclaimed
        // intervals for all touched reads.
        const auto profT0g = steady_clock::now();
        for(const uint32_t oidValue : touchedOrientedReads) {
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            ++candidateGeneration[uint64_t(oid.getReadId())];
            pushCurrentUnclaimedIntervals(oid);
        }
        if(profileClaiming) profRepush += seconds(steady_clock::now() - profT0g);

        claimedAnchors += window.claimedAnchorCount;
        totalReadIntervals += window.readIntervals.size();
        anchorWindows.push_back(std::move(window));
    };

    // Initialize the heap: all strand-0 reads, ordered by length (longest first).
    {
        const auto profT0h = steady_clock::now();
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
        if(profileClaiming) profInitSeed += seconds(steady_clock::now() - profT0h);
    }

    // Process the heap. For each candidate (longest base span first):
    // - If the full journey is unclaimed, create one window for the whole journey.
    // - Otherwise (only when fragments are allowed), find contiguous unclaimed
    //   intervals within the journey and create a separate window for each.
    //
    // The loop body is factored into runPqPass so it can be run twice: a first
    // pass that seeds pristine full-journey cores only, and an optional second
    // pass (see computeAnchorWindowsUnclaimed) that tiles the leftover
    // unclaimed intervals into fragment windows, again longest-span first.
    uint64_t fullJourneyWindows = 0;
    uint64_t fragmentWindows = 0;

    auto runPqPass = [&](bool allowFragments) {
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
            const auto profT0i = steady_clock::now();
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
            if(profileClaiming) profUnclaimedScan += seconds(steady_clock::now() - profT0i);

            if(unclaimedIntervals.empty()) continue;

            // Create a window for each unclaimed interval that passes the
            // base span threshold.
            const bool isFullJourney = (unclaimedIntervals.size() == 1 &&
                unclaimedIntervals[0].first == 0 &&
                unclaimedIntervals[0].second == uint32_t(journey.size()));

            // First pass (allowFragments == false): only seed reads whose entire
            // journey is still unclaimed, producing pristine disjoint cores. The
            // second pass (allowFragments == true) tiles the leftover unclaimed
            // runs into fragment windows.
            if(!allowFragments && !isFullJourney) continue;

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

                const auto profT0j = steady_clock::now();
                createWindow(candidate.backboneOrientedReadId, intervalBegin, intervalEnd);
                if(profileClaiming) profCreateWindowTotal += seconds(steady_clock::now() - profT0j);

                if(isFullJourney) {
                    ++fullJourneyWindows;
                } else {
                    ++fragmentWindows;
                }
            }
        }
    };

    // Pure priority-queue mode (legacy/original behavior). When true, run a
    // single PQ pass with fragments allowed from the start: each contiguous
    // unclaimed run of a read's journey seeds its own window, longest base span
    // first, with no preceding full-journey-cores pass. Experimental toggle;
    // when false, the two-pass scheme (cores, then optional tiling) runs.
    //
    // Disabled (false): run the two-pass scheme so windows are seeded only from
    // full read journeys. Pass 1 (runPqPass(false)) seeds a window only from a
    // read whose entire journey is still unclaimed, producing pristine disjoint
    // cores. Pass 2 (fragment tiling) runs only if tileUnclaimedIntervals is
    // set by the caller; with it off, reads whose journey overlaps an already-
    // claimed core contribute no window (their anchors stay unclaimed). This is
    // the "windows from full journeys only" behavior.
    constexpr bool purePriorityQueue = false;
    if(purePriorityQueue) {
        runPqPass(true);
    } else {
        // First pass: pristine full-journey cores only.
        runPqPass(false);

        // Optional second pass: tile the leftover unclaimed intervals into
        // fragment windows, longest-span first. Re-seed the heap from every
        // read's current unclaimed runs (the first pass drained it), then run
        // with fragments allowed. Disabled by default so the single-pass
        // behavior is unchanged.
        if(tileUnclaimedIntervals) {
            const uint64_t windowsBefore = anchorWindows.size();
            const auto profT0k = steady_clock::now();
            for(const ReadId readId : readIdsSortedByLength) {
                const OrientedReadId backboneOid(readId, 0);
                if(backboneOid.getValue() >= shasta2Journeys->size()) continue;
                const auto journey = (*shasta2Journeys)[backboneOid];
                if(journey.empty()) continue;
                ++candidateGeneration[uint64_t(readId)];  // invalidate stale.
                pushCurrentUnclaimedIntervals(backboneOid);
            }
            if(profileClaiming) profReseed += seconds(steady_clock::now() - profT0k);
            runPqPass(true);
            cout << timestamp << "Unclaimed-interval pass added "
                 << (anchorWindows.size() - windowsBefore)
                 << " fragment window(s)." << endl;
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

    if(profileClaiming) {
        cout << std::fixed << std::setprecision(2)
             << "  [claim-profile] initSeed=" << profInitSeed
             << "s unclaimedScan=" << profUnclaimedScan
             << "s reseed=" << profReseed
             << "s createWindowTotal=" << profCreateWindowTotal
             << "s (of which filteredPosBuild=" << profFilteredPosBuild
             << "s backboneMapBuild=" << profBackboneMapBuild
             << "s backboneDup=" << profBackboneDup
             << "s touchedDiscovery=" << profTouchedDiscovery
             << "s directCis=" << profDirectCis
             << "s dupCacheLookup=" << profDupCacheLookup
             << "s [hits=" << profDupCacheHits << " hitTime=" << profDupCacheHitTime
             << "s misses=" << profDupCacheMisses << " missTime=" << profDupCacheMissTime
             << "s avgMissUs=" << (profDupCacheMisses ? (profDupCacheMissTime * 1e6 / double(profDupCacheMisses)) : 0.0)
             << " buildKmers=" << profFdkaBuildKmers
             << "s seenSet=" << profFdkaSeenSet
             << "s markDup=" << profFdkaMarkDup
             << "s]"
             << "s sharedScan=" << profSharedScan
             << "s lis=" << profLis
             << "s claim=" << profClaim
             << "s repush=" << profRepush << "s)"
             << " sum=" << (profInitSeed + profUnclaimedScan + profReseed + profCreateWindowTotal)
             << "s / total=" << elapsedSeconds << "s"
             << std::defaultfloat << endl;
    }
}



// Convenience wrapper: disjoint full-journey cores, then a second pass that
// tiles the leftover unclaimed intervals into fragment windows (longest first).
void Assembler::computeAnchorWindowsWithUnclaimed(
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
    computeAnchorWindowsClean(
        shasta2Anchors,
        shasta2Journeys,
        readIdsSortedByLength,
        anchorWindows,
        threadCount,
        minCommonForBackbone,
        maxSkipForBackbone,
        minWindowBaseSpan,
        anchorDovetailWindow,
        /* tileUnclaimedIntervals = */ true);
}
