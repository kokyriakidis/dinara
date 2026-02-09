/**
 * @file AssemblerInvertedIndex.cpp
 * @brief High-performance alignment candidate discovery using an Inverted Index.
 *
 * This file implements a Hifiasm-compatible chaining algorithm for finding
 * potential read overlaps. The core data flow is:
 *
 *   1. **Inverted Index Construction**: Build a hash table mapping each
 *      canonical K-mer to a list of (ReadId, Position) occurrences.
 *   2. **Parallel Radix Sort**: Sort occurrences by K-mer ID using an O(N)
 *      LSD radix sort to group identical K-mers.
 *   3. **Hash Table Build**: Populate a power-of-2 sized hash table for O(1)
 *      K-mer lookups during the search phase.
 *   4. **Parallel Candidate Search**: For each read, query the index to find
 *      matching K-mers in other reads. Use DP chaining to score overlaps.
 *
 * Key performance features:
 *   - Structure-of-Arrays (SoA) layout for cache-efficient DP.
 *   - Early K-mer weighting based on occurrence frequency (Hifiasm-compatible).
 *   - Prefix quick-check to skip quadratic DP on monotonic hit runs.
 *   - Shared mcopy-fast style endpoint extraction for repetitive regions.
 *
 * @note All scoring and tie-breaking rules are strictly Hifiasm-compatible.
 */

#include "Assembler.hpp"
#include "hifiasmCoordinateTransforms.hpp"
#include "performanceLog.hpp"
#include "OrientedReadPair.hpp"
#include "timestamp.hpp"
#include "Reads.hpp"
#include <algorithm>
#include <barrier>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include <thread>

using namespace std;

#include "MultithreadedObject.tpp"
#include "Alignment.hpp"
#include "AlignmentCanonicalization.hpp"

namespace dinara {

/**
 * @brief Temporary hit structure used during K-mer matching.
 *
 * Represents a single K-mer match between the query read (Read A) and a
 * partner read. Hits are first collected in Array-of-Structures (AoS) format
 * for sorting, then converted to Structure-of-Arrays (SoA) for DP.
 *
 * @note Sorting is by (partnerReadId, posA) to group hits by read pair and
 *       ensure monotonically increasing query positions for the DP.
 */
struct InvertedIndexTempHit {
    ReadId partnerReadId;  ///< ID of the matching read (Read B).
    uint32_t posA;         ///< Base position of the K-mer in Read A.
    uint32_t posB;         ///< Base position of the K-mer in Read B.
    uint32_t ordinalA;     ///< Marker ordinal in Read A (for alignment output).
    uint8_t weight;        ///< Frequency-based weight (Hifiasm compatible).
    
    /// Comparison operator for sorting by (partnerReadId, posA).
    bool operator<(const InvertedIndexTempHit& other) const {
        if (partnerReadId != other.partnerReadId) return partnerReadId < other.partnerReadId;
        return posA < other.posA;
    }
};

/**
 * @brief Classifies an overlap based on Hifiasm's overlap type heuristic.
 *
 * This is used for Hifiasm's candidate limit feature, which caps the number
 * of reported overlaps per type to avoid excessive output in repetitive regions.
 *
 * @param start  Starting coordinate of the alignment on Read A.
 * @param end    Ending coordinate of the alignment on Read A.
 * @param readLen Total length of Read A.
 * @return 0 = Left overhang, 1 = Right overhang, 2 = Contained, 3 = Containing.
 */
// ============================================================================
// OVERLAP TYPE CLASSIFICATION (Hifiasm COV_W grouping)
// ============================================================================
// Reference: Hifiasm anchor.cpp overlap type classification for COV_W control
//
// Classifies overlaps into 4 types based on how the alignment spans the read:
// - Type 0: Left overhang  (alignment starts at position 0, ends before read end)
// - Type 1: Right overhang (alignment starts after position 0, ends at read end)
// - Type 2: Contained      (alignment spans entire read, start=0 and end=readLen-1)
// - Type 3: Containing     (alignment is internal, doesn't reach either end)
//
// Purpose: Different overlap types have different characteristics:
// - Type 2 (contained): Read is completely covered, may need special handling
// - Type 3 (containing): Internal alignment, subject to COV_W window control
// - Types 0/1 (overhangs): Typical overlap extensions
//
// This classification is used for:
// 1. Per-type max_n_chain filtering (different thresholds per type)
// 2. COV_W (Coverage Window) overload control for type 3 overlaps
// ============================================================================
static int getOverlapType(uint32_t start, uint32_t end, uint32_t readLen) {
    if (start == 0 && end >= readLen - 1) return 2; // Contained
    else if (start > 0 && end < readLen - 1) return 3; // Containing
    else return (start == 0) ? 0 : 1; // Left or Right overhang
}

// Hifiasm weak-overlap suppression constants (anchor.cpp).
static constexpr double HIFIASM_OFL = 0.95;
static constexpr uint32_t HIFIASM_CH_OCC = 4;
static constexpr uint32_t HIFIASM_CH_SC = 16;

struct HifiasmLchainDpOptions {
    int32_t maxSkip = 25;
    int32_t maxIter = 5000;
    int32_t maxDist = 5000;
    double chnPenGap = 0.0;
    double chnPenSkip = 0.0;
    bool quickCheck = false;
};

// Hifiasm parity for anchor.cpp:set_lchain_dp_op.
// ============================================================================
// HIFIASM CHAINING PARAMETER CONFIGURATION
// ============================================================================
// Configure dynamic programming parameters for the chaining algorithm.
// Reference: Hifiasm anchor.cpp:2272 set_lchain_dp_op()
//
// Parameters control the DP search space and penalty structure:
// - quickCheck: Enable fast prefix checking optimization (1 for ONT)
// - chnPenGap: Linear penalty coefficient for indels/gaps
// - chnPenSkip: Penalty for skipping bases between consecutive anchors
//
// The penalties use exponential decay based on k-mer size:
//   penalty = base * exp(-div * k)
// Where:
//   - div = 0.01 for ONT (accurate), 0.1 for HiFi
//   - Larger k → more specific → lower penalty (trust gaps more)
//   - Smaller k → less specific → higher penalty (penalize gaps more)
//
static inline HifiasmLchainDpOptions getHifiasmLchainDpOptions(
    const bool isAccurate,        // true=ONT, false=HiFi
    const uint32_t markerK)       // K-mer size
{
    HifiasmLchainDpOptions o;
    o.quickCheck = isAccurate;    // Enable fast path for ONT

    // Exponential decay factor: ONT uses slower decay (0.01 vs 0.1)
    const double div = isAccurate ? 0.01 : 0.1;

    // Base penalties before k-mer adjustment
    const double penGap = 0.5;      // Gap penalty base
    const double penSkip = 0.0005;  // Skip penalty base

    // Apply exponential decay based on k-mer size
    // Example for k=34, ONT: exp(-0.01*34) ≈ 0.712
    const double tmp = std::exp(-div * double(markerK));
    o.chnPenGap = penGap * tmp;     // Final gap penalty
    o.chnPenSkip = penSkip * tmp;   // Final skip penalty

    return o;
}

/**
 * @brief Per-thread scratchpad for high-performance DP chaining.
 *
 * This struct uses a Structure-of-Arrays (SoA) layout for the DP arrays,
 * which is critical for cache efficiency during the O(N^2) chaining loop.
 * The scratchpad is recycled across reads to avoid allocation overhead.
 *
 * Key arrays:
 *   - `dpSame`, `dpDiff`: DP scores for same-strand and opposite-strand chains.
 *   - `parentSame`, `parentDiff`: Backtrack pointers for chain reconstruction.
 *   - `backtrackVisit*`: marks nodes visited by descendants for max-skip pruning.
 *   - `chainOccurrences*`: chain anchor counts, used by mcopy/weak suppression.
 */
struct ThreadScratchpad {
    struct WeakFilterMeta {
        uint32_t qs = 0;
        uint32_t qe = 0;   // Half-open [qs, qe)
        uint32_t occ = 0;  // Number of anchors in the chain.
        int32_t score = 0;
    };

    // AoS for hit collection and sorting
    vector<InvertedIndexTempHit> flatHits;
    
    // Structure of Arrays (SoA) for cache-efficient DP scans
    vector<uint32_t> hitPosA, hitPosB, hitOrdinalA, hitOrdinalB;
    vector<uint32_t> hitOrderByPosB;
    vector<uint8_t> hitWeights;

    // DP score and backtrack arrays (int32_t since scores can be negative)
    vector<int32_t> dpSame, dpDiff;
    vector<int32_t> parentSame, parentDiff;
    vector<int32_t> backtrackVisitSame, backtrackVisitDiff;
    vector<uint32_t> chainOccurrencesSame, chainOccurrencesDiff;

    // Post-DP candidate extraction
    struct ChainCandidate {
        int32_t score;  // Changed to int32_t to match dpSame/dpDiff
        uint64_t chainLen; // For deterministic tie-breaking
        int32_t endK;
        bool isDiff;
        bool operator<(const ChainCandidate& other) const {
            if (score != other.score) return score > other.score;
            return chainLen < other.chainLen;
        }
    };
    vector<ChainCandidate> chainCandidates;
    vector<ChainCandidate> filteredCandidates;
    vector<WeakFilterMeta> weakMetas;
    vector<size_t> weakIdxWorkspace;
    vector<size_t> strongIdxWorkspace;
    vector<uint8_t> suppressWorkspace;
    vector<uint64_t> strongAnchorBegin;
    vector<uint64_t> strongAnchorEnd;
    vector<uint32_t> strongAnchorStartsFlat;
    vector<uint8_t> mcopyNodeUsed;
    vector<int32_t> mcopyPathNodes;
    vector<ChainCandidate> mcopySelectedCandidates;

    struct ChainInterval { uint32_t qs; uint32_t qe; };
    vector<ChainInterval> acceptedIntervalsSame;
    vector<ChainInterval> acceptedIntervalsDiff;
    vector<uint32_t> currentChainPath;

    void clear() {
        flatHits.clear();
        hitPosA.clear(); hitPosB.clear(); hitOrdinalA.clear(); hitOrdinalB.clear(); hitOrderByPosB.clear(); hitWeights.clear();
        dpSame.clear(); dpDiff.clear();
        parentSame.clear(); parentDiff.clear();
        backtrackVisitSame.clear(); backtrackVisitDiff.clear();
        chainOccurrencesSame.clear(); chainOccurrencesDiff.clear();
        chainCandidates.clear();
        filteredCandidates.clear();
        weakMetas.clear();
        weakIdxWorkspace.clear();
        strongIdxWorkspace.clear();
        suppressWorkspace.clear();
        strongAnchorBegin.clear();
        strongAnchorEnd.clear();
        strongAnchorStartsFlat.clear();
        mcopyNodeUsed.clear();
        mcopyPathNodes.clear();
        mcopySelectedCandidates.clear();
        acceptedIntervalsSame.clear();
        acceptedIntervalsDiff.clear();
        currentChainPath.clear();
    }
};

/**
 * @brief Fast reverse complement for K-mers using bit reversal.
 */
inline KmerId getRcKmerId(KmerId id, uint64_t k) {
    const KmerId mask = (KmerId(1) << k) - 1;
    const KmerId lsb = id & mask;
    const KmerId msb = (id >> k) & mask;
    
    auto reverseBits = [&](KmerId x) -> KmerId {
        if (sizeof(KmerId) <= 8) return bitReversal(uint64_t(x)) >> (64 - k);
        else return bitReversal((__uint128_t)x) >> (128 - k);
    };
    
    KmerId rc_lsb = (~reverseBits(lsb)) & mask;
    KmerId rc_msb = (~reverseBits(msb)) & mask;
    return (rc_msb << k) | rc_lsb;
}

/**
 * @brief Thread-safe hash function for K-mers.
 */
static inline uint64_t hashKmer(KmerId k) {
    const uint64_t* p = reinterpret_cast<const uint64_t*>(&k);
    uint64_t k1 = p[0];
    uint64_t k2 = sizeof(KmerId) > 8 ? p[1] : 0; 
    return k1 ^ (k2 + 0x9e3779b9 + (k1<<6) + (k1>>2));
}

struct QuickLinearChainResult {
    bool fullySolved = false;
    size_t solvedPrefix = 0;
    int32_t maxScore = 0;
    int32_t bestEnd = -1;
};

// ============================================================================
// QUICK LINEAR CHAINING OPTIMIZATION (Hifiasm quick_ck_lchain)
// ============================================================================
// Reference: Hifiasm anchor.cpp:1425-1470 (quick_ck_lchain)
//
// Purpose: Fast O(n) prefix chaining for strictly monotonic hit sequences
// - Many read pairs have perfectly collinear markers at the start
// - For these "clean" prefixes, we can skip the expensive O(n²) DP
// - Process hits one-by-one until monotonicity breaks
// - Return the solved prefix length for full DP to continue from
//
// Algorithm:
// 1. Start with first hit (i=0) as base
// 2. For each subsequent hit i:
//    a. Check if positions are strictly increasing on both reads:
//       - Same-strand: posA[i] > posA[i-1] AND posB[i] > posB[i-1]
//       - Diff-strand: posA[i] > posA[i-1] AND posB[i] < posB[i-1] (RC)
//    b. If monotonic: chain i → i-1 with score calculation
//    c. If NOT monotonic: STOP and return solved prefix
// 3. If entire sequence is monotonic: mark fullySolved=true
//
// Benefits:
// - Reduces DP from O(n²) to O(n) for clean alignments
// - Typical in high-quality data where many read pairs align perfectly
// - Saves significant CPU time on common cases
//
// Return:
// - solvedPrefix: Index where full DP should start
// - fullySolved: true if entire sequence was monotonic (skip full DP)
// - maxScore, bestEnd: Best chain endpoint found in prefix
// ============================================================================
static inline QuickLinearChainResult runQuickLinearChainPrefix(
    const bool isDiff,
    const bool useEcScoring,
    const vector<uint32_t>& hitPosA,
    const vector<uint32_t>& hitPosB,
    const vector<uint8_t>& hitWeights,
    const uint32_t kmerLen,
    const double bwRate,
    const double chnPenGap,
    const double chnPenSkip,
    const uint64_t readLenA,
    const uint64_t readLenB,
    vector<int32_t>& dp,
    vector<int32_t>& parent,
    vector<uint32_t>& chainOccurrences)
{
    QuickLinearChainResult result;
    const size_t n = hitPosA.size();
    if(n == 0) {
        return result;
    }

    const uint8_t span0 = uint8_t(std::min<uint32_t>(kmerLen, 255));
    dp[0] = int32_t(span0);
    parent[0] = -1;
    chainOccurrences[0] = 1;
    result.maxScore = dp[0];
    result.bestEnd = 0;
    result.solvedPrefix = 1;

    for(size_t i = 1; i < n; ++i) {
        if(hitPosA[i] <= hitPosA[i - 1]) {
            break;
        }

        const uint32_t posAi = hitPosA[i];
        const uint32_t posAj = hitPosA[i - 1];
        uint32_t posBi = 0;
        uint32_t posBj = 0;
        if(!isDiff) {
            if(hitPosB[i] <= hitPosB[i - 1]) {
                break;
            }
            posBi = hitPosB[i];
            posBj = hitPosB[i - 1];
        } else {
            if(hitPosB[i] >= hitPosB[i - 1]) {
                break;
            }
            posBi = uint32_t(readLenB - 1 - uint64_t(hitPosB[i]));
            posBj = uint32_t(readLenB - 1 - uint64_t(hitPosB[i - 1]));
        }

        const uint8_t spanI = uint8_t(std::min<uint32_t>(kmerLen, 255));
        int32_t sc = useEcScoring ?
            hifiasm_comput_sc_ch_ec(
                posAi,
                posBi,
                posAj,
                posBj,
                hitWeights[i],
                spanI,
                bwRate,
                chnPenGap,
                chnPenSkip,
                readLenA,
                readLenB) :
            hifiasm_comput_sc_ch(
                posAi,
                posBi,
                posAj,
                posBj,
                hitWeights[i],
                spanI,
                bwRate,
                chnPenGap,
                chnPenSkip,
                readLenA,
                readLenB);
        if(sc == INT32_MIN) {
            break;
        }
        sc += dp[i - 1];

        if(sc < int32_t(spanI)) {
            break;
        }

        dp[i] = sc;
        parent[i] = int32_t(i - 1);
        chainOccurrences[i] = chainOccurrences[i - 1] + 1U;
        if(sc > result.maxScore) {
            result.maxScore = sc;
            result.bestEnd = int32_t(i);
        }
        result.solvedPrefix = i + 1;
    }

    result.fullySolved = (result.solvedPrefix == n);
    return result;
}


// ============================================================================
// MULTI-COPY (MCOPY) CHAIN SELECTION
// ============================================================================
// Reference: Hifiasm ecovlp.cpp:3274, anchor.cpp:1797-1835
//
// Purpose: Extract multiple independent high-quality chains from the same read pair
// - Needed for repetitive regions, segmental duplications, and paralogous sequences
// - Ensures we don't miss alternative alignments when a read has multiple valid mappings
//
// Algorithm:
// 1. Select the best chain and mark all its nodes as "used"
// 2. Iterate through remaining candidates in score-descending order
// 3. For each candidate:
//    - Check if score >= bestScore * mcopyRate (default 0.7)
//    - Check if chain has >= mcopyKhitCutoff markers (default 32)
//    - Trace chain backward until hitting a "used" node
//    - If the independent portion has high enough score, accept it
//    - Mark all nodes in accepted chain as "used"
// 4. Stop when we've selected mcopyNum chains (default 3) or exhausted candidates
//
// Key Parameters (from hifiasm):
// - mcopyNum: Maximum chains to extract (hifiasm: 3)
// - mcopyRate: Minimum score ratio relative to best (hifiasm: 0.7)
// - mcopyKhitCutoff: Minimum markers required to enable mcopy (hifiasm: 32)
//
// This allows capturing biologically real multi-mappings while filtering noise.
// ============================================================================
static inline void applyMcopyFastSelection(
    vector<ThreadScratchpad::ChainCandidate>& chainCandidates,
    const vector<int32_t>& parentSame,
    const vector<int32_t>& parentDiff,
    const vector<int32_t>& dpSame,
    const vector<int32_t>& dpDiff,
    const vector<uint32_t>& chainOccurrencesSame,
    const vector<uint32_t>& chainOccurrencesDiff,
    const uint32_t mcopyNum,
    const double mcopyRate,
    const uint32_t mcopyKhitCutoff,
    vector<uint8_t>& nodeUsed,
    vector<int32_t>& pathNodes,
    vector<ThreadScratchpad::ChainCandidate>& selectedCandidates)
{
    if(chainCandidates.empty() || mcopyNum <= 1) {
        return;
    }

    const auto getOcc = [&](const ThreadScratchpad::ChainCandidate& c) -> uint32_t {
        return c.isDiff ? chainOccurrencesDiff[size_t(c.endK)] : chainOccurrencesSame[size_t(c.endK)];
    };
    const auto& bestCand = chainCandidates.front();
    if(getOcc(bestCand) < mcopyKhitCutoff) {
        return;
    }

    const size_t nodeCount = std::max(dpSame.size(), dpDiff.size());
    nodeUsed.assign(nodeCount, uint8_t(0));
    selectedCandidates.clear();
    selectedCandidates.reserve(std::min<size_t>(chainCandidates.size(), size_t(mcopyNum)));

    auto markChain = [&](const ThreadScratchpad::ChainCandidate& c) {
        const auto& parentArr = c.isDiff ? parentDiff : parentSame;
        for(int32_t k = c.endK; k != -1; k = parentArr[size_t(k)]) {
            nodeUsed[size_t(k)] = uint8_t(1);
        }
    };

    markChain(bestCand);
    selectedCandidates.push_back(bestCand);

    int32_t minScore = int32_t(double(bestCand.score) * mcopyRate);
    if(minScore < 1) {
        minScore = 1;
    }

    for(size_t ci = 1; ci < chainCandidates.size(); ++ci) {
        if(selectedCandidates.size() >= size_t(mcopyNum)) {
            break;
        }
        const auto& cand = chainCandidates[ci];
        if(cand.score < minScore) {
            break;
        }
        if(getOcc(cand) < mcopyKhitCutoff) {
            continue;
        }

        const auto& parentArr = cand.isDiff ? parentDiff : parentSame;
        const auto& dpArr = cand.isDiff ? dpDiff : dpSame;
        pathNodes.clear();
        int32_t ancestor = cand.endK;
        while(ancestor >= 0 && !nodeUsed[size_t(ancestor)]) {
            pathNodes.push_back(ancestor);
            ancestor = parentArr[size_t(ancestor)];
        }
        if(pathNodes.empty()) {
            continue;
        }

        int32_t adjustedScore = cand.score;
        if(ancestor >= 0) {
            adjustedScore -= dpArr[size_t(ancestor)];
        }
        if(adjustedScore < minScore) {
            continue;
        }
        // Keep tiny chains only for the best chain (hifiasm mcopy_fast behavior).
        if(pathNodes.size() <= 1) {
            continue;
        }

        for(const int32_t k : pathNodes) {
            nodeUsed[size_t(k)] = uint8_t(1);
        }
        auto selected = cand;
        selected.score = adjustedScore;
        selectedCandidates.push_back(selected);
    }

    if(selectedCandidates.size() < chainCandidates.size()) {
        chainCandidates = selectedCandidates;
    }
}

// Resolve hit target positions to marker ordinals in O(n log n + m) using:
// 1) one sort of hit indices by target position, then
// 2) one forward scan over marker positions.
// This avoids O(n log m) repeated binary searches on the hot path.
template<class MarkerContainer>
static inline bool mapHitPositionsToMarkerOrdinals(
    const vector<uint32_t>& hitPosB,
    const MarkerContainer& markersB,
    vector<uint32_t>& hitOrdinalB,
    vector<uint32_t>& orderByPosB)
{
    const size_t n = hitPosB.size();
    if (n == 0) {
        return true;
    }
    if (markersB.empty()) {
        return false;
    }

    orderByPosB.resize(n);
    std::iota(orderByPosB.begin(), orderByPosB.end(), uint32_t(0));
    std::sort(orderByPosB.begin(), orderByPosB.end(),
        [&](const uint32_t a, const uint32_t b) {
            if (hitPosB[a] != hitPosB[b]) {
                return hitPosB[a] < hitPosB[b];
            }
            return a < b;
        });

    size_t markerIdx = 0;
    for (const uint32_t hitIdx : orderByPosB) {
        const uint32_t pos = hitPosB[hitIdx];
        while (markerIdx < markersB.size() && markersB[markerIdx].position < pos) {
            ++markerIdx;
        }
        if (markerIdx >= markersB.size() || markersB[markerIdx].position != pos) {
            return false;
        }
        hitOrdinalB[hitIdx] = uint32_t(markerIdx);
    }
    return true;
}

// ============================================================================
// MARKER FREQUENCY WEIGHTING (Hifiasm K-mer Frequency Stratification)
// ============================================================================
// Reference: Hifiasm anchor.cpp:11, anchor.cpp:99
//
// Purpose: Weight shared markers (k-mers) based on their genome-wide frequency
// to balance specificity vs coverage in overlap detection.
//
// Three Frequency Tiers:
// 1. **Low-frequency (Rare/Informative)**: count <= lowFreqThreshold
//    - These are highly specific markers (e.g., unique or near-unique in genome)
//    - Weight: rareKmerWeight = 2 (hifiasm default)
//    - Threshold: coveragePeak * 0.333 (HA_KMER_GOOD_RATIO)
//
// 2. **Normal-frequency**: lowFreqThreshold < count < highFreqThreshold
//    - These are standard, reliable markers
//    - Weight: 1 (baseline)
//    - Middle of the coverage distribution
//
// 3. **High-frequency (Repetitive)**: count >= highFreqThreshold
//    - These are repetitive elements (e.g., transposons, tandem repeats)
//    - Weight: pow(1 + (count / highFreqWeightUnit), 1.1)
//    - Threshold: coveragePeak * 1.667 (hom_cov * (2.0 - HA_KMER_GOOD_RATIO))
//    - Normalized by highFreqWeightUnit = highFreqThreshold * 2
//    - Higher exponent (1.1) allows better discrimination in repeat-dense regions
//
// Rationale:
// - Rare markers get bonus weight because they're highly informative
// - Normal markers get standard weight (1) as baseline
// - Repetitive markers get frequency-dependent weight to maintain signal
//   in repeat-rich regions while not overwhelming the alignment
//
// This implements hifiasm's HIFIASM_NORMAL_W macro behavior for robust
// overlap detection across diverse genomic contexts.
// ============================================================================
static inline uint8_t computeInvertedIndexHitWeight(
    const uint32_t count,
    const uint64_t lowFreqThreshold,
    const uint64_t highFreqThreshold,
    const uint64_t highFreqWeightUnit,
    const uint32_t rareKmerWeight,
    const vector<uint8_t>& weightLut,
    const double weightExponent)
{
    if (count <= lowFreqThreshold) {
        return uint8_t(std::min<uint32_t>(255U, rareKmerWeight));
    }
    if (count >= highFreqThreshold) {
        const uint32_t w = 1U + uint32_t((uint64_t(count) + highFreqWeightUnit - 1ULL) / highFreqWeightUnit);
        if (w < 512U) {
            return weightLut[w];
        }
        return uint8_t(std::min<uint32_t>(255U, uint32_t(std::pow(double(w), weightExponent))));
    }
    return uint8_t(1);
}

// ============================================================================
// CHAINING CONFIGURATION INITIALIZATION
// ============================================================================
// Populate chaining parameters that are shared by discovery and PAF paths.
// Keeping this in one place avoids accidental drift between the two entry points.
//
// This function transfers all inverted-index-related options from the global
// OverlapCandidatesOptions into the InvertedIndexData structure that will be
// used by worker threads during parallel chaining.
//
// Parameters configured:
// - Frequency weighting: lowFreqMultiplier, highFreqMultiplier, weightExponent, rareKmerWeight
// - High-frequency downsampling: downsampleHighFrequencyMarkers, highFrequencySampleDistance, maxHighFrequencyPerStreak
// - Chain selection: highFactor, minNChain (for max_n_chain calculation)
// - DP chaining: lchainIsAccurate, useEcScoring
// - Mcopy extraction: enableMcopyFast, mcopyNum, mcopyRate, mcopyKhitCutoff
// - COV_W control: mcopyOcvWindow, mcopyOcvWeakKeepRatio
// - Overlap validation: minOverlapLength, maxEndFuzz, nonRedundantOverlapFraction
// ============================================================================
template<class InvertedIndexData>
static inline void configureInvertedIndexDataForChaining(
    InvertedIndexData& data,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
    const uint64_t coveragePeak,
    const double maxDriftRate)
{
    data.maxDriftRate = maxDriftRate;
    data.coveragePeak = coveragePeak;
    data.weightExponent = overlapCandidatesOptions.invertedIndexWeightExponent;
    data.lowFreqMultiplier = overlapCandidatesOptions.invertedIndexLowFreqMultiplier;
    data.highFreqMultiplier = overlapCandidatesOptions.invertedIndexHighFreqMultiplier;
    data.rareKmerWeight = overlapCandidatesOptions.invertedIndexRareKmerWeight;
    data.downsampleHighFrequencyMarkers = overlapCandidatesOptions.invertedIndexDownsampleHighFrequencyMarkers;
    data.highFrequencySampleDistance = overlapCandidatesOptions.invertedIndexHighFrequencySampleDistance;
    data.maxHighFrequencyPerStreak = overlapCandidatesOptions.invertedIndexMaxHighFrequencyPerStreak;
    data.highFactor = overlapCandidatesOptions.invertedIndexHighFactor;
    data.minNChain = overlapCandidatesOptions.invertedIndexMinNChain;
    data.nonRedundantOverlapFraction = overlapCandidatesOptions.invertedIndexNonRedundantOverlapFraction;
    data.lchainIsAccurate = overlapCandidatesOptions.invertedIndexLchainIsAccurate;
    data.useEcScoring = overlapCandidatesOptions.invertedIndexUseEcScoring;
    data.enableMcopyFast = overlapCandidatesOptions.invertedIndexEnableMcopyFast;
    data.mcopyNum = overlapCandidatesOptions.invertedIndexMcopyNum;
    data.mcopyRate = overlapCandidatesOptions.invertedIndexMcopyRate;
    data.mcopyKhitCutoff = overlapCandidatesOptions.invertedIndexMcopyKhitCutoff;
    data.mcopyOcvWindow = overlapCandidatesOptions.invertedIndexMcopyOcvWindow;
    data.mcopyOcvWeakKeepRatio = overlapCandidatesOptions.invertedIndexMcopyOcvWeakKeepRatio;
    data.minOverlapLength = overlapCandidatesOptions.minOverlapLength;
    data.maxEndFuzz = overlapCandidatesOptions.maxEndFuzz;
}

// ============================================================================
// WEIGHT LOOKUP TABLE (LUT) CONSTRUCTION
// ============================================================================
// Build a small (512-entry) lookup table for fast weight computation.
//
// Purpose: Precompute pow(i, weightExponent) for i ∈ [0, 511]
// - weightExponent is typically 1.1 (hifiasm default)
// - Used for high-frequency marker weighting: weight = pow(normalized_count, 1.1)
// - LUT avoids expensive pow() calls in the hot path during hit collection
//
// Coverage: The LUT covers normalized counts up to 511
// - For higher counts, we fall back to direct pow() computation
// - In practice, normalized counts rarely exceed 512
//
// Memory: 512 bytes (negligible compared to the inverted index itself)
// ============================================================================
template<class InvertedIndexData>
static inline void rebuildWeightLut(InvertedIndexData& data)
{
    data.weightLut.resize(512);
    for(size_t i = 0; i < data.weightLut.size(); i++) {
        data.weightLut[i] = uint8_t(std::min(255.0, std::pow(double(i), data.weightExponent)));
    }
}

// Release high-memory transient index buffers that are not needed after chaining.
template<class InvertedIndexData>
static inline void clearInvertedIndexTransientData(InvertedIndexData& data)
{
    data.compactOccurrences.clear();
    data.compactOccurrences.shrink_to_fit();
    data.hashTable.clear();
    data.hashTable.shrink_to_fit();
    data.strand0CanonicalKmerIds.clear();
    data.strand0CanonicalKmerIds.shrink_to_fit();
    data.strand0CanonicalOffsets.clear();
    data.strand0CanonicalOffsets.shrink_to_fit();
}

// Private class to encapsulate parallel logic (Codebase Pattern).
class InvertedIndexFinder : public MultithreadedObject<InvertedIndexFinder> {
public:
    InvertedIndexFinder(
        const Reads& reads,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MemoryMapped::VectorOfVectors<KmerId, uint64_t>& markerKmerIds,
        const Assembler::AlignmentCandidatesInvertedIndexData& invertedIndexData,
        MemoryMapped::Vector<OrientedReadPair>& candidates,
        MemoryMapped::Vector<Alignment>& precomputedAlignments,
        uint64_t maxChainLimit,
        uint32_t minChainedMarkerCount,
        uint64_t threadCount
    ) : 
        MultithreadedObject(*this),
        reads(reads),
        markers(markers),
        markerKmerIds(markerKmerIds),
        invertedIndexData(invertedIndexData),
        maxChainLimit(maxChainLimit),
        minChainedMarkerCount(minChainedMarkerCount)
    {
        const ReadId readCount = ReadId(markers.size() / 2); // Indexed by strand 0
        const size_t perThreadReserve = std::max<size_t>(
            10000,
            (size_t(readCount) * 8 + size_t(threadCount) - 1) / size_t(threadCount));

        // 1. Setup per-thread accumulation buffers and scratchpads
        threadCandidates.resize(threadCount);
        for(auto& v : threadCandidates) v.reserve(perThreadReserve);
        threadAlignments.resize(threadCount);
        for(auto& v : threadAlignments) v.reserve(perThreadReserve);
        threadScratchpads.resize(threadCount);

        // 2. Setup parallel workload partitioning
        setupLoadBalancing(readCount, 100);

        // 3. Kick off parallel worker threads
        runThreads(&InvertedIndexFinder::threadFunction, threadCount);

        // 4. Consolidate results: Resize global vectors once and copy results back
        size_t totalCandidatesFound = 0;
        for(const auto& v : threadCandidates) totalCandidatesFound += v.size();
        
        cout << "Discovery search complete. Merging " << totalCandidatesFound << " candidates." << endl;

        size_t candidateWritePos = candidates.size();
        candidates.resize(candidateWritePos + totalCandidatesFound);
        size_t alignmentWritePos = precomputedAlignments.size();
        precomputedAlignments.resize(alignmentWritePos + totalCandidatesFound);
        
        for(size_t i = 0; i < threadCount; i++) {
            const auto& v = threadCandidates[i];
            const auto& a = threadAlignments[i];
            if(!v.empty()) {
                std::copy(v.begin(), v.end(), candidates.begin() + candidateWritePos);
                candidateWritePos += v.size();
                std::copy(a.begin(), a.end(), precomputedAlignments.begin() + alignmentWritePos);
                alignmentWritePos += a.size();
            }
        }
    }

private:
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MemoryMapped::VectorOfVectors<KmerId, uint64_t>& markerKmerIds;
    const Assembler::AlignmentCandidatesInvertedIndexData& invertedIndexData;
    uint64_t maxChainLimit;
    uint32_t minChainedMarkerCount;

    vector<vector<OrientedReadPair>> threadCandidates;
    vector<vector<Alignment>> threadAlignments;
    vector<ThreadScratchpad> threadScratchpads;

    void threadFunction(size_t threadId) {
        vector<OrientedReadPair>& localCandidates = threadCandidates[threadId];
        vector<Alignment>& localAlignments = threadAlignments[threadId];
        ThreadScratchpad& scratch = threadScratchpads[threadId];
        struct EmittedChainedCandidate {
            OrientedReadPair candidate;
            Alignment alignment;
            int32_t score = 0;
            uint8_t overlapType = 0;
        };
        vector<EmittedChainedCandidate> emittedForRead;
        vector<EmittedChainedCandidate> filteredForRead;
        struct PendingHighFrequencyMarker {
            uint64_t startIdx = 0;
            uint32_t count = 0;
            uint32_t posA = 0;
            uint32_t ordinalA = 0;
            uint8_t weight = 1;
        };
        vector<PendingHighFrequencyMarker> highFrequencyStreak;
        highFrequencyStreak.reserve(64);
        
        const uint64_t hashMask = invertedIndexData.hashTable.size() - 1;
        const auto* hashTablePtr = invertedIndexData.hashTable.data();
        const uint64_t kmerLen = invertedIndexData.k;
        const double maxDriftRate = invertedIndexData.maxDriftRate;
        const uint64_t coveragePeak = invertedIndexData.coveragePeak;
        const double weightExponent = invertedIndexData.weightExponent;
        const double lowFreqMultiplier = invertedIndexData.lowFreqMultiplier;
        const double highFreqMultiplier = invertedIndexData.highFreqMultiplier;
        const uint32_t rareKmerWeight = invertedIndexData.rareKmerWeight;
        const uint64_t lowFreqThreshold = std::max<uint64_t>(2ULL, uint64_t(double(coveragePeak) * lowFreqMultiplier));
        const uint64_t highFreqThreshold = std::max<uint64_t>(3ULL, uint64_t(double(coveragePeak) * highFreqMultiplier));
        const uint64_t highFreqWeightUnit = std::max<uint64_t>(1ULL, highFreqThreshold * 2ULL);
        const bool downsampleHighFrequencyMarkers =
            invertedIndexData.downsampleHighFrequencyMarkers &&
            invertedIndexData.highFrequencySampleDistance > 0 &&
            invertedIndexData.maxHighFrequencyPerStreak > 0;
        const uint32_t highFrequencySampleDistance = std::max<uint32_t>(1U, invertedIndexData.highFrequencySampleDistance);
        const uint32_t maxHighFrequencyPerStreak = std::max<uint32_t>(1U, invertedIndexData.maxHighFrequencyPerStreak);
        const double nonRedundantOverlapFraction = invertedIndexData.nonRedundantOverlapFraction;
        const bool lchainIsAccurate = invertedIndexData.lchainIsAccurate;
        const bool useEcScoring = invertedIndexData.useEcScoring;
        const bool enableMcopyFast = invertedIndexData.enableMcopyFast;
        const uint32_t mcopyNum = std::max<uint32_t>(1U, invertedIndexData.mcopyNum);
        const double mcopyRate = std::max<double>(0.0, std::min<double>(1.0, invertedIndexData.mcopyRate));
        const uint32_t mcopyKhitCutoff = std::max<uint32_t>(1U, invertedIndexData.mcopyKhitCutoff);
        const uint32_t mcopyOcvWindow = std::max<uint32_t>(1U, invertedIndexData.mcopyOcvWindow);
        const double mcopyOcvWeakKeepRatio = std::max<double>(0.0, std::min<double>(1.0, invertedIndexData.mcopyOcvWeakKeepRatio));


        uint64_t startBatch, endBatch;
        while(getNextBatch(startBatch, endBatch)) {
            for(ReadId readIdA = ReadId(startBatch); readIdA != ReadId(endBatch); ++readIdA) {
                
                const OrientedReadId orientedReadIdA(readIdA, 0);
                const auto& markersA = markers[orientedReadIdA.getValue()];
                const bool haveCanonicalCache =
                    (size_t(readIdA) + 1 < invertedIndexData.strand0CanonicalOffsets.size());
                const KmerId* canonicalIdsA = nullptr;
                size_t canonicalCountA = 0;
                if(haveCanonicalCache) {
                    const uint64_t b = invertedIndexData.strand0CanonicalOffsets[size_t(readIdA)];
                    const uint64_t e = invertedIndexData.strand0CanonicalOffsets[size_t(readIdA) + 1];
                    if(e >= b && e <= invertedIndexData.strand0CanonicalKmerIds.size()) {
                        canonicalIdsA = invertedIndexData.strand0CanonicalKmerIds.data() + b;
                        canonicalCountA = size_t(e - b);
                    }
                }
                const auto& kmerIdsA = markerKmerIds[orientedReadIdA.getValue()];
                const size_t numMarkersA = canonicalIdsA ?
                    std::min(markersA.size(), canonicalCountA) :
                    std::min(markersA.size(), kmerIdsA.size());
                const uint64_t readLenA = reads.getReadRawSequenceLength(readIdA);
                
                scratch.clear();
                scratch.flatHits.reserve(numMarkersA * 2);
                emittedForRead.clear();
                filteredForRead.clear();
                highFrequencyStreak.clear();
                int64_t lastNonHighBoundaryPos = -1;

                auto computeHitWeight = [&](const uint32_t count) -> uint8_t {
                    return computeInvertedIndexHitWeight(
                        count,
                        lowFreqThreshold,
                        highFreqThreshold,
                        highFreqWeightUnit,
                        rareKmerWeight,
                        invertedIndexData.weightLut,
                        weightExponent);
                };

                auto appendMarkerHits = [&](const PendingHighFrequencyMarker& markerInfo) {
                    const auto* compactOccs = &invertedIndexData.compactOccurrences[markerInfo.startIdx];
                    for (uint32_t j = 0; j < markerInfo.count; ++j) {
                        if (compactOccs[j].readId != readIdA) {
                            scratch.flatHits.push_back(
                                {compactOccs[j].readId, markerInfo.posA, compactOccs[j].position, markerInfo.ordinalA, markerInfo.weight});
                        }
                    }
                };

                // ================================================================
                // HIGH-FREQUENCY MARKER DOWNSAMPLING (Hifiasm Repeat Handling)
                // ================================================================
                // Reference: Hifiasm sketch.cpp:12 (MAX_MAX_HIGH_OCC)
                //
                // Purpose: Prevent memory explosion and computational overhead from
                // consecutive high-frequency (repetitive) markers by intelligently
                // downsampling long stretches while preserving alignment signal.
                //
                // Strategy:
                // 1. When we encounter a "streak" of consecutive high-frequency markers,
                //    we buffer them instead of immediately processing all occurrences
                // 2. When the streak ends (encounter non-high-freq or end-of-read),
                //    we downsample the buffered markers:
                //    - Calculate span = rightBoundary - leftBoundary
                //    - keep = span / highFrequencySampleDistance (typically 500bp)
                //    - Cap at maxHighFrequencyPerStreak (hifiasm: MAX_MAX_HIGH_OCC = 16)
                // 3. Select 'keep' markers evenly distributed across the streak
                //    (e.g., for 100 markers keeping 5: indices 0, 25, 50, 75, 100)
                //
                // Parameters (from hifiasm):
                // - highFrequencySampleDistance: 500bp (sample density)
                // - maxHighFrequencyPerStreak: 16 (MAX_MAX_HIGH_OCC, absolute cap)
                //
                // Why this works:
                // - Preserves positional signal (markers are evenly spaced)
                // - Drastically reduces memory for tandem repeats and transposons
                // - Still maintains enough signal for proper alignment in repeat regions
                //
                // Example: A 5kb tandem repeat with 100 high-freq markers
                // → Sample 5000/500 = 10 markers → Keep 10 evenly-spaced markers
                // ================================================================
                auto flushHighFrequencyStreak = [&](const uint32_t rightBoundaryPos) {
                    if (highFrequencyStreak.empty()) {
                        return;
                    }
                    const uint32_t leftBoundaryPos = (lastNonHighBoundaryPos >= 0) ? uint32_t(lastNonHighBoundaryPos) : 0U;
                    const uint32_t span = (rightBoundaryPos > leftBoundaryPos) ? (rightBoundaryPos - leftBoundaryPos) : 0U;
                    uint32_t keep = uint32_t(double(span) / double(highFrequencySampleDistance) + 0.499);
                    if (keep > maxHighFrequencyPerStreak) {
                        keep = maxHighFrequencyPerStreak;
                    }
                    if (keep == 0) {
                        highFrequencyStreak.clear();
                        return;
                    }

                    const size_t selectedCount = std::min<size_t>(highFrequencyStreak.size(), keep);
                    for (size_t s = 0; s < selectedCount; ++s) {
                        const size_t idx = (uint64_t(s) * highFrequencyStreak.size()) / selectedCount;
                        appendMarkerHits(highFrequencyStreak[idx]);
                    }
                    highFrequencyStreak.clear();
                };

                // ================================================================
                // STEP 1: HIT COLLECTION & FREQUENCY-BASED WEIGHTING
                // ================================================================
                // Reference: Hifiasm anchor.cpp:1476-1540 (get_candidates)
                //
                // Purpose: Find all shared markers between read A and the genome-wide index
                // - Query the inverted index with each k-mer from read A
                // - Collect all partner reads that share this k-mer
                // - Apply frequency-based weighting to balance specificity vs coverage
                // - Downsample high-frequency (repetitive) markers to control memory
                //
                // Process:
                // 1. For each marker in read A:
                //    - Canonicalize the k-mer (min of forward and RC)
                //    - Look up in hash table to find all occurrences across reads
                //    - Count = genome-wide frequency of this k-mer
                // 2. Compute weight based on frequency tier:
                //    - Rare (count ≤ lowFreqThreshold): weight = 2 (informative)
                //    - Normal (lowFreq < count < highFreq): weight = 1 (baseline)
                //    - Repetitive (count ≥ highFreqThreshold): weight = pow(normalized_count, 1.1)
                // 3. If high-frequency marker:
                //    - Buffer it for potential downsampling (prevent memory explosion)
                // 4. Otherwise:
                //    - Flush any buffered streak and append this marker's hits
                //
                // Output: flatHits array containing (readIdB, posA, posB, ordinalA, weight)
                // for all shared markers, ready for DP chaining.
                // ================================================================
                for(size_t i = 0; i < numMarkersA; ++i) {
                    KmerId canonicalKId;
                    if(canonicalIdsA) {
                        canonicalKId = canonicalIdsA[i];
                    } else {
                        KmerId currentKId = kmerIdsA[i];
                        KmerId rcKId = getRcKmerId(currentKId, kmerLen);
                        canonicalKId = (currentKId < rcKId) ? currentKId : rcKId;
                    }

                    const uint32_t posA = markersA[i].position;
                    uint64_t slotIdx = hashKmer(canonicalKId) & hashMask;
                    uint64_t startIdx = 0;
                    uint32_t count = 0;
                    bool found = false;

                    // Search for the K-mer in the direct-addressing hash table
                    while(!hashTablePtr[slotIdx].empty) {
                        if(hashTablePtr[slotIdx].key == canonicalKId) {
                            startIdx = hashTablePtr[slotIdx].start;
                            count = hashTablePtr[slotIdx].count;
                            found = true;
                            break;
                        }
                        slotIdx = (slotIdx + 1) & hashMask;
                    }

                    if(!found) {
                        if(downsampleHighFrequencyMarkers) {
                            flushHighFrequencyStreak(posA);
                            lastNonHighBoundaryPos = posA;
                        }
                        continue;
                    }

                    const uint8_t hitWeight = computeHitWeight(count);
                    if(downsampleHighFrequencyMarkers && count >= highFreqThreshold) {
                        highFrequencyStreak.push_back({startIdx, count, posA, uint32_t(i), hitWeight});
                        continue;
                    }

                    if(downsampleHighFrequencyMarkers) {
                        flushHighFrequencyStreak(posA);
                    }
                    appendMarkerHits({startIdx, count, posA, uint32_t(i), hitWeight});
                    lastNonHighBoundaryPos = posA;
                }
                if(downsampleHighFrequencyMarkers && !highFrequencyStreak.empty()) {
                    const uint32_t readLenABoundary = uint32_t(std::min<uint64_t>(
                        readLenA, uint64_t(std::numeric_limits<uint32_t>::max())));
                    flushHighFrequencyStreak(readLenABoundary);
                }
                
                if(scratch.flatHits.empty()) continue;
                std::sort(scratch.flatHits.begin(), scratch.flatHits.end());

                // ================================================================
                // STEP 2: DP CHAINING PER READ PAIR
                // ================================================================
                // Reference: Hifiasm anchor.cpp:1542-1835 (lchain_qdp + chain selection)
                //
                // Purpose: For each partner read B that shares markers with read A:
                // - Run DP chaining to find optimal collinear marker chains
                // - Process both same-strand and diff-strand (RC) orientations
                // - Extract multiple high-quality chains per pair (mcopy)
                // - Store chain candidates for later per-read filtering
                //
                // Process Flow:
                // 1. Group flatHits by partner readId (hits are sorted by readIdB, posA)
                // 2. For each read pair (A, B):
                //    a. Skip if readIdB <= readIdA (avoid duplicates and self-comparisons)
                //    b. Check minimum anchor count threshold (minChainedMarkerCount)
                //    c. Transfer hit data to SoA (Struct-of-Arrays) for cache efficiency
                //    d. Map hit positions to marker ordinals on read B
                //    e. Run quick linear chaining for monotonic prefixes (optimization)
                //    f. Run full quadratic DP chaining for remaining hits
                //    g. Backtrack to extract chain endpoints
                //    h. Apply mcopy selection to extract multiple independent chains
                //    i. Store chain candidates (NOT filtered by score here!)
                // 3. NO per-pair limiting at this stage (collect all positive-score chains)
                //
                // Key Point: Score-based filtering happens later at PER-READ level (Step 4)
                // ================================================================
                size_t hitIter = 0;
                while(hitIter < scratch.flatHits.size()) {
                    const ReadId readIdB = scratch.flatHits[hitIter].partnerReadId;
                    
                    // Skip mirrored pairs and self-comparisons
                    if(readIdB <= readIdA) {
                        while(hitIter < scratch.flatHits.size() && scratch.flatHits[hitIter].partnerReadId == readIdB) hitIter++;
                        continue;
                    }

                    const size_t startInFlat = hitIter;
                    while(hitIter < scratch.flatHits.size() && scratch.flatHits[hitIter].partnerReadId == readIdB) hitIter++;
                    const size_t numHits = hitIter - startInFlat;
                    if(numHits == 0) continue;
                    // Hifiasm parity: chain_cutoff is a >= threshold.
                    // A pair with exactly minChainedMarkerCount anchors is still eligible.
                    if(minChainedMarkerCount > 0 && numHits < size_t(minChainedMarkerCount)) {
                        continue;
                    }

                    const uint64_t readLenB = reads.getReadRawSequenceLength(readIdB);
                    const OrientedReadId orientedReadB(readIdB, 0);
                    const auto& mB = markers[orientedReadB.getValue()];

                    // 2.1 Struct-of-Arrays (SoA) Transfer for cache-efficient DP access
                    scratch.hitPosA.resize(numHits);
                    scratch.hitPosB.resize(numHits);
                    scratch.hitOrdinalA.resize(numHits);
                    scratch.hitOrdinalB.assign(numHits, std::numeric_limits<uint32_t>::max());
                    scratch.hitWeights.resize(numHits);
                    for(size_t k = 0; k < numHits; ++k) {
                        const auto& h = scratch.flatHits[startInFlat + k];
                        scratch.hitPosA[k] = h.posA; scratch.hitPosB[k] = h.posB;
                        scratch.hitOrdinalA[k] = h.ordinalA; scratch.hitWeights[k] = h.weight;
                    }
                    if(!mapHitPositionsToMarkerOrdinals(
                        scratch.hitPosB, mB, scratch.hitOrdinalB, scratch.hitOrderByPosB)) {
                        continue;
                    }

                    // Pre-allocate/Reset DP work arrays
                    scratch.dpSame.assign(numHits, 0); scratch.dpDiff.assign(numHits, 0);
                    scratch.parentSame.assign(numHits, -1); scratch.parentDiff.assign(numHits, -1);
                    scratch.backtrackVisitSame.assign(numHits, -1); scratch.backtrackVisitDiff.assign(numHits, -1);
                    scratch.chainOccurrencesSame.assign(numHits, 1); scratch.chainOccurrencesDiff.assign(numHits, 1);

                    int32_t maxScSame = 0, maxScDiff = 0;
                    int32_t bestEndIdxSame = -1, bestEndIdxDiff = -1;
                    uint64_t maxExtSame = UINT64_MAX, maxExtDiff = UINT64_MAX;

                    const auto dpOptions = getHifiasmLchainDpOptions(
                        lchainIsAccurate,
                        uint32_t(kmerLen));
                    const int32_t MAX_ITER = dpOptions.maxIter;
                    const int32_t MAX_SKIP = dpOptions.maxSkip;
                    const int32_t MAX_DIST_X = dpOptions.maxDist;
                    const int32_t MAX_DIST_Y = dpOptions.maxDist;
                    const double CHN_PEN_GAP = dpOptions.chnPenGap;
                    const double CHN_PEN_SKIP = dpOptions.chnPenSkip;
                    const double BW_RATE = maxDriftRate;

                    // Utility: Calculate expected coordinate extension (length) for a hit
                    auto getAlignmentLength = [&](uint32_t pA, uint32_t pB) -> uint64_t {
                        uint32_t xB = pA, yB = pB; if (xB <= yB) { yB -= xB; xB = 0; } else { xB -= yB; yB = 0; }
                        uint64_t xR = readLenA - pA - 1, yR = readLenB - pB - 1;
                        uint32_t xE = (xR <= yR) ? (uint32_t)(readLenA - 1) : pA + (uint32_t)yR;
                        return (uint64_t)(xE - xB + 1);
                    };

                    // Hifiasm quick_ck_lchain parity:
                    // fast-chain strict monotonic prefixes and run quadratic DP on the rest.
                    QuickLinearChainResult quickSameResult;
                    QuickLinearChainResult quickDiffResult;
                    if(dpOptions.quickCheck) {
                        quickSameResult = runQuickLinearChainPrefix(
                            false,
                            useEcScoring,
                            scratch.hitPosA,
                            scratch.hitPosB,
                            scratch.hitWeights,
                            uint32_t(kmerLen),
                            BW_RATE,
                            CHN_PEN_GAP,
                            CHN_PEN_SKIP,
                            readLenA,
                            readLenB,
                            scratch.dpSame,
                            scratch.parentSame,
                            scratch.chainOccurrencesSame);
                        quickDiffResult = runQuickLinearChainPrefix(
                            true,
                            useEcScoring,
                            scratch.hitPosA,
                            scratch.hitPosB,
                            scratch.hitWeights,
                            uint32_t(kmerLen),
                            BW_RATE,
                            CHN_PEN_GAP,
                            CHN_PEN_SKIP,
                            readLenA,
                            readLenB,
                            scratch.dpDiff,
                            scratch.parentDiff,
                            scratch.chainOccurrencesDiff);
                        if(quickSameResult.bestEnd >= 0) {
                            maxScSame = quickSameResult.maxScore;
                            bestEndIdxSame = quickSameResult.bestEnd;
                            maxExtSame = getAlignmentLength(
                                scratch.hitPosA[size_t(bestEndIdxSame)],
                                scratch.hitPosB[size_t(bestEndIdxSame)]);
                        }
                        if(quickDiffResult.bestEnd >= 0) {
                            maxScDiff = quickDiffResult.maxScore;
                            bestEndIdxDiff = quickDiffResult.bestEnd;
                            maxExtDiff = getAlignmentLength(
                                scratch.hitPosA[size_t(bestEndIdxDiff)],
                                uint32_t(readLenB - 1 - scratch.hitPosB[size_t(bestEndIdxDiff)]));
                        }
                    }
                    const bool quickSame = quickSameResult.fullySolved;
                    const bool quickDiff = quickDiffResult.fullySolved;
                    const size_t dpStartSame = quickSame ? numHits : quickSameResult.solvedPrefix;
                    const size_t dpStartDiff = quickDiff ? numHits : quickDiffResult.solvedPrefix;
                    const size_t dpStart = std::min(dpStartSame, dpStartDiff);

                    // ================================================================
                    // STEP 2.2: DYNAMIC PROGRAMMING CHAINING (Core Algorithm)
                    // ================================================================
                    // Reference: Hifiasm anchor.cpp:1542-1718 (lchain_qdp)
                    //
                    // Purpose: Find optimal chains of collinear markers using DP scoring
                    // - Process both same-strand and diff-strand (RC) orientations
                    // - Build chains that maximize weighted alignment score
                    // - Apply gap penalties and enforce bandwidth constraints
                    //
                    // Algorithm Overview (for each hit i):
                    // 1. Look back at previous hits j < i within MAX_ITER window
                    // 2. For each valid predecessor j:
                    //    - Compute transition score using hifiasm_comput_sc_ch_ec()
                    //    - Check if positions are colinear within BW_RATE bandwidth
                    //    - Apply gap penalties (CHN_PEN_GAP, CHN_PEN_SKIP)
                    //    - Track best predecessor: dp[i] = max(dp[j] + score(j→i))
                    // 3. Maintain MAX_SKIP pruning: stop if no improvement after MAX_SKIP attempts
                    // 4. Apply rescue mechanism: if we've moved too far on target, rescan for best score
                    //
                    // Key Constraints:
                    // - MAX_ITER: Maximum lookback window (number of hits to consider)
                    // - MAX_SKIP: Early termination if no improvement
                    // - MAX_DIST_X/Y: Distance thresholds for rescue mechanism
                    // - BW_RATE: Bandwidth for diagonal drift tolerance
                    //
                    // Output: For each chain endpoint, stores:
                    // - dp[i]: Maximum score ending at hit i
                    // - parent[i]: Best predecessor for backtracking
                    // - chainOccurrences[i]: Number of markers in chain
                    // ================================================================
                    int32_t st_same = 0, st_diff = 0;
                    int32_t max_ii_same = -1, max_ii_diff = -1;
                    // quick_ck_lchain-style short-circuit:
                    // if both orientations are fully solved by the linear pass, skip quadratic DP.
                    const bool skipFullDp = quickSame && quickDiff;
                    for(size_t i = dpStart; i < numHits && !skipFullDp; ++i) {
                        const bool processSame = (!quickSame && i >= dpStartSame);
                        const bool processDiff = (!quickDiff && i >= dpStartDiff);
                        const uint32_t posAi = scratch.hitPosA[i], posBi = scratch.hitPosB[i];
                        const uint8_t spanI = (uint8_t)std::min((uint32_t)kmerLen, (uint32_t)255);

                        // Initialize with self score (raw span, NOT normalized - hifiasm line 1589)
                        // Weight normalization only applies when scoring between pairs
                        int32_t max_f_same = (int32_t)spanI;
                        int32_t max_f_diff = max_f_same;
                        int32_t max_j_same = -1, max_j_diff = -1;
                        int32_t n_skip_same = 0, n_skip_diff = 0;

                        // Update start position for same-strand (pure index-based, hifiasm line 1591)
                        if(processSame && (int32_t)i - st_same > MAX_ITER) st_same = (int32_t)i - MAX_ITER;

                        // Same-strand DP
                        int32_t end_j_same = st_same;
                        for(int32_t j = (int32_t)i - 1; processSame && j >= st_same; --j) {
                            // Use CURRENT hit's weight (scratch.hitWeights[i]), not previous hit's weight!
                            int32_t sc = useEcScoring ?
                                hifiasm_comput_sc_ch_ec(posAi, posBi, scratch.hitPosA[j], scratch.hitPosB[j],
                                    scratch.hitWeights[i], spanI,
                                    BW_RATE, CHN_PEN_GAP, CHN_PEN_SKIP,
                                    readLenA, readLenB) :
                                hifiasm_comput_sc_ch(posAi, posBi, scratch.hitPosA[j], scratch.hitPosB[j],
                                    scratch.hitWeights[i], spanI,
                                    BW_RATE, CHN_PEN_GAP, CHN_PEN_SKIP,
                                    readLenA, readLenB);
                            if(sc == INT32_MIN) continue;
                            sc += scratch.dpSame[j];
                            if(sc > max_f_same) {
                                max_f_same = sc;
                                max_j_same = j;
                                if(n_skip_same > 0) --n_skip_same;
                            } else if(scratch.backtrackVisitSame[j] == (int32_t)i) {
                                if(++n_skip_same > MAX_SKIP) break;
                            }
                            if(scratch.parentSame[j] >= 0) scratch.backtrackVisitSame[scratch.parentSame[j]] = (int32_t)i;
                            end_j_same = j;
                        }

                        // Rescue mechanism for same-strand (hifiasm lines 846-860)
                        // Track by TARGET position (hitPosB) as hifiasm does
                        if(processSame && (max_ii_same < 0 || (scratch.hitPosB[i] - scratch.hitPosB[max_ii_same] > (uint32_t)MAX_DIST_Y))) {
                            int32_t max_val = INT32_MIN;
                            max_ii_same = -1;
                            for(int32_t j = (int32_t)i - 1; j >= st_same; --j) {
                                if(max_val < scratch.dpSame[j]) {
                                    max_val = scratch.dpSame[j];
                                    max_ii_same = j;
                                }
                            }
                        }
                        if(processSame && max_ii_same >= 0 && max_ii_same < end_j_same) {
                            int32_t tmp = useEcScoring ?
                                hifiasm_comput_sc_ch_ec(posAi, posBi, scratch.hitPosA[max_ii_same], scratch.hitPosB[max_ii_same],
                                    scratch.hitWeights[i], spanI,
                                    BW_RATE, CHN_PEN_GAP, CHN_PEN_SKIP,
                                    readLenA, readLenB) :
                                hifiasm_comput_sc_ch(posAi, posBi, scratch.hitPosA[max_ii_same], scratch.hitPosB[max_ii_same],
                                    scratch.hitWeights[i], spanI,
                                    BW_RATE, CHN_PEN_GAP, CHN_PEN_SKIP,
                                    readLenA, readLenB);
                            if(tmp != INT32_MIN && max_f_same < tmp + scratch.dpSame[max_ii_same]) {
                                max_f_same = tmp + scratch.dpSame[max_ii_same];
                                max_j_same = max_ii_same;
                            }
                        }

                        // Store same-strand result
                        if(processSame) {
                            scratch.dpSame[i] = max_f_same;
                            scratch.parentSame[i] = max_j_same;
                            if(max_j_same >= 0) {
                                scratch.chainOccurrencesSame[i] = scratch.chainOccurrencesSame[max_j_same] + 1;
                            }
                        }
                        // Update max_ii tracker based on target position
                        if(processSame && (max_ii_same < 0 || (scratch.hitPosB[i] - scratch.hitPosB[max_ii_same] <= (uint32_t)MAX_DIST_Y &&
                                                 scratch.dpSame[max_ii_same] < scratch.dpSame[i]))) {
                            max_ii_same = (int32_t)i;
                        }

                        // Update start position for diff-strand (pure index-based)
                        if(processDiff && (int32_t)i - st_diff > MAX_ITER) st_diff = (int32_t)i - MAX_ITER;

                        // Diff-strand DP
                        int32_t end_j_diff = st_diff;
                        for(int32_t j = (int32_t)i - 1; processDiff && j >= st_diff; --j) {
                            // For diff-strand, check if B is decreasing
                            if(scratch.hitPosB[j] <= posBi) continue;

                            int32_t dA = (int32_t)posAi - (int32_t)scratch.hitPosA[j];
                            int32_t dB = (int32_t)scratch.hitPosB[j] - (int32_t)posBi; // Reversed for diff-strand

                            if(dA <= 0 || dA > MAX_DIST_X || dB <= 0 || dB > MAX_DIST_Y) continue;

                            // Gap size and minimum distance
                            int32_t dd = std::abs(dA - dB);
                            int32_t dg = std::min(dA, dB);

                            // Dynamic bandwidth check (only if dd > 16)
                            if(dd > 16) {
                                // For diff-strand, use reversed B coordinates for bandwidth calculation
                                uint32_t posBi_fwd = (uint32_t)(readLenB - 1 - posBi);
                                uint32_t posBj_fwd = (uint32_t)(readLenB - 1 - scratch.hitPosB[j]);
                                int32_t bw_diff = hifiasm_cal_bw(posAi, posBi_fwd, scratch.hitPosA[j], posBj_fwd,
                                                                  BW_RATE, readLenA, readLenB);
                                if(dd > bw_diff) continue;
                            }

                            // Base score using CURRENT hit's weight (comput_sc_ch style)
                            int32_t q_span = (int32_t)spanI;
                            int32_t sc = (q_span < dg) ? q_span : dg;
                            sc = HIFIASM_NORMAL_W(sc, (int32_t)scratch.hitWeights[i]);

                            // Apply gap penalty (comput_sc_ch style)
                            if(dd > 0 || (dg > q_span && dg > 0)) {
                                double lin_pen = CHN_PEN_GAP * (double)dd;
                                double a_pen = ((double)sc) * (((double)dd) / ((double)dg)) / BW_RATE;
                                if(useEcScoring && dd >= 4) {
                                    if(lin_pen < a_pen) {
                                        lin_pen = a_pen;
                                    }
                                } else {
                                    if(lin_pen > a_pen) {
                                        lin_pen = a_pen;
                                    }
                                }
                                lin_pen += CHN_PEN_SKIP * (double)dg;
                                sc -= (int32_t)lin_pen;
                            }

                            sc += scratch.dpDiff[j];
                            if(sc > max_f_diff) {
                                max_f_diff = sc;
                                max_j_diff = j;
                                if(n_skip_diff > 0) --n_skip_diff;
                            } else if(scratch.backtrackVisitDiff[j] == (int32_t)i) {
                                if(++n_skip_diff > MAX_SKIP) break;
                            }
                            if(scratch.parentDiff[j] >= 0) scratch.backtrackVisitDiff[scratch.parentDiff[j]] = (int32_t)i;
                            end_j_diff = j;
                        }

                        // Rescue mechanism for diff-strand
                        // For diff-strand, posB decreases, so check: posB[max_ii] - posB[i]
                        if(processDiff && (max_ii_diff < 0 || (scratch.hitPosB[max_ii_diff] - scratch.hitPosB[i] > (uint32_t)MAX_DIST_Y))) {
                            int32_t max_val = INT32_MIN;
                            max_ii_diff = -1;
                            for(int32_t j = (int32_t)i - 1; j >= st_diff; --j) {
                                if(max_val < scratch.dpDiff[j]) {
                                    max_val = scratch.dpDiff[j];
                                    max_ii_diff = j;
                                }
                            }
                        }
                        if(processDiff && max_ii_diff >= 0 && max_ii_diff < end_j_diff && scratch.hitPosB[max_ii_diff] > posBi) {
                            int32_t dA = (int32_t)posAi - (int32_t)scratch.hitPosA[max_ii_diff];
                            int32_t dB = (int32_t)scratch.hitPosB[max_ii_diff] - (int32_t)posBi;
                            if(dA > 0 && dB > 0) {
                                int32_t dd = std::abs(dA - dB);
                                int32_t dg = std::min(dA, dB);

                                // Check bandwidth (comput_sc_ch style)
                                bool bw_ok = (dd <= 16);
                                if(!bw_ok && dd > 16) {
                                    uint32_t posBi_fwd = (uint32_t)(readLenB - 1 - posBi);
                                    uint32_t posBmax_fwd = (uint32_t)(readLenB - 1 - scratch.hitPosB[max_ii_diff]);
                                    int32_t bw_diff = hifiasm_cal_bw(posAi, posBi_fwd, scratch.hitPosA[max_ii_diff], posBmax_fwd,
                                                                      BW_RATE, readLenA, readLenB);
                                    bw_ok = (dd <= bw_diff);
                                }

                                if(bw_ok) {
                                    // Base score using CURRENT hit's weight
                                    int32_t q_span = (int32_t)spanI;
                                    int32_t tmp = (q_span < dg) ? q_span : dg;
                                    tmp = HIFIASM_NORMAL_W(tmp, (int32_t)scratch.hitWeights[i]);

                                    // Apply gap penalty (comput_sc_ch style)
                                    if(dd > 0 || (dg > q_span && dg > 0)) {
                                        double lin_pen = CHN_PEN_GAP * (double)dd;
                                        double a_pen = ((double)tmp) * (((double)dd) / ((double)dg)) / BW_RATE;
                                        if(useEcScoring && dd >= 4) {
                                            if(lin_pen < a_pen) {
                                                lin_pen = a_pen;
                                            }
                                        } else {
                                            if(lin_pen > a_pen) {
                                                lin_pen = a_pen;
                                            }
                                        }
                                        lin_pen += CHN_PEN_SKIP * (double)dg;
                                        tmp -= (int32_t)lin_pen;
                                    }

                                    tmp += scratch.dpDiff[max_ii_diff];
                                    if(max_f_diff < tmp) {
                                        max_f_diff = tmp;
                                        max_j_diff = max_ii_diff;
                                    }
                                }
                            }
                        }

                        // Store diff-strand result
                        if(processDiff) {
                            scratch.dpDiff[i] = max_f_diff;
                            scratch.parentDiff[i] = max_j_diff;
                            if(max_j_diff >= 0) {
                                scratch.chainOccurrencesDiff[i] = scratch.chainOccurrencesDiff[max_j_diff] + 1;
                            }
                        }
                        // Update max_ii tracker for diff-strand (posB decreases)
                        if(processDiff && (max_ii_diff < 0 || (scratch.hitPosB[max_ii_diff] - scratch.hitPosB[i] <= (uint32_t)MAX_DIST_Y &&
                                                 scratch.dpDiff[max_ii_diff] < scratch.dpDiff[i]))) {
                            max_ii_diff = (int32_t)i;
                        }

                        // Update best scores
                        if(processSame && scratch.dpSame[i] > maxScSame) {
                            maxScSame = scratch.dpSame[i];
                            bestEndIdxSame = (int32_t)i;
                            maxExtSame = getAlignmentLength(posAi, posBi);
                        } else if(processSame && scratch.dpSame[i] == maxScSame && bestEndIdxSame >= 0) {
                            uint64_t ext = getAlignmentLength(posAi, posBi);
                            if(ext < maxExtSame) {
                                bestEndIdxSame = (int32_t)i;
                                maxExtSame = ext;
                            }
                        }

                        if(processDiff && scratch.dpDiff[i] > maxScDiff) {
                            maxScDiff = scratch.dpDiff[i];
                            bestEndIdxDiff = (int32_t)i;
                            maxExtDiff = getAlignmentLength(posAi, (uint32_t)(readLenB - 1 - posBi));
                        } else if(processDiff && scratch.dpDiff[i] == maxScDiff && bestEndIdxDiff >= 0) {
                            uint64_t ext = getAlignmentLength(posAi, (uint32_t)(readLenB - 1 - posBi));
                            if(ext < maxExtDiff) {
                                bestEndIdxDiff = (int32_t)i;
                                maxExtDiff = ext;
                            }
                        }
                    }

                    // ================================================================
                    // STEP 3: CHAIN CANDIDATE COLLECTION (Hifiasm Strategy)
                    // ================================================================
                    // Reference: Hifiasm anchor.cpp:1771-1835
                    //
                    // Key Point: NO per-pair filtering during DP chaining!
                    // - Collect ALL chain endpoints with positive score
                    // - Both same-strand and diff-strand chains are kept
                    // - Limiting happens later at per-READ level (line 1550)
                    //
                    // This differs from alternatives that filter by:
                    //   score >= max(minScore, bestScore * ratio)
                    //
                    // Why no per-pair filtering?
                    // - Each pair typically generates 1-2 chains
                    // - Accumulating across all partners creates large pool
                    // - More efficient to filter once at read level by quality
                    //

                    scratch.chainCandidates.clear();
                    for (size_t k = 0; k < numHits; ++k) {
                        // Collect same-strand chain endpoints
                        if (scratch.dpSame[k] > 0) {
                            scratch.chainCandidates.push_back({
                                scratch.dpSame[k],                          // Chain score
                                getAlignmentLength(scratch.hitPosA[k], scratch.hitPosB[k]),  // Alignment length
                                (int32_t)k,                                 // Endpoint index
                                false});                                    // Same-strand flag
                        }
                        // Collect diff-strand (reverse-complement) chain endpoints
                        if (scratch.dpDiff[k] > 0) {
                            scratch.chainCandidates.push_back({
                                scratch.dpDiff[k],
                                getAlignmentLength(scratch.hitPosA[k], (uint32_t)(readLenB - 1 - scratch.hitPosB[k])),
                                (int32_t)k,
                                true});                                     // Diff-strand flag
                        }
                    }

                    // Sort by score descending (ChainCandidate::operator< implements this)
                    // Higher scores = better quality chains placed first
                    std::sort(scratch.chainCandidates.begin(), scratch.chainCandidates.end());

                    if(enableMcopyFast &&
                        mcopyNum > 1 &&
                        numHits >= size_t(mcopyKhitCutoff))
                    {
                        applyMcopyFastSelection(
                            scratch.chainCandidates,
                            scratch.parentSame,
                            scratch.parentDiff,
                            scratch.dpSame,
                            scratch.dpDiff,
                            scratch.chainOccurrencesSame,
                            scratch.chainOccurrencesDiff,
                            mcopyNum,
                            mcopyRate,
                            mcopyKhitCutoff,
                            scratch.mcopyNodeUsed,
                            scratch.mcopyPathNodes,
                            scratch.mcopySelectedCandidates);
                    }

                    // ================================================================
                    // WEAK-CHAIN SUPPRESSION (Hifiasm Chain Redundancy Removal)
                    // ================================================================
                    // Reference: Hifiasm anchor.cpp:1840 (chain_cutoff logic)
                    //
                    // Purpose: Remove low-quality chains that are redundant with higher-quality chains
                    // - Prevents spurious overlaps from fragmentary/noisy alignments
                    // - Keeps strong chains, removes weak chains that overlap with strong ones
                    // - This is hifiasm's chain_cutoff = 2 behavior
                    //
                    // Classification:
                    // - **Strong chain**: Has >= minChainedMarkerCount anchors (e.g., ≥ 2)
                    // - **Weak chain**: Has < minChainedMarkerCount anchors
                    //
                    // Suppression Strategy:
                    // 1. Check if we have both weak AND strong chains (only run if both present)
                    // 2. For each weak chain:
                    //    - Find all strong chains that overlap it on query (read A)
                    //    - Check if the weak chain is DOMINATED by a strong chain:
                    //      a. Strong chain must have high overlap fraction (OFL = 95%)
                    //      b. Strong chain must be substantially better:
                    //         - Either >= CH_OCC (4) times more anchors
                    //         - Or >= CH_SC (16) higher score
                    //    - If dominated, mark weak chain for suppression
                    // 3. Keep all strong chains + non-dominated weak chains
                    //
                    // Constants (from hifiasm anchor.cpp):
                    // - HIFIASM_OFL = 0.95: Minimum overlap fraction to consider dominance
                    // - HIFIASM_CH_OCC = 4: Anchor count ratio for dominance
                    // - HIFIASM_CH_SC = 16: Score difference threshold for dominance
                    //
                    // This ensures we don't lose real overlaps while filtering noise.
                    // ================================================================
                    if(minChainedMarkerCount >= 2 && !scratch.chainCandidates.empty()) {
                        // Hifiasm lch-style gate: run the expensive weak-overlap suppression
                        // only when we have a mix of weak and strong chains for this pair.
                        bool hasWeakChain = false;
                        bool hasStrongChain = false;
                        for(const auto& cand : scratch.chainCandidates) {
                            const uint32_t occ = cand.isDiff ?
                                scratch.chainOccurrencesDiff[cand.endK] :
                                scratch.chainOccurrencesSame[cand.endK];
                            if(occ < minChainedMarkerCount) {
                                hasWeakChain = true;
                            } else {
                                hasStrongChain = true;
                            }
                            if(hasWeakChain && hasStrongChain) {
                                break;
                            }
                        }
                        const bool runWeakSuppression = hasWeakChain && hasStrongChain;

                        if(runWeakSuppression) {
                            const size_t candidateCount = scratch.chainCandidates.size();
                            auto& metas = scratch.weakMetas;
                            metas.resize(candidateCount);
                            auto& weakIdx = scratch.weakIdxWorkspace;
                            auto& strongIdx = scratch.strongIdxWorkspace;
                            weakIdx.clear();
                            strongIdx.clear();
                            weakIdx.reserve(candidateCount);
                            strongIdx.reserve(candidateCount);

                            for(size_t ci = 0; ci < candidateCount; ++ci) {
                                const auto& cand = scratch.chainCandidates[ci];
                                const auto& parentArr = cand.isDiff ? scratch.parentDiff : scratch.parentSame;
                                int32_t rootK = cand.endK;
                                while(parentArr[rootK] != -1) {
                                    rootK = parentArr[rootK];
                                }

                                uint32_t qs = markersA[scratch.hitOrdinalA[rootK]].position;
                                uint32_t qe = markersA[scratch.hitOrdinalA[cand.endK]].position + uint32_t(kmerLen);
                                if(qe < qs) {
                                    std::swap(qs, qe);
                                }

                                const uint32_t occ = cand.isDiff ?
                                    scratch.chainOccurrencesDiff[cand.endK] :
                                    scratch.chainOccurrencesSame[cand.endK];
                                metas[ci] = ThreadScratchpad::WeakFilterMeta{qs, qe, occ, cand.score};

                                if(occ < minChainedMarkerCount) {
                                    weakIdx.push_back(ci);
                                } else {
                                    strongIdx.push_back(ci);
                                }
                            }

                            if(!weakIdx.empty() && !strongIdx.empty()) {
                                std::sort(strongIdx.begin(), strongIdx.end(),
                                    [&](const size_t a, const size_t b) {
                                        if(metas[a].qs != metas[b].qs) {
                                            return metas[a].qs < metas[b].qs;
                                        }
                                        return metas[a].qe < metas[b].qe;
                                    });

                                auto& strongAnchorBegin = scratch.strongAnchorBegin;
                                auto& strongAnchorEnd = scratch.strongAnchorEnd;
                                auto& strongAnchorStartsFlat = scratch.strongAnchorStartsFlat;
                                strongAnchorBegin.assign(candidateCount, 0);
                                strongAnchorEnd.assign(candidateCount, 0);
                                uint64_t totalStrongAnchors = 0;
                                for(const size_t si : strongIdx) {
                                    totalStrongAnchors += metas[si].occ;
                                }
                                strongAnchorStartsFlat.clear();
                                strongAnchorStartsFlat.reserve(totalStrongAnchors);
                                for(const size_t si : strongIdx) {
                                    const auto& cand = scratch.chainCandidates[si];
                                    const auto& parentArr = cand.isDiff ? scratch.parentDiff : scratch.parentSame;
                                    const uint64_t begin = strongAnchorStartsFlat.size();
                                    scratch.currentChainPath.clear();
                                    for(int32_t k = cand.endK; k != -1; k = parentArr[k]) {
                                        scratch.currentChainPath.push_back(scratch.hitPosA[size_t(k)]);
                                    }
                                    std::reverse(scratch.currentChainPath.begin(), scratch.currentChainPath.end());
                                    strongAnchorStartsFlat.insert(
                                        strongAnchorStartsFlat.end(),
                                        scratch.currentChainPath.begin(),
                                        scratch.currentChainPath.end());
                                    scratch.currentChainPath.clear();
                                    strongAnchorBegin[si] = begin;
                                    strongAnchorEnd[si] = strongAnchorStartsFlat.size();
                                }

                                auto& suppress = scratch.suppressWorkspace;
                                suppress.assign(candidateCount, uint8_t(0));
                                for(const size_t wi : weakIdx) {
                                    const auto& w = metas[wi];
                                    const uint32_t weakSpan = (w.qe > w.qs) ? (w.qe - w.qs) : 0;
                                    const uint32_t minOverlap =
                                        std::max<uint32_t>(16U, uint32_t(double(weakSpan) * HIFIASM_OFL));
                                    const uint64_t minStrongOcc = uint64_t(w.occ) << HIFIASM_CH_OCC; // ocn
                                    const int64_t minStrongScore = int64_t(w.score) * int64_t(HIFIASM_CH_SC);

                                    for(const size_t si : strongIdx) {
                                        const auto& s = metas[si];
                                        if(s.qe <= w.qs) {
                                            continue;
                                        }
                                        if(s.qs >= w.qe) {
                                            break;
                                        }
                                        if(uint64_t(s.occ) < minStrongOcc) {
                                            continue;
                                        }
                                        if(int64_t(s.score) < minStrongScore) {
                                            continue;
                                        }

                                        const uint32_t os = std::max(w.qs, s.qs);
                                        const uint32_t oe = std::min(w.qe, s.qe);
                                        if(oe <= os) {
                                            continue;
                                        }
                                        const uint32_t overlap = oe - os;
                                        if(overlap < minOverlap) {
                                            continue;
                                        }

                                        // Hifiasm r485 guard: count strong-chain anchors whose full span
                                        // lies within [os, oe) <=> anchor start in [os, oe-kmerLen].
                                        const uint64_t begin = strongAnchorBegin[si];
                                        const uint64_t end = strongAnchorEnd[si];
                                        if(begin == end || oe < uint32_t(kmerLen)) {
                                            continue;
                                        }
                                        const uint32_t maxAnchorStart = oe - uint32_t(kmerLen);
                                        if(maxAnchorStart < os) {
                                            continue;
                                        }
                                        auto itBegin = std::lower_bound(
                                            strongAnchorStartsFlat.begin() + begin,
                                            strongAnchorStartsFlat.begin() + end,
                                            os);
                                        auto itEnd = std::upper_bound(
                                            itBegin,
                                            strongAnchorStartsFlat.begin() + end,
                                            maxAnchorStart);
                                        const uint64_t overlapAnchorCount =
                                            uint64_t(std::distance(itBegin, itEnd));
                                        if(overlapAnchorCount < minStrongOcc) {
                                            continue;
                                        }

                                        suppress[wi] = uint8_t(1);
                                        break;
                                    }
                                }

                                bool anySuppressed = false;
                                for(const size_t wi : weakIdx) {
                                    if(suppress[wi]) {
                                        anySuppressed = true;
                                        break;
                                    }
                                }
                                if(anySuppressed) {
                                    scratch.filteredCandidates.clear();
                                    scratch.filteredCandidates.reserve(candidateCount);
                                    for(size_t ci = 0; ci < candidateCount; ++ci) {
                                        if(!suppress[ci]) {
                                            scratch.filteredCandidates.push_back(scratch.chainCandidates[ci]);
                                        }
                                    }
                                    scratch.chainCandidates = std::move(scratch.filteredCandidates);
                                }
                            } // weak/strong candidates present
                        }
                    }

                    // ================================================================
                    // STEP 4: INTERVAL NON-REDUNDANCY & FINAL ALIGNMENT EXTRACTION
                    // ================================================================
                    // Reference: Hifiasm anchor.cpp interval redundancy filtering
                    //
                    // Purpose: Prevent emitting multiple chains that cover nearly identical
                    // query intervals (would create redundant/duplicate overlaps)
                    //
                    // Strategy:
                    // - Track accepted query intervals separately for same-strand and diff-strand
                    // - For each chain candidate (in score-descending order):
                    //   1. Backtrack from endpoint to get full chain path
                    //   2. Extract query interval [qOrdS, qOrdE] (in marker ordinals)
                    //   3. Check if this interval overlaps > 50% with any accepted interval
                    //   4. If redundant: SKIP (don't emit)
                    //   5. If novel: ACCEPT and add to accepted intervals
                    //
                    // This ensures we keep diverse overlaps while avoiding near-duplicates.
                    // ================================================================
                    scratch.acceptedIntervalsSame.clear(); scratch.acceptedIntervalsDiff.clear();

                    // Check if a candidate interval overlaps > nonRedundantOverlapFraction (default 0.5)
                    // with any already-accepted interval on the same strand
                    auto isOverlapLarge = [&](uint32_t qs, uint32_t qe, const vector<ThreadScratchpad::ChainInterval>& accepted) {
                         uint32_t len = qe - qs + 1;
                         for (const auto& iv : accepted) {
                             uint32_t oS = std::max(qs, iv.qs), oE = std::min(qe, iv.qe);
                             if (oE >= oS && (oE - oS + 1) > (uint32_t)(nonRedundantOverlapFraction * double(len))) return true;
                         }
                         return false;
                    };

                    for (const auto& cand : scratch.chainCandidates) {
                        int32_t currK = cand.endK;
                        const auto& parentArr = cand.isDiff ? scratch.parentDiff : scratch.parentSame;
                        scratch.currentChainPath.clear();
                        while(currK != -1) { scratch.currentChainPath.push_back(currK); currK = parentArr[currK]; }
                        std::reverse(scratch.currentChainPath.begin(), scratch.currentChainPath.end());
                        if (scratch.currentChainPath.empty()) continue;

                        uint32_t qOrdS = scratch.hitOrdinalA[scratch.currentChainPath.front()];
                        uint32_t qOrdE = scratch.hitOrdinalA[scratch.currentChainPath.back()];
                        if (cand.isDiff) { if (isOverlapLarge(qOrdS, qOrdE, scratch.acceptedIntervalsDiff)) continue; }
                        else { if (isOverlapLarge(qOrdS, qOrdE, scratch.acceptedIntervalsSame)) continue; }

                        // Extract alignment path
                        Alignment al; 
                        al.ordinals.reserve(scratch.currentChainPath.size());
                        for(uint32_t idxK : scratch.currentChainPath) {
                            const uint32_t ordB = scratch.hitOrdinalB[idxK];
                            if(ordB == std::numeric_limits<uint32_t>::max()) {
                                al.ordinals.clear();
                                break;
                            }
                            al.ordinals.push_back({scratch.hitOrdinalA[idxK], ordB});
                        }

                        if(!al.ordinals.empty()) {
                             uint32_t qPstart = markersA[scratch.hitOrdinalA[scratch.currentChainPath.front()]].position;
                             uint32_t qPend = markersA[scratch.hitOrdinalA[scratch.currentChainPath.back()]].position + (uint32_t)kmerLen;
                             uint32_t bPfirst = mB[al.ordinals.front()[1]].position;
                             uint32_t bPlast = mB[al.ordinals.back()[1]].position;
                             
                             uint32_t tS, tE;
                             if (cand.isDiff) {
                                 uint32_t minB = std::min(bPfirst, bPlast), maxB = std::max(bPfirst, bPlast) + (uint32_t)kmerLen;
                                 tS = (uint32_t)readLenB - maxB; tE = (uint32_t)readLenB - minB;
                             } else {
                                 tS = bPfirst; tE = bPlast + (uint32_t)kmerLen;
                             }

                             // Optional: minimum overlap length in bases (pre-extension span).
                             if(invertedIndexData.minOverlapLength > 0) {
                                 const uint64_t qSpan = uint64_t(qPend) - uint64_t(qPstart);
                                 const uint64_t tSpan = uint64_t(tE) - uint64_t(tS);
                                 if(std::min(qSpan, tSpan) < uint64_t(invertedIndexData.minOverlapLength)) {
                                     continue;
                                 }
                             }

                             // Optional: discard internal overlaps that would require large end extension.
                             if(invertedIndexData.maxEndFuzz > 0) {
                                 const uint32_t leftNeed = std::min(qPstart, tS);
                                 const int64_t qRight = int64_t(readLenA) - int64_t(qPend);
                                 const int64_t tRight = int64_t(readLenB) - int64_t(tE);
                                 const uint32_t rightNeed = uint32_t(std::min<int64_t>(std::max<int64_t>(qRight, 0), std::max<int64_t>(tRight, 0)));
                                 if(leftNeed > invertedIndexData.maxEndFuzz || rightNeed > invertedIndexData.maxEndFuzz) {
                                     continue;
                                 }
                             }

                             uint32_t fQs = qPstart, fQe = qPend, fTs = tS, fTe = tE;
                             if (fQs <= fTs) { fTs -= fQs; fQs = 0; } else { fQs -= fTs; fTs = 0; }
                             int64_t remQ = (int64_t)readLenA - fQe, remT = (int64_t)readLenB - fTe;
                             if (remQ <= remT) { fQe = (uint32_t)readLenA; fTe += (uint32_t)remQ; } else { fTe = (uint32_t)readLenB; fQe += (uint32_t)remT; } 

                             al.qs = fQs; al.qe = fQe;
                             const uint32_t markerCountA = uint32_t(markersA.size());
                             const uint32_t markerCountB = uint32_t(mB.size());

                             bool isSameStrand = true;
                             if (cand.isDiff) {
                                 isSameStrand = false;
                                 const auto p = dinara::rcIntervalToForward(uint32_t(readLenB), fTs, fTe);
                                 al.ts = p.first;
                                 al.te = p.second;
                                 const uint32_t numMB = markerCountB;
                                 for(auto& p : al.ordinals) {
                                     p[1] = numMB - 1 - p[1];
                                 }
                                 scratch.acceptedIntervalsDiff.push_back({qOrdS, qOrdE});
                             } else {
                                 isSameStrand = true;
                                 al.ts = fTs;
                                 al.te = fTe;
                                 scratch.acceptedIntervalsSame.push_back({qOrdS, qOrdE});
                             }

                             // Canonicalize candidate so readIds[0] < readIds[1], and keep alignment consistent.
                             ReadId cand0 = readIdA;
                             ReadId cand1 = readIdB;
                             canonicalizeCandidateAndAlignment(cand0, cand1, isSameStrand, al, markerCountA, markerCountB);
                             const uint8_t overlapType = uint8_t(getOverlapType(
                                 qPstart,
                                 qPend,
                                 uint32_t(readLenA)));
                             emittedForRead.push_back(EmittedChainedCandidate{
                                 OrientedReadPair(cand0, cand1, isSameStrand),
                                 std::move(al),
                                 cand.score,
                                 overlapType});
                        }
                    }
                }

                // ============================================================
                // HIFIASM MAX_N_CHAIN PER-READ FILTERING
                // ============================================================
                // Reference: Hifiasm anchor.cpp:1804-1833, 191-220
                //
                // Apply quality-based filtering at the PER-READ level:
                // - This happens AFTER chaining ALL read pairs
                // - NOT during per-pair DP chaining
                //
                // Algorithm:
                // 1. Sort all overlaps for this read by score (descending)
                // 2. Group by overlap type (0-3): Internal, Left, Right, Contained
                // 3. For each type, if count exceeds maxChainLimit:
                //    - Record threshold score at position maxChainLimit
                //    - Keep only overlaps with score >= threshold
                //
                // Example: maxChainLimit=150, coverage=30x
                // - Type 0 (internal): 200 overlaps → keep top 150 by score
                // - Type 1 (left): 50 overlaps → keep all (< 150)
                // - Type 2 (right): 180 overlaps → keep top 150 by score
                // - Type 3 (contained): 300 overlaps → keep top 150 + weak rescue
                //
                if(maxChainLimit > 0 && emittedForRead.size() > maxChainLimit) {
                    // Sort by: score (desc), overlap type, read IDs
                    std::sort(emittedForRead.begin(), emittedForRead.end(),
                        [](const EmittedChainedCandidate& a, const EmittedChainedCandidate& b) {
                            if(a.score != b.score) {
                                return a.score > b.score;  // Higher score first
                            }
                            if(a.overlapType != b.overlapType) {
                                return a.overlapType < b.overlapType;  // Group by type
                            }
                            if(a.candidate.readIds[0] != b.candidate.readIds[0]) {
                                return a.candidate.readIds[0] < b.candidate.readIds[0];
                            }
                            return a.candidate.readIds[1] < b.candidate.readIds[1];
                        });

                    // Count overlaps per type and record threshold score
                    uint64_t countByType[4] = {0, 0, 0, 0};      // Counts per overlap type
                    int32_t thresholdByType[4] = {0, 0, 0, 0};   // Score threshold per type
                    for(const auto& e : emittedForRead) {
                        const uint32_t type = std::min<uint32_t>(3, e.overlapType);
                        countByType[type]++;
                        // When we hit maxChainLimit for this type, record the score threshold
                        if(countByType[type] == maxChainLimit) {
                            thresholdByType[type] = e.score;
                        }
                    }

                    // ========================================================
                    // HIFIASM COV_W OVERLOAD CONTROL (For Contained Overlaps)
                    // ========================================================
                    // Reference: Hifiasm anchor.cpp:1966-2057
                    //
                    // Special handling for TYPE 3 (containing/contained) overlaps:
                    //
                    // Problem: Type 3 overlaps can saturate the maxChainLimit,
                    // potentially excluding important coverage in under-represented regions
                    //
                    // Solution: Window-based weak overlap rescue
                    // - Divide read into windows of size COV_W (default 3072bp)
                    // - Track coverage usage per window
                    // - Allow weak overlaps (score < threshold) if they add coverage
                    //   to under-saturated windows
                    //
                    // Logic:
                    // - Each window has capacity = window_size * (maxChainLimit / 2)
                    // - Weak overlap kept if: good_coverage / total_coverage >= 0.7
                    // - "good_coverage" = coverage in windows with usage < capacity
                    //
                    // Example: Read length 10kb, COV_W=3072, maxChainLimit=150
                    // - Windows: [0-3072), [3072-6144), [6144-9216), [9216-10000)
                    // - Each window capacity = 3072 * (150/2) = 230,400
                    // - Weak overlap accepted if >= 70% falls in non-saturated windows
                    //
                    const uint32_t readLenA32 = uint32_t(std::min<uint64_t>(
                        readLenA,
                        uint64_t(std::numeric_limits<uint32_t>::max())));
                    const bool useOcvWeakKeep =
                        (mcopyOcvWindow > 0) &&                  // COV_W enabled
                        (readLenA32 >= mcopyOcvWindow) &&        // Read long enough for windowing
                        (countByType[3] >= maxChainLimit);       // Type 3 exceeds limit
                    vector<uint64_t> ocvWindowUsage;
                    if(useOcvWeakKeep) {
                        const uint32_t windowCount = (readLenA32 + mcopyOcvWindow - 1U) / mcopyOcvWindow;
                        ocvWindowUsage.assign(windowCount, uint64_t(0));
                        for(uint32_t w = 0; w < windowCount; ++w) {
                            const uint64_t ws = uint64_t(w) * uint64_t(mcopyOcvWindow);
                            const uint64_t we = std::min<uint64_t>(ws + uint64_t(mcopyOcvWindow), uint64_t(readLenA32));
                            uint64_t cap = (we - ws) * uint64_t(maxChainLimit >> 1);
                            if(cap > uint64_t(std::numeric_limits<uint32_t>::max())) {
                                cap = uint64_t(std::numeric_limits<uint32_t>::max());
                            }
                            ocvWindowUsage[w] = (cap << 32); // high32=cap, low32=used
                        }
                    }

                    // Update window usage counters after accepting an overlap
                    // For each window spanned by the overlap, add the overlap length to used coverage
                    auto updateOcvWindows = [&](const EmittedChainedCandidate& e) {
                        if(!useOcvWeakKeep) {
                            return;
                        }
                        const uint32_t qs = std::min<uint32_t>(readLenA32, e.alignment.qs);
                        const uint32_t qe = std::min<uint32_t>(readLenA32, e.alignment.qe);
                        if(qe <= qs) {
                            return;
                        }
                        for(uint32_t w = qs / mcopyOcvWindow; w < ocvWindowUsage.size(); ++w) {
                            const uint32_t ws = w * mcopyOcvWindow;
                            const uint32_t we = std::min<uint32_t>(ws + mcopyOcvWindow, readLenA32);
                            const uint32_t os = std::max<uint32_t>(qs, ws);
                            const uint32_t oe = std::min<uint32_t>(qe, we);
                            if(oe <= os) {
                                if(ws >= qe) {
                                    break;
                                }
                                continue;
                            }
                            const uint32_t overlap = oe - os;
                            uint64_t v = ocvWindowUsage[w];
                            const uint32_t cap = uint32_t(v >> 32);
                            uint32_t used = uint32_t(v & uint64_t(std::numeric_limits<uint32_t>::max()));
                            const uint64_t sum = uint64_t(used) + uint64_t(overlap);
                            used = (sum > uint64_t(std::numeric_limits<uint32_t>::max())) ?
                                std::numeric_limits<uint32_t>::max() : uint32_t(sum);
                            ocvWindowUsage[w] = (uint64_t(cap) << 32) | uint64_t(used);
                            if(we >= qe) {
                                break;
                            }
                        }
                    };

                    // Decide if a weak (below-threshold) overlap should be rescued
                    // Keep if >= 70% of its coverage falls in non-saturated windows
                    // "good" = coverage in windows with room, "bad" = coverage in saturated windows
                    auto evaluateOcvWeakKeep = [&](const EmittedChainedCandidate& e) -> bool {
                        if(!useOcvWeakKeep) {
                            return false;
                        }
                        const uint32_t qs = std::min<uint32_t>(readLenA32, e.alignment.qs);
                        const uint32_t qe = std::min<uint32_t>(readLenA32, e.alignment.qe);
                        if(qe <= qs) {
                            return false;
                        }
                        uint64_t good = 0;
                        uint64_t bad = 0;
                        for(uint32_t w = qs / mcopyOcvWindow; w < ocvWindowUsage.size(); ++w) {
                            const uint32_t ws = w * mcopyOcvWindow;
                            const uint32_t we = std::min<uint32_t>(ws + mcopyOcvWindow, readLenA32);
                            const uint32_t os = std::max<uint32_t>(qs, ws);
                            const uint32_t oe = std::min<uint32_t>(qe, we);
                            if(oe <= os) {
                                if(ws >= qe) {
                                    break;
                                }
                                continue;
                            }
                            const uint32_t overlap = oe - os;
                            const uint64_t v = ocvWindowUsage[w];
                            const uint32_t cap = uint32_t(v >> 32);
                            const uint32_t used = uint32_t(v & uint64_t(std::numeric_limits<uint32_t>::max()));
                            if(uint64_t(used) + uint64_t(overlap) >= uint64_t(cap)) {
                                bad += uint64_t(overlap);
                            } else {
                                good += uint64_t(overlap);
                            }
                            if(we >= qe) {
                                break;
                            }
                        }
                        const uint64_t total = good + bad;
                        if(total == 0) {
                            return false;
                        }
                        return double(good) >= (double(total) * mcopyOcvWeakKeepRatio);
                    };

                    filteredForRead.clear();
                    filteredForRead.reserve(emittedForRead.size());
                    for(auto& e : emittedForRead) {
                        const uint32_t type = std::min<uint32_t>(3, e.overlapType);
                        bool keep = (e.score >= thresholdByType[type]);
                        if(!keep && type == 3U) {
                            keep = evaluateOcvWeakKeep(e);
                        }
                        if(keep) {
                            updateOcvWindows(e);
                            filteredForRead.push_back(std::move(e));
                        }
                    }
                    emittedForRead.swap(filteredForRead);
                }

                for(auto& e : emittedForRead) {
                    localCandidates.push_back(e.candidate);
                    localAlignments.push_back(std::move(e.alignment));
                }
            }
        }
    }
};
// Explicit template instantiation for the MultithreadedObject base class.
template class MultithreadedObject<InvertedIndexFinder>;

} // namespace dinara
using namespace dinara;

// =============================================================================
// Phase 1-4: Build the inverted index for overlap candidate discovery.
// =============================================================================
void Assembler::buildInvertedIndex(uint64_t threadCount) {
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    performanceLog << timestamp << "Building Inverted Index..." << endl;

    // =========================================================================
    // Phase 1: Inverted Index Construction
    // =========================================================================
    // We build an index mapping each canonical K-mer to all of its occurrences
    // (ReadId, Position). We only index Strand 0 because canonical K-mers are
    // strand-symmetric, halving memory usage. RC matches are detected during
    // the DP chaining phase by observing decreasing target positions.
    // =========================================================================
    checkMarkersAreOpen();
    if(!markerKmerIds->isOpen()) {
        throw runtime_error("Marker KmerIds not available for Inverted Index.");
    }
    invertedIndexData.k = assemblerInfo->k;
    const uint64_t totalMarkers = markers->totalSize();
    cout << "Building Inverted Index for " << totalMarkers << " markers." << endl;

    // We use static load balancing (not dynamic getNextBatch) to ensure that
    // each thread processes the exact same reads in both Pass 1 (Count) and
    // Pass 2 (Fill). This prevents buffer overflows from mismatched counts.
    const ReadId readCount = ReadId(markers->size() / 2);
    vector<uint64_t> threadMarkerCounts(threadCount, 0);
    vector<size_t> threadOffsets(threadCount, 0);
    vector<uint64_t> readMarkerCounts(size_t(readCount), 0);

    // Pass 1: Count markers per thread for memory allocation
    auto countFunction = [&](size_t threadId) {
        ReadId startRead = (ReadId)((uint64_t)readCount * threadId / threadCount);
        ReadId endRead   = (ReadId)((uint64_t)readCount * (threadId + 1) / threadCount);

        uint64_t count = 0;
        for(ReadId rId = startRead; rId != endRead; ++rId) {
            const auto& rMarkers = (*markers)[size_t(rId) << 1];
            const auto& rKmerIds = (*markerKmerIds)[size_t(rId) << 1];
            const size_t n = std::min<size_t>(rMarkers.size(), rKmerIds.size());
            readMarkerCounts[size_t(rId)] = uint64_t(n);
            count += n;
        }
        threadMarkerCounts[threadId] = count;
    };

    vector<std::thread> threads;
    for(size_t i = 0; i < threadCount; i++) threads.emplace_back(countFunction, i);
    for(auto& t : threads) t.join();
    threads.clear();

    // Compute global offsets
    size_t totalMarkersFound = 0;
    for(size_t i = 0; i < threadCount; i++) {
        threadOffsets[i] = totalMarkersFound;
        totalMarkersFound += threadMarkerCounts[i];
    }

    invertedIndexData.strand0CanonicalOffsets.resize(size_t(readCount) + 1, 0);
    for(size_t r = 0; r < size_t(readCount); ++r) {
        invertedIndexData.strand0CanonicalOffsets[r + 1] =
            invertedIndexData.strand0CanonicalOffsets[r] + readMarkerCounts[r];
    }
    if(invertedIndexData.strand0CanonicalOffsets.back() != totalMarkersFound) {
        throw runtime_error("buildInvertedIndex: inconsistent marker counts for canonical cache.");
    }
    
    invertedIndexData.occurrences.resize(totalMarkersFound);
    invertedIndexData.strand0CanonicalKmerIds.resize(totalMarkersFound);
    cout << "Index allocated for " << totalMarkersFound << " occurrences (Strand 0 only)." << endl;

    // Pass 2: Fill the occurrence array
    auto fillFunction = [&](size_t threadId) {
        ReadId startRead = (ReadId)((uint64_t)readCount * threadId / threadCount);
        ReadId endRead   = (ReadId)((uint64_t)readCount * (threadId + 1) / threadCount);
        size_t writeOffset = threadOffsets[threadId];

        for(ReadId rId = startRead; rId != endRead; ++rId) {
            const auto& rMarkers = (*markers)[size_t(rId) << 1];
            const auto& rKmerIds = (*markerKmerIds)[size_t(rId) << 1];
            const size_t n = std::min<size_t>(rMarkers.size(), rKmerIds.size());
            size_t canonicalWriteOffset = size_t(invertedIndexData.strand0CanonicalOffsets[size_t(rId)]);
            for(size_t i = 0; i < n; ++i) {
                KmerId kId = rKmerIds[i];
                KmerId rcKId = getRcKmerId(kId, invertedIndexData.k);
                KmerId canonicalKId = (kId < rcKId) ? kId : rcKId;

                invertedIndexData.occurrences[writeOffset++] = {
                    canonicalKId,
                    rId,
                    rMarkers[i].position
                };
                invertedIndexData.strand0CanonicalKmerIds[canonicalWriteOffset++] = canonicalKId;
            }
        }
    };

    // Run filling process in parallel
    for(size_t i = 0; i < threadCount; i++) {
        threads.emplace_back(fillFunction, i);
    }
    for(auto& t : threads) t.join();
    threads.clear();

    // =========================================================================
    // Phase 2: Parallel LSD Radix Sort
    // =========================================================================
    // We sort the occurrences by K-mer ID to group all markers with the same
    // K-mer together. This uses a parallel Least-Significant-Digit (LSD) radix
    // sort, which is O(N * numBytes) and stable. The parallelization strategy:
    //   1. Each thread computes a local histogram for its data partition.
    //   2. Histograms are combined to compute global bucket offsets.
    //   3. Each thread scatters its data to the correct output positions.
    // =========================================================================
    cout << "Sorting " << invertedIndexData.occurrences.size() << " occurrences..." << endl;
    
    if(!invertedIndexData.occurrences.empty()) {
        const size_t n = invertedIndexData.occurrences.size();
        // KmerId stores 2*k bits; higher bytes are always zero. Sorting only
        // active bytes avoids redundant radix passes for wide KmerId types.
        const size_t numBytes = std::max<size_t>(1, (2 * size_t(invertedIndexData.k) + 7) / 8);
        constexpr size_t bucketCount = 256;
        
        vector<InvertedIndexOccurrence> buffer(n);
        vector<InvertedIndexOccurrence>* src = &invertedIndexData.occurrences;
        vector<InvertedIndexOccurrence>* dst = &buffer;
        vector<size_t> threadHistograms(threadCount * bucketCount);
        vector<size_t> globalBucketCounts(bucketCount);
        vector<size_t> globalBucketOffsets(bucketCount);
        vector<size_t> writeOffsets(threadCount * bucketCount);
        uint8_t radixStage = 0; // 0=histogram, 1=scatter
        bool stopWorkers = false;
        size_t currentByteShift = 0;
        std::barrier syncPoint(ptrdiff_t(threadCount + 1));
        vector<thread> sortThreads;
        sortThreads.reserve(threadCount);
        for(size_t tid = 0; tid < threadCount; ++tid) {
            sortThreads.emplace_back([&, tid]() {
                const size_t start = (n * tid) / threadCount;
                const size_t end = (n * (tid + 1)) / threadCount;
                while(true) {
                    syncPoint.arrive_and_wait();
                    if(stopWorkers) {
                        break;
                    }
                    if(radixStage == 0) {
                        size_t* localHistogram = threadHistograms.data() + tid * bucketCount;
                        std::fill_n(localHistogram, bucketCount, size_t(0));
                        for(size_t i = start; i < end; ++i) {
                            const uint8_t byteVal = uint8_t(((*src)[i].kmerId >> currentByteShift) & 0xFF);
                            localHistogram[byteVal]++;
                        }
                    } else {
                        size_t* localWriteOffsets = writeOffsets.data() + tid * bucketCount;
                        for(size_t i = start; i < end; ++i) {
                            const uint8_t byteVal = uint8_t(((*src)[i].kmerId >> currentByteShift) & 0xFF);
                            (*dst)[localWriteOffsets[byteVal]++] = (*src)[i];
                        }
                    }
                    syncPoint.arrive_and_wait();
                }
            });
        }
        
        for (size_t byteIdx = 0; byteIdx < numBytes; ++byteIdx) {
            currentByteShift = byteIdx * 8;
            
            // 2.1 Parallel Histogram Calculation
            radixStage = 0;
            syncPoint.arrive_and_wait();
            syncPoint.arrive_and_wait();
            
            // 2.2 Global Offset Calculation
            std::fill(globalBucketCounts.begin(), globalBucketCounts.end(), 0);
            for(size_t b = 0; b < bucketCount; ++b) {
                for(size_t tid = 0; tid < threadCount; ++tid) {
                    globalBucketCounts[b] += threadHistograms[tid * bucketCount + b];
                }
            }
            globalBucketOffsets[0] = 0;
            for(size_t b = 1; b < bucketCount; ++b) {
                globalBucketOffsets[b] = globalBucketOffsets[b - 1] + globalBucketCounts[b - 1];
            }
            
            // 2.3 Individual Thread Starting Offsets
            for(size_t b = 0; b < bucketCount; ++b) {
                size_t current = globalBucketOffsets[b];
                for(size_t tid = 0; tid < threadCount; tid++) {
                    writeOffsets[tid * bucketCount + b] = current;
                    current += threadHistograms[tid * bucketCount + b];
                }
            }
            
            // 2.4 Parallel Data Scattering
            radixStage = 1;
            syncPoint.arrive_and_wait();
            syncPoint.arrive_and_wait();

            std::swap(src, dst);
        }

        stopWorkers = true;
        syncPoint.arrive_and_wait();
        for(auto& t : sortThreads) {
            t.join();
        }
        
        if (src != &invertedIndexData.occurrences) {
             invertedIndexData.occurrences = std::move(*src);
        }
    }

    // =========================================================================
    // Phase 3: Hash Table Construction (Open Addressing)
    // =========================================================================
    // Build a power-of-2 sized hash table for O(1) K-mer lookups. Each entry
    // stores the starting index and count of occurrences for a single K-mer.
    // We use linear probing for collision resolution. The load factor is ~50%
    // to balance memory usage and probe length.
    // =========================================================================
    cout << "Allocating Direct Addressing Hash Table..." << endl;
    
    uint64_t numDistinctKmers = 0;
    if(!invertedIndexData.occurrences.empty()) {
        numDistinctKmers = 1;
        KmerId lastKId = invertedIndexData.occurrences[0].kmerId;
        for(size_t i = 1; i < invertedIndexData.occurrences.size(); ++i) {
            if(invertedIndexData.occurrences[i].kmerId != lastKId) {
                numDistinctKmers++;
                lastKId = invertedIndexData.occurrences[i].kmerId;
            }
        }
    }
    cout << "Distinct K-mers found: " << numDistinctKmers << endl;

    uint64_t tableSize = 1;
    while(tableSize < numDistinctKmers * 2) tableSize *= 2; 
    invertedIndexData.hashTable.resize(tableSize); 
    
    if(!invertedIndexData.occurrences.empty()) {
        KmerId currentKmer = invertedIndexData.occurrences[0].kmerId;
        uint64_t blockStart = 0;
        uint64_t mask = tableSize - 1;

        auto insertIntoTable = [&](KmerId key, uint64_t startIdx, uint32_t count) {
            uint64_t slotIdx = hashKmer(key) & mask;
            while(!invertedIndexData.hashTable[slotIdx].empty) {
                slotIdx = (slotIdx + 1) & mask;
            }
            invertedIndexData.hashTable[slotIdx] = {key, startIdx, count, false};
        };

        for(uint64_t i = 1; i < invertedIndexData.occurrences.size(); ++i) {
            if(invertedIndexData.occurrences[i].kmerId != currentKmer) {
                insertIntoTable(currentKmer, blockStart, (uint32_t)(i - blockStart));
                currentKmer = invertedIndexData.occurrences[i].kmerId;
                blockStart = i;
            }
        }
        insertIntoTable(currentKmer, blockStart, (uint32_t)(invertedIndexData.occurrences.size() - blockStart));
    }

    // =========================================================================
    // Phase 4: Index Compaction
    // =========================================================================
    // The full occurrence struct contains the K-mer ID, but after building the
    // hash table, we only need (ReadId, Position) for each occurrence. This
    // compaction reduces the memory footprint by ~3x, which is critical for
    // large datasets. We parallelize this trivially across all occurrences.
    // =========================================================================
    const size_t nMarkersTotal = invertedIndexData.occurrences.size();
    invertedIndexData.compactOccurrences.resize(nMarkersTotal);
    
    auto compactFunction = [&](size_t tid) {
        size_t start = (nMarkersTotal * tid) / threadCount;
        size_t end = (nMarkersTotal * (tid + 1)) / threadCount;
        for(size_t i = start; i < end; ++i) {
            invertedIndexData.compactOccurrences[i] = {
                invertedIndexData.occurrences[i].readId,
                invertedIndexData.occurrences[i].position
            };
        }
    };
    vector<thread> compactThreads;
    for(size_t i = 0; i < threadCount; i++) compactThreads.emplace_back(compactFunction, i);
    for(auto& t : compactThreads) t.join();
    
    // Free high-memory intermediate vector
    invertedIndexData.occurrences.clear();
    invertedIndexData.occurrences.shrink_to_fit();
    cout << "Index construction complete." << endl;
}

// =============================================================================
// Phase 5: Run DP chaining on the built index to find alignment candidates.
// =============================================================================
void Assembler::chainAlignmentCandidates(
    double maxDriftRate,
    uint64_t maxChainLimit,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
    uint32_t minChainedMarkerCount,
    uint64_t threadCount
) {
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const auto startTime = std::chrono::steady_clock::now();
    performanceLog << timestamp << "Starting DP chaining for alignment candidates." << endl;

    // Check that buildInvertedIndex has been called
    if(invertedIndexData.compactOccurrences.empty()) {
        throw runtime_error("chainAlignmentCandidates: buildInvertedIndex must be called first.");
    }

    // =========================================================================
    // Parallel Candidate Search
    // =========================================================================
    // The main search loop. For each Read A (strand 0), we:
    //   1. Look up each of its K-mers in the hash table.
    //   2. Collect all (PartnerReadId, PositionA, PositionB) hits.
    //   3. Sort hits by PartnerReadId, then by PositionA.
    //   4. For each (A, B) pair, run DP chaining to find the best chain(s).
    //   5. Extract alignment coordinates and store the candidate pair.
    // The DP chaining logic is Hifiasm-compatible, using cumulative drift
    // constraints and a penalty based on gap rate.
    // =========================================================================
    alignmentCandidates.candidates.createNew(largeDataName("AlignmentCandidates"), largeDataPageSize);
    alignmentCandidatesAlignmentsData.alignments.createNew(largeDataName("AlignmentCandidatesInvertedIndex"), largeDataPageSize); 
    
    // Safety reserve for performance
    const ReadId readCount = ReadId(markers->size() / 2);
    alignmentCandidates.candidates.reserve(size_t(readCount) * 50);
    alignmentCandidatesAlignmentsData.alignments.reserve(size_t(readCount) * 50); 

    // Keep all tuning parameters synchronized with the PAF chaining path.
    configureInvertedIndexDataForChaining(
        invertedIndexData,
        overlapCandidatesOptions,
        assemblerInfo->kmerDistributionInfo.coveragePeak,
        maxDriftRate);
    rebuildWeightLut(invertedIndexData);

    // Launch the finder across all threads
    InvertedIndexFinder finder(
        getReads(),
        *markers,
        *markerKmerIds,
        invertedIndexData,
        alignmentCandidates.candidates,
        alignmentCandidatesAlignmentsData.alignments,
        maxChainLimit,
        minChainedMarkerCount,
        threadCount
    );

    alignmentCandidates.candidates.unreserve();
    alignmentCandidatesAlignmentsData.alignments.unreserve();
    
    // --- Final Cleanup ---
    // The index is a one-shot structure for this pass; release it aggressively.
    clearInvertedIndexTransientData(invertedIndexData);

    const auto endTime = std::chrono::steady_clock::now();
    const double totalSeconds = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime)).count());
    cout << timestamp << "Alignment discovery completed in " << totalSeconds << " s." << endl;
}

// =============================================================================
// Convenience wrapper that calls both buildInvertedIndex and chainAlignmentCandidates.
// =============================================================================
void Assembler::findAlignmentCandidatesInvertedIndex(
    double maxDriftRate,
    uint64_t maxChainLimit,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
    uint32_t minChainedMarkerCount,
    uint64_t threadCount
) {
    buildInvertedIndex(threadCount);
    chainAlignmentCandidates(maxDriftRate, maxChainLimit, overlapCandidatesOptions, minChainedMarkerCount, threadCount);
}

// =============================================================================
// Chain pre-imported PAF candidates using the inverted index.
// This assumes buildInvertedIndex has been called and alignmentCandidates.candidates
// has been populated by importAlignmentCandidatesFromPaf.
// =============================================================================
void Assembler::chainPafCandidates(
    double maxDriftRate,
    uint64_t maxChainLimit,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
    uint32_t minChainedMarkerCount,
    uint64_t threadCount
) {
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const auto startTime = std::chrono::steady_clock::now();
    performanceLog << timestamp << "Starting DP chaining for PAF-imported candidates." << endl;

    // Check that buildInvertedIndex has been called
    if(invertedIndexData.compactOccurrences.empty() || invertedIndexData.hashTable.empty()) {
        throw runtime_error("chainPafCandidates: buildInvertedIndex must be called first.");
    }

    // Check that PAF candidates have been imported
    if(!alignmentCandidates.candidates.isOpen || alignmentCandidates.candidates.size() == 0) {
        throw runtime_error("chainPafCandidates: No candidates imported. Call importAlignmentCandidatesFromPaf first.");
    }

    cout << timestamp << "Chaining " << alignmentCandidates.candidates.size() << " PAF-imported candidates..." << endl;

    // Match discovery-path configuration so PAF and discovery chaining are comparable.
    configureInvertedIndexDataForChaining(
        invertedIndexData,
        overlapCandidatesOptions,
        assemblerInfo->kmerDistributionInfo.coveragePeak,
        maxDriftRate);
    rebuildWeightLut(invertedIndexData);

    // PAF chaining pipeline:
    //   1) Snapshot imported candidate pairs so output vectors can be rebuilt in-place.
    //   2) For each pair, recollect marker hits from the index (same machinery as discovery).
    //   3) Run Hifiasm-parity DP, keeping only the orientation requested by the PAF pair.
    //   4) Reconstruct a single best chain + alignment per pair.
    //   5) Apply read-level max_n_chain style filtering and publish survivors.
    //
    // Even though input pairs come from PAF, the scoring/filtering code path is shared
    // conceptually with discovery to keep behavior as close as possible.
    // Create output storage for chained alignments
    alignmentCandidatesAlignmentsData.alignments.createNew(largeDataName("AlignmentCandidatesInvertedIndex"), largeDataPageSize);
    alignmentCandidatesAlignmentsData.alignments.reserve(alignmentCandidates.candidates.size());

    // Store original candidates before we modify them
    vector<OrientedReadPair> originalCandidates;
    originalCandidates.reserve(alignmentCandidates.candidates.size());
    for(size_t i = 0; i < alignmentCandidates.candidates.size(); i++) {
        originalCandidates.push_back(alignmentCandidates.candidates[i]);
    }

    // Clear and prepare to repopulate with chained results
    alignmentCandidates.candidates.clear();
    alignmentCandidates.candidates.reserve(originalCandidates.size());

    // Hoist all immutable chaining parameters to local constants before launching threads.
    // This avoids repeated pointer chasing into `invertedIndexData` in inner loops and
    // makes per-thread lambdas easier for the compiler to optimize.
    // Setup threading
    const uint64_t hashMask = invertedIndexData.hashTable.size() - 1;
    const auto* hashTablePtr = invertedIndexData.hashTable.data();
    const uint64_t kmerLen = invertedIndexData.k;
    const double maxDriftRateLocal = invertedIndexData.maxDriftRate;
    const uint64_t coveragePeak = invertedIndexData.coveragePeak;
    const uint64_t lowFreqThreshold = std::max<uint64_t>(
        2ULL, uint64_t(double(coveragePeak) * invertedIndexData.lowFreqMultiplier));
    const uint64_t highFreqThreshold = std::max<uint64_t>(
        3ULL, uint64_t(double(coveragePeak) * invertedIndexData.highFreqMultiplier));
    const uint64_t highFreqWeightUnit = std::max<uint64_t>(1ULL, highFreqThreshold * 2ULL);
    const bool downsampleHighFrequencyMarkers =
        invertedIndexData.downsampleHighFrequencyMarkers &&
        invertedIndexData.highFrequencySampleDistance > 0 &&
        invertedIndexData.maxHighFrequencyPerStreak > 0;
    const uint32_t highFrequencySampleDistance = std::max<uint32_t>(
        1U, invertedIndexData.highFrequencySampleDistance);
    const uint32_t maxHighFrequencyPerStreak = std::max<uint32_t>(
        1U, invertedIndexData.maxHighFrequencyPerStreak);
    const bool lchainIsAccurate = invertedIndexData.lchainIsAccurate;
    const bool useEcScoring = invertedIndexData.useEcScoring;
    const bool enableMcopyFast = invertedIndexData.enableMcopyFast;
    const uint32_t mcopyNum = std::max<uint32_t>(1U, invertedIndexData.mcopyNum);
    const double mcopyRate = std::max<double>(0.0, std::min<double>(1.0, invertedIndexData.mcopyRate));
    const uint32_t mcopyKhitCutoff = std::max<uint32_t>(1U, invertedIndexData.mcopyKhitCutoff);
    const uint32_t mcopyOcvWindow = std::max<uint32_t>(1U, invertedIndexData.mcopyOcvWindow);
    const double mcopyOcvWeakKeepRatio = std::max<double>(0.0, std::min<double>(1.0, invertedIndexData.mcopyOcvWeakKeepRatio));

    struct PafChainedCandidate {
        OrientedReadPair candidate;
        Alignment alignment;
        int32_t score = 0;
        uint8_t overlapType = 0;
        ReadId queryReadId = 0;
    };

    // Thread-local results
    vector<vector<PafChainedCandidate>> threadResults(threadCount);
    const size_t perThreadResultReserve =
        (originalCandidates.size() + threadCount - 1) / threadCount;
    for(auto& v : threadResults) {
        v.reserve(perThreadResultReserve);
    }
    
    // Per-thread processing:
    // each worker owns a scratchpad and writes only to its own result vector.
    // This keeps the hot path lock-free and avoids false sharing.
    const size_t batchSize = std::max(size_t(1), originalCandidates.size() / (threadCount * 10));
    setupLoadBalancing(originalCandidates.size(), batchSize);

    vector<std::thread> threads;
    for(size_t tid = 0; tid < threadCount; tid++) {
        threads.emplace_back([&, tid]() {
            ThreadScratchpad scratch;
            vector<PafChainedCandidate>& localResults = threadResults[tid];
            struct PendingHighFrequencyMarker {
                uint64_t startIdx = 0;
                uint32_t count = 0;
                uint32_t posA = 0;
                uint32_t ordinalA = 0;
                uint8_t weight = 1;
            };
            vector<PendingHighFrequencyMarker> highFrequencyStreak;
            highFrequencyStreak.reserve(64);

            uint64_t startBatch, endBatch;
            while(getNextBatch(startBatch, endBatch)) {
                for(size_t idx = startBatch; idx < endBatch; idx++) {
                    // Each iteration handles one imported PAF pair:
                    // hit recollection -> DP chaining -> chain extraction -> alignment emit.
                    const OrientedReadPair& pair = originalCandidates[idx];
                    const ReadId readIdA = pair.readIds[0];
                    const ReadId readIdB = pair.readIds[1];
                    const bool pafSameStrand = pair.isSameStrand;

                    const OrientedReadId orientedReadIdA(readIdA, 0);
                    const OrientedReadId orientedReadIdB(readIdB, 0);
                    const auto& markersA = (*markers)[orientedReadIdA.getValue()];
                    const auto& markersB = (*markers)[orientedReadIdB.getValue()];
                    const bool haveCanonicalCache =
                        (size_t(readIdA) + 1 < invertedIndexData.strand0CanonicalOffsets.size());
                    const KmerId* canonicalIdsA = nullptr;
                    size_t canonicalCountA = 0;
                    if(haveCanonicalCache) {
                        const uint64_t b = invertedIndexData.strand0CanonicalOffsets[size_t(readIdA)];
                        const uint64_t e = invertedIndexData.strand0CanonicalOffsets[size_t(readIdA) + 1];
                        if(e >= b && e <= invertedIndexData.strand0CanonicalKmerIds.size()) {
                            canonicalIdsA = invertedIndexData.strand0CanonicalKmerIds.data() + b;
                            canonicalCountA = size_t(e - b);
                        }
                    }
                    const auto& kmerIdsA = (*markerKmerIds)[orientedReadIdA.getValue()];
                    
                    const size_t numMarkersA = canonicalIdsA ?
                        std::min(markersA.size(), canonicalCountA) :
                        std::min(markersA.size(), kmerIdsA.size());
                    const size_t numMarkersB = markersB.size();
                    
                    // Without markers on either side, no anchor chain is possible.
                    if(numMarkersA == 0 || numMarkersB == 0) continue;

                    const uint64_t readLenA = reads->getReadRawSequenceLength(readIdA);
                    const uint64_t readLenB = reads->getReadRawSequenceLength(readIdB);

                    // Collect k-mer matches between the pair
                    scratch.clear();
                    scratch.flatHits.reserve(numMarkersA);
                    highFrequencyStreak.clear();
                    int64_t lastNonHighBoundaryPos = -1;

                    auto computeHitWeight = [&](const uint32_t count) -> uint8_t {
                        return computeInvertedIndexHitWeight(
                            count,
                            lowFreqThreshold,
                            highFreqThreshold,
                            highFreqWeightUnit,
                            invertedIndexData.rareKmerWeight,
                            invertedIndexData.weightLut,
                            invertedIndexData.weightExponent);
                    };

                    auto appendMarkerHits = [&](const PendingHighFrequencyMarker& markerInfo) {
                        const auto* compactOccs = &invertedIndexData.compactOccurrences[markerInfo.startIdx];
                        for (uint32_t j = 0; j < markerInfo.count; ++j) {
                            if (compactOccs[j].readId == readIdB) {
                                scratch.flatHits.push_back(
                                    {readIdB, markerInfo.posA, compactOccs[j].position, markerInfo.ordinalA, markerInfo.weight});
                            }
                        }
                    };

                    auto flushHighFrequencyStreak = [&](const uint32_t rightBoundaryPos) {
                        if (highFrequencyStreak.empty()) {
                            return;
                        }
                        const uint32_t leftBoundaryPos = (lastNonHighBoundaryPos >= 0) ? uint32_t(lastNonHighBoundaryPos) : 0U;
                        const uint32_t span = (rightBoundaryPos > leftBoundaryPos) ? (rightBoundaryPos - leftBoundaryPos) : 0U;
                        uint32_t keep = uint32_t(double(span) / double(highFrequencySampleDistance) + 0.499);
                        if (keep > maxHighFrequencyPerStreak) {
                            keep = maxHighFrequencyPerStreak;
                        }
                        if (keep == 0) {
                            highFrequencyStreak.clear();
                            return;
                        }

                        const size_t selectedCount = std::min<size_t>(highFrequencyStreak.size(), keep);
                        for (size_t s = 0; s < selectedCount; ++s) {
                            const size_t idx = (uint64_t(s) * highFrequencyStreak.size()) / selectedCount;
                            appendMarkerHits(highFrequencyStreak[idx]);
                        }
                        highFrequencyStreak.clear();
                    };

                    for(size_t i = 0; i < numMarkersA; i++) {
                        KmerId canonicalKId;
                        if(canonicalIdsA) {
                            canonicalKId = canonicalIdsA[i];
                        } else {
                            KmerId currentKId = kmerIdsA[i];
                            KmerId rcKId = getRcKmerId(currentKId, kmerLen);
                            canonicalKId = (currentKId < rcKId) ? currentKId : rcKId;
                        }
                        
                        const uint32_t posA = markersA[i].position;
                        uint64_t slotIdx = hashKmer(canonicalKId) & hashMask;
                        uint64_t startIdx = 0;
                        uint32_t count = 0;
                        bool found = false;

                        // Search hash table for this k-mer
                        while(!hashTablePtr[slotIdx].empty) {
                            if(hashTablePtr[slotIdx].key == canonicalKId) {
                                startIdx = hashTablePtr[slotIdx].start;
                                count = hashTablePtr[slotIdx].count;
                                found = true;
                                break;
                            }
                            slotIdx = (slotIdx + 1) & hashMask;
                        }

                        if(!found) {
                            if(downsampleHighFrequencyMarkers) {
                                flushHighFrequencyStreak(posA);
                                lastNonHighBoundaryPos = posA;
                            }
                            continue;
                        }

                        const uint8_t hitWeight = computeHitWeight(count);
                        if(downsampleHighFrequencyMarkers && count >= highFreqThreshold) {
                            highFrequencyStreak.push_back({startIdx, count, posA, uint32_t(i), hitWeight});
                            continue;
                        }

                        if(downsampleHighFrequencyMarkers) {
                            flushHighFrequencyStreak(posA);
                        }
                        appendMarkerHits({startIdx, count, posA, uint32_t(i), hitWeight});
                        lastNonHighBoundaryPos = posA;
                    }
                    if(downsampleHighFrequencyMarkers && !highFrequencyStreak.empty()) {
                        const uint32_t readLenABoundary = uint32_t(std::min<uint64_t>(
                            readLenA, uint64_t(std::numeric_limits<uint32_t>::max())));
                        flushHighFrequencyStreak(readLenABoundary);
                    }

                    if(scratch.flatHits.empty()) {
                        // No shared marker k-mers found, so there is no chain to score.
                        // Do not emit a candidate with an empty chain.
                        continue;
                    }

                    // Sort hits by position
                    std::sort(scratch.flatHits.begin(), scratch.flatHits.end(), 
                        [](const InvertedIndexTempHit& a, const InvertedIndexTempHit& b) {
                            return a.posA < b.posA;
                        });

                    const size_t numHits = scratch.flatHits.size();
                    // Hifiasm parity: chain_cutoff is a >= threshold.
                    // A pair with exactly minChainedMarkerCount anchors is still eligible.
                    if(minChainedMarkerCount > 0 && numHits < size_t(minChainedMarkerCount)) {
                        continue;
                    }

                    // Transfer hits from temporary AoS layout to SoA vectors.
                    // SoA layout improves cache locality in the O(N^2) DP loops.
                    scratch.hitPosA.resize(numHits);
                    scratch.hitPosB.resize(numHits);
                    scratch.hitOrdinalA.resize(numHits);
                    scratch.hitOrdinalB.assign(numHits, std::numeric_limits<uint32_t>::max());
                    scratch.hitWeights.resize(numHits);
                    for(size_t k = 0; k < numHits; k++) {
                        const auto& h = scratch.flatHits[k];
                        scratch.hitPosA[k] = h.posA;
                        scratch.hitPosB[k] = h.posB;
                        scratch.hitOrdinalA[k] = h.ordinalA;
                        scratch.hitWeights[k] = h.weight;
                    }
                    if(!mapHitPositionsToMarkerOrdinals(
                        scratch.hitPosB, markersB, scratch.hitOrdinalB, scratch.hitOrderByPosB)) {
                        continue;
                    }

                    // Hifiasm-style DP chaining for this fixed pair.
                    // Unlike discovery, PAF path enforces a single expected orientation
                    // (same/diff) based on the imported pair metadata.
                    scratch.dpSame.assign(numHits, 0);
                    scratch.dpDiff.assign(numHits, 0);
                    scratch.parentSame.assign(numHits, -1);
                    scratch.parentDiff.assign(numHits, -1);
                    scratch.backtrackVisitSame.assign(numHits, -1);
                    scratch.backtrackVisitDiff.assign(numHits, -1);
                    scratch.chainOccurrencesSame.assign(numHits, 1);
                    scratch.chainOccurrencesDiff.assign(numHits, 1);

                    const auto dpOptions = getHifiasmLchainDpOptions(
                        lchainIsAccurate,
                        uint32_t(kmerLen));
                    const int32_t MAX_ITER = dpOptions.maxIter;
                    const int32_t MAX_SKIP = dpOptions.maxSkip;
                    const int32_t MAX_DIST_X = dpOptions.maxDist;
                    const int32_t MAX_DIST_Y = dpOptions.maxDist;
                    const double CHN_PEN_GAP = dpOptions.chnPenGap;
                    const double CHN_PEN_SKIP = dpOptions.chnPenSkip;
                    const double BW_RATE = maxDriftRateLocal;

                    int32_t maxScSame = 0, maxScDiff = 0;
                    const bool runSame = pafSameStrand;
                    const bool runDiff = !pafSameStrand;
                    QuickLinearChainResult quickSameResult;
                    QuickLinearChainResult quickDiffResult;
                    if(dpOptions.quickCheck) {
                        if(runSame) {
                            quickSameResult = runQuickLinearChainPrefix(
                                false,
                                useEcScoring,
                                scratch.hitPosA,
                                scratch.hitPosB,
                                scratch.hitWeights,
                                uint32_t(kmerLen),
                                BW_RATE,
                                CHN_PEN_GAP,
                                CHN_PEN_SKIP,
                                readLenA,
                                readLenB,
                                scratch.dpSame,
                                scratch.parentSame,
                                scratch.chainOccurrencesSame);
                            if(quickSameResult.bestEnd >= 0) {
                                maxScSame = quickSameResult.maxScore;
                            }
                        }
                        if(runDiff) {
                            quickDiffResult = runQuickLinearChainPrefix(
                                true,
                                useEcScoring,
                                scratch.hitPosA,
                                scratch.hitPosB,
                                scratch.hitWeights,
                                uint32_t(kmerLen),
                                BW_RATE,
                                CHN_PEN_GAP,
                                CHN_PEN_SKIP,
                                readLenA,
                                readLenB,
                                scratch.dpDiff,
                                scratch.parentDiff,
                                scratch.chainOccurrencesDiff);
                            if(quickDiffResult.bestEnd >= 0) {
                                maxScDiff = quickDiffResult.maxScore;
                            }
                        }
                    }
                    const bool quickSame = runSame && quickSameResult.fullySolved;
                    const bool quickDiff = runDiff && quickDiffResult.fullySolved;
                    const size_t dpStartSame = runSame ? (quickSame ? numHits : quickSameResult.solvedPrefix) : numHits;
                    const size_t dpStartDiff = runDiff ? (quickDiff ? numHits : quickDiffResult.solvedPrefix) : numHits;
                    const size_t dpStart = std::min(dpStartSame, dpStartDiff);
                    int32_t st_same = 0, st_diff = 0;
                    int32_t max_ii_same = -1, max_ii_diff = -1;
                    // quick_ck_lchain-style short-circuit for the required orientation(s).
                    const bool skipFullDp =
                        (!runSame || quickSame) &&
                        (!runDiff || quickDiff);
                    for(size_t i = dpStart; i < numHits && !skipFullDp; ++i) {
                        const bool processSame = runSame && !quickSame && i >= dpStartSame;
                        const bool processDiff = runDiff && !quickDiff && i >= dpStartDiff;
                        const uint32_t posAi = scratch.hitPosA[i], posBi = scratch.hitPosB[i];
                        const uint8_t spanI = (uint8_t)std::min((uint32_t)kmerLen, (uint32_t)255);

                        // Initialize with self score (raw span, NOT normalized - hifiasm line 1589)
                        // Weight normalization only applies when scoring between pairs
                        int32_t max_f_same = (int32_t)spanI;
                        int32_t max_f_diff = max_f_same;
                        int32_t max_j_same = -1, max_j_diff = -1;
                        int32_t n_skip_same = 0, n_skip_diff = 0;

                        // Update start position for same-strand (pure index-based, hifiasm line 1591)
                        if(processSame && (int32_t)i - st_same > MAX_ITER) st_same = (int32_t)i - MAX_ITER;

                        // Same-strand DP
                        int32_t end_j_same = st_same;
                        for(int32_t j = (int32_t)i - 1; processSame && j >= st_same; --j) {
                            // Use CURRENT hit's weight (scratch.hitWeights[i]), not previous hit's weight!
                            int32_t sc = useEcScoring ?
                                hifiasm_comput_sc_ch_ec(posAi, posBi, scratch.hitPosA[j], scratch.hitPosB[j],
                                    scratch.hitWeights[i], spanI,
                                    BW_RATE, CHN_PEN_GAP, CHN_PEN_SKIP,
                                    readLenA, readLenB) :
                                hifiasm_comput_sc_ch(posAi, posBi, scratch.hitPosA[j], scratch.hitPosB[j],
                                    scratch.hitWeights[i], spanI,
                                    BW_RATE, CHN_PEN_GAP, CHN_PEN_SKIP,
                                    readLenA, readLenB);
                            if(sc == INT32_MIN) continue;
                            sc += scratch.dpSame[j];
                            if(sc > max_f_same) {
                                max_f_same = sc;
                                max_j_same = j;
                                if(n_skip_same > 0) --n_skip_same;
                            } else if(scratch.backtrackVisitSame[j] == (int32_t)i) {
                                if(++n_skip_same > MAX_SKIP) break;
                            }
                            if(scratch.parentSame[j] >= 0) scratch.backtrackVisitSame[scratch.parentSame[j]] = (int32_t)i;
                            end_j_same = j;
                        }

                        // Rescue for same-strand
                        // Track by TARGET position (hitPosB) as hifiasm does
                        if(processSame &&
                            (max_ii_same < 0 || (scratch.hitPosB[i] - scratch.hitPosB[max_ii_same] > (uint32_t)MAX_DIST_Y))) {
                            int32_t max_val = INT32_MIN;
                            max_ii_same = -1;
                            for(int32_t j = (int32_t)i - 1; j >= st_same; --j) {
                                if(max_val < scratch.dpSame[j]) {
                                    max_val = scratch.dpSame[j];
                                    max_ii_same = j;
                                }
                            }
                        }
                        if(processSame && max_ii_same >= 0 && max_ii_same < end_j_same) {
                            int32_t tmp = useEcScoring ?
                                hifiasm_comput_sc_ch_ec(posAi, posBi, scratch.hitPosA[max_ii_same], scratch.hitPosB[max_ii_same],
                                    scratch.hitWeights[i], spanI,
                                    BW_RATE, CHN_PEN_GAP, CHN_PEN_SKIP,
                                    readLenA, readLenB) :
                                hifiasm_comput_sc_ch(posAi, posBi, scratch.hitPosA[max_ii_same], scratch.hitPosB[max_ii_same],
                                    scratch.hitWeights[i], spanI,
                                    BW_RATE, CHN_PEN_GAP, CHN_PEN_SKIP,
                                    readLenA, readLenB);
                            if(tmp != INT32_MIN && max_f_same < tmp + scratch.dpSame[max_ii_same]) {
                                max_f_same = tmp + scratch.dpSame[max_ii_same];
                                max_j_same = max_ii_same;
                            }
                        }

                        if(processSame) {
                            scratch.dpSame[i] = max_f_same;
                            scratch.parentSame[i] = max_j_same;
                            if(max_j_same >= 0) {
                                scratch.chainOccurrencesSame[i] = scratch.chainOccurrencesSame[size_t(max_j_same)] + 1U;
                            }
                        }
                        // Update max_ii tracker based on target position
                        if(processSame &&
                            (max_ii_same < 0 || (scratch.hitPosB[i] - scratch.hitPosB[max_ii_same] <= (uint32_t)MAX_DIST_Y &&
                                                 scratch.dpSame[max_ii_same] < scratch.dpSame[i]))) {
                            max_ii_same = (int32_t)i;
                        }

                        // Diff-strand DP (pure index-based, lchain_dp style)
                        if(processDiff && (int32_t)i - st_diff > MAX_ITER) st_diff = (int32_t)i - MAX_ITER;

                        int32_t end_j_diff = st_diff;
                        for(int32_t j = (int32_t)i - 1; processDiff && j >= st_diff; --j) {
                            if(scratch.hitPosB[j] <= posBi) continue;

                            int32_t dA = (int32_t)posAi - (int32_t)scratch.hitPosA[j];
                            int32_t dB = (int32_t)scratch.hitPosB[j] - (int32_t)posBi;

                            if(dA <= 0 || dA > MAX_DIST_X || dB <= 0 || dB > MAX_DIST_Y) continue;

                            // Gap size and minimum distance
                            int32_t dd = std::abs(dA - dB);
                            int32_t dg = std::min(dA, dB);

                            // Dynamic bandwidth check (only if dd > 16) - chainPafCandidates version
                            if(dd > 16) {
                                uint32_t posBi_fwd = (uint32_t)(readLenB - 1 - posBi);
                                uint32_t posBj_fwd = (uint32_t)(readLenB - 1 - scratch.hitPosB[j]);
                                int32_t bw_diff = hifiasm_cal_bw(posAi, posBi_fwd, scratch.hitPosA[j], posBj_fwd,
                                                                  BW_RATE, readLenA, readLenB);
                                if(dd > bw_diff) continue;
                            }

                            // Base score using CURRENT hit's weight (comput_sc_ch style)
                            int32_t q_span = (int32_t)spanI;
                            int32_t sc = (q_span < dg) ? q_span : dg;
                            sc = HIFIASM_NORMAL_W(sc, (int32_t)scratch.hitWeights[i]);

                            // Apply gap penalty (comput_sc_ch style)
                            if(dd > 0 || (dg > q_span && dg > 0)) {
                                double lin_pen = CHN_PEN_GAP * (double)dd;
                                double a_pen = ((double)sc) * (((double)dd) / ((double)dg)) / BW_RATE;
                                if(useEcScoring && dd >= 4) {
                                    if(lin_pen < a_pen) {
                                        lin_pen = a_pen;
                                    }
                                } else {
                                    if(lin_pen > a_pen) {
                                        lin_pen = a_pen;
                                    }
                                }
                                lin_pen += CHN_PEN_SKIP * (double)dg;
                                sc -= (int32_t)lin_pen;
                            }

                            sc += scratch.dpDiff[j];
                            if(sc > max_f_diff) {
                                max_f_diff = sc;
                                max_j_diff = j;
                                if(n_skip_diff > 0) --n_skip_diff;
                            } else if(scratch.backtrackVisitDiff[j] == (int32_t)i) {
                                if(++n_skip_diff > MAX_SKIP) break;
                            }
                            if(scratch.parentDiff[j] >= 0) scratch.backtrackVisitDiff[scratch.parentDiff[j]] = (int32_t)i;
                            end_j_diff = j;
                        }

                        // Rescue for diff-strand
                        // For diff-strand, posB decreases, so check: posB[max_ii] - posB[i]
                        if(processDiff &&
                            (max_ii_diff < 0 || (scratch.hitPosB[max_ii_diff] - scratch.hitPosB[i] > (uint32_t)MAX_DIST_Y))) {
                            int32_t max_val = INT32_MIN;
                            max_ii_diff = -1;
                            for(int32_t j = (int32_t)i - 1; j >= st_diff; --j) {
                                if(max_val < scratch.dpDiff[j]) {
                                    max_val = scratch.dpDiff[j];
                                    max_ii_diff = j;
                                }
                            }
                        }
                        if(processDiff &&
                            max_ii_diff >= 0 && max_ii_diff < end_j_diff && scratch.hitPosB[max_ii_diff] > posBi) {
                            int32_t dA = (int32_t)posAi - (int32_t)scratch.hitPosA[max_ii_diff];
                            int32_t dB = (int32_t)scratch.hitPosB[max_ii_diff] - (int32_t)posBi;
                            if(dA > 0 && dB > 0) {
                                int32_t dd = std::abs(dA - dB);
                                int32_t dg = std::min(dA, dB);

                                // Check bandwidth (comput_sc_ch style)
                                bool bw_ok = (dd <= 16);
                                if(!bw_ok && dd > 16) {
                                    uint32_t posBi_fwd = (uint32_t)(readLenB - 1 - posBi);
                                    uint32_t posBmax_fwd = (uint32_t)(readLenB - 1 - scratch.hitPosB[max_ii_diff]);
                                    int32_t bw_diff = hifiasm_cal_bw(posAi, posBi_fwd, scratch.hitPosA[max_ii_diff], posBmax_fwd,
                                                                      BW_RATE, readLenA, readLenB);
                                    bw_ok = (dd <= bw_diff);
                                }

                                if(bw_ok) {
                                    // Base score using CURRENT hit's weight
                                    int32_t q_span = (int32_t)spanI;
                                    int32_t tmp = (q_span < dg) ? q_span : dg;
                                    tmp = HIFIASM_NORMAL_W(tmp, (int32_t)scratch.hitWeights[i]);

                                    // Apply gap penalty (comput_sc_ch style)
                                    if(dd > 0 || (dg > q_span && dg > 0)) {
                                        double lin_pen = CHN_PEN_GAP * (double)dd;
                                        double a_pen = ((double)tmp) * (((double)dd) / ((double)dg)) / BW_RATE;
                                        if(useEcScoring && dd >= 4) {
                                            if(lin_pen < a_pen) {
                                                lin_pen = a_pen;
                                            }
                                        } else {
                                            if(lin_pen > a_pen) {
                                                lin_pen = a_pen;
                                            }
                                        }
                                        lin_pen += CHN_PEN_SKIP * (double)dg;
                                        tmp -= (int32_t)lin_pen;
                                    }

                                    tmp += scratch.dpDiff[max_ii_diff];
                                    if(max_f_diff < tmp) {
                                        max_f_diff = tmp;
                                        max_j_diff = max_ii_diff;
                                    }
                                }
                            }
                        }

                        if(processDiff) {
                            scratch.dpDiff[i] = max_f_diff;
                            scratch.parentDiff[i] = max_j_diff;
                            if(max_j_diff >= 0) {
                                scratch.chainOccurrencesDiff[i] = scratch.chainOccurrencesDiff[size_t(max_j_diff)] + 1U;
                            }
                        }
                        // Update max_ii tracker for diff-strand (posB decreases)
                        if(processDiff &&
                            (max_ii_diff < 0 || (scratch.hitPosB[max_ii_diff] - scratch.hitPosB[i] <= (uint32_t)MAX_DIST_Y &&
                                                 scratch.dpDiff[max_ii_diff] < scratch.dpDiff[i]))) {
                            max_ii_diff = (int32_t)i;
                        }

                        if(processSame && scratch.dpSame[i] > maxScSame) {
                            maxScSame = scratch.dpSame[i];
                        }
                        if(processDiff && scratch.dpDiff[i] > maxScDiff) {
                            maxScDiff = scratch.dpDiff[i];
                        }
                    }

                    // Extract best chain based on PAF strand info.
                    // Enforce the orientation given by the PAF record and
                    // apply the same weak-chain cleanup used in discovery chaining.
                    // Hifiasm: collect ALL chains (no per-pair limiting)
                    const bool useSameStrand = pafSameStrand;
                    const bool wantDiff = !useSameStrand;

                    auto getAlignmentLength = [&](uint32_t pA, uint32_t pB) -> uint64_t {
                        uint32_t xB = pA, yB = pB;
                        if (xB <= yB) {
                            yB -= xB;
                            xB = 0;
                        } else {
                            xB -= yB;
                            yB = 0;
                        }
                        uint64_t xR = readLenA - pA - 1, yR = readLenB - pB - 1;
                        uint32_t xE = (xR <= yR) ? (uint32_t)(readLenA - 1) : pA + (uint32_t)yR;
                        return (uint64_t)(xE - xB + 1);
                    };

                    scratch.chainCandidates.clear();
                    for(size_t k = 0; k < numHits; ++k) {
                        // Collect all chains with positive score, matching PAF orientation
                        if(useSameStrand && scratch.dpSame[k] > 0) {
                            scratch.chainCandidates.push_back({
                                scratch.dpSame[k],
                                uint64_t(0), // Filled below.
                                int32_t(k),
                                false});
                        } else if(wantDiff && scratch.dpDiff[k] > 0) {
                            scratch.chainCandidates.push_back({
                                scratch.dpDiff[k],
                                uint64_t(0), // Filled below.
                                int32_t(k),
                                true});
                        }
                    }
                    if(scratch.chainCandidates.empty()) continue;

                    // Fill in alignment lengths
                    for(auto& cand : scratch.chainCandidates) {
                        if(cand.isDiff) {
                            cand.chainLen = getAlignmentLength(
                                scratch.hitPosA[size_t(cand.endK)],
                                uint32_t(readLenB - 1 - scratch.hitPosB[size_t(cand.endK)]));
                        } else {
                            cand.chainLen = getAlignmentLength(
                                scratch.hitPosA[size_t(cand.endK)],
                                scratch.hitPosB[size_t(cand.endK)]);
                        }
                    }

                    // Sort by score descending
                    std::sort(scratch.chainCandidates.begin(), scratch.chainCandidates.end());

                    if(enableMcopyFast &&
                        mcopyNum > 1 &&
                        numHits >= size_t(mcopyKhitCutoff))
                    {
                        applyMcopyFastSelection(
                            scratch.chainCandidates,
                            scratch.parentSame,
                            scratch.parentDiff,
                            scratch.dpSame,
                            scratch.dpDiff,
                            scratch.chainOccurrencesSame,
                            scratch.chainOccurrencesDiff,
                            mcopyNum,
                            mcopyRate,
                            mcopyKhitCutoff,
                            scratch.mcopyNodeUsed,
                            scratch.mcopyPathNodes,
                            scratch.mcopySelectedCandidates);
                    }

                    // Hifiasm-style weak-chain suppression (anchor.cpp: OFL/CH_OCC/CH_SC logic).
                    if(minChainedMarkerCount >= 2 && !scratch.chainCandidates.empty()) {
                        // Hifiasm lch-style gate: run weak suppression only when both
                        // weak and strong chains are present for this read pair.
                        bool hasWeakChain = false;
                        bool hasStrongChain = false;
                        for(const auto& cand : scratch.chainCandidates) {
                            const uint32_t occ = cand.isDiff ?
                                scratch.chainOccurrencesDiff[cand.endK] :
                                scratch.chainOccurrencesSame[cand.endK];
                            if(occ < minChainedMarkerCount) {
                                hasWeakChain = true;
                            } else {
                                hasStrongChain = true;
                            }
                            if(hasWeakChain && hasStrongChain) {
                                break;
                            }
                        }
                        const bool runWeakSuppression = hasWeakChain && hasStrongChain;

                        if(runWeakSuppression) {
                            const size_t candidateCount = scratch.chainCandidates.size();
                            auto& metas = scratch.weakMetas;
                            metas.resize(candidateCount);
                            auto& weakIdx = scratch.weakIdxWorkspace;
                            auto& strongIdx = scratch.strongIdxWorkspace;
                            weakIdx.clear();
                            strongIdx.clear();
                            weakIdx.reserve(candidateCount);
                            strongIdx.reserve(candidateCount);

                            for(size_t ci = 0; ci < candidateCount; ++ci) {
                                const auto& cand = scratch.chainCandidates[ci];
                                const auto& parentArrTmp = cand.isDiff ? scratch.parentDiff : scratch.parentSame;
                                int32_t rootK = cand.endK;
                                while(parentArrTmp[rootK] != -1) {
                                    rootK = parentArrTmp[rootK];
                                }

                                uint32_t qs = markersA[scratch.hitOrdinalA[rootK]].position;
                                uint32_t qe = markersA[scratch.hitOrdinalA[cand.endK]].position + uint32_t(kmerLen);
                                if(qe < qs) {
                                    std::swap(qs, qe);
                                }

                                const uint32_t occ = cand.isDiff ?
                                    scratch.chainOccurrencesDiff[cand.endK] :
                                    scratch.chainOccurrencesSame[cand.endK];
                                metas[ci] = ThreadScratchpad::WeakFilterMeta{qs, qe, occ, cand.score};

                                if(occ < minChainedMarkerCount) {
                                    weakIdx.push_back(ci);
                                } else {
                                    strongIdx.push_back(ci);
                                }
                            }

                            if(!weakIdx.empty() && !strongIdx.empty()) {
                                std::sort(strongIdx.begin(), strongIdx.end(),
                                    [&](const size_t a, const size_t b) {
                                        if(metas[a].qs != metas[b].qs) {
                                            return metas[a].qs < metas[b].qs;
                                        }
                                        return metas[a].qe < metas[b].qe;
                                    });

                                auto& strongAnchorBegin = scratch.strongAnchorBegin;
                                auto& strongAnchorEnd = scratch.strongAnchorEnd;
                                auto& strongAnchorStartsFlat = scratch.strongAnchorStartsFlat;
                                strongAnchorBegin.assign(candidateCount, 0);
                                strongAnchorEnd.assign(candidateCount, 0);
                                uint64_t totalStrongAnchors = 0;
                                for(const size_t si : strongIdx) {
                                    totalStrongAnchors += metas[si].occ;
                                }
                                strongAnchorStartsFlat.clear();
                                strongAnchorStartsFlat.reserve(totalStrongAnchors);
                                for(const size_t si : strongIdx) {
                                    const auto& cand = scratch.chainCandidates[si];
                                    const auto& parentArrTmp = cand.isDiff ? scratch.parentDiff : scratch.parentSame;
                                    const uint64_t begin = strongAnchorStartsFlat.size();
                                    scratch.currentChainPath.clear();
                                    for(int32_t k = cand.endK; k != -1; k = parentArrTmp[k]) {
                                        scratch.currentChainPath.push_back(scratch.hitPosA[size_t(k)]);
                                    }
                                    std::reverse(scratch.currentChainPath.begin(), scratch.currentChainPath.end());
                                    strongAnchorStartsFlat.insert(
                                        strongAnchorStartsFlat.end(),
                                        scratch.currentChainPath.begin(),
                                        scratch.currentChainPath.end());
                                    scratch.currentChainPath.clear();
                                    strongAnchorBegin[si] = begin;
                                    strongAnchorEnd[si] = strongAnchorStartsFlat.size();
                                }

                                auto& suppress = scratch.suppressWorkspace;
                                suppress.assign(candidateCount, uint8_t(0));
                                for(const size_t wi : weakIdx) {
                                    const auto& w = metas[wi];
                                    const uint32_t weakSpan = (w.qe > w.qs) ? (w.qe - w.qs) : 0;
                                    const uint32_t minOverlap =
                                        std::max<uint32_t>(16U, uint32_t(double(weakSpan) * HIFIASM_OFL));
                                    const uint64_t minStrongOcc = uint64_t(w.occ) << HIFIASM_CH_OCC; // ocn
                                    const int64_t minStrongScore = int64_t(w.score) * int64_t(HIFIASM_CH_SC);

                                    for(const size_t si : strongIdx) {
                                        const auto& s = metas[si];
                                        if(s.qe <= w.qs) {
                                            continue;
                                        }
                                        if(s.qs >= w.qe) {
                                            break;
                                        }
                                        if(uint64_t(s.occ) < minStrongOcc) {
                                            continue;
                                        }
                                        if(int64_t(s.score) < minStrongScore) {
                                            continue;
                                        }

                                        const uint32_t os = std::max(w.qs, s.qs);
                                        const uint32_t oe = std::min(w.qe, s.qe);
                                        if(oe <= os) {
                                            continue;
                                        }
                                        const uint32_t overlap = oe - os;
                                        if(overlap < minOverlap) {
                                            continue;
                                        }

                                        const uint64_t begin = strongAnchorBegin[si];
                                        const uint64_t end = strongAnchorEnd[si];
                                        if(begin == end || oe < uint32_t(kmerLen)) {
                                            continue;
                                        }
                                        const uint32_t maxAnchorStart = oe - uint32_t(kmerLen);
                                        if(maxAnchorStart < os) {
                                            continue;
                                        }
                                        auto itBegin = std::lower_bound(
                                            strongAnchorStartsFlat.begin() + begin,
                                            strongAnchorStartsFlat.begin() + end,
                                            os);
                                        auto itEnd = std::upper_bound(
                                            itBegin,
                                            strongAnchorStartsFlat.begin() + end,
                                            maxAnchorStart);
                                        const uint64_t overlapAnchorCount =
                                            uint64_t(std::distance(itBegin, itEnd));
                                        if(overlapAnchorCount < minStrongOcc) {
                                            continue;
                                        }

                                        suppress[wi] = uint8_t(1);
                                        break;
                                    }
                                }

                                bool anySuppressed = false;
                                for(const size_t wi : weakIdx) {
                                    if(suppress[wi]) {
                                        anySuppressed = true;
                                        break;
                                    }
                                }
                                if(anySuppressed) {
                                    scratch.filteredCandidates.clear();
                                    scratch.filteredCandidates.reserve(candidateCount);
                                    for(size_t ci = 0; ci < candidateCount; ++ci) {
                                        if(!suppress[ci]) {
                                            scratch.filteredCandidates.push_back(scratch.chainCandidates[ci]);
                                        }
                                    }
                                    scratch.chainCandidates = std::move(scratch.filteredCandidates);
                                }
                            }
                        }
                    }

                    int32_t bestEndIdx = -1;
                    int32_t selectedScore = 0;
                    for(const auto& cand : scratch.chainCandidates) {
                        if(cand.isDiff != wantDiff) {
                            continue;
                        }
                        if(minChainedMarkerCount > 0) {
                            const uint32_t occ = cand.isDiff ?
                                scratch.chainOccurrencesDiff[cand.endK] :
                                scratch.chainOccurrencesSame[cand.endK];
                            if(occ < minChainedMarkerCount) {
                                continue;
                            }
                        }
                        bestEndIdx = cand.endK;
                        selectedScore = cand.score;
                        break;
                    }
                    if(bestEndIdx < 0) {
                        continue;
                    }
                    const auto& parentArr = useSameStrand ? scratch.parentSame : scratch.parentDiff;

                    // Backtrack to get chain
                    scratch.currentChainPath.clear();
                    int32_t currK = bestEndIdx;
                    while(currK != -1) {
                        scratch.currentChainPath.push_back(currK);
                        currK = parentArr[currK];
                    }
                    std::reverse(scratch.currentChainPath.begin(), scratch.currentChainPath.end());

                    // Hifiasm parity: keep chains with length >= chain_cutoff.
                    if(minChainedMarkerCount > 0 &&
                        scratch.currentChainPath.size() < size_t(minChainedMarkerCount)) {
                        continue;
                    }

                    // Build alignment
                    Alignment al;
                    al.ordinals.reserve(scratch.currentChainPath.size());

                    for(uint32_t idxK : scratch.currentChainPath) {
                        const uint32_t ordB = scratch.hitOrdinalB[idxK];
                        if(ordB == std::numeric_limits<uint32_t>::max()) {
                            al.ordinals.clear();
                            break;
                        }
                        al.ordinals.push_back({scratch.hitOrdinalA[idxK], ordB});
                    }

                    if(!al.ordinals.empty()) {
                        // Set alignment coordinates.
                        // Note: Alignment stores ts/te in forward-read coordinates, even for reverse overlaps.
                        uint32_t qPstart = markersA[al.ordinals.front()[0]].position;
                        uint32_t qPend = markersA[al.ordinals.back()[0]].position + (uint32_t)kmerLen;
                        uint32_t tPstart = markersB[al.ordinals.front()[1]].position;
                        uint32_t tPend = markersB[al.ordinals.back()[1]].position + (uint32_t)kmerLen;

                        al.qs = qPstart;
                        al.qe = qPend;
                        if(useSameStrand) {
                            al.ts = tPstart;
                            al.te = tPend;
                        } else {
                            // In the diff-strand chain, the target marker positions are decreasing.
                            // Convert the (strand-0) forward interval to the canonical forward interval.
                            const uint32_t minB = std::min(tPstart, tPend >= uint32_t(kmerLen) ? (tPend - uint32_t(kmerLen)) : tPend);
                            const uint32_t maxB = std::max(tPstart, tPend >= uint32_t(kmerLen) ? (tPend - uint32_t(kmerLen)) : tPend) + uint32_t(kmerLen);
                            const uint32_t tsExt = uint32_t(readLenB) - maxB;
                            const uint32_t teExt = uint32_t(readLenB) - minB;
                            const auto p = dinara::rcIntervalToForward(uint32_t(readLenB), tsExt, teExt);
                            al.ts = p.first;
                            al.te = p.second;
                        }

                        // Optional: minimum overlap length in bases (pre-extension span).
                        if(invertedIndexData.minOverlapLength > 0) {
                            const uint64_t qSpan = uint64_t(al.qe) - uint64_t(al.qs);
                            const uint64_t tSpan = uint64_t(al.te) - uint64_t(al.ts);
                            if(std::min(qSpan, tSpan) < uint64_t(invertedIndexData.minOverlapLength)) {
                                continue;
                            }
                        }

                        // Optional: discard internal overlaps that would require large end extension.
                        if(invertedIndexData.maxEndFuzz > 0) {
                            // For reverse overlaps, Alignment stores ts/te in forward coordinates (ts < te),
                            // but the "end extension" heuristic must be computed in the overlap orientation.
                            // In the reverse-complement coordinate system, the target interval becomes:
                            // [lenB - te, lenB - ts], so the left overhang is (lenB - te) and the right overhang is ts.
                            const int64_t qRight = int64_t(readLenA) - int64_t(al.qe);
                            const uint32_t qRightNeed = uint32_t(std::max<int64_t>(qRight, 0));

                            uint32_t tLeftNeed = 0;
                            uint32_t tRightNeed = 0;
                            if(useSameStrand) {
                                tLeftNeed = al.ts;
                                const int64_t tRight = int64_t(readLenB) - int64_t(al.te);
                                tRightNeed = uint32_t(std::max<int64_t>(tRight, 0));
                            } else {
                                tLeftNeed = uint32_t(readLenB) - al.te;
                                tRightNeed = al.ts;
                            }

                            const uint32_t leftNeed = std::min(al.qs, tLeftNeed);
                            const uint32_t rightNeed = std::min(qRightNeed, tRightNeed);
                            if(leftNeed > invertedIndexData.maxEndFuzz || rightNeed > invertedIndexData.maxEndFuzz) {
                                continue;
                            }
                        }

                        // Flip ordinals for opposite strand
                        if(!useSameStrand) {
                            uint32_t numMB = (uint32_t)markersB.size();
                            for(auto& p : al.ordinals) p[1] = numMB - 1 - p[1];
                        }

                        // Canonicalize candidate so readIds[0] < readIds[1], and keep alignment consistent.
                        const uint32_t markerCountA = uint32_t(markersA.size());
                        const uint32_t markerCountB = uint32_t(markersB.size());
                        ReadId cand0 = readIdA;
                        ReadId cand1 = readIdB;
                        bool isSameStrand = useSameStrand;
                        const uint32_t readLenA32 = uint32_t(std::min<uint64_t>(
                            readLenA,
                            uint64_t(std::numeric_limits<uint32_t>::max())));
                        const uint8_t overlapType = uint8_t(getOverlapType(qPstart, qPend, readLenA32));
                        canonicalizeCandidateAndAlignment(cand0, cand1, isSameStrand, al, markerCountA, markerCountB);

                        localResults.push_back(PafChainedCandidate{
                            OrientedReadPair(cand0, cand1, isSameStrand),
                            std::move(al),
                            selectedScore,
                            overlapType,
                            readIdA
                        });
                    }
                }
            }
        });
    }

    for(auto& t : threads) t.join();

    // Merge thread-local results into a single vector for read-level filtering.
    // This is done once after all workers finish to keep worker code write-only/local.
    vector<PafChainedCandidate> mergedResults;
    {
        size_t total = 0;
        for(const auto& v : threadResults) {
            total += v.size();
        }
        mergedResults.reserve(total);
        for(auto& v : threadResults) {
            for(auto& r : v) {
                mergedResults.push_back(std::move(r));
            }
        }
    }

    // Apply Hifiasm-style max_n_chain filtering in the PAF path.
    // We perform this per query read, using the same per-overlap-type thresholds
    // and COV_W-style weak keep policy used by discovery chaining.
    if(maxChainLimit > 0 && !mergedResults.empty()) {
        std::sort(mergedResults.begin(), mergedResults.end(),
            [](const PafChainedCandidate& a, const PafChainedCandidate& b) {
                if(a.queryReadId != b.queryReadId) {
                    return a.queryReadId < b.queryReadId;
                }
                if(a.score != b.score) {
                    return a.score > b.score;
                }
                if(a.overlapType != b.overlapType) {
                    return a.overlapType < b.overlapType;
                }
                if(a.candidate.readIds[0] != b.candidate.readIds[0]) {
                    return a.candidate.readIds[0] < b.candidate.readIds[0];
                }
                return a.candidate.readIds[1] < b.candidate.readIds[1];
            });

        vector<PafChainedCandidate> filteredResults;
        filteredResults.reserve(mergedResults.size());
        size_t begin = 0;
        while(begin < mergedResults.size()) {
            size_t end = begin + 1;
            while(end < mergedResults.size() &&
                mergedResults[end].queryReadId == mergedResults[begin].queryReadId) {
                ++end;
            }

            if(end - begin <= maxChainLimit) {
                for(size_t i = begin; i < end; ++i) {
                    filteredResults.push_back(std::move(mergedResults[i]));
                }
                begin = end;
                continue;
            }

            uint64_t countByType[4] = {0, 0, 0, 0};
            int32_t thresholdByType[4] = {0, 0, 0, 0};
            for(size_t i = begin; i < end; ++i) {
                const uint32_t type = std::min<uint32_t>(3, mergedResults[i].overlapType);
                countByType[type]++;
                if(countByType[type] == maxChainLimit) {
                    thresholdByType[type] = mergedResults[i].score;
                }
            }

            const uint64_t readLenA = reads->getReadRawSequenceLength(mergedResults[begin].queryReadId);
            const uint32_t readLenA32 = uint32_t(std::min<uint64_t>(
                readLenA,
                uint64_t(std::numeric_limits<uint32_t>::max())));
            const bool useOcvWeakKeep =
                (mcopyOcvWindow > 0) &&
                (readLenA32 >= mcopyOcvWindow) &&
                (countByType[3] >= maxChainLimit);
            vector<uint64_t> ocvWindowUsage;
            if(useOcvWeakKeep) {
                const uint32_t windowCount = (readLenA32 + mcopyOcvWindow - 1U) / mcopyOcvWindow;
                ocvWindowUsage.assign(windowCount, uint64_t(0));
                for(uint32_t w = 0; w < windowCount; ++w) {
                    const uint64_t ws = uint64_t(w) * uint64_t(mcopyOcvWindow);
                    const uint64_t we = std::min<uint64_t>(ws + uint64_t(mcopyOcvWindow), uint64_t(readLenA32));
                    uint64_t cap = (we - ws) * uint64_t(maxChainLimit >> 1);
                    if(cap > uint64_t(std::numeric_limits<uint32_t>::max())) {
                        cap = uint64_t(std::numeric_limits<uint32_t>::max());
                    }
                    ocvWindowUsage[w] = (cap << 32);
                }
            }

            auto updateOcvWindows = [&](const PafChainedCandidate& e) {
                if(!useOcvWeakKeep) {
                    return;
                }
                const uint32_t qs = std::min<uint32_t>(readLenA32, e.alignment.qs);
                const uint32_t qe = std::min<uint32_t>(readLenA32, e.alignment.qe);
                if(qe <= qs) {
                    return;
                }
                for(uint32_t w = qs / mcopyOcvWindow; w < ocvWindowUsage.size(); ++w) {
                    const uint32_t ws = w * mcopyOcvWindow;
                    const uint32_t we = std::min<uint32_t>(ws + mcopyOcvWindow, readLenA32);
                    const uint32_t os = std::max<uint32_t>(qs, ws);
                    const uint32_t oe = std::min<uint32_t>(qe, we);
                    if(oe <= os) {
                        if(ws >= qe) {
                            break;
                        }
                        continue;
                    }
                    const uint32_t overlap = oe - os;
                    uint64_t v = ocvWindowUsage[w];
                    const uint32_t cap = uint32_t(v >> 32);
                    uint32_t used = uint32_t(v & uint64_t(std::numeric_limits<uint32_t>::max()));
                    const uint64_t sum = uint64_t(used) + uint64_t(overlap);
                    used = (sum > uint64_t(std::numeric_limits<uint32_t>::max())) ?
                        std::numeric_limits<uint32_t>::max() : uint32_t(sum);
                    ocvWindowUsage[w] = (uint64_t(cap) << 32) | uint64_t(used);
                    if(we >= qe) {
                        break;
                    }
                }
            };

            auto evaluateOcvWeakKeep = [&](const PafChainedCandidate& e) -> bool {
                if(!useOcvWeakKeep) {
                    return false;
                }
                const uint32_t qs = std::min<uint32_t>(readLenA32, e.alignment.qs);
                const uint32_t qe = std::min<uint32_t>(readLenA32, e.alignment.qe);
                if(qe <= qs) {
                    return false;
                }
                uint64_t good = 0;
                uint64_t bad = 0;
                for(uint32_t w = qs / mcopyOcvWindow; w < ocvWindowUsage.size(); ++w) {
                    const uint32_t ws = w * mcopyOcvWindow;
                    const uint32_t we = std::min<uint32_t>(ws + mcopyOcvWindow, readLenA32);
                    const uint32_t os = std::max<uint32_t>(qs, ws);
                    const uint32_t oe = std::min<uint32_t>(qe, we);
                    if(oe <= os) {
                        if(ws >= qe) {
                            break;
                        }
                        continue;
                    }
                    const uint32_t overlap = oe - os;
                    const uint64_t v = ocvWindowUsage[w];
                    const uint32_t cap = uint32_t(v >> 32);
                    const uint32_t used = uint32_t(v & uint64_t(std::numeric_limits<uint32_t>::max()));
                    if(uint64_t(used) + uint64_t(overlap) >= uint64_t(cap)) {
                        bad += uint64_t(overlap);
                    } else {
                        good += uint64_t(overlap);
                    }
                    if(we >= qe) {
                        break;
                    }
                }
                const uint64_t total = good + bad;
                if(total == 0) {
                    return false;
                }
                return double(good) >= (double(total) * mcopyOcvWeakKeepRatio);
            };

            for(size_t i = begin; i < end; ++i) {
                const uint32_t type = std::min<uint32_t>(3, mergedResults[i].overlapType);
                bool keep = (mergedResults[i].score >= thresholdByType[type]);
                if(!keep && type == 3U) {
                    keep = evaluateOcvWeakKeep(mergedResults[i]);
                }
                if(keep) {
                    updateOcvWindows(mergedResults[i]);
                    filteredResults.push_back(std::move(mergedResults[i]));
                }
            }

            begin = end;
        }
        mergedResults.swap(filteredResults);
    }

    for(auto& r : mergedResults) {
        alignmentCandidates.candidates.push_back(r.candidate);
        alignmentCandidatesAlignmentsData.alignments.push_back(std::move(r.alignment));
    }

    alignmentCandidates.candidates.unreserve();
    alignmentCandidatesAlignmentsData.alignments.unreserve();

    // Cleanup inverted index data
    clearInvertedIndexTransientData(invertedIndexData);

    const auto endTime = std::chrono::steady_clock::now();
    const double totalSeconds = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime)).count());
    cout << timestamp << "PAF candidate chaining completed in " << totalSeconds << " s." << endl;
    cout << timestamp << "Chained " << alignmentCandidates.candidates.size() << " candidates." << endl;
}
