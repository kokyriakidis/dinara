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
#include <array>
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
    uint32_t weight;       ///< Frequency-based weight (hifiasm k_mer_hit.cnt semantics).
    // Orientation bits relative to the canonical k-mer id:
    // - isRcA=1 means the observed k-mer on A is the RC of the canonical id.
    // - isRcB=1 means the observed k-mer on B is the RC of the canonical id.
    // In hifiasm terms, these correspond to z->rev and y->rev, and the per-hit
    // strand bucket is rev = isRcA ^ isRcB.
    uint8_t isRcA = 0;
    uint8_t isRcB = 0;
    
    /// Comparison operator for sorting by (partnerReadId, posA).
    bool operator<(const InvertedIndexTempHit& other) const {
        if (partnerReadId != other.partnerReadId) return partnerReadId < other.partnerReadId;
        return posA < other.posA;
    }
};

// Radix-sort `flatHits` by (partnerReadId, posA) using a packed 64-bit key.
// This replaces `std::sort` in the hot path; hifiasm uses radix sorts for the same reason.
static inline uint64_t flatHitPackedKey(const InvertedIndexTempHit& h)
{
    return (uint64_t(h.partnerReadId) << 32) | uint64_t(h.posA);
}

static void radixSortFlatHitsByPartnerReadIdAndPosA(
    vector<InvertedIndexTempHit>& hits,
    vector<InvertedIndexTempHit>& tmp)
{
    if(hits.size() <= 1) {
        return;
    }
    tmp.resize(hits.size());

    // LSD radix sort, 8 passes over the 64-bit packed key.
    // Stable counting-sort per byte ensures correct lexicographic order.
    for(int pass = 0; pass < 8; ++pass) {
        const uint32_t shift = uint32_t(pass * 8);
        std::array<size_t, 256> count{};
        for(const auto& h : hits) {
            const uint8_t b = uint8_t((flatHitPackedKey(h) >> shift) & 0xffULL);
            ++count[b];
        }
        size_t sum = 0;
        for(size_t i = 0; i < count.size(); ++i) {
            const size_t c = count[i];
            count[i] = sum;
            sum += c;
        }
        for(const auto& h : hits) {
            const uint8_t b = uint8_t((flatHitPackedKey(h) >> shift) & 0xffULL);
            tmp[count[b]++] = h;
        }
        hits.swap(tmp);
    }
}

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
// Hifiasm's `ha_ov_type()` classifies using inclusive coordinates (x_pos_s, x_pos_e) and len-1.
// Dinara stores alignment coordinates as half-open intervals [start, end), so we match semantics by
// comparing against 0 and readLen (end-exclusive).
static int getOverlapType(uint32_t start, uint32_t end, uint32_t readLen) {
    if (start == 0 && end >= readLen) return 2;        // contained in a longer read
    else if (start > 0 && end < readLen) return 3;     // containing a shorter read
    else return (start == 0) ? 0 : 1;                  // left or right overhang
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
	    vector<InvertedIndexTempHit> flatHitsTmp;
    
    // Structure of Arrays (SoA) for cache-efficient DP scans
    vector<uint32_t> hitPosA, hitPosB, hitOrdinalA, hitOrdinalB;
    vector<uint32_t> hitOrderByPosB;
    vector<uint32_t> hitWeights;

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
	        flatHitsTmp.clear();
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

// hifiasm's yak_hash64_64 (htab.h:150) used during minimizer sketching.
// We use it as a deterministic tie-break key when downsampling high-frequency marker streaks.
static inline uint64_t hifiasmYakHash64_64(uint64_t key)
{
    key = ~key + (key << 21);
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8);
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4);
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

static inline uint64_t foldKmerIdToUint64(const KmerId k)
{
    // Stable 128/256-bit -> 64-bit fold for ordering/tie-breaking only.
    const uint64_t lo = uint64_t(k);
    const uint64_t hi = uint64_t(k >> 64);
    return lo ^ hi;
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
    const vector<uint32_t>& hitWeights,
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

// ============================================================================
// STRICT HIFIASM PORT: lchain_qdp_mcopy_fast + quick_ck_lchain
// ============================================================================
// Source of truth: hifiasm 0.25.0-r726 (ec9a8b2)
// - Hash_Table.cpp:2007-2095 (quick_ck_lchain)
// - Hash_Table.cpp:2097-2284 (lchain_qdp_mcopy_fast)
// - Hash_Table.cpp:1475-1541 (cal_bw + comput_sc_ch_ec)
// - Hash_Table.cpp:779-809 (get_chainLen)
// - Hash_Table.cpp:1752-1780 (push_ovlp_chain_qgen)
//
// Notes:
// - This port is intentionally "mechanical": variable names and control flow match hifiasm
//   so behavior is step-parity given identical hit sets.
// - Dinara stores intervals as half-open, but hifiasm overlap_region uses inclusive end coords.
//   The port keeps inclusive ends internally and converts when building `Alignment`.
// ============================================================================
namespace {

static inline int32_t hifiasm_normal_w(int32_t x, int32_t y)
{
    // Hash_Table.cpp:20
    return (x >= y) ? (x / y) : 1;
}

struct HifiasmKmerHit {
    // Minimal k_mer_hit fields used by lchain DP.
    uint32_t readID = 0;       // partner read id (y_id)
    uint32_t self_offset = 0;  // x (query) "end" coordinate
    uint32_t offset = 0;       // y (target) "end" coordinate (already transformed for strand=1)
    uint32_t cnt = 0;          // (weight<<8) | span (span in low 8 bits)
    uint8_t strand = 0;        // 0=same, 1=diff (rev)

    // Extra payload used to build `Alignment` (ignored by DP).
    uint32_t ordinalA = 0;
    uint32_t ordinalB = 0;
    // Global index into a caller-owned array so chain-hit lists can outlive the per-pair `a` vector.
    // This mirrors how hifiasm uses `non_homopolymer_errors` as an offset into a global hit buffer.
    uint32_t globalIndex = 0;
};

struct HifiasmOverlapRegion {
    // Mirrors overlap_region fields used downstream in anchor.cpp.
    uint32_t x_id = 0;
    uint32_t y_id = 0;
    uint8_t x_pos_strand = 0;
    uint8_t y_pos_strand = 0;
    uint32_t x_pos_s = 0;  // inclusive
    uint32_t x_pos_e = 0;  // inclusive
    uint32_t y_pos_s = 0;  // inclusive (RC-oriented coordinates if y_pos_strand==1)
    uint32_t y_pos_e = 0;  // inclusive
    int32_t shared_seed = 0; // chaining score
    uint32_t align_length = 0; // number of hits in chain

    // Used by hifiasm as an offset into the chain-hit array; Dinara uses it as an offset
    // into `chainHitIndexFlat` (see below).
    uint32_t non_homopolymer_errors = 0;

    // Dinara-only convenience: unextended endpoints (half-open) for optional maxEndFuzz filtering.
    uint32_t raw_xs = 0, raw_xe = 0;
    uint32_t raw_ys = 0, raw_ye = 0;
};

struct HifiasmChainDataScratch {
    // Mirrors Chain_Data fields used by lchain_qdp_mcopy_fast.
    vector<int64_t> pre;   // dp->pre (backtrack pointers)
    vector<int64_t> tmp;   // dp->tmp (visit marks + reuse as heap/keys)
    vector<int32_t> score; // dp->score
    vector<int32_t> occ;   // dp->occ (reused as index workspace)

    void resize(size_t n)
    {
        pre.resize(n);
        tmp.resize(n);
        score.resize(n);
        occ.resize(n);
    }
};

static inline int32_t hifiasm_cal_bw(
    const HifiasmKmerHit* ai,
    const HifiasmKmerHit* aj,
    double bw_rate,
    int64_t sf_l,
    int64_t ot_l)
{
    // Hash_Table.cpp:1475-1488
    int64_t sf_s = int64_t(aj->self_offset);
    int64_t sf_e = int64_t(ai->self_offset) + 1;
    int64_t ot_s = int64_t(aj->offset);
    int64_t ot_e = int64_t(ai->offset) + 1;
    int64_t sf_r = sf_l - sf_e;
    int64_t ot_r = ot_l - ot_e;
    if (sf_s <= ot_s) sf_s = 0;
    else sf_s -= ot_s;
    if (sf_r <= ot_r) sf_e = sf_l;
    else sf_e += ot_r;
    return int32_t(double(sf_e - sf_s) * bw_rate);
}

static inline int32_t hifiasm_comput_sc_ch_ec(
    const HifiasmKmerHit* ai,
    const HifiasmKmerHit* aj,
    double bw_rate,
    double chn_pen_gap,
    double chn_pen_skip,
    int64_t sl,
    int64_t ol)
{
    // Hash_Table.cpp:1515-1541
    int32_t dq = int32_t(int64_t(ai->self_offset) - int64_t(aj->self_offset));
    if (dq <= 0) return INT32_MIN;
    int32_t dr = int32_t(int64_t(ai->offset) - int64_t(aj->offset));
    if (dr <= 0) return INT32_MIN;
    int32_t dd = (dr > dq) ? (dr - dq) : (dq - dr);
    if ((dd > 16) && (dd > hifiasm_cal_bw(ai, aj, bw_rate, sl, ol))) return INT32_MIN;
    int32_t dg = (dr < dq) ? dr : dq;
    int32_t q_span = int32_t(ai->cnt & 0xffu);
    int32_t sc = (q_span < dg) ? q_span : dg;
    sc = hifiasm_normal_w(sc, int32_t(ai->cnt >> 8));
    if (dd || (dg > q_span && dg > 0)) {
        double lin_pen = chn_pen_gap * double(dd);
        double a_pen = double(sc) * (double(dd) / double(dg)) / bw_rate;
        if (dd < 4) lin_pen = (lin_pen > a_pen) ? a_pen : lin_pen;
        else lin_pen = (lin_pen < a_pen) ? a_pen : lin_pen;
        lin_pen += chn_pen_skip * double(dg);
        sc -= int32_t(lin_pen);
    }
    return sc;
}

static inline int64_t hifiasm_get_chainLen(
    int64_t x_beg, int64_t x_end, int64_t xLen,
    int64_t y_beg, int64_t y_end, int64_t yLen)
{
    // Hash_Table.cpp:779-809
    if (x_beg <= y_beg) { y_beg = y_beg - x_beg; x_beg = 0; }
    else { x_beg = x_beg - y_beg; y_beg = 0; }

    int64_t x_right_length = xLen - x_end - 1;
    int64_t y_right_length = yLen - y_end - 1;
    if (x_right_length <= y_right_length) {
        x_end = xLen - 1;
        y_end = y_end + x_right_length;
    } else {
        x_end = x_end + y_right_length;
        y_end = yLen - 1;
    }
    return x_end - x_beg + 1;
}

static inline void hifiasm_push_ovlp_chain_qgen(
    HifiasmOverlapRegion& o,
    uint32_t xid,
    int64_t xl,
    int64_t yl,
    int64_t sc,
    const HifiasmKmerHit* beg,
    const HifiasmKmerHit* end)
{
    // Hash_Table.cpp:1752-1780
    o.x_id = xid;
    o.y_id = beg->readID;
    o.x_pos_strand = 0;
    o.y_pos_strand = beg->strand;
    o.x_pos_s = beg->self_offset;
    o.y_pos_s = beg->offset;
    o.x_pos_e = end->self_offset;
    o.y_pos_e = end->offset;

    // Save unextended endpoints for optional maxEndFuzz filtering.
    o.raw_xs = o.x_pos_s;
    o.raw_xe = o.x_pos_e + 1U;
    o.raw_ys = o.y_pos_s;
    o.raw_ye = o.y_pos_e + 1U;

    if (o.x_pos_s <= o.y_pos_s) { o.y_pos_s -= o.x_pos_s; o.x_pos_s = 0; }
    else { o.x_pos_s -= o.y_pos_s; o.y_pos_s = 0; }

    const int64_t xr = xl - int64_t(o.x_pos_e) - 1;
    const int64_t yr = yl - int64_t(o.y_pos_e) - 1;
    if (xr <= yr) { o.x_pos_e = uint32_t(xl - 1); o.y_pos_e += uint32_t(xr); }
    else { o.y_pos_e = uint32_t(yl - 1); o.x_pos_e += uint32_t(yr); }

    o.shared_seed = int32_t(sc);
    o.align_length = 0;
    o.non_homopolymer_errors = 0;
}

static inline void hifiasm_quick_ck_lchain(
    HifiasmKmerHit* a,
    int64_t a_n,
    int64_t xl,
    int64_t yl,
    double chn_pen_gap,
    double chn_pen_skip,
    double bw_rate,
    int64_t* p,
    int64_t* t,
    int32_t* f,
    int32_t* ii,
    int64_t* plus,
    int64_t* msc,
    int64_t* msc_i,
    int64_t* movl,
    int64_t* si,
    int64_t* ei)
{
    // Hash_Table.cpp:2007-2095
    if (a_n <= 0) return;
    int64_t l, k, is_srt = 1, z;
    HifiasmKmerHit *ai, *aj;
    int64_t dq, dr, dd, dg, q_span, sc, csc, ddt;
    int64_t plus0, msc0, msc_i0, movl0;
    double lin_pen, a_pen;

    *plus = 0;
    *msc = *msc_i = INT32_MIN;
    *movl = INT32_MAX;
    *si = 0;
    *ei = a_n;

    for (k = 1, l = 0; k <= a_n; k++) {
        if (k == a_n || a[k].strand != a[l].strand) {
            t[k - 1] = 0;
            ii[k - 1] = 0;
            if (is_srt) {
                plus0 = 0;
                msc0 = msc_i0 = INT32_MIN;
                movl0 = INT32_MAX;
                ddt = 0;

                p[l] = -1;
                f[l] = int32_t(a[l].cnt & 0xffu);
                if (f[l] >= msc0) { msc0 = f[l]; msc_i0 = l; }
                if (f[l] < plus0) plus0 = f[l];

                for (z = l + 1; z < k; z++) {
                    ai = &a[z];
                    aj = &a[z - 1];
                    dq = int64_t(ai->self_offset) - int64_t(aj->self_offset);
                    if (dq <= 0) break;
                    dr = int64_t(ai->offset) - int64_t(aj->offset);
                    if (dr <= 0) break;
                    dd = (dr > dq) ? (dr - dq) : (dq - dr);
                    if ((dd > 16) && (dd > hifiasm_cal_bw(&a[z], &a[z - 1], bw_rate, xl, yl))) break;
                    dg = (dr < dq) ? dr : dq;
                    q_span = int64_t(ai->cnt & 0xffu);
                    sc = (q_span < dg) ? q_span : dg;
                    sc = hifiasm_normal_w(int32_t(sc), int32_t(ai->cnt >> 8));
                    if (dd || (dg > q_span && dg > 0)) {
                        lin_pen = chn_pen_gap * double(dd);
                        a_pen = double(sc) * (double(dd) / double(dg)) / bw_rate;
                        if (dd < 4) lin_pen = (lin_pen > a_pen) ? a_pen : lin_pen;
                        else lin_pen = (lin_pen < a_pen) ? a_pen : lin_pen;
                        lin_pen += chn_pen_skip * double(dg);
                        sc -= int32_t(lin_pen);
                    }

                    sc += f[z - 1];
                    csc = int64_t(a[z].cnt & 0xffu);
                    if (sc < csc) break;
                    p[z] = z - 1;
                    f[z] = int32_t(sc);
                    ddt += dd;
                    if (f[z] >= msc0) { msc0 = f[z]; msc_i0 = z; }
                    if (f[z] < plus0) plus0 = f[z];
                }

                if ((z >= k) && (msc_i0 == (k - 1))) {
                    if ((k - l >= 2) && (ddt > 16) && (ddt > hifiasm_cal_bw(&a[k - 1], &a[l], bw_rate, xl, yl))) {
                        msc_i0 = INT32_MIN;
                    }
                    if (msc_i0 == (k - 1)) {
                        if (msc0 >= (*msc)) {
                            movl0 = hifiasm_get_chainLen(
                                a[msc_i0].self_offset, a[msc_i0].self_offset, xl,
                                a[msc_i0].offset, a[msc_i0].offset, yl);
                            if (msc0 > (*msc) || movl0 < (*movl)) {
                                *msc = msc0;
                                *msc_i = msc_i0;
                                *movl = movl0;
                            }
                        }
                        if (plus0 < (*plus)) *plus = plus0;
                        if ((*ei) > k) (*si) = k;
                        else (*ei) = l;
                    }
                }
            }
            l = k;
            is_srt = 1;
        } else {
            if ((a[k].self_offset <= a[k - 1].self_offset) || (a[k].offset <= a[k - 1].offset)) is_srt = 0;
            t[k - 1] = 0;
            ii[k - 1] = 0;
        }
    }
}

// Port of Hash_Table.cpp:2097-2284 (lchain_qdp_mcopy_fast), adapted to:
// - take `a` as a vector
// - emit overlaps + chain-hit indices into flat storage
static inline void hifiasm_lchain_qdp_mcopy_fast(
    vector<HifiasmKmerHit>& a,
    HifiasmChainDataScratch& dp,
    vector<HifiasmOverlapRegion>& res,
    vector<uint32_t>& chainHitIndexFlat,
    int64_t max_skip,
    int64_t max_iter,
    int64_t max_dis,
    double chn_pen_gap,
    double chn_pen_skip,
    double bw_rate,
    uint32_t xid,
    int64_t xl,
    int64_t yl,
    int64_t quick_check,
    int64_t mcopy_num,
    double mcopy_rate,
    int64_t mcopy_khit_cutoff)
{
    const int64_t a_n = int64_t(a.size());
    if (a_n <= 0) return;
    dp.resize(size_t(a_n));

    int64_t* p = dp.pre.data();
    int64_t* t = dp.tmp.data();
    int32_t* f = dp.score.data();
    int32_t* ii = dp.occ.data();

    int64_t max_f, n_skip, st, max_j, end_j, sc, msc, msc_i, max_ii, ovl, movl, plus = 0, min_sc, ch_n, si, ei;
    int32_t max, tmp;
    int64_t i, k, j, cL = 0;

	    if (quick_check) {
	        hifiasm_quick_ck_lchain(a.data(), a_n, xl, yl, chn_pen_gap, chn_pen_skip, bw_rate,
	            p, t, f, ii, &plus, &msc, &msc_i, &movl, &si, &ei);
	    } else {
	        msc = msc_i = INT32_MIN;
	        movl = INT32_MAX;
	        plus = 0;
	        si = 0;
	        ei = a_n;
	        std::fill(dp.tmp.begin(), dp.tmp.end(), int64_t(0));
	    }

	    // Safety (parity-preserving): if quick_ck_lchain didn't produce a usable best chain, fall back
	    // to full DP over the entire hit range so we always have a valid (msc_i, p[]) chain to emit.
	    // In hifiasm, msc_i should never stay <0 when a_n>0; if it does, it indicates a mismatch
	    // between our hit ordering/coordinate conventions and quick_ck_lchain's assumptions.
	    if (msc_i < 0 || si < 0 || ei < 0 || si > a_n || ei > a_n) {
	        msc = msc_i = INT32_MIN;
	        movl = INT32_MAX;
	        plus = 0;
	        si = 0;
	        ei = a_n;
	        std::fill(dp.tmp.begin(), dp.tmp.end(), int64_t(0));
	    }

    for (i = st = si, max_ii = -1; i < ei; ++i) {
        max_f = int64_t(a[size_t(i)].cnt & 0xffu);
        n_skip = 0;
        max_j = end_j = -1;
        if ((i - st) > max_iter) st = i - max_iter;
        while (a[size_t(i)].strand != a[size_t(st)].strand) ++st;

        for (j = i - 1; j >= st; --j) {
            sc = hifiasm_comput_sc_ch_ec(&a[size_t(i)], &a[size_t(j)], bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
            if (sc == INT32_MIN) continue;
            sc += f[size_t(j)];
            if (sc > max_f) {
                max_f = sc;
                max_j = j;
                if (n_skip > 0) --n_skip;
            } else if (t[size_t(j)] == (int32_t)i) {
                if (++n_skip > max_skip) break;
            }
            if (p[size_t(j)] >= 0) t[size_t(p[size_t(j)])] = i;
        }
        end_j = j;

        if ((max_ii < 0) ||
            (a[size_t(i)].self_offset > a[size_t(max_ii)].self_offset + max_dis) ||
            (a[size_t(i)].strand != a[size_t(max_ii)].strand)) {
            max = INT32_MIN;
            max_ii = -1;
            for (j = i - 1;
                 (j >= st) &&
                 (a[size_t(i)].self_offset <= max_dis + a[size_t(j)].self_offset) &&
                 (a[size_t(i)].strand == a[size_t(j)].strand);
                 --j) {
                if (max < f[size_t(j)]) {
                    max = f[size_t(j)];
                    max_ii = j;
                }
            }
        }

        if ((max_ii >= 0) && (max_ii < end_j) && (a[size_t(i)].strand == a[size_t(max_ii)].strand)) {
            tmp = hifiasm_comput_sc_ch_ec(&a[size_t(i)], &a[size_t(max_ii)], bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
            if (tmp != INT32_MIN && max_f < tmp + f[size_t(max_ii)]) {
                max_f = tmp + f[size_t(max_ii)];
                max_j = max_ii;
            }
        }
        f[size_t(i)] = int32_t(max_f);
        p[size_t(i)] = max_j;
        if ((max_ii < 0) ||
            ((a[size_t(i)].self_offset <= max_dis + a[size_t(max_ii)].self_offset) &&
             (a[size_t(i)].strand == a[size_t(max_ii)].strand) &&
             (f[size_t(max_ii)] < f[size_t(i)]))) {
            max_ii = i;
        }
        if (f[size_t(i)] >= msc) {
            ovl = hifiasm_get_chainLen(
                a[size_t(i)].self_offset, a[size_t(i)].self_offset, xl,
                a[size_t(i)].offset, a[size_t(i)].offset, yl);
            if (f[size_t(i)] > msc || ovl < movl) {
                msc = f[size_t(i)];
                msc_i = i;
                movl = ovl;
            }
        }
        if (f[size_t(i)] < plus) plus = f[size_t(i)];
        ii[size_t(i)] = 0;
    }

	    for (i = msc_i, cL = 0; i >= 0; i = p[size_t(i)]) {
	        ii[size_t(i)] = 1;
	        t[size_t(cL++)] = i;
	    }

    // mcopy selection (Hash_Table.cpp:2180+)
    if (mcopy_num > 1) {
        if (cL >= mcopy_khit_cutoff) {
            msc -= plus;
            min_sc = int64_t(double(msc) * mcopy_rate);
            ii[size_t(msc_i)] = 0;
	            for (i = ch_n = 0; i < a_n; ++i) {
	                f[size_t(i)] -= int32_t(plus);
	                if (i >= ch_n) t[size_t(i)] = 0;
	                if ((!ii[size_t(i)]) && (int64_t(f[size_t(i)]) >= min_sc)) {
	                    t[size_t(ch_n)] = (int64_t(uint64_t(uint32_t(f[size_t(i)])) << 32) | (uint64_t(i) << 1));
	                    ch_n++;
	                }
	            }
	            // If we ended up with <=1 candidate, the selection loop above overwrote `t[]` (which used
	            // to hold the best-chain backtrack indices). Hifiasm's downstream code assumes `t[]`
	            // still contains the best chain, so restore it here to avoid out-of-bounds access.
	            if (ch_n <= 1) {
	                msc += plus;
	                i = msc_i;
	                cL = 0;
	                while (i >= 0) { t[size_t(cL++)] = i; i = p[size_t(i)]; }
	            }
	            if (ch_n > 1) {
	                int64_t n_v, n_v0, ni, n_u, n_u0 = int64_t(res.size());
	                std::sort(t, t + ch_n,
	                    [](const int64_t a, const int64_t b) {
	                        return uint64_t(a) < uint64_t(b);
	                    });
	                for (k = ch_n - 1, n_v = n_u = 0; k >= 0 && n_u < mcopy_num; --k) {
	                    n_v0 = n_v;
	                    for (i = (int64_t)(((uint32_t)t[size_t(k)]) >> 1); i >= 0 && (t[size_t(i)] & 1) == 0; ) {
	                        ii[size_t(n_v++)] = int32_t(i);
	                        t[size_t(i)] |= 1;
	                        i = p[size_t(i)];
	                    }
	                    if (n_v0 == n_v) continue;
	                    const uint64_t key = uint64_t(t[size_t(k)]);
	                    const int64_t top = int64_t(uint32_t(key >> 32));
	                    sc = (i < 0) ? top : (top - int64_t(f[size_t(i)]));
	                    if (sc >= min_sc) {
                        HifiasmOverlapRegion z{};
                        hifiasm_push_ovlp_chain_qgen(z, xid, xl, yl, sc + plus,
                            &a[size_t(ii[size_t(n_v - 1)])], &a[size_t(ii[size_t(n_v0)])]);
                        if ((!n_u) || (n_v - n_v0 > 1)) {
                            z.align_length = uint32_t(n_v - n_v0);
                            z.x_id = uint32_t(n_v0); // temp: start in ii[]
                            // Store chain-hit indices in increasing self_offset order.
                            z.non_homopolymer_errors = uint32_t(chainHitIndexFlat.size());
                            ni = z.align_length;
                            for (j = 0; j < ni; ++j) {
                                const int64_t local = ii[size_t(n_v0 + (ni - j - 1))];
                                chainHitIndexFlat.push_back(a[size_t(local)].globalIndex);
                            }
                            res.push_back(z);
                            n_u++;
                        } else {
                            // non-best is tiny
                            n_v = n_v0;
                        }
                    } else {
                        n_v = n_v0;
                    }
                }
                if (res.size() > size_t(n_u0)) {
                    // If we emitted any mcopy overlaps, we are done.
                    return;
                }

                // Fall back to best chain if nothing survived.
                msc += plus;
                i = msc_i;
                cL = 0;
                while (i >= 0) { t[size_t(cL++)] = i; i = p[size_t(i)]; }
            }
        }
    }

    // Best chain emission (Hash_Table.cpp:2274+)
    HifiasmOverlapRegion z{};
    hifiasm_push_ovlp_chain_qgen(z, xid, xl, yl, msc, &a[size_t(t[size_t(cL - 1)])], &a[size_t(t[0])]);
    z.align_length = uint32_t(cL);
    z.non_homopolymer_errors = uint32_t(chainHitIndexFlat.size());
    for (i = 0; i < cL; ++i) {
        const int64_t local = t[size_t(cL - i - 1)];
        chainHitIndexFlat.push_back(a[size_t(local)].globalIndex);
    }
    res.push_back(z);
}

static inline int hifiasm_ha_ov_type(const HifiasmOverlapRegion& r, const uint32_t len)
{
    // anchor.cpp:86
    if (r.x_pos_s == 0 && r.x_pos_e == len - 1U) return 2;        // contained in a longer read
    else if (r.x_pos_s > 0 && r.x_pos_e < len - 1U) return 3;     // containing a shorter read
    else return (r.x_pos_s == 0) ? 0 : 1;
}

// Strict port of hifiasm anchor.cpp:lchain_qgen_mcopy_fast max_n_chain + ocv_w + r485 logic,
// adapted to Dinara's in-memory vectors.
static inline void hifiasm_lchain_qgen_mcopy_fast_postfilter(
    vector<HifiasmOverlapRegion>& ol,
    const uint64_t max_n_chain,
    const uint32_t chain_cutoff,
    const uint64_t ocv_w,
    const uint32_t rl,
    const vector<uint32_t>& chainHitIndexFlat,
    const vector<HifiasmKmerHit>& allHits)
{
    if (ol.empty()) {
        return;
    }

    uint64_t lch = 0;

    if (max_n_chain > 0 && ol.size() > max_n_chain) {
        // ks_introsort_or_ss: sort only by shared_seed descending.
        std::sort(ol.begin(), ol.end(),
            [](const HifiasmOverlapRegion& a, const HifiasmOverlapRegion& b) {
                return a.shared_seed > b.shared_seed;
            });

        int32_t n[4] = {0, 0, 0, 0};
        int32_t s[4] = {0, 0, 0, 0};
        for (size_t i = 0; i < ol.size(); ++i) {
            const int w = hifiasm_ha_ov_type(ol[i], rl);
            ++n[w];
            if (uint64_t(n[w]) == max_n_chain) s[w] = ol[i].shared_seed;
        }

        if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
            // ocv_w windows (COV_W) for type-3 overload control.
            vector<uint64_t> cc;
            uint64_t cwn = 0;
            if (ocv_w > 0 && (uint64_t(n[3]) >= max_n_chain) && (uint64_t(rl) >= ocv_w)) {
                cwn = (uint64_t(rl) / ocv_w) + ((uint64_t(rl) % ocv_w) ? 1ULL : 0ULL);
                cc.resize(size_t(cwn), uint64_t(0));
                for (uint64_t i = 0, cws = 0, cwe = 0; i < cwn; ++i) {
                    cwe = cws + ocv_w;
                    if (cwe > uint64_t(rl)) cwe = uint64_t(rl);
                    const uint64_t winLen = cwe - cws;
                    uint64_t cap = winLen * (max_n_chain >> 1);
                    if (cap > uint64_t(UINT32_MAX)) cap = uint64_t(UINT32_MAX);
                    cc[size_t(i)] = (cap << 32); // high32=cap, low32=used
                    cws += ocv_w;
                }
            }

            auto update_cc = [&](const HifiasmOverlapRegion& r) {
                if (cwn == 0) return;
                uint64_t m = uint64_t(r.x_pos_s) / ocv_w;
                const uint64_t rs = uint64_t(r.x_pos_s);
                const uint64_t re = uint64_t(r.x_pos_e) + 1ULL;
                for (uint64_t cws = m * ocv_w; m < cwn; ++m, cws += ocv_w) {
                    uint64_t cwe = cws + ocv_w;
                    if (cwe > uint64_t(rl)) cwe = uint64_t(rl);
                    const uint64_t os = (rs >= cws) ? rs : cws;
                    const uint64_t oe = (re <= cwe) ? re : cwe;
                    if (oe <= os) break;
                    const uint32_t add = uint32_t(oe - os);
                    const uint32_t used = uint32_t(cc[size_t(m)] & 0xffffffffULL);
                    if (uint64_t(used) + uint64_t(add) < uint64_t(UINT32_MAX)) {
                        cc[size_t(m)] += uint64_t(add);
                    } else {
                        cc[size_t(m)] = (cc[size_t(m)] & 0xffffffff00000000ULL) | uint64_t(UINT32_MAX);
                    }
                }
            };

            auto should_rescue_type3 = [&](const HifiasmOverlapRegion& r) -> bool {
                if (cwn == 0) return false;
                uint64_t m = uint64_t(r.x_pos_s) / ocv_w;
                uint64_t cw0 = 0, cw1 = 0;
                const uint64_t rs = uint64_t(r.x_pos_s);
                const uint64_t re = uint64_t(r.x_pos_e) + 1ULL;
                for (uint64_t cws = m * ocv_w; m < cwn; ++m, cws += ocv_w) {
                    uint64_t cwe = cws + ocv_w;
                    if (cwe > uint64_t(rl)) cwe = uint64_t(rl);
                    const uint64_t os = (rs >= cws) ? rs : cws;
                    const uint64_t oe = (re <= cwe) ? re : cwe;
                    if (oe <= os) break;
                    const uint64_t overlap = oe - os;
                    const uint64_t cap = (cc[size_t(m)] >> 32);
                    const uint64_t used = uint64_t(uint32_t(cc[size_t(m)] & 0xffffffffULL));
                    if (used + overlap >= cap) cw1 += overlap;
                    else cw0 += overlap;
                }
                const uint64_t total = cw0 + cw1;
                if (total == 0) return false;
                // anchor.cpp:2027 uses a hard-coded 0.7 threshold.
                return double(cw0) >= (double(total) * 0.7);
            };

            size_t k = 0;
            for (size_t i = 0; i < ol.size(); ++i) {
                const int w = hifiasm_ha_ov_type(ol[i], rl);
                bool keep = (ol[i].shared_seed >= s[w]);
                if (!keep && w == 3 && cwn > 0) {
                    keep = should_rescue_type3(ol[i]);
                }
                if (keep) {
                    update_cc(ol[i]);
                    if (chain_cutoff >= 2 && ol[i].align_length < chain_cutoff) lch = 1;
                    if (k != i) std::swap(ol[k], ol[i]);
                    ++k;
                }
            }
            ol.resize(k);
        }
    }

    // ks_introsort_or_xs: sort by x_pos_s (and x_pos_e for determinism).
    std::sort(ol.begin(), ol.end(),
        [](const HifiasmOverlapRegion& a, const HifiasmOverlapRegion& b) {
            if (a.x_pos_s != b.x_pos_s) return a.x_pos_s < b.x_pos_s;
            if (a.x_pos_e != b.x_pos_e) return a.x_pos_e < b.x_pos_e;
            if (a.y_id != b.y_id) return a.y_id < b.y_id;
            if (a.y_pos_strand != b.y_pos_strand) return a.y_pos_strand < b.y_pos_strand;
            if (a.y_pos_s != b.y_pos_s) return a.y_pos_s < b.y_pos_s;
            return a.y_pos_e < b.y_pos_e;
        });

    if (lch) {
        // r485 weak-overlap suppression block (anchor.cpp:2061-2096), adapted to use
        // each overlap's explicit chain-hit list in `chainHitIndexFlat`.
        size_t l = 0;
        for (size_t i = 0; i < ol.size(); ++i) {
            if (ol[i].align_length < chain_cutoff) {
                const uint64_t zs = uint64_t(ol[i].x_pos_s);
                const uint64_t ze = uint64_t(ol[i].x_pos_e) + 1ULL;
                uint64_t ob = uint64_t(double(ze - zs) * HIFIASM_OFL);
                if (ob < 16) ob = 16;
                const int64_t osc = int64_t(ol[i].shared_seed) * int64_t(HIFIASM_CH_SC);
                const uint64_t ocn = uint64_t(ol[i].align_length) << HIFIASM_CH_OCC;

                size_t k = 0;
                for (; k < ol.size() && ze > uint64_t(ol[k].x_pos_s); ++k) {
                    if (ol[k].align_length < chain_cutoff) continue;
                    if (uint64_t(ol[k].align_length) < ocn) continue;
                    if (int64_t(ol[k].shared_seed) < osc) continue;

                    const uint64_t rs = uint64_t(ol[k].x_pos_s);
                    const uint64_t re = uint64_t(ol[k].x_pos_e) + 1ULL;
                    const uint64_t os = (rs >= zs) ? rs : zs;
                    const uint64_t oe = (re <= ze) ? re : ze;
                    if (!((oe > os) && (oe - os) >= ob)) continue;

                    uint64_t kn = 0;
                    const uint64_t off = uint64_t(ol[k].non_homopolymer_errors);
                    const uint64_t n_hit = uint64_t(ol[k].align_length);
                    for (uint64_t j = 0; j < n_hit && kn < ocn; ++j) {
                        const uint64_t g = uint64_t(chainHitIndexFlat[size_t(off + j)]);
                        if (g >= allHits.size()) continue;
                        const uint64_t me = uint64_t(allHits[size_t(g)].self_offset);
                        const uint64_t span = uint64_t(allHits[size_t(g)].cnt & 0xffu);
                        // Match hifiasm's unsigned wrap semantics: ms = me - span.
                        const uint64_t ms = me - span;
                        if (ms >= os && me <= oe) ++kn;
                    }
                    if (kn >= ocn) break;
                }
                if (k < ol.size() && ze > uint64_t(ol[k].x_pos_s)) {
                    // Suppress this weak overlap.
                    continue;
                }
            }

            if (l != i) std::swap(ol[l], ol[i]);
            ++l;
        }
        ol.resize(l);
    }
}

} // anonymous namespace

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
static inline uint32_t computeInvertedIndexHitWeight(
    const uint32_t count,
    const uint64_t lowFreqThreshold,
    const uint64_t highFreqThreshold,
    const uint64_t highFreqWeightUnit,
    const uint32_t rareKmerWeight,
    const vector<uint32_t>& weightLut,
    const double weightExponent)
{
    if (count <= lowFreqThreshold) {
        return rareKmerWeight;
    }
    if (count >= highFreqThreshold) {
        const uint32_t w = 1U + uint32_t((uint64_t(count) + highFreqWeightUnit - 1ULL) / highFreqWeightUnit);
        if (w < weightLut.size()) {
            return weightLut[w];
        }
        return uint32_t(std::pow(double(w), weightExponent));
    }
    return uint32_t(1);
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
        data.weightLut[i] = uint32_t(std::pow(double(i), data.weightExponent));
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
    data.strand0CanonicalIsRc.clear();
    data.strand0CanonicalIsRc.shrink_to_fit();
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
        maxChainLimit(maxChainLimit)
    {
        // Kept for API compatibility. In strict hifiasm ONT parity mode, `chain_cutoff` is fixed
        // (2) and we do not apply an additional per-read-pair minimum here.
        (void)minChainedMarkerCount;

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
            ReadId partnerReadId = invalidReadId;
            uint32_t querySpan = 0; // pre-extension span on the query read (bases)
        };
	        vector<EmittedChainedCandidate> emittedForRead;
	        vector<EmittedChainedCandidate> filteredForRead;
		        struct PendingHighFrequencyMarker {
		            uint64_t startIdx = 0;
		            uint32_t count = 0;
		            uint64_t hashKey = 0; // tie-break key for high-frequency streak selection (hifiasm-like)
		            uint32_t posA = 0;
		            uint32_t ordinalA = 0;
		            uint32_t weight = 1;
		            uint8_t isRcA = 0;
		        };
	        vector<PendingHighFrequencyMarker> highFrequencyStreak;
	        highFrequencyStreak.reserve(64);
	        vector<size_t> highFrequencyStreakWorkspace;
	        highFrequencyStreakWorkspace.reserve(64);

        // Scratch storage for the strict hifiasm lchain_qdp_mcopy_fast port.
        HifiasmChainDataScratch hifiasmChainDpScratch;
        vector<HifiasmKmerHit> hifiasmAllHits;
        vector<HifiasmKmerHit> hifiasmPairHits;
        vector<HifiasmOverlapRegion> hifiasmOverlapRegions;
        vector<uint32_t> hifiasmChainHitIndexFlat;
        
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
	        const uint64_t highFreqThreshold = std::max<uint64_t>(1ULL, uint64_t(double(coveragePeak) * highFreqMultiplier));
	        // Hifiasm's high-occ minimizer suppression is driven by its filter table + minimizer selection,
	        // not by this (hom_cov-derived) high-occ threshold directly. In very-low-coverage/synthetic
	        // situations `highFreqThreshold` can be 1, which would incorrectly treat almost every kmer
	        // as "high-frequency" and cause aggressive downsampling. We gate downsampling behind a small
	        // absolute minimum so only genuinely repetitive markers are downsampled.
	        const uint64_t highFreqDownsampleThreshold = std::max<uint64_t>(3ULL, highFreqThreshold);
	        const uint64_t highFreqWeightUnit = std::max<uint64_t>(1ULL, highFreqThreshold * 2ULL);
        const bool downsampleHighFrequencyMarkers =
            invertedIndexData.downsampleHighFrequencyMarkers &&
            invertedIndexData.highFrequencySampleDistance > 0 &&
            invertedIndexData.maxHighFrequencyPerStreak > 0;
        const uint32_t highFrequencySampleDistance = std::max<uint32_t>(1U, invertedIndexData.highFrequencySampleDistance);
        const uint32_t maxHighFrequencyPerStreak = std::max<uint32_t>(1U, invertedIndexData.maxHighFrequencyPerStreak);
        const bool lchainIsAccurate = invertedIndexData.lchainIsAccurate;
        const bool enableMcopyFast = invertedIndexData.enableMcopyFast;
        const uint32_t mcopyNum = std::max<uint32_t>(1U, invertedIndexData.mcopyNum);
        const double mcopyRate = std::max<double>(0.0, std::min<double>(1.0, invertedIndexData.mcopyRate));
	        const uint32_t mcopyKhitCutoff = std::max<uint32_t>(1U, invertedIndexData.mcopyKhitCutoff);
	        const uint32_t mcopyOcvWindow = std::max<uint32_t>(1U, invertedIndexData.mcopyOcvWindow);
	        // Dinara only supports the strict hifiasm ONT lchain+mcopy path for inverted-index chaining.


        uint64_t startBatch, endBatch;
        while(getNextBatch(startBatch, endBatch)) {
            for(ReadId readIdA = ReadId(startBatch); readIdA != ReadId(endBatch); ++readIdA) {
                
                const OrientedReadId orientedReadIdA(readIdA, 0);
                const auto& markersA = markers[orientedReadIdA.getValue()];
	                const bool haveCanonicalCache =
	                    (size_t(readIdA) + 1 < invertedIndexData.strand0CanonicalOffsets.size());
	                const KmerId* canonicalIdsA = nullptr;
	                const uint8_t* canonicalIsRcA = nullptr;
	                size_t canonicalCountA = 0;
	                if(haveCanonicalCache) {
	                    const uint64_t b = invertedIndexData.strand0CanonicalOffsets[size_t(readIdA)];
	                    const uint64_t e = invertedIndexData.strand0CanonicalOffsets[size_t(readIdA) + 1];
	                    if(e >= b &&
	                       e <= invertedIndexData.strand0CanonicalKmerIds.size() &&
	                       e <= invertedIndexData.strand0CanonicalIsRc.size()) {
	                        canonicalIdsA = invertedIndexData.strand0CanonicalKmerIds.data() + b;
	                        canonicalIsRcA = invertedIndexData.strand0CanonicalIsRc.data() + b;
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

	                auto computeHitWeight = [&](const uint32_t count) -> uint32_t {
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
	                            const uint32_t posBEncoded = compactOccs[j].position;
	                            const uint32_t posB = posBEncoded & 0x7fffffffU;
	                            const uint8_t isRcB = uint8_t(posBEncoded >> 31);
	                            scratch.flatHits.push_back(
	                                {compactOccs[j].readId, markerInfo.posA, posB, markerInfo.ordinalA, markerInfo.weight, markerInfo.isRcA, isRcB});
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
                // 3. Select 'keep' markers from the streak using hifiasm's ordering:
                //    - Prefer smallest genome-wide occurrence count (more informative)
                //    - Tie-break by a deterministic hash key
                //    Then emit the selected markers in increasing posA order.
                //
                // Hifiasm detail (sketch.cpp hf_select):
                // It only selects a marker if `count < (pe-ps)` (where pe-ps is the streak span),
                // which effectively discards extremely high-occ kmers that are "too repetitive" even
                // relative to the local span. We replicate that guard as `count < span`.
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
                // Hifiasm sketch.cpp: for a streak of high-occ minimizers, choose up to
                // `max_high_occ = round((pe-ps)/sample_dist)` (capped) using the smallest (count, hash) keys.
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
                        // Hifiasm sketch.cpp: if `max_high_occ = round((pe-ps)/sample_dist)` is 0,
                        // it treats this as "no high-frequency streak worth downsampling".
                        // In that case, do not drop evidence; keep all markers in this streak.
                        for(const auto& m : highFrequencyStreak) {
                            appendMarkerHits(m);
                        }
                        highFrequencyStreak.clear();
                        return;
                    }
                        auto appendIfUseful = [&](const PendingHighFrequencyMarker& m) {
                            // Hifiasm guard: keep only if `count < span`.
                            if (span > 0 && m.count >= span) {
                                return;
                            }
                            appendMarkerHits(m);
                        };

                        const size_t selectedCount = std::min<size_t>(highFrequencyStreak.size(), keep);
                        if(selectedCount >= highFrequencyStreak.size()) {
                            for(const auto& m : highFrequencyStreak) {
                                appendIfUseful(m);
                            }
                            highFrequencyStreak.clear();
                            return;
                        }

                    highFrequencyStreakWorkspace.clear();
                    highFrequencyStreakWorkspace.reserve(highFrequencyStreak.size());
                    for(size_t i = 0; i < highFrequencyStreak.size(); ++i) {
                        highFrequencyStreakWorkspace.push_back(i);
                    }

                    auto better = [&](const size_t a, const size_t b) {
                        const auto& ma = highFrequencyStreak[a];
                        const auto& mb = highFrequencyStreak[b];
                        if(ma.count != mb.count) {
                            return ma.count < mb.count;
                        }
                        return ma.hashKey < mb.hashKey;
                    };
                    std::nth_element(
                        highFrequencyStreakWorkspace.begin(),
                        highFrequencyStreakWorkspace.begin() + selectedCount,
                        highFrequencyStreakWorkspace.end(),
                        better);
                    highFrequencyStreakWorkspace.resize(selectedCount);
                        std::sort(
                            highFrequencyStreakWorkspace.begin(),
                            highFrequencyStreakWorkspace.end(),
                            [&](const size_t a, const size_t b) {
                                return highFrequencyStreak[a].posA < highFrequencyStreak[b].posA;
                            });
                        for(const size_t idx : highFrequencyStreakWorkspace) {
                            appendIfUseful(highFrequencyStreak[idx]);
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
	                    uint8_t isRcA = 0;
	                    if(canonicalIdsA) {
	                        canonicalKId = canonicalIdsA[i];
	                        isRcA = canonicalIsRcA ? canonicalIsRcA[i] : uint8_t(0);
	                    } else {
	                        KmerId currentKId = kmerIdsA[i];
	                        KmerId rcKId = getRcKmerId(currentKId, kmerLen);
	                        canonicalKId = (currentKId < rcKId) ? currentKId : rcKId;
	                        isRcA = uint8_t(currentKId > rcKId);
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

		                    const uint32_t hitWeight = computeHitWeight(count);
		                    const uint64_t kmerHashKey = hifiasmYakHash64_64(foldKmerIdToUint64(canonicalKId));
		                    if(downsampleHighFrequencyMarkers && count > highFreqDownsampleThreshold) {
		                        highFrequencyStreak.push_back({startIdx, count, kmerHashKey, posA, uint32_t(i), hitWeight, isRcA});
		                        continue;
		                    }

                    if(downsampleHighFrequencyMarkers) {
                        flushHighFrequencyStreak(posA);
                    }
		                    appendMarkerHits({startIdx, count, kmerHashKey, posA, uint32_t(i), hitWeight, isRcA});
	                    lastNonHighBoundaryPos = posA;
	                }
                if(downsampleHighFrequencyMarkers && !highFrequencyStreak.empty()) {
                    const uint32_t readLenABoundary = uint32_t(std::min<uint64_t>(
                        readLenA, uint64_t(std::numeric_limits<uint32_t>::max())));
                    flushHighFrequencyStreak(readLenABoundary);
                }
                
	                if(scratch.flatHits.empty()) continue;
	                radixSortFlatHitsByPartnerReadIdAndPosA(scratch.flatHits, scratch.flatHitsTmp);

	                {
	                    // ================================================================
	                    // STEP 2 (HIFIASM ONT PARITY): lchain_qdp_mcopy_fast + mcopy + ocv_w
	                    // ================================================================
	                    // Reference:
	                    // - anchor.cpp:1920-2100 lchain_qgen_mcopy_fast
	                    // - Hash_Table.cpp:2007-2095 quick_ck_lchain
	                    // - Hash_Table.cpp:2097-2284 lchain_qdp_mcopy_fast
                    //
                    // Strict parity choices:
                    // - Build a hifiasm-like `k_mer_hit` array per partner read B that includes both strands.
                    // - Sort by (strand, self_offset, offset) like minimizers_qgen0.
                    // - Run the strict lchain_qdp_mcopy_fast port to emit best chain + optional mcopy alternates.
                    // - Apply strict max_n_chain + ocv_w rescue + r485 suppression on the aggregated overlap list.
                    // - Then convert surviving overlaps into Dinara `Alignment` objects.

                    const HifiasmLchainDpOptions dpOpt =
                        getHifiasmLchainDpOptions(lchainIsAccurate, uint32_t(kmerLen));
                    const uint8_t span = uint8_t(std::min<uint64_t>(kmerLen, 255ULL));

                    const uint32_t readLenA32 = uint32_t(std::min<uint64_t>(
                        readLenA, uint64_t(std::numeric_limits<uint32_t>::max())));

                    hifiasmChainDpScratch.resize(0);
                    hifiasmAllHits.clear();
                    hifiasmPairHits.clear();
                    hifiasmOverlapRegions.clear();
                    hifiasmChainHitIndexFlat.clear();

                    // Iterate per partner readIdB (flatHits are sorted by (readIdB, posA)).
                    size_t hitIter = 0;
                    while(hitIter < scratch.flatHits.size()) {
                        const ReadId readIdB = scratch.flatHits[hitIter].partnerReadId;

                        // Skip mirrored pairs and self-comparisons.
                        if(readIdB <= readIdA) {
                            while(hitIter < scratch.flatHits.size() && scratch.flatHits[hitIter].partnerReadId == readIdB) {
                                ++hitIter;
                            }
                            continue;
                        }

                        const size_t startInFlat = hitIter;
                        while(hitIter < scratch.flatHits.size() && scratch.flatHits[hitIter].partnerReadId == readIdB) {
                            ++hitIter;
                        }
                        const size_t numHits = hitIter - startInFlat;
                        if(numHits == 0) {
                            continue;
                        }

                        // Hifiasm ONT EC parity: `chain_cutoff` is passed as a constant 2
                        // (ecovlp.cpp:3274) and is applied inside minimizers_qgen_input(...)
                        // (anchor.cpp:1489) as a per-(readIdB, rev) minimum bucket size.
                        // In other words, a pair can be kept if (rev==0 has >=2 hits) OR (rev==1 has >=2 hits),
                        // and the smaller bucket is dropped entirely.
                        static constexpr uint32_t hifiasmChainCutoff = 2;
                        size_t revCount[2] = {0, 0};
                        for(size_t k = 0; k < numHits; ++k) {
                            const auto& h = scratch.flatHits[startInFlat + k];
                            const uint8_t rev = uint8_t(h.isRcA ^ h.isRcB); // z->rev ^ y->rev
                            ++revCount[rev ? 1 : 0];
                        }
                        const bool keepRev0 = (revCount[0] >= size_t(hifiasmChainCutoff));
                        const bool keepRev1 = (revCount[1] >= size_t(hifiasmChainCutoff));
                        if(!keepRev0 && !keepRev1) {
                            continue;
                        }

                        const uint64_t readLenB = reads.getReadRawSequenceLength(readIdB);
                        const OrientedReadId orientedReadB(readIdB, 0);
                        const auto& mB = markers[orientedReadB.getValue()];
                        // Map posB -> ordinalB for this (readA, readB) group.
                        scratch.hitPosB.resize(numHits);
                        scratch.hitOrdinalB.assign(numHits, std::numeric_limits<uint32_t>::max());
                        for(size_t k = 0; k < numHits; ++k) {
                            scratch.hitPosB[k] = scratch.flatHits[startInFlat + k].posB;
                        }
                        if(!mapHitPositionsToMarkerOrdinals(
                            scratch.hitPosB, mB, scratch.hitOrdinalB, scratch.hitOrderByPosB)) {
                            continue;
                        }

                        hifiasmPairHits.clear();
                        hifiasmPairHits.reserve(numHits);

                        for(size_t k = 0; k < numHits; ++k) {
                            const auto& h = scratch.flatHits[startInFlat + k];
                            const uint8_t rev = uint8_t(h.isRcA ^ h.isRcB); // z->rev ^ y->rev
                            if((rev == 0 && !keepRev0) || (rev == 1 && !keepRev1)) {
                                continue;
                            }

                            // Hifiasm uses minimizer end positions; Dinara stores marker start positions.
                            // Convert start -> end by adding span-1.
                            const uint32_t seedSpan = uint32_t(span);
                            const uint32_t selfOff = h.posA + (seedSpan - 1U);
                            const uint32_t offSame = h.posB + (seedSpan - 1U);
                            // anchor.cpp:1062-1064 simplifies to offset = tl - start - 1 for strand==1.
                            const uint32_t offDiff = uint32_t(readLenB - 1ULL - uint64_t(h.posB));

                            HifiasmKmerHit kh{};
                            kh.readID = uint32_t(readIdB);
                            kh.strand = rev;
                            kh.self_offset = selfOff;
                            kh.offset = (rev == 0) ? offSame : offDiff;

                            const uint32_t w = (h.weight > 0xffffffu) ? 0xffffffu : h.weight;
                            const uint32_t span8 = std::min<uint32_t>(seedSpan, 0xffu);
                            kh.cnt = (w << 8) | span8;

                            kh.ordinalA = h.ordinalA;
                            kh.ordinalB = scratch.hitOrdinalB[k];
                            kh.globalIndex = uint32_t(hifiasmAllHits.size());

                            hifiasmAllHits.push_back(kh);
                            hifiasmPairHits.push_back(kh);
                        }

                        // Hifiasm requires strand-groups sorted by self_offset and offset.
                        std::sort(hifiasmPairHits.begin(), hifiasmPairHits.end(),
                            [](const HifiasmKmerHit& a, const HifiasmKmerHit& b) {
                                if(a.strand != b.strand) return a.strand < b.strand;
                                if(a.self_offset != b.self_offset) return a.self_offset < b.self_offset;
                                return a.offset < b.offset;
                            });

                        const int64_t mcopy_num = enableMcopyFast ? int64_t(mcopyNum) : 1;
                        const double mcopy_rate = enableMcopyFast ? mcopyRate : 0.0;

                        hifiasm_lchain_qdp_mcopy_fast(
                            hifiasmPairHits,
                            hifiasmChainDpScratch,
                            hifiasmOverlapRegions,
                            hifiasmChainHitIndexFlat,
                            dpOpt.maxSkip,
                            dpOpt.maxIter,
                            dpOpt.maxDist,
                            dpOpt.chnPenGap,
                            dpOpt.chnPenSkip,
                            maxDriftRate,
                            uint32_t(readIdA),
                            int64_t(readLenA),
                            int64_t(readLenB),
                            dpOpt.quickCheck ? 1 : 0,
                            mcopy_num,
                            mcopy_rate,
                            int64_t(mcopyKhitCutoff));
                    }

                    // Dinara-only optional pruning before hifiasm max_n_chain/ocv_w logic.
                    if((invertedIndexData.minOverlapLength > 0 || invertedIndexData.maxEndFuzz > 0) &&
                        !hifiasmOverlapRegions.empty()) {
                        vector<HifiasmOverlapRegion> tmp;
                        tmp.reserve(hifiasmOverlapRegions.size());
                        for(const auto& r : hifiasmOverlapRegions) {
                            const uint64_t qSpan = uint64_t(r.raw_xe) - uint64_t(r.raw_xs);

                            uint32_t tS = r.raw_ys;
                            uint32_t tE = r.raw_ye;
                            const uint64_t readLenB = reads.getReadRawSequenceLength(ReadId(r.y_id));
                            if(r.y_pos_strand) {
                                const auto p = dinara::rcIntervalToForward(uint32_t(readLenB), tS, tE);
                                tS = p.first;
                                tE = p.second;
                            }
                            const uint64_t tSpan = uint64_t(tE) - uint64_t(tS);

                            if(invertedIndexData.minOverlapLength > 0) {
                                if(std::min(qSpan, tSpan) < uint64_t(invertedIndexData.minOverlapLength)) {
                                    continue;
                                }
                            }

                            if(invertedIndexData.maxEndFuzz > 0) {
                                const uint32_t qPstart = r.raw_xs;
                                const uint32_t qPend = r.raw_xe;
                                const uint32_t leftNeed = std::min(qPstart, tS);
                                const int64_t qRight = int64_t(readLenA) - int64_t(qPend);
                                const int64_t tRight = int64_t(readLenB) - int64_t(tE);
                                const uint32_t rightNeed = uint32_t(std::min<int64_t>(
                                    std::max<int64_t>(qRight, 0), std::max<int64_t>(tRight, 0)));
                                if(leftNeed > invertedIndexData.maxEndFuzz || rightNeed > invertedIndexData.maxEndFuzz) {
                                    continue;
                                }
                            }

                            tmp.push_back(r);
                        }
                        hifiasmOverlapRegions.swap(tmp);
                    }

                    // Strict hifiasm max_n_chain + ocv_w rescue + r485 suppression.
                    // Note: In hifiasm ONT EC (`ecovlp.cpp:3274`), `chain_cutoff` is passed as a constant 2.
                    // This cutoff controls when the r485 weak-overlap suppression block is enabled
                    // (it flags that at least one overlap has `align_length < chain_cutoff`).
                    // It is NOT the same as Dinara's `minChainedMarkerCount`.
                    static constexpr uint32_t hifiasmChainCutoffForR485 = 2;
                    hifiasm_lchain_qgen_mcopy_fast_postfilter(
                        hifiasmOverlapRegions,
                        maxChainLimit,
                        hifiasmChainCutoffForR485,
                        uint64_t(mcopyOcvWindow),
                        readLenA32,
                        hifiasmChainHitIndexFlat,
                        hifiasmAllHits);

                    // Convert surviving overlaps to Dinara candidates/alignments.
                    for(const auto& r : hifiasmOverlapRegions) {
                        const ReadId readIdB = ReadId(r.y_id);
                        const uint64_t readLenB = reads.getReadRawSequenceLength(readIdB);
                        const OrientedReadId orientedReadB(readIdB, 0);
                        const auto& mB = markers[orientedReadB.getValue()];
                        const uint32_t markerCountA = uint32_t(markersA.size());
                        const uint32_t markerCountB = uint32_t(mB.size());

                        const uint32_t qS = r.x_pos_s;
                        const uint32_t qE = r.x_pos_e + 1U;

                        uint32_t tS = r.y_pos_s;
                        uint32_t tE = r.y_pos_e + 1U;
                        bool isSameStrand = (r.y_pos_strand == 0);
                        if(!isSameStrand) {
                            const auto p = dinara::rcIntervalToForward(uint32_t(readLenB), tS, tE);
                            tS = p.first;
                            tE = p.second;
                        }

                        Alignment al;
                        const uint64_t off = uint64_t(r.non_homopolymer_errors);
                        const uint64_t n_hit = uint64_t(r.align_length);
                        al.ordinals.reserve(n_hit);
                        for(uint64_t j = 0; j < n_hit; ++j) {
                            const uint64_t g = uint64_t(hifiasmChainHitIndexFlat[size_t(off + j)]);
                            if(g >= hifiasmAllHits.size()) {
                                continue;
                            }
                            const auto& h = hifiasmAllHits[size_t(g)];
                            uint32_t ordB = h.ordinalB;
                            if(!isSameStrand) {
                                ordB = markerCountB - 1U - ordB;
                            }
                            al.ordinals.push_back({h.ordinalA, ordB});
                        }
                        if(al.ordinals.empty()) {
                            continue;
                        }

                        al.qs = qS;
                        al.qe = qE;
                        al.ts = tS;
                        al.te = tE;

                        // Canonicalize candidate so readIds[0] < readIds[1], and keep alignment consistent.
                        ReadId cand0 = readIdA;
                        ReadId cand1 = readIdB;
                        canonicalizeCandidateAndAlignment(cand0, cand1, isSameStrand, al, markerCountA, markerCountB);

                        const uint8_t overlapType =
                            uint8_t(getOverlapType(qS, qE, uint32_t(readLenA)));
                        const uint64_t preSpan64 = uint64_t(r.raw_xe) - uint64_t(r.raw_xs);
                        const uint32_t preSpan = (preSpan64 > uint64_t(std::numeric_limits<uint32_t>::max())) ?
                            std::numeric_limits<uint32_t>::max() : uint32_t(preSpan64);

                        emittedForRead.push_back(EmittedChainedCandidate{
                            OrientedReadPair(cand0, cand1, isSameStrand),
                            std::move(al),
                            r.shared_seed,
                            overlapType,
                            (cand0 == readIdA) ? cand1 : cand0,
                            preSpan});
	                    }

		                } // end lchain path

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
	                // Note: max_n_chain + ocv_w rescue + r485 suppression were already applied in the lchain stage
	                // (hifiasm_lchain_qgen_mcopy_fast_postfilter), so do not filter again here.

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
    invertedIndexData.strand0CanonicalIsRc.resize(totalMarkersFound);
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
                const uint8_t isRc = uint8_t(kId > rcKId);

                // Encode the observed-orientation bit into the top bit of the stored position.
                // Marker positions fit in 24 bits (Uint24), so this does not collide with real coordinates.
                const uint32_t encodedPosition = uint32_t(rMarkers[i].position) | (uint32_t(isRc) << 31);
                invertedIndexData.occurrences[writeOffset++] = {
                    canonicalKId,
                    rId,
                    encodedPosition
                };
                invertedIndexData.strand0CanonicalKmerIds[canonicalWriteOffset++] = canonicalKId;
                invertedIndexData.strand0CanonicalIsRc[canonicalWriteOffset - 1] = isRc;
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
	        1ULL, uint64_t(double(coveragePeak) * invertedIndexData.highFreqMultiplier));
	    const uint64_t highFreqDownsampleThreshold = std::max<uint64_t>(3ULL, highFreqThreshold);
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
	                uint64_t hashKey = 0; // hifiasm sketch.cpp: tie-break key (count, hash)
	                uint32_t posA = 0;
	                uint32_t ordinalA = 0;
	                uint32_t weight = 1;
	                uint8_t isRcA = 0;
	            };
	            vector<PendingHighFrequencyMarker> highFrequencyStreak;
	            highFrequencyStreak.reserve(64);
	            vector<size_t> highFrequencyStreakWorkspace;
	            highFrequencyStreakWorkspace.reserve(64);

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
	                    const uint8_t* canonicalIsRcA = nullptr;
	                    size_t canonicalCountA = 0;
	                    if(haveCanonicalCache) {
	                        const uint64_t b = invertedIndexData.strand0CanonicalOffsets[size_t(readIdA)];
	                        const uint64_t e = invertedIndexData.strand0CanonicalOffsets[size_t(readIdA) + 1];
	                        if(e >= b &&
	                           e <= invertedIndexData.strand0CanonicalKmerIds.size() &&
	                           e <= invertedIndexData.strand0CanonicalIsRc.size()) {
	                            canonicalIdsA = invertedIndexData.strand0CanonicalKmerIds.data() + b;
	                            canonicalIsRcA = invertedIndexData.strand0CanonicalIsRc.data() + b;
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

	                    auto computeHitWeight = [&](const uint32_t count) -> uint32_t {
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
		                                const uint32_t posBEncoded = compactOccs[j].position;
		                                const uint32_t posB = posBEncoded & 0x7fffffffU;
		                                const uint8_t isRcB = uint8_t(posBEncoded >> 31);
		                                scratch.flatHits.push_back(
		                                    {readIdB, markerInfo.posA, posB, markerInfo.ordinalA, markerInfo.weight, markerInfo.isRcA, isRcB});
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
	                            for(const auto& m : highFrequencyStreak) {
	                                appendMarkerHits(m);
	                            }
	                            highFrequencyStreak.clear();
	                            return;
	                        }
	                        auto appendIfUseful = [&](const PendingHighFrequencyMarker& m) {
	                            if (span > 0 && m.count >= span) {
	                                return;
	                            }
	                            appendMarkerHits(m);
	                        };

	                        const size_t selectedCount = std::min<size_t>(highFrequencyStreak.size(), keep);
	                        if(selectedCount >= highFrequencyStreak.size()) {
	                            for(const auto& m : highFrequencyStreak) {
	                                appendIfUseful(m);
	                            }
	                            highFrequencyStreak.clear();
	                            return;
	                        }

	                        highFrequencyStreakWorkspace.clear();
	                        highFrequencyStreakWorkspace.reserve(highFrequencyStreak.size());
	                        for(size_t i = 0; i < highFrequencyStreak.size(); ++i) {
	                            highFrequencyStreakWorkspace.push_back(i);
	                        }

	                        auto better = [&](const size_t a, const size_t b) {
	                            const auto& ma = highFrequencyStreak[a];
	                            const auto& mb = highFrequencyStreak[b];
	                            if(ma.count != mb.count) {
	                                return ma.count < mb.count;
	                            }
	                            return ma.hashKey < mb.hashKey;
	                        };
	                        std::nth_element(
	                            highFrequencyStreakWorkspace.begin(),
	                            highFrequencyStreakWorkspace.begin() + selectedCount,
	                            highFrequencyStreakWorkspace.end(),
	                            better);
	                        highFrequencyStreakWorkspace.resize(selectedCount);
	                        std::sort(
	                            highFrequencyStreakWorkspace.begin(),
	                            highFrequencyStreakWorkspace.end(),
	                            [&](const size_t a, const size_t b) {
	                                return highFrequencyStreak[a].posA < highFrequencyStreak[b].posA;
	                            });
	                        for(const size_t idx : highFrequencyStreakWorkspace) {
	                            appendIfUseful(highFrequencyStreak[idx]);
	                        }
	                        highFrequencyStreak.clear();
	                    };

	                    for(size_t i = 0; i < numMarkersA; i++) {
	                        KmerId canonicalKId;
	                        uint8_t isRcA = 0;
	                        if(canonicalIdsA) {
	                            canonicalKId = canonicalIdsA[i];
	                            isRcA = canonicalIsRcA ? canonicalIsRcA[i] : uint8_t(0);
	                        } else {
	                            KmerId currentKId = kmerIdsA[i];
	                            KmerId rcKId = getRcKmerId(currentKId, kmerLen);
	                            canonicalKId = (currentKId < rcKId) ? currentKId : rcKId;
	                            isRcA = uint8_t(currentKId > rcKId);
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

	                        const uint32_t hitWeight = computeHitWeight(count);
	                        const uint64_t kmerHashKey = hifiasmYakHash64_64(foldKmerIdToUint64(canonicalKId));
		                        if(downsampleHighFrequencyMarkers && count > highFreqDownsampleThreshold) {
		                            highFrequencyStreak.push_back({startIdx, count, kmerHashKey, posA, uint32_t(i), hitWeight, isRcA});
		                            continue;
		                        }

	                        if(downsampleHighFrequencyMarkers) {
	                            flushHighFrequencyStreak(posA);
	                        }
	                        appendMarkerHits({startIdx, count, kmerHashKey, posA, uint32_t(i), hitWeight, isRcA});
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
                            const int32_t wI =
                                (scratch.hitWeights[i] > uint32_t(std::numeric_limits<int32_t>::max())) ?
                                std::numeric_limits<int32_t>::max() : int32_t(scratch.hitWeights[i]);
                            sc = HIFIASM_NORMAL_W(sc, wI);

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
                                    const int32_t wI =
                                        (scratch.hitWeights[i] > uint32_t(std::numeric_limits<int32_t>::max())) ?
                                        std::numeric_limits<int32_t>::max() : int32_t(scratch.hitWeights[i]);
                                    tmp = HIFIASM_NORMAL_W(tmp, wI);

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

    // ============================================================
    // HIFIASM "ONE BEST OVERLAP PER PARTNER READ (y_id)" (PAF PATH)
    // ============================================================
    // Mirror the same pre-EC behavior as discovery chaining: per query read,
    // keep only the best overlap to each partner read by (score desc, overlapLen asc).
    //
    // This is intentionally done BEFORE max_n_chain filtering so the read-level
    // cap is not wasted on duplicate chains to the same partner.
    if(mergedResults.size() > 1) {
        std::sort(mergedResults.begin(), mergedResults.end(),
            [&](const PafChainedCandidate& a, const PafChainedCandidate& b) {
                if(a.queryReadId != b.queryReadId) {
                    return a.queryReadId < b.queryReadId;
                }
                const ReadId partnerA =
                    (a.candidate.readIds[0] == a.queryReadId) ? a.candidate.readIds[1] : a.candidate.readIds[0];
                const ReadId partnerB =
                    (b.candidate.readIds[0] == b.queryReadId) ? b.candidate.readIds[1] : b.candidate.readIds[0];
                if(partnerA != partnerB) {
                    return partnerA < partnerB;
                }
                if(a.score != b.score) {
                    return a.score > b.score;
                }
                const uint64_t spanA =
                    (a.candidate.readIds[0] == a.queryReadId) ? (uint64_t(a.alignment.qe) - uint64_t(a.alignment.qs))
                                                             : (uint64_t(a.alignment.te) - uint64_t(a.alignment.ts));
                const uint64_t spanB =
                    (b.candidate.readIds[0] == b.queryReadId) ? (uint64_t(b.alignment.qe) - uint64_t(b.alignment.qs))
                                                             : (uint64_t(b.alignment.te) - uint64_t(b.alignment.ts));
                if(spanA != spanB) {
                    return spanA < spanB;
                }
                if(a.candidate.isSameStrand != b.candidate.isSameStrand) {
                    return a.candidate.isSameStrand < b.candidate.isSameStrand;
                }
                if(a.candidate.readIds[0] != b.candidate.readIds[0]) {
                    return a.candidate.readIds[0] < b.candidate.readIds[0];
                }
                return a.candidate.readIds[1] < b.candidate.readIds[1];
            });

        vector<PafChainedCandidate> deduped;
        deduped.reserve(mergedResults.size());
        size_t i = 0;
        while(i < mergedResults.size()) {
            size_t j = i + 1;
            const ReadId qid = mergedResults[i].queryReadId;
            const ReadId partnerI =
                (mergedResults[i].candidate.readIds[0] == qid) ? mergedResults[i].candidate.readIds[1] : mergedResults[i].candidate.readIds[0];
            while(j < mergedResults.size() && mergedResults[j].queryReadId == qid) {
                const ReadId partnerJ =
                    (mergedResults[j].candidate.readIds[0] == qid) ? mergedResults[j].candidate.readIds[1] : mergedResults[j].candidate.readIds[0];
                if(partnerJ != partnerI) {
                    break;
                }
                ++j;
            }
            // Because of the sort order, mergedResults[i] is the best for this (query, partner).
            deduped.push_back(std::move(mergedResults[i]));
            i = j;
        }
        mergedResults.swap(deduped);
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
