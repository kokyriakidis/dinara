/**
 * @file AssemblerInvertedIndex.cpp
 * @brief Overlap candidate discovery via inverted index + hifiasm-parity chaining.
 *
 * Pipeline: index construction -> hit collection -> frequency-based weighting ->
 * high-freq downsampling -> per-partner DP chaining -> mcopy extraction ->
 * postfilter (max_n_chain + COV_W rescue + R485 suppression).
 *
 * Scoring, tie-breaking, and filtering match hifiasm v0.25.0-r726 (ec9a8b2).
 *
 * Two entry points:
 *   - Discovery path: findAlignmentCandidatesInvertedIndex (all-vs-all)
 *   - PAF path: chainPafCandidates (re-chain imported pairs)
 */

#include "Assembler.hpp"
#include "InvertedIndexBuilder.hpp"
#include "hifiasmCoordinateTransforms.hpp"
#include "performanceLog.hpp"
#include "OrientedReadPair.hpp"
#include "ProjectedAlignment.hpp"
#include "timestamp.hpp"
#include "Reads.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <mutex>
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

// ============================================================================
// LSD RADIX SORT FOR HIT ARRAYS
// ============================================================================
// Sort hits by (partnerReadId, posA) so that all shared k-mers between the
// same read pair are contiguous and in monotonically increasing query order.
// This packed-key approach merges the two 32-bit fields into a single 64-bit
// value and radix-sorts on it byte-by-byte.
//
// Why radix sort over std::sort?
//   - O(N) complexity vs O(N log N) for comparison-based sorts
//   - Hifiasm uses radix sorts in the same hot path for the same reason
//   - With millions of hits per read, the constant factor matters
//
// Key packing:  key = (partnerReadId << 32) | posA
//   - High 32 bits: groups all hits for the same partner together
//   - Low 32 bits:  orders hits by query position within each group
// ============================================================================

/**
 * @brief Pack the two sort keys into a single 64-bit value for radix sort.
 * @param h  The hit to pack.
 * @return Packed key with partnerReadId in bits [63:32] and posA in bits [31:0].
 */
static inline uint64_t flatHitPackedKey(const InvertedIndexTempHit& h)
{
    return (uint64_t(h.partnerReadId) << 32) | uint64_t(h.posA);
}

/**
 * @brief LSD (Least Significant Digit) radix sort for hit arrays.
 *
 * Performs 8 passes of counting sort, one per byte of the 64-bit packed key,
 * starting from the least significant byte. Each pass is stable, so after all
 * 8 passes the hits are in correct lexicographic order by (partnerReadId, posA).
 *
 * @param hits  The array to sort (modified in place).
 * @param tmp   Auxiliary buffer of the same size (avoids repeated allocation).
 *
 * @complexity O(N) time, O(N) space for the auxiliary buffer.
 * @note Uses hits.swap(tmp) on each pass, so the final sorted result is
 *       always in `hits` after an even number of passes (8 passes = even).
 */
static void radixSortFlatHitsByPartnerReadIdAndPosA(
    vector<InvertedIndexTempHit>& hits,
    vector<InvertedIndexTempHit>& tmp)
{
    if(hits.size() <= 1) {
        return;
    }
    tmp.resize(hits.size());

    // LSD radix sort, 8 passes over the 64-bit packed key.
    // Pass i sorts by byte i (bits [8i+7 : 8i]), from least to most significant.
    // Stable counting-sort per byte ensures correct lexicographic order.
    for(int pass = 0; pass < 8; ++pass) {
        const uint32_t shift = uint32_t(pass * 8);

        // Step 1: Count occurrences of each byte value (histogram).
        std::array<size_t, 256> count{};
        for(const auto& h : hits) {
            const uint8_t b = uint8_t((flatHitPackedKey(h) >> shift) & 0xffULL);
            ++count[b];
        }

        // Step 2: Convert histogram to prefix sums (exclusive scan).
        // After this, count[b] = starting output index for items with byte value b.
        size_t sum = 0;
        for(size_t i = 0; i < count.size(); ++i) {
            const size_t c = count[i];
            count[i] = sum;
            sum += c;
        }

        // Step 3: Scatter items into their sorted positions in the tmp buffer.
        // Post-increment count[b] so the next item with the same byte goes one slot later.
        for(const auto& h : hits) {
            const uint8_t b = uint8_t((flatHitPackedKey(h) >> shift) & 0xffULL);
            tmp[count[b]++] = h;
        }

        // Step 4: Swap so the sorted result is in `hits` for next pass.
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

// ============================================================================
// HIFIASM DP CHAINING OPTIONS
// ============================================================================
/**
 * @brief Dynamic programming parameters for anchor chaining (hifiasm parity).
 *
 * These parameters control the DP search space and scoring penalties used
 * in the linear-time chaining algorithm. They are configured differently
 * for ONT vs HiFi reads to account for their error profiles.
 *
 * Reference: hifiasm anchor.cpp:2272 (set_lchain_dp_op)
 */
struct HifiasmLchainDpOptions {
    int32_t maxSkip = 25;       // Max anchors to skip before pruning (prevents poor chains)
    int32_t maxIter = 5000;     // Max lookback window for DP (limits O(N²) to O(N*maxIter))
    int32_t maxDist = 5000;     // Max coordinate distance for max_ii tracking
    double chnPenGap = 0.0;     // Gap penalty coefficient (indels/mismatches)
    double chnPenSkip = 0.0;    // Skip penalty (bases between consecutive anchors)
    bool quickCheck = false;    // Enable fast O(N) prefix validation (ONT=true, HiFi=false)
};

// ============================================================================
// HIFIASM CHAINING PARAMETER CONFIGURATION
// ============================================================================
/**
 * @brief Configure DP parameters for anchor chaining based on sequencing technology.
 *
 * This function implements hifiasm's set_lchain_dp_op logic (anchor.cpp:2272),
 * which adjusts chaining penalties based on:
 *   1. Sequencing technology (ONT vs HiFi)
 *   2. K-mer size (larger k-mers → more specific → trust gaps more)
 *
 * ## Penalty Formula:
 *   penalty = base * exp(-div * k)
 *
 * Where:
 *   - div: Exponential decay rate (0.01 for ONT, 0.1 for HiFi)
 *   - k: K-mer size
 *   - base: Penalty coefficient (penGap or penSkip)
 *
 * ## Intuition:
 *   Larger k-mers provide stronger evidence of true homology, so we penalize
 *   gaps less. Smaller k-mers are noisier, requiring stricter gap penalties.
 *
 * ## ONT vs HiFi:
 *   - **ONT** (isAccurate=true): div=0.01 → slower decay → higher penalties
 *     * More conservative chaining for noisy reads
 *     * Enables quickCheck optimization (O(N) fast path for collinear hits)
 *     * Example (k=34): chnPenGap ≈ 0.091, quickCheck=true
 *
 *   - **HiFi** (isAccurate=false): div=0.1 → faster decay → lower penalties
 *     * Less conservative for high-quality reads
 *     * Full DP without quick check
 *     * Example (k=34): chnPenGap ≈ 0.300, quickCheck=false
 *
 * @param isAccurate  True for ONT (more accurate chaining), false for HiFi
 * @param markerK     K-mer size used for minimizer extraction
 * @return Configured DP options
 *
 * @complexity O(1) - Single exp() call
 * @reference Hifiasm anchor.cpp:2272-2285
 */
static inline HifiasmLchainDpOptions getHifiasmLchainDpOptions(
    const bool isAccurate,
    const uint32_t markerK) noexcept
{
    HifiasmLchainDpOptions opt;

    // ONT uses quick check for O(N) fast path on collinear hits
    opt.quickCheck = isAccurate;

    // Exponential decay configuration (ONT is more conservative)
    constexpr double DIV_ONT = 0.01;   // Slower decay for noisy reads
    constexpr double DIV_HIFI = 0.1;   // Faster decay for accurate reads
    constexpr double PEN_GAP_BASE = 0.5;      // Gap penalty coefficient
    constexpr double PEN_SKIP_BASE = 0.0005;  // Skip penalty coefficient

    const double div = isAccurate ? DIV_ONT : DIV_HIFI;
    const double decayFactor = std::exp(-div * static_cast<double>(markerK));

    // Final penalties with exponential decay
    opt.chnPenGap = PEN_GAP_BASE * decayFactor;
    opt.chnPenSkip = PEN_SKIP_BASE * decayFactor;

    return opt;
}

/**
 * @brief Per-thread scratchpad for high-performance DP chaining.
 *
 * This struct uses a Structure-of-Arrays (SoA) layout for the DP arrays,
 * which is critical for cache efficiency during the O(N^2) chaining loop.
 * The scratchpad is recycled across reads to avoid allocation overhead.
 *
 * ## Memory Layout Strategy (AoS → SoA conversion):
 *
 *   Hits are first collected in AoS format (InvertedIndexTempHit) during
 *   the k-mer matching phase for convenient per-hit field access. Before
 *   the DP loop, they are unpacked into separate SoA vectors so the inner
 *   loop touches only the position arrays, keeping cache lines full of
 *   useful data.
 *
 *   AoS (collection):                 SoA (DP):
 *   ┌─────────────────────┐          hitPosA:    [p0, p1, p2, ...]
 *   │ hit0: readB posA .. │   ──►    hitPosB:    [q0, q1, q2, ...]
 *   │ hit1: readB posA .. │          hitWeights:  [w0, w1, w2, ...]
 *   │ hit2: readB posA .. │          dpSame:      [s0, s1, s2, ...]
 *   └─────────────────────┘          parentSame:  [-1, 0,  1, ...]
 *
 * ## Key Array Groups:
 *
 *   **Hit coordinates** (read from index, used in DP scoring):
 *   - `hitPosA`, `hitPosB`: Base positions of each anchor on reads A and B.
 *   - `hitOrdinalA`, `hitOrdinalB`: Marker ordinals for alignment output.
 *   - `hitWeights`: Frequency-based weights for score normalization.
 *
 *   **DP state** (populated during the chaining loop):
 *   - `dpSame`, `dpDiff`: Best DP score ending at each anchor (same/diff strand).
 *   - `parentSame`, `parentDiff`: Backtrack pointers for chain reconstruction.
 *   - `backtrackVisit*`: Visit marks for hifiasm's max_skip pruning heuristic.
 *   - `chainOccurrences*`: Anchor counts per chain, used for weak-chain filtering.
 *
 *   **Post-DP workspaces** (used during chain extraction and filtering):
 *   - `chainCandidates`: All chains sorted by (score desc, chainLen asc).
 *   - `mcopy*`: Workspace for multi-copy chain extraction.
 *   - `weakMetas`, `suppress*`: Workspace for R485 weak-chain suppression.
 *   - `strongAnchor*`: Flattened anchor positions for overlap-based filtering.
 */
struct ThreadScratchpad {
    /// @brief Metadata for a single chain during weak-chain suppression.
    /// Stores the query-coordinate span [qs, qe), anchor count, and score
    /// so we can quickly test overlap and dominance conditions.
    struct WeakFilterMeta {
        uint32_t qs = 0;
        uint32_t qe = 0;   ///< Half-open interval [qs, qe) on the query read.
        uint32_t occ = 0;  ///< Number of anchors in the chain.
        int32_t score = 0;  ///< DP score of the chain.
    };

    // ---- Stage 1: Hit Collection (AoS format for sorting) ----
    vector<InvertedIndexTempHit> flatHits;     ///< Collected k-mer hits before DP.
    vector<InvertedIndexTempHit> flatHitsTmp;   ///< Auxiliary buffer for radix sort.

    // ---- Stage 2: SoA Unpacking (cache-efficient DP layout) ----
    vector<uint32_t> hitPosA, hitPosB, hitOrdinalA, hitOrdinalB;
    vector<uint32_t> hitOrderByPosB;   ///< Permutation array: indices sorted by posB.
    vector<uint32_t> hitWeights;        ///< Frequency-based weight per anchor.

    // ---- Stage 3: DP State (dual-strand chaining) ----
    vector<int32_t> dpSame, dpDiff;                 ///< Best DP score ending at anchor i.
    vector<int32_t> parentSame, parentDiff;         ///< Predecessor for backtracking.
    vector<int32_t> backtrackVisitSame, backtrackVisitDiff; ///< Visit marks for max_skip.
    vector<uint32_t> chainOccurrencesSame, chainOccurrencesDiff; ///< Anchor count per chain.

    // ---- Stage 4: Chain Candidate Extraction ----
    /// @brief A candidate chain endpoint with score and metadata.
    /// Sorted by (score descending, chainLen ascending) for hifiasm-compatible
    /// tie-breaking: among equal-score chains, prefer shorter alignment length
    /// (i.e., tighter, more compact overlaps).
    struct ChainCandidate {
        int32_t score;      ///< DP chain score.
        uint64_t chainLen;  ///< Normalized alignment length for tie-breaking.
        int32_t endK;       ///< Index of the chain's terminal anchor in the SoA arrays.
        bool isDiff;        ///< True if this chain is on the reverse-complement strand.
        /// @brief Comparison: higher score first, then shorter chainLen.
        bool operator<(const ChainCandidate& other) const {
            if (score != other.score) return score > other.score;
            return chainLen < other.chainLen;
        }
    };
    vector<ChainCandidate> chainCandidates;       ///< All chain endpoints, sorted by score.
    vector<ChainCandidate> filteredCandidates;    ///< Surviving chains after suppression.

    // ---- Stage 5: Weak-Chain Suppression Workspace ----
    // Hifiasm's R485 algorithm suppresses low-anchor chains that overlap with
    // high-anchor chains on the query read. These arrays partition candidates
    // into weak (few anchors) and strong (many anchors) groups, then test
    // whether each weak chain is "dominated" by a strong one.
    vector<WeakFilterMeta> weakMetas;              ///< Per-candidate [qs,qe) + occ + score.
    vector<size_t> weakIdxWorkspace;               ///< Indices of weak candidates.
    vector<size_t> strongIdxWorkspace;             ///< Indices of strong candidates.
    vector<uint8_t> suppressWorkspace;             ///< 1 = suppress this candidate.
    vector<uint64_t> strongAnchorBegin;            ///< Start offset in flat anchor array.
    vector<uint64_t> strongAnchorEnd;              ///< End offset in flat anchor array.
    vector<uint32_t> strongAnchorStartsFlat;       ///< Sorted posA values for binary search.

    // ---- Stage 6: Multi-Copy Extraction Workspace ----
    // Used by applyMcopyFastSelection to extract multiple independent chains.
    vector<uint8_t> mcopyNodeUsed;                 ///< 1 = this anchor is already claimed.
    vector<int32_t> mcopyPathNodes;                ///< Temporary backtrack path buffer.
    vector<ChainCandidate> mcopySelectedCandidates; ///< Accepted secondary chains.

    // ---- Miscellaneous scratch ----
    struct ChainInterval { uint32_t qs; uint32_t qe; }; ///< Query-span of an accepted chain.
    vector<ChainInterval> acceptedIntervalsSame;   ///< Used in gap-based interval analysis.
    vector<ChainInterval> acceptedIntervalsDiff;
    vector<uint32_t> currentChainPath;             ///< Backtracked anchor indices for current chain.

    /// @brief Reset all vectors without releasing memory (amortized O(1) per read).
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

// ============================================================================
// HIGH-FREQUENCY MARKER DOWNSAMPLING SUPPORT
// ============================================================================
// Shared data structures and logic for downsampling consecutive high-frequency
// markers. Used by both the discovery path (all-vs-all) and the PAF path
// (pair-specific) to avoid duplicating the streak selection algorithm.
//
// Hifiasm reference: sketch.cpp hf_select + select_mz_h
//
// The algorithm buffers consecutive high-frequency markers into a "streak",
// then selects up to `keep = round(span / sampleDistance)` markers (capped at
// MAX_MAX_HIGH_OCC=16) using a (count, hash) ordering that prefers rarer,
// more informative markers. Markers whose occurrence count exceeds the streak
// span are discarded entirely (hifiasm's `rid < pe - ps` guard).
// ============================================================================

/// Buffered high-frequency marker awaiting streak-level downsampling.
struct PendingHighFrequencyMarker {
    uint64_t startIdx = 0;   ///< Start offset in the compact occurrence array.
    uint32_t count = 0;      ///< Genome-wide occurrence count of this k-mer.
    uint64_t hashKey = 0;    ///< Deterministic tie-break key (yak_hash64_64).
    uint32_t posA = 0;       ///< Base position on the query read.
    uint32_t ordinalA = 0;   ///< Marker ordinal on the query read.
    uint32_t weight = 1;     ///< Frequency-based weight (packed into cnt>>8).
    uint8_t isRcA = 0;       ///< 1 if the observed k-mer is the RC of canonical.
};

/// @brief Downsample a streak of consecutive high-frequency markers.
///
/// Selects up to `keep = round(span / sampleDistance)` markers from the streak,
/// preferring those with the lowest occurrence count (ties broken by hash).
/// Markers with `count >= span` are discarded (too repetitive for the local
/// context). When `keep == 0` (streak too short), all markers are discarded
/// (hifiasm parity: select_mz_h skips hf_select when max_high_occ==0).
///
/// @param streak           Buffered high-frequency markers.
/// @param workspace        Reusable index buffer (avoids allocation).
/// @param lastNonHighPos   Position of the last non-high-freq marker (-1 if none).
/// @param rightBoundaryPos Right boundary of the streak (next non-high pos or read end).
/// @param sampleDistance    Spacing between retained markers (hifiasm: 500bp).
/// @param maxPerStreak     Hard cap on retained markers (hifiasm: MAX_MAX_HIGH_OCC=16).
/// @param emitFn           Callback invoked for each retained marker.
template<typename EmitFn>
static inline void flushHighFrequencyStreak(
    vector<PendingHighFrequencyMarker>& streak,
    vector<size_t>& workspace,
    int64_t lastNonHighPos,
    uint32_t rightBoundaryPos,
    uint32_t sampleDistance,
    uint32_t maxPerStreak,
    EmitFn&& emitFn)
{
    if (streak.empty()) return;

    const uint32_t leftPos = (lastNonHighPos >= 0) ? uint32_t(lastNonHighPos) : 0U;
    const uint32_t span = (rightBoundaryPos > leftPos) ? (rightBoundaryPos - leftPos) : 0U;
    uint32_t keep = uint32_t(double(span) / double(sampleDistance) + 0.499);
    if (keep > maxPerStreak) keep = maxPerStreak;

    // Hifiasm parity: when keep==0 (streak shorter than ~sampleDistance/2),
    // all high-freq markers in this streak are discarded. In hifiasm's
    // select_mz_h, max_high_occ==0 means hf_select is never called, so
    // the markers keep their non-zero rid and get squeezed out.
    if (keep == 0) {
        streak.clear();
        return;
    }

    const size_t selectedCount = std::min<size_t>(streak.size(), keep);
    if (selectedCount >= streak.size()) {
        for (const auto& m : streak) emitFn(m, span);
        streak.clear();
        return;
    }

    // Select top-k by (count asc, hashKey asc) using nth_element (O(n) average).
    // This matches hifiasm's heap-based selection in hf_select.
    workspace.clear();
    workspace.reserve(streak.size());
    for (size_t i = 0; i < streak.size(); ++i) workspace.push_back(i);

    std::nth_element(
        workspace.begin(),
        workspace.begin() + static_cast<ptrdiff_t>(selectedCount),
        workspace.end(),
        [&](size_t a, size_t b) {
            if (streak[a].count != streak[b].count) return streak[a].count < streak[b].count;
            return streak[a].hashKey < streak[b].hashKey;
        });
    workspace.resize(selectedCount);

    // Emit in posA order so downstream radix sort sees monotonic query positions.
    std::sort(workspace.begin(), workspace.end(),
        [&](size_t a, size_t b) { return streak[a].posA < streak[b].posA; });

    for (const size_t idx : workspace) emitFn(streak[idx], span);
    streak.clear();
}

// ============================================================================
// K-MER UTILITY FUNCTIONS
// K-mer canonicalization and hashing: defined in InvertedIndexBuilder.hpp.
using dinara::getRcKmerId;
using dinara::hashKmer;

/**
 * @brief Hifiasm's yak_hash64_64 bijective integer hash (htab.h:150).
 *
 * A Thomas Wang-style hash used in hifiasm's minimizer sketching.
 * In Dinara, we use it as a deterministic tie-breaking key when
 * downsampling high-frequency marker streaks: among k-mers with
 * the same occurrence count, we keep the one whose hash is smallest
 * to ensure both Dinara and hifiasm select the same representative.
 *
 * Properties: bijective (no collisions), fast, good avalanche.
 */
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

/**
 * @brief Fold a wide KmerId (128-bit or larger) into a single uint64_t.
 *
 * XOR-folds the upper and lower 64-bit halves together. This is used
 * only for ordering and tie-breaking (not for hash-table probing),
 * so collision probability is acceptable.
 */
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

// Strict port of hifiasm 0.25.0-r726 (ec9a8b2) chaining pipeline.
// Variable names and control flow match hifiasm for step-parity verification.
// Coordinates use inclusive ends internally (hifiasm convention).
//
// Source functions:
//   Hash_Table.cpp:1475  cal_bw
//   Hash_Table.cpp:1515  comput_sc_ch_ec
//   Hash_Table.cpp:779   get_chainLen
//   Hash_Table.cpp:1752  push_ovlp_chain_qgen
//   Hash_Table.cpp:2007  quick_ck_lchain
//   Hash_Table.cpp:2097  lchain_qdp_mcopy_fast
//   anchor.cpp:86        ha_ov_type
//   anchor.cpp:1920      lchain_qgen_mcopy_fast (postfilter)
namespace {

/// Downweight score by k-mer occurrence frequency.
/// @reference Hifiasm Hash_Table.cpp:20 (normal_w macro)
static inline int32_t hifiasm_normal_w(int32_t x, int32_t y) noexcept
{
    return (x >= y) ? (x / y) : 1;
}

// ============================================================================
// HIFIASM DATA STRUCTURES
// ============================================================================

/**
 * @brief K-mer anchor hit for DP chaining (hifiasm k_mer_hit parity).
 *
 * Represents a shared k-mer between query and target sequences. The DP
 * chaining algorithm connects these anchors to form high-scoring chains.
 *
 * ## Bit-Packing Scheme (cnt field):
 *   cnt = (occurrence_weight << 8) | k-mer_span
 *   - Low 8 bits: k-mer span (typically k-mer length)
 *   - High 24 bits: occurrence weight (for downweighting repeats)
 *
 * ## Coordinate System:
 *   - self_offset: Query position (always forward strand)
 *   - offset: Target position (RC-transformed if strand=1)
 *   - All positions are inclusive endpoints
 *
 * @reference Hifiasm Hash_Table.h:117-121 (k_mer_hit)
 */
struct HifiasmKmerHit {
    // --- Core DP fields (must match hifiasm layout) ---
    uint32_t readID = 0;       ///< Target read ID (y_id in hifiasm)
    uint32_t self_offset = 0;  ///< Query position (x coordinate)
    uint32_t offset = 0;       ///< Target position (y coordinate, RC if strand=1)
    uint32_t cnt = 0;          ///< Bit-packed: (weight<<8) | span
    uint8_t strand = 0;        ///< Strand: 0=same, 1=opposite (reverse-complement)

    // --- Dinara extensions (not used by DP, for Alignment construction) ---
    uint32_t ordinalA = 0;     ///< Query marker ordinal (for Dinara output)
    uint32_t ordinalB = 0;     ///< Target marker ordinal (for Dinara output)
    uint32_t globalIndex = 0;  ///< Index into global hit array (chain persistence)

    // Helper accessors for bit-packed cnt field
    [[nodiscard]] uint8_t span() const noexcept { return static_cast<uint8_t>(cnt & 0xFFu); }
    [[nodiscard]] uint32_t weight() const noexcept { return cnt >> 8; }
};

/**
 * @brief Sort k-mer hits by (self_offset, offset) for DP chaining.
 *
 * The DP chaining algorithm requires hits to be sorted in ascending order
 * by (self_offset, offset). Since hits are already sorted by self_offset
 * (inherited from posA sorting), we only need to sort ties by offset.
 *
 * ## Algorithm:
 *   1. Find runs of hits with identical self_offset
 *   2. Within each run, sort by offset (target position)
 *   3. Move to next run
 *
 * ## Complexity:
 *   - Best case: O(n) if no ties (already sorted)
 *   - Worst case: O(n log k) where k is max run length
 *   - Typical: O(n) since ties are rare with large genomes
 *
 * @param[in,out] hits  K-mer hit array (modified in-place)
 * @complexity O(n) expected, O(n log k) worst case
 */
static inline void sortHifiasmHitsBySelfOffsetThenOffsetRuns(
    vector<HifiasmKmerHit>& hits)
{
    size_t begin = 0;
    const size_t n = hits.size();

    while (begin < n) {
        // Find end of run with same self_offset
        const uint32_t self_offset_val = hits[begin].self_offset;
        size_t end = begin + 1;

        while (end < n && hits[end].self_offset == self_offset_val) {
            ++end;
        }

        // Sort run by offset if there are ties (run length > 1)
        if (end - begin > 1) {
            std::sort(
                hits.begin() + static_cast<ptrdiff_t>(begin),
                hits.begin() + static_cast<ptrdiff_t>(end),
                [](const HifiasmKmerHit& a, const HifiasmKmerHit& b) noexcept {
                    return a.offset < b.offset;
                });
        }

        begin = end;
    }
}

/**
 * @brief Overlap region representing a chained alignment (hifiasm overlap_region parity).
 *
 * Stores the result of DP chaining: a high-scoring path through anchor hits
 * connecting two reads. Coordinates are normalized (left-aligned, right-extended)
 * for consistent downstream processing.
 *
 * ## Coordinate System:
 *   - **Normalized coordinates** (x_pos_*, y_pos_*):
 *     * Left-normalized: one sequence starts at 0
 *     * Right-extended: extends to sequence boundaries
 *     * Inclusive endpoints [start, end]
 *
 * ## Strand Convention:
 *   - x_pos_strand: Always 0 (query is forward)
 *   - y_pos_strand: 0=same strand, 1=reverse-complement
 *
 * @reference Hifiasm Hash_Table.h:78-106 (overlap_region)
 */
struct HifiasmOverlapRegion {
    // --- Core overlap metadata ---
    uint32_t x_id = 0;          ///< Query read ID
    uint32_t y_id = 0;          ///< Target read ID
    uint8_t x_pos_strand = 0;   ///< Query strand (always 0)
    uint8_t y_pos_strand = 0;   ///< Target strand (0=fwd, 1=rev)

    // --- Normalized overlap coordinates (inclusive, extended) ---
    uint32_t x_pos_s = 0;       ///< Query start (left-normalized)
    uint32_t x_pos_e = 0;       ///< Query end (right-extended)
    uint32_t y_pos_s = 0;       ///< Target start (left-normalized, RC if strand=1)
    uint32_t y_pos_e = 0;       ///< Target end (right-extended, RC if strand=1)

    // --- Scoring and chain metadata ---
    int32_t shared_seed = 0;    ///< Chain DP score (sum of anchor scores)
    uint32_t align_length = 0;  ///< Number of anchors in chain

    /**
     * Dual-purpose field (hifiasm compatibility):
     *   - Hifiasm: Offset into k_mer_hit array (chain anchor indices)
     *   - Dinara: Offset into chainHitIndexFlat vector
     * Combined with align_length, defines the chain hit subarray.
     */
    uint32_t non_homopolymer_errors = 0;

};

/**
 * @brief Reusable scratch space for DP chaining computation (hifiasm Chain_Data parity).
 *
 * This structure holds temporary arrays for the chaining DP algorithm.
 * Reusing the same scratch space across reads avoids repeated allocations
 * and improves cache locality.
 *
 * ## Array Purposes:
 *   - **score (f[])**: DP scores - best score ending at anchor i
 *   - **pre (p[])**: Predecessor indices - backtrack pointers for chain reconstruction
 *   - **tmp (t[])**: Multi-purpose workspace:
 *     * During DP: Visit marks for max_skip pruning
 *     * During mcopy: Heap keys for chain extraction
 *   - **occ (ii[])**: Multi-purpose workspace:
 *     * During backtrack: Chain membership flags
 *     * During mcopy: Index buffer for chain anchors
 *
 * ## Memory Layout:
 *   All arrays use Structure-of-Arrays (SoA) layout for cache efficiency.
 *   The DP inner loop accesses these arrays sequentially, benefiting from
 *   hardware prefetching and cache line alignment.
 *
 * @reference Hifiasm Hash_Table.h:122-131 (Chain_Data)
 */
struct HifiasmChainDataScratch {
    vector<int32_t> score; ///< f[i]: Best DP score ending at anchor i
    vector<int64_t> pre;   ///< p[i]: Predecessor index for backtracking
    vector<int64_t> tmp;   ///< t[i]: Visit marks (DP) or heap keys (mcopy)
    vector<int32_t> occ;   ///< ii[i]: Chain flags (backtrack) or indices (mcopy)

    /**
     * @brief Resize all arrays to accommodate n anchors.
     * @complexity O(n) if reallocation needed, O(1) amortized
     */
    void resize(const size_t n) {
        score.resize(n);
        pre.resize(n);
        tmp.resize(n);
        occ.resize(n);
    }
};

/// Compute adaptive bandwidth: the maximum allowed diagonal deviation between
/// two anchors, scaled by the alignable span of the overlap region.
/// Near sequence ends the bandwidth is larger (more room for indels).
/// @reference Hifiasm Hash_Table.cpp:1475-1488 (cal_bw)
static inline int32_t hifiasm_cal_bw(
    const HifiasmKmerHit* ai,
    const HifiasmKmerHit* aj,
    double bw_rate,
    int64_t sf_l,
    int64_t ot_l) noexcept
{
    int64_t sf_s = static_cast<int64_t>(aj->self_offset);
    int64_t sf_e = static_cast<int64_t>(ai->self_offset) + 1;
    int64_t ot_s = static_cast<int64_t>(aj->offset);
    int64_t ot_e = static_cast<int64_t>(ai->offset) + 1;
    int64_t sf_r = sf_l - sf_e;
    int64_t ot_r = ot_l - ot_e;

    // Left-normalize: shift query start to 0 or offset by target start.
    if (sf_s <= ot_s) sf_s = 0;
    else sf_s -= ot_s;

    // Right-extend: extend query end to sequence boundary or by target remainder.
    if (sf_r <= ot_r) sf_e = sf_l;
    else sf_e += ot_r;

    return static_cast<int32_t>(static_cast<double>(sf_e - sf_s) * bw_rate);
}

/// Compute chaining score between two anchors with gap penalties.
///
/// Returns INT32_MIN if the connection is invalid (monotonicity violation or
/// bandwidth exceeded). Otherwise returns a score where higher is better:
///   score = min(span, min(dq, dr)) / weight - gap_penalty - skip_penalty
///
/// The penalty uses a split strategy: small gaps (dd < 4) take the milder of
/// linear vs adaptive penalty; large gaps take the harsher one.
///
/// @reference Hifiasm Hash_Table.cpp:1515-1541 (comput_sc_ch_ec)
static inline int32_t hifiasm_comput_sc_ch_ec(
    const HifiasmKmerHit* ai,
    const HifiasmKmerHit* aj,
    double bw_rate,
    double chn_pen_gap,
    double chn_pen_skip,
    int64_t sl,
    int64_t ol) noexcept
{
    // Monotonicity: both query and target must advance.
    const int32_t dq = static_cast<int32_t>(
        static_cast<int64_t>(ai->self_offset) - static_cast<int64_t>(aj->self_offset));
    if (dq <= 0) return INT32_MIN;

    const int32_t dr = static_cast<int32_t>(
        static_cast<int64_t>(ai->offset) - static_cast<int64_t>(aj->offset));
    if (dr <= 0) return INT32_MIN;

    // Diagonal deviation (gap difference).
    const int32_t dd = std::abs(dr - dq);
    if ((dd > 16) && (dd > hifiasm_cal_bw(ai, aj, bw_rate, sl, ol)))
        return INT32_MIN;

    // Base score: min of k-mer span and effective step, normalized by weight.
    const int32_t dg = std::min(dr, dq);
    const int32_t q_span = static_cast<int32_t>(ai->cnt & 0xFFu);
    int32_t sc = std::min(q_span, dg);
    sc = hifiasm_normal_w(sc, static_cast<int32_t>(ai->cnt >> 8));

    // Gap penalty: linear + adaptive, with split strategy by gap size.
    if (dd || (dg > q_span && dg > 0)) {
        double lin_pen = chn_pen_gap * static_cast<double>(dd);
        const double a_pen =
            static_cast<double>(sc) * (static_cast<double>(dd) / static_cast<double>(dg)) / bw_rate;

        if (dd < 4) lin_pen = (lin_pen > a_pen) ? a_pen : lin_pen;
        else         lin_pen = (lin_pen < a_pen) ? a_pen : lin_pen;

        lin_pen += chn_pen_skip * static_cast<double>(dg);
        sc -= static_cast<int32_t>(lin_pen);
    }

    return sc;
}

/// Minimap2-style chaining score between two anchors.
///
/// Uses fixed bandwidth and logarithmic gap penalty instead of hifiasm's
/// adaptive bandwidth and linear/adaptive penalty.
///
/// Score = min(q_span, dg) - gap_penalty
/// gap_penalty = 0.01 * span * dd + 0.5 * log2(dd + 1)
///
/// Returns INT32_MIN if connection is invalid (monotonicity violation,
/// bandwidth exceeded, or max gap exceeded).
///
/// @reference minimap2 lchain.c:mg_chain_dp_score (simplified for sr preset)
static inline int32_t minimap2_comput_sc(
    const HifiasmKmerHit* ai,
    const HifiasmKmerHit* aj,
    int32_t bw,             // Fixed bandwidth (minimap2 -r).
    int32_t maxGap,         // Max gap (minimap2 -g).
    double chn_pen_gap,     // Gap penalty coefficient.
    double chn_pen_skip     // Skip penalty coefficient.
) noexcept
{
    // Monotonicity: both query and target must advance.
    const int64_t dq = static_cast<int64_t>(ai->self_offset) - static_cast<int64_t>(aj->self_offset);
    if(dq <= 0) return INT32_MIN;

    const int64_t dr = static_cast<int64_t>(ai->offset) - static_cast<int64_t>(aj->offset);
    if(dr <= 0) return INT32_MIN;

    // Max gap check.
    if(dq > maxGap || dr > maxGap) return INT32_MIN;

    // Diagonal deviation.
    const int64_t dd = (dr > dq) ? (dr - dq) : (dq - dr);

    // Fixed bandwidth check.
    if(dd > bw) return INT32_MIN;

    // Base score.
    const int64_t dg = std::min(dq, dr);
    const int32_t q_span = static_cast<int32_t>(ai->cnt & 0xFFu);
    int32_t sc = std::min(q_span, static_cast<int32_t>(dg));
    sc = hifiasm_normal_w(sc, static_cast<int32_t>(ai->cnt >> 8));

    // Minimap2 gap penalty: linear + logarithmic.
    // minimap2 uses: pen = 0.01 * avg_qspan * dd + log2(dd) / 2
    // chn_pen_gap is used as the linear coefficient (0.01 * avg_qspan).
    if(dd > 0) {
        const double lin_pen = chn_pen_gap * static_cast<double>(dd);
        const double log_pen = 0.5 * log2(static_cast<double>(dd) + 1.0);
        sc -= static_cast<int32_t>(lin_pen + log_pen);
    }

    return sc;
}

/// Compute effective overlap length after left-normalizing and right-extending
/// the chain coordinates to sequence boundaries. Used for tie-breaking when
/// chains have identical scores (prefer shorter effective length = tighter fit).
/// @reference Hifiasm Hash_Table.cpp:779-809 (get_chainLen)
static inline int64_t hifiasm_get_chainLen(
    int64_t x_beg, int64_t x_end, int64_t xLen,
    int64_t y_beg, int64_t y_end, int64_t yLen) noexcept
{
    if (x_beg <= y_beg) { y_beg -= x_beg; x_beg = 0; }
    else                { x_beg -= y_beg; y_beg = 0; }

    const int64_t xr = xLen - x_end - 1;
    const int64_t yr = yLen - y_end - 1;
    if (xr <= yr) { x_end = xLen - 1; y_end += xr; }
    else          { x_end += yr; y_end = yLen - 1; }

    return x_end - x_beg + 1;
}

/// Build a normalized overlap region from the first and last anchors of a chain.
/// Coordinates are left-normalized (one start becomes 0) and right-extended
/// (one end reaches its sequence boundary).
/// @reference Hifiasm Hash_Table.cpp:1752-1780 (push_ovlp_chain_qgen)
static inline void hifiasm_push_ovlp_chain_qgen(
    HifiasmOverlapRegion& o,
    uint32_t xid,
    int64_t xl,
    int64_t yl,
    int64_t sc,
    const HifiasmKmerHit* beg,
    const HifiasmKmerHit* end) noexcept
{
    o.x_id = xid;
    o.y_id = beg->readID;
    o.x_pos_strand = 0;
    o.y_pos_strand = beg->strand;

    o.x_pos_s = beg->self_offset;
    o.y_pos_s = beg->offset;
    o.x_pos_e = end->self_offset;
    o.y_pos_e = end->offset;

    // Left-normalize: shift so the smaller start becomes 0.
    if (o.x_pos_s <= o.y_pos_s) { o.y_pos_s -= o.x_pos_s; o.x_pos_s = 0; }
    else                        { o.x_pos_s -= o.y_pos_s; o.y_pos_s = 0; }

    // Right-extend: extend until one sequence reaches its boundary.
    const int64_t xr = xl - static_cast<int64_t>(o.x_pos_e) - 1;
    const int64_t yr = yl - static_cast<int64_t>(o.y_pos_e) - 1;
    if (xr <= yr) { o.x_pos_e = static_cast<uint32_t>(xl - 1); o.y_pos_e += static_cast<uint32_t>(xr); }
    else          { o.y_pos_e = static_cast<uint32_t>(yl - 1); o.x_pos_e += static_cast<uint32_t>(yr); }

    o.shared_seed = static_cast<int32_t>(sc);
    o.align_length = 0;
    o.non_homopolymer_errors = 0;
}

/// O(N) quick-check: greedily chain monotonic same-strand segments.
///
/// Scans the sorted anchor array for maximal contiguous runs where both
/// coordinates strictly increase. Each such run is scored greedily (each
/// anchor chains to its immediate predecessor). If a run spans an entire
/// strand-homogeneous segment and passes cumulative bandwidth validation,
/// it narrows the [si, ei) range that the full DP must cover.
///
/// @param[out] msc, msc_i  Best chain score and its terminal anchor index.
/// @param[out] movl        Effective overlap length of best chain (tie-breaking).
/// @param[out] plus        Minimum accumulated score (quality floor for mcopy).
/// @param[out] si, ei      Unsolved DP range — full DP only runs over [si, ei).
/// @reference Hifiasm Hash_Table.cpp:2007-2095 (quick_ck_lchain)
static inline void hifiasm_quick_ck_lchain(
    HifiasmKmerHit* a, int64_t a_n,
    int64_t xl, int64_t yl,
    double chn_pen_gap, double chn_pen_skip, double bw_rate,
    int64_t* p, int64_t* t, int32_t* f, int32_t* ii,
    int64_t* plus, int64_t* msc, int64_t* msc_i, int64_t* movl,
    int64_t* si, int64_t* ei)
{
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

    // Scan for strand-homogeneous segments. At each boundary (strand change or
    // end-of-array), attempt greedy chaining if the segment was monotonic.
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
                f[l] = static_cast<int32_t>(a[l].cnt & 0xffu);
                if (f[l] >= msc0) { msc0 = f[l]; msc_i0 = l; }
                if (f[l] < plus0) plus0 = f[l];

                // Greedy forward scan: link each anchor to its immediate predecessor.
                for (z = l + 1; z < k; z++) {
                    ai = &a[z];
                    aj = &a[z - 1];

                    dq = static_cast<int64_t>(ai->self_offset) - static_cast<int64_t>(aj->self_offset);
                    if (dq <= 0) break;
                    dr = static_cast<int64_t>(ai->offset) - static_cast<int64_t>(aj->offset);
                    if (dr <= 0) break;

                    dd = (dr > dq) ? (dr - dq) : (dq - dr);
                    if ((dd > 16) && (dd > hifiasm_cal_bw(ai, aj, bw_rate, xl, yl))) break;

                    dg = (dr < dq) ? dr : dq;
                    q_span = static_cast<int64_t>(ai->cnt & 0xffu);
                    sc = (q_span < dg) ? q_span : dg;
                    sc = hifiasm_normal_w(static_cast<int32_t>(sc), static_cast<int32_t>(ai->cnt >> 8));

                    if (dd || (dg > q_span && dg > 0)) {
                        lin_pen = chn_pen_gap * static_cast<double>(dd);
                        a_pen = static_cast<double>(sc) * (static_cast<double>(dd) / static_cast<double>(dg)) / bw_rate;
                        if (dd < 4) lin_pen = (lin_pen > a_pen) ? a_pen : lin_pen;
                        else        lin_pen = (lin_pen < a_pen) ? a_pen : lin_pen;
                        lin_pen += chn_pen_skip * static_cast<double>(dg);
                        sc -= static_cast<int32_t>(lin_pen);
                    }

                    sc += f[z - 1];
                    csc = static_cast<int64_t>(a[z].cnt & 0xffu);
                    if (sc < csc) break;  // Starting fresh is better.

                    p[z] = z - 1;
                    f[z] = static_cast<int32_t>(sc);
                    ddt += dd;
                    if (f[z] >= msc0) { msc0 = f[z]; msc_i0 = z; }
                    if (f[z] < plus0) plus0 = f[z];
                }

                // Segment fully chained and best anchor is the last one.
                if ((z >= k) && (msc_i0 == (k - 1))) {
                    // Invalidate if cumulative diagonal drift exceeds bandwidth.
                    if ((k - l >= 2) && (ddt > 16) && (ddt > hifiasm_cal_bw(&a[k - 1], &a[l], bw_rate, xl, yl)))
                        msc_i0 = INT32_MIN;

                    if (msc_i0 == (k - 1)) {
                        if (msc0 >= (*msc)) {
                            movl0 = hifiasm_get_chainLen(
                                a[msc_i0].self_offset, a[msc_i0].self_offset, xl,
                                a[msc_i0].offset, a[msc_i0].offset, yl);
                            if (msc0 > (*msc) || movl0 < (*movl)) {
                                *msc = msc0; *msc_i = msc_i0; *movl = movl0;
                            }
                        }
                        if (plus0 < (*plus)) *plus = plus0;

                        // Narrow the unsolved DP range.
                        if ((*ei) > k) (*si) = k;
                        else           (*ei) = l;
                    }
                }
            }
            l = k;
            is_srt = 1;
        } else {
            if ((a[k].self_offset <= a[k - 1].self_offset) || (a[k].offset <= a[k - 1].offset))
                is_srt = 0;
            t[k - 1] = 0;
            ii[k - 1] = 0;
        }
    }
}

/// Main O(N * max_iter) DP loop for anchor chaining.
///
/// Computes f[i] = max_{j < i, same_strand} (score(i,j) + f[j]) with three
/// pruning strategies: lookback window (max_iter), max-skip (max_skip), and
/// sliding distance window (max_ii). Updates msc/msc_i/movl/plus in place.
///
/// chainingMode: 0 = hifiasm scoring, 1 = minimap2-sr scoring.
///
/// @reference Hifiasm Hash_Table.cpp:2097-2170 (main DP loop in lchain_qdp_mcopy_fast)
static inline void run_main_dp_loop(
    HifiasmKmerHit* a, int64_t si, int64_t ei,
    int64_t max_skip, int64_t max_iter, int64_t max_dis,
    double chn_pen_gap, double chn_pen_skip, double bw_rate,
    int64_t xl, int64_t yl,
    int64_t* p, int32_t* f, int64_t* t, int32_t* ii,
    int64_t* msc, int64_t* msc_i, int64_t* movl, int64_t* plus,
    int chainingMode = 0,
    int32_t minimap2Bw = 100,
    int32_t minimap2MaxGap = 100) noexcept
{
    int64_t max_f, n_skip, st, max_j, end_j, sc, max_ii = -1, ovl;
    int32_t max, tmp;

    // Lambda to compute score between two anchors using the selected mode.
    auto computeScore = [&](const HifiasmKmerHit* ai, const HifiasmKmerHit* aj) -> int32_t {
        if(chainingMode == 1) {
            return minimap2_comput_sc(ai, aj, minimap2Bw, minimap2MaxGap, chn_pen_gap, chn_pen_skip);
        } else {
            return hifiasm_comput_sc_ch_ec(ai, aj, bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);
        }
    };

    int64_t i, j;
    for (i = st = si; i < ei; ++i) {
        max_f = a[i].cnt & 0xFFu;
        n_skip = 0;
        max_j = end_j = -1;

        if ((i - st) > max_iter) st = i - max_iter;
        while (a[i].strand != a[st].strand) ++st;

        // Inner loop: find best predecessor j for anchor i.
        for (j = i - 1; j >= st; --j) {
            sc = computeScore(&a[i], &a[j]);
            if (sc == INT32_MIN) continue;
            sc += f[j];

            if (sc > max_f) {
                max_f = sc; max_j = j;
                if (n_skip > 0) --n_skip;
            } else if (t[j] == static_cast<int32_t>(i)) {
                if (++n_skip > max_skip) break;
            }
            if (p[j] >= 0) t[p[j]] = i;
        }
        end_j = j;

        // Sliding distance window: cache best score within max_dis.
        if ((max_ii < 0) ||
            (a[i].self_offset > a[max_ii].self_offset + max_dis) ||
            (a[i].strand != a[max_ii].strand)) {
            max = INT32_MIN; max_ii = -1;
            for (int64_t j = i - 1;
                 (j >= st) && (a[i].self_offset <= max_dis + a[j].self_offset) && (a[i].strand == a[j].strand);
                 --j) {
                if (max < f[j]) { max = f[j]; max_ii = j; }
            }
        }

        // Fallback: try max_ii if it was skipped by max-skip pruning.
        if ((max_ii >= 0) && (max_ii < end_j) && (a[i].strand == a[max_ii].strand)) {
            tmp = computeScore(&a[i], &a[max_ii]);
            if (tmp != INT32_MIN && max_f < tmp + f[max_ii]) {
                max_f = tmp + f[max_ii]; max_j = max_ii;
            }
        }

        f[i] = static_cast<int32_t>(max_f);
        p[i] = max_j;

        if ((max_ii < 0) ||
            ((a[i].self_offset <= max_dis + a[max_ii].self_offset) &&
             (a[i].strand == a[max_ii].strand) && (f[max_ii] < f[i])))
            max_ii = i;

        // Track global best chain (tie-break by overlap length).
        if (f[i] >= (*msc)) {
            ovl = hifiasm_get_chainLen(a[i].self_offset, a[i].self_offset, xl, a[i].offset, a[i].offset, yl);
            if (f[i] > (*msc) || ovl < (*movl)) {
                *msc = f[i]; *msc_i = i; *movl = ovl;
            }
        }
        if (f[i] < (*plus)) *plus = f[i];
        ii[i] = 0;
    }
}

/// Follow predecessor pointers from msc_i to reconstruct the best chain.
/// Stores indices in t[] in reverse order (t[0] = last anchor, t[cL-1] = first).
/// Marks chain membership in ii[] (1 = in chain).
/// @return Chain length cL.
static inline int64_t backtrack_best_chain(
    int64_t msc_i, int64_t* p, int64_t* t, int32_t* ii) noexcept
{
    int64_t cL = 0;
    for (int64_t i = msc_i; i >= 0; i = p[i]) {
        ii[i] = 1;
        t[cL++] = i;
    }
    return cL;
}

/// Emit the best chain as a HifiasmOverlapRegion with normalized coordinates.
/// t[0..cL-1] holds chain indices in reverse order (from backtrack_best_chain).
/// Chain hit global indices are appended to chainHitIndexFlat in forward order.
static inline void emit_best_chain_as_overlap(
    HifiasmKmerHit* a, int64_t* t, int64_t cL, int64_t msc,
    uint32_t xid, int64_t xl, int64_t yl,
    vector<HifiasmOverlapRegion>& res,
    vector<uint32_t>& chainHitIndexFlat) noexcept
{
    if (cL <= 0) return;

    HifiasmOverlapRegion z{};
    hifiasm_push_ovlp_chain_qgen(z, xid, xl, yl, msc, &a[t[cL - 1]], &a[t[0]]);
    z.align_length = static_cast<uint32_t>(cL);
    z.non_homopolymer_errors = static_cast<uint32_t>(chainHitIndexFlat.size());

    for (int64_t i = 0; i < cL; ++i)
        chainHitIndexFlat.push_back(a[t[cL - i - 1]].globalIndex);

    res.push_back(z);
}

/// Z-drop backtracking: walk back from a chain endpoint, stop when the
/// score drops by more than max_drop from the running maximum.
/// Returns the anchor index where the chain should be cut.
///
/// @reference minimap2 lchain.c:mg_chain_bk_end
static inline int64_t zdrop_chain_bk_end(
    int32_t max_drop,
    const int32_t* f,       // DP scores.
    const int64_t* p,       // Predecessor pointers.
    int32_t* t,             // Temporary marker array (0 = unvisited).
    int32_t k_score,        // Score at the chain endpoint.
    int64_t k_idx           // Anchor index of the chain endpoint.
) noexcept
{
    int64_t i = k_idx, end_i = -1, max_i = i;
    int32_t max_s = 0;
    if(i < 0 || t[i] != 0) return i;
    do {
        t[i] = 2;
        end_i = i = p[i];
        const int32_t s = (i < 0) ? k_score : k_score - f[i];
        if(s > max_s) { max_s = s; max_i = i; }
        else if(max_s - s > max_drop) break;
    } while(i >= 0 && t[i] == 0);
    // Reset temporary marks.
    for(i = k_idx; i >= 0 && i != end_i; i = p[i])
        t[i] = 0;
    return max_i;
}

/// Z-drop chain extraction: extract multiple chains from the DP table by
/// processing endpoints in descending score order. Each chain is cut at
/// Z-drop boundaries, naturally splitting at SV breakpoints.
///
/// @reference minimap2 lchain.c:mg_chain_backtrack
static inline void zdrop_chain_backtrack(
    HifiasmKmerHit* a,
    int64_t a_n,
    const int32_t* f,       // DP scores.
    const int64_t* p,       // Predecessor pointers.
    int32_t max_drop,       // Z-drop threshold (= bandwidth for minimap2-sr).
    int32_t min_cnt,        // Min anchors per chain (minimap2 -n).
    int32_t min_sc,         // Min chain score (minimap2 -m).
    int32_t max_chains,     // Max chains to extract (minimap2 -N). 0 = unlimited.
    uint32_t xid,
    int64_t xl, int64_t yl,
    vector<HifiasmOverlapRegion>& res,
    vector<uint32_t>& chainHitIndexFlat
)
{
    // Collect endpoints with score >= min_sc, sorted by score ascending.
    struct ScoreIdx { int32_t score; int64_t idx; };
    vector<ScoreIdx> z;
    z.reserve(size_t(a_n));
    for(int64_t i = 0; i < a_n; ++i) {
        if(f[i] >= min_sc) {
            z.push_back({f[i], i});
        }
    }
    if(z.empty()) return;

    // Sort ascending by score (we process from the end = descending).
    std::sort(z.begin(), z.end(), [](const ScoreIdx& a, const ScoreIdx& b) {
        return a.score < b.score;
    });

    // Temporary marker: 0 = available, 1 = claimed by a chain.
    vector<int32_t> t(size_t(a_n), 0);

    // Temporary storage for chain anchor indices (reversed order).
    vector<int32_t> v;
    v.reserve(size_t(a_n));

    // Process endpoints in descending score order.
    int32_t n_chains = 0;
    for(int64_t k = int64_t(z.size()) - 1; k >= 0; --k) {
        if(max_chains > 0 && n_chains >= max_chains) break;
        if(t[z[k].idx] != 0) continue;  // Already claimed.

        // Find where to cut this chain (Z-drop).
        int64_t end_i = zdrop_chain_bk_end(
            max_drop, f, p, t.data(), int32_t(z[k].score), z[k].idx);

        // Collect anchors from endpoint to cut point.
        size_t n_v0 = v.size();
        for(int64_t i = z[k].idx; i != end_i && i >= 0; i = p[i]) {
            if(t[i] != 0) break;  // Hit a claimed anchor.
            v.push_back(int32_t(i));
            t[i] = 1;
        }

        size_t chain_len = v.size() - n_v0;
        if(chain_len == 0) continue;

        // Compute chain score (score at endpoint minus score at cut point).
        int32_t sc;
        int64_t last_i = v.back();
        int64_t pred_of_last = p[last_i];
        sc = (pred_of_last < 0 || t[pred_of_last] != 0)
            ? f[z[k].idx]
            : f[z[k].idx] - f[pred_of_last];

        // Filter by min score and min anchor count.
        if(sc >= min_sc && int32_t(chain_len) >= min_cnt) {
            // Emit chain as overlap region.
            // Anchors in v are in reverse order (endpoint first).
            int64_t first_anchor = v[n_v0 + chain_len - 1];
            int64_t last_anchor = v[n_v0];

            HifiasmOverlapRegion region{};
            hifiasm_push_ovlp_chain_qgen(
                region, xid, xl, yl, int64_t(sc),
                &a[first_anchor], &a[last_anchor]);

            region.align_length = uint32_t(chain_len);
            region.non_homopolymer_errors = uint32_t(chainHitIndexFlat.size());

            // Store hit indices in forward order.
            for(size_t j = 0; j < chain_len; ++j) {
                chainHitIndexFlat.push_back(
                    a[v[n_v0 + (chain_len - 1 - j)]].globalIndex);
            }
            res.push_back(region);
            ++n_chains;
        } else {
            // Reject chain — unmark anchors.
            for(size_t j = n_v0; j < v.size(); ++j) {
                t[v[j]] = 0;
            }
            v.resize(n_v0);
        }
    }
}

/// Full DP chaining pipeline with optional quick-check and multi-copy extraction.
///
/// Phase 1: Optional O(N) quick-check for collinear hits (narrows DP range).
/// Phase 2: O(N * max_iter) DP with max-skip and distance-window pruning.
/// Phase 3: Backtrack best chain from predecessor pointers.
/// Phase 4: Optional mcopy extraction (up to mcopy_num secondary chains).
///          In minimap2-sr mode (chainingMode=1), uses Z-drop backtracking instead.
/// Phase 5: Emit best chain as overlap region (fallback if mcopy disabled/failed).
///
/// @reference Hifiasm Hash_Table.cpp:2097-2284 (lchain_qdp_mcopy_fast)
static inline void hifiasm_lchain_qdp_mcopy_fast(
    vector<HifiasmKmerHit>& a,
    HifiasmChainDataScratch& dp,
    vector<HifiasmOverlapRegion>& res,
    vector<uint32_t>& chainHitIndexFlat,
    int64_t max_skip, int64_t max_iter, int64_t max_dis,
    double chn_pen_gap, double chn_pen_skip, double bw_rate,
    uint32_t xid, int64_t xl, int64_t yl,
    int64_t quick_check,
    int64_t mcopy_num, double mcopy_rate, int64_t mcopy_khit_cutoff,
    int chainingMode = 0,
    int32_t minimap2Bw = 100,
    int32_t minimap2MaxGap = 100,
    int32_t minimap2MinChainScore = 25)
{
    const int64_t a_n = static_cast<int64_t>(a.size());
    if (a_n <= 0) return;

    dp.resize(static_cast<size_t>(a_n));
    int64_t* p  = dp.pre.data();
    int64_t* t  = dp.tmp.data();
    int32_t* f  = dp.score.data();
    int32_t* ii = dp.occ.data();

    int64_t msc, msc_i, movl, plus = 0;
    int64_t min_sc, ch_n, si, ei;
    int64_t i, k, cL = 0, sc;

    // Phase 1: Quick check (hifiasm mode only — uses hifiasm-specific scoring).
    if (quick_check && chainingMode == 0) {
        hifiasm_quick_ck_lchain(
            a.data(), a_n, xl, yl, chn_pen_gap, chn_pen_skip, bw_rate,
            p, t, f, ii, &plus, &msc, &msc_i, &movl, &si, &ei);
    } else {
        msc = msc_i = INT32_MIN;
        movl = INT32_MAX;
        plus = 0;
        si = 0; ei = a_n;
        std::fill(t, t + a_n, int64_t(0));
    }

    // Phase 2: Main DP.
    run_main_dp_loop(
        a.data(), si, ei, max_skip, max_iter, max_dis,
        chn_pen_gap, chn_pen_skip, bw_rate, xl, yl,
        p, f, t, ii, &msc, &msc_i, &movl, &plus,
        chainingMode, minimap2Bw, minimap2MaxGap);

    // minimap2-sr mode: use Z-drop backtracking instead of mcopy.
    if(chainingMode == 1) {
        zdrop_chain_backtrack(
            a.data(), a_n, f, p,
            minimap2Bw,             // max_drop = bandwidth (minimap2 convention).
            int32_t(mcopy_khit_cutoff > 0 ? mcopy_khit_cutoff : 2),  // min_cnt.
            minimap2MinChainScore,  // min_sc (minimap2 -m25).
            int32_t(mcopy_num),     // max_chains (minimap2 -N).
            xid, xl, yl,
            res, chainHitIndexFlat);
        return;
    }

    // Phase 3: Backtrack best chain.
    cL = backtrack_best_chain(msc_i, p, t, ii);

    // Phase 4: Multi-copy extraction.
    if (mcopy_num > 1) {
        if (cL >= mcopy_khit_cutoff) {
            // Normalize scores to non-negative range for threshold comparison.
            msc -= plus;
            min_sc = static_cast<int64_t>(static_cast<double>(msc) * mcopy_rate);
            ii[msc_i] = 0;

            // Collect candidates: pack (score << 32 | index << 1) into t[].
            for (i = ch_n = 0; i < a_n; ++i) {
                f[i] -= static_cast<int32_t>(plus);
                if (i >= ch_n) t[i] = 0;
                if ((!ii[i]) && (f[i] >= min_sc)) {
                    t[ch_n] = (static_cast<int64_t>(static_cast<uint64_t>(static_cast<uint32_t>(f[i])) << 32)
                             | (static_cast<uint64_t>(i) << 1));
                    ch_n++;
                }
            }

            if (ch_n <= 1) {
                // Too few candidates — restore best chain for fallback emit.
                msc += plus;
                cL = 0;
                for (i = msc_i; i >= 0; i = p[i]) t[cL++] = i;
            }

            if (ch_n > 1) {
                int64_t n_v, n_v0, ni, n_u;

                std::sort(t, t + ch_n, [](int64_t a, int64_t b) noexcept {
                    return static_cast<uint64_t>(a) < static_cast<uint64_t>(b);
                });

                // Extract top mcopy_num chains in descending score order.
                for (k = ch_n - 1, n_v = n_u = 0; k >= 0 && n_u < mcopy_num; --k) {
                    n_v0 = n_v;
                    for (i = static_cast<int64_t>(static_cast<uint32_t>(t[k]) >> 1);
                         i >= 0 && (t[i] & 1) == 0; ) {
                        ii[n_v++] = static_cast<int32_t>(i);
                        t[i] |= 1;
                        i = p[i];
                    }
                    if (n_v0 == n_v) continue;

                    sc = (i < 0) ? static_cast<int64_t>(static_cast<uint32_t>(static_cast<uint64_t>(t[k]) >> 32))
                                 : static_cast<int64_t>(static_cast<uint32_t>(static_cast<uint64_t>(t[k]) >> 32)) - f[i];

                    if (sc >= min_sc) {
                        HifiasmOverlapRegion z{};
                        hifiasm_push_ovlp_chain_qgen(z, xid, xl, yl, sc + plus, &a[ii[n_v - 1]], &a[ii[n_v0]]);

                        if ((!n_u) || (n_v - n_v0 > 1)) {
                            z.align_length = static_cast<uint32_t>(n_v - n_v0);
                            z.non_homopolymer_errors = static_cast<uint32_t>(chainHitIndexFlat.size());
                            ni = z.align_length;
                            for (uint32_t j = 0; j < ni; ++j)
                                chainHitIndexFlat.push_back(a[ii[n_v0 + (ni - j - 1)]].globalIndex);
                            res.push_back(z);
                            n_u++;
                        } else {
                            n_v = n_v0;
                        }
                    } else {
                        n_v = n_v0;
                    }
                }
                return;  // hifiasm parity: always return after ch_n > 1 block.
            }
        }
    }

    // Phase 5: Emit best chain (fallback).
    emit_best_chain_as_overlap(a.data(), t, cL, msc, xid, xl, yl, res, chainHitIndexFlat);
}

/// Classify overlap type by query coordinate coverage.
///   0 = prefix  [0, ...]       1 = suffix  [..., len-1]
///   2 = contained [0, len-1]   3 = internal [>0, <len-1]
/// @reference Hifiasm anchor.cpp:86 (ha_ov_type)
static inline int hifiasm_ha_ov_type(const HifiasmOverlapRegion& r, uint32_t len) noexcept
{
    if (r.x_pos_s == 0 && r.x_pos_e == len - 1U) return 2;
    if (r.x_pos_s > 0 && r.x_pos_e < len - 1U)   return 3;
    return (r.x_pos_s == 0) ? 0 : 1;
}

/// Three-stage postfilter applied after DP chaining:
///
/// Stage 1 (max_n_chain): Per-overlap-type score cap. Sort by score descending,
///   record the score at the max_n_chain-th position for each type, reject below.
///
/// Stage 2 (COV_W): Coverage-window rescue for type-3 overlaps. Divides the query
///   into fixed-width windows; a below-threshold type-3 overlap is rescued if ≥70%
///   of its span falls in under-capacity windows.
///
/// Stage 3 (R485): Weak-chain suppression. Chains with < chain_cutoff anchors are
///   suppressed if a strong chain overlaps ≥95% of their span with ≥16× more
///   anchors and ≥16× higher score.
///
/// @reference Hifiasm anchor.cpp:1920-2100 (lchain_qgen_mcopy_fast postfilter)
static inline void hifiasm_lchain_qgen_mcopy_fast_postfilter(
    vector<HifiasmOverlapRegion>& ol,
    uint64_t max_n_chain, uint32_t chain_cutoff, uint64_t ocv_w, uint32_t rl,
    const vector<uint32_t>& chainHitIndexFlat,
    const vector<HifiasmKmerHit>& allHits)
{
    if (ol.empty()) return;

    uint64_t lch = 0;

    // Stage 1: max_n_chain per-type cap.
    if (max_n_chain > 0 && ol.size() > max_n_chain) {
        std::sort(ol.begin(), ol.end(), [](const HifiasmOverlapRegion& a, const HifiasmOverlapRegion& b) {
            return a.shared_seed > b.shared_seed;
        });

        int32_t n[4] = {}, s[4] = {};
        for (size_t i = 0; i < ol.size(); ++i) {
            const int w = hifiasm_ha_ov_type(ol[i], rl);
            if (static_cast<uint64_t>(++n[w]) == max_n_chain) s[w] = ol[i].shared_seed;
        }

        if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
            // Stage 2: COV_W coverage-window rescue for type-3 overlaps.
            // Window state: high 32 bits = capacity, low 32 bits = used.
            vector<uint64_t> cc;
            uint64_t cwn = 0;
            if (ocv_w > 0 && uint64_t(n[3]) >= max_n_chain && uint64_t(rl) >= ocv_w) {
                cwn = uint64_t(rl) / ocv_w + (uint64_t(rl) % ocv_w ? 1ULL : 0ULL);
                cc.resize(cwn, 0ULL);
                for (uint64_t i = 0, cws = 0; i < cwn; ++i, cws += ocv_w) {
                    uint64_t cwe = std::min(cws + ocv_w, uint64_t(rl));
                    uint64_t cap = std::min((cwe - cws) * (max_n_chain >> 1), uint64_t(UINT32_MAX));
                    cc[i] = cap << 32;
                }
            }

            auto update_cc = [&](const HifiasmOverlapRegion& r) {
                if (!cwn) return;
                const uint64_t rs = r.x_pos_s, re = uint64_t(r.x_pos_e) + 1;
                for (uint64_t m = rs / ocv_w, cws = m * ocv_w; m < cwn; ++m, cws += ocv_w) {
                    uint64_t cwe = std::min(cws + ocv_w, uint64_t(rl));
                    uint64_t os = std::max(rs, cws), oe = std::min(re, cwe);
                    if (oe <= os) break;
                    uint32_t add = static_cast<uint32_t>(std::min(oe - os, uint64_t(UINT32_MAX)));
                    uint32_t used = static_cast<uint32_t>(cc[m]);
                    if (uint64_t(used) + add < UINT32_MAX) cc[m] += add;
                    else cc[m] = (cc[m] & 0xffffffff00000000ULL) | UINT32_MAX;
                }
            };

            auto should_rescue_type3 = [&](const HifiasmOverlapRegion& r) -> bool {
                if (!cwn) return false;
                uint64_t cw0 = 0, cw1 = 0;
                const uint64_t rs = r.x_pos_s, re = uint64_t(r.x_pos_e) + 1;
                for (uint64_t m = rs / ocv_w, cws = m * ocv_w; m < cwn; ++m, cws += ocv_w) {
                    uint64_t cwe = std::min(cws + ocv_w, uint64_t(rl));
                    uint64_t os = std::max(rs, cws), oe = std::min(re, cwe);
                    if (oe <= os) break;
                    uint64_t span = oe - os;
                    uint64_t cap = cc[m] >> 32, used = static_cast<uint32_t>(cc[m]);
                    if (used + span >= cap) cw1 += span;
                    else                    cw0 += span;
                }
                return double(cw0) >= double(cw0 + cw1) * 0.7;
            };

            size_t k = 0;
            for (size_t i = 0; i < ol.size(); ++i) {
                const int w = hifiasm_ha_ov_type(ol[i], rl);
                bool keep = (ol[i].shared_seed >= s[w]);
                if (!keep && w == 3 && cwn > 0) keep = should_rescue_type3(ol[i]);
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

    // If max_n_chain filter was skipped (ol.size() <= max_n_chain), lch was never
    // set above. In hifiasm, lch comes from the per-target chaining loop. We
    // replicate that by scanning surviving overlaps for weak chains.
    if (!lch && chain_cutoff >= 2) {
        for (size_t i = 0; i < ol.size() && !lch; ++i) {
            if (ol[i].align_length < chain_cutoff) lch = 1;
        }
    }

    // Stage 3: R485 weak-chain suppression.
    std::sort(ol.begin(), ol.end(), [](const HifiasmOverlapRegion& a, const HifiasmOverlapRegion& b) {
        return (uint64_t(a.x_pos_s) << 32 | a.x_pos_e) < (uint64_t(b.x_pos_s) << 32 | b.x_pos_e);
    });

    if (lch) {
        size_t l = 0;
        for (size_t i = 0; i < ol.size(); ++i) {
            if (ol[i].align_length < chain_cutoff) {
                const uint64_t zs = ol[i].x_pos_s, ze = uint64_t(ol[i].x_pos_e) + 1;
                uint64_t ob = static_cast<uint64_t>(double(ze - zs) * HIFIASM_OFL);
                if (ob < 16) ob = 16;
                const int64_t osc = int64_t(ol[i].shared_seed) * HIFIASM_CH_SC;
                const uint64_t ocn = uint64_t(ol[i].align_length) << HIFIASM_CH_OCC;

                size_t k = 0;
                for (; k < ol.size() && ze > ol[k].x_pos_s; ++k) {
                    if (ol[k].align_length < chain_cutoff) continue;
                    if (ol[k].align_length < ocn) continue;
                    if (int64_t(ol[k].shared_seed) < osc) continue;

                    uint64_t rs = ol[k].x_pos_s, re = uint64_t(ol[k].x_pos_e) + 1;
                    uint64_t os = std::max(rs, zs), oe = std::min(re, ze);
                    if (!(oe > os && oe - os >= ob)) continue;

                    // Count strong-chain anchors within the overlap zone.
                    uint64_t kn = 0;
                    const uint64_t off = ol[k].non_homopolymer_errors;
                    const uint64_t n_hit = ol[k].align_length;
                    if (off + n_hit > chainHitIndexFlat.size()) continue;

                    for (uint64_t j = 0; j < n_hit && kn < ocn; ++j) {
                        uint64_t g = chainHitIndexFlat[off + j];
                        if (g >= allHits.size()) continue;
                        uint64_t me = allHits[g].self_offset;
                        uint64_t ms = me - (allHits[g].cnt & 0xffu);
                        if (ms >= os && me <= oe) ++kn;
                    }
                    if (kn >= ocn) break;
                }
                if (k < ol.size() && ze > ol[k].x_pos_s) continue;
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
/// Map base positions on Read B to marker ordinals via a linear tandem scan.
/// Hits are sorted by posB so markers can be walked once. Returns false if
/// any hit position has no corresponding marker.
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
        [&](uint32_t a, uint32_t b) {
            return (hitPosB[a] != hitPosB[b]) ? hitPosB[a] < hitPosB[b] : a < b;
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

/// Frequency-based hit weight (hifiasm anchor.cpp:11 HA_KMER_GOOD_RATIO).
/// Three tiers: rare (weight=2), normal (weight=1), repetitive (pow-scaled).
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
        // Check LUT bounds with explicit type for clarity
        if (w < static_cast<uint32_t>(weightLut.size())) {
            return weightLut[w];
        }
        // Fallback: w beyond LUT range, compute dynamically
        return uint32_t(std::pow(double(w), weightExponent));
    }
    return uint32_t(1);
}

/// Transfer chaining parameters from OverlapCandidatesOptions into the
/// InvertedIndexData used by worker threads. Shared by discovery and PAF paths.
template<class InvertedIndexData>
static inline void configureInvertedIndexDataForChaining(
    InvertedIndexData& data,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
    const uint64_t coverageHet,
    const uint64_t coverageHom,
    const double maxDriftRate)
{
    data.maxDriftRate = maxDriftRate;
    data.coverageHet = coverageHet;
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
    // Effective chaining cutoff: min(coverageHom * 5, maxChainingFreq).
    // Matches hifiasm's min(hom_cov * 5, max_kmer_cnt) logic.
    data.maxChainingFreq = uint32_t(std::min(coverageHom * 5,
        uint64_t(overlapCandidatesOptions.maxChainingFreq)));
    data.chainingMode = overlapCandidatesOptions.chainingMode;
    data.minimap2Bw = overlapCandidatesOptions.minimap2Bw;
    data.minimap2MaxGap = overlapCandidatesOptions.minimap2MaxGap;
    data.minimap2MinChainScore = overlapCandidatesOptions.minimap2MinChainScore;
    data.referenceReadCount = overlapCandidatesOptions.referenceReadCount;
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

// ============================================================================
// PARALLEL OVERLAP DISCOVERY WORKER
// ============================================================================
// The InvertedIndexFinder class parallelizes the overlap discovery process.
// Each thread processes a range of reads ("read batches") and produces
// thread-local candidate pairs and alignments. After all threads finish,
// results are merged into the global output vectors.
//
// ## Threading Model:
//   - Reads are divided into batches (typically ~256 reads each).
//   - Each thread picks batches from a shared atomic counter (dynamic scheduling).
//   - Each thread has its own ThreadScratchpad to avoid false sharing.
//   - Results are accumulated in thread-local vectors, then merged once at the end.
//
// ## Per-Read Processing Pipeline (in threadFunction):
//   1. For each marker k-mer in Read A, query the inverted index for hits.
//   2. Collect hits into flatHits (AoS format), applying frequency weighting.
//   3. Optionally downsample high-frequency markers (hifiasm parity).
//   4. Radix sort flatHits by (partnerReadId, posA).
//   5. For each partner Read B, run DP chaining (same + diff strand).
//   6. Extract best chain → build Alignment → emit candidate.
//   7. Apply post-filters (max_n_chain, COV_W, R485).
// ============================================================================
class InvertedIndexFinder : public MultithreadedObject<InvertedIndexFinder> {
public:
    InvertedIndexFinder(
        const Reads& reads,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MemoryMapped::VectorOfVectors<KmerId, uint64_t>& markerKmerIds,
        const Assembler::AlignmentCandidatesInvertedIndexData& invertedIndexData,
        MemoryMapped::Vector<OrientedReadPair>& candidates,
        MemoryMapped::Vector<Alignment>& precomputedAlignments,
        MemoryMapped::Vector<int32_t>& precomputedSharedSeedScores,
        uint64_t maxChainLimit,
        uint64_t threadCount
    ) :
        MultithreadedObject(*this),
        reads(reads),
        markers(markers),
        markerKmerIds(markerKmerIds),
        invertedIndexData(invertedIndexData),
        maxChainLimit(maxChainLimit)
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
        threadSharedSeedScores.resize(threadCount);
        for(auto& v : threadSharedSeedScores) v.reserve(perThreadReserve);
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
        size_t sharedSeedWritePos = precomputedSharedSeedScores.size();
        precomputedSharedSeedScores.resize(sharedSeedWritePos + totalCandidatesFound);
        
        for(size_t i = 0; i < threadCount; i++) {
            const auto& v = threadCandidates[i];
            const auto& a = threadAlignments[i];
            const auto& s = threadSharedSeedScores[i];
            if(!v.empty()) {
                std::copy(v.begin(), v.end(), candidates.begin() + candidateWritePos);
                candidateWritePos += v.size();
                std::copy(a.begin(), a.end(), precomputedAlignments.begin() + alignmentWritePos);
                alignmentWritePos += a.size();
                std::copy(s.begin(), s.end(), precomputedSharedSeedScores.begin() + sharedSeedWritePos);
                sharedSeedWritePos += s.size();
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
    vector<vector<int32_t>> threadSharedSeedScores;
    vector<ThreadScratchpad> threadScratchpads;

    void threadFunction(size_t threadId) {
        vector<OrientedReadPair>& localCandidates = threadCandidates[threadId];
        vector<Alignment>& localAlignments = threadAlignments[threadId];
        vector<int32_t>& localSharedSeedScores = threadSharedSeedScores[threadId];
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
        vector<PendingHighFrequencyMarker> highFrequencyStreak;
        highFrequencyStreak.reserve(64);
        vector<size_t> highFrequencyStreakWorkspace;
        highFrequencyStreakWorkspace.reserve(64);

        // Scratch storage for the strict hifiasm lchain_qdp_mcopy_fast port.
        HifiasmChainDataScratch hifiasmChainDpScratch;
        vector<HifiasmKmerHit> hifiasmAllHits;
        vector<HifiasmKmerHit> hifiasmPairHits;
        vector<HifiasmKmerHit> hifiasmPairHitsStrand0;
        vector<HifiasmKmerHit> hifiasmPairHitsStrand1;
        vector<HifiasmOverlapRegion> hifiasmOverlapRegions;
        vector<uint32_t> hifiasmChainHitIndexFlat;
        
        const uint64_t hashMask = invertedIndexData.hashTable.size() - 1;
        const auto* hashTablePtr = invertedIndexData.hashTable.data();
        const uint64_t kmerLen = invertedIndexData.k;
        const double maxDriftRate = invertedIndexData.maxDriftRate;
        const uint64_t coverageHet = invertedIndexData.coverageHet;
        const double weightExponent = invertedIndexData.weightExponent;
        const double lowFreqMultiplier = invertedIndexData.lowFreqMultiplier;
        const double highFreqMultiplier = invertedIndexData.highFreqMultiplier;
        const uint32_t rareKmerWeight = invertedIndexData.rareKmerWeight;
        const uint64_t lowFreqThreshold = std::max<uint64_t>(2ULL, uint64_t(double(coverageHet) * lowFreqMultiplier));
        const uint64_t highFreqThreshold = std::max<uint64_t>(1ULL, uint64_t(double(coverageHet) * highFreqMultiplier));
        // Gate downsampling behind an absolute minimum (3) to avoid treating
        // nearly all k-mers as high-frequency in low-coverage datasets.
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
        const int chainingMode = invertedIndexData.chainingMode;
        const int32_t minimap2Bw = invertedIndexData.minimap2Bw;
        const int32_t minimap2MaxGap = invertedIndexData.minimap2MaxGap;
        const int32_t minimap2MinChainScore = invertedIndexData.minimap2MinChainScore;
        const uint64_t referenceReadCount = invertedIndexData.referenceReadCount;
        const uint32_t maxChainingFreq = invertedIndexData.maxChainingFreq;

        uint64_t startBatch, endBatch;
        while(getNextBatch(startBatch, endBatch)) {
            for(ReadId readIdA = ReadId(startBatch); readIdA != ReadId(endBatch); ++readIdA) {

                // Skip palindromic reads — they would produce
                // spurious self-overlaps on the opposite strand.
                if(reads.getFlags(readIdA).isPalindromic) continue;

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

                auto computeHitWeight = [&](uint32_t count) -> uint32_t {
                    return computeInvertedIndexHitWeight(
                        count, lowFreqThreshold, highFreqThreshold,
                        highFreqWeightUnit, rareKmerWeight,
                        invertedIndexData.weightLut, weightExponent);
                };

                // Discovery path: emit hits for all partners (readId != readIdA).
                auto appendMarkerHits = [&](const PendingHighFrequencyMarker& markerInfo) {
                    const auto* compactOccs = &invertedIndexData.compactOccurrences[markerInfo.startIdx];
                    for (uint32_t j = 0; j < markerInfo.count; ++j) {
                        if (compactOccs[j].readId != readIdA) {
                            const uint32_t posBEncoded = compactOccs[j].position;
                            const uint32_t posB = posBEncoded & 0x7fffffffU;
                            const uint8_t isRcB = uint8_t(posBEncoded >> 31);
                            scratch.flatHits.push_back(
                                {compactOccs[j].readId, markerInfo.posA, posB,
                                 markerInfo.ordinalA, markerInfo.weight,
                                 markerInfo.isRcA, isRcB});
                        }
                    }
                };

                // Flush buffered high-frequency markers through the shared
                // downsampling logic, emitting surviving hits via appendMarkerHits.
                auto doFlushStreak = [&](uint32_t rightBoundaryPos) {
                    flushHighFrequencyStreak(
                        highFrequencyStreak, highFrequencyStreakWorkspace,
                        lastNonHighBoundaryPos, rightBoundaryPos,
                        highFrequencySampleDistance, maxHighFrequencyPerStreak,
                        [&](const PendingHighFrequencyMarker& m, uint32_t span) {
                            // Hifiasm guard: discard markers whose count >= streak span.
                            if (span > 0 && m.count >= span) return;
                            appendMarkerHits(m);
                        });
                };

                // Hit collection: query each marker against the inverted index,
                // apply frequency-based weighting, and buffer high-freq streaks.
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

                    while(!hashTablePtr[slotIdx].empty) {
                        if(hashTablePtr[slotIdx].key == canonicalKId) {
                            startIdx = hashTablePtr[slotIdx].start;
                            count = hashTablePtr[slotIdx].count;
                            found = true;
                            break;
                        }
                        slotIdx = (slotIdx + 1) & hashMask;
                    }

                    if(!found || (maxChainingFreq > 0 && count > maxChainingFreq)) {
                        if(downsampleHighFrequencyMarkers) {
                            doFlushStreak(posA);
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
                        doFlushStreak(posA);
                    }
                    appendMarkerHits({startIdx, count, kmerHashKey, posA, uint32_t(i), hitWeight, isRcA});
                    lastNonHighBoundaryPos = posA;
                }
                if(downsampleHighFrequencyMarkers && !highFrequencyStreak.empty()) {
                    const uint32_t readLenABoundary = uint32_t(std::min<uint64_t>(
                        readLenA, uint64_t(std::numeric_limits<uint32_t>::max())));
                    doFlushStreak(readLenABoundary);
                }

                if(scratch.flatHits.empty()) continue;
                radixSortFlatHitsByPartnerReadIdAndPosA(scratch.flatHits, scratch.flatHitsTmp);

                // Per-partner chaining: build k_mer_hit array, run DP, apply postfilter.
                // Ref: anchor.cpp:1920 lchain_qgen_mcopy_fast
                    HifiasmLchainDpOptions dpOpt =
                        getHifiasmLchainDpOptions(lchainIsAccurate, uint32_t(kmerLen));

                    // In minimap2-sr mode, override gap/skip penalties to match
                    // minimap2's formula: chn_pen_gap = 0.01 * avg_qspan.
                    // chn_pen_skip is unused in minimap2_comput_sc (skip penalty
                    // is handled by the n_skip counter in the outer DP loop).
                    if(chainingMode == 1) {
                        dpOpt.chnPenGap = 0.01 * static_cast<double>(kmerLen);
                        dpOpt.chnPenSkip = 0.01 * static_cast<double>(kmerLen);
                        dpOpt.quickCheck = false;  // Not applicable in minimap2 mode.
                    }

                    const uint8_t span = uint8_t(std::min<uint64_t>(kmerLen, 255ULL));

                    const uint32_t readLenA32 = uint32_t(std::min<uint64_t>(
                        readLenA, uint64_t(std::numeric_limits<uint32_t>::max())));

                    hifiasmChainDpScratch.resize(0);
                    hifiasmAllHits.clear();
                    hifiasmPairHits.clear();
                    hifiasmOverlapRegions.clear();
                    hifiasmChainHitIndexFlat.clear();

                    // Optional debug: set DINARA_OVERLAP_DEBUG_READ to trace
                    // candidate discovery for a specific read.
                    // Parsed once per thread on first use.
                    static thread_local ReadId overlapDebugReadA = []() -> ReadId {
                        const char* s = std::getenv("DINARA_OVERLAP_DEBUG_READ");
                        if (s && *s) {
                            char* end = nullptr;
                            const unsigned long v = std::strtoul(s, &end, 10);
                            if (end && end != s && *end == 0
                                && v <= std::numeric_limits<uint32_t>::max()) {
                                return ReadId(uint32_t(v));
                            }
                        }
                        return invalidReadId;
                    }();
                    const bool isOverlapDebugRead =
                        (overlapDebugReadA != invalidReadId && readIdA == overlapDebugReadA);

                    static std::mutex overlapDebugMutex;
                    if (isOverlapDebugRead) {
                        std::lock_guard<std::mutex> lock(overlapDebugMutex);
                        cout << timestamp << "[OVERLAP-DBG] read=" << readIdA
                             << " flatHits=" << scratch.flatHits.size() << endl;
                    }

                    // Iterate per partner readIdB (flatHits are sorted by (readIdB, posA)).
                    size_t hitIter = 0;
                    while(hitIter < scratch.flatHits.size()) {
                        const ReadId readIdB = scratch.flatHits[hitIter].partnerReadId;

                        // Skip mirrored pairs, self-comparisons, and palindromic partners.
                        if(readIdB <= readIdA || reads.getFlags(readIdB).isPalindromic) {
                            while(hitIter < scratch.flatHits.size() && scratch.flatHits[hitIter].partnerReadId == readIdB) {
                                ++hitIter;
                            }
                            continue;
                        }

                        // In reference mode, skip read-vs-read pairs (neither is a reference).
                        if(referenceReadCount > 0
                           && readIdA >= ReadId(referenceReadCount)
                           && readIdB >= ReadId(referenceReadCount)) {
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

                        // Per-strand minimum hit count (hifiasm ecovlp.cpp:3274 chain_cutoff=2).
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
                            if (isOverlapDebugRead) {
                                std::lock_guard<std::mutex> lock(overlapDebugMutex);
                                cout << timestamp << "[OVERLAP-DBG]   partner=" << readIdB
                                     << " hits=" << numHits
                                     << " fwd=" << revCount[0]
                                     << " rev=" << revCount[1]
                                     << " -> REJECTED (chain_cutoff)" << endl;
                            }
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
                        hifiasmPairHitsStrand0.clear();
                        hifiasmPairHitsStrand1.clear();
                        hifiasmPairHitsStrand0.reserve(revCount[0]);
                        hifiasmPairHitsStrand1.reserve(revCount[1]);

                        for(size_t k = 0; k < numHits; ++k) {
                            const auto& h = scratch.flatHits[startInFlat + k];
                            const uint8_t rev = uint8_t(h.isRcA ^ h.isRcB); // z->rev ^ y->rev
                            if((rev == 0 && !keepRev0) || (rev == 1 && !keepRev1)) {
                                continue;
                            }

                            // Convert start positions to end positions (hifiasm convention).
                            const uint32_t seedSpan = uint32_t(span);
                            const uint32_t selfOff = h.posA + (seedSpan - 1U);
                            const uint32_t offSame = h.posB + (seedSpan - 1U);
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
                            if(rev == 0) {
                                hifiasmPairHitsStrand0.push_back(kh);
                            } else {
                                hifiasmPairHitsStrand1.push_back(kh);
                            }
                        }

                        // `self_offset` is already non-decreasing within each strand bucket (flatHits are posA-sorted).
                        // Match hifiasm's `(strand, self_offset, offset)` order by sorting only ties on `offset`.
                        sortHifiasmHitsBySelfOffsetThenOffsetRuns(hifiasmPairHitsStrand0);
                        sortHifiasmHitsBySelfOffsetThenOffsetRuns(hifiasmPairHitsStrand1);
                        hifiasmPairHits.insert(
                            hifiasmPairHits.end(),
                            hifiasmPairHitsStrand0.begin(),
                            hifiasmPairHitsStrand0.end());
                        hifiasmPairHits.insert(
                            hifiasmPairHits.end(),
                            hifiasmPairHitsStrand1.begin(),
                            hifiasmPairHitsStrand1.end());

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
                            int64_t(mcopyKhitCutoff),
                            chainingMode,
                            minimap2Bw,
                            minimap2MaxGap,
                            minimap2MinChainScore);
                    }

                    if (isOverlapDebugRead) {
                        std::lock_guard<std::mutex> lock(overlapDebugMutex);
                        cout << timestamp << "[OVERLAP-DBG] pre-postfilter chains="
                             << hifiasmOverlapRegions.size() << endl;
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

                    if (isOverlapDebugRead) {
                        std::lock_guard<std::mutex> lock(overlapDebugMutex);
                        cout << timestamp << "[OVERLAP-DBG] post-postfilter chains="
                             << hifiasmOverlapRegions.size() << endl;
                        for (const auto& r : hifiasmOverlapRegions) {
                            cout << timestamp << "[OVERLAP-DBG]   chain partner="
                                 << r.y_id
                                 << " strand=" << int(r.y_pos_strand)
                                 << " qSpan=" << (r.x_pos_e - r.x_pos_s)
                                 << " anchors=" << r.align_length << endl;
                        }
                    }

                    // Convert surviving overlaps to Dinara candidates/alignments.
                    for(const auto& r : hifiasmOverlapRegions) {
                        const ReadId readIdB = ReadId(r.y_id);
                        const uint64_t readLenB = reads.getReadRawSequenceLength(readIdB);
                        const OrientedReadId orientedReadB(readIdB, 0);
                        const auto& mB = markers[orientedReadB.getValue()];
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

                        // CRITICAL: Check bounds before accessing chain hit indices
                        if(off + n_hit > hifiasmChainHitIndexFlat.size()) {
                            continue;  // Skip this overlap - hit indices out of bounds
                        }

                        // Compute raw (pre-extension) marker chain span from the
                        // first and last chain hits. The extended coordinates
                        // (x_pos_s/x_pos_e) are left-normalized and right-extended
                        // to read boundaries, so they overestimate the actual
                        // marker-supported overlap region.
                        uint32_t rawQmin = UINT32_MAX, rawQmax = 0;
                        uint32_t rawTmin = UINT32_MAX, rawTmax = 0;

                        al.ordinals.reserve(n_hit);
                        for(uint64_t j = 0; j < n_hit; ++j) {
                            const uint64_t g = uint64_t(hifiasmChainHitIndexFlat[size_t(off + j)]);
                            if(g >= hifiasmAllHits.size()) {
                                continue;
                            }
                            const auto& h = hifiasmAllHits[size_t(g)];
                            uint32_t ordB = h.ordinalB;

                            // CRITICAL: Validate ordinal before RC transformation
                            // ordinalB may be UINT32_MAX if mapping failed
                            if(ordB == std::numeric_limits<uint32_t>::max()) {
                                continue;  // Skip this hit - ordinal not mapped
                            }

                            // Track raw base positions for chain span computation.
                            rawQmin = std::min(rawQmin, h.self_offset);
                            rawQmax = std::max(rawQmax, h.self_offset);
                            rawTmin = std::min(rawTmin, h.offset);
                            rawTmax = std::max(rawTmax, h.offset);

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

                        // Reject candidates whose raw marker chain span is below
                        // the threshold. Uses pre-extension base positions from
                        // the actual chain hits, not the extended coordinates.
                        if(invertedIndexData.minOverlapLength > 0) {
                            const uint32_t rawQspan = (rawQmax >= rawQmin) ? (rawQmax - rawQmin + 1) : 0;
                            const uint32_t rawTspan = (rawTmax >= rawTmin) ? (rawTmax - rawTmin + 1) : 0;
                            if(std::min(rawQspan, rawTspan) < invertedIndexData.minOverlapLength) {
                                continue;
                            }
                        }

                        // Optional: reject overlaps with large softclips.
                        // Two complementary checks:
                        // 1. Per-read: reject if a read is "internal" (both its clips exceed maxEndFuzz).
                        // 2. Joint: reject if neither read reaches a given side (min of both clips exceeds maxEndFuzz).
                        if(invertedIndexData.maxEndFuzz > 0) {
                            const uint32_t qRightClip = uint32_t(std::max<int64_t>(int64_t(readLenA) - int64_t(qE), 0));
                            const uint32_t tRightClip = uint32_t(std::max<int64_t>(int64_t(readLenB) - int64_t(tE), 0));
                            // if(qS > invertedIndexData.maxEndFuzz && qRightClip > invertedIndexData.maxEndFuzz) {
                            //     continue;
                            // }
                            // if(tS > invertedIndexData.maxEndFuzz && tRightClip > invertedIndexData.maxEndFuzz) {
                            //     continue;
                            // }
                            const uint32_t leftNeed = std::min(qS, tS);
                            const uint32_t rightNeed = std::min(qRightClip, tRightClip);
                            if(leftNeed > invertedIndexData.maxEndFuzz || rightNeed > invertedIndexData.maxEndFuzz) {
                                continue;
                            }
                        }

                        // readIdA < readIdB is guaranteed by the skip guard above,
                        // so no canonicalization is needed.
                        const uint8_t overlapType =
                            uint8_t(getOverlapType(qS, qE, uint32_t(readLenA)));

                        emittedForRead.push_back(EmittedChainedCandidate{
                            OrientedReadPair(readIdA, readIdB, isSameStrand),
                            std::move(al),
                            r.shared_seed,
                            overlapType,
                            readIdB,
                            0});
                        }

                // Postfilter (max_n_chain + ocv_w rescue + r485 suppression) was
                // already applied in hifiasm_lchain_qgen_mcopy_fast_postfilter.
                for(auto& e : emittedForRead) {
                    localCandidates.push_back(e.candidate);
                    localAlignments.push_back(std::move(e.alignment));
                    localSharedSeedScores.push_back(e.score);
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
// Uses count-then-scatter (hifiasm-style) to avoid the 24-byte intermediate
// array and radix sort. See InvertedIndexBuilder.hpp for implementation.
// =============================================================================
void Assembler::buildInvertedIndex(uint64_t threadCount) {
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    performanceLog << timestamp << "Building Inverted Index..." << endl;

    checkMarkersAreOpen();
    if(!markerKmerIds->isOpen()) {
        throw runtime_error("Marker KmerIds not available for Inverted Index.");
    }

    inverted_index_builder::build(
        invertedIndexData,
        *markers,
        *markerKmerIds,
        assemblerInfo->k,
        threadCount);
}

// =============================================================================
// Phase 5: Run DP chaining on the built index to find alignment candidates.
// =============================================================================
void Assembler::chainAlignmentCandidates(
    double maxDriftRate,
    uint64_t maxChainLimit,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
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
    alignmentCandidatesAlignmentsData.sharedSeedScores.createNew(largeDataName("AlignmentCandidatesInvertedIndexSharedSeed"), largeDataPageSize);
    
    // Safety reserve for performance
    const ReadId readCount = ReadId(markers->size() / 2);
    alignmentCandidates.candidates.reserve(size_t(readCount) * 50);
    alignmentCandidatesAlignmentsData.alignments.reserve(size_t(readCount) * 50); 
    alignmentCandidatesAlignmentsData.sharedSeedScores.reserve(size_t(readCount) * 50);

    // Keep all tuning parameters synchronized with the PAF chaining path.
    configureInvertedIndexDataForChaining(
        invertedIndexData,
        overlapCandidatesOptions,
        assemblerInfo->kmerDistributionInfo.coverageHet,
        assemblerInfo->kmerDistributionInfo.coverageHom,
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
        alignmentCandidatesAlignmentsData.sharedSeedScores,
        maxChainLimit,
        threadCount
    );

    alignmentCandidates.candidates.unreserve();
    alignmentCandidatesAlignmentsData.alignments.unreserve();
    alignmentCandidatesAlignmentsData.sharedSeedScores.unreserve();
    
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
    uint64_t threadCount
) {
    buildInvertedIndex(threadCount);
    chainAlignmentCandidates(maxDriftRate, maxChainLimit, overlapCandidatesOptions, threadCount);
}

// =============================================================================
// PAF CANDIDATE CHAINING PATH
// =============================================================================
// This is the second major entry point into the chaining pipeline. Unlike the
// discovery path (which builds the index and discovers ALL candidate overlaps),
// the PAF path receives pre-determined read pairs from an external PAF file
// and re-chains them using the same DP scoring.
//
// ## Why a Separate Path?
//
//   When an external aligner (e.g., minimap2) has already identified candidate
//   pairs, we skip the discovery phase entirely. We still need to chain them
//   with Dinara's hifiasm-compatible DP to produce the same scoring and overlap
//   coordinates that the assembler expects downstream.
//
// ## Key Differences from Discovery:
//
//   1. **Fixed orientation**: PAF records specify same-strand or diff-strand.
//      The DP only runs for the requested orientation (not both).
//   2. **Single-pair focus**: Each PAF pair is processed independently.
//      There is no radix-sort grouping by partner — we already know the pair.
//   3. **Post-filtering**: The same max_n_chain + COV_W + R485 filters are
//      applied, plus a deduplication pass to keep only the best overlap per
//      (query, partner) pair.
//
// ## Pipeline Stages:
//
//   ┌─────────────────────────────────────────────────────┐
//   │  1. Snapshot imported PAF pairs                      │
//   │  2. For each pair, collect shared k-mer hits         │
//   │  3. Run DP chaining (single orientation)             │
//   │  4. Apply mcopy + weak-chain suppression             │
//   │  5. Deduplicate (best overlap per partner)           │
//   │  6. Apply max_n_chain + COV_W per-read filtering     │
//   │  7. Emit survivors as Alignment candidates           │
//   └─────────────────────────────────────────────────────┘
//
// =============================================================================
void Assembler::chainPafCandidates(
    double maxDriftRate,
    uint64_t maxChainLimit,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
    uint64_t threadCount
) {
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    const uint32_t minChainedMarkerCount =
        uint32_t(std::max(0, overlapCandidatesOptions.minChainMarkerCount));

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
        assemblerInfo->kmerDistributionInfo.coverageHet,
        assemblerInfo->kmerDistributionInfo.coverageHom,
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
    alignmentCandidatesAlignmentsData.sharedSeedScores.createNew(largeDataName("AlignmentCandidatesInvertedIndexSharedSeed"), largeDataPageSize);
    alignmentCandidatesAlignmentsData.sharedSeedScores.reserve(alignmentCandidates.candidates.size());

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
    const uint64_t coverageHet = invertedIndexData.coverageHet;
    const uint64_t lowFreqThreshold = std::max<uint64_t>(
        2ULL, uint64_t(double(coverageHet) * invertedIndexData.lowFreqMultiplier));
    const uint64_t highFreqThreshold = std::max<uint64_t>(
        1ULL, uint64_t(double(coverageHet) * invertedIndexData.highFreqMultiplier));
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
    const uint32_t maxChainingFreq = invertedIndexData.maxChainingFreq;

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

                    // Skip pairs involving palindromic reads.
                    if(reads->getFlags(readIdA).isPalindromic ||
                       reads->getFlags(readIdB).isPalindromic) continue;

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

                    if(numMarkersA == 0 || numMarkersB == 0) continue;

                    const uint64_t readLenA = reads->getReadRawSequenceLength(readIdA);
                    const uint64_t readLenB = reads->getReadRawSequenceLength(readIdB);

                    scratch.clear();
                    scratch.flatHits.reserve(numMarkersA);
                    highFrequencyStreak.clear();
                    int64_t lastNonHighBoundaryPos = -1;

                    auto computeHitWeight = [&](uint32_t count) -> uint32_t {
                        return computeInvertedIndexHitWeight(
                            count, lowFreqThreshold, highFreqThreshold,
                            highFreqWeightUnit, invertedIndexData.rareKmerWeight,
                            invertedIndexData.weightLut, invertedIndexData.weightExponent);
                    };

                    // PAF path: emit only hits matching the specific partner readIdB.
                    auto appendMarkerHits = [&](const PendingHighFrequencyMarker& markerInfo) {
                        const auto* compactOccs = &invertedIndexData.compactOccurrences[markerInfo.startIdx];
                        for (uint32_t j = 0; j < markerInfo.count; ++j) {
                            if (compactOccs[j].readId == readIdB) {
                                const uint32_t posBEncoded = compactOccs[j].position;
                                const uint32_t posB = posBEncoded & 0x7fffffffU;
                                const uint8_t isRcB = uint8_t(posBEncoded >> 31);
                                scratch.flatHits.push_back(
                                    {readIdB, markerInfo.posA, posB, markerInfo.ordinalA,
                                     markerInfo.weight, markerInfo.isRcA, isRcB});
                            }
                        }
                    };

                    auto doFlushStreak = [&](uint32_t rightBoundaryPos) {
                        flushHighFrequencyStreak(
                            highFrequencyStreak, highFrequencyStreakWorkspace,
                            lastNonHighBoundaryPos, rightBoundaryPos,
                            highFrequencySampleDistance, maxHighFrequencyPerStreak,
                            [&](const PendingHighFrequencyMarker& m, uint32_t span) {
                                if (span > 0 && m.count >= span) return;
                                appendMarkerHits(m);
                            });
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

                        while(!hashTablePtr[slotIdx].empty) {
                            if(hashTablePtr[slotIdx].key == canonicalKId) {
                                startIdx = hashTablePtr[slotIdx].start;
                                count = hashTablePtr[slotIdx].count;
                                found = true;
                                break;
                            }
                            slotIdx = (slotIdx + 1) & hashMask;
                        }

                        if(!found || (maxChainingFreq > 0 && count > maxChainingFreq)) {
                            if(downsampleHighFrequencyMarkers) {
                                doFlushStreak(posA);
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
                            doFlushStreak(posA);
                        }
                        appendMarkerHits({startIdx, count, kmerHashKey, posA, uint32_t(i), hitWeight, isRcA});
                        lastNonHighBoundaryPos = posA;
                    }
                    if(downsampleHighFrequencyMarkers && !highFrequencyStreak.empty()) {
                        const uint32_t readLenABoundary = uint32_t(std::min<uint64_t>(
                            readLenA, uint64_t(std::numeric_limits<uint32_t>::max())));
                        doFlushStreak(readLenABoundary);
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

                        // ---------------------------------------------------------
                        // DIFF-STRAND DP (manually inlined, lchain_dp style)
                        // ---------------------------------------------------------
                        // For diff-strand (reverse-complement) overlaps, the target
                        // positions DECREASE as the query positions INCREASE.
                        // Hits are stored in original strand-0 coordinates, so for
                        // diff-strand we look for posB[j] > posB[i] (reverse order).
                        //
                        // The scoring uses the same comput_sc_ch formula, but with
                        // an explicit coordinate flip: posB is converted to forward
                        // coordinates (readLenB - 1 - posB) before computing the
                        // dynamic bandwidth threshold via cal_bw.
                        //
                        // This section is manually inlined (not calling run_main_dp_loop)
                        // because the PAF path runs only one orientation per pair,
                        // and the diff-strand distance checks differ materially from
                        // same-strand (reversed comparisons, flipped bandwidth).
                        // ---------------------------------------------------------
                        if(processDiff && (int32_t)i - st_diff > MAX_ITER) st_diff = (int32_t)i - MAX_ITER;

                        int32_t end_j_diff = st_diff;
                        for(int32_t j = (int32_t)i - 1; processDiff && j >= st_diff; --j) {
                            // For diff-strand: target position must be LARGER than current
                            // (positions are stored in strand-0 coords, so RC pairs have
                            // posB decreasing as posA increases).
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

                        // Flip ordinals for opposite strand
                        if(!useSameStrand) {
                            uint32_t numMB = (uint32_t)markersB.size();
                            for(auto& p : al.ordinals) p[1] = numMB - 1 - p[1];
                        }

                        // Optional: reject candidates with small pre-extension overlap span.
                        if(invertedIndexData.minOverlapLength > 0) {
                            const uint32_t qSpan = al.qe - al.qs;
                            const uint32_t tSpan = al.te - al.ts;
                            if(std::min(qSpan, tSpan) < invertedIndexData.minOverlapLength) {
                                continue;
                            }
                        }

                        // Optional: reject overlaps with large softclips.
                        // Two complementary checks:
                        // 1. Per-read: reject if a read is "internal" (both its clips exceed maxEndFuzz).
                        // 2. Joint: reject if neither read reaches a given side (min of both clips exceeds maxEndFuzz).
                        if(invertedIndexData.maxEndFuzz > 0) {
                            const uint32_t qRightClip = uint32_t(std::max<int64_t>(int64_t(readLenA) - int64_t(al.qe), 0));
                            uint32_t tLeftClip, tRightClip;
                            if(useSameStrand) {
                                tLeftClip = al.ts;
                                tRightClip = uint32_t(std::max<int64_t>(int64_t(readLenB) - int64_t(al.te), 0));
                            } else {
                                tLeftClip = uint32_t(readLenB) - al.te;
                                tRightClip = al.ts;
                            }
                            if(al.qs > invertedIndexData.maxEndFuzz && qRightClip > invertedIndexData.maxEndFuzz) {
                                continue;
                            }
                            if(tLeftClip > invertedIndexData.maxEndFuzz && tRightClip > invertedIndexData.maxEndFuzz) {
                                continue;
                            }
                            const uint32_t leftNeed = std::min(al.qs, tLeftClip);
                            const uint32_t rightNeed = std::min(qRightClip, tRightClip);
                            if(leftNeed > invertedIndexData.maxEndFuzz || rightNeed > invertedIndexData.maxEndFuzz) {
                                continue;
                            }
                        }

                        // readIdA < readIdB is guaranteed: PAF pairs come from the
                        // candidate table which was populated with canonical ordering.
                        const uint32_t readLenA32 = uint32_t(std::min<uint64_t>(
                            readLenA,
                            uint64_t(std::numeric_limits<uint32_t>::max())));
                        const uint8_t overlapType = uint8_t(getOverlapType(qPstart, qPend, readLenA32));

                        localResults.push_back(PafChainedCandidate{
                            OrientedReadPair(readIdA, readIdB, useSameStrand),
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

    // ============================================================
    // MERGE THREAD-LOCAL RESULTS
    // ============================================================
    // All worker threads have finished. Merge their independent result
    // vectors into a single vector for read-level post-processing.
    // This is done once after all workers finish to keep worker code
    // entirely write-only to thread-local storage (no synchronization).
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

    // ============================================================
    // EMIT FINAL RESULTS
    // ============================================================
    // Transfer surviving candidates and their alignments into Dinara's
    // memory-mapped output vectors. After this, the temporary in-memory
    // structures are no longer needed and the inverted index transient
    // data can be freed.
    for(auto& r : mergedResults) {
        alignmentCandidates.candidates.push_back(r.candidate);
        alignmentCandidatesAlignmentsData.alignments.push_back(std::move(r.alignment));
        alignmentCandidatesAlignmentsData.sharedSeedScores.push_back(r.score);
    }

    // Release over-reserved capacity now that final size is known.
    alignmentCandidates.candidates.unreserve();
    alignmentCandidatesAlignmentsData.alignments.unreserve();
    alignmentCandidatesAlignmentsData.sharedSeedScores.unreserve();

    // Free high-memory transient buffers (hash table, compact occurrences, etc.)
    // that are no longer needed after chaining is complete.
    clearInvertedIndexTransientData(invertedIndexData);

    const auto endTime = std::chrono::steady_clock::now();
    const double totalSeconds = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime)).count());
    cout << timestamp << "PAF candidate chaining completed in " << totalSeconds << " s." << endl;
    cout << timestamp << "Chained " << alignmentCandidates.candidates.size() << " candidates." << endl;
}


// ============================================================================
// Palindromic read detection via inverted index + DP chaining.
// Uses the same pipeline as overlap discovery but collects self-hits
// (strand 0 vs strand 1) instead of partner hits, then refines with
// ProjectedAlignment (astarpa) for base-level identity.
// Requires buildInvertedIndex to have been called first.
// ============================================================================

void Assembler::flagPalindromicReads(
    double maxDriftRate,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
    double alignedFractionThreshold,
    double maxErrorRate,
    uint64_t threadCount)
{
    performanceLog << timestamp
        << "Finding palindromic reads via inverted index chaining." << endl;

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    DINARA_ASSERT(!invertedIndexData.compactOccurrences.empty());

    configureInvertedIndexDataForChaining(
        invertedIndexData, overlapCandidatesOptions,
        assemblerInfo->kmerDistributionInfo.coverageHet,
        assemblerInfo->kmerDistributionInfo.coverageHom,
        maxDriftRate);
    rebuildWeightLut(invertedIndexData);

    const ReadId readCount = reads->readCount();
    const auto& allMarkers = *markers;
    const uint64_t kmerLen = assemblerInfo->k;

    const int64_t paMatch = 1;
    const int64_t paMismatch = -5;
    const int64_t paGapOpen = -4;
    const int64_t paGapExtend = -2;

    reads->assertReadsAndFlagsOfSameSize();
    for(ReadId readId = 0; readId < readCount; readId++) {
        reads->setPalindromicFlag(readId, false);
    }

    const uint64_t coverageHet = invertedIndexData.coverageHet;
    const uint32_t lowFreqThreshold = uint32_t(
        coverageHet * invertedIndexData.lowFreqMultiplier);
    const uint32_t highFreqThreshold = uint32_t(
        coverageHet * invertedIndexData.highFreqMultiplier);
    const uint32_t highFreqWeightUnit = std::max(1u, uint32_t(coverageHet));
    const uint32_t rareKmerWeight = invertedIndexData.rareKmerWeight;
    const double weightExponent = invertedIndexData.weightExponent;

    const uint64_t hashMask = invertedIndexData.hashTable.size() - 1;
    const auto* hashTablePtr = invertedIndexData.hashTable.data();

    std::atomic<uint64_t> palindromicCount{0};
    const uint64_t batchSize = 100;
    std::atomic<uint64_t> nextBatch{0};

    auto workerFn = [&]() {
        vector<HifiasmKmerHit> selfHits;
        HifiasmChainDataScratch dpScratch;
        vector<HifiasmOverlapRegion> overlapRegions;
        vector<uint32_t> chainHitIndexFlat;

        while(true) {
            const uint64_t batch = nextBatch.fetch_add(batchSize);
            if(batch >= readCount) break;
            const uint64_t batchEnd = std::min(
                batch + batchSize, uint64_t(readCount));

            for(ReadId readId = ReadId(batch);
                readId < ReadId(batchEnd); readId++) {

                const OrientedReadId oid0(readId, 0);
                const uint64_t readLen =
                    reads->getReadRawSequenceLength(readId);
                if(readLen == 0) continue;

                const auto& markersA = allMarkers[oid0.getValue()];
                const uint32_t numMarkersA = uint32_t(markersA.size());
                if(numMarkersA < 2) continue;

                // Get precomputed canonical kmer IDs for strand 0.
                const KmerId* canonicalIdsA = nullptr;
                const uint8_t* canonicalIsRcA = nullptr;
                bool hasCanonical =
                    (size_t(readId) + 1 <
                     invertedIndexData.strand0CanonicalOffsets.size());
                if(hasCanonical) {
                    const uint64_t b =
                        invertedIndexData.strand0CanonicalOffsets[
                            size_t(readId)];
                    const uint64_t e =
                        invertedIndexData.strand0CanonicalOffsets[
                            size_t(readId) + 1];
                    if(e - b == numMarkersA &&
                       e <= invertedIndexData
                            .strand0CanonicalKmerIds.size() &&
                       e <= invertedIndexData
                            .strand0CanonicalIsRc.size()) {
                        canonicalIdsA =
                            invertedIndexData
                                .strand0CanonicalKmerIds.data() + b;
                        canonicalIsRcA =
                            invertedIndexData
                                .strand0CanonicalIsRc.data() + b;
                    } else {
                        hasCanonical = false;
                    }
                }

                // Fall back to markerKmerIds if no canonical data.
                const auto kmerIdsA = hasCanonical
                    ? span<const KmerId>()
                    : (*markerKmerIds)[oid0.getValue()];

                selfHits.clear();

                for(uint32_t ordA = 0; ordA < numMarkersA; ordA++) {
                    KmerId currentKId;
                    uint8_t isRcA;
                    if(hasCanonical) {
                        currentKId = canonicalIdsA[ordA];
                        isRcA = canonicalIsRcA[ordA];
                    } else {
                        currentKId = kmerIdsA[ordA];
                        KmerId rcKId = getRcKmerId(currentKId, kmerLen);
                        if(rcKId < currentKId) {
                            currentKId = rcKId;
                            isRcA = 1;
                        } else {
                            isRcA = 0;
                        }
                    }

                    // Hash table lookup.
                    uint64_t slotIdx =
                        hashKmer(currentKId) & hashMask;
                    uint64_t startIdx = 0;
                    uint32_t count = 0;
                    bool found = false;
                    while(!hashTablePtr[slotIdx].empty) {
                        if(hashTablePtr[slotIdx].key == currentKId) {
                            startIdx = hashTablePtr[slotIdx].start;
                            count = hashTablePtr[slotIdx].count;
                            found = true;
                            break;
                        }
                        slotIdx = (slotIdx + 1) & hashMask;
                    }
                    if(!found) continue;

                    const uint32_t w = computeInvertedIndexHitWeight(
                        count, lowFreqThreshold, highFreqThreshold,
                        highFreqWeightUnit, rareKmerWeight,
                        invertedIndexData.weightLut, weightExponent);
                    if(w == 0) continue;

                    const uint32_t posA = markersA[ordA].position;
                    const uint32_t seedSpan = uint32_t(
                        std::min<uint64_t>(kmerLen, 255ULL));

                    for(uint64_t idx = startIdx; idx < startIdx + count; idx++) {
                        const auto& occ =
                            invertedIndexData.compactOccurrences[idx];
                        if(occ.readId != readId) continue;

                        const uint32_t posBEncoded = occ.position;
                        const uint32_t posB = posBEncoded & 0x7fffffffU;
                        const uint8_t isRcB =
                            uint8_t(posBEncoded >> 31);

                        // Strand of hit: isRcA ^ isRcB.
                        // For palindrome we want opposite strand (rev=1).
                        const uint8_t rev = isRcA ^ isRcB;
                        if(rev == 0) continue;

                        const uint32_t selfOff = posA + (seedSpan - 1U);
                        const uint32_t offDiff = uint32_t(
                            readLen - 1ULL - uint64_t(posB));

                        HifiasmKmerHit kh{};
                        kh.readID = readId;
                        kh.strand = 1;
                        kh.self_offset = selfOff;
                        kh.offset = offDiff;
                        kh.cnt = (std::min(w, 0xffffffu) << 8)
                                 | seedSpan;
                        kh.ordinalA = ordA;
                        kh.ordinalB =
                            std::numeric_limits<uint32_t>::max();
                        kh.globalIndex = uint32_t(selfHits.size());
                        selfHits.push_back(kh);
                    }
                }

                if(selfHits.size() < 2) continue;

                // Map ordinalB for each hit via binary search.
                // Markers are sorted by position so this is O(n log m)
                // with no heap allocations.
                {
                    bool allMapped = true;
                    for(auto& kh : selfHits) {
                        const uint32_t posB = uint32_t(
                            readLen - 1ULL - uint64_t(kh.offset));
                        // Binary search for marker at posB.
                        uint32_t lo = 0, hi = numMarkersA;
                        bool found = false;
                        while(lo < hi) {
                            const uint32_t mid = lo + (hi - lo) / 2;
                            const uint32_t mpos = markersA[mid].position;
                            if(mpos == posB) {
                                kh.ordinalB = numMarkersA - 1U - mid;
                                found = true;
                                break;
                            } else if(mpos < posB) {
                                lo = mid + 1;
                            } else {
                                hi = mid;
                            }
                        }
                        if(!found) { allMapped = false; break; }
                    }
                    if(!allMapped) continue;
                }

                sortHifiasmHitsBySelfOffsetThenOffsetRuns(selfHits);

                // Reassign globalIndex after sort so chain indices
                // map back to the correct elements.
                for(uint32_t i = 0; i < uint32_t(selfHits.size()); i++) {
                    selfHits[i].globalIndex = i;
                }

                overlapRegions.clear();
                chainHitIndexFlat.clear();

                hifiasm_lchain_qdp_mcopy_fast(
                    selfHits, dpScratch,
                    overlapRegions, chainHitIndexFlat,
                    25, 5000, 500,
                    1.0, 0.05,
                    maxDriftRate,
                    uint32_t(readId),
                    int64_t(readLen), int64_t(readLen),
                    1, 1, 0.0, 0);

                if(overlapRegions.empty()) continue;

                // Pick the chain with the largest raw anchor span
                // (self_offset of last hit minus first hit).
                // We cannot use x_pos_s/x_pos_e because push_ovlp_chain_qgen
                // left-normalizes and right-extends them to sequence boundaries.
                size_t bestIdx = 0;
                uint64_t bestRawSpan = 0;
                for(size_t ri = 0; ri < overlapRegions.size(); ++ri) {
                    const auto& r = overlapRegions[ri];
                    const uint64_t rOff = uint64_t(r.non_homopolymer_errors);
                    const uint64_t rN = uint64_t(r.align_length);
                    if(rN < 2 || rOff + rN > chainHitIndexFlat.size()) continue;
                    const uint32_t firstSelf = selfHits[chainHitIndexFlat[rOff]].self_offset;
                    const uint32_t lastSelf = selfHits[chainHitIndexFlat[rOff + rN - 1]].self_offset;
                    const uint64_t span = (lastSelf >= firstSelf)
                        ? uint64_t(lastSelf - firstSelf + 1)
                        : uint64_t(firstSelf - lastSelf + 1);
                    if(span > bestRawSpan) {
                        bestRawSpan = span;
                        bestIdx = ri;
                    }
                }
                if(bestRawSpan == 0) continue;

                // Chain span gate: the chain must cover at least
                // alignedFractionThreshold of the read length.
                const double chainSpanFrac = double(bestRawSpan) / double(readLen);
                if(chainSpanFrac < alignedFractionThreshold) continue;

                const auto& best = overlapRegions[bestIdx];
                const uint64_t off =
                    uint64_t(best.non_homopolymer_errors);
                const uint64_t nHit =
                    uint64_t(best.align_length);

                Alignment alignment;
                alignment.ordinals.reserve(nHit);
                for(uint64_t j = 0; j < nHit; j++) {
                    const auto& h = selfHits[
                        chainHitIndexFlat[size_t(off + j)]];
                    alignment.ordinals.push_back(
                        {h.ordinalA, h.ordinalB});
                }

                if(alignment.ordinals.size() < 2) continue;

                const array<OrientedReadId, 2> oids = {
                    OrientedReadId(readId, 0),
                    OrientedReadId(readId, 1)
                };

                ProjectedAlignment projectedAlignment(
                    *this, oids, alignment,
                    ProjectedAlignment::Method::QuickRawSparse,
                    paMatch, paMismatch, paGapOpen, paGapExtend);

                const uint64_t totalBases =
                    projectedAlignment.totalLength[0];
                if(totalBases == 0) continue;

                const double errorRate =
                    projectedAlignment.errorRate();

                if(errorRate <= maxErrorRate) {
                    reads->setPalindromicFlag(readId, true);
                    palindromicCount.fetch_add(1);
                }
            }
        }
    };

    vector<std::thread> threads;
    for(uint64_t t = 0; t < threadCount; t++) {
        threads.emplace_back(workerFn);
    }
    for(auto& t : threads) {
        t.join();
    }

    const size_t count = palindromicCount.load();
    assemblerInfo->palindromicReadCount = count;
    cout << "Flagged " << count
         << " reads as palindromic out of "
         << readCount << " total." << endl;
    if(count > 0) {
        cout << "Palindromic reads:";
        for(ReadId readId = ReadId(0); readId < ReadId(readCount); ++readId) {
            if(reads->getFlags(readId).isPalindromic) {
                cout << " " << readId;
            }
        }
        cout << endl;
    }
    cout << "Palindromic fraction is "
         << double(count) / double(readCount) << endl;
}
