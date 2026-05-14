/**
 * @file AssemblerInvertedIndex.cpp
 * @brief High-performance alignment candidate discovery using an Inverted Index.
 *
 * This file implements a Hifiasm-compatible chaining algorithm for finding
 * potential read overlaps. It is Dinara's largest algorithmic module (~4900 lines),
 * serving as the heart of the overlap-detection pipeline. Every scoring formula,
 * tie-break rule, and pruning heuristic is deliberately kept at step-parity with
 * hifiasm v0.25.0-r726 (ec9a8b2) so that the two programs produce identical
 * overlap graphs given equivalent input.
 *
 * ┌─────────────────────────────────────────────────────────────────────┐
 * │                  HIGH-LEVEL DATA-FLOW DIAGRAM                      │
 * │                                                                     │
 * │  Reads (FASTA/Q)                                                    │
 * │       │                                                             │
 * │       ▼                                                             │
 * │  ┌──────────────────────────────────────────────┐                   │
 * │  │ Phase 1: INVERTED INDEX CONSTRUCTION         │                   │
 * │  │  • Extract marker k-mers from every read     │                   │
 * │  │  • Emit (kmerId, readId, position) triples   │                   │
 * │  │  • Build parallel occurrence table            │                   │
 * │  └─────────────────────┬────────────────────────┘                   │
 * │                        ▼                                            │
 * │  ┌──────────────────────────────────────────────┐                   │
 * │  │ Phase 2: LSD RADIX SORT                      │                   │
 * │  │  • O(N) sort on 64-bit packed keys            │                   │
 * │  │  • Groups identical k-mers contiguously       │                   │
 * │  └─────────────────────┬────────────────────────┘                   │
 * │                        ▼                                            │
 * │  ┌──────────────────────────────────────────────┐                   │
 * │  │ Phase 3: HASH TABLE BUILD                    │                   │
 * │  │  • Power-of-2 open-addressing table           │                   │
 * │  │  • O(1) lookup of any k-mer's occurrence list │                   │
 * │  └─────────────────────┬────────────────────────┘                   │
 * │                        ▼                                            │
 * │  ┌──────────────────────────────────────────────┐                   │
 * │  │ Phase 4: COMPACTION & WEIGHT PRE-COMPUTATION │                   │
 * │  │  • Remove self-pairs and singletons           │                   │
 * │  │  • Compute frequency-based hit weights (LUT)  │                   │
 * │  └─────────────────────┬────────────────────────┘                   │
 * │                        ▼                                            │
 * │  ┌──────────────────────────────────────────────┐                   │
 * │  │ Phase 5: PARALLEL DP CHAINING (per read)     │                   │
 * │  │  • Collect shared k-mer hits vs every partner │                   │
 * │  │  • High-freq marker downsampling              │                   │
 * │  │  • Quick-check O(N) prefix on collinear runs  │                   │
 * │  │  • Full O(N·maxIter) DP with gap+skip penalty │                   │
 * │  │  • Backtrack → chain → overlap region         │                   │
 * │  │  • Multi-copy (mcopy) secondary extraction    │                   │
 * │  │  • Post-filter: max_n_chain + COV_W + R485    │                   │
 * │  └─────────────────────┬────────────────────────┘                   │
 * │                        ▼                                            │
 * │  Output: AlignmentCandidates (scored overlap graph)                 │
 * └─────────────────────────────────────────────────────────────────────┘
 *
 * ## Two Entry Points
 *
 * The module supports two distinct usage paths:
 *
 *   1. **Discovery Path** (`findAlignmentCandidatesInvertedIndex`):
 *      Builds the inverted index from scratch and discovers all candidate
 *      overlaps. This is the primary entry point during assembly.
 *
 *   2. **PAF Path** (`chainPafCandidates`):
 *      Receives pre-computed read pairs (from an external PAF file) and
 *      re-chains them using the same DP scoring. Skips index construction
 *      but runs identical chaining + post-filtering.
 *
 * ## Key Design Principles
 *
 *   - **Hifiasm Parity**: Scoring, tie-breaking, and filtering are step-for-step
 *     identical to hifiasm. Variable names in the "strict port" sections
 *     intentionally mirror the original C source for auditability.
 *   - **SoA Layout**: Per-thread scratch uses Structure-of-Arrays for cache
 *     efficiency during the O(N²) DP inner loop.
 *   - **Two Scoring Modes**: `comput_sc_ch` (HiFi / default) and
 *     `comput_sc_ch_ec` (ONT error-correction) differ only in how they
 *     choose between linear and adaptive gap penalties.
 *   - **Multi-Copy Extraction**: Repetitive genomic regions produce multiple
 *     valid chains per read pair; mcopy-fast extracts up to `mcopyNum`
 *     independent chains by greedy node exclusion.
 *   - **Post-Filtering Stack**: Three successive filters mirror hifiasm:
 *     (a) `max_n_chain` per-overlap-type cap, (b) COV_W window saturation
 *      for type-3 overlaps, (c) R485 weak-chain suppression.
 *
 * ## Section Index
 *
 *   Lines ~59-80 .... InvertedIndexTempHit (AoS hit for sorting)
 *   Lines ~82-117 ... radixSortFlatHitsByPartnerReadIdAndPosA (LSD radix sort)
 *   Lines ~119-166 .. getOverlapType (overlap classification for COV_W)
 *   Lines ~167-250 .. HifiasmLchainDpOptions (DP configuration)
 *   Lines ~251-340 .. ThreadScratchpad (per-thread SoA scratch)
 *   Lines ~342-390 .. Hash / reverse-complement utility functions
 *   Lines ~392-530 .. runQuickLinearChainPrefix (O(N) fast path)
 *   Lines ~533-656 .. applyMcopyFastSelection (multi-copy extraction)
 *   Lines ~656-698 .. Hifiasm scoring constants + normal_w
 *   Lines ~700-785 .. HifiasmKmerHit + sortHifiasmHitsBySelfOffsetThenOffsetRuns
 *   Lines ~787-875 .. HifiasmOverlapRegion + HifiasmChainDataScratch
 *   Lines ~875-1040 . hifiasm_cal_bw + hifiasm_comput_sc_ch_ec (scoring)
 *   Lines ~1042-1098  hifiasm_get_chainLen (chain length normalization)
 *   Lines ~1098-1190  hifiasm_push_ovlp_chain_qgen (overlap construction)
 *   Lines ~1192-1310  hifiasm_quick_ck_lchain (strict port)
 *   Lines ~1310-1460  lchain_qdp_mcopy_fast docstring + run_main_dp_loop
 *   Lines ~1461-1700  backtrack_best_chain + emit_best_chain_as_overlap
 *   Lines ~1700-2060  hifiasm_lchain_qdp_mcopy_fast (full pipeline)
 *   Lines ~2060-2290  hifiasm_ha_ov_type + postfilter
 *   Lines ~2290-2470  mapHitPositionsToMarkerOrdinals + weight/config
 *   Lines ~2470-3180  InvertedIndexFinder (parallel worker)
 *   Lines ~3186-3582  buildInvertedIndex + chainAlignmentCandidates
 *   Lines ~3583-4902  chainPafCandidates (PAF chaining path)
 *
 * @note All scoring and tie-breaking rules are strictly Hifiasm-compatible.
 */

#include "Assembler.hpp"
#include "InvertedIndexBuilder.hpp"
#include "hifiasmCoordinateTransforms.hpp"
#include "performanceLog.hpp"
#include "OrientedReadPair.hpp"
#include "timestamp.hpp"
#include "Reads.hpp"
#include <algorithm>
#include <array>

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
    /// Called at the start of each new read pair to prepare for reuse.
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

// ============================================================================
// HIFIASM SCORING HELPER FUNCTIONS
// ============================================================================

/**
 * @brief Normalize score by k-mer occurrence frequency (hifiasm weight normalization).
 *
 * Downweights common k-mers to prevent repetitive regions from dominating
 * the chaining score. If the k-mer occurs x times and has weight y:
 *   - If x >= y: score = x/y (downweight)
 *   - Otherwise: score = 1   (minimum score)
 *
 * @param x  Raw score (e.g., span or count)
 * @param y  K-mer occurrence/frequency weight
 * @return Normalized score (≥1)
 *
 * @complexity O(1)
 * @reference Hifiasm Hash_Table.cpp:20 (normal_w macro)
 */
static inline int32_t hifiasm_normal_w(const int32_t x, const int32_t y) noexcept
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

/**
 * @brief Calculate adaptive bandwidth constraint for anchor chaining.
 *
 * This function computes the maximum allowed coordinate deviation (bandwidth)
 * between two anchors ai and aj. The bandwidth adapts based on the "alignable
 * region" - the portion of both sequences that can participate in the alignment
 * given the current anchor positions.
 *
 * ## Algorithm:
 *   1. Compute the alignable region by normalizing start/end coordinates
 *   2. The alignable span is the min of query and target remaining lengths
 *   3. Bandwidth = alignable_span * bw_rate
 *
 * ## Intuition:
 *   Near sequence ends, we allow more flexibility (larger bandwidth) because
 *   there's more opportunity for valid gaps. In the middle, we're stricter.
 *
 * @param ai      Later anchor (query position i)
 * @param aj      Earlier anchor (query position j < i)
 * @param bw_rate Bandwidth rate (typically 0.3-0.5, i.e., 30-50% deviation)
 * @param sf_l    Query sequence length (self)
 * @param ot_l    Target sequence length (other)
 * @return Maximum allowed bandwidth (coordinate deviation)
 *
 * @complexity O(1)
 * @reference Hifiasm Hash_Table.cpp:1475-1488 (cal_bw)
 */
static inline int32_t hifiasm_cal_bw(
    const HifiasmKmerHit* ai,
    const HifiasmKmerHit* aj,
    const double bw_rate,
    const int64_t sf_l,
    const int64_t ot_l) noexcept
{
    // Extract anchor coordinates (half-open intervals)
    int64_t query_start = static_cast<int64_t>(aj->self_offset);
    int64_t query_end   = static_cast<int64_t>(ai->self_offset) + 1;
    int64_t target_start = static_cast<int64_t>(aj->offset);
    int64_t target_end   = static_cast<int64_t>(ai->offset) + 1;

    // Remaining lengths after current anchors
    const int64_t query_remaining = sf_l - query_end;
    const int64_t target_remaining = ot_l - target_end;

    // Normalize starts: offset by the minimum to align coordinates
    if (query_start <= target_start) {
        query_start = 0;
    } else {
        query_start -= target_start;
    }

    // Normalize ends: extend to full length if one sequence is shorter
    if (query_remaining <= target_remaining) {
        query_end = sf_l;
    } else {
        query_end += target_remaining;
    }

    // Alignable span determines bandwidth
    const int64_t alignable_span = query_end - query_start;
    return static_cast<int32_t>(static_cast<double>(alignable_span) * bw_rate);
}

/**
 * @brief Compute chaining score between two anchors with gap penalties.
 *
 * This is the core scoring function for the DP chaining algorithm. It computes
 * the cost of connecting anchor aj to anchor ai, accounting for:
 *   1. Coordinate differences (query and target gaps)
 *   2. Bandwidth constraints (alignment feasibility)
 *   3. Gap penalties (indels and skipped bases)
 *
 * ## Scoring Formula:
 *   base_score = min(k-mer_span, min(query_gap, target_gap))
 *   normalized_score = base_score / occurrence_weight
 *   penalty = chn_pen_gap * |query_gap - target_gap| + chn_pen_skip * min_gap
 *   final_score = normalized_score - penalty
 *
 * ## Return Values:
 *   - INT32_MIN: Invalid connection (monotonicity violation or bandwidth exceeded)
 *   - ≥0: Valid score (higher is better)
 *
 * @param ai          Later anchor (position i in DP)
 * @param aj          Earlier anchor (position j < i in DP)
 * @param bw_rate     Bandwidth rate for adaptive penalty (typically 0.3-0.5)
 * @param chn_pen_gap Gap penalty coefficient (penalizes indels)
 * @param chn_pen_skip Skip penalty (penalizes large gaps between anchors)
 * @param sl          Query sequence length (self)
 * @param ol          Target sequence length (other)
 * @return Chaining score or INT32_MIN if invalid
 *
 * @complexity O(1) - Single cal_bw call
 * @reference Hifiasm Hash_Table.cpp:1515-1541 (comput_sc_ch_ec)
 */
static inline int32_t hifiasm_comput_sc_ch_ec(
    const HifiasmKmerHit* ai,
    const HifiasmKmerHit* aj,
    const double bw_rate,
    const double chn_pen_gap,
    const double chn_pen_skip,
    const int64_t sl,
    const int64_t ol) noexcept
{
    // -------------------------------------------------------------------------
    // STEP 1: Check monotonicity (query and target must both advance)
    // -------------------------------------------------------------------------
    const int32_t dq = static_cast<int32_t>(
        static_cast<int64_t>(ai->self_offset) - static_cast<int64_t>(aj->self_offset));
    if (dq <= 0) return INT32_MIN;  // Query must advance

    const int32_t dr = static_cast<int32_t>(
        static_cast<int64_t>(ai->offset) - static_cast<int64_t>(aj->offset));
    if (dr <= 0) return INT32_MIN;  // Target must advance

    // -------------------------------------------------------------------------
    // STEP 2: Check bandwidth constraint (gap difference must be reasonable)
    // -------------------------------------------------------------------------
    const int32_t dd = std::abs(dr - dq);  // Gap difference (diagonal deviation)

    // Reject if deviation exceeds bandwidth (with small constant allowance)
    if ((dd > 16) && (dd > hifiasm_cal_bw(ai, aj, bw_rate, sl, ol))) {
        return INT32_MIN;
    }

    // -------------------------------------------------------------------------
    // STEP 3: Compute base score from k-mer span and gap length
    // -------------------------------------------------------------------------
    const int32_t dg = std::min(dr, dq);  // Min gap (effective overlap)
    const int32_t q_span = static_cast<int32_t>(ai->cnt & 0xFFu);  // K-mer span (low 8 bits)

    // Base score is limited by both k-mer span and gap
    int32_t score = std::min(q_span, dg);

    // Normalize by occurrence weight (downweight common k-mers)
    const int32_t occurrence_weight = static_cast<int32_t>(ai->cnt >> 8);
    score = hifiasm_normal_w(score, occurrence_weight);

    // -------------------------------------------------------------------------
    // STEP 4: Apply gap penalties if there's a mismatch or large gap
    // -------------------------------------------------------------------------
    if (dd > 0 || (dg > q_span && dg > 0)) {
        // Linear penalty: proportional to gap difference
        double linear_penalty = chn_pen_gap * static_cast<double>(dd);

        // Adaptive penalty: scales with score and gap ratio
        // CRITICAL: Only compute adaptive penalty if dg > 0 and bw_rate > 0 to avoid division by zero
        if (dg > 0 && bw_rate > 0.0) {
            const double adaptive_penalty =
                static_cast<double>(score) * (static_cast<double>(dd) / static_cast<double>(dg)) / bw_rate;

            // Choose penalty based on gap size (small vs large gaps)
            if (dd < 4) {
                linear_penalty = std::min(linear_penalty, adaptive_penalty);
            } else {
                linear_penalty = std::max(linear_penalty, adaptive_penalty);
            }
        }

        // Add skip penalty (penalizes bases between consecutive anchors)
        linear_penalty += chn_pen_skip * static_cast<double>(dg);

        score -= static_cast<int32_t>(linear_penalty);
    }

    return score;
}

/**
 * @brief Calculate effective chain length after coordinate normalization.
 *
 * This function extends the chain coordinates to sequence boundaries and
 * computes the resulting alignment span. The extension ensures that:
 *   1. Both sequences start at coordinate 0 (left-normalized)
 *   2. Both sequences extend as far as possible to the right
 *
 * ## Algorithm:
 *   1. Left-normalize: Offset both sequences so the smaller start becomes 0
 *   2. Right-extend: Extend both ends until one reaches its sequence boundary
 *   3. Return: Final span of the normalized alignment
 *
 * This is used for tie-breaking when chains have identical scores: prefer
 * the chain with better coverage (longer effective length).
 *
 * @param x_beg Query start position
 * @param x_end Query end position
 * @param xLen  Query sequence length
 * @param y_beg Target start position
 * @param y_end Target end position
 * @param yLen  Target sequence length
 * @return Effective chain length after normalization
 *
 * @complexity O(1)
 * @reference Hifiasm Hash_Table.cpp:779-809 (get_chainLen)
 */
static inline int64_t hifiasm_get_chainLen(
    int64_t x_beg, int64_t x_end, const int64_t xLen,
    int64_t y_beg, int64_t y_end, const int64_t yLen) noexcept
{
    // Left-normalize: shift both sequences so the smaller start becomes 0
    if (x_beg <= y_beg) {
        y_beg -= x_beg;
        x_beg = 0;
    } else {
        x_beg -= y_beg;
        y_beg = 0;
    }

    // Right-extend: extend until one sequence reaches its boundary
    const int64_t x_remaining = xLen - x_end - 1;
    const int64_t y_remaining = yLen - y_end - 1;

    if (x_remaining <= y_remaining) {
        x_end = xLen - 1;
        y_end += x_remaining;
    } else {
        x_end += y_remaining;
        y_end = yLen - 1;
    }

    // Return effective span
    return x_end - x_beg + 1;
}

/**
 * @brief Convert a chained anchor path into an overlap region structure.
 *
 * This function takes the first and last anchors of a DP chain and constructs
 * a normalized overlap region, with coordinates extended to sequence boundaries.
 *
 * ## Coordinate Transformation:
 *   1. **Raw coordinates**: Exact anchor positions (saved for filtering)
 *   2. **Left-normalize**: Shift both sequences so smaller start becomes 0
 *   3. **Right-extend**: Extend to sequence ends (maximizes overlap span)
 *
 * This normalization mirrors hifiasm's overlap representation and ensures
 * consistent coordinate systems for downstream filtering and alignment.
 *
 * @param[out] o    Output overlap region structure
 * @param xid       Query read ID
 * @param xl        Query sequence length
 * @param yl        Target sequence length
 * @param sc        Chain score (shared_seed)
 * @param beg       First anchor in chain (earliest query position)
 * @param end       Last anchor in chain (latest query position)
 *
 * @complexity O(1)
 * @reference Hifiasm Hash_Table.cpp:1752-1780 (push_ovlp_chain_qgen)
 */
static inline void hifiasm_push_ovlp_chain_qgen(
    HifiasmOverlapRegion& o,
    const uint32_t xid,
    const int64_t xl,
    const int64_t yl,
    const int64_t sc,
    const HifiasmKmerHit* beg,
    const HifiasmKmerHit* end) noexcept
{
    // -------------------------------------------------------------------------
    // STEP 1: Initialize basic overlap metadata
    // -------------------------------------------------------------------------
    o.x_id = xid;
    o.y_id = beg->readID;
    o.x_pos_strand = 0;           // Query is always forward
    o.y_pos_strand = beg->strand; // Target strand from chain

    // -------------------------------------------------------------------------
    // STEP 2: Extract raw anchor coordinates (inclusive)
    // -------------------------------------------------------------------------
    o.x_pos_s = beg->self_offset;  // Query start
    o.y_pos_s = beg->offset;       // Target start
    o.x_pos_e = end->self_offset;  // Query end
    o.y_pos_e = end->offset;       // Target end


    // -------------------------------------------------------------------------
    // STEP 3: Left-normalize coordinates (shift to origin)
    // -------------------------------------------------------------------------
    if (o.x_pos_s <= o.y_pos_s) {
        o.y_pos_s -= o.x_pos_s;
        o.x_pos_s = 0;
    } else {
        o.x_pos_s -= o.y_pos_s;
        o.y_pos_s = 0;
    }

    // -------------------------------------------------------------------------
    // STEP 4: Right-extend to sequence boundaries
    // -------------------------------------------------------------------------
    const int64_t x_remaining = xl - static_cast<int64_t>(o.x_pos_e) - 1;
    const int64_t y_remaining = yl - static_cast<int64_t>(o.y_pos_e) - 1;

    // CRITICAL: Only extend if remaining distances are positive (avoid underflow)
    if (x_remaining > 0 && y_remaining > 0) {
        if (x_remaining <= y_remaining) {
            o.x_pos_e = static_cast<uint32_t>(xl - 1);
            o.y_pos_e += static_cast<uint32_t>(x_remaining);
        } else {
            o.y_pos_e = static_cast<uint32_t>(yl - 1);
            o.x_pos_e += static_cast<uint32_t>(y_remaining);
        }
    } else if (x_remaining > 0) {
        // Only x can be extended
        o.x_pos_e = static_cast<uint32_t>(xl - 1);
    } else if (y_remaining > 0) {
        // Only y can be extended
        o.y_pos_e = static_cast<uint32_t>(yl - 1);
    }
    // else: both already at or beyond boundaries, no extension needed

    // -------------------------------------------------------------------------
    // STEP 5: Set scoring and metadata
    // -------------------------------------------------------------------------
    o.shared_seed = static_cast<int32_t>(sc);  // Chain DP score
    o.align_length = 0;                         // Set by caller (num anchors)
    o.non_homopolymer_errors = 0;               // Set by caller (chain hit index)
}

// ============================================================================
// HIFIASM QUICK CHECK LINEAR CHAIN (strict port)
// ============================================================================
// Reference: Hifiasm Hash_Table.cpp:2007-2095 (quick_ck_lchain)
//
// This function performs a fast O(N) pre-scan of the sorted anchor array,
// looking for "runs" of consecutive same-strand hits that are perfectly
// collinear (both query and target positions monotonically increase).
//
// ## How It Works:
//
//   The anchor array `a[0..a_n-1]` is sorted by (self_offset, offset) within
//   each strand group. This function scans left-to-right and identifies maximal
//   contiguous runs where:
//     - All hits share the same strand (a[k].strand == a[l].strand)
//     - Both coordinates are strictly increasing: a[k].self_offset > a[k-1].self_offset
//       AND a[k].offset > a[k-1].offset
//
//   For each such monotonic run, it chains the hits greedily (each hit chains
//   to its immediate predecessor) using the same comput_sc_ch scoring as the
//   full DP. If the entire run qualifies, it updates the best known score/index.
//
// ## Outputs (passed by pointer, hifiasm-style interface):
//
//   - *msc:   Best chain score found by quick-check
//   - *msc_i: Index of the terminal anchor of the best chain
//   - *movl:  Normalized chain length (for tie-breaking)
//   - *plus:  Minimum accumulated score (used as a "quality floor")
//   - *si:    Start index for full DP (everything before was solved)
//   - *ei:    End index for full DP (everything after was solved)
//
//   If quick-check solves the entire array, the full DP range [si, ei) may
//   be empty or already complete, allowing the caller to skip expensive DP.
//
// ## Bandwidth Validation:
//
//   After building a candidate chain, the function checks whether its total
//   diagonal drift (sum of |dr-dq| across all consecutive pairs) exceeds the
//   bandwidth threshold. If so, the chain is invalidated to prevent chaining
//   across structurally unrelated regions.
//
// ## Variable Naming (hifiasm parity):
//
//   - `a`: anchor array (HifiasmKmerHit*)
//   - `f[]`: DP score array (int32_t)
//   - `p[]`: predecessor array (int64_t)
//   - `t[]`: workspace / visit marks (int64_t)
//   - `ii[]`: workspace / chain flags (int32_t)
//   - `l`, `k`: segment boundary and scanner
//   - `is_srt`: flag indicating the current run is sorted
//   - `z`: inner loop variable scanning within a strand-homogeneous segment
// ============================================================================
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
    // Strict port of Hash_Table.cpp:2007-2095
    if (a_n <= 0) return;

    // --- Local variables (hifiasm naming preserved) ---
    int64_t l, k, is_srt = 1, z;  // l=segment start, k=scanner, is_srt=monotonic flag
    HifiasmKmerHit *ai, *aj;       // Pointer pair for scoring
    int64_t dq, dr, dd, dg, q_span, sc, csc, ddt;  // Scoring temporaries
    int64_t plus0, msc0, msc_i0, movl0;  // Per-segment best-chain tracking
    double lin_pen, a_pen;  // Gap penalty components

    // Initialize output accumulators to sentinel values.
    // *msc = best score seen so far (start with worst possible).
    // *movl = best chain length (start with largest to prefer shorter).
    // *si/*ei = DP range that still needs full DP [si, ei).
    *plus = 0;
    *msc = *msc_i = INT32_MIN;
    *movl = INT32_MAX;
    *si = 0;       // Start of unsolved range (grows as quick-check solves from left)
    *ei = a_n;     // End of unsolved range (shrinks as quick-check solves from right)

    // Outer loop: scan the anchor array and identify strand-homogeneous segments.
    // A "segment" is a maximal contiguous range [l, k) where all anchors share
    // the same strand. At segment boundaries (strand change or end-of-array),
    // we attempt to score the segment as a single greedy chain.
    //
    // CRITICAL: Loop goes to k == a_n so we flush the final segment.
    // When k == a_n, we do NOT access a[k] — the boundary check is k < a_n.
    for (k = 1, l = 0; k <= a_n; k++) {
        // Detect segment boundary: strand changes or end of array.
        if (k == a_n || (k < a_n && a[k].strand != a[l].strand)) {
            // Clear workspace entries for end of segment.
            t[k - 1] = 0;
            ii[k - 1] = 0;

            // Only attempt greedy chaining if the segment is monotonically sorted.
            // `is_srt` is cleared to 0 if any inversion was detected (see bottom of loop).
            if (is_srt) {
                // Per-segment accumulators.
                plus0 = 0;             // Minimum score in this segment's chain.
                msc0 = msc_i0 = INT32_MIN;  // Best score and its index.
                movl0 = INT32_MAX;     // Best chain length (for tie-breaking).
                ddt = 0;  // Accumulated diagonal drift across the chain.

                // Initialize the first anchor in the segment as a standalone chain.
                p[l] = -1;  // No predecessor for the first anchor.
                f[l] = int32_t(a[l].cnt & 0xffu);  // Self-score = k-mer span.
                if (f[l] >= msc0) { msc0 = f[l]; msc_i0 = l; }
                if (f[l] < plus0) plus0 = f[l];

                // Greedy forward scan: chain each anchor to its immediate predecessor.
                // Unlike full DP (which considers all previous anchors), quick-check
                // only links z → z-1, exploiting the fact that the run is sorted.
                for (z = l + 1; z < k; z++) {
                    ai = &a[z];
                    aj = &a[z - 1];

                    // Check strict monotonicity on the query axis.
                    dq = int64_t(ai->self_offset) - int64_t(aj->self_offset);
                    if (dq <= 0) break;  // Query must advance — monotonicity broken.

                    // Check strict monotonicity on the target axis.
                    dr = int64_t(ai->offset) - int64_t(aj->offset);
                    if (dr <= 0) break;  // Target must advance — monotonicity broken.

                    // Diagonal deviation (how far off-diagonal this step is).
                    dd = (dr > dq) ? (dr - dq) : (dq - dr);

                    // Bandwidth check: reject if diagonal drift is too large.
                    if ((dd > 16) && (dd > hifiasm_cal_bw(&a[z], &a[z - 1], bw_rate, xl, yl))) break;

                    // Score computation (identical to comput_sc_ch).
                    dg = (dr < dq) ? dr : dq;  // Effective step size.
                    q_span = int64_t(ai->cnt & 0xffu);  // K-mer span.
                    sc = (q_span < dg) ? q_span : dg;    // Base score = min(span, step).
                    sc = hifiasm_normal_w(int32_t(sc), int32_t(ai->cnt >> 8));  // Weight normalize.

                    // Apply gap penalty if there's any diagonal deviation or a large gap.
                    if (dd || (dg > q_span && dg > 0)) {
                        lin_pen = chn_pen_gap * double(dd);
                        a_pen = double(sc) * (double(dd) / double(dg)) / bw_rate;
                        // Small gaps (dd<4): use smaller penalty. Large gaps: use larger.
                        if (dd < 4) lin_pen = (lin_pen > a_pen) ? a_pen : lin_pen;
                        else lin_pen = (lin_pen < a_pen) ? a_pen : lin_pen;
                        lin_pen += chn_pen_skip * double(dg);  // Skip penalty for gap length.
                        sc -= int32_t(lin_pen);
                    }

                    // Add predecessor's score to get cumulative chain score.
                    sc += f[z - 1];

                    // Bail out if chaining through predecessor is worse than starting fresh.
                    csc = int64_t(a[z].cnt & 0xffu);  // Self-score of current anchor.
                    if (sc < csc) break;  // Better to start a new chain here.

                    // Accept this link: record predecessor and cumulative score.
                    p[z] = z - 1;
                    f[z] = int32_t(sc);
                    ddt += dd;  // Accumulate total diagonal drift.
                    if (f[z] >= msc0) { msc0 = f[z]; msc_i0 = z; }
                    if (f[z] < plus0) plus0 = f[z];
                }

                // If the entire segment was chained AND the best anchor is the last one,
                // this segment is fully solved — update the global outputs.
                if ((z >= k) && (msc_i0 == (k - 1))) {
                    // Validate that total diagonal drift doesn't exceed bandwidth.
                    // This catches pathological cases where each step is small but
                    // the cumulative drift across the whole chain is too large.
                    bool chain_valid = true;
                    if ((k - l >= 2) && (ddt > 16) && (ddt > hifiasm_cal_bw(&a[k - 1], &a[l], bw_rate, xl, yl))) {
                        chain_valid = false;  // Invalidate chain - exceeds bandwidth
                        msc_i0 = INT32_MIN;   // Sentinel value for invalid
                    }

                    // Only update global best if this chain passed validation.
                    if (chain_valid && msc_i0 >= 0) {
                        if (msc0 >= (*msc)) {
                            // Compute normalized chain length for tie-breaking.
                            movl0 = hifiasm_get_chainLen(
                                a[msc_i0].self_offset, a[msc_i0].self_offset, xl,
                                a[msc_i0].offset, a[msc_i0].offset, yl);
                            // Update best: prefer higher score, then shorter overlap.
                            if (msc0 > (*msc) || movl0 < (*movl)) {
                                *msc = msc0;
                                *msc_i = msc_i0;
                                *movl = movl0;
                            }
                        }
                        // Track minimum accumulated score for quality floor.
                        if (plus0 < (*plus)) *plus = plus0;

                        // Narrow the DP range: if this solved segment is at the left
                        // edge, advance *si past it; if at the right edge, retract *ei.
                        if ((*ei) > k) (*si) = k;
                        else (*ei) = l;
                    }
                }
            }
            // Advance segment start to current position for the next segment.
            l = k;
            is_srt = 1;  // Reset monotonicity flag for the new segment.
        } else {
            // Within a same-strand segment: check if monotonicity still holds.
            // If either coordinate fails to increase, mark the segment as unsorted.
            if ((a[k].self_offset <= a[k - 1].self_offset) || (a[k].offset <= a[k - 1].offset)) is_srt = 0;
            t[k - 1] = 0;
            ii[k - 1] = 0;
        }
    }
}

// ============================================================================
// CORE DP CHAINING ALGORITHM WITH MULTI-COPY EXTRACTION
// ============================================================================

/**
 * @brief Fast DP chaining with optional quick-check optimization and multi-copy extraction.
 *
 * This implements hifiasm's lchain_qdp_mcopy_fast algorithm (Hash_Table.cpp:2097-2284),
 * which finds optimal anchor chains using dynamic programming, with optional extraction
 * of secondary chains for repetitive regions.
 *
 * ## Algorithm Overview:
 *
 *   **Phase 1: Optional Quick Check** (if quick_check=true)
 *     - O(N) fast path for perfectly collinear hits
 *     - Identifies promising score range [si, ei] for full DP
 *     - Computes initial best score (msc, msc_i)
 *
 *   **Phase 2: Main DP** (O(N * max_iter) with pruning)
 *     ```
 *     For each anchor i:
 *       f[i] = max_{j < i} (score(i, j) + f[j])
 *       where j is within lookback window and passes bandwidth constraint
 *     ```
 *     - Pruning optimizations:
 *       * max_iter: Limit lookback window
 *       * max_skip: Stop after skipping too many ancestors
 *       * bandwidth: Reject connections with excessive coordinate deviation
 *
 *   **Phase 3: Backtracking**
 *     - Follow predecessor pointers from best anchor (msc_i)
 *     - Mark chain membership in ii[] array
 *
 *   **Phase 4: Multi-Copy Extraction** (optional, if mcopy_num > 1)
 *     - Extract up to mcopy_num high-scoring chains
 *     - Threshold: score >= mcopy_rate * best_score
 *     - Use radix sort for efficient top-k selection
 *     - Emit each chain as separate overlap region
 *
 * ## Complexity:
 *   - Time: O(N * max_iter) where N = number of anchors
 *     * Quick check adds O(N) preprocessing
 *     * Multi-copy adds O(N log N) sorting
 *   - Space: O(N) for DP arrays (reused scratch space)
 *
 * ## Key Optimizations:
 *   1. **Max-skip pruning**: Stop after skipping max_skip predecessors
 *   2. **Lookback window**: Only consider last max_iter anchors
 *   3. **Max_ii tracking**: Cache best score in distance window
 *   4. **Quick check**: O(N) fast path for collinear hits (ONT)
 *
 * ## Multi-Copy Extraction:
 *   For repetitive genomic regions, mcopy extracts multiple non-overlapping
 *   chains to capture alternative alignments. Useful for:
 *   - Segmental duplications
 *   - Tandem repeats
 *   - Paralogous genes
 *
 * @param[in]     a                   K-mer hit array (sorted by self_offset, offset)
 * @param[in,out] dp                  Reusable DP scratch space (resized internally)
 * @param[out]    res                 Output overlap regions (chains)
 * @param[out]    chainHitIndexFlat   Flattened chain hit indices (for Alignment construction)
 * @param max_skip       Max consecutive skipped predecessors before pruning
 * @param max_iter       Max lookback window (limits O(N²) to O(N*max_iter))
 * @param max_dis        Max coordinate distance for max_ii tracking
 * @param chn_pen_gap    Gap penalty coefficient (indel penalty)
 * @param chn_pen_skip   Skip penalty (bases between anchors)
 * @param bw_rate        Bandwidth rate (adaptive penalty, typically 0.3-0.5)
 * @param xid            Query read ID
 * @param xl             Query sequence length
 * @param yl             Target sequence length
 * @param quick_check    Enable O(N) quick check preprocessing (1=enable, 0=disable)
 * @param mcopy_num      Number of chains to extract (1=best only, >1=mcopy)
 * @param mcopy_rate     Score threshold for mcopy (fraction of best score, e.g., 0.2)
 * @param mcopy_khit_cutoff Minimum anchors required for mcopy extraction
 *
 * @complexity O(N * max_iter) + O(N log N) for mcopy, where N = number of anchors
 * @reference Hifiasm Hash_Table.cpp:2097-2284 (lchain_qdp_mcopy_fast)
 * @see hifiasm_quick_ck_lchain for quick check details
 * @see hifiasm_comput_sc_ch_ec for scoring function
 */

// ============================================================================
// EXTRACTED HELPER FUNCTIONS FOR HIFIASM CHAINING
// ============================================================================
// These functions extract distinct phases from the main chaining algorithm
// to improve readability, testability, and maintainability. All functions
// are marked 'inline' to ensure zero performance overhead.

/**
 * @brief Run the main O(N*max_iter) dynamic programming loop to find optimal chains.
 *
 * This implements the core chaining DP with multiple optimizations to achieve
 * practical O(N*max_iter) complexity instead of O(N²).
 *
 * ## DP Recurrence:
 * ```
 * f[i] = max_{j < i, same_strand} (score(i, j) + f[j])
 * ```
 * Where:
 * - f[i]: Best chain score ending at anchor i
 * - score(i, j): Transition score from anchor j to i, accounting for:
 *   1. Anchor span (k-mer length)
 *   2. Gap penalties (linear and skip)
 *   3. Bandwidth constraints
 *   4. Occurrence weight normalization
 * - p[i]: Predecessor pointer (index j that maximized f[i])
 *
 * ## Optimizations:
 *
 * ### 1. Lookback Window (max_iter)
 * - Only consider last max_iter predecessors (typically 50-100)
 * - Prevents O(N²) blowup on long reads
 * - Trade-off: May miss distant optimal predecessors, but rare in practice
 *
 * ### 2. Max-Skip Pruning
 * - Stop searching after max_skip (typically 25) consecutive skipped anchors
 * - "Skip" = anchor j doesn't improve f[i] (score(i,j) + f[j] ≤ f[i])
 * - Implemented via t[] array marking visited predecessors
 *
 * ### 3. Sliding Distance Window (max_ii caching)
 * - Cache best anchor within distance threshold (max_dis, typically 5000bp)
 * - Exploits locality: nearby anchors more likely to be good predecessors
 * - Provides fallback after max-skip pruning
 *
 * ### 4. Strand Filtering
 * - Only connect same-strand anchors (can't chain + and - strand)
 * - Implemented via while loop advancing start pointer
 *
 * @param a Anchor array
 * @param si Start index (from quick check or 0)
 * @param ei End index (from quick check or a_n)
 * @param max_skip Max consecutive skipped predecessors before pruning
 * @param max_iter Max lookback window size
 * @param max_dis Max distance for max_ii caching
 * @param chn_pen_gap Gap penalty coefficient
 * @param chn_pen_skip Skip penalty coefficient
 * @param bw_rate Bandwidth rate
 * @param xl Query length
 * @param yl Target length
 * @param p[out] Predecessor pointers
 * @param f[out] DP scores
 * @param t[inout] Temporary workspace (for skip tracking)
 * @param ii[out] Chain membership flags (cleared)
 * @param msc[inout] Best chain score
 * @param msc_i[inout] Best chain ending index
 * @param movl[inout] Overlap length (tie-breaking)
 * @param plus[inout] Minimum score (for normalization)
 *
 * @complexity O(N * max_iter) amortized
 */
static inline void run_main_dp_loop(
    HifiasmKmerHit* a,
    const int64_t si,
    const int64_t ei,
    const int64_t max_skip,
    const int64_t max_iter,
    const int64_t max_dis,
    const double chn_pen_gap,
    const double chn_pen_skip,
    const double bw_rate,
    const int64_t xl,
    const int64_t yl,
    int64_t* p,
    int32_t* f,
    int64_t* t,
    int32_t* ii,
    int64_t* msc,
    int64_t* msc_i,
    int64_t* movl,
    int64_t* plus) noexcept
{
    int64_t max_f, n_skip, st, max_j, end_j, sc;
    int64_t max_ii = -1;
    int64_t ovl;
    int32_t max, tmp;
    int64_t i, j;

    // Main DP loop: compute best chain ending at each anchor
    for (i = st = si; i < ei; ++i) {
        // Initialize: base score is k-mer span (no predecessors yet)
        max_f = static_cast<int64_t>(a[static_cast<size_t>(i)].cnt & 0xFFu);
        n_skip = 0;
        max_j = end_j = -1;

        // -----------------------------------------------------------------
        // Optimization 1: Limit lookback window to max_iter
        // -----------------------------------------------------------------
        if ((i - st) > max_iter) {
            st = i - max_iter;
        }

        // -----------------------------------------------------------------
        // Optimization 4: Ensure st points to same-strand anchor
        // -----------------------------------------------------------------
        // CRITICAL: Add bounds check to prevent infinite loop
        while (st < i && a[static_cast<size_t>(i)].strand != a[static_cast<size_t>(st)].strand) {
            ++st;
        }

        // If no matching strand found in lookback window, skip DP for this anchor
        if (st >= i || a[static_cast<size_t>(i)].strand != a[static_cast<size_t>(st)].strand) {
            // No valid predecessor with same strand - anchor stands alone
            max_f = static_cast<int64_t>(a[static_cast<size_t>(i)].cnt & 0xFFu);
            p[static_cast<size_t>(i)] = -1;
            f[static_cast<size_t>(i)] = static_cast<int32_t>(max_f);
            continue;
        }

        // -----------------------------------------------------------------
        // Inner DP Loop: Find best predecessor j for anchor i
        // -----------------------------------------------------------------
        for (j = i - 1; j >= st; --j) {
            // Compute transition score from j to i
            sc = hifiasm_comput_sc_ch_ec(
                &a[static_cast<size_t>(i)], &a[static_cast<size_t>(j)],
                bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);

            // INT32_MIN means invalid connection (bandwidth or monotonicity violation)
            if (sc == INT32_MIN) continue;

            // DP recurrence: total score = transition + predecessor score
            sc += f[static_cast<size_t>(j)];

            // Update best score and predecessor
            if (sc > max_f) {
                max_f = sc;
                max_j = j;
                if (n_skip > 0) --n_skip;  // Reset skip counter on improvement
            }
            // -----------------------------------------------------------------
            // Optimization 2: Max-skip pruning
            // -----------------------------------------------------------------
            // Stop if we've skipped too many consecutive ancestors
            else if (t[static_cast<size_t>(j)] == static_cast<int32_t>(i)) {
                if (++n_skip > max_skip) break;
            }

            // Mark that j's predecessor was visited (for max-skip tracking)
            if (p[static_cast<size_t>(j)] >= 0) {
                t[static_cast<size_t>(p[static_cast<size_t>(j)])] = i;
            }
        }
        end_j = j;  // Remember where we stopped

        // -----------------------------------------------------------------
        // Optimization 3: Sliding distance window (max_ii caching)
        // -----------------------------------------------------------------
        // max_ii tracks the anchor with the highest score within max_dis
        // distance. This allows recovery from skipped predecessors (max-skip
        // may have terminated early, but max_ii provides a fallback).

        // Recompute max_ii if:
        //   1. Not yet initialized (max_ii < 0)
        //   2. Out of distance window (> max_dis)
        //   3. Different strand (invalid connection)
        if ((max_ii < 0) ||
            (a[static_cast<size_t>(i)].self_offset >
             a[static_cast<size_t>(max_ii)].self_offset + max_dis) ||
            (a[static_cast<size_t>(i)].strand != a[static_cast<size_t>(max_ii)].strand)) {

            // Scan backwards to find best score within distance window
            max = INT32_MIN;
            max_ii = -1;
            for (j = i - 1;
                 (j >= st) &&
                 (a[static_cast<size_t>(i)].self_offset <=
                  max_dis + a[static_cast<size_t>(j)].self_offset) &&
                 (a[static_cast<size_t>(i)].strand == a[static_cast<size_t>(j)].strand);
                 --j) {
                if (max < f[static_cast<size_t>(j)]) {
                    max = f[static_cast<size_t>(j)];
                    max_ii = j;
                }
            }
        }

        // If max_ii is valid and was skipped by max-skip pruning (j < end_j),
        // try connecting to it as a fallback predecessor
        if ((max_ii >= 0) && (max_ii < end_j) &&
            (a[static_cast<size_t>(i)].strand == a[static_cast<size_t>(max_ii)].strand)) {

            tmp = hifiasm_comput_sc_ch_ec(
                &a[static_cast<size_t>(i)], &a[static_cast<size_t>(max_ii)],
                bw_rate, chn_pen_gap, chn_pen_skip, xl, yl);

            if (tmp != INT32_MIN && max_f < tmp + f[static_cast<size_t>(max_ii)]) {
                max_f = tmp + f[static_cast<size_t>(max_ii)];
                max_j = max_ii;
            }
        }

        // -----------------------------------------------------------------
        // Finalize DP for anchor i
        // -----------------------------------------------------------------
        f[static_cast<size_t>(i)] = static_cast<int32_t>(max_f);
        p[static_cast<size_t>(i)] = max_j;

        // Update max_ii for next iteration if current anchor is better
        if ((max_ii < 0) ||
            ((a[static_cast<size_t>(i)].self_offset <=
              max_dis + a[static_cast<size_t>(max_ii)].self_offset) &&
             (a[static_cast<size_t>(i)].strand == a[static_cast<size_t>(max_ii)].strand) &&
             (f[static_cast<size_t>(max_ii)] < f[static_cast<size_t>(i)]))) {
            max_ii = i;
        }

        // Track best overall chain (msc = max score, msc_i = anchor index)
        // Tie-break by overlap length (prefer longer coverage)
        if (f[static_cast<size_t>(i)] >= (*msc)) {
            ovl = hifiasm_get_chainLen(
                a[static_cast<size_t>(i)].self_offset, a[static_cast<size_t>(i)].self_offset, xl,
                a[static_cast<size_t>(i)].offset, a[static_cast<size_t>(i)].offset, yl);

            if (f[static_cast<size_t>(i)] > (*msc) || ovl < (*movl)) {
                (*msc) = f[static_cast<size_t>(i)];
                (*msc_i) = i;
                (*movl) = ovl;
            }
        }

        // Track minimum score (for mcopy score normalization)
        if (f[static_cast<size_t>(i)] < (*plus)) {
            (*plus) = f[static_cast<size_t>(i)];
        }

        // Clear chain membership flag (used in backtracking)
        ii[static_cast<size_t>(i)] = 0;
    }
}

/**
 * @brief Backtrack from best chain anchor to extract the full chain.
 *
 * After the DP loop, we have:
 * - f[]: Best chain scores ending at each anchor
 * - p[]: Predecessor pointers for each anchor
 * - msc_i: Index of anchor with best chain score
 *
 * This function follows the predecessor pointers from msc_i backwards
 * to reconstruct the optimal chain.
 *
 * ## Algorithm:
 * ```
 * Start at best anchor (msc_i)
 * While current anchor has a predecessor:
 *   1. Mark anchor as part of best chain (ii[i] = 1)
 *   2. Store anchor index in t[] (reverse order)
 *   3. Follow predecessor: i = p[i]
 *   4. Increment chain length: cL++
 * ```
 *
 * ## Output Format:
 * - t[0..cL-1]: Chain anchor indices in **reverse order**
 *   - t[0] = last anchor (highest query position)
 *   - t[cL-1] = first anchor (lowest query position)
 * - ii[i] = 1: Anchor i is in the best chain
 * - ii[i] = 0: Anchor i is not in the best chain
 *
 * @param msc_i Best chain ending index
 * @param p Predecessor pointers
 * @param t[out] Chain indices in reverse order
 * @param ii[out] Chain membership flags
 *
 * @return Chain length cL (number of anchors in the chain)
 *
 * @complexity O(cL) where cL = chain length (typically << N)
 */
static inline int64_t backtrack_best_chain(
    const int64_t msc_i,
    int64_t* p,
    int64_t* t,
    int32_t* ii) noexcept
{
    int64_t cL = 0;
    int64_t i;

    // Follow predecessor pointers backwards
    // Store chain in reverse order: t[0] = last anchor, t[cL-1] = first anchor
    for (i = msc_i; i >= 0; i = p[static_cast<size_t>(i)]) {
        ii[static_cast<size_t>(i)] = 1;  // Mark as part of best chain
        t[static_cast<size_t>(cL++)] = i; // Store index (reverse order)
    }

    return cL;
}

/**
 * @brief Convert the backtracked chain into a HifiasmOverlapRegion and emit it.
 *
 * This function takes the chain anchor indices (from backtracking) and:
 * 1. Computes overlap region coordinates (query and target spans)
 * 2. Normalizes coordinates to sequence boundaries
 * 3. Stores chain hit indices for alignment construction
 * 4. Appends the overlap region to the result vector
 *
 * ## Input Format (from backtracking):
 * - t[0..cL-1]: Chain anchor indices in **reverse order**
 *   - t[cL-1] = first anchor (beg)
 *   - t[0] = last anchor (end)
 * - cL: Chain length (number of anchors)
 *
 * ## Coordinate System:
 * - Query (x): Always forward strand
 * - Target (y): May be forward (strand=0) or reverse-complement (strand=1)
 * - Coordinates are normalized to sequence boundaries (left-extended, right-extended)
 *
 * @param a Anchor array
 * @param t Chain anchor indices (reverse order, size cL)
 * @param cL Chain length
 * @param msc Chain score
 * @param xid Query read ID
 * @param xl Query length
 * @param yl Target length
 * @param res[out] Result vector (overlap regions)
 * @param chainHitIndexFlat[out] Flattened hit indices (for alignment)
 *
 * @complexity O(cL) where cL = chain length
 */
static inline void emit_best_chain_as_overlap(
    HifiasmKmerHit* a,
    int64_t* t,
    const int64_t cL,
    const int64_t msc,
    const uint32_t xid,
    const int64_t xl,
    const int64_t yl,
    vector<HifiasmOverlapRegion>& res,
    vector<uint32_t>& chainHitIndexFlat) noexcept
{
    // CRITICAL: Only emit if we found a valid chain (cL > 0)
    // If no chain was found (msc_i < 0), cL == 0 and t[] is empty
    if (cL > 0) {
        HifiasmOverlapRegion z{};

        // -----------------------------------------------------------------
        // Compute overlap region from first and last anchor in chain
        // -----------------------------------------------------------------
        // Note: t[] is in reverse order, so t[cL-1] is first, t[0] is last
        hifiasm_push_ovlp_chain_qgen(
            z, xid, xl, yl, msc,
            &a[static_cast<size_t>(t[static_cast<size_t>(cL - 1)])],  // First anchor (earliest)
            &a[static_cast<size_t>(t[0])]);                             // Last anchor (latest)

        // -----------------------------------------------------------------
        // Set chain metadata
        // -----------------------------------------------------------------
        z.align_length = static_cast<uint32_t>(cL);  // Number of anchors in chain

        // Store starting index in chainHitIndexFlat (for Alignment construction)
        z.non_homopolymer_errors = static_cast<uint32_t>(chainHitIndexFlat.size());

        // -----------------------------------------------------------------
        // Append chain hit indices (for Alignment construction)
        // -----------------------------------------------------------------
        // Copy anchor indices in increasing self_offset order
        // (t[] is in reverse order, so iterate backwards)
        for (int64_t i = 0; i < cL; ++i) {
            const int64_t local = t[static_cast<size_t>(cL - i - 1)];
            chainHitIndexFlat.push_back(a[static_cast<size_t>(local)].globalIndex);
        }

        res.push_back(z);
    }
}

static inline void hifiasm_lchain_qdp_mcopy_fast(
    vector<HifiasmKmerHit>& a,
    HifiasmChainDataScratch& dp,
    vector<HifiasmOverlapRegion>& res,
    vector<uint32_t>& chainHitIndexFlat,
    const int64_t max_skip,
    const int64_t max_iter,
    const int64_t max_dis,
    const double chn_pen_gap,
    const double chn_pen_skip,
    const double bw_rate,
    const uint32_t xid,
    const int64_t xl,
    const int64_t yl,
    const int64_t quick_check,
    const int64_t mcopy_num,
    const double mcopy_rate,
    const int64_t mcopy_khit_cutoff)
{
    // =========================================================================
    // INITIALIZATION
    // =========================================================================
    const int64_t a_n = static_cast<int64_t>(a.size());
    if (a_n <= 0) return;  // No anchors to chain

    // Resize DP scratch space to accommodate all anchors
    dp.resize(static_cast<size_t>(a_n));

    // Raw array pointers for cache-efficient access
    int64_t* p = dp.pre.data();   // Predecessor pointers (backtracking)
    int64_t* t = dp.tmp.data();   // Visit marks (DP) or heap keys (mcopy)
    int32_t* f = dp.score.data(); // DP scores f[i]
    int32_t* ii = dp.occ.data();  // Chain membership flags

    // DP state variables (many locals moved to helper functions for clarity)
    int64_t msc, msc_i, movl, plus = 0;
    int64_t min_sc, ch_n, si, ei;
    int64_t i, k, cL = 0;
    int64_t sc;  // Still used in multi-copy extraction phase

    // =========================================================================
    // PHASE 1: OPTIONAL QUICK CHECK (O(N) fast path for collinear hits)
    // =========================================================================
    // Quick check attempts O(N) linear chaining for perfectly collinear hits.
    // If successful, it:
    //   1. Computes initial best score (msc, msc_i)
    //   2. Identifies promising DP range [si, ei]
    //   3. Initializes predecessor pointers p[] and scores f[]
    // If quick check fails or is disabled, fall back to full DP.

    if (quick_check) {
        hifiasm_quick_ck_lchain(
            a.data(), a_n, xl, yl, chn_pen_gap, chn_pen_skip, bw_rate,
            p, t, f, ii, &plus, &msc, &msc_i, &movl, &si, &ei);
    } else {
        // No quick check: initialize for full DP over entire range
        msc = msc_i = INT32_MIN;
        movl = INT32_MAX;
        plus = 0;
        si = 0;
        ei = a_n;
        std::fill(dp.tmp.begin(), dp.tmp.end(), int64_t(0));
    }

    // Safety: Validate quick check results. If invalid (e.g., due to coordinate
    // convention mismatch), fall back to full DP to ensure we always have a valid chain.
    if (msc_i < 0 || si < 0 || ei < 0 || si > a_n || ei > a_n) {
        msc = msc_i = INT32_MIN;
        movl = INT32_MAX;
        plus = 0;
        si = 0;
        ei = a_n;
        std::fill(dp.tmp.begin(), dp.tmp.end(), int64_t(0));
    }

    // =========================================================================
    // PHASE 2: MAIN DP LOOP (O(N * max_iter) with pruning)
    // =========================================================================
    // Run the main dynamic programming loop to find optimal chains.
    // Extracted to run_main_dp_loop() for readability and testability.
    run_main_dp_loop(
        a.data(), si, ei,
        max_skip, max_iter, max_dis,
        chn_pen_gap, chn_pen_skip, bw_rate,
        xl, yl,
        p, f, t, ii,
        &msc, &msc_i, &movl, &plus);

    // =========================================================================
    // PHASE 3: BACKTRACKING (Extract best chain)
    // =========================================================================
    // Follow predecessor pointers from best anchor (msc_i) to chain root.
    // Extracted to backtrack_best_chain() for readability and testability.
    cL = backtrack_best_chain(msc_i, p, t, ii);

    // =========================================================================
    // PHASE 4: MULTI-COPY EXTRACTION (Optional, for repetitive regions)
    // =========================================================================
    // Extract up to mcopy_num high-scoring chains for repetitive regions.
    // This is essential for:
    //   - Segmental duplications
    //   - Tandem repeats
    //   - Paralogous genes
    //
    // Algorithm:
    //   1. Normalize scores: f[i] -= plus (shift to non-negative)
    //   2. Collect candidates: anchors with score >= mcopy_rate * best_score
    //   3. Sort by score (descending) using radix sort
    //   4. Extract top mcopy_num non-overlapping chains
    //   5. For each chain, emit overlap region + chain hit indices
    if (mcopy_num > 1) {
        // Only attempt mcopy if best chain has enough anchors (quality threshold)
        if (cL >= mcopy_khit_cutoff) {
            // -----------------------------------------------------------------
            // Step 1: Normalize scores to non-negative range
            // -----------------------------------------------------------------
            // Shift all scores by -plus (minimum score) to make them ≥ 0
            // This is required for the radix sort and threshold comparison
            msc -= plus;

            // Safety: Ensure msc is non-negative after normalization
            if (msc < 0) msc = 0;

            // Compute threshold with proper rounding (floor for safe truncation)
            min_sc = static_cast<int64_t>(std::floor(static_cast<double>(msc) * mcopy_rate));

            // Clear best chain flag (allow it to compete with other chains)
            ii[static_cast<size_t>(msc_i)] = 0;

            // -----------------------------------------------------------------
            // Step 2: Collect candidate chains above threshold
            // -----------------------------------------------------------------
            // Build heap keys: (score << 32) | (index << 1)
            // The left-shift by 32 puts score in high bits for radix sort
            for (i = ch_n = 0; i < a_n; ++i) {
                // Normalize score
                f[static_cast<size_t>(i)] -= static_cast<int32_t>(plus);

                // Clear t[] workspace for candidate collection
                if (i >= ch_n) t[static_cast<size_t>(i)] = 0;

                // Collect anchors not in best chain and above threshold
                if ((!ii[static_cast<size_t>(i)]) &&
                    (static_cast<int64_t>(f[static_cast<size_t>(i)]) >= min_sc)) {

                    // Pack (score, index) into 64-bit key for sorting
                    t[static_cast<size_t>(ch_n)] =
                        (static_cast<int64_t>(static_cast<uint64_t>(static_cast<uint32_t>(f[static_cast<size_t>(i)])) << 32) |
                         (static_cast<uint64_t>(i) << 1));
                    ch_n++;
                }
            }

            // -----------------------------------------------------------------
            // Safety: Restore best chain if mcopy collection failed
            // -----------------------------------------------------------------
            // If we found ≤1 candidate, the loop above overwrote t[] (which
            // held the best-chain backtrack). Restore it to preserve hifiasm parity.
            if (ch_n <= 1) {
                msc += plus;  // Restore original score
                i = msc_i;
                cL = 0;
                // Rebuild best chain in t[]
                while (i >= 0) {
                    t[static_cast<size_t>(cL++)] = i;
                    i = p[static_cast<size_t>(i)];
                }
            }
            // -----------------------------------------------------------------
            // Step 3: Sort candidates by score (descending)
            // -----------------------------------------------------------------
            if (ch_n > 1) {
                int64_t n_v, n_v0, ni, n_u;
                const int64_t n_u0 = static_cast<int64_t>(res.size());

                // Sort by packed key: higher score → higher key → later in array
                // We iterate backwards (k = ch_n-1 down to 0) to extract in score order
                std::sort(t, t + ch_n,
                    [](const int64_t a, const int64_t b) noexcept {
                        return static_cast<uint64_t>(a) < static_cast<uint64_t>(b);
                    });

                // -----------------------------------------------------------------
                // Step 4: Extract top mcopy_num chains
                // -----------------------------------------------------------------
                // Iterate backwards through sorted candidates (highest score first)
                for (k = ch_n - 1, n_v = n_u = 0; k >= 0 && n_u < mcopy_num; --k) {
                    n_v0 = n_v;  // Remember start of current chain in ii[]

                    // Extract chain starting from anchor index in t[k]
                    // (index is in bits 1-31, bit 0 is used as visited flag)
                    i = static_cast<int64_t>((static_cast<uint32_t>(t[static_cast<size_t>(k)]) >> 1));

                    // Follow predecessor pointers until we hit a visited anchor
                    for (; i >= 0 && (t[static_cast<size_t>(i)] & 1) == 0; ) {
                        ii[static_cast<size_t>(n_v++)] = static_cast<int32_t>(i);
                        t[static_cast<size_t>(i)] |= 1;  // Mark as visited
                        i = p[static_cast<size_t>(i)];
                    }

                    // Skip if chain is empty (fully overlaps with previous chains)
                    if (n_v0 == n_v) continue;

                    // CRITICAL: Verify chain has valid bounds before accessing
                    // After continue above, we're guaranteed n_v > n_v0, so n_v >= n_v0 + 1
                    DINARA_ASSERT(n_v > n_v0);
                    DINARA_ASSERT(n_v > 0);

                    // Extract score from packed key and adjust for overlap
                    const uint64_t key = static_cast<uint64_t>(t[static_cast<size_t>(k)]);
                    const int64_t top = static_cast<int64_t>(static_cast<uint32_t>(key >> 32));

                    // If we stopped at an existing chain (i >= 0), subtract its score
                    sc = (i < 0) ? top : (top - static_cast<int64_t>(f[static_cast<size_t>(i)]));

                    // Emit chain if it passes threshold
                    if (sc >= min_sc) {
                        HifiasmOverlapRegion z{};

                        // Populate overlap region (coordinates, strand, etc.)
                        hifiasm_push_ovlp_chain_qgen(
                            z, xid, xl, yl, sc + plus,
                            &a[static_cast<size_t>(ii[static_cast<size_t>(n_v - 1)])],
                            &a[static_cast<size_t>(ii[static_cast<size_t>(n_v0)])]);

                        // Emit chain if it's the first or has >1 anchor
                        if ((!n_u) || (n_v - n_v0 > 1)) {
                            z.align_length = static_cast<uint32_t>(n_v - n_v0);
                            z.x_id = static_cast<uint32_t>(n_v0);  // Temp: start index in ii[]

                            // Store chain hit indices in increasing self_offset order
                            z.non_homopolymer_errors = static_cast<uint32_t>(chainHitIndexFlat.size());
                            ni = z.align_length;

                            // Copy indices in reverse order (ii[] is in reverse)
                            // Use uint32_t for iteration to prevent underflow in (ni - j - 1)
                            for (uint32_t j = 0; j < ni; ++j) {
                                const int64_t local = ii[static_cast<size_t>(n_v0 + (ni - j - 1))];
                                chainHitIndexFlat.push_back(a[static_cast<size_t>(local)].globalIndex);
                            }

                            res.push_back(z);
                            n_u++;
                        } else {
                            // Chain is too short (only 1 anchor), discard
                            n_v = n_v0;
                        }
                    } else {
                        // Score below threshold, discard
                        n_v = n_v0;
                    }
                }

                // If we successfully emitted mcopy chains, we're done
                if (res.size() > static_cast<size_t>(n_u0)) {
                    return;
                }

                // -----------------------------------------------------------------
                // Fallback: Restore best chain if mcopy extraction failed
                // -----------------------------------------------------------------
                msc += plus;  // Restore original score
                i = msc_i;
                cL = 0;
                while (i >= 0) {
                    t[static_cast<size_t>(cL++)] = i;
                    i = p[static_cast<size_t>(i)];
                }
            }
        }
    }

    // =========================================================================
    // PHASE 5: EMIT BEST CHAIN (Always executed, fallback if mcopy failed)
    // =========================================================================
    // Emit the best chain (msc, msc_i) as an overlap region.
    // Extracted to emit_best_chain_as_overlap() for readability and testability.
    emit_best_chain_as_overlap(a.data(), t, cL, msc, xid, xl, yl, res, chainHitIndexFlat);
}

/**
 * @brief Classify overlap type for filtering and statistics (hifiasm parity).
 *
 * Determines how the query read aligns relative to its full length:
 *   - Type 0: Query prefix (overhang at end)
 *   - Type 1: Query suffix (overhang at start)
 *   - Type 2: Query fully contained (covered end-to-end)
 *   - Type 3: Query contains target (internal alignment)
 *
 * This classification is used downstream for:
 *   - Filtering overlaps by topology
 *   - Computing containment relationships
 *   - Prioritizing overlap types during graph construction
 *
 * @param r    Overlap region to classify
 * @param len  Query read length
 * @return Overlap type (0=prefix, 1=suffix, 2=contained, 3=contains)
 *
 * @complexity O(1)
 * @reference Hifiasm anchor.cpp:86 (ha_ov_type)
 */
static inline int hifiasm_ha_ov_type(
    const HifiasmOverlapRegion& r,
    const uint32_t len) noexcept
{
    // Type 2: Query is fully covered [0, len-1]
    if (r.x_pos_s == 0 && r.x_pos_e == len - 1U) {
        return 2;  // Contained
    }
    // Type 3: Query alignment is internal (not touching either end)
    else if (r.x_pos_s > 0 && r.x_pos_e < len - 1U) {
        return 3;  // Contains
    }
    // Type 0 or 1: Query touches one end
    else {
        return (r.x_pos_s == 0) ? 0 : 1;  // Prefix : Suffix
    }
}

// ============================================================================
// POST-FILTER: MAX_N_CHAIN + COV_W + R485 WEAK-CHAIN SUPPRESSION
// ============================================================================
// Reference: Hifiasm anchor.cpp lchain_qgen_mcopy_fast max_n_chain + ocv_w + r485
//
// After the DP chaining phase produces overlap regions, this post-filter applies
// three successive pruning stages to control output volume and quality:
//
//   ┌─────────────────────────────────────────────────────────────────┐
//   │  Stage 1: MAX_N_CHAIN (per-overlap-type cap)                   │
//   │    • Sort overlaps by score descending                          │
//   │    • For each of the 4 overlap types (prefix/suffix/cont/int), │
//   │      record the score at the max_n_chain-th position.           │
//   │    • Reject overlaps below their type's threshold.              │
//   │                                                                 │
//   │  Stage 2: COV_W (Coverage Window overload control)              │
//   │    • Only applies to type-3 (internal/containing) overlaps      │
//   │    • Divides the query read into fixed-width windows            │
//   │    • Each window has a capacity = window_size * (max_n_chain/2) │
//   │    • A type-3 overlap below threshold is "rescued" if ≥70% of   │
//   │      its span falls in under-capacity windows                    │
//   │    • This prevents repetitive regions from drowning out signal   │
//   │                                                                 │
//   │  Stage 3: R485 (Weak-chain suppression)                         │
//   │    • Chains with few anchors (< chain_cutoff) are "weak"        │
//   │    • If a weak chain's query span overlaps a strong chain's     │
//   │      span by ≥ 95%, and the strong chain has ≥ 16× more anchors │
//   │      and ≥ 16× higher score, suppress the weak chain            │
//   │    • Prevents short false-positive chains from shadowing real    │
//   │      overlaps in repetitive regions                              │
//   └─────────────────────────────────────────────────────────────────┘
//
// @param ol             Vector of overlap regions (modified in place, may shrink)
// @param max_n_chain    Per-type cap on number of overlaps to keep
// @param chain_cutoff   Minimum anchor count for "strong" chains in R485
// @param ocv_w          COV_W window size (0 to disable)
// @param rl             Query read length
// @param chainHitIndexFlat  Flattened anchor indices for each chain
// @param allHits        Full anchor array (for R485 position lookups)
// ============================================================================
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

    uint64_t lch = 0;  // Counts overlaps that pass all filters.

    // =====================================================================
    // STAGE 1: MAX_N_CHAIN — Per-type score threshold
    // =====================================================================
    // Sort all overlaps by score descending. Then make a single pass to find,
    // for each of the 4 overlap types, the score at position max_n_chain.
    // Any overlap below its type's threshold is rejected (unless rescued).
    if (max_n_chain > 0 && ol.size() > max_n_chain) {
        std::sort(ol.begin(), ol.end(),
            [](const HifiasmOverlapRegion& a, const HifiasmOverlapRegion& b) {
                return a.shared_seed > b.shared_seed;
            });

        // n[w] = count of overlaps of type w seen so far.
        // s[w] = score at the max_n_chain-th overlap of type w (threshold).
        int32_t n[4] = {0, 0, 0, 0};
        int32_t s[4] = {0, 0, 0, 0};
        for (size_t i = 0; i < ol.size(); ++i) {
            const int w = hifiasm_ha_ov_type(ol[i], rl);
            ++n[w];
            if (uint64_t(n[w]) == max_n_chain) s[w] = ol[i].shared_seed;
        }

        // Only proceed if at least one type hit its cap.
        if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
            // =============================================================
            // STAGE 2: COV_W — Coverage Window rescue for type-3 overlaps
            // =============================================================
            // Divide the query read into fixed-width windows (ocv_w bases each).
            // Each window tracks how much alignment span has been consumed.
            // A below-threshold type-3 overlap is "rescued" if at least 70%
            // of its span falls in windows that still have capacity.
            //
            // Window state packing (uint64_t):
            //   High 32 bits = capacity (window_size * max_n_chain/2)
            //   Low 32 bits  = used (sum of overlap spans assigned so far)
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

                    // Safely compute overlap span with overflow protection
                    const uint64_t overlap_span = oe - os;
                    const uint32_t add = (overlap_span > uint64_t(UINT32_MAX)) ?
                        UINT32_MAX : uint32_t(overlap_span);

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

            // Apply max_n_chain + COV_W filtering in a single pass.
            // Overlaps above their type's threshold are always kept.
            // Type-3 overlaps below threshold get a COV_W rescue check.
            // Kept overlaps that have < chain_cutoff anchors mark lch=1
            // to trigger R485 weak-chain suppression below.
            size_t k = 0;
            for (size_t i = 0; i < ol.size(); ++i) {
                const int w = hifiasm_ha_ov_type(ol[i], rl);
                bool keep = (ol[i].shared_seed >= s[w]);
                if (!keep && w == 3 && cwn > 0) {
                    keep = should_rescue_type3(ol[i]);  // COV_W rescue.
                }
                if (keep) {
                    update_cc(ol[i]);  // Update window usage for future rescue checks.
                    if (chain_cutoff >= 2 && ol[i].align_length < chain_cutoff) lch = 1;
                    if (k != i) std::swap(ol[k], ol[i]);
                    ++k;
                }
            }
            ol.resize(k);
        }
    }

    // =================================================================
    // STAGE 3: R485 WEAK-CHAIN SUPPRESSION
    // =================================================================
    // Sort surviving overlaps by query start position so we can do
    // efficient interval overlap queries between weak and strong chains.
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
        // R485 weak-chain suppression (anchor.cpp:2061-2096).
        //
        // For each weak chain (align_length < chain_cutoff), check if there exists
        // any strong chain whose query span overlaps by ≥ OFL (95%), has ≥ CH_OCC
        // (16×) more anchors, has ≥ CH_SC (16×) higher score, AND has enough actual
        // anchor positions within the overlap zone. If so, suppress the weak chain.
        //
        // This prevents short, spurious chains (often caused by repetitive k-mers)
        // from diluting the overlap graph with false positives.
        size_t l = 0;
        for (size_t i = 0; i < ol.size(); ++i) {
            if (ol[i].align_length < chain_cutoff) {
                // Weak chain: compute its query span [zs, ze) and thresholds.
                const uint64_t zs = uint64_t(ol[i].x_pos_s);
                const uint64_t ze = uint64_t(ol[i].x_pos_e) + 1ULL;
                uint64_t ob = uint64_t(double(ze - zs) * HIFIASM_OFL);  // 95% overlap needed.
                if (ob < 16) ob = 16;  // Minimum overlap of 16 bases.
                const int64_t osc = int64_t(ol[i].shared_seed) * int64_t(HIFIASM_CH_SC);  // Score threshold.
                const uint64_t ocn = uint64_t(ol[i].align_length) << HIFIASM_CH_OCC;  // Anchor count threshold.

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

                    // CRITICAL: Check bounds before accessing chain hit indices
                    if(off + n_hit > chainHitIndexFlat.size()) {
                        continue;  // Skip this overlap - hit indices out of bounds
                    }

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
// ============================================================================
// MARKER ORDINAL RESOLUTION
// ============================================================================
// The inverted index stores hits by base position, but Dinara's Alignment
// structure needs marker ordinals (sequential index into the marker array).
// This function bridges the two representations by binary-searching the
// sorted marker array for each hit's base position on Read B.
//
// Why not resolve during hit collection?
//   Read A ordinals are known at collection time (we iterate A's markers).
//   But Read B ordinals require looking up each B-position in B's marker
//   array, which is cheaper to batch after sorting hits by posB.
//
// Algorithm:
//   1. Sort hit indices by posB → enables a single linear scan.
//   2. Walk through sorted markers and sorted hits in tandem.
//   3. For each hit, advance the marker cursor until position matches.
//   4. Record the matching marker ordinal.
//
// Returns false if any hit position has no corresponding marker (data error).
// ============================================================================
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

    // Step 1: Create an index permutation sorted by posB.
    // This lets us scan markers linearly instead of binary-searching each hit.
    orderByPosB.resize(n);
    std::iota(orderByPosB.begin(), orderByPosB.end(), uint32_t(0));
    std::sort(orderByPosB.begin(), orderByPosB.end(),
        [&](const uint32_t a, const uint32_t b) {
            if (hitPosB[a] != hitPosB[b]) {
                return hitPosB[a] < hitPosB[b];
            }
            return a < b;
        });

    // Step 2: Linear tandem scan — advance marker cursor and hit cursor together.
    size_t markerIdx = 0;
    for (const uint32_t hitIdx : orderByPosB) {
        const uint32_t pos = hitPosB[hitIdx];

        // Advance marker cursor past positions smaller than this hit's position.
        while (markerIdx < markersB.size() && markersB[markerIdx].position < pos) {
            ++markerIdx;
        }

        // The marker at markerIdx should now have position == pos.
        // If not, the hit references a position that doesn't exist as a marker.
        if (markerIdx >= markersB.size() || markersB[markerIdx].position != pos) {
            return false;  // Data integrity error — position not found.
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
        // Check LUT bounds with explicit type for clarity
        if (w < static_cast<uint32_t>(weightLut.size())) {
            return weightLut[w];
        }
        // Fallback: w beyond LUT range, compute dynamically
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
// - Overlap validation: nonRedundantOverlapFraction
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
        vector<HifiasmKmerHit> hifiasmPairHitsStrand0;
        vector<HifiasmKmerHit> hifiasmPairHitsStrand1;
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

                    // Safety: Validate sample distance to prevent division by zero
                    if (highFrequencySampleDistance == 0) {
                        // Should never happen (protected by max(1U, ...) at initialization)
                        // But if config is corrupt, keep all markers as fallback
                        for(const auto& m : highFrequencyStreak) {
                            appendMarkerHits(m);
                        }
                        highFrequencyStreak.clear();
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
                            int64_t(mcopyKhitCutoff));
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

                        // CRITICAL: Check bounds before accessing chain hit indices
                        if(off + n_hit > hifiasmChainHitIndexFlat.size()) {
                            continue;  // Skip this overlap - hit indices out of bounds
                        }

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

                        // Optional: reject candidates with small pre-extension overlap span.
                        if(invertedIndexData.minOverlapLength > 0) {
                            const uint32_t qSpan = qE - qS;
                            const uint32_t tSpan = tE - tS;
                            if(std::min(qSpan, tSpan) < invertedIndexData.minOverlapLength) {
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

                        // Canonicalize candidate so readIds[0] < readIds[1], and keep alignment consistent.
                        ReadId cand0 = readIdA;
                        ReadId cand1 = readIdB;
                        canonicalizeCandidateAndAlignment(cand0, cand1, isSameStrand, al, markerCountA, markerCountB);

                        const uint8_t overlapType =
                            uint8_t(getOverlapType(qS, qE, uint32_t(readLenA)));

                        emittedForRead.push_back(EmittedChainedCandidate{
                            OrientedReadPair(cand0, cand1, isSameStrand),
                            std::move(al),
                            r.shared_seed,
                            overlapType,
                            (cand0 == readIdA) ? cand1 : cand0,
                            0});
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

	                        // Safety: Validate sample distance to prevent division by zero
	                        if (highFrequencySampleDistance == 0) {
	                            // Should never happen (protected by max(1U, ...) at initialization)
	                            // But if config is corrupt, keep all markers as fallback
	                            for(const auto& m : highFrequencyStreak) {
	                                appendMarkerHits(m);
	                            }
	                            highFrequencyStreak.clear();
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
