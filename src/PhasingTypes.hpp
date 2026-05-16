#ifndef DINARA_PHASING_TYPES_HPP
#define DINARA_PHASING_TYPES_HPP

/// @file PhasingTypes.hpp
/// @brief Data structures for hifiasm-parity ONT overlap phasing.
///
/// Ports hifiasm's rphase_hc pipeline (Correct.cpp) to operate on
/// Dinara's OverlapCigarStore instead of hifiasm's in-memory CIGAR arrays.
///
/// The phasing pipeline detects heterozygous SNP sites from overlap CIGARs,
/// builds an overlap × site matrix, runs DP-based phasing to find consistent
/// haplotype chains, and labels each overlap as cis (same haplotype) or
/// trans (different haplotype).

#include "cstdint.hpp"
#include "vector.hpp"

namespace dinara {

// ============================================================================
// Constants (from hifiasm Correct.cpp / ecovlp.cpp)
// ============================================================================

/// Sliding window size for SNP detection.
/// Hifiasm: WINDOW_MAX_SIZE = 375 + (int)(1.0/0.02) + 3 = 428.
static constexpr uint32_t PHASING_WINDOW_SIZE = 428;

/// Minimum mismatch vote count to flag a position as candidate SNP.
/// A position needs ≥ (occ_thres + 1) = 2 overlaps mismatching to be flagged.
static constexpr uint32_t PHASING_OCC_THRES = 1;

/// Minimum match count (reference allele) for a confirmed SNP site.
static constexpr uint32_t PHASING_S_HAP_COV = 3;

/// Minimum alt allele count for a confirmed SNP site.
static constexpr uint32_t PHASING_INFOR_COV = 3;

/// Minimum indel size for SV-based phasing (hifiasm rphase_lidel).
static constexpr uint32_t PHASING_SV_MIN_LEN = 16;


// ============================================================================
// Per-overlap phasing state
// ============================================================================

/// Lightweight view of one overlap for phasing purposes.
/// Populated once per query read from AlignmentData + AlignmentInfo.
struct PhasingOverlap {
    uint32_t alignmentId;     ///< Index into assembler.alignmentData
    uint32_t targetReadId;    ///< Partner read ID
    uint32_t qs;              ///< Query start (forward-strand base coords)
    uint32_t qe;              ///< Query end
    uint32_t ts;              ///< Target start (forward-strand base coords)
    uint32_t te;              ///< Target end
    uint32_t cigarOffset;     ///< Offset into OverlapCigarStore token arena
    uint32_t cigarTokenCount; ///< Number of CIGAR tokens
    uint32_t errorCount;      ///< Non-homopolymer error count (for dedup scoring)
    uint8_t  isRev;           ///< 1 if target is reverse-complemented
    uint8_t  isMatch;         ///< 0=unclassified, 1=cis, 2=trans
    uint8_t  strong;          ///< 1=confidently phased
    uint8_t  queryIsRead0;    ///< 1 if query == readIds[0] in AlignmentData

    /// Number of confirmed SNP sites this overlap covers (for scoring).
    uint32_t confirmedSiteCount = 0;

    /// Number of mismatches at confirmed SNP sites.
    uint32_t confirmedMismatchCount = 0;
};

// ============================================================================
// SNP evidence
// ============================================================================

/// One observation at a candidate SNP site from one overlap.
/// Produced during pass 2 of the sliding window.
struct PhasingEvidence {
    uint32_t site;        ///< Position on query read (forward-strand base coords)
    uint32_t overlapIdx;  ///< Index into PhasingScratchpad::overlaps
    uint32_t siteIdx;     ///< Index into PhasingScratchpad::sites (hifiasm overlapSite)
    uint8_t  base;        ///< 2-bit base observed at this position (0=A,1=C,2=G,3=T)
    uint8_t  isAlt;       ///< 1 if this observation differs from the query base
    uint8_t  isHpc;       ///< 1 if position is in a periodic repeat (hpc_mask_ff)
};

/// Confirmed heterozygous SNP site on the query read.
///
/// Hifiasm reuses SnpStats fields (occ_0, overlap_num, score) across pipeline
/// stages with different semantics. We use separate fields for clarity:
///
///   matchCount / altCount / fwdStrandCount — immutable after buildSnpMatrix.
///   dpChainId — output of runDpPhasing (-1 = rejected, >=0 = chain ID).
///   labelMatchCount / labelFwdStrandCount — mutable copies decremented during
///       greedy labeling (hifiasm mutates occ_0/overlap_num in-place).
///
/// Hifiasm's `score` field is reused across labelCisTrans steps with three
/// different meanings. We split it into three separate flags:
///
///   transConfirmed — Step B: set when a trans overlap has a mismatch here.
///   cisReset       — Step D: set when a cis overlap has a mismatch here
///                    (overrides transConfirmed for multi_check filtering).
///   promoted       — Step E: set when multi_check promotes this weak site.
///
/// The combined check `isLabelConfirmed()` returns true when the site should
/// be treated as confirmed for Step F (final loop): transConfirmed and not
/// cisReset, OR promoted.
struct PhasingSite {
    uint32_t site;            ///< Position on query read
    uint8_t  queryBase;       ///< Query (reference) allele (2-bit)
    uint8_t  altBase;         ///< Alternative allele (2-bit)
    uint8_t  isHpc;           ///< Homopolymer context flag
    int32_t  dpChainId;       ///< DP chain assignment (-1 = unassigned)

    // --- Immutable counts (set by buildSnpMatrix) ---
    uint32_t matchCount;      ///< Overlaps matching query allele (occ_0)
    uint32_t altCount;        ///< Overlaps with alt allele (occ_1)
    uint32_t fwdStrandCount;  ///< Forward-strand ref-matching overlaps (+1 for query)

    // --- Mutable copies for greedy labeling (labelCisTrans) ---
    uint32_t labelMatchCount;    ///< Decremented as trans overlaps are processed
    uint32_t labelFwdStrandCount;///< Decremented for fwd-strand matches at trans overlaps

    // --- Per-step labeling flags (replaces hifiasm's reused `score`) ---
    uint8_t  transConfirmed;  ///< Step B: mismatch at this site in a trans overlap
    uint8_t  cisReset;        ///< Step D: mismatch at this site in a cis overlap
    uint8_t  promoted;        ///< Step E: multi_check promoted this weak site

    /// Combined check: is this site confirmed for Step F?
    /// Hifiasm equivalent: score == 1 at the point Step F runs.
    bool isLabelConfirmed() const {
        return (transConfirmed && !cisReset) || promoted;
    }

    /// Range into sorted evidence array for this site.
    uint32_t evidenceBegin;
    uint32_t evidenceEnd;
};

// ============================================================================
// Large indel (SV) evidence
// ============================================================================

/// A contiguous error region ≥ SV_MIN_LEN detected from a single overlap's CIGAR.
struct PhasingSvEvent {
    uint32_t overlapIdx;  ///< Index into PhasingScratchpad::overlaps
    uint32_t queryPos;    ///< Start position on query read
    uint32_t queryEnd;    ///< End position on query read
    uint32_t errorBases;  ///< Total insertion + deletion bases in this region
};

/// A cluster of SV events at similar positions across overlaps.
struct PhasingSvCluster {
    uint32_t consensusPos;    ///< Consensus start position
    uint32_t consensusEnd;    ///< Consensus end position
    uint32_t eventCount;      ///< Number of overlaps contributing
    uint32_t eventBegin;      ///< Range into sorted svEvents array
    uint32_t eventEnd;
};

// ============================================================================
// Thread-local scratchpad
// ============================================================================

/// Reusable workspace for phasing one query read.
/// Allocated once per thread, cleared between reads via clear().
struct PhasingScratchpad {

    // --- Per-query-read overlap list ---
    vector<PhasingOverlap> overlaps;

    // --- SNP detection (sliding window) ---
    vector<uint8_t> flag;             ///< Per-position mismatch vote count (window-local)
    vector<uint8_t> queryBases;       ///< Unpacked query sequence (2-bit → uint8_t)
    vector<uint32_t> candidatePositions; ///< Sorted candidate SNP positions

    /// Mismatch record: query position + target position from CIGAR walk.
    struct MismatchRecord {
        uint32_t qpos;
        uint32_t tpos;
    };
    vector<MismatchRecord> mismatchRecords;  ///< All mismatch positions per overlap
    vector<uint32_t> mismatchRangeBegin;     ///< Per-overlap range into mismatchRecords

    /// Indel record: query position range [qpos, tpos) affected by indel.
    /// (tpos field reused as end position for compactness.)
    struct IndelRecord {
        uint32_t qpos;  ///< Start of query range
        uint32_t tpos;  ///< End of query range (exclusive)
    };
    vector<IndelRecord> indelRecords;    ///< All indel ranges per overlap
    vector<uint32_t> indelRangeBegin;    ///< Per-overlap range into indelRecords

    vector<uint8_t> mismatchBitmap;  ///< [oi * numCand + ci] = 1 if overlap oi mismatches at candidate ci
    vector<uint8_t> indelBitmap;     ///< [oi * numCand + ci] = 1 if candidate ci is in indel for overlap oi

    // --- SNP evidence ---
    vector<PhasingEvidence> evidence; ///< All evidence entries (sorted by site, overlapIdx)
    vector<PhasingSite> sites;        ///< Confirmed SNP sites

    // --- DP phasing ---
    vector<int32_t> dpScore;          ///< DP score per site
    vector<int32_t> dpParent;         ///< DP predecessor per site
    vector<int32_t> dpChainId;        ///< Chain assignment per site
    vector<bool>    dpIsEndpoint;     ///< Endpoint flags for chain extraction
    vector<pair<int32_t, uint32_t>> dpEndpoints; ///< (score, index) for chain extraction
    vector<uint32_t> dpChain;         ///< Backtrack buffer for chain extraction

    // --- Allele grouping ---
    vector<uint32_t> sortedOverlapIndices; ///< Overlaps sorted by mismatch count

    // --- Sorting scratch buffers ---
    vector<PhasingEvidence> evidenceTmp;  ///< Output / scatter buffer for counting sort
    vector<PhasingEvidence> sortTmp;      ///< Per-site sort scatter buffer
    vector<uint32_t> countBuf;            ///< Reusable counting sort histogram

    // --- Large indel phasing ---
    vector<PhasingSvEvent> svEvents;
    vector<PhasingSvCluster> svClusters;

    /// Reset all vectors without releasing memory.
    void clear() {
        overlaps.clear();
        flag.clear();
        queryBases.clear();
        candidatePositions.clear();
        mismatchRecords.clear();
        mismatchRangeBegin.clear();
        indelRecords.clear();
        indelRangeBegin.clear();
        mismatchBitmap.clear();
        indelBitmap.clear();
        evidence.clear();
        sites.clear();
        dpScore.clear();
        dpParent.clear();
        dpChainId.clear();
        dpIsEndpoint.clear();
        dpEndpoints.clear();
        dpChain.clear();
        sortedOverlapIndices.clear();
        evidenceTmp.clear();
        sortTmp.clear();
        countBuf.clear();
        svEvents.clear();
        svClusters.clear();
    }
};

} // namespace dinara

#endif // DINARA_PHASING_TYPES_HPP
