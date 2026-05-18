#ifndef DINARA_PHASING_KMEANS_TYPES_HPP
#define DINARA_PHASING_KMEANS_TYPES_HPP

/// @file PhasingKmeansTypes.hpp
/// @brief Data structures for k-means overlap phasing.
///
/// Adapts pgphase/longcallD's iterative k-means read-haplotype clustering
/// to operate on dinara's OverlapCigarStore pairwise CIGARs.
///
/// Key difference from pgphase: there is no reference genome. Each backbone
/// read acts as its own reference. Variant sites are detected from pairwise
/// overlap CIGARs projected onto the backbone read's coordinate system.

#include "cstdint.hpp"
#include "vector.hpp"

namespace dinara {

// ============================================================================
// Variant categories (adapted from pgphase/longcallD)
// ============================================================================

/// Classification of a candidate variant site on the backbone read.
enum class KmVariantCategory : uint8_t {
    LowCoverage,       ///< Too few reads cover this position.
    LowAlleleFraction, ///< Alt allele fraction below threshold.
    StrandBias,        ///< Strand-biased alt observations.
    CleanHetSnp,       ///< Clean heterozygous SNP (balanced AF, not in repeat).
    CleanHetIndel,     ///< Clean heterozygous indel (balanced AF, not in repeat).
    CleanHom,          ///< Homozygous variant (AF > max_af).
    RepeatHetIndel,    ///< Indel in homopolymer/tandem repeat context.
    NoisyCandHet,      ///< Noisy-region het (from MSA recall, phase 2).
    NoisyCandHom,      ///< Noisy-region hom (from MSA recall, phase 2).
    NonVariant         ///< Demoted: inside noisy region or filtered.
};

/// Bitmask flags for category filtering in k-means (mirrors longcallD).
static constexpr uint32_t KM_CLEAN_HET_SNP   = 0x004u;
static constexpr uint32_t KM_CLEAN_HET_INDEL = 0x008u;
static constexpr uint32_t KM_REP_HET_VAR     = 0x010u;
static constexpr uint32_t KM_CLEAN_HOM        = 0x080u;
static constexpr uint32_t KM_NOISY_CAND_HET  = 0x100u;
static constexpr uint32_t KM_NOISY_CAND_HOM  = 0x200u;
static constexpr uint32_t KM_NON_VAR          = 0x800u;

/// Clean germline categories (used for initial k-means pass).
static constexpr uint32_t KM_GERMLINE_CLEAN =
    KM_CLEAN_HET_SNP | KM_CLEAN_HET_INDEL | KM_CLEAN_HOM;

/// All germline categories including noisy (used after MSA refinement).
static constexpr uint32_t KM_GERMLINE_ALL =
    KM_GERMLINE_CLEAN | KM_NOISY_CAND_HET | KM_NOISY_CAND_HOM;

inline uint32_t kmCategoryToFlag(KmVariantCategory c) {
    switch (c) {
        case KmVariantCategory::CleanHetSnp:   return KM_CLEAN_HET_SNP;
        case KmVariantCategory::CleanHetIndel: return KM_CLEAN_HET_INDEL;
        case KmVariantCategory::CleanHom:       return KM_CLEAN_HOM;
        case KmVariantCategory::RepeatHetIndel: return KM_REP_HET_VAR;
        case KmVariantCategory::NoisyCandHet:  return KM_NOISY_CAND_HET;
        case KmVariantCategory::NoisyCandHom:  return KM_NOISY_CAND_HOM;
        case KmVariantCategory::NonVariant:     return KM_NON_VAR;
        default: return 0;
    }
}

// ============================================================================
// Candidate variant site
// ============================================================================

/// A candidate variant site on the backbone read.
struct KmCandidate {
    uint32_t pos;              ///< Position on backbone read (base coords).
    uint8_t  refBase;          ///< Backbone base at this position (0-3=ACGT).
    uint8_t  type;             ///< 0=SNP, 1=insertion, 2=deletion.
    uint16_t refLen;           ///< Reference length consumed (1 for SNP, >1 for del).

    // Allele counts.
    uint32_t totalCov;         ///< Total overlaps covering this position.
    uint32_t refCov;           ///< Overlaps matching backbone allele.
    uint32_t altCov;           ///< Overlaps with alt allele (primary alt).
    uint32_t lowQualCov;       ///< Low-quality alt observations.
    uint32_t fwdRef;           ///< Forward-strand ref-matching overlaps.
    uint32_t revRef;           ///< Reverse-strand ref-matching overlaps.
    uint32_t fwdAlt;           ///< Forward-strand alt overlaps.
    uint32_t revAlt;           ///< Reverse-strand alt overlaps.

    /// Multi-allelic: per-allele coverage. Index 0 = ref, 1 = primary alt, 2+ = other alts.
    /// Empty when biallelic (use refCov/altCov directly).
    vector<int> alleCovs;
    int nUniqAlles = 2;

    double alleleFraction;     ///< altCov / totalCov.
    KmVariantCategory category;
    uint32_t categoryFlag;     ///< Bitmask flag for k-means filtering.

    bool isHomopolymerIndel = false;

    // K-means phasing state.
    /// Per-haplotype allele profile: alleCovs-sized vector per hap (0=unassigned, 1, 2).
    vector<int> hapAlleProfile[3];
    /// Consensus allele per hap (-1 = unknown).
    int hapConsAlle[3] = {-1, -1, -1};
    int32_t phaseSet = 0;
    int hapAlt = 0;            ///< Which hap carries the alt allele (1 or 2, 0=unresolved).
    int hapRef = 0;            ///< Which hap carries the ref allele.
};

// ============================================================================
// Per-overlap allele profile
// ============================================================================

/// Sparse allele observations for one overlap across candidate sites.
struct KmOverlapProfile {
    uint32_t overlapIdx;       ///< Index into KmScratchpad::overlaps.
    int startVarIdx = -1;      ///< First candidate index this overlap covers.
    int endVarIdx = -1;        ///< Last candidate index this overlap covers.
    /// Allele at each candidate in [startVarIdx, endVarIdx]:
    ///   0 = ref, 1 = primary alt, -1 = no observation, -2 = low-quality alt.
    vector<int> alleles;
};

// ============================================================================
// Noisy region interval
// ============================================================================

/// A contiguous region on the backbone read where CIGAR evidence is ambiguous.
struct KmNoisyRegion {
    uint32_t start;            ///< Start position on backbone read.
    uint32_t end;              ///< End position (exclusive).
    int label = 0;             ///< Metadata (variant density, merge info).
    bool done = false;         ///< True after MSA processing.
};

// ============================================================================
// Per-overlap phasing state (reused from PhasingTypes.hpp pattern)
// ============================================================================

/// Lightweight overlap descriptor for k-means phasing.
struct KmOverlap {
    uint32_t alignmentId;      ///< Index into assembler.alignmentData.
    uint32_t targetReadId;     ///< Partner read ID.
    uint32_t qs, qe;          ///< Query (backbone) start/end in base coords.
    uint32_t ts, te;          ///< Target start/end in base coords.
    uint32_t cigarOffset;
    uint32_t cigarTokenCount;
    uint8_t  isRev;            ///< 1 if target is reverse-complemented.
    uint8_t  queryIsRead0;     ///< 1 if backbone == readIds[0] in AlignmentData.
    int      hap = 0;          ///< Assigned haplotype: 0=unassigned, 1=hap1, 2=hap2.
    int32_t  phaseSet = -1;    ///< Phase set assignment.
    uint8_t  strong = 0;       ///< 1 = confidently phased.
};

// ============================================================================
// Thread-local scratchpad
// ============================================================================

/// Reusable workspace for k-means phasing of one backbone read.
struct KmScratchpad {
    // Overlaps for this backbone read.
    vector<KmOverlap> overlaps;

    // Unpacked backbone sequence.
    vector<uint8_t> backboneBases;

    // Candidate variant sites (sorted by position).
    vector<KmCandidate> candidates;

    // Per-overlap allele profiles.
    vector<KmOverlapProfile> overlapProfiles;

    // Noisy regions detected from classification.
    vector<KmNoisyRegion> noisyRegions;

    // Mismatch vote counts per backbone position (for candidate detection).
    vector<uint8_t> mismatchVotes;

    // Indel event positions per backbone position.
    vector<uint8_t> indelVotes;

    // Per-overlap mismatch records for allele counting.
    struct MismatchRecord {
        uint32_t backbonePos;  ///< Position on backbone read.
        uint32_t targetPos;    ///< Position on target read.
        uint32_t overlapIdx;   ///< Which overlap produced this.
    };
    vector<MismatchRecord> mismatchRecords;

    // Per-overlap indel records.
    struct IndelRecord {
        uint32_t backboneStart; ///< Start on backbone.
        uint32_t backboneEnd;   ///< End on backbone (exclusive).
        uint32_t overlapIdx;
        uint8_t  type;          ///< 2=insertion, 3=deletion.
        uint16_t len;           ///< Length of indel.
    };
    vector<IndelRecord> indelRecords;

    // Sorted overlap indices for k-means pivot selection.
    vector<uint32_t> validVarIdx;

    void clear() {
        overlaps.clear();
        backboneBases.clear();
        candidates.clear();
        overlapProfiles.clear();
        noisyRegions.clear();
        mismatchVotes.clear();
        indelVotes.clear();
        mismatchRecords.clear();
        indelRecords.clear();
        validVarIdx.clear();
    }
};

// ============================================================================
// Thresholds
// ============================================================================

/// Configuration for k-means phasing pipeline.
struct KmPhasingOptions {
    uint32_t minDepth = 5;
    uint32_t minAltDepth = 2;
    double   minAf = 0.20;
    double   maxAf = 0.80;
    double   strandBiasPval = 0.01;
    uint32_t noisyRegMergeDis = 500;
    uint32_t maxKmeansIter = 10;
    uint32_t minSpanningReads = 2; ///< Min reads spanning adjacent het sites for phase-set continuity.
};

} // namespace dinara

#endif // DINARA_PHASING_KMEANS_TYPES_HPP
