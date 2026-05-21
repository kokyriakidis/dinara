#ifndef DINARA_PHASING_KMEANS_TYPES_HPP
#define DINARA_PHASING_KMEANS_TYPES_HPP

/// @file PhasingKmeansTypes.hpp
/// @brief Data structures for k-means overlap phasing (pgphase/longcallD-style).

#include "cstdint.hpp"
#include "vector.hpp"
#include <array>
#include <cstring>

namespace dinara {

// ============================================================================
// Variant categories
// ============================================================================

enum class KmVariantCategory : uint8_t {
    LowCoverage, LowAlleleFraction, StrandBias,
    CleanHetSnp, CleanHetIndel, CleanHom,
    RepeatHetIndel, NoisyCandHet, NoisyCandHom, NonVariant
};

static constexpr uint32_t KM_CLEAN_HET_SNP   = 0x004u;
static constexpr uint32_t KM_CLEAN_HET_INDEL = 0x008u;
static constexpr uint32_t KM_REP_HET_VAR     = 0x010u;
static constexpr uint32_t KM_CLEAN_HOM        = 0x080u;
static constexpr uint32_t KM_NOISY_CAND_HET  = 0x100u;
static constexpr uint32_t KM_NOISY_CAND_HOM  = 0x200u;
static constexpr uint32_t KM_NON_VAR          = 0x800u;

static constexpr uint32_t KM_GERMLINE_CLEAN =
    KM_CLEAN_HET_SNP | KM_CLEAN_HET_INDEL | KM_CLEAN_HOM;
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
// Variant type (mirrors pgphase VariantType)
// ============================================================================

enum class KmVarType : uint8_t { Snp = 0, Insertion = 1, Deletion = 2 };

// ============================================================================
// DigarOp: per-read variant event (mirrors pgphase DigarOp)
// ============================================================================

/// One variant event on the backbone read, parsed from a pairwise CIGAR.
/// Equivalent to pgphase's DigarOp. Sorted by backbone position within each overlap.
struct KmDigarOp {
    uint32_t pos;          ///< Backbone position (0-based).
    KmVarType type;        ///< SNP, insertion, or deletion.
    uint8_t  altBase;      ///< For SNPs: alt base (0-3=ACGT). For indels: unused.
    uint16_t len;          ///< For indels: length. For SNPs: 1.

    /// Sort key: pos first, then type, then altBase, then len.
    bool operator<(const KmDigarOp& o) const {
        if (pos != o.pos) return pos < o.pos;
        if (type != o.type) return type < o.type;
        if (altBase != o.altBase) return altBase < o.altBase;
        return len < o.len;
    }
    bool operator==(const KmDigarOp& o) const {
        return pos == o.pos && type == o.type && altBase == o.altBase && len == o.len;
    }
};

// ============================================================================
// Variant key: candidate identity (mirrors pgphase VariantKey)
// ============================================================================

/// Unique identity of a candidate variant site.
/// Two events with the same key are the same variant.
struct KmVarKey {
    uint32_t pos;
    KmVarType type;
    uint8_t  altBase;      ///< For SNPs: 0-3. For indels: unused (0).
    uint16_t refLen;       ///< SNP=1, ins=0, del=length.
    uint16_t altLen;       ///< SNP=1, ins=length, del=0.

    bool operator<(const KmVarKey& o) const {
        // Sort by sort_pos (indels use pos-1 like pgphase), then type, then refLen, then altLen, then altBase.
        uint32_t sp1 = (type == KmVarType::Snp) ? pos : (pos > 0 ? pos - 1 : 0);
        uint32_t sp2 = (o.type == KmVarType::Snp) ? o.pos : (o.pos > 0 ? o.pos - 1 : 0);
        if (sp1 != sp2) return sp1 < sp2;
        if (type != o.type) return type < o.type;
        if (refLen != o.refLen) return refLen < o.refLen;
        if (altLen != o.altLen) return altLen < o.altLen;
        return altBase < o.altBase;
    }
    bool operator==(const KmVarKey& o) const {
        return pos == o.pos && type == o.type && altBase == o.altBase
            && refLen == o.refLen && altLen == o.altLen;
    }
};

// ============================================================================
// Candidate variant site (mirrors pgphase CandidateVariant)
// ============================================================================

struct KmCandidate {
    KmVarKey key;

    // Allele counts.
    int totalCov = 0;
    int refCov = 0;
    int altCov = 0;
    int lowQualCov = 0;
    int fwdRef = 0, revRef = 0;
    int fwdAlt = 0, revAlt = 0;

    /// Per-allele coverage. Index 0 = ref, 1 = primary alt, 2+ = other alts.
    vector<int> alleCovs;
    int nUniqAlles = 2;

    double alleleFraction = 0.0;
    KmVariantCategory category = KmVariantCategory::LowCoverage;
    uint32_t categoryFlag = 0;
    bool isHomopolymerIndel = false;

    // K-means phasing state.
    std::array<vector<int>, 3> hapAlleProfile;
    std::array<int, 3> hapConsAlle = {-1, -1, -1};
    int32_t phaseSet = 0;
    int hapAlt = 0;
    int hapRef = 0;
};

// ============================================================================
// Per-overlap allele profile (mirrors pgphase ReadVariantProfile)
// ============================================================================

struct KmOverlapProfile {
    uint32_t overlapIdx = 0;
    int startVarIdx = -1;
    int endVarIdx = -1;
    /// Allele at each candidate in [startVarIdx, endVarIdx]:
    ///   0 = ref, 1 = primary alt, -1 = no observation, -2 = low-quality alt.
    vector<int> alleles;
};

// ============================================================================
// Noisy region interval
// ============================================================================

struct KmNoisyRegion {
    uint32_t start;
    uint32_t end;
    int label = 0;
    bool done = false;
};

// ============================================================================
// Per-overlap phasing state
// ============================================================================

struct KmOverlap {
    uint32_t alignmentId;
    uint32_t targetReadId;
    uint32_t qs, qe;
    uint32_t ts, te;
    uint32_t cigarOffset;
    uint32_t cigarTokenCount;
    uint8_t  isRev;
    uint8_t  queryIsRead0;
    int      hap = 0;
    int32_t  phaseSet = -1;
    uint8_t  strong = 0;
};

// ============================================================================
// Thread-local scratchpad
// ============================================================================

struct KmScratchpad {
    vector<KmOverlap> overlaps;
    vector<uint8_t> backboneBases;

    // Per-overlap parsed digars. Flat buffer with per-overlap ranges.
    vector<KmDigarOp> digars;
    vector<uint32_t> digarBegin; ///< digars[digarBegin[oi]..digarEnd[oi])
    vector<uint32_t> digarEnd;

    // Per-overlap noisy regions (pgphase NoisyRegionBuilder output).
    // Flat buffer with per-overlap ranges, like digars.
    vector<KmNoisyRegion> overlapNoisyRegions;
    vector<uint32_t> overlapNoisyBegin;
    vector<uint32_t> overlapNoisyEnd;

    // Low-complexity intervals on backbone (sdust output, sorted).
    vector<std::pair<uint32_t,uint32_t>> lowComplexity;

    // Chunk-level noisy regions (merged from per-overlap, after pre_process_noisy_regs).
    vector<KmNoisyRegion> noisyRegions;

    // Candidate table (sorted by KmVarKey).
    vector<KmCandidate> candidates;

    // Per-overlap allele profiles.
    vector<KmOverlapProfile> overlapProfiles;

    // K-means pivot selection scratch.
    vector<uint32_t> validVarIdx;

    void clear() {
        overlaps.clear();
        backboneBases.clear();
        digars.clear();
        digarBegin.clear();
        digarEnd.clear();
        overlapNoisyRegions.clear();
        overlapNoisyBegin.clear();
        overlapNoisyEnd.clear();
        lowComplexity.clear();
        noisyRegions.clear();
        candidates.clear();
        overlapProfiles.clear();
        validVarIdx.clear();
    }
};

// ============================================================================
// Thresholds
// ============================================================================

struct KmPhasingOptions {
    uint32_t minDepth = 5;
    uint32_t minAltDepth = 2;
    double   minAf = 0.20;
    double   maxAf = 0.80;
    bool     isOnt = false;           ///< ONT mode: enables Fisher exact strand bias test.
    double   strandBiasPval = 0.01;   ///< p-value threshold for ONT strand bias (pgphase default).
    uint32_t noisyRegMaxXgaps = 5;    ///< Max indel span for repeat check (pgphase kDefaultNoisyRegMaxXgaps).
    uint32_t noisyRegMergeDis = 500;  ///< Merge distance for noisy region intervals.
    uint32_t noisyRegFlankLen = 10;   ///< Flank extension for noisy regions (pgphase kNoisyRegFlankLen).
    // minNoisyRegTotalDepth removed — longcallD doesn't have this guard.
    uint32_t sdustThreshold = 5;      ///< SDUST complexity threshold (pgphase kSdustThreshold).
    uint32_t sdustWindow = 20;        ///< SDUST window size (pgphase kSdustWindow).
    uint32_t maxKmeansIter = 10;
    uint32_t minSpanningReads = 2;
    int      minSvLen = 30;           ///< Insertions >= this length use fuzzy length-ratio collapsing.
};

// ============================================================================
// K-means runner (defined in AssemblerPhasingKmeans.cpp)
// ============================================================================

/// Run k-means phasing on candidates matching the given category flags.
void kmRunKmeans(KmScratchpad& scratch, const KmPhasingOptions& opts, uint32_t flags);

} // namespace dinara

#endif // DINARA_PHASING_KMEANS_TYPES_HPP
