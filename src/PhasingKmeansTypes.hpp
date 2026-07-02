#ifndef DINARA_PHASING_KMEANS_TYPES_HPP
#define DINARA_PHASING_KMEANS_TYPES_HPP

/// @file PhasingKmeansTypes.hpp
/// @brief Data structures for k-means overlap phasing (pgphase/longcallD-style).

#include "ReadId.hpp"
#include "cstdint.hpp"
#include "vector.hpp"
#include <array>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>

namespace dinara {

class Assembler;

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
    std::string altSeq;    ///< For insertions: inserted bases (ACGT chars). Empty otherwise.

    /// Sort key: pos first, then type, then altBase, then len, then altSeq.
    bool operator<(const KmDigarOp& o) const {
        if (pos != o.pos) return pos < o.pos;
        if (type != o.type) return type < o.type;
        if (altBase != o.altBase) return altBase < o.altBase;
        if (len != o.len) return len < o.len;
        return altSeq < o.altSeq;
    }
    bool operator==(const KmDigarOp& o) const {
        return pos == o.pos && type == o.type && altBase == o.altBase && len == o.len && altSeq == o.altSeq;
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
    std::string altSeq;    ///< Insertion sequence (ACGT chars). Empty for SNPs/dels.

    bool operator<(const KmVarKey& o) const {
        // Sort by sort_pos (indels use pos-1 like pgphase), then type, then
        // refLen, then altLen, then altBase, then altSeq.
        // Matches pgphase exact_comp_var_site ordering.
        uint32_t sp1 = (type == KmVarType::Snp) ? pos : (pos > 0 ? pos - 1 : 0);
        uint32_t sp2 = (o.type == KmVarType::Snp) ? o.pos : (o.pos > 0 ? o.pos - 1 : 0);
        if (sp1 != sp2) return sp1 < sp2;
        if (type != o.type) return type < o.type;
        if (refLen != o.refLen) return refLen < o.refLen;
        if (altLen != o.altLen) return altLen < o.altLen;
        if (altBase != o.altBase) return altBase < o.altBase;
        return altSeq < o.altSeq;
    }
    bool operator==(const KmVarKey& o) const {
        return pos == o.pos && type == o.type && altBase == o.altBase
            && refLen == o.refLen && altLen == o.altLen
            && altSeq == o.altSeq;
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
    uint32_t minDepth = 6;
    uint32_t minAltDepth = 3;
    double   minAf = 0.20;
    double   maxAf = 0.8;
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
// Fisher exact test (pgphase / longcallD fisher_exact).
// Used for strand bias detection in both k-means and MSA phasing.
// ============================================================================

inline double kmLogHypergeom(int a, int b, int c, int d)
{
    const int n1 = a + b;
    const int n2 = c + d;
    const int m1 = a + c;
    const int m2 = b + d;
    const int N  = n1 + n2;
    if (N <= 0) return -std::numeric_limits<double>::infinity();
    if (n1 > n2) return kmLogHypergeom(c, d, a, b);
    if (m1 > m2) return kmLogHypergeom(b, a, d, c);
    return std::lgamma(double(n1 + 1)) + std::lgamma(double(n2 + 1)) +
           std::lgamma(double(m1 + 1)) + std::lgamma(double(m2 + 1)) -
          (std::lgamma(double(a  + 1)) + std::lgamma(double(b  + 1)) +
           std::lgamma(double(c  + 1)) + std::lgamma(double(d  + 1)) +
           std::lgamma(double(N  + 1)));
}

inline double kmFisherExactTwoTail(int a, int b, int c, int d)
{
    if (a + b + c + d <= 0) return 1.0;
    const double p_observed = std::exp(kmLogHypergeom(a, b, c, d));
    double total_p = 0.0;
    int min_a = (0 > (a + c) - (b + d)) ? 0 : (a + c) - (b + d);
    const int max_a = (a + b) < (a + c) ? (a + b) : (a + c);
    const int denom = a + b + c + d;
    const int mode_a = denom > 0
        ? int((double(a + b) * double(a + c)) / double(denom))
        : 0;
    for (int delta = 0; delta <= max_a - min_a; ++delta) {
        int cur_a = mode_a + delta;
        if (cur_a <= max_a) {
            int cur_b = (a + b) - cur_a;
            int cur_c = (a + c) - cur_a;
            int cur_d = (b + d) - cur_b;
            if (cur_b >= 0 && cur_c >= 0 && cur_d >= 0) {
                const double p = std::exp(kmLogHypergeom(cur_a, cur_b, cur_c, cur_d));
                if (p <= p_observed + DBL_EPSILON) total_p += p;
            }
        }
        if (delta > 0) {
            cur_a = mode_a - delta;
            if (cur_a >= min_a) {
                int cur_b = (a + b) - cur_a;
                int cur_c = (a + c) - cur_a;
                int cur_d = (b + d) - cur_b;
                if (cur_b >= 0 && cur_c >= 0 && cur_d >= 0) {
                    const double p = std::exp(kmLogHypergeom(cur_a, cur_b, cur_c, cur_d));
                    if (p <= p_observed + DBL_EPSILON) total_p += p;
                }
            }
        }
    }
    return total_p;
}

// ============================================================================
// Repeat / homopolymer detection (pgphase var_is_homopolymer_pg,
// var_is_repeat_region_pg). Used by both k-means and MSA phasing.
// Backbone sequence is numeric: 0=A, 1=C, 2=G, 3=T.
// ============================================================================

/// Check if a variant sits in a homopolymer or short tandem repeat context.
/// Looks for a repeat unit of length 1..6 with ≥3 copies flanking the variant.
// Core repeat-context test, restricted to repeat unit lengths in
// [minUnitLen, maxUnitLen]. A variant is "in a repeat" if the reference flank
// (forward from endPos, or backward from startPos) consists of >=3 copies of a
// unit whose length is in the requested range. Callers use this to distinguish
// homopolymer context (unit length 1) from short-tandem-repeat context
// (unit length 2..6), which are otherwise identical tests over different ranges.
inline bool kmIsRepeatUnitRange(const uint8_t* seq, uint32_t seqLen,
                                const KmVarKey& key, int xid,
                                int minUnitLen, int maxUnitLen)
{
    if (seqLen == 0) return false;
    const int64_t sn = int64_t(seqLen);
    int64_t startPos, endPos;
    if (key.type == KmVarType::Snp) {
        startPos = int64_t(key.pos) - 1;
        endPos   = int64_t(key.pos) + 1;
    } else if (key.type == KmVarType::Insertion) {
        if (int(key.altLen) > xid) return false;
        startPos = int64_t(key.pos) - 1;
        endPos   = int64_t(key.pos);
    } else { // Deletion
        if (int(key.refLen) > xid) return false;
        startPos = int64_t(key.pos) + int64_t(key.refLen) - 1;
        endPos   = int64_t(key.pos);
    }
    constexpr int nCheckCopyNum = 3;
    auto safeBase = [&](int64_t p) -> uint8_t {
        if (p < 0 || p >= sn) return 4;
        return seq[p];
    };
    // Check forward from endPos.
    uint8_t refBases[6];
    for (int i = 0; i < 6; i++) refBases[i] = safeBase(endPos + i);
    for (int r = minUnitLen; r <= maxUnitLen; r++) {
        bool isHp = true;
        for (int i = 1; i < nCheckCopyNum && isHp; i++)
            for (int j = 0; j < r && isHp; j++)
                if (safeBase(endPos + i * r + j) != refBases[j]) isHp = false;
        if (isHp) return true;
    }
    // Check backward from startPos.
    for (int i = 0; i < 6; i++) refBases[i] = safeBase(startPos - i);
    for (int r = minUnitLen; r <= maxUnitLen; r++) {
        bool isHp = true;
        for (int i = 1; i < nCheckCopyNum && isHp; i++)
            for (int j = 0; j < r && isHp; j++)
                if (safeBase(startPos - i * r - j) != refBases[j]) isHp = false;
        if (isHp) return true;
    }
    return false;
}

// Repeat context over unit lengths 1..6 (homopolymers AND short tandem
// repeats). Kept for existing callers; equivalent to the union of the
// homopolymer (unit 1) and STR (unit 2..6) tests.
inline bool kmIsHomopolymer(const uint8_t* seq, uint32_t seqLen,
                            const KmVarKey& key, int xid)
{
    return kmIsRepeatUnitRange(seq, seqLen, key, xid, 1, 6);
}

/// Check if the deleted/inserted motif is a tandem repeat of the flanking reference.
inline bool kmIsRepeatRegion(const uint8_t* seq, uint32_t seqLen,
                             const KmVarKey& key, int xid)
{
    if (seqLen == 0) return false;
    const int64_t sn = int64_t(seqLen);
    if (key.type == KmVarType::Deletion) {
        const int delLen = int(key.refLen);
        if (delLen > xid) return false;
        const int len = delLen * 3;
        const int64_t pos = int64_t(key.pos);
        if (pos < 0 || pos + delLen + len > sn) return false;
        for (int i = 0; i < len; i++)
            if (seq[pos + i] != seq[pos + delLen + i]) return false;
        return true;
    }
    if (key.type == KmVarType::Insertion) {
        if (key.altSeq.empty()) return false;
        const int insLen = int(key.altSeq.size());
        if (insLen > xid) return false;
        const int len = insLen * 3;
        const int64_t pos = int64_t(key.pos);
        if (pos < 0 || pos + len > sn) return false;
        std::vector<uint8_t> alt_b(static_cast<size_t>(len));
        for (int j = 0; j < len; j++)
            alt_b[j] = seq[pos + j];
        for (int j = insLen; j < len; j++)
            alt_b[j] = alt_b[j - insLen];
        for (int j = 0; j < insLen; j++) {
            const char c = key.altSeq[j];
            uint8_t b = 4;
            if (c == 'A' || c == 'a') b = 0;
            else if (c == 'C' || c == 'c') b = 1;
            else if (c == 'G' || c == 'g') b = 2;
            else if (c == 'T' || c == 't') b = 3;
            alt_b[j] = b;
        }
        for (int j = 0; j < len; j++)
            if (seq[pos + j] != alt_b[j]) return false;
        return true;
    }
    return false;
}

// ============================================================================
// K-means runner (defined in AssemblerPhasingKmeans.cpp)
// ============================================================================

/// Run k-means phasing on candidates matching the given category flags.
void kmRunKmeans(KmScratchpad& scratch, const KmPhasingOptions& opts, uint32_t flags);

/// Write phasing results from scratchpad to AlignmentData.
void kmWriteResults(Assembler& assembler, ReadId backboneReadId, const KmScratchpad& scratch);

/// Cis refinement: detect cisDifferentCopy within the cis set.
void kmRefineCis(Assembler& assembler, ReadId backboneReadId,
    KmScratchpad& scratch, const KmPhasingOptions& opts,
    uint32_t bbLen, bool debug);

} // namespace dinara

#endif // DINARA_PHASING_KMEANS_TYPES_HPP
