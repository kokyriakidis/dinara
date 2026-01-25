/*
Dinara hifiasm-parity error-correction (Parity EC).

This translation unit implements a per-read pipeline modeled after hifiasm’s phasing-based
error-correction logic (Correct.cpp, rphase_hc + gen_rphase_dp + generate_haplotypes_*).

The implementation here intentionally avoids base-quality (QV) logic. “Parity” in this file
means we match hifiasm’s decision structure and state transitions using Dinara’s data model:
alignment-derived evidence streams instead of hifiasm’s overlap_region_alloc and internal
haplotype_evdience lists.

Per query read, the pipeline is:

  (A) Candidate gathering
      - Build a local list of overlaps (candidates[]) for a single query read.
      - Candidates are stored in the query’s coordinate system.
      - Each candidate carries a hifiasm-like state flag (is_match) that is updated across
        stages:
           is_match == 1 : kept/cis
           is_match == 2 : filtered/trans (for HiFi parity we never flip 2 back to 1)

  (B) SNP evidence collection (detectHetSites)
      - Read mismatch evidence from AlignedEvidenceStore for each overlap.
      - Project mismatch positions to query coordinates, then group by query position.
      - Compute coverage at each potential SNP site using a sweep-line style pass over
        candidate intervals.
      - Model each (site, alt-base) as a separate SnpStats row, and model per-overlap allele
        observations as HaplotypeEvidence entries.
      - Tag each emitted SNP row with a query-side homopolymer flag (if_is_homopolymer_strict),
        but do not drop rows solely due to query homopolymers at this stage.
      - Build “coverage holes” from the indel stream (query-coordinate gaps where the overlap
        has no aligned base) and subtract those overlaps from per-site coverage/ref support.

  (C) DP chaining (gen_rphase_dp)
      - Build compact bit-matrices per SNP row to support very fast linkage tests between sites
        (Ref bits, Alt bits, and AnyOther bits).
      - Run a dynamic program that chains sites when there is both ref-ref and alt-alt support
        across overlaps and no conflicting “other allele” evidence.
      - Apply hifiasm-style cc filtering based on coveragePeak and n_hap=2 (cut_rate=0.7,
        cut_bd=6).
      - Compact SNP rows down to the DP-retained set and reset score semantics to match hifiasm:
        DP selects rows, but downstream trans-closure starts with score = -1 again.

  (D) HiFi trans-closure (generate_haplotypes_naive_HiFi)
      - Seed trans overlaps using validated SNP rows, then apply “not real allele” decrements:
        trans overlaps reduce ref support occ_0 at sites they cover but do not carry as alt.
      - Optionally run multi_check recovery for dense weak-variant patterns.
      - Mark overlaps as trans (is_match=2), then reflect the decision into AlignmentData flags.

  (E) SV recovery (detectSVSites + generate_haplotypes_sv)
      - Cluster large indel events from the indel evidence stream.
      - Perform an analogous overlap marking stage for SV sites using the same trans-closure
        structure as hifiasm’s generate_haplotypes_sv.

The design is performance-oriented:
  - HifiasmECScratchPad holds all large vectors and is reused per thread to avoid per-read
    allocations.
  - Evidence is processed in query-sorted order to enable linear scans and lower_bound jumps.
  - Most “per-site” expensive work is done via bitsets in DP, not via nested loops.

Testing:
  - When compiled with -DDINARA_TESTING (integration tests), we expose small hooks/counters to
    validate DP behavior without changing production behavior.
*/

#include "Assembler.hpp"
#include "AlignedEvidenceStore.hpp"
#include "LongBaseSequence.hpp"
#include "Reads.hpp"
#include "chrono.hpp"
#include "timestamp.hpp"
#include <algorithm>
#include <atomic>
#include <cstring>
#include <immintrin.h>
#include <iostream>
#include <limits>
#include <mutex>
#include <vector>

using namespace dinara;
using namespace std;

inline int64_t comput_sc_rphase_strict(
    const uint64_t* flatBits,
    const uint64_t* flatAnyBits,
    size_t siteI, size_t siteJ,
    size_t nWords
);

#ifdef DINARA_TESTING
namespace {
    std::atomic<uint64_t> gDpSameSiteComparisons{0};
}

namespace dinara::hifiasmEcTestHooks {
    int64_t scoreLinkStrict(
        const std::vector<uint64_t>& flatBits,
        const std::vector<uint64_t>& flatAnyBits,
        size_t siteI,
        size_t siteJ,
        size_t nWords)
    {
        return comput_sc_rphase_strict(flatBits.data(), flatAnyBits.data(), siteI, siteJ, nWords);
    }

    void resetDpSameSiteComparisons()
    {
        gDpSameSiteComparisons.store(0, std::memory_order_relaxed);
    }

    uint64_t getDpSameSiteComparisons()
    {
        return gDpSameSiteComparisons.load(std::memory_order_relaxed);
    }
}
#endif

/*
Candidate overlap descriptor for a single query read.

All coordinates are stored in the query’s coordinate system:
  - qs/qe: query interval covered by this overlap (half-open).
  - ts/te: target interval in the other read (half-open).

The flags mirror hifiasm’s overlap_region bookkeeping:
  - is_match == 1 means the overlap is currently kept (cis).
  - is_match == 2 means the overlap is marked trans and should be filtered for this query.
  - strong mirrors hifiasm’s “strong” flag (kept for parity naming/telemetry).
*/
struct CandidateEC {
    uint32_t alignmentId;
    uint64_t qs, qe, ts, te;
    uint32_t targetId;
    bool isRev;
    uint8_t is_match = 1;
    uint8_t strong = 0;
};

/*
Per-allele statistics for a single SNP row.

Hifiasm models multi-allelic sites as multiple “rows” sharing the same genomic position,
one row per alternative base. We match that shape because downstream DP and trans-closure
logic is written in “row space”, not “site space”.

The occ_0/occ_1 counters are defined on the query read:
  - occ_1 is the number of overlaps supporting this specific alternative base at this site.
  - occ_0 is the number of overlaps supporting the query base (ref) at this site, plus 1 for
    the query read itself.
*/
struct SnpStats {
    uint32_t site;
    uint32_t occ_0;
    uint32_t occ_1;
    uint32_t fwd_ref_cov;
    char refBase; 
    char altBase; 
    uint8_t is_homopolymer;
    int score;
    int dpScore;
    
    SnpStats() {
        site = (uint32_t)-1;
        occ_0 = occ_1 = 0;
        fwd_ref_cov = 0;
        refBase = altBase = 0;
        is_homopolymer = 0;
        score = 0;
        dpScore = 0;
    }
};

/*
Per-overlap evidence entry for a single site.

This is the compact representation used by downstream phases:
  - overlapID indexes into the per-query candidates[] array.
  - site is the query-coordinate position of the event.
  - overlapSite is the SNP-row index (row in snpStats/svStats) for quick lookups.
  - type encodes the evidence kind:
      0 -> ref/supports query base at this site
      1 -> alt/mismatch base at this site
      2 -> SV/indel support (used in SV stage)
  - misBase stores the alt base for type==1.
  - hp is a homopolymer-suspect flag (parity with hifiasm’s hh_hp tagging).
*/
struct HaplotypeEvidence {
    uint32_t overlapID;
    uint32_t site;
    uint32_t overlapSite;
    uint8_t type;
    uint8_t misBase;
    bool hp;

    bool operator<(const HaplotypeEvidence& other) const {
        if(site != other.site) return site < other.site;
        return overlapID < other.overlapID;
    }
};

/*
Sweep-line event for coverage computation.

We compute per-site coverage across potentially many overlaps by turning overlap intervals into
start/end events in site-index space. Sorting these events lets us maintain an active set size
while iterating sites in increasing order.
*/
struct SweepEvent {
    size_t siteIdx;
    uint32_t candIdx;
    bool isEnd;
    bool operator<(const SweepEvent& other) const {
        if(siteIdx != other.siteIdx) return siteIdx < other.siteIdx;
        return isEnd < other.isEnd;
    }
};

struct RawSV {
    uint32_t overlapID;
    uint32_t site; 
    int64_t size; 
};

struct RphaseDpTiming {
    double buildAltAnyBits = 0.;
    double buildRefBits = 0.;
    double buildHpBits = 0.;
    double transitions = 0.;
    double extractPaths = 0.;
};

/*
Thread-local scratchpad for the parity EC pipeline.

The pipeline creates many temporary arrays (sites, events, bitsets, per-overlap groupings).
Reallocating them per read is expensive, so each worker thread keeps one scratchpad and reuses
storage across reads.

Most vectors are multi-purpose buffers that are “reinterpreted” at different phases. This is
intentional to reduce peak memory and allocation churn.
*/
struct HifiasmECScratchPad {
    struct GapInterval {
        uint32_t begin;
        uint32_t end;
    };

    vector<CandidateEC> candidates;
    vector<SnpStats> snpStats;
    vector<HaplotypeEvidence> hapEvidence;
    vector<SnpStats> svStats;
    vector<HaplotypeEvidence> svEvidence;

    vector<uint32_t> uniqueSites;
    vector<int32_t> diffTotal;
    vector<int32_t> diffFwd;
    vector<uint32_t> siteTotalCov;
    vector<uint32_t> siteFwdCov;

    vector<uint32_t> insertionOffsets;
    vector<GapInterval> insertionIntervals;
    vector<uint32_t> insertionBaseCount;

    vector<int> validIndices;
    vector<int64_t> f;
    vector<int> p;
    vector<int> indexMap;
    vector<uint64_t> flatBits;
    vector<uint64_t> flatAnyBits;
    vector<uint64_t> flatHpBits;
    vector<SweepEvent> events;
    vector<SweepEvent> gapEvents;
    vector<uint64_t> active;
    vector<uint64_t> gapped;

    vector<RawSV> rawSVs;
    vector<size_t> svIndices;
    vector<uint8_t> unpackedRead;
    vector<uint8_t> covered;
    vector<uint32_t> path;
    vector<uint64_t> supportBits;
    vector<uint64_t> conflictBits;

    void clear() {
        candidates.clear();
        snpStats.clear();
        hapEvidence.clear();
        svStats.clear();
        svEvidence.clear();
        uniqueSites.clear();
        diffTotal.clear();
        diffFwd.clear();
        siteTotalCov.clear();
        siteFwdCov.clear();
        insertionOffsets.clear();
        insertionIntervals.clear();
        insertionBaseCount.clear();
        validIndices.clear();
        f.clear();
        p.clear();
        indexMap.clear();
        flatBits.clear();
        flatAnyBits.clear();
        flatHpBits.clear();
        events.clear();
        gapEvents.clear();
        active.clear();
        gapped.clear();
        rawSVs.clear();
        svIndices.clear();
        unpackedRead.clear();
        covered.clear();
        path.clear();
        supportBits.clear();
        conflictBits.clear();
    }
};

/*
Return the base character at position pos in a read, with bounds checking.

This helper is used by homopolymer detection routines. Hot paths in the parity-EC pipeline use
pre-decoded 2-bit bases (see unpack2BitBases) to avoid repeated read accessor overhead.
*/
inline char getBase(const LongBaseSequenceView& read, uint64_t pos) {
    if(pos >= read.baseCount) return 'N';
    return read[pos].character(); 
}

/*
Convert an ASCII base into an integer in the range [0..4].

0=A, 1=C, 2=G, 3=T, 4=other. The value 4 is used as a safe “unknown” bucket.
*/
inline int base2int(char c) {
    switch(c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

/*
Decode Dinara’s packed LongBaseSequence into a flat 2-bit-per-base array.

The output array stores A=0,C=1,G=2,T=3. This is used to make homopolymer checks fast and
predictable (tight loops over a contiguous byte array).
*/
static void unpack2BitBases(const LongBaseSequenceView& seq, vector<uint8_t>& out)
{
    out.resize(seq.baseCount);
    const uint64_t blocks = (seq.baseCount + 63ULL) / 64ULL;
    for (uint64_t block = 0; block < blocks; ++block) {
        const uint64_t baseBegin = block * 64ULL;
        const uint64_t baseEnd = std::min(baseBegin + uint64_t(64), seq.baseCount);
        const uint64_t word0 = seq.begin[2ULL * block];
        const uint64_t word1 = seq.begin[2ULL * block + 1ULL];
        for (uint64_t j = 0; baseBegin + j < baseEnd; ++j) {
            const uint64_t bitIndex = 63ULL - j;
            const uint8_t bit0 = (uint8_t)((word0 >> bitIndex) & 1ULL);
            const uint8_t bit1 = (uint8_t)((word1 >> bitIndex) & 1ULL);
            out[baseBegin + j] = (uint8_t)((bit1 << 1U) | bit0);
        }
    }
}

/*
Hifiasm-style low-complexity/homopolymer mask (hpc_mask_ff).

This function checks whether the sequence around pos exhibits a strong periodic run pattern.
It is a direct adaptation of hifiasm’s masking logic and is used to tag “HP-suspect” evidence.

Implementation notes:
  - The scan considers multiple small periods r (1..hpc_min) and checks whether the maximal run
    length for that period exceeds a threshold scaled by r.
  - The check is performed in multiple variants (“including pos” vs “excluding pos”, forward and
    reverse) to mirror hifiasm’s behavior.
  - SeqType is either a raw byte pointer (unpacked query bases) or a vector-like object.
*/
template<typename SeqType>
static bool hpc_mask_ff(const SeqType& seq, int64_t len, uint32_t pos) {
    const int64_t hpc_len = 12;
    const int64_t hpc_min = 4;
    const int64_t hpc_cut = 2;
    
    auto getBase = [&](int64_t idx) -> uint8_t {
        if constexpr (std::is_pointer_v<SeqType>) return seq[idx];
        else return seq[idx].value;
    };

    const int64_t p = (int64_t)pos;
    const int64_t e = (p + hpc_len <= len) ? (p + hpc_len) : len;
    const int64_t s = (p >= hpc_len) ? (p - hpc_len) : 0;

    for (int64_t r = 1; r <= hpc_min; r++) {
        const int64_t rc = r * hpc_cut;
        
        int64_t k;
        for (k = p + r; (k < e) && ((k-r) >= s) && (getBase(k) == getBase(k-r)); k++);
        int64_t ze = k;
        for (k = p - 1; (k >= s) && ((k+r) < e) && (getBase(k) == getBase(k+r)); k--);
        int64_t zs = k + 1;
        if ((ze - zs > r) && (ze - zs >= rc)) return true;
        
        for (k = p + r + 1; (k < e) && ((k-r) >= s) && (getBase(k) == getBase(k-r)); k++); 
        ze = k; zs = p + 1;
        if ((ze - zs > r) && (ze - zs >= rc)) return true;
        
        for (k = p - r; (k >= s) && ((k+r) < e) && (getBase(k) == getBase(k+r)); k--);
        zs = k + 1;
        for (k = p + 1; (k < e) && ((k-r) >= s) && (getBase(k) == getBase(k-r)); k++);
        ze = k;
        if ((ze - zs > r) && (ze - zs >= rc)) return true;

        for (k = p - r - 1; (k >= s) && ((k+r) < e) && (getBase(k) == getBase(k+r)); k--);
        zs = k + 1; ze = p;
        if ((ze - zs > r) && (ze - zs >= rc)) return true;
    }
    return false;
}

/*
Hifiasm’s strict query-side homopolymer predicate (if_is_homopolymer_strict).

This function only inspects the query sequence. It is used to annotate SNP rows with a flag
that is consulted later in the pipeline (it is not, by itself, a reason to drop a site during
SNP emission).
*/
template<typename SeqType>
static bool if_is_homopolymer_strict(const SeqType& seq, int64_t len, uint32_t pos)
{
    auto getBase = [&](int64_t idx) -> uint8_t {
        if constexpr (std::is_pointer_v<SeqType>) return seq[idx];
        else return seq[idx].value;
    };

    const int64_t threshold = 3;
    const int64_t site = (int64_t)pos;

    int64_t beg = site - threshold;
    if (beg < 0) beg = 0;

    int64_t end = site + threshold;
    if (end >= len) end = len - 1;

    const uint8_t center = getBase(site);

    bool fSet = false;
    uint8_t fHomCh = 0;
    int64_t fHomLen = 0;
    for (int64_t i = site + 1; i <= end; ++i) {
        const uint8_t b = getBase(i);
        if (!fSet) {
            fHomCh = b;
            fHomLen = 1;
            fSet = true;
        } else {
            if (b != fHomCh) break;
            fHomLen++;
        }
    }

    bool bSet = false;
    uint8_t bHomCh = 0;
    int64_t bHomLen = 0;
    for (int64_t i = site - 1; i >= beg; --i) {
        const uint8_t b = getBase(i);
        if (!bSet) {
            bHomCh = b;
            bHomLen = 1;
            bSet = true;
        } else {
            if (b != bHomCh) break;
            bHomLen++;
        }
    }

    if (fSet && fHomCh == center) {
        fHomLen++;
    } else if (bSet && bHomCh == center) {
        bHomLen++;
    }

    if (fHomLen >= threshold || bHomLen >= threshold) {
        return true;
    }

    if (fSet && bSet &&
        center == fHomCh &&
        bHomCh == fHomCh &&
        (fHomLen + bHomLen >= threshold)) {
        return true;
    }

    return false;
}



/*
detectHetSites: build SNP candidate rows and per-overlap mismatch evidence for one query read.

Conceptually this stage corresponds to the “site aggregation” and “push_info” steps in hifiasm:
we turn raw per-alignment mismatch evidence into:
  - a set of candidate SNP rows (snpStats), including multi-allelic sites (one row per alt base)
  - a per-overlap list of mismatches at those sites (hapEvidence)

We also compute per-overlap query-coordinate “hole” intervals from the indel evidence stream.
These holes represent query positions where the overlap has no aligned target base, so the overlap
must not contribute coverage or ref-support for SNP calling inside those intervals.
*/
template<typename AlignmentContainer>
static void detectHetSites(
    Assembler& assembler,
    const Reads& reads,
    uint64_t queryReadId,
    const AlignmentContainer& alignmentData,
    HifiasmECScratchPad& scratch
) {
    const LongBaseSequenceView queryRead = reads.getRead(ReadId(queryReadId));
    auto& candidates = scratch.candidates;
    auto& hapEvidence = scratch.hapEvidence;
    auto& snpStats = scratch.snpStats;
    auto& insertionOffsets = scratch.insertionOffsets;
    auto& insertionIntervals = scratch.insertionIntervals;
    auto& insertionBaseCount = scratch.insertionBaseCount;

    /*
    Step 1: decode the query read into 2-bit bases.

    Several downstream checks compare mismatch bases against the query base and run homopolymer
    predicates. Doing these operations over a flat byte array avoids repeated read-view access
    overhead in tight loops.
    */
    auto& unpacked = scratch.unpackedRead;
    unpack2BitBases(queryRead, unpacked);

    insertionIntervals.clear();
    insertionOffsets.assign(candidates.size() + 1, 0);
    insertionBaseCount.assign(candidates.size(), 0);

    /*
    Step 2: per-overlap evidence extraction.

    For each overlap candidate, we:
      - build query-coordinate “hole” intervals from the indel stream (deletions in the query
        coordinate sense: query bases without an aligned partner)
      - append mismatch evidence from the SNP stream within the overlap span

    We do not emit “ref” evidence entries here. Ref support is derived later from coverage bitsets:
    “covered and not (alt or anyOther)”.
    */
    for (size_t k = 0; k < candidates.size(); ++k) {
        insertionOffsets[k] = (uint32_t)insertionIntervals.size();
        insertionBaseCount[k] = 0;

        const auto& cand = candidates[k];
        const auto& ad = alignmentData[cand.alignmentId];
        const size_t evidenceId = ad.info.alignmentId;
        
        if (evidenceId == invalid<size_t>) continue;

        const uint32_t bd = 0; 
        const uint32_t c_qs = cand.qs + bd;
        const uint32_t c_qe = (cand.qe > bd) ? (cand.qe - bd) : cand.qs;
        if (c_qe < c_qs) continue;

        uint32_t insertedBases = 0;
        {
            span<const IndelEvidence> indels;
            if (ad.readIds[1] == queryReadId) {
                indels = assembler.alignedEvidenceStore.getIndels0(uint32_t(evidenceId));
            } else if (ad.readIds[0] == queryReadId) {
                indels = assembler.alignedEvidenceStore.getIndels1(uint32_t(evidenceId));
            }

            const size_t baseOffset = insertionIntervals.size();
            auto appendInterval = [&](uint32_t begin, uint32_t end) {
                if (begin >= end) return;
                if (insertionIntervals.size() > baseOffset) {
                    auto& last = insertionIntervals.back();
                    if (begin <= last.end) {
                        insertedBases -= (last.end - last.begin);
                        last.end = std::max(last.end, end);
                        insertedBases += (last.end - last.begin);
                        return;
                    }
                }
                insertionIntervals.push_back({begin, end});
                insertedBases += (end - begin);
            };

            if (!indels.empty() && c_qe > c_qs) {
                auto it = std::lower_bound(
                    indels.begin(), indels.end(), c_qs,
                    [](const IndelEvidence& e, uint32_t value) { return e.pos() < value; }
                );
                for (; it != indels.end(); ++it) {
                    const uint32_t pos = it->pos();
                    if (pos >= c_qe) break;
                    if (!it->isDeletion()) continue;

                    uint32_t begin = pos;
                    uint32_t end = pos + it->len();
                    if (end <= c_qs || begin >= c_qe) continue;
                    begin = std::max(begin, c_qs);
                    end = std::min(end, c_qe);
                    appendInterval(begin, end);
                }
            }
        }
        insertionBaseCount[k] = insertedBases;

        const uint32_t evidenceId32 = uint32_t(evidenceId);
        auto addEv = [&](uint32_t currentPos, uint8_t base) {
            if (currentPos < unpacked.size() && base == unpacked[currentPos]) {
                return;
            }
            HaplotypeEvidence ev;
            ev.overlapID = (uint32_t)k;
            ev.site = currentPos;
            ev.overlapSite = invalid<uint32_t>;
            ev.type = 1;
            ev.misBase = base;
            ev.hp = false;
            hapEvidence.push_back(ev);
        };

        if (ad.readIds[1] == queryReadId) {
            assembler.alignedEvidenceStore.forEachSnp0InRange(evidenceId32, c_qs, c_qe, addEv);
        } else if (ad.readIds[0] == queryReadId) {
            assembler.alignedEvidenceStore.forEachSnp1InRange(evidenceId32, c_qs, c_qe, addEv);
        } else {
            continue;
        }
    }
    insertionOffsets[candidates.size()] = (uint32_t)insertionIntervals.size();
    
    std::sort(hapEvidence.begin(), hapEvidence.end());
    if (hapEvidence.empty()) return;

    /*
    Step 3: build the sorted list of unique query positions that have at least one mismatch.

    This list is later used for:
      - sweep-line coverage computation (O(#candidates + #sites))
      - mapping “hole” intervals into site-index ranges via lower_bound
    */
    auto& uniqueSites = scratch.uniqueSites;
    uniqueSites.clear();
    for (size_t i = 0; i < hapEvidence.size(); ) {
        uniqueSites.push_back(hapEvidence[i].site);
        uint32_t current = hapEvidence[i].site;
        while (i < hapEvidence.size() && hapEvidence[i].site == current) i++;
    }

    /*
    Step 4: compute coverage at each candidate site.

    For each unique query position in uniqueSites we need:
      - total number of overlaps that cover that position (siteTotalCov)
      - number of those overlaps on the “forward” orientation in our per-query coordinate system
        (siteFwdCov), used for strand-bias checks

    Naively testing every overlap against every site would be O(#overlaps * #sites). Instead we:
      - convert each overlap into a [start,end) range in site-index space using lower_bound
      - add that range to a difference array, then prefix-sum once

    We subtract each overlap’s “hole” intervals so that positions with no aligned base do not
    contribute to coverage or to derived ref counts.
    */
    const size_t nSites = uniqueSites.size();
    auto& diffTotal = scratch.diffTotal;
    auto& diffFwd = scratch.diffFwd;
    diffTotal.assign(nSites + 1, 0);
    diffFwd.assign(nSites + 1, 0);
    const uint32_t bd = 0; 
    
    for (size_t candIdx = 0; candIdx < candidates.size(); ++candIdx) {
        const auto& cand = candidates[candIdx];
        const uint32_t s = cand.qs + bd;
        const uint32_t e = (cand.qe > bd) ? (cand.qe - bd) : cand.qs;
        if (e <= s) continue;

        auto itStart = std::lower_bound(uniqueSites.begin(), uniqueSites.end(), s);
        auto itEnd = std::lower_bound(uniqueSites.begin(), uniqueSites.end(), e);
        
        size_t idxS = std::distance(uniqueSites.begin(), itStart);
        size_t idxE = std::distance(uniqueSites.begin(), itEnd);

        if (idxS < idxE) {
            diffTotal[idxS]++;
            diffTotal[idxE]--;
            if (!cand.isRev) {
                diffFwd[idxS]++;
                diffFwd[idxE]--;
            }
        }

        if (candIdx + 1 < insertionOffsets.size()) {
            const uint32_t offBegin = insertionOffsets[candIdx];
            const uint32_t offEnd = insertionOffsets[candIdx + 1];
            for (uint32_t ii = offBegin; ii < offEnd; ++ii) {
                const auto& interval = insertionIntervals[ii];
                auto itHoleS = std::lower_bound(uniqueSites.begin(), uniqueSites.end(), interval.begin);
                auto itHoleE = std::lower_bound(uniqueSites.begin(), uniqueSites.end(), interval.end);
                const size_t hS = std::distance(uniqueSites.begin(), itHoleS);
                const size_t hE = std::distance(uniqueSites.begin(), itHoleE);
                if (hS < hE) {
                    diffTotal[hS]--;
                    diffTotal[hE]++;
                    if (!cand.isRev) {
                        diffFwd[hS]--;
                        diffFwd[hE]++;
                    }
                }
            }
        }
    }

    auto& siteTotalCov = scratch.siteTotalCov;
    auto& siteFwdCov = scratch.siteFwdCov;
    siteTotalCov.resize(nSites);
    siteFwdCov.resize(nSites);
    int32_t curTotal = 0, curFwd = 0;
    for (size_t k = 0; k < nSites; ++k) {
        curTotal += diffTotal[k];
        curFwd += diffFwd[k];
        siteTotalCov[k] = (uint32_t)curTotal;
        siteFwdCov[k] = (uint32_t)curFwd;
    }

    /*
    Step 5: turn per-site mismatch counts + per-site coverage into SNP rows.

    For each query position in uniqueSites we:
      - count mismatch support per base (and per strand)
      - combine with coverage to estimate ref support at the site
      - apply minimal evidence and strand-bias filters
      - emit one SnpStats row per alternative base with enough support (multi-allelic modeled as
        multiple rows sharing the same site)
      - annotate the row with the query-side strict homopolymer flag
      - map each mismatch evidence entry at this site to the correct emitted row (overlapSite)
    */
    for (size_t siteIdx = 0, evidenceIdx = 0; siteIdx < nSites; ++siteIdx) {
        const uint32_t site = uniqueSites[siteIdx];
        const uint32_t startIdx = (uint32_t)evidenceIdx;
        
        /*
        Gather mismatch support at this site.

        hapEvidence is sorted by site, so we can consume a contiguous block of entries for the
        current site in one pass.
        */
        uint32_t misCountPerBase[4] = {0, 0, 0, 0};
        uint32_t misCountFwdPerBase[4] = {0, 0, 0, 0};
        uint32_t totalMisCount = 0;
        uint32_t totalAltReadsFwd = 0;
        while (evidenceIdx < hapEvidence.size() && hapEvidence[evidenceIdx].site == site) {
            const auto& ev = hapEvidence[evidenceIdx];
            const uint8_t base = ev.misBase;
            if (base < 4) {
                misCountPerBase[base]++;
                totalMisCount++;
                const uint32_t ovId = ev.overlapID;
                if (ovId < candidates.size() && !candidates[ovId].isRev) {
                    misCountFwdPerBase[base]++;
                    totalAltReadsFwd++;
                }
            }
            evidenceIdx++;
        }

        /*
        Hifiasm’s push_info only creates SNP rows for alternative alleles with at least 2 supports.
        We apply the same rule here: if no alt base reaches 2 overlaps, the site is ignored.
        */
        uint32_t maxMisCount = 0;
        for (uint8_t b = 0; b < 4; ++b) {
            maxMisCount = std::max(maxMisCount, misCountPerBase[b]);
        }
        if (maxMisCount < 2) {
            continue;
        }
        
        /*
        Combine mismatch counts with per-site coverage.

        totalCov and totalReadsFwd were computed in Step 4 with “hole” subtraction. refCov is
        derived as totalCov - totalMisCount (clamped at 0 for safety).
        */
        const uint32_t totalCov = siteTotalCov[siteIdx];
        const uint32_t totalReadsFwd = siteFwdCov[siteIdx];
        const uint32_t refCov = (totalCov >= totalMisCount) ? (totalCov - totalMisCount) : 0;

        /*
        Minimum total coverage filter.

        This is a pragmatic guard used in parity EC to avoid spending time on extremely low-support
        sites that are dominated by sequencing noise. It also matches the behavior expected by the
        downstream overlap marking logic.
        */
        if (totalCov < 5) continue; 
        
        /*
        Strand-bias filter on ref support.

        Hifiasm’s is_st_bs rejects extreme cases where essentially all ref support is on one strand
        when refCov is sufficiently large. We replicate that behavior using integer math.
        */
        if (refCov > 2) {
            uint32_t refReadsFwd = (totalReadsFwd >= totalAltReadsFwd) ? (totalReadsFwd - totalAltReadsFwd) : 0;
            if ((refReadsFwd + 2 >= refCov) && (uint64_t(refReadsFwd) * 100ULL >= uint64_t(refCov) * 95ULL)) continue;
            if ((refReadsFwd <= 2) && (uint64_t(refReadsFwd) * 100ULL <= uint64_t(refCov) * 5ULL)) continue;
        }

        /*
        Decide which alternative alleles to emit.

        We emit one SNP row per alt base that:
          - has at least 2 supporting overlaps at this site
          - passes an alt-strand-bias filter (again modeled after hifiasm push_info behavior)
        */
        uint32_t rowIndexForBase[4] = {
            invalid<uint32_t>, invalid<uint32_t>, invalid<uint32_t>, invalid<uint32_t>
        };
        uint32_t validAltCount = 0;

        for (uint8_t altBaseIdx = 0; altBaseIdx < 4; ++altBaseIdx) {
            const uint32_t misCount = misCountPerBase[altBaseIdx];
            if (misCount < 2) continue;

            const uint32_t altReadsFwd = misCountFwdPerBase[altBaseIdx];

            if (misCount > 2) { 
                if ((altReadsFwd + 2 >= misCount) && (uint64_t(altReadsFwd) * 100ULL >= uint64_t(misCount) * 95ULL)) continue;
                if ((altReadsFwd <= 2) && (uint64_t(altReadsFwd) * 100ULL <= uint64_t(misCount) * 5ULL)) continue;
            }
            rowIndexForBase[altBaseIdx] = 0;
            validAltCount++;
        }

        if (validAltCount == 0) {
            continue;
        }

        /*
        Record query-side strict homopolymer status.

        This flag is part of the per-row metadata and is consulted later; we do not drop the site
        here based solely on query-side homopolymers.
        */
        const uint8_t isHomopolymer = if_is_homopolymer_strict(
            unpacked.data(), (int64_t)queryRead.baseCount, site) ? 1 : 0;

        const uint32_t refReadsFwdSite = (totalReadsFwd >= totalAltReadsFwd) ? (totalReadsFwd - totalAltReadsFwd) : 0;

        /*
        Emit SnpStats rows.

        occ_0 includes +1 for the query read itself, which mirrors hifiasm’s modeling where the
        query contributes one unit of reference support at its own base.
        */
        for (uint8_t altBaseIdx = 0; altBaseIdx < 4; ++altBaseIdx) {
            if (rowIndexForBase[altBaseIdx] == invalid<uint32_t>) {
                continue;
            }

            const uint32_t misCount = misCountPerBase[altBaseIdx];

            SnpStats stat;
            stat.site = site;
            stat.occ_1 = misCount; 
            stat.occ_0 = refCov + 1;
            stat.fwd_ref_cov = refReadsFwdSite + 1;
            stat.refBase = queryRead[site].character();
            stat.altBase = Base::fromInteger(altBaseIdx).character();
            stat.is_homopolymer = isHomopolymer;
            stat.score = -1;

            rowIndexForBase[altBaseIdx] = (uint32_t)snpStats.size();
            snpStats.push_back(stat);
        }

        /*
        Attach each mismatch evidence entry at this site to the emitted SNP row.

        If the mismatch base did not pass the emission filters (e.g. weak third allele),
        overlapSite stays invalid and the evidence will be ignored by later stages.
        */
        for (size_t j = startIdx; j < evidenceIdx; ++j) {
            const uint8_t b = hapEvidence[j].misBase;
            if (b < 4) {
                const uint32_t row = rowIndexForBase[b];
                if (row != invalid<uint32_t>) {
                    hapEvidence[j].overlapSite = row;
                }
            }
        }
    }
}

/*
Bit-parallel approximation of hifiasm’s is_hpc_vec check.

Hifiasm has an additional robustness check that can discount HP-suspect observations and then
re-evaluate whether a site still meets minimal allele-support thresholds.

In Dinara parity-EC we keep the implementation available, but it is typically disabled for
performance: doing per-overlap target-side homopolymer masking can dominate runtime at scale.
When enabled, flatHpBits marks overlaps that are HP-suspect at a site and we subtract those
counts from occ_0/occ_1 to decide if the site would become non-informative.
*/
inline bool isHpcVecMaskedBitParallel(
    const SnpStats& stat,
    size_t siteIdx,
    const uint64_t* flatBits,
    const uint64_t* flatHpBits,
    size_t nWords,
    uint32_t s_hap_cov,
    uint32_t infor_cov
) {
    const uint64_t* rowBits = &flatBits[siteIdx * 2 * nWords];
    const uint64_t* rowHp = &flatHpBits[siteIdx * nWords];
    
    /*
    Count how many homopolymer-marked overlaps contribute to the ref and alt bitsets for this row.

    flatBits encodes ref/alt membership; flatHpBits encodes “HP-suspect” overlaps at this site.
    Taking the intersection gives the number of observations that would be discounted by the
    robustness filter.
    */
    uint64_t n0_hp = 0, n1_hp = 0;
    for (size_t w = 0; w < nWords; ++w) {
        n0_hp += __builtin_popcountll(rowBits[2*w] & rowHp[w]);
        n1_hp += __builtin_popcountll(rowBits[2*w+1] & rowHp[w]);
    }
    
    uint32_t n0_robust = (stat.occ_0 >= n0_hp) ? (uint32_t)(stat.occ_0 - n0_hp) : 0;
    uint32_t n1_robust = (stat.occ_1 >= n1_hp) ? (uint32_t)(stat.occ_1 - n1_hp) : 0;
    
    return (n0_robust < 2 || n1_robust < 2 || n0_robust < s_hap_cov || n1_robust < infor_cov);
}

/*
Strict phasing link test between two SNP rows (comput_sc_rphase).

We decide whether siteI and siteJ can be chained together by looking for evidence of both
haplotypes across overlaps:
  - at least one overlap supports ref at both sites (ref-ref)
  - at least one overlap supports alt at both sites (alt-alt)

We must also reject contradictions:
  - if an overlap has a usable allele state at one site but “unknown/other” at the other site,
    hifiasm treats the link as invalid (it indicates a third-allele pattern).
  - there is a special case in hifiasm: if the same overlap shows “other” at both sites, it is
    treated as ref-ref support (“rareRef”), not a contradiction.

The return value is INT64_MIN for “cannot link”, or a positive value (1) for “can link”.
*/
inline int64_t comput_sc_rphase_strict(
    const uint64_t* flatBits,
    const uint64_t* flatAnyBits,
    size_t siteI, size_t siteJ,
    size_t nWords
) {
    const uint64_t* rowI = &flatBits[siteI * 2 * nWords];
    const uint64_t* rowJ = &flatBits[siteJ * 2 * nWords];
    const uint64_t* anyI = &flatAnyBits[siteI * nWords];
    const uint64_t* anyJ = &flatAnyBits[siteJ * nWords];

    uint64_t p0 = 0, p1 = 0;
    for (size_t k = 0; k < nWords; ++k) {
        const uint64_t rI = rowI[2*k], aI = rowI[2*k+1], oI = anyI[k];
        const uint64_t rJ = rowJ[2*k], aJ = rowJ[2*k+1], oJ = anyJ[k];

        const uint64_t rareRef = oI & oJ;

        const uint64_t effRefI = rI | rareRef;
        const uint64_t effRefJ = rJ | rareRef;

        /*
        Reject “one-sided other” contradictions.

        If an overlap has an allele state at both sites (covered), but one site is “other” while the
        other site is ref/alt, hifiasm discards the link entirely. This enforces a strict biallelic
        model for link scoring. The special case “other at both sites” is handled via rareRef.
        */
        const uint64_t soloOther = (oI ^ oJ) & (rI | aI | oI) & (rJ | aJ | oJ);
        if (soloOther) return INT64_MIN;

        /*
        Reject ref-alt and alt-ref contradictions.

        A single overlap cannot support ref at one site and alt at the other site if we want to
        interpret the pair as part of a consistent phased chain.
        */
        if (effRefI & aJ) return INT64_MIN;
        if (aI & effRefJ) return INT64_MIN;

        const uint64_t rr = effRefI & effRefJ;
        const uint64_t aa = aI & aJ;
        p0 |= rr;
        p1 |= aa;
    }
    
    if (!p0 || !p1) return INT64_MIN;
    return 1;
}

/*
gen_rphase_dp: dynamic programming chaining of SNP rows for one query read.

This stage is the core of parity-EC. It answers: “which SNP rows are mutually consistent enough
to be considered phased/validated for this query read?”.

The implementation mirrors the shape of hifiasm’s gen_rphase_dp + gen_rphase_dp0_single_path:
  - Build a compact matrix representation of allele states per overlap per SNP row.
  - Score links between SNP rows with a strict biallelic consistency test.
  - Run a DP to find chains of mutually consistent rows.
  - Apply hifiasm-style cc filtering based on coveragePeak (hom_cov), n_hap=2, cut_rate=0.7,
    cut_bd=6.
  - Compact snpStats/hapEvidence down to the DP-retained rows and reset row score semantics to
    match hifiasm’s downstream expectations.

Important representation details:
  - nCands overlaps are represented as bitsets with nWords=(nCands+63)/64 64-bit words.
  - For each SNP row we store two bitsets interleaved in flatBits:
        rowRefBits[w] at flatBits[row*2*nWords + 2*w]
        rowAltBits[w] at flatBits[row*2*nWords + 2*w + 1]
    and one bitset flatAnyBits[row*nWords + w] for “other allele” observations at that site.
*/
static void gen_rphase_dp(
    Assembler& assembler,
    HifiasmECScratchPad& scratch,
    RphaseDpTiming* timing = nullptr
) {
    auto& snpStats = scratch.snpStats;
    auto& hapEvidence = scratch.hapEvidence;
    auto& candidates = scratch.candidates;

    if (snpStats.empty() || hapEvidence.empty()) return;

    const size_t nCands = candidates.size();
    const size_t nSites = snpStats.size();
    const size_t nWords = (nCands + 63) / 64;

    /*
    Step 1: build Alt and AnyOther bitsets from the mismatch evidence stream.

    hapEvidence only contains mismatch entries (type==1). For each SNP row we want:
      - which overlaps support the chosen alt base (Alt bitset)
      - which overlaps have a mismatch but not that chosen alt base (AnyOther bitset)

    The “AnyOther” bitset is what allows strict third-allele rejection in comput_sc_rphase_strict.
    */
    auto& flatBits = scratch.flatBits;
    auto& flatAnyBits = scratch.flatAnyBits;
    const auto tBuildAltAnyBegin = timing ? steady_clock::now() : steady_clock::time_point{};
    flatBits.assign(nSites * 2 * nWords, 0);
    flatAnyBits.assign(nSites * nWords, 0);
    
    auto& allMisBits = scratch.supportBits;
    size_t evIdx = 0;
    for (size_t i = 0; i < nSites; ) {
        const uint32_t s = snpStats[i].site;

        /*
        snpStats can contain multiple rows with the same site (multi-allelic). We process the whole
        [i, j) block for this site together so we only scan the hapEvidence block once.
        */
        size_t j = i + 1;
        while (j < nSites && snpStats[j].site == s) ++j;

        /*
        Build a small lookup table from alt base to row index for this site.
        This lets us place mismatch evidence into the correct SNP-row Alt bitset.
        */
        int rowForBase[4] = {-1, -1, -1, -1};
        for (size_t r = i; r < j; ++r) {
            const int b = base2int(snpStats[r].altBase);
            if (b >= 0 && b < 4) rowForBase[b] = int(r);
        }

        allMisBits.assign(nWords, 0);

        /*
        Consume the contiguous hapEvidence block for this site and:
          - set allMisBits for any overlap with any mismatch at this site
          - set per-row Alt bits for overlaps matching each row’s alt base
        */
        while (evIdx < hapEvidence.size() && hapEvidence[evIdx].site < s) ++evIdx;
        size_t nextEv = evIdx;
        while (nextEv < hapEvidence.size() && hapEvidence[nextEv].site == s) {
            const uint32_t ovId = hapEvidence[nextEv].overlapID;
            if (ovId < nCands) {
                const size_t w = ovId >> 6;
                const uint64_t mask = 1ULL << (ovId & 63);

                allMisBits[w] |= mask;

                const uint8_t b = hapEvidence[nextEv].misBase;
                if (b < 4) {
                    const int rowIndex = rowForBase[b];
                    if (rowIndex >= 0) {
                        uint64_t* row = &flatBits[size_t(rowIndex) * 2 * nWords];
                        row[2 * w + 1] |= mask;
                    }
                }
            }
            ++nextEv;
        }

        /*
        AnyOther bits for a row are mismatches at this site that are not the row’s chosen alt base.
        This is computed as (all mismatches) AND NOT (row alt mismatches).
        */
        for (size_t r = i; r < j; ++r) {
            const uint64_t* row = &flatBits[r * 2 * nWords];
            uint64_t* rowAny = &flatAnyBits[r * nWords];
            for (size_t w = 0; w < nWords; ++w) {
                rowAny[w] = allMisBits[w] & ~row[2 * w + 1];
            }
        }

        i = j;
    }

    if (timing) {
        timing->buildAltAnyBits += seconds(steady_clock::now() - tBuildAltAnyBegin);
    }

    /*
    Step 2: build Ref bitsets from coverage intervals.

    A candidate overlap is considered ref-supporting at a SNP row if:
      - it covers the row’s query position (within [qs,qe))
      - it is not inside a query-coordinate “hole” interval at that position
      - it does not have a mismatch at that position (neither the row’s alt base nor AnyOther)

    Coverage is computed with a sweep-line over site-indexed events, maintaining a bitset of
    currently active overlaps. A second sweep-line maintains overlaps that are “gapped” (holes)
    at the current site.
    */
    const auto tBuildRefBegin = timing ? steady_clock::now() : steady_clock::time_point{};
    auto& uniqueSites = scratch.uniqueSites;
    uniqueSites.clear();
    for (const auto& s : snpStats) uniqueSites.push_back(s.site);

    const bool haveInsertionHoles =
        (scratch.insertionOffsets.size() == nCands + 1) && !scratch.insertionIntervals.empty();

    auto& events = scratch.events;
    events.clear();
    for (size_t candIdx = 0; candIdx < nCands; ++candIdx) {
        const auto& cand = candidates[candIdx];
        if (cand.qe <= cand.qs) continue;

        auto itStart = std::lower_bound(uniqueSites.begin(), uniqueSites.end(), cand.qs);
        auto itEnd = std::lower_bound(uniqueSites.begin(), uniqueSites.end(), cand.qe);
        
        size_t idxS = (size_t)std::distance(uniqueSites.begin(), itStart);
        size_t idxE = (size_t)std::distance(uniqueSites.begin(), itEnd);

        if (idxS < idxE) {
            events.push_back({(uint32_t)idxS, (uint32_t)candIdx, false});
            events.push_back({(uint32_t)idxE, (uint32_t)candIdx, true});
        }
    }
    std::sort(events.begin(), events.end());

    /*
    Build events for “hole” intervals so we can maintain a gapped bitset in parallel with active.
    */
    auto& gapEvents = scratch.gapEvents;
    gapEvents.clear();
    if (haveInsertionHoles) {
        for (size_t candIdx = 0; candIdx < nCands; ++candIdx) {
            const uint32_t offBegin = scratch.insertionOffsets[candIdx];
            const uint32_t offEnd = scratch.insertionOffsets[candIdx + 1];
            for (uint32_t ii = offBegin; ii < offEnd; ++ii) {
                const auto& interval = scratch.insertionIntervals[ii];
                auto itStart = std::lower_bound(uniqueSites.begin(), uniqueSites.end(), interval.begin);
                auto itEnd = std::lower_bound(uniqueSites.begin(), uniqueSites.end(), interval.end);
                const size_t idxS = (size_t)std::distance(uniqueSites.begin(), itStart);
                const size_t idxE = (size_t)std::distance(uniqueSites.begin(), itEnd);
                if (idxS < idxE) {
                    gapEvents.push_back({idxS, (uint32_t)candIdx, false});
                    gapEvents.push_back({idxE, (uint32_t)candIdx, true});
                }
            }
        }
        std::sort(gapEvents.begin(), gapEvents.end());
    }

    auto& active = scratch.active;
    auto& gapped = scratch.gapped;
    active.assign(nWords, 0);
    gapped.assign(nWords, 0);

    size_t eventIdx = 0;
    size_t gapEventIdx = 0;
    for (size_t i = 0; i < nSites; ++i) {
        while (eventIdx < events.size() && events[eventIdx].siteIdx == i) {
            uint32_t c = events[eventIdx].candIdx;
            if (events[eventIdx].isEnd) active[c >> 6] &= ~(1ULL << (c & 63));
            else active[c >> 6] |= (1ULL << (c & 63));
            eventIdx++;
        }

        while (gapEventIdx < gapEvents.size() && gapEvents[gapEventIdx].siteIdx == i) {
            uint32_t c = gapEvents[gapEventIdx].candIdx;
            if (gapEvents[gapEventIdx].isEnd) gapped[c >> 6] &= ~(1ULL << (c & 63));
            else gapped[c >> 6] |= (1ULL << (c & 63));
            gapEventIdx++;
        }
        
        uint64_t* row = &flatBits[i * 2 * nWords];
        const uint64_t* rowAny = &flatAnyBits[i * nWords]; 
        for (size_t w = 0; w < nWords; ++w) {
            /*
            Ref-supporting overlaps are those that are currently active, not currently gapped,
            and not mismatching at this site.
            */
            uint64_t allMis = rowAny[w] | row[2*w+1];
            row[2*w] = active[w] & ~gapped[w] & ~allMis;
        }
    }
    if (timing) {
        timing->buildRefBits += seconds(steady_clock::now() - tBuildRefBegin);
    }

    /*
    Step 3: DP chaining and per-row scoring.

    At this point we have bit-matrices that allow O(nWords) linkage tests between any pair of SNP
    rows. We now:
      - select which rows are eligible for DP (must already be “real alleles”: >=s_hap_cov/infor_cov)
      - run DP transitions using comput_sc_rphase_strict
      - extract disjoint best paths and assign per-row score
      - finally compact to score==1 rows and reset score semantics for trans-closure
    */
    const uint32_t s_hap_cov = 3;
    const uint32_t infor_cov = 3;

    /*
    Compute hifiasm-style cc cutoff.

    Hifiasm defines:
      cc = ((het_cov > 0)? het_cov : (hom_cov / n_hap));
      cc *= cut_rate;
      cc = max(cc, cut_bd);

    In this parity implementation:
      - het_cov is not used (we use hom_cov derived from coveragePeak).
      - n_hap is fixed at 2.
      - cut_rate is 0.7 (represented as 7/10).
      - cut_bd is 6.

    cc is applied when turning DP paths into score==1 rows (including singletons).
    */
    uint32_t cc = 6;
    {
        constexpr uint64_t cut_bd = 6;
        constexpr uint64_t cut_rate_num = 7;
        constexpr uint64_t cut_rate_den = 10;
        constexpr uint64_t n_hap = 2;

        const uint64_t hom_cov = assembler.assemblerInfo.isOpen ?
            assembler.assemblerInfo->kmerDistributionInfo.coveragePeak :
            invalid<uint64_t>;

        uint64_t base = 0;
        if (hom_cov != invalid<uint64_t> && hom_cov > 0) {
            base = hom_cov / n_hap;
        }

        uint64_t cc64 = (base * cut_rate_num) / cut_rate_den;
        if (cc64 < cut_bd) cc64 = cut_bd;
        if (cc64 > uint64_t(std::numeric_limits<uint32_t>::max())) cc64 = uint64_t(std::numeric_limits<uint32_t>::max());
        cc = uint32_t(cc64);
    }

    auto& validIndices = scratch.validIndices;
    validIndices.clear();
    /*
    Decide which SNP rows participate in DP.

    SNP rows are created in detectHetSites when occ_1>=2 (push_info-like emission), but hifiasm
    only feeds “real allele” rows into DP: both alleles must meet s_hap_cov/infor_cov (default 3/3).
    */
    for (size_t i = 0; i < nSites; ++i) {
        if (snpStats[i].occ_0 >= s_hap_cov && snpStats[i].occ_1 >= infor_cov) {
            validIndices.push_back((int)i);
        }
        snpStats[i].score = -1;
    }

    if (validIndices.empty()) return;

    const auto tBuildHpBegin = timing ? steady_clock::now() : steady_clock::time_point{};
    /*
    Optional target-side homopolymer masking.

    Hifiasm has a robustness mechanism (is_hpc_vec) that can remove HP-suspect observations on the
    target side and then re-check whether a site still meets support thresholds.

    In Dinara parity-EC this is typically disabled because computing target-side HP masks requires
    heavy random access into many target reads per site and can dominate runtime. The bitset is
    still present so the logic can be enabled for experiments.
    */
    auto& flatHpBits = scratch.flatHpBits;
    flatHpBits.assign(nSites * nWords, 0);
    const bool useOverlapHpMask = false;
    if (timing) {
        timing->buildHpBits += seconds(steady_clock::now() - tBuildHpBegin);
    }

    /*
    DP arrays:
      - f[v] is the best chain score ending at validIndices[v]
      - p[v] stores the predecessor index in that best chain
    */
    const int nV = (int)validIndices.size();
    auto& f = scratch.f; f.assign(nV, 0);
    auto& p = scratch.p; p.assign(nV, -1);
    auto& rowIsValid = scratch.covered;
    rowIsValid.assign(nV, 0);

    /*
    validIndices can contain multiple rows with the same query position (multi-allelic). Linking
    rows at the same position is meaningless and wastes work. We therefore skip transitions that
    would compare within the current “same-site” block in validIndices.
    */
    auto& blockStart = scratch.indexMap;
    blockStart.resize(nV);
    int currBlockStart = 0;
    for (int i = 0; i < nV; ++i) {
        if (i == 0) {
            currBlockStart = 0;
        } else {
            const uint32_t prevSite = snpStats[validIndices[i - 1]].site;
            const uint32_t thisSite = snpStats[validIndices[i]].site;
            if (thisSite != prevSite) currBlockStart = i;
        }
        blockStart[i] = currBlockStart;
    }

    const auto tTransitionsBegin = timing ? steady_clock::now() : steady_clock::time_point{};
    for (int i = 0; i < nV; ++i) {
        int siteI = validIndices[i];
        
        /*
        Optional robustness gate (target-side HP discounting) before allowing this node to
        participate in DP. Disabled by default for performance.
        */
        if (useOverlapHpMask) {
            if (isHpcVecMaskedBitParallel(snpStats[siteI], siteI, flatBits.data(), flatHpBits.data(), nWords, s_hap_cov, infor_cov)) continue;
        }

        f[i] = 1;
        const int jEnd = blockStart[i];
        for (int j = 0; j < jEnd; ++j) {
            if (f[j] == 0) continue;
            int siteJ = validIndices[j];
#ifdef DINARA_TESTING
            if (snpStats[siteI].site == snpStats[siteJ].site) {
                gDpSameSiteComparisons.fetch_add(1, std::memory_order_relaxed);
            }
#endif
            /*
            Link test between the two SNP rows.

            comput_sc_rphase_strict implements the same “must have ref-ref and alt-alt support,
            and must have no contradictions” rule as hifiasm’s comput_sc_rphase.
            */
            int64_t sc = comput_sc_rphase_strict(flatBits.data(), flatAnyBits.data(), siteI, siteJ, nWords);
            if (sc == INT64_MIN) continue;
            
            if (f[j] + (int)sc > f[i]) {
                f[i] = f[j] + (int)sc;
                p[i] = j;
            }
        }
    }
    if (timing) {
        timing->transitions += seconds(steady_clock::now() - tTransitionsBegin);
    }

    /*
    Step 4: extract disjoint best paths from the DP predecessor graph.

    Hifiasm’s single-path DP assigns scores and then extracts paths in descending score order,
    marking visited nodes to avoid reusing the same SNP row in multiple chains.
    */
    const auto tExtractBegin = timing ? steady_clock::now() : steady_clock::time_point{};
    vector<pair<int64_t, int>> scoreIdx;
    for(int i=0; i<nV; ++i) scoreIdx.push_back({f[i], i});
    std::sort(scoreIdx.rbegin(), scoreIdx.rend());

    std::vector<int> pathNodes;
    for(auto& pair : scoreIdx) {
        int curr = pair.second;
        if(f[curr] < 1 || rowIsValid[curr]) continue;

        pathNodes.clear();
        while(curr != -1 && !rowIsValid[curr]) {
            rowIsValid[curr] = 1;
            pathNodes.push_back(curr);
            curr = p[curr];
        }

        /*
        Determine the “plus” score for this extracted path.

        In hifiasm’s no-QV branch:
          - any multi-site path is accepted (plus=1)
          - a singleton path is only accepted if it is not HP-suspect and has occ_0 >= cc
        We mirror that behavior here.
        */
        int plus = -1;
        if (pathNodes.size() > 1) {
            plus = 1;
        } else if (pathNodes.size() == 1) {
            const int siteV = validIndices[pathNodes[0]];
            if (snpStats[siteV].occ_0 >= cc) {
                bool hpMasked = false;
                {
                    const uint32_t querySite = snpStats[siteV].site;
                    const uint64_t* rowBits = &flatBits[siteV * 2 * nWords];
                    auto isGappedAtSite = [&](size_t candIdx, uint32_t sitePos) -> bool {
                        if (!haveInsertionHoles) return false;
                        const uint32_t offBegin = scratch.insertionOffsets[candIdx];
                        const uint32_t offEnd = scratch.insertionOffsets[candIdx + 1];
                        if (offBegin >= offEnd) return false;
                        const auto* beginIt = scratch.insertionIntervals.data() + offBegin;
                        const auto* endIt = scratch.insertionIntervals.data() + offEnd;
                        auto it = std::upper_bound(
                            beginIt, endIt, sitePos,
                            [](uint32_t value, const HifiasmECScratchPad::GapInterval& iv) { return value < iv.begin; }
                        );
                        if (it == beginIt) return false;
                        --it;
                        return (it->begin <= sitePos && sitePos < it->end);
                    };

                    uint64_t hpRef = 0, hpAlt = 0;
                    for (size_t candIdx = 0; candIdx < nCands; ++candIdx) {
                        const auto& cand = candidates[candIdx];
                        if (querySite < cand.qs || querySite >= cand.qe) continue;
                        if (isGappedAtSite(candIdx, querySite)) continue;

                        const auto& view = assembler.getReads().getRead(cand.targetId);
                        const uint32_t targetSite = cand.isRev ?
                            (uint32_t)(cand.te - (querySite - cand.qs) - 1) :
                            (uint32_t)(cand.ts + (querySite - cand.qs));
                        if (targetSite >= view.baseCount) continue;

                        if (!hpc_mask_ff(view, (int64_t)view.baseCount, targetSite)) continue;

                        const size_t w = candIdx >> 6;
                        const uint64_t mask = 1ULL << (candIdx & 63);
                        hpRef += __builtin_popcountll(rowBits[2 * w] & mask);
                        hpAlt += __builtin_popcountll(rowBits[2 * w + 1] & mask);
                    }

                    const uint32_t n0_robust = (snpStats[siteV].occ_0 >= hpRef) ? (uint32_t)(snpStats[siteV].occ_0 - hpRef) : 0;
                    const uint32_t n1_robust = (snpStats[siteV].occ_1 >= hpAlt) ? (uint32_t)(snpStats[siteV].occ_1 - hpAlt) : 0;
                    hpMasked = (n0_robust < 2 || n1_robust < 2 || n0_robust < s_hap_cov || n1_robust < infor_cov);
                }
                if (!hpMasked) plus = 1;
            }
        }

        for (int v : pathNodes) {
            const int snpIdx = validIndices[v];
            snpStats[snpIdx].score = (snpStats[snpIdx].occ_0 >= cc) ? plus : -1;
        }
    }
    if (timing) {
        timing->extractPaths += seconds(steady_clock::now() - tExtractBegin);
    }

    /*
    Step 5: compact the SNP row list to DP-retained rows and remap evidence.

    Hifiasm compacts snp_stat after DP so only retained rows continue into trans-closure. We do
    the same by keeping only rows with score==1, then:
      - build an oldRowIndex -> newRowIndex map
      - keep hapEvidence entries only for query positions that still exist in the compacted list
      - remap overlapSite indices through the map (or invalidate them if the row was dropped)

    Parity detail: after compaction, hifiasm resets snp_stat.score back to -1 before
    generate_haplotypes_naive_HiFi. We preserve the DP decision in dpScore and then reset score.
    */
    {
        const size_t oldNSites = snpStats.size();
        auto& indexMap = scratch.indexMap;
        indexMap.assign(oldNSites, -1);

        auto& keptSites = scratch.uniqueSites;
        keptSites.clear();
        keptSites.reserve(oldNSites);

        size_t write = 0;
        for (size_t i = 0; i < oldNSites; ++i) {
            if (snpStats[i].score != 1) continue;
            indexMap[i] = int(write);
            keptSites.push_back(snpStats[i].site);
            SnpStats kept = snpStats[i];
            kept.dpScore = kept.score;
            kept.score = -1;
            snpStats[write] = kept;
            ++write;
        }
        snpStats.resize(write);

        if (snpStats.empty()) {
            hapEvidence.clear();
            return;
        }

        std::sort(keptSites.begin(), keptSites.end());
        keptSites.erase(std::unique(keptSites.begin(), keptSites.end()), keptSites.end());

        size_t g = 0;
        size_t evWrite = 0;
        for (size_t i = 0; i < hapEvidence.size(); ++i) {
            const uint32_t site = hapEvidence[i].site;
            while (g < keptSites.size() && keptSites[g] < site) ++g;
            if (g == keptSites.size() || keptSites[g] != site) continue;

            auto ev = hapEvidence[i];
            if (ev.overlapSite < oldNSites) {
                const int mapped = indexMap[ev.overlapSite];
                ev.overlapSite = (mapped >= 0) ? uint32_t(mapped) : invalid<uint32_t>;
            } else {
                ev.overlapSite = invalid<uint32_t>;
            }
            hapEvidence[evWrite++] = ev;
        }
        hapEvidence.resize(evWrite);
    }
}

/*
generate_haplotypes_naive_HiFi: hifiasm-style trans-closure for SNP sites.

This stage consumes the post-DP SNP row list and marks overlaps as trans (is_match=2) when they
carry alt alleles at validated sites. It also applies hifiasm’s “not real allele” adjustment:
when an overlap is marked trans, it should no longer contribute ref support to other sites, so
we decrement occ_0 for sites covered by that overlap where it does not carry the alt allele.

The high-level flow matches hifiasm’s Correct.cpp:
  0) Drop adjacent SNP sites (distance 1) and rewrite evidence row indices.
  1) Sort evidence by overlap and count how many informative alleles each overlap carries.
  2) Seed trans overlaps from the most-informative overlaps and promote supported sites.
  3) Second pass: any remaining overlap that hits a promoted site becomes trans.
  4) Optional multi_check: rescue dense weak patterns and re-run final marking.

This function updates candidates[].is_match, and updates snpStats[].occ_0 as part of the not-real
allele decrements. The final application of deletions happens later in performHifiasmECParity.
*/
static void generate_haplotypes_naive_HiFi(
    Assembler& /* assembler */,
    HifiasmECScratchPad& scratch
) {
    auto& candidates = scratch.candidates;
    auto& snpStats = scratch.snpStats;
    auto& hapEvidence = scratch.hapEvidence;
    const auto& insertionOffsets = scratch.insertionOffsets;
    const auto& insertionIntervals = scratch.insertionIntervals;
    const auto& insertionBaseCount = scratch.insertionBaseCount;

    const size_t nCands = candidates.size();
    size_t nSites = snpStats.size();

    /*
    Reset overlap states for this stage.

    Earlier stages may have used is_match for intermediate bookkeeping. In hifiasm, SNP-stage
    trans-closure starts from a clean “all overlaps are cis” state.
    */
    for (auto& cand : candidates) cand.is_match = 1;

    if (nSites == 0 || nCands == 0 || hapEvidence.empty()) {
        return;
    }

    /*
    Constants.

    These are the correction-mode defaults used by hifiasm for the HiFi error-correction stage.
    */
    constexpr uint32_t s_hap_cov = 3;
    constexpr uint32_t infor_cov = 3;
    constexpr uint32_t multi_check_distance = 32;
    constexpr uint64_t up_num = 4;
    constexpr uint64_t up_den = 100;
    constexpr bool enable_multi_check = true;

    /*
    Step 0: drop adjacent SNP sites.

    Hifiasm removes SNP sites that are immediately adjacent (distance 1). This is a conservative
    filter to avoid dense local noise patterns. Because snpStats is row-based (multi-allelic rows),
    we first build “site groups” so we can remove an entire position (all rows at that position)
    in one operation and remap row indices in hapEvidence.
    */
    auto& siteGroups = scratch.rawSVs;
    siteGroups.clear();
    for (uint32_t i = 0; i < nSites; ) {
        const uint32_t site = snpStats[i].site;
        uint32_t j = i + 1;
        while (j < nSites && snpStats[j].site == site) ++j;
        RawSV g;
        g.overlapID = i;
        g.site = site;
        g.size = j;
        siteGroups.push_back(g);
        i = j;
    }

    /*
    Mark which site groups are kept.

    A group is dropped if its site is adjacent to either the previous or the next site group.
    */
    auto& keepGroup = scratch.covered;
    keepGroup.assign(siteGroups.size(), 1);
    for (size_t g = 0; g < siteGroups.size(); ++g) {
        const uint32_t site = siteGroups[g].site;
        const uint32_t prevSite = (g > 0) ? siteGroups[g - 1].site : invalid<uint32_t>;
        const uint32_t nextSite = (g + 1 < siteGroups.size()) ? siteGroups[g + 1].site : invalid<uint32_t>;
        const bool adjacent =
            (prevSite != invalid<uint32_t> && prevSite + 1 == site) ||
            (nextSite != invalid<uint32_t> && site + 1 == nextSite);
        if (adjacent) keepGroup[g] = 0;
    }

    /*
    Build a mapping from old SNP-row indices to new compacted indices.

    Each kept site group can contain multiple SNP rows (one per alt base). The mapping is applied
    to hapEvidence.overlapSite for evidence entries that remain after site dropping.
    */
    auto& rowMap = scratch.indexMap;
    rowMap.assign((int)nSites, -1);
    size_t newN = 0;
    for (size_t g = 0; g < siteGroups.size(); ++g) {
        if (!keepGroup[g]) continue;
        const uint32_t begin = siteGroups[g].overlapID;
        const uint32_t end = (uint32_t)siteGroups[g].size;
        for (uint32_t r = begin; r < end; ++r) {
            rowMap[(int)r] = (int)newN++;
        }
    }

    if (newN == 0) {
        snpStats.clear();
        hapEvidence.clear();
        return;
    }

    /*
    Compact snpStats in-place to only the rows belonging to kept site groups.
    */
    {
        size_t write = 0;
        for (size_t g = 0; g < siteGroups.size(); ++g) {
            if (!keepGroup[g]) continue;
            const uint32_t begin = siteGroups[g].overlapID;
            const uint32_t end = (uint32_t)siteGroups[g].size;
            for (uint32_t r = begin; r < end; ++r) {
                snpStats[write++] = snpStats[r];
            }
        }
        snpStats.resize(newN);
    }

    /*
    Compact hapEvidence:
      - drop evidence entries for sites that were removed as adjacent
      - remap overlapSite indices through rowMap for entries that remain
    */
    {
        size_t g = 0;
        size_t write = 0;
        for (size_t i = 0; i < hapEvidence.size(); ) {
            const uint32_t site = hapEvidence[i].site;
            while (g < siteGroups.size() && siteGroups[g].site < site) ++g;
            const bool keepSite = (g < siteGroups.size() && siteGroups[g].site == site && keepGroup[g]);

            size_t j = i + 1;
            while (j < hapEvidence.size() && hapEvidence[j].site == site) ++j;

            if (keepSite) {
                for (size_t k = i; k < j; ++k) {
                    auto ev = hapEvidence[k];
                    if (ev.overlapSite != invalid<uint32_t> && ev.overlapSite < nSites) {
                        const int mapped = rowMap[(int)ev.overlapSite];
                        ev.overlapSite = (mapped >= 0) ? (uint32_t)mapped : invalid<uint32_t>;
                    }
                    hapEvidence[write++] = ev;
                }
            }
            i = j;
        }
        hapEvidence.resize(write);
    }

    nSites = snpStats.size();
    if (nSites == 0 || hapEvidence.empty()) return;

    /*
    Rebuild siteGroups on the compacted snpStats.

    After adjacency dropping, later parts of trans-closure need per-site group boundaries for:
      - “not real allele” occ_0 decrements
      - multi_check’s neighbor-distance filters
    */
    siteGroups.clear();
    for (uint32_t i = 0; i < nSites; ) {
        const uint32_t site = snpStats[i].site;
        uint32_t j = i + 1;
        while (j < nSites && snpStats[j].site == site) ++j;
        RawSV g;
        g.overlapID = i;
        g.site = site;
        g.size = j;
        siteGroups.push_back(g);
        i = j;
    }

    /*
    Step 1: group evidence entries by overlap.

    Many downstream operations are “per overlap”: count how many informative sites a particular
    overlap hits, and scan all mismatch entries for that overlap. hapEvidence is sorted by site,
    so we build an overlapID -> [begin,end) index range by counting and then filling a permutation
    array (perm) into the existing hapEvidence vector.
    */
    auto& ovOffsets = scratch.siteTotalCov;
    auto& ovCursor = scratch.siteFwdCov;
    auto& perm = scratch.svIndices;
    ovOffsets.assign(nCands + 1, 0);

    for (const auto& ev : hapEvidence) {
        if (ev.overlapID >= nCands) continue;
        ++ovOffsets[ev.overlapID + 1];
    }
    for (size_t i = 1; i < ovOffsets.size(); ++i) {
        ovOffsets[i] += ovOffsets[i - 1];
    }
    ovCursor = ovOffsets;
    perm.resize(hapEvidence.size());
    for (size_t evIndex = 0; evIndex < hapEvidence.size(); ++evIndex) {
        const auto& ev = hapEvidence[evIndex];
        if (ev.overlapID >= nCands) continue;
        const size_t out = ovCursor[ev.overlapID]++;
        perm[out] = evIndex;
    }

    /*
    Step 2: compute an overlap ordering by “informativeness”.

    Hifiasm seeds trans overlaps by examining overlaps with many informative SNP alleles first.
    We compute a count per overlap (o) and sort overlaps by descending o.
    */
    auto& overlapSort = scratch.supportBits;
    overlapSort.clear();
    overlapSort.reserve(nCands);

    auto countInformativeMismatches = [&](size_t c) -> uint32_t {
        const size_t begin = ovOffsets[c];
        const size_t end = ovOffsets[c + 1];
        uint32_t o = 0;
        for (size_t k = begin; k < end; ++k) {
            const auto& ev = hapEvidence[perm[k]];
            const uint32_t row = ev.overlapSite;
            if (row == invalid<uint32_t> || row >= nSites) continue;
            const auto& s = snpStats[row];
            if (s.occ_0 < 2 || s.occ_1 < 2) continue;
            if (s.occ_0 >= s_hap_cov && s.occ_1 >= infor_cov) ++o;
        }
        return o;
    };

    for (size_t c = 0; c < nCands; ++c) {
        const uint32_t o = countInformativeMismatches(c);
        if (o == 0) continue;
        const uint64_t key = (uint64_t(uint32_t(-1) - o) << 32) | uint64_t(uint32_t(c));
        overlapSort.push_back(key);
    }
    std::sort(overlapSort.begin(), overlapSort.end());

    /*
    Helper: apply “not real allele” decrements for a trans overlap.

    In hifiasm’s evidence model, each overlap contributes explicit “ref” entries (hh_tp==0) at sites
    it covers. When an overlap is marked trans, those ref entries should no longer be counted as
    support for the query base at other sites. Hifiasm decrements occ_0 at those sites accordingly.

    Dinara does not store explicit ref entries. Instead we infer “ref at site” for this overlap as:
      - site is covered by the overlap interval [qs,qe)
      - site is not inside a query-coordinate hole interval for this overlap
      - the overlap does not have a mismatch at this site (no hapEvidence entry at that position)
    For any such covered-and-not-mismatching site, we decrement occ_0 for all SNP rows at that site.
    */
    auto decrementOcc0ForTransOverlap = [&](size_t c) {
        const auto& cand = candidates[c];
        const uint32_t qs = (uint32_t)cand.qs;
        const uint32_t qe = (uint32_t)cand.qe;
        if (qe <= qs) return;

        const bool haveInsertionHoles =
            (insertionOffsets.size() == nCands + 1) && !insertionIntervals.empty();
        uint32_t holeIdx = 0;
        uint32_t holeEnd = 0;
        if (haveInsertionHoles) {
            holeIdx = insertionOffsets[c];
            holeEnd = insertionOffsets[c + 1];
        }

        /*
        Collect all mismatch sites for this overlap.

        We collect by query position, not by row index, because “not real allele” decrements are
        applied per site group.
        */
        auto& mismatchSites = scratch.uniqueSites;
        mismatchSites.clear();
        {
            const size_t begin = ovOffsets[c];
            const size_t end = ovOffsets[c + 1];
            mismatchSites.reserve(end - begin);
            for (size_t k = begin; k < end; ++k) {
                mismatchSites.push_back(hapEvidence[perm[k]].site);
            }
            std::sort(mismatchSites.begin(), mismatchSites.end());
            mismatchSites.erase(std::unique(mismatchSites.begin(), mismatchSites.end()), mismatchSites.end());
        }

        /*
        Iterate site groups that fall within this overlap’s query span.
        */
        size_t g0 = 0;
        while (g0 < siteGroups.size() && siteGroups[g0].site < qs) ++g0;
        size_t m = 0;
        for (size_t g = g0; g < siteGroups.size(); ++g) {
            const uint32_t site = siteGroups[g].site;
            if (site >= qe) break;

            if (haveInsertionHoles) {
                while (holeIdx < holeEnd && insertionIntervals[holeIdx].end <= site) ++holeIdx;
                if (holeIdx < holeEnd) {
                    const auto& interval = insertionIntervals[holeIdx];
                    if (interval.begin <= site && site < interval.end) continue;
                }
            }

            while (m < mismatchSites.size() && mismatchSites[m] < site) ++m;
            if (m < mismatchSites.size() && mismatchSites[m] == site) continue;

            const uint32_t beginRow = siteGroups[g].overlapID;
            const uint32_t endRow = (uint32_t)siteGroups[g].size;
            for (uint32_t r = beginRow; r < endRow; ++r) {
                uint32_t& occ0 = snpStats[r].occ_0;
                occ0 = (occ0 > 1) ? (occ0 - 1) : 1U;
            }
        }
    };

    /*
    Step 3: seed trans overlaps in sorted order.

    For each overlap in overlapSort:
      - mark it trans (is_match=2)
      - promote any SNP rows it hits by setting snpStats[row].score = 1
      - apply not-real-allele decrements immediately so later overlaps see updated occ_0
    */
    for (uint64_t key : overlapSort) {
        const size_t c = size_t(uint32_t(key));
        const uint32_t o = countInformativeMismatches(c);
        if (o == 0) continue;

        candidates[c].is_match = 2;

        const size_t begin = ovOffsets[c];
        const size_t end = ovOffsets[c + 1];
        for (size_t k = begin; k < end; ++k) {
            const uint32_t row = hapEvidence[perm[k]].overlapSite;
            if (row == invalid<uint32_t> || row >= nSites) continue;
            snpStats[row].score = 1;
        }

        decrementOcc0ForTransOverlap(c);
    }

    /*
    Step 4: closure pass.

    Any overlap that is still cis (is_match==1) but has a mismatch at a promoted site becomes trans.
    */
    for (size_t c = 0; c < nCands; ++c) {
        if (candidates[c].is_match == 2) continue;
        const size_t begin = ovOffsets[c];
        const size_t end = ovOffsets[c + 1];
        for (size_t k = begin; k < end; ++k) {
            const uint32_t row = hapEvidence[perm[k]].overlapSite;
            if (row == invalid<uint32_t> || row >= nSites) continue;
            if (snpStats[row].score == 1) {
                candidates[c].is_match = 2;
                break;
            }
        }
    }

    /*
    Step 5: reset per-site temporary score marks for overlaps that stayed cis.

    Hifiasm resets snp_stat.score back to -1 for sites visited by overlaps that were not marked
    trans. This prevents “score==1” from leaking across overlaps in later logic.
    */
    for (size_t c = 0; c < nCands; ++c) {
        if (candidates[c].is_match == 2) continue;
        const size_t begin = ovOffsets[c];
        const size_t end = ovOffsets[c + 1];
        for (size_t k = begin; k < end; ++k) {
            const uint32_t row = hapEvidence[perm[k]].overlapSite;
            if (row == invalid<uint32_t> || row >= nSites) continue;
            snpStats[row].score = -1;
        }
    }

    /*
    Step 6: multi_check recovery (optional).

    Hifiasm’s multi_check is a recovery path for cases where an overlap carries many weak SNP
    alleles that do not meet the “real allele” threshold (occ_0/occ_1 >= s_hap_cov/infor_cov),
    but the density pattern is informative enough to treat some of them as valid.

    The logic is:
      - consider only overlaps that are still cis after the first trans-closure pass
      - collect candidate “weak rows” in that overlap (occ_0/occ_1>=2 but below real-allele threshold)
      - require enough weak rows relative to overlap alignment length (o/alignLen >= up)
      - apply a neighbor-distance filter so only isolated weak sites are rescued (no +/-32bp neighbors)
      - if a weak row is rescued by at least two overlaps, set its score=1
    */
    if (enable_multi_check) {
        auto& rescued = scratch.path;
        rescued.clear();

        auto& localRows = scratch.uniqueSites;
        auto& keptRows = scratch.diffTotal;
        localRows.clear();
        keptRows.clear();

        for (size_t c = 0; c < nCands; ++c) {
            if (candidates[c].is_match == 2) continue;
            const size_t begin = ovOffsets[c];
            const size_t end = ovOffsets[c + 1];
            if (begin >= end) continue;

            const uint64_t rawLen = (candidates[c].qe > candidates[c].qs) ? (candidates[c].qe - candidates[c].qs) : 0;
            const uint64_t holeLen = (c < insertionBaseCount.size()) ? insertionBaseCount[c] : 0;
            const uint64_t alignLen = (rawLen > holeLen) ? (rawLen - holeLen) : 0;
            if (alignLen == 0) continue;

            localRows.clear();
            for (size_t k = begin; k < end; ++k) {
                const auto& ev = hapEvidence[perm[k]];
                const uint32_t row = ev.overlapSite;
                if (row == invalid<uint32_t> || row >= nSites) continue;
                const auto& s = snpStats[row];
                if (s.occ_0 < 2 || s.occ_1 < 2) continue;
                if (s.occ_0 >= s_hap_cov && s.occ_1 >= infor_cov) continue;
                if (s.score == 1) continue;
                localRows.push_back(row);
            }

            const uint64_t o = localRows.size();
            if (o < 2) continue;
            if (o * up_den < alignLen * up_num) continue;

            std::sort(localRows.begin(), localRows.end());

            keptRows.clear();
            keptRows.reserve(localRows.size());
            for (size_t i = 0; i < localRows.size(); ++i) {
                const uint32_t row = localRows[i];
                const uint32_t site = snpStats[row].site;
                if (i > 0) {
                    const uint32_t prev = snpStats[localRows[i - 1]].site;
                    if (prev + multi_check_distance > site) continue;
                }
                if (i + 1 < localRows.size()) {
                    const uint32_t next = snpStats[localRows[i + 1]].site;
                    if (site + multi_check_distance > next) continue;
                }
                keptRows.push_back(int32_t(row));
            }

            if (keptRows.size() >= 2) {
                for (int32_t r : keptRows) rescued.push_back(uint32_t(r));
            }
        }

        if (!rescued.empty()) {
            std::sort(rescued.begin(), rescued.end());
            for (size_t i = 0, j = 0; i < rescued.size(); i = j) {
                const uint32_t row = rescued[i];
                j = i + 1;
                while (j < rescued.size() && rescued[j] == row) ++j;
                if (j - i >= 2) snpStats[row].score = 1;
            }
        }
    }

    /*
    Step 7: final trans marking after multi_check.

    After multi_check potentially marks additional rows score=1, we perform one more pass:
    any overlap that is still cis but has a mismatch at a score==1 row becomes trans.
    */
    for (size_t c = 0; c < nCands; ++c) {
        if (candidates[c].is_match == 2) continue;
        const size_t begin = ovOffsets[c];
        const size_t end = ovOffsets[c + 1];
        for (size_t k = begin; k < end; ++k) {
            const uint32_t row = hapEvidence[perm[k]].overlapSite;
            if (row == invalid<uint32_t> || row >= nSites) continue;
            const auto& s = snpStats[row];
            if (s.score != 1) continue;
            if (s.occ_0 < 2 || s.occ_1 < 2) continue;
            candidates[c].is_match = 2;
            break;
        }
    }
}

/*
compactPhasedSites: keep only validated SNP rows and remap evidence.

After the SNP pipeline finishes, many candidate rows remain with score==-1. Downstream stages
(SV marking and final “informative read” classification) only care about validated rows.

This helper performs a pure bookkeeping transform:
  1) snpStats: keep only rows with score==1, densely packed at the front.
  2) hapEvidence: remap overlapSite to the new compacted row index; invalidate evidence that
     pointed to discarded rows.

No validation/scoring decisions are made here; semantics are preserved.
*/
static void compactPhasedSites(HifiasmECScratchPad& scratch) {
    auto& snpStats = scratch.snpStats;
    auto& hapEvidence = scratch.hapEvidence;
    const size_t nSites = snpStats.size();
    auto& indexMap = scratch.indexMap;
    indexMap.assign(nSites, -1);
    
    size_t writeIdx = 0;
    for (size_t i = 0; i < nSites; ++i) {
        if (snpStats[i].score == 1) {
            indexMap[i] = (int)writeIdx;
            snpStats[writeIdx] = snpStats[i];
            writeIdx++;
        }
    }
    snpStats.resize(writeIdx);

    for (auto& ev : hapEvidence) {
        if (ev.overlapSite < nSites) {
            int newIdx = indexMap[ev.overlapSite];
            ev.overlapSite = (newIdx >= 0) ? (uint32_t)newIdx : invalid<uint32_t>;
        } else {
            ev.overlapSite = invalid<uint32_t>;
        }
    }
}

/*
detectSVSites: build per-read SV “sites” (large indels) and per-overlap SV evidence.

Dinara stores SNP mismatches and indels in separate evidence streams. This stage consumes the indel
stream, clusters large events into SV “site” rows, and emits:
  - svStats: per-site rows with allele support counts (occ_0/occ_1) in an SNP-like schema
  - svEvidence: per-overlap entries that indicate an overlap supports the SV allele at a site

Downstream SV marking reuses the same “overlap -> trans” state machine style as the SNP path.
*/
template<typename AlignmentContainer>
static void detectSVSites(
    Assembler& assembler,
    const Reads& reads, 
    uint64_t queryReadId,
    const AlignmentContainer& alignmentData, 
    HifiasmECScratchPad& scratch
) {
    auto& candidates = scratch.candidates;
    auto& svEvidence = scratch.svEvidence;
    auto& svStats = scratch.svStats;
    auto& rawSVs = scratch.rawSVs;
    auto& seen = scratch.covered;
    rawSVs.clear();
    svEvidence.clear();
    svStats.clear();

    /*
    SV clustering parameters.

    These are chosen to match the “large indel recovery” behavior used by hifiasm:
      - SV_MIN_LEN: minimum indel length to be considered an SV allele candidate
      - SV_WINDOW: clustering window on the query coordinate
      - SV_SIZE_RATIO: relative size tolerance within a cluster
      - SV_MERGE_GAP: when projected alignment fragments a long indel, merge nearby pieces
    */
    const int32_t SV_MIN_LEN = 20;
    const int32_t SV_WINDOW = 50;
    const double SV_SIZE_RATIO = 0.20;
    const uint32_t SV_MERGE_GAP = 2;

    /*
    Step 1: collect raw large-indel events from the per-overlap indel evidence stream.

    We only consider overlaps that are still kept after SNP trans-closure (is_match==1), mirroring
    hifiasm’s rphase_lidel which filters on overlap_list->list[ii].is_match == 1.

    For each overlap we:
      - scan indel evidence within its [qs,qe) span
      - merge adjacent same-type indel pieces within the overlap
      - record one RawSV entry per merged run if its length is >= SV_MIN_LEN
    */
    for(size_t k = 0; k < candidates.size(); ++k) {
        if (candidates[k].is_match != 1) continue;

        const auto& cand = candidates[k];
        const auto& ad = alignmentData[cand.alignmentId];
        size_t evidenceId = ad.info.alignmentId;
        
        if(evidenceId == invalid<size_t>) continue;

        span<const IndelEvidence> indels;
        if(ad.readIds[1] == queryReadId) {
            indels = assembler.alignedEvidenceStore.getIndels0(evidenceId);
        } else if(ad.readIds[0] == queryReadId) {
            indels = assembler.alignedEvidenceStore.getIndels1(evidenceId);
        } else continue;

        const uint32_t c_qs = (uint32_t)cand.qs;
        const uint32_t c_qe = (uint32_t)cand.qe;
        if (c_qe <= c_qs || indels.empty()) continue;

        /*
        IndelEvidence slices are stored in non-decreasing query-coordinate order, so we can
        lower_bound into the region and then break once we pass c_qe.
        */
        auto it = std::lower_bound(
            indels.begin(), indels.end(), c_qs,
            [](const IndelEvidence& e, uint32_t value) { return e.pos() < value; }
        );

        /*
        Merge split gap events within an overlap.

        Dinara’s projected alignment can fragment a single biological indel into multiple adjacent
        events (e.g., at segment boundaries). We merge those back together before clustering.
        */
        bool haveRun = false;
        uint8_t runType = 0;
        uint32_t runPos = 0;
        uint32_t runLenSum = 0;
        uint32_t runEnd = 0;
        auto flushRun = [&] {
            if (!haveRun) return;
            if (runLenSum >= (uint32_t)SV_MIN_LEN) {
                RawSV sv;
                sv.overlapID = (uint32_t)k;
                sv.site = runPos;
                sv.size = (runType == 0) ? (int64_t)runLenSum : -(int64_t)runLenSum;
                rawSVs.push_back(sv);
            }
            haveRun = false;
            runLenSum = 0;
            runEnd = 0;
        };

        for (; it != indels.end(); ++it) {
            const uint32_t pos = it->pos();
            if (pos >= c_qe) break;
            const uint32_t len = it->len();
            const uint8_t type = it->type();
            const uint32_t endPos = pos + len;

            if (!haveRun) {
                haveRun = true;
                runType = type;
                runPos = pos;
                runLenSum = len;
                runEnd = endPos;
                continue;
            }

            if (type == runType && pos <= runEnd + SV_MERGE_GAP) {
                runLenSum += len;
                runEnd = endPos;
            } else {
                flushRun();
                haveRun = true;
                runType = type;
                runPos = pos;
                runLenSum = len;
                runEnd = endPos;
            }
        }
        flushRun();
    }
    
    if(rawSVs.empty()) return;
    
    /*
    Step 2: cluster raw events by query position and validate frequency/bias.

    We cluster by query coordinate because all later parity logic is performed in query space.
    */
    std::sort(rawSVs.begin(), rawSVs.end(), [](const RawSV& a, const RawSV& b){
        return a.site < b.site;
    });
    
    for(size_t i = 0; i < rawSVs.size(); ) {
        uint32_t startPos = rawSVs[i].site;
        int64_t refSize = rawSVs[i].size; 
        
        size_t j = i;
        uint32_t supportCount = 0;
        uint32_t svFwdCount = 0;
        auto& supportOverlaps = scratch.path;
        supportOverlaps.clear();
        
        /*
        Gather unique supporting overlaps for this seed within the cluster window.
        */
        while(j < rawSVs.size() && rawSVs[j].site < startPos + SV_WINDOW) {
            int64_t sz = rawSVs[j].size;
            if ((sz > 0) == (refSize > 0)) {
                int64_t diff = std::abs(sz - refSize);
                if (diff <= std::abs(refSize) * SV_SIZE_RATIO) {
                    const uint32_t ov = rawSVs[j].overlapID;
                    if (ov < candidates.size()) {
                        if (seen.size() < candidates.size()) seen.assign(candidates.size(), 0);
                        if (!seen[ov]) {
                            seen[ov] = 1;
                            supportOverlaps.push_back(ov);
                            ++supportCount;
                            if (!candidates[ov].isRev) ++svFwdCount;
                        }
                    }
                }
            }
            j++;
        }
        
        /*
        Minimum support and strand-balance filter.

        This mirrors the intent of hifiasm’s SV recovery path: require at least 3 supporting overlaps
        and reject extreme strand bias at low counts.
        */
        if(supportCount >= 3) {
            if (supportCount > 5) {
                if (svFwdCount < 1 || svFwdCount >= supportCount) {
                     if (svFwdCount < supportCount * 0.05 || svFwdCount > supportCount * 0.95) {
                         i++; continue; 
                     }
                }
            } else if (svFwdCount == 0 || svFwdCount == supportCount) {
                 i++; continue;
            }

            SnpStats stat;
            stat.site = startPos;
            uint32_t totalCov = 0;
            for (const auto& cand2 : candidates) {
                if (cand2.is_match != 1) continue;
                if (startPos >= cand2.qs && startPos < cand2.qe) ++totalCov;
            }
            const uint32_t refCov = (totalCov >= supportCount) ? (totalCov - supportCount) : 0;
            stat.occ_1 = supportCount;
            stat.occ_0 = refCov + 1;
            stat.fwd_ref_cov = 0;
            stat.refBase = 'N';
            stat.altBase = 'N';
            stat.is_homopolymer = 0;
            stat.score = -1;
            
            /*
            Emit one per-overlap evidence entry for this SV site.

            svStats holds only aggregate support counts; svEvidence stores which overlaps support
            the SV allele and ties them back to the site row (overlapSite).
            */
            for (uint32_t ov : supportOverlaps) {
                HaplotypeEvidence ev;
                ev.overlapID = ov;
                ev.site = startPos;
                ev.type = 2;
                ev.overlapSite = (uint32_t)svStats.size();
                svEvidence.push_back(ev);
            }
            svStats.push_back(stat);
            i = j;
        } else {
            i++;
        }

        /*
        Reset per-overlap “seen” marks so the next cluster can deduplicate independently without
        allocating a fresh array.
        */
        if (!supportOverlaps.empty()) {
            for (uint32_t ov : supportOverlaps) {
                if (ov < seen.size()) seen[ov] = 0;
            }
        }
    }
}

/*
generate_haplotypes_sv: mark trans overlaps based on SV evidence (hifiasm parity).

This stage is the SV analogue of generate_haplotypes_naive_HiFi. It consumes:
  - svStats: SV site rows with support counts and score state
  - svEvidence: per-overlap SV allele observations (overlapID, site, overlapSite)
  - candidates: overlap state machine fields (is_match) that survive the SNP stage

The algorithm mirrors hifiasm’s overlap-driven marking:
  1) Sort overlaps by number of informative SV alleles (descending).
  2) Seed trans overlaps (is_match=2) from the sorted list.
  3) Promote SV rows touched by trans overlaps to score==1.
  4) Apply the “not real allele” decrement: for a trans overlap, decrement occ_0 at SV sites it
     covers but does not carry an SV event (discounts ref-support for that overlap).
  5) Propagate trans status: any remaining cis overlap carrying an SV at any promoted row becomes
     trans as well.
  6) Reset: rows that remain only supported by overlaps still labeled cis are set back to score==-1.
*/
static void generate_haplotypes_sv(
    Assembler& /* assembler */,
    HifiasmECScratchPad& scratch
) {
    auto& svEvidence = scratch.svEvidence;
    auto& candidates = scratch.candidates;
    auto& svStats = scratch.svStats;

    if (svEvidence.empty() || svStats.empty() || candidates.empty()) return;

    /*
    Minimal-support thresholds for “informative allele” counting when sorting overlaps.
    */
    constexpr uint32_t s_hap_cov = 3;
    constexpr uint32_t infor_cov = 3;

    /*
    Initialize site scores to “not validated” before the seeding pass.
    */
    for (auto& s : svStats) s.score = -1;

    const size_t nCands = candidates.size();
    const size_t nSites = svStats.size();

    /*
    Group svEvidence by overlapID using a counting partition (no global sort).

    After this:
      - ovOffsets[c]..ovOffsets[c+1] indexes perm entries for overlap c
      - perm[k] is an index into svEvidence

    We reuse scratch buffers (siteTotalCov/siteFwdCov/svIndices) to avoid per-read allocations.
    */
    auto& ovOffsets = scratch.siteTotalCov;
    auto& ovCursor = scratch.siteFwdCov;
    auto& perm = scratch.svIndices;
    ovOffsets.assign(nCands + 1, 0);
    for (const auto& ev : svEvidence) {
        if (ev.overlapID < nCands) ++ovOffsets[ev.overlapID + 1];
    }
    for (size_t i = 1; i < ovOffsets.size(); ++i) ovOffsets[i] += ovOffsets[i - 1];
    ovCursor = ovOffsets;
    perm.resize(svEvidence.size());
    for (size_t evIndex = 0; evIndex < svEvidence.size(); ++evIndex) {
        const auto& ev = svEvidence[evIndex];
        if (ev.overlapID >= nCands) continue;
        perm[ovCursor[ev.overlapID]++] = evIndex;
    }

    auto countInformativeAllelesInOverlap = [&](size_t c) -> uint32_t {
        const size_t begin = ovOffsets[c];
        const size_t end = ovOffsets[c + 1];
        uint32_t o = 0;
        for (size_t k = begin; k < end; ++k) {
            const auto& ev = svEvidence[perm[k]];
            const uint32_t row = ev.overlapSite;
            if (row == invalid<uint32_t> || row >= nSites) continue;
            const auto& s = svStats[row];
            if (s.occ_0 < 2 || s.occ_1 < 2) continue;
            if (s.occ_0 >= s_hap_cov && s.occ_1 >= infor_cov) ++o;
        }
        return o;
    };

    /*
    Sort overlaps by number of informative SV alleles (descending).
    */
    auto& overlapSort = scratch.supportBits;
    overlapSort.clear();
    overlapSort.reserve(nCands);
    for (size_t c = 0; c < nCands; ++c) {
        /*
        Only overlaps that remain “cis” after the SNP stage can participate here: SV marking can
        only further filter overlaps; it does not undo SNP-based decisions.
        */
        if (candidates[c].is_match != 1) continue;
        const uint32_t o = countInformativeAllelesInOverlap(c);
        if (o == 0) continue;
        const uint64_t key = (uint64_t(uint32_t(-1) - o) << 32) | uint64_t(uint32_t(c));
        overlapSort.push_back(key);
    }
    std::sort(overlapSort.begin(), overlapSort.end());

    /*
    Decrement occ_0 for “not real alleles” contributed by a trans overlap.

    For SVs, we infer “ref-support” as: the overlap covers the site in query coordinates but does
    not have an SV event at that site (no svEvidence entry). This mirrors hifiasm’s adjustment for
    trans overlaps, which can reduce which rows remain informative after seeding.
    */
    auto decrementOcc0ForTransOverlap = [&](size_t c) {
        const auto& cand = candidates[c];
        const uint32_t qs = (uint32_t)cand.qs;
        const uint32_t qe = (uint32_t)cand.qe;
        if (qe <= qs) return;

        /*
        Collect the set of SV sites where this overlap has an SV event (sorted unique).
        */
        auto& mismatchSites = scratch.uniqueSites;
        mismatchSites.clear();
        {
            const size_t begin = ovOffsets[c];
            const size_t end = ovOffsets[c + 1];
            mismatchSites.reserve(end - begin);
            for (size_t k = begin; k < end; ++k) {
                mismatchSites.push_back(svEvidence[perm[k]].site);
            }
            std::sort(mismatchSites.begin(), mismatchSites.end());
            mismatchSites.erase(std::unique(mismatchSites.begin(), mismatchSites.end()), mismatchSites.end());
        }

        /*
        svStats are emitted in increasing site order; scan the overlap’s covered range [qs, qe).

        Any covered SV site that does not appear in mismatchSites gets its occ_0 decremented
        (clamped to at least 1 to preserve the implicit query-read contribution).
        */
        auto it = std::lower_bound(
            svStats.begin(), svStats.end(), qs,
            [](const SnpStats& s, uint32_t value) { return s.site < value; }
        );
        size_t m = 0;
        for (; it != svStats.end(); ++it) {
            const uint32_t site = it->site;
            if (site >= qe) break;
            while (m < mismatchSites.size() && mismatchSites[m] < site) ++m;
            if (m < mismatchSites.size() && mismatchSites[m] == site) continue;
            uint32_t& occ0 = it->occ_0;
            occ0 = (occ0 > 1) ? (occ0 - 1) : 1U;
        }
    };

    /*
    Pass 1: seed trans overlaps, promote SV sites, and apply occ_0 decrements.
    */
    for (uint64_t key : overlapSort) {
        const size_t c = size_t(uint32_t(key));
        const uint32_t o = countInformativeAllelesInOverlap(c);
        if (o == 0) continue;
        if (candidates[c].is_match != 1) continue;

        /* Mark overlap as trans/filtered. */
        candidates[c].is_match = 2;

        /* Promote SV sites supported by this overlap. */
        const size_t begin = ovOffsets[c];
        const size_t end = ovOffsets[c + 1];
        for (size_t k = begin; k < end; ++k) {
            const uint32_t row = svEvidence[perm[k]].overlapSite;
            if (row == invalid<uint32_t> || row >= nSites) continue;
            svStats[row].score = 1;
        }

        /* Remove ref-support at other SV sites contributed by this trans overlap. */
        decrementOcc0ForTransOverlap(c);
    }

    /*
    Pass 2: propagate trans status.

    If a remaining cis overlap carries an SV at any promoted row (score==1), it becomes trans too.
    */
    for (uint64_t key : overlapSort) {
        const size_t c = size_t(uint32_t(key));
        if (candidates[c].is_match != 1) continue;
        const size_t begin = ovOffsets[c];
        const size_t end = ovOffsets[c + 1];
        uint32_t o = 0;
        for (size_t k = begin; k < end; ++k) {
            const uint32_t row = svEvidence[perm[k]].overlapSite;
            if (row == invalid<uint32_t> || row >= nSites) continue;
            const auto& s = svStats[row];
            if (s.occ_0 < 2 || s.occ_1 < 2) continue;
            if (s.score == 1) ++o;
        }
        if (o > 0) candidates[c].is_match = 2;
    }

    /*
    Pass 3: reset block.

    Rows still supported only by overlaps remaining cis are reset to score==-1 (hifiasm staging).
    */
    for (uint64_t key : overlapSort) {
        const size_t c = size_t(uint32_t(key));
        if (candidates[c].is_match != 1) continue;
        const size_t begin = ovOffsets[c];
        const size_t end = ovOffsets[c + 1];
        for (size_t k = begin; k < end; ++k) {
            const uint32_t row = svEvidence[perm[k]].overlapSite;
            if (row == invalid<uint32_t> || row >= nSites) continue;
            svStats[row].score = -1;
        }
    }
}



/*
performHifiasmECParity: round-1 parity-EC filtering across all reads.

Reads are processed in parallel. For each query read, we:
  - gather its overlap candidates and normalize coordinates
  - detect candidate SNP sites in query coordinates
  - run DP chaining to select consistent informative SNP rows
  - mark overlaps trans/cis with hifiasm-style closure on SNP evidence
  - compact validated SNP rows and remap evidence
  - detect SV (large indel) sites and apply an SV analogue of overlap marking
  - decide which overlaps to keep in the read graph and update alignment flags

Timings are aggregated per phase (summed over threads) to guide optimization work.
*/
void Assembler::performHifiasmECParity(uint64_t threadCount)
{ 
    cout << timestamp << "=== Hifiasm Parity EC Pipeline (Round 1) ===" << endl;

    const uint64_t readCount = reads->readCount();
    const auto tBeginAll = steady_clock::now();

    static std::mutex svFlipPrintMutex;
        
    /*
    Byte-addressable keep vector (thread-safe for independent writes).

    The current pipeline applies keep/deletion decisions directly to per-side deletion flags in
    alignmentData, but this vector is kept for compatibility with earlier parity variants.
    */
    vector<uint8_t> keepAlignment(alignmentData.size(), 0); 
    
    /* Parallel loop over reads (static block scheduling). */
    vector<thread> threads;
    uint64_t chunkSize = readCount / threadCount;
    if(chunkSize == 0) chunkSize = 1;

    struct alignas(64) PhaseTiming {
        double gatherCandidates = 0.;
        double snpDetect = 0.;
        double dp = 0.;
        double dpBuildAltAnyBits = 0.;
        double dpBuildRefBits = 0.;
        double dpBuildHpBits = 0.;
        double dpTransitions = 0.;
        double dpExtractPaths = 0.;
        double validateSnv = 0.;
        double compact = 0.;
        double svDetect = 0.;
        double validateSv = 0.;
        double finalizeFlags = 0.;
        uint64_t readsVisited = 0;
        uint64_t readsWithAlignments = 0;
        uint64_t readsWithCandidates = 0;
    };
    vector<PhaseTiming> timings(threadCount);

    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            uint64_t start = t * chunkSize;
            uint64_t end = (t == threadCount - 1) ? readCount : (t + 1) * chunkSize;

            /* Thread-local scratchpad to eliminate per-read allocations. */
            HifiasmECScratchPad scratch;
            PhaseTiming& timing = timings[t];

            for(uint64_t readId = start; readId < end; readId++) {
                timing.readsVisited++;
                /*
                Process query read in forward orientation and normalize every candidate overlap so:
                  - candidate.qs/qe are query-forward coordinates
                  - candidate.ts/te are target-forward coordinates (even if the alignment is reverse)
                */
                uint32_t strand = 0;
                OrientedReadId orientedReadId(dinara::ReadId(readId), strand);
                if(orientedReadId.getValue() >= alignmentTable.size()) continue;

                const auto& alignments = alignmentTable[orientedReadId.getValue()];
                if(alignments.empty()) continue;
                timing.readsWithAlignments++;

	                    /* Clear scratchpad for this query read. */
	                    scratch.clear();
	                    auto& candidates = scratch.candidates;
	                    candidates.reserve(alignments.size());
                    
                    const auto tGatherBegin = steady_clock::now();
                    for(uint32_t alignmentId : alignments) {
                        const auto& thisAlignmentData = alignmentData[alignmentId];
                        
	                        /* Skip alignments already deleted by earlier graph filtering. */
	                        if(thisAlignmentData.isDeleted()) continue;

                        CandidateEC candidate;
                        candidate.alignmentId = alignmentId;

	                        /*
	                        Use marker-based (non-extended) coordinates for coverage-sensitive logic.

	                        This avoids counting unaligned read tips (important for contained reads) and
	                        aligns with hifiasm’s intended definition of “covered region”.
	                        */
	                        const OrientedReadId queryOriented(ReadId(readId), 0);
                        OrientedReadId o0(thisAlignmentData.readIds[0], 0);
                        OrientedReadId o1(thisAlignmentData.readIds[1], thisAlignmentData.isSameStrand ? 0 : 1);
                        AlignmentInfo orientedInfo = thisAlignmentData.info;

                        if (o0.getReadId() != queryOriented.getReadId()) {
                            std::swap(o0, o1);
                            orientedInfo.swap();
                        }
                        DINARA_ASSERT(o0.getReadId() == queryOriented.getReadId());

                        if (o0.getStrand() != queryOriented.getStrand()) {
                            o0.flipStrand();
                            o1.flipStrand();
                            orientedInfo.reverseComplement();
                        }
                        DINARA_ASSERT(o0 == queryOriented);

	                        /* Compute aligned marker bounds in query-forward coordinates. */
                        uint32_t qsCore = thisAlignmentData.qs;
                        uint32_t qeCore = thisAlignmentData.qe;
                        uint32_t tsCoreOriented = thisAlignmentData.ts;
                        uint32_t teCoreOriented = thisAlignmentData.te;
                        const uint32_t kmerLen = assemblerInfo.isOpen ? uint32_t(assemblerInfo->k) : 0U;
                        if (markers && assemblerInfo.isOpen) {
                            const auto m0 = (*markers)[o0.getValue()];
                            const auto m1 = (*markers)[o1.getValue()];
                            if (!m0.empty() && !m1.empty()) {
                                qsCore = m0[orientedInfo.data[0].firstOrdinal].position;
                                qeCore = m0[orientedInfo.data[0].lastOrdinal].position + kmerLen;
                                tsCoreOriented = m1[orientedInfo.data[1].firstOrdinal].position;
                                teCoreOriented = m1[orientedInfo.data[1].lastOrdinal].position + kmerLen;
                            }
                        }

	                        /*
	                        Convert target bounds into the target’s forward coordinates.

	                        This ensures that later target-base extraction is a simple forward-indexed
	                        access regardless of alignment orientation.
	                        */
                        const uint32_t tLen = uint32_t(reads->getRead(o1.getReadId()).baseCount);
                        uint32_t tsFwd = tsCoreOriented;
                        uint32_t teFwd = teCoreOriented;
                        if (o1.getStrand() != 0) {
                            tsFwd = tLen - teCoreOriented;
                            teFwd = tLen - tsCoreOriented;
                        }

                        candidate.qs = qsCore;
                        candidate.qe = qeCore;
                        candidate.ts = tsFwd;
                        candidate.te = teFwd;
                        candidate.targetId = uint32_t(o1.getReadId());
                        candidate.isRev = (o1.getStrand() != 0);
                        candidate.is_match = 1;
                        candidate.strong = 0;

                        candidates.push_back(candidate);
                    }
                    timing.gatherCandidates += seconds(steady_clock::now() - tGatherBegin);

                    if(candidates.empty()) continue;
                    timing.readsWithCandidates++;

                    /*
                    Stage 1: SNP site detection on the query coordinate space.
                    */
                    const auto tSnpBegin = steady_clock::now();
	                    detectHetSites(*this, *reads, readId, alignmentData, scratch);
	                    timing.snpDetect += seconds(steady_clock::now() - tSnpBegin);
	                    
	                    /*
	                    Stage 2: DP chaining on SNP rows.
	                    */
	                    const auto tDpBegin = steady_clock::now();
	                    RphaseDpTiming dpDetail;
	                    gen_rphase_dp(*this, scratch, &dpDetail);
                    timing.dp += seconds(steady_clock::now() - tDpBegin);
                    timing.dpBuildAltAnyBits += dpDetail.buildAltAnyBits;
                    timing.dpBuildRefBits += dpDetail.buildRefBits;
                    timing.dpBuildHpBits += dpDetail.buildHpBits;
                    timing.dpTransitions += dpDetail.transitions;
                    timing.dpExtractPaths += dpDetail.extractPaths;

	                    /*
	                    Stage 3: SNP-based overlap marking/validation (hifiasm parity).
	                    */
                    const auto tValidateSnvBegin = steady_clock::now();
                    generate_haplotypes_naive_HiFi(*this, scratch);
                    timing.validateSnv += seconds(steady_clock::now() - tValidateSnvBegin);

                    const auto tCompactBegin = steady_clock::now();
                    compactPhasedSites(scratch);
                    timing.compact += seconds(steady_clock::now() - tCompactBegin);

	                    /*
	                    Stage 4: SV (large indel) detection + overlap marking.
	                    */
	                    const auto tSvDetectBegin = steady_clock::now();
                    detectSVSites(*this, *reads, readId, alignmentData, scratch);
                    timing.svDetect += seconds(steady_clock::now() - tSvDetectBegin);

                    /*
                    Capture overlap state before SV marking so we can report which overlaps are
                    newly flipped to trans by generate_haplotypes_sv.
                    */
                    auto& wasTransBeforeSv = scratch.covered;
                    wasTransBeforeSv.assign(candidates.size(), 0);
                    for (size_t c = 0; c < candidates.size(); ++c) {
                        wasTransBeforeSv[c] = (candidates[c].is_match == 2);
                    }

                    const auto tValidateSvBegin = steady_clock::now();
                    generate_haplotypes_sv(*this, scratch);
                    timing.validateSv += seconds(steady_clock::now() - tValidateSvBegin);

                    /*
                    Report overlaps newly flipped to trans due to SV evidence.

                    This is intentionally independent of SNP filtering: the message shows only
                    overlaps whose is_match changed from 1->2 in generate_haplotypes_sv.
                    */
                    size_t svFlipCount = 0;
                    for (size_t c = 0; c < candidates.size(); ++c) {
                        if (candidates[c].is_match == 2 && !wasTransBeforeSv[c]) ++svFlipCount;
                    }
                    if (svFlipCount) {
                        std::lock_guard<std::mutex> lock(svFlipPrintMutex);
                        cout << timestamp << "[SV] read " << readId << " flipped " << svFlipCount
                             << " overlaps to trans" << endl;
                        size_t printed = 0;
                        constexpr size_t maxToPrint = 32;
                        for (size_t c = 0; c < candidates.size() && printed < maxToPrint; ++c) {
                            if (candidates[c].is_match != 2 || wasTransBeforeSv[c]) continue;
                            const auto& cand = candidates[c];
                            cout << timestamp << "  target " << cand.targetId
                                 << " strand " << (cand.isRev ? 'R' : 'F')
                                 << " alignmentId " << cand.alignmentId << endl;
                            ++printed;
                        }
                        if (svFlipCount > printed) {
                            cout << timestamp << "  ... " << (svFlipCount - printed) << " more" << endl;
                        }
                    }

                    /*
                    Informative read heuristic (matches hifiasm intent).

                    If at least one SNP or SV site survives validation, we treat the read as
	                    informative and keep only cis overlaps. Otherwise, we keep all overlaps because
	                    there is no reliable heterozygous signal to separate cis/trans.
                    */
                    const bool isInformativeRead = !scratch.snpStats.empty() || !scratch.svStats.empty();

	                    /*
	                    Precompute which overlaps touch any informative site.

	                    This lets us set AlignmentData::coversHetSite in O(|evidence|) rather than an
	                    O(|overlaps| * |sites|) nested scan.
	                    */
	                    vector<uint8_t> candCoversInfo(candidates.size(), 0);
                    for (const auto& ev : scratch.hapEvidence) {
                        if (ev.overlapSite != invalid<uint32_t> && ev.overlapID < candidates.size()) {
                            candCoversInfo[ev.overlapID] = 1;
                        }
                    }
                    for (const auto& ev : scratch.svEvidence) {
                         if (ev.overlapID < candidates.size()) {
                             candCoversInfo[ev.overlapID] = 1;
                         }
                    }

                    const auto tFinalizeBegin = steady_clock::now();
                    for(size_t c = 0; c < candidates.size(); ++c) {
                        auto& cand = candidates[c];
                        auto& ad = alignmentData[cand.alignmentId];
                        
	                        bool keep = !isInformativeRead || (cand.is_match == 1);

                        if (keep) ad.info.isInReadGraph = 1;

	                        /*
	                        Set directional deletion reasons (phasing).

	                        We record EC/phasing decisions as a bit in the per-side delete reason masks.
	                        This ensures later filters (chemical/chimeric/hang/contained/...) can add
	                        additional reasons without being overwritten by phasing reruns or rescues.
	                        */
	                        if (keep) {
	                            ad.clearDeleteReasonsFromReadPerspective(ReadId(readId), AlignmentData::DeleteReasonPhase);
	                        } else {
	                            ad.addDeleteReasonsFromReadPerspective(ReadId(readId), AlignmentData::DeleteReasonPhase);
	                        }

	                        /*
	                        Update assembly flags.

	                        These are monotonic (only ever set), so concurrent writes are safe.
	                        */
	                        if (keep) ad.info.isInReadGraph = 1;
	                        if (candCoversInfo[c]) ad.coversHetSite = true;
	                    }
                    timing.finalizeFlags += seconds(steady_clock::now() - tFinalizeBegin);
                }
        });
    }


    for(auto& th : threads) th.join();

    PhaseTiming total;
    for(const auto& t : timings) {
        total.gatherCandidates += t.gatherCandidates;
        total.snpDetect += t.snpDetect;
        total.dp += t.dp;
        total.dpBuildAltAnyBits += t.dpBuildAltAnyBits;
        total.dpBuildRefBits += t.dpBuildRefBits;
        total.dpBuildHpBits += t.dpBuildHpBits;
        total.dpTransitions += t.dpTransitions;
        total.dpExtractPaths += t.dpExtractPaths;
        total.validateSnv += t.validateSnv;
        total.compact += t.compact;
        total.svDetect += t.svDetect;
        total.validateSv += t.validateSv;
        total.finalizeFlags += t.finalizeFlags;
        total.readsVisited += t.readsVisited;
        total.readsWithAlignments += t.readsWithAlignments;
        total.readsWithCandidates += t.readsWithCandidates;
    }

    const double tAll = seconds(steady_clock::now() - tBeginAll);
    cout << timestamp << "Parity EC Round 1 timings (sum over threads):" << endl;
    cout << timestamp << "  gather candidates: " << total.gatherCandidates << " s" << endl;
    cout << timestamp << "  SNP detect:        " << total.snpDetect << " s" << endl;
    cout << timestamp << "  DP phase:          " << total.dp << " s" << endl;
    cout << timestamp << "    build alt/any:    " << total.dpBuildAltAnyBits << " s" << endl;
    cout << timestamp << "    build ref:        " << total.dpBuildRefBits << " s" << endl;
    cout << timestamp << "    build HP:         " << total.dpBuildHpBits << " s" << endl;
    cout << timestamp << "    transitions:      " << total.dpTransitions << " s" << endl;
    cout << timestamp << "    extract paths:    " << total.dpExtractPaths << " s" << endl;
    cout << timestamp << "  SNV validate:      " << total.validateSnv << " s" << endl;
    cout << timestamp << "  compact sites:     " << total.compact << " s" << endl;
    cout << timestamp << "  SV detect:         " << total.svDetect << " s" << endl;
    cout << timestamp << "  SV validate:       " << total.validateSv << " s" << endl;
    cout << timestamp << "  finalize flags:    " << total.finalizeFlags << " s" << endl;
    cout << timestamp << "Parity EC Round 1 wall time: " << tAll << " s"
         << " (reads=" << total.readsVisited
         << ", withAlignments=" << total.readsWithAlignments
         << ", withCandidates=" << total.readsWithCandidates << ")" << endl;

    cout << timestamp << "Parity EC Round 1 Complete." << endl;
}

void Assembler::performHifiasmECFinalFilteringParity(uint64_t /* threadCount */)
{
    cout << timestamp << "=== Hifiasm Parity EC Final Filtering (ha_ec_ff) ===" << endl;
}
