/*
Dinara hifiasm-parity error-correction (Parity EC).

This translation unit implements a per-read pipeline modeled after hifiasm’s phasing-based
error-correction logic (Correct.cpp, rphase_hc + gen_rphase_dp + generate_haplotypes_*).

The implementation here intentionally avoids base-quality (QV) logic. “Parity” in this file
means we match hifiasm’s decision structure and state transitions using Dinara’s data model:
alignment-derived evidence streams instead of hifiasm’s overlap_region_alloc and internal
haplotype_evdience lists.

IMPORTANT: HIFIASM's ONT PATH ACTUALLY USES QV
  In hifiasm, the ONT mode (--ont) sets is_sc=1, which means base-quality arrays
  are loaded and the qv parameter passed to gen_rphase_dp is non-NULL. This causes
  gen_rphase_dp0_single_path to take the QV branch (using get_hq_value scoring).
  Dinara deliberately omits QV scoring and instead adapts the no-QV branch logic,
  supplemented with select hardening features borrowed from the QV branch (notably
  the +8bp consecutive-site demotion). See the detailed note at +8bp demotion below.

HIFIASM ONT EC CALL CHAIN (Correct.cpp, v0.25.0-r726):
  ecovlp.cpp:3301 -> rphase_hc(... std_bs=is_ont=1, dp=&chainDP, q8=&v8q ...)
    +-- Correct.cpp:20191  rphase_hc()        ONT entry point
         |-- push_info()                      SNP row emission (occ_0, occ_1, overlap_num)
         |-- gen_rphase_dp(... st_rate=0.05, st_max=2, qv=q8 ...)
         |    |-- Site QV filtering            (qv != NULL -> QV branch)
         |    |-- fill_incom()                 Fill incomplete overlap evidence
         |    +-- gen_rphase_dp0_single_path(... qual_a=qv->a ...)
         |         |-- DP transitions          comput_sc_rphase
         |         |-- Path extraction          Greedy by descending f[]
         |         +-- Path scoring             QV branch: get_hq_value + +8bp demotion
         |                                    no-QV branch: HP masking + cc threshold
         |                                    (Dinara uses no-QV + borrowed +8bp demotion)
         |-- Post-DP compaction               Keep score==1, recount occ_0/overlap_num
         |-- generate_haplotypes_naive_HiFi(... multi_check=0, st_rate=0.05, st_max=2 ...)
         |    |-- Adjacent site filter          Remove distance-1 sites
         |    |-- Two-pass trans-closure        Seed -> propagate via is_st_bs
         |    +-- Final trans marking           is_st_bs gated
         |-- rphase_lidel()                   (not implemented in Dinara)
         +-- generate_haplotypes_sv()         SV-based trans marking (no is_st_bs)

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
#include "hifiasmECInternals.hpp"
#include "timestamp.hpp"
#include <algorithm>
#include <atomic>
#include <cstdlib>
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

// CandidateEC, SnpStats, HaplotypeEvidence, SweepEvent, RawSV, RphaseDpTiming,
// and HifiasmECScratchPad are defined in hifiasmECInternals.hpp (included above).

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
================================================================================
detectHetSites: SNP Detection and Evidence Extraction (Stage 1 of EC Pipeline)
================================================================================

PURPOSE:
  Build candidate heterozygous SNP sites and per-overlap mismatch evidence for one
  query read. This is the first stage of hifiasm-parity error correction and
  corresponds to hifiasm's "site aggregation" and "push_info" steps.

INPUT:
  - queryReadId: The read being error-corrected
  - candidates: List of overlaps (already filtered and scored) covering this read
  - alignmentData: Container with alignment evidence (SNPs, indels) from overlap computation

OUTPUT (via scratch pad):
  - snpStats: Candidate SNP rows (one per alternative allele, multi-allelic modeled as
              multiple rows at the same query position)
  - hapEvidence: Per-overlap mismatch observations linked to SNP rows
  - insertionOffsets/insertionIntervals: Per-overlap "hole" intervals where the overlap
              has no aligned base (deletions in query coordinate sense)

ALGORITHM OVERVIEW:
  Step 1: Unpack query read into 2-bit array for fast base lookups
  Step 2: Extract mismatch evidence from alignment store and build hole intervals
  Step 3: Build sorted list of unique query positions with mismatches
  Step 4: Compute per-site coverage using sweep-line algorithm with difference arrays
  Step 5: Aggregate mismatch counts per base, apply filters, emit SNP rows

KEY HIFIASM PARITY DETAILS:
  - Only emit SNP rows when alt allele has >=2 supporting overlaps (push_info rule)
  - Apply strand-bias filters to both ref and alt alleles (is_st_bs equivalent)
  - Query read contributes +1 to ref support (occ_0 includes query itself)
  - Coverage computation subtracts "hole" intervals so gapped positions don't count
  - Multi-allelic sites create multiple rows, each with independent alt allele

PERFORMANCE NOTES:
  - Sweep-line coverage: O(nOverlaps + nSites) instead of O(nOverlaps * nSites)
  - Binary search for hole interval mapping: O(log nHoles) per site per overlap
  - Difference arrays avoid repeated iteration over overlap ranges

WHY THIS MATTERS:
  Accurate SNP detection is critical for downstream phasing (DP chaining) and
  transitive closure. Over-calling creates spurious trans overlaps; under-calling
  misses true heterozygosity. The filters here (minimum support, strand bias,
  coverage threshold) balance sensitivity and specificity to match hifiasm behavior.
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
            if (currentPos >= unpacked.size() || base == unpacked[currentPos]) {
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
    =============================================================================
    Step 4: Compute Coverage at Each Candidate Site Using Difference Arrays
    =============================================================================

    GOAL:
      For each unique query position in uniqueSites, compute:
        - siteTotalCov[i]: Total number of overlaps covering that position
        - siteFwdCov[i]: Number of forward-orientation overlaps covering that position

      These coverage values are later used to derive ref support counts:
        refCov = totalCov - totalMismatches

      Forward coverage is used for strand-bias filtering (is_st_bs equivalent).

    WHY DIFFERENCE ARRAYS:
      Naively checking every overlap against every site would be O(#overlaps * #sites).
      For a read with 100 overlaps and 500 sites, that's 50,000 tests.

      Instead, we use a difference array approach:
        1. Map each overlap's [qs, qe) range to site indices using binary search
        2. Mark "coverage starts" at index idxS: diffTotal[idxS]++
        3. Mark "coverage ends" at index idxE: diffTotal[idxE]--
        4. Prefix-sum the difference array once: O(#sites)

      This reduces complexity to O(#overlaps * log(#sites) + #sites).
      For 100 overlaps and 500 sites: ~1,200 operations instead of 50,000.

    HOLE INTERVAL SUBTRACTION:
      Overlaps can have "holes" (deletions in query coordinate sense) where there's
      no aligned target base. These positions should NOT contribute to coverage.

      We subtract hole intervals using the same difference array technique:
        - For each hole [begin, end) within an overlap
        - Map to site indices [hS, hE)
        - Subtract from coverage: diffTotal[hS]--, diffTotal[hE]++

      This ensures ref support is only counted at positions with actual aligned bases.

    PERFORMANCE IMPACT:
      This optimization is critical for reads with many overlaps. Without it,
      coverage computation would dominate runtime for high-coverage datasets.
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
    =============================================================================
    Step 5: Aggregate Mismatch Counts and Emit SNP Rows with Filtering
    =============================================================================

    GOAL:
      Transform per-site mismatch counts and coverage into validated SNP rows.
      Each row represents one candidate heterozygous site with:
        - Query position (site)
        - Reference base and one alternative base
        - Ref support count (occ_0) and alt support count (occ_1)
        - Strand-bias metadata (fwd_ref_cov)
        - Query-side homopolymer flag (is_homopolymer)

    ALGORITHM:
      For each unique query position:
        1. Count mismatch support per base (A/C/G/T) and per strand
        2. Derive ref support: refCov = totalCov - totalMismatches
        3. Apply filters (minimum coverage, minimum alt support, strand bias)
        4. Emit one SNP row per alternative base that passes filters
        5. Link mismatch evidence entries to emitted rows

    FILTERING STAGES:

      A. Minimum Alt Support (push_info rule):
         - Require at least 2 overlaps supporting each alt allele
         - This is the PRIMARY filter from hifiasm's push_info()
         - Prevents spurious single-mismatch noise from creating SNP rows

      B. Minimum Total Coverage:
         - Require totalCov >= 5
         - Pragmatic guard against low-quality sites dominated by noise
         - Matches hifiasm behavior and downstream overlap marking expectations

      C. Ref Strand-Bias Filter (is_st_bs equivalent):
         - When refCov > 2, check if ref support is extremely biased to one strand
         - Reject if >=95% forward or <=5% forward
         - Prevents overcalling sites where ref is strand-specific artifact

      D. Alt Strand-Bias Filter:
         - When misCount > 2, same bias check on alt allele
         - Ensures alt evidence isn't purely from one strand

    MULTI-ALLELIC MODELING:
      Sites with multiple alternative alleles (e.g., A->C and A->T) create
      MULTIPLE SNP rows, one per alt allele. Each row is independently:
        - Filtered for strand bias and minimum support
        - Scored in DP chaining (can link different alt alleles)
        - Subject to transitive closure decisions

    QUERY +1 CONTRIBUTION:
      occ_0 includes +1 for the query read itself:
        occ_0 = refCov + 1

      This mirrors hifiasm's modeling where the query contributes one unit of
      reference support at its own base. Without this, hom-ref sites would
      appear zero-coverage.

    HOMOPOLYMER ANNOTATION:
      The is_homopolymer flag marks sites in query-side HP runs. This is
      used later for robustness checks (singleton HP validation), but does
      NOT cause sites to be dropped at this stage.

    EVIDENCE LINKING:
      Each mismatch evidence entry in hapEvidence gets its overlapSite field
      set to the SNP row index it supports. If the mismatch base didn't pass
      filters (e.g., weak third allele), overlapSite stays invalid and the
      evidence is ignored by downstream stages.
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
================================================================================
gen_rphase_dp: DP Chaining and Phasing Validation (Stage 2 of EC Pipeline)
================================================================================

PURPOSE:
  Determine which SNP rows are mutually consistent and should be validated
  for error correction. This is the CORE of hifiasm-parity EC and directly
  corresponds to hifiasm's gen_rphase_dp + gen_rphase_dp0_single_path.

  The central question: "Which heterozygous sites can be reliably phased
  together based on consistent allele patterns across overlapping reads?"

INPUT (via scratch pad):
  - snpStats: Candidate SNP rows from detectHetSites
  - hapEvidence: Per-overlap mismatch observations
  - candidates: List of overlaps covering the query read
  - insertionOffsets/insertionIntervals: Hole intervals for coverage computation

OUTPUT (via scratch pad):
  - snpStats: COMPACTED to only DP-retained rows with dpScore=1
  - hapEvidence: COMPACTED and remapped to new row indices
  - Final score field reset to -1 (ready for transitive closure stage)

ALGORITHM OVERVIEW:
  Step 1: Build bit-matrices for allele states (Ref, Alt, AnyOther per SNP row)
  Step 2: Build Ref bitsets using sweep-line coverage computation
  Step 3: Run DP chaining with strict biallelic linkage test
  Step 4: Extract disjoint best paths and apply cc filtering
  Step 5: Compact SNP rows and evidence to only DP-retained entries

KEY CONCEPTS:

  A. BITSET REPRESENTATION:
     Overlaps are represented as 64-bit bitsets (nWords = ceil(nCands/64)).
     For each SNP row i, we store THREE bitsets:
       - Ref[i]: Overlaps supporting reference allele at this site
       - Alt[i]: Overlaps supporting alternative allele at this site
       - AnyOther[i]: Overlaps with mismatch but not the chosen alt base

     These bitsets enable O(nWords) linkage tests between any pair of rows.

  B. STRICT BIALLELIC LINKAGE (comput_sc_rphase_strict):
     Two SNP rows can link if:
       - At least one overlap supports ref-ref across both sites
       - At least one overlap supports alt-alt across both sites
       - No contradictions: overlaps with "other" at one site and known at
         the other are treated as invalid (third-allele indicator)

     Special case: overlaps with "other" at BOTH sites count as ref-ref support
     (hifiasm's "rareRef" pattern).

  C. DP CHAINING:
     Let f[v] = best chain score ending at validIndices[v]
     Let p[v] = predecessor in that chain

     For each SNP row v eligible for DP (occ_0 >= s_hap_cov, occ_1 >= infor_cov):
       f[v] = 1  // Base case: singleton path
       For each potential predecessor u (site[u] < site[v]):
         If linkage test passes and f[u] + linkScore > f[v]:
           f[v] = f[u] + linkScore
           p[v] = u

     Rows at the same query position (multi-allelic) don't link to each other.

  D. PATH EXTRACTION:
     Extract disjoint paths in descending f[] order:
       - Multi-site paths (length > 1) always accepted: score = 1
       - Singleton paths require:
           1. occ_0 >= cc (coverage cutoff)
           2. NOT HP-masked after target-side homopolymer discounting

  E. CC CUTOFF (Coverage Cutoff):
     Hifiasm defines cc = max(cut_bd, (hom_cov / n_hap) * cut_rate)
     Where:
       - hom_cov = coveragePeak (from k-mer distribution)
       - n_hap = 2 (diploid assumption)
       - cut_rate = 0.7 (70% of expected haploid coverage)
       - cut_bd = 6 (minimum absolute cutoff)

     This filters low-coverage singleton paths that might be noise.

  F. SINGLETON HP VALIDATION:
     For singleton paths ONLY, check if the site would remain informative
     after discounting HP-suspect observations on the target side:
       - Count overlaps where target position is in a homopolymer run
       - Subtract from occ_0 and occ_1 to get robust counts
       - Reject if robust counts fall below thresholds (s_hap_cov=3, infor_cov=3)

     This matches hifiasm's robustness check for isolated heterozygous sites.

HIFIASM PARITY DETAILS:
  - Uses same DP formulation as gen_rphase_dp0_single_path
  - Identical cc calculation (hom_cov / 2 * 0.7, min 6)
  - Same multi-allelic handling (rows at same position don't link)
  - Same score semantics: score=1 for retained, score=-1 for dropped
  - Compaction before transitive closure matches hifiasm's pipeline flow

PERFORMANCE NOTES:
  - Bitset operations: O(nWords) per linkage test, very fast
  - DP complexity: O(nValidRows^2 * nWords) worst case
    - In practice, much better due to same-site blocking
    - For 100 rows: ~5,000 comparisons * ~2 words = ~10K bitwise ops
  - Sweep-line coverage: O(nCands + nSites) instead of O(nCands * nSites)
  - Singleton HP validation: O(nCands) target read accesses per singleton
    - This is the only expensive part (disabled in main DP for performance)

WHY THIS MATTERS:
  DP chaining is what separates true heterozygosity from sequencing noise.
  Multi-site phased paths have strong evidence across multiple positions.
  Singleton paths require extra validation (cc threshold, HP check) to avoid
  false positives. Getting this right determines error correction accuracy.
*/
void gen_rphase_dp(
    Assembler& assembler,
    HifiasmECScratchPad& scratch,
    RphaseDpTiming* timing
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
    =============================================================================
    Step 3: DP Chaining and Per-Row Scoring (gen_rphase_dp parity)
    =============================================================================

    PURPOSE:
      Chain validated SNP rows into phased haplotype paths using dynamic programming,
      then assign per-row acceptance scores. This corresponds to hifiasm's
      gen_rphase_dp -> gen_rphase_dp0_single_path pipeline.

    AT THIS POINT WE HAVE:
      - snpStats[]: SNP rows from detectHetSites (push_info-like emission)
      - flatBits[]: Bit-matrices encoding per-overlap ref/alt/any allele observations
        that allow O(nWords) linkage tests between any pair of SNP rows

    WHAT THIS STEP DOES:
      1. Select eligible rows ("real alleles"): occ_0 >= s_hap_cov AND occ_1 >= infor_cov,
         with is_st_bs strand-bias filtering for ONT
      2. Run DP transitions using comput_sc_rphase_strict (bit-parallel linkage)
      3. Extract disjoint best paths via greedy extraction by descending f[]
      4. Score paths: multi-site auto-accept, singleton HP+cc validation
      5. Consecutive-position contamination check (from QV branch): if any node in a
         run of adjacent-bp sites was rejected, reject the entire run
      6. Compact to score==1 rows and reset score semantics for trans-closure

    HIFIASM NO-QV vs QV BRANCH:
      In hifiasm, gen_rphase_dp0_single_path has two branches controlled by the
      qual_a parameter (Correct.cpp ~line 9472):

      if(!qual_a) {   <-- NO-QV BRANCH (Dinara's primary model)
        Multi-site chains: auto-accept (plus=1)
        Singletons: require !is_hpc_vec && occ_0>=cc
        No +8bp demotion check

      } else {        <-- QV BRANCH (hifiasm's actual ONT path, since is_sc=1 -> qv!=NULL)
        Uses get_hq_value() to split base-quality into b0l/b0h, b1l/b1h
        Has +8bp consecutive-site demotion (krn=1 if all gaps <= 8bp)
        Multi-site chains: per-node QV scoring
        Singletons/demoted: stricter QV + 70% confidence + cc threshold
        Consecutive-position contamination check on QV-scored runs

      Dinara uses the no-QV branch structure but BORROWS two elements from
      the QV branch as hardening measures:
        1. +8bp consecutive-site demotion (see detailed note below)
        2. Consecutive-position contamination check: if any node in a run of
           adjacent-position sites (gap==1bp) was rejected, all nodes in the
           run are forced to rejected. See detailed note below.
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
    /*
    ONT strand-bias constants for is_st_bs filtering (hifiasm parity).

    Hifiasm's is_st_bs macro rejects sites where ref support is extremely strand-biased:
      #define is_st_bs(s, rr, mm) (((mm) != ((uint64_t)-1)) &&
            (((s).overlap_num + mm) >= ((s).occ_0)) &&
            ((((s).occ_0*(rr) + (s).overlap_num)) >= ((s).occ_0)))

    For ONT (std_bs=1): st_rate=0.05, st_max=2
    For HiFi (std_bs=0): st_rate=0.0, st_max=(uint64_t)-1 (disabled)

    Dinara maps: fwd_ref_cov ↔ hifiasm overlap_num (forward-strand ref count + 1 for query)
    */
    constexpr double st_rate = 0.05;
    constexpr uint64_t st_max = 2;

    auto isStrandBiased = [](const SnpStats& s, double stRate, uint64_t stMax) -> bool {
        if (stMax == UINT64_MAX) return false;
        if (uint64_t(s.fwd_ref_cov) + stMax < uint64_t(s.occ_0)) return false;
        if (double(s.occ_0) * stRate + double(s.fwd_ref_cov) < double(s.occ_0)) return false;
        return true;
    };

    for (size_t i = 0; i < nSites; ++i) {
        if (snpStats[i].occ_0 < 2 || snpStats[i].occ_1 < 2) { snpStats[i].score = -1; continue; }
        if (isStrandBiased(snpStats[i], st_rate, st_max)) { snpStats[i].score = -1; continue; }
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
    =============================================================================
    Step 4: Extract Disjoint Best Paths from DP Predecessor Graph
    =============================================================================

    PURPOSE:
      Extract the best set of non-overlapping (disjoint) paths from the DP
      predecessor graph built in Step 3. Each path represents a phased set
      of heterozygous sites that will be validated for error correction.

    ALGORITHM:
      1. Sort all DP nodes by their best chain score f[v] in descending order
      2. Greedily extract paths starting from highest-scoring unvisited nodes:
         a) For each node v with f[v] >= 1 and not yet visited:
            - Walk backwards through predecessors: v → p[v] → p[p[v]] → ...
            - Mark all nodes in this path as visited
            - Apply scoring rules to determine if path is accepted (plus = 1 or -1)
         b) Assign final score to all nodes in the path based on plus value

    PATH SCORING RULES (matches hifiasm's no-QV branch):
      - Multi-site paths (length > 1): ALWAYS accepted (plus = 1)
      - Singleton paths (length == 1): Conditional acceptance based on:
          1. Coverage threshold: occ_0 >= cc
          2. HP robustness: Not HP-masked after target-side discounting
        If both conditions met: plus = 1, otherwise plus = -1

      Final assignment:
        For each node v in path:
          score = (occ_0 >= cc) ? plus : -1

      This means even within an accepted path, individual nodes can be rejected
      if they don't meet the cc coverage threshold.

    WHY DISJOINT PATHS:
      Each SNP row can only belong to ONE path. Once a node is visited and
      assigned to a path, it cannot be used in another path. This prevents:
        - Double-counting evidence from the same site
        - Conflicting phasing assignments
        - Overlapping paths that would create inconsistent haplotypes

    GREEDY EXTRACTION ORDER:
      Processing in descending f[] order ensures that:
        - Longest, highest-scoring chains are extracted first
        - Remaining nodes may form shorter paths or singletons
        - Maximum total phasing information is preserved

    HIFIASM PARITY:
      This exactly matches hifiasm's gen_rphase_dp0_single_path extraction:
        - Same greedy order (descending f[])
        - Same disjoint constraint (rowIsValid marking)
        - Same multi-site vs singleton scoring rules
        - Same cc threshold application per node

    PERFORMANCE:
      - Sorting: O(nV log nV) where nV = number of DP-eligible rows
      - Path extraction: O(nV) total (each node visited once)
      - For typical read: nV ~ 50-200 rows → negligible cost
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
        =====================================================================================
        SINGLETON PATH VALIDATION WITH TARGET-SIDE HOMOPOLYMER MASKING
        =====================================================================================

        PURPOSE:
          Determine the "plus" score for this extracted DP path, which controls whether
          the SNP sites in the path are retained for error correction.

        HIFIASM PARITY LOGIC (no-QV branch):
          1. Multi-site paths (length > 1): ALWAYS accepted (plus = 1)
             - Strong evidence from phasing across multiple sites
             - Unlikely to be sequencing noise

          2. Singleton paths (length == 1): Require BOTH conditions:
             a) Coverage threshold: occ_0 >= cc (typically 6-12 based on coverage peak)
             b) HP robustness: Site remains informative after discounting HP-suspect observations

        WHY SINGLETON PATHS NEED EXTRA VALIDATION:
          Multi-site paths have mutual support across positions, making them robust.
          Singleton paths are isolated heterozygous sites with no phasing context,
          making them vulnerable to:
            - Sequencing errors in homopolymer runs
            - Low-coverage noise
            - Systematic sequencing artifacts

        ALGORITHM:
          For each singleton path:
            Step A: Check coverage threshold (occ_0 >= cc)
            Step B: Count HP-suspect observations on TARGET side:
                    - For each overlap covering this site
                    - Map query coordinate to target coordinate
                    - Check if target position is in a homopolymer run (hpc_mask_ff)
                    - Count HP-suspect overlaps separately for ref and alt alleles
            Step C: Compute robust counts after HP discounting:
                    - n0_robust = occ_0 - hpRef
                    - n1_robust = occ_1 - hpAlt
            Step D: Reject if robust counts fall below thresholds:
                    - Must have n0_robust >= 2 AND n1_robust >= 2 (minimum evidence)
                    - Must have n0_robust >= s_hap_cov=3 AND n1_robust >= infor_cov=3
            Step E: Accept only if NOT HP-masked (plus = 1)

        TARGET-SIDE HP MASKING EXPLAINED:
          Homopolymer runs (e.g., AAAAA) are error-prone in sequencing. When the
          TARGET read has a HP run at the aligned position, that observation is
          considered unreliable (HP-suspect) even if the overlap supports ref or alt.

          Example:
            Query:  ...GCTA...  (query position 100)
            Target: ...AAAAAAA... (target position 500, inside HP run)
            Even if this overlap shows a mismatch at query position 100, we discount
            it because the target position is HP-suspect.

        COORDINATE MAPPING (CRITICAL FOR CORRECTNESS):
          Mapping query coordinate to target coordinate requires careful arithmetic:

          1. Compute query offset: queryOffset = querySite - cand.qs
          2. Map to target coordinate:
             - Forward overlap: targetSite = cand.ts + queryOffset
             - Reverse overlap: targetSite = cand.te - queryOffset - 1

          OVERFLOW/UNDERFLOW PROTECTION (Bug fix 2026-Feb-09):
            Original buggy code:
              const uint32_t targetSite = cand.isRev ?
                  (uint32_t)(cand.te - (querySite - cand.qs) - 1) :
                  (uint32_t)(cand.ts + (querySite - cand.qs));
              if (targetSite >= view.baseCount) continue;

            Problem: Cast to uint32_t BEFORE bounds check can wrap large uint64_t values:
              - If cand.ts + queryOffset > 2^32, cast wraps to small value
              - If cand.te <= queryOffset, subtraction underflows, cast wraps
              - Wrapped value passes bounds check, causes out-of-bounds access → SEGFAULT

            Fix: Perform arithmetic and bounds checking in uint64_t, THEN cast:
              const uint64_t queryOffset = uint64_t(querySite) - cand.qs;
              uint64_t targetSite64;
              if (cand.isRev) {
                  if (cand.te <= queryOffset) continue;  // Underflow check
                  targetSite64 = cand.te - queryOffset - 1;
              } else {
                  targetSite64 = cand.ts + queryOffset;
              }
              if (targetSite64 >= view.baseCount) continue;  // Bounds check in uint64_t
              const uint32_t targetSite = (uint32_t)targetSite64;  // Safe cast

        GAP CHECKING:
          Overlaps can have "holes" (deletions in query coordinate sense) where there's
          no aligned target base. The isGappedAtSite helper uses binary search on the
          sorted gap intervals to check if the query position falls in a hole.

          If gapped, the overlap is skipped for HP masking (no aligned base to check).

        PERFORMANCE NOTES:
          This is the ONLY place in the DP pipeline where target-side HP masking happens.
          The main DP loop (useOverlapHpMask flag at line 1529) has it DISABLED for
          performance, but singleton validation needs this robustness check.

          Cost: O(nCands * nSingletons) target read accesses, but typically only a few
          singleton paths exist, so total cost is acceptable (~1-5% overhead).

        WHY THIS MATTERS:
          Singleton HP validation prevents false positive heterozygous calls in HP regions.
          Without this check, sequencing errors in HP runs would be called as SNPs and
          incorrectly mark overlaps as trans, damaging assembly contiguity.

          Multi-site paths don't need this check because phasing across positions provides
          strong evidence independent of HP-specific artifacts.
        */
        /*
        Lambda: Check if overlap is gapped at this site (binary search).
        Defined once, reused for all singleton and demoted-chain node evaluations.
        */
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

        /*
        Lambda: Per-node singleton HP-masking check.
        Returns true if the node should be REJECTED (HP-masked).
        This is the no-QV equivalent of hifiasm's is_hpc_vec check:
        discount target-side HP-suspect observations and re-check thresholds.
        */
        auto isNodeHpMasked = [&](int siteIdx) -> bool {
            const uint32_t querySite = snpStats[siteIdx].site;
            const uint64_t* rowBits = &flatBits[siteIdx * 2 * nWords];

            uint64_t hpRef = 0, hpAlt = 0;
            for (size_t candIdx = 0; candIdx < nCands; ++candIdx) {
                const auto& cand = candidates[candIdx];
                if (querySite < cand.qs || querySite >= cand.qe) continue;
                if (isGappedAtSite(candIdx, querySite)) continue;

                const uint64_t queryOffset = uint64_t(querySite) - cand.qs;
                uint64_t targetSite64;
                if (cand.isRev) {
                    if (cand.te <= queryOffset) continue;
                    targetSite64 = cand.te - queryOffset - 1;
                } else {
                    targetSite64 = cand.ts + queryOffset;
                }

                const auto& view = assembler.getReads().getRead(cand.targetId);
                if (targetSite64 >= view.baseCount) continue;
                const uint32_t targetSite = (uint32_t)targetSite64;

                if (!hpc_mask_ff(view, (int64_t)view.baseCount, targetSite)) continue;

                const size_t w = candIdx >> 6;
                const uint64_t mask = 1ULL << (candIdx & 63);
                hpRef += __builtin_popcountll(rowBits[2 * w] & mask);
                hpAlt += __builtin_popcountll(rowBits[2 * w + 1] & mask);
            }

            const uint32_t n0_robust = (snpStats[siteIdx].occ_0 >= hpRef) ? (uint32_t)(snpStats[siteIdx].occ_0 - hpRef) : 0;
            const uint32_t n1_robust = (snpStats[siteIdx].occ_1 >= hpAlt) ? (uint32_t)(snpStats[siteIdx].occ_1 - hpAlt) : 0;

            return (n0_robust < 2 || n1_robust < 2 || n0_robust < s_hap_cov || n1_robust < infor_cov);
        };

        /*
        =====================================================================================
        +8BP CONSECUTIVE-SITE DEMOTION — DELIBERATELY BORROWED FROM HIFIASM'S QV BRANCH
        =====================================================================================

        WHAT THIS DOES:
          If ALL consecutive sites in a multi-site DP chain are within +8bp of each
          other, the entire chain is "demoted" to singleton-level validation. Instead
          of auto-accepting (plus=1), each node must independently pass the singleton
          HP-masking + cc-threshold checks.

        WHERE THIS COMES FROM IN HIFIASM:
          This logic exists ONLY in hifiasm's QV branch of gen_rphase_dp0_single_path
          (Correct.cpp ~line 9554-9555):

            for (i = 1; (i < rn) && ((a[res->a[rn0+i]].site+8) >= a[res->a[rn0+i-1]].site); i++);
            if(i >= rn) krn = 1;  // demote: treat as singleton-length chain

          The no-QV branch (Correct.cpp ~line 9480-9510) does NOT have this check.
          It simply auto-accepts all multi-site chains unconditionally (plus=1).

        WHY DINARA INCLUDES IT DESPITE USING THE NO-QV BRANCH:
          Hifiasm's real ONT path (is_sc=1 → qv!=NULL) takes the QV branch, which
          DOES have this demotion. Since Dinara omits per-base QV scoring entirely,
          borrowing the +8bp demotion provides a useful safeguard: it prevents
          tightly clustered sequencing noise (common in ONT homopolymer regions)
          from being auto-accepted as a multi-site chain, without requiring the
          full QV scoring machinery.

          In effect, this makes Dinara's no-QV adaptation MORE conservative than
          the pure no-QV branch and CLOSER to the behavior of hifiasm's actual ONT
          code path (the QV branch), which is the desired outcome.

        ALGORITHM:
          pathNodes are stored in extraction order (end→start of the chain), so
          pathNodes[0] has the highest site position and pathNodes[N-1] the lowest.
          We check: for every consecutive pair, is siteCur + 8 >= sitePrev?
          If ALL pairs satisfy this, every gap is ≤ 8bp → demote.
          If ANY pair has a gap > 8bp, the chain has real spread → keep as multi-site.

        DEMOTION EFFECT:
          Demoted chains go through the singleton validation path below, where each
          node is independently checked:
            - Coverage: occ_0 >= cc
            - HP robustness: !isNodeHpMasked (target-side HP discounting)
          This is much stricter than the auto-accept that multi-site chains get.
        */
        bool demotedToSingleton = false;
        if (pathNodes.size() > 1) {
            demotedToSingleton = true;
            for (size_t pi = 1; pi < pathNodes.size(); ++pi) {
                const uint32_t siteCur = snpStats[validIndices[pathNodes[pi]]].site;
                const uint32_t sitePrev = snpStats[validIndices[pathNodes[pi - 1]]].site;
                // pathNodes are end→start, so sitePrev >= siteCur; check gap <= 8
                if (siteCur + 8 < sitePrev) {
                    demotedToSingleton = false;
                    break;
                }
            }
        }

        if (pathNodes.size() > 1 && !demotedToSingleton) {
            /*
            TRUE MULTI-SITE PATH — AUTO-ACCEPT (hifiasm no-QV: plus=1)

            Hifiasm reference (Correct.cpp ~line 9482):
              if(rn > 1) { plus = 1; }   ← unconditional accept for multi-site

            Then per-node (Correct.cpp ~line 9509):
              if(a[...].occ_0 >= cc) { score = plus; } else { score = -1; }

            So even within an accepted multi-site chain, individual nodes that
            don't meet the cc coverage threshold are still rejected. This prevents
            weak nodes from riding the coattails of strong chains.
            */
            for (int v : pathNodes) {
                const int snpIdx = validIndices[v];
                snpStats[snpIdx].score = (snpStats[snpIdx].occ_0 >= cc) ? 1 : -1;
            }
        } else {
            /*
            SINGLETON OR DEMOTED CHAIN — PER-NODE VALIDATION

            Hifiasm no-QV reference (Correct.cpp ~line 9490-9493):
              if(rn == 1) {
                if((!is_hpc_vec(&(a[...]), ...)) && (a[...].occ_0 >= cc)) plus = 1;
              }

            Two conditions must BOTH be met for acceptance:
              1. Coverage threshold:  occ_0 >= cc
                 cc = max(hom_cov/n_hap * cut_rate, cut_bd)
                 Typically 6-12 depending on sequencing depth.

              2. HP robustness:  !isNodeHpMasked()
                 Discount target-side homopolymer-suspect observations and
                 re-check that robust counts still exceed thresholds
                 (n0_robust >= s_hap_cov=3 AND n1_robust >= infor_cov=3).
                 This is the no-QV equivalent of hifiasm's is_hpc_vec.

            For demoted chains (all sites within +8bp, see demotion note above),
            each node goes through this same validation independently. This is
            stricter than the QV branch's demotion behavior (which uses QV-based
            scoring), but provides a reasonable conservative approximation.
            */
            for (int v : pathNodes) {
                const int siteV = validIndices[v];
                int nodeScore = -1;
                if (snpStats[siteV].occ_0 >= cc) {
                    if (!isNodeHpMasked(siteV)) {
                        nodeScore = 1;
                    }
                }
                snpStats[siteV].score = nodeScore;
            }
        }

        /*
        =====================================================================================
        CONSECUTIVE-POSITION CONTAMINATION CHECK — BORROWED FROM HIFIASM'S QV BRANCH
        =====================================================================================

        WHAT THIS DOES:
          After scoring all nodes in the extracted path (whether multi-site, demoted,
          or singleton), scan for runs of consecutive genomic positions (gap == 1bp)
          within the path. If ANY node in a consecutive run was rejected (score == -1),
          ALL nodes in that run are forced to score == -1.

        WHERE THIS COMES FROM IN HIFIASM:
          This logic lives in the QV branch of gen_rphase_dp0_single_path
          (Correct.cpp ~lines 9579-9591), immediately after per-node QV scoring:

            for (i = j = 0; i < rn; i = j) {
                plus = a[res->a[rn0 + i]].score;
                for (j = i + 1; (j < rn) && ((site[i] - site[j]) == (j - i)); j++) {
                    if(a[res->a[rn0 + j]].score == -1) plus = -1;
                }
                if(plus == -1) { for (; i < j; i++) a[res->a[rn0+i]].score = plus; }
            }

          It applies ONLY in the QV branch (else clause of `if(!qual_a)`).
          The no-QV branch that Dinara primarily follows does NOT have this check.

        WHY DINARA INCLUDES IT:
          Hifiasm's real ONT path (is_sc=1 → qv!=NULL) takes the QV branch, which
          DOES apply this contamination check. Since Dinara adapts the no-QV
          branch structure, borrowing this check provides a useful safeguard:

          Adjacent-base-position SNPs (gap == 1bp) are highly suspicious for
          systematic ONT sequencing artifacts, especially around homopolymer
          boundaries. If one site in a consecutive run fails validation (HP
          masking, coverage threshold, or QV scoring in hifiasm's case), the
          neighboring sites at adjacent positions are very likely to be artifacts
          of the same underlying sequencing error.

          This check requires NO quality values — it is pure structural logic
          operating on position coordinates and previously-assigned scores.

        ALGORITHM:
          pathNodes are stored in extraction order (end→start), so pathNodes[0]
          has the highest site position and pathNodes[N-1] the lowest.
          Sites are in DESCENDING order within pathNodes.

          We scan pathNodes sequentially, grouping runs where:
            site[pathNodes[i]] - site[pathNodes[j]] == (j - i)
          i.e., each successive path node has a site exactly 1bp lower.

          Within each such consecutive run:
            - If ANY member has score == -1, set ALL members to score == -1
            - If all members have score == 1, leave them unchanged

          This matches hifiasm's exact semantics: the outer loop advances by
          runs (i = j after each run), and the inner loop extends j while the
          gap is exactly 1bp per step.

        EXAMPLE:
          Path nodes with sites: [105, 104, 103, 100, 50]
                                  ^^^^^^^^^^^  ^^^  ^^
                                  run1 (gap=1) |    singleton
                                               singleton

          If site 104 has score==-1 but 105 and 103 have score==1:
            → All three (105, 104, 103) get score = -1
          Site 100 and 50 are unaffected (not part of a consecutive run).

        APPLICABILITY:
          This check applies to ALL path types:
          - Multi-site chains: nodes scored via auto-accept + cc gate
          - Demoted chains (+8bp): nodes scored via singleton HP validation
          - True singletons: only 1 node, so no consecutive run possible
          For singletons, pathNodes.size()==1 means the outer loop runs once
          with a run of length 1, which is a no-op.

        PERFORMANCE:
          O(pathNodes.size()) per path — negligible cost.
        */
        if (pathNodes.size() > 1) {
            for (size_t pi = 0, pj = 0; pi < pathNodes.size(); pi = pj) {
                const uint32_t siteI = snpStats[validIndices[pathNodes[pi]]].site;
                bool anyRejected = (snpStats[validIndices[pathNodes[pi]]].score == -1);

                // Extend run: each next node must have site exactly 1bp lower
                for (pj = pi + 1; pj < pathNodes.size(); ++pj) {
                    const uint32_t siteJ = snpStats[validIndices[pathNodes[pj]]].site;
                    // pathNodes are descending, so siteI > siteJ for consecutive
                    if (siteI - siteJ != (uint32_t)(pj - pi)) break;
                    if (snpStats[validIndices[pathNodes[pj]]].score == -1) anyRejected = true;
                }

                // If any member of the consecutive run was rejected, reject all
                if (anyRejected) {
                    for (size_t pk = pi; pk < pj; ++pk) {
                        snpStats[validIndices[pathNodes[pk]]].score = -1;
                    }
                }
            }
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
================================================================================
generate_haplotypes_naive_HiFi: Transitive Closure (Stage 3 of EC Pipeline)
================================================================================

PURPOSE:
  Mark overlaps as "trans" (is_match=2) based on validated heterozygous sites
  from DP chaining. This is the final stage of hifiasm-parity error correction
  and directly corresponds to hifiasm's generate_haplotypes_naive_HiFi.

  The central question: "Which overlaps come from the alternate haplotype
  (trans) versus the same haplotype (cis)?"

TRANSITIVE CLOSURE CONCEPT:
  Once we've validated heterozygous sites through DP (Stage 2), we use those
  sites to partition overlaps into cis vs trans:
    - CIS overlaps (is_match=1): Match the query haplotype, should be retained
    - TRANS overlaps (is_match=2): Match the alternate haplotype, should be removed

  The "transitive" aspect: if an overlap carries alt alleles at multiple
  validated sites, it's very likely trans. Other overlaps covering those
  same sites are also likely trans (transitivity through shared sites).

INPUT (via scratch pad):
  - snpStats: DP-validated SNP rows (dpScore=1, compacted)
  - hapEvidence: Per-overlap mismatch observations linked to SNP rows
  - candidates: List of overlaps covering the query read

OUTPUT:
  - candidates[].is_match: Updated to 1 (cis) or 2 (trans)
  - snpStats[].occ_0: Decremented for "not real allele" adjustments
  - snpStats[].score: Set to inform_check value (1 or 0) for promoted sites

ALGORITHM OVERVIEW:
  Step 0: Drop adjacent SNP sites (distance 1 filter)
  Step 1: Sort evidence by overlap, count informative alleles per overlap
  Step 2: Seed trans overlaps from most-informative overlaps
  Step 3: Second pass: mark remaining overlaps hitting promoted sites as trans
  Step 4: Optional multi_check: rescue dense weak patterns

KEY CONCEPTS:

  A. INFORMATIVE ALLELE COUNTING:
     For each overlap, count how many validated SNP sites it carries alt alleles.
     An overlap with many alt alleles is likely trans; one with few/zero is likely cis.

     Thresholds (hifiasm defaults):
       - infor_cov = 3: Minimum alt alleles to be "informative enough" for seeding
       - s_hap_cov = 3: Minimum support required for sites to be "real alleles"

  B. SITE PROMOTION:
     A SNP site becomes "promoted" (score=1) if:
       - It's supported by at least s_hap_cov overlaps that have been marked trans
       - OR it's rescued in multi_check pass

     Promoted sites are considered highly reliable indicators of trans overlaps.

  C. TWO-PASS MARKING:
     Pass 1: Seed trans from overlaps with >= infor_cov informative alleles
             - These are the most obvious trans overlaps
             - Promote sites supported by these seeded trans overlaps

     Pass 2: Transitivity propagation
             - Any overlap hitting a promoted site is marked trans
             - This captures borderline overlaps that might have fewer alt alleles
             - Promotes additional sites supported by these newly marked trans overlaps

  D. NOT REAL ALLELE ADJUSTMENT:
     When an overlap is marked trans, it should no longer contribute ref support
     to OTHER sites (where it doesn't carry the alt allele).

     Reason: A trans overlap represents a different haplotype. At sites where it
     matches the query, that's NOT evidence for "ref" — it's just a coincidental
     match between two haplotypes.

     Implementation:
       For each overlap marked trans:
         For each SNP site covered by this overlap:
           If overlap does NOT carry the alt allele at this site:
             Decrement occ_0 (ref support count)

     This adjustment refines the allele support estimates for future iterations.

  E. MULTI_CHECK (Dense Weak Pattern Rescue):
     Some reads have many weakly-supported heterozygous sites clustered together.
     Each individual site might not pass infor_cov, but collectively they indicate
     trans overlaps.

     The multi_check pass looks for:
       - Dense clusters of sites within multi_check_distance=32 bp
       - Collective support from the same overlaps
       - If cluster strength >= up_num/up_den threshold (4/100 = 4%), promote sites

     This prevents loss of valid trans overlaps in low-coverage or noisy regions.

  F. ADJACENT SITE FILTERING (Step 0):
     Hifiasm removes SNP sites that are immediately adjacent (distance 1).
     Rationale: Adjacent mismatches often indicate local sequencing errors or
     alignment artifacts rather than true heterozygosity.

     This conservative filter prevents overcalling trans in noisy regions.

HIFIASM PARITY DETAILS:
  - Same two-pass seeding + transitivity algorithm
  - Same infor_cov=3, s_hap_cov=3 thresholds
  - Same multi_check_distance=32, up_num=4, up_den=100
  - Same "not real allele" decrement logic
  - Same adjacent site (distance 1) filtering

STATE TRANSITIONS:
  Initial:    All overlaps start as is_match=1 (cis)
  After Pass 1: High-confidence trans overlaps marked is_match=2
  After Pass 2: Borderline trans overlaps marked is_match=2 via transitivity
  After multi_check: Additional trans overlaps rescued from dense weak patterns

  Final: is_match values are used in performHifiasmECParity to:
    - Mark trans overlaps for deletion (DeleteReasonPhase)
    - Update alignment data deleteReasons0/1 fields
    - Drive assembly graph construction and error correction decisions

PERFORMANCE NOTES:
  - Evidence sorting: O(nEvidence log nEvidence), typically small
  - Two-pass marking: O(nEvidence), linear in evidence size
  - Multi_check: O(nSites * nCands), can be expensive for dense patterns
  - Overall: Fast compared to DP chaining, usually <10% of EC time

WHY THIS MATTERS:
  Transitive closure is what actually APPLIES the error correction decisions.
  Without accurate cis/trans classification:
    - False trans: Deletes valid overlaps, fragments assembly
    - False cis: Retains errors, reduces assembly accuracy

  The multi-pass approach with site promotion and multi_check provides robust
  classification even in challenging regions with low coverage or high noise.

DOWNSTREAM IMPACT:
  The is_match values set here control:
    - Which overlaps are deleted from the assembly graph
    - Which bases are corrected in the reads
    - Assembly contiguity and accuracy
    - Final contig N50 and correctness

  This is the final and most impactful stage of the EC pipeline.
*/
void generate_haplotypes_naive_HiFi(
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
    /*
    Hifiasm disables multi_check for the ONT path (multi_check=0 in rphase_hc when std_bs=1).
    Multi_check is a HiFi-specific recovery mechanism for dense weak heterozygous patterns;
    it is not used in the ONT error-correction pipeline.
    */
    constexpr bool enable_multi_check = false;

    /*
    ONT strand-bias constants for is_st_bs filtering (hifiasm parity).

    Hifiasm's is_st_bs macro rejects sites where ref support is extremely strand-biased:
      #define is_st_bs(s, rr, mm) (((mm) != ((uint64_t)-1)) &&
            (((s).overlap_num + mm) >= ((s).occ_0)) &&
            ((((s).occ_0*(rr) + (s).overlap_num)) >= ((s).occ_0)))

    For ONT (std_bs=1): st_rate=0.05, st_max=2
    Dinara maps: fwd_ref_cov <-> hifiasm overlap_num (forward-strand ref count + 1 for query)
    */
    constexpr double st_rate_tc = 0.05;
    constexpr uint64_t st_max_tc = 2;

    auto isStrandBiasedTC = [](const SnpStats& s, double stRate, uint64_t stMax) -> bool {
        if (stMax == UINT64_MAX) return false;
        if (uint64_t(s.fwd_ref_cov) + stMax < uint64_t(s.occ_0)) return false;
        if (double(s.occ_0) * stRate + double(s.fwd_ref_cov) < double(s.occ_0)) return false;
        return true;
    };

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
            if (isStrandBiasedTC(s, st_rate_tc, st_max_tc)) continue;
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
                /*
                Hifiasm also decrements overlap_num (fwd_ref_cov) when the trans overlap
                is forward-strand. This keeps strand-bias tracking accurate as overlaps
                are marked trans. (Correct.cpp ~line 9006)
                */
                if (!candidates[c].isRev) {
                    uint32_t& frc = snpStats[r].fwd_ref_cov;
                    frc = (frc > 1) ? (frc - 1) : 1U;
                }
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

    Any overlap that is still cis (is_match==1) but has a mismatch at a promoted site becomes
    trans, provided the site still passes basic thresholds and strand-bias filtering.
    Hifiasm applies is_st_bs here too (Correct.cpp ~line 9033).
    */
    for (size_t c = 0; c < nCands; ++c) {
        if (candidates[c].is_match == 2) continue;
        const size_t begin = ovOffsets[c];
        const size_t end = ovOffsets[c + 1];
        uint32_t o = 0;
        for (size_t k = begin; k < end; ++k) {
            const uint32_t row = hapEvidence[perm[k]].overlapSite;
            if (row == invalid<uint32_t> || row >= nSites) continue;
            const auto& s = snpStats[row];
            if (s.occ_0 < 2 || s.occ_1 < 2) continue;
            if (isStrandBiasedTC(s, st_rate_tc, st_max_tc)) continue;
            if (s.score == 1) ++o;
        }
        if (o > 0) candidates[c].is_match = 2;
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
    Step 7: Final trans marking (hifiasm Correct.cpp ~line 9085-9110).

    After multi_check (disabled for ONT) potentially marks additional rows score=1,
    we perform one final pass: any overlap that is still cis (is_match==1) but has a
    mismatch at a score==1 row becomes trans (is_match=2).

    Hifiasm reference (generate_haplotypes_naive_HiFi, Correct.cpp ~line 9093-9099):
      if(s->score == 1 && (!(s->occ_0 < 2 || s->occ_1 < 2)) && (!(is_st_bs((*s), st_rate, st_max)))) {
          overlap_list->list[ii].is_match = 2;
      }

    Note: this final marking has THREE conditions, all of which we replicate:
      1. score == 1          (site has been promoted by trans-closure)
      2. occ_0 >= 2 AND occ_1 >= 2  (minimum allele support)
      3. !is_st_bs           (not strand-biased — ONT-specific guard)

    Contrast with generate_haplotypes_sv's final marking (Correct.cpp ~line 9230-9237),
    which does NOT check is_st_bs — only checks score==1 && occ_0>=2 && occ_1>=2.
    This is correct because SV evidence is less susceptible to strand bias artifacts.
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
            if (isStrandBiasedTC(s, st_rate_tc, st_max_tc)) continue;
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

    // Single mutex for all cout output emitted by this EC function.
    // This prevents interleaving between [SV] and [EC-DBG] messages across threads.
    static std::mutex ecCoutMutex;

    struct EcDebugReadSummary {
        bool present = false;
        ReadId readId = invalid<ReadId>;
        uint32_t readLength = 0;
        uint32_t candidateCount = 0;
        vector<ReadId> targetReads; // unique, sorted
    };
    struct EcDebugPairOverlapSummary {
        bool present = false;
        ReadId readId = invalid<ReadId>;
        uint32_t qs = 0;
        uint32_t qe = 0;
        vector<uint32_t> informativeSnpSitesInOverlap; // unique, sorted, after hole subtraction
    };

    // Optional debug hook: set DINARA_EC_DEBUG_ALIGNMENT_ID to an alignmentId (uint32)
    // to print detailed accounting of informative-site counts for that overlap from each read's view.
    uint32_t ecDebugAlignmentId = invalid<uint32_t>;
    array<ReadId, 2> ecDebugReadIds{invalid<ReadId>, invalid<ReadId>};
    bool ecDebugByReadPair = false;
    array<ReadId, 2> ecDebugReadPair{invalid<ReadId>, invalid<ReadId>};
    array<EcDebugReadSummary, 2> ecDebugReadSummaries;
    array<EcDebugPairOverlapSummary, 2> ecDebugPairOverlapSummaries;
    bool ecDebugPairIsSameStrand = true;
    {
        auto parseEnvU32 = [](const char* name, uint32_t& out) -> bool {
            const char* s = std::getenv(name);
            if(!s || !*s) return false;
            char* end = nullptr;
            const unsigned long v = std::strtoul(s, &end, 10);
            if(!(end && end != s && *end == 0 && v <= std::numeric_limits<uint32_t>::max())) {
                return false;
            }
            out = uint32_t(v);
            return true;
        };

        // Optional: debug by read pair (matches any alignmentId for those two reads).
        // Set DINARA_EC_DEBUG_READ0 and DINARA_EC_DEBUG_READ1 (uint32).
        uint32_t d0 = invalid<uint32_t>;
        uint32_t d1 = invalid<uint32_t>;
        if(parseEnvU32("DINARA_EC_DEBUG_READ0", d0) && parseEnvU32("DINARA_EC_DEBUG_READ1", d1)) {
            const ReadId r0 = ReadId(std::min(d0, d1));
            const ReadId r1 = ReadId(std::max(d0, d1));
            ecDebugByReadPair = true;
            ecDebugReadPair = {r0, r1};
            ecDebugReadSummaries[0].readId = r0;
            ecDebugReadSummaries[1].readId = r1;
            ecDebugPairOverlapSummaries[0].readId = r0;
            ecDebugPairOverlapSummaries[1].readId = r1;
            std::lock_guard<std::mutex> lock(ecCoutMutex);
            cout << timestamp << "[EC-DBG] Enabled for readPair=(" << r0 << "," << r1 << ")" << endl;
        }

        const char* s = std::getenv("DINARA_EC_DEBUG_ALIGNMENT_ID");
        if(s && *s) {
            char* end = nullptr;
            const unsigned long v = std::strtoul(s, &end, 10);
            if(end && end != s && *end == 0 && v <= std::numeric_limits<uint32_t>::max()) {
                ecDebugAlignmentId = uint32_t(v);
                if(ecDebugAlignmentId < alignmentData.size()) {
                    ecDebugReadIds = alignmentData[ecDebugAlignmentId].readIds;
                    std::lock_guard<std::mutex> lock(ecCoutMutex);
                    cout << timestamp << "[EC-DBG] Enabled for alignmentId=" << ecDebugAlignmentId
                         << " reads=(" << ecDebugReadIds[0] << "," << ecDebugReadIds[1] << ")" << endl;
                } else {
                    std::lock_guard<std::mutex> lock(ecCoutMutex);
                    cout << timestamp << "[EC-DBG] Requested alignmentId=" << ecDebugAlignmentId
                         << " but alignmentData.size()=" << alignmentData.size() << " (ignoring)" << endl;
                    ecDebugAlignmentId = invalid<uint32_t>;
                }
            }
        }
    }
        
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
	                        
		                        /*
		                        Skip alignments deleted on either side. Only use overlaps kept by both
		                        reads (keptByBothSides).
		                        */
		                        if(!thisAlignmentData.keptByBothSides()) continue;

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
                        uint32_t qsCore = 0;
                        uint32_t qeCore = 0;
                        uint32_t tsCoreOriented = 0;
                        uint32_t teCoreOriented = 0;
                        bool coordinatesFromMarkers = false;
                        const ReadId queryReadId = ReadId(readId);

                        const uint32_t kmerLen = assemblerInfo.isOpen ? uint32_t(assemblerInfo->k) : 0U;
                        if (markers && assemblerInfo.isOpen) {
                            const auto m0 = (*markers)[o0.getValue()];
                            const auto m1 = (*markers)[o1.getValue()];
                            if (!m0.empty() && !m1.empty()) {
                                qsCore = m0[orientedInfo.data[0].firstOrdinal].position;
                                qeCore = m0[orientedInfo.data[0].lastOrdinal].position + kmerLen;
                                tsCoreOriented = m1[orientedInfo.data[1].firstOrdinal].position;
                                teCoreOriented = m1[orientedInfo.data[1].lastOrdinal].position + kmerLen;
                                coordinatesFromMarkers = true;
                            }
                        }
                        if(!coordinatesFromMarkers) {
                            // AlignmentData stores base coordinates in forward read coordinates:
                            //   - qs/qe for readIds[0]
                            //   - ts/te for readIds[1] (already converted to forward, even for reverse overlaps).
                            // When the query read is readIds[1], we must swap accordingly.
                            if(thisAlignmentData.readIds[0] == queryReadId) {
                                qsCore = thisAlignmentData.qs;
                                qeCore = thisAlignmentData.qe;
                                tsCoreOriented = thisAlignmentData.ts;
                                teCoreOriented = thisAlignmentData.te;
                            } else {
                                qsCore = thisAlignmentData.ts;
                                qeCore = thisAlignmentData.te;
                                tsCoreOriented = thisAlignmentData.qs;
                                teCoreOriented = thisAlignmentData.qe;
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
                        // Marker positions for strand 1 are in reverse-complement coordinates and must be
                        // converted back to forward read coordinates. Base coordinates stored in AlignmentData
                        // are already in forward read coordinates and must NOT be converted again.
                        if (coordinatesFromMarkers && o1.getStrand() != 0) {
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

                    // Debug: capture per-read overlap neighborhood for the requested read pair.
                    if(ecDebugByReadPair) {
                        const ReadId queryReadId = ReadId(readId);
                        int idx = -1;
                        if(queryReadId == ecDebugReadPair[0]) {
                            idx = 0;
                        } else if(queryReadId == ecDebugReadPair[1]) {
                            idx = 1;
                        }
                        if(idx >= 0) {
                            vector<ReadId> targets;
                            targets.reserve(candidates.size());
                            for(const auto& cand : candidates) {
                                targets.push_back(ReadId(cand.targetId));
                            }
                            std::sort(targets.begin(), targets.end());
                            targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
                            const uint32_t qLen = uint32_t(reads->getRead(queryReadId).baseCount);

                            std::lock_guard<std::mutex> lock(ecCoutMutex);
                            if(!ecDebugReadSummaries[idx].present) {
                                ecDebugReadSummaries[idx].present = true;
                                ecDebugReadSummaries[idx].readLength = qLen;
                                ecDebugReadSummaries[idx].candidateCount = uint32_t(candidates.size());
                                ecDebugReadSummaries[idx].targetReads = std::move(targets);
                            }
                        }
                    }

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

                    if(readId == 0 && strand == 0) {
                        std::lock_guard<std::mutex> lock(ecCoutMutex);
                        cout << timestamp << "[EC-DBG] Surviving SNP sites for read 0-0 after compact: "
                             << scratch.snpStats.size() << " row(s)" << endl;
                        if(scratch.snpStats.empty()) {
                            cout << timestamp << "[EC-DBG]   (none)" << endl;
                        } else {
                            for(size_t snpRow = 0; snpRow < scratch.snpStats.size(); ++snpRow) {
                                const auto& s = scratch.snpStats[snpRow];
                                cout << timestamp << "[EC-DBG]   row=" << snpRow
                                     << " site=" << s.site
                                     << " ref=" << s.refBase
                                     << " alt=" << s.altBase
                                     << " occ0=" << s.occ_0
                                     << " occ1=" << s.occ_1
                                     << " fwdRef=" << s.fwd_ref_cov
                                     << " hp=" << int(s.is_homopolymer)
                                     << " score=" << s.score
                                     << " dpScore=" << s.dpScore
                                     << endl;
                            }
                        }
                    }

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

//                     /*
//                     Report overlaps newly flipped to trans due to SV evidence.
//
//                     This is intentionally independent of SNP filtering: the message shows only
//                     overlaps whose is_match changed from 1->2 in generate_haplotypes_sv.
//                     */
//                     size_t svFlipCount = 0;
//                     for (size_t c = 0; c < candidates.size(); ++c) {
//                         if (candidates[c].is_match == 2 && !wasTransBeforeSv[c]) ++svFlipCount;
//                     }
//                     if (svFlipCount) {
//                         std::lock_guard<std::mutex> lock(ecCoutMutex);
//                         cout << timestamp << "[SV] read " << readId << " flipped " << svFlipCount
//                              << " overlaps to trans" << endl;
//                         size_t printed = 0;
//                         constexpr size_t maxToPrint = 32;
//                         for (size_t c = 0; c < candidates.size() && printed < maxToPrint; ++c) {
//                             if (candidates[c].is_match != 2 || wasTransBeforeSv[c]) continue;
//                             const auto& cand = candidates[c];
//                             cout << timestamp << "  target " << cand.targetId
//                                  << " strand " << (cand.isRev ? 'R' : 'F')
//                                  << " alignmentId " << cand.alignmentId << endl;
//                             ++printed;
//                         }
//                         if (svFlipCount > printed) {
//                             cout << timestamp << "  ... " << (svFlipCount - printed) << " more" << endl;
//                         }
//                     }

                    /*
                    Informative read heuristic (matches hifiasm intent).

                    If at least one SNP or SV site survives validation, we treat the read as
	                    informative and keep only cis overlaps. Otherwise, we keep all overlaps because
	                    there is no reliable heterozygous signal to separate cis/trans.
                    */
                    const bool isInformativeRead = !scratch.snpStats.empty() || !scratch.svStats.empty();

	                    /*
	                    Precompute how many informative sites each overlap covers (query view).

	                    We count unique SNP/SV site positions retained by the EC pipeline that fall
	                    inside the overlap's [qs,qe) interval. For SNP sites we also exclude query
	                    "hole" intervals (positions without an aligned target base) for that overlap.
	                    */
	                    vector<uint32_t> informativeSnpSites;
	                    vector<uint32_t> candInformativeSiteCount(candidates.size(), 0);

	                    informativeSnpSites.clear();
	                    if(!scratch.snpStats.empty()) {
	                        informativeSnpSites.reserve(scratch.snpStats.size());
	                        for(const auto& s : scratch.snpStats) {
	                            informativeSnpSites.push_back(s.site);
	                        }
                        std::sort(informativeSnpSites.begin(), informativeSnpSites.end());
                        informativeSnpSites.erase(
                            std::unique(informativeSnpSites.begin(), informativeSnpSites.end()),
                            informativeSnpSites.end());
                    }

	                    // Note: informativeHetSiteCount{0,1} are currently computed using SNP sites only.
	                    // SV sites (large indels) are still used for cis/trans decisions, but are not used
	                    // for overlap ranking, because SNP sites are easier to interpret and compare across reads.

                    auto countSitesInRange = [](const vector<uint32_t>& sites, uint32_t begin, uint32_t end) -> uint32_t {
                        if(end <= begin || sites.empty()) {
                            return 0;
                        }
                        auto it0 = std::lower_bound(sites.begin(), sites.end(), begin);
                        auto it1 = std::lower_bound(sites.begin(), sites.end(), end);
                        return uint32_t(it1 - it0);
                    };

                    const bool haveInsertionHoles =
                        (scratch.insertionOffsets.size() == candidates.size() + 1) &&
                        !scratch.insertionIntervals.empty();

	                    for(size_t c = 0; c < candidates.size(); ++c) {
	                        const auto& cand = candidates[c];
	                        const uint32_t qs = uint32_t(cand.qs);
	                        const uint32_t qe = uint32_t(cand.qe);
	                        if(qe <= qs) {
	                            continue;
	                        }

	                        uint32_t count = 0;
	                        uint32_t snpCountInRange = 0;
	                        uint32_t snpCountSubtracted = 0;
	                        // SNP sites: subtract query "holes" for this overlap.
	                        if(!informativeSnpSites.empty()) {
	                            snpCountInRange = countSitesInRange(informativeSnpSites, qs, qe);
	                            count += snpCountInRange;
	                            if(haveInsertionHoles) {
	                                const uint32_t offBegin = scratch.insertionOffsets[c];
	                                const uint32_t offEnd = scratch.insertionOffsets[c + 1];
	                                for(uint32_t ii = offBegin; ii < offEnd; ++ii) {
	                                    const auto& hole = scratch.insertionIntervals[ii];
	                                    const uint32_t a = std::max(qs, hole.begin);
	                                    const uint32_t b = std::min(qe, hole.end);
	                                    if(b > a) {
	                                        const uint32_t sub = countSitesInRange(informativeSnpSites, a, b);
	                                        snpCountSubtracted += sub;
	                                        count -= sub;
	                                    }
	                                }
	                            }
	                        }

	                        candInformativeSiteCount[c] = count;

                            // Debug: for the requested read pair, capture the exact SNP sites that contribute
                            // to the informative count for the overlap between the two reads (query view).
                            if(ecDebugByReadPair) {
                                const ReadId queryReadId = ReadId(readId);
                                int idx = -1;
                                if(queryReadId == ecDebugReadPair[0]) {
                                    idx = 0;
                                } else if(queryReadId == ecDebugReadPair[1]) {
                                    idx = 1;
                                }
                                if(idx >= 0 && cand.targetId == uint32_t(ecDebugReadPair[1 - idx])) {
                                    std::lock_guard<std::mutex> lock(ecCoutMutex);
                                    if(!ecDebugPairOverlapSummaries[idx].present) {
                                        vector<uint32_t> sitesInOverlap;
                                        if(!informativeSnpSites.empty() && qe > qs) {
                                            auto it0 = std::lower_bound(informativeSnpSites.begin(), informativeSnpSites.end(), qs);
                                            auto it1 = std::lower_bound(informativeSnpSites.begin(), informativeSnpSites.end(), qe);
                                            sitesInOverlap.assign(it0, it1);

                                            // Subtract query holes for this overlap (same as the counter).
                                            if(haveInsertionHoles) {
                                                const uint32_t offBegin = scratch.insertionOffsets[c];
                                                const uint32_t offEnd = scratch.insertionOffsets[c + 1];
                                                if(offBegin < offEnd && !sitesInOverlap.empty()) {
                                                    vector<uint8_t> keepMask(sitesInOverlap.size(), 1);
                                                    for(uint32_t ii = offBegin; ii < offEnd; ++ii) {
                                                        const auto& hole = scratch.insertionIntervals[ii];
                                                        const uint32_t a = std::max(qs, hole.begin);
                                                        const uint32_t b = std::min(qe, hole.end);
                                                        if(b <= a) continue;
                                                        auto h0 = std::lower_bound(sitesInOverlap.begin(), sitesInOverlap.end(), a);
                                                        auto h1 = std::lower_bound(sitesInOverlap.begin(), sitesInOverlap.end(), b);
                                                        for(auto it = h0; it != h1; ++it) {
                                                            keepMask[size_t(it - sitesInOverlap.begin())] = 0;
                                                        }
                                                    }
                                                    size_t write = 0;
                                                    for(size_t i = 0; i < sitesInOverlap.size(); ++i) {
                                                        if(keepMask[i]) sitesInOverlap[write++] = sitesInOverlap[i];
                                                    }
                                                    sitesInOverlap.resize(write);
                                                }
                                            }
                                        }

                                        ecDebugPairOverlapSummaries[idx].present = true;
                                        ecDebugPairOverlapSummaries[idx].qs = qs;
                                        ecDebugPairOverlapSummaries[idx].qe = qe;
                                        ecDebugPairOverlapSummaries[idx].informativeSnpSitesInOverlap = std::move(sitesInOverlap);
                                        ecDebugPairIsSameStrand = alignmentData[cand.alignmentId].isSameStrand;
                                    }
                                }
                            }

	                        // Debug: print exact accounting for a specific overlap alignmentId.
	                        const bool debugThisOverlap =
	                            (ecDebugAlignmentId != invalid<uint32_t> && cand.alignmentId == ecDebugAlignmentId) ||
	                            (ecDebugByReadPair &&
	                             alignmentData[cand.alignmentId].readIds[0] == ecDebugReadPair[0] &&
	                             alignmentData[cand.alignmentId].readIds[1] == ecDebugReadPair[1]);
	                        if(debugThisOverlap) {
	                            std::lock_guard<std::mutex> lock(ecCoutMutex);
	                            const ReadId queryReadId = ReadId(readId);
	                            const AlignmentData& adDbg = alignmentData[cand.alignmentId];
	                            const uint32_t qLen = uint32_t(reads->getRead(queryReadId).baseCount);
	                            cout << timestamp << "[EC-DBG] QueryRead=" << queryReadId
	                                 << " qLen=" << qLen
	                                 << " candidates=" << candidates.size()
	                                 << " overlapReads=(" << adDbg.readIds[0] << "," << adDbg.readIds[1] << ")"
	                                 << " qs=" << qs << " qe=" << qe << " span=" << (qe - qs)
	                                 << " isInformativeRead=" << (isInformativeRead ? 1 : 0)
	                                 << " cand.is_match=" << int(cand.is_match)
	                                 << " keep=" << ((!isInformativeRead || cand.is_match == 1) ? 1 : 0)
	                                 << endl;
	                            cout << timestamp << "[EC-DBG] SNP rows=" << scratch.snpStats.size()
	                                 << " SNP uniqueSites=" << informativeSnpSites.size()
	                                 << " SV rows=" << scratch.svStats.size()
	                                 << endl;
	                            cout << timestamp << "[EC-DBG] Count: snpInRange=" << snpCountInRange
	                                 << " snpSubtracted=" << snpCountSubtracted
	                                 << " snpFinal=" << (snpCountInRange - snpCountSubtracted)
	                                 << " total=" << count
	                                 << endl;

	                            // Print a few SNP site positions in-range (after DP+multi_check+compact+unique).
	                            if(!informativeSnpSites.empty()) {
	                                auto it0 = std::lower_bound(informativeSnpSites.begin(), informativeSnpSites.end(), qs);
	                                auto it1 = std::lower_bound(informativeSnpSites.begin(), informativeSnpSites.end(), qe);
	                                const size_t inRange = size_t(it1 - it0);
	                                cout << timestamp << "[EC-DBG] SNP sites in-range=" << inRange << " first=";
	                                size_t printed = 0;
	                                for(auto it = it0; it != it1 && printed < 20; ++it, ++printed) {
	                                    if(printed) cout << ",";
	                                    cout << *it;
	                                }
	                                if(inRange > printed) cout << ",...";
	                                cout << endl;
	                            }

	                            // Print hole intervals for this overlap and how many SNP sites each subtracts.
	                            if(haveInsertionHoles) {
	                                const uint32_t offBegin = scratch.insertionOffsets[c];
	                                const uint32_t offEnd = scratch.insertionOffsets[c + 1];
	                                cout << timestamp << "[EC-DBG] HoleIntervals count=" << (offEnd - offBegin) << " ";
	                                uint32_t sumSub = 0;
	                                for(uint32_t ii = offBegin; ii < offEnd; ++ii) {
	                                    const auto& hole = scratch.insertionIntervals[ii];
	                                    const uint32_t a = std::max(qs, hole.begin);
	                                    const uint32_t b = std::min(qe, hole.end);
	                                    if(b <= a) continue;
	                                    const uint32_t sub = informativeSnpSites.empty() ? 0 : countSitesInRange(informativeSnpSites, a, b);
	                                    sumSub += sub;
	                                    cout << "[" << a << "," << b << ")->" << sub << " ";
	                                }
	                                cout << "sumSub=" << sumSub << endl;
	                            }

	                            // Sanity check: debug recomputation must match the stored candidate count.
	                            DINARA_ASSERT(count == candInformativeSiteCount[c]);
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
                            ad.setHifiasmEcMatchStateFromReadPerspective(ReadId(readId), cand.is_match);

	                        /*
	                        Update assembly flags.

	                        These are monotonic (only ever set), so concurrent writes are safe.
	                        */
                        if (keep) ad.info.isInReadGraph = 1;
	                        const ReadId queryReadId = ReadId(readId);
	                        const uint32_t informativeCount = candInformativeSiteCount[c];
	                        if(ad.readIds[0] == queryReadId) {
	                            ad.informativeHetSiteCount0 = informativeCount;
	                        } else if(ad.readIds[1] == queryReadId) {
	                            ad.informativeHetSiteCount1 = informativeCount;
	                        }

	                        const bool debugThisOverlap =
	                            (ecDebugAlignmentId != invalid<uint32_t> && cand.alignmentId == ecDebugAlignmentId) ||
	                            (ecDebugByReadPair &&
	                             ad.readIds[0] == ecDebugReadPair[0] &&
	                             ad.readIds[1] == ecDebugReadPair[1]);
	                        if(debugThisOverlap) {
	                            std::lock_guard<std::mutex> lock(ecCoutMutex);
	                            const uint32_t informativeCount = candInformativeSiteCount[c];
	                            const int side =
	                                (ad.readIds[0] == queryReadId) ? 0 :
	                                (ad.readIds[1] == queryReadId) ? 1 : -1;
	                            cout << timestamp << "[EC-DBG] Stored for queryRead=" << queryReadId
	                                 << " alignmentId=" << cand.alignmentId
	                                 << " side=" << side
	                                 << " informativeCount=" << informativeCount
	                                 << " informative0=" << ad.informativeHetSiteCount0
	                                 << " informative1=" << ad.informativeHetSiteCount1
	                                 << endl;
	                        }
                    }
                    timing.finalizeFlags += seconds(steady_clock::now() - tFinalizeBegin);
                }
        });
    }


    for(auto& th : threads) th.join();

    // Finalize informative-site overlap scores now that both read perspectives have been populated.
    for(AlignmentData& ad : alignmentData) {
        ad.updateInformativeHetSiteScore();
    }

    // If EC debugging is enabled, print the final counters for the requested overlap now that
    // both read views have (typically) been processed.
    if(ecDebugAlignmentId != invalid<uint32_t>) {
        std::lock_guard<std::mutex> lock(ecCoutMutex);
        const AlignmentData& ad = alignmentData[ecDebugAlignmentId];
        cout << timestamp << "[EC-DBG] Final alignmentId=" << ecDebugAlignmentId
             << " reads=(" << ad.readIds[0] << "," << ad.readIds[1] << ")"
             << " informative0=" << ad.informativeHetSiteCount0
             << " informative1=" << ad.informativeHetSiteCount1
             << " score=" << ad.informativeHetSiteScore
             << " cisByBoth=" << (ad.isCisByBothSides() ? 1 : 0)
             << " keptByBoth=" << (ad.keptByBothSides() ? 1 : 0)
             << " deleteReasons0=0x" << std::hex << ad.deleteReasons0
             << " deleteReasons1=0x" << std::hex << ad.deleteReasons1
             << std::dec
             << endl;
    }
    if(ecDebugByReadPair) {
        std::lock_guard<std::mutex> lock(ecCoutMutex);
        for(uint32_t alignmentId = 0; alignmentId < alignmentData.size(); ++alignmentId) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if(ad.readIds[0] != ecDebugReadPair[0] || ad.readIds[1] != ecDebugReadPair[1]) continue;
            cout << timestamp << "[EC-DBG] Final alignmentId=" << alignmentId
                 << " reads=(" << ad.readIds[0] << "," << ad.readIds[1] << ")"
                 << " informative0=" << ad.informativeHetSiteCount0
                 << " informative1=" << ad.informativeHetSiteCount1
                 << " score=" << ad.informativeHetSiteScore
                 << " cisByBoth=" << (ad.isCisByBothSides() ? 1 : 0)
                 << " keptByBoth=" << (ad.keptByBothSides() ? 1 : 0)
                 << " deleteReasons0=0x" << std::hex << ad.deleteReasons0
                 << " deleteReasons1=0x" << std::hex << ad.deleteReasons1
                 << std::dec
                 << endl;
        }

        auto printReadSummary = [&](const EcDebugReadSummary& s) {
            if(!s.present) {
                cout << timestamp << "[EC-DBG] ReadSummary read=" << s.readId << " missing" << endl;
                return;
            }
            cout << timestamp << "[EC-DBG] ReadSummary read=" << s.readId
                 << " len=" << s.readLength
                 << " candidates=" << s.candidateCount
                 << " uniqueTargets=" << s.targetReads.size()
                 << endl;
        };
        printReadSummary(ecDebugReadSummaries[0]);
        printReadSummary(ecDebugReadSummaries[1]);

        if(ecDebugReadSummaries[0].present && ecDebugReadSummaries[1].present) {
            const auto& a = ecDebugReadSummaries[0].targetReads;
            const auto& b = ecDebugReadSummaries[1].targetReads;
            size_t i = 0, j = 0, inter = 0;
            while(i < a.size() && j < b.size()) {
                if(a[i] == b[j]) { ++inter; ++i; ++j; }
                else if(a[i] < b[j]) { ++i; }
                else { ++j; }
            }
            const size_t uni = a.size() + b.size() - inter;
            const double jaccard = uni ? double(inter) / double(uni) : 0.0;
            cout << timestamp << "[EC-DBG] NeighborIntersection sharedTargets=" << inter
                 << " unionTargets=" << uni
                 << " jaccard=" << jaccard
                 << endl;
        }

        const auto& o0 = ecDebugPairOverlapSummaries[0];
        const auto& o1 = ecDebugPairOverlapSummaries[1];
        if(o0.present && o1.present) {
            cout << timestamp << "[EC-DBG] PairOverlap read=" << o0.readId
                 << " qs=" << o0.qs << " qe=" << o0.qe
                 << " informativeSnpSites=" << o0.informativeSnpSitesInOverlap.size()
                 << endl;
            cout << timestamp << "[EC-DBG] PairOverlap read=" << o1.readId
                 << " qs=" << o1.qs << " qe=" << o1.qe
                 << " informativeSnpSites=" << o1.informativeSnpSitesInOverlap.size()
                 << endl;

            const uint32_t span0 = (o0.qe > o0.qs) ? (o0.qe - o0.qs) : 0;
            const uint32_t span1 = (o1.qe > o1.qs) ? (o1.qe - o1.qs) : 0;
            if(span0 && span1) {
                auto countApproxMappedMatches = [&](const EcDebugPairOverlapSummary& from,
                                                    const EcDebugPairOverlapSummary& to,
                                                    bool sameStrand,
                                                    uint32_t toSpan) -> size_t {
                    constexpr uint32_t window = 5;
                    size_t matches = 0;
                    if(from.informativeSnpSitesInOverlap.empty() || to.informativeSnpSitesInOverlap.empty()) return 0;
                    const auto& toSites = to.informativeSnpSitesInOverlap;
                    for(const uint32_t p : from.informativeSnpSitesInOverlap) {
                        if(p < from.qs || p >= from.qe) continue;
                        const double t = double(p - from.qs) / double(from.qe - from.qs);
                        const double t2 = sameStrand ? t : (1.0 - t);
                        const uint32_t mapped = to.qs + uint32_t(t2 * double(toSpan) + 0.5);
                        const uint32_t lo = (mapped > window) ? (mapped - window) : 0;
                        const uint32_t hi = mapped + window + 1;
                        auto it0 = std::lower_bound(toSites.begin(), toSites.end(), lo);
                        if(it0 != toSites.end() && *it0 < hi) {
                            ++matches;
                        }
                    }
                    return matches;
                };

                const size_t m01 = countApproxMappedMatches(o0, o1, ecDebugPairIsSameStrand, span1);
                const size_t m10 = countApproxMappedMatches(o1, o0, ecDebugPairIsSameStrand, span0);
                cout << timestamp << "[EC-DBG] PairOverlapApproxShared sameStrand=" << (ecDebugPairIsSameStrand ? 1 : 0)
                     << " mappedMatches0to1=" << m01
                     << " mappedMatches1to0=" << m10
                     << endl;
            }
        } else {
            cout << timestamp << "[EC-DBG] PairOverlap missing: present0=" << (o0.present ? 1 : 0)
                 << " present1=" << (o1.present ? 1 : 0) << endl;
        }
    }

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
    uint64_t transMarked = 0;
    uint64_t keptByBoth = 0;
    for (const auto& ad : alignmentData) {
        if ((ad.deleteReasons0 & AlignmentData::DeleteReasonPhase) ||
            (ad.deleteReasons1 & AlignmentData::DeleteReasonPhase)) {
            ++transMarked;
        }
        if (ad.keptByBothSides()) {
            ++keptByBoth;
        }
    }

    cout << timestamp << "Parity EC Round 1 wall time: " << tAll << " s"
         << " (reads=" << total.readsVisited
         << ", withAlignments=" << total.readsWithAlignments
         << ", withCandidates=" << total.readsWithCandidates
         << ", transMarked=" << transMarked
         << ", keptByBoth=" << keptByBoth << ")" << endl;

    cout << timestamp << "Parity EC Round 1 Complete." << endl;
}

void Assembler::performHifiasmECFinalFilteringParity(uint64_t /* threadCount */)
{
    cout << timestamp << "=== Hifiasm Parity EC Final Filtering (ha_ec_ff) ===" << endl;
}



// Experimental EC parity using induced alignments through marker graph vertices.
//
// Instead of the SNP/SV detection + DP chaining pipeline, this function uses the
// *ordering consistency* of shared marker graph vertices between each pair of
// overlapping reads as the cis/trans phasing signal.
//
// For each read R (forward strand = strand 0):
//   For each overlapping candidate read R2:
//     Induced alignment = sequence of (R_ordinal_i, R2_ordinal_j) pairs at each
//     marker graph vertex shared between R and R2.
//   A candidate is "consistent" (cis) if the R2 ordinals are monotonically
//   non-decreasing when sorted by R ordinals — i.e., the shared vertices appear
//   in the same relative order in both reads.
//   A candidate is "inconsistent" (trans) if there is at least one ordinal inversion.
//
// A read is "informative" (analogous to having het SNP sites) if ANY of its
// candidates show an ordering inconsistency.  When informative, we keep only
// consistent candidates that also meet the minimum shared-vertex threshold.
// When not informative, all candidates are kept (no phasing signal available).
//
// Prerequisite: marker graph vertices must have been created before calling this.
// The function is a drop-in replacement for performHifiasmECParity in main.cpp.
void Assembler::performHifiasmECParityWithMarkerGraph(uint64_t threadCount)
{
    cout << timestamp << "=== Marker Graph Projected EC Parity ===" << endl;

    DINARA_ASSERT(markers && markers->isOpen());
    DINARA_ASSERT(markerGraph.vertices().isOpen());

    const uint64_t readCount = reads->readCount();
    const auto tBeginAll = steady_clock::now();

    // Minimum shared vertices required to label a candidate as cis.
    // Candidates with fewer shared vertices are kept regardless (insufficient evidence).
    const uint32_t minSharedVertices = 3;

    // Per-candidate bookkeeping.
    struct MgCandidate {
        uint32_t alignmentId  = 0;
        uint32_t targetReadId = 0;
        bool     isRev        = false;   // true = different-strand (RC) overlap
        // Induced alignment: (R_ordinal_i, target_ordinal_j) at shared vertices.
        // Populated during Stage 2; sorted by R_ordinal in Stage 3.
        vector<pair<uint32_t, uint32_t>> induced;
        bool     consistent   = true;    // no ordinal inversion found
        uint32_t sharedCount  = 0;       // |induced| after dedup
    };

    atomic<uint64_t> totalReads{0};
    atomic<uint64_t> readsWithCandidates{0};
    atomic<uint64_t> readsInformative{0};
    atomic<uint64_t> transMarkedNew{0};

    vector<thread> threads;
    const uint64_t chunkSize = max(uint64_t(1), readCount / threadCount);

    for (uint64_t t = 0; t < threadCount; ++t) {
        threads.emplace_back([&, t]() {
            const uint64_t tStart = t * chunkSize;
            const uint64_t tEnd   = (t == threadCount - 1) ? readCount : (t + 1) * chunkSize;

            // Thread-local scratch to avoid per-read heap allocations.
            vector<MgCandidate> candidates;
            // Sorted (key, candidateIdx) for O(log n) lookup during vertex traversal.
            // key = (readId << 1) | strand, where strand encodes isRev.
            vector<pair<uint64_t, uint32_t>> candLookup;

            uint64_t localReads = 0, localWithCandidates = 0;
            uint64_t localInformative = 0, localTransMarked = 0;

            for (uint64_t readId = tStart; readId < tEnd; ++readId) {
                ++localReads;
                const OrientedReadId orientedReadId(ReadId(readId), 0);
                if (orientedReadId.getValue() >= alignmentTable.size()) continue;

                const auto& alignments = alignmentTable[orientedReadId.getValue()];
                if (alignments.empty()) continue;

                // ----------------------------------------------------------------
                // Stage 1: gather non-deleted candidates.
                // Normalize so R is always in the strand-0 (forward) frame.
                // Same normalization logic as performHifiasmECParity.
                // ----------------------------------------------------------------
                candidates.clear();
                candLookup.clear();

                for (const uint32_t alignmentId : alignments) {
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if (!ad.keptByBothSides()) continue;

                    OrientedReadId o0(ad.readIds[0], 0);
                    OrientedReadId o1(ad.readIds[1], ad.isSameStrand ? 0 : 1);

                    // Swap so that o0 is always R.
                    if (o0.getReadId() != ReadId(readId)) {
                        swap(o0, o1);
                    }
                    DINARA_ASSERT(o0.getReadId() == ReadId(readId));

                    // Flip if R ended up in strand-1 after the swap.
                    if (o0.getStrand() != 0) {
                        o0.flipStrand();
                        o1.flipStrand();
                    }
                    DINARA_ASSERT(o0.getStrand() == 0);

                    MgCandidate cand;
                    cand.alignmentId  = alignmentId;
                    cand.targetReadId = uint32_t(o1.getReadId());
                    cand.isRev        = (o1.getStrand() != 0);

                    // Encode (readId, strand) into a uint64 key.
                    // isRev=false -> strand=0 ; isRev=true -> strand=1.
                    const uint64_t key =
                        (uint64_t(cand.targetReadId) << 1) | (cand.isRev ? 1U : 0U);
                    candLookup.emplace_back(key, uint32_t(candidates.size()));
                    candidates.push_back(move(cand));
                }

                if (candidates.empty()) continue;
                ++localWithCandidates;

                sort(candLookup.begin(), candLookup.end());

                // ----------------------------------------------------------------
                // Stage 2: build induced alignments via shared marker graph vertices.
                //
                // For each marker ordinal i in R (strand=0):
                //   vi = vertexTable[markerId]
                //   for each other marker at vi → (otherOriented, otherOrdinal):
                //     strand convention:
                //       isRev=false candidate → we want otherOriented.strand == 0
                //       isRev=true  candidate → we want otherOriented.strand == 1
                //     if found in candLookup: record (i, otherOrdinal) for that candidate.
                //
                // Ordinal monotonicity holds in both same-strand and RC overlaps:
                //   as R's ordinal i increases, the target ordinal j also increases
                //   (both iterate left-to-right in their respective oriented frames).
                // ----------------------------------------------------------------
                const MarkerId firstMkR =
                    MarkerId(markers->begin(orientedReadId.getValue()) - markers->begin());
                const uint32_t mkCountR =
                    uint32_t(markers->size(orientedReadId.getValue()));

                for (uint32_t i = 0; i < mkCountR; ++i) {
                    const MarkerId mkId = firstMkR + i;
                    const auto vi = markerGraph.vertexTable[mkId];
                    if (vi == MarkerGraph::invalidCompressedVertexId) continue;

                    for (const MarkerId otherMkId : markerGraph.getVertexMarkerIds(vi)) {
                        if (otherMkId == mkId) continue;   // skip R's own marker

                        OrientedReadId otherOriented;
                        uint32_t otherOrdinal;
                        tie(otherOriented, otherOrdinal) = findMarkerId(otherMkId);

                        // Lookup by (readId << 1 | strand).
                        const uint64_t key =
                            (uint64_t(otherOriented.getReadId()) << 1) |
                            uint64_t(otherOriented.getStrand());
                        const auto lo = lower_bound(
                            candLookup.begin(), candLookup.end(),
                            make_pair(key, uint32_t(0)));
                        if (lo == candLookup.end() || lo->first != key) continue;

                        candidates[lo->second].induced.emplace_back(i, otherOrdinal);
                        ++candidates[lo->second].sharedCount;
                    }
                }

                // ----------------------------------------------------------------
                // Stage 3: check ordering consistency for each candidate.
                //   Sort induced alignment by R ordinal.
                //   Deduplicate ties at the same R ordinal (keep first).
                //   Check if target ordinals are non-decreasing → consistent.
                // ----------------------------------------------------------------
                bool anyInconsistency = false;

                for (auto& cand : candidates) {
                    if (cand.sharedCount < 2) continue;   // need ≥2 points to check order

                    sort(cand.induced.begin(), cand.induced.end());   // sort by R ordinal

                    // Deduplicate: at equal R ordinals keep only the first entry.
                    {
                        uint32_t w = 0;
                        for (uint32_t r = 0; r < uint32_t(cand.induced.size()); ++r) {
                            if (w == 0 || cand.induced[r].first != cand.induced[w - 1].first) {
                                cand.induced[w++] = cand.induced[r];
                            }
                        }
                        cand.induced.resize(w);
                        cand.sharedCount = w;
                    }

                    if (cand.sharedCount < 2) continue;

                    // Monotone check: target ordinals non-decreasing as R ordinal increases.
                    bool mono = true;
                    for (uint32_t k = 1; k < uint32_t(cand.induced.size()); ++k) {
                        if (cand.induced[k].second < cand.induced[k - 1].second) {
                            mono = false;
                            break;
                        }
                    }
                    cand.consistent = mono;
                    if (!mono) anyInconsistency = true;
                }

                // "Informative" = at least one candidate has an ordering violation,
                // analogous to having validated het SNP/SV sites.
                const bool isInformative = anyInconsistency;
                if (isInformative) ++localInformative;

                // ----------------------------------------------------------------
                // Stage 4: phase decision and flag write-back.
                //   If informative : keep only consistent candidates with sufficient
                //                    shared-vertex evidence.
                //   If not informative : keep all (no signal to filter on).
                //
                // Same flag semantics as performHifiasmECParity:
                //   keep  → clear DeleteReasonPhase for this read's side
                //   trans → set   DeleteReasonPhase for this read's side
                // ----------------------------------------------------------------
                for (const auto& cand : candidates) {
                    AlignmentData& ad = alignmentData[cand.alignmentId];

                    bool keep;
                    if (!isInformative) {
                        keep = true;
                    } else {
                        keep = cand.consistent && (cand.sharedCount >= minSharedVertices);
                    }

                    if (keep) {
                        ad.clearDeleteReasonsFromReadPerspective(
                            ReadId(readId), AlignmentData::DeleteReasonPhase);
                        ad.info.isInReadGraph = 1;
                    } else {
                        ad.addDeleteReasonsFromReadPerspective(
                            ReadId(readId), AlignmentData::DeleteReasonPhase);
                        ++localTransMarked;
                    }
                }
            }

            totalReads          += localReads;
            readsWithCandidates += localWithCandidates;
            readsInformative    += localInformative;
            transMarkedNew      += localTransMarked;
        });
    }

    for (auto& th : threads) th.join();

    // Recompute the combined informative-site overlap score (now vertex-based).
    for (AlignmentData& ad : alignmentData) {
        ad.updateInformativeHetSiteScore();
    }

    const double tAll = seconds(steady_clock::now() - tBeginAll);
    uint64_t keptByBoth = 0;
    for (const AlignmentData& ad : alignmentData) {
        if (ad.keptByBothSides()) ++keptByBoth;
    }

    cout << timestamp
         << "Marker Graph EC wall time: " << tAll << " s"
         << "  reads=" << totalReads
         << "  withCandidates=" << readsWithCandidates
         << "  informative=" << readsInformative
         << "  transMarked=" << transMarkedNew
         << "  keptByBoth=" << keptByBoth
         << endl;
    cout << timestamp << "Marker Graph EC Complete." << endl;
}



// Debug helper: run het-site detection for a single read and print all SNP/SV sites
// together with the per-overlap evidence, so they can be inspected interactively.
//
// Call this after computeBaseAlignmentsAndStore() (which populates alignmentTable,
// alignmentData and alignedEvidenceStore).
void Assembler::debugPrintHetSitesForRead(uint64_t readId)
{
    using std::setw;
    using std::left;
    using std::right;

    const OrientedReadId orientedReadId(ReadId(readId), 0);
    if (orientedReadId.getValue() >= alignmentTable.size()) {
        cout << "[HetDebug] readId=" << readId << " not in alignmentTable (size="
             << alignmentTable.size() << ")\n";
        return;
    }
    const auto& alignments = alignmentTable[orientedReadId.getValue()];
    const uint32_t readLen = uint32_t(reads->getRead(ReadId(readId)).baseCount);

    cout << "\n";
    cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    cout << "║  Het-site debug  read=" << readId << "-0   len=" << readLen << " bp\n";
    cout << "╚══════════════════════════════════════════════════════════════════╝\n";
    cout << "  Alignments in table : " << alignments.size() << "\n";

    // -----------------------------------------------------------------------
    // Gather candidates (exact copy of the logic in performHifiasmECParity).
    // -----------------------------------------------------------------------
    HifiasmECScratchPad scratch;
    scratch.clear();
    auto& candidates = scratch.candidates;
    candidates.reserve(alignments.size());

    for (uint32_t alignmentId : alignments) {
        const AlignmentData& ad = alignmentData[alignmentId];
        // Include reads deleted only for phase reasons (trans) so that
        // detectHetSites can still see both allele groups.  Only exclude reads
        // deleted for other reasons (containment, chemistry, etc.).
        constexpr AlignmentData::DeleteReasonMask nonPhase =
            ~AlignmentData::DeleteReasonPhase;
        if ((ad.deleteReasons0 & nonPhase) || (ad.deleteReasons1 & nonPhase)) continue;

        CandidateEC candidate;
        candidate.alignmentId = alignmentId;

        const OrientedReadId queryOriented(ReadId(readId), 0);
        OrientedReadId o0(ad.readIds[0], 0);
        OrientedReadId o1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
        AlignmentInfo orientedInfo = ad.info;

        if (o0.getReadId() != queryOriented.getReadId()) {
            swap(o0, o1);
            orientedInfo.swap();
        }
        if (o0.getStrand() != queryOriented.getStrand()) {
            o0.flipStrand();
            o1.flipStrand();
            orientedInfo.reverseComplement();
        }

        uint32_t qsCore = 0, qeCore = 0, tsCoreOriented = 0, teCoreOriented = 0;
        bool coordinatesFromMarkers = false;
        const uint32_t kmerLen = assemblerInfo.isOpen ? uint32_t(assemblerInfo->k) : 0U;
        if (markers && assemblerInfo.isOpen) {
            const auto m0 = (*markers)[o0.getValue()];
            const auto m1 = (*markers)[o1.getValue()];
            if (!m0.empty() && !m1.empty()) {
                qsCore = m0[orientedInfo.data[0].firstOrdinal].position;
                qeCore = m0[orientedInfo.data[0].lastOrdinal].position + kmerLen;
                tsCoreOriented = m1[orientedInfo.data[1].firstOrdinal].position;
                teCoreOriented = m1[orientedInfo.data[1].lastOrdinal].position + kmerLen;
                coordinatesFromMarkers = true;
            }
        }
        if (!coordinatesFromMarkers) {
            if (ad.readIds[0] == ReadId(readId)) {
                qsCore = ad.qs; qeCore = ad.qe;
                tsCoreOriented = ad.ts; teCoreOriented = ad.te;
            } else {
                qsCore = ad.ts; qeCore = ad.te;
                tsCoreOriented = ad.qs; teCoreOriented = ad.qe;
            }
        }
        const uint32_t tLen = uint32_t(reads->getRead(o1.getReadId()).baseCount);
        uint32_t tsFwd = tsCoreOriented, teFwd = teCoreOriented;
        if (coordinatesFromMarkers && o1.getStrand() != 0) {
            tsFwd = tLen - teCoreOriented;
            teFwd = tLen - tsCoreOriented;
        }

        candidate.qs       = qsCore;
        candidate.qe       = qeCore;
        candidate.ts       = tsFwd;
        candidate.te       = teFwd;
        candidate.targetId = uint32_t(o1.getReadId());
        candidate.isRev    = (o1.getStrand() != 0);
        candidate.is_match = 1;
        candidate.strong   = 0;
        candidates.push_back(candidate);
    }

    cout << "  Kept candidates     : " << candidates.size() << "\n\n";
    if (candidates.empty()) {
        cout << "  (no candidates after keptByBothSides filter)\n";
        return;
    }

    // -----------------------------------------------------------------------
    // Run SNP detection.
    // -----------------------------------------------------------------------
    detectHetSites(*this, *reads, readId, alignmentData, scratch);

    // -----------------------------------------------------------------------
    // Print candidate summary.
    // -----------------------------------------------------------------------
    cout << "  Candidates (qs..qe in query-forward coords):\n";
    cout << "  " << left << setw(6) << "#"
         << setw(8) << "tgtRead"
         << setw(5) << "str"
         << setw(8) << "qs"
         << setw(8) << "qe"
         << setw(8) << "span"
         << setw(8) << "ts"
         << setw(8) << "te"
         << "alignId\n";
    cout << "  " << std::string(59, '-') << "\n";
    for (size_t ci = 0; ci < candidates.size(); ++ci) {
        const auto& c = candidates[ci];
        cout << "  " << left << setw(6) << ci
             << setw(8) << c.targetId
             << setw(5) << (c.isRev ? "RC" : "FW")
             << setw(8) << c.qs
             << setw(8) << c.qe
             << setw(8) << (c.qe - c.qs)
             << setw(8) << c.ts
             << setw(8) << c.te
             << c.alignmentId << "\n";
    }

    // -----------------------------------------------------------------------
    // Print SNP sites table.
    // -----------------------------------------------------------------------
    const auto& snpStats = scratch.snpStats;
    // Count how many sites pass the alt>=3 filter.
    size_t snpAlt3Count = 0;
    for (const auto& s : snpStats) { if (s.occ_1 >= 3) ++snpAlt3Count; }
    cout << "\n  SNP sites with alt>=3 (" << snpAlt3Count
         << " of " << snpStats.size() << " total):\n";

    // Build row → [overlapIdx list] from hapEvidence (needed for printing and
    // for the alt-supporter set used to filter the cis read list below).
    vector<vector<uint32_t>> rowToOverlaps(snpStats.size());
    for (const auto& ev : scratch.hapEvidence) {
        if (ev.type == 1 && ev.overlapSite != invalid<uint32_t>
                && ev.overlapSite < snpStats.size()) {
            rowToOverlaps[ev.overlapSite].push_back(ev.overlapID);
        }
    }

    // Collect reads that support the alt allele at any SNP with occ_1 >= 3.
    // These are expected to be on the opposite haplotype from readId.
    std::unordered_set<uint32_t> altSupporters;
    for (size_t ri = 0; ri < snpStats.size(); ++ri) {
        if (snpStats[ri].occ_1 < 3) continue;
        for (const uint32_t oi : rowToOverlaps[ri]) {
            if (oi < candidates.size())
                altSupporters.insert(candidates[oi].targetId);
        }
    }

    if (!snpStats.empty()) {
        cout << "  " << left
             << setw(5)  << "row"
             << setw(9)  << "pos"
             << setw(5)  << "ref"
             << setw(5)  << "alt"
             << setw(7)  << "occ0"
             << setw(7)  << "occ1"
             << setw(6)  << "HP"
             << setw(8)  << "fwdRef"
             << "score\n";
        cout << "  " << std::string(52, '-') << "\n";
        for (size_t ri = 0; ri < snpStats.size(); ++ri) {
            const auto& s = snpStats[ri];
            if (s.occ_1 < 3) continue;
            cout << "  " << left
                 << setw(5)  << ri
                 << setw(9)  << s.site
                 << setw(5)  << s.refBase
                 << setw(5)  << s.altBase
                 << setw(7)  << s.occ_0
                 << setw(7)  << s.occ_1
                 << setw(6)  << (s.is_homopolymer ? "yes" : "no")
                 << setw(8)  << s.fwd_ref_cov
                 << s.score << "\n";
        }

        // Per-overlap support breakdown for each SNP row.
        cout << "\n  Per-overlap alt evidence (row → overlap, alt>=3 only):\n";
        for (size_t ri = 0; ri < snpStats.size(); ++ri) {
            const auto& s = snpStats[ri];
            if (s.occ_1 < 3) continue;
            const auto& ovList = rowToOverlaps[ri];
            cout << "    row " << ri << " pos=" << s.site
                 << " " << s.refBase << "->" << s.altBase
                 << " | alt supported by " << ovList.size() << " overlap(s): ";
            for (uint32_t oi : ovList) {
                if (oi < candidates.size()) {
                    cout << "read" << candidates[oi].targetId
                         << "(" << (candidates[oi].isRev ? "RC" : "FW") << ") ";
                }
            }
            cout << "\n";
        }
    }

    // -----------------------------------------------------------------------
    // Run and print SV (large indel) sites.
    // -----------------------------------------------------------------------
    detectSVSites(*this, *reads, readId, alignmentData, scratch);

    const auto& svStats = scratch.svStats;
    cout << "\n  SV sites (" << svStats.size() << " rows after emission filters):\n";
    if (!svStats.empty()) {
        cout << "  " << left
             << setw(5)  << "row"
             << setw(9)  << "pos"
             << setw(7)  << "occ0"
             << setw(7)  << "occ1" << "\n";
        cout << "  " << std::string(28, '-') << "\n";
        for (size_t ri = 0; ri < svStats.size(); ++ri) {
            const auto& s = svStats[ri];
            cout << "  " << left
                 << setw(5)  << ri
                 << setw(9)  << s.site
                 << setw(7)  << s.occ_0
                 << setw(7)  << s.occ_1 << "\n";
        }
        // Per-overlap SV support.
        cout << "\n  Per-overlap SV evidence (row → overlap):\n";
        vector<vector<uint32_t>> svRowToOverlaps(svStats.size());
        for (const auto& ev : scratch.svEvidence) {
            if (ev.overlapSite != invalid<uint32_t>
                    && ev.overlapSite < svStats.size()) {
                svRowToOverlaps[ev.overlapSite].push_back(ev.overlapID);
            }
        }
        for (size_t ri = 0; ri < svStats.size(); ++ri) {
            const auto& s = svStats[ri];
            const auto& ovList = svRowToOverlaps[ri];
            cout << "    row " << ri << " pos=" << s.site
                 << " | supported by " << ovList.size() << " overlap(s): ";
            for (uint32_t oi : ovList) {
                if (oi < candidates.size()) {
                    cout << "read" << candidates[oi].targetId
                         << "(" << (candidates[oi].isRev ? "RC" : "FW") << ") ";
                }
            }
            cout << "\n";
        }
    }

    // -----------------------------------------------------------------------
    // Print cis reads: neighbours of readId whose side has no DeleteReasonPhase
    // AND are not alt-supporters at any confident SNP (occ_1 >= 3).
    // -----------------------------------------------------------------------
    cout << "\n  Cis reads (no DeleteReasonPhase, not an alt-supporter at occ_1>=3 SNP):\n";
    {
        vector<uint32_t> cisReads;
        for (const uint32_t alnId : alignments) {
            const AlignmentData& ad = alignmentData[alnId];
            const bool isR0 = (ad.readIds[0] == ReadId(readId));
            const AlignmentData::DeleteReasonMask ourReasons =
                isR0 ? ad.deleteReasons0 : ad.deleteReasons1;
            if (ourReasons & AlignmentData::DeleteReasonPhase) continue;
            const uint32_t partner = isR0 ? uint32_t(ad.readIds[1])
                                          : uint32_t(ad.readIds[0]);
            if (altSupporters.count(partner)) continue;
            cisReads.push_back(partner);
        }
        std::sort(cisReads.begin(), cisReads.end());
        cisReads.erase(std::unique(cisReads.begin(), cisReads.end()), cisReads.end());
        cout << "  count: " << cisReads.size() << "\n";
        for (const uint32_t r : cisReads) {
            cout << "    read" << r << "\n";
        }

        // Subset: cis reads that also appear as a candidate (valid coordinate entry).
        std::unordered_set<uint32_t> candidateTargets;
        for (const CandidateEC& c : candidates) {
            candidateTargets.insert(c.targetId);
        }
        cout << "\n  Cis direct-overlap reads (also in candidates):\n";
        uint32_t directCisCount = 0;
        for (const uint32_t r : cisReads) {
            if (candidateTargets.count(r)) {
                cout << "    read" << r << "\n";
                ++directCisCount;
            }
        }
        cout << "  count: " << directCisCount << "\n";
    }

    cout << "\n[HetDebug] done for read " << readId << "-0.\n\n";
}



void Assembler::debugDumpAlignedEvidenceForRead(uint64_t readId)
{
    const OrientedReadId orientedReadId(ReadId(readId), 0);
    if (orientedReadId.getValue() >= alignmentTable.size()) {
        cout << "[EvidenceDump] readId=" << readId << " not in alignmentTable (size="
             << alignmentTable.size() << ")\n";
        return;
    }
    if (alignedEvidenceStore.index.empty()) {
        cout << "[EvidenceDump] alignedEvidenceStore is empty. "
             << "Run computeBaseAlignmentsAndStore() first.\n";
        return;
    }

    const auto& alignments = alignmentTable[orientedReadId.getValue()];
    const uint32_t readLen = uint32_t(reads->getRead(ReadId(readId)).baseCount);
    constexpr uint32_t svMinLen = 20;

    cout << "\n";
    cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    cout << "║  Evidence dump  read=" << readId << "-0   len=" << readLen << " bp\n";
    cout << "╚══════════════════════════════════════════════════════════════════╝\n";
    cout << "  Alignments in table : " << alignments.size() << "\n";

    uint64_t totalSnps = 0;
    uint64_t totalIndels = 0;
    uint64_t totalSvCandidates = 0;

    for (const uint32_t alignmentId : alignments) {
        if (alignmentId >= alignmentData.size()) {
            continue;
        }
        const AlignmentData& ad = alignmentData[alignmentId];
        const size_t evidenceId = ad.info.alignmentId;
        if (evidenceId == invalid<size_t> || evidenceId >= alignedEvidenceStore.index.size()) {
            continue;
        }

        const bool queryIsRead0 = (ad.readIds[0] == ReadId(readId));
        const uint32_t partnerReadId = queryIsRead0 ? uint32_t(ad.readIds[1]) : uint32_t(ad.readIds[0]);
        const bool partnerIsRev = queryIsRead0 ? (!ad.isSameStrand) : false;
        const uint32_t qs = queryIsRead0 ? ad.qs : ad.ts;
        const uint32_t qe = queryIsRead0 ? ad.qe : ad.te;
        const uint32_t ts = queryIsRead0 ? ad.ts : ad.qs;
        const uint32_t te = queryIsRead0 ? ad.te : ad.qe;

        vector<pair<uint32_t, uint8_t> > snps;
        snps.reserve(32);
        auto collectSnp = [&](uint32_t pos, uint8_t base) {
            snps.emplace_back(pos, base);
        };
        const uint32_t evidenceId32 = uint32_t(evidenceId);
        if (queryIsRead0) {
            alignedEvidenceStore.forEachSnp1InRange(evidenceId32, qs, qe, collectSnp);
        } else {
            alignedEvidenceStore.forEachSnp0InRange(evidenceId32, qs, qe, collectSnp);
        }

        span<const IndelEvidence> indels =
            queryIsRead0 ?
            alignedEvidenceStore.getIndels1(evidenceId32) :
            alignedEvidenceStore.getIndels0(evidenceId32);

        vector<IndelEvidence> localIndels;
        localIndels.reserve(indels.size());
        for (const IndelEvidence& indel : indels) {
            if (indel.pos() >= qs && indel.pos() < qe) {
                localIndels.push_back(indel);
            }
        }

        uint64_t localSvCount = 0;
        for (const IndelEvidence& indel : localIndels) {
            if (indel.len() >= svMinLen) {
                ++localSvCount;
            }
        }

        totalSnps += snps.size();
        totalIndels += localIndels.size();
        totalSvCandidates += localSvCount;

        cout << "\n";
        cout << "[EvidenceDump] alignmentId=" << alignmentId
             << " evidenceId=" << evidenceId
             << " partner=read" << partnerReadId
             << "-" << (partnerIsRev ? 1 : 0)
             << " qs=" << qs
             << " qe=" << qe
             << " ts=" << ts
             << " te=" << te
             << " kept=" << (ad.keptByBothSides() ? "yes" : "no")
             << "\n";

        cout << "  SNPs (" << snps.size() << "):";
        if (snps.empty()) {
            cout << " none\n";
        } else {
            cout << "\n";
            for (const auto& [pos, base] : snps) {
                cout << "    pos=" << pos
                     << " ref=" << reads->getOrientedReadBase(orientedReadId, pos).character()
                     << " alt=" << Base::fromInteger(base).character()
                     << "\n";
            }
        }

        cout << "  Indels (" << localIndels.size() << "):";
        if (localIndels.empty()) {
            cout << " none\n";
        } else {
            cout << "\n";
            for (const IndelEvidence& indel : localIndels) {
                cout << "    pos=" << indel.pos()
                     << " type=" << (indel.isInsertion() ? "INS" : "DEL")
                     << " len=" << indel.len();
                if (indel.len() >= svMinLen) {
                    cout << "  [SV]";
                }
                cout << "\n";
            }
        }
    }

    cout << "\n";
    cout << "[EvidenceDump] summary for read " << readId << "-0:\n";
    cout << "  total SNP observations   : " << totalSnps << "\n";
    cout << "  total indel observations : " << totalIndels << "\n";
    cout << "  total SV candidates      : " << totalSvCandidates
         << " (indel len >= " << svMinLen << ")\n\n";
}

// ---------------------------------------------------------------------------
// debugDumpSnpSitesForRead
//
// Aggregates SNP evidence from all alignments involving the focal read
// (strand 0) and prints one line per (position, alt-base) pair where both
// the ref allele and the alt allele have at least minSupport supporting reads.
//
// ref_support at pos = alignments covering pos that show no mismatch there
// alt_support at pos = alignments covering pos that show a specific alt base
// ---------------------------------------------------------------------------
void Assembler::debugDumpSnpSitesForRead(uint64_t readId, uint32_t minSupport)
{
    const OrientedReadId orientedReadId(ReadId(readId), 0);
    if (orientedReadId.getValue() >= alignmentTable.size()) {
        cout << "[SnpSites] readId=" << readId << " not in alignmentTable\n";
        return;
    }
    if (alignedEvidenceStore.index.empty()) {
        cout << "[SnpSites] alignedEvidenceStore is empty. "
             << "Run computeBaseAlignmentsAndStore() first.\n";
        return;
    }

    const auto& alignments = alignmentTable[orientedReadId.getValue()];

    constexpr uint32_t svMinLen = 20;

    // Per alignment: span [qs, qe) on focal read + SNPs + SVs observed there.
    struct AlignSpan {
        uint32_t qs, qe;
        vector<pair<uint32_t, uint8_t>> snps; // {pos on focal read, alt base index}
        vector<IndelEvidence>           svs;  // indels with len >= svMinLen
    };
    vector<AlignSpan> spans;
    spans.reserve(alignments.size());

    for (const uint32_t alignmentId : alignments) {
        if (alignmentId >= alignmentData.size()) continue;
        const AlignmentData& ad = alignmentData[alignmentId];
        const size_t evidenceId = ad.info.alignmentId;
        if (evidenceId == invalid<size_t> || evidenceId >= alignedEvidenceStore.index.size()) continue;

        const bool queryIsRead0 = (ad.readIds[0] == ReadId(readId));
        const uint32_t qs = queryIsRead0 ? ad.qs : ad.ts;
        const uint32_t qe = queryIsRead0 ? ad.qe : ad.te;
        const uint32_t evidenceId32 = uint32_t(evidenceId);

        AlignSpan span;
        span.qs = qs;
        span.qe = qe;

        auto collectSnp = [&](uint32_t pos, uint8_t base) {
            span.snps.emplace_back(pos, base);
        };
        if (queryIsRead0)
            alignedEvidenceStore.forEachSnp1InRange(evidenceId32, qs, qe, collectSnp);
        else
            alignedEvidenceStore.forEachSnp0InRange(evidenceId32, qs, qe, collectSnp);

        const dinara::span<const IndelEvidence> indels =
            queryIsRead0 ?
            alignedEvidenceStore.getIndels1(evidenceId32) :
            alignedEvidenceStore.getIndels0(evidenceId32);
        for (const IndelEvidence& indel : indels)
            if (indel.pos() >= qs && indel.pos() < qe && indel.len() >= svMinLen)
                span.svs.push_back(indel);

        spans.push_back(std::move(span));
    }

    // Aggregate: for each (pos, altBase) count alt-supporting alignments.
    // Also build a set of SNP positions so we can compute ref support.
    map<uint32_t, map<uint8_t, uint32_t>> altCount; // pos → baseIdx → count
    for (const auto& sp : spans)
        for (const auto& [pos, base] : sp.snps)
            altCount[pos][base]++;

    // Ref support at pos = alignments covering pos with no SNP at pos.
    // Build a fast lookup: per alignment, which positions have a SNP.
    // Then for each candidate position, count spans that cover it without a SNP.
    auto refSupport = [&](uint32_t pos) -> uint32_t {
        uint32_t n = 0;
        for (const auto& sp : spans) {
            if (pos < sp.qs || pos >= sp.qe) continue;
            bool hasMismatch = false;
            for (const auto& [spos, sbase] : sp.snps)
                if (spos == pos) { hasMismatch = true; break; }
            if (!hasMismatch) ++n;
        }
        return n;
    };

    // Print header.
    cout << "\n[SnpSites] read=" << readId << "-0"
         << "  alignments=" << spans.size()
         << "  minSupport=" << minSupport << "\n";
    cout << "  " << string(60, '-') << "\n";
    cout << "  pos       ref  alt  refSupport  altSupport\n";
    cout << "  " << string(60, '-') << "\n";

    using std::setw;
    uint32_t printed = 0;
    for (const auto& [pos, baseCounts] : altCount) {
        const char refBase = reads->getOrientedReadBase(orientedReadId, pos).character();
        const uint32_t refSup = refSupport(pos);
        for (const auto& [baseIdx, altSup] : baseCounts) {
            if (altSup < minSupport || refSup < minSupport) continue;
            const char altBase = Base::fromInteger(baseIdx).character();
            cout << "  " << setw(8) << pos
                 << "  " << refBase
                 << "    " << altBase
                 << "    " << setw(10) << refSup
                 << "  " << setw(10) << altSup
                 << "\n";
            ++printed;
        }
    }
    if (printed == 0)
        cout << "  (no SNP sites with both alleles >= " << minSupport << " support)\n";
    cout << "  " << string(60, '-') << "\n";
    cout << "  Total SNP sites printed: " << printed << "\n";

    // ---- SV aggregation ----
    // Key: (pos, type, len); value: support count (alt-supporting alignments).
    // ref_support at pos = alignments covering pos with no SV of the same type+len there.
    struct SvKey {
        uint32_t pos;
        uint8_t  type; // 0=INS, 1=DEL
        uint32_t len;
        bool operator<(const SvKey& o) const {
            if (pos  != o.pos)  return pos  < o.pos;
            if (type != o.type) return type < o.type;
            return len < o.len;
        }
    };
    map<SvKey, uint32_t> svAltCount;
    for (const auto& sp : spans)
        for (const IndelEvidence& sv : sp.svs)
            svAltCount[{sv.pos(), sv.type(), sv.len()}]++;

    auto svRefSupport = [&](uint32_t pos, uint8_t type, uint32_t len) -> uint32_t {
        uint32_t n = 0;
        for (const auto& sp : spans) {
            if (pos < sp.qs || pos >= sp.qe) continue;
            bool hasSv = false;
            for (const IndelEvidence& sv : sp.svs)
                if (sv.pos() == pos && sv.type() == type && sv.len() == len)
                    { hasSv = true; break; }
            if (!hasSv) ++n;
        }
        return n;
    };

    cout << "\n";
    cout << "  " << string(60, '-') << "\n";
    cout << "  SV sites (indel len >= " << svMinLen << ")\n";
    cout << "  " << string(60, '-') << "\n";
    cout << "  pos       type  len   refSupport  altSupport\n";
    cout << "  " << string(60, '-') << "\n";

    uint32_t svPrinted = 0;
    for (const auto& [key, altSup] : svAltCount) {
        const uint32_t refSup = svRefSupport(key.pos, key.type, key.len);
        cout << "  " << setw(8) << key.pos
             << "  " << (key.type == 0 ? "INS" : "DEL")
             << "  " << setw(5) << key.len
             << "  " << setw(10) << refSup
             << "  " << setw(10) << altSup
             << "\n";
        ++svPrinted;
    }
    if (svPrinted == 0)
        cout << "  (no SV sites found)\n";
    cout << "  " << string(60, '-') << "\n";
    cout << "  Total SV sites printed: " << svPrinted << "\n\n";
}
