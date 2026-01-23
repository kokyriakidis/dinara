
#include "Assembler.hpp"
#include "LongBaseSequence.hpp" 
#include "Reads.hpp"            
#include "timestamp.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <cstring> // For memset
#include <immintrin.h> // For AVX2
#include "AlignedEvidenceStore.hpp" 
#include "chrono.hpp"

using namespace dinara;
using namespace std;

// Helper to sort candidates
struct CandidateEC {
    uint32_t alignmentId;
    int64_t score;
    uint64_t qs, qe, ts, te;
    
    // DP and Phasing fields
    uint32_t targetId;  // ReadId of target
    bool isRev;         // Orientation
};

// Hifiasm SnpStats equivalent
struct SnpStats {
    uint32_t site;    // Position on Query Read
    uint32_t occ_0;   // Ref support
    uint32_t occ_1;   // Alt support
    uint32_t fwd_ref_cov; // Ref Count on Strand 0 (Forward)
    char refBase; 
    char altBase; 
    int score;        // Output of DP
    
    SnpStats() {
        site = (uint32_t)-1;
        occ_0 = occ_1 = 0;
        fwd_ref_cov = 0;
        refBase = altBase = 0;
        score = 0;
    }
};

struct HaplotypeEvidence {
    uint32_t overlapID; // Index in candidates vector
    uint32_t site;      // Position on Query
    uint32_t overlapSite; // Index into SnpStats array
    uint8_t type;       // 0: ref (match), 1: alt (mismatch), 2: indel?
    uint8_t misBase;    // Base if mismatch
    bool hp;            // is homopolymer suspect (parity with hh_hp)
    
    // Sorting helper
    bool operator<(const HaplotypeEvidence& other) const {
        if(site != other.site) return site < other.site;
        return overlapID < other.overlapID;
    }
};

struct SweepEvent {
    size_t siteIdx;
    uint32_t candIdx;
    bool isEnd;
    bool operator<(const SweepEvent& other) const {
        if(siteIdx != other.siteIdx) return siteIdx < other.siteIdx;
        return isEnd < other.isEnd; // Starts before Ends at same site
    }
};

struct RawSV {
    uint32_t overlapID;
    uint32_t site; 
    int64_t size; 
};

/**
 * @brief Thread-local scratchpad to eliminate per-read allocations in Hifiasm EC.
 * Hoists all vectors to avoid the penalty of repeated mallocs in hot loops.
 */
struct HifiasmECScratchPad {
    vector<CandidateEC> candidates;
    vector<SnpStats> snpStats;
    vector<HaplotypeEvidence> hapEvidence;
    vector<SnpStats> svStats;
    vector<HaplotypeEvidence> svEvidence;
    
    // Internal helper buffers for detectHetSites
    vector<uint32_t> uniqueSites;
    vector<int32_t> diffTotal;
    vector<int32_t> diffFwd;
    vector<uint32_t> siteTotalCov;
    vector<uint32_t> siteFwdCov;
    
    // Internal helper buffers for gen_rphase_dp
    vector<int> validIndices;
    vector<int64_t> f;
    vector<int> p;
    vector<int> indexMap;
    vector<uint64_t> flatBits; // Interleaved: Ref, Alt
    vector<uint64_t> flatAnyBits; // Other mismatches at site
    vector<uint64_t> flatHpBits;
    vector<SweepEvent> events;
    vector<uint64_t> active;

    // Internal helper buffers for detectSVSites
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
        validIndices.clear();
        f.clear();
        p.clear();
        indexMap.clear();
        flatBits.clear();
        flatAnyBits.clear();
        flatHpBits.clear();
        events.clear();
        active.clear();
        rawSVs.clear();
        svIndices.clear();
        unpackedRead.clear();
        covered.clear();
        path.clear();
        supportBits.clear();
        conflictBits.clear();
    }
};

// Helper: Get Sequence Base
inline char getBase(const LongBaseSequenceView& read, uint64_t pos) {
    if(pos >= read.baseCount) return 'N';
    return read[pos].character(); 
}

// Map base to 0-4
inline int base2int(char c) {
    switch(c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

/**
 * @brief Checks if a site is part of a homopolymer run or low-complexity region.
 * Templated to support both raw uint8_t arrays (query) and direct read views (targets).
 */
template<typename SeqType>
static bool isHomopolymerMasked(const SeqType& seq, int64_t len, uint32_t pos) {
    const int64_t hpc_len = 10;   // ASM_HPC_LEN default
    const int64_t hpc_min = 4;    // HPC_RR default
    const int64_t hpc_cut = 6;    // HPC_CC default
    
    // Transparently handle different sequence accessors
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
        // Case 1: Including p, forward check
        for (k = p + r; (k < e) && ((k-r) >= s) && (getBase(k) == getBase(k-r)); k++);
        int64_t ze = k;
        for (k = p - 1; (k >= s) && ((k+r) < e) && (getBase(k) == getBase(k+r)); k--);
        int64_t zs = k + 1;
        if ((ze - zs > r) && (ze - zs >= rc)) return true;
        
        // Case 2: Excluding p, forward check
        for (k = p + r + 1; (k < e) && ((k-r) >= s) && (getBase(k) == getBase(k-r)); k++); 
        ze = k; zs = p + 1;
        if ((ze - zs > r) && (ze - zs >= rc)) return true;
        
        // Case 3: Including p, reverse check
        for (k = p - r; (k >= s) && ((k+r) < e) && (getBase(k) == getBase(k+r)); k--);
        zs = k + 1;
        for (k = p + 1; (k < e) && ((k-r) >= s) && (getBase(k) == getBase(k-r)); k++);
        ze = k;
        if ((ze - zs > r) && (ze - zs >= rc)) return true;

        // Case 4: Excluding p, reverse check
        for (k = p - r - 1; (k >= s) && ((k+r) < e) && (getBase(k) == getBase(k+r)); k--);
        zs = k + 1; ze = p;
        if ((ze - zs > r) && (ze - zs >= rc)) return true;
    }
    return false;
}



/**
 * @brief Identifies potential heterozygous sites (SNPs) by aggregating alignment evidence.
 * Strict parity with Hifiasm's informativeness filters.
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

    // --- Phase 1: Evidence Collection ---
    // In this phase, we iterate through all candidate alignments for the query read.
    // We look for mismatches (SNPs) that have been pre-computed and stored in alignedEvidenceStore.
    // These mismatches are collected as "HaplotypeEvidence".
    
    // Pre-unpack sequence for fast Masking access
    auto& unpacked = scratch.unpackedRead;
    unpacked.resize(queryRead.baseCount);
    for (size_t i = 0; i < queryRead.baseCount; ++i) unpacked[i] = (uint8_t)base2int(queryRead[i].character());

    for (size_t k = 0; k < candidates.size(); ++k) {
        const auto& cand = candidates[k];
        const auto& ad = alignmentData[cand.alignmentId];
        const size_t evidenceId = ad.info.alignmentId;
        
        if (evidenceId == invalid<size_t>) continue;

        // Apply boundary adjustments (typically 0 for ONT).
        const uint32_t bd = 0; 
        const uint32_t c_qs = cand.qs + bd;
        const uint32_t c_qe = (cand.qe > bd) ? (cand.qe - bd) : cand.qs;
        if (c_qe < c_qs) continue;

        // Map SNP coordinates from the alignment to the query read and iterate in range.
        const uint32_t evidenceId32 = uint32_t(evidenceId);
        auto addEv = [&](uint32_t currentPos, uint8_t base) {
            if (currentPos < unpacked.size() && base == unpacked[currentPos]) {
                return;
            }
            hapEvidence.push_back({
                .overlapID = (uint32_t)k,
                .site = currentPos,
                .overlapSite = invalid<uint32_t>, // Assigned during filtering
                .type = 1,        // Mismatch
                .misBase = base,
                .hp = false       // Assigned in Phase 3
            });
        };

        if (ad.readIds[1] == queryReadId) {
            assembler.alignedEvidenceStore.forEachSnp0InRange(evidenceId32, c_qs, c_qe, addEv);
        } else if (ad.readIds[0] == queryReadId) {
            assembler.alignedEvidenceStore.forEachSnp1InRange(evidenceId32, c_qs, c_qe, addEv);
        } else {
            continue;
        }
    }
    
    // Sort all evidence by position on the query read.
    std::sort(hapEvidence.begin(), hapEvidence.end());
    if (hapEvidence.empty()) return;

    // Identify unique genomic sites (positions on the query read) where mismatches were seen.
    auto& uniqueSites = scratch.uniqueSites;
    uniqueSites.clear();
    for (size_t i = 0; i < hapEvidence.size(); ) {
        uniqueSites.push_back(hapEvidence[i].site);
        uint32_t current = hapEvidence[i].site;
        while (i < hapEvidence.size() && hapEvidence[i].site == current) i++;
    }
    
    // --- Phase 2: Site-Space Coverage Aggregation ---
    // We need to know how many reads cover each site to calculate allele frequencies.
    // Instead of iterating O(Reads * Sites), we use a sweep-line algorithm (diffTotal/diffFwd)
    // to calculate coverage at every unique site in O(Reads + Sites).
    const size_t nSites = uniqueSites.size();
    auto& diffTotal = scratch.diffTotal;
    auto& diffFwd = scratch.diffFwd;
    diffTotal.assign(nSites + 1, 0);
    diffFwd.assign(nSites + 1, 0);
    const uint32_t bd = 0; 
    
    for (const auto& cand : candidates) {
        const uint32_t s = cand.qs + bd;
        const uint32_t e = (cand.qe > bd) ? (cand.qe - bd) : cand.qs;
        if (e <= s) continue;

        // Find the range of unique sites covered by this candidate.
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
    }

    // Prefix sum to finalize coverage counts.
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

    // --- Phase 3: Site Processing & Filtering ---
    // Now we decide which sites are likely heterozygous (SNVs).
    // Filters:
    // 1. Coverage > threshold (min 5 reads).
    // 2. Homopolymer masking: Sites in HP runs are filtered to avoid noise.
    // 3. Strand bias: Sites where alleles only appear on one strand are risky.
    for (size_t siteIdx = 0, evidenceIdx = 0; siteIdx < nSites; ++siteIdx) {
        const uint32_t site = uniqueSites[siteIdx];
        const uint32_t startIdx = (uint32_t)evidenceIdx;
        
        // Collect all evidence for this specific site into base-specific counts.
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

        // Quick reject: no allele has enough support to become an informative site.
        uint32_t maxMisCount = 0;
        for (uint8_t b = 0; b < 4; ++b) {
            maxMisCount = std::max(maxMisCount, misCountPerBase[b]);
        }
        if (maxMisCount < 2) {
            continue;
        }
        
        // O(1) Lookup for Coverage
        const uint32_t totalCov = siteTotalCov[siteIdx];
        const uint32_t totalReadsFwd = siteFwdCov[siteIdx];
        const uint32_t refCov = (totalCov >= totalMisCount) ? (totalCov - totalMisCount) : 0;

        // Selection Filters (Parity with Hifiasm push_info)
        if (totalCov < 5) continue; 
        
        // Reference Strand Bias (All Ref on Forward).
        // Hifiasm Parity: is_st_bs (Correct.cpp:9234). 
        // Checks Ref Strand Bias (Extreme case only for Parity).
        if (refCov > 2) {
            uint32_t refReadsFwd = (totalReadsFwd >= totalAltReadsFwd) ? (totalReadsFwd - totalAltReadsFwd) : 0;
            // Integer math: x >= 0.95*y  <=>  100*x >= 95*y
            if ((refReadsFwd + 2 >= refCov) && (uint64_t(refReadsFwd) * 100ULL >= uint64_t(refCov) * 95ULL)) continue;
            // Integer math: x <= 0.05*y  <=>  100*x <= 5*y
            if ((refReadsFwd <= 2) && (uint64_t(refReadsFwd) * 100ULL <= uint64_t(refCov) * 5ULL)) continue;
        }

        // Determine which alternative alleles are worth emitting for this site.
        uint32_t rowIndexForBase[4] = {
            invalid<uint32_t>, invalid<uint32_t>, invalid<uint32_t>, invalid<uint32_t>
        };
        uint32_t validAltCount = 0;

        // Important: Hifiasm creates a separate SnpStats for EACH alternative base that passes threshold.
        for (uint8_t altBaseIdx = 0; altBaseIdx < 4; ++altBaseIdx) {
            const uint32_t misCount = misCountPerBase[altBaseIdx];
            if (misCount < 2) continue; // Initial detection threshold (DP later requires 3)

            // Alt Strand Bias Check (Hifiasm parity in push_info)
            const uint32_t altReadsFwd = misCountFwdPerBase[altBaseIdx];

            // --- Hifiasm Parity: Strand Bias Filtering ---
            // Filter Alt Allele Bias
            if (misCount > 2) { 
                if ((altReadsFwd + 2 >= misCount) && (uint64_t(altReadsFwd) * 100ULL >= uint64_t(misCount) * 95ULL)) continue;
                if ((altReadsFwd <= 2) && (uint64_t(altReadsFwd) * 100ULL <= uint64_t(misCount) * 5ULL)) continue;
            }

            // Note: Site-level Ref Bias already checked above (lines 426-430).
            // No need for a second check here unless we have extremely high multiallelic complexity.

            // Defer emission until after homopolymer masking for this site.
            rowIndexForBase[altBaseIdx] = 0;
            validAltCount++;
        }

        if (validAltCount == 0) {
            continue;
        }

        // Homopolymer masking is relatively expensive. Do it only after cheaper filters passed.
        if (isHomopolymerMasked(unpacked.data(), (int64_t)queryRead.baseCount, site)) {
            continue;
        }

        const uint32_t refReadsFwdSite = (totalReadsFwd >= totalAltReadsFwd) ? (totalReadsFwd - totalAltReadsFwd) : 0;

        // Emit Valid SNP Stats for each surviving [Site, AltBase] pair and record its row index.
        for (uint8_t altBaseIdx = 0; altBaseIdx < 4; ++altBaseIdx) {
            if (rowIndexForBase[altBaseIdx] == invalid<uint32_t>) {
                continue;
            }

            const uint32_t misCount = misCountPerBase[altBaseIdx];

            // Emit Valid SNP Stats for this specific [Site, AltBase] pair
            SnpStats stat;
            stat.site = site;
            stat.occ_1 = misCount; 
            stat.occ_0 = refCov + 1; // +1 for query support
            stat.fwd_ref_cov = refReadsFwdSite + 1; // +1 for query support (Forward by convention)
            stat.refBase = queryRead[site].character();
            stat.altBase = Base::fromInteger(altBaseIdx).character();

            rowIndexForBase[altBaseIdx] = (uint32_t)snpStats.size();
            snpStats.push_back(stat);
        }

        // Assign overlapSite for evidence that matches an emitted alt base at this site.
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

/**
 * @brief Bit-parallel version of Hifiasm's is_hpc_vec.
 * Checks if a site remains informative after removing HP-suspect evidence.
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
    
    uint64_t n0_hp = 0, n1_hp = 0;
    size_t w = 0;
    // Vectorized popcount accumulation (if nWords > 4)
    for (; w < nWords; ++w) {
        n0_hp += __builtin_popcountll(rowBits[2*w] & rowHp[w]);
        n1_hp += __builtin_popcountll(rowBits[2*w+1] & rowHp[w]);
    }
    
    uint32_t n0_robust = (stat.occ_0 >= n0_hp) ? (uint32_t)(stat.occ_0 - n0_hp) : 0;
    uint32_t n1_robust = (stat.occ_1 >= n1_hp) ? (uint32_t)(stat.occ_1 - n1_hp) : 0;
    
    return (n0_robust < 2 || n1_robust < 2 || n0_robust < s_hap_cov || n1_robust < infor_cov);
}

/**
 * @brief Strictly calculates the phasing link score between two SNPs. 
 * Parity with Hifiasm: link is invalid if any read in overlap has a "Third Allele" contradiction.
 * Exception (Hifiasm Correct.cpp:9270): if both sites see a third allele, it's treated as Ref-Ref support.
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

        // Reads in overlap seeing 'Other' at I and J simultaneously (Hifiasm parity)
        const uint64_t rareRef = oI & oJ;
        
        // Corrected Alt/Other/Ref status
        const uint64_t effRefI = rI | rareRef;
        const uint64_t effRefJ = rJ | rareRef;

        // Condition 1: Strict Allele Consistency.
        // A phased link must be "Clean". If a read covers both sites but shows
        // a third allele at exactly one site (invalidating the biallelic assumption),
        // the entire link is discarded (Parity with Hifiasm Correct.cpp:9270).
        const uint64_t soloOther = (oI ^ oJ) & (rI | aI | oI) & (rJ | aJ | oJ);
        if (soloOther) return INT64_MIN;

        // Condition 2: No Ref-Alt or Alt-Ref contradictions.
        if (effRefI & aJ) return INT64_MIN;
        if (aI & effRefJ) return INT64_MIN;

        p0 |= (effRefI & effRefJ);
        p1 |= (aI & aJ);
    }
    
    // Parity check: Link MUST be supported on both haplotypes
    if (!p0 || !p1) return INT64_MIN;

    // Hifiasm Weighted Scoring: Sum of popcounts of Ref-Ref and Alt-Alt support.
    // This favors links with higher overlap depth.
    uint64_t score = 0;
    for (size_t k = 0; k < nWords; ++k) {
        const uint64_t rI = rowI[2*k], aI = rowI[2*k+1], oI = anyI[k];
        const uint64_t rJ = rowJ[2*k], aJ = rowJ[2*k+1], oJ = anyJ[k];
        const uint64_t rareRef = oI & oJ;
        score += __builtin_popcountll((rI | rareRef) & (rJ | rareRef));
        score += __builtin_popcountll(aI & aJ);
    }
    return (int64_t)score;
}

// --------------------------------------------------------
// Function: gen_rphase_dp
// --------------------------------------------------------
/**
 * @brief Dynamic Programming for phasing heterozygous sites. 
 * Replicates Hifiasm's gen_rphase_dp0_single_path with bit-parallel speedups.
 */
static void gen_rphase_dp(
    Assembler& assembler,
    HifiasmECScratchPad& scratch
) {
    auto& snpStats = scratch.snpStats;
    auto& hapEvidence = scratch.hapEvidence;
    auto& candidates = scratch.candidates;

    if (snpStats.empty() || hapEvidence.empty()) return;

    const size_t nCands = candidates.size();
    const size_t nSites = snpStats.size();
    const size_t nWords = (nCands + 63) / 64;

    // --- Phase 1: Matrix Construction (Interleaved for Locality) ---
    // We represent the presence of Ref, Alt, or 'Other' alleles in each read
    // using bitsets (64 reads per uint64_t). This allows O(N/64) phasing checks.
    auto& flatBits = scratch.flatBits;
    auto& flatAnyBits = scratch.flatAnyBits;
    flatBits.assign(nSites * 2 * nWords, 0); 
    flatAnyBits.assign(nSites * nWords, 0);
    
    // Populate Alt and AnyOther Bits in a single pass.
    // flatBits[2*w] = Ref alleles, flatBits[2*w+1] = Alt alleles.
    // flatAnyBits[w] = Any bases that are neither Ref nor Alt (noise/third alleles).
    auto& allMisBits = scratch.supportBits; // Temporary; overwritten later.
    size_t evIdx = 0;
    for (size_t i = 0; i < nSites; ) {
        const uint32_t s = snpStats[i].site;

        // Group all [site, alt] rows for this site to avoid re-scanning hapEvidence
        // once per alt allele.
        size_t j = i + 1;
        while (j < nSites && snpStats[j].site == s) ++j;

        // Map alt base (0..3) -> row index in [i, j) if present.
        int rowForBase[4] = {-1, -1, -1, -1};
        for (size_t r = i; r < j; ++r) {
            const int b = base2int(snpStats[r].altBase);
            if (b >= 0 && b < 4) rowForBase[b] = int(r);
        }

        allMisBits.assign(nWords, 0);

        // Find evidence at this site.
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

        // For each [site, alt] row, "AnyOther" are mismatches not equal to that alt.
        for (size_t r = i; r < j; ++r) {
            const uint64_t* row = &flatBits[r * 2 * nWords];
            uint64_t* rowAny = &flatAnyBits[r * nWords];
            for (size_t w = 0; w < nWords; ++w) {
                rowAny[w] = allMisBits[w] & ~row[2 * w + 1];
            }
        }

        i = j;
    }

    // Populate Ref Bits using a Sweep-line approach.
    // A read is 'Ref' at a site if it covers the site AND has no mismatches there.
    auto& uniqueSites = scratch.uniqueSites;
    uniqueSites.clear();
    for (const auto& s : snpStats) uniqueSites.push_back(s.site);

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

    auto& active = scratch.active;
    active.assign(nWords, 0);
    size_t eventIdx = 0;
    for (size_t i = 0; i < nSites; ++i) {
        while (eventIdx < events.size() && events[eventIdx].siteIdx == i) {
            uint32_t c = events[eventIdx].candIdx;
            if (events[eventIdx].isEnd) active[c >> 6] &= ~(1ULL << (c & 63));
            else active[c >> 6] |= (1ULL << (c & 63));
            eventIdx++;
        }
        
        uint64_t* row = &flatBits[i * 2 * nWords];
        const uint64_t* rowAny = &flatAnyBits[i * nWords]; 
        for (size_t w = 0; w < nWords; ++w) {
            // Interleaved Ref = Covered (Active) AND NOT (PrimaryAlt OR AnyOther)
            uint64_t allMis = rowAny[w] | row[2*w+1];
            row[2*w] = active[w] & ~allMis;
        }
    }

    // Populate flatHpBits (Bit-parallel HP mask)
    // Optimized: No target read unpacking, binary search for covered sites.
    auto& flatHpBits = scratch.flatHpBits;
    flatHpBits.assign(nSites * nWords, 0);
    
    for (size_t candIdx = 0; candIdx < nCands; ++candIdx) {
        const auto& cand = candidates[candIdx];
        const auto& view = assembler.getReads().getRead(cand.targetId);
        
        // Only iterate over SNPs that this candidate actually covers
        auto itStart = std::lower_bound(uniqueSites.begin(), uniqueSites.end(), cand.qs);
        auto itEnd = std::lower_bound(uniqueSites.begin(), uniqueSites.end(), cand.qe);
        size_t idxS = std::distance(uniqueSites.begin(), itStart);
        size_t idxE = std::distance(uniqueSites.begin(), itEnd);

        for (size_t i = idxS; i < idxE; ++i) {
            uint32_t site = uniqueSites[i];
            uint32_t targetSite = cand.isRev ? (uint32_t)(cand.te - (site - cand.qs) - 1) : (uint32_t)(cand.ts + (site - cand.qs));
            // Zero-unpacking HP check
            if (isHomopolymerMasked(view, (int64_t)view.baseCount, targetSite)) {
                flatHpBits[i * nWords + (candIdx >> 6)] |= (1ULL << (candIdx & 63));
            }
        }
    }

    // --- Phase 2: DP Filtering & Scoring ---
    // This is the core phasing logic. We use dynamic programming to find a continuous chain
    // of SNVs that are mutually consistent across the supporting reads.
    const uint32_t s_hap_cov = 3;
    const uint32_t infor_cov = 3;
    const uint32_t cc = 2; 

    auto& validIndices = scratch.validIndices;
    validIndices.clear();
    for (size_t i = 0; i < nSites; ++i) {
        if (snpStats[i].occ_0 >= s_hap_cov && snpStats[i].occ_1 >= infor_cov) {
            validIndices.push_back((int)i);
        }
        snpStats[i].score = -1;
    }

    if (validIndices.empty()) return;

    // DP formulation: f[i] = max length of a phased chain ending at SNV i.
    const int nV = (int)validIndices.size();
    auto& f = scratch.f; f.assign(nV, 0);
    auto& p = scratch.p; p.assign(nV, -1);
    auto& covered = scratch.covered; covered.assign(nV, 0);

    for (int i = 0; i < nV; ++i) {
        int siteI = validIndices[i];
        if (snpStats[siteI].occ_0 < cc) continue; 
        
        // Secondary Robustness Check
        if (isHpcVecMaskedBitParallel(snpStats[siteI], siteI, flatBits.data(), flatHpBits.data(), nWords, s_hap_cov, infor_cov)) continue;

        f[i] = 1;
        for (int j = 0; j < i; ++j) {
            int siteJ = validIndices[j];
            // Check if SNV i and SNV j can belong to the same haplotype.
            int64_t sc = comput_sc_rphase_strict(flatBits.data(), flatAnyBits.data(), siteI, siteJ, nWords);
            if (sc == INT64_MIN) continue;
            
            if (f[j] + (int)sc > f[i]) {
                f[i] = f[j] + (int)sc;
                p[i] = j;
            }
        }
    }

    // Greedy extract disjoint paths (Harvest uncovered segments)
    vector<pair<int64_t, int>> scoreIdx;
    for(int i=0; i<nV; ++i) scoreIdx.push_back({f[i], i});
    std::sort(scoreIdx.rbegin(), scoreIdx.rend());

    for(auto& pair : scoreIdx) {
        int curr = pair.second;
        if(f[curr] < 1 || covered[curr]) continue;
        
        // Trace back and harvest only the uncovered prefix
        while(curr != -1 && !covered[curr]) {
            covered[curr] = 1;
            int snpIdx = validIndices[curr];
            if (snpStats[snpIdx].occ_0 >= cc) {
                snpStats[snpIdx].score = 1;
            }
            curr = p[curr];
        }
    }
}

// --------------------------------------------------------
// Function: generate_haplotypes_naive_HiFi
// --------------------------------------------------------
/**
 * @brief High-performance bit-parallel alignment validation.
 * Keeps an alignment if it consistently supports exactly one haplotype
 * across all validated DP-chain SNVs, with zero noise alleles.
 */
static void generate_haplotypes_naive_HiFi(
    Assembler& /* assembler */,
    HifiasmECScratchPad& scratch
) {
    const size_t nCands = scratch.candidates.size();
    const size_t nWords = (nCands + 63) / 64;
    const size_t nSites = scratch.snpStats.size();

    // Reset scores
    for (auto& cand : scratch.candidates) cand.score = 0;

    auto& supportRef = scratch.supportBits;
    auto& supportAlt = scratch.conflictBits;
    auto& globalOther = scratch.active; // Reuse scratchpad buffer
    supportRef.assign(nWords, 0);
    supportAlt.assign(nWords, 0);
    globalOther.assign(nWords, 0);

    // Phase 1: Global support/conflict aggregation (O(ValidatedSites * Words))
    for (size_t i = 0; i < nSites; ++i) {
        if (scratch.snpStats[i].score != 1) continue;
        
        const uint64_t* row = &scratch.flatBits[i * 2 * nWords];
        const uint64_t* rowAny = &scratch.flatAnyBits[i * nWords];
        
        for (size_t w = 0; w < nWords; ++w) {
            supportRef[w] |= row[2 * w];
            supportAlt[w] |= row[2 * w + 1];
            globalOther[w] |= rowAny[w];
        }
    }

    // Phase 2: Alignment-level keep decision (O(Candidates))
    for (size_t c = 0; c < nCands; ++c) {
        size_t w = c >> 6;
        uint64_t mask = (1ULL << (c & 63));
        bool supportsRef = supportRef[w] & mask;
        bool supportsAlt = supportAlt[w] & mask;
        bool hasOther = globalOther[w] & mask;

        // --- Hifiasm Parity Filter ---
        // An alignment is KEPT if it does not contradict the phased haplotypes.
        // Contradictions are:
        // 1. Bridges both haplotypes: (supportsRef AND supportsAlt)
        // 2. Contains noise/third alleles at validated sites: (hasOther)
        
        bool isConflict = (supportsRef && supportsAlt) || hasOther;
        
        if (!isConflict) {
            scratch.candidates[c].score = 1;
        }
    }
}

/**
 * @brief Finalizes the SNV set and evidence mappings after validation.
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

// --------------------------------------------------------
// Helper: detectSVSites
// --------------------------------------------------------
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
    rawSVs.clear();
    svEvidence.clear();
    svStats.clear();

    const int32_t SV_MIN_LEN = 20;     // Hifiasm standard for SV-based recovery
    const int32_t SV_WINDOW = 50;      // Window for clustering nearby indels
    const double SV_SIZE_RATIO = 0.20; // 20% size variation allowed in a cluster

    // --- Phase 1: SV Collection ---
    for(size_t k = 0; k < candidates.size(); ++k) {
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

        // Evidence is expected to be non-decreasing within each per-alignment slice.
        auto it = std::lower_bound(
            indels.begin(), indels.end(), c_qs,
            [](const IndelEvidence& e, uint32_t value) { return e.pos() < value; }
        );
        for (; it != indels.end(); ++it) {
            const uint32_t pos = it->pos();
            if (pos >= c_qe) break;
            const uint32_t len = it->len();
            if (len < (uint32_t)SV_MIN_LEN) continue;

            RawSV sv;
            sv.overlapID = (uint32_t)k;
            sv.site = pos;
            // Use sign to distinguish Ins (+) from Del (-)
            sv.size = it->isInsertion() ? (int64_t)len : -(int64_t)len;
            rawSVs.push_back(sv);
        }
    }
    
    if(rawSVs.empty()) return;
    
    // Sort by site for windowed clustering
    std::sort(rawSVs.begin(), rawSVs.end(), [](const RawSV& a, const RawSV& b){
        return a.site < b.site;
    });
    
    // --- Phase 2: Windowed Clustering & Frequency Validation ---
    for(size_t i = 0; i < rawSVs.size(); ) {
        uint32_t startPos = rawSVs[i].site;
        int64_t refSize = rawSVs[i].size; 
        
        size_t j = i;
        int supportCount = 0;
        
        // Peek ahead to find support within window and size ratio
        while(j < rawSVs.size() && rawSVs[j].site < startPos + SV_WINDOW) {
            int64_t sz = rawSVs[j].size;
            // Check if same type (sign) and similar size
            if ((sz > 0) == (refSize > 0)) {
                int64_t diff = std::abs(sz - refSize);
                if (diff <= std::abs(refSize) * SV_SIZE_RATIO) {
                    supportCount++;
                }
            }
            j++;
        }
        
        // Hifiasm Parity: Recover alignment if SV has >= 3 supporting reads and strand balance
        if(supportCount >= 3) {
            // Check Strand Bias for SV (Parity with ha_ec_st_bs for SVs)
            int svFwdCount = 0;
            for(size_t k = i; k < j; ++k) {
                int64_t sz = rawSVs[k].size;
                if ((sz > 0) == (refSize > 0) && std::abs(sz - refSize) <= std::abs(refSize) * SV_SIZE_RATIO) {
                    if (!candidates[rawSVs[k].overlapID].isRev) svFwdCount++;
                }
            }
            
            // Refined bias thresholds matching ha_ec_st_bs
            if (supportCount > 5) {
                if (svFwdCount < 1 || svFwdCount >= supportCount) { // High confidence mode
                     if (svFwdCount < supportCount * 0.05 || svFwdCount > supportCount * 0.95) {
                         i++; continue; 
                     }
                }
            } else if (svFwdCount == 0 || svFwdCount == supportCount) {
                 // Low coverage SVs MUST be seen on both strands
                 i++; continue;
            }

            SnpStats stat;
            stat.site = startPos;
            stat.score = 1; // Mark as "Validated SV"
            
            // Re-scan to populate evidence for this specific validated SV cluster
            for(size_t k = i; k < j; ++k) {
                int64_t sz = rawSVs[k].size;
                if ((sz > 0) == (refSize > 0) && std::abs(sz - refSize) <= std::abs(refSize) * SV_SIZE_RATIO) {
                    HaplotypeEvidence ev;
                    ev.overlapID = rawSVs[k].overlapID;
                    ev.site = startPos;
                    ev.type = 2; // Type SV
                    ev.overlapSite = (uint32_t)svStats.size();
                    svEvidence.push_back(ev);
                }
            }
            svStats.push_back(stat);
            i = j; // Move to next cluster
        } else {
            i++; // Move to next potential seed
        }
    }
}

// --------------------------------------------------------
// Function: generate_haplotypes_sv
// --------------------------------------------------------
static void generate_haplotypes_sv(
    Assembler& /* assembler */,
    HifiasmECScratchPad& scratch
) {
    auto& svEvidence = scratch.svEvidence;
    auto& candidates = scratch.candidates;
    
    // Simple recovery: If a read supports a validated, frequent SV, we keep it.
    // This acts as a fallback for reads where SNV phasing was too sparse or fragmented.
    for(const auto& ev : svEvidence) {
        if(ev.overlapID < candidates.size()) {
            candidates[ev.overlapID].score = 1; 
        }
    }
}



// --------------------------------------------------------
// Main Function: performHifiasmECParity
// --------------------------------------------------------
void Assembler::performHifiasmECParity(uint64_t threadCount)
{ 
    static std::atomic<int> debugInformativeCount{0};
    cout << timestamp << "=== Hifiasm Parity EC Pipeline (Round 1) ===" << endl;

    const uint64_t readCount = reads->readCount();
    const auto tBeginAll = steady_clock::now();
        
    // Use uint8_t for thread-safe byte addressing
    vector<uint8_t> keepAlignment(alignmentData.size(), 0); 
    
    // Parallel loop over reads
    vector<thread> threads;
    uint64_t chunkSize = readCount / threadCount;
    if(chunkSize == 0) chunkSize = 1;

    struct alignas(64) PhaseTiming {
        double gatherCandidates = 0.;
        double snpDetect = 0.;
        double dp = 0.;
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

            // Thread-local scratchpad to eliminate per-read allocations
            HifiasmECScratchPad scratch;
            PhaseTiming& timing = timings[t];

            for(uint64_t readId = start; readId < end; readId++) {
                timing.readsVisited++;
                // Unified Pipeline (TASSD Accelerated)
                uint32_t strand = 0;
                OrientedReadId orientedReadId(dinara::ReadId(readId), strand);
                if(orientedReadId.getValue() >= alignmentTable.size()) continue;

                const auto& alignments = alignmentTable[orientedReadId.getValue()];
                if(alignments.empty()) continue;
                timing.readsWithAlignments++;

                    // Clear scratchpad for next read
                    scratch.clear();
                    auto& candidates = scratch.candidates;
                    candidates.reserve(alignments.size());
                    
                    const auto tGatherBegin = steady_clock::now();
                    for(uint32_t alignmentId : alignments) {
                        const auto& thisAlignmentData = alignmentData[alignmentId];
                        
                        // Skip deleted/filtered alignments
                        if(thisAlignmentData.isDeleted()) continue;

                        CandidateEC candidate;
                        candidate.alignmentId = alignmentId;
                        
                        if(readId == thisAlignmentData.readIds[0]) {
                            candidate.targetId = thisAlignmentData.readIds[1];
                            candidate.qs = thisAlignmentData.qs; 
                            candidate.qe = thisAlignmentData.qe;
                            candidate.ts = thisAlignmentData.ts;
                            candidate.te = thisAlignmentData.te;
                            candidate.isRev = !thisAlignmentData.isSameStrand; 
                        } else {
                            candidate.targetId = thisAlignmentData.readIds[0];
                            candidate.qs = thisAlignmentData.ts; 
                            candidate.qe = thisAlignmentData.te; 
                            candidate.ts = thisAlignmentData.qs;
                            candidate.te = thisAlignmentData.qe; 
                            candidate.isRev = !thisAlignmentData.isSameStrand;
                        }
                        candidate.score = 0; // Reset for EC keep decision
                        
                        candidates.push_back(candidate);
                    }
                    timing.gatherCandidates += seconds(steady_clock::now() - tGatherBegin);

                    if(candidates.empty()) continue;
                    timing.readsWithCandidates++;
                    
                    bool isDebugRead = false;
                    
                    // 1. SNP Detection Phase
                    const auto tSnpBegin = steady_clock::now();
                    detectHetSites(*this, *reads, readId, alignmentData, scratch);
                    timing.snpDetect += seconds(steady_clock::now() - tSnpBegin);
                    
                    // 2. Phasing Phase (DP)
                    const auto tDpBegin = steady_clock::now();
                    gen_rphase_dp(*this, scratch);
                    timing.dp += seconds(steady_clock::now() - tDpBegin);
                    
                    // 3. Alignment Validation
                    const auto tValidateSnvBegin = steady_clock::now();
                    generate_haplotypes_naive_HiFi(*this, scratch);
                    timing.validateSnv += seconds(steady_clock::now() - tValidateSnvBegin);

                    // Pick first 2 informative reads for debug tracking
                    if (debugInformativeCount < 2) {
                        bool hasSnvChain = false;
                        for (const auto& s : scratch.snpStats) if (s.score == 1) { hasSnvChain = true; break; }
                        if (hasSnvChain) {
                             if (debugInformativeCount.fetch_add(1) < 2) isDebugRead = true;
                        }
                    }
                    
                    vector<size_t> filteredBySnv;
                    if (isDebugRead) {
                        cout << "[DEBUG] Transitions for Informative Read " << readId << ":" << endl;
                    }
                    
                    if (isDebugRead) {
                        cout << "  - Filtered by SNV Phasing (no valid chain support):" << endl;
                        for (size_t c = 0; c < candidates.size(); ++c) {
                            if (candidates[c].score == 0) {
                                cout << "    * Filtered: " << candidates[c].targetId << ", Strand: " << (candidates[c].isRev ? "R" : "F") << endl;
                                filteredBySnv.push_back(c);
                            }
                        }
                        if (filteredBySnv.empty()) cout << "    * (None)" << endl;
                    }

                    const auto tCompactBegin = steady_clock::now();
                    compactPhasedSites(scratch);
                    timing.compact += seconds(steady_clock::now() - tCompactBegin);

                    // 4. Structural Variant (SV) Phase
                    const auto tSvDetectBegin = steady_clock::now();
                    detectSVSites(*this, *reads, readId, alignmentData, scratch);
                    timing.svDetect += seconds(steady_clock::now() - tSvDetectBegin);

                    const auto tValidateSvBegin = steady_clock::now();
                    generate_haplotypes_sv(*this, scratch);
                    timing.validateSv += seconds(steady_clock::now() - tValidateSvBegin);

                    if (isDebugRead) {
                        cout << "  - Recovered by SV Detection (backbone support):" << endl;
                        bool anyRecovered = false;
                        for (size_t c : filteredBySnv) {
                            if (candidates[c].score > 0) {
                                cout << "    * Recovered: " << candidates[c].targetId << ", Strand: " << (candidates[c].isRev ? "R" : "F") << endl;
                                anyRecovered = true;
                            }
                        }
                        if (!anyRecovered) cout << "    * (None)" << endl;
                    }
                    
                    // --- Connectivity Heuristic ---
                    // A read is "informative" if it has at least one SNV or SV 
                    // that survived filtering and was validated (score == 1).
                    // compactPhasedSites already pruned snpStats to score == 1 only.
                    const bool isInformativeRead = !scratch.snpStats.empty() || !scratch.svStats.empty();

                    // Pre-calculate which candidates cover at least one informative site
                    // using a local bitset for O(N_Evidence) instead of O(N_Cands * N_Evidence).
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

                    // Final Keep Decision & Flag Management
                    if (isDebugRead) {
                        cout << "  - Final Filtered Set (Informative=" << (isInformativeRead?"Yes":"No") << "):" << endl;
                        bool anyFinalDeletions = false;
                        for (const auto& cand : candidates) {
                            bool keep = !isInformativeRead || (cand.score > 0);
                            if (!keep) {
                                cout << "    * DELETED: " << cand.targetId << ", Strand: " << (cand.isRev ? "R" : "F") << endl;
                                anyFinalDeletions = true;
                            }
                        }
                        if (!anyFinalDeletions) cout << "    * (None)" << endl;
                    }

                    const auto tFinalizeBegin = steady_clock::now();
                    for(size_t c = 0; c < candidates.size(); ++c) {
                        auto& cand = candidates[c];
                        auto& ad = alignmentData[cand.alignmentId];
                        
                        // Decision Logic:
                        // 1. If the read is NOT informative (homozygous/repetitive), we KEEP all overlaps.
                        // 2. If the read IS informative, we only KEEP overlaps that don't conflict (cand.score == 1).
                        bool keep = !isInformativeRead || (cand.score > 0);

                        if (keep) ad.info.isInReadGraph = 1;

                        // Set Directional Deletion Flags (readIds[0] <= readIds[1] invariant)
                        if (readId == ad.readIds[0]) {
                            ad.isDeleted0 = !keep;
                        } else {
                            ad.isDeleted1 = !keep;
                        }

                        // Set Assembly Flags (shared across threads, logical OR)
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
