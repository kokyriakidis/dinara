
#include "Assembler.hpp"
#include "LongBaseSequence.hpp" 
#include "Reads.hpp"            
#include "timestamp.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <mutex>
#include <cstring> // For memset
#include <immintrin.h> // For AVX2
#include "AlignedEvidenceStore.hpp" 

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
    uint32_t site; // Position on Query Read
    uint32_t occ[4]; // A, C, G, T counts
    uint32_t occ_0;
    uint32_t occ_1; 
    uint32_t fwd_ref_cov; // Ref Count on Strand 0 (Forward). For Strand Bias check.
    char refBase; 
    char altBase; 
    
    int score; // Output of DP
    
    SnpStats() {
        memset(occ, 0, sizeof(occ));
        occ_0 = occ_1 = 0;
        overlap_num = 0;
        refBase = 0; altBase = 0;
        score = 0;
        site = (uint32_t)-1;
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
    vector<uint32_t> siteAnyCov;
    vector<uint64_t> active;

    // Internal helper buffers for detectSVSites
    vector<RawSV> rawSVs;
    vector<size_t> svIndices;
    vector<uint8_t> unpackedRead;
    vector<uint8_t> covered;
    vector<uint32_t> path;

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
        siteAnyCov.clear();
        active.clear();
        rawSVs.clear();
        svIndices.clear();
        unpackedRead.clear();
        covered.clear();
        path.clear();
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
 * Uses the pre-unpacked sequence for maximum speed.
 */
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
 * @brief Robustness filter parity with is_hpc_vec (Correct.cpp:9383).
 * Checks if a site remains informative after removing suspected homopolymer evidence.
 */
static bool isHpcVecMasked(SnpStats& ai, uint32_t id, const vector<HaplotypeEvidence>& hapEvidence, uint32_t s_hap_cov, uint32_t infor_cov) {
    uint32_t n0 = ai.occ_0;
    uint32_t n1 = ai.occ_1;
    
    // Find evidence for this site
    auto it = std::lower_bound(hapEvidence.begin(), hapEvidence.end(), ai.site, [](const HaplotypeEvidence& ev, uint32_t s){
        return ev.site < s;
    });

    while(it != hapEvidence.end() && it->site == ai.site) {
        if(it->hp) {
            if(it->type == 1 && it->overlapSite == id) {
                if(n1 > 0) n1--; 
            }
            // Note: Dinara's hapEvidence currently only stores mismatches (Type 1).
            // Hifiasm's is_hpc_vec also subtracts Ref-HP coverage.
            // However, our mismatch-based occ_0 usually remains stable.
        }
        it++;
    }
    
    if(n0 < 2 || n1 < 2 || n0 < s_hap_cov || n1 < infor_cov) return true;
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

        span<const SnpEvidence> snps;
        uint32_t currentPos = 0; 
        
        // Map SNP coordinates from the alignment to the query read.
        if (ad.readIds[1] == queryReadId) {
            snps = assembler.alignedEvidenceStore.getSnps0(evidenceId);
            currentPos = ad.ts;
        } else if (ad.readIds[0] == queryReadId) {
            snps = assembler.alignedEvidenceStore.getSnps1(evidenceId);
            currentPos = ad.qs;
        } else {
            continue;
        }

        // Apply boundary adjustments (typically 0 for ONT).
        const uint32_t bd = 0; 
        const uint32_t c_qs = cand.qs + bd;
        const uint32_t c_qe = (cand.qe > bd) ? (cand.qe - bd) : cand.qs;
        if (c_qe < c_qs) continue;

        for (const auto& ev : snps) {
            // Jump directly to the next coordinate where a mismatch exists
            currentPos += ev.delta();
            if (currentPos >= c_qs && currentPos < c_qe) {
                hapEvidence.push_back({
                    .overlapID = (uint32_t)k,
                    .site = currentPos,
                    .overlapSite = 0, // Assigned during filtering
                    .type = 1,        // Mismatch
                    .misBase = ev.base(),
                    .hp = false       // Assigned in Phase 3
                });
            }
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
        uint32_t totalMisCount = 0;
        while (evidenceIdx < hapEvidence.size() && hapEvidence[evidenceIdx].site == site) {
            uint8_t base = hapEvidence[evidenceIdx].misBase;
            if (base < 4) {
                misCountPerBase[base]++;
                totalMisCount++;
            }
            evidenceIdx++;
        }
        
        // Homopolymer Masking (Lazy check, once per site)
        bool hp = isHomopolymerMasked(unpacked.data(), (int64_t)queryRead.baseCount, site);
        // Ensure all evidence at this site inherits the HP status for DP filtering
        for (size_t j = startIdx; j < evidenceIdx; ++j) {
            hapEvidence[j].hp = hp;
        }
        if (hp) continue;

        // O(1) Lookup for Coverage
        const uint32_t totalCov = siteTotalCov[siteIdx];
        const uint32_t totalReadsFwd = siteFwdCov[siteIdx];
        const uint32_t refCov = (totalCov >= totalMisCount) ? (totalCov - totalMisCount) : 0;

        // Selection Filters (Parity with Hifiasm push_info)
        if (totalCov < 5) continue; 
        
        // Reference Strand Bias (All Ref on Forward).
        // Calculate altReadsFwd for the WHOLE site to derive refReadsFwd.
        uint32_t totalAltReadsFwd = 0;
        for (size_t j = startIdx; j < evidenceIdx; ++j) {
             const uint32_t ovId = hapEvidence[j].overlapID;
             if (ovId < candidates.size() && !candidates[ovId].isRev) {
                 totalAltReadsFwd++;
             }
        }
        // Hifiasm Parity: is_st_bs (Correct.cpp:9234). 
        // Checks Ref Strand Bias (Extreme case only for Parity).
        if (refCov > 2) {
            uint32_t refReadsFwd = (totalReadsFwd >= totalAltReadsFwd) ? (totalReadsFwd - totalAltReadsFwd) : 0;
            if (refReadsFwd + 2 >= refCov && refReadsFwd >= refCov * 0.95) continue;
            if (refReadsFwd <= 2 && refReadsFwd <= refCov * 0.05) continue;
        }

        // Important: Hifiasm creates a separate SnpStats for EACH alternative base that passes threshold.
        for (uint8_t altBaseIdx = 0; altBaseIdx < 4; ++altBaseIdx) {
            const uint32_t misCount = misCountPerBase[altBaseIdx];
            if (misCount <= 1) continue; // Hifiasm occ_thres = 1

            // Alt Strand Bias Check (Hifiasm parity in push_info)
            uint32_t altReadsFwd = 0;
            for (size_t j = startIdx; j < evidenceIdx; ++j) {
                if (hapEvidence[j].misBase == altBaseIdx) {
                    const uint32_t ovId = hapEvidence[j].overlapID;
                    if (ovId < candidates.size() && !candidates[ovId].isRev) altReadsFwd++;
                }
            }
            if (misCount > 2) { 
                if (altReadsFwd + 2 >= misCount && altReadsFwd >= misCount * 0.95) continue;
                if (altReadsFwd <= 2 && altReadsFwd <= misCount * 0.05) continue;
            }

            // Emit Valid SNP Stats for this specific [Site, AltBase] pair
            SnpStats stat;
            stat.site = site;
            stat.occ_1 = misCount; 
            stat.occ_0 = refCov + 1; // +1 for query support
            stat.fwd_ref_cov = 1; // Simplified for parity
            stat.refBase = queryRead[site].character();
            stat.altBase = Base::fromInteger(altBaseIdx).character();
            
            // Assign overlapSite only for evidence matching this specific alt base
            for (size_t j = startIdx; j < evidenceIdx; ++j) {
                if (hapEvidence[j].misBase == altBaseIdx) {
                    hapEvidence[j].overlapSite = (uint32_t)snpStats.size();
                }
            }
            snpStats.push_back(stat);
        }
    }
}

/**
 * @brief Bit-parallel version of Hifiasm's is_hpc_vec.
 * Checks if a site remains informative after removing HP-suspect evidence.
 */
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
    
    // Parity: Must have at least one supporting read for both alleles (nn[0]>0 && nn[1]>0)
    if (!p0 || !p1) return INT64_MIN;
    return 1; // Hifiasm unweighted DP
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
    size_t evIdx = 0;
    for (size_t i = 0; i < nSites; ++i) {
        const uint32_t s = snpStats[i].site;
        const char alt_c = snpStats[i].altBase;
        uint64_t* row = &flatBits[i * 2 * nWords];
        uint64_t* rowAny = &flatAnyBits[i * nWords];

        while (evIdx < hapEvidence.size() && hapEvidence[evIdx].site < s) evIdx++;
        size_t nextEv = evIdx;
        while (nextEv < hapEvidence.size() && hapEvidence[nextEv].site == s) {
            uint32_t ovId = hapEvidence[nextEv].overlapID;
            if (ovId < nCands) {
                if (Base::fromInteger(hapEvidence[nextEv].misBase).character() == alt_c) {
                    row[2 * (ovId >> 6) + 1] |= (1ULL << (ovId & 63));
                } else {
                    rowAny[ovId >> 6] |= (1ULL << (ovId & 63));
                }
            }
            nextEv++;
        }
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
        const auto& view = assembler.reads->getRead(cand.targetId);
        
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

    // --- Phase 3: Compaction & Index Remapping ---
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
            ev.overlapSite = (newIdx >= 0) ? (uint32_t)newIdx : invalid<uint32_t>();
        } else {
            ev.overlapSite = invalid<uint32_t>();
        }
    }
}

// --------------------------------------------------------
// Function: generate_haplotypes_naive_HiFi
// --------------------------------------------------------
static void generate_haplotypes_naive_HiFi(
    Assembler& assembler,
    HifiasmECScratchPad& scratch
) {
    auto& snpStats = scratch.snpStats;
    auto& hapEvidence = scratch.hapEvidence;
    auto& candidates = scratch.candidates;

    if (snpStats.empty()) return;

    // Reset scores for all candidates
    for (auto& cand : candidates) cand.score = 0;

    // Mark candidates that support any validated SNP in the DP chain.
    for (const auto& ev : hapEvidence) {
        if (ev.overlapSite < snpStats.size()) {
            const auto& s = snpStats[ev.overlapSite];
            // If the SNP at this evidence was validated by DP.
            if (s.score == 1 && ev.overlapID < candidates.size()) {
                candidates[ev.overlapID].score = 1; 
            }
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

    for(size_t k=0; k<candidates.size(); k++) {
        const auto& cand = candidates[k];
        const auto& ad = alignmentData[cand.alignmentId];
        size_t evidenceId = ad.info.alignmentId;
        
        if(evidenceId == invalid<size_t>) continue;

        span<const IndelEvidence> indels;
        
        if(ad.readIds[1] == queryReadId) {
            indels = assembler.alignedEvidenceStore.getIndels0(evidenceId);
        } else if(ad.readIds[0] == queryReadId) {
            indels = assembler.alignedEvidenceStore.getIndels1(evidenceId);
        } else {
            continue;
        }

        for(const auto& ev : indels) {
            if(ev.type() == 0 && ev.len() >= 16) {
                if(ev.pos() >= cand.qs && ev.pos() < cand.qe) {
                     RawSV sv;
                     sv.overlapID = (uint32_t)k;
                     sv.site = ev.pos() + (ev.len() / 2);
                     sv.size = -(int64_t)ev.len(); 
                     rawSVs.push_back(sv);
                }
            }
        }
    }
    
    if(rawSVs.empty()) return;
    
    std::sort(rawSVs.begin(), rawSVs.end(), [](const RawSV& a, const RawSV& b){
        return a.site < b.site;
    });
    
    for(size_t i=0; i<rawSVs.size(); ) {
        uint32_t site = rawSVs[i].site;
        int64_t refSize = rawSVs[i].size; 
        
        uint32_t endWin = site + 50;
        auto& indices = scratch.svIndices;
        indices.clear();
        while(i < rawSVs.size() && rawSVs[i].site < endWin) {
            indices.push_back(i);
            i++;
        }
        
        int validCount = 0;
        for(size_t idx : indices) {
             if(std::abs((double)(rawSVs[idx].size - refSize)) < std::abs((double)refSize)*0.25) { 
                 validCount++;
             }
        }
        
        if(validCount >= 3) {
            SnpStats stat;
            stat.site = site;
            stat.occ_1 = validCount;
            stat.occ_0 = (uint32_t)(candidates.size() - validCount);
            stat.score = -1;
            
            for(size_t idx : indices) {
                int64_t sz = rawSVs[idx].size;
                 if(std::abs((double)(sz - refSize)) < std::abs((double)refSize)*0.25) { 
                     HaplotypeEvidence ev;
                     ev.overlapID = rawSVs[idx].overlapID;
                     ev.site = site;
                     ev.type = 2; // SV
                     ev.overlapSite = (uint32_t)svStats.size();
                     svEvidence.push_back(ev); 
                 }
            }
            svStats.push_back(stat);
        }
    }
}

// --------------------------------------------------------
// Function: generate_haplotypes_sv
// --------------------------------------------------------
static void generate_haplotypes_sv(
    Assembler& assembler,
    HifiasmECScratchPad& scratch
) {
    auto& svEvidence = scratch.svEvidence;
    auto& candidates = scratch.candidates;
    if(svEvidence.empty()) return;
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
    cout << timestamp << "=== Hifiasm Parity EC Pipeline (Round 1) ===" << endl;

    const uint64_t readCount = reads->readCount();
        
    // Use uint8_t for thread-safe byte addressing
    vector<uint8_t> keepAlignment(alignmentData.size(), 0); 
    
    // Parallel loop over reads
    vector<thread> threads;
    uint64_t chunkSize = readCount / threadCount;
    if(chunkSize == 0) chunkSize = 1;

    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            uint64_t start = t * chunkSize;
            uint64_t end = (t == threadCount - 1) ? readCount : (t + 1) * chunkSize;

            // Thread-local scratchpad to eliminate per-read allocations
            HifiasmECScratchPad scratch;

            for(uint64_t readId = start; readId < end; readId++) {
                // Unified Pipeline (TASSD Accelerated)
                uint32_t strand = 0;
                OrientedReadId orientedReadId(dinara::ReadId(readId), strand);
                if(orientedReadId.getValue() >= alignmentTable.size()) continue;

                const auto& alignments = alignmentTable[orientedReadId.getValue()];
                if(alignments.empty()) continue;

                    // Clear scratchpad for next read
                    scratch.clear();
                    auto& candidates = scratch.candidates;
                    candidates.reserve(alignments.size());
                    
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

                    if(candidates.empty()) continue;
                    
                    // 1. SNP Detection Phase
                    // Collect all mismatches from aligned reads and filter them 
                    // based on coverage, strand bias, and homopolymer proximity.
                    detectHetSites(*this, *reads, readId, alignmentData, scratch);
                    
                    // 2. Phasing Phase (DP)
                    // Use Dynamic Programming to find the longest consistent chain of SNPs.
                    // This separates heterozygous sites from sequencing errors.
                    gen_rphase_dp(*this, scratch);
                    
                    // 3. Alignment Validation
                    // Reads that support the validated SNV chain are marked as 'Keep'.
                    generate_haplotypes_naive_HiFi(*this, scratch);

                    // 4. Structural Variant (SV) Phase
                    // Perform a secondary check for large Indels (SVs) to recover links
                    // that SNV-based phasing might miss.
                    detectSVSites(*this, *reads, readId, alignmentData, scratch);
                    generate_haplotypes_sv(*this, scratch);
                    
                    // Final Keep Decision
                    // If an alignment supports either an SNV chain or an SV chain, it is retained.
                    if(!candidates.empty()) {
                         for(const auto& cand : candidates) {
                              if (cand.score > 0) {
                                   keepAlignment[cand.alignmentId] = 1;
                              }
                         }
                    }
                }
        });
    }

    for(auto& t : threads) t.join();

    // Prune Alignment Data
    uint64_t deletedCount = 0;
    for(size_t i=0; i<alignmentData.size(); i++) {
        if(!keepAlignment[i]) {
            if(!alignmentData[i].isDeleted()) {
                 alignmentData[i].setDeleted(true);
                 deletedCount++;
            }
        }
    }
    cout << timestamp << "Parity EC Round 1 Complete. Deleted: " << deletedCount << endl;
}

void Assembler::performHifiasmECFinalFilteringParity(uint64_t threadCount)
{
    cout << timestamp << "=== Hifiasm Parity EC Final Filtering (ha_ec_ff) ===" << endl;
}
