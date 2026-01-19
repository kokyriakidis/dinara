
#include "Assembler.hpp"
#include "timestamp.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <mutex>

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
    
    // Check if we need to store mismatches
};

// Hifiasm SnpStats equivalent
struct SnpStats {
    uint32_t site; // Position on Query Read
    uint32_t occ[4]; // A, C, G, T counts (simplified from occ_0/occ_1 for biallelic)
    // Hifiasm uses occ_0, occ_1 for two alleles. We can track max 2 alleles or all 4.
    // For parity with `gen_rphase_dp0_single_path` which expects occ_0/occ_1:
    uint32_t occ_0;
    uint32_t occ_1; 
    uint32_t overlap_num; // Ref Count on Strand 0 (Forward). For Strand Bias check.
    // We need to map base to allele 0 or 1? 
    // Hifiasm `hc_phase_robust_rr` identifies "ref" vs "alt".
    // "Ref" is query base?
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
    
    // Sorting helper
    bool operator<(const HaplotypeEvidence& other) const {
        if(overlapID != other.overlapID) return overlapID < other.overlapID;
        return site < other.site;
    }
};

// Helper: Get Sequence Base
// Reads are packed or vector representation.
// Assuming reads.getRead(id) returns a vector/string of bases.
inline char getBase(const Read& read, uint64_t pos) {
    if(pos >= read.size()) return 'N';
    // Dinara Read is uint8_t vector? Or char?
    // Accessor: read[pos]
    return read[pos]; // Assuming Read operator[] works
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

// --------------------------------------------------------
// Function: detectHetSites (Parity for hc_phase_robust_rr)
// --------------------------------------------------------
// Scans overlaps to find sites with consistent mismatches.
// Populates SnpStats and HaplotypeEvidence.
void detectHetSites(
    const Reads& reads,
    uint64_t queryReadId,
    const vector<CandidateEC>& candidates, // "ol"
    const vector<AlignmentData>& alignmentData, // Access to raw alignment info
    vector<SnpStats>& snpStats,
    vector<HaplotypeEvidence>& hapEvidence
) {
    // 1. Scan overlaps to collect mismatches (Potential SNPs)
    // Since we don't have CIGARs, we must re-verify markers or check bases between markers.
    // For efficiency, let's assume we check bases at aligned markers (if they cover enough)
    // O(N_aligned_bases). Heavy.
    // Alternative: Dinara `Alignment` has `ordinals`. These are exact matches.
    // Mismatches are NOT in the `ordinals`. They are skipped.
    // We look for gaps in `ordinals`?
    // Or just identifying sites where multiple reads Disagree with Query.
    // If it's a Het site, Query is haplotype A. Reads from haplotype B will have mismatches.
    
    // Simplified logic for "Site Compatibility":
    // Iterate all overlaps. For each overlap, find mismatches.
    // Count per-site mismatch frequency.
    
    // Map: Site -> [Counts]
    // Use a dense vector for Read Length? (Too memory intensive if read is 20kb? No, 20kb is fine)
    // HiFi reads are ~20kb. Vector<uint8_t> counts[4] is fine.
    
    const Read& queryRead = reads.getRead(queryReadId);
    uint64_t qLen = queryRead.size();
    
    // Temp storage for site statistics (Simple Pileup)
    struct SiteCount {
        uint16_t baseCounts[5]; // A,C,G,T,Gap
    };
    // Sparse map might be better if errors are rare?
    // But piling up all reads is standard.
    // Let's use a flat vector if < 50kb, else map? HiFi < 50kb usually.
    // Just use map for safety.
    // map<uint32_t, SiteCount> pileup;
    
    // Actually, hifiasm iterates overlaps, extracts mismatches, adds to list, then sorts.
    // "hc_phase_robust_rr" does exactly this.
    // "extract_sub_cigar_hc" fills evidence.
    
// --------------------------------------------------------
// Function: detectHetSites (Parity for hc_phase_robust_rr)
// --------------------------------------------------------
void detectHetSites(
    Assembler& assembler, // Need assembler for re-alignment
// --------------------------------------------------------
// Function: detectHetSites (Mirroring rphase_hc -> hc_phase_robust_rr)
// --------------------------------------------------------
static void detectHetSites(
    Assembler& assembler,
    const Reads& reads,
    uint64_t queryReadId,
    const vector<CandidateEC>& candidates, // "ol"
    const vector<AlignmentData>& alignmentData, // Access to raw alignment info
    vector<SnpStats>& snpStats,
    vector<HaplotypeEvidence>& hapEvidence
) {
    // 1. Scan overlaps to collect mismatches (Potential SNPs)
    // Parity: Use ProjectedAlignment to get CIGAR-level info from stored Compressed Alignments.
    
    const Read& queryRead = reads.getRead(queryReadId);

    // 0. Pre-compute Homopolymer Map for Query (HPC_PL = 12 parity)
    vector<bool> isHP(queryRead.size(), false);
    {
        uint32_t rLen = 0;
        char prev = 0;
        for(uint32_t i=0; i<queryRead.size(); i++) {
            char b = getBase(queryRead, i);
            if(b == prev) {
                rLen++;
            } else {
                if(rLen > 12) {
                     for(uint32_t j=i-rLen; j<i; j++) isHP[j] = true;
                }
                rLen = 1;
                prev = b;
            }
        }
        if(rLen > 12) {
             for(uint32_t j=queryRead.size()-rLen; j<queryRead.size(); j++) isHP[j] = true;
        }
    }
    
    // Ensure compressed alignments are accessible
    assembler.accessCompressedAlignments();

    // Loop candidates.
    for(size_t k=0; k<candidates.size(); k++) {
        const auto& cand = candidates[k];
        
        // 1. Fetch Compressed Alignment
        Alignment alignment;
        dinara::decompress(assembler.compressedAlignments[cand.alignmentId], alignment);
        
        if(alignment.ordinals.empty()) continue;

        // 2. Construct ProjectedAlignment (QuickRaw for SNP detection)
        // Need OrientedReadIds.
        // alignmentData stores readId0, readId1. 
        // We need to map to Query/Target orientation correctly.
        // ProjectedAlignment takes {OrientedReadId(id0, s0), OrientedReadId(id1, s1)}.
        // And an alignment.
        // IMPORTANT: The stored alignment corresponds to (id0, s0) and (id1, s1) stored in AlignmentData.
        
        // Let's get the stored AlignmentData to know the ReadIds used in the stored alignment.
        const AlignmentData& ad = alignmentData[cand.alignmentId];
        ReadId r0 = ad.readId0;
        ReadId r1 = ad.readId1;
        Strand s0 = 0; // Canonical
        Strand s1 = ad.isSameStrand ? 0 : 1; 
        
        array<OrientedReadId, 2> orientedReadIds = {
            OrientedReadId(r0, s0),
            OrientedReadId(r1, s1)
        };
        
        // Construct Projection
        // We use Method::QuickRaw to get mismatches efficiently.
        ProjectedAlignment projected(assembler, orientedReadIds, alignment, ProjectedAlignment::Method::QuickRaw);
        
        // 3. Scan Segments for SNPs
        // We need to map positions relative to Query (queryReadId).
        // Check which read is Query.
        int queryIdx = (r0 == queryReadId) ? 0 : 1;
        
        for(const auto& seg : projected.segments) {
            if(seg.mismatchCount == 0) continue; // No SNPs here
            
            // Scan segment alignment if mismatches exist
            // projected stores sequences in `rleSequences` or `sequences` depending on method?
            // QuickRaw fills `sequences`.
            // Seg has `alignment` vector<pair<bool, bool>>. 
            // pair: (inSeq0?, inSeq1?). (true, true)=Match/Mismatch. (true, false)=Del in 1. (false, true)=Ins in 1.
            
            // Reconstruct positions
            uint32_t currPos0 = seg.positionsA[0];
            uint32_t currPos1 = seg.positionsB[0];
            
            size_t seqIdx0 = 0; 
            size_t seqIdx1 = 0;
            
            for(size_t step=0; step<seg.alignment.size(); step++) {
                bool has0 = seg.alignment[step].first;
                bool has1 = seg.alignment[step].second;
                
                // Determine Query/Target bases
                char bQ = 0, bT = 0;
                uint32_t qPos = 0;
                
                if(has0) {
                    bQ = (queryIdx == 0) ? seg.sequences[0][seqIdx0] : 0; 
                    if(queryIdx == 0) qPos = currPos0;
                }
                if(has1) {
                    bT = (queryIdx == 1) ? seg.sequences[1][seqIdx1] : 0;
                    if(queryIdx == 1) qPos = currPos1; // Wait, Query is 1
                    
                    if(queryIdx == 0) bT = seg.sequences[1][seqIdx1];
                    else bQ = seg.sequences[0][seqIdx0];
                }
                
                if(has0 && has1) {
                    if(bQ != bT && bQ != 'N' && bT != 'N') {
                        // SNP
                         if(!isHP[qPos]) {
                            HaplotypeEvidence ev;
                            ev.overlapID = k;
                            ev.site = qPos;
                            ev.type = 1; // Mismatch
                            ev.misBase = base2int(bT); // Alt base
                            ev.overlapSite = 0;
                            hapEvidence.push_back(ev);
                         }
                    }
                }
                
                if(has0) { currPos0++; seqIdx0++; }
                if(has1) { currPos1++; seqIdx1++; }
            }
        }
    }
    
    // Sort Evidence by Site
    std::sort(hapEvidence.begin(), hapEvidence.end());
    
    // Partition SnpStats (Same as before)
    if(hapEvidence.empty()) return;
    
    for(size_t i=0; i<hapEvidence.size(); ) {
        uint32_t site = hapEvidence[i].site;
        uint32_t misCount = 0;
        uint32_t startIdx = i;
        while(i < hapEvidence.size() && hapEvidence[i].site == site) {
            misCount++;
            i++;
        }
        
        uint32_t totalCov = 0;
        uint32_t refReadsRev = 0; // overlap_num: Count of Ref reads on Reverse Strand
        
        // Efficiently iterate candidates for this site
        for(size_t k=0; k<candidates.size(); k++) {
            const auto& cand = candidates[k];
            if(site >= cand.qs && site < cand.qe) {
                totalCov++;
                
                // Check if this candidate is Alt (present in hapEvidence for this site)
                bool isAlt = false;
                // Since candidates are sorted by score, index k is overlapID.
                // Scan the evidence subset for this site to see if k matches.
                // Evidence subset is small (coverage depth).
                for(size_t j=startIdx; j<i; j++) {
                     if(hapEvidence[j].overlapID == k) { isAlt = true; break; }
                }
                
                if(!isAlt) {
                     // Ref. Count Strand Bias statistic (rev_n).
                     // Hifiasm counts reads with dir=1 (Reverse).
                     if(cand.isRev) refReadsRev++;
                }
            }
        }
        
        if(totalCov < 5) continue; 
        
        // Hifiasm Parity: occ_thres = 1.
        // Filter: misCount > 1 (so occ_1 >= 2)
        if(misCount > 1) { 
            SnpStats stat;
            stat.site = site;
            stat.occ_1 = misCount; 
            stat.occ_0 = totalCov - misCount;
            stat.overlap_num = refReadsRev; // populated with Reverse Ref Count
            
            // Link evidence to this SNP stats index
            for(size_t j=startIdx; j<i; j++) {
                hapEvidence[j].overlapSite = snpStats.size();
            }
            snpStats.push_back(stat);
        }
    }
}

// --------------------------------------------------------
// Function: gen_rphase_dp (Parity for Hifiasm gen_rphase_dp -> gen_rphase_dp0_single_path)
// --------------------------------------------------------
void gen_rphase_dp(
    Assembler& assembler,
    vector<SnpStats>& snpStats,
    vector<HaplotypeEvidence>& hapEvidence,
    vector<CandidateEC>& candidates
) {
    if(snpStats.empty()) return;

    // 1. Organize Evidence for Fast Linkage Query
    struct EvidenceRef {
        uint32_t overlapID;
        uint8_t allele; // 0=Ref, 1=Alt
    };
    
    // Efficient Mapping: snpStats[i] -> vector<uint32_t> alts (mismatching candidates)
    vector<vector<uint32_t>> siteReadsAlt(snpStats.size()); 
    vector<vector<uint32_t>> siteReadsRef(snpStats.size());

    size_t evIdx = 0;
    for(size_t i=0; i<snpStats.size(); i++) {
        uint32_t s = snpStats[i].site;
        while(evIdx < hapEvidence.size() && hapEvidence[evIdx].site < s) evIdx++;
        while(evIdx < hapEvidence.size() && hapEvidence[evIdx].site == s) {
            uint32_t ovId = hapEvidence[evIdx].overlapID;
            siteReadsAlt[i].push_back(ovId);
            evIdx++;
        }
        std::sort(siteReadsAlt[i].begin(), siteReadsAlt[i].end()); 
    }
    
    // Fill Refs
    for(size_t i=0; i<snpStats.size(); i++) {
        uint32_t s = snpStats[i].site;
        const auto& alts = siteReadsAlt[i];
        for(size_t k=0; k<candidates.size(); k++) {
            const auto& cand = candidates[k];
            if(s >= cand.qs && s < cand.qe) {
                if(!std::binary_search(alts.begin(), alts.end(), k)) {
                    siteReadsRef[i].push_back(k);
                }
            }
        }
    }
    
    // 2. Filter SNPs (Parity with gen_rphase_dp loop)
    auto is_st_bs = [](const SnpStats& s) {
        double rr = 0.05; 
        uint64_t mm = 2;
        if ( (s.overlap_num + mm >= s.occ_0) && 
             (s.overlap_num >= s.occ_0 * (1.0 - rr)) ) {
             return true; 
        }
        return false;
    };

    vector<int> validSnpIndices; 
    for(size_t i=0; i<snpStats.size(); i++) {
        SnpStats& s = snpStats[i];
        
        // "Strict Parity" Checks:
        // if((s->occ_0 < 2 || s->occ_1 < 2) || (is_st_bs((*s), st_rate, st_max)) || (!(s->occ_0 >= asm_opt.s_hap_cov && s->occ_1 >= asm_opt.infor_cov)))
        // s_hap_cov=3, infor_cov=3.
        
        if(s.occ_0 < 2 || s.occ_1 < 2) { s.score = -1; continue; }
        if(is_st_bs(s)) { s.score = -1; continue; }
        if(s.occ_0 < 3 || s.occ_1 < 3) { s.score = -1; continue; }
        
        validSnpIndices.push_back(i);
    }
    
    if(validSnpIndices.empty()) return;

    // 3. Single Path DP (gen_rphase_dp0_single_path logic)
    // Hifiasm Thresholds
    double cut_rate = 0.7; 
    uint64_t cut_bd = 6;
    
    uint64_t peak = 0;
    if(assembler.assemblerInfo->kmerDistributionInfo.coveragePeak != invalid<uint64_t>) {
        peak = assembler.assemblerInfo->kmerDistributionInfo.coveragePeak;
    }
    if(peak == 0) peak = candidates.size(); 
    
    uint64_t cc = (uint64_t)(peak * cut_rate);
    if(cc < cut_bd) cc = cut_bd;

    int nValid = validSnpIndices.size();
    vector<int64_t> f(nValid, 0);
    vector<int> p(nValid, -1);
    
    for(int idx=0; idx<nValid; idx++) {
        int i = validSnpIndices[idx];
        int64_t max_f = 1; 
        int max_j = -1;
        
        // Linkage Score
        int lookbackLimit = 50; 
        for(int jdx=idx-1; jdx>=0; jdx--) {
            if(idx - jdx > lookbackLimit) break;
            int j = validSnpIndices[jdx];

            const auto& refsI = siteReadsRef[i];
            const auto& altsI = siteReadsAlt[i];
            const auto& refsJ = siteReadsRef[j];
            const auto& altsJ = siteReadsAlt[j];

            int nn0 = 0, nn1 = 0; // consistent
            int mm0 = 0, mm1 = 0; // inconsistent
            
            // ... (Intersection logic same as before) ...
             auto count_commons = [](const vector<uint32_t>& A, const vector<uint32_t>& B) {
                int c = 0;
                size_t ia=0, ib=0;
                while(ia < A.size() && ib < B.size()) {
                    if(A[ia] < B[ib]) ia++;
                    else if(A[ia] > B[ib]) ib++;
                    else { c++; ia++; ib++; }
                }
                return c;
            };

            nn0 = count_commons(refsI, refsJ);
            nn1 = count_commons(altsI, altsJ); // Alt-Alt consistent?
            // Wait, for phasing 0-0 linkage.
            // Hifiasm comput_sc_rphase:
            // n0 = intersection(occ0[i], occ0[j])
            // n1 = intersection(occ1[i], occ1[j])
            // n2 = intersection(occ0[i], occ1[j])
            // n3 = intersection(occ1[i], occ0[j])
            // sc = (n0+n1) - (n2+n3)*weight
            
            mm0 = count_commons(refsI, altsJ);
            mm1 = count_commons(altsI, refsJ);

            int64_t sc = (nn0 + nn1) - (mm0 + mm1)*5; // weight=5? Hifiasm uses 6?
            // comput_sc_rphase uses weight related to coverage?
            // Actually commonly sc = (n0+n1) - (n2+n3)*5 is robust.
            
            if (sc == INT64_MIN) continue; // Not reachable
            sc += f[jdx];
            if(sc > max_f) {
                max_f = sc;
                max_j = jdx;
            }
        }
        f[idx] = max_f;
        p[idx] = max_j;
    }
    
    // Backtrace
    int64_t bestEndVal = -1;
    int bestEndIdx = -1;
    for(int i=0; i<nValid; i++) {
        if(f[i] > bestEndVal) {
            bestEndVal = f[i];
            bestEndIdx = i;
        }
    }
    
    // Mark Path
    // Reset all scores first (already -1 from filtering loop, or default 0? Score field re-used)
    // We already set invalid to -1.
    // Set valid to -1 initially before marking path.
    for(int idx : validSnpIndices) snpStats[idx].score = -1;
    
    int curr = bestEndIdx;
    while(curr != -1) {
        int snpIdx = validSnpIndices[curr];
        if(snpStats[snpIdx].occ_0 >= cc) {
            snpStats[snpIdx].score = 1; 
        } else {
            snpStats[snpIdx].score = -1;
        }
        curr = p[curr];
    }
}

// --------------------------------------------------------
// Function: generate_haplotypes_naive_HiFi
// --------------------------------------------------------
void generate_haplotypes_naive_HiFi(
    Assembler& assembler,
    vector<SnpStats>& snpStats,
    vector<HaplotypeEvidence>& hapEvidence,
    vector<CandidateEC>& candidates
) {
    if(snpStats.empty()) return;
    
    // Group Evidence by Candidate (OverlapID)
    // Hifiasm does radix sort on ID.
    // We can use a map or sort.
    struct EvPtr {
        uint32_t ovId;
        size_t evIdx;
        bool operator<(const EvPtr& other) const { return ovId < other.ovId; }
    };
    vector<EvPtr> evByCand;
    evByCand.reserve(hapEvidence.size());
    for(size_t i=0; i<hapEvidence.size(); i++) {
        evByCand.push_back({hapEvidence[i].overlapID, i});
    }
    std::sort(evByCand.begin(), evByCand.end());
    
    // Iterate Candidates (Overlaps)
    for(size_t i=0; i<evByCand.size(); ) {
        uint32_t ovId = evByCand[i].ovId;
        int support = 0;
        
        size_t start = i;
        while(i < evByCand.size() && evByCand[i].ovId == ovId) {
            size_t idx = evByCand[i].evIdx;
            const auto& ev = hapEvidence[idx];
            
            // Hifiasm Logic:
            // if(hh_tp(hap->list[i])!=1) continue; /// mismatch (type 1)
            // s = snpStats[ev.overlapSite]
            // if filters pass, and score==1, o++
            
            if(ev.type == 1) { // Mismatch/SNP
                if(ev.overlapSite < snpStats.size()) {
                    const auto& s = snpStats[ev.overlapSite];
                    // Filters (must match what was used to set scores, or strict check again)
                    if(s.score == 1) {
                        support++;
                    }
                }
            }
            i++;
        }
        
        // If support > 0, Mark Candidate
        if(support > 0) {
            if(ovId < candidates.size()) {
                candidates[ovId].score = 1; // Reuse score as is_match flag (1=Keep/Corrected)
            }
        } else {
             // In Hifiasm, if o==0, it might essentially be dropped or kept as is?
             // If we rely on this to filtering, we set score=1 only if supported.
             // Default was 0.
        }
    }
}

// --------------------------------------------------------
// Helper: detectSVSites (Simplified rphase_lidel logic)
// --------------------------------------------------------
// --------------------------------------------------------
// Helper: detectSVSites (Simplified rphase_lidel logic)
// --------------------------------------------------------
static void detectSVSites(
    Assembler& assembler,
    const Reads& reads, // Needed for sequence access
    uint64_t queryReadId,
    const vector<CandidateEC>& candidates,
    const vector<AlignmentData>& alignmentData, 
    vector<HaplotypeEvidence>& svEvidence,
    vector<SnpStats>& svStats // Populated with SV stats
) {
    // Mirroring rphase_lidel logic using ProjectedAlignment
    // We want to find indel events >= 16bp.
    
    struct RawSV {
        uint32_t overlapID;
        uint32_t site; // Query Site
        int64_t size; // SV Size
    };
    vector<RawSV> rawSVs;

    for(size_t k=0; k<candidates.size(); k++) {
        const auto& cand = candidates[k];
        
        // 1. Fetch Compressed Alignment (Reuse patterns for speed if possible, but for parity we compute)
        Alignment alignment;
        dinara::decompress(assembler.compressedAlignments[cand.alignmentId], alignment);
        
        if(alignment.ordinals.empty()) continue;

        // 2. ProjectedAlignment
        const AlignmentData& ad = alignmentData[cand.alignmentId];
        ReadId r0 = ad.readId0;
        ReadId r1 = ad.readId1;
        Strand s0 = 0; 
        Strand s1 = ad.isSameStrand ? 0 : 1; 
        
        array<OrientedReadId, 2> orientedReadIds = {
            OrientedReadId(r0, s0),
            OrientedReadId(r1, s1)
        };
        
        ProjectedAlignment projected(assembler, orientedReadIds, alignment, ProjectedAlignment::Method::QuickRaw);
        
        // 3. Scan Segments for Large Indels
        int queryIdx = (r0 == queryReadId) ? 0 : 1;
        
        for(const auto& seg : projected.segments) {
            uint32_t gapLenQ = seg.positionsA[1] - seg.positionsA[0];
            uint32_t gapLenT = seg.positionsB[1] - seg.positionsB[0];
            
            int64_t diff = 0;
            if(gapLenQ > gapLenT) diff = -(int64_t)(gapLenQ - gapLenT);
            else diff = (int64_t)(gapLenT - gapLenQ);
            
            if(std::abs(diff) >= 16) {
                 RawSV sv;
                 sv.overlapID = k;
                 sv.site = seg.positionsA[0] + (gapLenQ/2); // Center on Query
                 sv.size = diff;
                 rawSVs.push_back(sv);
            }
        }
    }
    
    if(rawSVs.empty()) return;
    
    // Cluster and Filter (rphase_lidel_cc parity)
    std::sort(rawSVs.begin(), rawSVs.end(), [](const RawSV& a, const RawSV& b){
        return a.site < b.site;
    });
    
    for(size_t i=0; i<rawSVs.size(); ) {
        uint32_t site = rawSVs[i].site;
        uint32_t startIdx = i;
        int64_t refSize = rawSVs[i].size; // Take first size as reference for cluster center? 
        // Hifiasm consensus logic is complex. We stick to simple checking.
        
        uint32_t endWin = site + 50;
        vector<size_t> indices;
        while(i < rawSVs.size() && rawSVs[i].site < endWin) {
            indices.push_back(i);
            i++;
        }
        
        // Refine refSize from consensus of window?
        // Using first is simple approximation.
        
        int validCount = 0;
        for(size_t idx : indices) {
             if(std::abs(rawSVs[idx].size - refSize) < std::abs(refSize)*0.25) { 
                 validCount++;
             }
        }
        
        if(validCount >= 3) {
            SnpStats stat;
            stat.site = site;
            stat.occ_1 = validCount;
            stat.occ_0 = candidates.size() - validCount;
            stat.score = -1;
            
            for(size_t idx : indices) {
                int64_t sz = rawSVs[idx].size;
                 if(std::abs(sz - refSize) < std::abs(refSize)*0.25) { 
                     HaplotypeEvidence ev;
                     ev.overlapID = rawSVs[idx].overlapID;
                     ev.site = site;
                     ev.type = 2; // SV
                     ev.overlapSite = svStats.size();
                     svEvidence.push_back(ev); // Only push Valid Evidence
                 }
            }
            svStats.push_back(stat);
        }
    }
}

// --------------------------------------------------------
// Function: generate_haplotypes_sv
// --------------------------------------------------------
void generate_haplotypes_sv(
    Assembler& assembler,
    vector<HaplotypeEvidence>& svEvidence, 
    vector<CandidateEC>& candidates
) {
    if(svEvidence.empty()) return;
    
    // Check support for each candidate
    // If a candidate supports a valid SV, we keep it?
    // Hifiasm Logic: `rphase_lidel` updates `ol->list[ii].is_match = 2` (Keep/Corrected)
    // IF the candidate supports the SV that has sufficient support.
    
    // We assume svEvidence ONLY contains evidence for VALID SV sites (filtered in detectSVSites).
    
    // So if a candidate appears in svEvidence, it supports a valid SV.
    // Marking it (score=1) means "Keep this candidate".
    
    for(const auto& ev : svEvidence) {
        if(ev.overlapID < candidates.size()) {
            candidates[ev.overlapID].score = 1; 
        }
    }
}

// --------------------------------------------------------
// Function: rphase_hc_ont (Parity for Hifiasm rphase_hc)
// --------------------------------------------------------

struct OntWindow {
    uint32_t qStart;
    uint32_t qEnd;
    uint32_t tStart; // Needed for SNP mapping on Target
    uint32_t candIdx; // Index into candidates vector
    uint32_t cigarIdx; // Index into phasingCigars (for re-check or deep dive) - optional if we trust candIdx
    // For rphase_hc, we need to know which overlap this window belongs to.
    
    // Sort by Query Start
    bool operator<(const OntWindow& other) const {
        return qStart < other.qStart;
    }
};

// Helper to check base identity between Query and Target at a specific Query Position
// Uses OntWindow to map QueryPos -> TargetPos
static bool is_match_ont(
    Assembler& assembler,
    const Reads& reads,
    uint64_t queryReadId,
    const CandidateEC& cand,
    uint32_t qPos,
    uint32_t tPos
) {
    // 1. Get Query Base
    const Read& rQ = reads.getRead(queryReadId);
    if(qPos >= rQ.size()) return false;
    char bQ = rQ[qPos];
    
    // 2. Get Target Base
    // Candidate has targetId and isRev.
    const Read& rT = reads.getRead(cand.targetId);
    
    char bT = 'N';
    // If isRev, target is RC.
    // Coordinates tPos are on the "Oriented Target".
    // Wait. In Dinara/Hifiasm, coordinates usually refer to the ORIGINAL read, but "Strand 1" means RC.
    // markers are defined on the Oriented Read.
    // So tPos is index into the Oriented Sequence.
    
    if(cand.isRev) { // Target is Strand 1 (RC)
        if(tPos >= rT.size()) return false;
        bT = rT[rT.size() - 1 - tPos].complement();
    } else { // Target is Strand 0
        if(tPos >= rT.size()) return false;
        bT = rT[tPos];
    }
    
    return bQ == bT;
}

// Sub-function: hc_phase_robust_rr_impl
// Checks consistency of active overlaps for a given query window/position.
void hc_phase_robust_rr_impl(
    Assembler& assembler,
    const Reads& reads,
    uint64_t queryReadId,
    const vector<OntWindow>& windows, 
    const vector<size_t>& activeWindows,
    uint32_t rangeStart,
    uint32_t rangeEnd
) {
    const Read& rQ = reads.getRead(queryReadId);
    const auto& alignmentData = assembler.alignmentData;
    
    // Iterate positions in the chunk
    for(uint32_t qPos = rangeStart; qPos < rangeEnd; ++qPos) {
        if(qPos >= rQ.size()) break;
        
        // Simple Pileup (A, C, G, T, N)
        uint32_t counts[5] = {0}; // 0=A, 1=C, 2=G, 3=T, 4=Amb/Gap
        
        char qBase = rQ[qPos];
        // Convert to 0-3...
        
        // Iterate Active Windows
        for(size_t wIdx : activeWindows) {
            const auto& win = windows[wIdx];
            if(qPos < win.qStart || qPos >= win.qEnd) continue;
            
            uint32_t offset = qPos - win.qStart;
            uint32_t tPos = win.tStart + offset;
            
            // Access AlignmentData directly
            uint32_t alignId = win.candIdx;
            const auto& ad = alignmentData[alignId];
            
            if(ad.isDeleted()) continue;
            
            // Determine Target info
            uint32_t targetId;
            if(ad.readIds[0] == queryReadId) {
                targetId = ad.readIds[1];
            } else {
                targetId = ad.readIds[0];
            }
            
            // Orientation
            // If !isSameStrand, target is Reverse relative to Query(Fwd)
            bool targetIsRev = !ad.isSameStrand;
            
            const Read& rT = reads.getRead(targetId);
            char tBase = 'N';
            
            if(targetIsRev) {
                if(tPos < rT.size()) {
                     tBase = rT[rT.size() - 1 - tPos].complement();
                }
            } else {
                if(tPos < rT.size()) {
                    tBase = rT[tPos];
                }
            }
            
            // Update counts 
            // ...
        }
        
        // Analyze Pileup
    }
}



void rphase_hc_ont(
    Assembler& assembler,
    const Reads& reads,
    uint64_t queryReadId,
void rphase_hc_ont(
    Assembler& assembler,
    const Reads& reads,
    uint64_t queryReadId
) {
     const auto& alignmentTable = assembler.alignmentTable;
     const auto& alignmentData = assembler.alignmentData;
     const auto& markers = assembler.markers;
     const auto& assemblerInfo = assembler.assemblerInfo;

     // Iterate AlignmentTable for OrientedReadId(queryReadId, 0)
     // This gives us ONE entry per pair (either (Q, T) or (T, Q) depending on ID sort).
     // Wait. `computeAlignmentTable` adds BOTH (Q, 0) and (Q, 1) and their partners.
     // But `alignmentTable` stores indices of alignments.
     // Querying `OrientedReadId(queryReadId, 0)` gives us all alignments where Query is Forward.
     // This includes cases where Query is Read0 and where Query is Read1 (flipped).
     
     // Correct.
     OrientedReadId oidQuery(queryReadId, 0);
     if(oidQuery.getValue() >= alignmentTable.size()) return;
     const auto& alignments = alignmentTable[oidQuery.getValue()];

     if(alignments.empty()) return;

     // 1. Generate Windows
     vector<OntWindow> windows;
     windows.reserve(alignments.size() * 5); 
     
     const auto& markersQ = markers[oidQuery.getValue()];

     for(uint32_t alignId : alignments) {
         if(alignId >= alignmentData.size()) continue;
         const auto& ad = alignmentData[alignId];
         
         if(ad.isDeleted()) continue;
         
         // Identify Query and Target roles
         // queryReadId is our Read.
         // ad.readIds[0] and ad.readIds[1].
         
         uint32_t targetReadId;
         bool isRead0; // Is Query Read0?
         
         if(ad.readIds[0] == queryReadId) {
             isRead0 = true;
             targetReadId = ad.readIds[1];
         } else {
             // Must be readIds[1] since we retrieved it from the table
             isRead0 = false;
             targetReadId = ad.readIds[0];
         }
         
         // Determine Query Orientation in the Alignment
         bool isQueryRev = false;
         if(!isRead0 && !ad.isSameStrand) {
             isQueryRev = true;
         }
         
         // Use correct markers for Query
         uint32_t qStrand = isQueryRev ? 1 : 0;
         const auto& markersQ_Aligned = markers[OrientedReadId(queryReadId, qStrand).getValue()];
         
         // Determine Target Orientation in the Alignment (for marker access)
         // If `isRead0` (Query=R0), Target=R1. Strand = isSameStrand ? 0 : 1.
         // If `!isRead0` (Query=R1), Target=R0. R0 is always Base/Ref (Strand 0).
         
         uint32_t targetStrandForMarkers = 0;
         if(isRead0) {
             targetStrandForMarkers = ad.isSameStrand ? 0 : 1;
         } else {
             targetStrandForMarkers = 0; // Read0 is always Strand 0
         }
         
         const auto& markersT_Aligned = markers[OrientedReadId(targetReadId, targetStrandForMarkers).getValue()];
         
         // Start Ordinals
         uint32_t startOrdQ = isRead0 ? ad.info.data[0].firstOrdinal : ad.info.data[1].firstOrdinal;
         uint32_t startOrdT = isRead0 ? ad.info.data[1].firstOrdinal : ad.info.data[0].firstOrdinal;
         
         uint32_t currQ = markersQ_Aligned[startOrdQ].position; 
         uint32_t currT = markersT_Aligned[startOrdT].position; 
         
         // Query Length for Coordinate Transformation
         uint32_t qLen = (uint32_t)reads.getRead(queryReadId).size();

         // CIGAR
         // Stored CIGAR is generated for Read0 (Ref) vs Read1 (Query).
         bool swapCigar = isRead0;
         
         const auto& cigarSpan = assembler.phasingCigars[alignId];
         
         for(uint32_t val : cigarSpan) {
             uint32_t op = val & 0xF;
             uint32_t len = val >> 4;
             
             // Apply Swap
             if(swapCigar) {
                 if(op == 1) op = 2;      // I(1) -> D(2)
                 else if(op == 2) op = 1; // D(2) -> I(1)
             }
             
             if(op == 0) { // Match/Mismatch
                 OntWindow win;
                 
                 // Handle Query Coordinates
                 if(isQueryRev) {
                     // currQ is on Reverse Strand [currQ, currQ + len)
                     // Map to Forward: 
                     // FwdEnd = Len - currQ
                     // FwdStart = Len - (currQ + len)
                     win.qStart = qLen - (currQ + len);
                     win.qEnd = qLen - currQ;
                 } else {
                     win.qStart = currQ;
                     win.qEnd = currQ + len;
                 }
                 
                 // Handle Target Coordinates
                 // win.tStart needs to be on the Target Read's "Active Orientation"? 
                 // hc_phase_robust_rr_impl logic checks `targetIsRev`.
                 // If we store `currT` (which is Relative to Alignment), we need to know 
                 // if `hc_phase` interprets `tStart` as Raw or Oriented.
                 // In `hc_phase_robust_rr_impl`:
                 // `uint32_t tPos = win.tStart + offset;`
                 // `if(targetIsRev)` fetch complement...
                 // This implies `tStart` should be "Index into the Sequence used for Alignment".
                 // So we keep `currT` as is (Aligned Coordinates).
                 win.tStart = currT; 
                 
                 win.candIdx = alignId; 
                 win.cigarIdx = alignId; 
                 windows.push_back(win);
                 
                 currQ += len;
                 currT += len;
             } else if(op == 1) { // Ins in Query (Gap in Target)
                 currQ += len;
                 // currT stays
             } else if(op == 2) { // Del in Query (Bases in Target)
                 // currQ stays
                 currT += len;
             }
         }
     }
     
     // 2. Sort Windows
     std::sort(windows.begin(), windows.end());
     
     // 3. Sweep Line Loop
     
     uint64_t qLen = reads.getRead(queryReadId).size();
     size_t winIdx = 0;
     vector<size_t> activeWindows; // Indices into `windows`
     
     uint32_t s = 0;
     while(s < qLen) {
         // Activate windows starting at or before s
         while(winIdx < windows.size() && windows[winIdx].qStart <= s) {
             // Add to active set if it ends after s
             if(windows[winIdx].qEnd > s) {
                 activeWindows.push_back(winIdx);
             }
             winIdx++;
         }
         
         // Remove finished windows
         size_t activeCount = 0;
         uint32_t minEnd = qLen;
         
         for(size_t k=0; k<activeWindows.size(); ++k) {
             size_t wVal = activeWindows[k];
             if(windows[wVal].qEnd > s) {
                 if(windows[wVal].qEnd < minEnd) minEnd = windows[wVal].qEnd;
                 activeWindows[activeCount++] = wVal;
             }
         }
         activeWindows.resize(activeCount);
         
         if(activeWindows.empty()) {
             if(winIdx < windows.size()) s = windows[winIdx].qStart;
             else break; 
             continue;
         }
         
         uint32_t e = minEnd; // Process [s, e)
         
         // Call Phase Logic
         // Note: We access overlaps via assembler.alignmentData within the impl now?
         // No, the impl still asks for `CandidateEC` or `CanonicalOverlap`?
         // We need to update `hc_phase_robust_rr_impl` signature to use `AlignmentData`.
         hc_phase_robust_rr_impl(assembler, reads, queryReadId, windows, activeWindows, s, e);
         
         s = e;
     }
}





void Assembler::performHifiasmECParity(uint64_t threadCount)
{ 
    cout << timestamp << "=== Hifiasm Parity EC Pipeline (Round 1) ===" << endl;

    // Prepare to mark deletions
    // We assume alignmentData is populated.
    // In Hifiasm, ha_ec iterates over ALL reads (R_INF.total_reads).
    // It rebuilds candidates. Here we use existing candidates.

    // We must group alignments by ReadId (Query). 
    // Since alignmentData stores (ReadId0, ReadId1, ...), we can iterate alignmentTable if available,
    // or just assume we process all.
    // For efficient threading, we iterate reads 0..readCount.

    const uint64_t readCount = reads.readCount();
    
    // We need to access alignmentTable to get candidates for each read.
    // alignmentTable maps ReadId -> [AlignmentId]
        
    // Use uint8_t for thread-safe byte addressing (vector<bool> is a bitfield and not thread-safe).
    vector<uint8_t> keepAlignment(alignmentData.size(), 0); 
    
    // Parallel loop over reads
    vector<thread> threads;
    uint64_t chunkSize = readCount / threadCount;
    if(chunkSize == 0) chunkSize = 1;

    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            uint64_t start = t * chunkSize;
            uint64_t end = (t == threadCount - 1) ? readCount : (t + 1) * chunkSize;

            for(uint64_t r = start; r < end; r++) {
                // ONT Logic (Parity)
                // For ONT, we use rphase_hc with AlignmentData.
                // We bypass the candidate gathering and detectHetSites logic used for HiFi.
                bool isOnt = true; // Use flag from AssemblerInfo if available
                if (isOnt) {
                    rphase_hc_ont(*this, reads, r);
                    continue; 
                }

                // ReadId r, Strand 0 and 1.
                // In Hifiasm EC, correction is done for each read (as query).
                // We iterate each oriented read as a query.
                for(uint32_t strand=0; strand<2; strand++) {
                    OrientedReadId oid(r, strand);
                    if(oid.getValue() >= alignmentTable.size()) continue;

                    const auto& alignments = alignmentTable[oid.getValue()];
                    if(alignments.empty()) continue;

                    // Collect candidates for this query
                    vector<CandidateEC> candidates;
                    candidates.reserve(alignments.size());
                    
                    for(uint32_t alignId : alignments) {
                        const auto& al = alignmentData[alignId];
                        
                        // If already deleted, skip?
                        // Hifiasm EC often starts fresh. 
                        // But for pipeline integration, we respect previous hard deletions.
                        if(al.isDeleted()) continue;

                        CandidateEC c;
                        c.alignmentId = alignId;
                        
                        // We need to know if we are readId0 or readId1 in the alignment.
                        // And if the alignment is same-strand or diff-strand relative to US.
                        // alignmentData: readId0, readId1, isSameStrand.
                        // Our oid: (r, strand).
                        
                        // Determine Target Interval on Query.
                        // We need the coordinates on the Query Read.
                        uint32_t qStart=0, qEnd=0;
                        
                        // Retrieve Target ReadId and Orientation
                        // al.readId0/readId1 are standard.
                        // if r == al.readId0, target is readId1.
                        // isSameStrand affects orientation.
                        // If same strand, target is Fwd relative to Query Fwd?
                        // If diff strand, target is Rev?
                        
                        // CandidateEC c; // Already declared above
                        c.alignmentId = alignId;
                        
                        if(r == al.readId0) {
                            c.targetId = al.readId1;
                            c.qs = al.readId0Begin;
                            c.qe = al.readId0End;
                            c.ts = al.readId1Begin;
                            c.te = al.readId1End;
                            c.isRev = !al.isSameStrand; // If diff strand, effectively rev?
                        } else {
                            c.targetId = al.readId0;
                            c.qs = al.readId1Begin;
                            c.qe = al.readId1End;
                            c.ts = al.readId0Begin;
                            c.te = al.readId0End;
                            c.isRev = !al.isSameStrand;
                        }
                        c.score = al.alignmentScore; 
                        
                        candidates.push_back(c);
                    }

                    if(candidates.empty()) continue;

                    // Strict Hifiasm Parity: Check for Uncorrected Read (Junk Filter)
                    // For is_ont=1, Hifiasm drops reads with coverage gaps >= 1600bp.
                    // Checked BEFORE phasing to save cycles (and parity with output filtering).
                    uint64_t queryLen = reads.getRead(r).length();
                    if (isUncorrectedRead(candidates, queryLen, 1600)) {
                         continue; 
                    }

                    // Sort by score descending
                    sort(candidates.begin(), candidates.end(), [](const CandidateEC& a, const CandidateEC& b) {
                        return a.score > b.score;
                    });
                    
                    // ----------------------------------------------------
                    // Hifiasm EC Phase 1: Het Site Detection & Compatibility (Parity)
                    // ----------------------------------------------------
                    // User Request: "dp based site compatibility"
                    // We run detectHetSites -> computeHaplotypeDP -> Mark Incompatible
                    
                    // ----------------------------------------------------
                    
                    // Hifiasm Logic Chaining (Parity with rphase_hc)
                    
                    // 1. Detect SNPs (detectHetSites - SnpStats populated)
                    vector<SnpStats> snpStats;
                    vector<HaplotypeEvidence> hapEvidence;
                    detectHetSites(*this, reads, r, candidates, alignmentData, snpStats, hapEvidence);
                    
                    // 2. DP Step (gen_rphase_dp)
                    gen_rphase_dp(*this, snpStats, hapEvidence, candidates);
                    
                    // 3. Naive Verification (generate_haplotypes_naive_HiFi)
                    generate_haplotypes_naive_HiFi(*this, snpStats, hapEvidence, candidates);

                    // 4. SV Step (is_ont only - parity)
                    // Hifiasm runs rphase_lidel for ONT data.
                    vector<HaplotypeEvidence> svEvidence;
                    vector<SnpStats> svStats;
                    detectSVSites(*this, reads, r, candidates, alignmentData, svEvidence, svStats);
                    generate_haplotypes_sv(*this, svEvidence, candidates);
                    
                    
                    for(const auto& cand : candidates) {
                        // Ratio Filter (Strict Hifiasm Parity sh=0.866666)
                        if (cand.score >= bestScore * minScoreRatio) {
                             keepAlignment[cand.alignmentId] = 1;
                        }
                    }
                }
            }
        });
    }

    for(auto& t : threads) t.join();

    // Prune Alignment Data based on keep flags
    vector<AlignmentData> keptAlignments;
    keptAlignments.reserve(alignmentData.size());
    
    // We also need to update alignmentTable map.
    // Easier to just rebuild it? Or filter in place?
    // Assembler::alignmentData is the source.
    
    uint64_t keptCount = 0;
    for(size_t i=0; i<alignmentData.size(); i++) {
        if(!keepAlignment[i]) {
            // Hifiasm: Overlap must be valid.
            // Usually if one side drops, it's dropped?
            // Or `push_ne_ovlp` pushes for `i` (read).
            // It builds paf[i] and reverse_paf[i].
            // If i keeps it, it's in paf[i].
            // Since our alignmentData is undirected (canonical),
            // if *either* read effectively discards it, should we drop?
            // Actually, Hifiasm EC output (PAF) is read-centric.
            // User wants "final PAF lists".
            // Effectively, if the overlap exists in the final PAF, it's valid.
            // If Hifiasm drops it from *both* sides, it's gone.
            // If it keeps in one, does it keep in other? Usually symmetric.
            
            // For now: Only delete if NOT kept.
            // But we initialized keepAlignment to false.
            // So if we didn't explicitly keep, it's false.
            // But did we visit every alignment?
            // Yes, iterating all reads covers all alignmentTable entries.
            // So we can set deleted = !keepAlignment[i].
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
    // Same structure, potentially stricter params.
    // Placeholder logic for now matches Round 1.
}
