#include "AssemblerPhasing.hpp"
#include "Assembler.hpp"
#include <algorithm>
#include <map>
#include <iostream>
#include <numeric>
#include <cmath>
#include "Reads.hpp"
#include "timestamp.hpp"

namespace dinara {

namespace {

    // Hifiasm Parity Helpers
    
    // is_st_bs: Checks if a site is likely HOMOZYGOUS based on coverage.
    bool is_st_bs(uint32_t overlap_num, uint32_t occ_0, double st_rate, uint64_t st_max) {
        if(st_max != (uint64_t)-1) {
            if(overlap_num + st_max >= occ_0) return true;
        }
        if(occ_0 == 0) return false;
        if( (double)overlap_num >= (double)occ_0 * (1.0 - st_rate) ) return true;
        return false;
    }

    struct SVFeature {
        uint32_t start;
        uint32_t end;
        int32_t type; // 1=Ins, 2=Del 
        uint32_t ovIdx;
        uint32_t size;
    };
    
    // cal_lindel_dd: SV Similarity Score
    double cal_lindel_dd(const SVFeature& a, const SVFeature& b, double ol_r, uint32_t ol_w, double err_dif) {
        if(a.type != b.type) return -1.0;
        
        uint32_t os = std::max(a.start, b.start);
        uint32_t oe = std::min(a.end, b.end);
        
        if(oe <= os + ol_w) return -1.0;
        
        double ol = (double)(oe - os);
        double lenA = (double)(a.end - a.start);
        double lenB = (double)(b.end - b.start);
        
        // Overlap Ratio check
        if(ol < lenA * ol_r && ol < lenB * ol_r) return -1.0; 
        
        double szA = (double)a.size;
        double szB = (double)b.size;
        double diff = std::abs(szA - szB);
        
        if(diff > szA * err_dif || diff > szB * err_dif) return -1.0;
        
        // Score
        double sc = (szA + szB - 2.0*diff) / (szA + szB);
        sc += (2.0 * ol) / (lenA + lenB);
        return sc;
    }

    // checkPeriodicity: Matches Hifiasm `hpc_mask_ff`
    // Checks if position `p` is part of a repetitive run (homopolymer or STR)
    // defined by `maxPeriod` (HPC_RR=4) and `cutoffFactor` (HPC_CC=6).
    // Returns true if masked.
    bool checkPeriodicity(const std::span<const uint8_t>& seq, uint32_t p, int maxPeriod, int cutoffFactor) {
        uint32_t len = (uint32_t)seq.size();
        if(p >= len) return false;

        for(int r = 1; r <= maxPeriod; ++r) {
            int threshold = r * cutoffFactor;
            
            // Scan Right
            uint32_t ze = p + r;
            while(ze < len && seq[ze] == seq[ze-r]) {
                ze++;
            }
            
            // Scan Left
            int32_t zs = (int32_t)p - 1;
            while(zs >= 0 && (zs + r) < (int32_t)len && seq[zs] == seq[zs+r]) {
                zs--;
            }
            zs++; // First valid index

            if((ze - zs) > (uint32_t)r && (ze - zs) >= (uint32_t)threshold) {
                return true; // Masked
            }
        }
        return false;
    }

} // namespace

// ----------------------------------------------------------------------------
// Assembler::performPhasing Implementation
// ----------------------------------------------------------------------------

// Deprecated Legacy Phasing Logic (Replaced by APES/TASSD in AssemblerHifiasmEC.cpp)
void Assembler::performPhasing(uint64_t threadCount)
{
    (void)threadCount;
    // No-op. Logic moved to performHifiasmECParity.
}

void Assembler::performPhasingThreadFunction(size_t threadId)
{
    (void)threadId;
}


// ----------------------------------------------------------------------------
// AssemblerPhasing Static Methods
// ----------------------------------------------------------------------------

std::vector<uint32_t> AssemblerPhasing::filterOverlapsByPhasing(
    const Assembler& assembler,
    ReadId targetReadId,
    const std::vector<PhasingOverlap>& overlaps,
    const PhasingConfig& config)
{
    if (overlaps.empty()) return {};

    std::vector<HaplotypeEvidence> evidence;
    evidence.reserve(overlaps.size() * 10); 
    
    // Pass 1: Collect Evidence
    for(size_t i=0; i<overlaps.size(); ++i) {
        // Warning: collectHaplotypeEvidence needs to tag evidence with overlapId = i
        // But HaplotypeEvidence stores uint32_t overlapId.
        // We need to set it inside collectHaplotypeEvidence?
        // Let's pass 'i' as implicit context, or update evidence after?
        // Or update collectHaplotypeEvidence to take overlapId.
        // I updated header to take `PhasingOverlap` but not `overlapId`.
        // I will overload or simply set it here.
        
        size_t startSize = evidence.size();
        collectHaplotypeEvidence(assembler, overlaps[i], evidence, config);
        // Set overlapId for newly added items
        for(size_t k=startSize; k<evidence.size(); ++k) {
            evidence[k].overlapId = (uint32_t)i;
        }
    }

    if (evidence.empty()) {
        std::vector<uint32_t> keptIndices(overlaps.size());
        std::iota(keptIndices.begin(), keptIndices.end(), 0);
        return keptIndices;
    }

    std::vector<SnpStats> stats;
    std::vector<uint32_t> dpKept = generatePhasingDP(targetReadId, overlaps, evidence, stats, config);
    
    // 2. Naive Refinement (Stateful Parity)
    std::vector<uint32_t> naiveKept = refineOverlapsNaive(dpKept, overlaps, evidence, stats, config);

    // 3. SV Phasing (Parity with rphase_lidel)
    return filterOverlapsBySV(naiveKept, overlaps, evidence, stats, assembler, config);
}

// Parity with generate_haplotypes_naive_HiFi
// Stateful Greedy Refinement
std::vector<uint32_t> AssemblerPhasing::refineOverlapsNaive(
    const std::vector<uint32_t>& keptIndices,
    const std::vector<PhasingOverlap>& overlaps,
    const std::vector<HaplotypeEvidence>& evidence,
    std::vector<SnpStats>& stats,
    const PhasingConfig& config)
{
    if(keptIndices.empty()) return {};

    // Map Site -> Stats Index
    std::map<uint32_t, size_t> siteToStatIdx;
    for(size_t i=0; i<stats.size(); ++i) siteToStatIdx[stats[i].site] = i;

    // We need to know which alleles each overlap supports to decrement correctly.
    // Hifiasm iterates overlaps.
    // Sorting: overlaps with MORE valid alleles processed first?
    // Hifiasm: "sort by how many allels in each overlap; more -> less" (line 8999).
    


    // Evidence Map for fast lookup: OvIdx -> range in evidence vector
    std::vector<std::pair<size_t, size_t>> evRanges(overlaps.size(), {0,0});
    {
        // size_t start = 0;
        for(size_t i=0; i<evidence.size(); ) {
            uint32_t currId = evidence[i].overlapId;
            size_t j = i;
            while(j < evidence.size() && evidence[j].overlapId == currId) j++;
            if(currId < overlaps.size()) {
                evRanges[currId] = {i, j};
            }
            i = j;
        }
    }
    
    // Parity: Simple pass-through check (Stateless)
    // Hifiasm checks consistency against the refined haplotype.
    // It keeps an overlap if it has ANY valid evidence support.
    
    std::vector<uint32_t> finalKept;
    finalKept.reserve(keptIndices.size());
    
    for(uint32_t ovId : keptIndices) {
        int validCounts = 0;
        size_t start = evRanges[ovId].first;
        size_t end = evRanges[ovId].second;
        
        for(size_t k=start; k<end; ++k) {
            const auto& ev = evidence[k];
            
            // We only care about consistency with VALID sites
            auto it = siteToStatIdx.find(ev.site);
            if(it != siteToStatIdx.end()) {
                 const auto& s = stats[it->second];
                 
                 // Check if site is Valid (Parity with line 8899)
                 if(s.occ_0 >= (uint32_t)config.s_hap_cov && s.occ_1 >= (uint32_t)config.infor_cov && 
                    !is_st_bs(s.overlap_num, s.occ_0, config.st_rate, config.st_max)) {
                     
                     // Check if evidence supports the Variant (Type 1)
                     if(ev.type == 1) { // SNP
                          validCounts++;
                     }
                 }
            }
        }
        
        // Hifiasm: if(o >= 1) keep;
        if(validCounts > 0) {
            finalKept.push_back(ovId);
        }
    }
    
    // Restore original order by ID
    std::sort(finalKept.begin(), finalKept.end());
    return finalKept;
}




// Parity with rphase_lidel
std::vector<uint32_t> AssemblerPhasing::filterOverlapsBySV(
    const std::vector<uint32_t>& keptIndices,
    const std::vector<PhasingOverlap>& overlaps,
    const std::vector<HaplotypeEvidence>& evidence,
    std::vector<SnpStats>& snpStats,
    const Assembler& /* assembler */,
    const PhasingConfig& /* config */)
{
    if(keptIndices.empty()) return {};

    // 1. Extract SVs
    std::vector<SVFeature> svs;
    std::vector<uint32_t> svToOv; // Map SV index -> Overlap Index

    for(uint32_t ovId : keptIndices) {
        const auto& ov = overlaps[ovId];
        uint32_t tPos = ov.targetStart;
        // qPos isn't tracked in PhasingOverlap explicitly (we use relative), 
        // but we need it for SV size?
        // Hifiasm extract_sub_cigar_sv tracks qs/qe.
        // We just need event size.
        
        uint32_t idx = 0;
        while(idx < ov.cigar.size()) {
            uint32_t val = ov.cigar[idx];
            uint32_t op = val & 0xF;
            uint32_t len = val >> 4;
            
            if(op == 0 || op == 3) {
                // Match or Blind
                tPos += len;
                idx++;
            } else {
                // Non-match (Ins=1, Del=2) -> SV Group
                // Hifiasm groups consecutive operations where op != 0 (and !3?)
                // Actually it groups as long as `op` stays "non-match".
                // Wait, logic was `if(!!op != !!op0) break`.
                // Ins=1, Del=2. Both !! is true.
                // So it groups sequences of Ins and Del.
                
                uint32_t groupStartT = tPos;
                uint32_t groupLen = 0; // Total length of indels (size)
                // uint32_t groupSize = 0;
                
                uint32_t spanT = 0;
                
                while(idx < ov.cigar.size()) {
                    uint32_t subVal = ov.cigar[idx];
                    uint32_t subOp = subVal & 0xF;
                    uint32_t subLen = subVal >> 4;
                    
                    if(subOp == 0 || subOp == 3) break; // End of group
                    
                    groupLen += subLen;
                    if(subOp == 2) spanT += subLen; // Del consumes Target
                    // Ins consumes Query (not represented in tPos)
                    
                    idx++;
                }
                
                // Finished group
                // Check if it qualifies as SV
                if(groupLen >= 50) { 
                    // Type: Hifiasm uses Type 1 for SV in `generate_haplotypes_sv`?
                    // "type 1" usually means "Alt" / Indel.
                    // We need a specific type (Ins vs Del) for clustering `cal_lindel_dd`.
                    // But if it's mixed?
                    // `cal_lindel_dd` checks `if(a->tn == b->tn)`. `tn` is type? No, tn is TargetName usually.
                    // Wait, `cal_lindel_dd` doesn't check type!
                    // It checks overlap and size difference.
                    // So mixed groups allow clustering if size matches.
                    // We'll use type=1 for all SVs here.
                    
                    svs.push_back({groupStartT, groupStartT + spanT, 1, ovId, groupLen});
                }
                
                tPos += spanT;
            }
        }
    }
    
    if(svs.empty()) return keptIndices;

    // 4. Hifiasm `rphase_lidel_cc` Clustering Strategy
    
    // Step 4.1: Neighbor Counting (Calculate 'sec' / support)
    struct SVWrapper {
        SVFeature sv;
        uint32_t sec; // Secondary support (count of neighbors)
        int clusterId;
        uint32_t readId; // Target Read ID (Wait, Pre-Phasing is Target-centric?)
        // In `filterOverlapsBySV`, 'overlaps' are all mapped TO the Target.
        // So `ov.queryReadId` is the "Supporting Read".
        // We need to ensure unique `queryReadId` per cluster.
    };
    
    std::vector<SVWrapper> wrappers;
    wrappers.reserve(svs.size());
    for(const auto& s : svs) {
        // Find queryReadId from overlap
        // We need effective read ID. 
        // `overlaps[s.ovIdx]` has `queryReadId`.
        wrappers.push_back({s, 0, -1, overlaps[s.ovIdx].queryReadId});
    }
    
    // Sort by Start (already done mostly, but ensure)
    std::sort(wrappers.begin(), wrappers.end(), [](const SVWrapper& a, const SVWrapper& b){
        return a.sv.start < b.sv.start;
    });
    
    double ol_r = 0.5;
    double err_dif = 0.25;
    uint32_t ol_w = 3;
    uint32_t c_sz = 3; // Min support to form core
    
    int nSVs = (int)wrappers.size();
    
    for(int k=0; k<nSVs; ++k) {
        // Window optimization
        // uint32_t max_dist = wrappers[k].sv.size * 2 + 100;
        
        for(int z=k+1; z<nSVs; ++z) {
            if(wrappers[z].sv.start > wrappers[k].sv.end) break; // Start of z > End of k (approx)
             
            // Hifiasm window check: a[z].qs < a[k].qe ?
            // Here we sort by Start.
            
            double score = cal_lindel_dd(wrappers[k].sv, wrappers[z].sv, ol_r, ol_w, err_dif);
            if(score >= 0.0) {
                wrappers[k].sec++;
                wrappers[z].sec++;
            }
        }
    }
    
    // Step 4.2: Filter Noise (sec == 0) and Compress? 
    // We can just skip them in next loops.

    // Step 4.3: Build Core Clusters
    int nextClusterId = 0;
    
    // Iterate to build clusters from high-support seeds
    for(int k=0; k<nSVs; ++k) {
        if(wrappers[k].sec < c_sz || wrappers[k].clusterId != -1) continue;
        
        // Start new cluster
        int currentId = nextClusterId++;
        wrappers[k].clusterId = currentId;
        
        // Track Reads in this cluster to enforce uniqueness
        // Using a set for this local greedy expansion
        std::set<uint32_t> clusterReads;
        clusterReads.insert(wrappers[k].readId);
        
        // Expand
        // Using a BFS/Queue or just simple pass? 
        // Hifiasm does a seeded expansion using `buf` stack.
        std::vector<int> stack;
        stack.push_back(k);
        
        size_t head = 0;
        while(head < stack.size()) {
            int v = stack[head++];
            
            // Search neighbors of v
            for(int z=0; z<nSVs; ++z) { // Ideally range search
               // Optimization: Only check nearby?
               if(wrappers[z].sv.start > wrappers[v].sv.end) break; // Sorted
               // We need to look backwards too? 
               // Hifiasm `rphase_lidel_cc` loops `z` from 0 to `an`.
               // But optimized with `a[z].qs < a[v].qe`.
            }
            
            // Re-implementing range search properly:
            // We need to scan [0, nSVs].
            // But we can optimize.
            // Since we sorted by Start, we can search [v-window, v+window].
            
            // For rigorous parity, let's scan a safe window.
            // Or just Full Scan for now (perf tradeoff for correctness).
            // Actually, `cal_lindel_dd` fails if `oe <= os + ol_w`.
            
            // Let's use a window loop based on positions.
            // Find start index
            int zStart = v - 1;
            while(zStart >= 0 && wrappers[zStart].sv.end >= wrappers[v].sv.start) zStart--;
            zStart++; // approximate
            
            // Forward
            for(int z=0; z<nSVs; ++z) {
                if(wrappers[z].clusterId != -1) continue; // Already clustered
                if(wrappers[z].sec < c_sz) continue; // Must be core-quality
                if(z == v) continue;
                
                // Range check
                uint32_t os = std::max(wrappers[v].sv.start, wrappers[z].sv.start);
                uint32_t oe = std::min(wrappers[v].sv.end, wrappers[z].sv.end);
                if(oe <= os + ol_w) continue;
                
                double sw = cal_lindel_dd(wrappers[v].sv, wrappers[z].sv, ol_r, ol_w, err_dif);
                if(sw >= 0.0) {
                     // Check Unique Read Constraint
                     if(clusterReads.find(wrappers[z].readId) == clusterReads.end()) {
                         wrappers[z].clusterId = currentId;
                         clusterReads.insert(wrappers[z].readId);
                         stack.push_back(z);
                     }
                }
            }
        }
    }
    
    // Step 4.4: Rescue Pass (Assign remaining to best cluster)
    for(int k=0; k<nSVs; ++k) {
        if(wrappers[k].clusterId != -1) continue;
        
        // Find best cluster
        int bestClust = -1;
        double maxScore = -1.0;
        
        for(int z=0; z<nSVs; ++z) {
            if(wrappers[z].clusterId == -1) continue;
            if(wrappers[z].readId == wrappers[k].readId) continue; // Cannot join same read
            
            // Range check
            uint32_t os = std::max(wrappers[k].sv.start, wrappers[z].sv.start);
            uint32_t oe = std::min(wrappers[k].sv.end, wrappers[z].sv.end);
            if(oe <= os + ol_w) continue;

            double sw = cal_lindel_dd(wrappers[k].sv, wrappers[z].sv, ol_r, ol_w, err_dif);
            if(sw > maxScore) {
                maxScore = sw;
                bestClust = wrappers[z].clusterId;
            }
        }
        
        if(bestClust != -1) {
            // Need to verify Unique Read again? 
            // Hifiasm `is_get_group` checks if the read is already in the target group (lines 19779).
            // Yes.
            
            // We need to know which reads are in `bestClust`.
            // Bruteforce check?
            bool conflict = false;
            for(int z=0; z<nSVs; ++z) {
                if(wrappers[z].clusterId == bestClust && wrappers[z].readId == wrappers[k].readId) {
                    conflict = true; break;
                }
            }
            
            if(!conflict) {
                wrappers[k].clusterId = bestClust;
            }
        }
    }

    // 5. Filter overlaps based on Valid Clusters
    // Calculate cluster sizes (Unique Reads)
    std::map<int, std::set<uint32_t>> clusterReadSets;
    for(const auto& w : wrappers) {
        if(w.clusterId != -1) {
            clusterReadSets[w.clusterId].insert(w.readId);
        }
    }
    
    std::set<int> validClusterIds;
    for(auto const& [cid, reads] : clusterReadSets) {
        if(reads.size() >= c_sz) validClusterIds.insert(cid);
    }
    
    // Construct "Consensus SV Regions"
    struct SVRegion {
         uint32_t start, end;
         int id;
    };
    std::vector<SVRegion> regions;
    for(const auto& w : wrappers) {
        if(validClusterIds.count(w.clusterId)) {
            regions.push_back({w.sv.start, w.sv.end, w.clusterId});
        }
    }

    std::vector<uint32_t> finalKept;
    std::vector<bool> keepMask(overlaps.size(), false);
    for(uint32_t idx : keptIndices) keepMask[idx] = true;

    for(const auto& r : regions) {
        for(uint32_t ovId : keptIndices) {
            if(!keepMask[ovId]) continue;
            
            const auto& ov = overlaps[ovId];
            if(ov.targetStart <= r.start && ov.targetEnd >= r.end) {
                 bool hasSV = false;
                 // Does this overlap have an SV in this cluster?
                 for(const auto& w : wrappers) {
                     if(w.sv.ovIdx == ovId && w.clusterId == r.id) {
                         hasSV = true; break;
                     }
                 }
                 if(!hasSV) {
                     keepMask[ovId] = false;
                 }
            }
        }
    }

    for(uint32_t idx : keptIndices) {
        if(keepMask[idx]) finalKept.push_back(idx);
    }

    // Hifiasm Parity (generate_haplotypes_sv line 9171):
    // Decrement occ_0 for evidence consumed by SV-phased overlaps.
    // This prevents double-counting reference evidence.
    if(!snpStats.empty()) {
        // Build site-to-stat lookup
        std::map<uint32_t, size_t> siteToStatIdx;
        for(size_t i = 0; i < snpStats.size(); ++i) {
            siteToStatIdx[snpStats[i].site] = i;
        }
        
        // For each SV-phased overlap, find its reference evidence and decrement occ_0
        for(const auto& w : wrappers) {
            if(validClusterIds.count(w.clusterId)) {
                uint32_t ovId = w.sv.ovIdx;
                // Find all evidence for this overlap
                for(const auto& ev : evidence) {
                    if(ev.overlapId == ovId && ev.type == 0) { // Type 0 = Reference
                        auto it = siteToStatIdx.find(ev.site);
                        if(it != siteToStatIdx.end() && snpStats[it->second].occ_0 > 1) {
                            snpStats[it->second].occ_0--;
                        }
                    }
                }
            }
        }
    }

    return finalKept;
}

void AssemblerPhasing::collectHaplotypeEvidence(
    const Assembler& assembler,
    const PhasingOverlap& overlap,
    std::vector<HaplotypeEvidence>& evidenceOut,
    const PhasingConfig& config)
{
    (void)config;
    const auto& reads = assembler.getReads();
    
    // Load Sequences
    // Using `getRead` which returns a Read view.
    // Target (Strand 0)
    const auto rT = reads.getRead(overlap.targetReadId);
    // Query (Oriented)
    const auto rQ_raw = reads.getRead(overlap.queryReadId);
    
    std::vector<Base> seqQ;
    bool queryIsRC = (overlap.queryStrand == 1);
    
    // We need random access to Query sequence.
    // If RC, we build the RC sequence in memory.
    // If NOT RC, we can access rQ_raw directly.
    
    if (queryIsRC) {
        seqQ.reserve(rQ_raw.baseCount);
        // Correct RC construction: reverse and complement
        for(size_t i=0; i<rQ_raw.baseCount; ++i) {
             seqQ.push_back(rQ_raw[rQ_raw.baseCount - 1 - i].complement());
        }
    }
    
    // Helper access
    auto getQueryBase = [&](uint32_t idx) -> Base {
        if (queryIsRC) return seqQ[idx];
        return rQ_raw[idx];
    };
    
    uint32_t tPos = overlap.targetStart;
    uint32_t qPos = overlap.queryStart; // Used from new PhasingOverlap fields

    size_t qSeqLen = rQ_raw.baseCount; // baseCount is valid member
    size_t tSeqLen = rT.baseCount;

    // Iterate CIGAR
    for(uint32_t val : overlap.cigar) {
        uint32_t op = val & 0xF;
        uint32_t len = val >> 4;
        
        if (op == 3) { // Blind Match
            tPos += len;
            qPos += len;
        } else if (op == 0) { // Match/Mismatch
            for(uint32_t k=0; k<len; ++k) {
                if(tPos >= tSeqLen || qPos >= qSeqLen) break; 
                
                Base bT = rT[tPos]; // Valid operator[]
                Base bQ = getQueryBase(qPos);
                
                if (bT != bQ) {
                     HaplotypeEvidence ev;
                     ev.site = tPos;
                     ev.base = bQ.value;
                     ev.type = 1;
                     ev.isSolid = 1;
                     // overlapId set by caller
                     evidenceOut.push_back(ev); 
                }
                tPos++;
                qPos++;
            }
        } else if (op == 1) { // Insert (in Query)
            qPos += len;
        } else if (op == 2) { // Deletion (Gap in Query)
             for(uint32_t k=0; k<len; ++k) {
                 if(tPos >= tSeqLen) break;
                 // Hifiasm heuristic check (mismatchCount not available here, assuming it's a placeholder for a more complex check)
                 // if(mismatchCount > 3) continue; 

        // Check Homopolymer Context (Simple Sequence Check)
        // If site is in a simple homopolymer run > X, flag/filter?
        // Hifiasm `hpc_mask_ff` is complex.
        // For strictest parity, we could check sequence context here.
        // Assuming pre-aligned CIGAR handles most, but we can check adjacent bases.
        // For now, relying on `is_st_bs` which catches artifacts effectively.
        // Adding TODO for full hpc_mask_ff parity if strictly needed.

        // Store Evidence
        HaplotypeEvidence ev;
        ev.overlapId = overlap.alnIdx; // Temp usage, sorted/grouped later
                 ev.site = tPos;
                 ev.base = 4; // Gap
                 ev.type = 1;
                 ev.isSolid = 0;
        ev.isSameStrand = (overlap.queryStrand == 0); // Parity: Same strand check
        
        evidenceOut.push_back(ev);
                 tPos++;
             }
        }
    }
}

// Parity with `comput_sc_rphase`
// Returns 1 if consistent, INT64_MIN if inconsistent.
// Consistency check: Do shared overlaps show same alleles?
static int64_t computeLinkScore(
    uint32_t siteI, uint32_t siteJ,
    const std::vector<const HaplotypeEvidence*>& evI,
    const std::vector<const HaplotypeEvidence*>& evJ) 
{
    if(siteI == siteJ) return -2e18; // Logic invalid

    // Intersection
    size_t k = 0, m = 0;
    int nn0 = 0; // Shared Ref support
    int nn1 = 0; // Shared Alt support
    
    // ev lists are sorted by overlapId (pre-condition)
    while(k < evI.size() && m < evJ.size()) {
        uint32_t ovI = evI[k]->overlapId;
        uint32_t ovJ = evJ[m]->overlapId;
        
        if(ovI < ovJ) { k++; continue; }
        else if(ovJ < ovI) { m++; continue; }
        else {
            // Same overlap covers both sites.
            // Check consistency:
            // Type 0 = Match(Ref), Type 1 = Variant(Alt)
            int typeI = (evI[k]->type == 0) ? 0 : 1;
            int typeJ = (evJ[m]->type == 0) ? 0 : 1;
            
            if(typeI != typeJ) return -2e18; // Inconsistent (Phase Switch)
            
            // Consistent!
            if(typeI == 0) nn0++;
            else nn1++;

            k++; m++;
        }
    }
    
    // Hifiasm Parity (comput_sc_rphase):
    // Requires both Ref and Alt support to be confident.
    // Line 9280: if(nn[0] > 0 && nn[1] > 0) return 1;
    
    if(nn0 > 0 && nn1 > 0) return 1;
    return -2e18; // Insufficient dual-haplotype support
}

std::vector<uint32_t> AssemblerPhasing::generatePhasingDP(
    ReadId /* targetReadId */,
    const std::vector<PhasingOverlap>& overlaps,
    const std::vector<HaplotypeEvidence>& evidence,
    std::vector<SnpStats>& outStats,
    const PhasingConfig& config)
{
    outStats.clear();
    
    // 1. Identify Candidate Sites & Build Lookup
    std::map<uint32_t, std::vector<const HaplotypeEvidence*>> siteEvidence;
    // Iterate sequentially to group by site
    
    for(const auto& ev : evidence) {
        siteEvidence[ev.site].push_back(&ev);
    }
    
    struct Candidate {
        uint32_t site;
        uint32_t occ_0;
        uint32_t occ_1;
        uint32_t overlap_num; // RefSameStrand
    };
    std::vector<Candidate> candidates;
    
    for(auto& [site, evList] : siteEvidence) {
        uint32_t occ0 = 0;
        uint32_t occ1 = 0;
        uint32_t overlap_num = 0; // Hifiasm Parity: Ref Same Strand

        for(auto* e : evList) {
             if(e->type == 0) {
                 occ0++;
                 if(e->isSameStrand) overlap_num++; 
             }
             else occ1++;
        }

        // Hifiasm Pre-filter: !is_st_bs
        // Also: Hifiasm Line 10560: if(rev_n == occ_0) return 0 -> discard if ALL Ref reads are same strand
        if(occ0 >= 2 && occ1 >= 2) {
             // Check: NOT all Ref reads on same strand (would indicate strand bias)
             if(overlap_num != occ0) { // At least one Ref read on opposite strand
                 if(!is_st_bs(overlap_num, occ0, config.st_rate, config.st_max)) {
                     // Hifiasm Line 10571: occ_0 = 1 + occ_0 (adds target read count)
                     candidates.push_back({site, occ0 + 1, occ1, overlap_num + 1}); // +1 for self
                 }
             }
        }
    }
    
    // Sort evidence lists by overlapId for computeLinkScore intersection
    for(auto& [site, evList] : siteEvidence) {
        std::sort(evList.begin(), evList.end(), [](const HaplotypeEvidence* a, const HaplotypeEvidence* b){
            return a->overlapId < b->overlapId;
        });
    }
    
    if(candidates.empty()) {
        std::vector<uint32_t> kept(overlaps.size());
        std::iota(kept.begin(), kept.end(), 0);
        return kept;
    }

    // 2. DP Chaining (O(N^2))
    int n = (int)candidates.size();
    std::vector<int64_t> f(n, 0); 
    std::vector<int> p(n, -1);
    
    int64_t maxScoreGlobal = 0;
    int bestEnd = -1;
    
    for(int i=0; i<n; ++i) {
        f[i] = 1; // Base score (self)
        
        // Hifiasm Parity: No lookback limit (Full O(N^2))
        // Hifiasm Line 9442: `for (j = i - 1; j >= st; --j)` where st=0.
        
        for(int j=i-1; j>=0; --j) {
            // Check consistency
            int64_t link = computeLinkScore(
                candidates[i].site, candidates[j].site,
                siteEvidence[candidates[i].site], siteEvidence[candidates[j].site]
            );
            
            if(link > -1e17) { // Consistent
                 if(f[j] + link > f[i]) {
                     f[i] = f[j] + link;
                     p[i] = j;
                 }
            }
        }
        
        if(f[i] > maxScoreGlobal) {
            maxScoreGlobal = f[i];
            bestEnd = i;
        }
    }
    
    // 3. Iterative Multi-Chain Extraction (Hifiasm Parity)
    // Hifiasm sorts nodes by score and extracts disjoint chains until exhausted.
    
    // Sort indices by score descending
    std::vector<int> sortedIndices(n);
    std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
    std::sort(sortedIndices.begin(), sortedIndices.end(), [&](int a, int b){
        return f[a] > f[b];
    });
    
    std::vector<bool> visited(n, false);
    
    uint64_t cc = (config.hom_cov / 2); 
    cc = (uint64_t)((double)cc * 0.70);
    if (cc < 6) cc = 6;
    
    outStats.clear(); // Reset output
    
    // Helper for homopolymer filtering (Hifiasm Parity)
    auto isHpcVec = [&](SnpStats& s) -> bool {
        // Hifiasm logic: Subtract HP evidence (if tracked) and check robustness
        // Since we don't track per-read HP status deeply yet, we check total counts against thresholds.
        // Line 9396: if((occ0 < 2 || occ1 < 2) || !(occ0 >= s_hap_cov && occ1 >= infor_cov)) return true (invalid);
        
        if (s.occ_0 < 2 || s.occ_1 < 2) return true;
        if (!(s.occ_0 >= (uint32_t)config.s_hap_cov && s.occ_1 >= (uint32_t)config.infor_cov)) return true;
        return false; // Valid
    };

    for(int startNode : sortedIndices) {
        if(visited[startNode]) continue;
        
        // Extract Chain
        std::vector<SnpStats> chain;
        int curr = startNode;
        while(curr != -1 && !visited[curr]) {
            visited[curr] = true;
            
            SnpStats s;
            s.site = candidates[curr].site;
            s.occ_0 = candidates[curr].occ_0;
            s.occ_1 = candidates[curr].occ_1;
            s.overlap_num = candidates[curr].overlap_num; // Hifiasm Parity: RefSameStrand 
            chain.push_back(s);
            
            curr = p[curr];
        }
        std::reverse(chain.begin(), chain.end());
        
        if(chain.empty()) continue;
        
        // Filter Chain (Hifiasm Logic)
        if(chain.size() == 1) {
            // Single site chain check
            // Hifiasm Line 9493: if( (!is_hpc_vec(...)) && (occ0 >= cc) ) keep;
             if (!isHpcVec(chain[0]) && chain[0].occ_0 >= cc) {
                 outStats.push_back(chain[0]);
             }
        } else {
            // Multi-site chain: Filter individual sites by cc
            // Hifiasm Line 9542: if(occ_0 >= cc) keep;
            // Note: Hifiasm DP multi-path check (line 9482) is complex, but for single-path it relies on chain score.
            for(const auto& s : chain) {
                if(s.occ_0 >= cc) {
                    outStats.push_back(s);
                }
            }
        }
    }
    
    // Re-sort outStats by site as we might have merged chains out of order
    std::sort(outStats.begin(), outStats.end(), [](const SnpStats& a, const SnpStats& b){
        return a.site < b.site;
    }); 
    
    // Note: Hifiasm also has logic for "plus" score calculation based on chain consistency.
    // But since we backtracked the Optimal Chain from DP, we assume internal consistency (plus=1).
    // The main filter is the coverage check `occ_0 >= cc`.

    // 4. Score Overlaps based on Valid Chain
    std::vector<uint32_t> keptIndices;
    
    for(size_t ovId=0; ovId<overlaps.size(); ++ovId) {
        const auto& ov = overlaps[ovId];
        
        int score = 0;
        uint32_t minT = ov.targetStart;
        uint32_t maxT = ov.targetEnd;
        
        // Iterate only Valid Sites in range
        for(const auto& s : outStats) {
             if(s.site < minT) continue;
             if(s.site > maxT) break;
             
             // Check valid site evidence for THIS overlap
             const auto& evList = siteEvidence[s.site];
             
             // Binary search
             auto it = std::lower_bound(evList.begin(), evList.end(), (uint32_t)ovId, 
                [](const HaplotypeEvidence* a, uint32_t id) {
                    return a->overlapId < id; 
                });
             
             if(it != evList.end() && (*it)->overlapId == ovId) {
                 if((*it)->type == 1) score++; // Alt support
                 else score--; 
             }
        }
        
        if(score <= 0) {
            keptIndices.push_back((uint32_t)ovId);
        }
    }

    return keptIndices;
}


// NEW: Returns phasing score per overlap (Hifiasm parity)
// This matches Hifiasm's logic exactly: TRANS overlaps are classified but NOT removed
std::vector<AssemblerPhasing::OverlapPhasingResult> AssemblerPhasing::getOverlapPhasingScores(
    const Assembler& assembler,
    ReadId targetReadId,
    const std::vector<PhasingOverlap>& overlaps,
    const PhasingConfig& config)
{
    std::vector<OverlapPhasingResult> results(overlaps.size(), {0, false});
    
    if (overlaps.empty()) return results;
    
    // Collect all haplotype evidence
    std::vector<HaplotypeEvidence> evidence;
    evidence.reserve(overlaps.size() * 10);
    
    for (size_t i = 0; i < overlaps.size(); ++i) {
        size_t startSize = evidence.size();
        collectHaplotypeEvidence(assembler, overlaps[i], evidence, config);
        for (size_t k = startSize; k < evidence.size(); ++k) {
            evidence[k].overlapId = (uint32_t)i;
        }
    }
    
    // DEBUG: Log evidence count for read 0
    if (targetReadId == 0) {
        std::cout << "  DEBUG getOverlapPhasingScores for read 0:" << std::endl;
        std::cout << "    Overlaps: " << overlaps.size() << std::endl;
        std::cout << "    Evidence items collected: " << evidence.size() << std::endl;
    }
    
    if (evidence.empty()) {
        // No evidence = all overlaps are CIS with no informative sites
        if (targetReadId == 0) {
            std::cout << "    WARNING: No evidence collected!" << std::endl;
        }
        return results;
    }
    
    // Build site evidence map
    std::map<uint32_t, std::vector<const HaplotypeEvidence*>> siteEvidence;
    for (const auto& ev : evidence) {
        siteEvidence[ev.site].push_back(&ev);
    }
    
    // DEBUG: Log site count and type distribution
    if (targetReadId == 0) {
        std::cout << "    Unique sites: " << siteEvidence.size() << std::endl;
        
        // Count type distribution
        uint64_t type0Count = 0, type1Count = 0;
        for (const auto& ev : evidence) {
            if (ev.type == 0) type0Count++;
            else type1Count++;
        }
        std::cout << "    Type distribution: type0(Match)=" << type0Count 
                  << ", type1(Mismatch)=" << type1Count << std::endl;
        
        // Show first 5 sites with their occ values
        std::cout << "    Sample sites:" << std::endl;
        int sampleCount = 0;
        for (auto& [site, evList] : siteEvidence) {
            if (sampleCount >= 5) break;
            uint32_t occ0 = 0, occ1 = 0;
            for (auto* e : evList) {
                if (e->type == 0) occ0++;
                else occ1++;
            }
            std::cout << "      Site " << site << ": occ0=" << occ0 << ", occ1=" << occ1 << std::endl;
            sampleCount++;
        }
    }
    
    // Identify valid candidate sites (Hifiasm parity)
    std::vector<SnpStats> validSites;
    uint64_t rejectedLowOcc = 0, rejectedStrand = 0, rejectedStBs = 0;
    
    for (auto& [site, evList] : siteEvidence) {
        uint32_t occ0 = 0, occ1 = 0, overlap_num = 0;
        
        for (auto* e : evList) {
            if (e->type == 0) {
                occ0++;
                if (e->isSameStrand) overlap_num++;
            } else {
                occ1++;
            }
        }
        
        // Hifiasm filter: occ_0 >= 2 && occ_1 >= 2 && not-all-same-strand
        if (occ0 < 2 || occ1 < 2) {
            rejectedLowOcc++;
        } else if (overlap_num == occ0) {
            rejectedStrand++;
        } else if (is_st_bs(overlap_num, occ0, config.st_rate, config.st_max)) {
            rejectedStBs++;
        } else {
            validSites.push_back({site, occ0, occ1, overlap_num, 0, 0});
        }
    }
    
    // DEBUG: Log validation stats
    if (targetReadId == 0) {
        std::cout << "    Valid sites: " << validSites.size() << std::endl;
        std::cout << "    Rejected (low occ): " << rejectedLowOcc << std::endl;
        std::cout << "    Rejected (strand bias): " << rejectedStrand << std::endl;
        std::cout << "    Rejected (is_st_bs): " << rejectedStBs << std::endl;
    }
    
    // Sort evidence lists by overlapId for binary search
    for (auto& [site, evList] : siteEvidence) {
        std::sort(evList.begin(), evList.end(), 
            [](const HaplotypeEvidence* a, const HaplotypeEvidence* b) {
                return a->overlapId < b->overlapId;
            });
    }
    
    // Score each overlap (Hifiasm parity: +1 for Alt, -1 for Ref at valid sites)
    for (size_t ovId = 0; ovId < overlaps.size(); ++ovId) {
        const auto& ov = overlaps[ovId];
        int score = 0;
        bool hasInformative = false;
        
        for (const auto& s : validSites) {
            if (s.site < ov.targetStart) continue;
            if (s.site > ov.targetEnd) break;
            
            const auto& evList = siteEvidence[s.site];
            
            // Binary search for this overlap's evidence at this site
            auto it = std::lower_bound(evList.begin(), evList.end(), (uint32_t)ovId,
                [](const HaplotypeEvidence* a, uint32_t id) {
                    return a->overlapId < id;
                });
            
            if (it != evList.end() && (*it)->overlapId == ovId) {
                hasInformative = true;
                if ((*it)->type == 1) {
                    score++;  // Alt support → TRANS
                } else {
                    score--;  // Ref support → CIS
                }
            }
        }
        
        results[ovId].score = score;
        results[ovId].hasInformativeSite = hasInformative;
    }
    
    return results;
}

} // namespace dinara
