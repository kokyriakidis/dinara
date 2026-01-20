
#include "Assembler.hpp"
#include "LongBaseSequence.hpp" 
#include "Reads.hpp"            
#include "timestamp.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <mutex>
#include <cstring> // For memset
#include "AlignedEvidenceStore.hpp" // Ensure full definition is available

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
    uint32_t overlap_num; // Ref Count on Strand 0 (Forward). For Strand Bias check.
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

// --------------------------------------------------------
// Function: detectHetSites (Parity for hc_phase_robust_rr)
// --------------------------------------------------------
template<typename AlignmentContainer>
static void detectHetSites(
    Assembler& assembler,
    const Reads& reads,
    uint64_t queryReadId,
    const vector<CandidateEC>& candidates, // "ol"
    const AlignmentContainer& alignmentData, // Access to raw alignment info
    vector<SnpStats>& snpStats,
    vector<HaplotypeEvidence>& hapEvidence
) {
    const LongBaseSequenceView queryRead = reads.getRead(ReadId(queryReadId));

    // 0. Pre-compute Homopolymer Map for Query (HPC_PL = 12 parity)
    vector<bool> isHP(queryRead.baseCount, false);
    {
        uint32_t rLen = 0;
        char prev = 0;
        for(uint32_t i=0; i<queryRead.baseCount; i++) {
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
             for(uint32_t j=queryRead.baseCount-rLen; j<queryRead.baseCount; j++) isHP[j] = true;
        }
    }
    
    // Loop candidates.
    for(size_t k=0; k<candidates.size(); k++) {
        const auto& cand = candidates[k];
        
        // Access Alignment Data
        const auto& ad = alignmentData[cand.alignmentId];
        size_t evidenceId = ad.info.alignmentId;
        
        if(evidenceId == invalid<size_t>) continue;

        // Logic:
        span<const SnpEvidence> snps;
        uint32_t currentPos = 0; 
        
        // Corrected ReadId vs Alignment ReadIds
        // alignmentData typically stores readIds[0] and readIds[1]
        // We need to check which one is our queryReadId
        
        if(ad.readIds[1] == queryReadId) {
            // Corrected Read is Target (Stream 0 relative to Evidence Store structure?)
            // Actually, AlignedEvidenceStore stores (read0 -> read1).
            // snps0 is deviations of read1 from read0? Or vice versa?
            // Assuming getSnps0 corresponds to Read 0? Needs verification from APES definition.
            // Based on usage: "Corrected Read is Target (Stream 0)" implies we want variants ON Target.
            snps = assembler.alignedEvidenceStore.getSnps0(evidenceId);
            currentPos = ad.ts;
        } else if(ad.readIds[0] == queryReadId) {
            // Corrected Read is Query (Stream 1)
            snps = assembler.alignedEvidenceStore.getSnps1(evidenceId);
            currentPos = ad.qs;
        } else {
            continue;
        }

        // Iterate Stream (Branchless Accumulation Logic)
        for(const auto& ev : snps) {
            currentPos += ev.delta();
            
            if(currentPos >= cand.qs && currentPos < cand.qe) {
                if(!isHP[currentPos]) {
                    HaplotypeEvidence he;
                    he.overlapID = (uint32_t)k;
                    he.site = currentPos;
                    he.type = 1; // Mismatch
                    he.misBase = ev.base();
                    he.overlapSite = 0; // Assigned later
                    hapEvidence.push_back(he);
                }
            }
        }
    }
    
    // Sort Evidence by Site
    std::sort(hapEvidence.begin(), hapEvidence.end());
    
    // Partition SnpStats
    if(hapEvidence.empty()) return;
    
    for(size_t i=0; i<hapEvidence.size(); ) {
        uint32_t site = hapEvidence[i].site;
        uint32_t misCount = 0;
        uint32_t startIdx = (uint32_t)i;
        while(i < hapEvidence.size() && hapEvidence[i].site == site) {
            misCount++;
            i++;
        }
        
        uint32_t totalCov = 0;
        uint32_t refReadsRev = 0; 
        
        for(size_t k=0; k<candidates.size(); k++) {
            const auto& cand = candidates[k];
            if(site >= cand.qs && site < cand.qe) {
                totalCov++;
                
                bool isAlt = false;
                for(size_t j=startIdx; j<i; j++) {
                     if(hapEvidence[j].overlapID == k) { isAlt = true; break; }
                }
                
                if(!isAlt) {
                     if(cand.isRev) refReadsRev++;
                }
            }
        }
        
        if(totalCov < 5) continue; 
        
        // Hifiasm Parity: occ_thres = 1.
        if(misCount > 1) { 
            SnpStats stat;
            stat.site = site;
            stat.occ_1 = misCount; 
            stat.occ_0 = totalCov - misCount;
            stat.overlap_num = refReadsRev; 
            
            for(size_t j=startIdx; j<i; j++) {
                hapEvidence[j].overlapSite = (uint32_t)snpStats.size();
            }
            snpStats.push_back(stat);
        }
    }
}

// --------------------------------------------------------
// Function: gen_rphase_dp
// --------------------------------------------------------
static void gen_rphase_dp(
    Assembler& assembler,
    vector<SnpStats>& snpStats,
    vector<HaplotypeEvidence>& hapEvidence,
    vector<CandidateEC>& candidates
) {
    if(snpStats.empty()) return;

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
                if(!std::binary_search(alts.begin(), alts.end(), (uint32_t)k)) {
                    siteReadsRef[i].push_back((uint32_t)k);
                }
            }
        }
    }
    
    // 2. Filter SNPs
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
        
        if(s.occ_0 < 2 || s.occ_1 < 2) { s.score = -1; continue; }
        if(is_st_bs(s)) { s.score = -1; continue; }
        if(s.occ_0 < 3 || s.occ_1 < 3) { s.score = -1; continue; }
        
        validSnpIndices.push_back((int)i);
    }
    
    if(validSnpIndices.empty()) return;

    // 3. Single Path DP
    double cut_rate = 0.7; 
    uint64_t cut_bd = 6;
    
    uint64_t peak = 0;
    if(assembler.assemblerInfo->kmerDistributionInfo.coveragePeak != invalid<uint64_t>) {
        peak = assembler.assemblerInfo->kmerDistributionInfo.coveragePeak;
    }
    if(peak == 0) peak = candidates.size(); 
    
    uint64_t cc = (uint64_t)(peak * cut_rate);
    if(cc < cut_bd) cc = cut_bd;

    int nValid = (int)validSnpIndices.size();
    vector<int64_t> f(nValid, 0);
    vector<int> p(nValid, -1);
    
    for(int idx=0; idx<nValid; idx++) {
        int i = validSnpIndices[idx];
        int64_t max_f = 1; 
        int max_j = -1;
        
        int lookbackLimit = 50; 
        for(int jdx=idx-1; jdx>=0; jdx--) {
            if(idx - jdx > lookbackLimit) break;
            int j = validSnpIndices[jdx];

            const auto& refsI = siteReadsRef[i];
            const auto& altsI = siteReadsAlt[i];
            const auto& refsJ = siteReadsRef[j];
            const auto& altsJ = siteReadsAlt[j];

            int nn0 = 0, nn1 = 0;
            int mm0 = 0, mm1 = 0;

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
            nn1 = count_commons(altsI, altsJ); 
            mm0 = count_commons(refsI, altsJ);
            mm1 = count_commons(altsI, refsJ);

            int64_t sc = (nn0 + nn1) - (mm0 + mm1)*5; 
            
            if (sc == INT64_MIN) continue; 
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
static void generate_haplotypes_naive_HiFi(
    Assembler& assembler,
    vector<SnpStats>& snpStats,
    vector<HaplotypeEvidence>& hapEvidence,
    vector<CandidateEC>& candidates
) {
    if(snpStats.empty()) return;
    
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
    
    for(size_t i=0; i<evByCand.size(); ) {
        uint32_t ovId = evByCand[i].ovId;
        int support = 0;
        
        // Unused `start` variable removed
        while(i < evByCand.size() && evByCand[i].ovId == ovId) {
            size_t idx = evByCand[i].evIdx;
            const auto& ev = hapEvidence[idx];
            
            if(ev.type == 1) { // Mismatch/SNP
                if(ev.overlapSite < snpStats.size()) {
                    const auto& s = snpStats[ev.overlapSite];
                    if(s.score == 1) {
                        support++;
                    }
                }
            }
            i++;
        }
        
        if(support > 0) {
            if(ovId < candidates.size()) {
                candidates[ovId].score = 1; 
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
    const vector<CandidateEC>& candidates,
    const AlignmentContainer& alignmentData, 
    vector<HaplotypeEvidence>& svEvidence,
    vector<SnpStats>& svStats 
) {
    // Parity SV Struct
    struct RawSV {
        uint32_t overlapID;
        uint32_t site; 
        int64_t size; 
    };
    vector<RawSV> rawSVs;

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
        vector<size_t> indices;
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
    vector<HaplotypeEvidence>& svEvidence, 
    vector<CandidateEC>& candidates
) {
    if(svEvidence.empty()) return;
    for(const auto& ev : svEvidence) {
        if(ev.overlapID < candidates.size()) {
            candidates[ev.overlapID].score = 1; 
        }
    }
}

// --------------------------------------------------------
// Helper to identify uncorrected reads
// --------------------------------------------------------
static bool isUncorrectedRead(const vector<CandidateEC>& candidates, uint64_t queryLen, int threshold) {
    if (candidates.empty()) return true;
    if (queryLen < (uint64_t)threshold) return true;
    return false;
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

            for(uint64_t r = start; r < end; r++) {
                // Unified Pipeline (TASSD Accelerated)
                uint32_t strand = 0;
                OrientedReadId oid(dinara::ReadId(r), strand);
                if(oid.getValue() >= alignmentTable.size()) continue;

                const auto& alignments = alignmentTable[oid.getValue()];
                if(alignments.empty()) continue;

                    // Collect candidates for this query
                    vector<CandidateEC> candidates;
                    candidates.reserve(alignments.size());
                    
                    for(uint32_t alignId : alignments) {
                        const auto& al = alignmentData[alignId];
                        
                        if(al.isDeleted()) continue;

                        CandidateEC c;
                        c.alignmentId = alignId;
                        
                        // Fix: Use readIds array
                        if(r == al.readIds[0]) {
                            c.targetId = al.readIds[1];
                            c.qs = al.qs; 
                            c.qe = al.qe;
                            c.ts = al.ts;
                            c.te = al.te;
                            c.isRev = !al.isSameStrand; 
                        } else {
                            c.targetId = al.readIds[0];
                            c.qs = al.ts; 
                            c.qe = al.te; 
                            c.ts = al.qs;
                            c.te = al.qe; 
                            c.isRev = !al.isSameStrand;
                        }
                        c.score = al.info.markerCount; 
                        
                        candidates.push_back(c);
                    }

                    if(candidates.empty()) continue;

                    // Strict Hifiasm Parity: Check for Uncorrected Read
                    uint64_t queryLen = reads->getRead(dinara::ReadId(r)).baseCount;
                    if (isUncorrectedRead(candidates, queryLen, 1600)) {
                         continue; 
                    }

                    // Sort by score descending
                    sort(candidates.begin(), candidates.end(), [](const CandidateEC& a, const CandidateEC& b) {
                        return a.score > b.score;
                    });
                    
                    // 1. Detect SNPs
                    vector<SnpStats> snpStats;
                    vector<HaplotypeEvidence> hapEvidence;
                    detectHetSites(*this, *reads, r, candidates, alignmentData, snpStats, hapEvidence);
                    
                    // 2. DP Step
                    gen_rphase_dp(*this, snpStats, hapEvidence, candidates);
                    
                    // 3. Naive Verification
                    generate_haplotypes_naive_HiFi(*this, snpStats, hapEvidence, candidates);

                    // 4. SV Step
                    vector<HaplotypeEvidence> svEvidence;
                    vector<SnpStats> svStats;
                    detectSVSites(*this, *reads, r, candidates, alignmentData, svEvidence, svStats);
                    generate_haplotypes_sv(*this, svEvidence, candidates);
                    
                    // Final Keep Decision
                    if(!candidates.empty()) {
                         int64_t bestScore = candidates[0].score;
                         double minScoreRatio = 0.5; 
                         
                         for(const auto& cand : candidates) {
                             if (cand.score >= bestScore * minScoreRatio) {
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
