#include "Assembler.hpp"
#include "ProjectedAlignment.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "Alignment.hpp"
#include "timestamp.hpp"
#include <algorithm>
#include <vector>

using namespace dinara;
using namespace std;

// ----------------------------------------------------------------------------
// Chimeric Read Detection
// ----------------------------------------------------------------------------
// Mirrors hifiasm's logic (gen_chemical_arc_rf -> cal_chemical_r_adv)
// to detect coverage gaps in "Cis" alignments.

void Assembler::detectChimericReads(uint64_t threadCount)
{
    cout << timestamp << "Detecting chimeric reads..." << endl;
    
    // Resize the chimeric flag vector
    isChimericRead.createNew(largeDataName("IsChimericRead"), largeDataPageSize);
    isChimericRead.resize(reads->readCount());
    std::fill(isChimericRead.begin(), isChimericRead.end(), false);

    // Access necessary data
    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();
    
    // Ensure alignmentTable is accessible
    if (!alignmentTable.isOpen()) {
        alignmentTable.accessExistingReadOnly(largeDataName("AlignmentTable"));
    }

    const uint64_t batchSize = 1;
    setupLoadBalancing(reads->readCount(), batchSize);
    runThreads(&Assembler::detectChimericReadsThreadFunction, threadCount);

    // Count chimeric reads
    uint64_t chimericCount = 0;
    for(size_t i=0; i<isChimericRead.size(); i++) {
        if(isChimericRead[i]) chimericCount++;
    }
    cout << timestamp << "Detected " << chimericCount << " chimeric reads." << endl;
}

void Assembler::detectChimericReadsThreadFunction(size_t /* threadId */)
{
    // Loop over batches of reads
    uint64_t readIdBegin, readIdEnd;
    while(getNextBatch(readIdBegin, readIdEnd)) {
        for(ReadId readId = ReadId(readIdBegin); readId < ReadId(readIdEnd); readId++) {
            
            // We only process one orientation (Strand 0) and assume it covers the read.
            OrientedReadId orientedReadId(readId, 0);
            const uint64_t readLen = reads->getReadRawSequenceLength(readId);
            
            // Gather valid intervals (Coverage Events)
            vector<pair<int32_t, int>> events;
            
            // Access alignments for this read
            // alignmentTable indexes alignments by OrientedReadId.
            // CAUTION: Check if alignmentTable covers this ID.
            if (orientedReadId.getValue() >= alignmentTable.size()) continue;
            
            // Using range-based iteration if supported, or size/indexing
            const size_t n = alignmentTable.size(orientedReadId.getValue());
            for(size_t i=0; i<n; i++) {
                const uint64_t alignmentId = alignmentTable[orientedReadId.getValue()][i];
                const AlignmentData& ad = alignmentData[alignmentId];
                
                // Only consider "Cis" (kept) alignments. 
                // If status is Unknown, we fallback to isInReadGraph.
                bool isCis = (bool)(ad.cisTransStatus == CisTransStatus::Cis);
                if (ad.cisTransStatus == CisTransStatus::Unknown) {
                    if (ad.info.isInReadGraph) isCis = true;
                }
                if (!isCis) continue;
                
                // Determine coordinates on THIS read (readId)
                // AlignmentData usually stores canonical pair or oriented pair.
                // We need to deal with the coordinate system properly.
                // Dinara stores alignment between oriented reads.
                // The alignment intervals are in `ad.info` usually (markers), 
                // but we added explicit `qs, qe, ts, te` (bases).
                
                uint32_t start = 0;
                uint32_t end = 0;
                ReadId otherReadId;
                
                if (ad.readIds[0] == readId) {
                    // Coordinates on Read 0 are usually Forward Strand relative?
                    // ad.qs and ad.qe are valid for Read 0.
                    start = ad.qs;
                    end = ad.qe;
                    otherReadId = ad.readIds[1];
                } else if (ad.readIds[1] == readId) {
                    // Read 1
                    start = ad.ts;
                    end = ad.te;
                    otherReadId = ad.readIds[0];
                } else {
                    continue; 
                }

                const uint64_t otherReadLen = reads->getReadRawSequenceLength(otherReadId);

                // --- DUPLICATE FILTER (Hifiasm-style) ---
                // Filter out full-length overlaps that look like duplicates.
                // If the overlap is nearly full length for both reads AND lengths are similar (PCR dup).
                // They support the chimera across the junction artificially.
                
                // Length diff
                uint64_t lenDiff = (readLen > otherReadLen) ? (readLen - otherReadLen) : (otherReadLen - readLen);
                
                // Overlap length on this read
                uint64_t overlapLen = end - start;
                
                // Overlap length on other read (approximation using same overlap length, 
                // or we could use the data fields if we knew orientation). 
                // Detailed check:
                // Hifiasm checks: if (diff <= len*0.02 && diff <= otherLen*0.02)
                // AND overlap covers nearly everything (unaligned part <= len*0.02).
                
                const double dupRate = 0.02;
                if (lenDiff <= uint64_t(double(readLen)*dupRate) && 
                    lenDiff <= uint64_t(double(otherReadLen)*dupRate)) {
                    
                    // Check if it's a full overlap
                    // Unaligned on this read: readLen - overlapLen
                    // We check if unaligned part is small.
                    uint64_t unaligned = readLen - overlapLen;
                    if (unaligned <= uint64_t(double(readLen)*dupRate)) {
                        // Likely a duplicate. Skip it.
                        continue; 
                    }
                }
                
                // --- SHRINK LOGIC (Hifiasm-style: Preserve Tips) ---
                // We shrink the valid coverage interval by flank (256bp), but only if the alignment doesn't touch the read end.
                // Why: Coverage naturally drops to zero at the tips of any read. We don't want to flag a read as chimeric just 
                // because it has low coverage at the very start or end. We only care about "holes" that appear deep inside the read body.
                // "if (s0 > 0) s0 += cut_len" -> If start is NOT at tip, shrink it (internal gap potential).
                // "if (e0 < len) e0 -= cut_len" -> If end is NOT at tip, shrink it.
                // If it IS at tip, keep it at tip (to provide coverage there).
                
                int32_t flank = 256; 
                int32_t s = int32_t(start);
                int32_t e = int32_t(end);
                
                if (s > 0) s += flank;
                if (e < int32_t(readLen)) e -= flank;
                
                if (s < e) {
                    events.push_back({s, 1});
                    events.push_back({e, -1});
                }
            }
            
            // --- SWEEP LOGIC ---
            // We sweep from 0 to readLen.
            // If min coverage < threshold anywhere in [0, readLen), it's chimeric?
            // Hifiasm tracks min_cov among "regions between events".
            // AND Hifiasm starts tracking `st=0`.
            
            if (events.empty()) {
                // No valid intervals -> Chimeric
                isChimericRead[readId] = true;
            } else {
                std::sort(events.begin(), events.end());
                
                int coverage = 0;
                const int minCoverageThreshold = 2; 
                bool foundGap = false;
                
                int prevPos = 0; 
                
                // Hifiasm Logic: Check ALL intervals [prevPos, pos)
                for(const auto& ev : events) {
                    int pos = ev.first;
                    
                    if (pos > prevPos) {
                        // Check coverage in [prevPos, pos)
                        if (coverage < minCoverageThreshold) {
                            foundGap = true; 
                            break;
                        }
                    }
                    
                    // Update coverage AFTER checking the interval
                    coverage += ev.second;
                    prevPos = pos;
                }
                
                // Check tail [lastEvent, readLen)
                if (!foundGap && prevPos < int32_t(readLen)) {
                     if (coverage < minCoverageThreshold) {
                        foundGap = true;
                     }
                }
                
                if (foundGap) {
                    isChimericRead[readId] = true;
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Chimeric Read Rescue (Symmetry Restoration)
// ----------------------------------------------------------------------------
// Mirrors hifiasm's `try_rescue_overlaps`.
// Logic:
// 1. Iterate over reads flagged as Chimeric.
// 2. Collect "Safe Neighbors" (Non-chimeric, Valid Alignment).
// 3. STRICT CHECK: Do these neighbors support the SAME region?
//    - Hifiasm performs a sweep-line on the neighbors' overlaps.
//    - It finds the maximum pileup depth (max_dp).
// 4. If max_dp >= 4 (RescueThreshold), un-flag the read.

void Assembler::rescueChimericReads(uint64_t threadCount)
{
    cout << timestamp << "Rescuing chimeric reads..." << endl;
    
    // We need to know who is chimeric to check neighbors. 
    setupLoadBalancing(reads->readCount(), 1);
    runThreads(&Assembler::rescueChimericReadsThreadFunction, threadCount);

    // Recount
    uint64_t chimericCount = 0;
    for(size_t i=0; i<isChimericRead.size(); i++) {
        if(isChimericRead[i]) chimericCount++;
    }
    cout << timestamp << "Chimeric reads after rescue: " << chimericCount << "." << endl;
}

void Assembler::rescueChimericReadsThreadFunction(size_t /* threadId */)
{
    const uint32_t rescueThreshold = 4; // Hifiasm default
    
    uint64_t readIdBegin, readIdEnd;
    while(getNextBatch(readIdBegin, readIdEnd)) {
        for(ReadId readId = ReadId(readIdBegin); readId < ReadId(readIdEnd); readId++) {
            
            // Only process if currently flagged as Chimeric
            if (!isChimericRead[readId]) continue;
            
            OrientedReadId orientedReadId(readId, 0);
            
            // Check neighbors in alignment table
            if (orientedReadId.getValue() >= alignmentTable.size()) continue;
            
            vector<pair<uint32_t, int>> events;
            
            const size_t n = alignmentTable.size(orientedReadId.getValue());
            
            for(size_t i=0; i<n; i++) {
                const uint64_t alignmentId = alignmentTable[orientedReadId.getValue()][i];
                const AlignmentData& ad = alignmentData[alignmentId];
                
                // Identify the neighbor
                ReadId neighborId;
                uint32_t qs, qe; // Coordinates on CURRENT read (readId)
                
                if (ad.readIds[0] == readId) {
                    neighborId = ad.readIds[1];
                    qs = ad.qs; qe = ad.qe;
                } else {
                    neighborId = ad.readIds[0];
                    qs = ad.ts; qe = ad.te;
                }
                
                // 1. Is Neighbor Chimeric?
                if (isChimericRead[neighborId]) continue;
                
                // 2. Is Alignment Valid?
                bool isValid = (ad.cisTransStatus == CisTransStatus::Cis);
                if (ad.cisTransStatus == CisTransStatus::Unknown && ad.info.isInReadGraph) {
                    isValid = true;
                }
                
                if (isValid) {
                    // Add Interval [qs, qe) to sweep line events
                    // Use Type 1=Start, 2=End to ensure Start < End at same position (Start processed first) - Matches Hifiasm
                    if (qs < qe) {
                        events.push_back({qs, 1});
                        events.push_back({qe, 2});
                    }
                }
            }
            
            // 3. Sweep Line to check Consensus Depth (max_dp)
            // We need at least 'rescueThreshold' reads agreeing on the SAME region.
            if (events.size() / 2 >= rescueThreshold) {
                // Sort events:
                // Hifiasm sorts `qs<<1` (Start, LSB=0) and `qe<<1|1` (End, LSB=1).
                // This means at the same coordinate, START (even) comes before END (odd).
                // We use TYPE_START=1, TYPE_END=2.
                // pair compare: first(pos) then second(type).
                // So {pos, 1} < {pos, 2}. Matches Hifiasm.
                std::sort(events.begin(), events.end());
                
                int max_dp = 0;
                int dp = 0;
                uint32_t start = 0;
                uint32_t max_interval_s = 0;
                uint32_t max_interval_e = 0;
                int old_dp = 0;

                for(const auto& ev : events) {
                    old_dp = dp;
                    
                    if (ev.second == 1) { // Start
                        dp++;
                    } else { // End
                        dp--;
                    }
                    
                    // Logic mirroring Hifiasm's max_interval tracking
                    if (old_dp < dp) { // Increasing (Start processed)
                        // "if(dp >= max_dp) { start = b.a[j]>>1; max_dp = dp; }"
                        if (dp >= max_dp) {
                            start = ev.first;
                            max_dp = dp;
                        }
                    } else if (old_dp > dp) { // Decreasing (End processed)
                        // "if(old_dp == max_dp) { max_interval.s = start; max_interval.e = b.a[j]>>1; }"
                        if (old_dp == max_dp) {
                            max_interval_s = start;
                            max_interval_e = ev.first;
                        }
                    }
                }
                
                // Hifiasm rescues if max_dp >= threshold.
                // It then filters edges to only those spanning max_interval.
                // "if(qs <= max_interval.s && qe >= max_interval.e)"
                
                if (max_dp >= (int)rescueThreshold) {
                    // Rescue!
                    isChimericRead[readId] = false;
                    
                    // Second Pass: Filter edges (Enforce Consensus)
                    // Hifiasm effectively drops edges that don't support the consensus.
                    // We mark them as deleted.
                    
                    for(size_t i=0; i<n; i++) {
                        const uint64_t alignmentId = alignmentTable[orientedReadId.getValue()][i];
                        // We need access to writable AlignmentData
                        AlignmentData& ad = alignmentData[alignmentId];
                        
                        // Identify neighbor & coordinates again
                        ReadId neighborId;
                        uint32_t qs, qe;
                        if (ad.readIds[0] == readId) {
                            neighborId = ad.readIds[1];
                            qs = ad.qs; qe = ad.qe;
                        } else {
                            neighborId = ad.readIds[0];
                            qs = ad.ts; qe = ad.te;
                        }
                        
                         // Same validity checks as pass 1
                        if (isChimericRead[neighborId]) continue; // Chimeric neighbors likely shouldn't constrain us, but Hifiasm ignores them in "evidence".
                        // Wait, Hifiasm only loops over "evi" (the collected valid neighbors) for filtering.
                        // So we should only filter the *valid* edges? 
                        // If an edge was skipped in pass 1 (e.g. invalid status), it remains whatever it was.
                        // But if we rescue the read, we keep "all" edges unless we delete them.
                        // So we MUST process all edges and delete those that fail the check?
                        // Hifiasm loop: "for (j = 0; j < evi.n; j++)" -> It only adds valid edges to paf. 
                        // The implicit implication is that edges NOT in evi are NOT added (remain deleted).
                        // So yes, we should probably delete *everything* that doesn't pass.
                        
                        bool isValidType = (ad.cisTransStatus == CisTransStatus::Cis);
                        if (ad.cisTransStatus == CisTransStatus::Unknown && ad.info.isInReadGraph) isValidType = true;
                        
                        // Strict filter:
                        // 1. Must be Valid Type
                        // 2. Must span max_interval
                        bool keep = false;
                        if (isValidType) {
                            if (qs <= max_interval_s && qe >= max_interval_e) {
                                keep = true;
                            }
                        }
                        
                        if (!keep) {
                            ad.isDeleted = true;
                        }
                    }
                }
            }
        }
    }
}


