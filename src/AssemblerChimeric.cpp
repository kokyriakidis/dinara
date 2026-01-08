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


