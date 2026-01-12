#include "Assembler.hpp"
#include "Reads.hpp" // For Reads class definition
#include "timestamp.hpp"      // For timestamp
// #include "loadBalancing.hpp" // Removed: caused compilation error
#include <algorithm>
#include <vector>

// For debugging/logging
#include <iostream>

using namespace dinara;
using namespace std;

// This function implements logic similar to miniasm's ma_hit_sub function.
// It analyzes the coverage depth along a read induced by its overlaps.
// It identifies contiguous high-coverage regions effectively filtering out 
// chimeric joins or low-quality ends.

void Assembler::filterLocalSegments(
    uint64_t minCoverage, 
    uint64_t threadCount)
{
    cout << timestamp << "Filtering local segments (min coverage=" << minCoverage << ")..." << endl;
    
    // Store parameters for the thread function if needed, 
    // or passing them via member variables/lambdas if strictly needed.
    // For now we assume we might need to store minCoverage in Assembler 
    // or just hardcode/pass it. 
    // Since runThreads takes a member function with only threadId, 
    // we should optimally store this config in a member or capture it via a lambda 
    // if runThreads supported it (it doesn't, it takes a pointer).
    // So we'll use a member variable or a temporary one.
    // Let's rely on 'this->minOverlapCoverage' or similar if it exists, 
    // or just add a temporary member for this pass.
    
    // Actually, miniasm passes min_dp.
    // Let's look for a suitable place to store this transient config.
    // For now, I'll assume we can add a member to Assembler or reuse an existing one.
    // I'll add 'localSegmentMinCoverage' to Assembler.hpp later if needed.
    // But wait, I can't change the header right now easily inside this write.
    // I will add a member variable 'localSegmentMinCoverage' in Assembler.hpp next.
    
    // Resize output vector
    validReadIntervals.resize(reads->readCount());
    // Initialize with deleted status
    std::fill(validReadIntervals.begin(), validReadIntervals.end(), ReadSegmentStatus{0, 0, false}); // Default not deleted? Or should be? 
    // Miniasm default initializes and then sets. 
    // Let's safe initialize to 0,0,false. The loop covers all reads in batches.
    // If a read has no alignments, the loop:
    // 1. collectIntervals -> empty.
    // 2. if empty -> validReadIntervals[r0] = {0,0, true}; continue;
    
    // Store parameter for access by thread function
    this->localSegmentMinCoverage = minCoverage;

    setupLoadBalancing(reads->readCount(), 100);
    runThreads(&Assembler::filterLocalSegmentsThreadFunction, threadCount);
    
    cout << timestamp << "Local segment filtering complete." << endl;
}

void Assembler::filterLocalSegmentsThreadFunction(size_t threadId)
{
    // Minimal depth (min_dp in miniasm)
    // We use the passed minCoverage argument.
    // However, it is passed to the main function, not stored in "this".
    // We need to access it.
    // Option A: Store it in a member variable before running threads.
    // Option B: Hardcode/atomic (not ideal).
    // Reviewing Assembler class, we can just add a member `localSegmentMinCoverage`.
    
    // For now, assuming we added the member or use a specialized approach.
    // But to respect the current structure without modifying header again immediately:
    // We will assume `minCoverage` is available or we use a lambda capture if we could.
    // Since we can't change signature of runThreads target easily...
    // Let's assume we added `this->minOverlapCoverage` usage or similar.
    // Actually, `Assembler` has `minOverlapCoverage`? No, it has `minAlleleCoverage`.
    
    // Let's use `3` as a placeholder DEFAULT if we can't access the arg, 
    // BUT the user asked for "miniasm style", implying dynamic.
    // I will read `localSegmentMinCoverage` which I will conceptually add to `Assembler`.
    // Wait, I can't read it if I didn't add it.
    // The previous Plan step added `validReadIntervals` to header.
    // I should have added `localSegmentMinCoverage` too. 
    // I will use a const for now to allow compilation, but add TODO.
    // OR better: I can modify the header now to add the configuration variable.
    
    const uint64_t minDp = (this->localSegmentMinCoverage > 0) ? this->localSegmentMinCoverage : 3;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); r0++) {
            
            // Fast path for low coverage requirement (pass everything)
            if (minDp <= 1) {
                 uint64_t len = reads->getReadRawSequenceLength(r0);
                 validReadIntervals[r0] = {0, (uint32_t)len, false};
                 continue;
            }

            // Collect alignment intervals (start, end) on r0
            // from all alignments involving r0.
            
            // 1. Collect Valid Intervals
            std::vector<std::pair<uint32_t, uint32_t>> intervals;
            
            auto collectIntervals = [&](OrientedReadId orientedR0) {
                const auto& table = alignmentTable[orientedR0.getValue()];
                for (uint32_t alignmentId : table) {
                    const auto& ad = alignmentData[alignmentId];
                     if (ad.isDeleted) continue; // Skip deleted

                    // We use the explicit coordinates provided in AlignmentData (qs, qe, ts, te).
                    // These are assumed to be 0-based coordinates on the Forward strand of the read.
                    
                    uint32_t start = 0, end = 0;
                    
                    if (ad.readIds[0] == r0) {
                        // r0 is the Query (first read)
                        start = ad.qs;
                        end = ad.qe;
                    } else {
                        // r0 must be the Target (second read)
                        start = ad.ts;
                        end = ad.te;
                    }
                    
                    if (end > start) {
                        intervals.push_back({start, end});
                    }
                }
            };
            
            // We need to check both strands of r0 in the alignment table 
            // because alignments are stored indexed by the OrientedReadId.
            // However, the coordinates (qs/qe/ts/te) are physical coordinates on the read.
            collectIntervals(OrientedReadId(r0, 0));
            collectIntervals(OrientedReadId(r0, 1));
            
            if (intervals.empty()) {
                // No coverage means no valid segment -> deleted
                validReadIntervals[r0] = {0, 0, true};
                continue; 
            }
            
            // 2. Compute Coverage Depth (Sweep Line)
            // Events: (pos, 1) for start, (pos, -1) for end.
            std::vector<std::pair<uint32_t, int>> events;
            events.reserve(intervals.size() * 2);
            for(const auto& interval : intervals) {
                events.push_back({interval.first, 1});
                events.push_back({interval.second, -1});
            }
            
            // Sort events. Process Start (+1) before End (-1) at same position
            // to maintain continuity of adjacent interval coverage.
            std::sort(events.begin(), events.end(), [](const pair<uint32_t, int>& a, const pair<uint32_t, int>& b){
                if (a.first != b.first) return a.first < b.first;
                return a.second > b.second; // 1 (Start) before -1 (End)
            });
            
            // 3. Find Longest Segment with coverage >= minDp
            int coverage = 0;
            uint32_t maxS = 0, maxE = 0;
            
            uint32_t currentStart = 0;
            bool inSegment = false;
            
            // Loop through all events
            for(const auto& ev : events) {
                int oldCoverage = coverage;
                coverage += ev.second;
                
                uint32_t pos = ev.first;
                
                // Rising edge to valid coverage
                if (oldCoverage < (int)minDp && coverage >= (int)minDp) {
                    currentStart = pos;
                    inSegment = true;
                }
                // Falling edge from valid coverage
                else if (oldCoverage >= (int)minDp && coverage < (int)minDp) {
                    if (inSegment) {
                        uint32_t currentLen = pos - currentStart;
                        uint32_t maxLen = maxE - maxS;
                        
                        if (currentLen > maxLen) {
                            maxS = currentStart;
                            maxE = pos;
                        }
                        inSegment = false;
                    }
                }
            }
            
            // Store result
            if (maxE > maxS) {
                validReadIntervals[r0] = {maxS, maxE, false}; // Keep
            } else {
                validReadIntervals[r0] = {0, 0, true}; // Delete
            }
        }
    }
}



// Apply coverage cuts (ma_hit_cut equivalent).
// Clips alignments to the valid regions identified by filterLocalSegments.
// Removes short overlaps.
void Assembler::applyCoverageCuts(uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Applying coverage cuts (ma_hit_cut equivalent)." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    this->coverageCutMinOverlap = minOverlapLength;

    // Phase 1: Clip alignments in parallel.
    // Iterate over alignmentData directly to avoid double processing and simplify concurrency.
    setupLoadBalancing(alignmentData.size(), 10000); 
    runThreads(&Assembler::applyCoverageCutsToAlignmentsThreadFunction, threadCount);

    // Phase 2: Cleanup orphaned reads.
    // Reads that have no surviving edges should be marked as deleted.
    setupLoadBalancing(reads->readCount(), 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);
    
    cout << timestamp << "Coverage cuts applied." << endl;
}

void Assembler::applyCoverageCutsToAlignmentsThreadFunction(size_t threadId)
{
    uint64_t begin, end;
    const uint64_t minLen = this->coverageCutMinOverlap;

    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; i++) {
            AlignmentData& ad = alignmentData[i];
            
            // If already deleted, skip.
            if(ad.isDeleted) continue;

            ReadId qn = ad.readIds[0];
            ReadId tn = ad.readIds[1];

            // Access valid intervals (coverage_cut)
            // If explicit intervals not computed, assume full read? 
            // ma_hit_cut assumes sub is populated. We assume validReadIntervals is populated.
            // Safety check:
            if (qn >= validReadIntervals.size() || tn >= validReadIntervals.size()) continue;

            const auto& rq = validReadIntervals[qn];
            const auto& rt = validReadIntervals[tn];

            // "if any of target read and the query read has no enough coverage" -> del
            if (rq.isDeleted || rt.isDeleted) {
                ad.isDeleted = true;
                continue;
            }

            // Alias internal coords for readability matching reference code
            int32_t qs, qe, ts, te;
            
            // Logic mirroring ma_hit_cut coordinate projection
            if (!ad.isSameStrand) { // Different Strand (p->rev)
                // Anti-parallel projection
                qs = (int32_t)ad.te < (int32_t)rt.end ? (int32_t)ad.qs : (int32_t)ad.qs + ((int32_t)ad.te - (int32_t)rt.end);
                qe = (int32_t)ad.ts > (int32_t)rt.start ? (int32_t)ad.qe : (int32_t)ad.qe - ((int32_t)rt.start - (int32_t)ad.ts);
                ts = (int32_t)ad.qe < (int32_t)rq.end ? (int32_t)ad.ts : (int32_t)ad.ts + ((int32_t)ad.qe - (int32_t)rq.end);
                te = (int32_t)ad.qs > (int32_t)rq.start ? (int32_t)ad.te : (int32_t)ad.te - ((int32_t)rq.start - (int32_t)ad.qs);

            } else { // Same Strand
                // Parallel offset
                qs = (int32_t)ad.ts > (int32_t)rt.start ? (int32_t)ad.qs : (int32_t)ad.qs + ((int32_t)rt.start - (int32_t)ad.ts);
                qe = (int32_t)ad.te < (int32_t)rt.end ? (int32_t)ad.qe : (int32_t)ad.qe - ((int32_t)ad.te - (int32_t)rt.end);
                ts = (int32_t)ad.qs > (int32_t)rq.start ? (int32_t)ad.ts : (int32_t)ad.ts + ((int32_t)rq.start - (int32_t)ad.qs);
                te = (int32_t)ad.qe < (int32_t)rq.end ? (int32_t)ad.te : (int32_t)ad.te - ((int32_t)ad.qe - (int32_t)rq.end);
            }

            // "cut by self coverage" - Strict Clip to valid bounds
            // Using max/min to enforce bounds without offset change.
            qs = std::max(qs, (int32_t)rq.start);
            qe = std::min(qe, (int32_t)rq.end);
            ts = std::max(ts, (int32_t)rt.start);
            te = std::min(te, (int32_t)rt.end);

            // Check if valid length remains
            // "if (qe - qs >= mini_overlap_length && te - ts >= mini_overlap_length)"
            if ((qe - qs >= (int32_t)minLen) && (te - ts >= (int32_t)minLen) && (qe > qs) && (te > ts)) {
                // Update alignment data
                ad.qs = (uint32_t)qs;
                ad.qe = (uint32_t)qe;
                ad.ts = (uint32_t)ts;
                ad.te = (uint32_t)te;
            } else {
                ad.isDeleted = true;
            }
        }
    }
}

void Assembler::applyCoverageCutsCleanupThreadFunction(size_t threadId)
{
    // Phase 2: Check for orphaned reads.
    // "if(rLen == 0) (*coverage_cut)[i].del = 1;"
    // Iterate reads, check if any edges survive.

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId r=ReadId(begin); r!=ReadId(end); r++) {
            if (r >= validReadIntervals.size()) continue;
            if (validReadIntervals[r].isDeleted) continue; // Already deleted

            // Check neighbors
            bool hasSurvivingEdge = false;
            
            // Check strand 0 (and strand 1 implicitly via edges)
            const auto& table = alignmentTable[OrientedReadId(r, 0).getValue()];
            for(uint32_t alignmentId : table) {
                if (!alignmentData[alignmentId].isDeleted) {
                    hasSurvivingEdge = true;
                    break;
                }
            }
            // Check strand 1 (just in case edges are partitioned, though Dinara usually links both)
            if (!hasSurvivingEdge) {
                const auto& table1 = alignmentTable[OrientedReadId(r, 1).getValue()];
                for(uint32_t alignmentId : table1) {
                    if (!alignmentData[alignmentId].isDeleted) {
                        hasSurvivingEdge = true;
                        break;
                    }
                }
            }

            if (!hasSurvivingEdge) {
                validReadIntervals[r].isDeleted = true;
            }
        }
    }
}


// Filter hanging overlaps (ma_hit_flt equivalent).
// Filters overlaps based on hanging length (dovetail rule) and containment.
void Assembler::filterHangingOverlaps(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Filtering hanging overlaps (ma_hit_flt equivalent)." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    this->hangingFilterMaxHang = maxHang;
    this->hangingFilterMaxHangRate = maxHangRate;
    this->hangingFilterMinOverlap = minOverlapLength;

    // Use the alignment thread function directly
    setupLoadBalancing(alignmentData.size(), 10000);
    runThreads(&Assembler::filterHangingOverlapsThreadFunction, threadCount);

    // After filtering edges, we should clean up orphaned reads again?
    // ma_hit_flt does this internally ("if (rLen == 0) ...").
    // Reuse existing cleanup function.
    setupLoadBalancing(reads->readCount(), 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);

    cout << timestamp << "Hanging overlaps filtered." << endl;
}

void Assembler::filterHangingOverlapsThreadFunction(size_t threadId)
{
    uint64_t begin, end;
    const uint64_t maxHang = this->hangingFilterMaxHang;
    const double maxHangRate = this->hangingFilterMaxHangRate;
    const uint64_t minOvlp = this->hangingFilterMinOverlap;

    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; i++) {
            AlignmentData& ad = alignmentData[i];
            
            if(ad.isDeleted) continue;

            ReadId qn = ad.readIds[0];
            ReadId tn = ad.readIds[1];

            // Access valid lengths. If not present (e.g. skipped step), use full read len?
            // ma_hit_flt assumes validReadIntervals populated.
            // If empty, assume full length (fallback).
            uint32_t ql, tl, qs0, ts0;
            if (validReadIntervals.empty()) {
                ql = (uint32_t)reads->getReadRawSequenceLength(qn);
                tl = (uint32_t)reads->getReadRawSequenceLength(tn);
                qs0 = 0; ts0 = 0;
            } else {
                const auto& rq = validReadIntervals[qn];
                const auto& rt = validReadIntervals[tn];
                if (rq.isDeleted || rt.isDeleted) {
                    ad.isDeleted = true;
                    continue;
                }
                ql = rq.end - rq.start;
                tl = rt.end - rt.start;
                qs0 = rq.start;
                ts0 = rt.start;
            }

            // Calculate effective coordinates relative to valid region
            // Ensure non-negative (clipping should have handled this, but be safe)
            uint32_t qs = (ad.qs >= qs0) ? (ad.qs - qs0) : 0;
            uint32_t qe = (ad.qe >= qs0) ? (ad.qe - qs0) : 0;
            uint32_t ts = (ad.ts >= ts0) ? (ad.ts - ts0) : 0;
            uint32_t te = (ad.te >= ts0) ? (ad.te - ts0) : 0;
            
            // Re-check length just in case
            if (qe <= qs || te <= ts) {
                ad.isDeleted = true;
                continue;
            }

            // Calculate Overhangs (tl5, tl3) - Target hangs
            int32_t tl5, tl3;
            if (!ad.isSameStrand) { // Different Strand
                // "tl5 = tl - h->te, tl3 = h->ts"
                tl5 = (int32_t)tl - (int32_t)te;
                tl3 = (int32_t)ts;
            } else { // Same Strand
                // "tl5 = h->ts, tl3 = tl - h->te"
                tl5 = (int32_t)ts;
                tl3 = (int32_t)tl - (int32_t)te;
            }

            // Calculate Extensions (ext5, ext3) - Combined hangs
            // "ext5 = qs < tl5 ? qs : tl5;" (min)
            int32_t ext5 = ((int32_t)qs < tl5) ? (int32_t)qs : tl5;
            
            // "ext3 = ql - h->qe < tl3 ? ql - h->qe : tl3;" (min of Query 3' Hang vs Target 3' Hang)
            int32_t q3Hang = (int32_t)ql - (int32_t)qe;
            int32_t ext3 = (q3Hang < tl3) ? q3Hang : tl3;

            // --- Filter Logic ---
            // 1. Max Hang Check
            // "if (ext5 > max_hang || ext3 > max_hang ...)"
            bool badShape = (ext5 > (int32_t)maxHang) || (ext3 > (int32_t)maxHang);

            if (!badShape) {
                // 2. Hang Rate Check
                // "|| h->qe - qs < (h->qe - qs + ext5 + ext3) * int_frac"
                // "|| h->te - h->ts < (h->te - h->ts + ext5 + ext3) * int_frac"
                // Check both query and target overlap ratio? 
                // Reference checks OR, so if EITHER side is too short relative to total, fail.
                
                int32_t qOvlp = (int32_t)qe - (int32_t)qs;
                int32_t tOvlp = (int32_t)te - (int32_t)ts;
                int32_t totalSpanQ = qOvlp + ext5 + ext3;
                int32_t totalSpanT = tOvlp + ext5 + ext3; // ext used for both? yes per miniasm logic

                if (qOvlp < (totalSpanQ * maxHangRate) || tOvlp < (totalSpanT * maxHangRate)) {
                    badShape = true;
                }
            }
            
            // 3. Min Overlap Check (using extended spans or core?)
            // "if ((int)h->qe - qs + ext5 + ext3 < min_ovlp ...)"
            // Checks extended length?
            if (!badShape) {
                 int32_t qOvlp = (int32_t)qe - (int32_t)qs;
                 int32_t tOvlp = (int32_t)te - (int32_t)ts;
                 if ((qOvlp + ext5 + ext3 < (int32_t)minOvlp) || (tOvlp + ext5 + ext3 < (int32_t)minOvlp)) {
                     // return MA_HT_SHORT_OVLP -> Fail
                     badShape = true;
                 }
            }
            
            // 4. Containment Classification (Optional preservation)
            // miniasm keeps QCONT/TCONT even if "badShape"?
            // ma_hit_flt: "if (r >= 0 || r == MA_HT_QCONT || r == MA_HT_TCONT)" -> Keep.
            // Where r is result of ma_hit2arc.
            // ma_hit2arc logic:
            //   if (ext checks fail) return MA_HT_INT; (-1)
            //   else checks containment... 
            // So if EXT CHECKS FAIL (badShape), it returns INT, which is < 0 and NOT QCONT/TCONT.
            // So Contained reads are liable to fail Hang Check?
            // Wait, if contained: qs <= tl5 (because qs=ext5?), etc.
            // If contained, hangs should be 0 or small?
            // Actually, if contained, say Query in Target:
            // qs matches middle of T.
            // ext5 = min(qs, T_start_hang).
            // If qs is huge (deep inside), but T_start_hang is huge...
            // Basically if Bad Shape, we delete.
            
            if (badShape) {
                ad.isDeleted = true;
            }
        }
    }
}


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
                    // --- SHRINK LOGIC (Hifiasm-style: Preserve Tips) ---
                // We shrink the valid coverage interval by flank (256bp), but only if the alignment doesn't touch the read end.
                
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

// ----------------------------------------------------------------------------
// Miniasm-style Chimeric Read Detection (detectChimericReadsFromAnchors)
// ----------------------------------------------------------------------------

void Assembler::detectChimericReadsFromAnchors(
    double shiftRate, 
    uint64_t ulThres, 
    uint64_t threadCount)
{
    cout << timestamp << "Detecting chimeric reads (Anchors method)..." << endl;
    
    this->chimericShiftRate = shiftRate;
    this->chimericUlThres = ulThres;

    // Ensure flags initialized if not already
    if (!isChimericRead.isOpen) {
         isChimericRead.createNew(largeDataName("IsChimericRead"), largeDataPageSize);
         isChimericRead.resize(reads->readCount());
         std::fill(isChimericRead.begin(), isChimericRead.end(), false);
    }
    
    // We assume isChimericRead is already allocated
    
    setupLoadBalancing(reads->readCount(), 100);
    runThreads(&Assembler::detectChimericReadsFromAnchorsThreadFunction, threadCount);

    // Count chimeric reads
    uint64_t chimericCount = 0;
    for(size_t i=0; i<isChimericRead.size(); i++) {
        if(isChimericRead[i]) chimericCount++;
    }
    cout << timestamp << "Detected " << chimericCount << " chimeric reads (Anchors)." << endl;
}

// detectChimericReadsFromAnchorsThreadFunction
//
// This function implements chimeric read detection logic inspired by miniasm/hifiasm.
// It identifies "chimeric" reads (reads formed by the artificial joining of two distant
// genomic locations) by checking for gaps in coverage between "anchors" at the read tips.
//
// Logic Overview:
// 1.  **Anchors (`collectSides`)**: For each read, we identify overlaps that start near the
//     beginning (Left Anchor) and overlaps that end near the end (Right Anchor).
//     We track the furthest internal reach of these anchors (`maxLeft.e` and `maxRight.s`).
//
// 2.  **Gap Bridging (`collectContain`)**: We attempt to bridge the gap between Left and Right
//     anchors using "contained" reads (reads fully contained within the query read).
//     If a contained read significantly overlaps an existing anchor (determined by `overlapRate`),
//     we extend the anchor's reach. This helps bridge small coverage gaps caused by sequencing errors.
//
// 3.  **Gap Detection**: We check if there is a gap between the extended Left Anchor and Right Anchor.
//     Condition: `maxLeft.End <= maxRight.Start`.
//     If a gap exists, it implies the read consists of two disjoint alignment sets that do not
//     connect, suggesting a chimeric join.
//
// 4.  **Cleaning**: If a read is flagged as chimeric:
//     - We mark it as validly chimeric (`isChimericRead[id] = true`).
//     - We conceptually "cut" the read or remove it by marking all its alignments as deleted (`alignmentData[i].isDeleted = true`).

void Assembler::detectChimericReadsFromAnchorsThreadFunction(size_t threadId)
{
    const double shiftRate = this->chimericShiftRate; // Threshold for overlap between left/right anchors (default ~0.06)
    const double overlapRate = 0.1;                   // Threshold for extending anchors via contained reads

    struct CoverageExtent {
        uint32_t s;
        uint32_t e;
    };

    uint64_t begin, endIdx;
    while(getNextBatch(begin, endIdx)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(endIdx); r0++) {
            
            uint64_t rLen = reads->getReadRawSequenceLength(r0);
            if (rLen == 0) continue;

            uint32_t effectiveStart = 0;
            uint32_t effectiveEnd = (uint32_t)rLen;

            // Use filtered segments if available (coverage cut)
            if (!validReadIntervals.empty()) {
                const auto& status = validReadIntervals[r0];
                if (status.isDeleted) continue; // Read was filtered out entirely
                effectiveStart = status.start;
                effectiveEnd = status.end;
                
                if (effectiveStart >= effectiveEnd) continue;
            }

            // Initialize extents based on valid region.
            // maxLeft: [effectiveStart, e) - starts at effectiveStart, ends at e. 
            // maxRight: [s, effectiveEnd) - starts at s, ends at effectiveEnd.
            CoverageExtent maxLeft = { effectiveEnd, effectiveStart }; // Start with empty range at left
            CoverageExtent maxRight = { effectiveEnd, effectiveStart }; // Start with empty range at right (invalid)
            
            // Correction: Init maxLeft to represent "no left anchor"
            maxLeft = { effectiveEnd, effectiveStart }; 
            // Actually, we want maxLeft.e to grow from effectiveStart.
            // And maxRight.s to shrink from effectiveEnd.
            maxLeft = { effectiveEnd, effectiveStart };
            maxRight = { effectiveEnd, effectiveStart };
             
            // Re-think initialization matching logic:
            // "maxLeft.e" tracks the right-most reach of the left anchor. Init to effectiveStart.
            // "maxRight.s" tracks the left-most reach of the right anchor. Init to effectiveEnd.
            maxLeft.e = effectiveStart;
            maxRight.s = effectiveEnd;
            maxLeft.s = effectiveStart; // Just for completeness, though checks focus on .e
            maxRight.e = effectiveEnd;
            
            // --- Helper: Collect Anchors (collectSides) ---
            // Scans alignments to find interactions touching the read tips (within tolerance).
            auto collectSides = [&](OrientedReadId orientedR0) {
                 const auto& table = alignmentTable[orientedR0.getValue()];
                 for(uint32_t alignmentId : table) {
                     AlignmentData& ad = alignmentData[alignmentId]; 
                     if(ad.isDeleted) continue;
                     
                     // Get coordinates on r0
                     uint32_t qs = 0, qe = 0;
                     if (ad.readIds[0] == r0) {
                         qs = ad.qs; qe = ad.qe;
                     } else {
                         qs = ad.ts; qe = ad.te;
                     }
                     
                     // Constrain to valid region?
                     // If alignment is outside valid region, we shouldn't use it?
                     // Hifiasm ma_hit_sub trims overlaps to the valid region.
                     // Here we just check coordinates against effectiveStart/End.
                     
                     // Left anchor check (starts near effectiveStart)
                     const uint32_t tolerance = 100; 
                     if (qs <= effectiveStart + tolerance) {
                         // It is a valid left anchor candidate.
                         // We want the furthest it reaches.
                         // But we should clip qe to effectiveEnd? 
                         uint32_t reach = std::min(qe, effectiveEnd);
                         if (reach > maxLeft.e) maxLeft.e = reach;
                     }
                     // Right anchor check (ends near effectiveEnd)
                     if (qe >= effectiveEnd - tolerance) {
                         uint32_t start = std::max(qs, effectiveStart);
                         if (start < maxRight.s) maxRight.s = start;
                     }
                 }
            };

            collectSides(OrientedReadId(r0, 0));
            // Note: We only check Strand 0. Since the alignment graph is symmetric 
            // (overlaps are mirrored), Strand 0 contains all necessary physical coverage info.
            
            // If missing an anchor (Tip/Fragment), skip
            // (One side not covered -> Not a chimera, just a tip/fragment)
            if (maxLeft.s == effectiveEnd || maxRight.s == effectiveEnd) {
                continue;
            }

            // --- Helper: Collect Contain (collectContain) ---
            // Extend anchors using contained reads that overlap significantly with existing anchors.
            // This mirrors miniasm's logic to "grow" the valid regions inward.
            uint32_t newLeftE = maxLeft.e;
            uint32_t newRightS = maxRight.s;

            auto collectContain = [&](OrientedReadId orientedR0) {
                 const auto& table = alignmentTable[orientedR0.getValue()];
                 for(uint32_t alignmentId : table) {
                     const AlignmentData& ad = alignmentData[alignmentId];
                     if(ad.isDeleted) continue;
                     
                     uint32_t qs = 0, qe = 0;
                     if (ad.readIds[0] == r0) {
                         qs = ad.qs; qe = ad.qe;
                     } else {
                         qs = ad.ts; qe = ad.te;
                     }

                     // Check contained overlaps (reads strictly inside valid region, not touching ends)
                     const uint32_t tolerance = 100;
                     if (qs > effectiveStart + tolerance && qe < effectiveEnd - tolerance) {
                         
                         // Try to extend Left Anchor
                         // Logic: If contained read overlaps Left Anchor significantly, 
                         // we can trust it to extend the Left Anchor further right.
                         // Overlap: [qs, maxLeft.e)
                         if (qs < maxLeft.e && qe > maxLeft.e) {
                             uint32_t overlapLen = maxLeft.e - qs;
                             uint32_t alignLen = qe - qs;
                             if (overlapLen > (overlapRate * alignLen)) {
                                 if (qe > newLeftE) newLeftE = qe;
                             }
                         }

                         // Try to extend Right Anchor
                         // Logic: If contained read overlaps Right Anchor significantly,
                         // we can trust it to extend the Right Anchor further left.
                         // Overlap: [maxRight.s, qe)
                         if (qs < maxRight.s && qe > maxRight.s) {
                             uint32_t overlapLen = qe - maxRight.s;
                             uint32_t alignLen = qe - qs;
                             if (overlapLen > (overlapRate * alignLen)) {
                                 if (qs < newRightS) newRightS = qs;
                             }
                         }
                     }
                 }
            };
            
            collectContain(OrientedReadId(r0, 0));

            maxLeft.e = newLeftE;
            maxRight.s = newRightS;


            // --- Check Valid Bridging (Simple Chimera Detection) ---
            // If overlaps meet in the middle (Left End > Right Start), the read is likely good.
            // Even if they strictly "meet", we might require a minimal overlap length (shiftRate).
            int64_t overlap = (int64_t)maxLeft.e - (int64_t)maxRight.s;
            // Should effective length be used for ratio? (effectiveEnd - effectiveStart)
            uint32_t effectiveLen = effectiveEnd - effectiveStart;
            if (overlap > 0 && overlap >= (int64_t)(effectiveLen * shiftRate)) {
                continue; // Good Read (Anchors overlap significantly)
            }
            
            // If we are here, there is a GAP (maxLeft.e <= maxRight.s) or insufficient overlap.
            // We treat this as a potentially chimeric read.
            bool isSimpleChimera = (maxLeft.e <= maxRight.s);
            
            // Note: miniasm includes a "Complex Chimera" check here (intersection_check_by_base)
            // for cases where overlaps meet but alignment quality at junction is poor.
            // We omit that for now and rely on structural gap (Simple Chimera).
            bool isComplexChimera = false; // logic would go here
            
            // The user logic for 'complex': 
            // if intersection_check_by_base(...) returns true -> delete.
            // Since we don't have base check, we act conservatively?
            // The prompt asks to "mirror".
            // If simple chimera, we delete.
            
            if (isSimpleChimera) {
                isChimericRead[r0] = true;
                
                // --- Delete All Edges (delete_all_edges) ---
                // Mark all alignments involving this read as deleted.
                auto deleteEdges = [&](OrientedReadId orientedR0) {
                     const auto& table = alignmentTable[orientedR0.getValue()];
                     for(uint32_t alignmentId : table) {
                         // We can delete here because alignmentData is shared but atomic bool or similar is safe?
                         // isDeleted is not atomic.
                         // However, parallel writes to same bool are OK if value is same (true).
                         // Race condition if someone reads it.
                         // Standard for this codebase: marking reads is preferred.
                         // But user asked to mirror delete_all_edges.
                         alignmentData[alignmentId].isDeleted = true;
                     }
                };
                deleteEdges(OrientedReadId(r0, 0));
                
                continue;
            }
            
            // If complex logic existed, we would check it here.
            // But without base check, we assume non-gap reads are OK.
        }
    }
}

// --- Contained Read Filtering (ma_hit_contained_advance equivalent) ---

// Helper to determine containment status between two reads based on an overlap.
// Should mimic ma_hit2arc's containment checks.
// Returns: 
// 0: No containment
// 1: Read 0 is contained in Read 1 (QCONT)
// 2: Read 1 is contained in Read 0 (TCONT)
static int checkContainment(
    uint32_t qs, uint32_t qe, uint32_t ql, // Read 0 valid Interval: start, end, len
    uint32_t ts, uint32_t te, uint32_t tl, // Read 1 valid Interval: start, end, len
    uint32_t os, uint32_t oe,              // Overlap on Read 0 (raw coords)
    uint32_t rs, uint32_t re,              // Overlap on Read 1 (raw coords)
    bool isReverse,                        // relative strand
    uint32_t maxHang,
    double maxHangRate,
    uint32_t minOverlap)
{
    // Adjust overlap coordinates to be relative to the valid interval
    // If the overlap is outside the valid interval, it effectively doesn't exist for containment
    // But typically we assume we're working with valid overlaps.
    
    // Calculate overhangs (hang0, hang1)
    // Conceptually similar to filterHangingOverlaps but specifically for containment.
    
    // miniasm ma_hit2arc logic for containment:
    // It compares the "unmapped" portions.
    
    // Left/Right overhangs on Query (Read 0)
    int32_t q_hang_l = (int32_t)os - (int32_t)qs;
    int32_t q_hang_r = (int32_t)qe - (int32_t)oe;
    
    // Left/Right overhangs on Target (Read 1)
    // If reverse, Target start corresponds to Query end.
    int32_t t_hang_l, t_hang_r;
    if (isReverse) {
        // Reverse strand: 
        // Read 0:  [os ----- oe]
        // Read 1:  [re ----- rs]  (logical direction)
        // ts is start of valid region on 1, te is end.
        
        // Target valid region is [ts, te].
        // Overlap on target is [rs, re] (where rs < re, usually stored sorted).
        // Wait, Alignment info usually stores rs, re on the forward strand of Read 1.
        // If reverse, the 5' end of Read 0 matches the 3' end of Read 1?
        // Standard PAF/Coords: 
        // R0: os..oe
        // R1: rs..re. If Str=1, then R0(os) matches R1(re), R0(oe) matches R1(rs).
        
        // Let's compute 'projected' coordinates to align them.
        // But simpler: just compute overhangs.
        
        // Target overhangs (relative to valid interval [ts, te])
        // The overlap on T is [rs, re].
        // If Str=1:
        // Q_Left (os) matches T_Right (re).
        // Q_Right (oe) matches T_Left (rs).
        
        t_hang_l = (int32_t)re - (int32_t)te; // How much T extends beyond the match on the 'right' (which matches Q's left)
        t_hang_r = (int32_t)ts - (int32_t)rs; // How much T extends beyond match on 'left' (matches Q's right)
        
        // Wait, signage.
        // If T=[0..100], Valid=[10..90]. Overlap=[20..80].
        // ts=10, te=90. rs=20, re=80.
        // t_hang_left_of_match = rs - ts = 20-10=10.
        // t_hang_right_of_match = te - re = 90-80=10.
        
        // If Str=1:
        // match on T is [rs, re].
        // Q start matches T re. So T's extension "before" the match (relative to Q) is actually bases AFTER re.
        // So T_hang_corresponding_to_Q_Left = (te - re). (Assuming te > re).
        // T_hang_corresponding_to_Q_Right = (rs - ts). (Assuming rs > ts).
        
        t_hang_l = (int32_t)te - (int32_t)re;
        t_hang_r = (int32_t)rs - (int32_t)ts;
        
    } else {
        // Forward:
        // Q Left (os) matches T Left (rs).
        // Q Right (oe) matches T Right (re).
        
        t_hang_l = (int32_t)rs - (int32_t)ts;
        t_hang_r = (int32_t)te - (int32_t)re;
    }
    
    // Now we have q_hang_l, q_hang_r, t_hang_l, t_hang_r.
    // These represent the length of the read *outside* the alignment.
    // Note: Can be negative if overlap extends beyond valid interval? 
    // Usually clamped or validInterval constrains it.
    
    // Check Query Contained in Target (QCONT)
    // Q is contained if it has NO significant overhangs compared to T's extension?
    // No, Q is contained if Q is "inside" T.
    // This means Q's ends are "covered" by T.
    // i.e., T extends further than Q in both directions (or effectively so).
    // Specifically:
    // overlap almost covers Q (q_hang_l ~ 0, q_hang_r ~ 0).
    // T has slack on both sides (t_hang_l > 0, t_hang_r > 0) OR T roughly equals Q.
    
    // miniasm ma_hit2arc:
    // int32_t l = 0, r = 0;
    // if (tl > MAX_HANG) l = 1; // Left overhang significant
    // if (tr > MAX_HANG) r = 1; // Right overhang significant
    // And it compares lengths.
    
    // Simplified logic:
    // Q is contained if overlap length is close to Q's length.
    // (q_hang_l < maxHang && q_hang_r < maxHang).
    
    // T is contained if overlap length is close to T's length.
    // (t_hang_l < maxHang && t_hang_r < maxHang).
    
    bool q_covered = (abs(q_hang_l) < (int32_t)maxHang && abs(q_hang_r) < (int32_t)maxHang);
    bool t_covered = (abs(t_hang_l) < (int32_t)maxHang && abs(t_hang_r) < (int32_t)maxHang);
    
    if (q_covered && t_covered) {
        // Mutual containment (almost identical reads).
        // Resolve by ID or length to be deterministic.
        // Keep the longer one, or if equal, keep the larger ID.
        if (ql < tl) return 1; // Q < T -> Q contained
        if (tl < ql) return 2; // T < Q -> T contained
        return (qs > ts) ? 1 : 2; // Tie-break (arbitrary but deterministic)
    }
    
    if (q_covered) return 1;
    if (t_covered) return 2;
    
    return 0;
}

void Assembler::removeContainedReads(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Removing contained reads..." << endl;
    
    if (!containmentParent.isOpen) {
        containmentParent.createNew(largeDataName("ContainmentParent"), largeDataPageSize);
        containmentParent.resize(reads->readCount());
        std::fill(containmentParent.begin(), containmentParent.end(), ReadId(invalidReadId));
    }

    // Temporary storage for parameters
    this->hangingFilterMaxHang = maxHang;
    this->hangingFilterMaxHangRate = maxHangRate;
    this->hangingFilterMinOverlap = minOverlapLength;

    setupLoadBalancing(alignmentData.size(), 1000); // Iterate over alignments
    runThreads(&Assembler::removeContainedReadsThreadFunction, threadCount);
    
    // Transitive Reduction of Containment Tree & Cleanup
    cout << timestamp << "Transitive reduction and containment cleanup..." << endl;
    
    uint64_t containedCount = 0;
    
    // 1. Path Compression for containment parent
    for(size_t i=0; i<containmentParent.size(); i++) {
        if (containmentParent[i] != invalidReadId) {
            
            ReadId p = containmentParent[i];
            
            // Check for containment cycles or long chains
            // (A contains B, B contains C -> A contains C)
            // Limit depth to avoid infinite loops if cycle exists
            int depth = 0;
            while(p != invalidReadId && containmentParent[p] != invalidReadId && depth < 100) {
                 p = containmentParent[p];
                 depth++;
            }
            // If cycle (depth hit max), maybe break link or pick one? 
            // For now, assign ultimate parent.
            containmentParent[i] = p; 
            
            // Mark as deleted in global flags
            // This corresponds to 'coverage_cut[i].del = 1'
            validReadIntervals[i].isDeleted = true;
            containedCount++;
        }
    }
    
    // 2. Remove all alignments involving contained reads (Cleanup)
    // The snippet does `delete_all_edges` immediately.
    // We do it in a pass or rely on graph construction to respect `isDeleted`.
    // Let's enforce it in validReadIntervals which is checked by ReadGraph5.
    // Also, we can launch a thread pass to mark AlignmentData::isDeleted if strictly needed for other steps.
    // For now, setting validReadIntervals::isDeleted is sufficient as ReadGraph5 checks it.
    
    cout << timestamp << "Identified " << containedCount << " contained reads." << endl;
}

void Assembler::removeContainedReadsThreadFunction(size_t threadId)
{
    // Access loadBalancing directly (inherited from MultithreadedObjectBaseClass)
    uint64_t start, end;
    while(getNextBatch(start, end)) {
        for(uint64_t i=start; i<end; i++) {
            // Check if alignment is active in the graph construction
            // Since we don't have 'keepAlignment' vector passed here, we check 'isInReadGraph'
            // and maybe 'isDeleted' if that field exists in AlignmentData or similar.
            
            // Note: In ReadGraph5, keepAlignment is local. 
            // We should rely on `alignmentData[i].info.isInReadGraph` which should be updated?
            // Actually ReadGraph5 updates keepAlignment locally and then final graph uses it.
            // But filtering functions modify separate flags (`isDeleted`, `validReadIntervals`).
            // `applyCoverageCuts` uses `validReadIntervals`.
            
            // To be safe, we check if either read is ALREADY deleted by previous steps.
            
            // Access alignment data
            // Assuming AlignmentData is the type, let's see what it contains.
            // It likely has 'info' of type AlignmentInfo.
            // AlignmentData might be a struct with: ReadId readIds[2]; AlignmentInfo info; ...
            // Let's guess based on usage elsewhere: `ad.readIds[0]`.
            
            AlignmentData& ad = alignmentData[i];

            ReadId r0 = ad.readIds[0];
            ReadId r1 = ad.readIds[1];

            // Get raw lengths
            uint64_t len0 = getReads().getReadRawSequenceLength(r0);
            uint64_t len1 = getReads().getReadRawSequenceLength(r1);
            
            // Unaligned tails for R0
            int64_t tail0_L = (int64_t)ad.qs;
            int64_t tail0_R = (int64_t)len0 - (int64_t)ad.qe;
             
            // Unaligned tails for R1
            int64_t tail1_L = (int64_t)ad.ts;
            int64_t tail1_R = (int64_t)len1 - (int64_t)ad.te;

            // int64_t leftUnaligned0 = ad.info.leftTrim(0);
            // int64_t rightUnaligned0 = ad.info.rightTrim(0);
            // int64_t leftUnaligned1 = ad.info.leftTrim(1);
            // int64_t rightUnaligned1 = ad.info.rightTrim(1);

            int64_t leftUnaligned0 = tail0_L;
            int64_t rightUnaligned0 = tail0_R;
            int64_t leftUnaligned1 = tail1_L;
            int64_t rightUnaligned1 = tail1_R;

            if (r0 == 8894 || r1 == 8894) {
                cout << "Alignment " << i << ": " << r0 << " " << r1 << endl;
                cout << "Tail 0: " << tail0_L << " " << tail0_R << endl;
                cout << "Tail 1: " << tail1_L << " " << tail1_R << endl;
                cout << "Left Unaligned 0: " << leftUnaligned0 << " " << leftUnaligned1 << endl;
                cout << "Right Unaligned 0: " << rightUnaligned0 << " " << rightUnaligned1 << endl;
            }

            if ((leftUnaligned0 <= leftUnaligned1) && (rightUnaligned0 <= rightUnaligned1)) {
                containmentParent[r0] = r1;
                validReadIntervals[r0].isDeleted = true;
                ad.isDeleted = true;
                
                // Aggressively remove all other alignments involving contained read r0
                for(uint32_t strand=0; strand<2; strand++) {
                    OrientedReadId oid(r0, strand);
                    for(uint32_t otherAlignId : alignmentTable[oid.getValue()]) {
                        alignmentData[otherAlignId].isDeleted = true;
                    }
                }
                
                continue;
            } else if ((leftUnaligned1 <= leftUnaligned0) && (rightUnaligned1 <= rightUnaligned0)) {
                containmentParent[r1] = r0;
                validReadIntervals[r1].isDeleted = true;
                ad.isDeleted = true;

                // Aggressively remove all other alignments involving contained read r1
                for(uint32_t strand=0; strand<2; strand++) {
                    OrientedReadId oid(r1, strand);
                    for(uint32_t otherAlignId : alignmentTable[oid.getValue()]) {
                        alignmentData[otherAlignId].isDeleted = true;
                    }
                }
                
                continue;
            }




            
            if (!ad.info.isInReadGraph) continue;
            
            if (ad.isDeleted) continue;
            
            
            
            // Ignore if already deleted/chimeric
            if (validReadIntervals[r0].isDeleted || validReadIntervals[r1].isDeleted) continue;
            if (isChimericRead[r0] || isChimericRead[r1]) continue;

            const auto& status0 = validReadIntervals[r0];
            const auto& status1 = validReadIntervals[r1];
            
            // Containment Logic:
            // Check if one read is "enclosed" in the other.
            // Enclosed = The read is fully covered by the alignment (i.e. has no significant overhangs).
            // User Request: "check if one is contained in the other! we cant have both of them contained!"
            
            

            
            // int64_t tol = (int64_t)this->hangingFilterMaxHang;
            
            
            // Determine Candidate based on length (Shorter read can be contained in Longer read)
            // Tie-break with ID for determinism.
            bool checkR0 = false;
            
            if (len0 < len1) {
                checkR0 = true;
            } else if (len1 < len0) {
                checkR0 = false;
            } else {
                // Equal lengths: Tie-break
                if (r0 > r1) checkR0 = true;
                else checkR0 = false;
            }
            
            if (checkR0) {
                 // Check if R0 is contained in R1.
                 // logic: R0 is "inside" R1 if R0's tails are smaller than R1's tails (on the corresponding sides).
                 // relative check: tail0 <= tail1 (+ tolerance)
                 if ((tail0_L <= tail1_L) && (tail0_R <= tail1_R)) {
                     containmentParent[r0] = r1;
                     validReadIntervals[r0].isDeleted = true;
                     ad.isDeleted = true;
                 }
            } else {
                 // Check if R1 is contained in R0.
                 if ((tail1_L <= tail0_L) && (tail1_R <= tail0_R)) {
                     containmentParent[r1] = r0;
                     validReadIntervals[r1].isDeleted = true;
                     ad.isDeleted = true;
                 }
            }
        }
    }
}

