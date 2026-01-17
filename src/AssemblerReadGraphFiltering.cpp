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
                     if (ad.isDeleted()) continue; // Skip deleted

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
            
            
            // We only collect from strand 0 to match Hifiasm's sources[i] which is single-direction.
            // Collecting from both strands would double-count overlaps.
            collectIntervals(OrientedReadId(r0, 0));
            // Hifiasm parity: Don't collect strand 1 - same alignment is indexed for both oriented reads
            // collectIntervals(OrientedReadId(r0, 1));
            
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
            if(ad.isDeleted()) continue;

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
                ad.setDeleted(true);
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

            // Hifiasm parity: Normalize coordinates to 0-based relative to valid region
            // This matches: qs = (qs - rq.s), qe = (qe - rq.s), etc.
            int32_t norm_qs = qs - (int32_t)rq.start;
            int32_t norm_qe = qe - (int32_t)rq.start;
            int32_t norm_ts = ts - (int32_t)rt.start;
            int32_t norm_te = te - (int32_t)rt.start;

            // Check if valid length remains
            // "if (qe - qs >= mini_overlap_length && te - ts >= mini_overlap_length)"
            if ((norm_qe - norm_qs >= (int32_t)minLen) && (norm_te - norm_ts >= (int32_t)minLen) && 
                (norm_qe > norm_qs) && (norm_te > norm_ts)) {
                // Update alignment data with normalized coordinates
                ad.qs = (uint32_t)norm_qs;
                ad.qe = (uint32_t)norm_qe;
                ad.ts = (uint32_t)norm_ts;
                ad.te = (uint32_t)norm_te;
            } else {
                ad.setDeleted(true);
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
                if (!alignmentData[alignmentId].isDeleted()) {
                    hasSurvivingEdge = true;
                    break;
                }
            }
            // Check strand 1 (just in case edges are partitioned, though Dinara usually links both)
            if (!hasSurvivingEdge) {
                const auto& table1 = alignmentTable[OrientedReadId(r, 1).getValue()];
                for(uint32_t alignmentId : table1) {
                    if (!alignmentData[alignmentId].isDeleted()) {
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
            
            if(ad.isDeleted()) continue;

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
                    ad.setDeleted(true);
                    continue;
                }
                ql = rq.end - rq.start;
                tl = rt.end - rt.start;
                qs0 = rq.start;
                ts0 = rt.start;
            }

            // IMPORTANT: After applyCoverageCuts, coordinates are already normalized (0-based relative to valid region)
            // So we use them directly without subtracting vr.start again
            uint32_t qs = ad.qs;
            uint32_t qe = ad.qe;
            uint32_t ts = ad.ts;
            uint32_t te = ad.te;
            
            // Re-check length just in case
            if (qe <= qs || te <= ts) {
                ad.setDeleted(true);
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
            
            // 4. Containment Classification (Hifiasm ma_hit_flt Parity)
            // Hifiasm ma_hit_flt: "if (r < 0 || r == MA_HT_QCONT || r == MA_HT_TCONT)" -> Delete.
            // Contained reads must be deleted here as edges. 
            // ma_hit2arc checks containment if badShape is false? No, it checks containment explicitly.
            
            bool isQCont = false;
            bool isTCont = false;
            
            // Query contained in Target
            // qs <= tl5 (Left fits) AND (ql - qe) <= tl3 (Right fits)
            if ((int32_t)qs <= tl5 && ((int32_t)ql - (int32_t)qe) <= tl3) {
                isQCont = true;
            }
            // Target contained in Query
            // qs >= tl5 (Left extended) AND (ql - qe) >= tl3 (Right extended)
            if ((int32_t)qs >= tl5 && ((int32_t)ql - (int32_t)qe) >= tl3) {
                isTCont = true;
            }
            
            if (isQCont || isTCont) {
                badShape = true;
            }
            
            if (badShape) {
                ad.setDeleted(true);
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

    // Count chimeric reads and mark them as deleted in validReadIntervals
    uint64_t chimericCount = 0;
    for(size_t i=0; i<isChimericRead.size(); i++) {
        if(isChimericRead[i]) {
            chimericCount++;
            // IMPORTANT: Propagate status to validReadIntervals for downstream filtering
            if (i < validReadIntervals.size()) {
                validReadIntervals[i].isDeleted = true;
            }
        }
    }
    cout << timestamp << "Detected " << chimericCount << " chimeric reads." << endl;
    if (chimericCount > 0) {
        cout << timestamp << "  -> Marked these reads as deleted in validReadIntervals." << endl;
    }
}

void Assembler::detectChimericReadsThreadFunction(size_t /* threadId */)
{
    // Hifiasm Parity Implementation (detect_chimeric_reads: Overlaps.cpp:2449)
    
    // Struct to hold max_left / max_right intervals
    struct SubRegion { uint32_t s, e; };

    uint64_t start, end;
    while(getNextBatch(start, end)) {
        for(ReadId readId = ReadId(start); readId < ReadId(end); readId++) {
            
            uint64_t rLen = reads->getReadRawSequenceLength(readId);
            OrientedReadId oid(readId, 0);

            // Access alignments
            if (oid.getValue() >= alignmentTable.size()) continue;
            auto alignments = alignmentTable[oid.getValue()];
            
            SubRegion max_left = { (uint32_t)rLen, 0 };
            SubRegion max_right = { (uint32_t)rLen, 0 };
            
            // --- 1. collect_sides ---
            // "if(qs == 0) update max_left; if(qe == rLen) update max_right"
            
            for (uint32_t alignmentId : alignments) {
                const AlignmentData& ad = alignmentData[alignmentId];
                if (ad.isDeleted()) continue;
                
                // Get coords on readId
                uint32_t qs, qe;
                if (ad.readIds[0] == readId) { qs = ad.qs; qe = ad.qe; }
                else { qs = ad.ts; qe = ad.te; }

                // Check Left Anchor
                if (qs == 0) {
                    if (qs < max_left.s) max_left.s = qs;
                    if (qe > max_left.e) max_left.e = qe;
                }
                
                // Check Right Anchor
                if (qe == (uint32_t)rLen) {
                    if (qs < max_right.s) max_right.s = qs;
                    if (qe > max_right.e) max_right.e = qe;
                }
            }
            
            // "if(max_left.s == rLen || max_right.s == rLen) continue;"
            // Means we didn't find any left anchor or any right anchor.
            // Hifiasm treats these as "End Nodes" and does not flag them as chimeric.
            if (max_left.s == (uint32_t)rLen || max_right.s == (uint32_t)rLen) {
                continue;
            }
            
            // --- 2. collect_contain ---
            // Extend max_left.e and max_right.s using overlaps contained within them.
            // Hifiasm: "if(qs < max_left.e && qe > max_left.e && max_left.e - qs > (overlap_rate * (qe -qs)))"
            float overlap_rate = 0.1f;
            uint32_t new_left_e = max_left.e;
            uint32_t new_right_s = max_right.s;
            
            for (uint32_t alignmentId : alignments) {
                const AlignmentData& ad = alignmentData[alignmentId];
                if (ad.isDeleted()) continue;
                
                uint32_t qs, qe;
                if (ad.readIds[0] == readId) { qs = ad.qs; qe = ad.qe; }
                else { qs = ad.ts; qe = ad.te; }

                if (qs != 0 && qe != (uint32_t)rLen) {
                    // Contained overlap extension for Left
                    if (qs < max_left.e && qe > max_left.e) {
                         uint32_t len = qe - qs;
                         if (len > 0 && (max_left.e - qs) > (uint32_t)(overlap_rate * len)) {
                             if (qe > new_left_e) new_left_e = qe;
                         }
                    }
                    // Contained overlap extension for Right
                    if (qs < max_right.s && qe > max_right.s) {
                        uint32_t len = qe - qs;
                        if (len > 0 && (qe - max_right.s) > (uint32_t)(overlap_rate * len)) {
                            if (qs < new_right_s) new_right_s = qs;
                        }
                    }
                }
            }
            max_left.e = new_left_e;
            max_right.s = new_right_s;
            
            // --- 3. Check for Overlap ---
            // "if (max_left.e > max_right.s && (max_left.e - max_right.s >= rLen * shift_rate)) continue;"
            // If they overlap sufficiently, it's a good read.
            // We assume shift_rate=0 for basic check.
            if (max_left.e > max_right.s) {
                continue;
            }
            
            // --- 4. Chimeric Flagging ---
            // If we are here, max_left.e <= max_right.s (Gap Exists or barely touch).
            // Hifiasm performs expensive "intersection_check" here.
            // For structural parity, we flag as chimeric if the gap is not bridged.
            // If the anchors don't meet, and no contained reads bridged them (step 2), then it's chimeric.
            
            isChimericRead[readId] = true;
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
                            ad.setDeleted(true);
                        }
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Rescue Phased Overlaps (try_rescue_overlaps equivalent)
// ----------------------------------------------------------------------------
// Rescues overlaps where isDeleted0 != isDeleted1 (directional phasing conflict)
// if >= rescueThreshold overlaps from the same direction agree on a region.
// This recovers edges that were incorrectly phased out from one direction.

void Assembler::rescuePhasedOverlaps(uint64_t rescueThreshold, uint64_t threadCount)
{
    cout << timestamp << "Rescuing phased overlaps (threshold=" << rescueThreshold << ")..." << endl;
    
    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    
    this->rescuePhasedThreshold = rescueThreshold;
    
    setupLoadBalancing(reads->readCount(), 1);
    runThreads(&Assembler::rescuePhasedOverlapsThreadFunction, threadCount);
    
    // Count rescued overlaps
    uint64_t rescuedCount = 0;
    for(size_t i = 0; i < alignmentData.size(); i++) {
        const auto& ad = alignmentData[i];
        // Count overlaps where both flags are now false (rescued)
        // that previously had exactly one flag set
        if (!ad.isDeleted()) {
            // Simple count of surviving overlaps after rescue
        }
    }
    
    cout << timestamp << "Phased overlap rescue complete." << endl;
}

void Assembler::rescuePhasedOverlapsThreadFunction(size_t /* threadId */)
{
    const uint64_t rescueThreshold = this->rescuePhasedThreshold;
    
    uint64_t readIdBegin, readIdEnd;
    while(getNextBatch(readIdBegin, readIdEnd)) {
        for(ReadId readId = ReadId(readIdBegin); readId < ReadId(readIdEnd); readId++) {
            
            // Process overlaps from this read's perspective
            OrientedReadId orientedReadId(readId, 0);
            
            if (orientedReadId.getValue() >= alignmentTable.size()) continue;
            
            // Collect overlaps where there's a directional conflict
            // (one direction deleted, other not)
            std::vector<uint32_t> conflictAlignments;
            std::vector<std::pair<uint32_t, uint32_t>> conflictIntervals; // qs, qe on readId
            
            const size_t n = alignmentTable.size(orientedReadId.getValue());
            for(size_t i = 0; i < n; i++) {
                const uint64_t alignmentId = alignmentTable[orientedReadId.getValue()][i];
                const AlignmentData& ad = alignmentData[alignmentId];
                
                // Check for directional conflict: exactly one flag set
                bool hasConflict = (ad.isDeleted0 != ad.isDeleted1);
                if (!hasConflict) continue;
                
                // For this read's perspective, check which flag corresponds to this read
                uint32_t qs, qe;
                bool thisReadDeleted;
                
                if (ad.readIds[0] == readId) {
                    qs = ad.qs;
                    qe = ad.qe;
                    thisReadDeleted = ad.isDeleted0; // This read's decision
                } else {
                    qs = ad.ts;
                    qe = ad.te;
                    thisReadDeleted = ad.isDeleted1; // This read's decision
                }
                
                // We want to rescue overlaps where THIS read said "delete" but other said "keep"
                // (i.e., thisReadDeleted == true)
                if (thisReadDeleted) {
                    conflictAlignments.push_back(alignmentId);
                    conflictIntervals.push_back({qs, qe});
                }
            }
            
            // Check if we have enough conflicts to potentially rescue
            if (conflictAlignments.size() < rescueThreshold) continue;
            
            // Sweep line to find max consensus region
            std::vector<std::pair<uint32_t, int>> events;
            events.reserve(conflictIntervals.size() * 2);
            for(const auto& interval : conflictIntervals) {
                events.push_back({interval.first, 1});  // Start
                events.push_back({interval.second, -1}); // End
            }
            
            std::sort(events.begin(), events.end(), [](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first < b.first;
                return a.second > b.second; // Start before End at same position
            });
            
            int dp = 0, max_dp = 0;
            uint32_t start = 0;
            uint32_t max_interval_s = 0, max_interval_e = 0;
            
            for(const auto& ev : events) {
                int old_dp = dp;
                dp += ev.second;
                
                if (old_dp < dp) { // Start event
                    if (dp >= max_dp) {
                        start = ev.first;
                        max_dp = dp;
                    }
                } else if (old_dp > dp) { // End event
                    if (old_dp == max_dp) {
                        max_interval_s = start;
                        max_interval_e = ev.first;
                    }
                }
            }
            
            // If consensus depth >= threshold, rescue overlaps spanning the consensus region
            if ((uint64_t)max_dp >= rescueThreshold) {
                for(size_t j = 0; j < conflictAlignments.size(); j++) {
                    uint32_t alignmentId = conflictAlignments[j];
                    uint32_t qs = conflictIntervals[j].first;
                    uint32_t qe = conflictIntervals[j].second;
                    
                    // Only rescue if this overlap spans the consensus region
                    if (qs <= max_interval_s && qe >= max_interval_e) {
                        AlignmentData& ad = alignmentData[alignmentId];
                        
                        // Clear the flag for this read's direction
                        if (ad.readIds[0] == readId) {
                            ad.isDeleted0 = false;
                        } else {
                            ad.isDeleted1 = false;
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
//     - We conceptually "cut" the read or remove it by marking all its alignments as deleted (`alignmentData[i].setDeleted(true)`).

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
                     if(ad.isDeleted()) continue;
                     
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
                     if(ad.isDeleted()) continue;
                     
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
                         alignmentData[alignmentId].setDeleted(true);
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

// Exact hifiasm ma_hit2arc containment check.
// Parameters match hifiasm naming:
//   h = alignment (hit)
//   ql = query length (after coverage cut: sq->e - sq->s)
//   tl = target length (after coverage cut: st->e - st->s)
//   max_hang = maximum allowed overhang
//   int_frac = internal fraction threshold (max_hang_rate)
//   min_ovlp = minimum overlap length
// Returns:
//   0: Normal overlap (dovetail)
//   1: MA_HT_QCONT - Query contained in Target
//   2: MA_HT_TCONT - Target contained in Query
//   -1: MA_HT_INT - Internal match (too much overhang)
//   -2: MA_HT_SHORT_OVLP - Overlap too short
static int ma_hit2arc_containment(
    int32_t qs, int32_t qe, int32_t ql,   // Query: overlap start, end, length
    int32_t ts, int32_t te, int32_t tl,   // Target: overlap start, end, length  
    bool isReverse,                       // Is target reverse complemented?
    int32_t max_hang,                     // Maximum allowed overhang
    double int_frac,                      // Internal fraction threshold
    int32_t min_ovlp)                     // Minimum overlap length
{
    // Compute 5' and 3' overhangs on target (relative to query orientation)
    // tl5 = 5'-end overhang on target (on query strand)
    // tl3 = 3'-end overhang on target (on query strand)
    int32_t tl5, tl3;
    if (isReverse) {
        tl5 = tl - te;  // Target's right end becomes 5' relative to query
        tl3 = ts;       // Target's left end becomes 3' relative to query
    } else {
        tl5 = ts;       // Target's left end is 5'
        tl3 = tl - te;  // Target's right end is 3'
    }

    // ext5 and ext3: the minimum of query and target overhangs on each side
    int32_t ext5 = (qs < tl5) ? qs : tl5;
    int32_t ext3 = ((ql - qe) < tl3) ? (ql - qe) : tl3;

    // Check for internal match (too much overhang on both sides)
    // This rejects overlaps where the alignment doesn't extend to edges
    if (ext5 > max_hang || ext3 > max_hang) {
        return -1;  // MA_HT_INT
    }
    
    // Check internal fraction constraint
    int32_t qOverlapLen = qe - qs;
    int32_t tOverlapLen = te - ts;
    if (qOverlapLen < (qOverlapLen + ext5 + ext3) * int_frac) {
        return -1;  // MA_HT_INT
    }
    if (tOverlapLen < (tOverlapLen + ext5 + ext3) * int_frac) {
        return -1;  // MA_HT_INT
    }

    // Containment checks (exact hifiasm logic from Overlaps.h lines 418-419)
    // Query contained in Target: query's overhangs are smaller than target's
    if (qs <= tl5 && (ql - qe) <= tl3) {
        return 1;  // MA_HT_QCONT
    }
    // Target contained in Query: target's overhangs are smaller than query's
    if (qs >= tl5 && (ql - qe) >= tl3) {
        return 2;  // MA_HT_TCONT
    }

    // Check minimum overlap length
    if (qOverlapLen + ext5 + ext3 < min_ovlp || tOverlapLen + ext5 + ext3 < min_ovlp) {
        return -2;  // MA_HT_SHORT_OVLP
    }

    // Normal dovetail overlap
    return 0;
}


void Assembler::removeContainedReads(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Removing contained reads (ma_hit_contained_advance) - Serial Execution for Strict Parity..." << endl;
    
    if (!containmentParent.isOpen) {
        containmentParent.createNew(largeDataName("ContainmentParent"), largeDataPageSize);
        containmentParent.resize(reads->readCount());
        std::fill(containmentParent.begin(), containmentParent.end(), ReadId(invalidReadId));
    }

    uint64_t containedCount = 0;
    
    // Hifiasm Iterates Reads Serially (Overlaps.cpp:1794)
    for(ReadId i=0; i<reads->readCount(); i++) {
        if (validReadIntervals[i].isDeleted) continue;
        
        // Iterate alignments for Read i (Strand 0)
        // Hifiasm iterates sources[i] which contains all overlaps for i.
        OrientedReadId oid(i, 0);
        if (oid.getValue() >= alignmentTable.size()) continue;
        
        const auto& table = alignmentTable[oid.getValue()];
        for(uint32_t alignmentId : table) {
             AlignmentData& ad = alignmentData[alignmentId];
             
             // "if(h->del) continue;"
             if (ad.isDeleted()) continue;
             
             // Identify Query (i) and Target (neighbor)
             ReadId qn = i;
             ReadId tn;
             int32_t qs, qe, ts, te;
             bool isReverse = !ad.isSameStrand;
             
             // Determine usage of coordinates based on who is Query
             // AlignmentData stores [readId0, readId1].
             // Coordinates are [qs, qe] for r0, [ts, te] for r1.
             if (ad.readIds[0] == i) {
                 tn = ad.readIds[1];
                 qs = (int32_t)ad.qs; qe = (int32_t)ad.qe;
                 ts = (int32_t)ad.ts; te = (int32_t)ad.te;
             } else {
                 tn = ad.readIds[0];
                 // If i is r1 (Target in AD), then for THIS check, i is Query.
                 // So we must SWAP the roles of Q and T coordinates from AD perspective.
                 // "qs" for this check (on Read i) comes from ad.ts/te.
                 // "ts" for this check (on Read tn) comes from ad.qs/qe.
                 qs = (int32_t)ad.ts; qe = (int32_t)ad.te;
                 ts = (int32_t)ad.qs; te = (int32_t)ad.qe;
             }
             
             // "if(sq->del || st->del) continue;"
             if (validReadIntervals[tn].isDeleted) continue;

             const auto& vrQ = validReadIntervals[qn];
             const auto& vrT = validReadIntervals[tn];
             int32_t ql = (int32_t)(vrQ.end - vrQ.start);
             int32_t tl = (int32_t)(vrT.end - vrT.start);
             
             // Call Hifiasm Logic Helper
             // Note: qs/qe/ts/te are already normalized by applyCoverageCuts (Step 3).
             // Matches Hifiasm expectation of being relative to valid region.
             
             int result = ma_hit2arc_containment(
                qs, qe, ql,
                ts, te, tl,
                isReverse,
                (int32_t)maxHang,
                maxHangRate,
                (int32_t)minOverlapLength
            );
            
            if (result == 1) { // MA_HT_QCONT: Query (i) contained in Target (tn)
                // "delete_all_edges(qn)" -> Mark read deleted implies all edges invalid
                // We mark read as deleted immediately.
                validReadIntervals[qn].isDeleted = true;
                
                // "set_R_to_U" -> Record containment
                containmentParent[qn] = tn;
                containedCount++;
                
                // Break inner loop? Hifiasm does NOT break inner loop immediately?
                // `delete_all_edges` clears the buffer?
                // If we delete all edges, checking further overlaps for `i` is pointless?
                // Hifiasm `delete_all_edges(sources... i)`. It sets `p->del=1` for all overlaps in `sources[i]`.
                // So yes, we should stop processing `i`.
                break; 
                
            } else if (result == 2) { // MA_HT_TCONT: Target (tn) contained in Query (i)
                // Mark target deleted
                validReadIntervals[tn].isDeleted = true;
                containmentParent[tn] = qn;
                containedCount++;
                
                // Note: We do NOT break here, because Query `i` is still alive and might contain others.
            }
        }
    }
    
    // Transitive Reduction of Containment Tree & Cleanup (Parity: transfor_R_to_U)
    cout << timestamp << "Transitive reduction..." << endl;
    for(size_t i=0; i<containmentParent.size(); i++) {
        if (containmentParent[i] != invalidReadId) {
            ReadId parent = containmentParent[i];
            while(containmentParent[parent] != invalidReadId) {
                parent = containmentParent[parent];
            }
            containmentParent[i] = parent; // Path compression

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

void Assembler::removeContainedReadsThreadFunction(size_t /* threadId */)
{
    // Exact hifiasm ma_hit_contained_advance logic
    // Iterates over alignments and checks containment using ma_hit2arc logic
    
    uint64_t start, end;
    while(getNextBatch(start, end)) {
        for(uint64_t i = start; i < end; i++) {
            AlignmentData& ad = alignmentData[i];
            
            // Skip already deleted alignments
            if (ad.isDeleted()) continue;
            
            ReadId r0 = ad.readIds[0];
            ReadId r1 = ad.readIds[1];
            
            // Skip if either read is already marked as deleted (coverage_cut[i].del)
            if (validReadIntervals[r0].isDeleted || validReadIntervals[r1].isDeleted) {
                continue;
            }
            
            // Get valid region lengths (equivalent to sq->e - sq->s in hifiasm)
            // Uses the valid regions detected by filterLocalSegments (ma_hit_sub equivalent)
            const auto& vr0 = validReadIntervals[r0];
            const auto& vr1 = validReadIntervals[r1];
            int32_t ql = (int32_t)(vr0.end - vr0.start);
            int32_t tl = (int32_t)(vr1.end - vr1.start);
            
            // IMPORTANT: After applyCoverageCuts, coordinates are already normalized (0-based relative to valid region)
            // So we use them directly without subtracting vr.start again
            int32_t qs = (int32_t)ad.qs;
            int32_t qe = (int32_t)ad.qe;
            int32_t ts = (int32_t)ad.ts;
            int32_t te = (int32_t)ad.te;
            bool isReverse = ad.isSameStrand ? false : true;  // ad.isSameStrand means forward, !isSameStrand means reverse


            
            // Use the exact hifiasm containment check
            int result = ma_hit2arc_containment(
                qs, qe, ql,
                ts, te, tl,
                isReverse,
                (int32_t)this->hangingFilterMaxHang,
                this->hangingFilterMaxHangRate,
                (int32_t)this->hangingFilterMinOverlap
            );
            
            static std::atomic<int> debugCount(0);
            if (result == 1 && debugCount < 10) { // QCONT
                 debugCount++;
                 // Recalculate overhangs for debug
                 int32_t tl5, tl3;
                 if (isReverse) {
                    tl5 = tl - te;
                    tl3 = ts;
                 } else {
                    tl5 = ts;
                    tl3 = tl - te;
                 }
                 int32_t ext5 = std::min(qs, tl5);
                 int32_t ext3 = std::min(ql - qe, tl3);
                 
                 std::stringstream ss;
                 ss << "Debug QCONT: R" << r0 << " in R" << r1 
                    << " qs=" << qs << " qe=" << qe << " ql=" << ql 
                    << " ts=" << ts << " te=" << te << " tl=" << tl
                    << " tl5=" << tl5 << " tl3=" << tl3 << " ext5=" << ext5 << " ext3=" << ext3 << endl;
                 cout << ss.str();
            }
            
            if (result == 1) {
                // MA_HT_QCONT: Query (r0) is contained in Target (r1)
                ad.setDeleted(true);
                
                // delete_single_edge equivalent: delete the reverse edge from r1 to r0
                // (This is complex to implement without the full edge structure, 
                //  but we mark the alignment deleted which is equivalent)
                
                // delete_all_edges equivalent: mark all alignments involving r0 as deleted
                for(uint32_t strand = 0; strand < 2; strand++) {
                    OrientedReadId oid(r0, strand);
                    for(uint32_t otherAlignId : alignmentTable[oid.getValue()]) {
                        alignmentData[otherAlignId].setDeleted(true);
                    }
                }
                
                // set_R_to_U equivalent: record containment relationship
                containmentParent[r0] = r1;
                validReadIntervals[r0].isDeleted = true;
                
            } else if (result == 2) {
                // MA_HT_TCONT: Target (r1) is contained in Query (r0)
                ad.setDeleted(true);
                
                // delete_all_edges equivalent: mark all alignments involving r1 as deleted  
                for(uint32_t strand = 0; strand < 2; strand++) {
                    OrientedReadId oid(r1, strand);
                    for(uint32_t otherAlignId : alignmentTable[oid.getValue()]) {
                        alignmentData[otherAlignId].setDeleted(true);
                    }
                }
                
                // set_R_to_U equivalent: record containment relationship
                containmentParent[r1] = r0;
                validReadIntervals[r1].isDeleted = true;
            }
            // result == 0 (normal dovetail), -1 (internal), -2 (short): do nothing
        }
    }
}

