// AssemblerReadGraph6.cpp
// Read graph creation using CIGAR-based phasing (Hifiasm parity)

#include "Assembler.hpp"
#include "timestamp.hpp"
#include "Reads.hpp"
#include <vector>
#include <thread>

using namespace dinara;

void Assembler::createReadGraph6()
{
    cout << timestamp << "createReadGraph6 begins." << endl;

    // Check required data
    checkAlignmentDataAreOpen();
    checkPhasingCigarsAreOpen();

    const uint64_t threadCount = std::thread::hardware_concurrency();
    const uint64_t totalAlignments = alignmentData.size();
    
    // Helper to count active alignments (not deleted)
    auto countActiveAlignments = [this]() -> uint64_t {
        uint64_t active = 0;
        for(uint64_t i = 0; i < alignmentData.size(); i++) {
            if(!alignmentData[i].isDeleted()) active++;
        }
        return active;
    };
    
    // Helper to count reads with isDeleted0 or isDeleted1 set
    auto countPhasingFlags = [this]() -> std::pair<uint64_t, uint64_t> {
        uint64_t del0 = 0, del1 = 0;
        for(uint64_t i = 0; i < alignmentData.size(); i++) {
            if(alignmentData[i].isDeleted0) del0++;
            if(alignmentData[i].isDeleted1) del1++;
        }
        return {del0, del1};
    };
    
    cout << timestamp << "[DIAG] Total alignments: " << totalAlignments << endl;
    
    // Check initial state
    auto [initDel0, initDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] Initial state: isDeleted0=" << initDel0 
         << ", isDeleted1=" << initDel1 
         << ", active=" << countActiveAlignments() << endl;
    
    // Hifiasm parity: Full pre-graph filtering pipeline
    // Order matches Hifiasm's clean_graph / gen_init_sg flow
    
    // Step 1: Rescue phased overlaps with directional conflicts (try_rescue_overlaps equivalent)
    rescuePhasedOverlaps(4, threadCount);
    auto [afterRescueDel0, afterRescueDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] After rescuePhasedOverlaps: isDeleted0=" << afterRescueDel0 
         << ", isDeleted1=" << afterRescueDel1 
         << ", active=" << countActiveAlignments() << endl;
    
    // Step 2: Find valid read segments (ma_hit_sub equivalent)
    // Finds longest contiguous region with coverage >= minCoverage
    const uint64_t minCoverage = 3;  // Hifiasm default min_dp
    filterLocalSegments(minCoverage, threadCount);
    
    // Count deleted reads
    uint64_t deletedReads = 0;
    for(uint64_t r = 0; r < validReadIntervals.size(); r++) {
        if(validReadIntervals[r].isDeleted) deletedReads++;
    }
    cout << timestamp << "[DIAG] After filterLocalSegments: " << deletedReads << "/" << reads->readCount() << " reads marked deleted" << endl;
    cout << timestamp << "[DIAG] After filterLocalSegments: active alignments=" << countActiveAlignments() << endl;
    
    // Step 3: Clip overlaps to valid regions (ma_hit_cut equivalent)  
    // Also normalizes coordinates to 0-based relative to valid region
    const uint64_t minOverlapLength = 50;  // Hifiasm default mini_overlap_length
    applyCoverageCuts(minOverlapLength, threadCount);
    auto [afterCutDel0, afterCutDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] After applyCoverageCuts: isDeleted0=" << afterCutDel0 
         << ", isDeleted1=" << afterCutDel1 
         << ", active=" << countActiveAlignments() << endl;
    
    // Step 4: Filter hanging overlaps (ma_hit_flt equivalent)
    // Removes overlaps with excessive overhangs
    const uint64_t maxHang = 1000;
    const double maxHangRate = 0.8;
    filterHangingOverlaps(maxHang, maxHangRate, minOverlapLength, threadCount);
    auto [afterHangDel0, afterHangDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] After filterHangingOverlaps: isDeleted0=" << afterHangDel0 
         << ", isDeleted1=" << afterHangDel1 
         << ", active=" << countActiveAlignments() << endl;
    
    // Step 5: Remove contained reads (ma_hit_contained_advance equivalent)
    // Marks fully contained reads and removes their overlaps
    removeContainedReads(maxHang, maxHangRate, minOverlapLength, threadCount);
    auto [afterContainDel0, afterContainDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] After removeContainedReads: isDeleted0=" << afterContainDel0 
         << ", isDeleted1=" << afterContainDel1 
         << ", active=" << countActiveAlignments() << endl;
    
    // Step 6: Chimeric read detection (gen_chemical_arc_rf equivalent)
    // Detects reads with coverage gaps
    detectChimericReads(threadCount);
    rescueChimericReads(threadCount);
    auto [afterChimericDel0, afterChimericDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] After chimeric detection: isDeleted0=" << afterChimericDel0 
         << ", isDeleted1=" << afterChimericDel1 
         << ", active=" << countActiveAlignments() << endl;

    // Step 7: Final filtering pass - apply phasing decisions
    const uint64_t alignmentCount = alignmentData.size();
    std::vector<bool> keepAlignment(alignmentCount, true);
    
    uint64_t keptCount = 0;
    uint64_t phasedOutCount = 0;
    uint64_t palindromicCount = 0;
    uint64_t chimericCount = 0;
    uint64_t containedCount = 0;

    #pragma omp parallel for reduction(+:keptCount, phasedOutCount, palindromicCount, chimericCount, containedCount)
    for(uint64_t i = 0; i < alignmentCount; i++) {
        auto& ad = alignmentData[i];
        
        // 1. Check dual phasing flags (conservative AND: keep only if BOTH agree)
        if(ad.isDeleted()) {
            keepAlignment[i] = false;
            ad.info.isInReadGraph = 0;
            phasedOutCount++;
            continue;
        }
        
        // 2. Check if reads are marked as deleted in validReadIntervals
        if (validReadIntervals.size() > 0) {
            if (validReadIntervals[ad.readIds[0]].isDeleted || 
                validReadIntervals[ad.readIds[1]].isDeleted) {
                keepAlignment[i] = false;
                ad.info.isInReadGraph = 0;
                containedCount++;
                continue;
            }
        }
        
        // 3. Check palindromic reads
        if (reads->getFlags(ad.readIds[0]).isPalindromic || 
            reads->getFlags(ad.readIds[1]).isPalindromic) {
            keepAlignment[i] = false;
            ad.info.isInReadGraph = 0;
            palindromicCount++;
            continue;
        }
        
        // 4. Check chimeric reads
        if(isChimericRead.size() > 0) {
            if(isChimericRead[ad.readIds[0]] || isChimericRead[ad.readIds[1]]) {
                keepAlignment[i] = false;
                ad.info.isInReadGraph = 0;
                chimericCount++;
                continue;
            }
        }
        
        // Alignment passes all filters
        ad.info.isInReadGraph = 1;
        keptCount++;
    }

    cout << timestamp << "Phasing removed " << phasedOutCount << " alignments." << endl;
    cout << timestamp << "Contained/deleted reads removed " << containedCount << " alignments." << endl;
    cout << timestamp << "Palindromic filter removed " << palindromicCount << " alignments." << endl;
    cout << timestamp << "Chimeric filter removed " << chimericCount << " alignments." << endl;
    cout << timestamp << "Kept " << keptCount << " / " << alignmentCount << " alignments." << endl;

    // Step 8: Create read graph from kept alignments
    createReadGraphUsingSelectedAlignments(keepAlignment);
    
    cout << timestamp << "createReadGraph6 completed." << endl;
}
