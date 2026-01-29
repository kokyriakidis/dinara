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
    createReadGraph6(std::thread::hardware_concurrency());
}

void Assembler::createReadGraph6(uint64_t threadCount)
{
    cout << timestamp << "createReadGraph6 begins." << endl;

    // Check required data
    checkAlignmentDataAreOpen();
    // checkPhasingCigarsAreOpen(); // Removed: PhasingCigars replaced by AlignedEvidenceStore

    const uint64_t totalAlignments = alignmentData.size();
    
    // Helper to count active alignments (not deleted)
    auto countActiveAlignments = [this]() -> uint64_t {
        uint64_t active = 0;
        for(uint64_t i = 0; i < alignmentData.size(); i++) {
            if(alignmentData[i].keptByBothSides()) active++;
        }
        return active;
    };
    
    // Helper to count reads with isDeleted0 or isDeleted1 set
    auto countPhasingFlags = [this]() -> std::pair<uint64_t, uint64_t> {
        uint64_t del0 = 0, del1 = 0;
        for(uint64_t i = 0; i < alignmentData.size(); i++) {
            if(alignmentData[i].isDeleted0()) del0++;
            if(alignmentData[i].isDeleted1()) del1++;
        }
        return {del0, del1};
    };
    
    cout << timestamp << "[DIAG] Total alignments: " << totalAlignments << endl;
    
    // Check initial state
    auto [initDel0, initDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] Initial state: isDeleted0=" << initDel0 
         << ", isDeleted1=" << initDel1 
         << ", active=" << countActiveAlignments() << endl;
    
    // ONT-only parity step: chemical arc masking runs before overlap rescue and clean_graph.
    // It uses all overlaps (independent of EC/phasing) and annotates affected overlaps with
    // DeleteReasonChemical so later steps can avoid them.
    applyOntChemicalArcMask(threadCount);

    // Hifiasm parity: Full pre-graph filtering pipeline
    // Order matches Hifiasm's clean_graph / gen_init_sg flow
    
    // Step 1: Rescue phased overlaps with directional conflicts (try_rescue_overlaps equivalent)
    rescuePhasedOverlaps(4, threadCount);
    auto [afterRescueDel0, afterRescueDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] After rescuePhasedOverlaps: isDeleted0=" << afterRescueDel0 
         << ", isDeleted1=" << afterRescueDel1 
         << ", active=" << countActiveAlignments() << endl;
    
    // Step 2: Find valid read segments (ma_hit_sub equivalent)
    // Hifiasm default min_dp is asm_opt.min_overlap_coverage (default 0).
    const uint64_t minCoverage = 0;
    filterLocalSegments(minCoverage, threadCount);
    
    // Count deleted reads
    uint64_t deletedReads = 0;
    for(uint64_t r = 0; r < validReadIntervals.size(); r++) {
        if(validReadIntervals[r].isDeleted) deletedReads++;
    }
    cout << timestamp << "[DIAG] After filterLocalSegments: " << deletedReads << "/" << reads->readCount() << " reads marked deleted" << endl;
    cout << timestamp << "[DIAG] After filterLocalSegments: active alignments=" << countActiveAlignments() << endl;
    
    // Step 3: Detect chimeric reads (detect_chimeric_reads)
    // Runs after ma_hit_sub and before ma_hit_cut in hifiasm.
    detectChimericReads(threadCount);
    auto [afterChimericDel0, afterChimericDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] After detectChimericReads: isDeleted0=" << afterChimericDel0
         << ", isDeleted1=" << afterChimericDel1
         << ", active=" << countActiveAlignments() << endl;

    // Step 4: Clip overlaps to valid regions (ma_hit_cut equivalent)
    // Also normalizes coordinates to 0-based relative to valid region.
    const uint64_t minOverlapLength = 50;  // Hifiasm default mini_overlap_length
    applyCoverageCuts(minOverlapLength, threadCount);
    auto [afterCutDel0, afterCutDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] After applyCoverageCuts: isDeleted0=" << afterCutDel0 
         << ", isDeleted1=" << afterCutDel1 
         << ", active=" << countActiveAlignments() << endl;
    
    // Step 5: Filter hanging overlaps (ma_hit_flt equivalent)
    // Removes overlaps with excessive overhangs
    const uint64_t maxHang = 1000;
    const double maxHangRate = 0.8;
    filterHangingOverlaps(maxHang, maxHangRate, minOverlapLength, threadCount);
    auto [afterHangDel0, afterHangDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] After filterHangingOverlaps: isDeleted0=" << afterHangDel0 
         << ", isDeleted1=" << afterHangDel1 
         << ", active=" << countActiveAlignments() << endl;
    
    // // Step 6: Remove contained reads (ma_hit_contained_advance equivalent)
    // // Marks fully contained reads and removes their overlaps
    // removeContainedReads(maxHang, maxHangRate, minOverlapLength, threadCount);
    // auto [afterContainDel0, afterContainDel1] = countPhasingFlags();
    // cout << timestamp << "[DIAG] After removeContainedReads: isDeleted0=" << afterContainDel0
    //      << ", isDeleted1=" << afterContainDel1
    //      << ", active=" << countActiveAlignments() << endl;

    // Step 6b: Diagnostic-only contained read detection (does not remove overlaps).
    // Useful when experimenting with phasing/anchors without changing the overlap set.
    flagContainedReads(maxHang, maxHangRate, minOverlapLength, threadCount);
    uint64_t containedFlagCount = 0;
    for (ReadId r = 0; r < reads->readCount(); ++r) {
        if (reads->getFlags(r).isContained) {
            ++containedFlagCount;
        }
    }
    cout << timestamp << "[DIAG] After flagContainedReads: containedReads=" << containedFlagCount << endl;

    // Step 6c: For each contained read, keep only one best overlap (by dpScore) and prune all others.
    // This is a diagnostic/experimental alternative to removing contained reads entirely.
    pruneContainedReadsToOneBestOverlapByDpScore(threadCount);





    
    // Step 7: Final filtering pass - apply phasing decisions
    const uint64_t alignmentCount = alignmentData.size();
    std::vector<bool> keepAlignment(alignmentCount, true);
    
    uint64_t keptCount = 0;
    uint64_t phasedOutCount = 0;
    uint64_t filteredByPhase = 0;
    uint64_t filteredBySecondary = 0;
    uint64_t filteredByChemical = 0;
    uint64_t filteredByLocalSegment = 0;
    uint64_t filteredByCoverageCut = 0;
    uint64_t filteredByHanging = 0;
    uint64_t filteredByContained = 0;
    uint64_t palindromicCount = 0;
    uint64_t chimericCount = 0;
    uint64_t containedCount = 0;

    auto classifyFilterReason = [&](const AlignmentData& ad) {
        const AlignmentData::DeleteReasonMask reasons = ad.deleteReasons0 | ad.deleteReasons1;
        if (reasons & AlignmentData::DeleteReasonPhase) ++filteredByPhase;
        if (reasons & AlignmentData::DeleteReasonSecondary) ++filteredBySecondary;
        if (reasons & AlignmentData::DeleteReasonChemical) ++filteredByChemical;
        if (reasons & AlignmentData::DeleteReasonLocal) ++filteredByLocalSegment;
        if (reasons & AlignmentData::DeleteReasonCoverageCut) ++filteredByCoverageCut;
        if (reasons & AlignmentData::DeleteReasonHanging) ++filteredByHanging;
        if (reasons & AlignmentData::DeleteReasonContained) ++filteredByContained;
    };

    #pragma omp parallel for reduction(+:keptCount, phasedOutCount, filteredByPhase, filteredBySecondary, filteredByChemical, filteredByLocalSegment, filteredByCoverageCut, filteredByHanging, filteredByContained, palindromicCount, chimericCount, containedCount)
    for(uint64_t i = 0; i < alignmentCount; i++) {
        auto& ad = alignmentData[i];
        
        // 1. Conservative AND: keep an overlap only if BOTH reads keep it.
        if(!ad.keptByBothSides()) {
            keepAlignment[i] = false;
            ad.info.isInReadGraph = 0;
            phasedOutCount++;
            classifyFilterReason(ad);
            continue;
        }
        
        // 2. Check if reads are marked as deleted in validReadIntervals
        if (validReadIntervals.size() > 0) {
            if (validReadIntervals[ad.readIds[0]].isDeleted || 
                validReadIntervals[ad.readIds[1]].isDeleted) {
                const AlignmentData::DeleteReasonMask already =
                    ad.deleteReasons0 | ad.deleteReasons1;
                const AlignmentData::DeleteReasonMask majorReadLevel =
                    AlignmentData::DeleteReasonContained | AlignmentData::DeleteReasonChimeric;
                if ((already & majorReadLevel) == 0) {
                    ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonLocal);
                }
                keepAlignment[i] = false;
                ad.info.isInReadGraph = 0;
                containedCount++;
                continue;
            }
        }
        
        // 3. Check palindromic reads
        if (reads->getFlags(ad.readIds[0]).isPalindromic || 
            reads->getFlags(ad.readIds[1]).isPalindromic) {
            ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonPalindromic);
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
    if (phasedOutCount) {
        cout << timestamp << "  reasons among not-kept overlaps:"
             << " phase=" << filteredByPhase
             << " secondary=" << filteredBySecondary
             << " chemical=" << filteredByChemical
             << " ma_hit_sub(read)=" << filteredByLocalSegment
             << " ma_hit_cut=" << filteredByCoverageCut
             << " ma_hit_flt=" << filteredByHanging
             << " ma_hit_contained=" << filteredByContained
             << endl;
    }
    cout << timestamp << "Contained/deleted reads removed " << containedCount << " alignments." << endl;
    cout << timestamp << "Palindromic filter removed " << palindromicCount << " alignments." << endl;
    cout << timestamp << "Chimeric filter removed " << chimericCount << " alignments." << endl;
    cout << timestamp << "Kept " << keptCount << " / " << alignmentCount << " alignments." << endl;

    // Step 8: Create read graph from kept alignments
    createReadGraphUsingSelectedAlignments(keepAlignment);
    // createDirectedReadGraphUsingSelectedAlignments(keepAlignment);

    // // Step 9: Create directed string graph arcs (hifiasm-style suffix->prefix)
    // createStringGraphUsingSelectedAlignments(keepAlignment);
    //
    // // Step 10: Initial string-graph cleaning matching hifiasm defaults.
    // // Hifiasm runs `asg_arc_del_trans(gap_fuzz)` inside `gen_init_sg`, then cuts short tips.
    // // Defaults: gap_fuzz=1000, max_short_tip=3 (tip unitigs of <=3 reads).
    // cleanStringGraphInitialHifiasm(/*gapFuzz*/1000, /*maxShortTipReads*/3);
    //
    // // Step 11: Next-stage topology clean (hifiasm pre_clean intent): remove very small bubbles + re-cut tips.
    // cleanStringGraphPreCleanHifiasm(/*maxShortTipReads*/3);
    //
    // // Step 12: Additional cleanup rounds (approximate hifiasm clean_graph rounds):
    // // progressively drop short overlaps by overlap-length ratio, with intermittent pre-clean.
    // // Hifiasm defaults: clean_round=4, min_drop_rate=0.2, max_drop_rate=0.8, finalMinOverlapLen=2000.
    // cleanStringGraphDropOverlapRoundsHifiasm(
    //     /*cleanRounds*/4,
    //     /*minDropRate*/0.2,
    //     /*maxDropRate*/0.8,
    //     /*maxShortTipReads*/3,
    //     /*finalMinOverlapLen*/2000);
    //
    // // Next hifiasm stage: compress the cleaned string graph into unitigs.
    // createUnitigGraphFromStringGraph();
    //
    // // Additional hifiasm-like cleaning at the unitig level (topology-only).
    // cleanUnitigGraphInitialHifiasm(/*gapFuzz*/1000, /*maxShortTipUnitigs*/3);
    // cleanUnitigGraphPreCleanHifiasm(/*maxShortTipUnitigs*/3);
    // cleanUnitigGraphDropOverlapRoundsHifiasm(
    //     /*cleanRounds*/4,
    //     /*minDropRate*/0.2,
    //     /*maxDropRate*/0.8,
    //     /*maxShortTipUnitigs*/3,
    //     /*finalMinOverlapLen*/2000);
    
    cout << timestamp << "createReadGraph6 completed." << endl;
}

void Assembler::createReadGraphFromFilteredAlignments()
{
    cout << timestamp << "createReadGraphFromFilteredAlignments begins." << endl;
    checkAlignmentDataAreOpen();

    const uint64_t alignmentCount = alignmentData.size();
    std::vector<uint8_t> keepAlignmentByte(alignmentCount, 0);
    uint64_t keptCount = 0;
    uint64_t filteredCount = 0;

    #pragma omp parallel for reduction(+:keptCount, filteredCount)
    for(uint64_t i = 0; i < alignmentCount; i++) {
        AlignmentData& ad = alignmentData[i];
        // Hifiasm-style conservative keep rule: keep an overlap only if BOTH reads keep it.
        const bool keptByBothSides = ad.keptByBothSides();
        if (!keptByBothSides) {
            keepAlignmentByte[i] = 0;
            ad.info.isInReadGraph = 0;
            filteredCount++;
            continue;
        }

        // Also check if reads are marked deleted globally (e.g. from chimeras/containment).
        if (validReadIntervals.size() > 0) {
            if (validReadIntervals[ad.readIds[0]].isDeleted ||
                validReadIntervals[ad.readIds[1]].isDeleted) {
                keepAlignmentByte[i] = 0;
                ad.info.isInReadGraph = 0;
                filteredCount++;
                continue;
            }
        }

        keepAlignmentByte[i] = 1;
        ad.info.isInReadGraph = 1;
        keptCount++;
    }

    cout << timestamp << "Filtered out " << filteredCount << " alignments." << endl;
    cout << timestamp << "Kept " << keptCount << " / " << alignmentCount << " alignments for read graph." << endl;

    std::vector<bool> keepAlignment(alignmentCount, false);
    for (uint64_t i = 0; i < alignmentCount; ++i) {
        keepAlignment[i] = (keepAlignmentByte[i] != 0);
    }
    createReadGraphUsingSelectedAlignments(keepAlignment);
    createDirectedReadGraphUsingSelectedAlignments(keepAlignment);
    cout << timestamp << "createReadGraphFromFilteredAlignments completed." << endl;
}
