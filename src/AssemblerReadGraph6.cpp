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

void Assembler::createReadGraphFromEcParityCisOverlaps()
{
    createReadGraphFromEcParityCisOverlaps(std::thread::hardware_concurrency(), /*rebuildDirectedReadGraph*/ false);
}

void Assembler::createReadGraphFromEcParityCisOverlaps(
    uint64_t /*threadCount*/,
    bool rebuildDirectedReadGraph)
{
    cout << timestamp << "createReadGraphFromEcParityCisOverlaps begins." << endl;
    checkAlignmentDataAreOpen();

    const uint64_t alignmentCount = alignmentData.size();

    // Avoid vector<bool> in parallel loops (it is bit-packed and can race).
    vector<uint8_t> keepAlignmentByte(alignmentCount, 0);
    uint64_t keptCount = 0;
    uint64_t filteredCount = 0;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:keptCount, filteredCount)
#endif
    for(uint64_t i = 0; i < alignmentCount; i++) {
        const AlignmentData& ad = alignmentData[i];
        const bool cis0 = ((ad.deleteReasons0 & AlignmentData::DeleteReasonPhase) == 0);
        const bool cis1 = ((ad.deleteReasons1 & AlignmentData::DeleteReasonPhase) == 0);
        const bool keep = cis0 && cis1;
        keepAlignmentByte[i] = keep ? 1 : 0;
        if(keep) {
            ++keptCount;
        } else {
            ++filteredCount;
        }
    }

    cout << timestamp << "Kept " << keptCount << " / " << alignmentCount
         << " alignments (phasing cis only), filtered " << filteredCount << "." << endl;

    vector<bool> keepAlignment(alignmentCount, false);
    for(uint64_t i = 0; i < alignmentCount; ++i) {
        keepAlignment[i] = (keepAlignmentByte[i] != 0);
    }

    rebuildReadGraphUsingSelectedAlignments(std::move(keepAlignment), rebuildDirectedReadGraph);

    cout << timestamp << "createReadGraphFromEcParityCisOverlaps completed." << endl;
}

void Assembler::createReadGraphFromEcParityCisOverlapsCoveringInformativeSites()
{
    createReadGraphFromEcParityCisOverlapsCoveringInformativeSites(
        std::thread::hardware_concurrency(),
        /*rebuildDirectedReadGraph*/ false);
}

void Assembler::createReadGraphFromEcParityCisOverlapsCoveringInformativeSites(
    uint64_t /*threadCount*/,
    bool rebuildDirectedReadGraph)
{
    cout << timestamp << "createReadGraphFromEcParityCisOverlapsCoveringInformativeSites begins." << endl;
    checkAlignmentDataAreOpen();

    const uint64_t alignmentCount = alignmentData.size();
    static constexpr uint32_t minInformativeSiteCount = 2;

    // Avoid vector<bool> in parallel loops (it is bit-packed and can race).
    vector<uint8_t> keepAlignmentByte(alignmentCount, 0);
    uint64_t keptCount = 0;
    uint64_t filteredByPhase = 0;
    uint64_t filteredByNoInformativeSite = 0;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:keptCount, filteredByPhase, filteredByNoInformativeSite)
#endif
    for(uint64_t i = 0; i < alignmentCount; i++) {
        const AlignmentData& ad = alignmentData[i];
        const bool cis0 = ((ad.deleteReasons0 & AlignmentData::DeleteReasonPhase) == 0);
        const bool cis1 = ((ad.deleteReasons1 & AlignmentData::DeleteReasonPhase) == 0);
        if(!(cis0 && cis1)) {
            keepAlignmentByte[i] = 0;
            ++filteredByPhase;
            continue;
        }
        if(!ad.coversHetSiteAtLeast(minInformativeSiteCount)) {
            keepAlignmentByte[i] = 0;
            ++filteredByNoInformativeSite;
            continue;
        }
        keepAlignmentByte[i] = 1;
        ++keptCount;
    }

    cout << timestamp << "Kept " << keptCount << " / " << alignmentCount
         << " alignments (cis + covers >= " << minInformativeSiteCount << " informative sites)." << endl;
    cout << timestamp << "Filtered: phase=" << filteredByPhase
         << ", noInformativeSite=" << filteredByNoInformativeSite << endl;

    vector<bool> keepAlignment(alignmentCount, false);
    for(uint64_t i = 0; i < alignmentCount; ++i) {
        keepAlignment[i] = (keepAlignmentByte[i] != 0);
    }

    rebuildReadGraphUsingSelectedAlignments(std::move(keepAlignment), rebuildDirectedReadGraph);

    cout << timestamp << "createReadGraphFromEcParityCisOverlapsCoveringInformativeSites completed." << endl;
}

void Assembler::createReadGraph6(uint64_t threadCount)
{
    cout << timestamp << "createReadGraph6 begins." << endl;

    // Preconditions: overlap/alignment records must be loaded.
    // This routine applies overlap-level and read-level filters, then rebuilds the read graph
    // from the surviving overlaps.
    checkAlignmentDataAreOpen();
    // checkPhasingCigarsAreOpen(); // Removed: PhasingCigars replaced by AlignedEvidenceStore

    const uint64_t totalAlignments = alignmentData.size();

    // Count per-side overlap deletion flags. These are overlap-local diagnostics and are
    // different from read-level deletion (validReadIntervals[r].isDeleted).
    auto countPhasingFlags = [this]() -> std::pair<uint64_t, uint64_t> {
        uint64_t del0 = 0, del1 = 0;
        for(uint64_t i = 0; i < alignmentData.size(); i++) {
            if(alignmentData[i].isDeleted0()) del0++;
            if(alignmentData[i].isDeleted1()) del1++;
        }
        return {del0, del1};
    };

    // Count reads currently disabled at the read level (ma_hit_sub/chimera/contained pipeline state).
    auto countDeletedReads = [this]() -> uint64_t {
        if(validReadIntervals.empty()) {
            return 0;
        }
        uint64_t deleted = 0;
        for(const auto& interval : validReadIntervals) {
            if(interval.isDeleted) {
                ++deleted;
            }
        }
        return deleted;
    };
    
    cout << timestamp << "[DIAG] Total alignments: " << totalAlignments << endl;
    
    // Initial overlap-level state before read-level filtering.
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
    
    // Step 1: Find valid read segments (ma_hit_sub equivalent).
    // Parity note:
    // - hifiasm default min_overlap_coverage is 0
    // - ma_hit_sub treats min_dp <= 1 as "full read is valid"
    // So with minCoverage=0 we keep [0, readLen) for every read, but still populate
    // validReadIntervals, which downstream stages use as canonical read-level state.
    const uint64_t minCoverage = 0;
    filterLocalSegments(minCoverage, threadCount);

    const uint64_t deletedReads = countDeletedReads();
    cout << timestamp << "[DIAG] After filterLocalSegments: " << deletedReads << "/" << reads->readCount() << " reads marked deleted" << endl;
    cout << timestamp << "[DIAG] After filterLocalSegments: active alignments=" << countActiveAlignments() << endl;

    // Step 2: Detect chimeric reads (detect_chimeric_reads)
    // Runs after ma_hit_sub and before ma_hit_cut in hifiasm. This stage marks chimeric
    // reads at read-level state (validReadIntervals/isChimericRead) and propagates overlap
    // deletion reasons so downstream graph construction can prune incident edges.
    detectChimericReads(threadCount);
    auto [afterChimericDel0, afterChimericDel1] = countPhasingFlags();
    cout << timestamp << "[DIAG] After detectChimericReads: isDeleted0=" << afterChimericDel0
         << ", isDeleted1=" << afterChimericDel1
         << ", active=" << countActiveAlignments() << endl;
    cout << timestamp << "[DIAG] After detectChimericReads: deletedReads="
         << countDeletedReads() << "/" << reads->readCount() << endl;

    // Step 3: Clip overlaps to valid regions (ma_hit_cut equivalent)
    // Also normalizes coordinates to 0-based relative to valid region.
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


    





    // // Step 6a: Contained read detection (does not remove overlaps).
    // flagContainedReads(maxHang, maxHangRate, minOverlapLength, threadCount);
    // uint64_t containedFlagCount = 0;
    // for (ReadId r = 0; r < reads->readCount(); ++r) {
    //     if (reads->getFlags(r).isContained) {
    //         ++containedFlagCount;
    //     }
    // }
    // cout << timestamp << "[DIAG] After flagContainedReads: containedReads=" << containedFlagCount << endl;

    // // Step 6c: For each contained read, keep only one best overlap (by dpScore) and prune all others.
    // // This is a diagnostic/experimental alternative to removing contained reads entirely.
    // pruneContainedReadsToOneBestOverlapByDpScore(threadCount);





    
    // Step 7: Apply final AND-semantics filter
    // By this point, all filtering steps (2-6) have set delete reasons on individual alignment sides.
    // This step applies the conservative AND rule: keep an alignment ONLY if both sides keep it.
    const uint64_t alignmentCount = alignmentData.size();
    std::vector<bool> keepAlignment(alignmentCount, true);

    uint64_t keptCount = 0;
    uint64_t filteredCount = 0;

    // Breakdown of filter reasons for diagnostics (alignments can have multiple reasons)
    uint64_t filteredByPhase = 0;
    uint64_t filteredBySecondary = 0;
    uint64_t filteredByChemical = 0;
    uint64_t filteredByLocalSegment = 0;
    uint64_t filteredByCoverageCut = 0;
    uint64_t filteredByHanging = 0;
    uint64_t filteredByContained = 0;

    for(uint64_t i = 0; i < alignmentCount; i++) {
        auto& ad = alignmentData[i];

        // AND semantics: keep only if BOTH sides keep it (no delete reasons on either side)
        if(!ad.keptByBothSides()) {
            keepAlignment[i] = false;
            ad.info.isInReadGraph = 0;
            filteredCount++;

            // Classify the reason(s) - alignments can have multiple delete reasons
            const AlignmentData::DeleteReasonMask reasons = ad.deleteReasons0 | ad.deleteReasons1;
            if (reasons & AlignmentData::DeleteReasonPhase) ++filteredByPhase;
            if (reasons & AlignmentData::DeleteReasonSecondary) ++filteredBySecondary;
            if (reasons & AlignmentData::DeleteReasonChemical) ++filteredByChemical;
            if (reasons & AlignmentData::DeleteReasonLocal) ++filteredByLocalSegment;
            if (reasons & AlignmentData::DeleteReasonCoverageCut) ++filteredByCoverageCut;
            if (reasons & AlignmentData::DeleteReasonHanging) ++filteredByHanging;
            if (reasons & AlignmentData::DeleteReasonContained) ++filteredByContained;
            continue;
        }

        // Alignment passes all filters - mark as in read graph
        ad.info.isInReadGraph = 1;
        keptCount++;
    }

    cout << timestamp << "Read graph filtering complete: " << filteredCount << " alignments filtered, "
         << keptCount << " kept (" << (100.0 * keptCount / alignmentCount) << "%)." << endl;
    if (filteredCount > 0) {
        cout << timestamp << "  Filter reason breakdown (alignments can have multiple reasons):" << endl;
        cout << timestamp << "    phase=" << filteredByPhase
             << " secondary=" << filteredBySecondary
             << " chemical=" << filteredByChemical
             << " local(ma_hit_sub)=" << filteredByLocalSegment
             << " cut(ma_hit_cut)=" << filteredByCoverageCut
             << " hanging(ma_hit_flt)=" << filteredByHanging
             << " contained(ma_hit_contained)=" << filteredByContained
             << endl;
    }

    // Step 8: Create read graph from kept alignments.
    // Directed read graph is optional for createReadGraph6 (it is expensive and not
    // required by all pipelines). Enable if needed.
    rebuildReadGraphUsingSelectedAlignments(keepAlignment, /*rebuildDirectedReadGraph*/false);

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

#ifdef _OPENMP
    #pragma omp parallel for reduction(+:keptCount, filteredCount)
#endif
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
        if (!validReadIntervals.empty()) {
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
    rebuildReadGraphUsingSelectedAlignments(keepAlignment, /*rebuildDirectedReadGraph*/true);
    cout << timestamp << "createReadGraphFromFilteredAlignments completed." << endl;
}
