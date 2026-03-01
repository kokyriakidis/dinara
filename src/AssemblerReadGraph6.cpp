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

        // Propagate read-level deletions (ma_hit_sub/isolated) to the alignment reasons.
        // Hifiasm ma_sg_gen parity: if(sq->del || st->del) continue;
        if (!validReadIntervals.empty()) {
            if (validReadIntervals[ad.readIds[0]].isDeleted ||
                validReadIntervals[ad.readIds[1]].isDeleted) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonLocal);
            }
        }

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

    // Step 9 (hifiasm `ma_sg_gen`): create directed string graph arcs from kept overlaps.
    createStringGraphUsingSelectedAlignments(keepAlignment);

    stringGraph.writeGfa("ReadGraph6-StringGraph.gfa", *reads);

    // Step 9a: make readGraph match the just-built stringGraph immediately.
    rebuildReadGraphFromCurrentStringGraph(/*rebuildDirectedReadGraph*/false);

    // Step 10 (hifiasm `asg_arc_del_trans`): pure transitive reduction (no tip cutting here).
    reduceStringGraphTransitiveHifiasm(/*gapFuzz*/1000);

    stringGraph.writeGfa("ReadGraph6-StringGraph-TransitiveReduction.gfa", *reads);

    // Step 10a: Reflect string-graph deletions in the read graph so changes are visible there.
    rebuildReadGraphFromCurrentStringGraph(/*rebuildDirectedReadGraph*/false);

    // // Step 10b (hifiasm `ul_clean_gfa` ONT): cut weak arcs before the initial tip cut.
    // // Hifiasm does this only when `asm_opt.is_ont`.
    // // Note: Dinara currently has no explicit ONT/HiFi switch for createReadGraph6, so this is always run.
    // // // TODO(hifiasm-parity): gate this behind an explicit ONT mode flag.
    // // cutStringGraphWeakArcsOntHifiasm(/*maxExtReads*/3, /*lenRatio*/0.975, /*minDiff*/16);

    // Step 11 (hifiasm `ul_clean_gfa` initial): cut short tips before main cleaning loop.
    // Hifiasm calls `asg_arc_cut_tips(sg, max_tip, ...)` as first operation in ul_clean_gfa.
    // Default max_short_tip=3 means unitigs with ≤3 reads are considered tips.
    // TODO(hifiasm-parity): hifiasm's `asg_arc_cut_tips(..., telo_end_t* te)` skips cutting if the tip start
    // or any vertex along the candidate tip-unitig is telomere-marked (`te->hh[readId]`), to avoid trimming
    // true chromosome ends. Dinara's `StringGraph` currently has only `readDeleted` (no telomere state), so
    // `cutStringGraphTips()` cannot implement this "don't cut if telomeric" rule yet.
    // TODO(hifiasm-parity): UL/ONT mode in hifiasm uses `asg_arc_cut_tips(..., is_ou, R_to_U* ru, ...)` and
    // consults per-arc `asg_arc_t.ou` along the candidate path (tracking `mm_ou = MIN(ou)` and adjusting the
    // effective extension length before deciding if it is "short"), plus an additional deletion pass using
    // `R_to_U` when `ru && is_ou`. Dinara's `StringGraphArc` does not store `ou` and the pipeline does not
    // provide `is_ou/ru`, so `cutStringGraphTips()` currently matches only the non-OU/non-ru logic.
    cutStringGraphTips(/*maxShortTipReads*/3);

    stringGraph.writeGfa("ReadGraph6-StringGraph-CutTips.gfa", *reads);

    // Step 11a: Sync read graph after tip cutting.
    rebuildReadGraphFromCurrentStringGraph(/*rebuildDirectedReadGraph*/false);

    // Step 12 (hifiasm `ul_clean_gfa` main loop): iterative cleaning with 4 rounds.
    // Run each sub-step explicitly so we can observe its effect in isolation.
    // Round drop ratios: 0.2 → 0.4 → 0.6 → 0.8 (hifiasm defaults).
    {
        const uint32_t cleanRounds = 0;
        const double minDropRate = 0.2;
        const double maxDropRate = 0.8;
        const uint32_t maxShortTipReads = 3;
        const uint32_t maxCycleReads = 100;

        const double step = (cleanRounds == 1) ? 0. : (maxDropRate - minDropRate) / double(cleanRounds - 1);
        double dropRatio = minDropRate;

        for (uint32_t round = 0; round < cleanRounds; ++round, dropRatio += step) {
            if (dropRatio > maxDropRate) dropRatio = maxDropRate;
            const uint32_t dropPct = uint32_t(dropRatio * 100.0 + 0.5);
            const string roundPrefix = "ReadGraph6-StringGraph-Round" + to_string(round + 1) + "-";

            // asg_iterative_semi_circ (hifiasm ul_clean_gfa line 3078).
            cleanStringGraphBreakShortCycles(maxCycleReads, maxShortTipReads);
            stringGraph.writeGfa(roundPrefix + "1BreakCycles.gfa", *reads);

            // asg_arc_cut_chimeric (near-NOP for clean HiFi; el=1 for all arcs).
            cleanStringGraphChimericReads();
            stringGraph.writeGfa(roundPrefix + "2Chimeric.gfa", *reads);

            // asg_arc_cut_tips (hifiasm ul_clean_gfa line 3082).
            cutStringGraphTips(maxShortTipReads);
            stringGraph.writeGfa(roundPrefix + "3CutTips.gfa", *reads);

            // asg_arc_cut_inexact (near-NOP for clean HiFi; el=1 for all arcs).
            cleanStringGraphInexactOverlaps(maxShortTipReads, /*minDiff*/0);
            stringGraph.writeGfa(roundPrefix + "4InexactOverlaps.gfa", *reads);

            // asg_arc_cut_tips (hifiasm ul_clean_gfa line 3087).
            cutStringGraphTips(maxShortTipReads);
            stringGraph.writeGfa(roundPrefix + "5CutTips.gfa", *reads);

            // asg_arc_cut_length (drop short overlaps by ratio).
            cleanStringGraphDropShortOverlaps(dropRatio, /*minOverlapLen*/0, maxShortTipReads);
            stringGraph.writeGfa(roundPrefix + "6DropShortOverlaps-" + to_string(dropPct) + "pct.gfa", *reads);

            // asg_arc_cut_tips (hifiasm ul_clean_gfa line 3092).
            cutStringGraphTips(maxShortTipReads);
            stringGraph.writeGfa(roundPrefix + "7CutTips.gfa", *reads);

            // asg_arc_cut_bub_links (HARD_OL_DROP=0.6, HARD_OL_SEC_DROP=0.85).
            cleanStringGraphBubbleLinks(/*lenRat*/0.6, /*secLenRat*/0.85, maxShortTipReads);
            stringGraph.writeGfa(roundPrefix + "8BubbleLinks.gfa", *reads);

            // asg_arc_cut_complex_bub_links (HARD_OL_DROP=0.6).
            cleanStringGraphComplexBubbleLinks(/*lenRat*/0.6);
            stringGraph.writeGfa(roundPrefix + "9ComplexBubbleLinks.gfa", *reads);

            // asg_arc_cut_tips (hifiasm ul_clean_gfa line 3107).
            cutStringGraphTips(maxShortTipReads);
            stringGraph.writeGfa(roundPrefix + "10CutTips.gfa", *reads);
        }

        // // Post-loop cleaning (hifiasm ul_clean_gfa lines 3144–3174).

        // // Post-loop: one more semi-circ pass (hifiasm line 3144).
        // cleanStringGraphBreakShortCycles(maxCycleReads, maxShortTipReads);
        // stringGraph.writeGfa("ReadGraph6-StringGraph-PostLoop-1BreakCycles.gfa", *reads);

        // // Post-loop: large-indel arc cut (hifiasm line 3149; NOP for clean HiFi).
        // cleanStringGraphLargeIndelArcs(maxShortTipReads, /*minDiff*/0);
        // cutStringGraphTips(maxShortTipReads);
        // stringGraph.writeGfa("ReadGraph6-StringGraph-PostLoop-2LargeIndel.gfa", *reads);

        // // Post-loop: asg_arc_cut_length at HARD_ORTHOLOGY_DROP=0.4 (hifiasm line 3157).
        // cleanStringGraphDropShortOverlaps(/*HARD_ORTHOLOGY_DROP*/0.4, /*minOverlapLen*/0, maxShortTipReads);
        // cutStringGraphTips(maxShortTipReads);
        // stringGraph.writeGfa("ReadGraph6-StringGraph-PostLoop-3DropHardOrthology.gfa", *reads);

        // // Post-loop: asg_arc_cut_length at min_ovlp_drop_ratio=0.2, l_drop=2000 (hifiasm line 3163).
        // cleanStringGraphDropShortOverlaps(/*minDropRate*/0.2, /*finalMinOverlapLen*/2000, maxShortTipReads);
        // cutStringGraphTips(maxShortTipReads);
        // stringGraph.writeGfa("ReadGraph6-StringGraph-PostLoop-4DropFinal.gfa", *reads);

        // // Post-loop: asg_cut_semi_circ standalone pass (hifiasm line 3174, LIM_LEN=100).
        // cutStringGraphSemiCircular(/*limLen*/100);
        // stringGraph.writeGfa("ReadGraph6-StringGraph-PostLoop-5SemiCircular.gfa", *reads);

        // // post_rescue: SKIPPED — requires coverage/R_to_U infrastructure not in dinara.
    }

    // Step 12a: Sync read graph after iterative cleaning.
    rebuildReadGraphFromCurrentStringGraph(/*rebuildDirectedReadGraph*/false);

    // // Compress the cleaned string graph into unitigs (hifiasm parity: output_unitig_graph / ma_ug_gen).
    // // No separate unitig-graph cleaning: hifiasm cleans the string graph in ul_clean_gfa and then
    // // generates unitigs as a final output step — there is no post-generation unitig cleaning phase.
    // createUnitigGraphFromStringGraph();
    // unitigGraph.writeGfa("ReadGraph6-UnitigGraph.gfa");

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
