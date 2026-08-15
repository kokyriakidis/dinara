// AssemblerReadGraph6.cpp
// Read graph creation using CIGAR-based phasing (Hifiasm parity)

#include "Assembler.hpp"
#include "overlapClassification.hpp"
#include "timestamp.hpp"
#include "Reads.hpp"
#include <vector>
#include <thread>

using namespace dinara;

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

void Assembler::createReadGraphFromPhasingCisOverlaps()
{
    createReadGraphFromPhasingCisOverlaps(
        std::thread::hardware_concurrency(),
        /*rebuildDirectedReadGraph*/ false);
}

void Assembler::createReadGraphFromPhasingCisOverlaps(
    uint64_t /*threadCount*/,
    bool rebuildDirectedReadGraph)
{
    cout << timestamp << "createReadGraphFromPhasingCisOverlaps begins." << endl;
    checkAlignmentDataAreOpen();

    const uint64_t alignmentCount = alignmentData.size();

    vector<uint8_t> keepAlignmentByte(alignmentCount, 0);
    uint64_t keptCount = 0;
    uint64_t filteredDeleted = 0;
    uint64_t filteredTrans = 0;
    uint64_t unlabeledCount = 0;

    for(uint64_t i = 0; i < alignmentCount; i++) {
        const AlignmentData& ad = alignmentData[i];

        // Skip overlaps deleted by earlier stages (chimeric, weak-contradicted,
        // internal, short, contained).  Matches hifiasm's ma_sg_gen which
        // skips overlaps with del=1.
        if(!ad.keptByBothSides()) {
            ++filteredDeleted;
            continue;
        }

        // hifiasmEcMatchState: 0=unlabeled, 1=cis, 2=trans, 3=cisDifferentCopy.
        // Keep if neither side is trans or different-copy.
        // Unlabeled (0) overlaps that survived earlier stages are kept —
        // matches hifiasm where sources[] contains both cis (ml=1) and
        // unlabeled (ml=0) overlaps.
        const bool exclude0 = (ad.hifiasmEcMatchState0 == 2 || ad.hifiasmEcMatchState0 == 3);
        const bool exclude1 = (ad.hifiasmEcMatchState1 == 2 || ad.hifiasmEcMatchState1 == 3);
        if(exclude0 || exclude1) {
            ++filteredTrans;
        } else {
            keepAlignmentByte[i] = 1;
            ++keptCount;
            if(ad.hifiasmEcMatchState0 == 0 || ad.hifiasmEcMatchState1 == 0) {
                ++unlabeledCount;
            }
        }
    }

    cout << timestamp << "Kept " << keptCount << " / " << alignmentCount
         << " alignments (filtered " << filteredDeleted << " deleted, "
         << filteredTrans << " trans)." << endl;
    cout << timestamp << "Of kept: " << unlabeledCount
         << " have at least one unlabeled side." << endl;

    vector<bool> keepAlignment(alignmentCount, false);
    for(uint64_t i = 0; i < alignmentCount; ++i) {
        keepAlignment[i] = (keepAlignmentByte[i] != 0);
    }

    rebuildReadGraphUsingSelectedAlignments(std::move(keepAlignment), rebuildDirectedReadGraph);

    cout << timestamp << "createReadGraphFromPhasingCisOverlaps completed." << endl;
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


/*
asg_arc_del_trans (hifiasm) — exact port for the ReadGraph.

Transitive reduction removes redundant edges from the overlap graph.
An edge v->x is transitive (redundant) if there exists an intermediate
vertex w such that v->w->x, with len(v->w) + len(w->x) <= len(v->x) + fuzz.

Builds proper string graph arcs using ma_hit2arc_full to compute the
correct arc direction and node length (unique extension of the source
read). Each alignment produces up to 2 directed arcs (one from each
read's perspective), matching hifiasm's ma_sg_gen which processes each
overlap from both reads' sources arrays.

Arcs are sorted by arc length (ascending) per vertex. The fuzz parameter
(default 5000bp) allows for small length discrepancies from inexact
overlaps.

After marking transitive arcs, the corresponding alignments are removed
and the read graph is rebuilt.
*/
uint64_t Assembler::transitiveReductionOnReadGraph(int32_t fuzz)
{
    cout << timestamp << "transitiveReductionOnReadGraph begins (fuzz=" << fuzz << ")." << endl;

    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();

    const uint32_t vertexCount = uint32_t(reads->readCount() * 2);
    const uint64_t edgeCount = readGraph.edges.size();

    if (edgeCount == 0) {
        cout << timestamp << "transitiveReductionOnReadGraph: no edges." << endl;
        return 0;
    }

    // Build per-vertex directed arc lists using ma_hit2arc_full to compute
    // correct string graph arc directions and node lengths.
    struct Arc {
        uint32_t targetVertex;  // readId<<1 | strand
        uint32_t arcLen;        // node length (unique extension of source read)
        uint32_t alignmentId;
    };

    vector<vector<Arc>> vertexArcs(vertexCount);
    uint64_t skippedNonDovetail = 0;

    // Process each alignment once. ma_hit2arc_full is called twice per
    // alignment (once with each read as query), producing up to 2 directed
    // arcs. This matches hifiasm's ma_sg_gen which iterates sources[i] for
    // all reads, where each overlap appears in both reads' source arrays.

    // Collect the set of alignment IDs present in the read graph.
    const uint64_t alignmentCount = alignmentData.size();
    vector<bool> inReadGraph(alignmentCount, false);
    for (uint64_t edgeId = 0; edgeId < edgeCount; edgeId += 2) {
        const uint32_t alignmentId = uint32_t(readGraph.edges[edgeId].alignmentId);
        if (alignmentId < alignmentCount) {
            inReadGraph[alignmentId] = true;
        }
    }

    for (uint32_t alignmentId = 0; alignmentId < alignmentCount; alignmentId++) {
        if (!inReadGraph[alignmentId]) continue;
        const AlignmentData& ad = alignmentData[alignmentId];

        const ReadId queryId = ad.readIds[0];
        const ReadId targetId = ad.readIds[1];
        const int32_t ql = int32_t(reads->getReadRawSequenceLength(queryId));
        const int32_t tl = int32_t(reads->getReadRawSequenceLength(targetId));
        const bool rev = !ad.isSameStrand;

        // Arc from query's perspective (matches hifiasm sources[queryId]).
        {
            DovetailArc sgArc;
            const int r = ma_hit2arc_full(
                queryId, targetId,
                int32_t(ad.qs), int32_t(ad.qe), ql,
                int32_t(ad.ts), int32_t(ad.te), tl,
                rev,
                /*max_hang*/ 1000, /*int_frac*/ 0.8f, /*min_ovlp*/ 50,
                sgArc);
            if (r >= 0) {
                vertexArcs[sgArc.sourceVertex].push_back({
                    sgArc.targetVertex, sgArc.arcLen, alignmentId});
            } else {
                ++skippedNonDovetail;
            }
        }

        // Arc from target's perspective (matches hifiasm sources[targetId]).
        // Swap query/target coordinates; rev stays the same.
        {
            DovetailArc sgArc;
            const int r = ma_hit2arc_full(
                targetId, queryId,
                int32_t(ad.ts), int32_t(ad.te), tl,
                int32_t(ad.qs), int32_t(ad.qe), ql,
                rev,
                /*max_hang*/ 1000, /*int_frac*/ 0.8f, /*min_ovlp*/ 50,
                sgArc);
            if (r >= 0) {
                vertexArcs[sgArc.sourceVertex].push_back({
                    sgArc.targetVertex, sgArc.arcLen, alignmentId});
            } else {
                ++skippedNonDovetail;
            }
        }
    }

    // Sort each vertex's arcs by arc length (ascending = longest overlap first).
    for (uint32_t v = 0; v < vertexCount; v++) {
        sort(vertexArcs[v].begin(), vertexArcs[v].end(),
            [](const Arc& a, const Arc& b) { return a.arcLen < b.arcLen; });
    }

    // Transitive reduction (asg_arc_del_trans).
    // mark[w]: 0=vacant, 1=neighbor of v, 2=reachable via two-hop (transitive).
    vector<uint8_t> mark(vertexCount, 0);
    vector<bool> alignmentRemoved(alignmentCount, false);
    uint64_t reducedCount = 0;

    for (uint32_t v = 0; v < vertexCount; v++) {
        const auto& arcsV = vertexArcs[v];
        const uint32_t nv = uint32_t(arcsV.size());
        if (nv == 0) continue;

        for (uint32_t i = 0; i < nv; i++) {
            mark[arcsV[i].targetVertex] = 1;
        }

        const uint32_t maxLen = arcsV[nv - 1].arcLen + uint32_t(fuzz);

        for (uint32_t i = 0; i < nv; i++) {
            const uint32_t w = arcsV[i].targetVertex;
            if (mark[w] != 1) continue;

            const auto& arcsW = vertexArcs[w];
            for (uint32_t j = 0; j < uint32_t(arcsW.size()); j++) {
                if (arcsW[j].arcLen + arcsV[i].arcLen > maxLen) break;
                if (mark[arcsW[j].targetVertex]) {
                    mark[arcsW[j].targetVertex] = 2;
                }
            }
        }

        for (uint32_t i = 0; i < nv; i++) {
            if (mark[arcsV[i].targetVertex] == 2) {
                alignmentRemoved[arcsV[i].alignmentId] = true;
                ++reducedCount;
            }
            mark[arcsV[i].targetVertex] = 0;
        }
    }

    // Rebuild the read graph without removed alignments.
    vector<bool> keepAlignment(alignmentCount, false);
    uint64_t removedAlignments = 0;
    for (uint32_t alignmentId = 0; alignmentId < alignmentCount; alignmentId++) {
        if (inReadGraph[alignmentId] && !alignmentRemoved[alignmentId]) {
            keepAlignment[alignmentId] = true;
        }
    }
    for (uint32_t alignmentId = 0; alignmentId < alignmentCount; alignmentId++) {
        if (inReadGraph[alignmentId] && alignmentRemoved[alignmentId]) {
            ++removedAlignments;
        }
    }

    rebuildReadGraphUsingSelectedAlignments(keepAlignment, false);

    cout << timestamp << "transitiveReductionOnReadGraph: removed "
         << removedAlignments << " alignments (" << reducedCount
         << " directed arcs reduced, " << skippedNonDovetail
         << " non-dovetail arcs skipped)." << endl;

    return removedAlignments;
}


