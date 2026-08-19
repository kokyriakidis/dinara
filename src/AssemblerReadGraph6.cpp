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
