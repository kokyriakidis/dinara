// Dinara.
#include "Assembler.hpp"
#include "Reads.hpp"
#include "overlapClassification.hpp"
#include "timestamp.hpp"
using namespace dinara;

// Standard library.
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>



void Assembler::filterOverlapsByRegionalCliques(
    uint64_t minIntervalOverlap,
    uint64_t minRegionSize,
    double minCliqueFraction,
    uint64_t /* threadCount */)
{
    const auto tBegin = std::chrono::steady_clock::now();
    cout << timestamp << "filterOverlapsByRegionalCliques begins "
         << "(minIntervalOverlap=" << minIntervalOverlap
         << ", minRegionSize=" << minRegionSize
         << ", minCliqueFraction=" << minCliqueFraction << ")." << endl;

    checkAlignmentDataAreOpen();
    DINARA_ASSERT(alignmentTable.isOpen());

    const ReadId readCount = reads->readCount();
    const uint64_t alignmentCount = alignmentData.size();

    // Sort reads by length descending for cascading processing.
    vector<ReadId> readOrder(readCount);
    std::iota(readOrder.begin(), readOrder.end(), ReadId(0));
    std::sort(readOrder.begin(), readOrder.end(),
        [this](ReadId a, ReadId b) {
            return reads->getReadRawSequenceLength(a) > reads->getReadRawSequenceLength(b);
        });

    // Helper: check if two intervals [s1,e1) and [s2,e2) overlap by at least minOverlap bases.
    auto intervalsOverlap = [](uint32_t s1, uint32_t e1, uint32_t s2, uint32_t e2,
                               uint64_t minOverlap) -> bool {
        if (s1 >= e2 || s2 >= e1) return false;
        const uint32_t overlapStart = std::max(s1, s2);
        const uint32_t overlapEnd = std::min(e1, e2);
        return (overlapEnd > overlapStart) && (uint64_t(overlapEnd - overlapStart) >= minOverlap);
    };

    // Per-candidate info for processing one read.
    struct CandidateInfo {
        uint32_t alignmentId;     // index into alignmentData
        ReadId partnerId;         // physical read id of the partner
        uint32_t startOnR;        // interval start on R
        uint32_t endOnR;          // interval end on R
        int64_t dpScore;          // affine gap DP score
        bool isDovetail;          // dovetail vs contained/internal

        // Sort key: dovetail first, then by dpScore descending.
        bool operator<(const CandidateInfo& o) const {
            if (isDovetail != o.isDovetail) return isDovetail > o.isDovetail; // dovetail first
            return dpScore > o.dpScore; // higher score first
        }
    };

    // Accepted read entry in the interval set.
    struct AcceptedEntry {
        uint32_t startOnR;
        uint32_t endOnR;
        ReadId readId;
    };

    // Helper: check if readA and readB have a non-deleted overlap via alignmentTable.
    // Linear scan of readA's alignmentTable section (typically 50-200 entries).
    auto haveNonDeletedOverlap = [this](ReadId readA, ReadId readB) -> bool {
        const OrientedReadId orientedReadIdA(readA, 0);
        const auto section = alignmentTable[orientedReadIdA.getValue()];
        for (const uint32_t alignmentIdx : section) {
            const AlignmentData& ad = alignmentData[alignmentIdx];
            const ReadId otherRead = (ad.readIds[0] == readA) ? ad.readIds[1] : ad.readIds[0];
            if (otherRead == readB) {
                return ad.keptByBothSides();
            }
        }
        return false;
    };

    // Reusable buffers.
    vector<CandidateInfo> candidates;
    vector<AcceptedEntry> acceptedSet;
    vector<uint32_t> overlappingIndices; // indices into acceptedSet

    uint64_t totalFlagged = 0;
    uint64_t totalAccepted = 0;
    uint64_t readsProcessed = 0;

    // Parameters for ma_hit2arc_containment (same defaults as filterHangingOverlaps).
    const int32_t maxHang = 1000;
    const double intFrac = 0.8;
    const int32_t minOvlp = 50;

    for (const ReadId readId : readOrder) {
        const uint64_t rLen = reads->getReadRawSequenceLength(readId);
        if (rLen == 0) continue;

        // Gather all non-deleted overlaps for this read.
        const OrientedReadId orientedReadId(readId, 0);
        const auto section = alignmentTable[orientedReadId.getValue()];

        candidates.clear();
        for (const uint32_t alignmentIdx : section) {
            AlignmentData& ad = alignmentData[alignmentIdx];
            if (!ad.keptByBothSides()) continue;

            // Determine which read is R and which is the partner.
            const bool isRead0 = (ad.readIds[0] == readId);
            const ReadId partnerId = isRead0 ? ad.readIds[1] : ad.readIds[0];

            // Compute interval on R.
            // AlignmentData stores: readIds[0] on strand 0, qs/qe for readIds[0], ts/te for readIds[1].
            uint32_t startOnR, endOnR;
            if (isRead0) {
                startOnR = ad.qs;
                endOnR = ad.qe;
            } else {
                startOnR = ad.ts;
                endOnR = ad.te;
            }

            // Skip degenerate intervals.
            if (endOnR <= startOnR) continue;

            // Classify using the same hifiasm ma_hit2arc_containment function.
            // The function expects query=R, target=partner.
            // When R is readIds[0]: qs/qe are R's coords, ts/te are partner's coords.
            // When R is readIds[1]: swap so R is query.
            int32_t classQs, classQe, classQl, classTs, classTe, classTl;
            bool classIsReverse;
            if (isRead0) {
                classQs = (int32_t)ad.qs;
                classQe = (int32_t)ad.qe;
                classQl = (int32_t)rLen;
                classTs = (int32_t)ad.ts;
                classTe = (int32_t)ad.te;
                classTl = (int32_t)reads->getReadRawSequenceLength(partnerId);
                classIsReverse = !ad.isSameStrand;
            } else {
                // R is readIds[1], partner is readIds[0].
                // Swap query/target so R is the query.
                classQs = (int32_t)ad.ts;
                classQe = (int32_t)ad.te;
                classQl = (int32_t)rLen;
                classTs = (int32_t)ad.qs;
                classTe = (int32_t)ad.qe;
                classTl = (int32_t)reads->getReadRawSequenceLength(partnerId);
                classIsReverse = !ad.isSameStrand;
            }

            const int classification = ma_hit2arc_containment(
                classQs, classQe, classQl,
                classTs, classTe, classTl,
                classIsReverse, maxHang, intFrac, minOvlp);

            // result 0 = dovetail, 1 = R contained, 2 = partner contained
            // result < 0 = internal/too-short (should already be filtered, but handle gracefully)
            const bool isDovetail = (classification == 0);

            candidates.push_back({
                alignmentIdx,
                partnerId,
                startOnR,
                endOnR,
                ad.info.dpScore,
                isDovetail
            });
        }

        // Sort: dovetail first, then by dpScore descending.
        std::sort(candidates.begin(), candidates.end());

        // Greedy interval-based clique building.
        acceptedSet.clear();

        for (auto& cand : candidates) {
            // Find accepted reads whose intervals overlap this candidate's interval
            // by at least minIntervalOverlap.
            overlappingIndices.clear();
            for (uint32_t i = 0; i < acceptedSet.size(); i++) {
                if (intervalsOverlap(acceptedSet[i].startOnR, acceptedSet[i].endOnR,
                                     cand.startOnR, cand.endOnR, minIntervalOverlap)) {
                    overlappingIndices.push_back(i);
                }
            }

            if (overlappingIndices.empty()) {
                // New region -- accept unconditionally.
                acceptedSet.push_back({cand.startOnR, cand.endOnR, cand.partnerId});
                totalAccepted++;
                continue;
            }

            if (overlappingIndices.size() < minRegionSize) {
                // Not enough evidence to reject -- accept.
                acceptedSet.push_back({cand.startOnR, cand.endOnR, cand.partnerId});
                totalAccepted++;
                continue;
            }

            // Check clique property: candidate must overlap at least minCliqueFraction
            // of the accepted reads covering this region. This allows diploid/polyploid
            // regions where a read from one haplotype won't overlap reads from the other.
            const uint64_t required = (uint64_t)std::ceil(overlappingIndices.size() * minCliqueFraction);
            uint64_t cliqueHits = 0;
            for (const uint32_t idx : overlappingIndices) {
                if (haveNonDeletedOverlap(cand.partnerId, acceptedSet[idx].readId)) {
                    cliqueHits++;
                }
            }
            const bool passesClique = (cliqueHits >= required);

            if (passesClique) {
                acceptedSet.push_back({cand.startOnR, cand.endOnR, cand.partnerId});
                totalAccepted++;
            } else {
                // Flag this overlap as spurious from both perspectives.
                alignmentData[cand.alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonClique);
                totalFlagged++;
            }
        }

        readsProcessed++;
        if (readsProcessed % 100000 == 0) {
            cout << timestamp << "  Processed " << readsProcessed << " / " << readCount
                 << " reads. Accepted: " << totalAccepted << ", Flagged: " << totalFlagged << endl;
        }
    }

    const auto tEnd = std::chrono::steady_clock::now();
    const double tElapsed = std::chrono::duration<double>(tEnd - tBegin).count();

    cout << timestamp << "filterOverlapsByRegionalCliques completed in "
         << std::fixed << std::setprecision(1) << tElapsed << " s." << endl;
    cout << timestamp << "  Reads processed: " << readsProcessed << endl;
    cout << timestamp << "  Overlaps accepted: " << totalAccepted << endl;
    cout << timestamp << "  Overlaps flagged: " << totalFlagged
         << " / " << alignmentCount << " total alignments" << endl;
}
