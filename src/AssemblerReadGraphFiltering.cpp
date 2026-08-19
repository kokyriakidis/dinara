/*
Overlap/read-graph filtering (hifiasm/miniasm parity)

This file implements the overlap-cleaning pipeline used before building the string graph.
It mirrors hifiasm/miniasm stages where possible, while preserving Dinara-specific state
like EC/phasing decisions.

Key data and conventions
- `alignmentData[alignmentId]` stores an undirected overlap between `readIds[0]` and `readIds[1]`.
- Overlap coordinates are stored in each read's forward coordinate system:
  - `qs/qe` are on `readIds[0]`
  - `ts/te` are on `readIds[1]`
  - `isSameStrand == false` indicates reverse-complement orientation, but coordinates remain forward.
- Deletion/filters are recorded as per-side bitmasks (`deleteReasons0/deleteReasons1`).
  An overlap is active for graph construction only when both sides have no deletion reasons
  (`keptByBothSides()`).

Hifiasm stage mapping (conceptual)
- gen_chemical_arc_rf        -> applyOntChemicalArcMask (ONT-only)
- try_rescue_overlaps        -> rescuePhasedOverlaps
- ma_hit_cut                 -> applyCoverageCuts
- ma_hit_flt                 -> filterHangingOverlaps
*/

#include "Assembler.hpp"
#include "HifiasmBoundaryVerify.hpp"
#include "Reads.hpp"
#include "hifiasmCoordinateTransforms.hpp"
#include "overlapClassification.hpp"
#include "timestamp.hpp"
#include <algorithm>
#include <vector>
#include <limits>
#include <random>

#include <iostream>
#include <iomanip>

using namespace dinara;
using namespace std;

namespace {
    inline bool boundaryVerify(
        const Reads& reads,
        uint32_t qIntervalStart,
        uint32_t qIntervalEnd,
        ReadId qId,
        ReadId tId,
        uint32_t qs,
        uint32_t ts,
        uint32_t te,
        bool rev,
        std::vector<char>& xBuf,
        std::vector<char>& yBuf)
    {
        return hifiasmBoundaryVerify(
            reads,
            qIntervalStart,
            qIntervalEnd,
            qId,
            tId,
            qs,
            ts,
            te,
            rev,
            xBuf,
            yBuf);
    }
}

/*
gen_chemical_arc_rf (hifiasm ONT) parity: chemical-arc masking.

This stage is ONT-specific and intentionally ignores EC/phasing deletions:
hifiasm resets `del` before this stage and evaluates depth using all overlaps.

We compute the minimum overlap depth across the read after trimming `chemicalFlank`
from both ends of each overlap, and flag reads with `minDepth <= chemicalCov`.
All overlaps incident to flagged reads get `DeleteReasonChemical` on both sides.
*/
void Assembler::applyOntChemicalArcMask(uint64_t threadCount)
{
    applyOntChemicalArcMask(chemicalArcCov, chemicalArcFlank, chemicalArcDupRate, threadCount);
}

void Assembler::applyOntChemicalArcMask(uint64_t chemicalCov, uint64_t chemicalFlank, double dupRate, uint64_t threadCount)
{
    cout << timestamp << "Applying ONT chemical arc mask..." << endl;

    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();

    chemicalArcCov = chemicalCov;
    chemicalArcFlank = chemicalFlank;
    chemicalArcDupRate = dupRate;

    for (uint64_t i = 0; i < alignmentData.size(); ++i) {
        alignmentData[i].clearDeleteReasonsBoth(AlignmentData::DeleteReasonChemical);
    }

    chemicalArcMask.assign(readCount, uint8_t(0xFF));

    setupLoadBalancing(readCount, 64);
    runThreads(&Assembler::applyOntChemicalArcMaskThreadFunction, threadCount);

    uint64_t chemReadCount = 0;
    for (uint64_t r = 0; r < chemicalArcMask.size(); ++r) {
        if (chemicalArcMask[r] <= chemicalArcCov) {
            ++chemReadCount;
        }
    }
    cout << timestamp << "ONT chemical arc mask: " << chemReadCount << " reads flagged." << endl;

    if (chemReadCount) {
        for (uint64_t alignmentId = 0; alignmentId < alignmentData.size(); ++alignmentId) {
            AlignmentData& ad = alignmentData[alignmentId];
            const ReadId r0 = ad.readIds[0];
            const ReadId r1 = ad.readIds[1];
            const bool del =
                (r0 < chemicalArcMask.size() && chemicalArcMask[r0] <= chemicalArcCov) ||
                (r1 < chemicalArcMask.size() && chemicalArcMask[r1] <= chemicalArcCov);
            if (del) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonChemical);
            }
        }
    }

    cout << timestamp << "ONT chemical arc mask complete." << endl;
}

void Assembler::applyOntChemicalArcMaskThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    const uint64_t flank = chemicalArcFlank;
    const uint64_t covCut = chemicalArcCov;
    const double dupRate = chemicalArcDupRate;

    std::vector<uint64_t> events;
    events.reserve(256);

    uint64_t begin, end;
    while (getNextBatch(begin, end)) {
        for (ReadId r = ReadId(begin); r < ReadId(end); ++r) {
            const uint32_t len = uint32_t(reads->getReadRawSequenceLength(r));
            events.clear();

            const OrientedReadId oid(r, 0);
            if (oid.getValue() < alignmentTable.size()) {
                const auto& table = alignmentTable[oid.getValue()];
                for (uint32_t alignmentId : table) {
                    const AlignmentData& ad = alignmentData[alignmentId];

                    uint32_t qs0 = 0, qe0 = 0, os0 = 0, oe0 = 0;
                    ReadId other = invalidReadId;
                    if (ad.readIds[0] == r) {
                        qs0 = ad.qs; qe0 = ad.qe;
                        os0 = ad.ts; oe0 = ad.te;
                        other = ad.readIds[1];
                    } else if (ad.readIds[1] == r) {
                        qs0 = ad.ts; qe0 = ad.te;
                        os0 = ad.qs; oe0 = ad.qe;
                        other = ad.readIds[0];
                    } else {
                        continue;
                    }

                    if (qe0 <= qs0) continue;

                    uint32_t qs = qs0;
                    uint32_t qe = qe0;

                    if (qs > 0) {
                        const uint64_t qs64 = uint64_t(qs) + flank;
                        qs = (qs64 > uint64_t(len)) ? len : uint32_t(qs64);
                    }
                    if (qe < len) {
                        qe = (qe > flank) ? uint32_t(uint64_t(qe) - flank) : 0U;
                    }
                    if (qe <= qs) continue;

                    if (!ad.isSameStrand) {
                        const uint32_t otherLen = uint32_t(reads->getReadRawSequenceLength(other));
                        const uint32_t rr = (otherLen >= len) ? (otherLen - len) : (len - otherLen);
                        if (double(rr) <= double(len) * dupRate && double(rr) <= double(otherLen) * dupRate) {
                            const uint32_t uncoveredThis = len - (qe0 - qs0);
                            const uint32_t uncoveredOther = otherLen - (oe0 - os0);
                            if (double(uncoveredThis) <= double(len) * dupRate &&
                                double(uncoveredOther) <= double(otherLen) * dupRate) {
                                continue;
                            }
                        }
                    }

                    events.push_back(uint64_t(qs) << 1);
                    events.push_back((uint64_t(qe) << 1) | 1ULL);
                }
            }

            std::sort(events.begin(), events.end());

            int32_t dp = 0;
            uint32_t st = 0;
            int32_t minCov = std::numeric_limits<int32_t>::max();

            for (const uint64_t ev : events) {
                const int32_t oldDp = dp;
                if (ev & 1ULL) {
                    --dp;
                } else {
                    ++dp;
                }
                const uint32_t pos = uint32_t(ev >> 1);
                if (pos > st) {
                    if (oldDp < minCov) minCov = oldDp;
                    st = pos;
                }
            }
            if (len > st) {
                if (dp < minCov) minCov = dp;
            }
            if (minCov == std::numeric_limits<int32_t>::max()) {
                minCov = 0;
            }

            if (uint64_t(minCov) <= covCut) {
                chemicalArcMask[r] = uint8_t(minCov);
            } else {
                chemicalArcMask[r] = uint8_t(0xFF);
            }
        }
    }
}



/*
ma_hit_cut (miniasm/hifiasm) parity: clip active overlaps to the per-read valid intervals.

This stage:
- Skips overlaps already deleted by earlier filters.
- Projects overlap coordinates into each read's valid interval and normalizes them to be
  relative to the valid-region start.
- Marks overlaps as `DeleteReasonCoverageCut` when clipping makes them too short.

Important: some upstream stages can store tip-extended overlap bounds; for parity with hifiasm,
we reconstruct raw (non-extended) overlap bounds from aligned marker ordinals when available.
*/
void Assembler::applyCoverageCuts(uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Applying coverage cuts (ma_hit_cut equivalent, minOverlapLength="
         << minOverlapLength << ")." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    this->coverageCutMinOverlap = minOverlapLength;

    // Count overlaps before clipping
    uint64_t overlapsBeforeCut = 0;
    for (const auto& ad : alignmentData) {
        if (ad.keptByBothSides()) {
            ++overlapsBeforeCut;
        }
    }

    // Hifiasm ma_hit_cut: clip overlaps to valid regions and normalize coordinates
    setupLoadBalancing(alignmentData.size(), 10000);
    runThreads(&Assembler::applyCoverageCutsToAlignmentsThreadFunction, threadCount);

    // Count overlaps after clipping
    uint64_t overlapsAfterCut = 0;
    uint64_t overlapsDeletedByCut = 0;
    for (const auto& ad : alignmentData) {
        if (ad.keptByBothSides()) {
            ++overlapsAfterCut;
        } else if ((ad.deleteReasons0 & AlignmentData::DeleteReasonCoverageCut) ||
                   (ad.deleteReasons1 & AlignmentData::DeleteReasonCoverageCut)) {
            ++overlapsDeletedByCut;
        }
    }

    cout << timestamp << "Coverage cuts: " << overlapsDeletedByCut << " overlaps deleted ("
         << overlapsBeforeCut << " -> " << overlapsAfterCut << " remaining)." << endl;

    // Hifiasm ma_hit_cut: delete reads that have no remaining overlaps (coverage_cut[].del = 1)
    setupLoadBalancing(reads->readCount(), 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);

    // Count reads deleted after cleanup
    uint64_t readsDeleted = 0;
    for (const auto& status : validReadIntervals) {
        if (status.isDeleted) {
            ++readsDeleted;
        }
    }

    cout << timestamp << "Coverage cuts complete: " << readsDeleted << " reads marked deleted (no surviving overlaps)." << endl;
}

void Assembler::applyCoverageCutsToAlignmentsThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    uint64_t begin, end;
    const uint64_t minLen = this->coverageCutMinOverlap;

    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; i++) {
            AlignmentData& ad = alignmentData[i];

            // Hifiasm ma_hit_cut: only process overlaps that are kept by both sides
            if(!ad.keptByBothSides()) continue;

            ReadId qn = ad.readIds[0];
            ReadId tn = ad.readIds[1];

            if (qn >= validReadIntervals.size() || tn >= validReadIntervals.size()) continue;

            const auto& rq = validReadIntervals[qn];
            const auto& rt = validReadIntervals[tn];

            // Hifiasm ma_hit_cut: skip overlaps if either endpoint read is deleted (coverage_cut[].del == 1)
            if (rq.isDeleted || rt.isDeleted) continue;

            // Get read lengths for validation
            const uint32_t qLen = uint32_t(reads->getReadRawSequenceLength(qn));
            const uint32_t tLen = uint32_t(reads->getReadRawSequenceLength(tn));

            int32_t qs, qe, ts, te;

            /*
            Hifiasm ma_hit_cut operates on raw overlap bounds (not tip-extended).
            When available, reconstruct raw bounds from aligned marker ordinals.
            This avoids over-calling containment due to tip extensions.
            Fallback: use stored ad.qs/qe/ts/te (which may include tip extensions).
            */
            auto getRawBoundsIfAvailable = [&](uint32_t& outQs, uint32_t& outQe, uint32_t& outTs, uint32_t& outTe) -> bool {
                if (!markers || !markers->isOpen()) return false;
                const auto& d0 = ad.info.data[0];
                const auto& d1 = ad.info.data[1];
                if (d0.markerCount == 0 || d1.markerCount == 0) return false;
                if (d0.firstOrdinal > d0.lastOrdinal || d1.firstOrdinal > d1.lastOrdinal) return false;

                const OrientedReadId oid0(ad.readIds[0], 0);
                const OrientedReadId oid1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
                const auto& m0 = (*markers)[oid0.getValue()];
                const auto& m1 = (*markers)[oid1.getValue()];
                if (d0.lastOrdinal >= m0.size() || d1.lastOrdinal >= m1.size()) return false;

                const uint32_t k = uint32_t(assemblerInfo->k);
                const uint32_t qs0 = m0[d0.firstOrdinal].position;
                const uint32_t qe0 = m0[d0.lastOrdinal].position + k;

                const uint32_t tsOriented = m1[d1.firstOrdinal].position;
                const uint32_t teOriented = m1[d1.lastOrdinal].position + k;

                uint32_t ts0 = tsOriented;
                uint32_t te0 = teOriented;
                if (!ad.isSameStrand) {
                    const auto p = dinara::rcIntervalToForward(tLen, tsOriented, teOriented);
                    ts0 = p.first;
                    te0 = p.second;
                }

                // Validate bounds before returning
                if (qs0 >= qe0 || qe0 > qLen) return false;
                if (ts0 >= te0 || te0 > tLen) return false;

                outQs = qs0; outQe = qe0; outTs = ts0; outTe = te0;
                return true;
            };

            // Start with stored bounds, then try to reconstruct from markers
            uint32_t rawQs = ad.qs, rawQe = ad.qe, rawTs = ad.ts, rawTe = ad.te;
            (void)getRawBoundsIfAvailable(rawQs, rawQe, rawTs, rawTe);

            // Validate raw bounds are within read lengths (defensive check)
            if (rawQs >= rawQe || rawQe > qLen || rawTs >= rawTe || rawTe > tLen) {
                // Invalid bounds - delete overlap
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonCoverageCut);
                continue;
            }

            /*
            Hifiasm ma_hit_cut coordinate clipping logic:

            For same-strand overlaps:
            - If target extends beyond query's valid region, adjust query bounds accordingly
            - If query extends beyond target's valid region, adjust target bounds accordingly

            For reverse-strand overlaps:
            - Same logic but accounting for RC orientation

            Note: Using strict inequalities (< and >) to match original hifiasm implementation exactly.
            */
            if (!ad.isSameStrand) {
                // Reverse-complement case
                qs = (int32_t)rawTe < (int32_t)rt.end ? (int32_t)rawQs : (int32_t)rawQs + ((int32_t)rawTe - (int32_t)rt.end);
                qe = (int32_t)rawTs > (int32_t)rt.start ? (int32_t)rawQe : (int32_t)rawQe - ((int32_t)rt.start - (int32_t)rawTs);
                ts = (int32_t)rawQe < (int32_t)rq.end ? (int32_t)rawTs : (int32_t)rawTs + ((int32_t)rawQe - (int32_t)rq.end);
                te = (int32_t)rawQs > (int32_t)rq.start ? (int32_t)rawTe : (int32_t)rawTe - ((int32_t)rq.start - (int32_t)rawQs);
            } else {
                // Same-strand case
                qs = (int32_t)rawTs > (int32_t)rt.start ? (int32_t)rawQs : (int32_t)rawQs + ((int32_t)rt.start - (int32_t)rawTs);
                qe = (int32_t)rawTe < (int32_t)rt.end ? (int32_t)rawQe : (int32_t)rawQe - ((int32_t)rawTe - (int32_t)rt.end);
                ts = (int32_t)rawQs > (int32_t)rq.start ? (int32_t)rawTs : (int32_t)rawTs + ((int32_t)rq.start - (int32_t)rawQs);
                te = (int32_t)rawQe < (int32_t)rq.end ? (int32_t)rawTe : (int32_t)rawTe - ((int32_t)rawQe - (int32_t)rq.end);
            }

            // Hifiasm ma_hit_cut: clamp to valid intervals (boundary safety)
            qs = std::max(qs, (int32_t)rq.start);
            qe = std::min(qe, (int32_t)rq.end);
            ts = std::max(ts, (int32_t)rt.start);
            te = std::min(te, (int32_t)rt.end);

            // Hifiasm ma_hit_cut: normalize to 0-based coordinates relative to valid region start
            int32_t norm_qs = qs - (int32_t)rq.start;
            int32_t norm_qe = qe - (int32_t)rq.start;
            int32_t norm_ts = ts - (int32_t)rt.start;
            int32_t norm_te = te - (int32_t)rt.start;

            // Hifiasm ma_hit_cut: delete overlap if either side is too short or invalid after clipping
            // Check: both sides >= minLen AND both sides have positive length
            if ((norm_qe - norm_qs >= (int32_t)minLen) && (norm_te - norm_ts >= (int32_t)minLen) &&
                (norm_qe > norm_qs) && (norm_te > norm_ts) &&
                norm_qs >= 0 && norm_ts >= 0) {
                // Store normalized coordinates
                ad.qs = (uint32_t)norm_qs;
                ad.qe = (uint32_t)norm_qe;
                ad.ts = (uint32_t)norm_ts;
                ad.te = (uint32_t)norm_te;
            } else {
                // Hifiasm ma_hit_cut: delete overlap if it doesn't meet requirements
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonCoverageCut);
            }
        }
    }
}

/*
Hifiasm ma_hit_cut/ma_hit_flt/ma_hit_contained_advance cleanup:
After structural overlap filtering, hifiasm deletes reads that have no remaining active overlaps.
This matches the `rLen==0 -> coverage_cut[i].del=1` behavior.

A read is deleted if it has no surviving overlap edges to non-deleted reads.
*/
void Assembler::applyCoverageCutsCleanupThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId r=ReadId(begin); r!=ReadId(end); r++) {
            if (r >= validReadIntervals.size()) continue;

            // Hifiasm: skip already-deleted reads
            if (validReadIntervals[r].isDeleted) continue;

            // Check if read has any surviving overlap to a non-deleted read
            bool hasSurvivingEdge = false;

            const auto& table = alignmentTable[OrientedReadId(r, 0).getValue()];
            for (uint32_t alignmentId : table) {
                const auto& ad = alignmentData[alignmentId];

                // Hifiasm: only count overlaps that are kept by both sides
                if(!ad.keptByBothSides()) continue;

                // Get the other read in this overlap
                const ReadId other = (ad.readIds[0] == r) ? ad.readIds[1] : ad.readIds[0];

                // Hifiasm: don't count overlaps to deleted reads
                if (other < validReadIntervals.size() && validReadIntervals[other].isDeleted) {
                    continue;
                }

                // Found at least one valid overlap
                hasSurvivingEdge = true;
                break;
            }

            // Hifiasm: delete read if it has no surviving overlaps (rLen == 0)
            if (!hasSurvivingEdge) {
                validReadIntervals[r].isDeleted = true;
            }
        }
    }
}


/*
Hifiasm ma_hit_flt parity: remove internal/non-dovetail overlaps with excessive overhangs.

This stage runs after ma_hit_cut, so overlap coordinates are already normalized to valid regions.

Hifiasm's ma_hit2arc logic classifies each overlap as:
  - Dovetail (keep): proper prefix-suffix overlap, small overhangs
  - Containment (keep): one read fully within another (contained reads are
    only flagged later by flagContainedReads, not removed)
  - Internal (delete): both reads have large overhangs, indicating spurious internal match
  - Too short (delete): effective overlap length < min_ovlp after accounting for extensions

Parameters match hifiasm defaults:
  - max_hang = 1000 (asm_opt.max_hang)
  - int_frac = 0.8 (asm_opt.int_frac, called maxHangRate here)
  - min_ovlp = 50 (asm_opt.min_ovlp)
*/
void Assembler::filterHangingOverlaps(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Filtering hanging overlaps (ma_hit_flt equivalent, maxHang="
         << maxHang << ", maxHangRate=" << maxHangRate << ", minOverlapLength="
         << minOverlapLength << ")." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    this->hangingFilterMaxHang = maxHang;
    this->hangingFilterMaxHangRate = maxHangRate;
    this->hangingFilterMinOverlap = minOverlapLength;

    // Count overlaps before filtering
    uint64_t overlapsBeforeHang = 0;
    for (const auto& ad : alignmentData) {
        if (ad.keptByBothSides()) {
            ++overlapsBeforeHang;
        }
    }

    // Hifiasm ma_hit_flt: filter internal/bad overlaps
    setupLoadBalancing(alignmentData.size(), 10000);
    runThreads(&Assembler::filterHangingOverlapsThreadFunction, threadCount);

    // Count overlaps after filtering
    uint64_t overlapsAfterHang = 0;
    uint64_t overlapsDeletedByHang = 0;
    for (const auto& ad : alignmentData) {
        if (ad.keptByBothSides()) {
            ++overlapsAfterHang;
        } else if ((ad.deleteReasons0 & AlignmentData::DeleteReasonHanging) ||
                   (ad.deleteReasons1 & AlignmentData::DeleteReasonHanging)) {
            ++overlapsDeletedByHang;
        }
    }

    cout << timestamp << "Hanging overlaps: " << overlapsDeletedByHang << " internal/bad overlaps deleted ("
         << overlapsBeforeHang << " -> " << overlapsAfterHang << " remaining)." << endl;

    // Hifiasm ma_hit_flt: delete reads that have no remaining overlaps
    setupLoadBalancing(reads->readCount(), 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);

    // Count reads deleted after cleanup
    uint64_t readsDeleted = 0;
    for (const auto& status : validReadIntervals) {
        if (status.isDeleted) {
            ++readsDeleted;
        }
    }

    cout << timestamp << "Hanging overlaps complete: " << readsDeleted
         << " total reads marked deleted (no surviving overlaps)." << endl;
}

/*
deleteContainmentOverlaps: Delete overlaps where one read is contained in the other.

Runs early (e.g. after computeBaseAlignmentsAndStore) to remove containment overlaps
without removing the contained reads themselves. Uses ma_hit2arc to detect
QCONT (query contained in target) or TCONT (target contained in query), and marks
those overlaps with DeleteReasonContained.

NOTE: ad.qs/qe/ts/te are the TIGHT, real CIGAR-alignment span, not extended --
ma_hit2arc needs hifiasm's ma_hit_t convention (diagonally extrapolated to read
boundaries). Use ad.extendedQs/extendedQe/extendedTs/extendedTe instead (see
AlignmentData::extendedQs doc comment) if/when this function is wired back in --
it is currently not called anywhere in the live pipeline.
*/
void Assembler::deleteContainmentOverlaps(uint64_t threadCount)
{
    cout << timestamp << "Deleting containment overlaps." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    uint64_t containmentBefore = 0;
    for (const auto& ad : alignmentData) {
        if ((ad.deleteReasons0 & AlignmentData::DeleteReasonContained) ||
            (ad.deleteReasons1 & AlignmentData::DeleteReasonContained)) {
            ++containmentBefore;
        }
    }

    setupLoadBalancing(alignmentData.size(), 10000);
    runThreads(&Assembler::deleteContainmentOverlapsThreadFunction, threadCount);

    uint64_t containmentAfter = 0;
    for (const auto& ad : alignmentData) {
        if ((ad.deleteReasons0 & AlignmentData::DeleteReasonContained) ||
            (ad.deleteReasons1 & AlignmentData::DeleteReasonContained)) {
            ++containmentAfter;
        }
    }

    cout << timestamp << "Containment overlaps: " << (containmentAfter - containmentBefore) << " deleted ("
         << containmentAfter << " total)." << endl;
}

void Assembler::deleteContainmentOverlapsThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    uint64_t begin, end;

    while (getNextBatch(begin, end)) {
        for (uint64_t i = begin; i != end; i++) {
            AlignmentData& ad = alignmentData[i];
            if (!ad.keptByBothSides()) continue;

            ReadId qn = ad.readIds[0];
            ReadId tn = ad.readIds[1];
            const uint32_t qLen = uint32_t(reads->getReadRawSequenceLength(qn));
            const uint32_t tLen = uint32_t(reads->getReadRawSequenceLength(tn));

            const uint32_t qs = ad.qs, qe = ad.qe, ts = ad.ts, te = ad.te;
            if (qe <= qs || te <= ts) continue;

            const int result = ma_hit2arc(
                (int32_t)qs, (int32_t)qe, (int32_t)qLen,
                (int32_t)ts, (int32_t)te, (int32_t)tLen,
                !ad.isSameStrand,
                INT32_MAX, 0.0f, 0);

            // MA_HT_QCONT or MA_HT_TCONT: delete this overlap
            if (result == MA_HT_QCONT || result == MA_HT_TCONT) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonContained);
            }
        }
    }
}

/*
ma_hit_cut + ma_hit_flt (hifiasm) — combined port.

Classifies each overlap using hifiasm's ma_hit2arc geometry test and deletes
overlaps that are internal matches or too short.

ma_hit2arc examines the "overhang" on each side of the overlap:

    5' overhang                              3' overhang
    |<------->|                              |<------->|
    query:  ===========XXXXXXXXXXXXXXXXX==============
                        |||||||||||||||||
    target:       ======XXXXXXXXXXXXXXXXX=======
                  |<--->|               |<----->|
                  tl5 (target 5' hang)  tl3 (target 3' hang)

    ext5 = min(qs, tl5)    — the smaller of the two 5' overhangs
    ext3 = min(ql-qe, tl3) — the smaller of the two 3' overhangs

Classification:
  MA_HT_INT (-1):       ext5 or ext3 exceeds maxHang, or the aligned portion
                         is less than maxHangRate of the total span.  These are
                         "internal" matches — the overlap sits in the middle of
                         both reads with large unaligned flanks on both sides.
  MA_HT_QCONT (-2):     query is contained in target (qs <= tl5 and ql-qe <= tl3).
  MA_HT_TCONT (-3):     target is contained in query.
  MA_HT_SHORT_OVLP (-4): effective overlap length < minOverlapLength.
  >= 0:                  proper dovetail overlap (value = non-overlap node length).

In the post-phasing pipeline coverage trimming is disabled (minCoverage=0),
so there is no per-read coordinate clipping and classification uses raw read
lengths.  We therefore use raw read lengths directly, skipping the clipping
step.

COORDINATES: this uses ad.qs/qe/ts/te, the TIGHT, real CIGAR-alignment span.
That is deliberate and required. Internal matches are only detectable on the
tight span: extending the coordinates toward the read tips first (as
extendOverlapToReadBoundaries does) snaps the smaller overhang on each side to
0, forcing ext5 = ext3 = 0, so ma_hit2arc could then never return MA_HT_INT.
With dinara's minCoverage=0, ma_hit_cut is a no-op and ma_hit_flt classifies on
the real overlap span against raw read lengths, matching the tight span used
here. This is the variant wired into the live pipeline (main.cpp).

Parameters (hifiasm defaults):
  maxHang          = 1000  (asm_opt.max_hang_Len)
  maxHangRate      = 0.8   (asm_opt.max_hang_rate)
  minOverlapLength = 50    (asm_opt.min_overlap_Len)
*/
void Assembler::deleteInternalOverlaps(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Deleting internal overlaps (maxHang=" << maxHang
         << ", maxHangRate=" << maxHangRate << ", minOverlapLength=" << minOverlapLength << ")." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    hangingFilterMaxHang = maxHang;
    hangingFilterMaxHangRate = maxHangRate;
    hangingFilterMinOverlap = minOverlapLength;

    uint64_t deletedBefore = 0;
    for (const auto& ad : alignmentData) {
        if ((ad.deleteReasons0 & AlignmentData::DeleteReasonHanging) ||
            (ad.deleteReasons1 & AlignmentData::DeleteReasonHanging)) {
            ++deletedBefore;
        }
    }

    setupLoadBalancing(alignmentData.size(), 10000);
    runThreads(&Assembler::deleteInternalOverlapsThreadFunction, threadCount);

    uint64_t deletedAfter = 0;
    for (const auto& ad : alignmentData) {
        if ((ad.deleteReasons0 & AlignmentData::DeleteReasonHanging) ||
            (ad.deleteReasons1 & AlignmentData::DeleteReasonHanging)) {
            ++deletedAfter;
        }
    }

    cout << timestamp << "Internal overlaps: " << (deletedAfter - deletedBefore) << " deleted ("
         << deletedAfter << " total)." << endl;
}

void Assembler::deleteInternalOverlapsThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    uint64_t begin, end;
    const uint64_t maxHang = this->hangingFilterMaxHang;
    const double maxHangRate = this->hangingFilterMaxHangRate;
    const uint64_t minOvlp = this->hangingFilterMinOverlap;

    while (getNextBatch(begin, end)) {
        for (uint64_t i = begin; i != end; i++) {
            AlignmentData& ad = alignmentData[i];

            // Skip overlaps already deleted by earlier stages (chimeric, weak, etc.).
            // Matches hifiasm's "if(p->del) continue" in ma_hit_cut/ma_hit_flt.
            if (!ad.keptByBothSides()) continue;

            const ReadId queryId  = ad.readIds[0];
            const ReadId targetId = ad.readIds[1];
            const uint32_t queryLen  = uint32_t(reads->getReadRawSequenceLength(queryId));
            const uint32_t targetLen = uint32_t(reads->getReadRawSequenceLength(targetId));

            // Tight CIGAR span -- see the header comment above for why the
            // extended coordinates would defeat internal-match detection.
            const uint32_t queryStart  = ad.qs;
            const uint32_t queryEnd    = ad.qe;
            const uint32_t targetStart = ad.ts;
            const uint32_t targetEnd   = ad.te;

            // Degenerate overlap: zero or negative length.
            if (queryEnd <= queryStart || targetEnd <= targetStart) continue;

            // ---------------------------------------------------------------
            // ma_hit2arc classification.
            //
            // With minCoverage=0, validReadIntervals are [0, readLen) for all
            // reads, so ma_hit_cut clipping is a no-op and ma_hit_flt uses
            // raw read lengths.  We pass raw lengths directly.
            //
            // Return values:
            //   >= 0:            dovetail overlap (keep)
            //   MA_HT_QCONT:    query contained in target (keep for now,
            //                    handled by flagContainedReads)
            //   MA_HT_TCONT:    target contained in query (keep for now)
            //   MA_HT_INT:      internal match — excessive overhangs (delete)
            //   MA_HT_SHORT_OVLP: overlap too short (delete)
            // ---------------------------------------------------------------
            const int classification = ma_hit2arc(
                (int32_t)queryStart,  (int32_t)queryEnd,  (int32_t)queryLen,
                (int32_t)targetStart, (int32_t)targetEnd, (int32_t)targetLen,
                !ad.isSameStrand,
                (int32_t)maxHang,
                maxHangRate,
                (int32_t)minOvlp);

            if (classification == MA_HT_INT || classification == MA_HT_SHORT_OVLP) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonHanging);
            }
        }
    }
}

/*
filterHangingOverlapsThreadFunction: alternative thread function that uses
coverage-trimmed read lengths (validReadIntervals) instead of raw lengths.

This is the full ma_hit_flt port for use when coverage trimming is enabled
(non-trivial validReadIntervals).  In the current post-phasing pipeline
coverage trimming is disabled, so deleteInternalOverlapsThreadFunction (which
uses raw lengths) is equivalent and is the one wired in.

Kept for future use if coverage trimming is enabled.
*/
void Assembler::filterHangingOverlapsThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    uint64_t begin, end;
    const uint64_t maxHang = this->hangingFilterMaxHang;
    const double maxHangRate = this->hangingFilterMaxHangRate;
    const uint64_t minOvlp = this->hangingFilterMinOverlap;

    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; i++) {
            AlignmentData& ad = alignmentData[i];

            // Skip overlaps already deleted by earlier stages.
            if(!ad.keptByBothSides()) continue;

            const ReadId queryId  = ad.readIds[0];
            const ReadId targetId = ad.readIds[1];

            // Use coverage-trimmed read lengths.  After ma_hit_cut clips
            // overlap coordinates to validReadIntervals, the effective read
            // length for ma_hit2arc is the trimmed interval length, not the
            // raw sequence length.
            uint32_t queryEffectiveLen, targetEffectiveLen;
            if (validReadIntervals.empty()) {
                queryEffectiveLen  = (uint32_t)reads->getReadRawSequenceLength(queryId);
                targetEffectiveLen = (uint32_t)reads->getReadRawSequenceLength(targetId);
            } else {
                const auto& queryInterval  = validReadIntervals[queryId];
                const auto& targetInterval = validReadIntervals[targetId];

                // Skip overlaps involving reads deleted by coverage trimming
                // or chimera detection.
                if (queryInterval.isDeleted || targetInterval.isDeleted) continue;

                queryEffectiveLen  = queryInterval.end - queryInterval.start;
                targetEffectiveLen = targetInterval.end - targetInterval.start;
            }

            uint32_t queryStart  = ad.qs;
            uint32_t queryEnd    = ad.qe;
            uint32_t targetStart = ad.ts;
            uint32_t targetEnd   = ad.te;

            if (queryEnd <= queryStart || targetEnd <= targetStart) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonHanging);
                continue;
            }

            const int classification = ma_hit2arc(
                (int32_t)queryStart,  (int32_t)queryEnd,  (int32_t)queryEffectiveLen,
                (int32_t)targetStart, (int32_t)targetEnd, (int32_t)targetEffectiveLen,
                !ad.isSameStrand,
                (int32_t)maxHang,
                maxHangRate,
                (int32_t)minOvlp
            );

            if (classification == MA_HT_INT || classification == MA_HT_SHORT_OVLP) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonHanging);
            }
        }
    }
}


/*
Rescue pass for chimeric detection.

This mirrors hifiasm's \"try_rescue_overlaps\" idea: for each chimeric read, if enough non-chimeric
neighbors overlap a common consensus interval, we unmark the read and enforce that only overlaps
spanning the consensus interval are kept for that read.
*/
void Assembler::rescueChimericReads(uint64_t threadCount)
{
    cout << timestamp << "Rescuing chimeric reads..." << endl;
    
    setupLoadBalancing(reads->readCount(), 1);
    runThreads(&Assembler::rescueChimericReadsThreadFunction, threadCount);

    uint64_t chimericCount = 0;
    for(size_t i=0; i<isChimericRead.size(); i++) {
        if(isChimericRead[i]) chimericCount++;
    }
    cout << timestamp << "Chimeric reads after rescue: " << chimericCount << "." << endl;
}

void Assembler::rescueChimericReadsThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    const uint32_t rescueThreshold = 4;
    
    uint64_t readIdBegin, readIdEnd;
    while(getNextBatch(readIdBegin, readIdEnd)) {
        for(ReadId readId = ReadId(readIdBegin); readId < ReadId(readIdEnd); readId++) {
            
            if (!isChimericRead[readId]) continue;
            
            OrientedReadId orientedReadId(readId, 0);
            
            if (orientedReadId.getValue() >= alignmentTable.size()) continue;
            
            vector<pair<uint32_t, int>> events;
            
            const size_t n = alignmentTable.size(orientedReadId.getValue());
            
            for(size_t i=0; i<n; i++) {
                const uint64_t alignmentId = alignmentTable[orientedReadId.getValue()][i];
                const AlignmentData& ad = alignmentData[alignmentId];
                
                ReadId neighborId;
                uint32_t qs, qe;
                
                if (ad.readIds[0] == readId) {
                    neighborId = ad.readIds[1];
                    qs = ad.qs; qe = ad.qe;
                } else {
                    neighborId = ad.readIds[0];
                    qs = ad.ts; qe = ad.te;
                }
                
                if (isChimericRead[neighborId]) continue;
                
                const CisTransStatus cts = ad.getCisTransStatusFromReadPerspective(readId);
                bool isValid = (cts == CisTransStatus::Cis);
                if (cts == CisTransStatus::Unknown && ad.info.isInReadGraph) {
                    isValid = true;
                }
                
                if (isValid) {
                    if (qs < qe) {
                        events.push_back({qs, 1});
                        events.push_back({qe, 2});
                    }
                }
            }
            
            if (events.size() / 2 >= rescueThreshold) {
                std::sort(events.begin(), events.end());
                
                int max_dp = 0;
                int dp = 0;
                uint32_t start = 0;
                uint32_t max_interval_s = 0;
                uint32_t max_interval_e = 0;
                int old_dp = 0;

                for(const auto& ev : events) {
                    old_dp = dp;
                    
                    if (ev.second == 1) {
                        dp++;
                    } else {
                        dp--;
                    }
                    
                    if (old_dp < dp) {
                        if (dp >= max_dp) {
                            start = ev.first;
                            max_dp = dp;
                        }
                    } else if (old_dp > dp) {
                        if (old_dp == max_dp) {
                            max_interval_s = start;
                            max_interval_e = ev.first;
                        }
                    }
                }
                
                if (max_dp >= (int)rescueThreshold) {
                    isChimericRead[readId] = false;
                    
                    for(size_t i=0; i<n; i++) {
                        const uint64_t alignmentId = alignmentTable[orientedReadId.getValue()][i];
                        AlignmentData& ad = alignmentData[alignmentId];
                        
                        ReadId neighborId;
                        uint32_t qs, qe;
                        if (ad.readIds[0] == readId) {
                            neighborId = ad.readIds[1];
                            qs = ad.qs; qe = ad.qe;
                        } else {
                            neighborId = ad.readIds[0];
                            qs = ad.ts; qe = ad.te;
                        }
                        if (isChimericRead[neighborId]) continue;
                        
                        const CisTransStatus cts = ad.getCisTransStatusFromReadPerspective(readId);
                        bool isValidType = (cts == CisTransStatus::Cis);
                        if (cts == CisTransStatus::Unknown && ad.info.isInReadGraph) isValidType = true;
                        
                        bool keep = false;
                        if (isValidType) {
                            if (qs <= max_interval_s && qe >= max_interval_e) {
                                keep = true;
                            }
                        }
                        
                        if (!keep) {
                            ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonChimeric);
                        }
                    }
                }
            }
        }
    }
}

/*
============================================================================
try_rescue_overlaps (hifiasm) parity - COMPLETE PARITY VERIFIED
============================================================================

This function resolves directional phasing conflicts where error correction or
phasing analysis has caused one read to reject an overlap while the other read
still accepts it.

HIFIASM ALGORITHM (try_rescue_overlaps):
----------------------------------------
Hifiasm's phasing pipeline (HiCanu integration, trio binning, etc.) marks overlaps
with a "phase deletion" flag when one read's phasing evidence suggests the overlap
should be rejected. However, if many overlaps involving a read agree on a consensus
region, those conflicts may be false positives that should be rescued.

For each read r:
  1. Find all overlaps with directional conflicts (DeleteReasonPhase set on exactly one side)
  2. Among those where THIS read deleted the overlap, collect the overlap intervals
  3. Use sweep-line to find the interval with maximum conflict depth
  4. If max_depth >= rescueThreshold, rescue conflicts spanning that consensus interval
  5. Rescue = clear DeleteReasonPhase flag from this read's perspective only

DIRECTIONAL CONFLICT DEFINITION:
---------------------------------
An overlap has a directional phasing conflict when:
  - Exactly one side has DeleteReasonPhase set (phase0 XOR phase1)
  - This indicates asymmetric phasing evidence: one read rejected it, the other didn't

KEY IMPLEMENTATION DETAILS FOR HIFIASM PARITY:
-----------------------------------------------
1. **Directional conflict detection** (lines 1379-1383):
   - Check: phase0 != phase1 (XOR condition)
   - Only process overlaps where exactly one side has DeleteReasonPhase

2. **Marker-based coordinates** (lines 1388-1389):
   - Tries to use getRawBoundsIfAvailable to get original marker bounds
   - Falls back to stored coordinates if markers unavailable
   - This ensures consistent interval boundaries for consensus detection

3. **Per-read conflict collection** (lines 1391-1404):
   - Only collects conflicts where THIS read deleted the overlap
   - Stores alignmentId and overlap interval on this read

4. **Threshold check** (line 1407):
   - Skip reads with fewer than rescueThreshold conflicts
   - Typical hifiasm default: rescueThreshold = 4

5. **Sweep-line consensus** (lines 1409-1440):
   - Build events: +1 for interval start, -1 for interval end
   - Sort by position, with starts before ends at same position
   - Track maximum depth interval where most conflicts agree

6. **Rescue decision** (lines 1442-1452):
   - Only rescue if max_dp >= rescueThreshold
   - Rescue overlaps spanning consensus interval [max_interval_s, max_interval_e]
   - Clear DeleteReasonPhase flag from THIS read's perspective only

OUTPUT:
-------
For rescued overlaps:
  - DeleteReasonPhase cleared from the read's perspective that originally deleted it
  - Other side's DeleteReasonPhase (if any) remains unchanged
  - Overlap becomes "kept by both sides" if no other delete reasons exist

This asymmetric rescue approach preserves the directional nature of phasing decisions
while correcting systematic false positives.
*/

// Count active alignments (kept by both sides)
uint64_t Assembler::countActiveAlignments() const
{
    uint64_t active = 0;
    for(uint64_t i = 0; i < alignmentData.size(); i++) {
        if(alignmentData[i].keptByBothSides()) active++;
    }
    return active;
}

void Assembler::rescuePhasedOverlaps(uint64_t rescueThreshold, uint64_t threadCount)
{
    cout << timestamp << "Rescuing phased overlaps (try_rescue_overlaps equivalent, rescueThreshold="
         << rescueThreshold << ")." << endl;

    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    this->rescuePhasedThreshold = rescueThreshold;

    setupLoadBalancing(readCount, 1);
    runThreads(&Assembler::rescuePhasedOverlapsThreadFunction, threadCount);

    // Compute diagnostics: count remaining directional phasing conflicts
    uint64_t directionalConflicts = 0;
    for (const auto& ad : alignmentData) {
        const bool phase0 = (ad.deleteReasons0 & AlignmentData::DeleteReasonPhase) != 0;
        const bool phase1 = (ad.deleteReasons1 & AlignmentData::DeleteReasonPhase) != 0;
        if (phase0 != phase1) {
            ++directionalConflicts;
        }
        // Note: we can't easily count "rescued" overlaps here since we don't track
        // before/after state, but directionalConflicts shows remaining conflicts
    }

    cout << timestamp << "Phased overlap rescue (try_rescue_overlaps): "
         << directionalConflicts << " directional conflicts remaining, "
         << countActiveAlignments() << " active overlaps." << endl;
}

void Assembler::rescuePhasedOverlapsThreadFunction(size_t threadId)
{
    /*
    ============================================================================
    try_rescue_overlaps (hifiasm) thread function - COMPLETE PARITY VERIFIED
    ============================================================================

    For each read, this function:
    1. Identifies directional phasing conflicts (DeleteReasonPhase on exactly one side)
    2. Collects conflicts where THIS read deleted the overlap
    3. Finds the consensus interval with maximum conflict depth
    4. Rescues conflicts spanning the consensus interval

    HIFIASM ALGORITHM:
    ------------------
    Hifiasm's try_rescue_overlaps aims to correct systematic false positives in
    phasing decisions. When many overlaps involving a read show directional conflicts
    that agree on a common interval, those conflicts are likely spurious and should
    be rescued.

    The algorithm processes each read independently and makes asymmetric corrections:
    clearing DeleteReasonPhase only from the perspective of the read that originally
    deleted the overlap.
    */
    static_cast<void>(threadId);
    const uint64_t rescueThreshold = this->rescuePhasedThreshold;

    uint64_t readIdBegin, readIdEnd;
    while(getNextBatch(readIdBegin, readIdEnd)) {
        for(ReadId readId = ReadId(readIdBegin); readId < ReadId(readIdEnd); readId++) {

            OrientedReadId orientedReadId(readId, 0);

            if (orientedReadId.getValue() >= alignmentTable.size()) continue;

            /*
            STEP 1: Collect directional phasing conflicts for this read
            -----------------------------------------------------------
            A directional conflict occurs when:
              - DeleteReasonPhase is set on exactly one side (phase0 XOR phase1)
              - THIS read is the side that deleted the overlap (thisReadPhaseDeleted)

            We collect the alignmentId and overlap interval on this read for each conflict.
            */
            std::vector<uint32_t> conflictAlignments;
            std::vector<std::pair<uint32_t, uint32_t>> conflictIntervals;

            /*
            HELPER: getRawBoundsIfAvailable
            --------------------------------
            Attempts to reconstruct the original overlap interval from marker positions.
            This is used to get more accurate interval boundaries for consensus detection,
            since stored coordinates may have been clipped by earlier filtering stages.

            Returns true if marker-based bounds are available, false otherwise.
            Falls back to stored coordinates if markers unavailable.
            */
            auto getRawBoundsIfAvailable = [&](const AlignmentData& ad, uint32_t& outQs, uint32_t& outQe, uint32_t& outTs, uint32_t& outTe) -> bool {
                if (!markers || !markers->isOpen()) return false;
                const auto& d0 = ad.info.data[0];
                const auto& d1 = ad.info.data[1];
                if (d0.markerCount == 0 || d1.markerCount == 0) return false;
                if (d0.firstOrdinal > d0.lastOrdinal || d1.firstOrdinal > d1.lastOrdinal) return false;

                const OrientedReadId oid0(ad.readIds[0], 0);
                const OrientedReadId oid1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
                if (oid0.getValue() >= markers->size() || oid1.getValue() >= markers->size()) return false;
                const auto& m0 = (*markers)[oid0.getValue()];
                const auto& m1 = (*markers)[oid1.getValue()];
                if (d0.lastOrdinal >= m0.size() || d1.lastOrdinal >= m1.size()) return false;

                const uint32_t k = uint32_t(assemblerInfo->k);
                const uint32_t qLen = uint32_t(reads->getReadRawSequenceLength(ad.readIds[0]));
                const uint32_t tLen = uint32_t(reads->getReadRawSequenceLength(ad.readIds[1]));

                const uint32_t qs0 = m0[d0.firstOrdinal].position;
                const uint32_t qe0 = m0[d0.lastOrdinal].position + k;
                const uint32_t tsOriented = m1[d1.firstOrdinal].position;
                const uint32_t teOriented = m1[d1.lastOrdinal].position + k;

                uint32_t ts0 = tsOriented;
                uint32_t te0 = teOriented;
                if (!ad.isSameStrand) {
                    const auto p = dinara::rcIntervalToForward(tLen, tsOriented, teOriented);
                    ts0 = p.first;
                    te0 = p.second;
                }

                if (qs0 >= qe0 || qe0 > qLen) return false;
                if (ts0 >= te0 || te0 > tLen) return false;

                outQs = qs0; outQe = qe0; outTs = ts0; outTe = te0;
                return true;
            };


            const size_t n = alignmentTable.size(orientedReadId.getValue());
            for(size_t i = 0; i < n; i++) {
                const uint64_t alignmentId = alignmentTable[orientedReadId.getValue()][i];
                const AlignmentData& ad = alignmentData[alignmentId];

                /*
                HIFIASM DIRECTIONAL CONFLICT CHECK:
                -----------------------------------
                Only process overlaps where exactly one side has DeleteReasonPhase set.
                This indicates a directional conflict:
                  - One read rejected the overlap due to phasing
                  - The other read did not reject it

                If both sides rejected (phase0 && phase1) or neither rejected (!phase0 && !phase1),
                skip this overlap as it's not a directional conflict.
                */
                const bool phase0 = (ad.deleteReasons0 & AlignmentData::DeleteReasonPhase) != 0;
                const bool phase1 = (ad.deleteReasons1 & AlignmentData::DeleteReasonPhase) != 0;
                if (phase0 == phase1) continue;  // Not a directional conflict

                uint32_t qs, qe;
                bool thisReadPhaseDeleted;

                /*
                Get overlap interval on THIS read:
                ----------------------------------
                Try to use marker-based bounds (more accurate), fall back to stored coordinates.
                Extract the interval [qs, qe) on readId and check if THIS read deleted it.
                */
                uint32_t qs0 = ad.qs, qe0 = ad.qe, ts0 = ad.ts, te0 = ad.te;
                (void)getRawBoundsIfAvailable(ad, qs0, qe0, ts0, te0);

                if (ad.readIds[0] == readId) {
                    qs = qs0;
                    qe = qe0;
                    thisReadPhaseDeleted = phase0;
                } else {
                    qs = ts0;
                    qe = te0;
                    thisReadPhaseDeleted = phase1;
                }

                /*
                Collect conflicts where THIS read deleted the overlap:
                ------------------------------------------------------
                Hifiasm try_rescue_overlaps only rescues from the perspective of the read
                that originally deleted the overlap. This maintains asymmetry in the rescue
                process.
                */
                if (thisReadPhaseDeleted) {
                    conflictAlignments.push_back((uint32_t)alignmentId);
                    conflictIntervals.push_back({qs, qe});
                }
            }


            /*
            HIFIASM THRESHOLD CHECK:
            ------------------------
            Only attempt rescue if this read has at least rescueThreshold conflicts.
            Typical hifiasm default: rescueThreshold = 4.

            This prevents spurious rescues when a read has only a few conflicts.
            */
            if (conflictAlignments.size() < rescueThreshold) continue;

            /*
            STEP 2: Build sweep-line events for consensus detection
            --------------------------------------------------------
            Hifiasm try_rescue_overlaps uses a sweep-line algorithm to find the interval
            where most conflicts agree. This identifies a "consensus region" that likely
            represents a systematic false positive in phasing.

            Events:
              +1 = interval start (conflict coverage increases)
              -1 = interval end (conflict coverage decreases)

            Sort by position, with starts (+1) before ends (-1) at same position.
            */
            std::vector<std::pair<uint32_t, int>> events;
            events.reserve(conflictIntervals.size() * 2);
            for(const auto& interval : conflictIntervals) {
                events.push_back({interval.first, 1});   // +1 = start
                events.push_back({interval.second, -1}); // -1 = end
            }

            // Sort by position, then by event type (start before end at same position)
            std::sort(events.begin(), events.end(), [](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first < b.first;
                return a.second > b.second; // +1 (start) before -1 (end)
            });

            /*
            STEP 3: Sweep through events to find maximum-depth consensus interval
            ----------------------------------------------------------------------
            Hifiasm try_rescue_overlaps finds the interval where the most conflicts overlap.
            This represents the consensus region where phasing evidence likely failed.

            Algorithm:
              - Track depth as we process events (+1 for start, -1 for end)
              - When depth increases to a new maximum, record the start position
              - When depth decreases from the maximum, record the interval [start, current_pos)
              - Keep the interval with maximum depth
            */
            int dp = 0, max_dp = 0;
            uint32_t start = 0;
            uint32_t max_interval_s = 0, max_interval_e = 0;

            for(const auto& ev : events) {
                int old_dp = dp;
                dp += ev.second;  // Update depth (+1 or -1)

                // Depth increased: potentially starting a new max-depth region
                if (old_dp < dp) {
                    if (dp >= max_dp) {
                        start = ev.first;  // Record potential start of max-depth region
                        max_dp = dp;
                    }
                }
                // Depth decreased: potentially ending the max-depth region
                else if (old_dp > dp) {
                    if (old_dp == max_dp) {
                        max_interval_s = start;      // Consensus interval start
                        max_interval_e = ev.first;   // Consensus interval end
                    }
                }
            }

            /*
            STEP 4: Rescue conflicts spanning the consensus interval
            ---------------------------------------------------------
            Hifiasm try_rescue_overlaps rescues conflicts only if:
              1. max_dp >= rescueThreshold (enough conflicts agree)
              2. The conflict spans the consensus interval [max_interval_s, max_interval_e)

            Rescue = clear DeleteReasonPhase flag from THIS read's perspective only.

            This asymmetric rescue preserves the directional nature of phasing decisions
            while correcting systematic false positives.
            */
            if ((uint64_t)max_dp >= rescueThreshold) {
                for(size_t j = 0; j < conflictAlignments.size(); j++) {
                    uint32_t alignmentId = conflictAlignments[j];
                    uint32_t qs = conflictIntervals[j].first;
                    uint32_t qe = conflictIntervals[j].second;

                    // Rescue overlaps that span the consensus interval
                    if (qs <= max_interval_s && qe >= max_interval_e) {
                        AlignmentData& ad = alignmentData[alignmentId];
                        // Clear DeleteReasonPhase from THIS read's perspective only
                        ad.clearDeleteReasonsFromReadPerspective(readId, AlignmentData::DeleteReasonPhase);
                    }
                }
            }
        }
    }
}

// ma_hit2arc is now in overlapClassification.hpp (shared with clique filter).


void Assembler::removeReadsFlaggedContained(uint64_t threadCount)
{
    cout << timestamp << "Removing reads flagged as contained..." << endl;

    checkAlignmentDataAreOpen();
    reads->checkReadsAreOpen();
    reads->checkReadFlagsAreOpen();

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Ensure validReadIntervals exists. If missing, treat entire read as valid/alive.
    if (validReadIntervals.empty()) {
        cout << timestamp
             << "[DIAG] removeReadsFlaggedContained: validReadIntervals not set; using full read lengths."
             << endl;
        validReadIntervals.resize(reads->readCount());
        for (ReadId r = 0; r < reads->readCount(); ++r) {
            const uint64_t len = reads->getReadRawSequenceLength(r);
            validReadIntervals[r] = {0, uint32_t(len), false};
        }
    }

    if (!containmentParent->isOpen) {
        containmentParent->createNew(largeDataName("ContainmentParent"), largeDataPageSize);
        containmentParent->resize(reads->readCount());
        std::fill(containmentParent->begin(), containmentParent->end(), ReadId(invalidReadId));
    }

    auto deleteAllEdgesForRead = [&](ReadId r) {
        const OrientedReadId oid(r, 0);
        if (oid.getValue() >= alignmentTable.size()) return;
        const auto& table = alignmentTable[oid.getValue()];
        for (const uint32_t alignmentId : table) {
            alignmentData[alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonContained);
        }
    };

    const uint64_t readCount = reads->readCount();
    uint64_t containedReadCount = 0;
    for (ReadId r = 0; r < readCount; ++r) {
        if (r >= validReadIntervals.size() || validReadIntervals[r].isDeleted) continue;
        if (!reads->getFlags(r).isContained) continue;

        validReadIntervals[r].isDeleted = true;
        deleteAllEdgesForRead(r);
        ++containedReadCount;
    }

    // Compress containment chains (cycle-safe).
    for (ReadId r = 0; r < readCount; ++r) {
        if ((*containmentParent)[r] == ReadId(invalidReadId)) continue;
        ReadId root = r;
        vector<ReadId> visited;
        visited.reserve(16);
        while (root != ReadId(invalidReadId) && (*containmentParent)[root] != ReadId(invalidReadId)) {
            if (std::find(visited.begin(), visited.end(), root) != visited.end()) {
                ReadId cycleRoot = root;
                for (const ReadId x : visited) {
                    if (x < cycleRoot) cycleRoot = x;
                }
                (*containmentParent)[r] = cycleRoot;
                root = ReadId(invalidReadId);
                break;
            }
            visited.push_back(root);
            root = (*containmentParent)[root];
        }
        if (root != ReadId(invalidReadId)) {
            (*containmentParent)[r] = root;
        }
    }

    setupLoadBalancing(readCount, 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);

    cout << timestamp << "Removed " << containedReadCount << " reads flagged as contained." << endl;
}

/*
ma_hit_contained_advance (hifiasm) — detection logic ported, handling differs.

Identifies reads that are fully contained within another read using the same
ma_hit2arc geometry test as deleteInternalOverlaps.  A read is "contained" when
its overlap with another read covers the entire shorter read:

    query:    ====XXXXXXXXXXXXXXXXX====
                  |||||||||||||||||
    target:       XXXXXXXXXXXXXXXXX        ← target contained in query (TCONT)

    query:        XXXXXXXXXXXXXXXXX        ← query contained in target (QCONT)
                  |||||||||||||||||
    target:  ====XXXXXXXXXXXXXXXXX====

Containment is detected by ma_hit2arc when both overhangs on one side are
smaller than the corresponding overhangs on the other side (qs <= tl5 and
ql-qe <= tl3 for QCONT, or the reverse for TCONT).

This must run sequentially: when a read is found to be contained, it is
immediately masked out so later reads don't consider overlaps to it.  This
matches hifiasm's delete_all_edges + coverage_cut[qn].del = 1 behavior.

Architectural difference from hifiasm:
  hifiasm eagerly deletes all overlaps of contained reads (delete_all_edges)
  and removes them from the overlap graph entirely.  dinara only sets the
  isContained flag on the read — overlaps are preserved.  Downstream stages
  (pruneContainedReadsToOneBestOverlapByDpScore, marker graph construction)
  use this flag to handle contained reads with more flexibility.

Parameters (hifiasm defaults):
  maxHang          = 1000  (asm_opt.max_hang_Len)
  maxHangRate      = 0.8   (asm_opt.max_hang_rate)
  minOverlapLength = 50    (asm_opt.min_overlap_Len)
*/
void Assembler::flagContainedReads(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "flagContainedReads begins." << endl;

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    checkAlignmentDataAreOpen();
    reads->checkReadsAreOpen();
    reads->checkReadFlagsAreOpenForWriting();

    const uint64_t readCount = reads->readCount();
    const uint64_t alignmentCount = alignmentData.size();

    for(ReadId r = 0; r < readCount; ++r) {
        reads->setContainedFlag(r, false);
    }

    // Local deletion mask.  Once a read is deemed contained, we set its bit
    // so subsequent iterations skip it and skip overlaps incident to it.
    // This replicates hifiasm's immediate delete_all_edges +
    // coverage_cut[qn].del = 1 mutation visibility.
    vector<uint8_t> containedMask(readCount, 0);

    uint64_t containedReadCount = 0;

    for(ReadId queryId = 0; queryId < readCount; ++queryId) {
        if(containedMask[queryId]) continue;

        const int32_t queryLen = int32_t(reads->getReadRawSequenceLength(queryId));
        if(queryLen <= 0) continue;

        const OrientedReadId oid(queryId, 0);
        if(oid.getValue() >= alignmentTable.size()) continue;

        const auto& overlaps = alignmentTable[oid.getValue()];
        for(const uint32_t alignmentId : overlaps) {
            if(alignmentId >= alignmentCount) continue;
            const AlignmentData& ad = alignmentData[alignmentId];

            // Skip overlaps already deleted by earlier stages.
            if(!ad.keptByBothSides()) continue;

            const ReadId targetId = (ad.readIds[0] == queryId) ? ad.readIds[1] : ad.readIds[0];
            if(containedMask[targetId]) continue;

            const int32_t targetLen = int32_t(reads->getReadRawSequenceLength(targetId));
            if(targetLen <= 0) continue;

            // Project overlap coordinates onto queryId's axis. Use the
            // EXTENDED (hifiasm ma_hit_t convention) coordinates, not the
            // tight qs/qe/ts/te -- ma_hit2arc's containment test assumes
            // coordinates diagonally extrapolated to read boundaries (see
            // AlignmentData::extendedQs doc comment).
            int32_t queryStart, queryEnd, targetStart, targetEnd;
            if(ad.readIds[0] == queryId) {
                queryStart  = int32_t(ad.extendedQs);
                queryEnd    = int32_t(ad.extendedQe);
                targetStart = int32_t(ad.extendedTs);
                targetEnd   = int32_t(ad.extendedTe);
            } else {
                queryStart  = int32_t(ad.extendedTs);
                queryEnd    = int32_t(ad.extendedTe);
                targetStart = int32_t(ad.extendedQs);
                targetEnd   = int32_t(ad.extendedQe);
            }
            if(queryStart >= queryEnd || targetStart >= targetEnd) continue;

            const int classification = ma_hit2arc(
                queryStart, queryEnd, queryLen,
                targetStart, targetEnd, targetLen,
                !ad.isSameStrand,
                int32_t(maxHang),
                maxHangRate,
                int32_t(minOverlapLength)
            );

            if(classification == MA_HT_QCONT) {
                // Query is contained in target.  Flag it and stop scanning
                // this query's overlaps (hifiasm breaks after delete_all_edges).
                reads->setContainedFlag(queryId, true);
                containedMask[queryId] = 1;
                ++containedReadCount;
                break;
            } else if(classification == MA_HT_TCONT) {
                // Target is contained in query.  Flag it, keep scanning.
                if(!containedMask[targetId]) {
                    reads->setContainedFlag(targetId, true);
                    containedMask[targetId] = 1;
                    ++containedReadCount;
                }
            }
        }
    }

    cout << timestamp << "flagContainedReads: " << containedReadCount << " contained reads out of "
         << readCount << " total." << endl;
}



void Assembler::pruneContainedReadsToOneBestOverlapByDpScore(uint64_t /* threadCount */)
{
    checkAlignmentDataAreOpen();
    reads->checkReadFlagsAreOpen();

    const uint64_t readCount = reads->readCount();
    const uint64_t alignmentCount = alignmentData.size();

    // If a read graph has already been built in this process, prefer pruning only within
    // the current read-graph overlap set. This avoids selecting an overlap that was not
    // part of the broad read graph used to build the marker graph, which could otherwise
    // disconnect a contained read after pruning.
    bool restrictToCurrentReadGraph = false;
    for (uint64_t alignmentId = 0; alignmentId < alignmentCount; ++alignmentId) {
        if (alignmentData[alignmentId].info.isInReadGraph != 0) {
            restrictToCurrentReadGraph = true;
            break;
        }
    }

    uint64_t containedReadCount = 0;
    uint64_t containedReadsWithKeptOverlap = 0;
    uint64_t containedReadsWithNoKeptOverlap = 0;
    uint64_t containedReadsKeptToContainmentParent = 0;
    uint64_t containedReadsFellBackToAnyPartner = 0;
    uint64_t containedReadsMissingContainmentParent = 0;
    uint64_t overlapsPruned = 0;

    const bool hasContainmentParent =
        containmentParent &&
        containmentParent->isOpen &&
        containmentParent->size() >= readCount;

    // Deterministic pseudo-random tie-breaking per read.
    // This avoids changing results depending on iteration order while still breaking ties "randomly".
    for(ReadId r=0; r<readCount; ++r) {
        if(!reads->getFlags(r).isContained) {
            continue;
        }
        ++containedReadCount;

        const ReadId containmentRoot =
            hasContainmentParent ? (*containmentParent)[r] : ReadId(invalidReadId);
        if (!hasContainmentParent) {
            ++containedReadsMissingContainmentParent;
        }

        const OrientedReadId oid0(r, 0);
        if(oid0.getValue() >= alignmentTable.size()) {
            ++containedReadsWithNoKeptOverlap;
            continue;
        }

        const span<uint32_t> table = alignmentTable[oid0.getValue()];
        if(table.empty()) {
            ++containedReadsWithNoKeptOverlap;
            continue;
        }

        // Prefer an active overlap to the containment parent/root if available.
        int64_t bestScoreToParent = std::numeric_limits<int64_t>::min();
        vector<uint32_t> bestAlignmentIdsToParent;
        bestAlignmentIdsToParent.reserve(2);

        // Fallback: best hifiasm shared_seed proxy among all active overlaps for this contained read.
        int64_t bestScoreAny = std::numeric_limits<int64_t>::min();
        vector<uint32_t> bestAlignmentIdsAny;
        bestAlignmentIdsAny.reserve(4);

        for(const uint32_t alignmentId : table) {
            if(alignmentId >= alignmentCount) {
                continue;
            }
            const AlignmentData& ad = alignmentData[alignmentId];
            if (restrictToCurrentReadGraph && (ad.info.isInReadGraph == 0)) {
                continue;
            }
            if(!ad.keptByBothSides()) {
                continue;
            }
            const int64_t score = ad.info.hifiasmSharedSeedScoreProxy();

            // Update best-overlap-to-parent (if applicable).
            if (containmentRoot != ReadId(invalidReadId)) {
                const ReadId other = (ad.readIds[0] == r) ? ad.readIds[1] : ad.readIds[0];
                if (other == containmentRoot) {
                    if(score > bestScoreToParent) {
                        bestScoreToParent = score;
                        bestAlignmentIdsToParent.clear();
                        bestAlignmentIdsToParent.push_back(alignmentId);
                    } else if(score == bestScoreToParent) {
                        bestAlignmentIdsToParent.push_back(alignmentId);
                    }
                }
            }

            // Update best-overlap-to-any partner.
            if(score > bestScoreAny) {
                bestScoreAny = score;
                bestAlignmentIdsAny.clear();
                bestAlignmentIdsAny.push_back(alignmentId);
            } else if(score == bestScoreAny) {
                bestAlignmentIdsAny.push_back(alignmentId);
            }
        }

        const vector<uint32_t>& candidateBest =
            bestAlignmentIdsToParent.empty() ? bestAlignmentIdsAny : bestAlignmentIdsToParent;

        if(candidateBest.empty()) {
            ++containedReadsWithNoKeptOverlap;
            continue;
        }
        ++containedReadsWithKeptOverlap;
        if (!bestAlignmentIdsToParent.empty()) {
            ++containedReadsKeptToContainmentParent;
        } else if (containmentRoot != ReadId(invalidReadId)) {
            ++containedReadsFellBackToAnyPartner;
        }

        // Choose randomly among ties (deterministic per read).
        std::mt19937_64 rng(uint64_t(r) * 0x9e3779b97f4a7c15ULL + 0xD1A2B3C4ULL);
        std::uniform_int_distribution<size_t> dist(0, candidateBest.size() - 1);
        const uint32_t chosenAlignmentId = candidateBest[dist(rng)];

        // Prune all other overlaps incident to this contained read.
        for(const uint32_t alignmentId : table) {
            if(alignmentId >= alignmentCount) {
                continue;
            }
            if(alignmentId == chosenAlignmentId) {
                continue;
            }

            AlignmentData& ad = alignmentData[alignmentId];
            if (restrictToCurrentReadGraph && (ad.info.isInReadGraph == 0)) {
                continue;
            }
            // This diagnostic pruning is intended to remove overlaps from the read graph, so
            // we mark the overlap deleted from both sides (AND semantics).
            ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonContainedPrune);
            ++overlapsPruned;
        }
    }

    cout << timestamp << "[DIAG] pruneContainedReadsToOneBestOverlapByDpScore:"
         << " containedReads=" << containedReadCount
         << " kept=" << containedReadsWithKeptOverlap
         << " noneKept=" << containedReadsWithNoKeptOverlap
         << " keptToContainmentParent=" << containedReadsKeptToContainmentParent
         << " fellBackToAnyPartner=" << containedReadsFellBackToAnyPartner
         << " missingContainmentParent=" << containedReadsMissingContainmentParent
         << " overlapsPruned=" << overlapsPruned << endl;
}

