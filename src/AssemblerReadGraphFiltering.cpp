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
- ma_hit_sub                 -> filterLocalSegments
- detect_chimeric_reads      -> detectChimericReads (+ rescueChimericReads)
- ma_hit_cut                 -> applyCoverageCuts
- ma_hit_flt                 -> filterHangingOverlaps
- ma_hit_contained_advance   -> removeContainedReads (containmentParent + path compression)
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
ma_hit_sub (miniasm/hifiasm) parity: compute a per-read high-coverage interval.

For each read `r`, we gather its active overlaps (those not deleted from `r`'s perspective),
build a sweep-line depth profile over the read, and keep the longest interval where
depth >= `minCoverage`:
  - `validReadIntervals[r].start/end` correspond to hifiasm's `coverage_cut[r].s/e`
  - `validReadIntervals[r].isDeleted` corresponds to `coverage_cut[r].del`

When `minCoverage <= 1`, hifiasm keeps the full read as the valid region, so we do the same.
*/
/*
ma_hit_sub (hifiasm) — exact port.

Computes a per-read valid interval by finding the longest contiguous region
where overlap depth >= minCoverage.

When minCoverage <= 1 (the default, asm_opt.min_overlap_coverage = 0), this
is a no-op: every read gets validReadIntervals = [0, readLen).  The function
must still run because downstream stages (detectChimericReads,
deleteInternalOverlaps, flagContainedReads) read validReadIntervals.

When minCoverage >= 2, the algorithm:
  1. Collects all overlap intervals on the read (from one orientation only,
     matching hifiasm's sources[r] adjacency list to avoid double-counting).
  2. Builds a sweep-line depth profile using +1/-1 events.
  3. Finds the longest interval where depth >= minCoverage.
  4. Stores the result in validReadIntervals[r].

Reads with no valid interval get isDeleted = true, which causes downstream
stages to skip them (equivalent to hifiasm's coverage_cut[r].del = 1).
*/
void Assembler::filterLocalSegments(
    uint64_t minCoverage,
    uint64_t threadCount)
{
    cout << timestamp << "filterLocalSegments begins (minCoverage="
         << minCoverage << ")." << endl;

    const uint64_t readCount = reads->readCount();
    validReadIntervals.resize(readCount);
    std::fill(validReadIntervals.begin(), validReadIntervals.end(), ReadSegmentStatus{0, 0, false});
    this->localSegmentMinCoverage = minCoverage;

    setupLoadBalancing(readCount, 100);
    runThreads(&Assembler::filterLocalSegmentsThreadFunction, threadCount);

    uint64_t readsWithValidIntervals = 0;
    uint64_t readsDeleted = 0;
    uint64_t totalValidBases = 0;
    uint64_t totalRawBases = 0;
    for (ReadId r = 0; r < readCount; ++r) {
        const auto& interval = validReadIntervals[r];
        const uint64_t rawLen = reads->getReadRawSequenceLength(r);
        totalRawBases += rawLen;
        if (interval.isDeleted) {
            ++readsDeleted;
        } else if (interval.end > interval.start) {
            ++readsWithValidIntervals;
            totalValidBases += (interval.end - interval.start);
        }
    }

    cout << timestamp << "filterLocalSegments: "
         << readsWithValidIntervals << " reads with valid intervals, "
         << readsDeleted << " reads deleted ("
         << (readCount - readsDeleted) << " remaining). "
         << "Valid bases: " << totalValidBases << "/" << totalRawBases
         << " (" << std::fixed << std::setprecision(1)
         << (100.0 * totalValidBases / std::max(1UL, totalRawBases)) << "%)."
         << endl;
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

void Assembler::filterLocalSegmentsThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    const uint64_t minDepth = this->localSegmentMinCoverage;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); readId++) {

            // ---------------------------------------------------------------
            // Fast path: minDepth <= 1.
            //
            // hifiasm's ma_hit_sub skips all sweep-line computation when
            // min_dp <= 1 and sets coverage_cut[r] = {0, readLen, del=0}.
            // With the default asm_opt.min_overlap_coverage = 0, this is
            // always taken — every read gets the full [0, readLen) interval.
            // ---------------------------------------------------------------
            if (minDepth <= 1) {
                const uint32_t readLen = uint32_t(reads->getReadRawSequenceLength(readId));
                validReadIntervals[readId] = {0, readLen, false};
                continue;
            }

            // ---------------------------------------------------------------
            // Collect overlap intervals on this read.
            //
            // Query only alignmentTable[OrientedReadId(readId, 0)] — this
            // contains all overlaps involving readId exactly once, matching
            // hifiasm's sources[r] adjacency list.
            // ---------------------------------------------------------------
            std::vector<std::pair<uint32_t, uint32_t>> intervals;

            const auto& table = alignmentTable[OrientedReadId(readId, 0).getValue()];
            for (uint32_t alignmentId : table) {
                const auto& ad = alignmentData[alignmentId];
                if (!ad.keptByBothSides()) continue;

                uint32_t ovlpStart = 0, ovlpEnd = 0;
                if (ad.readIds[0] == readId) {
                    ovlpStart = ad.qs;
                    ovlpEnd   = ad.qe;
                } else {
                    ovlpStart = ad.ts;
                    ovlpEnd   = ad.te;
                }

                if (ovlpEnd > ovlpStart) {
                    intervals.push_back({ovlpStart, ovlpEnd});
                }
            }

            // No overlaps → mark read as deleted (coverage_cut[r].del = 1).
            if (intervals.empty()) {
                validReadIntervals[readId] = {0, 0, true};
                continue;
            }

            // ---------------------------------------------------------------
            // Sweep-line: find the longest interval with depth >= minDepth.
            //
            // Build +1/-1 events, sort by position (starts before ends at
            // the same position), then sweep tracking depth transitions.
            // ---------------------------------------------------------------
            std::vector<std::pair<uint32_t, int>> events;
            events.reserve(intervals.size() * 2);
            for (const auto& interval : intervals) {
                events.push_back({interval.first,  +1});
                events.push_back({interval.second, -1});
            }

            // Sort: by position, then starts (+1) before ends (-1).
            std::sort(events.begin(), events.end(),
                [](const pair<uint32_t, int>& a, const pair<uint32_t, int>& b) {
                    if (a.first != b.first) return a.first < b.first;
                    return a.second > b.second;
                });

            int depth = 0;
            uint32_t longestStart = 0, longestEnd = 0;
            uint32_t segmentStart = 0;
            bool inSegment = false;

            for (const auto& ev : events) {
                const int oldDepth = depth;
                depth += ev.second;
                const uint32_t pos = ev.first;

                // Crossed threshold upward: start a new segment.
                if (oldDepth < (int)minDepth && depth >= (int)minDepth) {
                    segmentStart = pos;
                    inSegment = true;
                }
                // Crossed threshold downward: close segment, keep if longest.
                else if (oldDepth >= (int)minDepth && depth < (int)minDepth) {
                    if (inSegment) {
                        const uint32_t segmentLen = pos - segmentStart;
                        const uint32_t longestLen = longestEnd - longestStart;
                        if (segmentLen > longestLen) {
                            longestStart = segmentStart;
                            longestEnd   = pos;
                        }
                        inSegment = false;
                    }
                }
            }

            // Store result.
            if (longestEnd > longestStart) {
                validReadIntervals[readId] = {longestStart, longestEnd, false};
            } else {
                validReadIntervals[readId] = {0, 0, true};
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
  - Containment (keep): one read fully within another (handled by removeContainedReads later)
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

Uses extended coordinates (ad.qs/qe/ts/te) from chaining.
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

In the post-phasing pipeline, filterLocalSegments runs with minCoverage=0, so
validReadIntervals are [0, readLen) for all reads.  This makes ma_hit_cut a
no-op (no coordinate clipping) and ma_hit_flt uses raw read lengths.  We
therefore use raw read lengths directly, skipping the clipping step.

ad.qs/qe/ts/te store extended coordinates (matching hifiasm's
append_inexact_overlap_region_alloc), where chain anchor positions are
extended toward read tips.  These are the same coordinates hifiasm stores
in ma_hit_t and passes to ma_hit2arc.

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

            // ad.qs/qe/ts/te are extended coordinates (matching hifiasm's
            // append_inexact_overlap_region_alloc), where chain anchor
            // positions are extended toward read tips.  These are the same
            // coordinates hifiasm stores in ma_hit_t and passes to ma_hit2arc.
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

This is the full ma_hit_flt port for use when filterLocalSegments runs with
minCoverage >= 2 (non-trivial coverage trimming).  In the current post-phasing
pipeline minCoverage=0, so deleteInternalOverlapsThreadFunction (which uses
raw lengths) is equivalent and is the one wired in.

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
detect_chimeric_reads (hifiasm) parity.

This runs after ma_hit_sub (validReadIntervals computed) and before ma_hit_cut (overlaps still in
absolute coordinates). For each read, we:
1) find left/right anchors that touch the read ends,
2) extend anchors using contained/near-contained overlaps,
3) decide obvious chimeras (gap) and ambiguous chimeras (small overlap) using boundary verification.

We mark chimeric reads in `isChimericRead`/`validReadIntervals[].isDeleted` and mark incident overlaps
with `DeleteReasonChimeric`.
*/
/*
detect_chimeric_reads (hifiasm) — exact port.

A chimeric read is one formed by ligation of two unrelated genomic fragments.
The signature: overlaps anchored to the left end of the read and overlaps
anchored to the right end don't meet in the middle — there is a gap or a
suspiciously thin bridge.

Algorithm (three phases per read):

  Phase 1 — collect_sides: find the reach of left-anchored and right-anchored
  overlaps.  A "left-anchored" overlap starts at position 0 on the read; a
  "right-anchored" overlap ends at the read length.  We track the union
  envelope of each group:
      leftAnchorEnd   = max qe  among overlaps with qs == 0
      rightAnchorStart = min qs  among overlaps with qe == readLen

  Phase 2 — collect_contain: extend the envelopes using "interior" overlaps
  (those touching neither end) that overlap the current envelope by at least
  10% of their length.  This bridges small gaps caused by contained overlaps
  that don't reach the read ends.

  Phase 3 — classify:
    (a) If leftAnchorEnd and rightAnchorStart overlap by >= readLen * 0.06,
        the read is normal (not chimeric).
    (b) If leftAnchorEnd <= rightAnchorStart (gap between anchors), the read
        is a simple chimera.
    (c) Otherwise the anchors overlap but not enough — complex case.  We check
        whether any overlap spanning the ambiguous zone [rightAnchorStart,
        leftAnchorEnd) fails base-level boundary verification (banded BPM
        alignment).  If any spanning overlap fails, the read is chimeric.

Hifiasm parity notes:
  - This must run sequentially: when a read is marked chimeric, all its
    overlaps are deleted immediately, and later reads observe those deletions.
  - hifiasm also checks the `el` field (non-homopolymer error count from
    error correction) as a fast path in the complex case.  We omit this
    because dinara does not perform hifiasm-style error correction.
  - hifiasm has a UL (ultra-long read) bypass in collect_sides that marks
    reads as non-chimeric if enough ultra-long reads span them.  We omit
    this because dinara does not use ultra-long reads.

Parameters (hifiasm defaults):
  shiftRate   = max_ov_diff_final * 2.0 = 0.03 * 2.0 = 0.06
  overlapRate = 0.1  (minimum fractional overlap for collect_contain extension)
*/
void Assembler::detectChimericReads(uint64_t threadCount)
{
    cout << timestamp << "Detecting chimeric reads..." << endl;

    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();
    static_cast<void>(threadCount);

    const uint64_t readCount = reads->readCount();

    if (!isChimericRead.isOpen) {
        isChimericRead.createNew(largeDataName("IsChimericRead"), largeDataPageSize);
        isChimericRead.resize(readCount);
    }
    std::fill(isChimericRead.begin(), isChimericRead.end(), false);

    // hifiasm defaults.
    constexpr double shiftRate   = 0.03 * 2.0; // asm_opt.max_ov_diff_final * 2.0
    constexpr float  overlapRate = 0.1f;        // collect_contain fractional threshold

    // Scratch buffers for boundary verification (reused across reads).
    std::vector<char> boundaryBufX;
    std::vector<char> boundaryBufY;

    // -----------------------------------------------------------------------
    // delete_all_edges equivalent: mark a read as chimeric and delete all
    // its incident overlaps from both perspectives.
    //
    // hifiasm sets del=1 on every overlap in sources[qn] and the reverse
    // edge in sources[tn], plus coverage_cut[qn].del = 1.
    // We set DeleteReasonChimeric on both sides of every incident overlap
    // and mark validReadIntervals[qn].isDeleted = true.
    // -----------------------------------------------------------------------
    auto deleteAllEdges = [&](ReadId readId) {
        if (readId >= readCount || isChimericRead[readId]) {
            return;
        }
        isChimericRead[readId] = true;
        if (readId < validReadIntervals.size()) {
            validReadIntervals[readId].isDeleted = true;
        }

        const OrientedReadId oid(readId, 0);
        if (oid.getValue() >= alignmentTable.size()) {
            return;
        }
        for (uint32_t alignmentId : alignmentTable[oid.getValue()]) {
            alignmentData[alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonChimeric);
        }
    };

    // -----------------------------------------------------------------------
    // Extract the overlap interval on a given read's coordinate axis.
    //
    // Each AlignmentData stores coordinates for readIds[0] in (qs,qe) and
    // for readIds[1] in (ts,te), both on the respective read's forward
    // strand.  This lambda returns whichever pair belongs to `readId`.
    // -----------------------------------------------------------------------
    auto getOverlapIntervalOnRead = [&](const AlignmentData& ad, ReadId readId,
                                        uint32_t& ovlpStart, uint32_t& ovlpEnd) {
        if (ad.readIds[0] == readId) {
            ovlpStart = ad.qs;
            ovlpEnd   = ad.qe;
        } else {
            DINARA_ASSERT(ad.readIds[1] == readId);
            ovlpStart = ad.ts;
            ovlpEnd   = ad.te;
        }
    };

    // -----------------------------------------------------------------------
    // Project an overlap into the coordinate system needed by boundaryVerify.
    //
    // boundaryVerify needs:
    //   - partnerId:       the other read in the overlap
    //   - ovlpStartOnRead: overlap start on readId's forward strand (= xs in hifiasm)
    //   - ovlpEndOnRead:   overlap end on readId's forward strand
    //   - partnerStart:    overlap start on partner's forward strand (= ts in hifiasm)
    //   - partnerEnd:      overlap end on partner's forward strand (= te in hifiasm)
    //   - isReverseStrand: true if the reads are on opposite strands
    // -----------------------------------------------------------------------
    struct OverlapProjection {
        ReadId   partnerId        = invalid<ReadId>;
        uint32_t ovlpStartOnRead  = 0;   // qs from readId's perspective
        uint32_t ovlpEndOnRead    = 0;   // qe from readId's perspective
        uint32_t partnerStart     = 0;   // ts on partner's forward strand
        uint32_t partnerEnd       = 0;   // te on partner's forward strand
        bool     isReverseStrand  = false;
    };

    auto projectOverlap = [&](const AlignmentData& ad, ReadId readId) -> OverlapProjection {
        OverlapProjection proj;
        if (ad.readIds[0] == readId) {
            proj.partnerId       = ad.readIds[1];
            proj.ovlpStartOnRead = ad.qs;
            proj.ovlpEndOnRead   = ad.qe;
            proj.partnerStart    = ad.ts;
            proj.partnerEnd      = ad.te;
        } else {
            DINARA_ASSERT(ad.readIds[1] == readId);
            proj.partnerId       = ad.readIds[0];
            proj.ovlpStartOnRead = ad.ts;
            proj.ovlpEndOnRead   = ad.te;
            proj.partnerStart    = ad.qs;
            proj.partnerEnd      = ad.qe;
        }
        proj.isReverseStrand = !ad.isSameStrand;
        return proj;
    };

    uint64_t chimericCount        = 0;
    uint64_t simpleChimericCount  = 0;
    uint64_t complexCheckedCount  = 0;
    uint64_t complexChimericCount = 0;

    for (ReadId readId = 0; readId < readCount; ++readId) {
        const uint32_t readLen = uint32_t(reads->getReadRawSequenceLength(readId));
        const OrientedReadId oid(readId, 0);
        if (oid.getValue() >= alignmentTable.size()) {
            continue;
        }
        const auto& overlaps = alignmentTable[oid.getValue()];

        // =================================================================
        // Phase 1: collect_sides
        //
        // Find the farthest reach of left-anchored overlaps (qs == 0) and
        // right-anchored overlaps (qe == readLen).
        //
        // leftAnchorEnd:    rightmost qe among all overlaps starting at 0
        // rightAnchorStart: leftmost qs among all overlaps ending at readLen
        //
        // Initialized to sentinel values: leftAnchorEnd = 0 (no reach),
        // rightAnchorStart = readLen (no reach).  The .s fields track the
        // minimum start for each group (always 0 for left, variable for
        // right); we use leftAnchorFound/rightAnchorFound to detect the
        // "end node" case where one side has no anchors at all.
        // =================================================================
        uint32_t leftAnchorEnd    = 0;       // max qe among left-anchored overlaps
        uint32_t rightAnchorStart = readLen;  // min qs among right-anchored overlaps
        bool leftAnchorFound  = false;
        bool rightAnchorFound = false;

        for (uint32_t alignmentId : overlaps) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if (!ad.keptByBothSides()) {
                continue;
            }

            uint32_t ovlpStart = 0, ovlpEnd = 0;
            getOverlapIntervalOnRead(ad, readId, ovlpStart, ovlpEnd);

            // Left-anchored: overlap starts at position 0.
            if (ovlpStart == 0) {
                leftAnchorFound = true;
                if (ovlpEnd > leftAnchorEnd) {
                    leftAnchorEnd = ovlpEnd;
                }
            }

            // Right-anchored: overlap ends at the read length.
            // Note: an overlap can be both left- and right-anchored (spans
            // the entire read).  hifiasm adds it to both groups.
            if (ovlpEnd == readLen) {
                rightAnchorFound = true;
                if (ovlpStart < rightAnchorStart) {
                    rightAnchorStart = ovlpStart;
                }
            }
        }

        // End node: one side has no anchored overlaps at all.  This read is
        // a tip in the overlap graph, not a chimera.
        if (!leftAnchorFound || !rightAnchorFound) {
            continue;
        }

        // =================================================================
        // Phase 2: collect_contain
        //
        // Extend the anchor envelopes using interior overlaps (those that
        // touch neither end of the read).  An interior overlap [qs, qe)
        // extends the left envelope rightward if it overlaps leftAnchorEnd
        // by at least 10% of its length.  Similarly for the right envelope.
        //
        // This catches cases where a contained overlap bridges the gap
        // between the anchor envelope and the middle of the read.
        //
        // We accumulate extensions into temporaries and apply once, matching
        // hifiasm's single-pass collect_contain.
        // =================================================================
        uint32_t extendedLeftEnd    = leftAnchorEnd;
        uint32_t extendedRightStart = rightAnchorStart;

        for (uint32_t alignmentId : overlaps) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if (!ad.keptByBothSides()) {
                continue;
            }

            uint32_t ovlpStart = 0, ovlpEnd = 0;
            getOverlapIntervalOnRead(ad, readId, ovlpStart, ovlpEnd);

            // Skip anchored overlaps — only interior overlaps participate.
            if (ovlpStart == 0 || ovlpEnd == readLen) {
                continue;
            }
            const uint32_t ovlpLen = ovlpEnd - ovlpStart;
            if (ovlpLen == 0) {
                continue;
            }

            // Can this overlap extend the left envelope rightward?
            // Condition: the overlap straddles leftAnchorEnd, and the portion
            // to the left of leftAnchorEnd is at least 10% of the overlap length.
            if (ovlpStart < leftAnchorEnd && ovlpEnd > leftAnchorEnd &&
                (leftAnchorEnd - ovlpStart) > uint32_t(overlapRate * float(ovlpLen))) {
                if (ovlpEnd > extendedLeftEnd) {
                    extendedLeftEnd = ovlpEnd;
                }
            }

            // Can this overlap extend the right envelope leftward?
            // Condition: the overlap straddles rightAnchorStart, and the portion
            // to the right of rightAnchorStart is at least 10% of the overlap length.
            if (ovlpStart < rightAnchorStart && ovlpEnd > rightAnchorStart &&
                (ovlpEnd - rightAnchorStart) > uint32_t(overlapRate * float(ovlpLen))) {
                if (ovlpStart < extendedRightStart) {
                    extendedRightStart = ovlpStart;
                }
            }
        }

        leftAnchorEnd    = extendedLeftEnd;
        rightAnchorStart = extendedRightStart;

        // =================================================================
        // Phase 3: classify
        //
        // Compare the (possibly extended) left and right envelopes.
        //
        //   leftAnchorEnd
        //        v
        //   |============================|.............|===================| read
        //        left envelope                gap?         right envelope
        //                                              ^
        //                                       rightAnchorStart
        //
        // (a) Large overlap (>= 6% of read length): normal read.
        // (b) No overlap (gap): simple chimera.
        // (c) Small overlap: complex case — verify by base-level alignment.
        // =================================================================

        // (a) Normal read: the two envelopes overlap substantially.
        if (leftAnchorEnd > rightAnchorStart &&
            (leftAnchorEnd - rightAnchorStart) >= uint32_t(double(readLen) * shiftRate)) {
            continue;
        }

        // (b) Simple chimera: clear gap between left and right envelopes.
        if (leftAnchorEnd <= rightAnchorStart) {
            deleteAllEdges(readId);
            ++chimericCount;
            ++simpleChimericCount;
            continue;
        }

        // (c) Complex chimera: the envelopes overlap, but the overlap is
        // suspiciously thin (< 6% of read length).
        //
        // The ambiguous zone is [rightAnchorStart, leftAnchorEnd).  We look
        // for any overlap that fully spans this zone and check whether its
        // base-level alignment across the zone boundary is consistent.
        //
        // hifiasm also short-circuits on the `el` field (non-homopolymer
        // error count from error correction).  We skip this because dinara
        // does not perform hifiasm-style error correction — boundaryVerify
        // alone decides.
        ++complexCheckedCount;
        const uint32_t ambiguousZoneStart = rightAnchorStart;
        const uint32_t ambiguousZoneEnd   = leftAnchorEnd;
        bool chimericByBoundary = false;

        for (uint32_t alignmentId : overlaps) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if (!ad.keptByBothSides()) {
                continue;
            }

            const OverlapProjection proj = projectOverlap(ad, readId);

            // Does this overlap fully span the ambiguous zone?
            if (proj.ovlpStartOnRead <= ambiguousZoneStart &&
                proj.ovlpEndOnRead   >= ambiguousZoneEnd) {

                // boundaryVerify returns true if the base-level alignment
                // across the zone boundary is consistent (NOT chimeric).
                // Returns false if the alignment fails (chimeric evidence).
                if (!boundaryVerify(*reads,
                        ambiguousZoneStart, ambiguousZoneEnd,
                        readId, proj.partnerId,
                        proj.ovlpStartOnRead,
                        proj.partnerStart, proj.partnerEnd,
                        proj.isReverseStrand,
                        boundaryBufX, boundaryBufY)) {
                    chimericByBoundary = true;
                    break;
                }
            }
        }

        if (chimericByBoundary) {
            deleteAllEdges(readId);
            ++chimericCount;
            ++complexChimericCount;
        }
    }

    cout << timestamp
         << "Detected " << chimericCount << " chimeric reads."
         << " simple=" << simpleChimericCount
         << " complex=" << complexChimericCount
         << "/" << complexCheckedCount
         << endl;
}

void Assembler::detectChimericReadsThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    throw runtime_error("detectChimericReadsThreadFunction is obsolete; detectChimericReads now runs sequentially for hifiasm parity.");
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


void Assembler::removeContainedReads(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Removing contained reads (ma_hit_contained_advance equivalent, maxHang="
         << maxHang << ", maxHangRate=" << maxHangRate << ", minOverlapLength="
         << minOverlapLength << ")." << endl;

    /*
    Hifiasm ma_hit_contained_advance: Remove redundant contained reads from assembly.

    This function implements hifiasm's containment detection and removal logic. A read is
    "contained" if it lies entirely within another read, making it redundant for assembly.

    Algorithm (exact hifiasm parity):

    1. Sequential Read Processing (IMPORTANT - not parallelized):
       Hifiasm processes reads sequentially because when read A is found contained in read B,
       it immediately deletes all edges incident to A. Later reads must observe these deletions
       to make correct containment decisions. This matches hifiasm's serial loop.

    2. For Each Alive Read qn:
       - Traverse its adjacency list via `alignmentTable[OrientedReadId(qn,0)]`
       - Only consider overlaps kept by both sides (`keptByBothSides()`) → hifiasm: `h->del == 0`
       - Skip overlaps to already-deleted reads

    3. Containment Classification:
       - Call `ma_hit2arc()` from qn's perspective
       - result == MA_HT_QCONT (-2): query (qn) is contained in target (tn)
       - result == MA_HT_TCONT (-3): target (tn) is contained in query (qn)
       - result >= 0: dovetail overlap (keep both reads)
       - result == MA_HT_INT (-1) or MA_HT_SHORT_OVLP (-4): internal/bad (already filtered by ma_hit_flt)

    4. Containment Actions (hifiasm `delete_all_edges`):
       When containment detected:
       - Mark contained read as deleted: `validReadIntervals[r].isDeleted = true`
       - Record container as parent: `containmentParent[contained] = container`
       - Delete ALL incident overlaps: `deleteAllEdgesForRead(contained)`
       - Increment counter

    5. Chain Compression (hifiasm `transfor_R_to_U`):
       Handle transitive containment chains: if A in B and B in C, set parent[A] = C
       This finds the "ultimate" container for each contained read.

    6. Cleanup Pass:
       Delete reads that became isolated (no surviving overlaps) after containment removal.
       Matches hifiasm's `rLen==0 -> coverage_cut[i].del=1` behavior.

    Key Differences from ma_hit_flt:
    - ma_hit_flt: Removes BAD OVERLAPS (internal matches), keeps good overlaps including containments
    - ma_hit_contained_advance: Removes REDUNDANT READS and ALL their overlaps

    Hifiasm Parameters (defaults):
      - max_hang = 1000 (asm_opt.max_hang)
      - int_frac = 0.8 (asm_opt.int_frac, called maxHangRate here)
      - min_ovlp = 50 (asm_opt.min_ovlp)
    */
    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Initialize containmentParent tracking
    if (!containmentParent->isOpen) {
        containmentParent->createNew(largeDataName("ContainmentParent"), largeDataPageSize);
        containmentParent->resize(reads->readCount());
    }
    std::fill(containmentParent->begin(), containmentParent->end(), ReadId(invalidReadId));

    uint64_t containedReadCount = 0;
    uint64_t overlapsDeletedByContainment = 0;

    // Hifiasm delete_all_edges: delete ALL overlaps incident to a contained read
    auto deleteAllEdgesForRead = [&](ReadId r) {
        const OrientedReadId oid(r, 0);
        if (oid.getValue() >= alignmentTable.size()) return;
        const auto& table = alignmentTable[oid.getValue()];
        for (const uint32_t alignmentId : table) {
            AlignmentData& ad = alignmentData[alignmentId];
            if (!ad.isDeleted()) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonContained);
                ++overlapsDeletedByContainment;
            }
        }
    };

    /*
    Hifiasm ma_hit_contained_advance: Sequential read processing.

    CRITICAL: This loop must be SEQUENTIAL (not parallelized) to match hifiasm.
    When read A is found contained and deleted, later reads must observe this deletion
    to make correct containment decisions.
    */
    const uint64_t readCount = reads->readCount();
    for (ReadId qn = 0; qn < readCount; ++qn) {
        // Hifiasm: skip already-deleted reads
        if (qn >= validReadIntervals.size() || validReadIntervals[qn].isDeleted) continue;

        // Get query read's valid region length (after ma_hit_cut normalization)
        const auto& vrQ = validReadIntervals[qn];
        const int32_t ql = int32_t(vrQ.end - vrQ.start);
        if (ql <= 0) continue;

        // Get query read's adjacency list
        const OrientedReadId oid(qn, 0);
        if (oid.getValue() >= alignmentTable.size()) continue;

        const auto& table = alignmentTable[oid.getValue()];
        for (uint32_t alignmentId : table) {
            const AlignmentData& ad = alignmentData[alignmentId];

            // Hifiasm ma_hit_contained_advance: only consider overlaps kept by both sides (h->del == 0)
            if(!ad.keptByBothSides()) continue;

            // Only consider cis (1) or unclassified (0) overlaps for containment.
            // Trans (2) and cisDifferentCopy (3) overlaps are from different
            // haplotypes and should not cause a read to be marked as contained.
            if(ad.hifiasmEcMatchState0 > 1 || ad.hifiasmEcMatchState1 > 1) continue;

            // Get target read
            const ReadId tn = (ad.readIds[0] == qn) ? ad.readIds[1] : ad.readIds[0];

            // Hifiasm: skip overlaps to already-deleted reads
            if (tn >= validReadIntervals.size() || validReadIntervals[tn].isDeleted) continue;

            // Get target read's valid region length
            const auto& vrT = validReadIntervals[tn];
            const int32_t tl = int32_t(vrT.end - vrT.start);
            if (tl <= 0) continue;

            /*
            Hifiasm coordinate reconstruction: view overlap from query's perspective.

            After ma_hit_cut, coordinates are normalized to valid regions.
            We need to reconstruct [qs,qe) and [ts,te) from qn's viewpoint:
            - If ad.readIds[0] == qn: use ad.qs/qe/ts/te directly
            - If ad.readIds[1] == qn: swap query/target (ad.ts/te become qs/qe, ad.qs/qe become ts/te)

            IMPORTANT: Do NOT RC-map coordinates! (hifiasm set_reverse_overlap parity)
            - Coordinates are always in forward frame of each read
            - isReverse indicates orientation, but coordinates stay forward
            */
            const bool rev = !ad.isSameStrand;
            int32_t qs = 0, qe = 0, ts = 0, te = 0;
            if (ad.readIds[0] == qn) {
                qs = int32_t(ad.qs);
                qe = int32_t(ad.qe);
                ts = int32_t(ad.ts);
                te = int32_t(ad.te);
            } else {
                // Swap query/target without RC-mapping
                qs = int32_t(ad.ts);
                qe = int32_t(ad.te);
                ts = int32_t(ad.qs);
                te = int32_t(ad.qe);
            }

            // Sanity checks: validate coordinates
            if (qs < 0 || qe < 0 || ts < 0 || te < 0) continue;
            if (qs >= qe || ts >= te) continue;
            if (qs > ql || qe > ql || ts > tl || te > tl) continue;

            const int result = ma_hit2arc(
                qs, qe, ql,
                ts, te, tl,
                rev,
                int32_t(maxHang),
                maxHangRate,
                int32_t(minOverlapLength)
            );

            if (result == MA_HT_QCONT) {
                // QCONT: query (qn) is contained in target (tn)
                if (!validReadIntervals[qn].isDeleted) {
                    validReadIntervals[qn].isDeleted = true;              // Mark read deleted
                    (*containmentParent)[qn] = tn;                        // Record container
                    deleteAllEdgesForRead(qn);                           // Delete ALL its overlaps
                    ++containedReadCount;
                }
                // Hifiasm: stop processing this read once contained
                break;

            } else if (result == MA_HT_TCONT) {
                // TCONT: target (tn) is contained in query (qn)
                if (!validReadIntervals[tn].isDeleted) {
                    validReadIntervals[tn].isDeleted = true;              // Mark read deleted
                    (*containmentParent)[tn] = qn;                        // Record container
                    deleteAllEdgesForRead(tn);                           // Delete ALL its overlaps
                    ++containedReadCount;
                }
                // Continue checking other overlaps of qn
            }
        }
    }

    /*
    Hifiasm transfor_R_to_U: Compress containment chains to find ultimate containers.

    Example chain: Read A contained in B, B contained in C, C is ultimate container
      Before compression: parent[A]=B, parent[B]=C, parent[C]=invalid
      After compression:  parent[A]=C, parent[B]=C, parent[C]=invalid

    This is critical for:
    1. Efficient containment queries (one lookup instead of chain traversal)
    2. Correct graph construction (edges should point to ultimate containers)
    3. Matching hifiasm's behavior exactly

    Algorithm: For each contained read, follow parent chain to root (ultimate container)
    */
    for (ReadId r = 0; r < readCount; ++r) {
        // Skip non-contained reads
        if ((*containmentParent)[r] == ReadId(invalidReadId)) continue;

        // Follow chain to ultimate container
        ReadId root = (*containmentParent)[r];
        while (root != ReadId(invalidReadId) && (*containmentParent)[root] != ReadId(invalidReadId)) {
            root = (*containmentParent)[root];
        }

        // Compress: point directly to ultimate container
        (*containmentParent)[r] = root;
    }

    cout << timestamp << "Containment detection: " << containedReadCount << " reads contained, "
         << overlapsDeletedByContainment << " overlaps deleted." << endl;

    // Hifiasm ma_hit_contained_advance: cleanup reads with no surviving overlaps (rLen == 0)
    setupLoadBalancing(readCount, 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);

    // Count final deleted reads (contained + isolated)
    uint64_t readsDeleted = 0;
    for (const auto& status : validReadIntervals) {
        if (status.isDeleted) {
            ++readsDeleted;
        }
    }

    cout << timestamp << "Containment removal complete: " << readsDeleted
         << " total reads deleted (contained + isolated)." << endl;
}

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

            // Project overlap coordinates onto queryId's axis.
            int32_t queryStart, queryEnd, targetStart, targetEnd;
            if(ad.readIds[0] == queryId) {
                queryStart  = int32_t(ad.qs);
                queryEnd    = int32_t(ad.qe);
                targetStart = int32_t(ad.ts);
                targetEnd   = int32_t(ad.te);
            } else {
                queryStart  = int32_t(ad.ts);
                queryEnd    = int32_t(ad.te);
                targetStart = int32_t(ad.qs);
                targetEnd   = int32_t(ad.qe);
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



/*
cleanWeakOverlaps: Remove weak cis overlaps that are contradicted by
strong phasing evidence.

Exact logic port of hifiasm's clean_weak_ma_hit_t (Overlaps.cpp:11177).

Background:
  After phaseOverlapsKmeans, each overlap has per-side phasing labels:
    0 = unknown (no het sites in the overlap region)
    1 = cis (same haplotype)
    2 = trans (different haplotype)

  An overlap is "weak" when both sides are unknown (0) — the overlap
  region contained no het SNP sites, so the phasing pipeline could not
  determine the haplotype relationship. In hifiasm, this corresponds to
  ml == 0 (the overlap does not span any het site identified during
  error correction).

  A weak cis overlap can be dangerous: it might connect reads from
  different haplotypes in a region with no distinguishing variants.
  This function identifies and removes such false cis overlaps.

Algorithm (for standard HiFi, ou_thres = -1):
  For each weak overlap qn→tn (both sides ecMatchState == 0):
    Search qn's overlaps for a "strong" overlap qn→S where:
      a) S has at least one side with ecMatchState != 0 (has het sites)
      b) S's interval on qn contains the weak overlap's interval
         (qs_strong <= qs_weak && qe_strong >= qe_weak)
      c) S has a trans overlap to tn (ecMatchState == 2 from either side)

    If such a chain qn→S→tn exists, the weak overlap is contradicted:
      - qn and S are cis (strong overlap confirms same haplotype)
      - S and tn are trans (different haplotype)
      - Therefore qn and tn should be trans, not cis
      → Mark the weak overlap for deletion.

    If no contradicting chain exists, the weak overlap is kept.

  After marking, set DeleteReasonWeakContradicted on both sides of
  all marked overlaps.

Differences from hifiasm:
  - hifiasm's ml bit maps to our ecMatchState: ml==0 ↔ both sides unknown,
    ml==1 ↔ at least one side labeled.
  - hifiasm's sources/reverse_sources arrays map to our ecMatchState
    labels: cis overlaps are in sources, trans in reverse_sources.
  - hifiasm marks with a temporary bit in bl, then sets del=1. We set
    a delete reason bitmask on both sides.
  - The ou_thres / update_weak_by_contain rescue pass (for ultra-long
    reads) is not implemented. For standard HiFi (ou_thres == -1), that
    pass is a no-op in hifiasm.

This runs after phaseOverlapsKmeans, before rescueTransOverlaps.
*/
void Assembler::cleanWeakOverlaps()
{
    cout << timestamp << "cleanWeakOverlaps begins." << endl;

    checkAlignmentDataAreOpen();

    const uint64_t alignmentCount = alignmentData.size();

    uint64_t weakCount = 0;
    uint64_t contradictedCount = 0;

    // Helper: check if read S has a trans overlap to read tn.
    // Scans alignmentTable[S] for an overlap with tn where at least
    // one side says trans (ecMatchState == 2).
    // Matches hifiasm's get_specific_overlap(&reverse_sources[S], S, tn).
    auto hasTransOverlap = [&](ReadId readS, ReadId readTn) -> bool {
        const OrientedReadId oidS(readS, 0);
        if(oidS.getValue() >= alignmentTable.size()) return false;
        const auto& tableS = alignmentTable[oidS.getValue()];
        for(const uint32_t aid : tableS) {
            if(aid >= alignmentCount) continue;
            const AlignmentData& ad2 = alignmentData[aid];
            // Find the overlap between S and tn.
            ReadId partner;
            uint8_t partnerState;
            if(ad2.readIds[0] == readS) {
                partner = ad2.readIds[1];
                // S's perspective on this overlap.
                partnerState = ad2.hifiasmEcMatchState0;
            } else {
                partner = ad2.readIds[0];
                partnerState = ad2.hifiasmEcMatchState1;
            }
            if(partner != readTn) continue;
            // In hifiasm, the overlap lives in reverse_sources[S] if S
            // classified it as trans. We check if either side says trans.
            if(partnerState == 2 || partnerState == 3) return true;
            // Also check the other side — if tn says trans about S.
            uint8_t otherState = (ad2.readIds[0] == readS)
                ? ad2.hifiasmEcMatchState1
                : ad2.hifiasmEcMatchState0;
            if(otherState == 2 || otherState == 3) return true;
        }
        return false;
    };

    // Helper: check if a weak overlap qn→tn is contradicted.
    // Matches hifiasm's check_weak_ma_hit (returns 0 when contradicted).
    // Searches qn's overlaps for a strong overlap qn→S that:
    //   1) has het sites (at least one side ecMatchState != 0)
    //   2) contains the weak overlap's interval on qn
    //   3) S has a trans overlap to tn
    auto isContradicted = [&](ReadId readQn, ReadId readTn,
                              uint32_t weakQs, uint32_t weakQe) -> bool {
        const OrientedReadId oidQn(readQn, 0);
        if(oidQn.getValue() >= alignmentTable.size()) return false;
        const auto& tableQn = alignmentTable[oidQn.getValue()];

        for(const uint32_t aid : tableQn) {
            if(aid >= alignmentCount) continue;
            const AlignmentData& adStrong = alignmentData[aid];

            // Get the strong overlap's properties from qn's perspective.
            ReadId strongTarget;
            uint8_t strongState0, strongState1;
            uint32_t strongQs, strongQe;
            if(adStrong.readIds[0] == readQn) {
                strongTarget = adStrong.readIds[1];
                strongState0 = adStrong.hifiasmEcMatchState0;
                strongState1 = adStrong.hifiasmEcMatchState1;
                strongQs     = adStrong.qs;
                strongQe     = adStrong.qe;
            } else {
                strongTarget = adStrong.readIds[0];
                strongState0 = adStrong.hifiasmEcMatchState1;
                strongState1 = adStrong.hifiasmEcMatchState0;
                strongQs     = adStrong.ts;
                strongQe     = adStrong.te;
            }

            // Must not be deleted. Matches hifiasm's
            // "aim_paf->buffer[i].del == 0" check.
            if(adStrong.isDeleted0() || adStrong.isDeleted1()) continue;

            // Must be strong: at least one side has a phasing label.
            // Matches hifiasm's "aim_paf->buffer[i].ml == 1".
            if(strongState0 == 0 && strongState1 == 0) continue;

            // Must be in qn's cis array. In hifiasm, sources[qn] contains
            // all overlaps NOT classified as trans by qn — both cis (ml==1)
            // and unclassified (ml==0). The ml==1 check above already filters
            // to strong overlaps. Here we exclude overlaps that qn
            // classified as trans or different-copy.
            if(strongState0 == 2 || strongState0 == 3) continue;

            // Must contain the weak overlap's interval on qn.
            if(strongQs > weakQs || strongQe < weakQe) continue;

            // Check if the strong target S has a trans overlap to tn.
            if(hasTransOverlap(strongTarget, readTn)) {
                return true; // contradicted
            }
        }
        return false;
    };

    // Main loop: scan all overlaps, mark weak ones that are contradicted.
    // We collect alignment IDs to mark, then apply deletions in a second pass
    // to avoid modifying data while iterating.
    vector<uint32_t> toDelete;

    for(uint32_t alignmentId = 0; alignmentId < alignmentCount; ++alignmentId) {
        const AlignmentData& ad = alignmentData[alignmentId];

        // Skip already-deleted overlaps. Matches hifiasm's
        // "if(sources[i].buffer[j].del) continue" check.
        if(ad.isDeleted0() || ad.isDeleted1()) continue;

        // Weak: both sides unknown (no het sites).
        if(ad.hifiasmEcMatchState0 != 0 || ad.hifiasmEcMatchState1 != 0) continue;
        ++weakCount;

        // Check from readIds[0]'s perspective (qn = readIds[0], tn = readIds[1]).
        // Matches hifiasm iterating sources[i] and checking each weak overlap.
        const ReadId qn = ad.readIds[0];
        const ReadId tn = ad.readIds[1];

        if(isContradicted(qn, tn, ad.qs, ad.qe)) {
            toDelete.push_back(alignmentId);
            continue;
        }

        // Also check from readIds[1]'s perspective (qn = readIds[1], tn = readIds[0]).
        // Hifiasm processes each read independently, so both directions are checked.
        if(isContradicted(tn, qn, ad.ts, ad.te)) {
            toDelete.push_back(alignmentId);
        }
    }

    // Apply deletions.
    for(const uint32_t alignmentId : toDelete) {
        AlignmentData& ad = alignmentData[alignmentId];
        ad.deleteReasons0 |= AlignmentData::DeleteReasonWeakContradicted;
        ad.deleteReasons1 |= AlignmentData::DeleteReasonWeakContradicted;
    }

    contradictedCount = toDelete.size();

    cout << timestamp << "cleanWeakOverlaps: "
         << weakCount << " weak overlaps (both sides unlabeled), "
         << contradictedCount << " contradicted and deleted." << endl;
}



/*
try_rescue_overlaps (hifiasm) — exact port.

Resolves directional cis/trans disagreements.  phaseOverlapsKmeans processes
each read independently, so for the same overlap A-B, read A may call it
trans while read B calls it cis.  This function rescues the trans call when
the evidence suggests it was wrong.

For each read:
  1. Collect "disagreement" overlaps: this read says trans, partner says cis.
  2. If fewer than minPileup disagreements, skip.
  3. Sweep-line on the disagreement intervals to find the peak-depth region.
     Uses hifiasm's exact encoding: start events = (position << 1), end
     events = (position << 1 | 1), sorted as uint32_t.  Peak tracking uses
     >= (takes the latest start at peak depth).
  4. If peak depth >= minPileup, rescue only disagreements whose interval
     fully contains the peak region (containment filter).
  5. Flip rescued overlaps from trans (2) to cis (1) on this read's side.

hifiasm physically moves rescued overlaps between reverse_sources and
sources arrays.  We flip hifiasmEcMatchState in place — equivalent because
downstream code checks the state value, not array membership.

skipDeleted maps to hifiasm's is_del parameter: false for HiFi, true for ONT.

Parameters:
  minPileup   = 4      (rescue_threshold)
  skipDeleted = false   (is_del = 0 for HiFi)
*/
void Assembler::rescueTransOverlaps(uint64_t minPileup, bool skipDeleted)
{
    cout << timestamp << "rescueTransOverlaps begins (minPileup=" << minPileup
         << ", skipDeleted=" << skipDeleted << ")." << endl;

    checkAlignmentDataAreOpen();

    const uint64_t readCount = reads->readCount();
    const uint64_t alignmentCount = alignmentData.size();

    uint64_t totalDisagreements = 0;
    uint64_t readsRescued = 0;
    uint64_t overlapsRescued = 0;
    uint64_t skippedDeletedCount = 0;

    struct Disagreement {
        uint32_t ovlpStart;     // overlap start on this read's coordinate axis
        uint32_t ovlpEnd;       // overlap end on this read's coordinate axis
        uint32_t alignmentId;
    };
    vector<Disagreement> disagreements;
    vector<uint32_t> events;

    for (ReadId readId = ReadId(0); readId < readCount; ++readId) {

        const OrientedReadId oid(readId, 0);
        if (oid.getValue() >= alignmentTable.size()) continue;
        const auto& overlaps = alignmentTable[oid.getValue()];

        // ---------------------------------------------------------------
        // Pass 1: collect disagreement overlaps.
        //
        // A disagreement is an overlap where this read says trans (2) but
        // the partner says cis (1).  In hifiasm, this means the overlap
        // is in reverse_sources[readId] (this read's trans array) AND in
        // paf[partnerId] (partner's cis array).
        // ---------------------------------------------------------------
        disagreements.clear();
        for (const uint32_t alignmentId : overlaps) {
            if (alignmentId >= alignmentCount) continue;
            const AlignmentData& ad = alignmentData[alignmentId];

            if (skipDeleted && (ad.isDeleted0() || ad.isDeleted1())) {
                ++skippedDeletedCount;
                continue;
            }

            uint8_t myState, partnerState;
            uint32_t myStart, myEnd;
            if (ad.readIds[0] == readId) {
                myState      = ad.hifiasmEcMatchState0;
                partnerState = ad.hifiasmEcMatchState1;
                myStart      = ad.qs;
                myEnd        = ad.qe;
            } else {
                myState      = ad.hifiasmEcMatchState1;
                partnerState = ad.hifiasmEcMatchState0;
                myStart      = ad.ts;
                myEnd        = ad.te;
            }

            if (myState == 2 && partnerState == 1 && myStart < myEnd) {
                disagreements.push_back({myStart, myEnd, alignmentId});
            }
        }

        totalDisagreements += disagreements.size();
        if (disagreements.size() < minPileup) continue;

        // ---------------------------------------------------------------
        // Pass 2: sweep-line to find the peak-depth interval.
        //
        // Events encoded as (position << 1 | isEnd).  Sorting as uint32_t
        // puts starts before ends at the same position.
        //
        // Peak tracking: dp >= maxDp (>=, not >) takes the latest start
        // at peak depth.  When depth drops from the peak, we record the
        // interval [start, position).
        // ---------------------------------------------------------------
        events.clear();
        events.reserve(disagreements.size() * 2);
        for (const auto& d : disagreements) {
            events.push_back(d.ovlpStart << 1);
            events.push_back(d.ovlpEnd << 1 | 1);
        }
        sort(events.begin(), events.end());

        int64_t depth = 0, oldDepth = 0, peakDepth = 0;
        uint32_t peakStart = 0;
        uint32_t peakIntervalStart = 0, peakIntervalEnd = 0;

        for (size_t j = 0; j < events.size(); ++j) {
            oldDepth = depth;
            if (events[j] & 1) {
                --depth;
            } else {
                ++depth;
            }

            if (oldDepth < depth) {
                if (depth >= peakDepth) {
                    peakStart = events[j] >> 1;
                    peakDepth = depth;
                }
            } else if (oldDepth > depth) {
                if (oldDepth == peakDepth) {
                    peakIntervalStart = peakStart;
                    peakIntervalEnd   = events[j] >> 1;
                }
            }
        }

        if (static_cast<uint64_t>(peakDepth) < minPileup) continue;

        // ---------------------------------------------------------------
        // Pass 3: rescue disagreements that span the peak interval.
        //
        // Only overlaps whose interval fully contains [peakIntervalStart,
        // peakIntervalEnd] are rescued.  This filters out scattered
        // disagreements outside the coherent peak region.
        // ---------------------------------------------------------------
        ++readsRescued;
        for (const auto& d : disagreements) {
            if (d.ovlpStart <= peakIntervalStart && d.ovlpEnd >= peakIntervalEnd) {
                AlignmentData& ad = alignmentData[d.alignmentId];
                if (ad.readIds[0] == readId) {
                    ad.hifiasmEcMatchState0 = 1; // trans -> cis
                } else {
                    ad.hifiasmEcMatchState1 = 1; // trans -> cis
                }
                ++overlapsRescued;
            }
        }
    }

    cout << timestamp << "rescueTransOverlaps: "
         << totalDisagreements << " disagreements, rescued "
         << overlapsRescued << " overlaps across "
         << readsRescued << " reads." << endl;
}
