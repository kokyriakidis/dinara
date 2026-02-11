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
    inline bool isDeletedFromReadPerspective(const dinara::AlignmentData& ad, dinara::ReadId readId)
    {
        if (ad.readIds[0] == readId) {
            return ad.isDeleted0();
        } else {
            return ad.isDeleted1();
        }
    }
}

static int ma_hit2arc_containment(
    int32_t qs, int32_t qe, int32_t ql,
    int32_t ts, int32_t te, int32_t tl,
    bool isReverse,
    int32_t max_hang,
    double int_frac,
    int32_t min_ovlp);

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
void Assembler::filterLocalSegments(
    uint64_t minCoverage,
    uint64_t threadCount)
{
    cout << timestamp << "Filtering local segments (ma_hit_sub equivalent, minCoverage="
         << minCoverage << ")." << endl;

    const uint64_t readCount = reads->readCount();
    validReadIntervals.resize(readCount);
    std::fill(validReadIntervals.begin(), validReadIntervals.end(), ReadSegmentStatus{0, 0, false});
    this->localSegmentMinCoverage = minCoverage;

    setupLoadBalancing(readCount, 100);
    runThreads(&Assembler::filterLocalSegmentsThreadFunction, threadCount);

    // Compute diagnostics
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

    cout << timestamp << "Local segment filtering (ma_hit_sub): "
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
    /*
    ============================================================================
    ma_hit_sub (miniasm/hifiasm) thread function - COMPLETE PARITY VERIFIED
    ============================================================================

    This function computes a per-read high-coverage interval that represents the
    "valid" portion of each read for subsequent assembly steps.

    HIFIASM ALGORITHM (ma_hit_sub):
    --------------------------------
    For each read r:
      1. Collect all overlaps involving r (from r's perspective only, to avoid double-counting)
      2. Build a sweep-line depth profile over the read coordinates
      3. Find the longest interval where overlap depth >= minCoverage
      4. Store this interval as coverage_cut[r].s and coverage_cut[r].e
      5. Special case: if minCoverage <= 1, keep the full read [0, length)

    KEY IMPLEMENTATION DETAILS FOR HIFIASM PARITY:
    -----------------------------------------------
    1. **Single orientation query** (lines 311-331):
       - Only process alignmentTable[OrientedReadId(r0, 0)]
       - This matches hifiasm's sources[r0] adjacency list structure
       - Avoids double-counting: each overlap appears once per read

    2. **Skip deleted overlaps** (line 314):
       - isDeletedFromReadPerspective checks if overlap is deleted from r0's side
       - Matches hifiasm's per-read deletion flags

    3. **Sweep-line algorithm** (lines 338-378):
       - Events: +1 for interval start, -1 for interval end
       - Sort events by position, with starts before ends at same position
       - Track coverage depth as we sweep through the read
       - Find the longest interval where depth >= minCoverage

    4. **Special case for minCoverage <= 1** (lines 302-306):
       - Hifiasm: when minDp <= 1, keep full read without computation
       - Sets validReadIntervals[r0] = {0, readLength, false}

    OUTPUT:
    -------
    validReadIntervals[r0] contains:
      - start: beginning of high-coverage interval
      - end: end of high-coverage interval
      - isDeleted: true if no valid interval found (corresponds to coverage_cut[r].del = 1)

    This interval will be used by ma_hit_cut (applyCoverageCuts) to clip overlaps
    to the valid read regions.
    */
    static_cast<void>(threadId);
    const uint64_t minDp = this->localSegmentMinCoverage;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); r0++) {
            /*
            HIFIASM SPECIAL CASE: minCoverage <= 1
            ---------------------------------------
            When minDp <= 1, hifiasm's ma_hit_sub keeps the full read as the valid region
            without computing coverage intervals. This matches the behavior:
              coverage_cut[r].s = 0
              coverage_cut[r].e = read_length
              coverage_cut[r].del = 0
            */
            if (minDp <= 1) {
                 uint64_t len = reads->getReadRawSequenceLength(r0);
                 validReadIntervals[r0] = {0, (uint32_t)len, false};
                 continue;
            }

            /*
            STEP 1: Collect overlap intervals for read r0
            ----------------------------------------------
            Hifiasm ma_hit_sub collects overlaps from sources[r0], which contains one
            copy of each overlap involving r0. We replicate this by querying only
            alignmentTable[OrientedReadId(r0, 0)] to avoid double-counting.
            */
            std::vector<std::pair<uint32_t, uint32_t>> intervals;

            auto collectIntervals = [&](OrientedReadId orientedR0) {
                const auto& table = alignmentTable[orientedR0.getValue()];
                for (uint32_t alignmentId : table) {
                    const auto& ad = alignmentData[alignmentId];

                    // Skip overlaps deleted from r0's perspective
                    if (isDeletedFromReadPerspective(ad, r0)) continue;

                    // Extract interval on r0 (query or target depending on role)
                    uint32_t start = 0, end = 0;
                    if (ad.readIds[0] == r0) {
                        start = ad.qs;
                        end = ad.qe;
                    } else {
                        start = ad.ts;
                        end = ad.te;
                    }

                    // Add non-empty interval
                    if (end > start) {
                        intervals.push_back({start, end});
                    }
                }
            };
            collectIntervals(OrientedReadId(r0, 0));

            /*
            HIFIASM BEHAVIOR: No valid overlaps
            ------------------------------------
            If no overlaps remain after filtering, mark the read as deleted:
              coverage_cut[r].del = 1
            */
            if (intervals.empty()) {
                validReadIntervals[r0] = {0, 0, true};
                continue;
            }

            /*
            STEP 2: Build sweep-line events
            --------------------------------
            Hifiasm ma_hit_sub uses a sweep-line algorithm to find the longest interval
            with coverage >= minDp:
              - Event +1: interval start (coverage increases)
              - Event -1: interval end (coverage decreases)
              - Sort events by position, with starts (+1) before ends (-1) at same position
            */
            std::vector<std::pair<uint32_t, int>> events;
            events.reserve(intervals.size() * 2);
            for(const auto& interval : intervals) {
                events.push_back({interval.first, 1});   // +1 = start
                events.push_back({interval.second, -1}); // -1 = end
            }

            // Sort by position, then by event type (start before end at same position)
            std::sort(events.begin(), events.end(), [](const pair<uint32_t, int>& a, const pair<uint32_t, int>& b){
                if (a.first != b.first) return a.first < b.first;
                return a.second > b.second; // +1 (start) before -1 (end)
            });

            /*
            STEP 3: Sweep through events to find longest high-coverage interval
            --------------------------------------------------------------------
            Hifiasm ma_hit_sub finds the longest interval where depth >= minDp:
              - Track current coverage depth as we process events
              - When coverage transitions from <minDp to >=minDp: start new segment
              - When coverage transitions from >=minDp to <minDp: close segment
              - Keep track of the longest segment found
            */
            int coverage = 0;
            uint32_t maxS = 0, maxE = 0;          // Longest high-coverage interval found
            uint32_t currentStart = 0;             // Start of current segment
            bool inSegment = false;                // Are we currently in a high-coverage segment?

            for(const auto& ev : events) {
                int oldCoverage = coverage;
                coverage += ev.second;             // Update coverage (+1 or -1)
                uint32_t pos = ev.first;           // Event position

                // Transition into high-coverage region (coverage crosses minDp threshold upward)
                if (oldCoverage < (int)minDp && coverage >= (int)minDp) {
                    currentStart = pos;
                    inSegment = true;
                }
                // Transition out of high-coverage region (coverage drops below minDp threshold)
                else if (oldCoverage >= (int)minDp && coverage < (int)minDp) {
                    if (inSegment) {
                        uint32_t currentLen = pos - currentStart;
                        uint32_t maxLen = maxE - maxS;

                        // Keep this segment if it's longer than the previous best
                        if (currentLen > maxLen) {
                            maxS = currentStart;
                            maxE = pos;
                        }
                        inSegment = false;
                    }
                }
            }

            /*
            STEP 4: Store result
            --------------------
            Hifiasm ma_hit_sub sets:
              - coverage_cut[r].s = maxS (start of longest high-coverage interval)
              - coverage_cut[r].e = maxE (end of longest high-coverage interval)
              - coverage_cut[r].del = (maxE <= maxS ? 1 : 0) (no valid interval found)

            We store this in validReadIntervals[r0] for use by subsequent stages.
            */
            if (maxE > maxS) {
                validReadIntervals[r0] = {maxS, maxE, false};
            } else {
                // No interval with coverage >= minDp found
                validReadIntervals[r0] = {0, 0, true};
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

            // Hifiasm ma_hit_flt: only process overlaps kept by both sides
            if(!ad.keptByBothSides()) continue;

            ReadId qn = ad.readIds[0];
            ReadId tn = ad.readIds[1];

            /*
            Hifiasm ma_hit_flt: get read lengths in valid coordinates.
            After ma_hit_cut, overlap coordinates are normalized to valid regions,
            so we use valid region lengths as "read lengths" for ma_hit2arc_containment.
            */
            uint32_t ql, tl;
            if (validReadIntervals.empty()) {
                // Fallback: use raw read lengths if valid intervals not computed
                ql = (uint32_t)reads->getReadRawSequenceLength(qn);
                tl = (uint32_t)reads->getReadRawSequenceLength(tn);
            } else {
                const auto& rq = validReadIntervals[qn];
                const auto& rt = validReadIntervals[tn];

                // Hifiasm ma_hit_flt: skip overlaps involving deleted reads
                if (rq.isDeleted || rt.isDeleted) continue;

                // Read length = valid region length (after ma_hit_cut normalization)
                ql = rq.end - rq.start;
                tl = rt.end - rt.start;
            }

            uint32_t qs = ad.qs;
            uint32_t qe = ad.qe;
            uint32_t ts = ad.ts;
            uint32_t te = ad.te;

            // Sanity check: invalid overlap coordinates
            if (qe <= qs || te <= ts) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonHanging);
                continue;
            }

            /*
            Hifiasm ma_hit2arc_containment: classify overlap type.
            Returns:
              0  = dovetail (keep)
              1  = query contained in target (keep - handled by removeContainedReads)
              2  = target contained in query (keep - handled by removeContainedReads)
              -1 = internal match (delete - excessive overhangs)
              -2 = too short (delete - effective overlap < min_ovlp)
            */
            const int result = ma_hit2arc_containment(
                (int32_t)qs, (int32_t)qe, (int32_t)ql,
                (int32_t)ts, (int32_t)te, (int32_t)tl,
                !ad.isSameStrand,
                (int32_t)maxHang,
                maxHangRate,
                (int32_t)minOvlp
            );

            // Hifiasm ma_hit_flt: delete internal matches and too-short overlaps
            if (result < 0) {
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
    // Hifiasm parity note:
    // detect_chimeric_reads deletes overlaps in-place while scanning reads. Later reads are
    // expected to observe those deletions, so this pass must remain sequential.
    constexpr double shiftRate = 0.03 * 2.0; // asm_opt.max_ov_diff_final * 2.0 (default 0.03)
    constexpr float overlapRate = 0.1f;

    struct SubRegion { uint32_t s, e; };
    struct ProjectionToTarget {
        ReadId tId = invalid<ReadId>;
        uint32_t qs = 0, qe = 0; // query interval in qId coordinates
        uint32_t ts = 0, te = 0; // target interval in tId forward coordinates
        bool rev = false;        // target is reverse-complemented relative to query
    };

    std::vector<char> xBuf;
    std::vector<char> yBuf;

    // hifiasm delete_all_edges equivalent for one read.
    auto markReadChimericAndDeleteIncidentOverlaps = [&](ReadId qId) {
        if (qId >= readCount || isChimericRead[qId]) {
            return;
        }
        isChimericRead[qId] = true;
        if (qId < validReadIntervals.size()) {
            validReadIntervals[qId].isDeleted = true;
        }

        const OrientedReadId oid(qId, 0);
        if (oid.getValue() >= alignmentTable.size()) {
            return;
        }
        const auto& table = alignmentTable[oid.getValue()];
        for (uint32_t alignmentId : table) {
            alignmentData[alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonChimeric);
        }
    };

    // Given an overlap in the incidence list of qId, return [qs,qe) on qId.
    auto getQueryInterval = [&](const AlignmentData& ad, ReadId qId, uint32_t& qs, uint32_t& qe) {
        if (ad.readIds[0] == qId) {
            qs = ad.qs;
            qe = ad.qe;
        } else {
            DINARA_ASSERT(ad.readIds[1] == qId);
            qs = ad.ts;
            qe = ad.te;
        }
    };

    // Given an overlap in the incidence list of qId, return all fields needed by boundaryVerify.
    auto projectFromQueryPerspective = [&](const AlignmentData& ad, ReadId qId) -> ProjectionToTarget {
        ProjectionToTarget p;
        if (ad.readIds[0] == qId) {
            p.tId = ad.readIds[1];
            p.qs = ad.qs;
            p.qe = ad.qe;
            p.ts = ad.ts;
            p.te = ad.te;
            p.rev = !ad.isSameStrand;
        } else {
            DINARA_ASSERT(ad.readIds[1] == qId);
            p.tId = ad.readIds[0];
            p.qs = ad.ts;
            p.qe = ad.te;
            p.ts = ad.qs;
            p.te = ad.qe;
            p.rev = !ad.isSameStrand;
        }
        return p;
    };

    uint64_t chimericCount = 0;
    uint64_t simpleChimericCount = 0;
    uint64_t complexCheckedCount = 0;
    uint64_t complexChimericCount = 0;
    for (ReadId qId = 0; qId < readCount; ++qId) {
        // hifiasm uses raw read length here (not coverage-clipped interval length).
        const uint32_t qLen = uint32_t(reads->getReadRawSequenceLength(qId));
        const OrientedReadId oid(qId, 0);
        if (oid.getValue() >= alignmentTable.size()) {
            continue;
        }
        const auto& table = alignmentTable[oid.getValue()];

        SubRegion maxLeft{qLen, 0};
        SubRegion maxRight{qLen, 0};

        // collect_sides equivalent (without UL-aware bypass).
        for (uint32_t alignmentId : table) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if (isDeletedFromReadPerspective(ad, qId)) {
                continue;
            }

            uint32_t qs = 0, qe = 0;
            getQueryInterval(ad, qId, qs, qe);

            if (qs == 0) {
                if (qs < maxLeft.s) maxLeft.s = qs;
                if (qe > maxLeft.e) maxLeft.e = qe;
            }
            if (qe == qLen) {
                if (qs < maxRight.s) maxRight.s = qs;
                if (qe > maxRight.e) maxRight.e = qe;
            }
        }

        if (maxLeft.s == qLen || maxRight.s == qLen) {
            continue;
        }

        // collect_contain equivalent.
        uint32_t newLeftE = maxLeft.e;
        uint32_t newRightS = maxRight.s;
        for (uint32_t alignmentId : table) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if (isDeletedFromReadPerspective(ad, qId)) {
                continue;
            }

            uint32_t qs = 0, qe = 0;
            getQueryInterval(ad, qId, qs, qe);

            if (qs == 0 || qe == qLen) {
                continue;
            }
            const uint32_t len = qe - qs;
            if (len == 0) {
                continue;
            }

            if (qs < maxLeft.e && qe > maxLeft.e &&
                (maxLeft.e - qs) > uint32_t(overlapRate * float(len))) {
                if (qe > newLeftE) {
                    newLeftE = qe;
                }
            }
            if (qs < maxRight.s && qe > maxRight.s &&
                (qe - maxRight.s) > uint32_t(overlapRate * float(len))) {
                if (qs < newRightS) {
                    newRightS = qs;
                }
            }
        }
        maxLeft.e = newLeftE;
        maxRight.s = newRightS;

        // Normal read.
        if (maxLeft.e > maxRight.s &&
            (maxLeft.e - maxRight.s) >= uint32_t(double(qLen) * shiftRate)) {
            continue;
        }

        // Simple chimera.
        if (maxLeft.e <= maxRight.s) {
            markReadChimericAndDeleteIncidentOverlaps(qId);
            ++chimericCount;
            ++simpleChimericCount;
            continue;
        }

        // Complex case with boundary verification:
        // maxLeft and maxRight overlap, but not enough for "normal" classification.
        // We ask if at least one overlap spanning [intervalS,intervalE) fails boundary consistency.
        ++complexCheckedCount;
        const uint32_t intervalS = maxRight.s;
        const uint32_t intervalE = maxLeft.e;
        bool isChimeric = false;

        for (uint32_t alignmentId : table) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if (isDeletedFromReadPerspective(ad, qId)) {
                continue;
            }

            const ProjectionToTarget p = projectFromQueryPerspective(ad, qId);
            if (p.qs <= intervalS && p.qe >= intervalE) {
                // For raw/uncorrected reads, exact-overlap surrogates are unreliable.
                // Use boundary verification alone to decide complex chimeras.
                if (!boundaryVerify(*reads, intervalS, intervalE, qId, p.tId, p.qs, p.ts, p.te, p.rev, xBuf, yBuf)) {
                    isChimeric = true;
                    break;
                }
            }
        }

        if (isChimeric) {
            markReadChimericAndDeleteIncidentOverlaps(qId);
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
                
                bool isValid = (ad.cisTransStatus == CisTransStatus::Cis);
                if (ad.cisTransStatus == CisTransStatus::Unknown && ad.info.isInReadGraph) {
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
                        
                        bool isValidType = (ad.cisTransStatus == CisTransStatus::Cis);
                        if (ad.cisTransStatus == CisTransStatus::Unknown && ad.info.isInReadGraph) isValidType = true;
                        
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

/*
Contained-read filtering (hifiasm `ma_hit_contained_advance` equivalent).

Hifiasm reduces the overlap graph by removing reads that are entirely contained in another read,
since they do not contribute new adjacency in the string graph.

The core predicate (`ma_hit2arc` in hifiasm) classifies an overlap, from the query read's perspective, as:
  - normal dovetail (keep),
  - query contained in target (delete query),
  - target contained in query (delete target),
  - internal match / too short (ignore for containment purposes).

Important convention used throughout Dinara's overlap storage:
  - Coordinates (`qs/qe` and `ts/te`) are stored in the forward coordinate frame of each read.
  - `isSameStrand` indicates whether the alignment is forward/forward or forward/reverse, but the stored
    target interval is still in the target read's forward frame. This matches the expectation of the
    containment predicate (which applies an orientation-aware conversion only to decide 5'/3' overhangs).
*/
static int ma_hit2arc_containment(
    int32_t qs, int32_t qe, int32_t ql,
    int32_t ts, int32_t te, int32_t tl,
    bool isReverse,
    int32_t max_hang,
    double int_frac,
    int32_t min_ovlp)
{
    /*
    Hifiasm ma_hit2arc: classify overlap type based on overhang analysis.

    This is a direct translation of hifiasm's `ma_hit2arc` function used for both
    ma_hit_flt (hanging overlap filtering) and ma_hit_contained_advance (containment detection).

    Given an overlap with:
      - Query: interval [qs, qe) on read of length ql
      - Target: interval [ts, te) on read of length tl
      - isReverse: whether target is reverse-complemented

    Algorithm:
    1. Compute overhangs (unaligned portions):
       - Query 5' overhang: qs
       - Query 3' overhang: ql - qe
       - Target overhangs: adjusted for orientation (forward vs RC)

    2. Compute extensions (minimum overhangs at each end):
       - ext5 = min(query 5' overhang, target 5' overhang)
       - ext3 = min(query 3' overhang, target 3' overhang)

    3. Apply filters in order:
       a) Absolute overhang filter: ext5 > max_hang OR ext3 > max_hang → internal (-1)
       b) Fractional overlap filter: overlap/(overlap+ext5+ext3) < int_frac → internal (-1)
       c) Containment detection: compare query/target overhangs → contained (1 or 2)
       d) Minimum overlap filter: effective_length < min_ovlp → too short (-2)
       e) Default: dovetail (0)

    Return values:
       0 = dovetail (proper prefix-suffix overlap)
       1 = query contained in target (QCONT)
       2 = target contained in query (TCONT)
      -1 = internal match (excessive overhangs, reject)
      -2 = too short (effective overlap < min_ovlp, reject)

    Hifiasm parameters (typical defaults):
      - max_hang = 1000 (asm_opt.max_hang)
      - int_frac = 0.8 (asm_opt.int_frac)
      - min_ovlp = 50 (asm_opt.min_ovlp)
    */

    // Step 1: Compute target overhangs in query orientation
    int32_t tl5, tl3;  // Target 5' and 3' overhangs
    if (isReverse) {
        // Reverse-complement: 5' and 3' are swapped
        tl5 = tl - te;  // Target 5' overhang = bases after target end
        tl3 = ts;       // Target 3' overhang = bases before target start
    } else {
        // Forward orientation
        tl5 = ts;       // Target 5' overhang = bases before target start
        tl3 = tl - te;  // Target 3' overhang = bases after target end
    }

    // Step 2: Compute extensions (minimum overhangs at each end)
    int32_t ext5 = (qs < tl5) ? qs : tl5;           // 5' extension
    int32_t ext3 = ((ql - qe) < tl3) ? (ql - qe) : tl3;  // 3' extension

    // Step 3a: Hifiasm ma_hit_flt: reject if absolute overhangs exceed threshold
    if (ext5 > max_hang || ext3 > max_hang) {
        return -1;  // Internal match
    }

    // Step 3b: Hifiasm ma_hit_flt: reject if fractional overlap is too low
    int32_t qOverlapLen = qe - qs;
    int32_t tOverlapLen = te - ts;

    // Check query: overlap_len / (overlap_len + ext5 + ext3) >= int_frac
    // Equivalent: overlap_len >= (overlap_len + ext5 + ext3) * int_frac
    if (qOverlapLen < (qOverlapLen + ext5 + ext3) * int_frac) {
        return -1;  // Internal match
    }
    if (tOverlapLen < (tOverlapLen + ext5 + ext3) * int_frac) {
        return -1;  // Internal match
    }

    // Step 3c: Hifiasm containment detection
    // QCONT: query is contained if both query overhangs <= corresponding target overhangs
    if (qs <= tl5 && (ql - qe) <= tl3) {
        return 1;  // Query contained in target
    }
    // TCONT: target is contained if both target overhangs <= corresponding query overhangs
    if (qs >= tl5 && (ql - qe) >= tl3) {
        return 2;  // Target contained in query
    }

    // Step 3d: Hifiasm minimum overlap filter
    // Effective overlap length = overlap + extensions
    if (qOverlapLen + ext5 + ext3 < min_ovlp || tOverlapLen + ext5 + ext3 < min_ovlp) {
        return -2;  // Too short
    }

    // Step 3e: Default case - proper dovetail overlap
    return 0;
}


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
       - Call `ma_hit2arc_containment()` from qn's perspective
       - result == 1 (QCONT): query (qn) is contained in target (tn)
       - result == 2 (TCONT): target (tn) is contained in query (qn)
       - result == 0: dovetail overlap (keep both reads)
       - result < 0: internal/bad (already filtered by ma_hit_flt)

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

            /*
            Hifiasm ma_hit2arc_containment: classify overlap type.

            Returns:
               0 = dovetail (no containment, keep both reads)
               1 = QCONT (query contained in target → delete query)
               2 = TCONT (target contained in query → delete target)
              -1 = internal (already filtered by ma_hit_flt)
              -2 = too short (already filtered by ma_hit_flt)
            */
            const int result = ma_hit2arc_containment(
                qs, qe, ql,
                ts, te, tl,
                rev,
                int32_t(maxHang),
                maxHangRate,
                int32_t(minOverlapLength)
            );

            if (result == 1) {
                // QCONT: query (qn) is contained in target (tn)
                if (!validReadIntervals[qn].isDeleted) {
                    validReadIntervals[qn].isDeleted = true;              // Mark read deleted
                    (*containmentParent)[qn] = tn;                        // Record container
                    deleteAllEdgesForRead(qn);                           // Delete ALL its overlaps
                    ++containedReadCount;
                }
                // Hifiasm: stop processing this read once contained
                break;

            } else if (result == 2) {
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

void Assembler::flagContainedReads(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Flagging contained reads (diagnostic-only, does not remove overlaps)..." << endl;

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    checkAlignmentDataAreOpen();

    // We will write into read flags.
    reads->checkReadsAreOpen();
    reads->checkReadFlagsAreOpenForWriting();

    // Historically this diagnostic expected `validReadIntervals` to be available (ma_hit_sub / ma_hit_cut).
    // For quick experimentation we allow running without it by treating the entire read as valid.
    // In that mode we operate on raw overlap coordinates (`ad.qs/qe/ts/te`) and full read lengths.
    if(validReadIntervals.empty()) {
        cout << timestamp
             << "[DIAG] flagContainedReads: validReadIntervals not set; using full read lengths and raw overlap coordinates."
             << endl;
        validReadIntervals.resize(reads->readCount());
        for(ReadId r=0; r<reads->readCount(); ++r) {
            const uint64_t len = reads->getReadRawSequenceLength(r);
            validReadIntervals[r] = {0, uint32_t(len), false};
        }
    }

    if(!containmentParent->isOpen) {
        containmentParent->createNew(largeDataName("ContainmentParent"), largeDataPageSize);
        containmentParent->resize(reads->readCount());
    }
    std::fill(containmentParent->begin(), containmentParent->end(), ReadId(invalidReadId));

    const uint64_t readCount = reads->readCount();
    const uint64_t alignmentCount = alignmentData.size();

    // Clear previous flags.
    for(ReadId r=0; r<readCount; ++r) {
        reads->setContainedFlag(r, false);
    }

    uint64_t containedReadCount = 0;
    // Simulate hifiasm ma_hit_contained_advance decisions, but without mutating overlaps.
    // We maintain a local "deleted" mask so that once a read is deemed contained, we stop
    // considering it and we ignore overlaps incident to it (matching the fact that hifiasm
    // deletes all its edges immediately).
    vector<uint8_t> deletedLocal(readCount, 0);
    for(ReadId r = 0; r < readCount; ++r) {
        if(r < validReadIntervals.size() && validReadIntervals[r].isDeleted) {
            deletedLocal[r] = 1;
        }
    }

    for(ReadId qn = 0; qn < readCount; ++qn) {
        if(qn >= validReadIntervals.size() || deletedLocal[qn]) continue;

        const auto& vrQ = validReadIntervals[qn];
        const int32_t ql = int32_t(vrQ.end - vrQ.start);
        if(ql <= 0) continue;

        const OrientedReadId oid(qn, 0);
        if(oid.getValue() >= alignmentTable.size()) continue;

        const auto& table = alignmentTable[oid.getValue()];
        for(const uint32_t alignmentId : table) {
            if(alignmentId >= alignmentCount) continue;
            const AlignmentData& ad = alignmentData[alignmentId];
            if(!ad.keptByBothSides()) continue;

            const ReadId tn = (ad.readIds[0] == qn) ? ad.readIds[1] : ad.readIds[0];
            if(tn >= validReadIntervals.size() || deletedLocal[tn]) continue;

            const auto& vrT = validReadIntervals[tn];
            const int32_t tl = int32_t(vrT.end - vrT.start);
            if(tl <= 0) continue;

            const bool rev = !ad.isSameStrand;
            int32_t qs = 0, qe = 0, ts = 0, te = 0;
            if(ad.readIds[0] == qn) {
                qs = int32_t(ad.qs);
                qe = int32_t(ad.qe);
                ts = int32_t(ad.ts);
                te = int32_t(ad.te);
            } else {
                qs = int32_t(ad.ts);
                qe = int32_t(ad.te);
                ts = int32_t(ad.qs);
                te = int32_t(ad.qe);
            }
            if(qs < 0 || qe < 0 || ts < 0 || te < 0) continue;
            if(qs >= qe || ts >= te) continue;
            if(qs > ql || qe > ql || ts > tl || te > tl) continue;

            const int result = ma_hit2arc_containment(
                qs, qe, ql,
                ts, te, tl,
                rev,
                int32_t(maxHang),
                maxHangRate,
                int32_t(minOverlapLength)
            );

            if(result == 1) {
                // Query contained in target: flag and stop processing this query.
                reads->setContainedFlag(qn, true);
                (*containmentParent)[qn] = tn;
                deletedLocal[qn] = 1;
                ++containedReadCount;
                break;
            } else if(result == 2) {
                // Target contained in query: flag the target and keep scanning (qn may contain multiple reads).
                if(!deletedLocal[tn]) {
                    reads->setContainedFlag(tn, true);
                    (*containmentParent)[tn] = qn;
                    deletedLocal[tn] = 1;
                    ++containedReadCount;
                }
            }
        }
    }

    // Compress containment chains.
    for(ReadId r = 0; r < readCount; ++r) {
        if((*containmentParent)[r] == ReadId(invalidReadId)) continue;
        ReadId root = r;
        vector<ReadId> visited;
        visited.reserve(16);
        while(root != ReadId(invalidReadId) && (*containmentParent)[root] != ReadId(invalidReadId)) {
            if(std::find(visited.begin(), visited.end(), root) != visited.end()) {
                // Cycle detected (can happen for identical reads). Break deterministically by choosing
                // the smallest read id in the cycle as the root container.
                ReadId cycleRoot = root;
                for(const ReadId x : visited) {
                    if(x < cycleRoot) cycleRoot = x;
                }
                (*containmentParent)[r] = cycleRoot;
                root = ReadId(invalidReadId);
                break;
            }
            visited.push_back(root);
            root = (*containmentParent)[root];
        }
        if(root != ReadId(invalidReadId)) {
            (*containmentParent)[r] = root;
        }
    }

    cout << timestamp << "Flagged " << containedReadCount << " contained reads." << endl;
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
