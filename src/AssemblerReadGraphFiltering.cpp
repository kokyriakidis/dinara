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
#include "Reads.hpp"
#include "hifiasmCoordinateTransforms.hpp"
#include "timestamp.hpp"
#include <algorithm>
#include <vector>
#include <limits>

#include <iostream>

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
    constexpr uint32_t HifiasmWindow = 375;
    constexpr uint32_t HifiasmThresholdMaxSize = 31;
    constexpr double HifiasmMaxOvDiffEc = 0.04;

    inline char complementBase(char c)
    {
        switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        default: return 'N';
        }
    }

    inline void extractReadSubstring(
        const Reads& reads,
        ReadId readId,
        uint32_t start,
        uint32_t length,
        bool reverseComplement,
        std::vector<char>& out)
    {
        out.resize(length);
        const auto& read = reads.getRead(readId);
        const uint32_t readLen = uint32_t(read.baseCount);
        if (length == 0) return;

        if (!reverseComplement) {
            for (uint32_t i = 0; i < length; ++i) {
                const uint32_t pos = start + i;
                out[i] = (pos < readLen) ? read[pos].character() : 'N';
            }
        } else {
            for (uint32_t i = 0; i < length; ++i) {
                const uint32_t pos = start + (length - 1U - i);
                const char b = (pos < readLen) ? read[pos].character() : 'N';
                out[i] = complementBase(b);
            }
        }
    }

    inline bool determineOverlapRegion(
        int threshold,
        int64_t yStart,
        int64_t yLen,
        int64_t windowLen,
        int& extraBegin,
        int& extraEnd,
        int64_t& clippedYStart,
        int64_t& clippedYLen)
    {
        if (yStart < 0 || yStart >= yLen ||
            (yLen - yStart + 2 * threshold + HifiasmThresholdMaxSize) < windowLen) {
            return false;
        }

        extraBegin = 0;
        extraEnd = 0;
        clippedYStart = yStart - threshold;
        clippedYLen = std::min<int64_t>(windowLen, yLen - clippedYStart);
        extraEnd = int(windowLen - clippedYLen);
        if (clippedYStart < 0) {
            extraBegin = int(-clippedYStart);
            clippedYStart = 0;
            clippedYLen -= extraBegin;
        }
        return clippedYLen > 0;
    }

    inline int minEditDistancePatternToAnySubstring(const char* pattern, int m, const char* text, int n)
    {
        std::vector<int16_t> prev(n + 1);
        std::vector<int16_t> curr(n + 1);
        std::fill(prev.begin(), prev.end(), 0);

        for (int i = 1; i <= m; ++i) {
            curr[0] = int16_t(i);
            for (int j = 1; j <= n; ++j) {
                const int cost = (pattern[i - 1] == text[j - 1]) ? 0 : 1;
                const int a = prev[j] + 1;
                const int b = curr[j - 1] + 1;
                const int c = prev[j - 1] + cost;
                curr[j] = int16_t(std::min({a, b, c}));
            }
            prev.swap(curr);
        }

        int best = prev[0];
        for (int j = 1; j <= n; ++j) best = std::min(best, int(prev[j]));
        return best;
    }

    inline bool verifySingleWindow(
        const Reads& reads,
        uint32_t xStart,
        uint32_t xEndInclusive,
        uint32_t overlapXs,
        uint32_t overlapYs,
        ReadId xId,
        ReadId yId,
        bool xStrand,
        double maxOvDiffEc,
        std::vector<char>& xBuf,
        std::vector<char>& yBuf)
    {
        if (xEndInclusive < xStart) return false;
        const uint32_t xLen = xEndInclusive - xStart + 1;
        if (xLen == 0) return false;

        int threshold = int(double(xLen) * maxOvDiffEc);
        if (threshold == 0 && xLen >= 4) threshold = 1;
        const int64_t windowLen = int64_t(xLen) + (int64_t(threshold) << 1);

        const int64_t yReadLen = int64_t(reads.getReadRawSequenceLength(yId));
        const int64_t yStart = int64_t(xStart) - int64_t(overlapXs) + int64_t(overlapYs);

        int extraBegin = 0;
        int extraEnd = 0;
        int64_t clippedYStart = 0;
        int64_t clippedYLen = 0;
        if (!determineOverlapRegion(threshold, yStart, yReadLen, windowLen, extraBegin, extraEnd, clippedYStart, clippedYLen)) {
            return false;
        }

        extractReadSubstring(reads, xId, xStart, xLen, xStrand, xBuf);

        yBuf.assign(size_t(windowLen), 'N');
        std::vector<char> yCore;
        extractReadSubstring(reads, yId, uint32_t(clippedYStart), uint32_t(clippedYLen), false, yCore);
        std::copy(yCore.begin(), yCore.end(), yBuf.begin() + extraBegin);

        const int dist = minEditDistancePatternToAnySubstring(xBuf.data(), int(xBuf.size()), yBuf.data(), int(yBuf.size()));
        return dist <= threshold;
    }

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
        if (qIntervalEnd <= qIntervalStart) return false;
        const uint32_t intervalLen = qIntervalEnd - qIntervalStart;

        const uint32_t tLen = uint32_t(reads.getReadRawSequenceLength(tId));
        const uint32_t xs = qs;
        const uint32_t ys = rev ? (tLen - te) : ts;

        uint32_t tIntervalStart = (qIntervalStart - xs) + ys;
        if (tIntervalStart >= tLen) return false;
        uint32_t tIntervalEndInclusive = tIntervalStart + intervalLen - 1;
        if (tIntervalEndInclusive >= tLen) tIntervalEndInclusive = tLen - 1;
        if (tIntervalEndInclusive < tIntervalStart) return false;

        const uint32_t tIntervalLen = tIntervalEndInclusive - tIntervalStart + 1;
        if (tIntervalLen <= HifiasmWindow) {
            return verifySingleWindow(
                reads,
                tIntervalStart,
                tIntervalEndInclusive,
                ys,
                xs,
                tId,
                qId,
                rev,
                HifiasmMaxOvDiffEc,
                xBuf,
                yBuf);
        }

        if (!verifySingleWindow(
                reads,
                tIntervalStart,
                tIntervalStart + HifiasmWindow - 1,
                ys,
                xs,
                tId,
                qId,
                rev,
                HifiasmMaxOvDiffEc,
                xBuf,
                yBuf)) {
            return false;
        }
        if (!verifySingleWindow(
                reads,
                tIntervalEndInclusive - HifiasmWindow + 1,
                tIntervalEndInclusive,
                ys,
                xs,
                tId,
                qId,
                rev,
                HifiasmMaxOvDiffEc,
                xBuf,
                yBuf)) {
            return false;
        }
        return true;
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
    cout << timestamp << "Filtering local segments (min coverage=" << minCoverage << ")..." << endl;

    validReadIntervals.resize(reads->readCount());
    std::fill(validReadIntervals.begin(), validReadIntervals.end(), ReadSegmentStatus{0, 0, false});
    this->localSegmentMinCoverage = minCoverage;

    setupLoadBalancing(reads->readCount(), 100);
    runThreads(&Assembler::filterLocalSegmentsThreadFunction, threadCount);
    
    cout << timestamp << "Local segment filtering complete." << endl;
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
    Thread function for filterLocalSegments (ma_hit_sub).

    For each read r0, we look only at overlaps in `alignmentTable[OrientedReadId(r0,0)]` and
    keep those that are not deleted from r0's perspective. This matches hifiasm/miniasm's
    per-read adjacency (`sources[r0]`) and avoids double-counting the same overlap.
    */
    static_cast<void>(threadId);
    const uint64_t minDp = this->localSegmentMinCoverage;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); r0++) {
            if (minDp <= 1) {
                 uint64_t len = reads->getReadRawSequenceLength(r0);
                 validReadIntervals[r0] = {0, (uint32_t)len, false};
                 continue;
            }

            std::vector<std::pair<uint32_t, uint32_t>> intervals;
            
            auto collectIntervals = [&](OrientedReadId orientedR0) {
                const auto& table = alignmentTable[orientedR0.getValue()];
                for (uint32_t alignmentId : table) {
                    const auto& ad = alignmentData[alignmentId];
                    if (isDeletedFromReadPerspective(ad, r0)) continue;

                    uint32_t start = 0, end = 0;
                    
                    if (ad.readIds[0] == r0) {
                        start = ad.qs;
                        end = ad.qe;
                    } else {
                        start = ad.ts;
                        end = ad.te;
                    }
                    
                    if (end > start) {
                        intervals.push_back({start, end});
                    }
                }
            };
            collectIntervals(OrientedReadId(r0, 0));
            
            if (intervals.empty()) {
                validReadIntervals[r0] = {0, 0, true};
                continue; 
            }
            
            std::vector<std::pair<uint32_t, int>> events;
            events.reserve(intervals.size() * 2);
            for(const auto& interval : intervals) {
                events.push_back({interval.first, 1});
                events.push_back({interval.second, -1});
            }
            
            std::sort(events.begin(), events.end(), [](const pair<uint32_t, int>& a, const pair<uint32_t, int>& b){
                if (a.first != b.first) return a.first < b.first;
                return a.second > b.second;
            });
            
            int coverage = 0;
            uint32_t maxS = 0, maxE = 0;
            
            uint32_t currentStart = 0;
            bool inSegment = false;
            
            for(const auto& ev : events) {
                int oldCoverage = coverage;
                coverage += ev.second;
                
                uint32_t pos = ev.first;
                
                if (oldCoverage < (int)minDp && coverage >= (int)minDp) {
                    currentStart = pos;
                    inSegment = true;
                }
                else if (oldCoverage >= (int)minDp && coverage < (int)minDp) {
                    if (inSegment) {
                        uint32_t currentLen = pos - currentStart;
                        uint32_t maxLen = maxE - maxS;
                        
                        if (currentLen > maxLen) {
                            maxS = currentStart;
                            maxE = pos;
                        }
                        inSegment = false;
                    }
                }
            }
            
            if (maxE > maxS) {
                validReadIntervals[r0] = {maxS, maxE, false};
            } else {
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
    cout << timestamp << "Applying coverage cuts (ma_hit_cut equivalent)." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    this->coverageCutMinOverlap = minOverlapLength;

    setupLoadBalancing(alignmentData.size(), 10000); 
    runThreads(&Assembler::applyCoverageCutsToAlignmentsThreadFunction, threadCount);

    setupLoadBalancing(reads->readCount(), 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);
    
    cout << timestamp << "Coverage cuts applied." << endl;
}

void Assembler::applyCoverageCutsToAlignmentsThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    uint64_t begin, end;
    const uint64_t minLen = this->coverageCutMinOverlap;

    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; i++) {
            AlignmentData& ad = alignmentData[i];
            
            if(!ad.keptByBothSides()) continue;

            ReadId qn = ad.readIds[0];
            ReadId tn = ad.readIds[1];

            if (qn >= validReadIntervals.size() || tn >= validReadIntervals.size()) continue;

            const auto& rq = validReadIntervals[qn];
            const auto& rt = validReadIntervals[tn];

            if (rq.isDeleted || rt.isDeleted) continue;

            int32_t qs, qe, ts, te;

            /*
            When available, reconstruct raw (non-extended) overlap bounds from aligned marker ordinals.
            This avoids feeding tip-extended bounds into ma_hit_cut, which can over-call containment.
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
                const uint32_t tLen = uint32_t(reads->getReadRawSequenceLength(ad.readIds[1]));

                uint32_t ts0 = tsOriented;
                uint32_t te0 = teOriented;
                if (!ad.isSameStrand) {
                    const auto p = dinara::rcIntervalToForward(tLen, tsOriented, teOriented);
                    ts0 = p.first;
                    te0 = p.second;
                }

                const uint32_t qLen = uint32_t(reads->getReadRawSequenceLength(ad.readIds[0]));
                if (qs0 >= qe0 || qe0 > qLen) return false;
                if (ts0 >= te0 || te0 > tLen) return false;

                outQs = qs0; outQe = qe0; outTs = ts0; outTe = te0;
                return true;
            };

            uint32_t rawQs = ad.qs, rawQe = ad.qe, rawTs = ad.ts, rawTe = ad.te;
            (void)getRawBoundsIfAvailable(rawQs, rawQe, rawTs, rawTe);

            if (!ad.isSameStrand) {
                qs = (int32_t)rawTe < (int32_t)rt.end ? (int32_t)rawQs : (int32_t)rawQs + ((int32_t)rawTe - (int32_t)rt.end);
                qe = (int32_t)rawTs > (int32_t)rt.start ? (int32_t)rawQe : (int32_t)rawQe - ((int32_t)rt.start - (int32_t)rawTs);
                ts = (int32_t)rawQe < (int32_t)rq.end ? (int32_t)rawTs : (int32_t)rawTs + ((int32_t)rawQe - (int32_t)rq.end);
                te = (int32_t)rawQs > (int32_t)rq.start ? (int32_t)rawTe : (int32_t)rawTe - ((int32_t)rq.start - (int32_t)rawQs);

            } else {
                qs = (int32_t)rawTs > (int32_t)rt.start ? (int32_t)rawQs : (int32_t)rawQs + ((int32_t)rt.start - (int32_t)rawTs);
                qe = (int32_t)rawTe < (int32_t)rt.end ? (int32_t)rawQe : (int32_t)rawQe - ((int32_t)rawTe - (int32_t)rt.end);
                ts = (int32_t)rawQs > (int32_t)rq.start ? (int32_t)rawTs : (int32_t)rawTs + ((int32_t)rq.start - (int32_t)rawQs);
                te = (int32_t)rawQe < (int32_t)rq.end ? (int32_t)rawTe : (int32_t)rawTe - ((int32_t)rawQe - (int32_t)rq.end);
            }

            qs = std::max(qs, (int32_t)rq.start);
            qe = std::min(qe, (int32_t)rq.end);
            ts = std::max(ts, (int32_t)rt.start);
            te = std::min(te, (int32_t)rt.end);

            int32_t norm_qs = qs - (int32_t)rq.start;
            int32_t norm_qe = qe - (int32_t)rq.start;
            int32_t norm_ts = ts - (int32_t)rt.start;
            int32_t norm_te = te - (int32_t)rt.start;

            if ((norm_qe - norm_qs >= (int32_t)minLen) && (norm_te - norm_ts >= (int32_t)minLen) && 
                (norm_qe > norm_qs) && (norm_te > norm_ts)) {
                ad.qs = (uint32_t)norm_qs;
                ad.qe = (uint32_t)norm_qe;
                ad.ts = (uint32_t)norm_ts;
                ad.te = (uint32_t)norm_te;
            } else {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonCoverageCut);
            }
        }
    }
}

/*
After structural overlap filtering, hifiasm deletes reads that have no remaining overlaps.
This matches the `rLen==0 -> coverage_cut[i].del=1` behavior in ma_hit_cut/ma_hit_flt/ma_hit_contained_advance.
*/
void Assembler::applyCoverageCutsCleanupThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId r=ReadId(begin); r!=ReadId(end); r++) {
            if (r >= validReadIntervals.size()) continue;
            if (validReadIntervals[r].isDeleted) continue;
            bool hasSurvivingEdge = false;
            
            const auto& table = alignmentTable[OrientedReadId(r, 0).getValue()];
            for (uint32_t alignmentId : table) {
                const auto& ad = alignmentData[alignmentId];
                if(!ad.keptByBothSides()) continue;
                const ReadId other = (ad.readIds[0] == r) ? ad.readIds[1] : ad.readIds[0];
                if (other < validReadIntervals.size() && validReadIntervals[other].isDeleted) continue;
                hasSurvivingEdge = true;
                break;
            }

            if (!hasSurvivingEdge) {
                validReadIntervals[r].isDeleted = true;
            }
        }
    }
}


/*
ma_hit_flt (miniasm/hifiasm) parity: remove internal/non-dovetail overlaps with excessive hangs.
This stage runs after ma_hit_cut, so overlap coordinates are already normalized to valid regions.
*/
void Assembler::filterHangingOverlaps(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Filtering hanging overlaps (ma_hit_flt equivalent)." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    this->hangingFilterMaxHang = maxHang;
    this->hangingFilterMaxHangRate = maxHangRate;
    this->hangingFilterMinOverlap = minOverlapLength;

    setupLoadBalancing(alignmentData.size(), 10000);
    runThreads(&Assembler::filterHangingOverlapsThreadFunction, threadCount);

    setupLoadBalancing(reads->readCount(), 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);

    cout << timestamp << "Hanging overlaps filtered." << endl;
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
            
            if(!ad.keptByBothSides()) continue;

            ReadId qn = ad.readIds[0];
            ReadId tn = ad.readIds[1];

            uint32_t ql, tl;
            if (validReadIntervals.empty()) {
                ql = (uint32_t)reads->getReadRawSequenceLength(qn);
                tl = (uint32_t)reads->getReadRawSequenceLength(tn);
            } else {
                const auto& rq = validReadIntervals[qn];
                const auto& rt = validReadIntervals[tn];
                if (rq.isDeleted || rt.isDeleted) continue;
                ql = rq.end - rq.start;
                tl = rt.end - rt.start;
            }

            uint32_t qs = ad.qs;
            uint32_t qe = ad.qe;
            uint32_t ts = ad.ts;
            uint32_t te = ad.te;
            
            if (qe <= qs || te <= ts) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonHanging);
                continue;
            }

            int32_t tl5, tl3;
            if (!ad.isSameStrand) {
                tl5 = (int32_t)tl - (int32_t)te;
                tl3 = (int32_t)ts;
            } else {
                tl5 = (int32_t)ts;
                tl3 = (int32_t)tl - (int32_t)te;
            }

            int32_t ext5 = ((int32_t)qs < tl5) ? (int32_t)qs : tl5;
            
            int32_t q3Hang = (int32_t)ql - (int32_t)qe;
            int32_t ext3 = (q3Hang < tl3) ? q3Hang : tl3;

            const int result = ma_hit2arc_containment(
                (int32_t)qs, (int32_t)qe, (int32_t)ql,
                (int32_t)ts, (int32_t)te, (int32_t)tl,
                !ad.isSameStrand,
                (int32_t)maxHang,
                maxHangRate,
                (int32_t)minOvlp
            );

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

    const uint64_t readCount = reads->readCount();

    if (!isChimericRead.isOpen) {
        isChimericRead.createNew(largeDataName("IsChimericRead"), largeDataPageSize);
        isChimericRead.resize(readCount);
    }
    std::fill(isChimericRead.begin(), isChimericRead.end(), false);

    chimericReadTmp.assign(readCount, 0);

    setupLoadBalancing(readCount, 64);
    runThreads(&Assembler::detectChimericReadsThreadFunction, threadCount);

    uint64_t chimericCount = 0;
    for (ReadId r = 0; r < readCount; ++r) {
        if (!chimericReadTmp[r]) continue;
        ++chimericCount;
        isChimericRead[r] = true;
        if (r < validReadIntervals.size()) {
            validReadIntervals[r].isDeleted = true;
        }
    }

    if (chimericCount) {
        for (uint64_t alignmentId = 0; alignmentId < alignmentData.size(); ++alignmentId) {
            AlignmentData& ad = alignmentData[alignmentId];
            const ReadId r0 = ad.readIds[0];
            const ReadId r1 = ad.readIds[1];
            if (r0 < readCount && chimericReadTmp[r0]) { ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonChimeric); continue; }
            if (r1 < readCount && chimericReadTmp[r1]) { ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonChimeric); continue; }
        }
    }

    cout << timestamp << "Detected " << chimericCount << " chimeric reads." << endl;
}

void Assembler::detectChimericReadsThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    constexpr double shiftRate = 0.06;
    constexpr float overlapRate = 0.1f;

    struct SubRegion { uint32_t s, e; };
    std::vector<char> xBuf;
    std::vector<char> yBuf;

    uint64_t begin, end;
    while (getNextBatch(begin, end)) {
        for (ReadId qId = ReadId(begin); qId < ReadId(end); ++qId) {
            const uint32_t qLen = uint32_t(reads->getReadRawSequenceLength(qId));
            const OrientedReadId oid(qId, 0);
            if (oid.getValue() >= alignmentTable.size()) continue;

            SubRegion maxLeft{qLen, 0};
            SubRegion maxRight{qLen, 0};

            const auto& table = alignmentTable[oid.getValue()];
            for (uint32_t alignmentId : table) {
                const AlignmentData& ad = alignmentData[alignmentId];
                if (ad.isDeleted()) continue;
                if (isDeletedFromReadPerspective(ad, qId)) continue;

                uint32_t qs, qe;
                if (ad.readIds[0] == qId) { qs = ad.qs; qe = ad.qe; }
                else { qs = ad.ts; qe = ad.te; }

                if (qs == 0) {
                    if (qs < maxLeft.s) maxLeft.s = qs;
                    if (qe > maxLeft.e) maxLeft.e = qe;
                }
                if (qe == qLen) {
                    if (qs < maxRight.s) maxRight.s = qs;
                    if (qe > maxRight.e) maxRight.e = qe;
                }
            }

            if (maxLeft.s == qLen || maxRight.s == qLen) continue;

            uint32_t newLeftE = maxLeft.e;
            uint32_t newRightS = maxRight.s;
            for (uint32_t alignmentId : table) {
                const AlignmentData& ad = alignmentData[alignmentId];
                if (ad.isDeleted()) continue;
                if (isDeletedFromReadPerspective(ad, qId)) continue;

                uint32_t qs, qe;
                if (ad.readIds[0] == qId) { qs = ad.qs; qe = ad.qe; }
                else { qs = ad.ts; qe = ad.te; }

                if (qs == 0 || qe == qLen) continue;
                const uint32_t len = qe - qs;
                if (len == 0) continue;

                if (qs < maxLeft.e && qe > maxLeft.e) {
                    if ((maxLeft.e - qs) > uint32_t(overlapRate * float(len))) {
                        if (qe > maxLeft.e && qe > newLeftE) newLeftE = qe;
                    }
                }
                if (qs < maxRight.s && qe > maxRight.s) {
                    if ((qe - maxRight.s) > uint32_t(overlapRate * float(len))) {
                        if (qs < maxRight.s && qs < newRightS) newRightS = qs;
                    }
                }
            }
            maxLeft.e = newLeftE;
            maxRight.s = newRightS;

            if (maxLeft.e > maxRight.s && (maxLeft.e - maxRight.s) >= uint32_t(double(qLen) * shiftRate)) {
                continue;
            }

            if (maxLeft.e <= maxRight.s) {
                chimericReadTmp[qId] = 1;
                continue;
            }

            const uint32_t intervalS = maxRight.s;
            const uint32_t intervalE = maxLeft.e;

            bool isChimeric = false;
            for (uint32_t alignmentId : table) {
                const AlignmentData& ad = alignmentData[alignmentId];
                if (ad.isDeleted()) continue;
                if (isDeletedFromReadPerspective(ad, qId)) continue;

                ReadId tId;
                uint32_t qs, qe, ts, te;
                bool rev;
                if (ad.readIds[0] == qId) {
                    tId = ad.readIds[1];
                    qs = ad.qs; qe = ad.qe;
                    ts = ad.ts; te = ad.te;
                    rev = !ad.isSameStrand;
                } else {
                    tId = ad.readIds[0];
                    qs = ad.ts; qe = ad.te;
                    ts = ad.qs; te = ad.qe;
                    rev = !ad.isSameStrand;
                }

                if (qs <= intervalS && qe >= intervalE) {
                    if (!boundaryVerify(*reads, intervalS, intervalE, qId, tId, qs, ts, te, rev, xBuf, yBuf)) {
                        isChimeric = true;
                        break;
                    }
                }
            }

            if (isChimeric) {
                chimericReadTmp[qId] = 1;
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
try_rescue_overlaps (hifiasm) parity: resolve directional phasing conflicts.

If `isDeleted0 != isDeleted1`, one read rejected the overlap while the other kept it.
For each read, we collect its conflict overlaps where that read deleted the overlap, find a
maximum-depth consensus interval, and rescue conflicts spanning that interval by clearing the
phase deletion bit only from that read's perspective.
*/
void Assembler::rescuePhasedOverlaps(uint64_t rescueThreshold, uint64_t threadCount)
{
    cout << timestamp << "Rescuing phased overlaps (threshold=" << rescueThreshold << ")..." << endl;
    
    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    
    this->rescuePhasedThreshold = rescueThreshold;
    
    setupLoadBalancing(reads->readCount(), 1);
    runThreads(&Assembler::rescuePhasedOverlapsThreadFunction, threadCount);
    
    cout << timestamp << "Phased overlap rescue complete." << endl;
}

void Assembler::rescuePhasedOverlapsThreadFunction(size_t threadId)
{
    static_cast<void>(threadId);
    const uint64_t rescueThreshold = this->rescuePhasedThreshold;
    
    uint64_t readIdBegin, readIdEnd;
    while(getNextBatch(readIdBegin, readIdEnd)) {
        for(ReadId readId = ReadId(readIdBegin); readId < ReadId(readIdEnd); readId++) {
            
            OrientedReadId orientedReadId(readId, 0);
            
            if (orientedReadId.getValue() >= alignmentTable.size()) continue;
            
            std::vector<uint32_t> conflictAlignments;
            std::vector<std::pair<uint32_t, uint32_t>> conflictIntervals;

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
                
                bool hasConflict = (ad.isDeleted0() != ad.isDeleted1());
                if (!hasConflict) continue;
                
                uint32_t qs, qe;
                bool thisReadDeleted;

                uint32_t qs0 = ad.qs, qe0 = ad.qe, ts0 = ad.ts, te0 = ad.te;
                (void)getRawBoundsIfAvailable(ad, qs0, qe0, ts0, te0);
                
                if (ad.readIds[0] == readId) {
                    qs = qs0;
                    qe = qe0;
                    thisReadDeleted = ad.isDeleted0();
                } else {
                    qs = ts0;
                    qe = te0;
                    thisReadDeleted = ad.isDeleted1();
                }
                
                if (thisReadDeleted) {
                    conflictAlignments.push_back((uint32_t)alignmentId);
                    conflictIntervals.push_back({qs, qe});
                }
            }
            
            if (conflictAlignments.size() < rescueThreshold) continue;
            
            std::vector<std::pair<uint32_t, int>> events;
            events.reserve(conflictIntervals.size() * 2);
            for(const auto& interval : conflictIntervals) {
                events.push_back({interval.first, 1});
                events.push_back({interval.second, -1});
            }
            
            std::sort(events.begin(), events.end(), [](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first < b.first;
                return a.second > b.second;
            });
            
            int dp = 0, max_dp = 0;
            uint32_t start = 0;
            uint32_t max_interval_s = 0, max_interval_e = 0;
            
            for(const auto& ev : events) {
                int old_dp = dp;
                dp += ev.second;
                
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
            
            if ((uint64_t)max_dp >= rescueThreshold) {
                for(size_t j = 0; j < conflictAlignments.size(); j++) {
                    uint32_t alignmentId = conflictAlignments[j];
                    uint32_t qs = conflictIntervals[j].first;
                    uint32_t qe = conflictIntervals[j].second;
                    
                    if (qs <= max_interval_s && qe >= max_interval_e) {
                        AlignmentData& ad = alignmentData[alignmentId];
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
    This is a direct translation of hifiasm's `ma_hit2arc` containment checks:

    Given a query interval [qs, qe) on a query of length ql, and a target interval [ts, te) on a target
    of length tl, compute:
      - query 5' overhang:  qs
      - query 3' overhang:  ql - qe
      - target 5'/3' overhangs in query orientation (depends on whether the target is reverse-complemented)

    We then:
      1) reject internal matches where both sides have too-large overhangs (absolute or fractional),
      2) detect containment (query-in-target or target-in-query),
      3) ensure the effective overlap length meets the minimum requirement.

    Return values match the local callers:
      0 = dovetail, 1 = query contained, 2 = target contained, -1 = internal, -2 = too short.
    */
    int32_t tl5, tl3;
    if (isReverse) {
        tl5 = tl - te;
        tl3 = ts;
    } else {
        tl5 = ts;
        tl3 = tl - te;
    }

    int32_t ext5 = (qs < tl5) ? qs : tl5;
    int32_t ext3 = ((ql - qe) < tl3) ? (ql - qe) : tl3;

    if (ext5 > max_hang || ext3 > max_hang) {
        return -1;
    }
    
    int32_t qOverlapLen = qe - qs;
    int32_t tOverlapLen = te - ts;
    if (qOverlapLen < (qOverlapLen + ext5 + ext3) * int_frac) {
        return -1;
    }
    if (tOverlapLen < (tOverlapLen + ext5 + ext3) * int_frac) {
        return -1;
    }

    if (qs <= tl5 && (ql - qe) <= tl3) {
        return 1;
    }
    if (qs >= tl5 && (ql - qe) >= tl3) {
        return 2;
    }

    if (qOverlapLen + ext5 + ext3 < min_ovlp || tOverlapLen + ext5 + ext3 < min_ovlp) {
        return -2;
    }

    return 0;
}


void Assembler::removeContainedReads(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Removing contained reads (ma_hit_contained_advance equivalent)..." << endl;

    /*
    High-level algorithm (mirrors hifiasm, adapted to Dinara's storage):

    1) Work only on reads that are still alive after prior filters (coverage cut + other overlap filters).
    2) For each read qn, traverse its adjacency list via `alignmentTable[OrientedReadId(qn,0)]`.
       We only consider overlaps that are still kept from both ends (`keptByBothSides()`), matching hifiasm's
       `h->del == 0` constraint.
    3) Reconstruct the overlap coordinates from qn's perspective and call `ma_hit2arc_containment`.
       If it reports containment, mark the contained read deleted in `validReadIntervals`, record its parent
       in `containmentParent`, and immediately delete all its incident overlaps (hifiasm `delete_all_edges`).
    4) After the sweep, compress `containmentParent` chains to the final container (hifiasm `transfor_R_to_U`).
    5) Finally, drop reads that became isolated (no surviving overlaps) using the existing cleanup worker.
    */
    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    if (!containmentParent->isOpen) {
        containmentParent->createNew(largeDataName("ContainmentParent"), largeDataPageSize);
        containmentParent->resize(reads->readCount());
    }
    std::fill(containmentParent->begin(), containmentParent->end(), ReadId(invalidReadId));

    uint64_t containedReadCount = 0;

    auto deleteAllEdgesForRead = [&](ReadId r) {
        const OrientedReadId oid(r, 0);
        if (oid.getValue() >= alignmentTable.size()) return;
        const auto& table = alignmentTable[oid.getValue()];
        for (const uint32_t alignmentId : table) {
            alignmentData[alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonContained);
        }
    };

    const uint64_t readCount = reads->readCount();
    for (ReadId qn = 0; qn < readCount; ++qn) {
        if (qn >= validReadIntervals.size() || validReadIntervals[qn].isDeleted) continue;

        const auto& vrQ = validReadIntervals[qn];
        const int32_t ql = int32_t(vrQ.end - vrQ.start);
        if (ql <= 0) continue;

        const OrientedReadId oid(qn, 0);
        if (oid.getValue() >= alignmentTable.size()) continue;

        const auto& table = alignmentTable[oid.getValue()];
        for (uint32_t alignmentId : table) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if(!ad.keptByBothSides()) continue;

            const ReadId tn = (ad.readIds[0] == qn) ? ad.readIds[1] : ad.readIds[0];
            if (tn >= validReadIntervals.size() || validReadIntervals[tn].isDeleted) continue;

            const auto& vrT = validReadIntervals[tn];
            const int32_t tl = int32_t(vrT.end - vrT.start);
            if (tl <= 0) continue;

            const bool rev = !ad.isSameStrand;
            int32_t qs = 0, qe = 0, ts = 0, te = 0;
            if (ad.readIds[0] == qn) {
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
            if (qs < 0 || qe < 0 || ts < 0 || te < 0) continue;
            if (qs >= qe || ts >= te) continue;
            if (qs > ql || qe > ql || ts > tl || te > tl) continue;

            const int result = ma_hit2arc_containment(
                qs, qe, ql,
                ts, te, tl,
                rev,
                int32_t(maxHang),
                maxHangRate,
                int32_t(minOverlapLength)
            );

            if (result == 1) {
                if (!validReadIntervals[qn].isDeleted) {
                    validReadIntervals[qn].isDeleted = true;
                    (*containmentParent)[qn] = tn;
                    deleteAllEdgesForRead(qn);
                    ++containedReadCount;
                }
                break;
            } else if (result == 2) {
                if (!validReadIntervals[tn].isDeleted) {
                    validReadIntervals[tn].isDeleted = true;
                    (*containmentParent)[tn] = qn;
                    deleteAllEdgesForRead(tn);
                    ++containedReadCount;
                }
            }
        }
    }

    for (ReadId r = 0; r < readCount; ++r) {
        if ((*containmentParent)[r] == ReadId(invalidReadId)) continue;
        ReadId root = (*containmentParent)[r];
        while (root != ReadId(invalidReadId) && (*containmentParent)[root] != ReadId(invalidReadId)) {
            root = (*containmentParent)[root];
        }
        (*containmentParent)[r] = root;
    }

    setupLoadBalancing(readCount, 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);

    cout << timestamp << "Identified " << containedReadCount << " contained reads." << endl;
}

void Assembler::flagContainedReads(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Flagging contained reads (diagnostic-only, does not remove overlaps)..." << endl;

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // We will write into read flags.
    reads->checkReadFlagsAreOpenForWriting();

    if(validReadIntervals.empty()) {
        throw runtime_error("flagContainedReads requires validReadIntervals (run filterLocalSegments/applyCoverageCuts first).");
    }

    if(!containmentParent->isOpen) {
        containmentParent->createNew(largeDataName("ContainmentParent"), largeDataPageSize);
        containmentParent->resize(reads->readCount());
    }
    std::fill(containmentParent->begin(), containmentParent->end(), ReadId(invalidReadId));

    const uint64_t readCount = reads->readCount();

    // Clear previous flags.
    for(ReadId r=0; r<readCount; ++r) {
        reads->setContainedFlag(r, false);
    }

    uint64_t containedReadCount = 0;

    // Same containment test as removeContainedReads, but we only record the result.
    for(ReadId qn = 0; qn < readCount; ++qn) {
        if(qn >= validReadIntervals.size() || validReadIntervals[qn].isDeleted) continue;
        if(reads->getFlags(qn).isContained) continue;

        const auto& vrQ = validReadIntervals[qn];
        const int32_t ql = int32_t(vrQ.end - vrQ.start);
        if(ql <= 0) continue;

        const OrientedReadId oid(qn, 0);
        if(oid.getValue() >= alignmentTable.size()) continue;

        const auto& table = alignmentTable[oid.getValue()];
        for(uint32_t alignmentId : table) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if(!ad.keptByBothSides()) continue;

            const ReadId tn = (ad.readIds[0] == qn) ? ad.readIds[1] : ad.readIds[0];
            if(tn >= validReadIntervals.size() || validReadIntervals[tn].isDeleted) continue;
            if(reads->getFlags(tn).isContained) continue;

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
                if(!reads->getFlags(qn).isContained) {
                    reads->setContainedFlag(qn, true);
                    (*containmentParent)[qn] = tn;
                    ++containedReadCount;
                }
                break;
            } else if(result == 2) {
                if(!reads->getFlags(tn).isContained) {
                    reads->setContainedFlag(tn, true);
                    (*containmentParent)[tn] = qn;
                    ++containedReadCount;
                }
            }
        }
    }

    // Compress containment chains.
    for(ReadId r = 0; r < readCount; ++r) {
        if((*containmentParent)[r] == ReadId(invalidReadId)) continue;
        ReadId root = (*containmentParent)[r];
        while(root != ReadId(invalidReadId) && (*containmentParent)[root] != ReadId(invalidReadId)) {
            root = (*containmentParent)[root];
        }
        (*containmentParent)[r] = root;
    }

    cout << timestamp << "Flagged " << containedReadCount << " contained reads." << endl;
}
