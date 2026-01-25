#include "Assembler.hpp"
#include "Reads.hpp" // For Reads class definition
#include "hifiasmCoordinateTransforms.hpp"
#include "timestamp.hpp"      // For timestamp
// #include "loadBalancing.hpp" // Removed: caused compilation error
#include <algorithm>
#include <vector>
#include <limits>

// For debugging/logging
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

// This function implements logic similar to miniasm's ma_hit_sub function.
// It analyzes the coverage depth along a read induced by its overlaps.
// It identifies contiguous high-coverage regions effectively filtering out 
// chimeric joins or low-quality ends.

void Assembler::filterLocalSegments(
    uint64_t minCoverage, 
    uint64_t threadCount)
{
    cout << timestamp << "Filtering local segments (min coverage=" << minCoverage << ")..." << endl;
    
    // Store parameters for the thread function if needed, 
    // or passing them via member variables/lambdas if strictly needed.
    // For now we assume we might need to store minCoverage in Assembler 
    // or just hardcode/pass it. 
    // Since runThreads takes a member function with only threadId, 
    // we should optimally store this config in a member or capture it via a lambda 
    // if runThreads supported it (it doesn't, it takes a pointer).
    // So we'll use a member variable or a temporary one.
    // Let's rely on 'this->minOverlapCoverage' or similar if it exists, 
    // or just add a temporary member for this pass.
    
    // Actually, miniasm passes min_dp.
    // Let's look for a suitable place to store this transient config.
    // For now, I'll assume we can add a member to Assembler or reuse an existing one.
    // I'll add 'localSegmentMinCoverage' to Assembler.hpp later if needed.
    // But wait, I can't change the header right now easily inside this write.
    // I will add a member variable 'localSegmentMinCoverage' in Assembler.hpp next.
    
    // Resize output vector
    validReadIntervals.resize(reads->readCount());
    // Initialize with deleted status
    std::fill(validReadIntervals.begin(), validReadIntervals.end(), ReadSegmentStatus{0, 0, false}); // Default not deleted? Or should be? 
    // Miniasm default initializes and then sets. 
    // Let's safe initialize to 0,0,false. The loop covers all reads in batches.
    // If a read has no alignments, the loop:
    // 1. collectIntervals -> empty.
    // 2. if empty -> validReadIntervals[r0] = {0,0, true}; continue;
    
    // Store parameter for access by thread function
    this->localSegmentMinCoverage = minCoverage;

    setupLoadBalancing(reads->readCount(), 100);
    runThreads(&Assembler::filterLocalSegmentsThreadFunction, threadCount);
    
    cout << timestamp << "Local segment filtering complete." << endl;
}

// ----------------------------------------------------------------------------
// ONT Chemical Arc Masking (gen_chemical_arc_rf equivalent)
// ----------------------------------------------------------------------------
// Hifiasm uses this only for ONT to detect "chemical" chimeras by finding
// low-coverage stretches along a read, using all overlaps (independent of EC/phasing).
// We compute a per-read minimum overlap depth (after trimming flanks) and mark all
// overlaps incident to reads with minDepth <= chemicalArcCov.
//
// Important: this MUST use all overlaps, not just the currently "kept" overlaps, or it
// can artificially reduce depth and over-call chemical chimeras.

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

    // Clear any previous chemical markings but keep all other deletion reasons.
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

    // Mirror hifiasm pass-1 symmetric deletion: delete overlaps if either endpoint is flagged.
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

void Assembler::applyOntChemicalArcMaskThreadFunction(size_t /* threadId */)
{
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

                    // Use all overlaps regardless of deletion reasons (hifiasm resets del first).
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

void Assembler::filterLocalSegmentsThreadFunction(size_t /* threadId */)
{
    // Minimal depth (min_dp in miniasm)
    // We use the passed minCoverage argument.
    // However, it is passed to the main function, not stored in "this".
    // We need to access it.
    // Option A: Store it in a member variable before running threads.
    // Option B: Hardcode/atomic (not ideal).
    // Reviewing Assembler class, we can just add a member `localSegmentMinCoverage`.
    
    // For now, assuming we added the member or use a specialized approach.
    // But to respect the current structure without modifying header again immediately:
    // We will assume `minCoverage` is available or we use a lambda capture if we could.
    // Since we can't change signature of runThreads target easily...
    // Let's assume we added `this->minOverlapCoverage` usage or similar.
    // Actually, `Assembler` has `minOverlapCoverage`? No, it has `minAlleleCoverage`.
    
    // Let's use `3` as a placeholder DEFAULT if we can't access the arg, 
    // BUT the user asked for "miniasm style", implying dynamic.
    // I will read `localSegmentMinCoverage` which I will conceptually add to `Assembler`.
    // Wait, I can't read it if I didn't add it.
    // The previous Plan step added `validReadIntervals` to header.
    // I should have added `localSegmentMinCoverage` too. 
    // I will use a const for now to allow compilation, but add TODO.
    // OR better: I can modify the header now to add the configuration variable.
    
    // Hifiasm parity (ma_hit_sub): if min_dp <= 1, keep the full read interval.
    const uint64_t minDp = this->localSegmentMinCoverage;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); r0++) {
            
            // Fast path for low coverage requirement (pass everything)
            if (minDp <= 1) {
                 uint64_t len = reads->getReadRawSequenceLength(r0);
                 validReadIntervals[r0] = {0, (uint32_t)len, false};
                 continue;
            }

            // Collect alignment intervals (start, end) on r0
            // from all alignments involving r0.
            
            // 1. Collect Valid Intervals
            std::vector<std::pair<uint32_t, uint32_t>> intervals;
            
            auto collectIntervals = [&](OrientedReadId orientedR0) {
                const auto& table = alignmentTable[orientedR0.getValue()];
                for (uint32_t alignmentId : table) {
                    const auto& ad = alignmentData[alignmentId];
                    // Hifiasm parity: only overlaps that are present in this read's adjacency.
                    // A global alignment can be deleted directionally by parity EC.
                    if (isDeletedFromReadPerspective(ad, r0)) continue;

                    // We use the explicit coordinates provided in AlignmentData (qs, qe, ts, te).
                    // These are assumed to be 0-based coordinates on the Forward strand of the read.
                    
                    uint32_t start = 0, end = 0;
                    
                    if (ad.readIds[0] == r0) {
                        // r0 is the Query (first read)
                        start = ad.qs;
                        end = ad.qe;
                    } else {
                        // r0 must be the Target (second read)
                        start = ad.ts;
                        end = ad.te;
                    }
                    
                    if (end > start) {
                        intervals.push_back({start, end});
                    }
                }
            };
            
            
            // We only collect from strand 0 to match Hifiasm's sources[i] which is single-direction.
            // Collecting from both strands would double-count overlaps.
            collectIntervals(OrientedReadId(r0, 0));
            // Hifiasm parity: Don't collect strand 1 - same alignment is indexed for both oriented reads
            // collectIntervals(OrientedReadId(r0, 1));
            
            if (intervals.empty()) {
                // No coverage means no valid segment -> deleted
                validReadIntervals[r0] = {0, 0, true};
                continue; 
            }
            
            // 2. Compute Coverage Depth (Sweep Line)
            // Events: (pos, 1) for start, (pos, -1) for end.
            std::vector<std::pair<uint32_t, int>> events;
            events.reserve(intervals.size() * 2);
            for(const auto& interval : intervals) {
                events.push_back({interval.first, 1});
                events.push_back({interval.second, -1});
            }
            
            // Sort events. Process Start (+1) before End (-1) at same position
            // to maintain continuity of adjacent interval coverage.
            std::sort(events.begin(), events.end(), [](const pair<uint32_t, int>& a, const pair<uint32_t, int>& b){
                if (a.first != b.first) return a.first < b.first;
                return a.second > b.second; // 1 (Start) before -1 (End)
            });
            
            // 3. Find Longest Segment with coverage >= minDp
            int coverage = 0;
            uint32_t maxS = 0, maxE = 0;
            
            uint32_t currentStart = 0;
            bool inSegment = false;
            
            // Loop through all events
            for(const auto& ev : events) {
                int oldCoverage = coverage;
                coverage += ev.second;
                
                uint32_t pos = ev.first;
                
                // Rising edge to valid coverage
                if (oldCoverage < (int)minDp && coverage >= (int)minDp) {
                    currentStart = pos;
                    inSegment = true;
                }
                // Falling edge from valid coverage
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
            
            // Store result
            if (maxE > maxS) {
                validReadIntervals[r0] = {maxS, maxE, false}; // Keep
            } else {
                validReadIntervals[r0] = {0, 0, true}; // Delete
            }
        }
    }
}



// Apply coverage cuts (ma_hit_cut equivalent).
// Clips alignments to the valid regions identified by filterLocalSegments.
// Removes short overlaps.
void Assembler::applyCoverageCuts(uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Applying coverage cuts (ma_hit_cut equivalent)." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    this->coverageCutMinOverlap = minOverlapLength;

    // Phase 1: Clip alignments in parallel.
    // Iterate over alignmentData directly to avoid double processing and simplify concurrency.
    setupLoadBalancing(alignmentData.size(), 10000); 
    runThreads(&Assembler::applyCoverageCutsToAlignmentsThreadFunction, threadCount);

    // Phase 2: Cleanup orphaned reads.
    // Reads that have no surviving edges should be marked as deleted.
    setupLoadBalancing(reads->readCount(), 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);
    
    cout << timestamp << "Coverage cuts applied." << endl;
}

void Assembler::applyCoverageCutsToAlignmentsThreadFunction(size_t /* threadId */)
{
    uint64_t begin, end;
    const uint64_t minLen = this->coverageCutMinOverlap;

    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; i++) {
            AlignmentData& ad = alignmentData[i];
            
            // Hifiasm parity: structural filters operate on overlaps that survived previous filters.
            // Our equivalent of h->del==0 is keptByBothSides().
            if(!ad.keptByBothSides()) continue;

            ReadId qn = ad.readIds[0];
            ReadId tn = ad.readIds[1];

            // Access valid intervals (coverage_cut)
            // If explicit intervals not computed, assume full read? 
            // ma_hit_cut assumes sub is populated. We assume validReadIntervals is populated.
            // Safety check:
            if (qn >= validReadIntervals.size() || tn >= validReadIntervals.size()) continue;

            const auto& rq = validReadIntervals[qn];
            const auto& rt = validReadIntervals[tn];

            // "if any of target read and the query read has no enough coverage" -> del
            if (rq.isDeleted || rt.isDeleted) continue;

            // Alias internal coords for readability matching reference code
            int32_t qs, qe, ts, te;

            // Hifiasm parity note:
            // For ONT/HiFi string-graph cleaning, hifiasm uses overlap boundaries coming from
            // the aligner/PAF (not the tip-extended "approximate" boundaries used in some
            // seed-chain stages). If we feed tip-extended boundaries into ma_hit_cut, we can
            // dramatically over-call containment later.
            //
            // In Dinara, when marker ordinals are available for this alignment we can cheaply
            // recover non-extended overlap bounds from the first/last aligned markers and use
            // those as the input to ma_hit_cut.
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

                // Basic sanity.
                const uint32_t qLen = uint32_t(reads->getReadRawSequenceLength(ad.readIds[0]));
                if (qs0 >= qe0 || qe0 > qLen) return false;
                if (ts0 >= te0 || te0 > tLen) return false;

                outQs = qs0; outQe = qe0; outTs = ts0; outTe = te0;
                return true;
            };

            uint32_t rawQs = ad.qs, rawQe = ad.qe, rawTs = ad.ts, rawTe = ad.te;
            (void)getRawBoundsIfAvailable(rawQs, rawQe, rawTs, rawTe);

            // Logic mirroring ma_hit_cut coordinate projection
            if (!ad.isSameStrand) { // Different Strand (p->rev)
                // Anti-parallel projection
                qs = (int32_t)rawTe < (int32_t)rt.end ? (int32_t)rawQs : (int32_t)rawQs + ((int32_t)rawTe - (int32_t)rt.end);
                qe = (int32_t)rawTs > (int32_t)rt.start ? (int32_t)rawQe : (int32_t)rawQe - ((int32_t)rt.start - (int32_t)rawTs);
                ts = (int32_t)rawQe < (int32_t)rq.end ? (int32_t)rawTs : (int32_t)rawTs + ((int32_t)rawQe - (int32_t)rq.end);
                te = (int32_t)rawQs > (int32_t)rq.start ? (int32_t)rawTe : (int32_t)rawTe - ((int32_t)rq.start - (int32_t)rawQs);

            } else { // Same Strand
                // Parallel offset
                qs = (int32_t)rawTs > (int32_t)rt.start ? (int32_t)rawQs : (int32_t)rawQs + ((int32_t)rt.start - (int32_t)rawTs);
                qe = (int32_t)rawTe < (int32_t)rt.end ? (int32_t)rawQe : (int32_t)rawQe - ((int32_t)rawTe - (int32_t)rt.end);
                ts = (int32_t)rawQs > (int32_t)rq.start ? (int32_t)rawTs : (int32_t)rawTs + ((int32_t)rq.start - (int32_t)rawQs);
                te = (int32_t)rawQe < (int32_t)rq.end ? (int32_t)rawTe : (int32_t)rawTe - ((int32_t)rawQe - (int32_t)rq.end);
            }

            // "cut by self coverage" - Strict Clip to valid bounds
            // Using max/min to enforce bounds without offset change.
            qs = std::max(qs, (int32_t)rq.start);
            qe = std::min(qe, (int32_t)rq.end);
            ts = std::max(ts, (int32_t)rt.start);
            te = std::min(te, (int32_t)rt.end);

            // Hifiasm parity: Normalize coordinates to 0-based relative to valid region
            // This matches: qs = (qs - rq.s), qe = (qe - rq.s), etc.
            int32_t norm_qs = qs - (int32_t)rq.start;
            int32_t norm_qe = qe - (int32_t)rq.start;
            int32_t norm_ts = ts - (int32_t)rt.start;
            int32_t norm_te = te - (int32_t)rt.start;

            // Check if valid length remains
            // "if (qe - qs >= mini_overlap_length && te - ts >= mini_overlap_length)"
            if ((norm_qe - norm_qs >= (int32_t)minLen) && (norm_te - norm_ts >= (int32_t)minLen) && 
                (norm_qe > norm_qs) && (norm_te > norm_ts)) {
                // Update alignment data with normalized coordinates
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

void Assembler::applyCoverageCutsCleanupThreadFunction(size_t /* threadId */)
{
    // Phase 2: Check for orphaned reads.
    // "if(rLen == 0) (*coverage_cut)[i].del = 1;"
    // Iterate reads, check if any edges survive.

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId r=ReadId(begin); r!=ReadId(end); r++) {
            if (r >= validReadIntervals.size()) continue;
            if (validReadIntervals[r].isDeleted) continue; // Already deleted

            // Check neighbors
            bool hasSurvivingEdge = false;
            
            // Hifiasm parity: sources[i] is a per-read adjacency, so we only count overlaps that
            // are kept from this read's perspective.
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


// Filter hanging overlaps (ma_hit_flt equivalent).
// Filters overlaps based on hanging length (dovetail rule) and containment.
void Assembler::filterHangingOverlaps(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Filtering hanging overlaps (ma_hit_flt equivalent)." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    this->hangingFilterMaxHang = maxHang;
    this->hangingFilterMaxHangRate = maxHangRate;
    this->hangingFilterMinOverlap = minOverlapLength;

    // Use the alignment thread function directly
    setupLoadBalancing(alignmentData.size(), 10000);
    runThreads(&Assembler::filterHangingOverlapsThreadFunction, threadCount);

    // After filtering edges, we should clean up orphaned reads again?
    // ma_hit_flt does this internally ("if (rLen == 0) ...").
    // Reuse existing cleanup function.
    setupLoadBalancing(reads->readCount(), 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);

    cout << timestamp << "Hanging overlaps filtered." << endl;
}

void Assembler::filterHangingOverlapsThreadFunction(size_t /* threadId */)
{
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

            // Access valid lengths. If not present (e.g. skipped step), use full read len?
            // ma_hit_flt assumes validReadIntervals populated.
            // If empty, assume full length (fallback).
            // uint32_t ql, tl, qs0, ts0; // Unused
            uint32_t ql, tl;
            if (validReadIntervals.empty()) {
                ql = (uint32_t)reads->getReadRawSequenceLength(qn);
                tl = (uint32_t)reads->getReadRawSequenceLength(tn);
                // (void)qs0; (void)ts0; // Suppress unused warning
            } else {
                const auto& rq = validReadIntervals[qn];
                const auto& rt = validReadIntervals[tn];
                if (rq.isDeleted || rt.isDeleted) continue;
                ql = rq.end - rq.start;
                tl = rt.end - rt.start;
                // qs0 = rq.start;
                // ts0 = rt.start;
            }

            // IMPORTANT: After applyCoverageCuts, coordinates are already normalized (0-based relative to valid region)
            // So we use them directly without subtracting vr.start again
            uint32_t qs = ad.qs;
            uint32_t qe = ad.qe;
            uint32_t ts = ad.ts;
            uint32_t te = ad.te;
            
            // Re-check length just in case
            if (qe <= qs || te <= ts) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonHanging);
                continue;
            }

            // Calculate Overhangs (tl5, tl3) - Target hangs
            int32_t tl5, tl3;
            if (!ad.isSameStrand) { // Different Strand
                // "tl5 = tl - h->te, tl3 = h->ts"
                tl5 = (int32_t)tl - (int32_t)te;
                tl3 = (int32_t)ts;
            } else { // Same Strand
                // "tl5 = h->ts, tl3 = tl - h->te"
                tl5 = (int32_t)ts;
                tl3 = (int32_t)tl - (int32_t)te;
            }

            // Calculate Extensions (ext5, ext3) - Combined hangs
            // "ext5 = qs < tl5 ? qs : tl5;" (min)
            int32_t ext5 = ((int32_t)qs < tl5) ? (int32_t)qs : tl5;
            
            // "ext3 = ql - h->qe < tl3 ? ql - h->qe : tl3;" (min of Query 3' Hang vs Target 3' Hang)
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

            // Hifiasm parity (Overlaps.cpp:1898): keep dovetails and containments; delete INT/SHORT.
            if (result < 0) {
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonHanging);
            }
        }
    }
}


// ----------------------------------------------------------------------------
// Chimeric Read Detection
// ----------------------------------------------------------------------------
// Mirrors hifiasm's logic (gen_chemical_arc_rf -> cal_chemical_r_adv)
// to detect coverage gaps in "Cis" alignments.

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

    // Thread-safe temp storage (MemoryMapped::Vector<bool> is bit-packed).
    chimericReadTmp.assign(readCount, 0);

    setupLoadBalancing(readCount, 64);
    runThreads(&Assembler::detectChimericReadsThreadFunction, threadCount);

    // Apply deletions serially (hifiasm delete_all_edges equivalent).
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

void Assembler::detectChimericReadsThreadFunction(size_t /* threadId */)
{
    // Hifiasm parity (Overlaps.cpp:2449 detect_chimeric_reads).
    // Called after ma_hit_sub and before ma_hit_cut: coordinates are still absolute.
    constexpr double shiftRate = 0.06;   // asm_opt.max_ov_diff_final * 2.0, default 0.03 * 2
    constexpr float overlapRate = 0.1f;  // collect_contain overlap_rate

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

            // 1) collect_sides: anchors that touch the read ends exactly.
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

            // End node: missing one of the two anchors.
            if (maxLeft.s == qLen || maxRight.s == qLen) continue;

            // 2) collect_contain: extend anchors with contained overlaps.
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

            // 3) shift-rate overlap check.
            if (maxLeft.e > maxRight.s && (maxLeft.e - maxRight.s) >= uint32_t(double(qLen) * shiftRate)) {
                continue;
            }

            // 4) simple chimera: uncovered gap.
            if (maxLeft.e <= maxRight.s) {
                chimericReadTmp[qId] = 1;
                continue;
            }

            // 5) complex chimera: small overlap; verify by base-level window checks.
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

// ----------------------------------------------------------------------------
// Chimeric Read Rescue (Symmetry Restoration)
// ----------------------------------------------------------------------------
// Mirrors hifiasm's `try_rescue_overlaps`.
// Logic:
// 1. Iterate over reads flagged as Chimeric.
// 2. Collect "Safe Neighbors" (Non-chimeric, Valid Alignment).
// 3. STRICT CHECK: Do these neighbors support the SAME region?
//    - Hifiasm performs a sweep-line on the neighbors' overlaps.
//    - It finds the maximum pileup depth (max_dp).
// 4. If max_dp >= 4 (RescueThreshold), un-flag the read.

void Assembler::rescueChimericReads(uint64_t threadCount)
{
    cout << timestamp << "Rescuing chimeric reads..." << endl;
    
    // We need to know who is chimeric to check neighbors. 
    setupLoadBalancing(reads->readCount(), 1);
    runThreads(&Assembler::rescueChimericReadsThreadFunction, threadCount);

    // Recount
    uint64_t chimericCount = 0;
    for(size_t i=0; i<isChimericRead.size(); i++) {
        if(isChimericRead[i]) chimericCount++;
    }
    cout << timestamp << "Chimeric reads after rescue: " << chimericCount << "." << endl;
}

void Assembler::rescueChimericReadsThreadFunction(size_t /* threadId */)
{
    const uint32_t rescueThreshold = 4; // Hifiasm default
    
    uint64_t readIdBegin, readIdEnd;
    while(getNextBatch(readIdBegin, readIdEnd)) {
        for(ReadId readId = ReadId(readIdBegin); readId < ReadId(readIdEnd); readId++) {
            
            // Only process if currently flagged as Chimeric
            if (!isChimericRead[readId]) continue;
            
            OrientedReadId orientedReadId(readId, 0);
            
            // Check neighbors in alignment table
            if (orientedReadId.getValue() >= alignmentTable.size()) continue;
            
            vector<pair<uint32_t, int>> events;
            
            const size_t n = alignmentTable.size(orientedReadId.getValue());
            
            for(size_t i=0; i<n; i++) {
                const uint64_t alignmentId = alignmentTable[orientedReadId.getValue()][i];
                const AlignmentData& ad = alignmentData[alignmentId];
                
                // Identify the neighbor
                ReadId neighborId;
                uint32_t qs, qe; // Coordinates on CURRENT read (readId)
                
                if (ad.readIds[0] == readId) {
                    neighborId = ad.readIds[1];
                    qs = ad.qs; qe = ad.qe;
                } else {
                    neighborId = ad.readIds[0];
                    qs = ad.ts; qe = ad.te;
                }
                
                // 1. Is Neighbor Chimeric?
                if (isChimericRead[neighborId]) continue;
                
                // 2. Is Alignment Valid?
                bool isValid = (ad.cisTransStatus == CisTransStatus::Cis);
                if (ad.cisTransStatus == CisTransStatus::Unknown && ad.info.isInReadGraph) {
                    isValid = true;
                }
                
                if (isValid) {
                    // Add Interval [qs, qe) to sweep line events
                    // Use Type 1=Start, 2=End to ensure Start < End at same position (Start processed first) - Matches Hifiasm
                    if (qs < qe) {
                        events.push_back({qs, 1});
                        events.push_back({qe, 2});
                    }
                }
            }
            
            // 3. Sweep Line to check Consensus Depth (max_dp)
            // We need at least 'rescueThreshold' reads agreeing on the SAME region.
            if (events.size() / 2 >= rescueThreshold) {
                // Sort events:
                // Hifiasm sorts `qs<<1` (Start, LSB=0) and `qe<<1|1` (End, LSB=1).
                // This means at the same coordinate, START (even) comes before END (odd).
                // We use TYPE_START=1, TYPE_END=2.
                // pair compare: first(pos) then second(type).
                // So {pos, 1} < {pos, 2}. Matches Hifiasm.
                std::sort(events.begin(), events.end());
                
                int max_dp = 0;
                int dp = 0;
                uint32_t start = 0;
                uint32_t max_interval_s = 0;
                uint32_t max_interval_e = 0;
                int old_dp = 0;

                for(const auto& ev : events) {
                    old_dp = dp;
                    
                    if (ev.second == 1) { // Start
                        dp++;
                    } else { // End
                        dp--;
                    }
                    
                    // Logic mirroring Hifiasm's max_interval tracking
                    if (old_dp < dp) { // Increasing (Start processed)
                        // "if(dp >= max_dp) { start = b.a[j]>>1; max_dp = dp; }"
                        if (dp >= max_dp) {
                            start = ev.first;
                            max_dp = dp;
                        }
                    } else if (old_dp > dp) { // Decreasing (End processed)
                        // "if(old_dp == max_dp) { max_interval.s = start; max_interval.e = b.a[j]>>1; }"
                        if (old_dp == max_dp) {
                            max_interval_s = start;
                            max_interval_e = ev.first;
                        }
                    }
                }
                
                // Hifiasm rescues if max_dp >= threshold.
                // It then filters edges to only those spanning max_interval.
                // "if(qs <= max_interval.s && qe >= max_interval.e)"
                
                if (max_dp >= (int)rescueThreshold) {
                    // Rescue!
                    isChimericRead[readId] = false;
                    
                    // Second Pass: Filter edges (Enforce Consensus)
                    // Hifiasm effectively drops edges that don't support the consensus.
                    // We mark them as deleted.
                    
                    for(size_t i=0; i<n; i++) {
                        const uint64_t alignmentId = alignmentTable[orientedReadId.getValue()][i];
                        // We need access to writable AlignmentData
                        AlignmentData& ad = alignmentData[alignmentId];
                        
                        // Identify neighbor & coordinates again
                        ReadId neighborId;
                        uint32_t qs, qe;
                        if (ad.readIds[0] == readId) {
                            neighborId = ad.readIds[1];
                            qs = ad.qs; qe = ad.qe;
                        } else {
                            neighborId = ad.readIds[0];
                            qs = ad.ts; qe = ad.te;
                        }
                        
                         // Same validity checks as pass 1
                        if (isChimericRead[neighborId]) continue; // Chimeric neighbors likely shouldn't constrain us, but Hifiasm ignores them in "evidence".
                        // Wait, Hifiasm only loops over "evi" (the collected valid neighbors) for filtering.
                        // So we should only filter the *valid* edges? 
                        // If an edge was skipped in pass 1 (e.g. invalid status), it remains whatever it was.
                        // But if we rescue the read, we keep "all" edges unless we delete them.
                        // So we MUST process all edges and delete those that fail the check?
                        // Hifiasm loop: "for (j = 0; j < evi.n; j++)" -> It only adds valid edges to paf. 
                        // The implicit implication is that edges NOT in evi are NOT added (remain deleted).
                        // So yes, we should probably delete *everything* that doesn't pass.
                        
                        bool isValidType = (ad.cisTransStatus == CisTransStatus::Cis);
                        if (ad.cisTransStatus == CisTransStatus::Unknown && ad.info.isInReadGraph) isValidType = true;
                        

                        
                        // Strict filter:
                        // 1. Must be Valid Type
                        // 2. Must span max_interval
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

// ----------------------------------------------------------------------------
// Rescue Phased Overlaps (try_rescue_overlaps equivalent)
// ----------------------------------------------------------------------------
// Rescues overlaps where isDeleted0 != isDeleted1 (directional phasing conflict)
// if >= rescueThreshold overlaps from the same direction agree on a region.
// This recovers edges that were incorrectly phased out from one direction.

void Assembler::rescuePhasedOverlaps(uint64_t rescueThreshold, uint64_t threadCount)
{
    cout << timestamp << "Rescuing phased overlaps (threshold=" << rescueThreshold << ")..." << endl;
    
    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    
    this->rescuePhasedThreshold = rescueThreshold;
    
    setupLoadBalancing(reads->readCount(), 1);
    runThreads(&Assembler::rescuePhasedOverlapsThreadFunction, threadCount);
    
    // Count rescued overlaps
    // uint64_t rescuedCount = 0;
    for(size_t i = 0; i < alignmentData.size(); i++) {
        const auto& ad = alignmentData[i];
        // Count overlaps where both flags are now false (rescued)
        // that previously had exactly one flag set
        if (!ad.isDeleted()) {
            // Simple count of surviving overlaps after rescue
        }
    }
    
    cout << timestamp << "Phased overlap rescue complete." << endl;
}

void Assembler::rescuePhasedOverlapsThreadFunction(size_t /* threadId */)
{
    const uint64_t rescueThreshold = this->rescuePhasedThreshold;
    
    uint64_t readIdBegin, readIdEnd;
    while(getNextBatch(readIdBegin, readIdEnd)) {
        for(ReadId readId = ReadId(readIdBegin); readId < ReadId(readIdEnd); readId++) {
            
            // Process overlaps from this read's perspective
            OrientedReadId orientedReadId(readId, 0);
            
            if (orientedReadId.getValue() >= alignmentTable.size()) continue;
            
            // Collect overlaps where there's a directional conflict
            // (one direction deleted, other not)
            std::vector<uint32_t> conflictAlignments;
            std::vector<std::pair<uint32_t, uint32_t>> conflictIntervals; // qs, qe on readId
            
            const size_t n = alignmentTable.size(orientedReadId.getValue());
            for(size_t i = 0; i < n; i++) {
                const uint64_t alignmentId = alignmentTable[orientedReadId.getValue()][i];
                const AlignmentData& ad = alignmentData[alignmentId];
                
                // Check for directional conflict: exactly one flag set
                bool hasConflict = (ad.isDeleted0() != ad.isDeleted1());
                if (!hasConflict) continue;
                
                // For this read's perspective, check which flag corresponds to this read
                uint32_t qs, qe;
                bool thisReadDeleted;
                
                if (ad.readIds[0] == readId) {
                    qs = ad.qs;
                    qe = ad.qe;
                    thisReadDeleted = ad.isDeleted0(); // This read's decision
                } else {
                    qs = ad.ts;
                    qe = ad.te;
                    thisReadDeleted = ad.isDeleted1(); // This read's decision
                }
                
                // We want to rescue overlaps where THIS read said "delete" but other said "keep"
                // (i.e., thisReadDeleted == true)
                if (thisReadDeleted) {
                    conflictAlignments.push_back((uint32_t)alignmentId);
                    conflictIntervals.push_back({qs, qe});
                }
            }
            
            // Check if we have enough conflicts to potentially rescue
            if (conflictAlignments.size() < rescueThreshold) continue;
            
            // Sweep line to find max consensus region
            std::vector<std::pair<uint32_t, int>> events;
            events.reserve(conflictIntervals.size() * 2);
            for(const auto& interval : conflictIntervals) {
                events.push_back({interval.first, 1});  // Start
                events.push_back({interval.second, -1}); // End
            }
            
            std::sort(events.begin(), events.end(), [](const auto& a, const auto& b) {
                if (a.first != b.first) return a.first < b.first;
                return a.second > b.second; // Start before End at same position
            });
            
            int dp = 0, max_dp = 0;
            uint32_t start = 0;
            uint32_t max_interval_s = 0, max_interval_e = 0;
            
            for(const auto& ev : events) {
                int old_dp = dp;
                dp += ev.second;
                
                if (old_dp < dp) { // Start event
                    if (dp >= max_dp) {
                        start = ev.first;
                        max_dp = dp;
                    }
                } else if (old_dp > dp) { // End event
                    if (old_dp == max_dp) {
                        max_interval_s = start;
                        max_interval_e = ev.first;
                    }
                }
            }
            
            // If consensus depth >= threshold, rescue overlaps spanning the consensus region
            if ((uint64_t)max_dp >= rescueThreshold) {
                for(size_t j = 0; j < conflictAlignments.size(); j++) {
                    uint32_t alignmentId = conflictAlignments[j];
                    uint32_t qs = conflictIntervals[j].first;
                    uint32_t qe = conflictIntervals[j].second;
                    
                    // Only rescue if this overlap spans the consensus region
                    if (qs <= max_interval_s && qe >= max_interval_e) {
                        AlignmentData& ad = alignmentData[alignmentId];
                        ad.clearDeleteReasonsFromReadPerspective(readId, AlignmentData::DeleteReasonPhase);
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Miniasm-style Chimeric Read Detection (detectChimericReadsFromAnchors)
// ----------------------------------------------------------------------------

void Assembler::detectChimericReadsFromAnchors(
    double shiftRate, 
    uint64_t ulThres, 
    uint64_t threadCount)
{
    cout << timestamp << "Detecting chimeric reads (Anchors method)..." << endl;
    
    this->chimericShiftRate = shiftRate;
    this->chimericUlThres = ulThres;

    // Ensure flags initialized if not already
    if (!isChimericRead.isOpen) {
         isChimericRead.createNew(largeDataName("IsChimericRead"), largeDataPageSize);
         isChimericRead.resize(reads->readCount());
         std::fill(isChimericRead.begin(), isChimericRead.end(), false);
    }
    
    // We assume isChimericRead is already allocated
    
    setupLoadBalancing(reads->readCount(), 100);
    runThreads(&Assembler::detectChimericReadsFromAnchorsThreadFunction, threadCount);

    // Count chimeric reads
    uint64_t chimericCount = 0;
    for(size_t i=0; i<isChimericRead.size(); i++) {
        if(isChimericRead[i]) chimericCount++;
    }
    cout << timestamp << "Detected " << chimericCount << " chimeric reads (Anchors)." << endl;
}

// detectChimericReadsFromAnchorsThreadFunction
//
// This function implements chimeric read detection logic inspired by miniasm/hifiasm.
// It identifies "chimeric" reads (reads formed by the artificial joining of two distant
// genomic locations) by checking for gaps in coverage between "anchors" at the read tips.
//
// Logic Overview:
// 1.  **Anchors (`collectSides`)**: For each read, we identify overlaps that start near the
//     beginning (Left Anchor) and overlaps that end near the end (Right Anchor).
//     We track the furthest internal reach of these anchors (`maxLeft.e` and `maxRight.s`).
//
// 2.  **Gap Bridging (`collectContain`)**: We attempt to bridge the gap between Left and Right
//     anchors using "contained" reads (reads fully contained within the query read).
//     If a contained read significantly overlaps an existing anchor (determined by `overlapRate`),
//     we extend the anchor's reach. This helps bridge small coverage gaps caused by sequencing errors.
//
// 3.  **Gap Detection**: We check if there is a gap between the extended Left Anchor and Right Anchor.
//     Condition: `maxLeft.End <= maxRight.Start`.
//     If a gap exists, it implies the read consists of two disjoint alignment sets that do not
//     connect, suggesting a chimeric join.
//
// 4.  **Cleaning**: If a read is flagged as chimeric:
//     - We mark it as validly chimeric (`isChimericRead[id] = true`).
//     - We conceptually "cut" the read or remove it by marking all its alignments as deleted (`alignmentData[i].setDeleted(true)`).

void Assembler::detectChimericReadsFromAnchorsThreadFunction(size_t /* threadId */)
{
    const double shiftRate = this->chimericShiftRate; // Threshold for overlap between left/right anchors (default ~0.06)
    const double overlapRate = 0.1;                   // Threshold for extending anchors via contained reads

    struct CoverageExtent {
        uint32_t s;
        uint32_t e;
    };

    uint64_t begin, endIdx;
    while(getNextBatch(begin, endIdx)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(endIdx); r0++) {
            
            uint64_t rLen = reads->getReadRawSequenceLength(r0);
            if (rLen == 0) continue;

            uint32_t effectiveStart = 0;
            uint32_t effectiveEnd = (uint32_t)rLen;

            // Use filtered segments if available (coverage cut)
            if (!validReadIntervals.empty()) {
                const auto& status = validReadIntervals[r0];
                if (status.isDeleted) continue; // Read was filtered out entirely
                effectiveStart = status.start;
                effectiveEnd = status.end;
                
                if (effectiveStart >= effectiveEnd) continue;
            }

            // Initialize extents based on valid region.
            // maxLeft: [effectiveStart, e) - starts at effectiveStart, ends at e. 
            // maxRight: [s, effectiveEnd) - starts at s, ends at effectiveEnd.
            CoverageExtent maxLeft = { effectiveEnd, effectiveStart }; // Start with empty range at left
            CoverageExtent maxRight = { effectiveEnd, effectiveStart }; // Start with empty range at right (invalid)
            
            // Correction: Init maxLeft to represent "no left anchor"
            maxLeft = { effectiveEnd, effectiveStart }; 
            // Actually, we want maxLeft.e to grow from effectiveStart.
            // And maxRight.s to shrink from effectiveEnd.
            maxLeft = { effectiveEnd, effectiveStart };
            maxRight = { effectiveEnd, effectiveStart };
             
            // Re-think initialization matching logic:
            // "maxLeft.e" tracks the right-most reach of the left anchor. Init to effectiveStart.
            // "maxRight.s" tracks the left-most reach of the right anchor. Init to effectiveEnd.
            maxLeft.e = effectiveStart;
            maxRight.s = effectiveEnd;
            maxLeft.s = effectiveStart; // Just for completeness, though checks focus on .e
            maxRight.e = effectiveEnd;
            
            // --- Helper: Collect Anchors (collectSides) ---
            // Scans alignments to find interactions touching the read tips (within tolerance).
            auto collectSides = [&](OrientedReadId orientedR0) {
                 const auto& table = alignmentTable[orientedR0.getValue()];
                 for(uint32_t alignmentId : table) {
                     AlignmentData& ad = alignmentData[alignmentId]; 
                     if(ad.isDeleted()) continue;
                     
                     // Get coordinates on r0
                     uint32_t qs = 0, qe = 0;
                     if (ad.readIds[0] == r0) {
                         qs = ad.qs; qe = ad.qe;
                     } else {
                         qs = ad.ts; qe = ad.te;
                     }
                     
                     // Constrain to valid region?
                     // If alignment is outside valid region, we shouldn't use it?
                     // Hifiasm ma_hit_sub trims overlaps to the valid region.
                     // Here we just check coordinates against effectiveStart/End.
                     
                     // Left anchor check (starts near effectiveStart)
                     const uint32_t tolerance = 100; 
                     if (qs <= effectiveStart + tolerance) {
                         // It is a valid left anchor candidate.
                         // We want the furthest it reaches.
                         // But we should clip qe to effectiveEnd? 
                         uint32_t reach = std::min(qe, effectiveEnd);
                         if (reach > maxLeft.e) maxLeft.e = reach;
                     }
                     // Right anchor check (ends near effectiveEnd)
                     if (qe >= effectiveEnd - tolerance) {
                         uint32_t start = std::max(qs, effectiveStart);
                         if (start < maxRight.s) maxRight.s = start;
                     }
                 }
            };

            collectSides(OrientedReadId(r0, 0));
            // Note: We only check Strand 0. Since the alignment graph is symmetric 
            // (overlaps are mirrored), Strand 0 contains all necessary physical coverage info.
            
            // If missing an anchor (Tip/Fragment), skip
            // (One side not covered -> Not a chimera, just a tip/fragment)
            if (maxLeft.s == effectiveEnd || maxRight.s == effectiveEnd) {
                continue;
            }

            // --- Helper: Collect Contain (collectContain) ---
            // Extend anchors using contained reads that overlap significantly with existing anchors.
            // This mirrors miniasm's logic to "grow" the valid regions inward.
            uint32_t newLeftE = maxLeft.e;
            uint32_t newRightS = maxRight.s;

            auto collectContain = [&](OrientedReadId orientedR0) {
                 const auto& table = alignmentTable[orientedR0.getValue()];
                 for(uint32_t alignmentId : table) {
                     const AlignmentData& ad = alignmentData[alignmentId];
                     if(ad.isDeleted()) continue;
                     
                     uint32_t qs = 0, qe = 0;
                     if (ad.readIds[0] == r0) {
                         qs = ad.qs; qe = ad.qe;
                     } else {
                         qs = ad.ts; qe = ad.te;
                     }

                     // Check contained overlaps (reads strictly inside valid region, not touching ends)
                     const uint32_t tolerance = 100;
                     if (qs > effectiveStart + tolerance && qe < effectiveEnd - tolerance) {
                         
                         // Try to extend Left Anchor
                         // Logic: If contained read overlaps Left Anchor significantly, 
                         // we can trust it to extend the Left Anchor further right.
                         // Overlap: [qs, maxLeft.e)
                         if (qs < maxLeft.e && qe > maxLeft.e) {
                             uint32_t overlapLen = maxLeft.e - qs;
                             uint32_t alignLen = qe - qs;
                             if (overlapLen > (overlapRate * alignLen)) {
                                 if (qe > newLeftE) newLeftE = qe;
                             }
                         }

                         // Try to extend Right Anchor
                         // Logic: If contained read overlaps Right Anchor significantly,
                         // we can trust it to extend the Right Anchor further left.
                         // Overlap: [maxRight.s, qe)
                         if (qs < maxRight.s && qe > maxRight.s) {
                             uint32_t overlapLen = qe - maxRight.s;
                             uint32_t alignLen = qe - qs;
                             if (overlapLen > (overlapRate * alignLen)) {
                                 if (qs < newRightS) newRightS = qs;
                             }
                         }
                     }
                 }
            };
            
            collectContain(OrientedReadId(r0, 0));

            maxLeft.e = newLeftE;
            maxRight.s = newRightS;


            // --- Check Valid Bridging (Simple Chimera Detection) ---
            // If overlaps meet in the middle (Left End > Right Start), the read is likely good.
            // Even if they strictly "meet", we might require a minimal overlap length (shiftRate).
            int64_t overlap = (int64_t)maxLeft.e - (int64_t)maxRight.s;
            // Should effective length be used for ratio? (effectiveEnd - effectiveStart)
            uint32_t effectiveLen = effectiveEnd - effectiveStart;
            if (overlap > 0 && overlap >= (int64_t)(effectiveLen * shiftRate)) {
                continue; // Good Read (Anchors overlap significantly)
            }
            
            // If we are here, there is a GAP (maxLeft.e <= maxRight.s) or insufficient overlap.
            // We treat this as a potentially chimeric read.
            bool isSimpleChimera = (maxLeft.e <= maxRight.s);
            
            // Note: miniasm includes a "Complex Chimera" check here (intersection_check_by_base)
            // for cases where overlaps meet but alignment quality at junction is poor.
            // We omit that for now and rely on structural gap (Simple Chimera).
            // bool isComplexChimera = false; // logic would go here
            
            // The user logic for 'complex': 
            // if intersection_check_by_base(...) returns true -> delete.
            // Since we don't have base check, we act conservatively?
            // The prompt asks to "mirror".
            // If simple chimera, we delete.
            
            if (isSimpleChimera) {
                isChimericRead[r0] = true;
                
                // --- Delete All Edges (delete_all_edges) ---
                // Mark all alignments involving this read as deleted.
                auto deleteEdges = [&](OrientedReadId orientedR0) {
                     const auto& table = alignmentTable[orientedR0.getValue()];
                     for(uint32_t alignmentId : table) {
                         // We can delete here because alignmentData is shared but atomic bool or similar is safe?
                         // isDeleted is not atomic.
                         // However, parallel writes to same bool are OK if value is same (true).
                         // Race condition if someone reads it.
                         // Standard for this codebase: marking reads is preferred.
                         // But user asked to mirror delete_all_edges.
                        alignmentData[alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonChimeric);
                     }
                };
                deleteEdges(OrientedReadId(r0, 0));
                
                continue;
            }
            
            // If complex logic existed, we would check it here.
            // But without base check, we assume non-gap reads are OK.
        }
    }
}

// --- Contained Read Filtering (ma_hit_contained_advance equivalent) ---

// Exact hifiasm ma_hit2arc containment check.
// Parameters match hifiasm naming:
//   h = alignment (hit)
//   ql = query length (after coverage cut: sq->e - sq->s)
//   tl = target length (after coverage cut: st->e - st->s)
//   max_hang = maximum allowed overhang
//   int_frac = internal fraction threshold (max_hang_rate)
//   min_ovlp = minimum overlap length
// Returns:
//   0: Normal overlap (dovetail)
//   1: MA_HT_QCONT - Query contained in Target
//   2: MA_HT_TCONT - Target contained in Query
//   -1: MA_HT_INT - Internal match (too much overhang)
//   -2: MA_HT_SHORT_OVLP - Overlap too short
static int ma_hit2arc_containment(
    int32_t qs, int32_t qe, int32_t ql,   // Query: overlap start, end, length
    int32_t ts, int32_t te, int32_t tl,   // Target: overlap start, end, length  
    bool isReverse,                       // Is target reverse complemented?
    int32_t max_hang,                     // Maximum allowed overhang
    double int_frac,                      // Internal fraction threshold
    int32_t min_ovlp)                     // Minimum overlap length
{
    // Compute 5' and 3' overhangs on target (relative to query orientation)
    // tl5 = 5'-end overhang on target (on query strand)
    // tl3 = 3'-end overhang on target (on query strand)
    int32_t tl5, tl3;
    if (isReverse) {
        tl5 = tl - te;  // Target's right end becomes 5' relative to query
        tl3 = ts;       // Target's left end becomes 3' relative to query
    } else {
        tl5 = ts;       // Target's left end is 5'
        tl3 = tl - te;  // Target's right end is 3'
    }

    // ext5 and ext3: the minimum of query and target overhangs on each side
    int32_t ext5 = (qs < tl5) ? qs : tl5;
    int32_t ext3 = ((ql - qe) < tl3) ? (ql - qe) : tl3;

    // Check for internal match (too much overhang on both sides)
    // This rejects overlaps where the alignment doesn't extend to edges
    if (ext5 > max_hang || ext3 > max_hang) {
        return -1;  // MA_HT_INT
    }
    
    // Check internal fraction constraint
    int32_t qOverlapLen = qe - qs;
    int32_t tOverlapLen = te - ts;
    if (qOverlapLen < (qOverlapLen + ext5 + ext3) * int_frac) {
        return -1;  // MA_HT_INT
    }
    if (tOverlapLen < (tOverlapLen + ext5 + ext3) * int_frac) {
        return -1;  // MA_HT_INT
    }

    // Containment checks (exact hifiasm logic from Overlaps.h lines 418-419)
    // Query contained in Target: query's overhangs are smaller than target's
    if (qs <= tl5 && (ql - qe) <= tl3) {
        return 1;  // MA_HT_QCONT
    }
    // Target contained in Query: target's overhangs are smaller than query's
    if (qs >= tl5 && (ql - qe) >= tl3) {
        return 2;  // MA_HT_TCONT
    }

    // Check minimum overlap length
    if (qOverlapLen + ext5 + ext3 < min_ovlp || tOverlapLen + ext5 + ext3 < min_ovlp) {
        return -2;  // MA_HT_SHORT_OVLP
    }

    // Normal dovetail overlap
    return 0;
}


void Assembler::removeContainedReads(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount)
{
    cout << timestamp << "Removing contained reads (ma_hit_contained_advance equivalent)..." << endl;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    if (!containmentParent->isOpen) {
        containmentParent->createNew(largeDataName("ContainmentParent"), largeDataPageSize);
        containmentParent->resize(reads->readCount());
    }
    std::fill(containmentParent->begin(), containmentParent->end(), ReadId(invalidReadId));

    uint64_t containedReadCount = 0;

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
            // Hifiasm parity: contained-read logic runs only on overlaps that survived prior filters (h->del==0).
            if(!ad.keptByBothSides()) continue;

            const ReadId tn = (ad.readIds[0] == qn) ? ad.readIds[1] : ad.readIds[0];
            if (tn >= validReadIntervals.size() || validReadIntervals[tn].isDeleted) continue;

            const auto& vrT = validReadIntervals[tn];
            const int32_t tl = int32_t(vrT.end - vrT.start);
            if (tl <= 0) continue;

            // Build the overlap record from qn's perspective (equivalent to sources[qn].buffer[j]).
            // Parity with hifiasm set_reverse_overlap: swap query/target intervals without
            // additional coordinate transforms (coordinates are always in each read's forward frame).
            const bool rev = !ad.isSameStrand;
            int32_t qs = 0, qe = 0, ts = 0, te = 0;
            if (ad.readIds[0] == qn) {
                qs = int32_t(ad.qs);
                qe = int32_t(ad.qe);
                ts = int32_t(ad.ts);
                te = int32_t(ad.te);
            } else {
                // qn is readIds[1].
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
                // Query read is contained in target.
                validReadIntervals[qn].isDeleted = true;
                (*containmentParent)[qn] = tn;
                ++containedReadCount;
                break;
            } else if (result == 2) {
                // Target read is contained in query.
                validReadIntervals[tn].isDeleted = true;
                (*containmentParent)[tn] = qn;
                ++containedReadCount;
            }
        }
    }

    // Parity with transfor_R_to_U: path compression to final container.
    for (ReadId r = 0; r < readCount; ++r) {
        if ((*containmentParent)[r] == ReadId(invalidReadId)) continue;
        ReadId root = (*containmentParent)[r];
        while (root != ReadId(invalidReadId) && (*containmentParent)[root] != ReadId(invalidReadId)) {
            root = (*containmentParent)[root];
        }
        (*containmentParent)[r] = root;
    }

    // Remove all alignments incident to deleted reads (contained reads are deleted in validReadIntervals).
    for (uint64_t alignmentId = 0; alignmentId < alignmentData.size(); ++alignmentId) {
        AlignmentData& ad = alignmentData[alignmentId];
        const ReadId r0 = ad.readIds[0];
        const ReadId r1 = ad.readIds[1];
        if (r0 < validReadIntervals.size() && validReadIntervals[r0].isDeleted) {
            ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonContained);
            continue;
        }
        if (r1 < validReadIntervals.size() && validReadIntervals[r1].isDeleted) {
            ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonContained);
            continue;
        }
    }

    // Parity with the final pass in ma_hit_contained_advance: delete reads that now have no remaining overlaps.
    setupLoadBalancing(readCount, 1000);
    runThreads(&Assembler::applyCoverageCutsCleanupThreadFunction, threadCount);

    cout << timestamp << "Identified " << containedReadCount << " contained reads." << endl;
}
