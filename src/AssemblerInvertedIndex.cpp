/**
 * @file AssemblerInvertedIndex.cpp
 * @brief High-performance alignment candidate discovery using an Inverted Index.
 *
 * This file implements a Hifiasm-compatible chaining algorithm for finding
 * potential read overlaps. The core data flow is:
 *
 *   1. **Inverted Index Construction**: Build a hash table mapping each
 *      canonical K-mer to a list of (ReadId, Position) occurrences.
 *   2. **Parallel Radix Sort**: Sort occurrences by K-mer ID using an O(N)
 *      LSD radix sort to group identical K-mers.
 *   3. **Hash Table Build**: Populate a power-of-2 sized hash table for O(1)
 *      K-mer lookups during the search phase.
 *   4. **Parallel Candidate Search**: For each read, query the index to find
 *      matching K-mers in other reads. Use DP chaining (with optional AVX2
 *      SIMD pre-filtering) to score potential alignments.
 *
 * Key performance features:
 *   - Structure-of-Arrays (SoA) layout for cache-efficient DP.
 *   - Early K-mer weighting based on occurrence frequency (Hifiasm-compatible).
 *   - AVX2 SIMD pre-filter to skip non-viable DP predecessors (optional).
 *   - Collinear fast-path to skip O(N^2) DP for trivially monotonic hit sets.
 *
 * @note All scoring and tie-breaking rules are strictly Hifiasm-compatible.
 */

#include "Assembler.hpp"
#include "hifiasmCoordinateTransforms.hpp"
#include "performanceLog.hpp"
#include "OrientedReadPair.hpp"
#include "timestamp.hpp"
#include "Reads.hpp"
#include <algorithm>
#include <map>
#include <mutex>
#include <cmath>
#include <vector>
#include <thread>
#include <functional>

// AVX2 intrinsics are only available on x86_64 architectures.
// The code includes a scalar fallback for other platforms.
#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

using namespace std;

#include "MultithreadedObject.tpp"
#include "Alignment.hpp"
#include "AlignmentCanonicalization.hpp"

namespace dinara {

/**
 * @brief Temporary hit structure used during K-mer matching.
 *
 * Represents a single K-mer match between the query read (Read A) and a
 * partner read. Hits are first collected in Array-of-Structures (AoS) format
 * for sorting, then converted to Structure-of-Arrays (SoA) for DP.
 *
 * @note Sorting is by (partnerReadId, posA) to group hits by read pair and
 *       ensure monotonically increasing query positions for the DP.
 */
struct InvertedIndexTempHit {
    ReadId partnerReadId;  ///< ID of the matching read (Read B).
    uint32_t posA;         ///< Base position of the K-mer in Read A.
    uint32_t posB;         ///< Base position of the K-mer in Read B.
    uint32_t ordinalA;     ///< Marker ordinal in Read A (for alignment output).
    uint8_t weight;        ///< Frequency-based weight (Hifiasm compatible).
    
    /// Comparison operator for sorting by (partnerReadId, posA).
    bool operator<(const InvertedIndexTempHit& other) const {
        if (partnerReadId != other.partnerReadId) return partnerReadId < other.partnerReadId;
        return posA < other.posA;
    }
};

/**
 * @brief Classifies an overlap based on Hifiasm's overlap type heuristic.
 *
 * This is used for Hifiasm's candidate limit feature, which caps the number
 * of reported overlaps per type to avoid excessive output in repetitive regions.
 *
 * @param start  Starting coordinate of the alignment on Read A.
 * @param end    Ending coordinate of the alignment on Read A.
 * @param readLen Total length of Read A.
 * @return 0 = Left overhang, 1 = Right overhang, 2 = Contained, 3 = Containing.
 */
static int getOverlapType(uint32_t start, uint32_t end, uint32_t readLen) {
    if (start == 0 && end >= readLen - 1) return 2; // Contained
    else if (start > 0 && end < readLen - 1) return 3; // Containing
    else return (start == 0) ? 0 : 1; // Left or Right overhang
}

/**
 * @brief Per-thread scratchpad for high-performance DP chaining.
 *
 * This struct uses a Structure-of-Arrays (SoA) layout for the DP arrays,
 * which is critical for cache efficiency during the O(N^2) chaining loop.
 * The scratchpad is recycled across reads to avoid allocation overhead.
 *
 * Key arrays:
 *   - `dpSame`, `dpDiff`: DP scores for same-strand and opposite-strand chains.
 *   - `parentSame`, `parentDiff`: Backtrack pointers for chain reconstruction.
 *   - `cumDriftSame`, `cumDriftDiff`: Cumulative indel count for drift constraint.
 *   - `cumLenSame`, `cumLenDiff`: Cumulative alignment length for gap rate.
 */
struct ThreadScratchpad {
    // AoS for hit collection and sorting
    vector<InvertedIndexTempHit> flatHits;
    
    // Structure of Arrays (SoA) for cache-efficient DP scans
    vector<uint32_t> hitPosA, hitPosB, hitOrdinalA;
    vector<uint8_t> hitWeights;

    // DP score and backtrack arrays
    vector<uint32_t> dpSame, dpDiff;
    vector<int32_t> parentSame, parentDiff;
    vector<uint32_t> cumDriftSame, cumDriftDiff;
    vector<uint32_t> cumLenSame, cumLenDiff;
    vector<int32_t> backtrackVisitSame, backtrackVisitDiff;
    vector<uint32_t> chainOccurrencesSame, chainOccurrencesDiff;

    // Post-DP candidate extraction
    struct ChainCandidate {
        uint32_t score;
        uint64_t chainLen; // For deterministic tie-breaking
        int32_t endK;
        bool isDiff; 
        bool operator<(const ChainCandidate& other) const {
            if (score != other.score) return score > other.score;
            return chainLen < other.chainLen;
        }
    };
    vector<ChainCandidate> chainCandidates;
    vector<int> candidateTypes;
    vector<ChainCandidate> filteredCandidates;

    struct ChainInterval { uint32_t qs; uint32_t qe; };
    vector<ChainInterval> acceptedIntervalsSame;
    vector<ChainInterval> acceptedIntervalsDiff;
    vector<uint32_t> currentChainPath;

    void clear() {
        flatHits.clear();
        hitPosA.clear(); hitPosB.clear(); hitOrdinalA.clear(); hitWeights.clear();
        dpSame.clear(); dpDiff.clear();
        parentSame.clear(); parentDiff.clear();
        cumDriftSame.clear(); cumDriftDiff.clear();
        cumLenSame.clear(); cumLenDiff.clear();
        backtrackVisitSame.clear(); backtrackVisitDiff.clear();
        chainOccurrencesSame.clear(); chainOccurrencesDiff.clear();
        chainCandidates.clear();
        candidateTypes.clear();
        filteredCandidates.clear();
        acceptedIntervalsSame.clear();
        acceptedIntervalsDiff.clear();
        currentChainPath.clear();
    }
};

/**
 * @brief Fast reverse complement for K-mers using bit reversal.
 */
inline KmerId getRcKmerId(KmerId id, uint64_t k) {
    const KmerId mask = (KmerId(1) << k) - 1;
    const KmerId lsb = id & mask;
    const KmerId msb = (id >> k) & mask;
    
    auto reverseBits = [&](KmerId x) -> KmerId {
        if (sizeof(KmerId) <= 8) return bitReversal(uint64_t(x)) >> (64 - k);
        else return bitReversal((__uint128_t)x) >> (128 - k);
    };
    
    KmerId rc_lsb = (~reverseBits(lsb)) & mask;
    KmerId rc_msb = (~reverseBits(msb)) & mask;
    return (rc_msb << k) | rc_lsb;
}

/**
 * @brief Thread-safe hash function for K-mers.
 */
static inline uint64_t hashKmer(KmerId k) {
    const uint64_t* p = reinterpret_cast<const uint64_t*>(&k);
    uint64_t k1 = p[0];
    uint64_t k2 = sizeof(KmerId) > 8 ? p[1] : 0; 
    return k1 ^ (k2 + 0x9e3779b9 + (k1<<6) + (k1>>2));
}

// Private class to encapsulate parallel logic (Codebase Pattern).
class InvertedIndexFinder : public MultithreadedObject<InvertedIndexFinder> {
public:
    InvertedIndexFinder(
        const Reads& reads,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MemoryMapped::VectorOfVectors<KmerId, uint64_t>& markerKmerIds,
        const Assembler::AlignmentCandidatesInvertedIndexData& invertedIndexData,
        MemoryMapped::Vector<OrientedReadPair>& candidates,
        MemoryMapped::Vector<Alignment>& precomputedAlignments,
        uint64_t maxChainLimit,
        uint32_t minChainedMarkerCount,
        uint64_t threadCount
    ) : 
        MultithreadedObject(*this),
        reads(reads),
        markers(markers),
        markerKmerIds(markerKmerIds),
        invertedIndexData(invertedIndexData),
        candidates(candidates),
        precomputedAlignments(precomputedAlignments),
        maxChainLimit(maxChainLimit),
        minChainedMarkerCount(minChainedMarkerCount),
        threadCount(threadCount)
    {
        // 1. Setup per-thread accumulation buffers and scratchpads
        threadCandidates.resize(threadCount);
        for(auto& v : threadCandidates) v.reserve(10000);
        threadAlignments.resize(threadCount);
        for(auto& v : threadAlignments) v.reserve(10000);
        threadScratchpads.resize(threadCount);

        // 2. Setup parallel workload partitioning
        const ReadId readCount = ReadId(markers.size() / 2); // Indexed by strand 0
        setupLoadBalancing(readCount, 100);

        // 3. Kick off parallel worker threads
        runThreads(&InvertedIndexFinder::threadFunction, threadCount);

        // 4. Consolidate results: Resize global vectors once and copy results back
        size_t totalCandidatesFound = 0;
        for(const auto& v : threadCandidates) totalCandidatesFound += v.size();
        
        cout << "Discovery search complete. Merging " << totalCandidatesFound << " candidates." << endl;

        size_t candidateWritePos = candidates.size();
        candidates.resize(candidateWritePos + totalCandidatesFound);
        size_t alignmentWritePos = precomputedAlignments.size();
        precomputedAlignments.resize(alignmentWritePos + totalCandidatesFound);
        
        for(size_t i = 0; i < threadCount; i++) {
            const auto& v = threadCandidates[i];
            const auto& a = threadAlignments[i];
            if(!v.empty()) {
                std::copy(v.begin(), v.end(), candidates.begin() + candidateWritePos);
                candidateWritePos += v.size();
                std::copy(a.begin(), a.end(), precomputedAlignments.begin() + alignmentWritePos);
                alignmentWritePos += a.size();
            }
        }
    }

private:
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MemoryMapped::VectorOfVectors<KmerId, uint64_t>& markerKmerIds;
    const Assembler::AlignmentCandidatesInvertedIndexData& invertedIndexData;
    [[maybe_unused]] MemoryMapped::Vector<OrientedReadPair>& candidates;
    [[maybe_unused]] MemoryMapped::Vector<Alignment>& precomputedAlignments;
    uint64_t maxChainLimit;
    uint32_t minChainedMarkerCount;
    [[maybe_unused]] uint64_t threadCount;

    vector<vector<OrientedReadPair>> threadCandidates;
    vector<vector<Alignment>> threadAlignments;
    vector<ThreadScratchpad> threadScratchpads;

    void threadFunction(size_t threadId) {
        vector<OrientedReadPair>& localCandidates = threadCandidates[threadId];
        vector<Alignment>& localAlignments = threadAlignments[threadId];
        ThreadScratchpad& scratch = threadScratchpads[threadId];
        
        const uint64_t hashMask = invertedIndexData.hashTable.size() - 1;
        const auto* hashTablePtr = invertedIndexData.hashTable.data();
        const uint64_t kmerLen = invertedIndexData.k;
        const double maxDriftRate = invertedIndexData.maxDriftRate;
        const uint64_t coveragePeak = invertedIndexData.coveragePeak;
        const uint8_t* weightLUT = invertedIndexData.weightLut.data();
        const double weightExponent = invertedIndexData.weightExponent;
        const double lowFreqMultiplier = invertedIndexData.lowFreqMultiplier;
        const double highFreqMultiplier = invertedIndexData.highFreqMultiplier;
        const uint32_t rareKmerWeight = invertedIndexData.rareKmerWeight;
        const double chainFilterRatio = invertedIndexData.chainFilterRatio;
        const uint32_t chainFilterMinScore = invertedIndexData.chainFilterMinScore;
        const double nonRedundantOverlapFraction = invertedIndexData.nonRedundantOverlapFraction;


        uint64_t startBatch, endBatch;
        while(getNextBatch(startBatch, endBatch)) {
            for(ReadId readIdA = ReadId(startBatch); readIdA != ReadId(endBatch); ++readIdA) {
                
                const OrientedReadId orientedReadIdA(readIdA, 0);
                const auto& markersA = markers[orientedReadIdA.getValue()];
                const auto& kmerIdsA = markerKmerIds[orientedReadIdA.getValue()];
                const size_t numMarkersA = std::min(markersA.size(), kmerIdsA.size());
                
                scratch.clear();
                scratch.flatHits.reserve(numMarkersA * 2);
                
                // --- Step 1: Hit Collection & Early Weighting ---
                // We scan markers in Read A and find matches in the Inverted Index.
                for(size_t i = 0; i < numMarkersA; ++i) {
                    KmerId currentKId = kmerIdsA[i];
                    KmerId rcKId = getRcKmerId(currentKId, kmerLen);
                    KmerId canonicalKId = (currentKId < rcKId) ? currentKId : rcKId;

                    uint64_t slotIdx = hashKmer(canonicalKId) & hashMask;
                    const uint32_t posA = markersA[i].position;

                    // Search for the K-mer in the direct-addressing hash table
                    while(!hashTablePtr[slotIdx].empty) {
                        if(hashTablePtr[slotIdx].key == canonicalKId) {
                            const uint64_t startIdx = hashTablePtr[slotIdx].start;
                            const uint32_t count = hashTablePtr[slotIdx].count;
                             
                            // Early Weighting: Constant weight for all hits of this k-mer.
                            // Frequent/Repetitive k-mers are penalized using a pre-computed LUT.
                            // The LUT is generated once at startup using pow(w, exponent) for Hifiasm compatibility.
                            uint8_t hitWeight = 1;
                            const uint64_t lowFreq = (uint64_t)(double(coveragePeak) * lowFreqMultiplier);
                            const uint64_t highFreq = (uint64_t)(double(coveragePeak) * highFreqMultiplier);
                            
                            if (count <= std::max(2UL, lowFreq)) {
                                hitWeight = uint8_t(std::min<uint32_t>(255U, rareKmerWeight)); // Rare/Informative k-mer (High value)
                            } else if (count >= std::max(3UL, highFreq)) {
                                uint32_t w = 1 + (uint32_t)((count + (highFreq * 2) - 1) / (highFreq * 2 == 0 ? 1 : highFreq * 2));
                                if(w < 512) {
                                    hitWeight = weightLUT[w];
                                } else {
                                    hitWeight = (uint8_t)std::min(255U, (uint32_t)pow((double)w, weightExponent));
                                }
                            }

                            const auto* compactOccs = &invertedIndexData.compactOccurrences[startIdx];
                            for(uint32_t j = 0; j < count; ++j) {
                                if(compactOccs[j].readId != readIdA) {
                                    scratch.flatHits.push_back({compactOccs[j].readId, posA, compactOccs[j].position, (uint32_t)i, hitWeight}); 
                                }
                            }
                            break;
                        }
                        slotIdx = (slotIdx + 1) & hashMask;
                    }
                }
                
                if(scratch.flatHits.empty()) continue;
                std::sort(scratch.flatHits.begin(), scratch.flatHits.end());

                // --- Step 2: DP Chaining per Read Pair ---
                const uint64_t readLenA = reads.getReadRawSequenceLength(readIdA);
                size_t hitIter = 0;
                while(hitIter < scratch.flatHits.size()) {
                    const ReadId readIdB = scratch.flatHits[hitIter].partnerReadId;
                    
                    // Skip mirrored pairs and self-comparisons
                    if(readIdB <= readIdA) {
                        while(hitIter < scratch.flatHits.size() && scratch.flatHits[hitIter].partnerReadId == readIdB) hitIter++;
                        continue;
                    }

                    const size_t startInFlat = hitIter;
                    while(hitIter < scratch.flatHits.size() && scratch.flatHits[hitIter].partnerReadId == readIdB) hitIter++;
                    const size_t numHits = hitIter - startInFlat;
                    if(numHits == 0) continue;

                    const uint64_t readLenB = reads.getReadRawSequenceLength(readIdB);

                    // 2.1 Struct-of-Arrays (SoA) Transfer for cache-efficient DP access
                    scratch.hitPosA.assign(numHits, 0); scratch.hitPosB.assign(numHits, 0);
                    scratch.hitOrdinalA.assign(numHits, 0); scratch.hitWeights.assign(numHits, 0);
                    for(size_t k = 0; k < numHits; ++k) {
                        const auto& h = scratch.flatHits[startInFlat + k];
                        scratch.hitPosA[k] = h.posA; scratch.hitPosB[k] = h.posB;
                        scratch.hitOrdinalA[k] = h.ordinalA; scratch.hitWeights[k] = h.weight;
                    }

                    // Pre-allocate/Reset DP work arrays
                    scratch.dpSame.assign(numHits, 0); scratch.dpDiff.assign(numHits, 0);
                    scratch.parentSame.assign(numHits, -1); scratch.parentDiff.assign(numHits, -1);
                    scratch.cumDriftSame.assign(numHits, 0); scratch.cumDriftDiff.assign(numHits, 0);
                    scratch.cumLenSame.assign(numHits, 0); scratch.cumLenDiff.assign(numHits, 0);
                    scratch.backtrackVisitSame.assign(numHits, -1); scratch.backtrackVisitDiff.assign(numHits, -1);
                    scratch.chainOccurrencesSame.assign(numHits, 1); scratch.chainOccurrencesDiff.assign(numHits, 1);

                    uint32_t maxScSame = 0, maxScDiff = 0; 
                    int32_t bestEndIdxSame = -1, bestEndIdxDiff = -1; 
                    uint64_t maxExtSame = UINT64_MAX, maxExtDiff = UINT64_MAX;
                    
                    const int64_t dRscaled = (int64_t)(maxDriftRate * 1024.0);

                    // Utility: Calculate expected coordinate extension (length) for a hit
                    auto getAlignmentLength = [&](uint32_t pA, uint32_t pB) -> uint64_t {
                        uint32_t xB = pA, yB = pB; if (xB <= yB) { yB -= xB; xB = 0; } else { xB -= yB; yB = 0; }
                        uint64_t xR = readLenA - pA - 1, yR = readLenB - pB - 1;
                        uint32_t xE = (xR <= yR) ? (uint32_t)(readLenA - 1) : pA + (uint32_t)yR;
                        return (uint64_t)(xE - xB + 1);
                    };

                    // [DISABLED] Collinear Fast-Path: The O(N) pre-check rarely pays off on
                    // noisy long-read data. Kept here for reference/future benchmarking.
                    #if 0
                    bool isStrictlyCollinear = true; 
                    if (numHits > 1) { 
                        for (size_t k = 1; k < numHits; ++k) {
                            int32_t dx = (int32_t)scratch.hitPosA[k] - (int32_t)scratch.hitPosA[k-1];
                            int32_t dy = (int32_t)scratch.hitPosB[k] - (int32_t)scratch.hitPosB[k-1];
                            if (dx <= 0 || dy <= 0 || std::abs(dy - dx) > (int32_t)(maxDriftRate * dx)) {
                                isStrictlyCollinear = false; break;
                            }
                        } 
                    }

                    if (isStrictlyCollinear && numHits > 0) {
                        uint32_t baseSc = scratch.hitWeights[0] > 1 ? (uint32_t)kmerLen / scratch.hitWeights[0] : (uint32_t)kmerLen;
                        scratch.dpSame[0] = std::max(1U, baseSc); 
                        int32_t driftS = 0, lenS = 0; bool validS = true;
                        for (size_t k = 1; k < numHits; ++k) {
                            int32_t dx = (int32_t)scratch.hitPosA[k] - (int32_t)scratch.hitPosA[k-1];
                            int32_t dy = (int32_t)scratch.hitPosB[k] - (int32_t)scratch.hitPosB[k-1];
                            int32_t dd = std::abs(dx - dy); driftS += dd; lenS += dy;
                            if (dy <= 0 || driftS > lenS * maxDriftRate || (dd > 31 && dd > std::min(dx, dy) * maxDriftRate)) { 
                                validS = false; break; 
                            }
                            uint32_t wS = scratch.hitWeights[k] > 1 ? std::min((uint32_t)std::min(dx, dy), (uint32_t)kmerLen) / scratch.hitWeights[k] : std::min((uint32_t)std::min(dx, dy), (uint32_t)kmerLen);
                            int32_t pnlty = (lenS > 0) ? (int32_t)(((int64_t)driftS * wS * 1024) / ((int64_t)lenS * dRscaled)) : 0;
                            scratch.dpSame[k] = scratch.dpSame[k-1] + std::max(1, (int32_t)wS - pnlty);
                            scratch.parentSame[k] = (int32_t)(k - 1); 
                            scratch.cumDriftSame[k] = driftS; scratch.cumLenSame[k] = lenS; 
                            scratch.chainOccurrencesSame[k] = (uint32_t)(k + 1);
                        }
                        if (validS) { maxScSame = scratch.dpSame[numHits - 1]; bestEndIdxSame = (int32_t)(numHits - 1); goto end_dp_chaining; }
                    }
                    #endif

                    // General DP Chaining with SIMD Pre-filtering
                    for(size_t k = 0; k < numHits; ++k) {
                        const uint32_t posAk = scratch.hitPosA[k], posBk = scratch.hitPosB[k];
                        const uint32_t initSc = std::max(1U, scratch.hitWeights[k] > 1 ? (uint32_t)kmerLen / scratch.hitWeights[k] : (uint32_t)kmerLen);
                        scratch.dpSame[k] = scratch.dpDiff[k] = initSc;
                        
                        size_t streakS = 0, streakD = 0, gapStreakS = 0, gapStreakD = 0; 
                        int32_t j = (int32_t)k - 1;

                        for(; j >= 0; --j) {
// [DISABLED] SIMD Pre-filter: On non-repetitive data, the overhead of loading
// 8 predecessors into AVX2 registers often exceeds the benefit of skipping.
// Enable via -DUSE_SIMD_DP_PREFILTER for benchmarking in repetitive regions.
#if defined(__AVX2__) && defined(USE_SIMD_DP_PREFILTER)
                            // SIMD Block: Quick-scan 8 predecessors for strand and drift compatibility
                            if (j >= 7) {
                                __m256i vPA = _mm256_loadu_si256((const __m256i*)&scratch.hitPosA[j-7]);
                                __m256i vPB = _mm256_loadu_si256((const __m256i*)&scratch.hitPosB[j-7]);
                                __m256i vDA = _mm256_sub_epi32(_mm256_set1_epi32(posAk), vPA);
                                __m256i vDBS = _mm256_sub_epi32(_mm256_set1_epi32(posBk), vPB);
                                __m256i vDBD = _mm256_sub_epi32(vPB, _mm256_set1_epi32(posBk));

                                __m256i vMdA = _mm256_cmpgt_epi32(vDA, _mm256_setzero_si256());
                                __m256i vMS = _mm256_and_si256(vMdA, _mm256_cmpgt_epi32(vDBS, _mm256_setzero_si256()));
                                __m256i vMD = _mm256_and_si256(vMdA, _mm256_cmpgt_epi32(vDBD, _mm256_setzero_si256()));
                                
                                __m256i vL = _mm256_mullo_epi32(vDA, _mm256_set1_epi32((int)dRscaled));
                                // Note: AVX2 has no _mm256_cmplt_epi32, so we use _mm256_cmpgt_epi32 with swapped operands.
                                __m256i vDriftValS = _mm256_slli_epi32(_mm256_abs_epi32(_mm256_sub_epi32(vDBS, vDA)), 10);
                                __m256i vDriftValD = _mm256_slli_epi32(_mm256_abs_epi32(_mm256_sub_epi32(vDBD, vDA)), 10);
                                __m256i vDriftS = _mm256_and_si256(vMS, _mm256_cmpgt_epi32(vL, vDriftValS));
                                __m256i vDriftD = _mm256_and_si256(vMD, _mm256_cmpgt_epi32(vL, vDriftValD));

                                if (_mm256_testz_si256(vDriftS, vDriftS) && _mm256_testz_si256(vDriftD, vDriftD)) {
                                    j -= 7; continue; 
                                }
                            }
#endif
                            int32_t dA = (int32_t)posAk - (int32_t)scratch.hitPosA[j];
                            int32_t dB = (int32_t)posBk - (int32_t)scratch.hitPosB[j];
                            
                            // Case 1: Hits on the same strand
                            if(dB > 0 && dA > 0) { 
                                // Hifiasm's chaining constraint uses the query/self-axis as the length accumulator.
                                // Here that is Read A (dA).
                                uint32_t nDr = scratch.cumDriftSame[j] + std::abs(dB - dA);
                                uint32_t nLn = scratch.cumLenSame[j] + uint32_t(dA);
                                if((int64_t)((uint64_t)nDr << 10) <= (int64_t)dRscaled * (int64_t)nLn) {
                                    uint32_t wS = scratch.hitWeights[j] > 1 ? std::min((uint32_t)std::min(dA, dB), (uint32_t)kmerLen) / scratch.hitWeights[j] : std::min((uint32_t)std::min(dA, dB), (uint32_t)kmerLen);
                                    int32_t pnlty = (nLn > 0) ? (int32_t)(((int64_t)nDr * wS * 1024) / ((int64_t)nLn * dRscaled)) : 0;
                                    int32_t sc = (int32_t)scratch.dpSame[j] + (int32_t)wS - pnlty;
                                    if(sc > (int32_t)scratch.dpSame[k]) { 
                                        scratch.dpSame[k] = std::max(1, sc); scratch.parentSame[k] = j; 
                                        scratch.cumDriftSame[k] = nDr; scratch.cumLenSame[k] = nLn; 
                                        scratch.chainOccurrencesSame[k] = scratch.chainOccurrencesSame[j] + 1; 
                                        streakS = 0; if (gapStreakS > 0) --gapStreakS; 
                                    } else { 
                                        if(++streakS > 25) break; 
                                        if (scratch.backtrackVisitSame[j] == (int32_t)k && ++gapStreakS > 25) break; 
                                    }
                                    if (scratch.parentSame[j] >= 0) scratch.backtrackVisitSame[scratch.parentSame[j]] = (int32_t)k;
                                }
                            } 
                            // Case 2: Hits on opposite strands
                            else if (dB < 0 && dA > 0) { 
                                int32_t aB = -dB; 
                                uint32_t nDr = scratch.cumDriftDiff[j] + std::abs(aB - dA);
                                uint32_t nLn = scratch.cumLenDiff[j] + uint32_t(dA);
                                if((int64_t)((uint64_t)nDr << 10) <= (int64_t)dRscaled * (int64_t)nLn) {
                                    uint32_t wS = scratch.hitWeights[j] > 1 ? std::min((uint32_t)std::min(dA, aB), (uint32_t)kmerLen) / scratch.hitWeights[j] : std::min((uint32_t)std::min(dA, aB), (uint32_t)kmerLen);
                                    int32_t pnlty = (nLn > 0) ? (int32_t)(((int64_t)nDr * wS * 1024) / ((int64_t)nLn * dRscaled)) : 0;
                                    int32_t sc = (int32_t)scratch.dpDiff[j] + (int32_t)wS - pnlty;
                                    if(sc > (int32_t)scratch.dpDiff[k]) { 
                                        scratch.dpDiff[k] = std::max(1, sc); scratch.parentDiff[k] = j; 
                                        scratch.cumDriftDiff[k] = nDr; scratch.cumLenDiff[k] = nLn; 
                                        scratch.chainOccurrencesDiff[k] = scratch.chainOccurrencesDiff[j] + 1; 
                                        streakD = 0; if (gapStreakD > 0) --gapStreakD; 
                                    } else { 
                                        if(++streakD > 25) break; 
                                        if (scratch.backtrackVisitDiff[j] == (int32_t)k && ++gapStreakD > 25) break; 
                                    }
                                    if (scratch.parentDiff[j] >= 0) scratch.backtrackVisitDiff[scratch.parentDiff[j]] = (int32_t)k;
                                }
                            }
                        }
                        
                        // Update best score and extension for current pair A-B
                        if(scratch.dpSame[k] > maxScSame) { 
                            maxScSame = scratch.dpSame[k]; bestEndIdxSame = (int32_t)k; maxExtSame = getAlignmentLength(posAk, posBk); 
                        } else if (scratch.dpSame[k] == maxScSame && bestEndIdxSame >= 0) { 
                            uint64_t ext = getAlignmentLength(posAk, posBk); 
                            if (ext < maxExtSame) { bestEndIdxSame = (int32_t)k; maxExtSame = ext; } 
                        }
                        
                        if(scratch.dpDiff[k] > maxScDiff) { 
                            maxScDiff = scratch.dpDiff[k]; bestEndIdxDiff = (int32_t)k; maxExtDiff = getAlignmentLength(posAk, (uint32_t)(readLenB - 1 - posBk)); 
                        } else if (scratch.dpDiff[k] == maxScDiff && bestEndIdxDiff >= 0) { 
                            uint64_t ext = getAlignmentLength(posAk, (uint32_t)(readLenB - 1 - posBk)); 
                            if (ext < maxExtDiff) { bestEndIdxDiff = (int32_t)k; maxExtDiff = ext; } 
                        }
                    }

                    // --- Step 3: Filtering and Candidate Generation ---
                    uint32_t bestScPair = std::max(maxScSame, maxScDiff);
                    if (bestScPair < 2) continue;
                    uint32_t filterThresh = std::max<uint32_t>(
                        chainFilterMinScore,
                        uint32_t(double(bestScPair) * chainFilterRatio));

                    scratch.chainCandidates.clear();
                    for (size_t k = 0; k < numHits; ++k) {
                        if (scratch.dpSame[k] >= filterThresh) {
                            scratch.chainCandidates.push_back({scratch.dpSame[k], getAlignmentLength(scratch.hitPosA[k], scratch.hitPosB[k]), (int32_t)k, false});
                        }
                        if (scratch.dpDiff[k] >= filterThresh) {
                            scratch.chainCandidates.push_back({scratch.dpDiff[k], getAlignmentLength(scratch.hitPosA[k], (uint32_t)(readLenB - 1 - scratch.hitPosB[k])), (int32_t)k, true});
                        }
                    }
                    std::sort(scratch.chainCandidates.begin(), scratch.chainCandidates.end());

                    // Apply Hifiasm's candidate limit heuristic if needed
                    if (scratch.chainCandidates.size() > maxChainLimit && maxChainLimit > 0) {
                        scratch.candidateTypes.assign(scratch.chainCandidates.size(), 0);
                        int32_t countByType[4] = {0}, threshByType[4] = {0};
                        for (size_t ci = 0; ci < scratch.chainCandidates.size(); ++ci) {
                            const auto& cand = scratch.chainCandidates[ci];
                            int32_t rootK = cand.endK;
                            const auto& parentArr = cand.isDiff ? scratch.parentDiff : scratch.parentSame;
                            while (parentArr[rootK] != -1) rootK = parentArr[rootK];
                            
                            uint32_t qs = markersA[scratch.hitOrdinalA[rootK]].position;
                            uint32_t qe = markersA[scratch.hitOrdinalA[cand.endK]].position + (uint32_t)kmerLen;
                            int type = getOverlapType(qs, qe, (uint32_t)readLenA);
                            scratch.candidateTypes[ci] = type;
                            countByType[type]++;
                            if ((uint64_t)countByType[type] == maxChainLimit) threshByType[type] = cand.score;
                        }
                        
                        if (threshByType[0] || threshByType[1] || threshByType[2] || threshByType[3]) {
                            scratch.filteredCandidates.clear();
                            for (size_t ci = 0; ci < scratch.chainCandidates.size(); ++ci) {
                                if (scratch.chainCandidates[ci].score >= (uint32_t)threshByType[scratch.candidateTypes[ci]]) {
                                    scratch.filteredCandidates.push_back(scratch.chainCandidates[ci]);
                                }
                            }
                            scratch.chainCandidates = std::move(scratch.filteredCandidates);
                        }
                    }

                    // --- Step 4: Interval Non-redundancy & Final Alignment Extraction ---
                    scratch.acceptedIntervalsSame.clear(); scratch.acceptedIntervalsDiff.clear();
                    auto isOverlapLarge = [&](uint32_t qs, uint32_t qe, const vector<ThreadScratchpad::ChainInterval>& accepted) {
                         uint32_t len = qe - qs + 1;
                         for (const auto& iv : accepted) {
                             uint32_t oS = std::max(qs, iv.qs), oE = std::min(qe, iv.qe);
                             if (oE >= oS && (oE - oS + 1) > (uint32_t)(nonRedundantOverlapFraction * double(len))) return true;
                         }
                         return false;
                    };

                    for (const auto& cand : scratch.chainCandidates) {
                        // Filter low-support candidates before doing any backtracking/extraction.
                        // We keep only candidates with strictly more than minChainedMarkerCount marker hits
                        // so the candidate table (and downstream work) stays high-confidence.
                        if(minChainedMarkerCount > 0) {
                            const uint32_t occ = cand.isDiff ?
                                scratch.chainOccurrencesDiff[cand.endK] :
                                scratch.chainOccurrencesSame[cand.endK];
                            if(occ <= minChainedMarkerCount) {
                                continue;
                            }
                        }

                        int32_t currK = cand.endK;
                        const auto& parentArr = cand.isDiff ? scratch.parentDiff : scratch.parentSame;
                        scratch.currentChainPath.clear();
                        while(currK != -1) { scratch.currentChainPath.push_back(currK); currK = parentArr[currK]; }
                        std::reverse(scratch.currentChainPath.begin(), scratch.currentChainPath.end());
                        if (scratch.currentChainPath.empty()) continue;

                        uint32_t qOrdS = scratch.hitOrdinalA[scratch.currentChainPath.front()];
                        uint32_t qOrdE = scratch.hitOrdinalA[scratch.currentChainPath.back()];
                        if (cand.isDiff) { if (isOverlapLarge(qOrdS, qOrdE, scratch.acceptedIntervalsDiff)) continue; }
                        else { if (isOverlapLarge(qOrdS, qOrdE, scratch.acceptedIntervalsSame)) continue; }

                        // Extract alignment path
                        Alignment al; 
                        al.ordinals.reserve(scratch.currentChainPath.size());
                        const OrientedReadId orientedReadB(readIdB, 0);
                        const auto& mB = markers[orientedReadB.getValue()];
                        bool validChain = true;

                        for(uint32_t idxK : scratch.currentChainPath) {
                            uint32_t tPos = scratch.hitPosB[idxK];
                            auto itB = std::lower_bound(mB.begin(), mB.end(), tPos, [](const CompressedMarker& m, uint32_t v){ return m.position < v; });
                            if(itB != mB.end() && itB->position == tPos) {
                                al.ordinals.push_back({scratch.hitOrdinalA[idxK], (uint32_t)(itB - mB.begin())});
                            } else { validChain = false; break; }
                        }

                        if(validChain) {
                             uint32_t qPstart = markersA[scratch.hitOrdinalA[scratch.currentChainPath.front()]].position;
                             uint32_t qPend = markersA[scratch.hitOrdinalA[scratch.currentChainPath.back()]].position + (uint32_t)kmerLen;
                             uint32_t bPfirst = mB[al.ordinals.front()[1]].position;
                             uint32_t bPlast = mB[al.ordinals.back()[1]].position;
                             
                             uint32_t tS, tE;
                             if (cand.isDiff) {
                                 uint32_t minB = std::min(bPfirst, bPlast), maxB = std::max(bPfirst, bPlast) + (uint32_t)kmerLen;
                                 tS = (uint32_t)readLenB - maxB; tE = (uint32_t)readLenB - minB;
                             } else {
                                 tS = bPfirst; tE = bPlast + (uint32_t)kmerLen;
                             }

                             // Optional: minimum overlap length in bases (pre-extension span).
                             if(invertedIndexData.minOverlapLength > 0) {
                                 const uint64_t qSpan = uint64_t(qPend) - uint64_t(qPstart);
                                 const uint64_t tSpan = uint64_t(tE) - uint64_t(tS);
                                 if(std::min(qSpan, tSpan) < uint64_t(invertedIndexData.minOverlapLength)) {
                                     continue;
                                 }
                             }

                             // Optional: discard internal overlaps that would require large end extension.
                             if(invertedIndexData.maxEndFuzz > 0) {
                                 const uint32_t leftNeed = std::min(qPstart, tS);
                                 const int64_t qRight = int64_t(readLenA) - int64_t(qPend);
                                 const int64_t tRight = int64_t(readLenB) - int64_t(tE);
                                 const uint32_t rightNeed = uint32_t(std::min<int64_t>(std::max<int64_t>(qRight, 0), std::max<int64_t>(tRight, 0)));
                                 if(leftNeed > invertedIndexData.maxEndFuzz || rightNeed > invertedIndexData.maxEndFuzz) {
                                     continue;
                                 }
                             }

                             uint32_t fQs = qPstart, fQe = qPend, fTs = tS, fTe = tE;
                             if (fQs <= fTs) { fTs -= fQs; fQs = 0; } else { fQs -= fTs; fTs = 0; }
                             int64_t remQ = (int64_t)readLenA - fQe, remT = (int64_t)readLenB - fTe;
                             if (remQ <= remT) { fQe = (uint32_t)readLenA; fTe += (uint32_t)remQ; } else { fTe = (uint32_t)readLenB; fQe += (uint32_t)remT; } 

                             al.qs = fQs; al.qe = fQe;
                             const uint32_t markerCountA = uint32_t(markersA.size());
                             const uint32_t markerCountB = uint32_t(mB.size());

                             bool isSameStrand = true;
                             if (cand.isDiff) {
                                 isSameStrand = false;
                                 const auto p = dinara::rcIntervalToForward(uint32_t(readLenB), fTs, fTe);
                                 al.ts = p.first;
                                 al.te = p.second;
                                 const uint32_t numMB = markerCountB;
                                 for(auto& p : al.ordinals) {
                                     p[1] = numMB - 1 - p[1];
                                 }
                                 scratch.acceptedIntervalsDiff.push_back({qOrdS, qOrdE});
                             } else {
                                 isSameStrand = true;
                                 al.ts = fTs;
                                 al.te = fTe;
                                 scratch.acceptedIntervalsSame.push_back({qOrdS, qOrdE});
                             }

                             // Canonicalize candidate so readIds[0] < readIds[1], and keep alignment consistent.
                             ReadId cand0 = readIdA;
                             ReadId cand1 = readIdB;
                             canonicalizeCandidateAndAlignment(cand0, cand1, isSameStrand, al, markerCountA, markerCountB);
                             localCandidates.push_back(OrientedReadPair(cand0, cand1, isSameStrand));
                             localAlignments.push_back(std::move(al));
                        }
                    }
                }
            }
        }
    }
};
// Explicit template instantiation for the MultithreadedObject base class.
template class MultithreadedObject<InvertedIndexFinder>;

} // namespace dinara
using namespace dinara;

// =============================================================================
// Phase 1-4: Build the inverted index for overlap candidate discovery.
// =============================================================================
void Assembler::buildInvertedIndex(uint64_t threadCount) {
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    performanceLog << timestamp << "Building Inverted Index..." << endl;

    // =========================================================================
    // Phase 1: Inverted Index Construction
    // =========================================================================
    // We build an index mapping each canonical K-mer to all of its occurrences
    // (ReadId, Position). We only index Strand 0 because canonical K-mers are
    // strand-symmetric, halving memory usage. RC matches are detected during
    // the DP chaining phase by observing decreasing target positions.
    // =========================================================================
    checkMarkersAreOpen();
    if(!markerKmerIds->isOpen()) {
        throw runtime_error("Marker KmerIds not available for Inverted Index.");
    }
    invertedIndexData.k = assemblerInfo->k;
    const uint64_t totalMarkers = markers->totalSize();
    cout << "Building Inverted Index for " << totalMarkers << " markers." << endl;

    // We use static load balancing (not dynamic getNextBatch) to ensure that
    // each thread processes the exact same reads in both Pass 1 (Count) and
    // Pass 2 (Fill). This prevents buffer overflows from mismatched counts.
    checkMarkersAreOpen();
    const ReadId readCount = ReadId(markers->size() / 2);
    vector<uint64_t> threadMarkerCounts(threadCount, 0);
    vector<size_t> threadOffsets(threadCount, 0);

    // Pass 1: Count markers per thread for memory allocation
    auto countFunction = [&](size_t threadId) {
        ReadId startRead = (ReadId)((uint64_t)readCount * threadId / threadCount);
        ReadId endRead   = (ReadId)((uint64_t)readCount * (threadId + 1) / threadCount);

        uint64_t count = 0;
        for(ReadId rId = startRead; rId != endRead; ++rId) {
            count += (*markers)[size_t(rId) << 1].size();
        }
        threadMarkerCounts[threadId] = count;
    };

    vector<std::thread> threads;
    for(size_t i = 0; i < threadCount; i++) threads.emplace_back(countFunction, i);
    for(auto& t : threads) t.join();
    threads.clear();

    // Compute global offsets
    size_t totalMarkersFound = 0;
    for(size_t i = 0; i < threadCount; i++) {
        threadOffsets[i] = totalMarkersFound;
        totalMarkersFound += threadMarkerCounts[i];
    }
    
    invertedIndexData.occurrences.resize(totalMarkersFound);
    cout << "Index allocated for " << totalMarkersFound << " occurrences (Strand 0 only)." << endl;

    // Pass 2: Fill the occurrence array
    auto fillFunction = [&](size_t threadId) {
        ReadId startRead = (ReadId)((uint64_t)readCount * threadId / threadCount);
        ReadId endRead   = (ReadId)((uint64_t)readCount * (threadId + 1) / threadCount);
        size_t writeOffset = threadOffsets[threadId];

        for(ReadId rId = startRead; rId != endRead; ++rId) {
            const auto& rMarkers = (*markers)[size_t(rId) << 1];
            const auto& rKmerIds = (*markerKmerIds)[size_t(rId) << 1];
            if(rMarkers.size() != rKmerIds.size()) continue; 

            for(size_t i = 0; i < rMarkers.size(); ++i) {
                KmerId kId = rKmerIds[i];
                KmerId rcKId = getRcKmerId(kId, invertedIndexData.k);
                KmerId canonicalKId = (kId < rcKId) ? kId : rcKId;

                invertedIndexData.occurrences[writeOffset++] = {
                    canonicalKId,
                    rId,
                    rMarkers[i].position
                };
            }
        }
    };

    // Run filling process in parallel
    for(size_t i = 0; i < threadCount; i++) {
        threads.emplace_back(fillFunction, i);
    }
    for(auto& t : threads) t.join();
    threads.clear();

    // =========================================================================
    // Phase 2: Parallel LSD Radix Sort
    // =========================================================================
    // We sort the occurrences by K-mer ID to group all markers with the same
    // K-mer together. This uses a parallel Least-Significant-Digit (LSD) radix
    // sort, which is O(N * numBytes) and stable. The parallelization strategy:
    //   1. Each thread computes a local histogram for its data partition.
    //   2. Histograms are combined to compute global bucket offsets.
    //   3. Each thread scatters its data to the correct output positions.
    // =========================================================================
    cout << "Sorting " << invertedIndexData.occurrences.size() << " occurrences..." << endl;
    
    if(!invertedIndexData.occurrences.empty()) {
        const size_t n = invertedIndexData.occurrences.size();
        const size_t numBytes = sizeof(KmerId);
        
        vector<InvertedIndexOccurrence> buffer(n);
        vector<InvertedIndexOccurrence>* src = &invertedIndexData.occurrences;
        vector<InvertedIndexOccurrence>* dst = &buffer;
        
        for (size_t byteIdx = 0; byteIdx < numBytes; ++byteIdx) {
            
            // 2.1 Parallel Histogram Calculation
            vector<vector<size_t>> threadHistograms(threadCount, vector<size_t>(256, 0));
            auto countHistFunc = [&](size_t tid) {
                 size_t start = (n * tid) / threadCount;
                 size_t end = (n * (tid + 1)) / threadCount;
                 for(size_t i = start; i < end; ++i) {
                     uint8_t byteVal = (uint8_t)(((*src)[i].kmerId >> (byteIdx * 8)) & 0xFF);
                     threadHistograms[tid][byteVal]++;
                 }
            };
            vector<thread> sortThreads;
            for(size_t i = 0; i < threadCount; i++) sortThreads.emplace_back(countHistFunc, i);
            for(auto& t : sortThreads) t.join();
            sortThreads.clear();
            
            // 2.2 Global Offset Calculation
            size_t globalBucketCounts[256] = {0};
            for(size_t b = 0; b < 256; ++b) {
                for(size_t tid = 0; tid < threadCount; ++tid) {
                    globalBucketCounts[b] += threadHistograms[tid][b];
                }
            }
            size_t globalBucketOffsets[256];
            globalBucketOffsets[0] = 0;
            for(size_t b = 1; b < 256; ++b) {
                globalBucketOffsets[b] = globalBucketOffsets[b - 1] + globalBucketCounts[b - 1];
            }
            
            // 2.3 Individual Thread Starting Offsets
            vector<vector<size_t>> writeOffsets(threadCount, vector<size_t>(256));
            for(size_t b = 0; b < 256; ++b) {
                size_t current = globalBucketOffsets[b];
                for(size_t tid = 0; tid < threadCount; tid++) {
                    writeOffsets[tid][b] = current;
                    current += threadHistograms[tid][b];
                }
            }
            
            // 2.4 Parallel Data Scattering
            auto scatterDataFunc = [&](size_t tid) {
                 size_t start = (n * tid) / threadCount;
                 size_t end = (n * (tid + 1)) / threadCount;
                 for(size_t i = start; i < end; ++i) {
                     uint8_t byteVal = (uint8_t)(((*src)[i].kmerId >> (byteIdx * 8)) & 0xFF);
                     (*dst)[writeOffsets[tid][byteVal]++] = (*src)[i];
                 }
            };
            for(size_t i = 0; i < threadCount; i++) sortThreads.emplace_back(scatterDataFunc, i);
            for(auto& t : sortThreads) t.join();
            sortThreads.clear();

            std::swap(src, dst);
        }
        
        if (src != &invertedIndexData.occurrences) {
             invertedIndexData.occurrences = std::move(*src);
        }
    }

    // =========================================================================
    // Phase 3: Hash Table Construction (Open Addressing)
    // =========================================================================
    // Build a power-of-2 sized hash table for O(1) K-mer lookups. Each entry
    // stores the starting index and count of occurrences for a single K-mer.
    // We use linear probing for collision resolution. The load factor is ~50%
    // to balance memory usage and probe length.
    // =========================================================================
    cout << "Allocating Direct Addressing Hash Table..." << endl;
    
    uint64_t numDistinctKmers = 0;
    if(!invertedIndexData.occurrences.empty()) {
        numDistinctKmers = 1;
        KmerId lastKId = invertedIndexData.occurrences[0].kmerId;
        for(size_t i = 1; i < invertedIndexData.occurrences.size(); ++i) {
            if(invertedIndexData.occurrences[i].kmerId != lastKId) {
                numDistinctKmers++;
                lastKId = invertedIndexData.occurrences[i].kmerId;
            }
        }
    }
    cout << "Distinct K-mers found: " << numDistinctKmers << endl;

    uint64_t tableSize = 1;
    while(tableSize < numDistinctKmers * 2) tableSize *= 2; 
    invertedIndexData.hashTable.resize(tableSize); 
    
    if(!invertedIndexData.occurrences.empty()) {
        KmerId currentKmer = invertedIndexData.occurrences[0].kmerId;
        uint64_t blockStart = 0;
        uint64_t mask = tableSize - 1;

        auto insertIntoTable = [&](KmerId key, uint64_t startIdx, uint32_t count) {
            uint64_t slotIdx = hashKmer(key) & mask;
            while(!invertedIndexData.hashTable[slotIdx].empty) {
                slotIdx = (slotIdx + 1) & mask;
            }
            invertedIndexData.hashTable[slotIdx] = {key, startIdx, count, false};
        };

        for(uint64_t i = 1; i < invertedIndexData.occurrences.size(); ++i) {
            if(invertedIndexData.occurrences[i].kmerId != currentKmer) {
                insertIntoTable(currentKmer, blockStart, (uint32_t)(i - blockStart));
                currentKmer = invertedIndexData.occurrences[i].kmerId;
                blockStart = i;
            }
        }
        insertIntoTable(currentKmer, blockStart, (uint32_t)(invertedIndexData.occurrences.size() - blockStart));
    }

    // =========================================================================
    // Phase 4: Index Compaction
    // =========================================================================
    // The full occurrence struct contains the K-mer ID, but after building the
    // hash table, we only need (ReadId, Position) for each occurrence. This
    // compaction reduces the memory footprint by ~3x, which is critical for
    // large datasets. We parallelize this trivially across all occurrences.
    // =========================================================================
    const size_t nMarkersTotal = invertedIndexData.occurrences.size();
    invertedIndexData.compactOccurrences.resize(nMarkersTotal);
    
    auto compactFunction = [&](size_t tid) {
        size_t start = (nMarkersTotal * tid) / threadCount;
        size_t end = (nMarkersTotal * (tid + 1)) / threadCount;
        for(size_t i = start; i < end; ++i) {
            invertedIndexData.compactOccurrences[i] = {
                invertedIndexData.occurrences[i].readId,
                invertedIndexData.occurrences[i].position
            };
        }
    };
    vector<thread> compactThreads;
    for(size_t i = 0; i < threadCount; i++) compactThreads.emplace_back(compactFunction, i);
    for(auto& t : compactThreads) t.join();
    
    // Free high-memory intermediate vector
    invertedIndexData.occurrences.clear();
    invertedIndexData.occurrences.shrink_to_fit();
    cout << "Index construction complete." << endl;
}

// =============================================================================
// Phase 5: Run DP chaining on the built index to find alignment candidates.
// =============================================================================
void Assembler::chainAlignmentCandidates(
    double maxDriftRate,
    uint64_t maxChainLimit,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
    uint32_t minChainedMarkerCount,
    uint64_t threadCount
) {
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const auto startTime = std::chrono::steady_clock::now();
    performanceLog << timestamp << "Starting DP chaining for alignment candidates." << endl;

    // Check that buildInvertedIndex has been called
    if(invertedIndexData.compactOccurrences.empty()) {
        throw runtime_error("chainAlignmentCandidates: buildInvertedIndex must be called first.");
    }

    // =========================================================================
    // Parallel Candidate Search
    // =========================================================================
    // The main search loop. For each Read A (strand 0), we:
    //   1. Look up each of its K-mers in the hash table.
    //   2. Collect all (PartnerReadId, PositionA, PositionB) hits.
    //   3. Sort hits by PartnerReadId, then by PositionA.
    //   4. For each (A, B) pair, run DP chaining to find the best chain(s).
    //   5. Extract alignment coordinates and store the candidate pair.
    // The DP chaining logic is Hifiasm-compatible, using cumulative drift
    // constraints and a penalty based on gap rate.
    // =========================================================================
    alignmentCandidates.candidates.createNew(largeDataName("AlignmentCandidates"), largeDataPageSize);
    alignmentCandidatesAlignmentsData.alignments.createNew(largeDataName("AlignmentCandidatesInvertedIndex"), largeDataPageSize); 
    
    // Safety reserve for performance
    const ReadId readCount = ReadId(markers->size() / 2);
    alignmentCandidates.candidates.reserve(size_t(readCount) * 50);
    alignmentCandidatesAlignmentsData.alignments.reserve(size_t(readCount) * 50); 

    invertedIndexData.maxDriftRate = maxDriftRate;
    invertedIndexData.coveragePeak = assemblerInfo->kmerDistributionInfo.coveragePeak;
    invertedIndexData.weightExponent = overlapCandidatesOptions.invertedIndexWeightExponent;
    invertedIndexData.lowFreqMultiplier = overlapCandidatesOptions.invertedIndexLowFreqMultiplier;
    invertedIndexData.highFreqMultiplier = overlapCandidatesOptions.invertedIndexHighFreqMultiplier;
    invertedIndexData.rareKmerWeight = overlapCandidatesOptions.invertedIndexRareKmerWeight;
    invertedIndexData.chainFilterRatio = overlapCandidatesOptions.invertedIndexChainFilterRatio;
    invertedIndexData.chainFilterMinScore = overlapCandidatesOptions.invertedIndexChainFilterMinScore;
    invertedIndexData.nonRedundantOverlapFraction = overlapCandidatesOptions.invertedIndexNonRedundantOverlapFraction;
    invertedIndexData.minOverlapLength = overlapCandidatesOptions.minOverlapLength;
    invertedIndexData.maxEndFuzz = overlapCandidatesOptions.maxEndFuzz;
    invertedIndexData.weightLut.resize(512);
    for(size_t i=0; i<invertedIndexData.weightLut.size(); i++) {
        invertedIndexData.weightLut[i] = (uint8_t)std::min(255.0, std::pow((double)i, invertedIndexData.weightExponent));
    }

    // Launch the finder across all threads
    InvertedIndexFinder finder(
        getReads(),
        *markers,
        *markerKmerIds,
        invertedIndexData,
        alignmentCandidates.candidates,
        alignmentCandidatesAlignmentsData.alignments,
        maxChainLimit,
        minChainedMarkerCount,
        threadCount
    );

    alignmentCandidates.candidates.unreserve();
    alignmentCandidatesAlignmentsData.alignments.unreserve();
    
    // --- Final Cleanup ---
    invertedIndexData.compactOccurrences.clear();
    invertedIndexData.compactOccurrences.shrink_to_fit();
    invertedIndexData.hashTable.clear();
    invertedIndexData.hashTable.shrink_to_fit();

    const auto endTime = std::chrono::steady_clock::now();
    const double totalSeconds = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime)).count());
    cout << timestamp << "Alignment discovery completed in " << totalSeconds << " s." << endl;
}

// =============================================================================
// Convenience wrapper that calls both buildInvertedIndex and chainAlignmentCandidates.
// =============================================================================
void Assembler::findAlignmentCandidatesInvertedIndex(
    double maxDriftRate,
    uint64_t maxChainLimit,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
    uint32_t minChainedMarkerCount,
    uint64_t threadCount
) {
    buildInvertedIndex(threadCount);
    chainAlignmentCandidates(maxDriftRate, maxChainLimit, overlapCandidatesOptions, minChainedMarkerCount, threadCount);
}

// =============================================================================
// Chain pre-imported PAF candidates using the inverted index.
// This assumes buildInvertedIndex has been called and alignmentCandidates.candidates
// has been populated by importAlignmentCandidatesFromPaf.
// =============================================================================
void Assembler::chainPafCandidates(
    double maxDriftRate,
    uint64_t maxChainLimit,
    const OverlapCandidatesOptions& overlapCandidatesOptions,
    uint32_t minChainedMarkerCount,
    uint64_t threadCount
) {
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const auto startTime = std::chrono::steady_clock::now();
    performanceLog << timestamp << "Starting DP chaining for PAF-imported candidates." << endl;

    // Check that buildInvertedIndex has been called
    if(invertedIndexData.compactOccurrences.empty() && invertedIndexData.hashTable.empty()) {
        throw runtime_error("chainPafCandidates: buildInvertedIndex must be called first.");
    }

    // Check that PAF candidates have been imported
    if(!alignmentCandidates.candidates.isOpen || alignmentCandidates.candidates.size() == 0) {
        throw runtime_error("chainPafCandidates: No candidates imported. Call importAlignmentCandidatesFromPaf first.");
    }

    cout << timestamp << "Chaining " << alignmentCandidates.candidates.size() << " PAF-imported candidates..." << endl;

    // Store parameters for chaining
    invertedIndexData.maxDriftRate = maxDriftRate;
    invertedIndexData.coveragePeak = assemblerInfo->kmerDistributionInfo.coveragePeak;
    invertedIndexData.weightExponent = overlapCandidatesOptions.invertedIndexWeightExponent;
    invertedIndexData.lowFreqMultiplier = overlapCandidatesOptions.invertedIndexLowFreqMultiplier;
    invertedIndexData.highFreqMultiplier = overlapCandidatesOptions.invertedIndexHighFreqMultiplier;
    invertedIndexData.rareKmerWeight = overlapCandidatesOptions.invertedIndexRareKmerWeight;
    invertedIndexData.chainFilterRatio = overlapCandidatesOptions.invertedIndexChainFilterRatio;
    invertedIndexData.chainFilterMinScore = overlapCandidatesOptions.invertedIndexChainFilterMinScore;
    invertedIndexData.nonRedundantOverlapFraction = overlapCandidatesOptions.invertedIndexNonRedundantOverlapFraction;
    invertedIndexData.minOverlapLength = overlapCandidatesOptions.minOverlapLength;
    invertedIndexData.maxEndFuzz = overlapCandidatesOptions.maxEndFuzz;
    invertedIndexData.weightLut.resize(512);
    for(size_t i=0; i<invertedIndexData.weightLut.size(); i++) {
        invertedIndexData.weightLut[i] = (uint8_t)std::min(255.0, std::pow((double)i, invertedIndexData.weightExponent));
    }

    // Create output storage for chained alignments
    alignmentCandidatesAlignmentsData.alignments.createNew(largeDataName("AlignmentCandidatesInvertedIndex"), largeDataPageSize);
    alignmentCandidatesAlignmentsData.alignments.reserve(alignmentCandidates.candidates.size());

    // Store original candidates before we modify them
    vector<OrientedReadPair> originalCandidates;
    originalCandidates.reserve(alignmentCandidates.candidates.size());
    for(size_t i = 0; i < alignmentCandidates.candidates.size(); i++) {
        originalCandidates.push_back(alignmentCandidates.candidates[i]);
    }

    // Clear and prepare to repopulate with chained results
    alignmentCandidates.candidates.clear();
    alignmentCandidates.candidates.reserve(originalCandidates.size());

    // Setup threading
    const uint64_t hashMask = invertedIndexData.hashTable.size() - 1;
    const auto* hashTablePtr = invertedIndexData.hashTable.data();
    const uint64_t kmerLen = invertedIndexData.k;
    const double maxDriftRateLocal = invertedIndexData.maxDriftRate;
    const uint64_t coveragePeak = invertedIndexData.coveragePeak;

    // Thread-local results
    vector<vector<OrientedReadPair>> threadCandidates(threadCount);
    vector<vector<Alignment>> threadAlignments(threadCount);
    
    // Per-thread processing
    const size_t batchSize = std::max(size_t(1), originalCandidates.size() / (threadCount * 10));
    setupLoadBalancing(originalCandidates.size(), batchSize);

    vector<std::thread> threads;
    for(size_t tid = 0; tid < threadCount; tid++) {
        threads.emplace_back([&, tid]() {
            ThreadScratchpad scratch;
            vector<OrientedReadPair>& localCandidates = threadCandidates[tid];
            vector<Alignment>& localAlignments = threadAlignments[tid];
            
            const int64_t dRscaled = (int64_t)(maxDriftRateLocal * 1024.0);
            
            uint64_t startBatch, endBatch;
            while(getNextBatch(startBatch, endBatch)) {
                for(size_t idx = startBatch; idx < endBatch; idx++) {
                    const OrientedReadPair& pair = originalCandidates[idx];
                    const ReadId readIdA = pair.readIds[0];
                    const ReadId readIdB = pair.readIds[1];
                    const bool pafSameStrand = pair.isSameStrand;

                    const OrientedReadId orientedReadIdA(readIdA, 0);
                    const OrientedReadId orientedReadIdB(readIdB, 0);
                    const auto& markersA = (*markers)[orientedReadIdA.getValue()];
                    const auto& kmerIdsA = (*markerKmerIds)[orientedReadIdA.getValue()];
                    const auto& markersB = (*markers)[orientedReadIdB.getValue()];
                    const auto& kmerIdsB = (*markerKmerIds)[orientedReadIdB.getValue()];
                    
                    const size_t numMarkersA = std::min(markersA.size(), kmerIdsA.size());
                    const size_t numMarkersB = std::min(markersB.size(), kmerIdsB.size());
                    
                    if(numMarkersA == 0 || numMarkersB == 0) continue;

                    const uint64_t readLenA = reads->getReadRawSequenceLength(readIdA);
                    const uint64_t readLenB = reads->getReadRawSequenceLength(readIdB);

                    // Build position lookup for read B
                    std::map<uint32_t, uint32_t> posBToOrdinalB;
                    for(size_t j = 0; j < numMarkersB; j++) {
                        posBToOrdinalB[markersB[j].position] = (uint32_t)j;
                    }

                    // Collect k-mer matches between the pair
                    scratch.clear();
                    scratch.flatHits.reserve(numMarkersA);

                    for(size_t i = 0; i < numMarkersA; i++) {
                        KmerId currentKId = kmerIdsA[i];
                        KmerId rcKId = getRcKmerId(currentKId, kmerLen);
                        KmerId canonicalKId = (currentKId < rcKId) ? currentKId : rcKId;
                        
                        uint64_t slotIdx = hashKmer(canonicalKId) & hashMask;
                        const uint32_t posA = markersA[i].position;

                        // Search hash table for this k-mer
                        while(!hashTablePtr[slotIdx].empty) {
                            if(hashTablePtr[slotIdx].key == canonicalKId) {
                                const uint64_t startIdx = hashTablePtr[slotIdx].start;
                                const uint32_t count = hashTablePtr[slotIdx].count;
                                
                                uint8_t hitWeight = 1;
                                const uint64_t lowFreq = (uint64_t)(double(coveragePeak) * invertedIndexData.lowFreqMultiplier);
                                const uint64_t highFreq = (uint64_t)(double(coveragePeak) * invertedIndexData.highFreqMultiplier);
                                
                                if (count <= std::max(2UL, lowFreq)) {
                                    hitWeight = uint8_t(std::min<uint32_t>(255U, invertedIndexData.rareKmerWeight));
                                } else if (count >= std::max(3UL, highFreq)) {
                                    uint32_t w = 1 + (uint32_t)((count + (highFreq * 2) - 1) / (highFreq * 2 == 0 ? 1 : highFreq * 2));
                                    hitWeight = (w < 512) ?
                                        invertedIndexData.weightLut[w] :
                                        (uint8_t)std::min(255U, (uint32_t)pow((double)w, invertedIndexData.weightExponent));
                                }

                                // Find occurrences in read B only
                                const auto* compactOccs = &invertedIndexData.compactOccurrences[startIdx];
                                for(uint32_t j = 0; j < count; j++) {
                                    if(compactOccs[j].readId == readIdB) {
                                        scratch.flatHits.push_back({readIdB, posA, compactOccs[j].position, (uint32_t)i, hitWeight});
                                    }
                                }
                                break;
                            }
                            slotIdx = (slotIdx + 1) & hashMask;
                        }
                    }

                    if(scratch.flatHits.empty()) {
                        // No shared marker k-mers found, so there is no chain to score.
                        // Do not emit a candidate with an empty chain.
                        continue;
                    }

                    // Sort hits by position
                    std::sort(scratch.flatHits.begin(), scratch.flatHits.end(), 
                        [](const InvertedIndexTempHit& a, const InvertedIndexTempHit& b) {
                            return a.posA < b.posA;
                        });

                    const size_t numHits = scratch.flatHits.size();

                    // Transfer to SoA
                    scratch.hitPosA.assign(numHits, 0);
                    scratch.hitPosB.assign(numHits, 0);
                    scratch.hitOrdinalA.assign(numHits, 0);
                    scratch.hitWeights.assign(numHits, 0);
                    for(size_t k = 0; k < numHits; k++) {
                        scratch.hitPosA[k] = scratch.flatHits[k].posA;
                        scratch.hitPosB[k] = scratch.flatHits[k].posB;
                        scratch.hitOrdinalA[k] = scratch.flatHits[k].ordinalA;
                        scratch.hitWeights[k] = scratch.flatHits[k].weight;
                    }

                    // Run DP chaining (same-strand only if PAF says same strand, else opposite)
                    scratch.dpSame.assign(numHits, 0);
                    scratch.dpDiff.assign(numHits, 0);
                    scratch.parentSame.assign(numHits, -1);
                    scratch.parentDiff.assign(numHits, -1);
                    scratch.cumDriftSame.assign(numHits, 0);
                    scratch.cumDriftDiff.assign(numHits, 0);
                    scratch.cumLenSame.assign(numHits, 0);
                    scratch.cumLenDiff.assign(numHits, 0);

                    uint32_t maxScSame = 0, maxScDiff = 0;
                    int32_t bestEndIdxSame = -1, bestEndIdxDiff = -1;

                    for(size_t k = 0; k < numHits; k++) {
                        const uint32_t posAk = scratch.hitPosA[k], posBk = scratch.hitPosB[k];
                        const uint32_t initSc = std::max(1U, scratch.hitWeights[k] > 1 ? (uint32_t)kmerLen / scratch.hitWeights[k] : (uint32_t)kmerLen);
                        scratch.dpSame[k] = scratch.dpDiff[k] = initSc;

                        for(int32_t j = (int32_t)k - 1; j >= 0; --j) {
                            int32_t dA = (int32_t)posAk - (int32_t)scratch.hitPosA[j];
                            int32_t dB = (int32_t)posBk - (int32_t)scratch.hitPosB[j];

                            // Same strand case
                            if(dB > 0 && dA > 0) {
                                uint32_t nDr = scratch.cumDriftSame[j] + std::abs(dB - dA);
                                uint32_t nLn = scratch.cumLenSame[j] + uint32_t(dA);
                                if((int64_t)((uint64_t)nDr << 10) <= (int64_t)dRscaled * (int64_t)nLn) {
                                    uint32_t wS = std::min((uint32_t)std::min(dA, dB), (uint32_t)kmerLen);
                                    int32_t sc = (int32_t)scratch.dpSame[j] + (int32_t)wS;
                                    if(sc > (int32_t)scratch.dpSame[k]) {
                                        scratch.dpSame[k] = std::max(1, sc);
                                        scratch.parentSame[k] = j;
                                        scratch.cumDriftSame[k] = nDr;
                                        scratch.cumLenSame[k] = nLn;
                                    }
                                }
                            }
                            // Opposite strand case
                            else if(dB < 0 && dA > 0) {
                                int32_t aB = -dB;
                                uint32_t nDr = scratch.cumDriftDiff[j] + std::abs(aB - dA);
                                uint32_t nLn = scratch.cumLenDiff[j] + uint32_t(dA);
                                if((int64_t)((uint64_t)nDr << 10) <= (int64_t)dRscaled * (int64_t)nLn) {
                                    uint32_t wS = std::min((uint32_t)std::min(dA, aB), (uint32_t)kmerLen);
                                    int32_t sc = (int32_t)scratch.dpDiff[j] + (int32_t)wS;
                                    if(sc > (int32_t)scratch.dpDiff[k]) {
                                        scratch.dpDiff[k] = std::max(1, sc);
                                        scratch.parentDiff[k] = j;
                                        scratch.cumDriftDiff[k] = nDr;
                                        scratch.cumLenDiff[k] = nLn;
                                    }
                                }
                            }
                        }

                        if(scratch.dpSame[k] > maxScSame) { maxScSame = scratch.dpSame[k]; bestEndIdxSame = (int32_t)k; }
                        if(scratch.dpDiff[k] > maxScDiff) { maxScDiff = scratch.dpDiff[k]; bestEndIdxDiff = (int32_t)k; }
                    }

                    // Extract best chain based on PAF strand info
                    // Enforce the orientation given by the PAF record.
                    // This avoids ambiguous cases where the same/diff DP scores are similar due to
                    // canonical k-mer matching (which does not encode strand).
                    bool useSameStrand = pafSameStrand;
                    int32_t bestEndIdx = useSameStrand ? bestEndIdxSame : bestEndIdxDiff;
                    const auto& parentArr = useSameStrand ? scratch.parentSame : scratch.parentDiff;

                    if(bestEndIdx < 0) {
                        continue;
                    }

                    // Backtrack to get chain
                    scratch.currentChainPath.clear();
                    int32_t currK = bestEndIdx;
                    while(currK != -1) {
                        scratch.currentChainPath.push_back(currK);
                        currK = parentArr[currK];
                    }
                    std::reverse(scratch.currentChainPath.begin(), scratch.currentChainPath.end());

                    if(minChainedMarkerCount > 0 &&
                        scratch.currentChainPath.size() <= size_t(minChainedMarkerCount)) {
                        continue;
                    }

                    // Build alignment
                    Alignment al;
                    al.ordinals.reserve(scratch.currentChainPath.size());
                    bool validChain = true;

                    for(uint32_t idxK : scratch.currentChainPath) {
                        uint32_t tPos = scratch.hitPosB[idxK];
                        auto it = posBToOrdinalB.find(tPos);
                        if(it != posBToOrdinalB.end()) {
                            al.ordinals.push_back({scratch.hitOrdinalA[idxK], it->second});
                        } else {
                            validChain = false;
                            break;
                        }
                    }

                    if(validChain && !al.ordinals.empty()) {
                        // Set alignment coordinates.
                        // Note: Alignment stores ts/te in forward-read coordinates, even for reverse overlaps.
                        uint32_t qPstart = markersA[al.ordinals.front()[0]].position;
                        uint32_t qPend = markersA[al.ordinals.back()[0]].position + (uint32_t)kmerLen;
                        uint32_t tPstart = markersB[al.ordinals.front()[1]].position;
                        uint32_t tPend = markersB[al.ordinals.back()[1]].position + (uint32_t)kmerLen;

                        al.qs = qPstart;
                        al.qe = qPend;
                        if(useSameStrand) {
                            al.ts = tPstart;
                            al.te = tPend;
                        } else {
                            // In the diff-strand chain, the target marker positions are decreasing.
                            // Convert the (strand-0) forward interval to the canonical forward interval.
                            const uint32_t minB = std::min(tPstart, tPend >= uint32_t(kmerLen) ? (tPend - uint32_t(kmerLen)) : tPend);
                            const uint32_t maxB = std::max(tPstart, tPend >= uint32_t(kmerLen) ? (tPend - uint32_t(kmerLen)) : tPend) + uint32_t(kmerLen);
                            const uint32_t tsExt = uint32_t(readLenB) - maxB;
                            const uint32_t teExt = uint32_t(readLenB) - minB;
                            const auto p = dinara::rcIntervalToForward(uint32_t(readLenB), tsExt, teExt);
                            al.ts = p.first;
                            al.te = p.second;
                        }

                        // Optional: minimum overlap length in bases (pre-extension span).
                        if(invertedIndexData.minOverlapLength > 0) {
                            const uint64_t qSpan = uint64_t(al.qe) - uint64_t(al.qs);
                            const uint64_t tSpan = uint64_t(al.te) - uint64_t(al.ts);
                            if(std::min(qSpan, tSpan) < uint64_t(invertedIndexData.minOverlapLength)) {
                                continue;
                            }
                        }

                        // Optional: discard internal overlaps that would require large end extension.
                        if(invertedIndexData.maxEndFuzz > 0) {
                            // For reverse overlaps, Alignment stores ts/te in forward coordinates (ts < te),
                            // but the "end extension" heuristic must be computed in the overlap orientation.
                            // In the reverse-complement coordinate system, the target interval becomes:
                            // [lenB - te, lenB - ts], so the left overhang is (lenB - te) and the right overhang is ts.
                            const int64_t qRight = int64_t(readLenA) - int64_t(al.qe);
                            const uint32_t qRightNeed = uint32_t(std::max<int64_t>(qRight, 0));

                            uint32_t tLeftNeed = 0;
                            uint32_t tRightNeed = 0;
                            if(useSameStrand) {
                                tLeftNeed = al.ts;
                                const int64_t tRight = int64_t(readLenB) - int64_t(al.te);
                                tRightNeed = uint32_t(std::max<int64_t>(tRight, 0));
                            } else {
                                tLeftNeed = uint32_t(readLenB) - al.te;
                                tRightNeed = al.ts;
                            }

                            const uint32_t leftNeed = std::min(al.qs, tLeftNeed);
                            const uint32_t rightNeed = std::min(qRightNeed, tRightNeed);
                            if(leftNeed > invertedIndexData.maxEndFuzz || rightNeed > invertedIndexData.maxEndFuzz) {
                                continue;
                            }
                        }

                        // Flip ordinals for opposite strand
                        if(!useSameStrand) {
                            uint32_t numMB = (uint32_t)markersB.size();
                            for(auto& p : al.ordinals) p[1] = numMB - 1 - p[1];
                        }

                        // Canonicalize candidate so readIds[0] < readIds[1], and keep alignment consistent.
                        const uint32_t markerCountA = uint32_t(markersA.size());
                        const uint32_t markerCountB = uint32_t(markersB.size());
                        ReadId cand0 = readIdA;
                        ReadId cand1 = readIdB;
                        bool isSameStrand = useSameStrand;
                        canonicalizeCandidateAndAlignment(cand0, cand1, isSameStrand, al, markerCountA, markerCountB);

                        localCandidates.push_back(OrientedReadPair(cand0, cand1, isSameStrand));
                        localAlignments.push_back(std::move(al));
                    }
                }
            }
        });
    }

    for(auto& t : threads) t.join();

    // Merge results
    for(size_t tid = 0; tid < threadCount; tid++) {
        for(const auto& c : threadCandidates[tid]) {
            alignmentCandidates.candidates.push_back(c);
        }
        for(const auto& a : threadAlignments[tid]) {
            alignmentCandidatesAlignmentsData.alignments.push_back(a);
        }
    }

    alignmentCandidates.candidates.unreserve();
    alignmentCandidatesAlignmentsData.alignments.unreserve();

    // Cleanup inverted index data
    invertedIndexData.compactOccurrences.clear();
    invertedIndexData.compactOccurrences.shrink_to_fit();
    invertedIndexData.hashTable.clear();
    invertedIndexData.hashTable.shrink_to_fit();

    const auto endTime = std::chrono::steady_clock::now();
    const double totalSeconds = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime)).count());
    cout << timestamp << "PAF candidate chaining completed in " << totalSeconds << " s." << endl;
    cout << timestamp << "Chained " << alignmentCandidates.candidates.size() << " candidates." << endl;
}
