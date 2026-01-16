#include "Assembler.hpp"
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

using namespace std;
// using namespace dinara; // Removed to avoid ambiguity, we use explicit namespace block.

// Include Template Implementation.
#include "MultithreadedObject.tpp"
#include "Alignment.hpp"

namespace dinara {

// Helper struct (Global scope).
struct InvertedIndexTempHit {
    ReadId partnerReadId;
    uint32_t posA;
    uint32_t posB;
    uint32_t ordinalA;
    uint8_t weight;
    
    bool operator<(const InvertedIndexTempHit& other) const {
        if (partnerReadId != other.partnerReadId) return partnerReadId < other.partnerReadId;
        return posA < other.posA;
    }
};

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
        uint64_t threadCount
    ) : 
        MultithreadedObject(*this),
        reads(reads),
        markers(markers),
        markerKmerIds(markerKmerIds),
        invertedIndexData(invertedIndexData),
        candidates(candidates),
        precomputedAlignments(precomputedAlignments),
        threadCount(threadCount)
    {
        // 1. Setup Work Areas.
        threadCandidates.resize(threadCount);
        for(auto& v : threadCandidates) v.reserve(10000);
        threadAlignments.resize(threadCount);
        for(auto& v : threadAlignments) v.reserve(10000);

        // 2. Setup Load Balancing.
        // Safety: Ensure we use markers.size()/2 for read count to match mapped file.
        const ReadId readCount = ReadId(markers.size() / 2);
        setupLoadBalancing(readCount, 100);

        // 3. Run Threads.
        runThreads(&InvertedIndexFinder::threadFunction, threadCount);

        // 4. Merge Results.
        size_t totalCandidates = 0;
        for(const auto& v : threadCandidates) totalCandidates += v.size();
        
        cout << "Deep Parity: Found " << totalCandidates << " candidates." << endl;

        size_t candidateWriteOffset = candidates.size();
        candidates.resize(candidateWriteOffset + totalCandidates);
        size_t alignmentWriteOffset = precomputedAlignments.size();
        precomputedAlignments.resize(alignmentWriteOffset + totalCandidates);
        
        for(size_t i=0; i<threadCount; i++) {
            const auto& v = threadCandidates[i];
            const auto& a = threadAlignments[i];
            if(!v.empty()) {
                std::copy(v.begin(), v.end(), candidates.begin() + candidateWriteOffset);
                candidateWriteOffset += v.size();
                std::copy(a.begin(), a.end(), precomputedAlignments.begin() + alignmentWriteOffset);
                alignmentWriteOffset += a.size();
            }
        }
        
        // 5. Cleanup.
        threadCandidates.clear(); // Free memory
        threadAlignments.clear();
    }

private:
    const Reads& reads;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MemoryMapped::VectorOfVectors<KmerId, uint64_t>& markerKmerIds;
    const Assembler::AlignmentCandidatesInvertedIndexData& invertedIndexData;
    MemoryMapped::Vector<OrientedReadPair>& candidates;
    MemoryMapped::Vector<Alignment>& precomputedAlignments;
    uint64_t threadCount;

    vector<vector<OrientedReadPair>> threadCandidates;
    vector<vector<Alignment>> threadAlignments;

    void threadFunction(size_t threadId) {
        vector<OrientedReadPair>& localCandidates = threadCandidates[threadId];
        vector<Alignment>& localAlignments = threadAlignments[threadId];
        
        uint64_t mask = invertedIndexData.hashTable.size() - 1;
        const auto* hashTablePtr = invertedIndexData.hashTable.data();
        
        auto hashKmer = [&](KmerId k) -> uint64_t {
             const uint64_t* p = reinterpret_cast<const uint64_t*>(&k);
             uint64_t k1 = p[0];
             uint64_t k2 = sizeof(KmerId) > 8 ? p[1] : 0; 
             return k1 ^ (k2 + 0x9e3779b9 + (k1<<6) + (k1>>2));
        };

        // Reuse vectors inside loop (defensive).
        uint64_t begin, end;
        while(getNextBatch(begin, end)) {
            for(ReadId readIdA=ReadId(begin); readIdA!=ReadId(end); ++readIdA) {
                
                const OrientedReadId orientedReadIdA(readIdA, 0);
                const auto& markersA = markers[orientedReadIdA.getValue()];
                const auto& kmerIdsA = markerKmerIds[orientedReadIdA.getValue()]; // Safe Mapped access?
                const size_t numMarkers = std::min(markersA.size(), kmerIdsA.size()); // Defensive Clamp
                
                // Local vectors.
                vector<InvertedIndexTempHit> flatHits;
                flatHits.reserve(numMarkers * 2);
                
                vector<uint32_t> dpSame;
                vector<uint32_t> dpDiff;
                vector<int32_t> parentSame;
                vector<int32_t> parentDiff;
                vector<uint32_t> cumulativeDriftSame;
                vector<uint32_t> cumulativeDriftDiff;
                vector<uint32_t> cumulativeLengthSame;
                vector<uint32_t> cumulativeLengthDiff;

                // 1. Collect Hits.
                for(size_t i=0; i<numMarkers; ++i) {
                    // Canonical Query
                    Kmer kmer(kmerIdsA[i], invertedIndexData.k);
                    KmerId rcKmerId = KmerId(kmer.reverseComplement(invertedIndexData.k).id(invertedIndexData.k));
                    KmerId canonicalKmerId = kmerIdsA[i] < rcKmerId ? kmerIdsA[i] : rcKmerId;

                    uint64_t h = hashKmer(canonicalKmerId) & mask;
                    const KmerId kmerId = canonicalKmerId;
                    const uint32_t posA = markersA[i].position;

                    uint64_t idx = h;
                    while(!hashTablePtr[idx].empty) {
                        if(hashTablePtr[idx].key == kmerId) {
                            uint64_t start = hashTablePtr[idx].start;
                            uint32_t count = hashTablePtr[idx].count;
                            
                            // Iterate occurrences (compact).
                             const auto* occPtr = &invertedIndexData.compactOccurrences[start];
                             // Avoid bounds check in inner loop for speed, assuming valid Index build.
                             for(uint32_t j=0; j<count; ++j) {
                                 const auto& occurrence = occPtr[j];
                                 if(occurrence.readId != readIdA) {
                                      uint32_t val = count > 255 ? 255 : count;
                                      flatHits.push_back({occurrence.readId, posA, occurrence.position, (uint32_t)i, (uint8_t)val}); // Store count as weight
                                 }
                             }
                            break;
                        }
                        idx = (idx + 1) & mask;
                    }
                }
                
                // 2. Sort.
                if(flatHits.empty()) continue;
                std::sort(flatHits.begin(), flatHits.end());

                // 3. DP.
                size_t i = 0;
                while(i < flatHits.size()) {
                    ReadId readIdB = flatHits[i].partnerReadId;
                    
                    if(readIdB <= readIdA) {
                        while(i < flatHits.size() && flatHits[i].partnerReadId == readIdB) i++;
                        continue;
                    }

                    // Define group.
                    size_t start = i;
                    while(i < flatHits.size() && flatHits[i].partnerReadId == readIdB) i++;
                    size_t endGroup = i;
                    size_t numHits = endGroup - start;
                    
                    if(numHits < invertedIndexData.minMarkerCount) continue;

                    // Init DP.
                    dpSame.assign(numHits, 0); 
                    dpDiff.assign(numHits, 0);
                    parentSame.assign(numHits, -1);
                    parentDiff.assign(numHits, -1);
                    cumulativeDriftSame.assign(numHits, 0);
                    cumulativeDriftDiff.assign(numHits, 0);
                    cumulativeLengthSame.assign(numHits, 0); 
                    cumulativeLengthDiff.assign(numHits, 0); 

                    uint32_t maxChainSame = 0;
                    uint32_t maxChainDiff = 0;
                    int32_t bestIdxSame = -1;
                    int32_t bestIdxDiff = -1;
                    
                    const uint32_t kmerLength = (uint32_t)invertedIndexData.k;
                    const double bandwidthPenaltyFactor = 1.0; 
                    const int64_t driftRateInt = (int64_t)(invertedIndexData.maxDriftRate * 1024.0);
                    
                    for(size_t k=0; k<numHits; ++k) {
                        size_t idx = start + k;
                        
                        // Seed.
                        uint32_t count = flatHits[idx].weight; // We stored count in 'weight' field
                        uint32_t score = count > 1 ? kmerLength / count : kmerLength; 
                        if (score == 0) score = 1;

                        dpSame[k] = score;
                        dpDiff[k] = score;
                        cumulativeLengthSame[k] = kmerLength;
                        cumulativeLengthDiff[k] = kmerLength;

                        // Chain.
                        const size_t MAX_SKIP = 25;
                        size_t n_skip_same = 0;
                        size_t n_skip_diff = 0;
                        
                        for(size_t localJ=k-1; localJ!=SIZE_MAX; --localJ) {
                            size_t idxJ = start + localJ;
                            int32_t deltaA = (int32_t)flatHits[idx].posA - (int32_t)flatHits[idxJ].posA;
                            int32_t deltaB = (int32_t)flatHits[idx].posB - (int32_t)flatHits[idxJ].posB;
                            
                            if(deltaB > 0) { // Same Strand
                                if (deltaA == 0) continue; // Ensure strictly increasing posA
                                
                                int32_t drift = std::abs(deltaB - deltaA);
                                uint32_t newCumulativeDrift = cumulativeDriftSame[localJ] + drift;
                                uint32_t newCumulativeLength = cumulativeLengthSame[localJ] + deltaA;
                                
                                int32_t distance_min = std::min(deltaA, std::abs(deltaB));
                                uint32_t baseScore = std::min((uint32_t)distance_min, kmerLength);
                                
                                // Apply Hifiasm Weighting: score = base_score / count
                                uint32_t weightedScore = count > 1 ? baseScore / count : baseScore;
                                if (weightedScore == 0 && baseScore > 0) weightedScore = 1;
                                
                                double gap_rate = (newCumulativeLength > 0) ? ((double)newCumulativeDrift / (double)newCumulativeLength) : 0.0;
                                int32_t penalty = (int32_t)(gap_rate * weightedScore * bandwidthPenaltyFactor / 1024.0);
                                
                                if((newCumulativeDrift << 10) <= driftRateInt * newCumulativeLength) {
                                    int32_t candidateScore = (int32_t)dpSame[localJ] + (int32_t)weightedScore - penalty;
                                    if(candidateScore > (int32_t)dpSame[k]) {
                                        dpSame[k] = (uint32_t)std::max(1, candidateScore);
                                        parentSame[k] = (int32_t)localJ; // Update Parent
                                        cumulativeDriftSame[k] = newCumulativeDrift;
                                        cumulativeLengthSame[k] = newCumulativeLength;
                                        n_skip_same = 0;
                                    } else {
                                        if(++n_skip_same > MAX_SKIP) break;
                                    }
                                }
                            } else if (deltaB < 0) { // Diff Strand
                                if (deltaA == 0) continue; // Ensure strictly increasing posA

                                int32_t absDeltaB = -deltaB;
                                int32_t drift = std::abs(absDeltaB - deltaA);
                                uint32_t newCumulativeDrift = cumulativeDriftDiff[localJ] + drift;
                                uint32_t newCumulativeLength = cumulativeLengthDiff[localJ] + deltaA;
                                
                                int32_t distance_min = std::min(deltaA, absDeltaB);
                                uint32_t baseScore = std::min((uint32_t)distance_min, kmerLength);
                                
                                // Apply Hifiasm Weighting
                                uint32_t weightedScore = count > 1 ? baseScore / count : baseScore;
                                if (weightedScore == 0 && baseScore > 0) weightedScore = 1;

                                double gap_rate = (newCumulativeLength > 0) ? ((double)newCumulativeDrift / (double)newCumulativeLength) : 0.0;
                                int32_t penalty = (int32_t)(gap_rate * weightedScore * bandwidthPenaltyFactor / 1024.0);

                                if((newCumulativeDrift << 10) <= driftRateInt * newCumulativeLength) {
                                    int32_t candidateScore = (int32_t)dpDiff[localJ] + (int32_t)weightedScore - penalty;
                                    if(candidateScore > (int32_t)dpDiff[k]) {
                                        dpDiff[k] = (uint32_t)std::max(1, candidateScore);
                                        parentDiff[k] = (int32_t)localJ; // Update Parent
                                        cumulativeDriftDiff[k] = newCumulativeDrift;
                                        cumulativeLengthDiff[k] = newCumulativeLength;
                                        n_skip_diff = 0;
                                    } else {
                                        if(++n_skip_diff > MAX_SKIP) break;
                                    }
                                }
                            }
                        }
                        
                        if(dpSame[k] > maxChainSame) {
                            maxChainSame = dpSame[k];
                            bestIdxSame = (int32_t)k;
                        }
                        if(dpDiff[k] > maxChainDiff) {
                            maxChainDiff = dpDiff[k];
                            bestIdxDiff = (int32_t)k;
                        }
                    }

                    // Hifiasm-style Iterative Chaining:
                    // 1. Identify "peaks" (endpoints) with score >= 0.8 * bestScore.
                    // 2. Sort by score descending.
                    // 3. Greedily keep chains that don't overlap > 50% with accepted chains on Query (Read A).

                    uint32_t bestScore = std::max(maxChainSame, maxChainDiff);
                    if (bestScore < invertedIndexData.minMarkerCount) continue;

                    uint32_t threshold = (uint32_t)(bestScore * 0.80);
                    if (threshold < (uint32_t)invertedIndexData.minMarkerCount) threshold = (uint32_t)invertedIndexData.minMarkerCount;

                    struct ChainCandidate {
                        uint32_t score;
                        int32_t endK;
                        bool isDiff; // true if Diff, false if Same
                        bool operator<(const ChainCandidate& other) const {
                            return score > other.score; // Descending
                        }
                    };
                    vector<ChainCandidate> candidates;

                    // Collect candidates from Same
                    for (size_t k = 0; k < numHits; ++k) {
                        if (dpSame[k] >= threshold) {
                            // Only add local maxima (peaks) to reduce redundancy? 
                            // Hifiasm checks if extending improves score.
                            // Here, we take endpoints. If k extends k-1, k has higher score usually.
                            // Simple approach: Take all above threshold, sort will prioritize ends.
                            candidates.push_back({dpSame[k], (int32_t)k, false});
                        }
                    }

                    // Collect candidates from Diff
                    for (size_t k = 0; k < numHits; ++k) {
                        if (dpDiff[k] >= threshold) {
                            candidates.push_back({dpDiff[k], (int32_t)k, true});
                        }
                    }

                    std::sort(candidates.begin(), candidates.end());

                    // Accepted chains (Query Interval [qs, qe])
                    struct Interval { uint32_t qs; uint32_t qe; };
                    vector<Interval> acceptedIntervalsSame;
                    vector<Interval> acceptedIntervalsDiff;

                    // Hifiasm Overlap Threshold: 50% of the shorter chain (or new chain)
                    auto isOverlapping = [&](uint32_t qs, uint32_t qe, const vector<Interval>& existing) {
                         uint32_t len = qe - qs + 1;
                         for (const auto& iv : existing) {
                             uint32_t overlapStart = std::max(qs, iv.qs);
                             uint32_t overlapEnd = std::min(qe, iv.qe);
                             if (overlapEnd >= overlapStart) {
                                  uint32_t overlapLen = overlapEnd - overlapStart + 1;
                                  if (overlapLen > (uint32_t)(0.5 * len)) return true; // >50% overlap of NEW chain rejected
                                  // Note: Hifiasm might check both, but shielding new (shorter/worse) chains is the goal.
                             }
                         }
                         return false;
                    };


                    for (const auto& cand : candidates) {
                         // Reconstruct chain to find qs
                         int32_t curr = cand.endK;
                         const auto& parents = cand.isDiff ? parentDiff : parentSame;
                         
                         vector<uint32_t> chainIndices;
                         while(curr != -1) {
                             chainIndices.push_back(curr);
                             curr = parents[curr];
                         }
                         std::reverse(chainIndices.begin(), chainIndices.end());
                         
                         if (chainIndices.empty()) continue;

                         // Get Query Range
                         uint32_t qs = flatHits[start + chainIndices.front()].ordinalA;
                         uint32_t qe = flatHits[start + chainIndices.back()].ordinalA;
                         
                         // Check overlap
                         if (cand.isDiff) {
                             if (isOverlapping(qs, qe, acceptedIntervalsDiff)) continue;
                         } else {
                             if (isOverlapping(qs, qe, acceptedIntervalsSame)) continue;
                         }

                         // Valid non-overlapping chain! Construct Alignment.
                         Alignment alignment;
                         alignment.ordinals.reserve(chainIndices.size());
                         const OrientedReadId orientedReadIdB(readIdB, 0); 
                         const auto& markersB = markers[orientedReadIdB.getValue()];

                         bool chainValid = true;
                         for(size_t idxK_i=0; idxK_i < chainIndices.size(); ++idxK_i) {
                              int32_t idxK = chainIndices[idxK_i];
                              const auto& hit = flatHits[start + idxK];
                              uint32_t ordA = hit.ordinalA;
                              uint32_t targetPosB = hit.posB;

                              auto itB = std::lower_bound(markersB.begin(), markersB.end(), targetPosB, 
                                  [](const CompressedMarker& m, uint32_t val){ return m.position < val; });
                              
                              if(itB != markersB.end() && itB->position == targetPosB) {
                                  uint32_t ordB = (uint32_t)(itB - markersB.begin());
                                  alignment.ordinals.push_back({ordA, ordB});
                              } else {
                                  chainValid = false; break;
                              }
                         }

                         if(chainValid) {
                             if (cand.isDiff) {
                                 // Adjust for RC
                                 uint32_t totalMarkersB = (uint32_t)markersB.size();
                                 for(auto& pair : alignment.ordinals) {
                                     pair[1] = totalMarkersB - 1 - pair[1];
                                 }
                                 acceptedIntervalsDiff.push_back({qs, qe});
                                 localCandidates.push_back(OrientedReadPair(readIdA, readIdB, false)); // Diff
                             } else {
                                 acceptedIntervalsSame.push_back({qs, qe});
                                 localCandidates.push_back(OrientedReadPair(readIdA, readIdB, true)); // Same
                             }
                             localAlignments.push_back(std::move(alignment));
                         }
                    }
                }
            }
        }
    }
};
// Explicit instantiation.
template class MultithreadedObject<InvertedIndexFinder>;

} // namespace dinara
using namespace dinara; // For Assembler access below

// Optimized Packed Occurrence Structure.
// 24 bytes (16 for KmerId + 4 ReadId + 4 Pos).
// Standard ReadId is uint32_t. Position is uint32_t.
// So we can store this efficiently. 
// We separate KmerId (Key) from Payload (Value) for sorting if we use Structure of Arrays?
// No, sorting requires moving Key and Value together.
// Using Array of Structures (AoS) is cache efficient for sorting.

// Wait. markerKmerIds stores Kmers separate from markers (Pos).
// We are merging them here.

// Hifiasm uses Binning to Parallelize sort.
// We implement a simplified Parallel Fill -> Serial Sort (std::sort).
// Given std::sort is highly optimized, it is sufficient.

void Assembler::findAlignmentCandidatesInvertedIndex(
    uint64_t minMarkerCount, 
    double maxDriftRate,     
    uint64_t threadCount
) {
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const auto tBegin = std::chrono::steady_clock::now();
    performanceLog << timestamp << "Finding alignment candidates using Inverted Index." << endl;

    // 1. Build the Inverted Index.
    // We assume markerKmerIds and markers are populated and consistent.
    checkMarkersAreOpen();
    if(!markerKmerIds.isOpen()) {
        throw runtime_error("Marker KmerIds not available for Inverted Index.");
    }
    
    // Initialize k
    invertedIndexData.k = assemblerInfo->k;

    // Allocate huge vector.
    // Total markers = markers.totalSize().
    const uint64_t totalMarkers = markers.totalSize();
    cout << "Building Inverted Index for " << totalMarkers << " markers." << endl;
    // invertedIndexData.occurrences.resize(totalMarkers); // Initial resize removed, now exact resize after counting.

    // Populate occurrences.
    // Can be parallelized per read batch. But vector is contiguous.
    // Serial populate is fast (sequential memory). Parallel populate requires pre-calculating offsets.
    // markers.size() gives vector of vectors size (2*readCount).
    // markers stores ReadId implicitly by index?
    // markers[i] is for OrientedRead i?
    // OrderedReadId(readId, 0).
    // We only need one strand?
    // If we match ReadA (Strand0) to ReadB (Strand0) -> Same Strand.
    // If we match ReadA (Strand0) to ReadB (Strand1) -> Opposite Strand?
    // markers stores both strands.
    // Hifiasm usually indexes Canonical K-mers?
    // KmerId is Canonical.
    // If KmerId is Canonical, occurrences from Strand0 and Strand1 are mathematically related.
    // Storing occurrences from ONLY Strand 0 is sufficient if we handle RC?
    // MarkerKmerIds stores Canonical KmerId.
    // A marker at Pos X on Strand 0 with KmerId K.
    // At Pos L-X-1 on Strand 1, it has the SAME KmerId K (because K is canonical).
    // So storing occurrences for BOTH strands is redundant?
    // If we index Strand 0 only:
    // Query Read A (Strand 0).
    // Match Read B (Strand 0) -> Same Strand.
    // Read C(0) sequence: K3, K2, K1. -> We detect "Reverse Match".
    // So "Same Strand" vs "Opposite Strand" is detected by order of matches (Chain).
    // (Increasing vs Decreasing positions).
    // So we ONLY need to store occurrences from Strand 0!
    // This halves the index size.
    // markers[i] where i is even (0, 2, ...) is Strand 0.

    // Parallel Fill.
    // We divide reads among threads.
    // Each thread writes to a separate chunk of the huge vector?
    // No, we don't know the exact count per thread beforehand unless we perform a prefix sum.
    // Step 1: Count markers per thread.
    // Step 2: Calculate offsets.
    // Step 3: Fill.

    // Using `setupLoadBalancing` logic.
    // Ensure readCount matches markers size to avoid OOB access in MemoryMapped vector.
    // -------------------------------------------------------------------------
    // Static Load Balancing for 2-Pass Algorithm
    // We MUST use static partitioning to ensure Thread X handles the EXACT same
    // set of reads in Pass 1 (Count) and Pass 2 (Fill).
    // Dynamic load balancing (getNextBatch) would cause buffer overflows/corruption
    // because a thread might get fewer markers in Pass 1 and more in Pass 2.
    // -------------------------------------------------------------------------

    checkMarkersAreOpen();
    const ReadId readCount = ReadId(markers.size() / 2);
    vector<uint64_t> threadMarkerCounts(threadCount, 0);
    vector<size_t> threadOffsets(threadCount, 0);

    // Pass 1: Count Markers
    auto countMarkers = [&](size_t threadId) {
        // Static Partitioning
        ReadId start = (ReadId)((uint64_t)readCount * threadId / threadCount);
        ReadId end   = (ReadId)((uint64_t)readCount * (threadId + 1) / threadCount);

        uint64_t count = 0;
        for(ReadId readId=start; readId!=end; ++readId) {
            // Only Strand 0. markers stores oriented reads.
            // 2 * readid = Strand 0.
            count += markers[size_t(readId) << 1].size();
        }
        threadMarkerCounts[threadId] = count;
    };

    // Run Pass 1
    vector<std::thread> threads;
    for(size_t i=0; i<threadCount; i++) {
        threads.emplace_back(countMarkers, i);
    }
    for(auto& t : threads) t.join();
    threads.clear();

    // Prefix sum
    size_t currentOffset = 0;
    for(size_t i=0; i<threadCount; i++) {
        threadOffsets[i] = currentOffset;
        currentOffset += threadMarkerCounts[i];
    }
    
    // Resize exact
    invertedIndexData.occurrences.resize(currentOffset);
    cout << "Index contains " << invertedIndexData.occurrences.size() << " occurrences (Strand 0 only)." << endl;

    // Pass 2: Fill Markers
    auto fillMarkers = [&](size_t threadId) {
        // Same Partitioning
        ReadId start = (ReadId)((uint64_t)readCount * threadId / threadCount);
        ReadId end   = (ReadId)((uint64_t)readCount * (threadId + 1) / threadCount);
        
        size_t offset = threadOffsets[threadId];

        for(ReadId readId=start; readId!=end; ++readId) {
            const auto& readMarkers = markers[size_t(readId) << 1];
            const auto& readKmerIds = markerKmerIds[size_t(readId) << 1];
            
            if(readMarkers.size() != readKmerIds.size()) {
                 continue; 
            }

            for(size_t i=0; i<readMarkers.size(); ++i) {
                // Bounds check optimization: verify offset < occurrences.size()?
                // Should be safe due to deterministic counting.
                // Canonical Index Fill
                Kmer kmer(readKmerIds[i], invertedIndexData.k);
                KmerId rcKmerId = KmerId(kmer.reverseComplement(invertedIndexData.k).id(invertedIndexData.k));
                KmerId canonicalKmerId = readKmerIds[i] < rcKmerId ? readKmerIds[i] : rcKmerId;

                invertedIndexData.occurrences[offset++] = {
                    canonicalKmerId,
                    readId,
                    readMarkers[i].position
                };
            }
        }
    };

    // Run Pass 2
    for(size_t i=0; i<threadCount; i++) {
        threads.emplace_back(fillMarkers, i);
    }
    for(auto& t : threads) t.join();
    threads.clear();


    // Sort.
    cout << "Sorting Inverted Index (Parallel Radix Sort)..." << endl;
    
    if(!invertedIndexData.occurrences.empty()) {
        const size_t n = invertedIndexData.occurrences.size();
        const size_t numBytes = sizeof(KmerId);
        
        // 1. Parallel Radix Sort.
        // Divide into chunks for parallel counting?
        // Or simpler: Just Parallel MSB Radix Sort (Binning).
        // Since we want O(N) global stability, LSD is easier but harder to parallelize without barrier.
        // Let's implement Parallel LSD Radix Sort with Thread-Local Histograms.
        
        vector<InvertedIndexOccurrence> buffer(n);
        vector<InvertedIndexOccurrence>* src = &invertedIndexData.occurrences;
        vector<InvertedIndexOccurrence>* dst = &buffer;
        
        for (size_t byteIdx = 0; byteIdx < numBytes; ++byteIdx) {
            
            // Step 1: Parallel Count.
            vector<vector<size_t>> histograms(threadCount, vector<size_t>(256, 0));
            auto countFunc = [&](size_t tid) {
                 size_t start = (n * tid) / threadCount;
                 size_t end = (n * (tid + 1)) / threadCount;
                 for(size_t i=start; i<end; ++i) {
                     uint8_t b = (uint8_t)(((*src)[i].kmerId >> (byteIdx * 8)) & 0xFF);
                     histograms[tid][b]++;
                 }
            };
            vector<thread> sortThreads;
            for(size_t i=0; i<threadCount; i++) sortThreads.emplace_back(countFunc, i);
            for(auto& t : sortThreads) t.join();
            sortThreads.clear();
            
            // Step 2: Global Offsets.
            size_t globalCounts[256] = {0};
            for(size_t b=0; b<256; ++b) {
                for(size_t tid=0; tid<threadCount; ++tid) {
                    globalCounts[b] += histograms[tid][b];
                }
            }
            size_t globalOffsets[256];
            globalOffsets[0] = 0;
            for(size_t b=1; b<256; ++b) globalOffsets[b] = globalOffsets[b-1] + globalCounts[b-1];
            
            // Calculate thread-specific start offsets for each bucket.
            vector<vector<size_t>> threadOffsets(threadCount, vector<size_t>(256));
            for(size_t b=0; b<256; ++b) {
                size_t current = globalOffsets[b];
                for(size_t tid=0; tid<threadCount; ++tid) {
                    threadOffsets[tid][b] = current;
                    current += histograms[tid][b];
                }
            }
            
            // Step 3: Parallel Scatter.
            auto scatterFunc = [&](size_t tid) {
                 size_t start = (n * tid) / threadCount;
                 size_t end = (n * (tid + 1)) / threadCount;
                 for(size_t i=start; i<end; ++i) {
                     uint8_t b = (uint8_t)(((*src)[i].kmerId >> (byteIdx * 8)) & 0xFF);
                     (*dst)[threadOffsets[tid][b]++] = (*src)[i];
                 }
            };
            for(size_t i=0; i<threadCount; i++) sortThreads.emplace_back(scatterFunc, i);
            for(auto& t : sortThreads) t.join();
            sortThreads.clear();

            std::swap(src, dst);
        }
        
        if (src != &invertedIndexData.occurrences) {
             invertedIndexData.occurrences = std::move(*src);
        }
    }

    // Build Lookup Map (Open Addressing Hash Table).
    cout << "Building Linear Probing Hash Table..." << endl;
    
    uint64_t numDistinct = 0;
    if(!invertedIndexData.occurrences.empty()) {
        numDistinct = 1;
        KmerId last = invertedIndexData.occurrences[0].kmerId;
        for(size_t i=1; i<invertedIndexData.occurrences.size(); ++i) {
            if(invertedIndexData.occurrences[i].kmerId != last) {
                numDistinct++;
                last = invertedIndexData.occurrences[i].kmerId;
            }
        }
    }
    cout << "Distinct K-mers: " << numDistinct << endl;

    uint64_t tableSize = 1;
    while(tableSize < numDistinct * 2) tableSize *= 2; 
    cout << "Table Size: " << tableSize << endl;
    
    invertedIndexData.hashTable.resize(tableSize); 
    
    auto hashKmer = [&](KmerId k) -> uint64_t {
        uint64_t h = 0;
        const uint64_t* p = reinterpret_cast<const uint64_t*>(&k);
        uint64_t k1 = p[0];
        uint64_t k2 = sizeof(KmerId) > 8 ? p[1] : 0;
        h = k1 ^ (k2 + 0x9e3779b9 + (k1<<6) + (k1>>2)); 
        return h;
    };

    if(!invertedIndexData.occurrences.empty()) {
        KmerId currentKmer = invertedIndexData.occurrences[0].kmerId;
        uint64_t start = 0;
        uint64_t mask = tableSize - 1;

        auto insert = [&](KmerId key, uint64_t startIdx, uint32_t count) {
            uint64_t idx = hashKmer(key) & mask;
            while(!invertedIndexData.hashTable[idx].empty) {
                idx = (idx + 1) & mask;
            }
            invertedIndexData.hashTable[idx] = {key, startIdx, count, false};
        };

        for(uint64_t i=1; i<invertedIndexData.occurrences.size(); ++i) {
            if(invertedIndexData.occurrences[i].kmerId != currentKmer) {
                uint32_t count = (uint32_t)(i - start);
                insert(currentKmer, start, count);
                
                currentKmer = invertedIndexData.occurrences[i].kmerId;
                start = i;
            }
        }
        uint32_t count = (uint32_t)(invertedIndexData.occurrences.size() - start);
        insert(currentKmer, start, count);
    }
    cout << "Hash Table built." << endl;

    // OPTIMIZATION: Convert to Compact Index (SoAish).
    cout << "Compacting Index..." << endl;
    invertedIndexData.compactOccurrences.resize(invertedIndexData.occurrences.size());
    for(size_t i=0; i<invertedIndexData.occurrences.size(); ++i) {
        invertedIndexData.compactOccurrences[i] = {
            invertedIndexData.occurrences[i].readId,
            invertedIndexData.occurrences[i].position
        };
    }
    
    // Free heavy vector.
    invertedIndexData.occurrences.clear();
    invertedIndexData.occurrences.shrink_to_fit();
    cout << "Index Compacted." << endl;



    // 2. Compute Candidates (Parallel).
    alignmentCandidates.candidates.createNew(largeDataName("AlignmentCandidates"), largeDataPageSize);
    alignmentCandidatesAlignmentsData.alignments.createNew(largeDataName("AlignmentCandidatesInvertedIndex"), largeDataPageSize); // Added
    
    // Reserve estimation.
    checkMarkersAreOpen();
    
    alignmentCandidates.candidates.reserve(size_t(readCount) * 50);
    alignmentCandidatesAlignmentsData.alignments.reserve(size_t(readCount) * 50); // Added

    invertedIndexData.minMarkerCount = minMarkerCount;
    invertedIndexData.maxDriftRate = maxDriftRate;

    // Use InvertedIndexFinder pattern.
    InvertedIndexFinder finder(
        getReads(),
        markers,
        markerKmerIds,
        invertedIndexData,
        alignmentCandidates.candidates,
        alignmentCandidatesAlignmentsData.alignments, // Added
        threadCount
    );

    alignmentCandidates.candidates.unreserve();
    alignmentCandidatesAlignmentsData.alignments.unreserve(); // Added
    
    // Cleanup Index.
    invertedIndexData.compactOccurrences.clear();
    invertedIndexData.compactOccurrences.shrink_to_fit();
    invertedIndexData.hashTable.clear();
    invertedIndexData.hashTable.shrink_to_fit();

    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    performanceLog << timestamp << "Inverted Index Candidate Finding (Optimized) completed in " << tTotal << " s." << endl;
}
