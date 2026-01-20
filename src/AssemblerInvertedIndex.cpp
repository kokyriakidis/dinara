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

// Hifiasm parity: Overlap type classification for max_n_chain filtering.
// Type 0: Left overhang (query starts at position 0)
// Type 1: Right overhang (query ends at read length)
// Type 2: Contained (query fully spans read)
// Type 3: Containing (neither end at boundary)
static int getOverlapType(uint32_t qs, uint32_t qe, uint32_t queryLen) {
    if (qs == 0 && qe >= queryLen - 1) return 2;  // Contained
    else if (qs > 0 && qe < queryLen - 1) return 3;  // Containing
    else return (qs == 0) ? 0 : 1;  // Left or Right overhang
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
    uint64_t maxChainLimit;
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
                const uint64_t readLenA = reads.getReadRawSequenceLength(readIdA);
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
                    
                    // Hifiasm Parity: usage of minMarkerCount is replaced by MinScore = K.
                    // We must allow small hit counts (e.g. 1) if they generate a high score (Unique Marker).
                    if(numHits == 0) continue;

                    // --- Hifiasm Parity: Gradient Scoring Logic (ha_get_new_candidates) ---
                    // Calculate dynamic thresholds based on Coverage Peak.
                    // HA_KMER_GOOD_RATIO = 0.333 (anchor.cpp:11)
                    uint64_t peak = invertedIndexData.coveragePeak;
                    uint64_t low_occ = (uint64_t)(peak * 0.333);
                    if (low_occ < 2) low_occ = 2;
                    uint64_t high_occ = (uint64_t)(peak * 1.667); // 2.0 - 0.333
                    if (high_occ < low_occ) high_occ = low_occ + 1; // Safety

                    for(size_t k=0; k<numHits; ++k) {
                        size_t idx = start + k;
                        // flatHits.weight currently holds the Global Count (from Index).
                        uint32_t cnt = flatHits[idx].weight;
                        
                        if (cnt > low_occ && cnt < high_occ) {
                            flatHits[idx].weight = 1; // Good (Unique-ish)
                        } else if (cnt <= low_occ) {
                            flatHits[idx].weight = 2; // Low (Rare, slight penalty? Hifiasm logic)
                        } else {
                            // High (Repetitive Penalty)
                            // Hifiasm: cnt = 1 + ((cnt + 2*high - 1) / (2*high))
                            // Then pow(cnt, 1.1)
                            // Normal_w divides score by this weight. Higher weight = Lower Score.
                            uint64_t twoHigh = high_occ * 2;
                            if (twoHigh == 0) twoHigh = 1; // Safety
                            uint32_t w = 1 + (uint32_t)((cnt + twoHigh - 1) / twoHigh);
                            flatHits[idx].weight = (uint32_t)pow((double)w, 1.1);
                        }
                    }

                    // Init DP.
                    dpSame.assign(numHits, 0); 
                    dpDiff.assign(numHits, 0);
                    parentSame.assign(numHits, -1);
                    parentDiff.assign(numHits, -1);
                    cumulativeDriftSame.assign(numHits, 0);
                    cumulativeDriftDiff.assign(numHits, 0);
                    cumulativeLengthSame.assign(numHits, 0); 
                    cumulativeLengthDiff.assign(numHits, 0);

                    // Hifiasm Parity: tmp arrays for chain skip tracking
                    vector<int32_t> tmpSame(numHits, -1);
                    vector<int32_t> tmpDiff(numHits, -1);

                    // Hifiasm Parity: occ arrays for chain length counting (dp->occ)
                    vector<uint32_t> occSame(numHits, 1);
                    vector<uint32_t> occDiff(numHits, 1); 

                    uint32_t maxChainSame = 0;
                    uint32_t maxChainDiff = 0;
                    int32_t bestIdxSame = -1;
                    int32_t bestIdxDiff = -1;
                    // Hifiasm Parity: tie-breaking by chainLen for determinism
                    uint64_t minChainLenSame = UINT64_MAX;
                    uint64_t minChainLenDiff = UINT64_MAX;
                    
                    // Hifiasm Parity: get_chainLen - computes overlap length after coordinate extension
                    // Used for tie-breaking when scores are equal (prefer shorter chainLen)
                    // readLenA is hoisted outside loop
                    const uint64_t readLenB = reads.getReadRawSequenceLength(readIdB);
                    auto get_chainLen = [&](uint32_t posA, uint32_t posB) -> uint64_t {
                        // Extend start: whoever has shorter overhang extends to 0
                        uint32_t x_beg = posA;
                        uint32_t y_beg = posB;
                        if (x_beg <= y_beg) {
                            y_beg = y_beg - x_beg;
                            x_beg = 0;
                        } else {
                            x_beg = x_beg - y_beg;
                            y_beg = 0;
                        }
                        // Extend end: whoever has shorter remaining extends to read length
                        uint64_t x_right = readLenA - posA - 1;
                        uint64_t y_right = readLenB - posB - 1;
                        uint32_t x_end, y_end;
                        if (x_right <= y_right) {
                            x_end = (uint32_t)(readLenA - 1);
                            y_end = posB + (uint32_t)x_right;
                        } else {
                            x_end = posA + (uint32_t)y_right;
                            y_end = (uint32_t)(readLenB - 1);
                        }
                        return x_end - x_beg + 1;
                    };
                    
                    const uint32_t kmerLength = (uint32_t)invertedIndexData.k;
                    // Hifiasm Parity: bandwidth_penalty = 1 / band_width_threshold
                    const double bandwidthPenalty = 1.0 / invertedIndexData.maxDriftRate;
                    const int64_t driftRateInt = (int64_t)(invertedIndexData.maxDriftRate * 1024.0);

                    // Hifiasm Parity: ha_chain_check fast path
                    // If all hits are perfectly collinear, skip full DP.
                    bool canSkipDP = true;
                    if (numHits > 1) {
                        for (size_t checkK = 1; checkK < numHits && canSkipDP; ++checkK) {
                            int32_t dA = (int32_t)flatHits[start + checkK].posA - (int32_t)flatHits[start + checkK - 1].posA;
                            int32_t dB = (int32_t)flatHits[start + checkK].posB - (int32_t)flatHits[start + checkK - 1].posB;
                            if (dA <= 0 || dB <= 0 || std::abs(dB - dA) > (int32_t)(invertedIndexData.maxDriftRate * dA)) {
                                canSkipDP = false;
                            }
                        }
                    }

                    if (canSkipDP && numHits > 0) {
                        // Hifiasm Parity: ha_chain_check - compute actual scores with weighting and penalty
                        // First hit
                        uint32_t count0 = flatHits[start].weight;
                        uint32_t score0 = count0 > 1 ? kmerLength / count0 : kmerLength;
                        if (score0 == 0) score0 = 1;
                        dpSame[0] = score0;
                        parentSame[0] = -1;
                        cumulativeDriftSame[0] = 0;
                        cumulativeLengthSame[0] = 0;
                        occSame[0] = 1;
                        
                        int32_t tot_indel = 0;
                        int32_t tot_len = 0;
                        const int32_t THRESHOLD_MAX_SIZE = 31;  // Hifiasm constant
                        bool fastPathValid = true;
                        
                        for (size_t checkK = 1; checkK < numHits && fastPathValid; ++checkK) {
                            int32_t dx = (int32_t)flatHits[start + checkK].posA - (int32_t)flatHits[start + checkK - 1].posA;
                            int32_t dy = (int32_t)flatHits[start + checkK].posB - (int32_t)flatHits[start + checkK - 1].posB;
                            int32_t dd = std::abs(dx - dy);
                            tot_indel += dd;
                            tot_len += dy;
                            
                            // Hifiasm Parity: cumulative drift check (line 876)
                            if (dy <= 0 || tot_indel > tot_len * invertedIndexData.maxDriftRate) {
                                fastPathValid = false;
                                break;
                            }
                            
                            int32_t dg = std::min(dx, dy);
                            
                            // Hifiasm Parity: THRESHOLD_MAX_SIZE check (line 878)
                            if (dd > THRESHOLD_MAX_SIZE && dd > dg * invertedIndexData.maxDriftRate) {
                                fastPathValid = false;
                                break;
                            }
                            
                            int32_t baseScore = std::min((uint32_t)dg, kmerLength);
                            // Hifiasm uses a[i].cnt (Local Count)
                            uint32_t countK = flatHits[start + checkK].weight;
                            uint32_t weightedScore = countK > 1 ? baseScore / countK : baseScore;
                            if (weightedScore == 0 && baseScore > 0) weightedScore = 1;
                            
                            double gap_rate = (tot_len > 0) ? ((double)tot_indel / (double)tot_len) : 0.0;
                            int32_t penalty = (int32_t)(gap_rate * weightedScore * bandwidthPenalty);
                            
                            dpSame[checkK] = dpSame[checkK - 1] + weightedScore - penalty;
                            parentSame[checkK] = (int32_t)(checkK - 1);
                            cumulativeDriftSame[checkK] = tot_indel;
                            cumulativeLengthSame[checkK] = tot_len;
                            occSame[checkK] = (uint32_t)(checkK + 1);
                        }
                        
                        if (!fastPathValid) {
                            // Fast path failed, fall through to full DP
                            goto full_dp;
                        }
                        
                        // Set best for Same strand
                        maxChainSame = dpSame[numHits - 1];
                        bestIdxSame = (int32_t)(numHits - 1);
                        
                        // Skip full DP, go directly to chain extraction
                        goto skip_dp_same;
                    }
                    
                    full_dp:  // Label for fallback when fast path fails
                    for(size_t k=0; k<numHits; ++k) {
                        size_t idx = start + k;
                        
                        // Seed.
                        // Hifiasm Parity: Local Count
                        uint32_t count = flatHits[idx].weight;
                        uint32_t score = count > 1 ? kmerLength / count : kmerLength; 
                        if (score == 0) score = 1;

                        dpSame[k] = score;
                        dpDiff[k] = score;
                        // Hifiasm Parity: self_length starts at 0, not kmerLength!
                        // dp->self_length[i] = max_self_length (which is 0 initially)
                        cumulativeLengthSame[k] = 0;
                        cumulativeLengthDiff[k] = 0;

                        // Chain.
                        // Hifiasm Parity: calculate_overlap_region_by_chaining (Hash_Table.cpp:1140) uses 25.
                        // calculate_ug_chaining uses 50, but we are doing overlap generation.
                        const size_t MAX_SKIP = 25;
                        size_t n_max_skip_same = 0;
                        size_t n_max_skip_diff = 0;
                        // Hifiasm Parity: chain skip counters
                        size_t n_chn_skip_same = 0;
                        size_t n_chn_skip_diff = 0;
                        
                        for(size_t localJ=k-1; localJ!=SIZE_MAX; --localJ) {
                            size_t idxJ = start + localJ;
                            int32_t deltaA = (int32_t)flatHits[idx].posA - (int32_t)flatHits[idxJ].posA;
                            int32_t deltaB = (int32_t)flatHits[idx].posB - (int32_t)flatHits[idxJ].posB;
                            
                            if(deltaB > 0) { // Same Strand
                                // Hifiasm Parity: complete skip condition
                                if (deltaB == 0 || deltaA <= 0) continue;
                                
                                int32_t drift = std::abs(deltaB - deltaA);
                                uint32_t newCumulativeDrift = cumulativeDriftSame[localJ] + drift;
                                // Hifiasm Parity: use deltaB (distance_self_pos), not deltaA
                                uint32_t newCumulativeLength = cumulativeLengthSame[localJ] + std::abs(deltaB);
                                
                                int32_t distance_min = std::min(deltaA, std::abs(deltaB));
                                uint32_t baseScore = std::min((uint32_t)distance_min, kmerLength);
                                
                                // Hifiasm Parity: score = base_score / a[j].cnt (PREDECESSOR's count, not current!)
                                // Use Local Count
                                uint32_t predecessorCount = flatHits[idxJ].weight;
                                uint32_t weightedScore = predecessorCount > 1 ? baseScore / predecessorCount : baseScore;
                                if (weightedScore == 0 && baseScore > 0) weightedScore = 1;
                                
                                // Hifiasm Parity: penalty = gap_rate * score * (1/band_width)
                                double gap_rate = (newCumulativeLength > 0) ? ((double)newCumulativeDrift / (double)newCumulativeLength) : 0.0;
                                int32_t penalty = (int32_t)(gap_rate * weightedScore * bandwidthPenalty);
                                
                                if(((uint64_t)newCumulativeDrift << 10) <= driftRateInt * newCumulativeLength) {
                                    int32_t candidateScore = (int32_t)dpSame[localJ] + (int32_t)weightedScore - penalty;
                                    if(candidateScore > (int32_t)dpSame[k]) {
                                        dpSame[k] = (uint32_t)std::max(1, candidateScore);
                                        parentSame[k] = (int32_t)localJ;
                                        cumulativeDriftSame[k] = newCumulativeDrift;
                                        cumulativeLengthSame[k] = newCumulativeLength;
                                        // Hifiasm Parity: update occ (chain length)
                                        occSame[k] = occSame[localJ] + 1;
                                        n_max_skip_same = 0;
                                        if (n_chn_skip_same > 0) --n_chn_skip_same; // Hifiasm Parity
                                    } else {
                                        if(++n_max_skip_same > MAX_SKIP) break;
                                        // Hifiasm Parity: chain skip check
                                        if (tmpSame[localJ] == (int32_t)k) {
                                            if (++n_chn_skip_same > MAX_SKIP) break;
                                        }
                                    }
                                    // Hifiasm Parity: propagate tmp for chain tracking
                                    if (parentSame[localJ] >= 0) tmpSame[parentSame[localJ]] = (int32_t)k;
                                }
                            } else if (deltaB < 0) { // Diff Strand
                                // Hifiasm Parity: complete skip condition
                                if (deltaA <= 0) continue;

                                int32_t absDeltaB = -deltaB;
                                int32_t drift = std::abs(absDeltaB - deltaA);
                                uint32_t newCumulativeDrift = cumulativeDriftDiff[localJ] + drift;
                                // Hifiasm Parity: use absDeltaB (distance_self_pos), not deltaA
                                uint32_t newCumulativeLength = cumulativeLengthDiff[localJ] + absDeltaB;
                                
                                int32_t distance_min = std::min(deltaA, absDeltaB);
                                uint32_t baseScore = std::min((uint32_t)distance_min, kmerLength);
                                
                                // Hifiasm Parity: score = base_score / a[j].cnt (PREDECESSOR's count, not current!)
                                // Use Local Count
                                uint32_t predecessorCount = flatHits[idxJ].weight;
                                uint32_t weightedScore = predecessorCount > 1 ? baseScore / predecessorCount : baseScore;
                                if (weightedScore == 0 && baseScore > 0) weightedScore = 1;

                                // Hifiasm Parity: penalty = gap_rate * score * (1/band_width)
                                double gap_rate = (newCumulativeLength > 0) ? ((double)newCumulativeDrift / (double)newCumulativeLength) : 0.0;
                                int32_t penalty = (int32_t)(gap_rate * weightedScore * bandwidthPenalty);

                                if(((uint64_t)newCumulativeDrift << 10) <= driftRateInt * newCumulativeLength) {
                                    int32_t candidateScore = (int32_t)dpDiff[localJ] + (int32_t)weightedScore - penalty;
                                    if(candidateScore > (int32_t)dpDiff[k]) {
                                        dpDiff[k] = (uint32_t)std::max(1, candidateScore);
                                        parentDiff[k] = (int32_t)localJ;
                                        cumulativeDriftDiff[k] = newCumulativeDrift;
                                        cumulativeLengthDiff[k] = newCumulativeLength;
                                        // Hifiasm Parity: update occ (chain length)
                                        occDiff[k] = occDiff[localJ] + 1;
                                        n_max_skip_diff = 0;
                                        if (n_chn_skip_diff > 0) --n_chn_skip_diff; // Hifiasm Parity
                                    } else {
                                        if(++n_max_skip_diff > MAX_SKIP) break;
                                        // Hifiasm Parity: chain skip check
                                        if (tmpDiff[localJ] == (int32_t)k) {
                                            if (++n_chn_skip_diff > MAX_SKIP) break;
                                        }
                                    }
                                    // Hifiasm Parity: propagate tmp for chain tracking
                                    if (parentDiff[localJ] >= 0) tmpDiff[parentDiff[localJ]] = (int32_t)k;
                                }
                            }
                        }
                        
                        // Hifiasm Parity: tie-breaking by chainLen for determinism
                        // When scores are equal, prefer shorter chainLen
                        if(dpSame[k] > maxChainSame) {
                            maxChainSame = dpSame[k];
                            bestIdxSame = (int32_t)k;
                            minChainLenSame = get_chainLen(flatHits[start + k].posA, flatHits[start + k].posB);
                        } else if (dpSame[k] == maxChainSame && bestIdxSame >= 0) {
                            uint64_t thisChainLen = get_chainLen(flatHits[start + k].posA, flatHits[start + k].posB);
                            if (thisChainLen < minChainLenSame) {
                                bestIdxSame = (int32_t)k;
                                minChainLenSame = thisChainLen;
                            }
                        }
                        if(dpDiff[k] > maxChainDiff) {
                            maxChainDiff = dpDiff[k];
                            bestIdxDiff = (int32_t)k;
                            // For diff strand, flip posB to simulate parallel alignment (Hifiasm parity)
                            uint32_t posB_flipped = (uint32_t)(readLenB - 1 - flatHits[start + k].posB);
                            minChainLenDiff = get_chainLen(flatHits[start + k].posA, posB_flipped);
                        } else if (dpDiff[k] == maxChainDiff && bestIdxDiff >= 0) {
                            uint32_t posB_flipped = (uint32_t)(readLenB - 1 - flatHits[start + k].posB);
                            uint64_t thisChainLen = get_chainLen(flatHits[start + k].posA, posB_flipped);
                            if (thisChainLen < minChainLenDiff) {
                                bestIdxDiff = (int32_t)k;
                                minChainLenDiff = thisChainLen;
                            }
                        }
                    }

                    skip_dp_same:  // Hifiasm Parity: label for fast path

                    // Hifiasm-style Iterative Chaining:
                    // 1. Identify "peaks" (endpoints) with score >= 0.8 * bestScore.
                    // 2. Sort by score descending.
                    // 3. Greedily keep chains that don't overlap > 50% with accepted chains on Query (Read A).

                    // --- Hifiasm Parity: Exact Filtering Logic ---
                    // 1. Chains must have Score >= kmerLength (equivalent to 1 unique marker or ~16 repetitive markers).
                    // 2. We do NOT filter by minMarkerCount anymore, as Hifiasm accepts single unique markers.
                    
                    uint32_t bestScore = std::max(maxChainSame, maxChainDiff);
                    // Hifiasm Parity: chain_cutoff = 2 
                    // Condition: if (score < 2) drop; => Keeps score >= 2.
                    if (bestScore < 2) continue; 

                    uint32_t threshold = (uint32_t)(bestScore * 0.80);
                    // Ensure threshold is at least >2 (so 3)
                    if (threshold <= 2) threshold = 3;


                    struct ChainCandidate {
                        uint32_t score;
                        uint64_t chainLen;  // Hifiasm Parity: for deterministic tie-breaking
                        int32_t endK;
                        bool isDiff; // true if Diff, false if Same
                        bool operator<(const ChainCandidate& other) const {
                            if (score != other.score) return score > other.score; // Descending by score
                            return chainLen < other.chainLen; // Tie-break by shorter chainLen
                        }
                    };
                    vector<ChainCandidate> candidates;

                    // Collect candidates from Same
                    for (size_t k = 0; k < numHits; ++k) {
                        if (dpSame[k] >= threshold) {
                            uint64_t cLen = get_chainLen(flatHits[start + k].posA, flatHits[start + k].posB);
                            candidates.push_back({dpSame[k], cLen, (int32_t)k, false});
                        }
                    }

                    // Collect candidates from Diff
                    for (size_t k = 0; k < numHits; ++k) {
                        if (dpDiff[k] >= threshold) {
                            // For diff strand, flip posB to simulate parallel alignment
                            uint32_t posB_flipped = (uint32_t)(readLenB - 1 - flatHits[start + k].posB);
                            uint64_t cLen = get_chainLen(flatHits[start + k].posA, posB_flipped);
                            candidates.push_back({dpDiff[k], cLen, (int32_t)k, true});
                        }
                    }

                    std::sort(candidates.begin(), candidates.end());

                    // --- Hifiasm Parity: max_n_chain filtering per overlap type ---
                    // If too many candidates, filter by category to keep top N per type.
                    // This prevents quadratic blowup on repetitive reads.
                    // Algorithm: (anchor.cpp:508-537)
                    // 1. Already sorted by score descending
                    // 2. Count per type, record score threshold at N-th position
                    // 3. Keep only candidates with score >= threshold for their type
                    
                    if (candidates.size() > maxChainLimit && maxChainLimit > 0) {
                        const uint64_t queryLen = readLenA;
                        const uint64_t kLen = invertedIndexData.k;
                        
                        // Cache overlap type per candidate (avoid double reconstruction)
                        vector<int> candidateTypes(candidates.size());
                        
                        // First pass: Determine type and count per category
                        int32_t n[4] = {0, 0, 0, 0};
                        uint32_t s[4] = {0, 0, 0, 0};
                        
                        for (size_t ci = 0; ci < candidates.size(); ++ci) {
                            const auto& cand = candidates[ci];
                            
                            // Reconstruct chain to get qs_base/qe_base
                            int32_t first = cand.endK;
                            const auto& parents = cand.isDiff ? parentDiff : parentSame;
                            while (parents[first] != -1) first = parents[first];
                            
                            uint32_t qs_base = markersA[flatHits[start + first].ordinalA].position;
                            uint32_t qe_base = markersA[flatHits[start + cand.endK].ordinalA].position + (uint32_t)kLen;
                            
                            int w = getOverlapType(qs_base, qe_base, (uint32_t)queryLen);
                            candidateTypes[ci] = w;
                            
                            n[w]++;
                            if ((uint64_t)n[w] == maxChainLimit) s[w] = cand.score;
                        }
                        
                        // Second pass: Filter using cached types (no reconstruction needed)
                        if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
                            vector<ChainCandidate> filteredCandidates;
                            filteredCandidates.reserve(maxChainLimit * 4);
                            
                            for (size_t ci = 0; ci < candidates.size(); ++ci) {
                                int w = candidateTypes[ci];
                                if (candidates[ci].score >= s[w]) {
                                    filteredCandidates.push_back(candidates[ci]);
                                }
                            }
                            candidates = std::move(filteredCandidates);
                        }
                    }

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
                             // --- Coordinate Extension (Hifiasm Parity) ---
                             // Compute extended coordinates here to avoid recomputation in AssemblerAlign.cpp
                             // Logic matches AssemblerAlign.cpp exactly.
                             {
                                  // 1. Get Read Lengths (Reuse values)
                                  const uint64_t len0 = readLenA;
                                  const uint64_t len1 = readLenB;
                                  const uint32_t kVal = invertedIndexData.k;

                                  // 2. Get Query Coords (Always Strand 0 / Forward)
                                  uint32_t qs_marker = markersA[flatHits[start + chainIndices.front()].ordinalA].position;
                                  uint32_t qe_marker = markersA[flatHits[start + chainIndices.back()].ordinalA].position + kVal;
                                  
                                  // 3. Get Target Coords (Strand 0 / Forward)
                                  uint32_t firstOrdB = alignment.ordinals.front()[1];
                                  uint32_t lastOrdB = alignment.ordinals.back()[1];
                                  uint32_t posB_first = markersB[firstOrdB].position;
                                  uint32_t posB_last = markersB[lastOrdB].position;
                                  
                                  uint32_t ts_marker, te_marker;
                                  
                                  if (cand.isDiff) {
                                       // Diff Strand: Work in Reverse Frame to match Hifiasm extension logic.
                                       uint32_t min_fwd = std::min(posB_first, posB_last);
                                       uint32_t max_fwd_end = std::max(posB_first, posB_last) + kVal;
                                       
                                       // Convert to Reverse Frame (Open)
                                       ts_marker = (uint32_t)len1 - max_fwd_end;
                                       te_marker = (uint32_t)len1 - min_fwd;
                                  } else {
                                       // Same Strand: Work in Forward Frame
                                       ts_marker = posB_first;
                                       te_marker = posB_last + kVal;
                                  }
                                  
                                  // 4. Extension Logic (Hifiasm / AssemblerAlign Parity)
                                  uint32_t qs_ext = qs_marker;
                                  uint32_t qe_ext = qe_marker;
                                  uint32_t ts_ext = ts_marker;
                                  uint32_t te_ext = te_marker;

                                  // Extend start
                                  if (qs_ext <= ts_ext) {
                                      ts_ext = ts_ext - qs_ext;
                                      qs_ext = 0;
                                  } else {
                                      qs_ext = qs_ext - ts_ext;
                                      ts_ext = 0;
                                  }

                                  // Extend end
                                  int64_t q_right = (int64_t)len0 - (int64_t)qe_ext;
                                  int64_t t_right = (int64_t)len1 - (int64_t)te_ext;

                                  if (q_right <= t_right) {
                                      qe_ext = (uint32_t)len0;
                                      te_ext = te_ext + (uint32_t)q_right;
                                  } else {
                                      te_ext = (uint32_t)len1;
                                      qe_ext = qe_ext + (uint32_t)t_right;
                                  }

                                  // 5. Store Results
                                  alignment.qs = qs_ext;
                                  alignment.qe = qe_ext;

                                  if (cand.isDiff) {
                                      // Convert Target back to Forward Strand (Storage Convention)
                                      alignment.ts = (uint32_t)len1 - te_ext - 1;
                                      alignment.te = (uint32_t)len1 - ts_ext - 1;
                                  } else {
                                      alignment.ts = ts_ext;
                                      alignment.te = te_ext;
                                  }
                             }

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
    double maxDriftRate,
    uint64_t maxChainLimit,
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
    if(!markerKmerIds->isOpen()) {
        throw runtime_error("Marker KmerIds not available for Inverted Index.");
    }
    
    // Initialize k
    invertedIndexData.k = assemblerInfo->k;

    // Allocate huge vector.
    // Total markers = markers.totalSize().
    const uint64_t totalMarkers = markers->totalSize();
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
    const ReadId readCount = ReadId(markers->size() / 2);
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
            count += (*markers)[size_t(readId) << 1].size();
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
            const auto& readMarkers = (*markers)[size_t(readId) << 1];
            const auto& readKmerIds = (*markerKmerIds)[size_t(readId) << 1];
            
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

    invertedIndexData.maxDriftRate = maxDriftRate;
    invertedIndexData.coveragePeak = assemblerInfo->kmerDistributionInfo.coveragePeak; // Parity

    // Use InvertedIndexFinder pattern.
    InvertedIndexFinder finder(
        getReads(),
        *markers,
        *markerKmerIds,
        invertedIndexData,
        alignmentCandidates.candidates,
        alignmentCandidatesAlignmentsData.alignments, // Added
        maxChainLimit,
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
