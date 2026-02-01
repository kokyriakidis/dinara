// PngImage.hpp must be included first because of png issues on Ubuntu 16.04.
#include "PngImage.hpp"

// Dinara.
#include "Assembler.hpp"
#include "Alignment.hpp"
#include "AlignmentGraph.hpp"
#include "Align4.hpp"
#include "AssemblerOptions.hpp"
#include "compressAlignment.hpp"
#include "performanceLog.hpp"
#include "ProjectedAlignment.hpp"
#include "AlignedEvidenceStore.hpp"
#include "Reads.hpp"
#include "span.hpp"
#include "timestamp.hpp"

// Standard libraries.
#include "chrono.hpp"
#include <cmath>
#include <filesystem>
#include "iterator.hpp"
#include "tuple.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>

using namespace dinara;
using namespace std;

void Assembler::computeAlignmentsWithEvidence(
    const AlignOptions& alignOptions,
    uint64_t threadCount
) {
    const auto tBegin = steady_clock::now();

    cout << timestamp << "Begin computing alignments with evidence for ";
    cout << alignmentCandidates.candidates.size() << " alignment candidates." << endl;

    // Check that we have what we need.
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Store parameters so they are accessible to the threads.
    auto& data = computeAlignmentsData;
    data.alignOptions = &alignOptions;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Pick the batch size for computing alignments.
    // Use the same logic as legacy to ensure similar load balancing.
    size_t batchSize = 10;
    if(batchSize > alignmentCandidates.candidates.size()/threadCount) {
        batchSize = alignmentCandidates.candidates.size()/threadCount;
    }
    if(batchSize == 0) {
        batchSize = 1;
    }

    // Prepare data structures for each thread.
    data.threadAlignmentData.assign(threadCount, {});
    data.threadCompressedAlignments.resize(threadCount);
    data.threadEvidenceStores.assign(threadCount, {});
    
    // Always resize these to avoid SIGSEGV during aggregation.
    data.threadProjectedAlignmentTime.assign(threadCount, 0.0);
    data.threadCollectionTime.assign(threadCount, 0.0);
    data.threadFilteredByErrorRate.assign(threadCount, 0);
    data.threadFilteredByErrorRateGap.assign(threadCount, 0);
    data.threadFilteredByGapCount.assign(threadCount, 0);

    // If needed for variant clustering, resize this as well.
    if (assemblerInfo->readGraphCreationMethod == 5) {
        data.threadVariantClusteringPositionPairs.resize(threadCount);
    }
    
    performanceLog << timestamp << "Alignment computation begins (Unified Chaining & Evidence path)." << endl;
    cout << timestamp << "Alignment computation begins (Unified Chaining & Evidence path)." << endl;
    
    setupLoadBalancing(alignmentCandidates.candidates.size(), batchSize);
    runThreads(&Assembler::computeAlignmentsWithEvidenceThreadFunction, threadCount);
    
    performanceLog << timestamp << "Alignment computation completed." << endl;
    cout << timestamp << "Alignment computation completed." << endl;

    // Store the alignments found by each thread into global shared memory.
    performanceLog << timestamp << "Storing the alignment found by each thread." << endl;

    // Pre-size global containers to avoid repeated remaps/reallocations.
    uint64_t totalAlignments = 0;
    uint64_t totalEvidenceIndex = 0;
    uint64_t totalSnp0 = 0;
    uint64_t totalIndel0 = 0;
    uint64_t totalSnp1 = 0;
    uint64_t totalIndel1 = 0;
    uint64_t totalSnpCheckpoints0 = 0;
    uint64_t totalSnpCheckpoints1 = 0;
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        totalAlignments += data.threadAlignmentData[threadId].size();
        const AlignedEvidenceStore& localStore = data.threadEvidenceStores[threadId];
        totalEvidenceIndex += localStore.index.size();
        totalSnp0 += localStore.snpStream0.size();
        totalIndel0 += localStore.indelStream0.size();
        totalSnp1 += localStore.snpStream1.size();
        totalIndel1 += localStore.indelStream1.size();
        totalSnpCheckpoints0 += localStore.snpCheckpoints0.size();
        totalSnpCheckpoints1 += localStore.snpCheckpoints1.size();
    }

    alignmentData.createNew(largeDataName("AlignmentData"), largeDataPageSize, 0, totalAlignments);
    compressedAlignments.createNew(largeDataName("CompressedAlignments"), largeDataPageSize);

    alignedEvidenceStore.clear();
    alignedEvidenceStore.reserve(
        totalEvidenceIndex,
        totalSnp0,
        totalSnpCheckpoints0,
        totalIndel0,
        totalSnp1,
        totalSnpCheckpoints1,
        totalIndel1
    );

    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
        const size_t idShift = alignmentData.size(); // Current global count serves as offset
        for(const AlignmentData& ad: threadAlignmentData) {
            alignmentData.push_back(ad);
            alignmentData.back().info.alignmentId += idShift;
        }

        const auto threadCompressedAlignments = data.threadCompressedAlignments[threadId];
        const auto size = threadCompressedAlignments->size();
        for(size_t i=0; i<size; i++) {
            compressedAlignments.appendVector(
                (*threadCompressedAlignments)[i].begin(),
                (*threadCompressedAlignments)[i].end()
            );
        }

        // Merge AlignedEvidenceStore (APES/TASSD)
        AlignedEvidenceStore& localStore = data.threadEvidenceStores[threadId];
        
        // Append Indexes (Adjusting offsets based on global stream sizes)
        uint64_t globalSnpOffset0 = alignedEvidenceStore.snpStream0.size();
        uint64_t globalSnpCheckpointOffset0 = alignedEvidenceStore.snpCheckpoints0.size();
        uint64_t globalIndelOffset0 = alignedEvidenceStore.indelStream0.size();
        uint64_t globalSnpOffset1 = alignedEvidenceStore.snpStream1.size();
        uint64_t globalSnpCheckpointOffset1 = alignedEvidenceStore.snpCheckpoints1.size();
        uint64_t globalIndelOffset1 = alignedEvidenceStore.indelStream1.size();
        
        for (auto& entry : localStore.index) {
            entry.snpOffset0 += globalSnpOffset0;
            entry.indelOffset0 += globalIndelOffset0;
            entry.snpCheckpointOffset0 += globalSnpCheckpointOffset0;
            entry.snpOffset1 += globalSnpOffset1;
            entry.indelOffset1 += globalIndelOffset1;
            entry.snpCheckpointOffset1 += globalSnpCheckpointOffset1;
            alignedEvidenceStore.index.push_back(entry);
        }

        // Append Streams
        alignedEvidenceStore.snpStream0.insert(
            alignedEvidenceStore.snpStream0.end(),
            localStore.snpStream0.begin(),
            localStore.snpStream0.end()
        );
        alignedEvidenceStore.snpCheckpoints0.insert(
            alignedEvidenceStore.snpCheckpoints0.end(),
            localStore.snpCheckpoints0.begin(),
            localStore.snpCheckpoints0.end()
        );
        alignedEvidenceStore.indelStream0.insert(
            alignedEvidenceStore.indelStream0.end(),
            localStore.indelStream0.begin(),
            localStore.indelStream0.end()
        );
        
        alignedEvidenceStore.snpStream1.insert(alignedEvidenceStore.snpStream1.end(), localStore.snpStream1.begin(), localStore.snpStream1.end());
        alignedEvidenceStore.snpCheckpoints1.insert(alignedEvidenceStore.snpCheckpoints1.end(), localStore.snpCheckpoints1.begin(), localStore.snpCheckpoints1.end());
        alignedEvidenceStore.indelStream1.insert(alignedEvidenceStore.indelStream1.end(), localStore.indelStream1.begin(), localStore.indelStream1.end());

        localStore.clear();
    }

    // Release unused allocated memory.
    alignmentData.unreserve();
    compressedAlignments.unreserve();

    // Create the alignment table (mapping from read to its alignments).
    // This is required for subsequent filtering phases.
    performanceLog << timestamp << "Creating alignment table." << endl;
    computeAlignmentTable();

    const auto tEnd = steady_clock::now();
    const double elapsedSeconds = seconds(tEnd - tBegin);
    performanceLog << timestamp << "Done computing alignments. Elapsed time: " << elapsedSeconds << " s." << endl;
    cout << timestamp << "Done computing alignments. Elapsed time: " << elapsedSeconds << " s." << endl;
}

void Assembler::computeAlignmentsWithEvidenceThreadFunction(size_t threadId) {
    auto& data = computeAlignmentsData;
    const AlignOptions& alignOptions = *data.alignOptions;
    auto& threadAlignmentData = data.threadAlignmentData[threadId];
    
    // Initialize compressed alignment storage for this thread.
    data.threadCompressedAlignments[threadId] = make_shared<MemoryMapped::VectorOfVectors<char, uint64_t>>();
    auto& thisThreadCompressedAlignments = *data.threadCompressedAlignments[threadId];
    thisThreadCompressedAlignments.createNew(
        largeDataName("Thread" + to_string(threadId) + "-CompressedAlignments"),
        largeDataPageSize);

    AlignedEvidenceStore& store = data.threadEvidenceStores[threadId];
    string compressedAlignment;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t candidateIndex = begin; candidateIndex != end; candidateIndex++) {
            const OrientedReadPair& candidate = alignmentCandidates.candidates[candidateIndex];
            
            // This refactored flow EXCLUSIVELY uses precomputed chains.
            // Chaining is now performed upfront during candidate generation/PAF import.
            const Alignment& alignment = alignmentCandidatesAlignmentsData.alignments[candidateIndex];
            if(alignment.ordinals.empty()) {
                continue;
            }
            // Skip low-support candidates early.
            // This avoids spending time in projected alignment construction for pairs
            // that cannot possibly meet Align.minAlignedMarkerCount.
            if(alignOptions.minAlignedMarkerCount > 0 &&
                alignment.ordinals.size() < size_t(alignOptions.minAlignedMarkerCount)) {
                continue;
            }
            array<OrientedReadId, 2> orientedReadIds = {
                OrientedReadId(candidate.readIds[0], 0),
                OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1)
            };

            // Compute projected alignment metrics and sparse diffs (mismatches/indels).
            const auto tProjStart = steady_clock::now();
            const ProjectedAlignment projectedAlignment(
                *this,
                orientedReadIds,
                alignment,
                ProjectedAlignment::Method::QuickRawSparse,
                alignOptions.overlapDpMatchScore,
                alignOptions.overlapDpMismatchScore,
                alignOptions.overlapDpGapOpen1,
                alignOptions.overlapDpGapExtend1,
                alignOptions.overlapDpGapOpen2,
                alignOptions.overlapDpGapExtend2);
            const auto tProjEnd = steady_clock::now();
            data.threadProjectedAlignmentTime[threadId] += seconds(tProjEnd - tProjStart);
            
            // Error rate filtering.
            if (projectedAlignment.errorRate() > alignOptions.maxErrorRate) {
                data.threadFilteredByErrorRate[threadId]++;
                continue;
            }

            // Create alignment info summary.
            AlignmentInfo alignmentInfo(alignment, 
                uint32_t((*markers).size(orientedReadIds[0].getValue())),
                uint32_t((*markers).size(orientedReadIds[1].getValue()))
            );
            
            alignmentInfo.errorRate = float(projectedAlignment.errorRate());
            alignmentInfo.mismatchCount = uint32_t(projectedAlignment.mismatchCount);
            alignmentInfo.errorRateGaps = float(projectedAlignment.errorRateGaps());
            alignmentInfo.gapCount = uint32_t(projectedAlignment.totalDeletionCount);
            alignmentInfo.gapEventCount = uint32_t(projectedAlignment.totalGapEventCount);
            alignmentInfo.dpScore = projectedAlignment.totalDpScore;

            AlignmentData thisAlignmentData(candidate, alignmentInfo);
            
            // Coordinates were extended to read boundaries during chaining.
            thisAlignmentData.qs = alignment.qs;
            thisAlignmentData.qe = alignment.qe;
            thisAlignmentData.ts = alignment.ts;
            thisAlignmentData.te = alignment.te;
            
            thisAlignmentData.hasLargeIndel = projectedAlignment.hasLargeIndel;
            thisAlignmentData.cisTransStatus = CisTransStatus::Unknown;
            thisAlignmentData.coversHetSite = false;
            thisAlignmentData.deleteReasons0 = AlignmentData::DeleteReasonNone;
            thisAlignmentData.deleteReasons1 = AlignmentData::DeleteReasonNone;

            // --- Populate AlignedEvidenceStore (APES/TASSD) ---
            // Evidence is stored in dual streams (Target-View and Query-View)
            // ensuring Canonical Coordinate Monotonicity.
            thisAlignmentData.info.alignmentId = store.beginAlignment();

            const LongBaseSequenceView tView = reads->getRead(orientedReadIds[1].getReadId());
            const bool tRev = orientedReadIds[1].getStrand();
            DINARA_ASSERT(tView.baseCount <= uint64_t(SnpEvidence::POS_MASK) + 1ULL);
            const uint32_t tRawLen = uint32_t(tView.baseCount);

            static const uint8_t complementBase[4] = {3, 2, 1, 0};

            // Stream 1 (Query-view): positions are in read0 forward coordinates.
            {
                for(const auto& m : projectedAlignment.sparseMismatches) {
                    store.addSnp1(m.position0, m.base1);
                }

                for(const auto& indel : projectedAlignment.sparseIndels) {
                    if(indel.op == 'I') {
                        store.addIndel1(indel.position0, indel.length, 0);
                    } else if(indel.op == 'D') {
                        store.addIndel1(indel.position0, indel.length, 1);
                    } else {
                        DINARA_ASSERT(0);
                    }
                }
            }

            // Stream 0 (Target-view): positions are in read1 forward coordinates.
            {
                if(!tRev) {
                    for(const auto& m : projectedAlignment.sparseMismatches) {
                        store.addSnp0(m.position1, m.base0);
                    }

                    for(const auto& indel : projectedAlignment.sparseIndels) {
                        if(indel.op == 'I') {
                            store.addIndel0(indel.position1, indel.length, 1);
                        } else if(indel.op == 'D') {
                            store.addIndel0(indel.position1, indel.length, 0);
                        } else {
                            DINARA_ASSERT(0);
                        }
                    }

                } else {
                    // Reverse in the alignment: emit in increasing canonical coordinates.
                    for(auto it = projectedAlignment.sparseMismatches.rbegin();
                        it != projectedAlignment.sparseMismatches.rend(); ++it) {

                        const uint32_t posOriented = it->position1;
                        DINARA_ASSERT(posOriented < tRawLen);
                        const uint32_t pos = (tRawLen - 1U) - posOriented;
                        DINARA_ASSERT(it->base0 < 4);
                        // Complement read0 base into read1's forward frame.
                        store.addSnp0(pos, complementBase[it->base0]);
                    }

                    for(auto it = projectedAlignment.sparseIndels.rbegin();
                        it != projectedAlignment.sparseIndels.rend(); ++it) {

                        const uint32_t posOriented = it->position1;
                        DINARA_ASSERT(posOriented < tRawLen);
                        if(it->op == 'I') {
                            const uint32_t pos = tRawLen - (posOriented + it->length);
                            store.addIndel0(pos, it->length, 1);
                        } else if(it->op == 'D') {
                            const uint32_t pos = (tRawLen - 1U) - posOriented;
                            store.addIndel0(pos, it->length, 0);
                        } else {
                            DINARA_ASSERT(0);
                        }
                    }
                }
            }

            threadAlignmentData.push_back(thisAlignmentData);
            
            dinara::compress(alignment, compressedAlignment);
            thisThreadCompressedAlignments.appendVector(
                compressedAlignment.begin(),
                compressedAlignment.end()
            );
        }
    }
    
    thisThreadCompressedAlignments.unreserve();
}



void Assembler::filterSecondaryAlignmentsPerReadPair(
    uint64_t threadCount,
    bool requireNonRedundantOnBothReads)
{
    // --- Hifiasm-style per-read-pair redundancy filtering ---
    // For each read r0, consider all overlaps to a given partner read (r1).
    // Keep the strongest chains and delete redundant/secondary ones:
    // - score ratio filter (hifiasm -p 0.8)
    // - overlap-on-r0 filter (hifiasm -M 0.5)
    // Optionally also require non-redundancy on the partner read interval (ts/te).

    const auto tBegin = steady_clock::now();
    cout << timestamp << "Filtering secondary alignments per read pair begins." << endl;

    removedSecondaryAlignmentCount = 0;
    removedSecondaryAlignmentBySymmetryOnlyCount = 0;
    filterSecondaryRequireNonRedundantOnBothReads = requireNonRedundantOnBothReads;

    // Run threads.
    const uint64_t readCount = reads->readCount();
    setupLoadBalancing(readCount, 100); // Batch size heuristic.
    runThreads(&Assembler::filterSecondaryAlignmentsPerReadPairThreadFunction, threadCount);

    const auto tEnd = steady_clock::now();
    const double tSeconds = seconds(tEnd - tBegin);
    cout << timestamp << "Filtering secondary alignments per read pair ends in " << tSeconds << " s." << endl;
    cout << timestamp << "Removed " << removedSecondaryAlignmentCount << " redundant alignments (marked isDeleted).";
    if(filterSecondaryRequireNonRedundantOnBothReads) {
        cout << " symmetryOnlyRemoved=" << removedSecondaryAlignmentBySymmetryOnlyCount;
    }
    cout << endl;
}



void Assembler::filterSecondaryAlignmentsPerReadPairThreadFunction(size_t)
{
    uint64_t begin, end;
    uint64_t localRemovedCount = 0;
    uint64_t localSymmetryOnlyRemovedCount = 0;

    // Reuse buffers across reads to avoid repeated allocations.
    struct CandidateInfo {
        uint32_t alignmentId;
        int64_t score;
        uint64_t qs;
        uint64_t qe;
        uint64_t ts;
        uint64_t te;
        bool keptOld = false; // kept by legacy (r0-only) logic
        bool keptNew = false; // kept by current logic (with optional symmetry)
    };
    vector<CandidateInfo> group;
    vector<CandidateInfo> kept;
    group.reserve(16);
    kept.reserve(16);

    // Default batch size from setupLoadBalancing.
    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); r0++) {
            const OrientedReadId orientedR0(r0, 0);
            const auto& table = alignmentTable[orientedR0.getValue()];
            if(table.empty()) {
                continue;
            }

            group.clear();
            kept.clear();

            ReadId currentTarget = invalidReadId;

            auto flushGroup = [&]() {
                if(group.empty()) {
                    return;
                }
                if(group.size() == 1) {
                    // A single alignment for this read pair cannot be redundant.
                    return;
                }

                // Sort group by score descending.
                sort(group.begin(), group.end(), [](const CandidateInfo& a, const CandidateInfo& b) {
                    return a.score > b.score;
                });

                const int64_t bestScore = group[0].score;
                auto failsScoreRatio = [&](const CandidateInfo& cand) -> bool {
                    // 1) Score ratio filter (hifiasm -p 0.8), done with integer math:
                    // cand.score < bestScore * 0.8  <=>  cand.score * 5 < bestScore * 4
                    return cand.score * 5 < bestScore * 4;
                };
                auto overlapsMoreThanHalf = [&](uint64_t s0, uint64_t e0, uint64_t s1, uint64_t e1, uint64_t len) -> bool {
                    if(len == 0) {
                        return false;
                    }
                    const uint64_t oStart = std::max(s0, s1);
                    const uint64_t oEnd = std::min(e0, e1);
                    if(oEnd <= oStart) {
                        return false;
                    }
                    const uint64_t oLen = oEnd - oStart;
                    return oLen * 2 > len;
                };
                auto overlapsKeptOnR0 = [&](const CandidateInfo& cand, const vector<CandidateInfo>& keptIntervals) -> bool {
                    const uint64_t len = cand.qe - cand.qs;
                    for(const auto& existing : keptIntervals) {
                        if(overlapsMoreThanHalf(cand.qs, cand.qe, existing.qs, existing.qe, len)) {
                            return true;
                        }
                    }
                    return false;
                };
                auto overlapsKeptOnR1 = [&](const CandidateInfo& cand, const vector<CandidateInfo>& keptIntervals) -> bool {
                    const uint64_t len = cand.te - cand.ts;
                    for(const auto& existing : keptIntervals) {
                        if(overlapsMoreThanHalf(cand.ts, cand.te, existing.ts, existing.te, len)) {
                            return true;
                        }
                    }
                    return false;
                };

                // Fast path: legacy behavior (r0-only redundancy).
                if(!filterSecondaryRequireNonRedundantOnBothReads) {
                    for(auto& cand : group) {
                        cand.keptNew = false;
                    }
                    kept.clear();
                    for(auto& cand : group) {
                        if(failsScoreRatio(cand)) {
                            continue;
                        }
                        if(overlapsKeptOnR0(cand, kept)) {
                            continue;
                        }
                        cand.keptNew = true;
                        kept.push_back(cand);
                    }
                    for(const auto& cand : group) {
                        if(!cand.keptNew) {
                            alignmentData[cand.alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
                            localRemovedCount++;
                        }
                    }
                    return;
                }

                // First pass: legacy behavior (r0-only redundancy).
                for(auto& cand : group) {
                    cand.keptOld = false;
                }
                kept.clear();
                for(auto& cand : group) {
                    if(failsScoreRatio(cand)) {
                        continue;
                    }
                    if(overlapsKeptOnR0(cand, kept)) {
                        continue;
                    }
                    cand.keptOld = true;
                    kept.push_back(cand);
                }

                // Second pass: actual behavior (optionally symmetric).
                for(auto& cand : group) {
                    cand.keptNew = false;
                }
                kept.clear();
                for(auto& cand : group) {
                    if(failsScoreRatio(cand)) {
                        continue;
                    }
                    if(overlapsKeptOnR0(cand, kept)) {
                        continue;
                    }
                    if(filterSecondaryRequireNonRedundantOnBothReads && overlapsKeptOnR1(cand, kept)) {
                        continue;
                    }
                    cand.keptNew = true;
                    kept.push_back(cand);
                }

                // Mark deletions and count symmetry-only removals.
                for(const auto& cand : group) {
                    if(!cand.keptNew) {
                        alignmentData[cand.alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
                        localRemovedCount++;
                        if(cand.keptOld) {
                            localSymmetryOnlyRemovedCount++;
                        }
                    }
                }
            };

            // alignmentTable[OrientedReadId(r0,0)] is already sorted by partner OrientedReadId,
            // so partner ReadIds are contiguous (strand 0/1 are adjacent for the same read).
            // Skip the prefix where the partner oriented read id is < orientedR0 (partner read id < r0).
            // This avoids scanning alignments that would be skipped by canonical ordering checks.
            size_t firstIndex = 0;
            {
                size_t lo = 0;
                size_t hi = table.size();
                while(lo < hi) {
                    const size_t mid = (lo + hi) / 2;
                    const uint32_t alignmentId = table[mid];
                    const OrientedReadId other = alignmentData[alignmentId].getOther(orientedR0);
                    if(other < orientedR0) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                firstIndex = lo;
            }

            for(size_t tableIndex = firstIndex; tableIndex < table.size(); ++tableIndex) {
                const uint32_t alignmentId = table[tableIndex];
                const auto& ad = alignmentData[alignmentId];

                // Mark self-overlaps as deleted immediately.
                if(ad.readIds[0] == ad.readIds[1]) {
                    alignmentData[alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
                    localRemovedCount++;
                    continue;
                }

                // Canonical ordering: only process from the smaller ReadId side.
                // AlignmentData is expected to be stored with readIds[0] < readIds[1].
                // So for this r0, we only act when r0 == readIds[0].
                if(ad.readIds[0] != r0) {
                    continue;
                }
                const ReadId target = ad.readIds[1];

                if(target != currentTarget) {
                    flushGroup();
                    group.clear();
                    currentTarget = target;
                }

                const auto& info = ad.info;

                // Overlap length on r0.
                const uint64_t overlapLenOnR0 = uint64_t(ad.qe) - uint64_t(ad.qs);

                const int64_t score = info.hifiasmDpScoreOrApprox(overlapLenOnR0);

                // Overlap interval on r0 (forward coordinates).
                const uint64_t qs = ad.qs;
                const uint64_t qe = ad.qe;
                const uint64_t ts = ad.ts;
                const uint64_t te = ad.te;

                group.push_back(CandidateInfo{alignmentId, score, qs, qe, ts, te});
            }

            flushGroup();
        }
    }

    removedSecondaryAlignmentCount += localRemovedCount;
    removedSecondaryAlignmentBySymmetryOnlyCount += localSymmetryOnlyRemovedCount;
}
