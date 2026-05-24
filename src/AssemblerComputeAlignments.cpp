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

void Assembler::computeBaseAlignmentsAndStore(
    const AlignOptions& alignOptions,
    uint64_t threadCount
) {
    const auto tBegin = steady_clock::now();
    const size_t candidateCount = alignmentCandidates.candidates.size();

    cout << timestamp << "Begin computing alignments with evidence for ";
    cout << candidateCount << " alignment candidates." << endl;

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
    if(batchSize > candidateCount/threadCount) {
        batchSize = candidateCount/threadCount;
    }
    if(batchSize == 0) {
        batchSize = 1;
    }

    // Prepare data structures for each thread.
    data.threadAlignmentData.resize(threadCount);
    for(auto& v : data.threadAlignmentData) {
        v.clear();
    }
    data.threadCompressedAlignments.resize(threadCount);
    data.threadEvidenceStores.resize(threadCount);
    for(auto& store : data.threadEvidenceStores) {
        store.clear();
    }
    data.threadCigarStores.resize(threadCount);
    for(auto& store : data.threadCigarStores) {
        store.clear();
    }
    
    // Always resize these to avoid SIGSEGV during aggregation.
    data.threadProjectedAlignmentTime.assign(threadCount, 0.0);
    data.threadCollectionTime.assign(threadCount, 0.0);
    data.threadFilteredByErrorRate.assign(threadCount, 0);
    data.threadFilteredByErrorRateGap.assign(threadCount, 0);
    data.threadFilteredByGapCount.assign(threadCount, 0);

    // If needed for variant clustering, resize this as well.
    if (assemblerInfo->readGraphCreationMethod == 5) {
        data.threadVariantClusteringPositionPairs.resize(threadCount);
    } else {
        data.threadVariantClusteringPositionPairs.clear();
    }
    
    performanceLog << timestamp << "Alignment computation begins (Unified Chaining & Evidence path)." << endl;
    cout << timestamp << "Alignment computation begins (Unified Chaining & Evidence path)." << endl;
    
    setupLoadBalancing(candidateCount, batchSize);
    runThreads(&Assembler::computeBaseAlignmentsAndStoreThreadFunction, threadCount);
    
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
    uint64_t totalCigarTokens = 0;
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
        const OverlapCigarStore& localCigarStore = data.threadCigarStores[threadId];
        totalCigarTokens += localCigarStore.tokenCount();
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

    overlapCigarStore.clear();
    overlapCigarStore.reserve(totalCigarTokens);

    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
        const size_t idShift = alignmentData.size(); // Current global count serves as offset.
        if(!threadAlignmentData.empty()) {
            // Merge thread-local packed CIGAR arena first to get the offset shift.
            const uint32_t cigarOffsetShift = overlapCigarStore.merge(data.threadCigarStores[threadId]);
            data.threadCigarStores[threadId].clear();

            alignmentData.resize(idShift + threadAlignmentData.size());
            std::copy(
                threadAlignmentData.begin(),
                threadAlignmentData.end(),
                alignmentData.begin() + idShift);
            for(size_t i = 0; i < threadAlignmentData.size(); ++i) {
                alignmentData[idShift + i].info.alignmentId += idShift;
                if(alignmentData[idShift + i].info.cigarOffset != uint32_t(-1)) {
                    alignmentData[idShift + i].info.cigarOffset += cigarOffsetShift;
                }
            }
        }

        const auto& threadCompressedAlignments = *data.threadCompressedAlignments[threadId];
        const auto size = threadCompressedAlignments.size();
        for(size_t i=0; i<size; i++) {
            compressedAlignments.appendVector(
                threadCompressedAlignments[i].begin(),
                threadCompressedAlignments[i].end()
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
        
        const size_t localIndexSize = localStore.index.size();
        const size_t globalIndexBegin = alignedEvidenceStore.index.size();
        alignedEvidenceStore.index.resize(globalIndexBegin + localIndexSize);
        for(size_t i = 0; i < localIndexSize; ++i) {
            auto entry = localStore.index[i];
            entry.snpOffset0 += globalSnpOffset0;
            entry.indelOffset0 += globalIndelOffset0;
            entry.snpCheckpointOffset0 += globalSnpCheckpointOffset0;
            entry.snpOffset1 += globalSnpOffset1;
            entry.indelOffset1 += globalIndelOffset1;
            entry.snpCheckpointOffset1 += globalSnpCheckpointOffset1;
            alignedEvidenceStore.index[globalIndexBegin + i] = entry;
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

    // Report memory usage for both storage mechanisms.
    {
        const size_t cigarTotalBytes = overlapCigarStore.memoryUsage();

        const size_t evidSnp0 = alignedEvidenceStore.snpStream0.size() * sizeof(SnpEvidence);
        const size_t evidSnp1 = alignedEvidenceStore.snpStream1.size() * sizeof(SnpEvidence);
        const size_t evidIndel0 = alignedEvidenceStore.indelStream0.size() * sizeof(IndelEvidence);
        const size_t evidIndel1 = alignedEvidenceStore.indelStream1.size() * sizeof(IndelEvidence);
        const size_t evidCkpt0 = alignedEvidenceStore.snpCheckpoints0.size() * sizeof(SnpCheckpoint);
        const size_t evidCkpt1 = alignedEvidenceStore.snpCheckpoints1.size() * sizeof(SnpCheckpoint);
        const size_t evidIdx = alignedEvidenceStore.index.size() * sizeof(AlignedEvidenceStore::IndexEntry);
        const size_t evidTotalBytes = evidSnp0 + evidSnp1 + evidIndel0 + evidIndel1 + evidCkpt0 + evidCkpt1 + evidIdx;

        auto kb = [](size_t b) { return b / 1024.0; };
        auto mb = [](size_t b) { return b / (1024.0 * 1024.0); };

        cout << timestamp << "AlignedEvidenceStore: "
             << alignedEvidenceStore.index.size() << " alignments, "
             << mb(evidTotalBytes) << " MB total"
             << " (snp0=" << kb(evidSnp0) << " KB"
             << ", snp1=" << kb(evidSnp1) << " KB"
             << ", indel0=" << kb(evidIndel0) << " KB"
             << ", indel1=" << kb(evidIndel1) << " KB"
             << ", ckpt0=" << kb(evidCkpt0) << " KB"
             << ", ckpt1=" << kb(evidCkpt1) << " KB"
             << ", index=" << kb(evidIdx) << " KB)." << endl;

        cout << timestamp << "OverlapCigarStore: "
             << overlapCigarStore.tokenCount() << " tokens, "
             << mb(cigarTotalBytes) << " MB." << endl;

        cout << timestamp << "Memory ratio (CigarStore / EvidenceStore): "
             << (evidTotalBytes > 0 ? double(cigarTotalBytes) / double(evidTotalBytes) : 0.0)
             << "x." << endl;
    }

    const auto tEnd = steady_clock::now();
    const double elapsedSeconds = seconds(tEnd - tBegin);
    performanceLog << timestamp << "Done computing alignments. Elapsed time: " << elapsedSeconds << " s." << endl;
    cout << timestamp << "Done computing alignments. Elapsed time: " << elapsedSeconds << " s." << endl;
}



void Assembler::computeAlignmentDataFromChainedCandidatesOnly(
    const AlignOptions& alignOptions,
    uint64_t threadCount)
{
    const auto tBegin = steady_clock::now();
    const size_t candidateCount = alignmentCandidates.candidates.size();

    cout << timestamp << "Begin lightweight alignmentData materialization for "
         << candidateCount << " chained candidates." << endl;

    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    (void)threadCount; // This path is intentionally simple; no base DP is run.

    const auto& candidates = alignmentCandidates.candidates;
    const auto& precomputedAlignments = alignmentCandidatesAlignmentsData.alignments;
    const auto& precomputedSharedSeedScores = alignmentCandidatesAlignmentsData.sharedSeedScores;
    const size_t minAlignedMarkerCount = (alignOptions.minAlignedMarkerCount > 0) ?
        size_t(alignOptions.minAlignedMarkerCount) : 0;

    alignmentData.createNew(
        largeDataName("AlignmentData"),
        largeDataPageSize,
        0,
        candidateCount);
    compressedAlignments.createNew(largeDataName("CompressedAlignments"), largeDataPageSize);
    alignedEvidenceStore.clear();

    string compressedAlignment;
    uint64_t skippedEmpty = 0;
    uint64_t skippedShort = 0;

    for(uint64_t candidateIndex=0; candidateIndex<candidateCount; candidateIndex++) {
        const Alignment& alignment = precomputedAlignments[candidateIndex];
        if(alignment.ordinals.empty()) {
            ++skippedEmpty;
            continue;
        }
        if(minAlignedMarkerCount > 0 &&
            alignment.ordinals.size() < minAlignedMarkerCount) {
            ++skippedShort;
            continue;
        }

        const OrientedReadPair& candidate = candidates[candidateIndex];
        const array<OrientedReadId, 2> orientedReadIds = {
            OrientedReadId(candidate.readIds[0], 0),
            OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1)
        };
        const array<span<const CompressedMarker>, 2> markerSpans = {
            (*markers)[orientedReadIds[0].getValue()],
            (*markers)[orientedReadIds[1].getValue()]
        };

        AlignmentInfo alignmentInfo(
            alignment,
            uint32_t(markerSpans[0].size()),
            uint32_t(markerSpans[1].size()));
        if(candidateIndex < precomputedSharedSeedScores.size()) {
            alignmentInfo.sharedSeedScore = precomputedSharedSeedScores[candidateIndex];
        }

        AlignmentData thisAlignmentData(candidate, alignmentInfo);
        thisAlignmentData.qs = alignment.qs;
        thisAlignmentData.qe = alignment.qe;
        thisAlignmentData.ts = alignment.ts;
        thisAlignmentData.te = alignment.te;
        thisAlignmentData.hasLargeIndel = false;
        thisAlignmentData.informativeHetSiteCount0 = 0;
        thisAlignmentData.informativeHetSiteCount1 = 0;
        thisAlignmentData.informativeHetSiteScore = 0;
        thisAlignmentData.deleteReasons0 = AlignmentData::DeleteReasonNone;
        thisAlignmentData.deleteReasons1 = AlignmentData::DeleteReasonNone;

        alignmentData.push_back(thisAlignmentData);

        dinara::compress(alignment, compressedAlignment);
        compressedAlignments.appendVector(
            compressedAlignment.begin(),
            compressedAlignment.end());
    }

    alignmentData.unreserve();
    compressedAlignments.unreserve();

    performanceLog << timestamp << "Creating alignment table." << endl;
    computeAlignmentTable();

    const auto tEnd = steady_clock::now();
    const double elapsedSeconds = seconds(tEnd - tBegin);
    cout << timestamp << "Done lightweight alignmentData materialization. "
         << "kept=" << alignmentData.size()
         << " skippedEmpty=" << skippedEmpty
         << " skippedShort=" << skippedShort
         << " elapsed=" << elapsedSeconds << " s." << endl;
}



void Assembler::computeBaseAlignmentsAndStoreThreadFunction(size_t threadId) {
    auto& data = computeAlignmentsData;
    const AlignOptions& alignOptions = *data.alignOptions;
    auto& threadAlignmentData = data.threadAlignmentData[threadId];
    const auto& candidates = alignmentCandidates.candidates;
    const auto& precomputedAlignments = alignmentCandidatesAlignmentsData.alignments;
    const auto& precomputedSharedSeedScores = alignmentCandidatesAlignmentsData.sharedSeedScores;
    const uint32_t markerK = uint32_t(assemblerInfo->k);
    const bool collectProjectedTiming = (assemblerInfo->readGraphCreationMethod == 5);
    const size_t minAlignedMarkerCount = (alignOptions.minAlignedMarkerCount > 0) ?
        size_t(alignOptions.minAlignedMarkerCount) : 0;
    const double maxErrorRate = alignOptions.maxErrorRate;
    const int64_t dpMatchScore = alignOptions.overlapDpMatchScore;
    const int64_t dpMismatchScore = alignOptions.overlapDpMismatchScore;
    const int64_t dpGapOpen1 = alignOptions.overlapDpGapOpen1;
    const int64_t dpGapExtend1 = alignOptions.overlapDpGapExtend1;
    
    // Initialize compressed alignment storage for this thread.
    data.threadCompressedAlignments[threadId] = make_shared<MemoryMapped::VectorOfVectors<char, uint64_t>>();
    auto& thisThreadCompressedAlignments = *data.threadCompressedAlignments[threadId];
    thisThreadCompressedAlignments.createNew(
        largeDataName("Thread" + to_string(threadId) + "-CompressedAlignments"),
        largeDataPageSize);

    AlignedEvidenceStore& store = data.threadEvidenceStores[threadId];
    OverlapCigarStore& cigarStore = data.threadCigarStores[threadId];
    string compressedAlignment;
    array<OrientedReadId, 2> orientedReadIds;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t candidateIndex = begin; candidateIndex != end; candidateIndex++) {
            const OrientedReadPair& candidate = candidates[candidateIndex];
            
            // This refactored flow EXCLUSIVELY uses precomputed chains.
            // Chaining is now performed upfront during candidate generation/PAF import.
            const Alignment& alignment = precomputedAlignments[candidateIndex];
            if(alignment.ordinals.empty()) {
                continue;
            }
            // Skip low-support candidates early.
            // This avoids spending time in projected alignment construction for pairs
            // that cannot possibly meet Align.minAlignedMarkerCount.
            if(minAlignedMarkerCount > 0 &&
                alignment.ordinals.size() < minAlignedMarkerCount) {
                continue;
            }
            orientedReadIds[0] = OrientedReadId(candidate.readIds[0], 0);
            orientedReadIds[1] = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);
            const array<LongBaseSequenceView, 2> sequenceViews = {
                reads->getRead(orientedReadIds[0].getReadId()),
                reads->getRead(orientedReadIds[1].getReadId())
            };
            const array<span<const CompressedMarker>, 2> markerSpans = {
                (*markers)[orientedReadIds[0].getValue()],
                (*markers)[orientedReadIds[1].getValue()]
            };

            // Compute projected alignment metrics and sparse diffs (mismatches/indels).
            steady_clock::time_point tProjStart;
            if(collectProjectedTiming) {
                tProjStart = steady_clock::now();
            }
            const ProjectedAlignment projectedAlignment(
                markerK,
                orientedReadIds,
                sequenceViews,
                alignment,
                markerSpans,
                ProjectedAlignment::Method::QuickRawSparse,
                dpMatchScore,
                dpMismatchScore,
                dpGapOpen1,
                dpGapExtend1,
                &cigarStore);
            if(collectProjectedTiming) {
                data.threadProjectedAlignmentTime[threadId] += seconds(steady_clock::now() - tProjStart);
            }
            
            // Error rate filtering.
            const double projectedErrorRate = projectedAlignment.errorRate();
            if(projectedErrorRate > maxErrorRate) {
                data.threadFilteredByErrorRate[threadId]++;
                continue;
            }

            // Create alignment info summary.
            AlignmentInfo alignmentInfo(alignment, 
                uint32_t(markerSpans[0].size()),
                uint32_t(markerSpans[1].size())
            );
            if(candidateIndex < precomputedSharedSeedScores.size()) {
                alignmentInfo.sharedSeedScore = precomputedSharedSeedScores[candidateIndex];
            }
            
            alignmentInfo.errorRate = float(projectedErrorRate);
            alignmentInfo.mismatchCount = uint32_t(projectedAlignment.mismatchCount);
            alignmentInfo.nonHomopolymerErrorCount = uint32_t(projectedAlignment.nonHomopolymerErrorCount);
            const double projectedGapErrorRate = projectedAlignment.errorRateGaps();
            alignmentInfo.errorRateGaps = float(projectedGapErrorRate);
            alignmentInfo.gapCount = uint32_t(projectedAlignment.totalIndelBaseCount);
            alignmentInfo.gapEventCount = uint32_t(projectedAlignment.totalGapEventCount);
            alignmentInfo.dpScore = projectedAlignment.totalDpScore;

            AlignmentData thisAlignmentData(candidate, alignmentInfo);
            
            // Store CIGAR boundary positions (marker-based, not extended).
            // read0 is always strand 0 → forward coordinates.
            thisAlignmentData.qs = projectedAlignment.cigarRead0Start;
            thisAlignmentData.qe = projectedAlignment.cigarRead0End;
            // read1 may be strand 1 → convert oriented coords to forward.
            if (candidate.isSameStrand) {
                thisAlignmentData.ts = projectedAlignment.cigarRead1Start;
                thisAlignmentData.te = projectedAlignment.cigarRead1End;
            } else {
                const uint32_t tLen = uint32_t(sequenceViews[1].baseCount);
                thisAlignmentData.ts = tLen - projectedAlignment.cigarRead1End;
                thisAlignmentData.te = tLen - projectedAlignment.cigarRead1Start;
            }
            
            thisAlignmentData.hasLargeIndel = projectedAlignment.hasLargeIndel;
            thisAlignmentData.informativeHetSiteCount0 = 0;
            thisAlignmentData.informativeHetSiteCount1 = 0;
            thisAlignmentData.informativeHetSiteScore = 0;
            thisAlignmentData.deleteReasons0 = AlignmentData::DeleteReasonNone;
            thisAlignmentData.deleteReasons1 = AlignmentData::DeleteReasonNone;

            // --- Populate AlignedEvidenceStore (APES/TASSD) ---
            // Evidence is stored in dual streams (Target-View and Query-View)
            // ensuring Canonical Coordinate Monotonicity.
            const LongBaseSequenceView tView = sequenceViews[1];
            const bool tRev = orientedReadIds[1].getStrand();
            DINARA_ASSERT(tView.baseCount <= uint64_t(SnpEvidence::POS_MASK) + 1ULL);
            const uint32_t tRawLen = uint32_t(tView.baseCount);

            thisAlignmentData.info.alignmentId = store.beginAlignment();
            thisAlignmentData.info.cigarOffset     = projectedAlignment.cigarOffset;
            thisAlignmentData.info.cigarTokenCount = projectedAlignment.cigarTokenCount;

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
                        if(it->op == 'I') {
                            const uint32_t pos = tRawLen - (posOriented + it->length);
                            store.addIndel0(pos, it->length, 1);
                        } else if(it->op == 'D') {
                            // Gap in sequence1 (read1): sparseIndel.position1 is a boundary
                            // in read1 oriented coordinates. Convert boundary b -> len-b in
                            // read1 forward coordinates.
                            DINARA_ASSERT(posOriented <= tRawLen);
                            const uint32_t pos = tRawLen - posOriented;
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
                // Hifiasm tie-break when shared_seed ties: prefer smaller overlapLen.
                sort(group.begin(), group.end(), [](const CandidateInfo& a, const CandidateInfo& b) {
                    if(a.score != b.score) {
                        return a.score > b.score;
                    }
                    const uint64_t lenA = a.qe - a.qs;
                    const uint64_t lenB = b.qe - b.qs;
                    if(lenA != lenB) {
                        return lenA < lenB;
                    }
                    return a.alignmentId < b.alignmentId;
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

                // This filter is intended to run after phasing/EC parity:
                // only consider overlaps that are cis on BOTH reads.
                // (Trans overlaps already carry DeleteReasonPhase and should not compete for "best" here.)
                if(!ad.isCisByBothSides()) {
                    continue;
                }

                // Skip alignments already deleted for other reasons. We only want to mark *additional*
                // redundancies among alignments that are still eligible for the graph.
                if(!ad.keptByBothSides()) {
                    continue;
                }

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

                // Hifiasm keeps the best chain(s) by overlap_region.shared_seed (chaining DP score),
                // not by a base-level alignment DP score. We use a proxy score derived from the
                // marker alignment strength (see AlignmentInfo::hifiasmSharedSeedScoreProxy()).
                (void)overlapLenOnR0;
                const int64_t score = info.hifiasmSharedSeedScoreProxy();

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



void Assembler::keepOnlyBestAlignmentPerReadPairByDpScore(uint64_t threadCount)
{
    // For each read pair (r0,r1), multiple alignment chains can exist.
    // This keeps only the single best-scoring chain (hifiasm shared_seed proxy score),
    // and deletes the rest as redundant/secondary.
    //
    // This is intentionally aggressive. It can remove legitimate alternative placements
    // (repeats / multi-mapping) and should be used as an experimental clean-up step.
    const auto tBegin = steady_clock::now();
    cout << timestamp << "Keeping only best alignment per read pair by DP score begins." << endl;

    removedSecondaryAlignmentCount = 0;
    removedSecondaryAlignmentBySymmetryOnlyCount = 0;

    const uint64_t readCount = reads->readCount();
    setupLoadBalancing(readCount, 100); // Batch size heuristic.
    runThreads(&Assembler::keepOnlyBestAlignmentPerReadPairByDpScoreThreadFunction, threadCount);

    const auto tEnd = steady_clock::now();
    const double tSeconds = seconds(tEnd - tBegin);
    cout << timestamp << "Keeping only best alignment per read pair by DP score ends in " << tSeconds << " s." << endl;
    cout << timestamp << "Removed " << removedSecondaryAlignmentCount << " alignments (marked DeleteReasonSecondary)." << endl;
}



void Assembler::keepOnlyBestAlignmentPerReadPairByDpScoreThreadFunction(size_t)
{
    uint64_t begin = 0, end = 0;
    uint64_t localRemovedCount = 0;

    // Reuse buffers across reads to avoid repeated allocations.
    struct CandidateInfo {
        uint32_t alignmentId = 0;
        int64_t score = 0;
        ReadId target = invalidReadId;
        uint64_t overlapLen = 0;
    };
    vector<CandidateInfo> group;
    group.reserve(16);

    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); r0++) {
            const OrientedReadId orientedR0(r0, 0);
            const auto& table = alignmentTable[orientedR0.getValue()];
            if(table.empty()) {
                continue;
            }

            group.clear();
            ReadId currentTarget = invalidReadId;

            auto flushGroup = [&]() {
                if(group.size() <= 1) {
                    return;
                }
                // Keep the single best by score.
                size_t bestIndex = 0;
                for(size_t i = 1; i < group.size(); ++i) {
                    if(group[i].score > group[bestIndex].score) {
                        bestIndex = i;
                    } else if(group[i].score == group[bestIndex].score && group[i].overlapLen < group[bestIndex].overlapLen) {
                        // Hifiasm tie-break when shared_seed ties: prefer smaller overlapLen.
                        bestIndex = i;
                    }
                }
                for(size_t i = 0; i < group.size(); ++i) {
                    if(i == bestIndex) continue;
                    alignmentData[group[i].alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
                    ++localRemovedCount;
                }
            };

            // Same canonical prefix skip as filterSecondaryAlignmentsPerReadPair to reduce scanning.
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

                // Only consider overlaps that passed phasing (cis on both) and are still eligible.
                if(!ad.isCisByBothSides()) {
                    continue;
                }
                if(!ad.keptByBothSides()) {
                    continue;
                }

                // Mark self-overlaps as deleted immediately.
                if(ad.readIds[0] == ad.readIds[1]) {
                    alignmentData[alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
                    ++localRemovedCount;
                    continue;
                }

                // Canonical ordering: only process from the smaller ReadId side.
                if(ad.readIds[0] != r0) {
                    continue;
                }
                const ReadId target = ad.readIds[1];

                if(target != currentTarget) {
                    flushGroup();
                    group.clear();
                    currentTarget = target;
                }

                const uint64_t overlapLenOnR0 = uint64_t(ad.qe) - uint64_t(ad.qs);
                const int64_t score = ad.info.hifiasmSharedSeedScoreProxy();
                group.push_back(CandidateInfo{alignmentId, score, target, overlapLenOnR0});
            }

            flushGroup();
        }
    }

    removedSecondaryAlignmentCount += localRemovedCount;
}



/**
 * @brief Deduplicate ONT overlap chains to keep one overlap per partner read (hifiasm parity).
 *
 * ## Purpose
 * After lchain+mcopy chaining and base-level alignment, multiple overlap chains may exist
 * for the same (query, partner) read pair. This function implements hifiasm's dedup_chains
 * logic (ecovlp.cpp:2984) to collapse these to a single best overlap per partner.
 *
 * ## Hifiasm Reference: ecovlp.cpp:2984-3030 `dedup_chains(overlap_region_alloc* ol)`
 *
 * ### Hifiasm Algorithm:
 * ```c
 * 1. Sort overlaps by y_id (partner read ID)
 * 2. Group overlaps by y_id
 * 3. For each group with >1 entry:
 *    a. Find best overlap based on:
 *       - Priority 1: Minimize is_match (0=perfect, 1=cis, 2=trans)
 *       - Priority 2: Maximize score = (x_pos_e + 1 - x_pos_s) - 12*non_homopolymer_errors
 *       - Priority 3: Maximize span = x_pos_e + 1 - x_pos_s (tie-breaker)
 *    b. Swap best to position m, increment m
 * 4. Trim ol->length to m (remove non-best in-place)
 * ```
 *
 * ### Dinara Implementation:
 * ```cpp
 * 1. For each query read r0:
 *    a. Process alignments where r0 is readIds[0] (canonical ordering)
 *    b. Group alignments by partner readIds[1]
 *    c. Skip alignments with non-phase deletion reasons
 * 2. For each partner group:
 *    a. Find best alignment based on:
 *       - Priority 1: Minimize isMatchRank (1=cis, 2=trans)
 *       - Priority 2: Maximize score = span - 12*errorCount
 *       - Priority 3: Maximize span (tie-breaker)
 *       - Priority 4: Minimize alignmentId (determinism tie-breaker)
 *    b. Mark non-best as DeleteReasonSecondary (lazy deletion)
 * ```
 *
 * ## Algorithm Equivalence
 *
 * Close to hifiasm's dedup_chains logic, but not byte-for-byte identical.
 *
 * | Aspect | Hifiasm | Dinara | Equivalent? |
 * |--------|---------|--------|-------------|
 * | **Grouping** | By y_id | By readIds[1] (partner) | ✅ Same |
 * | **Priority 1** | Minimize is_match | Minimize stored EC match state | Near-parity |
 * | **Priority 2** | Maximize (span - 12*errors) | Maximize (span - 12*errors) | Approximate |
 * | **Priority 3** | Maximize span | Maximize span | ✅ Same |
 * | **Deletion** | In-place compaction | Mark DeleteReasonSecondary | Behaviorally similar |
 *
 * ### Field Mappings:
 *
 * | Hifiasm Field | Dinara Field | Notes |
 * |---------------|--------------|-------|
 * | `y_id` | `ad.readIds[1]` | Partner read ID |
 * | `x_pos_s, x_pos_e` | `ad.qs, ad.qe` | Query span (when readIds[0] = r0) |
 * | `is_match` | `hifiasmEcMatchState{0,1}` | Stored per read perspective |
 * | `non_homopolymer_errors` | `nonHomopolymerErrorCount` | Stored from ProjectedAlignment when available |
 *
 * ### Current limitation
 *
 * Dinara now stores `nonHomopolymerErrorCount` from ProjectedAlignment's
 * hifiasm-style homopolymer-aware CIGAR accounting. If that field is missing,
 * we fall back to `mismatchCount + gapCount`.
 *
 * ## Why Deduplication is Necessary
 *
 * The lchain+mcopy algorithm (with mcopy_num > 1) can extract multiple high-scoring
 * chains between the same read pair. This is useful during chaining to:
 * 1. Handle repetitive genomic regions (segmental duplications)
 * 2. Capture alternative alignments for ambiguous regions
 * 3. Allow phasing to select the best cis overlap
 *
 * However, after phasing marks cis/trans overlaps, we want only **one** overlap per
 * partner for downstream assembly. This function implements that final deduplication.
 *
 * ## Execution Context
 *
 * This function MUST run after:
 * 1. ✅ `computeBaseAlignmentsAndStore()` - base-level alignment computed
 * 2. ✅ `performHifiasmECParity()` - phasing marks cis/trans (DeleteReasonPhase)
 * 3. ✅ Before final alignment filtering - allows best overlap to survive
 *
 * ## Thread Safety
 *
 * - Parallelized per query read (no conflicts between threads)
 * - Each thread processes disjoint sets of query reads
 * - Atomic counter for removed overlap count
 *
 * @param threadCount Number of threads for parallel processing
 *
 * @complexity O(N*M) where N=read count, M=avg alignments per read
 *             Amortized O(N log M) per read due to binary search for firstIndex
 *
 * @see hifiasm ecovlp.cpp:2984 dedup_chains()
 * @see performHifiasmECParity() for cis/trans phasing
 * @see AlignmentData::DeleteReasonSecondary for lazy deletion marker
 */
void Assembler::deduplicateOntChainsPerPartnerReadHifiasmLike(uint64_t threadCount)
{
    const auto tBegin = steady_clock::now();
    cout << timestamp << "ONT-style dedup (one overlap per partner) begins." << endl;

    removedOntDeduplicatedChainCount = 0;

    const uint64_t readCount = reads->readCount();
    setupLoadBalancing(readCount, 100); // Same heuristic as other per-read filters.
    runThreads(&Assembler::deduplicateOntChainsPerPartnerReadHifiasmLikeThreadFunction, threadCount);

    const auto tEnd = steady_clock::now();
    const double tSeconds = seconds(tEnd - tBegin);
    cout << timestamp << "ONT-style dedup ends in " << tSeconds << " s. Removed "
         << removedOntDeduplicatedChainCount.load() << " redundant overlaps (marked DeleteReasonSecondary)." << endl;
}



/**
 * @brief Worker thread function for ONT chain deduplication per query read.
 *
 * This implements the per-read logic equivalent to hifiasm's dedup_chains.
 * Each thread processes a disjoint subset of query reads.
 *
 * @see deduplicateOntChainsPerPartnerReadHifiasmLike() for algorithm overview
 */
void Assembler::deduplicateOntChainsPerPartnerReadHifiasmLikeThreadFunction(size_t)
{
    uint64_t begin = 0, end = 0;
    uint64_t localRemoved = 0;

    // =========================================================================
    // Per-Alignment Metadata for Group Selection
    // =========================================================================
    struct CandidateInfo {
        uint32_t alignmentId = 0;
        ReadId partner = invalidReadId;
        uint32_t errorCount = 0; // hifiasm-style non-homopolymer error count
        uint32_t span = 0;       // qe - qs (query span)
        uint8_t isMatchRank = 1; // Hifiasm-style: 1=cis, 2=trans
        int64_t score = std::numeric_limits<int64_t>::min(); // span - 12*errorCount
    };
    vector<CandidateInfo> group;
    group.reserve(8);  // Typical case: few overlaps per partner

    // =========================================================================
    // Lambda: Select Best Overlap Per Partner Group (Hifiasm Logic)
    // =========================================================================
    // Implements hifiasm ecovlp.cpp:2994-3014 selection logic:
    //   - Priority 1: Minimize stored hifiasm-style match state
    //   - Priority 2: Maximize score (span - 12*errors)
    //   - Priority 3: Maximize span (tie-breaker)
    //   - Priority 4: Minimize alignmentId (Dinara's determinism tie-breaker)
    auto flushGroup = [&]() {
        if(group.size() <= 1) {
            return;  // No deduplication needed for single overlap
        }

        // -----------------------------------------------------------------
        // Find Best Overlap in Group (Hifiasm Priority Order)
        // -----------------------------------------------------------------
        size_t best = 0;
        for(size_t i = 1; i < group.size(); ++i) {
            const auto& a = group[i];
            const auto& b = group[best];

            // Priority 1: Prefer smaller hifiasm-style match state.
            if(a.isMatchRank != b.isMatchRank) {
                if(a.isMatchRank < b.isMatchRank) best = i;
                continue;
            }

            // Priority 2: Maximize score = span - 12*errors
            // Matches hifiasm: sc = plus - minus; if(sc > mm_sc) sf = 1;
            if(a.score != b.score) {
                if(a.score > b.score) best = i;
                continue;
            }

            // Priority 3: Maximize span (tie-breaker)
            // Matches hifiasm: if((sc == mm_sc) && ((z->x_pos_e+1-z->x_pos_s) > ...)) sf = 1;
            if(a.span != b.span) {
                if(a.span > b.span) best = i;
                continue;
            }

            // Priority 4: Minimize alignmentId (Dinara determinism tie-breaker;
            // hifiasm leaves perfect ties iteration-order dependent).
            if(a.alignmentId < b.alignmentId) {
                best = i;
            }
        }

        // -----------------------------------------------------------------
        // Mark Non-Best as Secondary (Lazy Deletion)
        // -----------------------------------------------------------------
        // Hifiasm swaps best to position m and trims array (in-place deletion).
        // Dinara marks as DeleteReasonSecondary (lazy deletion, removed later).
        for(size_t i = 0; i < group.size(); ++i) {
            if(i == best) continue;  // Keep best overlap
            alignmentData[group[i].alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
            ++localRemoved;
        }
    };

    // =========================================================================
    // Main Loop: Process Each Query Read in Thread's Batch
    // =========================================================================
    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); ++r0) {
            const OrientedReadId orientedR0(r0, 0);
            const auto& table = alignmentTable[orientedR0.getValue()];
            if(table.empty()) {
                continue;  // No alignments for this read
            }

            group.clear();
            ReadId currentPartner = invalidReadId;

            // -----------------------------------------------------------------
            // Binary Search: Skip Prefix Where Partner < r0 (Avoid Duplicates)
            // -----------------------------------------------------------------
            // Dinara stores alignments bidirectionally (both (r0, r1) and (r1, r0)).
            // To avoid processing the same pair twice, we only process pairs where
            // r0 < partner (canonical ordering). Binary search finds the first such entry.
            size_t firstIndex = 0;
            {
                size_t lo = 0;
                size_t hi = table.size();
                while(lo < hi) {
                    const size_t mid = (lo + hi) / 2;
                    const uint32_t alignmentId = table[mid];
                    const OrientedReadId other = alignmentData[alignmentId].getOther(orientedR0);
                    if(other < orientedR0) {
                        lo = mid + 1;  // Partner < r0, skip
                    } else {
                        hi = mid;      // Partner >= r0, could be first
                    }
                }
                firstIndex = lo;
            }

            // -----------------------------------------------------------------
            // Group Alignments by Partner (Hifiasm Groups by y_id)
            // -----------------------------------------------------------------
            for(size_t tableIndex = firstIndex; tableIndex < table.size(); ++tableIndex) {
                const uint32_t alignmentId = table[tableIndex];
                const auto& ad = alignmentData[alignmentId];

                // Canonical ordering: only process groups where r0 is readIds[0]
                // so ad.qs/qe represent span on this query read.
                if(ad.readIds[0] != r0) {
                    continue;
                }
                const ReadId partner = ad.readIds[1];

                // -----------------------------------------------------------------
                // Filter: Skip Overlaps with Non-Phase Deletion Reasons
                // -----------------------------------------------------------------
                // Match hifiasm's ONT behavior: dedup runs after phasing marks
                // cis/trans (is_match), but before overlaps are hard-removed.
                // Therefore we allow phase-deleted (trans) overlaps to participate,
                // but skip overlaps with OTHER deletion reasons (e.g., low quality).
                const auto nonPhase0 = ad.deleteReasons0 & ~AlignmentData::DeleteReasonPhase;
                const auto nonPhase1 = ad.deleteReasons1 & ~AlignmentData::DeleteReasonPhase;
                if(nonPhase0 != AlignmentData::DeleteReasonNone ||
                   nonPhase1 != AlignmentData::DeleteReasonNone) {
                    continue;  // Already filtered, skip
                }

                // -----------------------------------------------------------------
                // New Partner Group: Flush Previous Group
                // -----------------------------------------------------------------
                if(partner != currentPartner) {
                    flushGroup();  // Deduplicate previous partner's group
                    group.clear();
                    currentPartner = partner;
                }

                // -----------------------------------------------------------------
                // Compute Overlap Metrics (Hifiasm Fields)
                // -----------------------------------------------------------------
                const uint32_t span = ad.qe - ad.qs;  // Matches hifiasm: x_pos_e + 1 - x_pos_s

                const uint32_t err =
                    (ad.info.nonHomopolymerErrorCount != invalid<uint32_t>) ?
                    ad.info.nonHomopolymerErrorCount :
                    ((ad.info.mismatchCount == invalid<uint32_t>) ? 0U : ad.info.mismatchCount) +
                    ((ad.info.gapCount == invalid<uint32_t>) ? 0U : ad.info.gapCount);

                // Hifiasm-style EC match state: 1=cis, 2=trans.
                // Lower value (cis) is preferred in dedup.
                const uint8_t isMatchRank = ad.getHifiasmEcMatchStateFromReadPerspective(r0);

                // Score formula: span - 12*errors
                // Matches hifiasm: sc = plus - minus (ecovlp.cpp:2997)
                const int64_t score = int64_t(span) - 12 * int64_t(err);

                group.push_back(CandidateInfo{alignmentId, partner, err, span, isMatchRank, score});
            }

            // Flush final partner group for this query read
            flushGroup();
        }
    }

    // Update global removed count (atomic)
    removedOntDeduplicatedChainCount += localRemoved;
}



// ---------------------------------------------------------------------------
// Pre-phasing chain deduplication (hifiasm dedup_chains port)
// ---------------------------------------------------------------------------
// For each (readIds[0], readIds[1]) pair, keep only the single best chain.
// Runs before marker graph construction, so no phasing info is available.
//
// Scoring (matches hifiasm ecovlp.cpp:2984):
//   1. Maximize score = span - 12 * nonHomopolymerErrors
//   2. Maximize span (tiebreaker)
//   3. Minimize alignmentId (determinism tiebreaker)
//
// Unlike deduplicateOntChainsPerPartnerReadHifiasmLike, this version:
//   - Does not filter by deletion reasons (nothing is deleted yet)
//   - Does not use is_match / phasing state (phasing hasn't run)
//   - Processes ALL alignments regardless of cis/trans status

void Assembler::dedupChainsPrePhasing(uint64_t threadCount)
{
    const auto tBegin = steady_clock::now();
    cout << timestamp << "Pre-phasing chain dedup (one chain per read pair) begins." << endl;

    removedPrePhasingDedupCount = 0;

    const uint64_t readCount = reads->readCount();
    setupLoadBalancing(readCount, 100);
    runThreads(&Assembler::dedupChainsPrePhasingThreadFunction, threadCount);

    const auto tEnd = steady_clock::now();
    const double tSeconds = seconds(tEnd - tBegin);
    cout << timestamp << "Pre-phasing chain dedup ends in " << tSeconds << " s. Removed "
         << removedPrePhasingDedupCount.load() << " secondary chains." << endl;
}



// Per-read-pair chain dedup.
//
// Dinara's mcopy chaining (hifiasm_lchain_qdp_mcopy_fast) can produce
// multiple chains between the same read pair on the same strand.
// Hifiasm's chain_DP only produces one chain per (query, target, strand),
// so it doesn't need this step.
//
// For each (readA, readB) pair we keep at most one chain per strand:
//   - group[0]: same-strand overlaps
//   - group[1]: opposite-strand overlaps
// This matches hifiasm's separate paf[] / reverse_paf[] storage where
// same-strand and opposite-strand overlaps never compete.
//
// Within each strand group, the chain with the highest sharedSeedScore
// (hifiasm's shared_seed = DP chain score) wins. All others are marked
// DeleteReasonSecondary.
//
// Note: this is separate from the max_n_chain per-read global cap,
// which is already applied during chaining.
void Assembler::dedupChainsPrePhasingThreadFunction(size_t)
{
    uint64_t begin = 0, end = 0;
    uint64_t localRemoved = 0;

    struct CandidateInfo {
        uint32_t alignmentId;
        int32_t sharedSeed;
    };

    vector<CandidateInfo> group[2];
    group[0].reserve(4);
    group[1].reserve(4);

    // Process accumulated groups for the current partner.
    // For each strand bucket, keep only the best chain.
    auto flushGroups = [&]() {
        for(int s = 0; s < 2; s++) {
            if(group[s].size() <= 1) {
                group[s].clear();
                continue;
            }

            size_t best = 0;
            for(size_t i = 1; i < group[s].size(); ++i) {
                if(group[s][i].sharedSeed > group[s][best].sharedSeed) {
                    best = i;
                }
            }

            for(size_t i = 0; i < group[s].size(); ++i) {
                if(i == best) continue;
                alignmentData[group[s][i].alignmentId].addDeleteReasonsBoth(
                    AlignmentData::DeleteReasonSecondary);
                ++localRemoved;
            }
            group[s].clear();
        }
    };

    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); ++r0) {
            const OrientedReadId orientedR0(r0, 0);
            const auto& table = alignmentTable[orientedR0.getValue()];
            if(table.empty()) {
                continue;
            }

            group[0].clear();
            group[1].clear();
            ReadId currentPartner = invalidReadId;

            // The alignment table is sorted by partner. Binary search
            // to skip entries where partner < r0 (those are handled
            // when the other read is the query).
            size_t firstIndex = 0;
            {
                size_t lo = 0;
                size_t hi = table.size();
                while(lo < hi) {
                    const size_t mid = (lo + hi) / 2;
                    const uint32_t alignmentId = table[mid];
                    const OrientedReadId other =
                        alignmentData[alignmentId].getOther(orientedR0);
                    if(other < orientedR0) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                firstIndex = lo;
            }

            for(size_t tableIndex = firstIndex; tableIndex < table.size();
                ++tableIndex) {
                const uint32_t alignmentId = table[tableIndex];
                const auto& ad = alignmentData[alignmentId];

                // Only process each pair once: r0 must be readIds[0].
                if(ad.readIds[0] != r0) {
                    continue;
                }

                // Only dedup cis (1) or unclassified (0) chains.
                // Trans (2) and cisDifferentCopy (3) chains represent
                // different structural relationships and should not
                // compete with cis chains.
                if(ad.hifiasmEcMatchState0 > 1 || ad.hifiasmEcMatchState1 > 1) {
                    continue;
                }

                // Self-overlaps should not exist; remove if found.
                if(ad.readIds[0] == ad.readIds[1]) {
                    alignmentData[alignmentId].addDeleteReasonsBoth(
                        AlignmentData::DeleteReasonSecondary);
                    ++localRemoved;
                    continue;
                }

                const ReadId partner = ad.readIds[1];

                // New partner — flush the previous partner's groups.
                if(partner != currentPartner) {
                    flushGroups();
                    currentPartner = partner;
                }

                const int32_t sharedSeed =
                    (ad.info.sharedSeedScore != invalid<int32_t>) ?
                    ad.info.sharedSeedScore : 0;

                const int strandIdx = ad.isSameStrand ? 0 : 1;
                group[strandIdx].push_back(CandidateInfo{alignmentId, sharedSeed});
            }

            // Flush the last partner's groups.
            flushGroups();
        }
    }

    removedPrePhasingDedupCount += localRemoved;
}
