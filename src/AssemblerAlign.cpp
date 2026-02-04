// PngImage.hpp must be included first because of png issues on Ubuntu 16.04.
#include "PngImage.hpp"

// Dinara.
#include "Assembler.hpp"
#include "hifiasmCoordinateTransforms.hpp"
#include "Alignment.hpp"
#include "AlignmentGraph.hpp"
#include "Align4.hpp"
#include "Align6.hpp"
#include "AssemblerOptions.hpp"
#include "compressAlignment.hpp"
#include "performanceLog.hpp"
#include "ProjectedAlignment.hpp"
#include "AlignedEvidenceStore.hpp"
#include "Reads.hpp"
#include "span.hpp"
#include "timestamp.hpp"
using namespace dinara;

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



// Compute a marker alignment of two oriented reads.
void Assembler::alignOrientedReads(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    size_t maxSkip,     // Maximum ordinal skip allowed.
    size_t maxDrift,    // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency
)
{
    alignOrientedReads(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        maxSkip, maxDrift, maxMarkerFrequency
        );
}
void Assembler::alignOrientedReads(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    size_t maxSkip, // Maximum ordinal skip allowed.
    size_t maxDrift,    // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency
)
{
    reads->checkReadsAreOpen();
    reads->checkReadNamesAreOpen();
    checkMarkersAreOpen();


    // Get the markers sorted by kmerId.
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    vector<MarkerWithOrdinal> markers1SortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId[0]);
    getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId[1]);

    // Call the lower level function.
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    const bool debug = true;
    alignOrientedReads(
        markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);

    // Compute the AlignmentInfo.
    uint32_t leftTrim;
    uint32_t rightTrim;
    tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();
    cout << orientedReadId0 << " has " << reads->getRead(orientedReadId0.getReadId()).baseCount;
    cout << " bases and " << markersSortedByKmerId[0].size() << " markers." << endl;
    cout << orientedReadId1 << " has " << reads->getRead(orientedReadId1.getReadId()).baseCount;
    cout << " bases and " << markersSortedByKmerId[1].size() << " markers." << endl;
    cout << "The alignment has " << alignmentInfo.markerCount;
    cout << " markers. Left trim " << leftTrim;
    cout << " markers, right trim " << rightTrim << " markers." << endl;

#if 0
    // For convenience, also write the two oriented reads.
    ofstream fasta("AlignedOrientedReads.fasta");
    writeOrientedRead(orientedReadId0, fasta);
    writeOrientedRead(orientedReadId1, fasta);
#endif
}



// This lower level version takes as input vectors of
// markers already sorted by kmerId.
void Assembler::alignOrientedReads(
    const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
    size_t maxSkip,             // Maximum ordinal skip allowed.
    size_t maxDrift,            // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency
)
{
    // Compute the alignment.
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    const bool debug = true;
    alignOrientedReads(
        markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);
}



void Assembler::alignOrientedReads(
    const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
    size_t maxSkip,             // Maximum ordinal skip allowed.
    size_t maxDrift,            // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency,
    bool debug,
    AlignmentGraph& graph,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo
)
{
    align(markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);
}


bool Assembler::computeAlignmentParity(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    Alignment& alignment
)
{
    // Parity Alignment Reconstruction
    // Uses actual assembler parameters to reconstruct exact alignment path.
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();

    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId[0]);
    getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId[1]);

    AlignmentGraph graph;
    AlignmentInfo alignmentInfo;
    const bool debug = false;
    
    // Use stored alignment parameters
    size_t maxSkip = assemblerInfo->actualMaxSkip > 0 ? assemblerInfo->actualMaxSkip : 30; // Default fallback
    size_t maxDrift = assemblerInfo->actualMaxDrift > 0 ? assemblerInfo->actualMaxDrift : 15;
    uint32_t maxMarkerFrequency = 100; // Standard

    align(
        markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);
        
    return alignmentInfo.markerCount > 0;
}



// Compute marker alignments of an oriented read with all reads
// for which we have an alignment.
void Assembler::alignOverlappingOrientedReads(
    ReadId readId, Strand strand,
    size_t maxSkip,                 // Maximum ordinal skip allowed.
    size_t maxDrift,                // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency,
    size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
    size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
)
{
    alignOverlappingOrientedReads(
        OrientedReadId(readId, strand),
        maxSkip, maxDrift, maxMarkerFrequency, minAlignedMarkerCount, maxTrim);
}



void Assembler::alignOverlappingOrientedReads(
    OrientedReadId orientedReadId0,
    size_t maxSkip,                 // Maximum ordinal skip allowed.
    size_t maxDrift,                // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency,
    size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
    size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
)
{
    // Check that we have what we need.
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Get the markers for orientedReadId0.
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId[0]);

    // Loop over all alignments involving this oriented read.
    size_t goodAlignmentCount = 0;
    for(const uint64_t i: alignmentTable[orientedReadId0.getValue()]) {
        const AlignmentData& ad = alignmentData[i];

        // Get the other oriented read involved in this alignment.
        const OrientedReadId orientedReadId1 = ad.getOther(orientedReadId0);

        // Get the markers for orientedReadId1.
        getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId[1]);

        // Compute the alignment.
        AlignmentGraph graph;
        Alignment alignment;
        AlignmentInfo alignmentInfo;
        const bool debug = false;
        alignOrientedReads(
            markersSortedByKmerId,
            maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);

        uint32_t leftTrim;
        uint32_t rightTrim;
        tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();

        cout << orientedReadId0 << " " << orientedReadId1 << " " << alignmentInfo.markerCount;
        if(alignmentInfo.markerCount) {
            cout << " " << leftTrim << " " << rightTrim;
            if(alignmentInfo.markerCount >= minAlignedMarkerCount && leftTrim<=maxTrim && rightTrim<=maxTrim) {
                cout << " good";
                ++goodAlignmentCount;
            }
        }
        cout << endl;

    }
    cout << "Found " << goodAlignmentCount << " alignments out of ";
    cout << alignmentTable[orientedReadId0.getValue()].size() << "." << endl;

}



// Minimal struct to store overlap info for transitive inference
struct PafOverlap {
    ReadId otherId;
    bool sameStrand;
    uint32_t start;
    uint32_t end;
};

void Assembler::importAlignmentCandidatesFromPaf(const string& pafFilePath)
{
    if (!std::filesystem::exists(pafFilePath)) {
        throw runtime_error("PAF file not found: " + pafFilePath);
    }

    cout << timestamp << "Loading alignment candidates from " << pafFilePath << " with transitive expansion..." << endl;
    alignmentCandidates.candidates.createNew(largeDataName("AlignmentCandidates"), largeDataPageSize);
    
    // Store overlaps per read for geometric checking: overlaps[readId] -> list of intervals
    vector<vector<PafOverlap>> overlaps(reads->readCount());
    
    // Vectors to store direct candidates for later sorting
    // We'll create thread-local storage for parsing if needed, but parsing is usually I/O bound.
    // Let's do single-threaded parse, but it's fast.
    
    std::ifstream pafFile(pafFilePath);
    string line;
    uint64_t lineCount = 0;
    
    while (std::getline(pafFile, line)) {
        std::stringstream ss(line);
        string qName, tName, strandStr;
        uint64_t qLen, qStart, qEnd, tLen, tStart, tEnd, mapQ, alignLen;
        
        // PAF columns:
        // 1: Query name
        // 2: Query length
        // 3: Query start
        // 4: Query end
        // 5: Strand (+/-)
        // 6: Target name
        // 7: Target length
        // 8: Target start
        // 9: Target end
        // 10: Number of residue matches
        // 11: Alignment block length
        
        if (!(ss >> qName >> qLen >> qStart >> qEnd >> strandStr >> tName >> tLen >> tStart >> tEnd >> mapQ >> alignLen)) {
            continue; 
        }

        if (alignLen < 200) {
            continue;
        }

        try {
            ReadId readId0 = reads->getReadId(qName);
            ReadId readId1 = reads->getReadId(tName);

            if (readId0 == invalid<ReadId> || readId1 == invalid<ReadId>) continue;
            // Ignore palindromic for now
            if (reads->getFlags(readId0).isPalindromic || reads->getFlags(readId1).isPalindromic) continue;
            if (readId0 == readId1) continue; 

            bool isSameStrand = (strandStr == "+");
            
            // Store for transitive inference (Geometric Graph)
            // Q aligns to T.
            // On Q, the interval is [qStart, qEnd].
            // On T, the interval is [tStart, tEnd].
            
            // Add neighbor T to Q
            overlaps[readId0].push_back({readId1, isSameStrand, (uint32_t)qStart, (uint32_t)qEnd});
            // Add neighbor Q to T
            overlaps[readId1].push_back({readId0, isSameStrand, (uint32_t)tStart, (uint32_t)tEnd});
            
            // Store direct candidate
            if (readId0 > readId1) swap(readId0, readId1);
            alignmentCandidates.candidates.push_back(OrientedReadPair(readId0, readId1, isSameStrand));
            
            lineCount++;
            
        } catch (...) {
            continue;
        }
    }
    cout << timestamp << "Parsed " << lineCount << " PAF lines." << endl;

    // Transitive Expansion with Geometric Filter
    cout << timestamp << "Generating transitive candidates (Distance 2 with overlap check)..." << endl;
    
    uint64_t threadCount = std::thread::hardware_concurrency();
    vector<vector<OrientedReadPair>> threadNewCandidates(threadCount);
    
    const uint64_t batchSize = 200;
    setupLoadBalancing(reads->readCount(), batchSize);

    vector<std::thread> threads;
    for (size_t t = 0; t < threadCount; t++) {
        threads.emplace_back([&, t]() {
            uint64_t start, end;
            while (getNextBatch(start, end)) {
                for (uint64_t r = start; r < end; r++) {
                    auto& neighbors = overlaps[r]; // Mutable reference, safe as r is exclusive to thread
                    if (neighbors.size() < 2) continue; 
                    
                    // OPTIMIZATION: Sort by start position to enable sliding window
                    std::sort(neighbors.begin(), neighbors.end(), 
                        [](const PafOverlap& a, const PafOverlap& b) {
                            return a.start < b.start;
                        });
                    
                    // Check overlapping pairs using sliding window
                    for (size_t i = 0; i < neighbors.size(); i++) {
                        const auto& n1 = neighbors[i];
                        
                        // We need overlap > 200.
                        // Implies: min(n1.end, n2.end) - max(n1.start, n2.start) > 200
                        // Since sorted by start, n2.start >= n1.start.
                        // So max(starts) = n2.start.
                        // Overlap = min(n1.end, n2.end) - n2.start > 200
                        // We need min(n1.end, n2.end) > n2.start + 200
                        // Necessary condition: n1.end > n2.start + 200  (since min(...) <= n1.end)
                        // So if n2.start >= n1.end - 200, we can stop, as n2.start will only increase.
                        
                        uint32_t limit = (n1.end > 200) ? n1.end - 200 : 0;

                        for (size_t j = i + 1; j < neighbors.size(); j++) {
                            const auto& n2 = neighbors[j];
                            
                            // Optimization: Early exit
                            if (n2.start >= limit) break; 
                            
                            // Check if distinct neighbors
                            if (n1.otherId == n2.otherId) continue;
                            
                            // Strict overlap calculation
                            uint32_t overlapStart = n2.start; // Since sorted
                            uint32_t overlapEnd = std::min(n1.end, n2.end);
                            
                            if (overlapEnd > overlapStart + 200) {
                                // Overlap exists! Transitive candidate identified.
                                bool impliedSame = (n1.sameStrand == n2.sameStrand);
                                
                                ReadId rA = n1.otherId;
                                ReadId rC = n2.otherId;
                                if (rA > rC) swap(rA, rC);
                                
                                threadNewCandidates[t].push_back(OrientedReadPair(rA, rC, impliedSame));
                            }
                        }
                    }
                }
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }
    
    // Merge new candidates
    for (const auto& vec : threadNewCandidates) {
        for (const auto& c : vec) {
            alignmentCandidates.candidates.push_back(c);
        }
    }

    // Sort and remove duplicates.
    std::sort(alignmentCandidates.candidates.begin(), alignmentCandidates.candidates.end(),
        [](const OrientedReadPair& a, const OrientedReadPair& b) {
            return tie(a.readIds[0], a.readIds[1], a.isSameStrand) <
                   tie(b.readIds[0], b.readIds[1], b.isSameStrand);
        });
    
    auto it = std::unique(alignmentCandidates.candidates.begin(), alignmentCandidates.candidates.end(),
        [](const OrientedReadPair& a, const OrientedReadPair& b) {
            return a.readIds[0] == b.readIds[0] &&
                   a.readIds[1] == b.readIds[1] &&
                   a.isSameStrand == b.isSameStrand;
        });
    
    alignmentCandidates.candidates.resize(it - alignmentCandidates.candidates.begin());

    alignmentCandidates.unreserve();
    cout << timestamp << "Total unique candidates (Direct + Transitive): " << alignmentCandidates.candidates.size() << endl;
}


// Compute alignments.
void Assembler::computeAlignments(

    const AlignOptions& alignOptions,

    // Number of threads. If zero, a number of threads equal to
    // the number of virtual processors is used.
    uint64_t threadCount
)
{

    const auto tBegin = steady_clock::now();


    performanceLog << timestamp << "Begin computing alignments for ";
    performanceLog << alignmentCandidates.candidates.size() << " alignment candidates." << endl;

    // Check that we have what we need.
    reads->checkReadsAreOpen();
    // Note: kmerChecker is not required when using SIMD minimizers for marker generation.
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Store parameters so they are accessible to the threads.
    auto& data = computeAlignmentsData;
    data.alignOptions = &alignOptions;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // For alignment method 4, compute sorted markers.
    if(alignOptions.alignMethod == 4) {
        cout << timestamp << "Computing sorted markers." << endl;
        computeSortedMarkers(threadCount);
    }

    // For alignment method 5, compute low frequency markers.
    if(alignOptions.alignMethod == 5) {
        cout << timestamp << "Computing unique markers." << endl;
        computeLowFrequencyMarkers(1, threadCount);
    }

    // For alignment method 6, compute Align6Markers.
    if(alignOptions.alignMethod == 6) {
        performanceLog << timestamp << "Computing Align6Markers." << endl;
        computeAlign6Markers(threadCount);
        performanceLog << timestamp << "Done computing Align6Markers." << endl;
    }

    // Pick the batch size for computing alignments.
    size_t batchSize = 10;
    if(batchSize > alignmentCandidates.candidates.size()/threadCount) {
        batchSize = alignmentCandidates.candidates.size()/threadCount;
    }
    if(batchSize == 0) {
        batchSize = 1;
    }


    // Compute the alignments.
    data.threadAlignmentData.resize(threadCount);
    data.threadCompressedAlignments.resize(threadCount);
    data.threadEvidenceStores.resize(threadCount);
    
    // Always resize these to avoid SIGSEGV during aggregation, even if not used.
    data.threadProjectedAlignmentTime.assign(threadCount, 0.0);
    data.threadCollectionTime.assign(threadCount, 0.0);
    data.threadFilteredByErrorRate.assign(threadCount, 0);
    data.threadFilteredByErrorRateGap.assign(threadCount, 0);
    data.threadFilteredByGapCount.assign(threadCount, 0);

    if (assemblerInfo->readGraphCreationMethod == 5) {
        data.threadVariantClusteringPositionPairs.resize(threadCount);
    }
    
    performanceLog << timestamp << "Alignment computation begins." << endl;
    cout << timestamp << "Alignment computation begins." << endl;
    setupLoadBalancing(alignmentCandidates.candidates.size(), batchSize);
    runThreads(&Assembler::computeAlignmentsThreadFunction, threadCount);
    performanceLog << timestamp << "Alignment computation completed." << endl;
    cout << timestamp << "Alignment computation completed." << endl;

    // Store the alignments found by each thread.
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
        
        // Append Indexes (Adjusting offsets)
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

        // Clear local store to free memory immediately
        localStore.clear();
    }

    // Release unused allocated memory.
    alignmentData.unreserve();
    compressedAlignments.unreserve();





    // Aggregate timing statistics from all threads
    double totalProjectedAlignmentTime = 0.0;
    double totalCollectionTime = 0.0;
    double estimatedProjectedWallTime = 0.0;
    double estimatedCollectionWallTime = 0.0;
    
    if (assemblerInfo->readGraphCreationMethod == 5) {
        for (size_t i = 0; i < threadCount; i++) {
            totalProjectedAlignmentTime += data.threadProjectedAlignmentTime[i];
            totalCollectionTime += data.threadCollectionTime[i];
        }
    
        // Estimate wall-clock time by dividing summed thread time by thread count
        estimatedProjectedWallTime = totalProjectedAlignmentTime / threadCount;
        estimatedCollectionWallTime = totalCollectionTime / threadCount;
    }

    if (assemblerInfo->readGraphCreationMethod == 5) {
        uint64_t totalFilteredByErrorRate = 0;
        uint64_t totalFilteredByErrorRateGap = 0;
        uint64_t totalFilteredByGapCount = 0;
        for (size_t i = 0; i < threadCount; i++) {
            totalFilteredByErrorRate += data.threadFilteredByErrorRate[i];
            totalFilteredByErrorRateGap += data.threadFilteredByErrorRateGap[i];
            totalFilteredByGapCount += data.threadFilteredByGapCount[i];
        }
        cout << "Time spent in projected alignment construction (all threads): " << totalProjectedAlignmentTime << " s" << endl;
        cout << "Time spent collecting variant clustering position pairs (all threads): " << totalCollectionTime << " s" << endl;
        cout << "Alignments filtered by error rate (> 0.07): " << totalFilteredByErrorRate << endl;
        cout << "Alignments filtered by gap error rate (> 0.006): " << totalFilteredByErrorRateGap << endl;
        cout << "Alignments filtered by gap count (> 64): " << totalFilteredByGapCount << endl;
    }
    if (assemblerInfo->readGraphCreationMethod == 5) {
#if DINARA_ENABLE_VARIANT_CLUSTERING
        variantClusteringProjectedAlignmentTime = estimatedProjectedWallTime;
        variantClusteringCollectionTime = estimatedCollectionWallTime;

        cout << "\nVariant clustering collection timing (estimated wall-clock):" << endl;
        cout << "  Projected alignment: ~" << estimatedProjectedWallTime << " s" << endl;
        cout << "  Position pair collection: ~" << estimatedCollectionWallTime << " s" << endl;
        performanceLog << timestamp << "Projected alignment (estimated wall-clock): " << estimatedProjectedWallTime << " s" << endl;
        performanceLog << timestamp << "Collection (estimated wall-clock): " << estimatedCollectionWallTime << " s" << endl;

        // Store position pairs collected by each thread
        performanceLog << timestamp << "Storing position pairs for variant clustering." << endl;
        cout << timestamp << "Storing position pairs for variant clustering." << endl;
        storeVariantClusteringPositionPairs(threadCount, data);
        performanceLog << timestamp << "Done storing position pairs for variant clustering." << endl;
        cout << timestamp << "Done storing position pairs for variant clustering." << endl;
#else
        throw runtime_error("readGraphCreationMethod=5 requires DINARA_ENABLE_VARIANT_CLUSTERING=ON.");
#endif
    }



    // Cleanup.
    if(alignOptions.alignMethod == 4) {
        sortedMarkers->remove();
    }
    if(alignOptions.alignMethod == 5) {
        lowFrequencyMarkers->remove();
    }
    if(alignOptions.alignMethod == 6) {
        align6Markers->remove();
    }

    cout << "Found and stored " << alignmentData.size() << " good alignments." << endl;
    performanceLog << timestamp << "Creating alignment table." << endl;
    computeAlignmentTable();

    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    performanceLog << timestamp << "Computation of alignments ";
    performanceLog << "completed in " << tTotal << " s." << endl;

    performanceLog << timestamp;
}



void Assembler::computeAlignmentsThreadFunction(size_t threadId)
{

    array<OrientedReadId, 2> orientedReadIds;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    string compressedAlignment;

    const bool debug = false;
    auto& data = computeAlignmentsData;
    const size_t alignmentMethod = data.alignOptions->alignMethod;
    const uint32_t maxMarkerFrequency = data.alignOptions->maxMarkerFrequency;
    const size_t maxSkip = data.alignOptions->maxSkip;
    const size_t maxDrift = data.alignOptions->maxDrift;
    const size_t minAlignedMarkerCount = data.alignOptions->minAlignedMarkerCount;
    const double minAlignedFraction = data.alignOptions->minAlignedFraction;
    const size_t maxTrim = data.alignOptions->maxTrim;
    const int matchScore = data.alignOptions->matchScore;
    const int mismatchScore = data.alignOptions->mismatchScore;
    const int gapScore = data.alignOptions->gapScore;
    const double downsamplingFactor = data.alignOptions->downsamplingFactor;
    const int bandExtend = data.alignOptions->bandExtend;
    const int maxBand = data.alignOptions->maxBand;
    const bool suppressContainments = data.alignOptions->suppressContainments;
    const double align5DriftRateTolerance = data.alignOptions->align5DriftRateTolerance;
    const uint64_t align5MinBandExtend = data.alignOptions->align5MinBandExtend;
    const bool computeProjectedAlignmentMetrics = true;


    // Align4-specific items.
    Align4::Options align4Options;
    MemoryMapped::ByteAllocator byteAllocator;
    if(alignmentMethod == 4) {
        align4Options.deltaX = data.alignOptions->align4DeltaX;
        align4Options.deltaY = data.alignOptions->align4DeltaY;
        align4Options.minEntryCountPerCell = data.alignOptions->align4MinEntryCountPerCell;
        align4Options.maxDistanceFromBoundary = data.alignOptions->align4MaxDistanceFromBoundary;
        align4Options.minAlignedMarkerCount = minAlignedMarkerCount;
        align4Options.minAlignedFraction = minAlignedFraction;
        align4Options.maxSkip = maxSkip;
        align4Options.maxDrift = maxDrift;
        align4Options.maxTrim = maxTrim;
        align4Options.maxBand = maxBand;
        align4Options.matchScore = matchScore;
        align4Options.mismatchScore = mismatchScore;
        align4Options.gapScore = gapScore;
        byteAllocator.createNew(
            largeDataName("tmp-ByteAllocator-" + to_string(threadId)),
            largeDataPageSize, 2ULL * 1024 * 1024 * 1024);
    }

    ofstream nullStream;
    Align6 align6(
        assemblerInfo->k,
        maxSkip,
        maxDrift,
        data.alignOptions->align6Options,
        assemblerInfo->kmerDistributionInfo,
        nullStream);
    if((threadId == 0) and (alignmentMethod == 6)) {
        std::lock_guard<std::mutex> lock(mutex);
        align6.writeGlobalFrequencyCriteria(cout);
    }

    vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
    
    shared_ptr< MemoryMapped::VectorOfVectors<char, uint64_t> > thisThreadCompressedAlignmentsPointer =
        make_shared< MemoryMapped::VectorOfVectors<char, uint64_t> >();
    data.threadCompressedAlignments[threadId] = thisThreadCompressedAlignmentsPointer;
    auto& thisThreadCompressedAlignments = *thisThreadCompressedAlignmentsPointer;
    thisThreadCompressedAlignments.createNew(
        largeDataName("tmp-ThreadGlobalCompressedAlignments-" + to_string(threadId)),
        largeDataPageSize);

    if (assemblerInfo->readGraphCreationMethod == 5) {
        // Create the vector to contain the position pairs collected by this thread.
        shared_ptr< MemoryMapped::Vector< pair<OrientedReadId, uint32_t> > > thisThreadVariantClusteringPositionPairsPointer =
            make_shared< MemoryMapped::Vector< pair<OrientedReadId, uint32_t> > >();
        data.threadVariantClusteringPositionPairs[threadId] = thisThreadVariantClusteringPositionPairsPointer;
        auto& thisThreadVariantClusteringPositionPairs = *thisThreadVariantClusteringPositionPairsPointer;
        thisThreadVariantClusteringPositionPairs.createNew(
            largeDataName("tmp-ThreadVariantClusteringPositionPairs-" + to_string(threadId)),
            largeDataPageSize);
    }



#if 0
    // A vector to store the time taken to compute each alignment.
    vector< pair<uint64_t, double> > elapsedTime;
#endif

    const uint64_t messageFrequency = max(1UL, min(1000000UL, alignmentCandidates.candidates.size()/20));

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        for(size_t i=begin; i!=end; i++) {
#if 0
            const auto steadyClock0 = std::chrono::steady_clock::now();
#endif

            const OrientedReadPair& candidate = alignmentCandidates.candidates[i];
            DINARA_ASSERT(candidate.readIds[0] < candidate.readIds[1]);

            // Get the oriented read ids, with the first one on strand 0.
            orientedReadIds[0] = OrientedReadId(candidate.readIds[0], 0);
            orientedReadIds[1] = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);

            if((i % messageFrequency) == 0){
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << "Working on alignment " << i;
                performanceLog << " of " << alignmentCandidates.candidates.size();
                // performanceLog << ": " << orientedReadIds[0] << " " << orientedReadIds[1];
                performanceLog << endl;
            }


            // Compute the alignment.
            bool precomputedUsed = false;
            try {
                if(alignmentMethod == 0) {
                    // ...
                }
                
                // Direct Chain Propagation Check
                if(!alignmentCandidatesAlignmentsData.alignments.empty() && 
                   (alignmentMethod == 5 || alignmentMethod == 6)) { // Only supported for these methods/InvertedIndex
                    if(i < alignmentCandidatesAlignmentsData.alignments.size()) {
                        // Use precomputed alignment.
                        alignment = alignmentCandidatesAlignmentsData.alignments[i];
                        uint32_t mCount0 = uint32_t((*markers)[orientedReadIds[0].getValue()].size());
                        uint32_t mCount1 = uint32_t((*markers)[orientedReadIds[1].getValue()].size());
                        alignmentInfo.create(alignment, mCount0, mCount1);
                        precomputedUsed = true;
                    }
                }

                if(precomputedUsed) {
                    // Skip alignment computation.
                } else if(alignmentMethod == 0) {

                    // Get the markers for the two oriented reads in this candidate.
                    for(size_t j=0; j<2; j++) {
                        getMarkersSortedByKmerId(orientedReadIds[j], markersSortedByKmerId[j]);
                    }

                    // Compute the Alignment.
                    alignOrientedReads(
                        markersSortedByKmerId,
                        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);

                } else if(alignmentMethod == 1) {
                    alignOrientedReads1(orientedReadIds[0], orientedReadIds[1],
                        matchScore, mismatchScore, gapScore,
                        alignment, alignmentInfo);
                } else if(alignmentMethod == 3) {
                    alignOrientedReads3(orientedReadIds[0], orientedReadIds[1],
                        matchScore, mismatchScore, gapScore,
                        downsamplingFactor, bandExtend, maxBand,
                        alignment, alignmentInfo);
                } else if(alignmentMethod == 4) {
                    alignOrientedReads4(orientedReadIds[0], orientedReadIds[1],
                        align4Options,
                        byteAllocator,
                        alignment, alignmentInfo,
                        false);
                    DINARA_ASSERT(byteAllocator.isEmpty());
                } else if(alignmentMethod == 5) {
                    ofstream nullStream;
                    alignOrientedReads5(orientedReadIds[0], orientedReadIds[1],
                        matchScore, mismatchScore, gapScore,
                        align5DriftRateTolerance, align5MinBandExtend,
                        alignment, alignmentInfo,
                        nullStream);
                } else if(alignmentMethod == 6) {
                    ofstream nullStream;
                    alignOrientedReads6(orientedReadIds[0], orientedReadIds[1],
                        alignment, alignmentInfo, align6);
                } else {
                    DINARA_ASSERT(0);
                }
            } catch (std::exception& e) {
                std::lock_guard<std::mutex> lock(mutex);
                cout <<
                    "An error occurred while computing a marker alignment "
                    " of oriented reads " << orientedReadIds[0] << " and " << orientedReadIds[1] <<
                    ". This alignment candidate will be skipped. Error description is: " <<
                    e.what() << endl;
                continue;
            } catch(...) {
                std::lock_guard<std::mutex> lock(mutex);
                cout <<
                    "An error occurred while computing a marker alignment "
                    " of oriented reads " << orientedReadIds[0] << " and " << orientedReadIds[1] <<
                    ". This alignment candidate will be skipped. " << endl;
                continue;
            }

#if 0
            const auto steadyClock1 = std::chrono::steady_clock::now();
            const auto deltaT = 1.e-9 * double(std::chrono::duration_cast<std::chrono::nanoseconds>(steadyClock1 - steadyClock0).count());
            elapsedTime.push_back({i, deltaT});
#endif

            // If the alignment has too few markers, skip it.
            if(alignment.ordinals.size() < minAlignedMarkerCount) {
                continue;
            }

            // If the aligned fraction is too small, skip it.
            if(min(alignmentInfo.alignedFraction(0), alignmentInfo.alignedFraction(1)) < minAlignedFraction) {
                continue;
            }

            // If the alignment has too much trim, skip it.
            uint32_t leftTrim;
            uint32_t rightTrim;
            tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();
            if(leftTrim>maxTrim || rightTrim>maxTrim) {
                continue;
            }

            // For alignment methods other than method 0, we also need to check for
            // maxSip and maxDrift. Method 0 does that automatically.
            if(alignmentMethod != 0) {
                if(alignment.maxSkip() > maxSkip) {
                    continue;
                }
                if(alignment.maxDrift() > maxDrift) {
                    continue;
                }
            }

            // If getting here, this is a good alignment.

            // Compute projected alignment metrics.
            bool hasLargeIndel = false;
            const auto tProjStart = steady_clock::now();
            const ProjectedAlignment projectedAlignment(
                *this,
                orientedReadIds,
                alignment,
                ProjectedAlignment::Method::QuickRawSparse,
                data.alignOptions->overlapDpMatchScore,
                data.alignOptions->overlapDpMismatchScore,
                data.alignOptions->overlapDpGapOpen1,
                data.alignOptions->overlapDpGapExtend1,
                data.alignOptions->overlapDpGapOpen2,
                data.alignOptions->overlapDpGapExtend2);
            const auto tProjEnd = steady_clock::now();
            
            alignmentInfo.errorRate = float(projectedAlignment.errorRate());
            alignmentInfo.mismatchCount = uint32_t(projectedAlignment.mismatchCount);
            alignmentInfo.errorRateGaps = float(projectedAlignment.errorRateGaps());
            alignmentInfo.gapCount = uint32_t(projectedAlignment.totalDeletionCount);
            alignmentInfo.gapEventCount = uint32_t(projectedAlignment.totalGapEventCount); // Transfer gap events
            alignmentInfo.dpScore = projectedAlignment.totalDpScore;
            

            data.threadProjectedAlignmentTime[threadId] += seconds(tProjEnd - tProjStart);

            const uint64_t ql = projectedAlignment.totalLength[0];
            const uint64_t tl = projectedAlignment.totalLength[1];
            
            const double errorRateThreshold = 0.07;
            const uint64_t totalErrors = uint64_t(projectedAlignment.mismatchCount) + uint64_t(projectedAlignment.totalDeletionCount);
            if ((totalErrors > (ql * errorRateThreshold)) || (totalErrors > (tl * errorRateThreshold))) {
                data.threadFilteredByErrorRate[threadId]++;
                continue;
            }

            // const double gapRateThreshold = 0.006;
            // const uint64_t totalGapCount = uint64_t(projectedAlignment.totalDeletionCount);
            // if ((totalGapCount > (ql * gapRateThreshold)) || (totalGapCount > (tl * gapRateThreshold))) {
            //     data.threadFilteredByErrorRateGap[threadId]++;
            //     continue;
            // }

            // // Skip alignments with any single indel >= 64 bases.
            // if (projectedAlignment.maxIndelSize > 64) {
            //     data.threadFilteredByGapCount[threadId]++;
            //     continue;
            // }

            // // Collect position pairs for variant clustering
            // // Only collect those with SNP differences (No indels)
            // const auto tCollectStart = steady_clock::now();
            // collectVariantClusteringPositionPairs(
            //     projectedAlignment,
            //     orientedReadIds,
            //     *data.threadVariantClusteringPositionPairs[threadId]
            // );
            // const auto tCollectEnd = steady_clock::now();
            // data.threadCollectionTime[threadId] += seconds(tCollectEnd - tCollectStart);

            // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " good." << endl;
            AlignmentData thisAlignmentData(candidate, alignmentInfo);
            
            // Calculate coordinates efficiently using direct marker access
            if (precomputedUsed) {
                // Option A: Use stored coordinates from Alignment (computed during InvertedIndex chaining)
                // This avoids re-fetching markers and re-computing extension logic.
                // Parity logic (including -1 for Diff) was already applied in InvertedIndexFinder.
                thisAlignmentData.qs = alignment.qs;
                thisAlignmentData.qe = alignment.qe;
                thisAlignmentData.ts = alignment.ts;
                thisAlignmentData.te = alignment.te;
            } else {
                // Fallback: Compute coordinates from markers (Legacy path for other aligners)
                const auto markers0 = (*markers)[orientedReadIds[0].getValue()];
                uint32_t qs_marker = markers0[alignmentInfo.data[0].firstOrdinal].position;
                uint32_t qe_marker = markers0[alignmentInfo.data[0].lastOrdinal].position + uint32_t(assemblerInfo->k);
    
                const auto markers1 = (*markers)[orientedReadIds[1].getValue()];
                uint32_t ts_marker = markers1[alignmentInfo.data[1].firstOrdinal].position;
                uint32_t te_marker = markers1[alignmentInfo.data[1].lastOrdinal].position + uint32_t(assemblerInfo->k);
    
                // Read lengths for coordinate extension
                const uint64_t len0 = reads->getReadRawSequenceLength(orientedReadIds[0].getReadId());
                const uint64_t len1 = reads->getReadRawSequenceLength(orientedReadIds[1].getReadId());
                
                // --- Hifiasm Parity: Extend coordinates to read boundaries ---
                // This matches append_inexact_overlap_region_alloc in Hash_Table.cpp:374-398
                // The logic extends alignment coordinates so they reach read tips symmetrically.
                
                uint32_t qs_ext = qs_marker;
                uint32_t qe_ext = qe_marker;
                uint32_t ts_ext = ts_marker;
                uint32_t te_ext = te_marker;
                
                // Extend start: whoever has shorter overhang extends to 0, other adjusts proportionally
                if (qs_ext <= ts_ext) {
                    ts_ext = ts_ext - qs_ext;  // Reduce ts by qs amount
                    qs_ext = 0;                 // qs goes to 0
                } else {
                    qs_ext = qs_ext - ts_ext;  // Reduce qs by ts amount
                    ts_ext = 0;                 // ts goes to 0
                }
                
                // Extend end: whoever has shorter remaining distance extends to read length
                // Hifiasm uses closed coordinates: x_right_length = xLen - x_pos_e - 1
                // Dinara uses half-open [start, end): equivalent is len0 - qe_ext
                // But we need to be careful about the semantics:
                // Hifiasm's x_pos_e is inclusive (last base), Dinara's qe_ext is exclusive (one past last)
                // So: Hifiasm's "remaining" = xLen - x_pos_e - 1 = positions after x_pos_e
                // Dinara's "remaining" = len0 - qe_ext = same thing (qe_ext points one past last aligned base)
                
                // For closed-to-half-open conversion: closed_end + 1 = half_open_end
                // So if Hifiasm has x_pos_e as last aligned base, Dinara has qe_ext = x_pos_e + 1
                // Thus: len0 - qe_ext = len0 - (x_pos_e + 1) = (len0 - 1) - x_pos_e = x_right_length ✓
                
                int64_t q_right = (int64_t)len0 - (int64_t)qe_ext;  // Safe signed arithmetic
                int64_t t_right = (int64_t)len1 - (int64_t)te_ext;
                
                if (q_right <= t_right) {
                    qe_ext = (uint32_t)len0;           // qe extends to read end (half-open)
                    te_ext = te_ext + (uint32_t)q_right;  // te extends by same amount
                } else {
                    te_ext = (uint32_t)len1;           // te extends to read end (half-open)
                    qe_ext = qe_ext + (uint32_t)t_right;  // qe extends by same amount
                }
                
                thisAlignmentData.qs = qs_ext;
                thisAlignmentData.qe = qe_ext;
                
                // --- Convert target coordinates to FORWARD STRAND (hifiasm convention) ---
                // When isSameStrand=false (reverse), the ts/te are on the reverse-complement strand.
                // We need to flip them to represent positions on the forward strand.
                if (!candidate.isSameStrand) {
                    // Convert reverse-complement half-open interval [ts_ext, te_ext) to forward [ts, te).
                    const auto p = dinara::rcIntervalToForward(uint32_t(len1), ts_ext, te_ext);
                    thisAlignmentData.ts = p.first;
                    thisAlignmentData.te = p.second;
                } else {
                    thisAlignmentData.ts = ts_ext;
                    thisAlignmentData.te = te_ext;
                }
            }

            // Check for large indels (>= 6 bases)
            thisAlignmentData.hasLargeIndel = projectedAlignment.hasLargeIndel;

            // Cis/Trans Status
            thisAlignmentData.cisTransStatus = CisTransStatus::Unknown;
            
            // Alignment covers an informative het site (counted per read perspective during EC parity)
            thisAlignmentData.informativeHetSiteCount0 = 0;
            thisAlignmentData.informativeHetSiteCount1 = 0;
            thisAlignmentData.informativeHetSiteScore = 0;
            
            // Alignment deletion reasons (none by default)
            thisAlignmentData.deleteReasons0 = AlignmentData::DeleteReasonNone;
            thisAlignmentData.deleteReasons1 = AlignmentData::DeleteReasonNone;

            // --- Populate AlignedEvidenceStore (APES/TASSD) ---
            // Store sparse mismatch/indel evidence (no per-base trace scanning).
            {
                AlignedEvidenceStore& store = data.threadEvidenceStores[threadId];
                thisAlignmentData.info.alignmentId = store.beginAlignment();

                const LongBaseSequenceView tView = reads->getRead(orientedReadIds[1].getReadId());
                const bool tRev = orientedReadIds[1].getStrand();
                DINARA_ASSERT(tView.baseCount <= uint64_t(SnpEvidence::POS_MASK) + 1ULL);
                const uint32_t tRawLen = uint32_t(tView.baseCount);

                static const uint8_t complementBase[4] = {3, 2, 1, 0};

                // Stream 1 (read0/query-view): store read1 base in the oriented frame.
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

                // Stream 0 (read1/target-view): positions are in read1 forward coordinates.
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
                    // Opposite strand: emit in increasing canonical coordinates and
                    // complement read0 base into read1's forward frame.
                    for(auto it = projectedAlignment.sparseMismatches.rbegin();
                        it != projectedAlignment.sparseMismatches.rend(); ++it) {

                        const uint32_t posOriented = it->position1;
                        DINARA_ASSERT(posOriented < tRawLen);
                        const uint32_t pos = (tRawLen - 1U) - posOriented;
                        DINARA_ASSERT(it->base0 < 4);
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

            // Store AlignmentData and the corresponding compressed alignment (same order).
            threadAlignmentData.push_back(thisAlignmentData);

            // Store the alignment in compressed form.
            dinara::compress(alignment, compressedAlignment);
            thisThreadCompressedAlignments.appendVector(
                compressedAlignment.c_str(),
                compressedAlignment.c_str() + compressedAlignment.size()
            );
            
#if 0
            threadAlignmentData.push_back(thisAlignmentData);



            // --- Populate AlignedEvidenceStore (APES/TASSD) ---
            // The AlignedEvidenceStore uses a dual-stream architecture (Target-View vs Query-View)
            // to allow O(1) lookups of evidence during phasing without re-projecting.
            // We populate both streams simultaneously to avoid duplicate computation.

            AlignedEvidenceStore& store = data.threadEvidenceStores[threadId];
            thisAlignmentData.info.alignmentId = store.beginAlignment(); // Returns ID, implicitly syncs with alignmentData index

            const LongBaseSequenceView qView = reads->getRead(orientedReadIds[0].getReadId());
            const LongBaseSequenceView tView = reads->getRead(orientedReadIds[1].getReadId());
            const bool qRev = orientedReadIds[0].getStrand();
            const bool tRev = orientedReadIds[1].getStrand();
            const uint64_t qRawLen = qView.baseCount;
            const uint64_t tRawLen = tView.baseCount;

            // Helper lambdas to access bases in their RAW (Oriented) coordinate system.
            auto getQ = [&](uint64_t p) { return qRev ? qView[qRawLen - 1 - p].complement() : qView[p]; };
            auto getT = [&](uint64_t p) { return tRev ? tView[tRawLen - 1 - p].complement() : tView[p]; };

            // =================================================================================
            // UNIFIED PHASING LOGIC (Forward & Backward Passes)
            // =================================================================================
            // We populate streams based on Strand Orientation to ensure Monotonic Canonical Coordinates.
            // - Forward Reads (!Rev): Populated in Forward Pass (Low->High Oriented).
            // - Reverse Reads (Rev): Populated in Backward Pass (High->Low Oriented).
            // This handles all cases: SameStrand (F/F, R/R) and DiffStrand (F/R, R/F).

            const auto& segments = projectedAlignment.segments;

            // --- PASS 1: FORWARD ITERATION (Handles Forward Reads) ---
            if (!qRev || !tRev) {
                // Initialize pointers to Alignment Start (Oriented)
                uint32_t currentQ = thisAlignmentData.qs;
                uint32_t currentT = thisAlignmentData.ts;
                
                uint32_t lastSnpQ = thisAlignmentData.qs;
                uint32_t lastSnpT = thisAlignmentData.ts;

                if (!segments.empty()) {
                    // Left Tail Handling
                    // Implicit match from start to first segment
                    uint32_t firstQ = segments.front().positionsA[0];
                    if (firstQ > currentQ) {
                        uint32_t gapLen = firstQ - currentQ;
                        uint32_t gapT   = segments.front().positionsA[1] - currentT;

                        if (!qRev) { // Stream 1
                            uint32_t adv = 0;
                            while(adv + SnpEvidence::MAX_DELTA <= gapLen) {
                                store.addSnp1(SnpEvidence::MAX_DELTA, getQ(lastSnpQ + SnpEvidence::MAX_DELTA).value);
                                lastSnpQ += SnpEvidence::MAX_DELTA; adv += SnpEvidence::MAX_DELTA;
                            }
                        }
                        if (!tRev) { // Stream 0
                            uint32_t adv = 0;
                            while(adv + SnpEvidence::MAX_DELTA <= gapT) {
                                store.addSnp0(SnpEvidence::MAX_DELTA, getT(lastSnpT + SnpEvidence::MAX_DELTA).value);
                                lastSnpT += SnpEvidence::MAX_DELTA; adv += SnpEvidence::MAX_DELTA;
                            }
                        }
                    }
                    // Sync up
                    currentQ = segments.front().positionsA[0];
                    currentT = segments.front().positionsA[1];

                    for (const auto& segment : segments) {
                        // 1. Implicit Gap (between previous segment end and current segment start)
                        uint32_t gapQ = segment.positionsA[0] - currentQ;
                        uint32_t gapT = segment.positionsA[1] - currentT;

                        if (!qRev) {
                            uint32_t adv = 0;
                            while(adv + SnpEvidence::MAX_DELTA <= gapQ) {
                                store.addSnp1(SnpEvidence::MAX_DELTA, getQ(lastSnpQ + SnpEvidence::MAX_DELTA).value);
                                lastSnpQ += SnpEvidence::MAX_DELTA; adv += SnpEvidence::MAX_DELTA;
                            }
                        }
                        if (!tRev) {
                            uint32_t adv = 0;
                            while(adv + SnpEvidence::MAX_DELTA <= gapT) {
                                store.addSnp0(SnpEvidence::MAX_DELTA, getT(lastSnpT + SnpEvidence::MAX_DELTA).value);
                                lastSnpT += SnpEvidence::MAX_DELTA; adv += SnpEvidence::MAX_DELTA;
                            }
                        }

                        currentQ = segment.positionsA[0];
                        currentT = segment.positionsA[1];

                        // 2. Steps
                        for(auto it = segment.alignment.begin(); it != segment.alignment.end(); ) {
                            bool advQ = it->first;
                            bool advT = it->second;
                            
                            if (advQ && advT) { // Match/Mismatch
                                if (getQ(currentQ) != getT(currentT)) {
                                    // Mismatch: Store Partner Base
                                    // Use getQ/getT directly (already returns Aligned/Canonical Bases)
                                    if (!qRev) { // Stream 1 expects Target Base
                                        uint32_t d = currentQ - lastSnpQ;
                                        while(d > SnpEvidence::MAX_DELTA) {
                                            store.addSnp1(SnpEvidence::MAX_DELTA, getQ(lastSnpQ + SnpEvidence::MAX_DELTA).value);
                                            lastSnpQ += SnpEvidence::MAX_DELTA; d -= SnpEvidence::MAX_DELTA;
                                        }
                                        store.addSnp1((uint16_t)d, getT(currentT).value);
                                        lastSnpQ = currentQ;
                                    }
                                    if (!tRev) { // Stream 0 expects Query Base
                                        uint32_t d = currentT - lastSnpT;
                                        while(d > SnpEvidence::MAX_DELTA) {
                                            store.addSnp0(SnpEvidence::MAX_DELTA, getT(lastSnpT + SnpEvidence::MAX_DELTA).value);
                                            lastSnpT += SnpEvidence::MAX_DELTA; d -= SnpEvidence::MAX_DELTA;
                                        }
                                        store.addSnp0((uint16_t)d, getQ(currentQ).value);
                                        lastSnpT = currentT;
                                    }
                                }
                                currentQ++; currentT++;
                                ++it;
                            } else {
                                uint32_t len = 0;
                                auto it2 = it;
                                while(it2 != segment.alignment.end() && it2->first == advQ && it2->second == advT) {
                                    len++;
                                    ++it2;
                                }
                                if (!advQ && advT) { // Gap in Q
                                    if (!qRev) store.addIndel1(currentQ, len, 0); // Ins in T relative to Q
                                    if (!tRev) store.addIndel0(currentT, len, 1); // Del in Q relative to T
                                    currentT += len;
                                } else if (advQ && !advT) { // Gap in T
                                    if (!qRev) store.addIndel1(currentQ, len, 1); // Del in T relative to Q
                                    if (!tRev) store.addIndel0(currentT, len, 0); // Ins in Q relative to T
                                    currentQ += len;
                                }
                                it = it2;
                            }
                        }
                        currentQ = segment.positionsB[0];
                        currentT = segment.positionsB[1];
                    }
                    
                    // Right Tail Handling
                    uint32_t endQ = thisAlignmentData.qe;
                    uint32_t endT = thisAlignmentData.te;
                    
                    if (!qRev && endQ > currentQ) {
                        uint32_t gap = endQ - currentQ;
                        uint32_t adv = 0;
                        while(adv + SnpEvidence::MAX_DELTA <= gap) {
                            store.addSnp1(SnpEvidence::MAX_DELTA, getQ(lastSnpQ + SnpEvidence::MAX_DELTA).value);
                            lastSnpQ += SnpEvidence::MAX_DELTA; adv += SnpEvidence::MAX_DELTA;
                        }
                    }
                    if (!tRev && endT > currentT) {
                        uint32_t gap = endT - currentT;
                        uint32_t adv = 0;
                        while(adv + SnpEvidence::MAX_DELTA <= gap) {
                            store.addSnp0(SnpEvidence::MAX_DELTA, getT(lastSnpT + SnpEvidence::MAX_DELTA).value);
                            lastSnpT += SnpEvidence::MAX_DELTA; adv += SnpEvidence::MAX_DELTA;
                        }
                    }
                }
            } // End Pass 1

            // --- PASS 2: BACKWARD ITERATION (Handles Reverse Reads) ---
            if (qRev || tRev) {
                // Initialize pointers to Oriented High
                // Query: qs/qe are Oriented. High = qe.
                // Target: ts/te are Canonical. High = Len - ts.
                uint32_t currentQ = thisAlignmentData.qe;
                uint32_t currentT = tRawLen - thisAlignmentData.ts;

                // For Delta Encoding on Reverse Reads:
                // We iterate Logic High Oriented -> Logic Low Oriented.
                // This corresponds to Canonical Low -> Canonical High.
                // Query: Oriented High qe maps to Canonical Low (Len - qe).
                // Target: Oriented High (Len - ts) maps to Canonical Low ts.
                uint32_t lastSnpQ = qRawLen - thisAlignmentData.qe;
                uint32_t lastSnpT = thisAlignmentData.ts;

                if (!segments.empty()) {
                    // Iterate segments in Reverse
                    for (auto it = segments.rbegin(); it != segments.rend(); ++it) {
                        const auto& segment = *it;

                        // 1. Implicit Gap (High Oriented side)
                        // Gap between current (High) and segment End (High).
                        // In Canonical space, this is the gap BEFORE the segment (Low->High).
                        uint32_t gapQ = currentQ - segment.positionsB[0];
                        uint32_t gapT = currentT - segment.positionsB[1];
                        
                        if (qRev) {
                            uint32_t adv = 0;
                            while(adv + SnpEvidence::MAX_DELTA <= gapQ) {
                                // Dummy at Canonical `lastSnpQ + SnpEvidence::MAX_DELTA`.
                                // To get base, we Map Canonical -> Oriented.
                                // Oriented = qRawLen - 1 - (lastSnpQ + SnpEvidence::MAX_DELTA).
                                store.addSnp1(SnpEvidence::MAX_DELTA, getQ(qRawLen - 1 - (lastSnpQ + SnpEvidence::MAX_DELTA)).value);
                                lastSnpQ += SnpEvidence::MAX_DELTA; adv += SnpEvidence::MAX_DELTA;
                            }
                        }
                        if (tRev) {
                            uint32_t adv = 0;
                            while(adv + SnpEvidence::MAX_DELTA <= gapT) {
                                store.addSnp0(SnpEvidence::MAX_DELTA, getT(tRawLen - 1 - (lastSnpT + SnpEvidence::MAX_DELTA)).value);
                                lastSnpT += SnpEvidence::MAX_DELTA; adv += SnpEvidence::MAX_DELTA;
                            }
                        }

                        currentQ = segment.positionsB[0];
                        currentT = segment.positionsB[1];
                        
                        // 2. Steps (Reverse)
                        for (auto sIt = segment.alignment.rbegin(); sIt != segment.alignment.rend(); ) {
                            bool advQ = sIt->first;
                            bool advT = sIt->second;
                            
                            if (advQ && advT) { // Match/Mismatch
                                currentQ--; currentT--;
                                if (getQ(currentQ) != getT(currentT)) {
                                    // Canonical coords
                                    uint32_t canQ = qRawLen - 1 - currentQ;
                                    uint32_t canT = tRawLen - 1 - currentT;

                                    if (qRev) {
                                        uint32_t d = canQ - lastSnpQ;
                                        while(d > SnpEvidence::MAX_DELTA) {
                                            store.addSnp1(SnpEvidence::MAX_DELTA, getQ(qRawLen - 1 - (lastSnpQ + SnpEvidence::MAX_DELTA)).value);
                                            lastSnpQ += SnpEvidence::MAX_DELTA; d -= SnpEvidence::MAX_DELTA;
                                        }
                                        store.addSnp1((uint16_t)d, getT(currentT).value); // Partner Base
                                        lastSnpQ = canQ;
                                    }
                                    if (tRev) {
                                        uint32_t d = canT - lastSnpT;
                                        while(d > SnpEvidence::MAX_DELTA) {
                                            store.addSnp0(SnpEvidence::MAX_DELTA, getT(tRawLen - 1 - (lastSnpT + SnpEvidence::MAX_DELTA)).value);
                                            lastSnpT += SnpEvidence::MAX_DELTA; d -= SnpEvidence::MAX_DELTA;
                                        }
                                        store.addSnp0((uint16_t)d, getQ(currentQ).value); // Partner Base
                                        lastSnpT = canT;
                                    }
                                }
                                ++sIt;
                            } else {
                                uint32_t len = 0;
                                auto sIt2 = sIt;
                                while(sIt2 != segment.alignment.rend() && sIt2->first == advQ && sIt2->second == advT) {
                                    len++;
                                    ++sIt2;
                                }
                                if (!advQ && advT) { // Gap in Q -> Dec T
                                    // Canonical T: tRawLen - 1 - newest currentT
                                    // newest currentT = currentT - len
                                    // Canonical start pos = tRawLen - 1 - (currentT - 1)
                                    if (tRev) store.addIndel0(tRawLen - 1 - (currentT - 1), len, 1); // Del in Q
                                    if (qRev) store.addIndel1(qRawLen - 1 - currentQ, len, 0); // Ins in T
                                    currentT -= len;
                                } else if (advQ && !advT) { // Gap in T -> Dec Q
                                    if (qRev) store.addIndel1(qRawLen - 1 - (currentQ - 1), len, 1); // Del in T
                                    if (tRev) store.addIndel0(tRawLen - 1 - currentT, len, 0); // Ins in Q
                                    currentQ -= len;
                                }
                                sIt = sIt2;
                            }
                        }
                        currentQ = segment.positionsA[0];
                        currentT = segment.positionsA[1];
                    }

                    // Tail Handling (Low Oriented / High Canonical)
                    // Corresponds to "Right Tail" of Canonical. 
                    // Oriented Low Target: tRawLen - 1 - te.
                    // Oriented Low Query: qRawLen - 1 - qe.
                    
                    uint32_t endQ_Oriented = thisAlignmentData.qs;
                    uint32_t endT_Oriented = tRawLen - thisAlignmentData.te;
                    
                    if (qRev && currentQ > endQ_Oriented) {
                        uint32_t gap = currentQ - endQ_Oriented;
                        uint32_t adv = 0;
                        while(adv + SnpEvidence::MAX_DELTA <= gap) {
                            store.addSnp1(SnpEvidence::MAX_DELTA, getQ(qRawLen - 1 - (lastSnpQ + SnpEvidence::MAX_DELTA)).value);
                            lastSnpQ += SnpEvidence::MAX_DELTA; adv += SnpEvidence::MAX_DELTA;
                        }
                    }
                    if (tRev && currentT > endT_Oriented) {
                        uint32_t gap = currentT - endT_Oriented;
                        uint32_t adv = 0;
                        while(adv + SnpEvidence::MAX_DELTA <= gap) {
                            store.addSnp0(SnpEvidence::MAX_DELTA, getT(tRawLen - 1 - (lastSnpT + SnpEvidence::MAX_DELTA)).value);
                            lastSnpT += SnpEvidence::MAX_DELTA; adv += SnpEvidence::MAX_DELTA;
                        }
                    }
                }
            } // End Pass 2
#endif
        }
    }

    if(alignmentMethod == 4) {

        std::lock_guard<std::mutex> lock(mutex);
        cout << "Thread " << threadId << " byte allocator: " <<
            byteAllocator.getMaxAllocatedByteCount() << "/" <<
            2ULL * 1024 * 1024 * 1024 << endl;
    }

    thisThreadCompressedAlignments.unreserve();

#if 0
    // Write the elapsed time taken by each alignment.
    ofstream csv("AlignTime-Thread" + to_string(threadId) + ".csv");
    for(const auto& p: elapsedTime) {
        const uint64_t i = p.first;
        const OrientedReadPair& candidate = alignmentCandidates.candidates[i];
        DINARA_ASSERT(candidate.readIds[0] < candidate.readIds[1]);
        const OrientedReadId orientedReadId0 = OrientedReadId(candidate.readIds[0], 0);
        const OrientedReadId orientedReadId1 = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);
        const double t = p.second;
        csv << i << ",";
        csv << orientedReadId0 << ",";
        csv << orientedReadId1 << ",";
        csv << t << "\n";
    }
#endif
}



void Assembler::accessCompressedAlignments()
{
    compressedAlignments.accessExistingReadOnly(
        largeDataName("CompressedAlignments"));
}



// Compute alignmentTable from alignmentData.
// This could be made multithreaded if it becomes a bottleneck.
void Assembler::computeAlignmentTable()
{
    // Avoid rebuilding if already computed in this process.
    // Some pipelines call alignment computation routines more than once during experimentation;
    // rebuilding the alignment table is expensive and also fails if the memory-mapped object
    // is still open.
    if(alignmentTable.isOpen()) {
        performanceLog << timestamp << "Alignment table already exists - skipping recomputation." << endl;
        return;
    }

    alignmentTable.createNew(largeDataName("AlignmentTable"), largeDataPageSize);
    alignmentTable.beginPass1(ReadId(2 * reads->readCount()));
    for(const AlignmentData& ad: alignmentData) {
        const auto& readIds = ad.readIds;
        OrientedReadId orientedReadId0(readIds[0], 0);
        OrientedReadId orientedReadId1(readIds[1], ad.isSameStrand ? 0 : 1);
        alignmentTable.incrementCount(orientedReadId0.getValue());
        alignmentTable.incrementCount(orientedReadId1.getValue());
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        alignmentTable.incrementCount(orientedReadId0.getValue());
        alignmentTable.incrementCount(orientedReadId1.getValue());
    }
    alignmentTable.beginPass2();
    for(uint32_t i=0; i<alignmentData.size(); i++) {
        const AlignmentData& ad = alignmentData[i];
        const auto& readIds = ad.readIds;
        OrientedReadId orientedReadId0(readIds[0], 0);
        OrientedReadId orientedReadId1(readIds[1], ad.isSameStrand ? 0 : 1);
        alignmentTable.store(orientedReadId0.getValue(), i);
        alignmentTable.store(orientedReadId1.getValue(), i);
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        alignmentTable.store(orientedReadId0.getValue(), i);
        alignmentTable.store(orientedReadId1.getValue(), i);
    }
    alignmentTable.endPass2();



    // Sort each section of the alignment table by OrientedReadId.
    vector< pair<OrientedReadId, uint32_t> > v;
    for(ReadId readId0=0; readId0<reads->readCount(); readId0++) {
        for(Strand strand0=0; strand0<2; strand0++) {
            const OrientedReadId orientedReadId0(readId0, strand0);

            // Access the section of the alignment table for this oriented read.
            const span<uint32_t> alignmentTableSection =
                alignmentTable[orientedReadId0.getValue()];

            // Store pairs(OrientedReadId, alignmentIndex).
            v.clear();
            for(uint32_t alignmentIndex: alignmentTableSection) {
                const AlignmentData& alignment = alignmentData[alignmentIndex];
                const OrientedReadId orientedReadId1 = alignment.getOther(orientedReadId0);
                v.push_back(make_pair(orientedReadId1, alignmentIndex));
            }

            // Sort.
            sort(v.begin(), v.end());

            // Store the sorted alignmentIndex.
            for(size_t i=0; i<v.size(); i++) {
                alignmentTableSection[i] = v[i].second;
            }
        }
    }

    alignmentTable.unreserve();

}



void Assembler::accessAlignmentData()
{
    alignmentData.accessExistingReadOnly(largeDataName("AlignmentData"));
    alignmentTable.accessExistingReadOnly(largeDataName("AlignmentTable"));
}
void Assembler::accessAlignmentDataReadWrite()
{
    alignmentData.accessExistingReadWrite(largeDataName("AlignmentData"));
    alignmentTable.accessExistingReadWrite(largeDataName("AlignmentTable"));
}



void Assembler::checkAlignmentDataAreOpen() const
{
    if(!alignmentData.isOpen || !alignmentTable.isOpen()) {
        throw runtime_error("Alignment data are not accessible.");
    }
}





// Find in the alignment table the alignments involving
// a given oriented read, and return them with the correct
// orientation (this may involve a swap and/or reverse complement
// of the AlignmentInfo stored in the alignmentTable).
vector< pair<OrientedReadId, dinara::AlignmentInfo> >
    Assembler::findOrientedAlignments(
        OrientedReadId orientedReadId0Argument,
        bool inReadGraphOnly) const
{
    const ReadId readId0 = orientedReadId0Argument.getReadId();
    const ReadId strand0 = orientedReadId0Argument.getStrand();

    vector< pair<OrientedReadId, AlignmentInfo> > result;

    // Loop over alignment involving this read, as stored in the
    // alignment table.
    const auto alignmentTable0 = alignmentTable[orientedReadId0Argument.getValue()];
    for(const auto i: alignmentTable0) {
        const AlignmentData& ad = alignmentData[i];

        // Get the oriented read ids that the AlignmentData refers to.
        OrientedReadId orientedReadId0(ad.readIds[0], 0);
        OrientedReadId orientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
        AlignmentInfo alignmentInfo = ad.info;

        // Skip it if it is not in the read graph and only alignments
        // in the read graph were requested.
        if(inReadGraphOnly and (not alignmentInfo.isInReadGraph)) {
            continue;
        }

        // Swap oriented reads, if necessary.
        if(orientedReadId0.getReadId() != readId0) {
            swap(orientedReadId0, orientedReadId1);
            alignmentInfo.swap();
        }
        DINARA_ASSERT(orientedReadId0.getReadId() == readId0);

        // Reverse complement, if necessary.
        if(orientedReadId0.getStrand() != strand0) {
            orientedReadId0.flipStrand();
            orientedReadId1.flipStrand();
            alignmentInfo.reverseComplement();
        }
        DINARA_ASSERT(orientedReadId0.getStrand() == strand0);
        DINARA_ASSERT(orientedReadId0 == orientedReadId0Argument);

        result.push_back(make_pair(orientedReadId1, alignmentInfo));
    }
    return result;
}



void Assembler::getAlignmentIdsSortedByInformativeSites(
    OrientedReadId orientedReadId0,
    vector<uint32_t>& alignmentIds,
    bool inReadGraphOnly) const
{
    checkAlignmentDataAreOpen();
    DINARA_ASSERT(alignmentTable.isOpen());

    const ReadId readId0 = orientedReadId0.getReadId();
    const span<const uint32_t> table = alignmentTable[orientedReadId0.getValue()];

    alignmentIds.clear();
    alignmentIds.reserve(table.size());
    for(const uint32_t alignmentId : table) {
        const AlignmentData& ad = alignmentData[alignmentId];
        if(inReadGraphOnly && (ad.info.isInReadGraph == 0)) {
            continue;
        }
        alignmentIds.push_back(alignmentId);
    }

    auto overlapLenOnRead = [&](const AlignmentData& ad) -> uint64_t {
        if(ad.readIds[0] == readId0) {
            return uint64_t(ad.qe) - uint64_t(ad.qs);
        } else {
            return uint64_t(ad.te) - uint64_t(ad.ts);
        }
    };

    sort(
        alignmentIds.begin(),
        alignmentIds.end(),
        [&](uint32_t aId, uint32_t bId) {
            const AlignmentData& a = alignmentData[aId];
            const AlignmentData& b = alignmentData[bId];
            const uint32_t aCount = a.getInformativeHetSiteCountFromReadPerspective(readId0);
            const uint32_t bCount = b.getInformativeHetSiteCountFromReadPerspective(readId0);
            if(aCount != bCount) {
                return aCount > bCount;
            }
            const uint64_t aLen = overlapLenOnRead(a);
            const uint64_t bLen = overlapLenOnRead(b);
            if(aLen != bLen) {
                return aLen > bLen;
            }
            return aId < bId;
        });
}



void Assembler::getCisAlignmentIdsSortedByInformativeSites(
    OrientedReadId orientedReadId0,
    vector<uint32_t>& alignmentIds,
    bool keptByBothSidesOnly) const
{
    checkAlignmentDataAreOpen();
    DINARA_ASSERT(alignmentTable.isOpen());

    const ReadId readId0 = orientedReadId0.getReadId();
    const span<const uint32_t> table = alignmentTable[orientedReadId0.getValue()];

    alignmentIds.clear();
    alignmentIds.reserve(table.size());
    for(const uint32_t alignmentId : table) {
        const AlignmentData& ad = alignmentData[alignmentId];
        if(!ad.isCisByBothSides()) {
            continue;
        }
        if(keptByBothSidesOnly && !ad.keptByBothSides()) {
            continue;
        }
        alignmentIds.push_back(alignmentId);
    }

    auto overlapLenOnRead = [&](const AlignmentData& ad) -> uint64_t {
        if(ad.readIds[0] == readId0) {
            return uint64_t(ad.qe) - uint64_t(ad.qs);
        } else {
            return uint64_t(ad.te) - uint64_t(ad.ts);
        }
    };

    sort(
        alignmentIds.begin(),
        alignmentIds.end(),
        [&](uint32_t aId, uint32_t bId) {
            const AlignmentData& a = alignmentData[aId];
            const AlignmentData& b = alignmentData[bId];
            if(a.informativeHetSiteScore != b.informativeHetSiteScore) {
                return a.informativeHetSiteScore > b.informativeHetSiteScore;
            }
            const uint64_t aLen = overlapLenOnRead(a);
            const uint64_t bLen = overlapLenOnRead(b);
            if(aLen != bLen) {
                return aLen > bLen;
            }
            return aId < bId;
        });
}

void Assembler::getAllCisAlignmentIdsSortedByInformativeSites(
    vector<uint32_t>& alignmentIds,
    bool keptByBothSidesOnly) const
{
    checkAlignmentDataAreOpen();

    const uint64_t alignmentCount = alignmentData.size();

    alignmentIds.clear();
    alignmentIds.reserve(alignmentCount);
    for(uint32_t alignmentId = 0; alignmentId < alignmentCount; ++alignmentId) {
        const AlignmentData& ad = alignmentData[alignmentId];
        if(!ad.isCisByBothSides()) {
            continue;
        }
        if(keptByBothSidesOnly && !ad.keptByBothSides()) {
            continue;
        }
        alignmentIds.push_back(alignmentId);
    }

    sort(
        alignmentIds.begin(),
        alignmentIds.end(),
        [&](uint32_t aId, uint32_t bId) {
            const AlignmentData& a = alignmentData[aId];
            const AlignmentData& b = alignmentData[bId];
            if(a.informativeHetSiteScore != b.informativeHetSiteScore) {
                return a.informativeHetSiteScore > b.informativeHetSiteScore;
            }
            return aId < bId;
        });
}



// Flag palindromic reads.
void Assembler::flagPalindromicReads(
    uint32_t maxSkip,
    uint32_t maxDrift,
    uint32_t maxMarkerFrequency,
    double alignedFractionThreshold,
    double nearDiagonalFractionThreshold,
    uint32_t deltaThreshold,
    uint64_t threadCount)
{
    performanceLog << timestamp << "Finding palindromic reads." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store the parameters so all threads can see them.
    flagPalindromicReadsData.maxSkip = maxSkip;
    flagPalindromicReadsData.maxDrift = maxDrift;
    flagPalindromicReadsData.maxMarkerFrequency = maxMarkerFrequency;
    flagPalindromicReadsData.alignedFractionThreshold = alignedFractionThreshold;
    flagPalindromicReadsData.nearDiagonalFractionThreshold = nearDiagonalFractionThreshold;
    flagPalindromicReadsData.deltaThreshold = deltaThreshold;

    // Reset all palindromic flags.
    reads->assertReadsAndFlagsOfSameSize();
    const ReadId readCount = reads->readCount();
    for(ReadId readId=0; readId<readCount; readId++) {
        reads->setPalindromicFlag(readId, false);
    }

    // Do it in parallel.
    setupLoadBalancing(readCount, 1000);
    runThreads(&Assembler::flagPalindromicReadsThreadFunction, threadCount);

    // Count the reads flagged as palindromic.
    size_t palindromicReadCount = 0;
    for(ReadId readId=0; readId<readCount; readId++) {
        if(reads->getFlags(readId).isPalindromic) {
            ++palindromicReadCount;
        }
    }
    assemblerInfo->palindromicReadCount = palindromicReadCount;
    cout << "Flagged " << palindromicReadCount <<
        " reads as palindromic out of " << readCount << " total." << endl;
    cout << "Palindromic fraction is " <<
        double(palindromicReadCount)/double(readCount) << endl;

}



void Assembler::flagPalindromicReadsThreadFunction(uint64_t)
{

    // Work areas used inside the loop and defined here
    // to reduce memory allocation activity.
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;

    // Make local copies of the parameters.
    // const uint32_t maxSkip = flagPalindromicReadsData.maxSkip;
    // const uint32_t maxDrift = flagPalindromicReadsData.maxDrift;
    // const uint32_t maxMarkerFrequency = flagPalindromicReadsData.maxMarkerFrequency;
    const double alignedFractionThreshold = flagPalindromicReadsData.alignedFractionThreshold;
    const double nearDiagonalFractionThreshold = flagPalindromicReadsData.nearDiagonalFractionThreshold;
    const uint32_t deltaThreshold = flagPalindromicReadsData.deltaThreshold;
    const int maxUncoveredBases = flagPalindromicReadsData.maxUncoveredBases;

    // Default values for alignOrientedReads5 (LowHash), matching AssemblerOptions defaults
    const int matchScore = 6;
    const int mismatchScore = -1;
    const int gapScore = -1;
    const double align5DriftRateTolerance = 0.05;
    const uint64_t align5MinBandExtend = 10;
    ofstream nullStream;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    reads->assertReadsAndFlagsOfSameSize();

    while(getNextBatch(begin, end)) {

        // Loop over all reads in this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            alignOrientedReads5(OrientedReadId(readId, 0), OrientedReadId(readId, 1),
                        matchScore, mismatchScore, gapScore,
                        align5DriftRateTolerance, align5MinBandExtend,
                        alignment, alignmentInfo,
                        nullStream);

            // Calculate Metrics
            const uint32_t alignedMarkerCount = alignmentInfo.markerCount;
            const uint32_t totalMarkerCount = uint32_t((*markers)[OrientedReadId(readId, 0).getValue()].size());
            const double alignedFraction = alignmentInfo.alignedFraction(0);
            const uint32_t alignmentRange = alignmentInfo.range(0);

            // Metric 1: Extended Span (for full-length palindromes with gaps/adapters)
            bool isSpanPass = false;
            uint64_t unaligned = 0; // Declare outside to be visible for debug

            // We require *some* markers to calculate span, but low count is acceptable for 'gappy' alignments.
            if (alignedMarkerCount > 30) { 

                uint64_t readLength = reads->getReadRawSequenceLength(readId);
                
                uint64_t unaligned = (readLength > alignmentRange) ? (readLength - alignmentRange) : 0;
                 

                if (int(unaligned) <= maxUncoveredBases) {
                    isSpanPass = true;
                }
            }

            // Metric 2: Aligned Fraction (Legacy/Density Check)
            // If the span is good, we don't care about the fraction (handles large gaps).
            bool isEnoughMarkers = isSpanPass || (alignedFraction >= alignedFractionThreshold);

            if(!isEnoughMarkers) {
                // Debug Logging for targeted read if it FAILs
                if (readId == 1180) {
                    {
                        cout << "DEBUG Palindrome 1180 FAILED: alignedMarkerCount=" << alignedMarkerCount 
                            << " alignedFraction=" << alignedFraction
                            << " unalignedBases=" << unaligned << endl;
                    }
                }
                continue;
            }

            // If the alignment has too few markers near the diagonal, skip it.
            size_t nearDiagonalMarkerCount = 0;
            for(size_t i=0; i<alignment.ordinals.size(); i++) {
                const array<uint32_t, 2>& ordinals = alignment.ordinals[i];
                const int32_t ordinal0 = int32_t(ordinals[0]);
                const int32_t ordinal1 = int32_t(ordinals[1]);
                const uint32_t delta = abs(ordinal0 - ordinal1);
                if(delta < deltaThreshold) {
                    nearDiagonalMarkerCount++;
                }
            }
            const double nearDiagonalFraction = double(nearDiagonalMarkerCount)/double(totalMarkerCount);
            if(nearDiagonalFraction < nearDiagonalFractionThreshold) {
                if (readId == 1180) {
                    {
                        cout << "DEBUG Palindrome 1180 FAILED: nearDiagonalFraction=" << nearDiagonalFraction << endl;
                    }
                }
                continue;
            }

            // If we got here, mark the read as palindromic.
            reads->setPalindromicFlag(readId, true);

        }
    }
}



void Assembler::analyzeAlignmentMatrix(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1)
{
    // Get the oriented reads.
    const OrientedReadId orientedReadId0(readId0, strand0);
    const OrientedReadId orientedReadId1(readId1, strand1);

    // Get the markers sorted by kmerId.
    vector<MarkerWithOrdinal> markers0;
    vector<MarkerWithOrdinal> markers1;
    getMarkersSortedByKmerId(orientedReadId0, markers0);
    getMarkersSortedByKmerId(orientedReadId1, markers1);

    // Some iterators we will need.
    using MarkerIterator = vector<MarkerWithOrdinal>::const_iterator;
    const MarkerIterator begin0 = markers0.begin();
    const MarkerIterator end0   = markers0.end();
    const MarkerIterator begin1 = markers1.begin();
    const MarkerIterator end1   = markers1.end();

    // The number of markers in each oriented read.
    const int64_t n0 = end0 - begin0;
    const int64_t n1 = end1 - begin1;

    // We will use coordinates
    // x = ordinal0 + ordinal1
    // y = ordinal0 - ordinal1 (offset)
    // In these coordinates, diagonals in the alignment matrix
    // are lines of constant y and so they become horizontal.
    const int64_t xMin = 0;
    const int64_t xMax = n0 + n1 - 2;
    const int64_t yMin = -n1;
    const int64_t yMax = n0 - 1;
    const int64_t nx= xMax -xMin + 1;
    const int64_t ny= yMax -yMin + 1;

    // Create a histogram in cells of size (dx, dy).
    const int64_t dx = 100;
    const int64_t dy = 20;
    const int64_t nxCells = (nx-1)/dx + 1;
    const int64_t nyCells = (ny-1)/dy + 1;
    vector< vector<uint64_t> > histogram(nxCells, vector<uint64_t>(nyCells, 0));
    cout << "nxCells " << nxCells << endl;
    cout << "nyCells " << nyCells << endl;



    // Joint loop over the markers, looking for common k-mer ids.
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->kmerId < it1->kmerId) {
            ++it0;
        } else if(it1->kmerId < it0->kmerId) {
            ++it1;
        } else {

            // We found a common k-mer id.
            const KmerId kmerId = it0->kmerId;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in kmers0 and kmers1.
            MarkerIterator it0Begin = it0;
            MarkerIterator it1Begin = it1;
            MarkerIterator it0End = it0Begin;
            MarkerIterator it1End = it1Begin;
            while(it0End!=end0 && it0End->kmerId==kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->kmerId==kmerId) {
                ++it1End;
            }


            // Loop over pairs in the streaks.
            for(MarkerIterator jt0=it0Begin; jt0!=it0End; ++jt0) {
                const int64_t ordinal0 = jt0->ordinal;
                for(MarkerIterator jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const int64_t ordinal1 = int64_t(jt1->ordinal);

                    const int64_t x = ordinal0 + ordinal1;
                    const int64_t y = ordinal0 - ordinal1;
                    const int64_t ix = (x-xMin) / dx;
                    const int64_t iy = (y-yMin) / dy;
                    DINARA_ASSERT(ix >= 0);
                    DINARA_ASSERT(iy >= 0);
                    DINARA_ASSERT(ix < nxCells);
                    DINARA_ASSERT(iy < nyCells);

                    ++histogram[ix][iy];
                }
            }


            // Continue joint loop over k-mers.
            it0 = it0End;
            it1 = it1End;
        }
    }

    ofstream csv("Histogram.csv");
    PngImage image = PngImage(int(nxCells), int(nyCells));
    uint64_t activeCellCount = 0;
    for(int64_t iy=0; iy<nyCells; iy++) {
        for(int64_t ix=0; ix<nxCells; ix++) {
            const int64_t frequency = histogram[ix][iy];
            if(frequency > 0) {
                csv << frequency;
            }
            if(frequency >= 10) {
                ++activeCellCount;
            }
            // const int r = (frequency <= 10) ? 0 : min(int(255), int(10*frequency));
            const int r = (frequency>=10) ? 255 : 0;
            // const int r = min(int(255), int(10*frequency));
            const int g = r;
            const int b = r;
            image.setPixel(int(ix), int(iy), r, g, b);
            csv << ",";
        }
        csv << "\n";
    }
    image.write("Histogram.png");
    cout << activeCellCount << " active cells out of " << nxCells*nyCells << endl;


}


// Count the common marker near a given ordinal offset for
// two oriented reads. This can be used to check
// whether an alignmnent near the specified ordinal offset exists.
uint32_t Assembler::countCommonMarkersNearOffset(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int32_t offset,
    int32_t offsetTolerance
)
{
    const int32_t minOffset = offset - offsetTolerance;
    const int32_t maxOffset = offset + offsetTolerance;
    return countCommonMarkersWithOffsetIn(
        orientedReadId0, orientedReadId1,
        minOffset, maxOffset);
}
uint32_t Assembler::countCommonMarkersWithOffsetIn(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int32_t minOffset,
    int32_t maxOffset
)
{
#if 1
    const bool debug = false;
#else
    const bool debug =
        orientedReadId0.getReadId() == 5 and
        orientedReadId0.getStrand() == 0 and
        orientedReadId1.getReadId() == 32 and
        orientedReadId1.getStrand() == 1;
#endif
    if(debug) {
        cout << "countCommonMarkersWithOffsetIn" << endl;
    }

    // Get the markers sorted by kmerId.
    checkMarkersAreOpen();
    vector<MarkerWithOrdinal> markers0;
    vector<MarkerWithOrdinal> markers1;
    getMarkersSortedByKmerId(orientedReadId0, markers0);
    getMarkersSortedByKmerId(orientedReadId1, markers1);

    // Some iterators we will need.
    using MarkerIterator = vector<MarkerWithOrdinal>::const_iterator;
    const MarkerIterator begin0 = markers0.begin();
    const MarkerIterator end0   = markers0.end();
    const MarkerIterator begin1 = markers1.begin();
    const MarkerIterator end1   = markers1.end();



    // Main loop, looking for common markers.
    uint32_t count = 0;
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->kmerId < it1->kmerId) {
            ++it0;
        } else if(it1->kmerId < it0->kmerId) {
            ++it1;
        } else {

            // We found a common marker.
            const KmerId kmerId = it0->kmerId;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in markers0 and markers1.
            MarkerIterator it0Begin = it0;
            MarkerIterator it1Begin = it1;
            MarkerIterator it0End = it0Begin;
            MarkerIterator it1End = it1Begin;
            while(it0End!=end0 && it0End->kmerId==kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->kmerId==kmerId) {
                ++it1End;
            }

            // Loop over pairs of markers in the streaks.
            for(MarkerIterator jt0=it0Begin; jt0!=it0End; ++jt0) {
                for(MarkerIterator jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const int32_t ordinalOffset = int32_t(jt0->ordinal) - int32_t(jt1->ordinal);
                    if(debug) {
                        cout << jt0->ordinal << " " <<
                            jt1->ordinal << " " <<
                            ordinalOffset<< endl;
                    }
                    if(ordinalOffset >= minOffset and ordinalOffset <= maxOffset) {
                        ++count;
                    }

                }
            }
            // Continue main loop.
            it0 = it0End;
            it1 = it1End;
        }

    }

    return count;
}



// Check if an alignment between two reads should be suppressed,
// bases on the setting of command line option
// --Align.sameChannelReadAlignment.suppressDeltaThreshold.
bool Assembler::suppressAlignment(
    ReadId readId0,
    ReadId readId1,
    uint64_t delta)
{

    // If the ch meta data fields of the two reads are missing or different,
    // don't suppress the alignment.
    // Check the channel first for efficiency,
    // so we can return faster in most cases.
    const auto ch0 = reads->getMetaData(readId0, "ch");
    if(ch0.empty()) {
        return false;
    }
    const auto ch1 = reads->getMetaData(readId1, "ch");
    if(ch1.empty()) {
        return false;
    }
    if(ch0 not_eq ch1) {
        return false;
    }



    // If the sampleid meta data fields of the two reads are missing or different,
    // don't suppress the alignment.
    const auto sampleid0 = reads->getMetaData(readId0, "sampleid");
    if(sampleid0.empty()) {
        return false;
    }
    const auto sampleid1 = reads->getMetaData(readId1, "sampleid");
    if(sampleid1.empty()) {
        return false;
    }
    if(sampleid0 not_eq sampleid1) {
        return false;
    }



    // If the runid meta data fields of the two reads are missing or different,
    // don't suppress the alignment.
    const auto runid0 = reads->getMetaData(readId0, "runid");
    if(runid0.empty()) {
        return false;
    }
    const auto runid1 = reads->getMetaData(readId1, "runid");
    if(runid1.empty()) {
        return false;
    }
    if(runid0 not_eq runid1) {
        return false;
    }

    // cout << "Checking " << readId0 << " " << readId1 << endl;


    // If the read meta data fields of the two reads are missing,
    // don't suppress the alignment.
    const auto read0 = reads->getMetaData(readId0, "read");
    if(read0.empty()) {
        return false;
    }
    const auto read1 = reads->getMetaData(readId1, "read");
    if(read1.empty()) {
        return false;
    }
    // cout << read0 << " " << read1 << endl;




    // Convert the read meta data fields to integers.
    // Keep in mind the span<char> is not null-terminated.
    const int64_t r0 = int64_t(atoul(read0));
    const int64_t r1 = int64_t(atoul(read1));
    // cout << r0 << " " << r1 << endl;



    // Suppress the alignment if the absolute difference of the
    // read meta data fields is less than delta.
    return std::abs(r0 - r1) < int64_t(delta);

}



// Remove all alignment candidates for which suppressAlignment
// returns false.
void Assembler::suppressAlignmentCandidates(
    uint64_t delta,
    uint64_t threadCount)
{
    performanceLog << timestamp << "Suppressing alignment candidates." << endl;

    // Allocate memory for flags to keep track of which alignments
    // should be suppressed.
    suppressAlignmentCandidatesData.suppress.createNew(
        largeDataName("tmp-suppressAlignmentCandidates"), largeDataPageSize);
    const uint64_t candidateCount = alignmentCandidates.candidates.size();
    suppressAlignmentCandidatesData.suppress.resize(candidateCount);

    // Figure out which candidates should be suppressed.
    suppressAlignmentCandidatesData.delta = delta;
    const uint64_t batchSize = 10000;
    setupLoadBalancing(alignmentCandidates.candidates.size(), batchSize);
    runThreads(&Assembler::suppressAlignmentCandidatesThreadFunction, threadCount);

    ofstream csv("SuppressedAlignmentCandidates.csv");
    csv << "ReadId0,ReadId1,SameStrand,Name0,Name1,MetaData0,MetaData1" << endl;

    // Suppress the alignment candidates we flagged.
    cout << "Number of alignment candidates before suppression is " << candidateCount << endl;
    uint64_t j = 0;
    uint64_t suppressCount = 0;
    for(uint64_t i=0; i<candidateCount; i++) {
        if(suppressAlignmentCandidatesData.suppress[i]) {
            ++suppressCount;
            const ReadId readId0 = alignmentCandidates.candidates[i].readIds[0];
            const ReadId readId1 = alignmentCandidates.candidates[i].readIds[1];
            csv << readId0 << "," << readId1 << ","
                << (alignmentCandidates.candidates[i].isSameStrand ? "Yes" : "No") << ","
                << reads->getReadName(readId0) << "," << reads->getReadName(readId1) << ","
                << reads->getReadMetaData(readId0) << "," << reads->getReadMetaData(readId1) << endl;
        } else {
            alignmentCandidates.candidates[j++] =
                alignmentCandidates.candidates[i];
        }
    }
    DINARA_ASSERT(j + suppressCount == candidateCount);
    alignmentCandidates.candidates.resize(j);
    cout << "Suppressed " << suppressCount << " alignment candidates." << endl;
    cout << "Number of alignment candidates after suppression is " << j << endl;


    // Clean up.
    suppressAlignmentCandidatesData.suppress.remove();

    performanceLog << timestamp << "Done suppressing alignment candidates." << endl;
}



void Assembler::suppressAlignmentCandidatesThreadFunction(uint64_t)
{
    const uint64_t delta = suppressAlignmentCandidatesData.delta;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over candidate alignments in this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const OrientedReadPair& p = alignmentCandidates.candidates[i];
            // cout << "Checking " << p.readIds[0] << " " << p.readIds[1] <<  " " << int(p.isSameStrand) << endl;
            suppressAlignmentCandidatesData.suppress[i] =
                suppressAlignment(p.readIds[0], p.readIds[1], delta);
        }

    }
}
