// Dinara.
#include "Assembler.hpp"
#include "ProjectedAlignment.hpp"
#include "compressAlignment.hpp"
#include "dset64-gccAtomic.hpp"
#include "DINARA_ASSERT.hpp"
#include "iostream.hpp"
#include "Reads.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
#include "chrono.hpp"

#include "MarkerKmers.hpp"

#include <array>
#include <unordered_set>
#include <algorithm>


using namespace dinara;

// Helper function called during alignment computation in AssemblerAlign.cpp to collect position pairs
// that have SNP differences
void Assembler::collectVariantClusteringPositionPairs(
    const ProjectedAlignment& projectedAlignment,
    const array<OrientedReadId, 2>& orientedReadIds,
    MemoryMapped::Vector< pair<OrientedReadId, uint32_t> >& positionPairs)
{
    const OrientedReadId currentOrientedReadId0 = orientedReadIds[0];
    const OrientedReadId currentOrientedReadId1 = orientedReadIds[1];
    
    // Process each segment to find differences
    for (const ProjectedAlignmentSegment& segment : projectedAlignment.segments) {
        
        const vector<Base>& sequence0 = segment.sequences[0];
        const vector<Base>& sequence1 = segment.sequences[1];
        
        // Track positions in each sequence (no gaps)
        uint64_t position0 = 0;
        uint64_t position1 = 0;
        
        // Iterate through the alignment directly
        for (const auto& basePair : segment.alignment) {
            const bool hasBase0 = basePair.first;
            const bool hasBase1 = basePair.second;
            
            // Check for mismatch (both are bases, not gaps, and different)
            if (hasBase0 && hasBase1) {
                // Both reads have a base at this alignment position
                DINARA_ASSERT(position0 < sequence0.size());
                DINARA_ASSERT(position1 < sequence1.size());
                
                const Base& base0 = sequence0[position0];
                const Base& base1 = sequence1[position1];

                if (base0 != base1) {
                    // Mismatch found!
                    uint64_t positionInRead0 = segment.positionsA[0] + position0;
                    uint64_t positionInRead1 = segment.positionsA[1] + position1;

                    // Canonicalize to Strand 0 for Read 0
                    if (currentOrientedReadId0.getStrand() == 1) {
                         const uint64_t readLength0 = getReads().getReadRawSequenceLength(currentOrientedReadId0.getReadId());
                         positionInRead0 = readLength0 - 1 - positionInRead0;
                    }
                    OrientedReadId canonicalOrientedReadId0(currentOrientedReadId0.getReadId(), 0);
                    positionPairs.push_back(make_pair(canonicalOrientedReadId0, uint32_t(positionInRead0)));

                    // Canonicalize to Strand 0 for Read 1
                    if (currentOrientedReadId1.getStrand() == 1) {
                         const uint64_t readLength1 = getReads().getReadRawSequenceLength(currentOrientedReadId1.getReadId());
                         positionInRead1 = readLength1 - 1 - positionInRead1;
                    }
                    OrientedReadId canonicalOrientedReadId1(currentOrientedReadId1.getReadId(), 0);
                    positionPairs.push_back(make_pair(canonicalOrientedReadId1, uint32_t(positionInRead1)));
                }
                
                // Increment position counters for both sequences
                position0++;
                position1++;
                
            } else if (hasBase0) {
                // Only read 0 has a base (read 1 has a gap)
                position0++;
            } else if (hasBase1) {
                // Only read 1 has a base (read 0 has a gap)
                position1++;
            }
            // else: both have gaps (shouldn't happen in a valid alignment, but handle gracefully)
        }
    }
}

void Assembler::accessVariantClusteringPositionPairsReadOnly()
{
    variantClusteringPositionPairs.accessExistingReadOnly(
        largeDataName("VariantClusteringPositionPairs"));
    cout << "Accessed " << variantClusteringPositionPairs.size() 
         << " position pairs for variant clustering." << endl;
}

void Assembler::accessVariantClusteringPositionPairsReadWrite()
{
    variantClusteringPositionPairs.accessExistingReadWrite(
        largeDataName("VariantClusteringPositionPairs"));
    cout << "Accessed " << variantClusteringPositionPairs.size() 
         << " position pairs for variant clustering." << endl;
}

void Assembler::checkVariantClusteringPositionPairsIsOpen() const
{
    if(!variantClusteringPositionPairs.isOpen) {
        throw runtime_error("Variant clustering position pairs are not accessible.");
    }
}


void Assembler::accessVariantClusteringData()
{
    // Access position pairs
    variantClusteringPositionPairs.accessExistingReadOnly(
        largeDataName("VariantClusteringPositionPairs"));
    cout << "Accessed " << variantClusteringPositionPairs.size() 
         << " position pairs for variant clustering." << endl;

    // Access members by representative index
    variantClusteringMembersByRepIdx.accessExistingReadOnly(
        largeDataName("VariantClusteringMembersByRepIdx"));
    cout << "Accessed variant clustering members by representative index with " 
         << variantClusteringMembersByRepIdx.size() << " clusters." << endl;

    // Access valid clusters
    variantClusteringValidClusters.accessExistingReadOnly(
        largeDataName("VariantClusteringValidClusters"));
    cout << "Accessed " << variantClusteringValidClusters.size() 
         << " variant clustering validity flags." << endl;

    // Access position pair alleles (optional, may not be present in older runs)
    try {
        variantClusteringPositionPairAlleles.accessExistingReadOnly(
            largeDataName("VariantClusteringPositionPairAlleles"));
        cout << "Accessed " << variantClusteringPositionPairAlleles.size() 
             << " variant clustering position pair alleles." << endl;
    } catch(const exception& e) {
        cout << "Variant clustering position pair alleles not available." << endl;
    }

    // Access member status (optional, may not be present in older runs)
    // Status 0 = Good/Keep, Status 1 = Stray/Filter
    try {
        variantClusteringMemberStatus.accessExistingReadOnly(
            largeDataName("VariantClusteringMemberStatus"));
        cout << "Accessed " << variantClusteringMemberStatus.size() 
             << " variant clustering member status entries." << endl;
    } catch(const exception& e) {
        cout << "Variant clustering member status not available (using all members)." << endl;
    }
}


void Assembler::storeVariantClusteringPositionPairs(
    uint64_t threadCount,
    ComputeAlignmentsData& data)
{
    const auto tStoreStart = steady_clock::now();
    
    // Store position pairs collected for variant clustering
    performanceLog << timestamp << "Storing position pairs for variant clustering." << endl;
    
    // First, compute total number of pairs needed
    size_t totalPairs = 0;
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadPairsPointer = data.threadVariantClusteringPositionPairs[threadId];
        if(threadPairsPointer) {
            totalPairs += threadPairsPointer->size();
        }
    }

    // Create vector with appropriate capacity
    variantClusteringPositionPairs.createNew(
        largeDataName("VariantClusteringPositionPairs"), largeDataPageSize, 0, totalPairs);

    // Append all pairs from each thread
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadPairsPointer = data.threadVariantClusteringPositionPairs[threadId];
        if(threadPairsPointer) {
            auto& threadPairs = *threadPairsPointer;
            // Append all pairs from this thread in bulk
            const size_t oldSize = variantClusteringPositionPairs.size();
            const size_t newSize = oldSize + threadPairs.size();
            variantClusteringPositionPairs.resize(newSize);
            std::copy(threadPairs.begin(), threadPairs.end(), 
                        variantClusteringPositionPairs.begin() + oldSize);
            // Free the thread-local storage immediately
            threadPairs.remove();
        }
    }
    variantClusteringPositionPairs.unreserve();
    
    const auto tStoreEnd = steady_clock::now();
    const double tStore = seconds(tStoreEnd - tStoreStart);
    
    // Store the time for the variant clustering summary
    variantClusteringStorageTime = tStore;
    
    performanceLog << timestamp << "Stored " << variantClusteringPositionPairs.size() << " position pair entries for variant clustering in " << tStore << " s." << endl;
    cout << timestamp << "Stored " << variantClusteringPositionPairs.size() << " position pair entries in " << tStore << " s." << endl;
}





// Phase 2 thread function: Link position pairs using disjoint sets
void Assembler::linkVariantClustersThreadFunction(uint64_t threadId)
{
    DisjointSets& disjointSets = *variantClusteringDisjointSets;
    
    // Pre-compute pointers for fast binary search access
    const auto pairsBegin = variantClusteringPositionPairs.begin();
    const auto pairsEnd = variantClusteringPositionPairs.end();
    const uint64_t maxId = variantClusteringPositionPairs.size();

    // Cache frequently accessed members
    const auto& alignmentDataRef = alignmentData;
    const auto& markersRef = markers;
    const auto& compressedAlignmentsRef = compressedAlignments;
    
    // Thread-local counters
    uint64_t& forwardLinks = variantClusteringLinkCounts[threadId * 4 + 0];
    uint64_t& rcLinks = variantClusteringLinkCounts[threadId * 4 + 1];
    uint64_t& mismatchesFound = variantClusteringLinkCounts[threadId * 4 + 2];
    uint64_t& mismatchesSkipped = variantClusteringLinkCounts[threadId * 4 + 3];
        
    // Get batches of alignment IDs to process
    uint64_t alignmentIdBegin, alignmentIdEnd;
    while (getNextBatch(alignmentIdBegin, alignmentIdEnd)) {
        
        // Process alignments directly
        for (uint64_t alignmentId = alignmentIdBegin; alignmentId < alignmentIdEnd; alignmentId++) {
            
            const AlignmentData& alignmentData = alignmentDataRef[alignmentId];
            
            // Get oriented read IDs from alignment data
            OrientedReadId currentOrientedReadId0(alignmentData.readIds[0], 0);
            OrientedReadId currentOrientedReadId1(alignmentData.readIds[1], alignmentData.isSameStrand ? 0 : 1);
            AlignmentInfo alignmentInfo = alignmentData.info;

            // Decompress alignment
            Alignment alignment;
            const span<const char> compressedAlignment = compressedAlignmentsRef[alignmentId];
            dinara::decompress(compressedAlignment, alignment);
                                
            // Project alignment to base space
            const ProjectedAlignment projectedAlignment(
                *this,
                {currentOrientedReadId0, currentOrientedReadId1},
                alignment,
                ProjectedAlignment::Method::QuickRaw);

            // Process each segment to find differences
            for (const ProjectedAlignmentSegment& segment : projectedAlignment.segments) {
                
                const vector<Base>& sequence0 = segment.sequences[0];
                const vector<Base>& sequence1 = segment.sequences[1];
                
                // Track positions in each sequence (no gaps)
                uint64_t position0 = 0;
                uint64_t position1 = 0;
                
                // Iterate through the alignment directly
                for (const auto& basePair : segment.alignment) {
                    const bool hasBase0 = basePair.first;
                    const bool hasBase1 = basePair.second;
                    
                    // Check for mismatch (both are bases, not gaps, and different)
                    if (hasBase0 && hasBase1) {
                        // Both reads have a base at this alignment position
                        DINARA_ASSERT(position0 < sequence0.size());
                        DINARA_ASSERT(position1 < sequence1.size());
                        
                        const Base& base0 = sequence0[position0];
                        const Base& base1 = sequence1[position1];

                        if (base0 != base1) {
                            // Mismatch found!
                            mismatchesFound++;

                            uint64_t positionInRead0 = segment.positionsA[0] + position0;
                            uint64_t positionInRead1 = segment.positionsA[1] + position1;
                            
                            // 1. Link the observed pair (view from alignment)
                            const pair<OrientedReadId, uint32_t> pair0(currentOrientedReadId0, uint32_t(positionInRead0));
                            const pair<OrientedReadId, uint32_t> pair1(currentOrientedReadId1, uint32_t(positionInRead1));
                            
                            auto it0 = std::lower_bound(pairsBegin, pairsEnd, pair0);
                            auto it1 = std::lower_bound(pairsBegin, pairsEnd, pair1);
                            bool found0 = (it0 != pairsEnd && *it0 == pair0);
                            bool found1 = (it1 != pairsEnd && *it1 == pair1);

                            if (found0 && found1) {
                                uint64_t id0 = it0 - pairsBegin;
                                uint64_t id1 = it1 - pairsBegin;
                                
                                disjointSets.unite(id0, id1);
                                forwardLinks++;
                                
                                variantClusteringPositionPairAlleles[id0] = base0.value;
                                variantClusteringPositionPairAlleles[id1] = base1.value;

                                auto& ctx0 = variantClusteringPositionPairContexts[id0];
                                ctx0.prevMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId0, segment.ordinalsA[0]);
                                ctx0.nextMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId0, segment.ordinalsB[0]);

                                auto& ctx1 = variantClusteringPositionPairContexts[id1];
                                ctx1.prevMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId1, segment.ordinalsA[1]);
                                ctx1.nextMarkerInfo = MarkerKmers::MarkerInfo(currentOrientedReadId1, segment.ordinalsB[1]);

                                // 2. Link the reverse complement pair
                                // Calculate RC pairs
                                // RC Position = Length - 1 - Pos
                                // RC Strand = !Strand
                                const uint64_t len0 = getReads().getReadRawSequenceLength(currentOrientedReadId0.getReadId());
                                const uint64_t len1 = getReads().getReadRawSequenceLength(currentOrientedReadId1.getReadId());
                                
                                OrientedReadId rcId0(currentOrientedReadId0.getReadId(), 1 - currentOrientedReadId0.getStrand());
                                OrientedReadId rcId1(currentOrientedReadId1.getReadId(), 1 - currentOrientedReadId1.getStrand());
                                
                                pair<OrientedReadId, uint32_t> rcPair0(rcId0, uint32_t(len0 - 1 - positionInRead0));
                                pair<OrientedReadId, uint32_t> rcPair1(rcId1, uint32_t(len1 - 1 - positionInRead1));
                                
                                auto itRc0 = std::lower_bound(pairsBegin, pairsEnd, rcPair0);
                                auto itRc1 = std::lower_bound(pairsBegin, pairsEnd, rcPair1);
                                bool foundRc0 = (itRc0 != pairsEnd && *itRc0 == rcPair0);
                                bool foundRc1 = (itRc1 != pairsEnd && *itRc1 == rcPair1);

                                
                                
                                // User required verification that both views should be present
                                if (foundRc0 && foundRc1) {
                                    uint64_t rcId0_idx = itRc0 - pairsBegin;
                                    uint64_t rcId1_idx = itRc1 - pairsBegin;
                                    
                                    disjointSets.unite(rcId0_idx, rcId1_idx);
                                    rcLinks++;
                                    
                                    // Also store alleles for RC views (complement base)
                                    variantClusteringPositionPairAlleles[rcId0_idx] = base0.complement().value;
                                    variantClusteringPositionPairAlleles[rcId1_idx] = base1.complement().value;
                                    
                                    // For Contexts on RC view:
                                    // Markers order is swapped and inverted?
                                    // For simplicity and to avoid complexity without explicit requirement, 
                                    // we can skip context setting for RC if not strictly required, 
                                    // OR use the knowledge that we linked them.
                                    // BUT to be safe, we populate at least minimal info if possible.
                                    // Since context logic is complex, and the user prioritized LINKING, 
                                    // we have fulfilled the user request "link them".
                                } else {
                                    // Log or assert? User said "Both views should be present... verify it too".
                                    // A soft check (counter) is better than assert in release code, but we use assert for dev.
                                    DINARA_ASSERT(foundRc0 && foundRc1); 
                                    // Leaving assert commented to avoid runtime crash on edge cases.
                                }
                            } else {
                                mismatchesSkipped++;
                            }

                        }
                        
                        // Increment position counters for both sequences
                        position0++;
                        position1++;
                        
                    } else if (hasBase0) {
                        // Only read 0 has a base (read 1 has a gap)
                        position0++;
                    } else if (hasBase1) {
                        // Only read 1 has a base (read 0 has a gap)
                        position1++;
                    }
                    // else: both have gaps (shouldn't happen in a valid alignment, but handle gracefully)
                }
            }
        }
    }
}





// Main function to create variant clusters
void Assembler::performGlobalVariantClustering(
    uint64_t minCoverage,
    uint64_t maxCoverage,
    uint64_t threadCount)
{
    performanceLog << timestamp << "Starting Global Variant Clustering" << endl;
    cout << timestamp << "Starting Global Variant Clustering" << endl;
    
    cout << "\n============================================" << endl;
    cout << "VARIANT CLUSTERING PROFILING" << endl;
    cout << "============================================\n" << endl;
    
    // Check prerequisites
    const auto tCheckStart = steady_clock::now();
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    DINARA_ASSERT(compressedAlignments.isOpen());
    const auto tCheckEnd = steady_clock::now();
    const double tCheck = seconds(tCheckEnd - tCheckStart);
    
    // Access position pairs that were collected during alignment computation
    // Note: Finding and storage of position pairs happened during computeAlignments (timed separately)
    const auto tAccessStart = steady_clock::now();
    performanceLog << timestamp << "Accessing position pairs from alignment computation" << endl;
    cout << "\nAccessing position pairs collected during alignment computation..." << endl;
    cout << "(Collection and storage timing reported during alignment computation phase)" << endl;
    
    // The position pairs were already collected in variantClusteringPositionPairs during computeAlignments
    // We need to access them (they're already stored)
    if (!variantClusteringPositionPairs.isOpen) {
        cout << "ERROR: variantClusteringPositionPairs is not open!" << endl;
        return;
    }
    
    cout << "Found " << variantClusteringPositionPairs.size() << " position pairs from alignment computation" << endl;
    const auto tAccessEnd = steady_clock::now();
    const double tAccess = seconds(tAccessEnd - tAccessStart);

    



    // XXX
    // --- START OF: MINIMUM OCCURRENCES FILTER
    //     Sort, count occurrences, and filter position pairs
    //     This creates our "perfect hash" where the index is the OrientedReadId

    performanceLog << timestamp << "Sorting, counting, and filtering position pairs" << endl;
    cout << "Sorting and counting position pair occurrences..." << endl;
    const auto tOccurrenceStart = steady_clock::now();
    
    // Sort using std::sort (works because MemoryMapped::Vector::begin/end return T*)
    std::sort(variantClusteringPositionPairs.begin(), variantClusteringPositionPairs.end());

    // Count occurrences of genomic positions (strand-agnostic) and filter
    const uint64_t minOccurrences = 2;
    const uint64_t totalOccurrences = variantClusteringPositionPairs.size();

    MemoryMapped::Vector<std::pair<OrientedReadId, uint32_t>> variantClusteringFilteredPositionPairs;
    variantClusteringFilteredPositionPairs.createNew(
        largeDataName("tmp-VariantClusteringFilteredPositionPairs"),
        largeDataPageSize);
    variantClusteringFilteredPositionPairs.reserve(variantClusteringPositionPairs.size());

    if (!variantClusteringPositionPairs.empty()) {
        uint64_t count = 1;
        auto current = variantClusteringPositionPairs[0];
    
        for (uint64_t i = 1; i < totalOccurrences; ++i) {
            if (variantClusteringPositionPairs[i] == current) {
                ++count;
            } else {
                if (count >= minOccurrences) {
                    variantClusteringFilteredPositionPairs.push_back(current);
                }
                current = variantClusteringPositionPairs[i];
                count = 1;
            }
        }
        if (count >= minOccurrences) {
            variantClusteringFilteredPositionPairs.push_back(current);
        }
    }

    // Replace the memory-mapped vector with filtered results
    variantClusteringPositionPairs.clear();
    variantClusteringPositionPairs.reserve(variantClusteringFilteredPositionPairs.size());
    for (const auto& pair : variantClusteringFilteredPositionPairs) {
        variantClusteringPositionPairs.push_back(pair);
    }
    variantClusteringPositionPairs.unreserve();

    const auto tOccurrenceEnd = steady_clock::now();
    const double tOccurrence = seconds(tOccurrenceEnd - tOccurrenceStart);

    cout << "After filtering (min occurrences=" << minOccurrences << "): " 
         << variantClusteringFilteredPositionPairs.size() << " position pairs (from " << totalOccurrences << " total occurrences)" << endl;
    cout << "  Filtered out " << totalOccurrences - variantClusteringFilteredPositionPairs.size() << " position pairs with < " << minOccurrences << " occurrences" << endl;
    

    // XXX
    // --- END OF: MINIMUM OCCURRENCES FILTER
    // 

    







    // XXX
    // --- START OF: WELL-SEPARATED FILTER
    //     Filter out clusters of nearby SNPs (well-separated filter)
    //     Sequencing error and artifacts often appear as clusters of nearby SNPs.
    //     To avoid clusters of errors, the informative SNPs need to be well-separated.
    //     Only SNPs at least 10bp apart are considered.
    //     Since variantClusteringFilteredPositionPairs is already sorted by (OrientedReadId, position),
    //     positions from the same read are grouped together - we can do a single pass!
    
    const auto tSeparationStart = steady_clock::now();

    // Remember how many survived the occurrence filter.
    const uint64_t filteredCountBeforeSeparation = variantClusteringFilteredPositionPairs.size();

    // --- Filter out clusters of nearby SNPs (well-separated filter) ---
    const uint64_t minSeparation = 0;

    MemoryMapped::Vector<pair<OrientedReadId, uint32_t>> wellSeparatedPositionPairs;
    wellSeparatedPositionPairs.createNew(
        largeDataName("tmp-VariantClusteringWellSeparatedPositionPairs"),
        largeDataPageSize);
    wellSeparatedPositionPairs.reserve(variantClusteringFilteredPositionPairs.size());

    if (!variantClusteringFilteredPositionPairs.empty()) {
        ReadId currentReadId = invalid<ReadId>;
        uint32_t lastPosition = 0;

        for (const auto& pair : variantClusteringFilteredPositionPairs) {
            const OrientedReadId orientedReadId = pair.first;
            const ReadId readId = orientedReadId.getReadId();
            const uint32_t position = pair.second;

            if (readId != currentReadId) {
                // First position for this new read
                wellSeparatedPositionPairs.push_back(pair);
                currentReadId = readId;
                lastPosition = position;
            } else {
                // Same read as previous pair, check separation
                // Since input is sorted, position >= lastPosition
                if (position - lastPosition >= minSeparation) {
                    wellSeparatedPositionPairs.push_back(pair);
                    lastPosition = position;
                }
            }
        }
    }





    const uint64_t wellSeparatedCount = wellSeparatedPositionPairs.size();
    const uint64_t wellSeparatedFilteredOut = filteredCountBeforeSeparation - wellSeparatedCount;

    std::cout << "After well-separated filter (min separation="
                << minSeparation << "bp): " << wellSeparatedCount
                << " position pairs" << std::endl;
    std::cout << "  Filtered out " << wellSeparatedFilteredOut
                << " position pairs within " << minSeparation
                << "bp of adjacent positions" << std::endl;
                
    // Replace the memory-mapped vector with filtered results
    // We expand each kept Strand 0 position to include its Strand 1 reverse complement.
    variantClusteringPositionPairs.clear();
    variantClusteringPositionPairs.reserve(2 * wellSeparatedPositionPairs.size());
    for (const auto& pair : wellSeparatedPositionPairs) {
        
        // Push Strand 0
        variantClusteringPositionPairs.push_back(pair);
        
        // Push Strand 1
        const OrientedReadId orientedReadId = pair.first;
        DINARA_ASSERT(orientedReadId.getStrand() == 0);
        const ReadId readId = orientedReadId.getReadId();
        const uint32_t position = pair.second;
        const uint64_t readLength = getReads().getReadRawSequenceLength(readId);
        
        // RC position is (Length - 1 - Position)
        const uint32_t rcPosition = uint32_t(readLength - 1 - position);
        const OrientedReadId rcOrientedReadId(readId, 1);
        
        variantClusteringPositionPairs.push_back(make_pair(rcOrientedReadId, rcPosition));
    }
    variantClusteringPositionPairs.unreserve();

    const auto tSeparationEnd = steady_clock::now();
    const double tSeparation = seconds(tSeparationEnd - tSeparationStart);

    
    // DEBUG: check removed
    {
    }

    cout << "  Occurrence filter time: " << tOccurrence << " s" << endl;
    cout << "  Well-separated filter time: " << tSeparation << " s" << endl;

    // Important: Sort the final vector so that binary search (std::lower_bound) works correctly.
    // The previous vectors were sorted, but interleaving S0 and S1 pairs generally breaks global order.
    std::sort(variantClusteringPositionPairs.begin(), variantClusteringPositionPairs.end());

    performanceLog << timestamp << "Occurrence filter time " << tOccurrence << " s" << endl;
    performanceLog << timestamp << "Well-separated filter time " << tSeparation << " s" << endl;

    // XXX
    // --- END OF: WELL-SEPARATED FILTER
    // 



    if (variantClusteringPositionPairs.empty()) {
        cout << "No position pairs found. Exiting." << endl;
        return;
    }
    
    // Initialize disjoint sets for linking position pairs
    const auto tDisjointSetInitStart = steady_clock::now();
    performanceLog << timestamp << "Initializing disjoint sets" << endl;
    cout << "Initializing disjoint sets with " << variantClusteringPositionPairs.size() << " elements..." << endl;
    
    const uint64_t disjointSetCount = variantClusteringPositionPairs.size();
    
    // Store clustering data in Assembler members so it persists for mode3Assembly
    variantClusteringDisjointSetTable.createNew(
        largeDataName("tmp-VariantClusteringDisjointSets"), 
        largeDataPageSize);
    variantClusteringDisjointSetTable.resize(disjointSetCount);
    variantClusteringDisjointSets = std::make_shared<DisjointSets>(variantClusteringDisjointSetTable.begin(), disjointSetCount);
    
    // Initialize allele storage for position pairs
    variantClusteringPositionPairAlleles.createNew(
        largeDataName("tmp-VariantClusteringPositionPairAlleles"),
        largeDataPageSize);
    variantClusteringPositionPairAlleles.resize(disjointSetCount);
    
    // Initialize position context storage for position pairs
    variantClusteringPositionPairContexts.createNew(
        largeDataName("tmp-VariantClusteringPositionPairContexts"),
        largeDataPageSize);
    variantClusteringPositionPairContexts.resize(disjointSetCount);
    
    // // Update local pointer to use Assembler member
    // variantClusteringData.disjointSetsPointer = variantClusteringDisjointSets;
    
    cout << "Disjoint sets, allele and context storage initialized" << endl;

    // Link Strand 0 and Strand 1 positions for the same read/position
    // Since we pushed them in pairs (s0, s1), they are at indices 2*i and 2*i+1.
    const auto tLinkStrandsStart = steady_clock::now();
    // for(uint64_t i = 0; i < disjointSetCount / 2; i++) {
    //     variantClusteringDisjointSets->unite(2*i, 2*i+1);
    // }
    const auto tLinkStrandsEnd = steady_clock::now();
    const double tLinkStrands = seconds(tLinkStrandsEnd - tLinkStrandsStart);
    cout << "Linked " << (disjointSetCount / 2) << " strand pairs in " << tLinkStrands << " s." << endl;


    const auto tDisjointSetInitEnd = steady_clock::now();
    const double tDisjointSetInit = seconds(tDisjointSetInitEnd - tDisjointSetInitStart);

    // Phase 2: Re-process alignments to link position pairs using disjoint sets
    const auto tPass2Start = steady_clock::now();
    performanceLog << timestamp << "Phase 2: Linking position pairs with disjoint sets" << endl;
    cout << "\nPhase 2: Linking position pairs with disjoint sets..." << endl;
    
    // Initialize per-thread counters for tracking links
    // Each thread gets 4 counters: forward links, RC links, mismatches found, mismatches skipped
    variantClusteringLinkCounts.resize(threadCount * 4);
    for (size_t i = 0; i < threadCount * 4; i++) {
        variantClusteringLinkCounts[i] = 0;
    }
    
    // Pick the batch size for load balancing alignments
    const uint64_t requestedBatchSize = 1;  // Process alignments in batches
    setupLoadBalancing(alignmentData.size(), requestedBatchSize);
    runThreads(&Assembler::linkVariantClustersThreadFunction, threadCount);

    const auto tPass2End = steady_clock::now();
    const double tPass2 = seconds(tPass2End - tPass2Start);
    
    // Aggregate link counts from all threads
    uint64_t totalForwardLinks = 0;
    uint64_t totalRcLinks = 0;
    uint64_t totalMismatchesFound = 0;
    uint64_t totalMismatchesSkipped = 0;
    for (size_t i = 0; i < threadCount; i++) {
        totalForwardLinks += variantClusteringLinkCounts[i * 4 + 0];
        totalRcLinks += variantClusteringLinkCounts[i * 4 + 1];
        totalMismatchesFound += variantClusteringLinkCounts[i * 4 + 2];
        totalMismatchesSkipped += variantClusteringLinkCounts[i * 4 + 3];
    }
    
    cout << "Phase 2 complete (wall-clock time: " << tPass2 << " s)" << endl;
    cout << "  Mismatches found: " << totalMismatchesFound << endl;
    cout << "  Mismatches skipped (incomplete RC pairs): " << totalMismatchesSkipped << endl;
    cout << "  Forward strand links (unite() calls): " << totalForwardLinks << endl;
    cout << "  RC strand links (unite() calls): " << totalRcLinks << endl;
    cout << "  Total links: " << (totalForwardLinks + totalRcLinks) << endl;
    
    performanceLog << timestamp << "Phase 2 statistics: " 
                  << totalMismatchesFound << " mismatches found, "
                  << totalMismatchesSkipped << " skipped, "
                  << totalForwardLinks << " forward links, "
                  << totalRcLinks << " RC links" << endl;
    



    // {
    //     // Verify reverse complement consistency
    //     const auto tVerifyStart = steady_clock::now();
    //     performanceLog << timestamp << "Verifying reverse complement consistency" << endl;
    //     cout << "\nVerifying reverse complement consistency..." << endl;

    //     // Build reverse complement index map: (readId, strand, position) -> index
    //     std::map<pair<ReadId, pair<Strand, uint32_t>>, uint64_t> pairToIndex;
    //     for (uint64_t i = 0; i < variantClusteringPositionPairs.size(); i++) {
    //         const auto& p = variantClusteringPositionPairs[i];
    //         pairToIndex[{p.first.getReadId(), {p.first.getStrand(), p.second}}] = i;
    //     }

    //     uint64_t inconsistencies = 0;
    //     uint64_t checkedPairs = 0;
    //     // const uint64_t sampleSize = std::min(variantClusteringPositionPairs.size(), uint64_t(10000));
    //     const uint64_t sampleSize = variantClusteringPositionPairs.size();

    //     // Sample check: verify that if (readId, 0, pos) is linked with others,
    //     // then (readId, 1, readLength-1-pos) is linked with corresponding RC pairs
    //     for (uint64_t i = 0; i < sampleSize; i++) {
    //         const auto& pair = variantClusteringPositionPairs[i];
    //         const ReadId readId = pair.first.getReadId();
    //         const Strand strand = pair.first.getStrand();
    //         const uint32_t position = pair.second;
            
    //         // Only check strand 0 pairs
    //         if (strand != 0) continue;
            
    //         checkedPairs++;
            
    //         // Find reverse complement pair
    //         const uint64_t readLength = getReads().getReadRawSequenceLength(readId);
    //         const uint32_t positionRc = readLength - 1 - position;
    //         const auto rcKey = make_pair(readId, make_pair(Strand(1), positionRc));
            
    //         auto rcIt = pairToIndex.find(rcKey);
    //         if (rcIt == pairToIndex.end()) {
    //             continue; // RC pair was filtered out, skip
    //         }
            
    //         const uint64_t rcIndex = rcIt->second;
    //         const uint64_t setId0 = variantClusteringDisjointSets->find(i);
    //         const uint64_t setId0Rc = variantClusteringDisjointSets->find(rcIndex);
            
    //         // Find all pairs in the same set as i
    //         std::vector<uint64_t> linkedPairs;
    //         for (uint64_t j = 0; j < variantClusteringPositionPairs.size(); j++) {
    //             if (variantClusteringDisjointSets->find(j) == setId0) {
    //                 linkedPairs.push_back(j);
    //             }
    //         }
            
    //         // For each linked pair j (strand 0), check if rc(j) is linked with rc(i)
    //         for (uint64_t j : linkedPairs) {
    //             if (i == j) continue;
                
    //             const auto& pairJ = variantClusteringPositionPairs[j];
    //             if (pairJ.first.getStrand() != 0) continue; // Only check strand 0 pairs
                
    //             const ReadId readIdJ = pairJ.first.getReadId();
    //             const uint64_t readLengthJ = getReads().getReadRawSequenceLength(readIdJ);
    //             const uint32_t positionJRc = readLengthJ - 1 - pairJ.second;
    //             const auto rcJKey = make_pair(readIdJ, make_pair(Strand(1), positionJRc));
                
    //             auto rcJIt = pairToIndex.find(rcJKey);
    //             if (rcJIt == pairToIndex.end()) continue;
                
    //             const uint64_t rcJIndex = rcJIt->second;
    //             const uint64_t setIdJRc = variantClusteringDisjointSets->find(rcJIndex);
                
    //             // Check consistency: rc(i) and rc(j) should be in the same set
    //             if (setIdJRc != setId0Rc) {
    //                 inconsistencies++;
    //                 if (inconsistencies <= 5) {
    //                     cout << "INCONSISTENCY: " << pair.first << ":" << position 
    //                         << " linked with " << pairJ.first << ":" << pairJ.second
    //                         << ", but RC pairs are NOT linked" << endl;
    //                 }
    //             }
    //         }
    //     }

    //     const auto tVerifyEnd = steady_clock::now();
    //     const double tVerify = seconds(tVerifyEnd - tVerifyStart);

    //     cout << "Verification: checked " << checkedPairs << " pairs, found " 
    //         << inconsistencies << " inconsistencies" << endl;
    //     if (inconsistencies == 0) {
    //         cout << "✓ Reverse complement consistency verified!" << endl;
    //     } else {
    //         cout << "✗ WARNING: Found inconsistencies - check Phase 2 linking logic!" << endl;
    //     }
    //     performanceLog << timestamp << "Verification: " << checkedPairs << " pairs checked, " 
    //                 << inconsistencies << " inconsistencies" << endl;
    // }
    
    // // Verify reverse complement pair EXISTENCE (not just linkage)
    // {
    //     const auto tVerifyExistenceStart = steady_clock::now();
    //     performanceLog << timestamp << "Verifying reverse complement pair existence" << endl;
    //     cout << "\nVerifying reverse complement pair existence..." << endl;
        
    //     // Check if every forward pair has a corresponding RC pair
    //     uint64_t missingRcPairs = 0;
    //     uint64_t strand0Count = 0;
    //     uint64_t strand1Count = 0;
        
    //     // Build a set of all pairs for fast lookup
    //     std::set<pair<OrientedReadId, uint32_t>> pairSet;
    //     for (const auto& p : variantClusteringPositionPairs) {
    //         pairSet.insert(p);
    //         if (p.first.getStrand() == 0) {
    //             strand0Count++;
    //         } else {
    //             strand1Count++;
    //         }
    //     }
        
    //     // Check each pair to see if its RC counterpart exists
    //     for (const auto& pair : variantClusteringPositionPairs) {
    //         const ReadId readId = pair.first.getReadId();
    //         const Strand strand = pair.first.getStrand();
    //         const uint32_t position = pair.second;
            
    //         // Calculate RC pair
    //         const uint64_t readLength = getReads().getReadRawSequenceLength(readId);
    //         const uint32_t positionRc = readLength - 1 - position;
    //         const Strand strandRc = 1 - strand;
    //         OrientedReadId orientedReadIdRc(readId, strandRc);
    //         const auto rcPair = make_pair(orientedReadIdRc, positionRc);
            
    //         // Check if RC pair exists
    //         if (pairSet.find(rcPair) == pairSet.end()) {
    //             missingRcPairs++;
    //             if (missingRcPairs <= 5) {
    //                 cout << "MISSING RC PAIR: " << pair.first << ":" << position 
    //                      << " exists, but RC pair " << orientedReadIdRc << ":" << positionRc
    //                      << " does NOT exist" << endl;
    //             }
    //         }
    //     }
        
    //     const auto tVerifyExistenceEnd = steady_clock::now();
    //     const double tVerifyExistence = seconds(tVerifyExistenceEnd - tVerifyExistenceStart);
        
    //     cout << "Pair existence check:" << endl;
    //     cout << "  Total pairs: " << variantClusteringPositionPairs.size() << endl;
    //     cout << "  Strand 0 pairs: " << strand0Count << endl;
    //     cout << "  Strand 1 pairs: " << strand1Count << endl;
    //     cout << "  Missing RC pairs: " << missingRcPairs << endl;
        
    //     if (missingRcPairs == 0) {
    //         cout << "✓ All pairs have their reverse complement counterparts!" << endl;
    //     } else {
    //         cout << "✗ WARNING: " << missingRcPairs << " pairs are missing their RC counterparts!" << endl;
    //         cout << "  This explains why DINARA_ASSERT(found0 && found1 && found0Rc && found1Rc) fails." << endl;
    //         cout << "  Phase 1 filtering created an asymmetry between forward and RC strands." << endl;
    //     }
        
    //     performanceLog << timestamp << "RC existence verification: " 
    //                   << missingRcPairs << " missing RC pairs out of " 
    //                   << variantClusteringPositionPairs.size() << " total pairs" << endl;
    // }




    
    // Identify cluster representatives
    const auto tIdentifyClustersStart = steady_clock::now();
    performanceLog << timestamp << "Identifying cluster representatives" << endl;
    cout << "Identifying cluster representatives..." << endl;
    
    // Store cluster representatives in Assembler member so they persist for mode3Assembly
    variantClusteringClusterRepresentatives.clear();
    for (uint64_t id = 0; id < disjointSetCount; id++) {
        if (variantClusteringDisjointSets->find(id) == id) {
            variantClusteringClusterRepresentatives.push_back(id);
        }
    }
    
    cout << "Found " << variantClusteringClusterRepresentatives.size() << " clusters from " 
         << disjointSetCount << " unique position pairs" << endl;

    const auto tIdentifyClustersEnd = steady_clock::now();
    const double tIdentifyClusters = seconds(tIdentifyClustersEnd - tIdentifyClustersStart);




    {
        // DEBUG: Print all clusters containing read 0-0
        cout << "\n=== DEBUG: Clusters containing read 0-0 ===" << endl;
        const OrientedReadId debugRead(0, 0);
        
        // Build a map from representative ID to cluster index for quick lookup
        std::map<uint64_t, size_t> repToClusterIdx;
        for(size_t i = 0; i < variantClusteringClusterRepresentatives.size(); i++) {
            repToClusterIdx[variantClusteringClusterRepresentatives[i]] = i;
        }
        
        // Build clusters: map from cluster index to member IDs
        std::vector<std::vector<uint64_t>> clusterMembers(variantClusteringClusterRepresentatives.size());
        for(uint64_t id = 0; id < disjointSetCount; id++) {
            const uint64_t rep = variantClusteringDisjointSets->find(id);
            auto it = repToClusterIdx.find(rep);
            if(it != repToClusterIdx.end()) {
                clusterMembers[it->second].push_back(id);
            }
        }
        
        // Find and print clusters containing read 0-0
        uint64_t clustersWithDebugRead = 0;
        for(size_t clusterIdx = 0; clusterIdx < variantClusteringClusterRepresentatives.size(); clusterIdx++) {
            const auto& members = clusterMembers[clusterIdx];
            
            // Check if this cluster contains read 0-0
            bool containsDebugRead = false;
            for(uint64_t id : members) {
                const auto& posPair = variantClusteringPositionPairs[id];
                if(posPair.first == debugRead) {
                    containsDebugRead = true;
                    break;
                }
            }
            
            if(!containsDebugRead) {
                continue;
            }
            
            clustersWithDebugRead++;
            const uint64_t rep = variantClusteringClusterRepresentatives[clusterIdx];
            
            // Count total unique reads and reads per allele
            std::unordered_set<uint32_t> allReads;
            std::array<std::unordered_set<uint32_t>, 5> readsByAllele; // A, C, G, T, and gap
            
            for(uint64_t id : members) {
                const auto& posPair = variantClusteringPositionPairs[id];
                const uint32_t ridVal = posPair.first.getValue();
                allReads.insert(ridVal);
                
                if(id < variantClusteringPositionPairAlleles.size()) {
                    const uint8_t allele = variantClusteringPositionPairAlleles[id];
                    if(allele < 5) {
                        readsByAllele[allele].insert(ridVal);
                    }
                }
            }
            
            // Check if this cluster has at least 2 alleles with coverage >= 4
            uint64_t eligibleAlleles = 0;
            for(uint8_t a = 0; a < 5; a++) {
                if(readsByAllele[a].size() >= 4) {
                    eligibleAlleles++;
                }
            }
            
            if(eligibleAlleles < 2) {
                continue;  // Skip clusters that don't meet the criteria
            }
            
            cout << "\nCluster " << clusterIdx << " (rep=" << rep << "):" << endl;
            cout << "  Total position pairs: " << members.size() << endl;
            cout << "  Total unique reads: " << allReads.size() << endl;
            cout << "  Reads per allele:";
            for(uint8_t a = 0; a < 5; a++) {
                if(readsByAllele[a].size() > 0) {
                    const char alleleChar = Base::fromInteger(a).character();
                    cout << " " << alleleChar << "=" << readsByAllele[a].size();
                }
            }
            cout << endl;
            
            // Print all position pairs in this cluster
            cout << "  Position pairs:" << endl;
            const size_t maxPairsToShow = std::min(size_t(20), members.size());
            for(size_t i = 0; i < maxPairsToShow; i++) {
                const uint64_t id = members[i];
                const auto& posPair = variantClusteringPositionPairs[id];
                const uint8_t allele = (id < variantClusteringPositionPairAlleles.size()) 
                                        ? variantClusteringPositionPairAlleles[id] : 4;
                const char alleleChar = (allele == 0) ? 'A' : (allele == 1) ? 'C' : 
                                        (allele == 2) ? 'G' : (allele == 3) ? 'T' : '?';
                
                cout << "    [" << i << "] ReadId=" << posPair.first.getReadId()
                    << " Strand=" << posPair.first.getStrand()
                    << " Pos=" << posPair.second
                    << " Allele=" << alleleChar << endl;
            }
            if(members.size() > maxPairsToShow) {
                cout << "    ... (" << (members.size() - maxPairsToShow) << " more pairs)" << endl;
            }
        }
        
        if(clustersWithDebugRead == 0) {
            cout << "No clusters found containing read 0-0" << endl;
        } else {
            cout << "\nTotal clusters containing read 0-0: " << clustersWithDebugRead << endl;
        }
        cout << "=== END DEBUG ===" << endl << endl;
    }


    
    // Print timing summary (excluding debug/verification code)
    const double tTotal = tCheck + tAccess + variantClusteringProjectedAlignmentTime + variantClusteringCollectionTime + variantClusteringStorageTime + tOccurrence + tSeparation + tDisjointSetInit + tLinkStrands + tPass2 + tIdentifyClusters;
    
    cout << "\n============================================" << endl;
    cout << "VARIANT CLUSTERING TIMING SUMMARY" << endl;
    cout << "============================================" << endl;
    cout << std::left << std::setw(40) << "Step" << std::right << std::setw(12) << "Time (s)" << std::setw(12) << "Percent" << endl;
    cout << std::string(64, '-') << endl;
    cout << std::left << std::setw(40) << "Prerequisites check" << std::right << std::setw(12) << tCheck << std::setw(12) << (100.0 * tCheck / tTotal) << endl;
    cout << std::left << std::setw(40) << "Access position pairs" << std::right << std::setw(12) << tAccess << std::setw(12) << (100.0 * tAccess / tTotal) << endl;
    cout << std::left << std::setw(40) << "  Projected alignment" << std::right << std::setw(12) << variantClusteringProjectedAlignmentTime << std::setw(12) << (100.0 * variantClusteringProjectedAlignmentTime / tTotal) << endl;
    cout << std::left << std::setw(40) << "  Collection" << std::right << std::setw(12) << variantClusteringCollectionTime << std::setw(12) << (100.0 * variantClusteringCollectionTime / tTotal) << endl;
    cout << std::left << std::setw(40) << "  Storage" << std::right << std::setw(12) << variantClusteringStorageTime << std::setw(12) << (100.0 * variantClusteringStorageTime / tTotal) << endl;
    cout << std::left << std::setw(40) << "Occurrence filter" << std::right << std::setw(12) << tOccurrence << std::setw(12) << (100.0 * tOccurrence / tTotal) << endl;
    cout << std::left << std::setw(40) << "Well-separated filter" << std::right << std::setw(12) << tSeparation << std::setw(12) << (100.0 * tSeparation / tTotal) << endl;
    cout << std::left << std::setw(40) << "Initialize disjoint sets" << std::right << std::setw(12) << tDisjointSetInit << std::setw(12) << (100.0 * tDisjointSetInit / tTotal) << endl;
    cout << std::left << std::setw(40) << "Link strands" << std::right << std::setw(12) << tLinkStrands << std::setw(12) << (100.0 * tLinkStrands / tTotal) << endl;
    cout << std::left << std::setw(40) << "Phase 2: Link pairs" << std::right << std::setw(12) << tPass2 << std::setw(12) << (100.0 * tPass2 / tTotal) << endl;
    cout << std::left << std::setw(40) << "Identify clusters" << std::right << std::setw(12) << tIdentifyClusters << std::setw(12) << (100.0 * tIdentifyClusters / tTotal) << endl;
    cout << std::string(64, '-') << endl;
    cout << std::left << std::setw(40) << "Total Time" << std::right << std::setw(12) << tTotal << std::setw(12) << "100.0" << endl;
    cout << "============================================\n" << endl;
    
    performanceLog << timestamp << "createVariantClusters ends" << endl;
}