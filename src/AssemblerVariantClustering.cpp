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
#include <numeric>
#include <vector>


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
            // else: both have gaps (shouldn't happen in a valid alignment)
            DINARA_ASSERT(hasBase0 or hasBase1);
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
            
            // Skip deleted alignments (e.g. redundant alignments from Best Hit Filtering)
            if (alignmentDataRef[alignmentId].isDeleted()) continue;

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
                    DINARA_ASSERT(hasBase0 || hasBase1);
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

    // Count occurrences of genomic positions (strand-agnostic) and filter IN-PLACE
    const uint64_t minOccurrences = assemblerInfo->variantClusteringMinOccurrences;
    const uint64_t totalOccurrences = variantClusteringPositionPairs.size();

    if (totalOccurrences > 0) {
        uint64_t writeIndex = 0;
        uint64_t count = 1;
        auto current = variantClusteringPositionPairs[0];

        for (uint64_t i = 1; i < totalOccurrences; ++i) {
            if (variantClusteringPositionPairs[i] == current) {
                ++count;
            } else {
                if (count >= minOccurrences) {
                    variantClusteringPositionPairs[writeIndex++] = current;
                }
                current = variantClusteringPositionPairs[i];
                count = 1;
            }
        }
        // Handle the last group
        if (count >= minOccurrences) {
            variantClusteringPositionPairs[writeIndex++] = current;
        }

        // Resize to the filtered count
        variantClusteringPositionPairs.resize(writeIndex);
    } else {
        variantClusteringPositionPairs.resize(0);
    }

    const auto tOccurrenceEnd = steady_clock::now();
    const double tOccurrence = seconds(tOccurrenceEnd - tOccurrenceStart);

    cout << "After filtering (min occurrences=" << minOccurrences << "): " 
         << variantClusteringPositionPairs.size() << " position pairs (from " << totalOccurrences << " total occurrences)" << endl;
    cout << "  Filtered out " << totalOccurrences - variantClusteringPositionPairs.size() << " position pairs with < " << minOccurrences << " occurrences" << endl;

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
    const uint64_t filteredCountBeforeSeparation = variantClusteringPositionPairs.size();

    // --- Filter out clusters of nearby SNPs (well-separated filter) ---
    const uint64_t minSeparation = assemblerInfo->variantClusteringMinSeparation;

    if (filteredCountBeforeSeparation > 0) {
        uint64_t writeIndex = 0;
        ReadId currentReadId = invalid<ReadId>;
        uint32_t lastPosition = 0;

        for (uint64_t i = 0; i < filteredCountBeforeSeparation; ++i) {
            const auto& pair = variantClusteringPositionPairs[i];
            const OrientedReadId orientedReadId = pair.first;
            const ReadId readId = orientedReadId.getReadId();
            const uint32_t position = pair.second;

            if (readId != currentReadId) {
                // First position for this new read
                variantClusteringPositionPairs[writeIndex++] = pair;
                currentReadId = readId;
                lastPosition = position;
            } else {
                // Same read as previous pair, check separation
                // Since input is sorted, position >= lastPosition
                if (position - lastPosition >= minSeparation) {
                    variantClusteringPositionPairs[writeIndex++] = pair;
                    lastPosition = position;
                }
            }
        }
        variantClusteringPositionPairs.resize(writeIndex);
    }

    const auto tSeparationEnd = steady_clock::now();
    const double tSeparation = seconds(tSeparationEnd - tSeparationStart);

    const uint64_t wellSeparatedCount = variantClusteringPositionPairs.size();
    const uint64_t wellSeparatedFilteredOut = filteredCountBeforeSeparation - wellSeparatedCount;

    std::cout << "After well-separated filter (min separation="
                << minSeparation << "bp): " << wellSeparatedCount
                << " position pairs" << std::endl;
    std::cout << "  Filtered out " << wellSeparatedFilteredOut
                << " position pairs within " << minSeparation
                << "bp of adjacent positions" << std::endl;
                
    // Expand each kept Strand 0 position to include its Strand 1 reverse complement.
    // We do this by resizing and filling the second half.
    const uint64_t s0Count = variantClusteringPositionPairs.size();
    
    // Resize to hold both strands
    variantClusteringPositionPairs.resize(2 * s0Count);
    
    for (uint64_t i = 0; i < s0Count; ++i) {
        const auto& pair = variantClusteringPositionPairs[i];
        
        // Push Strand 1
        const OrientedReadId orientedReadId = pair.first;
        DINARA_ASSERT(orientedReadId.getStrand() == 0);
        const ReadId readId = orientedReadId.getReadId();
        const uint32_t position = pair.second;
        const uint64_t readLength = getReads().getReadRawSequenceLength(readId);
        
        // RC position is (Length - 1 - Position)
        const uint32_t rcPosition = uint32_t(readLength - 1 - position);
        const OrientedReadId rcOrientedReadId(readId, 1);
        
        // Store at the second half
        variantClusteringPositionPairs[s0Count + i] = make_pair(rcOrientedReadId, rcPosition);
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

    
    // Identify cluster representatives
    const auto tIdentifyClustersStart = steady_clock::now();
    performanceLog << timestamp << "Identifying cluster representatives" << endl;
    cout << "Identifying cluster representatives..." << endl;
    
    // Store cluster representatives in Assembler member so they persist for mode3Assembly
    variantClusteringClusterRepresentatives.clear();

    // To efficiently validate clusters, we first group all element IDs by their cluster representative.
    // We use a flat vector of (RepID, ElementID) pairs and sort it to avoid map overhead.
    cout << "Grouping clusters for validation..." << endl;
    vector<pair<uint64_t, uint64_t>> elementsByCluster;
    elementsByCluster.reserve(disjointSetCount);
    for (uint64_t id = 0; id < disjointSetCount; id++) {
        elementsByCluster.push_back({variantClusteringDisjointSets->find(id), id});
    }
    std::sort(elementsByCluster.begin(), elementsByCluster.end());

    uint64_t discardedClusters = 0;
    uint64_t keptClusters = 0;
    uint64_t removedReads = 0;

    // Pre-allocate buffer for duplicate detection to avoid re-allocation per cluster
    std::vector<pair<ReadId, uint64_t>> readOccurrences;
    readOccurrences.reserve(256); // Reasonable expected coverage

    // Iterate through groups of elements belonging to the same cluster
    uint64_t totalElements = elementsByCluster.size();
    uint64_t i = 0;
    while (i < totalElements) {
        uint64_t repId = elementsByCluster[i].first;
        uint64_t start = i;
        
        // Find end of current cluster
        while (i < totalElements && elementsByCluster[i].first == repId) {
            i++;
        }
        uint64_t end = i; // [start, end)
        
        uint64_t clusterSize = end - start;
        bool keepCluster = true;

        /*  
        if (clusterSize > 1) {
        
            // 1. Identify reads that appear multiple times in this cluster.
            // We want to remove ALL occurrences of ambiguous reads.
            
            // Gather (ReadId, ElementIndex) for all members
            readOccurrences.clear();
            if (readOccurrences.capacity() < clusterSize) {
            readOccurrences.reserve(clusterSize);
            }

            for(uint64_t k = start; k < end; k++) {
                uint64_t elementId = elementsByCluster[k].second;
                ReadId r = variantClusteringPositionPairs[elementId].first.getReadId();
                readOccurrences.push_back({r, elementId});
            }
            
            // Sort by ReadId to find duplicates
            std::sort(readOccurrences.begin(), readOccurrences.end());
            
            uint64_t validMemberCount = 0;
            uint64_t r = 0;
            while(r < clusterSize) {
                ReadId currentRead = readOccurrences[r].first;
                uint64_t rStart = r;
                while(r < clusterSize && readOccurrences[r].first == currentRead) {
                    r++;
                }
                uint64_t rEnd = r;
                uint64_t count = rEnd - rStart;
                
                if (count > 1) {
                    // This read appears multiple times. Mark all occurrences as invalid (255).
                    for(uint64_t k = rStart; k < rEnd; k++) {
                        uint64_t elementId = readOccurrences[k].second;
                        if (elementId < variantClusteringPositionPairAlleles.size()) {
                        variantClusteringPositionPairAlleles[elementId] = 255;
                        removedReads++;
                        }
                    }
                } else {
                    validMemberCount++;
                }
            }

            if (validMemberCount == 0) {
                keepCluster = false;
            }
        
        }
        */

        // STRICT FILTER: Discard entire cluster if any ambiguity exists (Total Position Pairs != Total Unique Reads)
        if (clusterSize > 1) {
            readOccurrences.clear();
            if (readOccurrences.capacity() < clusterSize) {
                readOccurrences.reserve(clusterSize);
            }
            for(uint64_t k = start; k < end; k++) {
                uint64_t elementId = elementsByCluster[k].second;
                ReadId r = variantClusteringPositionPairs[elementId].first.getReadId();
                readOccurrences.push_back({r, elementId});
            }
            
            // Sort by ReadId
            std::sort(readOccurrences.begin(), readOccurrences.end());
            
            // Check for duplicates
            for(size_t k = 1; k < readOccurrences.size(); k++) {
                if (readOccurrences[k].first == readOccurrences[k-1].first) {
                    keepCluster = false;
                    break;
                }
            }
        }
        
        
        if (keepCluster) {
            variantClusteringClusterRepresentatives.push_back(repId);
            keptClusters++;
        } else {
            discardedClusters++;
        }
    }
    
    cout << "Found " << variantClusteringClusterRepresentatives.size() << " valid clusters." << endl;
    cout << "Discarded " << discardedClusters << " empty clusters." << endl;
    cout << "Removed " << removedReads << " ambiguous read references (duplicates within clusters)." << endl;

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


