#include "Assembler.hpp"
#include "timestamp.hpp"
#include "deduplicate.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <mutex>
#include <set>
#include "Reads.hpp"
#include <numeric>
#include <algorithm>
#include <vector>

using namespace dinara;



bool isSiteInHomopolymerRegion(
    uint64_t sitePositionInRead,
    const dinara::LongBaseSequenceView& readSequence,
    uint64_t homopolymerContextLength = 12,  // HPC_PL (homopolymer flank length)
    uint64_t maxRepeatLength = 4,            // HPC_RR (homopolymer repeat range)
    uint64_t homopolymerCutoff = 2           // HPC_CC (homopolymer cutoff)
) {
    const uint64_t readLength = readSequence.baseCount;
    DINARA_ASSERT(sitePositionInRead < readLength);

    // Convert to int64_t for compatibility with hifiasm logic
    int64_t p = static_cast<int64_t>(sitePositionInRead);
    int64_t sn = static_cast<int64_t>(readLength);
    int64_t hpc_flk = static_cast<int64_t>(homopolymerContextLength);
    int64_t hpc_rr = static_cast<int64_t>(maxRepeatLength);
    int64_t hpc_cutoff = static_cast<int64_t>(homopolymerCutoff);

    // Calculate boundaries (same as hifiasm)
    int64_t s = ((p >= hpc_flk) ? (p - hpc_flk) : 0);
    int64_t e = (((p + hpc_flk) <= sn) ? (p + hpc_flk) : sn);

    // Safety checks
    if (s < 0) s = 0;
    if (e > sn) e = sn;
    DINARA_ASSERT(s < e);

    for (int64_t r = 1; r <= hpc_rr; r++) {
        int64_t rc = r * hpc_cutoff;
        
        // Pattern 1: Including site (forward and backward)
        int64_t k, zs, ze;
        
        // Forward scan from position p + r
        for (k = p + r; (k < e) && ((k - r) >= s) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k - r)].value); k++);
        ze = k;
        if (ze > e) ze = e;
        
        // Backward scan from position p - 1
        for (k = p - 1; (k >= s) && ((k + r) < e) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k + r)].value); k--);
        zs = k + 1;
        if (zs < s) zs = s;
        
        if (((ze - zs) > r) && ((ze - zs) >= rc)) {
            return true;
        }

        // Pattern 2: Excluding site (forward only)
        for (k = p + r + 1; (k < e) && ((k - r) >= s) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k - r)].value); k++);
        zs = p + 1;
        if (zs < s) zs = s;
        ze = k;
        if (ze > e) ze = e;
        
        if (((ze - zs) > r) && ((ze - zs) >= rc)) {
            return true;
        }

        // Pattern 3: Including site (backward and forward)
        for (k = p - r; (k >= s) && ((k + r) < e) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k + r)].value); k--);
        zs = k + 1;
        if (zs < s) zs = s;
        
        for (k = p + 1; (k < e) && ((k - r) >= s) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k - r)].value); k++);
        ze = k;
        if (ze > e) ze = e;
        
        if (((ze - zs) > r) && ((ze - zs) >= rc)) {
            return true;
        }

        // Pattern 4: Excluding site (backward only)
        for (k = p - r - 1; (k >= s) && ((k + r) < e) && 
             (readSequence[static_cast<uint64_t>(k)].value == readSequence[static_cast<uint64_t>(k + r)].value); k--);
        zs = k + 1;
        if (zs < s) zs = s;
        ze = p;
        if (ze > e) ze = e;
        
        if (((ze - zs) > r) && ((ze - zs) >= rc)) {
            return true;
        }
    }
    
    return false;
}



// -----------------------------------------------------------------------------
// Overload: isSiteInHomopolymerRegion accepting a vector<Base> sequence
// -----------------------------------------------------------------------------
inline bool isSiteInHomopolymerRegion(
    uint64_t sitePositionInRead,
    const std::vector<dinara::Base>& readSequence,
    uint64_t homopolymerContextLength = 12,  // HPC_PL (homopolymer flank length)
    uint64_t maxRepeatLength        = 4,     // HPC_RR (homopolymer repeat range)
    uint64_t homopolymerCutoff      = 2)     // HPC_CC (homopolymer cutoff)
{
    const uint64_t readLength = readSequence.size();
    DINARA_ASSERT(sitePositionInRead < readLength);

    // Convert to int64_t for computation identical to the other overload.
    int64_t p          = static_cast<int64_t>(sitePositionInRead);
    int64_t sn         = static_cast<int64_t>(readLength);
    int64_t hpc_flk    = static_cast<int64_t>(homopolymerContextLength);
    int64_t hpc_rr     = static_cast<int64_t>(maxRepeatLength);
    int64_t hpc_cutoff = static_cast<int64_t>(homopolymerCutoff);

    int64_t s = (p >= hpc_flk) ? (p - hpc_flk) : 0;
    int64_t e = ((p + hpc_flk) <= sn) ? (p + hpc_flk) : sn;
    if(s < 0) s = 0;
    if(e > sn) e = sn;
    DINARA_ASSERT(s < e);

    for(int64_t r = 1; r <= hpc_rr; ++r) {
        int64_t rc = r * hpc_cutoff;
        int64_t k, zs, ze;

        // Pattern 1: Including site (forward & backward)
        for(k = p + r; (k < e) && ((k - r) >= s) && (readSequence[k].value == readSequence[k - r].value); ++k) {}
        ze = (k > e) ? e : k;
        for(k = p - 1; (k >= s) && ((k + r) < e) && (readSequence[k].value == readSequence[k + r].value); --k) {}
        zs = k + 1;
        if(zs < s) zs = s;
        if(((ze - zs) > r) && ((ze - zs) >= rc)) return true;

        // Pattern 2: Excluding site (forward only)
        for(k = p + r + 1; (k < e) && ((k - r) >= s) && (readSequence[k].value == readSequence[k - r].value); ++k) {}
        zs = p + 1;
        if(zs < s) zs = s;
        ze = (k > e) ? e : k;
        if(((ze - zs) > r) && ((ze - zs) >= rc)) return true;

        // Pattern 3: Including site (backward & forward)
        for(k = p - r; (k >= s) && ((k + r) < e) && (readSequence[k].value == readSequence[k + r].value); --k) {}
        zs = k + 1;
        if(zs < s) zs = s;
        for(k = p + 1; (k < e) && ((k - r) >= s) && (readSequence[k].value == readSequence[k - r].value); ++k) {}
        ze = (k > e) ? e : k;
        if(((ze - zs) > r) && ((ze - zs) >= rc)) return true;

        // Pattern 4: Excluding site (backward only)
        for(k = p - r - 1; (k >= s) && ((k + r) < e) && (readSequence[k].value == readSequence[k + r].value); --k) {}
        zs = k + 1;
        if(zs < s) zs = s;
        ze = p;
        if(ze > e) ze = e;
        if(((ze - zs) > r) && ((ze - zs) >= rc)) return true;
    }

    return false;
}


// Function to check for a specific strand bias pattern
// Returns true if the strand bias is detected, and populates outDominantStrandHet1 and outDominantStrandHet2.
// dominantStrand values: 0 (all on strand 0), 1 (all on strand 1), 2 (no dominance or not enough reads).
bool hasSiteStrandBiasReadGraph5(
    const std::set<OrientedReadId>& hetBase1OrientedReads,
    const std::set<OrientedReadId>& hetBase2OrientedReads,
    uint64_t& outDominantStrandHet1,
    uint64_t& outDominantStrandHet2)
{
    // For this specific bias pattern, we need at least one read supporting each allele's strand pattern.
    const uint64_t minReadsPerAlleleForThisBias = 1;

    outDominantStrandHet1 = 2; // Initialize to no dominance
    if (hetBase1OrientedReads.size() >= minReadsPerAlleleForThisBias) {
        uint64_t hetBase1CountStrand0 = 0;
        uint64_t hetBase1CountStrand1 = 0;
        for (const auto& hetBase1OrientedRead : hetBase1OrientedReads) {
            if (hetBase1OrientedRead.getStrand() == 0) {
                hetBase1CountStrand0++;
            } else {
                hetBase1CountStrand1++;
            }
        }
        if (hetBase1CountStrand0 > 0 && hetBase1CountStrand1 == 0) {
            outDominantStrandHet1 = 0; // All reads for hetBase1 are on strand 0
        } else if (hetBase1CountStrand1 > 0 && hetBase1CountStrand0 == 0) {
            outDominantStrandHet1 = 1; // All reads for hetBase1 are on strand 1
        }
    }

    outDominantStrandHet2 = 2; // Initialize to no dominance
    if (hetBase2OrientedReads.size() >= minReadsPerAlleleForThisBias) {
        uint64_t hetBase2CountStrand0 = 0;
        uint64_t hetBase2CountStrand1 = 0;
        for (const auto& hetBase2OrientedRead : hetBase2OrientedReads) {
            if (hetBase2OrientedRead.getStrand() == 0) {
                hetBase2CountStrand0++;
            } else {
                hetBase2CountStrand1++;
            }
        }
        if (hetBase2CountStrand0 > 0 && hetBase2CountStrand1 == 0) {
            outDominantStrandHet2 = 0; // All reads for hetBase2 are on strand 0
        } else if (hetBase2CountStrand1 > 0 && hetBase2CountStrand0 == 0) {
            outDominantStrandHet2 = 1; // All reads for hetBase2 are on strand 1
        }
    }

    // Check for the specific bias pattern:
    // Both alleles must show complete strand dominance, and on opposite strands.
    if (outDominantStrandHet1 != 2 && outDominantStrandHet2 != 2 && outDominantStrandHet1 != outDominantStrandHet2) {
        return true; // Bias detected
    }

    return false; // No bias detected
}







void Assembler::createReadGraph5()
{
    cout << timestamp << "createReadGraph5 begins." << endl;

    // Check that we have the necessary data.
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    checkVariantClusteringPositionPairsIsOpen();

    if(!variantClusteringDisjointSets) {
        throw runtime_error("Variant clustering disjoint sets not available.");
    }

    // 1. Build membersByRepIdx (Cluster -> PositionPairs map)
    // We use the 2-pass approach for MemoryMapped::VectorOfVectors to be memory efficient
    // and avoid creating a huge temporary vector<vector> in RAM.
    cout << timestamp << "Building cluster to reads map." << endl;
    
    const uint64_t positionPairCount = variantClusteringPositionPairs.size();
    const uint64_t clusterCount = variantClusteringDisjointSets->size();

    variantClusteringMembersByRepIdx.createNew(
        largeDataName("VariantClusteringMembersByRepIdx"),
        largeDataPageSize
    );

    // Pass 1: Count the number of members in each cluster (representative).
    // This prepares the "Table of Contents" (TOC).
    variantClusteringMembersByRepIdx.beginPass1(clusterCount);
    for(uint64_t i=0; i<positionPairCount; i++) {
        const uint64_t clusterId = variantClusteringDisjointSets->find(i);
        variantClusteringMembersByRepIdx.incrementCount(clusterId);
    }

    // Pass 2: Store the member IDs.
    // This fills the data vector.
    variantClusteringMembersByRepIdx.beginPass2();
    for(uint64_t i=0; i<positionPairCount; i++) {
        const uint64_t clusterId = variantClusteringDisjointSets->find(i);
        variantClusteringMembersByRepIdx.store(clusterId, i);
    }
    variantClusteringMembersByRepIdx.endPass2();
    
    cout << timestamp << "Map built. " << clusterCount << " potential clusters." << endl;

    // 2. Initialize output data structures.
    // It will be used as a lookup table (or boolean array) to mark valid clusters.
    variantClusteringValidClusters.createNew(
        largeDataName("VariantClusteringValidClusters"),
        largeDataPageSize);
    variantClusteringValidClusters.resize(clusterCount);
    // Initialize to 0 (invalid)
    std::fill(variantClusteringValidClusters.begin(), variantClusteringValidClusters.end(), 0);

    // 3. Run threads to perform global clustervalidity check.
    // A cluster is valid if it has at least 2 alleles with coverage >= minAlleleCoverage.
    cout << timestamp << "Running global cluster validity checks." << endl;
    const uint32_t threadCount = std::thread::hardware_concurrency();
    const uint64_t minAlleleCoverage = 5;
    setupLoadBalancing(clusterCount, 1);
    runThreads(&Assembler::computeClusterValidityThreadFunction, threadCount);
    cout << timestamp << "Global cluster validity checks completed." << endl;


    // 4. Prepare for Haplotype Voting
    variantClusteringValidClustersCompatible.createNew(
        largeDataName("VariantClusteringValidClustersCompatible"),
        largeDataPageSize
    );
    variantClusteringValidClustersCompatible.resize(clusterCount);
    std::fill(variantClusteringValidClustersCompatible.begin(), variantClusteringValidClustersCompatible.end(), 0);

    // 4. Initialize Global Haplotype Graph
    // Calculate total oriented reads (2 * readCount)
    const uint64_t orientedReadCount = 2 * reads->readCount();
    
    // Create the graph with orientedReadCount vertices
    cout << timestamp << "Initializing Global Haplotype Graph with " << orientedReadCount << " vertices." << endl;
    globalHaplotypeGraph = std::make_shared<HaplotypeGraph>(orientedReadCount);

    // 5. Run threads to generate votes and add edges directly
    cout << timestamp << "Running compatibility checks and building haplotype graph." << endl;
    setupLoadBalancing(getReads().readCount(), 1); 
    runThreads(&Assembler::createReadGraph5ThreadFunction, threadCount);
    cout << timestamp << "Haplotype graph construction completed." << endl;
    
    cout << timestamp << "Global Haplotype Graph built with " 
         << boost::num_vertices(*globalHaplotypeGraph) << " vertices and "
         << boost::num_edges(*globalHaplotypeGraph) << " edges." << endl;

    // 7. Refine Clusters using the Graph
    cout << timestamp << "Refining clusters using haplotype graph." << endl;
    
    // This vector is used to track the quality or "membership status" of each read within a variant cluster.
    // Here is the breakdown:
    // 1. Initialization: It creates a new memory-mapped vector sized to positionPairCount (the total number of read-variant pairs).
    // 2. Default State: It fills the vector with 0, which represents a "Good" or "Keep" status.
    // 3. Usage: Later in the refineClustersThreadFunction this vector is updated. If a read is identified 
    // as "stray" (meaning it connects more strongly to reads of a different allele in the haplotype graph than to its own), 
    // then its status is changed to 1 (Filter/Stray).
    // Essentially, this prepares the "report card" for every read in every cluster, assuming they are all valid until 
    // proven otherwise by the graph refinement step.
    variantClusteringMemberStatus.createNew(
        largeDataName("VariantClusteringMemberStatus"),
        largeDataPageSize
    );
    variantClusteringMemberStatus.resize(positionPairCount);
    std::fill(variantClusteringMemberStatus.begin(), variantClusteringMemberStatus.end(), 0); // Default 0 = Good

    setupLoadBalancing(clusterCount, 1);
    runThreads(&Assembler::refineClustersThreadFunction, threadCount);
    
    cout << timestamp << "Cluster refinement completed." << endl;

    // Print how many clusters are valid and how many are compatible.
    uint64_t validClusterCount = 0;
    uint64_t compatibleClusterCount = 0;
    for(uint64_t i=0; i<clusterCount; i++) {
        if(variantClusteringValidClusters[i]) validClusterCount++;
        if(variantClusteringValidClustersCompatible[i]) compatibleClusterCount++;
    }
    cout << timestamp << "Valid clusters: " << validClusterCount << endl;
    cout << timestamp << "Compatible clusters: " << compatibleClusterCount << endl;

    // Filter alignments based on the Final Haplotype Graph
    cout << timestamp << "Filtering alignments using Global Haplotype Graph compatibility." << endl;

    // Get the total number of stored alignments.
    const uint64_t alignmentCount = alignmentData.size();
    DINARA_ASSERT(compressedAlignments.size() == alignmentCount);

    // Flag alignments to be kept.
    // "Initially keep all alignments and for each read start removing its alignments that do not exist in the globalHaplotypeGraph"
    // "If the orientedReadIds ... are not in the haplotypeGraph (meaning they did not have any het site for phasing) ... include those alignments"
    // Logic: VETO. If a read is Phased (has out-edges), it provides an exclusive whitelist of neighbors.
    // If an alignment is not in that whitelist, it is removed.
    
    vector<bool> keepAlignment(alignmentCount, true);
    uint64_t keptAlignmentCount = 0;

    const auto& graph = *globalHaplotypeGraph;
    const size_t numVertices = boost::num_vertices(graph);

    // Parallelize filtering for efficiency
    #pragma omp parallel for reduction(+:keptAlignmentCount)
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
         const auto& ad = alignmentData[alignmentId];
         
         bool veto = false;

         // Check Veto from Read 0 Perspective
         // We only check Strand 0 because we only populated out-edges for Strand 0 in createReadGraph5ThreadFunction.
         OrientedReadId::Int u1 = OrientedReadId(ad.readIds[0], 0).getValue();
         if (u1 < numVertices && boost::out_degree(u1, graph) > 0) {
             // Read 0 is phased. It MUST have an edge to Read 1 with correct orientation.
             OrientedReadId::Int v_target = OrientedReadId(ad.readIds[1], ad.isSameStrand ? 0 : 1).getValue();
             
             bool found = false;
             auto outEdges = boost::out_edges(u1, graph);
             for(auto it = outEdges.first; it != outEdges.second; ++it) {
                 if (boost::target(*it, graph) == v_target) {
                     found = true;
                     break;
                 }
             }
             if (!found) veto = true;
         }

         // Check Veto from Read 1 Perspective (if not already vetoed)
         if (!veto) {
             OrientedReadId::Int u2 = OrientedReadId(ad.readIds[1], 0).getValue();
             if (u2 < numVertices && boost::out_degree(u2, graph) > 0) {
                 // Read 1 is phased. It MUST have an edge to Read 0 with correct orientation.
                 // If ad.isSameStrand: R1(0) -> R0(0)
                 // If !ad.isSameStrand: R1(0) -> R0(1)
                 OrientedReadId::Int v_target = OrientedReadId(ad.readIds[0], ad.isSameStrand ? 0 : 1).getValue();
                 
                 bool found = false;
                 auto outEdges = boost::out_edges(u2, graph);
                 for(auto it = outEdges.first; it != outEdges.second; ++it) {
                     if (boost::target(*it, graph) == v_target) {
                         found = true;
                         break;
                     }
                 }
                 if (!found) veto = true;
             }
         }

         if (!veto) {
             alignmentData[alignmentId].info.isInReadGraph = 1;
             keptAlignmentCount++;
         } else {
             keepAlignment[alignmentId] = false;
             alignmentData[alignmentId].info.isInReadGraph = 0;
         }
    }

    cout << timestamp << "Kept " << keptAlignmentCount << " / " << alignmentCount << " alignments after haplotype filtering." << endl;

    // Create the read graph using FILTERED alignments.
    createReadGraphUsingSelectedAlignments(keepAlignment);
}


// Global Cluster Validity Check Function to check if a cluster is valid.
// A cluster is valid if:
// 1. It has at least 2 alleles with coverage >= minAlleleCoverage
// 2. It is not in a homopolymer region (checked on multiple reads)
// 3. It does not exhibit strand bias (following hifiasm's strand balance criteria)
//
// Strand Bias Filter (hifiasm-style):
// - Each allele must have reads from BOTH strands (at least 2 from each)
// - The minority strand must be >= 35% of total coverage for that allele
//   (relaxed to 24% if minority strand has >= 5 reads)
// - Total coverage for the site must be >= 10 reads
void Assembler::computeClusterValidityThreadFunction(uint64_t threadId) {
    const uint64_t currentMinAlleleCoverage = this->minAlleleCoverage;
    
    // === Strand Bias Thresholds (hifiasm-style) ===
    const uint32_t minReadsPerStrand = 2;           // Minimum reads from each strand per allele
    const uint32_t minStrand0Coverage = 4;          // Minimum reads from strand 0 per allele
    const double strandImbalanceThreshold = 0.35;   // Minority strand must be >= 35% of total
    const double relaxedImbalanceThreshold = 0.24;  // Relaxed to 24% when min strand >= 5 reads
    const uint32_t relaxedMinStrandCoverage = 5;    // Threshold for using relaxed imbalance
    const uint32_t minTotalCoverage = 10;           // Minimum total reads per allele
    
    // === Homopolymer Thresholds ===
    const double homopolymerFractionThreshold = 0.5;  // If >50% of reads show homopolymer context, filter
    const uint64_t minReadsForHomopolymerCheck = 3;   // Need at least 3 reads to reliably check
    const uint64_t maxHomopolymerChecks = 10;         // Max reads to check for homopolymer (for efficiency)
    
    // Loop over a batch of cluster IDs assigned to this thread
    uint64_t clusterIdBegin, clusterIdEnd;
    while (getNextBatch(clusterIdBegin, clusterIdEnd)) {
        for (uint64_t clusterId = clusterIdBegin; clusterId < clusterIdEnd; clusterId++) {
            const auto& members = variantClusteringMembersByRepIdx[clusterId];
            
            if (members.empty()) continue;
            
            // === SINGLE PASS: Collect all data in one iteration ===
            // Allele counts
            std::array<uint32_t, 5> alleleCounts = {0};
            // Per-allele strand counts (for strand bias detection)
            std::array<uint32_t, 5> strand0Counts = {0};
            std::array<uint32_t, 5> strand1Counts = {0};
            
            // Fixed-size array for homopolymer check positions (avoid heap allocation)
            std::array<std::pair<OrientedReadId, uint32_t>, 10> homopolymerCheckPositions;
            uint64_t homopolymerCheckCount = 0;

            // Single pass through all members
            for (uint64_t memberIdx : members) {
                const auto& pp = variantClusteringPositionPairs[memberIdx];
                const OrientedReadId orientedReadId = pp.first;
                const uint32_t position = pp.second;
                
                const uint8_t allele = variantClusteringPositionPairAlleles[memberIdx];
                if (allele < 5) {
                    alleleCounts[allele]++;
                    // Track strand counts per allele
                    if (orientedReadId.getStrand() == 0) {
                        strand0Counts[allele]++;
                    } else {
                        strand1Counts[allele]++;
                    }
                }
                
                // Collect first few positions for homopolymer check
                if (homopolymerCheckCount < maxHomopolymerChecks) {
                    homopolymerCheckPositions[homopolymerCheckCount++] = {orientedReadId, position};
                }
            }

            // === STEP 2: Check allele coverage threshold ===
            std::array<uint8_t, 5> significantAlleles;
            uint8_t significantCount = 0;
            for (uint8_t a = 0; a < 5; a++) {
                if (alleleCounts[a] >= currentMinAlleleCoverage) {
                    significantAlleles[significantCount++] = a;
                }
            }
            
            // Must have at least 2 valid alleles to be a het site
            if (significantCount < 2) {
                continue;
            }

            // === STEP 3: Check for strand bias (hifiasm-style) ===
            // TODO: Re-enable strand bias filtering after testing
            /*
            // Find the two alleles with highest coverage
            uint8_t allele1 = significantAlleles[0];
            uint8_t allele2 = significantAlleles[1];
            
            if (significantCount > 2) {
                // Find top 2 by coverage using simple comparison
                for (uint8_t i = 2; i < significantCount; i++) {
                    uint8_t candidate = significantAlleles[i];
                    if (alleleCounts[candidate] > alleleCounts[allele1]) {
                        allele2 = allele1;
                        allele1 = candidate;
                    } else if (alleleCounts[candidate] > alleleCounts[allele2]) {
                        allele2 = candidate;
                    }
                }
                // Ensure allele1 has higher coverage
                if (alleleCounts[allele2] > alleleCounts[allele1]) {
                    std::swap(allele1, allele2);
                }
            }
            
            // Lambda to check strand bias for a single allele (hifiasm filter_one_snp logic)
            auto passesStrandFilter = [&](uint8_t allele) -> bool {
                uint32_t occ_0 = strand0Counts[allele];
                uint32_t occ_1 = strand1Counts[allele];
                uint32_t total = occ_0 + occ_1;
                
                // Must have reads from both strands
                if (occ_0 < minReadsPerStrand || occ_1 < minReadsPerStrand) {
                    return false;
                }
                
                // Minimum total coverage
                if (total < minTotalCoverage) {
                    return false;
                }
                
                // Calculate strand balance
                uint32_t minStrand = std::min(occ_0, occ_1);
                double available = static_cast<double>(minStrand) / static_cast<double>(total);
                
                // Use relaxed threshold if minority strand has enough coverage
                if (minStrand >= relaxedMinStrandCoverage) {
                    // Relaxed threshold (24%)
                    if (available < relaxedImbalanceThreshold) {
                        return false;
                    }
                } else {
                    // Standard threshold (35%) + minimum strand 0 coverage
                    if (available < strandImbalanceThreshold || occ_0 < minStrand0Coverage) {
                        return false;
                    }
                }
                
                return true;
            };
            
            // Check strand bias for both top alleles
            bool hasStrandBias = !passesStrandFilter(allele1) || !passesStrandFilter(allele2);
            
            if (hasStrandBias) {
                continue;  // Filter out strand-biased sites
            }
            */

            // === STEP 4: Check for homopolymer context ===
            // TODO: Re-enable homopolymer filtering after testing strand bias filter
            /*
            bool isHomopolymerSite = false;
            
            if (homopolymerCheckCount >= minReadsForHomopolymerCheck) {
                uint64_t homopolymerHitCount = 0;
                uint64_t checkedCount = 0;
                
                for (uint64_t i = 0; i < homopolymerCheckCount; i++) {
                    const OrientedReadId orientedReadId = homopolymerCheckPositions[i].first;
                    const uint32_t position = homopolymerCheckPositions[i].second;
                    
                    // Get the read length
                    const uint64_t readLength = getReads().getReadRawSequenceLength(orientedReadId.getReadId());
                    
                    // Make sure position is valid
                    if (position >= readLength) continue;
                    
                    // Get read sequence for homopolymer check
                    LongBaseSequenceView readSequence = getReads().getRead(orientedReadId.getReadId());
                    
                    // For strand 1, we need to check the reverse complement position
                    uint64_t checkPosition = position;
                    if (orientedReadId.getStrand() == 1) {
                        checkPosition = readLength - 1 - position;
                    }
                    
                    if (checkPosition < readLength) {
                        if (isSiteInHomopolymerRegion(checkPosition, readSequence)) {
                            homopolymerHitCount++;
                        }
                        checkedCount++;
                    }
                }
                
                // If majority of reads show homopolymer context, mark as homopolymer site
                if (checkedCount >= minReadsForHomopolymerCheck) {
                    double homopolymerFraction = static_cast<double>(homopolymerHitCount) / static_cast<double>(checkedCount);
                    if (homopolymerFraction > homopolymerFractionThreshold) {
                        isHomopolymerSite = true;
                    }
                }
            }
            
            if (isHomopolymerSite) {
                continue;  // Filter out homopolymer sites
            }
            */

            // === STEP 5: Cluster passes all filters, mark as valid ===
            variantClusteringValidClusters[clusterId] = 1;
        }
    }
}


// Helper function for Strand Bias Filtering
// Returns TRUE if the allele is BIASED towards Strand 0 (Forward).
//
// Arguments:
//   refAlleleCoverage:          Total number of reads with the reference allele.
//   refAlleleCoverageOnStrand0: Number of reads with the reference allele on Strand 0 (Forward).
//   mm:                         Max reads allowed on the "minor" strand (default 2).
//   rr:                         Minimum ratio of "major" strand reads (default 0.05 => 95%).
static bool isStrandBiased(uint64_t refAlleleCoverage, uint64_t refAlleleCoverageOnStrand0, uint64_t mm = 2, double rr = 0.05) {
    
    // Site is biased if all reads with the reference allele are on one strand
    if (refAlleleCoverage == refAlleleCoverageOnStrand0) return true;

    // Condition 1: Few minor strand reads (major + mm >= total)
    bool cond1 = (refAlleleCoverageOnStrand0 + mm >= refAlleleCoverage);

    // Condition 2: High major strand percentage
    bool cond2 = (double(refAlleleCoverage) * rr + double(refAlleleCoverageOnStrand0) >= double(refAlleleCoverage));

    // We only check for bias towards Strand 0 because the reference read will be on Strand 0.
    // If the allele is heavily biased towards Strand 0, it is likely a systematic error 
    // rather than a true variant supported by both strands.
    if (cond1 && cond2) {
        return true;
    }

    return false;
}

void Assembler::createReadGraph5ThreadFunction(uint64_t threadId)
{
    // Access data
    const auto& positionPairs = variantClusteringPositionPairs;
    auto& disjointSets = *variantClusteringDisjointSets;
    const auto& membersByRepIdx = variantClusteringMembersByRepIdx;

    uint64_t readIdBegin, readIdEnd;
    while (getNextBatch(readIdBegin, readIdEnd)) {
        for (ReadId readId = readIdBegin; readId < readIdEnd; readId++) {
            
            // Only process strand 0
            // The process is symmetric for the two strands.
            OrientedReadId currentOrientedReadId0(readId, 0);
            
            // Find clusters involving this read on strand 0.
            // The position pairs are sorted by read id and strand.
            auto it = std::lower_bound(positionPairs.begin(), positionPairs.end(), make_pair(currentOrientedReadId0, 0u));
            
            struct ClusterOnRead {
                uint64_t clusterId;
                uint32_t position;
                uint64_t indexInPositionPairs;
            };
            vector<ClusterOnRead> clusters;

            // Find clusters that involve this read on strand 0.
            while(it != positionPairs.end() && it->first == currentOrientedReadId0) {
                uint64_t index = it - positionPairs.begin();
                uint64_t clusterId = disjointSets.find(index);
                clusters.push_back({clusterId, it->second, index});
                ++it;
            }

            // If no clusters that involve this read on strand 0, skip it
            if (clusters.empty()) {
                continue;
            }

            vector<ClusterOnRead> filteredClusters;
            for(const auto& cluster : clusters) {
                // O(1) Lookup
                if (variantClusteringValidClusters[cluster.clusterId]) {
                    filteredClusters.push_back(cluster);
                }
            }
            
            // If no valid clusters after filtering, skip this read
            if (filteredClusters.empty()) {
                continue;
            }
            
            // Use filtered clusters for the rest of the analysis
            clusters = std::move(filteredClusters);




            // --------------------------------------------------------------
            // Compatibility Check with Multi-Allelic Splitting
            // --------------------------------------------------------------

            // Check if we have any clusters to process
            if (clusters.empty()) {
                continue;
            }

            // --- Step 1: Generate Virtual Clusters (Nodes for DP) ---
            // A physical cluster with alleles {Target, A, B} becomes 
            // two virtual clusters: (Target vs A) and (Target vs B).
            
            struct VirtualCluster {
                uint64_t clusterId;
                uint8_t targetAllele;
                uint8_t altAllele;
                // ID in the original 'clusters' vector to retrieve members
                uint32_t originalIndex; 
            };
            vector<VirtualCluster> dpNodes;
            dpNodes.reserve(clusters.size() * 2); // Heuristic reserve

            const uint64_t currentMinAlleleCoverage = this->minAlleleCoverage; 

            for(uint32_t i=0; i<clusters.size(); i++) {
                const auto& cluster = clusters[i];
                
                // Identify Target Allele (The allele supported by the current read we are phasing)
                uint8_t targetAllele = 255;
                if (cluster.indexInPositionPairs < variantClusteringPositionPairAlleles.size()) {
                    targetAllele = variantClusteringPositionPairAlleles[cluster.indexInPositionPairs];
                }
                DINARA_ASSERT(targetAllele < 5); // Invalid target allele

                // Count alleles in this cluster
                // ALSO count Strand-0 (Forward) occurrences for bias checking.
                const auto& members = membersByRepIdx[cluster.clusterId];
                std::array<uint32_t, 5> alleleCounts = {0};       // Total
                std::array<uint32_t, 5> alleleStrand0Counts = {0}; // Forward (Strand 0)
                
                for(uint64_t memberIdx : members) {
                    if(memberIdx < variantClusteringPositionPairAlleles.size()) {
                        uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                        if(a < 5) {
                            alleleCounts[a]++;
                            if (this->variantClusteringPositionPairs[memberIdx].first.getStrand() == 0) {
                                alleleStrand0Counts[a]++;
                            }
                        }
                    }
                }

                // Check Strand Bias for the Reference Allele
                // If the Reference allele is strand biased, we filter this physical cluster entirely.
                if (isStrandBiased(alleleCounts[targetAllele], alleleStrand0Counts[targetAllele])) {
                    continue;
                }

                // Check if the Reference Allele has sufficient coverage
                if (alleleCounts[targetAllele] < currentMinAlleleCoverage) {
                    continue;
                }

                // Create a Virtual Cluster for each significant ALT allele
                // It is guaranteed that there is at least one significant ALT allele.
                // The point is to split multiallelic clusters into binary het clusters.
                for(uint8_t a=0; a<5; a++) {
                    if (a == targetAllele) continue;
                    if (alleleCounts[a] >= currentMinAlleleCoverage) {
                        dpNodes.push_back({cluster.clusterId, targetAllele, a, i});
                    }
                }
            }

            const size_t N = dpNodes.size();
            if (N == 0) continue;

            // --- Step 2: Build Read Phase Vectors for each DP Node ---
            // nodeReads[i] contains list of {ReadId, isPhase1} for dpNodes[i]
            // Phase 0 = Matches Target. Phase 1 = Matches Alt. 

            struct ReadPhase {
                OrientedReadId::Int orientedReadIdValue; // Using the raw integer value
                bool isPhase1; // false=Target, true=Alt
            };
            vector<vector<ReadPhase>> nodeReads(N);

            for(size_t i=0; i<N; i++) {
                const auto& node = dpNodes[i];
                const auto& members = membersByRepIdx[node.clusterId];
                nodeReads[i].reserve(members.size());

                for(uint64_t memberIdx : members) {
                    const auto& pp = positionPairs[memberIdx];
                    
                    OrientedReadId orientedReadId = pp.first;
                    
                    // Skip target oriented read itself (given that the clusters involve this read on strand 0).
                    if(orientedReadId == currentOrientedReadId0) continue;

                    if(memberIdx < variantClusteringPositionPairAlleles.size()) {
                        uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                        
                        if (a == node.targetAllele) {
                            nodeReads[i].push_back({orientedReadId.getValue(), false}); 
                        } else if (a == node.altAllele) {
                            nodeReads[i].push_back({orientedReadId.getValue(), true});  
                        }
                    }
                }
                
                // Sort by OrientedReadId value for fast intersection
                std::sort(nodeReads[i].begin(), nodeReads[i].end(), 
                    [](const ReadPhase& a, const ReadPhase& b) { 
                        return a.orientedReadIdValue < b.orientedReadIdValue; 
                    });
            }

            // --- Step 3: DP for Longest Compatible Group (LCG) ---
            vector<int> LCG(N, 1);
            vector<int> parent(N, -1);

            for(int i = 0; i < N; i++) {
                for(int j = 0; j < i; j++) {
                    // Don't link two virtual nodes coming from the SAME physical cluster.
                    // This is to avoid linking multiallelic sites coming from the same physical cluster.
                    // We force the DP chain to pick at most one virtual node per physical genomic site, 
                    // ensuring we don't double-count the same site or create invalid paths.
                    // This enforces the constraint: "One genomic location = One node in the path."
                    if (dpNodes[i].clusterId == dpNodes[j].clusterId) continue;

                    // Check compatibility(i, j)
                    // Compatible if NO read conflicts.
                    // Conflict: Read r present in both, but has Phase 0 in one and Phase 1 in other.
                    
                    bool compatible = true;
                    const auto& orientedReadsI = nodeReads[i];
                    const auto& orientedReadsJ = nodeReads[j];
                    
                    auto itI = orientedReadsI.begin();
                    auto itJ = orientedReadsJ.begin();
                    
                    int supportPhase0 = 0; // Reads linking Target -> Target
                    int supportPhase1 = 0; // Reads linking Alt -> Alt

                    // Linear intersection scan to compare sorted sets.
                    while(itI != orientedReadsI.end() && itJ != orientedReadsJ.end()) {
                        if(itI->orientedReadIdValue < itJ->orientedReadIdValue) {
                            ++itI;
                        } else if(itJ->orientedReadIdValue < itI->orientedReadIdValue) {
                            ++itJ;
                        } else {
                            // Same OrientedRead covers both.
                            if(itI->isPhase1 != itJ->isPhase1) {
                                // CONFLICT found! (Phase 0 -> 1 or 1 -> 0).
                                compatible = false; 
                                break; 
                            } else {
                                // Consistent!
                                if (itI->isPhase1 == false) supportPhase0++;
                                else supportPhase1++;
                            }
                            ++itI;
                            ++itJ;
                        }
                    }

                    // COMPATIBILITY CRITERIA:
                    // 1. Must be compatible (no conflicts)
                    // 2. Must have ACTIVE support (shared reads)
                    
                    if (compatible) {
                        // "Strong Linkage" (Support for BOTH alleles)
                        // This ensures we are tracking two real haplotypes, not just linking noise.
                        if (supportPhase0 > 0 && supportPhase1 > 0) {
                            if(LCG[j] + 1 > LCG[i]) {
                                LCG[i] = LCG[j] + 1; // Longest compatible group length
                                parent[i] = j; // Parent node in the DP chain
                            }
                        }
                    }
                }
            }

            // --- Step 4: Traceback and Chain Identification (Multi-Chain) ---
            // This is the core of the DP algorithm.
            // It identifies the longest compatible group of virtual nodes.
            // It also identifies isolated sites that are supported by a sufficiently high number of reads.
            // It then marks the clusters as valid and collects the reads that support the Target and Alt alleles.
 
            vector<bool> isAssigned(N, false);
            vector<pair<int, int>> sortedLCGIndices; // (Length, Index)
            vector<pair<int, int>> isolatedSitesIndices;

            for(int i = 0; i < N; i++) {
                if (LCG[i] > 1) {
                    sortedLCGIndices.push_back({LCG[i], i});
                } else {
                    isolatedSitesIndices.push_back({LCG[i], i});
                }
            }

            // Sort chains by length descending
            std::sort(sortedLCGIndices.rbegin(), sortedLCGIndices.rend());
            vector<vector<int>> validChains;

            const uint64_t minMultiChainCoverage = 6; // Stricter threshold for multi-chain sites
            
            // 1. Process Multi-Node Chains
            for (const auto& p : sortedLCGIndices) {
                int endNode = p.second;
                if (!isAssigned[endNode]) {
                    vector<int> chain;
                    int curr = endNode;
                    while (curr != -1 && !isAssigned[curr]) {
                        // Calculate target allele support
                        uint64_t targetAlleleSupport = 0;
                        for (const auto& rp : nodeReads[curr]) {
                            if (!rp.isPhase1) targetAlleleSupport++;
                        }

                        // Filter by Target Allele Support
                        if (targetAlleleSupport >= minMultiChainCoverage) {
                            chain.push_back(curr);
                        }
                        isAssigned[curr] = true;
                        curr = parent[curr];
                    }
                    
                    // No need to reverse if we just iterate, but for correctness:
                    // Check if it is empty chain
                    if (chain.size() > 0) {
                        std::reverse(chain.begin(), chain.end());
                        validChains.push_back(chain);
                    }
                }
            }

            // -------------------------------------------------------------------------
            // Filter Helper for Isolated Sites in Homopolymer Regions
            // -------------------------------------------------------------------------
            auto isHomopolymerDominated = [&](const VirtualCluster& node) -> bool {
                const auto& members = membersByRepIdx[node.clusterId];
                
                int64_t n0 = 0;
                int64_t n1 = 0;

                // 1. Initial counts
                for(uint64_t memberIdx : members) {
                     if(memberIdx < variantClusteringPositionPairAlleles.size()) {
                        uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                        if (a == node.targetAllele) n0++;
                        else if (a == node.altAllele) n1++;
                     }
                }
                
                int64_t occ_0 = n0;
                int64_t occ_1 = n1;

                // 2. Remove votes from reads in homopolymer regions
                // We access the reads object once outside loop if possible, but inside lambda we need 'this'.
                const auto& reads = this->getReads();

                for(uint64_t memberIdx : members) {
                    if(memberIdx < variantClusteringPositionPairAlleles.size()) {
                        const auto& pp = positionPairs[memberIdx];
                        ReadId rId = pp.first.getReadId();
                                                
                        uint32_t pos = pp.second; 

                        // Retrieve the raw sequence (LongBaseSequenceView)
                        dinara::LongBaseSequenceView readSequence = reads.getRead(rId);

                        // Check if this read is in a homopolymer region at this site
                        // Note: readSequence is always Strand 0 (Forward).
                        // If the read is on Strand 1 (Reverse), the position 'pos' is relative to the RC read.
                        // We need to convert 'pos' to the corresponding coordinate on Strand 0.
                        uint64_t checkPos = pos;
                        if (pp.first.getStrand() == 1) {
                            uint64_t readLen = readSequence.baseCount;
                            DINARA_ASSERT(pos < readLen);
                            checkPos = readLen - 1 - pos;
                        }
                        
                        // "Keep the default values for the other parameters"
                        if (isSiteInHomopolymerRegion(checkPos, readSequence)) {
                            uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                            if (a == node.targetAllele) {
                                occ_0--;
                            } else if (a == node.altAllele) {
                                occ_1--;
                            }
                        }
                    }
                }

                // 3. Check thresholds
                if (occ_0 < 2 || occ_1 < 2 || !(occ_0 >= currentMinAlleleCoverage && occ_1 >= currentMinAlleleCoverage)) return true;

                return false;
            };

            // 2. Process Isolated Sites
            // An isolated site is considered informative only if the reference read allele 
            // is supported by a sufficiently high number of reads AND is not HP-dominated.
            const uint64_t minIsolatedCoverage = 14; // Stricter threshold for isolated sites
            
            for (const auto& p : isolatedSitesIndices) {
                int nodeIdx = p.second;
                if (!isAssigned[nodeIdx]) {
                    // Check support for this specific Virtual Cluster
                    
                    // Count reads supporting the reference read allele (Target Allele / Phase 0)
                    uint64_t targetAlleleSupport = 0;
                    for (const auto& rp : nodeReads[nodeIdx]) {
                        if (!rp.isPhase1) { // isPhase1=false means Target Allele
                            targetAlleleSupport++;
                        }
                    }

                    if (targetAlleleSupport >= minIsolatedCoverage) {
                        // Check Homopolymer Filter
                        if (!isHomopolymerDominated(dpNodes[nodeIdx])) {
                            validChains.push_back({nodeIdx});
                            isAssigned[nodeIdx] = true;
                        }
                    }
                }
            }

            // --- Step 5: Haplotype Assignment (Naive HiFi) ---

            // Data structures for the algorithm
            struct SiteStat {
                uint32_t dpNodeIdx;
                uint32_t position;
                uint32_t occ_0; // Target allele count
                uint32_t occ_1; // Alt allele count
                uint32_t score; // 0=Unconfirmed, 1=Confirmed Het
                bool is_filtered;
            };
            
            struct ReadSiteInfo {
                uint32_t siteIndex; // Index in 'sites' vector
                bool isPhase1;      // false=Target, true=Alt
            };

            struct ReadStat {
                OrientedReadId::Int orientedReadIdValue;
                std::vector<ReadSiteInfo> coveredSites;
                uint32_t informative_score; // 'o' in user algo
                uint8_t is_trans;           // 0=Unknown, 1=CIS, 2=TRANS
            };

            // 1. Data Preparation
            // Collect unique sites from all valid chains
            std::vector<SiteStat> sites;
            std::vector<ReadStat> readStats;
            std::unordered_map<OrientedReadId::Int, size_t> readMap; // ReadValue -> Index in readStats

            // We need to flatten the chains into a list of unique sites
            std::set<int> usedNodeIndices;
            for (const auto& chain : validChains) {
                for (int nodeIdx : chain) {
                    usedNodeIndices.insert(nodeIdx);
                }
            }
            
            // Mark Compatible Clusters
            for (int nodeIdx : usedNodeIndices) {
                 uint64_t originalClusterId = dpNodes[nodeIdx].clusterId;

                // Mark compatible atomically
                __sync_bool_compare_and_swap(&variantClusteringValidClustersCompatible[originalClusterId], 0, 1);

                // Mark RC cluster
                {
                    const auto& members = membersByRepIdx[originalClusterId];
                    if (!members.empty()) {
                        uint64_t memberIdx = members[0];
                        const auto& pp = positionPairs[memberIdx];
                        ReadId readId = pp.first.getReadId();
                        Strand strand = pp.first.getStrand();
                        uint32_t position = pp.second;
                        
                        Strand rcStrand = 1 - strand;
                        OrientedReadId rcOrientedReadId(readId, rcStrand);
                        uint64_t readLength = getReads().getReadRawSequenceLength(readId);
                        uint32_t rcPosition = uint32_t(readLength - 1 - position);
                        
                        auto searchKey = std::make_pair(rcOrientedReadId, rcPosition);
                        auto itRc = std::lower_bound(positionPairs.begin(), positionPairs.end(), searchKey);
                        if (itRc != positionPairs.end() && *itRc == searchKey) {
                            uint64_t rcMemberIdx = itRc - positionPairs.begin();
                            uint64_t rcClusterId = disjointSets.find(rcMemberIdx);
                            if (rcClusterId != originalClusterId) {
                                __sync_bool_compare_and_swap(&variantClusteringValidClustersCompatible[rcClusterId], 0, 1);
                            }
                        }
                    }
                }
            }

            // Build Sites and Reads
            for (int nodeIdx : usedNodeIndices) {
                // Get Position
                uint32_t pos = 0;
                uint64_t cId = dpNodes[nodeIdx].clusterId;
                if (!membersByRepIdx[cId].empty()) {
                    pos = positionPairs[membersByRepIdx[cId][0]].second;
                }

                // Initial Counts (Target=0, Alt=0 - will fill from reads)
                SiteStat s = {(uint32_t)nodeIdx, pos, 0, 0, 0, false};
                
                uint32_t currentSiteIdx = (uint32_t)sites.size();
                sites.push_back(s);

                // Process Reads in this Node
                for (const auto& rp : nodeReads[nodeIdx]) {
                    // Find or Create ReadStat
                    if (readMap.find(rp.orientedReadIdValue) == readMap.end()) {
                        size_t idx = readStats.size();
                        readStats.push_back({rp.orientedReadIdValue, {}, 0, 0});
                        readMap[rp.orientedReadIdValue] = idx;
                    }
                    size_t rIdx = readMap[rp.orientedReadIdValue];
                    
                    // Add site info to read
                    readStats[rIdx].coveredSites.push_back({currentSiteIdx, rp.isPhase1});

                    // Update Site Counts (temporarily, will refine later)
                    if (!rp.isPhase1) sites.back().occ_0++;
                    else sites.back().occ_1++;
                }
            }

            // 2.1 Adjacent Site Filter
            // Sort sites by position
            // Note: 'sites' vector order corresponds to indices in ReadSiteInfo. 
            // If we sort 'sites', we break those indices.
            // Instead of sorting 'sites' and breaking indices, let's create a sorted view or just iterate carefully.
            // Actually, sorting sites by position is needed for the 1bp check. 
            // We can check 1bp diff by comparing all pairs O(N^2) or sorting. N is small (chain nodes).
            // Let's sort an index vector.
            std::vector<uint32_t> sitesByPos(sites.size());
            std::iota(sitesByPos.begin(), sitesByPos.end(), 0);
            std::sort(sitesByPos.begin(), sitesByPos.end(), [&](uint32_t a, uint32_t b) {
                return sites[a].position < sites[b].position;
            });

            for (size_t i = 0; i < sites.size(); ++i) {
                // Check neighbors in sorted list
                uint32_t siteIdx = sitesByPos[i];
                uint32_t myPos = sites[siteIdx].position;
                
                // Check Prev
                if (i > 0) {
                    uint32_t prevIdx = sitesByPos[i-1];
                    if (sites[prevIdx].position + 1 == myPos) {
                        sites[siteIdx].is_filtered = true;
                    }
                }
                // Check Next
                if (i + 1 < sites.size()) {
                    uint32_t nextIdx = sitesByPos[i+1];
                    if (sites[nextIdx].position == myPos + 1) {
                         sites[siteIdx].is_filtered = true;
                    }
                }
            }

            // 2.2 Count Informative Sites Per Read
            for (auto& r : readStats) {
                r.informative_score = 0;
                for (const auto& info : r.coveredSites) {
                     SiteStat& s = sites[info.siteIndex];
                     
                     if (s.is_filtered) continue;
                     
                     // Must be mismatch (Phase 1 / Alt) to be "informative" for TRANS detection?
                     // User code: "if (hh_tp(hap->list[i]) != 1) continue; // Must be mismatch"
                     // Wait, hh_tp != 1 means Mismatch? Or Match?
                     // In user code: 
                     //   Line 8973: if (overlap_list->list[ii].is_match == 1) overlap_list->list[ii].is_match = 2; (Mark TRANS)
                     //   Inside loop: if (hh_tp == 1) s->score=1 (MISMATCH pos)
                     // So hh_tp=1 is Mismatch/Alt.
                     // User 2.2 says: "if (hh_tp != 1) continue; // Must be mismatch"
                     // So yes, we only count Mismatches (Alt alleles) as "informative" for *identifying* trans overlaps?
                     // But later: "Only overlaps with o > 0 are added".
                     // Let's follow user: info.isPhase1 == true (Alt)
                     if (!info.isPhase1) continue; 

                     // Coverage check
                     if (s.occ_0 < 2 || s.occ_1 < 2) continue;

                     // Strand bias check (Skip for now, or assume handled upstream)
                     
                     // "if (s->occ_0 >= 3 && s->occ_1 >= 3) o++"
                     if (s.occ_0 >= 3 && s.occ_1 >= 3) {
                         r.informative_score++;
                     }
                }
            }

            // 2.3 Sort Reads (Overlaps)
            // Sort by informative_score descending
            std::vector<uint32_t> readIndices(readStats.size());
            std::iota(readIndices.begin(), readIndices.end(), 0);
            std::sort(readIndices.begin(), readIndices.end(), [&](uint32_t a, uint32_t b) {
                return readStats[a].informative_score > readStats[b].informative_score;
            });

            // 2.4 First Loop: Assign Trans and Update Coverage
            for (uint32_t rIdx : readIndices) {
                ReadStat& r = readStats[rIdx];
                
                // If informative score > 0 (from 2.2 logic: overlap has mismatch at strong site)
                // "Only overlaps with o > 0 are added to processing queue" implied filter?
                // Or just proceed. User says "For each overlap...".
                // User code: "if (overlap_list->list[ii].is_match == 1) is_match = 2"
                // This implies ALL reads in this sorted list are candidates for TRANS if they have mismatches?
                // Actually, the sorting key `o` is based on mismatches. If `o > 0`, it has mismatches at strong sites.
                // So likely these are the ones we assume are TRANS first?
                // Wait, "2.4 ... For each overlap (most informative first): ... mark as TRANS"
                // It seems it aggressively marks them as TRANS if they have mismatches?
                // Let's implement the loop actions.
                
                // Check if it's candidate for TRANS
                // The user logic seems to assume we are processing a list of potential *TRANS* overlaps?
                // "Only overlaps with o > 0 (at least one mismatch at informative site) are added to the processing queue."
                // So if informative_score == 0, we skip this logic (it remains default/CIS or UNKNOWN).
                if (r.informative_score == 0) continue;

                // Mark as TRANS
                r.is_trans = 2; // TRANS

                for (const auto& info : r.coveredSites) {
                    SiteStat& s = sites[info.siteIndex];
                    if (s.is_filtered) continue;

                    if (info.isPhase1) { // Mismatch / Alt
                        // "MISMATCH position -> mark site as confirmed het"
                        s.score = 1; 
                    } else { // Match / Target
                        // "MATCH position -> remove this overlap's ref vote"
                        if (s.occ_0 > 0) s.occ_0--;
                    }
                }
            }

            // 2.7 Second Loop: Catch-Up Assignment
            // "For overlaps not yet assigned" -> i.e. informative_score == 0 or not visited?
            // "or overlaps not yet assigned, check against confirmed het sites"
            for (uint32_t rIdx : readIndices) {
                ReadStat& r = readStats[rIdx];
                
                // If already TRANS, skip
                if (r.is_trans == 2) continue;

                int o = 0;
                for (const auto& info : r.coveredSites) {
                    SiteStat& s = sites[info.siteIndex];
                    // "if (hh_tp != 1) continue" -> Mismatch only
                    if (!info.isPhase1) continue;
                    
                    // "if (s->occ_0 < 2 || s->occ_1 < 2) continue" -> Updated counts!
                    if (s.occ_0 < 2 || s.occ_1 < 2) continue;
                    
                    // "if (s->score == 1) o++" -> Confirmed het
                    if (s.score == 1) o++;
                }

                // Mark as TRANS if mismatch at confirmed het site
                if (o > 0) {
                    r.is_trans = 2;
                }
            }

            // 2.8 Output Generation
            // Default: if is_trans != 2, it is CIS (is_trans=1 or 0)
            
            vector<OrientedReadId::Int> finalHapOrientedReadIds;
            // Add Target Read itself (always in phase)
            finalHapOrientedReadIds.push_back(currentOrientedReadId0.getValue());

            for (const auto& r : readStats) {
                if (r.is_trans != 2) {
                    // In-Phase (CIS)
                     finalHapOrientedReadIds.push_back(r.orientedReadIdValue);
                } 
                // else Out-Of-Phase (TRANS) - dropped from final haplotype helper
            }

            // Sort & Unique
            std::sort(finalHapOrientedReadIds.begin(), finalHapOrientedReadIds.end());
            finalHapOrientedReadIds.erase(std::unique(finalHapOrientedReadIds.begin(), finalHapOrientedReadIds.end()), finalHapOrientedReadIds.end());

            // Print finalHapOrientedReadIds with one orientedReadId per line if readId is 0 (debug)
            if(currentOrientedReadId0.getReadId() == 87 and currentOrientedReadId0.getStrand() == 0) {
                cout << "Final Haplotype Oriented Reads: " << endl;
                for(OrientedReadId::Int orientedReadIdValue : finalHapOrientedReadIds) {
                    cout << OrientedReadId::fromValue(orientedReadIdValue).getString() << endl;
                }
            }

            // CAST VOTES DIRECTLY TO GRAPH
            // Calculate total weight (total number of nodes in all valid chains)
            // Note: The concept of "nodes" might be less relevant now that we filter sites,
            // but we can still use validChains size or just count valid sites?
            // User code didn't specify weight changes. Let's keep original totalWeight logic.
            uint32_t totalWeight = 0;
             for (const auto& chain : validChains) {
                totalWeight += (uint32_t)chain.size();
            }

            {
                std::lock_guard<std::mutex> lock(haplotypeGraphMutex);
                for(OrientedReadId::Int otherOrientedReadIdValue : finalHapOrientedReadIds) {
                    if(otherOrientedReadIdValue != currentOrientedReadId0.getValue()) {
                        // Try to add edge with totalWeight
                        auto result = boost::add_edge(currentOrientedReadId0.getValue(), otherOrientedReadIdValue, totalWeight, *globalHaplotypeGraph);
                    }
                }
            }










            // KEEP THAT FOR REFERENCE - DO NOT DELETE
            // // --- Step 5: Mark Valid Clusters & Identify Haplotype Reads ---
            
            // // Collect all In-Phase reads (supporting Target Allele) 
            // // and Out-Of-Phase reads (supporting Alt Allele) from VALID chains.
            // vector<uint64_t> inPhaseReads;
            // vector<uint64_t> outOfPhaseReads;

            // for (const auto& chain : validChains) {
            //     for (int nodeIdx : chain) {

            //         // Mark cluster as valid in global array
            //         uint64_t originalClusterId = dpNodes[nodeIdx].clusterId;

            //         // Benign race: multiple threads may write 1 simultaneously.
            //         // Mark compatible atomically to avoid data races
            //         __sync_bool_compare_and_swap(&variantClusteringValidClustersCompatible[originalClusterId], 0, 1);

            //         // Also mark the RC cluster as compatible to maintain symmetry
            //         // Since we only process strand 0 reads, we need to find and mark the RC cluster
            //         {
            //             // Get a member from this cluster to find its RC
            //             const auto& members = membersByRepIdx[originalClusterId];
            //             if (!members.empty()) {
            //                 uint64_t memberIdx = members[0];
            //                 const auto& pp = positionPairs[memberIdx];
            //                 ReadId readId = pp.first.getReadId();
            //                 Strand strand = pp.first.getStrand();
            //                 uint32_t position = pp.second;
                            
            //                 // Find the RC: same read, opposite strand, complementary position
            //                 Strand rcStrand = 1 - strand;
            //                 OrientedReadId rcOrientedReadId(readId, rcStrand);
            //                 uint64_t readLength = getReads().getReadRawSequenceLength(readId);
            //                 uint32_t rcPosition = uint32_t(readLength - 1 - position);
                            
            //                 // Binary search for the RC position pair
            //                 auto searchKey = std::make_pair(rcOrientedReadId, rcPosition);
            //                 auto itRc = std::lower_bound(positionPairs.begin(), positionPairs.end(), searchKey);
                            
            //                 if (itRc != positionPairs.end() && *itRc == searchKey) {
            //                     uint64_t rcMemberIdx = itRc - positionPairs.begin();
            //                     uint64_t rcClusterId = disjointSets.find(rcMemberIdx);
                                
            //                     // Mark RC cluster as compatible too
            //                     if (rcClusterId != originalClusterId) {
            //                         __sync_bool_compare_and_swap(&variantClusteringValidClustersCompatible[rcClusterId], 0, 1);
            //                     }
            //                 }
            //             }
            //         }

            //         // Collect reads
            //         for (const auto& rp : nodeReads[nodeIdx]) {
            //             if (rp.isPhase1 == false) {
            //                 // Supported Target Allele -> In Phase
            //                 inPhaseReads.push_back(rp.orientedReadIdValue);
            //             } else {
            //                 // Supported Alt Allele -> Out of Phase
            //                 outOfPhaseReads.push_back(rp.orientedReadIdValue);
            //             }
            //         }
            //     }
            // }

            // // Sort & Unique for set operations
            // std::sort(inPhaseReads.begin(), inPhaseReads.end());
            // inPhaseReads.erase(std::unique(inPhaseReads.begin(), inPhaseReads.end()), inPhaseReads.end());

            // std::sort(outOfPhaseReads.begin(), outOfPhaseReads.end());
            // outOfPhaseReads.erase(std::unique(outOfPhaseReads.begin(), outOfPhaseReads.end()), outOfPhaseReads.end());

            // // Final Set: InPhase MINUS OutOfPhase
            // // A read is only trusted if it NEVER supported the Alt allele in any valid cluster.
            



            // // Calculate total weight (total number of nodes in all valid chains)
            // // This is the number of sites that were compatible and were used to cast votes
            // uint32_t totalWeight = 0;
            // for (const auto& chain : validChains) {
            //     totalWeight += (uint32_t)chain.size();
            // }

            // for (const auto& chain : validChains) {
            //     for (int nodeIdx : chain) {
            //         // Mark cluster as valid... (already done above)

            //         // Collect reads
            //         for (const auto& rp : nodeReads[nodeIdx]) {
            //             if (rp.isPhase1 == false) {
            //                 // Supported Target Allele -> In Phase
            //                 inPhaseReads.push_back(rp.orientedReadIdValue);
            //             } else {
            //                 // Supported Alt Allele -> Out of Phase
            //                 outOfPhaseReads.push_back(rp.orientedReadIdValue);
            //             }
            //         }
            //     }
            // }

            // // Sort & Unique for set operations
            // std::sort(inPhaseReads.begin(), inPhaseReads.end());
            // inPhaseReads.erase(std::unique(inPhaseReads.begin(), inPhaseReads.end()), inPhaseReads.end());

            // std::sort(outOfPhaseReads.begin(), outOfPhaseReads.end());
            // outOfPhaseReads.erase(std::unique(outOfPhaseReads.begin(), outOfPhaseReads.end()), outOfPhaseReads.end());

            // // Final Set: InPhase MINUS OutOfPhase
            // vector<OrientedReadId::Int> finalHapOrientedReadIds;
            // // Add Target Read itself (always in phase)
            // finalHapOrientedReadIds.push_back(currentOrientedReadId0.getValue());
            // std::set_difference(
            //     inPhaseReads.begin(), inPhaseReads.end(),
            //     outOfPhaseReads.begin(), outOfPhaseReads.end(),
            //     std::back_inserter(finalHapOrientedReadIds)
            // );

            // // Sort again to include the target read
            // std::sort(finalHapOrientedReadIds.begin(), finalHapOrientedReadIds.end());
            // finalHapOrientedReadIds.erase(std::unique(finalHapOrientedReadIds.begin(), finalHapOrientedReadIds.end()), finalHapOrientedReadIds.end());
            

            // // Print finalHapOrientedReadIds with one orientedReadId per line if readId is 0
            // if(currentOrientedReadId0.getReadId() == 87 and currentOrientedReadId0.getStrand() == 0) {
            //     cout << "Final Haplotype Oriented Reads: " << endl;
            //     for(OrientedReadId::Int orientedReadIdValue : finalHapOrientedReadIds) {
            //         cout << OrientedReadId::fromValue(orientedReadIdValue).getString() << endl;
            //     }
            // }


            // // CAST VOTES DIRECTLY TO GRAPH
            // // Add votes: currentOrientedReadId -> otherOrientedReadId
            // {
            //     std::lock_guard<std::mutex> lock(haplotypeGraphMutex);
            //     for(OrientedReadId::Int otherOrientedReadIdValue : finalHapOrientedReadIds) {
            //         if(otherOrientedReadIdValue != currentOrientedReadId0.getValue()) {
            //             // Try to add edge with totalWeight
            //             auto result = boost::add_edge(currentOrientedReadId0.getValue(), otherOrientedReadIdValue, totalWeight, *globalHaplotypeGraph);
            //         }
            //     }
            // }


        }
    }
}


void Assembler::refineClustersThreadFunction(uint64_t threadId) {
    const uint64_t currentMinAlleleCoverage = this->minAlleleCoverage;
    const auto& graph = *globalHaplotypeGraph;
    
    uint64_t clusterIdBegin, clusterIdEnd;
    while (getNextBatch(clusterIdBegin, clusterIdEnd)) {
        for (uint64_t clusterId = clusterIdBegin; clusterId < clusterIdEnd; clusterId++) {
            
            // Only process valid/compatible clusters
            if (!variantClusteringValidClustersCompatible[clusterId]) {
                continue;
            }

            const auto& members = variantClusteringMembersByRepIdx[clusterId];
            
            // 1. Identify Valid Alleles (those with sufficient coverage)
            std::array<uint32_t, 5> alleleCounts = {0};
            for (uint64_t memberIdx : members) {
                if (memberIdx < variantClusteringPositionPairAlleles.size()) {
                    uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                    if (a < 5) alleleCounts[a]++;
                }
            }
            
            std::vector<uint8_t> validAlleles;
            for(uint8_t a=0; a<5; a++) {
                if (alleleCounts[a] >= currentMinAlleleCoverage) {
                    validAlleles.push_back(a);
                }
            }
            
            // If fewer than 2 valid alleles, no phasing to refine
            if (validAlleles.size() < 2) continue;

            // 2. Partition Reads by Allele using oriented read IDs
            std::array<std::unordered_set<uint32_t>, 5> groups;
            for (uint64_t memberIdx : members) {
                if (memberIdx < variantClusteringPositionPairAlleles.size()) {
                    uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                    if (a < 5) {
                        const auto& pp = variantClusteringPositionPairs[memberIdx];
                        groups[a].insert(pp.first.getValue());
                    }
                }
            }

            // 3. Process each allele group separately
            for (uint8_t allele = 0; allele < 5; allele++) {
                const auto& group = groups[allele];
                
                // If group is too small, we can't really do intersection of 2 reads + others.
                // If size < 3, we just keep them (default status 0).
                if(group.size() < 3) {
                    continue;
                }

                // Store (ReadVal, Weight) pairs
                std::vector<std::pair<uint32_t, uint32_t>> readWeights;
                readWeights.reserve(group.size());

                // Access weight map
                auto weightMap = boost::get(boost::edge_weight, graph);

                for(uint32_t readVal : group) {
                    uint32_t weight = 0;
                    // Check if vertex exists in graph
                    if (readVal >= boost::num_vertices(graph)) {
                        readWeights.push_back({readVal, 0}); 
                    } else {
                        // Get weight from first outgoing edge
                        auto outEdges = boost::out_edges(readVal, graph);
                        if(outEdges.first != outEdges.second) {
                            auto edge = *outEdges.first;
                            weight = weightMap[edge];
                        }
                    }
                    readWeights.push_back({readVal, weight});
                }

                // Sort by Weight Descending
                std::sort(readWeights.begin(), readWeights.end(), 
                    [](const std::pair<uint32_t, uint32_t>& a, const std::pair<uint32_t, uint32_t>& b) {
                        return a.second > b.second;
                    });

                // If the top 1 has 0 weight (no edges), keep all.
                if (readWeights[0].second == 0) {
                    continue;
                }

                // Select Top 1
                uint32_t r1 = readWeights[0].first;
                
                // Filter
                for(const auto& rw : readWeights) {
                    uint32_t readVal = rw.first;
                    
                    // Keep Top 1
                    if (readVal == r1) {
                        continue; 
                    }

                    // Check connectivity to R1
                    bool connectedToR1 = false;

                    if (r1 < boost::num_vertices(graph) && readVal < boost::num_vertices(graph)) {
                        
                        // Check R1 neighbors
                        auto neighborsR1 = boost::adjacent_vertices(r1, graph);
                        for(auto it = neighborsR1.first; it != neighborsR1.second; ++it) {
                            if(*it == readVal) { connectedToR1 = true; break; }
                        }
                        if (!connectedToR1) {
                            auto inEdgesR1 = boost::in_edges(r1, graph);
                            for(auto it = inEdgesR1.first; it != inEdgesR1.second; ++it) {
                                if(boost::source(*it, graph) == readVal) { connectedToR1 = true; break; }
                            }
                        }
                    }

                    if (connectedToR1) {
                        // Keep
                    } else {
                        // Mark as Stray
                        // Need to find memberIdx
                        for(uint64_t memberIdx : members) {
                            if (memberIdx < variantClusteringPositionPairAlleles.size()) {
                                if (variantClusteringPositionPairAlleles[memberIdx] == allele && 
                                    variantClusteringPositionPairs[memberIdx].first.getValue() == readVal) {
                                    variantClusteringMemberStatus[memberIdx] = 1; // Mark as Stray
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
