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
#include <cstdint>

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
        
        // Helper lambda to check if a read is phased and voting for the other
        auto checkPhasing = [&](ReadId sourceRead, ReadId targetRead, bool sameStrand) -> std::pair<bool, bool> {
            OrientedReadId::Int u = OrientedReadId(sourceRead, 0).getValue();
            if (u >= numVertices || boost::out_degree(u, graph) == 0) return {false, false}; // Not Phased
            
            // Is Phased, check for edge with Weight 1 (CIS)
            OrientedReadId::Int v_target = OrientedReadId(targetRead, sameStrand ? 0 : 1).getValue();
            auto outEdges = boost::out_edges(u, graph);
            for(auto it = outEdges.first; it != outEdges.second; ++it) {
                // We check if the target matches AND the weight is 1 (CIS).
                // Edges with weight 2 (trans) must NOT count for keeping the alignment.
                if (boost::target(*it, graph) == v_target) {
                    uint32_t w = boost::get(boost::edge_weight, graph, *it);
                    if (w == 1) return {true, true}; // Phased & Found (CIS)
                }
            }
            return {true, false}; // Phased & Not Found (or Found TRANS only, which counts as Not Found for keeping)
        };

        // Check Read 0 -> Read 1
        std::pair<bool, bool> r0 = checkPhasing(ad.readIds[0], ad.readIds[1], ad.isSameStrand);
        // Check Read 1 -> Read 0
        std::pair<bool, bool> r1 = checkPhasing(ad.readIds[1], ad.readIds[0], ad.isSameStrand);

        bool keep = true;
        
        // If BOTH are phased, we require at least one confirmation.
        // If one is unphased (no het sites), we keep the overlap.
        if (r0.first && r1.first) {
            keep = (r0.second || r1.second);
        }
        // Else (neither phased) -> Keep (default)

        // Filter out alignments involving palindromic reads
        if (reads->getFlags(ad.readIds[0]).isPalindromic || reads->getFlags(ad.readIds[1]).isPalindromic) {
            keep = false;
        }

        if (keep) {
            alignmentData[alignmentId].info.isInReadGraph = 1;
            keptAlignmentCount++;
        } else {
            keepAlignment[alignmentId] = false;
            alignmentData[alignmentId].info.isInReadGraph = 0;
        }
    }

    cout << timestamp << "Kept " << keptAlignmentCount << " / " << alignmentCount << " alignments after haplotype filtering." << endl;

    // --- Old Chimeric Detection (Kept per user request) ---
    detectChimericReads(threadCount);
    rescueChimericReads(threadCount);

    // --- New Filtering Pipeline (Mirrors ma_hit_sub -> detect_chimeric -> ma_hit_cut -> ma_hit_flt) ---
    cout << timestamp << "Running read graph filtering pipeline..." << endl;

    // 1. ma_hit_sub: Filter local segments (coverage > 3)
    // Mirrors: ma_hit_sub(min_dp, src, n_read, readLen, mini_overlap_length, cov);
    // We use minCoverage = 3 as standard default.
    filterLocalSegments(3, threadCount);

    // 2. detect_chimeric_reads: Identify and filter chimeric reads using anchor logic
    // Mirrors: detect_chimeric_reads(src, n_read, readLen, *cov, asm_opt.max_ov_diff_final*2.0, ul, UL_COV_THRES);
    // shiftRate = 0.05 (approx max_ov_diff_final * 2), ulThres = 0 (unlimited)
    detectChimericReadsFromAnchors(0.05, 0, threadCount);

    // 3. ma_hit_cut: Clip alignments to valid regions
    // Mirrors: ma_hit_cut(src, n_read, readLen, mini_overlap_length, cov);
    // minOverlapLength = 0 (or a small value like 100 bp if desired)
    applyCoverageCuts(0, threadCount);

    // 4. ma_hit_flt: Filter hanging overlaps (dovetail check)
    // Mirrors: ma_hit_flt(src, n_read, *cov, max_hang_length, mini_overlap_length);
    // maxHang = 1000, maxHangRate = 0.8, minOverlap = 0
    filterHangingOverlaps(1000, 0.8, 0, threadCount);

    // // 5. ma_hit_contained_advance: Remove contained reads
    // // Mirrors: ma_hit_contained_advance(...);
    // // Uses the same maxHang (1000) and minOverlap (0).
    // removeContainedReads(1000, 0.8, 0, threadCount);

    // Sync keepAlignment status with isDeleted flags
    // The functions above mark alignments as `isDeleted`. We need to respect that.
    uint64_t filteredCount = 0;
    #pragma omp parallel for reduction(+:filteredCount)
    for(uint64_t alignmentId=0; alignmentId<alignmentCount; alignmentId++) {
        if (!keepAlignment[alignmentId]) continue; 

        const auto& ad = alignmentData[alignmentId];
        // If marked deleted by any stage of the pipeline, drop it.
        // Also check if read is chimeric (though detectChimericReadsFromAnchors handles deleting edges, safe to double check).
        if (ad.isDeleted || isChimericRead[ad.readIds[0]] || isChimericRead[ad.readIds[1]]) {
             keepAlignment[alignmentId] = false;
             alignmentData[alignmentId].info.isInReadGraph = 0;
             filteredCount++;
        }
    }
    cout << timestamp << "Pipeline removed " << filteredCount << " additional alignments." << endl;

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

            // Single pass through all members - WITH DEDUPLICATION (Best Hit Policy)
            // To prevent double counting from multiple alignments (e.g. both strands),
            // we ensure each ReadId contributes at most once.
            // Since we can't easily access alignment score here, we take the "first" one encountered
            // (which is arbitrary but prevents bias inflation).
            
            // Fixed-size optimization: Use small vector for seen ReadIds to avoid set overhead
            // Most clusters have < 100 coverage.
            std::vector<ReadId> seenReads;
            seenReads.reserve(members.size());

            for (uint64_t memberIdx : members) {
                const auto& pp = variantClusteringPositionPairs[memberIdx];
                const OrientedReadId orientedReadId = pp.first;
                
                // Deduplication check
                bool seen = false;
                for(ReadId r : seenReads) {
                    if (r == orientedReadId.getReadId()) {
                        seen = true;
                        break;
                    }
                }
                if (seen) continue; // Skip duplicate contribution from same read
                seenReads.push_back(orientedReadId.getReadId());

                const uint32_t position = pp.second;
                
                const uint8_t allele = variantClusteringPositionPairAlleles[memberIdx];
                if (allele < 5) {
                    // Check for palindromic reads and skip them for coverage counting
                    // This "reduces the coverage" as requested for checking potential het sites.
                    if (!reads->getFlags(orientedReadId.getReadId()).isPalindromic) {
                        alleleCounts[allele]++;
                        // Track strand counts per allele
                        if (orientedReadId.getStrand() == 0) {
                            strand0Counts[allele]++;
                        } else {
                            strand1Counts[allele]++;
                        }
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
            
            // --- Step 1: Generate Virtual Clusters & Populate Reads (Nodes for DP) ---
            // A physical cluster with alleles {Target, A, B} becomes 
            // two virtual clusters: (Target vs A) and (Target vs B).
            // Each VirtualCluster now self-contains its evidence (reads).

            struct ReadPhase {
                OrientedReadId::Int orientedReadIdValue; // Using the raw integer value
                bool isPhase1; // false=Target, true=Alt
            };

            struct VirtualCluster {
                uint64_t clusterId;
                uint8_t targetAllele;
                uint8_t altAllele;
                std::vector<ReadPhase> reads; // Evidence
                
                // Haplotype Assignment (Step 5) Fields
                uint32_t position = 0;
                uint32_t occ_0 = 0; 
                uint32_t occ_1 = 0;
                int32_t score = 0; // 0=Unconfirmed, 1=Confirmed Het, -1=Invalidated
                bool is_filtered = false;
            };
            vector<VirtualCluster> dpNodes;
            
            const uint64_t currentMinAlleleCoverage = this->minAlleleCoverage; 

            // Find clusters that involve this read on strand 0.
            while(it != positionPairs.end() && it->first == currentOrientedReadId0) {
                uint64_t indexInPositionPairs = it - positionPairs.begin();
                uint64_t clusterId = disjointSets.find(indexInPositionPairs);
                
                // Identify Target Allele (The allele supported by the current read we are phasing)
                uint8_t targetAllele = 255;
                if (indexInPositionPairs < variantClusteringPositionPairAlleles.size()) {
                    targetAllele = variantClusteringPositionPairAlleles[indexInPositionPairs];
                }
                
                // If allele is invalid (e.g. filtered out duplicate), skip this position.
                if (targetAllele >= 5) {
                    it++;
                    continue;
                }

                // Count alleles in this cluster
                // ALSO count Strand-0 (Forward) occurrences for bias checking.
                const auto& members = membersByRepIdx[clusterId];
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

                // Check Strand Bias & Coverage for the Reference Allele
                // If the Reference allele fails, we cannot use it as a pivot.
                if (!isStrandBiased(alleleCounts[targetAllele], alleleStrand0Counts[targetAllele])) {
                     // Check if the Reference Allele has sufficient coverage
                    if (alleleCounts[targetAllele] >= currentMinAlleleCoverage) {
                        
                        // 4. Identify Valid Chains (Graph Nodes)
                        std::vector<std::vector<int>> validChains;
                        std::vector<int> validChainsDumpster; // Captured sites that failed filters
                        std::vector<bool> visited(dpNodes.size(), false);
                        // Identify valid Alt alleles
                        std::vector<uint8_t> validAlts;
                        for(uint8_t a=0; a<5; a++) {
                            if (a == targetAllele) continue;

                            // Check Alt Coverage
                            if (alleleCounts[a] < currentMinAlleleCoverage) continue;

                            // // Check Alt Strand Bias
                            // if (isStrandBiased(alleleCounts[a], alleleStrand0Counts[a])) continue;

                            validAlts.push_back(a);
                        }

                        if (!validAlts.empty()) {
                            // Create Virtual Clusters
                            size_t startIdx = dpNodes.size();
                            for (uint8_t alt : validAlts) {
                                VirtualCluster vc;
                                vc.clusterId = clusterId;
                                vc.targetAllele = targetAllele;
                                vc.altAllele = alt;
                                vc.reads.reserve(members.size()); // Upper bound
                                dpNodes.push_back(std::move(vc));
                            }

                            // POPULATE READS for these new virtual clusters
                            // We iterate members once and distribute to relevant VCs
                            for(uint64_t memberIdx : members) {
                                const auto& pp = positionPairs[memberIdx];
                                OrientedReadId orientedReadId = pp.first;
                                
                                // Skip target oriented read itself
                                if(orientedReadId == currentOrientedReadId0) continue;

                                if(memberIdx < variantClusteringPositionPairAlleles.size()) {
                                    uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                                    
                                    if (a == targetAllele) {
                                        // Supported Target: Add (Phase 0) to ALL VCs derived from this physical cluster
                                        for(size_t k=0; k<validAlts.size(); k++) {
                                            dpNodes[startIdx + k].reads.push_back({orientedReadId.getValue(), false});
                                        }
                                    } else {
                                        // Supported Alt: Check if it matches one of our valid Alts
                                        for(size_t k=0; k<validAlts.size(); k++) {
                                            if (a == validAlts[k]) {
                                                dpNodes[startIdx + k].reads.push_back({orientedReadId.getValue(), true});
                                                break; // Found the matching alt VC
                                            }
                                        }
                                    }
                                }
                            }

                            // SORT READS for these new virtual clusters
                            for(size_t k=0; k<validAlts.size(); k++) {
                                auto& reads = dpNodes[startIdx + k].reads;
                                std::sort(reads.begin(), reads.end(), 
                                    [](const ReadPhase& a, const ReadPhase& b) { 
                                        return a.orientedReadIdValue < b.orientedReadIdValue; 
                                    });
                            }
                        }
                    }
                }

                ++it;
            }

            // Get the number of Virtual Clusters (DP nodes)
            const size_t N = dpNodes.size();

            // If no DP nodes, skip this read
            if (N == 0) continue;

            // --- Step 2: Strict DP for Phasing Consistency (gen_rphase_dp) ---
            // Replaces previous LCG logic with strict "comput_sc_rphase" consistency check.
            
            // DP definition:
            // dp[i] = max score of a consistent chain ending at node i
            // parent[i] = predecessor node index
            // Initialize with 1 because every node is a valid chain of length 1 by itself.
            std::vector<int64_t> dpScore(N, 1); 
            std::vector<int> parent(N, -1);

            // DP Loop
            for(int i = 0; i < N; i++) {
                int64_t currentMaxScore = 1; // Base score (chain of length 1)
                int bestPredecessorIndex = -1;

                // Optimization: Iterate backwards to find closest predecessor first.
                for(int j = i - 1; j >= 0; j--) {
                    // Optimization: If the best possible score from j doesn't improve our current max, skip.
                    // Since we iterate backwards, we likely found a good/closest predecessor early.
                    if (dpScore[j] + 1 <= currentMaxScore) continue;

                    // 1. Physical Cluster constraint
                    // Only allow chains between different clusters (don't chain multiallelic clusters).
                    if (dpNodes[i].clusterId == dpNodes[j].clusterId) continue;

                    // 2. Strict Consistency Check (comput_sc_rphase)
                    // Check if the reads support the same phase.
                    const auto& readsI = dpNodes[i].reads;
                    const auto& readsJ = dpNodes[j].reads;

                    int supportTarget = 0; // Support matching Ref-Ref
                    int supportAlt = 0; // Support matching Alt-Alt
                    bool inconsistent = false;

                    // Sorted intersection
                    auto itI = readsI.begin();
                    auto itJ = readsJ.begin();
                    auto endI = readsI.end();
                    auto endJ = readsJ.end();

                    // We compare dpNodes[i].reads and dpNodes[j].reads.
                    // Since they are sorted, we use a single linear pass (like merging two sorted lists).
                    // Conflict: If any orientedReadId is shared but has different phases (Ref vs Alt), we flag inconsistent = true and abort.
                    // Support: If phases match, we count support for Ref (supportTarget) and Alt (supportAlt).
                    while(itI != endI && itJ != endJ) {
                        if(itI->orientedReadIdValue < itJ->orientedReadIdValue) {
                            ++itI;
                        } else if(itJ->orientedReadIdValue < itI->orientedReadIdValue) {
                            ++itJ;
                        } else {
                            // Shared Read
                            if (itI->isPhase1 != itJ->isPhase1) {
                                // Conflict! (One says Ref, one says Alt)
                                inconsistent = true;
                                break;
                            } else {
                                // Consistent
                                if (itI->isPhase1 == false) supportTarget++;
                                else supportAlt++;
                            }
                            ++itI;
                            ++itJ;
                        }
                    }

                    if (inconsistent) continue; // Score = INT64_MIN -> skip

                    // Check for positive evidence (Strict Requirement)
                    // Ensure that we have at least one read supporting the Reference haplotype (supportTarget) 
                    // AND at least one read supporting the Alternative haplotype (supportAlt) connecting these two sites.
                    // Just because two sites don't conflict (e.g., no reads jumping from Ref to Alt) doesn't mean 
                    // they are truly linked. They might just be far apart with no shared reads, or only one haplotype might be covered.
                    // This eliminates "hollow" chains where one haplotype is actually missing or unsupported.
                    if (supportTarget > 0 && supportAlt > 0) {
                        // Consistent and Linked!
                        int64_t candidateScore = 1 + dpScore[j];
                        // Maximization: We check if this new path coming from $j$ is longer than 
                        // any best path we found from other predecessors so far (currentMaxScore).
                        // If it is longer, we update currentMaxScore (new best score) and record bestPredecessorIndex (best predecessor) 
                        // to reconstruct the path later.
                        if (candidateScore > currentMaxScore) {
                            currentMaxScore = candidateScore;
                            bestPredecessorIndex = j;
                        }
                    }
                }
                dpScore[i] = currentMaxScore;
                parent[i] = bestPredecessorIndex;
            }

            // --- Step 3: Traceback & Selection ---
            // Extract chains starting from best scores.
            
            std::vector<int> sortedNodeIndices(N);
            std::iota(sortedNodeIndices.begin(), sortedNodeIndices.end(), 0);
            std::sort(sortedNodeIndices.begin(), sortedNodeIndices.end(), [&](int a, int b) {
                return dpScore[a] > dpScore[b]; // Descending score
            });

            vector<bool> isAssigned(N, false);
            vector<vector<int>> validChains;
            vector<int> validChainsDumpster; // Captured sites that failed filters

            // Filter Helper Definition: Checks if the site is dominated by homopolymer errors.
            // Returns true if valid coverage (non-HP reads) is too low for either allele.
            auto isHomopolymerDominated = [&](const VirtualCluster& node) -> bool {
                int64_t validCoverageTarget = 0;
                int64_t validCoverageAlt = 0;
                
                const auto& members = membersByRepIdx[node.clusterId];
                const auto& reads = this->getReads();

                for(uint64_t memberIdx : members) {
                    if(memberIdx >= variantClusteringPositionPairAlleles.size()) continue;

                    uint8_t a = variantClusteringPositionPairAlleles[memberIdx];
                    bool isTarget = (a == node.targetAllele);
                    bool isAlt = (a == node.altAllele);

                    if (!isTarget && !isAlt) continue; // Skip noise/irrelevant alleles

                    // Check Homopolymer Status
                    const auto& pp = positionPairs[memberIdx];
                    ReadId rId = pp.first.getReadId();
                    uint32_t pos = pp.second; 
                    
                    dinara::LongBaseSequenceView readSequence = reads.getRead(rId);
                    uint64_t checkPos = pos;
                    if (pp.first.getStrand() == 1) checkPos = readSequence.baseCount - 1 - pos;
                    
                    if (isSiteInHomopolymerRegion(checkPos, readSequence)) {
                        continue; // Ignored (HP region)
                    }

                    if (isTarget) validCoverageTarget++;
                    else validCoverageAlt++;
                }
                
                // Check thresholds (Min 3 non-HP reads for both)
                if (validCoverageTarget < 3 || validCoverageAlt < 3) return true; 
                return false; 
            };

            for (int idx : sortedNodeIndices) {
                if (isAssigned[idx]) continue;

                // Traceback to reconstruct the chain
                std::vector<int> chain;
                int curr = idx;
                while (curr != -1 && !isAssigned[curr]) {
                    isAssigned[curr] = true;
                    chain.push_back(curr);
                    curr = parent[curr];
                }
                if (chain.empty()) continue;

                // Chain is [End, ..., Start]. Reverse to [Start, ..., End]
                std::reverse(chain.begin(), chain.end());

                // Evaluate Chain Quality
                if (chain.size() > 1) {
                    // MULTI-NODE CHAIN:
                    // Only keep sites that meet coverage requirements within the chain
                    std::vector<int> filteredChain;
                    for(int nodeIdx : chain) {
                        // Check Target Support coverage.
                        // We use the pre-filtered reads vector for efficiency.
                        int supportTarget = 0;
                        for(const auto& rp : dpNodes[nodeIdx].reads) if(!rp.isPhase1) supportTarget++;
                        
                        if (supportTarget >= minMultiNodeChainSupport) {
                            filteredChain.push_back(nodeIdx);
                        } else {
                            // DUMPSTER: Low support sites in chains
                            validChainsDumpster.push_back(nodeIdx);
                        }
                    }
                    if (!filteredChain.empty()) validChains.push_back(filteredChain);

                } else {
                    // ISOLATED SITE:
                    int nodeIdx = chain[0];
                    int supportTarget = 0;
                    for(const auto& rp : dpNodes[nodeIdx].reads) if(!rp.isPhase1) supportTarget++;
                    
                    if (supportTarget >= minIsolatedSiteSupport) {
                        // Strict Homopolymer Check for isolated sites
                        if (!isHomopolymerDominated(dpNodes[nodeIdx])) {
                            validChains.push_back({nodeIdx});
                        } else {
                            // DUMPSTER: Homopolymer failure
                            validChainsDumpster.push_back(nodeIdx); 
                        }
                    } else {
                        // DUMPSTER: Low support isolated
                        validChainsDumpster.push_back(nodeIdx);
                    }
                }
            }


            // --- Step 5: Haplotype Assignment (Naive HiFi) ---
            
            // Data structures
            // 'SiteStat' is now merged into 'VirtualCluster' (dpNodes)
            // We only need Read structures locally.

            struct ReadSiteInfo {
                uint32_t dpNodeIdx; // Index in 'dpNodes'
                bool isPhase1;      // false=Ref/Match, true=Alt/Mismatch
            };

            struct ReadStat {
                OrientedReadId::Int orientedReadIdValue;
                std::vector<ReadSiteInfo> coveredSites;
                uint32_t informative_score; // 'o'
                uint8_t is_trans;           // 0=Unknown, 1=CIS, 2=TRANS
                bool strong;
            };

            // 1. Data Preparation
            std::vector<ReadStat> readStats;
            std::unordered_map<OrientedReadId::Int, size_t> readMap;

            // Collect unique sites
            std::set<int> usedNodeIndices;
            for (const auto& chain : validChains) {
                for (int nodeIdx : chain) usedNodeIndices.insert(nodeIdx);
            }

            // Mark Compatible Clusters
            for (int nodeIdx : usedNodeIndices) {
                uint64_t cId = dpNodes[nodeIdx].clusterId;
                __sync_bool_compare_and_swap(&variantClusteringValidClustersCompatible[cId], 0, 1);

                // Mark RC
                if (!membersByRepIdx[cId].empty()) {
                    const auto& pp = positionPairs[membersByRepIdx[cId][0]];
                    ReadId readId = pp.first.getReadId();
                    Strand strand = pp.first.getStrand();
                    uint32_t position = pp.second;
                    Strand rcStrand = 1 - strand;
                    OrientedReadId rcId(readId, rcStrand);
                    uint32_t rcPos = uint32_t(getReads().getReadRawSequenceLength(readId) - 1 - position);
                    auto searchKey = std::make_pair(rcId, rcPos);
                    auto itRc = std::lower_bound(positionPairs.begin(), positionPairs.end(), searchKey);
                    if (itRc != positionPairs.end() && *itRc == searchKey) {
                        uint64_t rcCId = disjointSets.find(itRc - positionPairs.begin());
                        if (rcCId != cId) __sync_bool_compare_and_swap(&variantClusteringValidClustersCompatible[rcCId], 0, 1);
                    }
                }
            }

            // Init Sites (in dpNodes) & Reads
            for (int nodeIdx : usedNodeIndices) {
                // Initialize VirtualCluster fields
                dpNodes[nodeIdx].occ_0 = 0;
                dpNodes[nodeIdx].occ_1 = 0;
                dpNodes[nodeIdx].score = 0;
                dpNodes[nodeIdx].is_filtered = false;
                
                uint64_t cId = dpNodes[nodeIdx].clusterId;
                if (!membersByRepIdx[cId].empty()) {
                    dpNodes[nodeIdx].position = positionPairs[membersByRepIdx[cId][0]].second;
                } else {
                    dpNodes[nodeIdx].position = 0; 
                }

                // Process Reads
                for (const auto& rp : dpNodes[nodeIdx].reads) {
                    if (readMap.find(rp.orientedReadIdValue) == readMap.end()) {
                        size_t idx = readStats.size();
                        readStats.push_back(ReadStat{rp.orientedReadIdValue, {}, 0, 1, false}); // Default is_trans=1 (CIS)
                        readMap[rp.orientedReadIdValue] = idx;
                    }
                    size_t rIdx = readMap[rp.orientedReadIdValue];
                    readStats[rIdx].coveredSites.push_back(ReadSiteInfo{(uint32_t)nodeIdx, rp.isPhase1});
                    
                    if (!rp.isPhase1) dpNodes[nodeIdx].occ_0++;
                    else dpNodes[nodeIdx].occ_1++;
                }
            }

            // 2.1 Filter Adjacent SNPs (Hifiasm "Block-based" Logic)
            // Goal: Remove ANY site that is exactly 1bp away from another site.
            // Why Block-Based? 
            // - We may have multiple nodes at the same position (multi-allelic sites or duplicates).
            // - These nodes form a "block" [l, k) in the sorted list.
            // - We must check if the ENTIRE block is adjacent to the *previous* unique position or *next* unique position.
            // - If adjacency is found, ALL nodes in the block are filtered out.

            std::vector<int> sortedNodeIndicesByPos(usedNodeIndices.begin(), usedNodeIndices.end());
            std::sort(sortedNodeIndicesByPos.begin(), sortedNodeIndicesByPos.end(), [&](int a, int b){
                return dpNodes[a].position < dpNodes[b].position;
            });

            size_t n = sortedNodeIndicesByPos.size();
            size_t l = 0; // Start index of current block
            for (size_t k = 1; k <= n; ++k) {
                // Determine if we reached the end of a block of identical positions:
                // k == n: End of list
                // pos[k] != pos[l]: New unique position found
                if (k == n || dpNodes[sortedNodeIndicesByPos[k]].position != dpNodes[sortedNodeIndicesByPos[l]].position) {
                    bool filterBlock = false;
                    uint32_t currentPos = dpNodes[sortedNodeIndicesByPos[l]].position;

                    // 1. Check Previous Block Adjacency
                    // If this isn't the first block (l > 0), look at the last element of the PREVIOUS block (l-1)
                    if (l > 0) {
                        uint32_t prevPos = dpNodes[sortedNodeIndicesByPos[l-1]].position;
                        if (currentPos == prevPos + 1) {
                            filterBlock = true; // Block is 1bp after previous block
                        }
                    }

                    // 2. Check Next Block Adjacency
                    // If this isn't the last block (k < n), look at the first element of the NEXT block (k)
                    if (k < n) {
                        uint32_t nextPos = dpNodes[sortedNodeIndicesByPos[k]].position;
                        if (currentPos + 1 == nextPos) {
                            filterBlock = true; // Block is 1bp before next block
                        }
                    }

                    // 3. Apply Filter
                    // If either neighbor was exactly 1bp away, invalidate EVERY node in this current block.
                    if (filterBlock) {
                        for (size_t i = l; i < k; ++i) {
                            dpNodes[sortedNodeIndicesByPos[i]].is_filtered = true;
                        }
                    }

                    // Advance start pointer 'l' to 'k' (start of the next block)
                    l = k;
                }
            }

            // 2.2 Calculate Informative Scores (Hifiasm "o" metric)
            // Goal: Score each read based on how many "strong" mismatches it provides.
            // This score ('o') is used to process the most informative reads first.
            for (auto& r : readStats) {
                r.informative_score = 0;
                for (const auto& info : r.coveredSites) {
                    VirtualCluster& s = dpNodes[info.dpNodeIdx];
                    if (s.is_filtered) continue;
                    
                    // Criteria 1: Must be a Mismatch (Phase 1 / Alt)
                    // Matches (Reference allele) do not help distinguish the *other* haplotype in this logic.
                    if (!info.isPhase1) continue; 

                    // Criteria 2: Site must have minimum coverage (occ >= 2)
                    // If a site is poorly covered on either allele, it's too noisy to trust.
                    if (s.occ_0 < 2 || s.occ_1 < 2) continue;
                    
                    // Criteria 3: High-Quality / "Strong" Site (occ >= 3)
                    // We only increment the informative score if the site coverage is robust (>= 3).
                    // This creates the 'o' score used for sorting.
                    if (s.occ_0 >= 3 && s.occ_1 >= 3) {
                        r.informative_score++;
                    }
                }
            }

            // 2.3 Sort reads by their informative_score in descending order.
            // Reads with higher informative_score will be processed first.
            // We want to process the reads with the strongest evidence of being "Other Haplotype" (TRANS) first. 
            // Identifying these correctly clarifies the picture for ambiguous reads later.
            std::vector<uint32_t> readIndices(readStats.size());
            std::iota(readIndices.begin(), readIndices.end(), 0);
            std::sort(readIndices.begin(), readIndices.end(), [&](uint32_t a, uint32_t b) {
                return readStats[a].informative_score > readStats[b].informative_score;
            });

            // 2.4 First Loop: Greedy TRANS Assignment & Coverage Update
            // Goals:
            // 1. Identify definitive TRANS (Other Haplotype) reads based on their score.
            // 2. "Lock in" the sites they mismatch as Confirmed Heterozygous (score=1).
            // 3. Remove their support from Reference alleles (since we know they aren't Reference).
            for (uint32_t rIdx : readIndices) {
                ReadStat& r = readStats[rIdx];
                
                // Hifiasm Logic: We re-calculate 'o' with CURRENT site counts.
                // Site counts change as we process reads, so a read that looked informative initially 
                // might lose that status if its supporting sites get "downgraded" to not informative.
                // We verify if the read *still* has strong support after previous decrements.
                int o = 0;
                for (const auto& info : r.coveredSites) {
                    VirtualCluster& s = dpNodes[info.dpNodeIdx];
                    if (s.is_filtered) continue;
                    if (!info.isPhase1) continue; // Only count Mismatches
                    if (s.occ_0 < 2 || s.occ_1 < 2) continue; // Weak sites ignored
                    if (s.occ_0 >= 3 && s.occ_1 >= 3) o++;    // Strong mismatch count
                }

                // If read has no strong mismatches anymore, skip it (leave as potential CIS/Ambiguous).
                if (o == 0) continue;

                // Action 1: Mark as TRANS (Other Haplotype)
                r.is_trans = 2; // 2 = Definitely TRANS

                // Action 2: Update Sites based on this read's assignment
                for (const auto& info : r.coveredSites) {
                    VirtualCluster& s = dpNodes[info.dpNodeIdx];
                    if (s.is_filtered) continue;

                    if (info.isPhase1) { 
                        // Mismatch -> This TRANS read confirms this site is a real Difference.
                        // Mark site as "Confirmed Heterozygous" (Anchor).
                        s.score = 1; 
                    } else { 
                        // Match -> This TRANS read happens to match Reference here (homozygosity/noise).
                        // Since this read is TRANS, its vote for Reference is invalid.
                        // Decrement occ_0 to remove this false support.
                        if (s.occ_0 > 0) s.occ_0--;
                    }
                }
            }

            // 2.5 Second Loop: "Catch-up" Pass
            // Goal: Rescue ambiguous reads using the anchors established in the First Pass.
            // Why? Some reads might not have had enough strong sites to be called TRANS initially ('informative' count too low).
            // However, if they conflict with a site that is NOW "Confirmed Heterozygous" (score=1), 
            // that is strong evidence they belong to the other haplotype.
            // Logic: "If you clash with a trusted anchor, you are TRANS."
            for (uint32_t rIdx : readIndices) {
                ReadStat& r = readStats[rIdx];
                
                // Skip reads already definitively assigned as TRANS
                if (r.is_trans == 2) continue;

                int o = 0;
                for (const auto& info : r.coveredSites) {
                    VirtualCluster& s = dpNodes[info.dpNodeIdx];
                    if (s.is_filtered) continue;
                    if (!info.isPhase1) continue; // Mismatch only
                    if (s.occ_0 < 2 || s.occ_1 < 2) continue; // Min coverage check
                    
                    // Critical Check: Does this mismatch happen at a CONFIRMED Anchor?
                    // s.score was set to 1 in the First Pass.
                    if (s.score == 1) o++;
                }

                // If at least one Confirmed Mismatch exists -> Mark as TRANS
                if (o > 0) {
                    r.is_trans = 2; 
                }
            }

            // 2.6 Reset/Invalidate Contradictory Sites in CIS Reads
            // Goal: Clean up sites that show contradictory evidence in the reads we believe are CIS.
            // Logic:
            // - At this point, any read with is_trans != 2 is tentatively considered "CIS" (In-Phase).
            // - If a CIS read has a Mismatch (Phase 1) at a site, that is a contradiction.
            //   (A CIS read should technically Match indices, not Mismatch).
            // - Instead of discarding the read, we blame the site: strict logic says this site is unreliable.
            // - Action: Invalidate the site (score = -1) so it cannot determine future assignments.
            for (auto& r : readStats) {
                if (r.is_trans == 1) { // Currently considered CIS
                    for (const auto& info : r.coveredSites) {
                        if (info.isPhase1) { // Mismatch at this site
                            // Contradiction found -> Invalidate site
                            dpNodes[info.dpNodeIdx].score = -1; 
                        }
                    }
                }
            }

            // 2.7 Final Assignment & Strong Marking
            // Goal: Finalize TRANS/CIS status and identify "Strong" (high-confidence) reads.
            // Why? Even after previous loops, some CIS reads might still technically conflict 
            // with sites that are definitely Heterozygous anchors. We need a final safety check.
            for (auto& r : readStats) {
                if (r.is_trans == 2) {
                    // TRANS reads are already confirmed strong by definition (they drove the assignment).
                    r.strong = true;
                }
                else if (r.is_trans == 1) { // Tentative CIS
                    // Safety Check: Does this CIS read actually Mismatch a trusted Anchor?
                    // Anchor definition here: score=1 (Confirmed Het) AND robust coverage (>=2) on both alleles.
                    for (const auto& info : r.coveredSites) {
                        VirtualCluster& s = dpNodes[info.dpNodeIdx];
                        if (s.is_filtered) continue;
                        
                        if (s.score == 1 && s.occ_0 >= 2 && s.occ_1 >= 2) {
                            if (info.isPhase1) { // Mismatch!
                                // It contradicts a trusted anchor -> It MUST be TRANS.
                                r.is_trans = 2;
                                r.strong = true;
                                break;
                            }
                        }
                    }
                    
                    // Final "Strong" Check for surviving CIS reads:
                    // A CIS read is only "Strong" if it explicitly MATCHES at least one trusted anchor.
                    // (i.e., it's not just "not conflicting", it's positively confirmed).
                    if (r.is_trans == 1) {
                        for (const auto& info : r.coveredSites) {
                            VirtualCluster& s = dpNodes[info.dpNodeIdx];
                            if (s.score == 1 && s.occ_0 >= 2 && s.occ_1 >= 2) {
                                // Implicit: Since we passed the loop above, we know info.isPhase0 (Match).
                                r.strong = true;
                            }
                        }
                    }
                }
            }

            // 2.8 Output Generation
            // Default: if is_trans != 2, it is CIS (is_trans=1 or 0)
            
            vector<OrientedReadId::Int> finalHapOrientedReadIds;
            // Add Target Read itself (always in phase)
            finalHapOrientedReadIds.push_back(currentOrientedReadId0.getValue());

            // 2.7.5 Dumpster Diving Validation
            // Goal: Rescue evidence from sites that were discarded during chain filtering (validChainsDumpster).
            // Logic: Same consensus check. If a "Dumpster" site is effectively validated by Strong TRANS reads,
            // we can use it to flip Weak CIS reads even though the site itself was considered low quality/noisy.
            
            size_t dumpsterFlippedCount = 0;
            
            /* 
            // We need to initialize the "Dumpster" sites first (calculate position, occ, etc.)
            // Note: These sites were NOT in usedNodeIndices, so they have raw data.
            for (int dNodeIdx : validChainsDumpster) {
                // Calculation matches Init Sites loop (lines ~1069)
                dpNodes[dNodeIdx].occ_0 = 0;
                dpNodes[dNodeIdx].occ_1 = 0;
                // dpNodes[dNodeIdx].reads is already populated
                for (const auto& rp : dpNodes[dNodeIdx].reads) {
                     if (!rp.isPhase1) dpNodes[dNodeIdx].occ_0++;
                     else dpNodes[dNodeIdx].occ_1++;
                }
                // No need to map position or extensive filtering for this heuristic
            }

            for (auto& r : readStats) {
                if (r.is_trans == 1 && !r.strong) { // Weak CIS
                    // We need to check if this read covers any Dumpster sites.
                    // This is inefficient (O(Reads * DumpsterSites)), but correct for testing.
                    // Ideally, we would have added Dumpster sites to r.coveredSites, but that would pollute the main logic.
                    // Instead, we iterate the Dumpster list and check if 'r' is present in the site's read list.
                    
                    for (int dNodeIdx : validChainsDumpster) {
                        VirtualCluster& s = dpNodes[dNodeIdx];
                        // Basic quality check for the dumpster site itself (e.g. at least SOME coverage)
                        if (s.occ_0 < 2 || s.occ_1 < 2) continue; 

                        // Does 'r' mismatch here?
                        bool rMismatchesHere = false;
                        for (const auto& rp : s.reads) {
                            if (rp.orientedReadIdValue == r.orientedReadIdValue && rp.isPhase1) {
                                rMismatchesHere = true;
                        break;
                            }
                        }
                        if (!rMismatchesHere) continue;

                        // Validate this site using OTHER reads
                        bool supportedByStrongTrans = false;
                        bool conflictWithStrongReads = false;

                        for (const auto& otherRp : s.reads) {
                            if (otherRp.orientedReadIdValue == r.orientedReadIdValue) continue; // Skip self
                            if (readMap.find(otherRp.orientedReadIdValue) == readMap.end()) continue;
                            
                            size_t otherIdx = readMap[otherRp.orientedReadIdValue];
                            ReadStat& otherR = readStats[otherIdx];

                            // Check Strong TRANS opinions
                            if (otherR.is_trans == 2 && otherR.strong) {
                                if (otherRp.isPhase1) supportedByStrongTrans = true; // Confirms TRANS
                                else {
                                    conflictWithStrongReads = true; // Contradiction: Strong TRANS matches Ref
                                    break;
                                }
                            }
                            // Check Strong CIS opinions
                            else if (otherR.is_trans == 1 && otherR.strong) {
                                if (otherRp.isPhase1) {
                                    conflictWithStrongReads = true; // Contradiction: Strong CIS matches Alt
                                    break;
                                }
                            }
                        }

                        if (supportedByStrongTrans && !conflictWithStrongReads) {
                            r.is_trans = 2; // Flip to TRANS
                            dumpsterFlippedCount++;
                            break; // One valid site is enough
                        }
                    }
                }
            }
            if (dumpsterFlippedCount > 0) {
                 std::cout << "Dumpster Diving Validation: Flipped " << dumpsterFlippedCount << " Weak CIS reads to TRANS." << std::endl;
            }
            */

            // 2.8 Output Generation
            vector<std::pair<OrientedReadId::Int, int>> finalGraphEdges; // Pair<ReadVal, Weight>
            
            for (const auto& r : readStats) {
                if (r.is_trans != 2) {
                    // In-Phase (CIS) - Includes both Strong CIS and surviving Weak CIS
                    finalGraphEdges.push_back({r.orientedReadIdValue, 1}); // Weight 1 for CIS
                    finalHapOrientedReadIds.push_back(r.orientedReadIdValue);
                } 
                else if (r.is_trans == 2) {
                    // Out-Of-Phase (TRANS)
                    finalGraphEdges.push_back({r.orientedReadIdValue, 2}); // Weight 2 for TRANS
                    // Do NOT add to finalHapOrientedReadIds (helper set usually implies Same Haplotype)
                }
            }
            finalGraphEdges.push_back({currentOrientedReadId0.getValue(), 1}); // Self is CIS (Weight 1)

            // Sort & Unique (using Read ID as primary key, but we need weight consistency)
            // Just sorting by ID is enough if we trust uniqueness logic.
            std::sort(finalGraphEdges.begin(), finalGraphEdges.end());
            finalGraphEdges.erase(std::unique(finalGraphEdges.begin(), finalGraphEdges.end()), finalGraphEdges.end());

            // Also maintain the CIS-only list for debug printing if needed
            std::sort(finalHapOrientedReadIds.begin(), finalHapOrientedReadIds.end());
            finalHapOrientedReadIds.erase(std::unique(finalHapOrientedReadIds.begin(), finalHapOrientedReadIds.end()), finalHapOrientedReadIds.end());

            // Print finalHapOrientedReadIds with one orientedReadId per line if readId is 0 (debug)
            if(currentOrientedReadId0.getReadId() == 8820 and currentOrientedReadId0.getStrand() == 0) {
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
                for(const auto& edge : finalGraphEdges) {
                    OrientedReadId::Int otherVal = edge.first;
                    int weightType = edge.second; // 1=CIS, 2=TRANS
                    
                    if(otherVal != currentOrientedReadId0.getValue()) {
                        // Add edge with weight = weightType (1 or 2)
                        // Note: totalWeight (chain size) is ignored as per user instruction to use specific codes.
                        boost::add_edge(currentOrientedReadId0.getValue(), otherVal, weightType, *globalHaplotypeGraph);
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
