#ifndef DINARA_ASSEMBLER_HPP
#define DINARA_ASSEMBLER_HPP

#ifndef DINARA_ENABLE_VARIANT_CLUSTERING
#define DINARA_ENABLE_VARIANT_CLUSTERING 0
#endif

// Dinara.
#include "Alignment.hpp"
#include "AlignmentCandidates.hpp"
#include "Align6Marker.hpp"
#include "AssemblerOptions.hpp"
#include "AssemblyGraph2Statistics.hpp"
#include "HttpServer.hpp"
#include "Kmer.hpp"
#include "KmerDistributionInfo.hpp"
#include "MappedMemoryOwner.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "MarkerGraphEdgePairInfo.hpp"
#include "MemoryMappedObject.hpp"
#include "MultithreadedObject.hpp"
#include "BidirectionalReadGraph.hpp"
#include "ReadGraph.hpp"
#include "StringGraph.hpp"
#include "UnitigGraph.hpp"
#include "OrientedUnitigId.hpp"

#include "ReadId.hpp"
#include "AlignedEvidenceStore.hpp"
#include "OverlapCigarStore.hpp"
#include "dinaraTypes.hpp"
#include "MarkerKmers.hpp"
#include "mode3-Anchor.hpp"
#include "PafImport.hpp"

// Standard library.
#include "memory.hpp"
#include "string.hpp"
#include "utility.hpp"

// Boost.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include <atomic>
#include <mutex>
#include <unordered_map>

namespace dinara {

    class Assembler;
    class InvertedIndexFinder;
    class AssemblerInfo;
    class Alignment;
    class AlignmentData;
    class AlignmentGraph;
    class AlignmentInfo;
    class AlignOptions;
    class Align6Options;
    class Align6;
    class AssemblerOptions;
    class AssembledSegment;
    class AssemblyGraph2;
    class ClusterGraph;
    class CompressedAssemblyGraph;
    class ConsensusCaller;
    class Histogram2;
    class InducedAlignment;
    class KmerChecker;
    class KmerCounter;
    class KmersOptions;
    class LocalAssemblyGraph;
    class LocalAlignmentCandidateGraph;
    class LocalAlignmentGraph;
    class LocalMarkerGraph0;
    class LocalReadGraph;
    class LocalReadGraphTriangles;
    class LocalStringGraph;
    class LocalUnitigGraph;
    class LocalMarkerGraph0RequestParameters;
    class LongBaseSequences;
    class MarkerConnectivityGraph;
    class MarkerConnectivityGraphVertexMap;
    class MarkerKmers;
    class Mode2AssemblyOptions;
    class Mode3AssemblyOptions;
    class Mode3Assembler;
    class OrientedReadPair;
    class Reads;
    class ReferenceOverlapMap;
    class ProjectedAlignment;
    struct AnchorWindow;
    class Shasta2Anchors;
    class Shasta2Journeys;
    class Shasta2AnchorGraph;
    class Shasta2AssemblyGraph;
    class Shasta2AssemblyGraphOptions;
    class Shasta2AssemblyGraphPostprocessor;

    namespace mode0 {
        class AssemblyGraph;
    }

    namespace mode3 {
        class Anchors;
        class AnchorGraph;
        class BidirectedAnchors;
        class DirectedAnchorGraph;
    }

    namespace MemoryMapped {
        class ByteAllocator;
    }

    // Write an html form to select strand.
    void writeStrandSelection(
        ostream&,               // The html stream to write the form to.
        const string& name,     // The selection name.
        bool select0,           // Whether strand 0 is selected.
        bool select1);          // Whether strand 1 is selected.

    namespace Align4 {
        class MatrixEntry;
        class Options;
    }

    extern template class MultithreadedObject<Assembler>;
}

namespace spoa {
    class AlignmentEngine;
    class Graph;
}

class DisjointSets;


// Sanity check that we are compiling on x86_64 or aarch64
#if !__x86_64__ && !__aarch64__
#error "Dinara can only be built on an x86_64 machine (64-bit Intel/AMD) or an ARM64 machine. "
#endif



// Class used to store various pieces of assembler information in shared memory.
class dinara::AssemblerInfo {
public:

    // The read representation used: 0 = raw sequence, 1 = RLE sequence
    uint64_t readRepresentation;

    // The length of k-mers used to define markers.
    size_t k;

    // The method used to generate kmers (--Kmers.generationMethod).
    uint64_t kmerGenerationMethod;

    // The page size in use for this run.
    size_t largeDataPageSize;

    // Assembly mode (0=haploid, 2=phased).
    uint64_t assemblyMode;

    // Read graph creation method.
    uint64_t readGraphCreationMethod;

    // Variant Clustering options
    uint64_t variantClusteringMinOccurrences;
    uint64_t variantClusteringMinSeparation;


    // Statistics on the number of reads discarded on input.
    // These are incremented during each call to addReadsFromFasta.

    // The number of reads and raw bases discarded because the read
    // contained invalid bases.
    uint64_t discardedInvalidBaseReadCount = 0;
    uint64_t discardedInvalidBaseBaseCount = 0; // Only counts the valid bases in those reads.

    // The number of reads and raw bases discarded because the read length
    // was less than minReadLength.
    uint64_t discardedShortReadReadCount = 0;
    uint64_t discardedShortReadBaseCount = 0;

    // The number of reads and raw bases discarded because the read
    // contained repeat counts greater than 255.
    uint64_t discardedBadRepeatCountReadCount = 0;
    uint64_t discardedBadRepeatCountBaseCount = 0;



    // Statistics for the reads kept in the assembly
    // and not discarded on input.
    // These are computed and stored by histogramReadLength.
    size_t readCount = 0;
    size_t baseCount = 0;
    size_t readN50 = 0;
    uint64_t minReadLength = 0;

    // Other read statistics.
    size_t palindromicReadCount = 0;
    size_t chimericReadCount = 0;
    uint64_t isolatedReadCount = 0;
    uint64_t isolatedReadBaseCount = 0;

    // The coverage distribution of marker k-mers.
    // Only filled in for --Align.alignMethod 6.
    KmerDistributionInfo kmerDistributionInfo;

    // Alignment criteria actually used.
    // For readGraph creation method 0, they are the values specified
    // by the command line options and/or configuration.
    // For readGraph creation method 2, they are dynamically selected
    // based on alignments statistics.
    double actualMinAlignedFraction = 0;
    uint64_t actualMinAlignedMarkerCount = 0;
    uint64_t actualMaxDrift = 0;
    uint64_t actualMaxSkip = 0;
    uint64_t actualMaxTrim = 0;

    // Marker graph statistics.
    size_t markerGraphVerticesNotIsolatedCount = 0;
    size_t markerGraphEdgesNotRemovedCount = 0;
    uint64_t markerGraphMinCoverageUsed = 0;

    // Assembly graph statistics.
    size_t assemblyGraphAssembledEdgeCount = 0;
    size_t totalAssembledSegmentLength = 0;
    size_t longestAssembledSegmentLength = 0;
    size_t assembledSegmentN50 = 0;

    // Mode 2 assembly statistics.
    AssemblyGraph2Statistics assemblyGraph2Statistics;

    // Performance information.
    double assemblyElapsedTimeSeconds = 0.;
    double averageCpuUtilization;
    uint64_t peakMemoryUsage = 0ULL;
    uint64_t threadCount = 0;
    uint64_t virtualCpuCount = 0;
    uint64_t totalAvailableMemory = 0;

    inline string peakMemoryUsageForSummaryStats() {
        return peakMemoryUsage > 0 ? to_string(peakMemoryUsage) : "Not determined.";
    }
};

class dinara::Assembler :
    public MultithreadedObject<Assembler>,
    public MappedMemoryOwner,
    public HttpServer {
public:
    friend class InvertedIndexFinder;


    /***************************************************************************

    The constructors specify the file name prefix for binary data files.
    If this is a directory name, it must include the final "/".

    The constructor also specifies the page size for binary data files.
    Typically, for a large run binary data files will reside in a huge page
    file system backed by 2MB pages.
    1GB huge pages are also supported.
    The page sizes specified here must be equal to, or be an exact multiple of,
    the actual size of the pages backing the data.

    ***************************************************************************/

    // Constructor to be called one to create a new run.
    Assembler(
        const string& largeDataFileNamePrefix,
        bool createNew,
        uint64_t readRepresentation, // 0 = raw sequence, 1 = RLE sequence. Only used if createNew.
        size_t largeDataPageSize);

    // Add reads.
    // The reads in the specified file are added to those already previously present.
    void addReads(
        const string& fileName,
        uint64_t minReadLength,
        bool noCache,
        uint64_t threadCount);

    // Create a histogram of read lengths.
    void histogramReadLength(const string& fileName);


    // Functions related to markers.
    // See the beginning of Marker.hpp for more information.
    void findMarkers(uint64_t threadCount);
    void findMarkersSimdClosedSyncmers(uint64_t threadCount, int k, int s);
    // hifiasmFilter (a hifiasm_filter_t*, passed as void* to avoid the C bridge
    // header here) and hifiasmSampleDist enable hifiasm's overlap-path minimizer
    // filters on the hifiasm path; both are ignored when useHifiasm is false or
    // hifiasmFilter is null. The caller owns hifiasmFilter.
    void findMarkersSimdMinimizers(uint64_t threadCount, int k, int w,
        bool useHifiasm = true,
        const void* hifiasmFilter = nullptr, int hifiasmSampleDist = 0);
    void accessMarkers();
    void writeMarkers(ReadId, Strand, const string& fileName);

    // Write the reads that overlap a given read.
    void writeOverlappingReads(ReadId, Strand, const string& fileName);

    // Compute a marker alignment of two oriented reads.
    void alignOrientedReads(
        ReadId, Strand,
        ReadId, Strand,
        size_t maxSkip,  // Maximum ordinal skip allowed.
        size_t maxDrift, // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency
    );

    // Compute marker alignments of an oriented read with all reads
    // for which we have an Overlap.
    void alignOverlappingOrientedReads(
        ReadId, Strand,
        size_t maxSkip,   // Maximum ordinal skip allowed.
        size_t maxDrift,  // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency,
        size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
        size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
    );



    // Compute alignment and return the object (for EC parity)
    bool computeAlignmentParity(
        OrientedReadId id0,
        OrientedReadId id1,
        Alignment& outAlignment
    );

    // Compute an alignment for each alignment candidate.
    // Store summary information for the ones that are good enough,
    // without storing details of the alignment.
    void computeAlignments(
        const AlignOptions&,
        // Number of threads. If zero, a number of threads equal to
        // the number of virtual processors is used.
        uint64_t threadCount
    );

    // New unified alignment flow with evidence storage.
    void computeBaseAlignmentsAndStore(
        const AlignOptions&,
        uint64_t threadCount
    );

    // Lightweight marker-chain materialization for marker-graph prototypes.
    // Populates alignmentData, compressedAlignments, and alignmentTable from
    // precomputed chained marker alignments, but skips projected banded/base
    // alignment and sparse evidence generation.
    void computeAlignmentDataFromChainedCandidatesOnly(
        const AlignOptions&,
        uint64_t threadCount
    );

    class GlobalMismatchSiteClusters {
    public:
        // Unique nodes, each a (readId, position) pair in read forward coordinates.
        vector< pair<ReadId, uint32_t> > nodes;

        // Cluster representatives (indices into nodes).
        vector<uint64_t> clusterRepresentatives;

        // Cluster membership in CSR-like form:
        // members for cluster k are clusterMembers[clusterMemberOffsets[k] .. clusterMemberOffsets[k+1]).
        vector<uint64_t> clusterMemberOffsets;
        vector<uint64_t> clusterMembers;

        // Allele counts for each cluster, in Base value order (A,C,G,T => 0..3).
        vector< array<uint32_t, 4> > alleleCounts;

        // If alleleCounts have been updated using transitive (match-including) coverage.
        // Same length as alleleCounts.
        vector<uint8_t> alleleCountsAreTransitive;

        // If alleleCountsAreTransitive[clusterId]=1, this stores the number of reads/positions
        // that contributed to the transitive site count for that cluster.
        // Same length as alleleCounts.
        vector<uint32_t> transitiveSiteMemberCounts;
    };

    // Cluster mismatching (SNP) read positions into transitive "global sites".
    // This computes ProjectedAlignment for each stored overlap and unions the paired mismatch
    // positions (read0,pos0) <-> (read1,pos1) into connected components.
    GlobalMismatchSiteClusters clusterMismatchingPositionsIntoGlobalHetSites(
        const AlignOptions& alignOptions,
        uint64_t threadCount,
        bool includeDeletedAlignments = false,
        bool readGraphOnly = false
    ) const;

    // Like clusterMismatchingPositionsIntoGlobalHetSites, but only explores the connected component
    // reachable from the given seed read using alignmentTable, optionally with exploration limits.
    // Intended for fast debugging/printing (for example, sites involving read 0).
    GlobalMismatchSiteClusters clusterMismatchingPositionsIntoGlobalHetSitesReachableFromRead(
        ReadId seedReadId,
        const AlignOptions& alignOptions,
        uint64_t threadCount,
        uint64_t maxReadsToProcess = 0,
        uint64_t maxAlignmentsToProcess = 0,
        bool includeDeletedAlignments = false,
        bool readGraphOnly = false
    ) const;

    class TransitiveHetSiteCoverage {
    public:
        // Members as (readId, position) in forward coordinates. This can include reads
        // that match at the site (no mismatch token), discovered transitively using overlaps.
        vector< pair<ReadId, uint32_t> > members;

        // Allele counts over members, in Base value order (A,C,G,T => 0..3).
        array<uint32_t, 4> alleleCounts{0, 0, 0, 0};

        // Debug/performance counters.
        uint64_t alignmentsScanned = 0;
        uint64_t mappingHoles = 0;       // positions falling in deletions
        uint64_t mappingConflicts = 0;   // same read mapped to different positions
        bool hitNodeLimit = false;
        bool hitAlignmentLimit = false;
    };

    // Starting from a (readId, position) site seed, traverse overlaps transitively and
    // map the position across reads using only the sparse indel evidence (plus stored anchors).
    // This can be used to accumulate per-site allele counts across reads that do not directly
    // mismatch each other (simulating reference mapping in a de novo setting).
    TransitiveHetSiteCoverage gatherTransitiveHetSiteCoverage(
        ReadId seedReadId,
        uint32_t seedPosition,
        const AlignOptions& alignOptions,
        uint64_t maxNodesToVisit = 0,
        uint64_t maxAlignmentsToScan = 0,
        bool includeDeletedAlignments = false
    ) const;

    class GlobalHetSiteAlleleMembers {
    public:
        class Member {
        public:
            ReadId readId = invalidReadId;
            uint32_t position = 0;
        };

        // All members, grouped by (siteId, allele) where allele is in Base order A,C,G,T.
        vector<Member> members;

        // Offsets per site:
        // offsets[siteId] = {beginA, beginC, beginG, beginT, end}.
        vector< array<uint64_t, 5> > offsets;

        // Debug/performance counters.
        uint64_t propagatedAssignments = 0; // total unique (siteId,readId) assignments
        uint64_t mappingHoles = 0;          // positions falling in deletions
        uint64_t mappingConflicts = 0;      // same (siteId,readId) mapped to different positions
    };

    class GlobalHetSiteOrientedAlleleMembers {
    public:
        class Member {
        public:
            OrientedReadId orientedReadId;
            uint32_t position = 0; // Position on orientedReadId, 0-based.
        };

        // All members, grouped by (siteId, allele) where allele is in Base order A,C,G,T
        // in the oriented coordinate frame.
        vector<Member> members;

        // Offsets per site:
        // offsets[siteId] = {beginA, beginC, beginG, beginT, end}.
        vector< array<uint64_t, 5> > offsets;
    };

    // Compute full per-allele member lists for mismatch-defined global sites using only read-graph overlaps.
    // This starts from mismatch cluster members as seeds and propagates site positions transitively across
    // readGraph edges (alignmentData.info.isInReadGraph), mapping positions with sparse indels + marker anchors.
    // If seedReadId is specified, only sites containing that read are seeded and each such site is seeded from
    // a single position on that read (minimum observed position in the mismatch cluster).
    GlobalHetSiteAlleleMembers computeGlobalHetSiteAlleleMembersUsingReadGraph(
        const GlobalMismatchSiteClusters& clusters,
        const AlignOptions& alignOptions,
        uint64_t maxPendingTasks = 0,
        bool includeDeletedAlignments = false,
        ReadId seedReadId = invalidReadId
    ) const;

    class GlobalHetSitePositionVerificationStats {
    public:
        uint32_t siteId = 0;
        uint64_t expectedMembers = 0;
        uint64_t reachedMembers = 0;
        uint64_t checkedMappings = 0;
        uint64_t mismatchedPositions = 0;
        uint64_t mappingHoles = 0;
        uint64_t mappingFailures = 0;
        bool hitNodeLimit = false;
        bool hitAlignmentLimit = false;
    };

    // Debug/diagnostic: verify that the per-read positions assigned to a given siteId are reachable
    // and consistent via readGraph-only overlaps, starting from the mismatch seeds of that site.
    GlobalHetSitePositionVerificationStats debugVerifyGlobalHetSitePositionsUsingReadGraph(
        const GlobalMismatchSiteClusters& clusters,
        const GlobalHetSiteAlleleMembers& members,
        uint32_t siteId,
        const AlignOptions& alignOptions,
        uint64_t maxNodesToVisit = 0,
        uint64_t maxAlignmentsToScan = 0,
        bool includeDeletedAlignments = false
    ) const;

    class GlobalHetSiteReadIndex {
    public:
        class ReadSite {
        public:
            uint32_t siteId = 0;
            uint32_t readPosition = 0;
            uint8_t allele = 0; // Base value: A,C,G,T => 0..3
        };

        // Per-site propagated support counts in forward coordinates.
        // Indexed by siteId, in Base order A,C,G,T.
        vector< array<uint32_t, 4> > siteAlleleCounts;

        // Per-site pass flag after filtering.
        vector<uint8_t> sitePassesFilter;

        // Per-read, sorted by readPosition then siteId.
        // Contains only filtered, per-read-consistent site memberships.
        vector< vector<ReadSite> > sitesByRead;

        // Summary counters.
        uint64_t keptSiteCount = 0;
        uint64_t droppedAmbiguousReadSiteCount = 0;
        uint64_t droppedInvalidReadSiteCount = 0;
    };

    // Build a read-centric index of filtered global het sites from propagated members.
    // Site filter: keep sites with at least minAllelesWithMinSupport alleles having support
    // >= minAlleleSupport. Read-site filter: for a given (read,site), keep only if all
    // candidate entries agree on both position and allele; otherwise drop that read-site.
    GlobalHetSiteReadIndex buildFilteredGlobalHetSiteReadIndex(
        const GlobalHetSiteAlleleMembers& members,
        uint32_t minAlleleSupport = 3,
        uint32_t minAllelesWithMinSupport = 2
    ) const;

    // Return a strand assignment (0/1) for reads reachable from a seed read using readGraph-only overlaps.
    // For each kept overlap between reads r0 and r1, we enforce:
    //   strand[r1] = strand[r0] XOR (isSameStrand ? 0 : 1).
    // Unreachable reads get -1. Conflicts are counted when an already-assigned read is implied to
    // have the opposite strand via a different path (can happen in repetitive regions).
    vector<int8_t> computeReadGraphStrandsFromSeed(
        ReadId seedReadId,
        uint64_t& conflicts,
        bool includeDeletedAlignments = false
    ) const;

    // Convert forward-coordinate allele members to an oriented coordinate frame using a read->strand assignment.
    // For strand==1, positions are flipped as (len-1)-pos and alleles are complemented.
    GlobalHetSiteOrientedAlleleMembers orientGlobalHetSiteAlleleMembers(
        const GlobalHetSiteAlleleMembers& forwardMembers,
        const vector<int8_t>& strandByRead
    ) const;

    // Update per-cluster allele counts using transitive site coverage (includes matching reads).
    // This does not change cluster membership (clusters are still defined by mismatch connectivity),
    // it only updates clusters.alleleCounts based on the transitive mapper.
    void updateGlobalMismatchSiteClusterAlleleCountsWithTransitiveCoverage(
        GlobalMismatchSiteClusters& clusters,
        const AlignOptions& alignOptions,
        uint64_t maxNodesToVisitPerSite = 0,
        uint64_t maxAlignmentsToScanPerSite = 0,
        ReadId restrictToClustersInvolvingRead = invalidReadId,
        uint64_t maxClustersToUpdate = 0,
        bool includeDeletedAlignments = false
    ) const;

    // Old Phasing Logic Stub (for AssemblerPhasing.cpp compatibility)
    void performPhasing(uint64_t threadCount);
    void accessAlignmentData();
    void accessAlignmentDataReadWrite();


    // Loop over all alignments in the read graph
    // to create vertices of the global marker graph.
    // Throw away vertices with coverage (number of markers)
    // less than minCoverage or more than maxCoverage.
    // Also throw away "bad" vertices - that is, vertices
    // with more than one marker on the same oriented read.
    void createMarkerGraphVertices(

        // Minimum coverage (number of markers) for a vertex
        // of the marker graph to be kept.
        size_t minCoverage,

        // Maximum coverage (number of markers) for a vertex
        // of the marker graph to be kept.
        size_t maxCoverage,

        // Minimum coverage per strand (number of markers required
        // on each strand) for a vertex of the marker graph to be kept.
        uint64_t minCoveragePerStrand,

        // Flag that specifies whether to allow more than one marker on the
        // same oriented read id on a single marker graph vertex.
        bool allowDuplicateMarkers,

        // These two are used by PeakFinder in the automatic selection
        // of minCoverage when minCoverage is set to 0.
        double peakFinderMinAreaFraction,
        uint64_t peakFinderAreaStartIndex,

        // Number of threads. If zero, a number of threads equal to
        // the number of virtual processors is used.
        uint64_t threadCount
    );

    // BRG-aware variant of createMarkerGraphVertices.
    // Uses BidirectionalReadGraph edges (skipping isDeleted) for the
    // disjoint-set union step instead of ReadGraph edges.
    // The marker graph reflects the cleaned BRG overlap set.
    void createMarkerGraphVerticesFromBrg(
        size_t minCoverage,
        size_t maxCoverage,
        uint64_t minCoveragePerStrand,
        bool allowDuplicateMarkers,
        double peakFinderMinAreaFraction,
        uint64_t peakFinderAreaStartIndex,
        uint64_t threadCount
    );

    // Filter marker graph vertices whose marker k-mer is a short-period exact repeat
    // (including homopolymers). This removes vertices that tend to generate unreliable
    // anchors and artifacts in repetitive regions.
    // Must be called after createMarkerGraphVertices and before reverse-complement vertices/edges.
    void filterMarkerGraphVerticesByRepeatKmers(uint64_t threadCount);

    // Filter marker graph vertices whose marker k-mer has low sequence complexity,
    // assessed by counting distinct sub-k-mers of lengths 1, 2, 3, ...
    // Uses the same criterion as Shasta2's --min-anchor-distinct-subkmer-count option.
    // Must be called after createMarkerGraphVertices and before reverse-complement vertices/edges.
    void filterMarkerGraphVerticesByDistinctSubkmerCount(uint64_t threadCount);

    // Filter marker graph vertices where reads were grouped by transitive collapse
    // at k-mer positions outside their chaining range.
    // Must be called after createMarkerGraphVertices and computeCandidateTable.
    void filterMarkerGraphVerticesByChainConsistency(uint64_t threadCount);

    // Diagnostic (does not remove anything yet): detect the same class of
    // transitive-collapse false merges as filterMarkerGraphVerticesByChainConsistency
    // above, but via cross-read order agreement on consecutive journey vertices
    // instead of alignment-ordinal-range reconstruction. Must be called after
    // createMarkerGraphVertices.
    void filterMarkerGraphVerticesByJourneyOrderConsistency(uint64_t threadCount);

    // Create mode3 anchors from a subset of marker graph vertices selected by a sweep-line over
    // overlap start/end events on each oriented read (using read-graph overlaps).
    // This produces fewer anchors than using all marker graph vertices, while preserving
    // marker-graph semantics for each anchor.
    shared_ptr<mode3::Anchors> createAnchorsFromMarkerGraphVerticesAtOverlapEvents(
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t threadCount);

    // Create mode3 anchors from marker graph vertices selected as follows:
    // - Use a sweep-line over overlap start/end events for each oriented read (from readGraph overlaps).
    // - For each maximal interval where the active overlap count is >0, scan all marker ordinals in the interval
    //   and select the marker graph vertex (canonicalized by RC) with maximum vertex coverage in the requested range.
    // This produces significantly fewer, stronger anchors than using all marker graph vertices.
    shared_ptr<mode3::Anchors> createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval(
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t threadCount,
        bool enableColinearityPeeling = false,
        double minDominantFractionToPeel = 0.9);

    // Like createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval, but each selected marker graph vertex
    // is validated/split using the filtered readGraph overlaps among the oriented reads present in the vertex.
    // This avoids a single bridging read collapsing two unrelated regions into one anchor.
    shared_ptr<mode3::Anchors> createAnchorsFromMarkerGraphVerticesBestPerOverlapIntervalDecomposed(
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t threadCount);

    // Create mode3 anchors from *all* marker graph vertices, but split each marker graph vertex
    // into multiple anchors if the oriented reads inside the vertex fall into multiple clusters
    // based on surviving readGraph overlaps. This mitigates DSU transitive-collapse caused by
    // chimeric/bridging reads.
    shared_ptr<mode3::Anchors> createAnchorsFromMarkerGraphVerticesSplitUsingReadGraph(
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        const Mode3AssemblyOptions& mode3Options,
        uint64_t threadCount);

    // Create mode3 anchors directly from filtered overlaps without using markerGraph vertices.
    // For each oriented read, sweep overlap start/end events to find maximal intervals with active overlaps,
    // then select one marker ordinal per interval and gather matching marker ordinals on overlapping reads
    // using AlignmentInfo ordinal-offset bounds and k-mer validation.
    shared_ptr<mode3::Anchors> createAnchorsFromOverlapsBestPerOverlapInterval(
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t threadCount);

    // BidirectionalReadGraph-aware variant of createAnchorsFromOverlapsBestPerOverlapInterval.
    // Uses orientation-aware traversal via edge.traverse() instead of the strand-doubled
    // ReadGraph, preserving cross-strand overlaps at inversion/segdup boundaries.
    shared_ptr<mode3::Anchors> createAnchorsFromOverlapsBestPerOverlapIntervalBidirectional(
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t threadCount);

    // Create anchors from BRG-aware (self-RC) marker graph vertices.
    // Handles self-RC vertices by extracting only strand-0 markers for the
    // forward anchor and deriving the RC anchor by flipping.
    shared_ptr<mode3::Anchors> createAnchorsFromBrgMarkerGraphVertices(
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t threadCount);

    // Create BRG-native anchors (no RC doubling) from self-RC marker
    // graph vertices.  Returns a BidirectedAnchors object.
    shared_ptr<mode3::BidirectedAnchors> createBidirectedAnchors(
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t minEdgeCoverage,
        uint64_t threadCount);

    /// Save the bidirected anchor graph to GFA. Call after createBidirectedAnchors
    /// and after each modification (transitive reduction, unitigify, etc.) to persist
    /// the current state. Requires assembler.bidirectedAnchors to be set.
    void saveAnchorGraph(const std::string& fileName, bool includePaths = false) const;

    // Create and run Verkko-style directed anchor graph resolution.
    void runDirectedAnchorGraphResolution();


    // Find the vertex of the global marker graph that contains a given marker.
    // The marker is specified by the ReadId and Strand of the oriented read
    // it belongs to, plus the ordinal of the marker in the oriented read.
    MarkerGraphVertexId getGlobalMarkerGraphVertex(
        ReadId,
        Strand,
        uint32_t ordinal) const;

    // Find the markers contained in a given vertex of the global marker graph.
    // Returns the markers as tuples(read id, strand, ordinal).
    vector< tuple<ReadId, Strand, uint32_t> >
        getGlobalMarkerGraphVertexMarkers(MarkerGraph::VertexId) const;



    // Approximate transitive reduction of the marker graph.
    // This does the following, in this order:
    // - All edges with coverage less than or equal to lowCoverageThreshold
    //   are marked wasRemovedByTransitiveReduction.
    // - All edges with coverage 1 and a marker skip
    //   greater than edgeMarkerSkipThreshold
    //   are marked wasRemovedByTransitiveReduction.
    // - Edges with coverage greater than lowCoverageThreshold
    //   and less then highCoverageThreshold are processed in
    //   ordered of increasing coverage:
    //   * For each such edge A->B, we look for a path of length
    //     at most maxDistance starting at A and ending at B  that does not use
    //     edge A->B and also does not use any
    //     edges already marked wasRemovedByTransitiveReduction.
    //   * If such a path is found, the edge is marked
    //     wasRemovedByTransitiveReduction.
    // - Edges with coverage highCoverageThreshold or greater
    //   are left untouched.
    // The marker graph is guaranteed to be strand symmetric
    // when this begins, and we have to guarantee that it remains
    // strand symmetric when this ends.
    // To achieve this, we always process the two edges
    // in a reverse complemented pair together.
    void transitiveReduction(
        size_t lowCoverageThreshold,
        size_t highCoverageThreshold,
        size_t maxDistance,
        size_t edgeMarkerSkipThreshold);



    // Various pieces of assembler information stored in shared memory.
    // See class AssemblerInfo for more information.
public:
    MemoryMapped::Object<AssemblerInfo> assemblerInfo;
    uint64_t getMarkerGraphMinCoverageUsed() const
    {
        return assemblerInfo->markerGraphMinCoverageUsed;
    }
    
    // Filtering parameters for phased chains (passed from ReadGraphOptions)
    int minMultiNodeChainSupport = 6;
    int minIsolatedSiteSupport = 6;

private:



    // Reads in RLE representation.
    shared_ptr<Reads> reads;
public:
    const Reads& getReads() const {
        DINARA_ASSERT(reads);
        return *reads;
    }

    uint64_t adjustCoverageAndGetNewMinReadLength(uint64_t desiredCoverage);

    // Write a csv file with summary information for each read.
public:
    void writeReadsSummary();
    
    // ONT overlap phasing (hifiasm rphase_hc parity).
    void phaseOverlaps(uint64_t threadCount);

    // K-means overlap phasing (pgphase/longcallD-style iterative clustering).
    // Detects clean het sites from CIGAR walks, phases overlaps via k-means,
    // identifies noisy regions for future targeted MSA refinement.
    // isOnt: enables ONT-specific Fisher exact strand bias filter.
    void phaseOverlapsKmeans(uint64_t threadCount, bool isOnt = false, bool useEvidenceStore = false);

    // POA-based overlap phasing using Theseus MSA on anchor windows.
    // Requires shasta2Anchors and shasta2Journeys to be populated.
    // Coexists with phaseOverlapsKmeans — both write hifiasmEcMatchState.
    void phaseOverlapsMSA(uint64_t threadCount);

    // Hifiasm Error Correction
    void performHifiasmECParity(uint64_t threadCount);
    // Experimental: EC parity using induced alignments through marker graph vertices.
    // Requires marker graph vertices to be already created.
    // Replaces the SNP/SV detection pipeline with vertex-ordering consistency.
    void performHifiasmECParityWithMarkerGraph(uint64_t threadCount);
    // Debug: run het-site detection for one read and print all SNP/SV sites.
    // Call after computeBaseAlignmentsAndStore().
    void debugPrintHetSitesForRead(uint64_t readId);
    // Debug: dump all raw SNP/indel evidence for one read across its overlaps.
    // Call after computeBaseAlignmentsAndStore().
    void debugDumpAlignedEvidenceForRead(uint64_t readId);
    // Debug: aggregate SNP sites for one read; print only positions where both
    // ref and alt alleles have at least minSupport supporting alignments.
    void debugDumpSnpSitesForRead(uint64_t readId, uint32_t minSupport = 3);
    // Experimental global-site based phasing/EC pass.
    void performGlobalSiteECParity(uint64_t threadCount);


    void computeReadIdsSortedByName();

    // Find duplicate reads, as determined by name (not sequence).
    // This also sets the isDuplicate and discardDueToDuplicates read flags
    // and summarizes what it found Duplicates.csv.
    void findDuplicateReads(const string& handleDuplicates);


private:




    // The KmerChecker can find out if a given KmerId is a marker.
    shared_ptr<KmerChecker> kmerChecker;
    public:
    void createKmerChecker(
        const KmersOptions& kmersOptions,
        uint64_t threadCount);
    void setKmerChecker(shared_ptr<KmerChecker> inputKmerChecker)
    {
        kmerChecker = inputKmerChecker;
    }
    const shared_ptr<KmerChecker>& getKmerChecker() const
    {
        return kmerChecker;
    }
    void accessKmerChecker();

    // This one should eventually go away, but there are several scripts
    // that depend on it.
    void accessKmers()
    {
        accessKmerChecker();
    }


private:


    // Hash a KmerId in such a way that it has the same hash as its reverse
    // complement. This is used by alignment method 3 to downsample markers.
    uint32_t hashKmerId(KmerId) const;

    // The markers on all oriented reads. Indexed by OrientedReadId::getValue().
public:
    shared_ptr<MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>> markers;
private:
    void checkMarkersAreOpen() const;

    // Get markers sorted by KmerId for a given OrientedReadId.
    void getMarkersSortedByKmerId(
        OrientedReadId,
        vector<MarkerWithOrdinal>&) const;

    // Given a marker by its OrientedReadId and ordinal,
    // return the corresponding global marker id.
public:
    MarkerId getMarkerId(OrientedReadId, uint32_t ordinal) const;
private:
    MarkerId getReverseComplementMarkerId(OrientedReadId, uint32_t ordinal) const;
    MarkerId getMarkerId(const MarkerDescriptor& m) const
    {
        return getMarkerId(m.first, m.second);
    }
    MarkerId getReverseComplementMarkerId(const MarkerDescriptor& m) const
    {
        return getReverseComplementMarkerId(m.first, m.second);
    }

    // Inverse of the above: given a global marker id,
    // return its OrientedReadId and ordinal.
    // This requires a binary search in the markers toc.
    // This could be avoided, at the cost of storing
    // an additional 4 bytes per marker.
public:
    pair<OrientedReadId, uint32_t> findMarkerId(MarkerId) const;


    // KmerIds for all markers. Indexed by OrientedReadId::getValue().
    // Only stored during alignment computation, and then freed.
    shared_ptr<MemoryMapped::VectorOfVectors<KmerId, uint64_t>> markerKmerIds;
    void computeMarkerKmerIds(uint64_t threadCount);
    void cleanupMarkerKmerIds();
private:
    void computeMarkerKmerIdsThreadFunction(size_t threadId);


    // Pairs (KmerId, ordinal), sorted by KmerId, for each oriented read.
    // Indexed by orientedReadId.getValue().
    // Used by alignment method 4.
public:
    shared_ptr<MemoryMapped::VectorOfVectors< pair<KmerId, uint32_t>, uint64_t>> sortedMarkers;
    void computeSortedMarkers(uint64_t threadCount);
    bool accessSortedMarkers();
private:
    void computeSortedMarkersThreadFunction(size_t threadId);
    // void computeSortedMarkersThreadFunction1(size_t threadId);
    // void computeSortedMarkersThreadFunction2(size_t threadId);



    // Low frequency markers for each oriented read.
    // This stores, for each oriented read, the ordinals corresponding
    // to marker with low frequency (up to maxMarkerFrequency), sorted by KmerId.
    // Used by alignment method 5. It is only stored durign alignment
    // computation.
public:
    shared_ptr<MemoryMapped::VectorOfVectors<uint32_t, uint64_t>> lowFrequencyMarkers;
    void computeLowFrequencyMarkers(uint64_t maxMarkerFrequency, uint64_t threadCount);
    void computeLowFrequencyMarkers(
        const span<const KmerId>&,  // The marker k-mers for the oriented reads (sorted by ordinal)
        uint64_t maxMarkerFrequency,
        vector<uint32_t>&);         // The ordinals of the low frequency markers, sorted by KmerId.
private:
    void computeLowFrequencyMarkersThreadFunctionPass1(uint64_t threadId);
    void computeLowFrequencyMarkersThreadFunctionPass2(uint64_t threadId);
    void computeLowFrequencyMarkersThreadFunctionPass12(uint64_t pass);
    class ComputeLowFrequencyMarkersData {
    public:
        uint64_t maxMarkerFrequency;
    };
    ComputeLowFrequencyMarkersData computeLowFrequencyMarkersData;



    // Align6Markers, sorted by KmerId, for each oriented read.
    // Indexed by orientedReadId.getValue().
    // Used by alignment method 6.
public:
    shared_ptr<MemoryMapped::VectorOfVectors<Align6Marker, uint64_t>> align6Markers;
    void computeAlign6Markers(uint64_t threadCount);
    void accessAlign6Markers();
private:
    void computeAlign6MarkersThreadFunction(size_t threadId);



    // The alignment associated to each alignment candidate.
    // Computed by findAlignmentCandidatesInvertedIndex (Direct Chain Propagation).
    // If empty, computed by computeAlignments.
    class AlignmentCandidatesAlignmentsData {
    public:
        MemoryMapped::Vector<Alignment> alignments;
        MemoryMapped::Vector<int32_t> sharedSeedScores;
    };
public:
    AlignmentCandidatesAlignmentsData alignmentCandidatesAlignmentsData;
private:
    
    // Low level functions to get marker Kmers/KmerIds of an oriented read.
    // They are obtained from the reads and not from CompressedMarker::kmerId,
    // which will soon go away.

    // Get the marker Kmer for an oriented read and ordinal.
    Kmer getOrientedReadMarkerKmer(OrientedReadId, uint32_t ordinal) const;
    Kmer getOrientedReadMarkerKmerStrand0(ReadId, uint32_t ordinal) const;
    Kmer getOrientedReadMarkerKmerStrand1(ReadId, uint32_t ordinal) const;

    // Get the marker KmerId for an oriented read and ordinal.
    KmerId getOrientedReadMarkerKmerId(OrientedReadId, uint32_t ordinal) const;

    // Get all marker Kmers for an oriented read.
    void getOrientedReadMarkerKmers(OrientedReadId, const span<Kmer>&) const;
    void getOrientedReadMarkerKmersStrand0(ReadId, const span<Kmer>&) const;
    void getOrientedReadMarkerKmersStrand1(ReadId, const span<Kmer>&) const;

    // Get all marker KmerIds for an oriented read.
    void getOrientedReadMarkerKmerIds(OrientedReadId, const span<KmerId>&) const;
    void getOrientedReadMarkerKmerIdsStrand0(ReadId, const span<KmerId>&) const;
    void getOrientedReadMarkerKmerIdsStrand1(ReadId, const span<KmerId>&) const;

    // Get all MarkerWithOrdinals for an oriented read (includes position, KmerId, and ordinal).
    void getOrientedReadMarkers(OrientedReadId, const span<MarkerWithOrdinal>&) const;
    void getOrientedReadMarkersStrand0(ReadId, const span<MarkerWithOrdinal>&) const;
    void getOrientedReadMarkersStrand1(ReadId, const span<MarkerWithOrdinal>&) const;

    // Get all Align6::Markers for an oriented read (includes KmerId, and ordinal, and global frequency).
    // Sorted by KmerId.
    void getOrientedReadAlign6Markers(OrientedReadId, const span<Align6Marker>&) const;
    void getOrientedReadAlign6MarkersStrand0(ReadId, const span<Align6Marker>&) const;
    void getOrientedReadAlign6MarkersStrand1(ReadId, const span<Align6Marker>&) const;

    // Get all marker Kmers/KmerIds for a read in both orientations.
    void getReadMarkerKmers(
        ReadId,
        const span<Kmer>& Kmers0,
        const span<Kmer>& Kmers1) const;
    void getReadMarkerKmerIds(
        ReadId,
        const span<KmerId>& kmerIds0,
        const span<KmerId>& kmerIds1) const;

    // Get the Kmer/KmerId for an oriented read at a given marker ordinal.
    Kmer getOrientedReadMarkerKmer(OrientedReadId, uint64_t ordinal) const;
    KmerId getOrientedReadMarkerKmerId(OrientedReadId, uint64_t ordinal) const;



    // Given a MarkerId, compute the MarkerId of the
    // reverse complemented marker.
    MarkerId findReverseComplement(MarkerId) const;

    // Counting of marker Kmers.
public:
    shared_ptr<KmerCounter> kmerCounter;
    // Count k-mers in markers.
    // If the k-mer counts were already computed, this recomputes them.
    // The previous k-mer counts are lost.
    void countKmers(
        uint64_t threadCount,
        const string& globalFrequencyOverrideDirectory
    );
    
    // Optimized variant using markerKmerIds.
    void countKmersFromMarkerKmerIds(uint64_t threadCount);

    // This creates various histograms of k-mer frequencies:
    // 1. One line per k-mer. Frequencies of this k-mer in reads and in the marker graph.
    void accessKmerCounts();



    // The MarkerKmers keep track of the locations in the oriented reads
    // where each marker k-mer appears. It is only used for alignment-free assembly
    // (--Assembly.mode 3 --Assembly.mode3.anchorCreationMethod FromMarkerKmers)
    // but can also be created using CreateMarkerKmers.py.
    shared_ptr<MarkerKmers> markerKmers;
    void createMarkerKmers(uint64_t threadCount);
    void accessMarkerKmers();



    // Flag palindromic reads.
    // Uses inverted index + DP chaining (same pipeline as overlap discovery)
    // to find self-alignments (strand 0 vs strand 1), then ProjectedAlignment
    // with astarpa for base-level identity.
    // Requires buildInvertedIndex to have been called first.
    void flagPalindromicReads(
        double maxDriftRate,
        const OverlapCandidatesOptions& overlapCandidatesOptions,
        double alignedFractionThreshold,
        double maxErrorRate,
        uint64_t threadCount);

	    // Filter secondary/redundant alignments per read pair (hifiasm-style).
		    void filterSecondaryAlignmentsPerReadPairThreadFunction(size_t threadId);
		    void keepOnlyBestAlignmentPerReadPairByDpScoreThreadFunction(size_t threadId);
		    void deduplicateOntChainsPerPartnerReadHifiasmLikeThreadFunction(size_t threadId);
		    std::atomic<uint64_t> removedSecondaryAlignmentCount;
		    std::atomic<uint64_t> removedSecondaryAlignmentBySymmetryOnlyCount;
		    std::atomic<uint64_t> removedOntDeduplicatedChainCount;
		    bool filterSecondaryRequireNonRedundantOnBothReads = false;
		public:
		    void filterSecondaryAlignmentsPerReadPair(
		        uint64_t threadCount,
		        bool requireNonRedundantOnBothReads = false);
	public:
	    // ONT EC parity helper:
	    // In hifiasm's ONT pipeline, `h_ec_lchain` can emit multiple chains per (query, partner).
	    // Later, after base-level alignment refinement and phasing (`gen_hc_r_alin_ea`, `rphase_hc`),
	    // hifiasm collapses duplicates with `dedup_chains` (ecovlp.cpp:2984), using base-level error
	    // statistics and cis/trans state to select one best overlap per partner.
	    //
	    // Dinara's inverted-index lchain+mcopy path mirrors the "emit multiple chains" behavior.
	    // This function mirrors the "dedup later" stage and marks all but the best overlap per
	    // (readIds[0], readIds[1]) as `DeleteReasonSecondary`.
	    void deduplicateOntChainsPerPartnerReadHifiasmLike(uint64_t threadCount);
	public:
	    // More aggressive than filterSecondaryAlignmentsPerReadPair:
	    // for each (read0, read1) pair, keep only the single best overlap by DP score
	    // among overlaps that are cis on both reads and not already deleted.
	    void keepOnlyBestAlignmentPerReadPairByDpScore(uint64_t threadCount);
public:
    // Pre-phasing chain deduplication: keep one best chain per read pair.
    // Port of hifiasm's dedup_chains (ecovlp.cpp:2984) for use before
    // marker graph construction, when phasing info is not yet available.
    // Uses score = span - 12*errors, then span as tiebreaker.
    void dedupChainsPrePhasing(uint64_t threadCount);
private:
    void dedupChainsPrePhasingThreadFunction(size_t threadId);
    std::atomic<uint64_t> removedPrePhasingDedupCount;
public:

    // Remove all chains for read pairs that have multiple chains on the
    // same strand. Such multi-chain pairs are likely repeat-induced.
    void removeMultiChainAlignments(uint64_t threadCount);
private:
    void removeMultiChainAlignmentsThreadFunction(size_t threadId);
    std::atomic<uint64_t> removedMultiChainCount;
public:
    // Hifiasm-style filtering methods (called from main.cpp)
    void filterLocalSegments(uint64_t minCoverage, uint64_t threadCount);
    void applyCoverageCuts(uint64_t minOverlapLength, uint64_t threadCount);
    void filterHangingOverlaps(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount);
    /// Delete overlaps where one read is contained in the other (ma_hit2arc containment).
    /// Runs early (e.g. after computeBaseAlignmentsAndStore) to remove containment overlaps.
    void deleteContainmentOverlaps(uint64_t threadCount);
    /// Delete internal overlaps (ma_hit2arc MA_HT_INT/MA_HT_SHORT_OVLP: excessive overhangs or too short).
    /// Runs early (e.g. after deleteContainmentOverlaps) to remove spurious internal matches.
    /// Uses stored CIGAR boundary coordinates directly.
    void deleteInternalOverlaps(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount);
    /// Same as deleteInternalOverlaps but extends CIGAR boundary coordinates
    /// toward read tips before classification (matching hifiasm's
    /// append_inexact_overlap_region_alloc). Less aggressive than the
    /// non-extended version.
    void deleteInternalOverlapsExtended(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount);
    void filterOverlapsByRegionalCliques(uint64_t minIntervalOverlap, uint64_t minRegionSize, double minCliqueFraction, uint64_t threadCount);
    void removeContainedReads(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount);
    void removeReadsFlaggedContained(uint64_t threadCount);
    void flagContainedReads(uint64_t maxHang, double maxHangRate, uint64_t minOverlapLength, uint64_t threadCount);
    void pruneContainedReadsToOneBestOverlapByDpScore(uint64_t threadCount);

    // Rescue overlaps where one side says trans but the partner says cis.
    // Adapted from hifiasm's try_rescue_overlaps. If ≥ minPileup disagreement
    // overlaps pile up spatially on a read, that read's trans calls are flipped
    // to cis. Call after phaseOverlapsKmeans, before createReadGraphFromPhasingCisOverlaps.
    // skipDeleted: when true, skip overlaps with any deleteReason set by earlier
    // pipeline stages (chimeric, contained, weak, etc.). Maps to hifiasm's
    // is_del parameter, which is true for ONT and false for HiFi.
    void rescueTransOverlaps(uint64_t minPileup = 4, bool skipDeleted = false);

    void applyOntChemicalArcMask(uint64_t threadCount);
    void applyOntChemicalArcMask(uint64_t chemicalCov, uint64_t chemicalFlank, double dupRate, uint64_t threadCount);

private:
    // Filter local segments (coverage based) - thread function
    void filterLocalSegmentsThreadFunction(size_t threadId);
    
    // Structure for local segment filtering results
    struct ReadSegmentStatus {
        uint32_t start = 0;
        uint32_t end = 0;
        bool isDeleted = false;
    };

    // Store valid intervals and status for each read.
    std::vector<ReadSegmentStatus> validReadIntervals;

    // Temporary storage for filtering parameter, accessible by thread function.
    uint64_t localSegmentMinCoverage = 0;

    // Apply coverage cuts (ma_hit_cut equivalent) - thread functions
    void applyCoverageCutsToAlignmentsThreadFunction(size_t threadId);
    void applyCoverageCutsCleanupThreadFunction(size_t threadId);
    uint64_t coverageCutMinOverlap = 0;

    // Filter hanging overlaps (ma_hit_flt equivalent) - thread function
    void filterHangingOverlapsThreadFunction(size_t threadId);
    // Delete containment overlaps - thread function
    void deleteContainmentOverlapsThreadFunction(size_t threadId);
    // Delete internal overlaps - thread functions
    void deleteInternalOverlapsThreadFunction(size_t threadId);
    void deleteInternalOverlapsExtendedThreadFunction(size_t threadId);
    uint64_t hangingFilterMaxHang = 1000;
    double hangingFilterMaxHangRate = 0.8;
    uint64_t hangingFilterMinOverlap = 0;

    // ONT chemical arc masking (hifiasm gen_chemical_arc_rf equivalent).
    void applyOntChemicalArcMaskThreadFunction(size_t threadId);
    uint64_t chemicalArcCov = 1;
    uint64_t chemicalArcFlank = 256;
    double chemicalArcDupRate = 0.02;
    std::vector<uint8_t> chemicalArcMask;


    // Mapping from contained read to its container (parent).
    // Initialized to invalidReadId.
    shared_ptr<MemoryMapped::Vector<ReadId>> containmentParent;


    // Check if an alignment between two reads should be suppressed,
    // bases on the setting of command line option
    // --Align.sameChannelReadAlignment.suppressDeltaThreshold.
    bool suppressAlignment(ReadId, ReadId, uint64_t delta);

    // Remove all alignment candidates for which suppressAlignment
    // returns false.
public:
    void suppressAlignmentCandidates(uint64_t delta, uint64_t threadCount);
private:
    class SuppressAlignmentCandidatesData {
    public:
        uint64_t delta;
        MemoryMapped::Vector<bool> suppress; // For each alignment candidate.
    };
    SuppressAlignmentCandidatesData suppressAlignmentCandidatesData;
    void suppressAlignmentCandidatesThreadFunction(size_t threadId);



    // Alignment candidates found by the LowHash algorithm.
    // They all have readId0<readId1.
    // They are interpreted with readId0 on strand 0.
public:
    AlignmentCandidates alignmentCandidates;
private:

    // PAF-imported overlap intervals, keyed by canonical read pair (readId0<readId1).
    // Populated by importAlignmentCandidatesFromPaf and consumed by chainPafCandidates
    // to constrain shared-minimizer collection to the interval hifiasm already agreed on.
    // Each pair can hold both a same-strand and a reverse overlap (PafPairIntervals),
    // since inverted repeats can produce both. PafCandidateInterval/PafPairIntervals are
    // defined in PafImport.hpp; coordinates are half-open base positions on each read,
    // target coordinates forward-strand (Alignment::ts/te).
    // Key packs (readId0<<32 | readId1) with readId0<readId1.
public:
    std::unordered_map<uint64_t, PafPairIntervals> pafCandidateIntervals;
private:

public:
    void writeAlignmentCandidates(bool useReadName=false, bool verbose=false) const;
private:


    // Use the LowHash (modified MinHash) algorithm to find candidate alignments.
    // Use as features sequences of m consecutive special k-mers.
public:
    void findAlignmentCandidatesLowHash0(
        size_t m,                       // Number of consecutive k-mers that define a feature.
        double hashFraction,            // Low hash threshold.
        // Iteration control. See MinHashOptions for details.
        size_t minHashIterationCount,
        double alignmentCandidatesPerRead,

        size_t log2MinHashBucketCount,  // Base 2 log of number of buckets for lowHash.
        size_t minBucketSize,           // The minimum size for a bucket to be used.
        size_t maxBucketSize,           // The maximum size for a bucket to be used.
        size_t minFrequency,            // Minimum number of lowHash hits for a pair to become a candidate.
        uint64_t threadCount
    );
    void markAlignmentCandidatesAllPairs();
    void accessAlignmentCandidates();
    void accessAlignmentCandidateTable();
    vector<OrientedReadPair> getAlignmentCandidates() const;
    void computeCandidateTable();
    // threadCount == 0 means "use all hardware threads".
    void importAlignmentCandidatesFromPaf(const string& pafFilePath, uint64_t threadCount = 0);
    
    // Chain pre-imported PAF candidates using the inverted index.
    // buildInvertedIndex must be called before this.
    void chainPafCandidates(
        double maxDriftRate,
        uint64_t maxChainLimit,
        const OverlapCandidatesOptions& overlapCandidatesOptions,
        uint64_t threadCount
    );

private:
    void checkAlignmentCandidatesAreOpen() const;


    // LowHash statistics for read.
    // For each read we count the number of times a low hash
    // puts the read in:
    // - A sparse bucket (bucketSize < minFrequency), index 0.
    // - A good bucket (minFrequency <= bucketSize <= maxFrequency), index 1.
    // - A crowded bucket (bucketSize > maxFrequency), index 2.
    MemoryMapped::Vector< array<uint64_t, 3> > readLowHashStatistics;
    void accessReadLowHashStatistics();

    bool createLocalAlignmentCandidateGraph(
        vector<OrientedReadId>& starts,
        uint32_t maxDistance,               // How far to go from starting oriented read.
        bool allowChimericReads,
        double timeout,                     // Or 0 for no timeout.
        bool inGoodAlignmentsRequired,      // Only add an edge to the local graph if it's in the "good" alignments
        bool inReadgraphRequired,           // Only add an edge to the local graph if it's in the ReadGraph
        LocalAlignmentCandidateGraph& graph
    );

    // This method is used as an alternative to createLocalAlignmentCandidateGraph, in the case that the user
    // wants to see only the edges that are inferred from the PAF, and none others. Coloring/labelling w.r.t.
    // the different subgroups still applies (candidate, good alignment, read graph)
    bool createLocalReferenceGraph(
        vector<OrientedReadId>& starts,
        uint32_t maxDistance,           // How far to go from starting oriented read.
        bool allowChimericReads,
        double timeout,                 // Or 0 for no timeout.
        LocalAlignmentCandidateGraph& graph
    );
public:
    // Construct a subgraph around a read in the candidate graph
    // and write all subgraph reads to a Fasta file.
    void writeLocalAlignmentCandidateReads(
            ReadId readId,
            Strand strand,
            uint32_t maxDistance,
            bool allowChimericReads,
            bool allowCrossStrandEdges,
            bool allowInconsistentAlignmentEdges
    );
private:
    // Compute a marker alignment of two oriented reads.
    void alignOrientedReads(
        OrientedReadId,
        OrientedReadId,
        size_t maxSkip,     // Maximum ordinal skip allowed.
        size_t maxDrift,    // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency
    );
    // This lower level version takes as input vectors of
    // markers already sorted by kmerId.
    void alignOrientedReads(
        const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
        size_t maxSkip,  // Maximum ordinal skip allowed.
        size_t maxDrift, // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency
    );
    // This version allows reusing the AlignmentGraph and Alignment
    void alignOrientedReads(
        const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
        size_t maxSkip,             // Maximum ordinal skip allowed.
        size_t maxDrift,            // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency,
        bool debug,
        AlignmentGraph&,
        Alignment&,
        AlignmentInfo&
    );
public:
    void analyzeAlignmentMatrix(ReadId, Strand, ReadId, Strand);
private:


    // Alternative alignment functions with 1 suffix (SeqAn).
    void alignOrientedReads1(
        OrientedReadId,
        OrientedReadId,
        int matchScore,
        int mismatchScore,
        int gapScore);
    void alignOrientedReads1(
        OrientedReadId,
        OrientedReadId,
        int matchScore,
        int mismatchScore,
        int gapScore,
        Alignment&,
        AlignmentInfo&);
public:
    void alignOrientedReads1(
        ReadId, Strand,
        ReadId, Strand,
        int matchScore,
        int mismatchScore,
        int gapScore);
private:


    // Alternative alignment function with 3 suffix (SeqAn, banded).
    void alignOrientedReads3(
        OrientedReadId,
        OrientedReadId,
        int matchScore,
        int mismatchScore,
        int gapScore,
        double downsamplingFactor,
        int bandExtend,
        int maxBand,
        Alignment&,
        AlignmentInfo&);



    // Chimeric Read Detection
public:
    void detectChimericReads(uint64_t threadCount);
    void rescueChimericReads(uint64_t threadCount);
    
    // Rescue phased overlaps (try_rescue_overlaps equivalent)
    // Rescues overlaps with directional phasing conflicts if sufficient neighbor support exists.
    void rescuePhasedOverlaps(uint64_t rescueThreshold, uint64_t threadCount);

    // Count active alignments (kept by both sides)
    uint64_t countActiveAlignments() const;

private:
    void detectChimericReadsThreadFunction(size_t threadId);
    void rescueChimericReadsThreadFunction(size_t threadId);
    void rescuePhasedOverlapsThreadFunction(size_t threadId);
    uint64_t rescuePhasedThreshold = 4; // Default threshold from Hifiasm

    // Data for chimeric detection (if needed across threads)
    class ChimericDetectionData {
    public:
         // Thread-local or shared data
    };
    ChimericDetectionData chimericDetectionData;
    
    // Store reads identified as chimeric
    MemoryMapped::Vector<bool> isChimericRead;
    // Thread-safe temporary buffer used when detecting chimeric reads in parallel.
    // We avoid writing to MemoryMapped::Vector<bool> concurrently because it is bit-packed.
    std::vector<uint8_t> chimericReadTmp;

public:

    // Python-callable version.
    void alignOrientedReads4(
        ReadId, Strand,
        ReadId, Strand,
        uint64_t deltaX,
        uint64_t deltaY,
        uint64_t minEntryCountPerCell,
        uint64_t maxDistanceFromBoundary,
        uint64_t minAlignedMarkerCount,
        double minAlignedFraction,
        uint64_t maxSkip,
        uint64_t maxDrift,
        uint64_t maxTrim,
        uint64_t maxBand,
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore);

    // Align two reads using alignment method 4.
    // If debug is true, detailed output to html is produced.
    // Otherwise, html is not used.
    void alignOrientedReads4(
        OrientedReadId,
        OrientedReadId,
        const Align4::Options&,
        MemoryMapped::ByteAllocator&,
        Alignment&,
        AlignmentInfo&,
        bool debug);

    // Intermediate level version used by the http server.
    void alignOrientedReads4(
        OrientedReadId,
        OrientedReadId,
        uint64_t deltaX,
        uint64_t deltaY,
        uint64_t minEntryCountPerCell,
        uint64_t maxDistanceFromBoundary,
        uint64_t minAlignedMarkerCount,
        double minAlignedFraction,
        uint64_t maxSkip,
        uint64_t maxDrift,
        uint64_t maxTrim,
        uint64_t maxBand,
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore,
        Alignment&,
        AlignmentInfo&
        );

    // Alignment method 5.
    void alignOrientedReads5(
        OrientedReadId,
        OrientedReadId,
        int matchScore,
        int mismatchScore,
        int gapScore,
        double driftRateTolerance,
        uint64_t minBandExtend,
        Alignment&,
        AlignmentInfo&,
        ostream& html);

    // Alignment method 6.
    void alignOrientedReads6(
        OrientedReadId,
        OrientedReadId,
        Alignment&,
        AlignmentInfo&,
        Align6&);

private:


    // Create a local alignment graph starting from a given oriented read
    // and walking out a given distance on the global alignment graph.
    // An alignment graph is an undirected graph in which each vertex
    // represents an oriented read. Two vertices are joined by an
    // undirected edge if we have a found a good alignment,
    // and stored it, between the corresponding oriented reads.
    bool createLocalAlignmentGraph(
        OrientedReadId,
        size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
        size_t maxTrim,                 // Maximum left/right trim (expressed in bases).
        uint32_t distance,              // How far to go from starting oriented read.
        double timeout,                 // Or 0 for no timeout.
        LocalAlignmentGraph&
    );

    // Compute marker alignments of an oriented read with all reads
    // for which we have an Overlap.
    void alignOverlappingOrientedReads(
        OrientedReadId,
        size_t maxSkip,                 // Maximum ordinal skip allowed.
        size_t maxDrift,                // Maximum ordinal drift allowed.
        uint32_t maxMarkerFrequency,
        size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
        size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
    );



    // Count the common marker near a given ordinal offset for
    // two oriented reads. This can be used to check
    // whether an alignmnent near the specified ordinal offset exists.
public:
    uint32_t countCommonMarkersNearOffset(
        OrientedReadId,
        OrientedReadId,
        int32_t offset,
        int32_t offsetTolerance
    );
    uint32_t countCommonMarkersWithOffsetIn(
        OrientedReadId,
        OrientedReadId,
        int32_t minOffset,
        int32_t maxOffset
    );
private:



    // The good alignments we found.
    // They are stored with readId0<readId1 and with strand0==0.
    // The order in compressedAlignments matches that in alignmentData.
public:
    MemoryMapped::Vector<AlignmentData> alignmentData;
    const MemoryMapped::VectorOfVectors<uint32_t, uint32_t>& getAlignmentTable() const
    {
        return alignmentTable;
    }
    const OverlapCigarStore& getOverlapCigarStore() const
    {
        return overlapCigarStore;
    }
#ifdef DINARA_TESTING
public:
    // Test-only hook: allows integration tests to construct a consistent alignmentTable
    // from manually-populated alignmentData.
    void computeAlignmentTableForTesting()
    {
        computeAlignmentTable();
    }

    uint64_t getRemovedSecondaryAlignmentBySymmetryOnlyCountForTesting() const
    {
        return removedSecondaryAlignmentBySymmetryOnlyCount.load();
    }

    struct ValidReadIntervalForTesting {
        uint32_t start = 0;
        uint32_t end = 0;
        bool isDeleted = false;
    };

    ValidReadIntervalForTesting getValidReadIntervalForTesting(ReadId readId) const
    {
        if (readId >= validReadIntervals.size()) {
            return {};
        }
        const auto& v = validReadIntervals[readId];
        return {v.start, v.end, v.isDeleted};
    }

    ReadId getContainmentRootForTesting(ReadId readId) const
    {
        if (!containmentParent || !containmentParent->isOpen) {
            return ReadId(invalidReadId);
        }
        if (readId >= containmentParent->size()) {
            return ReadId(invalidReadId);
        }
        return (*containmentParent)[readId];
    }
#endif
private:
    MemoryMapped::VectorOfVectors<char, uint64_t> compressedAlignments;

    void checkAlignmentDataAreOpen() const;
public:
    void accessCompressedAlignments();
private:



    // Get the stored alignments involving a given oriented read.
    // This performs swaps and reverse complementing as necessary,
    // to return alignments in which the first oriented read is
    // the one specified as the argument.
    class StoredAlignmentInformation {
    public:
        uint64_t alignmentId;
        OrientedReadId orientedReadId;
        Alignment alignment;
    };
    void getStoredAlignments(
        OrientedReadId,
        vector<StoredAlignmentInformation>&) const;
    void getStoredAlignments(
        OrientedReadId,
        const vector<OrientedReadId>&,
        vector<StoredAlignmentInformation>&) const;



    // The alignment table stores the AlignmentData that each oriented read is involved in.
    // Stores, for each OrientedReadId, a vector of indexes into the alignmentData vector.
    // Indexed by OrientedReadId::getValue(),
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> alignmentTable;
    void computeAlignmentTable();





    // Private functions and data used by computeAlignments.
    void computeAlignmentsThreadFunction(size_t threadId);
    void computeBaseAlignmentsAndStoreThreadFunction(size_t threadId);
    class ComputeAlignmentsData {
    public:

        // Not owned.
        const AlignOptions* alignOptions = 0;

        // The AlignmentInfo found by each thread.
        vector< vector<AlignmentData> > threadAlignmentData;

        // Compressed alignments corresponding to the AlignmentInfo found by each thread.
        vector< shared_ptr< MemoryMapped::VectorOfVectors<char, uint64_t> > > threadCompressedAlignments;

        // Variant clustering position pairs collected by each thread
        // Each pair is (OrientedReadId, position)
        vector< shared_ptr< MemoryMapped::Vector< pair<OrientedReadId, uint32_t> > > > threadVariantClusteringPositionPairs;
        
        // Timing accumulators for variant clustering (per thread, in seconds)
        vector<double> threadProjectedAlignmentTime;  // Time spent in ProjectedAlignment construction
        vector<double> threadCollectionTime;          // Time spent in collectVariantClusteringPositionPairs
        vector<uint64_t> threadFilteredByErrorRate; 
        vector<uint64_t> threadFilteredByErrorRateGap; 
        vector<uint64_t> threadFilteredByGapCount;

        // Thread-local packed CIGAR stores (hifiasm-style uint16_t tokens).
        // One store per thread to avoid locking.
        vector<OverlapCigarStore> threadCigarStores;

        // Thread-local Evidence Stores (APES/TASSD)
        // One store per thread to avoid locking.
        vector<AlignedEvidenceStore> threadEvidenceStores;
    };
    ComputeAlignmentsData computeAlignmentsData;



    // Find in the alignment table the alignments involving
    // a given oriented read, and return them with the correct
    // orientation (this may involve a swap and/or reverse complement
    // of the AlignmentInfo stored in the alignmentTable).
    vector< pair<OrientedReadId, AlignmentInfo> >
        findOrientedAlignments(OrientedReadId, bool inReadGraphOnly,
            vector<uint32_t>* alignmentIds = nullptr) const;

public:
    // Like findOrientedAlignments, but returns alignment ids sorted by the number of informative
    // het/SV sites covered by the overlap (from the query read's perspective), descending.
    // This does not change the ordering in alignmentTable (which is sorted by partner OrientedReadId).
    void getAlignmentIdsSortedByInformativeSites(
        OrientedReadId,
        vector<uint32_t>& alignmentIds,
        bool inReadGraphOnly = false) const;

    // Return alignment ids involving this oriented read that are cis in both views
    // (DeleteReasonPhase not set on either side) and sorted by informativeHetSiteScore, descending.
    // If keptByBothSidesOnly is true, also require no deletion reasons on either side.
    void getCisAlignmentIdsSortedByInformativeSites(
        OrientedReadId,
        vector<uint32_t>& alignmentIds,
        bool keptByBothSidesOnly = false) const;

    // Like getCisAlignmentIdsSortedByInformativeSites, but returns all alignments
    // (not limited to a specific read).
    void getAllCisAlignmentIdsSortedByInformativeSites(
        vector<uint32_t>& alignmentIds,
        bool keptByBothSidesOnly = false) const;



    // Analyze the stored alignments involving a given oriented read.
private:
    void analyzeAlignments1(ReadId, Strand) const;


    // Read graph and related functions and data.
    // For more information, see comments in ReadGraph.hpp.
	public:
	    ReadGraph readGraph;
	    ReadGraph readGraphAllAlignments;
	    DirectedReadGraph directedReadGraph;
	    BidirectionalReadGraph bidirectionalReadGraph;
	    StringGraph stringGraph;
	    UnitigGraph unitigGraph;
    void createReadGraph(
        uint32_t maxAlignmentCount,
        bool preferAlignedFraction);

    void createReadGraph2(
        uint32_t maxAlignmentCount,
        double markerCountPercentile,
        double alignedFractionPercentile,
        double maxSkipPercentile,
        double maxDriftPercentile,
        double maxTrimPercentile);

    void setReadGraph2Criteria(
            double markerCountPercentile,
            double alignedFractionPercentile,
            double maxSkipPercentile,
            double maxDriftPercentile,
            double maxTrimPercentile);

    bool passesReadGraph2Criteria(const AlignmentInfo& info) const;
    void accessReadGraph();
    void accessReadGraphReadWrite();
    void checkReadGraphIsOpen() const;
    void accessDirectedReadGraph();
    void checkDirectedReadGraphIsOpen() const;
    // BidirectionalReadGraph: one vertex per read, one edge per alignment.
    void createBidirectionalReadGraph();
    void createBidirectionalReadGraphFromSelectedAlignments(const vector<bool>& keepAlignment);
    void accessBidirectionalReadGraph();
    void accessBidirectionalReadGraphReadWrite();
    void checkBidirectionalReadGraphIsOpen() const;
    void removeBidirectionalReadGraph();

    // BidirectionalReadGraph cleaning (string-graph-style operations on BRG).
    // Individual operations (each builds a temporary directed view):
    uint64_t reduceBidirectionalReadGraphTransitive(uint32_t gapFuzz = 1000);
    uint64_t cutBidirectionalReadGraphTips(uint32_t maxShortTipReads = 3);
    uint64_t cutBidirectionalReadGraphWeakArcs(
        uint32_t maxExtReads = 3,
        double lenRatio = 0.975,
        uint32_t minDiff = 16);
    // Combined cleaning passes:
    void cleanBidirectionalReadGraphInitial(
        uint32_t gapFuzz = 1000,
        uint32_t maxShortTipReads = 3);
    void cleanBidirectionalReadGraphIterative(
        uint32_t cleanRounds = 4,
        double minDropRate = 0.2,
        double maxDropRate = 0.8,
        uint32_t maxShortTipReads = 3);
    void accessStringGraph();
    void accessStringGraphReadWrite();
    void checkStringGraphIsOpen() const;
    void accessUnitigGraph();
    void accessUnitigGraphReadWrite();
    void checkUnitigGraphIsOpen() const;
    void removeReadGraphBridges(uint64_t maxDistance);
    void analyzeReadGraph();
    void readGraphClustering();
    void writeReadGraphEdges(bool useReadName=false) const;

    void createReadGraph3(uint64_t maxAlignmentCount);
    void createReadGraph4(uint64_t maxAlignmentCount);

    void createReadGraph4AllAlignments(uint32_t maxAlignmentCount);
    void createReadGraph4withStrandSeparation(
        uint64_t maxAlignmentCount,
        double epsilon,
        double delta,
        double WThreshold,
        double WThresholdForBreaks
        );
    void removeReadGraph();

    void accessReadGraphAllAlignments();
    void accessReadGraphAllAlignmentsReadWrite();
    void checkReadGraphAllAlignmentsIsOpen() const;

    // Read graph creation method 6: Uses CIGAR-based phasing (isDeleted0/isDeleted1 flags)
    // instead of variant clustering. Provides Hifiasm-parity for ONT/HiFi data.
    void createReadGraph6();
    void createReadGraph6(uint64_t threadCount);

    // Create a read graph keeping all alignments without any filtering.
    // Use this together with marker graph vertex coverage thresholds
    // (minCoverage/maxCoverage in createMarkerGraphVertices) to filter
    // instead of pre-filtering per read.
    // Build a read graph from all alignments.
    // If pruneContained is true, contained reads (flagged by flagContainedReads)
    // keep only their single best alignment to a non-contained read.
    void createReadGraphAllAlignments(bool pruneContained = false);

    // Create a read graph using only the cis/trans (phasing) decisions produced by
    // performHifiasmECParity. This ignores all non-phasing deletion reasons.
    // An overlap is kept iff neither read marked it with DeleteReasonPhase.
    void createReadGraphFromEcParityCisOverlaps();
    void createReadGraphFromEcParityCisOverlaps(uint64_t threadCount, bool rebuildDirectedReadGraph);

    // Read graph from phaseOverlaps labels (hifiasmEcMatchState).
    // Keeps overlaps where neither side is TRANS (state 2).
    // Unlabeled overlaps (state 0) are kept.
    void createReadGraphFromPhasingCisOverlaps();
    void createReadGraphFromPhasingCisOverlaps(uint64_t threadCount, bool rebuildDirectedReadGraph);

    // asg_arc_del_trans port: remove transitive edges from the read graph.
    // An edge v→x is transitive if there exists v→w→x with
    // len(v→w) + len(w→x) <= len(v→x) + fuzz.
    uint64_t transitiveReductionOnReadGraph(int32_t fuzz = 5000);

    // Like createReadGraphFromEcParityCisOverlaps, but only keep cis overlaps that
    // cover at least one informative site (as recorded by AlignmentData::coversHetSite()
    // / informativeHetSiteCount{0,1} during performHifiasmECParity).
    void createReadGraphFromEcParityCisOverlapsCoveringInformativeSites();
    void createReadGraphFromEcParityCisOverlapsCoveringInformativeSites(uint64_t threadCount, bool rebuildDirectedReadGraph);
    
    // Canonical per-Read overlap storage  // Convert alignmentData to OverlapIndex
    
    // Create read graph directly from OverlapIndex (Option A: direct use, no alignmentData mapping)
    
	    // Create read graph from alignments after filterSecondaryAlignmentsPerReadPair (no phasing)
	    void createReadGraphFromFilteredAlignments();

	    // Rebuild read graph (and optionally directed read graph) from a provided keep vector.
	    // This removes the existing read graph data structures before recreating them.
	    void rebuildReadGraphUsingSelectedAlignments(vector<bool> keepAlignment, bool rebuildDirectedReadGraph = false);

    // Triangle and least square analysis of the read graph
    // to flag inconsistent alignments.
    void flagInconsistentAlignments(
        uint64_t triangleErrorThreshold,
        uint64_t leastSquareErrorThreshold,
        uint64_t leastSquareMaxDistance,
        uint64_t threadCount);
private:
    void flagInconsistentAlignmentsThreadFunction1(size_t threadId);
    void flagInconsistentAlignmentsThreadFunction2(size_t threadId);
    class FlagInconsistentAlignmentsData {
    public:

        // Arguments of flagInconsistentAlignments, stored here
        // to make them visible to the threads.
        uint64_t triangleErrorThreshold;
        uint64_t leastSquareErrorThreshold;
        uint64_t leastSquareMaxDistance;

        // The alignment offset for each edge of the read graph,
        // oriented with the lowest OrientedReadId first.
        MemoryMapped::Vector<int32_t> edgeOffset;

        // The inconsistent read graph edge ids found by each thread.
        vector< vector<uint64_t> > threadEdgeIds;
    };
    FlagInconsistentAlignmentsData flagInconsistentAlignmentsData;
public:



    // Functions and data used with read creation for iterative assembly.
    void createReadGraphUsingPseudoPaths(
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore,
        double mismatchSquareFactor,
        double minScore,
        uint64_t maxAlignmentCount,
        uint64_t threadCount);
    class CreateReadGraphsingPseudoPathsAlignmentData {
    public:
        uint64_t alignedMarkerCount = 0;

        // Pseudo-path alignment information.
        uint64_t weakMatchCount = 0;
        uint64_t strongMatchCount = 0;
        uint64_t mismatchCount = 0;
    };
    class CreateReadGraphUsingPseudoPathsData {
    public:
        int64_t matchScore;
        int64_t mismatchScore;
        int64_t gapScore;

        // The pseudopaths of all oriented reads.
        // Indexed by OrientedReadId::getValue().
        vector< vector<AssemblyGraphEdgeId> > pseudoPaths;

        // Vector to store information about each alignment.
        vector<CreateReadGraphsingPseudoPathsAlignmentData> alignmentInfos;
    };
    CreateReadGraphUsingPseudoPathsData createReadGraphUsingPseudoPathsData;

    // Thread function used to compute pseudoPaths.
    void createReadGraphUsingPseudoPathsThreadFunction1(size_t threadId);
    // Thread functions used to align pseudopaths.
    void createReadGraphUsingPseudoPathsThreadFunction2(size_t threadId);



#if 0
    // Functions and data for the version that uses mini-assemblies.
private:
    void createReadGraph2ThreadFunction(size_t threadId);
    void createReadGraph2LowLevel(ReadId);
    class CreateReadGraph2Data {
    public:
        vector<bool> keepAlignment;
    };
    CreateReadGraph2Data createReadGraph2Data;
public:
#endif



    // Approximate strand separation in the read graph.
    void flagCrossStrandReadGraphEdges1(int maxDistance, uint64_t threadCount);
private:
    void flagCrossStrandReadGraphEdges1ThreadFunction(size_t threadId);
    class FlagCrossStrandReadGraphEdges1Data {
    public:
        size_t maxDistance;
        vector<bool> isNearStrandJump;
    };
    FlagCrossStrandReadGraphEdges1Data flagCrossStrandReadGraphEdges1Data;
public:


    // Strict strand separation in the read graph.
    void flagCrossStrandReadGraphEdges2();
    void flagCrossStrandReadGraphEdges4();
    void flagCrossStrandReadGraphEdges5();



	    // Create the ReadGraph given a bool vector that specifies which
	    // alignments should be used in the read graph.
	    void createReadGraphUsingSelectedAlignments(vector<bool>& keepAlignment);
	    void createDirectedReadGraphUsingSelectedAlignments(vector<bool>& keepAlignment);
	    void createStringGraphUsingSelectedAlignments(const vector<bool>& keepAlignment);
    void createUnitigGraphFromStringGraph();
    void createReadGraphUsingAllAlignments(vector<bool>& keepAlignment);
    void rebuildReadGraphFromCurrentStringGraph(bool rebuildDirectedReadGraph = false);
    uint64_t reduceStringGraphTransitiveHifiasm(uint32_t gapFuzz);
    uint64_t cutStringGraphTips(uint32_t maxShortTipReads);
    void cleanStringGraphInitialHifiasm(uint32_t gapFuzz, uint32_t maxShortTipReads);
    void cleanStringGraphIterativeHifiasm(
        uint32_t cleanRounds,
        double minDropRate,
        double maxDropRate,
        uint32_t maxShortTipReads);
    void cleanStringGraphPreCleanHifiasm(uint32_t maxShortTipReads);
    void cleanStringGraphDropShortOverlaps(double lenRatio, uint32_t minOverlapLen, uint32_t maxShortTipReads);
	    void cleanStringGraphDropOverlapRoundsHifiasm(
	        uint32_t cleanRounds,
	        double minDropRate,
	        double maxDropRate,
	        uint32_t maxShortTipReads,
	        uint32_t finalMinOverlapLen);

		    // ONT-only hifiasm parity: weak arc cutting (ul_clean_gfa: asg_arc_cut_weak) on the StringGraph.
		    uint64_t cutStringGraphWeakArcsOntHifiasm(uint32_t maxExtReads, double lenRatio, uint32_t minDiff);
		    // Hifiasm parity: `asg_iterative_semi_circ` (semi-circular edge cutting, plus optional chimeric-bubble cut).
		    // `limLen` corresponds to hifiasm's LIM_LEN (100).
		    // `normalLen` corresponds to hifiasm's `normal_len` / `max_tip` (typically maxShortTipReads).
		    uint64_t cleanStringGraphBreakShortCycles(uint32_t limLen);
		    uint64_t cleanStringGraphBreakShortCycles(uint32_t limLen, uint32_t normalLen);
		    uint64_t cleanStringGraphRemoveSingleNodeBubbles(uint32_t maxShortTipReads);
		    uint64_t cleanStringGraphChimericReads();
		    uint64_t cleanStringGraphInexactOverlaps(uint32_t maxShortTipReads, uint32_t minDiff);
		    uint64_t cleanStringGraphBubbleLinks(double lenRat, double secLenRat, uint32_t maxShortTipReads);
		    uint64_t cleanStringGraphComplexBubbleLinks(double lenRat);
		    uint64_t cleanStringGraphLargeIndelArcs(uint32_t maxShortTipReads, uint32_t minDiff);
		    uint64_t cutStringGraphSemiCircular(uint32_t limLen);

	    uint64_t transitiveReduceUnitigGraph(uint32_t gapFuzz);
	    uint64_t cutUnitigGraphTips(uint32_t maxShortTipUnitigs);
	    uint64_t removeUnitigGraphOneStepBubbles();
	    void cleanUnitigGraphInitialHifiasm(uint32_t gapFuzz, uint32_t maxShortTipUnitigs);
	    void cleanUnitigGraphPreCleanHifiasm(uint32_t maxShortTipUnitigs);
	    void cleanUnitigGraphDropShortOverlaps(double dropRatio, uint32_t minOverlapLen);
	    void cleanUnitigGraphDropOverlapRoundsHifiasm(
        uint32_t cleanRounds,
        double minDropRate,
        double maxDropRate,
        uint32_t maxShortTipUnitigs,
        uint32_t finalMinOverlapLen);
    uint64_t cleanUnitigGraphBreakShortCycles(uint32_t maxCycleUnitigs);



public:
    // Use the read graph to flag chimeric reads.
    void flagChimericReads(size_t maxDistance, uint64_t threadCount);
private:
    class FlagChimericReadsData {
    public:
        size_t maxDistance;
    };
    FlagChimericReadsData flagChimericReadsData;
    void flagChimericReadsThreadFunction(size_t threadId);


	    // Create a local subgraph of the global read graph,
	    // starting at a given vertex and extending out to a specified
	    // distance (number of edges).
	    bool createLocalReadGraph(
	            OrientedReadId start,
	            uint32_t maxDistance,   // How far to go from starting oriented read.
	            bool allowChimericReads,
	            bool allowCrossStrandEdges,
	            bool allowInconsistentAlignmentEdges,
	            double timeout,         // Or 0 for no timeout.
	            LocalReadGraph&);

	    // Create a local subgraph of the global read graph,
	    // starting at any number of  given vertexes and extending out to a specified
	    // distance (number of edges).
	    bool createLocalReadGraph(
	        const vector<OrientedReadId>& starts,
	        uint32_t maxDistance,   // How far to go from starting oriented read.
	        bool allowChimericReads,
	        bool allowCrossStrandEdges,
	        bool allowInconsistentAlignmentEdges,
	        double timeout,         // Or 0 for no timeout.
	        LocalReadGraph&);

        // Create a local subgraph of the directed read graph,
        // starting at a given vertex and extending out to a specified
        // distance (number of arcs).
        // This is rendered as a LocalReadGraph and uses readGraph edge ids
        // to preserve the same visualization/analysis pipeline as exploreReadGraph.
        bool createLocalDirectedReadGraph(
                OrientedReadId start,
                uint32_t maxDistance,
                bool allowChimericReads,
                bool allowCrossStrandEdges,
                bool allowInconsistentAlignmentEdges,
                double timeout,
                LocalReadGraph&);

        bool createLocalDirectedReadGraph(
            const vector<OrientedReadId>& starts,
            uint32_t maxDistance,
            bool allowChimericReads,
            bool allowCrossStrandEdges,
            bool allowInconsistentAlignmentEdges,
            double timeout,
            LocalReadGraph&);

        // Create a local subgraph of the BidirectionalReadGraph.
        // Uses orientation-aware BFS via edge.traverse().
        // Produces a LocalReadGraph with derived OrientedReadId vertices.
        bool createLocalBidirectionalReadGraph(
            const vector<OrientedReadId>& starts,
            uint32_t maxDistance,
            bool allowChimericReads,
            bool allowInconsistentAlignmentEdges,
            double timeout,
            LocalReadGraph&);

	    bool createLocalStringGraph(
	        const vector<OrientedReadId>& starts,
	        uint32_t maxDistance,
	        bool allowChimericReads,
        bool followOutgoing,
        bool followIncoming,
        double timeout,
        LocalStringGraph&);

    bool createLocalUnitigGraph(
        const vector<OrientedUnitigId>& starts,
        uint32_t maxDistance,
        bool followOutgoing,
        bool followIncoming,
        double timeout,
        LocalUnitigGraph&);

    // Triangle analysis of the local read graph.
    // Returns a vector of triangles and their alignment residuals,
    // sorted by decreasing residual.
    void triangleAnalysis(
        LocalReadGraph&,
        LocalReadGraphTriangles&) const;

    // Singular value decomposition analysis of the local read graph.
    void leastSquareAnalysis(
        LocalReadGraph&,
        vector<double>& singularValues) const;

public:



    // Write a FASTA file containing all reads that appear in
    // the local read graph.
    void writeLocalReadGraphReads(
        ReadId,
        Strand,
        uint32_t maxDistance,
        bool allowChimericReads,
        bool allowCrossStrandEdges,
        bool allowInconsistentAlignmentEdges);



    AlignedEvidenceStore alignedEvidenceStore;

    // Per-overlap packed CIGARs (hifiasm-style uint16_t tokens).
    // Populated during computeBaseAlignmentsAndStore.
    // Token arena; each overlap's (cigarOffset, cigarTokenCount) in AlignmentInfo
    // points directly into this arena.
    OverlapCigarStore overlapCigarStore;

    void performHifiasmECFinalFilteringParity(uint64_t threadCount);

    // Old function (to be removed/replaced)
    void performHifiasmECFiltering(uint64_t threadCount);
    void performPhasingThreadFunction(uint64_t threadId);
    
    // Phasing using canonical OverlapIndex (sets is_match/strong flags)
    void performPhasingCanonical(uint64_t threadCount);
    void performPhasingCanonicalThreadFunction(uint64_t threadId);



    // Compute connected components of the read graph.
    // This just writes a csv file and has no other side effects
    // (nothing is stored).
    void computeReadGraphConnectedComponents() const;



    // Private functions and data used by createMarkerGraphVertices.
private:
    void createMarkerGraphVerticesThreadFunction1(size_t threadId);
    void createMarkerGraphVerticesThreadFunction2(size_t threadId);
    void createMarkerGraphVerticesThreadFunction21(size_t threadId);
    void createMarkerGraphVerticesThreadFunction3(size_t threadId);
    void createMarkerGraphVerticesThreadFunction4(size_t threadId);
    void createMarkerGraphVerticesThreadFunction5(size_t threadId);
    void createMarkerGraphVerticesThreadFunction45(int);
    void createMarkerGraphVerticesThreadFunction6(size_t threadId);
    void createMarkerGraphVerticesThreadFunction7(size_t threadId);
    void createMarkerGraphVerticesDebug1(uint64_t stage);
    class CreateMarkerGraphVerticesData {
    public:

        // Parameters.
        uint64_t minCoveragePerStrand;
        bool allowDuplicateMarkers;

        // The total number of oriented markers.
        uint64_t orientedMarkerCount;

        // Disjoint sets data structures.
        shared_ptr<DisjointSets> disjointSetsPointer;

        // The disjoint set that each oriented marker was assigned to.
        // See createMarkerGraphVertices for details.
        MemoryMapped::Vector<MarkerGraph::VertexId> disjointSetTable;

        // Work area used for multiple purposes.
        // See createMarkerGraphVertices for details.
        MemoryMapped::Vector<MarkerGraph::VertexId> workArea;

        // The markers in each disjoint set with coverage in the requested range.
        MemoryMapped::VectorOfVectors<MarkerId, MarkerGraph::VertexId> disjointSetMarkers;

        // Flag disjoint sets that contain more than one marker on the same oriented read.
        MemoryMapped::Vector<bool> isBadDisjointSet;

    };
    CreateMarkerGraphVerticesData createMarkerGraphVerticesData;



    void checkMarkerGraphVerticesAreAvailable() const;

    // Check for consistency of globalMarkerGraphVertex and globalMarkerGraphVertices.
    void checkMarkerGraphVertices(
        size_t minCoverage,
        size_t maxCoverage);



    // Marker graph.
public:
    MarkerGraph markerGraph;
    void removeMarkerGraph()
    {
        markerGraph.remove();
    }

    // Find the reverse complement of each marker graph vertex.
    void findMarkerGraphReverseComplementVertices(uint64_t threadCount);
    void accessMarkerGraphVertices(bool readWriteAccess = false);
    void accessMarkerGraphReverseComplementVertex(bool readWriteAccess = false);
    void removeMarkerGraphVertices();
    void accessDisjointSetsHistogram();
private:
    void findMarkerGraphReverseComplementVerticesThreadFunction1(size_t threadId);
    void findMarkerGraphReverseComplementVerticesThreadFunction2(size_t threadId);



    // Given a marker graph vertex, follow all of the contributing oriented
    // reads to their next vertex, but without moving forward more than
    // maxSkip markers.
    // In the returned vector, each entry correspond to a marker in the given vertex
    // (in the same order) and gives the next VertexId for that oriented read.
    // The next VertexId can be invalidVertexId if the oriented read has no vertices
    // past the starting VertexId.
    void findNextMarkerGraphVertices(
        MarkerGraphVertexId,
        uint32_t maxSkip,
        vector<MarkerGraphVertexId>&) const;

    // Find the common KmerId for all the markers of a marker graph vertex.
    KmerId getMarkerGraphVertexKmerId(MarkerGraphVertexId) const;

    // Clean up marker graph vertices that have duplicate markers
    // (more than one marker on the same oriented reads).
    // Such vertices are only generated when using --MarkerGraph.allowDuplicateMarkers.
public:
    void cleanupDuplicateMarkers(
        uint64_t threadCount,
        uint64_t minCoverage,
        uint64_t minCoveragePerStrand,
        double pattern1Threshold,
        bool pattern1CreateNewVertices,
        bool pattern2CreateNewVertices);
private:
    void cleanupDuplicateMarkersThreadFunction(size_t threadId);
    void cleanupDuplicateMarkersPattern1(
        MarkerGraph::VertexId,
        uint64_t minCoverage,
        uint64_t minCoveragePerStrand,
        bool createNewVertices,
        vector<MarkerDescriptor>&,
        vector<bool>& isDuplicateOrientedReadId,
        bool debug,
        ostream& out);
    void cleanupDuplicateMarkersPattern2(
        MarkerGraph::VertexId,
        uint64_t minCoverage,
        uint64_t minCoveragePerStrand,
        bool createNewVertices,
        vector<MarkerDescriptor>&,
        vector<bool>& isDuplicateOrientedReadId,
        bool debug,
        ostream& out);
    class CleanupDuplicateMarkersData {
    public:
        uint64_t minCoverage;
        uint64_t minCoveragePerStrand;
        double pattern1Threshold;
        bool pattern1CreateNewVertices;
        bool pattern2CreateNewVertices;
        uint64_t badVertexCount;    // Total number of vertices with duplicate markers.
        uint64_t pattern1Count;
        uint64_t pattern2Count;

        MarkerGraph::VertexId nextVertexId;

        // Get the next vertex id, then increment it in thread safe way.
        MarkerGraph::VertexId getAndIncrementNextVertexId()
        {
            return __sync_fetch_and_add(&nextVertexId, 1);
        }

    };
    CleanupDuplicateMarkersData cleanupDuplicateMarkersData;



    // Create marker graph edges.
public:
    void createMarkerGraphEdges(uint64_t threadCount);
    void accessMarkerGraphEdges(bool accessEdgesReadWrite, bool accessConnectivityReadWrite = false);
    void accessMarkerGraphEdgeMarkerIntervals();
    void checkMarkerGraphEdgesIsOpen() const;
    void accessMarkerGraphConsensus();
private:
    void createMarkerGraphEdgesThreadFunction0(size_t threadId);
    void createMarkerGraphEdgesThreadFunction1(size_t threadId);
    void createMarkerGraphEdgesThreadFunction2(size_t threadId);
    void createMarkerGraphEdgesThreadFunction12(size_t threadId, size_t pass);
    void createMarkerGraphEdgesBySourceAndTarget(uint64_t threadCount);
    class CreateMarkerGraphEdgesData {
    public:
        vector< shared_ptr< MemoryMapped::Vector<MarkerGraph::Edge> > > threadEdges;
        vector< shared_ptr< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> > > threadEdgeMarkerIntervals;
    };
    CreateMarkerGraphEdgesData createMarkerGraphEdgesData;



    // "Strict" version of createMarkerGraphEdges.
    // Differences from createMarkerGraphEdges:
    // - Will only create edges in which all contributing oriented reads have
    //   exactly the same RLE sequence. If more than one distinct RLE sequence
    //   is present, the edge is split into two parallel edges.
    // - Enforces minEdgeCoverage and minEdgeCoveragePerStrand.
    //   An edge is not generated if the total number of oriented
    //   reads on the edge is less than minEdgeCoverage,
    //   of it the number of oriented reads on each strand is less
    //   than minEdgeCoveragePerStrand.
    // - The main loop is written differently - it loops over reads
    //   rather than marker graph vertices.
    // Because of these strict criteria, this version generates frequent breaks
    // in contiguity that must later be fixed by other means.
public:
    void createMarkerGraphEdgesStrict(
        uint64_t minEdgeCoverage,
        uint64_t minEdgeCoveragePerStrand,
        uint64_t threadCount);
private:
    void createMarkerGraphEdgesStrictPass1(size_t threadId);
    void createMarkerGraphEdgesStrictPass2(size_t threadId);
    void createMarkerGraphEdgesStrictPass12(size_t threadId, uint64_t pass);
    void createMarkerGraphEdgesStrictPass3(size_t threadId);
    class CreateMarkerGraphEdgesStrictData {
    public:
        uint64_t minEdgeCoverage;
        uint64_t minEdgeCoveragePerStrand;

        // Marker intervals, stored in a VectorOfVectors
        // indexed by the source vertex id, vertexId0.
        // Each interval stores the target vertex, vertexId1.
        class MarkerIntervalInfo {
        public:
            MarkerGraph::VertexId vertexId1;
            OrientedReadId orientedReadId;
            uint32_t ordinal0;
            uint32_t ordinal1;
            bool operator<(const MarkerIntervalInfo& that) const
            {
                return tie(vertexId1, orientedReadId, ordinal0, ordinal1) <
                    tie(that.vertexId1, that.orientedReadId, that.ordinal0, that.ordinal1);
            }
        };
        MemoryMapped::VectorOfVectors<MarkerIntervalInfo, uint64_t> markerIntervalInfos;

        // The edges and corresponding marker intervals found by each thread.
        vector< shared_ptr< MemoryMapped::Vector<MarkerGraph::Edge> > > threadEdges;
        vector< shared_ptr< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> > > threadEdgeMarkerIntervals;

        // Class used in createMarkerGraphEdgesStrictPass3.
        class MarkerIntervalInfo3 {
        public:
            MarkerInterval markerInterval;
            uint32_t overlap;       // The number of overlapping bases between the markers.
            vector<Base> sequence;  // RLE
            bool operator<(const MarkerIntervalInfo3& that) const
            {
                return tie(overlap, sequence, markerInterval) <
                    tie(that.overlap, that.sequence, that.markerInterval);
            }
        };
    };
    CreateMarkerGraphEdgesStrictData createMarkerGraphEdgesStrictData;



    // Edge creation for Mode 3 assembly.
    // - Like createMarkerGraphEdgesStrict, this
    //   will only create edges in which all contributing oriented reads have
    //   exactly the same sequence. If more than one distinct sequence
    //   is present, the edge is split into two parallel edges.
    // - This only generates primary marker graph edges, defined as follows:
    //   * Edge coverage is >= minPrimaryCoverage and <= maxPrimaryCoverage.
    //   * Both its vertices have no duplicate oriented read.
    //   * The edge marker interval has no duplicate oriented read.
    // - This will only create the MarkerGraph::edgeMarkerIntervals
    //   and MarkerGraph::edgeSequence and nothing else.
public:
    void createPrimaryMarkerGraphEdges(
        uint64_t minPrimaryCoverage,
        uint64_t maxPrimaryCoverage,
        uint64_t threadCount);
private:
    void createPrimaryMarkerGraphEdgesThreadFunction(uint64_t threadId);
    class CreatePrimaryMarkerGraphEdgesData {
    public:
        uint64_t minPrimaryCoverage;
        uint64_t maxPrimaryCoverage;

        // The marker intervals of the edges found by each thread.
        vector< shared_ptr< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> > > threadMarkerIntervals;

        // The corresponding sequences
        vector< shared_ptr< MemoryMapped::VectorOfVectors<Base, uint64_t> > > threadSequences;
    };
    CreatePrimaryMarkerGraphEdgesData createPrimaryMarkerGraphEdgesData;



    // Write out the sets of parallel marker graph edges.
    // Only createMarkerGraphedgesStrict can create parallel edges.
public:
    void writeParallelMarkerGraphEdges() const;



    // Analyze and compare the read compositions of two marker graph edges.
    // This can only be done if the two edges have no duplicate OrientedReadIds
    // in the markers. In that case, each OrientedReadId of an edge
    // corresponds to one and only one markerInterval for each edge.
    bool analyzeMarkerGraphEdgePair(
        MarkerGraphEdgeId,
        MarkerGraphEdgeId,
        MarkerGraphEdgePairInfo&
        ) const;
    void writeHtmlMarkerGraphEdgePairInfo(
        ostream& html,
        MarkerGraphEdgeId,
        MarkerGraphEdgeId,
        const MarkerGraphEdgePairInfo&
        ) const;

    // Count the number of common oriented reads between two marker graph edges.
    // This assumes, WITHOUT CHECKING, that each of the two edges has no duplicate
    // oriented reads. This assumption is satisfied for primary marker graph edges
    // in Mode 3 assembly.
    uint64_t countCommonOrientedReadsUnsafe(MarkerGraphEdgeId, MarkerGraphEdgeId) const;

    // Estimate the offset, in bases, between two marker graph edges.
    // This assumes, WITHOUT CHECKING, that each of the two edges has no duplicate
    // oriented reads. This assumption is satisfied for primary marker graph edges
    // in Mode 3 assembly.
    // If there are common oriented reads between the two edges, this uses
    // countCommonOrientedReadsUnsafe.
    // This can fail, in which case it returns invalid<uint64_t>.
    uint64_t estimateBaseOffsetUnsafe(MarkerGraphEdgeId, MarkerGraphEdgeId) const;


    // Function createMarkerGraphSecondaryEdges can be called after createMarkerGraphEdgesStrict
    // to create a minimal amount of additional non-strict edges (secondary edges)
    // sufficient to restore contiguity.
    void createMarkerGraphSecondaryEdges(
        uint32_t secondaryEdgeMaxSkip,
        uint64_t threadCount);
private:
    void createMarkerGraphSecondaryEdges(
        uint32_t secondaryEdgeMaxSkip,
        bool aggressive,
        uint64_t threadCount);
public:



    // Cluster the oriented reads on a marker graph edge based on their sequence.
    // This returns a vector of connected components.
    // Each connected component is an index into the marker intervals for the edge.
    vector< vector<uint64_t> > clusterMarkerGraphEdgeOrientedReads(
        MarkerGraphEdgeId,
        double errorRateThreshold,
        bool debug) const;



    // Use clusterMarkerGraphEdgeOrientedReads to split secondary marker graph edges
    // where necessary.
    void splitMarkerGraphSecondaryEdges(
        double errorRateThreshold,
        uint64_t minCoverage,
        uint64_t threadCount);
    void splitMarkerGraphSecondaryEdgesThreadFunction(size_t threadId);
    class SplitMarkerGraphSecondaryEdgesData {
    public:
        double errorRateThreshold;
        uint64_t minCoverage;
        uint64_t initialSecondaryCount;
        uint64_t splitCount;
        uint64_t createdCount;

        // The new edges that were created by each thread.
        class Edge {
        public:
            MarkerGraphVertexId source;
            MarkerGraphVertexId target;
            vector<MarkerInterval> markerIntervals;
        };
        vector< vector<Edge> > threadEdges;
    };
    SplitMarkerGraphSecondaryEdgesData splitMarkerGraphSecondaryEdgesData;



    // Set marker graph edge flags to specified values for all marker graph edges.
    // Specify any value other than 0 or 1 leaves that flag unchanged.
    // Only useful for debugging.
    void setMarkerGraphEdgeFlags(
        uint8_t wasRemovedByTransitiveReduction,
        uint8_t wasPruned,
        uint8_t isSuperBubbleEdge,
        uint8_t isLowCoverageCrossEdge,
        uint8_t wasAssembled);



    // Find the reverse complement of each marker graph edge.
public:
    void findMarkerGraphReverseComplementEdges(uint64_t threadCount);
    void accessMarkerGraphReverseComplementEdge();
private:
    void findMarkerGraphReverseComplementEdgesThreadFunction1(size_t threadId);
    void findMarkerGraphReverseComplementEdgesThreadFunction2(size_t threadId);


    // Check that the marker graph is strand symmetric.
    // This can only be called after both findMarkerGraphReverseComplementVertices
    // and findMarkerGraphReverseComplementEdges have been called.
public:
    void checkMarkerGraphIsStrandSymmetric(uint64_t threadCount = 0);
private:
    void checkMarkerGraphIsStrandSymmetricThreadFunction1(size_t threadId);
    void checkMarkerGraphIsStrandSymmetricThreadFunction2(size_t threadId);



public:

    // Prune leaves from the strong subgraph of the global marker graph.
    void pruneMarkerGraphStrongSubgraph(size_t iterationCount);

private:



    // Private access functions for the global marker graph.
    // See the public section for some more that are callable from Python.

    // Find the vertex of the global marker graph that contains a given marker.
    // The marker is specified by the ReadId and Strand of the oriented read
    // it belongs to, plus the ordinal
    // The thread function runs on one thread at a time.
    void findMarkersSimdClosedSyncmersPass1(size_t threadId);
    void findMarkersSimdClosedSyncmersPass2(size_t threadId);
    class FindMarkersSimdClosedSyncmersData {
    public:
         int k;
         int w; // interpreted as s for syncmers
    };
    FindMarkersSimdClosedSyncmersData findMarkersSimdClosedSyncmersData;

    void findMarkersSimdMinimizersPass1(size_t threadId);
    void findMarkersSimdMinimizersPass2(size_t threadId);
    class FindMarkersSimdMinimizersData {
    public:
         int k;
         int w; // window size for minimizers
         bool useHifiasm = false; // use hifiasm's no-HPC sketcher as position source
         // Optional hifiasm overlap-path minimizer filter. When non-null, the
         // hifiasm sketch applies this high-occurrence k-mer filter plus
         // distance subsampling so markers match hifiasm's overlap seeds. Held
         // as void* to keep the C bridge header out of this C++ header; the
         // marker .cpp casts it back to hifiasm_filter_t*. Not owned here (the
         // driver in findMarkersSimdMinimizers owns and frees it).
         const void* hifiasmFilter = nullptr;
         int hifiasmSampleDist = 0; // subsampling distance; <=0 disables it
         // Per-read staged markers (forward strand: (position, kmerId)), filled
         // by the single sketch pass and consumed by the store pass so each read
         // is sketched exactly once. Indexed by ReadId; strand 1 is derived from
         // strand 0 during the store pass. Entries are freed as they are stored.
         std::vector<std::vector<std::pair<uint32_t, KmerId>>> stagedMarkers;
    };
    FindMarkersSimdMinimizersData findMarkersSimdMinimizersData;

public:
    // Prune existing markers based on KmerCounter frequencies.
    // filterRepeatKmers: also drop markers whose k-mer is a short-period tandem
    //   repeat (periods 1-6, thresholds {6,4,4,4,4,4}), the same predicate used
    //   by filterMarkerGraphVerticesByRepeatKmers.
    // filterLowComplexity: also drop markers whose k-mer is low-complexity by
    //   distinct sub-k-mer count (lengths 1-3, thresholds {4,12,24}), the same
    //   predicate used by filterMarkerGraphVerticesByDistinctSubkmerCount.
    // Applying these at the minimizer stage removes repeat/low-complexity
    // minimizers before marker-graph construction, so they never seed vertices.
    void applyKmerCountFilter(
        uint64_t minFreq, uint64_t maxFreq, uint64_t threadCount,
        bool filterPalindromes = true,
        bool filterRepeatKmers = false,
        bool filterLowComplexity = false);

    // Remove all markers from reads whose marker span covers less than
    // minSpanFraction of the read length. Span = lastMarkerPos + k - firstMarkerPos.
    void filterReadsByMarkerSpanCoverage(double minSpanFraction, uint64_t threadCount);

    // Remove markers whose k-mer appears more than once in the reference reads
    // (read IDs 0..referenceReadCount-1). Non-unique reference k-mers are
    // ambiguous anchors and cannot be used for read-to-reference alignment.
    // Remove reference k-mers that appear more than maxRefKmerFreq
    // times. Default 1 removes all non-unique k-mers. Higher values
    // retain low-copy repeat markers for VNTR SV detection.
    void removeNonUniqueReferenceMarkers(
        uint64_t referenceReadCount,
        uint64_t threadCount,
        uint64_t maxRefKmerFreq = 1);

    // Alignment candidates using Inverted Index (modular pipeline).
    // Phase 1-4: Build the inverted index for overlap candidate discovery.
    // buildCanonicalCache=true (default) stores a per-marker canonical k-mer
    // cache that lets chaining skip reverse-complement recomputation, at a cost
    // of ~17 bytes/marker of persistent RAM. Set false to save that memory on
    // low-RAM machines; the query phase then recomputes canonicalization inline.
    void buildInvertedIndex(uint64_t threadCount, bool buildCanonicalCache = true);

    // Phase 5: Run DP chaining on the built index to find alignment candidates.
    void chainAlignmentCandidates(
        double maxDriftRate,
        uint64_t maxChainLimit,
        const OverlapCandidatesOptions& overlapCandidatesOptions,
        uint64_t threadCount
    );

    // Convenience wrapper that calls both buildInvertedIndex and chainAlignmentCandidates.
    void findAlignmentCandidatesInvertedIndex(
        double maxDriftRate,
        uint64_t maxChainLimit,
        const OverlapCandidatesOptions& overlapCandidatesOptions,
        uint64_t threadCount
    );
public:
    // Compact Structure for Query Phase (8 bytes).
    struct CompactOccurrence {
        ReadId readId;
        uint32_t position;
    };

	    class AlignmentCandidatesInvertedIndexData {
	    public:
	         double maxDriftRate;
	         uint64_t k; // k-mer length for canonicalization
	         uint64_t coverageHet; // Added for Hifiasm Parity (Gradient Scoring)
             // InvertedIndex chaining configuration (see OverlapCandidatesOptions for meaning).
             double weightExponent = 1.1;
             double lowFreqMultiplier = 0.333;
             double highFreqMultiplier = 1.667;
             uint32_t rareKmerWeight = 2;
             bool downsampleHighFrequencyMarkers = true;
             uint32_t highFrequencySampleDistance = 500;
             uint32_t maxHighFrequencyPerStreak = 16;
             double highFactor = 5.0;        // Hifiasm high_factor: max_n_chain = max(hom_cov * high_factor, min_n_chain)
             uint32_t minNChain = 100;       // Hifiasm MIN_N_CHAIN: minimum max_n_chain value
	             double nonRedundantOverlapFraction = 0.5;
		             bool lchainIsAccurate = true;
		             bool useEcScoring = true;
	             bool enableMcopyFast = true;
             // Chaining scoring mode: 0 = hifiasm, 1 = minimap2-sr.
             int chainingMode = 0;
             // When > 0, only chain pairs where at least one read is a
             // reference (readId < referenceReadCount). Skips read-vs-read.
             uint64_t referenceReadCount = 0;
             int32_t minimap2Bw = 100;
             int32_t minimap2MaxGap = 100;
             int32_t minimap2MinChainScore = 25;
             uint32_t mcopyNum = 3;
             double mcopyRate = 0.70;
             uint32_t mcopyKhitCutoff = 32;
             uint32_t mcopyOcvWindow = 3072;
             double mcopyOcvWeakKeepRatio = 0.70;
             uint32_t minOverlapLength = 0;  // If >0, reject candidates whose min(qSpan, tSpan) < threshold.
             uint32_t maxEndFuzz = 0;        // If >0, reject candidates needing more extension to read ends.
             uint32_t maxChainingFreq = 1000;  // Skip kmers with frequency above this during hit collection.
	             vector<uint32_t> weightLut; // size 512 (pow(weightBase, weightExponent) truncated)
	         
         // Compact vector for Query (8 bytes/hit).
         vector<CompactOccurrence> compactOccurrences;

	         // Canonical k-mer ids for strand-0 markers, laid out read-contiguously.
	         // offsets[r]..offsets[r+1]-1 corresponds to read r strand 0.
	         // This cache lets chaining avoid per-marker reverse-complement work.
	         vector<uint64_t> strand0CanonicalOffsets;
	         vector<KmerId> strand0CanonicalKmerIds;
	         // For each entry in strand0CanonicalKmerIds, stores whether the observed k-mer on strand 0
	         // was the reverse complement of the canonical k-mer (1) or the canonical itself (0).
	         // This is the hifiasm z->rev equivalent needed to compute per-hit rev = z->rev ^ y->rev.
	         vector<uint8_t> strand0CanonicalIsRc;


         // Open Addressing Hash Table (Linear Probing).
         // Key = KmerId. Value = {Start, Count}.
         // Stored as a flat vector. size must be Power of 2.
         struct HashEntry {
             KmerId key;
             uint64_t start;
             uint32_t count;
             bool empty = true;
         };
         vector<HashEntry> hashTable;
    };

    AlignmentCandidatesInvertedIndexData invertedIndexData;

private:
    void applyKmerCountFilterThreadFunctionPass1(size_t threadId);
    void applyKmerCountFilterThreadFunctionPass2(size_t threadId);
    class ApplyKmerCountFilterData {
    public:
        uint64_t minFreq;
        uint64_t maxFreq;
        bool filterPalindromes;
        bool filterRepeatKmers;
        bool filterLowComplexity;
        shared_ptr<MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>> oldMarkers;
        shared_ptr<MemoryMapped::VectorOfVectors<KmerId, uint64_t>> oldMarkerKmerIds;
        
        // Optimization: Cache validity bits per read to avoid redundant work in Pass2.
        // Each vector<uint8_t> is a packed bitset (8 markers per byte).
        std::vector<std::vector<uint8_t>> markerValidity;
        
        // Cache read lengths to avoid accessing Reads in Pass2.
        std::vector<uint64_t> readLengths;
        
        // Cache reverse complement KmerIds per read, computed in Pass1.
        // Only valid markers have their rcKmerId stored (sparse storage by index).
        // This eliminates the reverseComplement() call in Pass2.
        std::vector<std::vector<KmerId>> rcKmerIds;
    };
    ApplyKmerCountFilterData applyKmerCountFilterData;

    // Use a k-mer checker to select markers.
    // If not used, kmerChecker is empty (null pointer).
    // If the marker is not contained in any vertex, return
    // MarkerGraph::invalidVertexId.
public:
    MarkerGraph::VertexId getGlobalMarkerGraphVertex(
        OrientedReadId,
        uint32_t ordinal) const;
private:

    // Get pairs (ordinal, marker graph vertex id) for all markers of an oriented read.
    // The pairs are returned sorted by ordinal.
    void getMarkerGraphVertices(
        OrientedReadId,
        vector< pair<uint32_t, MarkerGraph::VertexId> >&);

    // Find the markers contained in a given vertex of the global marker graph.
    // The markers are stored as pairs(oriented read id, ordinal).
    void getGlobalMarkerGraphVertexMarkers(
        MarkerGraph::VertexId,
        vector< pair<OrientedReadId, uint32_t> >&) const;

    void getGlobalMarkerGraphVertexChildren(
        MarkerGraphVertexId,
        vector< pair<MarkerGraphVertexId, vector<MarkerInterval> > >&,
        vector< pair<MarkerGraphVertexId, MarkerInterval> >& workArea
        ) const;

    // Given two marker graph vertices, get the marker intervals
    // that a possible edge between the two vertices would have.
    void getMarkerIntervals(
        MarkerGraphVertexId,
        MarkerGraphVertexId,
        vector<MarkerInterval>&
        ) const;

    // Return true if a vertex of the global marker graph has more than
    // one marker for at least one oriented read id.
    bool isBadMarkerGraphVertex(MarkerGraph::VertexId) const;

    // Write csv files with detailed marker graph information.
    void debugWriteMarkerGraph(const string& fileNamePrefix = "") const;

    // Write a csv file with information on all marker graph vertices for which
    // isBadMarkerGraphVertex returns true.
public:
    void writeBadMarkerGraphVertices() const;
private:

    // Find out if a vertex is a forward or backward leaf of the pruned
    // strong subgraph of the marker graph.
    // A forward leaf is a vertex with out-degree 0.
    // A backward leaf is a vertex with in-degree 0.
    bool isForwardLeafOfMarkerGraphPrunedStrongSubgraph(MarkerGraph::VertexId) const;
    bool isBackwardLeafOfMarkerGraphPrunedStrongSubgraph(MarkerGraph::VertexId) const;

    // Given an edge of the pruned strong subgraph of the marker graph,
    // return the next/previous edge in the linear chain the edge belongs to.
    // If the edge is the last/first edge in its linear chain, return MarkerGraph::invalidEdgeId.
    MarkerGraphEdgeId nextEdgeInMarkerGraphPrunedStrongSubgraphChain(MarkerGraphEdgeId) const;
    MarkerGraphEdgeId previousEdgeInMarkerGraphPrunedStrongSubgraphChain(MarkerGraphEdgeId) const;

    // Return the out-degree or in-degree (number of outgoing/incoming edges)
    // of a vertex of the pruned strong subgraph of the marker graph.
    size_t markerGraphPrunedStrongSubgraphOutDegree(MarkerGraph::VertexId) const;
    size_t markerGraphPrunedStrongSubgraphInDegree (MarkerGraph::VertexId) const;

    // Return true if an edge disconnects the local subgraph.
    bool markerGraphEdgeDisconnectsLocalStrongSubgraph(
        MarkerGraphEdgeId edgeId,
        size_t maxDistance,

        // Work areas, to reduce memory allocation activity.

        // Each of these two must be sized maxDistance+1.
        array<vector< vector<MarkerGraphEdgeId> >, 2>& verticesByDistance,

        // Each of these two must be sized globalMarkerGraphVertices.size()
        // and set to all false on entry.
        // It is left set to all false on exit, so it can be reused.
        array<vector<bool>, 2>& vertexFlags
        ) const;



    // Each oriented read corresponds to a path in the marker graph.
    // This function computes a subset of that path
    // covering the specified range of marker ordinals for the given
    // oriented read.
    void computeOrientedReadMarkerGraphPath(
        OrientedReadId,
        uint32_t firstOrdinal,
        uint32_t lastOrdinal,
        vector<MarkerGraphEdgeId>& path,
        vector< pair<uint32_t, uint32_t> >& pathOrdinals
        ) const;

    // Create the marker connectivity graph starting with a given marker.
    void createMarkerConnectivityGraph(
        OrientedReadId,
        uint32_t ordinal,
        bool useReadGraphAlignmentsOnly,
        MarkerConnectivityGraph&) const;
    void createMarkerConnectivityGraph(
        OrientedReadId,
        uint32_t ordinal,
        bool useReadGraphAlignmentsOnly,
        MarkerConnectivityGraph&,
        MarkerConnectivityGraphVertexMap&) const;

    // Compute an alignment between two oriented reads
    // induced by the marker graph. See InducedAlignment.hpp for more
    // information.
    void computeInducedAlignment(
        OrientedReadId,
        OrientedReadId,
        InducedAlignment&
    );

    // Compute induced alignments between an oriented read orientedReadId0
    // and the oriented reads stored sorted in orientedReadIds1.
    void computeInducedAlignments(
        OrientedReadId orientedReadId0,
        const vector<OrientedReadId>& orientedReadIds1,
        vector<InducedAlignment>& inducedAlignments);

    // Fill in compressed ordinals of an InducedAlignment.
    void fillCompressedOrdinals(
        OrientedReadId,
        OrientedReadId,
        InducedAlignment&);

    // Find the markers aligned to a given marker.
    // This is slow and cannot be used during assembly.
    void findAlignedMarkers(
        OrientedReadId, uint32_t ordinal,
        bool useReadGraphAlignmentsOnly,
        vector< pair<OrientedReadId, uint32_t> >&) const;



    // Extract a local subgraph of the global marker graph.
    bool extractLocalMarkerGraph(
        OrientedReadId,
        uint32_t ordinal,
        uint64_t distance,
        int timeout,                 // Or 0 for no timeout.
        uint64_t minVertexCoverage,
        uint64_t minEdgeCoverage,
        bool useWeakEdges,
        bool usePrunedEdges,
        bool useSuperBubbleEdges,
        bool useLowCoverageCrossEdges,
        bool useRemovedSecondaryEdges,
        LocalMarkerGraph0&
        );
    bool extractLocalMarkerGraph(
        MarkerGraph::VertexId,
        uint64_t distance,
        int timeout,                 // Or 0 for no timeout.
        uint64_t minVertexCoverage,
        uint64_t minEdgeCoverage,
        bool useWeakEdges,
        bool usePrunedEdges,
        bool useSuperBubbleEdges,
        bool useLowCoverageCrossEdges,
        bool useRemovedSecondaryEdges,
        LocalMarkerGraph0&
        );

    // Compute consensus sequence for a vertex of the marker graph.
    void computeMarkerGraphVertexConsensusSequence(
        MarkerGraph::VertexId,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts
        );



    // Class used to store spoa details needed by the http server.
    // See computeMarkerGraphEdgeConsensusSequenceUsingSpoa for details.
    class ComputeMarkerGraphEdgeConsensusSequenceUsingSpoaDetail {
    public:

        // If there is a very long marker interval,
        // the shortest sequence is used as consensus.
        // In that case, this flag is set and nothing else is stored.
        bool hasLongMarkerInterval;

        // Assembly mode: 1=overlapping bases, 2=intervening bases.
        int assemblyMode;   // 1 or 2.

        // Data stored when hasLongMarkerInterval is set.
        size_t iShortest;

        // Data stored for assembly mode 1.



        // Data stored for assembly mode 2.

        // The alignment for each distinct sequence.
        // Indexed by distinct sequence index.
        vector<string> msa;

        // The consensus, including gap bases.
        vector<AlignedBase> alignedConsensus;
        vector<uint8_t> alignedRepeatCounts;

        // The indexes of oriented reads that have each of the distinct sequences.
        // Indexed by distinct sequence index (same as the index used
        // for the msa vector above).
        vector< vector<size_t> > distinctSequenceOccurrences;

        // The alignment row corresponding to each oriented read.
        vector<int> alignmentRow;
    };



    // Use spoa to compute consensus sequence for an edge of the marker graph.
    // This does not include the bases corresponding to the flanking markers.
    void computeMarkerGraphEdgeConsensusSequenceUsingSpoa(
        MarkerGraphEdgeId,
        uint32_t markerGraphEdgeLengthThresholdForConsensus,
        const std::unique_ptr<spoa::AlignmentEngine>& spoaAlignmentEngine,
        spoa::Graph& spoaAlignmentGraph,
        vector<Base>& sequence,
        vector<uint32_t>& repeatCounts,
        uint8_t& overlappingBaseCount,
        ComputeMarkerGraphEdgeConsensusSequenceUsingSpoaDetail&,
        vector< pair<uint32_t, CompressedCoverageData> >* coverageData // Optional
        );



    // Simplify the marker graph.
    // The first argument is a number of marker graph edges.
    // See the code for detail on its meaning and how it is used.
public:
    void simplifyMarkerGraph(
        const vector<size_t>& maxLength, // One value for each iteration.
        bool debug);
private:
    void simplifyMarkerGraphIterationPart1(
        size_t iteration,
        size_t maxLength,
        bool debug);
    void simplifyMarkerGraphIterationPart2(
        size_t iteration,
        size_t maxLength,
        bool debug);



    // Create a coverage histogram for vertices and edges of the
    // marker graph. This counts all vertices that are not isolated
    // (are connected to no edges that are not marked removed)
    // and all edges that are not marked as removed.
    // Output is to csv files.
public:
    void computeMarkerGraphCoverageHistogram();


    // In the assembly graph, each vertex corresponds to a linear chain
    // of edges in the pruned strong subgraph of the marker graph.
    // A directed vertex A->B is created if the last marker graph vertex
    // of the edge chain corresponding to A coincides with the
    // first marker graph vertex of the edge chain corresponding to B.
    shared_ptr<mode0::AssemblyGraph> assemblyGraphPointer;
    void removeAssemblyGraph()
    {
        assemblyGraphPointer.reset();
    }
    void createAssemblyGraphVertices();
    void accessAssemblyGraphVertices();
    void createAssemblyGraphEdges();
    void accessAssemblyGraphEdgeLists();
    void accessAssemblyGraphEdges();
    void accessAssemblyGraphOrientedReadsByEdge();
    void writeAssemblyGraph(const string& fileName) const;
    void pruneAssemblyGraph(uint64_t pruneLength);

    // Gather and write out all reads that contributed to
    // each assembly graph edge.
    void gatherOrientedReadsByAssemblyGraphEdge(uint64_t threadCount);
    void writeOrientedReadsByAssemblyGraphEdge();
private:
    void gatherOrientedReadsByAssemblyGraphEdgePass1(size_t threadId);
    void gatherOrientedReadsByAssemblyGraphEdgePass2(size_t threadId);
    void gatherOrientedReadsByAssemblyGraphEdgePass(int);

    // Extract a local assembly graph from the global assembly graph.
    // This returns false if the timeout was exceeded.
    bool extractLocalAssemblyGraph(
        AssemblyGraphEdgeId,
        int distance,
        double timeout,
        LocalAssemblyGraph&) const;
public:
    void colorGfaBySimilarityToSegment(
        AssemblyGraphEdgeId,
        uint64_t minVertexCount,
        uint64_t minEdgeCount);


    // Compute consensus repeat counts for each vertex of the marker graph.
    void assembleMarkerGraphVertices(uint64_t threadCount);
    void accessMarkerGraphVertexRepeatCounts();
private:
    void assembleMarkerGraphVerticesThreadFunction(size_t threadId);
public:



    // Optional computation of coverage data for marker graph vertices.
    // This is only called if Assembly.storeCoverageData in dinara.conf is True.
    void computeMarkerGraphVerticesCoverageData(uint64_t threadCount);
private:
    void computeMarkerGraphVerticesCoverageDataThreadFunction(size_t threadId);
    class ComputeMarkerGraphVerticesCoverageDataData {
    public:

        // The results computed by each thread.
        // For each threadId:
        // - threadVertexIds[threadId] contains the vertex ids processed by each thread.
        // - threadVertexCoverageData[threadId] contains the coverage data for those vertices.
        vector< shared_ptr<
            MemoryMapped::Vector<MarkerGraph::VertexId> > > threadVertexIds;
        vector< shared_ptr<
            MemoryMapped::VectorOfVectors<pair<uint32_t, CompressedCoverageData>, uint64_t> > >
            threadVertexCoverageData;
    };
    ComputeMarkerGraphVerticesCoverageDataData computeMarkerGraphVerticesCoverageDataData;



    // Find the set of assembly graph edges encountered on a set
    // of edges in the marker graph. The given marker graph edges
    // could form a path, but don't have to.
    void findAssemblyGraphEdges(
        const vector<MarkerGraphEdgeId>& markerGraphEdges,
        vector<AssemblyGraphEdgeId>& assemblyGraphEdges
        ) const;



    // Pseudo-paths.
    // An oriented read corresponds to a path (sequence of adjacent edges)
    // in the marker graph, which
    // can be computed via computeOrientedReadMarkerGraphPath.
    // That path encounters a sequence of assembly graph edges,AssemblyGraph::EdgeId
    // which is not necessarily a path in the assembly graph
    // because not all marker graph edges belong to an assembly graph edge.
    // We call this sequence the pseudo-path of an oriented read in the assembly graph.
public:
    class PseudoPathEntry {
    public:
        AssemblyGraphEdgeId segmentId;

        // The first and last ordinal on the oriented read
        // where this assembly graph edge (segment) is encountered.
        uint32_t firstOrdinal;
        uint32_t lastOrdinal;

        // The first and last position in the assembly graph edge
        // (segment) where this oriented read is encountered.
        uint32_t firstPosition;
        uint32_t lastPosition;

        // The number of marker graph edges on the oriented read
        // where this assembly graph edge (segment) is encountered.
        uint32_t markerGraphEdgeCount;
    };
    using PseudoPath = vector<PseudoPathEntry>;
    void computePseudoPath(
        OrientedReadId,

        // The marker graph path computed using computeOrientedReadMarkerGraphPath.
        // This is computed by this function - it does not neet to be filled in
        // in advance.
        vector<MarkerGraphEdgeId>& path,
        vector< pair<uint32_t, uint32_t> >& pathOrdinals,

        // The pseudo-path computed by this function.
        PseudoPath&) const;
    void writePseudoPath(ReadId, Strand) const;
    static void getPseudoPathSegments(const PseudoPath&, vector<AssemblyGraphEdgeId>&);



    // Detangle the AssemblyGraph.
    void detangle();    // detangleMethod 1
    void detangle2(     // detangleMethod 2
        uint64_t diagonalReadCountMin,
        uint64_t offDiagonalReadCountMax,
        double detangleOffDiagonalRatio
         );



    // CompressedAssemblyGraph.
    // Note that we have no persistent version of this.
    // It must be created from scratch each time.
    void createCompressedAssemblyGraph();
    shared_ptr<CompressedAssemblyGraph> compressedAssemblyGraph;
    void colorCompressedAssemblyGraph(const string&);


public:
    // Mark as isLowCoverageCrossEdge all low coverage cross edges
    // of the assembly graph and the corresponding marker graph edges.
    // These edges are then considered removed.
    // An edge v0->v1 of the assembly graph is a cross edge if:
    // - in-degree(v0)=1, out-degree(v0)>1
    // - in-degree(v1)>1, out-degree(v1)=1
    // A cross edge is marked as isCrossEdge if its average edge coverage
    // is <= crossEdgeCoverageThreshold.
    void removeLowCoverageCrossEdges(uint32_t crossEdgeCoverageThreshold);



    // Assemble consensus sequence and repeat counts for each marker graph edge.
    void assembleMarkerGraphEdges(
        uint64_t threadCount,

        // This controls when we give up trying to compute consensus for long edges.
        uint32_t markerGraphEdgeLengthThresholdForConsensus,

        // Request storing detailed coverage information in binary format.
        bool storeCoverageData,

        // Request assembling all edges (used by Mode 2 assembly)
        bool assembleAllEdges
        );
private:
    void assembleMarkerGraphEdgesThreadFunction(size_t threadId);
    class AssembleMarkerGraphEdgesData {
    public:

        // The arguments to assembleMarkerGraphEdges, stored here so
        // they are accessible to the threads.
        uint32_t markerGraphEdgeLengthThresholdForConsensus;
        bool storeCoverageData;
        bool assembleAllEdges;

        // The results computed by each thread.
        // For each threadId:
        // threadEdgeIds[threadId] contains the edge ids processed by each thread.
        // threadEdgeConsensusSequence[threadId]  and
        // threadEdgeConsensusOverlappingBaseCount[threadId] contains the corresponding
        // consensus sequence and repeat counts.
        // These are temporary data which are eventually gathered into
        // MarkerGraph::edgeConsensus and MarkerGraph::edgeConsensusOverlappingBaseCount
        // before assembleMarkerGraphEdges completes.
        // See their definition for more details about their meaning.
        vector< shared_ptr< MemoryMapped::Vector<MarkerGraphEdgeId> > > threadEdgeIds;
        vector< shared_ptr< MemoryMapped::VectorOfVectors<pair<Base, uint8_t>, uint64_t> > > threadEdgeConsensus;
        vector< shared_ptr< MemoryMapped::Vector<uint8_t> > > threadEdgeConsensusOverlappingBaseCount;

        vector< shared_ptr<
            MemoryMapped::VectorOfVectors<pair<uint32_t, CompressedCoverageData>, uint64_t> > >
            threadEdgeCoverageData;
    };
    AssembleMarkerGraphEdgesData assembleMarkerGraphEdgesData;

    // Access coverage data for vertices and edges of the marker graph.
    // This is only available if the run had Assembly.storeCoverageData set to True
    // in dinara.conf.
public:
    void accessMarkerGraphCoverageData();
private:



    // Assemble Mode 3 sequence for all marker graph edges.
    // See the comments before MarkerGraph::edgeSequence for more information.
    // For now this is done sequentially.
public:
    void assembleMarkerGraphEdgesMode3();



    // Assemble sequence for an edge of the assembly graph.
private:
    void assembleAssemblyGraphEdge(
        AssemblyGraphEdgeId,
        bool storeCoverageData,
        AssembledSegment&);
    // Lower level version that works on a generic marker graph path.
    void assembleAssemblyGraphEdge(
        const span<const MarkerGraphEdgeId>&,
        bool storeCoverageData,
        AssembledSegment&);
public:
    AssembledSegment assembleAssemblyGraphEdge(
        AssemblyGraphEdgeId,
        bool storeCoverageData);


    // Assemble sequence for all edges of the assembly graph.
    void assemble(
        uint64_t threadCount,
        uint32_t storeCoverageDataCsvLengthThreshold);
    void accessAssemblyGraphSequences();
    void computeAssemblyStatistics();
private:
    class AssembleData {
    public:
        uint32_t storeCoverageDataCsvLengthThreshold;

        // The results created by each thread.
        // All indexed by threadId.
        vector< vector<AssemblyGraphEdgeId> > edges;
        vector< shared_ptr<LongBaseSequences> > sequences;
        vector< shared_ptr<MemoryMapped::VectorOfVectors<uint8_t, uint64_t> > > repeatCounts;
        void allocate(uint64_t threadCount);
        void free();
    };
    AssembleData assembleData;
    void assembleThreadFunction(size_t threadId);



    // Write the assembly graph in GFA 1.0 format defined here:
    // https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
public:
    void writeGfa1(const string& fileName);
    void writeGfa1BothStrands(const string& fileName);
    void writeGfa1BothStrandsNoSequence(const string& fileName);
private:
    // Construct the CIGAR string given two vectors of repeat counts.
    // Used by writeGfa1.
    static void constructCigarString(
        const span<uint8_t>& repeatCounts0,
        const span<uint8_t>& repeatCounts1,
        string&
        );

public:

    // Write assembled sequences in FASTA format.
    void writeFasta(const string& fileName);



    // Write a csv file that can be used to color the double-stranded GFA
    // in Bandage based on the presence of two oriented reads
    // on each assembly graph edge.
    // Red    =  only oriented read id 0 is present
    // Blue   =  only oriented read id 1 is present
    // Purple =  both oriented read id 0 and oriented read id 1 are present
    // Grey   =  neither oriented read id 0 nor oriented read id 1 are present
    void colorGfaWithTwoReads(
        ReadId readId0, Strand strand0,
        ReadId readId1, Strand strand1,
        const string& fileName
        ) const;



    // Color key segments in the gfa file.
    // A segment (assembly graph edge) v0->v1 is a key segment if in-degree(v0)<2 and
    // out_degree(v)<2, that is, there is no uncertainty on what preceeds
    // and follows the segment.
    void colorGfaKeySegments(const string& fileName) const;



    // Write a csv file describing the marker graph path corresponding to an
    // oriented read and the corresponding pseudo-path on the assembly graph.
    void writeOrientedReadPath(ReadId, Strand, const string& fileName) const;



    // Data and functions used for the http server.
    // This function puts the server into an endless loop
    // of processing requests.
    void writeHtmlBegin(ostream&) const;
    void writeHtmlEnd(ostream&) const;
    void writeAssemblySummary(ostream&);
    void writeAssemblySummaryBody(ostream&);
    void writeAssemblySummaryJson(ostream&);
    void writeAssemblyIndex(ostream&) const;
    static void writeStyle(ostream& html);


    void writeNavigation(ostream&) const;
    void writeNavigation(
        ostream& html,
        const string& title,
        const vector<pair <string, string> >&) const;

    static void writePngToHtml(
        ostream& html,
        const string& pngFileName,
        const string useMap = ""
        );
    static void writeGnuPlotPngToHtml(
        ostream& html,
        int width,
        int height,
        const string& gnuplotCommands);

    void fillServerFunctionTable();
    void processRequest(
        const vector<string>& request,
        ostream&,
        const BrowserInformation&) override;
    void exploreSummary(const vector<string>&, ostream&);
    void exploreRead(const vector<string>&, ostream&);
    void exploreReadRaw(const vector<string>&, ostream&);
    void exploreReadRle(const vector<string>&, ostream&);
    void exploreLookupRead(const vector<string>&, ostream&);
    void exploreReadSequence(const vector<string>&, ostream&);
    void exploreReadMarkers(const vector<string>&, ostream&);
    void exploreMarkerKmers(const vector<string>&, ostream&);
    void blastRead(const vector<string>&, ostream&);
    void exploreAlignmentCandidateGraph(const vector<string>& request, ostream& html);
    void exploreAlignments(const vector<string>&, ostream&);
    void exploreAlignmentCoverage(const vector<string>&, ostream&);
    void exploreAlignment(const vector<string>&, ostream&);
    void alignSequencesInBaseRepresentation(const vector<string>&, ostream&);
    void exploreAlignmentGraph(const vector<string>&, ostream&);
    void exploreReadGraph(const vector<string>&, ostream&);
    void exploreBidirectionalReadGraph(const vector<string>&, ostream&);
    void exploreStringGraph(const vector<string>&, ostream&);
    void exploreUnitigGraph(const vector<string>&, ostream&);
    void exploreUndirectedReadGraph(const vector<string>&, ostream&);
    void exploreDirectedReadGraph(const vector<string>&, ostream&);
    void exploreCompressedAssemblyGraph(const vector<string>&, ostream&);
    static bool parseCommaSeparatedReadIDs(string& commaSeparatedReadIds, vector<OrientedReadId>& readIds, ostream& html);
    static void addScaleSvgButtons(ostream&, uint64_t sizePixels);
    class HttpServerData {
    public:
        shared_ptr<LocalAlignmentCandidateGraph> referenceOverlapGraph;

        using ServerFunction = void (Assembler::*) (
            const vector<string>& request,
            ostream&);
        std::map<string, ServerFunction> functionTable;
        string docsDirectory;
        string referenceFastaFileName = "reference.fa";

        const AssemblerOptions* assemblerOptions = 0;

        void createGraphEdgesFromOverlapMap(const ReferenceOverlapMap& overlapMap);

    };
    HttpServerData httpServerData;

    // For the display of the alignment candidate graph, we can optionally
    // specify a PAF file containing alignments of reads to the reference.
    // Persistent data structures from loading the PAF are stored as
    // members of HttpServerData
    void loadAlignmentsPafFile(const string& alignmentsPafFileAbsolutePath);

    // Display alignments in an html table.
    void displayAlignments(
        OrientedReadId,
        const vector< pair<OrientedReadId, AlignmentInfo> >&,
        bool showIsInReadGraphFlag,
        ostream&,
        const vector<uint32_t>* alignmentIds = nullptr) const;
    void displayAlignment(
        OrientedReadId orientedReadId0,
        OrientedReadId orientedReadId1,
        const AlignmentInfo& alignment,
        ostream&) const;


    // Functions and data used by the http server
    // for display of the local marker graph.
    void exploreMarkerGraph0(const vector<string>&, ostream&);
    void exploreMarkerGraph1(const vector<string>&, ostream&);
    void getLocalMarkerGraph0RequestParameters(
        const vector<string>&,
        LocalMarkerGraph0RequestParameters&) const;
    void exploreMarkerGraphVertex(const vector<string>&, ostream&);
    void exploreMarkerGraphEdge(const vector<string>&, ostream&);
    void exploreMarkerGraphEdgePair(const vector<string>&, ostream&);
    void exploreMarkerCoverage(const vector<string>&, ostream&);
    void exploreMarkerGraphInducedAlignment(const vector<string>&, ostream&);
    void followReadInMarkerGraph(const vector<string>&, ostream&);
    void exploreMarkerConnectivity(const vector<string>&, ostream&);
    void renderEditableAlignmentConfig(
        const int method,
        const uint64_t maxSkip,
        const uint64_t maxDrift,
        const uint32_t maxMarkerFrequency,
        const uint64_t minAlignedMarkerCount,
        const double minAlignedFraction,
        const uint64_t maxTrim,
        const int matchScore,
        const int mismatchScore,
        const int gapScore,
        const double downsamplingFactor,
        int bandExtend,
        int maxBand,
        uint64_t align4DeltaX,
        uint64_t align4DeltaY,
        uint64_t align4MinEntryCountPerCell,
        uint64_t align4MaxDistanceFromBoundary,
        double align5DriftRateTolerance,
        uint64_t align5MinBandExtend,
        const Align6Options&,
        ostream& html
    );
    void writeColorPicker(ostream& html, string svgId);

    void writeMakeAllTablesCopyable(ostream&) const;

    // Do bulk sampling of reads and accumulate stats about their alignments
    void assessAlignments(const vector<string>& request, ostream& html);
    void sampleReads(vector<OrientedReadId>& sample, uint64_t n);
    void sampleReads(vector<OrientedReadId>& sample, uint64_t n, uint64_t minLength, uint64_t maxLength);
    void sampleReadsFromDeadEnds(
            vector<OrientedReadId>& sample,
            vector<bool>& isLeftEnd,
            uint64_t n);

    void sampleReadsFromDeadEnds(
            vector<OrientedReadId>& sample,
            vector<bool>& isLeftEnd,
            uint64_t n,
            uint64_t minLength,
            uint64_t maxLength);

    void countDeadEndOverhangs(
            const vector<pair<OrientedReadId, AlignmentInfo> >& allAlignmentInfo,
            const vector<bool>& isLeftEnd,
            Histogram2& overhangLengths,
            uint32_t minOverhang);

    // Compute all alignments for a given read.
    // This can be slow for large assemblies,
    // and therefore the computation is multithreaded.
    void computeAllAlignments(const vector<string>&, ostream&);
    void computeAllAlignmentsThreadFunction(size_t threadId);
    class ComputeAllAlignmentsData {
    public:
        OrientedReadId orientedReadId0;
        size_t minMarkerCount;
        size_t maxSkip;
        size_t maxDrift;
        uint32_t maxMarkerFrequency;
        size_t minAlignedMarkerCount;
        double minAlignedFraction;
        size_t maxTrim;
        int method;
        int matchScore;
        int mismatchScore;
        int gapScore;
        double downsamplingFactor;
        int bandExtend;
        int maxBand;
        uint64_t align4DeltaX;
        uint64_t align4DeltaY;
        uint64_t align4MinEntryCountPerCell;
        uint64_t align4MaxDistanceFromBoundary;
        double align5DriftRateTolerance;
        uint64_t align5MinBandExtend;
        Align6Options align6Options;
        // The alignments found by each thread.
        vector< vector< pair<OrientedReadId, AlignmentInfo> > > threadAlignments;
    };
    ComputeAllAlignmentsData computeAllAlignmentsData;


    // Access all available assembly data, without thorwing an exception
    // on failures.
public:
    void accessAllSoft();

    // Store assembly time.
    void storeAssemblyTime(
        double elapsedTimeSeconds,
        double averageCpuUtilization);

    void storePeakMemoryUsage(uint64_t peakMemoryUsage);

    // Functions and data used by the http server
    // for display of the local assembly graph.
private:
    void exploreAssemblyGraph(const vector<string>&, ostream&);
    class LocalAssemblyGraphRequestParameters {
    public:
        AssemblyGraphEdgeId edgeId;
        bool edgeIdIsPresent;
        uint32_t maxDistance;
        bool maxDistanceIsPresent;
        bool useDotLayout;
        bool showVertexLabels;
        bool showEdgeLabels;
        uint32_t sizePixels;
        bool sizePixelsIsPresent;
        double timeout;
        bool timeoutIsPresent;
        void writeForm(ostream&, AssemblyGraphEdgeId edgeCount) const;
        bool hasMissingRequiredParameters() const;
    };
    void getLocalAssemblyGraphRequestParameters(
        const vector<string>&,
        LocalAssemblyGraphRequestParameters&) const;
    void exploreAssemblyGraphEdge(const vector<string>&, ostream&);
    void exploreAssemblyGraphEdgesSupport(const vector<string>&, ostream&);


public:
    // Functions and data to save binary data.
    // This happens under control of --saveBinaryData
    // and is only implemented for Mode 3 assembly.

    // Binary data will be saved to this directory,
    // stored including the final slash:
    // Data/ if memoryMode is anonymous.
    // DataOnDisk if memoryMode is filesystem.
    string saveBinaryDataDirectory;
    void createSaveBinaryDataDirectory(const string& memoryMode);

    using SaveBinaryDataFunction = void (Assembler::*)() const;
    void initiateSaveBinaryData(SaveBinaryDataFunction);
    vector<std::thread> saveBinaryDataThreads;
    void waitForSaveBinaryDataThreads();

    void saveMarkers() const;


private:
    // Set up the ConsensusCaller used to compute the "best"
    // base and repeat count at each assembly position.
    // The argument to setupConsensusCaller specifies
    // the consensus caller to be used.
    // It can be one of the following:
    // - Modal
    //   Selects the SimpleConsensusCaller.
    // - Median
    //   Selects the MedianConsensusCaller.
    // - Bayesian:fileName
    //   Selects the SimpleBayesianConsensusCaller,
    //   using fileName as the configuration file.
    //   Filename must be an absolute path (it must begin with "/").
public:
    void setupConsensusCaller(const string&);
private:
    shared_ptr<ConsensusCaller> consensusCaller;
public:



    // Assembly graph for mode 2 assembly.
    shared_ptr<AssemblyGraph2> assemblyGraph2Pointer;
    void createAssemblyGraph2(
        uint64_t pruneLength,
        const Mode2AssemblyOptions&,
        uint64_t threadCount,
        bool debug);



    // Mode 3 assembly.
    shared_ptr<Mode3Assembler> mode3Assembler;
    void accessMode3Assembler();

    // BRG-native anchors (stored for HTTP server visualization).
    shared_ptr<mode3::BidirectedAnchors> bidirectedAnchors;
    void accessBidirectedAnchors();

    // Shasta2-style Anchors for Mode 3.
    std::shared_ptr<Shasta2Anchors> shasta2Anchors;
    std::shared_ptr<Shasta2Journeys> shasta2Journeys;

    // Early structural layer before het detection / phasing: per-strand-0-read CSR of
    // oriented reads that co-occur on any Shasta2 journey anchor (marker-graph collapse;
    // can exceed alignmentTable neighbors). Use to scope windows, anchor separation, evidence.
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> strand0JourneyCoReads;
    void computeStrand0JourneyCoReadsTable();

    // Non-overlapping journeys derived from shasta2Journeys.
    // Indexed by OrientedReadId.getValue(). Computed after journey creation.
    // AnchorId is uint64_t (shasta2::AnchorId).
    std::vector<std::vector<uint64_t>> shasta2LinearJourneys;

    // MSA het-site variant events (AssemblerMSAHetSites.cpp).
    // Produced by vg deconstruct on per-anchor-pair POA GFAs with embedded read paths.
    // Each event corresponds to one VCF record (one snarl, one alt allele).
    // Nested snarls are linked via parentSiteId / level fields.

    enum class MSAVariantType : uint8_t { SNP, INSERTION, DELETION, MNP, COMPLEX };

    struct VariantEvent {
        // Anchor-pair context
        uint32_t    segmentIndex;   // jPos in linear journey
        uint32_t    focalPosStart;  // base position of anchor[jPos]   in focal read
        uint32_t    focalPosEnd;    // base position of anchor[jPos+1] in focal read + k

        // VCF coordinates (relative to focal read as reference path)
        uint32_t    vcfPos;         // 1-based POS in focal read
        std::string vcfId;          // snarl ID from vg (e.g. ">2>5")

        // Alleles
        std::string refAllele;      // REF (focal-read allele)
        std::string altAllele;      // ALT (one alt; multi-allelic sites produce multiple events)
        MSAVariantType varType;

        // Snarl-tree linkage (from vg deconstruct -n LV/PS/PA/RS/RD tags)
        int         level;          // 0 = top-level; >0 = nested
        std::string parentSiteId;   // PS: ID of parent snarl (empty if level==0)
        int         parentAllele;   // PA: which allele of parent contains this site's ref path.
                                    //     0 = focal read passes through this nested site
                                    //     (it is inside the parent's REF branch).
                                    //     >0 = focal read does NOT enter this site
                                    //     (it is inside a parent ALT branch).
                                    //     -1 = not set (level==0).
        uint32_t    topLevelPosStart; // RS: start of top-level containing site on focal read
        uint32_t    topLevelPosEnd;   // RD: end   of top-level containing site on focal read

        // Spanning reads assigned by vg deconstruct GT column.
        // Star-allele reads (GT=*) are excluded — they carry no evidence at this level.
        std::vector<OrientedReadId> refReads;  // GT == 0
        std::vector<OrientedReadId> altReads;  // GT == altIdx for this event
    };

    // All variant events for a given oriented read, across all anchor-pair segments.
    // Indexed by OrientedReadId.getValue().
    // Ordered by (segmentIndex, vcfPos, level) — top-level events before children.
    // This is the primary input for het-site phasing.
    std::vector<std::vector<VariantEvent>> shasta2VariantEvents;

    // Per-anchor-pair POA-MSA for het-site detection (AssemblerMSAHetSites.cpp).
    // Writes one GFA per inter-anchor segment of the given read's linear journey.
    void computeMSAHetSites(ReadId focalReadId, uint32_t strand);
    void computeTheseusMarkerGraphMSAPrototype(
        uint64_t maxAnchorPairs,
        uint64_t maxReadsPerPair,
        uint64_t threadCount);
    void computeTheseusTargetBackboneMSAPrototype(
        uint64_t maxReads,
        uint64_t threadCount);
    void computeTheseusReadWindowMSAPrototype(
        uint64_t threadCount);
    void computeTheseusReadWindowMSAPrototype(
        shared_ptr<Shasta2Anchors> shasta2Anchors,
        shared_ptr<Shasta2Journeys> shasta2Journeys,
        uint64_t threadCount);
    // Partition anchor journeys into disjoint windows.
    // Clean version: for each touched read, keeps only anchors shared with
    // the backbone (by anchor ID), enforces backbone order via LIS, and
    // discards non-shared and out-of-order anchors. Guarantees strand-consistent
    // edges in the resulting AnchorGraph.
    // anchorDovetailWindow (optional output): when non-null, receives the
    // forward-oriented anchorId -> owning windowId map for claimed dovetail
    // anchors (empty when whole-journey claiming is disabled). Lets the anchor
    // graph treat dovetails as part of the window for path-finding.
    // tileUnclaimedIntervals (optional): when true, after the first pass that
    // seeds pristine full-journey cores, run a second priority-queue pass that
    // tiles the leftover unclaimed intervals into fragment windows (longest
    // base span first). When false (default), only full-journey cores are
    // created and the seams between them are left unclaimed.
    // threadCount is accepted but unused: the claiming sweep is a single
    // priority queue over shared anchor ownership with a strict largest-first
    // order, which is inherently serial (not just unparallelized) -- see the
    // definition for why. Kept for parity with sibling mode3-pipeline calls.
    void computeAnchorWindowsClean(
        shared_ptr<Shasta2Anchors> shasta2Anchors,
        shared_ptr<Shasta2Journeys> shasta2Journeys,
        const vector<ReadId>& readIdsSortedByLength,
        vector<AnchorWindow>& anchorWindows,
        uint64_t threadCount,
        uint64_t minCommonForBackbone = 2,
        uint64_t maxSkipForBackbone = 10,
        uint64_t minWindowBaseSpan = 4000,
        vector<uint32_t>* anchorDovetailWindow = nullptr,
        bool tileUnclaimedIntervals = false);

    // Convenience wrapper: compute disjoint full-journey cores, then tile the
    // leftover unclaimed intervals into fragment windows (longest first).
    // Equivalent to computeAnchorWindowsClean(..., tileUnclaimedIntervals=true).
    void computeAnchorWindowsWithUnclaimed(
        shared_ptr<Shasta2Anchors> shasta2Anchors,
        shared_ptr<Shasta2Journeys> shasta2Journeys,
        const vector<ReadId>& readIdsSortedByLength,
        vector<AnchorWindow>& anchorWindows,
        uint64_t threadCount,
        uint64_t minCommonForBackbone = 2,
        uint64_t maxSkipForBackbone = 10,
        uint64_t minWindowBaseSpan = 4000,
        vector<uint32_t>* anchorDovetailWindow = nullptr);

    // Per-window progressive abPOA.
    // For each window, seed an abPOA graph with the backbone read sequence,
    // then progressively align every other member read to the subgraph
    // spanned by the backbone anchors it shares with the backbone (its
    // anchor interval). Writes one GFA per window for diagnostics; the
    // abPOA graph / MSA is the substrate for later het-site detection.
    void computeWindowAbpoaGraphs(
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys,
        const string& outputPrefix,
        uint64_t threadCount) const;

    // Leaf-snarl (multi-column MNP/het-block) detection per anchor window,
    // from independent pairwise ksw2 alignments (no shared multi-sequence
    // graph, no per-interval fragmentation): one banded ksw2 global alignment
    // per inter-anchor segment per member, FULLSPAN (every column the member
    // spans gets a real aligned entry). Verification only -- writes nothing
    // to AnchorWindow::hetBubbles.
    void computeWindowKsw2LeafSnarls(
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys,
        const AlignOptions& alignOptions,
        uint64_t threadCount) const;

    // Leaf-snarl detection using ProjectedAlignment (the same fast pairwise
    // aligner computeBaseAlignmentsAndStore uses for the whole-dataset
    // overlap graph) instead of ksw2's per-segment banded DP: ONE
    // ProjectedAlignment call per member (covering all its inter-anchor
    // segments internally), not one per segment. Verification only.
    void computeWindowProjectedAlignmentLeafSnarls(
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys,
        const AlignOptions& alignOptions,
        uint64_t threadCount) const;

    // PRODUCTION het-bubble detection using ProjectedAlignment: same sparse
    // evidence/pinning machinery as computeWindowProjectedAlignmentLeafSnarls
    // above, but emits real AnchorWindow::hetBubbles (window.hetBubbles is
    // cleared and repopulated) instead of only logging -- a drop-in
    // alternative engine, selected the same way as ksw2/abpoa/intervalpoa.
    // Multi-column MNP snarls are skipped (HetAnchor::alleleBase is a single
    // base and can't represent them without extending the anchor format).
    // Returns total bubbles; also reports hetWindows and totalBubbles counts.
    uint32_t projAlnDetectHetBubblesAllWindows(
        vector<AnchorWindow>& windows,
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys,
        const AlignOptions& alignOptions,
        double hetMinVaf,
        uint64_t hetMinSupport,
        bool hetDropHomopolymer,
        bool hetDropRepeat,
        uint64_t threadCount,
        uint64_t& hetWindowsOut,
        uint64_t& totalBubblesOut,
        const vector<bool>* skipWindow = nullptr) const;

    // Detect clean het SNPs in an anchor window using Theseus MSA.
    // Returns the number of SNPs passing strand bias and repeat filtering.
    uint32_t msaDetectSnpsInWindow(
        const AnchorWindow& window,
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys) const;

    // Detect clean het SNPs in an anchor window using CIGAR-based variant parsing.
    // Uses pairwise CIGARs from OverlapCigarStore instead of Theseus MSA.
    uint32_t cigarDetectSnpsInWindow(
        AnchorWindow& window,
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys) const;

    // Detect clean het SNPs in an anchor window using banded ksw2 (2-piece
    // affine) alignment of each member's inter-anchor segments against the
    // backbone. Self-contained: uses persisted shared-anchor pins, so it works
    // for transitive members and needs no global alignment table. SNPs only.
    // noisyRegSlideWin / noisyRegMaxXgaps control the per-read CIGAR-density
    // noise filter (pgphase-style): within a window of noisyRegSlideWin backbone
    // bases, if the summed mismatch+indel size exceeds noisyRegMaxXgaps the span
    // is flagged noisy and its SNP votes are excluded. Tune per read technology
    // (HiFi ~100/5, ONT ~25/5).
    uint32_t ksw2DetectSnpsInWindow(
        AnchorWindow& window,
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys,
        const AlignOptions& alignOptions,
        int noisyRegSlideWin = 100,
        int noisyRegMaxXgaps = 5) const;

    // Detect clean het BUBBLES in an anchor window using the same banded ksw2
    // member pileup as ksw2DetectSnpsInWindow, but emit AnchorWindow::hetBubbles
    // in the exact format produced by testAbpoaMultiSegmentMSA (arms + leadHom +
    // hom, rawPositions = absolute oriented read base positions at the k=2
    // anchor columns, backboneOffset = bbPos - windowBbBegin). This is a drop-in
    // alternative to the abPOA het path: the downstream plan/append/stage passes
    // consume the bubbles unchanged. Returns the number of bubbles emitted.
    uint32_t ksw2DetectHetBubblesInWindow(
        AnchorWindow& window,
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys,
        const AlignOptions& alignOptions,
        double hetMinVaf,
        uint64_t hetMinSupport,
        bool hetDropHomopolymer,
        bool hetDropRepeat,
        int noisyRegSlideWin = 100,
        int noisyRegMaxXgaps = 5) const;

    // Shasta2 LocalAssembly7-style het detection: instead of one POA over the
    // whole window (testAbpoaMultiSegmentMSA, which builds a single growing
    // graph per window and serializes on one thread), this tiles the window at
    // its backbone anchors and runs an INDEPENDENT small POA (spoa) per
    // consecutive-anchor interval. Each interval's reads are aligned to each
    // other (true POA quality, unlike ksw2's star alignment), the graphs stay
    // tiny (one inter-anchor gap, not the whole 100kb window), and the intervals
    // are embarrassingly parallel. Per-interval MSA columns are mapped back to
    // absolute backbone offsets and merged into per-read KwMemberProfile rows,
    // then the shared emitHetBubblesFromProfiles tail emits the identical
    // AnchorWindow::hetBubbles. Returns the number of bubbles emitted.
    uint32_t intervalPoaDetectHetBubblesInWindow(
        AnchorWindow& window,
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys,
        const AlignOptions& alignOptions,
        double hetMinVaf,
        uint64_t hetMinSupport,
        bool hetDropHomopolymer,
        bool hetDropRepeat) const;

    // All-windows per-interval POA het detection with global interval load
    // balancing (shasta2 assembleChainsMultithreaded model): every interval of
    // every window is one work unit, flattened into a single list, sorted by
    // descending cost, and run batch=1 across threads. Each interval writes its
    // own fragment (no shared state); a serial-per-window merge then emits the
    // identical AnchorWindow::hetBubbles as the single-window path. This keeps
    // all threads busy even when a few large windows hold most intervals.
    // Returns total bubbles; also reports hetWindows and totalBubbles counts.
    uint32_t intervalPoaDetectHetBubblesAllWindows(
        std::vector<AnchorWindow>& windows,
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys,
        double hetMinVaf,
        uint64_t hetMinSupport,
        bool hetDropHomopolymer,
        bool hetDropRepeat,
        uint64_t threadCount,
        uint64_t& hetWindowsOut,
        uint64_t& totalBubblesOut,
        // Optional per-window skip mask (size == windows.size()). When set and
        // skipWindow[w] is true, window w is left homozygous (no het detection).
        // Used to suppress het calls in highly connected tangle windows.
        const std::vector<bool>* skipWindow = nullptr) const;

    // Test computeAnchorWindowsClean on the longest read, then build
    // a restricted anchor graph from the kept anchors and write GFA.
    void testAnchorWindowsCleanLongestRead(
        uint64_t threadCount,
        uint64_t minInterWindowCoverage);

    // Write AnchorWindowsClean.gfa and .csv from pre-computed windows.
    // Uses cleanHetSnpCount to gate alternate-path output.
    void writeAnchorWindowsCleanGfa(
        const vector<AnchorWindow>& anchorWindows,
        uint64_t minInterWindowCoverage);

    // Original version: claims all unclaimed anchors in the touched range.
    void computeAnchorWindows(
        shared_ptr<Shasta2Anchors> shasta2Anchors,
        shared_ptr<Shasta2Journeys> shasta2Journeys,
        const vector<ReadId>& readIdsSortedByLength,
        vector<AnchorWindow>& anchorWindows,
        uint64_t threadCount);

    // Per-window all-reads multi-segment MSA (abPOA engine). Builds one MSA per
    // anchor window from all oriented reads sharing >=2 of the window's backbone
    // anchors, folding read pieces into an abPOA partial-order graph via
    // abpoa_align_sequence_to_subgraph. Writes
    // testAbpoaMultiSegmentMSA_window<N>.{fasta,gfa}.
    void testAbpoaMultiSegmentMSA(
        const shared_ptr<Shasta2Anchors>& shasta2Anchors,
        const shared_ptr<Shasta2Journeys>& shasta2Journeys,
        vector<AnchorWindow>& anchorWindows,
        uint64_t threadCount,
        double hetMinVaf,
        uint64_t hetMinSupport,
        bool hetDropHomopolymer,
        bool hetDropRepeat);

    // Per-thread state and the thread function for the parallel per-window MSA
    // loop (one window per work item, dynamic load balancing). Each worker owns
    // a distinct AnchorWindow, so windows are processed independently; only the
    // shared counters and cout are serialized (atomics + mutex).
    class AbpoaMultiSegmentMSAData {
    public:
        const shared_ptr<Shasta2Anchors>* shasta2Anchors = nullptr;
        const shared_ptr<Shasta2Journeys>* shasta2Journeys = nullptr;
        vector<AnchorWindow>* anchorWindows = nullptr;
        uint64_t windowEnd = 0;                 // process windows in [0, windowEnd)
        std::atomic<uint64_t> processed{0};
        std::atomic<uint64_t> produced{0};
        // Het-SNP detection tunables (from Assembly.mode3.het* options).
        double hetMinVaf = 0.12;
        uint64_t hetMinSupport = 0;             // 0 = auto-derive from coverage
        bool hetDropHomopolymer = false;        // drop unit-length-1 context
        bool hetDropRepeat = false;             // drop unit-length-2..6 context
    };
    AbpoaMultiSegmentMSAData abpoaMultiSegmentMSAData;
    void testAbpoaMultiSegmentMSAThreadFunction(size_t threadId);

    // isHet detection for one window. Diagnostic output is written to `out`
    // (buffered by the caller and flushed atomically) instead of cout, so it can
    // run in parallel without interleaving.
    bool runOneWindowAbpoaMultiSegmentMSA(
        const shared_ptr<Shasta2Anchors>& shasta2Anchors,
        const shared_ptr<Shasta2Journeys>& shasta2Journeys,
        AnchorWindow& window,
        std::ostream& out,
        double hetMinVaf,
        uint64_t hetMinSupport,
        bool hetDropHomopolymer,
        bool hetDropRepeat);

    // Build a single multi-segment Theseus MSA for one focal read using
    // all its direct overlaps from alignmentTable. Evaluates feasibility
    // of per-read MSA for het-site detection.
    void testDirectOverlapMSA(
        const shared_ptr<Shasta2Anchors>& shasta2Anchors,
        const shared_ptr<Shasta2Journeys>& shasta2Journeys,
        ReadId focalReadId = ReadId(0));

    std::shared_ptr<Shasta2AnchorGraph> shasta2AnchorGraph;
    std::shared_ptr<Shasta2AssemblyGraph> shasta2AssemblyGraph;
    std::map<string, std::shared_ptr<Shasta2AssemblyGraphPostprocessor> > shasta2AssemblyGraphTable;

    // Use a separate mapped-memory prefix for Shasta2 objects so we can
    // keep upstream shasta2 object names ("Journeys", "AnchorGraph", ...)
    // without colliding with other Dinara pipelines.
    MappedMemoryOwner shasta2MappedMemoryOwner() const
    {
        MappedMemoryOwner owner(*this);
        if(!owner.largeDataFileNamePrefix.empty()) {
            owner.largeDataFileNamePrefix += "Shasta2-";
        }
        return owner;
    }

    void exploreShasta2AnchorGraph(const vector<string>& request, ostream& html);
    void exploreShasta2Anchor(const vector<string>& request, ostream& html);
    void exploreShasta2AnchorPair(const vector<string>& request, ostream& html);
    void exploreShasta2AnchorPair2(const vector<string>& request, ostream& html);
    void exploreShasta2Journey(const vector<string>& request, ostream& html);
    void exploreShasta2LocalAnchorGraph(const vector<string>& request, ostream& html);
    void exploreShasta2LocalReadAnchorGraph(const vector<string>& request, ostream& html);
    void exploreShasta2LocalReadGraph(const vector<string>& request, ostream& html);

    void exploreShasta2Segments(const vector<string>& request, ostream& html);
    void exploreShasta2SegmentSequence(const vector<string>& request, ostream& html);
    void exploreShasta2SegmentSteps(const vector<string>& request, ostream& html);
    void exploreShasta2SegmentStepSupport(const vector<string>& request, ostream& html);
    void exploreShasta2SegmentStep(const vector<string>& request, ostream& html);
    void exploreShasta2TangleMatrix(const vector<string>& request, ostream& html);
    void exploreShasta2SegmentPair(const vector<string>& request, ostream& html);
    void exploreShasta2SimilarSequences(const vector<string>& request, ostream& html);

    Shasta2AssemblyGraphPostprocessor& getShasta2AssemblyGraph(
        const string& assemblyStage,
        const Shasta2AssemblyGraphOptions&);

    // Global AnchorGraph (created for all anchors, not per-component).
    shared_ptr<mode3::AnchorGraph> anchorGraph;

    // Verkko-style directed anchor graph (built from BRG anchors).
    shared_ptr<mode3::DirectedAnchorGraph> directedAnchorGraph;

    // If the coverage range for primary marker graph edges is not
    // specified, this uses the disjoint sets histogram to compute reasonable values.
    pair<uint64_t, uint64_t> getPrimaryCoverageRange();

    // Assemble sequence between two primary edges.
    void fillMode3AssemblyPathStep(const vector<string>&, ostream&);

    // Top level function for Mode 3 assembly.
    void mode3Assembly(
        uint64_t threadCount,
        shared_ptr<mode3::Anchors>,
        const Mode3AssemblyOptions&,
        bool debug
    );
    // Same, but use existing Anchors. Python callable.
    void mode3Reassembly(
        uint64_t threadCount,
        const Mode3AssemblyOptions&,
        bool debug
    );

    // Alignment-free version of mode 3 assembly.
    void alignmentFreeAssembly(
        const Mode3AssemblyOptions&,
        const vector<string>& anchorFileAbsolutePaths,
        uint64_t threadCount);

    // Http server functions related to Mode 3 assembly.
    void exploreAnchor(const vector<string>&, ostream&);
    void exploreAnchorPair(const vector<string>&, ostream&);
    void exploreJourney(const vector<string>&, ostream&);
    void exploreReadFollowing(const vector<string>&, ostream&);
    void exploreLocalAssembly(const vector<string>&, ostream&);
    void exploreLocalAnchorGraph(const vector<string>&, ostream&);
    void exploreMode3AssemblyGraph(const vector<string>&, ostream&);
    void exploreSegment(const vector<string>&, ostream&);
    void exploreReadFollowingAssemblyGraph(const vector<string>&, ostream&);

    // Http server functions for BRG-native anchors.
    void exploreBidirectedAnchor(const vector<string>&, ostream&);
    void exploreBidirectedJourney(const vector<string>&, ostream&);
    void exploreBidirectedAnchorGraph(const vector<string>&, ostream&);
    void exploreBidirectedAnchorGraphNode(const vector<string>&, ostream&);
    void exploreBidirectedAnchorGraphPath(const vector<string>&, ostream&);

    // Http server function for the global AnchorGraph.
    void exploreAnchorGraph(const vector<string>&, ostream&);

    // Http server functions for the directed anchor graph (Verkko-style).
    void exploreDirectedAnchorGraph(const vector<string>&, ostream&);
    void exploreDirectedAnchorGraphNode(const vector<string>&, ostream&);
    void exploreDirectedAnchorGraphPath(const vector<string>&, ostream&);


public:
    void test();
};

#if DINARA_TESTING
namespace dinara::testing {
    // Test-only hook: apply the same overlap-support splitting (including bridge removal
    // and quasi-clique peeling) used by marker-vertex anchor decomposition.
    vector<vector<uint32_t>> splitVertexByOverlapSupportForTesting(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage);

    // Test-only hook: run Markov Clustering (MCL) on an (undirected) adjacency list.
    // The adjacency does not need to be symmetric; it is symmetrized internally.
    vector<vector<uint32_t>> mclClusterForTesting(
        const vector<vector<uint32_t>>& adj,
        double inflation,
        uint32_t maxIterations);

    // Test-only hook: split using a provided core mask (see Mode3 vertexSplit.useNonContainedCores).
    vector<vector<uint32_t>> splitVertexByOverlapSupportWithCoreMaskForTesting(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        const vector<uint8_t>& isCore,
        uint32_t coreMinSize,
        uint32_t attachMinSupport,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage);

    // Test-only hook: run the same "auto" logic as the vertex-based anchor splitter:
    // - optional non-contained core splitting + attachment
    // - optional MCL secondary splitting (with suspicious-vertex checks)
    // - quasi-clique peeling
    vector<vector<uint32_t>> autoSplitVertexForTesting(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        const vector<uint8_t>& isCore,
        bool useNonContainedCores,
        uint32_t coreMinSize,
        uint32_t attachMinSupport,
        bool useMclSecondary,
        uint32_t mclMinVertexSize,
        double mclInflation,
        uint32_t mclMaxIterations,
        double suspiciousMaxDensity,
        double suspiciousMaxAverageClustering,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage);

    // Same as autoSplitVertexForTesting, but also returns whether the MCL branch was attempted.
    std::pair<vector<vector<uint32_t>>, bool> autoSplitVertexWithMclTriedFlagForTesting(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        const vector<uint8_t>& isCore,
        bool useNonContainedCores,
        uint32_t coreMinSize,
        uint32_t attachMinSupport,
        bool useMclSecondary,
        uint32_t mclMinVertexSize,
        double mclInflation,
        uint32_t mclMaxIterations,
        double suspiciousMaxDensity,
        double suspiciousMaxAverageClustering,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage);

    // Test-only hook: clique-cover splitter (maximal cliques on core reads + attachment + peeling).
    vector<vector<uint32_t>> splitVertexByCliqueCoverForTesting(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        const vector<uint8_t>& isCore,
        uint32_t attachMinSupport,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage);
}
#endif

#endif
