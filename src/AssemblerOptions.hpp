#ifndef DINARA_ASSEMBLER_OPTIONS_HPP
#define DINARA_ASSEMBLER_OPTIONS_HPP


/*******************************************************************************

CLASSES DESCRIBING ASSEMBLER OPTIONS

There are two types of options:
- Options that can be used both on the command line and in a configuration file
  ("configurable options").
- Options that can only be used in a configuration file
  ("non-configurable options").

Configuration files are divided in sections formatted like this:

[SectionName]
optionName = optionValue

The command line syntax corresponding to the above is
--SectionName.optionName optionValue

Each SectionName corresponds to a class defined below. For example,
section [Align] corresonds to class AlignOptions.

If the option is a Boolean switch, use True or False as the optionValue.



ADDING A NEW CONFIGURABLE OPTION

1. Add the option to the class corresponding to the desired section.
2. Modify the write function to that class to also write the newly added option.
3. Modify AssemblerOptions::addConfigurableOptions to reflect the new option,
   making sure to include a default value and at least a minimal help message.
4. Document the option in dinara/docs/CommandLineOptions.html.
5. If the option requires validation add it at the appropriate place in
   dinara/srcMain/main.cpp.

For options not ready for end users, it is fine in steps 3 4 5 to use a
comment just saying "Experimental - leave at default value"
or something to that effect.



ADDING A NEW NON-CONFIGURABLE OPTION

1. Add the option to class CommandLineOnlyOptions.
2. Modify AssemblerOptions::addCommandLineOnlyOptions to reflect the new option,
   making sure to include a default value and at least a minimal help message.
3. Document the option in dinara/docs/CommandLineOptions.html.
4. If the option requires validation add it at the appropriate place in
   dinara/srcMain/main.cpp.

For options not ready for end users, it is fine in steps 3 4 5 to use a
comment just saying "Experimental - leave at default value"
or something to that effect.

*******************************************************************************/

// Boost libraries.
#include <boost/program_options.hpp>

// Standard library.
#include "iostream.hpp"
#include "string.hpp"
#include "vector.hpp"

namespace dinara {
    class AlignOptions;
    class Align6Options;
    class AssemblerOptions;
    class AssemblyOptions;
    class CommandLineOnlyOptions;
    class KmersOptions;
    class MarkerGraphOptions;
    class MinHashOptions;
    class Mode3AssemblyOptions;
    class OverlapCandidatesOptions;
    class PalindromicReadOptions;
    class ReadsOptions;
    class ReadGraphOptions;
    class VariantClusteringOptions;

    // Function to convert a bool to True or False for better
    // compatibility with Python scripts.
    string convertBoolToPythonString(bool);
}



// Options only allowed on the command line and not in the configuration file.
class dinara::CommandLineOnlyOptions {
public:
    string configName;
    vector <string> inputFileNames;
    vector <string> anchorFileNames;
    string assemblyDirectory;
    string shasta2OutputDirectory;
    string command;
    string memoryMode;
    string memoryBacking;
    uint32_t threadCount;
    bool suppressStdoutLog;
    string exploreAccess;
    uint16_t port;
    string alignmentsPafFile;
    bool saveBinaryData;
};


class dinara::PalindromicReadOptions {
public:
    bool skipFlagging;
    double alignedFractionThreshold;
    double maxErrorRate;
    void write(ostream&) const;
};


// Options in the [Reads] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Reads.".
class dinara::ReadsOptions {
public:
    uint64_t representation;    // 0 = Raw, 1=RLE
    int minReadLength;
    bool noCache;
    string desiredCoverageString;
    uint64_t desiredCoverage;

    // String to control handling of duplicate reads.
    // Can be one of:
    // useAllCopies
    // useOneCopy
    // useNone
    // forbid
    // See ReadFlags.hpp for the meaning of each option.
    string handleDuplicates;

    PalindromicReadOptions palindromicReads;

    void write(ostream&) const;

    void parseDesiredCoverageString();
};



// Options in the [Kmers] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Kmers.".
class dinara::KmersOptions {
public:
    int generationMethod;
    int k;
    double probability;
    double enrichmentThreshold;
    uint64_t distanceThreshold;
    string file;
    string globalFrequencyOverrideDirectory;

    // Minimum fraction of read length that must be covered by the marker span
    // (lastMarkerPos + k - firstMarkerPos). Reads below this are discarded.
    // 0 disables the filter.
    double minMarkerSpanFraction;

    void write(ostream&) const;
};



// Options in the [MinHash] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "MinHash.".
class dinara::MinHashOptions {
public:
    int version;
    int m;
    double hashFraction;
    int minHashIterationCount;
    double alignmentCandidatesPerRead;
    int minBucketSize;
    int maxBucketSize;
    int minFrequency;
    bool allPairs;
    void write(ostream&) const;
};


// Options in the [OverlapCandidates] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "OverlapCandidates.".
class dinara::OverlapCandidatesOptions {
public:
    string method;              // "MinHash" or "InvertedIndex".
    double driftRateTolerance;  // Drift rate tolerance for chaining. Hifiasm: 0.05 for ONT, 0.02 for HiFi.
    int minChainMarkerCount = 2;    // Minimum marker count required for a chained overlap candidate. Hifiasm: min_lc_cnt = 2 (inter.cpp:496).
    uint32_t minOverlapLength = 1000; // Minimum overlap span (bases) for a candidate to be kept (min of query/target spans); 0 disables.
    uint32_t maxEndFuzz = 0;           // Deprecated and IGNORED: dovetail/internal filtering moved to deleteInternalOverlaps (ma_hit2arc). Kept so old command lines still parse.
    uint32_t maxChainingFreq = 1000;  // Skip kmers with frequency above this during chaining (markers still kept for journeys).

    // When > 0, only chain pairs where at least one read has
    // readId < referenceReadCount. Skips read-vs-read pairs.
    uint64_t referenceReadCount = 0;

    void write(ostream&) const;
};



class dinara::Align6Options {
public:
    uint64_t maxLocalFrequency;
    uint64_t minGlobalFrequency;
    uint64_t maxGlobalFrequency;
    double maxGlobalFrequencyMultiplier;
    uint64_t minLowFrequencyCount;
    double driftRateTolerance;
    uint64_t maxInBandCount;
    double maxInBandRatio;

    void write(ostream&) const;
};



// Options in the [Align] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Align.".
class dinara::AlignOptions {
public:
    int alignMethod;
    int maxSkip;
    int maxDrift;
    int maxTrim;
    int maxMarkerFrequency;
    int minAlignedMarkerCount;
    double minAlignedFraction;
    int matchScore;
    int mismatchScore;
    int gapScore;
    double downsamplingFactor;
    int bandExtend;
    int maxBand;
    int sameChannelReadAlignmentSuppressDeltaThreshold;
    bool suppressContainments;
    double maxErrorRate;

    // When true, compute a base-level alignment (per-segment A*PA2 CIGAR) for
    // each overlap and store it, deriving base statistics (edit distance,
    // mismatches, indels, error rate). When false (default), skip base
    // alignment entirely: only the marker-ordinal chain and its span are kept.
    // Downstream marker/read-graph construction reads the ordinal chain, not the
    // base CIGAR, so the CIGAR is wasted work unless a base-level consumer (e.g.
    // CIGAR-based phasing/MSA) is enabled. The A*PA2 implementation is retained;
    // this flag only gates whether it runs.
    bool computeBaseAlignmentCigar = false;

    // Overlap/base DP scoring parameters (used to compute AlignmentInfo::dpScore from a base-level CIGAR).
    // These should be configured to match hifiasm's overlap-alignment scoring model.
    // Current hifiasm overlap scoring is single-affine: gapCost(k) = O1 + k*E1.
    // O2/E2 are retained for backward-compatible configuration plumbing but are ignored.
    int64_t overlapDpMatchScore;
    int64_t overlapDpMismatchScore;
    int64_t overlapDpGapOpen1;
    int64_t overlapDpGapExtend1;
    int64_t overlapDpGapOpen2;
    int64_t overlapDpGapExtend2;

    // Align4.
    uint64_t align4DeltaX;
    uint64_t align4DeltaY;
    uint64_t align4MinEntryCountPerCell;
    uint64_t align4MaxDistanceFromBoundary;

    // Align5.
    double align5DriftRateTolerance;
    uint64_t align5MinBandExtend;

    // Align6.
    Align6Options align6Options;

    void write(ostream&) const;
};



// Options in the [ReadGraph] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "ReadGraph.".
class dinara::ReadGraphOptions {
public:
    int creationMethod;
    int maxAlignmentCount;
    bool preferAlignedFraction;
    int maxChimericReadDistance;
    uint64_t strandSeparationMethod;
    int crossStrandMaxDistance;
    bool removeConflicts;
    double markerCountPercentile;
    double alignedFractionPercentile;
    double maxSkipPercentile;
    double maxDriftPercentile;
    double maxTrimPercentile;
    bool flagInconsistentAlignments;
    uint64_t flagInconsistentAlignmentsTriangleErrorThreshold;
    uint64_t flagInconsistentAlignmentsLeastSquareErrorThreshold;
    uint64_t flagInconsistentAlignmentsLeastSquareMaxDistance;
    // New readGraph4withStrandSeparation options
    double epsilon;
    double delta;
    double WThreshold;
    double WThresholdForBreaks;
    
    // Cluster graph options.
    uint64_t clusterGraphMinEdgeCoverage;
    
    // Filtering options for phased chains/sites
    int minMultiNodeChainSupport;
    int minIsolatedSiteSupport;

    // If set, hifiasm-style "secondary overlap" filtering requires non-redundancy
    // on both reads (query and target), not just on the query read.
    bool filterSecondaryRequireNonRedundantOnBothReads = true;

    void write(ostream& ) const;
};



// Options in the [MarkerGraph] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "MarkerGraph.".
class dinara::MarkerGraphOptions {
public:
    int minCoverage;
    int maxCoverage;
    int minCoveragePerStrand;
    uint64_t minEdgeCoverage;
    uint64_t minEdgeCoveragePerStrand;
    bool allowDuplicateMarkers;
    int lowCoverageThreshold;
    int highCoverageThreshold;
    int maxDistance;
    int edgeMarkerSkipThreshold;
    double peakFinderMinAreaFraction;
    uint64_t peakFinderAreaStartIndex;
    bool alwaysSave;

    // Optional diagnostics.
    bool writeVertexCoverageHistogram;
    string vertexCoverageHistogramFileName;
    bool vertexCoverageHistogramCanonicalOnly;

    void write(ostream&) const;
};
// Options for variant clustering.
class dinara::VariantClusteringOptions {
public:
    uint64_t minOccurrences;
    uint64_t minSeparation;
    void write(ostream&) const;
};




// Assembly options that are specific to Mode 3 assembly.
// See source code in the mode3 namespace
// (source files with a mode3- prefix) for more information
class dinara::Mode3AssemblyOptions {
public:

    string anchorCreationMethod;

    uint64_t minAnchorCoverage;
    uint64_t maxAnchorCoverage;
    double minAnchorCoverageMultiplier;
    double maxAnchorCoverageMultiplier;

    // Options used by anchor creation methods that split marker graph vertices
    // using readGraph overlap support (for example: FromMarkerGraphVerticesSplitUsingReadGraph).
    class VertexSplitOptions {
    public:
        // If true, run Markov Clustering (MCL) as a secondary splitter for "suspicious"
        // vertices that remain a single cluster after the default bridge-removal + peeling logic.
        bool useMclSecondary;

        // Only consider MCL for vertices with at least this many oriented reads.
        uint32_t mclMinVertexSize;

        // MCL inflation parameter (controls cluster granularity). Typical values: 1.4-2.5.
        double mclInflation;

        // Maximum MCL iterations.
        uint32_t mclMaxIterations;

        // Trigger MCL only when the overlap-support graph looks "non-clique-like".
        // MCL is attempted only if both conditions hold:
        // - density <= suspiciousMaxDensity
        // - averageClustering <= suspiciousMaxAverageClustering
        double suspiciousMaxDensity;
        double suspiciousMaxAverageClustering;

        // If true, attempt to split vertices using only non-contained reads
        // as "core" evidence, then attach contained reads to exactly one core cluster.
        // This helps when contained reads act as bridges between unrelated regions/strands.
        bool useNonContainedCores;

        // Minimum number of core (non-contained) oriented reads required to attempt core-based splitting.
        uint32_t coreMinSize;

        // Minimum number of edges required to attach a non-core (typically contained) read to a cluster.
        uint32_t attachMinSupport;

        void write(ostream&) const;
    };
    VertexSplitOptions vertexSplitOptions;

    // Options used to clean up the PrimaryGraph.
    class PrimaryGraphOptions {
    public:

        // Parameter to control removal of weak edges.
        double maxLoss;

        // Parameters to control removal of cross edges.
        uint64_t crossEdgesLowCoverageThreshold;
        uint64_t crossEdgesHighCoverageThreshold;

        void write(ostream&) const;
    };
    PrimaryGraphOptions primaryGraphOptions;



    class AssemblyGraphOptions {
    public:

        // Detangle tolerances.
        uint64_t detangleToleranceLow;
        uint64_t detangleToleranceHigh;

        bool suppressBubbleCleanup;

        // Bayesian model.
        double epsilon;
        double minLogP;

        // Other thresholds used by the mode3::AssemblyGraph
        uint64_t longBubbleThreshold;
        double phaseErrorThreshold;
        double bubbleErrorThreshold;
        uint64_t bubbleCleanupMaxOffset;
        uint64_t chainTerminalCommonThreshold;
        uint64_t superbubbleLengthThreshold1;
        uint64_t superbubbleLengthThreshold2;
        uint64_t superbubbleLengthThreshold3;
        uint64_t superbubbleLengthThreshold4;
        uint64_t pruneLength;

        void write(ostream&) const;
    };
    AssemblyGraphOptions assemblyGraphOptions;



    // Minimum shared read count for an inter-window edge to be created.
    // Candidates with fewer shared reads are discarded.
    uint64_t minInterWindowCoverage;

    // Minimum anchor pair coverage (common reads) for an inter-window edge.
    // Edges whose anchor pair has fewer common reads are discarded.
    uint64_t minInterWindowEdgeCoverage;

    // Minimum common reads between consecutive backbone anchors.
    // Backbone journeys are filtered to keep the longest subsequence
    // where every consecutive pair meets this threshold.
    uint64_t minCommonForBackbone;

    // Maximum number of positions to look back when filtering
    // backbone journeys for well-supported consecutive pairs.
    uint64_t maxSkipForBackbone;

    // Minimum anchor-pair coverage (common two-sided reads) for an
    // anchor-graph edge to be considered for per-edge MSA het detection
    // (experimental, DINARA_HET_ON_GRAPH=1). Edges below this are skipped.
    uint64_t minCommonForHet;

    // Minimum base span (first anchor to last anchor) for a read's
    // journey to be accepted as a window backbone.
    uint64_t minWindowBaseSpan;

    // --- Per-window abPOA het-SNP detection tunables ---
    // Minimum variant allele fraction for an alt allele to be accepted.
    double hetMinVaf;

    // Minimum per-allele read support. If 0, it is auto-derived from the
    // k-mer coverage histogram (peak/2 * 0.7, floored at 6), matching the
    // hifiasm het-site rule.
    uint64_t hetMinSupport;

    // If true, drop SNPs whose backbone context is a homopolymer run
    // (repeat unit length 1). Default false.
    bool hetDropHomopolymer;

    // If true, drop SNPs whose backbone context is a short-tandem-repeat run
    // (repeat unit length 2..6). Default false.
    // Both default false: the flank-linearity test already requires a clean
    // (bubble-free) homozygous base on each side, so repeat-context SNPs that
    // pass it are real het sites; dropping them discarded far more true SNPs
    // than it kept.
    bool hetDropRepeat;

    // Skip het-anchor detection in windows that are highly connected on both
    // sides -- i.e. windows with at least this many distinct incoming AND this
    // many distinct outgoing inter-window neighbors. Such windows sit at
    // tangles/repeats where per-window het calls are unreliable. A window is
    // skipped only when BOTH its in-degree and out-degree meet the threshold.
    // 0 disables the gate (default), so behavior is unchanged unless set.
    uint64_t hetMaxWindowInDegree;
    uint64_t hetMaxWindowOutDegree;

    // Options used by class mode3::LocalAssembly
    class LocalAssemblyOptions {
    public:

        // The estimated offset gets extended by this ratio to
        // decide how much to extend reads that only appear in edgeIdA or edgeIdB.
        double estimatedOffsetRatio;

        // Vertex sampling rate, used to set minVertexCoverage.
        // Only used if minVertexCoverage is 0 on input to
        // mode3::LocalAssembly constructor.
        double vertexSamplingRate;

        // Alignment parameters.
        int64_t matchScore;
        int64_t mismatchScore;
        int64_t gapScore;

        // Number of bases (not markers) that can be skipped by an alignment.
        uint64_t maxSkipBases;

        // The maximum tolerated length drift of each read.
        // Used to compute the band for banded alignments.
        double maxDrift;

        // Minimum half band, in markers.
        uint64_t minHalfBand;

        // Minimum ration of scorew to best possible score for
        // an alignment to be used.
        double minScoreRatio;

        // The maximum length of an MSA alignment we are willing to compute.
        uint64_t maxMsaLength;

        void write(ostream&) const;
    };
    LocalAssemblyOptions localAssemblyOptions;

    void write(ostream&) const;
};



// Options in the [Assembly] section of the configuration file.
// Can also be entered on the command line with option names
// beginning with "Assembly.".
class dinara::AssemblyOptions {
public:
    int markerGraphEdgeLengthThresholdForConsensus;
    string consensusCallerString;
    string consensusCaller;
    bool storeCoverageData;
    int storeCoverageDataCsvLengthThreshold;
    bool writeReadsByAssembledSegment;
    uint64_t pruneLength;

    // Options that control detangling.
    int detangleMethod;
    uint64_t detangleDiagonalReadCountMin;
    uint64_t detangleOffDiagonalReadCountMax;
    double detangleOffDiagonalRatio;

    // Options that control iterative assembly.
    bool iterative;
    uint64_t iterativeIterationCount;
    int64_t iterativePseudoPathAlignMatchScore;
    int64_t iterativePseudoPathAlignMismatchScore;
    int64_t iterativePseudoPathAlignGapScore;
    double iterativeMismatchSquareFactor;
    double iterativeMinScore;
    uint64_t iterativeMaxAlignmentCount;
    uint64_t iterativeBridgeRemovalIterationCount;
    uint64_t iterativeBridgeRemovalMaxDistance;

    // Mode 3 assembly options.
    Mode3AssemblyOptions mode3Options;

    void write(ostream&) const;

    // If a relative path is provided for a Bayesian consensus caller
    // replace it with its absolute path.
    void parseConsensusCallerString();
};



class dinara::AssemblerOptions {
public:

    // Object containing the options.
    CommandLineOnlyOptions commandLineOnlyOptions;
    ReadsOptions readsOptions;
    KmersOptions kmersOptions;
    MinHashOptions minHashOptions;
    OverlapCandidatesOptions overlapCandidatesOptions;
    AlignOptions alignOptions;
    ReadGraphOptions readGraphOptions;
    VariantClusteringOptions variantClusteringOptions;
    MarkerGraphOptions markerGraphOptions;
    AssemblyOptions assemblyOptions;

    // Constructor from a command line.
    // If the command line includes a --config option,
    // the specified built-in configuration or configuration file
    // is used to fill the AssemblyOptions,
    // but values specified on the command line take precedence.
    AssemblerOptions(int argumentCount, const char** arguments);

    // Constructor from a configuration file.
    // This only fills in the configurable options specified in
    // the given configuration file. Command line only options
    // are left at their defaults.
    AssemblerOptions(const string& fileName);

    // Add configurable options to the Boost option description object.
    void addCommandLineOnlyOptions();
    void addConfigurableOptions();

    // Write the options as a config file.
    void write(ostream&) const;

    // Boost program_options library objects.
    boost::program_options::options_description commandLineOnlyOptionsDescription;
    boost::program_options::options_description configurableOptionsDescription;
    boost::program_options::options_description allOptionsDescription;

    // This one is the same as allOptionsDescription, with
    // "--invalidOption" added to capture invalid positional options.
    vector<string> invalidPositionalOptions;
    boost::program_options::options_description allOptionsIncludingInvalidDescription;

};

#endif
