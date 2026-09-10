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

    // Run hifiasm's aligned overlap path (candidate detection + base-level
    // alignment + filter) rather than the raw pre-alignment candidate set. True
    // by default: it is the only source of hifiasm's per-overlap CIGAR, which
    // spans the full overlap box including the flanks that chaining extends
    // past the outermost anchor -- a span dinara's own A*PA2 layer cannot reach,
    // since it only aligns between consecutive chain anchors. See the long
    // comment at the hifiOpt setup in srcMain/main.cpp for the measured costs
    // (slower overlap detection, a smaller candidate set) and why the CIGARs and
    // the filtering cannot be taken separately.
    bool useHifiasmBaseAlignment = true;

    // Run hifiasm's overlapper with its ONT preset (--ont). True by default.
    //
    // This was previously never set at all: hifiOpt is zero-initialised, so the
    // overlapper ran hifiasm's HiFi preset on ONT reads. One flag controls
    // several things inside hifiasm at once:
    //     max_ov_diff_ec   0.04 -> 0.07   (max overlap error rate)
    //     alignment window  775 -> 375    (WINDOW_HC -> WINDOW_OHC)
    //     chaining band    0.02 -> 0.05   (h_ec_lchain bw_thres)
    //     rl_cut / sc_cut  -1 / 1 -> ONT defaults, is_sc -> 1
    // A 4% error ceiling and a 0.02 band are tight for raw ONT divergence and
    // indel drift, so the HiFi preset was discarding real overlaps.
    //
    // Measured on the GIAB fixture: overlaps 59614 -> 68514, anchors
    // 50630 -> 52550, het sites after filtering 425 -> 662, and the
    // chain-consistency rejection rate -- our measure of marker-graph
    // corruption -- FELL from 6.8% to 4.9%. More data and cleaner data at
    // once, which is what you see when a ceiling was cutting signal rather
    // than noise.
    //
    // Note Align.alignmentWindowLength defaults to 375 to match WINDOW_OHC;
    // that default was only correct once this became true.
    //
    // Set false for HiFi input.
    bool hifiasmIsOnt = true;

    // Windowed overlap acceptance, a port of hifiasm's align_hc_ed_post_extz +
    // pass_qovlp rule (see ProjectedAlignment::alignedWindowFraction). Unlike
    // maxErrorRate, which is a single average over the whole overlap, this is a
    // LOCAL criterion: an overlap that averages well but contains one badly
    // diverging stretch fails it. Both filters are applied; they reject
    // different alignments.
    //
    // Inactive unless computeBaseAlignmentCigar is true (there are no base-level
    // errors to window without it), and disabled outright by a fraction of 0.
    uint32_t alignmentWindowLength = 375;           // hifiasm WINDOW_OHC (ONT)
    double minAlignedWindowFraction = 0.9;          // hifiasm OVERLAP_THRESHOLD_HIFI_FILTER

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

    // Minimum per-edge coverage for the journey anchor graph: an edge between
    // consecutive journey anchors is kept only if at least this many reads
    // traverse the adjacency (posB == posA+1). Default 0 keeps every consecutive
    // pair. Values > 0 threshold adjacency coverage, which can
    // isolate well-supported anchors whose reads reach a neighbor through
    // intermediate anchors.
    uint64_t minJourneyEdgeCoverage;

    // Candidate SNP-site detection from the imported hifiasm CIGARs
    // (Assembler::detectCigarSnpSites). Detection and reporting only -- it
    // creates no anchors and changes nothing downstream; it exists to measure
    // whether the disagreement-fraction signal separates real het sites from
    // sequencing error before anything is built on it. Needs
    // Align.useHifiasmBaseAlignment (the default), since it reads the CIGARs.
    bool detectSnpSites = false;
    // A position is a candidate when at least this many covering partners
    // disagree, and (optionally) the disagreeing share falls in [min, max].
    //
    // The defaults are hifiasm's criterion exactly: >= 2 disagreeing partners
    // and NO fraction constraint (Correct.cpp uses snp_threshold = 1, tested as
    // flag[i] > snp_threshold, with nothing else). Detection is deliberately
    // permissive; deciding whether a minority allele is real belongs to the
    // statistical and context filters downstream, not to a frequency cutoff
    // here.
    //
    // Measured on the GIAB fixture: a [0.2, 0.8] window cut 3662 sites to 518,
    // which does remove real junk (2129 monoallelic sites down to 49) but also
    // discards 1027 of 1493 BIALLELIC sites -- 69% of them. Those sit in the
    // tails: a second allele carried by only a few of ~34 reads, or the owning
    // read itself carrying the rare allele. A frequency cutoff cannot separate
    // "rare real allele" from "recurrent error"; a binomial test against the
    // assumed error rate can, which is why hifiasm has no such cutoff.
    //
    // The fraction bounds are kept as a diagnostic, defaulted off.
    uint64_t snpSiteMinDisagree = 2;
    double snpSiteMinFraction = 0.0;
    double snpSiteMaxFraction = 1.0;
    // Fisher exact p-value below which a site's minority allele is judged
    // strand-biased and dropped. Same test and form as the longcallD-derived
    // one already in the tree (kmFisherExactTwoTail on fwd/rev alt counts
    // against a balanced expectation). A real variant is seen from both
    // directions; one that is not is a strand-specific systematic error.
    double snpSiteStrandBiasPValue = 0.01;
    // hifiasm's biallelic-reduction gates (Correct.cpp, before
    // InsertSNPVector). It does not split a multi-allelic site into several
    // binary ones -- it picks the single strongest alternate and requires the
    // site to be effectively biallelic, dropping it if two alternates tie, if
    // reference+alternate do not account for snpSiteMinPurity of the pileup, or
    // if the alternate is not snpSiteMinAltDominance of all the disagreement.
    // Note hifiasm's own source marks both thresholds "Fix-attention: looks
    // definitely wrong", so they are values to revisit, not gospel.
    double snpSiteMinPurity = 0.95;
    double snpSiteMinAltDominance = 0.70;
    // When two independent anchors land on the same base of one read, shasta2
    // forbids keeping both (positionOffsetAB asserts strictly increasing
    // positions), so one must lose. True keeps the het anchor, false the
    // primary. True by default: preferring primary made het anchors lose EVERY
    // collision they were in (11.1% of all het-anchor occurrences), while the
    // primary that loses instead carries ~2x the coverage and pays
    // proportionally far less.
    bool journeyTiePreferHet = true;
    // Minimum minor-allele fraction for a het site. Measured against the
    // published HG002 mat-vs-pat het track: false positives cluster at VAF ~0.1
    // and real hets at 0.4-0.5, so 0.25 separates them cleanly (precision
    // 75.4% -> 97.5% at a 0.9-point recall cost).
    double snpSiteMinAlleleFraction = 0.25;
    // Sequence-context gates. filterStr (repeat unit 2..6) helps; the
    // homopolymer gate is measured to be net harmful -- it was the largest
    // single cause of missed real variants (49 of 58 filter-rejected truth
    // SNVs) while the VAF floor already removes the false positives it was
    // meant to catch -- so it defaults OFF. See AssemblerCigarSnpSites.cpp.
    bool snpSiteFilterHomopolymer = false;
    bool snpSiteFilterStr = true;
    // hifiasm's `cc`, from the LIVE phasing path (gen_rphase_dp0_single_path):
    //     cc = ((het_cov > 0) ? het_cov : (hom_cov / ploidy));
    //     cc *= cut_rate;  if (cc < cut_bd) cc = cut_bd;
    // with cut_rate = 0.7 and cut_bd = 6 at hifiasm's call site. It gates the
    // MAJOR allele's read count, so it asks whether a site's dominant allele
    // carries the support a real haplotype would. Because each haplotype at a
    // true het site gets about het_cov reads, this adapts to the dataset's own
    // coverage instead of being a guessed constant -- on the GIAB fixture
    // hifiasm's peaks are het 22 / hom 43, giving 0.7 x 22 = 15.
    double snpSiteAlleleCoverageRate = 0.7;
    uint64_t snpSiteAlleleCoverageFloor = 6;
    uint64_t snpSitePloidy = 2;
    // Turn the surviving sites into het anchors: one per allele arm, then the
    // caller rebuilds journeys and the anchor graph from scratch. No surgery on
    // an existing graph. Off by default -- detection alone changes nothing.
    bool createSnpSiteAnchors = false;

    // EXPERIMENTAL: before building anchors, re-check each surviving site with
    // an MSA over the interval between two anchors every member shares. The
    // pairwise CIGARs can agree that two bases differ while disagreeing about
    // whether they are the same locus; an MSA answers that in one frame. Off by
    // default -- it can only remove sites, so it trades recall for precision
    // and the trade has to be measured before it is trusted.
    bool msaVerifySnpSites = false;
    // Half-width of the window each member read contributes to the MSA,
    // centred on its own copy of the site. Against the HG002 truth track on
    // E821, 50 and 150 are indistinguishable (4392 true variants kept either
    // way); 150 is the default for margin against nearby indels.
    uint64_t msaVerifyFlankBases = 150;
    // Fraction of placed members that must agree, applied both to the MSA
    // column their SNP bases land in and to the arm/allele partition. Nearly
    // flat between 0.6 and 0.9 on E821; only 1.0 changes much.
    double msaVerifyMinAgreement = 0.9;


    // Assumed per-read sequencing error rate used by the per-edge MSA het
    // detector's allele significance test: a one-sided binomial test asks
    // whether an allele's read count is explainable as errors misreading the
    // run's dominant allele at this rate, or is significant enough (p <=
    // 0.05) to be a real second haplotype. Same idea and default as the
    // myloasm SNPmer caller's minor-allele test. Lower values make the test
    // MORE PERMISSIVE (a given minor-allele count looks more surprising
    // against a lower assumed noise floor, so fewer reads are needed to pass);
    // higher values make it STRICTER (more minor-allele reads are needed to
    // stand out above a higher assumed noise floor).
    double hetErrorRate;

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
