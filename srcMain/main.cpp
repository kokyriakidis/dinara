// The static executable provides
// basic functionality and reduced performance.
// For full functionality use the shared library built
// under directory src.

// Dinara.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "buildId.hpp"
#if DINARA_ENABLE_VARIANT_CLUSTERING
#include "ClusterGraph.hpp"
#endif
#include "filesystem.hpp"
#include "mode3-Anchor.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "Tee.hpp"
#include "timestamp.hpp"
#include "platformDependent.hpp"


using namespace dinara;

// Boost libraries.
#include <boost/program_options.hpp>
#include  <boost/chrono/process_cpu_clocks.hpp>

//  Linux.
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

// Standard library.
#include "chrono.hpp"
#include <filesystem>
#include "iostream.hpp"

#include "stdexcept.hpp"


// Shasta 2 Integration
#include "AssemblerShasta2Anchors.hpp"
#include <thread>
#include <vector>

// Minimizer filtering
#include "MinimizerChecker.hpp"
#include "MarkerKmers.hpp"
#include "KmerCounter.hpp"



namespace dinara {
    namespace main {

        void main(int argumentCount, const char** arguments);

        void setupRunDirectory(
            const string& memoryMode,
            const string& memoryBacking,
            size_t& pageSize,
            string& dataDirectory
            );

        void setupHugePages();
        void segmentFaultHandler(int);
        void setupSegmentFaultHandler();

        // Functions that implement --command keywords
        void assemble(const AssemblerOptions&, int argumentCount, const char** arguments);
        void saveBinaryData(const AssemblerOptions&);
        void cleanupBinaryData(const AssemblerOptions&);
        void explore(const AssemblerOptions&);
        void listCommands();

        const std::set<string> commands = {
            "assemble",
            "saveBinaryData",
            "cleanupBinaryData",
            "explore",
            "listCommands"};



        void assemble(
            Assembler&,
            const AssemblerOptions&,
            vector<string> inputNames);

    }

    // This is used to duplicate cout output to stdout.log.
    Tee tee;
    ofstream dinaraLog;
}




int main(int argumentCount, const char** arguments)
{
    try {

        dinara::main::main(argumentCount, arguments);

    } catch(const boost::program_options::error_with_option_name& e) {
        cout << "Invalid option: " << e.what() << endl;
        return 1;
    } catch (const runtime_error& e) {
        cout << timestamp << e.what() << endl;
        return 2;
    } catch (const std::bad_alloc& e) {
        cout << timestamp << e.what() << endl;
        cout << "Memory allocation failure." << endl;
        cout << "This assembly requires more memory than available." << endl;
        cout << "Rerun on a larger machine." << endl;
        return 2;
    } catch (const exception& e) {
        cout << timestamp << e.what() << endl;
        return 3;
    } catch (...) {
        cout << timestamp << "Terminated after catching a non-standard exception." << endl;
        return 4;
    }
    return 0;
}



void dinara::main::segmentFaultHandler(int)
{
    char message[] = "\nA segment fault occurred. Please report it by filing an "
        "issue on the Dinara repository and attaching the entire log output. "
        "To file an issue, point your browser to https://github.com/kokyriakidis/dinara/issues\n";
    ::write(fileno(stderr), message, sizeof(message));
    ::_exit(1);
}

void dinara::main::setupSegmentFaultHandler()
{
    struct sigaction action;
    ::memset(&action, 0, sizeof(action));
    action.sa_handler = &segmentFaultHandler;
    sigaction(SIGSEGV, &action, 0);
}


void dinara::main::main(int argumentCount, const char** arguments)
{
    setupSegmentFaultHandler();

    // Parse command line options and the configuration file, if one was specified.
    AssemblerOptions assemblerOptions(argumentCount, arguments);

    // Check that we have a valid command.
    auto it = commands.find(assemblerOptions.commandLineOnlyOptions.command);
    if(it ==commands.end()) {
        const string message = "Invalid command " + assemblerOptions.commandLineOnlyOptions.command;
        listCommands();
        throw runtime_error(message);
    }

    // Execute the requested command.
    if(assemblerOptions.commandLineOnlyOptions.command == "assemble") {
        assemble(assemblerOptions, argumentCount, arguments);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "cleanupBinaryData") {
        cleanupBinaryData(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "saveBinaryData") {
        saveBinaryData(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "explore") {
        explore(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "listCommands") {
        listCommands();
        return;
    }

    // We already checked for a valid command above, so if we get here
    // the above logic is missing code for one of the valid commands.
    DINARA_ASSERT(0);

}




// Implementation of --command assemble.
void dinara::main::assemble(
    const AssemblerOptions& assemblerOptions,
    int argumentCount, const char** arguments)
{
    DINARA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "assemble");


    // Various checks for option validity.

    if(assemblerOptions.kmersOptions.k > 62 or assemblerOptions.kmersOptions.k < 6) {
        throw runtime_error("Invalid value specified for --Kmers.k. Must be between 6 and 62.");
    }

    if((assemblerOptions.kmersOptions.k % 2) == 1) {
        throw runtime_error("Invalid value specified for --Kmers.k. Must be even.");
    }

    // Check that we have at least one input file.
    if(assemblerOptions.commandLineOnlyOptions.inputFileNames.empty()) {
        throw runtime_error("Specify at least one input file "
            "using command line option --input.");
    }

    if( assemblerOptions.alignOptions.alignMethod <  0 or
        assemblerOptions.alignOptions.alignMethod == 2 or
        assemblerOptions.alignOptions.alignMethod >  6) {
        throw runtime_error("Align method " + to_string(assemblerOptions.alignOptions.alignMethod) +
            " is not valid. Valid options are 0 through 6 except 2.");
    }

    // Find absolute paths of the input files.
    // We will use them below after changing directory to the output directory.
    vector<string> inputFileAbsolutePaths;
    for(const string& inputFileName: assemblerOptions.commandLineOnlyOptions.inputFileNames) {
        if(!std::filesystem::exists(inputFileName)) {
            throw runtime_error("Input file not found: " + inputFileName);
        }
        if(!std::filesystem::is_regular_file(inputFileName)) {
            throw runtime_error("Input file is not a regular file: " + inputFileName);
        }
        inputFileAbsolutePaths.push_back(filesystem::getAbsolutePath(inputFileName));
    }

    // Create the assembly directory. If it exists, stop.
    bool exists = std::filesystem::exists(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
    if (exists) {
        throw runtime_error(
            assemblerOptions.commandLineOnlyOptions.assemblyDirectory +
            " already exists. Remove it first \n"
            "or use --assemblyDirectory to specify a different assembly directory."
        );
    } else {
        DINARA_ASSERT(std::filesystem::create_directory(assemblerOptions.commandLineOnlyOptions.assemblyDirectory));
    }

    // Make the assembly directory current.
    std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);

    // Open the performance log.
    openPerformanceLog("performance.log");
    performanceLog << timestamp << "Assembly begins." << endl;

    // Open stdout.log and "tee" (duplicate) stdout to it.
    if(not assemblerOptions.commandLineOnlyOptions.suppressStdoutLog) {
        dinaraLog.open("stdout.log");
        tee.duplicate(cout, dinaraLog);
    }

    // Echo out the command line options.
    cout << timestamp << "Assembly begins.\nCommand line:" << endl;
    for(int i=0; i<argumentCount; i++) {
        cout << arguments[i] << " ";
    }
    cout << endl;

    // Set up the run directory as required by the memoryMode and memoryBacking options.
    size_t pageSize = 0;
    string dataDirectory;
    setupRunDirectory(
        assemblerOptions.commandLineOnlyOptions.memoryMode,
        assemblerOptions.commandLineOnlyOptions.memoryBacking,
        pageSize,
        dataDirectory);

    // Write out the option in effect to dinara.conf.
    {
        ofstream configurationFile("dinara.conf");
        assemblerOptions.write(configurationFile);
    }
    cout << "For options in use for this assembly, see dinara.conf in the assembly directory." << endl;

    // Create the Assembler.
    Assembler assembler(dataDirectory, true, assemblerOptions.readsOptions.representation, pageSize);
    assembler.assemblerInfo->readGraphCreationMethod = assemblerOptions.readGraphOptions.creationMethod;
    assembler.assemblerInfo->assemblyMode = assemblerOptions.assemblyOptions.mode;
    assembler.minMultiNodeChainSupport = assemblerOptions.readGraphOptions.minMultiNodeChainSupport;
    assembler.minIsolatedSiteSupport = assemblerOptions.readGraphOptions.minIsolatedSiteSupport;
    assembler.assemblerInfo->variantClusteringMinOccurrences = assemblerOptions.variantClusteringOptions.minOccurrences;
    assembler.assemblerInfo->variantClusteringMinSeparation = assemblerOptions.variantClusteringOptions.minSeparation;

    // Run the assembly.
    assemble(assembler, assemblerOptions, inputFileAbsolutePaths);

    cout << timestamp << "Assembly ends." << endl;
    performanceLog << timestamp << "Assembly ends." << endl;
}



// Set up the run directory as required by the memoryMode and memoryBacking options.
void dinara::main::setupRunDirectory(
    const string& memoryMode,
    const string& memoryBacking,
    size_t& pageSize,
    string& dataDirectory
    )
{

    if(memoryMode == "anonymous") {

        if(memoryBacking == "disk") {

            // This combination is meaningless.
            throw runtime_error("\"--memoryMode anonymous\" is not allowed in combination "
                "with \"--memoryBacking disk\".");

        } else if(memoryBacking == "4K") {

            // Anonymous memory on 4KB pages.
            // This combination is the default.
            // It does not require root privilege.
            dataDirectory = "";
            pageSize = 4096;

        } else if(memoryBacking == "2M") {

            // Anonymous memory on 2MB pages.
            // This may require root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            // Root privilege is not required if 2M pages have already
            // been set up as required.
#ifdef __APPLE__
            throw runtime_error("Option --memoryBacking 2M is not supported on macOS.");
#else
            setupHugePages();
            pageSize = 2 * 1024 * 1024;
#endif

        } else {
            throw runtime_error("Invalid value specified for --memoryBacking: " + memoryBacking +
                "\nValid values are: disk, 4K, 2M.");
        }

    } else if(memoryMode == "filesystem") {

        if(memoryBacking == "disk") {

            // Binary files on disk.
            // This does not require root privilege.
            DINARA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 4096;

        } else if(memoryBacking == "4K") {

            // Binary files on the tmpfs filesystem
            // (filesystem in memory backed by 4K pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
#ifdef __APPLE__
            throw runtime_error("Option --memoryMode filesystem --memoryBacking 4K is not supported on macOS.");
#else
            DINARA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 4096;
            const string command = "sudo mount -t tmpfs -o size=0 tmpfs Data";
            const int errorCode = ::system(command.c_str());
            if(errorCode != 0) {
                throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                    " running command: " + command);
            }
#endif

        } else if(memoryBacking == "2M") {

            // Binary files on the hugetlbfs filesystem
            // (filesystem in memory backed by 2M pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
#ifdef __APPLE__
            throw runtime_error("Option --memoryMode filesystem --memoryBacking 2M is not supported on macOS.");
#else
            setupHugePages();
            DINARA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 2 * 1024 * 1024;
            const uid_t userId = ::getuid();
            const gid_t groupId = ::getgid();
            const string command = "sudo mount -t hugetlbfs -o pagesize=2M"
                ",uid=" + to_string(userId) +
                ",gid=" + to_string(groupId) +
                " none Data";
            const int errorCode = ::system(command.c_str());
            if(errorCode != 0) {
                throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                    " running command: " + command);
            }
#endif

        } else {
            throw runtime_error("Invalid value specified for --memoryBacking: " + memoryBacking +
                "\nValid values are: disk, 4K, 2M.");
        }

    } else {
        throw runtime_error("Invalid value specified for --memoryMode: " + memoryMode +
            "\nValid values are: anonymous, filesystem.");
    }
}



// This runs the entire assembly, under the following assumptions:
// - The current directory is the run directory.
// - The Data directory has already been created and set up, if necessary.
// - The input file names are either absolute,
//   or relative to the run directory, which is the current directory.
void dinara::main::assemble(
    Assembler& assembler,
    const AssemblerOptions& assemblerOptions,
    vector<string> inputFileNames)
{
    const auto steadyClock0 = std::chrono::steady_clock::now();
    const auto userClock0 = boost::chrono::process_user_cpu_clock::now();
    const auto systemClock0 = boost::chrono::process_system_cpu_clock::now();

    // Adjust the number of threads, if necessary.
    uint64_t threadCount = assemblerOptions.commandLineOnlyOptions.threadCount;
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "This assembly will use " << threadCount << " threads." << endl;

    // If --saveBinaryData was requested,
    // create the directory where binary data will be saved.
    if (assemblerOptions.commandLineOnlyOptions.saveBinaryData) {
        assembler.createSaveBinaryDataDirectory(assemblerOptions.commandLineOnlyOptions.memoryMode);
    }


    // Add reads from the specified input files.
    performanceLog << timestamp << "Begin loading reads from " << inputFileNames.size() << " files." << endl;
    const auto t0 = steady_clock::now();
    for(const string& inputFileName: inputFileNames) {

        assembler.addReads(
            inputFileName,
            assemblerOptions.readsOptions.minReadLength,
            assemblerOptions.readsOptions.noCache,
            threadCount);
    }

    if(assembler.getReads().readCount() == 0) {
        throw runtime_error("There are no input reads.");
    }



    // If requested, increase the read length cutoff
    // to reduce coverage to the specified amount.
    if (assemblerOptions.readsOptions.desiredCoverage > 0) {
        // Write out the read length histogram using provided minReadLength.
        assembler.histogramReadLength("ExtendedReadLengthHistogram.csv");

        const auto newMinReadLength = assembler.adjustCoverageAndGetNewMinReadLength(
            assemblerOptions.readsOptions.desiredCoverage);

        const auto oldMinReadLength = uint64_t(assemblerOptions.readsOptions.minReadLength);

        if (newMinReadLength == 0ULL) {
            throw runtime_error(
                "With Reads.minReadLength " +
                to_string(assemblerOptions.readsOptions.minReadLength) +
                ", total available coverage is " +
                to_string(assembler.getReads().getTotalBaseCount()) +
                ", less than desired coverage " +
                to_string(assemblerOptions.readsOptions.desiredCoverage) +
                ". Try reducing Reads.minReadLength if appropriate or get more coverage."
            );
        }

        // Adjusting coverage should only ever reduce coverage if necessary.
        DINARA_ASSERT(newMinReadLength >= oldMinReadLength);
    }

    assembler.computeReadIdsSortedByName();
    assembler.histogramReadLength("ReadLengthHistogram.csv");

    const auto t1 = steady_clock::now();
    performanceLog << timestamp << "Done loading reads from " << inputFileNames.size() << " files." << endl;
    performanceLog << "Read loading took " << seconds(t1-t0) << "s." << endl;

    // Find duplicate reads and handle them according to the setting
    // of --Reads.handleDuplicates. The default option is "useOneCopy".
    assembler.findDuplicateReads(assemblerOptions.readsOptions.handleDuplicates);

    // Find markers using either SIMD closed syncmers or the default k-mer based method.
    if(assemblerOptions.kmersOptions.useSimdClosedSyncmers) {
        // Use SIMD-accelerated closed syncmers for initial marker generation (no filtering).
        // This generates a superset of the markers we eventually want.
        assembler.findMarkersSimdClosedSyncmers(
            threadCount,
            assemblerOptions.kmersOptions.k,
            assemblerOptions.kmersOptions.syncmerS);

        // Compute histogram using the pre-calculated KmerIds.
        // This avoids accessing the Reads data structure (Cache Misses).
        assembler.countKmersFromMarkerKmerIds(threadCount);
        
        // Retrieve peak and set thresholds.
        const uint64_t coveragePeak = assembler.assemblerInfo->kmerDistributionInfo.coveragePeak;
        const uint64_t minFreq = 4;
        const uint64_t maxFreq = 5 * coveragePeak;
        const uint64_t distinctKmerCount = assembler.kmerCounter->kmerIdFrequencies.size();

        cout << "Analyzing " << distinctKmerCount << " distinct minimizer k-mers." << endl;
        cout << "Filtering minimizers: Peak coverage is " << coveragePeak << "." << endl;
        cout << "Keeping k-mers with frequency [" << minFreq << ", " << maxFreq << "]." << endl;
             
        // Prune the existing markers in-place using the KmerCounter and markerKmerIds.
        assembler.applyKmerCountFilter(minFreq, maxFreq, threadCount);

        // Initialize KmerChecker for HttpServer diagnostics (optional).
        cout << "Initializing KmerChecker for diagnostics." << endl;
        assembler.createKmerChecker(assemblerOptions.kmersOptions, threadCount);
            
    } else {
        // Use the default k-mer based method.
        // Initialize the KmerChecker, which has the information needed
        // to decide if a k-mer is a marker.
        assembler.createKmerChecker(assemblerOptions.kmersOptions, threadCount);

        // Find the markers in the reads.
        assembler.findMarkers(threadCount);

        // Compute marker KmerIds (required for LowHash and alignment).
        // The SIMD path already creates these, but findMarkers does not.
        assembler.computeMarkerKmerIds(threadCount);
    }
    assembler.initiateSaveBinaryData(&Assembler::saveMarkers);



    // // Flag palindromic reads.
    // // These will be excluded from further processing.
    // if(!assemblerOptions.readsOptions.palindromicReads.skipFlagging) {
    //     assembler.palindromicMinAlignedMarkerCount = assemblerOptions.readsOptions.palindromicReads.minAlignedMarkerCount;
    //     assembler.palindromicMaxUncoveredBases = assemblerOptions.readsOptions.palindromicReads.maxUncoveredBases;
    //     assembler.flagPalindromicReads(
    //         assemblerOptions.readsOptions.palindromicReads.maxSkip,
    //         assemblerOptions.readsOptions.palindromicReads.maxDrift,
    //         assemblerOptions.readsOptions.palindromicReads.maxMarkerFrequency,
    //         assemblerOptions.readsOptions.palindromicReads.alignedFractionThreshold,
    //         assemblerOptions.readsOptions.palindromicReads.nearDiagonalFractionThreshold,
    //         assemblerOptions.readsOptions.palindromicReads.deltaThreshold,
    //         threadCount);
    // }


    // Compute maxChainLimit from coverage (Hifiasm parity: max(100, hom_cov * 5))
    const uint64_t coveragePeak = assembler.assemblerInfo->kmerDistributionInfo.coveragePeak;
    const uint64_t maxChainLimit = std::max(100UL, coveragePeak * 5);

    // Always build the inverted index for k-mer lookups (needed by both paths)
    assembler.buildInvertedIndex(threadCount);

    // Find and chain alignment candidates.
    if(!assemblerOptions.commandLineOnlyOptions.overlapsFromPafFile.empty()) {
        // PAF path: Import candidate pairs from PAF, then chain them using the inverted index.
        assembler.importAlignmentCandidatesFromPaf(assemblerOptions.commandLineOnlyOptions.overlapsFromPafFile);
        assembler.chainPafCandidates(
            assemblerOptions.overlapCandidatesOptions.driftRateTolerance,
            maxChainLimit,
            assemblerOptions.overlapCandidatesOptions,
            uint32_t(std::max(0, assemblerOptions.overlapCandidatesOptions.minChainMarkerCount)),
            threadCount
        );
    } else {
        // Inverted Index path: Discover candidate pairs via k-mer matches and chain them.
        assembler.chainAlignmentCandidates(
            assemblerOptions.overlapCandidatesOptions.driftRateTolerance,
            maxChainLimit,
            assemblerOptions.overlapCandidatesOptions,
            uint32_t(std::max(0, assemblerOptions.overlapCandidatesOptions.minChainMarkerCount)),
            threadCount
        );
    }

    // Compute alignments with variant evidence storage.
    assembler.computeAlignmentsWithEvidence(
        assemblerOptions.alignOptions,
        threadCount);

    // For http server and debugging/development purposes, generate an exhaustive table of candidates.
    // This can be done after alignment computation (it depends only on the candidate list).
    assembler.computeCandidateTable();

    // // Filter secondary/redundant alignments per read pair (Hifiasm Parity)
    // assembler.filterSecondaryAlignmentsPerReadPair(
    //     threadCount,
    //     assemblerOptions.readGraphOptions.filterSecondaryRequireNonRedundantOnBothReads);

    // // =========================================================================
    // // Overlap Filtering + Clean ReadGraph
    // // =========================================================================
    // // Build global mismatch-site clusters before EC parity paths.
    // // This guarantees global-site construction runs first for both parity modes.
    // cout << timestamp << "Precomputing global mismatch-site clusters before EC parity..." << endl;
    // const auto preEcGlobalHetClusters = assembler.clusterMismatchingPositionsIntoGlobalHetSites(
    //     assemblerOptions.alignOptions,
    //     threadCount,
    //     false,  // includeDeletedAlignments
    //     false   // readGraphOnly
    // );
    // cout << timestamp << "Precomputed global mismatch-site clusters: clusters="
    //      << preEcGlobalHetClusters.clusterRepresentatives.size()
    //      << " nodes=" << preEcGlobalHetClusters.nodes.size()
    //      << endl;

    // Default path: Hifiasm-style overlap filtering/parity (ha_ec + ha_ec_ff semantics).
    // Optional path: experimental global-site phasing/parity.
    const bool useGlobalSiteEcParity = (::getenv("DINARA_USE_GLOBAL_SITE_EC") != nullptr);
    if (useGlobalSiteEcParity) {
        cout << timestamp << "Using experimental global-site EC parity path." << endl;
        assembler.performGlobalSiteECParity(threadCount);
    } else {
        assembler.performHifiasmECParity(threadCount);
    }

    // assembler.performHifiasmECFinalFilteringParity(threadCount);
    // Clean overlap filtering (ma_hit_sub/cut/flt/contained + chimera detection) and read graph creation.
    // This uses conservative AND parity semantics (both reads must keep the overlap).
    assembler.createReadGraph6(threadCount);

    // Global mismatch-site diagnostics and export are expensive and intended for debugging.
    // Keep them off by default in production runs to preserve assembly throughput.
    const bool runGlobalHetDiagnostics = (::getenv("DINARA_ENABLE_GLOBAL_HET_DEBUG") != nullptr);
    if (runGlobalHetDiagnostics) {
    // Global mismatch sites + full per-allele member lists using only readGraph overlaps.
    // This is the fastest way to approximate "pileup across all reads" without a reference.
    {
        const ReadId focalReadId = ReadId(3);
        const Reads& reads = assembler.getReads();

        const auto clusters = assembler.clusterMismatchingPositionsIntoGlobalHetSitesReachableFromRead(
            focalReadId,
            assemblerOptions.alignOptions,
            threadCount,
            0,      // maxReadsToProcess
            0,      // maxAlignmentsToProcess
            false,  // includeDeletedAlignments
            true    // readGraphOnly
        );
        const uint32_t clusterSiteCount = uint32_t(
            clusters.clusterMemberOffsets.empty() ? 0 : (clusters.clusterMemberOffsets.size() - 1));

        static const char baseToAscii[] = {'A', 'C', 'G', 'T'};

        // Collect mismatch-defined sites that explicitly involve focalReadId (it has a mismatch at that site).
        struct FocalMismatchSite {
            uint32_t focalPos = 0;
            uint32_t siteId = 0;
        };
        vector<FocalMismatchSite> focalMismatchSites;
        focalMismatchSites.reserve(clusterSiteCount);
        for (size_t siteId = 0; siteId + 1 < clusters.clusterMemberOffsets.size(); siteId++) {
            const uint64_t begin = clusters.clusterMemberOffsets[siteId];
            const uint64_t end = clusters.clusterMemberOffsets[siteId + 1];
            uint32_t bestPos = std::numeric_limits<uint32_t>::max();
            for (uint64_t i = begin; i < end; i++) {
                const auto& node = clusters.nodes[clusters.clusterMembers[i]];
                if (node.first == focalReadId) {
                    bestPos = std::min(bestPos, node.second);
                }
            }
            if (bestPos != std::numeric_limits<uint32_t>::max()) {
                focalMismatchSites.push_back(FocalMismatchSite{bestPos, uint32_t(siteId)});
            }
        }
        sort(focalMismatchSites.begin(), focalMismatchSites.end(),
            [](const FocalMismatchSite& a, const FocalMismatchSite& b) {
                if (a.focalPos != b.focalPos) {
                    return a.focalPos < b.focalPos;
                }
                return a.siteId < b.siteId;
            });
        cout << timestamp << "Read" << focalReadId
             << " mismatch sites (readGraph clusters): " << focalMismatchSites.size() << endl;

        const auto propagationStart = std::chrono::steady_clock::now();
        cout << timestamp << "GlobalHetSite member propagation (readGraph): starting..." << endl;
        const auto members = assembler.computeGlobalHetSiteAlleleMembersUsingReadGraph(
            clusters,
            assemblerOptions.alignOptions,
            0,      // maxPendingTasks
            false,  // includeDeletedAlignments
            focalReadId
        );
        const auto propagationSeconds = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - propagationStart).count();

        cout << timestamp << "GlobalHetSite member propagation (readGraph): sites=" << clusterSiteCount
             << " propagatedAssignments=" << members.propagatedAssignments
             << " mappingHoles=" << members.mappingHoles
             << " mappingConflicts=" << members.mappingConflicts
             << " elapsedSec=" << propagationSeconds
             << endl;

        // Compute a consistent strand assignment for reads reachable from the focal read in the read graph.
        // This lets us export member positions (and alleles) in a single, focal-oriented coordinate frame.
        uint64_t strandConflicts = 0;
        vector<int8_t> strandByRead = assembler.computeReadGraphStrandsFromSeed(
            focalReadId,
            strandConflicts,
            false // includeDeletedAlignments
        );
        {
            uint64_t assigned = 0;
            for (const int8_t v : strandByRead) {
                if (v != -1) {
                    assigned++;
                }
            }
            cout << timestamp << "ReadGraph strand assignment: assigned=" << assigned
                 << " conflicts=" << strandConflicts << endl;
        }

        static const uint8_t complementBase[4] = {3, 2, 1, 0};
        const auto orientedMembers = assembler.orientGlobalHetSiteAlleleMembers(members, strandByRead);

        // Spot-check: verify that a few sites involving the focal read have self-consistent positions
        // under a readGraph-only multi-source traversal seeded from the mismatch members.
        {
            const size_t checkCount = std::min<size_t>(3, focalMismatchSites.size());
            for (size_t i = 0; i < checkCount; i++) {
                const uint32_t siteId = focalMismatchSites[i].siteId;
                const auto stats = assembler.debugVerifyGlobalHetSitePositionsUsingReadGraph(
                    clusters,
                    members,
                    siteId,
                    assemblerOptions.alignOptions,
                    20000,   // maxNodesToVisit
                    200000,  // maxAlignmentsToScan
                    false    // includeDeletedAlignments
                );
                cout << timestamp << "GlobalHetSite verify: siteId=" << siteId
                     << " expected=" << stats.expectedMembers
                     << " reached=" << stats.reachedMembers
                     << " checked=" << stats.checkedMappings
                     << " mismatched=" << stats.mismatchedPositions
                     << " holes=" << stats.mappingHoles
                     << " fails=" << stats.mappingFailures
                     << " hitNodeLimit=" << stats.hitNodeLimit
                     << " hitAlignmentLimit=" << stats.hitAlignmentLimit
                     << endl;
            }
        }

        const uint32_t siteCount = uint32_t(members.offsets.size());

        // Precompute mismatch-member counts and per-allele mismatch counts once per site.
        vector<array<uint32_t, 4> > mismatchCountsForward(siteCount, array<uint32_t, 4>{0, 0, 0, 0});
        vector<array<uint32_t, 4> > mismatchCountsOriented(siteCount, array<uint32_t, 4>{0, 0, 0, 0});
        vector<uint64_t> mismatchMembersBySite(siteCount, 0);
        for (uint32_t siteId = 0; siteId < siteCount && (siteId + 1) < clusters.clusterMemberOffsets.size(); siteId++) {
            const uint64_t begin = clusters.clusterMemberOffsets[siteId];
            const uint64_t end = clusters.clusterMemberOffsets[siteId + 1];
            mismatchMembersBySite[siteId] = end - begin;
            for (uint64_t j = begin; j < end; j++) {
                const auto& node = clusters.nodes[clusters.clusterMembers[j]];
                uint8_t b = reads.getOrientedReadBase(OrientedReadId(node.first, 0), node.second).value;
                if (b >= 4) {
                    continue;
                }
                mismatchCountsForward[siteId][b]++;
                const int8_t s = (uint64_t(node.first) < strandByRead.size()) ? strandByRead[uint64_t(node.first)] : int8_t(-1);
                if (s == 1) {
                    b = complementBase[b];
                }
                mismatchCountsOriented[siteId][b]++;
            }
        }

        // Precompute oriented support counts and total members per site once.
        // These are reused in filtering, printing, and export.
        vector<array<uint64_t, 4> > orientedSiteCounts(siteCount, array<uint64_t, 4>{0, 0, 0, 0});
        vector<uint64_t> orientedSiteMembers(siteCount, 0);
        const uint32_t orientedCount = uint32_t(orientedMembers.offsets.size());
        for (uint32_t siteId = 0; siteId < siteCount && siteId < orientedCount; siteId++) {
            const auto& off = orientedMembers.offsets[siteId];
            for (int allele = 0; allele < 4; allele++) {
                orientedSiteCounts[siteId][allele] = off[allele + 1] - off[allele];
                orientedSiteMembers[siteId] += orientedSiteCounts[siteId][allele];
            }
        }

        // Keep only robust multiallelic sites: at least 2 alleles with support >= 3.
        static constexpr uint32_t minAlleleSupportForExport = 3;
        static constexpr uint32_t minAlleleCountForExport = 2;
        vector<uint8_t> sitePassesMultiallelic(siteCount, 0);
        for (uint32_t siteId = 0; siteId < siteCount; siteId++) {
            uint32_t supportedAlleles = 0;
            for (int allele = 0; allele < 4; allele++) {
                if (orientedSiteCounts[siteId][allele] >= minAlleleSupportForExport) {
                    supportedAlleles++;
                }
            }
            sitePassesMultiallelic[siteId] = uint8_t(supportedAlleles >= minAlleleCountForExport);
        }

        const auto readIndex = assembler.buildFilteredGlobalHetSiteReadIndex(
            members,
            minAlleleSupportForExport,
            minAlleleCountForExport
        );
        const uint32_t invalidPos = std::numeric_limits<uint32_t>::max();
        const vector<Assembler::GlobalHetSiteReadIndex::ReadSite> emptyFocalReadSites;
        const auto& focalReadSites =
            (uint64_t(focalReadId) < readIndex.sitesByRead.size()) ?
            readIndex.sitesByRead[uint64_t(focalReadId)] :
            emptyFocalReadSites;
        vector<uint32_t> focalReadPosBySite(siteCount, invalidPos);
        vector<char> focalReadAlleleBySite(siteCount, '?');
        for (const auto& s : focalReadSites) {
            if (s.siteId < siteCount) {
                focalReadPosBySite[s.siteId] = s.readPosition;
                focalReadAlleleBySite[s.siteId] = baseToAscii[s.allele];
            }
        }

        vector<Assembler::GlobalHetSiteReadIndex::ReadSite> filteredFocalReadSites;
        filteredFocalReadSites.reserve(focalReadSites.size());
        for (const auto& s : focalReadSites) {
            if (s.siteId < sitePassesMultiallelic.size() && sitePassesMultiallelic[s.siteId]) {
                filteredFocalReadSites.push_back(s);
            }
        }

        cout << timestamp << "Read" << focalReadId
             << " projected global het sites after multiallelic filter: "
             << filteredFocalReadSites.size() << " / " << focalReadSites.size()
             << " (need >= " << minAlleleCountForExport << " alleles with support >= "
             << minAlleleSupportForExport << ")" << endl;

        // Print 10 sites involving focalReadId.
        const size_t toPrint = std::min<size_t>(10, filteredFocalReadSites.size());
        for (size_t i = 0; i < toPrint; i++) {
            const uint32_t siteId = filteredFocalReadSites[i].siteId;
            const uint32_t focalPos = filteredFocalReadSites[i].readPosition;
            const uint64_t mismatchMembers =
                (siteId < mismatchMembersBySite.size()) ? mismatchMembersBySite[siteId] : 0;
            const auto mismatchCounts =
                (siteId < mismatchCountsForward.size()) ?
                mismatchCountsForward[siteId] :
                std::array<uint32_t, 4>{0, 0, 0, 0};
            const auto siteCounts =
                (siteId < orientedSiteCounts.size()) ?
                orientedSiteCounts[siteId] :
                std::array<uint64_t, 4>{0, 0, 0, 0};
            const uint64_t siteMembers = (siteId < orientedSiteMembers.size()) ? orientedSiteMembers[siteId] : 0;

            cout << timestamp
                 << "GlobalHetSite[" << i << "]"
                 << " read" << focalReadId << "Pos=" << focalPos
                 << " read" << focalReadId << "Allele=" << baseToAscii[filteredFocalReadSites[i].allele]
                 << " mismatchMembers=" << mismatchMembers
                 << " mismatchCounts(A,C,G,T)=(" << mismatchCounts[0] << "," << mismatchCounts[1] << "," << mismatchCounts[2] << "," << mismatchCounts[3] << ")"
                 << " siteMembers=" << siteMembers
                 << " siteCounts(A,C,G,T)=(" << siteCounts[0] << "," << siteCounts[1] << "," << siteCounts[2] << "," << siteCounts[3] << ")"
                 << " members={";

            // Show up to 8 members per allele.
            for (int allele = 0; allele < 4; allele++) {
                const uint64_t b0 = orientedMembers.offsets[siteId][allele];
                const uint64_t b1 = orientedMembers.offsets[siteId][allele + 1];
                const uint64_t show = std::min<uint64_t>(b1 - b0, 8);
                if (show == 0) {
                    continue;
                }
                cout << baseToAscii[allele] << ":{";
                for (uint64_t k = 0; k < show; k++) {
                    const auto& m = orientedMembers.members[b0 + k];
                    cout << m.orientedReadId.getReadId()
                         << (m.orientedReadId.getStrand() == 1 ? "rc" : "fw")
                         << "-" << m.position;
                    if (k + 1 < show) {
                        cout << ",";
                    }
                }
                if ((b1 - b0) > show) {
                    cout << ",...";
                }
                cout << "}";
            }
            cout << "}" << endl;
        }

        // Export all SNP sites involving read 0 (summary + full per-allele member list).
        {
            const string summaryFileName = "Read" + to_string(focalReadId) + "GlobalHetSitesSummary.tsv";
            const string membersFileName = "Read" + to_string(focalReadId) + "GlobalHetSitesMembers.tsv";
            std::ofstream summary(summaryFileName);
            std::ofstream membersOut(membersFileName);
            if (!summary || !membersOut) {
                cout << timestamp << "Failed to open export files for read " << focalReadId << "." << endl;
            } else {
                summary << "siteId\treadPos\tmismatchMembers\tmismatchA\tmismatchC\tmismatchG\tmismatchT"
                        << "\tsiteMembers\tsiteA\tsiteC\tsiteG\tsiteT\treadAllele\n";
                membersOut << "siteId\treadPos0\treadAllele\treadId\treadStrand\tposition0\tpositionForward0\treadLength\n";

                // Prefer the mismatch-defined sites for "SNP sites of read0".
                // If there are none, fall back to the propagated membership list.
                const bool useMismatchSites = !focalMismatchSites.empty();
                const size_t exportCount = useMismatchSites ? focalMismatchSites.size() : filteredFocalReadSites.size();
                vector<uint32_t> readLengths(reads.readCount(), 0);
                for (uint64_t iRead = 0; iRead < reads.readCount(); iRead++) {
                    const ReadId rid = ReadId(iRead);
                    readLengths[iRead] = uint32_t(reads.getRead(rid).baseCount);
                }
                size_t exportedCount = 0;
                size_t filteredOutCount = 0;
                for (size_t idx = 0; idx < exportCount; idx++) {
                    const uint32_t siteId = useMismatchSites ? focalMismatchSites[idx].siteId : filteredFocalReadSites[idx].siteId;
                    if (siteId >= readIndex.sitePassesFilter.size() || readIndex.sitePassesFilter[siteId] == 0) {
                        filteredOutCount++;
                        continue;
                    }
                    if (siteId >= sitePassesMultiallelic.size() || sitePassesMultiallelic[siteId] == 0) {
                        filteredOutCount++;
                        continue;
                    }
                    if (siteId >= focalReadPosBySite.size() || focalReadPosBySite[siteId] == invalidPos) {
                        // Keep per-read-consistency filtering strict for DP-ready exports.
                        filteredOutCount++;
                        continue;
                    }
                    const uint32_t readPos = focalReadPosBySite[siteId];
                    const char readAllele = focalReadAlleleBySite[siteId];

                    const uint64_t mismatchMembers =
                        (siteId < mismatchMembersBySite.size()) ? mismatchMembersBySite[siteId] : 0;
                    const auto mismatchCounts =
                        (siteId < mismatchCountsOriented.size()) ?
                        mismatchCountsOriented[siteId] :
                        std::array<uint32_t, 4>{0, 0, 0, 0};
                    const auto siteCounts =
                        (siteId < orientedSiteCounts.size()) ?
                        orientedSiteCounts[siteId] :
                        std::array<uint64_t, 4>{0, 0, 0, 0};

                    // Export members using the focal-oriented coordinate frame (strandByRead),
                    // plus the original forward coordinates for debugging.
                    const uint32_t readPos0 = readPos;
                    // Export members in oriented coordinates (position0/1) consistent with readStrand,
                    // plus forward positions for debugging.
                    for (int allele = 0; allele < 4; allele++) {
                        const uint64_t b0 = orientedMembers.offsets[siteId][allele];
                        const uint64_t b1 = orientedMembers.offsets[siteId][allele + 1];
                        for (uint64_t k = b0; k < b1; k++) {
                            const auto& om = orientedMembers.members[k];
                            const ReadId rid = om.orientedReadId.getReadId();
                            const Strand strand = om.orientedReadId.getStrand();
                            const uint32_t posOriented0 = om.position;
                            const uint32_t len = (uint64_t(rid) < readLengths.size()) ? readLengths[uint64_t(rid)] : 0;
                            if (len == 0 || posOriented0 >= len) {
                                continue;
                            }
                            const uint32_t posFwd0 = (strand == 1) ? ((len - 1U) - posOriented0) : posOriented0;

                            // Allele char is already in oriented frame by construction (bucketed by allele).
                            const char alleleChar = baseToAscii[allele];

                            membersOut << siteId << "\t" << readPos0
                                       << "\t" << alleleChar
                                       << "\t" << rid
                                       << "\t" << int(strand)
                                       << "\t" << posOriented0
                                       << "\t" << posFwd0
                                       << "\t" << len
                                       << "\n";
                        }
                    }

                    const uint64_t siteMembers = (siteId < orientedSiteMembers.size()) ? orientedSiteMembers[siteId] : 0;
                    summary << siteId << "\t" << readPos
                            << "\t" << mismatchMembers
                            << "\t" << mismatchCounts[0] << "\t" << mismatchCounts[1] << "\t" << mismatchCounts[2] << "\t" << mismatchCounts[3]
                            << "\t" << siteMembers
                            << "\t" << siteCounts[0] << "\t" << siteCounts[1] << "\t" << siteCounts[2] << "\t" << siteCounts[3]
                            << "\t" << readAllele
                            << "\n";
                    exportedCount++;
                }

                cout << timestamp << "Wrote read" << focalReadId << " global het sites to " << summaryFileName
                     << " and " << membersFileName
                     << " (sites=" << exportedCount
                     << ", filteredOut=" << filteredOutCount
                     << ", criteria: >= " << minAlleleCountForExport
                     << " alleles with support >= " << minAlleleSupportForExport
                     << ")." << endl;
            }
        }
    }
    } else {
        cout << timestamp << "Skipping global-het diagnostics/export. "
             << "Set DINARA_ENABLE_GLOBAL_HET_DEBUG=1 to enable." << endl;
    }

    // return;


    // vector<uint32_t> ids;
    //
    // // After performHifiasmECParity(...) (it sets DeleteReasonPhase + informative counts/scores).
    // assembler.getAllCisAlignmentIdsSortedByInformativeSites(ids);
    //
    // // Print the 10 first sorted alignmentId and the informative sites they share
    // const size_t n = std::min<size_t>(10, ids.size());
    // for(size_t i = 0; i < n; ++i) {
    //     const uint32_t alignmentId = ids[i];
    //     const auto& ad = assembler.alignmentData[alignmentId];
    //
    //     // NOTE: we do NOT currently compute the exact "shared" informative-site intersection.
    //     // We have per-side counts and an overlap score:
    //     //   ad.informativeHetSiteCount0, ad.informativeHetSiteCount1
    //     //   ad.informativeHetSiteScore = max(count0,count1)
    //     const uint32_t sharedLowerBound = std::min(ad.informativeHetSiteCount0, ad.informativeHetSiteCount1);
    //
    //     cout << "rank=" << i
    //          << " alignmentId=" << alignmentId
    //          << " reads=(" << ad.readIds[0] << "," << ad.readIds[1] << ")"
    //          << " informative0=" << ad.informativeHetSiteCount0
    //          << " informative1=" << ad.informativeHetSiteCount1
    //          << " score=" << ad.informativeHetSiteScore
    //          << " sharedLB=" << sharedLowerBound
    //          << "\n";
    // }
    //
    //
    //
    //
    // // If you also want to exclude anything with other delete reasons:
    // // assembler.getAllCisAlignmentIdsSortedByInformativeSites(ids, /*keptByBothSidesOnly=*/true);
    //
    // // `ids` is now: CIS in both views (no DeleteReasonPhase on either side),
    // // sorted by `alignmentData[id].informativeHetSiteScore` descending.
    //
    // // assembler.createReadGraphFromEcParityCisOverlaps(threadCount, /*rebuildDirectedReadGraph*/ false);
    // assembler.createReadGraphFromEcParityCisOverlapsCoveringInformativeSites(threadCount,  false);

    // Snapshot the broad keep-set used for marker-graph collapse.
    std::vector<bool> keepForMarkerGraph(assembler.alignmentData.size(), false);
    for (uint64_t i = 0; i < keepForMarkerGraph.size(); ++i) {
        keepForMarkerGraph[i] = (assembler.alignmentData[i].info.isInReadGraph != 0);
    }

    // Mode 3 assembly requires reads in raw representation (not RLE).
    DINARA_ASSERT(assemblerOptions.readsOptions.representation == 0);

    // The marker length must be even.
    DINARA_ASSERT((assembler.assemblerInfo->k %2) == 0);




    // To create a complete marker graph, generate all vertices
    // regardless of coverage, and allow duplicate markers on vertices.
    assembler.createMarkerGraphVertices(
        2,                                              // minVertexCoverage
        std::numeric_limits<uint64_t>::max(),           // maxVertexCoverage
        0,                                              // minVertexCoveragePerStrand
        false,                                           // allowDuplicateMarkers
        std::numeric_limits<double>::signaling_NaN(),   // For peak finder, unused because minVertexCoverage is not 0.
        invalid<uint64_t>,                              // For peak finder, unused because minVertexCoverage is not 0.
        threadCount);

    // Filter marker graph vertices whose marker k-mers are short-period repeats (including homopolymers).
    // This reduces unreliable anchors and artifacts in repetitive regions.
    assembler.filterMarkerGraphVerticesByRepeatKmers(threadCount);

    // Find the reverse complement of each marker graph vertex.
    // We need the reverse complement vertices to be populated for Mode 3 anchor generation.
    assembler.findMarkerGraphReverseComplementVertices(threadCount);

    // // Clean up of duplicate markers, if requested and necessary.
    // if(assemblerOptions.markerGraphOptions.allowDuplicateMarkers and
    //     assemblerOptions.markerGraphOptions.cleanupDuplicateMarkers) {
    //     assembler.cleanupDuplicateMarkers(
    //         threadCount,
    //         assembler.getMarkerGraphMinCoverageUsed(),    // Stored by createMarkerGraphVertices.
    //         assemblerOptions.markerGraphOptions.minCoveragePerStrand,
    //         assemblerOptions.markerGraphOptions.duplicateMarkersPattern1Threshold,
    //         true, true);
    //     }

    // Create edges of the marker graph.
    assembler.createMarkerGraphEdges(threadCount);
    assembler.findMarkerGraphReverseComplementEdges(threadCount);

    if(assemblerOptions.markerGraphOptions.writeVertexCoverageHistogram) {
        cout << timestamp << "Writing marker graph vertex coverage histogram to " <<
            assemblerOptions.markerGraphOptions.vertexCoverageHistogramFileName << "." << endl;
        assembler.markerGraph.writeVertexCoverageHistogram(
            assemblerOptions.markerGraphOptions.vertexCoverageHistogramFileName,
            assemblerOptions.markerGraphOptions.vertexCoverageHistogramCanonicalOnly);
    }

    // // Now that the marker graph is built from the broad read-graph overlap set,
    // // we can tighten the overlap set by pruning overlaps for contained reads.
    // // This changes the read graph (rebuilt below) but keeps the marker graph intact.
    // const uint64_t minOverlapLengthForContainment = 50;
    // const uint64_t maxHangForContainment = 1000;
    // const double maxHangRateForContainment = 0.8;
    // assembler.flagContainedReads(
    //     maxHangForContainment,
    //     maxHangRateForContainment,
    //     minOverlapLengthForContainment,
    //     threadCount);
    // {
    //     uint64_t containedFlagCount = 0;
    //     for (ReadId r = 0; r < assembler.getReads().readCount(); ++r) {
    //         if (assembler.getReads().getFlags(r).isContained) {
    //             ++containedFlagCount;
    //         }
    //     }
    //     cout << timestamp << "[DIAG] After flagContainedReads: containedReads=" << containedFlagCount << endl;
    // }



    // Declare anchors pointer here to avoid scope issues
    shared_ptr<mode3::Anchors> anchors;

    // Compute the coverage range for primary marker graph edges (anchors).
    // This is done BEFORE marker graph vertex creation so filtering happens at source.
    uint64_t minPrimaryCoverage = assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage;
    uint64_t maxPrimaryCoverage = assemblerOptions.assemblyOptions.mode3Options.maxAnchorCoverage;
    if((minPrimaryCoverage == 0) and (maxPrimaryCoverage == 0)) {
        tie(minPrimaryCoverage, maxPrimaryCoverage) = assembler.getPrimaryCoverageRange();
        cout << "Automatically determined: minAnchorCoverage = " << minPrimaryCoverage <<
            ", maxAnchorCoverage = " << maxPrimaryCoverage << endl;
        minPrimaryCoverage = uint64_t(std::round(
            double(minPrimaryCoverage) * assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverageMultiplier));
        maxPrimaryCoverage = uint64_t(std::round(
            double(maxPrimaryCoverage) * assemblerOptions.assemblyOptions.mode3Options.maxAnchorCoverageMultiplier));
        cout << "After applying specified multipliers: minAnchorCoverage = " << minPrimaryCoverage <<
            ", maxAnchorCoverage = " << maxPrimaryCoverage << endl;
    } else {
        cout << "Using minAnchorCoverage = " << minPrimaryCoverage <<
            ", maxAnchorCoverage = " << maxPrimaryCoverage << endl;
    }

    // // Robust vertex-based anchors: split marker graph vertices using surviving readGraph overlaps
    // // (bridge removal + quasi-clique peeling), then emit per-cluster anchors (+ RC anchors).
    // anchors = assembler.createAnchorsFromMarkerGraphVerticesSplitUsingReadGraph(
    //     minPrimaryCoverage, maxPrimaryCoverage, assemblerOptions.assemblyOptions.mode3Options, threadCount);


    anchors =
            make_shared<mode3::Anchors>(
                MappedMemoryOwner(assembler),
                assembler.getReads(),
                assembler.assemblerInfo->k,
                *assembler.markers,
                assembler.markerGraph,
                minPrimaryCoverage,
                maxPrimaryCoverage,
                threadCount,
                true); // createFromVertices


    // anchors = assembler.createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval(
    //     minPrimaryCoverage,
    //     maxPrimaryCoverage,
    //     threadCount,
    //     /*enableColinearityPeeling*/ false,
    //     /*minDominantFractionToPeel*/ 0.9);

    // anchors = assembler.createAnchorsFromMarkerGraphVerticesBestPerOverlapIntervalDecomposed(
    //     minPrimaryCoverage, maxPrimaryCoverage, threadCount);


















    // // Step 6c: For each contained read, keep only one best overlap (by dpScore) and prune all others.
    // // This is a diagnostic/experimental alternative to removing contained reads entirely.
    // assembler.pruneContainedReadsToOneBestOverlapByDpScore(threadCount);
    //
    //
    // std::vector<bool> keepForReadGraph(assembler.alignmentData.size(), false);
    // for (uint64_t i = 0; i < keepForReadGraph.size(); ++i) {
    //     if (!keepForMarkerGraph[i]) continue;
    //
    //     const auto& ad = assembler.alignmentData[i];
    //
    //     // Strict rule: only keep if both sides still keep it after your extra filtering.
    //     if (!ad.keptByBothSides()) continue;
    //
    //     keepForReadGraph[i] = true;
    // }
    //
    // // Rebuild the read graph using the tightened keep-set.
    // assembler.rebuildReadGraphUsingSelectedAlignments(
    //     std::move(keepForReadGraph),
    //     /*rebuildDirectedReadGraph*/false);














    // // Alternatives (disabled):
    // anchors = assembler.createAnchorsFromMarkerGraphVerticesBestPerOverlapIntervalDecomposed(
    //     minPrimaryCoverage, maxPrimaryCoverage, threadCount);
    // anchors = make_shared<mode3::Anchors>(
    //     MappedMemoryOwner(assembler),
    //     assembler.getReads(),
    //     assembler.assemblerInfo->k,
    //     *assembler.markers,
    //     assembler.markerGraph,
    //     minPrimaryCoverage,
    //     maxPrimaryCoverage,
    //     threadCount,
    //     true); // createFromVertices


    // // Construct the mode3::Anchors (for HTTP server visualization).
    // // This must be done BEFORE createShasta2Anchors.
    // if(assemblerOptions.assemblyOptions.mode3Options.anchorCreationMethod ==
    //     "FromMarkerGraphVerticesAtOverlapEvents") {
    //     anchors = assembler.createAnchorsFromMarkerGraphVerticesAtOverlapEvents(
    //         minPrimaryCoverage,
    //         maxPrimaryCoverage,
    //         threadCount);
    // } else if(assemblerOptions.assemblyOptions.mode3Options.anchorCreationMethod ==
    //     "FromMarkerGraphVerticesBestPerOverlapInterval") {
    //     anchors = assembler.createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval(
    //         minPrimaryCoverage,
    //         maxPrimaryCoverage,
    //         threadCount);
    // } else if(assemblerOptions.assemblyOptions.mode3Options.anchorCreationMethod ==
    //     "FromOverlapsBestPerOverlapInterval") {
    //     anchors = assembler.createAnchorsFromOverlapsBestPerOverlapInterval(
    //         minPrimaryCoverage,
    //         maxPrimaryCoverage,
    //         threadCount);
    // } else {
    //     anchors =
    //         make_shared<mode3::Anchors>(
    //             MappedMemoryOwner(assembler),
    //             assembler.getReads(),
    //             assembler.assemblerInfo->k,
    //             *assembler.markers,
    //             assembler.markerGraph,
    //             minPrimaryCoverage,
    //             maxPrimaryCoverage,
    //             threadCount,
    //             true); // createFromVertices
    // }
    

    // Compute oriented read journeys.
    anchors->computeJourneys(threadCount);

    // Run Mode 3 assembly (initializes mode3Assembler for HTTP server).
    assembler.mode3Assembly(threadCount, anchors, assemblerOptions.assemblyOptions.mode3Options, false);

    // Store elapsed time for assembly.
    const auto steadyClock1 = std::chrono::steady_clock::now();
    const auto userClock1 = boost::chrono::process_user_cpu_clock::now();
    const auto systemClock1 = boost::chrono::process_system_cpu_clock::now();
    const double elapsedTime = 1.e-9 * double((
        std::chrono::duration_cast<std::chrono::nanoseconds>(steadyClock1 - steadyClock0)).count());
    const double userTime = 1.e-9 * double((
        boost::chrono::duration_cast<boost::chrono::nanoseconds>(userClock1 - userClock0)).count());
    const double systemTime = 1.e-9 * double((
        boost::chrono::duration_cast<boost::chrono::nanoseconds>(systemClock1 - systemClock0)).count());
    const double averageCpuUtilization =
        (userTime + systemTime) / (double(std::thread::hardware_concurrency()) * elapsedTime);
    assembler.storeAssemblyTime(elapsedTime, averageCpuUtilization);

    // Store peak memory usage.
    uint64_t peakMemoryUsage = getPeakMemoryUsage();
    assembler.storePeakMemoryUsage(peakMemoryUsage);

    // Store other performance information.
    assembler.assemblerInfo->threadCount = threadCount;
    assembler.assemblerInfo->virtualCpuCount = std::thread::hardware_concurrency();
    assembler.assemblerInfo->totalAvailableMemory = getTotalPhysicalMemory();

    // Write a summary of read information.
    assembler.writeReadsSummary();

    // Write the assembly summary.
    ofstream html("AssemblySummary.html");
    assembler.writeAssemblySummary(html);
    ofstream json("AssemblySummary.json");
    assembler.writeAssemblySummaryJson(json);
    ofstream htmlIndex("index.html");
    assembler.writeAssemblyIndex(htmlIndex);

    // If --saveBinaryData was requested and Mode assembly is 3,
    // wait for save binary data threads to finish.
    if(not assembler.saveBinaryDataDirectory.empty()) {
        assembler.waitForSaveBinaryDataThreads();
    }

    performanceLog << timestamp << endl;
    performanceLog << "Assembly time statistics:\n"
        "    Elapsed seconds: " << elapsedTime << "\n"
        "    Elapsed minutes: " << elapsedTime/60. << "\n"
        "    Elapsed hours:   " << elapsedTime/3600. << "\n";
    performanceLog << "Average CPU utilization: " << averageCpuUtilization << endl;
    performanceLog << "Peak Memory usage: " << peakMemoryUsage << " bytes = " <<
        int(std::round(double(peakMemoryUsage) / (1024. * 1024. * 1024.)) ) << " GiB" << endl;


    // Create shasta2 anchors equivalent to the marker graph vertices.
    // This allows downstream processing using shasta2 tools.
    createShasta2Anchors(assembler, assemblerOptions, threadCount, anchors);

}



// This function sets nr_overcommit_hugepages for 2MB pages
// to a little below total memory.
// If the setting needs to be modified, it acquires
// root privilege via sudo. This may result in the
// user having to enter a password.
void dinara::main::setupHugePages()
{

    // Get the total memory size.
    const uint64_t totalMemoryBytes = sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);

    // Figure out how much memory we want to allow for 2MB pages.
    const uint64_t MB = 1024 * 1024;
    const uint64_t GB = MB * 1024;
    const uint64_t maximumHugePageMemoryBytes = totalMemoryBytes - 8 * GB;
    const uint64_t maximumHugePageMemoryHugePages = maximumHugePageMemoryBytes / (2 * MB);

    // Check what we have it set to.
    const string fileName = "/sys/kernel/mm/hugepages/hugepages-2048kB/nr_overcommit_hugepages";
    ifstream file(fileName);
    if(!file) {
        throw runtime_error("Error opening " + fileName + " for read.");
    }
    uint64_t currentValue = 0;
    file >> currentValue;
    file.close();

    // If it's set to at least what we want, don't do anything.
    // When this happens, root access is not required.
    if(currentValue >= maximumHugePageMemoryHugePages) {
        return;
    }

    // Use sudo to set.
    const string command =
        "sudo sh -c \"echo " +
        to_string(maximumHugePageMemoryHugePages) +
        " > " + fileName + "\"";
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " running command: " + command);
    }

}



// Implementation of --command saveBinaryData.
// This copies Data to DataOnDisk.
void dinara::main::saveBinaryData(
    const AssemblerOptions& assemblerOptions)
{
    DINARA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "saveBinaryData");

    // Locate the Data directory.
    const string dataDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/Data";
    if(!std::filesystem::exists(dataDirectory)) {
        throw runtime_error(dataDirectory + " does not exist, nothing done.");
    }

    // Check that the DataOnDisk directory does not exist.
    const string dataOnDiskDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/DataOnDisk";
    if(std::filesystem::exists(dataOnDiskDirectory)) {
        throw runtime_error(dataOnDiskDirectory + " already exists, nothing done.");
    }

    // Copy Data to DataOnDisk.
    const string command = "cp -rp " + dataDirectory + " " + dataOnDiskDirectory;
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " running command:\n" + command);
    }
    cout << "Binary data successfully saved." << endl;
}



// Implementation of --command cleanupBinaryData.
void dinara::main::cleanupBinaryData(
    const AssemblerOptions& assemblerOptions)
{
    DINARA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "cleanupBinaryData");

    // Locate the Data directory.
    const string dataDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/Data";
    if(!std::filesystem::exists(dataDirectory)) {
        cout << dataDirectory << " does not exist, nothing done." << endl;
        return;
    }

    // Unmount it and remove it.
    ::system(("sudo umount " + dataDirectory).c_str());
    const int errorCode = ::system(string("rm -rf " + dataDirectory).c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " removing " + dataDirectory);
    }
    cout << "Cleanup of " << dataDirectory << " successful." << endl;

    // If the DataOnDisk directory exists, create a symbolic link
    // Data->DataOnDisk.
    const string dataOnDiskDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/DataOnDisk";
    if(std::filesystem::exists(dataOnDiskDirectory)) {
        std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
        const string command = "ln -s DataOnDisk Data";
        ::system(command.c_str());
    }

}

// Implementation of --command explore.
void dinara::main::explore(
    const AssemblerOptions& assemblerOptions)
{
    // If a paf file was specified, find its absolute path
    // before we switch to the assembly directory.
    string alignmentsPafFileAbsolutePath;
    if(not assemblerOptions.commandLineOnlyOptions.alignmentsPafFile.empty()) {
        if(!std::filesystem::exists(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile)) {
            throw runtime_error(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile + " not found.");
        }
        if(!std::filesystem::is_regular_file(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile)) {
            throw runtime_error(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile + " is not a regular file.");
        }
        alignmentsPafFileAbsolutePath = filesystem::getAbsolutePath(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile);
    }

    // Go to the assembly directory.
    std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);

    // Check that we have the binary data.
    if(!std::filesystem::exists("Data")) {
        throw runtime_error("Binary directory \"Data\" not available "
        " in assembly directory " +
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory +
        ". Use \"--memoryMode filesystem\", possibly followed by "
        "\"--command saveBinaryData\" and \"--command cleanupBinaryData\" "
        "if you want to make sure the binary data are persistently available on disk. "
        "See the documentations are some of these options require root access."
        );
        return;
    }

    // Create the Assembler.
    Assembler assembler("Data/", false, 1, 0);

    // Set up the consensus caller.
    if(assembler.getReads().representation == 1) {
        cout << "Setting up consensus caller " <<
            assemblerOptions.assemblyOptions.consensusCaller << endl;
    }
    assembler.setupConsensusCaller(assemblerOptions.assemblyOptions.consensusCaller);

    // Access all available binary data.
    assembler.httpServerData.assemblerOptions = &assemblerOptions;
    assembler.accessAllSoft();

    string executablePath = filesystem::executablePath();
    // On Linux it will be something like - `/path/to/install_root/bin/dinara`

    string executableBinPath = executablePath.substr(0, executablePath.find_last_of('/'));
    string installRootPath = executableBinPath.substr(0, executableBinPath.find_last_of('/'));
    string docsPath = installRootPath + "/docs";

    if (std::filesystem::is_directory(docsPath)) {
        assembler.httpServerData.docsDirectory = docsPath;
    } else {
        cout << "Documentation is not available." << endl;
        assembler.httpServerData.docsDirectory = "";
    }

    // Load the paf file, if one was specified.
    if(not alignmentsPafFileAbsolutePath.empty()) {
        assembler.loadAlignmentsPafFile(alignmentsPafFileAbsolutePath);
    }

    // Start the http server.
    bool localOnly;
    bool sameUserOnly;
    if(assemblerOptions.commandLineOnlyOptions.exploreAccess == "user") {
        localOnly = true;
        sameUserOnly = true;
    } else if(assemblerOptions.commandLineOnlyOptions.exploreAccess == "local") {
        localOnly = true;
        sameUserOnly = false;
    } else if (assemblerOptions.commandLineOnlyOptions.exploreAccess == "unrestricted"){
        localOnly = false;
        sameUserOnly = false;
    } else {
        throw runtime_error("Invalid value specified for --exploreAccess. "
            "Only use this option if you understand its security implications."
        );
    }
    assembler.explore(
        assemblerOptions.commandLineOnlyOptions.port,
        localOnly,
        sameUserOnly);
}




void dinara::main::listCommands()
{
    cout << "Valid commands are:" << endl;
    for(const string& command: commands) {
        cout << command << endl;
    }
}
