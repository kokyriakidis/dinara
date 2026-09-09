// The static executable provides
// basic functionality and reduced performance.
// For full functionality use the shared library built
// under directory src.

// Dinara.
#include "Assembler.hpp"
#include "HetAnchorK.hpp"
#include "HetSitePreparation.hpp"
#include "platformDependent.hpp"
#include "AssemblerOptions.hpp"
#include "filesystem.hpp"

#include "Shasta2Anchors.hpp"
#include "Shasta2AnchorPair.hpp"
#include "Shasta2Journeys.hpp"
#include "Shasta2AnchorGraph.hpp"

#include "performanceLog.hpp"
#include "Reads.hpp"
#include "Tee.hpp"
#include "timestamp.hpp"

// hifiasm candidate-overlap detector (submodule). C API used to generate read
// overlaps (with aligned intervals and CIGARs) directly in memory.
#include "hifiasm_overlaps.h"


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
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <fstream>
#include <map>
#include <set>
#include <functional>
#include <unordered_map>
#include <unordered_set>



// Shasta 2 Integration
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

// Minimizer filtering



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

#ifdef DINARA_LONG_MARKERS
    // With capacity-128 Kmers (256-bit KmerId), k can be up to 126.
    constexpr uint64_t maxK = 126;
#else
    constexpr uint64_t maxK = 62;
#endif
    if(assemblerOptions.kmersOptions.k < 6 or
       uint64_t(assemblerOptions.kmersOptions.k) > maxK) {
        throw runtime_error("Invalid value specified for --Kmers.k. Must be between 6 and " +
            to_string(maxK) + ".");
    }

    // k must be even. Note this rules out matching hifiasm's own default of 51,
    // and it is not free: anchor positions are marker MIDPOINTS (marker start +
    // k/2), and a midpoint mirrors exactly under reverse complement only when
    // 2*(k/2) == k-1, i.e. only for ODD k. Measured at k=50 on the GIAB fixture:
    // midpoint1 - ((len-1) - midpoint0) == +1 for 200454/200454 markers, with no
    // other value occurring. That offset is uniform across every marker, so
    // relative positions, offsets and orderings are unaffected -- but code that
    // compares a strand-0 position against a strand-1 position expecting an
    // exact mirror will be off by one.
    //
    // The original reason for requiring even k is not recorded here and was not
    // rediscovered; the only concrete dependency found in the tree is the
    // dormant AssemblerAbpoaMultiSegmentMSA.cpp, which extends its backbone by
    // k/2 per side to cover a full k-mer. Do not relax this guard without
    // establishing what else assumes it.
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
    assembler.assemblerInfo->assemblyMode = 3;
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



namespace {

// ===========================================================================
// DESIGN NOTE: markers AND overlaps both run no-HPC at the same marker k.
//
// dinara MARKERS are exact-match anchors in RAW space: two reads "share a
// marker" only when the raw k-mer is byte-identical, so markers are sketched
// no-HPC. Overlap detection is run in the SAME space (hifiOpt.no_hpc = 1,
// k = w = Kmers.k) so hifiasm's native chain anchors land 1:1 on dinara marker
// positions and no re-chaining is needed. The overlap-path minimizer filter is
// likewise built no-HPC at the marker k and reused by the overlapper, so the
// frequency counts correspond exactly to the minimizers being filtered.
//
// Selecting at the marker k (not a longer hifiasm 51-mer) makes each marker the
// FULL canonical minimizer, which is required for strand symmetry: RC(marker)
// is then exactly the same k-mer on the other strand, so reverse-strand
// overlaps share KmerIds.
//
// Only revisit this if dinara markers move into HPC space (a large change
// touching read representation, coordinates, and alignment).
// ===========================================================================

// Shared-read-store bridge to hifiasm.
//
// dinara already holds every read in memory (2-bit packed, in class Reads).
// hifiasm's marker filter build and overlap detection each read the input
// FASTA/FASTQ files themselves, so running both re-reads the same bytes from
// disk twice. The store bridge (hifiasm_reads_store_load + the *_from_store
// entry points) instead loads the reads into hifiasm's read store ONCE from
// memory, so both stages run with no file I/O.
//
// The store carries raw ASCII bases. That matches dinara's in-memory reads
// only when Reads.representation == 0 (raw); with RLE (representation == 1)
// the in-memory bases are run-length-encoded and would NOT reproduce the
// sequences hifiasm expects, so the caller must fall back to the file path.
//
// This RAII helper converts dinara's 2-bit reads to ASCII, loads them, and
// releases the store on scope exit. hifiasm copies the bases into its own
// store during load, so the temporary ASCII buffers are freed immediately
// after each read is inserted (only one read's worth is live at a time).
class HifiasmReadStore {
public:
    HifiasmReadStore(const Reads& reads) : loaded(false)
    {
        const ReadId n = reads.readCount();
        readViews.resize(n);
        // One contiguous ASCII buffer per read, kept alive until after the
        // load call returns (hifiasm_read_t stores pointers into them).
        asciiSeqs.resize(n);
        for(ReadId i = 0; i < n; i++) {
            const LongBaseSequenceView seq = reads.getRead(i);
            const uint64_t len = seq.baseCount;
            std::string& s = asciiSeqs[i];
            s.resize(len);
            for(uint64_t j = 0; j < len; j++) {
                s[j] = seq[j].character();  // 2-bit base -> 'A'/'C'/'G'/'T'
            }
            const span<const char> name = reads.getReadName(i);

            hifiasm_read_t& r = readViews[i];
            r.seq = s.data();
            r.seq_len = len;
            r.name = name.data();
            r.name_len = uint32_t(name.size());
        }
        const int rc = hifiasm_reads_store_load(
            readViews.data(), uint64_t(readViews.size()));
        if(rc != 0) {
            throw runtime_error(
                "Failed to load reads into hifiasm store (code " +
                to_string(rc) + ").");
        }
        loaded = true;
        // hifiasm has copied the bases into its own 2-bit store; the ASCII
        // staging buffers are no longer needed.
        asciiSeqs.clear();
        asciiSeqs.shrink_to_fit();
    }

    ~HifiasmReadStore()
    {
        if(loaded) {
            hifiasm_reads_store_release();
        }
    }

    HifiasmReadStore(const HifiasmReadStore&) = delete;
    HifiasmReadStore& operator=(const HifiasmReadStore&) = delete;

private:
    bool loaded;
    vector<hifiasm_read_t> readViews;
    vector<std::string> asciiSeqs;
};

} // anonymous namespace


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
    const uint64_t averageReadLength =
        assembler.getReads().getTotalBaseCount() / assembler.getReads().readCount();
    cout << "Average read length: " << averageReadLength << " bp." << endl;



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

    // Load the reads into hifiasm's in-memory read store ONCE, so both the
    // marker filter build (no-HPC) and overlap detection run against the same
    // store with no file re-reads. See HifiasmReadStore above.
    //
    // The store carries raw ASCII bases, so it requires Reads.representation == 0
    // (raw): with RLE the in-memory bases are run-length-encoded and would not
    // reproduce the sequences hifiasm expects. hifiasm overlap detection is the
    // only supported path, and it produces the native chains markers are built
    // from, so RLE is rejected up front rather than silently taking a different
    // route. The store is held for the whole marker+overlap section and released
    // when hifiasmStore goes out of scope at the end of assemble().
    if(assemblerOptions.readsOptions.representation != 0) {
        throw runtime_error(
            "Reads.representation must be 0 (raw). The hifiasm overlap path "
            "requires raw bases in the in-memory read store; RLE is not "
            "supported.");
    }
    performanceLog << timestamp
        << "Loading " << assembler.getReads().readCount()
        << " reads into hifiasm in-memory store (no file re-reads)." << endl;
    HifiasmReadStore hifiasmStore(assembler.getReads());
    performanceLog << timestamp << "hifiasm read store loaded." << endl;

    // Markers and overlaps share one no-HPC minimizer space at the marker k
    // (Kmers.k), so hifiasm's native chain anchors land 1:1 on dinara markers.
    const int markerK = int(assemblerOptions.kmersOptions.k);
    if(markerK <= 0) {
        throw runtime_error("Marker k (Kmers.k) must be positive.");
    }

    // Build hifiasm's overlap-path minimizer filter (no-HPC high-occurrence
    // k-mer filter) over the loaded store. Overlap detection REUSES this exact
    // filter (see hifiOpt.filter below) so it skips its own filter-build pass
    // (one fewer yak count pass) and seeds against the identical k-mer set. The
    // filter is keyed by (k, w, HPC), so it MUST be built at the same marker k /
    // no-HPC the overlapper sketches at. Ownership is held here until overlaps
    // are generated, then destroyed.
    hifiasm_filter_opt_t filterOpt = {};
    filterOpt.threads = int(threadCount);
    filterOpt.k_mer_length = markerK;
    filterOpt.mz_win = markerK;
    filterOpt.is_hpc = 0;        // dinara markers are no-HPC
    filterOpt.min_read_len = -1; // keep all reads (match dinara's set)

    performanceLog << timestamp
        << "Building hifiasm overlap-path minimizer filter (no-HPC, k=w="
        << markerK << ") from in-memory read store." << endl;
    hifiasm_filter_t* overlapReuseFilter =
        hifiasm_build_filter_from_store(&filterOpt);
    if(overlapReuseFilter == nullptr) {
        throw runtime_error(
            "Failed to build hifiasm minimizer filter for markers.");
    }
    performanceLog << timestamp << "hifiasm minimizer filter built." << endl;

    // Markers are derived from hifiasm's native overlap chains AFTER detection
    // (createMarkersFromNativeChain), so no marker sketch runs here and marker
    // creation is deferred. The chain anchors are the only marker positions that
    // can become anchors; sketching a full marker set would just add
    // marker-graph singletons. Chain-derived markers cluster in aligned
    // interiors and omit a read's extreme minimizers, so the span-coverage read
    // filter is not meaningful and is not applied on this path.
    cout << "Deferring marker creation to hifiasm native chain." << endl;

    // Coverage distribution (coverageLow/Het/Hom) feeds phasing and EC.
    // hifiasm already computed the marker k-mer count histogram while building
    // the filter (ha_ft_gen), so we take its authoritative peaks directly
    // instead of a second dinara-side histogram pass over every marker.
    //
    // NOTE: the actual chaining-frequency cutoff (hifiasm's low_occ/high_occ)
    // is computed internally by the hifiasm submodule from its own hom_cov,
    // not from dinara's OverlapCandidates.maxChainingFreq -- neither
    // hifiasm_filter_opt_t nor hifiasm_ovlp_opt_t exposes a frequency-cutoff
    // field, so that option cannot currently be threaded through to it.
    {
        const int homCov = hifiasm_filter_hom_cov(overlapReuseFilter);
        const int hetCov = hifiasm_filter_het_cov(overlapReuseFilter);
        const int lowCov = hifiasm_filter_low_cov(overlapReuseFilter);
        auto& dist = assembler.assemblerInfo->kmerDistributionInfo;
        dist.coverageHom = (homCov > 0) ? uint64_t(homCov) : 0;
        dist.coverageHet = (hetCov > 0) ? uint64_t(hetCov) : dist.coverageHom;
        dist.coverageLow = (lowCov > 0) ? uint64_t(lowCov) : 2;
        cout << "Marker k-mer coverage from hifiasm peaks:"
                " low "  << dist.coverageLow <<
                ", het " << dist.coverageHet <<
                ", hom " << dist.coverageHom << endl;
    }

    // Initialize KmerChecker for HttpServer diagnostics (optional).
    cout << "Initializing KmerChecker for diagnostics." << endl;
    assembler.createKmerChecker(assemblerOptions.kmersOptions, threadCount);

    // Generate read overlaps with the bundled hifiasm library from the loaded
    // in-memory read store (no file re-reads) and import them. hifiasm exports,
    // per overlap, the candidate pair, its aligned interval, AND its native
    // dense chain (its own colinear-DP anchors); dinara maps those anchors
    // straight to marker ordinals rather than re-chaining. Overlaps come back in
    // memory keyed by read name.
    //
    // Overlap detection runs in dinara's marker space: no-HPC at the marker k,
    // reusing the prebuilt filter. This unifies the overlapper with marker
    // selection so (a) the reused filter's k-mer counts correspond to the
    // minimizers being filtered and (b) hifiasm's native chain anchors land 1:1
    // on dinara marker positions.
    // Align.useHifiasmBaseAlignment (default true) selects hifiasm's ALIGNED
    // overlap path: candidate detection followed by gen_hc_r_alin_ea's
    // base-level alignment and filter. That path is the default because it is
    // the only source of hifiasm's per-overlap CIGAR, and that CIGAR is the one
    // thing dinara cannot reconstruct for itself: it spans the full overlap BOX
    // (x_pos_s..x_pos_e), which chaining extends outward along the diagonal to a
    // read end (Hash_Table.cpp), whereas dinara's own A*PA2 layer can only align
    // BETWEEN consecutive chain anchors and so never covers the flanks. Measured
    // on the GIAB fixture: box 11337 bases vs chain span 9472. With a CIGAR over
    // the whole box, any position in it can be mapped between the two reads
    // (OverlapCigarStore::queryToTarget), which chain-only alignment cannot do.
    //
    // The costs this buys, both measured on the same fixture: hifiasm's
    // alignment adds ~1.09 s to overlap detection, and its filter drops the
    // candidate set from 47146 to 29851 (overlaps whose windows fail to align).
    // The CIGARs are not separable from the filter -- both come from the aligned
    // path -- so taking one means taking the other.
    //
    // Setting this false restores the raw pre-alignment candidate set (only
    // h_ec_lchain runs, no CIGAR): faster and more permissive, and the mode
    // tests/hifiasm_overlap_parity.sh verifies against real hifiasm's own
    // --dbg-ovec output. The one thing the raw mode changes that dinara must
    // compensate for itself: hifiasm's aligned path forward-adjusts a
    // reverse-strand overlap's t_start/t_end to natural-forward coordinates as
    // part of alignment; the raw path never does that (there is no alignment
    // step to do it in), so importAlignmentCandidatesFromMemory is told via
    // rawCandidates below to undo the missing adjustment itself -- see its
    // Assembler.hpp comment.
    hifiasm_ovlp_opt_t hifiOpt = {};
    hifiOpt.threads        = int(threadCount);
    hifiOpt.no_hpc         = 1;
    hifiOpt.k_mer_length   = markerK;
    hifiOpt.mz_win         = markerK;
    hifiOpt.filter         = overlapReuseFilter;
    hifiOpt.is_ont         = assemblerOptions.alignOptions.hifiasmIsOnt ? 1 : 0;
    hifiOpt.raw_candidates =
        assemblerOptions.alignOptions.useHifiasmBaseAlignment ? 0 : 1;
    performanceLog << timestamp
        << "Overlap detection: no-HPC, k=w=" << markerK
        << ", reusing prebuilt marker filter, "
        << (hifiOpt.raw_candidates ?
            "raw candidates only (no base-level alignment)." :
            "WITH hifiasm base-level alignment.") << endl;

    hifiasm_overlap_t* ov = nullptr;
    uint64_t nOv = 0;
    char* names = nullptr;
    uint64_t* nameOff = nullptr;
    uint64_t nReads = 0;
    // hifiasm's own already-computed base-level CIGAR (query-forward /
    // target-alignment-orientation, same convention as the native chain).
    // ha_detect_candidates_from_store() always computes this internally per
    // overlap (concatenated across that overlap's aligned windows) regardless
    // of whether it is requested here -- cigarOut only controls whether the
    // already-computed tokens get copied out to the caller or discarded, so
    // requesting it costs essentially nothing extra (measured: ~7.0s either
    // way on the GIAB HG002 chr1:15-15.4Mb fixture, 989 reads).
    // On by default; set DINARA_STORE_HIFIASM_CIGAR=0 to fall back to
    // interval-only import (dinara re-derives the per-segment CIGAR with
    // A*PA2 in computeBaseAlignmentsAndStore, as before this was wired up).
    const bool storeHifiasmCigar = [](){
        const char* env = std::getenv("DINARA_STORE_HIFIASM_CIGAR");
        return env == nullptr || env[0] != '0';
    }();
    uint16_t* cigar = nullptr;
    uint64_t cigarLen = 0;
    // Native dense chain anchors: hifiasm's own colinear-DP chain per overlap.
    uint64_t* chain = nullptr;
    uint64_t chainLen = 0;
    performanceLog << timestamp
        << "Generating overlaps with hifiasm (in memory) from the "
        << "loaded read store." << endl;
    const int rc = hifiasm_detect_overlaps_from_store(
        &hifiOpt, &ov, &nOv, &names, &nameOff, &nReads,
        storeHifiasmCigar ? &cigar : nullptr,
        storeHifiasmCigar ? &cigarLen : nullptr,
        &chain, &chainLen);
    if(rc != 0) {
        hifiasm_overlaps_mem_free(ov, names, nameOff, cigar);
        throw runtime_error("hifiasm in-memory overlap detection failed (code "
            + to_string(rc) + ").");
    }
    performanceLog << timestamp << "hifiasm produced " << nOv
        << " overlaps over " << nReads << " reads (in memory)." << endl;
    cout << timestamp << "hifiasm CIGAR: " << (storeHifiasmCigar ?
        (to_string(cigarLen) + " tokens requested and stored") :
        string("not requested (interval-only import)")) << endl;

    // dinara uses hifiasm's native dense chain anchors (its colinear-DP seeds,
    // exported per overlap), mapped directly to marker ordinals in
    // computeBaseAlignmentsAndStore. Because overlap detection runs no-HPC at
    // the marker k, each anchor lands on a marker 1:1, so no re-chaining is
    // needed. When storeHifiasmCigar is on, hifiasm's own base-level CIGAR is
    // ALSO stored per pair (HifiasmImportedCigarStore), available for reuse
    // instead of recomputing an equivalent alignment with A*PA2 -- see
    // ProjectedAlignment::constructFromHifiasmCigar (not yet wired into
    // computeBaseAlignmentsAndStoreThreadFunction as of this change; storage
    // only for now).
    assembler.importAlignmentCandidatesFromMemory(
        ov, nOv, names, nameOff, nReads,
        cigar, cigarLen,
        /*chain*/ chain, /*chainLen*/ chainLen, threadCount,
        assemblerOptions.overlapCandidatesOptions.minOverlapLength,
        /*rawCandidates*/ bool(hifiOpt.raw_candidates));

    // Deferred marker creation: build the marker set from hifiasm's native chain
    // anchors (the shared k-mer positions), then persist. Must happen while
    // ov/names/nameOff/chain are still alive and before
    // computeBaseAlignmentsAndStore (which maps chain anchors to marker ordinals).
    // createMarkersFromNativeChain builds BOTH markers and markerKmerIds (strand
    // 0 and the RC-mirrored strand 1), so computeMarkerKmerIds is not needed.
    assembler.assemblerInfo->k = markerK;
    assembler.createMarkersFromNativeChain(
        ov, nOv, names, nameOff, nReads,
        chain, chainLen, threadCount);
    assembler.initiateSaveBinaryData(&Assembler::saveMarkers);


    hifiasm_overlaps_mem_free(ov, names, nameOff, cigar);
    free(chain);  // native chain arena (plain uint64_t array; owned by caller)
    chain = nullptr;

    // The overlap-path filter (built above and reused for overlap detection) is
    // no longer needed once overlaps are generated.
    hifiasm_filter_destroy(overlapReuseFilter);
    overlapReuseFilter = nullptr;

    // Compute base-level pairwise alignments for all overlaps and store
    // the resulting CIGARs. These are used downstream for CIGAR-based
    // SNP/indel detection in the phasing windows.
    assembler.computeBaseAlignmentsAndStore(
        assemblerOptions.alignOptions,
        threadCount);

    // Build a vector of ReadIds sorted by read length (longest first).
    const Reads& reads = assembler.getReads();
    const ReadId readCount = reads.readCount();
    vector<ReadId> readIdsSortedByLength(readCount);
    iota(readIdsSortedByLength.begin(), readIdsSortedByLength.end(), ReadId(0));
    sort(readIdsSortedByLength.begin(), readIdsSortedByLength.end(),
        [&reads](ReadId a, ReadId b) {
            return reads.getReadRawSequenceLength(a) >
                   reads.getReadRawSequenceLength(b);
        });

    // For http server and debugging/development purposes, generate an exhaustive table of candidates.
    // This can be done after alignment computation (it depends only on the candidate list).
    assembler.computeCandidateTable();


    // assembler.phaseOverlaps(threadCount);
    // assembler.phaseOverlapsKmeans(threadCount);

    // assembler.performHifiasmECParity(threadCount);

    // ---- Post-phasing overlap cleaning ----

    // Delete internal (non-dovetail) overlaps: matches that dangle off both
    // reads instead of reaching a read end -- repeat-induced or spurious
    // overlaps that hifiasm's ma_hit2arc classifies as MA_HT_INT / MA_HT_SHORT.
    //
    // This drops ONLY internal / short overlaps. Containments (MA_HT_QCONT /
    // MA_HT_TCONT) are intentionally KEPT: contained reads are only *flagged*
    // later (flagContainedReads), never removed here, because dropping their
    // overlaps would fragment the graph.
    //
    // Coordinates: this uses the TIGHT CIGAR span (ad.qs/qe/ts/te). Internal
    // matches are only detectable on tight coordinates -- extending them to the
    // read tips (extendOverlapToReadBoundaries) snaps the smaller overhang to 0
    // on each side, forcing ext5 = ext3 = 0, so ma_hit2arc could never return
    // MA_HT_INT. With dinara's minCoverage=0, hifiasm's ma_hit_flt likewise
    // classifies on the real overlap span against raw read lengths, matching the
    // tight span here.
    // assembler.deleteInternalOverlaps(/* maxHang */ 1000, /* maxHangRate */ 0.8, /* minOverlapLength */ 50, threadCount);

    // Build the read graph used for marker graph vertex construction.
    // Includes only alignments that:
    // - Are not deleted (no deleteReasons on either side — filters out
    //   chimeric reads, multi-chain pairs, and other earlier removals).
    // - Are not trans (state 2) or cisDifferentCopy (state 3) on either side.
    // Cis (state 1) and unlabeled (state 0) alignments are kept.
    // The read graph edges drive the disjoint set merges in createMarkerGraphVertices.
    assembler.createReadGraphFromPhasingCisOverlaps();

    // Set min and max marker graph vertex coverage thresholds.
    // const uint64_t minAnchorCoverage = std::max((uint64_t)3, (uint64_t)(0.15 * double(coverageHet) / 2));
    // const uint64_t maxAnchorCoverage = (uint64_t)(1.5 * double(coverageHet));
    
    const uint64_t minVertexCoverage = assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage;
    const uint64_t maxVertexCoverage = std::numeric_limits<uint64_t>::max();

    // Build marker graph vertices by transitive alignment collapse.
    // Each alignment merges its aligned marker pairs into a disjoint set.
    // Connected components become vertices. Vertices are then filtered by:
    // - Coverage: must have [minVertexCoverage, maxVertexCoverage] markers.
    // - No duplicate reads: a vertex cannot contain two markers from the
    //   same readId (either strand). This prevents collapsing both strands
    //   of a read into the same vertex.
    assembler.createMarkerGraphVertices(
        minVertexCoverage,
        maxVertexCoverage,
        0,                                              // minVertexCoveragePerStrand (disabled)
        false,                                          // allowDuplicateMarkers
        std::numeric_limits<double>::signaling_NaN(),   // unused (minVertexCoverage != 0)
        invalid<uint64_t>,                              // unused (minVertexCoverage != 0)
        threadCount);

    // Repeat-kmer and low-complexity filtering now happens at the minimizer
    // stage (applyKmerCountFilter with filterRepeatKmers / filterLowComplexity),
    // so these marker-graph vertex filters are redundant -- the offending k-mers
    // never seed a vertex. Left commented out; re-enable if the minimizer-stage
    // filters are ever turned off.
    //
    // Remove vertices whose k-mer is a short-period tandem repeat (period 1-5,
    // including homopolymers). Thresholds: {6, 4, 4, 4, 4}. Removes ~18.5%.
    // assembler.filterMarkerGraphVerticesByRepeatKmers(threadCount);
    //
    // Remove vertices whose k-mer has low sequence complexity (distinct
    // sub-k-mers of lengths 1, 2, 3). Thresholds: {4, 12, 24}. Removes ~4.5%.
    // assembler.filterMarkerGraphVerticesByDistinctSubkmerCount(threadCount);

    // Remove vertices where the transitive collapse grouped reads at k-mer
    // positions outside their direct chaining range. For each pair of reads
    // in a vertex, check that the vertex ordinal falls within the chain
    // start/end for both reads. Vertices failing this check were created by
    // indirect transitive paths (A→B→C) where A and C have no direct
    // alignment support at that position — typically false merges in repeats.
    assembler.filterMarkerGraphVerticesByChainConsistency(threadCount);

    // Pair each marker graph vertex with its reverse complement vertex.
    // Required before anchor generation, which needs RC-consistent vertices.
    assembler.findMarkerGraphReverseComplementVertices(threadCount);


    const uint64_t minAnchorCoverage = assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage;
    const uint64_t maxAnchorCoverage = std::numeric_limits<uint64_t>::max();

    // const uint64_t minPrimaryCoverage = assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage;;
    // const uint64_t maxPrimaryCoverage = assemblerOptions.assemblyOptions.mode3Options.maxAnchorCoverage;;
    cout << "Using: minAnchorCoverage = " << minAnchorCoverage <<
        ", maxAnchorCoverage = " << maxAnchorCoverage << endl;


    const MappedMemoryOwner shasta2Owner = assembler.shasta2MappedMemoryOwner();
    
    assembler.shasta2Anchors = make_shared<Shasta2Anchors>(
        shasta2Owner,
        assembler.getReads(),
        assembler.assemblerInfo->k,
        *assembler.markers,
        assembler.markerGraph,
        threadCount,
        minAnchorCoverage,
        maxAnchorCoverage);
        auto& shasta2Anchors = assembler.shasta2Anchors;

    // External-anchor export is deferred until after per-window MSA het-anchor
    // generation (testAbpoaMultiSegmentMSA below), so that any anchors created
    // from detected het sites are included in the exported set. Journeys and
    // windows are built from shasta2Anchors but do not depend on the export,
    // so moving the export down is dependency-safe.
    const string externalAnchorsName =
        std::filesystem::absolute("Shasta2ExternalAnchors").string();

    // Compute journeys.
    cout << timestamp << "Creating Shasta2Journeys..." << endl;
    assembler.shasta2Journeys = make_shared<Shasta2Journeys>(
        2 * assembler.getReads().readCount(),
        shasta2Anchors,
        threadCount,
        shasta2Owner);
    auto& shasta2Journeys = assembler.shasta2Journeys;

    // Filter each read's journey to its longest well-supported anchor chain
    // (every consecutive pair sharing >= minCommonForBackbone reads, bounded
    // look-back maxSkipForBackbone). Runs per read independently and rewrites
    // the stored journeys + positionInJourney before any windowing decision, so
    // every downstream stage sees the cleaned chains.
    cout << timestamp << "Filtering journeys by anchor chaining..." << endl;
    shasta2Journeys->filterByAnchorChaining(
        assemblerOptions.assemblyOptions.mode3Options.minCommonForBackbone,
        assemblerOptions.assemblyOptions.mode3Options.maxSkipForBackbone,
        threadCount);

    // MSA-based overlap phasing — disabled, replaced by CIGAR-based window pipeline.
    // assembler.phaseOverlapsMSA(threadCount);

    // Flag contained reads so they can be excluded from inter-window edge discovery.
    cout << timestamp << "Flagging contained reads..." << endl;
    assembler.flagContainedReads(1000, 0.8, 0, threadCount);



    // Build the DETECTION anchor graph from the current (pre-het) journeys,
    // run POA-based het-site detection on its edges to append new het
    // anchors, rebuild journeys so the new anchors take their place in every
    // affected read's journey, and rebuild the anchor graph fresh from the
    // rebuilt journeys. This graph -- not the detection one -- is what
    // everything downstream (export, GFA, assembly graph) uses.
    //
    // This must run BEFORE journey tie resolution and external-anchor export
    // below, so shasta2Anchors holds every het anchor this pass creates before
    // either one runs -- otherwise anchors created here would silently never
    // reach the shasta2 export. See Shasta2AnchorGraphHetOnGraph.cpp's file
    // header for why detection only appends anchors and never touches a graph
    // directly.
    // The per-edge coverage threshold defaults to 0 (see
    // Assembly.mode3.minJourneyEdgeCoverage). filterByAnchorChaining has
    // already run (above), so the filtered journeys are the source of truth:
    // every consecutive pair surviving in a filtered journey becomes an edge.
    const uint64_t minEdgeCoverage =
        assemblerOptions.assemblyOptions.mode3Options.minJourneyEdgeCoverage;
    // Snapshot the anchor count before the detector appends anything,
    // so the journeys rebuild below (Shasta2Journeys::rebuildAfterNewAnchors)
    // can fold in exactly the new anchors it creates without also pulling in
    // any other het/hom anchors some other pass might append for a purpose
    // that doesn't include internal graph participation.
    const Shasta2AnchorId newAnchorsBegin = shasta2Anchors->size();

    // CIGAR-driven het anchors. This runs HERE, not next to the CIGAR import,
    // because it needs shasta2Anchors to append to and shasta2Journeys to be
    // rebuilt afterwards -- the imported CIGAR store lives for the whole run,
    // so detection is free to happen late.
    //
    // Same contract as the abPOA detector: detect, append one anchor per
    // allele arm, and let the caller rebuild journeys and then the anchor graph
    // from scratch. No surgery on an existing graph -- once a new anchor is
    // just another entry in each read's journey, the (independently verified)
    // journey->graph builder produces the right topology on its own, and reads
    // in neither arm take the direct flank-to-flank edge for free.
    if(assemblerOptions.assemblyOptions.mode3Options.detectSnpSites) {
        vector<Assembler::CigarSnpSite> snpSites;
        assembler.detectCigarSnpSites(
            assemblerOptions.assemblyOptions.mode3Options.snpSiteMinDisagree,
            assemblerOptions.assemblyOptions.mode3Options.snpSiteMinFraction,
            assemblerOptions.assemblyOptions.mode3Options.snpSiteMaxFraction,
            assemblerOptions.assemblyOptions.mode3Options.hetErrorRate,
            assemblerOptions.assemblyOptions.mode3Options.snpSiteStrandBiasPValue,
            assemblerOptions.assemblyOptions.mode3Options.snpSiteMinPurity,
            assemblerOptions.assemblyOptions.mode3Options.snpSiteMinAltDominance,
            assemblerOptions.assemblyOptions.mode3Options.snpSiteMinAlleleFraction,
            assemblerOptions.assemblyOptions.mode3Options.snpSiteFilterHomopolymer,
            assemblerOptions.assemblyOptions.mode3Options.snpSiteFilterStr,
            assemblerOptions.assemblyOptions.mode3Options.snpSiteAlleleCoverageRate,
            assemblerOptions.assemblyOptions.mode3Options.snpSiteAlleleCoverageFloor,
            assemblerOptions.assemblyOptions.mode3Options.snpSitePloidy,
            assemblerOptions.assemblyOptions.mode3Options.createSnpSiteAnchors ?
                &snpSites : nullptr,
            threadCount);

        if(assemblerOptions.assemblyOptions.mode3Options.createSnpSiteAnchors) {
            // Collapse sites that are the SAME LOCUS detected twice, then
            // drop arm members an existing anchor already holds. Both are pure
            // functions of their inputs and live in HetSitePreparation.cpp so
            // they can be tested; see that header for why each is needed.
            {
                const uint64_t duplicateSites = collapseDuplicateLoci(snpSites);
                cout << timestamp << "CIGAR het sites: dropped "
                     << duplicateSites << " duplicate detection(s) of a locus "
                        "already claimed by a better-supported site, leaving "
                     << snpSites.size() << "." << endl;

                const auto occupied = buildOccupiedPositions(*shasta2Anchors);
                const uint64_t alreadyAnchored = dropAlreadyAnchoredArmMembers(
                    snpSites, occupied, hetAnchorKHalf());
                cout << timestamp << "  " << alreadyAnchored
                     << " arm member(s) already sit in an anchor that isolates "
                        "their allele; het anchors are built from the rest." << endl;
            }

            // DINARA_HET_ARM_DUMP: one line per arm, "siteIndex armIndex" then
            // the member read names. The arms ARE the haplotype partition this
            // path exists to produce, so this is what lets a ground-truth check
            // ask whether the two arms of a site separate reads by haplotype.
            // Measured on the 989-read fixture against an independent minimap2
            // alignment: 82.2% of sites perfectly haplotype-pure, 99.08% mean
            // purity, 0.77% of members in the opposite arm.
            std::unique_ptr<ofstream> armDump;
            if(const char* path = std::getenv("DINARA_HET_ARM_DUMP")) {
                armDump = std::make_unique<ofstream>(path);
            }

            uint64_t created = 0, skippedThin = 0;
            for(uint64_t siteIndex = 0; siteIndex < snpSites.size(); siteIndex++) {
                const Assembler::CigarSnpSite& site = snpSites[siteIndex];
                uint64_t armIndex = 0;
                for(const auto& members: site.alleles) {
                    const uint64_t thisArm = armIndex++;
                    if(armDump) {
                        (*armDump) << siteIndex << '\t' << thisArm;
                        for(const auto& m: members) {
                            const auto nm =
                                assembler.getReads().getReadName(m.first.getReadId());
                            (*armDump) << '\t';
                            armDump->write(&*nm.begin(), std::streamsize(nm.size()));
                        }
                        (*armDump) << '\n';
                    }
                    // appendHetAnchorPair's own floor: an arm needs at least
                    // two members to be an anchor at all.
                    if(members.size() < 2) { ++skippedThin; continue; }
                    shasta2Anchors->appendHetAnchorPair(members);
                    ++created;
                }
            }
            cout << timestamp << "CIGAR het anchors: created " << created
                 << " from " << snpSites.size() << " sites ("
                 << skippedThin << " arms too thin)." << endl;
        }
    }

    // Rebuild journeys if the CIGAR-driven detector appended anchors. This is
    // the step that makes a new
    // anchor real: until a read's journey contains it, it is just an isolated
    // vertex the graph builder never links, so the assembly graph comes out
    // unchanged and the anchors do nothing. (That is exactly what happened when
    // this rebuild was still nested inside a since-removed branch.)
    if(shasta2Anchors->size() > newAnchorsBegin) {
        cout << timestamp << "Rebuilding journeys to include "
             << (shasta2Anchors->size() - newAnchorsBegin)
             << " new het anchors..." << endl;
        shasta2Journeys->journeyTiePreferHet =
            assemblerOptions.assemblyOptions.mode3Options.journeyTiePreferHet;
        shasta2Journeys->minAnchorCoverage =
            assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage;
        shasta2Journeys->rebuildAfterNewAnchors(newAnchorsBegin, threadCount);
    }
    cout << timestamp << "Creating Shasta2AnchorGraph from journeys "
         << "(final pass, minEdgeCoverage=" << minEdgeCoverage << ")..." << endl;
    assembler.shasta2AnchorGraph = make_shared<Shasta2AnchorGraph>(
        *shasta2Anchors,
        *shasta2Journeys,
        minEdgeCoverage,
        threadCount);
    // The export's drop map is DERIVED from the journeys (below), not resolved
    // a second time here. An earlier version did resolve it independently --
    // grouping each read's anchor occurrences by exported position and picking a
    // keeper -- which meant two implementations of one policy. They disagreed:
    // 446 drops agreed, 11 were derived-only and 27 resolved-only, and for a
    // while they even used opposite tie-break preferences, so the exported
    // anchors contradicted the exported graph. The journeys already carry the
    // answer, so the export reads it rather than recomputing it.
    Shasta2Anchors::ExternalAnchorDropMap journeyTieDropMap;

    // The drop map the export applies should be a FUNCTION OF THE JOURNEYS, not
    // a second opinion about them. The journeys were already rebuilt with ties
    // resolved; re-deriving that decision here on the anchor objects means two
    // implementations of one policy, which is exactly how they came to disagree
    // before (the rebuild kept het, the export kept primary, and the exported
    // anchors contradicted the exported graph).
    //
    // So recompute it directly: a member (anchor, read) is exported iff that
    // anchor actually occurs in that read's journey. Anything the rebuild
    // dropped is dropped here by construction, and nothing else is. The
    // independently-computed map is kept alongside only to report whether the
    // two agree -- if they ever diverge, the derivation is authoritative and the
    // difference is the bug.
    {
        // Straight from the rebuild's own record of what it dropped.
        uint64_t derivedDrops = 0;
        for(const auto& [canonicalId, readId]: shasta2Journeys->journeyTieDrops) {
            auto& vec = journeyTieDropMap[canonicalId];
            if(std::find(vec.begin(), vec.end(), readId) == vec.end()) {
                vec.push_back(readId);
                ++derivedDrops;
            }
        }
        cout << timestamp << "Export drop map taken from the journey rebuild: "
             << derivedDrops << " member(s) removed when journeys were rebuilt."
             << endl;
    }

    // Write external anchors. Deferred to here (after MSA het-anchor
    // generation) so newly generated het anchors are part of the exported set.
    // The drop map removes the members that would otherwise collide in a
    // per-read journey (see the tie resolution above).
    cout << timestamp << "Writing Shasta2 external anchors to "
         << externalAnchorsName << "..." << endl;
    const uint64_t exportedExternalAnchorCount =
        shasta2Anchors->writeExternalAnchors(
            externalAnchorsName, true, &journeyTieDropMap);
    cout << timestamp << "Wrote " << exportedExternalAnchorCount
         << " external anchors for Shasta2. Use --external-anchors-name "
         << externalAnchorsName << endl;
    // The export subtracts hetAnchorKHalf() uniformly from every stored midpoint
    // (see writeExternalAnchors), so shasta2 must be loaded with the MATCHING
    // --k: 2 by default, 0 for the experimental DINARA_HET_K=0 path. A mismatch
    // shifts every anchor by one base. Report it so the caller passes the right
    // value to the downstream shasta2 invocation.
    cout << timestamp << "Shasta2 must load these external anchors with --k "
         << hetAnchorK()
         << (hetAnchorK() == 0 ? " (EXPERIMENTAL DINARA_HET_K=0)." : ".") << endl;

    // Verify and finalize the anchor graph (already built above for the
    // journey path), then export it.
    {
        // VERIFICATION: independently recompute the journey-consecutive
        // (anchorA, anchorB) tally directly from the rebuilt journeys
        // (bypassing findChildren/Shasta2AnchorPair entirely) and cross-check
        // against every edge in the FINAL anchor graph.
        if(const char* env = std::getenv("DINARA_VERIFY_ANCHOR_GRAPH")) {
            if(env[0] != '0') {
                cout << timestamp << "Verifying anchor graph construction "
                        "against an independent journey walk..." << endl;
                std::unordered_map<uint64_t, uint64_t> independentTally;
                const uint64_t orientedReadCount = shasta2Journeys->size();
                for(uint64_t v = 0; v < orientedReadCount; v++) {
                    const OrientedReadId orientedReadId = OrientedReadId::fromValue(v);
                    const auto journey = (*shasta2Journeys)[orientedReadId];
                    for(uint64_t i = 0; i + 1 < journey.size(); i++) {
                        const uint64_t a = uint64_t(journey[i]);
                        const uint64_t b = uint64_t(journey[i + 1]);
                        independentTally[(a << 32) | b]++;
                    }
                }

                uint64_t edgesChecked = 0, edgesMismatched = 0, edgesUnsupported = 0;
                auto edgeRange = boost::edges(*assembler.shasta2AnchorGraph);
                for(auto it = edgeRange.first; it != edgeRange.second; ++it) {
                    const auto& edge = (*assembler.shasta2AnchorGraph)[*it];
                    const uint64_t a = edge.anchorIdA;
                    const uint64_t b = edge.anchorIdB;
                    const uint64_t declaredCoverage = edge.coverage();
                    const auto tallyIt = independentTally.find((a << 32) | b);
                    ++edgesChecked;
                    if(tallyIt == independentTally.end()) {
                        ++edgesUnsupported;
                        cout << timestamp << "  MISMATCH: edge " << a << "->" << b
                             << " coverage=" << declaredCoverage
                             << " but independent walk found ZERO occurrences." << endl;
                    } else if(tallyIt->second != declaredCoverage) {
                        ++edgesMismatched;
                        cout << timestamp << "  MISMATCH: edge " << a << "->" << b
                             << " coverage=" << declaredCoverage
                             << " independent tally=" << tallyIt->second << endl;
                    }
                }

                uint64_t missingEdges = 0;
                for(const auto& [key, count] : independentTally) {
                    if(count < minEdgeCoverage) continue;
                    const Shasta2AnchorId a = key >> 32;
                    const Shasta2AnchorId b = key & 0xffffffffULL;
                    const auto found = boost::edge(a, b, *assembler.shasta2AnchorGraph);
                    if(!found.second) {
                        ++missingEdges;
                        if(missingEdges <= 20) {
                            cout << timestamp << "  MISSING EDGE: " << a << "->" << b
                                 << " independent tally=" << count
                                 << " (>= minEdgeCoverage=" << minEdgeCoverage
                                 << ") but no graph edge exists." << endl;
                        }
                    }
                }

                cout << timestamp << "Anchor graph verification: " << edgesChecked
                     << " graph edges checked (" << edgesMismatched << " coverage mismatches, "
                     << edgesUnsupported << " with zero independent support), "
                     << independentTally.size() << " distinct journey-consecutive pairs found, "
                     << missingEdges << " missing edges (tally >= minEdgeCoverage but no graph edge)."
                     << endl;
            }
        }

        // Tip cleanup: a het/hom anchor can still end up one-sided (e.g. if
        // every one of its member reads happens to have it as the first/last
        // entry of its filtered journey); iterate until stable.
        for(;;) {
            const uint64_t hetTips =
                assembler.shasta2AnchorGraph->removeHetArmTips(*shasta2Anchors);
            if(hetTips == 0) break;
        }

        assembler.shasta2AnchorGraph->writeGfa("Shasta2AnchorGraph.gfa");
        assembler.shasta2AnchorGraph->writeCsv("Shasta2AnchorGraph.csv");
        cout << timestamp << "Wrote Shasta2AnchorGraph.gfa / .csv" << endl;

        // Persist the anchor graph as binary, not just the GFA/CSV views.
        //
        // saveForShasta2 writes the shasta2-compatible MemoryMapped format
        // (a boost archive of shasta2::AnchorGraph). This is the ONLY file
        // shasta2 can load via --external-anchor-graph-name; dinara's own
        // save() format is a different, incompatible archive and will segfault
        // shasta2 if passed there. Write it into the run's Data/ directory under
        // the same "Data/Shasta2-Shasta2AnchorGraph" name the binary-data path
        // uses, so the standard --external-anchor-graph-name works unchanged.
        string externalAnchorGraphName =
            assembler.shasta2MappedMemoryOwner().largeDataName("Shasta2AnchorGraph");
        if(externalAnchorGraphName.empty()) {
            // No binary-data directory (e.g. memory-mode anonymous): fall back to
            // a plain file in the working directory.
            externalAnchorGraphName = "Shasta2ExternalAnchorGraph";
        }
        externalAnchorGraphName =
            std::filesystem::absolute(externalAnchorGraphName).string();
        assembler.shasta2AnchorGraph->saveForShasta2(
            externalAnchorGraphName, *shasta2Anchors, &journeyTieDropMap);
        cout << timestamp << "Wrote shasta2 anchor graph. Use "
             << "--external-anchor-graph-name " << externalAnchorGraphName << endl;
    } 

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
    const uint64_t peakMemoryUsage = getPeakMemoryUsage();
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

    return;
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
