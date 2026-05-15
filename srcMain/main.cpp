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
#include "mode3-BidirectedAnchor.hpp"
#include "mode3-DirectedAnchors.hpp"
#include "mode3-DirectedAnchorGraph.hpp"
#include "mode3-AnchorGraph.hpp"
#include "mode3-AnchorGraphSuperbubbles.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2AnchorsFromSplitVertices.hpp"
#include "Shasta2Journeys.hpp"
#include "Shasta2AnchorGraph.hpp"
#include "Shasta2AssemblyGraph.hpp"
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
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include "iostream.hpp"
#include <set>
#include <unordered_set>

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


namespace {

using BubbleEndpointPair = std::pair<std::string, std::string>;

BubbleEndpointPair canonicalBubblePair(std::string a, std::string b)
{
    if(a <= b) {
        return {std::move(a), std::move(b)};
    }
    return {std::move(b), std::move(a)};
}


std::string normalizeBubbleEndpoint(std::string s)
{
    if(!s.empty() && (s.back() == '+' || s.back() == '-')) {
        s.pop_back();
    }
    return s;
}


std::string shellQuote(const std::string& s)
{
    std::string quoted = "'";
    for(const char c: s) {
        if(c == '\'') {
            quoted += "'\"'\"'";
        } else {
            quoted += c;
        }
    }
    quoted += "'";
    return quoted;
}


bool envFlagIsEnabled(const char* variableName)
{
    const char* value = ::getenv(variableName);
    if(value == nullptr) {
        return false;
    }
    std::string s(value);
    std::transform(
        s.begin(),
        s.end(),
        s.begin(),
        [](unsigned char c) { return char(std::tolower(c)); });
    return
        s == "1" ||
        s == "true" ||
        s == "yes" ||
        s == "on";
}


bool envFlagIsDisabled(const char* variableName)
{
    const char* value = ::getenv(variableName);
    if(value == nullptr) {
        return false;
    }
    std::string s(value);
    std::transform(
        s.begin(),
        s.end(),
        s.begin(),
        [](unsigned char c) { return char(std::tolower(c)); });
    return
        s == "0" ||
        s == "false" ||
        s == "no" ||
        s == "off";
}


uint64_t envUintOrDefault(const char* variableName, uint64_t defaultValue)
{
    const char* value = ::getenv(variableName);
    if(value == nullptr || *value == '\0') {
        return defaultValue;
    }
    try {
        const uint64_t parsed = std::stoull(value);
        return std::max<uint64_t>(1, parsed);
    } catch(...) {
        return defaultValue;
    }
}


void maybeRunBubbleFinderDirectedSuperbubbleComparison(
    const mode3::AnchorGraph& anchorGraph,
    const std::vector<std::pair<
        mode3::AnchorGraph::vertex_descriptor,
        mode3::AnchorGraph::vertex_descriptor> >& onoderaSuperbubbles,
    uint64_t defaultThreadCount)
{
    if(envFlagIsDisabled("DINARA_BUBBLEFINDER_COMPARE")) {
        cout << "[SuperbubbleDetection] BubbleFinder comparison disabled by "
            "DINARA_BUBBLEFINDER_COMPARE." << endl;
        return;
    }

    std::string bubbleFinderBinary = "BubbleFinder";
    if(const char* value = ::getenv("DINARA_BUBBLEFINDER_BIN")) {
        if(*value != '\0') {
            bubbleFinderBinary = value;
        }
    }

    const uint64_t bubbleFinderThreadCount =
        envUintOrDefault("DINARA_BUBBLEFINDER_THREADS", defaultThreadCount);
    const uint64_t bubbleFinderSampleCount =
        envUintOrDefault("DINARA_BUBBLEFINDER_SAMPLE_COUNT", 10);

    const std::string gfaFileName = "AnchorGraph-BubbleFinder.gfa";
    const std::string outputFileName =
        "AnchorGraph-BubbleFinder.directed-superbubbles.txt";
    const std::string stdoutLogFileName = "AnchorGraph-BubbleFinder.stdout.log";
    const std::string stderrLogFileName = "AnchorGraph-BubbleFinder.stderr.log";

    std::unordered_set<std::string> nodeLabels;
    nodeLabels.reserve(num_vertices(anchorGraph));
    for(auto vp = boost::vertices(anchorGraph); vp.first != vp.second; ++vp.first) {
        const auto v = *vp.first;
        nodeLabels.insert(std::to_string(anchorGraph.getAnchorId(v)));
    }

    std::vector<std::pair<std::string, std::string> > directedEdges;
    directedEdges.reserve(num_edges(anchorGraph));
    for(auto ep = boost::edges(anchorGraph); ep.first != ep.second; ++ep.first) {
        const auto e = *ep.first;
        const auto v0 = source(e, anchorGraph);
        const auto v1 = target(e, anchorGraph);
        directedEdges.emplace_back(
            std::to_string(anchorGraph.getAnchorId(v0)),
            std::to_string(anchorGraph.getAnchorId(v1)));
    }

    {
        std::ofstream gfaFile(gfaFileName);
        if(!gfaFile) {
            cout << "[SuperbubbleDetection] BubbleFinder comparison skipped: "
                "failed to open " << gfaFileName << " for writing." << endl;
            return;
        }

        gfaFile << "H\tVN:Z:1.0\n";
        for(const std::string& label: nodeLabels) {
            gfaFile << "S\t" << label << "\t*\n";
        }
        for(const auto& [u, v]: directedEdges) {
            gfaFile << "L\t" << u << "\t+\t" << v << "\t+\t0M\n";
        }
    }

    const std::string command =
        shellQuote(bubbleFinderBinary) +
        " directed-superbubbles --gfa-directed -g " +
        shellQuote(gfaFileName) +
        " -o " +
        shellQuote(outputFileName) +
        " -j " +
        std::to_string(bubbleFinderThreadCount) +
        " > " +
        shellQuote(stdoutLogFileName) +
        " 2> " +
        shellQuote(stderrLogFileName);

    cout << "[SuperbubbleDetection] Running BubbleFinder directed-superbubbles "
        "for AnchorGraph comparison." << endl;
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        cout << "[SuperbubbleDetection] BubbleFinder comparison failed with exit code "
            << errorCode << ". Command stderr: " << stderrLogFileName << endl;
        return;
    }

    std::set<BubbleEndpointPair> onoderaPairs;
    onoderaPairs.clear();
    for(const auto& p: onoderaSuperbubbles) {
        const auto entrance = std::to_string(anchorGraph.getAnchorId(p.first));
        const auto exit = std::to_string(anchorGraph.getAnchorId(p.second));
        onoderaPairs.insert(canonicalBubblePair(entrance, exit));
    }

    std::ifstream bubbleFinderOutput(outputFileName);
    if(!bubbleFinderOutput) {
        cout << "[SuperbubbleDetection] BubbleFinder comparison skipped: "
            "failed to open " << outputFileName << " for reading." << endl;
        return;
    }

    uint64_t declaredPairCount = 0;
    if(!(bubbleFinderOutput >> declaredPairCount)) {
        cout << "[SuperbubbleDetection] BubbleFinder comparison skipped: "
            "invalid output header in " << outputFileName << "." << endl;
        return;
    }

    std::set<BubbleEndpointPair> bubbleFinderPairs;
    std::vector<BubbleEndpointPair> bubbleFinderSamples;
    uint64_t parsedPairCount = 0;
    uint64_t ignoredAuxiliaryPairCount = 0;
    std::string u;
    std::string v;
    while(bubbleFinderOutput >> u >> v) {
        parsedPairCount++;
        u = normalizeBubbleEndpoint(u);
        v = normalizeBubbleEndpoint(v);
        if(nodeLabels.find(u) == nodeLabels.end() ||
           nodeLabels.find(v) == nodeLabels.end()) {
            ignoredAuxiliaryPairCount++;
            continue;
        }
        const auto pair = canonicalBubblePair(u, v);
        bubbleFinderPairs.insert(pair);
        if(bubbleFinderSamples.size() < bubbleFinderSampleCount) {
            bubbleFinderSamples.push_back(pair);
        }
    }

    uint64_t sharedPairCount = 0;
    for(const auto& pair: onoderaPairs) {
        if(bubbleFinderPairs.find(pair) != bubbleFinderPairs.end()) {
            sharedPairCount++;
        }
    }

    const uint64_t onoderaOnlyCount = onoderaPairs.size() - sharedPairCount;
    const uint64_t bubbleFinderOnlyCount = bubbleFinderPairs.size() - sharedPairCount;

    cout << "[SuperbubbleDetection] BubbleFinder comparison: "
        << "Onodera=" << onoderaPairs.size()
        << ", BubbleFinder=" << bubbleFinderPairs.size()
        << ", shared=" << sharedPairCount
        << ", onoderaOnly=" << onoderaOnlyCount
        << ", bubbleFinderOnly=" << bubbleFinderOnlyCount
        << ", ignoredAuxiliary=" << ignoredAuxiliaryPairCount;
    if(parsedPairCount != declaredPairCount) {
        cout << ", declaredPairs=" << declaredPairCount
            << ", parsedPairs=" << parsedPairCount;
    }
    cout << "." << endl;

    if(!bubbleFinderSamples.empty()) {
        cout << "[SuperbubbleDetection] BubbleFinder sample pairs:";
        for(const auto& p: bubbleFinderSamples) {
            cout << " " << p.first << "->" << p.second;
        }
        cout << endl;
    }
}

} // namespace




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
        // // Use SIMD-accelerated closed syncmers for initial marker generation (no filtering).
        // assembler.findMarkersSimdClosedSyncmers(
        //     threadCount,
        //     assemblerOptions.kmersOptions.k,
        //     assemblerOptions.kmersOptions.syncmerS);

        // Use SIMD-accelerated minimizers instead of closed syncmers.
        // For hifiasm-like behavior with k=w, use syncmerS parameter as window size.
        // Density ≈ 2/w (smaller w = denser sampling, larger w = sparser sampling)
        assembler.findMarkersSimdMinimizers(
            threadCount,
            assemblerOptions.kmersOptions.k,
            assemblerOptions.kmersOptions.k);  // Using kmer length as window size w

        // Compute histogram using the pre-calculated KmerIds.
        assembler.countKmersFromMarkerKmerIds(threadCount);
        
        // Retrieve peak and set thresholds.
        // Hifiasm hard-filters k-mers above max_kmer_cnt (default 2000) during
        // sketching. K-mers between highFreqThreshold (coveragePeak * 1.667)
        // and this cutoff are kept but downsampled per-streak during chaining.
        const uint64_t coveragePeak = assembler.assemblerInfo->kmerDistributionInfo.coveragePeak;
        const uint64_t minFreq = 2;
        const uint64_t maxFreq = uint64_t(
            assemblerOptions.overlapCandidatesOptions.invertedIndexMaxKmerCount);
        const bool removePalindromicKmers = true;
        uint64_t distinctKmerCount = 0;
        for(uint64_t bucketId=0; bucketId<assembler.kmerCounter->kmerIdFrequencies.size(); bucketId++) {
            distinctKmerCount += assembler.kmerCounter->kmerIdFrequencies[bucketId].size();
        }

        cout << "Analyzing " << distinctKmerCount << " distinct minimizer k-mers." << endl;
        cout << "Filtering minimizers: Peak coverage is " << coveragePeak << "." << endl;
        cout << "Keeping k-mers with frequency [" << minFreq << ", " << maxFreq << "]";
        if(removePalindromicKmers) {
            cout << " and excluding palindromic k-mers";
        }
        cout << "." << endl;

        [[maybe_unused]]
        auto writeReadMarkerGapDiagnostic = [&assembler](const string& label, ReadId readId) {
            const OrientedReadId oid(readId, 0);
            const auto read = assembler.getReads().getRead(readId);
            const auto readMarkers = (*assembler.markers)[oid.getValue()];
            const auto readKmerIds = (*assembler.markerKmerIds)[oid.getValue()];
            const uint64_t k = assembler.assemblerInfo->k;

            class GapInfo {
            public:
                string type;
                uint32_t begin = 0;
                uint32_t end = 0;
                uint64_t markerStartDistance = 0;
                uint64_t noMarkerBases = 0;
                uint64_t leftFrequency = 0;
                uint64_t rightFrequency = 0;
            };

            auto frequency = [&](KmerId kmerId) {
                const Kmer kmer(kmerId, k);
                const KmerId rcKmerId = kmer.reverseComplement(k).id(k);
                return assembler.kmerCounter->getFrequencyFast(min(kmerId, rcKmerId));
            };

            vector<GapInfo> gaps;
            if(readMarkers.empty()) {
                gaps.push_back(GapInfo{
                    "wholeRead",
                    0,
                    uint32_t(read.baseCount),
                    read.baseCount,
                    read.baseCount,
                    0,
                    0});
            } else {
                const uint32_t firstPosition = readMarkers.front().position;
                if(firstPosition != 0) {
                    gaps.push_back(GapInfo{
                        "prefix",
                        0,
                        firstPosition,
                        firstPosition,
                        firstPosition,
                        0,
                        frequency(readKmerIds.front())});
                }
                for(uint64_t i=1; i<readMarkers.size(); i++) {
                    const uint32_t previousPosition = readMarkers[i - 1].position;
                    const uint32_t nextPosition = readMarkers[i].position;
                    const uint64_t markerStartDistance = nextPosition > previousPosition ?
                        nextPosition - previousPosition : 0;
                    const uint64_t previousEnd = uint64_t(previousPosition) + k;
                    const uint64_t noMarkerBases = uint64_t(nextPosition) > previousEnd ?
                        uint64_t(nextPosition) - previousEnd : 0;
                    gaps.push_back(GapInfo{
                        "internal",
                        previousPosition,
                        nextPosition,
                        markerStartDistance,
                        noMarkerBases,
                        frequency(readKmerIds[i - 1]),
                        frequency(readKmerIds[i])});
                }
                const uint64_t lastEnd = uint64_t(readMarkers.back().position) + k;
                if(read.baseCount > lastEnd) {
                    gaps.push_back(GapInfo{
                        "suffix",
                        uint32_t(lastEnd),
                        uint32_t(read.baseCount),
                        read.baseCount - lastEnd,
                        read.baseCount - lastEnd,
                        frequency(readKmerIds.back()),
                        0});
                }
            }

            sort(gaps.begin(), gaps.end(),
                [](const GapInfo& a, const GapInfo& b) {
                    if(a.noMarkerBases != b.noMarkerBases) {
                        return a.noMarkerBases > b.noMarkerBases;
                    }
                    return a.markerStartDistance > b.markerStartDistance;
                });

            cout << timestamp << "[MarkerGapDiagnostic] " << label
                 << " readId=" << readId
                 << " readLength=" << read.baseCount
                 << " markerCount=" << readMarkers.size()
                 << endl;
            for(uint64_t i=0; i<min<uint64_t>(10, gaps.size()); i++) {
                const GapInfo& gap = gaps[i];
                cout << timestamp << "  rank=" << i
                     << " type=" << gap.type
                     << " begin=" << gap.begin
                     << " end=" << gap.end
                     << " markerStartDistance=" << gap.markerStartDistance
                     << " noMarkerBases=" << gap.noMarkerBases
                     << " leftFrequency=" << gap.leftFrequency
                     << " rightFrequency=" << gap.rightFrequency
                     << endl;
            }
        };

        // writeReadMarkerGapDiagnostic("beforeFrequencyFilter", ReadId(3729));
             
        // Prune the existing minimizer markers in-place.
        // applyKmerCountFilter keeps a marker only if:
        // - its canonical k-mer frequency is in the inclusive range [minFreq, maxFreq],
        // - and, by default, the k-mer is not palindromic/self-reverse-complementary.
        // The function rebuilds both markers and markerKmerIds from the pre-filtered
        // arrays, preserving only marker positions whose matching k-mer id passes.
        assembler.applyKmerCountFilter(minFreq, maxFreq, threadCount, removePalindromicKmers);
        
        // writeReadMarkerGapDiagnostic("afterFrequencyFilter", ReadId(3729));

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

        // Compute k-mer histogram to get coveragePeak (needed by phasing).
        // The SIMD path does this via countKmersFromMarkerKmerIds.
        assembler.countKmersFromMarkerKmerIds(threadCount);
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


    // ========================================================================
    // HIFIASM MAX_N_CHAIN CALCULATION (Per-Read Overlap Limiting)
    // ========================================================================
    // Reference: Hifiasm CommandLines.cpp:413, anchor.cpp:191-220
    //
    // Compute max_n_chain = max(MIN_N_CHAIN, hom_cov * high_factor)
    //
    // Purpose: Limit total number of overlaps kept per read across ALL partners
    // - Applied AFTER chaining all read pairs (not during DP chaining)
    // - Prevents memory explosion from reads with many partners
    // - Keeps top N overlaps by score, grouped by overlap type
    //
    // Hifiasm defaults:
    //   high_factor = 5.0 (CommandLines.cpp:271)
    //   MIN_N_CHAIN = 100 (CommandLines.h:28)
    //
    // Example: coverage = 30x → max(100, 30*5) = 150 overlaps per read
    //
    const uint64_t coveragePeak = assembler.assemblerInfo->kmerDistributionInfo.coveragePeak;
    const uint64_t maxChainLimit = std::max<uint64_t>(
        assemblerOptions.overlapCandidatesOptions.invertedIndexMinNChain,     // MIN_N_CHAIN = 100
        uint64_t(double(coveragePeak) * assemblerOptions.overlapCandidatesOptions.invertedIndexHighFactor + 0.499));  // hom_cov * 5.0

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
            threadCount
        );
    } else {
        // Inverted Index path: Discover candidate pairs via k-mer matches and chain them.
        assembler.chainAlignmentCandidates(
            assemblerOptions.overlapCandidatesOptions.driftRateTolerance,
            maxChainLimit,
            assemblerOptions.overlapCandidatesOptions,
            threadCount
        );
    }

    // // Experimental path: run FASTGA directly on the chained candidate spans
    // // and store sparse differences keyed by candidate index.
    // assembler.alignChainedCandidatesWithFastGA(threadCount);
    // return;

    // // Lightweight marker-chain materialization.
    // // The marker graph vertex builder needs alignmentData/compressedAlignments,
    // // but this prototype does not need projected banded/base alignments or evidence.
    // assembler.computeAlignmentDataFromChainedCandidatesOnly(
    //     assemblerOptions.alignOptions,
    //     threadCount);

    // Previous full evidence path. Use this instead of the lightweight path
    // when base-level projected alignments and SNP/indel evidence are needed.
    assembler.computeBaseAlignmentsAndStore(
        assemblerOptions.alignOptions,
        threadCount);

    // For http server and debugging/development purposes, generate an exhaustive table of candidates.
    // This can be done after alignment computation (it depends only on the candidate list).
    assembler.computeCandidateTable();

    assembler.phaseOverlaps(threadCount);

    assembler.createReadGraphFromPhasingCisOverlaps();
    

    // Set min and max marker graph vertex coverage thresholds.
    const uint64_t minAnchorCoverage = 8;
    const uint64_t maxAnchorCoverage = 5 * coveragePeak;

    // Build marker graph vertices using transitive alignments collapse.
    assembler.createMarkerGraphVertices(
        minAnchorCoverage,                                              // minVertexCoverage
        maxAnchorCoverage,                                              // maxVertexCoverage
        0,                                              // minVertexCoveragePerStrand
        false,                                          // allowDuplicateMarkers
        std::numeric_limits<double>::signaling_NaN(),   // unused (minVertexCoverage != 0)
        invalid<uint64_t>,                              // unused (minVertexCoverage != 0)
        threadCount);
    
    // // Filter marker graph vertices whose marker k-mers are short-period repeats (including homopolymers).
    // // This reduces unreliable anchors and artifacts in repetitive regions.
    // assembler.filterMarkerGraphVerticesByRepeatKmers(threadCount);

    // // Filter marker graph vertices whose marker k-mers have low sequence complexity
    // // (too few distinct sub-k-mers of lengths 1, 2, 3, ...).
    // assembler.filterMarkerGraphVerticesByDistinctSubkmerCount(threadCount);

    // Find the reverse complement of each marker graph vertex.
    // We need the reverse complement vertices to be populated for anchor generation.
    assembler.findMarkerGraphReverseComplementVertices(threadCount);


    // const uint64_t minPrimaryCoverage = assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage;;
    // const uint64_t maxPrimaryCoverage = assemblerOptions.assemblyOptions.mode3Options.maxAnchorCoverage;;
    cout << "Using: minAnchorCoverage = " << minAnchorCoverage <<
        ", maxAnchorCoverage = " << maxAnchorCoverage << endl;


    const MappedMemoryOwner shasta2Owner = assembler.shasta2MappedMemoryOwner();
    
    // Create Shasta2Anchors with vertex splitting.
    // Splits marker graph vertices whose reads were merged by transitive
    // closure but lack direct pairwise overlaps in the read graph.
    assembler.shasta2Anchors = createShasta2AnchorsFromSplitVertices(
            shasta2Owner,
            assembler.getReads(),
            assembler.assemblerInfo->k,
            *assembler.markers,
            assembler.markerGraph,
            assembler.readGraph,
            threadCount,
            minAnchorCoverage,
            maxAnchorCoverage);
    auto& shasta2Anchors = assembler.shasta2Anchors;

    // double kmerDensity = 1.0;
    // cout << timestamp << "Filtering Shasta2Anchors with Shasta2 hashed k-mer checker..." << endl;
    // // shasta2Anchors->filterByShasta2HashedKmerChecker(
    // //     assemblerOptions.kmersOptions.probability);
    // shasta2Anchors->filterByShasta2HashedKmerChecker(
    //     kmerDensity);

    const string externalAnchorsName =
        std::filesystem::absolute("Shasta2ExternalAnchors").string();
    cout << timestamp << "Writing Shasta2 external anchors to "
         << externalAnchorsName << "..." << endl;
    const uint64_t exportedExternalAnchorCount =
        shasta2Anchors->writeExternalAnchors(externalAnchorsName);
    cout << timestamp << "Wrote " << exportedExternalAnchorCount
         << " external anchors for Shasta2. Use --external-anchors-name "
         << externalAnchorsName << endl;

    // Compute journeys.
    cout << timestamp << "Creating Shasta2Journeys..." << endl;
    assembler.shasta2Journeys = make_shared<Shasta2Journeys>(
        2 * assembler.getReads().readCount(),
        shasta2Anchors,
        threadCount,
        shasta2Owner);
    auto& shasta2Journeys = assembler.shasta2Journeys;

    // Create the Shasta2AnchorGraph.
    const uint64_t minEdgeCoverage = 2;
    cout << timestamp << "Creating Shasta2AnchorGraph..." << endl;
    assembler.shasta2AnchorGraph = make_shared<Shasta2AnchorGraph>(
        *shasta2Anchors,
        *shasta2Journeys,
        minEdgeCoverage,
        threadCount);
    auto& shasta2AnchorGraph = assembler.shasta2AnchorGraph;

    // Save the pre-transitive-reduction Shasta2 anchor graph so the HTTP server
    // can load and visualize it even when we return before later assembly stages.
    shasta2AnchorGraph->saveAnchorGraph("Shasta2AnchorGraph");
    shasta2AnchorGraph->writeGfa("Shasta2AnchorGraph.gfa");

    return;













    

    // The marker graph vertex builder iterates readGraph edges. For this
    // diagnostic prototype, keep all chained alignments and let marker-graph
    // coverage/repeat/complexity filters do the pruning.
    assembler.createReadGraphAllAlignments();

    // // Delete overlaps where one read is contained in the other.
    // assembler.deleteContainmentOverlaps(threadCount);

    // // // Delete internal overlaps (excessive overhangs or too short).
    // // assembler.deleteInternalOverlaps(500, 0.8, 50, threadCount);

    

    // assembler.debugDumpSnpSitesForRead(0, 3);
    // return;


    // coveragePeak already defined above (line ~979).
    const uint64_t minFreq = 8;
    const uint64_t maxFreq = 5 * coveragePeak;


    // Build marker graph vertices needed by performHifiasmECParityWithMarkerGraph.
    assembler.createMarkerGraphVertices(
        minFreq,                                              // minVertexCoverage
        maxFreq,                                              // maxVertexCoverage
        0,                                              // minVertexCoveragePerStrand
        false,                                          // allowDuplicateMarkers
        std::numeric_limits<double>::signaling_NaN(),   // unused (minVertexCoverage != 0)
        invalid<uint64_t>,                              // unused (minVertexCoverage != 0)
        threadCount);
    // assembler.filterMarkerGraphVerticesByRepeatKmers(threadCount);
    // assembler.filterMarkerGraphVerticesByDistinctSubkmerCount(threadCount);
    assembler.findMarkerGraphReverseComplementVertices(threadCount);


    {
    // Same Shasta2 anchor/journey coverage as the main assembly path (see shasta2 block
    // later in this file): reuse these objects when wiring Theseus to Shasta2.
    const uint64_t minPrimaryCoverage = 8;
    const uint64_t maxPrimaryCoverage = 5 * coveragePeak;
    const MappedMemoryOwner shasta2OwnerEarly = assembler.shasta2MappedMemoryOwner();
    cout << timestamp << "Creating Shasta2Anchors for Theseus read-window prototype..." << endl;
    assembler.shasta2Anchors = make_shared<Shasta2Anchors>(
        shasta2OwnerEarly,
        assembler.getReads(),
        assembler.assemblerInfo->k,
        *assembler.markers,
        assembler.markerGraph,
        threadCount,
        minPrimaryCoverage,
        maxPrimaryCoverage);

    
    cout << timestamp << "Creating Shasta2Journeys for Theseus read-window prototype..." << endl;
    assembler.shasta2Journeys = make_shared<Shasta2Journeys>(
        2 * assembler.getReads().readCount(),
        assembler.shasta2Anchors,
        threadCount,
        shasta2OwnerEarly);

    // Structural scaffold before het detection / phasing: journey co-read CSR (marker graph).
    assembler.computeStrand0JourneyCoReadsTable();

    // Diagnostic prototype: partition Shasta2 anchor journeys into windows, then later
    // run one Theseus MSA per anchor-window interval.
    assembler.computeTheseusReadWindowMSAPrototype(
        assembler.shasta2Anchors,
        assembler.shasta2Journeys,
        threadCount);
    // assembler.computeTheseusMarkerGraphMSAPrototype(
    //     std::numeric_limits<uint64_t>::max(),    // maxAnchorPairs
    //     std::numeric_limits<uint64_t>::max(),    // maxReadsPerPair
    //     threadCount);
    // assembler.computeTheseusTargetBackboneMSAPrototype(
    //     12800,    // maxReads
    //     threadCount);
    return;
    }

    // Run marker-graph-projected EC parity (updates delete flags on alignments).
    assembler.performHifiasmECParityWithMarkerGraph(threadCount);

    // Print het sites for read 0-0 using the surviving (non-deleted) candidates.
    assembler.debugPrintHetSitesForRead(0);

    return;






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
    // const bool useGlobalSiteEcParity = (::getenv("DINARA_USE_GLOBAL_SITE_EC") != nullptr);
    const bool useGlobalSiteEcParity = false;
    if (useGlobalSiteEcParity) {
        cout << timestamp << "Using experimental global-site EC parity path." << endl;
        assembler.performGlobalSiteECParity(threadCount);
    } else {
        assembler.performHifiasmECParity(threadCount);
    }

    // return;



    // =========================================================================
    // Read graph construction — two alternatives:
    //   (A) All alignments: skip EC parity filtering entirely; defer to
    //       marker graph vertex coverage thresholds (minCoverage / maxCoverage).
    //   (B) Phased alignments: keep only overlaps where BOTH sides passed the
    //       phasing step of performHifiasmECParity (DeleteReasonPhase not set
    //       on either side).
    // =========================================================================
    // assembler.createReadGraphAllAlignments();       // (A)
    assembler.createReadGraphFromEcParityCisOverlaps(); // (B)

    

    // // Create the bidirectional read graph.
    // // This is built from the same alignments that are in the ReadGraph
    // // (isInReadGraph == 1), but stores one vertex per physical read and
    // // one edge per alignment (no strand doubling).
    // assembler.createBidirectionalReadGraph();

    // return;

    // // Clean the BRG using string-graph-style operations:
    // // Step 10: transitive reduction + Step 11: tip cutting.
    // // This removes redundant edges while operating on the undirected BRG
    // // via a temporary directed-arc view derived from alignment coordinates.
    // assembler.cleanBidirectionalReadGraphInitial(
    //     /*gapFuzz*/1000,
    //     /*maxShortTipReads*/3);

    // // Global mismatch-site diagnostics and export are expensive and intended for debugging.
    // // Keep them off by default in production runs to preserve assembly throughput.
    // const bool runGlobalHetDiagnostics = true;
    // if (runGlobalHetDiagnostics) {
    // // Global mismatch sites + full per-allele member lists using only readGraph overlaps.
    // // This is the fastest way to approximate "pileup across all reads" without a reference.
    // {
    //     const OrientedReadId focalOrientedReadId(ReadId(24347), 1);
    //     const ReadId focalReadId = focalOrientedReadId.getReadId();
    //     const Strand focalReadStrand = focalOrientedReadId.getStrand();
    //     const string focalReadLabel =
    //         to_string(uint64_t(focalReadId)) + "-" + to_string(uint64_t(focalReadStrand));
    //     const Reads& reads = assembler.getReads();
    //     const uint32_t focalReadLength = uint32_t(reads.getRead(focalReadId).baseCount);

    //     const auto clusters = assembler.clusterMismatchingPositionsIntoGlobalHetSitesReachableFromRead(
    //         focalReadId,
    //         assemblerOptions.alignOptions,
    //         threadCount,
    //         0,      // maxReadsToProcess
    //         0,      // maxAlignmentsToProcess
    //         false,  // includeDeletedAlignments
    //         true    // readGraphOnly
    //     );
    //     const uint32_t clusterSiteCount = uint32_t(
    //         clusters.clusterMemberOffsets.empty() ? 0 : (clusters.clusterMemberOffsets.size() - 1));

    //     static const char baseToAscii[] = {'A', 'C', 'G', 'T'};

    //     // Collect mismatch-defined sites that explicitly involve focalReadId (it has a mismatch at that site).
    //     struct FocalMismatchSite {
    //         uint32_t focalPos = 0;
    //         uint32_t siteId = 0;
    //     };
    //     vector<FocalMismatchSite> focalMismatchSites;
    //     focalMismatchSites.reserve(clusterSiteCount);
    //     for (size_t siteId = 0; siteId + 1 < clusters.clusterMemberOffsets.size(); siteId++) {
    //         const uint64_t begin = clusters.clusterMemberOffsets[siteId];
    //         const uint64_t end = clusters.clusterMemberOffsets[siteId + 1];
    //         uint32_t bestPos = std::numeric_limits<uint32_t>::max();
    //         for (uint64_t i = begin; i < end; i++) {
    //             const auto& node = clusters.nodes[clusters.clusterMembers[i]];
    //             if (node.first == focalReadId) {
    //                 bestPos = std::min(bestPos, node.second);
    //             }
    //         }
    //         if (bestPos != std::numeric_limits<uint32_t>::max()) {
    //             if (focalReadStrand == 1 && focalReadLength > 0 && bestPos < focalReadLength) {
    //                 bestPos = (focalReadLength - 1U) - bestPos;
    //             }
    //             focalMismatchSites.push_back(FocalMismatchSite{bestPos, uint32_t(siteId)});
    //         }
    //     }
    //     sort(focalMismatchSites.begin(), focalMismatchSites.end(),
    //         [](const FocalMismatchSite& a, const FocalMismatchSite& b) {
    //             if (a.focalPos != b.focalPos) {
    //                 return a.focalPos < b.focalPos;
    //             }
    //             return a.siteId < b.siteId;
    //         });
    //     cout << timestamp << "Read" << focalReadLabel
    //          << " mismatch sites (readGraph clusters): " << focalMismatchSites.size() << endl;

    //     const auto propagationStart = std::chrono::steady_clock::now();
    //     cout << timestamp << "GlobalHetSite member propagation (readGraph): starting..." << endl;
    //     const auto members = assembler.computeGlobalHetSiteAlleleMembersUsingReadGraph(
    //         clusters,
    //         assemblerOptions.alignOptions,
    //         0,      // maxPendingTasks
    //         false,  // includeDeletedAlignments
    //         focalReadId
    //     );
    //     const auto propagationSeconds = std::chrono::duration_cast<std::chrono::seconds>(
    //         std::chrono::steady_clock::now() - propagationStart).count();

    //     cout << timestamp << "GlobalHetSite member propagation (readGraph): sites=" << clusterSiteCount
    //          << " propagatedAssignments=" << members.propagatedAssignments
    //          << " mappingHoles=" << members.mappingHoles
    //          << " mappingConflicts=" << members.mappingConflicts
    //          << " elapsedSec=" << propagationSeconds
    //          << endl;

    //     // Compute a consistent strand assignment for reads reachable from the focal read in the read graph.
    //     // This lets us export member positions (and alleles) in a single, focal-oriented coordinate frame.
    //     uint64_t strandConflicts = 0;
    //     vector<int8_t> strandByRead = assembler.computeReadGraphStrandsFromSeed(
    //         focalReadId,
    //         strandConflicts,
    //         false // includeDeletedAlignments
    //     );
    //     if (focalReadStrand == 1) {
    //         for (int8_t& v : strandByRead) {
    //             if (v != -1) {
    //                 v = int8_t(v ^ 1);
    //             }
    //         }
    //     }
    //     {
    //         uint64_t assigned = 0;
    //         for (const int8_t v : strandByRead) {
    //             if (v != -1) {
    //                 assigned++;
    //             }
    //         }
    //         cout << timestamp << "ReadGraph strand assignment: assigned=" << assigned
    //              << " conflicts=" << strandConflicts << endl;
    //     }

    //     static const uint8_t complementBase[4] = {3, 2, 1, 0};
    //     const auto orientedMembers = assembler.orientGlobalHetSiteAlleleMembers(members, strandByRead);

    //     // Spot-check: verify that a few sites involving the focal read have self-consistent positions
    //     // under a readGraph-only multi-source traversal seeded from the mismatch members.
    //     {
    //         const size_t checkCount = std::min<size_t>(3, focalMismatchSites.size());
    //         for (size_t i = 0; i < checkCount; i++) {
    //             const uint32_t siteId = focalMismatchSites[i].siteId;
    //             const auto stats = assembler.debugVerifyGlobalHetSitePositionsUsingReadGraph(
    //                 clusters,
    //                 members,
    //                 siteId,
    //                 assemblerOptions.alignOptions,
    //                 20000,   // maxNodesToVisit
    //                 200000,  // maxAlignmentsToScan
    //                 false    // includeDeletedAlignments
    //             );
    //             cout << timestamp << "GlobalHetSite verify: siteId=" << siteId
    //                  << " expected=" << stats.expectedMembers
    //                  << " reached=" << stats.reachedMembers
    //                  << " checked=" << stats.checkedMappings
    //                  << " mismatched=" << stats.mismatchedPositions
    //                  << " holes=" << stats.mappingHoles
    //                  << " fails=" << stats.mappingFailures
    //                  << " hitNodeLimit=" << stats.hitNodeLimit
    //                  << " hitAlignmentLimit=" << stats.hitAlignmentLimit
    //                  << endl;
    //         }
    //     }

    //     const uint32_t siteCount = uint32_t(members.offsets.size());

    //     // Precompute mismatch-member counts and per-allele mismatch counts once per site.
    //     vector<array<uint32_t, 4> > mismatchCountsForward(siteCount, array<uint32_t, 4>{0, 0, 0, 0});
    //     vector<array<uint32_t, 4> > mismatchCountsOriented(siteCount, array<uint32_t, 4>{0, 0, 0, 0});
    //     vector<uint64_t> mismatchMembersBySite(siteCount, 0);
    //     for (uint32_t siteId = 0; siteId < siteCount && (siteId + 1) < clusters.clusterMemberOffsets.size(); siteId++) {
    //         const uint64_t begin = clusters.clusterMemberOffsets[siteId];
    //         const uint64_t end = clusters.clusterMemberOffsets[siteId + 1];
    //         mismatchMembersBySite[siteId] = end - begin;
    //         for (uint64_t j = begin; j < end; j++) {
    //             const auto& node = clusters.nodes[clusters.clusterMembers[j]];
    //             uint8_t b = reads.getOrientedReadBase(OrientedReadId(node.first, 0), node.second).value;
    //             if (b >= 4) {
    //                 continue;
    //             }
    //             mismatchCountsForward[siteId][b]++;
    //             const int8_t s = (uint64_t(node.first) < strandByRead.size()) ? strandByRead[uint64_t(node.first)] : int8_t(-1);
    //             if (s == 1) {
    //                 b = complementBase[b];
    //             }
    //             mismatchCountsOriented[siteId][b]++;
    //         }
    //     }

    //     // Precompute oriented support counts and total members per site once.
    //     // These are reused in filtering, printing, and export.
    //     vector<array<uint64_t, 4> > orientedSiteCounts(siteCount, array<uint64_t, 4>{0, 0, 0, 0});
    //     vector<uint64_t> orientedSiteMembers(siteCount, 0);
    //     const uint32_t orientedCount = uint32_t(orientedMembers.offsets.size());
    //     for (uint32_t siteId = 0; siteId < siteCount && siteId < orientedCount; siteId++) {
    //         const auto& off = orientedMembers.offsets[siteId];
    //         for (int allele = 0; allele < 4; allele++) {
    //             orientedSiteCounts[siteId][allele] = off[allele + 1] - off[allele];
    //             orientedSiteMembers[siteId] += orientedSiteCounts[siteId][allele];
    //         }
    //     }

    //     // Keep only robust multiallelic sites: at least 2 alleles with support >= 3.
    //     static constexpr uint32_t minAlleleSupportForExport = 3;
    //     static constexpr uint32_t minAlleleCountForExport = 2;
    //     vector<uint8_t> sitePassesMultiallelic(siteCount, 0);
    //     for (uint32_t siteId = 0; siteId < siteCount; siteId++) {
    //         uint32_t supportedAlleles = 0;
    //         for (int allele = 0; allele < 4; allele++) {
    //             if (orientedSiteCounts[siteId][allele] >= minAlleleSupportForExport) {
    //                 supportedAlleles++;
    //             }
    //         }
    //         sitePassesMultiallelic[siteId] = uint8_t(supportedAlleles >= minAlleleCountForExport);
    //     }

    //     const auto readIndex = assembler.buildFilteredGlobalHetSiteReadIndex(
    //         members,
    //         minAlleleSupportForExport,
    //         minAlleleCountForExport
    //     );
    //     const uint32_t invalidPos = std::numeric_limits<uint32_t>::max();
    //     const vector<Assembler::GlobalHetSiteReadIndex::ReadSite> emptyFocalReadSites;
    //     const auto& focalReadSites =
    //         (uint64_t(focalReadId) < readIndex.sitesByRead.size()) ?
    //         readIndex.sitesByRead[uint64_t(focalReadId)] :
    //         emptyFocalReadSites;
    //     vector<uint32_t> focalReadPosBySite(siteCount, invalidPos);
    //     vector<char> focalReadAlleleBySite(siteCount, '?');
    //     for (const auto& s : focalReadSites) {
    //         if (s.siteId < siteCount) {
    //             uint32_t pos = s.readPosition;
    //             uint8_t allele = s.allele;
    //             if (focalReadStrand == 1 && focalReadLength > 0 && pos < focalReadLength) {
    //                 pos = (focalReadLength - 1U) - pos;
    //                 allele = complementBase[allele];
    //             }
    //             focalReadPosBySite[s.siteId] = pos;
    //             focalReadAlleleBySite[s.siteId] = baseToAscii[allele];
    //         }
    //     }

    //     vector<Assembler::GlobalHetSiteReadIndex::ReadSite> filteredFocalReadSites;
    //     filteredFocalReadSites.reserve(focalReadSites.size());
    //     for (const auto& s : focalReadSites) {
    //         if (s.siteId < sitePassesMultiallelic.size() && sitePassesMultiallelic[s.siteId]) {
    //             filteredFocalReadSites.push_back(s);
    //         }
    //     }

    //     cout << timestamp << "Read" << focalReadLabel
    //          << " projected global het sites after multiallelic filter: "
    //          << filteredFocalReadSites.size() << " / " << focalReadSites.size()
    //          << " (need >= " << minAlleleCountForExport << " alleles with support >= "
    //          << minAlleleSupportForExport << ")" << endl;

    //     // Print 10 sites involving focalReadId.
    //     const size_t toPrint = std::min<size_t>(10, filteredFocalReadSites.size());
    //     for (size_t i = 0; i < toPrint; i++) {
    //         const uint32_t siteId = filteredFocalReadSites[i].siteId;
    //         const uint32_t focalPos =
    //             (siteId < focalReadPosBySite.size()) ? focalReadPosBySite[siteId] : invalidPos;
    //         const char focalAllele =
    //             (siteId < focalReadAlleleBySite.size()) ? focalReadAlleleBySite[siteId] : '?';
    //         const uint64_t mismatchMembers =
    //             (siteId < mismatchMembersBySite.size()) ? mismatchMembersBySite[siteId] : 0;
    //         const auto mismatchCounts =
    //             (siteId < mismatchCountsForward.size()) ?
    //             mismatchCountsForward[siteId] :
    //             std::array<uint32_t, 4>{0, 0, 0, 0};
    //         const auto siteCounts =
    //             (siteId < orientedSiteCounts.size()) ?
    //             orientedSiteCounts[siteId] :
    //             std::array<uint64_t, 4>{0, 0, 0, 0};
    //         const uint64_t siteMembers = (siteId < orientedSiteMembers.size()) ? orientedSiteMembers[siteId] : 0;

    //         cout << timestamp
    //              << "GlobalHetSite[" << i << "]"
    //              << " read" << focalReadLabel << "Pos=" << focalPos
    //              << " read" << focalReadLabel << "Allele=" << focalAllele
    //              << " mismatchMembers=" << mismatchMembers
    //              << " mismatchCounts(A,C,G,T)=(" << mismatchCounts[0] << "," << mismatchCounts[1] << "," << mismatchCounts[2] << "," << mismatchCounts[3] << ")"
    //              << " siteMembers=" << siteMembers
    //              << " siteCounts(A,C,G,T)=(" << siteCounts[0] << "," << siteCounts[1] << "," << siteCounts[2] << "," << siteCounts[3] << ")"
    //              << " members={";

    //         // Show up to 8 members per allele.
    //         for (int allele = 0; allele < 4; allele++) {
    //             const uint64_t b0 = orientedMembers.offsets[siteId][allele];
    //             const uint64_t b1 = orientedMembers.offsets[siteId][allele + 1];
    //             const uint64_t show = std::min<uint64_t>(b1 - b0, 8);
    //             if (show == 0) {
    //                 continue;
    //             }
    //             cout << baseToAscii[allele] << ":{";
    //             for (uint64_t k = 0; k < show; k++) {
    //                 const auto& m = orientedMembers.members[b0 + k];
    //                 cout << m.orientedReadId.getReadId()
    //                      << (m.orientedReadId.getStrand() == 1 ? "rc" : "fw")
    //                      << "-" << m.position;
    //                 if (k + 1 < show) {
    //                     cout << ",";
    //                 }
    //             }
    //             if ((b1 - b0) > show) {
    //                 cout << ",...";
    //             }
    //             cout << "}";
    //         }
    //         cout << "}" << endl;
    //     }

    //     // Export all SNP sites involving read 0 (summary + full per-allele member list).
    //     {
    //         const string summaryFileName = "Read" + focalReadLabel + "GlobalHetSitesSummary.tsv";
    //         const string membersFileName = "Read" + focalReadLabel + "GlobalHetSitesMembers.tsv";
    //         std::ofstream summary(summaryFileName);
    //         std::ofstream membersOut(membersFileName);
    //         if (!summary || !membersOut) {
    //             cout << timestamp << "Failed to open export files for read " << focalReadLabel << "." << endl;
    //         } else {
    //             summary << "siteId\treadPos\tmismatchMembers\tmismatchA\tmismatchC\tmismatchG\tmismatchT"
    //                     << "\tsiteMembers\tsiteA\tsiteC\tsiteG\tsiteT\treadAllele\n";
    //             membersOut << "siteId\treadPos0\treadAllele\treadId\treadStrand\tposition0\tpositionForward0\treadLength\n";

    //             // Prefer the mismatch-defined sites for "SNP sites of read0".
    //             // If there are none, fall back to the propagated membership list.
    //             const bool useMismatchSites = !focalMismatchSites.empty();
    //             const size_t exportCount = useMismatchSites ? focalMismatchSites.size() : filteredFocalReadSites.size();
    //             vector<uint32_t> readLengths(reads.readCount(), 0);
    //             for (uint64_t iRead = 0; iRead < reads.readCount(); iRead++) {
    //                 const ReadId rid = ReadId(iRead);
    //                 readLengths[iRead] = uint32_t(reads.getRead(rid).baseCount);
    //             }
    //             size_t exportedCount = 0;
    //             size_t filteredOutCount = 0;
    //             for (size_t idx = 0; idx < exportCount; idx++) {
    //                 const uint32_t siteId = useMismatchSites ? focalMismatchSites[idx].siteId : filteredFocalReadSites[idx].siteId;
    //                 if (siteId >= readIndex.sitePassesFilter.size() || readIndex.sitePassesFilter[siteId] == 0) {
    //                     filteredOutCount++;
    //                     continue;
    //                 }
    //                 if (siteId >= sitePassesMultiallelic.size() || sitePassesMultiallelic[siteId] == 0) {
    //                     filteredOutCount++;
    //                     continue;
    //                 }
    //                 if (siteId >= focalReadPosBySite.size() || focalReadPosBySite[siteId] == invalidPos) {
    //                     // Keep per-read-consistency filtering strict for DP-ready exports.
    //                     filteredOutCount++;
    //                     continue;
    //                 }
    //                 const uint32_t readPos = focalReadPosBySite[siteId];
    //                 const char readAllele = focalReadAlleleBySite[siteId];

    //                 const uint64_t mismatchMembers =
    //                     (siteId < mismatchMembersBySite.size()) ? mismatchMembersBySite[siteId] : 0;
    //                 const auto mismatchCounts =
    //                     (siteId < mismatchCountsOriented.size()) ?
    //                     mismatchCountsOriented[siteId] :
    //                     std::array<uint32_t, 4>{0, 0, 0, 0};
    //                 const auto siteCounts =
    //                     (siteId < orientedSiteCounts.size()) ?
    //                     orientedSiteCounts[siteId] :
    //                     std::array<uint64_t, 4>{0, 0, 0, 0};

    //                 // Export members using the focal-oriented coordinate frame (strandByRead),
    //                 // plus the original forward coordinates for debugging.
    //                 const uint32_t readPos0 = readPos;
    //                 // Export members in oriented coordinates (position0/1) consistent with readStrand,
    //                 // plus forward positions for debugging.
    //                 for (int allele = 0; allele < 4; allele++) {
    //                     const uint64_t b0 = orientedMembers.offsets[siteId][allele];
    //                     const uint64_t b1 = orientedMembers.offsets[siteId][allele + 1];
    //                     for (uint64_t k = b0; k < b1; k++) {
    //                         const auto& om = orientedMembers.members[k];
    //                         const ReadId rid = om.orientedReadId.getReadId();
    //                         const Strand strand = om.orientedReadId.getStrand();
    //                         const uint32_t posOriented0 = om.position;
    //                         const uint32_t len = (uint64_t(rid) < readLengths.size()) ? readLengths[uint64_t(rid)] : 0;
    //                         if (len == 0 || posOriented0 >= len) {
    //                             continue;
    //                         }
    //                         const uint32_t posFwd0 = (strand == 1) ? ((len - 1U) - posOriented0) : posOriented0;

    //                         // Allele char is already in oriented frame by construction (bucketed by allele).
    //                         const char alleleChar = baseToAscii[allele];

    //                         membersOut << siteId << "\t" << readPos0
    //                                    << "\t" << alleleChar
    //                                    << "\t" << rid
    //                                    << "\t" << int(strand)
    //                                    << "\t" << posOriented0
    //                                    << "\t" << posFwd0
    //                                    << "\t" << len
    //                                    << "\n";
    //                     }
    //                 }

    //                 const uint64_t siteMembers = (siteId < orientedSiteMembers.size()) ? orientedSiteMembers[siteId] : 0;
    //                 summary << siteId << "\t" << readPos
    //                         << "\t" << mismatchMembers
    //                         << "\t" << mismatchCounts[0] << "\t" << mismatchCounts[1] << "\t" << mismatchCounts[2] << "\t" << mismatchCounts[3]
    //                         << "\t" << siteMembers
    //                         << "\t" << siteCounts[0] << "\t" << siteCounts[1] << "\t" << siteCounts[2] << "\t" << siteCounts[3]
    //                         << "\t" << readAllele
    //                         << "\n";
    //                 exportedCount++;
    //             }

    //             cout << timestamp << "Wrote read" << focalReadLabel << " global het sites to " << summaryFileName
    //                  << " and " << membersFileName
    //                  << " (sites=" << exportedCount
    //                  << ", filteredOut=" << filteredOutCount
    //                  << ", criteria: >= " << minAlleleCountForExport
    //                  << " alleles with support >= " << minAlleleSupportForExport
    //                  << ")." << endl;
    //         }
    //     }
    // }
    // } else {
    //     cout << timestamp << "Skipping global-het diagnostics/export. "
    //          << "Set DINARA_ENABLE_GLOBAL_HET_DEBUG=1 to enable." << endl;
    // }

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

    // // Snapshot the broad keep-set used for marker-graph collapse.
    // std::vector<bool> keepForMarkerGraph(assembler.alignmentData.size(), false);
    // for (uint64_t i = 0; i < keepForMarkerGraph.size(); ++i) {
    //     keepForMarkerGraph[i] = (assembler.alignmentData[i].info.isInReadGraph != 0);
    // }

    // Mode 3 assembly requires reads in raw representation (not RLE).
    DINARA_ASSERT(assemblerOptions.readsOptions.representation == 0);

    // The marker length must be even.
    DINARA_ASSERT((assembler.assemblerInfo->k %2) == 0);


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

    // Filter marker graph vertices whose marker k-mers have low sequence complexity
    // (too few distinct sub-k-mers of lengths 1, 2, 3, ...).
    assembler.filterMarkerGraphVerticesByDistinctSubkmerCount(threadCount);

    // Find the reverse complement of each marker graph vertex.
    // We need the reverse complement vertices to be populated for Mode 3 anchor generation.
    assembler.findMarkerGraphReverseComplementVertices(threadCount);

    // // Create edges of the marker graph.
    // assembler.createMarkerGraphEdges(threadCount);
    // assembler.findMarkerGraphReverseComplementEdges(threadCount);

    if(assemblerOptions.markerGraphOptions.writeVertexCoverageHistogram) {
        cout << timestamp << "Writing marker graph vertex coverage histogram to " <<
            assemblerOptions.markerGraphOptions.vertexCoverageHistogramFileName << "." << endl;
        assembler.markerGraph.writeVertexCoverageHistogram(
            assemblerOptions.markerGraphOptions.vertexCoverageHistogramFileName,
            assemblerOptions.markerGraphOptions.vertexCoverageHistogramCanonicalOnly);
    }

    {
    const uint64_t minPrimaryCoverage = 2;
    const uint64_t maxPrimaryCoverage = std::numeric_limits<uint64_t>::max();
    // const uint64_t minPrimaryCoverage = assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage;;
    // const uint64_t maxPrimaryCoverage = assemblerOptions.assemblyOptions.mode3Options.maxAnchorCoverage;;
    cout << "Using: minAnchorCoverage = " << minPrimaryCoverage <<
        ", maxAnchorCoverage = " << maxPrimaryCoverage << endl;

    // // // Declare anchors pointer here to avoid scope issues
    // shared_ptr<mode3::Anchors> anchors;
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

    // // Compute oriented read journeys.
    // anchors->computeJourneys(threadCount);

    // assembler.mode3Assembly(threadCount, anchors, assemblerOptions.assemblyOptions.mode3Options, false);

    



    const MappedMemoryOwner shasta2Owner = assembler.shasta2MappedMemoryOwner();
    
    // Create Shasta2Anchors
    // We use the markerGraph structure to define anchors.
    assembler.shasta2Anchors = make_shared<Shasta2Anchors>(
            shasta2Owner,
            assembler.getReads(),
            assembler.assemblerInfo->k,
            *assembler.markers,
            assembler.markerGraph,
            threadCount,
            minPrimaryCoverage,
            maxPrimaryCoverage);
    auto& shasta2Anchors = assembler.shasta2Anchors;

    // double kmerDensity = 1.0;
    // cout << timestamp << "Filtering Shasta2Anchors with Shasta2 hashed k-mer checker..." << endl;
    // // shasta2Anchors->filterByShasta2HashedKmerChecker(
    // //     assemblerOptions.kmersOptions.probability);
    // shasta2Anchors->filterByShasta2HashedKmerChecker(
    //     kmerDensity);

    const string externalAnchorsName =
        std::filesystem::absolute("Shasta2ExternalAnchors").string();
    cout << timestamp << "Writing Shasta2 external anchors to "
         << externalAnchorsName << "..." << endl;
    const uint64_t exportedExternalAnchorCount =
        shasta2Anchors->writeExternalAnchors(externalAnchorsName);
    cout << timestamp << "Wrote " << exportedExternalAnchorCount
         << " external anchors for Shasta2. Use --external-anchors-name "
         << externalAnchorsName << endl;

    // Compute journeys.
    cout << timestamp << "Creating Shasta2Journeys..." << endl;
    assembler.shasta2Journeys = make_shared<Shasta2Journeys>(
        2 * assembler.getReads().readCount(),
        shasta2Anchors,
        threadCount,
        shasta2Owner);
    auto& shasta2Journeys = assembler.shasta2Journeys;

    // // --- Remove overlapping anchors from journeys ---
    // // Two consecutive anchors on the same oriented read overlap when the base
    // // position of anchor i+1 is less than position(anchor_i) + k, i.e. the
    // // k-mers share bases. We greedily keep the first anchor of any overlapping
    // // pair (anchors are already in ordinal order so the first has the smaller
    // // position). The result is stored as a plain vector so the memory-mapped
    // // journeys are not modified.
    // cout << timestamp << "Removing overlapping anchors from journeys..." << endl;
    // {
    //     const uint64_t k = assembler.assemblerInfo->k;
    //     const auto& mkrs = *assembler.markers;
    //     const uint64_t orientedReadCount = 2 * assembler.getReads().readCount();

    //     assembler.shasta2LinearJourneys.resize(orientedReadCount);
    //     uint64_t totalRemoved = 0;

    //     for (uint64_t i = 0; i < orientedReadCount; i++) {
    //         const OrientedReadId oid = OrientedReadId::fromValue(ReadId(i));
    //         const Shasta2Journey journey = (*shasta2Journeys)[oid];
    //         std::vector<Shasta2AnchorId>& linear = assembler.shasta2LinearJourneys[i];
    //         linear.clear();
    //         linear.reserve(journey.size());

    //         uint32_t prevEnd = 0; // end base position of the last kept anchor
    //         for (const Shasta2AnchorId anchorId : journey) {
    //             const uint32_t ordinal = shasta2Anchors->getOrdinal(anchorId, oid);
    //             const uint32_t pos = uint32_t(mkrs[i][ordinal].position);
    //             if (linear.empty() || pos >= prevEnd) {
    //                 linear.push_back(anchorId);
    //                 prevEnd = pos + uint32_t(k);
    //             } else {
    //                 ++totalRemoved;
    //             }
    //         }
    //     }
    //     cout << timestamp << "  Removed " << totalRemoved << " overlapping anchors." << endl;
    // }

    // Create the Shasta2AnchorGraph.
    const uint64_t minEdgeCoverage = 2;
    cout << timestamp << "Creating Shasta2AnchorGraph..." << endl;
    assembler.shasta2AnchorGraph = make_shared<Shasta2AnchorGraph>(
        *shasta2Anchors,
        *shasta2Journeys,
        minEdgeCoverage,
        threadCount);
    auto& shasta2AnchorGraph = assembler.shasta2AnchorGraph;

    // Save the pre-transitive-reduction Shasta2 anchor graph so the HTTP server
    // can load and visualize it even when we return before later assembly stages.
    shasta2AnchorGraph->saveAnchorGraph("Shasta2AnchorGraph");
    shasta2AnchorGraph->writeGfa("Shasta2AnchorGraph.gfa");
    shasta2AnchorGraph->writeBubbleFinderGraph("Shasta2AnchorGraph.graph");


    // Transitive reduction.
    // Shared parameters for local path-based simplification steps.
    const uint64_t transitiveReductionMaxEdgeCoverage = 10;
    const uint64_t transitiveReductionMaxDistance = 10;

    shasta2AnchorGraph->transitiveReduction(
        transitiveReductionMaxEdgeCoverage,
        transitiveReductionMaxDistance);

    shasta2AnchorGraph->writeGfa("Shasta2AnchorGraph-transitive-reduction.gfa");

    // Post-transitive-reduction cleanup:
    // cut low-read linear stalks that start at a tip and reach a branch point
    // before involving more than 3 distinct oriented reads across all anchors
    // in the traversed chain.
    const uint64_t maxWeakTipReadCount = 3;
    shasta2AnchorGraph->cutWeakStalksLeadingToBranch(
        *shasta2Anchors,
        maxWeakTipReadCount);

    // Save the final assembly-enabled state used by the HTTP server.
    shasta2AnchorGraph->saveAnchorGraph("Shasta2AnchorGraph");
    shasta2AnchorGraph->writeGfa("Shasta2AnchorGraph-transitive-reduction-weak-stalk-cut.gfa");
    // shasta2AnchorGraph->writeBubbleFinderGraph("Shasta2AnchorGraph-transitive-reduction-weak-stalk-cut.graph");

    // Next Shasta2 stage: create the AssemblyGraph, then simplify/assemble.
    cout << timestamp << "Creating Shasta2AssemblyGraph..." << endl;
    Shasta2AssemblyGraphOptions shasta2AssemblyGraphOptions;
    shasta2AssemblyGraphOptions.simplifyMaxIterationCount = 3;
    shasta2AssemblyGraphOptions.threadCount = threadCount;
    shasta2AssemblyGraphOptions.writeIntermediateAssemblyStages = true;
    assembler.shasta2AssemblyGraph = make_shared<Shasta2AssemblyGraph>(
        *shasta2Anchors,
        *shasta2Journeys,
        *shasta2AnchorGraph,
        shasta2AssemblyGraphOptions);
    auto& shasta2AssemblyGraph = assembler.shasta2AssemblyGraph;


    return;
    } // end of second assembly path block scope











    


    // // // Declare anchors pointer here to avoid scope issues
    // shared_ptr<mode3::Anchors> anchors;
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
    // assembler.mode3Assembly(threadCount, anchors, assemblerOptions.assemblyOptions.mode3Options, false);


    

    
    
    

    

    // cout << timestamp << "Simplifying and assembling Shasta2AssemblyGraph..." << endl;
    // shasta2AssemblyGraph->simplifyAndAssemble();

    // return;

    // dag.writeGfa("DirectedAnchorGraph-initial.gfa", true);

    // // Step 3: MBG-style cleaning phase
    // cout << timestamp << "Starting MBG-style cleaning phase..." << endl;
    // const uint64_t maxResolveLength = 500000;
    // const bool doRoundCleaning = true;
    // const bool doGuessworkCleaning = true;
    // const uint64_t maxUnconditionalResolveLength = 0;
    // const bool copycountFilterHeuristic = false;
    // const uint64_t maxLocalResolve = 0;
    // const bool resolvePalindromesGlobal = false;

    // // Step 3a: Remove low-coverage tips (MBG pre-resolve pass)
    // auto tipStats1 = dag.removeLowCoverageTips(3.0, 10.0, 10000);
    // if(tipStats1.nodesRemoved > 0) {
    //     cout << "  Removed " << tipStats1.nodesRemoved << " tip nodes, "
    //          << tipStats1.edgesRemoved << " edges." << endl;
    //     dag.unitigifyAll();
    //     cout << "  After tip removal: "
    //          << dag.nodeCount() << " nodes, "
    //          << dag.edgeCount() << " edges." << endl;
    // }

    // // Step 3b: Remove low-coverage crosslinks
    // auto crosslinkStats = dag.removeLowCoverageCrosslinks(2.0, 10);
    // if(crosslinkStats.edgesRemoved > 0) {
    //     cout << "  Removed " << crosslinkStats.edgesRemoved
    //          << " crosslink edges." << endl;
    //     dag.unitigifyAll();
    //     cout << "  After crosslink removal: "
    //          << dag.nodeCount() << " nodes, "
    //          << dag.edgeCount() << " edges." << endl;
    // }

    // // Step 3c: Copy-number cleaning (MBG pre-resolve, guesswork mode)
    // if(doGuessworkCleaning) {
    //     double totalCov = 0.0;
    //     uint64_t count = 0;
    //     for(uint64_t segId = 0; segId < dag.totalNodeCount(); ++segId) {
    //         if(!dag.nodeExists(segId)) continue;
    //         totalCov += dag.getPathCoverage(segId);
    //         count++;
    //     }
    //     const double avgCov = count > 0 ? totalCov / double(count) : 1.0;
    //     auto copyStats = dag.cleanComponentsByCopynumber(
    //         avgCov,
    //         50000,
    //         0,
    //         max(maxResolveLength, maxLocalResolve),
    //         {},
    //         0);
    //     if(copyStats.nodesRemoved > 0 || copyStats.edgesRemoved > 0) {
    //         cout << "  Copy-number cleaning removed "
    //              << copyStats.nodesRemoved << " nodes, "
    //              << copyStats.edgesRemoved << " edges." << endl;
    //         dag.unitigifyAll();
    //         cout << "  After copy-number cleaning: "
    //              << dag.nodeCount() << " nodes, "
    //              << dag.edgeCount() << " edges." << endl;
    //     }
    // }

    // cout << timestamp << "Cleaning phase complete." << endl;
    // dag.writeSummary(cout);
    // dag.writeGfa("DirectedAnchorGraph-After-Cleaning.gfa", true);

    

    // // Step 4: Resolution rounds (MBG two-pass: minCoverage, then 1)
    // const uint64_t initialMinEdgeSupport = 20;

    // // Step 4a: Resolve with minEdgeSupport = initial pass.
    // cout << timestamp
    //      << "Resolution step with minEdgeSupport="
    //      << initialMinEdgeSupport << endl;
    // dag.resolveRound(
    //     initialMinEdgeSupport,
    //     maxResolveLength,
    //     doRoundCleaning,
    //     doGuessworkCleaning,
    //     maxUnconditionalResolveLength,
    //     copycountFilterHeuristic,
    //     maxLocalResolve,
    //     resolvePalindromesGlobal);
    // dag.unitigifyAll();
    // cout << "  After resolution step " << initialMinEdgeSupport << ": "
    //      << dag.nodeCount() << " nodes, "
    //      << dag.edgeCount() << " edges, "
    //      << dag.pathCount() << " paths." << endl;

    // // Step 4b: MBG-style second pass with minimal support.
    // cout << timestamp << "Resolution final low-support pass (minEdgeSupport=1)" << endl;
    // dag.resolveRound(
    //     1,
    //     maxResolveLength,
    //     doRoundCleaning,
    //     doGuessworkCleaning,
    //     maxUnconditionalResolveLength,
    //     copycountFilterHeuristic,
    //     maxLocalResolve,
    //     resolvePalindromesGlobal);
    // dag.unitigifyAll();
    // cout << "  After final low-support pass: "
    //      << dag.nodeCount() << " nodes, "
    //      << dag.edgeCount() << " edges, "
    //      << dag.pathCount() << " paths." << endl;

    // // Step 5: Write final graph
    // dag.verifyEdgeConsistency();
    // dag.writeSummary(cout);
    // dag.writeGfa("DirectedAnchorGraph.gfa");
    // dag.writePaths("DirectedAnchorGraph.paths.gaf");

    // return;


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


    // // Create shasta2 anchors equivalent to the marker graph vertices.
    // // This allows downstream processing using shasta2 tools.
    // createShasta2Anchors(assembler, assemblerOptions, threadCount, anchors);

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
