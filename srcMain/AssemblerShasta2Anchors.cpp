#include "Reads.hpp"
#include "AssemblerShasta2Anchors.hpp"
#include "Assembler.hpp"
#include "mode3-Anchor.hpp"
#include "timestamp.hpp"
#include "shasta2/Anchor.hpp"
#include "shasta2/Reads.hpp"
#include "shasta2/Markers.hpp"
#include "shasta2/MarkerKmers.hpp"
#include "shasta2/MappedMemoryOwner.hpp"
#include "shasta2/MappedMemoryOwner.hpp"
#include "shasta2/Options.hpp"
#include "shasta2/Journeys.hpp"
#include "shasta2/AnchorGraph.hpp"
#include "shasta2/AssemblyGraph.hpp"
#include "shasta2/ReadSummary.hpp"
#include "shasta2/Tee.hpp"
#include "shasta2/performanceLog.hpp"

#include <thread>
#include <vector>
#include <memory>
#include <filesystem>
#include <iostream>
#include <algorithm>

using namespace std;
namespace fs = std::filesystem;

// Helper structs for safe memory aliasing (must be outside block scope).
namespace {
    // Layout of MemoryMapped::Vector<T>
    // Matches shasta2/MemoryMappedVector.hpp definition.
    template<class T> struct VectorLayout {
        void* header;
        T* data;
        bool isOpen;
        bool isOpenWithWriteAccess;
        string fileName;
    };

    // Layout of MemoryMapped::VectorOfVectors<T, Int>
    // Matches shasta2/MemoryMappedVectorOfVectors.hpp definition.
    template<class T, class Int> struct VecOfVecLayout {
        VectorLayout<Int> toc;
        VectorLayout<Int> count;
        VectorLayout<T> data;
        string name;
        size_t pageSize;
    };
    // Layout of LongBaseSequences
    // Matches shasta2/LongBaseSequence.hpp definition.
    struct LongBaseSequencesLayout {
        VectorLayout<uint64_t> baseCount;
        VecOfVecLayout<shasta2::Base, uint64_t> data;
        string name;
    };
} // parameter namespace

namespace shasta2 {
    // Define AssemblerInfo locally to avoid including "shasta2/Assembler.hpp"
    // which leads to header conflicts with Dinara's "Assembler.hpp".
    class AssemblerInfo {
    public:
        uint64_t k;
        double markerDensity;
        uint64_t largeDataPageSize;
    };
}

namespace dinara {

    // Helper function to perform the actual conversion (multithreaded).
    // This populates the shasta2::Anchors object from Dinara's MarkerGraph.
    static void populateAnchors(
        Assembler& assembler,
        shasta2::Anchors& anchors,
        uint64_t threadCount
    ) {
        if(threadCount == 0) threadCount = std::thread::hardware_concurrency();
        auto& vertices = assembler.markerGraph.vertices();
        uint64_t anchorCount = vertices.size();
        
        // Pass 1: Count the number of markers for each anchor.
        // We need to determine the storage required for each anchor in shasta2::Anchors.
        anchors.anchorMarkerInfos.beginPass1(anchorCount);
        
        std::vector<std::thread> threads;
        for(size_t t=0; t<threadCount; t++) {
            threads.emplace_back([&, t]() {
                size_t batchSize = anchorCount / threadCount;
                if(batchSize == 0) batchSize = 1;
                size_t begin = t * batchSize;
                size_t end = (t == threadCount - 1) ? anchorCount : begin + batchSize;
                if(begin >= anchorCount) return;
                
                // Reuse vector to avoid reallocation in inner loop.
                std::vector<dinara::OrientedReadId> uniqueReadIds;

                for(size_t i=begin; i<end; ++i) {
                    const auto& vertex = vertices[i];
                    uniqueReadIds.clear();
                    uniqueReadIds.reserve(vertex.size());

                    for(dinara::MarkerId markerId : vertex) {
                        dinara::OrientedReadId dReadId;
                        uint32_t ordinal;
                        std::tie(dReadId, ordinal) = assembler.findMarkerId(markerId);
                        uniqueReadIds.push_back(dReadId);
                    }

                    // Deduplicate read IDs. Shasta2 anchors enforce one marker per read per anchor.
                    std::sort(uniqueReadIds.begin(), uniqueReadIds.end());
                    auto last = std::unique(uniqueReadIds.begin(), uniqueReadIds.end());
                    size_t uniqueCount = std::distance(uniqueReadIds.begin(), last);

                    // Increment the counter for this anchor.
                    for(size_t k=0; k<uniqueCount; k++) {
                        anchors.anchorMarkerInfos.incrementCount(i);
                    }
                }
            });
        }
        for(auto& t : threads) t.join();
        threads.clear();
        
        // Pass 2: Populate the anchor data.
        anchors.anchorMarkerInfos.beginPass2();
        
        for(size_t t=0; t<threadCount; t++) {
            threads.emplace_back([&, t]() {
                size_t batchSize = anchorCount / threadCount;
                if(batchSize == 0) batchSize = 1;
                size_t begin = t * batchSize;
                size_t end = (t == threadCount - 1) ? anchorCount : begin + batchSize;
                if(begin >= anchorCount) return;
                
                // Reuse vector to avoid reallocation.
                std::vector<shasta2::AnchorMarkerInfo> amis;

                for(size_t i=begin; i<end; ++i) {
                    const auto& vertex = vertices[i];
                    amis.clear();
                    amis.reserve(vertex.size());

                    for(dinara::MarkerId markerId : vertex) {
                        dinara::OrientedReadId dReadId;
                        uint32_t ordinal;
                        std::tie(dReadId, ordinal) = assembler.findMarkerId(markerId);
                        
                        // Convert to shasta2 types.
                        shasta2::OrientedReadId sReadId = shasta2::OrientedReadId::fromValue(dReadId.getValue());
                        shasta2::AnchorMarkerInfo ami;
                        ami.orientedReadId = sReadId;
                        ami.positionInJourney = shasta2::invalid<uint32_t>;
                        amis.push_back(ami);
                    }

                    // Sort by OrientedReadId as required by Shasta2.
                    std::sort(amis.begin(), amis.end());

                    // Deduplicate logic: keep only unique orientedReadIds.
                    auto last = std::unique(amis.begin(), amis.end(), 
                        [](const shasta2::AnchorMarkerInfo& a, const shasta2::AnchorMarkerInfo& b){
                            return a.orientedReadId == b.orientedReadId;
                        });
                    amis.erase(last, amis.end());

                    // Reverse the vector before storing.
                    // Shasta's VectorOfVectors::store(i, value) typically fills data starting from the *end* 
                    // of the allocated block for vector i, moving backwards.
                    // To ensure the data ends up sorted [Small, Medium, Large] in memory, we must store 
                    // the elements in storage order, which is [Large, Medium, Small].
                    std::reverse(amis.begin(), amis.end());

                    // Store the data.
                    for(const auto& ami : amis) {
                        anchors.anchorMarkerInfos.store(i, ami);
                    }
                }
            });
        }
        for(auto& t : threads) t.join();
        threads.clear();
        
        anchors.anchorMarkerInfos.endPass2();
        
        // Initialize kmerIndex for AnchorInfos (used later in assembly).
        anchors.anchorInfos.resize(anchorCount);
        
        for(size_t t=0; t<threadCount; t++) {
            threads.emplace_back([&, t]() {
                size_t batchSize = anchorCount / threadCount;
                if(batchSize == 0) batchSize = 1;
                size_t begin = t * batchSize;
                size_t end = (t == threadCount - 1) ? anchorCount : begin + batchSize;
                if(begin >= anchorCount) return;

                for(size_t i=begin; i<end; ++i) {
                    anchors.anchorInfos[i].kmerIndex = shasta2::invalid<uint64_t>;
                }
            });
        }
        for(auto& t : threads) t.join();
    }

    // Populate shasta2::Anchors from an existing dinara::mode3::Anchors.
    // This is used to keep Shasta2 downstream assembly consistent with the same anchor
    // selection logic used by Dinara (for example "BestPerOverlapInterval").
    static void populateAnchorsFromMode3Anchors(
        const mode3::Anchors& dinaraAnchors,
        shasta2::Anchors& anchors,
        uint64_t threadCount)
    {
        if(threadCount == 0) threadCount = std::thread::hardware_concurrency();

        const uint64_t anchorCount = dinaraAnchors.size();

        anchors.anchorMarkerInfos.beginPass1(anchorCount);

        std::vector<std::thread> threads;
        threads.reserve(threadCount);
        for(uint64_t t=0; t<threadCount; t++) {
            threads.emplace_back([&, t]() {
                uint64_t batchSize = anchorCount / threadCount;
                if(batchSize == 0) batchSize = 1;
                const uint64_t begin = t * batchSize;
                const uint64_t end = (t == threadCount - 1) ? anchorCount : std::min(anchorCount, begin + batchSize);
                if(begin >= anchorCount) return;

                std::vector<shasta2::AnchorMarkerInfo> amis;
                for(uint64_t anchorId=begin; anchorId<end; anchorId++) {
                    const mode3::Anchor a = dinaraAnchors[anchorId];
                    amis.clear();
                    amis.reserve(a.size());
                    for(const auto& mi : a) {
                        shasta2::AnchorMarkerInfo ami;
                        ami.orientedReadId = shasta2::OrientedReadId::fromValue(mi.orientedReadId.getValue());
                        ami.positionInJourney = shasta2::invalid<uint32_t>;
                        amis.push_back(ami);
                    }

                    std::sort(amis.begin(), amis.end());
                    const auto last = std::unique(amis.begin(), amis.end(),
                        [](const shasta2::AnchorMarkerInfo& x, const shasta2::AnchorMarkerInfo& y) {
                            return x.orientedReadId == y.orientedReadId;
                        });
                    const uint64_t uniqueCount = uint64_t(std::distance(amis.begin(), last));
                    for(uint64_t k=0; k<uniqueCount; k++) {
                        anchors.anchorMarkerInfos.incrementCount(anchorId);
                    }
                }
            });
        }
        for(auto& th : threads) th.join();
        threads.clear();

        anchors.anchorMarkerInfos.beginPass2();

        for(uint64_t t=0; t<threadCount; t++) {
            threads.emplace_back([&, t]() {
                uint64_t batchSize = anchorCount / threadCount;
                if(batchSize == 0) batchSize = 1;
                const uint64_t begin = t * batchSize;
                const uint64_t end = (t == threadCount - 1) ? anchorCount : std::min(anchorCount, begin + batchSize);
                if(begin >= anchorCount) return;

                std::vector<shasta2::AnchorMarkerInfo> amis;
                for(uint64_t anchorId=begin; anchorId<end; anchorId++) {
                    const mode3::Anchor a = dinaraAnchors[anchorId];
                    amis.clear();
                    amis.reserve(a.size());
                    for(const auto& mi : a) {
                        shasta2::AnchorMarkerInfo ami;
                        ami.orientedReadId = shasta2::OrientedReadId::fromValue(mi.orientedReadId.getValue());
                        ami.positionInJourney = shasta2::invalid<uint32_t>;
                        amis.push_back(ami);
                    }

                    std::sort(amis.begin(), amis.end());
                    const auto last = std::unique(amis.begin(), amis.end(),
                        [](const shasta2::AnchorMarkerInfo& x, const shasta2::AnchorMarkerInfo& y) {
                            return x.orientedReadId == y.orientedReadId;
                        });
                    amis.erase(last, amis.end());
                    std::reverse(amis.begin(), amis.end());

                    for(const auto& ami : amis) {
                        anchors.anchorMarkerInfos.store(anchorId, ami);
                    }
                }
            });
        }
        for(auto& th : threads) th.join();
        threads.clear();

        anchors.anchorMarkerInfos.endPass2();

        anchors.anchorInfos.resize(anchorCount);
        for(uint64_t t=0; t<threadCount; t++) {
            threads.emplace_back([&, t]() {
                uint64_t batchSize = anchorCount / threadCount;
                if(batchSize == 0) batchSize = 1;
                const uint64_t begin = t * batchSize;
                const uint64_t end = (t == threadCount - 1) ? anchorCount : std::min(anchorCount, begin + batchSize);
                if(begin >= anchorCount) return;

                for(uint64_t i=begin; i<end; ++i) {
                    anchors.anchorInfos[i].kmerIndex = shasta2::invalid<uint64_t>;
                }
            });
        }
        for(auto& th : threads) th.join();
    }

    [[maybe_unused]] static void populateAnchorsFromMode3AnchorsWithShasta2KmerFilters(
        const mode3::Anchors& dinaraAnchors,
        shasta2::Anchors& anchors,
        const shasta2::Markers& markers,
        const shasta2::MarkerKmers& markerKmers,
        uint64_t k,
        const vector<uint64_t>& maxAnchorRepeatLength,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint64_t /*threadCount*/)
    {
        const uint64_t inputAnchorCount = dinaraAnchors.size();
        uint64_t skippedCoverage = 0;
        uint64_t skippedEmpty = 0;
        uint64_t skippedLowComplexity = 0;
        uint64_t skippedMissingKmer = 0;
        uint64_t skippedDupReadIdInKmer = 0;

        // Cache for "this canonical marker k-mer has duplicate ReadIds in MarkerKmers".
        // -1: unknown, 0: ok, 1: has duplicate ReadIds (repeat within a read -> skip).
        vector<int8_t> hasDuplicateReadIdByKmerIndex(markerKmers.size(), int8_t(-1));

        auto shouldSkipKmerDueToReadDuplicates = [&](uint64_t kmerIndex) -> bool {
            int8_t& cached = hasDuplicateReadIdByKmerIndex[kmerIndex];
            if(cached != int8_t(-1)) {
                return cached != 0;
            }
            const auto markerInfos = markerKmers[kmerIndex];
            bool hasDuplicate = false;
            for(size_t i=1; i<markerInfos.size(); i++) {
                if(markerInfos[i].orientedReadId.getReadId() == markerInfos[i-1].orientedReadId.getReadId()) {
                    hasDuplicate = true;
                    break;
                }
            }
            cached = hasDuplicate ? int8_t(1) : int8_t(0);
            return hasDuplicate;
        };

        auto shouldSkipKmerDueToLowComplexity = [&](const shasta2::Kmer& kmer) -> bool {
            for(uint64_t i=0; i<maxAnchorRepeatLength.size(); i++) {
                const uint64_t period = i + 1;
                const uint64_t maxAllowedCopyNumber = maxAnchorRepeatLength[i];
                if(kmer.countExactRepeatCopies(period, k) > maxAllowedCopyNumber) {
                    return true;
                }
            }
            return false;
        };

        vector<uint64_t> keptInputAnchorIds;
        keptInputAnchorIds.reserve(inputAnchorCount / 4);

        for(uint64_t inputAnchorId=0; inputAnchorId<inputAnchorCount; inputAnchorId++) {
            const mode3::Anchor a = dinaraAnchors[inputAnchorId];
            const uint64_t coverage = a.size();
            if(coverage < minAnchorCoverage || coverage > maxAnchorCoverage) {
                ++skippedCoverage;
                continue;
            }
            if(a.empty()) {
                ++skippedEmpty;
                continue;
            }

            const auto& first = a.front();
            const shasta2::OrientedReadId orientedReadId =
                shasta2::OrientedReadId::fromValue(first.orientedReadId.getValue());
            const shasta2::Kmer kmer0 = markers.getKmer(orientedReadId, first.ordinal0);

            if(shouldSkipKmerDueToLowComplexity(kmer0)) {
                ++skippedLowComplexity;
                continue;
            }

            const shasta2::Kmer kmerRc = kmer0.reverseComplement(k);
            const shasta2::Kmer canonicalKmer = (kmer0 <= kmerRc) ? kmer0 : kmerRc;
            const uint64_t kmerIndex = markerKmers.getGlobalIndex(canonicalKmer);
            if(kmerIndex == shasta2::invalid<uint64_t>) {
                ++skippedMissingKmer;
                continue;
            }
            if(shouldSkipKmerDueToReadDuplicates(kmerIndex)) {
                ++skippedDupReadIdInKmer;
                continue;
            }

            keptInputAnchorIds.push_back(inputAnchorId);
        }

        const uint64_t anchorCount = keptInputAnchorIds.size();
        cout << timestamp << "Shasta2 k-mer filters kept " << anchorCount << " / " << inputAnchorCount
             << " overlap-derived anchors."
             << " skipped: coverage=" << skippedCoverage
             << " empty=" << skippedEmpty
             << " lowComplexity=" << skippedLowComplexity
             << " missingKmer=" << skippedMissingKmer
             << " dupReadIdInKmer=" << skippedDupReadIdInKmer
             << endl;
        anchors.anchorMarkerInfos.beginPass1(anchorCount);

        // Pass 1: count marker infos for each kept anchor (distinct ReadIds).
        for(uint64_t newAnchorId=0; newAnchorId<anchorCount; newAnchorId++) {
            const mode3::Anchor a = dinaraAnchors[keptInputAnchorIds[newAnchorId]];
            ReadId prevReadId = invalid<ReadId>;
            uint64_t uniqueReadCount = 0;
            for(const auto& mi : a) {
                const ReadId readId = mi.orientedReadId.getReadId();
                if(readId != prevReadId) {
                    ++uniqueReadCount;
                    prevReadId = readId;
                }
            }
            for(uint64_t i=0; i<uniqueReadCount; i++) {
                anchors.anchorMarkerInfos.incrementCount(newAnchorId);
            }
        }

        anchors.anchorMarkerInfos.beginPass2();
        std::vector<shasta2::AnchorMarkerInfo> amis;
        for(uint64_t newAnchorId=0; newAnchorId<anchorCount; newAnchorId++) {
            const mode3::Anchor a = dinaraAnchors[keptInputAnchorIds[newAnchorId]];
            amis.clear();
            amis.reserve(a.size());
            for(const auto& mi : a) {
                shasta2::AnchorMarkerInfo ami;
                ami.orientedReadId = shasta2::OrientedReadId::fromValue(mi.orientedReadId.getValue());
                ami.positionInJourney = shasta2::invalid<uint32_t>;
                amis.push_back(ami);
            }

            std::sort(amis.begin(), amis.end());
            const auto last = std::unique(amis.begin(), amis.end(),
                [](const shasta2::AnchorMarkerInfo& x, const shasta2::AnchorMarkerInfo& y) {
                    return x.orientedReadId == y.orientedReadId;
                });
            amis.erase(last, amis.end());
            std::reverse(amis.begin(), amis.end());
            for(const auto& ami : amis) {
                anchors.anchorMarkerInfos.store(newAnchorId, ami);
            }
        }
        anchors.anchorMarkerInfos.endPass2();

        anchors.anchorInfos.resize(anchorCount);
        for(uint64_t i=0; i<anchorCount; i++) {
            anchors.anchorInfos[i].kmerIndex = shasta2::invalid<uint64_t>;
        }
    }

    shasta2::Tee shastaTee;
    ofstream shastaLog;

    void createShasta2Anchors(
        dinara::Assembler& assembler,
        const dinara::AssemblerOptions& dinaraOptions,
        uint64_t threadCount
    ) {
        createShasta2Anchors(assembler, dinaraOptions, threadCount, {});
    }

    void createShasta2Anchors(
        dinara::Assembler& assembler,
        const dinara::AssemblerOptions& dinaraOptions,
        uint64_t threadCount,
        const std::shared_ptr<const dinara::mode3::Anchors>& precomputedDinaraAnchors
    ) {
        cout << timestamp << "Creating Shasta2 Anchors..." << endl;

        const string& shastaOut = dinaraOptions.commandLineOnlyOptions.shasta2OutputDirectory;
        
        // 1. Ensure Dinara's data is accessible via absolute paths before we chdir.
        string originalDinaraPrefix = assembler.largeDataFileNamePrefix;
        if (not originalDinaraPrefix.empty() and not fs::path(originalDinaraPrefix).is_absolute()) {
             assembler.largeDataFileNamePrefix = fs::absolute(originalDinaraPrefix).string();
        }

        // 2. Setup Shasta2 Output Directory.
        if (not fs::exists(shastaOut)) {
            fs::create_directory(shastaOut);
        }
        
        // 3. chdir avoided. We use explicit paths.
        // fs::path originalCwd = fs::current_path();
        // fs::current_path(shastaOut);

        // 4. Initialize Shasta2 Options.
        // We use the same argc/argv trick to get default values, then override.
        {
        int argc = 1;
        char* name = (char*)"dinara";
        char* argv[] = {name};
        shasta2::Options options(argc, argv);
        
        options.assemblyDirectory = fs::absolute(shastaOut).string();
        options.threadCount = threadCount;
        options.memoryMode = dinaraOptions.commandLineOnlyOptions.memoryMode;
        options.memoryBacking = dinaraOptions.commandLineOnlyOptions.memoryBacking;
        options.writeIntermediateAssemblyStages = true;

        // Write shasta2.conf
        {
            const string fileName = (fs::path(shastaOut) / "shasta2.conf").string();
            ofstream configurationFile(fileName);
            options.write(configurationFile);
        }

        // Open logs
        shasta2::openPerformanceLog( (fs::path(shastaOut) / "performance.log").string() );
        shasta2::performanceLog << timestamp << "Shasta2 Assembly begins." << endl;
        shastaLog.open( (fs::path(shastaOut) / "stdout.log").string() );
        shastaTee.duplicate(cout, shastaLog);

        // Determine dataDirectory prefix for shastaOwner.
        string dataDirectory = "";
        size_t shastaPageSize = assembler.largeDataPageSize;
        if (options.memoryMode == "filesystem") {
            dataDirectory = (fs::path(shastaOut) / "Data/").string();
            if (not fs::exists(dataDirectory)) {
                fs::create_directory(dataDirectory);
            }
        }

        shasta2::MappedMemoryOwner shastaOwner(dataDirectory, shastaPageSize);

        // Options: Anchor creation.
        // Use the same coverage range used for Dinara mode3 anchor creation, so Shasta2 sees
        // the same anchors when we populate anchors from Dinara anchors.
        uint64_t minAnchorCoverageDinara = dinaraOptions.assemblyOptions.mode3Options.minAnchorCoverage;
        uint64_t maxAnchorCoverageDinara = dinaraOptions.assemblyOptions.mode3Options.maxAnchorCoverage;
        if((minAnchorCoverageDinara == 0) and (maxAnchorCoverageDinara == 0)) {
            tie(minAnchorCoverageDinara, maxAnchorCoverageDinara) = assembler.getPrimaryCoverageRange();
            minAnchorCoverageDinara = uint64_t(std::round(
                double(minAnchorCoverageDinara) * dinaraOptions.assemblyOptions.mode3Options.minAnchorCoverageMultiplier));
            maxAnchorCoverageDinara = uint64_t(std::round(
                double(maxAnchorCoverageDinara) * dinaraOptions.assemblyOptions.mode3Options.maxAnchorCoverageMultiplier));
        }

        options.minAnchorCoverage = uint32_t(minAnchorCoverageDinara);
        options.maxAnchorCoverage = uint32_t(maxAnchorCoverageDinara);
        options.maxAnchorRepeatLength = {8, 3, 3, 3, 3};

        // Options: Anchor graph.
        options.minAnchorGraphEdgeCoverage = 4;
        options.transitiveReductionMaxEdgeCoverage = 10;
        options.transitiveReductionMaxDistance = 20;

        // Options: Detangling.
        options.detangleMaxLogP = 50.;
        options.detangleMinLogPDelta = 30.;
        options.detangleEpsilon = 0.1;

        // Options: Read following.
        options.readFollowingMinCommonCount = 2;
        options.readFollowingMinCorrectedJaccard = 0.7;
        options.readFollowingSegmentLengthThreshold = 30000;
        options.readFollowingPruneLength = 10000;

        // =================================================================================
        // DANGER: DIRECT MEMORY ALIASING IMPLEMENTATION
        // =================================================================================
        // In anonymous memory mode, we cannot open files to initialize shasta2 objects 
        // because the data exists only in RAM. 
        //
        // Solution: We manually construct the shasta2 objects and "alias" their internal 
        // pointers to point to the existing data owned by the 'dinara::assembler' instance.
        //
        // SAFETY CRITICAL REQUIREMENTS:
        // 1. LIFETIME: The shasta2 objects must NOT free/unmap the memory when they are destroyed.
        //    The memory is owned by 'assembler' and must remain valid.
        //    We achieve this by zeroing out the shasta2 objects before destruction (via custom deleter).
        // 2. LAYOUT: We assume 'dinara' and 'shasta2' classes use binary-compatible 
        //    implementations for their member types (LongBaseSequences, VectorOfVectors, etc.).
        //    We perform member-wise copying to account for differences in class member order/presence.
        // =================================================================================

        // =================================================================================
        // DANGER: DIRECT MEMORY ALIASING IMPLEMENTATION
        // =================================================================================
        // ... (comments kept as is, but adding Info creation here)

        // 0. Create Shasta2 AssemblerInfo (Required for shasta2 executable to start).
        // This must be done manually since we are not using shasta2::Assembler.
        {
            shasta2::MemoryMapped::Object<shasta2::AssemblerInfo> info;
            info.createNew(shastaOwner.largeDataName("Info"), shastaPageSize);
            info->k = assembler.assemblerInfo->k;
            info->markerDensity = 0.1; // Default/Symbolic value (not strictly used by Dinara)
            info->largeDataPageSize = shastaPageSize;
        }

        // ---------------------------------------------------------------------------------
        // 1. Setup Shasta2 Reads
        // ---------------------------------------------------------------------------------
        
        // Custom deleter detects anonymous mode and zeroes out the object to prevent 'munmap' calls.
        shared_ptr<shasta2::Reads> shastaReadsPtr(new shasta2::Reads(), [&](shasta2::Reads* p) {
            if (p && assembler.largeDataFileNamePrefix.empty()) {
                struct ReadsLayout {
                    LongBaseSequencesLayout reads;
                    VecOfVecLayout<char, uint64_t> readNames;
                    VectorLayout<shasta2::ReadId> readIdsSortedByName;
                };
                auto* n = reinterpret_cast<ReadsLayout*>(p);
                
                // LongBaseSequences (reads)
                n->reads.baseCount.isOpen = false;
                n->reads.data.toc.isOpen = false;
                n->reads.data.count.isOpen = false;
                n->reads.data.data.isOpen = false;
                new (&n->reads.baseCount.fileName) string();
                new (&n->reads.data.toc.fileName) string();
                new (&n->reads.data.data.fileName) string();
                new (&n->reads.data.name) string();
                new (&n->reads.name) string();

                // VectorOfVectors (readNames)
                n->readNames.toc.isOpen = false;
                n->readNames.count.isOpen = false;
                n->readNames.data.isOpen = false;
                new (&n->readNames.toc.fileName) string();
                new (&n->readNames.data.fileName) string();
                new (&n->readNames.name) string();

                // Vector (readIdsSortedByName)
                n->readIdsSortedByName.isOpen = false;
                new (&n->readIdsSortedByName.fileName) string();
            }
            delete p;
        });

        if (assembler.largeDataFileNamePrefix.empty()) {
            cout << timestamp << "Anonymous memory detected. Aliasing Reads memory directly..." << endl;
             
            // Define the memory layout of shasta2::Reads (and Dinara::Reads prefix) to allow writing to its private members.
            struct Shasta2ReadsLayout {
                LongBaseSequencesLayout reads;
                VecOfVecLayout<char, uint64_t> readNames;
                VectorLayout<shasta2::ReadId> readIdsSortedByName;
            };
            
            // Cast the opaque pointer to our layout struct.
            auto* target = reinterpret_cast<Shasta2ReadsLayout*>(shastaReadsPtr.get());
            
            // Cast source (Dinara Reads) to the same layout (prefix matches).
            const auto* source = reinterpret_cast<const Shasta2ReadsLayout*>(&assembler.getReads());

            // Member 1: 'reads' (LongBaseSequences)
            std::memcpy((void*)&target->reads, (const void*)&source->reads, sizeof(target->reads));

            // Member 2: 'readNames' (VectorOfVectors<char>)
            std::memcpy((void*)&target->readNames, (const void*)&source->readNames, sizeof(target->readNames));

            // Member 3: 'readIdsSortedByName' (Vector<ReadId>)
            std::memcpy((void*)&target->readIdsSortedByName, (const void*)&source->readIdsSortedByName, sizeof(target->readIdsSortedByName));
            
            // Note: We SKIP 'readRepeatCounts', 'readMetaData', etc. which exist in dinara::Reads 
            // but not in shasta2::Reads or are not required here.

        } else {
            // Standard File-Based Access: Standalone Compatibility via Symlinks
            try {
                // Determine the Dinara Data directory
                fs::path dinaraDataPath = fs::path(assembler.largeDataFileNamePrefix);
                if (dinaraDataPath.has_filename()) {
                    dinaraDataPath = dinaraDataPath.parent_path();
                }

                // Selective Symlinking: Only link what Shasta2 actually needs.
                if (fs::exists(dinaraDataPath)) {
                    // List of files that Shasta2 needs to function using Dinara's data.
                    // List of files that Shasta2 needs to function using Dinara's data.
                    // List of files that Shasta2 needs to function using Dinara's data.
                    // We define {SourceCandidate, Target}. We list preferred/likely sources first.
                    // The loop below will try candidates in order and link the first one that exists.
                    const std::vector<std::pair<string, string>> requiredFiles = {
                        // Reads-BaseCount (Vector -> Single File typically, but check .toc too)
                        {"Reads-BaseCount", "Reads-BaseCount"},
                        {"Reads-BaseCount.toc", "Reads-BaseCount"}, // Fallback

                        // Reads sequence data: Try "Reads-Bases" (likely) then "Reads-Data" (legacy/mapped)
                        {"Reads-Bases.toc", "Reads-Bases.toc"},
                        {"Reads-Data.toc", "Reads-Bases.toc"},

                        {"Reads-Bases.count", "Reads-Bases.count"},
                        {"Reads-Data.count", "Reads-Bases.count"},

                        {"Reads-Bases.data", "Reads-Bases.data"},
                        {"Reads-Data.data", "Reads-Bases.data"},
                        
                        // ReadNames (VectorOfVectors)
                        {"ReadNames.toc", "ReadNames.toc"},
                        {"ReadNames.count", "ReadNames.count"},
                        {"ReadNames.data", "ReadNames.data"},

                        // Markers (VectorOfVectors)
                        {"Markers.toc", "Markers.toc"},
                        {"Markers.count", "Markers.count"},
                        {"Markers.data", "Markers.data"},

                        // ReadIdsSortedByName (Vector -> Single File typically)
                        {"ReadIdsSortedByName", "ReadIdsSortedByName"},
                        {"ReadIdsSortedByName.toc", "ReadIdsSortedByName"}
                    };

                    for (const auto& filePair : requiredFiles) {
                        fs::path sourcePath = dinaraDataPath / filePair.first;
                        fs::path targetPath = fs::path(dataDirectory) / filePair.second;
                        
                        if (fs::exists(sourcePath) && !fs::exists(targetPath)) {
                            fs::create_symlink(sourcePath, targetPath);
                        }
                    }
                }

                shastaReadsPtr->access(
                    dataDirectory + "Reads",
                    dataDirectory + "ReadNames",
                    dataDirectory + "ReadIdsSortedByName"
                );
            } catch (const std::exception& e) {
                 cout << timestamp << "Error accessing Shasta2 Reads: " << e.what() << endl;
                 throw;
            }
        }


        // ---------------------------------------------------------------------------------
        // 2. Setup Shasta2 Markers
        // ---------------------------------------------------------------------------------
        
        shasta2::Markers* shastaMarkersRawPtr = nullptr;
        
        if (assembler.largeDataFileNamePrefix.empty()) {
            cout << timestamp << "Anonymous memory detected. Aliasing Markers memory directly..." << endl;
             
            // Shasta2 Markers constructor (the one taking 'reads') UNCONDITIONALLY attempts to open 
            // the "Markers.toc" and "Markers.data" files via 'accessExistingReadOnly'.
            // In anonymous mode, MappedMemoryOwner returns "" for names, causing access to ".toc".
            // Even if we gave it a name, it would look for a file that doesn't exist.
            //
            // Workaround: We use a temporary owner with a prefix to create valid dummy files 
            // to satisfy the constructor's checks, then hijack the object.
            
            string tempPrefix = "Temp_Shasta2_Fix_";
            shasta2::MappedMemoryOwner tempOwner(tempPrefix, assembler.largeDataPageSize);
            
            // Create dummy files.
            {
                shasta2::MemoryMapped::Vector<uint64_t> dummyToc;
                dummyToc.createNew(tempOwner.largeDataName("Markers.toc"), assembler.largeDataPageSize);
                shasta2::MemoryMapped::Vector<shasta2::Marker> dummyData;
                dummyData.createNew(tempOwner.largeDataName("Markers.data"), assembler.largeDataPageSize);
            }
            
            try {
                // Now construct the object using the dummy files.
                // This initializes the 'reads' and 'k' references properly.
                shastaMarkersRawPtr = new shasta2::Markers(
                    tempOwner,
                    assembler.assemblerInfo->k,
                    shastaReadsPtr);
            } catch(...) {
                // Cleanup if construction fails.
                fs::remove(tempOwner.largeDataName("Markers.toc"));
                fs::remove(tempOwner.largeDataName("Markers.data"));
                throw;
            }
            
            // Cleanup dummy files from disk. 
            // The object currently holds open file descriptors to them.
            fs::remove(tempOwner.largeDataName("Markers.toc"));
            fs::remove(tempOwner.largeDataName("Markers.data"));


            // We need to alias the 'VectorOfVectors' base class part of Markers.
            using TargetVectorType = shasta2::MemoryMapped::VectorOfVectors<shasta2::Marker, uint64_t>;
            
            // Access the source vector from the assembler.
            auto* sourceVec = assembler.markers.get();

            if (!sourceVec) {
                throw runtime_error("Assembler markers not initialized.");
            }
            
            // Access the base class pointer on the target.
            auto* targetVec = static_cast<TargetVectorType*>(shastaMarkersRawPtr);

            // Close the dummy file handles to prevent leaks and prepare for overwrite.
            // VectorOfVectors::close() syncs and unmaps.
            targetVec->close();

            // Overwrite the VectorOfVectors state in 'target' with 'source' (Aliasing).
            // This copies the internal pointers (header, data) so both objects point to the same Dinara memory.
            if (sizeof(TargetVectorType) != sizeof(shasta2::MemoryMapped::VectorOfVectors<dinara::CompressedMarker, uint64_t>)) {
                 if (sizeof(dinara::CompressedMarker) != sizeof(shasta2::Marker)) {
                     throw runtime_error("Layout Mismatch: dinara::CompressedMarker size != shasta2::Marker size.");
                 }
                 throw runtime_error("Layout Mismatch: VectorOfVectors container size mismatch.");
            }
            std::memcpy(
                static_cast<void*>(targetVec),
                static_cast<const void*>(sourceVec),
                sizeof(TargetVectorType));
            
        } else {
            // Standard File-Based Access: Selected symlinks are already created above.
            // We just need to make sure the Markers object can access them.
            try {
                 shastaMarkersRawPtr = new shasta2::Markers(
                     shastaOwner,
                     assembler.assemblerInfo->k,
                     shastaReadsPtr);
            } catch(const std::exception& e) {
                 cout << timestamp << "Error accessing Shasta2 Markers: " << e.what() << endl;
                 throw;
            }
        }
        
        // Wrap the raw pointer in a shared_ptr with the safety deleter.
        shared_ptr<shasta2::Markers> shastaMarkersPtr(shastaMarkersRawPtr, [&](shasta2::Markers* p) {
            if (p) {
                if (assembler.largeDataFileNamePrefix.empty()) {
                    // Safety: Defuse destruction to prevent unmapping aliased memory.
                    // Access the VectorOfVectors private members using a layout struct.
                    
                    using TargetVectorType = shasta2::MemoryMapped::VectorOfVectors<shasta2::Marker, uint64_t>;
                    auto* vecBase = static_cast<TargetVectorType*>(p);
                    
                    // Cast to layout struct to access private members.
                    auto* layout = reinterpret_cast<VecOfVecLayout<shasta2::Marker, uint64_t>*>(vecBase);
                    
                    // Neutralize the vectors.
                    layout->toc.isOpen = false;
                    layout->count.isOpen = false;
                    layout->data.isOpen = false;
                    
                    // IMPORTANT: Reset strings to empty to avoid double-free in destructor.
                    // Since we memcpy'd from Dinara, these strings might point to the same buffers.
                    // Placement new (empty string) is safer than just assignment if the memory is already aliased.
                    new (&layout->toc.fileName) string();
                    new (&layout->count.fileName) string();
                    new (&layout->data.fileName) string();
                    new (&layout->name) string();
                }
                
                // For the manually allocated case, we used placement new (or operator new + constructor).
                // Standard 'delete' should be fine as it calls destructor then frees memory.
                delete p; 
            }
        });

        // Reference for convenience in existing code.
        shasta2::Markers& shastaMarkers = *shastaMarkersPtr;


        // ---------------------------------------------------------------------------------
        // 3. Compute MarkerKmers
        // ---------------------------------------------------------------------------------
        // We need dummy ReadSummaries because Dinara doesn't persist them, but Shasta2 API requires them.
        shasta2::MemoryMapped::Vector<shasta2::ReadSummary> readSummaries;
        readSummaries.createNew(
            shastaOwner.largeDataName("ReadSummaries"), 
            assembler.largeDataPageSize
        );
        readSummaries.resize(shastaReadsPtr->readCount()); 
        
        cout << timestamp << "Computing Shasta2 MarkerKmers..." << endl;
        shasta2::MarkerKmers shastaMarkerKmers(
             assembler.assemblerInfo->k,
             shastaOwner,      
             *shastaReadsPtr,  
             readSummaries,    
             shastaMarkers,    
             threadCount       
        );

        // 4. Create and Populate Shasta2 Anchors.
        cout << timestamp << "Creating and Populating Shasta2 Anchors..." << endl;
        auto shastaAnchors = make_shared<shasta2::Anchors>(
            "Anchors",
            shastaOwner,
            *shastaReadsPtr,
            assembler.assemblerInfo->k,
            shastaMarkerKmers,
            minAnchorCoverageDinara,
            maxAnchorCoverageDinara,
            options.maxAnchorRepeatLength,
            vector<uint64_t>{4, 12, 24},  // minAnchorDistinctSubkmerCount (shasta2 default)
            threadCount);

        if(precomputedDinaraAnchors) {
            cout << timestamp
                 << "Populating Shasta2 anchors from precomputed Dinara anchors (exact handoff), count="
                 << precomputedDinaraAnchors->size() << endl;
            populateAnchorsFromMode3Anchors(*precomputedDinaraAnchors, *shastaAnchors, threadCount);
        } else {
            const string& anchorMethod = dinaraOptions.assemblyOptions.mode3Options.anchorCreationMethod;
            if(anchorMethod == "FromMarkerGraphVerticesBestPerOverlapInterval") {
                cout << timestamp << "Populating Shasta2 anchors using Dinara BestPerOverlapInterval anchors..." << endl;
                const auto dinaraAnchors = assembler.createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval(
                    minAnchorCoverageDinara,
                    maxAnchorCoverageDinara,
                    threadCount);
                populateAnchorsFromMode3Anchors(*dinaraAnchors, *shastaAnchors, threadCount);
            } else if(anchorMethod == "FromMarkerGraphVerticesAtOverlapEvents") {
                cout << timestamp << "Populating Shasta2 anchors using Dinara OverlapEvents anchors..." << endl;
                const auto dinaraAnchors = assembler.createAnchorsFromMarkerGraphVerticesAtOverlapEvents(
                    minAnchorCoverageDinara,
                    maxAnchorCoverageDinara,
                    threadCount);
                populateAnchorsFromMode3Anchors(*dinaraAnchors, *shastaAnchors, threadCount);
            } else if(anchorMethod == "FromOverlapsBestPerOverlapInterval") {
                cout << timestamp << "Populating Shasta2 anchors using Dinara overlap-only BestPerOverlapInterval anchors..." << endl;
                const auto dinaraAnchors = assembler.createAnchorsFromOverlapsBestPerOverlapInterval(
                    minAnchorCoverageDinara,
                    maxAnchorCoverageDinara,
                    threadCount);
                // Use exactly Dinara's overlap-derived anchors as-is.
                // Dinara already applies overlap-based filtering/selection; additional Shasta2 k-mer filters
                // would change the anchor set and make downstream comparisons confusing.
                populateAnchorsFromMode3Anchors(*dinaraAnchors, *shastaAnchors, threadCount);
            } else {
                // Default: treat marker graph vertices as anchors.
                if((not assembler.markerGraph.verticesPointer) or
                    (not assembler.markerGraph.verticesPointer->isOpen())) {
                    throw runtime_error(
                        "createShasta2Anchors: anchorCreationMethod=" + anchorMethod +
                        " requires marker graph vertices, but markerGraph is not open/initialized. "
                        "Use --Assembly.mode3.anchorCreationMethod FromMarkerGraphVertices... or "
                        "FromOverlapsBestPerOverlapInterval.");
                }
                populateAnchors(assembler, *shastaAnchors, threadCount);
            }
        }
	        
	        cout << timestamp << "Shasta2 Anchors created (" << shastaAnchors->size() << " anchors)." << endl;

	        
	        // 5. Downstream Assembly Pipeline (Journeys -> AnchorGraph -> AssemblyGraph).
        cout << timestamp << "Proceeding with Downstream Shasta2 Assembly..." << endl;
        
        // Create Journeys.
        cout << timestamp << "Creating Journeys..." << endl;
        shasta2::Journeys journeys(
            2 * shastaReadsPtr->readCount(),
            shastaAnchors,
            threadCount,
            shastaOwner);
        

        // Store Anchor Gaps (Ported from Shasta2 Assembler.cpp).
        cout << timestamp << "Storing Anchor Gaps..." << endl;
        {
            const uint32_t kHalf = uint32_t(shastaAnchors->kHalf);
            // Loop over all ReadIds.
            for(shasta2::ReadId readId=0; readId<shastaReadsPtr->readCount(); readId++) {
                 shasta2::ReadSummary& readSummary = readSummaries[readId];
                 const uint32_t readLength = uint32_t(shastaReadsPtr->getReadSequenceLength(readId));
                 
                 // Put it on strand 0.
                 const shasta2::OrientedReadId orientedReadId(readId, 0);

                 // Get journey.
                 const auto journey = journeys[orientedReadId];

                 if(journey.empty()) {
                     readSummary.initialAnchorGap = readLength;
                     readSummary.middleAnchorGap = readLength;
                     readSummary.finalAnchorGap = readLength;
                     continue;
                 }

                 // Compute largest gap between adjacent anchors.
                 uint32_t maxGap = 0;
                 for(uint64_t i1=1; i1<journey.size(); i1++) {
                     const uint64_t i0 = i1 - 1;
                     const shasta2::AnchorId anchorId0 = journey[i0];
                     const shasta2::AnchorId anchorId1 = journey[i1];
                     
                     const uint32_t position0 = shastaAnchors->getPosition(anchorId0, orientedReadId);
                     const uint32_t position1 = shastaAnchors->getPosition(anchorId1, orientedReadId);
                     
                     if (position1 > position0) {
                        const uint32_t gap = position1 - position0;
                        maxGap = std::max(maxGap, gap);
                     }
                 }
                 readSummary.middleAnchorGap = maxGap;

                 // Initial gap.
                 const shasta2::AnchorId anchorId0 = journey.front();
                 const uint32_t position0 = shastaAnchors->getPosition(anchorId0, orientedReadId);
                 readSummary.initialAnchorGap = position0;

                 // Final gap.
                 const shasta2::AnchorId anchorId1 = journey.back();
                 const uint32_t position1 = shastaAnchors->getPosition(anchorId1, orientedReadId);
                 if (readLength > position1) {
                    readSummary.finalAnchorGap = readLength - position1;
                 } else {
                    readSummary.finalAnchorGap = 0;
                 }
            }
        }


        // Create AnchorGraph.
        cout << timestamp << "Creating AnchorGraph..." << endl;
        shasta2::AnchorGraph anchorGraph(
            *shastaAnchors,
            journeys,
            options.minAnchorGraphEdgeCoverage);
        
        // Transitive Reduction.
        cout << timestamp << "Performing AnchorGraph Transitive Reduction..." << endl;
        anchorGraph.transitiveReduction(
            options.transitiveReductionMaxEdgeCoverage,
            options.transitiveReductionMaxDistance);

        // Create AssemblyGraph.
        cout << timestamp << "Creating AssemblyGraph..." << endl;
        {
            shasta2::AssemblyGraph assemblyGraph(
                *shastaAnchors,
                journeys,
                anchorGraph,
                options);

            // Save graphs for inspection/server.
            anchorGraph.save("AnchorGraph");

            // Final Assembly Step.
            cout << timestamp << "Simplifying and Assembling..." << endl;
            assemblyGraph.simplifyAndAssemble();
        } 
        cout << timestamp << "AssemblyGraph destroyed." << endl;

        // Write Read Summaries.
        cout << timestamp << "Writing ReadSummaries.csv..." << endl;
        {
             const string fileName = (fs::path(shastaOut) / "ReadSummaries.csv").string();
             ofstream csv(fileName);
             csv << "ReadId,Use for assembly,Is palindromic,Has high error rare,Palindromic rate,Initial marker error rate,Marker error rate,Initial anchor gap,Middle anchor gap,Final anchor gap,\n";
             for(shasta2::ReadId readId=0; readId<readSummaries.size(); readId++) {
                 const shasta2::ReadSummary& readSummary = readSummaries[readId];
                 csv << readId << "," <<
                     (readSummary.isInUse() ? "Yes" : "No") << "," <<
                     (readSummary.isPalindromic ? "Yes" : "No") <<
                     (readSummary.hasHighErrorRate ? "Yes" : "No") << "," <<
                     readSummary.palindromicRate << "," <<
                     readSummary.initialMarkerErrorRate << "," <<
                     readSummary.markerErrorRate << "," <<
                     readSummary.initialAnchorGap << "," <<
                     readSummary.middleAnchorGap << "," <<
                     readSummary.finalAnchorGap << ",\n";
             }
        }

        }
        cout << timestamp << "Shasta2 Assembly Completed." << endl;
        shasta2::performanceLog << timestamp << "Shasta2 Assembly Completed." << endl;

        // Restore original directory and Dinara prefix.
        // fs::current_path(originalCwd);
        assembler.largeDataFileNamePrefix = originalDinaraPrefix;
    }

}
