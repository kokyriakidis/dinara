#include "AssemblerShasta2Anchors.hpp"
#include "Assembler.hpp"
#include "timestamp.hpp"
#include "iostream.hpp" // Dinara's iostream wrapper or std? main.cpp uses "iostream.hpp"

// Shasta 2 Integration
#include "shasta2/src/Anchor.hpp"
#include "shasta2/src/Reads.hpp"
#include "shasta2/src/Markers.hpp"
#include "shasta2/src/MarkerKmers.hpp"
#include "shasta2/src/MappedMemoryOwner.hpp"

#include <thread>
#include <vector>
#include <memory>
#include <filesystem>
#include <iostream>

#include "shasta2/src/Options.hpp"
#include "shasta2/src/Journeys.hpp"
#include "shasta2/src/AnchorGraph.hpp"
#include "shasta2/src/AssemblyGraph.hpp"
#include "shasta2/src/ReadSummary.hpp"

using namespace std;
namespace fs = std::filesystem;

namespace dinara {

    // Helper function to perform the actual conversion (multithreaded)
    static void populateAnchors(
        Assembler& assembler,
        shasta2::Anchors& anchors,
        uint64_t threadCount
    ) {
        if(threadCount == 0) threadCount = std::thread::hardware_concurrency();
        auto& vertices = assembler.markerGraph.vertices();
        uint64_t anchorCount = vertices.size();
        
        // Pass 1: Count sizes
        anchors.anchorMarkerInfos.beginPass1(anchorCount);
        
        std::vector<std::thread> threads;
        for(size_t t=0; t<threadCount; t++) {
            threads.emplace_back([&, t]() {
                size_t batchSize = anchorCount / threadCount;
                if(batchSize == 0) batchSize = 1;
                size_t begin = t * batchSize;
                size_t end = (t == threadCount - 1) ? anchorCount : begin + batchSize;
                if(begin >= anchorCount) return;
                
                for(size_t i=begin; i<end; ++i) {
                    const auto& vertex = vertices[i];
                    std::vector<dinara::OrientedReadId> uniqueReadIds;
                    uniqueReadIds.reserve(vertex.size());

                    for(dinara::MarkerId markerId : vertex) {
                        dinara::OrientedReadId dReadId;
                        uint32_t ordinal;
                        std::tie(dReadId, ordinal) = assembler.findMarkerId(markerId);
                        uniqueReadIds.push_back(dReadId);
                    }

                    // Sort and deduplicate to match logic in Pass 2
                    std::sort(uniqueReadIds.begin(), uniqueReadIds.end());
                    auto last = std::unique(uniqueReadIds.begin(), uniqueReadIds.end());
                    size_t uniqueCount = std::distance(uniqueReadIds.begin(), last);

                    // Increment allowed count
                    for(size_t k=0; k<uniqueCount; k++) {
                        anchors.anchorMarkerInfos.incrementCount(i);
                    }
                }
            });
        }
        for(auto& t : threads) t.join();
        threads.clear();
        
        // Pass 2: Populate Data
        anchors.anchorMarkerInfos.beginPass2();
        
        for(size_t t=0; t<threadCount; t++) {
            threads.emplace_back([&, t]() {
                size_t batchSize = anchorCount / threadCount;
                if(batchSize == 0) batchSize = 1;
                size_t begin = t * batchSize;
                size_t end = (t == threadCount - 1) ? anchorCount : begin + batchSize;
                if(begin >= anchorCount) return;
                
                for(size_t i=begin; i<end; ++i) {
                    const auto& vertex = vertices[i];
                    std::vector<shasta2::AnchorMarkerInfo> amis;
                    amis.reserve(vertex.size());

                    for(dinara::MarkerId markerId : vertex) {
                        dinara::OrientedReadId dReadId;
                        uint32_t ordinal;
                        // assembler.findMarkerId contains the logic to get ReadId/ordinal
                        std::tie(dReadId, ordinal) = assembler.findMarkerId(markerId);
                        
                        // Convert to shasta2 types using fromValue
                        shasta2::OrientedReadId sReadId = shasta2::OrientedReadId::fromValue(dReadId.getValue());
                        shasta2::AnchorMarkerInfo ami(sReadId, ordinal);
                        ami.positionInJourney = shasta2::invalid<uint32_t>;
                        amis.push_back(ami);
                    }

                    // Sort by OrientedReadId to satisfy Shasta2 requirements
                    std::sort(amis.begin(), amis.end());

                    // Deduplicate to satisfy strict strict inequality check in Anchor::check
                    // which implies only one marker per read per anchor.
                    auto last = std::unique(amis.begin(), amis.end(), 
                        [](const shasta2::AnchorMarkerInfo& a, const shasta2::AnchorMarkerInfo& b){
                            return a.orientedReadId == b.orientedReadId;
                        });
                    amis.erase(last, amis.end());

                    // Reverse AMIs because shasta2::VectorOfVectors::store fills BACKWARDS (from end to begin).
                    // We want them sorted [Small ... Large] in memory.
                    // Store fills: [ ... store(1) ... store(0) ] if we push 0 then 1?
                    // Hypothesis: store writes to (--count).
                    // If we want [Small, Medium, Large].
                    // Last slot filled first (with first store).
                    // So we must store LARGE first.
                    // So we must provide Descending order: Large, Medium, Small.
                    // Sorted amis are [Small, Medium, Large].
                    // Reversing gives [Large, Medium, Small].
                    std::reverse(amis.begin(), amis.end());

                    // Store sorted AMIs (now Descending)
                    for(const auto& ami : amis) {
                        anchors.anchorMarkerInfos.store(i, ami);
                    }
                }
            });
        }
        for(auto& t : threads) t.join();
        threads.clear();
        
        anchors.anchorMarkerInfos.endPass2();
        
        // Fill AnchorInfos
        // anchorInfos is already created by Anchors constructor.
        anchors.anchorInfos.resize(anchorCount);
        
        // Parallelize initialization of invalid kmers
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

    void createShasta2Anchors(
        dinara::Assembler& assembler,
        uint64_t threadCount
    ) {
        cout << timestamp << "Creating Shasta2 Anchors from Marker Graph..." << endl;
        
        // Use a local MappedMemoryOwner for shasta objects
        // Initialize with "./" prefix because empty prefix causes largeDataName to return empty string
        // (treated as anonymous mapping), which fails when accessing existing files like Markers.
        shasta2::MappedMemoryOwner shastaOwner("./", assembler.largeDataPageSize);

        // 1. Initialize shasta2 artifacts reusing dinara data paths
        auto shastaReadsPtr = make_shared<shasta2::Reads>();
        try {
            cout << "Accessing Shasta2 Reads with paths:" << endl;
            cout << "  Reads: " << assembler.largeDataName("Reads") << endl;
            cout << "  ReadNames: " << assembler.largeDataName("ReadNames") << endl;
            cout << "  ReadIdsSortedByName: " << assembler.largeDataName("ReadIdsSortedByName") << endl;

            shastaReadsPtr->access(
                assembler.largeDataName("Reads"),
                assembler.largeDataName("ReadNames"),
                assembler.largeDataName("ReadIdsSortedByName")
            );
        } catch (const std::exception& e) {
             cout << "Error accessing Shasta2 Reads: " << e.what() << endl;
             // If ReadIdsSortedByName is missing, we might cope or fail.
             // But for now, let's rethrow or exit.
             throw;
        }

        // Shasta2 Markers expects "Markers" in the current directory.
        // Create symlinks to dinara's Data/Markers.
        try {
            string markersPrefix = assembler.largeDataName("Markers");
            cout << "Symlinking Markers from " << markersPrefix << "..." << endl;
            if(fs::exists(markersPrefix + ".toc")) {
                 if(fs::exists("Markers.toc")) fs::remove("Markers.toc");
                 fs::create_symlink(markersPrefix + ".toc", "Markers.toc");
            }
            if(fs::exists(markersPrefix + ".data")) {
                 if(fs::exists("Markers.data")) fs::remove("Markers.data");
                 fs::create_symlink(markersPrefix + ".data", "Markers.data");
            }
        } catch(const std::exception& e) {
            cout << "Warning: Symlink creation for Markers failed: " << e.what() << endl;
        }


        // Shasta2 Markers
        // We use check if Markers.toc exists before accessing.
        cout << "Checking Markers files..." << endl;
        if(fs::exists("Markers.toc")) cout << "Markers.toc exists." << endl;
        else cout << "Markers.toc DOES NOT exist." << endl;

        cout << "Accessing Shasta2 Markers..." << endl;
        shasta2::Markers* shastaMarkersPtr = nullptr;
        try {
             shastaMarkersPtr = new shasta2::Markers(
                 shastaOwner, // Correct owner
                 assembler.assemblerInfo->k,
                 shastaReadsPtr);
        } catch(const std::exception& e) {
             cout << "Error accessing Shasta2 Markers: " << e.what() << endl;
             throw;
        }
        shasta2::Markers& shastaMarkers = *shastaMarkersPtr;


        // Shasta2 MarkerKmers
        // We must COMPUTE them because they likely don't exist as shasta artifacts.
        // To compute, we need ReadSummaries. Dinara doesn't keep them, so we create valid dummies.
        cout << "Creating dummy ReadSummaries for MarkerKmers computation..." << endl;
        shasta2::MemoryMapped::Vector<shasta2::ReadSummary> readSummaries;
        readSummaries.createNew(
            "tmp-ReadSummaries", 
            assembler.largeDataPageSize
        );
        readSummaries.resize(shastaReadsPtr->readCount()); 
        // Default constructor of ReadSummary sets isPalindromic=false, hasHighErrorRate=false.
        // So all reads are used.

        cout << "Computing Shasta2 MarkerKmers..." << endl;
        shasta2::MarkerKmers shastaMarkerKmers(
             assembler.assemblerInfo->k,
             shastaOwner,      // owner
             *shastaReadsPtr,  // reads
             readSummaries,    // readSummaries
             shastaMarkers,    // markers
             threadCount       // threadCount
        );
        readSummaries.remove(); // Clean up temp file

        // 2. Create Shasta2 Anchors
        // This will create new files Shasta2Anchors-AnchorMarkerInfos etc.
        cout << "Creating Shasta2 Anchors object..." << endl;
        auto shastaAnchors = make_shared<shasta2::Anchors>(
            "Shasta2Anchors", // baseName
            shastaOwner,
            *shastaReadsPtr,
            assembler.assemblerInfo->k,
            shastaMarkers,
            shastaMarkerKmers);

        // 3. Populate Anchors from Marker Graph
        // This logic remains similar to original plan
        populateAnchors(assembler, *shastaAnchors, threadCount);
        cout << timestamp << "Shasta2 Anchors created." << endl;
        
        // DEBUG: Print anchor statistics
        cout << "DEBUG: Anchor count: " << shastaAnchors->size() << endl;
        cout << "DEBUG: Total AnchorMarkerInfos: " << shastaAnchors->anchorMarkerInfos.totalSize() << endl;
        
        // Sample some anchors to verify content
        if (shastaAnchors->size() > 0) {
            size_t sampleIdx = 0;
            auto sampleAnchor = (*shastaAnchors)[sampleIdx];
            cout << "DEBUG: Anchor " << sampleIdx << " has " << sampleAnchor.size() << " marker infos" << endl;
            for (size_t j = 0; j < std::min((size_t)3, sampleAnchor.size()); j++) {
                const auto& ami = sampleAnchor[j];
                cout << "  MarkerInfo " << j << ": readId=" << ami.orientedReadId.getValue() 
                     << " ordinal=" << ami.ordinal << endl;
            }
        }
        
        // DEBUG: Anchor coverage histogram - write to file
        std::map<uint64_t, uint64_t> coverageHistogram;
        for (size_t i = 0; i < shastaAnchors->size(); i++) {
            uint64_t coverage = (*shastaAnchors)[i].size();
            coverageHistogram[coverage]++;
        }
        
        // Write to file
        std::ofstream histFile("AnchorCoverageHistogram.csv");
        histFile << "Coverage,AnchorCount\n";
        for (const auto& [cov, count] : coverageHistogram) {
            histFile << cov << "," << count << "\n";
        }
        histFile.close();
        cout << "DEBUG: Anchor coverage histogram written to AnchorCoverageHistogram.csv" << endl;
        
        // Also print summary to console
        cout << "DEBUG: Coverage distribution summary:" << endl;
        for (const auto& [cov, count] : coverageHistogram) {
            if (cov <= 10 || count > 100) {
                cout << "  Coverage " << cov << ": " << count << " anchors" << endl;
            }
        }


        // 4. Downstream Assembly
        cout << timestamp << "Proceeding with Downstream Shasta2 Assembly..." << endl;

        // Initialize Options
        int argc = 1;
        char* name = (char*)"dinara";
        char* argv[] = {name};
        shasta2::Options options(argc, argv);
        options.threadCount = threadCount;
        
        // Map relevant dinara options if necessary, e.g. min/max coverage
        // Lower minAnchorGraphEdgeCoverage to allow more edges (default is 6)
        options.minAnchorGraphEdgeCoverage = 2;
        
        // Create Journeys
        cout << "Creating Journeys..." << endl;
        shasta2::Journeys journeys(
            2 * shastaReadsPtr->readCount(),
            shastaAnchors,
            threadCount,
            shastaOwner);
        
        // DEBUG: Check journey statistics
        uint64_t nonEmptyJourneys = 0;
        uint64_t totalJourneyLength = 0;
        for (uint64_t i = 0; i < 2 * shastaReadsPtr->readCount(); i++) {
            auto journey = journeys[shasta2::OrientedReadId::fromValue(i)];
            if (journey.size() > 0) {
                nonEmptyJourneys++;
                totalJourneyLength += journey.size();
            }
        }
        cout << "DEBUG: Non-empty journeys: " << nonEmptyJourneys << " / " << 2 * shastaReadsPtr->readCount() << endl;
        cout << "DEBUG: Total journey length: " << totalJourneyLength << endl;
        if (nonEmptyJourneys > 0) {
            cout << "DEBUG: Average journey length: " << (double)totalJourneyLength / nonEmptyJourneys << endl;
        }

        // Create AnchorGraph
        cout << "Creating AnchorGraph..." << endl;
        shasta2::AnchorGraph anchorGraph(
            *shastaAnchors,
            journeys,
            options.minAnchorGraphEdgeCoverage);
        
        // Transitive Reduction of AnchorGraph
        // This is a critical step in shasta2 assembly pipeline.
        cout << "performing AnchorGraph Transitive Reduction..." << endl;
        anchorGraph.transitiveReduction(
            options.transitiveReductionMaxEdgeCoverage,
            options.transitiveReductionMaxDistance);

        // Create AssemblyGraph
        cout << "Creating AssemblyGraph..." << endl;
        shasta2::AssemblyGraph assemblyGraph(
            *shastaAnchors,
            journeys,
            anchorGraph,
            options);

        // Save AnchorGraph for HTTP server access
        cout << "Saving AnchorGraph..." << endl;
        anchorGraph.save("AnchorGraph");

        // Simplify and Assemble
        // This method handles the full assembly process including outputting GFA/Fasta.
        cout << "Simplifying and Assembling..." << endl;
        assemblyGraph.simplifyAndAssemble();

        cout << timestamp << "Shasta2 Assembly Completed." << endl;
    }

}

