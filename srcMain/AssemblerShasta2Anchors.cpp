#include "AssemblerShasta2Anchors.hpp"
#include "Assembler.hpp"
#include "timestamp.hpp"
#include "shasta2/Anchor.hpp"
#include "shasta2/Reads.hpp"
#include "shasta2/Markers.hpp"
#include "shasta2/MarkerKmers.hpp"
#include "shasta2/MappedMemoryOwner.hpp"
#include "shasta2/Options.hpp"
#include "shasta2/Journeys.hpp"
#include "shasta2/AnchorGraph.hpp"
#include "shasta2/AssemblyGraph.hpp"
#include "shasta2/ReadSummary.hpp"

#include <thread>
#include <vector>
#include <memory>
#include <filesystem>
#include <iostream>
#include <algorithm>

using namespace std;
namespace fs = std::filesystem;

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
                        ami.ordinal = ordinal;
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

    void createShasta2Anchors(
        dinara::Assembler& assembler,
        uint64_t threadCount
    ) {
        cout << timestamp << "Creating Shasta2 Anchors from Marker Graph..." << endl;
        
        // Use a local MappedMemoryOwner for new shasta objects.
        shasta2::MappedMemoryOwner shastaOwner("./", assembler.largeDataPageSize);

        // 1. Access existing Reads data using shasta2 structures.
        auto shastaReadsPtr = make_shared<shasta2::Reads>();
        try {
            shastaReadsPtr->access(
                assembler.largeDataName("Reads"),
                assembler.largeDataName("ReadNames"),
                assembler.largeDataName("ReadIdsSortedByName")
            );
        } catch (const std::exception& e) {
             cout << timestamp << "Error accessing Shasta2 Reads: " << e.what() << endl;
             throw;
        }

        // 2. Setup Markers.
        // Shasta2 Markers expects "Markers" files in the current directory (or specific path).
        // Since Dinara stores them in Data/, we symlink them here for Shasta2 to find.
        try {
            const string markersPrefix = assembler.largeDataName("Markers");
            if(fs::exists(markersPrefix + ".toc")) {
                 if(fs::exists("Markers.toc")) fs::remove("Markers.toc");
                 fs::create_symlink(markersPrefix + ".toc", "Markers.toc");
            }
            if(fs::exists(markersPrefix + ".data")) {
                 if(fs::exists("Markers.data")) fs::remove("Markers.data");
                 fs::create_symlink(markersPrefix + ".data", "Markers.data");
            }
        } catch(const std::exception& e) {
            cout << timestamp << "Warning: Symlink creation for Markers failed: " << e.what() << endl;
        }

        shasta2::Markers* shastaMarkersPtr = nullptr;
        try {
             shastaMarkersPtr = new shasta2::Markers(
                 shastaOwner,
                 assembler.assemblerInfo->k,
                 shastaReadsPtr);
        } catch(const std::exception& e) {
             cout << timestamp << "Error accessing Shasta2 Markers: " << e.what() << endl;
             throw;
        }
        shasta2::Markers& shastaMarkers = *shastaMarkersPtr;


        // 3. Compute MarkerKmers.
        // We need dummy ReadSummaries because Dinara doesn't persist them, but Shasta2 API requires them.
        shasta2::MemoryMapped::Vector<shasta2::ReadSummary> readSummaries;
        readSummaries.createNew(
            "ReadSummaries", 
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
            "Shasta2Anchors",
            shastaOwner,
            *shastaReadsPtr,
            assembler.assemblerInfo->k,
            shastaMarkers,
            shastaMarkerKmers);

        populateAnchors(assembler, *shastaAnchors, threadCount);
        
        cout << timestamp << "Shasta2 Anchors created (" << shastaAnchors->size() << " anchors)." << endl;

        
        // 5. Downstream Assembly Pipeline (Journeys -> AnchorGraph -> AssemblyGraph).
        cout << timestamp << "Proceeding with Downstream Shasta2 Assembly..." << endl;

        // Initialize Options.
        int argc = 1;
        char* name = (char*)"dinara";
        char* argv[] = {name};
        shasta2::Options options(argc, argv);
        options.threadCount = threadCount;
        
        // Options: Anchor creation.
        options.minAnchorCoverage = 10;
        options.maxAnchorCoverage = 40;
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
        
        
        options.writeIntermediateAssemblyStages = true;
        
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

                 // Get markers and journey.
                 const auto orientedReadMarkers = shastaMarkers[orientedReadId.getValue()];
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
                     
                     const uint32_t ordinal0 = shastaAnchors->getOrdinal(anchorId0, orientedReadId);
                     const uint32_t ordinal1 = shastaAnchors->getOrdinal(anchorId1, orientedReadId);

                     const uint32_t position0 = orientedReadMarkers[ordinal0].position + kHalf;
                     const uint32_t position1 = orientedReadMarkers[ordinal1].position + kHalf;
                     
                     if (position1 > position0) {
                        const uint32_t gap = position1 - position0;
                        maxGap = std::max(maxGap, gap);
                     }
                 }
                 readSummary.middleAnchorGap = maxGap;

                 // Initial gap.
                 const shasta2::AnchorId anchorId0 = journey.front();
                 const uint32_t ordinal0 = shastaAnchors->getOrdinal(anchorId0, orientedReadId);
                 readSummary.initialAnchorGap = orientedReadMarkers[ordinal0].position + kHalf;

                 // Final gap.
                 const shasta2::AnchorId anchorId1 = journey.back();
                 const uint32_t ordinal1 = shastaAnchors->getOrdinal(anchorId1, orientedReadId);
                 if (readLength > (orientedReadMarkers[ordinal1].position + kHalf)) {
                    readSummary.finalAnchorGap = readLength - orientedReadMarkers[ordinal1].position - kHalf;
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
             ofstream csv("ReadSummaries.csv");
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

        cout << timestamp << "Shasta2 Assembly Completed." << endl;
    }

}
