// Dinara.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "ConsensusCaller.hpp"
#include "compressAlignment.hpp"
#include "Coverage.hpp"
#include "dset64-gccAtomic.hpp"
#include "extractKmer.hpp"
#include "PeakFinder.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "chrono.hpp"
#include <map>
#include <queue>



// Loop over all alignments in the read graph
// to create vertices of the global marker graph.
// Throw away vertices with coverage (number of markers)
// less than minCoverage or more than maxCoverage.
// Also throw away "bad" vertices - that is, vertices
// with more than one marker on the same oriented read.
void Assembler::createMarkerGraphVertices(

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
)
{

    // Flag to control debug output.
    // Only turn on for a very small test run.
    const bool debug = false;

    // using VertexId = MarkerGraph::VertexId;
    // using CompressedVertexId = CompressedGlobalMarkerGraphVertexId;

    const auto tBegin = steady_clock::now();
    performanceLog << timestamp << "Begin computing marker graph vertices." << endl;

    // Check that we have what we need.
    reads->checkReadsAreOpen();
    reads->checkReadFlagsAreOpen();
    reads->checkReadFlagsAreOpen();
    // DINARA_ASSERT(kmerChecker); // Not needed for simple vertex creation
    checkMarkersAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    DINARA_ASSERT(compressedAlignments.isOpen());


    // Store parameters so they are accessible to the threads.
    auto& data = createMarkerGraphVerticesData;
    data.allowDuplicateMarkers = allowDuplicateMarkers;
    data.minCoveragePerStrand = minCoveragePerStrand;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Initialize computation of the global marker graph.
    data.orientedMarkerCount = markers->totalSize();

    data.disjointSetTable.createNew(
        largeDataName("tmp-DisjointSetTable"),
        largeDataPageSize);
    // DisjointSets data structure needs an additional 64 bits per entry, in order to implement
    // a lock-free, union-find operation. You can find more information in dset64-gccatomic.hpp.
    // Once the set representatives have been found, we have no need for these extra 64 bits per entry.
    //
    // We allocate twice as much space in data.disjointSetTable so that the underlying memory
    // can be used as an array of 128 bit integers of size data.orientedMarkerCount. This allows
    // us to compact the data in-place, there by reducing memory usage.
    data.disjointSetTable.reserveAndResize(data.orientedMarkerCount * 2);

    // Have DisjointSets use the memory allocated in and managed by data.disjointSetTable.
    data.disjointSetsPointer = std::make_shared<DisjointSets>(
        reinterpret_cast<DisjointSets::Aint*>(data.disjointSetTable.begin()),
        data.orientedMarkerCount
    );



    // Update the disjoint set data structure for each alignment
    // in the read graph.
    performanceLog << timestamp << "Disjoint set computation begins." << endl;
    size_t batchSize = 10000;
    setupLoadBalancing(readGraph.edges.size(), batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction1, threadCount);
    performanceLog << timestamp << "Disjoint set computation completed." << endl;



    // Find the disjoint set that each oriented marker was assigned to.
    // Iterate till each marker has its set representative populated in the parent (lower 64 bits)
    uint64_t pass = 1;
    do {
        (data.disjointSetsPointer)->parentUpdated = 0;
        performanceLog << "    " << timestamp << " Iteration  " << pass << endl;
        setupLoadBalancing(data.orientedMarkerCount, batchSize);
        runThreads(&Assembler::createMarkerGraphVerticesThreadFunction2, threadCount);
        performanceLog << "    " << timestamp << " Updated parent of - " << (data.disjointSetsPointer)->parentUpdated << " entries." << endl;
        pass++;
    } while ((data.disjointSetsPointer)->parentUpdated > 0 && pass <= 10);

    if (pass > 10) {
        // This should never happen. Even in a highly parallel environment (128 threads), convergence happens
        // in 2 or 3 passes. It's definitely worth investigating if this ever happens.
        string errorMsg = "DisjointSets parent information did not converge in " + to_string(pass) + " iterations.";
        throw runtime_error(errorMsg);
    }


    performanceLog << timestamp << "Verifying convergence of parent information." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction21, threadCount);
    performanceLog << timestamp << "Done verifying convergence of parent information." << endl;


    // data.disjointSetTable now has the correct set representative for entry N at location 2*N.
    // That's because DisjointSets stores parent information in the lower 64 bits of the 128 bits
    // it uses for each entry. Since we only care about these bits, we can compact data.disjointSetTable
    // and free up half the memory.
    // This bit seems tricky to parallelize. It's not worth the effort as this is pretty fast as is.
    performanceLog << timestamp << "Compacting the Disjoint Set data-structure." << endl;
    for(uint64_t i=0; i<data.orientedMarkerCount; i++) {
        data.disjointSetTable[i] = data.disjointSetTable[2*i];
    }
    data.disjointSetTable.resize(data.orientedMarkerCount);
    data.disjointSetTable.unreserve();
    performanceLog << timestamp << "Done compacting the Disjoint Set data-structure." << endl;

    // Don't need the DisjointSets data-structure any more.
    data.disjointSetsPointer = 0;



    // Debug output.
    if(debug) {
        createMarkerGraphVerticesDebug1(0);
    }



    // Count the number of markers in each disjoint set
    // and store it in data.workArea.
    // We don't want to combine this with the previous block
    // because it would significantly increase the peak memory usage.
    // This way, we allocate data.workArea only after compacting
    // data.disjointSetTable.
    performanceLog << timestamp << "Counting the number of markers in each disjoint set." << endl;
    data.workArea.createNew(
        largeDataName("tmp-WorkArea"),
        largeDataPageSize);
    data.workArea.reserveAndResize(data.orientedMarkerCount);
    fill(data.workArea.begin(), data.workArea.end(), 0ULL);
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction3, threadCount);



    // Debug output.
    if(debug) {
        ofstream out("WorkArea-initial-count.csv");
        for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
            out << markerId << "," << data.workArea[markerId] << "\n";
        }

    }



    // At this point, data.workArea contains the number of oriented markers in
    // each disjoint set.
    // Compute a histogram of this distribution and write it to a csv file.
    {
        vector<uint64_t> histogram;
        for(MarkerGraph::VertexId i=0; i<data.orientedMarkerCount; i++) {
            const MarkerGraph::VertexId markerCount = data.workArea[i];
            if(markerCount == 0) {
                continue;
            }
            if(markerCount >= histogram.size()) {
                histogram.resize(markerCount+1, 0);
            }
            ++histogram[markerCount];
        }

        // Store the disjoint sets histogram in a MemoryMapped::Vector.
        // This is used when flagging primary marker graph edges for Mode 3 assembly.
        // This stored pairs(coverage, frequency).
        // Only pairs where the frequency is not zero are stored.
        {
            markerGraph.disjointSetsHistogram.createNew(
                largeDataName("DisjointSetsHistogram"),
                largeDataPageSize);
            for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
                const uint64_t frequency = histogram[coverage];
                if(frequency) {
                    markerGraph.disjointSetsHistogram.push_back({coverage, frequency});
                }
            }
        }

        ofstream csv("DisjointSetsHistogram.csv");
        csv << "Coverage,Frequency\n";
        for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
            const uint64_t frequency = histogram[coverage];
            if(frequency) {
                csv << coverage << "," << frequency << "\n";
            }
        }

        if (minCoverage == 0) {
            try {
                dinara::PeakFinder p;
                p.findPeaks(histogram);
                minCoverage = p.findXCutoff(histogram, peakFinderMinAreaFraction, peakFinderAreaStartIndex);
                cout << "Automatically selected value of MarkerGraph.minCoverage "
                    "is " << minCoverage << endl;
            }
            catch (PeakFinderException& e){
                minCoverage = 5;
                cout <<
                    "Unable to automatically select MarkerGraph.minCoverage. "
                    "No significant cutoff found in disjoint sets size distribution. "
                    "Observed peak has percent total area of " << e.observedPercentArea << endl <<
                    "minPercentArea is " << e.minPercentArea << endl <<
                    "See DisjointSetsHistogram.csv."
                    "Using MarkerGraph.minCoverage = " << minCoverage << endl;
            }
        }
    }
    // Store the value of minCoverage actually used.
    assemblerInfo->markerGraphMinCoverageUsed = minCoverage;



    // At this point, data.workArea contains the number of oriented markers in
    // each disjoint set.
    // Replace it with a new numbering, counting only disjoint sets
    // with size not less than minCoverage and not greater than maxCoverage.
    // Note that this numbering is not yet the final vertex numbering,
    // as we will later remove "bad" vertices
    // (vertices with more than one marker on the same read).
    // This block is recursive and cannot be multithreaded.
    performanceLog << timestamp << "Renumbering the disjoint sets." << endl;
    MarkerGraph::VertexId newDisjointSetId = 0ULL;
    for(MarkerGraph::VertexId oldDisjointSetId=0;
        oldDisjointSetId<data.orientedMarkerCount; ++oldDisjointSetId) {
        auto& w = data.workArea[oldDisjointSetId];
        const MarkerGraph::VertexId markerCount = w;
        if(markerCount<minCoverage || markerCount>maxCoverage) {
            w = MarkerGraph::invalidVertexId;
        } else {
            w = newDisjointSetId;
            ++newDisjointSetId;
        }
    }
    const auto disjointSetCount = newDisjointSetId;
    cout << "Kept " << disjointSetCount << " disjoint sets with coverage in the requested range." << endl;



    // Debug output.
    if(debug) {
        ofstream out("WorkArea-initial-renumbering.csv");
        for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
            out << markerId << "," << data.workArea[markerId] << "\n";
        }
    }


    // Reassign vertices to disjoint sets using this new numbering.
    // Vertices assigned to no disjoint set will store MarkerGraph::invalidVertexId.
    // This could be multithreaded if necessary.
    performanceLog << timestamp << "Assigning vertices to renumbered disjoint sets." << endl;
    for(MarkerGraph::VertexId markerId=0;
        markerId<data.orientedMarkerCount; ++markerId) {
        auto& d = data.disjointSetTable[markerId];
        const auto oldId = d;
        const auto newId = data.workArea[oldId];
        d = newId;
    }
    // We no longer need the workArea.
    data.workArea.remove();



    // Debug output.
    if(debug) {
        createMarkerGraphVerticesDebug1(1);
    }



    // At this point, data.disjointSetTable contains, for each oriented marker,
    // the disjoint set that the oriented marker is assigned to,
    // and all disjoint sets have a number of markers in the requested range.
    // Vertices not assigned to any disjoint set store a disjoint set
    // equal to MarkerGraph::invalidVertexId.



    // Gather the markers in each disjoint set.
    data.disjointSetMarkers.createNew(
        largeDataName("tmp-DisjointSetMarkers"),
        largeDataPageSize);
    performanceLog << timestamp << "Gathering markers in disjoint sets, pass1." << endl;
    data.disjointSetMarkers.beginPass1(disjointSetCount);
    performanceLog << timestamp << "Processing " << data.orientedMarkerCount << " oriented markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction4, threadCount);
    performanceLog << timestamp << "Gathering markers in disjoint sets, pass2." << endl;
    data.disjointSetMarkers.beginPass2();
    performanceLog << timestamp << "Processing " << data.orientedMarkerCount << " oriented markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction5, threadCount);
    data.disjointSetMarkers.endPass2();



    // Sort the markers in each disjoint set.
    performanceLog << timestamp << "Sorting the markers in each disjoint set." << endl;
    setupLoadBalancing(disjointSetCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction6, threadCount);



    // Debug output.
    if(debug) {
        const uint64_t stage = 2;

        vector<uint64_t> histogram;
        for(uint64_t i=0; i<data.disjointSetMarkers.size(); i++) {
            const uint64_t size = data.disjointSetMarkers.size(i);
            if(histogram.size() <= size) {
                histogram.resize(size+1, 0);
            }
            ++histogram[size];
        }
        ofstream csv("Debug-DisjointSets-Histogram-" + to_string(stage) + ".csv");
        csv << "Size,Frequency\n";
        for(uint64_t size=0; size<histogram.size(); size++) {
            csv << size << "," << histogram[size] << "\n";
        }

    }



    // Flag "bad" disjoint sets for which we don't want to
    // create marker graph vertices. A disjoint set can be flagged
    // as bad for one of two reasons:
    // - It contains more than one marker on the same oriented read.
    // - It does not contain at least minCoveragePerStrand supporting
    //   oriented reads on each strand.
    data.isBadDisjointSet.createNew(
        largeDataName("tmp-IsBadDisjointSet"),
        largeDataPageSize);
    data.isBadDisjointSet.reserveAndResize(disjointSetCount);
    performanceLog << timestamp << "Flagging bad disjoint sets." << endl;
    setupLoadBalancing(disjointSetCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction7, threadCount);
    const size_t badDisjointSetCount = std::count(
        data.isBadDisjointSet.begin(), data.isBadDisjointSet.end(), true);
    cout << "Found " << badDisjointSetCount << " disjoint sets "
        "with more than one marker on a single oriented read "
        "or with less than " << minCoveragePerStrand <<
        " supporting oriented reads on each strand." << endl;



    // Renumber the disjoint sets again, this time without counting the ones marked as bad.
    performanceLog << timestamp << "Renumbering disjoint sets to remove the bad ones." << endl;
    data.workArea.createNew(
        largeDataName("tmp-WorkArea"),
        largeDataPageSize);
    data.workArea.reserveAndResize(disjointSetCount);
    newDisjointSetId = 0ULL;
    for(MarkerGraph::VertexId oldDisjointSetId=0;
        oldDisjointSetId<disjointSetCount; ++oldDisjointSetId) {
        auto& w = data.workArea[oldDisjointSetId];
        if(data.isBadDisjointSet[oldDisjointSetId]) {
            w = MarkerGraph::invalidVertexId;
        } else {
            w = newDisjointSetId;
            ++newDisjointSetId;
        }
    }
    DINARA_ASSERT(newDisjointSetId + badDisjointSetCount == disjointSetCount);



    // Debug output.
    if(debug) {
        ofstream out("WorkArea-final-renumbering.csv");
        for(MarkerId markerId=0; markerId<disjointSetCount; markerId++) {
            out << markerId << "," << data.workArea[markerId] << "\n";
        }
    }


    // Compute the final disjoint set number for each marker.
    // That becomes the vertex id assigned to that marker.
    // This could be multithreaded.
    performanceLog << timestamp << "Assigning vertex ids to markers." << endl;
    markerGraph.vertexTable.createNew(
        largeDataName("MarkerGraphVertexTable"),
        largeDataPageSize);
    markerGraph.vertexTable.reserveAndResize(data.orientedMarkerCount);
    for(MarkerGraph::VertexId markerId=0;
        markerId<data.orientedMarkerCount; ++markerId) {
        auto oldValue = data.disjointSetTable[markerId];
        if(oldValue == MarkerGraph::invalidVertexId) {
            markerGraph.vertexTable[markerId] = MarkerGraph::invalidCompressedVertexId;
        } else {
        	markerGraph.vertexTable[markerId] = data.workArea[oldValue];
        }
    }

    data.workArea.remove();
    data.disjointSetTable.remove();


    // Store the disjoint sets that are not marked bad.
    // Each corresponds to a vertex of the global marker graph.
    // This could be multithreaded.
    performanceLog << timestamp << "Gathering the markers of each vertex of the marker graph." << endl;
    markerGraph.constructVertices();
    markerGraph.vertices().createNew(
        largeDataName("MarkerGraphVertices"),
        largeDataPageSize);
    for(MarkerGraph::VertexId oldDisjointSetId=0;
        oldDisjointSetId<disjointSetCount; ++oldDisjointSetId) {
        if(data.isBadDisjointSet[oldDisjointSetId]) {
            continue;
        }
        markerGraph.vertices().appendVector();
        const auto markers = data.disjointSetMarkers[oldDisjointSetId];
        for(const MarkerId markerId: markers) {
            markerGraph.vertices().append(markerId);
        }
    }
    markerGraph.vertices().unreserve();

    data.isBadDisjointSet.remove();
    data.disjointSetMarkers.remove();


    // Check that the data structures we created are consistent with each other.
    // This could be expensive. Remove when we know this code works.
    // cout << timestamp << "Checking marker graph vertices." << endl;
    // checkMarkerGraphVertices(minCoverage, maxCoverage);




    // Debug output.
    if(debug) {
        const uint64_t stage = 3;
        ofstream csv1("Debug-Vertices-" + to_string(stage) + ".csv");
        for(uint64_t i=0 ;i<markerGraph.vertices().size(); i++) {
            const auto v = markerGraph.vertices()[i];
            csv1 << v.size() << ",";
            for(const MarkerId markerId: v) {
                OrientedReadId orientedReadId;
                uint32_t ordinal;
                tie(orientedReadId, ordinal) = findMarkerId(markerId);
                csv1 << orientedReadId << "-" << ordinal << ",";
            }
            csv1 << "\n";
        }

        // Also create a histogram.
        vector<uint64_t> histogram;
        for(uint64_t i=0 ;i<markerGraph.vertices().size(); i++) {
            const auto v = markerGraph.vertices()[i];
            const uint64_t size = v.size();
            if(histogram.size() <= size) {
                histogram.resize(size+1, 0);
            }
            ++histogram[size];
        }
        ofstream csv3("Debug-Vertices-Histogram-" + to_string(stage) + ".csv");
        csv3 << "Size,Frequency\n";
        for(uint64_t size=0; size<histogram.size(); size++) {
            csv3 << size << "," << histogram[size] << "\n";
        }
    }



    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    performanceLog << timestamp << "Computation of global marker graph vertices ";
    performanceLog << "completed in " << tTotal << " s." << endl;

}



void Assembler::createMarkerGraphVerticesThreadFunction1(uint64_t)
{

    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;

    auto& data = createMarkerGraphVerticesData;
    const std::shared_ptr<DisjointSets> disjointSetsPointer = data.disjointSetsPointer;

    const auto& storedAlignments = compressedAlignments;
    uint64_t alignmentId;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // We process read graph edges in pairs.
        // In each pair, the second edge is the reverse complement of the first.
        DINARA_ASSERT((begin%2) == 0);
        DINARA_ASSERT((end%2) == 0);

        for(size_t i=begin; i!=end; i+=2) {

            // Get the oriented read ids we want to align.
            array<OrientedReadId, 2> orientedReadIds;
            const ReadGraphEdge& readGraphEdge = readGraph.edges[i];
            alignmentId = readGraphEdge.alignmentId;

            // Check that the next edge is the reverse complement of
            // this edge.
            {
                const ReadGraphEdge& readGraphNextEdge = readGraph.edges[i + 1];
                array<OrientedReadId, 2> nextEdgeOrientedReadIds = readGraphNextEdge.orientedReadIds;
                nextEdgeOrientedReadIds[0].flipStrand();
                nextEdgeOrientedReadIds[1].flipStrand();
                DINARA_ASSERT(nextEdgeOrientedReadIds == readGraphEdge.orientedReadIds);
            }

            // If the edge is flagged as crossing strands, skip it.
            if(readGraphEdge.crossesStrands) {
                continue;
            }

            // If the edge has an alignment flagged as inconsistent, skip it.
            if(readGraphEdge.hasInconsistentAlignment) {
                continue;
            }

            orientedReadIds = readGraphEdge.orientedReadIds;
            DINARA_ASSERT(orientedReadIds[0] < orientedReadIds[1]);

            // If either of the reads is flagged chimeric, skip it.
            if( reads->getFlags(orientedReadIds[0].getReadId()).isChimeric ||
                reads->getFlags(orientedReadIds[1].getReadId()).isChimeric) {
                continue;
            }

            // Sanity check.
            DINARA_ASSERT(alignmentData[alignmentId].info.isInReadGraph);

            // Decompress this alignment.
            span<const char> compressedAlignment = storedAlignments[alignmentId];
            dinara::decompress(compressedAlignment, alignment);


            // In the global marker graph, merge pairs
            // of aligned markers.
            for(const auto& p: alignment.ordinals) {
                const uint32_t ordinal0 = p[0];
                const uint32_t ordinal1 = p[1];
                const MarkerId markerId0 = getMarkerId(orientedReadIds[0], ordinal0);
                const MarkerId markerId1 = getMarkerId(orientedReadIds[1], ordinal1);
                disjointSetsPointer->unite(markerId0, markerId1);

                // Also merge the reverse complemented markers.
                // This guarantees that the marker graph remains invariant
                // under strand swap.
                disjointSetsPointer->unite(
                    findReverseComplement(markerId0),
                    findReverseComplement(markerId1));
            }
        }
    }

}



void Assembler::createMarkerGraphVerticesThreadFunction2(uint64_t)
{
    DisjointSets& disjointSets = *createMarkerGraphVerticesData.disjointSetsPointer;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            // Update parent information.
            disjointSets.find(i, true);
        }
    }
}

void Assembler::createMarkerGraphVerticesThreadFunction21(uint64_t)
{
    DisjointSets& disjointSets = *createMarkerGraphVerticesData.disjointSetsPointer;
    const auto& disjointSetTable = createMarkerGraphVerticesData.disjointSetTable;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            // Verify that parent has been populated (I.e convergence happened after some iterations of
            // createMarkerGraphVerticesThreadFunction2)
            DINARA_ASSERT(disjointSets.parent(i) == disjointSets.find(i));
            // Verify that reinterpreting DisjointSets data as a vector of uint64_t will work as expected.
            DINARA_ASSERT(disjointSets.parent(i) == disjointSetTable[2*i]);
        }
    }
}


void Assembler::createMarkerGraphVerticesThreadFunction3(uint64_t)
{
    const auto& disjointSetTable = createMarkerGraphVerticesData.disjointSetTable;
    auto& workArea = createMarkerGraphVerticesData.workArea;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            const uint64_t disjointSetId = disjointSetTable[i];

            // Increment the set size in a thread-safe way.
            __sync_fetch_and_add(&workArea[disjointSetId], 1ULL);
        }
    }
}



void Assembler::createMarkerGraphVerticesThreadFunction4(uint64_t)
{
    createMarkerGraphVerticesThreadFunction45(4);
}



void Assembler::createMarkerGraphVerticesThreadFunction5(uint64_t)
{
    createMarkerGraphVerticesThreadFunction45(5);
}



void Assembler::createMarkerGraphVerticesThreadFunction6(uint64_t)
{
    auto& disjointSetMarkers = createMarkerGraphVerticesData.disjointSetMarkers;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerGraph::VertexId i=begin; i!=end; ++i) {
            auto markers = disjointSetMarkers[i];
            sort(markers.begin(), markers.end());
        }
    }
}



// Flag "bad" disjoint sets for which we don't want to
// create marker graph vertices. A disjoint set can be flagged
// as bad for one of two reasons:
// - It contains more than one marker on the same oriented read
//   (But this check if suppressed if allowDuplicateMarkers is set).
// - It does not contain at least minCoveragePerStrand supporting
//   oriented reads on each strand.
void Assembler::createMarkerGraphVerticesThreadFunction7(uint64_t)
{
    const auto& disjointSetMarkers = createMarkerGraphVerticesData.disjointSetMarkers;
    auto& isBadDisjointSet = createMarkerGraphVerticesData.isBadDisjointSet;
    const auto allowDuplicateMarkers = createMarkerGraphVerticesData.allowDuplicateMarkers;
    const auto minCoveragePerStrand = createMarkerGraphVerticesData.minCoveragePerStrand;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerGraph::VertexId disjointSetId=begin; disjointSetId!=end; ++disjointSetId) {
            auto markers = disjointSetMarkers[disjointSetId];
            const size_t markerCount = markers.size();
            DINARA_ASSERT(markerCount > 0);
            isBadDisjointSet[disjointSetId] = false;
            if(markerCount == 1) {
                if(1 < minCoveragePerStrand) {
                    isBadDisjointSet[disjointSetId] = true;
                }
                continue;
            }
            array<uint64_t, 2> countByStrand = {0, 0};
            for(size_t j=0; j<markerCount; j++) {
                const MarkerId& markerId = markers[j];
                OrientedReadId orientedReadId;
                tie(orientedReadId, ignore) = findMarkerId(markerId);
                ++countByStrand[orientedReadId.getStrand()];

                if((not allowDuplicateMarkers) and j > 0) {
                    const MarkerId& previousMarkerId = markers[j-1];
                    OrientedReadId previousOrientedReadId;
                    tie(previousOrientedReadId, ignore) = findMarkerId(previousMarkerId);
                    if(orientedReadId.getReadId() == previousOrientedReadId.getReadId()) {
                        isBadDisjointSet[disjointSetId] = true;
                        break;
                    }
                }
            }

            // If we did not flag it above due to more than one marker on the
            // same oriented read, check it for sufficient coverage on each strand.
            if(not isBadDisjointSet[disjointSetId]) {
                isBadDisjointSet[disjointSetId] =
                    (countByStrand[0] < minCoveragePerStrand) or
                    (countByStrand[1] < minCoveragePerStrand);
            }
        }
    }
}



void Assembler::createMarkerGraphVerticesThreadFunction45(int value)
{
    DINARA_ASSERT(value==4 || value==5);
    const auto& disjointSetTable = createMarkerGraphVerticesData.disjointSetTable;
    auto& disjointSetMarkers = createMarkerGraphVerticesData.disjointSetMarkers;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            const uint64_t disjointSetId = disjointSetTable[i];
            if(disjointSetId == MarkerGraph::invalidVertexId) {
                continue;
            }
            if(value == 4) {
                disjointSetMarkers.incrementCountMultithreaded(disjointSetId);
            } else {
                disjointSetMarkers.storeMultithreaded(disjointSetId, i);
            }
        }
   }
}



void Assembler::createMarkerGraphVerticesDebug1(uint64_t stage)
{
    const auto& data = createMarkerGraphVerticesData;

    // Dump the disjoint sets table.
    ofstream csv1("Debug-DisjointSets-Table-" + to_string(stage) + ".csv");
    csv1 << "MarkerId,ReadId,Strand,Ordinal,DisjointSet\n";
    MarkerId markerIdCheck = 0;
    const ReadId readCount = getReads().readCount();
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const uint64_t thisOrientedReadMarkerCount = markers->size(orientedReadId.getValue());
            for(uint32_t ordinal=0; ordinal<thisOrientedReadMarkerCount; ordinal++) {
                const MarkerId markerId = getMarkerId(orientedReadId, ordinal);
                DINARA_ASSERT(markerId == markerIdCheck++);
                csv1 <<
                    markerId << "," <<
                    readId << "," <<
                    strand << "," <<
                    ordinal << "," <<
                    data.disjointSetTable[markerId] << "\n";
            }
        }
    }

    // Gather the markers in disjoint sets and sort them so
    // the result is reproducible. This is expensive and can only be done like
    // this in debug code.
    std::map<MarkerId, vector<MarkerId> > m;
    for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
        const auto disjointSetId = data.disjointSetTable[markerId];
        if(disjointSetId != MarkerGraph::invalidVertexId) {
            m[disjointSetId].push_back(markerId);
        }
    }
    vector< vector<MarkerId> > v;
    for(const auto& p: m) {
        v.push_back(p.second);
    }
    sort(v.begin(), v.end());
    ofstream csv2("Debug-DisjointSets" + to_string(stage) + ".csv");
    for(const auto& s: v) {
        csv2 << v.size() << ",";
        for(const MarkerId markerId: s) {
            OrientedReadId orientedReadId;
            uint32_t ordinal;
            tie(orientedReadId, ordinal) = findMarkerId(markerId);
            csv2 << orientedReadId << "-" << ordinal << ",";
        }
        csv2 << "\n";
    }

    // Also create a histogram.
    vector<uint64_t> histogram;
    for(const auto& s: v) {
        const uint64_t size = s.size();
        if(histogram.size() <= size) {
            histogram.resize(size+1, 0);
        }
        ++histogram[size];
    }
    ofstream csv3("Debug-DisjointSets-Histogram-" + to_string(stage) + ".csv");
    csv3 << "Size,Frequency\n";
    for(uint64_t size=0; size<histogram.size(); size++) {
        csv3 << size << "," << histogram[size] << "\n";
    }

}






void Assembler::accessMarkerGraphVertices(bool readWriteAccess)
{
    markerGraph.vertexTable.accessExisting(
        largeDataName("MarkerGraphVertexTable"), readWriteAccess);

    markerGraph.constructVertices();
    markerGraph.vertices().accessExisting(
        largeDataName("MarkerGraphVertices"), readWriteAccess);
}



void Assembler::checkMarkerGraphVerticesAreAvailable() const
{
    if(!markerGraph.vertices().isOpen() || !markerGraph.vertexTable.isOpen) {
        throw runtime_error("Vertices of the marker graph are not accessible.");
    }
}



MarkerGraph::VertexId Assembler::getGlobalMarkerGraphVertex(
    OrientedReadId orientedReadId,
    uint32_t ordinal) const
{
    const MarkerId markerId =  getMarkerId(orientedReadId, ordinal);
    return markerGraph.vertexTable[markerId];
}



// Get pairs (ordinal, marker graph vertex id) for all markers of an oriented read.
// The pairs are returned sorted by ordinal.
void Assembler::getMarkerGraphVertices(
    OrientedReadId orientedReadId,
    vector< pair<uint32_t, MarkerGraph::VertexId> >& v)
{
    const uint32_t markerCount = uint32_t(markers->size(orientedReadId.getValue()));
    v.clear();
    for(uint32_t ordinal=0; ordinal<markerCount; ordinal++) {
        const MarkerGraph::VertexId vertexId =
            getGlobalMarkerGraphVertex(orientedReadId, ordinal);
        if(vertexId != MarkerGraph::invalidCompressedVertexId) {
            v.push_back(make_pair(ordinal, vertexId));
        }
    }
}






// Find the reverse complement of each marker graph vertex.
void Assembler::findMarkerGraphReverseComplementVertices(uint64_t threadCount)
{
    performanceLog << timestamp << "Begin findMarkerGraphReverseComplementVertices."
        << endl;

    // Check that we have what we need.
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Get the number of vertices in the marker graph.
    using VertexId = MarkerGraph::VertexId;
    const VertexId vertexCount = markerGraph.vertexCount();

    // Allocate the vector to hold the reverse complemented
    // vertex id for each vertex.
    if(not markerGraph.reverseComplementVertex.isOpen) {
        markerGraph.reverseComplementVertex.createNew(
            largeDataName("MarkerGraphReverseComplementVertex"),
            largeDataPageSize);
    }
    markerGraph.reverseComplementVertex.resize(vertexCount);

    // Check each vertex.
    setupLoadBalancing(vertexCount, 10000);
    runThreads(&Assembler::findMarkerGraphReverseComplementVerticesThreadFunction1,
        threadCount);

    // Check that the reverse complement of the reverse complement of a
    // vertex is the vertex itself.
    setupLoadBalancing(vertexCount, 10000);
    runThreads(&Assembler::findMarkerGraphReverseComplementVerticesThreadFunction2,
        threadCount);
    performanceLog << timestamp << "End findMarkerGraphReverseComplementVertices." << endl;

}



void Assembler::findMarkerGraphReverseComplementVerticesThreadFunction1(uint64_t)
{
    using VertexId = MarkerGraph::VertexId;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices in this batch.
        for (VertexId vertexId=begin; vertexId!=end; vertexId++) {

            // Get the markers of this vertex.
            const span<MarkerId> vertexMarkers = markerGraph.getVertexMarkerIds(vertexId);
            DINARA_ASSERT(vertexMarkers.size() > 0);

            // Get the first marker of this vertex.
            const MarkerId firstMarkerId = vertexMarkers[0];

            /// Find the reverse complemented marker.
            const MarkerId firstMarkerIdReverseComplement = findReverseComplement(firstMarkerId);

            // Find the corresponding vertex.
            const VertexId vertexIdReverseComplement = markerGraph.vertexTable[firstMarkerIdReverseComplement];
            DINARA_ASSERT(vertexIdReverseComplement != MarkerGraph::invalidCompressedVertexId);

            // Now check that we get the same reverse complement vertex for all markers.
            for(const MarkerId markerId: vertexMarkers) {
                const MarkerId markerIdReverseComplement = findReverseComplement(markerId);
                DINARA_ASSERT(markerGraph.vertexTable[markerIdReverseComplement] == vertexIdReverseComplement);
            }

            // Store the reverse complement vertex.
            markerGraph.reverseComplementVertex[vertexId] = vertexIdReverseComplement;
        }
    }
}



void Assembler::findMarkerGraphReverseComplementVerticesThreadFunction2(uint64_t)
{
    using VertexId = MarkerGraph::VertexId;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        for (VertexId vertexId=begin; vertexId!=end; vertexId++) {
            const VertexId vertexIdReverseComplement =
                markerGraph.reverseComplementVertex[vertexId];
            DINARA_ASSERT(
                markerGraph.reverseComplementVertex[vertexIdReverseComplement] == vertexId);
        }
    }
}



void Assembler::accessMarkerGraphReverseComplementVertex(bool readWriteAccess)
{
    markerGraph.reverseComplementVertex.accessExisting(
        largeDataName("MarkerGraphReverseComplementeVertex"),
        readWriteAccess);
}




void Assembler::accessMarkerGraphReverseComplementEdge()
{
    markerGraph.reverseComplementEdge.accessExistingReadOnly(
        largeDataName("MarkerGraphReverseComplementeEdge"));
}





void Assembler::accessMarkerGraphEdges(
    bool accessEdgesReadWrite,
    bool accessConnectivityReadWrite)
{
    if(accessEdgesReadWrite) {
        markerGraph.edgeMarkerIntervals.accessExistingReadWrite(
            largeDataName("GlobalMarkerGraphEdgeMarkerIntervals"));
        markerGraph.edges.accessExistingReadWrite(
            largeDataName("GlobalMarkerGraphEdges"));
    } else {
        markerGraph.edgeMarkerIntervals.accessExistingReadOnly(
            largeDataName("GlobalMarkerGraphEdgeMarkerIntervals"));
        markerGraph.edges.accessExistingReadOnly(
            largeDataName("GlobalMarkerGraphEdges"));
    }

    if(accessConnectivityReadWrite) {
        markerGraph.edgesBySource.accessExistingReadWrite(
            largeDataName("GlobalMarkerGraphEdgesBySource"));
        markerGraph.edgesByTarget.accessExistingReadWrite(
            largeDataName("GlobalMarkerGraphEdgesByTarget"));
    } else {
        markerGraph.edgesBySource.accessExistingReadOnly(
            largeDataName("GlobalMarkerGraphEdgesBySource"));
        markerGraph.edgesByTarget.accessExistingReadOnly(
            largeDataName("GlobalMarkerGraphEdgesByTarget"));
    }
}





void Assembler::accessMarkerGraphConsensus()
{
    if(assemblerInfo->assemblyMode == 3) {
        markerGraph.edgeSequence.accessExistingReadOnly(largeDataName("MarkerGraphEdgesSequence"));

    } else {
        if(assemblerInfo->readRepresentation == 1) {
            markerGraph.vertexRepeatCounts.accessExistingReadOnly(
                largeDataName("MarkerGraphVertexRepeatCounts"));
        }
        markerGraph.edgeConsensus.accessExistingReadOnly(
            largeDataName("MarkerGraphEdgesConsensus"));
        markerGraph.edgeConsensusOverlappingBaseCount.accessExistingReadOnly(
            largeDataName("MarkerGraphEdgesConsensusOverlappingBaseCount"));
    }

}






// Find the common KmerId for all the markers of a marker graph vertex.
KmerId Assembler::getMarkerGraphVertexKmerId(MarkerGraphVertexId vertexId) const
{
    return markerGraph.getVertexKmerId(
        vertexId,
        assemblerInfo->k,
        *reads,
        *markers);
}
