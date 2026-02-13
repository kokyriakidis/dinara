// ============================================================================
// BRG-aware marker graph vertex creation.
//
// Identical to createMarkerGraphVertices() except for the disjoint-set union
// step: instead of iterating ReadGraph edges (strand-doubled), we iterate
// BidirectionalReadGraph edges (one per alignment, undirected, physical-read
// indexed).  Only non-deleted BRG edges are used, so the marker graph
// reflects the cleaned overlap set produced by cleanBidirectionalReadGraph*.
//
// The output is a standard MarkerGraph (vertexTable, vertices, disjoint-set
// histogram).  All downstream consumers (filterMarkerGraphVerticesByRepeatKmers,
// findMarkerGraphReverseComplementVertices, createMarkerGraphEdges, anchor
// creation, etc.) work unchanged.
// ============================================================================

// Dinara.
#include "Assembler.hpp"
#include "chrono.hpp"
#include "compressAlignment.hpp"
#include "dset64-gccAtomic.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"

using namespace dinara;
using namespace std;



void Assembler::createMarkerGraphVerticesFromBrg(

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
    const auto tBegin = steady_clock::now();
    performanceLog << timestamp << "Begin computing marker graph vertices from BRG." << endl;

    // Check that we have what we need.
    reads->checkReadsAreOpen();
    reads->checkReadFlagsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    checkBidirectionalReadGraphIsOpen();
    DINARA_ASSERT(compressedAlignments.isOpen());

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store parameters so they are accessible to the threads.
    auto& data = createMarkerGraphVerticesData;
    data.allowDuplicateMarkers = allowDuplicateMarkers;
    data.minCoveragePerStrand = minCoveragePerStrand;

    // Initialize computation of the global marker graph.
    data.orientedMarkerCount = markers->totalSize();

    data.disjointSetTable.createNew(
        largeDataName("tmp-DisjointSetTable"),
        largeDataPageSize);
    // DisjointSets needs 128 bits per entry for lock-free union-find.
    // We allocate 2× so the memory can be used as 128-bit integers.
    data.disjointSetTable.reserveAndResize(data.orientedMarkerCount * 2);

    // Have DisjointSets use the memory allocated in data.disjointSetTable.
    data.disjointSetsPointer = std::make_shared<DisjointSets>(
        reinterpret_cast<DisjointSets::Aint*>(data.disjointSetTable.begin()),
        data.orientedMarkerCount
    );


    // =====================================================================
    // STEP 1: Disjoint-set union from BRG edges.
    //
    // This is the BRG-specific part.  Instead of iterating ReadGraph edges
    // (strand-doubled pairs), we iterate BRG edges directly.  Each BRG edge
    // corresponds to one alignment.  We derive OrientedReadIds via
    // edge.traverse(), decompress the alignment, and unite marker pairs.
    // =====================================================================
    performanceLog << timestamp << "BRG disjoint set computation begins." << endl;
    cout << timestamp << "BRG: building marker graph vertices from "
         << bidirectionalReadGraph.edges.size() << " BRG edges." << endl;

    {
        const uint64_t edgeCount = bidirectionalReadGraph.edges.size();
        const uint64_t chunk = std::max<uint64_t>(1, edgeCount / threadCount + 1);

        vector<thread> threads;
        threads.reserve(threadCount);

        for(uint64_t t = 0; t < threadCount; ++t) {
            threads.emplace_back([&, t]() {
                const uint64_t begin = t * chunk;
                const uint64_t end = std::min(edgeCount, (t + 1) * chunk);

                Alignment alignment;
                DisjointSets& disjointSets = *data.disjointSetsPointer;

                for(uint64_t edgeId = begin; edgeId < end; ++edgeId) {
                    const auto& edge = bidirectionalReadGraph.edges[edgeId];

                    // Skip deleted edges (BRG cleaning removes these).
                    if(edge.isDeleted) {
                        continue;
                    }

                    // Skip edges with inconsistent alignments.
                    if(edge.hasInconsistentAlignment) {
                        continue;
                    }

                    const uint64_t alignmentId = edge.alignmentId;

                    // Skip if not in read graph (safety check).
                    if(!alignmentData[alignmentId].info.isInReadGraph) {
                        continue;
                    }

                    // If either read is chimeric, skip.
                    if(reads->getFlags(edge.readIds[0]).isChimeric ||
                       reads->getFlags(edge.readIds[1]).isChimeric) {
                        continue;
                    }

                    // Derive OrientedReadIds from the BRG edge.
                    // Convention: readIds[0] enters at strand 0.
                    const ReadId readId0 = edge.readIds[0];
                    const OrientedReadId orientedReadId0(readId0, 0);
                    const auto [readId1, strand1] = edge.traverse(readId0, Strand(0));
                    const OrientedReadId orientedReadId1(readId1, strand1);

                    // Decompress the alignment.
                    span<const char> compressedAlignment = compressedAlignments[alignmentId];
                    dinara::decompress(compressedAlignment, alignment);

                    // The compressed alignment stores ordinals for the canonical
                    // (readIds[0] strand-0, readIds[1] strand-as-stored) pair.
                    // Since we constructed orientedReadId0 as (readIds[0], 0)
                    // and orientedReadId1 via traverse(), the ordinals are correct
                    // for the (orientedReadId0, orientedReadId1) pair.

                    // Unite aligned marker pairs, their reverse complements,
                    // and merge each marker with its own RC.
                    // This collapses forward and RC kmers into one vertex,
                    // so that a single marker graph vertex represents the
                    // same genomic position regardless of strand.
                    for(const auto& p : alignment.ordinals) {
                        const uint32_t ordinal0 = p[0];
                        const uint32_t ordinal1 = p[1];
                        const MarkerId markerId0 = getMarkerId(orientedReadId0, ordinal0);
                        const MarkerId markerId1 = getMarkerId(orientedReadId1, ordinal1);
                        const MarkerId markerId0Rc = findReverseComplement(markerId0);
                        const MarkerId markerId1Rc = findReverseComplement(markerId1);

                        // Unite the aligned marker pair.
                        disjointSets.unite(markerId0, markerId1);

                        // Merge each marker with its reverse complement.
                        // This ensures that a kmer and its RC kmer share
                        // the same vertex — each vertex is self-RC.
                        disjointSets.unite(markerId0, markerId0Rc);
                        disjointSets.unite(markerId1, markerId1Rc);
                    }
                }
            });
        }

        for(auto& th : threads) {
            th.join();
        }
    }

    performanceLog << timestamp << "BRG disjoint set computation completed." << endl;


    // =====================================================================
    // STEPS 2+: Find representatives, compact, renumber, gather, filter.
    // This is identical to the code in createMarkerGraphVertices after
    // the disjoint-set union step.
    // =====================================================================

    // Find the disjoint set that each oriented marker was assigned to.
    const size_t batchSize = 10000;
    uint64_t pass = 1;
    do {
        (data.disjointSetsPointer)->parentUpdated = 0;
        performanceLog << "    " << timestamp << " Iteration  " << pass << endl;
        setupLoadBalancing(data.orientedMarkerCount, batchSize);
        runThreads(&Assembler::createMarkerGraphVerticesThreadFunction2, threadCount);
        performanceLog << "    " << timestamp << " Updated parent of - "
                       << (data.disjointSetsPointer)->parentUpdated << " entries." << endl;
        pass++;
    } while ((data.disjointSetsPointer)->parentUpdated > 0 && pass <= 10);

    if (pass > 10) {
        string errorMsg = "DisjointSets parent information did not converge in "
                          + to_string(pass) + " iterations.";
        throw runtime_error(errorMsg);
    }

    performanceLog << timestamp << "Verifying convergence of parent information." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction21, threadCount);
    performanceLog << timestamp << "Done verifying convergence of parent information." << endl;

    // Compact: parent at 2*N → move to N.
    performanceLog << timestamp << "Compacting the Disjoint Set data-structure." << endl;
    for(uint64_t i = 0; i < data.orientedMarkerCount; i++) {
        data.disjointSetTable[i] = data.disjointSetTable[2 * i];
    }
    data.disjointSetTable.resize(data.orientedMarkerCount);
    data.disjointSetTable.unreserve();
    performanceLog << timestamp << "Done compacting the Disjoint Set data-structure." << endl;

    // Don't need DisjointSets any more.
    data.disjointSetsPointer = 0;


    // Count the number of markers in each disjoint set.
    performanceLog << timestamp << "Counting the number of markers in each disjoint set." << endl;
    data.workArea.createNew(
        largeDataName("tmp-WorkArea"),
        largeDataPageSize);
    data.workArea.reserveAndResize(data.orientedMarkerCount);
    fill(data.workArea.begin(), data.workArea.end(), 0ULL);
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction3, threadCount);


    // Compute histogram and store it.
    {
        vector<uint64_t> histogram;
        for(MarkerGraph::VertexId i = 0; i < data.orientedMarkerCount; i++) {
            const MarkerGraph::VertexId markerCount = data.workArea[i];
            if(markerCount == 0) {
                continue;
            }
            if(markerCount >= histogram.size()) {
                histogram.resize(markerCount + 1, 0);
            }
            ++histogram[markerCount];
        }

        // Store the disjoint sets histogram.
        {
            markerGraph.disjointSetsHistogram.createNew(
                largeDataName("DisjointSetsHistogram"),
                largeDataPageSize);
            for(uint64_t coverage = 0; coverage < histogram.size(); coverage++) {
                const uint64_t frequency = histogram[coverage];
                if(frequency) {
                    markerGraph.disjointSetsHistogram.push_back({coverage, frequency});
                }
            }
        }

        ofstream csv("DisjointSetsHistogram.csv");
        csv << "Coverage,Frequency\n";
        for(uint64_t coverage = 0; coverage < histogram.size(); coverage++) {
            const uint64_t frequency = histogram[coverage];
            if(frequency) {
                csv << coverage << "," << frequency << "\n";
            }
        }

        if(minCoverage == 0) {
            // Automatic selection not implemented for BRG path;
            // require explicit minCoverage.
            throw runtime_error(
                "createMarkerGraphVerticesFromBrg requires explicit minCoverage > 0.");
        }
    }
    // Store the value of minCoverage actually used.
    assemblerInfo->markerGraphMinCoverageUsed = minCoverage;


    // Renumber disjoint sets keeping only those with coverage in [minCoverage, maxCoverage].
    performanceLog << timestamp << "Renumbering the disjoint sets." << endl;
    MarkerGraph::VertexId newDisjointSetId = 0ULL;
    for(MarkerGraph::VertexId oldDisjointSetId = 0;
        oldDisjointSetId < data.orientedMarkerCount; ++oldDisjointSetId) {
        auto& w = data.workArea[oldDisjointSetId];
        const MarkerGraph::VertexId mc = w;
        if(mc < minCoverage || mc > maxCoverage) {
            w = MarkerGraph::invalidVertexId;
        } else {
            w = newDisjointSetId;
            ++newDisjointSetId;
        }
    }
    const auto disjointSetCount = newDisjointSetId;
    cout << "BRG: Kept " << disjointSetCount
         << " disjoint sets with coverage in the requested range." << endl;


    // Reassign vertices to disjoint sets using new numbering.
    performanceLog << timestamp << "Assigning vertices to renumbered disjoint sets." << endl;
    for(MarkerGraph::VertexId markerId = 0;
        markerId < data.orientedMarkerCount; ++markerId) {
        auto& d = data.disjointSetTable[markerId];
        const auto oldId = d;
        const auto newId = data.workArea[oldId];
        d = newId;
    }
    data.workArea.remove();


    // Gather the markers in each disjoint set.
    data.disjointSetMarkers.createNew(
        largeDataName("tmp-DisjointSetMarkers"),
        largeDataPageSize);
    performanceLog << timestamp << "Gathering markers in disjoint sets, pass1." << endl;
    data.disjointSetMarkers.beginPass1(disjointSetCount);
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction4, threadCount);
    performanceLog << timestamp << "Gathering markers in disjoint sets, pass2." << endl;
    data.disjointSetMarkers.beginPass2();
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction5, threadCount);
    data.disjointSetMarkers.endPass2();


    // Sort the markers in each disjoint set.
    performanceLog << timestamp << "Sorting the markers in each disjoint set." << endl;
    setupLoadBalancing(disjointSetCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction6, threadCount);


    // Flag "bad" disjoint sets.
    data.isBadDisjointSet.createNew(
        largeDataName("tmp-IsBadDisjointSet"),
        largeDataPageSize);
    data.isBadDisjointSet.reserveAndResize(disjointSetCount);
    performanceLog << timestamp << "Flagging bad disjoint sets." << endl;
    setupLoadBalancing(disjointSetCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction7, threadCount);
    const size_t badDisjointSetCount = std::count(
        data.isBadDisjointSet.begin(), data.isBadDisjointSet.end(), true);
    cout << "BRG: Found " << badDisjointSetCount << " bad disjoint sets "
        "(duplicate markers or insufficient per-strand coverage)." << endl;


    // Renumber again, excluding bad sets.
    performanceLog << timestamp << "Renumbering disjoint sets to remove the bad ones." << endl;
    data.workArea.createNew(
        largeDataName("tmp-WorkArea"),
        largeDataPageSize);
    data.workArea.reserveAndResize(disjointSetCount);
    newDisjointSetId = 0ULL;
    for(MarkerGraph::VertexId oldDisjointSetId = 0;
        oldDisjointSetId < disjointSetCount; ++oldDisjointSetId) {
        auto& w = data.workArea[oldDisjointSetId];
        if(data.isBadDisjointSet[oldDisjointSetId]) {
            w = MarkerGraph::invalidVertexId;
        } else {
            w = newDisjointSetId;
            ++newDisjointSetId;
        }
    }
    DINARA_ASSERT(newDisjointSetId + badDisjointSetCount == disjointSetCount);


    // Compute the final vertex id for each marker.
    performanceLog << timestamp << "Assigning vertex ids to markers." << endl;
    markerGraph.vertexTable.createNew(
        largeDataName("MarkerGraphVertexTable"),
        largeDataPageSize);
    markerGraph.vertexTable.reserveAndResize(data.orientedMarkerCount);
    for(MarkerGraph::VertexId markerId = 0;
        markerId < data.orientedMarkerCount; ++markerId) {
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
    performanceLog << timestamp << "Gathering the markers of each vertex of the marker graph." << endl;
    markerGraph.constructVertices();
    markerGraph.vertices().createNew(
        largeDataName("MarkerGraphVertices"),
        largeDataPageSize);
    for(MarkerGraph::VertexId oldDisjointSetId = 0;
        oldDisjointSetId < disjointSetCount; ++oldDisjointSetId) {
        if(data.isBadDisjointSet[oldDisjointSetId]) {
            continue;
        }
        markerGraph.vertices().appendVector();
        const auto dsMarkers = data.disjointSetMarkers[oldDisjointSetId];
        for(const MarkerId markerId : dsMarkers) {
            markerGraph.vertices().append(markerId);
        }
    }
    markerGraph.vertices().unreserve();

    data.isBadDisjointSet.remove();
    data.disjointSetMarkers.remove();


    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "BRG: Marker graph vertex creation complete: "
         << markerGraph.vertexCount() << " vertices in " << tTotal << " s." << endl;
    performanceLog << timestamp << "Computation of marker graph vertices from BRG "
                   << "completed in " << tTotal << " s." << endl;
}
