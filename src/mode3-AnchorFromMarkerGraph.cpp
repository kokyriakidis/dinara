/********************************************************************************

Creation of Anchors from the MarkerGraph.

Each Anchor corresponds to a "primary marker graph edge"
defined as follows:

- All contributing oriented reads have exactly the same sequence.
  If more than one distinct sequence is present, the edge (anchor)
  is split into two edges (anchors).
- Edge coverage is >= minPrimaryCoverage and <= maxPrimaryCoverage.
- Both vertices have no duplicate ReadIds, and as a result the
  resulting anchor has no duplicate ReadIds.

This uses as input the MarkerGraph vertices. MarkerGraph edges
are created implicitly and never stored. Instead, the corresponding
information is stored directly in the Anchors.

This code is equivalent to Assembler::createPrimaryMarkerGraphEdges.

********************************************************************************/

#include "mode3-Anchor.hpp"
#include "findMarkerId.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "MarkerGraphEdgePairInfo.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace dinara;
using namespace mode3;


Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage,
    uint64_t threadCount,
    bool createFromVertices) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    performanceLog << timestamp << "Anchor creation from the marker graph begins." << endl;

    // Sanity checks and get kHalf.
    DINARA_ASSERT((k %2) == 0);
    kHalf = k / 2;
    DINARA_ASSERT(reads.representation == 0);
    DINARA_ASSERT(markers.isOpen());
    DINARA_ASSERT(markerGraph.vertices().isOpen());
    DINARA_ASSERT(markerGraph.vertexTable.isOpen);

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store the arguments so the threads can see them.
    auto& data = constructFromMarkerGraphData;
    data.minPrimaryCoverage = minPrimaryCoverage;
    data.maxPrimaryCoverage = maxPrimaryCoverage;
    data.createFromVertices = createFromVertices;

    data.markerGraphPointer = &markerGraph;

    // Make space for the edges found by each thread.
    constructFromMarkerGraphData.threadMarkerIntervals.resize(threadCount);
    constructFromMarkerGraphData.threadSequences.resize(threadCount);

    // Parallelize over the MarkerGraph source vertex.
    // Each thread stores in a separate vector the Anchors it finds.
    uint64_t batchSize = 1000;
    setupLoadBalancing(markerGraph.vertices().size(), batchSize);
    if(createFromVertices) {
        runThreads(&Anchors::constructFromMarkerGraphVerticesThreadFunction, threadCount);
    } else {
        runThreads(&Anchors::constructFromMarkerGraphThreadFunction, threadCount);
    }

    // Gather the Anchors found by all threads.
    anchorMarkerIntervals.createNew(
            largeDataName("AnchorMarkerIntervals"),
            largeDataPageSize);
    anchorSequences.createNew(
        largeDataName("AnchorSequences"), largeDataPageSize);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadMarkerIntervals = *data.threadMarkerIntervals[threadId];
        auto& threadSequences = *data.threadSequences[threadId];
        DINARA_ASSERT(threadMarkerIntervals.size() == threadSequences.size());
        for(uint64_t i=0; i<threadMarkerIntervals.size(); i++) {
            const auto thisAnchorMarkerIntervals = threadMarkerIntervals[i];
            anchorMarkerIntervals.appendVector();
            for(const auto& threadMarkerInterval: thisAnchorMarkerIntervals) {
                anchorMarkerIntervals.append(
                    AnchorMarkerInterval(threadMarkerInterval.orientedReadId, threadMarkerInterval.ordinal0));
            }
            const span<Base> sequence = threadSequences[i];
            anchorSequences.appendVector(sequence.begin(), sequence.end());
        }
        threadMarkerIntervals.remove();
        threadSequences.remove();
        data.threadMarkerIntervals[threadId] = 0;
        data.threadSequences[threadId] = 0;
    }
    data.threadMarkerIntervals.clear();
    data.threadSequences.clear();
    cout << "Found " << anchorMarkerIntervals.size() << " anchors." << endl;

    // Initialize the AnchorInfos with the ordinal offset, which is always
    // 1 when creating anchors from the marker graph.
    // The componentId and localAnchorIdInComponent fields will be filled later,
    // when processing each component.
    anchorInfos.createNew(largeDataName("AnchorInfos"), largeDataPageSize);
    anchorInfos.resize(anchorMarkerIntervals.size());
    const uint32_t ordinalOffset = createFromVertices ? 0 : 1;
    for(AnchorInfo& anchorInfo: anchorInfos) {
        anchorInfo.ordinalOffset = ordinalOffset;
        anchorInfo.componentId = invalid<uint32_t>;
        anchorInfo.localAnchorIdInComponent = invalid<uint64_t>;
    }

    performanceLog << timestamp << "Anchor creation from the marker graph ends." << endl;
}



// Create Anchors from a selected subset of marker graph vertices.
Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    const vector<MarkerGraphVertexId>& selectedVertexIds,
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage,
    uint64_t threadCount) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    performanceLog << timestamp << "Anchor creation from selected marker graph vertices begins." << endl;

    // Sanity checks and get kHalf.
    DINARA_ASSERT((k % 2) == 0);
    kHalf = k / 2;
    DINARA_ASSERT(reads.representation == 0);
    DINARA_ASSERT(markers.isOpen());
    DINARA_ASSERT(markerGraph.vertices().isOpen());
    DINARA_ASSERT(markerGraph.vertexTable.isOpen);
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store the arguments so the threads can see them.
    auto& data = constructFromMarkerGraphData;
    data.minPrimaryCoverage = minPrimaryCoverage;
    data.maxPrimaryCoverage = maxPrimaryCoverage;
    data.createFromVertices = true;
    data.markerGraphPointer = &markerGraph;
    data.selectedVertexIdsPointer = &selectedVertexIds;

    // Make space for the edges found by each thread.
    data.threadMarkerIntervals.resize(threadCount);
    data.threadSequences.resize(threadCount);

    // Parallelize over the selected vertices.
    const uint64_t batchSize = 1000;
    setupLoadBalancing(selectedVertexIds.size(), batchSize);
    runThreads(&Anchors::constructFromSelectedMarkerGraphVerticesThreadFunction, threadCount);

    // Gather the Anchors found by all threads.
    anchorMarkerIntervals.createNew(
        largeDataName("AnchorMarkerIntervals"),
        largeDataPageSize);
    anchorSequences.createNew(
        largeDataName("AnchorSequences"), largeDataPageSize);
    for(uint64_t threadId=0; threadId<threadCount; threadId++) {
        auto& threadMarkerIntervals = *data.threadMarkerIntervals[threadId];
        auto& threadSequences = *data.threadSequences[threadId];
        DINARA_ASSERT(threadMarkerIntervals.size() == threadSequences.size());
        for(uint64_t i=0; i<threadMarkerIntervals.size(); i++) {
            const auto thisAnchorMarkerIntervals = threadMarkerIntervals[i];
            anchorMarkerIntervals.appendVector();
            for(const auto& threadMarkerInterval: thisAnchorMarkerIntervals) {
                anchorMarkerIntervals.append(
                    AnchorMarkerInterval(threadMarkerInterval.orientedReadId, threadMarkerInterval.ordinal0));
            }
            const span<Base> sequence = threadSequences[i];
            anchorSequences.appendVector(sequence.begin(), sequence.end());
        }
        threadMarkerIntervals.remove();
        threadSequences.remove();
        data.threadMarkerIntervals[threadId] = 0;
        data.threadSequences[threadId] = 0;
    }
    data.threadMarkerIntervals.clear();
    data.threadSequences.clear();
    data.selectedVertexIdsPointer = nullptr;

    cout << "Found " << anchorMarkerIntervals.size() << " anchors." << endl;

    // Initialize the AnchorInfos with ordinalOffset=0 for vertex anchors.
    anchorInfos.createNew(largeDataName("AnchorInfos"), largeDataPageSize);
    anchorInfos.resize(anchorMarkerIntervals.size());
    for(AnchorId anchorId=0; anchorId<anchorMarkerIntervals.size(); anchorId++) {
        anchorInfos[anchorId].ordinalOffset = 0;
        anchorInfos[anchorId].componentId = invalid<uint32_t>;
        anchorInfos[anchorId].localAnchorIdInComponent = invalid<uint64_t>;
    }

    performanceLog << timestamp << "Anchor creation from selected marker graph vertices ends." << endl;
}



void Anchors::constructFromMarkerGraphThreadFunction(uint64_t threadId)
{
    // Access the data set up by createPrimaryMarkerGraphEdges.
    using ThreadMarkerInterval = ConstructFromMarkerGraphData::ThreadMarkerInterval;
    auto& data = constructFromMarkerGraphData;

    // Get the primary coverage range.
    const uint64_t minPrimaryCoverage = data.minPrimaryCoverage;
    const uint64_t maxPrimaryCoverage = data.maxPrimaryCoverage;

    const MarkerGraph& markerGraph = *data.markerGraphPointer;
    const MemoryMapped::VectorOfVectors<MarkerId, MarkerGraph::CompressedVertexId>& markerGraphVertices =
        markerGraph.vertices();

    // Create the vector to contain the marker intervals for the Anchors found by this thread.
    shared_ptr< MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t> > markerIntervalsPointer =
        make_shared< MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t> >();
    data.threadMarkerIntervals[threadId] = markerIntervalsPointer;
    MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t>&
        markerIntervals = *markerIntervalsPointer;
    markerIntervals.createNew(
            largeDataName("tmp-ThreadAnchorMarkerIntervals-" + to_string(threadId)),
            largeDataPageSize);

    // Create the vector to contain the sequences for the Anchors found by this thread.
    shared_ptr< MemoryMapped::VectorOfVectors<Base, uint64_t> > sequencesPointer =
        make_shared< MemoryMapped::VectorOfVectors<Base, uint64_t> >();
    data.threadSequences[threadId] = sequencesPointer;
    MemoryMapped::VectorOfVectors<Base, uint64_t>& sequences = *sequencesPointer;
    sequences.createNew(
            largeDataName("tmp-ThreadAnchorSequences-" + to_string(threadId)),
            largeDataPageSize);

    // For each sequence and next vertex we encounter, we store the corresponding MarkerIntervals.
    class AnchorCandidate {
    public:
        vector<Base> sequence;
        MarkerGraphVertexId vertexId1;
        vector<ThreadMarkerInterval> markerIntervals;
        AnchorCandidate(
            const vector<Base>& sequence,
            MarkerGraphVertexId vertexId1,
            const ThreadMarkerInterval& markerInterval) :
            sequence(sequence),
            vertexId1(vertexId1),
            markerIntervals(1, markerInterval) {}
    };
    vector<AnchorCandidate> anchorCandidates;

    // Work vector used in the loop below.
    vector<Base> sequence;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over source vertex ids assigned to this batch.
        for(MarkerGraphVertexId vertexId0=begin; vertexId0!=end; vertexId0++) {

            anchorCandidates.clear();

            // Access the markers of this vertex.
            const auto vertexMarkerIds = markerGraphVertices[vertexId0];

            // If vertex coverage of this vertex is less than minPrimaryCoverage,
            // no primary edges start here for sure.
            if(vertexMarkerIds.size() < minPrimaryCoverage) {
                continue;
            }

            // If this vertex has duplicate ReadIds, no primary edge can begin here.
            if(markerGraph.vertexHasDuplicateReadIds(vertexId0, markers)) {
                continue;
            }

            // Loop over the markers of this vertex.
            for(const MarkerId markerId0: vertexMarkerIds) {
                OrientedReadId orientedReadId;
                uint32_t ordinal0;
                tie(orientedReadId, ordinal0) = dinara::findMarkerId(markerId0, markers);

                // Find the next marker graph vertex visited by this OrientedReadId.
                const uint32_t ordinal1 = ordinal0 + 1;
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];
                if(ordinal1 == orientedReadMarkers.size()) {
                    // There is no next vertex.
                    continue;
                }
                DINARA_ASSERT(ordinal1 < orientedReadMarkers.size());
                const MarkerId markerId1 = markerId0 + 1;
                const MarkerGraphVertexId vertexId1 = markerGraph.vertexTable[markerId1];

                // If vertexId1 has coverage less than minPrimaryCoverage, there
                // cannot be a primary edge vertexId0->vertexId1.
                if(markerGraphVertices.size(vertexId1) < minPrimaryCoverage) {
                    continue;
                }

                // Get the sequence between these two markers.
                const uint32_t position0 = uint32_t(markers.begin()[markerId0].position + kHalf);
                const uint32_t position1 = uint32_t(markers.begin()[markerId1].position + kHalf);
                sequence.clear();
                for(uint32_t position=position0; position<position1; position++) {
                    const Base base = reads.getOrientedReadBase(orientedReadId, position);
                    sequence.push_back(base);
                }

                // Construct this ThreadMarkerInterval.
                ThreadMarkerInterval markerInterval;
                markerInterval.orientedReadId = orientedReadId;
                markerInterval.ordinal0 = ordinal0;

                // If we already have an AnchorCandidate with this sequence and vertexId1, add id there.
                // Otherwise create a new Info.
                bool found = false;
                for(AnchorCandidate& anchorCandidate: anchorCandidates) {
                    if(anchorCandidate.sequence == sequence and anchorCandidate.vertexId1 == vertexId1) {
                        anchorCandidate.markerIntervals.push_back(markerInterval);
                        found = true;
                        break;
                    }
                }
                if(not found) {
                    anchorCandidates.push_back(AnchorCandidate(sequence, vertexId1, markerInterval));
                }
            }

            // Each of our Infos object can generate an Anchor, but we have to check that:
            // - Coverage is in the primary coverage range.
            // - vertexId1 does not have duplicate ReadIds.
            // In addition, we want to make sure that pairs of reverse complemented Anchors
            // are numbered consecutively.
            // To achieve this, we do the following:
            // - If the first OrientedReadId is on strand 1, we do nothing.
            // - Otherwise, we store this Anchor and its reverse complement.
            for(AnchorCandidate& anchorCandidate: anchorCandidates) {
                if(anchorCandidate.markerIntervals.size() < minPrimaryCoverage) {
                    continue;
                }
                if(anchorCandidate.markerIntervals.size() > maxPrimaryCoverage) {
                    continue;
                }
                if(markerGraph.vertexHasDuplicateReadIds(anchorCandidate.vertexId1, markers)) {
                    continue;
                }

                // If the first OrientedReadId is on strand 1, we do nothing.
                if(anchorCandidate.markerIntervals.front().orientedReadId.getStrand() == 1) {
                    continue;
                }

                // Store this Anchor.
                markerIntervals.appendVector(anchorCandidate.markerIntervals);
                sequences.appendVector(anchorCandidate.sequence);

                // Reverse complement the marker intervals.
                for(auto& markerInterval: anchorCandidate.markerIntervals) {
                    markerInterval.orientedReadId.flipStrand();
                    const uint64_t markerCount = markers.size(markerInterval.orientedReadId.getValue());
                    markerInterval.ordinal0 = uint32_t(markerCount) - 2 - markerInterval.ordinal0;
                }

                // Reverse complement the sequence.
                std::reverse(anchorCandidate.sequence.begin(), anchorCandidate.sequence.end());
                for(Base& b: anchorCandidate.sequence) {
                    b.complementInPlace();
                }

                // Store the reverse complement Anchor.
                markerIntervals.appendVector(anchorCandidate.markerIntervals);
                sequences.appendVector(anchorCandidate.sequence);
            }

        }
    }
}



void Anchors::constructFromMarkerGraphVerticesThreadFunction(uint64_t threadId)
{
    // Access the data set up by createPrimaryMarkerGraphEdges.
    using ThreadMarkerInterval = ConstructFromMarkerGraphData::ThreadMarkerInterval;
    auto& data = constructFromMarkerGraphData;

    // Get the primary coverage range.
    const uint64_t minPrimaryCoverage = data.minPrimaryCoverage;
    const uint64_t maxPrimaryCoverage = data.maxPrimaryCoverage;

    const MarkerGraph& markerGraph = *data.markerGraphPointer;
    const MemoryMapped::VectorOfVectors<MarkerId, MarkerGraph::CompressedVertexId>& markerGraphVertices =
        markerGraph.vertices();

    // Create the vector to contain the marker intervals for the Anchors found by this thread.
    shared_ptr< MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t> > markerIntervalsPointer =
        make_shared< MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t> >();
    data.threadMarkerIntervals[threadId] = markerIntervalsPointer;
    MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t>&
        markerIntervals = *markerIntervalsPointer;
    markerIntervals.createNew(
            largeDataName("tmp-ThreadAnchorMarkerIntervals-" + to_string(threadId)),
            largeDataPageSize);

    // Create the vector to contain the sequences for the Anchors found by this thread.
    // For vertex-based anchors, sequences are empty.
    shared_ptr< MemoryMapped::VectorOfVectors<Base, uint64_t> > sequencesPointer =
        make_shared< MemoryMapped::VectorOfVectors<Base, uint64_t> >();
    data.threadSequences[threadId] = sequencesPointer;
    MemoryMapped::VectorOfVectors<Base, uint64_t>& sequences = *sequencesPointer;
    sequences.createNew(
            largeDataName("tmp-ThreadAnchorSequences-" + to_string(threadId)),
            largeDataPageSize);

    // Temporary storage for current anchor markers.
    vector<ThreadMarkerInterval> currentAnchorMarkers;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices assigned to this batch.
        for(MarkerGraphVertexId vertexId=begin; vertexId!=end; vertexId++) {

            // To ensure canonical ordering and consecutive numbering of RC pairs,
            // we process a vertex pair only when we encounter the "canonical" one
            // (e.g. vertexId <= rcVertexId).
            // We then generate both the Forward and Reverse Complement anchors immediately.
            MarkerGraphVertexId rcVertexId = markerGraph.reverseComplementVertex[vertexId];
            if(vertexId > rcVertexId) {
                continue;
            }

            // Access the markers of this vertex.
            const auto vertexMarkerIds = markerGraphVertices[vertexId];

            // Check coverage.
            if(vertexMarkerIds.size() < minPrimaryCoverage || vertexMarkerIds.size() > maxPrimaryCoverage) {
                continue;
            }

            // Check for duplicate ReadIds.
            if(markerGraph.vertexHasDuplicateReadIds(vertexId, markers)) {
                continue;
            }

            // Gather markers for the Forward anchor.
            currentAnchorMarkers.clear();
            currentAnchorMarkers.reserve(vertexMarkerIds.size());

            // If the vertex is on strand 1 (canonical check heuristic),
            // the existing code skipped it.
            // But here we rely on vertexId <= rcVertexId.
            // However, we should check if the markers inside are strand 0 or 1?
            // MarkerGraph vertices usually separate strands (unless palindrome).
            // Let's gather them.

            for(const MarkerId markerId: vertexMarkerIds) {
                 OrientedReadId orientedReadId;
                 uint32_t ordinal0;
                 tie(orientedReadId, ordinal0) = dinara::findMarkerId(markerId, markers);

                 ThreadMarkerInterval interval;
                 interval.orientedReadId = orientedReadId;
                 interval.ordinal0 = ordinal0;
                 currentAnchorMarkers.push_back(interval);
            }

            // Sort and deduplicate by oriented read.
            // Many downstream routines (for example analyzeAnchorPair/countCommon) assume
            // anchor marker intervals are sorted by OrientedReadId and contain no duplicates.
            std::sort(currentAnchorMarkers.begin(), currentAnchorMarkers.end(),
                [](const ThreadMarkerInterval& a, const ThreadMarkerInterval& b) {
                    if(a.orientedReadId != b.orientedReadId) {
                        return a.orientedReadId < b.orientedReadId;
                    }
                    return a.ordinal0 < b.ordinal0;
                });
            currentAnchorMarkers.erase(
                std::unique(currentAnchorMarkers.begin(), currentAnchorMarkers.end(),
                    [](const ThreadMarkerInterval& a, const ThreadMarkerInterval& b) {
                        return a.orientedReadId == b.orientedReadId;
                    }),
                currentAnchorMarkers.end());

            // Re-check coverage after deduplication.
            if(currentAnchorMarkers.size() < minPrimaryCoverage || currentAnchorMarkers.size() > maxPrimaryCoverage) {
                continue;
            }

            // Check if the "representative" marker is Strand 1?
            // Existing edge-based logic: "if (front.strand == 1) continue".
            // If vertexId < rcVertexId, usually this corresponds to Strand 0.
            // But checking explicitly matches old logic style.
            if(!currentAnchorMarkers.empty()) {
                if(currentAnchorMarkers.front().orientedReadId.getStrand() == 1) {
                    // This creates a conflict if vertexId < rcVertexId but it contains Strand 1.
                    // However, markerGraph logic usually assigns LOWER vertexId to Strand 0.
                    // If we skip here, we might miss it if RC vertex (Strand 0) is > vertexId.
                    // Unlikely.
                    // We will trust vertexId <= rcVertexId check as primary.
                    // And duplicate checks ensuring we don't output twice.
                    // But if it's a palindrome (vertexId == rcVertexId), we process.
                }
            }

            // Store Forward Anchor.
            markerIntervals.appendVector(currentAnchorMarkers);
            sequences.appendVector(vector<Base>()); // Empty sequence

            // Store Reverse Complement Anchor.
            for(auto& interval : currentAnchorMarkers) {
                // Flip strand.
                interval.orientedReadId.flipStrand();
                // Re-calculate ordinal.
                // For a single marker M at ordinal k in Read with N markers:
                // RC(M) is at ordinal N - 1 - k.
                // We need N.
                const uint64_t markerCount = markers.size(interval.orientedReadId.getValue());
                interval.ordinal0 = uint32_t(markerCount) - 1 - interval.ordinal0;
            }
            std::sort(currentAnchorMarkers.begin(), currentAnchorMarkers.end(),
                [](const ThreadMarkerInterval& a, const ThreadMarkerInterval& b) {
                    if(a.orientedReadId != b.orientedReadId) {
                        return a.orientedReadId < b.orientedReadId;
                    }
                    return a.ordinal0 < b.ordinal0;
                });
            currentAnchorMarkers.erase(
                std::unique(currentAnchorMarkers.begin(), currentAnchorMarkers.end(),
                    [](const ThreadMarkerInterval& a, const ThreadMarkerInterval& b) {
                        return a.orientedReadId == b.orientedReadId;
                    }),
                currentAnchorMarkers.end());
            markerIntervals.appendVector(currentAnchorMarkers);
            sequences.appendVector(vector<Base>()); // Empty sequence
        }
    }
}


void Anchors::constructFromSelectedMarkerGraphVerticesThreadFunction(uint64_t threadId)
{
    using ThreadMarkerInterval = ConstructFromMarkerGraphData::ThreadMarkerInterval;
    auto& data = constructFromMarkerGraphData;
    DINARA_ASSERT(data.selectedVertexIdsPointer);
    const vector<MarkerGraphVertexId>& selectedVertexIds = *data.selectedVertexIdsPointer;

    const uint64_t minPrimaryCoverage = data.minPrimaryCoverage;
    const uint64_t maxPrimaryCoverage = data.maxPrimaryCoverage;

    const MarkerGraph& markerGraph = *data.markerGraphPointer;
    const auto& markerGraphVertices = markerGraph.vertices();

    // Create the vector to contain the marker intervals for the Anchors found by this thread.
    shared_ptr< MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t> > markerIntervalsPointer =
        make_shared< MemoryMapped::VectorOfVectors<ThreadMarkerInterval, uint64_t> >();
    data.threadMarkerIntervals[threadId] = markerIntervalsPointer;
    auto& markerIntervals = *markerIntervalsPointer;
    markerIntervals.createNew(
        largeDataName("tmp-ThreadAnchorMarkerIntervals-" + to_string(threadId)),
        largeDataPageSize);

    // Create the vector to contain the sequences for the Anchors found by this thread.
    // For vertex-based anchors, sequences are empty.
    shared_ptr< MemoryMapped::VectorOfVectors<Base, uint64_t> > sequencesPointer =
        make_shared< MemoryMapped::VectorOfVectors<Base, uint64_t> >();
    data.threadSequences[threadId] = sequencesPointer;
    auto& sequences = *sequencesPointer;
    sequences.createNew(
        largeDataName("tmp-ThreadAnchorSequences-" + to_string(threadId)),
        largeDataPageSize);

    vector<ThreadMarkerInterval> currentAnchorMarkers;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t i=begin; i!=end; ++i) {
            const MarkerGraphVertexId vertexId = selectedVertexIds[i];
            const MarkerGraphVertexId rcVertexId = markerGraph.reverseComplementVertex[vertexId];
            if(vertexId > rcVertexId) {
                continue; // not canonical
            }

            const auto vertexMarkerIds = markerGraphVertices[vertexId];
            if(vertexMarkerIds.size() < minPrimaryCoverage) {
                continue;
            }

            currentAnchorMarkers.clear();
            currentAnchorMarkers.reserve(vertexMarkerIds.size());
            for(const MarkerId markerId: vertexMarkerIds) {
                OrientedReadId orientedReadId;
                uint32_t ordinal0;
                tie(orientedReadId, ordinal0) = dinara::findMarkerId(markerId, markers);
                ThreadMarkerInterval interval;
                interval.orientedReadId = orientedReadId;
                interval.ordinal0 = ordinal0;
                currentAnchorMarkers.push_back(interval);
            }

            // Sort and deduplicate by oriented read.
            std::sort(currentAnchorMarkers.begin(), currentAnchorMarkers.end(),
                [](const ThreadMarkerInterval& a, const ThreadMarkerInterval& b) {
                    if(a.orientedReadId != b.orientedReadId) {
                        return a.orientedReadId < b.orientedReadId;
                    }
                    return a.ordinal0 < b.ordinal0;
                });
            currentAnchorMarkers.erase(
                std::unique(currentAnchorMarkers.begin(), currentAnchorMarkers.end(),
                    [](const ThreadMarkerInterval& a, const ThreadMarkerInterval& b) {
                        return a.orientedReadId == b.orientedReadId;
                    }),
                currentAnchorMarkers.end());

            if(currentAnchorMarkers.size() < minPrimaryCoverage || currentAnchorMarkers.size() > maxPrimaryCoverage) {
                continue;
            }

            markerIntervals.appendVector(currentAnchorMarkers);
            sequences.appendVector(vector<Base>());

            // Reverse complement anchor.
            for(auto& interval : currentAnchorMarkers) {
                interval.orientedReadId.flipStrand();
                const uint64_t markerCount = markers.size(interval.orientedReadId.getValue());
                interval.ordinal0 = uint32_t(markerCount) - 1 - interval.ordinal0;
            }
            std::sort(currentAnchorMarkers.begin(), currentAnchorMarkers.end(),
                [](const ThreadMarkerInterval& a, const ThreadMarkerInterval& b) {
                    if(a.orientedReadId != b.orientedReadId) {
                        return a.orientedReadId < b.orientedReadId;
                    }
                    return a.ordinal0 < b.ordinal0;
                });
            currentAnchorMarkers.erase(
                std::unique(currentAnchorMarkers.begin(), currentAnchorMarkers.end(),
                    [](const ThreadMarkerInterval& a, const ThreadMarkerInterval& b) {
                        return a.orientedReadId == b.orientedReadId;
                    }),
                currentAnchorMarkers.end());
            markerIntervals.appendVector(currentAnchorMarkers);
            sequences.appendVector(vector<Base>());
        }
    }
}
