// AssemblerBidirectionalReadGraph.cpp
//
// Construction and access methods for the BidirectionalReadGraph.
//
// The BidirectionalReadGraph stores one vertex per physical read and one edge
// per kept alignment (no strand doubling).  Every edge records the relative
// orientation (isSameStrand), allowing orientation-aware traversal without
// discarding cross-strand edges.
//
// Phase 1 of the BRG migration: populate from existing alignmentData[].

#include "Assembler.hpp"
#include "LocalReadGraph.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"

#include "chrono.hpp"
#include <queue>

using namespace dinara;



// ============================================================================
// createBidirectionalReadGraph
// ============================================================================
// Populates the BidirectionalReadGraph from alignmentData[].
//
// An alignment is included iff it is currently marked as in the read graph
// (alignmentData[i].info.isInReadGraph == 1).  This means the BRG should be
// built *after* whatever filtering pipeline has run (createReadGraph6,
// createReadGraphFromFilteredAlignments, etc.) so that isInReadGraph reflects
// the final kept-set.
//
// Unlike the strand-doubled ReadGraph which produces 2 edges per alignment,
// the BRG produces exactly 1 edge per alignment.
void Assembler::createBidirectionalReadGraph()
{
    cout << timestamp << "createBidirectionalReadGraph begins." << endl;
    checkAlignmentDataAreOpen();

    const uint64_t alignmentCount = alignmentData.size();
    const uint64_t nReads = reads->readCount();

    // Remove any existing BRG data structures to start clean.
    bidirectionalReadGraph.remove();

    // Create the edges vector.
    bidirectionalReadGraph.edges.createNew(
        largeDataName("BidirectionalReadGraphEdges"), largeDataPageSize);

    uint64_t keptCount = 0;
    uint64_t skippedCount = 0;

    for(uint64_t alignmentId = 0; alignmentId < alignmentCount; ++alignmentId) {
        const AlignmentData& ad = alignmentData[alignmentId];

        // Only include alignments that are in the read graph.
        if(!ad.info.isInReadGraph) {
            ++skippedCount;
            continue;
        }

        BidirectionalReadGraphEdge edge;
        edge.readIds[0] = ad.readIds[0];
        edge.readIds[1] = ad.readIds[1];
        edge.alignmentId = alignmentId & 0x1fff'ffff'ffff'ffff;  // 61-bit mask
        edge.isSameStrand = ad.isSameStrand ? 1 : 0;
        edge.isDeleted = 0;
        edge.hasInconsistentAlignment = 0;

        bidirectionalReadGraph.edges.push_back(edge);
        ++keptCount;
    }

    // Release over-allocated capacity.
    bidirectionalReadGraph.edges.unreserve();

    cout << timestamp << "BidirectionalReadGraph: " << keptCount
         << " edges from " << alignmentCount << " alignments ("
         << skippedCount << " skipped)." << endl;

    // Build adjacency index.
    bidirectionalReadGraph.connectivity.createNew(
        largeDataName("BidirectionalReadGraphConnectivity"), largeDataPageSize);
    bidirectionalReadGraph.buildConnectivity(nReads);

    // Count isolated reads (no BRG edges).
    uint64_t isolatedReadCount = 0;
    for(ReadId readId = 0; readId < nReads; ++readId) {
        bool hasEdge = false;
        for(const uint32_t edgeId : bidirectionalReadGraph.connectivity[readId]) {
            if(!bidirectionalReadGraph.edges[edgeId].isDeleted) {
                hasEdge = true;
                break;
            }
        }
        if(!hasEdge) {
            ++isolatedReadCount;
        }
    }

    // Summary statistics.
    cout << timestamp << "BidirectionalReadGraph: "
         << nReads << " reads, "
         << keptCount << " edges, "
         << isolatedReadCount << " isolated reads." << endl;

    // Comparison with the strand-doubled ReadGraph.
    if(readGraph.edges.isOpen) {
        const uint64_t strandDoubledEdges = readGraph.edges.size();
        cout << timestamp << "BidirectionalReadGraph vs ReadGraph: "
             << keptCount << " BRG edges, "
             << strandDoubledEdges << " strand-doubled edges ("
             << strandDoubledEdges / 2 << " unique alignments)." << endl;
    }

    cout << timestamp << "createBidirectionalReadGraph completed." << endl;
}



// ============================================================================
// createBidirectionalReadGraphFromSelectedAlignments
// ============================================================================
// Populates the BidirectionalReadGraph from a caller-provided keep vector.
// This allows building the BRG with a custom selection of alignments
// (independent of the isInReadGraph flag).
void Assembler::createBidirectionalReadGraphFromSelectedAlignments(
    const vector<bool>& keepAlignment)
{
    cout << timestamp << "createBidirectionalReadGraphFromSelectedAlignments begins." << endl;
    checkAlignmentDataAreOpen();

    const uint64_t alignmentCount = alignmentData.size();
    const uint64_t nReads = reads->readCount();
    DINARA_ASSERT(keepAlignment.size() == alignmentCount);

    // Remove any existing BRG data structures.
    bidirectionalReadGraph.remove();

    // Create edge vector.
    bidirectionalReadGraph.edges.createNew(
        largeDataName("BidirectionalReadGraphEdges"), largeDataPageSize);

    uint64_t keptCount = 0;

    for(uint64_t alignmentId = 0; alignmentId < alignmentCount; ++alignmentId) {
        if(!keepAlignment[alignmentId]) {
            continue;
        }

        const AlignmentData& ad = alignmentData[alignmentId];

        BidirectionalReadGraphEdge edge;
        edge.readIds[0] = ad.readIds[0];
        edge.readIds[1] = ad.readIds[1];
        edge.alignmentId = alignmentId & 0x1fff'ffff'ffff'ffff;  // 61-bit mask
        edge.isSameStrand = ad.isSameStrand ? 1 : 0;
        edge.isDeleted = 0;
        edge.hasInconsistentAlignment = 0;

        bidirectionalReadGraph.edges.push_back(edge);
        ++keptCount;
    }

    bidirectionalReadGraph.edges.unreserve();

    cout << timestamp << "BidirectionalReadGraph: " << keptCount
         << " edges from " << alignmentCount << " alignments." << endl;

    // Build adjacency index.
    bidirectionalReadGraph.connectivity.createNew(
        largeDataName("BidirectionalReadGraphConnectivity"), largeDataPageSize);
    bidirectionalReadGraph.buildConnectivity(nReads);

    cout << timestamp << "createBidirectionalReadGraphFromSelectedAlignments completed." << endl;
}



// ============================================================================
// Access methods
// ============================================================================
void Assembler::accessBidirectionalReadGraph()
{
    bidirectionalReadGraph.edges.accessExistingReadOnly(
        largeDataName("BidirectionalReadGraphEdges"));
    bidirectionalReadGraph.connectivity.accessExistingReadOnly(
        largeDataName("BidirectionalReadGraphConnectivity"));
}

void Assembler::accessBidirectionalReadGraphReadWrite()
{
    bidirectionalReadGraph.edges.accessExistingReadWrite(
        largeDataName("BidirectionalReadGraphEdges"));
    bidirectionalReadGraph.connectivity.accessExistingReadWrite(
        largeDataName("BidirectionalReadGraphConnectivity"));
}

void Assembler::checkBidirectionalReadGraphIsOpen() const
{
    if(!bidirectionalReadGraph.edges.isOpen) {
        throw runtime_error("Bidirectional read graph edges are not accessible.");
    }
    if(!bidirectionalReadGraph.connectivity.isOpen()) {
        throw runtime_error("Bidirectional read graph connectivity is not accessible.");
    }
}



// ============================================================================
// removeBidirectionalReadGraph
// ============================================================================
void Assembler::removeBidirectionalReadGraph()
{
    bidirectionalReadGraph.remove();
}



// ============================================================================
// createLocalBidirectionalReadGraph
// ============================================================================
// Create a local subgraph of the BidirectionalReadGraph for visualization.
//
// This performs a BFS on the BRG (one vertex per physical read), but emits
// a LocalReadGraph with OrientedReadId vertices so the result can be rendered
// with the same pipeline as the legacy ReadGraph explorer.
//
// Each physical read appears exactly once in the output graph.
// Visited status is tracked by ReadId, not OrientedReadId.
// The displayed strand is the "derived strand" from the BFS — i.e. the
// inferred orientation relative to the seed read.
//
// Orientation is derived via edge.traverse():
//   - Each start read enters with its user-supplied strand.
//   - Edges propagate strand through: toStrand = fromStrand ^ (!isSameStrand).
//
// The "crossesStrands" flag on LocalReadGraphEdge is set to !isSameStrand,
// which causes those edges to be drawn in orange (same coloring convention
// as the legacy ReadGraph viewer).
bool Assembler::createLocalBidirectionalReadGraph(
    const vector<OrientedReadId>& starts,
    uint32_t maxDistance,
    bool allowChimericReads,
    bool allowInconsistentAlignmentEdges,
    double timeout,
    LocalReadGraph& graph)
{
    using std::chrono::steady_clock;
    const auto startTime = steady_clock::now();

    checkBidirectionalReadGraphIsOpen();

    // Track which physical reads have been visited.
    // Maps ReadId → the OrientedReadId assigned in the LocalReadGraph.
    // This ensures each physical read appears exactly once, regardless of
    // how many paths reach it or which derived strand each path implies.
    std::map<ReadId, OrientedReadId> visitedReads;

    std::queue<OrientedReadId> q;

    for(const OrientedReadId& start : starts) {
        const ReadId readId = start.getReadId();
        const Strand strand = start.getStrand();

        // If the starting read is chimeric and we don't allow chimeric reads, skip.
        if(!allowChimericReads && reads->getFlags(readId).isChimeric) {
            continue;
        }

        // Skip if this physical read was already added by a prior start vertex.
        if(visitedReads.count(readId)) {
            continue;
        }

        const OrientedReadId orientedReadId(readId, strand);
        graph.addVertex(
            orientedReadId,
            uint32_t((*markers)[orientedReadId.getValue()].size()),
            reads->getFlags(readId).isChimeric,
            0);
        visitedReads[readId] = orientedReadId;
        q.push(orientedReadId);
    }

    // BFS.
    while(!q.empty()) {

        // Check timeout.
        if(timeout > 0. && (seconds(steady_clock::now() - startTime) > timeout)) {
            graph.clear();
            return false;
        }

        // Dequeue.
        const OrientedReadId orientedReadId0 = q.front();
        q.pop();
        const ReadId readId0 = orientedReadId0.getReadId();
        const Strand strand0 = orientedReadId0.getStrand();
        const uint32_t distance0 = graph.getDistance(orientedReadId0);

        // Loop over BRG edges incident to this read.
        for(const uint32_t edgeIndex : bidirectionalReadGraph.connectivity[readId0]) {
            DINARA_ASSERT(edgeIndex < bidirectionalReadGraph.edges.size());
            const BidirectionalReadGraphEdge& edge =
                bidirectionalReadGraph.edges[edgeIndex];

            // Skip deleted edges.
            if(edge.isDeleted) {
                continue;
            }
            if(!allowInconsistentAlignmentEdges && edge.hasInconsistentAlignment) {
                continue;
            }

            // Orientation-aware traversal: derive the neighbor's strand.
            const auto [toReadId, toStrand] = edge.traverse(readId0, strand0);

            // If this read is chimeric and we don't allow chimeric reads, skip.
            if(!allowChimericReads && reads->getFlags(toReadId).isChimeric) {
                continue;
            }

            // Get alignment information for markerCount.
            const AlignmentData& alignment = alignmentData[edge.alignmentId];
            OrientedReadId alignmentOrientedReadId0(alignment.readIds[0], 0);
            OrientedReadId alignmentOrientedReadId1(
                alignment.readIds[1], alignment.isSameStrand ? 0 : 1);
            AlignmentInfo alignmentInfo = alignment.info;
            if(alignmentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
                swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
                alignmentInfo.swap();
            }
            if(alignmentOrientedReadId0.getStrand() != orientedReadId0.getStrand()) {
                alignmentOrientedReadId0.flipStrand();
                alignmentOrientedReadId1.flipStrand();
                alignmentInfo.reverseComplement();
            }
            DINARA_ASSERT(alignmentOrientedReadId0 == orientedReadId0);
            const uint32_t markerCount = alignmentInfo.markerCount;

            // The "crossesStrands" visual flag is set for cross-strand BRG edges
            // so they are drawn in orange, matching the legacy coloring convention.
            const bool crossesStrands = !edge.isSameStrand;

            // BFS expansion — track by ReadId so each physical read appears once.
            const uint32_t distance1 = distance0 + 1;
            if(distance0 < maxDistance) {
                if(visitedReads.count(toReadId) == 0) {
                    // First time visiting this physical read.
                    const OrientedReadId orientedReadId1(toReadId, toStrand);
                    graph.addVertex(
                        orientedReadId1,
                        uint32_t((*markers)[orientedReadId1.getValue()].size()),
                        reads->getFlags(toReadId).isChimeric,
                        distance1);
                    visitedReads[toReadId] = orientedReadId1;
                    q.push(orientedReadId1);
                }
                // Use the canonical OrientedReadId assigned to this read.
                graph.addEdge(
                    orientedReadId0,
                    visitedReads[toReadId],
                    markerCount,
                    edgeIndex,
                    crossesStrands);
            } else {
                DINARA_ASSERT(distance0 == maxDistance);
                if(visitedReads.count(toReadId)) {
                    graph.addEdge(
                        orientedReadId0,
                        visitedReads[toReadId],
                        markerCount,
                        edgeIndex,
                        crossesStrands);
                }
            }
        }
    }

    return true;
}
