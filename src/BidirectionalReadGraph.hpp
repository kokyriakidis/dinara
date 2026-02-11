#ifndef DINARA_BIDIRECTIONAL_READ_GRAPH_HPP
#define DINARA_BIDIRECTIONAL_READ_GRAPH_HPP

/*******************************************************************************

BidirectionalReadGraph — A single-vertex-per-read overlap graph
===============================================================

Motivation
----------
The legacy ReadGraph is "strand-doubled": every physical read has two vertices
(strand 0 and strand 1), and every alignment produces two edges (one + its RC).
This works well in collinear regions but creates problems at strand-contact
boundaries (inversions, segmental duplications, palindromes):

  1. Strand-separation (flagCrossStrandReadGraphEdges) rejects biologically
     correct overlaps that happen to span an inversion breakpoint.
  2. Downstream BFS (computeReadGraphStrandsFromSeed) encounters conflicting
     strand assignments via different paths at inversions.
  3. Anchor creation skips cross-strand edges, causing coverage gaps and
     fragmented unitigs at inversion breakpoints.

Data model
----------
  Vertices : N  (one per physical read, indexed by ReadId)
  Edges    : A  (one per kept alignment — not 2A)
  Each edge stores:
    - readIds[2]   : the two reads (readIds[0] < readIds[1], same canonicalization as OrientedReadPair)
    - isSameStrand : relative orientation (true = forward/forward, false = forward/reverse-complement)
    - alignmentId  : back-pointer into alignmentData[]

Orientation-aware traversal
---------------------------
Every traversal carries a "running strand" s ∈ {0,1} that tracks the current
read's inferred orientation:

    Enter read u with derived strand s.
    For every edge (u, v, sameStrand):
        derived strand of v  =  s ^ (!sameStrand)

This naturally handles inversions: a same-strand edge preserves orientation,
a cross-strand edge flips it.  No edges need to be discarded.

Relationship to existing ReadGraph
----------------------------------
Phase 0 installs BidirectionalReadGraph as a parallel data structure alongside
the existing ReadGraph.  No existing code is modified.  Later phases will
migrate consumers from ReadGraph → BidirectionalReadGraph one at a time.

*******************************************************************************/

// Dinara.
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

// Standard library.
#include "cstdint.hpp"
#include "vector.hpp"
#include <limits>
#include <queue>
#include <map>


namespace dinara {
    class BidirectionalReadGraph;
    class BidirectionalReadGraphEdge;
}



// ============================================================================
// BidirectionalReadGraphEdge
// ============================================================================
// One per alignment.  Canonical form: readIds[0] < readIds[1].
// Same canonicalization as OrientedReadPair / AlignmentData.
class dinara::BidirectionalReadGraphEdge {
public:
    array<ReadId, 2> readIds;

    // The id of the alignment that corresponds to this edge.
    uint64_t alignmentId : 61;

    // Relative strand of the two reads.
    // 1 = same strand (forward/forward), 0 = different strands (forward/RC).
    uint64_t isSameStrand : 1;

    // Soft-delete flag for iterative cleaning without reallocating.
    uint64_t isDeleted : 1;

    // Additional flags carried over from the legacy graph for compatibility.
    uint64_t hasInconsistentAlignment : 1;

    BidirectionalReadGraphEdge()
    {
        readIds = {invalidReadId, invalidReadId};
        alignmentId = 0;
        isSameStrand = 1;
        isDeleted = 0;
        hasInconsistentAlignment = 0;
    }

    // -----------------------------------------------------------------------
    // Traversal helpers
    // -----------------------------------------------------------------------

    /// Given one of the two read IDs, return the other.
    ReadId getOther(ReadId readId) const
    {
        if(readId == readIds[0]) {
            return readIds[1];
        } else {
            DINARA_ASSERT(readId == readIds[1]);
            return readIds[0];
        }
    }

    /// Orientation-aware neighbour query.
    /// Given that we entered `fromReadId` with derived strand `fromStrand`,
    /// returns (toReadId, toStrand) where toStrand accounts for the edge's
    /// relative orientation.
    ///
    ///   toStrand = fromStrand                  if isSameStrand
    ///   toStrand = fromStrand ^ 1              if !isSameStrand
    ///
    std::pair<ReadId, Strand> traverse(ReadId fromReadId, Strand fromStrand) const
    {
        const ReadId toReadId = getOther(fromReadId);
        const Strand toStrand = isSameStrand ? fromStrand : (fromStrand ^ 1);
        return {toReadId, toStrand};
    }
};



// ============================================================================
// BidirectionalReadGraph
// ============================================================================
// Undirected overlap graph with one vertex per physical read.
//
// Edges are stored in a flat vector.  Adjacency is indexed by ReadId
// (size = readCount), not OrientedReadId::getValue() (size = 2*readCount).
//
// This halves memory for the connectivity index and edge storage compared
// to the strand-doubled ReadGraph.
class dinara::BidirectionalReadGraph {
public:

    // The flat edge array.
    MemoryMapped::Vector<BidirectionalReadGraphEdge> edges;

    // Adjacency lists: connectivity[readId] → span of edge indices incident to readId.
    // Indexed by ReadId (range [0, readCount)).
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> connectivity;

    // Number of (non-deleted) vertices.  Set during buildConnectivity().
    uint64_t readCount() const { return connectivity.size(); }

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Build the connectivity index from the current edges[] array.
    /// Must be called after all edges have been pushed.
    /// @param nReads  Total number of physical reads.
    void buildConnectivity(uint64_t nReads);

    // -----------------------------------------------------------------------
    // Traversal
    // -----------------------------------------------------------------------

    /// Find all neighbours of a read (distance 1), regardless of strand.
    /// Skips deleted edges.  Results are sorted by ReadId.
    void findNeighbors(
        ReadId readId,
        vector<ReadId>& neighbors) const;

    /// Find all neighbours within maxDistance hops.
    /// Skips deleted edges.  Results are sorted by ReadId.
    void findNeighbors(
        ReadId readId,
        uint64_t maxDistance,
        vector<ReadId>& neighbors) const;

    /// BFS shortest path from read0 to read1, skipping deleted edges.
    /// Returns edge IDs of the path, or an empty vector if unreachable
    /// within maxDistance.
    static const uint32_t infiniteDistance;
    void computeShortPath(
        ReadId readId0,
        ReadId readId1,
        size_t maxDistance,
        vector<uint32_t>& path,
        // Work areas (sized readCount, caller responsibility):
        vector<uint32_t>& distance,
        vector<ReadId>& reachedVertices,
        vector<uint32_t>& parentEdges) const;

    /// Orientation-aware BFS:  starting from (seedReadId, seedStrand=0),
    /// propagate derived strand to all reachable reads.
    /// Returns a vector indexed by ReadId:
    ///   -1 = unreachable, 0 or 1 = derived strand.
    /// Also returns the number of BFS conflicts (a read reached via two
    /// paths that imply different strands — indicates inversion boundary).
    struct StrandResult {
        vector<int8_t> strandByRead;   // indexed by ReadId
        uint64_t conflictCount = 0;
    };
    StrandResult propagateStrands(
        ReadId seedReadId,
        Strand seedStrand = 0) const;

    // -----------------------------------------------------------------------
    // Lifetime
    // -----------------------------------------------------------------------
    void unreserve();
    void remove();
};



#endif

