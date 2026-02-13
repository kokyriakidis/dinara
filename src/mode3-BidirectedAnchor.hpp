#ifndef DINARA_MODE3_BIDIRECTED_ANCHOR_HPP
#define DINARA_MODE3_BIDIRECTED_ANCHOR_HPP

/******************************************************************************

BRG-native anchor and journey data structures.

In the bidirected read graph (BRG), there is no separate forward and reverse
complement representation.  Reads are physical (one per ReadId) and overlaps
are undirected edges with a same-strand/diff-strand flag.

Traditional dinara anchors operate in OrientedReadId space: each anchor has a
forward copy and a reverse-complement copy (AnchorId XOR 1), and journeys are
per-OrientedReadId.

BidirectedAnchors operate in ReadId space:
  - Each anchor is a single entity (no RC doubling).
  - Marker intervals store ReadId + ordinal + the strand from which the marker
    was observed.  This "observation strand" serves as the direction: when
    a read traverses an anchor, the strand tells you which way.
  - Journeys are per-ReadId.  Each journey entry is (BidirectedAnchorId, strand),
    where strand encodes the direction of traversal:
      strand 0  →  the read's strand-0 markers hit this anchor
      strand 1  →  the read's strand-1 markers hit this anchor
    Consecutive entries in a journey define edges in the bidirected anchor
    graph.  The pair of strands at the endpoints encodes the edge type
    (same-strand = forward, diff-strand = inversion).

Primitives inspired by MBG (Verkko's graph builder):
  - OrientedBidirectedAnchor = (BidirectedAnchorId, Strand) — an anchor with direction.
  - reverse(o): flip the strand → access the "other end" of the node.
  - canon(from, to): canonical edge representation for deduplication.
  - revCompPath(path): reverse path + flip all strands.
  - BidirectedVectorWithDirection<T>: container indexed by (anchorId, strand),
    with separate forward/backward storage per node end.

This header defines the data structures.  The algorithms (creation from BRG
marker graph, journey computation, neighbor finding, edge computation) live
in mode3-BidirectedAnchor.cpp and AssemblerBidirectedAnchors.cpp.

******************************************************************************/

#include "dinaraTypes.hpp"
#include "invalid.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MemoryMappedVector.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

#include <algorithm>
#include <cassert>
#include <map>
#include <span>
#include <utility>
#include <vector>

namespace dinara {

    class CompressedMarker;
    class MarkerGraph;
    class Reads;

    namespace mode3 {

        using BidirectedAnchorId = uint64_t;

        class BidirectedAnchorMarkerInterval;
        class BidirectedAnchorInfo;
        class BidirectedJourneyEntry;
        class BidirectedAnchors;


        // ===============================================================
        // OrientedBidirectedAnchor — a directed node in the bidirected graph.
        // Equivalent to MBG's std::pair<size_t, bool>.
        // ===============================================================

        struct OrientedBidirectedAnchor {
            BidirectedAnchorId anchorId;
            Strand strand;  // 0 = forward, 1 = reverse

            OrientedBidirectedAnchor() :
                anchorId(invalid<BidirectedAnchorId>), strand(0) {}

            OrientedBidirectedAnchor(BidirectedAnchorId anchorId, Strand strand) :
                anchorId(anchorId), strand(strand) {}

            bool operator==(const OrientedBidirectedAnchor& o) const {
                return anchorId == o.anchorId && strand == o.strand;
            }
            bool operator!=(const OrientedBidirectedAnchor& o) const {
                return !(*this == o);
            }
            bool operator<(const OrientedBidirectedAnchor& o) const {
                if(anchorId != o.anchorId) return anchorId < o.anchorId;
                return strand < o.strand;
            }

            bool isValid() const { return anchorId != invalid<BidirectedAnchorId>; }
        };


        // ===============================================================
        // Core bidirected-graph primitives (MBG-style).
        // ===============================================================

        // Flip the orientation of a directed anchor.
        // Provides access to the "other end" of a bidirected node.
        inline OrientedBidirectedAnchor reverse(OrientedBidirectedAnchor o) {
            return OrientedBidirectedAnchor(o.anchorId, 1 - o.strand);
        }

        // Canonical representation of a directed edge.
        // Picks the lexicographic minimum of (from, to) vs (reverse(to), reverse(from)).
        // This ensures that an edge and its reverse-complement counterpart
        // map to the same canonical form, like MBG's canon().
        inline std::pair<OrientedBidirectedAnchor, OrientedBidirectedAnchor>
        canon(OrientedBidirectedAnchor from, OrientedBidirectedAnchor to) {
            OrientedBidirectedAnchor revTo = reverse(to);
            OrientedBidirectedAnchor revFrom = reverse(from);
            // Compare (from, to) with (reverse(to), reverse(from))
            if(revTo < from || (revTo == from && revFrom < to)) {
                return {revTo, revFrom};
            }
            return {from, to};
        }

        // Reverse-complement an entire path (sequence of oriented anchors).
        // Reverses the order and flips each orientation.
        inline std::vector<OrientedBidirectedAnchor>
        revCompPath(const std::vector<OrientedBidirectedAnchor>& path) {
            std::vector<OrientedBidirectedAnchor> result;
            result.reserve(path.size());
            for(auto it = path.rbegin(); it != path.rend(); ++it) {
                result.push_back(reverse(*it));
            }
            return result;
        }


        // ===============================================================
        // BidirectedVectorWithDirection<T> — a container indexed by
        // OrientedBidirectedAnchor (anchorId, strand), with separate storage
        // for forward (strand=0) and reverse (strand=1) entries.
        // Inspired by MBG's VectorWithDirection<T>.
        // ===============================================================

        template<class T>
        class BidirectedVectorWithDirection {
        public:
            void resize(uint64_t nodeCount) {
                forward.resize(nodeCount);
                backward.resize(nodeCount);
            }
            void resize(uint64_t nodeCount, const T& value) {
                forward.resize(nodeCount, value);
                backward.resize(nodeCount, value);
            }
            uint64_t size() const {
                return forward.size();
            }
            void clear() {
                forward.clear();
                backward.clear();
            }
            void emplace_back() {
                forward.emplace_back();
                backward.emplace_back();
            }
            void emplace_back(const T& value) {
                forward.emplace_back(value);
                backward.emplace_back(value);
            }

            T& operator[](OrientedBidirectedAnchor o) {
                return o.strand == 0 ? forward[o.anchorId] : backward[o.anchorId];
            }
            const T& operator[](OrientedBidirectedAnchor o) const {
                return o.strand == 0 ? forward[o.anchorId] : backward[o.anchorId];
            }

        private:
            std::vector<T> forward;
            std::vector<T> backward;
        };


        // ===============================================================
        // BidirectedEdge — a directed edge between two oriented anchors.
        // The edge encodes the bidirected connectivity: the pair of
        // strands at the two endpoints determines the edge type.
        // ===============================================================

        struct BidirectedEdge {
            OrientedBidirectedAnchor from;
            OrientedBidirectedAnchor to;
            uint64_t coverage = 0;  // number of reads supporting this edge

            BidirectedEdge() {}
            BidirectedEdge(OrientedBidirectedAnchor from, OrientedBidirectedAnchor to, uint64_t coverage = 0) :
                from(from), to(to), coverage(coverage) {}

            bool operator==(const BidirectedEdge& o) const {
                return from == o.from && to == o.to;
            }
            bool operator<(const BidirectedEdge& o) const {
                if(from != o.from) return from < o.from;
                return to < o.to;
            }
        };

    }
}


// A marker interval in a BRG anchor.
// Identifies a physical read, its marker ordinal, and the strand
// from which this marker participates in the anchor.
class dinara::mode3::BidirectedAnchorMarkerInterval {
public:
    ReadId readId;
    uint32_t ordinal;       // Marker ordinal on strand 0 of the read.
    Strand  strand;         // 0 or 1: which strand's marker is in this anchor.

    BidirectedAnchorMarkerInterval() {}
    BidirectedAnchorMarkerInterval(ReadId readId, uint32_t ordinal, Strand strand) :
        readId(readId), ordinal(ordinal), strand(strand) {}
};


// Per-anchor metadata.
class dinara::mode3::BidirectedAnchorInfo {
public:
    uint32_t componentId = invalid<uint32_t>;
    uint64_t localAnchorIdInComponent = invalid<uint64_t>;
};


// A single entry in a per-ReadId journey.
// Records which BRG anchor was visited and which strand the read uses there.
// This is equivalent to an OrientedBidirectedAnchor but stored as a journey entry.
class dinara::mode3::BidirectedJourneyEntry {
public:
    BidirectedAnchorId anchorId;
    Strand strand;          // direction of traversal (0 = forward, 1 = reverse)

    BidirectedJourneyEntry() {}
    BidirectedJourneyEntry(BidirectedAnchorId anchorId, Strand strand) :
        anchorId(anchorId), strand(strand) {}

    // Convert to/from OrientedBidirectedAnchor.
    BidirectedJourneyEntry(OrientedBidirectedAnchor o) :
        anchorId(o.anchorId), strand(o.strand) {}
    OrientedBidirectedAnchor oriented() const {
        return OrientedBidirectedAnchor(anchorId, strand);
    }
};


// The main BRG-native anchor container.
// Owns anchors, journeys, and a precomputed bidirected edge table.
class dinara::mode3::BidirectedAnchors :
    public MultithreadedObject<BidirectedAnchors>,
    public MappedMemoryOwner {
public:

    // ---- Construction ----

    // Create BidirectedAnchors from self-RC marker graph vertices produced by
    // createMarkerGraphVerticesFromBrg (with forward/RC kmer merging).
    // Each self-RC vertex becomes one BrgAnchor.
    // Non-self-RC vertices are also handled (should be rare after merging).
    BidirectedAnchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph& markerGraph,
        uint64_t minCoverage,
        uint64_t maxCoverage,
        uint64_t threadCount);

    // Access existing BidirectedAnchors from binary data (for HTTP server).
    // Loads anchorMarkerIntervals, anchorInfos, and journeys from disk,
    // then recomputes the in-memory edge table.
    BidirectedAnchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers);

    // ---- Accessors ----

    uint64_t size() const { return anchorCount; }

    // Get the marker intervals for a given anchor.
    std::span<const BidirectedAnchorMarkerInterval> operator[](BidirectedAnchorId anchorId) const;

    // Get the coverage (number of distinct physical reads) of an anchor.
    uint64_t coverage(BidirectedAnchorId anchorId) const;


    // ---- Journeys ----

    // Compute per-ReadId journeys.  Each journey is the ordered sequence of
    // (BidirectedAnchorId, strand) pairs along a physical read, sorted by ordinal.
    void computeJourneys(uint64_t threadCount);

    // Access the journey for a given ReadId.
    std::span<const BidirectedJourneyEntry> journey(ReadId readId) const;

    // Write journeys to CSV for diagnostics.
    void writeJourneys() const;


    // ---- Edge table ----

    // Precompute the bidirected edge table from journeys.
    // For every consecutive pair of journey entries (A,sA) → (B,sB),
    // an edge is recorded from OrientedBidirectedAnchor(A,sA) to OrientedBidirectedAnchor(B,sB),
    // AND the symmetric reverse edge from reverse(B,sB) to reverse(A,sA).
    // Edge coverage counts how many reads support the edge.
    void computeEdges(uint64_t threadCount);

    // Get all outgoing edges from a directed anchor (O(1) lookup).
    std::span<const BidirectedEdge> getOutEdges(OrientedBidirectedAnchor oa) const;

    // Get all incoming edges to a directed anchor (O(1) lookup via reverse).
    // Incoming edges to (a, s) = outgoing edges from reverse(a, s), with
    // source/target reversed.
    std::span<const BidirectedEdge> getInEdges(OrientedBidirectedAnchor oa) const;

    // Write edges to CSV for diagnostics.
    void writeEdges() const;


    // ---- Neighbor finding (from precomputed edges) ----

    // Convenience structures and functions using the edge table.
    struct Neighbor {
        BidirectedAnchorId anchorId;
        Strand strand;
        uint64_t count;
    };

    // Get forward neighbors of an oriented anchor from the edge table.
    void findForwardNeighbors(
        OrientedBidirectedAnchor oa,
        std::vector<Neighbor>& neighbors,
        uint64_t minCount = 0) const;

    // Get backward neighbors of an oriented anchor from the edge table.
    void findBackwardNeighbors(
        OrientedBidirectedAnchor oa,
        std::vector<Neighbor>& neighbors,
        uint64_t minCount = 0) const;

    // Legacy overloads for compatibility.
    void findForwardNeighbors(
        BidirectedAnchorId anchorId,
        Strand strandAtAnchor,
        std::vector<Neighbor>& neighbors,
        uint64_t minCount = 0) const {
        findForwardNeighbors(OrientedBidirectedAnchor(anchorId, strandAtAnchor), neighbors, minCount);
    }
    void findBackwardNeighbors(
        BidirectedAnchorId anchorId,
        Strand strandAtAnchor,
        std::vector<Neighbor>& neighbors,
        uint64_t minCount = 0) const {
        findBackwardNeighbors(OrientedBidirectedAnchor(anchorId, strandAtAnchor), neighbors, minCount);
    }


    // ---- Common read counting ----

    // Count common physical reads between two anchors.
    uint64_t countCommon(BidirectedAnchorId, BidirectedAnchorId) const;


    // ---- Data members ----

    const Reads& reads;
    uint64_t k;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;

private:
    uint64_t anchorCount = 0;

    // Marker intervals for each anchor.
    // Indexed by BidirectedAnchorId.
    MemoryMapped::VectorOfVectors<BidirectedAnchorMarkerInterval, uint64_t> anchorMarkerIntervals;

    // Per-anchor metadata.
    MemoryMapped::Vector<BidirectedAnchorInfo> anchorInfos;

    // Per-ReadId journeys (sorted by ordinal within each read).
    MemoryMapped::VectorOfVectors<BidirectedJourneyEntry, uint64_t> journeys;

    // ---- Edge table storage ----
    // Indexed by OrientedBidirectedAnchor via BidirectedVectorWithDirection.
    // Each entry is a vector of outgoing BidirectedEdges.
    // The symmetry invariant is enforced: if edge (A→B) exists,
    // then edge (reverse(B)→reverse(A)) also exists with same coverage.
    BidirectedVectorWithDirection<std::vector<BidirectedEdge>> outEdges;

    // Reverse edges — indexed by OrientedBidirectedAnchor, stores incoming edges.
    // inEdges[oa] contains all edges whose .to == oa.
    BidirectedVectorWithDirection<std::vector<BidirectedEdge>> inEdges;

    bool edgesComputed = false;

    // Temporary storage used during journey computation.
    struct JourneyWithOrdinal {
        BidirectedAnchorId anchorId;
        Strand strand;
        uint32_t ordinal;  // for sorting
    };
    MemoryMapped::VectorOfVectors<JourneyWithOrdinal, uint64_t> journeysWithOrdinals;

    // Thread functions for journey computation.
    void computeJourneysThreadFunction12(uint64_t pass);
    void computeJourneysThreadFunction1(uint64_t threadId);
    void computeJourneysThreadFunction2(uint64_t threadId);
    void computeJourneysThreadFunction3(uint64_t threadId);
    void computeJourneysThreadFunction4(uint64_t threadId);
};

#endif
