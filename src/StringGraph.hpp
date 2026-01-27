#ifndef DINARA_STRING_GRAPH_HPP
#define DINARA_STRING_GRAPH_HPP

/*
Hifiasm-style string graph (directed arcs between oriented reads).

This is intentionally distinct from `ReadGraph`:
- `ReadGraph` is an undirected connectivity structure used throughout Dinara.
- `StringGraph` is a directed overlap graph meant to match hifiasm/miniasm semantics:
    arc u->v means the suffix (3') of oriented read u overlaps the prefix (5') of oriented read v.
  Each arc has a reverse-complement twin at `arcId^1`: (v^1)->(u^1).
*/

#include "MemoryMappedVectorOfVectors.hpp"
#include "MemoryMappedVector.hpp"
#include "ReadId.hpp"

namespace dinara {
    class StringGraph;
    class StringGraphArc;
}


class dinara::StringGraphArc {
public:
    // Oriented vertex ids (OrientedReadId::getValue()).
    uint32_t from = uint32_t(invalidReadId);
    uint32_t to = uint32_t(invalidReadId);

    // `len` matches hifiasm's `asg_arc_len` (the non-overlap extension length along `from`).
    // `overlapLen` matches hifiasm's `p->ol` (overlap length).
    uint32_t len = 0;
    uint32_t overlapLen = 0;

    // Alignment id that generated this arc (for debugging / tracing).
    uint64_t alignmentId = std::numeric_limits<uint64_t>::max();

    // Reserved for later cleaning stages (transitive reduction, etc).
    uint8_t del = 0;
};


class dinara::StringGraph {
public:
    // Arcs are stored with reverse-complement twin pairs at consecutive positions.
    // Reverse complement arc id is `arcId ^ 1`.
    MemoryMapped::Vector<StringGraphArc> arcs;

    // Outgoing adjacency: for each oriented read (vertex), arc ids of arcs leaving that vertex.
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> outgoing;

    // Incoming adjacency: for each oriented read (vertex), arc ids of arcs entering that vertex.
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> incoming;

    // Per-read deletion flag used by string-graph cleaning (matches `asg_t::seq[].del` in hifiasm).
    // Size is `readCount` (not `2*readCount`).
    MemoryMapped::Vector<uint8_t> readDeleted;

    uint64_t getReverseComplementArcId(uint64_t arcId) const;
    void unreserve();
    void remove();
};

#endif
