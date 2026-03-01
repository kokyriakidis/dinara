#ifndef DINARA_UNITIG_GRAPH_HPP
#define DINARA_UNITIG_GRAPH_HPP

/*
Unitig graph built by compressing the cleaned StringGraph into maximal unitigs.

This follows hifiasm's `ma_ug_gen` construction:
- Unitigs are maximal paths of oriented reads where each internal oriented read has
  out-degree==1 on its suffix side and in-degree==1 on its prefix side (in string-graph sense).
- The unitig graph has 2*unitigCount oriented vertices, and directed arcs between them.
*/

#include "MemoryMappedVector.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

#include <iosfwd>
#include <string>

namespace dinara {
    class UnitigGraph;
    class UnitigGraphArc;
    class UnitigInfo;
}


class dinara::UnitigInfo {
public:
    uint32_t start = uint32_t(invalidReadId); // Oriented StringGraph vertex id.
    uint32_t end = uint32_t(invalidReadId);   // Oriented StringGraph vertex id.
    uint64_t len = 0;                         // Hifiasm-style unitig length (sum of ext + last node length).
    uint32_t n = 0;                           // Number of oriented reads in the unitig path.
    uint8_t circ = 0;                         // 1 if circular (no start/end).
};


class dinara::UnitigGraphArc {
public:
    // Oriented unitig vertex ids (2*unitigId + strand).
    uint32_t from = uint32_t(invalidReadId);
    uint32_t to = uint32_t(invalidReadId);

    // `len` matches hifiasm's unitig-graph arc length (unitigLen - overlapLen, min 1).
    uint32_t len = 0;
    uint32_t overlapLen = 0;

    // Source StringGraph arc id that generated this unitig arc (for debugging).
    uint64_t stringArcId = std::numeric_limits<uint64_t>::max();

    uint8_t del = 0;
};


class dinara::UnitigGraph {
public:
    // Per-unitig metadata and the unitig path (hifiasm ma_utg_t::a entries).
    MemoryMapped::Vector<UnitigInfo> unitigs;
    MemoryMapped::VectorOfVectors<uint64_t, uint32_t> unitigPaths;

    // Directed arcs between oriented unitigs, stored with reverse-complement twin pairs at consecutive positions.
    // Reverse complement arc id is `arcId ^ 1`.
    MemoryMapped::Vector<UnitigGraphArc> arcs;

    // Outgoing / incoming adjacency for oriented unitigs (size 2*unitigCount).
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> outgoing;
    MemoryMapped::VectorOfVectors<uint32_t, uint32_t> incoming;

    // Per-unitig deletion flag (matches hifiasm `asg_t::seq[].del` at unitig level).
    // Size is `unitigCount`.
    MemoryMapped::Vector<uint8_t> unitigDeleted;

    void remove();
    void unreserve();

    void writeGfa(const std::string& fileName) const;
    void writeGfa(std::ostream& gfa) const;
};

#endif

