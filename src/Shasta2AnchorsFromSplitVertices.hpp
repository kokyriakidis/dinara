#pragma once

// Shasta2AnchorsFromSplitVertices.hpp
//
// Creates Shasta2Anchors by splitting marker graph vertices whose reads
// were merged by transitive closure but lack direct pairwise overlaps.
// Each vertex is partitioned into overlap-connected components using the
// read graph, producing anchors that only contain mutually overlapping reads.

#include "Shasta2Anchors.hpp"

#include <memory>
#include <cstdint>

namespace dinara {
    class ReadGraph;

    // Create Shasta2Anchors from marker graph vertices, splitting vertices
    // whose reads are not mutually connected in the read graph.
    // Drop-in replacement for the standard Shasta2Anchors constructor.
    // The read graph must contain only cis overlaps (all edges are confirmed
    // cis from both perspectives after phasing).
    std::shared_ptr<Shasta2Anchors> createShasta2AnchorsFromSplitVertices(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph& markerGraph,
        const ReadGraph& readGraph,
        uint64_t threadCount,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage);
}
