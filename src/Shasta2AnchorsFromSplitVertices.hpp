#pragma once

// Shasta2AnchorsFromSplitVertices.hpp
//
// Creates Shasta2Anchors by splitting marker graph vertices whose reads
// were merged by transitive closure but lack direct pairwise overlaps.
// Each vertex is partitioned into overlap-connected components using the
// read graph, producing anchors that only contain mutually overlapping reads.
//
// When alignmentData and alignmentTable are provided, the split also
// validates that both reads' ordinals at the vertex fall within their
// pairwise alignment's ordinal range (prevents false grouping from
// transitive closure contamination).

#include "Shasta2Anchors.hpp"
#include "Alignment.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

#include <memory>
#include <cstdint>

namespace dinara {
    class ReadGraph;

    std::shared_ptr<Shasta2Anchors> createShasta2AnchorsFromSplitVertices(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph& markerGraph,
        const ReadGraph& readGraph,
        uint64_t threadCount,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        const MemoryMapped::Vector<AlignmentData>* alignmentData = nullptr,
        const MemoryMapped::VectorOfVectors<uint32_t, uint32_t>* alignmentTable = nullptr);
}
