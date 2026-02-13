#pragma once

#include "mode3-Anchor.hpp"

namespace dinara {
namespace mode3 {

// DirectedAnchors is a specialized version of Anchors that stores
// base positions directly in the journey data structure.
// This is used to optimize DirectedAnchorGraph construction by
// avoiding slow linear searches for anchor locations in reads.
// ============================================================================
// DirectedAnchors — An optimized version of mode3::Anchors.
//
// Unlike the legacy Anchors class which requires searching for occurrences,
// DirectedAnchors precomputes "journeys" that store base-pair positions
// directly. This allows the DirectedAnchorGraph to be built in near-linear
// time without slow coordinate lookups.
// ============================================================================
class DirectedAnchors : public Anchors {
public:
    // We inherit mostly everything from Anchors to keep compatibility,
    // but we use a different journey structure.
    // ============================================================================
// JourneyAnchor — stores an anchor ID plus its base-pair coordinates on a read.
// Using explicit start/end positions allows MBG-style O(1) overlap calculation.
// ============================================================================
    struct JourneyAnchor {
        AnchorId anchorId;
        uint32_t start; // Base position in oriented read
        uint32_t end;   // End position in oriented read
        
        // Comparison for sorting anchors within a journey by their read position.
        bool operator<(const JourneyAnchor& other) const {
            if(start != other.start) return start < other.start;
            return anchorId < other.anchorId;
        }
    };

    // Replace the default journeys with one that includes positions.
    MemoryMapped::VectorOfVectors<JourneyAnchor, uint64_t> journeysWithPositions;

    DirectedAnchors(
        const MappedMemoryOwner& mappedMemoryOwner,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph& markerGraph,
        uint64_t minPrimaryCoverage,
        uint64_t maxPrimaryCoverage,
        uint64_t threadCount,
        bool createFromVertices);

    void computeJourneysWithPositions(uint64_t threadCount);
};

} // namespace mode3
} // namespace dinara
