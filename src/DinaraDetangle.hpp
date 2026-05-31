#pragma once

#include "AnchorWindows.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"

#include <set>
#include <vector>

namespace dinara {

// A bypass edge created by detangling: connects the last anchor of a
// previous window directly to the first anchor of a next window,
// bypassing the tangled window in between.
struct DetangleBypassEdge {
    Shasta2AnchorId anchorIdA;  // Last anchor in prev window.
    Shasta2AnchorId anchorIdB;  // First anchor in next window.
};

// Detangle windows by bypassing all through-flows.
//
// For each window X with reads flowing A -> X -> B, creates bypass
// edges connecting A directly to B and removes the flow reads from
// X's backbone anchors. Every through-flow is bypassed regardless
// of read count.
//
// Returns the number of windows detangled.
// bypassEdges: edges to add to the anchor graph.
uint64_t detangleWindows(
    Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const std::vector<AnchorWindow>& anchorWindows,
    std::vector<DetangleBypassEdge>& bypassEdges);

} // namespace dinara
