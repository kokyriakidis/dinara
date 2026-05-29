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

// Detangle windows with multiple distinct through-flows.
//
// For each tangled window, identifies flows (triplets prev->current->next),
// creates bypass edges connecting prev directly to next, and removes
// flow reads from the current window's backbone anchors.
//
// Returns the number of windows detangled.
// bypassEdges: edges to add to the anchor graph.
uint64_t detangleWindows(
    Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const std::vector<AnchorWindow>& anchorWindows,
    uint64_t minFlowCoverage,
    std::vector<DetangleBypassEdge>& bypassEdges);

} // namespace dinara
