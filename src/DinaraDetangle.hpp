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

// Detangle windows using Verkko-style triplet resolution.
//
// For each window X with multiple predecessors and/or successors,
// checks if every significant edge is covered by a solid triplet.
// If so, creates bypass edges for each solid (prev, succ) pairing
// and removes the flow reads from X's backbone anchors.
// Linear windows (≤1 predecessor AND ≤1 successor) are never bypassed.
//
// Returns the number of windows detangled.
// bypassEdges: edges to add to the anchor graph.
uint64_t detangleWindows(
    Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const std::vector<AnchorWindow>& anchorWindows,
    std::vector<DetangleBypassEdge>& bypassEdges);

// Detangle windows using the G-test (log-likelihood ratio test).
//
// For each window X with multiple predecessors and/or successors,
// builds a tangle matrix from transitionReads and runs the G-test
// to find the best connectivity hypothesis. If the best hypothesis
// is injective (each entrance maps to at most one exit and vice versa),
// creates bypass edges for each connected (entrance, exit) pair
// and removes the flow reads from X's backbone anchors.
//
// epsilon: noise parameter for the G-test (default 0.1).
//          Higher values tolerate more off-diagonal noise.
//
// Returns the number of windows detangled.
uint64_t detangleWindowsGTest(
    Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const std::vector<AnchorWindow>& anchorWindows,
    std::vector<DetangleBypassEdge>& bypassEdges,
    double epsilon = 0.1);

} // namespace dinara
