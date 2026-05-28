#pragma once

#include "AnchorWindows.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"

#include <map>
#include <vector>

namespace dinara {

// Detangle windows with multiple distinct through-flows.
// Splits backbone anchors so each path gets its own copies.
//
// Returns the number of windows detangled.
// anchorSplitMap: originalAnchorId -> [newAnchorId per path].
uint64_t detangleWindows(
    Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const std::vector<AnchorWindow>& anchorWindows,
    uint64_t minFlowCoverage,
    std::map<Shasta2AnchorId, std::vector<Shasta2AnchorId>>& anchorSplitMap);

} // namespace dinara
