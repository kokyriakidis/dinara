#pragma once

// Compute per-window transition counts directly from journeys,
// without building the anchor graph. Populates each AnchorWindow's
// transitionReads map and per-read previousWindow/nextWindow fields.
//
// This decouples transition counting from graph construction so
// detangling can run before the first graph build.

#include "AnchorWindows.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"

#include <vector>

namespace dinara {

// Walk all journeys, map anchors to windows via backbone positions,
// and populate:
//   - AnchorWindow::transitionReads[(prev, next)] -> read list
//   - AnchorWindowReadInterval::previousWindow / nextWindow
//   - AnchorWindow::backbonePreviousWindow / backboneNextWindow
// anchorDovetailWindow (optional): forward-oriented anchorId -> owning windowId
// for claimed dovetail anchors. When provided and non-empty, dovetails are
// added to anchorToWindow (additive, backbone stays authoritative) so they are
// treated as part of their window for transition counting.
void computeWindowTransitions(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    std::vector<AnchorWindow>& anchorWindows,
    const std::vector<uint32_t>* anchorDovetailWindow = nullptr);

} // namespace dinara
