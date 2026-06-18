#pragma once

// Read-only diagnostic for the shared dovetail-halo model.
//
// Each window's "halo" is the set of forward-oriented anchors from its backbone
// plus the whole journeys of all reads touching the backbone (produced by
// computeAnchorWindowsClean). Halos are multi-owner: a read connecting two
// neighboring windows contributes its journey to BOTH halos, so the
// intersection of two windows' halos is the overlap that should register their
// connection.
//
// This diagnostic measures those overlaps without modifying any state, to
// validate the halo premise before building overlap-based connectivity:
//   - how many window pairs share halo anchors,
//   - the distribution of overlap sizes,
//   - per-window: how many other windows it overlaps (specific vs promiscuous).

#include <cstdint>
#include <vector>

namespace dinara {

// windowHalos[windowId] = sorted/unsorted forward anchor ids in that halo.
// anchorCount is the total number of anchors (used to size the inverted index).
void reportWindowHaloOverlaps(
    const std::vector<std::vector<uint32_t>>& windowHalos,
    uint64_t anchorCount);

} // namespace dinara
