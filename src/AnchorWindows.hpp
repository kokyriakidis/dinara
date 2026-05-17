#ifndef DINARA_ANCHOR_WINDOWS_HPP
#define DINARA_ANCHOR_WINDOWS_HPP

#include "ReadId.hpp"
#include "Shasta2Anchors.hpp"
#include "cstdint.hpp"

#include <vector>

namespace dinara {

// A read's interval within an anchor window, expressed as journey positions.
struct AnchorWindowReadInterval {
    OrientedReadId orientedReadId;
    uint32_t begin = 0; // Inclusive position in the oriented read journey.
    uint32_t end = 0;   // Exclusive position in the oriented read journey.
    uint32_t touchedAnchorCount = 0; // Anchors shared with the backbone in this interval.
};

// An anchor window: a contiguous interval on a backbone read's journey,
// plus the intervals on all other reads that share anchors with the backbone.
struct AnchorWindow {
    uint32_t windowId = 0;
    OrientedReadId backboneOrientedReadId;
    uint32_t backboneBegin = 0; // Inclusive journey position on backbone.
    uint32_t backboneEnd = 0;   // Exclusive journey position on backbone.
    uint32_t claimedAnchorCount = 0;
    std::vector<AnchorWindowReadInterval> readIntervals;
};

} // namespace dinara

#endif // DINARA_ANCHOR_WINDOWS_HPP
