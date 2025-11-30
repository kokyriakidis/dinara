#ifndef DINARA_MARKER_INTERVAL_HPP
#define DINARA_MARKER_INTERVAL_HPP

// Dinara.
#include "ReadId.hpp"

// Standard library.
#include "array.hpp"
#include "tuple.hpp"
#include "vector.hpp"

namespace dinara {
    class MarkerInterval;
    class MarkerIntervalWithRepeatCounts;
}


// Class to describe the interval between
// two markers on an oriented read.
// The two markers are not necessarily consecutive.
// However, the second marker has a higher ordinal
// than the first.
class dinara::MarkerInterval {
public:
    OrientedReadId orientedReadId;

    // The ordinals of the two markers.
    array<uint32_t, 2> ordinals;

    MarkerInterval() {}
    MarkerInterval(
        OrientedReadId orientedReadId,
        uint32_t ordinal0,
        uint32_t ordinal1) :
        orientedReadId(orientedReadId)
    {
        ordinals[0] = ordinal0;
        ordinals[1] = ordinal1;
    }

    bool operator==(const MarkerInterval& that) const
    {
        return
            tie(orientedReadId, ordinals[0], ordinals[1])
            ==
            tie(that.orientedReadId, that.ordinals[0], that.ordinals[1]);
    }
    bool operator<(const MarkerInterval& that) const
    {
        return
            tie(orientedReadId, ordinals[0], ordinals[1])
            <
            tie(that.orientedReadId, that.ordinals[0], that.ordinals[1]);
    }
};



class dinara::MarkerIntervalWithRepeatCounts :
    public MarkerInterval {
public:
    vector<uint8_t> repeatCounts;

    // The constructor does not fill in the repeat counts.
    MarkerIntervalWithRepeatCounts(const MarkerInterval& markerInterval) :
        MarkerInterval(markerInterval){}
};

#endif
