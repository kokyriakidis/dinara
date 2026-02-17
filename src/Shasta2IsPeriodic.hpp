#pragma once

#include "cstdint.hpp"

namespace dinara {
    template<class RandomIterator> inline bool shasta2IsPeriodic(RandomIterator begin, RandomIterator end, uint64_t period);
}

// Return true if the sequence defined by the given iterators is periodic with period.
template<class RandomIterator> inline bool dinara::shasta2IsPeriodic(
    RandomIterator begin,
    RandomIterator end,
    const uint64_t period)
{
    if((end - begin) % period) {
        return false;
    }

    for(RandomIterator itA=begin; ; ++itA) {
        const RandomIterator itB = itA + period;
        if(itB >= end) {
            break;
        }
        if(*itA != *itB) {
            return false;
        }
    }
    return true;
}
