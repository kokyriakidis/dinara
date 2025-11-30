#ifndef DINARA_CHRONO_HPP
#define DINARA_CHRONO_HPP

#include <chrono>

/*******************************************************************************

Usage pattern (from code in dinara namespace):

const auto t0 = steady_clock::now();
const auto t1 = steady_clock::now();
const double t01 = seconds(t1-t0);   // Can use auto instead of double.

*******************************************************************************/

namespace dinara {
    using std::chrono::steady_clock;

    template<class Duration> double seconds(Duration duration)
    {
        return std::chrono::duration<double>(duration).count();
    }
}

#endif
