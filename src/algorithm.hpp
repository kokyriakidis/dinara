#ifndef DINARA_ALGORITHM_HPP
#define DINARA_ALGORITHM_HPP

#include <algorithm>

// With C++20, it is best to use the algorithms in the std::ranges namespace.
// For more information see:
// https://en.cppreference.com/w/cpp/algorithm/ranges

namespace dinara {
    using std::ranges::copy;
    using std::fill;
    using std::find;
    using std::max;
    using std::min;
    using std::sort;
    using std::swap;
    using std::unique;
}

#endif

