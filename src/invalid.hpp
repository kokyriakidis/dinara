#ifndef DINARA_INVALID_HPP
#define DINARA_INVALID_HPP

// In many contexts, we use invalid<T>
// to indicate a value that is invalid, uninitialized, or unknown.

#include <numeric>

namespace dinara {
    template<class T> static const T invalid = std::numeric_limits<T>::max();
    template<class T> static const T unlimited = std::numeric_limits<T>::max();
}

#endif
