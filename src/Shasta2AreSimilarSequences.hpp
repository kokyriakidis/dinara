#pragma once

#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "vector.hpp"

namespace dinara {
    class Base;

    bool shasta2AreSimilarSequences(
        const vector<Base>& x,
        const vector<Base>& y,
        const vector<uint64_t>& minRepeatCount,
        ostream& html);
}
