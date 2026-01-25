#ifndef DINARA_HIFIASM_COORDINATE_TRANSFORMS_HPP
#define DINARA_HIFIASM_COORDINATE_TRANSFORMS_HPP

#include "DINARA_ASSERT.hpp"
#include "cstdint.hpp"
#include "utility.hpp"

namespace dinara {

// Convert a half-open interval [s,e) on the reverse-complement oriented read
// (coordinates on the RC sequence, 0 <= s <= e <= len) into forward-strand
// coordinates on the underlying read.
// For half-open intervals this mapping is:
//   [s,e)  ->  [len - e, len - s)
inline std::pair<uint32_t, uint32_t> rcIntervalToForward(uint32_t len, uint32_t s, uint32_t e)
{
    DINARA_ASSERT(s <= e);
    DINARA_ASSERT(e <= len);
    return {len - e, len - s};
}

} // namespace dinara

#endif

