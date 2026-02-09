#ifndef DINARA_HIFIASM_BOUNDARY_VERIFY_HPP
#define DINARA_HIFIASM_BOUNDARY_VERIFY_HPP

#include "ReadId.hpp"

#include "cstdint.hpp"
#include "vector.hpp"

namespace dinara {
    class Reads;

    // Hifiasm-parity boundary verification used by chimeric-read detection.
    // This mirrors hifiasm's boundary_verify/verify_single_window stack and
    // returns true if the projected boundary interval is sequence-consistent.
    bool hifiasmBoundaryVerify(
        const Reads& reads,
        uint32_t qIntervalStart,
        uint32_t qIntervalEnd,
        ReadId qId,
        ReadId tId,
        uint32_t qs,
        uint32_t ts,
        uint32_t te,
        bool rev,
        vector<char>& xBuffer,
        vector<char>& yBuffer);
}

#endif

