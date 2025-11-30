#ifndef DINARA_CSTDINT_HPP
#define DINARA_CSTDINT_HPP

#include <cstdint>

// Sanity check that we are compiling on x86_64 or aarch64
#if !__x86_64__ && !__aarch64__
#error "Dinara can only be built on an x86_64 machine (64-bit Intel/AMD) or an ARM64 machine. "
#endif


namespace dinara {

    using std::uint8_t;
    using std::uint16_t;
    using std::uint32_t;
    using std::uint64_t;

    using std::int8_t;
    using std::int16_t;
    using std::int32_t;
    using std::int64_t;

}

#endif
