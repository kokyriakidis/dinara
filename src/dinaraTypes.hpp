#pragma once

#ifdef DINARA_LONG_MARKERS
#include "Uint256.hpp"
#endif

#include <cstdint.hpp>

namespace dinara {

    // A KmerId requires twice the number of bits of the integer type
    // used to define ShortBaseSequence.
    using KmerId16 = uint32_t;
    using KmerId32 = uint64_t;
    using KmerId64 = __uint128_t;

#ifdef DINARA_LONG_MARKERS
    // With capacity-128 Kmers, KmerId needs 256 bits. Use a flat POD type
    // (four uint64_t words) rather than boost::multiprecision::uint256_t:
    // KmerId is memory-mapped, hashed by raw bytes, and used as a hash-table
    // key, all of which require a trivially copyable, flat byte layout.
    using KmerId128 = Uint256;
    using KmerId = KmerId128;
#else
    using KmerId = KmerId64;
#endif

    using ReadId = uint32_t;
    using Strand = ReadId;

    using MarkerId = uint64_t;
    using MarkerGraphVertexId = uint64_t;
    using MarkerGraphEdgeId = uint64_t;

    using AssemblyGraphVertexId = uint64_t;
    using AssemblyGraphEdgeId = uint64_t;
}

