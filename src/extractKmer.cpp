// Extract the kmer at a given position in a LongBaseSequenceView.
#include "extractKmer.hpp"
#include "LongBaseSequence.hpp"
#include "DINARA_ASSERT.hpp"
#include "ShortBaseSequence.hpp"
using namespace dinara;



// Extract n bits from x, starting at position xPosition,
// and store them in y, starting at position yPosition,
// leaving the remaining bits of y unchanged.
// Here, xPosition and yPosition are counted with
// 0 at the most significant bit and moving towards the least
// significant bit.
// THIS IS DONE SEPARATELY ON x[0], y[0] AND x[1], y[1].
template<class Int> inline void dinara::extractBits(
    const uint64_t* x,
    uint64_t xPosition, // 0 = MSB
    uint64_t n,
    Int* y,
    uint64_t yPosition  // 0 = MSB
    )
{
    using std::cout;
    using std::endl;

    uint64_t x0 = x[0];
    uint64_t x1 = x[1];
    Int& y0 = y[0];
    Int& y1 = y[1];


    DINARA_ASSERT(xPosition + n <= 64);

    // Shift x right so the n bits are the least significant.
    const uint64_t xShift = 64 - xPosition - n;
    x0 >>= xShift;
    x1 >>= xShift;

    // Copy to an Int.
    Int z0 = Int(x0);
    Int z1 = Int(x1);

    // Shift left so the n bits are in the right place.
    const uint64_t zShift = 8 * sizeof(Int) - yPosition - n;
    z0 <<= zShift;
    z1 <<= zShift;

    // Copy these n bits to y without changing the remaining bits.
    // Build the mask in a 128-bit intermediate: n is bounded by 64 (a single
    // read word), but zShift can reach 8*sizeof(Int)-1, which is up to 127 for
    // a capacity-128 Kmer. Doing the shift in 64-bit space (the previous
    // ((1UL<<n)-1UL)<<zShift) is undefined for zShift>=64 -- this is why no
    // ShortBaseSequence<__uint128_t> instantiation existed.
    const __uint128_t maskBits =
        ((__uint128_t(1) << n) - __uint128_t(1)) << zShift;
    const Int zMask = Int(maskBits);
    const Int yMask = ~zMask;
    y0 = (y0 & yMask) | (z0 & zMask);
    y1 = (y1 & yMask) | (z1 & zMask);

}



template<class Int> void dinara::extractKmer(
    const LongBaseSequenceView& v,
    uint64_t position,
    uint64_t length,
    ShortBaseSequence<Int>& s)
{

    // Sanity checks.
    DINARA_ASSERT(length <= s.capacity);
    DINARA_ASSERT(position + length <= v.baseCount);

    s.data[0] = 0;
    s.data[1] = 0;

    // Bases are stored two words per 64-base block. A k-mer can straddle
    // several blocks: up to 2 for a capacity-64 Kmer, but up to 3 for a
    // capacity-128 Kmer (e.g. k=70 starting at intra-block offset 63 spans
    // 1 + 64 + 5 bases). Iterate block by block, copying at most 64 bases per
    // step so each extractBits call satisfies xPosition + n <= 64.
    uint64_t remaining = length;
    uint64_t srcPosition = position;    // absolute base position in the read
    uint64_t dstPosition = 0;           // base position within the ShortBaseSequence
    while(remaining > 0) {
        const uint64_t blockWord0 = (srcPosition >> 6) << 1;
        array<uint64_t, 2> ww = {v.begin[blockWord0], v.begin[blockWord0 + 1]};

        const uint64_t offsetInBlock = srcPosition & 63;          // 0..63, MSB-based
        const uint64_t take = min(remaining, 64 - offsetInBlock); // <= 64 - offset

        extractBits(&(ww[0]), offsetInBlock, take, &(s.data[0]), dstPosition);

        srcPosition += take;
        dstPosition += take;
        remaining   -= take;
    }
}



// Explicit instantiations.
template void dinara::extractKmer(
    const LongBaseSequenceView&,
    uint64_t position,
    uint64_t length,
    ShortBaseSequence<uint8_t>&);

template void dinara::extractKmer(
    const LongBaseSequenceView&,
    uint64_t position,
    uint64_t length,
    ShortBaseSequence<uint16_t>&);

template void dinara::extractKmer(
    const LongBaseSequenceView&,
    uint64_t position,
    uint64_t length,
    ShortBaseSequence<uint32_t>&);

template void dinara::extractKmer(
    const LongBaseSequenceView&,
    uint64_t position,
    uint64_t length,
    ShortBaseSequence<uint64_t>&);

// Capacity-128 Kmer (Kmer128), used only when DINARA_LONG_MARKERS is enabled.
// Guarded because ShortBaseSequence<__uint128_t>::id() requires
// BitCounter<__uint128_t>::doubleSizeType, which is defined only in that build.
#ifdef DINARA_LONG_MARKERS
template void dinara::extractKmer(
    const LongBaseSequenceView&,
    uint64_t position,
    uint64_t length,
    ShortBaseSequence<__uint128_t>&);
#endif

