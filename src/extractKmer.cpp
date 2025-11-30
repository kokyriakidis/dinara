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
    const Int zMask = Int(((1UL << n) - 1UL) << zShift);
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

    // Access the first two words containing the k-mer we want.
    const uint64_t i0 = (position >> 6) << 1;
    const uint64_t i1 = i0 + 1;
    array<uint64_t, 2> ww01 = {v.begin[i0], v.begin[i1]};

    // The starting position of the k-mer in the first two words.
    const uint64_t position01 = position & 63;
    // cout << "position01 " << position01 << endl;

    // The number of k-mer bases in the first two words.
    const uint64_t length01 = min(length, 64 - position01);

    // Store these length01 bases at the beginning of the ShortBaseSequence.
    s.data[0] = 0;
    s.data[1] = 0;
    extractBits(&(ww01[0]), position01, length01, &(s.data[0]), 0);

    // If the k-mer is entirely contained in the first two words, we are done.
    if(length01 == length) {
        return;
    }

    // Get the remaining bases from the next two words.
    // The number of k-mer bases in the second two words
    const uint64_t length23 = length - length01;

    // Access the second two words.
    const uint64_t i2 = i1 + 1;
    const uint64_t i3 = i2 + 1;
    array<uint64_t, 2> ww23 = {v.begin[i2], v.begin[i3]};

    // Store the most significant length23 bits in the ShortBaseSequence,
    // following the ones we already stored.
    extractBits(&(ww23[0]), 0, length23, &(s.data[0]), length01);

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

