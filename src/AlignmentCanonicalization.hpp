#pragma once

#include "Alignment.hpp"

namespace dinara {

// Canonicalize a candidate (readId0/readId1 + orientation) so that readId0 < readId1.
// The associated Alignment is transformed to remain consistent with the candidate:
// - The alignment always refers to (readId0 on strand 0, readId1 on strand 0 if same-strand,
//   or strand 1 if diff-strand).
// - Coordinates (qs/qe on readId0, ts/te on readId1) are stored in forward-read coordinates
//   for each read (canonical Dinara convention).
inline void canonicalizeCandidateAndAlignment(
    ReadId& readId0,
    ReadId& readId1,
    bool& isSameStrand,
    Alignment& alignment,
    uint32_t markerCount0,
    uint32_t markerCount1)
{
    // Already canonical.
    if(readId0 < readId1) {
        return;
    }

    DINARA_ASSERT(readId0 != readId1);
    std::swap(readId0, readId1);

    // The alignment is currently for (old readId0 strand0, old readId1 strand{0/1}).
    // We need it for (new readId0 strand0, new readId1 strand{0/1}).
    if(isSameStrand) {
        // (A0,B0) -> (B0,A0)
        alignment.swap();
    } else {
        // (A0,B1) -> (B0,A1):
        // reverse-complement both reads to turn (A0,B1) into (A1,B0),
        // then swap the two columns.
        alignment.reverseComplement(markerCount0, markerCount1);
        alignment.swap();
    }

    // Coordinates are always stored in forward-read coordinates.
    // Swapping query/target just swaps their intervals.
    std::swap(alignment.qs, alignment.ts);
    std::swap(alignment.qe, alignment.te);
}

} // namespace dinara

