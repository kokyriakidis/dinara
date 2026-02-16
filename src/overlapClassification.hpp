#ifndef DINARA_OVERLAP_CLASSIFICATION_HPP
#define DINARA_OVERLAP_CLASSIFICATION_HPP

#include "cstdint.hpp"

namespace dinara {

/*
Hifiasm ma_hit2arc: classify overlap type based on overhang analysis.

This is a direct translation of hifiasm's `ma_hit2arc` function used for both
ma_hit_flt (hanging overlap filtering) and ma_hit_contained_advance (containment detection).

Given an overlap with:
  - Query: interval [qs, qe) on read of length ql
  - Target: interval [ts, te) on read of length tl
  - isReverse: whether target is reverse-complemented

Algorithm:
1. Compute overhangs (unaligned portions):
   - Query 5' overhang: qs
   - Query 3' overhang: ql - qe
   - Target overhangs: adjusted for orientation (forward vs RC)

2. Compute extensions (minimum overhangs at each end):
   - ext5 = min(query 5' overhang, target 5' overhang)
   - ext3 = min(query 3' overhang, target 3' overhang)

3. Apply filters in order:
   a) Absolute overhang filter: ext5 > max_hang OR ext3 > max_hang -> internal (-1)
   b) Fractional overlap filter: overlap/(overlap+ext5+ext3) < int_frac -> internal (-1)
   c) Containment detection: compare query/target overhangs -> contained (1 or 2)
   d) Minimum overlap filter: effective_length < min_ovlp -> too short (-2)
   e) Default: dovetail (0)

Return values:
   0 = dovetail (proper prefix-suffix overlap)
   1 = query contained in target (QCONT)
   2 = target contained in query (TCONT)
  -1 = internal match (excessive overhangs, reject)
  -2 = too short (effective overlap < min_ovlp, reject)

Hifiasm parameters (typical defaults):
  - max_hang = 1000 (asm_opt.max_hang)
  - int_frac = 0.8 (asm_opt.int_frac)
  - min_ovlp = 50 (asm_opt.min_ovlp)

Coordinates:
  - qs/qe are on the query read's forward strand.
  - ts/te are on the target read's forward strand (even for RC overlaps).
  - isSameStrand indicates whether the alignment is forward/forward or forward/reverse.
  - The function applies an orientation-aware conversion to determine 5'/3' overhangs.
*/
inline int ma_hit2arc_containment(
    int32_t qs, int32_t qe, int32_t ql,
    int32_t ts, int32_t te, int32_t tl,
    bool isReverse,
    int32_t max_hang,
    double int_frac,
    int32_t min_ovlp)
{
    // Step 1: Compute target overhangs in query orientation
    int32_t tl5, tl3;  // Target 5' and 3' overhangs
    if (isReverse) {
        // Reverse-complement: 5' and 3' are swapped
        tl5 = tl - te;  // Target 5' overhang = bases after target end
        tl3 = ts;       // Target 3' overhang = bases before target start
    } else {
        // Forward orientation
        tl5 = ts;       // Target 5' overhang = bases before target start
        tl3 = tl - te;  // Target 3' overhang = bases after target end
    }

    // Step 2: Compute extensions (minimum overhangs at each end)
    int32_t ext5 = (qs < tl5) ? qs : tl5;           // 5' extension
    int32_t ext3 = ((ql - qe) < tl3) ? (ql - qe) : tl3;  // 3' extension

    // Step 3a: Hifiasm ma_hit_flt: reject if absolute overhangs exceed threshold
    if (ext5 > max_hang || ext3 > max_hang) {
        return -1;  // Internal match
    }

    // Step 3b: Hifiasm ma_hit_flt: reject if fractional overlap is too low
    int32_t qOverlapLen = qe - qs;
    int32_t tOverlapLen = te - ts;

    // Check query: overlap_len / (overlap_len + ext5 + ext3) >= int_frac
    // Equivalent: overlap_len >= (overlap_len + ext5 + ext3) * int_frac
    if (qOverlapLen < (qOverlapLen + ext5 + ext3) * int_frac) {
        return -1;  // Internal match
    }
    if (tOverlapLen < (tOverlapLen + ext5 + ext3) * int_frac) {
        return -1;  // Internal match
    }

    // Step 3c: Hifiasm containment detection
    // QCONT: query is contained if both query overhangs <= corresponding target overhangs
    if (qs <= tl5 && (ql - qe) <= tl3) {
        return 1;  // Query contained in target
    }
    // TCONT: target is contained if both target overhangs <= corresponding query overhangs
    if (qs >= tl5 && (ql - qe) >= tl3) {
        return 2;  // Target contained in query
    }

    // Step 3d: Hifiasm minimum overlap filter
    // Effective overlap length = overlap + extensions
    if (qOverlapLen + ext5 + ext3 < min_ovlp || tOverlapLen + ext5 + ext3 < min_ovlp) {
        return -2;  // Too short
    }

    // Step 3e: Default case - proper dovetail overlap
    return 0;
}

} // namespace dinara

#endif // DINARA_OVERLAP_CLASSIFICATION_HPP
