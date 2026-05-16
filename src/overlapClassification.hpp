#ifndef DINARA_OVERLAP_CLASSIFICATION_HPP
#define DINARA_OVERLAP_CLASSIFICATION_HPP

#include "cstdint.hpp"

namespace dinara {

// Return value constants — identical to hifiasm (Overlaps.h).
#define MA_HT_INT        (-1)
#define MA_HT_QCONT      (-2)
#define MA_HT_TCONT      (-3)
#define MA_HT_SHORT_OVLP (-4)

/*
Exact port of hifiasm's ma_hit2arc (Overlaps.h).

Classifies an overlap between a query read and a target read.
In hifiasm the function also populates an asg_arc_t for dovetail overlaps;
we only need the classification, so the arc output is omitted.

Parameters (matching hifiasm):
  qs, qe  — query start/end on forward strand
  ql      — query length
  ts, te  — target start/end on forward strand
  tl      — target length
  rev     — true if target is reverse-complemented
  max_hang — maximum allowed extension (default 1000)
  int_frac — minimum overlap fraction (default 0.8f)
  min_ovlp — minimum effective overlap length (default 50)

Return values (matching hifiasm):
  >= 0            dovetail overlap (value is the non-overlap length l)
  MA_HT_INT       (-1)  internal match (excessive overhangs or low overlap fraction)
  MA_HT_QCONT     (-2)  query contained in target
  MA_HT_TCONT     (-3)  target contained in query
  MA_HT_SHORT_OVLP (-4)  overlap too short

///in default, max_hang = 1000, int_frac = 0.8, min_ovlp = 50
*/
inline int ma_hit2arc(
    int32_t qs, int32_t qe, int32_t ql,
    int32_t ts, int32_t te, int32_t tl,
    bool rev,
    int32_t max_hang,
    float int_frac,
    int32_t min_ovlp)
{
    int32_t tl5, tl3, ext5, ext3;
    uint32_t l;

    ///if query and target are in different strand
    if (rev) tl5 = tl - te, tl3 = ts; // tl5: 5'-end overhang (on the query strand); tl3: similar
    else tl5 = ts, tl3 = tl - te;

    ///ext5 and ext3 is the hang on left side and right side, respectively
    ext5 = qs < tl5? qs : tl5;
    ext3 = ql - qe < tl3? ql - qe : tl3;

    ///ext3 and ext5 should be always 0
    if (ext5 > max_hang || ext3 > max_hang
    || qe - qs < (qe - qs + ext5 + ext3) * int_frac
    || te - ts < (te - ts + ext5 + ext3) * int_frac)
    {
        return MA_HT_INT;
    }

    if (qs <= tl5 && ql - qe <= tl3) return MA_HT_QCONT; // query contained in target
    else if (qs >= tl5 && ql - qe >= tl3) return MA_HT_TCONT; // target contained in query
    else if (qs > tl5) l = qs - tl5; ///query-to-target overlap, l is the length of node in string graph
    else l = (ql - qe) - tl3; ///target-to-query overlap, l is the length of node in string graph
    if (qe - qs + ext5 + ext3 < min_ovlp || te - ts + ext5 + ext3 < min_ovlp) return MA_HT_SHORT_OVLP; // short overlap

    return (int)l;
}

} // namespace dinara

#endif // DINARA_OVERLAP_CLASSIFICATION_HPP
