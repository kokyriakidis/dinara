#pragma once

#include "cstdint.hpp"
#include <climits>
#include <limits>

namespace dinara {

// =============================================================================
// Hifiasm lchain_dp scoring functions (Hash_Table.cpp)
// Used for read-to-read overlap detection
// =============================================================================

// Weight normalization macro (Hash_Table.cpp:20)
// Returns max(1, x/y) if x >= y, otherwise returns 1
#define HIFIASM_NORMAL_W(x, y) ((x) >= (y) ? (x) / (y) : 1)

// Dynamic bandwidth calculation (Hash_Table.cpp:1475-1488)
// Computes bandwidth based on alignment boundaries and bw_rate
static inline int32_t hifiasm_cal_bw(
    uint32_t posAi, uint32_t posBi,     // Current hit positions
    uint32_t posAj, uint32_t posBj,     // Previous hit positions
    double bw_rate,                      // Bandwidth rate (maxDriftRate)
    uint64_t readLenA, uint64_t readLenB // Read lengths
) {
    // Calculate alignment span on both reads
    int64_t sf_s = (int64_t)posAj;
    int64_t sf_e = (int64_t)posAi + 1;
    int64_t ot_s = (int64_t)posBj;
    int64_t ot_e = (int64_t)posBi + 1;

    // Remaining bases after current position
    int64_t sf_r = (int64_t)readLenA - sf_e;
    int64_t ot_r = (int64_t)readLenB - ot_e;

    // Adjust start position
    if (sf_s <= ot_s) {
        sf_s = 0;
    } else {
        sf_s -= ot_s;
    }

    // Adjust end position
    if (sf_r <= ot_r) {
        sf_e = (int64_t)readLenA;
    } else {
        sf_e += ot_r;
    }

    return (int32_t)((sf_e - sf_s) * bw_rate);
}

// Hifiasm's lchain_dp scoring function (Hash_Table.cpp:1490-1513)
// This is the CORRECT function for read-to-read overlaps
static inline int32_t hifiasm_comput_sc_ch(
    uint32_t posAi, uint32_t posBi,     // Current hit positions
    uint32_t posAj, uint32_t posBj,     // Previous hit positions
    uint32_t weightI, uint8_t spanI,    // CURRENT hit's weight and span (not previous!)
    double bw_rate,                      // Bandwidth rate
    double chn_pen_gap,                  // Gap penalty coefficient
    double chn_pen_skip,                 // Skip penalty coefficient
    uint64_t readLenA, uint64_t readLenB // Read lengths
) {
    // Query and target distances
    int32_t dq = (int32_t)posAi - (int32_t)posAj;
    if (dq <= 0) return INT32_MIN;

    int32_t dr = (int32_t)posBi - (int32_t)posBj;
    if (dr <= 0) return INT32_MIN;

    // Gap size (indel) and minimum distance
    int32_t dd = (dr > dq) ? (dr - dq) : (dq - dr);  // Indel size
    int32_t dg = (dr < dq) ? dr : dq;                 // Min distance

    // Dynamic bandwidth check (only if dd > 16)
    if (dd > 16) {
        int32_t bw = hifiasm_cal_bw(posAi, posBi, posAj, posBj, bw_rate, readLenA, readLenB);
        if (dd > bw) return INT32_MIN;
    }

    // Base score: min(span, dg) normalized by weight
    int32_t q_span = (int32_t)spanI;
    int32_t sc = (q_span < dg) ? q_span : dg;
    const int32_t w =
        (weightI > uint32_t(std::numeric_limits<int32_t>::max())) ?
        std::numeric_limits<int32_t>::max() : int32_t(weightI);
    sc = HIFIASM_NORMAL_W(sc, w);

    // Apply gap penalty if there are indels or if dg exceeds span
    if (dd > 0 || (dg > q_span && dg > 0)) {
        double lin_pen = chn_pen_gap * (double)dd;
        double a_pen = ((double)sc) * (((double)dd) / ((double)dg)) / bw_rate;

        // Take minimum of linear and adaptive penalties
        if (lin_pen > a_pen) {
            lin_pen = a_pen;
        }

        // Add skip penalty proportional to distance
        lin_pen += chn_pen_skip * (double)dg;

        sc -= (int32_t)lin_pen;
    }

    return sc;
}

// Hifiasm EC chaining score (Hash_Table.cpp:1515-1541, comput_sc_ch_ec).
// Differs from comput_sc_ch by using a different long-gap penalty behavior:
//   - dd < 4  => use min(linear, adaptive)
//   - dd >= 4 => use max(linear, adaptive)
static inline int32_t hifiasm_comput_sc_ch_ec(
    uint32_t posAi, uint32_t posBi,
    uint32_t posAj, uint32_t posBj,
    uint32_t weightI, uint8_t spanI,
    double bw_rate,
    double chn_pen_gap,
    double chn_pen_skip,
    uint64_t readLenA, uint64_t readLenB
) {
    int32_t dq = (int32_t)posAi - (int32_t)posAj;
    if (dq <= 0) return INT32_MIN;

    int32_t dr = (int32_t)posBi - (int32_t)posBj;
    if (dr <= 0) return INT32_MIN;

    int32_t dd = (dr > dq) ? (dr - dq) : (dq - dr);
    int32_t dg = (dr < dq) ? dr : dq;

    if (dd > 16) {
        int32_t bw = hifiasm_cal_bw(posAi, posBi, posAj, posBj, bw_rate, readLenA, readLenB);
        if (dd > bw) return INT32_MIN;
    }

    int32_t q_span = (int32_t)spanI;
    int32_t sc = (q_span < dg) ? q_span : dg;
    const int32_t w =
        (weightI > uint32_t(std::numeric_limits<int32_t>::max())) ?
        std::numeric_limits<int32_t>::max() : int32_t(weightI);
    sc = HIFIASM_NORMAL_W(sc, w);

    if (dd > 0 || (dg > q_span && dg > 0)) {
        double lin_pen = chn_pen_gap * (double)dd;
        double a_pen = ((double)sc) * (((double)dd) / ((double)dg)) / bw_rate;
        if (dd < 4) {
            lin_pen = (lin_pen > a_pen) ? a_pen : lin_pen;
        } else {
            lin_pen = (lin_pen < a_pen) ? a_pen : lin_pen;
        }
        lin_pen += chn_pen_skip * (double)dg;
        sc -= (int32_t)lin_pen;
    }

    return sc;
}

// Hifiasm DP chaining parameters (from inter.h)
struct HifiasmDPParams {
    int32_t max_iter = 10000;        // Maximum lookback window (hifiasm default)
    int32_t max_skip = 25;           // Maximum consecutive skips before breaking
    int32_t bw = 500;                // Bandwidth for diagonal constraint
    float chn_pen_gap = 0.1f;        // Gap penalty coefficient (W_CHN_PEN_GAP in inter.h:28)
    int32_t max_dist_x = 50000;      // Max distance on query axis
    int32_t max_dist_y = 50000;      // Max distance on reference axis
    int32_t min_cnt = 2;             // Minimum chain length
    int32_t min_score = 30;          // Minimum chain score
};

} // namespace dinara
