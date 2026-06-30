#ifndef DINARA_MSA_HET_KERNEL_HPP
#define DINARA_MSA_HET_KERNEL_HPP

// Pure-matrix kernels for abPOA-MSA-based het-site detection and 2-base
// allele-anchor construction. These functions take an abPOA-style row-column
// MSA matrix and contain the correctness-critical coordinate logic, with NO
// dependency on abPOA, the Assembler, or any I/O — so they are unit-testable
// in isolation (see tests/test_msa_het_kernel.cpp).
//
// MSA matrix conventions (matching abPOA abc->msa_base):
//   - msa[row][col], row 0 = backbone (the spine), rows 1..n-1 = members.
//   - Base encoding: 0=A 1=C 2=G 3=T, 4=N. Gap is a caller-supplied value
//     (abPOA uses abpt->m, typically 5). Anything == gapValue is a gap;
//     anything >= 4 and != gapValue is treated as N (non-informative).
//
// THREE coordinate frames — never conflate:
//   - MSA column   c  in [0, msaLen): includes gap columns from member insertions.
//   - backbone base bb in [0, backboneLen): the spine, derived from row-0 non-gaps.
//   - read base position: index into the oriented read (what export needs).
//
// The window-local backbone base bb maps to a full-read position via a
// caller-supplied origin: fullReadPos = bbStartBase + bb.

#include <cstdint>
#include <vector>
#include <unordered_map>

namespace dinara {
namespace msahet {

constexpr int kMsaBaseN = 4; // encodings >= this (and != gap) are N / unknown.

// --------------------------------------------------------------------------
// Pass 1: column <-> backbone-base maps, by walking row 0 (the backbone).
// --------------------------------------------------------------------------
struct ColumnMaps {
    // colToBb[c] = backbone base index at column c, or -1 if c is a member-
    // insertion column (row 0 is a gap there).
    std::vector<int> colToBb;
    // bbToCol[bb] = the MSA column holding backbone base bb.
    std::vector<int> bbToCol;
    int backboneLen = 0; // number of non-gap entries in row 0.
};

// Build the column maps from row 0. msaRow0 has length msaLen.
inline ColumnMaps buildColumnMaps(const uint8_t* msaRow0, int msaLen, int gapValue) {
    ColumnMaps m;
    m.colToBb.assign(msaLen, -1);
    m.bbToCol.reserve(msaLen);
    int bb = 0;
    for (int c = 0; c < msaLen; ++c) {
        if (int(msaRow0[c]) != gapValue) {
            m.colToBb[c] = bb;
            m.bbToCol.push_back(c);
            ++bb;
        }
    }
    m.backboneLen = bb;
    return m;
}

// --------------------------------------------------------------------------
// Pass 2 helper: tally member alleles at one backbone column.
// --------------------------------------------------------------------------
struct ColumnTally {
    uint32_t count[4] = {0, 0, 0, 0};
    uint32_t fwd[4]   = {0, 0, 0, 0};
    uint32_t rev[4]   = {0, 0, 0, 0};
    uint32_t covered  = 0; // member rows with a real base (0..3) here.

    // Most frequent base != refBase. Returns 0xFF if none.
    uint8_t altBase(uint8_t refBase) const {
        int best = -1;
        uint32_t bestN = 0;
        for (int b = 0; b < 4; ++b) {
            if (uint8_t(b) == refBase) continue;
            if (count[b] > bestN) { bestN = count[b]; best = b; }
        }
        return best < 0 ? uint8_t(0xFF) : uint8_t(best);
    }
};

// Tally member rows (1..nSeq-1) at MSA column c. strandIsRev[row-1] gives the
// member's strand (true = reverse) for strand-bias accounting.
inline ColumnTally tallyColumn(
    const uint8_t* const* msa, int nSeq, int c, int gapValue,
    const std::vector<bool>& strandIsRev)
{
    ColumnTally t;
    for (int row = 1; row < nSeq; ++row) {
        const uint8_t v = msa[row][c];
        if (int(v) == gapValue) continue;
        if (v >= kMsaBaseN) continue; // N / unknown
        t.count[v]++;
        if (strandIsRev[size_t(row - 1)]) t.rev[v]++; else t.fwd[v]++;
        t.covered++;
    }
    return t;
}

// --------------------------------------------------------------------------
// Pass 3: read base position at a set of columns of interest, per member row.
// --------------------------------------------------------------------------
// The member's MSA row is built by concatenating contiguous read segments, so
// its non-gap entries are exactly the read's bases in order with no holes. The
// read position at the row's first non-gap is rowStartReadPos; each subsequent
// non-gap advances by 1. We record only the requested columns.
//
// Returns a map keyed (uint64_t(row) << 32 | uint32_t(col)) -> read position.
inline void recoverReadPositions(
    const uint8_t* rowBases, int msaLen, int gapValue,
    uint32_t rowStartReadPos,
    int row,
    const std::vector<int>& columnsOfInterestSorted,
    std::unordered_map<uint64_t, uint32_t>& out)
{
    if (columnsOfInterestSorted.empty()) return;
    size_t want = 0;
    uint32_t readPos = rowStartReadPos;
    for (int c = 0; c < msaLen; ++c) {
        const bool nonGap = int(rowBases[c]) != gapValue;
        if (want < columnsOfInterestSorted.size() && c == columnsOfInterestSorted[want]) {
            if (nonGap) {
                out[(uint64_t(uint32_t(row)) << 32) | uint32_t(c)] = readPos;
            }
            ++want;
        }
        if (nonGap) ++readPos;
    }
}

// --------------------------------------------------------------------------
// Agreed-neighbor selection for a 2-base allele anchor.
// --------------------------------------------------------------------------
// Given a het site at backbone base bb (column C0) and an allele group (the
// member rows carrying that allele, plus optionally the backbone), pick the
// better of the two CONSECUTIVE-backbone neighbor pairs:
//     right (bb, bb+1)  or  left (bb-1, bb).
// For each candidate the neighbor's consensus base is the group majority; a read
// matches if it is non-gap at both columns and equals (alleleBase, neighborCons)
// in backbone order. The candidate with more matching reads wins (tie -> right).
struct NeighborChoice {
    bool   valid = false;
    int    leftBb = -1;      // min(bb, neighborBb): anchor's first backbone base.
    int    leftCol = -1;     // MSA column of leftBb.
    uint8_t base0 = 0xFF;    // base at leftBb (backbone order).
    uint8_t base1 = 0xFF;    // base at leftBb+1 (backbone order).
    uint32_t matchCount = 0; // reads agreeing on the 2-mer.
};

// Majority base among group rows at column c (non-gap, 0..3 only). 0xFF if none.
inline uint8_t groupMajorityBase(
    const uint8_t* const* msa, int c, int gapValue,
    const std::vector<int>& groupRows)
{
    uint32_t cnt[4] = {0, 0, 0, 0};
    for (int row : groupRows) {
        const uint8_t v = msa[row][c];
        if (int(v) == gapValue || v >= kMsaBaseN) continue;
        cnt[v]++;
    }
    int best = -1; uint32_t bestN = 0;
    for (int b = 0; b < 4; ++b) if (cnt[b] > bestN) { bestN = cnt[b]; best = b; }
    return best < 0 ? uint8_t(0xFF) : uint8_t(best);
}

// Count group rows matching the 2-mer (baseLeft at colLeft, baseRight at colRight).
inline uint32_t countTwoMerMatches(
    const uint8_t* const* msa, int colLeft, int colRight, int gapValue,
    uint8_t baseLeft, uint8_t baseRight, const std::vector<int>& groupRows)
{
    uint32_t n = 0;
    for (int row : groupRows) {
        const uint8_t l = msa[row][colLeft];
        const uint8_t r = msa[row][colRight];
        if (int(l) == gapValue || int(r) == gapValue) continue;
        if (l == baseLeft && r == baseRight) ++n;
    }
    return n;
}

// alleleBase is this group's base at the het column C0 (bb).
inline NeighborChoice chooseNeighbor(
    const uint8_t* const* msa, int gapValue,
    const ColumnMaps& maps, int bb, uint8_t alleleBase,
    const std::vector<int>& groupRows)
{
    auto evalPair = [&](int leftBb) -> NeighborChoice {
        NeighborChoice nc;
        const int rightBb = leftBb + 1;
        if (leftBb < 0 || rightBb >= maps.backboneLen) return nc;
        const int colLeft  = maps.bbToCol[size_t(leftBb)];
        const int colRight = maps.bbToCol[size_t(rightBb)];

        // The neighbor column is the one that is NOT the het column.
        const bool hetIsLeft = (leftBb == bb);
        const int neighborCol = hetIsLeft ? colRight : colLeft;
        const uint8_t neighborCons = groupMajorityBase(msa, neighborCol, gapValue, groupRows);
        if (neighborCons == 0xFF) return nc;

        const uint8_t baseL = hetIsLeft ? alleleBase : neighborCons;
        const uint8_t baseR = hetIsLeft ? neighborCons : alleleBase;
        const uint32_t m = countTwoMerMatches(msa, colLeft, colRight, gapValue,
                                              baseL, baseR, groupRows);
        nc.valid = true;
        nc.leftBb = leftBb;
        nc.leftCol = colLeft;
        nc.base0 = baseL;
        nc.base1 = baseR;
        nc.matchCount = m;
        return nc;
    };

    const NeighborChoice right = evalPair(bb);     // (bb, bb+1)
    const NeighborChoice left  = evalPair(bb - 1); // (bb-1, bb)

    // Tie -> right.
    if (right.valid && (!left.valid || right.matchCount >= left.matchCount)) return right;
    if (left.valid) return left;
    return NeighborChoice{};
}

} // namespace msahet
} // namespace dinara

#endif // DINARA_MSA_HET_KERNEL_HPP
