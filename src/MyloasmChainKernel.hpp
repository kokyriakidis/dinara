// Port of myloasm's chaining DP (myloasm src/mapping.rs: dp_anchors_v2/dp_inner,
// https://github.com/bluenote-1577/myloasm, commit 9531fbb), wired in as
// chainingMode=2 so it can be run head-to-head against dinara's existing
// hifiasm/minimap2-style chaining (run_main_dp_loop, AssemblerInvertedIndex.cpp)
// on real data.
//
// This is a literal port, including real myloasm quirks a from-scratch
// reimplementation would not introduce (see myloasmChainAnchors below) -- the
// goal is to test myloasm's actual behavior, not an idealized version of it.
//
// Pure, Assembler-independent functions (like MsaHetKernel.hpp) so they can be
// unit-tested directly (tests/test_myloasm_chaining.cpp) without pulling in the
// rest of the assembler; production code in AssemblerInvertedIndex.cpp includes
// this header too.
//
// --- Anchor encoding note ---
// myloasm packs strand into pos1's sign bit (pos1 = (rel_strand << 31) |
// selfOffset, mapping.rs:291) and bitwise-NOTs pos2 for reverse-strand anchors
// (mapping.rs:576-577) so both strands' distance arithmetic (dist2 = s2 - e2)
// comes out uniformly positive-when-collinear. Dinara's existing
// HifiasmKmerHit::offset field already has the same monotonicity property via a
// different transform: for opposite-strand hits it stores
// readLenB - 1 - posB (AssemblerInvertedIndex.cpp, "offDiff") -- an
// order-reversing affine function of posB, exactly like myloasm's bitwise-NOT
// (M - posB). dp_inner only ever uses DIFFERENCES between two anchors' pos2
// values (dist2 = s2 - e2), never an absolute pos2 value, and the constant
// additive offset between "readLenB - 1 - posB" and "M - posB" cancels out of
// every such difference (and out of relative sort order, since both are
// order-reversing with the same sign). So dinara's `offset` field can be used
// directly as myloasm's pos2, unmodified, for both strands -- confirmed here by
// derivation, and empirically by a dedicated reverse-strand unit test.

#pragma once

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

namespace dinara {

struct MyloasmChainParams {
    int32_t matchScore = 11;   // myloasm's own minimizer-chaining default (mapping.rs:939-945, match_score = k, default k=11).
    int32_t gapCost = 1;
    int32_t maxGap = 200;      // myloasm MAX_GAP_CHAINING (constants.rs:6).
    int32_t doubleGap = 10000; // myloasm CompareTwinReadOptions::default() (types.rs:1703-1723).
    int64_t maxSkip = 10;
    int64_t maxIter = 50;      // myloasm's chaining "band" default (mapping.rs:823-827).
    int64_t minChainLength = 3;
};

// Struct-of-arrays input, mirroring the subset of HifiasmKmerHit fields this
// kernel needs (self_offset, offset, strand) so it stays independent of that
// (TU-local) type and is trivial to build by hand in unit tests.
struct MyloasmAnchorInput {
    std::vector<uint32_t> selfOffset; // Query (read A) position -- HifiasmKmerHit::self_offset.
    std::vector<uint32_t> offset;     // Target (read B) position -- HifiasmKmerHit::offset.
    std::vector<uint8_t>  strand;     // 0 = same strand, 1 = opposite -- HifiasmKmerHit::strand.

    size_t size() const noexcept { return selfOffset.size(); }
};

struct MyloasmChain {
    int32_t score = 0;
    std::vector<int64_t> originalIndices; // Into the MyloasmAnchorInput arrays, in forward (genomic) order.
};

namespace detail {

inline int32_t myloasmEncodePos1(uint32_t selfOffset, uint8_t strand) noexcept
{
    const uint32_t raw = (uint32_t(strand) << 31) | selfOffset;
    return static_cast<int32_t>(raw);
}

inline int32_t myloasmAbs(int32_t x) noexcept
{
    return x < 0 ? -x : x;
}

// Faithful port of myloasm's dp_inner::<false> (mapping.rs:687-810). myloasm only
// ever instantiates the REV=false version in production (its own comment at
// mapping.rs:595 says the REV const generic is no longer needed since the pos1/
// pos2 encoding handles both strands), so that is the only variant ported here.
inline void myloasmDpInner(
    const std::vector<int32_t>& pos1s,
    const std::vector<int32_t>& pos2s,
    std::vector<int32_t>& f,
    std::vector<int32_t>& p,
    std::vector<int32_t>& t,
    const MyloasmChainParams& params) noexcept
{
    const int64_t n = int64_t(pos1s.size());
    int64_t st = 0;
    int32_t maxIi = -1;

    for(int64_t i = 0; i < n; ++i) {
        const int32_t s1 = pos1s[size_t(i)];
        const int32_t s2 = pos2s[size_t(i)];

        // Advance the pos1 window left boundary (O(n) total across all i).
        while(st < i && s1 > params.doubleGap + pos1s[size_t(st)]) {
            ++st;
        }
        const int64_t lo = (i >= params.maxIter) ? std::max(i - params.maxIter, st) : st;

        int32_t maxF = params.matchScore;
        int32_t maxJ = -1;
        int64_t nSkip = 0;
        int64_t endJ = lo;
        const int32_t strandI = s1 >> 31;

        for(int64_t j = i - 1; j >= lo; --j) {
            endJ = j;
            const int32_t e1 = pos1s[size_t(j)];
            const int32_t e2 = pos2s[size_t(j)];

            const int32_t dist1 = s1 - e1;
            const int32_t dist2 = s2 - e2;
            const bool sameStrand = (strandI == (e1 >> 31));
            if(dist1 <= 0 || dist2 <= 0 || dist2 > params.doubleGap || !sameStrand) {
                continue;
            }

            const int32_t gapPenalty = myloasmAbs(dist1 - dist2);
            if(gapPenalty > params.maxGap) {
                continue;
            }

            const int32_t kmerOverlapScore = std::min(std::min(dist1, dist2), params.matchScore);
            const int32_t sc = f[size_t(j)] + kmerOverlapScore - params.gapCost * gapPenalty;

            if(sc > maxF) {
                maxF = sc;
                maxJ = int32_t(j);
                if(nSkip > 0) {
                    --nSkip;
                }
            } else if(t[size_t(j)] == int32_t(i)) {
                ++nSkip;
                if(nSkip > params.maxSkip) {
                    break;
                }
            }

            const int32_t pj = p[size_t(j)];
            if(pj >= 0) {
                t[size_t(pj)] = int32_t(i);
            }
        }

        // max_ii fallback: rescan if the cached best predecessor fell out of the pos1 window.
        const bool rescan = (maxIi < 0) || (s1 > params.doubleGap + pos1s[size_t(maxIi)]);
        if(rescan) {
            int32_t best = INT32_MIN;
            maxIi = -1;
            for(int64_t j = i - 1; j >= lo; --j) {
                if(f[size_t(j)] > best) {
                    best = f[size_t(j)];
                    maxIi = int32_t(j);
                }
            }
        }
        if(maxIi >= 0 && int64_t(maxIi) < endJ) {
            const int32_t e1 = pos1s[size_t(maxIi)];
            const int32_t e2 = pos2s[size_t(maxIi)];
            const int32_t dist1 = s1 - e1;
            const int32_t dist2 = s2 - e2;
            if(dist1 > 0 && dist2 > 0 && dist2 <= params.doubleGap) {
                const int32_t gapPenalty = myloasmAbs(dist1 - dist2);
                if(gapPenalty <= params.maxGap) {
                    const int32_t kmerOverlapScore = std::min(std::min(dist1, dist2), params.matchScore);
                    const int32_t sc = f[size_t(maxIi)] + kmerOverlapScore - params.gapCost * gapPenalty;
                    if(sc > maxF) {
                        maxF = sc;
                        maxJ = maxIi;
                    }
                }
            }
        }

        f[size_t(i)] = maxF;
        p[size_t(i)] = maxJ;

        // Note: myloasm's update here is a plain score comparison, with no
        // distance/strand gate -- unlike dinara's own run_main_dp_loop, which
        // gates this update the way minimap2 does. Ported as-is.
        if(maxIi < 0 || f[size_t(maxIi)] < f[size_t(i)]) {
            maxIi = int32_t(i);
        }
    }
}

} // namespace detail

// Faithful port of myloasm's dp_anchors_v2 (mapping.rs:527-682): sorts anchors by
// myloasm's own (pos1, pos2) key, runs the DP, then greedily extracts chains from
// every candidate endpoint above a score heuristic, in descending score order,
// claiming anchors as it goes.
//
// Ported INCLUDING a real myloasm quirk: when a candidate's traceback collides
// with an already-claimed anchor, the anchors visited before the collision are
// still marked claimed (burned) even though the resulting chain is discarded
// (good_chain=false) -- and if that discarded prefix's length happens to be
// below minChainLength, the ENTIRE remaining (lower-scoring) candidate list is
// abandoned too (myloasm breaks the outer loop there, not just skips this one
// candidate). This is myloasm's actual production behavior, not a bug
// introduced by this port -- see mapping.rs:639-675.
inline std::vector<MyloasmChain> myloasmChainAnchors(
    const MyloasmAnchorInput& anchors,
    const MyloasmChainParams& params)
{
    const int64_t n = int64_t(anchors.size());
    if(n == 0) {
        return {};
    }

    // Sort by (pos1, pos2) exactly as myloasm's dp_anchors_v2 does (mapping.rs:548-558).
    // See the file header comment: dinara's `offset` field is used as pos2
    // unmodified, for both strands.
    std::vector<int64_t> sortedIndices(static_cast<size_t>(n));
    for(int64_t i = 0; i < n; ++i) {
        sortedIndices[size_t(i)] = i;
    }
    std::sort(sortedIndices.begin(), sortedIndices.end(),
        [&](int64_t a, int64_t b) noexcept {
            const int32_t pa = detail::myloasmEncodePos1(anchors.selfOffset[size_t(a)], anchors.strand[size_t(a)]);
            const int32_t pb = detail::myloasmEncodePos1(anchors.selfOffset[size_t(b)], anchors.strand[size_t(b)]);
            if(pa != pb) {
                return pa < pb;
            }
            return int32_t(anchors.offset[size_t(a)]) < int32_t(anchors.offset[size_t(b)]);
        });

    std::vector<int32_t> pos1s(static_cast<size_t>(n));
    std::vector<int32_t> pos2s(static_cast<size_t>(n));
    for(int64_t i = 0; i < n; ++i) {
        const int64_t orig = sortedIndices[size_t(i)];
        pos1s[size_t(i)] = detail::myloasmEncodePos1(anchors.selfOffset[size_t(orig)], anchors.strand[size_t(orig)]);
        pos2s[size_t(i)] = int32_t(anchors.offset[size_t(orig)]);
    }

    std::vector<int32_t> f(size_t(n), params.matchScore);
    std::vector<int32_t> p(size_t(n), -1);
    std::vector<int32_t> t(size_t(n), -1);
    detail::myloasmDpInner(pos1s, pos2s, f, p, t, params);

    // Chain reconstruction (mapping.rs:621-675). Repurpose `t` as the claimed marker.
    std::fill(t.begin(), t.end(), 0);

    std::vector<std::pair<int32_t, int64_t>> bestIndicesOrdered; // (score, index).
    for(int64_t i = 0; i < n; ++i) {
        if(f[size_t(i)] > int32_t(params.minChainLength) * params.matchScore / 2) {
            bestIndicesOrdered.emplace_back(f[size_t(i)], i);
        }
    }
    std::stable_sort(bestIndicesOrdered.begin(), bestIndicesOrdered.end(),
        [](const std::pair<int32_t, int64_t>& a, const std::pair<int32_t, int64_t>& b) noexcept {
            return a.first > b.first;
        });

    std::vector<MyloasmChain> chains;
    for(const auto& scoreAndIndex : bestIndicesOrdered) {
        const int32_t score = scoreAndIndex.first;
        const int64_t bestIndex = scoreAndIndex.second;
        if(t[size_t(bestIndex)] != 0) {
            continue;
        }

        std::vector<int64_t> chainSorted; // Sorted-order indices, endpoint first.
        int64_t idx = bestIndex;
        bool goodChain = true;
        while(idx >= 0) {
            if(t[size_t(idx)] != 0) {
                goodChain = false;
                break;
            }
            t[size_t(idx)] = 1;
            chainSorted.push_back(idx);
            idx = p[size_t(idx)];
        }

        if(int64_t(chainSorted.size()) < params.minChainLength) {
            break; // Matches myloasm: abandon all remaining candidates too, not just this one.
        }

        if(goodChain) {
            std::reverse(chainSorted.begin(), chainSorted.end());
            MyloasmChain chain;
            chain.score = score;
            chain.originalIndices.reserve(chainSorted.size());
            for(const int64_t sortedIdx : chainSorted) {
                chain.originalIndices.push_back(sortedIndices[size_t(sortedIdx)]);
            }
            chains.push_back(std::move(chain));
        }
    }

    return chains;
}

} // namespace dinara
