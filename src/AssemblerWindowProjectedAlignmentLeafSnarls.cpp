/// @file AssemblerWindowProjectedAlignmentLeafSnarls.cpp
/// @brief Leaf-snarl (multi-column MNP/het-block) detection per anchor window,
///        from ProjectedAlignment (the SAME fast pairwise aligner
///        computeBaseAlignmentsAndStore uses for the whole-dataset overlap
///        graph) instead of ksw2's per-segment banded DP.
///
/// Motivation: the ksw2-pairwise leaf-snarl detector (AssemblerWindowKsw2LeafSnarls.cpp)
/// calls ksw_extd2_sse ONCE PER (member, inter-anchor segment) pair --
/// measured at 27015 calls for one test region, ~8x more than the
/// interval-abPOA engine's 3346 shared abpoa_msa calls for the SAME windows,
/// because ksw2 never shares work across members touching the same interval.
/// ProjectedAlignment::QuickRawSparse takes a member's ENTIRE pin list
/// (Alignment.ordinals) in ONE constructor call and loops over all
/// inter-anchor segments INTERNALLY -- so this detector pays ONE call per
/// member, not one per segment, AND skips the actual alignment DP entirely
/// for any segment whose sequences are byte-identical (the common case for
/// HiFi data), via constructQuickRawSparse's memcmp fast path
/// (ProjectedAlignment.cpp:362-366). This is the same mechanism that makes
/// computeBaseAlignmentsAndStore fast for the whole overlap graph (5.07s for
/// 104191 alignments on one test region) -- reused here for window-relative
/// per-member alignment instead of read-vs-read overlap alignment.
///
/// Produces SPARSE evidence (only mismatch/indel positions, not a dense
/// per-column walk), fed into findLeafSnarlsFromSparseEvidence
/// (WindowIntervalPoa.hpp), which defaults every covered position to
/// OnBackbone and applies the sparse exceptions -- the natural fit for a
/// diff-only aligner, as opposed to findLeafSnarlsFromProfiles's dense
/// alignedCols model (built for engines that walk every column explicitly).
///
/// computeWindowProjectedAlignmentLeafSnarls itself is verification only:
/// writes nothing to AnchorWindow::hetBubbles, just logs -- a comparison
/// point against the ksw2-pairwise, interval-abPOA, and abPOA-shared-graph
/// leaf-snarl detectors on the same windows. This file also defines
/// projAlnDetectHetBubblesAllWindows, the PRODUCTION counterpart that feeds
/// the same sparse evidence through buildHetBubbleFromLeafSnarl and writes
/// real AnchorWindow::hetBubbles, selectable as DINARA_HET_ENGINE=projaln.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "Alignment.hpp"
#include "AssemblerOptions.hpp"
#include "HetAnchorK.hpp"
#include "OverlapCigarStore.hpp"
#include "ProjectedAlignment.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "WindowIntervalPoa.hpp"
#include "invalid.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace dinara;

namespace {

// A shared anchor pin: backbone marker ordinal + position, member marker
// ordinal + position. Ordinals feed Alignment.ordinals (what ProjectedAlignment
// actually needs); positions are only used for the LIS colinearity filter and
// logging, mirroring the ksw2/abPOA detectors' own pin structs.
struct PaPin { uint32_t bbOrd; uint32_t bbPos; uint32_t cOrd; uint32_t cPos; };

// Keep the longest colinear (strictly diagonal-increasing) run of pins by
// backbone position: sort by backbone position, then LIS over read position.
// Identical in intent to lsnLisByCPos/lisByCPos in the ksw2 detectors.
void paLisByCPos(vector<PaPin>& pins) {
    const uint32_t n = uint32_t(pins.size());
    if (n < 2) return;
    vector<uint32_t> tails;
    vector<int32_t> pred(n, -1);
    for (uint32_t i = 0; i < n; i++) {
        uint32_t lo = 0, hi = uint32_t(tails.size());
        while (lo < hi) {
            const uint32_t mid = (lo + hi) / 2;
            if (pins[tails[mid]].cPos < pins[i].cPos) lo = mid + 1;
            else hi = mid;
        }
        if (lo > 0) pred[i] = int32_t(tails[lo - 1]);
        if (lo == uint32_t(tails.size())) tails.push_back(i);
        else tails[lo] = i;
    }
    vector<PaPin> kept;
    for (int32_t c = int32_t(tails.back()); c >= 0; c = pred[uint32_t(c)])
        kept.push_back(pins[uint32_t(c)]);
    reverse(kept.begin(), kept.end());
    pins.swap(kept);
}

// Builds a SparseMemberEvidence directly from a STORED pairwise CIGAR
// (OverlapCigarStore, from computeBaseAlignmentsAndStore) between the
// window's backbone and one member, instead of computing a fresh
// ProjectedAlignment -- avoiding redundant alignment work for the ~92% of
// window members whose pair already has a stored CIGAR. Falls back to fresh
// ProjectedAlignment (the caller's existing path) whenever this returns
// false: no stored CIGAR, or no coverage overlapping the window.
//
// Coordinate derivation. computeBaseAlignmentsAndStoreThreadFunction always
// builds the ProjectedAlignment that produced this CIGAR from
// orientedReadIds[0] = OrientedReadId(ad.readIds[0], 0) and
// orientedReadIds[1] = OrientedReadId(ad.readIds[1], ad.isSameStrand?0:1)
// (AssemblerComputeAlignments.cpp:432-433) -- i.e. the CIGAR's own (xk,yk)
// walk is ALWAYS in that fixed native frame, regardless of which strand the
// CALLER's bbOid/cOid happen to be. So for each side independently:
//   nativeIsStrand1 = (this side is "read1" AND !ad.isSameStrand)
//   flip = nativeIsStrand1 != (requested OrientedReadId's strand == 1)
// A raw single-base position p becomes `len-1-p` under flip; a raw
// half-open range [s,e) becomes [len-e,len-s). A base value gets
// complemented under flip too (a flip is a real reverse-complement, not
// just renumbering) -- needed for the member's altBase, not for the
// backbone (bbSeqVec already supplies the backbone's own base in bbOid's
// frame; the CIGAR only ever supplies bbPos as an index into it).
//
// Delta (readPos - bbPos) breakpoints: ONE per match/mismatch run, not
// accumulated by propagating +len/-len through insertions/deletions the
// way the ProjectedAlignment path does. Reason: when bbFlip is active, the
// CIGAR's native (xk,yk) walk order is the REVERSE of bbPos-ascending
// order, and an insertion's effect on delta going forward flips sign under
// that reversal (verified by hand with a worked numeric example while
// developing this) -- propagating a running "+len for insertion, -len for
// deletion" accumulator would silently pick up the wrong sign whenever
// bbFlip is set. Instead, each match/mismatch run's OWN delta is computed
// directly from ONE of its own base pairs (bbRaw offset 0 vs its
// corresponding tRaw), independently transformed to final (bbPos,readPos)
// via bbTransformPoint/targetTransformPoint -- never by propagating a
// signed step from a neighboring run. Since a run's bbRaw range and tRaw
// range have the SAME length and advance in lockstep, delta is constant
// across the whole run regardless of flip direction, so one sample per run
// suffices, and consecutive runs' independently-computed deltas already
// correctly reflect whatever indel separates them.
bool buildSparseEvidenceFromCigar(
    const Assembler& assembler,
    ReadId backboneReadId,
    uint32_t backboneLen,
    OrientedReadId bbOid,
    OrientedReadId cOid,
    uint32_t alignmentId,
    uint32_t windowBbBegin,
    uint32_t windowBbEnd,
    SparseMemberEvidence& outEv)
{
    const auto& cigarStore = assembler.getOverlapCigarStore();
    const auto& ad = assembler.alignmentData[alignmentId];
    const auto& info = ad.info;
    if (info.cigarOffset == uint32_t(-1) || info.cigarTokenCount == 0) return false;

    const bool bbIsRead0 = (ad.readIds[0] == backboneReadId);
    const ReadId targetReadId = bbIsRead0 ? ad.readIds[1] : ad.readIds[0];
    const uint32_t targetLen = uint32_t(assembler.getReads().getRead(targetReadId).baseCount);
    const bool isRev = !ad.isSameStrand; // true iff read1's native frame is strand1

    // needsRcBb/needsRcTarget: true iff THIS side is the one that's "read1"
    // (native frame strand1) when isRev. Exactly one can be true (mutually
    // exclusive on bbIsRead0), matching that exactly one side is "read1".
    const bool needsRcBb = !bbIsRead0 && isRev;
    const bool needsRcTarget = bbIsRead0 && isRev;
    const bool bbFlip = needsRcBb != (bbOid.getStrand() == 1);
    const bool targetFlip = needsRcTarget != (cOid.getStrand() == 1);

    // Reverse-complementing BOTH sides of a real overlap together preserves
    // the SAME alignment (RC(A) vs RC(B) is the mirror image of A vs B);
    // flipping only ONE side describes a DIFFERENT, unrelated relationship.
    // So a valid (bbOid, cOid) pair consistent with this specific `ad` must
    // have bbFlip == targetFlip -- if they differ, the requested pairing
    // does NOT correspond to the same overlap this CIGAR describes (e.g. a
    // read with an internal near-palindrome/duplication can have the anchor
    // graph settle on the opposite relative orientation from whatever
    // single candidate overlap survived into the CIGAR store for these two
    // raw reads). Verified on real data: a mismatch here produced a
    // monotonically-decreasing readPos as bbPos increased -- a genuinely
    // wrong result, not just a differently-placed-but-valid one. Bail out
    // to the fresh-ProjectedAlignment fallback instead.
    if (bbFlip != targetFlip) return false;

    uint64_t read1Start;
    if (ad.isSameStrand) {
        read1Start = ad.ts;
    } else {
        const uint32_t rlen = uint32_t(assembler.getReads().getRead(ad.readIds[1]).baseCount);
        read1Start = rlen - ad.te;
    }

    struct DeltaSample { uint32_t bbPos; int64_t delta; };
    vector<DeltaSample> deltaSamples;
    vector<SparseSnp> snps;
    vector<pair<uint32_t, uint32_t>> deletionRanges;
    uint32_t covMin = UINT32_MAX, covMax = 0;

    auto bbTransformPoint = [&](uint32_t rawP) {
        return bbFlip ? (backboneLen - 1 - rawP) : rawP;
    };
    auto bbTransformRange = [&](uint32_t rawLo, uint32_t rawHi) {
        return bbFlip ? std::make_pair(backboneLen - rawHi, backboneLen - rawLo)
                      : std::make_pair(rawLo, rawHi);
    };
    auto targetTransformPoint = [&](uint32_t rawP) {
        return targetFlip ? (targetLen - 1 - rawP) : rawP;
    };

    cigarStore.forEachOpWithPositions(
        info.cigarOffset, info.cigarTokenCount,
        ad.qs, read1Start,
        [&](uint8_t op, uint32_t len, uint64_t xk64, uint64_t yk64) {
            const uint32_t xk = uint32_t(xk64), yk = uint32_t(yk64);
            if (op == 0 || op == 1) {
                const uint32_t bbRawStart = bbIsRead0 ? xk : yk;
                const uint32_t tRawStart = bbIsRead0 ? yk : xk;
                const uint32_t bbRawEnd = std::min(bbRawStart + len, backboneLen);
                const auto [bbLo, bbHi] = bbTransformRange(bbRawStart, bbRawEnd);
                if (bbLo < covMin) covMin = bbLo;
                if (bbHi > covMax) covMax = bbHi;
                // One delta sample for this whole run, taken at whichever
                // raw index i0 maps to the run's SMALLEST final bbPos
                // (bbLo) -- NOT always i0=0. When bbFlip is set,
                // bbTransformPoint reverses order within the run, so
                // bbRawStart (i0=0) maps to the run's LARGEST final bbPos
                // instead; sampling there would attribute this run's delta
                // to the wrong end of its range once sorted with other
                // runs' breakpoints (caught via a downstream monotonicity
                // assertion on a real dataset -- verified by hand that
                // bbTransformPoint(bbRawStart) landed at bbHi-1, not bbLo,
                // whenever bbFlip was set).
                {
                    const uint32_t i0 = bbFlip ? (bbRawEnd - bbRawStart - 1) : 0;
                    const uint32_t tRawAtLo = tRawStart + i0;
                    deltaSamples.push_back({
                        bbLo,
                        int64_t(targetTransformPoint(tRawAtLo)) - int64_t(bbLo)});
                }
                if (op == 1) { // mismatch run: every base in it differs
                    for (uint32_t b = 0; b < len; b++) {
                        const uint32_t bbRaw = bbRawStart + b;
                        const uint32_t tRaw = tRawStart + b;
                        const uint32_t bbPos = bbTransformPoint(bbRaw);
                        if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;
                        // tRaw is a position in the CIGAR's own native walk
                        // frame, which is target-forward UNLESS
                        // needsRcTarget (then it's target-reverse, by
                        // construction of read1Start above) -- normalize to
                        // a real forward-storage index FIRST. targetFlip
                        // (needsRcTarget XOR cOid-is-reverse) is valid for
                        // POSITION math (targetTransformPoint below: two
                        // reversals with the same length cancel correctly),
                        // but is NOT the right condition for the base VALUE:
                        // using it directly to both index and complement
                        // silently read the wrong raw base whenever
                        // needsRcTarget was true (caught by hand: readPos
                        // matched exactly between this path and the trusted
                        // fresh-ProjectedAlignment reference, but the
                        // reported base didn't -- a real base-identity bug,
                        // not alignment-placement ambiguity).
                        const uint32_t forwardIdx = needsRcTarget ? (targetLen - 1 - tRaw) : tRaw;
                        uint8_t base = assembler.getReads().getRead(targetReadId)[forwardIdx].value;
                        if (cOid.getStrand() == 1) base = uint8_t((~base) & 3);
                        const uint32_t readPos = targetTransformPoint(tRaw);
                        snps.push_back(SparseSnp{bbPos, base, readPos});
                    }
                }
            } else if (op == CigarOpIns || op == CigarOpDel) {
                // Backbone read (xk if bbIsRead0, else yk) consumed by this indel
                // iff the op consumes that read's sequence → deletion on the
                // backbone. Otherwise it is an insertion relative to the backbone.
                const bool bbConsumed = bbIsRead0 ? opConsumesQuery(op)
                                                  : opConsumesTarget(op);
                if (bbConsumed) { // deletion on backbone: backbone-consuming
                    const uint32_t bbRawStart = bbIsRead0 ? xk : yk;
                    const uint32_t bbRawEnd = std::min(bbRawStart + len, backboneLen);
                    const auto [bbLo, bbHi] = bbTransformRange(bbRawStart, bbRawEnd);
                    if (bbLo < covMin) covMin = bbLo;
                    if (bbHi > covMax) covMax = bbHi;
                    deletionRanges.push_back({bbLo, bbHi});
                }
                // Insertions need no explicit delta sample: the match/
                // mismatch runs immediately before and after it already
                // each carry their own independently-correct delta.
            }
        });

    if (covMin >= covMax) return false;

    outEv.oid = cOid;
    outEv.bbCovBegin = std::max(covMin, windowBbBegin);
    outEv.bbCovEnd = std::min(covMax, windowBbEnd);
    if (outEv.bbCovBegin >= outEv.bbCovEnd) return false;

    outEv.snps.clear();
    for (const SparseSnp& s : snps) outEv.snps.push_back(s);
    sort(outEv.snps.begin(), outEv.snps.end(),
        [](const SparseSnp& a, const SparseSnp& b) { return a.bbPos < b.bbPos; });

    outEv.deletionRanges.clear();
    for (const auto& d : deletionRanges) {
        const uint32_t a = std::max(d.first, windowBbBegin);
        const uint32_t b = std::min(d.second, windowBbEnd);
        if (a < b) outEv.deletionRanges.push_back({a, b});
    }
    sort(outEv.deletionRanges.begin(), outEv.deletionRanges.end());

    // deltaSamples are keyed by their already-transformed bbPos (one per
    // match/mismatch run), so sorting is the only reordering needed --
    // regardless of whether bbFlip reversed the CIGAR's native walk order.
    // deltaAtStart = the smallest-bbPos sample's own delta (covers queries
    // at or before it; harmless overlap with deltaBreakpoints[0] itself).
    sort(deltaSamples.begin(), deltaSamples.end(),
        [](const DeltaSample& a, const DeltaSample& b) { return a.bbPos < b.bbPos; });
    outEv.deltaBreakpoints.clear();
    if (!deltaSamples.empty()) {
        outEv.deltaAtStart = deltaSamples.front().delta;
        for (const DeltaSample& s : deltaSamples) outEv.deltaBreakpoints.push_back({s.bbPos, s.delta});
    } else {
        outEv.deltaAtStart = 0;
    }

    if (getenv("DINARA_CIGAR_REUSE_SELFCHECK") != nullptr) {
        for (size_t i = 1; i < deltaSamples.size(); i++) {
            const int64_t rp0 = int64_t(deltaSamples[i - 1].bbPos) + deltaSamples[i - 1].delta;
            const int64_t rp1 = int64_t(deltaSamples[i].bbPos) + deltaSamples[i].delta;
            if (rp1 < rp0) {
                cout << "      selfCheckFAIL member=" << cOid
                     << " sample[" << (i - 1) << "]=(bb" << deltaSamples[i - 1].bbPos
                     << ",delta" << deltaSamples[i - 1].delta << ",rp" << rp0 << ")"
                     << " sample[" << i << "]=(bb" << deltaSamples[i].bbPos
                     << ",delta" << deltaSamples[i].delta << ",rp" << rp1 << ")"
                     << " bbFlip=" << bbFlip << " targetFlip=" << targetFlip
                     << " bbIsRead0=" << bbIsRead0 << " isSameStrand=" << ad.isSameStrand << endl;
            }
        }
    }

    outEv.noisyRanges.clear();
    return true;
}

// Builds a SparseMemberEvidence via a fresh ProjectedAlignment over this
// member's own shared-anchor pins -- the ORIGINAL (pre-CIGAR-reuse) path,
// kept as: (a) the fallback for any pair with no usable stored CIGAR (the
// anchor-graph-only gap), and (b) the trusted reference for cross-checking
// buildSparseEvidenceFromCigar's output (DINARA_CIGAR_REUSE_VALIDATE=1).
// pinsRaw is sorted and LIS-filtered in place. Returns false if fewer than
// 2 colinear pins survive, or the resulting alignment has no coverage
// inside the window.
bool buildSparseEvidenceFromProjectedAlignment(
    const Assembler& assembler,
    OrientedReadId bbOid,
    OrientedReadId cOid,
    vector<PaPin>& pinsRaw,
    uint32_t windowBbBegin,
    uint32_t windowBbEnd,
    int64_t dpMatchScore,
    int64_t dpMismatchScore,
    int64_t dpGapOpen1,
    int64_t dpGapExtend1,
    SparseMemberEvidence& outEv)
{
    if (pinsRaw.size() < 2) return false;

    sort(pinsRaw.begin(), pinsRaw.end(),
        [](const PaPin& a, const PaPin& b) {
            return a.bbPos < b.bbPos || (a.bbPos == b.bbPos && a.cPos < b.cPos);
        });
    paLisByCPos(pinsRaw);
    if (pinsRaw.size() < 2) return false;

    Alignment sparseAlign;
    sparseAlign.ordinals.reserve(pinsRaw.size());
    for (const PaPin& p : pinsRaw)
        sparseAlign.ordinals.push_back({p.bbOrd, p.cOrd});

    const ProjectedAlignment pa(
        assembler, {bbOid, cOid}, sparseAlign,
        ProjectedAlignment::Method::QuickRawSparse,
        dpMatchScore, dpMismatchScore, dpGapOpen1, dpGapExtend1);

    SparseMemberEvidence ev;
    ev.oid = cOid;
    ev.bbCovBegin = max(pa.cigarRead0Start, windowBbBegin);
    ev.bbCovEnd = min(pa.cigarRead0End, windowBbEnd);
    if (ev.bbCovBegin >= ev.bbCovEnd) return false;

    for (const auto& mm : pa.sparseMismatches) {
        if (mm.position0 < windowBbBegin || mm.position0 >= windowBbEnd) continue;
        ev.snps.push_back(SparseSnp{mm.position0, mm.base1, mm.position1});
    }
    ev.deltaAtStart = int64_t(pa.cigarRead1Start) - int64_t(pa.cigarRead0Start);
    ev.deltaBreakpoints.reserve(pa.sparseIndels.size());
    {
        int64_t delta = ev.deltaAtStart;
        for (const auto& ind : pa.sparseIndels) {
            if (ind.op == 'D') {
                delta -= int64_t(ind.length);
                ev.deltaBreakpoints.push_back({ind.position0 + ind.length, delta});
                const uint32_t a = max(ind.position0, windowBbBegin);
                const uint32_t b = min(ind.position0 + ind.length, windowBbEnd);
                if (a < b) ev.deletionRanges.push_back({a, b});
            } else { // 'I'
                delta += int64_t(ind.length);
                ev.deltaBreakpoints.push_back({ind.position0, delta});
            }
        }
    }
    sort(ev.snps.begin(), ev.snps.end(),
        [](const SparseSnp& a, const SparseSnp& b) { return a.bbPos < b.bbPos; });
    sort(ev.deletionRanges.begin(), ev.deletionRanges.end());

    outEv = std::move(ev);
    return true;
}

// Cross-check helper for DINARA_CIGAR_REUSE_VALIDATE=1: true if two
// SparseMemberEvidence instances for the SAME member agree on everything
// that feeds downstream classification/pinning -- coverage bounds, the
// exact snp set (bbPos+altBase+readPos), and the exact deletion-range set.
// deltaBreakpoints/deltaAtStart are NOT compared directly (the two builders
// are free to place breakpoints differently, e.g. one-per-CIGAR-run vs
// one-per-indel-event, while still agreeing on every actual query result);
// instead the caller probes sparseColAt at representative positions.
// Compares only within [lo,hi) = the two evidences' OVERLAPPING coverage --
// CIGAR reuse legitimately covers more of a member's true overlap than the
// window-local pin set (the pin-based path is bounded by whatever anchors
// happen to fall in THIS window), so requiring identical bbCovBegin/End
// would flag a real improvement as a mismatch.
bool sparseEvidenceCoreEqual(const SparseMemberEvidence& a, const SparseMemberEvidence& b) {
    const uint32_t lo = max(a.bbCovBegin, b.bbCovBegin);
    const uint32_t hi = min(a.bbCovEnd, b.bbCovEnd);
    if (lo >= hi) return false;
    auto within = [lo, hi](uint32_t p) { return p >= lo && p < hi; };

    vector<SparseSnp> snpsA, snpsB;
    for (const auto& s : a.snps) if (within(s.bbPos)) snpsA.push_back(s);
    for (const auto& s : b.snps) if (within(s.bbPos)) snpsB.push_back(s);
    sort(snpsA.begin(), snpsA.end(), [](const SparseSnp& x, const SparseSnp& y) { return x.bbPos < y.bbPos; });
    sort(snpsB.begin(), snpsB.end(), [](const SparseSnp& x, const SparseSnp& y) { return x.bbPos < y.bbPos; });
    if (snpsA.size() != snpsB.size()) return false;
    for (size_t i = 0; i < snpsA.size(); i++) {
        if (snpsA[i].bbPos != snpsB[i].bbPos) return false;
        if (snpsA[i].altBase != snpsB[i].altBase) return false;
        if (snpsA[i].readPos != snpsB[i].readPos) return false;
    }

    auto clipRanges = [lo, hi](const vector<pair<uint32_t, uint32_t>>& in) {
        vector<pair<uint32_t, uint32_t>> out;
        for (auto d : in) {
            d.first = max(d.first, lo);
            d.second = min(d.second, hi);
            if (d.first < d.second) out.push_back(d);
        }
        sort(out.begin(), out.end());
        return out;
    };
    return clipRanges(a.deletionRanges) == clipRanges(b.deletionRanges);
}

// Cross-check helper for DINARA_SPARSE_CLASSIFY_VALIDATE=1: true if two
// LeafSnarl lists are identical (same runs in the same order, same alleles
// with the same base keys and the same member SETS -- member order within an
// allele can legitimately differ only if the two callers iterated members in
// a different order, which neither does here, so this checks order too).
bool leafSnarlsEqual(const vector<LeafSnarl>& a, const vector<LeafSnarl>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i].start != b[i].start || a[i].end != b[i].end) return false;
        if (a[i].alleles.size() != b[i].alleles.size()) return false;
        for (size_t j = 0; j < a[i].alleles.size(); j++) {
            if (a[i].alleles[j].bases != b[i].alleles[j].bases) return false;
            vector<OrientedReadId> ma = a[i].alleles[j].members;
            vector<OrientedReadId> mb = b[i].alleles[j].members;
            sort(ma.begin(), ma.end());
            sort(mb.begin(), mb.end());
            if (ma != mb) return false;
        }
    }
    return true;
}

// Builds an AnchorWindow::HetBubble from a single-column leaf snarl (one
// divergent backbone column, bracketed by a leading/trailing hom), mirroring
// emitHetBubblesFromProfiles's per-site construction (WindowHetProfiles.hpp)
// exactly -- same ref/alt arm structure, same leadHom/hom bracketing, same
// post-pinning minSupport re-check, same "both homs must bracket every arm"
// requirement -- but pinning members via the SPARSE evidence model
// (sparsePinnedKmerCol/sparsePinnedPointCol, WindowIntervalPoa.hpp) instead
// of KwMemberProfile's dense colAt.
//
// Only single-column snarls apply: AnchorWindow::HetAnchor::alleleBase is a
// single uint8_t, so it cannot represent a multi-base MNP allele without
// extending the anchor format itself (out of scope here) -- the caller
// filters to snarl.end - snarl.start == 2 before calling this.
//
// Returns false (bubble left untouched) if the site doesn't survive pinning
// (either arm falls below minSupport after pinning, or the two homs don't
// bracket every surviving arm).
bool buildHetBubbleFromLeafSnarl(
    const LeafSnarl& snarl,
    const unordered_map<uint64_t, const SparseMemberEvidence*>& evidenceByOid,
    const vector<uint8_t>& bbSeqVec,
    OrientedReadId bbOid,
    uint32_t windowBbBegin,
    uint64_t minSupport,
    AnchorWindow::HetBubble& outBubble,
    string* failReason = nullptr)
{
    const uint32_t pos = snarl.start + 1; // the single divergent column
    const bool hetK0 = (hetAnchorK() == 0);
    const uint8_t predBase = bbSeqVec[pos - 1];

    AnchorWindow::HetBubble bubble;
    bubble.backboneOffset = pos - windowBbBegin;

    vector<const SparseMemberEvidence*> homMemberEv;

    auto lookup = [&](OrientedReadId oid) -> const SparseMemberEvidence* {
        auto it = evidenceByOid.find(oid.getValue());
        return it == evidenceByOid.end() ? nullptr : it->second;
    };

    // Reference arm (allele 0, per classifyLeafSnarls's seeding order --
    // WindowIntervalPoa.hpp always seeds the ref allele first).
    AnchorWindow::HetAnchor refArm;
    refArm.backboneOffset = bubble.backboneOffset;
    refArm.predBase = predBase;
    refArm.alleleBase = bbSeqVec[pos];
    refArm.isRef = true;
    refArm.members.push_back({bbOid, hetK0 ? pos : pos - 1});
    for (const OrientedReadId& oid : snarl.alleles[0].members) {
        if (oid == bbOid) continue;
        const SparseMemberEvidence* ev = lookup(oid);
        if (!ev) continue;
        uint32_t readPos = 0;
        const bool ok = hetK0
            ? sparsePinnedPointCol(*ev, pos, refArm.alleleBase, refArm.alleleBase, readPos)
            : sparsePinnedKmerCol(*ev, pos - 1, predBase, refArm.alleleBase,
                                  predBase, refArm.alleleBase, readPos);
        if (!ok) continue;
        refArm.members.push_back({oid, readPos});
        homMemberEv.push_back(ev);
    }
    // Post-pinning re-check (the bug we found and fixed in
    // AssemblerWindowKsw2HetSites.cpp: the site's raw support count is not
    // the same as the count that survives pinning).
    if (refArm.members.size() < size_t(minSupport)) {
        if (failReason) *failReason = "refArm too small after pinning: " +
            std::to_string(refArm.members.size()) + " (raw=" +
            std::to_string(snarl.alleles[0].members.size()) + ")";
        return false;
    }
    bubble.alleles.push_back(std::move(refArm));

    // Alternate arms.
    for (size_t ai = 1; ai < snarl.alleles.size(); ai++) {
        const LeafSnarlAllele& allele = snarl.alleles[ai];
        AnchorWindow::HetAnchor arm;
        arm.backboneOffset = bubble.backboneOffset;
        arm.predBase = predBase;
        arm.alleleBase = allele.bases.empty() ? 0 : allele.bases[0]; // single column
        arm.isRef = false;
        for (const OrientedReadId& oid : allele.members) {
            const SparseMemberEvidence* ev = lookup(oid);
            if (!ev) continue;
            uint32_t readPos = 0;
            const bool ok = hetK0
                ? sparsePinnedPointCol(*ev, pos, arm.alleleBase, refArm.alleleBase, readPos)
                : sparsePinnedKmerCol(*ev, pos - 1, predBase, arm.alleleBase,
                                      predBase, refArm.alleleBase, readPos);
            if (!ok) continue;
            arm.members.push_back({oid, readPos});
            homMemberEv.push_back(ev);
        }
        if (arm.members.size() < size_t(minSupport)) {
            if (failReason) *failReason = "altArm too small after pinning: " +
                std::to_string(arm.members.size()) + " (raw=" +
                std::to_string(allele.members.size()) + ")";
            continue;
        }
        bubble.alleles.push_back(std::move(arm));
    }
    if (bubble.alleles.size() < 2) {
        if (failReason && failReason->empty()) *failReason = "fewer than 2 arms survived pinning";
        return false;
    }

    // Leading hom [predPrevBase, predBase] at predPrev (column pos-2).
    {
        const uint32_t leadOff = (pos - 2) - windowBbBegin;
        AnchorWindow::HetAnchor leadHom;
        leadHom.backboneOffset = leadOff;
        leadHom.predBase = bbSeqVec[pos - 2];
        leadHom.alleleBase = predBase;
        leadHom.isRef = true;
        leadHom.members.push_back({bbOid, pos - 2});
        for (const SparseMemberEvidence* ev : homMemberEv) {
            uint32_t readPos = 0;
            const bool ok = hetK0
                ? sparsePinnedPointCol(*ev, pos - 2, leadHom.predBase, leadHom.predBase, readPos)
                : sparsePinnedKmerCol(*ev, pos - 2, leadHom.predBase, leadHom.alleleBase,
                                      leadHom.predBase, predBase, readPos);
            if (!ok) continue;
            leadHom.members.push_back({ev->oid, readPos});
        }
        if (leadHom.members.size() > 1) {
            bubble.predBackboneOffset = leadOff;
            bubble.leadHom = std::move(leadHom);
        }
    }

    // Trailing hom [succBase, nextBase] at commonSucc (column pos+1).
    {
        const uint32_t succOff = (pos + 1) - windowBbBegin;
        AnchorWindow::HetAnchor homAnchor;
        homAnchor.backboneOffset = succOff;
        homAnchor.predBase = bbSeqVec[pos + 1];
        homAnchor.alleleBase = bbSeqVec[pos + 2];
        homAnchor.isRef = true;
        homAnchor.members.push_back({bbOid, pos + 1});
        for (const SparseMemberEvidence* ev : homMemberEv) {
            uint32_t readPos = 0;
            const bool ok = hetK0
                ? sparsePinnedPointCol(*ev, pos + 1, homAnchor.predBase, homAnchor.predBase, readPos)
                : sparsePinnedKmerCol(*ev, pos + 1, homAnchor.predBase, homAnchor.alleleBase,
                                      homAnchor.predBase, homAnchor.alleleBase, readPos);
            if (!ok) continue;
            homAnchor.members.push_back({ev->oid, readPos});
        }
        if (homAnchor.members.size() > 1) {
            bubble.succBackboneOffset = succOff;
            bubble.hom = std::move(homAnchor);
        }
    }

    // Both homs must bracket EVERY allele arm (see the design note on
    // AnchorWindow::HetBubble, AnchorWindows.hpp, for why: a direct
    // backbone->arm edge has an empty read intersection since a minority
    // allele's reads diverge at the flanking k=50 anchor too).
    {
        auto oidSet = [](const AnchorWindow::HetAnchor& h) {
            vector<OrientedReadId> s;
            s.reserve(h.members.size());
            for (const auto& m : h.members) s.push_back(m.orientedReadId);
            sort(s.begin(), s.end());
            s.erase(unique(s.begin(), s.end()), s.end());
            return s;
        };
        auto shares = [](const vector<OrientedReadId>& a, const AnchorWindow::HetAnchor& arm) {
            for (const auto& m : arm.members)
                if (binary_search(a.begin(), a.end(), m.orientedReadId)) return true;
            return false;
        };
        const vector<OrientedReadId> leadOids = oidSet(bubble.leadHom);
        const vector<OrientedReadId> homOids = oidSet(bubble.hom);
        bool bracketed = (bubble.leadHom.members.size() > 1) && (bubble.hom.members.size() > 1);
        if (!bracketed && failReason) {
            *failReason = "leadHom/hom too small: leadHom=" +
                std::to_string(bubble.leadHom.members.size()) + " hom=" +
                std::to_string(bubble.hom.members.size());
        }
        if (bracketed) {
            for (const AnchorWindow::HetAnchor& arm : bubble.alleles) {
                if (!shares(leadOids, arm) || !shares(homOids, arm)) {
                    bracketed = false;
                    if (failReason) *failReason = "an arm shares no member with leadHom/hom";
                    break;
                }
            }
        }
        if (!bracketed) return false;
    }

    outBubble = std::move(bubble);
    return true;
}

} // anonymous namespace

void Assembler::computeWindowProjectedAlignmentLeafSnarls(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const AlignOptions& alignOptions,
    uint64_t threadCount) const
{
    if (!markers) {
        cout << timestamp << "computeWindowProjectedAlignmentLeafSnarls: markers not available. Skipping." << endl;
        return;
    }
    const Reads& rds = getReads();
    const auto& mkrs = *markers;
    const uint64_t k = assemblerInfo->k;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
        if (threadCount == 0) threadCount = 1;
    }

    cout << timestamp << "computeWindowProjectedAlignmentLeafSnarls: processing "
         << anchorWindows.size() << " windows with " << threadCount << " threads." << endl;

    // Same sign convention as computeBaseAlignmentsAndStore
    // (AssemblerComputeAlignments.cpp:397-400): these AlignOptions fields are
    // already signed addends (match positive, mismatch/gaps negative), used
    // directly with NO extra negation.
    const int64_t dpMatchScore = alignOptions.overlapDpMatchScore;
    const int64_t dpMismatchScore = alignOptions.overlapDpMismatchScore;
    const int64_t dpGapOpen1 = alignOptions.overlapDpGapOpen1;
    const int64_t dpGapExtend1 = alignOptions.overlapDpGapExtend1;

    std::atomic<uint64_t> nextWindow{0};
    std::atomic<uint64_t> windowsProcessed{0};
    std::atomic<uint64_t> membersAligned{0};
    std::atomic<uint64_t> projectedAlignmentCalls{0};
    std::atomic<uint64_t> snarlsFound{0};
    std::atomic<uint64_t> pinGatherNs{0};      // readPins loop + sort/LIS per member
    std::atomic<uint64_t> projAlignNs{0};      // ProjectedAlignment construction + evidence extraction
    std::atomic<uint64_t> classifyNs{0};       // findLeafSnarlsFromSparseEvidence
    std::atomic<uint64_t> slowClassifyNs{0};   // DINARA_SPARSE_CLASSIFY_VALIDATE's dense cross-check only

    auto worker = [&]() {
        for (;;) {
            const uint64_t wi = nextWindow.fetch_add(1);
            if (wi >= anchorWindows.size()) break;
            const AnchorWindow& window = anchorWindows[wi];

            const OrientedReadId bbOid = window.backboneOrientedReadId;
            const auto bbJ = journeys[bbOid];

            if (window.backboneEnd <= window.backboneBegin + 1) continue;
            const uint32_t firstAnchorJP = window.backboneBegin;
            const uint32_t lastAnchorJP = window.backboneEnd - 1;
            if (firstAnchorJP >= bbJ.size() || lastAnchorJP >= bbJ.size()) continue;

            const uint32_t firstOrd = anchors.getOrdinal(bbJ[firstAnchorJP], bbOid);
            const uint32_t lastOrd = anchors.getOrdinal(bbJ[lastAnchorJP], bbOid);
            if (firstOrd == invalid<uint32_t> || lastOrd == invalid<uint32_t>) continue;

            const uint32_t windowBbBegin = mkrs[bbOid.getValue()][firstOrd].position;
            const uint32_t windowBbEnd = mkrs[bbOid.getValue()][lastOrd].position + uint32_t(k);
            if (windowBbEnd <= windowBbBegin) continue;

            vector<uint8_t> bbSeqVec(windowBbEnd);
            for (uint32_t i = 0; i < windowBbEnd; i++)
                bbSeqVec[i] = rds.getOrientedReadBase(bbOid, i).value;

            // Live per-read pins, same source as the ksw2/abPOA detectors:
            // every anchor a read shares with the backbone inside this
            // window, read from the anchor's own live member list.
            const auto tPinGather0 = std::chrono::steady_clock::now();
            unordered_map<uint64_t, vector<PaPin>> readPins;
            for (uint32_t jp = window.backboneBegin; jp < window.backboneEnd; jp++) {
                if (jp >= bbJ.size()) break;
                const Shasta2AnchorId aid = bbJ[jp];
                const uint32_t bbOrd = anchors.getOrdinal(aid, bbOid);
                if (bbOrd == invalid<uint32_t> || bbOrd >= mkrs[bbOid.getValue()].size()) continue;
                const uint32_t bbPos = mkrs[bbOid.getValue()][bbOrd].position;
                if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;

                const auto anchor = anchors[aid];
                for (const auto& info : anchor) {
                    const OrientedReadId cOid = info.orientedReadId;
                    if (cOid == bbOid) continue;
                    if (cOid.getValue() >= journeys.size()) continue;
                    const uint32_t cOrd = info.ordinal;
                    if (cOrd >= mkrs[cOid.getValue()].size()) continue;
                    const uint32_t cPos = mkrs[cOid.getValue()][cOrd].position;
                    readPins[cOid.getValue()].push_back(PaPin{bbOrd, bbPos, cOrd, cPos});
                }
            }
            pinGatherNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - tPinGather0).count()), std::memory_order_relaxed);

            vector<SparseMemberEvidence> evidence;
            evidence.reserve(readPins.size());

            for (auto& [orientedValue, pinsRaw] : readPins) {
                const OrientedReadId cOid = OrientedReadId::fromValue(ReadId(orientedValue));
                if (pinsRaw.size() < 2) continue;

                const auto tProjAlign0 = std::chrono::steady_clock::now();
                sort(pinsRaw.begin(), pinsRaw.end(),
                    [](const PaPin& a, const PaPin& b) {
                        return a.bbPos < b.bbPos || (a.bbPos == b.bbPos && a.cPos < b.cPos);
                    });
                paLisByCPos(pinsRaw);
                if (pinsRaw.size() < 2) continue;

                // ONE Alignment (marker-ordinal correspondences) for this
                // member's WHOLE pin list -- ProjectedAlignment loops over
                // every inter-anchor gap internally in its own constructor,
                // so this is one call per member, not one per segment.
                Alignment sparseAlign;
                sparseAlign.ordinals.reserve(pinsRaw.size());
                for (const PaPin& p : pinsRaw)
                    sparseAlign.ordinals.push_back({p.bbOrd, p.cOrd});

                const ProjectedAlignment pa(
                    *this, {bbOid, cOid}, sparseAlign,
                    ProjectedAlignment::Method::QuickRawSparse,
                    dpMatchScore, dpMismatchScore, dpGapOpen1, dpGapExtend1);
                projectedAlignmentCalls.fetch_add(1, std::memory_order_relaxed);

                SparseMemberEvidence ev;
                ev.oid = cOid;
                ev.bbCovBegin = max(pa.cigarRead0Start, windowBbBegin);
                ev.bbCovEnd = min(pa.cigarRead0End, windowBbEnd);
                if (ev.bbCovBegin >= ev.bbCovEnd) continue;

                for (const auto& mm : pa.sparseMismatches) {
                    if (mm.position0 < windowBbBegin || mm.position0 >= windowBbEnd) continue;
                    ev.snps.push_back(SparseSnp{mm.position0, mm.base1, mm.position1});
                }
                // deltaAtStart/deltaBreakpoints: the EXACT piecewise (readPos
                // - bbPos) offset, built by replaying every indel event the
                // aligner recorded, in CIGAR order (already non-decreasing in
                // position0 by construction -- a single linear walk). This
                // is exact everywhere in [bbCovBegin, bbCovEnd) outside a
                // deletionRange, unlike interpolating from the nearest
                // shared-anchor pin (which can be many kb away and gets
                // blocked by ANY unrelated indel in that wide span).
                ev.deltaAtStart = int64_t(pa.cigarRead1Start) - int64_t(pa.cigarRead0Start);
                ev.deltaBreakpoints.reserve(pa.sparseIndels.size());
                {
                    int64_t delta = ev.deltaAtStart;
                    for (const auto& ind : pa.sparseIndels) {
                        if (ind.op == 'D') {
                            delta -= int64_t(ind.length);
                            ev.deltaBreakpoints.push_back({ind.position0 + ind.length, delta});
                            const uint32_t a = max(ind.position0, windowBbBegin);
                            const uint32_t b = min(ind.position0 + ind.length, windowBbEnd);
                            if (a < b) ev.deletionRanges.push_back({a, b});
                        } else { // 'I': consumes no backbone position, delta applies from here on
                            delta += int64_t(ind.length);
                            ev.deltaBreakpoints.push_back({ind.position0, delta});
                        }
                    }
                }
                sort(ev.snps.begin(), ev.snps.end(),
                    [](const SparseSnp& a, const SparseSnp& b) { return a.bbPos < b.bbPos; });
                sort(ev.deletionRanges.begin(), ev.deletionRanges.end());

                const char* dumpOidEnv = getenv("DINARA_PROJALN_DUMP_OID");
                if (dumpOidEnv != nullptr && cOid == OrientedReadId(string(dumpOidEnv))) {
                    cout << "      dumpEvidence window=" << window.windowId << " oid=" << cOid
                         << " cigarRead0=[" << pa.cigarRead0Start << "," << pa.cigarRead0End << ")"
                         << " cigarRead1=[" << pa.cigarRead1Start << "," << pa.cigarRead1End << ")"
                         << " nSnps=" << ev.snps.size() << " nDeletions=" << ev.deletionRanges.size()
                         << " nPins=" << pinsRaw.size() << endl;
                    cout << "        mismatches:";
                    for (const auto& mm : pa.sparseMismatches)
                        cout << " (bb" << mm.position0 << "/mb" << mm.position1
                             << " " << "ACGT"[mm.base0 < 4 ? mm.base0 : 0]
                             << "->" << "ACGT"[mm.base1 < 4 ? mm.base1 : 0] << ")";
                    cout << endl;
                    cout << "        indels:";
                    for (const auto& ind : pa.sparseIndels)
                        cout << " (bb" << ind.position0 << "/mb" << ind.position1
                             << " " << ind.op << ind.length << ")";
                    cout << endl;
                }

                // Noise gate: OFF BY DEFAULT here, unlike the ksw2 detector.
                // Verified directly on real data that enabling it (same
                // win=100/maxS=5 as ksw2's HiFi tuning) REGRESSES this
                // detector from 26/27 exact matches (+2 new) down to 9/27:
                // a member (240-0) whose mismatches trace exactly onto known
                // real reference SNP positions (bb33615, bb34173, bb34674,
                // bb35339, ...) got its votes at MULTIPLE real sites excluded,
                // because this test region has an unusually dense cluster of
                // TRUE heterozygous SNPs (~10 within 300bp) that trips a
                // local-density heuristic meant to catch alignment-ambiguity
                // artifacts, not real closely-spaced variation. That
                // ambiguity is specific to ksw2's BANDED heuristic DP (where
                // an indel's exact placement genuinely is uncertain
                // read-to-read); ProjectedAlignment's astarpa2_simple is an
                // EXACT aligner, so this detector doesn't have the same
                // per-read gap-placement smearing problem the noise gate was
                // built to solve in the first place. Opt in with
                // DINARA_PROJALN_ENABLE_NOISE_GATE=1 to re-test if this ever
                // needs revisiting (e.g. against denser/noisier data).
                if (getenv("DINARA_PROJALN_ENABLE_NOISE_GATE") != nullptr) {
                    struct NoiseEvent { uint32_t pos; uint32_t len; int count; };
                    vector<NoiseEvent> events;
                    events.reserve(pa.sparseMismatches.size() + pa.sparseIndels.size());
                    for (const auto& mm : pa.sparseMismatches)
                        events.push_back({mm.position0, 1, 1});
                    for (const auto& ind : pa.sparseIndels)
                        events.push_back({ind.position0, ind.op == 'D' ? ind.length : 0, int(ind.length)});
                    sort(events.begin(), events.end(),
                        [](const NoiseEvent& a, const NoiseEvent& b) { return a.pos < b.pos; });
                    IpoaNoiseTracker noise(100, 5, ev.noisyRanges);
                    for (const NoiseEvent& e : events) noise.observe(e.pos, e.len, e.count);
                    noise.finish();

                    if (getenv("DINARA_PROJALN_NOISE_DEBUG") != nullptr && !ev.noisyRanges.empty()) {
                        cout << "      noiseProbe window=" << window.windowId << " oid=" << cOid
                             << " snps=" << ev.snps.size()
                             << " deletionRanges=" << ev.deletionRanges.size()
                             << " noisyRangeCount=" << ev.noisyRanges.size();
                        for (const auto& r : ev.noisyRanges)
                            cout << " [" << r.first << "," << r.second << ")";
                        cout << endl;
                    }
                }

                projAlignNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - tProjAlign0).count()), std::memory_order_relaxed);
                evidence.push_back(std::move(ev));
                membersAligned.fetch_add(1);
            }

            if (evidence.size() < 2) continue;

            if (const char* probeEnv = getenv("DINARA_PROJALN_LEAFSNARL_PROBE")) {
                const uint32_t probePos = uint32_t(atol(probeEnv));
                uint32_t nCov = 0, nMatch = 0, nMismatch = 0, nDeleted = 0, nNoisy = 0;
                for (const auto& ev : evidence) {
                    if (probePos < ev.bbCovBegin || probePos >= ev.bbCovEnd) continue;
                    nCov++;
                    bool isSnp = false, isDel = false, isNoisy2 = false;
                    for (const auto& snp : ev.snps) if (snp.bbPos == probePos) { isSnp = true; break; }
                    for (const auto& d : ev.deletionRanges)
                        if (probePos >= d.first && probePos < d.second) { isDel = true; break; }
                    for (const auto& nr : ev.noisyRanges)
                        if (probePos >= nr.first && probePos < nr.second) { isNoisy2 = true; break; }
                    if (isNoisy2) nNoisy++;
                    else if (isDel) nDeleted++;
                    else if (isSnp) nMismatch++;
                    else nMatch++;
                }
                cout << "    projAlnLeafSnarlProbe window=" << window.windowId
                     << " pos=" << probePos << " nMembersTotal=" << evidence.size()
                     << " covering=" << nCov << " match=" << nMatch
                     << " mismatch=" << nMismatch << " deleted=" << nDeleted
                     << " noisy=" << nNoisy << endl;
            }

            const vector<uint8_t> bbWindowBases(
                bbSeqVec.begin() + windowBbBegin, bbSeqVec.begin() + windowBbEnd);
            const auto tClassify0 = std::chrono::steady_clock::now();
            const auto snarls = findLeafSnarlsFromSparseEvidenceFast(
                bbOid, evidence, bbWindowBases, windowBbBegin,
                /*minSupport=*/6, /*minVaf=*/0.12,
                /*dropHomopolymer=*/true, /*dropRepeat=*/true);
            classifyNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - tClassify0).count()), std::memory_order_relaxed);
            snarlsFound.fetch_add(snarls.size());

            // Cross-check the sparse-native classifier against the original
            // dense-array one on the SAME evidence -- off by default (the
            // dense path duplicates the classifyNs cost this was built to
            // avoid), opt in with DINARA_SPARSE_CLASSIFY_VALIDATE=1.
            if (getenv("DINARA_SPARSE_CLASSIFY_VALIDATE") != nullptr) {
                const auto tSlow0 = std::chrono::steady_clock::now();
                const auto slowSnarls = findLeafSnarlsFromSparseEvidence(
                    bbOid, evidence, bbWindowBases, windowBbBegin,
                    /*minSupport=*/6, /*minVaf=*/0.12,
                    /*dropHomopolymer=*/true, /*dropRepeat=*/true);
                slowClassifyNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - tSlow0).count()), std::memory_order_relaxed);
                const bool match = leafSnarlsEqual(snarls, slowSnarls);
                cout << "    sparseClassifyValidate window=" << window.windowId
                     << " bb=" << bbOid << " fast=" << snarls.size()
                     << " slow=" << slowSnarls.size()
                     << " match=" << (match ? "YES" : "NO") << endl;
                if (!match) {
                    auto dumpSpans = [](const char* tag, const vector<LeafSnarl>& v) {
                        cout << "      " << tag << ":";
                        for (const auto& s : v) cout << " [" << s.start << "," << s.end << ")";
                        cout << endl;
                    };
                    dumpSpans("fast", snarls);
                    dumpSpans("slow", slowSnarls);
                }
            }

            for (const LeafSnarl& s : snarls) {
                const int64_t localStart = int64_t(s.start) - int64_t(windowBbBegin);
                const int64_t ctxBegin = max<int64_t>(0, localStart - 10);
                const int64_t ctxEnd = min<int64_t>(int64_t(bbWindowBases.size()), localStart + 12);
                string ctx;
                for (int64_t p = ctxBegin; p < ctxEnd; p++) {
                    const uint8_t b = bbWindowBases[size_t(p)];
                    ctx += "ACGTN"[b < 4 ? b : 4];
                }
                cout << "    windowProjAlnLeafSnarl window=" << window.windowId
                     << " bb=" << bbOid
                     << " [" << s.start << "," << s.end << ")"
                     << " ctx=" << ctx
                     << " alleles=" << s.alleles.size();
                for (const LeafSnarlAllele& al : s.alleles) {
                    cout << " {n=" << al.members.size() << " bases=";
                    for (uint8_t b : al.bases) cout << "ACGTN"[b < 4 ? b : 4];
                    cout << "}";
                }
                cout << endl;
            }

            // Verification pass: try to build a real AnchorWindow::HetBubble
            // (ref/alt arms + bracketing leadHom/hom, all members pinned to
            // exact read positions) for every SINGLE-COLUMN snarl, the same
            // structure emitHetBubblesFromProfiles/ksw2DetectHetBubblesInWindow
            // produce. Multi-column MNPs are skipped: HetAnchor::alleleBase
            // is a single base, so they aren't representable in this anchor
            // format without extending it (out of scope here). Logs the
            // result; does NOT write to window.hetBubbles (whichever engine
            // is actually selected for this run owns that). Off by default;
            // opt in with DINARA_PROJALN_BUILD_HETBUBBLES=1.
            if (getenv("DINARA_PROJALN_BUILD_HETBUBBLES") != nullptr) {
                unordered_map<uint64_t, const SparseMemberEvidence*> evidenceByOid;
                evidenceByOid.reserve(evidence.size());
                for (const auto& ev : evidence) evidenceByOid[ev.oid.getValue()] = &ev;

                uint64_t multiColumnSkipped = 0, bubblesBuilt = 0, bubblesFailed = 0;
                for (const LeafSnarl& s : snarls) {
                    if (s.end - s.start != 2) { multiColumnSkipped++; continue; }
                    AnchorWindow::HetBubble bubble;
                    string failReason;
                    if (!buildHetBubbleFromLeafSnarl(
                            s, evidenceByOid, bbSeqVec, bbOid, windowBbBegin,
                            /*minSupport=*/6, bubble, &failReason)) {
                        bubblesFailed++;
                        cout << "    projAlnHetBubbleFailed window=" << window.windowId
                             << " bb=" << bbOid << " [" << s.start << "," << s.end << ")"
                             << " reason=\"" << failReason << "\"" << endl;
                        continue;
                    }
                    bubblesBuilt++;
                    cout << "    projAlnHetBubble window=" << window.windowId
                         << " bb=" << bbOid
                         << " backboneOffset=" << bubble.backboneOffset
                         << " leadHomMembers=" << bubble.leadHom.members.size()
                         << " homMembers=" << bubble.hom.members.size()
                         << " arms=" << bubble.alleles.size();
                    for (const auto& arm : bubble.alleles)
                        cout << " {isRef=" << arm.isRef << " n=" << arm.members.size()
                             << " base=" << "ACGT"[arm.alleleBase < 4 ? arm.alleleBase : 0] << "}";
                    cout << endl;
                }
                if (multiColumnSkipped > 0 || bubblesBuilt > 0 || bubblesFailed > 0) {
                    cout << "    projAlnHetBubbleSummary window=" << window.windowId
                         << " built=" << bubblesBuilt << " failed=" << bubblesFailed
                         << " multiColumnSkipped=" << multiColumnSkipped << endl;
                }
            }

            windowsProcessed.fetch_add(1);
        }
    };

    vector<std::thread> threads;
    threads.reserve(threadCount);
    for (uint64_t t = 0; t < threadCount; t++) threads.emplace_back(worker);
    for (auto& th : threads) th.join();

    cout << timestamp << "computeWindowProjectedAlignmentLeafSnarls: processed " << windowsProcessed.load()
         << " windows, aligned " << membersAligned.load() << " members via "
         << projectedAlignmentCalls.load() << " ProjectedAlignment calls, found "
         << snarlsFound.load() << " leaf snarls." << endl;
    cout << timestamp << "  timing (summed across all windows/threads): "
         << "pin gather=" << (double(pinGatherNs.load()) / 1e6) << " ms, "
         << "ProjectedAlignment construct+extract=" << (double(projAlignNs.load()) / 1e6) << " ms, "
         << "leaf-snarl classify (fast)=" << (double(classifyNs.load()) / 1e6) << " ms";
    if (slowClassifyNs.load() > 0) {
        cout << ", leaf-snarl classify (dense cross-check)="
             << (double(slowClassifyNs.load()) / 1e6) << " ms";
    }
    cout << "." << endl;
}

// PRODUCTION engine: same per-member ProjectedAlignment sparse-evidence
// construction as computeWindowProjectedAlignmentLeafSnarls above (pin
// gather, ONE ProjectedAlignment call per member, SparseMemberEvidence with
// the exact deltaBreakpoints offset), but feeds findLeafSnarlsFromSparseEvidence's
// output through buildHetBubbleFromLeafSnarl and writes the result into
// window.hetBubbles for real -- a drop-in alternative to
// intervalPoaDetectHetBubblesAllWindows/ksw2DetectHetBubblesInWindow,
// selected via DINARA_HET_ENGINE=projaln (main.cpp). No debug dumps/probes/
// noise-gate here (the noise gate regressed this detector on real data, see
// the design note above); this is the lean production path.
uint32_t Assembler::projAlnDetectHetBubblesAllWindows(
    vector<AnchorWindow>& windows,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const AlignOptions& alignOptions,
    double hetMinVaf,
    uint64_t hetMinSupport,
    bool hetDropHomopolymer,
    bool hetDropRepeat,
    uint64_t threadCount,
    uint64_t& hetWindowsOut,
    uint64_t& totalBubblesOut,
    const vector<bool>* skipWindow) const
{
    hetWindowsOut = 0;
    totalBubblesOut = 0;
    if (!markers) {
        cout << timestamp << "projAlnDetectHetBubblesAllWindows: markers not available. Skipping." << endl;
        return 0;
    }
    const Reads& rds = getReads();
    const auto& mkrs = *markers;
    const uint64_t k = assemblerInfo->k;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
        if (threadCount == 0) threadCount = 1;
    }

    const int64_t dpMatchScore = alignOptions.overlapDpMatchScore;
    const int64_t dpMismatchScore = alignOptions.overlapDpMismatchScore;
    const int64_t dpGapOpen1 = alignOptions.overlapDpGapOpen1;
    const int64_t dpGapExtend1 = alignOptions.overlapDpGapExtend1;

    // hetMinSupport==0 is the "auto" sentinel (AssemblerOptions.hpp): derive
    // from the k-mer coverage histogram with the SAME rule
    // ksw2DetectHetBubblesInWindow/testAbpoaMultiSegmentMSA use (cut_bd=6,
    // cut_rate=0.7, n_hap=2) -- without this, the caller's raw 0 reaches
    // findLeafSnarlsFromSparseEvidence/buildHetBubbleFromLeafSnarl as a
    // literal minSupport of 0, which is not what "auto" means.
    uint64_t effectiveMinSupport = hetMinSupport;
    if (effectiveMinSupport == 0) {
        constexpr uint64_t cut_bd = 6, cut_rate_num = 7, cut_rate_den = 10, n_hap = 2;
        const uint64_t coverageHet = assemblerInfo.isOpen ?
            assemblerInfo->kmerDistributionInfo.coverageHet : invalid<uint64_t>;
        uint64_t base = 0;
        if (coverageHet != invalid<uint64_t> && coverageHet > 0) base = coverageHet / n_hap;
        uint64_t cc = (base * cut_rate_num) / cut_rate_den;
        if (cc < cut_bd) cc = cut_bd;
        effectiveMinSupport = cc;
    }

    // Shared, read-only across all windows/threads: per-oriented-read list
    // of alignmentIds touching it, used below to avoid a fresh
    // ProjectedAlignment call for any member whose pair with the backbone
    // already has a stored CIGAR from computeBaseAlignmentsAndStore --
    // same lookup pattern cigarDetectSnpsInWindow uses (AssemblerWindowCigarMSA.cpp).
    const auto& alnTable = getAlignmentTable();
    std::atomic<uint64_t> cigarReuseCount{0};
    std::atomic<uint64_t> freshAlignCount{0};
    std::atomic<uint64_t> pinGatherNs{0};    // per-window readPins loop (anchor scan)
    std::atomic<uint64_t> alnTableNs{0};     // per-window partnerReadToAlnIds build
    std::atomic<uint64_t> cigarReuseNs{0};   // findBestAlnForWindow + buildSparseEvidenceFromCigar
    std::atomic<uint64_t> freshAlignNs{0};   // sort/LIS + buildSparseEvidenceFromProjectedAlignment
    std::atomic<uint64_t> classifyNs{0};     // findLeafSnarlsFromSparseEvidenceFast
    std::atomic<uint64_t> bubbleBuildNs{0};  // buildHetBubbleFromLeafSnarl loop
    // DINARA_CIGAR_REUSE_VALIDATE base-correctness breakdown: exactAgree =
    // same bbPos, same altBase, same readPos (the two builders fully agree);
    // sameBbPosDisagree = same bbPos but DIFFERENT altBase/readPos -- this
    // is the real-bug signal (a base fetched from the wrong place), and
    // should be ~zero; onlyCigar/onlyRef = a snp recorded by one side but
    // not the other at that exact position -- the expected signature of
    // benign indel-placement ambiguity (the site is real, just landed a
    // base or two to one side depending which anchor set resolved the tie).
    std::atomic<uint64_t> exactAgreeCount{0};
    std::atomic<uint64_t> sameBbPosDisagreeCount{0};
    std::atomic<uint64_t> onlyCigarCount{0};
    std::atomic<uint64_t> onlyRefCount{0};

    std::atomic<uint64_t> nextWindow{0};
    std::atomic<uint64_t> hetWindows{0};
    std::atomic<uint64_t> totalBubbles{0};

    auto worker = [&]() {
        for (;;) {
            const uint64_t wi = nextWindow.fetch_add(1);
            if (wi >= windows.size()) break;
            AnchorWindow& window = windows[wi];
            window.hetBubbles.clear();
            if (skipWindow != nullptr && (*skipWindow)[wi]) continue;

            const OrientedReadId bbOid = window.backboneOrientedReadId;
            const auto bbJ = journeys[bbOid];

            if (window.backboneEnd <= window.backboneBegin + 1) continue;
            const uint32_t firstAnchorJP = window.backboneBegin;
            const uint32_t lastAnchorJP = window.backboneEnd - 1;
            if (firstAnchorJP >= bbJ.size() || lastAnchorJP >= bbJ.size()) continue;

            const uint32_t firstOrd = anchors.getOrdinal(bbJ[firstAnchorJP], bbOid);
            const uint32_t lastOrd = anchors.getOrdinal(bbJ[lastAnchorJP], bbOid);
            if (firstOrd == invalid<uint32_t> || lastOrd == invalid<uint32_t>) continue;

            const uint32_t windowBbBegin = mkrs[bbOid.getValue()][firstOrd].position;
            const uint32_t windowBbEnd = mkrs[bbOid.getValue()][lastOrd].position + uint32_t(k);
            if (windowBbEnd <= windowBbBegin) continue;

            vector<uint8_t> bbSeqVec(windowBbEnd);
            for (uint32_t i = 0; i < windowBbEnd; i++)
                bbSeqVec[i] = rds.getOrientedReadBase(bbOid, i).value;

            // Per-window: partner ReadId -> alignmentIds against the
            // backbone, for the CIGAR-store reuse path below. Mirrors
            // cigarDetectSnpsInWindow's own lookup exactly.
            const auto tAlnTable0 = std::chrono::steady_clock::now();
            const ReadId bbReadId = bbOid.getReadId();
            const uint32_t bbLen = uint32_t(rds.getRead(bbReadId).baseCount);
            unordered_map<uint64_t, vector<uint32_t>> partnerReadToAlnIds;
            if (bbOid.getValue() < alnTable.size()) {
                const auto bbAlnIds = alnTable[bbOid.getValue()];
                for (size_t ai = 0; ai < bbAlnIds.size(); ai++) {
                    const uint32_t alnId = bbAlnIds[ai];
                    const auto& ad = alignmentData[alnId];
                    if (ad.isDeleted0() && ad.isDeleted1()) continue;
                    if (ad.info.cigarOffset == uint32_t(-1) || ad.info.cigarTokenCount == 0) continue;
                    const ReadId partnerId = (ad.readIds[0] == bbReadId) ? ad.readIds[1] : ad.readIds[0];
                    partnerReadToAlnIds[uint64_t(partnerId)].push_back(alnId);
                }
            }
            alnTableNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - tAlnTable0).count()), std::memory_order_relaxed);
            auto findBestAlnForWindow = [&](ReadId partnerId) -> uint32_t {
                auto it = partnerReadToAlnIds.find(uint64_t(partnerId));
                if (it == partnerReadToAlnIds.end()) return uint32_t(-1);
                for (uint32_t alnId : it->second) {
                    const auto& ad = alignmentData[alnId];
                    const bool bbIsR0 = (ad.readIds[0] == bbReadId);
                    const uint32_t bs = bbIsR0 ? ad.qs : ad.ts;
                    const uint32_t be = bbIsR0 ? ad.qe : ad.te;
                    if (be > windowBbBegin && bs < windowBbEnd) return alnId;
                }
                return uint32_t(-1);
            };

            const auto tPinGather0 = std::chrono::steady_clock::now();
            unordered_map<uint64_t, vector<PaPin>> readPins;
            for (uint32_t jp = window.backboneBegin; jp < window.backboneEnd; jp++) {
                if (jp >= bbJ.size()) break;
                const Shasta2AnchorId aid = bbJ[jp];
                const uint32_t bbOrd = anchors.getOrdinal(aid, bbOid);
                if (bbOrd == invalid<uint32_t> || bbOrd >= mkrs[bbOid.getValue()].size()) continue;
                const uint32_t bbPos = mkrs[bbOid.getValue()][bbOrd].position;
                if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;

                const auto anchor = anchors[aid];
                for (const auto& info : anchor) {
                    const OrientedReadId cOid = info.orientedReadId;
                    if (cOid == bbOid) continue;
                    if (cOid.getValue() >= journeys.size()) continue;
                    const uint32_t cOrd = info.ordinal;
                    if (cOrd >= mkrs[cOid.getValue()].size()) continue;
                    const uint32_t cPos = mkrs[cOid.getValue()][cOrd].position;
                    readPins[cOid.getValue()].push_back(PaPin{bbOrd, bbPos, cOrd, cPos});
                }
            }
            pinGatherNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - tPinGather0).count()), std::memory_order_relaxed);

            vector<SparseMemberEvidence> evidence;
            evidence.reserve(readPins.size());

            for (auto& [orientedValue, pinsRaw] : readPins) {
                const OrientedReadId cOid = OrientedReadId::fromValue(ReadId(orientedValue));

                // Try the stored CIGAR first: computeBaseAlignmentsAndStore
                // already aligned this pair for ~92% of members, so reusing
                // it avoids a fully redundant fresh ProjectedAlignment call.
                // ON BY DEFAULT as of this change; opt out with
                // DINARA_CIGAR_REUSE_DISABLE=1 to fall back to fresh
                // ProjectedAlignment for every member (the pre-reuse
                // behavior), e.g. for A/B comparison.
                //
                // History: two real bugs were found and fixed in the
                // coordinate transform (an orientation-flip mismatch that
                // silently reused the wrong stored alignment -- guarded by
                // the bbFlip==targetFlip check below -- and a
                // breakpoint-placement bug under bbFlip). A third apparent
                // "regression" (small region: one window's snarlsFound
                // dropped 32->3) turned out to be a real base-identity bug
                // too, not structural ambiguity as first assumed: the
                // mismatch base fetch reused the POSITION-only collapsed
                // bbFlip/targetFlip flag to also index+complement the base
                // VALUE, which is only valid when the CIGAR's native walk
                // frame is already target-forward. Fixed by normalizing the
                // raw CIGAR-walk index to forward storage via needsRcTarget
                // FIRST, then complementing based on cOid's own strand
                // ALONE -- two independent steps, not one collapsed flag.
                // After that fix, the previously-regressed window matches
                // the fresh-ProjectedAlignment-only baseline exactly
                // (32 snarls found, 29 emitted, both ways), and spot-checked
                // large-region-only sites are clean, balanced biallelic
                // splits (not scattered noise) recovered by CIGAR reuse's
                // wider per-member coverage.
                if (getenv("DINARA_CIGAR_REUSE_DISABLE") == nullptr) {
                    const auto tCigar0 = std::chrono::steady_clock::now();
                    const uint32_t alnId = findBestAlnForWindow(cOid.getReadId());
                    SparseMemberEvidence ev;
                    const bool cigarOk = (alnId != uint32_t(-1)) && buildSparseEvidenceFromCigar(
                        *this, bbReadId, bbLen, bbOid, cOid, alnId,
                        windowBbBegin, windowBbEnd, ev);
                    cigarReuseNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - tCigar0).count()), std::memory_order_relaxed);
                    if (cigarOk) {
                            cigarReuseCount.fetch_add(1, std::memory_order_relaxed);

                            // Off by default (doubles the alignment work
                            // this path exists to avoid); opt in with
                            // DINARA_CIGAR_REUSE_VALIDATE=1 to cross-check
                            // against the trusted fresh-ProjectedAlignment
                            // path for every CIGAR-reuse hit.
                            if (getenv("DINARA_CIGAR_REUSE_VALIDATE") != nullptr) {
                                vector<PaPin> pinsCopy = pinsRaw;
                                SparseMemberEvidence evRef;
                                const bool haveRef = buildSparseEvidenceFromProjectedAlignment(
                                    *this, bbOid, cOid, pinsCopy, windowBbBegin, windowBbEnd,
                                    dpMatchScore, dpMismatchScore, dpGapOpen1, dpGapExtend1, evRef);
                                bool ok = haveRef && sparseEvidenceCoreEqual(ev, evRef);
                                // Also probe sparseColAt agreement at every
                                // 20th covered position (deltaBreakpoints
                                // placement differs between the two
                                // builders even when they agree on every
                                // query -- this is what actually matters).
                                if (ok) {
                                    const uint32_t lo = max(ev.bbCovBegin, evRef.bbCovBegin);
                                    const uint32_t hi = min(ev.bbCovEnd, evRef.bbCovEnd);
                                    for (uint32_t p = lo; p < hi; p += 20) {
                                        SparseColResult r1, r2;
                                        const bool ok1 = sparseColAt(ev, p, bbSeqVec[p], r1);
                                        const bool ok2 = sparseColAt(evRef, p, bbSeqVec[p], r2);
                                        if (ok1 != ok2 || (ok1 && (r1.readPos != r2.readPos || r1.readBase != r2.readBase))) {
                                            ok = false;
                                            break;
                                        }
                                    }
                                }
                                if (haveRef) {
                                    const uint32_t lo = max(ev.bbCovBegin, evRef.bbCovBegin);
                                    const uint32_t hi = min(ev.bbCovEnd, evRef.bbCovEnd);
                                    unordered_map<uint32_t, SparseSnp> refByPos;
                                    for (const auto& s : evRef.snps)
                                        if (s.bbPos >= lo && s.bbPos < hi) refByPos[s.bbPos] = s;
                                    unordered_set<uint32_t> seenCigarPos;
                                    for (const auto& s : ev.snps) {
                                        if (s.bbPos < lo || s.bbPos >= hi) continue;
                                        seenCigarPos.insert(s.bbPos);
                                        auto it2b = refByPos.find(s.bbPos);
                                        if (it2b == refByPos.end()) { onlyCigarCount.fetch_add(1); continue; }
                                        if (it2b->second.altBase == s.altBase && it2b->second.readPos == s.readPos) {
                                            exactAgreeCount.fetch_add(1);
                                        } else {
                                            sameBbPosDisagreeCount.fetch_add(1);
                                            if (getenv("DINARA_CIGAR_REUSE_DISAGREE_DUMP") != nullptr) {
                                                cout << "      sameBbPosDisagree window=" << window.windowId
                                                     << " bb=" << bbOid << " member=" << cOid
                                                     << " bbPos=" << s.bbPos
                                                     << " cigar=(alt" << int(s.altBase) << ",rp" << s.readPos << ")"
                                                     << " ref=(alt" << int(it2b->second.altBase) << ",rp" << it2b->second.readPos << ")"
                                                     << endl;
                                            }
                                        }
                                    }
                                    for (const auto& [pos, s] : refByPos)
                                        if (!seenCigarPos.count(pos)) onlyRefCount.fetch_add(1);
                                }
                                cout << "    cigarReuseValidate window=" << window.windowId
                                     << " bb=" << bbOid << " member=" << cOid
                                     << " haveRef=" << haveRef << " match=" << (ok ? "YES" : "NO") << endl;
                                if (!ok && getenv("DINARA_CIGAR_REUSE_DUMP") != nullptr) {
                                    const auto& adDbg = alignmentData[alnId];
                                    const auto it2 = partnerReadToAlnIds.find(uint64_t(cOid.getReadId()));
                                    cout << "      alnId=" << alnId
                                         << " readIds=[" << adDbg.readIds[0] << "," << adDbg.readIds[1] << "]"
                                         << " isSameStrand=" << adDbg.isSameStrand
                                         << " qs/qe=[" << adDbg.qs << "," << adDbg.qe << ")"
                                         << " ts/te=[" << adDbg.ts << "," << adDbg.te << ")"
                                         << " candidateAlnIdsForPartner=" << (it2 != partnerReadToAlnIds.end() ? it2->second.size() : 0)
                                         << " rawPinCount=" << pinsRaw.size() << endl;
                                    cout << "      cigar: cov=[" << ev.bbCovBegin << "," << ev.bbCovEnd
                                         << ") nSnps=" << ev.snps.size() << " nDel=" << ev.deletionRanges.size()
                                         << " nBreak=" << ev.deltaBreakpoints.size()
                                         << " deltaAtStart=" << ev.deltaAtStart << endl;
                                    cout << "      ref  : cov=[" << evRef.bbCovBegin << "," << evRef.bbCovEnd
                                         << ") nSnps=" << evRef.snps.size() << " nDel=" << evRef.deletionRanges.size()
                                         << " nBreak=" << evRef.deltaBreakpoints.size()
                                         << " deltaAtStart=" << evRef.deltaAtStart << endl;
                                    cout << "      cigar snps:";
                                    for (const auto& s : ev.snps) cout << " (bb" << s.bbPos << "/rp" << s.readPos << " alt" << int(s.altBase) << ")";
                                    cout << endl << "      ref   snps:";
                                    for (const auto& s : evRef.snps) cout << " (bb" << s.bbPos << "/rp" << s.readPos << " alt" << int(s.altBase) << ")";
                                    cout << endl;
                                    cout << "      cigar dels:";
                                    for (const auto& d : ev.deletionRanges) cout << " [" << d.first << "," << d.second << ")";
                                    cout << endl << "      ref   dels:";
                                    for (const auto& d : evRef.deletionRanges) cout << " [" << d.first << "," << d.second << ")";
                                    cout << endl;
                                    // Raw sequence dump around the divergent region, for
                                    // manual inspection: is this a repeat/homopolymer
                                    // stretch (expected ambiguity) or a clean region
                                    // (would indicate a real bug somewhere)?
                                    const uint32_t dumpLo = ev.bbCovBegin, dumpHi = ev.bbCovEnd;
                                    string bbSeq;
                                    for (uint32_t p = dumpLo; p < dumpHi; p++)
                                        bbSeq += "ACGT"[bbSeqVec[p] < 4 ? bbSeqVec[p] : 0];
                                    cout << "      backbone[" << dumpLo << "," << dumpHi << "): " << bbSeq << endl;
                                    const auto memberSeq = getReads().getRead(cOid.getReadId());
                                    string mSeq;
                                    const uint32_t mLen = uint32_t(memberSeq.baseCount);
                                    for (uint32_t p = 0; p < mLen; p++) {
                                        const uint8_t b = cOid.getStrand() == 0
                                            ? memberSeq[p].value
                                            : uint8_t((~memberSeq[mLen - 1 - p].value) & 3);
                                        mSeq += "ACGT"[b < 4 ? b : 0];
                                    }
                                    cout << "      member  [0," << mLen << "): " << mSeq << endl;
                                }
                            }

                            evidence.push_back(std::move(ev));
                            continue;
                    }
                }

                const auto tFresh0 = std::chrono::steady_clock::now();
                SparseMemberEvidence ev;
                if (!buildSparseEvidenceFromProjectedAlignment(
                        *this, bbOid, cOid, pinsRaw, windowBbBegin, windowBbEnd,
                        dpMatchScore, dpMismatchScore, dpGapOpen1, dpGapExtend1, ev)) {
                    freshAlignNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - tFresh0).count()), std::memory_order_relaxed);
                    continue;
                }
                freshAlignNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - tFresh0).count()), std::memory_order_relaxed);

                freshAlignCount.fetch_add(1, std::memory_order_relaxed);
                evidence.push_back(std::move(ev));
            }

            // A het/hom anchor can never legitimately contain the same raw
            // read twice. It normally can't here either -- but the SAME raw
            // read occasionally appears in readPins under BOTH orientations
            // (an anchor-graph property: e.g. a read touching this window's
            // backbone via two separate anchors that ended up on opposite
            // strands, plausible near a duplication/near-palindrome).
            // Previously this was harmless in practice because usually only
            // one orientation's fresh-ProjectedAlignment pin set actually
            // survived (enough colinear pins, real coverage); CIGAR reuse
            // succeeds unconditionally whenever a correctly-oriented stored
            // alignment exists, so both orientations can now survive
            // together -- verified: this is what caused a real
            // "anchor contains the same ReadId on both strands" export
            // failure on the large test region. Keep only the
            // wider-coverage orientation per raw read.
            if (evidence.size() >= 2) {
                unordered_map<uint64_t, size_t> bestByReadId; // ReadId -> index into evidence
                for (size_t i = 0; i < evidence.size(); i++) {
                    const uint64_t rid = uint64_t(evidence[i].oid.getReadId());
                    auto it = bestByReadId.find(rid);
                    if (it == bestByReadId.end()) {
                        bestByReadId[rid] = i;
                    } else if ((evidence[i].bbCovEnd - evidence[i].bbCovBegin) >
                               (evidence[it->second].bbCovEnd - evidence[it->second].bbCovBegin)) {
                        it->second = i;
                    }
                }
                if (bestByReadId.size() < evidence.size()) {
                    vector<SparseMemberEvidence> deduped;
                    deduped.reserve(bestByReadId.size());
                    for (const auto& [rid, idx] : bestByReadId) deduped.push_back(std::move(evidence[idx]));
                    evidence = std::move(deduped);
                }
            }

            if (evidence.size() < 2) continue;

            if (const char* probeEnv = getenv("DINARA_CIGAR_REUSE_PROBE_POS")) {
                const uint32_t probePos = uint32_t(atol(probeEnv));
                if (probePos >= windowBbBegin && probePos < windowBbEnd) {
                    cout << "    cigarReuseProbe window=" << window.windowId
                         << " bb=" << bbOid << " pos=" << probePos << " bbBase="
                         << "ACGT"[bbSeqVec[probePos] < 4 ? bbSeqVec[probePos] : 0] << endl;
                    for (const auto& ev : evidence) {
                        if (probePos < ev.bbCovBegin || probePos >= ev.bbCovEnd) {
                            cout << "      " << ev.oid << " NotCovered" << endl;
                            continue;
                        }
                        bool isDel = false;
                        for (const auto& d : ev.deletionRanges)
                            if (probePos >= d.first && probePos < d.second) { isDel = true; break; }
                        if (isDel) { cout << "      " << ev.oid << " Deleted" << endl; continue; }
                        uint8_t altBase = 0xff;
                        for (const auto& s : ev.snps) if (s.bbPos == probePos) { altBase = s.altBase; break; }
                        SparseColResult r;
                        const bool ok = sparseColAt(ev, probePos, bbSeqVec[probePos], r);
                        cout << "      " << ev.oid
                             << " snpAlt=" << (altBase == 0xff ? string("-") : string(1, "ACGT"[altBase < 4 ? altBase : 0]))
                             << " colAt=" << (ok ? string(1, "ACGT"[r.readBase < 4 ? r.readBase : 0]) : string("FAIL"))
                             << " readPos=" << (ok ? int64_t(r.readPos) : -1) << endl;
                    }
                }
            }

            const vector<uint8_t> bbWindowBases(
                bbSeqVec.begin() + windowBbBegin, bbSeqVec.begin() + windowBbEnd);
            const auto tClassify0 = std::chrono::steady_clock::now();
            const auto snarls = findLeafSnarlsFromSparseEvidenceFast(
                bbOid, evidence, bbWindowBases, windowBbBegin,
                effectiveMinSupport, hetMinVaf, hetDropHomopolymer, hetDropRepeat);
            classifyNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - tClassify0).count()), std::memory_order_relaxed);
            if (snarls.empty()) continue;

            const auto tBubble0 = std::chrono::steady_clock::now();
            unordered_map<uint64_t, const SparseMemberEvidence*> evidenceByOid;
            evidenceByOid.reserve(evidence.size());
            for (const auto& ev : evidence) evidenceByOid[ev.oid.getValue()] = &ev;

            uint64_t emitted = 0;
            for (const LeafSnarl& s : snarls) {
                if (s.end - s.start != 2) continue; // multi-column MNP: not representable yet
                AnchorWindow::HetBubble bubble;
                if (!buildHetBubbleFromLeafSnarl(
                        s, evidenceByOid, bbSeqVec, bbOid, windowBbBegin,
                        effectiveMinSupport, bubble, nullptr)) continue;
                window.hetBubbles.push_back(std::move(bubble));
                emitted++;
            }
            bubbleBuildNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - tBubble0).count()), std::memory_order_relaxed);
            if (getenv("DINARA_CIGAR_REUSE_DEBUG") != nullptr) {
                cout << "    projAlnWindowBubbles window=" << window.windowId
                     << " bb=" << bbOid << " snarlsFound=" << snarls.size()
                     << " emitted=" << emitted << " evidenceCount=" << evidence.size() << ":";
                for (const auto& s : snarls) cout << " [" << s.start << "," << s.end << ")";
                cout << endl;
            }
            if (emitted > 0) {
                hetWindows.fetch_add(1);
                totalBubbles.fetch_add(emitted);
            }
        }
    };

    vector<std::thread> threads;
    threads.reserve(threadCount);
    for (uint64_t t = 0; t < threadCount; t++) threads.emplace_back(worker);
    for (auto& th : threads) th.join();

    hetWindowsOut = hetWindows.load();
    totalBubblesOut = totalBubbles.load();
    cout << timestamp << "projAlnDetectHetBubblesAllWindows: "
         << cigarReuseCount.load() << " members via stored CIGAR reuse, "
         << freshAlignCount.load() << " via fresh ProjectedAlignment." << endl;
    cout << timestamp << "  timing (summed across all windows/threads): "
         << "pin gather=" << (double(pinGatherNs.load()) / 1e6) << " ms, "
         << "alnTable build=" << (double(alnTableNs.load()) / 1e6) << " ms, "
         << "CIGAR reuse (lookup+build)=" << (double(cigarReuseNs.load()) / 1e6) << " ms, "
         << "fresh align (sort/LIS+ProjectedAlignment)=" << (double(freshAlignNs.load()) / 1e6) << " ms, "
         << "classify=" << (double(classifyNs.load()) / 1e6) << " ms, "
         << "bubble build=" << (double(bubbleBuildNs.load()) / 1e6) << " ms." << endl;
    if (getenv("DINARA_CIGAR_REUSE_VALIDATE") != nullptr) {
        const uint64_t exact = exactAgreeCount.load();
        const uint64_t disagree = sameBbPosDisagreeCount.load();
        const uint64_t onlyC = onlyCigarCount.load();
        const uint64_t onlyR = onlyRefCount.load();
        cout << timestamp << "  base-correctness breakdown: exactAgree=" << exact
             << " sameBbPosDISAGREE=" << disagree << " (real-bug signal, should be ~0)"
             << " onlyCigar=" << onlyC << " onlyRef=" << onlyR
             << " (both expected under benign indel-placement ambiguity)" << endl;
    }
    return uint32_t(totalBubblesOut);
}
