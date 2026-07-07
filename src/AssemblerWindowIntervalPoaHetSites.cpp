/// @file AssemblerWindowIntervalPoaHetSites.cpp
/// @brief Per-window het-SNP detection using one small POA per consecutive-anchor
///        interval (Shasta2 LocalAssembly7 model), instead of a single
///        whole-window POA.
///
/// Motivation
/// ----------
/// testAbpoaMultiSegmentMSA builds ONE abPOA graph per window and folds every
/// read into it. That graph grows as reads are added, so a read's segment DP
/// cost scales with the accumulated graph node span (not local depth), and the
/// whole window runs on a single thread (setupLoadBalancing(windowEnd, 1)). On a
/// 100kb / 200-read window this is the dominant cost of the assembler.
///
/// Shasta2's LocalAssembly7 avoids this: it never builds a whole-window graph.
/// It assembles between consecutive anchors, so each POA sees only the reads
/// spanning one short gap and a tiny graph. This engine ports that idea to
/// het detection:
///
///   1. Tile the window at its backbone anchors -> intervals [a_i, a_{i+1}).
///   2. For each interval, gather the backbone segment (row 0) + each spanning
///      read's bracketed subsequence, and run one abPOA MSA (abpoa_msa) with a
///      thread-local, per-interval-reset handle.
///   3. Map each interval's MSA columns to ABSOLUTE backbone offsets (row 0 is
///      the backbone) and append aligned columns / SNPs / deletions into that
///      read's window-global KwMemberProfile.
///   4. The shared emitHetBubblesFromProfiles tail emits the identical
///      AnchorWindow::hetBubbles.
///
/// vs the other engines:
///   ksw2    : star alignment (each read vs backbone) -> alt alleles smear
///             across columns near indels. This engine aligns reads to EACH
///             OTHER within an interval, so a shared variant lands in one column.
///   whole-window abPOA : same POA quality, but this keeps graphs tiny and the
///             intervals independent (parallelizable), so it does not pay the
///             graph-growth / single-thread cost.
///
/// The window/read enumeration (live-anchor pins + LIS) is identical to the
/// ksw2 path so all engines see the same member set.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "AssemblerOptions.hpp"
#include "Base.hpp"
#include "PhasingKmeansTypes.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "WindowHetProfiles.hpp"
#include "globalMsa.hpp"
#include "invalid.hpp"

#include <abpoa.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <vector>

using namespace std;
using namespace dinara;

// Detect het bubbles in one window with per-interval POA. Signature mirrors
// ksw2DetectHetBubblesInWindow so main.cpp can select engines
// interchangeably.
uint32_t Assembler::intervalPoaDetectHetBubblesInWindow(
    AnchorWindow& window,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const AlignOptions& /*alignOptions*/,
    double hetMinVaf,
    uint64_t hetMinSupport,
    bool hetDropHomopolymer,
    bool hetDropRepeat) const
{
    const Reads& rds = getReads();
    const auto& mkrs = *markers;
    const uint64_t k = assemblerInfo->k;
    const OrientedReadId bbOid = window.backboneOrientedReadId;
    const ReadId bbReadId = bbOid.getReadId();
    const uint32_t bbLen = uint32_t(rds.getRead(bbReadId).baseCount);
    const auto bbJ = journeys[bbOid];

    window.hetBubbles.clear();
    if (window.backboneEnd <= window.backboneBegin + 1) return 0;

    // Backbone coordinate frame -- identical to the abPOA/ksw2 paths:
    // START of the first anchor to END of the last anchor (marker.position, no
    // half-k), so backboneOffset = bbPos - windowBbBegin matches the abPOA frame
    // exactly.
    const uint32_t firstAnchorJP = window.backboneBegin;
    const uint32_t lastAnchorJP = window.backboneEnd - 1;
    if (firstAnchorJP >= bbJ.size() || lastAnchorJP >= bbJ.size()) return 0;

    const uint32_t firstOrd = anchors.getOrdinal(bbJ[firstAnchorJP], bbOid);
    const uint32_t lastOrd = anchors.getOrdinal(bbJ[lastAnchorJP], bbOid);
    if (firstOrd == invalid<uint32_t> || lastOrd == invalid<uint32_t>) return 0;

    const uint32_t windowBbBegin = mkrs[bbOid.getValue()][firstOrd].position;
    const uint32_t windowBbEnd = mkrs[bbOid.getValue()][lastOrd].position + uint32_t(k);
    if (windowBbEnd <= windowBbBegin) return 0;

    vector<uint8_t> bbSeqVec(bbLen);
    for (uint32_t i = 0; i < bbLen; i++)
        bbSeqVec[i] = rds.getOrientedReadBase(bbOid, i).value;

    const bool debug = (getenv("DINARA_HET_DEBUG") != nullptr);

    // --- Live-anchor pin enumeration (identical to ksw2 path) ------
    // Each backbone anchor a read shares inside the window becomes a pin
    // (backbone base position, read base position). The LIS over read positions
    // drops off-diagonal (repeat) pins. The kept pins define the read's segment
    // boundaries: between consecutive pins is one interval to POA.
    struct Pin { uint32_t bbPos; uint32_t cPos; };
    unordered_map<uint64_t, vector<Pin>> readPins;

    for (uint32_t jp = window.backboneBegin; jp < window.backboneEnd; jp++) {
        if (jp >= bbJ.size()) break;
        const Shasta2AnchorId aid = bbJ[jp];
        const uint32_t bbOrd = anchors.getOrdinal(aid, bbOid);
        if (bbOrd == invalid<uint32_t> || bbOrd >= mkrs[bbOid.getValue()].size())
            continue;
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
            readPins[cOid.getValue()].push_back(Pin{bbPos, cPos});
        }
    }

    auto lisByCPos = [](vector<Pin>& pins) {
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
        vector<Pin> kept;
        for (int32_t c = int32_t(tails.back()); c >= 0; c = pred[uint32_t(c)])
            kept.push_back(pins[uint32_t(c)]);
        reverse(kept.begin(), kept.end());
        pins.swap(kept);
    };

    // Collect the ordered, colinear pin list for every member that keeps >=2
    // pins. These define both the intervals and each member's segment endpoints.
    struct Member { OrientedReadId oid; vector<Pin> pins; };
    vector<Member> members;
    members.reserve(readPins.size());
    for (auto& [orientedValue, pinsRaw] : readPins) {
        const OrientedReadId cOid = OrientedReadId::fromValue(ReadId(orientedValue));
        if (pinsRaw.size() < 2) continue;
        sort(pinsRaw.begin(), pinsRaw.end(),
            [](const Pin& a, const Pin& b) {
                return a.bbPos < b.bbPos || (a.bbPos == b.bbPos && a.cPos < b.cPos);
            });
        lisByCPos(pinsRaw);
        if (pinsRaw.size() < 2) continue;
        members.push_back({cOid, std::move(pinsRaw)});
    }
    if (members.size() < 2) return 0;   // need >=2 reads for a het

    // --- Build the set of interval breakpoints ------------------------------
    // A breakpoint is a backbone position where some read has an anchor pin.
    // Sorted unique pin positions partition the window into intervals; the reads
    // spanning [b_i, b_{i+1}) are POA'd together. Using the union of all pins
    // (not just the backbone's own consecutive anchors) keeps intervals short
    // even where the backbone is sparsely anchored, and guarantees every
    // per-read segment endpoint coincides with an interval boundary.
    vector<uint32_t> breakpoints;
    for (const Member& m : members)
        for (const Pin& p : m.pins)
            breakpoints.push_back(p.bbPos);
    sort(breakpoints.begin(), breakpoints.end());
    breakpoints.erase(unique(breakpoints.begin(), breakpoints.end()), breakpoints.end());
    if (breakpoints.size() < 2) return 0;

    // Map a backbone position to its breakpoint index (only pin positions map).
    unordered_map<uint32_t, uint32_t> bpIndex;
    bpIndex.reserve(breakpoints.size() * 2);
    for (uint32_t i = 0; i < breakpoints.size(); i++) bpIndex[breakpoints[i]] = i;

    // Per-read profile accumulator, keyed by oid value. Filled interval by
    // interval, then flushed into the profiles vector for the emit tail.
    struct Accum {
        OrientedReadId oid;
        vector<KwAlignedCol> alignedCols;
        vector<KwSnp> snps;
        vector<pair<uint32_t, uint32_t>> deletionRanges;
        int64_t firstBb = -1, lastBb = -1;
    };
    unordered_map<uint64_t, Accum> accums;
    accums.reserve(members.size() * 2);

    // For each member, index its pins by breakpoint so a given interval can
    // pull the read's [cPos_i, cPos_{i+1}) subsequence directly. pinAtBp[oid]
    // maps breakpoint index -> read base position at that anchor.
    // We iterate intervals as the outer loop (so each interval is one POA and is
    // the natural unit of parallelism), and for each interval gather the members
    // that have BOTH endpoints pinned (i.e. actually span the interval).
    struct MemberBp { uint32_t oidValue; uint32_t cPos; };
    // intervalMembers[i] = members pinned at breakpoint i, with their read pos.
    vector<vector<MemberBp>> atBreakpoint(breakpoints.size());
    for (const Member& m : members) {
        for (const Pin& p : m.pins) {
            auto it = bpIndex.find(p.bbPos);
            if (it != bpIndex.end())
                atBreakpoint[it->second].push_back({m.oid.getValue(), p.cPos});
        }
    }
    for (auto& v : atBreakpoint)
        sort(v.begin(), v.end(),
            [](const MemberBp& a, const MemberBp& b) { return a.oidValue < b.oidValue; });

    // One-sided read support (semi-global). Reads pinned at only ONE interval
    // boundary are still placed into the interval POA; their bases are clipped to
    // the aligned footprint during column extraction so the free end does not
    // force a global end-gap. This recovers the ~19% of read placements the
    // strict both-anchors gate drops (up to 62% in the largest windows).
    //
    // ON BY DEFAULT. Without it, a read that spans a het locus but is pinned at
    // only one of the interval's two backbone anchors is dropped, leaving a
    // COVERAGE HOLE: the read's cov range still spans the SNP, so the emit tail
    // miscounts it as ref support, yet no aligned column exists to build the arm
    // from -- yielding coverage-1 ref arms and, in windows seeded by the
    // opposite haplotype, whole het sites with no recoverable ref allele. With
    // one-sided placement these reads get real columns, restoring valid
    // multi-member ref/alt arms. Set DINARA_IPOA_ONESIDED=0 to disable.
    const char* oneSidedEnv = getenv("DINARA_IPOA_ONESIDED");
    const bool oneSidedEnabled = (oneSidedEnv == nullptr) || (oneSidedEnv[0] != '0');

    uint32_t intervalsRun = 0;
    uint64_t twoSidedTotal = 0;    // reads pinned at both interval boundaries
    uint64_t oneSidedTotal = 0;    // reads pinned at exactly one boundary
    uint64_t oneSidedPlaced = 0;   // one-sided reads actually POA'd (when enabled)
    uint32_t intervalsWithOneSided = 0;

    // --- Per-interval POA ---------------------------------------------------
    // Each interval is an independent abPOA MSA. Windows are the unit of
    // parallelism (see main.cpp), so this loop runs serially within a window.
    for (uint32_t bi = 0; bi + 1 < breakpoints.size(); bi++) {
        const uint32_t bbBegin = breakpoints[bi];
        const uint32_t bbEnd = breakpoints[bi + 1];
        if (bbEnd <= bbBegin) continue;
        if (bbEnd <= windowBbBegin || bbBegin >= windowBbEnd) continue;

        // Members spanning this interval = those pinned at BOTH bi and bi+1.
        // Intersect the two sorted breakpoint member lists. The <-/>- branches
        // are reads present at only one boundary (one-sided); with oneSidedEnabled
        // they are captured for semi-global placement, else counted only.
        //   side: 0 = both, 1 = left-only (pinned at bi), 2 = right-only (bi+1).
        // For one-sided reads the free-side subsequence is bounded to ~segLen
        // bases past the pinned anchor (clamped to the read end) so the kOV
        // alignment stays cheap and local.
        struct SpanMember { OrientedReadId oid; uint32_t cBegin; uint32_t cEnd; uint8_t side; };
        vector<SpanMember> spanMembers;
        uint64_t oneSidedHere = 0;
        const uint32_t segLenApprox = bbEnd - bbBegin;
        auto readLenOf = [&](uint32_t oidValue) -> uint32_t {
            const OrientedReadId o = OrientedReadId::fromValue(ReadId(oidValue));
            return uint32_t(rds.getRead(o.getReadId()).baseCount);
        };
        {
            const auto& lo = atBreakpoint[bi];
            const auto& hi = atBreakpoint[bi + 1];
            size_t a = 0, b = 0;
            auto addLeftOnly = [&](const MemberBp& e) {
                oneSidedHere++;
                if (!oneSidedEnabled) return;
                // Pinned at left anchor (cPos), free on the right: take up to
                // segLenApprox bases forward, clamped to the read end.
                const uint32_t rl = readLenOf(e.oidValue);
                const uint32_t cB = e.cPos;
                if (cB >= rl) return;
                const uint32_t cE = min(rl, cB + max<uint32_t>(segLenApprox, 1));
                if (cE > cB)
                    spanMembers.push_back(
                        {OrientedReadId::fromValue(ReadId(e.oidValue)), cB, cE, 1});
            };
            auto addRightOnly = [&](const MemberBp& e) {
                oneSidedHere++;
                if (!oneSidedEnabled) return;
                // Pinned at right anchor (cPos = read base AT bbEnd), free on the
                // left: take up to segLenApprox bases backward, clamped to 0.
                const uint32_t cE = e.cPos;
                if (cE == 0) return;
                const uint32_t back = max<uint32_t>(segLenApprox, 1);
                const uint32_t cB = (cE > back) ? (cE - back) : 0;
                if (cE > cB)
                    spanMembers.push_back(
                        {OrientedReadId::fromValue(ReadId(e.oidValue)), cB, cE, 2});
            };
            while (a < lo.size() && b < hi.size()) {
                if (lo[a].oidValue < hi[b].oidValue) { addLeftOnly(lo[a]); a++; }
                else if (lo[a].oidValue > hi[b].oidValue) { addRightOnly(hi[b]); b++; }
                else {
                    const uint32_t cB = lo[a].cPos;
                    const uint32_t cE = hi[b].cPos;
                    if (cE > cB)
                        spanMembers.push_back(
                            {OrientedReadId::fromValue(ReadId(lo[a].oidValue)), cB, cE, 0});
                    a++; b++;
                }
            }
            while (a < lo.size()) { addLeftOnly(lo[a]); a++; }
            while (b < hi.size()) { addRightOnly(hi[b]); b++; }
        }
        uint64_t twoSidedHere = 0;
        for (const auto& sm : spanMembers) if (sm.side == 0) twoSidedHere++;
        twoSidedTotal += twoSidedHere;
        oneSidedTotal += oneSidedHere;
        oneSidedPlaced += (spanMembers.size() - twoSidedHere);
        if (oneSidedHere > 0) intervalsWithOneSided++;
        if (spanMembers.empty()) continue;

        // Interval backbone segment [bbBegin, bbEnd), aligned by the anchor
        // marker so member subsequences [cBegin, cEnd) hit the same locus.
        const uint32_t segLen = bbEnd - bbBegin;

        // Order rows: backbone (row 0), then two-sided reads, then one-sided.
        // abpoa_msa keeps input order (sort_input_seq=0), so the RC-MSA rows
        // come back in this same order and rowMembers[si-1] is the read for
        // row si. All rows are aligned in global mode; one-sided reads (kOV
        // semantics) are only used when DINARA_IPOA_ONESIDED is set and are
        // clipped to their aligned footprint during column extraction.
        vector<const SpanMember*> ordered;
        ordered.reserve(spanMembers.size());
        for (const SpanMember& sm : spanMembers) if (sm.side == 0) ordered.push_back(&sm);
        const uint32_t nTwoSided = uint32_t(ordered.size());
        for (const SpanMember& sm : spanMembers) if (sm.side != 0) ordered.push_back(&sm);
        if (nTwoSided < 2) continue;   // need >=2 two-sided reads to anchor the POA

        // Fill a 0..3 base-code buffer for one member's [cB, cE) subsequence.
        auto toCodes = [&](const OrientedReadId o, uint32_t cB, uint32_t cE,
                           vector<uint8_t>& out) {
            out.resize(cE - cB);
            for (uint32_t i = cB; i < cE; i++)
                out[i - cB] = uint8_t(rds.getOrientedReadBase(o, i).value & 3);
        };

        // Per-interval abPOA MSA. One graph per interval keeps each POA tiny and
        // independent (the whole point of the interval engine); abPOA is the same
        // aligner the whole-window path uses, with long-read-tuned affine scoring
        // (match=2, mismatch=4, gap 4/2) and adaptive banding (wb=10). The handle
        // is thread-local and reset internally by abpoa_msa on each call, so it
        // is reused across intervals on this thread with no per-interval alloc.
        // Row order follows input order (sort_input_seq=0): row 0 = backbone,
        // rows 1..N = ordered[] in the same order, matching the previous spoa path.
        struct AbHandle {
            abpoa_t* ab = nullptr;
            abpoa_para_t* abpt = nullptr;
            AbHandle() {
                ab = abpoa_init();
                abpt = abpoa_init_para();
                abpt->align_mode = ABPOA_GLOBAL_MODE;
                abpt->gap_mode = ABPOA_AFFINE_GAP;
                abpt->match = 2;
                abpt->mismatch = 4;
                abpt->gap_open1 = 4;
                abpt->gap_ext1 = 2;
                abpt->gap_open2 = 0;   // single-piece affine gap
                abpt->gap_ext2 = 0;
                abpt->wb = 10;         // adaptive band; <0 disables banding
                abpt->wf = 0.01;
                abpt->disable_seeding = 1;   // segmentation is driven by intervals
                abpt->progressive_poa = 0;
                abpt->sort_input_seq = 0;    // keep row 0 = backbone
                abpt->out_msa = 1;
                abpt->out_cons = 0;
                abpt->ret_cigar = 1;   // required: RC-MSA row extraction needs it
                abpoa_post_set_para(abpt);
            }
            ~AbHandle() { if (ab) abpoa_free(ab); if (abpt) abpoa_free_para(abpt); }
        };
        static thread_local AbHandle ah;

        // Assemble the input sequence set: row 0 = backbone segment, then each
        // ordered member's subsequence, all as 0..3 code arrays for abpoa_msa.
        // Reused thread-local scratch buffers avoid per-interval reallocation.
        const int nSeq = int(ordered.size()) + 1;
        static thread_local vector<vector<uint8_t>> seqStore;
        static thread_local vector<uint8_t*> seqPtrs;
        static thread_local vector<int> seqLens;
        seqStore.resize(nSeq);
        seqPtrs.resize(nSeq);
        seqLens.resize(nSeq);

        {
            auto& bb = seqStore[0];
            bb.resize(segLen);
            for (uint32_t i = 0; i < segLen; i++) bb[i] = uint8_t(bbSeqVec[bbBegin + i] & 3);
            seqPtrs[0] = bb.data();
            seqLens[0] = int(segLen);
        }
        bool addFailed = false;
        for (uint32_t oi = 0; oi < ordered.size(); oi++) {
            const SpanMember& sm = *ordered[oi];
            auto& codes = seqStore[oi + 1];
            toCodes(sm.oid, sm.cBegin, sm.cEnd, codes);
            if (codes.empty()) { addFailed = true; break; }
            seqPtrs[oi + 1] = codes.data();
            seqLens[oi + 1] = int(codes.size());
        }
        if (addFailed) continue;

        // The thread-local handle is reused across intervals, so reset it before
        // each MSA (abpoa_msa only auto-resets when the seq store is empty).
        int maxLen = 0;
        for (int r = 0; r < nSeq; r++) if (seqLens[r] > maxLen) maxLen = seqLens[r];
        abpoa_reset(ah.ab, ah.abpt, maxLen > 0 ? maxLen : 1);
        abpoa_msa(ah.ab, ah.abpt, nSeq, nullptr, seqLens.data(),
                  seqPtrs.data(), nullptr, nullptr);
        const abpoa_cons_t* abc = ah.ab->abc;
        if (abc == nullptr || abc->msa_len <= 0 || abc->n_seq != nSeq) continue;
        ++intervalsRun;

        // abPOA base codes: 0..3 = ACGT, 4 = N, 5 = gap (abpt->m). Map directly
        // to AlignedBase integer values: 0..3 stay, gap(5)/N(4) -> 4 (AlignedBase
        // gap). Row 0 = backbone, so downstream isGap()/SNP handling matches the
        // previous spoa path.
        const uint64_t ncols = uint64_t(abc->msa_len);
        vector<vector<AlignedBase>> alignment(nSeq);
        for (int r = 0; r < nSeq; r++) {
            alignment[r].resize(ncols);
            const uint8_t* row = abc->msa_base[r];
            for (uint64_t c = 0; c < ncols; c++) {
                const uint8_t code = row[c];
                alignment[r][c] = AlignedBase::fromInteger(uint8_t(code < 4 ? code : 4));
            }
        }

        // Column -> absolute backbone position. Row 0 is the backbone segment;
        // walk it, assigning bbBegin, bbBegin+1, ... to each non-gap column.
        // Backbone-gap columns (insertions) get bbPos = -1 and are ignored for
        // SNP calling (they cannot be substitution columns).
        vector<int64_t> colBbPos(ncols, -1);
        {
            uint32_t bbPos = bbBegin;
            for (uint64_t col = 0; col < ncols; col++) {
                if (!alignment[0][col].isGap()) {
                    colBbPos[col] = int64_t(bbPos);
                    bbPos++;
                }
            }
        }

        // Each member row -> aligned columns / SNPs / deletions, accumulated
        // into the read's window-global profile. Read positions are recovered by
        // walking the member's subsequence in parallel with the columns.
        //
        // Footprint clipping: for EVERY row we only record between the row's
        // first and last non-gap column. For two-sided reads this is the whole
        // interval; for one-sided (kOV) reads it restricts recording to where
        // the read actually aligned, so the free end produces no phantom
        // deletions or aligned columns beyond the read's real extent.
        for (uint32_t si = 1; si < alignment.size(); si++) {
            const SpanMember& sm = *ordered[si - 1];
            const auto& row = alignment[si];
            Accum& acc = accums[sm.oid.getValue()];
            acc.oid = sm.oid;

            // Find the aligned footprint [colFirst, colLast] of this row.
            int64_t colFirst = -1, colLast = -1;
            for (uint64_t col = 0; col < ncols; col++) {
                if (!row[col].isGap()) { if (colFirst < 0) colFirst = int64_t(col); colLast = int64_t(col); }
            }
            if (colFirst < 0) continue;   // empty row (should not happen)

            // Read base position at the start of the footprint: count non-gap
            // read bases before colFirst and add to cBegin.
            uint32_t readAbs = sm.cBegin;
            for (int64_t col = 0; col < colFirst; col++)
                if (!row[col].isGap()) readAbs++;

            int64_t pendingDelBegin = -1;
            bool entered = false;

            // True-pin bounds for this member. A member's readAbs at a backbone
            // column becomes the read's exported position for any het/hom anchor
            // built there, and the downstream monotonicity verifier requires that
            // position to lie strictly between the read's TRUE marker positions at
            // the bracketing backbone anchors. cBegin is a true pin only when the
            // member is pinned at the LEFT boundary (side 0 or 1); cEnd is true
            // only when pinned at the RIGHT boundary (side 0 or 2). A guessed
            // bound (segLenApprox on the free end of a one-sided read) gives no
            // valid coordinate, so columns positioned relative to it must be
            // dropped. Enforcing readAbs in [cBegin, cEnd) against the TRUE pins
            // guarantees every recorded column is monotonic with the backbone
            // anchors it will edge to -- eliminating backward edges from both
            // one-sided guesses and interior POA drift.
            const bool leftPinTrue  = (sm.side == 0 || sm.side == 1);
            const bool rightPinTrue = (sm.side == 0 || sm.side == 2);

            for (int64_t col = colFirst; col <= colLast; col++) {
                const AlignedBase mb = row[col];
                const int64_t bbPos = colBbPos[col];

                if (bbPos >= 0) {
                    // Backbone-bearing column.
                    if (!mb.isGap()) {
                        const uint8_t code = mb.value & 0xff;
                        // Drop columns whose walked readAbs violates a TRUE pin
                        // (POA misalignment or one-sided overshoot): they cannot
                        // be positioned consistently with the backbone anchors.
                        const bool posValid =
                            (!leftPinTrue  || readAbs >= sm.cBegin) &&
                            (!rightPinTrue || readAbs <  sm.cEnd);
                        if (!posValid) { entered = true; readAbs++; continue; }
                        acc.alignedCols.push_back(
                            KwAlignedCol{uint32_t(bbPos), readAbs, code});
                        if (acc.firstBb < 0) acc.firstBb = bbPos;
                        acc.lastBb = bbPos;
                        const uint8_t bbBase = bbSeqVec[uint32_t(bbPos)] & 3;
                        if (code < 4 && code != bbBase)
                            acc.snps.push_back(KwSnp{uint32_t(bbPos), code});  // altBase = code
                        if (pendingDelBegin >= 0) {
                            acc.deletionRanges.push_back(
                                {uint32_t(pendingDelBegin), uint32_t(bbPos)});
                            pendingDelBegin = -1;
                        }
                        entered = true;
                        readAbs++;
                    } else {
                        // Member gap at a backbone position => deletion, only
                        // after the read has entered its footprint (guaranteed
                        // here since we start at colFirst).
                        if (entered && pendingDelBegin < 0)
                            pendingDelBegin = bbPos;
                    }
                } else {
                    // Backbone-gap column (insertion). Advance read pos if the
                    // member has a base; not a SNP column.
                    if (!mb.isGap()) readAbs++;
                }
            }
            // A deletion run open within the footprint is closed at the column
            // after colLast's backbone position (i.e. the read's last aligned
            // backbone pos + 1). Since the footprint ends on a non-gap column,
            // pendingDelBegin is always resolved; nothing to flush past the end.
            (void)pendingDelBegin;
        }
    }

    // --- Flush accumulators into KwMemberProfile rows -----------------------
    vector<KwMemberProfile> profiles;
    profiles.reserve(accums.size());
    for (auto& [oidValue, acc] : accums) {
        (void)oidValue;
        if (acc.firstBb < 0 || acc.lastBb < 0) continue;
        if (acc.alignedCols.empty()) continue;

        KwMemberProfile prof;
        prof.oid = acc.oid;
        prof.bbCovBegin = uint32_t(acc.firstBb);
        prof.bbCovEnd = uint32_t(acc.lastBb) + 1;
        prof.snps = std::move(acc.snps);

        // alignedCols are produced interval-by-interval in increasing bbPos, but
        // two intervals share their boundary anchor position, so a column can be
        // emitted twice (once as an interval's right endpoint, once as the
        // next interval's left endpoint). Sort and dedup by bbPos (keep first)
        // so colAt()'s binary search is valid and support is not double-counted.
        prof.alignedCols = std::move(acc.alignedCols);
        sort(prof.alignedCols.begin(), prof.alignedCols.end(),
            [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos < b.bbPos; });
        prof.alignedCols.erase(
            unique(prof.alignedCols.begin(), prof.alignedCols.end(),
                [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos == b.bbPos; }),
            prof.alignedCols.end());

        // Enforce read-position monotonicity: as bbPos increases, readPos must
        // strictly increase (a read's bases advance along the backbone). The
        // one-sided (kOV) path can place a free-end base at a readPos that does
        // not fit the read's true coordinates (its subsequence was bounded by a
        // fixed base count, so indels drift it), which would later fail the
        // downstream monotonicity verifier. Drop any column that violates the
        // invariant and any SNP at a dropped position. No-op for two-sided-only
        // profiles (already monotone by construction).
        {
            vector<KwAlignedCol> kept;
            kept.reserve(prof.alignedCols.size());
            uint32_t lastReadPos = 0; bool have = false;
            unordered_map<uint32_t, uint8_t> keptBaseAt;   // bbPos -> readBase kept
            for (const KwAlignedCol& c : prof.alignedCols) {
                if (have && c.readPos <= lastReadPos) continue;  // backward/stall: drop
                kept.push_back(c);
                keptBaseAt[c.bbPos] = c.readBase;
                lastReadPos = c.readPos;
                have = true;
            }
            prof.alignedCols.swap(kept);
            // Drop SNPs whose column was removed (or whose base no longer matches
            // the kept column at that position).
            vector<KwSnp> keptSnps;
            keptSnps.reserve(prof.snps.size());
            for (const KwSnp& s : prof.snps) {
                auto it = keptBaseAt.find(s.bbPos);
                if (it != keptBaseAt.end() && it->second == s.altBase)
                    keptSnps.push_back(s);
            }
            prof.snps.swap(keptSnps);
        }

        // SNPs can likewise be duplicated at shared boundaries; dedup by
        // (bbPos, base).
        sort(prof.snps.begin(), prof.snps.end(),
            [](const KwSnp& a, const KwSnp& b) {
                return a.bbPos < b.bbPos || (a.bbPos == b.bbPos && a.altBase < b.altBase);
            });
        prof.snps.erase(
            unique(prof.snps.begin(), prof.snps.end(),
                [](const KwSnp& a, const KwSnp& b) {
                    return a.bbPos == b.bbPos && a.altBase == b.altBase;
                }),
            prof.snps.end());

        // Merge deletion ranges (they may abut across interval boundaries).
        if (!acc.deletionRanges.empty()) {
            sort(acc.deletionRanges.begin(), acc.deletionRanges.end());
            vector<pair<uint32_t, uint32_t>> merged;
            merged.reserve(acc.deletionRanges.size());
            for (const auto& r : acc.deletionRanges) {
                if (!merged.empty() && r.first <= merged.back().second)
                    merged.back().second = max(merged.back().second, r.second);
                else merged.push_back(r);
            }
            prof.deletionRanges = std::move(merged);
        }

        profiles.push_back(std::move(prof));
    }

    if (debug) {
        const double oneSidedFrac = (twoSidedTotal + oneSidedTotal) > 0
            ? double(oneSidedTotal) / double(twoSidedTotal + oneSidedTotal) : 0.0;
        cout << "    intervalPoaDetectHetBubbles bb=" << bbOid
             << " window=[" << windowBbBegin << "," << windowBbEnd << ")"
             << " members=" << members.size()
             << " intervals=" << intervalsRun
             << " twoSided=" << twoSidedTotal
             << " oneSided=" << oneSidedTotal
             << " oneSidedPlaced=" << oneSidedPlaced
             << " oneSidedFrac=" << std::fixed << std::setprecision(3) << oneSidedFrac
             << std::defaultfloat
             << " intervalsWithOneSided=" << intervalsWithOneSided
             << " profiles=" << profiles.size() << endl;
    }

    // Deduplicate profiles by PHYSICAL ReadId. A read can appear in the window
    // under both orientations (forward oid and its reverse complement) as two
    // separate members; the strict both-anchors gate rarely gives both enough
    // overlapping coverage to matter, but one-sided placement can extend both
    // far enough that the same ReadId lands in one anchor on both strands, which
    // the external-anchor export rejects. Keep, per ReadId, the profile with the
    // most aligned columns (strongest evidence) and drop the other.
    {
        unordered_map<uint32_t, uint32_t> bestByReadId;  // readId -> index in profiles
        vector<bool> drop(profiles.size(), false);
        for (uint32_t i = 0; i < profiles.size(); i++) {
            const uint32_t rid = profiles[i].oid.getReadId();
            auto it = bestByReadId.find(rid);
            if (it == bestByReadId.end()) {
                bestByReadId[rid] = i;
            } else {
                const uint32_t j = it->second;
                if (profiles[i].alignedCols.size() > profiles[j].alignedCols.size()) {
                    drop[j] = true; bestByReadId[rid] = i;
                } else {
                    drop[i] = true;
                }
            }
        }
        if (bestByReadId.size() != profiles.size()) {
            vector<KwMemberProfile> kept;
            kept.reserve(bestByReadId.size());
            for (uint32_t i = 0; i < profiles.size(); i++)
                if (!drop[i]) kept.push_back(std::move(profiles[i]));
            profiles.swap(kept);
        }
    }

    if (profiles.size() < 2) return 0;

    // Diploid het coverage estimate for the auto-derived minSupport rule (same
    // source as the ksw2/abPOA paths).
    const uint64_t coverageHet = assemblerInfo.isOpen ?
        assemblerInfo->kmerDistributionInfo.coverageHet : invalid<uint64_t>;

    return emitHetBubblesFromProfiles(
        window, profiles, bbSeqVec, windowBbBegin, windowBbEnd,
        bbOid, k, hetMinVaf, hetMinSupport,
        hetDropHomopolymer, hetDropRepeat, coverageHet, "intervalpoa");
}
