#ifndef DINARA_WINDOW_INTERVAL_POA_HPP
#define DINARA_WINDOW_INTERVAL_POA_HPP

/// @file WindowIntervalPoa.hpp
/// @brief Decomposed per-interval POA het detection, split so intervals from ALL
///        windows can be flattened into one global work list and load-balanced
///        across threads (shasta2 assembleChainsMultithreaded model).
///
/// The single-window entry point intervalPoaDetectHetBubblesInWindow and the
/// all-windows orchestrator intervalPoaDetectHetBubblesAllWindows both build on
/// three pieces defined here:
///
///   buildIpoaPlan(window)   -> IpoaPlan   (all per-window precomputed state:
///                                          backbone frame, member pins, interval
///                                          breakpoints, per-breakpoint member
///                                          index). Cheap; done once per window.
///   runIpoaInterval(plan,bi)-> IpoaFragment (ONE abPOA MSA for interval bi,
///                                          writing per-read column/SNP/deletion
///                                          contributions to its OWN fragment --
///                                          no shared state, so steps run in
///                                          parallel like shasta2 AssemblySteps).
///   mergeAndEmitIpoaWindow(plan, fragments) -> uint32_t (serial per window:
///                                          concatenate fragments into per-read
///                                          KwMemberProfiles -- the concat is
///                                          order-independent because cols/snps
///                                          are sorted+deduped and firstBb/lastBb
///                                          are min/max -- then run the shared
///                                          emitHetBubblesFromProfiles tail).
///
/// This mirrors shasta2::mode3::AssemblyGraph::assembleChainsMultithreaded:
/// each consecutive-anchor step is one unit of work, all steps are flattened
/// across the whole graph, sorted by descending cost, and pulled by threads
/// one at a time (batch=1). Each step writes a disjoint slot (there: a chain's
/// stepSequences[i]; here: a fragment), then a serial pass combines them.

#include "AnchorWindows.hpp"
#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "WindowHetProfiles.hpp"
#include "globalMsa.hpp"
#include "invalid.hpp"

#include <abpoa.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace dinara {

// --- Optional per-phase timing (DINARA_HET_TIMING=1) --------------------------
// Off by default and effectively free when off (one bool check per phase). Each
// worker thread accumulates into thread-local nanosecond counters, then flushes
// them into these global atomics once, so there is no per-interval atomic
// contention. The orchestrator reads them after joining and prints a breakdown.
struct IpoaTiming {
    std::atomic<std::uint64_t> planNs{0};        // buildIpoaPlan
    std::atomic<std::uint64_t> setupNs{0};       // member gather + code arrays
    std::atomic<std::uint64_t> msaNs{0};         // abpoa_reset + abpoa_msa
    std::atomic<std::uint64_t> extractNs{0};     // MSA copy + column/SNP walk
    std::atomic<std::uint64_t> accumNs{0};       // fold fragment into window
    std::atomic<std::uint64_t> mergeNs{0};       // mergeAndEmitIpoaWindow
    std::atomic<std::uint64_t> resetNs{0};       // abpoa_reset only
    std::atomic<std::uint64_t> intervals{0};     // intervals attempted
    std::atomic<std::uint64_t> msaRuns{0};       // intervals that ran an MSA
    std::atomic<std::uint64_t> skippedIntervals{0}; // skipped by member pre-filter
    std::atomic<std::uint64_t> skippedByLen{0};  // skipped: maxLen over DP budget
    std::atomic<std::uint64_t> seqSum{0};        // sum of nSeq over runs
    std::atomic<std::uint64_t> seqMax{0};        // max nSeq
    std::atomic<std::uint64_t> lenSum{0};        // sum of maxLen over runs
    std::atomic<std::uint64_t> lenMax{0};        // max maxLen
    std::atomic<std::uint64_t> colSum{0};        // sum of msa_len over runs
    std::atomic<std::uint64_t> colMax{0};        // max msa_len
    // Per-run (reset+msa) time distribution, to see outlier vs uniform cost.
    std::atomic<std::uint64_t> runMaxNs{0};      // slowest single reset+msa
    std::atomic<std::uint64_t> runOver100us{0};  // runs > 100 us
    std::atomic<std::uint64_t> runOver1ms{0};    // runs > 1 ms
    std::atomic<std::uint64_t> runOver10ms{0};   // runs > 10 ms
    std::atomic<std::uint64_t> nsOver1ms{0};     // total ns spent by >1ms runs
    void reset() {
        planNs = 0; setupNs = 0; msaNs = 0; extractNs = 0;
        accumNs = 0; mergeNs = 0; resetNs = 0;
        intervals = 0; msaRuns = 0; skippedIntervals = 0; skippedByLen = 0;
        seqSum = 0; seqMax = 0; lenSum = 0; lenMax = 0; colSum = 0; colMax = 0;
        runMaxNs = 0; runOver100us = 0; runOver1ms = 0; runOver10ms = 0;
        nsOver1ms = 0;
    }
};

inline void ipoaAtomicMax(std::atomic<std::uint64_t>& a, std::uint64_t v) {
    std::uint64_t cur = a.load(std::memory_order_relaxed);
    while (v > cur && !a.compare_exchange_weak(cur, v,
                                               std::memory_order_relaxed)) {}
}

inline IpoaTiming& ipoaTiming() { static IpoaTiming t; return t; }
inline bool ipoaTimingEnabled() {
    static const bool on = (std::getenv("DINARA_HET_TIMING") != nullptr);
    return on;
}

// Max interval span (backbone bp between consecutive breakpoints) that will be
// POA'd. 0 = NO CAP (default): every region is processed. This exists only as a
// diagnostic / profiling knob, NOT a correctness setting -- skipping wide
// intervals silently drops any het inside them, which is unacceptable for a
// variant caller. Wide anchor-less gaps (repeats / low coverage / SV) are slow
// (a global affine POA is O(span * graphNodes) over ~50 reads) and drift-prone
// (one-sided members inject a full span-length GUESSED sequence), but the fix
// is to SUB-TILE them into small POAs, not to skip them. Set
// DINARA_IPOA_MAX_INTERVAL_BP=N only to measure the cost of that tail.
inline std::uint32_t ipoaMaxIntervalBp() {
    static const std::uint32_t v = []() -> std::uint32_t {
        const char* e = std::getenv("DINARA_IPOA_MAX_INTERVAL_BP");
        if (e == nullptr) return 0u;   // no cap: process all regions
        char* endp = nullptr;
        const unsigned long p = std::strtoul(e, &endp, 10);
        return (endp != e) ? std::uint32_t(p) : 0u;
    }();
    return v;
}
using IpoaClock = std::chrono::steady_clock;
inline std::uint64_t ipoaNsSince(const IpoaClock::time_point& t0) {
    return std::uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
        IpoaClock::now() - t0).count());
}

// A member read pinned at a breakpoint, with its base position there.
struct IpoaMemberBp {
    std::uint32_t oidValue;
    std::uint32_t cPos;
};

// All per-window state needed to run any interval independently and to merge +
// emit afterwards. Built once per window by buildIpoaPlan.
struct IpoaPlan {
    bool valid = false;

    OrientedReadId bbOid;
    std::uint64_t k = 0;

    // Backbone base frame (marker START .. last marker START + k). backbone
    // offsets used everywhere are absolute oriented positions on the bb read.
    std::uint32_t windowBbBegin = 0;
    std::uint32_t windowBbEnd = 0;

    // Backbone read length. Backbone bases are read on demand from the Reads
    // store (runIpoaInterval) or materialized transiently at merge time
    // (buildBbSeqVec), NEVER held per-window across all windows -- that full-read
    // vector, times thousands of windows, is a memory blowup.
    std::uint32_t bbLen = 0;

    // Free the heavy per-window state after this window has been merged+emitted.
    // Keeps resident memory bounded to in-flight windows (like the old
    // per-window path) instead of all windows at once.
    void freeHeavy() {
        std::vector<std::uint32_t>().swap(breakpoints);
        std::vector<std::vector<IpoaMemberBp>>().swap(atBreakpoint);
    }

    // Sorted unique interval breakpoints (backbone positions where some read is
    // pinned). Intervals are [breakpoints[bi], breakpoints[bi+1]).
    std::vector<std::uint32_t> breakpoints;

    // atBreakpoint[i] = members pinned at breakpoints[i], sorted by oidValue.
    std::vector<std::vector<IpoaMemberBp>> atBreakpoint;

    // Member count (for debug / gating only).
    std::size_t memberCount = 0;

    std::uint32_t intervalCount() const {
        return breakpoints.size() >= 2 ?
            std::uint32_t(breakpoints.size() - 1) : 0u;
    }
};

// Per-read contribution from ONE interval. Concatenated (order-independent) into
// the window's accumulators during merge.
struct IpoaReadFrag {
    std::uint32_t oidValue = 0;
    std::vector<KwAlignedCol> alignedCols;
    std::vector<KwSnp> snps;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> deletionRanges;
    // Backbone positions with an insertion immediately BEFORE them (the aligned
    // base at this position is preceded by >=1 read base with no backbone
    // column). Mirrors deletionRanges' use of the trailing backbone position.
    std::vector<std::uint32_t> insertionSites;
    std::int64_t firstBb = -1;
    std::int64_t lastBb = -1;
};

// Output of one interval step.
struct IpoaFragment {
    std::vector<IpoaReadFrag> reads;
    // Instrumentation (summed over intervals during merge).
    std::uint64_t twoSided = 0;
    std::uint64_t oneSided = 0;
    std::uint64_t oneSidedPlaced = 0;
    bool ran = false;          // an MSA was actually performed
    bool hadOneSided = false;  // this interval had >=1 one-sided read
};

// Build the per-window plan. mkrs = *markers. Returns a plan with valid=false if
// the window cannot yield het bubbles (too small, missing ordinals, <2 members,
// <2 breakpoints); callers treat that as "0 bubbles".
inline IpoaPlan buildIpoaPlan(
    const AnchorWindow& window,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const MemoryMapped::VectorOfVectors<CompressedMarker, std::uint64_t>& mkrs,
    const Reads& rds,
    std::uint64_t k)
{
    using std::vector;
    using std::uint32_t;
    using std::uint64_t;
    using std::unordered_map;
    using std::sort;

    IpoaPlan plan;
    plan.bbOid = window.backboneOrientedReadId;
    plan.k = k;

    const OrientedReadId bbOid = plan.bbOid;
    const ReadId bbReadId = bbOid.getReadId();
    const uint32_t bbLen = uint32_t(rds.getRead(bbReadId).baseCount);
    const auto bbJ = journeys[bbOid];

    if (window.backboneEnd <= window.backboneBegin + 1) return plan;

    const uint32_t firstAnchorJP = window.backboneBegin;
    const uint32_t lastAnchorJP = window.backboneEnd - 1;
    if (firstAnchorJP >= bbJ.size() || lastAnchorJP >= bbJ.size()) return plan;

    const uint32_t firstOrd = anchors.getOrdinal(bbJ[firstAnchorJP], bbOid);
    const uint32_t lastOrd = anchors.getOrdinal(bbJ[lastAnchorJP], bbOid);
    if (firstOrd == invalid<uint32_t> || lastOrd == invalid<uint32_t>) return plan;

    plan.windowBbBegin = mkrs[bbOid.getValue()][firstOrd].position;
    plan.windowBbEnd = mkrs[bbOid.getValue()][lastOrd].position + uint32_t(k);
    if (plan.windowBbEnd <= plan.windowBbBegin) return plan;

    plan.bbLen = bbLen;

    // Live-anchor pin enumeration (identical to the ksw2 path). Each backbone
    // anchor a read shares in the window becomes a pin (bbPos, cPos).
    struct Pin { uint32_t bbPos; uint32_t cPos; };
    unordered_map<uint64_t, vector<Pin>> readPins;
    for (uint32_t jp = window.backboneBegin; jp < window.backboneEnd; jp++) {
        if (jp >= bbJ.size()) break;
        const Shasta2AnchorId aid = bbJ[jp];
        const uint32_t bbOrd = anchors.getOrdinal(aid, bbOid);
        if (bbOrd == invalid<uint32_t> || bbOrd >= mkrs[bbOid.getValue()].size())
            continue;
        const uint32_t bbPos = mkrs[bbOid.getValue()][bbOrd].position;
        if (bbPos < plan.windowBbBegin || bbPos >= plan.windowBbEnd) continue;

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

    // LIS over read positions to drop off-diagonal (repeat) pins.
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
    if (members.size() < 2) return plan;
    plan.memberCount = members.size();

    // Breakpoints = union of all members' pin positions.
    for (const Member& m : members)
        for (const Pin& p : m.pins)
            plan.breakpoints.push_back(p.bbPos);
    sort(plan.breakpoints.begin(), plan.breakpoints.end());
    plan.breakpoints.erase(
        unique(plan.breakpoints.begin(), plan.breakpoints.end()),
        plan.breakpoints.end());
    if (plan.breakpoints.size() < 2) return plan;

    unordered_map<uint32_t, uint32_t> bpIndex;
    bpIndex.reserve(plan.breakpoints.size() * 2);
    for (uint32_t i = 0; i < plan.breakpoints.size(); i++)
        bpIndex[plan.breakpoints[i]] = i;

    plan.atBreakpoint.assign(plan.breakpoints.size(), {});
    for (const Member& m : members) {
        for (const Pin& p : m.pins) {
            auto it = bpIndex.find(p.bbPos);
            if (it != bpIndex.end())
                plan.atBreakpoint[it->second].push_back({m.oid.getValue(), p.cPos});
        }
    }
    for (auto& v : plan.atBreakpoint)
        sort(v.begin(), v.end(),
            [](const IpoaMemberBp& a, const IpoaMemberBp& b) {
                return a.oidValue < b.oidValue;
            });

    plan.valid = true;
    return plan;
}

// Whether one-sided (semi-global) reads are placed. ON by default; opt out with
// DINARA_IPOA_ONESIDED=0. Read once by the orchestrator and passed in so every
// interval step sees the same value without repeated getenv.
inline bool ipoaOneSidedEnabled() {
    const char* e = std::getenv("DINARA_IPOA_ONESIDED");
    return (e == nullptr) || (e[0] != '0');
}

// abPOA handle wrapper. Thread-local; reset per interval by abpoa_reset. Long-
// read-tuned affine scoring, adaptive banding, minimizer seeding + partitioning
// (so wide gaps become many tiny window POAs, not one huge DP), row order
// preserved so row 0 is always the backbone.
struct IpoaAbHandle {
    abpoa_t* ab = nullptr;
    abpoa_para_t* abpt = nullptr;
    IpoaAbHandle() {
        ab = abpoa_init();
        abpt = abpoa_init_para();
        abpt->align_mode = ABPOA_GLOBAL_MODE;
        abpt->gap_mode = ABPOA_AFFINE_GAP;
        abpt->match = 2;
        abpt->mismatch = 4;
        abpt->gap_open1 = 4;
        abpt->gap_ext1 = 2;
        abpt->gap_open2 = 0;
        abpt->gap_ext2 = 0;
        abpt->wb = 10;
        abpt->wf = 0.01;
        // Minimizer-based seeding + partitioning (shasta2 LocalAssembly values).
        // With disable_seeding=1 AND progressive_poa=0, abpoa takes the
        // abpoa_poa() branch: ONE monolithic O(qlen * graphNodes) DP over the
        // whole interval. On wide anchor-less gaps (up to ~4.7 kb x ~50 reads)
        // that single DP was ~500 ms and dominated het time. Leaving seeding on
        // (as shasta2 does) makes abpoa take the abpoa_anchor_poa() branch: it
        // finds shared minimizers, partitions each sequence into small windows
        // at those anchors, and runs a tiny DP per window -- exactly the
        // sub-tiling we'd otherwise hand-roll, done inside abpoa. Row order is
        // still input order (abpoa_anchor_poa adds each read at its ORIGINAL
        // index read_id = read_id_map[_i], so msa_base[0] stays the backbone).
        abpt->disable_seeding = 0;
        abpt->w = 6;
        abpt->k = 9;
        abpt->min_w = 10;
        abpt->progressive_poa = 0;
        abpt->sort_input_seq = 0;   // keep row 0 = backbone
        abpt->out_msa = 1;
        abpt->out_cons = 0;
        abpt->ret_cigar = 1;
        abpoa_post_set_para(abpt);
    }
    ~IpoaAbHandle() { if (ab) abpoa_free(ab); if (abpt) abpoa_free_para(abpt); }
};

// A member's placement over a backbone range [bbBegin, bbEnd): the read's base
// interval [cBegin, cEnd) and which ends are TRUE pins. side: 0=both, 1=left-
// only, 2=right-only. Interval mode may fill the free end of a one-sided member
// with a guess; range mode always uses true pins (both ends), so side==0 there.
struct IpoaSpanMember {
    OrientedReadId oid;
    std::uint32_t cBegin;
    std::uint32_t cEnd;
    std::uint8_t side;
};

inline void runIpoaOnRows(
    const IpoaPlan& plan, std::uint32_t bbBegin, std::uint32_t bbEnd,
    std::vector<IpoaSpanMember>& spanMembers, const Reads& rds,
    IpoaAbHandle& ah, IpoaFragment& frag, bool timing,
    const IpoaClock::time_point& tSetup0, int minMembers);

// ============================================================================
// Leaf-snarl detection directly from the abpoa graph (substitution-only).
// ============================================================================
//
// emitHetBubblesFromProfiles (WindowHetProfiles.hpp) scans the flattened MSA
// column matrix (abc->msa_base[r][c]) one backbone column at a time, so a
// genuine multi-base block shared identically by the same reads is never
// recognized as one bubble -- it is only ever found (if at all) as several
// unrelated single-column sites. This works directly on abpoa's own graph
// instead, which additionally lets it tell apart "same value by coincidence"
// from "genuinely the same aligned position" when nearby structural
// complexity is involved.
//
// Per-read node membership is captured DIRECTLY from abpoa's own alignment
// result (qpos_to_node_id, filled in by abpoa_add_subgraph_alignment /
// abpoa_add_graph_alignment during construction) rather than reconstructed
// afterward by walking read_ids bitmasks on graph edges. That reconstruction
// was tried first and has a real bug for multi-segment members (a read added
// piece by piece, one abpoa_add_subgraph_alignment call per inter-anchor
// segment, matching AssemblerWindowAbpoaGraph.cpp): abpoa registers a read's
// bit on the OUTGOING edge of a shared anchor node only when a further
// segment continues past it (inc_both_ends=1 on the next segment's own first
// step); a member whose own rightmost pin IS that anchor -- i.e. exactly its
// own endpoint, which differs read by read -- has no such edge, even though
// it genuinely matches there, so edge-based reconstruction misclassifies
// every read as "diverged" exactly at its own endpoint. Reading
// qpos_to_node_id directly avoids this: it is abpoa's own record of exactly
// which node each of a read's bases aligned to, with no dependency on which
// edge happened to get a bit registered.
//
// Grouping rule: consecutive backbone columns belong to the SAME leaf snarl
// only while the exact SET of diverged (off-backbone) members is identical
// from one column to the next -- not merely "some divergence exists". Two
// reads that diverge independently at adjacent columns have DIFFERENT
// diverged-row sets at those two columns (e.g. {R2} vs {R3}), so they are
// kept as two separate single-column snarls; a genuine multi-base block
// shared identically by the same reads keeps the same set across all its
// columns and is correctly merged into one multi-base snarl. A column where
// the diverged set changes at all -- growing, shrinking, or swapping
// membership -- starts a new snarl there, ending the previous one at that
// boundary.
//
// Scoped to substitutions only: an anchor needs a real base position for
// every member, which a deletion or insertion doesn't have. A member with an
// Invalid (non-substitution) divergence anywhere in a candidate span is
// excluded from allele counting there -- same treatment as "not covered" --
// rather than disqualifying the whole site for every other member. (An
// earlier version disqualified the whole site, which in real noisy data --
// almost every multi-row site has at least one indel-affected read mixed in
// with otherwise-clean substitution support -- silently discarded
// everything.)

struct LeafSnarlAllele {
    std::vector<OrientedReadId> members;
    std::vector<std::uint8_t> bases;  // length == end - start - 1; one real base per column.
};

struct LeafSnarl {
    std::uint32_t start;  // absolute backbone position: last full-agreement column before.
    std::uint32_t end;    // absolute backbone position: first full-agreement column after.
    std::vector<LeafSnarlAllele> alleles;  // >= 2; includes the backbone-matching (ref) allele.
};

// NotCovered is distinct from Invalid: Invalid means the member covers this
// backbone position but has no substitution there (a deletion, or the
// alignment simply doesn't reach it -- ambiguous either way, so scoped out).
// NotCovered means this backbone position is outside the member's own
// coverage range entirely -- it was never going to align here, which is the
// common case for any window wider than a single read's overlap, and must
// NOT be treated as divergence (see the design note on IpoaMemberInfo).
enum class IpoaColState : std::uint8_t { OnBackbone, Substitution, Invalid, NotCovered };

// Small per-member metadata: identity plus its coverage bound.
// bbCovBegin/bbCovEnd (backbone position indices, half-open, same indexing as
// backbonePath/bbBases below) bound where this member actually has coverage
// -- normally its first and last shared anchor with the backbone. A window
// is typically far wider than any single member's own overlap with the
// backbone (a window chains together many reads, each covering only part of
// it), so most members do not cover most positions. Without this bound,
// "doesn't visit this node" is indistinguishable from "was never aligned
// here at all", and every uncovered position gets misclassified as a
// deletion for every member that simply isn't there -- in practice this
// showed up as ~100% of positions in a window appearing to "diverge".
struct IpoaMemberInfo {
    OrientedReadId oid;
    std::uint32_t bbCovBegin = 0;
    std::uint32_t bbCovEnd = 0;
};

// Which members visited each graph node, captured directly from abpoa's own
// alignment result (qpos_to_node_id, filled in by abpoa_add_subgraph_alignment
// / abpoa_add_graph_alignment during construction) rather than reconstructed
// afterward by walking read_ids bitmasks on graph edges. That reconstruction
// was tried first and has a real bug for multi-segment members (a read added
// piece by piece, one abpoa_add_subgraph_alignment call per inter-anchor
// segment, matching AssemblerWindowAbpoaGraph.cpp): abpoa registers a read's
// bit on the OUTGOING edge of a shared anchor node only when a further
// segment continues past it (inc_both_ends=1 on the next segment's own first
// step); a member whose own rightmost pin IS that anchor -- i.e. exactly its
// own endpoint, which differs read by read -- has no such edge, even though
// it genuinely matches there, so edge-based reconstruction misclassifies
// every read as "diverged" exactly at its own endpoint. Reading
// qpos_to_node_id directly avoids this: it is abpoa's own record of exactly
// which node each of a read's bases aligned to, with no dependency on which
// edge happened to get a bit registered.
//
// Stored as a CSR (offset + flat member-index array), built once from the
// flat list of (nodeId, memberIndex) pairs collected during construction --
// NOT as one vector<bool> per member sized to the whole graph. A member
// typically visits only a small fraction of a window's nodes, so a per-member
// bitset wastes memory proportional to graph size rather than actual
// coverage, and turns classification into "ask every member a question"
// instead of "look up who's already here".
struct IpoaNodeVisitors {
    std::vector<std::uint32_t> offset;  // size nodeCount + 1
    std::vector<int> memberIdx;         // flat, size offset.back()

    // Member indices that visited nodeId, as a raw (pointer, count) span.
    std::pair<const int*, std::uint32_t> at(int nodeId) const {
        if (nodeId < 0 || std::size_t(nodeId) + 1 >= offset.size()) return {nullptr, 0};
        const std::uint32_t begin = offset[std::size_t(nodeId)];
        const std::uint32_t end = offset[std::size_t(nodeId) + 1];
        return {memberIdx.data() + begin, end - begin};
    }

    // Build the CSR from a flat (possibly unsorted) list of (nodeId, member)
    // pairs via counting sort -- O(visits + nodeCount), no comparison sort.
    static IpoaNodeVisitors build(
        int nodeCount, const std::vector<std::pair<int, int>>& visits)
    {
        IpoaNodeVisitors nv;
        nv.offset.assign(std::size_t(nodeCount) + 1, 0);
        for (const auto& [nodeId, memberIdx] : visits) {
            (void)memberIdx;
            if (nodeId >= 0 && nodeId < nodeCount) nv.offset[std::size_t(nodeId) + 1]++;
        }
        for (std::size_t k = 1; k < nv.offset.size(); k++) nv.offset[k] += nv.offset[k - 1];
        nv.memberIdx.resize(nv.offset.back());
        std::vector<std::uint32_t> cursor(nv.offset.begin(), nv.offset.end() - 1);
        for (const auto& [nodeId, memberIdx] : visits) {
            if (nodeId >= 0 && nodeId < nodeCount) {
                nv.memberIdx[cursor[std::size_t(nodeId)]++] = memberIdx;
            }
        }
        return nv;
    }
};

// backbonePath[i] = the graph node id of the backbone's own base at position
// i, captured directly from qpos_to_node_id at seed time (not reconstructed
// via edge traversal, for the same reason as IpoaNodeVisitors above).
// bbBases/bbBeginAbs = the backbone's own base codes over the same span and
// the absolute backbone position of bbBases[0]. minSupport/minVaf/
// dropHomopolymer/dropRepeat mirror the identical-named gates in
// emitHetBubblesFromProfiles (WindowHetProfiles.hpp).
//
// Performance note: state/altBase are one flat, contiguous, row-major array
// each (index = i*nMembers+r), not vector<vector<T>> -- the latter is n
// separate heap allocations, and walks between them on every position-major
// access (which is the access pattern used throughout this function: outer
// loop over position, inner loop over member). diverged is similarly one
// flat CSR (offset + flat member-index array) instead of vector<vector<int>>,
// for the same reason.
// Shared classification+grouping core, independent of how state/altBase were
// populated (from an abPOA graph, from independent pairwise alignments, or
// any other source): given the fully-built per-position/per-member state,
// group runs of identically-diverged columns into leaf snarls. See
// findLeafSnarlsFromGraph and findLeafSnarlsFromPairwiseColumns for the two
// current front-ends that build state/altBase and delegate here.
inline std::vector<LeafSnarl> classifyLeafSnarls(
    std::uint32_t n,
    int nMembers,
    const std::vector<IpoaColState>& state,
    const std::vector<std::uint8_t>& altBase,
    const std::vector<IpoaMemberInfo>& members,
    OrientedReadId backboneOid,
    const std::vector<std::uint8_t>& bbBases,
    std::uint32_t bbBeginAbs,
    std::uint64_t minSupport,
    double minVaf,
    bool dropHomopolymer,
    bool dropRepeat)
{
    using std::uint8_t;
    using std::uint32_t;
    using std::uint64_t;
    using std::vector;

    vector<LeafSnarl> result;
    if (n == 0 || n != bbBases.size()) return result;

    const size_t nMembersSz = size_t(nMembers);
    auto idx = [nMembersSz](uint32_t i, int r) { return size_t(i) * nMembersSz + size_t(r); };

    // Per-column set of diverged members, as one flat CSR (offset + flat
    // member-index array) built in two passes over the now-known state array
    // -- counting then filling -- rather than n separate vector<int> rows.
    // Exact-equality comparison between consecutive columns' spans drives
    // grouping (see the file-level design note above for why exact-set
    // matching, not "any divergence", is required). NotCovered is excluded
    // here alongside OnBackbone: a member never aligned at this position
    // carries no information about it and must not count as divergence.
    auto isDiverged = [](IpoaColState s) {
        return s == IpoaColState::Substitution || s == IpoaColState::Invalid;
    };
    vector<uint32_t> divergedOffset(size_t(n) + 1, 0);
    for (uint32_t i = 0; i < n; i++) {
        uint32_t cnt = 0;
        for (int r = 0; r < nMembers; r++) if (isDiverged(state[idx(i, r)])) cnt++;
        divergedOffset[i + 1] = divergedOffset[i] + cnt;
    }
    vector<int> divergedFlat(divergedOffset[n]);
    {
        vector<uint32_t> cursor(divergedOffset.begin(), divergedOffset.end() - 1);
        for (uint32_t i = 0; i < n; i++) {
            for (int r = 0; r < nMembers; r++) {
                if (isDiverged(state[idx(i, r)])) divergedFlat[cursor[i]++] = r;
            }
        }
    }
    auto divergedSpan = [&](uint32_t i) {
        return std::pair<const int*, uint32_t>(
            divergedFlat.data() + divergedOffset[i], divergedOffset[i + 1] - divergedOffset[i]);
    };
    auto divergedEmpty = [&](uint32_t i) { return divergedOffset[i] == divergedOffset[i + 1]; };
    auto divergedEqual = [&](uint32_t a, uint32_t b) {
        const auto [pa, na] = divergedSpan(a);
        const auto [pb, nb] = divergedSpan(b);
        return na == nb && std::equal(pa, pa + na, pb);
    };

    // Flank-linearity gate, generalized from emitHetBubblesFromProfiles's
    // single-column flanksLinear (WindowHetProfiles.hpp) to a [runStart,
    // runEnd) span: requires >=4 clean columns on each side of the WHOLE
    // span (not just around one position), where "clean" means no other
    // divergence -- substitution or Invalid (deletion/insertion) -- reaches
    // minSupport there. This is what keeps a real, isolated site from being
    // called right next to a repeat/indel-prone region, same intent as the
    // single-column version, just checked against the span's own two edges
    // instead of a single position's immediate neighbors.
    auto flanksLinear = [&](uint32_t spanStart, uint32_t spanEnd) -> bool {
        if (spanStart < 4 || spanEnd + 4 > n) return false;
        for (const uint32_t c : {spanStart - 4, spanStart - 3, spanStart - 2, spanStart - 1,
                                  spanEnd,       spanEnd + 1,   spanEnd + 2,   spanEnd + 3}) {
            uint64_t invalidCount = 0;
            uint64_t altCounts[4] = {0, 0, 0, 0};
            for (int r = 0; r < nMembers; r++) {
                const size_t cell = idx(c, r);
                if (state[cell] == IpoaColState::Invalid) {
                    invalidCount++;
                } else if (state[cell] == IpoaColState::Substitution) {
                    const uint8_t b = altBase[cell];
                    if (b < 4) altCounts[b]++;
                }
            }
            if (invalidCount >= minSupport) return false;
            for (uint64_t cnt : altCounts) if (cnt >= minSupport) return false;
        }
        return true;
    };

    uint32_t i = 0;
    while (i < n) {
        if (divergedEmpty(i)) { i++; continue; }
        const uint32_t runStart = i;
        uint32_t runEnd = i + 1;
        while (runEnd < n && divergedEqual(runEnd, runStart)) runEnd++;

        if (!flanksLinear(runStart, runEnd)) { i = runEnd; continue; }

        {
            // A member with an Invalid or NotCovered state anywhere in the
            // span is excluded from allele counting for this span (see
            // design notes on IpoaColState/IpoaMemberInfo): Invalid means
            // real but non-substitution divergence there; NotCovered means
            // this member has no information about part of the span at all,
            // so it cannot contribute a full-span key either way.
            vector<bool> rowExcluded(nMembersSz, false);
            for (int r = 0; r < nMembers; r++) {
                for (uint32_t c = runStart; c < runEnd; c++) {
                    const IpoaColState s = state[idx(c, r)];
                    if (s == IpoaColState::Invalid || s == IpoaColState::NotCovered) {
                        rowExcluded[size_t(r)] = true;
                        break;
                    }
                }
            }

            auto keyFor = [&](int r) {
                vector<uint8_t> k;
                k.reserve(runEnd - runStart);
                for (uint32_t c = runStart; c < runEnd; c++) {
                    const size_t cell = idx(c, r);
                    k.push_back(state[cell] == IpoaColState::OnBackbone ? bbBases[c] : altBase[cell]);
                }
                return k;
            };

            // Seed with the ref allele (the backbone's own sequence over this
            // span) so it always exists even if every member is excluded.
            vector<vector<uint8_t>> alleleKeys;
            vector<LeafSnarlAllele> alleles;
            {
                const vector<uint8_t> refKey(bbBases.begin() + runStart, bbBases.begin() + runEnd);
                alleleKeys.push_back(refKey);
                LeafSnarlAllele refAllele;
                refAllele.members.push_back(backboneOid);
                refAllele.bases = refKey;
                alleles.push_back(std::move(refAllele));
            }
            const int refAlleleIndex = 0;

            for (int r = 0; r < nMembers; r++) {
                if (rowExcluded[size_t(r)]) continue;
                vector<uint8_t> k = keyFor(r);
                size_t ai = alleleKeys.size();
                for (size_t j = 0; j < alleleKeys.size(); j++) {
                    if (alleleKeys[j] == k) { ai = j; break; }
                }
                if (ai == alleleKeys.size()) {
                    alleleKeys.push_back(k);
                    LeafSnarlAllele allele;
                    allele.bases = k;
                    alleles.push_back(std::move(allele));
                }
                alleles[ai].members.push_back(members[size_t(r)].oid);
            }

            uint64_t totalMembers = 0;
            for (const auto& a : alleles) totalMembers += a.members.size();

            vector<LeafSnarlAllele> gated;
            for (int ai = 0; ai < int(alleles.size()); ai++) {
                LeafSnarlAllele& allele = alleles[size_t(ai)];
                if (allele.members.size() < minSupport) continue;
                if (ai != refAlleleIndex) {
                    const double af = double(allele.members.size()) / double(totalMembers);
                    if (af < minVaf) continue;
                }
                gated.push_back(std::move(allele));
            }

            bool passesRepeatGate = true;
            if (dropHomopolymer || dropRepeat) {
                KmVarKey vkey;
                vkey.type = KmVarType::Snp;
                vkey.refLen = 1; vkey.altLen = 1; vkey.altBase = 0;
                // Check both ends of the span (matches the single-column gate,
                // which checks the SNP's own single position via the same
                // pos-1/pos+1 flank logic inside kmIsRepeatUnitRange).
                for (uint32_t localPos : {runStart, runEnd - 1}) {
                    vkey.pos = localPos;
                    const bool inHomopolymer = kmIsRepeatUnitRange(
                        bbBases.data(), uint32_t(bbBases.size()), vkey, 0, 1, 1);
                    const bool inStr = kmIsRepeatUnitRange(
                        bbBases.data(), uint32_t(bbBases.size()), vkey, 0, 2, 6);
                    if (dropHomopolymer && inHomopolymer) passesRepeatGate = false;
                    if (dropRepeat && inStr) passesRepeatGate = false;
                }
            }

            if (passesRepeatGate && gated.size() >= 2) {
                LeafSnarl snarl;
                snarl.start = bbBeginAbs + runStart - 1;
                snarl.end = bbBeginAbs + runEnd;
                snarl.alleles = std::move(gated);
                result.push_back(std::move(snarl));
            }
        }

        i = runEnd;
    }

    return result;
}

// Populates state/altBase from an abPOA graph's own alignment record
// (qpos_to_node_id / aligned_node_id), then delegates the classification and
// grouping to classifyLeafSnarls. See the design note on IpoaNodeVisitors for
// why membership is read from qpos_to_node_id rather than reconstructed from
// read_ids edge bits.
inline std::vector<LeafSnarl> findLeafSnarlsFromGraph(
    const abpoa_graph_t* abg,
    OrientedReadId backboneOid,
    const std::vector<int>& backbonePath,
    const std::vector<IpoaMemberInfo>& members,
    const IpoaNodeVisitors& nodeVisitors,
    const std::vector<std::uint8_t>& bbBases,
    std::uint32_t bbBeginAbs,
    std::uint64_t minSupport,
    double minVaf,
    bool dropHomopolymer,
    bool dropRepeat)
{
    using std::uint8_t;
    using std::uint32_t;
    using std::vector;

    if (abg == nullptr) return {};

    const uint32_t n = uint32_t(backbonePath.size());
    if (n == 0 || n != bbBases.size()) return {};

    const int nMembers = int(members.size());
    const size_t nMembersSz = size_t(nMembers);

    // Flat, contiguous, row-major (position-major) storage: cell (i, r) is
    // state[size_t(i)*nMembersSz + r]. One allocation each, matching the
    // position-major access pattern used everywhere below.
    auto idx = [nMembersSz](uint32_t i, int r) { return size_t(i) * nMembersSz + size_t(r); };
    vector<IpoaColState> state(size_t(n) * nMembersSz, IpoaColState::NotCovered);
    vector<uint8_t> altBase(size_t(n) * nMembersSz, 0);

    for (uint32_t i = 0; i < n; i++) {
        const int bbNodeId = backbonePath[i];
        const abpoa_node_t& node = abg->node[bbNodeId];

        // Default per covered member: Invalid (pessimistic) unless upgraded
        // below by direct membership in a visitor list. NotCovered members
        // keep the array's global default and are never touched again.
        for (int r = 0; r < nMembers; r++) {
            if (i >= members[size_t(r)].bbCovBegin && i < members[size_t(r)].bbCovEnd) {
                state[idx(i, r)] = IpoaColState::Invalid;
            }
        }

        // Members actually at the backbone's own node: OnBackbone. Iterates
        // only the (typically small) list of members truly there, not every
        // member in the window.
        {
            const auto [ptr, cnt] = nodeVisitors.at(bbNodeId);
            for (uint32_t k = 0; k < cnt; k++) {
                const size_t cell = idx(i, ptr[k]);
                if (state[cell] == IpoaColState::Invalid) state[cell] = IpoaColState::OnBackbone;
            }
        }

        // Members at an aligned alternative (same column, different base):
        // Substitution. aligned_node_n is small (a handful of alternate
        // bases at most), so this stays cheap.
        for (int a = 0; a < node.aligned_node_n; a++) {
            const int altId = node.aligned_node_id[a];
            if (altId < 0) continue;
            const auto [ptr, cnt] = nodeVisitors.at(altId);
            const uint8_t altBaseVal = abg->node[altId].base;
            for (uint32_t k = 0; k < cnt; k++) {
                const size_t cell = idx(i, ptr[k]);
                if (state[cell] == IpoaColState::Invalid) {
                    state[cell] = IpoaColState::Substitution;
                    altBase[cell] = altBaseVal;
                }
            }
        }
    }

    return classifyLeafSnarls(
        n, nMembers, state, altBase, members, backboneOid, bbBases, bbBeginAbs,
        minSupport, minVaf, dropHomopolymer, dropRepeat);
}

// Populates state/altBase from KwMemberProfile (WindowHetProfiles.hpp) --
// the SAME engine-agnostic per-member representation emitHetBubblesFromProfiles
// already consumes -- then delegates to classifyLeafSnarls. This means this
// function plugs directly into whichever engine already built `profiles` for
// a window (the production per-interval abPOA engine, or the ksw2 pairwise
// engine) with ZERO new alignment work: both engines populate profiles by
// aligning each member independently (one small bounded POA/DP per interval
// or per member, no shared multi-sequence graph spanning the whole window),
// so per-member cost here is bounded purely by that member's own overlap
// span and is completely unaffected by how many other members were already
// processed for this window -- see the design notes in
// AssemblerWindowAbpoaGraph.cpp on the whole-window shared-graph blowup this
// avoids by construction.
//
// A profile's aligned column whose readBase differs from the backbone's own
// base there is a substitution, matching is OnBackbone. A position in
// [bbCovBegin, bbCovEnd) with no corresponding column (a deletion), or one
// flagged noisy (untrustworthy locally-dense mismatch/indel region -- same
// exclusion emitHetBubblesFromProfiles's memberCall applies), is Invalid;
// positions outside that range are NotCovered.
inline std::vector<LeafSnarl> findLeafSnarlsFromProfiles(
    OrientedReadId backboneOid,
    const std::vector<KwMemberProfile>& profiles,
    const std::vector<std::uint8_t>& bbBases,
    std::uint32_t bbBeginAbs,
    std::uint64_t minSupport,
    double minVaf,
    bool dropHomopolymer,
    bool dropRepeat)
{
    using std::uint8_t;
    using std::uint32_t;
    using std::vector;

    const uint32_t n = uint32_t(bbBases.size());
    if (n == 0) return {};

    const int nMembers = int(profiles.size());
    const size_t nMembersSz = size_t(nMembers);

    vector<IpoaMemberInfo> members(nMembersSz);
    for (int r = 0; r < nMembers; r++) {
        members[size_t(r)].oid = profiles[size_t(r)].oid;
        members[size_t(r)].bbCovBegin = profiles[size_t(r)].bbCovBegin;
        members[size_t(r)].bbCovEnd = profiles[size_t(r)].bbCovEnd;
    }

    auto idx = [nMembersSz](uint32_t i, int r) { return size_t(i) * nMembersSz + size_t(r); };
    vector<IpoaColState> state(size_t(n) * nMembersSz, IpoaColState::NotCovered);
    vector<uint8_t> altBase(size_t(n) * nMembersSz, 0);

    for (int r = 0; r < nMembers; r++) {
        const KwMemberProfile& prof = profiles[size_t(r)];
        const uint32_t covBegin = std::max(prof.bbCovBegin, bbBeginAbs);
        const uint32_t covEnd = std::min(prof.bbCovEnd, bbBeginAbs + n);
        for (uint32_t pos = covBegin; pos < covEnd; pos++) {
            state[idx(pos - bbBeginAbs, r)] = IpoaColState::Invalid;
        }
        for (const KwAlignedCol& col : prof.alignedCols) {
            if (col.bbPos < bbBeginAbs || col.bbPos >= bbBeginAbs + n) continue;
            if (prof.isNoisy(col.bbPos)) continue; // untrustworthy: leave as Invalid.
            const uint32_t i = col.bbPos - bbBeginAbs;
            const size_t cell = idx(i, r);
            if (col.readBase >= 4) continue; // N: leave as Invalid, not a substitution.
            if (col.readBase == bbBases[i]) {
                state[cell] = IpoaColState::OnBackbone;
            } else {
                state[cell] = IpoaColState::Substitution;
                altBase[cell] = col.readBase;
            }
        }
    }

    return classifyLeafSnarls(
        n, nMembers, state, altBase, members, backboneOid, bbBases, bbBeginAbs,
        minSupport, minVaf, dropHomopolymer, dropRepeat);
}

// Sliding-window CIGAR-density noise tracker: shared by any per-member
// pairwise alignment detector (ksw2, ProjectedAlignment, ...) that needs to
// flag locally dense mismatch/indel clusters as untrustworthy -- a port of
// pgphase's XidQueue, originally duplicated per-file in this codebase's
// per-engine detectors to avoid unrelated naming collisions, consolidated
// here once a second detector needed the identical logic. Each variant event
// (mismatch=1, insertion=len, deletion=len) is pushed in backbone-coordinate
// order; within a window of `win` backbone bases, if the summed event size
// exceeds `maxS`, the spanned backbone interval is flagged noisy. Contiguous/
// overlapping noisy spans are merged. Emitted ranges are half-open
// [begin, end) backbone positions. Without this gate, independent per-member
// alignment noise (dense local mismatch/indel clusters, each member solving
// its own optimal but locally ambiguous gap placement) produces spurious
// apparent divergence.
struct IpoaNoiseTracker {
    int win;
    int maxS;
    std::vector<std::uint32_t> pos;
    std::vector<std::uint32_t> len;
    std::vector<int> count;
    std::size_t front = 0;
    long total = 0;
    long curStart = -1;
    long curEnd = -1;
    std::vector<std::pair<std::uint32_t, std::uint32_t>>& out;

    IpoaNoiseTracker(int win_, int maxS_, std::vector<std::pair<std::uint32_t, std::uint32_t>>& out_)
        : win(win_), maxS(maxS_), out(out_) {}

    void observe(std::uint32_t p, std::uint32_t l, int c) {
        pos.push_back(p);
        len.push_back(l);
        count.push_back(c);
        total += c;
        const std::size_t rear = pos.size() - 1;

        while (front <= rear &&
               std::int64_t(pos[front]) + std::int64_t(len[front]) - 1 <= std::int64_t(p) - win) {
            total -= count[front];
            ++front;
        }

        if (c <= 0) return;
        if (total <= maxS) return;

        const long noisyStart = long(pos[front]);
        const long noisyEnd = long(pos[rear]) + long(len[rear]);

        if (curStart == -1) {
            curStart = noisyStart;
            curEnd = noisyEnd;
            return;
        }
        if (noisyStart <= curEnd) {
            curEnd = std::max(curEnd, noisyEnd);
            return;
        }
        out.push_back({std::uint32_t(curStart), std::uint32_t(curEnd)});
        curStart = noisyStart;
        curEnd = noisyEnd;
    }

    void finish() {
        if (curStart == -1) return;
        out.push_back({std::uint32_t(curStart), std::uint32_t(curEnd)});
        curStart = -1;
        curEnd = -1;
    }
};

// One sparse mismatch, with the member's exact read position at that
// backbone column (not just the alt base) -- needed to PIN a het/hom anchor
// (see sparseColAt/sparsePinnedKmerCol/sparsePinnedPointCol below), not just
// to classify divergence.
struct SparseSnp {
    std::uint32_t bbPos;
    std::uint8_t altBase;
    std::uint32_t readPos;
};

// A member's SPARSE evidence against the backbone (e.g. from
// ProjectedAlignment::sparseMismatches/sparseIndels): only the positions
// where it DIFFERS from the backbone, not a dense per-column list. Distinct
// from KwMemberProfile/findLeafSnarlsFromProfiles above, which requires a
// real KwAlignedCol entry at EVERY covered position to tell "matches
// backbone" apart from "not aligned here" -- appropriate for a dense
// aligner's per-column walk, but wasteful for a sparse-diff aligner where
// the vast majority of covered positions simply match and were never
// explicitly recorded as such.
struct SparseMemberEvidence {
    OrientedReadId oid;
    std::uint32_t bbCovBegin = 0;
    std::uint32_t bbCovEnd = 0;
    std::vector<SparseSnp> snps;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> deletionRanges; // [begin,end)
    // Piecewise-constant (readPos - bbPos) offset, as a sorted list of
    // (bbPos, deltaFromHere) breakpoints: "starting at this bbPos
    // (inclusive), delta = deltaFromHere, until the next breakpoint". Built
    // once from the FULL, exhaustive indel event list (every insertion/
    // deletion the aligner recorded across the member's whole covered
    // range) -- NOT from sparse shared-anchor pins. This matters: pins are
    // only as dense as the anchor graph (can be many kb apart), so "nearest
    // pin, reject if any indel lies between it and the target column" was
    // rejecting almost every column whenever an unrelated indel happened to
    // sit anywhere in that wide span, even though the aligner's own CIGAR
    // already gives the EXACT delta at every position with no ambiguity.
    // A deletion at bbPos=p, length L (backbone consumes L bases the member
    // doesn't have) shifts delta by -L starting at p+L (positions in
    // [p,p+L) are handled by deletionRanges instead, not by delta at all).
    // An insertion at bbPos=p, length L (member has L extra bases the
    // backbone doesn't) shifts delta by +L starting at p (no backbone
    // position is skipped, so the very next backbone column already sees
    // the new delta).
    std::vector<std::pair<std::uint32_t, std::int64_t>> deltaBreakpoints;
    std::int64_t deltaAtStart = 0; // delta for bbPos < first breakpoint (or all of bbCovBegin..bbCovEnd if none)
    // Locally dense mismatch/indel clusters (IpoaNoiseTracker output): this
    // member's votes inside these ranges are untrustworthy and excluded from
    // both OnBackbone and Substitution classification, matching the same
    // exclusion emitHetBubblesFromProfiles's memberCall applies to
    // KwMemberProfile::isNoisy.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> noisyRanges;
};

// Populates state/altBase from SPARSE per-member evidence: every position in
// [bbCovBegin, bbCovEnd) defaults to OnBackbone (the aligner didn't bother
// recording a match explicitly), EXCEPT positions inside a deletionRange
// (Invalid -- no base there) or listed in snps (Substitution). This is the
// mirror image of findLeafSnarlsFromProfiles's default (Invalid unless an
// explicit alignedCols entry says otherwise): correct here because a sparse
// diff-only aligner recording NOTHING at a position specifically means "it
// matched", not "unknown".
inline std::vector<LeafSnarl> findLeafSnarlsFromSparseEvidence(
    OrientedReadId backboneOid,
    const std::vector<SparseMemberEvidence>& evidence,
    const std::vector<std::uint8_t>& bbBases,
    std::uint32_t bbBeginAbs,
    std::uint64_t minSupport,
    double minVaf,
    bool dropHomopolymer,
    bool dropRepeat)
{
    using std::uint8_t;
    using std::uint32_t;
    using std::vector;

    const uint32_t n = uint32_t(bbBases.size());
    if (n == 0) return {};

    const int nMembers = int(evidence.size());
    const size_t nMembersSz = size_t(nMembers);

    vector<IpoaMemberInfo> members(nMembersSz);
    for (int r = 0; r < nMembers; r++) {
        members[size_t(r)].oid = evidence[size_t(r)].oid;
        members[size_t(r)].bbCovBegin = evidence[size_t(r)].bbCovBegin;
        members[size_t(r)].bbCovEnd = evidence[size_t(r)].bbCovEnd;
    }

    auto idx = [nMembersSz](uint32_t i, int r) { return size_t(i) * nMembersSz + size_t(r); };
    vector<IpoaColState> state(size_t(n) * nMembersSz, IpoaColState::NotCovered);
    vector<uint8_t> altBase(size_t(n) * nMembersSz, 0);

    for (int r = 0; r < nMembers; r++) {
        const SparseMemberEvidence& ev = evidence[size_t(r)];
        const uint32_t covBegin = std::max(ev.bbCovBegin, bbBeginAbs);
        const uint32_t covEnd = std::min(ev.bbCovEnd, bbBeginAbs + n);
        for (uint32_t pos = covBegin; pos < covEnd; pos++) {
            state[idx(pos - bbBeginAbs, r)] = IpoaColState::OnBackbone;
        }
        for (const auto& del : ev.deletionRanges) {
            const uint32_t a = std::max(del.first, covBegin);
            const uint32_t b = std::min(del.second, covEnd);
            for (uint32_t pos = a; pos < b; pos++) {
                state[idx(pos - bbBeginAbs, r)] = IpoaColState::Invalid;
            }
        }
        // Noisy ranges: untrustworthy regardless of what's underneath (a
        // match default OR a deletion), matching memberCall's precedence
        // (isNoisy checked unconditionally, before looking at snps).
        for (const auto& noisy : ev.noisyRanges) {
            const uint32_t a = std::max(noisy.first, covBegin);
            const uint32_t b = std::min(noisy.second, covEnd);
            for (uint32_t pos = a; pos < b; pos++) {
                state[idx(pos - bbBeginAbs, r)] = IpoaColState::Invalid;
            }
        }
        for (const SparseSnp& snp : ev.snps) {
            if (snp.bbPos < covBegin || snp.bbPos >= covEnd) continue;
            if (snp.altBase >= 4) continue; // N: leave whatever it was, not a substitution.
            const size_t cell = idx(snp.bbPos - bbBeginAbs, r);
            if (state[cell] == IpoaColState::Invalid) continue; // deleted or noisy: not a substitution.
            state[cell] = IpoaColState::Substitution;
            altBase[cell] = snp.altBase;
        }
    }

    return classifyLeafSnarls(
        n, nMembers, state, altBase, members, backboneOid, bbBases, bbBeginAbs,
        minSupport, minVaf, dropHomopolymer, dropRepeat);
}

// Sparse-native equivalent of findLeafSnarlsFromSparseEvidence: produces the
// IDENTICAL LeafSnarl list (validated by DINARA_SPARSE_CLASSIFY_VALIDATE=1,
// which runs both and asserts agreement) WITHOUT ever materializing the dense
// n*nMembers state/altBase array -- that array, and classifyLeafSnarls's own
// per-column scan of it, are the "~26% of runtime" a sparse-diff aligner like
// ProjectedAlignment doesn't actually need: the vast majority of positions in
// a wide window are OnBackbone for every member (nothing was ever recorded
// there), so walking every (position, member) pair to discover that is pure
// waste when the sparse evidence already tells you directly where the real
// divergence is.
//
// Key insight: NONE of the three things classifyLeafSnarls' dense state array
// is used for actually needs a per-column scan when the source is sparse:
//   1. Run boundaries (divergedSpan/divergedEqual) only depend on
//      Substitution/Invalid -- i.e. exactly the recorded snps/deletionRanges/
//      noisyRanges. A member's diverged-or-not status can only CHANGE at the
//      start/end of one of these recorded spans, so sweeping just those
//      breakpoints (not every column) finds the exact same maximal
//      constant-diverged-set runs classifyLeafSnarls's column-by-column
//      divergedEqual walk does.
//   2. rowExcluded's NotCovered half is a pure bounds check
//      ([runStart,runEnd) subset of [bbCovBegin,bbCovEnd)?), not a per-column
//      state read. Its Invalid half is answered by checking the (typically
//      tiny) deletionRanges/noisyRanges lists for overlap with the run --
//      also not a per-column read.
//   3. keyFor only needs a base per (run column, surviving member) -- and
//      every surviving member (rowExcluded already dropped the rest) is
//      either OnBackbone (no recorded snp at that column: use bbBases) or
//      Substitution (a recorded snp: use its altBase), decided by a lookup
//      against that member's own small snps list, not a dense array.
// flanksLinear (checking 8 flank columns) is the one place this still reads
// per-member per-column state, but only for the handful of runs actually
// found -- not for the window's silent majority.
inline std::vector<LeafSnarl> findLeafSnarlsFromSparseEvidenceFast(
    OrientedReadId backboneOid,
    const std::vector<SparseMemberEvidence>& evidence,
    const std::vector<std::uint8_t>& bbBases,
    std::uint32_t bbBeginAbs,
    std::uint64_t minSupport,
    double minVaf,
    bool dropHomopolymer,
    bool dropRepeat)
{
    using std::uint8_t;
    using std::uint32_t;
    using std::uint64_t;
    using std::vector;

    vector<LeafSnarl> result;
    const uint32_t n = uint32_t(bbBases.size());
    if (n == 0) return result;
    const int nMembers = int(evidence.size());

    // Sparse per-member state at an ABSOLUTE backbone position: only called
    // for flanksLinear's 8 flank columns and for keyFor's run columns, never
    // for the window's bulk. altBase is set only when the return is
    // Substitution.
    // Binary search, not linear scan: deletionRanges/noisyRanges/snps are
    // each sorted (and disjoint, for the range lists) by construction, but
    // are also now LONG for a CIGAR-reuse-sourced member (its wider
    // per-member coverage means 50-190+ entries per list is common) --
    // stateAt is called 8*nMembers times per candidate run inside
    // flanksLinear alone, so a linear scan here was the actual dominant
    // cost of the whole detector once CIGAR reuse widened coverage
    // (measured: ~60% of runtime), NOT the diverged-set bookkeeping.
    auto rangeContains = [](const vector<pair<uint32_t, uint32_t>>& ranges, uint32_t pos) {
        auto it = std::upper_bound(ranges.begin(), ranges.end(), pos,
            [](uint32_t p, const pair<uint32_t, uint32_t>& r) { return p < r.first; });
        if (it == ranges.begin()) return false;
        --it;
        return pos >= it->first && pos < it->second;
    };
    auto stateAt = [&](int r, uint32_t absPos, uint8_t& outAltBase) -> IpoaColState {
        const SparseMemberEvidence& ev = evidence[size_t(r)];
        if (absPos < ev.bbCovBegin || absPos >= ev.bbCovEnd) return IpoaColState::NotCovered;
        if (rangeContains(ev.deletionRanges, absPos)) return IpoaColState::Invalid;
        if (rangeContains(ev.noisyRanges, absPos)) return IpoaColState::Invalid;
        auto sit = std::lower_bound(ev.snps.begin(), ev.snps.end(), absPos,
            [](const SparseSnp& s, uint32_t p) { return s.bbPos < p; });
        if (sit != ev.snps.end() && sit->bbPos == absPos && sit->altBase < 4) {
            outAltBase = sit->altBase;
            return IpoaColState::Substitution;
        }
        return IpoaColState::OnBackbone;
    };

    auto flanksLinear = [&](uint32_t spanStartAbs, uint32_t spanEndAbs) -> bool {
        if (spanStartAbs < bbBeginAbs + 4 || spanEndAbs + 4 > bbBeginAbs + n) return false;
        for (const uint32_t c : {spanStartAbs - 4, spanStartAbs - 3, spanStartAbs - 2, spanStartAbs - 1,
                                  spanEndAbs,       spanEndAbs + 1,   spanEndAbs + 2,   spanEndAbs + 3}) {
            uint64_t invalidCount = 0;
            uint64_t altCounts[4] = {0, 0, 0, 0};
            for (int r = 0; r < nMembers; r++) {
                uint8_t ab = 0;
                const IpoaColState s = stateAt(r, c, ab);
                if (s == IpoaColState::Invalid) invalidCount++;
                else if (s == IpoaColState::Substitution) altCounts[ab]++;
            }
            if (invalidCount >= minSupport) return false;
            for (uint64_t cnt : altCounts) if (cnt >= minSupport) return false;
        }
        return true;
    };

    // Breakpoint events: (position, member, +1 at the start of a
    // Substitution/Invalid span for that member, -1 one past its end).
    // Absolute backbone coordinates, clipped to [bbBeginAbs, bbBeginAbs+n).
    struct BEv { uint32_t pos; int member; int delta; };
    vector<BEv> events;
    events.reserve(nMembers * 2);
    for (int r = 0; r < nMembers; r++) {
        const SparseMemberEvidence& ev = evidence[size_t(r)];
        for (const SparseSnp& s : ev.snps) {
            if (s.altBase >= 4) continue;
            if (s.bbPos < bbBeginAbs || s.bbPos >= bbBeginAbs + n) continue;
            events.push_back({s.bbPos, r, +1});
            events.push_back({s.bbPos + 1, r, -1});
        }
        for (const auto& d : ev.deletionRanges) {
            const uint32_t a = std::max(d.first, bbBeginAbs);
            const uint32_t b = std::min(d.second, bbBeginAbs + n);
            if (a < b) { events.push_back({a, r, +1}); events.push_back({b, r, -1}); }
        }
        for (const auto& nr : ev.noisyRanges) {
            const uint32_t a = std::max(nr.first, bbBeginAbs);
            const uint32_t b = std::min(nr.second, bbBeginAbs + n);
            if (a < b) { events.push_back({a, r, +1}); events.push_back({b, r, -1}); }
        }
    }
    if (events.empty()) return result;
    // Removals sort before additions at the SAME position (delta -1 < +1):
    // a member whose divergence span ends exactly where another of its own
    // spans begins (e.g. two adjacent single-column snps) must stay
    // continuously in divergedSet with no one-position gap. Applying the add
    // before the remove would incorrectly drop it for that one position
    // (addMember is a no-op when already present, so remove-after-add still
    // removes it), splitting or losing a run that should span both columns.
    std::sort(events.begin(), events.end(),
        [](const BEv& x, const BEv& y) {
            if (x.pos != y.pos) return x.pos < y.pos;
            return x.delta < y.delta;
        });

    // Live diverged-member set, updated as the sweep crosses each
    // breakpoint. Represented as a fixed-width bitset (vector<uint64_t>,
    // one bit per member) rather than a sorted vector<int>: with dense
    // per-member evidence (e.g. CIGAR-reuse's wider coverage feeding many
    // more events per member) and nMembers in the hundreds, the sorted
    // vector's O(nMembers) insert/erase/compare per breakpoint measurably
    // dominated runtime (findLeafSnarlsFromSparseEvidenceFast's classify
    // phase went from ~10% to ~60% of this detector's time once CIGAR
    // reuse widened per-member coverage). A bitset turns set/clear into
    // O(1) (single word, single bit) and equality/emptiness into O(words)
    // word-at-a-time comparisons -- no per-member work at all. Between two
    // consecutive distinct event positions this set is constant by
    // construction; a run (below) is a MAXIMAL span of such constant-set
    // positions, possibly crossing several breakpoints that net out to the
    // same set (see the note on runOpen).
    const size_t nWords = (size_t(nMembers) + 63) / 64;
    vector<uint64_t> divergedSet(nWords, 0);
    auto addMember = [&](int r) { divergedSet[size_t(r) / 64] |= (uint64_t(1) << (size_t(r) % 64)); };
    auto removeMember = [&](int r) { divergedSet[size_t(r) / 64] &= ~(uint64_t(1) << (size_t(r) % 64)); };
    auto bitsetPopcount = [](const vector<uint64_t>& b) {
        uint64_t c = 0;
        for (uint64_t w : b) c += uint64_t(__builtin_popcountll(w));
        return c;
    };

    // A run is [openRunStart, someLaterPos) with a FIXED diverged-member set
    // (openRunSet, snapshotted when the run opened). Breakpoints only close
    // the run when the resulting set actually DIFFERS from openRunSet -- two
    // breakpoints in a row that net out to the same set (e.g. one member's
    // divergence ends exactly where a DIFFERENT member's begins, touching
    // but not changing the overall set) must NOT split what
    // classifyLeafSnarls' column-by-column divergedEqual would see as one
    // unbroken run. Splitting it would additionally make each half look
    // divergent to the other's flanksLinear check, so both fragments would
    // spuriously fail flank-linearity and the whole real run would vanish.
    bool runOpen = false;
    uint32_t openRunStart = 0;
    vector<uint64_t> openRunSet;

    auto processRun = [&](uint32_t runStartAbs, uint32_t runEndAbs) {
        if (runStartAbs >= runEndAbs) return;
        if (!flanksLinear(runStartAbs, runEndAbs)) return;

        const uint32_t runStart = runStartAbs - bbBeginAbs; // local, for bbBases indexing
        const uint32_t runEnd = runEndAbs - bbBeginAbs;

        // rowExcluded: NotCovered is a pure coverage-bounds check; Invalid is
        // a small-list overlap check -- neither needs a per-column scan.
        vector<bool> rowExcluded(size_t(nMembers), false);
        for (int r = 0; r < nMembers; r++) {
            const SparseMemberEvidence& ev = evidence[size_t(r)];
            if (ev.bbCovBegin > runStartAbs || ev.bbCovEnd < runEndAbs) {
                rowExcluded[size_t(r)] = true;
                continue;
            }
            for (const auto& d : ev.deletionRanges)
                if (d.first < runEndAbs && d.second > runStartAbs) { rowExcluded[size_t(r)] = true; break; }
            if (rowExcluded[size_t(r)]) continue;
            for (const auto& nr : ev.noisyRanges)
                if (nr.first < runEndAbs && nr.second > runStartAbs) { rowExcluded[size_t(r)] = true; break; }
        }

        auto keyFor = [&](int r) {
            vector<uint8_t> k;
            k.reserve(runEnd - runStart);
            for (uint32_t c = runStart; c < runEnd; c++) {
                uint8_t ab = 0;
                const IpoaColState s = stateAt(r, bbBeginAbs + c, ab);
                k.push_back(s == IpoaColState::Substitution ? ab : bbBases[c]);
            }
            return k;
        };

        vector<vector<uint8_t>> alleleKeys;
        vector<LeafSnarlAllele> alleles;
        {
            const vector<uint8_t> refKey(bbBases.begin() + runStart, bbBases.begin() + runEnd);
            alleleKeys.push_back(refKey);
            LeafSnarlAllele refAllele;
            refAllele.members.push_back(backboneOid);
            refAllele.bases = refKey;
            alleles.push_back(std::move(refAllele));
        }
        const int refAlleleIndex = 0;

        for (int r = 0; r < nMembers; r++) {
            if (rowExcluded[size_t(r)]) continue;
            vector<uint8_t> k = keyFor(r);
            size_t ai = alleleKeys.size();
            for (size_t j = 0; j < alleleKeys.size(); j++) {
                if (alleleKeys[j] == k) { ai = j; break; }
            }
            if (ai == alleleKeys.size()) {
                alleleKeys.push_back(k);
                LeafSnarlAllele allele;
                allele.bases = k;
                alleles.push_back(std::move(allele));
            }
            alleles[ai].members.push_back(evidence[size_t(r)].oid);
        }

        uint64_t totalMembers = 0;
        for (const auto& a : alleles) totalMembers += a.members.size();

        vector<LeafSnarlAllele> gated;
        for (int ai = 0; ai < int(alleles.size()); ai++) {
            LeafSnarlAllele& allele = alleles[size_t(ai)];
            if (allele.members.size() < minSupport) continue;
            if (ai != refAlleleIndex) {
                const double af = double(allele.members.size()) / double(totalMembers);
                if (af < minVaf) continue;
            }
            gated.push_back(std::move(allele));
        }

        bool passesRepeatGate = true;
        if (dropHomopolymer || dropRepeat) {
            KmVarKey vkey;
            vkey.type = KmVarType::Snp;
            vkey.refLen = 1; vkey.altLen = 1; vkey.altBase = 0;
            for (uint32_t localPos : {runStart, runEnd - 1}) {
                vkey.pos = localPos;
                const bool inHomopolymer = kmIsRepeatUnitRange(
                    bbBases.data(), uint32_t(bbBases.size()), vkey, 0, 1, 1);
                const bool inStr = kmIsRepeatUnitRange(
                    bbBases.data(), uint32_t(bbBases.size()), vkey, 0, 2, 6);
                if (dropHomopolymer && inHomopolymer) passesRepeatGate = false;
                if (dropRepeat && inStr) passesRepeatGate = false;
            }
        }

        if (passesRepeatGate && gated.size() >= 2) {
            LeafSnarl snarl;
            snarl.start = runStartAbs - 1;
            snarl.end = runEndAbs;
            snarl.alleles = std::move(gated);
            result.push_back(std::move(snarl));
        }
    };

    size_t ei = 0;
    while (ei < events.size()) {
        const uint32_t pos = events[ei].pos;
        size_t ej = ei;
        while (ej < events.size() && events[ej].pos == pos) {
            if (events[ej].delta > 0) addMember(events[ej].member);
            else removeMember(events[ej].member);
            ej++;
        }
        ei = ej;

        if (runOpen && divergedSet == openRunSet) continue; // touching, unchanged set: keep the run open

        // Every alt allele is built EXCLUSIVELY from divergedSet members (a
        // recorded snp always differs from the backbone base, by
        // construction -- see ProjectedAlignment's mismatch-only recording
        // -- so a member outside divergedSet can only ever key into the ref
        // allele). An alt allele's member count is therefore bounded by
        // openRunSet.size(), so openRunSet.size() < minSupport makes gating
        // (>=2 alleles, each >=minSupport) unreachable regardless of what
        // flanksLinear/keyFor would say -- skip the whole run for free. This
        // matters: real HiFi windows carry plenty of ordinary per-base
        // sequencing noise (isolated single-member "divergence"), and
        // without this check EVERY such blip still pays flanksLinear's
        // 8*nMembers sparse lookups before minSupport ever gets to reject it.
        if (runOpen) {
            if (bitsetPopcount(openRunSet) >= minSupport) processRun(openRunStart, pos);
            runOpen = false;
        }
        if (bitsetPopcount(divergedSet) >= minSupport) {
            runOpen = true;
            openRunStart = pos;
            openRunSet = divergedSet;
        }
    }
    if (runOpen && bitsetPopcount(openRunSet) >= minSupport) processRun(openRunStart, bbBeginAbs + n);

    return result;
}

// Result of a sparse column lookup: the member's exact read position and
// base at a given backbone column. Mirrors KwAlignedCol's role for the dense
// aligners, but computed from sparse evidence instead of a stored list.
struct SparseColResult {
    std::uint32_t readPos;
    std::uint8_t readBase;
};

// The sparse equivalent of KwMemberProfile::colAt, needed to PIN a het/hom
// anchor (recover the member's exact read position at a specific backbone
// column) from sparse evidence. Three cases:
//   1. bbPos outside coverage, or inside a deletion: not resolvable (no base
//      there) -- returns false, matching colAt's nullptr.
//   2. bbPos is a recorded mismatch: exact, no interpolation needed.
//   3. bbPos is an ordinary (unrecorded) match column: readPos = bbPos +
//      delta, where delta is looked up from ev.deltaBreakpoints (the
//      piecewise-constant offset built from the aligner's OWN exhaustive
//      indel list -- see SparseMemberEvidence's comment for why this
//      replaces nearest-pin interpolation). No adjacency/indel-in-between
//      check is needed here: the breakpoints already encode the exact
//      cumulative effect of every indel over the whole covered range, so
//      the delta this returns is exact, not an approximation.
// bbBaseAtPos is the backbone's own base at bbPos, used only to fill
// out.readBase for the interpolated (match) case -- a match column's base is
// the backbone's own base there BY DEFINITION under the sparse "unrecorded
// = matches" model, so this keeps the function's contract identical to
// colAt's (always a real base, never a sentinel).
inline bool sparseColAt(
    const SparseMemberEvidence& ev,
    std::uint32_t bbPos,
    std::uint8_t bbBaseAtPos,
    SparseColResult& out)
{
    if (bbPos < ev.bbCovBegin || bbPos >= ev.bbCovEnd) return false;
    // Binary search, not linear scan: both lists are sorted (and
    // deletionRanges is disjoint) by construction, and can be long for a
    // CIGAR-reuse-sourced member (wider per-member coverage means more
    // recorded events) -- this is called once per member per anchor column
    // in buildHetBubbleFromLeafSnarl, so a linear scan here was a real,
    // measurable cost (same class of issue found and fixed in stateAt,
    // WindowIntervalPoa.hpp's classify path).
    {
        auto dit = std::upper_bound(ev.deletionRanges.begin(), ev.deletionRanges.end(), bbPos,
            [](std::uint32_t p, const std::pair<std::uint32_t, std::uint32_t>& r) { return p < r.first; });
        if (dit != ev.deletionRanges.begin()) {
            const auto& d = *std::prev(dit);
            if (bbPos >= d.first && bbPos < d.second) return false;
        }
    }

    {
        auto sit = std::lower_bound(ev.snps.begin(), ev.snps.end(), bbPos,
            [](const SparseSnp& s, std::uint32_t p) { return s.bbPos < p; });
        if (sit != ev.snps.end() && sit->bbPos == bbPos) {
            out.readPos = sit->readPos;
            out.readBase = sit->altBase;
            return true;
        }
    }

    // Ordinary match column: delta = the last breakpoint at or before bbPos
    // (or deltaAtStart if bbPos precedes every breakpoint).
    std::int64_t delta = ev.deltaAtStart;
    auto it = std::upper_bound(ev.deltaBreakpoints.begin(), ev.deltaBreakpoints.end(), bbPos,
        [](std::uint32_t pos, const std::pair<std::uint32_t, std::int64_t>& bp) {
            return pos < bp.first;
        });
    if (it != ev.deltaBreakpoints.begin()) delta = std::prev(it)->second;

    const std::int64_t readPos = std::int64_t(bbPos) + delta;
    if (readPos < 0) return false;
    out.readPos = std::uint32_t(readPos);
    out.readBase = bbBaseAtPos;
    return true;
}

// Sparse equivalent of pinnedKmerCol (WindowHetProfiles.hpp): valid only if
// the member spells the 2-mer [kmer0, kmer1] starting at predPos with no gap
// between the two bases (readPos advances by exactly 1). bbBaseAtPredPos/
// bbBaseAtPredPos1 are the backbone's own bases at predPos and predPos+1
// (passed through to sparseColAt).
inline bool sparsePinnedKmerCol(
    const SparseMemberEvidence& ev,
    std::uint32_t predPos, std::uint8_t kmer0, std::uint8_t kmer1,
    std::uint8_t bbBaseAtPredPos, std::uint8_t bbBaseAtPredPos1,
    std::uint32_t& outReadPos)
{
    SparseColResult c0;
    if (!sparseColAt(ev, predPos, bbBaseAtPredPos, c0)) return false;
    if (c0.readBase != kmer0) return false;

    SparseColResult c1;
    if (!sparseColAt(ev, predPos + 1, bbBaseAtPredPos1, c1)) return false;
    if (c1.readPos != c0.readPos + 1) return false; // insertion between them
    if (c1.readBase != kmer1) return false;

    outReadPos = c0.readPos;
    return true;
}

// Sparse equivalent of pinnedPointCol: valid only if the member carries
// `base` at atPos. bbBaseAtPos is the backbone's own base at atPos (passed
// through to sparseColAt).
inline bool sparsePinnedPointCol(
    const SparseMemberEvidence& ev,
    std::uint32_t atPos, std::uint8_t base, std::uint8_t bbBaseAtPos,
    std::uint32_t& outReadPos)
{
    SparseColResult c;
    if (!sparseColAt(ev, atPos, bbBaseAtPos, c)) return false;
    if (c.readBase != base) return false;
    outReadPos = c.readPos;
    return true;
}

// Shared MSA + emit for a backbone interval. Given `spanMembers` already placed
// over [bbBegin, bbEnd), run ONE abPOA MSA with row 0 = backbone segment and
// emit each member's per-column contribution into `frag`. Split out from
// runIpoaInterval so the member-selection and the MSA/extraction concerns are
// separable. Requires >=2 two-sided members. Safe to call from any thread.
//
// minMembers is a per-position support pre-filter: a het at a backbone position
// needs both a ref allele (>=minSupport members, counting the backbone) and an
// alt allele (>=minSupport members) over DISJOINT member sets, so a passing
// position needs >= 2*minSupport-1 covering members. Every member here spans the
// whole interval, so spanMembers.size() is an upper bound on per-position
// coverage: if it is below minMembers no position in this interval can ever
// pass, and the abPOA MSA is skipped. This never drops a real het (it is a
// strict upper bound) and never affects the backbone chain (a skipped interval
// simply holds no bubble -> a plain bbA_i->bbA_i+1 edge is staged). Pass
// minMembers<=0 to disable the pre-filter.
inline void runIpoaOnRows(
    const IpoaPlan& plan,
    std::uint32_t bbBegin,
    std::uint32_t bbEnd,
    std::vector<IpoaSpanMember>& spanMembers,
    const Reads& rds,
    IpoaAbHandle& ah,
    IpoaFragment& frag,
    bool timing,
    const IpoaClock::time_point& tSetup0,
    int minMembers)
{
    using std::vector;
    using std::uint8_t;
    using std::uint32_t;
    using std::uint64_t;
    using std::int64_t;

    if (spanMembers.empty()) return;

    // Per-position support pre-filter: skip the MSA when no position could pass.
    // Kill switch (DINARA_IPOA_NO_SKIP=1) forces every interval to run, for
    // measuring the filter's speedup and confirming it changes no output.
    static const bool noSkip = (std::getenv("DINARA_IPOA_NO_SKIP") != nullptr);
    if (!noSkip && minMembers > 0 && spanMembers.size() < std::size_t(minMembers)) {
        if (timing) ipoaTiming().skippedIntervals.fetch_add(
            1, std::memory_order_relaxed);
        return;
    }

    const uint32_t segLen = bbEnd - bbBegin;

    // Order rows: two-sided first (they anchor the POA), then one-sided.
    vector<const IpoaSpanMember*> ordered;
    ordered.reserve(spanMembers.size());
    for (const IpoaSpanMember& sm : spanMembers) if (sm.side == 0) ordered.push_back(&sm);
    const uint32_t nTwoSided = uint32_t(ordered.size());
    for (const IpoaSpanMember& sm : spanMembers) if (sm.side != 0) ordered.push_back(&sm);
    if (nTwoSided < 2) return;   // need >=2 two-sided reads to anchor the POA

    auto toCodes = [&](const OrientedReadId o, uint32_t cB, uint32_t cE,
                       vector<uint8_t>& out) {
        out.resize(cE - cB);
        for (uint32_t i = cB; i < cE; i++)
            out[i - cB] = uint8_t(rds.getOrientedReadBase(o, i).value & 3);
    };

    // Assemble the input sequence set: row 0 = backbone segment, rows 1..N =
    // ordered members. Thread-local scratch avoids per-interval reallocation.
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
        for (uint32_t i = 0; i < segLen; i++)
            bb[i] = uint8_t(rds.getOrientedReadBase(plan.bbOid, bbBegin + i).value & 3);
        seqPtrs[0] = bb.data();
        seqLens[0] = int(segLen);
    }
    bool addFailed = false;
    for (uint32_t oi = 0; oi < ordered.size(); oi++) {
        const IpoaSpanMember& sm = *ordered[oi];
        auto& codes = seqStore[oi + 1];
        toCodes(sm.oid, sm.cBegin, sm.cEnd, codes);
        if (codes.empty()) { addFailed = true; break; }
        seqPtrs[oi + 1] = codes.data();
        seqLens[oi + 1] = int(codes.size());
    }
    if (addFailed) return;

    int maxLen = 0;
    for (int r = 0; r < nSeq; r++) if (seqLens[r] > maxLen) maxLen = seqLens[r];

    // Length guard against abpoa's monolithic DP blowup. When reads over an
    // anchor-less span (a divergent repeat / SV / low-coverage gap) share no
    // minimizers, abpoa_anchor_poa cannot partition the interval and falls back
    // to one O(node_n * qlen) DP matrix; adaptive banding does not bound it. At
    // genome scale a single such interval demanded a 128 GiB SIMDMalloc and OOM-
    // killed the run. maxLen (the longest row) upper-bounds both DP dimensions,
    // so worst-case bytes ~= maxLen^2 * ~32. Cap maxLen so that stays within a
    // safe budget (default 16 kb -> ~8 GiB peak per thread); pathological long
    // intervals are skipped rather than crashing the whole assembly. Override
    // with DINARA_IPOA_MAX_LEN=0 to disable, or a custom cap in bases.
    static const int maxLenCap = [] {
        const char* e = std::getenv("DINARA_IPOA_MAX_LEN");
        if (e == nullptr) return 16384;
        const int v = std::atoi(e);
        return v;   // 0 disables the guard; negative treated as disabled below
    }();
    if (maxLenCap > 0 && maxLen > maxLenCap) {
        if (timing) ipoaTiming().skippedByLen.fetch_add(
            1, std::memory_order_relaxed);
        return;
    }

    if (timing) {
        ipoaTiming().setupNs.fetch_add(ipoaNsSince(tSetup0),
                                       std::memory_order_relaxed);
    }
    const IpoaClock::time_point tReset0 =
        timing ? IpoaClock::now() : IpoaClock::time_point{};

    abpoa_reset(ah.ab, ah.abpt, maxLen > 0 ? maxLen : 1);

    const IpoaClock::time_point tMsa0 =
        timing ? IpoaClock::now() : IpoaClock::time_point{};

    abpoa_msa(ah.ab, ah.abpt, nSeq, nullptr, seqLens.data(),
              seqPtrs.data(), nullptr, nullptr);

    if (timing) {
        IpoaTiming& tt = ipoaTiming();
        const uint64_t resetDt = ipoaNsSince(tReset0);   // includes msa below? no
        const uint64_t msaDt = ipoaNsSince(tMsa0);
        // resetDt above spans reset+msa; recompute reset-only as (tMsa0-tReset0).
        const uint64_t resetOnly = uint64_t(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                tMsa0 - tReset0).count());
        const uint64_t runDt = resetOnly + msaDt;        // reset + msa
        tt.resetNs.fetch_add(resetOnly, std::memory_order_relaxed);
        tt.msaNs.fetch_add(msaDt, std::memory_order_relaxed);
        tt.msaRuns.fetch_add(1, std::memory_order_relaxed);
        tt.seqSum.fetch_add(uint64_t(nSeq), std::memory_order_relaxed);
        ipoaAtomicMax(tt.seqMax, uint64_t(nSeq));
        tt.lenSum.fetch_add(uint64_t(maxLen), std::memory_order_relaxed);
        ipoaAtomicMax(tt.lenMax, uint64_t(maxLen));
        if (ah.ab->abc && ah.ab->abc->msa_len > 0) {
            tt.colSum.fetch_add(uint64_t(ah.ab->abc->msa_len),
                                std::memory_order_relaxed);
            ipoaAtomicMax(tt.colMax, uint64_t(ah.ab->abc->msa_len));
        }
        ipoaAtomicMax(tt.runMaxNs, runDt);
        if (runDt > 100000)   tt.runOver100us.fetch_add(1, std::memory_order_relaxed);
        if (runDt > 1000000) {
            tt.runOver1ms.fetch_add(1, std::memory_order_relaxed);
            tt.nsOver1ms.fetch_add(runDt, std::memory_order_relaxed);
        }
        if (runDt > 10000000) tt.runOver10ms.fetch_add(1, std::memory_order_relaxed);
        (void)resetDt;
    }
    const IpoaClock::time_point tExtract0 =
        timing ? IpoaClock::now() : IpoaClock::time_point{};

    const abpoa_cons_t* abc = ah.ab->abc;
    if (abc == nullptr || abc->msa_len <= 0 || abc->n_seq != nSeq) {
        if (timing) ipoaTiming().extractNs.fetch_add(ipoaNsSince(tExtract0),
                                                     std::memory_order_relaxed);
        return;
    }
    frag.ran = true;

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

    // Column -> absolute backbone position (row 0 = backbone).
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

    frag.reads.reserve(ordered.size());
    for (uint32_t si = 1; si < alignment.size(); si++) {
        const IpoaSpanMember& sm = *ordered[si - 1];
        const auto& row = alignment[si];

        IpoaReadFrag rf;
        rf.oidValue = sm.oid.getValue();

        int64_t colFirst = -1, colLast = -1;
        for (uint64_t col = 0; col < ncols; col++) {
            if (!row[col].isGap()) {
                if (colFirst < 0) colFirst = int64_t(col);
                colLast = int64_t(col);
            }
        }
        if (colFirst < 0) continue;

        uint32_t readAbs = sm.cBegin;
        for (int64_t col = 0; col < colFirst; col++)
            if (!row[col].isGap()) readAbs++;

        int64_t pendingDelBegin = -1;
        bool pendingInsertion = false;
        bool entered = false;

        // Bound the read position by the member's segment [cBegin, cEnd) on BOTH
        // sides, whether each boundary is a TRUE pin or a guessed one-sided end.
        //
        // On a one-sided read the free end (cEnd for left-only, cBegin for
        // right-only) is a base-count GUESS (cB + segLenApprox). Indels make the
        // read's true base count across the span differ from the backbone's, so
        // readAbs drifts past the read's true position at the next shared anchor.
        // A later het arm/hom anchor built from that drifted readPos then lands
        // beyond the read's real position at the neighbouring backbone anchor,
        // producing a backward intra-window edge (caught by monotonity
        // verification as an off-by-one, which is the overshoot minus the
        // primary anchor's k/2 frame offset).
        //
        // Clamping to the guessed boundary drops those uncertain overshoot
        // columns -- we lose a little coverage on the guessed end, never emit a
        // position we cannot trust. (This is what commit aeba3de intended; it
        // previously bounded only the true-pinned side, leaving the drift side
        // unchecked.)
        for (int64_t col = colFirst; col <= colLast; col++) {
            const AlignedBase mb = row[col];
            const int64_t bbPos = colBbPos[col];

            if (bbPos >= 0) {
                if (!mb.isGap()) {
                    const uint8_t code = mb.value & 0xff;
                    const bool posValid =
                        (readAbs >= sm.cBegin) && (readAbs < sm.cEnd);
                    if (!posValid) { entered = true; readAbs++; continue; }
                    rf.alignedCols.push_back(
                        KwAlignedCol{uint32_t(bbPos), readAbs, code});
                    if (rf.firstBb < 0) rf.firstBb = bbPos;
                    rf.lastBb = bbPos;
                    // Backbone base for the SNP check: index the interval's
                    // already-materialized backbone codes (row 0 = seqStore[0],
                    // built [bbBegin, bbEnd) and &3), NOT a per-column
                    // getOrientedReadBase call. That call in the innermost loop
                    // (per column, per read, per interval) was the dominant het
                    // cost after bbSeqVec was removed for the OOM fix.
                    const uint8_t bbBase = seqStore[0][uint32_t(bbPos) - bbBegin];
                    if (code < 4 && code != bbBase)
                        rf.snps.push_back(KwSnp{uint32_t(bbPos), code});
                    if (pendingDelBegin >= 0) {
                        rf.deletionRanges.push_back(
                            {uint32_t(pendingDelBegin), uint32_t(bbPos)});
                        pendingDelBegin = -1;
                    }
                    // An insertion accumulated in the backbone-gap columns just
                    // before this aligned base: attach it to this backbone pos.
                    if (pendingInsertion) {
                        rf.insertionSites.push_back(uint32_t(bbPos));
                        pendingInsertion = false;
                    }
                    entered = true;
                    readAbs++;
                } else {
                    if (entered && pendingDelBegin < 0)
                        pendingDelBegin = bbPos;
                }
            } else {
                if (!mb.isGap()) {
                    // Read base with no backbone column = inserted base. Flag it
                    // (only once we are inside the read's aligned span) so the
                    // next aligned backbone column records the insertion site.
                    if (entered) pendingInsertion = true;
                    readAbs++;
                }
            }
        }
        (void)pendingDelBegin;

        if (!rf.alignedCols.empty())
            frag.reads.push_back(std::move(rf));
    }

    if (timing) ipoaTiming().extractNs.fetch_add(ipoaNsSince(tExtract0),
                                                 std::memory_order_relaxed);
}

// --- Mode front-ends: choose spanMembers, then call runIpoaOnRows ------------

// INTERVAL MODE (baseline, unchanged behavior). One POA between adjacent
// breakpoints [bi, bi+1]. Two-sided members come from the endpoint intersection;
// one-sided members get a GUESSED free end (cB + segLenApprox). This is the
// proven path -- kept byte-for-byte so `mode=interval` reproduces prior output.
inline void runIpoaInterval(
    const IpoaPlan& plan,
    std::uint32_t bi,
    const Reads& rds,
    bool oneSidedEnabled,
    IpoaAbHandle& ah,
    IpoaFragment& frag,
    int minMembers = 0)
{
    using std::vector;
    using std::uint8_t;
    using std::uint32_t;
    using std::uint64_t;
    using std::min;
    using std::max;

    frag.reads.clear();
    frag.twoSided = 0;
    frag.oneSided = 0;
    frag.oneSidedPlaced = 0;
    frag.ran = false;
    frag.hadOneSided = false;

    const bool timing = ipoaTimingEnabled();
    if (timing) ipoaTiming().intervals.fetch_add(1, std::memory_order_relaxed);
    const IpoaClock::time_point tSetup0 =
        timing ? IpoaClock::now() : IpoaClock::time_point{};

    const uint32_t bbBegin = plan.breakpoints[bi];
    const uint32_t bbEnd = plan.breakpoints[bi + 1];
    if (bbEnd <= bbBegin) return;
    if (bbEnd <= plan.windowBbBegin || bbBegin >= plan.windowBbEnd) return;

    // Diagnostic-only interval span cap (default 0 = off; every region is
    // processed). Wide anchor-less gaps are slow and drift-prone, but the fix is
    // to sub-tile them, NOT to skip them -- skipping silently drops any het
    // inside the gap. This knob only exists to MEASURE that tail's cost.
    if (const uint32_t cap = ipoaMaxIntervalBp())
        if (bbEnd - bbBegin > cap) return;

    vector<IpoaSpanMember> spanMembers;
    uint64_t oneSidedHere = 0;
    const uint32_t segLenApprox = bbEnd - bbBegin;
    auto readLenOf = [&](uint32_t oidValue) -> uint32_t {
        const OrientedReadId o = OrientedReadId::fromValue(ReadId(oidValue));
        return uint32_t(rds.getRead(o.getReadId()).baseCount);
    };
    {
        const auto& lo = plan.atBreakpoint[bi];
        const auto& hi = plan.atBreakpoint[bi + 1];
        size_t a = 0, b = 0;
        auto addLeftOnly = [&](const IpoaMemberBp& e) {
            oneSidedHere++;
            if (!oneSidedEnabled) return;
            const uint32_t rl = readLenOf(e.oidValue);
            const uint32_t cB = e.cPos;
            if (cB >= rl) return;
            const uint32_t cE = min(rl, cB + max<uint32_t>(segLenApprox, 1));
            if (cE > cB)
                spanMembers.push_back(
                    {OrientedReadId::fromValue(ReadId(e.oidValue)), cB, cE, 1});
        };
        auto addRightOnly = [&](const IpoaMemberBp& e) {
            oneSidedHere++;
            if (!oneSidedEnabled) return;
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
    frag.twoSided = twoSidedHere;
    frag.oneSided = oneSidedHere;
    frag.oneSidedPlaced = spanMembers.size() - twoSidedHere;
    frag.hadOneSided = (oneSidedHere > 0);

    runIpoaOnRows(plan, bbBegin, bbEnd, spanMembers, rds, ah, frag,
                  timing, tSetup0, minMembers);
}

// Per-read accumulator. Same shape as the original in-function Accum. Filled by
// streaming each interval fragment in as it completes (accumulateIpoaFragment)
// so per-interval outputs never accumulate -- peak memory per window is just the
// final accumulators (~depth*span), matching the original per-window path.
struct IpoaAccum {
    OrientedReadId oid;
    std::vector<KwAlignedCol> alignedCols;
    std::vector<KwSnp> snps;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> deletionRanges;
    std::vector<std::uint32_t> insertionSites;
    std::int64_t firstBb = -1;
    std::int64_t lastBb = -1;
};

// Per-window streaming accumulator: the map plus its own mutex, so workers can
// fold interval fragments in concurrently without a whole-window serial pass.
struct IpoaWindowAccum {
    std::unordered_map<std::uint32_t, IpoaAccum> accums;
    std::mutex mutex;
    // Instrumentation, summed as fragments arrive.
    std::uint32_t intervalsRun = 0;
    std::uint32_t intervalsWithOneSided = 0;
    std::uint64_t twoSidedTotal = 0;
    std::uint64_t oneSidedTotal = 0;
    std::uint64_t oneSidedPlaced = 0;
};

// Lock-free fold of one interval fragment into a window accumulator. Use this
// when the caller owns `wa` exclusively (the per-window path: one thread runs
// all of a window's intervals serially into its own `wa`, so no other thread
// touches it). Avoids ~intervalCount() uncontended lock/unlock pairs per window.
inline void accumulateIpoaFragmentUnlocked(IpoaWindowAccum& wa, IpoaFragment& frag)
{
    if (frag.ran) wa.intervalsRun++;
    if (frag.hadOneSided) wa.intervalsWithOneSided++;
    wa.twoSidedTotal += frag.twoSided;
    wa.oneSidedTotal += frag.oneSided;
    wa.oneSidedPlaced += frag.oneSidedPlaced;
    for (IpoaReadFrag& rf : frag.reads) {
        IpoaAccum& acc = wa.accums[rf.oidValue];
        acc.oid = OrientedReadId::fromValue(ReadId(rf.oidValue));
        acc.alignedCols.insert(acc.alignedCols.end(),
            rf.alignedCols.begin(), rf.alignedCols.end());
        acc.snps.insert(acc.snps.end(), rf.snps.begin(), rf.snps.end());
        acc.deletionRanges.insert(acc.deletionRanges.end(),
            rf.deletionRanges.begin(), rf.deletionRanges.end());
        acc.insertionSites.insert(acc.insertionSites.end(),
            rf.insertionSites.begin(), rf.insertionSites.end());
        if (rf.firstBb >= 0 && (acc.firstBb < 0 || rf.firstBb < acc.firstBb))
            acc.firstBb = rf.firstBb;
        if (rf.lastBb > acc.lastBb) acc.lastBb = rf.lastBb;
    }
}

// Fold one completed interval fragment into a window accumulator, then leave the
// fragment ready to be cleared and reused. Thread-safe (locks the window); for
// callers that fold fragments from multiple threads into a shared `wa`.
inline void accumulateIpoaFragment(IpoaWindowAccum& wa, IpoaFragment& frag)
{
    std::lock_guard<std::mutex> lock(wa.mutex);
    if (frag.ran) wa.intervalsRun++;
    if (frag.hadOneSided) wa.intervalsWithOneSided++;
    wa.twoSidedTotal += frag.twoSided;
    wa.oneSidedTotal += frag.oneSided;
    wa.oneSidedPlaced += frag.oneSidedPlaced;
    for (IpoaReadFrag& rf : frag.reads) {
        IpoaAccum& acc = wa.accums[rf.oidValue];
        acc.oid = OrientedReadId::fromValue(ReadId(rf.oidValue));
        acc.alignedCols.insert(acc.alignedCols.end(),
            rf.alignedCols.begin(), rf.alignedCols.end());
        acc.snps.insert(acc.snps.end(), rf.snps.begin(), rf.snps.end());
        acc.deletionRanges.insert(acc.deletionRanges.end(),
            rf.deletionRanges.begin(), rf.deletionRanges.end());
        acc.insertionSites.insert(acc.insertionSites.end(),
            rf.insertionSites.begin(), rf.insertionSites.end());
        if (rf.firstBb >= 0 && (acc.firstBb < 0 || rf.firstBb < acc.firstBb))
            acc.firstBb = rf.firstBb;
        if (rf.lastBb > acc.lastBb) acc.lastBb = rf.lastBb;
    }
}

// Merge all interval fragments of one window into KwMemberProfiles, apply the
// same dedup/monotonicity/RC cleanups as the original single-window path, then
// run the shared emit tail. Serial; called once per window.
//
// The concatenation is order-independent: alignedCols and snps are sorted+
// deduped by position, deletionRanges are sorted+merged, and firstBb/lastBb are
// taken as min/max -- so it does not matter which interval (or thread) produced
// a given fragment or in what order fragments arrive.
inline std::uint32_t mergeAndEmitIpoaWindow(
    AnchorWindow& window,
    const IpoaPlan& plan,
    IpoaWindowAccum& wa,
    const Reads& rds,
    std::uint64_t coverageHet,
    double hetMinVaf,
    std::uint64_t hetMinSupport,
    bool hetDropHomopolymer,
    bool hetDropRepeat)
{
    using std::vector;
    using std::pair;
    using std::uint8_t;
    using std::uint32_t;
    using std::uint64_t;
    using std::unordered_map;
    using std::sort;
    using std::max;
    using std::min;
    using std::unique;

    window.hetBubbles.clear();

    // Fragments were streamed into wa.accums as intervals completed; consume it.
    unordered_map<uint32_t, IpoaAccum>& accums = wa.accums;

    // Flush accumulators into KwMemberProfile rows.
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

        prof.alignedCols = std::move(acc.alignedCols);
        sort(prof.alignedCols.begin(), prof.alignedCols.end(),
            [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos < b.bbPos; });
        prof.alignedCols.erase(
            unique(prof.alignedCols.begin(), prof.alignedCols.end(),
                [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos == b.bbPos; }),
            prof.alignedCols.end());

        // Read-position monotonicity: drop columns where readPos does not
        // strictly increase with bbPos, and any SNP at a dropped column.
        {
            vector<KwAlignedCol> kept;
            kept.reserve(prof.alignedCols.size());
            uint32_t lastReadPos = 0; bool have = false;
            unordered_map<uint32_t, uint8_t> keptBaseAt;
            for (const KwAlignedCol& c : prof.alignedCols) {
                if (have && c.readPos <= lastReadPos) continue;
                kept.push_back(c);
                keptBaseAt[c.bbPos] = c.readBase;
                lastReadPos = c.readPos;
                have = true;
            }
            prof.alignedCols.swap(kept);
            vector<KwSnp> keptSnps;
            keptSnps.reserve(prof.snps.size());
            for (const KwSnp& s : prof.snps) {
                auto it = keptBaseAt.find(s.bbPos);
                if (it != keptBaseAt.end() && it->second == s.altBase)
                    keptSnps.push_back(s);
            }
            prof.snps.swap(keptSnps);
        }

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

        if (!acc.insertionSites.empty()) {
            sort(acc.insertionSites.begin(), acc.insertionSites.end());
            acc.insertionSites.erase(
                unique(acc.insertionSites.begin(), acc.insertionSites.end()),
                acc.insertionSites.end());
            prof.insertionSites = std::move(acc.insertionSites);
        }

        profiles.push_back(std::move(prof));
    }

    // Deduplicate profiles by PHYSICAL ReadId (a read can appear under both
    // orientations); keep the one with the most aligned columns.
    {
        unordered_map<uint32_t, uint32_t> bestByReadId;
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

    // Materialize the full backbone base vector transiently for the emit tail
    // (it indexes absolute backbone positions, e.g. bbSeqVec[pos-2]). Built here
    // and freed on return, so at most `threadCount` of these exist at once --
    // never one per window across all windows.
    vector<uint8_t> bbSeqVec(plan.bbLen);
    for (uint32_t i = 0; i < plan.bbLen; i++)
        bbSeqVec[i] = rds.getOrientedReadBase(plan.bbOid, i).value;

    // Leaf-snarl detection (verification only; not fed into window.hetBubbles
    // or anything else -- production het-site output above is unaffected).
    // Reuses the SAME per-member profiles the interval-abPOA engine already
    // built for emitHetBubblesFromProfiles: no new alignment work, and no
    // shared-graph memory risk, since each profile came from an independent
    // small per-interval abPOA MSA (see the design note on
    // findLeafSnarlsFromProfiles for why this is architecturally safe for
    // large/dense windows, unlike the one-shared-graph-per-window approach in
    // AssemblerWindowAbpoaGraph.cpp).
    if (std::getenv("DINARA_INTERVALPOA_LEAFSNARL_DEBUG") != nullptr) {
        const std::vector<std::uint8_t> bbWindowBases(
            bbSeqVec.begin() + plan.windowBbBegin, bbSeqVec.begin() + plan.windowBbEnd);
        if (const char* probeEnv = std::getenv("DINARA_INTERVALPOA_LEAFSNARL_PROBE")) {
            const std::uint32_t probePos = std::uint32_t(std::atol(probeEnv));
            std::uint32_t nCov = 0, nAligned = 0, nMatch = 0, nMismatch = 0, nDeleted = 0;
            for (const auto& prof : profiles) {
                if (probePos < prof.bbCovBegin || probePos >= prof.bbCovEnd) continue;
                nCov++;
                const auto* c = prof.colAt(probePos);
                if (c) {
                    nAligned++;
                    if (c->readBase == bbWindowBases[probePos - plan.windowBbBegin]) nMatch++;
                    else nMismatch++;
                } else if (prof.isDeleted(probePos)) {
                    nDeleted++;
                }
            }
            std::cout << "    intervalPoaLeafSnarlProbe window=" << window.windowId
                 << " pos=" << probePos << " nProfilesTotal=" << profiles.size()
                 << " covering=" << nCov << " aligned=" << nAligned
                 << " match=" << nMatch << " mismatch=" << nMismatch
                 << " deleted=" << nDeleted << " neitherAlignedNorDeleted="
                 << (nCov - nAligned - nDeleted) << std::endl;
        }
        const auto snarls = findLeafSnarlsFromProfiles(
            plan.bbOid, profiles, bbWindowBases, plan.windowBbBegin,
            /*minSupport=*/6, /*minVaf=*/0.12,
            /*dropHomopolymer=*/hetDropHomopolymer, /*dropRepeat=*/hetDropRepeat);
        for (const LeafSnarl& s : snarls) {
            const std::int64_t localStart = std::int64_t(s.start) - std::int64_t(plan.windowBbBegin);
            const std::int64_t ctxBegin = std::max<std::int64_t>(0, localStart - 10);
            const std::int64_t ctxEnd = std::min<std::int64_t>(
                std::int64_t(bbWindowBases.size()), localStart + 12);
            std::string ctx;
            for (std::int64_t p = ctxBegin; p < ctxEnd; p++) {
                const std::uint8_t b = bbWindowBases[size_t(p)];
                ctx += "ACGTN"[b < 4 ? b : 4];
            }
            std::cout << "    intervalPoaLeafSnarl window=" << window.windowId
                 << " bb=" << plan.bbOid
                 << " [" << s.start << "," << s.end << ")"
                 << " ctx=" << ctx
                 << " alleles=" << s.alleles.size();
            for (const LeafSnarlAllele& al : s.alleles) {
                std::cout << " {n=" << al.members.size() << " bases=";
                for (std::uint8_t b : al.bases) std::cout << "ACGTN"[b < 4 ? b : 4];
                std::cout << "}";
            }
            std::cout << std::endl;
        }
    }

    return emitHetBubblesFromProfiles(
        window, profiles, bbSeqVec, plan.windowBbBegin, plan.windowBbEnd,
        plan.bbOid, plan.k, hetMinVaf, hetMinSupport,
        hetDropHomopolymer, hetDropRepeat, coverageHet, "intervalpoa");
}

}  // namespace dinara

#endif
