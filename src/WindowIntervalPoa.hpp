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
        intervals = 0; msaRuns = 0;
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
    const IpoaClock::time_point& tSetup0);

// Shared MSA + emit for a backbone interval. Given `spanMembers` already placed
// over [bbBegin, bbEnd), run ONE abPOA MSA with row 0 = backbone segment and
// emit each member's per-column contribution into `frag`. Split out from
// runIpoaInterval so the member-selection and the MSA/extraction concerns are
// separable. Requires >=2 two-sided members. Safe to call from any thread.
inline void runIpoaOnRows(
    const IpoaPlan& plan,
    std::uint32_t bbBegin,
    std::uint32_t bbEnd,
    std::vector<IpoaSpanMember>& spanMembers,
    const Reads& rds,
    IpoaAbHandle& ah,
    IpoaFragment& frag,
    bool timing,
    const IpoaClock::time_point& tSetup0)
{
    using std::vector;
    using std::uint8_t;
    using std::uint32_t;
    using std::uint64_t;
    using std::int64_t;

    if (spanMembers.empty()) return;

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
                    entered = true;
                    readAbs++;
                } else {
                    if (entered && pendingDelBegin < 0)
                        pendingDelBegin = bbPos;
                }
            } else {
                if (!mb.isGap()) readAbs++;
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
    IpoaFragment& frag)
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
                  timing, tSetup0);
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

    return emitHetBubblesFromProfiles(
        window, profiles, bbSeqVec, plan.windowBbBegin, plan.windowBbEnd,
        plan.bbOid, plan.k, hetMinVaf, hetMinSupport,
        hetDropHomopolymer, hetDropRepeat, coverageHet, "intervalpoa");
}

}  // namespace dinara

#endif
