/// @file AssemblerWindowIntervalPoaHetSites.cpp
/// @brief Per-window het-SNP detection using one small POA per consecutive-anchor
///        interval (Shasta2 LocalAssembly model), instead of a single
///        whole-window POA.
///
/// Motivation
/// ----------
/// A whole-window abPOA builds ONE graph per window and folds every read into it.
/// That graph grows as reads are added, so a read's segment DP cost scales with
/// the accumulated graph span (not local depth), and the whole window runs on a
/// single thread. On a 100kb / 200-read window this dominates assembler time.
///
/// Shasta2's LocalAssembly avoids this: it never builds a whole-window graph. It
/// assembles between consecutive anchors, so each POA sees only the reads
/// spanning one short gap and a tiny graph. This engine ports that idea to het
/// detection:
///
///   1. Tile the window at its backbone anchors -> intervals [a_i, a_{i+1}).
///   2. For each interval, gather the backbone segment (row 0) + each spanning
///      read's bracketed subsequence, and run one abPOA MSA.
///   3. Map each interval's MSA columns to ABSOLUTE backbone offsets and merge
///      into per-read KwMemberProfile rows.
///   4. The shared emitHetBubblesFromProfiles tail emits AnchorWindow::hetBubbles.
///
/// Parallelism (shasta2 assembleChainsMultithreaded model)
/// -------------------------------------------------------
/// The per-interval MSA is the unit of work. Rather than parallelize over
/// windows (which leaves threads idle when a few windows hold most intervals),
/// intervalPoaDetectHetBubblesAllWindows flattens EVERY interval from EVERY
/// window into one global list, sorts by descending cost, and lets threads pull
/// steps one at a time (batch=1) -- exactly like shasta2 flattens every
/// consecutive-anchor AssemblyStep across the whole assembly graph. Each step
/// writes its own IpoaFragment (no shared state), then a serial per-window pass
/// merges fragments and emits bubbles.
///
/// The heavy lifting (plan build, interval MSA, merge+emit) lives in
/// WindowIntervalPoa.hpp so both the single-window entry point and the
/// all-windows orchestrator share identical logic.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "AssemblerOptions.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "WindowIntervalPoa.hpp"
#include "invalid.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

using namespace std;
using namespace dinara;


namespace {

// Emit the DINARA_HET_DEBUG per-window instrumentation line from a window's
// streamed accumulator. Shared by both entry points.
void ipoaDebugLine(const IpoaPlan& plan, const IpoaWindowAccum& wa)
{
    const double oneSidedFrac = (wa.twoSidedTotal + wa.oneSidedTotal) > 0
        ? double(wa.oneSidedTotal) / double(wa.twoSidedTotal + wa.oneSidedTotal) : 0.0;
    cout << "    intervalPoaDetectHetBubbles bb=" << plan.bbOid
         << " window=[" << plan.windowBbBegin << "," << plan.windowBbEnd << ")"
         << " members=" << plan.memberCount
         << " intervals=" << wa.intervalsRun
         << " twoSided=" << wa.twoSidedTotal
         << " oneSided=" << wa.oneSidedTotal
         << " oneSidedPlaced=" << wa.oneSidedPlaced
         << " oneSidedFrac=" << std::fixed << std::setprecision(3) << oneSidedFrac
         << std::defaultfloat
         << " intervalsWithOneSided=" << wa.intervalsWithOneSided << endl;
}

}  // anonymous namespace


// Detect het bubbles in ONE window with per-interval POA. Kept for callers that
// want a single window; the all-windows orchestrator below is preferred for the
// main pipeline because it load-balances intervals across the whole run.
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

    window.hetBubbles.clear();

    const IpoaPlan plan = buildIpoaPlan(window, anchors, journeys, mkrs, rds, k);
    if (!plan.valid) return 0;

    const bool oneSidedEnabled = ipoaOneSidedEnabled();
    const bool debug = (getenv("DINARA_HET_DEBUG") != nullptr);

    // Stream each interval into the accumulator, reusing one scratch fragment so
    // per-interval outputs never accumulate.
    const uint32_t nIntervals = plan.intervalCount();
    IpoaWindowAccum wa;
    IpoaAbHandle ah;
    IpoaFragment scratch;
    for (uint32_t bi = 0; bi < nIntervals; bi++) {
        runIpoaInterval(plan, bi, rds, oneSidedEnabled, ah, scratch);
        accumulateIpoaFragmentUnlocked(wa, scratch);
    }

    const uint64_t coverageHet = assemblerInfo.isOpen ?
        assemblerInfo->kmerDistributionInfo.coverageHet : invalid<uint64_t>;

    const uint32_t n = mergeAndEmitIpoaWindow(
        window, plan, wa, rds, coverageHet,
        hetMinVaf, hetMinSupport, hetDropHomopolymer, hetDropRepeat);

    if (debug) ipoaDebugLine(plan, wa);
    return n;
}


// Detect het bubbles across ALL windows with global interval load balancing.
//
// Phase 1 (serial): build a plan per window and enumerate its intervals into a
//   flat step list, each step recording (window, interval, cost). Cost is a
//   cheap span*rows proxy so the biggest MSAs dispatch first (shasta2 sorts
//   steps by descending offsetInBases for the same reason).
// Phase 2 (parallel, batch=1): each thread pulls one step, runs its interval
//   MSA into that step's fragment slot. No locks on the hot path -- every slot
//   is written by exactly one thread.
// Phase 3 (parallel over windows): merge each window's fragments and emit
//   bubbles. Windows are independent so this parallelizes cleanly.
uint32_t Assembler::intervalPoaDetectHetBubblesAllWindows(
    vector<AnchorWindow>& windows,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    double hetMinVaf,
    uint64_t hetMinSupport,
    bool hetDropHomopolymer,
    bool hetDropRepeat,
    uint64_t threadCount,
    uint64_t& hetWindowsOut,
    uint64_t& totalBubblesOut) const
{
    const Reads& rds = getReads();
    const auto& mkrs = *markers;
    const uint64_t k = assemblerInfo->k;
    const bool oneSidedEnabled = ipoaOneSidedEnabled();
    const bool debug = (getenv("DINARA_HET_DEBUG") != nullptr);
    if (threadCount == 0) threadCount = std::thread::hardware_concurrency();

    const uint64_t coverageHet = assemblerInfo.isOpen ?
        assemblerInfo->kmerDistributionInfo.coverageHet : invalid<uint64_t>;

    std::atomic<uint64_t> hetWindows{0};
    std::atomic<uint64_t> totalBubbles{0};
    std::mutex debugMutex;

    // Clear all het-bubble slots up front (windows not processed / homozygous
    // must end up empty).
    for (AnchorWindow& w : windows) w.hetBubbles.clear();

    // --- Per-window work-stealing ------------------------------------------
    // One persistent thread pool; each thread pulls the next window and handles
    // it ENTIRELY on its own: build plan, run all its intervals serially into
    // its own accumulator, merge+emit, free. Full utilization whenever there are
    // more windows than threads (the normal case) with none of the overhead a
    // flatten-all-intervals scheme adds -- no per-chunk barriers, no accumulator
    // locking (one thread owns the window), one pool spawn for the whole run.
    //
    // Memory stays bounded: at most ~threadCount windows are resident (one per
    // worker), and each window streams its intervals so its peak is just the
    // accumulator (~depth*span) -- the same footprint as the original per-window
    // driver that worked.
    //
    // Windows are visited biggest-first (by backbone span) so a dominant window
    // starts early and does not end up as an end-of-run straggler running alone.
    vector<uint32_t> order(windows.size());
    for (uint32_t i = 0; i < order.size(); i++) order[i] = i;
    sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        return windows[a].baseSpan > windows[b].baseSpan;
    });

    {
        std::atomic<uint64_t> nextIdx{0};
        const uint64_t nWindows = windows.size();
        auto worker = [&]() {
            IpoaAbHandle ah;       // thread-local abPOA handle, reused
            IpoaFragment scratch;  // thread-local, reused across intervals
            for (;;) {
                const uint64_t idx = nextIdx.fetch_add(1);
                if (idx >= nWindows) break;
                const uint32_t w = order[idx];

                IpoaPlan plan =
                    buildIpoaPlan(windows[w], anchors, journeys, mkrs, rds, k);
                if (!plan.valid) continue;

                // All intervals of this window run serially into its own
                // accumulator -- single owner, so no locking is needed.
                IpoaWindowAccum wa;
                const uint32_t nI = plan.intervalCount();
                for (uint32_t bi = 0; bi < nI; bi++) {
                    runIpoaInterval(plan, bi, rds, oneSidedEnabled, ah, scratch);
                    accumulateIpoaFragmentUnlocked(wa, scratch);
                }

                const uint32_t n = mergeAndEmitIpoaWindow(
                    windows[w], plan, wa, rds, coverageHet,
                    hetMinVaf, hetMinSupport, hetDropHomopolymer, hetDropRepeat);
                if (n > 0) { hetWindows.fetch_add(1); totalBubbles.fetch_add(n); }
                if (debug) {
                    std::lock_guard<std::mutex> lock(debugMutex);
                    ipoaDebugLine(plan, wa);
                }
            }
        };
        vector<std::thread> pool;
        pool.reserve(threadCount);
        for (uint64_t t = 0; t < threadCount; t++) pool.emplace_back(worker);
        for (auto& th : pool) th.join();
    }

    hetWindowsOut = hetWindows.load();
    totalBubblesOut = totalBubbles.load();
    return uint32_t(totalBubbles.load());
}
