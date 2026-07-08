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
#include <chrono>
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

    const uint64_t coverageHet = assemblerInfo.isOpen ?
        assemblerInfo->kmerDistributionInfo.coverageHet : invalid<uint64_t>;

    // Per-position support pre-filter threshold: a het needs a ref allele and a
    // competing allele, each with >=minSupport members over disjoint sets, so a
    // passing position needs >= 2*minSupport-1 covering members. Every span
    // member covers the whole interval, so an interval with fewer members can
    // never yield a het and its MSA is skipped (see runIpoaOnRows).
    const int minSupport = resolveHetMinSupport(hetMinSupport, coverageHet);
    const int minMembers = 2 * minSupport - 1;

    // Stream each interval into the accumulator, reusing one scratch fragment so
    // per-interval outputs never accumulate.
    const uint32_t nIntervals = plan.intervalCount();
    IpoaWindowAccum wa;
    IpoaAbHandle ah;
    IpoaFragment scratch;
    for (uint32_t bi = 0; bi < nIntervals; bi++) {
        runIpoaInterval(plan, bi, rds, oneSidedEnabled, ah, scratch, minMembers);
        accumulateIpoaFragmentUnlocked(wa, scratch);
    }

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

    // Per-position support pre-filter threshold: a passing het position needs a
    // ref allele and a competing allele, each with >=minSupport members over
    // disjoint sets, hence >= 2*minSupport-1 covering members. Span members all
    // cover the whole interval, so an interval below this can never yield a het
    // and its MSA is skipped (see runIpoaOnRows).
    const int minSupport = resolveHetMinSupport(hetMinSupport, coverageHet);
    const int minMembers = 2 * minSupport - 1;

    const bool timing = ipoaTimingEnabled();
    if (timing) ipoaTiming().reset();
    const auto tAll0 = std::chrono::steady_clock::now();

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

                const auto tPlan0 = timing ?
                    std::chrono::steady_clock::now() :
                    std::chrono::steady_clock::time_point{};
                IpoaPlan plan =
                    buildIpoaPlan(windows[w], anchors, journeys, mkrs, rds, k);
                if (timing)
                    ipoaTiming().planNs.fetch_add(ipoaNsSince(tPlan0),
                                                  std::memory_order_relaxed);
                if (!plan.valid) continue;

                // All intervals of this window run serially into its own
                // accumulator -- single owner, so no locking is needed.
                IpoaWindowAccum wa;
                const uint32_t nI = plan.intervalCount();
                for (uint32_t bi = 0; bi < nI; bi++) {
                    runIpoaInterval(plan, bi, rds, oneSidedEnabled, ah, scratch,
                                    minMembers);
                    const auto tAcc0 = timing ?
                        std::chrono::steady_clock::now() :
                        std::chrono::steady_clock::time_point{};
                    accumulateIpoaFragmentUnlocked(wa, scratch);
                    if (timing)
                        ipoaTiming().accumNs.fetch_add(ipoaNsSince(tAcc0),
                                                       std::memory_order_relaxed);
                }

                const auto tMerge0 = timing ?
                    std::chrono::steady_clock::now() :
                    std::chrono::steady_clock::time_point{};
                const uint32_t n = mergeAndEmitIpoaWindow(
                    windows[w], plan, wa, rds, coverageHet,
                    hetMinVaf, hetMinSupport, hetDropHomopolymer, hetDropRepeat);
                if (timing)
                    ipoaTiming().mergeNs.fetch_add(ipoaNsSince(tMerge0),
                                                   std::memory_order_relaxed);
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

    if (timing) {
        const double wallSecs = std::chrono::duration_cast<
            std::chrono::duration<double>>(
                std::chrono::steady_clock::now() - tAll0).count();
        const IpoaTiming& t = ipoaTiming();
        // The phase counters are summed across threads (thread-busy time), so
        // they add up to ~threadCount * wall. ms() reports thread-seconds; the
        // percentage is share of total thread-busy time, which is what tells you
        // where the CPU actually went.
        const double planS    = double(t.planNs.load())    / 1e9;
        const double setupS   = double(t.setupNs.load())   / 1e9;
        const double resetS   = double(t.resetNs.load())   / 1e9;
        const double msaS     = double(t.msaNs.load())     / 1e9;
        const double extractS = double(t.extractNs.load()) / 1e9;
        const double accumS   = double(t.accumNs.load())   / 1e9;
        const double mergeS   = double(t.mergeNs.load())   / 1e9;
        const double busyS = planS + setupS + resetS + msaS + extractS +
                             accumS + mergeS;
        const double denom = busyS > 0 ? busyS : 1.0;
        const uint64_t nInt = t.intervals.load();
        const uint64_t nRun = t.msaRuns.load();
        const double nRunD = nRun ? double(nRun) : 1.0;
        auto pct = [&](double s) { return 100.0 * s / denom; };
        std::cout << std::fixed << std::setprecision(2)
                  << "[HetTiming] wall=" << wallSecs << "s"
                  << " threadBusy=" << busyS << "s"
                  << " intervals=" << nInt
                  << " msaRuns=" << nRun
                  << " skipped=" << t.skippedIntervals.load()
                  << " skippedByLen=" << t.skippedByLen.load() << "\n"
                  << "[HetTiming]   buildPlan   " << planS    << "s (" << pct(planS)    << "%)\n"
                  << "[HetTiming]   setup       " << setupS   << "s (" << pct(setupS)   << "%)  member-gather + code arrays\n"
                  << "[HetTiming]   abpoa_reset " << resetS   << "s (" << pct(resetS)   << "%)  [" << (1e6 * resetS / nRunD) << " us/run]\n"
                  << "[HetTiming]   abpoa_msa   " << msaS     << "s (" << pct(msaS)     << "%)  [" << (1e6 * msaS / nRunD) << " us/run]\n"
                  << "[HetTiming]   extract     " << extractS << "s (" << pct(extractS) << "%)  msa copy + column walk\n"
                  << "[HetTiming]   accumulate  " << accumS   << "s (" << pct(accumS)   << "%)  fold fragment\n"
                  << "[HetTiming]   mergeEmit   " << mergeS   << "s (" << pct(mergeS)   << "%)  merge + emit bubbles\n"
                  << "[HetTiming]   msa dims: nSeq avg=" << (double(t.seqSum.load()) / nRunD) << " max=" << t.seqMax.load()
                  << "  maxLen(bp) avg=" << (double(t.lenSum.load()) / nRunD) << " max=" << t.lenMax.load()
                  << "  msaCols avg=" << (double(t.colSum.load()) / nRunD) << " max=" << t.colMax.load() << "\n"
                  << "[HetTiming]   run dist: slowest=" << (double(t.runMaxNs.load()) / 1e6) << "ms"
                  << "  >100us=" << t.runOver100us.load()
                  << "  >1ms=" << t.runOver1ms.load()
                  << "  >10ms=" << t.runOver10ms.load()
                  << "  time in >1ms runs=" << (double(t.nsOver1ms.load()) / 1e9) << "s ("
                  << (100.0 * double(t.nsOver1ms.load()) / 1e9 / ((msaS + resetS) > 0 ? (msaS + resetS) : 1.0)) << "% of reset+msa)\n"
                  << std::defaultfloat << std::flush;
    }

    hetWindowsOut = hetWindows.load();
    totalBubblesOut = totalBubbles.load();
    return uint32_t(totalBubbles.load());
}
