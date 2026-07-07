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
        accumulateIpoaFragment(wa, scratch);
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

    // --- Chunked processing to bound memory --------------------------------
    // Building EVERY window's plan and holding EVERY window's fragments at once
    // is O(genome): each window's aligned-column fragments are ~depth*span, and
    // summed over thousands of windows that is many GB -> OOM. Process windows
    // in small chunks so only a chunk's plans+fragments are resident at a time.
    //
    // The chunk size is a WINDOW COUNT, defaulting to threadCount. This makes
    // resident memory match the OLD per-window driver exactly (it too kept
    // ~threadCount windows in flight, one per worker) -- so this path can never
    // use more het memory than the version that worked. Within a chunk,
    // intervals are still flattened and work-stolen across all threads, so when
    // a few big windows dominate their thousands of intervals fan out over every
    // thread -- the load-balancing win, at old-path memory. Override with
    // DINARA_IPOA_CHUNK_WINDOWS to trade memory for cross-window overlap.
    uint32_t chunkWindows = uint32_t(threadCount);
    if (const char* e = getenv("DINARA_IPOA_CHUNK_WINDOWS")) {
        const long v = atol(e);
        if (v > 0) chunkWindows = uint32_t(v);
    }
    if (chunkWindows == 0) chunkWindows = 1;

    struct Step {
        uint32_t windowIndex;    // absolute window index
        uint32_t localIndex;     // index within the current chunk
        uint32_t intervalIndex;
        uint64_t cost;           // descending-cost load balancing
        bool operator<(const Step& that) const { return cost > that.cost; }
    };

    uint32_t wBegin = 0;
    while (wBegin < windows.size()) {
        const uint32_t wEnd =
            std::min<uint32_t>(uint32_t(windows.size()), wBegin + chunkWindows);
        const uint32_t chunkN = wEnd - wBegin;

        // Build this chunk's plans and flatten its interval steps. Each window
        // has one streaming accumulator; interval outputs are folded in as they
        // complete and never stored, so peak memory per window is just its
        // accumulator (~depth*span) -- matching the original per-window path.
        vector<IpoaPlan> plans(chunkN);
        vector<IpoaWindowAccum> accums(chunkN);
        vector<std::atomic<uint32_t>> remaining(chunkN);
        vector<Step> steps;
        for (uint32_t li = 0; li < chunkN; li++) {
            const uint32_t w = wBegin + li;
            plans[li] = buildIpoaPlan(windows[w], anchors, journeys, mkrs, rds, k);
            const uint32_t nI = plans[li].valid ? plans[li].intervalCount() : 0u;
            remaining[li].store(nI);
            if (nI == 0) continue;
            const IpoaPlan& plan = plans[li];
            for (uint32_t bi = 0; bi < nI; bi++) {
                const uint32_t span = plan.breakpoints[bi + 1] - plan.breakpoints[bi];
                const uint64_t rows =
                    plan.atBreakpoint[bi].size() + plan.atBreakpoint[bi + 1].size();
                steps.push_back({w, li, bi, uint64_t(span) * (rows + 1)});
            }
        }

        // Biggest MSAs first so a straggler never runs alone at the end.
        sort(steps.begin(), steps.end());

        // Merge+emit a chunk window once its last interval finishes, then free
        // its heavy state immediately.
        auto finishWindow = [&](uint32_t li) {
            const uint32_t w = wBegin + li;
            const uint32_t n = mergeAndEmitIpoaWindow(
                windows[w], plans[li], accums[li], rds, coverageHet,
                hetMinVaf, hetMinSupport, hetDropHomopolymer, hetDropRepeat);
            if (n > 0) { hetWindows.fetch_add(1); totalBubbles.fetch_add(n); }
            if (debug) {
                std::lock_guard<std::mutex> lock(debugMutex);
                ipoaDebugLine(plans[li], accums[li]);
            }
            plans[li].freeHeavy();
        };

        // Run this chunk's interval MSAs (work-stealing, batch=1). Each thread
        // reuses one scratch fragment: run interval -> fold into the window
        // accumulator -> clear. No per-interval storage survives.
        {
            std::atomic<uint64_t> nextStep{0};
            const uint64_t nSteps = steps.size();
            auto worker = [&]() {
                IpoaAbHandle ah;       // thread-local abPOA handle, reused
                IpoaFragment scratch;  // thread-local, reused across intervals
                for (;;) {
                    const uint64_t i = nextStep.fetch_add(1);
                    if (i >= nSteps) break;
                    const Step& s = steps[i];
                    runIpoaInterval(plans[s.localIndex], s.intervalIndex,
                                    rds, oneSidedEnabled, ah, scratch);
                    accumulateIpoaFragment(accums[s.localIndex], scratch);
                    if (remaining[s.localIndex].fetch_sub(1) == 1)
                        finishWindow(s.localIndex);
                }
            };
            vector<std::thread> pool;
            pool.reserve(threadCount);
            for (uint64_t t = 0; t < threadCount; t++) pool.emplace_back(worker);
            for (auto& th : pool) th.join();
        }

        wBegin = wEnd;
    }

    hetWindowsOut = hetWindows.load();
    totalBubblesOut = totalBubbles.load();
    return uint32_t(totalBubbles.load());
}
