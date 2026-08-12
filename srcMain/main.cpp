// The static executable provides
// basic functionality and reduced performance.
// For full functionality use the shared library built
// under directory src.

// Dinara.
#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "InvertedIndexBuilder.hpp"
#include "AssemblerOptions.hpp"
#include "WindowHetProfiles.hpp"
#include "buildId.hpp"
#if DINARA_ENABLE_VARIANT_CLUSTERING
#include "ClusterGraph.hpp"
#endif
#include "filesystem.hpp"
#include "mode3-Anchor.hpp"
#include "mode3-BidirectedAnchor.hpp"
#include "mode3-DirectedAnchors.hpp"
#include "mode3-DirectedAnchorGraph.hpp"
#include "mode3-AnchorGraph.hpp"
#include "mode3-AnchorGraphSuperbubbles.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2AnchorPair.hpp"
#include "Shasta2AnchorsFromSplitVertices.hpp"
#include "Shasta2Journeys.hpp"
#include "Shasta2AnchorGraph.hpp"
#include "Shasta2AssemblyGraph.hpp"
#include "BidirectedAnchorGraph.hpp"
#include "DinaraDetangle.hpp"
#include "WindowTransitions.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "Tee.hpp"
#include "timestamp.hpp"
#include "platformDependent.hpp"

// hifiasm candidate-overlap detector (submodule). C API used by the
// --overlapsFromHifiasm path to generate read overlaps as PAF.
#include "hifiasm_overlaps.h"


using namespace dinara;

// Boost libraries.
#include <boost/program_options.hpp>
#include  <boost/chrono/process_cpu_clocks.hpp>

//  Linux.
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

// Standard library.
#include "chrono.hpp"
#include <algorithm>
#include <numeric>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include "iostream.hpp"
#include <set>
#include <functional>
#include <unordered_map>
#include <unordered_set>

#include "stdexcept.hpp"


// Shasta 2 Integration
#include "AssemblerShasta2Anchors.hpp"
#include "HetAnchorK.hpp"
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

// Minimizer filtering
#include "MinimizerChecker.hpp"
#include "MarkerKmers.hpp"
#include "KmerCounter.hpp"



namespace dinara {
    namespace main {

        void main(int argumentCount, const char** arguments);

        void setupRunDirectory(
            const string& memoryMode,
            const string& memoryBacking,
            size_t& pageSize,
            string& dataDirectory
            );

        void setupHugePages();
        void segmentFaultHandler(int);
        void setupSegmentFaultHandler();

        // Functions that implement --command keywords
        void assemble(const AssemblerOptions&, int argumentCount, const char** arguments);
        void saveBinaryData(const AssemblerOptions&);
        void cleanupBinaryData(const AssemblerOptions&);
        void explore(const AssemblerOptions&);
        void listCommands();

        const std::set<string> commands = {
            "assemble",
            "saveBinaryData",
            "cleanupBinaryData",
            "explore",
            "listCommands"};



        void assemble(
            Assembler&,
            const AssemblerOptions&,
            vector<string> inputNames);

    }

    // This is used to duplicate cout output to stdout.log.
    Tee tee;
    ofstream dinaraLog;
}


namespace {

using BubbleEndpointPair = std::pair<std::string, std::string>;

BubbleEndpointPair canonicalBubblePair(std::string a, std::string b)
{
    if(a <= b) {
        return {std::move(a), std::move(b)};
    }
    return {std::move(b), std::move(a)};
}


std::string normalizeBubbleEndpoint(std::string s)
{
    if(!s.empty() && (s.back() == '+' || s.back() == '-')) {
        s.pop_back();
    }
    return s;
}


std::string shellQuote(const std::string& s)
{
    std::string quoted = "'";
    for(const char c: s) {
        if(c == '\'') {
            quoted += "'\"'\"'";
        } else {
            quoted += c;
        }
    }
    quoted += "'";
    return quoted;
}


bool envFlagIsEnabled(const char* variableName)
{
    const char* value = ::getenv(variableName);
    if(value == nullptr) {
        return false;
    }
    std::string s(value);
    std::transform(
        s.begin(),
        s.end(),
        s.begin(),
        [](unsigned char c) { return char(std::tolower(c)); });
    return
        s == "1" ||
        s == "true" ||
        s == "yes" ||
        s == "on";
}


bool envFlagIsDisabled(const char* variableName)
{
    const char* value = ::getenv(variableName);
    if(value == nullptr) {
        return false;
    }
    std::string s(value);
    std::transform(
        s.begin(),
        s.end(),
        s.begin(),
        [](unsigned char c) { return char(std::tolower(c)); });
    return
        s == "0" ||
        s == "false" ||
        s == "no" ||
        s == "off";
}


uint64_t envUintOrDefault(const char* variableName, uint64_t defaultValue)
{
    const char* value = ::getenv(variableName);
    if(value == nullptr || *value == '\0') {
        return defaultValue;
    }
    try {
        const uint64_t parsed = std::stoull(value);
        return std::max<uint64_t>(1, parsed);
    } catch(...) {
        return defaultValue;
    }
}


// Run `body(windowIndex)` for every window across `threadCount` worker threads
// (work-stealing over a shared atomic counter, batch=1). The post-MSA het passes
// (plan, merge, stage, verify) are all window-local: they read shared read-only
// state (the anchor store, journeys, the frozen primary set) and mutate only
// their own window, so they parallelize cleanly over windows the same way the
// MSA phase already does. Only the anchor APPEND pass stays serial (it grows the
// memory-mapped store and assigns ids in order); it is not run through here.
//
// `body` must be safe to call concurrently for distinct window indices: use
// thread-local scratch (bbAnchors/bbOffset) and atomic counters. Any exception
// thrown by `body` (the verifiers throw on a structural violation) is captured
// and rethrown on the calling thread after all workers join, so a violation
// still aborts the run with its original message.
template<class Body>
void parallelForEachWindow(uint64_t windowCount, uint64_t threadCount, Body&& body)
{
    if(threadCount == 0) threadCount = std::thread::hardware_concurrency();
    if(threadCount == 0) threadCount = 1;
    threadCount = std::min<uint64_t>(threadCount, std::max<uint64_t>(1, windowCount));

    std::atomic<uint64_t> nextIdx{0};
    std::mutex errMutex;
    std::exception_ptr firstError;

    auto worker = [&]() {
        for(;;) {
            const uint64_t w = nextIdx.fetch_add(1);
            if(w >= windowCount) break;
            try {
                body(w);
            } catch(...) {
                std::lock_guard<std::mutex> lock(errMutex);
                if(!firstError) firstError = std::current_exception();
                // Fast-drain the rest so we stop promptly after a failure.
                nextIdx.store(windowCount, std::memory_order_relaxed);
                break;
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(threadCount);
    for(uint64_t t = 0; t < threadCount; t++) pool.emplace_back(worker);
    for(auto& th : pool) th.join();

    if(firstError) std::rethrow_exception(firstError);
}


// ============================================================================
// Intra-window het-bubble wiring helpers.
//
// After per-window MSA het detection stages het bubbles on each window, three
// serial passes turn those bubbles into anchor-graph structure:
//   1. planWindowHetBubbles   — assign each bubble to the backbone interval that
//                               strictly contains its flank span; drop the rest.
//   2. appendWindowHetAnchors — create the anchor pairs for planned bubbles.
//   3. stageWindowIntraEdges  — stage the backbone chain + het bubble edges.
// All three need the window's ordered backbone anchors and their base offsets,
// computed once by computeWindowBackbone.
// ============================================================================

// Ordered backbone anchors of a window plus their base offsets from the window
// start. Returns false (empty vectors) when the window has fewer than two
// backbone anchors, i.e. no interval to contain a bubble.
bool computeWindowBackbone(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const AnchorWindow& window,
    uint32_t kHalf,
    vector<Shasta2AnchorId>& bbAnchors,
    vector<uint32_t>& bbOffset)
{
    bbAnchors.clear();
    bbOffset.clear();

    const OrientedReadId backboneOid = window.backboneOrientedReadId;
    const auto backboneJourney = journeys[backboneOid];

    const auto& positions = window.filteredBackbonePositions;
    if(!positions.empty()) {
        for(const uint32_t pos : positions)
            bbAnchors.push_back(backboneJourney[pos]);
    } else {
        for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++)
            bbAnchors.push_back(backboneJourney[pos]);
    }
    if(bbAnchors.size() < 2) {
        bbAnchors.clear();
        return false;
    }

    const uint32_t backboneBeginPos =
        anchors.getPosition(bbAnchors.front(), backboneOid) - kHalf;
    bbOffset.resize(bbAnchors.size());
    for(size_t i = 0; i < bbAnchors.size(); i++)
        bbOffset[i] = anchors.getPosition(bbAnchors[i], backboneOid) - backboneBeginPos;
    return true;
}


// Per-window backbone (bbAnchors/bbOffset), computed ONCE and shared by all
// three het-anchor wiring passes below (plan, merge, stage). Nothing between
// window creation and the stage pass mutates a window's own backbone fields
// (filteredBackbonePositions/backboneBegin/backboneEnd) or moves an existing
// anchor's position -- the intervening het-anchor append only ADDS new
// anchors -- so the three passes would otherwise recompute byte-identical
// results. computeWindowBackbone's cost is O(backbone anchors x anchor
// coverage) per window (getPosition linearly scans an anchor's markers to
// find one read), so computing it once instead of three times removes 2/3 of
// that work across every window in the assembly.
class WindowBackboneCache {
public:
    WindowBackboneCache(
        const Shasta2Anchors& anchors,
        const Shasta2Journeys& journeys,
        const vector<AnchorWindow>& windows,
        uint32_t kHalf,
        uint64_t threadCount) :
        bbAnchors_(windows.size()),
        bbOffset_(windows.size())
    {
        parallelForEachWindow(windows.size(), threadCount, [&](uint64_t wi) {
            computeWindowBackbone(anchors, journeys, windows[wi], kHalf,
                bbAnchors_[wi], bbOffset_[wi]);
        });
    }

    // False when the window has fewer than two backbone anchors (no interval
    // to contain a bubble); bbAnchors(wi)/bbOffset(wi) are empty in that case.
    bool hasBackbone(uint64_t wi) const { return !bbAnchors_[wi].empty(); }
    const vector<Shasta2AnchorId>& bbAnchors(uint64_t wi) const { return bbAnchors_[wi]; }
    const vector<uint32_t>& bbOffset(uint64_t wi) const { return bbOffset_[wi]; }

private:
    vector<vector<Shasta2AnchorId>> bbAnchors_;
    vector<vector<uint32_t>> bbOffset_;
};


// Step 1 of the per-window het pipeline: the explicit "middle-2 backbone
// shift". Every k=50 backbone anchor is exported to shasta2 as the CENTERED
// 2-base clip of its 50-mer footprint: the exported raw position is
// (storedMidpoint - 1), i.e. bases [midpoint-1, midpoint]. Because the full
// 50-mer agrees across all member reads, the centered 2 bases agree too, so the
// clip is a valid k=2 anchor while keeping guaranteed-good flanking sequence on
// both sides. This is the SAME frame het detection runs against (both use the
// anchor-midpoint origin; see computeWindowBackbone and the abPOA backboneCodes
// span), and the same frame the monotonicity verifier checks. Making it an
// explicit, named step keeps the three coordinate frames (midpoint, exported,
// MSA) reconciled in one place.
//
// bbOffset[i] is anchor i's midpoint offset from the window start; the exported
// clip sits k/2 bases earlier, so the exported window-local offset is
// bbOffset[i] - k/2 (1 for k=2, 0 for k=0). The export shift is UNIFORM across
// anchor classes (shasta2 uses a single --k), so the backbone uses the same
// hetAnchorKHalf() the het anchors do. Returns those exported offsets.
// (bbOffset[0] >= kHalf by construction, so the subtraction never underflows.)
void computeWindowShiftedBackbone(
    const vector<uint32_t>& bbOffset,
    vector<uint32_t>& bbExportedOffset)
{
    const uint32_t hetKHalf = hetAnchorKHalf();
    bbExportedOffset.resize(bbOffset.size());
    for(size_t i = 0; i < bbOffset.size(); i++) {
        DINARA_ASSERT(bbOffset[i] >= hetKHalf);   // midpoint frame: >= kHalf
        bbExportedOffset[i] = bbOffset[i] - hetKHalf;
    }
}


// Set of (oriented-read, exported/shasta2 position) markers owned by primary
// anchors. Built once over the primary anchors BEFORE any het/hom append. Every
// anchor is exported at (storedMidpoint - k/2) with the SAME uniform export
// shift and shasta2 re-adds that same k/2 on load, so an anchor's shasta2 marker
// position equals its dinara stored midpoint regardless of k. Thus two anchors
// occupy the same shasta2 (read, position) marker iff their stored midpoints
// coincide on that read -- a k-independent test, so this set is compared in the
// stored-midpoint frame directly. Used during planning to drop any het bubble
// whose bracketing homs would duplicate an existing primary anchor's markers
// (shasta2 read-following would otherwise pair the two ids into a zero-length
// assembly step and assert).
using PrimaryMarkerSet = std::unordered_set<uint64_t>;

static inline uint64_t primaryMarkerKey(OrientedReadId oid, uint32_t position)
{
    return (uint64_t(oid.getValue()) << 32) | uint64_t(position);
}

PrimaryMarkerSet buildPrimaryMarkerSet(const Shasta2Anchors& anchors)
{
    // Called before any het/hom append, so anchors.size() is the primary count
    // and every stored marker belongs to a primary anchor.
    PrimaryMarkerSet set;
    // Without this, the set grows one insert at a time with no reserve --
    // for ~21M markers that means dozens of full rehashes (each touching every
    // element inserted so far), which profiling showed costing ~11s on a small
    // test region. anchorMarkerInfos.totalSize() is the exact final element
    // count (sum of every anchor's member count), available up front with no
    // extra pass, so reserve once and insert without ever rehashing.
    set.reserve(anchors.anchorMarkerInfos.totalSize());
    const uint64_t primaryCount = anchors.size();
    for(Shasta2AnchorId aid = 0; aid < primaryCount; aid++) {
        const Shasta2Anchor anchor = anchors[aid];
        for(const Shasta2AnchorMarkerInfo& mi : anchor)
            set.insert(primaryMarkerKey(mi.orientedReadId, mi.position));
    }
    return set;
}

// True if any member of this het/hom anchor lands on an existing primary anchor's
// (read, position) marker. A het/hom member's stored midpoint is rawPosition +
// hetAnchorKHalf() (1 for k=2, 0 for k=0), matching appendHetAnchorPair; the
// primary set is keyed in the same stored-midpoint frame (see buildPrimaryMarkerSet).
bool hetAnchorCollidesWithPrimary(
    const AnchorWindow::HetAnchor& a, const PrimaryMarkerSet& primarySet)
{
    // rawPosition + k/2 is the anchor's stored midpoint = the SNP read position
    // in both k=2 (predReadPos+1) and k=0 (snpReadPos+0). See HetAnchorK.hpp.
    const uint32_t hetKHalf = hetAnchorKHalf();
    for(const auto& m : a.members)
        if(primarySet.count(primaryMarkerKey(m.orientedReadId, m.rawPosition + hetKHalf)))
            return true;
    return false;
}


// Assign each het bubble in a window to the backbone interval [bbA_i, bbA_{i+1})
// that STRICTLY contains its full flank span [predBackboneOffset,
// succBackboneOffset] (leading hom .. trailing hom). A bubble whose flank span
// crosses/touches a backbone anchor, lacks either bracketing hom, or lacks any
// arm is contained by NO interval and is dropped (plannedInterval = -1).
// Accumulates counts into plannedBubbles/dropped.
//
// Because every bubble is self-bracketed by its own leading and trailing hom,
// there is no last-in-interval special case: consecutive bubbles in the same
// interval are chained hom-to-hom and every backbone-touching edge is a
// hom<->backbone edge with a non-empty read intersection.
// Drop the members of a bracketing hom whose EXPORTED position would create a
// backward edge against the flanking backbone anchor on that same read.
//
// Only the homs touch backbone anchors (arms connect solely to homs). A hom is
// exported at rawPosition + hetKHalf; a backbone anchor at markerPosition +
// k/2. The verifier compares these stored positions directly, requiring every
// shared read to step strictly forward:
//   bbA_i    -> leadHom : leadHom.pos  >  getPosition(bbA_i,   read)
//   hom      -> bbA_i+1 : getPosition(bbA_i+1, read)  >  hom.pos
//
// A member's rawPosition can drift (one-sided reads walk a guessed span through
// indels), landing on or past the read's TRUE position at a flanking backbone
// anchor -- and that anchor may not even be one of the read's POA breakpoints,
// so the interval walk never bounded against it. getPosition consults the read's
// actual marker position on the backbone anchor, so it catches the drift
// regardless. Drop only the offending members; the hom keeps the rest (it stays
// valid as long as it retains members, checked by the caller).
void dropBackwardHomMembers(
    const Shasta2Anchors& anchors,
    AnchorWindow::HetAnchor& hom,
    uint32_t hetKHalf,
    Shasta2AnchorId flankingBackbone,
    bool homIsLeading,             // true: bbA -> hom; false: hom -> bbA
    uint64_t& droppedMembers)
{
    vector<AnchorWindow::HetAnchorMember> kept;
    kept.reserve(hom.members.size());
    for(const auto& m : hom.members) {
        const uint32_t bbPos = anchors.getPosition(flankingBackbone, m.orientedReadId);
        if(bbPos == invalid<uint32_t>) {
            // Read is not on this backbone anchor: no edge on this read, keep it.
            kept.push_back(m);
            continue;
        }
        const uint32_t homPos = m.rawPosition + hetKHalf;
        const bool forward = homIsLeading ? (homPos > bbPos)   // bbA -> leadHom
                                          : (bbPos > homPos);   // hom  -> bbA
        if(forward) kept.push_back(m);
        else ++droppedMembers;
    }
    hom.members.swap(kept);
}

void planWindowHetBubbles(
    AnchorWindow& window,
    const Shasta2Anchors& anchors,
    const vector<Shasta2AnchorId>& bbAnchors,
    const vector<uint32_t>& bbOffset,
    uint32_t hetKHalf,
    const PrimaryMarkerSet& primarySet,
    uint64_t& plannedBubbles,
    uint64_t& dropped,
    uint64_t& droppedPrimaryCollision,
    uint64_t& droppedBackwardMembers)
{
    for(auto& b : window.hetBubbles) {
        b.plannedInterval = -1;
        const bool haveLead = !b.leadHom.members.empty();
        const bool haveTail = !b.hom.members.empty();
        bool hasArm = false;
        for(const auto& a : b.alleles) if(!a.members.empty()) { hasArm = true; break; }
        if(!(haveLead && haveTail && hasArm)) { ++dropped; continue; }

        // The span runs from the leading hom to the trailing hom. Find the
        // interval by the span's left edge (predBackboneOffset) and require the
        // right edge (succBackboneOffset) to fall strictly before the interval's
        // upper backbone anchor, so the whole chain nests inside one interval.
        const uint32_t leftOff = b.predBackboneOffset;
        const auto ub = std::upper_bound(bbOffset.begin(), bbOffset.end(), leftOff);
        if(ub == bbOffset.begin() || ub == bbOffset.end()) { ++dropped; continue; }
        const size_t i = size_t(ub - bbOffset.begin()) - 1;
        // The bracketing homs must not collide with the flanking backbone
        // anchors' EXPORTED positions. bbOffset[] is in the anchor-midpoint
        // frame, but every anchor is exported at position (midpoint - 1) (the
        // uniform k=2 half-length). So the backbone successor anchor bbOffset[i+1]
        // is exported at bbOffset[i+1]-1: a trailing hom at succBackboneOffset ==
        // bbOffset[i+1]-1 would map to the SAME read position as that backbone
        // anchor, giving the hom->backbone edge a zero-length (equal-position)
        // read that shasta2's LocalAssembly6 rejects (positionB > positionA).
        // Require succBackboneOffset + k/2 < bbOffset[i+1] on the right (k/2 = 1
        // for k=2, 0 for k=0). The left side needs no extra margin: the
        // predecessor anchor is likewise exported at bbOffset[i]-k/2, no closer
        // to the leading hom, so strict bbOffset[i] < predBackboneOffset already
        // avoids collision.
        const uint32_t hetKHalf = hetAnchorKHalf();
        if(!(bbOffset[i] < leftOff && b.succBackboneOffset + hetKHalf < bbOffset[i + 1])) {
            ++dropped;
            continue;
        }

        // Per-read backward-edge guard. The offset checks above use the backbone
        // MIDPOINT frame, which is shared across reads; they do not catch a
        // single member whose own rawPosition drifted onto/past that read's TRUE
        // position at a flanking backbone anchor. Drop such members from the
        // bracketing homs (only homs touch backbone anchors). bbAnchors[i] flanks
        // the leading hom (bbA_i -> leadHom); bbAnchors[i+1] flanks the trailing
        // hom (hom -> bbA_i+1).
        dropBackwardHomMembers(anchors, b.leadHom, hetKHalf,
            bbAnchors[i], /*homIsLeading=*/true, droppedBackwardMembers);
        dropBackwardHomMembers(anchors, b.hom, hetKHalf,
            bbAnchors[i + 1], /*homIsLeading=*/false, droppedBackwardMembers);
        // A hom must keep >=2 members to remain a valid anchor; if the drift
        // gutted one, drop the whole bubble.
        if(b.leadHom.members.size() < 2 || b.hom.members.size() < 2) {
            ++dropped;
            continue;
        }

        // Beyond the flanking backbone anchors, ANY of the bubble's k=2 anchors
        // (both bracketing homs AND the allele arms) can coincide with an
        // INTERVENING primary anchor: the backbone is a sparse subset of
        // primaries, so a primary can sit strictly inside the interval. A hom
        // column is homozygous, which is exactly what defines a primary anchor,
        // so homs are naturally drawn onto such positions. An allele ARM is not
        // exempt: heterozygosity constrains the allele BASE, not the read
        // POSITION, and a minority allele that is a coherent haplotype forms its
        // own k=50 primary anchor over exactly the same (read, position) markers.
        // Emitting the arm/hom as a second anchor id on a marker already owned by
        // a primary makes shasta2 read-following pair the two into a zero-length
        // assembly step and assert (the primary<->het "same read position"
        // collision). Drop the whole bubble when ANY of its anchors collides with
        // a primary, keeping only bubbles that sit cleanly inside the interval.
        bool primaryCollision =
            hetAnchorCollidesWithPrimary(b.leadHom, primarySet) ||
            hetAnchorCollidesWithPrimary(b.hom, primarySet);
        for(const auto& arm : b.alleles) {
            if(primaryCollision) break;
            if(hetAnchorCollidesWithPrimary(arm, primarySet)) primaryCollision = true;
        }
        if(primaryCollision) {
            ++dropped;
            ++droppedPrimaryCollision;
            continue;
        }

        b.plannedInterval = int32_t(i);
        ++plannedBubbles;
    }
}


// Pass 1.5 of the per-window het pipeline: merge coincident hom anchors.
//
// Why the collision arises. Adjacent bubbles N and N+1 in the same backbone
// interval are chained trailing-hom(N) -> leading-hom(N+1). The trailing hom
// sits at commonSucc_N (succBackboneOffset_N) and the leading hom at
// predPrev_{N+1} (predBackboneOffset_{N+1}). Flank linearity only guarantees >=2
// linear backbone bases between accepted SNPs, so when two SNPs are exactly 3 bp
// apart these two backbone offsets coincide: commonSucc_N == predPrev_{N+1}.
// That single backbone column is ONE abPOA node. Appending both homs as separate
// anchor ids builds two anchors on one node and stages a hom_N -> leadHom_{N+1}
// edge whose endpoints report equal positions on every shared read -- exactly
// the equal-position violation step 3 verifies and shasta2's LocalAssembly6
// asserts on (positionB > positionA).
//
// The two homs share the NODE but not their MEMBER SETS: cur.hom is built from
// SNP N's spanning reads, nxt.leadHom from SNP N+1's spanning reads, so a read
// spanning only one SNP appears in only one list. Reads present in BOTH are
// recovered at the same rawPosition (same node, same per-read position map), so
// the physically correct shared anchor is the UNION of the two member sets --
// every read passing through this backbone column belongs to this k=2 column.
//
// The fix. Union nxt.leadHom's members into cur.hom (dedup by OrientedReadId,
// asserting positions agree on the overlap), then mark N+1's leading hom as
// shared (sharedLeadFromBubble) so the append pass reuses N's trailing-hom
// anchor id instead of allocating a second one and the staging pass omits the
// redundant leadHom_{N+1} step. The chain becomes
//   ...arms_N -> sharedHom -> arms_{N+1}...
// with a single anchor on the shared node, keeping BOTH SNPs. This is the
// root-cause representation fix, not an export-time drop.
void mergeWindowCoincidentHoms(
    AnchorWindow& window,
    uint64_t& mergedHoms)
{
    if(window.hetBubbles.empty()) {
        return;
    }

    // Group planned bubbles by interval, then within each interval look for
    // adjacent pairs whose homs coincide. Previously this bucketed into
    // vector<vector<...>>(bbAnchors.size()) -- one bucket per BACKBONE anchor,
    // not per bubble -- so every window paid for an allocation and a scan
    // sized by its full backbone length (thousands of positions for a
    // full-journey window) even when it had zero or a handful of bubbles.
    // Sorting the actual bubbles by (plannedInterval, backboneOffset) gets the
    // same grouping and order, scaled by bubble count instead.
    vector<AnchorWindow::HetBubble*> planned;
    planned.reserve(window.hetBubbles.size());
    for(auto& b : window.hetBubbles) {
        if(b.plannedInterval >= 0) planned.push_back(&b);
    }
    if(planned.size() < 2) {
        return;
    }
    std::sort(planned.begin(), planned.end(),
        [](const AnchorWindow::HetBubble* a, const AnchorWindow::HetBubble* b) {
            if(a->plannedInterval != b->plannedInterval)
                return a->plannedInterval < b->plannedInterval;
            return a->backboneOffset < b->backboneOffset;
        });

    for(size_t bj = 0; bj + 1 < planned.size(); bj++) {
        AnchorWindow::HetBubble& cur = *planned[bj];
        AnchorWindow::HetBubble& nxt = *planned[bj + 1];
        if(cur.plannedInterval != nxt.plannedInterval) continue;
        // Coincident iff the two homs sit on the same backbone column, which
        // (nodeToBackboneOffset is a bijection) means the SAME POA node.
        if(cur.succBackboneOffset != nxt.predBackboneOffset) continue;

        // The two homs sit on one node but were built from DIFFERENT read
        // sets: cur.hom from SNP N's spanning reads, nxt.leadHom from SNP
        // N+1's spanning reads. A read spanning only one of the two SNPs is
        // in only one list, so the lists are not equal (and their allele-
        // grouped orderings differ). But because both lookups use the same
        // node and the same per-read position map, a read present in BOTH is
        // recovered at the SAME rawPosition. The physically correct shared
        // anchor is therefore the UNION of the two member sets (every read
        // through this backbone column belongs to this k=2 column); an
        // intersection would needlessly weaken the anchor. Merge nxt.leadHom
        // into cur.hom by OrientedReadId, deduplicating and asserting the
        // position invariant on the overlap.
        std::unordered_map<uint64_t, uint32_t> posByRead;
        posByRead.reserve(cur.hom.members.size() + nxt.leadHom.members.size());
        for(const auto& m : cur.hom.members)
            posByRead.emplace(m.orientedReadId.getValue(), m.rawPosition);
        for(const auto& m : nxt.leadHom.members) {
            const auto [it, inserted] =
                posByRead.emplace(m.orientedReadId.getValue(), m.rawPosition);
            // Read in both sets: same node + same per-read map => same
            // rawPosition. A mismatch would mean the coincidence assumption
            // is broken, so assert.
            if(!inserted) DINARA_ASSERT(it->second == m.rawPosition);
            else cur.hom.members.push_back(m);   // union: add the extra read
        }
        // Fold: bubble N+1 reuses bubble N's (now unioned) trailing hom as
        // its leading hom.
        nxt.sharedLeadFromBubble = int64_t(planned[bj] - window.hetBubbles.data());
        if(getenv("DINARA_APPEND_HET_DEBUG") != nullptr) {
            cout << "mergeDebug: window=" << window.windowId
                 << " cur.backboneOffset=" << cur.backboneOffset
                 << " cur.plannedInterval=" << cur.plannedInterval
                 << " nxt.backboneOffset=" << nxt.backboneOffset
                 << " nxt.plannedInterval=" << nxt.plannedInterval
                 << " storedIndex=" << nxt.sharedLeadFromBubble
                 << " window.hetBubbles.size()=" << window.hetBubbles.size() << endl;
        }
        ++mergedHoms;
    }
}


// Create the canonical/RC anchor pair for one het/hom anchor from its member
// list and record the assigned id on the anchor. No-op for empty member lists.
void appendHetAnchor(Shasta2Anchors& anchors, AnchorWindow::HetAnchor& a)
{
    if(a.members.empty()) return;
    if(a.members.size() < 2 && getenv("DINARA_APPEND_HET_DEBUG") != nullptr) {
        cout << "appendHetAnchorDebug: size=" << a.members.size()
             << " backboneOffset=" << a.backboneOffset
             << " isRef=" << a.isRef
             << " alleleBase=" << int(a.alleleBase)
             << " predBase=" << int(a.predBase) << endl;
    }
    vector<std::pair<OrientedReadId, uint32_t>> members;
    members.reserve(a.members.size());
    for(const auto& m : a.members)
        members.push_back({m.orientedReadId, m.rawPosition});
    a.anchorId = anchors.appendHetAnchorPair(members);
}


// Append anchor pairs for every planned bubble in a window: one per non-empty
// allele arm, plus BOTH bracketing hom anchors (leading + trailing). Only
// bubbles the planner accepted (plannedInterval >= 0) get anchors, and the
// planner already required both homs and >=1 arm, so every appended anchor is
// wired by the staging pass. Accumulates pair counts.
void appendWindowHetAnchors(
    Shasta2Anchors& anchors,
    AnchorWindow& window,
    uint64_t& hetAnchorPairs,
    uint64_t& homAnchorPairs)
{
    for(auto& bubble : window.hetBubbles) {
        if(bubble.plannedInterval < 0) continue;
        // Leading hom: brackets the bubble upstream (backbone -> leadHom edge).
        // If this bubble's leading hom was merged onto the preceding bubble's
        // trailing hom (coincident-hom, Pass 1.5), reuse that anchor id rather
        // than allocating a second anchor on the same POA node.
        if(bubble.sharedLeadFromBubble >= 0) {
            const AnchorWindow::HetBubble& prev =
                window.hetBubbles[size_t(bubble.sharedLeadFromBubble)];
            if(prev.plannedInterval < 0 && getenv("DINARA_APPEND_HET_DEBUG") != nullptr) {
                cout << "sharedLeadDebug: window=" << window.windowId
                     << " bubbleOff=" << bubble.backboneOffset
                     << " sharedLeadFromBubble=" << bubble.sharedLeadFromBubble
                     << " prev.backboneOffset=" << prev.backboneOffset
                     << " prev.plannedInterval=" << prev.plannedInterval
                     << " bubble.plannedInterval=" << bubble.plannedInterval << endl;
            }
            // Pass 1.5 only sets sharedLeadFromBubble for pairs where BOTH
            // bubbles are planned and coincident, so prev.hom was appended and
            // its anchorId is valid here.
            DINARA_ASSERT(prev.plannedInterval >= 0);
            DINARA_ASSERT(prev.hom.anchorId != invalid<Shasta2AnchorId>);
            bubble.leadHom.anchorId = prev.hom.anchorId;
        } else {
            appendHetAnchor(anchors, bubble.leadHom);
            ++homAnchorPairs;
        }
        for(auto& allele : bubble.alleles) {
            if(allele.members.empty()) continue;
            appendHetAnchor(anchors, allele);
            ++hetAnchorPairs;
        }
        // Trailing hom: brackets the bubble downstream (arms -> hom edge, then
        // hom -> next leadHom or hom -> backbone successor).
        appendHetAnchor(anchors, bubble.hom);
        ++homAnchorPairs;
    }
}


// Stage the window-local anchor-graph edges: for each backbone interval, a plain
// backbone edge when it holds no planned bubble, otherwise the hom-bracketed
// series chain
//   bbA_i -> leadHom_0 -> arms_0 -> hom_0
//         -> leadHom_1 -> arms_1 -> hom_1 -> ... -> hom_{n-1} -> bbA_{i+1}
// (the direct backbone edge is omitted, the chain replaces it). Each bubble is
// bracketed by its own leading and trailing hom so the two backbone-touching
// edges (bbA_i -> leadHom_0, hom_{n-1} -> bbA_{i+1}) are hom<->backbone edges
// whose read intersection is non-empty; minority allele arms never touch the
// backbone directly. Emits each forward edge plus its RC mirror. Accumulates
// edge/interval counts.
void stageWindowIntraEdges(
    AnchorWindow& window,
    uint64_t anchorCount,
    const vector<Shasta2AnchorId>& bbAnchors,
    const vector<uint32_t>& bbOffset,
    uint64_t& stagedBackbone,
    uint64_t& stagedHet,
    uint64_t& chainedIntervals)
{
    // Emit a forward edge plus its RC mirror (v^1 -> u^1) onto the window.
    // isHet controls how the constructor builds the edge: false = the normal
    // marker-based path (backbone-only endpoints), true = the direct
    // read-intersection path (any het/hom endpoint, which has no marker ordinal).
    auto emitEdge = [&](Shasta2AnchorId A, Shasta2AnchorId B,
                        uint32_t off, bool isHet, uint64_t& counter) {
        if(uint64_t(A) >= anchorCount || uint64_t(B) >= anchorCount) return;
        window.intraWindowEdges.push_back({A, B, off, isHet});
        ++counter;
        const Shasta2AnchorId Arc = Shasta2AnchorId(uint64_t(A) ^ 1ULL);
        const Shasta2AnchorId Brc = Shasta2AnchorId(uint64_t(B) ^ 1ULL);
        if(uint64_t(Arc) < anchorCount && uint64_t(Brc) < anchorCount) {
            window.intraWindowEdges.push_back({Brc, Arc, off, isHet});
            ++counter;
        }
    };

    // A step in an interval's series chain: either a SINGLE anchor (backbone
    // endpoint or hom separator) or a BUBBLE (>=2 allele arms).
    struct ChainStep {
        vector<Shasta2AnchorId> ids;   // 1 = single, >=2 = bubble arms
        uint32_t off = 0;              // representative backbone offset
    };

    // Group the planned bubbles by their assigned backbone interval. Bubbles
    // with plannedInterval < 0 were dropped and are ignored here.
    vector<vector<const AnchorWindow::HetBubble*>> byInterval(bbAnchors.size());
    for(const auto& bubble : window.hetBubbles) {
        if(bubble.plannedInterval < 0) continue;
        const size_t i = size_t(bubble.plannedInterval);
        if(i + 1 < bbAnchors.size()) byInterval[i].push_back(&bubble);
    }

    // For each backbone interval: no contained bubble -> plain backbone edge;
    // otherwise build the hom-bracketed series chain and OMIT the direct
    // backbone edge (the chain replaces the shortcut).
    for(uint64_t i = 0; i + 1 < bbAnchors.size(); i++) {
        auto& bubs = byInterval[i];
        if(bubs.empty()) {
            emitEdge(bbAnchors[i], bbAnchors[i + 1], 0, false, stagedBackbone);
            continue;
        }
        ++chainedIntervals;

        // Order the interval's bubbles by position.
        std::sort(bubs.begin(), bubs.end(),
            [](const AnchorWindow::HetBubble* a, const AnchorWindow::HetBubble* b) {
                return a->backboneOffset < b->backboneOffset;
            });

        // Build the ordered step list. Each bubble contributes THREE steps in
        // order: its leading hom, its allele arms, its trailing hom. The
        // flanking backbone anchors bracket the whole chain, and consecutive
        // bubbles chain trailing-hom -> next-leading-hom:
        //   bbA_i -> leadHom_0 -> arms_0 -> hom_0
        //         -> leadHom_1 -> arms_1 -> hom_1
        //         -> ... -> hom_{n-1} -> bbA_{i+1}
        // Every backbone-touching edge is now hom<->backbone (bbA_i -> leadHom_0
        // and hom_{n-1} -> bbA_{i+1}), so its read intersection is non-empty.
        // The append pass created both homs for every planned bubble, so their
        // anchorIds are valid here.
        vector<ChainStep> steps;
        steps.push_back({{bbAnchors[i]}, bbOffset[i]});
        for(size_t bj = 0; bj < bubs.size(); bj++) {
            const auto* b = bubs[bj];
            // Leading hom. When this bubble's leading hom was merged onto the
            // preceding bubble's trailing hom (coincident-hom, Pass 1.5), that
            // trailing hom is already in the step list as the previous bubble's
            // hom step, so omit the redundant leadHom step here: the chain runs
            // ...arms_prev -> sharedHom -> arms_this... through the single shared
            // anchor. (leadHom.anchorId still equals the shared id, set by the
            // append pass, but re-adding it would stage a self-edge.)
            if(b->sharedLeadFromBubble < 0 &&
               b->leadHom.anchorId != invalid<Shasta2AnchorId> &&
               uint64_t(b->leadHom.anchorId) < anchorCount) {
                steps.push_back({{b->leadHom.anchorId}, b->predBackboneOffset});
            }
            // Allele arms.
            ChainStep bubbleStep;
            bubbleStep.off = b->backboneOffset;
            for(const auto& allele : b->alleles) {
                if(allele.anchorId == invalid<Shasta2AnchorId>) continue;
                if(uint64_t(allele.anchorId) >= anchorCount) continue;
                bubbleStep.ids.push_back(allele.anchorId);
            }
            steps.push_back(std::move(bubbleStep));
            // Trailing hom.
            if(b->hom.anchorId != invalid<Shasta2AnchorId> &&
               uint64_t(b->hom.anchorId) < anchorCount) {
                steps.push_back({{b->hom.anchorId}, b->succBackboneOffset});
            }
        }
        steps.push_back({{bbAnchors[i + 1]}, bbOffset[i + 1]});

        // Connect consecutive steps. Every edge touches a het/hom anchor
        // (backbone->leadHom, leadHom->arms, arms->hom, hom->leadHom,
        // hom->backbone), so isHet=true. Offsets are strictly increasing by
        // construction (containment guarantees bbOffset[i] < predOff < snpOff <
        // succOff < bbOffset[i+1], bubbles sorted), so no fallback offset is used.
        for(size_t s = 0; s + 1 < steps.size(); s++) {
            const ChainStep& p = steps[s];
            const ChainStep& q = steps[s + 1];
            const uint32_t off = (q.off > p.off) ? (q.off - p.off) : 1;
            for(const Shasta2AnchorId A : p.ids)
                for(const Shasta2AnchorId B : q.ids)
                    emitEdge(A, B, off, true, stagedHet);
        }
    }
}


// Step 3 of the per-window het pipeline: monotonicity verification. Every staged
// intra-window edge A -> B must be strictly FORWARD on every read the two
// anchors share: the read's exported position at B must exceed its exported
// position at A. shasta2's LocalAssembly6 walks read-following in increasing
// position and rejects (asserts) a step whose target position is <= its source
// position, so a single backward edge here corrupts the whole export.
//
// The exported read position is (storedMidpoint - 1) for EVERY anchor class
// (primary and k=2 het/hom alike; see writeExternalAnchors and
// appendHetAnchorPair), so ordering by the stored midpoint position is identical
// to ordering by the exported position: the -1 cancels in the comparison.
// Backbone anchors already satisfy monotonicity by construction (journey order),
// and the shift (step 1) is uniform, so this pass is an assertion that the
// het-bubble wiring did not introduce a backward step. A violation is a bug in
// the planning/staging logic, so we throw rather than silently drop.
//
// Both anchor member lists are sorted by OrientedReadId, so shared reads are
// found by a linear merge. RC-mirror edges are covered because they are staged
// as their own edges and verified in the same pass.
void verifyWindowEdgeMonotonicity(
    const Shasta2Anchors& anchors,
    const AnchorWindow& window,
    uint64_t& checkedEdges,
    uint64_t& checkedReadSteps)
{
    for(const AnchorWindow::IntraWindowEdge& edge : window.intraWindowEdges) {
        const Shasta2Anchor anchorA = anchors[edge.anchorIdA];
        const Shasta2Anchor anchorB = anchors[edge.anchorIdB];
        ++checkedEdges;

        auto itA = anchorA.begin();
        auto itB = anchorB.begin();
        const auto endA = anchorA.end();
        const auto endB = anchorB.end();
        while(itA != endA && itB != endB) {
            if(itA->orientedReadId < itB->orientedReadId) { ++itA; continue; }
            if(itB->orientedReadId < itA->orientedReadId) { ++itB; continue; }
            // Shared read: exported ordering must be strictly forward. Both sides
            // subtract the same 1, so compare stored positions directly.
            if(!(itB->position > itA->position)) {
                throw runtime_error(
                    "Backward intra-window edge detected during monotonicity "
                    "verification: anchor " +
                    shasta2AnchorIdToString(edge.anchorIdA) + " -> anchor " +
                    shasta2AnchorIdToString(edge.anchorIdB) +
                    " on oriented read " +
                    to_string(itA->orientedReadId.getValue()) +
                    " has non-increasing position (" +
                    to_string(itA->position) + " -> " +
                    to_string(itB->position) + ").");
            }
            ++checkedReadSteps;
            ++itA;
            ++itB;
        }
    }
}


// Structural verification of a window's intra-window graph. Enforces the two
// invariants the staging pass (stageWindowIntraEdges) is supposed to guarantee:
//
//   (1) LINEAR: with each bubble's parallel allele arms condensed to a single
//       node, the window's anchor graph is a set of simple paths -- every
//       condensed node has in-degree <=1 and out-degree <=1, and every weakly
//       connected component is a path (equal source/sink counts, all nodes
//       reachable by walking sources). A backbone run with no bubble is a
//       straight chain; a bubble is a diamond (hom -> {arms} -> hom) that
//       condenses to one node. intraWindowEdges stores BOTH strands (each
//       forward edge plus its reversed RC mirror on id^1 anchors), so a correct
//       window is TWO node-disjoint paths (forward + reverse-complement). A node
//       with two distinct predecessors/successors means a fork/join that is not
//       a condensed bubble.
//
//   (2) HOM-FLANKED: every het allele arm is bracketed by hom/backbone anchors
//       on BOTH sides. Concretely, no edge connects one het allele arm directly
//       to another (arm->arm), and every arm has at least one incoming and one
//       outgoing edge (so it cannot dangle). Since staging only ever wires
//       arms as hom -> {arms} -> hom, any arm neighbor that is itself a het
//       allele arm (rather than a hom/backbone anchor) is a violation.
//
// hetFirst is the first het/hom anchor id; anchors below it are backbone/primary
// (treated as hom for flanking). het allele arms vs hom anchors are distinguished
// via the window's HetBubble records (allele anchorIds are arms; leadHom/hom
// anchorIds are homs). Throws on violation -- these are staging bugs.
void verifyWindowGraphStructure(
    const AnchorWindow& window,
    Shasta2AnchorId hetFirst,
    uint64_t& checkedWindows,
    uint64_t& checkedBubbles)
{
    using std::unordered_map;
    using std::unordered_set;
    using std::vector;

    if(window.intraWindowEdges.empty()) return;

    // Classify every het/hom anchor that appears in this window's bubbles:
    // allele arms are "arm" nodes; leadHom/hom are "hom" nodes. Anchors below
    // hetFirst are backbone/primary -> hom-like. Anything not recorded (should
    // not happen for het edges) is treated as hom-like so it cannot mask an
    // arm->arm edge.
    //
    // intraWindowEdges stores BOTH strands: each forward edge A->B plus its RC
    // mirror B^1 -> A^1 (reversed direction, strand-flipped ids, id^1). The two
    // strands are node-disjoint chains running in opposite directions, so we
    // must NOT merge them: node identity is the RAW anchor id (a forward-strand
    // path and an RC-strand path). The bubble records only carry one strand's
    // ids, so for each recorded id we register BOTH strands (id and id^1) as the
    // same class -- arms on the forward strand and their id^1 mirrors on the RC
    // strand are all arms; likewise for homs.
    auto flip = [](Shasta2AnchorId a) -> Shasta2AnchorId {
        return a ^ Shasta2AnchorId(1);
    };
    unordered_set<Shasta2AnchorId> armAnchors;   // raw ids, both strands
    unordered_set<Shasta2AnchorId> homAnchors;   // raw ids, both strands
    for(const AnchorWindow::HetBubble& b : window.hetBubbles) {
        bool any = false;
        for(const AnchorWindow::HetAnchor& a : b.alleles)
            if(a.anchorId != invalid<Shasta2AnchorId>) {
                armAnchors.insert(a.anchorId);
                armAnchors.insert(flip(a.anchorId));
                any = true;
            }
        if(b.leadHom.anchorId != invalid<Shasta2AnchorId>) {
            homAnchors.insert(b.leadHom.anchorId);
            homAnchors.insert(flip(b.leadHom.anchorId));
        }
        if(b.hom.anchorId != invalid<Shasta2AnchorId>) {
            homAnchors.insert(b.hom.anchorId);
            homAnchors.insert(flip(b.hom.anchorId));
        }
        if(any) ++checkedBubbles;
    }
    auto isArm = [&](Shasta2AnchorId a) {
        return armAnchors.find(a) != armAnchors.end();
    };

    // Assign condensed-node ids: node identity is the raw anchor id, so the two
    // strands stay separate. Each backbone/hom anchor is its own node; all arms
    // of the same bubble on the same strand collapse to one shared node id
    // (forward arms -> one node, their id^1 RC mirrors -> a second node).
    unordered_map<Shasta2AnchorId, uint64_t> node;   // raw anchorId -> node id
    uint64_t nextNode = 0;
    auto nodeOf = [&](Shasta2AnchorId a) -> uint64_t {
        auto it = node.find(a);
        if(it != node.end()) return it->second;
        const uint64_t id = nextNode++;
        node.emplace(a, id);
        return id;
    };
    // Collapse each bubble's arms onto one node id, per strand.
    for(const AnchorWindow::HetBubble& b : window.hetBubbles) {
        for(bool rc : {false, true}) {
            uint64_t shared = std::numeric_limits<uint64_t>::max();
            for(const AnchorWindow::HetAnchor& a : b.alleles) {
                if(a.anchorId == invalid<Shasta2AnchorId>) continue;
                const Shasta2AnchorId id = rc ? flip(a.anchorId) : a.anchorId;
                if(shared == std::numeric_limits<uint64_t>::max())
                    shared = nodeOf(id);
                else
                    node[id] = shared;
            }
        }
    }

    // Condensed directed graph.
    unordered_map<uint64_t, unordered_set<uint64_t>> succ, pred;
    unordered_set<uint64_t> nodes;
    for(const auto& e : window.intraWindowEdges) {
        // (2) HOM-FLANKED: reject a direct arm -> arm edge (two het alleles
        // wired together with no hom between them).
        if(isArm(e.anchorIdA) && isArm(e.anchorIdB)) {
            throw runtime_error(
                "Window " + to_string(window.windowId) +
                ": het allele arm connected directly to another het allele arm (" +
                shasta2AnchorIdToString(e.anchorIdA) + " -> " +
                shasta2AnchorIdToString(e.anchorIdB) +
                "); bubbles must be hom-flanked on both sides.");
        }
        const uint64_t na = nodeOf(e.anchorIdA);
        const uint64_t nb = nodeOf(e.anchorIdB);
        if(na == nb) continue;   // intra-bubble parallel arm edge (condensed away)
        nodes.insert(na);
        nodes.insert(nb);
        succ[na].insert(nb);
        pred[nb].insert(na);
    }

    // (2) HOM-FLANKED (dangle check): every arm node must have >=1 predecessor
    // and >=1 successor (a hom on each side).
    {
        unordered_set<uint64_t> armNodes;
        for(Shasta2AnchorId a : armAnchors) {
            auto it = node.find(a);
            if(it != node.end()) armNodes.insert(it->second);
        }
        for(uint64_t n : armNodes) {
            const bool hasPred = pred.count(n) && !pred[n].empty();
            const bool hasSucc = succ.count(n) && !succ[n].empty();
            if(!hasPred || !hasSucc) {
                throw runtime_error(
                    "Window " + to_string(window.windowId) +
                    ": het bubble arm is not hom-flanked on both sides (missing " +
                    string(!hasPred ? "leading" : "trailing") + " hom anchor).");
            }
        }
    }

    // (1) LINEAR: every weakly-connected component of the condensed graph must
    // be a simple path. Note intraWindowEdges stores BOTH strands (each forward
    // edge A->B plus its RC mirror Brc->Arc), so a correct window has TWO
    // node-disjoint paths (forward strand + reverse-complement); we therefore
    // check "each component is a path", not "one global path".
    //   - No node may have in-degree >1 or out-degree >1 (a fork or join means a
    //     branch that is not a condensed bubble -> non-linear).
    //   - #sources (in-deg 0) == #sinks (out-deg 0): each path has exactly one
    //     of each.
    //   - Walking forward from every source visits every node exactly once (no
    //     cycles, nothing stranded off a path).
    vector<uint64_t> sourceNodes;
    uint64_t sinks = 0;
    for(uint64_t n : nodes) {
        const size_t outDeg = succ.count(n) ? succ[n].size() : 0;
        const size_t inDeg  = pred.count(n) ? pred[n].size() : 0;
        if(outDeg > 1 || inDeg > 1) {
            // Diagnostic: identify the offending anchors and their class.
            string detail;
            for(const auto& e : window.intraWindowEdges) {
                const uint64_t na = node.count(e.anchorIdA) ? node[e.anchorIdA] : ~0ull;
                const uint64_t nb = node.count(e.anchorIdB) ? node[e.anchorIdB] : ~0ull;
                if(na != n && nb != n) continue;
                auto cls = [&](Shasta2AnchorId a) {
                    if(isArm(a)) return string("arm");
                    if(homAnchors.count(a)) return string("hom");
                    return (a >= hetFirst) ? string("het?") : string("bb");
                };
                detail += "\n    " + shasta2AnchorIdToString(e.anchorIdA) + "(" +
                          cls(e.anchorIdA) + ") -> " +
                          shasta2AnchorIdToString(e.anchorIdB) + "(" +
                          cls(e.anchorIdB) + ") isHet=" + to_string(e.isHet);
            }
            throw runtime_error(
                "Window " + to_string(window.windowId) +
                ": intra-window graph is not linear (condensed node has in-degree " +
                to_string(inDeg) + ", out-degree " + to_string(outDeg) +
                "); expected linear chains of backbone anchors and bubbles."
                " Edges touching this node:" + detail);
        }
        if(inDeg == 0)  sourceNodes.push_back(n);
        if(outDeg == 0) ++sinks;
    }
    if(sourceNodes.size() != sinks) {
        throw runtime_error(
            "Window " + to_string(window.windowId) +
            ": intra-window graph is not a set of paths (" +
            to_string(sourceNodes.size()) + " sources, " + to_string(sinks) +
            " sinks); expected equal counts.");
    }
    // Walk each source to its sink; total visited must equal node count.
    uint64_t visited = 0;
    unordered_set<uint64_t> seen;
    for(uint64_t src : sourceNodes) {
        uint64_t cur = src;
        while(true) {
            if(!seen.insert(cur).second) break;   // cycle / re-entry guard
            ++visited;
            if(!succ.count(cur) || succ[cur].empty()) break;
            cur = *succ[cur].begin();
        }
    }
    if(visited != nodes.size()) {
        throw runtime_error(
            "Window " + to_string(window.windowId) +
            ": intra-window graph has nodes off the backbone path(s) or a cycle (" +
            to_string(visited) + " of " + to_string(nodes.size()) +
            " condensed nodes reachable from sources).");
    }

    ++checkedWindows;
}


void maybeRunBubbleFinderDirectedSuperbubbleComparison(
    const mode3::AnchorGraph& anchorGraph,
    const std::vector<std::pair<
        mode3::AnchorGraph::vertex_descriptor,
        mode3::AnchorGraph::vertex_descriptor> >& onoderaSuperbubbles,
    uint64_t defaultThreadCount)
{
    if(envFlagIsDisabled("DINARA_BUBBLEFINDER_COMPARE")) {
        cout << "[SuperbubbleDetection] BubbleFinder comparison disabled by "
            "DINARA_BUBBLEFINDER_COMPARE." << endl;
        return;
    }

    std::string bubbleFinderBinary = "BubbleFinder";
    if(const char* value = ::getenv("DINARA_BUBBLEFINDER_BIN")) {
        if(*value != '\0') {
            bubbleFinderBinary = value;
        }
    }

    const uint64_t bubbleFinderThreadCount =
        envUintOrDefault("DINARA_BUBBLEFINDER_THREADS", defaultThreadCount);
    const uint64_t bubbleFinderSampleCount =
        envUintOrDefault("DINARA_BUBBLEFINDER_SAMPLE_COUNT", 10);

    const std::string gfaFileName = "AnchorGraph-BubbleFinder.gfa";
    const std::string outputFileName =
        "AnchorGraph-BubbleFinder.directed-superbubbles.txt";
    const std::string stdoutLogFileName = "AnchorGraph-BubbleFinder.stdout.log";
    const std::string stderrLogFileName = "AnchorGraph-BubbleFinder.stderr.log";

    std::unordered_set<std::string> nodeLabels;
    nodeLabels.reserve(num_vertices(anchorGraph));
    for(auto vp = boost::vertices(anchorGraph); vp.first != vp.second; ++vp.first) {
        const auto v = *vp.first;
        nodeLabels.insert(std::to_string(anchorGraph.getAnchorId(v)));
    }

    std::vector<std::pair<std::string, std::string> > directedEdges;
    directedEdges.reserve(num_edges(anchorGraph));
    for(auto ep = boost::edges(anchorGraph); ep.first != ep.second; ++ep.first) {
        const auto e = *ep.first;
        const auto v0 = source(e, anchorGraph);
        const auto v1 = target(e, anchorGraph);
        directedEdges.emplace_back(
            std::to_string(anchorGraph.getAnchorId(v0)),
            std::to_string(anchorGraph.getAnchorId(v1)));
    }

    {
        std::ofstream gfaFile(gfaFileName);
        if(!gfaFile) {
            cout << "[SuperbubbleDetection] BubbleFinder comparison skipped: "
                "failed to open " << gfaFileName << " for writing." << endl;
            return;
        }

        gfaFile << "H\tVN:Z:1.0\n";
        for(const std::string& label: nodeLabels) {
            gfaFile << "S\t" << label << "\t*\n";
        }
        for(const auto& [u, v]: directedEdges) {
            gfaFile << "L\t" << u << "\t+\t" << v << "\t+\t0M\n";
        }
    }

    const std::string command =
        shellQuote(bubbleFinderBinary) +
        " directed-superbubbles --gfa-directed -g " +
        shellQuote(gfaFileName) +
        " -o " +
        shellQuote(outputFileName) +
        " -j " +
        std::to_string(bubbleFinderThreadCount) +
        " > " +
        shellQuote(stdoutLogFileName) +
        " 2> " +
        shellQuote(stderrLogFileName);

    cout << "[SuperbubbleDetection] Running BubbleFinder directed-superbubbles "
        "for AnchorGraph comparison." << endl;
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        cout << "[SuperbubbleDetection] BubbleFinder comparison failed with exit code "
            << errorCode << ". Command stderr: " << stderrLogFileName << endl;
        return;
    }

    std::set<BubbleEndpointPair> onoderaPairs;
    onoderaPairs.clear();
    for(const auto& p: onoderaSuperbubbles) {
        const auto entrance = std::to_string(anchorGraph.getAnchorId(p.first));
        const auto exit = std::to_string(anchorGraph.getAnchorId(p.second));
        onoderaPairs.insert(canonicalBubblePair(entrance, exit));
    }

    std::ifstream bubbleFinderOutput(outputFileName);
    if(!bubbleFinderOutput) {
        cout << "[SuperbubbleDetection] BubbleFinder comparison skipped: "
            "failed to open " << outputFileName << " for reading." << endl;
        return;
    }

    uint64_t declaredPairCount = 0;
    if(!(bubbleFinderOutput >> declaredPairCount)) {
        cout << "[SuperbubbleDetection] BubbleFinder comparison skipped: "
            "invalid output header in " << outputFileName << "." << endl;
        return;
    }

    std::set<BubbleEndpointPair> bubbleFinderPairs;
    std::vector<BubbleEndpointPair> bubbleFinderSamples;
    uint64_t parsedPairCount = 0;
    uint64_t ignoredAuxiliaryPairCount = 0;
    std::string u;
    std::string v;
    while(bubbleFinderOutput >> u >> v) {
        parsedPairCount++;
        u = normalizeBubbleEndpoint(u);
        v = normalizeBubbleEndpoint(v);
        if(nodeLabels.find(u) == nodeLabels.end() ||
           nodeLabels.find(v) == nodeLabels.end()) {
            ignoredAuxiliaryPairCount++;
            continue;
        }
        const auto pair = canonicalBubblePair(u, v);
        bubbleFinderPairs.insert(pair);
        if(bubbleFinderSamples.size() < bubbleFinderSampleCount) {
            bubbleFinderSamples.push_back(pair);
        }
    }

    uint64_t sharedPairCount = 0;
    for(const auto& pair: onoderaPairs) {
        if(bubbleFinderPairs.find(pair) != bubbleFinderPairs.end()) {
            sharedPairCount++;
        }
    }

    const uint64_t onoderaOnlyCount = onoderaPairs.size() - sharedPairCount;
    const uint64_t bubbleFinderOnlyCount = bubbleFinderPairs.size() - sharedPairCount;

    cout << "[SuperbubbleDetection] BubbleFinder comparison: "
        << "Onodera=" << onoderaPairs.size()
        << ", BubbleFinder=" << bubbleFinderPairs.size()
        << ", shared=" << sharedPairCount
        << ", onoderaOnly=" << onoderaOnlyCount
        << ", bubbleFinderOnly=" << bubbleFinderOnlyCount
        << ", ignoredAuxiliary=" << ignoredAuxiliaryPairCount;
    if(parsedPairCount != declaredPairCount) {
        cout << ", declaredPairs=" << declaredPairCount
            << ", parsedPairs=" << parsedPairCount;
    }
    cout << "." << endl;

    if(!bubbleFinderSamples.empty()) {
        cout << "[SuperbubbleDetection] BubbleFinder sample pairs:";
        for(const auto& p: bubbleFinderSamples) {
            cout << " " << p.first << "->" << p.second;
        }
        cout << endl;
    }
}

} // namespace




int main(int argumentCount, const char** arguments)
{
    try {

        dinara::main::main(argumentCount, arguments);

    } catch(const boost::program_options::error_with_option_name& e) {
        cout << "Invalid option: " << e.what() << endl;
        return 1;
    } catch (const runtime_error& e) {
        cout << timestamp << e.what() << endl;
        return 2;
    } catch (const std::bad_alloc& e) {
        cout << timestamp << e.what() << endl;
        cout << "Memory allocation failure." << endl;
        cout << "This assembly requires more memory than available." << endl;
        cout << "Rerun on a larger machine." << endl;
        return 2;
    } catch (const exception& e) {
        cout << timestamp << e.what() << endl;
        return 3;
    } catch (...) {
        cout << timestamp << "Terminated after catching a non-standard exception." << endl;
        return 4;
    }
    return 0;
}



void dinara::main::segmentFaultHandler(int)
{
    char message[] = "\nA segment fault occurred. Please report it by filing an "
        "issue on the Dinara repository and attaching the entire log output. "
        "To file an issue, point your browser to https://github.com/kokyriakidis/dinara/issues\n";
    ::write(fileno(stderr), message, sizeof(message));
    ::_exit(1);
}

void dinara::main::setupSegmentFaultHandler()
{
    struct sigaction action;
    ::memset(&action, 0, sizeof(action));
    action.sa_handler = &segmentFaultHandler;
    sigaction(SIGSEGV, &action, 0);
}


void dinara::main::main(int argumentCount, const char** arguments)
{
    setupSegmentFaultHandler();

    // Parse command line options and the configuration file, if one was specified.
    AssemblerOptions assemblerOptions(argumentCount, arguments);

    // Check that we have a valid command.
    auto it = commands.find(assemblerOptions.commandLineOnlyOptions.command);
    if(it ==commands.end()) {
        const string message = "Invalid command " + assemblerOptions.commandLineOnlyOptions.command;
        listCommands();
        throw runtime_error(message);
    }

    // Execute the requested command.
    if(assemblerOptions.commandLineOnlyOptions.command == "assemble") {
        assemble(assemblerOptions, argumentCount, arguments);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "cleanupBinaryData") {
        cleanupBinaryData(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "saveBinaryData") {
        saveBinaryData(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "explore") {
        explore(assemblerOptions);
        return;
    } else if(assemblerOptions.commandLineOnlyOptions.command == "listCommands") {
        listCommands();
        return;
    }

    // We already checked for a valid command above, so if we get here
    // the above logic is missing code for one of the valid commands.
    DINARA_ASSERT(0);

}




// Implementation of --command assemble.
void dinara::main::assemble(
    const AssemblerOptions& assemblerOptions,
    int argumentCount, const char** arguments)
{
    DINARA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "assemble");


    // Various checks for option validity.

#ifdef DINARA_LONG_MARKERS
    // With capacity-128 Kmers (256-bit KmerId), k can be up to 126.
    constexpr uint64_t maxK = 126;
#else
    constexpr uint64_t maxK = 62;
#endif
    if(assemblerOptions.kmersOptions.k > maxK or assemblerOptions.kmersOptions.k < 6) {
        throw runtime_error("Invalid value specified for --Kmers.k. Must be between 6 and " +
            to_string(maxK) + ".");
    }

    if((assemblerOptions.kmersOptions.k % 2) == 1) {
        throw runtime_error("Invalid value specified for --Kmers.k. Must be even.");
    }

    // Check that we have at least one input file.
    if(assemblerOptions.commandLineOnlyOptions.inputFileNames.empty()) {
        throw runtime_error("Specify at least one input file "
            "using command line option --input.");
    }

    if( assemblerOptions.alignOptions.alignMethod <  0 or
        assemblerOptions.alignOptions.alignMethod == 2 or
        assemblerOptions.alignOptions.alignMethod >  6) {
        throw runtime_error("Align method " + to_string(assemblerOptions.alignOptions.alignMethod) +
            " is not valid. Valid options are 0 through 6 except 2.");
    }

    // Find absolute paths of the input files.
    // We will use them below after changing directory to the output directory.
    vector<string> inputFileAbsolutePaths;
    for(const string& inputFileName: assemblerOptions.commandLineOnlyOptions.inputFileNames) {
        if(!std::filesystem::exists(inputFileName)) {
            throw runtime_error("Input file not found: " + inputFileName);
        }
        if(!std::filesystem::is_regular_file(inputFileName)) {
            throw runtime_error("Input file is not a regular file: " + inputFileName);
        }
        inputFileAbsolutePaths.push_back(filesystem::getAbsolutePath(inputFileName));
    }

    // Create the assembly directory. If it exists, stop.
    bool exists = std::filesystem::exists(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
    if (exists) {
        throw runtime_error(
            assemblerOptions.commandLineOnlyOptions.assemblyDirectory +
            " already exists. Remove it first \n"
            "or use --assemblyDirectory to specify a different assembly directory."
        );
    } else {
        DINARA_ASSERT(std::filesystem::create_directory(assemblerOptions.commandLineOnlyOptions.assemblyDirectory));
    }

    // Make the assembly directory current.
    std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);

    // Open the performance log.
    openPerformanceLog("performance.log");
    performanceLog << timestamp << "Assembly begins." << endl;

    // Open stdout.log and "tee" (duplicate) stdout to it.
    if(not assemblerOptions.commandLineOnlyOptions.suppressStdoutLog) {
        dinaraLog.open("stdout.log");
        tee.duplicate(cout, dinaraLog);
    }

    // Echo out the command line options.
    cout << timestamp << "Assembly begins.\nCommand line:" << endl;
    for(int i=0; i<argumentCount; i++) {
        cout << arguments[i] << " ";
    }
    cout << endl;

    // Set up the run directory as required by the memoryMode and memoryBacking options.
    size_t pageSize = 0;
    string dataDirectory;
    setupRunDirectory(
        assemblerOptions.commandLineOnlyOptions.memoryMode,
        assemblerOptions.commandLineOnlyOptions.memoryBacking,
        pageSize,
        dataDirectory);

    // Write out the option in effect to dinara.conf.
    {
        ofstream configurationFile("dinara.conf");
        assemblerOptions.write(configurationFile);
    }
    cout << "For options in use for this assembly, see dinara.conf in the assembly directory." << endl;

    // Create the Assembler.
    Assembler assembler(dataDirectory, true, assemblerOptions.readsOptions.representation, pageSize);
    assembler.assemblerInfo->readGraphCreationMethod = assemblerOptions.readGraphOptions.creationMethod;
    assembler.assemblerInfo->assemblyMode = assemblerOptions.assemblyOptions.mode;
    assembler.minMultiNodeChainSupport = assemblerOptions.readGraphOptions.minMultiNodeChainSupport;
    assembler.minIsolatedSiteSupport = assemblerOptions.readGraphOptions.minIsolatedSiteSupport;
    assembler.assemblerInfo->variantClusteringMinOccurrences = assemblerOptions.variantClusteringOptions.minOccurrences;
    assembler.assemblerInfo->variantClusteringMinSeparation = assemblerOptions.variantClusteringOptions.minSeparation;

    // Run the assembly.
    assemble(assembler, assemblerOptions, inputFileAbsolutePaths);

    cout << timestamp << "Assembly ends." << endl;
    performanceLog << timestamp << "Assembly ends." << endl;
}



// Set up the run directory as required by the memoryMode and memoryBacking options.
void dinara::main::setupRunDirectory(
    const string& memoryMode,
    const string& memoryBacking,
    size_t& pageSize,
    string& dataDirectory
    )
{

    if(memoryMode == "anonymous") {

        if(memoryBacking == "disk") {

            // This combination is meaningless.
            throw runtime_error("\"--memoryMode anonymous\" is not allowed in combination "
                "with \"--memoryBacking disk\".");

        } else if(memoryBacking == "4K") {

            // Anonymous memory on 4KB pages.
            // This combination is the default.
            // It does not require root privilege.
            dataDirectory = "";
            pageSize = 4096;

        } else if(memoryBacking == "2M") {

            // Anonymous memory on 2MB pages.
            // This may require root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            // Root privilege is not required if 2M pages have already
            // been set up as required.
#ifdef __APPLE__
            throw runtime_error("Option --memoryBacking 2M is not supported on macOS.");
#else
            setupHugePages();
            pageSize = 2 * 1024 * 1024;
#endif

        } else {
            throw runtime_error("Invalid value specified for --memoryBacking: " + memoryBacking +
                "\nValid values are: disk, 4K, 2M.");
        }

    } else if(memoryMode == "filesystem") {

        if(memoryBacking == "disk") {

            // Binary files on disk.
            // This does not require root privilege.
            DINARA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 4096;

        } else if(memoryBacking == "4K") {

            // Binary files on the tmpfs filesystem
            // (filesystem in memory backed by 4K pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
#ifdef __APPLE__
            throw runtime_error("Option --memoryMode filesystem --memoryBacking 4K is not supported on macOS.");
#else
            DINARA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 4096;
            const string command = "sudo mount -t tmpfs -o size=0 tmpfs Data";
            const int errorCode = ::system(command.c_str());
            if(errorCode != 0) {
                throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                    " running command: " + command);
            }
#endif

        } else if(memoryBacking == "2M") {

            // Binary files on the hugetlbfs filesystem
            // (filesystem in memory backed by 2M pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
#ifdef __APPLE__
            throw runtime_error("Option --memoryMode filesystem --memoryBacking 2M is not supported on macOS.");
#else
            setupHugePages();
            DINARA_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 2 * 1024 * 1024;
            const uid_t userId = ::getuid();
            const gid_t groupId = ::getgid();
            const string command = "sudo mount -t hugetlbfs -o pagesize=2M"
                ",uid=" + to_string(userId) +
                ",gid=" + to_string(groupId) +
                " none Data";
            const int errorCode = ::system(command.c_str());
            if(errorCode != 0) {
                throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                    " running command: " + command);
            }
#endif

        } else {
            throw runtime_error("Invalid value specified for --memoryBacking: " + memoryBacking +
                "\nValid values are: disk, 4K, 2M.");
        }

    } else {
        throw runtime_error("Invalid value specified for --memoryMode: " + memoryMode +
            "\nValid values are: anonymous, filesystem.");
    }
}



// This runs the entire assembly, under the following assumptions:
// - The current directory is the run directory.
// - The Data directory has already been created and set up, if necessary.
// - The input file names are either absolute,
//   or relative to the run directory, which is the current directory.
void dinara::main::assemble(
    Assembler& assembler,
    const AssemblerOptions& assemblerOptions,
    vector<string> inputFileNames)
{
    const auto steadyClock0 = std::chrono::steady_clock::now();
    const auto userClock0 = boost::chrono::process_user_cpu_clock::now();
    const auto systemClock0 = boost::chrono::process_system_cpu_clock::now();

    // Adjust the number of threads, if necessary.
    uint64_t threadCount = assemblerOptions.commandLineOnlyOptions.threadCount;
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "This assembly will use " << threadCount << " threads." << endl;

    // If --saveBinaryData was requested,
    // create the directory where binary data will be saved.
    if (assemblerOptions.commandLineOnlyOptions.saveBinaryData) {
        assembler.createSaveBinaryDataDirectory(assemblerOptions.commandLineOnlyOptions.memoryMode);
    }


    // Add reads from the specified input files.
    performanceLog << timestamp << "Begin loading reads from " << inputFileNames.size() << " files." << endl;
    const auto t0 = steady_clock::now();
    for(const string& inputFileName: inputFileNames) {

        assembler.addReads(
            inputFileName,
            assemblerOptions.readsOptions.minReadLength,
            assemblerOptions.readsOptions.noCache,
            threadCount);
    }

    if(assembler.getReads().readCount() == 0) {
        throw runtime_error("There are no input reads.");
    }
    const uint64_t averageReadLength =
        assembler.getReads().getTotalBaseCount() / assembler.getReads().readCount();
    cout << "Average read length: " << averageReadLength << " bp." << endl;



    // If requested, increase the read length cutoff
    // to reduce coverage to the specified amount.
    if (assemblerOptions.readsOptions.desiredCoverage > 0) {
        // Write out the read length histogram using provided minReadLength.
        assembler.histogramReadLength("ExtendedReadLengthHistogram.csv");

        const auto newMinReadLength = assembler.adjustCoverageAndGetNewMinReadLength(
            assemblerOptions.readsOptions.desiredCoverage);

        const auto oldMinReadLength = uint64_t(assemblerOptions.readsOptions.minReadLength);

        if (newMinReadLength == 0ULL) {
            throw runtime_error(
                "With Reads.minReadLength " +
                to_string(assemblerOptions.readsOptions.minReadLength) +
                ", total available coverage is " +
                to_string(assembler.getReads().getTotalBaseCount()) +
                ", less than desired coverage " +
                to_string(assemblerOptions.readsOptions.desiredCoverage) +
                ". Try reducing Reads.minReadLength if appropriate or get more coverage."
            );
        }

        // Adjusting coverage should only ever reduce coverage if necessary.
        DINARA_ASSERT(newMinReadLength >= oldMinReadLength);
    }

    assembler.computeReadIdsSortedByName();
    assembler.histogramReadLength("ReadLengthHistogram.csv");

    const auto t1 = steady_clock::now();
    performanceLog << timestamp << "Done loading reads from " << inputFileNames.size() << " files." << endl;
    performanceLog << "Read loading took " << seconds(t1-t0) << "s." << endl;

    // Find duplicate reads and handle them according to the setting
    // of --Reads.handleDuplicates. The default option is "useOneCopy".
    assembler.findDuplicateReads(assemblerOptions.readsOptions.handleDuplicates);

    // Marker generation method selection.
    //
    // The SIMD minimizer path is the main path and is taken by default because
    // Kmers.useHifiasmMinimizers defaults to true: it generates markers from
    // hifiasm's own no-HPC sketcher plus the overlap-path minimizer filter, so
    // the marker seeds match the seeds hifiasm uses for overlaps. This is what
    // pairs with the PAF/hifiasm overlap path.
    //
    // Within the SIMD path the position source is:
    //   - hifiasm sketcher            when useHifiasmMinimizers (default)
    //   - simd-minimizers library     when only useSimdClosedSyncmers is set
    // The legacy k-mer marker method is used only when BOTH flags are false.
    const bool useSimdMinimizerPath =
        assemblerOptions.kmersOptions.useHifiasmMinimizers ||
        assemblerOptions.kmersOptions.useSimdClosedSyncmers;
    if(useSimdMinimizerPath) {
        // // Use SIMD-accelerated closed syncmers for initial marker generation (no filtering).
        // assembler.findMarkersSimdClosedSyncmers(
        //     threadCount,
        //     assemblerOptions.kmersOptions.k,
        //     assemblerOptions.kmersOptions.syncmerS);

        // Use SIMD-accelerated minimizers instead of closed syncmers.
        // For hifiasm-like behavior with k=w, use syncmerS parameter as window size.
        // Density ≈ 2/w (smaller w = denser sampling, larger w = sparser sampling)
        // Optionally build hifiasm's overlap-path minimizer filter over the
        // input reads so markers match the seeds hifiasm uses for overlaps: a
        // no-HPC high-occurrence k-mer filter plus distance subsampling. Built
        // once here and shared read-only across all marker threads; the handle
        // is a standalone hash of k-mer values (no read store), so it is safe
        // to share and outlives this call until we destroy it below.
        //
        // Only meaningful on the hifiasm minimizer path. When active, the
        // redundant downstream marker frequency prune (applyKmerCountFilter) is
        // skipped so the marker set is exactly the hf + sample_dist sketch
        // output (true parity with the overlap path).
        hifiasm_filter_t* hifiasmMarkerFilter = nullptr;
        const bool wantHifiasmMarkerFilter =
            assemblerOptions.kmersOptions.useHifiasmMinimizers &&
            !inputFileNames.empty();
        if(wantHifiasmMarkerFilter) {
            performanceLog << timestamp
                << "Building hifiasm overlap-path minimizer filter (no-HPC, k=w="
                << assemblerOptions.kmersOptions.k << ") over "
                << inputFileNames.size() << " input file(s)." << endl;

            vector<const char*> markerReadFiles;
            markerReadFiles.reserve(inputFileNames.size());
            for(const string& f: inputFileNames) {
                markerReadFiles.push_back(f.c_str());
            }

            hifiasm_filter_opt_t filterOpt = {};
            filterOpt.threads = int(threadCount);
            filterOpt.k_mer_length = assemblerOptions.kmersOptions.k;
            filterOpt.mz_win = assemblerOptions.kmersOptions.k; // w == k
            filterOpt.is_hpc = 0;        // dinara markers are no-HPC
            filterOpt.min_read_len = -1; // keep all reads (match dinara's set)

            hifiasmMarkerFilter = hifiasm_build_filter(
                markerReadFiles.data(), int(markerReadFiles.size()), &filterOpt);
            if(hifiasmMarkerFilter == nullptr) {
                throw runtime_error(
                    "Failed to build hifiasm minimizer filter for markers.");
            }
            performanceLog << timestamp
                << "hifiasm minimizer filter built." << endl;
        }

        assembler.findMarkersSimdMinimizers(
            threadCount,
            assemblerOptions.kmersOptions.k,
            assemblerOptions.kmersOptions.k,  // Using kmer length as window size w
            assemblerOptions.kmersOptions.useHifiasmMinimizers,
            hifiasmMarkerFilter,
            assemblerOptions.kmersOptions.hifiasmMarkerSampleDist);

        // Whether the overlap-path filter was actually applied to the markers.
        const bool hifiasmMarkerFilterApplied = (hifiasmMarkerFilter != nullptr);

        // Compute histogram using the pre-calculated KmerIds. Still needed even
        // when the hf filter is active: it fills coverageHet/coverageHom, which
        // phasing and the chaining-frequency cutoff below depend on.
        assembler.countKmersFromMarkerKmerIds(threadCount);
        
        // Retrieve peak and set thresholds.
        // Hifiasm hard-filters k-mers above max_kmer_cnt (default 2000) during
        // sketching. K-mers between highFreqThreshold (coverageHet * 1.667)
        // and this cutoff are kept but downsampled per-streak during chaining.
        const uint64_t coverageHet = assembler.assemblerInfo->kmerDistributionInfo.coverageHet;
        const uint64_t coverageHom = assembler.assemblerInfo->kmerDistributionInfo.coverageHom;
        const uint64_t minFreq = 2;
        // Keep all markers with freq >= 2 in the marker arrays (for complete
        // journeys/anchors). High-frequency kmers are excluded from chaining
        // via maxChainingFreq in the inverted index.
        const uint64_t maxFreq = std::numeric_limits<uint64_t>::max();
        const bool removePalindromicKmers = true;
        uint64_t distinctKmerCount = 0;
        for(uint64_t bucketId=0; bucketId<assembler.kmerCounter->kmerIdFrequencies.size(); bucketId++) {
            distinctKmerCount += assembler.kmerCounter->kmerIdFrequencies[bucketId].size();
        }

        cout << "Analyzing " << distinctKmerCount << " distinct minimizer k-mers." << endl;
        cout << "Filtering minimizers: het coverage=" << coverageHet
             << ", hom coverage=" << coverageHom << "." << endl;
        cout << "Keeping k-mers with frequency >= " << minFreq;
        if(removePalindromicKmers) {
            cout << " and excluding palindromic k-mers";
        }
        // Effective chaining cutoff: min(coverageHom * 5, maxChainingFreq).
        // Matches hifiasm's min(hom_cov * 5, max_kmer_cnt) logic.
        const uint64_t effectiveChainingFreq = min(coverageHom * 5,
            uint64_t(assemblerOptions.overlapCandidatesOptions.maxChainingFreq));
        cout << ". Chaining limited to frequency <= " << effectiveChainingFreq << "." << endl;
             
        // Prune the existing minimizer markers in-place.
        // applyKmerCountFilter keeps a marker only if:
        // - its canonical k-mer frequency is in the inclusive range [minFreq, maxFreq],
        // - the k-mer is not palindromic/self-reverse-complementary (default),
        // - the k-mer is not a short-period tandem repeat (filterRepeatKmers), and
        // - the k-mer is not low-complexity by distinct sub-k-mer count
        //   (filterLowComplexity).
        // The last two apply the SAME predicates the marker-graph vertex filters
        // (filterMarkerGraphVerticesByRepeatKmers /
        // filterMarkerGraphVerticesByDistinctSubkmerCount) use, but at the
        // minimizer stage, so repeat/low-complexity minimizers never seed a
        // marker or marker-graph vertex. The function rebuilds both markers and
        // markerKmerIds from the pre-filtered arrays, preserving only marker
        // positions whose matching k-mer id passes.
        const bool filterRepeatKmers = true;
        const bool filterLowComplexity = true;
        if(hifiasmMarkerFilterApplied) {
            // The hifiasm overlap-path filter (high-occurrence k-mer filter +
            // distance subsampling) already selected the markers, exactly as
            // hifiasm selects overlap seeds. Running applyKmerCountFilter on top
            // would prune a second, different way and break that parity, so it
            // is intentionally skipped here. coverageHet/coverageHom from
            // countKmersFromMarkerKmerIds above are still available for phasing
            // and the chaining-frequency cutoff.
            cout << "Skipping downstream marker frequency filter: markers already "
                 << "filtered by hifiasm's overlap-path minimizer filter." << endl;
        } else {
            assembler.applyKmerCountFilter(
                minFreq, maxFreq, threadCount, removePalindromicKmers,
                filterRepeatKmers, filterLowComplexity);
        }

        // writeReadMarkerGapDiagnostic("afterFrequencyFilter", ReadId(3729));

        // The marker filter handle is no longer needed once markers are built.
        if(hifiasmMarkerFilter != nullptr) {
            hifiasm_filter_destroy(hifiasmMarkerFilter);
            hifiasmMarkerFilter = nullptr;
        }

        // Initialize KmerChecker for HttpServer diagnostics (optional).
        cout << "Initializing KmerChecker for diagnostics." << endl;
        assembler.createKmerChecker(assemblerOptions.kmersOptions, threadCount);
            
    } else {
        // Use the default k-mer based method.
        // Initialize the KmerChecker, which has the information needed
        // to decide if a k-mer is a marker.
        assembler.createKmerChecker(assemblerOptions.kmersOptions, threadCount);

        // Find the markers in the reads.
        assembler.findMarkers(threadCount);

        // Compute marker KmerIds (required for LowHash and alignment).
        // The SIMD path already creates these, but findMarkers does not.
        assembler.computeMarkerKmerIds(threadCount);

        // Compute k-mer histogram to get coverageHet (needed by phasing).
        // The SIMD path does this via countKmersFromMarkerKmerIds.
        assembler.countKmersFromMarkerKmerIds(threadCount);
    }

    // Filter reads whose marker span covers less than the threshold fraction
    // of the read length. Reads with sparse markers contribute noise.
    if(assemblerOptions.kmersOptions.minMarkerSpanFraction > 0.0) {
        assembler.filterReadsByMarkerSpanCoverage(
            assemblerOptions.kmersOptions.minMarkerSpanFraction, threadCount);
    }

    assembler.initiateSaveBinaryData(&Assembler::saveMarkers);



    // ========================================================================
    // HIFIASM MAX_N_CHAIN CALCULATION (Per-Read Overlap Limiting)
    // ========================================================================
    // Reference: Hifiasm CommandLines.cpp:413, anchor.cpp:191-220
    //
    // Compute max_n_chain = max(MIN_N_CHAIN, hom_cov * high_factor)
    //
    // Purpose: Limit total number of overlaps kept per read across ALL partners
    // - Applied AFTER chaining all read pairs (not during DP chaining)
    // - Prevents memory explosion from reads with many partners
    // - Keeps top N overlaps by score, grouped by overlap type
    //
    // Hifiasm defaults:
    //   high_factor = 5.0 (CommandLines.cpp:271)
    //   MIN_N_CHAIN = 100 (CommandLines.h:28)
    //
    // Example: coverage = 30x → max(100, 30*5) = 150 overlaps per read
    //
    const uint64_t coverageHet = assembler.assemblerInfo->kmerDistributionInfo.coverageHet;
    const uint64_t maxChainLimit = std::max<uint64_t>(
        assemblerOptions.overlapCandidatesOptions.invertedIndexMinNChain,     // MIN_N_CHAIN = 100
        uint64_t(double(coverageHet) * assemblerOptions.overlapCandidatesOptions.invertedIndexHighFactor + 0.499));  // hom_cov * 5.0

    // Build the inverted index: for each k-mer, store which reads contain
    // it and at what position. Used for overlap candidate discovery and
    // for chaining marker matches into alignments.
    assembler.buildInvertedIndex(threadCount);

    // // Detect palindromic reads — reads whose reverse complement aligns well
    // // to themselves. These cause spurious overlaps because both strands map
    // // to the same genomic location. Flagged reads are excluded from overlap
    // // candidate discovery.
    // if(!assemblerOptions.readsOptions.palindromicReads.skipFlagging) {
    //     assembler.flagPalindromicReads(
    //         assemblerOptions.overlapCandidatesOptions.driftRateTolerance,
    //         assemblerOptions.overlapCandidatesOptions,
    //         assemblerOptions.readsOptions.palindromicReads.alignedFractionThreshold,
    //         assemblerOptions.readsOptions.palindromicReads.maxErrorRate,
    //         threadCount);
    // }

    // Discover overlapping read pairs and chain their shared markers into
    // alignments. Each chain represents a collinear sequence of marker
    // matches between two reads, scored by the number of shared seeds.
    // Two paths are supported:
    // - hifiasm: read overlaps are generated by the bundled hifiasm library
    //   (writing hifiasm.ovlp.paf), then imported via the PAF path.
    // - PAF: read pairs are imported from an external PAF file, then chained.
    // - Inverted index: read pairs are discovered by shared k-mer lookups,
    //   then chained. This is the default path.
    string overlapsPafFile = assemblerOptions.commandLineOnlyOptions.overlapsFromPafFile;
    if(assemblerOptions.commandLineOnlyOptions.overlapsFromHifiasm) {
        if(!overlapsPafFile.empty()) {
            throw runtime_error(
                "--overlapsFromHifiasm and --overlapsFromPafFile are mutually exclusive.");
        }
        if(inputFileNames.empty()) {
            throw runtime_error("--overlapsFromHifiasm requires at least one --input read file.");
        }

        // Run hifiasm candidate overlap detection on the same input reads.
        // The current directory is the run directory, so the PAF is written
        // there as hifiasm.ovlp.paf. hifiasm re-derives read ids by name; the
        // names match those dinara loaded from the same input file(s).
        performanceLog << timestamp << "Generating overlaps with hifiasm from "
            << inputFileNames.size() << " input file(s)." << endl;

        vector<const char*> readFiles;
        readFiles.reserve(inputFileNames.size());
        for(const string& f: inputFileNames) {
            readFiles.push_back(f.c_str());
        }

        hifiasm_ovlp_opt_t hifiOpt = {};
        hifiOpt.threads = int(threadCount);

        char* pafPathC = nullptr;
        const int rc = hifiasm_detect_overlaps(
            readFiles.data(), int(readFiles.size()),
            "hifiasm", &hifiOpt, &pafPathC);
        if(rc != 0 || pafPathC == nullptr) {
            if(pafPathC) {
                free(pafPathC);
            }
            throw runtime_error("hifiasm overlap detection failed (code "
                + to_string(rc) + ").");
        }
        overlapsPafFile = pafPathC;
        free(pafPathC);
        performanceLog << timestamp << "hifiasm overlaps written to "
            << overlapsPafFile << endl;
    }

    if(!overlapsPafFile.empty()) {
        assembler.importAlignmentCandidatesFromPaf(overlapsPafFile, threadCount);
        assembler.chainPafCandidates(
            assemblerOptions.overlapCandidatesOptions.driftRateTolerance,
            maxChainLimit,
            assemblerOptions.overlapCandidatesOptions,
            threadCount
        );
    } else {
        assembler.chainAlignmentCandidates(
            assemblerOptions.overlapCandidatesOptions.driftRateTolerance,
            maxChainLimit,
            assemblerOptions.overlapCandidatesOptions,
            threadCount
        );
    }

    // // Lightweight marker-chain materialization.
    // // The marker graph vertex builder needs alignmentData/compressedAlignments,
    // // but this prototype does not need projected banded/base alignments or evidence.
    // assembler.computeAlignmentDataFromChainedCandidatesOnly(
    //     assemblerOptions.alignOptions,
    //     threadCount);

    // Compute base-level pairwise alignments for all overlaps and store
    // the resulting CIGARs. These are used downstream for CIGAR-based
    // SNP/indel detection in the phasing windows.
    assembler.computeBaseAlignmentsAndStore(
        assemblerOptions.alignOptions,
        threadCount);

    // Build a vector of ReadIds sorted by read length (longest first).
    const Reads& reads = assembler.getReads();
    const ReadId readCount = reads.readCount();
    vector<ReadId> readIdsSortedByLength(readCount);
    iota(readIdsSortedByLength.begin(), readIdsSortedByLength.end(), ReadId(0));
    sort(readIdsSortedByLength.begin(), readIdsSortedByLength.end(),
        [&reads](ReadId a, ReadId b) {
            return reads.getReadRawSequenceLength(a) >
                   reads.getReadRawSequenceLength(b);
        });

    // For http server and debugging/development purposes, generate an exhaustive table of candidates.
    // This can be done after alignment computation (it depends only on the candidate list).
    assembler.computeCandidateTable();


    // assembler.phaseOverlaps(threadCount);
    // assembler.phaseOverlapsKmeans(threadCount);

    // // Keep one best chain per read pair per strand (hifiasm dedup port).
    // // Same-strand and opposite-strand overlaps are deduped independently,
    // // matching hifiasm's separate paf/reverse_paf storage.
    // // Runs after phasing so that only cis and unclassified overlaps compete.
    // assembler.dedupChainsPrePhasing(threadCount);

    // For each read pair, if there are multiple chains on the same strand,
    // delete all of them. A legitimate overlap produces one chain per strand.
    // Multiple chains indicate the reads overlap in a repeat region where
    // the chainer found multiple plausible paths. Keeping any of them risks
    // merging distinct repeat copies during transitive collapse.
    assembler.removeMultiChainAlignments(threadCount);

    // assembler.performHifiasmECParity(threadCount);

    // ---- Post-phasing overlap cleaning (hifiasm order: Overlaps.cpp:39390-39726) ----

    // Initialize per-read valid intervals. With minCoverage=0 this sets
    // every read's interval to [0, readLen). Must run before the next two
    // steps which read validReadIntervals to compute hangs and classify
    // overlaps.
    assembler.filterLocalSegments(/* minCoverage */ 0, threadCount);

    // // Flag chimeric reads — reads whose overlaps don't form a consistent
    // // linear arrangement. A chimeric read has left-side and right-side
    // // overlaps that come from different genomic locations (artifact of
    // // library prep joining unrelated fragments). All overlaps of flagged
    // // reads are deleted. The chimeric flag is also checked by
    // // createMarkerGraphVertices, which skips chimeric reads.
    // assembler.detectChimericReads(threadCount);

    // Delete internal overlaps — overlaps where both reads extend
    // significantly beyond the aligned region on the same side.
    //
    // NOTE (corrected): ad.qs/qe/ts/te are the TIGHT, real CIGAR-alignment
    // span (computeBaseAlignmentsAndStoreThreadFunction), not extended to
    // read tips. AlignmentData::extendedQs/extendedQe/extendedTs/extendedTe
    // are the hifiasm ma_hit_t-convention coordinates (diagonally
    // extrapolated to read boundaries via extendOverlapToReadBoundaries,
    // overlapClassification.hpp) that ma_hit2arc actually needs.
    //
    // Two variants:
    //   deleteInternalOverlaps:         feeds the TIGHT ad.qs/qe/ts/te
    //                                   straight to ma_hit2arc -- NOT hifiasm
    //                                   parity despite the name/old comment here.
    //   deleteInternalOverlapsExtended: correctly uses
    //                                   ad.extendedQs/extendedQe/extendedTs/extendedTe.
    //                                   This is the hifiasm-parity one.
    // assembler.deleteInternalOverlaps(/* maxHang */ 1000, /* maxHangRate */ 0.8, /* minOverlapLength */ 50, threadCount);
    // assembler.deleteInternalOverlapsExtended(/* maxHang */ 1000, /* maxHangRate */ 0.8, /* minOverlapLength */ 50, threadCount);

    // assembler.removeContainedReads(/* maxHang */ 1000, /* maxHangRate */ 0.8, /* minOverlapLength */ 50, threadCount);

    // assembler.rescueTransOverlaps(/* minPileup */ 4, /* skipDeleted */ true);

    // Build the read graph used for marker graph vertex construction.
    // Includes only alignments that:
    // - Are not deleted (no deleteReasons on either side — filters out
    //   chimeric reads, multi-chain pairs, and other earlier removals).
    // - Are not trans (state 2) or cisDifferentCopy (state 3) on either side.
    // Cis (state 1) and unlabeled (state 0) alignments are kept.
    // The read graph edges drive the disjoint set merges in createMarkerGraphVertices.
    assembler.createReadGraphFromPhasingCisOverlaps();

    // // 8. Transitive reduction: remove redundant edges where v→x can be
    // //    reached through v→w→x within fuzz tolerance.
    // //    Port of asg_arc_del_trans.
    // assembler.transitiveReductionOnReadGraph(/* fuzz */ 5000);

    // Set min and max marker graph vertex coverage thresholds.
    // const uint64_t minAnchorCoverage = std::max((uint64_t)3, (uint64_t)(0.15 * double(coverageHet) / 2));
    // const uint64_t maxAnchorCoverage = (uint64_t)(1.5 * double(coverageHet));
    
    const uint64_t minVertexCoverage = assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage;
    const uint64_t maxVertexCoverage = std::numeric_limits<uint64_t>::max();

    // Build marker graph vertices by transitive alignment collapse.
    // Each alignment merges its aligned marker pairs into a disjoint set.
    // Connected components become vertices. Vertices are then filtered by:
    // - Coverage: must have [minVertexCoverage, maxVertexCoverage] markers.
    // - No duplicate reads: a vertex cannot contain two markers from the
    //   same readId (either strand). This prevents collapsing both strands
    //   of a read into the same vertex.
    assembler.createMarkerGraphVertices(
        minVertexCoverage,
        maxVertexCoverage,
        0,                                              // minVertexCoveragePerStrand (disabled)
        false,                                          // allowDuplicateMarkers
        std::numeric_limits<double>::signaling_NaN(),   // unused (minVertexCoverage != 0)
        invalid<uint64_t>,                              // unused (minVertexCoverage != 0)
        threadCount);

    // Repeat-kmer and low-complexity filtering now happens at the minimizer
    // stage (applyKmerCountFilter with filterRepeatKmers / filterLowComplexity),
    // so these marker-graph vertex filters are redundant -- the offending k-mers
    // never seed a vertex. Left commented out; re-enable if the minimizer-stage
    // filters are ever turned off.
    //
    // Remove vertices whose k-mer is a short-period tandem repeat (period 1-5,
    // including homopolymers). Thresholds: {6, 4, 4, 4, 4}. Removes ~18.5%.
    // assembler.filterMarkerGraphVerticesByRepeatKmers(threadCount);
    //
    // Remove vertices whose k-mer has low sequence complexity (distinct
    // sub-k-mers of lengths 1, 2, 3). Thresholds: {4, 12, 24}. Removes ~4.5%.
    // assembler.filterMarkerGraphVerticesByDistinctSubkmerCount(threadCount);

    // Remove vertices created by transitive-collapse false merges, detected via
    // cross-read order agreement on consecutive journey vertices: if one read's
    // journey visits vertex A then B, and some other read visits the same pair
    // in the opposite order, that is a direct, cheap witness that A and/or B
    // cannot both be single, correctly-placed genomic loci. Complementary to,
    // not a subset of, the chain-consistency filter below: that filter only
    // checks read pairs that are direct alignment candidates, so it can miss
    // violations between reads that never became candidates -- on one test
    // region every flagged vertex here was also caught below, on another
    // ~30% were not. Low recall by design (a windowed sweep confirmed it's a
    // narrower, cheaper signal, not a stand-in for the broader filter below).
    // Runs first so the chain-consistency filter below sees a slightly
    // smaller, cheaply-pre-cleaned vertex set.
    assembler.filterMarkerGraphVerticesByJourneyOrderConsistency(threadCount);

    // Remove vertices where the transitive collapse grouped reads at k-mer
    // positions outside their direct chaining range. For each pair of reads
    // in a vertex, check that the vertex ordinal falls within the chain
    // start/end for both reads. Vertices failing this check were created by
    // indirect transitive paths (A→B→C) where A and C have no direct
    // alignment support at that position — typically false merges in repeats.
    assembler.filterMarkerGraphVerticesByChainConsistency(threadCount);

    // Pair each marker graph vertex with its reverse complement vertex.
    // Required before anchor generation, which needs RC-consistent vertices.
    assembler.findMarkerGraphReverseComplementVertices(threadCount);


    const uint64_t minAnchorCoverage = assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage;
    const uint64_t maxAnchorCoverage = std::numeric_limits<uint64_t>::max();

    // const uint64_t minPrimaryCoverage = assemblerOptions.assemblyOptions.mode3Options.minAnchorCoverage;;
    // const uint64_t maxPrimaryCoverage = assemblerOptions.assemblyOptions.mode3Options.maxAnchorCoverage;;
    cout << "Using: minAnchorCoverage = " << minAnchorCoverage <<
        ", maxAnchorCoverage = " << maxAnchorCoverage << endl;


    const MappedMemoryOwner shasta2Owner = assembler.shasta2MappedMemoryOwner();
    
    assembler.shasta2Anchors = make_shared<Shasta2Anchors>(
        shasta2Owner,
        assembler.getReads(),
        assembler.assemblerInfo->k,
        *assembler.markers,
        assembler.markerGraph,
        threadCount,
        minAnchorCoverage,
        maxAnchorCoverage);
        auto& shasta2Anchors = assembler.shasta2Anchors;

    // External-anchor export is deferred until after per-window MSA het-anchor
    // generation (testAbpoaMultiSegmentMSA below), so that any anchors created
    // from detected het sites are included in the exported set. Journeys and
    // windows are built from shasta2Anchors but do not depend on the export,
    // so moving the export down is dependency-safe.
    const string externalAnchorsName =
        std::filesystem::absolute("Shasta2ExternalAnchors").string();

    // Compute journeys.
    cout << timestamp << "Creating Shasta2Journeys..." << endl;
    assembler.shasta2Journeys = make_shared<Shasta2Journeys>(
        2 * assembler.getReads().readCount(),
        shasta2Anchors,
        threadCount,
        shasta2Owner);
    auto& shasta2Journeys = assembler.shasta2Journeys;

    // Filter each read's journey to its longest well-supported anchor chain
    // (every consecutive pair sharing >= minCommonForBackbone reads, bounded
    // look-back maxSkipForBackbone). Runs per read independently and rewrites
    // the stored journeys + positionInJourney before any windowing decision, so
    // every downstream stage sees the cleaned chains.
    cout << timestamp << "Filtering journeys by anchor chaining..." << endl;
    shasta2Journeys->filterByAnchorChaining(
        assemblerOptions.assemblyOptions.mode3Options.minCommonForBackbone,
        assemblerOptions.assemblyOptions.mode3Options.maxSkipForBackbone,
        threadCount);

    // MSA-based overlap phasing — disabled, replaced by CIGAR-based window pipeline.
    // assembler.phaseOverlapsMSA(threadCount);

    // Flag contained reads so they can be excluded from inter-window edge discovery.
    cout << timestamp << "Flagging contained reads..." << endl;
    assembler.flagContainedReads(1000, 0.8, 0, threadCount);

    // Compute anchor windows.
    cout << timestamp << "Computing anchor windows..." << endl;
    // ========================================================================
    // PHASE 1: Window creation.
    // Partitions anchors into disjoint windows. A first pass seeds pristine
    // full-journey cores (one window = one whole read journey); a second pass
    // then tiles the leftover unclaimed base spans into fragment windows.
    // Produces `anchorWindows` only — no inter-window edges, no transitions.
    // ========================================================================
    vector<AnchorWindow> anchorWindows;
    vector<uint32_t> anchorDovetailWindow;
    const uint64_t minCommonForBackbone =
        assemblerOptions.assemblyOptions.mode3Options.minCommonForBackbone;
    const uint64_t maxSkipForBackbone =
        assemblerOptions.assemblyOptions.mode3Options.maxSkipForBackbone;
    // Two-pass window creation. Pass 1 seeds a window only from a read whose
    // entire journey is still unclaimed (pristine disjoint full-journey core),
    // claiming its anchors so windows never overlap. Pass 2
    // (tileUnclaimedIntervals=true) then tiles the leftover unclaimed runs into
    // fragment windows (longest base span first), so reads whose journey overlaps
    // an already-claimed core still contribute windows from their remaining
    // contiguous unclaimed intervals rather than being dropped.
    assembler.computeAnchorWindowsClean(
        assembler.shasta2Anchors,
        assembler.shasta2Journeys,
        readIdsSortedByLength,
        anchorWindows,
        threadCount,
        minCommonForBackbone,
        maxSkipForBackbone,
        assemblerOptions.assemblyOptions.mode3Options.minWindowBaseSpan,
        &anchorDovetailWindow,
        /* tileUnclaimedIntervals = */ true);

    // High-connectivity het gate: optionally suppress het detection in windows
    // that sit at tangles/repeats -- i.e. windows with many distinct incoming
    // AND many distinct outgoing inter-window neighbors. Per-window het calls in
    // such windows are unreliable, so we leave them homozygous.
    //
    // In/out degree is derived from the same per-window transition map the
    // anchor-graph constructor consumes. computeWindowTransitions clears and
    // repopulates its output fields on entry, so calling it here (before het
    // detection) and again later (before the graph build) is safe and
    // idempotent for the same journeys/windows: it reads only backbone
    // journeys and each window's own readIntervals/backboneBegin/End/
    // filteredBackbonePositions, none of which the intervening het-bubble
    // plan/append/stage passes touch (they mutate hetBubbles/intraWindowEdges
    // only). So the two calls would produce byte-identical output -- windowTransitionsComputed
    // tracks whether the gate already paid for it, so the graph-build call
    // below can skip a second full recomputation.
    std::vector<bool> hetSkipWindow;
    bool windowTransitionsComputed = false;
    {
        const uint64_t hetMaxInDeg =
            assemblerOptions.assemblyOptions.mode3Options.hetMaxWindowInDegree;
        const uint64_t hetMaxOutDeg =
            assemblerOptions.assemblyOptions.mode3Options.hetMaxWindowOutDegree;
        // Gate is active only when BOTH thresholds are set (>0); either at 0
        // disables it (0 would otherwise match every window trivially).
        if(hetMaxInDeg > 0 && hetMaxOutDeg > 0) {
            windowTransitionsComputed = true;
            computeWindowTransitions(*shasta2Anchors, *shasta2Journeys,
                anchorWindows, &anchorDovetailWindow);

            hetSkipWindow.assign(anchorWindows.size(), false);
            const uint32_t noW = AnchorWindowReadInterval::noWindow;
            uint64_t gatedWindows = 0;
            for(uint64_t w = 0; w < anchorWindows.size(); w++) {
                // Count distinct incoming and outgoing neighbor windows from the
                // transition keys (previousWindow, nextWindow); noWindow (read
                // starts/ends here) is not a real neighbor.
                std::set<uint32_t> inNbrs;
                std::set<uint32_t> outNbrs;
                for(const auto& kv : anchorWindows[w].transitionReads) {
                    const uint32_t prev = kv.first.first;
                    const uint32_t next = kv.first.second;
                    if(prev != noW) inNbrs.insert(prev);
                    if(next != noW) outNbrs.insert(next);
                }
                if(inNbrs.size() >= hetMaxInDeg && outNbrs.size() >= hetMaxOutDeg) {
                    hetSkipWindow[w] = true;
                    gatedWindows++;
                }
            }
            cout << timestamp << "High-connectivity het gate active"
                 << " (minInDegree=" << hetMaxInDeg
                 << ", minOutDegree=" << hetMaxOutDeg
                 << "): gating " << gatedWindows << " of "
                 << anchorWindows.size() << " windows out of het detection."
                 << endl;
        }
    }

    // Per-window het-bubble detection. Interchangeable engines produce the
    // SAME output (AnchorWindow::hetBubbles), so the downstream plan/append/stage
    // passes are identical regardless of which one runs:
    //   - abpoa (default): all-reads multi-segment abPOA MSA per window, bubbles
    //     read off the partial-order graph.
    //   - ksw2: banded per-member ksw2 pileup against the backbone, bubbles
    //     synthesized from the member SNP pileup (ksw2DetectHetBubblesInWindow).
    // Engine selection via env DINARA_HET_ENGINE:
    //   unset / "intervalpoa" (DEFAULT) : per-interval parallel POA. ~5x faster
    //         than whole-window abPOA with equal-or-better het recall; each
    //         window's anchor intervals are POA'd independently and in parallel.
    //   "abpoa"   : legacy whole-window multi-segment abPOA (one growing graph
    //         per window, one window per thread). Capped by DINARA_MSA_MAX_WINDOWS
    //         (default 1, 0=all).
    //   "ksw2"    : star alignment against the backbone (no read-to-read POA).
    //   "projaln" : per-member ProjectedAlignment (computeBaseAlignmentsAndStore's
    //         own fast aligner) sparse-evidence pinning -- see
    //         projAlnDetectHetBubblesAllWindows, AssemblerWindowProjectedAlignmentLeafSnarls.cpp.
    {
        const char* hetEngineEnv = getenv("DINARA_HET_ENGINE");
        const string hetEngine = (hetEngineEnv != nullptr) ? string(hetEngineEnv) : "";
        const bool useKsw2HetEngine = (hetEngine == "ksw2");
        const bool useAbpoaHetEngine = (hetEngine == "abpoa");
        const bool useProjAlnHetEngine = (hetEngine == "projaln");
        // Default engine is per-interval POA: selected explicitly or when unset
        // and no other engine was requested.
        const bool useIntervalPoaHetEngine =
            (hetEngine == "intervalpoa") ||
            (!useKsw2HetEngine && !useAbpoaHetEngine && !useProjAlnHetEngine);

        if(useIntervalPoaHetEngine) {
            cout << timestamp << "Detecting het bubbles with the per-interval POA"
                    " engine on " << anchorWindows.size() << " windows on "
                 << threadCount << " threads..." << endl;
            const auto tHet0 = steady_clock::now();
            const double hetMinVaf = assemblerOptions.assemblyOptions.mode3Options.hetMinVaf;
            const uint64_t hetMinSupport = assemblerOptions.assemblyOptions.mode3Options.hetMinSupport;
            const bool hetDropHomopolymer = assemblerOptions.assemblyOptions.mode3Options.hetDropHomopolymer;
            const bool hetDropRepeat = assemblerOptions.assemblyOptions.mode3Options.hetDropRepeat;
            // Global interval load balancing (shasta2 assembleChainsMultithreaded
            // model): every interval of every window is one work unit, flattened
            // into a single list, sorted biggest-first, and run batch=1 across
            // threads. This keeps all threads busy even when a few large windows
            // hold most intervals -- the previous per-window scheduling left
            // threads idle once the small windows finished while one thread
            // ground through a window with thousands of serial intervals.
            uint64_t hetWindows = 0;
            uint64_t totalBubbles = 0;
            assembler.intervalPoaDetectHetBubblesAllWindows(
                anchorWindows, *shasta2Anchors, *shasta2Journeys,
                hetMinVaf, hetMinSupport, hetDropHomopolymer, hetDropRepeat,
                threadCount, hetWindows, totalBubbles,
                hetSkipWindow.empty() ? nullptr : &hetSkipWindow);
            const double hetSecs = seconds(steady_clock::now() - tHet0);
            dinaraPrintHetEmitTiming();
            // Verification-only: run the whole-window abpoa graph (one graph
            // per seed/backbone read, members added piece-by-piece via their
            // own shared anchors) and its leaf-snarl detector alongside the
            // live per-interval pipeline, logging results without feeding
            // them into window.hetBubbles or anything else live. Off by
            // default; opt in with DINARA_WINDOW_ABPOA_DEBUG=1.
            if (getenv("DINARA_WINDOW_ABPOA_DEBUG") != nullptr) {
                assembler.computeWindowAbpoaGraphs(
                    anchorWindows, *shasta2Anchors, *shasta2Journeys,
                    "windowAbpoaDebug.", threadCount);
            }
            // Verification-only: leaf-snarl detection from independent
            // pairwise ksw2 alignments (no shared graph, no per-interval
            // fragmentation) -- a comparison point against both the abPOA
            // shared-graph version above and the interval-abPOA-profile
            // version (DINARA_INTERVALPOA_LEAFSNARL_DEBUG=1, inside
            // intervalPoaDetectHetBubblesAllWindows itself). Off by default;
            // opt in with DINARA_KSW2_LEAFSNARL_DEBUG=1.
            if (getenv("DINARA_KSW2_LEAFSNARL_DEBUG") != nullptr) {
                assembler.computeWindowKsw2LeafSnarls(
                    anchorWindows, *shasta2Anchors, *shasta2Journeys,
                    assemblerOptions.alignOptions, threadCount);
            }
            // Verification-only: leaf-snarl detection using ProjectedAlignment
            // (computeBaseAlignmentsAndStore's own fast aligner) instead of
            // ksw2's per-segment banded DP. Off by default; opt in with
            // DINARA_PROJALN_LEAFSNARL_DEBUG=1.
            if (getenv("DINARA_PROJALN_LEAFSNARL_DEBUG") != nullptr) {
                assembler.computeWindowProjectedAlignmentLeafSnarls(
                    anchorWindows, *shasta2Anchors, *shasta2Journeys,
                    assemblerOptions.alignOptions, threadCount);
            }
            cout << timestamp << "per-interval POA het-bubble detection complete."
                 << " hetWindows=" << hetWindows
                 << " homWindows=" << (anchorWindows.size() - hetWindows)
                 << " totalHetBubbles=" << totalBubbles
                 << " seconds=" << std::fixed << std::setprecision(2) << hetSecs
                 << std::defaultfloat << endl;
        } else if(useKsw2HetEngine) {
            cout << timestamp << "Detecting het bubbles with the ksw2 pileup"
                    " engine on " << anchorWindows.size() << " windows..." << endl;
            const auto tHet0 = steady_clock::now();
            uint64_t hetWindows = 0, totalBubbles = 0;
            // CIGAR-density noise filter window/threshold (HiFi defaults 100/5).
            constexpr int noisyRegSlideWin = 100;
            constexpr int noisyRegMaxXgaps = 5;
            for(uint64_t wi = 0; wi < anchorWindows.size(); wi++) {
                if(!hetSkipWindow.empty() && hetSkipWindow[wi]) continue;
                AnchorWindow& window = anchorWindows[wi];
                const uint32_t n = assembler.ksw2DetectHetBubblesInWindow(
                    window, *shasta2Anchors, *shasta2Journeys,
                    assemblerOptions.alignOptions,
                    assemblerOptions.assemblyOptions.mode3Options.hetMinVaf,
                    assemblerOptions.assemblyOptions.mode3Options.hetMinSupport,
                    assemblerOptions.assemblyOptions.mode3Options.hetDropHomopolymer,
                    assemblerOptions.assemblyOptions.mode3Options.hetDropRepeat,
                    noisyRegSlideWin, noisyRegMaxXgaps);
                if(n > 0) { hetWindows++; totalBubbles += n; }
            }
            const double hetSecs = seconds(steady_clock::now() - tHet0);
            cout << timestamp << "ksw2 het-bubble detection complete."
                 << " hetWindows=" << hetWindows
                 << " homWindows=" << (anchorWindows.size() - hetWindows)
                 << " totalHetBubbles=" << totalBubbles
                 << " seconds=" << std::fixed << std::setprecision(2) << hetSecs
                 << std::defaultfloat << endl;
            dinaraPrintHetEmitTiming();
            { extern void dinaraPrintKswCallCount(); dinaraPrintKswCallCount(); }
        } else if(useProjAlnHetEngine) {
            cout << timestamp << "Detecting het bubbles with the ProjectedAlignment"
                    " engine on " << anchorWindows.size() << " windows on "
                 << threadCount << " threads..." << endl;
            const auto tHet0 = steady_clock::now();
            uint64_t hetWindows = 0, totalBubbles = 0;
            assembler.projAlnDetectHetBubblesAllWindows(
                anchorWindows, *shasta2Anchors, *shasta2Journeys,
                assemblerOptions.alignOptions,
                assemblerOptions.assemblyOptions.mode3Options.hetMinVaf,
                assemblerOptions.assemblyOptions.mode3Options.hetMinSupport,
                assemblerOptions.assemblyOptions.mode3Options.hetDropHomopolymer,
                assemblerOptions.assemblyOptions.mode3Options.hetDropRepeat,
                threadCount, hetWindows, totalBubbles,
                hetSkipWindow.empty() ? nullptr : &hetSkipWindow);
            const double hetSecs = seconds(steady_clock::now() - tHet0);
            cout << timestamp << "ProjectedAlignment het-bubble detection complete."
                 << " hetWindows=" << hetWindows
                 << " homWindows=" << (anchorWindows.size() - hetWindows)
                 << " totalHetBubbles=" << totalBubbles
                 << " seconds=" << std::fixed << std::setprecision(2) << hetSecs
                 << std::defaultfloat << endl;
        } else {  // useAbpoaHetEngine (DINARA_HET_ENGINE=abpoa)
            const auto tAbpoa0 = steady_clock::now();
            assembler.testAbpoaMultiSegmentMSA(
                assembler.shasta2Anchors,
                assembler.shasta2Journeys,
                anchorWindows,
                threadCount,
                assemblerOptions.assemblyOptions.mode3Options.hetMinVaf,
                assemblerOptions.assemblyOptions.mode3Options.hetMinSupport,
                assemblerOptions.assemblyOptions.mode3Options.hetDropHomopolymer,
                assemblerOptions.assemblyOptions.mode3Options.hetDropRepeat);
            const double abpoaSecs = seconds(steady_clock::now() - tAbpoa0);
            cout << timestamp << "abPOA het-bubble detection complete."
                 << " windows=" << anchorWindows.size()
                 << " threads=" << threadCount
                 << " seconds=" << std::fixed << std::setprecision(2) << abpoaSecs
                 << std::defaultfloat << endl;
        }
    }

    // Verification-only timing comparison: cigarDetectSnpsInWindow reuses
    // computeBaseAlignmentsAndStore's ALREADY-COMPUTED pairwise CIGARs
    // (overlapCigarStore) instead of doing fresh per-window realignment (ksw2
    // or interval-abPOA above) -- a comparison point for whether that
    // already-paid alignment cost, reused directly, beats redoing alignment
    // per window. Writes only window.cleanHetSnpCount/hetSnps (not
    // hetBubbles), so it does not disturb whichever engine ran above. Off by
    // default; opt in with DINARA_CIGAR_SNP_DEBUG=1.
    if (getenv("DINARA_CIGAR_SNP_DEBUG") != nullptr) {
        cout << timestamp << "Running CIGAR-reuse SNP detection on "
             << anchorWindows.size() << " windows (comparison only)..." << endl;
        const auto tCigar0 = steady_clock::now();
        uint64_t totalCigarSnps = 0;
        for(AnchorWindow& window : anchorWindows) {
            window.cleanHetSnpCount = assembler.cigarDetectSnpsInWindow(
                window, *shasta2Anchors, *shasta2Journeys);
            totalCigarSnps += window.cleanHetSnpCount;
        }
        const double cigarSecs = seconds(steady_clock::now() - tCigar0);
        cout << timestamp << "CIGAR-reuse SNP detection complete."
             << " totalSnps=" << totalCigarSnps
             << " seconds=" << std::fixed << std::setprecision(2) << cigarSecs
             << std::defaultfloat << endl;
    }

    // Turn the staged het bubbles into anchor-graph structure in three serial
    // passes (plan -> append -> stage edges). The passes are ordered so anchors
    // are created only for bubbles that will be wired: the planner marks which
    // bubbles are usable, the append pass creates exactly those anchors, and the
    // staging pass wires them. See the helpers in the anonymous namespace above.
    const uint32_t hetKHalf = uint32_t(shasta2Anchors->k / 2);

    // Backbone cache: computed once per window here, then shared read-only by
    // Pass 1, Pass 1.5, and Pass 3 below (see WindowBackboneCache for why this
    // is safe -- nothing between window creation and Pass 3 changes a
    // window's backbone fields or moves an existing anchor's position).
    cout << timestamp << "Computing window backbones..." << endl;
    const auto tBackboneCache0 = steady_clock::now();
    const WindowBackboneCache backboneCache(
        *shasta2Anchors, *shasta2Journeys, anchorWindows, hetKHalf, threadCount);
    cout << timestamp << "Window backbones computed in "
         << seconds(steady_clock::now() - tBackboneCache0) << "s." << endl;

    // Pass 1: plan. Assign each bubble to the backbone interval that strictly
    // contains its flank span; drop the rest (plannedInterval = -1). Also drop
    // any bubble whose bracketing homs would collide with an existing primary
    // anchor's (read, position) marker (see planWindowHetBubbles).
    uint64_t plannedBubbles = 0, droppedUncontained = 0, droppedPrimaryCollision = 0;
    uint64_t droppedBackwardMembers = 0;
    {
        // (read, position) markers owned by primary anchors, built once before
        // any het anchor is appended (so the store still holds only primaries).
        const auto tPrimarySet0 = steady_clock::now();
        const PrimaryMarkerSet primarySet = buildPrimaryMarkerSet(*shasta2Anchors);
        cout << timestamp << "buildPrimaryMarkerSet: " << primarySet.size()
             << " markers in " << seconds(steady_clock::now() - tPrimarySet0) << "s." << endl;
        const auto tPlanLoop0 = steady_clock::now();
        // Parallel over windows: planning is window-local (mutates only its own
        // window's bubbles) and reads shared state read-only (the store via
        // getPosition, the frozen primarySet, the frozen backboneCache).
        std::atomic<uint64_t> aPlanned{0}, aUncontained{0}, aPrimaryCollision{0};
        std::atomic<uint64_t> aBackwardMembers{0};
        parallelForEachWindow(anchorWindows.size(), threadCount, [&](uint64_t wi) {
            AnchorWindow& window = anchorWindows[wi];
            if(!backboneCache.hasBackbone(wi)) {
                // Fewer than two backbone anchors: no interval, drop all bubbles.
                for(auto& b : window.hetBubbles) { b.plannedInterval = -1; }
                aUncontained.fetch_add(window.hetBubbles.size(),
                                       std::memory_order_relaxed);
                return;
            }
            // NOTE: het members are stored at rawPosition + hetAnchorKHalf()
            // (1 for k=2, 0 for k=0; see appendHetAnchorPair), NOT the store's
            // k/2. Pass the HET half-length, not the k/2 used for the backbone
            // frame, so the backward-member comparison uses the same frame the
            // verifier sees.
            uint64_t planned = 0, uncontained = 0, primaryCollision = 0, backward = 0;
            planWindowHetBubbles(window, *shasta2Anchors,
                                 backboneCache.bbAnchors(wi), backboneCache.bbOffset(wi),
                                 /*hetKHalf=*/hetAnchorKHalf(), primarySet,
                                 planned, uncontained,
                                 primaryCollision, backward);
            aPlanned.fetch_add(planned, std::memory_order_relaxed);
            aUncontained.fetch_add(uncontained, std::memory_order_relaxed);
            aPrimaryCollision.fetch_add(primaryCollision, std::memory_order_relaxed);
            aBackwardMembers.fetch_add(backward, std::memory_order_relaxed);
        });
        plannedBubbles = aPlanned.load();
        droppedUncontained = aUncontained.load();
        droppedPrimaryCollision = aPrimaryCollision.load();
        droppedBackwardMembers = aBackwardMembers.load();
        cout << timestamp << "Planning loop: " << seconds(steady_clock::now() - tPlanLoop0) << "s." << endl;
        cout << timestamp << "Planned " << plannedBubbles
             << " contained het bubbles (" << droppedUncontained
             << " dropped, of which " << droppedPrimaryCollision
             << " for hom/primary anchor collision, "
             << droppedBackwardMembers << " backward hom members removed)." << endl;
    }

    // Pass 1.5: merge coincident hom anchors. When two adjacent bubbles in the
    // same interval have succBackboneOffset_N == predBackboneOffset_{N+1} (their
    // SNPs are exactly 3 bp apart), the trailing hom of bubble N and the leading
    // hom of bubble N+1 are the SAME abPOA node -- same members, same recovered
    // read positions. Appending them as two separate anchor ids would build two
    // anchors on one node and stage a hom_N -> leadHom_{N+1} edge that is
    // equal-position on every shared read (the crash step 3 catches). Instead we
    // fold bubble N+1's leading hom onto bubble N's trailing hom so the chain
    // becomes ...arms_N -> sharedHom -> arms_{N+1}..., keeping BOTH SNPs.
    {
        const auto tMergeLoop0 = steady_clock::now();
        // Parallel over windows: the merge is window-local (folds one window's
        // adjacent-bubble homs) and reads no mutable shared state (the
        // backboneCache is frozen after its construction above).
        std::atomic<uint64_t> aMergedHoms{0};
        parallelForEachWindow(anchorWindows.size(), threadCount, [&](uint64_t wi) {
            if(!backboneCache.hasBackbone(wi)) {
                return;
            }
            AnchorWindow& window = anchorWindows[wi];
            uint64_t merged = 0;
            mergeWindowCoincidentHoms(window, merged);
            aMergedHoms.fetch_add(merged, std::memory_order_relaxed);
        });
        cout << timestamp << "Merge loop: " << seconds(steady_clock::now() - tMergeLoop0) << "s." << endl;
        cout << timestamp << "Merged " << aMergedHoms.load()
             << " coincident hom anchors (adjacent SNPs 3 bp apart)." << endl;
    }

    // Pass 2: append. Create the canonical/RC anchor pairs for planned bubbles.
    // This grows the memory-mapped store and invalidates outstanding anchor
    // spans, so it must run after all window processing.
    {
        cout << timestamp << "Appending het anchors from "
             << anchorWindows.size() << " windows..." << endl;
        uint64_t hetAnchorPairs = 0, homAnchorPairs = 0;
        for(AnchorWindow& window : anchorWindows)
            appendWindowHetAnchors(*shasta2Anchors, window,
                                   hetAnchorPairs, homAnchorPairs);
        cout << timestamp << "  (" << homAnchorPairs
             << " hom separator anchor pairs)" << endl;
        cout << timestamp << "Appended " << hetAnchorPairs
             << " het anchor pairs (" << (2 * hetAnchorPairs)
             << " anchors); store now has " << shasta2Anchors->size()
             << " anchors." << endl;
    }

    // Pass 3: stage window-local (intra-window) anchor-graph edges on the
    // windows: the backbone chain (consecutive backbone anchors + RC mirror) and
    // the het bubbles (pred->allele->succ + RC mirror). These edges are fully
    // determined by a single window, so they live in the window structure; the
    // anchor graph constructor replays window.intraWindowEdges directly.
    // Inter-window edges are NOT staged here (they need multiple windows and are
    // discovered in a global pass). Runs after the het-anchor append so allele
    // anchor ids are assigned.
    {
        const uint64_t anchorCount = shasta2Anchors->size();
        // Parallel over windows: staging is window-local (mutates only its own
        // window.intraWindowEdges) and reads the store read-only. anchorCount is
        // fixed here (append is complete), so all endpoint ids are valid. The
        // backboneCache is frozen and unaffected by the append (see its comment).
        std::atomic<uint64_t> aBackbone{0}, aHet{0}, aChained{0};
        parallelForEachWindow(anchorWindows.size(), threadCount, [&](uint64_t wi) {
            if(!backboneCache.hasBackbone(wi)) {
                return;
            }
            AnchorWindow& window = anchorWindows[wi];
            const vector<Shasta2AnchorId>& bbAnchors = backboneCache.bbAnchors(wi);
            const vector<uint32_t>& bbOffset = backboneCache.bbOffset(wi);
            thread_local vector<uint32_t> bbExportedOffset;
            // Step 1: middle-2 backbone shift. Materialize the exported (k=2)
            // frame the edges and monotonicity check operate in. The shift is a
            // uniform -1 per anchor, so it preserves the backbone's strict
            // ordering; assert that here so a broken window backbone is caught
            // before it reaches the edge staging.
            computeWindowShiftedBackbone(bbOffset, bbExportedOffset);
            for(size_t i = 1; i < bbExportedOffset.size(); i++)
                DINARA_ASSERT(bbExportedOffset[i] > bbExportedOffset[i - 1]);
            uint64_t backbone = 0, het = 0, chained = 0;
            stageWindowIntraEdges(window, anchorCount, bbAnchors, bbOffset,
                                  backbone, het, chained);
            aBackbone.fetch_add(backbone, std::memory_order_relaxed);
            aHet.fetch_add(het, std::memory_order_relaxed);
            aChained.fetch_add(chained, std::memory_order_relaxed);
        });
        const uint64_t stagedBackbone = aBackbone.load();
        const uint64_t stagedHet = aHet.load();
        const uint64_t chainedIntervals = aChained.load();
        cout << timestamp << "Staged intra-window edges: "
             << stagedBackbone << " backbone, " << stagedHet << " het ("
             << chainedIntervals << " chained intervals)." << endl;
    }

    // Step 3: monotonicity verification. Every staged intra-window edge must be
    // strictly forward on the reads its endpoints share (exported position at B
    // > exported position at A). Runs after staging (edges present) and after
    // the het-anchor append (all endpoint ids valid). A backward edge is a bug
    // in planning/staging, so the verifier throws.
    {
        // Parallel over windows: both verifiers are read-only. They THROW on a
        // violation; parallelForEachWindow captures the first exception and
        // rethrows it on this thread after the workers join, so a violation
        // still aborts the run with its original message.
        std::atomic<uint64_t> aCheckedEdges{0}, aCheckedReadSteps{0};
        parallelForEachWindow(anchorWindows.size(), threadCount, [&](uint64_t wi) {
            uint64_t edges = 0, steps = 0;
            verifyWindowEdgeMonotonicity(*shasta2Anchors, anchorWindows[wi],
                                         edges, steps);
            aCheckedEdges.fetch_add(edges, std::memory_order_relaxed);
            aCheckedReadSteps.fetch_add(steps, std::memory_order_relaxed);
        });
        cout << timestamp << "Verified " << aCheckedEdges.load()
             << " intra-window edges forward-monotonic across "
             << aCheckedReadSteps.load() << " shared-read steps." << endl;

        // Structural check: each window's intra-window graph must be a single
        // linear path of backbone anchors and bubbles, and every het bubble arm
        // must be flanked by hom anchors on both sides.
        const Shasta2AnchorId hetFirst = shasta2Anchors->hetAnchorFirstId;
        std::atomic<uint64_t> aStructWindows{0}, aStructBubbles{0};
        parallelForEachWindow(anchorWindows.size(), threadCount, [&](uint64_t wi) {
            uint64_t windows = 0, bubbles = 0;
            verifyWindowGraphStructure(anchorWindows[wi], hetFirst,
                                       windows, bubbles);
            aStructWindows.fetch_add(windows, std::memory_order_relaxed);
            aStructBubbles.fetch_add(bubbles, std::memory_order_relaxed);
        });
        cout << timestamp << "Verified intra-window graph structure: "
             << aStructWindows.load() << " windows linear, "
             << aStructBubbles.load() << " het bubbles hom-flanked." << endl;
    }

    // Global per-read journey tie resolution. shasta2 builds, for every oriented
    // read, the ordered list of ALL anchors that read belongs to (sorted by the
    // read's exported position) and asserts they are STRICTLY INCREASING; two
    // anchors at the SAME exported position on one read triggers "Invalid
    // Journey ...". dinara's own guards are all LOCAL (per intra-window edge, per
    // exported edge, or het/hom-vs-primary only); none enforce this GLOBAL
    // per-read constraint across all het/hom anchors, which the k=0 membership
    // relaxation (pinnedPointCol drops the k=2 flank-adjacency guard) can now
    // violate: two het/hom anchors in different windows can pin the same read at
    // the same exported base position.
    //
    // Model. The exportable unit is (canonicalAnchorId, readId): shasta2 loads
    // only the canonical (even) anchors and REGENERATES the RC of each, so one
    // unit becomes TWO journey occurrences -- one on readId-0, one on readId-1.
    // dinara's store already holds each RC twin as the odd anchor
    // (appendHetAnchorPair), so iterating the FULL store (even + odd) reproduces
    // exactly the occurrence set shasta2 sees. An occurrence's unit is
    // (anchorId & ~1, orientedReadId.getReadId()); dropping that unit removes
    // BOTH the direct and the RC-induced occurrence for the read.
    //
    // Resolution. For each oriented read, group occurrences by exported position;
    // for every position with >1 unit keep one (primary over het/hom, then higher
    // coverage, then lower canonical id) and record the rest in the drop map.
    // Dropping only REMOVES occurrences, so a position group can only shrink
    // (n -> <=1) and no new tie can appear; a single pass therefore leaves every
    // (read, position) group with <=1 survivor regardless of mirror-strand
    // interactions. writeExternalAnchors omits the dropped members.
    Shasta2Anchors::ExternalAnchorDropMap journeyTieDropMap;
    {
        const uint32_t exportShift = hetAnchorKHalf();
        const Shasta2AnchorId hetFirst = shasta2Anchors->hetAnchorFirstId;
        auto isHet = [&](Shasta2AnchorId canonicalId) -> bool {
            return hetFirst != invalid<Shasta2AnchorId> && canonicalId >= hetFirst;
        };
        auto classOf = [&](Shasta2AnchorId canonicalId) -> const char* {
            return isHet(canonicalId) ? "het/hom" : "primary";
        };
        // Per oriented read: (exportedPosition, canonicalAnchorId). The unit's
        // readId is the read's own ReadId, identical for both strands.
        std::unordered_map<uint64_t, vector<std::pair<uint32_t, Shasta2AnchorId>>> byRead;
        const uint64_t anchorCount = shasta2Anchors->size();
        for(Shasta2AnchorId id = 0; id < anchorCount; id++) {
            const Shasta2AnchorId canonicalId = id & ~Shasta2AnchorId(1);
            const Shasta2Anchor anchor = (*shasta2Anchors)[id];
            for(const Shasta2AnchorMarkerInfo& mi : anchor) {
                byRead[mi.orientedReadId.getValue()].push_back(
                    {mi.position - exportShift, canonicalId});
            }
        }
        // Coverage of a canonical anchor (member count), used to pick the keeper.
        auto coverageOf = [&](Shasta2AnchorId canonicalId) -> uint64_t {
            return (*shasta2Anchors)[canonicalId].size();
        };
        // Prefer to KEEP: primary over het/hom, then higher coverage, then lower
        // canonical id. Returns true if a should be kept over b.
        auto keepAOverB = [&](Shasta2AnchorId a, Shasta2AnchorId b) -> bool {
            const bool aHet = isHet(a), bHet = isHet(b);
            if(aHet != bHet) return !aHet;                 // primary wins
            const uint64_t ca = coverageOf(a), cb = coverageOf(b);
            if(ca != cb) return ca > cb;                    // higher coverage wins
            return a < b;                                   // lower id wins
        };
        // Record a dropped unit (canonicalId, readId), de-duplicated.
        auto recordDrop = [&](Shasta2AnchorId canonicalId, ReadId readId) {
            auto& v = journeyTieDropMap[canonicalId];
            if(std::find(v.begin(), v.end(), readId) == v.end()) {
                v.push_back(readId);
            }
        };
        uint64_t tieGroups = 0, unitsDropped = 0, readsWithTie = 0;
        const uint64_t maxReport = 20;
        uint64_t reported = 0;
        for(auto& [oidValue, occ] : byRead) {
            std::sort(occ.begin(), occ.end());
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            const ReadId readId = oid.getReadId();
            bool readCounted = false;
            for(size_t i = 0; i < occ.size(); ) {
                // Advance over a run of equal positions [i, j).
                size_t j = i + 1;
                while(j < occ.size() && occ[j].first == occ[i].first) ++j;
                if(j - i > 1) {
                    ++tieGroups;
                    if(!readCounted) { ++readsWithTie; readCounted = true; }
                    // Choose the keeper among the tied canonical anchors.
                    Shasta2AnchorId keeper = occ[i].second;
                    for(size_t t = i + 1; t < j; t++) {
                        if(keepAOverB(occ[t].second, keeper)) keeper = occ[t].second;
                    }
                    // Drop every tied unit except the keeper.
                    for(size_t t = i; t < j; t++) {
                        const Shasta2AnchorId cId = occ[t].second;
                        if(cId == keeper) continue;
                        recordDrop(cId, readId);
                        ++unitsDropped;
                    }
                    if(reported < maxReport) {
                        ++reported;
                        cout << "  [journey-tie] read " << oid
                             << " pos " << occ[i].first << " keep anchor "
                             << keeper << " (" << classOf(keeper) << "), drop";
                        for(size_t t = i; t < j; t++) {
                            if(occ[t].second == keeper) continue;
                            cout << " " << occ[t].second
                                 << " (" << classOf(occ[t].second) << ")";
                        }
                        cout << endl;
                    }
                }
                i = j;
            }
        }
        cout << timestamp << "Journey tie resolution: " << tieGroups
             << " tied position group(s) on " << readsWithTie
             << " read(s); dropping " << unitsDropped
             << " unit(s) across " << journeyTieDropMap.size()
             << " anchor(s) (--k " << hetAnchorK() << ")." << endl;
    }

    // Write external anchors. Deferred to here (after MSA het-anchor
    // generation) so newly generated het anchors are part of the exported set.
    // The drop map removes the members that would otherwise collide in a
    // per-read journey (see the tie resolution above).
    cout << timestamp << "Writing Shasta2 external anchors to "
         << externalAnchorsName << "..." << endl;
    const uint64_t exportedExternalAnchorCount =
        shasta2Anchors->writeExternalAnchors(
            externalAnchorsName, true, &journeyTieDropMap);
    cout << timestamp << "Wrote " << exportedExternalAnchorCount
         << " external anchors for Shasta2. Use --external-anchors-name "
         << externalAnchorsName << endl;
    // The export subtracts hetAnchorKHalf() uniformly from every stored midpoint
    // (see writeExternalAnchors), so shasta2 must be loaded with the MATCHING
    // --k: 2 by default, 0 for the experimental DINARA_HET_K=0 path. A mismatch
    // shifts every anchor by one base. Report it so the caller passes the right
    // value to the downstream shasta2 invocation.
    cout << timestamp << "Shasta2 must load these external anchors with --k "
         << hetAnchorK()
         << (hetAnchorK() == 0 ? " (EXPERIMENTAL DINARA_HET_K=0)." : ".") << endl;

    // Build and export the anchor graph, including het-anchor bubble edges.
    // computeWindowTransitions fills the per-window transition fields the graph
    // constructor consumes; the constructor then wires backbone chains,
    // inter-window edges, and het bubbles. anchorDovetailWindow was populated
    // above (with the windows). This is the second end product alongside the
    // exported anchors. Skip if the high-connectivity het gate already
    // computed it above -- see the note there for why the result would be
    // identical either way.
    {
        if(!windowTransitionsComputed) {
            computeWindowTransitions(*shasta2Anchors, *shasta2Journeys, anchorWindows,
                &anchorDovetailWindow);
        }

        const uint64_t minInterWindowCoverage =
            assemblerOptions.assemblyOptions.mode3Options.minInterWindowCoverage;
        const uint64_t minInterWindowEdgeCoverage =
            assemblerOptions.assemblyOptions.mode3Options.minInterWindowEdgeCoverage;

        cout << timestamp << "Creating Shasta2AnchorGraph (with het bubbles) from "
             << anchorWindows.size() << " anchor windows..." << endl;
        assembler.shasta2AnchorGraph = make_shared<Shasta2AnchorGraph>(
            *shasta2Anchors,
            *shasta2Journeys,
            anchorWindows,
            minInterWindowCoverage,
            minInterWindowEdgeCoverage,
            threadCount,
            &assembler.getReads(),
            nullptr, // bypassEdges
            nullptr, // detourWindowPairs
            &anchorDovetailWindow);

        // Trim dangling backbone ends that extend past the outermost inter-
        // window connection points. Now that inter-window edges are created,
        // each window's backbone typically has unsupported overhangs beyond its
        // first/last inter-window link; trimming disables those dead-end edges
        // before the graph is compressed into the assembly graph. Called once:
        // the trim criterion (first/last anchor carrying any inter-window edge)
        // is unaffected by disabling backbone edges, so a second standalone call
        // would report the same count -- it only advances when interleaved with
        // other filters that change inter-window connectivity.
        {
            const uint64_t trimmed =
                assembler.shasta2AnchorGraph->trimBackbones(anchorWindows, *shasta2Journeys);
            cout << timestamp << "trimBackbones: disabled backbone edges for "
                 << trimmed << " overhang anchors." << endl;
        }

        // Tip cleanup before export. Two complementary passes:
        //   (B) removeHetArmTips fixes the tip SOURCE -- het/hom allele arms left
        //       one-sided because addHetEdge dropped one of the arm's two flank
        //       edges (empty read intersection / non-forward offset). It disables
        //       the arm's surviving edge and cascades to any stranded hom.
        //   (A) removeAnchorGraphTips is a general safety net that removes any
        //       remaining interior one-sided anchor (e.g. from inter-window
        //       coverage gating or journey-tie edge drops), while preserving real
        //       contig/telomere ends (anchors with an inter-window edge) and
        //       backbone window boundaries.
        // Iterate the pair until neither removes anything, since each can expose
        // a new tip for the other.
        for(uint64_t tipPass = 0; ; ++tipPass) {
            const uint64_t hetTips =
                assembler.shasta2AnchorGraph->removeHetArmTips(*shasta2Anchors);
            const uint64_t genTips =
                assembler.shasta2AnchorGraph->removeAnchorGraphTips(
                    *shasta2Anchors, anchorWindows, *shasta2Journeys);
            if(hetTips == 0 && genTips == 0) break;
        }

        assembler.shasta2AnchorGraph->writeGfa("Shasta2AnchorGraph.gfa", &anchorWindows);
        assembler.shasta2AnchorGraph->writeCsv("Shasta2AnchorGraph.csv");
        cout << timestamp << "Wrote Shasta2AnchorGraph.gfa / .csv" << endl;

        // Persist the anchor graph as binary, not just the GFA/CSV views.
        //
        // saveForShasta2 writes the shasta2-compatible MemoryMapped format
        // (a boost archive of shasta2::AnchorGraph). This is the ONLY file
        // shasta2 can load via --external-anchor-graph-name; dinara's own
        // save() format is a different, incompatible archive and will segfault
        // shasta2 if passed there. Write it into the run's Data/ directory under
        // the same "Data/Shasta2-Shasta2AnchorGraph" name the binary-data path
        // uses, so the standard --external-anchor-graph-name works unchanged.
        string externalAnchorGraphName =
            assembler.shasta2MappedMemoryOwner().largeDataName("Shasta2AnchorGraph");
        if(externalAnchorGraphName.empty()) {
            // No binary-data directory (e.g. memory-mode anonymous): fall back to
            // a plain file in the working directory.
            externalAnchorGraphName = "Shasta2ExternalAnchorGraph";
        }
        externalAnchorGraphName =
            std::filesystem::absolute(externalAnchorGraphName).string();
        assembler.shasta2AnchorGraph->saveForShasta2(
            externalAnchorGraphName, *shasta2Anchors, &journeyTieDropMap);
        cout << timestamp << "Wrote shasta2 anchor graph. Use "
             << "--external-anchor-graph-name " << externalAnchorGraphName << endl;
    }

    // Build the assembly graph from the anchor graph. This collapses each
    // maximal non-branching anchor chain into a single segment (compress), so
    // the graph is far smaller than the anchor graph while preserving its
    // topology. Iterative short-tip removal + compress then cleans dangling
    // ends, mirroring the cleanup used in the later (disabled) pipeline.
    {
        cout << timestamp << "Creating Shasta2AssemblyGraph (with het bubbles) from "
             << "the anchor graph..." << endl;
        Shasta2AssemblyGraphOptions shasta2AssemblyGraphOptions;
        assembler.shasta2AssemblyGraph = make_shared<Shasta2AssemblyGraph>(
            *shasta2Anchors,
            *shasta2Journeys,
            *assembler.shasta2AnchorGraph,
            anchorWindows,
            shasta2AssemblyGraphOptions);
        auto& shasta2AssemblyGraph = assembler.shasta2AssemblyGraph;
        shasta2AssemblyGraph->compress();
        shasta2AssemblyGraph->writeGfa("Shasta2AssemblyGraph.gfa");

        // Iterative tip removal + compress. Shorter tips are processed first so
        // their removal can expose longer ones; loop until nothing changes.
        {
            const uint32_t maxTipWindows = 3;
            const uint64_t maxTipLength = (maxTipWindows - 1) * averageReadLength;
            for(uint64_t cleanRound = 0; ; cleanRound++) {
                uint64_t changeCount = 0;
                changeCount += shasta2AssemblyGraph->removeShortTips(maxTipWindows, maxTipLength);
                shasta2AssemblyGraph->compress();
                if(changeCount == 0) break;
            }
        }
        shasta2AssemblyGraph->writeGfa("Shasta2AssemblyGraph-cleaned.gfa");
        cout << timestamp << "Wrote Shasta2AssemblyGraph.gfa / -cleaned.gfa" << endl;
    }

    // Store elapsed time for assembly.
    const auto steadyClock1 = std::chrono::steady_clock::now();
    const auto userClock1 = boost::chrono::process_user_cpu_clock::now();
    const auto systemClock1 = boost::chrono::process_system_cpu_clock::now();
    const double elapsedTime = 1.e-9 * double((
        std::chrono::duration_cast<std::chrono::nanoseconds>(steadyClock1 - steadyClock0)).count());
    const double userTime = 1.e-9 * double((
        boost::chrono::duration_cast<boost::chrono::nanoseconds>(userClock1 - userClock0)).count());
    const double systemTime = 1.e-9 * double((
        boost::chrono::duration_cast<boost::chrono::nanoseconds>(systemClock1 - systemClock0)).count());
    const double averageCpuUtilization =
        (userTime + systemTime) / (double(std::thread::hardware_concurrency()) * elapsedTime);
    assembler.storeAssemblyTime(elapsedTime, averageCpuUtilization);

    // Store peak memory usage.
    const uint64_t peakMemoryUsage = getPeakMemoryUsage();
    assembler.storePeakMemoryUsage(peakMemoryUsage);

    // Store other performance information.
    assembler.assemblerInfo->threadCount = threadCount;
    assembler.assemblerInfo->virtualCpuCount = std::thread::hardware_concurrency();
    assembler.assemblerInfo->totalAvailableMemory = getTotalPhysicalMemory();

    // Write a summary of read information.
    assembler.writeReadsSummary();

    // Write the assembly summary.
    ofstream html("AssemblySummary.html");
    assembler.writeAssemblySummary(html);
    ofstream json("AssemblySummary.json");
    assembler.writeAssemblySummaryJson(json);
    ofstream htmlIndex("index.html");
    assembler.writeAssemblyIndex(htmlIndex);

    if(not assembler.saveBinaryDataDirectory.empty()) {
        assembler.waitForSaveBinaryDataThreads();
    }

    performanceLog << timestamp << endl;
    performanceLog << "Assembly time statistics:\n"
        "    Elapsed seconds: " << elapsedTime << "\n"
        "    Elapsed minutes: " << elapsedTime/60. << "\n"
        "    Elapsed hours:   " << elapsedTime/3600. << "\n";
    performanceLog << "Average CPU utilization: " << averageCpuUtilization << endl;
    performanceLog << "Peak Memory usage: " << peakMemoryUsage << " bytes = " <<
        int(std::round(double(peakMemoryUsage) / (1024. * 1024. * 1024.)) ) << " GiB" << endl;

    return;
}



// This function sets nr_overcommit_hugepages for 2MB pages
// to a little below total memory.
// If the setting needs to be modified, it acquires
// root privilege via sudo. This may result in the
// user having to enter a password.
void dinara::main::setupHugePages()
{

    // Get the total memory size.
    const uint64_t totalMemoryBytes = sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);

    // Figure out how much memory we want to allow for 2MB pages.
    const uint64_t MB = 1024 * 1024;
    const uint64_t GB = MB * 1024;
    const uint64_t maximumHugePageMemoryBytes = totalMemoryBytes - 8 * GB;
    const uint64_t maximumHugePageMemoryHugePages = maximumHugePageMemoryBytes / (2 * MB);

    // Check what we have it set to.
    const string fileName = "/sys/kernel/mm/hugepages/hugepages-2048kB/nr_overcommit_hugepages";
    ifstream file(fileName);
    if(!file) {
        throw runtime_error("Error opening " + fileName + " for read.");
    }
    uint64_t currentValue = 0;
    file >> currentValue;
    file.close();

    // If it's set to at least what we want, don't do anything.
    // When this happens, root access is not required.
    if(currentValue >= maximumHugePageMemoryHugePages) {
        return;
    }

    // Use sudo to set.
    const string command =
        "sudo sh -c \"echo " +
        to_string(maximumHugePageMemoryHugePages) +
        " > " + fileName + "\"";
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " running command: " + command);
    }

}



// Implementation of --command saveBinaryData.
// This copies Data to DataOnDisk.
void dinara::main::saveBinaryData(
    const AssemblerOptions& assemblerOptions)
{
    DINARA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "saveBinaryData");

    // Locate the Data directory.
    const string dataDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/Data";
    if(!std::filesystem::exists(dataDirectory)) {
        throw runtime_error(dataDirectory + " does not exist, nothing done.");
    }

    // Check that the DataOnDisk directory does not exist.
    const string dataOnDiskDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/DataOnDisk";
    if(std::filesystem::exists(dataOnDiskDirectory)) {
        throw runtime_error(dataOnDiskDirectory + " already exists, nothing done.");
    }

    // Copy Data to DataOnDisk.
    const string command = "cp -rp " + dataDirectory + " " + dataOnDiskDirectory;
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " running command:\n" + command);
    }
    cout << "Binary data successfully saved." << endl;
}



// Implementation of --command cleanupBinaryData.
void dinara::main::cleanupBinaryData(
    const AssemblerOptions& assemblerOptions)
{
    DINARA_ASSERT(assemblerOptions.commandLineOnlyOptions.command == "cleanupBinaryData");

    // Locate the Data directory.
    const string dataDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/Data";
    if(!std::filesystem::exists(dataDirectory)) {
        cout << dataDirectory << " does not exist, nothing done." << endl;
        return;
    }

    // Unmount it and remove it.
    ::system(("sudo umount " + dataDirectory).c_str());
    const int errorCode = ::system(string("rm -rf " + dataDirectory).c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " removing " + dataDirectory);
    }
    cout << "Cleanup of " << dataDirectory << " successful." << endl;

    // If the DataOnDisk directory exists, create a symbolic link
    // Data->DataOnDisk.
    const string dataOnDiskDirectory =
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory + "/DataOnDisk";
    if(std::filesystem::exists(dataOnDiskDirectory)) {
        std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);
        const string command = "ln -s DataOnDisk Data";
        ::system(command.c_str());
    }

}

// Implementation of --command explore.
void dinara::main::explore(
    const AssemblerOptions& assemblerOptions)
{
    // If a paf file was specified, find its absolute path
    // before we switch to the assembly directory.
    string alignmentsPafFileAbsolutePath;
    if(not assemblerOptions.commandLineOnlyOptions.alignmentsPafFile.empty()) {
        if(!std::filesystem::exists(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile)) {
            throw runtime_error(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile + " not found.");
        }
        if(!std::filesystem::is_regular_file(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile)) {
            throw runtime_error(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile + " is not a regular file.");
        }
        alignmentsPafFileAbsolutePath = filesystem::getAbsolutePath(assemblerOptions.commandLineOnlyOptions.alignmentsPafFile);
    }

    // Go to the assembly directory.
    std::filesystem::current_path(assemblerOptions.commandLineOnlyOptions.assemblyDirectory);

    // Check that we have the binary data.
    if(!std::filesystem::exists("Data")) {
        throw runtime_error("Binary directory \"Data\" not available "
        " in assembly directory " +
        assemblerOptions.commandLineOnlyOptions.assemblyDirectory +
        ". Use \"--memoryMode filesystem\", possibly followed by "
        "\"--command saveBinaryData\" and \"--command cleanupBinaryData\" "
        "if you want to make sure the binary data are persistently available on disk. "
        "See the documentations are some of these options require root access."
        );
        return;
    }

    // Create the Assembler.
    Assembler assembler("Data/", false, 1, 0);

    // Set up the consensus caller.
    if(assembler.getReads().representation == 1) {
        cout << "Setting up consensus caller " <<
            assemblerOptions.assemblyOptions.consensusCaller << endl;
    }
    assembler.setupConsensusCaller(assemblerOptions.assemblyOptions.consensusCaller);

    // Access all available binary data.
    assembler.httpServerData.assemblerOptions = &assemblerOptions;
    assembler.accessAllSoft();

    string executablePath = filesystem::executablePath();
    // On Linux it will be something like - `/path/to/install_root/bin/dinara`

    string executableBinPath = executablePath.substr(0, executablePath.find_last_of('/'));
    string installRootPath = executableBinPath.substr(0, executableBinPath.find_last_of('/'));
    string docsPath = installRootPath + "/docs";

    if (std::filesystem::is_directory(docsPath)) {
        assembler.httpServerData.docsDirectory = docsPath;
    } else {
        cout << "Documentation is not available." << endl;
        assembler.httpServerData.docsDirectory = "";
    }

    // Load the paf file, if one was specified.
    if(not alignmentsPafFileAbsolutePath.empty()) {
        assembler.loadAlignmentsPafFile(alignmentsPafFileAbsolutePath);
    }

    // Start the http server.
    bool localOnly;
    bool sameUserOnly;
    if(assemblerOptions.commandLineOnlyOptions.exploreAccess == "user") {
        localOnly = true;
        sameUserOnly = true;
    } else if(assemblerOptions.commandLineOnlyOptions.exploreAccess == "local") {
        localOnly = true;
        sameUserOnly = false;
    } else if (assemblerOptions.commandLineOnlyOptions.exploreAccess == "unrestricted"){
        localOnly = false;
        sameUserOnly = false;
    } else {
        throw runtime_error("Invalid value specified for --exploreAccess. "
            "Only use this option if you understand its security implications."
        );
    }
    assembler.explore(
        assemblerOptions.commandLineOnlyOptions.port,
        localOnly,
        sameUserOnly);
}


void dinara::main::listCommands()
{
    cout << "Valid commands are:" << endl;
    for(const string& command: commands) {
        cout << command << endl;
    }
}
