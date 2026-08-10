// AssemblerWindowAbpoaGraph.cpp
//
// Per-window abPOA MSA construction: one graph per window, each member read
// aligned in a small number of capped-span batched align+add calls (grouping
// several consecutive inter-anchor segments per call) rather than one call
// per tiny inter-anchor segment.
//
// DESIGN NOTE — abPOA stays base-level; anchors are pins, not nodes.
//   abPOA is a Partial Order Alignment library: every graph node is exactly one
//   residue (one base). There is no mode in which a node holds a k=50 anchor
//   k-mer; the DP, the row-column msa_base matrix, and the consensus all assume
//   residue-level nodes. So the backbone is seeded as its whole BASE sequence
//   (one node per base, node id = bbPos + 2), NOT as anchor nodes.
//
//   The k=50 anchors are used only as PINS: AnchorPin{backbonePos, readOrdinal}
//   defines where each member's shared-anchor span is DP-aligned against the
//   base-level backbone subgraph (abpoa_align_sequence_to_subgraph with a node
//   range). They constrain placement; they never become graph content.
//
//   This base resolution is a requirement, not a limitation: het sites and the
//   2-base allele anchors minted from this MSA are SUB-anchor features (a single
//   SNP lives inside one 50-base anchor's span). An anchor-node graph — which is
//   what Shasta2AnchorGraph already provides for assembly topology — is blind to
//   that signal by construction. The correct hybrid is:
//       anchor graph (k=50, topology)  --pins-->  per-window abPOA (base-level)
//   resolving fine structure between anchors, from which finer 2-base allele
//   anchors are derived and handed back to the anchor-graph world.
//
//   Performance: per-base backbone is cheap here. Seeding an empty graph adds
//   nodes with no DP.

//
// For each anchor window we build ONE abPOA graph:
//   1. Seed the graph with the backbone read's base sequence spanning the
//      window's backbone anchors. The backbone is sequence 0 (the spine);
//      backbone base position p maps to graph node id p + 2. Seeding an empty
//      graph does no DP, so the backbone may be arbitrarily long.
//   2. For every other member read, find the backbone anchors it shares with
//      the backbone (its "pins"). Walk the pins left to right, greedily
//      grouping consecutive inter-anchor segments into one align+add call as
//      long as the group's combined backbone span stays under
//      maxBatchSpanBases -- not one call per inter-anchor segment (see the
//      "Why capped-span batching" note below for why this is both correct
//      and much faster). The member occupies a single MSA row.
//   3. After all members are added, extract the row-column MSA matrix
//      (abc->msa_base) and write it as a CSV per window.
//
// Why capped-span batching, not one call per inter-anchor segment: earlier
// code aligned each member's inter-anchor gaps as separate
// abpoa_align_sequence_to_subgraph/abpoa_add_subgraph_alignment calls,
// reasoning that abPOA's DP is O(graph_span x query_len) and a whole ~150kb
// window in one shot would make an unbanded span x span matrix (~128 GiB) and
// abort. Splitting at every inter-anchor gap avoided that, but bought nothing
// for DP cost while paying abpoa_add_subgraph_alignment's internal full-graph
// topological re-sort (unconditional on every call) once per gap -- measured
// at 497 members over 110563 segments in one real run, i.e. a call ~222x more
// often than needed, 470s of a ~475s total vs. 2.6s of actual DP.
//
// The natural fix -- one call per member, spanning its first pin to its last
// -- is NOT safe in general: abPOA's DP buffer is allocated as
// matrix_row_n x qlen regardless of the wb/wf band parameters (confirmed by
// reading abpoa_align_simd.c's simd_abpoa_var macro: dp_sn scales with the
// full qlen, matrix_row_n with the full node range; wb/wf only bound which
// cells get COMPUTED, not the buffer size). A member's own span is not
// reliably bounded by its own read length either -- a spurious or
// repeat-driven shared anchor can pair positions far apart -- so an
// unconditional one-call-per-member batch hit exactly the "unbanded
// span x span matrix" blowup the per-segment split was meant to avoid, just
// at member-span instead of window-span scale (this crashed on a larger,
// denser test region). Capping each batch's span at a small fixed constant
// (maxBatchSpanBases) bounds worst-case memory per call unconditionally,
// independent of read length, window size, or anchor placement, while still
// cutting the resort count by roughly maxBatchSpanBases / (anchor spacing).
//
// Intervals are computed per OrientedReadId (read + strand), never per ReadId;
// a read's two strands are distinct OrientedReadIds with distinct pins.
//
// Anchor/sequence handling mirrors AssemblerAnchorWindowsClean.cpp::createWindow
// and PhasingKmeansAlign.cpp::abpoaMsaRun.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "Base.hpp"
#include "Marker.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "DINARA_ASSERT.hpp"
#include "timestamp.hpp"
#include "WindowIntervalPoa.hpp"

#include "abpoa/abpoa.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// A shared backbone anchor for one oriented member read.
struct AnchorPin {
    int backbonePos;       // 0-based base position into the backbone sequence.
    uint32_t readOrdinal;  // Marker ordinal of this anchor on the member read.
};

// A member oriented read and its shared-anchor pins (in backbone order).
struct OrientedReadMember {
    OrientedReadId orientedReadId;
    vector<AnchorPin> pins;   // Sorted ascending by backbonePos.
    int backboneSpan = 0;     // pins.back - pins.front; used for add order.
};

inline uint8_t baseToInt(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4; // N
    }
}

// abPOA encodes bases 0..3 = ACGT, 4 = N, and the gap value equals abpt->m.
inline char intToBaseOrGap(uint8_t v, int gapValue) {
    if(int(v) == gapValue) return '-';
    switch(v) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

} // anonymous namespace

void Assembler::computeWindowAbpoaGraphs(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const string& outputPrefix,
    uint64_t threadCount) const
{
    if(!markers) {
        cout << timestamp << "computeWindowAbpoaGraphs: markers not available. Skipping." << endl;
        return;
    }
    const Reads& reads = getReads();
    const uint32_t kHalf = assemblerInfo.isOpen ? uint32_t(assemblerInfo->k / 2) : 0U;

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
        if(threadCount == 0) threadCount = 1;
    }

    cout << timestamp << "computeWindowAbpoaGraphs: processing "
         << anchorWindows.size() << " windows with " << threadCount
         << " threads (prefix=\"" << outputPrefix << "\")." << endl;

    std::atomic<uint64_t> nextWindow{0};
    std::atomic<uint64_t> windowsWritten{0};
    std::atomic<uint64_t> windowsSkipped{0};
    std::atomic<uint64_t> membersAligned{0};
    std::atomic<uint64_t> segmentsAligned{0};  // number of batched (capped-span) align+add calls, across all members
    std::atomic<uint64_t> snarlDetectNs{0};  // time in IpoaNodeVisitors::build + findLeafSnarlsFromGraph only
    std::atomic<uint64_t> generateRcMsaNs{0};  // time in abpoa_generate_rc_msa only
    std::atomic<uint64_t> csvWriteNs{0};       // time in the fprintf-per-cell CSV loop only
    std::atomic<uint64_t> alignToSubgraphNs{0};  // time in abpoa_align_sequence_to_subgraph only
    std::atomic<uint64_t> addSubgraphAlnNs{0};   // time in abpoa_add_subgraph_alignment only

    // Extract the base sequence of an oriented read between two marker ordinals
    // (midpoint to midpoint), encoded as 0123/4. Mirrors the extraction used in
    // PhasingKmeansAlign.cpp / createWindow.
    const auto extractSeq0123 = [&](OrientedReadId oid, uint32_t ordA, uint32_t ordB) -> vector<uint8_t> {
        vector<uint8_t> seq;
        if(ordA >= ordB) return seq;
        const auto readMarkers = (*markers)[oid.getValue()];
        if(ordB >= readMarkers.size()) return seq;
        const uint32_t beginPos = readMarkers[ordA].position + kHalf;
        const uint32_t endPos   = readMarkers[ordB].position + kHalf;
        if(endPos <= beginPos) return seq;
        seq.reserve(endPos - beginPos);
        for(uint32_t pos = beginPos; pos < endPos; pos++) {
            seq.push_back(baseToInt(reads.getOrientedReadBase(oid, pos).character()));
        }
        return seq;
    };

    auto worker = [&]() {
        // Each thread owns its own abPOA instances — no shared mutable state.
        for(;;) {
            const uint64_t wi = nextWindow.fetch_add(1);
            if(wi >= anchorWindows.size()) break;
            const AnchorWindow& window = anchorWindows[wi];

            const OrientedReadId backboneOid = window.backboneOrientedReadId;

            // Ordered backbone journey positions forming the spine.
            const auto backboneJourney = journeys[backboneOid];
            vector<uint32_t> bbPositions = window.filteredBackbonePositions;
            if(bbPositions.empty()) {
                for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                    bbPositions.push_back(pos);
                }
            }
            if(bbPositions.size() < 2) {
                windowsSkipped.fetch_add(1);
                continue;
            }

            // Backbone base coordinates: midpoint of first and last backbone
            // anchor. The backbone sequence spans [bbStartBase, bbEndBase).
            const uint32_t firstOrdinal = anchors.getOrdinal(
                backboneJourney[bbPositions.front()], backboneOid);
            const uint32_t lastOrdinal = anchors.getOrdinal(
                backboneJourney[bbPositions.back()], backboneOid);
            const auto bbMarkers = (*markers)[backboneOid.getValue()];
            if(lastOrdinal >= bbMarkers.size() || firstOrdinal >= lastOrdinal) {
                windowsSkipped.fetch_add(1);
                continue;
            }
            const uint32_t bbStartBase = bbMarkers[firstOrdinal].position + kHalf;
            const uint32_t bbEndBase   = bbMarkers[lastOrdinal].position + kHalf;
            if(bbEndBase <= bbStartBase) {
                windowsSkipped.fetch_add(1);
                continue;
            }

            // Backbone sequence (0123/4).
            vector<uint8_t> backboneSeq;
            backboneSeq.reserve(bbEndBase - bbStartBase);
            for(uint32_t pos = bbStartBase; pos < bbEndBase; pos++) {
                backboneSeq.push_back(baseToInt(reads.getOrientedReadBase(backboneOid, pos).character()));
            }
            const int backboneLen = int(backboneSeq.size());
            if(backboneLen < 2) {
                windowsSkipped.fetch_add(1);
                continue;
            }

            // Consume the shared-anchor pins persisted on each member during
            // window construction (readJourneyPos, backboneJourneyPos), and
            // convert them to (backbone base position, read marker ordinal).
            // No journey re-intersection here — that work was done once in
            // createWindow.
            vector<OrientedReadMember> members;
            members.reserve(window.readIntervals.size());

            for(const AnchorWindowReadInterval& ri : window.readIntervals) {
                const OrientedReadId oid = ri.orientedReadId;
                if(oid == backboneOid) continue;
                if(oid.getValue() >= journeys.size()) continue;
                if(ri.sharedPins.size() < 2) continue;
                const auto readMarkers = (*markers)[oid.getValue()];

                vector<AnchorPin> pins;
                pins.reserve(ri.sharedPins.size());
                for(const AnchorWindowSharedPin& sp : ri.sharedPins) {
                    // Backbone base position: midpoint of the anchor on the
                    // backbone, relative to backbone sequence start. Must lie
                    // within the seeded backbone span.
                    if(sp.backboneJourneyPos >= uint32_t(backboneJourney.size())) continue;
                    const Shasta2AnchorId aid = backboneJourney[sp.backboneJourneyPos];
                    const uint32_t bbOrd = anchors.getOrdinal(aid, backboneOid);
                    if(bbOrd >= bbMarkers.size()) continue;
                    const uint32_t bbMid = bbMarkers[bbOrd].position + kHalf;
                    if(bbMid < bbStartBase || bbMid >= bbEndBase) continue;
                    const int backbonePos = int(bbMid - bbStartBase);

                    // Read marker ordinal of the same anchor on this oriented read.
                    const uint32_t readOrd = anchors.getOrdinal(aid, oid);
                    if(readOrd >= readMarkers.size()) continue;

                    pins.push_back(AnchorPin{backbonePos, readOrd});
                }
                if(pins.size() < 2) continue;

                // Pins arrive in backbone order from construction; re-sort
                // defensively in case any were dropped above.
                std::sort(pins.begin(), pins.end(),
                    [](const AnchorPin& a, const AnchorPin& b) {
                        return a.backbonePos < b.backbonePos;
                    });

                OrientedReadMember m;
                m.orientedReadId = oid;
                m.pins = std::move(pins);
                m.backboneSpan = m.pins.back().backbonePos - m.pins.front().backbonePos;
                if(m.backboneSpan <= 0) continue;
                members.push_back(std::move(m));
            }

            // Add members from the biggest backbone span to the smallest, so the
            // longest-overlapping reads shape the graph spine before shorter
            // ones. Ties broken by oriented read id for determinism.
            std::sort(members.begin(), members.end(),
                [](const OrientedReadMember& a, const OrientedReadMember& b) {
                    if(a.backboneSpan != b.backboneSpan) return a.backboneSpan > b.backboneSpan;
                    return a.orientedReadId < b.orientedReadId;
                });

            const int nMembers = int(members.size());
            const int totalSeqs = 1 + nMembers; // backbone + members.

            // --- Build the per-window abPOA graph ---
            abpoa_t* ab = abpoa_init();
            abpoa_para_t* abpt = abpoa_init_para();
            abpt->align_mode = ABPOA_GLOBAL_MODE;
            // sub_aln: align reads to a subgraph (node-range restricted). Without
            // it abpoa_align_sequence_to_subgraph ignores the node bounds.
            abpt->sub_aln = 1;
            abpt->inc_path_score = 1;
            abpt->out_msa  = 1;   // populate abc->msa_base (row-column MSA).
            // out_cons MUST be 0 here. We add each member as multiple segments
            // sharing one read_id, so a node's coverage (summed read_id bits)
            // can exceed the per-cluster sequence count. abPOA's consensus path
            // (abpoa_cons_phred_score) then asserts n_cov > n_seq and aborts
            // ("unexpected n_cov/n_seq"). We only need the MSA matrix, so skip
            // consensus entirely.
            abpt->out_cons = 0;
            abpt->out_gfa  = 0;
            abpt->sort_input_seq = 0;   // keep backbone as seq 0 (the spine).
            abpt->progressive_poa = 0;
            abpt->max_n_cons = 1;
            abpt->use_qv = 0;
            abpoa_post_set_para(abpt);

            ab->abs->n_seq = totalSeqs;

            // Seed backbone as read_id 0. On an empty graph,
            // abpoa_add_subgraph_alignment creates nodes 2..L+1 directly (no DP),
            // so backbone base position p maps to node id p + 2. Also capture
            // that mapping directly via qpos_to_node_id (backbonePath[p] = node
            // id of backbone base p) for findLeafSnarlsFromGraph below, rather
            // than assuming the p+2 formula or re-deriving it from graph edges
            // afterward -- see the design note on findLeafSnarlsFromGraph
            // (WindowIntervalPoa.hpp) for why edge-based reconstruction is
            // unreliable for multi-segment members.
            vector<int> backbonePath(size_t(backboneLen), -1);
            {
                abpoa_res_t res{};
                res.graph_cigar = nullptr;
                res.n_cigar = 0;
                abpoa_align_sequence_to_subgraph(
                    ab, abpt, 0, 1, backboneSeq.data(), backboneLen, &res);
                abpoa_add_subgraph_alignment(
                    ab, abpt, 0, 1, backboneSeq.data(), nullptr,
                    backboneLen, backbonePath.data(), res, 0, totalSeqs, 0);
                if(res.n_cigar) free(res.graph_cigar);
            }

            // Add each member segment-by-segment between consecutive shared
            // anchors. All segments of one member use the same read_id, so the
            // member fills a single MSA row. Also accumulate the flat list of
            // (nodeId, memberIndex) pairs actually visited (captured directly
            // from each segment's own qpos_to_node_id output) for
            // findLeafSnarlsFromGraph below -- one growing flat array here,
            // turned into a proper CSR (IpoaNodeVisitors) once the final node
            // count is known, rather than one vector<bool> per member sized
            // to the whole graph (see the design note on IpoaNodeVisitors).
            const size_t nMembersSz = size_t(nMembers);
            vector<IpoaMemberInfo> memberInfos(nMembersSz);
            vector<std::pair<int, int>> nodeVisits;
            nodeVisits.reserve(size_t(backboneLen) * 2);  // rough amortization hint
            for(int i = 0; i < nMembers; i++) {
                const OrientedReadMember& m = members[size_t(i)];
                const int readId = i + 1;
                IpoaMemberInfo& info = memberInfos[size_t(i)];
                info.oid = m.orientedReadId;
                // Coverage bound: this member has no information outside its
                // own first/last shared anchor with the backbone (see the
                // design note on IpoaMemberInfo -- without this, every
                // position a member doesn't reach looks like a deletion).
                info.bbCovBegin = std::uint32_t(m.pins.front().backbonePos);
                info.bbCovEnd = std::uint32_t(m.pins.back().backbonePos) + 1;

                // Batch consecutive inter-anchor segments into chunks capped
                // at maxBatchSpanBases backbone bases each, instead of one
                // align+add call per tiny inter-anchor segment. abPOA's
                // subgraph-add unconditionally re-sorts the whole graph
                // topologically on every call (abpoa_graph.c:
                // abpoa_add_subgraph_alignment), so per-segment calls paid
                // that O(graph size) resort ~222 times per member on average
                // (110563 segments / 497 members measured); chunking cuts
                // that by roughly maxBatchSpanBases / (anchor spacing).
                //
                // The cap is essential, NOT an optional tuning knob: abPOA's
                // DP buffer is allocated as matrix_row_n x qlen regardless of
                // the wb/wf band parameters (confirmed by reading
                // abpoa_align_simd.c's simd_abpoa_var macro -- dp_sn scales
                // with the FULL qlen, matrix_row_n with the full node range;
                // wb/wf only bound which cells get COMPUTED, not the buffer
                // size). Batching one member's ENTIRE walk into a single call
                // (an earlier version of this code) hits exactly the
                // "unbanded span x span matrix" blowup the original
                // per-segment split existed to avoid -- just at
                // member-span scale instead of window scale -- and a member's
                // span is not reliably bounded by its own read length (a
                // spurious/repeat-driven shared anchor can pair positions far
                // apart). A fixed cap bounds worst-case memory per call
                // unconditionally, independent of read length, window size,
                // or anchor placement.
                //
                // The pins skipped over within one chunk are no longer forced
                // as exact hard pins (the DP is free to choose its own path
                // across the chunk), but findLeafSnarlsFromGraph only depends
                // on the OUTER coverage bounds (bbCovBegin/bbCovEnd, set above
                // from front()/back()), so this doesn't affect correctness
                // there. HiFi reads are accurate enough that the free
                // alignment should still pass through the same true-match
                // positions the anchors mark, except occasionally within
                // homopolymer/STR runs where the placement of an indel is
                // genuinely ambiguous -- exactly what the leaf-snarl repeat
                // gate already filters out.
                // Cap BOTH the backbone-side span and the read-side span of a
                // batch. Capping only the backbone side is not enough: a read
                // can carry a large insertion relative to the backbone (a
                // duplication/repeat expansion) inside a short backbone
                // range, which grows the query length independent of the
                // backbone span and can still blow up the row_n x qlen DP
                // buffer even with a tight backbone cap.
                //
                // Neither cap alone is actually sufficient, though: abPOA's
                // subgraph node range is NOT simply the nodes whose backbone
                // position falls in [incBeg,incEnd] -- abpoa_subgraph_nodes
                // calls abpoa_upstream_index/abpoa_downstream_index, which
                // walk in/out edges and EXPAND the range until it is
                // topologically closed. In a window with many members already
                // added, that expansion can pull in a much wider graph-index
                // range than the backbone-position distance would suggest
                // (every prior member's substitution/insertion nodes in that
                // area interleave into the index space). This was observed
                // directly: a 16 GiB single-allocation failure on a larger,
                // denser test region even with both position-based caps in
                // place. So after computing the actual subgraph range, check
                // its real cell count (matrix_row_n x qlen) against a hard
                // budget and shrink the chunk (fewer segments) until it fits
                // -- down to a single segment as the unconditional floor,
                // which matches the pre-batching per-segment behavior that
                // never had this problem (a single inter-anchor gap's own
                // qlen is small enough that even a very expanded row_n stays
                // within budget).
                constexpr int maxBatchSpanBases = 2000;
                constexpr int64_t maxBatchCells = 20'000'000; // ~200MB at ~10 bytes/cell
                for(size_t s = 0; s + 1 < m.pins.size(); ) {
                    size_t e = s;
                    while(e + 1 < m.pins.size() &&
                          m.pins[e + 1].backbonePos - m.pins[s].backbonePos <= maxBatchSpanBases &&
                          int(m.pins[e + 1].readOrdinal - m.pins[s].readOrdinal) <= maxBatchSpanBases) {
                        // The shared anchors between two oriented reads are
                        // guaranteed monotone by the underlying alignment:
                        // ordering pins by backbone position must also order
                        // them strictly ascending in read ordinal. A
                        // violation means a broken invariant upstream, not
                        // normal data, so assert.
                        DINARA_ASSERT(m.pins[e + 1].backbonePos > m.pins[e].backbonePos);
                        DINARA_ASSERT(m.pins[e + 1].readOrdinal > m.pins[e].readOrdinal);
                        e++;
                    }
                    if(e == s) e = s + 1; // guarantee progress even if one segment alone exceeds the cap.

                    // Shrink the chunk if the ACTUAL graph-index range and/or
                    // ACTUAL extracted read length (queried/built below) turns
                    // out larger than the position-based caps above assumed.
                    // The read length MUST be measured via extractSeq0123
                    // (real base positions, readMarkers[ord].position), not
                    // approximated from the readOrdinal difference: ordinals
                    // count MARKERS, not bases, and marker density can be
                    // sparse in low-complexity/repeat stretches, so a small
                    // ordinal gap can still span a huge base range (this was
                    // the actual cause of a 16 GiB allocation failure that
                    // survived an earlier version of this cap, which used the
                    // ordinal difference as a proxy for read length instead
                    // of the real extracted size).
                    int incBeg = 0, incEnd = 0, excBeg = 0, excEnd = 1;
                    vector<uint8_t> segSeq;
                    while(true) {
                        incBeg = m.pins[s].backbonePos + 2;
                        incEnd = m.pins[e].backbonePos + 2 - 1;
                        excBeg = 0; excEnd = 1;
                        abpoa_subgraph_nodes(ab, abpt, incBeg, incEnd, &excBeg, &excEnd);
                        segSeq = extractSeq0123(m.orientedReadId, m.pins[s].readOrdinal, m.pins[e].readOrdinal);
                        if(e == s + 1) break; // single segment: unconditional floor, proceed regardless of size.
                        const int64_t rowN = int64_t(ab->abg->node_id_to_index[excEnd])
                            - int64_t(ab->abg->node_id_to_index[excBeg]) + 1;
                        if(rowN * int64_t(segSeq.size()) <= maxBatchCells) break;
                        e--;
                    }

                    DINARA_ASSERT(incBeg <= incEnd);
                    DINARA_ASSERT(!segSeq.empty());

                    if(getenv("DINARA_ABPOA_BATCH_DEBUG") != nullptr) {
                        const int64_t rowNDbg = int64_t(ab->abg->node_id_to_index[excEnd])
                            - int64_t(ab->abg->node_id_to_index[excBeg]) + 1;
                        fprintf(stderr, "abpoaBatch window=%llu member=%d s=%zu e=%zu rowN=%lld qlen=%zu cells=%lld\n",
                            (unsigned long long)wi, i, s, e, (long long)rowNDbg, segSeq.size(),
                            (long long)rowNDbg * (long long)segSeq.size());
                        fflush(stderr);
                    }
                    abpoa_res_t res{};
                    res.graph_cigar = nullptr;
                    res.n_cigar = 0;
                    const auto tAlign0 = std::chrono::steady_clock::now();
                    abpoa_align_sequence_to_subgraph(
                        ab, abpt, excBeg, excEnd,
                        segSeq.data(), int(segSeq.size()), &res);
                    const auto tAlign1 = std::chrono::steady_clock::now();
                    vector<int> segQpos(segSeq.size(), -1);
                    abpoa_add_subgraph_alignment(
                        ab, abpt, excBeg, excEnd,
                        segSeq.data(), nullptr,
                        int(segSeq.size()), segQpos.data(), res, readId, totalSeqs, 1);
                    const auto tAlign2 = std::chrono::steady_clock::now();
                    alignToSubgraphNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        tAlign1 - tAlign0).count()), std::memory_order_relaxed);
                    addSubgraphAlnNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                        tAlign2 - tAlign1).count()), std::memory_order_relaxed);
                    if(res.n_cigar) free(res.graph_cigar);
                    segmentsAligned.fetch_add(1);

                    for(const int nodeId : segQpos) {
                        if(nodeId < 0) continue;
                        nodeVisits.push_back({nodeId, i});
                    }

                    s = e;
                }
                membersAligned.fetch_add(1);
            }

            // --- Leaf-snarl detection (verification only; not fed into
            // anything else -- window.hetBubbles/anchor construction still
            // come exclusively from the interval-POA + WindowHetProfiles
            // pipeline). Thresholds match the auto-derived defaults used
            // there (resolveHetMinSupport's floor and the standard het gates)
            // rather than pulling in coverage estimation here. ---
            {
                const auto tSnarl0 = std::chrono::steady_clock::now();

                // Build the CSR from the flat (nodeId, memberIndex) list --
                // O(V + node_n) counting sort, no comparison sort needed.
                const IpoaNodeVisitors nodeVisitors =
                    IpoaNodeVisitors::build(ab->abg->node_n, nodeVisits);

                // DINARA_SKIP_REPEAT_GATE=1 disables the homopolymer/STR gate,
                // for verification only -- to see by eye which additional
                // sites it is actually excluding.
                const bool skipRepeatGate = (getenv("DINARA_SKIP_REPEAT_GATE") != nullptr);
                const auto snarls = findLeafSnarlsFromGraph(
                    ab->abg, backboneOid, backbonePath, memberInfos, nodeVisitors,
                    backboneSeq, bbStartBase,
                    /*minSupport=*/6, /*minVaf=*/0.12,
                    /*dropHomopolymer=*/!skipRepeatGate, /*dropRepeat=*/!skipRepeatGate);

                snarlDetectNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - tSnarl0).count()), std::memory_order_relaxed);

                for(const LeafSnarl& s : snarls) {
                    // Sequence context (verification only): dump +/-10 backbone
                    // bases around the snarl so the homopolymer/STR gate's
                    // effect can be checked by eye against the actual sequence,
                    // not just trusted because dropHomopolymer/dropRepeat=true.
                    const int64_t localStart = int64_t(s.start) - int64_t(bbStartBase);
                    const int64_t ctxBegin = std::max<int64_t>(0, localStart - 10);
                    const int64_t ctxEnd = std::min<int64_t>(backboneLen, localStart + 12);
                    string ctx;
                    for(int64_t p = ctxBegin; p < ctxEnd; p++) {
                        ctx += intToBaseOrGap(backboneSeq[size_t(p)], -1);
                    }
                    cout << "    windowAbpoaLeafSnarl window=" << window.windowId
                         << " bb=" << backboneOid
                         << " [" << s.start << "," << s.end << ")"
                         << " ctx=" << ctx
                         << " alleles=" << s.alleles.size();
                    for(const LeafSnarlAllele& al : s.alleles) {
                        cout << " {n=" << al.members.size() << " bases=";
                        for(uint8_t b : al.bases) cout << "ACGTN"[b < 4 ? b : 4];
                        cout << "}";
                    }
                    cout << endl;
                }
            }

            // --- Extract the row-column MSA matrix and write CSV ---
            const auto tRcMsa0 = std::chrono::steady_clock::now();
            abpoa_generate_rc_msa(ab, abpt);
            generateRcMsaNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - tRcMsa0).count()), std::memory_order_relaxed);
            const abpoa_cons_t* abc = ab->abc;

            bool wrote = false;
            if(abc != nullptr && abc->msa_base != nullptr && abc->msa_len > 0 && abc->n_seq > 0) {
                const auto tCsv0 = std::chrono::steady_clock::now();
                const int gapValue = abpt->m; // abPOA stores gaps as value m (=5).
                const string csvPath = outputPrefix + "window" +
                    std::to_string(uint64_t(window.windowId)) + ".msa.csv";
                FILE* f = fopen(csvPath.c_str(), "w");
                if(f != nullptr) {
                    // Header: row label, then one column per MSA column.
                    fprintf(f, "row,orientedReadId,role");
                    for(int c = 0; c < abc->msa_len; c++) fprintf(f, ",c%d", c);
                    fprintf(f, "\n");

                    // Row 0 = backbone; rows 1..nMembers = members in add order.
                    for(int rowReadId = 0; rowReadId < abc->n_seq; rowReadId++) {
                        const char* role = (rowReadId == 0) ? "backbone" : "member";
                        uint64_t oidVal;
                        if(rowReadId == 0) {
                            oidVal = uint64_t(backboneOid.getValue());
                        } else {
                            oidVal = uint64_t(members[size_t(rowReadId - 1)].orientedReadId.getValue());
                        }
                        fprintf(f, "%d,%llu,%s", rowReadId, (unsigned long long)oidVal, role);
                        const uint8_t* rowBases = abc->msa_base[rowReadId];
                        for(int c = 0; c < abc->msa_len; c++) {
                            fprintf(f, ",%c", intToBaseOrGap(rowBases[c], gapValue));
                        }
                        fprintf(f, "\n");
                    }

                    // No consensus rows: out_cons is disabled (see above).
                    fclose(f);
                    wrote = true;
                }
                csvWriteNs.fetch_add(uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - tCsv0).count()), std::memory_order_relaxed);
            }
            if(wrote) windowsWritten.fetch_add(1);
            else windowsSkipped.fetch_add(1);

            abpoa_free(ab);
            abpoa_free_para(abpt);
        }
    };

    vector<std::thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t = 0; t < threadCount; t++) {
        threads.emplace_back(worker);
    }
    for(auto& th : threads) {
        th.join();
    }

    cout << timestamp << "computeWindowAbpoaGraphs: wrote " << windowsWritten.load()
         << " MSA CSVs, skipped " << windowsSkipped.load() << " windows, aligned "
         << membersAligned.load() << " members over " << segmentsAligned.load()
         << " capped-span batches." << endl;
    cout << timestamp << "  timing (summed across all windows/threads): "
         << "leaf-snarl detection=" << (double(snarlDetectNs.load()) / 1e6) << " ms, "
         << "abpoa_generate_rc_msa=" << (double(generateRcMsaNs.load()) / 1e6) << " ms, "
         << "CSV write=" << (double(csvWriteNs.load()) / 1e6) << " ms, "
         << "abpoa_align_sequence_to_subgraph=" << (double(alignToSubgraphNs.load()) / 1e6) << " ms, "
         << "abpoa_add_subgraph_alignment=" << (double(addSubgraphAlnNs.load()) / 1e6) << " ms." << endl;
}
