// Dinara.
#include "Assembler.hpp"
#include "timestamp.hpp"
#include "chrono.hpp"
#include "Reads.hpp"

// Standard library.
#include "algorithm.hpp"
#include <queue>

using namespace dinara;
using namespace std;

namespace {
    constexpr int ASG_ET_MERGEABLE = 0;
    constexpr int ASG_ET_TIP = 1;
    constexpr int ASG_ET_MULTI_OUT = 2;
    constexpr int ASG_ET_MULTI_NEI = 3;

    inline void symmetrizeArcDeletion(StringGraph& g)
    {
        DINARA_ASSERT((g.arcs.size() & 1ULL) == 0);
        for (uint64_t arcId = 0; arcId < g.arcs.size(); arcId += 2) {
            const uint8_t del = uint8_t(g.arcs[arcId].del | g.arcs[arcId ^ 1ULL].del);
            g.arcs[arcId].del = del;
            g.arcs[arcId ^ 1ULL].del = del;
        }
    }

    inline uint32_t countOutgoingNonDeleted(const StringGraph& g, uint32_t v)
    {
        uint32_t n = 0;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) ++n;
        }
        return n;
    }

    inline uint32_t firstOutgoingNonDeleted(const StringGraph& g, uint32_t v, uint32_t& arcIdOut)
    {
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) {
                arcIdOut = arcId;
                return 1;
            }
        }
        return 0;
    }

    // Hifiasm parity: implements `asg_is_utg_end` using our stored directed arcs.
    // Note: hifiasm uses `asg_arc_a(g, v^1)` here, which corresponds to outgoing arcs of `v^1`.
    inline int stringGraphIsUtgEnd(const StringGraph& g, uint32_t v, uint32_t* nextVertex)
    {
        const uint32_t v1 = v ^ 1U;
        const uint32_t nv = countOutgoingNonDeleted(g, v1);
        if (nv == 0) {
            return ASG_ET_TIP;
        }
        if (nv > 1) {
            return ASG_ET_MULTI_OUT;
        }

        uint32_t arcId = 0;
        const uint32_t found = firstOutgoingNonDeleted(g, v1, arcId);
        DINARA_ASSERT(found == 1);
        const uint32_t to = g.arcs[arcId].to;
        if (nextVertex) {
            *nextVertex = to;
        }

        const uint32_t w = to ^ 1U;
        const uint32_t nw = countOutgoingNonDeleted(g, w);
        if (nw != 1) {
            return ASG_ET_MULTI_NEI;
        }
        return ASG_ET_MERGEABLE;
    }

    // Hifiasm parity: implements `asg_extend`.
    inline int stringGraphExtend(const StringGraph& g, uint32_t v, int maxExtReads, vector<uint32_t>& path)
    {
        path.clear();
        path.push_back(v);
        int ret = ASG_ET_MERGEABLE;
        do {
            uint32_t nextVertex = 0;
            ret = stringGraphIsUtgEnd(g, v ^ 1U, &nextVertex);
            if (ret != ASG_ET_MERGEABLE) {
                break;
            }
            path.push_back(nextVertex);
            v = nextVertex;
        } while (--maxExtReads > 0);
        return ret;
    }

    // Hifiasm parity: implements `asg_check_unambi1` (unique outgoing neighbor, else invalid).
    inline uint32_t stringGraphCheckUnambiguous1(const StringGraph& g, uint32_t v)
    {
        uint32_t next = std::numeric_limits<uint32_t>::max();
        uint32_t count = 0;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (g.arcs[arcId].del) continue;
            next = g.arcs[arcId].to;
            ++count;
            if (count > 1) break;
        }
        return (count == 1) ? next : std::numeric_limits<uint32_t>::max();
    }

    // Hifiasm parity: implements `asg_topocut_aux`.
    inline int stringGraphTopocutAux(const StringGraph& g, uint32_t v, int maxExtReads)
    {
        int nExt = 1;
        for (; nExt < maxExtReads && v != std::numeric_limits<uint32_t>::max(); ++nExt) {
            if (stringGraphCheckUnambiguous1(g, v ^ 1U) == std::numeric_limits<uint32_t>::max()) {
                --nExt;
                break;
            }
            v = stringGraphCheckUnambiguous1(g, v);
        }
        return nExt;
    }

    // -----------------------------------------------------------------------
    // Port of hifiasm `asg_bub_pop1_label` (Overlaps.cpp) and
    // `bubble_identify_worker` + `asg_arc_identify_simple_bubbles_multi`
    // (gfa_ut.cpp).
    //
    // These are used to compute a `seqVis` mask that marks bubble source,
    // interior, and merge-twin vertices.  `asg_arc_cut_length` (Phase 1)
    // skips vertices with seq_vis != 0 so that arcs from bubble branch-points
    // are not collected as short-overlap deletion candidates.
    // -----------------------------------------------------------------------

    // Per-vertex BFS state (≡ binfo_s_t in hifiasm Overlaps.h).
    //   parent    ≡ binfo_s_t::p  — parent in the BFS tree
    //   dist      ≡ binfo_s_t::d  — shortest accumulated arc-len from v0
    //   remaining ≡ binfo_s_t::r  — incoming arcs of w not yet processed by BFS
    //   visited   ≡ binfo_s_t::s  — 0 = unseen, 1 = seen
    struct BubbleBfsState {
        uint32_t parent    = 0;
        uint32_t dist      = 0;
        uint32_t remaining = 0;
        uint8_t  visited   = 0;
    };

    // Port of `asg_bub_pop1_label`.
    //
    // BFS from oriented vertex v0, tracking n_pending (incoming arcs not yet
    // processed).  Succeeds when BFS converges to a single merge vertex.
    //
    // Preconditions:
    //   • g.outgoing[] contains ONLY non-deleted arcs (after rebuildStringGraphAdjacency).
    //   • bfsState[] is pre-allocated to size vertexCount and zero-initialised.
    //
    // Returns true if a proper bubble was found.
    // On success:
    //   ready[0]   is the merge/end vertex  (≡ buf_s_t::S.a[0])
    //   visited[]  are all BFS-visited vertices excluding v0, including end
    //              (≡ buf_s_t::b)
    // BFS state is ALWAYS fully reset before returning.
    bool stringGraphBubblePop1Label(
        const StringGraph&      g,
        uint32_t                v0,
        uint64_t                maxDist,
        vector<BubbleBfsState>& bfsState,
        vector<uint32_t>&       ready,    // ≡ buf_s_t::S
        vector<uint32_t>&       visited   // ≡ buf_s_t::b
    )
    {
        if (g.readDeleted.isOpen && g.readDeleted[v0 >> 1U]) return false;
        // g.outgoing[] contains only non-deleted arcs after rebuildStringGraphAdjacency.
        if (g.outgoing[v0].size() < 2) return false;

        ready.clear();
        visited.clear();
        bfsState[v0].dist = 0;
        ready.push_back(v0);

        uint32_t nPending = 0;
        uint32_t nTips    = 0;
        uint32_t tipEnd   = uint32_t(-1);
        bool     found    = false;

        do {
            const uint32_t v  = ready.back(); ready.pop_back();
            const uint32_t dv = bfsState[v].dist;

            bool tooFar = false;
            for (const uint32_t arcId : g.outgoing[v]) {
                const uint32_t w = g.arcs[arcId].to;    // ≡ av[i].v
                const uint32_t l = g.arcs[arcId].len;   // ≡ (uint32_t)av[i].ul (non-overlap extension)

                // Cycle back to the start read → abort entirely.
                if ((w >> 1U) == (v0 >> 1U)) goto bfs_done;

                // Distance limit: abort this vertex's expansion.
                if (uint64_t(dv) + uint64_t(l) > maxDist) { tooFar = true; break; }

                {
                    BubbleBfsState& tw = bfsState[w];
                    if (tw.visited == 0) {
                        // First visit: initialise per-vertex BFS state.
                        visited.push_back(w);
                        tw.parent    = v;
                        tw.visited   = 1;
                        tw.dist      = dv + l;
                        // Count incoming arcs of w = non-deleted outgoing arcs of w^1.
                        // (≡ get_real_length(g, w^1, NULL) in hifiasm)
                        tw.remaining = uint32_t(g.outgoing[w ^ 1U].size());
                        ++nPending;
                    } else {
                        // Already visited: update shortest-distance parent.
                        if (dv + l < tw.dist) { tw.dist = dv + l; tw.parent = v; }
                    }
                    // When all incoming arcs of w have been processed, w is ready.
                    if (--(tw.remaining) == 0) {
                        if (!g.outgoing[w].empty()) {
                            ready.push_back(w);
                        } else {
                            // Dead-end tip: tolerate at most one.
                            if (nTips != 0) goto bfs_done;
                            ++nTips;
                            tipEnd = w;
                        }
                        --nPending;
                    }
                }
            }

            if (tooFar) goto bfs_done;

            // One-tip case: promote tip to merge candidate if no other work remains.
            if (nTips == 1) {
                if (tipEnd != uint32_t(-1) && nPending == 0 && ready.empty()) {
                    ready.push_back(tipEnd);
                    break;
                }
                goto bfs_done;
            }

            if (ready.empty()) goto bfs_done;

        } while (ready.size() > 1 || nPending > 0);

        found = true;

    bfs_done:
        // Reset BFS state for all visited vertices (v0 is not in visited; its
        // dist is set explicitly at the top of the next call).
        for (const uint32_t u : visited) {
            bfsState[u].visited = 0;
            bfsState[u].dist    = 0;
        }
        return found;
    }

    // Port of `asg_arc_identify_simple_bubbles_multi` + `bubble_identify_worker`.
    //
    // Returns seqVis[v] for each oriented vertex v (size = vertexCount = 2*readCount):
    //   0  normal vertex
    //   1  bubble source, interior node, or merge-node reverse complement
    //
    // Note: cross-node detection (seqVis=2, via check_if_cross) is omitted.
    // It requires asg_is_single_edge + detect_bubble_end_with_bubbles which
    // are not yet ported.  The Phase 1 skip in asg_arc_cut_length uses
    // `seq_vis != 0`, so seqVis=1 alone captures the main divergence.
    vector<uint8_t> stringGraphMarkBubbleVertices(
        const StringGraph& g,
        uint32_t           vertexCount
    )
    {
        // maxDist: safety cut-off for BFS path length (sum of arc::len values).
        // Since dv and l are uint32_t, their sum fits in uint64_t without overflow.
        // Using UINT64_MAX/4 means the distance check never fires in practice
        // (real bubbles are bounded by the graph diameter); the BFS terminates by
        // convergence (n_pending=0, |ready|=1) or cycle detection instead.
        // This matches the spirit of get_s_bub_pop_max_dist_advance which computes
        // the empirical max bubble distance, but avoids porting dfs_subgraph_s_advance.
        const uint64_t maxDist = uint64_t(-1) / 4;

        vector<uint8_t>        seqVis(vertexCount, 0);
        vector<BubbleBfsState> bfsState(vertexCount);
        vector<uint32_t>       ready;
        vector<uint32_t>       visited;

        // Single-threaded port of bubble_identify_worker (hifiasm uses kt_for).
        for (uint32_t v = 0; v < vertexCount; ++v) {
            if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) continue;
            // asg_arc_n(g,v) < 2 || get_real_length(g,v) < 2
            // Both collapse to g.outgoing[v].size() < 2 after rebuildStringGraphAdjacency.
            if (g.outgoing[v].size() < 2) continue;
            if (seqVis[v] == 1) continue;  // already marked (≡ seq_vis[v] != 1 in hifiasm)

            if (stringGraphBubblePop1Label(g, v, maxDist, bfsState, ready, visited)) {
                // beg = v, end = ready[0].
                // visited[] includes end, does NOT include beg (matches hifiasm b.b).
                const uint32_t endV = ready[0];
                for (const uint32_t u : visited) {
                    if (u == endV) continue;  // end handled separately
                    seqVis[u       ] = 1;     // interior: mark both orientations
                    seqVis[u ^ 1U  ] = 1;
                }
                seqVis[v        ] = 1;   // bubble source (forward only)
                seqVis[endV ^ 1U] = 1;   // merge-end reverse complement (not endV itself)
            }
        }

        return seqVis;
    }
}



static void sortStringGraphOutgoingByLen(Assembler& assembler, uint64_t vertexCount)
{
    for (uint32_t v = 0; v < vertexCount; ++v) {
        uint32_t* b = assembler.stringGraph.outgoing.begin(v);
        uint32_t* e = assembler.stringGraph.outgoing.end(v);
        std::stable_sort(b, e, [&](uint32_t aId, uint32_t bId) {
            const auto& a = assembler.stringGraph.arcs[aId];
            const auto& bArc = assembler.stringGraph.arcs[bId];
            return a.len < bArc.len;
        });
    }
}



static void rebuildStringGraphAdjacency(Assembler& assembler, uint64_t vertexCount)
{
    // Rebuild adjacency lists to exclude deleted arcs (hifiasm's asg_cleanup + asg_symm effect).
    assembler.stringGraph.outgoing.clear();
    if (assembler.stringGraph.incoming.isOpen()) {
        assembler.stringGraph.incoming.clear();
    }

    assembler.stringGraph.outgoing.beginPass1(uint32_t(vertexCount));
    if (assembler.stringGraph.incoming.isOpen()) {
        assembler.stringGraph.incoming.beginPass1(uint32_t(vertexCount));
    }

    for (uint64_t arcId = 0; arcId < assembler.stringGraph.arcs.size(); ++arcId) {
        const auto& a = assembler.stringGraph.arcs[arcId];
        if (a.del) continue;
        assembler.stringGraph.outgoing.incrementCount(a.from);
        if (assembler.stringGraph.incoming.isOpen()) {
            assembler.stringGraph.incoming.incrementCount(a.to);
        }
    }

    assembler.stringGraph.outgoing.beginPass2();
    if (assembler.stringGraph.incoming.isOpen()) {
        assembler.stringGraph.incoming.beginPass2();
    }
    for (uint64_t arcId = 0; arcId < assembler.stringGraph.arcs.size(); ++arcId) {
        const auto& a = assembler.stringGraph.arcs[arcId];
        if (a.del) continue;
        assembler.stringGraph.outgoing.store(a.from, uint32_t(arcId));
        if (assembler.stringGraph.incoming.isOpen()) {
            assembler.stringGraph.incoming.store(a.to, uint32_t(arcId));
        }
    }
    assembler.stringGraph.outgoing.endPass2();
    if (assembler.stringGraph.incoming.isOpen()) {
        assembler.stringGraph.incoming.endPass2();
    }

    // Sort outgoing adjacency by `len` (required by transitive reduction).
    sortStringGraphOutgoingByLen(assembler, vertexCount);
}



static void cleanupStringGraphLikeHifiasm(Assembler& assembler, uint32_t vertexCount)
{
    // Hifiasm `asg_cleanup` effect:
    // - hard-remove arcs marked deleted
    // - hard-remove arcs incident to deleted reads
    // Here we keep arc storage stable and rebuild adjacency to include only surviving arcs.
    StringGraph& g = assembler.stringGraph;
    if (g.readDeleted.isOpen) {
        for (uint64_t arcId = 0; arcId < g.arcs.size(); ++arcId) {
            auto& a = g.arcs[arcId];
            if (a.del) {
                continue;
            }
            const uint32_t fromRead = a.from >> 1U;
            const uint32_t toRead = a.to >> 1U;
            if ((fromRead < g.readDeleted.size() && g.readDeleted[fromRead]) ||
                (toRead < g.readDeleted.size() && g.readDeleted[toRead])) {
                a.del = 1;
            }
        }
    }
    rebuildStringGraphAdjacency(assembler, vertexCount);
}



static uint64_t deleteStringGraphMultiArcsLikeHifiasm(Assembler& assembler, uint32_t vertexCount)
{
    // Hifiasm `asg_arc_del_multi`: for each v, keep one v->w arc and delete additional duplicates.
    StringGraph& g = assembler.stringGraph;
    vector<uint32_t> count(vertexCount, 0);
    uint64_t removed = 0;

    for (uint32_t v = 0; v < vertexCount; ++v) {
        const span<const uint32_t> out = g.outgoing[v];
        if (out.size() < 2) {
            continue;
        }

        for (uint64_t i = out.size(); i > 0; --i) {
            const uint32_t arcId = out[i - 1];
            ++count[g.arcs[arcId].to];
        }
        for (uint64_t i = out.size(); i > 0; --i) {
            const uint32_t arcId = out[i - 1];
            const uint32_t w = g.arcs[arcId].to;
            if (--count[w] != 0) {
                g.arcs[arcId].del = 1;
                ++removed;
            }
        }
    }

    if (removed) {
        cleanupStringGraphLikeHifiasm(assembler, vertexCount);
    }
    return removed;
}



static uint64_t deleteStringGraphAsymmetricArcsLikeHifiasm(Assembler& assembler, uint32_t vertexCount)
{
    // Hifiasm `asg_arc_del_asymm`: delete u->v if (v^1)->(u^1) is absent.
    StringGraph& g = assembler.stringGraph;
    uint64_t removed = 0;

    for (uint32_t u = 0; u < vertexCount; ++u) {
        const span<const uint32_t> out = g.outgoing[u];
        for (const uint32_t arcId : out) {
            const uint32_t v = g.arcs[arcId].to;
            const uint32_t rcFrom = v ^ 1U;
            const uint32_t rcTo = u ^ 1U;

            bool found = false;
            for (const uint32_t arcRcId : g.outgoing[rcFrom]) {
                if (g.arcs[arcRcId].to == rcTo) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                g.arcs[arcId].del = 1;
                ++removed;
            }
        }
    }

    if (removed) {
        cleanupStringGraphLikeHifiasm(assembler, vertexCount);
    }
    return removed;
}



static void symmetrizeStringGraphLikeHifiasm(Assembler& assembler, uint32_t vertexCount)
{
    // Hifiasm `asg_symm` = delete multi-arcs then delete asymmetric arcs.
    (void)deleteStringGraphMultiArcsLikeHifiasm(assembler, vertexCount);
    (void)deleteStringGraphAsymmetricArcsLikeHifiasm(assembler, vertexCount);
}



static uint64_t stringGraphTransitiveReduce(
    Assembler& assembler,
    uint32_t fuzz,
    uint32_t vertexCount,
    bool rebuildBefore)
{
    StringGraph& g = assembler.stringGraph;

    // Ensure adjacency exists and is sorted by `len`.
    if (rebuildBefore) {
        rebuildStringGraphAdjacency(assembler, vertexCount);
    } else {
        sortStringGraphOutgoingByLen(assembler, vertexCount);
    }

    vector<uint8_t> mark(vertexCount, 0);
    uint64_t reduced = 0;

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) {
            for (const uint32_t arcId : g.outgoing[v]) {
                ++reduced;
                g.arcs[arcId].del = 1;
            }
            continue;
        }

        const span<const uint32_t> out = g.outgoing[v];
        if (out.empty()) continue;

        for (const uint32_t arcId : out) {
            mark[g.arcs[arcId].to] = 1;
        }

        const uint32_t L = g.arcs[out.back()].len + fuzz;
        for (const uint32_t arcIdVw : out) {
            const uint32_t w = g.arcs[arcIdVw].to;
            if (mark[w] != 1) continue;

            const span<const uint32_t> outW = g.outgoing[w];
            for (const uint32_t arcIdWx : outW) {
                if (g.arcs[arcIdWx].len + g.arcs[arcIdVw].len > L) {
                    break;
                }
                const uint32_t x = g.arcs[arcIdWx].to;
                if (mark[x]) {
                    mark[x] = 2;
                }
            }
        }

        for (const uint32_t arcId : out) {
            const uint32_t w = g.arcs[arcId].to;
            if (mark[w] == 2) {
                g.arcs[arcId].del = 1;
                ++reduced;
            }
            mark[w] = 0;
        }
    }

    if (reduced) {
        // Hifiasm `asg_arc_del_trans` post-processing:
        // `asg_cleanup(g); asg_symm(g);`
        cleanupStringGraphLikeHifiasm(assembler, vertexCount);
        symmetrizeStringGraphLikeHifiasm(assembler, vertexCount);
    }
    return reduced;
}



static uint64_t stringGraphCutTips(Assembler& assembler, uint32_t maxShortTipReads, uint32_t vertexCount)
{
    // -----------------------------------------------------------------------
    // Port of hifiasm `asg_arc_cut_tips` (gfa_ut.cpp), non-ou/non-telomere path.
    //
    // A "tip" is a short dead-end chain that starts at a vertex with no
    // incoming arcs and terminates before reaching maxShortTipReads reads:
    //
    //   tip[0] ──► tip[1] ──► ... ──► tip[k-1] ──► (branch or dead-end)
    //      ▲
    //   no incoming arcs (tip[0]^1 has out-degree 0)
    //
    // A tip is "short" when the chain ends (non-MERGEABLE) within
    // maxShortTipReads extension steps.  Long tips are left intact.
    //
    // Algorithm (two-phase, matching hifiasm):
    //   Phase 1 — scan all vertices; for each tip start, walk the chain and
    //             record short tips as candidates keyed by
    //             (pathLength << 32 | tipStart) for length-ascending sort.
    //   Phase 2 — process candidates shortest-first; re-check each one
    //             because earlier deletions in this phase may have invalidated
    //             a candidate, then delete every read in the tip chain.
    //
    // Divergence from hifiasm:
    //   `asg_arc_cut_tips` also handles `is_ou` (unitig-overlap extension)
    //   and `te` (telomere-anchored tip protection).  Dinara's string graph
    //   encodes neither, so only the base tip-cutting path is implemented.
    // -----------------------------------------------------------------------

    StringGraph& g = assembler.stringGraph;

    if (maxShortTipReads == 0) return 0;

    // -----------------------------------------------------------------------
    // isTipStart — true when oriented vertex v has no incoming arcs.
    // Equivalently: v^1 has no non-deleted outgoing arcs.
    // Hifiasm: asg_arc_a(g, v^1) with nv == 0.
    // -----------------------------------------------------------------------
    auto isTipStart = [&](uint32_t v) -> bool {
        return countOutgoingNonDeleted(g, v ^ 1U) == 0;
    };

    // -----------------------------------------------------------------------
    // deleteRead — mark the read and all incident arcs (plus twins) deleted.
    // Hifiasm equivalent: asg_seq_del, with arc deletion deferred to asg_cleanup.
    // We mark twins eagerly so that re-checks in Phase 2 see a consistent graph.
    // -----------------------------------------------------------------------
    auto deleteRead = [&](uint32_t orientedVertex) {
        const ReadId readId = ReadId(orientedVertex >> 1U);
        if (g.readDeleted.isOpen && readId < g.readDeleted.size()) {
            g.readDeleted[readId] = 1;
        }
        const uint32_t v0 = uint32_t(readId) << 1U;
        const uint32_t v1 = v0 ^ 1U;
        for (const uint32_t arcId : g.outgoing[v0]) {
            g.arcs[arcId      ].del = 1;   // arc FROM v0
            g.arcs[arcId ^ 1U ].del = 1;   // twin arc INTO v0
        }
        for (const uint32_t arcId : g.outgoing[v1]) {
            g.arcs[arcId      ].del = 1;   // arc FROM v1
            g.arcs[arcId ^ 1U ].del = 1;   // twin arc INTO v1
        }
    };

    // -----------------------------------------------------------------------
    // Phase 1: Collect candidate tip-start vertices.
    //
    // Each candidate is packed as (pathLength << 32 | tipStart) so that
    // sorting by value yields shortest-tip-first order for Phase 2.
    // -----------------------------------------------------------------------
    vector<uint64_t> candidates;
    candidates.reserve(vertexCount / 8 + 16);

    for (uint32_t tipStart = 0; tipStart < vertexCount; ++tipStart) {
        if (g.readDeleted.isOpen && g.readDeleted[tipStart >> 1U]) continue;
        if (!isTipStart(tipStart)) continue;  // v^1 has outgoing arcs → not a tip

        // Walk the chain from tipStart, extending up to maxShortTipReads steps.
        // Each step follows the unique MERGEABLE continuation via stringGraphIsUtgEnd.
        uint32_t cursor   = tipStart;
        uint32_t pathLen  = 1;   // reads in chain including tipStart
        uint32_t extSteps = 0;   // successful extension steps taken
        for (; extSteps < maxShortTipReads; ++extSteps) {
            uint32_t nextVertex = 0;
            if (stringGraphIsUtgEnd(g, cursor ^ 1U, &nextVertex) != ASG_ET_MERGEABLE) break;
            cursor = nextVertex;
            ++pathLen;
        }

        // Short tip: the chain ended before maxShortTipReads steps.
        // (extSteps == maxShortTipReads means the tip is long enough to keep.)
        if (extSteps < maxShortTipReads) {
            candidates.push_back((uint64_t(pathLen) << 32) | uint64_t(tipStart));
        }
    }

    // Sort ascending: shorter tips are processed first so that deleting a
    // short tip can never inadvertently extend a longer candidate.
    std::sort(candidates.begin(), candidates.end());

    // -----------------------------------------------------------------------
    // Phase 2: Validate and cut tips, shortest first.
    //
    // Re-check each candidate: earlier deletions in this phase may have
    // changed the graph so that a previously identified tip is no longer
    // valid or is no longer short.
    // -----------------------------------------------------------------------
    uint64_t cutCount = 0;
    vector<uint32_t> tipPath;
    tipPath.reserve(maxShortTipReads + 2);

    for (const uint64_t candidate : candidates) {
        const uint32_t tipStart = uint32_t(candidate);  // lower 32 bits

        // Skip if this read was already deleted by an earlier tip cut.
        if (g.readDeleted.isOpen && g.readDeleted[tipStart >> 1U]) continue;

        // Re-verify: is it still a tip start (no incoming arcs)?
        if (!isTipStart(tipStart)) continue;

        // Re-walk the chain from this (still-valid) tip start.
        tipPath.clear();
        tipPath.push_back(tipStart);

        uint32_t cursor   = tipStart;
        uint32_t extSteps = 0;
        for (; extSteps < maxShortTipReads; ++extSteps) {
            uint32_t nextVertex = 0;
            if (stringGraphIsUtgEnd(g, cursor ^ 1U, &nextVertex) != ASG_ET_MERGEABLE) break;
            cursor = nextVertex;
            tipPath.push_back(cursor);
        }

        // Delete the tip if it is still short (may have shrunk since Phase 1).
        if (extSteps < maxShortTipReads) {
            for (const uint32_t v : tipPath) {
                deleteRead(v);
            }
            ++cutCount;
        }
    }

    // Propagate arc deletions to twins and rebuild adjacency (≡ asg_cleanup).
    if (cutCount > 0) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(assembler, vertexCount);
    }

    return cutCount;
}

static uint64_t stringGraphRemoveSingleNodeBubbles(Assembler& assembler, uint32_t maxShortTipReads, uint32_t vertexCount)
{
    // -----------------------------------------------------------------------
    // Port of hifiasm `asg_arc_del_single_node_directly` (called with
    // longLen_thres = asm_opt.max_short_tip = maxShortTipReads), which
    // delegates each candidate vertex to `test_single_node_bubble_directly`.
    //
    // Removes the weaker arm of a symmetric single-node bubble:
    //
    //         ┌─ arcVa ─► a ─► t
    //   v ───►┤                ▲
    //         └─ arcVb ─► b ──┘
    //
    // Conditions (corresponding to the Len[0]==1 && Len[1]==1 branch of
    // test_single_node_bubble_directly after check_small_bubble):
    //   • v has exactly 2 non-deleted outgoing arcs (to a and b).
    //   • a and b each have exactly 1 incoming arc, and it comes from v.
    //   • a and b each have exactly 1 outgoing arc, and both lead to the
    //     same merge vertex t.
    //
    // Deletion criterion divergence from hifiasm:
    //   Hifiasm uses `el` (exact-overlap) flags and `is_abnormal` status to
    //   decide which arm to delete: it only deletes when one arm has a
    //   non-exact overlap (el==0) or is flagged abnormal, and it skips the
    //   bubble entirely when both arms are fully exact.
    //   Dinara has no stored `el` or `is_abnormal` equivalents, so it
    //   always deletes the arm with lower total overlap length
    //   (overlapLen(v→arm) + overlapLen(arm→t)).  This makes dinara more
    //   aggressive: it may delete arms that hifiasm would protect.
    //
    // After all deletions, asg_cleanup is replicated by
    // symmetrizeArcDeletion + rebuildStringGraphAdjacency.
    // -----------------------------------------------------------------------

    StringGraph& g = assembler.stringGraph;
    if (!g.incoming.isOpen()) {
        return 0;
    }

    // -----------------------------------------------------------------------
    // inDegree — count non-deleted arcs pointing INTO vertex v.
    // -----------------------------------------------------------------------
    auto inDegree = [&](uint32_t v) -> uint32_t {
        uint32_t n = 0;
        for (const uint32_t arcId : g.incoming[v]) {
            if (!g.arcs[arcId].del) ++n;
        }
        return n;
    };

    // -----------------------------------------------------------------------
    // hasIncomingFrom — true if there is a non-deleted arc from `src` to `v`.
    // -----------------------------------------------------------------------
    auto hasIncomingFrom = [&](uint32_t v, uint32_t src) -> bool {
        for (const uint32_t arcId : g.incoming[v]) {
            const auto& a = g.arcs[arcId];
            if (!a.del && a.from == src) return true;
        }
        return false;
    };

    // -----------------------------------------------------------------------
    // singleOutgoingArc — true iff v has exactly one non-deleted outgoing arc.
    // Sets arcIdOut and target to that arc's id and destination.
    // Corresponds to hifiasm nv==1 check inside check_small_bubble.
    // -----------------------------------------------------------------------
    auto singleOutgoingArc = [&](uint32_t v, uint32_t& arcIdOut, uint32_t& target) -> bool {
        bool found = false;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (g.arcs[arcId].del) continue;
            if (found) return false;  // second non-deleted arc → not unique
            arcIdOut = arcId;
            target   = g.arcs[arcId].to;
            found    = true;
        }
        return found;
    };

    uint64_t deletedArcCount = 0;

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) continue;

        // Collect the exactly 2 non-deleted outgoing arcs of v.
        // (hifiasm: asg_arc_n(g,v) == 2, equivalent post-cleanup)
        array<uint32_t, 2> outArcIds{};
        uint32_t outCount = 0;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (g.arcs[arcId].del) continue;
            if (outCount < 2) outArcIds[outCount] = arcId;
            ++outCount;
        }
        if (outCount != 2) continue;

        const uint32_t arcVa = outArcIds[0];
        const uint32_t arcVb = outArcIds[1];
        const uint32_t a     = g.arcs[arcVa].to;  // arm 0
        const uint32_t b     = g.arcs[arcVb].to;  // arm 1
        if (a == b) continue;  // degenerate: both arcs go to the same vertex

        // Each arm must have exactly 1 incoming arc, coming from v.
        // (hifiasm: asg_is_single_edge(g, av[i].v, v>>1) == 1)
        if (inDegree(a) != 1 || !hasIncomingFrom(a, v)) continue;
        if (inDegree(b) != 1 || !hasIncomingFrom(b, v)) continue;

        // Each arm must have exactly 1 outgoing arc and both must lead to
        // the same merge vertex t.  (hifiasm: nv==1 && nw==1 in check_small_bubble,
        // plus path-length computation confirming both paths reach the same endNode.)
        uint32_t arcAtoT = 0, tA = 0;
        uint32_t arcBtoT = 0, tB = 0;
        if (!singleOutgoingArc(a, arcAtoT, tA)) continue;
        if (!singleOutgoingArc(b, arcBtoT, tB)) continue;
        if (tA != tB) continue;

        // The merge vertex t must have exactly 2 incoming arcs (one from each arm).
        // (hifiasm: asg_is_single_edge(g, arm_target, arm>>1) == 2 inside check_small_bubble)
        if (inDegree(tA) != 2) continue;

        // -----------------------------------------------------------------------
        // Deletion criterion: remove the arm with the lower total overlap length.
        // Hifiasm instead uses el-flags / is_abnormal — see comment at top.
        // -----------------------------------------------------------------------
        const uint32_t scoreA = g.arcs[arcVa].overlapLen + g.arcs[arcAtoT].overlapLen;
        const uint32_t scoreB = g.arcs[arcVb].overlapLen + g.arcs[arcBtoT].overlapLen;

        const bool removeB          = (scoreA >= scoreB);
        const uint32_t delArcVx     = removeB ? arcVb  : arcVa;   // v → deleted arm
        const uint32_t delArcXtoT   = removeB ? arcBtoT : arcAtoT; // deleted arm → t

        g.arcs[delArcVx    ].del = 1;
        g.arcs[delArcVx ^ 1U].del = 1;   // twin: t^1 → v^1 side of deleted arm
        g.arcs[delArcXtoT    ].del = 1;
        g.arcs[delArcXtoT ^ 1U].del = 1; // twin: t^1 → arm^1

        deletedArcCount += 4;
    }

    if (deletedArcCount > 0) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(assembler, vertexCount);
    }
    return deletedArcCount;
}

uint64_t Assembler::cleanStringGraphRemoveSingleNodeBubbles(uint32_t maxShortTipReads)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen() || !stringGraph.incoming.isOpen()) {
        throw runtime_error("cleanStringGraphRemoveSingleNodeBubbles: StringGraph arcs/outgoing/incoming must be open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph single-node bubble removal begins." << endl;

    const uint64_t removed = stringGraphRemoveSingleNodeBubbles(*this, maxShortTipReads, vertexCount);

    cout << timestamp << "  removed " << removed << " bubble arcs (maxShortTipReads=" << maxShortTipReads << ")" << endl;
    cout << timestamp << "String graph single-node bubble removal complete in "
         << seconds(steady_clock::now() - t0) << " s" << endl;
    return removed;
}



uint64_t Assembler::reduceStringGraphTransitiveHifiasm(uint32_t gapFuzz)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen) {
        throw runtime_error("reduceStringGraphTransitiveHifiasm: StringGraph is not open.");
    }
    if (!stringGraph.outgoing.isOpen()) {
        throw runtime_error("reduceStringGraphTransitiveHifiasm: StringGraph adjacency is not open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph transitive reduction (hifiasm asg_arc_del_trans) begins." << endl;

    // This entry point is used right after ma_sg_gen, so adjacency is already current.
    const uint64_t reduced = stringGraphTransitiveReduce(*this, gapFuzz, vertexCount, /*rebuildBefore*/false);
    cout << timestamp << "  transitively reduced " << reduced << " arcs (gapFuzz=" << gapFuzz << ")" << endl;
    cout << timestamp << "String graph transitive reduction complete in " << seconds(steady_clock::now() - t0) << " s" << endl;
    return reduced;
}



uint64_t Assembler::cutStringGraphTips(uint32_t maxShortTipReads)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen) {
        throw runtime_error("cutStringGraphTips: StringGraph is not open.");
    }
    if (!stringGraph.outgoing.isOpen()) {
        throw runtime_error("cutStringGraphTips: StringGraph adjacency is not open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph tip cutting (hifiasm asg_arc_cut_tips) begins." << endl;

    const uint64_t cut = stringGraphCutTips(*this, maxShortTipReads, vertexCount);
    cout << timestamp << "  cut " << cut << " tips (maxShortTipReads=" << maxShortTipReads << ")" << endl;
    cout << timestamp << "String graph tip cutting complete in " << seconds(steady_clock::now() - t0) << " s" << endl;
    return cut;
}

uint64_t Assembler::cutStringGraphWeakArcsOntHifiasm(uint32_t maxExtReads, double lenRatio, uint32_t minDiff)
{
    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();
    if (!stringGraph.arcs.isOpen) {
        throw runtime_error("cutStringGraphWeakArcsOntHifiasm: StringGraph is not open.");
    }
    if (!stringGraph.outgoing.isOpen()) {
        throw runtime_error("cutStringGraphWeakArcsOntHifiasm: StringGraph adjacency is not open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph weak-arc cutting (hifiasm asg_arc_cut_weak, ONT path) begins." << endl;

    StringGraph& g = stringGraph;
    if (maxExtReads == 0) return 0;

    auto hasAlignmentBetweenOrientedReads = [&](OrientedReadId a, OrientedReadId b) -> bool {
        const span<const uint32_t> section = alignmentTable[a.getValue()];
        auto it = std::lower_bound(
            section.begin(),
            section.end(),
            b,
            [&](uint32_t alignmentId, const OrientedReadId& target) {
                const OrientedReadId other = alignmentData[alignmentId].getOther(a);
                return other < target;
            });
        if (it == section.end()) return false;
        return alignmentData[*it].getOther(a) == b;
    };

    auto hasAnyAlignmentBetweenReads = [&](ReadId a, ReadId b) -> bool {
        for (Strand sa = 0; sa < 2; ++sa) {
            for (Strand sb = 0; sb < 2; ++sb) {
                if (hasAlignmentBetweenOrientedReads(OrientedReadId(a, sa), OrientedReadId(b, sb))) {
                    return true;
                }
            }
        }
        return false;
    };

    auto readIsDeleted = [&](uint32_t v) -> bool {
        return g.readDeleted.isOpen && (v >> 1U) < g.readDeleted.size() && g.readDeleted[v >> 1U];
    };

    // Approximate parity with hifiasm gfa_ut.cpp:asg_arc_cut_weak for the ONT call site used by ul_clean_gfa:
    //   asg_arc_cut_weak(sg, &bu, max_tip, 0.975, 0, is_ou, 0, 1, 16, UL_COV_THRES-1, 0, rev, NULL, NULL)
    // We intentionally do not implement trio/OU/R_to_U paths here (Dinara does not store `ou`).

    vector<uint64_t> candidates;
    candidates.reserve(1024);

    // Phase 1: gather weak-arc candidates.
    // We use "not the best-overlap outgoing arc" as the weak-arc proxy, and require that the best target
    // overlaps the weak target (dedup check), matching the intent of hifiasm's `is_dedup_weak_arc`.
    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (readIsDeleted(v)) continue;
        const span<const uint32_t> out = g.outgoing[v];

        uint32_t maxOl = 0;
        uint32_t outCount = 0;
        uint32_t maxArcId = std::numeric_limits<uint32_t>::max();
        for (const uint32_t arcId : out) {
            const auto& a = g.arcs[arcId];
            if (a.del) continue;
            ++outCount;
            if (a.overlapLen > maxOl) {
                maxOl = a.overlapLen;
                maxArcId = arcId;
            }
        }
        if (outCount < 2) continue;
        if (maxArcId == std::numeric_limits<uint32_t>::max()) continue;

        const ReadId strongRead = ReadId(g.arcs[maxArcId].to >> 1U);
        for (const uint32_t arcId : out) {
            const auto& a = g.arcs[arcId];
            if (a.del) continue;
            if (arcId == maxArcId) continue;

            if (a.overlapLen >= maxOl) continue;
            if (uint32_t(double(a.overlapLen)) + minDiff > maxOl) continue;
            if (double(a.overlapLen) > double(maxOl) * lenRatio) continue;

            const ReadId weakRead = ReadId(a.to >> 1U);
            if (!hasAnyAlignmentBetweenReads(strongRead, weakRead)) continue;

            candidates.push_back((uint64_t(a.overlapLen) << 32) | uint64_t(arcId));
        }
    }

    std::sort(candidates.begin(), candidates.end());

    // Phase 2: evaluate candidates in increasing overlap length and delete if supported by local criteria.
    uint64_t cut = 0;
    for (const uint64_t key : candidates) {
        const uint32_t arcId = uint32_t(key);
        if (arcId >= g.arcs.size()) continue;
        if (g.arcs[arcId].del) continue;

        const uint32_t v = g.arcs[arcId].from;
        const uint32_t to = g.arcs[arcId].to;
        const uint32_t w = to ^ 1U; // hifiasm: w = ve->v^1
        if (readIsDeleted(v) || readIsDeleted(w)) continue;

        const uint32_t kv = countOutgoingNonDeleted(g, v);
        const uint32_t kw = countOutgoingNonDeleted(g, w);
        if (kv <= 1 && kw <= 1) continue;

        const uint32_t twinId = arcId ^ 1U;
        if (twinId >= g.arcs.size()) continue;
        const auto& ve = g.arcs[arcId];
        const auto& we = g.arcs[twinId];

        const uint32_t mmOl = std::min(ve.overlapLen, we.overlapLen);

        // Find best alternative on v side that overlaps the weak target read.
        uint32_t vOlMax = 0;
        for (const uint32_t altId : g.outgoing[v]) {
            const auto& a = g.arcs[altId];
            if (a.del) continue;
            if (a.to == ve.to) continue;
            if (a.overlapLen <= ve.overlapLen) continue;
            if (!hasAnyAlignmentBetweenReads(ReadId(a.to >> 1U), ReadId(ve.to >> 1U))) continue;
            vOlMax = std::max(vOlMax, a.overlapLen);
        }
        if (kv >= 2) {
            if (vOlMax == 0) continue;
            if (double(mmOl) > double(vOlMax) * lenRatio) continue;
            if (mmOl + minDiff > vOlMax) continue;
        }

        // Find best alternative on w side that overlaps the weak source read (mirrors symmetric checks).
        uint32_t wOlMax = 0;
        for (const uint32_t altId : g.outgoing[w]) {
            const auto& a = g.arcs[altId];
            if (a.del) continue;
            if (a.to == we.to) continue;
            if (a.overlapLen <= we.overlapLen) continue;
            if (!hasAnyAlignmentBetweenReads(ReadId(a.to >> 1U), ReadId(we.to >> 1U))) continue;
            wOlMax = std::max(wOlMax, a.overlapLen);
        }
        if (kw >= 2) {
            if (wOlMax == 0) continue;
            if (double(mmOl) > double(wOlMax) * lenRatio) continue;
            if (mmOl + minDiff > wOlMax) continue;
        }

        // Topology criterion (hifiasm `is_topo=1` path).
        bool toDel = false;
        if (kv > 1 && kw > 1) {
            toDel = true;
        } else if (kw == 1) {
            if (stringGraphTopocutAux(g, w ^ 1U, int(maxExtReads)) < int(maxExtReads)) toDel = true;
        } else if (kv == 1) {
            if (stringGraphTopocutAux(g, v ^ 1U, int(maxExtReads)) < int(maxExtReads)) toDel = true;
        }

        if (toDel) {
            g.arcs[arcId].del = 1;
            g.arcs[twinId].del = 1;
            ++cut;
        }
    }

    if (cut) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(*this, vertexCount);
    }

    cout << timestamp << "  cut " << cut << " weak arcs (maxExtReads=" << maxExtReads
         << ", lenRatio=" << lenRatio << ", minDiff=" << minDiff << ")" << endl;
    cout << timestamp << "String graph weak-arc cutting complete in " << seconds(steady_clock::now() - t0) << " s" << endl;
    return cut;
}



void Assembler::cleanStringGraphInitialHifiasm(uint32_t gapFuzz, uint32_t maxShortTipReads)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen) {
        throw runtime_error("cleanStringGraphInitialHifiasm: StringGraph is not open.");
    }
    if (!stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphInitialHifiasm: StringGraph adjacency is not open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());

    const auto t0 = steady_clock::now();
    cout << timestamp << "Cleaning string graph (hifiasm initial clean): transitive reduction + tip cut" << endl;

    const uint64_t reduced = stringGraphTransitiveReduce(*this, gapFuzz, vertexCount, /*rebuildBefore*/true);
    cout << timestamp << "  transitively reduced " << reduced << " arcs (gapFuzz=" << gapFuzz << ")" << endl;

    if (maxShortTipReads > 0) {
        const uint64_t cut = stringGraphCutTips(*this, maxShortTipReads, vertexCount);
        cout << timestamp << "  cut " << cut << " tips (maxShortTipReads=" << maxShortTipReads << ")" << endl;
    }

    cout << timestamp << "String graph clean complete in " << seconds(steady_clock::now() - t0) << " s" << endl;
}



void Assembler::cleanStringGraphIterativeHifiasm(
    uint32_t cleanRounds,
    double minDropRate,
    double maxDropRate,
    uint32_t maxShortTipReads)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphIterativeHifiasm: StringGraph is not open.");
    }

    if (cleanRounds == 0) return;
    if (minDropRate <= 0. || maxDropRate <= 0. || minDropRate > maxDropRate || maxDropRate >= 1.) {
        throw runtime_error("cleanStringGraphIterativeHifiasm: invalid drop-rate range.");
    }

    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph iterative cleaning (hifiasm ul_clean_gfa) begins: "
         << cleanRounds << " rounds, drop ratio " << minDropRate << " → " << maxDropRate << endl;

    const double step = (cleanRounds == 1) ? 0. : (maxDropRate - minDropRate) / double(cleanRounds - 1);
    double dropRatio = minDropRate;

    for (uint32_t round = 0; round < cleanRounds; ++round, dropRatio += step) {
        if (dropRatio > maxDropRate) dropRatio = maxDropRate;

        cout << timestamp << "  Round " << (round + 1) << "/" << cleanRounds
             << " (dropRatio=" << dropRatio << ")" << endl;

        // Hifiasm ul_clean_gfa round structure:
        // 1. Semi-circular removal (approximated by cycle breaking)
        // 2. Bubble identification + chimeric cutting + tips
        // 3. Bubble identification + inexact overlap removal + tips
        // 4. Bubble identification + overlap-length-ratio cut + tips
        // 5. Bubble identification + bubble links + complex bubble links + tips

        // Phase 1: Topology cleaning (cycles, bubbles, tips)
        cleanStringGraphPreCleanHifiasm(maxShortTipReads);

        // Phase 2: Overlap-based cleaning (drop short overlaps by ratio)
        cleanStringGraphDropShortOverlaps(dropRatio, /*minOverlapLen*/0, maxShortTipReads);

        // Phase 3: Re-clean topology after overlap changes
        cleanStringGraphPreCleanHifiasm(maxShortTipReads);
    }

    cout << timestamp << "String graph iterative cleaning complete in "
         << seconds(steady_clock::now() - t0) << " s" << endl;
}



void Assembler::cleanStringGraphPreCleanHifiasm(uint32_t maxShortTipReads)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphPreCleanHifiasm: StringGraph is not open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    const auto t0 = steady_clock::now();
    cout << timestamp << "    String graph pre-clean (cycles + bubbles + tips) begins." << endl;

    uint64_t totalBubbleArcRemoved = 0;
    uint64_t totalTipsCut = 0;
    uint64_t totalCyclesBroken = 0;
    for (uint32_t iter = 0; iter < 10; ++iter) {
        // Break semi-circular edges and cut chimeric bubbles (hifiasm asg_iterative_semi_circ; te omitted).
        const uint64_t cyclesBroken = cleanStringGraphBreakShortCycles(/*limLen*/100, /*normalLen*/maxShortTipReads);
        totalCyclesBroken += cyclesBroken;

        const uint64_t bubbleRemoved = stringGraphRemoveSingleNodeBubbles(*this, maxShortTipReads, vertexCount);
        totalBubbleArcRemoved += bubbleRemoved;

        const uint64_t tipsCut = stringGraphCutTips(*this, maxShortTipReads, vertexCount);
        totalTipsCut += tipsCut;

        if (cyclesBroken == 0 && bubbleRemoved == 0 && tipsCut == 0) {
            break;
        }
    }

    cout << timestamp << "      broke " << totalCyclesBroken << " cycles, removed "
         << totalBubbleArcRemoved << " bubble arcs, cut " << totalTipsCut
         << " tips (" << seconds(steady_clock::now() - t0) << " s)" << endl;
}



void Assembler::cleanStringGraphDropShortOverlaps(double lenRatio, uint32_t minOverlapLen, uint32_t maxShortTipReads)
{
    // -----------------------------------------------------------------------
    // Port of hifiasm `asg_arc_cut_length` (gfa_ut.cpp) for the dinara case:
    //   is_ou=0, is_trio=0, is_topo=1, min_diff=0, rev=NULL, rI=NULL.
    //
    // Unlike the simpler `asg_arc_del_short` (per-vertex ratio filter), this
    // function is arc-centric: it evaluates each candidate arc against BOTH
    // of its endpoints before deciding to delete it.
    //
    // A candidate arc (v → dest, twin: dest^1 → v^1) is deleted when:
    //
    //   mm_ol = min(arc.overlapLen, twin.overlapLen)   // min of both directions
    //
    //   V-side: mm_ol ≤ olMaxV × lenRatio   (if v has ≥ 2 arcs)
    //   W-side: mm_ol ≤ olMaxW × lenRatio   (if w = dest^1 has ≥ 2 arcs)
    //   Bridge guard: NOT (kv ≤ 1 AND kw ≤ 1)         // skip bridge arcs
    //
    //   Topological guard (asg_topocut_aux / stringGraphTopocutAux):
    //     kv > 1 AND kw > 1  → delete (both endpoints have alternatives)
    //     kw == 1             → delete only if w is on a short linear chain
    //     kv == 1             → delete only if v is on a short linear chain
    //
    // The two-phase structure (collect → sort ascending → process) matches
    // hifiasm: weakest arcs are evaluated first so that early deletions update
    // the in-degree counts seen by later candidates in the same pass.
    //
    // Phase 1 skips bubble-marked vertices (seq_vis != 0 in hifiasm).
    // seqVis is computed by stringGraphMarkBubbleVertices, a port of
    // asg_arc_identify_simple_bubbles_multi + asg_bub_pop1_label.
    // Note: cross-node marking (seqVis=2, via check_if_cross) is omitted.
    //
    // Dinara-specific:
    //   `minOverlapLen > 0` bypasses the ratio check for arcs with
    //   mm_ol < minOverlapLen (they are always drop candidates, still subject
    //   to the topological guard).  This is not present in hifiasm and is
    //   only used by the final-pass call in cleanStringGraphDropOverlapRoundsHifiasm.
    //
    // After deletions, symmetrizeArcDeletion + rebuildStringGraphAdjacency
    // replicate hifiasm's asg_cleanup.
    // -----------------------------------------------------------------------

    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphDropShortOverlaps: StringGraph is not open.");
    }
    if (lenRatio <= 0. || lenRatio >= 1.) {
        throw runtime_error("cleanStringGraphDropShortOverlaps: lenRatio must be in (0,1).");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    const int      maxExt      = int(maxShortTipReads);  // asg_topocut_aux max_ext
    StringGraph&   g           = stringGraph;

    // Ensure adjacency lists are up to date before iterating.
    rebuildStringGraphAdjacency(*this, vertexCount);

    // Compute bubble-vertex markers (≡ hifiasm seq_vis set by
    // asg_arc_identify_simple_bubbles_multi before asg_arc_cut_length).
    // seqVis[v] = 1: bubble source, interior node, or merge-end reverse.
    // These vertices are skipped in Phase 1 so their arcs are not nominated
    // as short-overlap deletion candidates.
    const vector<uint8_t> seqVis = stringGraphMarkBubbleVertices(g, vertexCount);

    // -----------------------------------------------------------------------
    // Phase 1: Collect candidate arcs.
    //
    // For each non-deleted, non-bubble vertex with ≥ 2 non-deleted outgoing
    // arcs, push every non-deleted outgoing arc as (overlapLen << 32 | arcId).
    // Sorting ascending by this key gives weakest-arc-first order in Phase 2.
    // (hifiasm: radix_sort_srt64 ascending on (ol << 32 | arcIndex))
    // -----------------------------------------------------------------------
    vector<uint64_t> candidates;
    candidates.reserve(g.arcs.size() / 4 + 16);

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) continue;
        if (seqVis[v] != 0) continue;  // skip bubble nodes (≡ hifiasm seq_vis != 0 check)

        uint32_t kv = 0;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) ++kv;
        }
        if (kv < 2) continue;

        for (const uint32_t arcId : g.outgoing[v]) {
            if (g.arcs[arcId].del) continue;
            candidates.push_back((uint64_t(g.arcs[arcId].overlapLen) << 32) | uint64_t(arcId));
        }
    }

    std::sort(candidates.begin(), candidates.end());  // ascending overlap = weakest first

    // -----------------------------------------------------------------------
    // Phase 2: Evaluate and delete candidate arcs, weakest first.
    // -----------------------------------------------------------------------
    uint64_t deletedArcCount = 0;

    for (const uint64_t cand : candidates) {
        const uint32_t arcId  = uint32_t(cand);
        const uint32_t twinId = arcId ^ 1U;

        if (g.arcs[arcId].del) continue;  // already deleted by an earlier candidate

        const uint32_t v = g.arcs[arcId].from;       // source of the forward arc
        const uint32_t w = g.arcs[arcId].to ^ 1U;    // dest^1 = source of twin arc

        if (g.readDeleted.isOpen &&
            (g.readDeleted[v >> 1U] || g.readDeleted[w >> 1U])) continue;

        // mm_ol: the minimum overlap length from both directions of this arc.
        // (hifiasm: mm_ol = MIN(ve->ol, we->ol))
        const uint32_t mmOl = min(g.arcs[arcId ].overlapLen,
                                  g.arcs[twinId].overlapLen);

        // Ratio check helper: true when mm_ol passes the lenRatio filter
        // (meaning the arc is "close enough" to the best at that endpoint
        // and should be KEPT).  Also returns false (→ candidate survives the
        // check) when minOverlapLen > 0 and mm_ol is below the absolute floor.
        const bool belowFloor = (minOverlapLen > 0 && mmOl < minOverlapLen);

        // --- V-side ---
        uint32_t kv = 0, olMaxV = 0;
        for (const uint32_t aId : g.outgoing[v]) {
            if (g.arcs[aId].del) continue;
            ++kv;
            if (g.arcs[aId].overlapLen > olMaxV) olMaxV = g.arcs[aId].overlapLen;
        }
        if (kv < 1) continue;  // arc became isolated since Phase 1
        if (kv >= 2 && !belowFloor) {
            // Arc is within lenRatio of the best on v's side → keep.
            // Direct float comparison matching hifiasm: if (mm_ol > ol_max * len_rat) continue
            if (double(mmOl) > double(olMaxV) * lenRatio) continue;
            // min_diff=0: the guard (mm_ol + 0 > ol_max) is never true, skipped.
        }

        // --- W-side ---
        uint32_t kw = 0, olMaxW = 0;
        for (const uint32_t aId : g.outgoing[w]) {
            if (g.arcs[aId].del) continue;
            ++kw;
            if (g.arcs[aId].overlapLen > olMaxW) olMaxW = g.arcs[aId].overlapLen;
        }
        if (kw < 1) continue;
        if (kw >= 2 && !belowFloor) {
            if (double(mmOl) > double(olMaxW) * lenRatio) continue;
        }

        // Bridge guard: if this arc is the only one on both sides, keep it.
        // (hifiasm: if (kv <= 1 && kw <= 1) continue)
        if (kv <= 1 && kw <= 1) continue;

        // --- Topological guard (is_topo=1, asg_topocut_aux) ---
        // Prevents disconnecting subgraphs that are not short tips.
        bool toDel = false;
        if (kv > 1 && kw > 1) {
            toDel = true;  // both endpoints have alternatives — safe to delete
        } else if (kw == 1) {
            // w has only this arc; delete only if w is on a short linear chain.
            // (hifiasm: if (asg_topocut_aux(g, w^1, max_ext) < max_ext) to_del=1)
            if (stringGraphTopocutAux(g, w ^ 1U, maxExt) < maxExt) toDel = true;
        } else {
            // kv == 1: v has only this arc.
            if (stringGraphTopocutAux(g, v ^ 1U, maxExt) < maxExt) toDel = true;
        }

        if (toDel) {
            g.arcs[arcId ].del = 1;
            g.arcs[twinId].del = 1;
            deletedArcCount += 2;
        }
    }

    // Propagate deletions and rebuild adjacency (≡ asg_cleanup).
    if (deletedArcCount > 0) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(*this, vertexCount);
    }
}



void Assembler::cleanStringGraphDropOverlapRoundsHifiasm(
    uint32_t cleanRounds,
    double minDropRate,
    double maxDropRate,
    uint32_t maxShortTipReads,
    uint32_t finalMinOverlapLen)
{
    if (cleanRounds == 0) return;
    if (minDropRate <= 0. || maxDropRate <= 0. || minDropRate > maxDropRate || maxDropRate >= 1.) {
        throw runtime_error("cleanStringGraphDropOverlapRoundsHifiasm: invalid drop-rate range.");
    }

    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph overlap-drop rounds begin (cleanRounds=" << cleanRounds
         << ", minDropRate=" << minDropRate << ", maxDropRate=" << maxDropRate << ")." << endl;

    const double step = (cleanRounds == 1) ? 0. : (maxDropRate - minDropRate) / double(cleanRounds - 1);
    double drop = minDropRate;
    for (uint32_t i = 0; i < cleanRounds; ++i, drop += step) {
        if (drop > maxDropRate) drop = maxDropRate;
        cleanStringGraphDropShortOverlaps(drop, /*minOverlapLen*/0, maxShortTipReads);
        cleanStringGraphPreCleanHifiasm(maxShortTipReads);
    }

    // Final pass similar in spirit to `asg_arc_del_too_short_overlaps(sg, 2000, min_drop_rate, ...)`.
    if (finalMinOverlapLen > 0) {
        cleanStringGraphDropShortOverlaps(minDropRate, finalMinOverlapLen, maxShortTipReads);
        cleanStringGraphPreCleanHifiasm(maxShortTipReads);
    }

    cout << timestamp << "String graph overlap-drop rounds complete in " << seconds(steady_clock::now() - t0) << " s" << endl;
}



	uint64_t Assembler::cleanStringGraphBreakShortCycles(uint32_t limLen)
	{
	    // Backwards-compatible wrapper (hifiasm default max_tip is typically 3).
	    return cleanStringGraphBreakShortCycles(limLen, /*normalLen*/3);
	}

	uint64_t Assembler::cleanStringGraphBreakShortCycles(uint32_t limLen, uint32_t normalLen)
	{
	    reads->checkReadsAreOpen();
	    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen() || !stringGraph.incoming.isOpen()) {
	        throw runtime_error("cleanStringGraphBreakShortCycles: StringGraph arcs/outgoing/incoming must be open.");
	    }
	    if (normalLen > 0 && !alignmentTable.isOpen()) {
	        throw runtime_error("cleanStringGraphBreakShortCycles: alignmentTable must be open to run chimeric-bubble cutting.");
	    }

	    // -----------------------------------------------------------------------
	    // Port of hifiasm `asg_iterative_semi_circ` (ul_clean_gfa, pop_chimer=1,
	    // te=NULL — telomere guard not yet implemented).
	    //
	    // Two cleaning passes are alternated until the graph stops changing:
	    //
	    //  1. Semi-circular arc removal  (asg_cut_semi_circ,  LIM_LEN = 100)
	    //     Removes the "back arc" that completes a short semi-circle:
	    //
	    //         v ──────────────────► ... ──► e
	    //         │                             │
	    //         v^1 ◄─ (back arc, deleted) ─ e^1
	    //
	    //     Condition: v has exactly one outgoing arc; v^1 has ≥ 2; the
	    //     forward path from v is shorter than limLen and is not a dead-end
	    //     or loop.
	    //
	    //  2. Chimeric-bubble read deletion  (asg_cut_chimeric_bub, normalLen)
	    //     Deletes a read v that sits in a bubble where a longer alternative
	    //     path exists, indicating v is likely a chimeric artefact:
	    //
	    //       ┌─ arc0 ─► nghbr0 ─► [alt path, >normalLen reads] ─► altEnd ─►┐
	    //  v ──►┤                                                               merge
	    //  v^1 ─┤                                                               ├──► target
	    //       └─ arc1 ────────────────────────────────────────────────────────┘
	    //
	    //     Condition: v has exactly one outgoing arc and one incoming arc;
	    //     both neighbours are branch points (in-degree 2); the read passes
	    //     the overlap-based chimera test; and the alt path converges at the
	    //     same merge vertex as arc1.
	    //
	    // After convergence, asg_cleanup is replicated by symmetrizeArcDeletion
	    // + cleanupStringGraphLikeHifiasm.
	    // -----------------------------------------------------------------------

	    const uint32_t readCount   = uint32_t(reads->readCount());
	    const uint32_t vertexCount = 2 * readCount;  // each read → two oriented vertices
	    StringGraph& g = stringGraph;

	    // Rebuild adjacency so outgoing lists are consistent and sorted by length,
	    // matching hifiasm's pre-conditions for the cleaning passes.
	    rebuildStringGraphAdjacency(*this, vertexCount);

	    // -----------------------------------------------------------------------
	    // Return codes for followLimitPath — identical to hifiasm's constants.
	    // -----------------------------------------------------------------------
	    static constexpr uint32_t LONG_TIPS  = 0;  // path exceeded length limit
	    static constexpr uint32_t TWO_INPUT  = 1;  // vertex has 2 incoming arcs (merge point)
	    static constexpr uint32_t TWO_OUTPUT = 2;  // vertex has 2 outgoing arcs (branch point)
	    static constexpr uint32_t MUL_INPUT  = 3;  // vertex has >2 incoming arcs
	    static constexpr uint32_t MUL_OUTPUT = 4;  // vertex has >2 outgoing arcs
	    static constexpr uint32_t END_TIPS   = 5;  // dead end (0 outgoing arcs)
	    static constexpr uint32_t LOOP       = 7;  // path cycled back to start

	    // -----------------------------------------------------------------------
	    // isVertexDeleted — true if the read owning vertex v is soft-deleted.
	    // -----------------------------------------------------------------------
	    auto isVertexDeleted = [&](uint32_t v) -> bool {
	        return g.readDeleted.isOpen && g.readDeleted[v >> 1U];
	    };

	    // -----------------------------------------------------------------------
	    // getArcs — hifiasm: get_arcs.
	    // Count non-deleted outgoing arcs of vertex v.  If arcIdOut is non-null,
	    // fill up to arcIdOutCount entries with arc indices (useful for fetching
	    // the single arc when the caller already knows degree == 1).
	    // -----------------------------------------------------------------------
	    auto getArcs = [&](uint32_t v, uint32_t* arcIdOut, uint32_t arcIdOutCount) -> uint32_t {
	        uint32_t count  = 0;
	        uint32_t stored = 0;
	        for (const uint32_t arcId : g.outgoing[v]) {
	            if (g.arcs[arcId].del) continue;
	            if (arcIdOut && stored < arcIdOutCount) {
	                arcIdOut[stored++] = arcId;
	            }
	            ++count;
	        }
	        return count;
	    };

	    // -----------------------------------------------------------------------
	    // followLimitPath — hifiasm: follow_limit_path.
	    //
	    // Walk the unique outgoing path from `start`, stepping one arc at a time
	    // as long as the current vertex has exactly one outgoing arc and the next
	    // vertex has exactly one incoming arc.  Stops and returns a code when a
	    // topological event ends the walk.
	    //
	    //   *endVertex — last vertex visited before the stopping condition
	    //   *pathLen   — number of vertices stepped through (including endVertex)
	    // -----------------------------------------------------------------------
	    auto followLimitPath = [&](uint32_t  start,
	                               uint32_t* endVertex,
	                               uint32_t* pathLen,
	                               uint32_t  limit) -> uint32_t
	    {
	        uint32_t v = start;
	        *pathLen = 0;

	        while (true) {
	            ++(*pathLen);
	            uint32_t arcId   = 0;
	            const uint32_t outDeg = getArcs(v, &arcId, 1);
	            *endVertex = v;  // record before advancing, so callers see the last good vertex

	            if (outDeg == 0)          return END_TIPS;    // dead end
	            if (outDeg == 2)          return TWO_OUTPUT;  // fork reached
	            if (outDeg >  2)          return MUL_OUTPUT;
	            if (*pathLen > limit)     return LONG_TIPS;   // budget exceeded

	            // Exactly one outgoing arc: follow it to the next vertex w.
	            const uint32_t w     = g.arcs[arcId].to;
	            // Incoming degree of w = outgoing degree of w^1 (asg convention).
	            const uint32_t inDeg = getArcs(w ^ 1U, nullptr, 0);
	            v = w;

	            if (inDeg == 2) return TWO_INPUT;   // merge point reached
	            if (inDeg >  2) return MUL_INPUT;
	            if (v == start) return LOOP;         // cycled back to origin
	        }
	    };

	    // -----------------------------------------------------------------------
	    // cutSemiCircularOnce — hifiasm: asg_cut_semi_circ(g, LIM_LEN, is_clean=0).
	    //
	    // For each vertex v that has exactly one outgoing arc and more than one
	    // incoming arc, walk the outgoing path.  If the path is short and clean
	    // (not a dead-end or loop), delete the "back arc" from v^1 to endVertex^1
	    // that completes the semi-circle, together with its twin.
	    // -----------------------------------------------------------------------
	    auto cutSemiCircularOnce = [&](uint32_t limLenArg) -> uint64_t {
	        if (limLenArg == 0) return 0;
	        uint64_t deletedCount = 0;

	        for (uint32_t v = 0; v < vertexCount; ++v) {
	            if (isVertexDeleted(v)) continue;

	            // v must have >1 incoming arcs (= outdegree of v^1 > 1).
	            // Raw size is a cheap pre-filter before the O(degree) deleted scan.
	            if (g.outgoing[v ^ 1U].size() <= 1) continue;
	            if (getArcs(v ^ 1U, nullptr, 0) <= 1) continue;

	            // v must have exactly one outgoing arc.
	            if (g.outgoing[v].empty()) continue;
	            if (getArcs(v, nullptr, 0) != 1) continue;

	            // Walk the unique path forward from v.
	            uint32_t endVertex = 0, pathLen = 0;
	            const uint32_t reason = followLimitPath(v, &endVertex, &pathLen, limLenArg);

	            // Skip paths that are too long, loops, or dead-ends — those are not
	            // semi-circles (matching hifiasm's identical filter condition).
	            if (pathLen > limLenArg ||
	                reason == LONG_TIPS  ||
	                reason == LOOP       ||
	                reason == END_TIPS) {
	                continue;
	            }

	            // Delete the back arc  v^1 → endVertex^1  and its twin  endVertex → v.
	            const uint32_t backArcTarget = endVertex ^ 1U;
	            for (const uint32_t arcId : g.outgoing[v ^ 1U]) {
	                if (g.arcs[arcId].del) continue;
	                if (g.arcs[arcId].to == backArcTarget) {
	                    g.arcs[arcId      ].del = 1;  // forward:  v^1 → endVertex^1
	                    g.arcs[arcId ^ 1U ].del = 1;  // twin:     endVertex → v
	                    ++deletedCount;
	                }
	            }
	        }
	        return deletedCount;
	    };

	    // -----------------------------------------------------------------------
	    // deleteRead — mark a read and all arcs incident to it as deleted.
	    //
	    // Hifiasm's asg_seq_del only marks the read; arc removal is deferred to
	    // asg_cleanup after the whole loop.  We delete arcs eagerly so that
	    // getArcs() on neighbours returns the correct degree within the same pass,
	    // preventing false positives in subsequent iterations.
	    //
	    // outgoing[v0] and outgoing[v1] cover arcs FROM the read.  arcId^1 is the
	    // twin arc (TO the read from some other read's reverse vertex), so marking
	    // it here makes the deletion fully symmetric without waiting for the
	    // symmetrizeArcDeletion call at the end.
	    // -----------------------------------------------------------------------
	    auto deleteRead = [&](ReadId readId) {
	        if (g.readDeleted.isOpen && readId < g.readDeleted.size()) {
	            g.readDeleted[readId] = 1;
	        }
	        const uint32_t v0 = uint32_t(readId) << 1U;  // forward oriented vertex
	        const uint32_t v1 = v0 ^ 1U;                  // reverse oriented vertex

	        for (const uint32_t arcId : g.outgoing[v0]) {
	            g.arcs[arcId      ].del = 1;  // arc FROM v0
	            g.arcs[arcId ^ 1U ].del = 1;  // twin arc INTO v0 (from some other read's v^1)
	        }
	        for (const uint32_t arcId : g.outgoing[v1]) {
	            g.arcs[arcId      ].del = 1;  // arc FROM v1
	            g.arcs[arcId ^ 1U ].del = 1;  // twin arc INTO v1 (from some other read's v^1)
	        }
	    };

	    // -----------------------------------------------------------------------
	    // Scratch space and per-read cache for the chimera test below.
	    // -----------------------------------------------------------------------
	    vector<uint64_t> intervalScratch;
	    intervalScratch.reserve(256);
	    // chimericCache: -1 = not yet tested, 0 = not chimeric, 1 = chimeric.
	    vector<int8_t> chimericCache(readCount, int8_t(-1));

	    // -----------------------------------------------------------------------
	    // isExactOverlapForSupChimeric — hifiasm: if_exact / ma_hit_t::el.
	    //
	    // Hifiasm admits only "exact" overlaps (el=1) to the chimera interval
	    // test.  Dinara has no stored el-equivalent, so all overlaps are treated
	    // as exact.  If an el-equivalent is added later, restore the filter below.
	    // -----------------------------------------------------------------------
	    auto isExactOverlap = [&](const AlignmentData& ad) -> bool {
	        (void)ad;
	        return true;
	        // if (ad.hasLargeIndel) return false;
	        // if (ad.info.errorRate != invalid<float> && ad.info.errorRate > 0.01f) return false;
	        // return true;
	    };

	    // -----------------------------------------------------------------------
	    // readLenForChimeraTest — effective read length used for chimera intervals.
	    // Uses the trimmed valid interval when available; falls back to raw length.
	    // -----------------------------------------------------------------------
	    auto readLenForChimeraTest = [&](ReadId readId) -> uint32_t {
	        if (!validReadIntervals.empty() && readId < validReadIntervals.size()) {
	            const auto& vi = validReadIntervals[readId];
	            if (!vi.isDeleted && vi.end >= vi.start) {
	                return vi.end - vi.start;
	            }
	        }
	        return uint32_t(reads->getRead(readId).baseCount);
	    };

	    // -----------------------------------------------------------------------
	    // queryInterval — extract [qs, qe) of readId from one AlignmentData.
	    // Handles the two orientations (readId may be on either side of the
	    // alignment) and clamps qe to rLen.  Returns false for degenerate ranges.
	    // -----------------------------------------------------------------------
	    auto queryInterval = [&](const AlignmentData& ad, ReadId readId,
	                             uint32_t rLen,
	                             uint32_t& qs, uint32_t& qe) -> bool
	    {
	        if (ad.readIds[0] == readId) { qs = ad.qs; qe = ad.qe; }
	        else                         { qs = ad.ts; qe = ad.te; }
	        if (qe <= qs) return false;
	        if (qe > rLen) qe = rLen;
	        return true;
	    };

	    // -----------------------------------------------------------------------
	    // ifSupChimeric — hifiasm: if_sup_chimeric(src, rLen, b, if_exact=1).
	    //
	    // Returns true when the read's overlaps suggest it is chimeric, i.e. there
	    // is an uncovered gap between the region spanned by left-anchored overlaps
	    // (qs == 0) and the region spanned by right-anchored overlaps (qe == rLen).
	    //
	    // Two passes over the overlap set:
	    //
	    //   Pass 1 (fast):
	    //     Compute l1 = max right-end of left-anchored overlaps.
	    //     Compute r0 = min left-start of right-anchored overlaps.
	    //     • l1 > r0 → the two anchored regions already overlap → not chimeric.
	    //     • l1 == 0 or r1 == 0 → one side has no coverage at all → chimeric.
	    //
	    //   Pass 2 (interval depth sweep):
	    //     Seed with the merged anchor intervals, add all internal overlaps
	    //     (qs > 0 and qe < rLen), sort as (position << 1 | isEnd) events,
	    //     and sweep a depth counter to find maximal covered runs.  Record the
	    //     run touching position 0 (→ new l1) and the one touching rLen (→ new
	    //     r0).  If l1 ≤ r0 there is a gap → chimeric.
	    //
	    // Divergence from hifiasm: hifiasm passes the full ma_hit_t_alloc (all
	    // overlaps, filtered only by .del).  Dinara uses alignmentTable (read-graph
	    // overlaps, further filtered by isInReadGraph + keptByBothSides +
	    // isDeletedFromReadPerspective).  Overlaps absent from the read graph are
	    // invisible here, making this test slightly less sensitive than hifiasm's.
	    // -----------------------------------------------------------------------
	    auto ifSupChimeric = [&](ReadId readId) -> bool {
	        if (readId >= readCount) return false;

	        int8_t& cached = chimericCache[readId];
	        if (cached != int8_t(-1)) return cached != 0;

	        const uint32_t rLen = readLenForChimeraTest(readId);

	        // l1: rightmost end reached by left-anchored (qs == 0) overlaps.
	        // r0: leftmost start reached by right-anchored (qe == rLen) overlaps.
	        // r1: rLen when any right-anchored overlap exists, 0 otherwise.
	        uint32_t l1 = 0, r0 = rLen, r1 = 0;

	        const span<const uint32_t> overlaps =
	            alignmentTable[OrientedReadId(readId, 0).getValue()];

	        // Helper: accept only non-deleted, read-graph, exact overlaps.
	        auto accept = [&](const AlignmentData& ad) -> bool {
	            return ad.info.isInReadGraph
	                && ad.keptByBothSides()
	                && !ad.isDeletedFromReadPerspective(readId)
	                && isExactOverlap(ad);
	        };

	        // --- Pass 1: identify left-anchored and right-anchored coverage -------
	        for (const uint32_t idx : overlaps) {
	            const AlignmentData& ad = alignmentData[idx];
	            if (!accept(ad)) continue;

	            uint32_t qs = 0, qe = 0;
	            if (!queryInterval(ad, readId, rLen, qs, qe)) continue;

	            if (qs == 0)    { if (qe > l1) l1 = qe; }               // extends left coverage
	            if (qe == rLen) { if (qs < r0) r0 = qs; r1 = rLen; }   // extends right coverage
	        }

	        // Left and right anchored regions already overlap → not chimeric.
	        if (l1 > r0) { cached = 0; return false; }

	        // One side completely uncovered → definitely chimeric.
	        if (l1 == 0 || r1 == 0) { cached = 1; return true; }

	        // --- Pass 2: depth sweep with internal overlaps ----------------------
	        // Events are encoded as (position << 1) | isClosing so a single sort
	        // puts openings before closings at the same position.
	        const size_t sweepBase = intervalScratch.size();

	        // Seed with the anchor intervals found in Pass 1.
	        // (Both conditions are always true at this point, but we mirror
	        //  hifiasm's guards for exactness.)
	        if (l1 > 0) {
	            intervalScratch.push_back( uint64_t(0)    << 1);           // left anchor open
	            intervalScratch.push_back((uint64_t(l1)   << 1) | 1ULL);  // left anchor close
	        }
	        if (r1 > r0) {
	            intervalScratch.push_back( uint64_t(r0)   << 1);           // right anchor open
	            intervalScratch.push_back((uint64_t(rLen)  << 1) | 1ULL); // right anchor close
	        }

	        // Add internal overlaps (neither left- nor right-anchored).
	        for (const uint32_t idx : overlaps) {
	            const AlignmentData& ad = alignmentData[idx];
	            if (!accept(ad)) continue;

	            uint32_t qs = 0, qe = 0;
	            if (!queryInterval(ad, readId, rLen, qs, qe)) continue;
	            if (qs == 0 || qe == rLen) continue;  // anchored: already seeded above

	            intervalScratch.push_back( uint64_t(qs) << 1);
	            intervalScratch.push_back((uint64_t(qe) << 1) | 1ULL);
	        }

	        std::sort(intervalScratch.begin() + sweepBase, intervalScratch.end());

	        // Sweep: track runs where coverage depth ≥ 1.
	        // Record the run touching position 0 (→ l1) and the one touching rLen (→ r0).
	        l1 = 0; r0 = rLen; r1 = 0;
	        uint32_t runStart = 0;
	        int32_t  depth    = 0;

	        for (size_t i = sweepBase; i < intervalScratch.size(); ++i) {
	            const uint64_t event   = intervalScratch[i];
	            const uint32_t pos     = uint32_t(event >> 1);
	            const int32_t  prevDep = depth;

	            if (event & 1ULL) --depth;  // closing event
	            else               ++depth; // opening event

	            if (prevDep < 1 && depth >= 1) {
	                runStart = pos;              // covered run begins here
	            } else if (prevDep >= 1 && depth < 1) {
	                // Covered run just ended at pos.
	                if (runStart == 0) l1 = pos;              // run anchored at left edge
	                if (pos == rLen)  { r0 = runStart; r1 = pos; } // run anchored at right edge
	            }
	        }

	        intervalScratch.resize(sweepBase);  // release scratch for next call

	        // Chimeric when the left-coverage extent does not reach right-coverage start.
	        const bool isSuspectedChimeric = (l1 <= r0);
	        cached = isSuspectedChimeric ? 1 : 0;
	        return isSuspectedChimeric;
	    };

	    // -----------------------------------------------------------------------
	    // cutChimericBubblesOnce — hifiasm: asg_cut_chimeric_bub(g, src, in,
	    //                          normal_len, is_clean=0, te=NULL).
	    //
	    // Scans for reads that appear as the short arm of a bubble: the read v
	    // has exactly one outgoing arc (arc0) and one incoming arc (arc1), both
	    // neighbours are branch points, an alternative path of length >normalLen
	    // bypasses v and reconverges at the same merge vertex, and the read's
	    // overlap pattern is consistent with a chimeric join.  Such reads are
	    // deleted.
	    // -----------------------------------------------------------------------
	    auto cutChimericBubblesOnce = [&](uint32_t normalLenArg) -> uint64_t {
	        if (normalLenArg == 0) return 0;
	        uint64_t deletedCount = 0;

	        for (uint32_t v = 0; v < vertexCount; ++v) {
	            if (isVertexDeleted(v)) continue;

	            // v must have exactly one outgoing arc (arc0) and one incoming arc (arc1).
	            // In asg, incoming arcs of v are outgoing arcs of v^1.
	            uint32_t arc0 = 0, arc1 = 0;
	            if (getArcs(v,      &arc0, 1) != 1) continue;
	            if (getArcs(v ^ 1U, &arc1, 1) != 1) continue;

	            // Both neighbours must be branch points: each must have in-degree 2
	            // (= 2 outgoing arcs of their reverse vertex).
	            if (getArcs(g.arcs[arc0].to ^ 1U, nullptr, 0) != 2) continue;
	            if (getArcs(g.arcs[arc1].to ^ 1U, nullptr, 0) != 2) continue;

	            // The read must pass the overlap-based chimera test.
	            const ReadId readId = ReadId(v >> 1U);
	            if (!ifSupChimeric(readId)) continue;

	            // Among the two outgoing arcs of (arc0.to)^1, find the one that does
	            // NOT loop back to v^1 — that is the entry point of the alternative path.
	            const uint32_t branchVertex = g.arcs[arc0].to ^ 1U;
	            const uint32_t backToV1     = v ^ 1U;
	            uint32_t altPathStart = numeric_limits<uint32_t>::max();
	            for (const uint32_t arcId : g.outgoing[branchVertex]) {
	                if (g.arcs[arcId].del)         continue;
	                if (g.arcs[arcId].to == backToV1) continue;  // skip the arc back to v^1
	                altPathStart = g.arcs[arcId].to;
	                break;
	            }
	            if (altPathStart == numeric_limits<uint32_t>::max()) continue;

	            // Walk the alternative path until it hits a merge point (TWO_INPUT).
	            uint32_t altPathEnd = 0, altPathLen = 0;
	            if (followLimitPath(altPathStart, &altPathEnd, &altPathLen,
	                                numeric_limits<uint32_t>::max()) != TWO_INPUT) continue;

	            // The alternative path must be longer than normalLen reads —
	            // a short alternative might just be another artefact.
	            if (altPathLen <= normalLenArg) continue;

	            // The single outgoing arc from altPathEnd must lead to the same merge
	            // vertex that arc1 points to, confirming the bubble closes correctly.
	            uint32_t mergeArc = 0;
	            if (getArcs(altPathEnd, &mergeArc, 1) != 1) continue;
	            if (g.arcs[mergeArc].to != g.arcs[arc1].to) continue;

	            // All conditions met: v is a chimeric read inside a bubble.  Delete it.
	            deleteRead(readId);
	            ++deletedCount;
	        }
	        return deletedCount;
	    };

	    // -----------------------------------------------------------------------
	    // Convergence loop — hifiasm: asg_iterative_semi_circ.
	    // Alternate semi-circular cutting and chimeric-bubble deletion until
	    // neither pass removes anything in a full sweep of the graph.
	    // -----------------------------------------------------------------------
	    uint64_t totalDeleted = 0;
	    while (true) {
	        const uint64_t semiCircDeleted = cutSemiCircularOnce(limLen);
	        const uint64_t chimericDeleted = cutChimericBubblesOnce(normalLen);
	        const uint64_t roundTotal      = semiCircDeleted + chimericDeleted;
	        if (roundTotal == 0) break;
	        totalDeleted += roundTotal;
	    }

	    // Propagate any remaining asymmetric arc deletions and compact the graph
	    // (hifiasm: asg_cleanup called once after the loop).
	    if (totalDeleted > 0) {
	        symmetrizeArcDeletion(g);
	        cleanupStringGraphLikeHifiasm(*this, vertexCount);
	    }
	    return totalDeleted;
	}



// -----------------------------------------------------------------------
// Del-aware BFS bubble pop — like stringGraphBubblePop1Label but does NOT
// require outgoing[] to contain only non-deleted arcs.  Instead it checks
// arc.del inside the loop and uses countOutgoingNonDeleted() wherever the
// original used g.outgoing[...].size().
//
// Used by cleanStringGraphBubbleLinks, which temporarily marks arcs deleted
// without rebuilding adjacency between iterations.
//
// Returns true when a proper bubble was found.  On success ready[0] is the
// merge vertex (same as stringGraphBubblePop1Label).
// -----------------------------------------------------------------------
static bool stringGraphBubblePop1LabelDynamic(
    const StringGraph&      g,
    uint32_t                v0,
    uint64_t                maxDist,
    vector<BubbleBfsState>& bfsState,
    vector<uint32_t>&       ready,
    vector<uint32_t>&       visited
)
{
    if (g.readDeleted.isOpen && g.readDeleted[v0 >> 1U]) return false;
    // Need ≥2 non-deleted outgoing arcs from v0.
    if (countOutgoingNonDeleted(g, v0) < 2) return false;

    ready.clear();
    visited.clear();
    bfsState[v0].dist = 0;
    ready.push_back(v0);

    uint32_t nPending = 0;
    uint32_t nTips    = 0;
    uint32_t tipEnd   = uint32_t(-1);
    bool     found    = false;

    do {
        const uint32_t v  = ready.back(); ready.pop_back();
        const uint32_t dv = bfsState[v].dist;

        bool tooFar = false;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (g.arcs[arcId].del) continue;     // skip deleted arcs
            const uint32_t w = g.arcs[arcId].to;
            const uint32_t l = g.arcs[arcId].len;

            if ((w >> 1U) == (v0 >> 1U)) goto bfs_done_dyn;
            if (uint64_t(dv) + uint64_t(l) > maxDist) { tooFar = true; break; }

            {
                BubbleBfsState& tw = bfsState[w];
                if (tw.visited == 0) {
                    visited.push_back(w);
                    tw.parent    = v;
                    tw.visited   = 1;
                    tw.dist      = dv + l;
                    // Count non-deleted incoming arcs of w = non-deleted outgoing of w^1.
                    tw.remaining = countOutgoingNonDeleted(g, w ^ 1U);
                    ++nPending;
                } else {
                    if (dv + l < tw.dist) { tw.dist = dv + l; tw.parent = v; }
                }
                if (--(tw.remaining) == 0) {
                    if (countOutgoingNonDeleted(g, w) > 0) {
                        ready.push_back(w);
                    } else {
                        if (nTips != 0) goto bfs_done_dyn;
                        ++nTips;
                        tipEnd = w;
                    }
                    --nPending;
                }
            }
        }

        if (tooFar) goto bfs_done_dyn;

        if (nTips == 1) {
            if (tipEnd != uint32_t(-1) && nPending == 0 && ready.empty()) {
                ready.push_back(tipEnd);
                break;
            }
            goto bfs_done_dyn;
        }

        if (ready.empty()) goto bfs_done_dyn;

    } while (ready.size() > 1 || nPending > 0);

    found = true;

bfs_done_dyn:
    for (const uint32_t u : visited) {
        bfsState[u].visited = 0;
        bfsState[u].dist    = 0;
    }
    return found;
}



// -----------------------------------------------------------------------
// 1. cleanStringGraphChimericReads  ≡ hifiasm asg_arc_cut_chimeric
//
// Removes reads that are chimeric based on their overlap-coverage pattern.
// Only acts on reads where at least one incident arc has el==0 (hasLargeIndel).
// For clean HiFi data (all el=1) this is a NOP.
// -----------------------------------------------------------------------
uint64_t Assembler::cleanStringGraphChimericReads()
{
    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphChimericReads: StringGraph not open.");
    }

    const uint32_t readCount   = uint32_t(reads->readCount());
    const uint32_t vertexCount = 2 * readCount;
    StringGraph&   g           = stringGraph;
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph chimeric-read cutting (hifiasm asg_arc_cut_chimeric) begins." << endl;

    // Recompute bubble vertex marks (same as hifiasm does before asg_arc_cut_chimeric).
    rebuildStringGraphAdjacency(*this, vertexCount);
    const vector<uint8_t> seqVis = stringGraphMarkBubbleVertices(g, vertexCount);

    // Helper: effective read length for chimera test (matches cleanStringGraphBreakShortCycles).
    auto readLen = [&](ReadId id) -> uint32_t {
        if (!validReadIntervals.empty() && id < validReadIntervals.size()) {
            const auto& vi = validReadIntervals[id];
            if (!vi.isDeleted && vi.end >= vi.start) return vi.end - vi.start;
        }
        return uint32_t(reads->getRead(id).baseCount);
    };

    uint64_t deletedCount = 0;

    for (uint32_t v = 0; v < vertexCount; ++v) {
        const ReadId readId = ReadId(v >> 1U);
        if (g.readDeleted.isOpen && g.readDeleted[readId]) continue;
        if (v < seqVis.size() && seqVis[v] != 0) continue;

        // v must have exactly 1 non-deleted outgoing arc with el==0.
        uint32_t kv = 0;
        bool hasElZero = false;
        uint32_t arc0id = 0;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (g.arcs[arcId].del) continue;
            ++kv;
            if (!g.arcs[arcId].el) { hasElZero = true; arc0id = arcId; }
        }
        if (kv != 1 || !hasElZero) continue;

        // v^1 must have exactly 1 non-deleted outgoing arc.
        uint32_t kv1 = countOutgoingNonDeleted(g, v ^ 1U);
        if (kv1 != 1) continue;

        // Both targets must have in-degree ≥ 2 (= outgoing of target^1 ≥ 2).
        const uint32_t dest0 = g.arcs[arc0id].to;
        uint32_t arc1id = 0;
        for (const uint32_t arcId : g.outgoing[v ^ 1U]) {
            if (!g.arcs[arcId].del) { arc1id = arcId; break; }
        }
        const uint32_t dest1 = g.arcs[arc1id].to;
        if (countOutgoingNonDeleted(g, dest0 ^ 1U) < 2) continue;
        if (countOutgoingNonDeleted(g, dest1 ^ 1U) < 2) continue;

        // Phase-2 el_n check (≡ hifiasm asg_arc_cut_chimeric Phase 2):
        // dest0^1 must have at least one non-deleted el=1 arc not pointing back to v^1.
        // If all arcs from dest0^1 (other than the back-arc) are el=0, this isn't a true chimera signal.
        {
            const uint32_t wV      = dest0 ^ 1U;
            const uint32_t backV1  = v ^ 1U;
            bool elN = false;
            for (const uint32_t arcId : g.outgoing[wV]) {
                if (g.arcs[arcId].del) continue;
                if (g.arcs[arcId].to == backV1) continue;
                if (g.arcs[arcId].el) { elN = true; break; }
            }
            if (!elN) continue;
        }

        // Chimera check: scan overlaps of readId.
        const uint32_t rLen = readLen(readId);
        uint32_t l1 = 0, r0 = rLen, r1 = 0;

        const span<const uint32_t> alns = alignmentTable[OrientedReadId(readId, 0).getValue()];
        for (const uint32_t idx : alns) {
            const AlignmentData& ad = alignmentData[idx];
            if (!ad.info.isInReadGraph || !ad.keptByBothSides()) continue;
            if (ad.hasLargeIndel) continue;   // el==0 → skip (same as hifiasm)

            uint32_t qs, qe;
            if (ad.readIds[0] == readId) { qs = ad.qs; qe = ad.qe; }
            else                          { qs = ad.ts; qe = ad.te; }
            if (qe > rLen) qe = rLen;
            if (qe <= qs) continue;

            if (qs == 0)    { if (qe > l1) l1 = qe; }
            if (qe == rLen) { if (qs < r0) r0 = qs; r1 = rLen; }
        }

        // Chimeric: gap between left-coverage and right-coverage.
        const bool isChimeric = (r1 > 0 && l1 > 0 && l1 <= r0);
        if (!isChimeric) continue;

        // Delete the read and all incident arcs.
        if (g.readDeleted.isOpen && readId < g.readDeleted.size())
            g.readDeleted[readId] = 1;
        for (const uint32_t arcId : g.outgoing[v]) {
            g.arcs[arcId].del = 1; g.arcs[arcId ^ 1U].del = 1;
        }
        for (const uint32_t arcId : g.outgoing[v ^ 1U]) {
            g.arcs[arcId].del = 1; g.arcs[arcId ^ 1U].del = 1;
        }
        ++deletedCount;
    }

    if (deletedCount > 0) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(*this, vertexCount);
    }

    cout << timestamp << "  deleted " << deletedCount << " chimeric reads (NOP for clean HiFi el=1 data)" << endl;
    cout << timestamp << "String graph chimeric-read cutting complete in "
         << seconds(steady_clock::now() - t0) << " s" << endl;
    return deletedCount;
}



// -----------------------------------------------------------------------
// 2. cleanStringGraphInexactOverlaps  ≡ hifiasm asg_arc_cut_inexact
//
// Removes el==0 (hasLargeIndel) arcs that are weaker than the best overlap
// at their endpoints.  For clean HiFi data (all el=1) this is a NOP.
// -----------------------------------------------------------------------
uint64_t Assembler::cleanStringGraphInexactOverlaps(
    uint32_t maxShortTipReads,
    uint32_t minDiff)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphInexactOverlaps: StringGraph not open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    const int      maxExt      = int(maxShortTipReads);
    StringGraph&   g           = stringGraph;
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph inexact-overlap cutting (hifiasm asg_arc_cut_inexact) begins." << endl;

    // asg_arc_cut_inexact uses seq_vis[v]==0 filter (same as cut_length / cut_bub_links).
    rebuildStringGraphAdjacency(*this, vertexCount);
    const vector<uint8_t> seqVis = stringGraphMarkBubbleVertices(g, vertexCount);

    // Collect el==0 candidate arcs at non-bubble, multi-degree vertices, sorted by overlapLen ascending.
    vector<uint64_t> candidates;
    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) continue;
        if (v < seqVis.size() && seqVis[v] != 0) continue;   // skip bubble vertices
        const uint32_t kv = countOutgoingNonDeleted(g, v);
        if (kv < 2) continue;
        for (const uint32_t arcId : g.outgoing[v]) {
            const auto& a = g.arcs[arcId];
            if (a.del || a.el) continue;   // only el==0 arcs
            candidates.push_back((uint64_t(a.overlapLen) << 32) | uint64_t(arcId));
        }
    }
    std::sort(candidates.begin(), candidates.end());

    uint64_t deletedCount = 0;

    for (const uint64_t cand : candidates) {
        const uint32_t arcId  = uint32_t(cand);
        const uint32_t twinId = arcId ^ 1U;
        if (g.arcs[arcId].del) continue;

        const uint32_t v = g.arcs[arcId].from;
        const uint32_t w = g.arcs[arcId].to ^ 1U;  // dest^1

        if (g.readDeleted.isOpen &&
            (g.readDeleted[v >> 1U] || g.readDeleted[w >> 1U])) continue;

        const uint32_t mmOl = min(g.arcs[arcId].overlapLen, g.arcs[twinId].overlapLen);

        // V-side max overlap.
        uint32_t kv = 0, olMaxV = 0;
        for (const uint32_t aId : g.outgoing[v]) {
            if (g.arcs[aId].del) continue;
            ++kv;
            if (g.arcs[aId].overlapLen > olMaxV) olMaxV = g.arcs[aId].overlapLen;
        }

        // W-side max overlap.
        uint32_t kw = 0, olMaxW = 0;
        for (const uint32_t aId : g.outgoing[w]) {
            if (g.arcs[aId].del) continue;
            ++kw;
            if (g.arcs[aId].overlapLen > olMaxW) olMaxW = g.arcs[aId].overlapLen;
        }

        // V-side: skip if this arc is tied for best on v's outgoing (keep it).
        if (kv >= 2 && mmOl >= olMaxV) continue;
        // V-side: skip if arc is within minDiff of best (always false for minDiff=0).
        if (kv >= 2 && (mmOl + minDiff) > olMaxV) continue;

        // W-side: same independent checks.
        if (kw >= 2 && mmOl >= olMaxW) continue;
        if (kw >= 2 && (mmOl + minDiff) > olMaxW) continue;

        // Bridge guard: don't remove the only connection between two degree-1 vertices.
        if (kv <= 1 && kw <= 1) continue;

        bool toDel = false;
        if (kv > 1 && kw > 1) {
            toDel = true;
        } else if (kw == 1) {
            if (stringGraphTopocutAux(g, w ^ 1U, maxExt) < maxExt) toDel = true;
        } else if (kv == 1) {
            if (stringGraphTopocutAux(g, v ^ 1U, maxExt) < maxExt) toDel = true;
        }

        if (toDel) {
            g.arcs[arcId ].del = 1;
            g.arcs[twinId].del = 1;
            ++deletedCount;
        }
    }

    if (deletedCount > 0) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(*this, vertexCount);
    }

    cout << timestamp << "  deleted " << deletedCount << " inexact-overlap arcs (NOP for clean HiFi el=1 data)" << endl;
    cout << timestamp << "String graph inexact-overlap cutting complete in "
         << seconds(steady_clock::now() - t0) << " s" << endl;
    return deletedCount;
}



// -----------------------------------------------------------------------
// 3. cleanStringGraphBubbleLinks  ≡ hifiasm asg_arc_cut_bub_links
//
// Removes arcs that appear to be "false" multi-incoming links into a bubble.
// Uses del-aware bubble BFS (stringGraphBubblePop1LabelDynamic).
// The `sec_check` path (hifiasm `trans_path_check`) is implemented using
// alignmentData in place of hifiasm's rev/rI overlap arrays.
// -----------------------------------------------------------------------
uint64_t Assembler::cleanStringGraphBubbleLinks(
    double   lenRat,
    double   secLenRat,
    uint32_t maxShortTipReads)
{
    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphBubbleLinks: StringGraph not open.");
    }

    const uint32_t readCount   = uint32_t(reads->readCount());
    const uint32_t vertexCount = 2 * readCount;
    StringGraph&   g           = stringGraph;
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph bubble-link cutting (hifiasm asg_arc_cut_bub_links) begins." << endl;

    rebuildStringGraphAdjacency(*this, vertexCount);
    const vector<uint8_t> seqVis = stringGraphMarkBubbleVertices(g, vertexCount);

    // BFS scratch (reused per call to stringGraphBubblePop1LabelDynamic).
    vector<BubbleBfsState> bfsState(vertexCount);
    vector<uint32_t>       bfsReady;
    vector<uint32_t>       bfsVisited;
    const uint64_t         maxDist = uint64_t(-1) / 4;

    // -----------------------------------------------------------------------
    // Helpers for the sec_check path (≡ hifiasm trans_path_check).
    // Uses alignmentData in place of hifiasm's rev/rI overlap arrays.
    // -----------------------------------------------------------------------

    // Follow the unambiguous chain from oriented vertex `start`.
    // Stops at branches, tips, or when the next vertex has multiple incoming arcs.
    // Returns the read IDs along the path (≡ hifiasm follow_limit_path with lim=∞).
    auto followPath = [&](uint32_t start) -> vector<ReadId> {
        vector<ReadId> pathReads;
        uint32_t v = start;
        while (true) {
            const ReadId rid = ReadId(v >> 1U);
            if (g.readDeleted.isOpen && g.readDeleted[rid]) break;
            pathReads.push_back(rid);
            // Exactly one non-deleted outgoing arc required to continue.
            const uint32_t next = stringGraphCheckUnambiguous1(g, v);
            if (next == std::numeric_limits<uint32_t>::max()) break;
            // Loop detection.
            if ((next >> 1U) == (start >> 1U)) break;
            // Stop if next has multiple incoming arcs (outgoing of next^1 ≥ 2).
            if (countOutgoingNonDeleted(g, next ^ 1U) >= 2) break;
            v = next;
        }
        return pathReads;
    };

    // Persistent mark array for transPathCheck (avoids repeated allocation).
    vector<uint8_t> tpcMark(readCount, 0);

    // Implements hifiasm trans_path_check using alignmentData.
    // Follows paths from oriented vertices `a` and `b`, then checks whether
    // reads on the shorter path have >30% overlap connections to the longer path.
    // Returns true  ≡ trans_path_check == 1
    // Returns false ≡ trans_path_check == 0 or -1
    auto transPathCheck = [&](uint32_t a, uint32_t b) -> bool {
        if (a == b) return false;

        vector<ReadId> pathA = followPath(a);
        vector<ReadId> pathB = followPath(b);

        // Both paths must be longer than maxShortTipReads (≡ l[0] <= minLen → -1).
        if (pathA.size() <= size_t(maxShortTipReads) ||
            pathB.size() <= size_t(maxShortTipReads)) return false;

        // Orient so shorter = pathA, longer = pathB.
        const vector<ReadId>* shorter = &pathA;
        const vector<ReadId>* longer  = &pathB;
        if (pathA.size() > pathB.size()) std::swap(shorter, longer);

        // Mark reads in the longer path.
        for (const ReadId rid : *longer) {
            if (rid < readCount) tpcMark[rid] = 1;
        }

        // Count: minCount = all overlapping reads reachable from shorter;
        //        maxCount = those in the longer path.
        // Uses all alignments for each read (≡ hifiasm rev[qi], no isInReadGraph filter).
        double minCount = 0.0, maxCount = 0.0;
        for (const ReadId qi : *shorter) {
            const span<const uint32_t> alns =
                alignmentTable[OrientedReadId(qi, 0).getValue()];
            for (const uint32_t idx : alns) {
                const AlignmentData& ad = alignmentData[idx];
                const ReadId ti = (ad.readIds[0] == qi) ? ad.readIds[1] : ad.readIds[0];
                if (g.readDeleted.isOpen && ti < g.readDeleted.size() &&
                    g.readDeleted[ti]) continue;
                ++minCount;
                if (ti < readCount && tpcMark[ti]) ++maxCount;
            }
        }

        // Restore marks.
        for (const ReadId rid : *longer) {
            if (rid < readCount) tpcMark[rid] = 0;
        }

        if (minCount == 0.0) return false;
        return (maxCount / minCount) > 0.3;
    };

    // Collect candidate source vertices: seqVis==0, ≥2 non-deleted outgoing arcs.
    // Key is (sum_overlapLen, v) — smaller sum processed first (weakest bubble links).
    vector<pair<uint64_t, uint32_t>> candidates;
    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) continue;
        if (v < seqVis.size() && seqVis[v] != 0) continue;
        const uint32_t kv = countOutgoingNonDeleted(g, v);
        if (kv < 2) continue;
        uint64_t sumOl = 0;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) sumOl += g.arcs[arcId].overlapLen;
        }
        candidates.push_back({sumOl, v});
    }
    std::sort(candidates.begin(), candidates.end());

    uint64_t deletedCount = 0;

    for (const auto& [sumOl, v] : candidates) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) continue;

        // For each non-deleted outgoing arc from v, check the "bubble link" condition.
        // Collect non-del outgoing arcs.
        vector<uint32_t> vArcs;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) vArcs.push_back(arcId);
        }
        if (vArcs.size() < 2) continue;

        bool allPass = true;
        uint32_t secCheck = 0;
        // For each arm arc: v → target_v
        // w = target_v ^ 1 is the reverse side of the target read.
        // Check that kw >= 2 and the arc overlap is short relative to w's best incoming.
        vector<uint32_t> armIncoming;  // arcs of w that are NOT the back-arc to v
        armIncoming.clear();

        for (const uint32_t arcId : vArcs) {
            const uint32_t targetV = g.arcs[arcId].to;
            const uint32_t w = targetV ^ 1U;

            const uint32_t kw = countOutgoingNonDeleted(g, w);
            if (kw < 2) { allPass = false; break; }

            // me = min overlapLen of w's outgoing arcs, excluding the arc back to v
            //      (≡ hifiasm: initialized to (uint32_t)-1, updated with "if(ol < me) me = ol").
            //      Comment in hifiasm: "me is the shortest edge except aw[t]".
            uint32_t me = std::numeric_limits<uint32_t>::max();
            const uint32_t backTo = v ^ 1U;
            for (const uint32_t aid : g.outgoing[w]) {
                if (g.arcs[aid].del) continue;
                if (g.arcs[aid].to == backTo) continue;
                if (g.arcs[aid].overlapLen < me) me = g.arcs[aid].overlapLen;
                armIncoming.push_back(aid);
            }
            if (me == std::numeric_limits<uint32_t>::max()) { allPass = false; break; }  // no arm arcs

            const double ol = double(g.arcs[arcId].overlapLen);
            if (ol > double(me) * lenRat) {
                if (ol > double(me) * secLenRat) { allPass = false; break; }
                ++secCheck;
            }
        }

        if (!allPass) continue;

        if (secCheck > 0) {
            // ---------------------------------------------------------------
            // sec_check path: ≡ hifiasm trans_path_check on vArcs targets.
            //
            // Forward check: all pairs of vArcs[0].to vs vArcs[i].to (i>0)
            // must be transitively connected via alignmentData.
            // If forward check fails, try arm incoming check instead.
            // ---------------------------------------------------------------
            bool forwardOk = true;
            uint32_t firstTarget = std::numeric_limits<uint32_t>::max();
            for (const uint32_t arcId : vArcs) {
                if (g.arcs[arcId].del) continue;
                const uint32_t target = g.arcs[arcId].to;
                if (firstTarget == std::numeric_limits<uint32_t>::max()) {
                    firstTarget = target;
                } else if (!transPathCheck(firstTarget, target)) {
                    forwardOk = false;
                    break;
                }
            }

            if (!forwardOk) {
                // Arm incoming check: sources of arm arcs must have degree 2,
                // and their targets must be transitively connected pairwise.
                if (armIncoming.size() < 2) continue;
                bool armOk = true;
                uint32_t firstArmTarget = std::numeric_limits<uint32_t>::max();
                for (const uint32_t aid : armIncoming) {
                    if (g.arcs[aid].del) continue;
                    // Source (w) must have exactly 2 non-deleted outgoing arcs.
                    if (countOutgoingNonDeleted(g, g.arcs[aid].from) != 2) {
                        armOk = false;
                        break;
                    }
                    const uint32_t armTarget = g.arcs[aid].to;
                    if (firstArmTarget == std::numeric_limits<uint32_t>::max()) {
                        firstArmTarget = armTarget;
                    } else if (!transPathCheck(firstArmTarget, armTarget)) {
                        armOk = false;
                        break;
                    }
                }
                if (!armOk) continue;
            }
        }

        // Test whether removing armIncoming arcs reveals a clean double bubble.
        // Temporarily mark armIncoming arcs (and their twins) as deleted.
        for (const uint32_t aid : armIncoming) {
            g.arcs[aid      ].del = 1;
            g.arcs[aid ^ 1U ].del = 1;
        }

        bool isFalseBub = false;
        uint32_t mergeVertex = std::numeric_limits<uint32_t>::max();

        if (countOutgoingNonDeleted(g, v) >= 2) {
            // First bubble pop: from v.
            if (stringGraphBubblePop1LabelDynamic(g, v, maxDist, bfsState, bfsReady, bfsVisited)) {
                mergeVertex = bfsReady.empty() ? std::numeric_limits<uint32_t>::max() : bfsReady[0];
            }
        }

        if (mergeVertex != std::numeric_limits<uint32_t>::max()) {
            // Temporarily delete v's non-del outgoing arcs.
            vector<uint32_t> vArcsTmp;
            for (const uint32_t arcId : g.outgoing[v]) {
                if (!g.arcs[arcId].del) {
                    vArcsTmp.push_back(arcId);
                    g.arcs[arcId      ].del = 1;
                    g.arcs[arcId ^ 1U ].del = 1;
                }
            }

            // Restore arm incoming first (they were blocking the second pop).
            for (const uint32_t aid : armIncoming) {
                g.arcs[aid      ].del = 0;
                g.arcs[aid ^ 1U ].del = 0;
            }

            // Second bubble pop: from mergeVertex^1.
            if (stringGraphBubblePop1LabelDynamic(g, mergeVertex ^ 1U, maxDist, bfsState, bfsReady, bfsVisited)) {
                isFalseBub = true;
            }

            // Restore vArcsTmp if not a false bubble.
            if (!isFalseBub) {
                for (const uint32_t arcId : vArcsTmp) {
                    g.arcs[arcId      ].del = 0;
                    g.arcs[arcId ^ 1U ].del = 0;
                }
                // Re-delete armIncoming since we didn't commit.
                for (const uint32_t aid : armIncoming) {
                    g.arcs[aid      ].del = 1;
                    g.arcs[aid ^ 1U ].del = 1;
                }
            }
            // On success vArcsTmp stays deleted and armIncoming stays restored.
        }

        if (!isFalseBub) {
            // Restore everything.
            for (const uint32_t aid : armIncoming) {
                g.arcs[aid      ].del = 0;
                g.arcs[aid ^ 1U ].del = 0;
            }
        } else {
            deletedCount++;
        }
    }

    if (deletedCount > 0) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(*this, vertexCount);
    }

    cout << timestamp << "  deleted " << deletedCount << " bubble-link arc groups (asg_arc_cut_bub_links)" << endl;
    cout << timestamp << "String graph bubble-link cutting complete in "
         << seconds(steady_clock::now() - t0) << " s" << endl;
    return deletedCount;
}



// -----------------------------------------------------------------------
// 4. cleanStringGraphComplexBubbleLinks  ≡ hifiasm asg_arc_cut_complex_bub_links
//
// Removes all outgoing arcs of a vertex v if each target has multiple
// incoming arcs and all arc overlaps are short relative to the best
// alternative incoming at that target.  After deletion, verifies with
// stringGraphMarkBubbleVertices that deleted arcs are not in bubbles.
// -----------------------------------------------------------------------
uint64_t Assembler::cleanStringGraphComplexBubbleLinks(
    double lenRat)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphComplexBubbleLinks: StringGraph not open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    StringGraph&   g           = stringGraph;
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph complex-bubble-link cutting (hifiasm asg_arc_cut_complex_bub_links) begins." << endl;

    rebuildStringGraphAdjacency(*this, vertexCount);
    const vector<uint8_t> seqVis = stringGraphMarkBubbleVertices(g, vertexCount);

    // Collect candidates.
    vector<pair<uint64_t, uint32_t>> candidates;
    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) continue;
        if (v < seqVis.size() && seqVis[v] != 0) continue;
        if (countOutgoingNonDeleted(g, v) < 2) continue;
        uint64_t sumOl = 0;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) sumOl += g.arcs[arcId].overlapLen;
        }
        candidates.push_back({sumOl, v});
    }
    std::sort(candidates.begin(), candidates.end());

    vector<uint64_t> deletedArcIds;

    for (const auto& [sumOl, v] : candidates) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) continue;

        vector<uint32_t> vArcs;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) vArcs.push_back(arcId);
        }
        if (vArcs.size() < 2) continue;

        bool allPass = true;
        for (const uint32_t arcId : vArcs) {
            const uint32_t targetV = g.arcs[arcId].to;
            const uint32_t w = targetV ^ 1U;

            const uint32_t kw = countOutgoingNonDeleted(g, w);
            if (kw < 2) { allPass = false; break; }

            const uint32_t backTo = v ^ 1U;
            uint32_t me = std::numeric_limits<uint32_t>::max();  // min arm overlap (≡ hifiasm "shortest edge")
            for (const uint32_t aid : g.outgoing[w]) {
                if (g.arcs[aid].del) continue;
                if (g.arcs[aid].to == backTo) continue;
                if (g.arcs[aid].overlapLen < me) me = g.arcs[aid].overlapLen;
            }
            if (me == std::numeric_limits<uint32_t>::max()) { allPass = false; break; }  // no arm arcs

            if (double(g.arcs[arcId].overlapLen) > double(me) * lenRat)
                { allPass = false; break; }
        }

        if (!allPass) continue;

        // Delete all of v's non-deleted outgoing arcs.
        for (const uint32_t arcId : vArcs) {
            g.arcs[arcId      ].del = 1;
            g.arcs[arcId ^ 1U ].del = 1;
            deletedArcIds.push_back(uint64_t(arcId));
        }
    }

    uint64_t permanentlyDeleted = 0;

    if (!deletedArcIds.empty()) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(*this, vertexCount);

        // Re-compute bubble vertices after deletions (≡ asg_arc_identify_simple_bubbles_multi).
        // Hifiasm logic: if an arc's endpoint IS a bubble node → keep it deleted (true bubble link).
        //               if an arc's endpoint is NOT a bubble node → restore it (over-deleted).
        const vector<uint8_t> seqVisNew = stringGraphMarkBubbleVertices(g, vertexCount);
        bool anyRestored = false;
        for (const uint64_t aid64 : deletedArcIds) {
            const uint32_t arcId  = uint32_t(aid64);
            const uint32_t twinId = arcId ^ 1U;
            if (arcId >= g.arcs.size() || twinId >= g.arcs.size()) continue;
            const uint32_t fv = g.arcs[arcId ].from;
            const uint32_t tv = g.arcs[arcId ].to;
            const bool touchesBubble =
                (fv < seqVisNew.size() && seqVisNew[fv] != 0) ||
                (tv < seqVisNew.size() && seqVisNew[tv] != 0) ||
                ((fv ^ 1U) < seqVisNew.size() && seqVisNew[fv ^ 1U] != 0) ||
                ((tv ^ 1U) < seqVisNew.size() && seqVisNew[tv ^ 1U] != 0);
            if (touchesBubble) {
                // Endpoint is a bubble node → keep deleted (correctly identified bubble link).
                ++permanentlyDeleted;
            } else {
                // Endpoint is NOT a bubble node → restore (over-deleted, not a true bubble link).
                g.arcs[arcId ].del = 0;
                g.arcs[twinId].del = 0;
                anyRestored = true;
            }
        }
        if (anyRestored) {
            symmetrizeArcDeletion(g);
            rebuildStringGraphAdjacency(*this, vertexCount);
        }
    }
    cout << timestamp << "  deleted " << permanentlyDeleted << " complex-bubble-link arcs (asg_arc_cut_complex_bub_links)" << endl;
    cout << timestamp << "String graph complex-bubble-link cutting complete in "
         << seconds(steady_clock::now() - t0) << " s" << endl;
    return permanentlyDeleted;
}



// -----------------------------------------------------------------------
// 5. cleanStringGraphLargeIndelArcs  ≡ hifiasm asg_cut_large_indel
//
// Removes el==0 (hasLargeIndel) arcs that are weaker than the best
// non-indel overlap at their endpoints.  For clean HiFi data (all el=1)
// this is a NOP.
// -----------------------------------------------------------------------
uint64_t Assembler::cleanStringGraphLargeIndelArcs(
    uint32_t maxShortTipReads,
    uint32_t minDiff)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphLargeIndelArcs: StringGraph not open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    const int      maxExt      = int(maxShortTipReads);
    StringGraph&   g           = stringGraph;
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph large-indel arc cutting (hifiasm asg_cut_large_indel) begins." << endl;

    // seqVis needed to skip bubble nodes (same as asg_cut_large_indel).
    rebuildStringGraphAdjacency(*this, vertexCount);
    const vector<uint8_t> seqVis = stringGraphMarkBubbleVertices(g, vertexCount);

    vector<uint64_t> candidates;
    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) continue;
        if (v < seqVis.size() && seqVis[v] != 0) continue;
        const uint32_t kv = countOutgoingNonDeleted(g, v);
        if (kv < 2) continue;
        for (const uint32_t arcId : g.outgoing[v]) {
            const auto& a = g.arcs[arcId];
            if (a.del || a.el) continue;   // only el==0 (large indel) arcs
            candidates.push_back((uint64_t(a.overlapLen) << 32) | uint64_t(arcId));
        }
    }
    std::sort(candidates.begin(), candidates.end());

    uint64_t deletedCount = 0;

    for (const uint64_t cand : candidates) {
        const uint32_t arcId  = uint32_t(cand);
        const uint32_t twinId = arcId ^ 1U;
        if (g.arcs[arcId].del) continue;

        const uint32_t v = g.arcs[arcId].from;
        const uint32_t w = g.arcs[arcId].to ^ 1U;

        if (g.readDeleted.isOpen &&
            (g.readDeleted[v >> 1U] || g.readDeleted[w >> 1U])) continue;

        uint32_t kv = 0, olMaxV = 0;
        for (const uint32_t aId : g.outgoing[v]) {
            if (g.arcs[aId].del) continue;
            ++kv;
            if (g.arcs[aId].overlapLen > olMaxV) olMaxV = g.arcs[aId].overlapLen;
        }

        uint32_t kw = 0, olMaxW = 0;
        for (const uint32_t aId : g.outgoing[w]) {
            if (g.arcs[aId].del) continue;
            ++kw;
            if (g.arcs[aId].overlapLen > olMaxW) olMaxW = g.arcs[aId].overlapLen;
        }

        const uint32_t veOl = g.arcs[arcId].overlapLen;
        const uint32_t weOl = g.arcs[twinId].overlapLen;

        // V-side: skip if this arc's overlap exceeds the v-side best (can't improve).
        // For minDiff==0 this condition is always false (veOl <= olMaxV by construction).
        if (kv >= 2 && (veOl + minDiff) > olMaxV) continue;

        // W-side: skip similarly.
        if (kw >= 2 && (weOl + minDiff) > olMaxW) continue;

        // Bridge guard: both sides degree-1 means removing this arc disconnects the graph.
        if (kv <= 1 && kw <= 1) continue;

        bool toDel = false;
        if (kv > 1 && kw > 1) {
            toDel = true;
        } else if (kw == 1) {
            if (stringGraphTopocutAux(g, w ^ 1U, maxExt) < maxExt) toDel = true;
        } else if (kv == 1) {
            if (stringGraphTopocutAux(g, v ^ 1U, maxExt) < maxExt) toDel = true;
        }

        if (toDel) {
            g.arcs[arcId ].del = 1;
            g.arcs[twinId].del = 1;
            ++deletedCount;
        }
    }

    if (deletedCount > 0) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(*this, vertexCount);
    }

    cout << timestamp << "  deleted " << deletedCount << " large-indel arcs (NOP for clean HiFi el=1 data)" << endl;
    cout << timestamp << "String graph large-indel arc cutting complete in "
         << seconds(steady_clock::now() - t0) << " s" << endl;
    return deletedCount;
}



// -----------------------------------------------------------------------
// 6. cutStringGraphSemiCircular  ≡ hifiasm asg_cut_semi_circ(sg, limLen, 1)
//
// Standalone post-loop semi-circular arc cut (single pass, is_clean=1).
// Removes the back-arc completing a short semi-circle at each qualifying
// vertex v (v has >1 incoming, =1 outgoing, and the forward path is short).
// -----------------------------------------------------------------------
uint64_t Assembler::cutStringGraphSemiCircular(uint32_t limLen)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cutStringGraphSemiCircular: StringGraph not open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    StringGraph&   g           = stringGraph;
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph semi-circular arc cutting (hifiasm asg_cut_semi_circ) begins." << endl;

    // Rebuild adjacency first (is_clean=1 in hifiasm).
    rebuildStringGraphAdjacency(*this, vertexCount);

    uint64_t deletedCount = 0;

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) continue;

        // v must have >1 incoming arcs (= outgoing of v^1 > 1).
        if (g.outgoing[v ^ 1U].size() <= 1) continue;
        if (countOutgoingNonDeleted(g, v ^ 1U) <= 1) continue;

        // v must have exactly one outgoing arc.
        if (g.outgoing[v].empty()) continue;
        if (countOutgoingNonDeleted(g, v) != 1) continue;

        // Walk the unique path forward from v, up to limLen steps.
        // Mirrors hifiasm follow_limit_path (END_TIPS / TWO_INPUT / TWO_OUTPUT / LOOP).
        uint32_t cursor  = v;
        uint32_t steps   = 0;
        bool     isLoop  = false;
        bool     isEnd   = false;   // END_TIPS: dead-end, no deletion (hifiasm skips)

        while (steps < limLen) {
            // Get the single non-deleted outgoing arc.
            uint32_t nextArcId = std::numeric_limits<uint32_t>::max();
            for (const uint32_t arcId : g.outgoing[cursor]) {
                if (!g.arcs[arcId].del) { nextArcId = arcId; break; }
            }
            if (nextArcId == std::numeric_limits<uint32_t>::max()) { isEnd = true; break; }  // END_TIPS

            const uint32_t next = g.arcs[nextArcId].to;
            ++steps;

            // Loop detection (≡ hifiasm LOOP).
            if (next == v) { isLoop = true; break; }

            // TWO_INPUT: next has ≥2 incoming arcs (= outgoing of next^1 ≥ 2).
            // hifiasm sets e = cursor (the vertex BEFORE next) and returns TWO_INPUT.
            // Do NOT advance cursor — it stays at the vertex before the join point.
            if (countOutgoingNonDeleted(g, next ^ 1U) >= 2) break;

            cursor = next;

            // TWO_OUTPUT / MUL_OUTPUT: cursor has ≠1 outgoing arc.
            // cursor (= next) is the branch vertex; hifiasm returns TWO_OUTPUT with e = cursor.
            if (countOutgoingNonDeleted(g, cursor) != 1) break;
        }

        if (steps == 0 || isLoop || isEnd) continue;

        // cursor is the end vertex `e`.
        // Delete the arc v^1 → cursor^1 (if it exists) and its twin.
        const uint32_t backTarget = cursor ^ 1U;
        for (const uint32_t arcId : g.outgoing[v ^ 1U]) {
            if (g.arcs[arcId].del) continue;
            if (g.arcs[arcId].to == backTarget) {
                g.arcs[arcId      ].del = 1;
                g.arcs[arcId ^ 1U ].del = 1;
                ++deletedCount;
                break;
            }
        }
    }

    if (deletedCount > 0) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(*this, vertexCount);
    }

    cout << timestamp << "  deleted " << deletedCount << " semi-circular arcs (asg_cut_semi_circ)" << endl;
    cout << timestamp << "String graph semi-circular cutting complete in "
         << seconds(steady_clock::now() - t0) << " s" << endl;
    return deletedCount;
}
