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
    for (uint32_t v = 0; v < vertexCount; ++v) {
        uint32_t* b = assembler.stringGraph.outgoing.begin(v);
        uint32_t* e = assembler.stringGraph.outgoing.end(v);
        std::sort(b, e, [&](uint32_t aId, uint32_t bId) {
            const auto& a = assembler.stringGraph.arcs[aId];
            const auto& bArc = assembler.stringGraph.arcs[bId];
            if (a.len != bArc.len) return a.len < bArc.len;
            return a.to < bArc.to;
        });
    }
}



static uint64_t stringGraphTransitiveReduce(Assembler& assembler, uint32_t fuzz, uint32_t vertexCount)
{
    StringGraph& g = assembler.stringGraph;

    // Ensure adjacency exists and is sorted by `len`.
    rebuildStringGraphAdjacency(assembler, vertexCount);

    vector<uint8_t> mark(vertexCount, 0);
    uint64_t reduced = 0;

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) {
            for (const uint32_t arcId : g.outgoing[v]) {
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
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(assembler, vertexCount);
    }
    return reduced;
}



static uint64_t stringGraphCutTips(Assembler& assembler, uint32_t maxShortTipReads, uint32_t vertexCount)
{
    StringGraph& g = assembler.stringGraph;

    uint64_t cutCount = 0;
    vector<uint32_t> path;

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) {
            continue;
        }
        if (stringGraphIsUtgEnd(g, v, nullptr) != ASG_ET_TIP) {
            continue;
        }
        const int ret = stringGraphExtend(g, v, int(maxShortTipReads), path);
        if (ret == ASG_ET_MERGEABLE) {
            continue;
        }

        for (const uint32_t vv : path) {
            const ReadId readId = ReadId(vv >> 1U);
            if (g.readDeleted.isOpen && readId < g.readDeleted.size()) {
                g.readDeleted[readId] = 1;
            }
            const uint32_t v0 = uint32_t(readId) << 1U;
            const uint32_t v1 = v0 ^ 1U;
            for (const uint32_t arcId : g.outgoing[v0]) {
                g.arcs[arcId].del = 1;
            }
            for (const uint32_t arcId : g.outgoing[v1]) {
                g.arcs[arcId].del = 1;
            }
        }

        ++cutCount;
    }

    if (cutCount) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(assembler, vertexCount);
    }

    return cutCount;
}

static uint64_t stringGraphRemoveSingleNodeBubbles(Assembler& assembler, uint32_t maxShortTipReads, uint32_t vertexCount)
{
    // This is a simplified, topology-only version of hifiasm's `asg_arc_del_single_node_directly`.
    // It removes the weaker branch of a 2-branch bubble that reconverges in one step:
    //   v -> a -> t
    //   v -> b -> t
    // with a and b having in-degree=1 and out-degree=1, and v having out-degree=2.
    StringGraph& g = assembler.stringGraph;
    if (!g.incoming.isOpen()) {
        return 0;
    }

    uint64_t reduced = 0;

    auto getNonDeletedOutgoingArcs = [&](uint32_t v, array<uint32_t, 2>& outArcIds, uint32_t& count) {
        count = 0;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (g.arcs[arcId].del) continue;
            if (count < 2) outArcIds[count] = arcId;
            ++count;
        }
    };

    auto uniqueOutgoingTo = [&](uint32_t v, uint32_t& to) -> bool {
        uint32_t arcId = 0;
        if (firstOutgoingNonDeleted(g, v, arcId) != 1) return false;
        if (countOutgoingNonDeleted(g, v) != 1) return false;
        to = g.arcs[arcId].to;
        return true;
    };

    auto incomingNonDeletedCount = [&](uint32_t v) -> uint32_t {
        uint32_t n = 0;
        for (const uint32_t arcId : g.incoming[v]) {
            if (!g.arcs[arcId].del) ++n;
        }
        return n;
    };

    auto hasIncomingFrom = [&](uint32_t v, uint32_t from) -> bool {
        for (const uint32_t arcId : g.incoming[v]) {
            const auto& a = g.arcs[arcId];
            if (a.del) continue;
            if (a.from == from) return true;
        }
        return false;
    };

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) {
            continue;
        }

        array<uint32_t, 2> outArcIds{};
        uint32_t outCount = 0;
        getNonDeletedOutgoingArcs(v, outArcIds, outCount);
        if (outCount != 2) {
            continue;
        }

        const uint32_t arcVa = outArcIds[0];
        const uint32_t arcVb = outArcIds[1];
        const uint32_t a = g.arcs[arcVa].to;
        const uint32_t b = g.arcs[arcVb].to;
        if (a == b) {
            continue;
        }

        if (incomingNonDeletedCount(a) != 1 || incomingNonDeletedCount(b) != 1) {
            continue;
        }
        if (!hasIncomingFrom(a, v) || !hasIncomingFrom(b, v)) {
            continue;
        }

        uint32_t tA = 0;
        uint32_t tB = 0;
        if (!uniqueOutgoingTo(a, tA) || !uniqueOutgoingTo(b, tB) || tA != tB) {
            continue;
        }
        const uint32_t t = tA;

        // Heuristic scoring: prefer the branch with more total overlap.
        uint32_t arcAtoT = 0;
        uint32_t arcBtoT = 0;
        {
            uint32_t tmp = 0;
            firstOutgoingNonDeleted(g, a, tmp);
            arcAtoT = tmp;
            firstOutgoingNonDeleted(g, b, tmp);
            arcBtoT = tmp;
        }
        const uint32_t scoreA = g.arcs[arcVa].overlapLen + g.arcs[arcAtoT].overlapLen;
        const uint32_t scoreB = g.arcs[arcVb].overlapLen + g.arcs[arcBtoT].overlapLen;

        // Only remove if at least one branch is a short unitig (<= maxShortTipReads reads),
        // mirroring hifiasm's intent for this pass (very small bubbles).
        // Our bubble has exactly one intermediate node per branch => 2 reads per branch path (a or b),
        // so treat it as length 1 edge between v and t with one interior node.
        if (maxShortTipReads < 2) {
            continue;
        }

        const bool removeB = (scoreA >= scoreB);
        const uint32_t delArcVx = removeB ? arcVb : arcVa;
        const uint32_t delArcXtoT = removeB ? arcBtoT : arcAtoT;

        g.arcs[delArcVx].del = 1;
        g.arcs[delArcVx ^ 1U].del = 1;
        g.arcs[delArcXtoT].del = 1;
        g.arcs[delArcXtoT ^ 1U].del = 1;
        reduced += 4;
    }

    if (reduced) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(assembler, vertexCount);
    }
    return reduced;
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

    const uint64_t reduced = stringGraphTransitiveReduce(*this, gapFuzz, vertexCount);
    cout << timestamp << "  transitively reduced " << reduced << " arcs (gapFuzz=" << gapFuzz << ")" << endl;

    if (maxShortTipReads > 0) {
        const uint64_t cut = stringGraphCutTips(*this, maxShortTipReads, vertexCount);
        cout << timestamp << "  cut " << cut << " tips (maxShortTipReads=" << maxShortTipReads << ")" << endl;
    }

    cout << timestamp << "String graph clean complete in " << seconds(steady_clock::now() - t0) << " s" << endl;
}



void Assembler::cleanStringGraphPreCleanHifiasm(uint32_t maxShortTipReads)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphPreCleanHifiasm: StringGraph is not open.");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    const auto t0 = steady_clock::now();
    cout << timestamp << "String graph pre-clean (small bubbles + tip cut) begins." << endl;

    uint64_t totalBubbleArcRemoved = 0;
    uint64_t totalTipsCut = 0;
    uint64_t totalCyclesBroken = 0;
    for (uint32_t iter = 0; iter < 10; ++iter) {
        // Break simple short cycles (hifiasm pre_clean removes simple circles).
        const uint64_t cyclesBroken = cleanStringGraphBreakShortCycles(/*maxCycleReads*/100);
        totalCyclesBroken += cyclesBroken;

        const uint64_t bubbleRemoved = stringGraphRemoveSingleNodeBubbles(*this, maxShortTipReads, vertexCount);
        totalBubbleArcRemoved += bubbleRemoved;

        const uint64_t tipsCut = stringGraphCutTips(*this, maxShortTipReads, vertexCount);
        totalTipsCut += tipsCut;

        if (cyclesBroken == 0 && bubbleRemoved == 0 && tipsCut == 0) {
            break;
        }
    }

    cout << timestamp << "  broke " << totalCyclesBroken << " cycle arcs; removed " << totalBubbleArcRemoved << " bubble arcs; cut " << totalTipsCut << " tips" << endl;
    cout << timestamp << "String graph pre-clean complete in " << seconds(steady_clock::now() - t0) << " s" << endl;
}



void Assembler::cleanStringGraphDropShortOverlaps(double dropRatio, uint32_t minOverlapLen)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen()) {
        throw runtime_error("cleanStringGraphDropShortOverlaps: StringGraph is not open.");
    }

    if (dropRatio <= 0. || dropRatio >= 1.) {
        throw runtime_error("cleanStringGraphDropShortOverlaps: dropRatio must be in (0,1).");
    }

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    StringGraph& g = stringGraph;

    // Ensure adjacency lists are consistent and sorted by `len` (not required here, but keeps invariant).
    rebuildStringGraphAdjacency(*this, vertexCount);

    vector<uint32_t> arcIds;
    arcIds.reserve(64);

    uint64_t removed = 0;
    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.readDeleted.isOpen && g.readDeleted[v >> 1U]) {
            continue;
        }
        const span<const uint32_t> out = g.outgoing[v];
        if (out.size() < 2) continue;

        arcIds.clear();
        uint32_t maxOl = 0;
        for (const uint32_t arcId : out) {
            const auto& a = g.arcs[arcId];
            if (a.del) continue;
            arcIds.push_back(arcId);
            maxOl = max(maxOl, a.overlapLen);
        }
        if (arcIds.size() < 2) continue;

        // Hifiasm's asg_arc_del_short uses overlap length ratio threshold.
        const uint32_t thres = max<uint32_t>(uint32_t(double(maxOl) * dropRatio + 0.499), minOverlapLen);

        // Keep the best overlap and drop the rest below threshold.
        std::ranges::sort(arcIds, [&](uint32_t aId, uint32_t bId) {
            const auto& a = g.arcs[aId];
            const auto& b = g.arcs[bId];
            if (a.overlapLen != b.overlapLen) return a.overlapLen > b.overlapLen;
            if (a.len != b.len) return a.len < b.len;
            return a.to < b.to;
        });

        for (size_t i = 1; i < arcIds.size(); ++i) {
            const uint32_t arcId = arcIds[i];
            if (g.arcs[arcId].overlapLen < thres) {
                g.arcs[arcId].del = 1;
                g.arcs[arcId ^ 1U].del = 1;
                removed += 2;
            }
        }
    }

    if (removed) {
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
        cleanStringGraphDropShortOverlaps(drop, /*minOverlapLen*/0);
        cleanStringGraphPreCleanHifiasm(maxShortTipReads);
    }

    // Final pass similar in spirit to `asg_arc_del_too_short_overlaps(sg, 2000, min_drop_rate, ...)`.
    if (finalMinOverlapLen > 0) {
        cleanStringGraphDropShortOverlaps(minDropRate, finalMinOverlapLen);
        cleanStringGraphPreCleanHifiasm(maxShortTipReads);
    }

    cout << timestamp << "String graph overlap-drop rounds complete in " << seconds(steady_clock::now() - t0) << " s" << endl;
}



uint64_t Assembler::cleanStringGraphBreakShortCycles(uint32_t maxCycleReads)
{
    reads->checkReadsAreOpen();
    if (!stringGraph.arcs.isOpen || !stringGraph.outgoing.isOpen() || !stringGraph.incoming.isOpen()) {
        throw runtime_error("cleanStringGraphBreakShortCycles: StringGraph arcs/outgoing/incoming must be open.");
    }
    if (maxCycleReads == 0) return 0;

    const uint32_t vertexCount = uint32_t(2 * reads->readCount());
    StringGraph& g = stringGraph;

    // Ensure adjacency is consistent (also keeps outgoing sorted by len).
    rebuildStringGraphAdjacency(*this, vertexCount);

    vector<uint8_t> visited(vertexCount, 0);
    vector<uint32_t> path;
    path.reserve(maxCycleReads + 2);

    auto outDeg = [&](uint32_t v) -> uint32_t {
        uint32_t n = 0;
        for (const uint32_t arcId : g.outgoing[v]) if (!g.arcs[arcId].del) ++n;
        return n;
    };
    auto inDeg = [&](uint32_t v) -> uint32_t {
        uint32_t n = 0;
        for (const uint32_t arcId : g.incoming[v]) if (!g.arcs[arcId].del) ++n;
        return n;
    };
    auto uniqueOutArc = [&](uint32_t v, uint32_t& arcIdOut) -> bool {
        if (outDeg(v) != 1) return false;
        return firstOutgoingNonDeleted(g, v, arcIdOut) == 1;
    };

    uint64_t broken = 0;
    for (uint32_t start = 0; start < vertexCount; ++start) {
        if (visited[start]) continue;
        if (g.readDeleted.isOpen && g.readDeleted[start >> 1U]) {
            visited[start] = 1;
            continue;
        }
        // Only consider vertices that are part of a simple 1-in/1-out unitig chain.
        if (outDeg(start) != 1 || inDeg(start) != 1) {
            visited[start] = 1;
            continue;
        }

        path.clear();
        uint32_t v = start;
        uint32_t steps = 0;
        bool cycle = false;
        while (true) {
            if (steps++ > maxCycleReads) {
                break;
            }
            if (visited[v] && v != start) {
                break;
            }
            path.push_back(v);

            uint32_t arcId = 0;
            if (!uniqueOutArc(v, arcId)) {
                break;
            }
            const uint32_t next = g.arcs[arcId].to;

            if (next == start) {
                // Close the cycle if the next vertex also satisfies 1-in/1-out.
                cycle = true;
                path.push_back(next);
                break;
            }

            // Continue only along 1-in/1-out nodes to ensure the cycle is "simple".
            if (outDeg(next) != 1 || inDeg(next) != 1) {
                break;
            }
            v = next;
        }

        for (const uint32_t u : path) {
            visited[u] = 1;
        }

        if (!cycle) {
            continue;
        }
        // Path includes start again as last element; cycle length in vertices is (path.size()-1).
        if (path.size() <= 2) {
            continue;
        }
        if (path.size() - 1 > maxCycleReads) {
            continue;
        }

        // Break the cycle by deleting the last edge (prev -> start), like hifiasm removing a back edge.
        const uint32_t prev = path[path.size() - 2];
        uint32_t arcPrev = 0;
        if (!uniqueOutArc(prev, arcPrev)) {
            continue;
        }
        if (g.arcs[arcPrev].to != start) {
            continue;
        }
        g.arcs[arcPrev].del = 1;
        g.arcs[arcPrev ^ 1U].del = 1;
        broken += 2;
    }

    if (broken) {
        symmetrizeArcDeletion(g);
        rebuildStringGraphAdjacency(*this, vertexCount);
    }
    return broken;
}
