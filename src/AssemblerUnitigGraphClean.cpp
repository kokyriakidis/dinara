// Dinara.
#include "Assembler.hpp"
#include "timestamp.hpp"
#include "chrono.hpp"

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

    inline void symmetrizeArcDeletion(UnitigGraph& g)
    {
        DINARA_ASSERT((g.arcs.size() & 1ULL) == 0);
        for (uint64_t arcId = 0; arcId < g.arcs.size(); arcId += 2) {
            const uint8_t del = uint8_t(g.arcs[arcId].del | g.arcs[arcId ^ 1ULL].del);
            g.arcs[arcId].del = del;
            g.arcs[arcId ^ 1ULL].del = del;
        }
    }

    inline uint32_t countOutgoingNonDeleted(const UnitigGraph& g, uint32_t v)
    {
        uint32_t n = 0;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) ++n;
        }
        return n;
    }

    inline uint32_t firstOutgoingNonDeleted(const UnitigGraph& g, uint32_t v, uint32_t& arcIdOut)
    {
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) {
                arcIdOut = arcId;
                return 1;
            }
        }
        return 0;
    }

    inline int unitigGraphIsUtgEnd(const UnitigGraph& g, uint32_t v, uint32_t* nextVertex)
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

    inline int unitigGraphExtend(const UnitigGraph& g, uint32_t v, int maxExt, vector<uint32_t>& path)
    {
        path.clear();
        path.push_back(v);
        int ret = ASG_ET_MERGEABLE;
        do {
            uint32_t nextVertex = 0;
            ret = unitigGraphIsUtgEnd(g, v ^ 1U, &nextVertex);
            if (ret != ASG_ET_MERGEABLE) {
                break;
            }
            path.push_back(nextVertex);
            v = nextVertex;
        } while (--maxExt > 0);
        return ret;
    }
}



static void rebuildUnitigGraphAdjacency(Assembler& assembler, uint64_t vertexCount)
{
    assembler.unitigGraph.outgoing.clear();
    if (assembler.unitigGraph.incoming.isOpen()) {
        assembler.unitigGraph.incoming.clear();
    }

    assembler.unitigGraph.outgoing.beginPass1(uint32_t(vertexCount));
    if (assembler.unitigGraph.incoming.isOpen()) {
        assembler.unitigGraph.incoming.beginPass1(uint32_t(vertexCount));
    }

    for (uint64_t arcId = 0; arcId < assembler.unitigGraph.arcs.size(); ++arcId) {
        const auto& a = assembler.unitigGraph.arcs[arcId];
        if (a.del) continue;
        assembler.unitigGraph.outgoing.incrementCount(a.from);
        if (assembler.unitigGraph.incoming.isOpen()) {
            assembler.unitigGraph.incoming.incrementCount(a.to);
        }
    }

    assembler.unitigGraph.outgoing.beginPass2();
    if (assembler.unitigGraph.incoming.isOpen()) {
        assembler.unitigGraph.incoming.beginPass2();
    }
    for (uint64_t arcId = 0; arcId < assembler.unitigGraph.arcs.size(); ++arcId) {
        const auto& a = assembler.unitigGraph.arcs[arcId];
        if (a.del) continue;
        assembler.unitigGraph.outgoing.store(a.from, uint32_t(arcId));
        if (assembler.unitigGraph.incoming.isOpen()) {
            assembler.unitigGraph.incoming.store(a.to, uint32_t(arcId));
        }
    }
    assembler.unitigGraph.outgoing.endPass2();
    if (assembler.unitigGraph.incoming.isOpen()) {
        assembler.unitigGraph.incoming.endPass2();
    }

    // Sort outgoing adjacency by `len` (required by transitive reduction).
    for (uint32_t v = 0; v < vertexCount; ++v) {
        uint32_t* b = assembler.unitigGraph.outgoing.begin(v);
        uint32_t* e = assembler.unitigGraph.outgoing.end(v);
        std::sort(b, e, [&](uint32_t aId, uint32_t bId) {
            const auto& a = assembler.unitigGraph.arcs[aId];
            const auto& bArc = assembler.unitigGraph.arcs[bId];
            if (a.len != bArc.len) return a.len < bArc.len;
            return a.to < bArc.to;
        });
    }
}



static uint64_t unitigGraphTransitiveReduce(Assembler& assembler, uint32_t fuzz, uint32_t vertexCount)
{
    UnitigGraph& g = assembler.unitigGraph;
    rebuildUnitigGraphAdjacency(assembler, vertexCount);

    vector<uint8_t> mark(vertexCount, 0);
    uint64_t reduced = 0;

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.unitigDeleted.isOpen && g.unitigDeleted[v >> 1U]) {
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
        rebuildUnitigGraphAdjacency(assembler, vertexCount);
    }
    return reduced;
}



static uint64_t unitigGraphCutTips(Assembler& assembler, uint32_t maxShortTipUnitigs, uint32_t vertexCount)
{
    UnitigGraph& g = assembler.unitigGraph;
    if (maxShortTipUnitigs == 0) return 0;

    uint64_t cutCount = 0;
    vector<uint32_t> path;
    path.reserve(maxShortTipUnitigs + 4);

    rebuildUnitigGraphAdjacency(assembler, vertexCount);

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.unitigDeleted.isOpen && g.unitigDeleted[v >> 1U]) continue;

        const int ret = unitigGraphExtend(g, v, int(maxShortTipUnitigs), path);
        if (ret == ASG_ET_MERGEABLE) continue;
        if (uint32_t(path.size()) > maxShortTipUnitigs) continue;

        if (ret == ASG_ET_TIP || ret == ASG_ET_MULTI_OUT || ret == ASG_ET_MULTI_NEI) {
            for (size_t i = 0; i + 1 < path.size(); ++i) {
                const uint32_t from = path[i];
                const uint32_t to = path[i + 1];
                for (const uint32_t arcId : g.outgoing[from]) {
                    if (g.arcs[arcId].del) continue;
                    if (g.arcs[arcId].to == to) {
                        g.arcs[arcId].del = 1;
                        g.arcs[arcId ^ 1U].del = 1;
                        ++cutCount;
                        break;
                    }
                }
            }
        }
    }

    if (cutCount) {
        symmetrizeArcDeletion(g);
        rebuildUnitigGraphAdjacency(assembler, vertexCount);
    }
    return cutCount;
}



static uint64_t unitigGraphRemoveOneStepBubbles(Assembler& assembler, uint32_t vertexCount)
{
    UnitigGraph& g = assembler.unitigGraph;
    rebuildUnitigGraphAdjacency(assembler, vertexCount);

    vector<uint32_t> inDegree(vertexCount, 0);
    for (uint32_t v = 0; v < vertexCount; ++v) {
        for (const uint32_t arcId : g.outgoing[v]) {
            if (g.arcs[arcId].del) continue;
            ++inDegree[g.arcs[arcId].to];
        }
    }

    uint64_t removed = 0;
    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.unitigDeleted.isOpen && g.unitigDeleted[v >> 1U]) continue;

        vector<uint32_t> outs;
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) outs.push_back(arcId);
        }
        if (outs.size() != 2) continue;

        const uint32_t a = g.arcs[outs[0]].to;
        const uint32_t b = g.arcs[outs[1]].to;
        if (a == b) continue;
        if (inDegree[a] != 1 || inDegree[b] != 1) continue;

        if (countOutgoingNonDeleted(g, a) != 1) continue;
        if (countOutgoingNonDeleted(g, b) != 1) continue;

        uint32_t arcA = 0, arcB = 0;
        (void)firstOutgoingNonDeleted(g, a, arcA);
        (void)firstOutgoingNonDeleted(g, b, arcB);
        const uint32_t tA = g.arcs[arcA].to;
        const uint32_t tB = g.arcs[arcB].to;
        if (tA != tB) continue;

        // Pick the weaker branch by total overlap length.
        const uint64_t scoreA = uint64_t(g.arcs[outs[0]].overlapLen) + uint64_t(g.arcs[arcA].overlapLen);
        const uint64_t scoreB = uint64_t(g.arcs[outs[1]].overlapLen) + uint64_t(g.arcs[arcB].overlapLen);
        const uint32_t killFirst = (scoreA < scoreB) ? outs[0] : outs[1];
        const uint32_t killMid = (scoreA < scoreB) ? arcA : arcB;

        if (!g.arcs[killFirst].del) {
            g.arcs[killFirst].del = 1;
            g.arcs[killFirst ^ 1U].del = 1;
            ++removed;
        }
        if (!g.arcs[killMid].del) {
            g.arcs[killMid].del = 1;
            g.arcs[killMid ^ 1U].del = 1;
            ++removed;
        }
    }

    if (removed) {
        symmetrizeArcDeletion(g);
        rebuildUnitigGraphAdjacency(assembler, vertexCount);
    }
    return removed;
}



uint64_t Assembler::cleanUnitigGraphBreakShortCycles(uint32_t maxCycleUnitigs)
{
    checkUnitigGraphIsOpen();
    UnitigGraph& g = unitigGraph;
    const uint32_t vertexCount = uint32_t(2 * g.unitigs.size());

    rebuildUnitigGraphAdjacency(*this, vertexCount);

    vector<uint8_t> visited(vertexCount, 0);
    uint64_t removed = 0;

    for (uint32_t start = 0; start < vertexCount; ++start) {
        if (visited[start]) continue;
        if (g.unitigDeleted.isOpen && g.unitigDeleted[start >> 1U]) {
            visited[start] = 1;
            continue;
        }
        if (countOutgoingNonDeleted(g, start) != 1) {
            visited[start] = 1;
            continue;
        }
        if (countOutgoingNonDeleted(g, start ^ 1U) != 1) {
            visited[start] = 1;
            continue;
        }

        vector<uint32_t> path;
        path.reserve(maxCycleUnitigs + 2);
        unordered_map<uint32_t, size_t> index;
        index.reserve(maxCycleUnitigs + 2);

        uint32_t v = start;
        while (true) {
            if (index.contains(v)) {
                const size_t i0 = index[v];
                const size_t cycleSize = path.size() - i0;
                if (cycleSize > 0 && cycleSize <= maxCycleUnitigs) {
                    const uint32_t from = path.back();
                    uint32_t arcId = 0;
                    if (firstOutgoingNonDeleted(g, from, arcId)) {
                        g.arcs[arcId].del = 1;
                        g.arcs[arcId ^ 1U].del = 1;
                        removed += 2;
                    }
                }
                break;
            }
            if (path.size() > maxCycleUnitigs + 1) break;
            index[v] = path.size();
            path.push_back(v);

            uint32_t arcId = 0;
            if (!firstOutgoingNonDeleted(g, v, arcId)) break;
            const uint32_t next = g.arcs[arcId].to;

            if (countOutgoingNonDeleted(g, next) != 1) break;
            if (countOutgoingNonDeleted(g, next ^ 1U) != 1) break;
            v = next;
        }

        for (const uint32_t u : path) {
            visited[u] = 1;
        }
    }

    if (removed) {
        symmetrizeArcDeletion(g);
        rebuildUnitigGraphAdjacency(*this, vertexCount);
    }
    return removed;
}



void Assembler::cleanUnitigGraphInitialHifiasm(uint32_t gapFuzz, uint32_t maxShortTipUnitigs)
{
    checkUnitigGraphIsOpen();
    const uint32_t vertexCount = uint32_t(2 * unitigGraph.unitigs.size());

    cout << timestamp << "Unitig graph clean begins." << endl;
    const auto t0 = steady_clock::now();

    const uint64_t reduced = unitigGraphTransitiveReduce(*this, gapFuzz, vertexCount);
    const uint64_t cut = unitigGraphCutTips(*this, maxShortTipUnitigs, vertexCount);

    cout << timestamp << "Unitig graph clean complete in " << seconds(steady_clock::now() - t0)
         << " s (transReduced=" << reduced << ", cutTips=" << cut << ")" << endl;
}



void Assembler::cleanUnitigGraphPreCleanHifiasm(uint32_t maxShortTipUnitigs)
{
    checkUnitigGraphIsOpen();
    const uint32_t vertexCount = uint32_t(2 * unitigGraph.unitigs.size());

    cout << timestamp << "Unitig graph pre-clean begins." << endl;
    const auto t0 = steady_clock::now();

    uint64_t total = 0;
    for (uint32_t round = 0; round < 10; ++round) {
        const uint64_t broken = cleanUnitigGraphBreakShortCycles(10);
        const uint64_t bubbles = unitigGraphRemoveOneStepBubbles(*this, vertexCount);
        const uint64_t tips = unitigGraphCutTips(*this, maxShortTipUnitigs, vertexCount);
        const uint64_t changed = broken + bubbles + tips;
        total += changed;
        if (changed == 0) break;
    }

    cout << timestamp << "Unitig graph pre-clean complete in " << seconds(steady_clock::now() - t0)
         << " s (changes=" << total << ")" << endl;
}



void Assembler::cleanUnitigGraphDropShortOverlaps(double dropRatio, uint32_t minOverlapLen)
{
    checkUnitigGraphIsOpen();
    UnitigGraph& g = unitigGraph;
    const uint32_t vertexCount = uint32_t(2 * g.unitigs.size());

    rebuildUnitigGraphAdjacency(*this, vertexCount);

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (g.unitigDeleted.isOpen && g.unitigDeleted[v >> 1U]) continue;

        vector<uint32_t> out;
        out.reserve(g.outgoing.size(v));
        for (const uint32_t arcId : g.outgoing[v]) {
            if (!g.arcs[arcId].del) out.push_back(arcId);
        }
        if (out.size() < 2) continue;

        std::sort(out.begin(), out.end(), [&](uint32_t aId, uint32_t bId) {
            return g.arcs[aId].overlapLen > g.arcs[bId].overlapLen;
        });

        const uint32_t maxOl = g.arcs[out.front()].overlapLen;
        const uint32_t threshold = std::max<uint32_t>(uint32_t(double(maxOl) * dropRatio), minOverlapLen);
        for (size_t i = 1; i < out.size(); ++i) {
            const uint32_t arcId = out[i];
            if (g.arcs[arcId].overlapLen < threshold) {
                g.arcs[arcId].del = 1;
                g.arcs[arcId ^ 1U].del = 1;
            }
        }
    }

    symmetrizeArcDeletion(g);
    rebuildUnitigGraphAdjacency(*this, vertexCount);
}



void Assembler::cleanUnitigGraphDropOverlapRoundsHifiasm(
    uint32_t cleanRounds,
    double minDropRate,
    double maxDropRate,
    uint32_t maxShortTipUnitigs,
    uint32_t finalMinOverlapLen)
{
    checkUnitigGraphIsOpen();

    cout << timestamp << "Unitig graph overlap-drop rounds begin (cleanRounds=" << cleanRounds
         << ", minDropRate=" << minDropRate << ", maxDropRate=" << maxDropRate << ")." << endl;
    const auto t0 = steady_clock::now();

    for (uint32_t r = 0; r < cleanRounds; ++r) {
        const double alpha = (cleanRounds <= 1) ? 1.0 : double(r) / double(cleanRounds - 1);
        const double drop = minDropRate + (maxDropRate - minDropRate) * alpha;
        cleanUnitigGraphDropShortOverlaps(drop, 0);
        cleanUnitigGraphPreCleanHifiasm(maxShortTipUnitigs);
    }

    cleanUnitigGraphDropShortOverlaps(minDropRate, finalMinOverlapLen);
    cleanUnitigGraphPreCleanHifiasm(maxShortTipUnitigs);

    cout << timestamp << "Unitig graph overlap-drop rounds complete in " << seconds(steady_clock::now() - t0) << " s" << endl;
}

