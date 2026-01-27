// Dinara.
#include "Assembler.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"

using namespace dinara;
using namespace std;

namespace {
    inline uint32_t countOutgoing(const StringGraph& g, uint32_t v)
    {
        return uint32_t(g.outgoing.size(v));
    }

    inline uint32_t arcFirst(const StringGraph& g, uint32_t v)
    {
        DINARA_ASSERT(g.outgoing.size(v) == 1);
        return g.outgoing[v][0];
    }
}



void Assembler::createUnitigGraphFromStringGraph()
{
    checkStringGraphIsOpen();
    reads->checkReadsAreOpen();

    const uint64_t readCount = reads->readCount();
    const uint32_t stringVertexCount = uint32_t(2 * readCount);

    auto effectiveReadLength = [&](ReadId readId) -> uint32_t {
        if (validReadIntervals.empty()) {
            return uint32_t(reads->getReadRawSequenceLength(readId));
        }
        const auto& r = validReadIntervals[readId];
        return r.end > r.start ? (r.end - r.start) : 0U;
    };

    // Recreate from scratch.
    unitigGraph.remove();
    unitigGraph.unitigs.createNew(largeDataName("UnitigGraphUnitigs"), largeDataPageSize);
    unitigGraph.unitigPaths.createNew(largeDataName("UnitigGraphPaths"), largeDataPageSize);
    unitigGraph.arcs.createNew(largeDataName("UnitigGraphArcs"), largeDataPageSize);
    unitigGraph.outgoing.createNew(largeDataName("UnitigGraphOutgoing"), largeDataPageSize);
    unitigGraph.incoming.createNew(largeDataName("UnitigGraphIncoming"), largeDataPageSize);
    unitigGraph.unitigDeleted.createNew(largeDataName("UnitigGraphUnitigDeleted"), largeDataPageSize);

    vector<int32_t> mark(stringVertexCount, 0);
    vector<uint64_t> path;
    path.reserve(64);

    // Build unitigs (hifiasm `ma_ug_gen` logic adapted to our StringGraph representation).
    for (uint32_t v = 0; v < stringVertexCount; ++v) {
        if (stringGraph.readDeleted.isOpen && stringGraph.readDeleted[v >> 1U]) continue;
        if (mark[v]) continue;
        if (countOutgoing(stringGraph, v) == 0 && countOutgoing(stringGraph, v ^ 1U) != 0) continue;

        mark[v] = 1;
        path.clear();

        uint32_t start = v;
        uint32_t end = v ^ 1U;
        uint64_t len = 0;

        // Forward extension.
        uint32_t w = v;
        while (true) {
            if (countOutgoing(stringGraph, w) != 1) break;
            const uint32_t arcId = arcFirst(stringGraph, w);
            const uint32_t x = stringGraph.arcs[arcId].to;
            if (countOutgoing(stringGraph, x ^ 1U) != 1) break;

            mark[x] = 1;
            mark[w ^ 1U] = 1;

            const uint32_t l = stringGraph.arcs[arcId].len;
            path.push_back((uint64_t(w) << 32) | uint64_t(l));
            end = x ^ 1U;
            len += l;
            w = x;
            if (x == v) break; // circular
        }

        if (start != (end ^ 1U) || path.empty()) {
            // Linear unitig: add the last vertex (end^1) with its node length.
            const uint32_t lastVertex = end ^ 1U;
            const uint32_t l = effectiveReadLength(ReadId(lastVertex >> 1U));
            path.push_back((uint64_t(lastVertex) << 32) | uint64_t(l));
            len += l;
        } else {
            // Circular unitig: skip backward extension.
            start = end = std::numeric_limits<uint32_t>::max();
        }

        // Backward extension (only for linear unitigs).
        if (start != std::numeric_limits<uint32_t>::max()) {
            uint32_t x = v;
            while (true) {
                if (countOutgoing(stringGraph, x ^ 1U) != 1) break;
                const uint32_t arcIdX = arcFirst(stringGraph, x ^ 1U);
                const uint32_t pred = stringGraph.arcs[arcIdX].to ^ 1U; // w->x
                if (countOutgoing(stringGraph, pred) != 1) break;

                mark[x] = 1;
                mark[pred ^ 1U] = 1;

                const uint32_t arcIdPred = arcFirst(stringGraph, pred);
                const uint32_t l = stringGraph.arcs[arcIdPred].len;
                path.insert(path.begin(), (uint64_t(pred) << 32) | uint64_t(l));
                start = pred;
                len += l;
                x = pred;
            }
            mark[start] = 1;
            mark[end] = 1;
        }

        UnitigInfo info;
        info.start = start;
        info.end = end;
        info.len = len;
        info.n = uint32_t(path.size());
        info.circ = uint8_t(start == std::numeric_limits<uint32_t>::max());

        const uint32_t unitigId = uint32_t(unitigGraph.unitigs.size());
        unitigGraph.unitigs.push_back(info);
        unitigGraph.unitigPaths.appendVector(path);
        (void)unitigId;
    }

    const uint32_t unitigCount = uint32_t(unitigGraph.unitigs.size());
    unitigGraph.unitigDeleted.reserveAndResize(unitigCount);
    std::fill(unitigGraph.unitigDeleted.begin(), unitigGraph.unitigDeleted.end(), uint8_t(0));

    // Build unitig-graph arcs (hifiasm `ma_ug_gen` second phase).
    std::fill(mark.begin(), mark.end(), -1);
    for (uint32_t i = 0; i < unitigCount; ++i) {
        const UnitigInfo& u = unitigGraph.unitigs[i];
        if (u.circ) continue;
        if (u.start != std::numeric_limits<uint32_t>::max()) {
            mark[u.start] = int32_t((i << 1) | 0);
            mark[u.end] = int32_t((i << 1) | 1);
        }
    }

    // Each StringGraph arc yields at most one unitig-graph arc; we create its RC twin explicitly.
    for (uint64_t stringArcId = 0; stringArcId + 1 < stringGraph.arcs.size(); stringArcId += 2) {
        const auto& a = stringGraph.arcs[stringArcId];
        if (a.del) continue;

        const int32_t fromMark = mark[a.from ^ 1U];
        const int32_t toMark = mark[a.to];
        if (fromMark < 0 || toMark < 0) {
            continue;
        }

        const uint32_t u = uint32_t(fromMark) ^ 1U;
        const uint32_t v = uint32_t(toMark);
        const uint32_t unitigId = u >> 1U;
        const uint64_t unitigLen = unitigGraph.unitigs[unitigId].len;
        int64_t l = int64_t(unitigLen) - int64_t(a.overlapLen);
        if (l < 0) l = 1;

        UnitigGraphArc ua;
        ua.from = u;
        ua.to = v;
        ua.overlapLen = a.overlapLen;
        ua.len = uint32_t(l);
        ua.stringArcId = stringArcId;
        ua.del = 0;
        unitigGraph.arcs.push_back(ua);

        UnitigGraphArc ub = ua;
        ub.from = ua.to ^ 1U;
        ub.to = ua.from ^ 1U;
        ub.stringArcId = stringArcId ^ 1ULL;
        unitigGraph.arcs.push_back(ub);
    }

    const uint32_t unitigVertexCount = 2 * unitigCount;

    // Build adjacency.
    unitigGraph.outgoing.beginPass1(unitigVertexCount);
    unitigGraph.incoming.beginPass1(unitigVertexCount);
    for (uint64_t arcId = 0; arcId < unitigGraph.arcs.size(); ++arcId) {
        const UnitigGraphArc& arc = unitigGraph.arcs[arcId];
        if (arc.del) continue;
        unitigGraph.outgoing.incrementCount(arc.from);
        unitigGraph.incoming.incrementCount(arc.to);
    }
    unitigGraph.outgoing.beginPass2();
    unitigGraph.incoming.beginPass2();
    for (uint64_t arcId = 0; arcId < unitigGraph.arcs.size(); ++arcId) {
        const UnitigGraphArc& arc = unitigGraph.arcs[arcId];
        if (arc.del) continue;
        unitigGraph.outgoing.store(arc.from, uint32_t(arcId));
        unitigGraph.incoming.store(arc.to, uint32_t(arcId));
    }
    unitigGraph.outgoing.endPass2();
    unitigGraph.incoming.endPass2();

    unitigGraph.unreserve();

    cout << timestamp << "Unitig graph: " << unitigCount
         << " unitigs, " << unitigGraph.arcs.size() << " arcs." << endl;
}



void Assembler::accessUnitigGraph()
{
    unitigGraph.unitigs.accessExistingReadOnly(largeDataName("UnitigGraphUnitigs"));
    unitigGraph.unitigPaths.accessExistingReadOnly(largeDataName("UnitigGraphPaths"));
    unitigGraph.arcs.accessExistingReadOnly(largeDataName("UnitigGraphArcs"));
    unitigGraph.outgoing.accessExistingReadOnly(largeDataName("UnitigGraphOutgoing"));
    unitigGraph.incoming.accessExistingReadOnly(largeDataName("UnitigGraphIncoming"));
    unitigGraph.unitigDeleted.accessExistingReadOnly(largeDataName("UnitigGraphUnitigDeleted"));
}



void Assembler::accessUnitigGraphReadWrite()
{
    unitigGraph.unitigs.accessExistingReadWrite(largeDataName("UnitigGraphUnitigs"));
    unitigGraph.unitigPaths.accessExistingReadWrite(largeDataName("UnitigGraphPaths"));
    unitigGraph.arcs.accessExistingReadWrite(largeDataName("UnitigGraphArcs"));
    unitigGraph.outgoing.accessExistingReadWrite(largeDataName("UnitigGraphOutgoing"));
    unitigGraph.incoming.accessExistingReadWrite(largeDataName("UnitigGraphIncoming"));
    unitigGraph.unitigDeleted.accessExistingReadWrite(largeDataName("UnitigGraphUnitigDeleted"));
}



void Assembler::checkUnitigGraphIsOpen() const
{
    if (!unitigGraph.unitigs.isOpen) {
        throw runtime_error("Unitig graph is not accessible.");
    }
    if (!unitigGraph.arcs.isOpen) {
        throw runtime_error("Unitig graph arcs are not accessible.");
    }
    if (!unitigGraph.outgoing.isOpen()) {
        throw runtime_error("Unitig graph outgoing adjacency is not accessible.");
    }
}

