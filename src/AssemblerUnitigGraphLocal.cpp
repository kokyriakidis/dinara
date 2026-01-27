// Dinara.
#include "Assembler.hpp"
#include "LocalUnitigGraph.hpp"
#include "chrono.hpp"
using namespace dinara;
using namespace std;

// Standard library.
#include <queue>



bool Assembler::createLocalUnitigGraph(
    const vector<OrientedUnitigId>& starts,
    uint32_t maxDistance,
    bool followOutgoing,
    bool followIncoming,
    double timeout,
    LocalUnitigGraph& graph)
{
    checkUnitigGraphIsOpen();

    const auto startTime = steady_clock::now();
    std::queue<OrientedUnitigId> q;

    const auto unitigIsDeleted = [&](const OrientedUnitigId& u) -> bool {
        const uint32_t id = u.getUnitigId();
        return unitigGraph.unitigDeleted.isOpen &&
            id < unitigGraph.unitigDeleted.size() &&
            unitigGraph.unitigDeleted[id];
    };

    auto addStartIfAllowed = [&](OrientedUnitigId orientedUnitigId)
    {
        if (unitigIsDeleted(orientedUnitigId)) return;
        if (graph.vertexExists(orientedUnitigId)) return;
        const uint32_t unitigId = orientedUnitigId.getUnitigId();
        const bool isCircular = unitigGraph.unitigs[unitigId].circ != 0;
        const uint32_t readCount = unitigGraph.unitigs[unitigId].n;
        graph.addVertex(orientedUnitigId, readCount, isCircular, 0);
        q.push(orientedUnitigId);
    };

    for (const OrientedUnitigId start : starts) {
        addStartIfAllowed(start);
    }

    auto addNeighborIfAllowed = [&](OrientedUnitigId neighbor, uint32_t distance)
    {
        if (unitigIsDeleted(neighbor)) return false;
        if (graph.vertexExists(neighbor)) return false;
        const uint32_t unitigId = neighbor.getUnitigId();
        const bool isCircular = unitigGraph.unitigs[unitigId].circ != 0;
        const uint32_t readCount = unitigGraph.unitigs[unitigId].n;
        graph.addVertex(neighbor, readCount, isCircular, distance);
        q.push(neighbor);
        return true;
    };

    auto processArc = [&](uint32_t arcId, OrientedUnitigId current, uint32_t currentDistance)
    {
        const UnitigGraphArc& arc = unitigGraph.arcs[arcId];
        const OrientedUnitigId from = OrientedUnitigId::fromValue(arc.from);
        const OrientedUnitigId to = OrientedUnitigId::fromValue(arc.to);

        const uint32_t currentValue = current.getValue();
        const OrientedUnitigId neighbor = (arc.from == currentValue) ? to : from;

        if (!graph.vertexExists(neighbor)) {
            addNeighborIfAllowed(neighbor, currentDistance + 1);
        }

        if (graph.vertexExists(from) && graph.vertexExists(to)) {
            graph.addEdge(from, to, arc.overlapLen, arc.len, arcId);
        }
    };

    while (!q.empty()) {
        if (timeout > 0) {
            const double elapsed = seconds(steady_clock::now() - startTime);
            if (elapsed > timeout) {
                return false;
            }
        }

        const OrientedUnitigId current = q.front();
        q.pop();
        const uint32_t currentDistance = graph.getDistance(current);
        if (currentDistance >= maxDistance) continue;
        const uint32_t currentValue = current.getValue();
        if (followOutgoing) {
            for (const uint32_t arcId : unitigGraph.outgoing[currentValue]) {
                processArc(arcId, current, currentDistance);
            }
        }
        if (followIncoming) {
            for (const uint32_t arcId : unitigGraph.incoming[currentValue]) {
                processArc(arcId, current, currentDistance);
            }
        }
    }

    return true;
}
