// Dinara.
#include "Assembler.hpp"
#include "LocalStringGraph.hpp"
#include "chrono.hpp"
#include "Reads.hpp"
using namespace dinara;
using namespace std;

// Standard library.
#include <queue>



bool Assembler::createLocalStringGraph(
    const vector<OrientedReadId>& starts,
    uint32_t maxDistance,
    bool allowChimericReads,
    bool followOutgoing,
    bool followIncoming,
    double timeout,
    LocalStringGraph& graph)
{
    checkStringGraphIsOpen();
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();

    const auto startTime = steady_clock::now();
    std::queue<OrientedReadId> q;

    auto addStartIfAllowed = [&](OrientedReadId orientedReadId)
    {
        if (!allowChimericReads && reads->getFlags(orientedReadId.getReadId()).isChimeric) {
            return;
        }
        if (stringGraph.readDeleted.isOpen &&
            orientedReadId.getReadId() < stringGraph.readDeleted.size() &&
            stringGraph.readDeleted[orientedReadId.getReadId()]) {
            return;
        }
        if (graph.vertexExists(orientedReadId)) {
            return;
        }
        graph.addVertex(orientedReadId,
            uint32_t((*markers)[orientedReadId.getValue()].size()),
            reads->getFlags(orientedReadId.getReadId()).isChimeric,
            0);
        q.push(orientedReadId);
    };

    for (const OrientedReadId start : starts) {
        addStartIfAllowed(start);
    }

    auto addNeighborIfAllowed = [&](OrientedReadId neighbor, uint32_t distance)
    {
        if (!allowChimericReads && reads->getFlags(neighbor.getReadId()).isChimeric) {
            return false;
        }
        if (stringGraph.readDeleted.isOpen &&
            neighbor.getReadId() < stringGraph.readDeleted.size() &&
            stringGraph.readDeleted[neighbor.getReadId()]) {
            return false;
        }
        if (graph.vertexExists(neighbor)) {
            return false;
        }
        graph.addVertex(neighbor,
            uint32_t((*markers)[neighbor.getValue()].size()),
            reads->getFlags(neighbor.getReadId()).isChimeric,
            distance);
        q.push(neighbor);
        return true;
    };

    auto processArc = [&](uint32_t arcId, OrientedReadId current, uint32_t currentDistance)
    {
        const StringGraphArc& arc = stringGraph.arcs[arcId];
        const OrientedReadId from = OrientedReadId::fromValue(arc.from);
        const OrientedReadId to = OrientedReadId::fromValue(arc.to);

        // Determine which endpoint is the neighbor from the viewpoint of `current`.
        const uint32_t currentValue = uint32_t(current.getValue());
        const OrientedReadId neighbor = (arc.from == currentValue) ? to : from;

        if (!graph.vertexExists(neighbor)) {
            addNeighborIfAllowed(neighbor, currentDistance + 1);
        }

        // Add the directed arc to the local graph if both endpoints are present.
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

        const OrientedReadId current = q.front();
        q.pop();
        const uint32_t currentDistance = graph.getDistance(current);
        if (currentDistance >= maxDistance) {
            continue;
        }
        const uint32_t currentValue = uint32_t(current.getValue());

        if (followOutgoing) {
            for (const uint32_t arcId : stringGraph.outgoing[currentValue]) {
                processArc(arcId, current, currentDistance);
            }
        }
        if (followIncoming && stringGraph.incoming.isOpen()) {
            for (const uint32_t arcId : stringGraph.incoming[currentValue]) {
                processArc(arcId, current, currentDistance);
            }
        }
    }

    return true;
}
