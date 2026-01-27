// Dinara.
#include "StringGraph.hpp"
#include "DINARA_ASSERT.hpp"

using namespace dinara;

uint64_t StringGraph::getReverseComplementArcId(uint64_t arcId) const
{
    const uint64_t rcId = arcId ^ 1ULL;
    const StringGraphArc& a = arcs[arcId];
    const StringGraphArc& b = arcs[rcId];
    DINARA_ASSERT((a.from ^ 1U) == b.to);
    DINARA_ASSERT((a.to ^ 1U) == b.from);
    DINARA_ASSERT(a.len == b.len);
    DINARA_ASSERT(a.overlapLen == b.overlapLen);
    return rcId;
}

void StringGraph::unreserve()
{
    if (arcs.isOpenWithWriteAccess) {
        arcs.unreserve();
    }
    if (outgoing.isOpenWithWriteAccess()) {
        outgoing.unreserve();
    }
    if (incoming.isOpenWithWriteAccess()) {
        incoming.unreserve();
    }
    if (readDeleted.isOpenWithWriteAccess) {
        readDeleted.unreserve();
    }
}

void StringGraph::remove()
{
    if (arcs.isOpen) {
        arcs.remove();
    }
    if (outgoing.isOpen()) {
        outgoing.remove();
    }
    if (incoming.isOpen()) {
        incoming.remove();
    }
    if (readDeleted.isOpen) {
        readDeleted.remove();
    }
}
