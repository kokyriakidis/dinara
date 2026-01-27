// Dinara.
#include "UnitigGraph.hpp"

using namespace dinara;


void UnitigGraph::unreserve()
{
    if (unitigs.isOpenWithWriteAccess) {
        unitigs.unreserve();
    }
    if (unitigPaths.isOpenWithWriteAccess()) {
        unitigPaths.unreserve();
    }
    if (arcs.isOpenWithWriteAccess) {
        arcs.unreserve();
    }
    if (outgoing.isOpenWithWriteAccess()) {
        outgoing.unreserve();
    }
    if (incoming.isOpenWithWriteAccess()) {
        incoming.unreserve();
    }
    if (unitigDeleted.isOpenWithWriteAccess) {
        unitigDeleted.unreserve();
    }
}


void UnitigGraph::remove()
{
    if (unitigs.isOpen) {
        unitigs.remove();
    }
    if (unitigPaths.isOpen()) {
        unitigPaths.remove();
    }
    if (arcs.isOpen) {
        arcs.remove();
    }
    if (outgoing.isOpen()) {
        outgoing.remove();
    }
    if (incoming.isOpen()) {
        incoming.remove();
    }
    if (unitigDeleted.isOpen) {
        unitigDeleted.remove();
    }
}
