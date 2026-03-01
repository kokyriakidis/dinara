// Dinara.
#include "UnitigGraph.hpp"

#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>

using namespace dinara;


void UnitigGraph::writeGfa(const std::string& fileName) const
{
    std::ofstream gfa(fileName);
    if (!gfa) {
        throw std::runtime_error("Cannot open " + fileName + " for writing.");
    }
    writeGfa(gfa);
}



void UnitigGraph::writeGfa(std::ostream& gfa) const
{
    if (!unitigs.isOpen) {
        throw std::runtime_error("UnitigGraph::writeGfa: unitigs are not accessible.");
    }
    if (!arcs.isOpen) {
        throw std::runtime_error("UnitigGraph::writeGfa: arcs are not accessible.");
    }

    gfa << "H\tVN:Z:1.0\n";

    const uint32_t unitigCount = uint32_t(unitigs.size());

    // S lines — one per unitig.
    // Name: utg<i>  Sequence: *  Tags: LN:i:<len>  RC:i:<n>
    for (uint32_t i = 0; i < unitigCount; ++i) {
        if (unitigDeleted.isOpen && unitigDeleted[i]) continue;
        const UnitigInfo& u = unitigs[i];
        gfa << "S\tutg" << i << "\t*\tLN:i:" << u.len << "\tRC:i:" << u.n << "\n";
    }

    // L lines — one per canonical arc (even arcId).
    for (uint64_t arcId = 0; arcId + 1 < arcs.size(); arcId += 2) {
        const UnitigGraphArc& arc = arcs[arcId];
        if (arc.del) continue;

        const uint32_t fromUtig = arc.from >> 1U;
        const uint32_t toUtig   = arc.to   >> 1U;
        if (unitigDeleted.isOpen && (unitigDeleted[fromUtig] || unitigDeleted[toUtig])) continue;

        const char fromOrient = (arc.from & 1U) ? '-' : '+';
        const char toOrient   = (arc.to   & 1U) ? '-' : '+';

        gfa << "L\tutg" << fromUtig << "\t" << fromOrient
            << "\tutg" << toUtig   << "\t" << toOrient
            << "\t" << arc.overlapLen << "M\n";
    }
}



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
