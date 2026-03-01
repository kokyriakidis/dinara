// Dinara.
#include "StringGraph.hpp"
#include "DINARA_ASSERT.hpp"
#include "Reads.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <ostream>
#include <stdexcept>
#include <unordered_set>
#include <unordered_map>
#include <vector>

using namespace dinara;

uint64_t StringGraph::getReverseComplementArcId(uint64_t arcId) const
{
    const uint64_t rcId = arcId ^ 1ULL;
    const StringGraphArc& a = arcs[arcId];
    const StringGraphArc& b = arcs[rcId];
    DINARA_ASSERT((a.from ^ 1U) == b.to);
    DINARA_ASSERT((a.to ^ 1U) == b.from);
    return rcId;
}

void StringGraph::writeGfa(const std::string& fileName) const
{
    std::ofstream gfa(fileName);
    if (!gfa) {
        throw std::runtime_error("Cannot open " + fileName + " for writing.");
    }
    writeGfa(gfa);
}

void StringGraph::writeGfa(std::ostream& gfa) const
{
    if (!arcs.isOpen) {
        throw std::runtime_error("StringGraph::writeGfa: arcs are not accessible.");
    }

    gfa << "H\tVN:Z:1.0\n";

    // Collect read ids that appear in at least one non-deleted canonical arc (even arcId).
    std::unordered_set<uint32_t> readIds;
    readIds.reserve(size_t(arcs.size() / 2));

    for (uint64_t arcId = 0; arcId + 1 < arcs.size(); arcId += 2) {
        const StringGraphArc& arc = arcs[arcId];
        if (arc.del) {
            continue;
        }
        readIds.insert(arc.from >> 1U);
        readIds.insert(arc.to >> 1U);
    }

    // S lines — one per read id that appears in the graph.
    std::vector<uint32_t> readIdsSorted(readIds.begin(), readIds.end());
    std::sort(readIdsSorted.begin(), readIdsSorted.end());
    for (const uint32_t readId : readIdsSorted) {
        gfa << "S\t" << readId << "\t*\n";
    }

    // L lines — one per canonical arc (even arcId).
    for (uint64_t arcId = 0; arcId + 1 < arcs.size(); arcId += 2) {
        const StringGraphArc& arc = arcs[arcId];
        if (arc.del) {
            continue;
        }

        const uint32_t fromRead = arc.from >> 1U;
        const uint32_t toRead = arc.to >> 1U;
        const char fromOrient = (arc.from & 1U) ? '-' : '+';
        const char toOrient = (arc.to & 1U) ? '-' : '+';

        gfa << "L\t" << fromRead << "\t" << fromOrient
            << "\t" << toRead << "\t" << toOrient
            << "\t" << arc.overlapLen << "M\n";
    }
}

namespace {
    inline std::string sanitizeGfaName(span<const char> name)
    {
        std::string s;
        s.reserve(name.size());
        for (const char c : name) {
            const unsigned char uc = static_cast<unsigned char>(c);
            if (uc < 33 || uc == 127 || std::isspace(uc)) {
                s.push_back('_');
            } else {
                s.push_back(c);
            }
        }
        // Segment names cannot be empty or "*".
        if (s.empty() || s == "*") {
            s.clear();
        }
        return s;
    }
}

void StringGraph::writeGfa(const std::string& fileName, const Reads& reads) const
{
    std::ofstream gfa(fileName);
    if (!gfa) {
        throw std::runtime_error("Cannot open " + fileName + " for writing.");
    }
    writeGfa(gfa, reads);
}

void StringGraph::writeGfa(std::ostream& gfa, const Reads& reads) const
{
    if (!arcs.isOpen) {
        throw std::runtime_error("StringGraph::writeGfa: arcs are not accessible.");
    }

    gfa << "H\tVN:Z:1.0\n";

    // Collect read ids that appear in at least one non-deleted canonical arc (even arcId).
    std::unordered_set<uint32_t> readIds;
    readIds.reserve(size_t(arcs.size() / 2));
    for (uint64_t arcId = 0; arcId + 1 < arcs.size(); arcId += 2) {
        const StringGraphArc& arc = arcs[arcId];
        if (arc.del) {
            continue;
        }
        readIds.insert(arc.from >> 1U);
        readIds.insert(arc.to >> 1U);
    }

    // Stable, deterministic ordering.
    std::vector<uint32_t> readIdsSorted(readIds.begin(), readIds.end());
    std::sort(readIdsSorted.begin(), readIdsSorted.end());

    // Map readId -> segment name (sanitized and made unique).
    std::unordered_map<uint32_t, std::string> segmentName;
    segmentName.reserve(readIdsSorted.size() * 2);
    std::unordered_set<std::string> usedNames;
    usedNames.reserve(readIdsSorted.size() * 2);

    const uint32_t readCount = uint32_t(reads.readCount());
    for (const uint32_t readId : readIdsSorted) {
        std::string name;
        if (readId < readCount) {
            name = sanitizeGfaName(reads.getReadName(ReadId(readId)));
        }
        if (name.empty()) {
            name = std::to_string(readId);
        }
        if (!usedNames.insert(name).second) {
            name += "_" + std::to_string(readId);
            usedNames.insert(name);
        }
        segmentName.emplace(readId, std::move(name));
    }

    // S lines.
    for (const uint32_t readId : readIdsSorted) {
        gfa << "S\t" << segmentName.at(readId) << "\t*\n";
    }

    // L lines — one per canonical arc (even arcId).
    for (uint64_t arcId = 0; arcId + 1 < arcs.size(); arcId += 2) {
        const StringGraphArc& arc = arcs[arcId];
        if (arc.del) {
            continue;
        }

        const uint32_t fromRead = arc.from >> 1U;
        const uint32_t toRead = arc.to >> 1U;
        const char fromOrient = (arc.from & 1U) ? '-' : '+';
        const char toOrient = (arc.to & 1U) ? '-' : '+';

        gfa << "L\t" << segmentName.at(fromRead) << "\t" << fromOrient
            << "\t" << segmentName.at(toRead) << "\t" << toOrient
            << "\t" << arc.overlapLen << "M\n";
    }
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
