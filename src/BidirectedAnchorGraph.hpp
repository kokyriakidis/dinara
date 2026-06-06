#pragma once

// BidirectedAnchorGraph.hpp
//
// Verkko/MBG-style bidirected anchor graph. Each node is a
// BidirectedAnchorId (one per locus, collapsing forward and RC).
// Edges connect OrientedAnchors. Like Verkko/MBG, the caller must
// add both traversal directions explicitly.
//
// Edge properties are stored once per canonical link (via canonAnchor()).
// Directional fields are stored relative to the canonical direction.
// getEdgeProperties() adjusts them for the query direction.
//
// Unitigification collapses linear chains of nodes into single segments,
// producing a compressed graph where each segment tracks its constituent
// anchors, window membership, and aggregate coverage.

#include "BidirectedAnchorId.hpp"

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace dinara {

struct AnchorWindow;
class Shasta2Journeys;

class BidirectedAnchorGraph {
public:
    // Edge properties stored once per canonical link.
    // Directional fields are relative to the canonical direction.
    struct EdgeProperties {
        // Symmetric fields.
        uint64_t coverage = 0;
        uint64_t offset = 0;  // Base-level distance between the two anchors.
        uint32_t sharedReadCount = 0;
        bool isInterWindow = false;
        bool useForAssembly = true;

        // Directional fields (relative to canonical direction).
        // supportingSpanPrev = span at the canonical "from" anchor.
        // supportingSpanNext = span at the canonical "to" anchor.
        uint32_t supportingSpanPrev = 0;
        uint32_t supportingSpanNext = 0;

        EdgeProperties swapped() const {
            EdgeProperties result = *this;
            std::swap(result.supportingSpanPrev, result.supportingSpanNext);
            return result;
        }
    };

    struct NodeProperties {
        uint32_t windowId = UINT32_MAX;
        uint32_t backboneReadId = UINT32_MAX;  // ReadId of the backbone read that owns this window.
    };

    // A unitig (compressed segment) — a linear chain of oriented anchors.
    struct Unitig {
        // Ordered chain of oriented anchors in this unitig.
        // For a single-node unitig, this has one entry.
        std::vector<OrientedAnchor> chain;

        // Window IDs traversed by this unitig (consecutive duplicates removed,
        // noWindow entries omitted).
        std::vector<uint32_t> windowSequence;

        // Average edge coverage along the chain (0 for single-node unitigs).
        double averageCoverage = 0.;

        // Total base-level length (sum of edge offsets along the chain).
        uint64_t totalOffset = 0;

        uint64_t anchorCount() const { return chain.size(); }
        OrientedAnchor front() const { return chain.front(); }
        OrientedAnchor back() const { return chain.back(); }
    };

    uint64_t numNodes() const { return nodeCount; }
    uint64_t numEdges() const { return edgeProperties.size(); }

    void resize(uint64_t n) {
        nodeCount = n;
        adjacency.resize(n);
        nodeProps.resize(n);
    }

    void setNodeWindow(BidirectedAnchorId id, uint32_t windowId) {
        nodeProps[id.value()].windowId = windowId;
    }

    void setNodeBackboneRead(BidirectedAnchorId id, uint32_t readId) {
        nodeProps[id.value()].backboneReadId = readId;
    }

    // Add a single directed traversal and store properties under
    // the canonical key. Caller must add both directions.
    void addEdge(OrientedAnchor from, OrientedAnchor to,
                 const EdgeProperties& props)
    {
        adjacency[from.first.value()][from.second].insert(to);
        auto key = canonAnchor(from, to);
        edgeProperties[key] = props;
    }

    // Add traversal without setting properties.
    void addTraversal(OrientedAnchor from, OrientedAnchor to)
    {
        adjacency[from.first.value()][from.second].insert(to);
    }

    bool hasEdge(OrientedAnchor from, OrientedAnchor to) const {
        auto idx = from.first.value();
        if(idx >= adjacency.size()) return false;
        auto it = adjacency[idx].find(from.second);
        if(it == adjacency[idx].end()) return false;
        return it->second.count(to) > 0;
    }

    std::vector<OrientedAnchor> getNeighbors(OrientedAnchor node) const {
        std::vector<OrientedAnchor> result;
        auto idx = node.first.value();
        if(idx >= adjacency.size()) return result;
        auto it = adjacency[idx].find(node.second);
        if(it == adjacency[idx].end()) return result;
        result.assign(it->second.begin(), it->second.end());
        return result;
    }

    // Get edge properties adjusted for the query direction.
    // Returns false if the edge doesn't exist.
    bool getEdgeProperties(OrientedAnchor from, OrientedAnchor to,
                           EdgeProperties& result) const
    {
        auto key = canonAnchor(from, to);
        auto it = edgeProperties.find(key);
        if(it == edgeProperties.end()) return false;
        if(isCanonicalAnchorDirection(from, to)) {
            result = it->second;
        } else {
            result = it->second.swapped();
        }
        return true;
    }

    // Mutable pointer to canonical storage. Caller must be aware
    // that directional fields are in the canonical direction.
    EdgeProperties* getEdgePropertiesMutable(OrientedAnchor from, OrientedAnchor to) {
        auto key = canonAnchor(from, to);
        auto it = edgeProperties.find(key);
        if(it == edgeProperties.end()) return nullptr;
        return &it->second;
    }

    // Const pointer to canonical storage (for existence checks).
    const EdgeProperties* getEdgePropertiesCanonical(OrientedAnchor from, OrientedAnchor to) const {
        auto key = canonAnchor(from, to);
        auto it = edgeProperties.find(key);
        if(it == edgeProperties.end()) return nullptr;
        return &it->second;
    }

    void removeEdge(OrientedAnchor from, OrientedAnchor to) {
        auto idx = from.first.value();
        if(idx < adjacency.size()) {
            adjacency[idx][from.second].erase(to);
        }
        edgeProperties.erase(canonAnchor(from, to));
    }

    uint64_t degree(BidirectedAnchorId id) const {
        auto idx = id.value();
        if(idx >= adjacency.size()) return 0;
        uint64_t d = 0;
        auto it = adjacency[idx].find(true);
        if(it != adjacency[idx].end()) d += it->second.size();
        it = adjacency[idx].find(false);
        if(it != adjacency[idx].end()) d += it->second.size();
        return d;
    }

    // Number of neighbors reachable from a specific side of a node.
    // sideDegree({v, true}) = number of edges leaving v+.
    uint64_t sideDegree(OrientedAnchor oa) const {
        auto idx = oa.first.value();
        if(idx >= adjacency.size()) return 0;
        auto it = adjacency[idx].find(oa.second);
        if(it == adjacency[idx].end()) return 0;
        return it->second.size();
    }

    // Remove an edge in both traversal directions (forward and RC mirror).
    void removeEdgeBothDirections(OrientedAnchor from, OrientedAnchor to) {
        removeEdge(from, to);
        removeEdge(reverseAnchor(to), reverseAnchor(from));
    }

    void normalizeOrientations(const std::vector<bool>& swapOrientation);
    void writeGfa(const std::string& fileName) const;
    void writeCsv(const std::string& fileName, uint32_t windowCount) const;

    // Collapse linear chains into unitigs. Returns the list of unitigs.
    // A node is internal to a chain when both sides have exactly one
    // neighbor and those neighbors reciprocally have degree 1 on the
    // entry side.
    std::vector<Unitig> unitigify() const;

    // Write unitig GFA where each segment is a unitig.
    void writeUnitigGfa(const std::string& fileName,
                        const std::vector<Unitig>& unitigs,
                        uint32_t windowCount) const;

    // Write unitig CSV for Bandage coloring.
    void writeUnitigCsv(const std::string& fileName,
                        const std::vector<Unitig>& unitigs,
                        uint32_t windowCount) const;

    // Write CSV colored by backbone read ID (windows from the same read
    // get the same color).
    void writeCsvByRead(const std::string& fileName, uint32_t readCount) const;
    void writeUnitigCsvByRead(const std::string& fileName,
                              const std::vector<Unitig>& unitigs,
                              uint32_t readCount) const;

    // For each window pair, keep only the exit edge at the last backbone
    // position and the entry edge at the first backbone position.
    // Removes redundant inter-window edges that fragment backbones.
    uint64_t filterInterWindowEdges(
        const std::vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys);

    // Bypass detour filter: when a neighbor window X enters window W at
    // backbone position i and exits at position j > i, create a bypass
    // edge connecting X's anchors on either side and remove the
    // inter-window edges.
    uint64_t bypassDetourFilter(
        const std::vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys);

private:
    uint64_t nodeCount = 0;
    std::vector<std::map<bool, std::set<OrientedAnchor>>> adjacency;
    std::vector<NodeProperties> nodeProps;
    std::map<std::pair<OrientedAnchor, OrientedAnchor>, EdgeProperties> edgeProperties;
};

} // namespace dinara
