#pragma once

// BidirectedAnchorGraph.hpp
//
// Verkko/MBG-style bidirected anchor graph. Each node represents an
// anchor pair (anchorId / 2), collapsing forward and RC into one node.
// Each node has two ends (forward = true, RC = false). Edges connect
// oriented node ends. Adding an edge automatically adds its RC mirror,
// so both traversal directions are stored, but they represent one
// biological link.

#include "cstdint.hpp"

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace dinara {

class BidirectedAnchorGraph {
public:
    // An oriented node: (nodeId, forward).
    // forward=true means the node is traversed in the forward direction.
    // forward=false means the node is traversed in the RC direction.
    using OrientedNode = std::pair<uint64_t, bool>;

    // Reverse an oriented node (flip orientation).
    static OrientedNode reverse(OrientedNode node) {
        return {node.first, !node.second};
    }

    // Canonical form of a link. Ensures a consistent ordering so that
    // an edge and its RC mirror map to the same canonical key.
    // RC mirror of (A, fwA) -> (B, fwB) is (B, !fwB) -> (A, !fwA).
    static std::pair<OrientedNode, OrientedNode> canon(
        OrientedNode from, OrientedNode to)
    {
        auto fwd = std::make_pair(from, to);
        auto rev = std::make_pair(reverse(to), reverse(from));
        return std::min(fwd, rev);
    }

    // Edge properties stored once per canonical link.
    struct EdgeProperties {
        uint64_t coverage = 0;
        uint32_t supportingSpanPrev = 0;
        uint32_t supportingSpanNext = 0;
        uint32_t sharedReadCount = 0;
        bool isInterWindow = false;
        bool useForAssembly = true;
    };

    // Node properties.
    struct NodeProperties {
        uint32_t windowId = UINT32_MAX;  // Normalized window ID (noWindow if unassigned).
    };

    // Number of nodes.
    uint64_t numNodes() const { return nodeCount; }

    // Number of canonical edges (each link counted once).
    uint64_t numEdges() const { return edgeProperties.size(); }

    // Resize to hold the given number of nodes.
    void resize(uint64_t n) {
        nodeCount = n;
        adjacency.resize(n);
        nodeProps.resize(n);
    }

    // Set node properties.
    void setNodeWindow(uint64_t nodeId, uint32_t windowId) {
        nodeProps[nodeId].windowId = windowId;
    }

    // Add a bidirected edge. Stores both traversal directions and
    // edge properties under the canonical key.
    void addEdge(OrientedNode from, OrientedNode to,
                 const EdgeProperties& props)
    {
        // Store both traversal directions.
        adjacency[from.first][from.second].insert(to);
        adjacency[to.first][!to.second].insert(reverse(from));

        // Store properties under canonical key.
        auto key = canon(from, to);
        edgeProperties[key] = props;
    }

    // Check if an edge exists.
    bool hasEdge(OrientedNode from, OrientedNode to) const {
        if(from.first >= adjacency.size()) return false;
        const auto& dirMap = adjacency[from.first];
        auto it = dirMap.find(from.second);
        if(it == dirMap.end()) return false;
        return it->second.count(to) > 0;
    }

    // Get neighbors of an oriented node.
    std::vector<OrientedNode> getNeighbors(OrientedNode node) const {
        std::vector<OrientedNode> result;
        if(node.first >= adjacency.size()) return result;
        const auto& dirMap = adjacency[node.first];
        auto it = dirMap.find(node.second);
        if(it == dirMap.end()) return result;
        result.assign(it->second.begin(), it->second.end());
        return result;
    }

    // Get edge properties (by canonical key).
    const EdgeProperties* getEdgeProperties(OrientedNode from, OrientedNode to) const {
        auto key = canon(from, to);
        auto it = edgeProperties.find(key);
        if(it == edgeProperties.end()) return nullptr;
        return &it->second;
    }

    EdgeProperties* getEdgeProperties(OrientedNode from, OrientedNode to) {
        auto key = canon(from, to);
        auto it = edgeProperties.find(key);
        if(it == edgeProperties.end()) return nullptr;
        return &it->second;
    }

    // Remove an edge (both directions + properties).
    void removeEdge(OrientedNode from, OrientedNode to) {
        if(from.first < adjacency.size()) {
            adjacency[from.first][from.second].erase(to);
        }
        auto revFrom = reverse(to);
        auto revTo = reverse(from);
        if(revFrom.first < adjacency.size()) {
            adjacency[revFrom.first][revFrom.second].erase(revTo);
        }
        edgeProperties.erase(canon(from, to));
    }

    // In-degree + out-degree for a node (counting both orientations).
    uint64_t degree(uint64_t nodeId) const {
        uint64_t d = 0;
        auto it = adjacency[nodeId].find(true);
        if(it != adjacency[nodeId].end()) d += it->second.size();
        it = adjacency[nodeId].find(false);
        if(it != adjacency[nodeId].end()) d += it->second.size();
        return d;
    }

    // Write GFA with proper bidirected orientations.
    void writeGfa(const std::string& fileName) const;

    // Write CSV with node colors by window.
    void writeCsv(const std::string& fileName, uint32_t windowCount) const;

private:
    uint64_t nodeCount = 0;

    // Adjacency: nodeId -> {orientation -> set of oriented neighbors}.
    // adjacency[nodeId][true] = neighbors when traversing node forward.
    // adjacency[nodeId][false] = neighbors when traversing node in RC.
    std::vector<std::map<bool, std::set<OrientedNode>>> adjacency;

    // Node properties indexed by nodeId.
    std::vector<NodeProperties> nodeProps;

    // Edge properties indexed by canonical (from, to) pair.
    std::map<std::pair<OrientedNode, OrientedNode>, EdgeProperties> edgeProperties;
};

} // namespace dinara
