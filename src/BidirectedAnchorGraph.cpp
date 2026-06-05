// BidirectedAnchorGraph.cpp

#include "BidirectedAnchorGraph.hpp"

#include <cmath>
#include <iomanip>

using namespace dinara;
using namespace std;


void BidirectedAnchorGraph::normalizeOrientations(const vector<bool>& swapOrientation)
{
    if(swapOrientation.size() != adjacency.size()) {
        throw runtime_error("swapOrientation size mismatch");
    }

    // Collect all canonical edges with their properties.
    struct EdgeRecord {
        OrientedAnchor from, to;
        EdgeProperties props;
    };
    vector<EdgeRecord> allEdges;
    allEdges.reserve(edgeProperties.size());
    for(const auto& [key, props] : edgeProperties) {
        allEdges.push_back({key.first, key.second, props});
    }

    // Clear everything.
    for(auto& dirMap : adjacency) {
        dirMap.clear();
    }
    edgeProperties.clear();

    // Re-add all edges with adjusted orientations.
    for(const auto& rec : allEdges) {
        OrientedAnchor from = rec.from;
        OrientedAnchor to = rec.to;

        if(swapOrientation[from.first.value()]) from.second = !from.second;
        if(swapOrientation[to.first.value()]) to.second = !to.second;

        addEdge(from, to, rec.props);
        addTraversal(reverseAnchor(to), reverseAnchor(from));
    }

    uint64_t flippedCount = 0;
    for(bool b : swapOrientation) if(b) ++flippedCount;
    cout << "Normalized orientations: flipped " << flippedCount
         << " of " << swapOrientation.size() << " nodes." << endl;
}


void BidirectedAnchorGraph::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    gfa << "H\tVN:Z:1.0\n";

    // Collect links by walking adjacency (preserves actual traversal
    // orientations). Deduplicate by tracking emitted pairs.
    struct LinkRecord {
        OrientedAnchor from, to;
        EdgeProperties props;  // Direction-adjusted copy.
    };
    vector<LinkRecord> links;
    set<pair<OrientedAnchor, OrientedAnchor>> emitted;
    set<uint64_t> activeNodes;

    for(uint64_t nodeIdx = 0; nodeIdx < nodeCount; nodeIdx++) {
        for(bool orient : {true, false}) {
            auto dirIt = adjacency[nodeIdx].find(orient);
            if(dirIt == adjacency[nodeIdx].end()) continue;
            for(const auto& neighbor : dirIt->second) {
                OrientedAnchor from = {BidirectedAnchorId(nodeIdx), orient};
                OrientedAnchor to = neighbor;

                // Skip if RC mirror already emitted.
                auto rcPair = make_pair(reverseAnchor(to), reverseAnchor(from));
                if(emitted.count(rcPair)) continue;

                // Get direction-adjusted properties.
                EdgeProperties props;
                if(!getEdgeProperties(from, to, props)) continue;
                if(!props.useForAssembly) continue;

                emitted.insert(make_pair(from, to));
                links.push_back({from, to, props});
                activeNodes.insert(from.first.value());
                activeNodes.insert(to.first.value());
            }
        }
    }

    // Write segments.
    for(const auto nodeIdx : activeNodes) {
        gfa << "S\t" << nodeIdx << "\t*\tLN:i:1";
        if(nodeIdx < nodeProps.size() && nodeProps[nodeIdx].windowId != UINT32_MAX) {
            gfa << "\twn:i:" << nodeProps[nodeIdx].windowId;
        }
        gfa << "\n";
    }

    // Write links with direction-adjusted properties.
    for(const auto& link : links) {
        gfa << "L\t" << link.from.first.value()
            << "\t" << (link.from.second ? "+" : "-")
            << "\t" << link.to.first.value()
            << "\t" << (link.to.second ? "+" : "-")
            << "\t0M"
            << "\tRC:i:" << link.props.coverage;

        if(link.props.supportingSpanPrev > 0 || link.props.supportingSpanNext > 0) {
            gfa << "\tsp:i:" << link.props.supportingSpanPrev
                << "\tsn:i:" << link.props.supportingSpanNext;
        }
        if(link.props.sharedReadCount > 0) {
            gfa << "\tsr:i:" << link.props.sharedReadCount;
        }

        gfa << "\n";
    }

    cout << "Wrote bidirected GFA to " << fileName
         << ": " << activeNodes.size() << " segments, "
         << links.size() << " links." << endl;
}


void BidirectedAnchorGraph::writeCsv(const string& fileName, uint32_t windowCount) const
{
    ofstream csv(fileName);
    if(!csv) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    csv << "Name,Color\n";

    if(windowCount == 0) return;

    set<uint64_t> activeNodes;
    for(const auto& [key, props] : edgeProperties) {
        if(!props.useForAssembly) continue;
        activeNodes.insert(key.first.first.value());
        activeNodes.insert(key.second.first.value());
    }

    for(const auto nodeIdx : activeNodes) {
        if(nodeIdx >= nodeProps.size()) continue;
        const uint32_t wid = nodeProps[nodeIdx].windowId;
        if(wid == UINT32_MAX) continue;

        const double hue = (360.0 * wid) / windowCount;
        const double s = 0.7;
        const double l = 0.5;
        const double c = (1.0 - std::abs(2.0 * l - 1.0)) * s;
        const double x = c * (1.0 - std::abs(std::fmod(hue / 60.0, 2.0) - 1.0));
        const double m = l - c / 2.0;
        double r1, g1, b1;
        if(hue < 60)       { r1 = c; g1 = x; b1 = 0; }
        else if(hue < 120) { r1 = x; g1 = c; b1 = 0; }
        else if(hue < 180) { r1 = 0; g1 = c; b1 = x; }
        else if(hue < 240) { r1 = 0; g1 = x; b1 = c; }
        else if(hue < 300) { r1 = x; g1 = 0; b1 = c; }
        else               { r1 = c; g1 = 0; b1 = x; }
        const int r = int((r1 + m) * 255);
        const int g = int((g1 + m) * 255);
        const int b = int((b1 + m) * 255);

        csv << nodeIdx << ","
            << "#" << hex << setfill('0')
            << setw(2) << r << setw(2) << g << setw(2) << b
            << dec << "\n";
    }
}
