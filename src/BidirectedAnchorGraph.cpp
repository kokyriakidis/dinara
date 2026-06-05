// BidirectedAnchorGraph.cpp

#include "BidirectedAnchorGraph.hpp"

#include <cmath>
#include <iomanip>

using namespace dinara;
using namespace std;


void BidirectedAnchorGraph::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    gfa << "H\tVN:Z:1.0\n";

    // Collect active nodes (nodes with at least one edge).
    set<uint64_t> activeNodes;
    for(const auto& [key, props] : edgeProperties) {
        activeNodes.insert(key.first.first);
        activeNodes.insert(key.second.first);
    }

    // Write segments.
    for(const auto nodeId : activeNodes) {
        gfa << "S\t" << nodeId << "\t*\tLN:i:1";
        if(nodeId < nodeProps.size() && nodeProps[nodeId].windowId != UINT32_MAX) {
            gfa << "\twn:i:" << nodeProps[nodeId].windowId;
        }
        gfa << "\n";
    }

    // Write links. edgeProperties is already keyed by canonical form,
    // so each link is emitted exactly once.
    for(const auto& [key, props] : edgeProperties) {
        if(!props.useForAssembly) continue;

        const auto& from = key.first;
        const auto& to = key.second;

        gfa << "L\t" << from.first << "\t" << (from.second ? "+" : "-")
            << "\t" << to.first << "\t" << (to.second ? "+" : "-")
            << "\t0M"
            << "\tRC:i:" << props.coverage;

        if(props.supportingSpanPrev > 0 || props.supportingSpanNext > 0) {
            gfa << "\tsp:i:" << props.supportingSpanPrev
                << "\tsn:i:" << props.supportingSpanNext;
        }
        if(props.sharedReadCount > 0) {
            gfa << "\tsr:i:" << props.sharedReadCount;
        }

        gfa << "\n";
    }

    cout << "Wrote bidirected GFA to " << fileName
         << ": " << activeNodes.size() << " segments, "
         << numEdges() << " links." << endl;
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
        activeNodes.insert(key.first.first);
        activeNodes.insert(key.second.first);
    }

    for(const auto nodeId : activeNodes) {
        if(nodeId >= nodeProps.size()) continue;
        const uint32_t wid = nodeProps[nodeId].windowId;
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

        csv << nodeId << ","
            << "#" << hex << setfill('0')
            << setw(2) << r << setw(2) << g << setw(2) << b
            << dec << "\n";
    }
}
