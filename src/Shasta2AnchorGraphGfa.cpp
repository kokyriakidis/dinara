// Shasta2.
#include "Shasta2AnchorGraph.hpp"

#include <boost/graph/iteration_macros.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>

using namespace dinara;
using namespace std;

void Shasta2AnchorGraph::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    gfa << "H\tVN:Z:1.0\n";

    // Collect vertices that have at least one active edge.
    std::set<vertex_descriptor> activeVertices;
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        if(!(*this)[e].useForAssembly) continue;
        activeVertices.insert(source(e, *this));
        activeVertices.insert(target(e, *this));
    }

    // Write active vertices only.
    for(const auto v : activeVertices) {
        gfa << "S\t" << v << "\t*\tLN:i:1\n";
    }

    // Write edges. All edges are forward-to-forward.
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        const auto& edge = (*this)[e];
        if(!edge.useForAssembly) continue;
        gfa << "L\t" << source(e, *this) << "\t+\t"
            << target(e, *this) << "\t+\t0M"
            << "\tRC:i:" << edge.coverage() << "\n";
    }
}


void Shasta2AnchorGraph::writeCsv(const string& fileName) const
{
    ofstream csv(fileName);
    if(!csv) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    csv << "Name,Color\n";

    if(windowCount == 0 || anchorToWindow.empty()) return;

    // Build window adjacency graph from active inter-window edges,
    // then greedy-color so adjacent windows get maximally different colors.
    auto normalize = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    std::map<uint32_t, std::set<uint32_t>> windowNeighbors;
    const uint64_t anchorCount = anchorToWindow.size();
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        if(!(*this)[e].useForAssembly) continue;
        const uint64_t srcVal = uint64_t(source(e, *this));
        const uint64_t dstVal = uint64_t(target(e, *this));
        if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
        const uint32_t srcWin = normalize(anchorToWindow[srcVal]);
        const uint32_t dstWin = normalize(anchorToWindow[dstVal]);
        if(srcWin == noWindow || dstWin == noWindow) continue;
        if(srcWin == dstWin) continue;
        windowNeighbors[srcWin].insert(dstWin);
        windowNeighbors[dstWin].insert(srcWin);
    }

    // Use a fixed palette of 20 colors. Assign each window the color
    // index that is least used among its neighbors, breaking ties by
    // picking the lowest index. This maximizes local contrast.
    const uint32_t numColors = 20;
    std::map<uint32_t, uint32_t> windowColorIndex;
    for(uint32_t w = 0; w < windowCount; w++) {
        // Count how many neighbors use each color.
        std::vector<uint32_t> neighborColorCount(numColors, 0);
        if(windowNeighbors.count(w)) {
            for(const uint32_t nb : windowNeighbors[w]) {
                if(windowColorIndex.count(nb)) {
                    neighborColorCount[windowColorIndex[nb]]++;
                }
            }
        }
        // Pick the color with the lowest neighbor count.
        uint32_t bestColor = 0;
        for(uint32_t c = 1; c < numColors; c++) {
            if(neighborColorCount[c] < neighborColorCount[bestColor]) {
                bestColor = c;
            }
        }
        windowColorIndex[w] = bestColor;
    }

    // HSL-to-RGB helper.
    auto hslToRgb = [](double h, double s, double l,
                       int& rOut, int& gOut, int& bOut) {
        const double c = (1.0 - std::abs(2.0 * l - 1.0)) * s;
        const double x = c * (1.0 - std::abs(std::fmod(h / 60.0, 2.0) - 1.0));
        const double m = l - c / 2.0;
        double r1, g1, b1;
        if(h < 60)       { r1 = c; g1 = x; b1 = 0; }
        else if(h < 120) { r1 = x; g1 = c; b1 = 0; }
        else if(h < 180) { r1 = 0; g1 = c; b1 = x; }
        else if(h < 240) { r1 = 0; g1 = x; b1 = c; }
        else if(h < 300) { r1 = x; g1 = 0; b1 = c; }
        else              { r1 = c; g1 = 0; b1 = x; }
        rOut = int((r1 + m) * 255);
        gOut = int((g1 + m) * 255);
        bOut = int((b1 + m) * 255);
    };

    // Generate a palette: spread colors evenly in hue, vary S/L per band.
    const double goldenRatio = 0.618033988749895;
    std::vector<std::tuple<int,int,int>> palette(numColors);
    for(uint32_t i = 0; i < numColors; i++) {
        const double hue = std::fmod(i * goldenRatio * 360.0, 360.0);
        const double sat = 0.6 + 0.3 * ((i % 3) / 2.0);
        const double lit = 0.35 + 0.2 * (((i / 3) % 3) / 2.0);
        int r, g, b;
        hslToRgb(hue, sat, lit, r, g, b);
        palette[i] = {r, g, b};
    }

    BGL_FORALL_VERTICES(v, *this, Shasta2AnchorGraph) {
        if(in_degree(v, *this) == 0 && out_degree(v, *this) == 0) continue;
        if(uint64_t(v) >= anchorToWindow.size()) continue;

        const uint32_t wid = anchorToWindow[uint64_t(v)];
        if(wid == noWindow) continue;

        const uint32_t fwdWid = normalize(wid);
        const uint32_t colorIdx = windowColorIndex.count(fwdWid)
            ? windowColorIndex[fwdWid] % numColors : 0;
        const auto& [r, g, b] = palette[colorIdx];

        csv << v << ","
            << "#" << hex << setfill('0')
            << setw(2) << r << setw(2) << g << setw(2) << b
            << dec << "\n";
    }
}



void Shasta2AnchorGraph::writeBubbleFinderGraph(const string& fileName, bool useForAssemblyOnly) const
{
    ofstream graphFile(fileName);
    if(!graphFile) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }

    std::set<vertex_descriptor> verticesInOutput;
    uint64_t edgeCount = 0;
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        const auto& edge = (*this)[e];
        if(useForAssemblyOnly && !edge.useForAssembly) {
            continue;
        }
        verticesInOutput.insert(source(e, *this));
        verticesInOutput.insert(target(e, *this));
        ++edgeCount;
    }

    graphFile << verticesInOutput.size() << " " << edgeCount << "\n";
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        const auto& edge = (*this)[e];
        if(useForAssemblyOnly && !edge.useForAssembly) {
            continue;
        }
        graphFile << source(e, *this) << " " << target(e, *this) << "\n";
    }
}
