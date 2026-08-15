// Shasta2.
#include "Shasta2AnchorGraph.hpp"

#include <boost/graph/iteration_macros.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>

using namespace dinara;
using namespace std;

void Shasta2AnchorGraph::writeGfa(const string& fileName,
                                  const vector<AnchorWindow>* anchorWindows) const
{
    ofstream gfa(fileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    gfa << "H\tVN:Z:1.0\n";

    const uint64_t anchorCount = anchorToWindow.size();

    auto normalize = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // Collect vertices that have at least one active edge.
    std::set<vertex_descriptor> activeVertices;
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        if(!(*this)[e].useForAssembly) continue;
        activeVertices.insert(source(e, *this));
        activeVertices.insert(target(e, *this));
    }

    // Write active vertices with window ID and strand tags.
    for(const auto v : activeVertices) {
        gfa << "S\t" << v << "\t*\tLN:i:1";
        if(uint64_t(v) < anchorCount) {
            const uint32_t wid = anchorToWindow[uint64_t(v)];
            if(wid != noWindow) {
                const uint32_t fwdWid = normalize(wid);
                const bool isRc = (wid >= windowCount);
                gfa << "\twn:i:" << fwdWid
                    << "\tws:Z:" << (isRc ? "rc" : "fw");
            }
        }
        gfa << "\n";
    }

    // Write edges. All edges are forward-to-forward.
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        const auto& edge = (*this)[e];
        if(!edge.useForAssembly) continue;

        const auto src = source(e, *this);
        const auto dst = target(e, *this);

        gfa << "L\t" << src << "\t+\t"
            << dst << "\t+\t0M"
            << "\tRC:i:" << edge.coverage();

        // Inter-window edge attributes.
        if(edge.supportingSpanPrev > 0 || edge.supportingSpanNext > 0) {
            gfa << "\tsp:i:" << edge.supportingSpanPrev
                << "\tsn:i:" << edge.supportingSpanNext;
        }
        if(edge.sharedReadCount > 0) {
            gfa << "\tsr:i:" << edge.sharedReadCount;
        }

        gfa << "\n";
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

    auto normalize = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // Collect vertices that have at least one active edge.
    std::set<vertex_descriptor> activeVertices;
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        if(!(*this)[e].useForAssembly) continue;
        activeVertices.insert(source(e, *this));
        activeVertices.insert(target(e, *this));
    }

    // One unique color per window. Golden angle spacing for distinct colors.
    for(const auto v : activeVertices) {
        if(uint64_t(v) >= anchorToWindow.size()) continue;

        const uint32_t wid = anchorToWindow[uint64_t(v)];
        if(wid == noWindow) continue;

        const uint32_t fwdWid = normalize(wid);

        const double hue = fmod(double(fwdWid) * 137.508, 360.0);
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

        csv << v << ","
            << "#" << hex << setfill('0')
            << setw(2) << r << setw(2) << g << setw(2) << b
            << dec << "\n";
    }
}

