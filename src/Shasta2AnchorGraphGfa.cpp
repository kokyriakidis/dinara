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



BidirectedAnchorGraph Shasta2AnchorGraph::toBidirected(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Journeys& journeys) const
{
    const uint64_t anchorCount = anchorToWindow.size();
    const uint64_t bidirNodeCount = (anchorCount + 1) / 2;

    BidirectedAnchorGraph bg;
    bg.resize(bidirNodeCount);

    auto normalizeWindow = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // Set node window and backbone read assignments using the forward anchor (even).
    for(uint64_t anchorId = 0; anchorId < anchorCount; anchorId += 2) {
        const auto bid = BidirectedAnchorId::fromDirected(Shasta2AnchorId(anchorId));
        if(anchorId < anchorToWindow.size()) {
            const uint32_t wid = anchorToWindow[anchorId];
            if(wid != noWindow) {
                const uint32_t normWid = normalizeWindow(wid);
                bg.setNodeWindow(bid, normWid);
                if(normWid < anchorWindows.size()) {
                    bg.setNodeBackboneRead(bid,
                        uint32_t(anchorWindows[normWid].backboneOrientedReadId.getReadId()));
                }
            }
        }
    }

    // Build the set of valid consecutive backbone anchor pairs.
    // For intra-window edges, only these pairs are converted.
    set<pair<uint64_t, uint64_t>> backbonePairs;
    for(const auto& window : anchorWindows) {
        const auto& positions = window.filteredBackbonePositions;
        if(positions.size() < 2) continue;
        const auto journey = journeys[window.backboneOrientedReadId];
        for(uint32_t i = 0; i + 1 < uint32_t(positions.size()); i++) {
            const uint64_t a = uint64_t(journey[positions[i]]);
            const uint64_t b = uint64_t(journey[positions[i + 1]]);
            backbonePairs.insert({a, b});
            // Also insert RC mirror pair.
            backbonePairs.insert({b ^ 1ULL, a ^ 1ULL});
        }
    }

    // Convert edges. Each directed edge maps to a bidirected traversal.
    // For intra-window edges, only consecutive backbone pairs are kept.
    // Inter-window edges are always kept.
    uint64_t skippedNonBackbone = 0;
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        const auto& edge = (*this)[e];
        if(!edge.useForAssembly) continue;

        const uint64_t srcRaw = uint64_t(source(e, *this));
        const uint64_t dstRaw = uint64_t(target(e, *this));

        // Check if this is an intra-window edge.
        const uint32_t srcWin = (srcRaw < anchorToWindow.size()) ? anchorToWindow[srcRaw] : noWindow;
        const uint32_t dstWin = (dstRaw < anchorToWindow.size()) ? anchorToWindow[dstRaw] : noWindow;
        const bool isInterWindow = (srcWin == noWindow || dstWin == noWindow ||
                                    normalizeWindow(srcWin) != normalizeWindow(dstWin));

        // For intra-window edges, only keep consecutive backbone pairs.
        if(!isInterWindow) {
            if(backbonePairs.find({srcRaw, dstRaw}) == backbonePairs.end()) {
                ++skippedNonBackbone;
                continue;
            }
        }

        const Shasta2AnchorId srcId = Shasta2AnchorId(srcRaw);
        const Shasta2AnchorId dstId = Shasta2AnchorId(dstRaw);

        const OrientedAnchor from = toOrientedAnchor(srcId);
        const OrientedAnchor to   = toOrientedAnchor(dstId);

        BidirectedAnchorGraph::EdgeProperties props;
        props.coverage = edge.coverage();
        props.offset = edge.offset;
        props.supportingSpanPrev = uint32_t(edge.supportingSpanPrev);
        props.supportingSpanNext = uint32_t(edge.supportingSpanNext);
        props.sharedReadCount = uint32_t(edge.sharedReadCount);
        props.isInterWindow = isInterWindow;
        props.useForAssembly = true;

        if(isCanonicalAnchorDirection(from, to) ||
           bg.getEdgePropertiesCanonical(from, to) == nullptr) {
            bg.addEdge(from, to, props);
        } else {
            bg.addTraversal(from, to);
        }
    }

    if(skippedNonBackbone > 0) {
        cout << "toBidirected: skipped " << skippedNonBackbone
             << " non-backbone intra-window edges." << endl;
    }

    cout << "Converted to bidirected graph: "
         << bg.numNodes() << " nodes, "
         << bg.numEdges() << " edges." << endl;

    // Normalize orientations so backbone chains appear as +/+ links.
    // Flip nodes whose backbone anchor is odd (strand-1 reads).
    vector<bool> swapOrientation(bidirNodeCount, false);
    for(uint32_t windowId = 0; windowId < windowCount; windowId++) {
        const AnchorWindow& window = anchorWindows[windowId];
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = journeys[backboneOid];

        const auto& positions = window.filteredBackbonePositions.empty()
            ? [&]() -> const vector<uint32_t>& {
                static thread_local vector<uint32_t> allPositions;
                allPositions.clear();
                for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                    allPositions.push_back(pos);
                }
                return allPositions;
            }()
            : window.filteredBackbonePositions;

        for(const uint32_t pos : positions) {
            const Shasta2AnchorId aid = backboneJourney[pos];
            const auto bid = BidirectedAnchorId::fromDirected(aid);
            if(bid.value() >= bidirNodeCount) continue;
            if(uint64_t(aid) & 1) {
                swapOrientation[bid.value()] = true;
            }
        }
    }
    bg.normalizeOrientations(swapOrientation);

    return bg;
}
