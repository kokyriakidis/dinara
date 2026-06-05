// BidirectedAnchorGraph.cpp

#include "BidirectedAnchorGraph.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

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


// Unitigify: collapse linear chains of nodes into unitigs.
//
// In a bidirected graph, a node v is internal to a linear chain when:
// - sideDegree(v+) == 1 and sideDegree(v-) == 1
// - The single neighbor w reached from v+ has sideDegree == 1 on the
//   entry side (!w_orient)
// - Similarly for the neighbor reached from v-
//
// The algorithm walks from non-internal nodes along chains of internal
// nodes, collecting oriented anchors.
vector<BidirectedAnchorGraph::Unitig> BidirectedAnchorGraph::unitigify() const
{
    auto isInternal = [&](BidirectedAnchorId id) -> bool {
        if(sideDegree({id, true}) != 1) return false;
        if(sideDegree({id, false}) != 1) return false;

        auto neighborsPlus = getNeighbors({id, true});
        if(neighborsPlus.size() != 1) return false;
        const OrientedAnchor& wPlus = neighborsPlus[0];
        if(sideDegree({wPlus.first, !wPlus.second}) != 1) return false;

        auto neighborsMinus = getNeighbors({id, false});
        if(neighborsMinus.size() != 1) return false;
        const OrientedAnchor& wMinus = neighborsMinus[0];
        if(sideDegree({wMinus.first, !wMinus.second}) != 1) return false;

        return true;
    };

    // Track which nodes have been consumed as internal chain members.
    vector<bool> visited(nodeCount, false);

    vector<Unitig> unitigs;

    for(uint64_t nodeIdx = 0; nodeIdx < nodeCount; nodeIdx++) {
        if(visited[nodeIdx]) continue;

        BidirectedAnchorId id(nodeIdx);
        if(degree(id) == 0) continue;
        if(isInternal(id)) continue;  // Will be picked up by a chain walk.

        // Branch/endpoint node. Start a chain from each side that has
        // exactly one neighbor.
        for(bool startOrient : {true, false}) {
            if(sideDegree({id, startOrient}) != 1) continue;

            auto startNeighbors = getNeighbors({id, startOrient});
            if(startNeighbors.size() != 1) continue;

            const OrientedAnchor firstNeighbor = startNeighbors[0];

            // Don't start if the neighbor was already visited
            // (chain was already built from the other end).
            if(visited[firstNeighbor.first.value()]) continue;

            Unitig unitig;
            unitig.chain.push_back({id, startOrient});

            OrientedAnchor current = firstNeighbor;
            while(true) {
                unitig.chain.push_back(current);

                if(!isInternal(current.first)) break;

                // Only mark internal nodes as visited. Branch/endpoint
                // nodes at chain ends must remain unvisited so they can
                // participate as endpoints of other chains.
                visited[current.first.value()] = true;

                // Exit from current on the traversal side (current.second).
                auto nextNeighbors = getNeighbors({current.first, current.second});
                if(nextNeighbors.size() != 1) break;
                current = nextNeighbors[0];

                if(current.first == id) break;  // Circular safety.
            }

            // Compute window sequence (consecutive duplicates removed).
            static constexpr uint32_t noWindow = UINT32_MAX;
            for(const auto& oa : unitig.chain) {
                auto idx = oa.first.value();
                if(idx < nodeProps.size()) {
                    uint32_t wid = nodeProps[idx].windowId;
                    if(wid != noWindow) {
                        if(unitig.windowSequence.empty() ||
                           unitig.windowSequence.back() != wid) {
                            unitig.windowSequence.push_back(wid);
                        }
                    }
                }
            }

            // Average edge coverage and total offset along the chain.
            if(unitig.chain.size() >= 2) {
                double totalCov = 0.;
                uint64_t edgeCount = 0;
                uint64_t totalOff = 0;
                for(uint64_t i = 0; i + 1 < unitig.chain.size(); i++) {
                    EdgeProperties props;
                    if(getEdgeProperties(unitig.chain[i], unitig.chain[i+1], props)) {
                        totalCov += double(props.coverage);
                        totalOff += props.offset;
                        edgeCount++;
                    }
                }
                if(edgeCount > 0) {
                    unitig.averageCoverage = totalCov / double(edgeCount);
                }
                unitig.totalOffset = totalOff;
            }

            unitigs.push_back(std::move(unitig));
        }
    }

    cout << "Unitigified: " << unitigs.size() << " unitigs from "
         << nodeCount << " nodes." << endl;

    uint64_t totalAnchors = 0;
    uint64_t maxAnchors = 0;
    for(const auto& u : unitigs) {
        totalAnchors += u.anchorCount();
        maxAnchors = max(maxAnchors, u.anchorCount());
    }
    cout << "  Total anchors in unitigs: " << totalAnchors << endl;
    cout << "  Max anchors in a unitig: " << maxAnchors << endl;
    if(!unitigs.empty()) {
        cout << "  Average anchors per unitig: "
             << double(totalAnchors) / double(unitigs.size()) << endl;
    }

    return unitigs;
}


// HSL to hex color string.
static string hslToHex(double hue, double s, double l)
{
    const double c = (1.0 - abs(2.0 * l - 1.0)) * s;
    const double x = c * (1.0 - abs(fmod(hue / 60.0, 2.0) - 1.0));
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

    ostringstream oss;
    oss << "#" << hex << setfill('0')
        << setw(2) << r << setw(2) << g << setw(2) << b;
    return oss.str();
}


// Write unitig GFA: each unitig is a segment, links connect unitig endpoints.
void BidirectedAnchorGraph::writeUnitigGfa(
    const string& fileName,
    const vector<Unitig>& unitigs,
    uint32_t windowCount) const
{
    ofstream gfa(fileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    gfa << "H\tVN:Z:1.0\n";

    // Write segments.
    for(uint64_t i = 0; i < unitigs.size(); i++) {
        const auto& u = unitigs[i];
        gfa << "S\t" << i << "\t*"
            << "\tLN:i:" << u.totalOffset;

        if(u.averageCoverage > 0.) {
            gfa << "\tRC:i:" << uint64_t(u.averageCoverage * double(u.totalOffset));
            gfa << "\tdp:f:" << fixed << setprecision(1) << u.averageCoverage;
        }

        if(!u.windowSequence.empty()) {
            gfa << "\tws:Z:";
            for(uint64_t j = 0; j < u.windowSequence.size(); j++) {
                if(j > 0) gfa << ",";
                gfa << u.windowSequence[j];
            }
        }

        gfa << "\tan:i:" << u.anchorCount();

        if(!u.chain.empty()) {
            gfa << "\tfa:Z:" << u.chain.front().first.value()
                << (u.chain.front().second ? "+" : "-");
            gfa << "\tla:Z:" << u.chain.back().first.value()
                << (u.chain.back().second ? "+" : "-");
        }

        gfa << "\n";
    }

    // Map oriented anchors to unitig entry points.
    // If neighbor == chain[0] of unitig j, we enter j from its front (j+).
    // If neighbor == reverseAnchor(chain.back()) of unitig j, we enter j
    // from its back in reverse (j-).
    map<OrientedAnchor, uint64_t> frontEntry, backEntry;
    for(uint64_t i = 0; i < unitigs.size(); i++) {
        frontEntry[unitigs[i].chain.front()] = i;
        backEntry[reverseAnchor(unitigs[i].chain.back())] = i;
    }

    // Deduplicate links (each link and its RC mirror are the same).
    using LinkKey = pair<pair<uint64_t,bool>, pair<uint64_t,bool>>;
    set<LinkKey> emittedLinks;

    auto emitLink = [&](uint64_t fromUnitig, bool fromOrient,
                        uint64_t toUnitig, bool toOrient,
                        const EdgeProperties& props) {
        LinkKey link = {{fromUnitig, fromOrient}, {toUnitig, toOrient}};
        LinkKey rcLink = {{toUnitig, !toOrient}, {fromUnitig, !fromOrient}};
        if(emittedLinks.count(rcLink)) return;
        emittedLinks.insert(link);
        gfa << "L\t" << fromUnitig << "\t" << (fromOrient ? "+" : "-")
            << "\t" << toUnitig << "\t" << (toOrient ? "+" : "-")
            << "\t0M"
            << "\tRC:i:" << props.coverage << "\n";
    };

    for(uint64_t i = 0; i < unitigs.size(); i++) {
        const auto& u = unitigs[i];

        // Links from the back of unitig i (i+).
        {
            const OrientedAnchor exitAnchor = u.chain.back();
            auto neighbors = getNeighbors({exitAnchor.first, exitAnchor.second});
            for(const auto& neighbor : neighbors) {
                EdgeProperties props;
                if(!getEdgeProperties(exitAnchor, neighbor, props)) continue;
                if(!props.useForAssembly) continue;

                auto itFront = frontEntry.find(neighbor);
                if(itFront != frontEntry.end()) {
                    emitLink(i, true, itFront->second, true, props);
                    continue;
                }
                auto itBack = backEntry.find(neighbor);
                if(itBack != backEntry.end()) {
                    emitLink(i, true, itBack->second, false, props);
                }
            }
        }

        // Links from the front of unitig i (i-).
        {
            const OrientedAnchor exitAnchor = reverseAnchor(u.chain.front());
            auto neighbors = getNeighbors({exitAnchor.first, exitAnchor.second});
            for(const auto& neighbor : neighbors) {
                EdgeProperties props;
                if(!getEdgeProperties(exitAnchor, neighbor, props)) continue;
                if(!props.useForAssembly) continue;

                auto itFront = frontEntry.find(neighbor);
                if(itFront != frontEntry.end()) {
                    emitLink(i, false, itFront->second, true, props);
                    continue;
                }
                auto itBack = backEntry.find(neighbor);
                if(itBack != backEntry.end()) {
                    emitLink(i, false, itBack->second, false, props);
                }
            }
        }
    }

    cout << "Wrote unitig GFA to " << fileName
         << ": " << unitigs.size() << " segments, "
         << emittedLinks.size() << " links." << endl;
}


void BidirectedAnchorGraph::writeUnitigCsv(
    const string& fileName,
    const vector<Unitig>& unitigs,
    uint32_t windowCount) const
{
    ofstream csv(fileName);
    if(!csv) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    csv << "Name,Color\n";

    if(windowCount == 0) return;

    for(uint64_t i = 0; i < unitigs.size(); i++) {
        const auto& u = unitigs[i];
        if(u.windowSequence.empty()) continue;

        const uint32_t wid = u.windowSequence[0];
        const double hue = (360.0 * wid) / windowCount;
        csv << i << "," << hslToHex(hue, 0.7, 0.5) << "\n";
    }
}
