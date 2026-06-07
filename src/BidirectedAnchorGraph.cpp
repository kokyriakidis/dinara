// BidirectedAnchorGraph.cpp

#include "BidirectedAnchorGraph.hpp"
#include "AnchorWindows.hpp"
#include "Shasta2Journeys.hpp"

#include <boost/graph/adjacency_list.hpp>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

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

    // Store swap flags for post-normalization anchor lookups.
    wasSwapped = swapOrientation;

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

        const double hue = fmod(double(wid) * 137.508, 360.0);
        csv << nodeIdx << "," << hslToHex(hue, 0.7, 0.5) << "\n";
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
vector<BidirectedAnchorGraph::Unitig> BidirectedAnchorGraph::unitigify(bool quiet) const
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

    // Deduplicate: a unitig and its RC mirror are the same.
    // For each unitig, compute its canonical form (the lexicographically
    // smaller of the chain and its RC mirror). Skip duplicates.
    {
        // Canonicalize a chain: return the lexicographically smaller of
        // the chain and its RC mirror.
        auto canonicalizeChain = [](const vector<OrientedAnchor>& chain)
            -> vector<OrientedAnchor>
        {
            vector<OrientedAnchor> rc(chain.size());
            for(uint64_t i = 0; i < chain.size(); i++) {
                rc[chain.size() - 1 - i] = reverseAnchor(chain[i]);
            }
            return (chain <= rc) ? chain : rc;
        };

        set<vector<OrientedAnchor>> seen;
        vector<Unitig> deduped;

        for(auto& u : unitigs) {
            auto canon = canonicalizeChain(u.chain);
            if(!seen.insert(canon).second) continue;

            // If the canonical form differs from the original, flip.
            if(canon != u.chain) {
                u.chain = std::move(canon);
                // Recompute window sequence for the reversed chain.
                u.windowSequence.clear();
                static constexpr uint32_t noWindow = UINT32_MAX;
                for(const auto& oa : u.chain) {
                    auto idx = oa.first.value();
                    if(idx < nodeProps.size()) {
                        uint32_t wid = nodeProps[idx].windowId;
                        if(wid != noWindow) {
                            if(u.windowSequence.empty() ||
                               u.windowSequence.back() != wid) {
                                u.windowSequence.push_back(wid);
                            }
                        }
                    }
                }
            }

            deduped.push_back(std::move(u));
        }

        unitigs = std::move(deduped);
    }

    if(!quiet) {
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
    }

    return unitigs;
}


// Remove tip unitigs from the bidirected graph.
// Follows the same pattern as removeTipWindows: build unitig-level
// adjacency, walk linear chains from dead-ends, process shortest-first,
// re-walk after each removal.
uint64_t BidirectedAnchorGraph::removeTips(uint64_t maxTipLength)
{
    // Iterate: unitigify, collect all tip chains sorted shortest-first,
    // remove non-overlapping tips in batch, repeat until none remain.

    struct UnitigEntry {
        uint64_t unitigIndex;
        bool orientation;  // true = +, false = -
    };

    uint64_t totalChainsRemoved = 0;
    uint64_t totalUnitigsRemoved = 0;

    while(true) {
        auto unitigs = unitigify(/*quiet=*/true);

        // Build entry map.
        std::map<OrientedAnchor, UnitigEntry> entryMap;
        for(uint64_t i = 0; i < unitigs.size(); i++) {
            const auto& front = unitigs[i].chain.front();
            const auto& back = unitigs[i].chain.back();
            entryMap[front]               = {i, true};
            entryMap[reverseAnchor(back)]  = {i, false};
            entryMap[reverseAnchor(front)] = {i, false};
            entryMap[back]                 = {i, true};
        }

        auto getUnitigLinks = [&](uint64_t i, bool backSide) -> std::vector<UnitigEntry>
        {
            const auto& u = unitigs[i];
            OrientedAnchor exitAnchor = backSide
                ? u.chain.back()
                : reverseAnchor(u.chain.front());
            auto neighbors = getNeighbors({exitAnchor.first, exitAnchor.second});
            std::vector<UnitigEntry> result;
            for(const auto& nbr : neighbors) {
                auto it = entryMap.find(nbr);
                if(it != entryMap.end() && it->second.unitigIndex != i) {
                    result.push_back(it->second);
                }
            }
            return result;
        };

        auto sideNeighborCount = [&](uint64_t i, bool backSide) -> uint64_t {
            auto links = getUnitigLinks(i, backSide);
            std::set<uint64_t> distinct;
            for(const auto& e : links) distinct.insert(e.unitigIndex);
            return distinct.size();
        };

        // Collect all tip chains, sorted shortest-first.
        // Each entry: (totalLen, chain of unitig indices).
        std::vector<std::pair<uint64_t, std::vector<uint64_t>>> candidates;

        for(uint64_t i = 0; i < unitigs.size(); i++) {
            uint64_t frontCount = sideNeighborCount(i, false);
            uint64_t backCount = sideNeighborCount(i, true);
            if((frontCount == 0) == (backCount == 0)) continue;

            bool exitBack = (backCount > 0);
            std::vector<uint64_t> chain;
            chain.push_back(i);
            uint64_t current = i;
            bool currentExitBack = exitBack;
            uint64_t totalLen = unitigs[i].totalOffset;

            while(totalLen <= maxTipLength) {
                auto links = getUnitigLinks(current, currentExitBack);
                std::set<uint64_t> distinctNeighbors;
                UnitigEntry nextEntry = {0, false};
                for(const auto& e : links) {
                    if(distinctNeighbors.insert(e.unitigIndex).second) {
                        nextEntry = e;
                    }
                }
                if(distinctNeighbors.size() != 1) break;

                uint64_t next = nextEntry.unitigIndex;
                bool cycle = false;
                for(uint64_t c : chain) {
                    if(c == next) { cycle = true; break; }
                }
                if(cycle) break;

                bool nextExitBack = nextEntry.orientation;
                if(sideNeighborCount(next, !nextExitBack) != 1) break;

                chain.push_back(next);
                totalLen += unitigs[next].totalOffset;
                current = next;
                currentExitBack = nextExitBack;
            }

            if(totalLen <= maxTipLength) {
                candidates.push_back({totalLen, std::move(chain)});
            }
        }

        if(candidates.empty()) break;

        // Sort shortest-first.
        std::sort(candidates.begin(), candidates.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

        // Remove non-overlapping tips in batch.
        std::set<uint64_t> removedThisRound;
        uint64_t chainsThisRound = 0;

        for(const auto& [len, chain] : candidates) {
            // Skip if any unitig in this chain was already removed.
            bool overlap = false;
            for(uint64_t idx : chain) {
                if(removedThisRound.count(idx)) { overlap = true; break; }
            }
            if(overlap) continue;

            // Remove all edges of all anchors in the chain.
            for(uint64_t idx : chain) {
                removedThisRound.insert(idx);
                const auto& u = unitigs[idx];
                for(const auto& oa : u.chain) {
                    auto fwdNeighbors = getNeighbors(oa);
                    for(const auto& nbr : fwdNeighbors) {
                        removeEdgeBothDirections(oa, nbr);
                    }
                    auto revNeighbors = getNeighbors(reverseAnchor(oa));
                    for(const auto& nbr : revNeighbors) {
                        removeEdgeBothDirections(reverseAnchor(oa), nbr);
                    }
                }
            }
            ++chainsThisRound;
        }

        totalChainsRemoved += chainsThisRound;
        totalUnitigsRemoved += removedThisRound.size();
    }

    cout << "removeTips: removed " << totalChainsRemoved << " tip chains ("
         << totalUnitigsRemoved << " unitigs, maxTipLength="
         << maxTipLength << ")." << endl;
    return totalChainsRemoved;
}


uint64_t BidirectedAnchorGraph::popSuperbubbles(uint64_t maxBubbleSize)
{
    // Doubled directed graph: each unitig i becomes two vertices
    // 2*i (i+) and 2*i+1 (i-). Edges follow Verkko's convention.

    using DGraph = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::bidirectionalS>;
    using dvertex = DGraph::vertex_descriptor;
    using dedge = pair<uint64_t, uint64_t>;

    auto unitigs = unitigify(/*quiet=*/true);
    const uint64_t n = unitigs.size();

    if(n == 0) {
        cout << "popSuperbubbles: no unitigs." << endl;
        return 0;
    }

    // Build entry map (same as writeUnitigGfa).
    struct UnitigEntryPoint {
        uint64_t unitigIndex;
        bool orientation;  // true = +, false = -
    };
    map<OrientedAnchor, UnitigEntryPoint> entryMap;
    for(uint64_t i = 0; i < n; i++) {
        const auto& front = unitigs[i].chain.front();
        const auto& back = unitigs[i].chain.back();
        entryMap[front]               = {i, true};
        entryMap[reverseAnchor(back)]  = {i, false};
        entryMap[reverseAnchor(front)] = {i, false};
        entryMap[back]                 = {i, true};
    }

    // Directed vertex index: unitig i, orientation + -> 2*i, orientation - -> 2*i+1.
    auto dv = [](uint64_t unitigIdx, bool orient) -> uint64_t {
        return 2 * unitigIdx + (orient ? 0 : 1);
    };

    // Build directed graph.
    DGraph dg(2 * n);

    // Map directed edges back to bidirected anchor-level edges.
    struct AnchorEdge {
        OrientedAnchor exitAnchor;
        OrientedAnchor neighborAnchor;
    };
    map<dedge, AnchorEdge> directedToAnchor;
    map<dedge, double> directedEdgeCoverage;
    set<dedge> addedEdges;

    for(uint64_t i = 0; i < n; i++) {
        const auto& u = unitigs[i];

        // Links from back of unitig i (i+).
        {
            const OrientedAnchor exitAnchor = u.chain.back();
            auto neighbors = getNeighbors({exitAnchor.first, exitAnchor.second});
            for(const auto& nbr : neighbors) {
                EdgeProperties props;
                if(!getEdgeProperties(exitAnchor, nbr, props)) continue;
                if(!props.useForAssembly) continue;

                auto it = entryMap.find(nbr);
                if(it == entryMap.end()) continue;

                uint64_t from = dv(i, true);
                uint64_t to = dv(it->second.unitigIndex, it->second.orientation);
                if(addedEdges.insert({from, to}).second) {
                    boost::add_edge(from, to, dg);
                    directedToAnchor[{from, to}] = {exitAnchor, nbr};
                    directedEdgeCoverage[{from, to}] = double(props.coverage);
                }

                // RC mirror.
                uint64_t rcFrom = dv(it->second.unitigIndex, !it->second.orientation);
                uint64_t rcTo = dv(i, false);
                if(addedEdges.insert({rcFrom, rcTo}).second) {
                    boost::add_edge(rcFrom, rcTo, dg);
                    directedToAnchor[{rcFrom, rcTo}] = {
                        reverseAnchor(nbr), reverseAnchor(exitAnchor)};
                    directedEdgeCoverage[{rcFrom, rcTo}] = double(props.coverage);
                }
            }
        }

        // Links from front of unitig i (i-).
        {
            const OrientedAnchor exitAnchor = reverseAnchor(u.chain.front());
            auto neighbors = getNeighbors({exitAnchor.first, exitAnchor.second});
            for(const auto& nbr : neighbors) {
                EdgeProperties props;
                if(!getEdgeProperties(exitAnchor, nbr, props)) continue;
                if(!props.useForAssembly) continue;

                auto it = entryMap.find(nbr);
                if(it == entryMap.end()) continue;

                uint64_t from = dv(i, false);
                uint64_t to = dv(it->second.unitigIndex, it->second.orientation);
                if(addedEdges.insert({from, to}).second) {
                    boost::add_edge(from, to, dg);
                    directedToAnchor[{from, to}] = {exitAnchor, nbr};
                    directedEdgeCoverage[{from, to}] = double(props.coverage);
                }

                uint64_t rcFrom = dv(it->second.unitigIndex, !it->second.orientation);
                uint64_t rcTo = dv(i, true);
                if(addedEdges.insert({rcFrom, rcTo}).second) {
                    boost::add_edge(rcFrom, rcTo, dg);
                    directedToAnchor[{rcFrom, rcTo}] = {
                        reverseAnchor(nbr), reverseAnchor(exitAnchor)};
                    directedEdgeCoverage[{rcFrom, rcTo}] = double(props.coverage);
                }
            }
        }
    }

    cout << "popSuperbubbles: directed graph: " << num_vertices(dg)
         << " vertices, " << num_edges(dg) << " edges from "
         << n << " unitigs." << endl;

    // Onodera et al. 2013 superbubble finder (inline, no BGL macros).
    // Uses remainingIncoming counting like findSuperbubbleOnoderaGfa.
    // Onodera et al. 2013, matching Verkko's find_bubble exactly.
    // Uses has_nonvisited_parent check instead of remainingIncoming counting.
    auto findBubble = [&](dvertex vStart) -> dvertex {
        if(out_degree(vStart, dg) < 2) return DGraph::null_vertex();

        std::unordered_set<dvertex> visited;
        std::unordered_set<dvertex> seen;
        vector<dvertex> S;

        S.push_back(vStart);
        seen.insert(vStart);

        while(!S.empty()) {
            if(maxBubbleSize > 0 && (visited.size() + seen.size()) > maxBubbleSize) {
                return DGraph::null_vertex();
            }

            dvertex v = S.back();
            S.pop_back();
            seen.erase(v);
            visited.insert(v);

            if(out_degree(v, dg) == 0) return DGraph::null_vertex();

            auto [oeBegin, oeEnd] = out_edges(v, dg);
            for(auto oeIt = oeBegin; oeIt != oeEnd; ++oeIt) {
                dvertex u = boost::target(*oeIt, dg);

                // Self-edge.
                if(u == v) return DGraph::null_vertex();
                // RC of u already visited (bidirected check).
                // In the doubled directed graph, RC of vertex 2*i is 2*i+1 and vice versa.
                dvertex rcU = (u ^ 1);
                if(visited.count(rcU)) return DGraph::null_vertex();
                // Back to start.
                if(u == vStart) return DGraph::null_vertex();
                // Already visited (back-edge).
                if(visited.count(u)) return DGraph::null_vertex();

                seen.insert(u);

                // Check if all parents of u have been visited.
                bool hasNonvisitedParent = false;
                auto [ieBegin, ieEnd] = in_edges(u, dg);
                for(auto ieIt = ieBegin; ieIt != ieEnd; ++ieIt) {
                    dvertex parent = boost::source(*ieIt, dg);
                    if(!visited.count(parent)) {
                        hasNonvisitedParent = true;
                        break;
                    }
                }
                if(!hasNonvisitedParent) {
                    S.push_back(u);
                }
            }

            if(S.size() == 1 && seen.size() == 1 && S[0] == *seen.begin()) {
                dvertex t = S[0];
                // Reject if exit has edge back to start.
                auto [teBegin, teEnd] = out_edges(t, dg);
                for(auto teIt = teBegin; teIt != teEnd; ++teIt) {
                    if(boost::target(*teIt, dg) == vStart) {
                        return DGraph::null_vertex();
                    }
                }
                return t;
            }
        }

        return DGraph::null_vertex();
    };

    // Find all superbubbles.
    struct Superbubble {
        dvertex source;
        dvertex target;
    };
    vector<Superbubble> superbubbles;

    uint64_t branchCount = 0;
    for(dvertex v = 0; v < num_vertices(dg); v++) {
        if(out_degree(v, dg) < 2) continue;
        ++branchCount;
        dvertex t = findBubble(v);
        if(t != DGraph::null_vertex()) {
            superbubbles.push_back({v, t});
        }
    }

    cout << "popSuperbubbles: " << branchCount << " branch vertices, found "
         << superbubbles.size() << " superbubbles." << endl;

    if(superbubbles.empty()) return 0;

    // For each superbubble, enumerate paths, keep widest, mark non-kept edges.
    set<dedge> globalKeptEdges;
    set<dedge> candidatesForRemoval;

    for(const auto& bubble : superbubbles) {
        // Collect bubble vertices via BFS from source.
        set<dvertex> bubbleVertices;
        {
            vector<dvertex> stack = {bubble.source};
            while(!stack.empty()) {
                dvertex v = stack.back();
                stack.pop_back();
                if(!bubbleVertices.insert(v).second) continue;
                if(v == bubble.target) continue;
                auto [oeBegin, oeEnd] = out_edges(v, dg);
                for(auto oeIt = oeBegin; oeIt != oeEnd; ++oeIt) {
                    dvertex w = boost::target(*oeIt, dg);
                    if(!bubbleVertices.count(w)) {
                        stack.push_back(w);
                    }
                }
            }
        }

        // Enumerate all paths from source to target via DFS.
        vector<vector<dedge>> allPaths;
        {
            struct Frame {
                dvertex vertex;
                vector<dedge> path;
            };
            vector<Frame> stack;
            stack.push_back({bubble.source, {}});

            while(!stack.empty()) {
                auto [v, path] = std::move(stack.back());
                stack.pop_back();

                if(v == bubble.target) {
                    allPaths.push_back(std::move(path));
                    continue;
                }

                auto [oeBegin, oeEnd] = out_edges(v, dg);
                for(auto oeIt = oeBegin; oeIt != oeEnd; ++oeIt) {
                    dvertex w = boost::target(*oeIt, dg);
                    if(w == bubble.target || bubbleVertices.count(w)) {
                        vector<dedge> newPath = path;
                        newPath.push_back({uint64_t(v), uint64_t(w)});
                        stack.push_back({w, std::move(newPath)});
                    }
                }
            }
        }

        if(allPaths.empty()) continue;

        // Find widest path (highest minimum edge coverage).
        uint64_t bestIdx = 0;
        double bestMinCov = 0.;
        for(uint64_t i = 0; i < allPaths.size(); i++) {
            double minCov = std::numeric_limits<double>::max();
            for(const auto& de : allPaths[i]) {
                auto it = directedEdgeCoverage.find(de);
                if(it != directedEdgeCoverage.end()) {
                    minCov = std::min(minCov, it->second);
                } else {
                    minCov = 0.;
                }
            }
            if(minCov > bestMinCov) {
                bestMinCov = minCov;
                bestIdx = i;
            }
        }

        for(const auto& de : allPaths[bestIdx]) {
            globalKeptEdges.insert(de);
        }
        for(uint64_t i = 0; i < allPaths.size(); i++) {
            if(i == bestIdx) continue;
            for(const auto& de : allPaths[i]) {
                candidatesForRemoval.insert(de);
            }
        }
    }

    // Remove edges not on any best path.
    set<dedge> edgesToRemove;
    for(const auto& de : candidatesForRemoval) {
        if(!globalKeptEdges.count(de)) {
            edgesToRemove.insert(de);
        }
    }

    // Map back to bidirected graph and remove.
    uint64_t removedCount = 0;
    for(const auto& de : edgesToRemove) {
        auto it = directedToAnchor.find(de);
        if(it == directedToAnchor.end()) continue;
        const auto& ae = it->second;
        if(hasEdge(ae.exitAnchor, ae.neighborAnchor)) {
            removeEdgeBothDirections(ae.exitAnchor, ae.neighborAnchor);
            ++removedCount;
        }
    }

    cout << "popSuperbubbles: popped " << superbubbles.size()
         << " superbubbles, removed " << removedCount
         << " edges." << endl;
    return removedCount;
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
    // Since RC-mirror unitigs were deduplicated, we register four
    // entry points per unitig (the unitig's own front/back plus
    // the RC mirror's front/back):
    //
    //   front                → i+  (enter unitig forward from left)
    //   reverseAnchor(back)  → i-  (enter unitig reversed from right)
    //   reverseAnchor(front) → i-  (RC mirror's back entry)
    //   back                 → i+  (RC mirror's front entry)
    struct UnitigEntryPoint {
        uint64_t unitigIndex;
        bool orientation;  // true = +, false = -
    };
    map<OrientedAnchor, UnitigEntryPoint> entryMap;

    for(uint64_t i = 0; i < unitigs.size(); i++) {
        const auto& front = unitigs[i].chain.front();
        const auto& back = unitigs[i].chain.back();

        entryMap[front] = {i, true};
        entryMap[reverseAnchor(back)] = {i, false};
        entryMap[reverseAnchor(front)] = {i, false};
        entryMap[back] = {i, true};
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

                auto it = entryMap.find(neighbor);
                if(it != entryMap.end()) {
                    emitLink(i, true, it->second.unitigIndex,
                             it->second.orientation, props);
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

                auto it = entryMap.find(neighbor);
                if(it != entryMap.end()) {
                    emitLink(i, false, it->second.unitigIndex,
                             it->second.orientation, props);
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
        const double hue = fmod(double(wid) * 137.508, 360.0);
        csv << i << "," << hslToHex(hue, 0.7, 0.5) << "\n";
    }
}


void BidirectedAnchorGraph::writeCsvByRead(const string& fileName, uint32_t readCount) const
{
    ofstream csv(fileName);
    if(!csv) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    csv << "Name,Color\n";

    if(readCount == 0) return;

    set<uint64_t> activeNodes;
    for(const auto& [key, props] : edgeProperties) {
        if(!props.useForAssembly) continue;
        activeNodes.insert(key.first.first.value());
        activeNodes.insert(key.second.first.value());
    }

    for(const auto nodeIdx : activeNodes) {
        if(nodeIdx >= nodeProps.size()) continue;
        const uint32_t rid = nodeProps[nodeIdx].backboneReadId;
        if(rid == UINT32_MAX) continue;

        // Golden angle spacing for perceptually distinct colors.
        const double hue = fmod(double(rid) * 137.508, 360.0);
        csv << nodeIdx << "," << hslToHex(hue, 0.7, 0.5) << "\n";
    }
}


void BidirectedAnchorGraph::writeUnitigCsvByRead(
    const string& fileName,
    const vector<Unitig>& unitigs,
    uint32_t readCount) const
{
    ofstream csv(fileName);
    if(!csv) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    csv << "Name,Color\n";

    if(readCount == 0) return;

    for(uint64_t i = 0; i < unitigs.size(); i++) {
        const auto& u = unitigs[i];
        if(u.chain.empty()) continue;

        const auto idx = u.chain.front().first.value();
        if(idx >= nodeProps.size()) continue;
        const uint32_t rid = nodeProps[idx].backboneReadId;
        if(rid == UINT32_MAX) continue;

        const double hue = fmod(double(rid) * 137.508, 360.0);
        csv << i << "," << hslToHex(hue, 0.7, 0.5) << "\n";
    }
}



// Bypass detour filter for the bidirected anchor graph.
//
// For each window w, walk its backbone anchors in order. For each
// neighbor window x that has edges both into and out of w's backbone,
// find pairs where x enters at backbone position i and exits at
// position j > i. Create a bypass edge connecting x's anchors on
// either side of the detour, and remove the inter-window edges.
//
// In the bidirected graph, backbone[i] = (node_i, orient_i):
// - Exit side (orient_i): getNeighbors gives nodes reachable going
//   forward along the backbone. Inter-window neighbors here are
//   "exits from w to x".
// - Entry side (!orient_i): getNeighbors gives nodes reachable going
//   backward. These are RC mirrors of edges arriving from x.
//   The real x-side source is reverseAnchor(neighbor).
//
// Bypass edge: reverseAnchor(inE.xNeighbor) -> outE.xNeighbor
// Edges to remove:
//   Entry: {bb[i].node, !bb[i].orient} -> inE.xNeighbor (and RC mirror)
//   Exit:  {bb[j].node,  bb[j].orient} -> outE.xNeighbor (and RC mirror)
uint64_t BidirectedAnchorGraph::trimBackbones(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Journeys& journeys)
{
    const uint32_t windowCount = uint32_t(anchorWindows.size());
    static constexpr uint32_t noWindow = UINT32_MAX;

    uint64_t trimmedAnchors = 0;
    uint64_t trimmedWindows = 0;

    for(uint32_t w = 0; w < windowCount; w++) {
        const auto& window = anchorWindows[w];
        const auto& positions = window.filteredBackbonePositions;
        if(positions.size() < 2) continue;

        const auto journey = journeys[window.backboneOrientedReadId];

        vector<OrientedAnchor> backbone;
        for(const uint32_t pos : positions) {
            backbone.push_back(toNormalizedOrientedAnchor(journey[pos]));
        }

        // Find the first and last backbone positions that have an
        // inter-window neighbor (on either side).
        auto hasInterWindowEdge = [&](uint32_t i) -> bool {
            const auto& bb = backbone[i];
            // Check exit side.
            for(const auto& nbr : getNeighbors({bb.first, bb.second})) {
                auto nIdx = nbr.first.value();
                if(nIdx < nodeProps.size()) {
                    uint32_t nWin = nodeProps[nIdx].windowId;
                    if(nWin != noWindow && nWin != w) return true;
                }
            }
            // Check entry side.
            for(const auto& nbr : getNeighbors({bb.first, !bb.second})) {
                auto nIdx = nbr.first.value();
                if(nIdx < nodeProps.size()) {
                    uint32_t nWin = nodeProps[nIdx].windowId;
                    if(nWin != noWindow && nWin != w) return true;
                }
            }
            return false;
        };

        uint32_t firstIW = UINT32_MAX;
        uint32_t lastIW = 0;
        for(uint32_t i = 0; i < uint32_t(backbone.size()); i++) {
            if(hasInterWindowEdge(i)) {
                if(firstIW == UINT32_MAX) firstIW = i;
                lastIW = i;
            }
        }

        // No inter-window edges at all — remove all edges (isolated window).
        if(firstIW == UINT32_MAX) {
            for(uint32_t i = 0; i < uint32_t(backbone.size()); i++) {
                const auto& bb = backbone[i];
                auto exitNbrs = getNeighbors({bb.first, bb.second});
                for(const auto& nbr : exitNbrs) {
                    removeEdgeBothDirections({bb.first, bb.second}, nbr);
                }
                auto entryNbrs = getNeighbors({bb.first, !bb.second});
                for(const auto& nbr : entryNbrs) {
                    removeEdgeBothDirections({bb.first, !bb.second}, nbr);
                }
                ++trimmedAnchors;
            }
            ++trimmedWindows;
            continue;
        }

        // Remove all edges of backbone anchors before firstIW and after lastIW.
        bool trimmed = false;
        for(uint32_t i = 0; i < firstIW; i++) {
            const auto& bb = backbone[i];
            auto exitNbrs = getNeighbors({bb.first, bb.second});
            for(const auto& nbr : exitNbrs) {
                removeEdgeBothDirections({bb.first, bb.second}, nbr);
            }
            auto entryNbrs = getNeighbors({bb.first, !bb.second});
            for(const auto& nbr : entryNbrs) {
                removeEdgeBothDirections({bb.first, !bb.second}, nbr);
            }
            ++trimmedAnchors;
            trimmed = true;
        }
        for(uint32_t i = lastIW + 1; i < uint32_t(backbone.size()); i++) {
            const auto& bb = backbone[i];
            auto exitNbrs = getNeighbors({bb.first, bb.second});
            for(const auto& nbr : exitNbrs) {
                removeEdgeBothDirections({bb.first, bb.second}, nbr);
            }
            auto entryNbrs = getNeighbors({bb.first, !bb.second});
            for(const auto& nbr : entryNbrs) {
                removeEdgeBothDirections({bb.first, !bb.second}, nbr);
            }
            ++trimmedAnchors;
            trimmed = true;
        }
        if(trimmed) ++trimmedWindows;
    }

    cout << "Trim backbones (bidirected): " << trimmedWindows
         << " windows trimmed, " << trimmedAnchors
         << " anchors trimmed." << endl;

    return trimmedAnchors;
}


uint64_t BidirectedAnchorGraph::bypassDetourFilter(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Journeys& journeys)
{
    const uint32_t windowCount = uint32_t(anchorWindows.size());
    static constexpr uint32_t noWindow = UINT32_MAX;

    uint64_t detoursFixed = 0;

    for(uint32_t w = 0; w < windowCount; w++) {
        const auto& window = anchorWindows[w];
        const auto& positions = window.filteredBackbonePositions;
        if(positions.size() < 2) continue;

        const auto journey = journeys[window.backboneOrientedReadId];

        // Build ordered backbone as OrientedAnchors (post-normalization).
        vector<OrientedAnchor> backbone;
        for(const uint32_t pos : positions) {
            backbone.push_back(toNormalizedOrientedAnchor(journey[pos]));
        }

        // Collect inter-window edges at each backbone position.
        struct IWEdge {
            OrientedAnchor xNeighbor;  // neighbor as returned by getNeighbors
            OrientedAnchor wAnchor;    // the backbone anchor
            uint32_t bbIdx;
        };

        // Group by neighbor window.
        map<uint32_t, vector<IWEdge>> entryEdges;  // x enters w
        map<uint32_t, vector<IWEdge>> exitEdges;   // w exits to x

        for(uint32_t i = 0; i < uint32_t(backbone.size()); i++) {
            const auto& bb = backbone[i];

            // Exit side: forward neighbors in other windows.
            for(const auto& nbr : getNeighbors({bb.first, bb.second})) {
                auto nIdx = nbr.first.value();
                if(nIdx >= nodeProps.size()) continue;
                uint32_t nWin = nodeProps[nIdx].windowId;
                if(nWin == noWindow || nWin == w) continue;
                exitEdges[nWin].push_back({nbr, bb, i});
            }

            // Entry side: backward neighbors in other windows.
            // These are RC mirrors of edges arriving from x.
            for(const auto& nbr : getNeighbors({bb.first, !bb.second})) {
                auto nIdx = nbr.first.value();
                if(nIdx >= nodeProps.size()) continue;
                uint32_t nWin = nodeProps[nIdx].windowId;
                if(nWin == noWindow || nWin == w) continue;
                entryEdges[nWin].push_back({nbr, bb, i});
            }
        }

        // For each neighbor window x with both entry and exit edges,
        // find detour pairs (entry at i, exit at j > i).
        for(const auto& [xWin, entries] : entryEdges) {
            auto exitIt = exitEdges.find(xWin);
            if(exitIt == exitEdges.end()) continue;
            const auto& exits = exitIt->second;

            // Sort exits by backbone index so we can find the closest
            // exit after each entry.
            auto sortedExits = exits;
            sort(sortedExits.begin(), sortedExits.end(),
                 [](const IWEdge& a, const IWEdge& b) {
                     return a.bbIdx < b.bbIdx;
                 });

            for(const auto& inE : entries) {
                // Find the closest bbIdx after this entry, then process
                // all exits at that bbIdx (there may be multiple X nodes
                // connecting at the same backbone position).
                uint32_t closestBbIdx = UINT32_MAX;
                for(const auto& outE : sortedExits) {
                    if(outE.bbIdx <= inE.bbIdx) continue;
                    closestBbIdx = outE.bbIdx;
                    break;
                }
                if(closestBbIdx == UINT32_MAX) continue;

                OrientedAnchor bypassFrom = reverseAnchor(inE.xNeighbor);

                for(const auto& outE : sortedExits) {
                    if(outE.bbIdx != closestBbIdx) continue;

                    OrientedAnchor bypassTo = outE.xNeighbor;

                    // Don't create self-loops.
                    if(bypassFrom.first == bypassTo.first) continue;

                    // Create bypass edge (both directions).
                    if(!hasEdge(bypassFrom, bypassTo)) {
                        EdgeProperties props;
                        props.useForAssembly = true;
                        addEdge(bypassFrom, bypassTo, props);
                        addTraversal(reverseAnchor(bypassTo), reverseAnchor(bypassFrom));
                    }

                    // Remove the exit edge (idempotent).
                    removeEdgeBothDirections(
                        {outE.wAnchor.first, outE.wAnchor.second},
                        outE.xNeighbor);
                }

                // Remove the entry edge.
                removeEdgeBothDirections(
                    {inE.wAnchor.first, !inE.wAnchor.second}, inE.xNeighbor);

                ++detoursFixed;
            }
        }
    }

    cout << "Bypass detour filter (bidirected): " << detoursFixed
         << " detours fixed." << endl;

    return detoursFixed;
}
