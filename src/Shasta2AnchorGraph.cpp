// Shasta.
#include "Shasta2AnchorGraph.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2AnchorPair.hpp"
#include "findSuperbubbleOnodera.hpp"
#include "GTest.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "ReadId.hpp"
#include "deduplicate.hpp"
#include "timestamp.hpp"
using namespace dinara;

namespace {
string anchorIdToString(Shasta2AnchorId anchorId)
{
    return shasta2AnchorIdToString(anchorId);
}
}

// Boost libraries.
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/serialization/vector.hpp>

// Standard library.
#include <algorithm>
#include <cstdlib>
#include <functional>
#include "fstream.hpp"
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "tuple.hpp"

// Explicit instantiation.
#include "MultithreadedObject.tpp"
namespace dinara {
    template class MultithreadedObject<Shasta2AnchorGraph>;
}



// Construct the Shasta2AnchorGraph from the Shasta2Journeys using
// the same edge creation rule as mode3::AnchorGraph:
// for each anchor, call Shasta2Anchors::findChildren and create one edge
// per child that satisfies minEdgeCoverage.
// The threadCount parameter is accepted for API compatibility but is not used here.
Shasta2AnchorGraph::Shasta2AnchorGraph(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    uint64_t minEdgeCoverage,
    uint64_t threadCount) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<Shasta2AnchorGraph>(*this)
{
    Shasta2AnchorGraph& anchorGraph = *this;
    static_cast<void>(threadCount);

    // Create the vertices, one for each AnchorId.
    // In the AnchorGraph, vertex_descriptors are AnchorIds.
    const uint64_t anchorCount = anchors.size();
    for(Shasta2AnchorId anchorId=0; anchorId<anchorCount; anchorId++) {
        add_vertex(anchorGraph);
    }

    nextEdgeId = 0;
    vector<Shasta2AnchorId> children;
    vector<uint64_t> counts;
    for(Shasta2AnchorId anchorIdA=0; anchorIdA<anchorCount; anchorIdA++) {
        anchors.findChildren(journeys, anchorIdA, children, counts, minEdgeCoverage);
        DINARA_ASSERT(children.size() == counts.size());
        for(uint64_t i=0; i<children.size(); i++) {
            const Shasta2AnchorId anchorIdB = children[i];
            Shasta2AnchorPair anchorPair(anchors, anchorIdA, anchorIdB, true);
            anchorPair.assertNoNegativeOffsets(anchors);
            if(anchorPair.orientedReadIds.empty()) {
                continue;
            }
            DINARA_ASSERT(anchors.countCommon(anchorIdA, anchorIdB) > 0);
            anchorGraph.addEdge(
                anchorPair.anchorIdA,
                anchorPair.anchorIdB,
                anchorPair.orientedReadIds,
                anchorPair.getAverageOffset(anchors));
        }
    }

    cout << "The anchor graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;
}



// Constructor from binary data.
Shasta2AnchorGraph::Shasta2AnchorGraph(const MappedMemoryOwner& mappedMemoryOwner, const string& name) :
    MappedMemoryOwner(mappedMemoryOwner),
    MultithreadedObject<Shasta2AnchorGraph>(*this)
{
    load(name);
}



void Shasta2AnchorGraph::save(ostream& s) const
{
    boost::archive::binary_oarchive archive(s);
    archive << *this;
}



void Shasta2AnchorGraph::load(istream& s)
{
    boost::archive::binary_iarchive archive(s);
    archive >> *this;
}



void Shasta2AnchorGraph::save(const string& name) const
{
    // If not using persistent binary data, do nothing.
    if(largeDataFileNamePrefix.empty()) {
        return;
    }

    // First save to a string.
    std::ostringstream s;
    save(s);
    const string dataString = s.str();

    // Now save the string to binary data.
    MemoryMapped::Vector<char> data;
    data.createNew(largeDataName(name), largeDataPageSize);
    data.resize(dataString.size());
    const char* begin = dataString.data();
    const char* end = begin + dataString.size();
    copy(begin, end, data.begin());
}



void Shasta2AnchorGraph::load(const string& name)
{
    // Access the binary data.
    MemoryMapped::Vector<char> data;
    try {
        data.accessExistingReadOnly(largeDataName(name));
    } catch (std::exception&) {
        throw runtime_error(name + " is not available.");
    }
    const string dataString(data.begin(), data.size());

    // Load it from here.
    std::istringstream s(dataString);
    try {
        load(s);
    } catch(std::exception& e) {
        throw runtime_error(string("Error reading " + name + ": ") + e.what());
    }
}



void Shasta2AnchorGraph::transitiveReduction(
    uint64_t transitiveReductionMaxEdgeCoverage,
    uint64_t transitiveReductionMaxDistance)
{
    Shasta2AnchorGraph& anchorGraph = *this;
    cout << "AnchorGraph transitive reduction begins." << endl;

    // Loop over edge coverage.
    // At each iteration we only consider edges with this coverage.
    vector<edge_descriptor> edgesToProcess;
    vector<edge_descriptor> edgesToRemove;
    for(uint64_t edgeCoverage=1; edgeCoverage<=transitiveReductionMaxEdgeCoverage; edgeCoverage++) {

        // Gather edges with this coverage.
        edgesToProcess.clear();
        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
            if(anchorGraph[e].useForAssembly and anchorGraph[e].coverage() == edgeCoverage) {
                edgesToProcess.push_back(e);
            }
        }

        // If there are none, there is nothing to do.
        if(edgesToProcess.empty()) {
            continue;
        }

        // Loop over all edges with this coverage.
        // This can be multithreaded.
        edgesToRemove.clear();
        for(const edge_descriptor e: edgesToProcess) {
            if(transitiveReductionCanRemove(e, transitiveReductionMaxDistance)) {
                edgesToRemove.push_back(e);
            }
        }

        // Turn off the useForAssembly flag for edges removed at this iteration over coverage.
        for(const edge_descriptor e: edgesToRemove) {
            disableEdge(e);
        }
        cout << "Edge coverage " << edgeCoverage <<
            ": processed " << edgesToProcess.size() <<
            " edges and flagged " << edgesToRemove.size() << endl;
    }
    cout << "AnchorGraph transitive reduction ends." << endl;

    uint64_t useForAssemblyCount = 0;
    BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
        if(anchorGraph[e].useForAssembly) {
            ++useForAssemblyCount;
        }
    }
    cout << useForAssemblyCount << " flagged for use in assembly out of " <<
        num_edges(anchorGraph) << " total." << endl;

}



uint64_t Shasta2AnchorGraph::cutWeakStalksLeadingToBranch(
    const Shasta2Anchors& anchors,
    uint64_t maxTipReadCount)
{
    // -----------------------------------------------------------------------
    // Post-transitive-reduction weak-stalk cutting on the assembly subgraph.
    //
    // We operate only on edges currently marked useForAssembly=true.
    //
    // A candidate stalk must:
    //   1. start at a tip in either directed orientation:
    //        - source-tip orientation: in-degree 0, walked forward
    //        - sink-tip orientation:   out-degree 0, walked backward
    //   2. follow a linear chain in that orientation until one of three stop conditions:
    //        a. the traversed chain reaches a branch point,
    //        b. the chain reaches a dead end,
    //        c. the union of oriented reads across all anchors seen so far
    //           exceeds maxTipReadCount.
    //
    // Cut rule agreed with the user:
    //   - If we hit a branch point while the union of supporting reads across
    //     the traversed stalk (excluding the terminal branch/merge anchor)
    //     is still <= maxTipReadCount, cut the whole chain.
    //   - If we hit a dead end, do not cut.
    //   - If the read union exceeds maxTipReadCount before reaching a branch
    //     point, stop and do not cut.
    //
    // "Cut the whole chain" means: mark all assembly edges in the traversed
    // prefix useForAssembly=false. Anchors themselves are not deleted.
    //
    // Important detail:
    //   The terminal branch/merge anchor is intentionally excluded from the
    //   read-union threshold. Otherwise, a weak low-read stalk that attaches
    //   into a high-coverage branch anchor would almost never satisfy the
    //   <= maxTipReadCount rule.
    // -----------------------------------------------------------------------

    Shasta2AnchorGraph& graph = *this;
    cout << "AnchorGraph weak-stalk branch cutting begins." << endl;

    auto inAssemblyDegree = [&](vertex_descriptor v) -> uint64_t {
        uint64_t degree = 0;
        BGL_FORALL_INEDGES(v, e, graph, Shasta2AnchorGraph) {
            if(graph[e].useForAssembly) {
                ++degree;
            }
        }
        return degree;
    };

    auto outAssemblyEdges = [&](vertex_descriptor v, vector<edge_descriptor>& edges) {
        edges.clear();
        BGL_FORALL_OUTEDGES(v, e, graph, Shasta2AnchorGraph) {
            if(graph[e].useForAssembly) {
                edges.push_back(e);
            }
        }
    };

    auto inAssemblyEdges = [&](vertex_descriptor v, vector<edge_descriptor>& edges) {
        edges.clear();
        BGL_FORALL_INEDGES(v, e, graph, Shasta2AnchorGraph) {
            if(graph[e].useForAssembly) {
                edges.push_back(e);
            }
        }
    };

    auto readUnionWithinThreshold = [&](
        const vector<vertex_descriptor>& chainVertices,
        uint64_t threshold) -> bool {
        std::unordered_set<uint64_t> orientedReadValues;
        orientedReadValues.reserve(threshold + 1);
        for(const vertex_descriptor v: chainVertices) {
            const Shasta2Anchor anchor = anchors[Shasta2AnchorId(v)];
            for(const auto& markerInfo: anchor) {
                orientedReadValues.insert(markerInfo.orientedReadId.getValue());
                if(orientedReadValues.size() > threshold) {
                    return false;
                }
            }
        }
        return true;
    };

    auto tryCollectWeakStalk = [&](
        vertex_descriptor vStart,
        bool forward,
        vector<edge_descriptor>& candidateEdgesToCut)
    {
        vector<edge_descriptor> assemblyOutEdges;
        vector<edge_descriptor> assemblyInEdges;
        vector<vertex_descriptor> chainVertices;
        vector<edge_descriptor> chainEdges;

        if(forward) {
            if(inAssemblyDegree(vStart) != 0) {
                return;
            }
            outAssemblyEdges(vStart, assemblyOutEdges);
            if(assemblyOutEdges.size() != 1) {
                return;
            }
        } else {
            outAssemblyEdges(vStart, assemblyOutEdges);
            if(!assemblyOutEdges.empty()) {
                return;
            }
            inAssemblyEdges(vStart, assemblyInEdges);
            if(assemblyInEdges.size() != 1) {
                return;
            }
        }

        chainVertices.clear();
        chainEdges.clear();
        chainVertices.push_back(vStart);

        vertex_descriptor current = vStart;
        bool shouldCut = false;

        while(true) {
            if(forward) {
                outAssemblyEdges(current, assemblyOutEdges);

                if(assemblyOutEdges.empty()) {
                    break;
                }
                if(assemblyOutEdges.size() > 1) {
                    shouldCut = !chainEdges.empty() && readUnionWithinThreshold(chainVertices, maxTipReadCount);
                    break;
                }

                const edge_descriptor e = assemblyOutEdges.front();
                const vertex_descriptor next = target(e, graph);
                const uint64_t nextInDegree = inAssemblyDegree(next);
                if(nextInDegree > 1) {
                    chainEdges.push_back(e);
                    shouldCut = true;
                    break;
                }

                chainEdges.push_back(e);
                chainVertices.push_back(next);

                if(!readUnionWithinThreshold(chainVertices, maxTipReadCount)) {
                    shouldCut = false;
                    break;
                }

                current = next;
            } else {
                inAssemblyEdges(current, assemblyInEdges);

                if(assemblyInEdges.empty()) {
                    break;
                }
                if(assemblyInEdges.size() > 1) {
                    shouldCut = !chainEdges.empty() && readUnionWithinThreshold(chainVertices, maxTipReadCount);
                    break;
                }

                const edge_descriptor e = assemblyInEdges.front();
                const vertex_descriptor previous = source(e, graph);
                vector<edge_descriptor> previousOutEdges;
                outAssemblyEdges(previous, previousOutEdges);
                if(previousOutEdges.size() > 1) {
                    chainEdges.push_back(e);
                    shouldCut = true;
                    break;
                }

                chainEdges.push_back(e);
                chainVertices.push_back(previous);

                if(!readUnionWithinThreshold(chainVertices, maxTipReadCount)) {
                    shouldCut = false;
                    break;
                }

                current = previous;
            }
        }

        if(shouldCut) {
            candidateEdgesToCut.insert(
                candidateEdgesToCut.end(),
                chainEdges.begin(),
                chainEdges.end());
        }
    };

    vector<edge_descriptor> candidateEdgesToCut;

    BGL_FORALL_VERTICES(vStart, graph, Shasta2AnchorGraph) {
        tryCollectWeakStalk(vStart, true, candidateEdgesToCut);
        tryCollectWeakStalk(vStart, false, candidateEdgesToCut);
    }

    uint64_t cutCount = 0;
    for(const edge_descriptor e: candidateEdgesToCut) {
        if(graph[e].useForAssembly) {
            disableEdge(e);
            ++cutCount;
        }
    }

    cout << "AnchorGraph weak-stalk branch cutting ends. Cut "
         << cutCount
         << " assembly edges (maxTipReadCount=" << maxTipReadCount << ")." << endl;
    return cutCount;
}



bool Shasta2AnchorGraph::transitiveReductionCanRemove(
    edge_descriptor e,
    uint64_t transitiveReductionMaxDistance) const
{
    const Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t edgeCoverage = anchorGraph[e].coverage();

    const vertex_descriptor v0 = source(e, anchorGraph);
    const vertex_descriptor v1 = target(e, anchorGraph);

    const bool debug = ((anchorIdToString(v0) == "45549+") and (anchorIdToString(v1) == "78505-"));

    // Do a forward BFS starting at v0, using edges
    // still marked as "use for assembly"
    // with coverage greater than edgeCoverage
    // and with maximum distance (number of edges)
    // equal to transitiveReductionMaxDistance.
    // If we encounter v1, return true.
    std::queue<vertex_descriptor> q;
    q.push(v0);

    // A map to store vertices already encountered and their distance from v0.
    std::map<vertex_descriptor, uint64_t> m;
    m.insert(make_pair(v0, 0));



    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor vA = q.front();
        q.pop();
        const auto itA = m.find(vA);
        DINARA_ASSERT(itA != m.end());
        const uint64_t distanceA = itA->second;
        const uint64_t distanceB = distanceA + 1;

        // Loop over its out-edges still marked as useForAssembly
        // and with sufficient coverage.
        BGL_FORALL_OUTEDGES(vA, eAB, anchorGraph, Shasta2AnchorGraph) {
            const Shasta2AnchorGraphEdge& edgeAB = anchorGraph[eAB];
            if(not edgeAB.useForAssembly) {
                continue;
            }

            // Only use edges with higher coverage for the BFS,
            if(edgeAB.coverage() <= edgeCoverage) {
                continue;
            }

            // If we reached v1, return true;
            const vertex_descriptor vB = target(eAB, anchorGraph);
            if(vB == v1) {
                if(debug) {
                    cout << "Edge " << anchorIdToString(v0) << " " << anchorIdToString(v1) <<
                        " flagged by transitive reduction." << endl;
                }
                return true;
            }

            // If we already encountered vB, don't do anything.
            if(m.contains(vB)) {
                continue;
            }

            if(distanceB < transitiveReductionMaxDistance) {
                q.push(vB);
                m.insert(make_pair(vB, distanceB));
            }
        }
    }

    // If getting here we did not encounter v1 in the BFS loop.
    if(debug) {
        cout << "Edge " << anchorIdToString(v0) << " " << anchorIdToString(v1) <<
            " not flagged by transitive reduction." << endl;
    }
    return false;
}



void Shasta2AnchorGraph::disableEdge(edge_descriptor e)
{
    Shasta2AnchorGraph& anchorGraph = *this;
    if(!anchorGraph[e].useForAssembly) return;
    anchorGraph[e].useForAssembly = false;

    const uint64_t srcVal = uint64_t(source(e, anchorGraph));
    const uint64_t dstVal = uint64_t(target(e, anchorGraph));
    const uint64_t anchorCount = num_vertices(anchorGraph);
    const uint64_t rcSrc = dstVal ^ 1ULL;
    const uint64_t rcDst = srcVal ^ 1ULL;
    if(rcSrc < anchorCount && rcDst < anchorCount) {
        auto [eit, exists] = boost::edge(rcSrc, rcDst, anchorGraph);
        if(exists) {
            anchorGraph[eit].useForAssembly = false;
        }
    }
}


uint64_t Shasta2AnchorGraph::windowTransitiveReduction()
{
    Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = num_vertices(anchorGraph);

    auto normalizeW = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // Build a window-level adjacency list and edge index.
    using WindowPairKey = std::pair<uint32_t, uint32_t>;
    std::set<WindowPairKey> windowEdges;
    std::map<WindowPairKey, vector<edge_descriptor>> windowPairEdges;
    std::map<uint32_t, std::set<uint32_t>> windowAdj;

    BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraphBaseClass) {
        if(!anchorGraph[e].useForAssembly) continue;
        const uint64_t src = uint64_t(source(e, anchorGraph));
        const uint64_t dst = uint64_t(target(e, anchorGraph));
        if(src >= anchorCount || dst >= anchorCount) continue;
        const uint32_t srcWin = anchorToWindow[src];
        const uint32_t dstWin = anchorToWindow[dst];
        if(srcWin == noWindow || dstWin == noWindow) continue;
        const uint32_t srcNorm = normalizeW(srcWin);
        const uint32_t dstNorm = normalizeW(dstWin);
        if(srcNorm == dstNorm) continue;
        windowEdges.insert({srcNorm, dstNorm});
        windowPairEdges[{srcNorm, dstNorm}].push_back(e);
        windowAdj[srcNorm].insert(dstNorm);
    }

    // For each window edge A→C, check if there exists a window B
    // such that A→B and B→C both exist (two-hop path).
    // If so, A→C is a transitive edge and can be removed.
    std::set<WindowPairKey> redundant;

    for(const auto& [wA, wC] : windowEdges) {
        auto itA = windowAdj.find(wA);
        if(itA == windowAdj.end()) continue;

        for(const uint32_t wB : itA->second) {
            if(wB == wC) continue;
            // Check if B→C exists.
            auto itB = windowAdj.find(wB);
            if(itB != windowAdj.end() && itB->second.count(wC)) {
                redundant.insert({wA, wC});
                break;
            }
        }
    }

    // Disable all anchor-level edges for redundant window pairs.
    uint64_t removedCount = 0;
    for(const auto& key : redundant) {
        auto it = windowPairEdges.find(key);
        if(it != windowPairEdges.end()) {
            for(const edge_descriptor e : it->second) {
                if(anchorGraph[e].useForAssembly) {
                    disableEdge(e);
                    ++removedCount;
                }
            }
        }
    }

    cout << "windowTransitiveReduction: found " << redundant.size()
         << " redundant window pairs, removed "
         << removedCount << " edges." << endl;
    return removedCount;
}


uint64_t Shasta2AnchorGraph::removeRcWindowConnections()
{
    Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = num_vertices(anchorGraph);

    uint64_t removedCount = 0;
    BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraphBaseClass) {
        if(!anchorGraph[e].useForAssembly) continue;
        const uint64_t src = uint64_t(source(e, anchorGraph));
        const uint64_t dst = uint64_t(target(e, anchorGraph));
        if(src >= anchorCount || dst >= anchorCount) continue;
        const uint32_t srcWin = anchorToWindow[src];
        const uint32_t dstWin = anchorToWindow[dst];
        if(srcWin == noWindow || dstWin == noWindow) continue;
        if(srcWin == dstWin) continue;

        // Check if src and dst are in RC-paired windows.
        // Window w and w + windowCount are RC counterparts.
        const bool srcIsRc = (srcWin >= windowCount);
        const bool dstIsRc = (dstWin >= windowCount);
        const uint32_t srcNorm = srcIsRc ? (srcWin - windowCount) : srcWin;
        const uint32_t dstNorm = dstIsRc ? (dstWin - windowCount) : dstWin;
        if(srcNorm == dstNorm) {
            // src is in window W, dst is in window W' (or vice versa).
            disableEdge(e);
            ++removedCount;
        }
    }

    cout << "removeRcWindowConnections: removed " << removedCount
         << " edges between RC window pairs." << endl;
    return removedCount;
}



void Shasta2AnchorGraph::verifyAgainstJourneys(
    const Shasta2Journeys& journeys,
    uint64_t minEdgeCoverage) const
{
    const char* const env = std::getenv("DINARA_VERIFY_ANCHOR_GRAPH");
    if((env == nullptr) or (env[0] == '0')) {
        return;
    }

    cout << timestamp << "Verifying anchor graph construction "
            "against an independent journey walk..." << endl;

    // Tally every journey-consecutive (anchorA, anchorB) pair by walking the
    // journeys directly. This deliberately uses none of the machinery the
    // construction path uses, so a shared bug cannot hide from the check.
    std::unordered_map<uint64_t, uint64_t> independentTally;
    const uint64_t orientedReadCount = journeys.size();
    for(uint64_t v = 0; v < orientedReadCount; v++) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(v);
        const auto journey = journeys[orientedReadId];
        for(uint64_t i = 0; i + 1 < journey.size(); i++) {
            const uint64_t a = uint64_t(journey[i]);
            const uint64_t b = uint64_t(journey[i + 1]);
            independentTally[(a << 32) | b]++;
        }
    }

    // Every edge in the graph must match the tally.
    uint64_t edgesChecked = 0, edgesMismatched = 0, edgesUnsupported = 0;
    const auto edgeRange = boost::edges(*this);
    for(auto it = edgeRange.first; it != edgeRange.second; ++it) {
        const auto& edge = (*this)[*it];
        const uint64_t a = edge.anchorIdA;
        const uint64_t b = edge.anchorIdB;
        const uint64_t declaredCoverage = edge.coverage();
        const auto tallyIt = independentTally.find((a << 32) | b);
        ++edgesChecked;
        if(tallyIt == independentTally.end()) {
            ++edgesUnsupported;
            cout << timestamp << "  MISMATCH: edge " << a << "->" << b
                 << " coverage=" << declaredCoverage
                 << " but independent walk found ZERO occurrences." << endl;
        } else if(tallyIt->second != declaredCoverage) {
            ++edgesMismatched;
            cout << timestamp << "  MISMATCH: edge " << a << "->" << b
                 << " coverage=" << declaredCoverage
                 << " independent tally=" << tallyIt->second << endl;
        }
    }

    // ...and every adjacency that clears minEdgeCoverage must have an edge.
    uint64_t missingEdges = 0;
    for(const auto& [key, count] : independentTally) {
        if(count < minEdgeCoverage) {
            continue;
        }
        const Shasta2AnchorId a = key >> 32;
        const Shasta2AnchorId b = key & 0xffffffffULL;
        const auto found = boost::edge(a, b, *this);
        if(!found.second) {
            ++missingEdges;
            if(missingEdges <= 20) {
                cout << timestamp << "  MISSING EDGE: " << a << "->" << b
                     << " independent tally=" << count
                     << " (>= minEdgeCoverage=" << minEdgeCoverage
                     << ") but no graph edge exists." << endl;
            }
        }
    }

    cout << timestamp << "Anchor graph verification: " << edgesChecked
         << " graph edges checked (" << edgesMismatched << " coverage mismatches, "
         << edgesUnsupported << " with zero independent support), "
         << independentTally.size() << " distinct journey-consecutive pairs found, "
         << missingEdges << " missing edges (tally >= minEdgeCoverage but no graph edge)."
         << endl;
}


uint64_t Shasta2AnchorGraph::removeHetArmTips(const Shasta2Anchors& anchors)
{
    Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = num_vertices(anchorGraph);
    const Shasta2AnchorId hetFirst = anchors.hetAnchorFirstId;
    if(hetFirst == invalid<Shasta2AnchorId>) {
        return 0;  // No het/hom anchors: nothing to do.
    }

    // Active in/out degree of a vertex (edges with useForAssembly == true).
    auto activeInDegree = [&](uint64_t a) -> uint64_t {
        uint64_t n = 0;
        auto ie = boost::in_edges(a, anchorGraph);
        for(auto it = ie.first; it != ie.second; ++it) {
            if(anchorGraph[*it].useForAssembly) ++n;
        }
        return n;
    };
    auto activeOutDegree = [&](uint64_t a) -> uint64_t {
        uint64_t n = 0;
        auto oe = boost::out_edges(a, anchorGraph);
        for(auto it = oe.first; it != oe.second; ++it) {
            if(anchorGraph[*it].useForAssembly) ++n;
        }
        return n;
    };

    // A het/hom anchor is a graph-interior node: it must have active edges on
    // BOTH sides. One-sided (or fully isolated-but-was-connected) => tip.
    // Isolated (0/0) het/hom anchors are staging artifacts with no edges to
    // disable, so they are not tips for our purposes (nothing to remove).
    auto isHetArmTip = [&](uint64_t a) -> bool {
        if(a < hetFirst || a >= anchorCount) return false;
        const uint64_t in = activeInDegree(a);
        const uint64_t out = activeOutDegree(a);
        if(in == 0 && out == 0) return false;   // already fully disconnected
        return (in == 0) || (out == 0);          // connected on only one side
    };

    // Disable all active edges incident to a vertex. disableEdge() also
    // disables the RC twin edge, so the mirror strand stays consistent.
    auto disableAllEdges = [&](uint64_t a) {
        if(a >= anchorCount) return;
        auto oe = boost::out_edges(a, anchorGraph);
        for(auto it = oe.first; it != oe.second; ++it) {
            if(anchorGraph[*it].useForAssembly) disableEdge(*it);
        }
        auto ie = boost::in_edges(a, anchorGraph);
        for(auto it = ie.first; it != ie.second; ++it) {
            if(anchorGraph[*it].useForAssembly) disableEdge(*it);
        }
    };

    // Iterate to a fixpoint: disabling a tip arm's surviving edge can strand a
    // neighboring hom anchor, turning it into a new one-sided het/hom tip.
    // disableEdge() disables the RC twin edge too, so processing the forward
    // anchor keeps both strands consistent without a separate RC pass.
    uint64_t totalDisabled = 0;
    uint64_t armsRemoved = 0;
    for(;;) {
        uint64_t disabledThisPass = 0;
        for(uint64_t a = hetFirst; a < anchorCount; ++a) {
            if(!isHetArmTip(a)) continue;
            const uint64_t before =
                activeInDegree(a) + activeOutDegree(a);
            disableAllEdges(a);
            disabledThisPass += before;
            ++armsRemoved;
        }
        totalDisabled += disabledThisPass;
        if(disabledThisPass == 0) break;
    }

    cout << "Remove het-arm tips: disabled " << totalDisabled
         << " edge(s) across " << armsRemoved
         << " one-sided het/hom anchor(s)." << endl;
    return totalDisabled;
}

