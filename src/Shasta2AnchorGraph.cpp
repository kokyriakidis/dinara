// Shasta.
#include "Shasta2AnchorGraph.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2AnchorPair.hpp"
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
#include "fstream.hpp"
#include <limits>
#include <map>
#include <queue>
#include <set>
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
            anchorPair.removeNegativeOffsets(anchors);
            if(anchorPair.orientedReadIds.empty()) {
                continue;
            }
            DINARA_ASSERT(anchors.countCommon(anchorIdA, anchorIdB) > 0);
            edge_descriptor e;
            tie(e, ignore) = add_edge(
                anchorPair.anchorIdA,
                anchorPair.anchorIdB,
                Shasta2AnchorGraphEdge(anchorPair, anchorPair.getAverageOffset(anchors), nextEdgeId++),
                anchorGraph);
            anchorGraph[e].useForAssembly = true;
        }
    }

    cout << "The anchor graph has " << num_vertices(*this) <<
        " vertices and " << num_edges(*this) << " edges." << endl;
}



// Construct from anchor windows.
// Each window becomes a chain of its backbone anchors.
// Inter-window edges are discovered by walking read journeys.
Shasta2AnchorGraph::Shasta2AnchorGraph(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const vector<AnchorWindow>& anchorWindows,
    uint64_t minInterWindowCoverage,
    uint64_t threadCount,
    const Reads* reads,
    const vector<DetangleBypassEdge>* bypassEdges) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<Shasta2AnchorGraph>(*this)
{
    Shasta2AnchorGraph& anchorGraph = *this;
    static_cast<void>(threadCount);

    // Create vertices, one per anchor.
    const uint64_t anchorCount = anchors.size();
    for(Shasta2AnchorId anchorId = 0; anchorId < anchorCount; anchorId++) {
        add_vertex(anchorGraph);
    }

    nextEdgeId = 0;

    // Build anchorId -> windowId and anchorId -> position-in-backbone maps.
    // For each original window W (windowId), we also create a mirror RC window
    // (windowId + windowCount) whose backbone anchors are the RC (anchorId ^ 1)
    // of the original backbone anchors. This lets strand-1 reads discover
    // inter-window transitions through the RC windows.
    windowCount = uint32_t(anchorWindows.size());
    anchorToWindow.assign(anchorCount, noWindow);
    vector<uint32_t> anchorToBackbonePos(anchorCount, 0);
    for(uint32_t windowId = 0; windowId < windowCount; windowId++) {
        const AnchorWindow& window = anchorWindows[windowId];
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = journeys[backboneOid];

        // Use filtered backbone positions if available, otherwise all positions.
        const auto& positions = window.filteredBackbonePositions.empty()
            ? [&]() -> const vector<uint32_t>& {
                // Build a temporary vector of all positions (stored per-window).
                static thread_local vector<uint32_t> allPositions;
                allPositions.clear();
                for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                    allPositions.push_back(pos);
                }
                return allPositions;
            }()
            : window.filteredBackbonePositions;

        for(const uint32_t pos : positions) {
            const uint64_t aid = uint64_t(backboneJourney[pos]);
            anchorToWindow[aid] = windowId;
            anchorToBackbonePos[aid] = pos;
            // Mirror RC window: map the RC anchor to windowId + windowCount.
            const uint64_t rcAid = aid ^ 1ULL;
            if(rcAid < anchorCount) {
                anchorToWindow[rcAid] = windowId + windowCount;
                anchorToBackbonePos[rcAid] = pos;
            }
        }
    }

    // Helper to add an edge if the anchor pair has shared oriented reads.
    auto addEdgeIfValid = [&](Shasta2AnchorId anchorIdA, Shasta2AnchorId anchorIdB) -> bool {
        Shasta2AnchorPair anchorPair(anchors, anchorIdA, anchorIdB, false);
        anchorPair.removeNegativeOffsets(anchors);
        if(anchorPair.orientedReadIds.empty()) {
            return false;
        }
        DINARA_ASSERT(anchors.countCommon(anchorIdA, anchorIdB) > 0);
        edge_descriptor e;
        tie(e, ignore) = add_edge(
            anchorPair.anchorIdA,
            anchorPair.anchorIdB,
            Shasta2AnchorGraphEdge(anchorPair, anchorPair.getAverageOffset(anchors), nextEdgeId++),
            anchorGraph);
        anchorGraph[e].useForAssembly = true;
        return true;
    };

    // Intra-window edges: consecutive filtered backbone anchor pairs,
    // for both the original windows and their RC mirrors.
    for(const AnchorWindow& window : anchorWindows) {
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = journeys[backboneOid];

        // Collect the backbone anchor IDs for this window.
        vector<Shasta2AnchorId> backboneAnchors;
        const auto& positions = window.filteredBackbonePositions;
        if(!positions.empty()) {
            for(const uint32_t pos : positions) {
                backboneAnchors.push_back(backboneJourney[pos]);
            }
        } else {
            for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                backboneAnchors.push_back(backboneJourney[pos]);
            }
        }

        if(backboneAnchors.size() < 2) continue;

        for(uint64_t i = 0; i + 1 < backboneAnchors.size(); i++) {
            addEdgeIfValid(backboneAnchors[i], backboneAnchors[i + 1]);
            // RC mirror edge.
            const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(backboneAnchors[i]) ^ 1ULL);
            const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(backboneAnchors[i + 1]) ^ 1ULL);
            if(uint64_t(rcA) < anchorCount && uint64_t(rcB) < anchorCount) {
                addEdgeIfValid(rcB, rcA);
            }
        }
    }

    // Alternate path edges: for each het window, add a chain
    // anchorIdA -> intermediate[0] -> ... -> intermediate[N-1] -> anchorIdB.
    // These form parallel paths (bubbles) at het sites.
    // Only emit for windows with detected het SNPs.
    // Also create RC mirror edges.
    uint64_t alternatePathEdgeCount = 0;
    for(const AnchorWindow& window : anchorWindows) {
        if(window.cleanHetSnpCount == 0) continue;
        for(const AnchorWindowAlternatePath& altPath : window.alternatePaths) {
            // Build the forward chain: A -> intermediates -> B.
            Shasta2AnchorId prevAnchorId = altPath.anchorIdA;
            vector<Shasta2AnchorId> forwardChain;
            forwardChain.push_back(prevAnchorId);
            for(const Shasta2AnchorId midAnchorId : altPath.intermediateAnchorIds) {
                if(addEdgeIfValid(prevAnchorId, midAnchorId)) {
                    ++alternatePathEdgeCount;
                }
                forwardChain.push_back(midAnchorId);
                prevAnchorId = midAnchorId;
            }
            // Last edge: last intermediate -> anchorIdB.
            if(addEdgeIfValid(prevAnchorId, altPath.anchorIdB)) {
                ++alternatePathEdgeCount;
            }
            forwardChain.push_back(altPath.anchorIdB);

            // RC mirror: reverse the chain and flip each anchor ID.
            for(uint64_t i = forwardChain.size() - 1; i > 0; i--) {
                const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(forwardChain[i]) ^ 1ULL);
                const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(forwardChain[i-1]) ^ 1ULL);
                if(uint64_t(rcA) < anchorCount && uint64_t(rcB) < anchorCount) {
                    if(addEdgeIfValid(rcA, rcB)) {
                        ++alternatePathEdgeCount;
                    }
                }
            }
        }
    }

    // Per-read window transition tracking.
    // Walk each read's journey, collect the sequence of normalized windows,
    // then populate:
    //   - AnchorWindowReadInterval::previousWindow / nextWindow
    //   - AnchorWindow::transitionReads[(prev, next)] -> read list
    {
        auto normalize = [&](uint32_t w) -> uint32_t {
            return (w >= windowCount) ? (w - windowCount) : w;
        };
        const uint32_t noW = AnchorWindowReadInterval::noWindow;

        // Build lookup: (windowId, orientedReadId) -> index in readIntervals.
        std::map<std::pair<uint32_t, uint32_t>, uint64_t> windowReadIndex;
        for(uint32_t wid = 0; wid < windowCount; wid++) {
            auto& window = const_cast<AnchorWindow&>(anchorWindows[wid]);
            for(uint64_t ri = 0; ri < window.readIntervals.size(); ri++) {
                windowReadIndex[{wid, window.readIntervals[ri].orientedReadId.getValue()}] = ri;
            }
        }

        const uint64_t journeyCount = journeys.size();
        for(uint64_t oidValue = 0; oidValue < journeyCount; oidValue++) {
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            const auto journey = journeys[oid];
            if(journey.empty()) continue;

            // Collect the sequence of distinct normalized windows this read visits.
            std::vector<uint32_t> windowSequence;
            for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
                const Shasta2AnchorId anchorId = journey[pos];
                if(uint64_t(anchorId) >= anchorCount) continue;
                const uint32_t windowId = anchorToWindow[uint64_t(anchorId)];
                if(windowId == noWindow) continue;
                const uint32_t normW = normalize(windowId);
                if(windowSequence.empty() || windowSequence.back() != normW) {
                    windowSequence.push_back(normW);
                }
            }

            // For each window in the sequence, update per-read fields and transition map.
            for(uint64_t i = 0; i < windowSequence.size(); i++) {
                const uint32_t wid = windowSequence[i];
                const uint32_t prev = (i > 0) ? windowSequence[i - 1] : noW;
                const uint32_t next = (i + 1 < windowSequence.size()) ? windowSequence[i + 1] : noW;

                auto it = windowReadIndex.find({wid, oidValue});
                if(it != windowReadIndex.end()) {
                    auto& interval = const_cast<AnchorWindow&>(anchorWindows[wid]).readIntervals[it->second];
                    interval.previousWindow = prev;
                    interval.nextWindow = next;
                }

                // Always add to transitionReads (even if read isn't in readIntervals,
                // it still provides transition evidence).
                auto& w = const_cast<AnchorWindow&>(anchorWindows[wid]);
                w.transitionReads[{prev, next}].push_back(oid);


            }
        }
    }

    // Inter-window edges: walk each read's journey and collect candidate
    // anchor pairs (lastAnchorInWindowA, firstAnchorInWindowB) for each
    // ordered window pair. Pick the candidate with the highest shared
    // read count (commonReadCount).

    // For each window pair, collect all candidate (anchorA, anchorB) pairs.
    struct AnchorPairKey {
        Shasta2AnchorId anchorIdA;
        Shasta2AnchorId anchorIdB;
        bool operator<(const AnchorPairKey& o) const {
            if(anchorIdA != o.anchorIdA) return anchorIdA < o.anchorIdA;
            return anchorIdB < o.anchorIdB;
        }
    };
    std::map<std::pair<uint32_t, uint32_t>,
             std::map<AnchorPairKey, uint32_t>> windowPairCandidates;

    uint64_t containedSkipCount = 0;
    const uint64_t journeyCount = journeys.size();
    for(uint64_t oidValue = 0; oidValue < journeyCount; oidValue++) {
        const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
        const auto journey = journeys[oid];
        if(journey.empty()) continue;

        // Skip contained reads for inter-window edge discovery.
        if(reads) {
            const ReadId readId = oid.getReadId();
            if(readId < reads->readCount() && reads->getFlags(readId).isContained) {
                ++containedSkipCount;
                continue;
            }
        }

        uint32_t currentWindow = noWindow;
        Shasta2AnchorId lastAnchorInCurrentWindow = 0;

        for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
            const Shasta2AnchorId anchorId = journey[pos];
            if(uint64_t(anchorId) >= anchorCount) continue;
            const uint32_t windowId = anchorToWindow[uint64_t(anchorId)];
            if(windowId == noWindow) continue;

            if(windowId == currentWindow) {
                lastAnchorInCurrentWindow = anchorId;
            } else {
                if(currentWindow != noWindow) {
                    auto key = std::make_pair(currentWindow, windowId);
                    AnchorPairKey apk{lastAnchorInCurrentWindow, anchorId};
                    windowPairCandidates[key][apk]++;
                }
                currentWindow = windowId;
                lastAnchorInCurrentWindow = anchorId;
            }
        }
    }

    // Diagnostic: count window pairs and candidates.
    {
        uint64_t totalCandidateEdges = 0;
        uint64_t fwFwPairs = 0, fwRcPairs = 0, rcFwPairs = 0, rcRcPairs = 0;
        for(const auto& [wp, cands] : windowPairCandidates) {
            totalCandidateEdges += cands.size();
            const bool aIsRc = (wp.first >= windowCount);
            const bool bIsRc = (wp.second >= windowCount);
            if(!aIsRc && !bIsRc) ++fwFwPairs;
            else if(!aIsRc && bIsRc) ++fwRcPairs;
            else if(aIsRc && !bIsRc) ++rcFwPairs;
            else ++rcRcPairs;
        }
        cout << "Inter-window edge discovery: " << windowPairCandidates.size()
             << " window pairs (" << fwFwPairs << " fw-fw, "
             << rcRcPairs << " rc-rc, "
             << fwRcPairs << " fw-rc, "
             << rcFwPairs << " rc-fw) with "
             << totalCandidateEdges << " candidate anchor pairs." << endl;

        // Count how many anchors are mapped vs total.
        uint64_t mappedCount = 0;
        for(uint64_t i = 0; i < anchorCount; i++) {
            if(anchorToWindow[i] != noWindow) mappedCount++;
        }
        cout << "anchorToWindow: " << mappedCount << " of " << anchorCount
             << " anchors mapped to " << windowCount << " original + "
             << windowCount << " RC mirror windows." << endl;
        if(reads) {
            cout << "Skipped " << containedSkipCount
                 << " contained oriented reads for inter-window edge discovery." << endl;
        }
    }

    // Identify parallel windows: windows where the set of incoming neighbor
    // windows is identical to the set of outgoing neighbor windows.
    // These represent parallel connections (e.g., haplotype bubbles) and
    // should not have inter-window edges added.
    std::set<uint32_t> parallelWindows;
    {
        const uint32_t noW = AnchorWindowReadInterval::noWindow;
        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& window = anchorWindows[w];
            std::set<uint32_t> inNeighbors, outNeighbors;
            for(const auto& [key, reads] : window.transitionReads) {
                if(key.first != noW) inNeighbors.insert(key.first);
                if(key.second != noW) outNeighbors.insert(key.second);
            }
            if(inNeighbors.size() == 2 && inNeighbors == outNeighbors) {
                parallelWindows.insert(w);
            }
        }
        if(!parallelWindows.empty()) {
            cout << "Parallel windows (same in/out neighbors, edges skipped):";
            for(const uint32_t w : parallelWindows) {
                cout << " " << w;
            }
            cout << endl;
        }
    }

    // For each window pair, pick the candidate with the most shared reads.
    uint64_t interWindowZeroPairs = 0;
    uint64_t interWindowBelowCoverage = 0;
    uint64_t interWindowParallel = 0;
    uint64_t interWindowCreated = 0;
    // Track created inter-window edges: (windowPair, anchorIdA, anchorIdB, readCount).
    struct InterWindowEdgeInfo {
        std::pair<uint32_t, uint32_t> windowPair;
        Shasta2AnchorId anchorIdA;
        Shasta2AnchorId anchorIdB;
        uint64_t readCount;
    };
    std::vector<InterWindowEdgeInfo> createdEdges;
    for(const auto& [windowPair, candidates] : windowPairCandidates) {
        // Skip edges involving parallel windows.
        const uint32_t normSrc = (windowPair.first >= windowCount)
            ? (windowPair.first - windowCount) : windowPair.first;
        const uint32_t normDst = (windowPair.second >= windowCount)
            ? (windowPair.second - windowCount) : windowPair.second;
        if(parallelWindows.count(normSrc) || parallelWindows.count(normDst)) {
            ++interWindowParallel;
            continue;
        }

        Shasta2AnchorPair bestPair;
        uint64_t bestSize = 0;
        for(const auto& [apk, count] : candidates) {
            Shasta2AnchorPair anchorPair(anchors, apk.anchorIdA, apk.anchorIdB, false);
            anchorPair.removeNegativeOffsets(anchors);
            if(anchorPair.size() > bestSize) {
                bestSize = anchorPair.size();
                bestPair = std::move(anchorPair);
            }
        }
        if(bestSize == 0) {
            ++interWindowZeroPairs;
        } else if(bestSize < minInterWindowCoverage) {
            ++interWindowBelowCoverage;
        } else {
            DINARA_ASSERT(anchors.countCommon(bestPair.anchorIdA, bestPair.anchorIdB) > 0);

            edge_descriptor e;
            tie(e, ignore) = add_edge(
                bestPair.anchorIdA,
                bestPair.anchorIdB,
                Shasta2AnchorGraphEdge(bestPair, bestPair.getAverageOffset(anchors), nextEdgeId++),
                anchorGraph);
            anchorGraph[e].useForAssembly = true;
            createdEdges.push_back({windowPair, bestPair.anchorIdA, bestPair.anchorIdB, bestSize});
            ++interWindowCreated;
        }
    }
    cout << "Inter-window edges: " << interWindowCreated << " created, "
         << interWindowParallel << " rejected (parallel window), "
         << interWindowZeroPairs << " rejected (zero forward-flow reads), "
         << interWindowBelowCoverage << " rejected (below minInterWindowCoverage="
         << minInterWindowCoverage << ")." << endl;

    // Populate per-window outEdges/inEdges from createdEdges.
    for(const auto& edgeInfo : createdEdges) {
        const uint32_t srcW = edgeInfo.windowPair.first;
        const uint32_t dstW = edgeInfo.windowPair.second;
        // anchorWindows only has windowCount entries (forward windows).
        // For RC mirror windows (>= windowCount), skip — they don't have
        // AnchorWindow entries. The RC edges are already tracked via the
        // forward window's edges + anchor ^ 1.
        if(srcW < windowCount) {
            auto& w = const_cast<AnchorWindow&>(anchorWindows[srcW]);
            w.outEdges.push_back({dstW, edgeInfo.anchorIdA, edgeInfo.anchorIdB, edgeInfo.readCount});
        }
        if(dstW < windowCount) {
            auto& w = const_cast<AnchorWindow&>(anchorWindows[dstW]);
            w.inEdges.push_back({srcW, edgeInfo.anchorIdA, edgeInfo.anchorIdB, edgeInfo.readCount});
        }
    }

    // Per-window connectivity summary using transitionReads.
    {
        const uint32_t noW = AnchorWindow::noWindow;
        cout << "Per-window connectivity (per-read transitions):" << endl;
        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& window = anchorWindows[w];

            // Aggregate incoming and outgoing counts from transitionReads.
            std::map<uint32_t, uint64_t> incomingCounts, outgoingCounts;
            uint64_t totalTransitionReads = 0;
            for(const auto& [key, reads] : window.transitionReads) {
                totalTransitionReads += reads.size();
                if(key.first != noW) {
                    incomingCounts[key.first] += reads.size();
                }
                if(key.second != noW) {
                    outgoingCounts[key.second] += reads.size();
                }
            }

            cout << "  Window " << w
                 << " (" << totalTransitionReads << " reads): incoming=[";
            bool first = true;
            for(const auto& [fromW, count] : incomingCounts) {
                if(!first) cout << ", ";
                cout << fromW << ":" << count;
                first = false;
            }
            cout << "] outgoing=[";
            first = true;
            for(const auto& [toW, count] : outgoingCounts) {
                if(!first) cout << ", ";
                cout << toW << ":" << count;
                first = false;
            }
            cout << "]";

            // Show transition flows if there are multiple distinct patterns.
            if(window.transitionReads.size() > 1) {
                cout << " flows={";
                first = true;
                for(const auto& [key, reads] : window.transitionReads) {
                    if(!first) cout << ", ";
                    cout << "(";
                    if(key.first == noW) cout << "-"; else cout << key.first;
                    cout << "→";
                    if(key.second == noW) cout << "-"; else cout << key.second;
                    cout << "):" << reads.size();
                    first = false;
                }
                cout << "}";
            }
            cout << endl;
        }
    }

    // Bypass edges from detangling: direct connections that skip over
    // tangled windows whose backbone reads have been removed.
    uint64_t bypassEdgeCount = 0;
    if(bypassEdges) {
        for(const auto& be : *bypassEdges) {
            if(addEdgeIfValid(be.anchorIdA, be.anchorIdB)) {
                ++bypassEdgeCount;
            }
        }
        if(bypassEdgeCount > 0) {
            cout << "Bypass edges: " << bypassEdgeCount << " created from "
                 << bypassEdges->size() << " candidates." << endl;
        }
    }

    // ========================================================================
    // Post-construction graph cleanup rules.
    // After creating inter-window edges, apply rules to remove artifacts.
    // More rules can be added below as needed.
    //
    // Rule 1: Trim backbone outside the bounding inter-window connection span.
    //   For each forward window, find the journey positions of all inter-window
    //   connection anchors (outgoing anchorIdA, incoming anchorIdB). The useful
    //   backbone span is bounded by:
    //     - Start: min of all incoming connection positions (or backbone start
    //       if no incoming edges)
    //     - End: max of all outgoing connection positions (or backbone end
    //       if no outgoing edges)
    //   Trim backbone vertices outside [start, end].
    //   Only intra-window edges are removed to avoid fragmenting the graph.
    // ========================================================================
    {
        // Remove only intra-window edges of a vertex (edges where both
        // endpoints belong to the same window). Preserves inter-window edges.
        auto clearIntraWindowEdges = [&](uint64_t vid) {
            if(vid >= anchorCount) return;
            const uint32_t vWindow = anchorToWindow[vid];

            std::vector<edge_descriptor> toRemove;

            auto oe = boost::out_edges(vid, anchorGraph);
            for(auto it = oe.first; it != oe.second; ++it) {
                const uint64_t tgt = uint64_t(boost::target(*it, anchorGraph));
                if(tgt < anchorCount && anchorToWindow[tgt] == vWindow) {
                    toRemove.push_back(*it);
                }
            }
            auto ie = boost::in_edges(vid, anchorGraph);
            for(auto it = ie.first; it != ie.second; ++it) {
                const uint64_t src = uint64_t(boost::source(*it, anchorGraph));
                if(src < anchorCount && anchorToWindow[src] == vWindow) {
                    toRemove.push_back(*it);
                }
            }

            for(const auto& e : toRemove) {
                boost::remove_edge(e, anchorGraph);
            }
        };

        // Helper: find the journey position of an anchor in a window's
        // backbone journey span [backboneBegin, backboneEnd).
        auto findJourneyPos = [&](
            const auto& journey,
            const AnchorWindow& window,
            Shasta2AnchorId anchor) -> int64_t
        {
            for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                if(journey[pos] == anchor) return int64_t(pos);
            }
            return -1;
        };

        uint64_t trimmedVertexCount = 0;
        uint64_t trimmedWindowCount = 0;

        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& window = anchorWindows[w];
            const auto& positions = window.filteredBackbonePositions;
            if(positions.empty()) continue;
            if(window.outEdges.empty() && window.inEdges.empty()) continue;

            const auto journey = journeys[window.backboneOrientedReadId];

            // Find bounding journey positions from inter-window edges.
            // Start bound: min of incoming connection positions.
            // End bound: max of outgoing connection positions.
            int64_t startBound = -1;  // earliest incoming
            int64_t endBound = -1;    // latest outgoing

            for(const auto& ie : window.inEdges) {
                const int64_t pos = findJourneyPos(journey, window, ie.anchorIdB);
                if(pos >= 0) {
                    if(startBound < 0 || pos < startBound) startBound = pos;
                }
            }

            for(const auto& oe : window.outEdges) {
                const int64_t pos = findJourneyPos(journey, window, oe.anchorIdA);
                if(pos >= 0) {
                    if(endBound < 0 || pos > endBound) endBound = pos;
                }
            }

            // If no incoming edges, start at backbone start (no head trimming).
            if(startBound < 0) startBound = int64_t(positions.front());
            // If no outgoing edges, end at backbone end (no tail trimming).
            if(endBound < 0) endBound = int64_t(positions.back());

            // Sanity: if start > end, skip (shouldn't happen in well-formed graph).
            if(startBound > endBound) continue;

            // Find the first filtered position >= startBound (trim before it).
            int64_t keepFirst = -1;
            for(uint64_t i = 0; i < positions.size(); i++) {
                if(int64_t(positions[i]) >= startBound) {
                    keepFirst = int64_t(i);
                    break;
                }
            }

            // Find the last filtered position <= endBound (trim after it).
            int64_t keepLast = -1;
            for(int64_t i = int64_t(positions.size()) - 1; i >= 0; i--) {
                if(int64_t(positions[i]) <= endBound) {
                    keepLast = i;
                    break;
                }
            }

            if(keepFirst < 0 || keepLast < 0 || keepFirst > keepLast) continue;

            const uint64_t headTrim = uint64_t(keepFirst);
            const uint64_t tailTrim = positions.size() - 1 - uint64_t(keepLast);

            if(headTrim == 0 && tailTrim == 0) continue;

            ++trimmedWindowCount;

            cout << "  Rule1 W" << w
                 << " in=" << window.inEdges.size()
                 << " out=" << window.outEdges.size()
                 << " span=[" << startBound << "," << endBound << "]"
                 << " keep=[" << keepFirst << "," << keepLast << "]/" << positions.size()
                 << " head=" << headTrim << " tail=" << tailTrim
                 << endl;

            // Trim head (before keepFirst).
            for(uint64_t i = 0; i < uint64_t(keepFirst); i++) {
                const Shasta2AnchorId aid = journey[positions[i]];
                clearIntraWindowEdges(uint64_t(aid));
                const uint64_t rcAid = uint64_t(aid) ^ 1ULL;
                clearIntraWindowEdges(rcAid);
                ++trimmedVertexCount;
            }

            // Trim tail (after keepLast).
            for(uint64_t i = uint64_t(keepLast) + 1; i < positions.size(); i++) {
                const Shasta2AnchorId aid = journey[positions[i]];
                clearIntraWindowEdges(uint64_t(aid));
                const uint64_t rcAid = uint64_t(aid) ^ 1ULL;
                clearIntraWindowEdges(rcAid);
                ++trimmedVertexCount;
            }
        }

        cout << "Rule 1: " << trimmedWindowCount << " windows trimmed, "
             << trimmedVertexCount << " vertices trimmed." << endl;
    }

    // ========================================================================
    // Detangle Case 2: 2x2 tangle matrix for internal inter-window edges.
    // COMMENTED OUT pending validation on larger datasets.
    // ========================================================================
#if 0
    {
        // Build boundary anchor set.
        std::set<Shasta2AnchorId> boundaryAnchors;
        for(const AnchorWindow& window : anchorWindows) {
            const auto backboneJourney = journeys[window.backboneOrientedReadId];
            const auto& positions = window.filteredBackbonePositions;
            if(positions.empty()) continue;
            const Shasta2AnchorId firstAnchor = backboneJourney[positions.front()];
            const Shasta2AnchorId lastAnchor = backboneJourney[positions.back()];
            boundaryAnchors.insert(firstAnchor);
            boundaryAnchors.insert(lastAnchor);
            boundaryAnchors.insert(Shasta2AnchorId(uint64_t(firstAnchor) ^ 1ULL));
            boundaryAnchors.insert(Shasta2AnchorId(uint64_t(lastAnchor) ^ 1ULL));
        }

        // Build anchor → (windowIndex, backbone position index) map.
        std::map<Shasta2AnchorId, std::pair<uint32_t, uint64_t>> anchorPosition;
        for(uint32_t wIdx = 0; wIdx < anchorWindows.size(); wIdx++) {
            const auto& window = anchorWindows[wIdx];
            const auto backboneJourney = journeys[window.backboneOrientedReadId];
            const auto& positions = window.filteredBackbonePositions;
            for(uint64_t i = 0; i < positions.size(); i++) {
                const Shasta2AnchorId aid = backboneJourney[positions[i]];
                anchorPosition[aid] = {wIdx, i};
            }
        }

        // Collect reads from consecutive linear backbone anchors,
        // counting how many anchors (steps) each read appears on.
        // Walks from startPosIdx in the given direction (forward or
        // backward) and stops when it hits an anchor that has inter-window
        // edges (degree disturbance) or reaches the window boundary.
        // Returns a map: readId -> stepCount.
        auto collectLinearReadsWithSteps = [&](
            const AnchorWindow& window,
            uint64_t startPosIdx,
            bool walkForward) -> std::map<uint32_t, uint64_t>
        {
            std::map<uint32_t, uint64_t> readSteps;
            const auto journey = journeys[window.backboneOrientedReadId];
            const auto& positions = window.filteredBackbonePositions;

            if(walkForward) {
                for(uint64_t i = startPosIdx; i < positions.size(); i++) {
                    const Shasta2AnchorId aid = journey[positions[i]];
                    const uint64_t aidVal = uint64_t(aid);

                    // Stop if this anchor has inter-window edges
                    // (i.e., it connects to a different window).
                    if(i != startPosIdx) {
                        bool hasInterWindowEdge = false;
                        BGL_FORALL_OUTEDGES(aid, outE, anchorGraph, Shasta2AnchorGraph) {
                            const uint64_t tgt = uint64_t(target(outE, anchorGraph));
                            const uint32_t tgtWin = (tgt < anchorCount) ? anchorToWindow[tgt] : noWindow;
                            const uint32_t aidWin = (aidVal < anchorCount) ? anchorToWindow[aidVal] : noWindow;
                            if(tgtWin != noWindow && tgtWin != aidWin) {
                                hasInterWindowEdge = true;
                                break;
                            }
                        }
                        if(!hasInterWindowEdge) {
                            BGL_FORALL_INEDGES(aid, inE, anchorGraph, Shasta2AnchorGraph) {
                                const uint64_t src = uint64_t(source(inE, anchorGraph));
                                const uint32_t srcWin2 = (src < anchorCount) ? anchorToWindow[src] : noWindow;
                                const uint32_t aidWin = (aidVal < anchorCount) ? anchorToWindow[aidVal] : noWindow;
                                if(srcWin2 != noWindow && srcWin2 != aidWin) {
                                    hasInterWindowEdge = true;
                                    break;
                                }
                            }
                        }
                        if(hasInterWindowEdge) break;
                    }

                    const auto span = anchors[aid];
                    for(const auto& mi : span) {
                        readSteps[mi.orientedReadId.getValue()]++;
                    }
                }
            } else {
                // Walk backward.
                for(int64_t i = int64_t(startPosIdx); i >= 0; i--) {
                    const Shasta2AnchorId aid = journey[positions[uint64_t(i)]];
                    const uint64_t aidVal = uint64_t(aid);

                    if(uint64_t(i) != startPosIdx) {
                        bool hasInterWindowEdge = false;
                        BGL_FORALL_OUTEDGES(aid, outE, anchorGraph, Shasta2AnchorGraph) {
                            const uint64_t tgt = uint64_t(target(outE, anchorGraph));
                            const uint32_t tgtWin = (tgt < anchorCount) ? anchorToWindow[tgt] : noWindow;
                            const uint32_t aidWin = (aidVal < anchorCount) ? anchorToWindow[aidVal] : noWindow;
                            if(tgtWin != noWindow && tgtWin != aidWin) {
                                hasInterWindowEdge = true;
                                break;
                            }
                        }
                        if(!hasInterWindowEdge) {
                            BGL_FORALL_INEDGES(aid, inE, anchorGraph, Shasta2AnchorGraph) {
                                const uint64_t src = uint64_t(source(inE, anchorGraph));
                                const uint32_t srcWin2 = (src < anchorCount) ? anchorToWindow[src] : noWindow;
                                const uint32_t aidWin = (aidVal < anchorCount) ? anchorToWindow[aidVal] : noWindow;
                                if(srcWin2 != noWindow && srcWin2 != aidWin) {
                                    hasInterWindowEdge = true;
                                    break;
                                }
                            }
                        }
                        if(hasInterWindowEdge) break;
                    }

                    const auto span = anchors[aid];
                    for(const auto& mi : span) {
                        readSteps[mi.orientedReadId.getValue()]++;
                    }
                }
            }
            return readSteps;
        };

        // Compute tangle matrix entry: sum of min(stepsInEntrance, stepsInExit)
        // for each read that appears in both, following shasta2's approach.
        auto tangleEntry = [](const std::map<uint32_t, uint64_t>& entrance,
                              const std::map<uint32_t, uint64_t>& exit) -> uint64_t {
            uint64_t total = 0;
            for(const auto& [readId, entranceSteps] : entrance) {
                auto it = exit.find(readId);
                if(it != exit.end()) {
                    total += std::min(entranceSteps, it->second);
                }
            }
            return total;
        };

        // Find internal inter-window edges and check the 2x2 matrix.
        uint64_t case2RemovedCount = 0;
        bool changed = true;
        while(changed) {
            changed = false;
            BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
                const Shasta2AnchorId srcAid = Shasta2AnchorId(uint64_t(source(e, anchorGraph)));
                const Shasta2AnchorId dstAid = Shasta2AnchorId(uint64_t(target(e, anchorGraph)));

                const uint64_t srcVal = uint64_t(srcAid);
                const uint64_t dstVal = uint64_t(dstAid);
                const uint32_t srcWin = (srcVal < anchorCount) ? anchorToWindow[srcVal] : noWindow;
                const uint32_t dstWin = (dstVal < anchorCount) ? anchorToWindow[dstVal] : noWindow;
                if(srcWin == noWindow || dstWin == noWindow) continue;
                if(srcWin == dstWin) continue;

                // Both endpoints must be internal.
                if(boundaryAnchors.count(srcAid)) continue;
                if(boundaryAnchors.count(dstAid)) continue;

                // Find backbone positions of both endpoints.
                auto srcIt = anchorPosition.find(srcAid);
                auto dstIt = anchorPosition.find(dstAid);
                if(srcIt == anchorPosition.end() || dstIt == anchorPosition.end()) continue;

                const uint32_t srcWIdx = srcIt->second.first;
                const uint64_t srcPosIdx = srcIt->second.second;
                const uint32_t dstWIdx = dstIt->second.first;
                const uint64_t dstPosIdx = dstIt->second.second;

                const auto& srcWindow = anchorWindows[srcWIdx];
                const auto& dstWindow = anchorWindows[dstWIdx];
                const auto& srcPositions = srcWindow.filteredBackbonePositions;
                const auto& dstPositions = dstWindow.filteredBackbonePositions;

                // Need at least one anchor before and after in both windows.
                if(srcPosIdx == 0 || srcPosIdx + 1 >= srcPositions.size()) continue;
                if(dstPosIdx == 0 || dstPosIdx + 1 >= dstPositions.size()) continue;

                // Collect reads with step counts from consecutive linear
                // backbone anchors before and after the connection point.
                // Walk stops at anchors with inter-window edges.
                // A-entrance: walk backward from srcPosIdx-1.
                // A-exit:     walk forward from srcPosIdx+1.
                // B-entrance: walk backward from dstPosIdx-1.
                // B-exit:     walk forward from dstPosIdx+1.
                const auto aEntrance = collectLinearReadsWithSteps(srcWindow, srcPosIdx - 1, false);
                const auto aExit     = collectLinearReadsWithSteps(srcWindow, srcPosIdx + 1, true);
                const auto bEntrance = collectLinearReadsWithSteps(dstWindow, dstPosIdx - 1, false);
                const auto bExit     = collectLinearReadsWithSteps(dstWindow, dstPosIdx + 1, true);

                // Build 2x2 tangle matrix with step-weighted entries.
                const uint64_t m00 = tangleEntry(aEntrance, aExit);  // A→A
                const uint64_t m01 = tangleEntry(aEntrance, bExit);  // A→B
                const uint64_t m10 = tangleEntry(bEntrance, aExit);  // B→A
                const uint64_t m11 = tangleEntry(bEntrance, bExit);  // B→B

                vector<vector<uint64_t>> tangleMatrix = {{m00, m01}, {m10, m11}};
                GTest gtest(tangleMatrix, 0.1, false, false);

                if(gtest.success && !gtest.hypotheses.empty()) {
                    const auto& best = gtest.hypotheses[0];
                    // The connection is false if the best hypothesis
                    // is diagonal: A→A and B→B only, no cross-connections.
                    const bool isDiagonal =
                        best.connectivityMatrix[0][0] &&
                        !best.connectivityMatrix[0][1] &&
                        !best.connectivityMatrix[1][0] &&
                        best.connectivityMatrix[1][1];

                    if(isDiagonal) {
                        boost::remove_edge(e, anchorGraph);
                        ++case2RemovedCount;
                        changed = true;
                        break; // Iterator invalidated, restart.
                    }
                }
            }
        }

        if(case2RemovedCount > 0) {
            cout << "Detangle Case 2: removed " << case2RemovedCount
                 << " false internal inter-window edges." << endl;
        }
    }
#endif

    // Validate: check that every edge has shared oriented reads.
    {
        uint64_t emptyEdgeCount = 0;
        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
            if(anchorGraph[e].anchorPair.orientedReadIds.empty()) {
                ++emptyEdgeCount;
            }
        }
        if(emptyEdgeCount > 0) {
            cout << "WARNING: " << emptyEdgeCount
                 << " edges have no shared oriented reads." << endl;
        }
    }

    // Recount edge types after trimming (clear_vertex may have removed edges).
    uint64_t finalIntraCount = 0;
    uint64_t finalAltPathCount = 0;
    uint64_t finalInterCount = 0;
    BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
        if(!anchorGraph[e].useForAssembly) {
            ++finalAltPathCount;
        } else {
            // Check if both endpoints belong to the same window.
            const uint64_t srcAnchor = uint64_t(source(e, anchorGraph));
            const uint64_t dstAnchor = uint64_t(target(e, anchorGraph));
            const uint32_t srcWin = (srcAnchor < anchorCount) ? anchorToWindow[srcAnchor] : noWindow;
            const uint32_t dstWin = (dstAnchor < anchorCount) ? anchorToWindow[dstAnchor] : noWindow;
            if(srcWin != noWindow && srcWin == dstWin) {
                ++finalIntraCount;
            } else {
                ++finalInterCount;
            }
        }
    }

    // Count non-isolated vertices (vertices with at least one edge).
    uint64_t nonIsolatedVertexCount = 0;
    for(Shasta2AnchorId v = 0; v < anchorCount; v++) {
        if(in_degree(v, anchorGraph) > 0 || out_degree(v, anchorGraph) > 0) {
            ++nonIsolatedVertexCount;
        }
    }

    cout << "The anchor graph has " << nonIsolatedVertexCount
         << " non-isolated vertices (" << num_vertices(*this) << " total), "
         << num_edges(*this) << " edges"
         << " (" << finalIntraCount << " intra-window, "
         << finalAltPathCount << " alternate-path, "
         << finalInterCount << " inter-window)." << endl;

    // Verify all edges: for each orientedReadId, check its base position
    // on both anchors. Report edges where any read has negative base offset
    // or is missing from one of the anchors.
    {
        uint64_t backwardEdgeCount = 0;
        uint64_t missingReadEdgeCount = 0;
        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraphBaseClass) {
            const auto& dEdge = anchorGraph[e];
            const auto& ap = dEdge.anchorPair;
            if(ap.orientedReadIds.empty()) continue;

            const Shasta2Anchor anchorA = anchors[ap.anchorIdA];
            const Shasta2Anchor anchorB = anchors[ap.anchorIdB];

            uint64_t forwardCount = 0;
            uint64_t backwardCount = 0;
            uint64_t onlyACount = 0;
            uint64_t onlyBCount = 0;
            uint64_t neitherCount = 0;

            auto itA = anchorA.begin();
            auto itB = anchorB.begin();

            for(const OrientedReadId orientedReadId : ap.orientedReadIds) {
                // Find this read on anchor A.
                while(itA != anchorA.end() && itA->orientedReadId < orientedReadId) ++itA;
                const bool isOnA = (itA != anchorA.end() && itA->orientedReadId == orientedReadId);

                // Find this read on anchor B.
                while(itB != anchorB.end() && itB->orientedReadId < orientedReadId) ++itB;
                const bool isOnB = (itB != anchorB.end() && itB->orientedReadId == orientedReadId);

                if(isOnA && isOnB) {
                    if(itB->position >= itA->position) {
                        ++forwardCount;
                    } else {
                        ++backwardCount;
                    }
                } else if(isOnA) {
                    ++onlyACount;
                } else if(isOnB) {
                    ++onlyBCount;
                } else {
                    ++neitherCount;
                }
            }

            if(backwardCount > 0 || neitherCount > 0) {
                if(backwardCount > 0) ++backwardEdgeCount;
                if(neitherCount > 0) ++missingReadEdgeCount;
                if(backwardEdgeCount + missingReadEdgeCount <= 10) {
                    cout << "EDGE CHECK " << ap.anchorIdA << " -> " << ap.anchorIdB
                         << ": forward=" << forwardCount
                         << " backward=" << backwardCount
                         << " onlyA=" << onlyACount
                         << " onlyB=" << onlyBCount
                         << " neither=" << neitherCount
                         << " total=" << ap.orientedReadIds.size()
                         << " (offset=" << dEdge.offset << ")" << endl;
                }
            }
        }
        cout << "Edge verification: " << backwardEdgeCount << " edges with backward reads, "
             << missingReadEdgeCount << " edges with reads on neither anchor." << endl;
    }
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
            anchorGraph[e].useForAssembly = false;
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
            graph[e].useForAssembly = false;
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




