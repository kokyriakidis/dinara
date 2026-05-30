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

    // Populate backbonePreviousWindow / backboneNextWindow for each window.
    for(uint32_t wid = 0; wid < windowCount; wid++) {
        auto& window = const_cast<AnchorWindow&>(anchorWindows[wid]);
        const uint32_t bbOid = window.backboneOrientedReadId.getValue();
        for(const auto& ri : window.readIntervals) {
            if(ri.orientedReadId.getValue() == bbOid) {
                window.backbonePreviousWindow = ri.previousWindow;
                window.backboneNextWindow = ri.nextWindow;
                break;
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

    // Build the set of endpoint window pairs from backbone transitions.
    auto normalize = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };
    endpointWindowPairs.clear();
    for(uint32_t wid = 0; wid < windowCount; wid++) {
        const auto& window = anchorWindows[wid];
        const uint32_t noW = AnchorWindowReadInterval::noWindow;
        if(window.backbonePreviousWindow != noW) {
            uint32_t prev = window.backbonePreviousWindow;
            endpointWindowPairs.insert({std::min(prev, wid), std::max(prev, wid)});
        }
        if(window.backboneNextWindow != noW) {
            uint32_t next = window.backboneNextWindow;
            endpointWindowPairs.insert({std::min(wid, next), std::max(wid, next)});
        }
    }

    // For each window pair, pick the candidate with the most shared reads.
    // Two passes: endpoint edges first (to reserve their anchors), then
    // remaining edges (skipping reserved anchors).
    uint64_t interWindowZeroPairs = 0;
    uint64_t interWindowBelowCoverage = 0;
    uint64_t interWindowCreated = 0;
    uint64_t interWindowEndpointCreated = 0;
    uint64_t interWindowInternalCreated = 0;
    uint64_t interWindowInternalSkipped = 0;
    // Track created inter-window edges: (windowPair, anchorIdA, anchorIdB, readCount).
    struct InterWindowEdgeInfo {
        std::pair<uint32_t, uint32_t> windowPair;
        Shasta2AnchorId anchorIdA;
        Shasta2AnchorId anchorIdB;
        uint64_t readCount;
    };
    std::vector<InterWindowEdgeInfo> createdEdges;

    // Anchors reserved by endpoint edges — other edges cannot use these.
    std::set<uint64_t> reservedAnchors;

    // Helper: create an edge from a chosen anchor pair.
    auto createInterWindowEdge = [&](
        const std::pair<uint32_t, uint32_t>& windowPair,
        Shasta2AnchorPair& bestPair,
        uint64_t bestSize)
    {
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
    };

    // Pass 1: Create endpoint edges and reserve their anchors.
    // Each raw window pair is processed independently (fw-fw, rc-rc, etc.).
    for(const auto& [windowPair, candidates] : windowPairCandidates) {
        const uint32_t srcNorm = normalize(windowPair.first);
        const uint32_t dstNorm = normalize(windowPair.second);
        if(srcNorm == dstNorm) continue;
        if(!endpointWindowPairs.count({std::min(srcNorm, dstNorm), std::max(srcNorm, dstNorm)})) continue;

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
            reservedAnchors.insert(uint64_t(bestPair.anchorIdA));
            reservedAnchors.insert(uint64_t(bestPair.anchorIdB));
            reservedAnchors.insert(uint64_t(bestPair.anchorIdA) ^ 1ULL);
            reservedAnchors.insert(uint64_t(bestPair.anchorIdB) ^ 1ULL);
            createInterWindowEdge(windowPair, bestPair, bestSize);
            ++interWindowEndpointCreated;
        }
    }

    // Store endpoint anchors for GFA tagging.
    endpointAnchors = reservedAnchors;

    // Early trim: disable backbone anchors beyond the endpoint anchors.
    // For each window, find the backbone positions of the endpoint anchors
    // (from pass 1). Disable all intra-window edges on anchors before the
    // head endpoint or after the tail endpoint, and add those anchors to
    // reservedAnchors so pass 2 can't use them.
    {
        // Per-window endpoint positions on the backbone.
        struct EndpointPos {
            uint32_t headPos = 0;           // position of head endpoint anchor
            uint32_t tailPos = UINT32_MAX;  // position of tail endpoint anchor
            bool hasHead = false;
            bool hasTail = false;
        };
        std::vector<EndpointPos> endpointPos(windowCount);

        for(const auto& edgeInfo : createdEdges) {
            const uint32_t srcNorm = normalize(edgeInfo.windowPair.first);
            const uint32_t dstNorm = normalize(edgeInfo.windowPair.second);

            if(srcNorm < windowCount) {
                const uint32_t pos = anchorToBackbonePos[uint64_t(edgeInfo.anchorIdA)];
                const auto& w = anchorWindows[srcNorm];
                if(w.backboneNextWindow == dstNorm) {
                    if(!endpointPos[srcNorm].hasTail || pos < endpointPos[srcNorm].tailPos) {
                        endpointPos[srcNorm].tailPos = pos;
                    }
                    endpointPos[srcNorm].hasTail = true;
                }
                if(w.backbonePreviousWindow == dstNorm) {
                    if(!endpointPos[srcNorm].hasHead || pos > endpointPos[srcNorm].headPos) {
                        endpointPos[srcNorm].headPos = pos;
                    }
                    endpointPos[srcNorm].hasHead = true;
                }
            }

            if(dstNorm < windowCount) {
                const uint32_t pos = anchorToBackbonePos[uint64_t(edgeInfo.anchorIdB)];
                const auto& w = anchorWindows[dstNorm];
                if(w.backbonePreviousWindow == srcNorm) {
                    if(!endpointPos[dstNorm].hasHead || pos > endpointPos[dstNorm].headPos) {
                        endpointPos[dstNorm].headPos = pos;
                    }
                    endpointPos[dstNorm].hasHead = true;
                }
                if(w.backboneNextWindow == srcNorm) {
                    if(!endpointPos[dstNorm].hasTail || pos < endpointPos[dstNorm].tailPos) {
                        endpointPos[dstNorm].tailPos = pos;
                    }
                    endpointPos[dstNorm].hasTail = true;
                }
            }
        }

        uint64_t earlyTrimCount = 0;
        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& ep = endpointPos[w];
            if(!ep.hasHead && !ep.hasTail) continue;

            const auto& window = anchorWindows[w];
            const auto& positions = window.filteredBackbonePositions;
            const auto journey = journeys[window.backboneOrientedReadId];

            for(const uint32_t pos : positions) {
                bool outsideBounds = false;
                if(ep.hasHead && pos < ep.headPos) outsideBounds = true;
                if(ep.hasTail && pos > ep.tailPos) outsideBounds = true;

                if(outsideBounds) {
                    const uint64_t aid = uint64_t(journey[pos]);
                    // Disable all edges on this anchor and its RC mirror.
                    auto disableAll = [&](uint64_t a) {
                        if(a >= anchorCount) return;
                        auto oe = boost::out_edges(a, anchorGraph);
                        for(auto it = oe.first; it != oe.second; ++it) {
                            if(anchorGraph[*it].useForAssembly)
                                disableEdge(*it);
                        }
                        auto ie = boost::in_edges(a, anchorGraph);
                        for(auto it = ie.first; it != ie.second; ++it) {
                            if(anchorGraph[*it].useForAssembly)
                                disableEdge(*it);
                        }
                    };
                    disableAll(aid);
                    disableAll(aid ^ 1ULL);
                    reservedAnchors.insert(aid);
                    reservedAnchors.insert(aid ^ 1ULL);
                    ++earlyTrimCount;
                }
            }
        }
        cout << "Early trim after endpoint edges: " << earlyTrimCount
             << " backbone anchors disabled." << endl;
    }

    // Pass 2: Create remaining (non-endpoint) edges, skipping reserved anchors.
    for(const auto& [windowPair, candidates] : windowPairCandidates) {
        const uint32_t srcNorm = normalize(windowPair.first);
        const uint32_t dstNorm = normalize(windowPair.second);
        if(srcNorm == dstNorm) continue;
        if(endpointWindowPairs.count({std::min(srcNorm, dstNorm), std::max(srcNorm, dstNorm)})) continue;

        Shasta2AnchorPair bestPair;
        uint64_t bestSize = 0;
        bool hadCandidates = false;
        for(const auto& [apk, count] : candidates) {
            hadCandidates = true;
            if(reservedAnchors.count(uint64_t(apk.anchorIdA)) ||
               reservedAnchors.count(uint64_t(apk.anchorIdB))) {
                continue;
            }
            Shasta2AnchorPair anchorPair(anchors, apk.anchorIdA, apk.anchorIdB, false);
            anchorPair.removeNegativeOffsets(anchors);
            if(anchorPair.size() > bestSize) {
                bestSize = anchorPair.size();
                bestPair = std::move(anchorPair);
            }
        }

        if(bestSize == 0) {
            if(hadCandidates) {
                ++interWindowInternalSkipped;
            } else {
                ++interWindowZeroPairs;
            }
        } else if(bestSize < minInterWindowCoverage) {
            ++interWindowBelowCoverage;
        } else {
            createInterWindowEdge(windowPair, bestPair, bestSize);
            ++interWindowInternalCreated;
        }
    }

    cout << "Inter-window edges: " << interWindowCreated << " created ("
         << interWindowEndpointCreated << " endpoint, "
         << interWindowInternalCreated << " internal), "
         << interWindowInternalSkipped << " internal skipped (reserved anchors), "
         << interWindowZeroPairs << " rejected (zero forward-flow reads), "
         << interWindowBelowCoverage << " rejected (below minInterWindowCoverage="
         << minInterWindowCoverage << ")." << endl;

    // ========================================================================
    // Filter lambdas (defined here, called in order below).
    // ========================================================================

    // Use the member function disableEdge() for all edge disabling.
    // (Defined at the end of this file.)

    // Recompute backbonePreviousWindow / backboneNextWindow for each window
    // based on current active edges. Must be called after any filter that
    // changes edge state before filters that use these fields.
    auto recomputeBackboneEndpoints = [&]() {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        for(uint32_t wid = 0; wid < windowCount; wid++) {
            auto& window = const_cast<AnchorWindow&>(anchorWindows[wid]);
            const auto journey = journeys[window.backboneOrientedReadId];

            // Walk the backbone read's full journey and build the sequence
            // of distinct normalized windows. Use anchorToWindow directly
            // (not filtered by active edges) because the backbone read's
            // path through windows is a property of the read, not the
            // current edge state. Filtering by active edges would miss
            // windows whose connecting anchors were trimmed.
            std::vector<uint32_t> windowSequence;
            std::vector<uint32_t> posToSeqIdx(journey.size(), UINT32_MAX);
            for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
                const uint64_t aid = uint64_t(journey[pos]);
                if(aid >= anchorCount) continue;
                const uint32_t rawW = anchorToWindow[aid];
                if(rawW == noWindow) continue;
                const uint32_t normW = normalize(rawW);
                if(windowSequence.empty() || windowSequence.back() != normW) {
                    windowSequence.push_back(normW);
                }
                posToSeqIdx[pos] = uint32_t(windowSequence.size() - 1);
            }

            // Find the occurrence of wid in the sequence that corresponds
            // to the backbone positions. Use the first backbone position
            // to locate the correct occurrence (the backbone read may
            // re-enter wid after visiting other windows, and we need the
            // occurrence that matches the backbone span).
            window.backbonePreviousWindow = AnchorWindowReadInterval::noWindow;
            window.backboneNextWindow = AnchorWindowReadInterval::noWindow;

            const auto& positions = window.filteredBackbonePositions;
            uint32_t seqIdx = UINT32_MAX;
            for(const uint32_t pos : positions) {
                if(pos < posToSeqIdx.size() && posToSeqIdx[pos] != UINT32_MAX) {
                    seqIdx = posToSeqIdx[pos];
                    break;
                }
            }

            if(seqIdx != UINT32_MAX && seqIdx < windowSequence.size() &&
               windowSequence[seqIdx] == wid) {
                if(seqIdx > 0) {
                    window.backbonePreviousWindow = windowSequence[seqIdx - 1];
                }
                if(seqIdx + 1 < windowSequence.size()) {
                    window.backboneNextWindow = windowSequence[seqIdx + 1];
                }
            }
        }
    };

    // Case 2: Remove single internal-edge connections between windows.
    // If two windows are connected by exactly one inter-window edge
    // and that edge is not an endpoint connection (not matching
    // backbonePreviousWindow or backboneNextWindow of either window),
    // it's a spurious single-point connection — delete it.
    auto runSingleEdgeFilter = [&]() {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        // Count inter-window edges per (normalized) window pair.
        // Also collect the edge descriptors.
        std::map<std::pair<uint32_t, uint32_t>, std::vector<edge_descriptor>> pairEdges;
        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
            if(!anchorGraph[e].useForAssembly) continue;
            const uint64_t srcVal = uint64_t(source(e, anchorGraph));
            const uint64_t dstVal = uint64_t(target(e, anchorGraph));
            if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
            const uint32_t srcRaw = anchorToWindow[srcVal];
            const uint32_t dstRaw = anchorToWindow[dstVal];
            if(srcRaw == noWindow || dstRaw == noWindow) continue;
            const uint32_t s = normalize(srcRaw);
            const uint32_t d = normalize(dstRaw);
            if(s == d) continue;
            auto key = (s < d) ? std::make_pair(s, d) : std::make_pair(d, s);
            pairEdges[key].push_back(e);
        }

        uint64_t singleEdgeRemovedCount = 0;
        for(const auto& [pair, edges] : pairEdges) {
            if(edges.size() != 1) continue;
            const uint32_t wA = pair.first;
            const uint32_t wB = pair.second;

            // Skip if this is an endpoint connection for either window.
            const auto& windowA = anchorWindows[wA];
            const auto& windowB = anchorWindows[wB];
            if(windowA.backbonePreviousWindow == wB ||
               windowA.backboneNextWindow == wB ||
               windowB.backbonePreviousWindow == wA ||
               windowB.backboneNextWindow == wA) {
                continue;
            }

            disableEdge(edges[0]);
            ++singleEdgeRemovedCount;
        }

        if(singleEdgeRemovedCount > 0) {
            cout << "Single-edge filter: removed " << singleEdgeRemovedCount
                 << " single-point internal connections." << endl;
        }
    };

    // Case 1: Bypass detour through another window.
    // Walking window w's backbone, if we find an incoming edge from
    // window X at anchor a_i and an outgoing edge to the same window X
    // at a later anchor a_j, then X's path detours through w's backbone.
    // Create a bypass edge in X (connecting X's anchors on either side
    // of the detour) and remove the inter-window edges.
    //
    // Example:
    //   Before: Window X: ... → x1          x2 → ...
    //                           |           ↑
    //           Window w: a0 → a_i → ... → a_j → a4
    //
    //   After:  Window X: ... → x1 ------→ x2 → ...  (bypass edge)
    //           Window w: a0 → a_i → ... → a_j → a4  (backbone intact)
    //           Inter-window edges x1→a_i and a_j→x2 removed.
    //
    auto runBypassDetourFilter = [&]() {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        uint64_t detoursFixed = 0;

        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& window = anchorWindows[w];
            const auto& positions = window.filteredBackbonePositions;
            if(positions.size() < 2) continue;
            const auto journey = journeys[window.backboneOrientedReadId];

            // Walk backbone anchors and collect inter-window connections.
            // For each backbone anchor, record incoming and outgoing
            // inter-window edges grouped by neighbor window.
            struct InterWindowEdge {
                edge_descriptor e;
                uint64_t otherAnchor; // anchor in the other window
                uint32_t bbIdx;       // index in backbone
            };
            // incoming[neighborWindow] = edges coming INTO w from that window
            // outgoing[neighborWindow] = edges going OUT of w to that window
            std::map<uint32_t, std::vector<InterWindowEdge>> incoming;
            std::map<uint32_t, std::vector<InterWindowEdge>> outgoing;

            for(uint32_t pi = 0; pi < positions.size(); pi++) {
                const uint64_t aid = uint64_t(journey[positions[pi]]);
                if(aid >= anchorCount) continue;

                // Check incoming edges (from other windows into this anchor).
                auto ie = boost::in_edges(aid, anchorGraph);
                for(auto it = ie.first; it != ie.second; ++it) {
                    if(!anchorGraph[*it].useForAssembly) continue;
                    const uint64_t src = uint64_t(boost::source(*it, anchorGraph));
                    if(src >= anchorCount) continue;
                    const uint32_t srcRaw = anchorToWindow[src];
                    if(srcRaw == noWindow) continue;
                    const uint32_t srcW = normalize(srcRaw);
                    if(srcW == w) continue; // intra-window, skip
                    incoming[srcW].push_back({*it, src, pi});
                }

                // Check outgoing edges (from this anchor to other windows).
                auto oe = boost::out_edges(aid, anchorGraph);
                for(auto it = oe.first; it != oe.second; ++it) {
                    if(!anchorGraph[*it].useForAssembly) continue;
                    const uint64_t tgt = uint64_t(boost::target(*it, anchorGraph));
                    if(tgt >= anchorCount) continue;
                    const uint32_t tgtRaw = anchorToWindow[tgt];
                    if(tgtRaw == noWindow) continue;
                    const uint32_t tgtW = normalize(tgtRaw);
                    if(tgtW == w) continue;
                    outgoing[tgtW].push_back({*it, tgt, pi});
                }
            }

            // For each neighbor window X that has both incoming and outgoing
            // connections, find pairs where incoming is at an earlier backbone
            // position than outgoing (X enters w then leaves w).
            // Skip endpoint connections — X is a legitimate neighbor.
            const uint32_t prevW = window.backbonePreviousWindow;
            const uint32_t nextW = window.backboneNextWindow;
            for(const auto& [xWindow, inEdges] : incoming) {
                // Skip if X is an endpoint of w or w is an endpoint of X.
                if(xWindow == prevW || xWindow == nextW) continue;
                const auto& xWin = anchorWindows[xWindow];
                if(xWin.backbonePreviousWindow == w || xWin.backboneNextWindow == w) continue;
                auto outIt = outgoing.find(xWindow);
                if(outIt == outgoing.end()) continue;
                const auto& outEdges = outIt->second;

                // Find earliest incoming and latest outgoing.
                // incoming at a_i (earliest), outgoing at a_j (latest after a_i).
                for(const auto& inE : inEdges) {
                    for(const auto& outE : outEdges) {
                        if(outE.bbIdx <= inE.bbIdx) continue;

                        // X enters w at inE.bbIdx, leaves at outE.bbIdx.
                        // Create bypass edge in X: inE.otherAnchor → outE.otherAnchor
                        const Shasta2AnchorId bypassSrc(inE.otherAnchor);
                        const Shasta2AnchorId bypassDst(outE.otherAnchor);
                        const bool bypassCreated = addEdgeIfValid(bypassSrc, bypassDst);
                        if(!bypassCreated) continue;
                        // Also create RC mirror bypass edge.
                        const Shasta2AnchorId rcBypassSrc(outE.otherAnchor ^ 1ULL);
                        const Shasta2AnchorId rcBypassDst(inE.otherAnchor ^ 1ULL);
                        addEdgeIfValid(rcBypassSrc, rcBypassDst);

                        // Remove the inter-window edges.
                        disableEdge(inE.e);
                        disableEdge(outE.e);
                        ++detoursFixed;
                    }
                }
            }
        }

        if(detoursFixed > 0) {
            cout << "Bypass detour filter: fixed " << detoursFixed
                 << " detours through other windows." << endl;
        }
    };

    // Bubble popping filter (graph surgery): for each window, detect
    // bubbles where inter-window edges create alternate paths between
    // backbone anchors. When found, keep the inter-window path and
    // remove the intra-window backbone segment between the two connection
    // points. Intermediate backbone anchors become orphaned.
    //
    // Example:
    //   Before: a0 → a1 → a2 → a3 → a4  (backbone of w)
    //                 |              ↑
    //                 → [X] ------→ |    (inter-window path)
    //
    //   After:  a0 → a1             a3 → a4  (backbone of w, linearized)
    //                 |              ↑
    //                 → [X] ------→ |    (inter-window path kept)
    //           a2 orphaned (all edges disabled)
    //
    auto runBubblePopFilter = [&]() {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        uint64_t bubblesPopped = 0;
        uint64_t anchorsDisconnected = 0;
        const uint32_t maxWindowsTraversed = 3;

        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& window = anchorWindows[w];
            const auto& positions = window.filteredBackbonePositions;
            if(positions.size() < 3) continue;
            const auto journey = journeys[window.backboneOrientedReadId];

            // Build ordered backbone anchor list and index lookup.
            std::vector<uint64_t> bbAnchors;
            std::map<uint64_t, uint32_t> bbAnchorIndex;
            for(uint32_t pi = 0; pi < positions.size(); pi++) {
                const uint64_t aid = uint64_t(journey[positions[pi]]);
                bbAnchors.push_back(aid);
                bbAnchorIndex[aid] = pi;
            }
            std::set<uint64_t> bbAnchorSet(bbAnchors.begin(), bbAnchors.end());

            const uint32_t prevW = window.backbonePreviousWindow;
            const uint32_t nextW = window.backboneNextWindow;
            std::set<uint32_t> orphanedIndices;

            for(uint32_t startIdx = 0; startIdx < bbAnchors.size(); startIdx++) {
                if(orphanedIndices.count(startIdx)) continue;
                const uint64_t startAid = bbAnchors[startIdx];

                // Check for internal inter-window outgoing edges
                // (skip edges to endpoint windows).
                bool hasInternalInterWindow = false;
                auto oe = boost::out_edges(startAid, anchorGraph);
                for(auto it = oe.first; it != oe.second; ++it) {
                    if(!anchorGraph[*it].useForAssembly) continue;
                    const uint64_t tgt = uint64_t(boost::target(*it, anchorGraph));
                    if(tgt >= anchorCount) continue;
                    const uint32_t tgtRaw = anchorToWindow[tgt];
                    if(tgtRaw == noWindow) continue;
                    const uint32_t tgtW = normalize(tgtRaw);
                    if(tgtW == w) continue;
                    if(tgtW == prevW || tgtW == nextW) continue;
                    hasInternalInterWindow = true;
                    break;
                }
                if(!hasInternalInterWindow) continue;

                // BFS from startAid through inter-window edges.
                struct BfsEntry {
                    uint64_t aid;
                    uint32_t windowsCrossed;
                    uint32_t currentWindow;
                };
                std::set<uint64_t> visited;
                std::queue<BfsEntry> q;
                visited.insert(startAid);

                // Seed BFS with internal inter-window neighbors.
                auto seedBfs = [&](uint64_t aid) {
                    auto oe2 = boost::out_edges(aid, anchorGraph);
                    for(auto it = oe2.first; it != oe2.second; ++it) {
                        if(!anchorGraph[*it].useForAssembly) continue;
                        const uint64_t tgt = uint64_t(boost::target(*it, anchorGraph));
                        if(tgt >= anchorCount || visited.count(tgt)) continue;
                        const uint32_t tgtRaw = anchorToWindow[tgt];
                        if(tgtRaw == noWindow) continue;
                        const uint32_t tgtW = normalize(tgtRaw);
                        if(tgtW == w) continue;
                        if(tgtW == prevW || tgtW == nextW) continue;
                        visited.insert(tgt);
                        q.push({tgt, 1, tgtW});
                    }
                    auto ie = boost::in_edges(aid, anchorGraph);
                    for(auto it = ie.first; it != ie.second; ++it) {
                        if(!anchorGraph[*it].useForAssembly) continue;
                        const uint64_t src = uint64_t(boost::source(*it, anchorGraph));
                        if(src >= anchorCount || visited.count(src)) continue;
                        const uint32_t srcRaw = anchorToWindow[src];
                        if(srcRaw == noWindow) continue;
                        const uint32_t srcW = normalize(srcRaw);
                        if(srcW == w) continue;
                        if(srcW == prevW || srcW == nextW) continue;
                        visited.insert(src);
                        q.push({src, 1, srcW});
                    }
                };
                seedBfs(startAid);

                uint32_t bestSinkIdx = UINT32_MAX;

                while(!q.empty()) {
                    auto entry = q.front();
                    q.pop();

                    if(bbAnchorSet.count(entry.aid)) {
                        auto idxIt = bbAnchorIndex.find(entry.aid);
                        if(idxIt != bbAnchorIndex.end() &&
                           idxIt->second > startIdx &&
                           idxIt->second < bestSinkIdx) {
                            bestSinkIdx = idxIt->second;
                        }
                        continue;
                    }

                    if(entry.windowsCrossed >= maxWindowsTraversed) continue;

                    auto processNeighbor = [&](uint64_t nbrAid) {
                        if(nbrAid >= anchorCount || visited.count(nbrAid)) return;
                        const uint32_t nbrRaw = anchorToWindow[nbrAid];
                        const uint32_t nbrW = (nbrRaw != noWindow) ? normalize(nbrRaw) : noWindow;
                        if(nbrW == w && !bbAnchorSet.count(nbrAid)) return;
                        uint32_t newWC = entry.windowsCrossed;
                        if(nbrW != noWindow && nbrW != w && nbrW != entry.currentWindow) {
                            newWC++;
                        }
                        visited.insert(nbrAid);
                        q.push({nbrAid, newWC, nbrW});
                    };

                    auto oe3 = boost::out_edges(entry.aid, anchorGraph);
                    for(auto it = oe3.first; it != oe3.second; ++it) {
                        if(!anchorGraph[*it].useForAssembly) continue;
                        processNeighbor(uint64_t(boost::target(*it, anchorGraph)));
                    }
                    auto ie3 = boost::in_edges(entry.aid, anchorGraph);
                    for(auto it = ie3.first; it != ie3.second; ++it) {
                        if(!anchorGraph[*it].useForAssembly) continue;
                        processNeighbor(uint64_t(boost::source(*it, anchorGraph)));
                    }
                }

                if(bestSinkIdx == UINT32_MAX) continue;
                if(bestSinkIdx <= startIdx + 1) continue;

                for(uint32_t idx = startIdx + 1; idx < bestSinkIdx; idx++) {
                    if(orphanedIndices.count(idx)) continue;
                    const uint64_t aid = bbAnchors[idx];

                    auto disableIntraWindow = [&](uint64_t a) {
                        if(a >= anchorCount) return;
                        auto oe4 = boost::out_edges(a, anchorGraph);
                        for(auto it = oe4.first; it != oe4.second; ++it) {
                            if(!anchorGraph[*it].useForAssembly) continue;
                            const uint64_t tgt = uint64_t(boost::target(*it, anchorGraph));
                            if(tgt >= anchorCount) continue;
                            const uint32_t tgtRaw = anchorToWindow[tgt];
                            if(tgtRaw == noWindow) continue;
                            if(normalize(tgtRaw) == w) {
                                disableEdge(*it);
                            }
                        }
                        auto ie4 = boost::in_edges(a, anchorGraph);
                        for(auto it = ie4.first; it != ie4.second; ++it) {
                            if(!anchorGraph[*it].useForAssembly) continue;
                            const uint64_t src = uint64_t(boost::source(*it, anchorGraph));
                            if(src >= anchorCount) continue;
                            const uint32_t srcRaw = anchorToWindow[src];
                            if(srcRaw == noWindow) continue;
                            if(normalize(srcRaw) == w) {
                                disableEdge(*it);
                            }
                        }
                    };
                    disableIntraWindow(aid);

                    orphanedIndices.insert(idx);
                    ++anchorsDisconnected;
                }
                ++bubblesPopped;
            }
        }

        if(bubblesPopped > 0) {
            cout << "Bubble pop: popped " << bubblesPopped
                 << " bubbles, disconnected " << anchorsDisconnected
                 << " backbone anchors." << endl;
        }
    };

    // Shortcut filter: a window is a shortcut window if its backbone
    // read transitions between two windows (backbonePreviousWindow and
    // backboneNextWindow) that are connected to each other.
    // The window is redundant — it shortcuts a path that already exists.
    auto runShortcutFilter = [&]() {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        // Build per-window neighbor sets from active edges.
        std::map<uint32_t, std::set<uint32_t>> neighbors;
        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
            if(!anchorGraph[e].useForAssembly) continue;
            const uint64_t srcVal = uint64_t(source(e, anchorGraph));
            const uint64_t dstVal = uint64_t(target(e, anchorGraph));
            if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
            const uint32_t srcRaw = anchorToWindow[srcVal];
            const uint32_t dstRaw = anchorToWindow[dstVal];
            if(srcRaw == noWindow || dstRaw == noWindow) continue;
            const uint32_t s = normalize(srcRaw);
            const uint32_t d = normalize(dstRaw);
            if(s == d) continue;
            neighbors[s].insert(d);
            neighbors[d].insert(s);
        }

        uint64_t shortcutRemovedCount = 0;

        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& window = anchorWindows[w];
            const uint32_t prevW = window.backbonePreviousWindow;
            const uint32_t nextW = window.backboneNextWindow;

            // Both endpoints must exist and be different.
            if(prevW == noWindow || nextW == noWindow) continue;
            if(prevW == nextW) continue;

            // Check if prevW and nextW are connected to each other.
            const auto prevIt = neighbors.find(prevW);
            if(prevIt == neighbors.end() || !prevIt->second.count(nextW)) {
                // They are not connected — not a shortcut window.
                continue;
            }

            // prevW and nextW are connected. This is a shortcut window.
            // Disable only edges between w and prevW/nextW.
            // Edges between w and other windows are preserved.
            BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
                if(!anchorGraph[e].useForAssembly) continue;
                const uint64_t srcVal = uint64_t(source(e, anchorGraph));
                const uint64_t dstVal = uint64_t(target(e, anchorGraph));
                if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
                const uint32_t srcRaw = anchorToWindow[srcVal];
                const uint32_t dstRaw = anchorToWindow[dstVal];
                if(srcRaw == noWindow || dstRaw == noWindow) continue;
                const uint32_t s = normalize(srcRaw);
                const uint32_t d = normalize(dstRaw);
                if(s == d) continue;
                // Only remove edges between w and prevW or nextW.
                if((s == w && (d == prevW || d == nextW)) ||
                   (d == w && (s == prevW || s == nextW))) {
                    disableEdge(e);
                    ++shortcutRemovedCount;
                }
            }
        }

        if(shortcutRemovedCount > 0) {
            cout << "Shortcut filter: removed " << shortcutRemovedCount
                 << " redundant shortcut edges." << endl;
        }
    };

    // ========================================================================
    // trimBackbones: remove dangling backbone ends that extend beyond
    // the outermost inter-window connections.
    //
    // Each window's backbone is an ordered chain of anchors along a
    // single read (the backbone read). The chain always goes forward:
    //   position 0 → position 1 → ... → position N-1
    // Intra-window edges connect consecutive positions. Inter-window
    // edges connect backbone anchors to anchors in other windows.
    //
    // The backbone chain direction is fixed by construction:
    //   - Head (position 0 side): where reads enter the window.
    //     Inter-window edges arriving here are INCOMING (other → this).
    //   - Tail (position N-1 side): where reads leave the window.
    //     Inter-window edges departing here are OUTGOING (this → other).
    //
    // Trimming removes backbone anchors that are beyond the outermost
    // inter-window connections:
    //   - Head trim: removes anchors before the first one that has an
    //     incoming inter-window edge. If no anchor has an incoming
    //     inter-window edge, head trim is skipped entirely (the window
    //     has no predecessor connection, so nothing to trim toward).
    //   - Tail trim: removes anchors after the last one that has an
    //     outgoing inter-window edge. If no anchor has an outgoing
    //     inter-window edge, tail trim is skipped entirely.
    //
    // Example (positions 0..6, inter-window edges at positions 2 and 5):
    //
    //   A0 -- A1 -- A2 -- A3 -- A4 -- A5 -- A6
    //                ↑ incoming          ↓ outgoing
    //   [trimmed]    [kept .............]  [trimmed]
    //   headTrim=2                         tailTrim=1
    //
    // RC mirror handling:
    //   Each forward anchor aid has an RC mirror aid^1. The RC mirror
    //   belongs to the virtual RC window (w + windowCount). Inter-window
    //   edges may exist only on the RC anchor (discovered by RC-strand
    //   reads). We must check both aid and aid^1 when deciding whether
    //   to trim, because disableAllEdges(aid^1) will destroy the RC
    //   anchor's edges too.
    //
    //   Edge direction flips on the RC mirror: the RC of edge A→B is
    //   B^1→A^1. So:
    //     - An OUT-edge on aid^1 = an IN-edge on aid (forward orientation)
    //     - An IN-edge on aid^1 = an OUT-edge on aid (forward orientation)
    //
    // Deletion:
    //   Trimmed anchors have no inter-window edges by construction
    //   (otherwise the trim would have stopped before them). So
    //   disableAllEdges only destroys intra-window edges. We call it
    //   on both aid and aid^1 to symmetrically disconnect both the
    //   forward and RC anchors. disableEdge() ensures each disabled
    //   edge also disables its RC mirror edge (dst^1 → src^1).
    //
    // Safety checks:
    //   - If headTrim >= positions.size() (no anchor has an incoming
    //     inter-window edge), headTrim resets to 0 — no head trim.
    //     The window may be isolated or tail-only; handled elsewhere.
    //   - If headTrim + tailTrim >= positions.size(), tailTrim resets
    //     to 0 — prevents consuming all positions from both sides.
    //   - Windows with <= 1 backbone position are skipped.
    //
    // Called after each filter step to clean up backbone ends exposed
    // by edge removal. filteredBackbonePositions is NOT updated — the
    // trimmed anchors remain in the array but have all edges disabled.
    // On subsequent calls, they are re-examined and re-trimmed (no-op
    // since their edges are already disabled).
    // ========================================================================
    auto trimBackbones = [&]() {

        // Disable all active edges of a single anchor. Iterates both
        // out-edges and in-edges, calling disableEdge on each active
        // one. disableEdge sets useForAssembly=false on the edge AND
        // its RC mirror (dst^1 → src^1), so this also partially
        // disconnects the RC anchor aid^1.
        auto disableAllEdges = [&](uint64_t aid) {
            if(aid >= anchorCount) return;
            auto oe = boost::out_edges(aid, anchorGraph);
            for(auto it = oe.first; it != oe.second; ++it) {
                if(anchorGraph[*it].useForAssembly)
                    disableEdge(*it);
            }
            auto ie = boost::in_edges(aid, anchorGraph);
            for(auto it = ie.first; it != ie.second; ++it) {
                if(anchorGraph[*it].useForAssembly)
                    disableEdge(*it);
            }
        };

        // Normalize a window ID: RC windows (>= windowCount) map to
        // their forward counterpart. Forward windows map to themselves.
        auto normalizeW = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        // Check if anchor `a` has any active incoming edge from a
        // different normalized window. An incoming edge means some
        // other window's anchor has a directed edge INTO this anchor.
        auto anchorHasInterWindowInEdge = [&](uint64_t a) -> bool {
            if(a >= anchorCount) return false;
            const uint32_t aWin = anchorToWindow[a];
            if(aWin == noWindow) return false;
            const uint32_t aNorm = normalizeW(aWin);
            auto ie = boost::in_edges(a, anchorGraph);
            for(auto it = ie.first; it != ie.second; ++it) {
                if(!anchorGraph[*it].useForAssembly) continue;
                const uint64_t src = uint64_t(boost::source(*it, anchorGraph));
                if(src < anchorCount) {
                    const uint32_t srcWin = anchorToWindow[src];
                    if(srcWin != noWindow && normalizeW(srcWin) != aNorm)
                        return true;
                }
            }
            return false;
        };

        // Check if anchor `a` has any active outgoing edge to a
        // different normalized window. An outgoing edge means this
        // anchor has a directed edge TO some other window's anchor.
        auto anchorHasInterWindowOutEdge = [&](uint64_t a) -> bool {
            if(a >= anchorCount) return false;
            const uint32_t aWin = anchorToWindow[a];
            if(aWin == noWindow) return false;
            const uint32_t aNorm = normalizeW(aWin);
            auto oe = boost::out_edges(a, anchorGraph);
            for(auto it = oe.first; it != oe.second; ++it) {
                if(!anchorGraph[*it].useForAssembly) continue;
                const uint64_t tgt = uint64_t(boost::target(*it, anchorGraph));
                if(tgt < anchorCount) {
                    const uint32_t tgtWin = anchorToWindow[tgt];
                    if(tgtWin != noWindow && normalizeW(tgtWin) != aNorm)
                        return true;
                }
            }
            return false;
        };

        // Check if anchor `aid` (or its RC mirror) has an incoming
        // inter-window edge in forward orientation.
        //
        // For the forward anchor aid: check in-edges directly.
        // For the RC anchor aid^1: check OUT-edges, because the RC
        // of edge X → aid is aid^1 → X^1 (direction reverses).
        // So an out-edge on aid^1 to another window = an in-edge
        // on aid from that window in forward orientation.
        auto hasIncomingInterWindowEdge = [&](uint64_t aid) -> bool {
            if(anchorHasInterWindowInEdge(aid)) return true;
            const uint64_t rcAid = aid ^ 1ULL;
            return (rcAid < anchorCount && anchorHasInterWindowOutEdge(rcAid));
        };

        // Check if anchor `aid` (or its RC mirror) has an outgoing
        // inter-window edge in forward orientation.
        //
        // For the forward anchor aid: check out-edges directly.
        // For the RC anchor aid^1: check IN-edges, because the RC
        // of edge aid → Y is Y^1 → aid^1 (direction reverses).
        // So an in-edge on aid^1 from another window = an out-edge
        // on aid to that window in forward orientation.
        auto hasOutgoingInterWindowEdge = [&](uint64_t aid) -> bool {
            if(anchorHasInterWindowOutEdge(aid)) return true;
            const uint64_t rcAid = aid ^ 1ULL;
            return (rcAid < anchorCount && anchorHasInterWindowInEdge(rcAid));
        };

        uint64_t trimmedVertexCount = 0;
        uint64_t trimmedWindowCount = 0;

        // Iterate forward windows only. RC windows are virtual (no
        // AnchorWindow entries). RC anchors are handled via
        // disableAllEdges(aid^1) when their forward counterpart is
        // trimmed.
        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& window = anchorWindows[w];
            const auto& positions = window.filteredBackbonePositions;
            if(positions.size() <= 1) continue;
            const auto journey = journeys[window.backboneOrientedReadId];

            // --- Head trim ---
            // Walk from position 0 (head) inward. Count consecutive
            // anchors that have no incoming inter-window edge. Stop
            // at the first anchor that does — that's where the window's
            // predecessor connection lands.
            uint64_t headTrim = 0;
            for(uint64_t i = 0; i < positions.size(); i++) {
                const uint64_t aid = uint64_t(journey[positions[i]]);
                if(hasIncomingInterWindowEdge(aid)) break;
                ++headTrim;
            }
            // Safety: if no anchor has an incoming inter-window edge,
            // headTrim == positions.size(). Reset to 0 — the window
            // has no predecessor, so don't trim the head. (The window
            // may be a telomere start or isolated; handled elsewhere.)
            if(headTrim >= positions.size()) headTrim = 0;

            // --- Tail trim ---
            // Walk from position N-1 (tail) inward toward headTrim.
            // Count consecutive anchors that have no outgoing
            // inter-window edge. Stop at the first anchor that does —
            // that's where the window's successor connection departs.
            uint64_t tailTrim = 0;
            for(int64_t i = int64_t(positions.size()) - 1; i >= int64_t(headTrim); i--) {
                const uint64_t aid = uint64_t(journey[positions[uint64_t(i)]]);
                if(hasOutgoingInterWindowEdge(aid)) break;
                ++tailTrim;
            }
            // Safety: if headTrim + tailTrim >= positions.size(), the
            // head and tail trims overlap (no outgoing edge found in
            // the remaining range). Reset tailTrim to 0 to avoid
            // consuming all positions.
            if(headTrim + tailTrim >= positions.size()) tailTrim = 0;

            // Nothing to trim for this window.
            if(headTrim == 0 && tailTrim == 0) continue;
            ++trimmedWindowCount;

            // --- Disable edges of trimmed head anchors ---
            // These anchors are before the first incoming inter-window
            // edge, so they have no inter-window edges (incoming or
            // outgoing). All their edges are intra-window.
            // disableAllEdges(aid): disables all edges of the forward
            //   anchor. Each disableEdge call also disables the RC
            //   mirror edge on aid^1.
            // disableAllEdges(aid^1): disables any remaining edges on
            //   the RC anchor that weren't RC mirrors of aid's edges
            //   (e.g., RC-only intra-window edges created when
            //   addEdgeIfValid succeeded for RC but not forward).
            for(uint64_t i = 0; i < headTrim; i++) {
                const uint64_t aid = uint64_t(journey[positions[i]]);
                disableAllEdges(aid);
                disableAllEdges(aid ^ 1ULL);
                ++trimmedVertexCount;
            }

            // --- Disable edges of trimmed tail anchors ---
            // Same logic: these anchors are after the last outgoing
            // inter-window edge, so they have no inter-window edges.
            for(uint64_t i = positions.size() - tailTrim; i < positions.size(); i++) {
                const uint64_t aid = uint64_t(journey[positions[i]]);
                disableAllEdges(aid);
                disableAllEdges(aid ^ 1ULL);
                ++trimmedVertexCount;
            }
        }

        cout << "Trim backbones: " << trimmedWindowCount << " windows trimmed, "
             << trimmedVertexCount << " vertices trimmed." << endl;
        return trimmedVertexCount;
    };

    // Dangling window cleanup: remove inter-window edges for windows
    // with edges on only one side (only incoming or only outgoing),
    // provided removal won't make any neighbor dangling.
    // Defined as a lambda so it can be called after each filter pass.
    // Remove dangling windows: windows with inter-window edges on only
    // one side (only incoming or only outgoing). A dangling window is
    // removed only if every neighbor has OTHER connections on the SAME
    // side where the dangling window connects, and at least one of those
    // other connections is from a non-dangling window (has both in and out).
    //
    // Test cases:
    //  1. Simple noise tip: W_tip -> N <- W_prev(non-dangling), N -> W_next.
    //     N has other incoming from W_prev (non-dangling) -> remove W_tip.
    //  2. Telomere (sole source): W_telo -> N -> W_next.
    //     N has no other incoming -> preserve W_telo.
    //  3. Two dangling tips on same neighbor: W_tip1 -> N <- W_tip2, N -> W_next.
    //     Other source into N is W_tip2 (dangling) -> not non-dangling -> preserve both.
    //  4. Dangling + non-dangling on same neighbor: W_tip -> N <- W_prev(non-dangling).
    //     W_prev is non-dangling -> remove W_tip.
    //  5. Dangling tip with multiple neighbors: W_tip -> N1, W_tip -> N2.
    //     Both N1 and N2 must pass the safety check. If either fails -> preserve.
    //  6. Incoming-only dangling: ... -> W_prev -> N -> W_tip.
    //     Same side = N's outgoing. If N only outputs to W_tip -> preserve (telomere end).
    //     If N also outputs to W_next (non-dangling) -> remove W_tip.
    //  7. Chain of dangling windows: W1 -> N <- W2 -> M <- W3 (all dangling).
    //     Each tip's neighbor's other sources are also dangling -> preserve all.
    //  8. RC mirror symmetry: forward and RC raw windows are always dangling
    //     in complementary directions. Processing either one via disableEdge
    //     correctly removes both sides.
    //  9. Cascading removal: after removing W_tip, neighbor N may lose edges,
    //     making other windows dangling. The while(changed) loop re-checks.
    //     Count maps are updated by disableAndUpdateCounts before re-check.
    // 10. Empty edgeIndices: if a window appears dangling by raw counts but
    //     has no active edges (stale counts from backbone cleanup), skip it.
    // Remove dead-end spurs: windows that connect to exactly one other
    // window, and that connection is internal (not an endpoint connection
    // for either window). The spur adds nothing to the graph.
    auto removeDeadEndSpurs = [&]() {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        // Build per-window neighbor sets from active edges.
        std::map<uint32_t, std::set<uint32_t>> neighbors;
        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
            if(!anchorGraph[e].useForAssembly) continue;
            const uint64_t srcVal = uint64_t(source(e, anchorGraph));
            const uint64_t dstVal = uint64_t(target(e, anchorGraph));
            if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
            const uint32_t srcRaw = anchorToWindow[srcVal];
            const uint32_t dstRaw = anchorToWindow[dstVal];
            if(srcRaw == noWindow || dstRaw == noWindow) continue;
            const uint32_t s = normalize(srcRaw);
            const uint32_t d = normalize(dstRaw);
            if(s == d) continue;
            neighbors[s].insert(d);
            neighbors[d].insert(s);
        }

        uint64_t spursRemoved = 0;

        for(uint32_t w = 0; w < windowCount; w++) {
            auto it = neighbors.find(w);
            if(it == neighbors.end()) continue;
            if(it->second.size() != 1) continue;

            const uint32_t neighbor = *it->second.begin();
            const auto& window = anchorWindows[w];
            const auto& nbrWindow = anchorWindows[neighbor];

            // Check that the connection is internal for the neighbor.
            // We don't check the spur window's endpoints — a spur with
            // one neighbor always has that neighbor as its endpoint,
            // so the check would always skip it.
            if(nbrWindow.backbonePreviousWindow == w ||
               nbrWindow.backboneNextWindow == w) continue;

            // Dead-end spur. Disable all inter-window edges of w.
            BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
                if(!anchorGraph[e].useForAssembly) continue;
                const uint64_t srcVal = uint64_t(source(e, anchorGraph));
                const uint64_t dstVal = uint64_t(target(e, anchorGraph));
                if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
                const uint32_t srcRaw = anchorToWindow[srcVal];
                const uint32_t dstRaw = anchorToWindow[dstVal];
                if(srcRaw == noWindow || dstRaw == noWindow) continue;
                const uint32_t s = normalize(srcRaw);
                const uint32_t d = normalize(dstRaw);
                if(s == d) continue;
                if(s == w || d == w) {
                    disableEdge(e);
                }
            }
            ++spursRemoved;
        }

        if(spursRemoved > 0) {
            cout << "Dead-end spur removal: removed " << spursRemoved
                 << " spur windows." << endl;
        }
    };

    auto removeDanglingWindows = [&](const string& label) {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };
        // Build per-edge info using raw window IDs.
        struct InterWindowEdge {
            edge_descriptor e;
            uint32_t srcWin;  // raw
            uint32_t dstWin;  // raw
        };
        std::vector<InterWindowEdge> interWindowEdges;
        // Raw counts for dangling detection.
        std::map<uint32_t, uint64_t> rawInCount, rawOutCount;
        // Normalized counts for safety check.
        std::map<uint32_t, uint64_t> normInCount, normOutCount;

        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
            if(!anchorGraph[e].useForAssembly) continue;
            const uint64_t srcVal = uint64_t(source(e, anchorGraph));
            const uint64_t dstVal = uint64_t(target(e, anchorGraph));
            if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
            const uint32_t srcWinRaw = anchorToWindow[srcVal];
            const uint32_t dstWinRaw = anchorToWindow[dstVal];
            if(srcWinRaw == noWindow || dstWinRaw == noWindow) continue;
            if(srcWinRaw == dstWinRaw) continue;
            const uint32_t srcNorm = normalize(srcWinRaw);
            const uint32_t dstNorm = normalize(dstWinRaw);
            if(srcNorm == dstNorm) continue;  // skip fwd<->RC of same window
            interWindowEdges.push_back({e, srcWinRaw, dstWinRaw});
            rawOutCount[srcWinRaw]++;
            rawInCount[dstWinRaw]++;
            normOutCount[srcNorm]++;
            normInCount[dstNorm]++;
        }

        // Collect all raw window IDs.
        std::set<uint32_t> allRawWindows;
        for(const auto& [w, c] : rawInCount) if(c > 0) allRawWindows.insert(w);
        for(const auto& [w, c] : rawOutCount) if(c > 0) allRawWindows.insert(w);

        uint64_t danglingRemovedCount = 0;
        uint64_t danglingWindowCount = 0;

        // Disable an edge + RC mirror and update local count maps.
        auto disableAndUpdateCounts = [&](edge_descriptor e) {
            if(!anchorGraph[e].useForAssembly) return;

            // Decrement counts for this edge and its RC mirror before disabling.
            auto decrementCounts = [&](edge_descriptor ed) {
                const uint64_t sv = uint64_t(source(ed, anchorGraph));
                const uint64_t dv = uint64_t(target(ed, anchorGraph));
                if(sv < anchorCount && dv < anchorCount) {
                    const uint32_t sw = anchorToWindow[sv];
                    const uint32_t dw = anchorToWindow[dv];
                    if(sw != noWindow && dw != noWindow &&
                       sw != dw && normalize(sw) != normalize(dw)) {
                        rawOutCount[sw]--;
                        rawInCount[dw]--;
                        normOutCount[normalize(sw)]--;
                        normInCount[normalize(dw)]--;
                        ++danglingRemovedCount;
                    }
                }
            };

            // Find RC mirror before disabling (disableEdge will disable both).
            const uint64_t srcVal = uint64_t(source(e, anchorGraph));
            const uint64_t dstVal = uint64_t(target(e, anchorGraph));
            const uint64_t rcSrc = dstVal ^ 1ULL;
            const uint64_t rcDst = srcVal ^ 1ULL;
            edge_descriptor mirrorE;
            bool hasMirror = false;
            if(rcSrc < anchorCount && rcDst < anchorCount) {
                auto [eit, exists] = boost::edge(rcSrc, rcDst, anchorGraph);
                if(exists && anchorGraph[eit].useForAssembly) {
                    hasMirror = true;
                    mirrorE = eit;
                }
            }

            decrementCounts(e);
            if(hasMirror) decrementCounts(mirrorE);
            disableEdge(e);  // disables both edge and RC mirror
        };

        bool changed = true;
        while(changed) {
            changed = false;
            for(const uint32_t w : allRawWindows) {
                const bool hasIn = (rawInCount.count(w) && rawInCount[w] > 0);
                const bool hasOut = (rawOutCount.count(w) && rawOutCount[w] > 0);
                if(hasIn == hasOut) continue;
                ++danglingWindowCount;

                // Collect this window's active inter-window edges.
                std::vector<size_t> edgeIndices;
                for(size_t i = 0; i < interWindowEdges.size(); i++) {
                    const auto& iwe = interWindowEdges[i];
                    if(!anchorGraph[iwe.e].useForAssembly) continue;
                    if(iwe.srcWin == w || iwe.dstWin == w) {
                        edgeIndices.push_back(i);
                    }
                }

                if(edgeIndices.empty()) continue;

                // Safety: only remove if every neighbor has OTHER connections
                // on the SAME side where w connects, and at least one of
                // those other connections is from a non-dangling window.
                // This preserves telomeres and avoids removing a tip whose
                // only "backup" is another dangling window.
                const uint32_t wNorm = normalize(w);
                bool safeToRemove = true;

                // Collect unique normalized neighbors to avoid redundant checks.
                std::set<uint32_t> checkedNeighbors;

                for(const size_t idx : edgeIndices) {
                    const auto& iwe = interWindowEdges[idx];
                    if(!anchorGraph[iwe.e].useForAssembly) continue;
                    const uint32_t nRaw = (iwe.srcWin == w) ? iwe.dstWin : iwe.srcWin;
                    const uint32_t nNorm = normalize(nRaw);

                    if(!checkedNeighbors.insert(nNorm).second) continue;

                    const bool wIsSource = (iwe.srcWin == w);

                    // Find other windows on the same side of N (excluding w).
                    // At least one must be non-dangling (has edges on both sides).
                    bool hasNonDanglingOther = false;
                    for(const auto& iwe2 : interWindowEdges) {
                        if(!anchorGraph[iwe2.e].useForAssembly) continue;
                        const uint32_t src2Norm = normalize(iwe2.srcWin);
                        const uint32_t dst2Norm = normalize(iwe2.dstWin);

                        uint32_t otherNorm;
                        if(wIsSource) {
                            // w -> N: same side = N's incoming. Find other sources.
                            if(dst2Norm != nNorm) continue;
                            if(src2Norm == wNorm) continue; // back to w
                            otherNorm = src2Norm;
                        } else {
                            // N -> w: same side = N's outgoing. Find other destinations.
                            if(src2Norm != nNorm) continue;
                            if(dst2Norm == wNorm) continue; // back to w
                            otherNorm = dst2Norm;
                        }

                        // Check that this other window is not dangling
                        // (has edges on both sides).
                        const uint64_t oIn = normInCount.count(otherNorm) ? normInCount[otherNorm] : 0;
                        const uint64_t oOut = normOutCount.count(otherNorm) ? normOutCount[otherNorm] : 0;
                        if(oIn > 0 && oOut > 0) {
                            hasNonDanglingOther = true;
                            break;
                        }
                    }

                    if(!hasNonDanglingOther) {
                        safeToRemove = false;
                        break;
                    }
                }

                if(!safeToRemove) continue;

                for(const size_t idx : edgeIndices) {
                    const auto& iwe = interWindowEdges[idx];
                    disableAndUpdateCounts(iwe.e);
                }

                // Also disable all edges (including intra-window) of
                // every anchor in this window, so the backbone doesn't
                // remain as a disconnected chain in the output.
                // Use normalized window to find the forward backbone,
                // then handle both forward and RC anchors.
                const uint32_t normW = normalize(w);
                if(normW < anchorWindows.size()) {
                    const auto& dangleWindow = anchorWindows[normW];
                    const auto& danglePositions = dangleWindow.filteredBackbonePositions;
                    const auto dangleJourney = journeys[dangleWindow.backboneOrientedReadId];
                    for(const uint32_t pos : danglePositions) {
                        const uint64_t aid = uint64_t(dangleJourney[pos]);
                        // Disable all edges of forward anchor.
                        auto oe = boost::out_edges(aid, anchorGraph);
                        for(auto it = oe.first; it != oe.second; ++it) {
                            if(anchorGraph[*it].useForAssembly)
                                disableEdge(*it);
                        }
                        auto ie = boost::in_edges(aid, anchorGraph);
                        for(auto it = ie.first; it != ie.second; ++it) {
                            if(anchorGraph[*it].useForAssembly)
                                disableEdge(*it);
                        }
                        // RC anchor edges are handled by disableEdge.
                    }
                }

                changed = true;
            }
        }

        cout << "Dangling window cleanup (" << label << "): found "
             << danglingWindowCount << " dangling windows, removed "
             << danglingRemovedCount << " inter-window edges." << endl;
        return danglingRemovedCount;
    };

    // Wrapper that runs dangling cleanup iteratively until stable.
    auto removeDanglingWindowsIterative = [&](const string& label) {
        uint64_t totalRemoved = 0;
        uint64_t iteration = 0;
        while(true) {
            const string iterLabel = label + " pass " + to_string(++iteration);
            const uint64_t removed = removeDanglingWindows(iterLabel);
            if(removed == 0) break;
            totalRemoved += removed;
        }
    };

    // Remove small windows: windows with few backbone anchors (≤ maxSmallWindow)
    // that are not needed for connectivity. Uses the same neighbor safety check
    // as dangling removal — only remove if every neighbor has other connections
    // on the same side from a non-dangling window.
    auto removeSmallWindows = [&](uint64_t maxSmallWindow) {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        // Build normalized in/out counts from active inter-window edges.
        std::map<uint32_t, uint64_t> normInCount, normOutCount;
        // Collect inter-window edges for neighbor lookup.
        struct InterWindowEdge {
            edge_descriptor e;
            uint32_t srcNorm;
            uint32_t dstNorm;
        };
        std::vector<InterWindowEdge> interWindowEdges;

        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
            if(!anchorGraph[e].useForAssembly) continue;
            const uint64_t srcVal = uint64_t(source(e, anchorGraph));
            const uint64_t dstVal = uint64_t(target(e, anchorGraph));
            if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
            const uint32_t srcW = anchorToWindow[srcVal];
            const uint32_t dstW = anchorToWindow[dstVal];
            if(srcW == noWindow || dstW == noWindow) continue;
            const uint32_t srcNorm = normalize(srcW);
            const uint32_t dstNorm = normalize(dstW);
            if(srcNorm == dstNorm) continue;
            interWindowEdges.push_back({e, srcNorm, dstNorm});
            normOutCount[srcNorm]++;
            normInCount[dstNorm]++;
        }

        uint64_t removedCount = 0;

        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& window = anchorWindows[w];
            const auto& positions = window.filteredBackbonePositions;
            if(positions.size() > maxSmallWindow) continue;

            // Check if this window has any active inter-window edges.
            bool hasActiveEdges = false;
            for(const auto& iwe : interWindowEdges) {
                if(!anchorGraph[iwe.e].useForAssembly) continue;
                if(iwe.srcNorm == w || iwe.dstNorm == w) {
                    hasActiveEdges = true;
                    break;
                }
            }
            if(!hasActiveEdges) continue;

            // Safety: for each neighbor, check that it has other connections
            // on the same side from a non-dangling window.
            bool safeToRemove = true;
            std::set<uint32_t> checkedNeighbors;

            for(const auto& iwe : interWindowEdges) {
                if(!anchorGraph[iwe.e].useForAssembly) continue;
                uint32_t nNorm;
                bool wIsSource;
                if(iwe.srcNorm == w) {
                    nNorm = iwe.dstNorm;
                    wIsSource = true;
                } else if(iwe.dstNorm == w) {
                    nNorm = iwe.srcNorm;
                    wIsSource = false;
                } else {
                    continue;
                }

                if(!checkedNeighbors.insert(nNorm).second) continue;

                // Find other windows on the same side of N (excluding w).
                // At least one must be non-dangling.
                bool hasNonDanglingOther = false;
                for(const auto& iwe2 : interWindowEdges) {
                    if(!anchorGraph[iwe2.e].useForAssembly) continue;
                    uint32_t otherNorm;
                    if(wIsSource) {
                        // w -> N: same side = N's incoming.
                        if(iwe2.dstNorm != nNorm) continue;
                        if(iwe2.srcNorm == w) continue;
                        otherNorm = iwe2.srcNorm;
                    } else {
                        // N -> w: same side = N's outgoing.
                        if(iwe2.srcNorm != nNorm) continue;
                        if(iwe2.dstNorm == w) continue;
                        otherNorm = iwe2.dstNorm;
                    }

                    const uint64_t oIn = normInCount.count(otherNorm) ? normInCount[otherNorm] : 0;
                    const uint64_t oOut = normOutCount.count(otherNorm) ? normOutCount[otherNorm] : 0;
                    if(oIn > 0 && oOut > 0) {
                        hasNonDanglingOther = true;
                        break;
                    }
                }

                if(!hasNonDanglingOther) {
                    safeToRemove = false;
                    break;
                }
            }

            if(!safeToRemove) continue;

            // Disable all inter-window edges of this window.
            for(const auto& iwe : interWindowEdges) {
                if(!anchorGraph[iwe.e].useForAssembly) continue;
                if(iwe.srcNorm == w || iwe.dstNorm == w) {
                    // Update counts before disabling.
                    normOutCount[iwe.srcNorm]--;
                    normInCount[iwe.dstNorm]--;
                    disableEdge(iwe.e);
                }
            }

            // Disable all backbone edges.
            const auto journey = journeys[window.backboneOrientedReadId];
            for(const uint32_t pos : positions) {
                const uint64_t aid = uint64_t(journey[pos]);
                auto oe = boost::out_edges(aid, anchorGraph);
                for(auto it = oe.first; it != oe.second; ++it) {
                    if(anchorGraph[*it].useForAssembly)
                        disableEdge(*it);
                }
                auto ie = boost::in_edges(aid, anchorGraph);
                for(auto it = ie.first; it != ie.second; ++it) {
                    if(anchorGraph[*it].useForAssembly)
                        disableEdge(*it);
                }
            }

            ++removedCount;
        }

        if(removedCount > 0) {
            cout << "Small window removal (maxSize=" << maxSmallWindow
                 << "): removed " << removedCount << " windows." << endl;
        }
    };

    // Remove windows that have no active inter-window edges.
    // Their backbone is left floating after other filters removed
    // all their inter-window connections.
    auto removeIsolatedWindows = [&]() {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        // Find which normalized windows still have active inter-window edges.
        std::set<uint32_t> connectedWindows;
        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
            if(!anchorGraph[e].useForAssembly) continue;
            const uint64_t srcVal = uint64_t(source(e, anchorGraph));
            const uint64_t dstVal = uint64_t(target(e, anchorGraph));
            if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
            const uint32_t srcRaw = anchorToWindow[srcVal];
            const uint32_t dstRaw = anchorToWindow[dstVal];
            if(srcRaw == noWindow || dstRaw == noWindow) continue;
            const uint32_t srcNorm = normalize(srcRaw);
            const uint32_t dstNorm = normalize(dstRaw);
            if(srcNorm != dstNorm) {
                connectedWindows.insert(srcNorm);
                connectedWindows.insert(dstNorm);
            }
        }

        uint64_t removedCount = 0;
        for(uint32_t w = 0; w < windowCount; w++) {
            if(connectedWindows.count(w)) continue;
            const auto& window = anchorWindows[w];
            const auto& positions = window.filteredBackbonePositions;
            if(positions.empty()) continue;
            const auto journey = journeys[window.backboneOrientedReadId];

            bool hadActive = false;
            for(const uint32_t pos : positions) {
                // Disable all edges of both forward and RC anchors.
                const uint64_t fwdAid = uint64_t(journey[pos]);
                const uint64_t rcAid = fwdAid ^ 1ULL;
                for(const uint64_t aid : {fwdAid, rcAid}) {
                    if(aid >= anchorCount) continue;
                    auto oe = boost::out_edges(aid, anchorGraph);
                    for(auto it = oe.first; it != oe.second; ++it) {
                        if(anchorGraph[*it].useForAssembly) {
                            disableEdge(*it);
                            hadActive = true;
                        }
                    }
                    auto ie = boost::in_edges(aid, anchorGraph);
                    for(auto it = ie.first; it != ie.second; ++it) {
                        if(anchorGraph[*it].useForAssembly) {
                            disableEdge(*it);
                            hadActive = true;
                        }
                    }
                }
            }
            if(hadActive) ++removedCount;
        }

        if(removedCount > 0) {
            cout << "Isolated window removal: removed " << removedCount
                 << " windows with no inter-window edges." << endl;
        }
    };

    // Cross-window filter: a window is a cross-window if its backbone
    // read transitions between two windows (backbonePreviousWindow and
    // backboneNextWindow) that are not connected to each other.
    // This means the backbone read bridges unrelated regions.
    auto runCrossWindowFilter = [&]() {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        // Build per-window neighbor sets from active edges.
        std::map<uint32_t, std::set<uint32_t>> neighbors;
        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
            if(!anchorGraph[e].useForAssembly) continue;
            const uint64_t srcVal = uint64_t(source(e, anchorGraph));
            const uint64_t dstVal = uint64_t(target(e, anchorGraph));
            if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
            const uint32_t srcRaw = anchorToWindow[srcVal];
            const uint32_t dstRaw = anchorToWindow[dstVal];
            if(srcRaw == noWindow || dstRaw == noWindow) continue;
            const uint32_t s = normalize(srcRaw);
            const uint32_t d = normalize(dstRaw);
            if(s == d) continue;
            neighbors[s].insert(d);
            neighbors[d].insert(s);
        }

        uint64_t crossWindowRemovedCount = 0;

        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& window = anchorWindows[w];
            const uint32_t prevW = window.backbonePreviousWindow;
            const uint32_t nextW = window.backboneNextWindow;

            // Both endpoints must exist and be different.
            if(prevW == noWindow || nextW == noWindow) continue;
            if(prevW == nextW) continue;

            // Check if prevW and nextW are connected to each other.
            const auto prevIt = neighbors.find(prevW);
            if(prevIt != neighbors.end() && prevIt->second.count(nextW)) {
                // They are connected — not a cross-window.
                continue;
            }

            // prevW and nextW are not connected. This is a cross-window.
            // Disable all its inter-window edges.
            BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
                if(!anchorGraph[e].useForAssembly) continue;
                const uint64_t srcVal = uint64_t(source(e, anchorGraph));
                const uint64_t dstVal = uint64_t(target(e, anchorGraph));
                if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
                const uint32_t srcRaw = anchorToWindow[srcVal];
                const uint32_t dstRaw = anchorToWindow[dstVal];
                if(srcRaw == noWindow || dstRaw == noWindow) continue;
                const uint32_t s = normalize(srcRaw);
                const uint32_t d = normalize(dstRaw);
                if(s == d) continue;
                if(s == w || d == w) {
                    disableEdge(e);
                    ++crossWindowRemovedCount;
                }
            }
        }

        if(crossWindowRemovedCount > 0) {
            cout << "Cross-window filter: removed " << crossWindowRemovedCount
                 << " spurious bridge edges." << endl;
        }
    };

    // ========================================================================
    // Filter pipeline: trimBackbones runs after each filter to clean up
    // backbone ends exposed by edge removal. recomputeBackboneEndpoints
    // runs before filters that use backbonePreviousWindow/backboneNextWindow.
    // ========================================================================
    trimBackbones();
    recomputeBackboneEndpoints();
    removeDeadEndSpurs();           // Remove single-neighbor internal spurs
    trimBackbones();
    recomputeBackboneEndpoints();
    runSingleEdgeFilter();          // Case 2: remove single-point connections
    trimBackbones();
    recomputeBackboneEndpoints();
    runBypassDetourFilter();        // Case 1: bypass detours through other windows
    trimBackbones();
    recomputeBackboneEndpoints();
    runBubblePopFilter();           // Case 3: pop bubbles returning to same window
    trimBackbones();
    recomputeBackboneEndpoints();
    runShortcutFilter();
    trimBackbones();
    recomputeBackboneEndpoints();
    runCrossWindowFilter();
    trimBackbones();
    removeIsolatedWindows();
    trimBackbones();
    removeSmallWindows(2);
    trimBackbones();
    removeDanglingWindowsIterative("post-filter");

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
    // Bypass edges from detangling — disabled for now.
#if 0
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
#endif

    // ========================================================================
    // Detangle Case 2: 2x2 tangle matrix for internal inter-window edges.
    // COMMENTED OUT pending validation on larger datasets.
    // ========================================================================
#if 0  // DISABLED for incremental filter testing
    // ========================================================================
    //
    // When two windows A and B are connected by an internal inter-window
    // edge (both endpoints are mid-backbone), check whether the connection
    // is real or spurious using a tangle matrix approach similar to
    // shasta2's Shasta2TangleMatrix1.
    //
    // For each window, collect reads from all backbone anchors before the
    // connection point ("entrance" side) and all anchors after it ("exit"
    // side). Each read's contribution is weighted by the number of anchors
    // (steps) it appears on in each side, following shasta2's per-read
    // matrix approach.
    //
    // Build a 2x2 matrix:
    //   A-entrance → A-exit (diagonal: A continues on its own)
    //   B-entrance → B-exit (diagonal: B continues on its own)
    //   A-entrance → B-exit (off-diagonal: cross from A to B)
    //   B-entrance → A-exit (off-diagonal: cross from B to A)
    //
    // If both diagonal entries are stronger than both off-diagonal
    // entries, the internal connection is false and the edge is removed.
    // ========================================================================
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

            // Check if an anchor or its RC mirror has any active inter-window edge.
            // Uses normalized window IDs to avoid false positives between
            // forward and RC anchors of the same window.
            auto normalizeWalk = [&](uint32_t ww) -> uint32_t {
                return (ww >= windowCount) ? (ww - windowCount) : ww;
            };
            auto walkAnchorHasInterWindow = [&](uint64_t a) -> bool {
                if(a >= anchorCount) return false;
                const uint32_t aWin = anchorToWindow[a];
                if(aWin == noWindow) return false;
                const uint32_t aNorm = normalizeWalk(aWin);
                BGL_FORALL_OUTEDGES(a, outE, anchorGraph, Shasta2AnchorGraph) {
                    if(!anchorGraph[outE].useForAssembly) continue;
                    const uint64_t tgt = uint64_t(target(outE, anchorGraph));
                    if(tgt < anchorCount) {
                        const uint32_t tgtWin = anchorToWindow[tgt];
                        if(tgtWin != noWindow && normalizeWalk(tgtWin) != aNorm)
                            return true;
                    }
                }
                BGL_FORALL_INEDGES(a, inE, anchorGraph, Shasta2AnchorGraph) {
                    if(!anchorGraph[inE].useForAssembly) continue;
                    const uint64_t src = uint64_t(source(inE, anchorGraph));
                    if(src < anchorCount) {
                        const uint32_t srcWin = anchorToWindow[src];
                        if(srcWin != noWindow && normalizeWalk(srcWin) != aNorm)
                            return true;
                    }
                }
                return false;
            };
            auto walkHasInterWindow = [&](uint64_t a) -> bool {
                if(walkAnchorHasInterWindow(a)) return true;
                const uint64_t rcA = a ^ 1ULL;
                if(rcA < anchorCount && walkAnchorHasInterWindow(rcA)) return true;
                return false;
            };

            if(walkForward) {
                for(uint64_t i = startPosIdx; i < positions.size(); i++) {
                    const Shasta2AnchorId aid = journey[positions[i]];

                    if(i != startPosIdx) {
                        if(walkHasInterWindow(uint64_t(aid))) break;
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

                    if(uint64_t(i) != startPosIdx) {
                        if(walkHasInterWindow(uint64_t(aid))) break;
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
        {
            BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
                if(!anchorGraph[e].useForAssembly) continue;
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
                        disableEdge(e);
                        ++case2RemovedCount;
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

    //removeDanglingWindowsIterative("post-detangle");

    // Diagnostic: dump per-window inter-window connectivity after all cleanup.
    // Check both normalized and raw window IDs, and count noWindow edges.
    {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        // Normalized check.
        std::map<uint32_t, std::set<uint32_t>> outNeighbors, inNeighbors;
        // Raw (un-normalized) check.
        std::map<uint32_t, std::set<uint32_t>> rawOutNeighbors, rawInNeighbors;
        // Edges involving noWindow anchors.
        uint64_t noWindowEdgeCount = 0;
        std::map<uint32_t, uint64_t> noWindowOutCount, noWindowInCount;

        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
            if(!anchorGraph[e].useForAssembly) continue;
            const uint64_t srcVal = uint64_t(source(e, anchorGraph));
            const uint64_t dstVal = uint64_t(target(e, anchorGraph));
            if(srcVal >= anchorCount || dstVal >= anchorCount) continue;
            const uint32_t srcWinRaw = anchorToWindow[srcVal];
            const uint32_t dstWinRaw = anchorToWindow[dstVal];

            // Track noWindow edges.
            if(srcWinRaw == noWindow || dstWinRaw == noWindow) {
                ++noWindowEdgeCount;
                if(srcWinRaw != noWindow && dstWinRaw == noWindow)
                    noWindowOutCount[normalize(srcWinRaw)]++;
                if(dstWinRaw != noWindow && srcWinRaw == noWindow)
                    noWindowInCount[normalize(dstWinRaw)]++;
                continue;
            }
            if(srcWinRaw == dstWinRaw) continue;

            // Raw.
            rawOutNeighbors[srcWinRaw].insert(dstWinRaw);
            rawInNeighbors[dstWinRaw].insert(srcWinRaw);

            // Normalized.
            const uint32_t srcWin = normalize(srcWinRaw);
            const uint32_t dstWin = normalize(dstWinRaw);
            if(srcWin == dstWin) continue;
            outNeighbors[srcWin].insert(dstWin);
            inNeighbors[dstWin].insert(srcWin);
        }

        cout << "Post-cleanup dangling check (normalized):" << endl;
        std::set<uint32_t> allWindows;
        for(const auto& [w, _] : outNeighbors) allWindows.insert(w);
        for(const auto& [w, _] : inNeighbors) allWindows.insert(w);
        uint64_t survivingDangling = 0;
        for(const uint32_t w : allWindows) {
            const bool hasIn = inNeighbors.count(w) && !inNeighbors[w].empty();
            const bool hasOut = outNeighbors.count(w) && !outNeighbors[w].empty();
            if(hasIn != hasOut) {
                ++survivingDangling;
                cout << "  DANGLING(norm) window " << w
                     << " in={";
                if(hasIn) for(const auto n : inNeighbors[w]) cout << n << " ";
                cout << "} out={";
                if(hasOut) for(const auto n : outNeighbors[w]) cout << n << " ";
                cout << "}" << endl;
            }
        }

        // Raw check.
        std::set<uint32_t> allRawWindows;
        for(const auto& [w, _] : rawOutNeighbors) allRawWindows.insert(w);
        for(const auto& [w, _] : rawInNeighbors) allRawWindows.insert(w);
        uint64_t rawSurvivingDangling = 0;
        for(const uint32_t w : allRawWindows) {
            const bool hasIn = rawInNeighbors.count(w) && !rawInNeighbors[w].empty();
            const bool hasOut = rawOutNeighbors.count(w) && !rawOutNeighbors[w].empty();
            if(hasIn != hasOut) {
                ++rawSurvivingDangling;
                cout << "  DANGLING(raw) window " << w
                     << (w >= windowCount ? " (RC)" : " (fwd)")
                     << " in={";
                if(hasIn) for(const auto n : rawInNeighbors[w]) cout << n << " ";
                cout << "} out={";
                if(hasOut) for(const auto n : rawOutNeighbors[w]) cout << n << " ";
                cout << "}" << endl;
            }
        }

        cout << "  Normalized: " << survivingDangling << " dangling, Raw: "
             << rawSurvivingDangling << " dangling, noWindow edges: "
             << noWindowEdgeCount << endl;
    }

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


