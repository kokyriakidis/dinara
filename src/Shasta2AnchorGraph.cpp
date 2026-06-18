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
            anchorPair.assertNoNegativeOffsets(anchors);
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



// Find detour window pairs.
// A detour occurs when window X has inter-window edges entering window W
// at backbone position i and exiting at position j > i. Returns the set
// of (W, X) pairs where this pattern exists. Used to suppress W→X→W
// transitions per-read during journey walks.
std::set<std::pair<uint32_t, uint32_t>> Shasta2AnchorGraph::findDetourWindowPairs(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Journeys& journeys) const
{
    const Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = num_vertices(anchorGraph);
    std::set<std::pair<uint32_t, uint32_t>> detourPairs;

    auto normalize = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    for(uint32_t w = 0; w < windowCount; w++) {
        const auto& window = anchorWindows[w];
        const auto& positions = window.filteredBackbonePositions;
        if(positions.size() < 2) continue;
        const auto journey = journeys[window.backboneOrientedReadId];

        // Collect inter-window edges at each backbone position.
        std::map<uint32_t, std::vector<uint32_t>> incomingBbIdx;  // xWindow → [bbIdx...]
        std::map<uint32_t, std::vector<uint32_t>> outgoingBbIdx;  // xWindow → [bbIdx...]

        for(uint32_t pi = 0; pi < positions.size(); pi++) {
            const uint64_t aid = uint64_t(journey[positions[pi]]);
            if(aid >= anchorCount) continue;

            auto ie = boost::in_edges(aid, anchorGraph);
            for(auto it = ie.first; it != ie.second; ++it) {
                if(!anchorGraph[*it].useForAssembly) continue;
                const uint64_t src = uint64_t(boost::source(*it, anchorGraph));
                if(src >= anchorCount) continue;
                const uint32_t srcRaw = anchorToWindow[src];
                if(srcRaw == noWindow) continue;
                const uint32_t srcW = normalize(srcRaw);
                if(srcW == w) continue;
                incomingBbIdx[srcW].push_back(pi);
            }

            auto oe = boost::out_edges(aid, anchorGraph);
            for(auto it = oe.first; it != oe.second; ++it) {
                if(!anchorGraph[*it].useForAssembly) continue;
                const uint64_t tgt = uint64_t(boost::target(*it, anchorGraph));
                if(tgt >= anchorCount) continue;
                const uint32_t tgtRaw = anchorToWindow[tgt];
                if(tgtRaw == noWindow) continue;
                const uint32_t tgtW = normalize(tgtRaw);
                if(tgtW == w) continue;
                outgoingBbIdx[tgtW].push_back(pi);
            }
        }

        // For each neighbor window X with both entry and exit,
        // check if any entry is at an earlier position than any exit.
        for(const auto& [xWindow, inIdxs] : incomingBbIdx) {
            auto outIt = outgoingBbIdx.find(xWindow);
            if(outIt == outgoingBbIdx.end()) continue;
            const auto& outIdxs = outIt->second;

            uint32_t minIn = *std::min_element(inIdxs.begin(), inIdxs.end());
            uint32_t maxOut = *std::max_element(outIdxs.begin(), outIdxs.end());

            if(maxOut > minIn) {
                // X detours through W.
                detourPairs.insert({w, xWindow});
            }
        }
    }

    cout << "findDetourWindowPairs: found " << detourPairs.size()
         << " detour window pairs." << endl;
    return detourPairs;
}


// Construct from anchor windows.
// Each window becomes a chain of its backbone anchors.
// Inter-window edges are discovered by walking read journeys.
Shasta2AnchorGraph::Shasta2AnchorGraph(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const vector<AnchorWindow>& anchorWindows,
    uint64_t minInterWindowCoverage,
    uint64_t minInterWindowEdgeCoverage,
    uint64_t threadCount,
    const Reads* reads,
    const vector<DetangleBypassEdge>* bypassEdges,
    const std::set<std::pair<uint32_t, uint32_t>>* detourWindowPairs,
    const vector<uint32_t>* anchorDovetailWindow) :
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

    // Edgeless-windows toggle. When true, inter-window edges are not created:
    // each window remains a bare anchor set (its intra-window backbone chain
    // is still built, but windows are not connected to each other). Downstream
    // stages still run on the resulting graph. Set to false to restore the
    // Stage A length-weighted reciprocal-best inter-window backbone.
    constexpr bool edgelessWindows = true;

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

    // Whole-journey dovetail membership (additive). anchorDovetailWindow maps
    // forward-oriented claimed dovetail anchors to their owning window. The
    // backbone fill above is authoritative (it uses filteredBackbonePositions),
    // so dovetails are only added where the anchor is still noWindow. This makes
    // each window's left/right dovetails part of the window for path-finding,
    // with the RC twin mapped to the mirror window. anchorToBackbonePos is left
    // at 0 for dovetails (they have no backbone position); transition selection
    // (which uses it) is expected to shift, as topology intentionally changes.
    if(anchorDovetailWindow != nullptr && !anchorDovetailWindow->empty()) {
        DINARA_ASSERT(anchorDovetailWindow->size() == anchorCount);
        uint64_t dovetailMapped = 0;
        for(uint64_t aid = 0; aid < anchorCount; aid++) {
            const uint32_t windowId = (*anchorDovetailWindow)[aid];
            if(windowId == noWindow) continue;
            if(anchorToWindow[aid] == noWindow) {
                anchorToWindow[aid] = windowId;
                ++dovetailMapped;
            }
            const uint64_t rcAid = aid ^ 1ULL;
            if(rcAid < anchorCount && anchorToWindow[rcAid] == noWindow) {
                anchorToWindow[rcAid] = windowId + windowCount;
            }
        }
        cout << "Shasta2AnchorGraph: mapped " << dovetailMapped
             << " dovetail anchors into windows." << endl;
    }

    // Helper to add an edge if the anchor pair has shared oriented reads.
    auto addEdgeIfValid = [&](Shasta2AnchorId anchorIdA, Shasta2AnchorId anchorIdB) -> bool {
        Shasta2AnchorPair anchorPair(anchors, anchorIdA, anchorIdB, false);
        anchorPair.assertNoNegativeOffsets(anchors);
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



    // Per-read window transition tracking.
    // If computeWindowTransitions() was already called, transitionReads
    // are populated. Recompute unconditionally so the graph's own
    // anchorToWindow (which may differ after detangling) is authoritative.
    {
        auto normalize = [&](uint32_t w) -> uint32_t {
            return (w >= windowCount) ? (w - windowCount) : w;
        };
        const uint32_t noW = AnchorWindowReadInterval::noWindow;

        // Clear any pre-existing transition data to avoid duplicates.
        for(uint32_t wid = 0; wid < windowCount; wid++) {
            auto& w = const_cast<AnchorWindow&>(anchorWindows[wid]);
            w.transitionReads.clear();
            w.backbonePreviousWindow = noW;
            w.backboneNextWindow = noW;
            for(auto& ri : w.readIntervals) {
                ri.previousWindow = noW;
                ri.nextWindow = noW;
            }
        }

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

    // Per-read transition: which anchor pair the read uses at the boundary.
    struct ReadTransition {
        uint32_t oidValue;
        Shasta2AnchorId firstAnchorInA;
        Shasta2AnchorId lastAnchorInA;
        Shasta2AnchorId firstAnchorInB;
        uint32_t firstJourneyPosInA;  // journey position of firstAnchorInA
        uint64_t supportingSpanProduct;
        uint64_t supportingSpanA;
        uint64_t supportingSpanB;
    };

    // Per window pair: all read transitions.
    std::map<std::pair<uint32_t, uint32_t>,
             std::vector<ReadTransition>> windowPairTransitions;

    // Per-read per-window: first and last base positions visited.
    struct ReadWindowSpan {
        uint32_t firstBasePos;
        uint32_t lastBasePos;
    };
    std::map<uint64_t, std::map<uint32_t, ReadWindowSpan>> readWindowSpans;

    // windowReads and readWindows are class members, populated here.

    uint64_t containedSkipCount = 0;
    const uint64_t journeyCount = journeys.size();
    for(uint64_t oidValue = 0; oidValue < journeyCount; oidValue++) {
        const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
        const auto journey = journeys[oid];
        if(journey.empty()) continue;

        if(reads) {
            const ReadId readId = oid.getReadId();
            if(readId < reads->readCount() && reads->getFlags(readId).isContained) {
                ++containedSkipCount;
            }
        }

        // Per-read journey walk with detour suppression.
        // When detourWindowPairs is provided and a read transitions from
        // window W to window X where (W, X) is a detour pair, we buffer
        // the transition. If the read returns to W, we discard the buffer
        // (the read was detouring through X). Otherwise we commit it.

        // First pass: collect (windowId, anchorId) for claimed anchors.
        struct WindowAnchor {
            uint32_t windowId;
            Shasta2AnchorId anchorId;
            uint32_t journeyPos;
        };
        std::vector<WindowAnchor> claimedAnchors;
        for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
            const Shasta2AnchorId anchorId = journey[pos];
            if(uint64_t(anchorId) >= anchorCount) continue;
            const uint32_t windowId = anchorToWindow[uint64_t(anchorId)];
            if(windowId == noWindow) continue;
            claimedAnchors.push_back({windowId, anchorId, pos});
        }

        // Build the effective window sequence, suppressing detour visits.
        // Each entry: (windowId, lastAnchorInWindow).
        struct WindowVisit {
            uint32_t windowId;
            Shasta2AnchorId firstAnchor;
            Shasta2AnchorId lastAnchor;
            uint32_t firstJourneyPos;
        };
        std::vector<WindowVisit> windowVisits;

        auto normalizeW = [&](uint32_t w) -> uint32_t {
            return (w >= windowCount) ? (w - windowCount) : w;
        };

        for(const auto& ca : claimedAnchors) {
            const uint32_t normalizedWin = normalizeW(ca.windowId);
            windowReads[normalizedWin].insert(uint32_t(oidValue));

            const uint32_t basePos = anchors.getPosition(ca.anchorId, oid);
            if(basePos != invalid<uint32_t>) {
                auto& spanMap = readWindowSpans[oidValue];
                auto spanIt = spanMap.find(ca.windowId);
                if(spanIt == spanMap.end()) {
                    spanMap[ca.windowId] = {basePos, basePos};
                } else {
                    spanIt->second.firstBasePos = std::min(spanIt->second.firstBasePos, basePos);
                    spanIt->second.lastBasePos = std::max(spanIt->second.lastBasePos, basePos);
                }
            }

            if(!windowVisits.empty() && windowVisits.back().windowId == ca.windowId) {
                // Same window — extend.
                windowVisits.back().lastAnchor = ca.anchorId;
            } else {
                // New window.
                windowVisits.push_back({ca.windowId, ca.anchorId, ca.anchorId, ca.journeyPos});
            }
        }

        // Suppress detour visits: scan for W→X→W patterns where (W, X)
        // is a detour pair. Remove the X visit and merge the two W visits.
        if(detourWindowPairs && !detourWindowPairs->empty()) {
            bool changed = true;
            while(changed) {
                changed = false;
                for(uint64_t i = 0; i + 2 < windowVisits.size(); i++) {
                    const uint32_t wNorm = normalizeW(windowVisits[i].windowId);
                    const uint32_t xNorm = normalizeW(windowVisits[i + 1].windowId);
                    const uint32_t w2Norm = normalizeW(windowVisits[i + 2].windowId);
                    if(wNorm != w2Norm) continue;
                    if(!detourWindowPairs->count({wNorm, xNorm})) continue;

                    // Suppress: merge W visits, remove X.
                    windowVisits[i].lastAnchor = windowVisits[i + 2].lastAnchor;
                    windowVisits.erase(windowVisits.begin() + int64_t(i + 1),
                                       windowVisits.begin() + int64_t(i + 3));
                    changed = true;
                    break;  // restart scan
                }
            }
        }

        // Emit transitions from the (possibly filtered) window sequence.
        for(uint64_t i = 0; i < windowVisits.size(); i++) {
            const auto& visit = windowVisits[i];
            readWindows[uint32_t(oidValue)].push_back(visit.windowId);
            if(i > 0) {
                const auto& prev = windowVisits[i - 1];
                auto key = std::make_pair(prev.windowId, visit.windowId);
                windowPairTransitions[key].push_back(
                    {uint32_t(oidValue), prev.firstAnchor, prev.lastAnchor, visit.firstAnchor,
                     prev.firstJourneyPos, 0, 0, 0});
            }
        }
    }

    // Second pass: finalize supporting spans.
    for(auto& [windowPair, transitions] : windowPairTransitions) {
        for(auto& t : transitions) {
            auto readSpanIt = readWindowSpans.find(t.oidValue);
            if(readSpanIt == readWindowSpans.end()) continue;
            const auto& spans = readSpanIt->second;
            uint64_t spanA = 0, spanB = 0;
            auto itA = spans.find(windowPair.first);
            auto itB = spans.find(windowPair.second);
            if(itA != spans.end()) spanA = itA->second.lastBasePos - itA->second.firstBasePos;
            if(itB != spans.end()) spanB = itB->second.lastBasePos - itB->second.firstBasePos;
            t.supportingSpanProduct = spanA * spanB;
            t.supportingSpanA = spanA;
            t.supportingSpanB = spanB;
        }
    }

    // Merge RC-mirror window pair transitions.
    {
        auto rcWindow = [&](uint32_t w) -> uint32_t {
            return (w >= windowCount) ? (w - windowCount) : (w + windowCount);
        };
        std::vector<std::pair<uint32_t, uint32_t>> keysToMerge;
        for(const auto& [wp, transitions] : windowPairTransitions) {
            const auto mirrorKey = std::make_pair(rcWindow(wp.second), rcWindow(wp.first));
            if(wp > mirrorKey) {
                keysToMerge.push_back(wp);
            }
        }
        uint64_t mergedTransitions = 0;
        for(const auto& key : keysToMerge) {
            const auto canonicalKey = std::make_pair(rcWindow(key.second), rcWindow(key.first));
            auto& srcTransitions = windowPairTransitions[key];
            auto& dstTransitions = windowPairTransitions[canonicalKey];
            for(auto& t : srcTransitions) {
                const Shasta2AnchorId newFirstInA =
                    Shasta2AnchorId(uint64_t(t.firstAnchorInB) ^ 1ULL);
                const Shasta2AnchorId newLastInA =
                    Shasta2AnchorId(uint64_t(t.firstAnchorInB) ^ 1ULL);
                const Shasta2AnchorId newFirstInB =
                    Shasta2AnchorId(uint64_t(t.lastAnchorInA) ^ 1ULL);
                // RC mirror: journey position not directly available,
                // set to 0 (RC transitions don't use it for edge creation).
                dstTransitions.push_back({
                    t.oidValue, newFirstInA, newLastInA, newFirstInB,
                    0, t.supportingSpanProduct, t.supportingSpanB, t.supportingSpanA});
                ++mergedTransitions;
            }
            windowPairTransitions.erase(key);
        }
        cout << "Merged " << mergedTransitions << " transitions from "
             << keysToMerge.size() << " non-canonical window pairs." << endl;
    }

    // Count shared reads between window pairs.
    auto countSharedReads = [&](uint32_t windowA, uint32_t windowB) -> uint64_t {
        const uint32_t normA = (windowA >= windowCount) ? (windowA - windowCount) : windowA;
        const uint32_t normB = (windowB >= windowCount) ? (windowB - windowCount) : windowB;
        auto itA = windowReads.find(normA);
        auto itB = windowReads.find(normB);
        if(itA == windowReads.end() || itB == windowReads.end()) return 0;
        const auto& setA = itA->second;
        const auto& setB = itB->second;
        const auto& smaller = (setA.size() <= setB.size()) ? setA : setB;
        const auto& larger  = (setA.size() <= setB.size()) ? setB : setA;
        uint64_t count = 0;
        for(const uint32_t r : smaller) {
            if(larger.count(r)) ++count;
        }
        return count;
    };

    // Diagnostics.
    {
        auto rcWindow = [&](uint32_t w) -> uint32_t {
            return (w >= windowCount) ? (w - windowCount) : (w + windowCount);
        };
        uint64_t totalTransitions = 0;
        uint64_t fwFwPairs = 0, fwRcPairs = 0, rcFwPairs = 0, rcRcPairs = 0;
        for(const auto& [wp, transitions] : windowPairTransitions) {
            totalTransitions += transitions.size();
            const bool aIsRc = (wp.first >= windowCount);
            const bool bIsRc = (wp.second >= windowCount);
            if(!aIsRc && !bIsRc) ++fwFwPairs;
            else if(!aIsRc && bIsRc) ++fwRcPairs;
            else if(aIsRc && !bIsRc) ++rcFwPairs;
            else ++rcRcPairs;
            const auto mirrorKey = std::make_pair(rcWindow(wp.second), rcWindow(wp.first));
            DINARA_ASSERT(wp <= mirrorKey);
        }
        cout << "Inter-window edge discovery (post-merge): " << windowPairTransitions.size()
             << " window pairs (" << fwFwPairs << " fw-fw, "
             << rcRcPairs << " rc-rc, "
             << fwRcPairs << " fw-rc, "
             << rcFwPairs << " rc-fw) with "
             << totalTransitions << " read transitions." << endl;
        if(reads) {
            cout << "Found " << containedSkipCount
                 << " contained oriented reads (included in inter-window edge discovery)." << endl;
        }
    }

    cout << "Inter-window discovery: " << windowPairTransitions.size()
         << " window pairs found." << endl;

    // Edgeless-windows: skip all inter-window edge creation. Windows stay as
    // bare anchor sets (intra-window backbone chains above are kept).
    if(edgelessWindows) {
        cout << "Edgeless-windows: skipping inter-window edge creation; "
             << "windows remain disjoint anchor sets." << endl;
    }

    // Stage A: length-weighted reciprocal-best inter-window edges.
    //
    // Start from all windows as nodes with no inter-window edges, then connect
    // A->B only when, scoring each window pair by the lengths of the distinct
    // reads that span it, B is A's best outgoing neighbor AND A is B's best
    // incoming neighbor (reciprocal best). This yields a clean linear local
    // backbone (degree <= 1 per side) without destructively dropping branches:
    // non-reciprocal pairs are simply deferred to Stage B (long-read bridging
    // across the mess), not deleted. "No single read may seed a backbone" is
    // honored by requiring at least minInterWindowCoverage distinct reads.
    //
    // Scoring uses summed spanning-read length, so longer reads carry more
    // weight (trust reads by length). windowPairTransitions is canonicalized
    // (mirror pairs merged and erased) while edge creation emits both the
    // forward and the RC mirror edge, so best-in/best-out is accumulated over
    // both the forward pair (A,B) and its mirror (rcWindow(B), rcWindow(A)).
    //
    // Gated by edgelessWindows: when edgeless, this entire inter-window edge
    // construction is skipped (kept here, not deleted, for easy restoration).
    if(!edgelessWindows)
    {
        auto rcWindow = [&](uint32_t w) -> uint32_t {
            return (w >= windowCount) ? (w - windowCount) : (w + windowCount);
        };
        auto readLengthOf = [&](uint32_t oidValue) -> uint64_t {
            if(reads == nullptr) return 1;  // fall back to unit weight (=count)
            const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
            const ReadId readId = oid.getReadId();
            if(readId >= reads->readCount()) return 1;
            return reads->getReadRawSequenceLength(readId);
        };

        // Per-pair: distinct spanning physical reads and summed read length.
        struct PairScore {
            uint64_t readCount = 0;
            uint64_t lengthScore = 0;
        };
        std::map<std::pair<uint32_t, uint32_t>, PairScore> pairScore;

        for(const auto& [windowPair, transitions] : windowPairTransitions) {
            // Distinct physical reads and their summed length.
            std::set<uint32_t> physicalReads;
            uint64_t lengthSum = 0;
            for(const auto& t : transitions) {
                if(physicalReads.insert(t.oidValue / 2).second) {
                    lengthSum += readLengthOf(t.oidValue);
                }
            }
            if(physicalReads.size() < minInterWindowCoverage) continue;

            const PairScore s{physicalReads.size(), lengthSum};
            pairScore[windowPair] = s;
            // RC mirror pair carries the same score (mirror edge is created too).
            const std::pair<uint32_t, uint32_t> mirror(
                rcWindow(windowPair.second), rcWindow(windowPair.first));
            if(mirror != windowPair) {
                pairScore[mirror] = s;
            }
        }

        // Best outgoing per source window and best incoming per dest window,
        // by length score (ties broken by read count then smaller window id).
        std::map<uint32_t, std::pair<uint32_t, PairScore>> bestOut; // src -> (dst, score)
        std::map<uint32_t, std::pair<uint32_t, PairScore>> bestIn;  // dst -> (src, score)
        auto better = [](const PairScore& x, uint32_t xId,
                         const PairScore& y, uint32_t yId) -> bool {
            if(x.lengthScore != y.lengthScore) return x.lengthScore > y.lengthScore;
            if(x.readCount != y.readCount) return x.readCount > y.readCount;
            return xId < yId;
        };
        for(const auto& [pair, s] : pairScore) {
            const uint32_t src = pair.first, dst = pair.second;
            auto itO = bestOut.find(src);
            if(itO == bestOut.end() ||
               better(s, dst, itO->second.second, itO->second.first)) {
                bestOut[src] = {dst, s};
            }
            auto itI = bestIn.find(dst);
            if(itI == bestIn.end() ||
               better(s, src, itI->second.second, itI->second.first)) {
                bestIn[dst] = {src, s};
            }
        }

        uint64_t interWindowCreated = 0;
        uint64_t interWindowSkipped = 0;
        uint64_t interWindowNotReciprocal = 0;

        for(const auto& [windowPair, transitions] : windowPairTransitions) {
            std::set<uint32_t> physicalReads;
            for(const auto& t : transitions) {
                physicalReads.insert(t.oidValue / 2);
            }
            if(physicalReads.size() < minInterWindowCoverage) {
                ++interWindowSkipped;
                continue;
            }

            // Reciprocal-best gate: B is A's best out AND A is B's best in.
            const uint32_t A = windowPair.first, B = windowPair.second;
            auto itO = bestOut.find(A);
            auto itI = bestIn.find(B);
            const bool reciprocal =
                itO != bestOut.end() && itO->second.first == B &&
                itI != bestIn.end()  && itI->second.first == A;
            if(!reciprocal) {
                ++interWindowNotReciprocal;
                continue;
            }

            Shasta2AnchorPair bestPair;
            uint64_t bestCoverage = 0;
            for(const auto& t : transitions) {
                Shasta2AnchorPair candidatePair(
                    anchors, t.lastAnchorInA, t.firstAnchorInB, false);
                candidatePair.assertNoNegativeOffsets(anchors);
                if(candidatePair.size() > bestCoverage) {
                    bestCoverage = candidatePair.size();
                    bestPair = std::move(candidatePair);
                }
            }

            if(bestCoverage == 0) continue;
            if(bestCoverage < minInterWindowEdgeCoverage) continue;

            if(addEdgeIfValid(bestPair.anchorIdA, bestPair.anchorIdB)) {
                ++interWindowCreated;
                // RC mirror edge.
                const Shasta2AnchorId rcA =
                    Shasta2AnchorId(uint64_t(bestPair.anchorIdA) ^ 1ULL);
                const Shasta2AnchorId rcB =
                    Shasta2AnchorId(uint64_t(bestPair.anchorIdB) ^ 1ULL);
                if(uint64_t(rcA) < anchorCount && uint64_t(rcB) < anchorCount) {
                    addEdgeIfValid(rcB, rcA);
                }
            }
        }

        cout << "Inter-window edges (Stage A, length-weighted reciprocal-best): "
             << interWindowCreated << " created, "
             << interWindowSkipped << " skipped (< " << minInterWindowCoverage
             << " reads), " << interWindowNotReciprocal
             << " deferred (not reciprocal-best)." << endl;
    }

    // All-edge construction (disabled): created one edge per coverage-passing
    // window pair regardless of score. Replaced by the Stage A length-weighted
    // reciprocal-best block above to build a clean linear backbone first.
#if 0
    {
        uint64_t interWindowCreated = 0;
        uint64_t interWindowSkipped = 0;

        for(const auto& [windowPair, transitions] : windowPairTransitions) {
            // Count distinct physical reads (deduplicate by read ID).
            std::set<uint32_t> physicalReads;
            for(const auto& t : transitions) {
                physicalReads.insert(t.oidValue / 2);
            }
            if(physicalReads.size() < minInterWindowCoverage) {
                ++interWindowSkipped;
                continue;
            }

            Shasta2AnchorPair bestPair;
            uint64_t bestCoverage = 0;
            for(const auto& t : transitions) {
                Shasta2AnchorPair candidatePair(
                    anchors, t.lastAnchorInA, t.firstAnchorInB, false);
                candidatePair.assertNoNegativeOffsets(anchors);
                if(candidatePair.size() > bestCoverage) {
                    bestCoverage = candidatePair.size();
                    bestPair = std::move(candidatePair);
                }
            }

            if(bestCoverage == 0) continue;
            if(bestCoverage < minInterWindowEdgeCoverage) continue;

            if(addEdgeIfValid(bestPair.anchorIdA, bestPair.anchorIdB)) {
                ++interWindowCreated;
                // RC mirror edge.
                const Shasta2AnchorId rcA =
                    Shasta2AnchorId(uint64_t(bestPair.anchorIdA) ^ 1ULL);
                const Shasta2AnchorId rcB =
                    Shasta2AnchorId(uint64_t(bestPair.anchorIdB) ^ 1ULL);
                if(uint64_t(rcA) < anchorCount && uint64_t(rcB) < anchorCount) {
                    addEdgeIfValid(rcB, rcA);
                }
            }
        }

        cout << "Inter-window edges: " << interWindowCreated << " created, "
             << interWindowSkipped << " skipped (< " << minInterWindowCoverage
             << " reads)." << endl;
    }
#endif

    // Old inter-window edge creation (disabled).
#if 0
    uint64_t interWindowZeroPairs = 0;
    uint64_t interWindowLowEdgeCoverage = 0;
    uint64_t interWindowCreated = 0;
    struct InterWindowEdgeInfo {
        std::pair<uint32_t, uint32_t> windowPair;
        Shasta2AnchorId anchorIdA;
        Shasta2AnchorId anchorIdB;
        uint64_t readCount;
    };
    std::vector<InterWindowEdgeInfo> createdEdges;

    // Helper: create an inter-window edge and its RC mirror.
    auto createInterWindowEdge = [&](
        const std::pair<uint32_t, uint32_t>& windowPair,
        Shasta2AnchorPair& bestPair,
        uint64_t sharedReads,
        uint64_t spanPrev,
        uint64_t spanNext)
    {
        DINARA_ASSERT(anchors.countCommon(bestPair.anchorIdA, bestPair.anchorIdB) > 0);
        // Forward edge.
        edge_descriptor e;
        tie(e, ignore) = add_edge(
            bestPair.anchorIdA,
            bestPair.anchorIdB,
            Shasta2AnchorGraphEdge(bestPair, bestPair.getAverageOffset(anchors), nextEdgeId++),
            anchorGraph);
        anchorGraph[e].useForAssembly = true;
        anchorGraph[e].supportingSpanPrev = spanPrev;
        anchorGraph[e].supportingSpanNext = spanNext;
        anchorGraph[e].sharedReadCount = sharedReads;
        createdEdges.push_back({windowPair, bestPair.anchorIdA, bestPair.anchorIdB, sharedReads});
        ++interWindowCreated;

        // RC mirror edge: reverse the anchor pair and flip both anchor IDs.
        // Spans swap: the RC mirror's prev is the forward's next and vice versa.
        // Skip for self-RC-mirror window pairs (W, W+wc) where the RC mirror
        // edge would land on the same window pair, creating a duplicate.
        const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(bestPair.anchorIdA) ^ 1ULL);
        const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(bestPair.anchorIdB) ^ 1ULL);
        const uint32_t rcSrc = (windowPair.second < windowCount)
            ? windowPair.second + windowCount : windowPair.second - windowCount;
        const uint32_t rcDst = (windowPair.first < windowCount)
            ? windowPair.first + windowCount : windowPair.first - windowCount;
        const bool isSelfRcMirror =
            (rcSrc == windowPair.first && rcDst == windowPair.second);
        if(!isSelfRcMirror &&
           uint64_t(rcA) < anchorCount && uint64_t(rcB) < anchorCount) {
            Shasta2AnchorPair rcPair(anchors, rcB, rcA, false);
            rcPair.assertNoNegativeOffsets(anchors);
            DINARA_ASSERT(rcPair.size() > 0);
            edge_descriptor eRc;
            tie(eRc, ignore) = add_edge(
                rcB, rcA,
                Shasta2AnchorGraphEdge(rcPair, rcPair.getAverageOffset(anchors), nextEdgeId++),
                anchorGraph);
            anchorGraph[eRc].useForAssembly = true;
            anchorGraph[eRc].supportingSpanPrev = spanNext;
            anchorGraph[eRc].supportingSpanNext = spanPrev;
            anchorGraph[eRc].sharedReadCount = sharedReads;
            createdEdges.push_back({{rcSrc, rcDst}, rcB, rcA, sharedReads});
            ++interWindowCreated;
        }
    };

    // For each window pair (A, B), create one edge using:
    //   - The latest lastAnchorInA across all reads (highest backbone pos in A)
    //   - The earliest firstAnchorInB across all reads (lowest backbone pos in B)
    // This maximizes backbone retention on both sides after trimming.
    // Backbone positions always increase from backboneBegin to backboneEnd
    // regardless of forward/RC, so "latest" = highest pos and "earliest" =
    // lowest pos in all cases.
    uint64_t interWindowLowCoverage = 0;
    for(const auto& [windowPair, transitions] : windowPairTransitions) {

        // Skip window pairs with insufficient read support.
        if(transitions.size() < minInterWindowCoverage) {
            ++interWindowLowCoverage;
            continue;
        }

        // Pick the latest lastAnchorInA (highest backbone position in A).
        // This maximizes the number of reads in A that can reach the
        // inter-window edge. Among transitions sharing that latest A
        // anchor, pick the firstAnchorInB that yields the highest
        // AnchorPair coverage.
        uint32_t bestPosA = 0;
        for(const auto& t : transitions) {
            const uint32_t posA = anchorToBackbonePos[uint64_t(t.lastAnchorInA)];
            if(posA > bestPosA) {
                bestPosA = posA;
            }
        }

        // Collect all transitions at the best A position and pick the
        // one whose (lastInA, firstInB) pair has the highest coverage.
        Shasta2AnchorPair anchorPair;
        uint64_t bestCoverage = 0;
        for(const auto& t : transitions) {
            if(anchorToBackbonePos[uint64_t(t.lastAnchorInA)] != bestPosA) continue;
            Shasta2AnchorPair candidatePair(anchors, t.lastAnchorInA, t.firstAnchorInB, false);
            candidatePair.assertNoNegativeOffsets(anchors);
            if(candidatePair.size() > bestCoverage) {
                bestCoverage = candidatePair.size();
                anchorPair = std::move(candidatePair);
            }
        }

        if(bestCoverage == 0) {
            ++interWindowZeroPairs;
            continue;
        }

        // Skip edges with insufficient anchor pair coverage.
        if(anchorPair.size() < minInterWindowEdgeCoverage) {
            ++interWindowLowEdgeCoverage;
            continue;
        }

        const uint64_t sharedReads = countSharedReads(windowPair.first, windowPair.second);
        // Use average spans across all transitions for the edge attributes.
        uint64_t totalSpanA = 0, totalSpanB = 0;
        for(const auto& t : transitions) {
            totalSpanA += t.supportingSpanA;
            totalSpanB += t.supportingSpanB;
        }
        const uint64_t avgSpanA = totalSpanA / transitions.size();
        const uint64_t avgSpanB = totalSpanB / transitions.size();
        createInterWindowEdge(windowPair, anchorPair, sharedReads, avgSpanA, avgSpanB);
    }

    cout << "Inter-window edges: " << interWindowCreated << " created, "
         << interWindowLowCoverage << " rejected (< " << minInterWindowCoverage << " reads), "
         << interWindowLowEdgeCoverage << " rejected (< " << minInterWindowEdgeCoverage << " edge coverage), "
         << interWindowZeroPairs << " rejected (no valid anchor pair)." << endl;
#endif  // Inter-window edge creation disabled

    // ========================================================================
    // Filter lambdas (defined here, called in order below).
    // ========================================================================

    // Use the member function disableEdge() for all edge disabling.
    // (Defined at the end of this file.)

    // Case 2: Remove single internal-edge connections between windows.
    // If two windows are connected by exactly one inter-window edge
    // and that edge is not an endpoint connection (not matching
    // backbonePreviousWindow or backboneNextWindow of either window),
    // it's a spurious single-point connection — delete it.
    auto runSingleEdgeFilter = [&]() {
        auto normalize = [&](uint32_t w2) -> uint32_t {
            return (w2 >= windowCount) ? (w2 - windowCount) : w2;
        };

        // Count unique anchor pairs per (normalized) window pair.
        // fw-fw and rc-rc edges for the same logical connection share
        // the same normalized anchor pair (min, max), so they count as one.
        struct WindowPairInfo {
            std::set<std::pair<uint64_t, uint64_t>> uniqueAnchorPairs;
            std::vector<edge_descriptor> edges;
        };
        std::map<std::pair<uint32_t, uint32_t>, WindowPairInfo> pairInfo;
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
            auto& info = pairInfo[key];
            info.edges.push_back(e);
            // Normalize anchor pair: use (min, max) so fw-fw and rc-rc
            // edges with anchors (a, b) and (a^1, b^1) map to the same pair.
            const uint64_t normSrc = std::min(srcVal, uint64_t(srcVal ^ 1ULL));
            const uint64_t normDst = std::min(dstVal, uint64_t(dstVal ^ 1ULL));
            const auto ap = (normSrc < normDst)
                ? std::make_pair(normSrc, normDst)
                : std::make_pair(normDst, normSrc);
            info.uniqueAnchorPairs.insert(ap);
        }

        // Helper: check if an anchor has both an incoming and outgoing
        // active intra-window edge (same normalized window).
        auto hasBothIntraWindowEdges = [&](uint64_t aid) -> bool {
            if(aid >= anchorCount) return false;
            const uint32_t aidRaw = anchorToWindow[aid];
            if(aidRaw == noWindow) return false;
            const uint32_t aidW = normalize(aidRaw);
            bool hasOut = false, hasIn = false;
            auto oe = boost::out_edges(aid, anchorGraph);
            for(auto it = oe.first; it != oe.second; ++it) {
                if(!anchorGraph[*it].useForAssembly) continue;
                const uint64_t tgt = uint64_t(target(*it, anchorGraph));
                if(tgt >= anchorCount) continue;
                const uint32_t tgtRaw = anchorToWindow[tgt];
                if(tgtRaw == noWindow) continue;
                if(normalize(tgtRaw) == aidW) { hasOut = true; break; }
            }
            auto ie = boost::in_edges(aid, anchorGraph);
            for(auto it = ie.first; it != ie.second; ++it) {
                if(!anchorGraph[*it].useForAssembly) continue;
                const uint64_t src = uint64_t(source(*it, anchorGraph));
                if(src >= anchorCount) continue;
                const uint32_t srcRaw = anchorToWindow[src];
                if(srcRaw == noWindow) continue;
                if(normalize(srcRaw) == aidW) { hasIn = true; break; }
            }
            return hasOut && hasIn;
        };

        uint64_t singleEdgeRemovedCount = 0;
        for(const auto& [pair, info] : pairInfo) {
            // 1. Must have only one unique anchor pair (single-point connection).
            if(info.uniqueAnchorPairs.size() != 1) continue;

            // 2. Both anchors must have both incoming and outgoing intra-window
            // edges, so removing this inter-window edge won't create a dead end.
            const uint64_t srcVal = uint64_t(source(info.edges[0], anchorGraph));
            const uint64_t dstVal = uint64_t(target(info.edges[0], anchorGraph));
            if(!hasBothIntraWindowEdges(srcVal) || !hasBothIntraWindowEdges(dstVal)) continue;

            // Remove all edges for this single-connection pair.
            for(const auto& e : info.edges) {
                disableEdge(e);
            }
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
            for(const auto& [xWindow, inEdges] : incoming) {
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

            std::set<uint32_t> orphanedIndices;

            for(uint32_t startIdx = 0; startIdx < bbAnchors.size(); startIdx++) {
                if(orphanedIndices.count(startIdx)) continue;
                const uint64_t startAid = bbAnchors[startIdx];

                // Check for internal (non-endpoint) inter-window outgoing edges.
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

                // Seed BFS with internal (non-endpoint) inter-window neighbors.
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

                    // Disable all edges of the orphaned anchor (intra- and inter-window).
                    auto disableAllEdgesOf = [&](uint64_t a) {
                        if(a >= anchorCount) return;
                        auto oe4 = boost::out_edges(a, anchorGraph);
                        for(auto it = oe4.first; it != oe4.second; ++it) {
                            if(anchorGraph[*it].useForAssembly) {
                                disableEdge(*it);
                            }
                        }
                        auto ie4 = boost::in_edges(a, anchorGraph);
                        for(auto it = ie4.first; it != ie4.second; ++it) {
                            if(anchorGraph[*it].useForAssembly) {
                                disableEdge(*it);
                            }
                        }
                    };
                    disableAllEdgesOf(aid);
                    disableAllEdgesOf(aid ^ 1ULL);

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
    // trimBackbones is now a member function — called via this->trimBackbones().

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

            // Dead-end spur: window w has only one neighbor. Disable all
            // inter-window edges of w.
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
    // backbone ends exposed by edge removal.
    // backbonePreviousWindow/backboneNextWindow are set once during
    // construction and do not change (recomputeBackboneEndpoints is not
    // needed since it always returns the same result).
    // ========================================================================
    // All filtering moved to main.cpp — called externally after construction.
#if 0
    trimBackbones(anchorWindows, journeys);
    removeDeadEndSpurs();           // Remove single-neighbor internal spurs
    trimBackbones(anchorWindows, journeys);
    runSingleEdgeFilter();          // Case 2: remove single-point connections
    trimBackbones(anchorWindows, journeys);
    runBypassDetourFilter();        // Case 1: bypass detours through other windows
    trimBackbones(anchorWindows, journeys);
    runBubblePopFilter();           // Case 3: pop bubbles returning to same window
    trimBackbones(anchorWindows, journeys);
    runShortcutFilter();
    trimBackbones(anchorWindows, journeys);
    runCrossWindowFilter();
    trimBackbones(anchorWindows, journeys);
    removeIsolatedWindows();
    trimBackbones(anchorWindows, journeys);
    removeSmallWindows(2);
    trimBackbones(anchorWindows, journeys);
    removeDanglingWindowsIterative("post-filter");
#endif

    // Inter-window edge population and bypass edges disabled.
#if 0
    // Populate per-window outEdges/inEdges from createdEdges.
    for(const auto& edgeInfo : createdEdges) {
        const uint32_t srcW = edgeInfo.windowPair.first;
        const uint32_t dstW = edgeInfo.windowPair.second;
        if(srcW < windowCount) {
            auto& w = const_cast<AnchorWindow&>(anchorWindows[srcW]);
            w.outEdges.push_back({dstW, edgeInfo.anchorIdA, edgeInfo.anchorIdB, edgeInfo.readCount});
        }
        if(dstW < windowCount) {
            auto& w = const_cast<AnchorWindow&>(anchorWindows[dstW]);
            w.inEdges.push_back({srcW, edgeInfo.anchorIdA, edgeInfo.anchorIdB, edgeInfo.readCount});
        }
    }



    // Bypass edges from detangling: direct connections that skip over
    // tangled windows whose flow reads have been removed.
    // No coverage threshold — bypass edges are validated by the detangle
    // algorithm (journey walk, X-span matching, anchor pair selection).
    // Requiring minInterWindowCoverage here would leave tangles unresolved.
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

    // Recount edge types (only active edges).
    uint64_t finalIntraCount = 0;
    uint64_t finalInterCount = 0;
    uint64_t finalDisabledCount = 0;
    uint64_t finalActiveCount = 0;
    std::set<Shasta2AnchorId> activeVertices;
    BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraph) {
        if(!anchorGraph[e].useForAssembly) {
            ++finalDisabledCount;
            continue;
        }
        ++finalActiveCount;
        activeVertices.insert(source(e, anchorGraph));
        activeVertices.insert(target(e, anchorGraph));
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

    cout << "The anchor graph has " << activeVertices.size()
         << " active vertices (" << num_vertices(*this) << " total), "
         << finalActiveCount << " active edges"
         << " (" << finalIntraCount << " intra-window, "
         << finalInterCount << " inter-window, "
         << finalDisabledCount << " disabled)." << endl;

    // Verify all edges: for each orientedReadId, check its base position
    // on both anchors. Report edges where any read has negative base offset
    // or is missing from one of the anchors.
    {
        uint64_t backwardEdgeCount = 0;
        uint64_t missingReadEdgeCount = 0;
        BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraphBaseClass) {
            const auto& dEdge = anchorGraph[e];
            if(!dEdge.useForAssembly) continue;
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

    // Trim backbone ends that extend beyond inter-window connection points.
    // Disabled: the bidirected pipeline handles this differently.
    // trimBackbones(anchorWindows, journeys);
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


void Shasta2AnchorGraph::writeWindowConnectionStats(
    const Shasta2Anchors& anchors,
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Journeys& journeys) const
{
    const Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = num_vertices(anchorGraph);

    auto normalizeW = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // Find the biggest window by baseSpan.
    uint32_t biggestWindow = 0;
    uint64_t biggestSpan = 0;
    for(uint32_t w = 0; w < windowCount; w++) {
        if(anchorWindows[w].baseSpan > biggestSpan) {
            biggestSpan = anchorWindows[w].baseSpan;
            biggestWindow = w;
        }
    }

    const auto& window = anchorWindows[biggestWindow];
    const uint32_t bbBegin = window.backboneBegin;
    const uint32_t bbEnd = window.backboneEnd;
    const uint32_t bbAnchors = bbEnd - bbBegin;

    cout << "\n=== Window connection stats for window " << biggestWindow
         << " (baseSpan=" << biggestSpan
         << ", backboneAnchors=" << bbAnchors
         << ", backbone=" << window.backboneOrientedReadId
         << ", prevWindow=" << window.backbonePreviousWindow
         << ", nextWindow=" << window.backboneNextWindow
         << ") ===" << endl;

    // Collect all inter-window edges involving this window.
    struct NeighborEdgeInfo {
        uint32_t neighborNormWindow;
        bool isOutgoing; // true = this window is source
        edge_descriptor e;
        uint64_t coverage;
        uint64_t offset;
        uint64_t sharedReadCount;
        uint64_t spanInThis;  // supportingSpan on this window's side
        uint64_t spanInOther; // supportingSpan on neighbor's side
        Shasta2AnchorId anchorInThis;
        Shasta2AnchorId anchorInOther;
    };
    vector<NeighborEdgeInfo> edgeInfos;

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

        const auto& edge = anchorGraph[e];

        if(srcNorm == biggestWindow) {
            edgeInfos.push_back({
                dstNorm, true, e,
                edge.coverage(), edge.offset, edge.sharedReadCount,
                edge.supportingSpanPrev, edge.supportingSpanNext,
                Shasta2AnchorId(src), Shasta2AnchorId(dst)
            });
        } else if(dstNorm == biggestWindow) {
            edgeInfos.push_back({
                srcNorm, false, e,
                edge.coverage(), edge.offset, edge.sharedReadCount,
                edge.supportingSpanNext, edge.supportingSpanPrev,
                Shasta2AnchorId(dst), Shasta2AnchorId(src)
            });
        }
    }

    // Group by neighbor window.
    std::map<uint32_t, vector<const NeighborEdgeInfo*>> byNeighbor;
    for(const auto& info : edgeInfos) {
        byNeighbor[info.neighborNormWindow].push_back(&info);
    }

    cout << "Total inter-window edges: " << edgeInfos.size()
         << " to " << byNeighbor.size() << " neighbor windows." << endl;

    // For each neighbor, show details.
    for(const auto& [neighborW, infos] : byNeighbor) {
        const bool neighborIsBackbonePrev =
            (window.backbonePreviousWindow != AnchorWindowReadInterval::noWindow &&
             normalizeW(window.backbonePreviousWindow) == neighborW);
        const bool neighborIsBackboneNext =
            (window.backboneNextWindow != AnchorWindowReadInterval::noWindow &&
             normalizeW(window.backboneNextWindow) == neighborW);

        cout << "\n  Neighbor window " << neighborW
             << " (baseSpan=" << anchorWindows[neighborW].baseSpan << ")";
        if(neighborIsBackbonePrev) cout << " [BACKBONE-PREV]";
        if(neighborIsBackboneNext) cout << " [BACKBONE-NEXT]";
        cout << ":" << endl;

        for(const auto* info : infos) {
            const string dir = info->isOutgoing ? "OUT" : "IN";

            // Find backbone position of the anchor in this window.
            const uint64_t anchorInThisId = uint64_t(info->anchorInThis);
            const auto bbJourney = journeys[window.backboneOrientedReadId];
            uint32_t bbPos = std::numeric_limits<uint32_t>::max();
            for(uint32_t p = bbBegin; p < bbEnd; p++) {
                if(uint64_t(bbJourney[p]) == anchorInThisId) {
                    bbPos = p;
                    break;
                }
            }

            // Position relative to backbone: how far from start/end.
            const string posLabel = (bbPos != std::numeric_limits<uint32_t>::max())
                ? ("bbPos=" + std::to_string(bbPos) +
                   " fromStart=" + std::to_string(bbPos - bbBegin) +
                   " fromEnd=" + std::to_string(bbEnd - 1 - bbPos))
                : "bbPos=N/A";

            cout << "    " << dir
                 << " cov=" << info->coverage
                 << " offset=" << info->offset
                 << " sharedReads=" << info->sharedReadCount
                 << " spanInThis=" << info->spanInThis
                 << " spanInOther=" << info->spanInOther
                 << " anchorInThis=" << anchorInThisId
                 << " anchorInOther=" << uint64_t(info->anchorInOther)
                 << " " << posLabel
                 << endl;
        }
    }

    // Show transition reads summary.
    cout << "\n  Transition reads (from transitionReads map):" << endl;
    for(const auto& [key, reads] : window.transitionReads) {
        const uint32_t prev = key.first;
        const uint32_t next = key.second;
        const string prevStr = (prev == AnchorWindowReadInterval::noWindow)
            ? "START" : std::to_string(normalizeW(prev));
        const string nextStr = (next == AnchorWindowReadInterval::noWindow)
            ? "END" : std::to_string(normalizeW(next));
        cout << "    " << prevStr << " -> [" << biggestWindow << "] -> " << nextStr
             << " : " << reads.size() << " reads" << endl;
    }

    cout << "=== End window " << biggestWindow << " stats ===\n" << endl;
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


uint64_t Shasta2AnchorGraph::removeInternalConnections(
    const Shasta2Anchors& anchors,
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Journeys& journeys)
{
    Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = num_vertices(anchorGraph);

    auto normalizeW = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // Collect backbone oriented read IDs so we can skip them.
    std::set<uint64_t> backboneReadIds;
    for(uint32_t w = 0; w < windowCount; w++) {
        backboneReadIds.insert(anchorWindows[w].backboneOrientedReadId.getValue());
        backboneReadIds.insert(anchorWindows[w].backboneOrientedReadId.getValue() ^ 1ULL);
    }

    // Build an index of inter-window edges by (srcNormWindow, dstNormWindow).
    using WindowPairKey = std::pair<uint32_t, uint32_t>;
    std::map<WindowPairKey, vector<edge_descriptor>> windowPairEdges;
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
        windowPairEdges[{srcNorm, dstNorm}].push_back(e);
    }

    // Collect non-backbone reads.
    vector<OrientedReadId> nonBackboneReads;
    for(uint64_t readId = 0; readId < journeys.size(); readId++) {
        if(backboneReadIds.count(readId)) continue;
        const auto journey = journeys[OrientedReadId::fromValue(uint32_t(readId))];
        if(journey.size() < 2) continue;
        nonBackboneReads.push_back(OrientedReadId::fromValue(uint32_t(readId)));
    }

    cout << "removeInternalConnections: processing " << nonBackboneReads.size()
         << " non-backbone reads, looking for A->B->A triplets." << endl;

    // Track which window pairs have had edges removed or bypasses created.
    std::set<WindowPairKey> removedPairs;
    std::set<WindowPairKey> createdBypasses;

    uint64_t totalRemoved = 0;
    uint64_t totalBypasses = 0;
    uint64_t totalRepeats = 0;

    struct WindowVisit {
        uint32_t normWindow;
        Shasta2AnchorId firstAnchor;
        Shasta2AnchorId lastAnchor;
    };

    // Helper to disable all edges for a window pair.
    auto disableWindowPair = [&](uint32_t wFrom, uint32_t wTo) {
        if(removedPairs.count({wFrom, wTo})) return;
        removedPairs.insert({wFrom, wTo});
        auto it = windowPairEdges.find({wFrom, wTo});
        if(it != windowPairEdges.end()) {
            for(const edge_descriptor e : it->second) {
                if(anchorGraph[e].useForAssembly) {
                    disableEdge(e);
                    ++totalRemoved;
                }
            }
        }
    };

    for(const OrientedReadId& orientedReadId : nonBackboneReads) {
        const auto journey = journeys[orientedReadId];

        // Build the window sequence for this read.
        vector<WindowVisit> ws;
        uint32_t prevNormWindow = noWindow;
        for(uint64_t pos = 0; pos < journey.size(); pos++) {
            const uint64_t aid = uint64_t(journey[pos]);
            if(aid >= anchorCount) continue;
            const uint32_t w = anchorToWindow[aid];
            if(w == noWindow) continue;
            const uint32_t normW = normalizeW(w);
            if(normW != prevNormWindow) {
                ws.push_back({normW, journey[pos], journey[pos]});
                prevNormWindow = normW;
            } else if(!ws.empty()) {
                ws.back().lastAnchor = journey[pos];
            }
        }

        // Repeatedly find and collapse the shortest repeated window span.
        // This ensures inside-out resolution: inner patterns (e.g. B->C->B)
        // are resolved before outer ones (e.g. A->B->A), creating the
        // correct bypass edges at each level.
        while(ws.size() >= 3) {
            // Find the shortest span (i,j) where ws[i] == ws[j], j >= i+2.
            uint64_t bestI = 0, bestJ = 0, bestSpan = ws.size() + 1;
            for(uint64_t i = 0; i + 2 <= ws.size(); i++) {
                const uint32_t wA = ws[i].normWindow;
                for(uint64_t j = i + 2; j < ws.size(); j++) {
                    if(ws[j].normWindow == wA) {
                        const uint64_t span = j - i;
                        if(span < bestSpan) {
                            bestI = i;
                            bestJ = j;
                            bestSpan = span;
                        }
                        break; // nearest repeat of ws[i], no need to look further
                    }
                }
                if(bestSpan == 2) break; // can't do better than a triplet
            }

            if(bestSpan > ws.size()) break; // no repeated window found

            const uint64_t i = bestI;
            const uint64_t j = bestJ;
            const uint32_t wA = ws[i].normWindow;

            ++totalRepeats;

            // Disable all inter-window edges along the detour path.
            for(uint64_t k = i; k < j; k++) {
                const uint32_t wFrom = ws[k].normWindow;
                const uint32_t wTo = ws[k + 1].normWindow;
                disableWindowPair(wFrom, wTo);
                disableWindowPair(wTo, wFrom);
            }

            // Create a bypass edge from A (last anchor of first visit)
            // to A (first anchor of second visit).
            if(!createdBypasses.count({wA, wA})) {
                createdBypasses.insert({wA, wA});

                const Shasta2AnchorId bypassFrom = ws[i].lastAnchor;
                const Shasta2AnchorId bypassTo = ws[j].firstAnchor;

                if(uint64_t(bypassFrom) != uint64_t(bypassTo) &&
                   anchors.countCommon(bypassFrom, bypassTo) > 0) {
                    Shasta2AnchorPair bypassPair(anchors, bypassFrom, bypassTo, false);
                    bypassPair.assertNoNegativeOffsets(anchors);
                    if(bypassPair.size() > 0) {

                        // Create forward edge.
                        edge_descriptor eBypass;
                        tie(eBypass, ignore) = add_edge(
                            uint64_t(bypassFrom), uint64_t(bypassTo),
                            Shasta2AnchorGraphEdge(bypassPair,
                                bypassPair.getAverageOffset(anchors), nextEdgeId++),
                            anchorGraph);
                        anchorGraph[eBypass].useForAssembly = true;

                        // Create RC mirror edge.
                        const Shasta2AnchorId rcFrom = Shasta2AnchorId(uint64_t(bypassFrom) ^ 1ULL);
                        const Shasta2AnchorId rcTo = Shasta2AnchorId(uint64_t(bypassTo) ^ 1ULL);
                        if(uint64_t(rcFrom) < anchorCount && uint64_t(rcTo) < anchorCount) {
                            Shasta2AnchorPair rcPair(anchors, rcTo, rcFrom, false);
                            rcPair.assertNoNegativeOffsets(anchors);
                            if(rcPair.size() > 0) {
                                edge_descriptor eRc;
                                tie(eRc, ignore) = add_edge(
                                    uint64_t(rcTo), uint64_t(rcFrom),
                                    Shasta2AnchorGraphEdge(rcPair,
                                        rcPair.getAverageOffset(anchors), nextEdgeId++),
                                    anchorGraph);
                                anchorGraph[eRc].useForAssembly = true;
                            }
                        }
                        ++totalBypasses;
                    }
                }
            }

            // Collapse: merge the two A visits, remove intermediates.
            ws[i].lastAnchor = ws[j].lastAnchor;
            ws.erase(ws.begin() + int64_t(i + 1), ws.begin() + int64_t(j + 1));
        }
    }

    cout << "removeInternalConnections: found " << totalRepeats
         << " repeated window visits across " << removedPairs.size()
         << " window pairs, removed " << totalRemoved
         << " internal edges, created " << totalBypasses
         << " bypass edges." << endl;
    return totalRemoved + totalBypasses;
}


uint64_t Shasta2AnchorGraph::trimBackbones(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Journeys& journeys)
{
    Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = num_vertices(anchorGraph);

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

    auto normalizeW = [&](uint32_t w2) -> uint32_t {
        return (w2 >= windowCount) ? (w2 - windowCount) : w2;
    };

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

    auto hasIncomingInterWindowEdge = [&](uint64_t aid) -> bool {
        if(anchorHasInterWindowInEdge(aid)) return true;
        const uint64_t rcAid = aid ^ 1ULL;
        return (rcAid < anchorCount && anchorHasInterWindowOutEdge(rcAid));
    };

    auto hasOutgoingInterWindowEdge = [&](uint64_t aid) -> bool {
        if(anchorHasInterWindowOutEdge(aid)) return true;
        const uint64_t rcAid = aid ^ 1ULL;
        return (rcAid < anchorCount && anchorHasInterWindowInEdge(rcAid));
    };

    uint64_t trimmedVertexCount = 0;
    uint64_t trimmedWindowCount = 0;

    for(uint32_t w = 0; w < windowCount; w++) {
        const auto& window = anchorWindows[w];

        // Use filteredBackbonePositions if available, otherwise all positions.
        static thread_local vector<uint32_t> allPositions;
        const vector<uint32_t>* positionsPtr;
        if(!window.filteredBackbonePositions.empty()) {
            positionsPtr = &window.filteredBackbonePositions;
        } else {
            allPositions.clear();
            for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                allPositions.push_back(pos);
            }
            positionsPtr = &allPositions;
        }
        const auto& positions = *positionsPtr;
        if(positions.size() <= 1) continue;
        const auto journey = journeys[window.backboneOrientedReadId];

        // The trim criterion is intentionally direction-agnostic: it stops at
        // the first/last anchor carrying ANY inter-window edge, not just an
        // incoming (head) or outgoing (tail) one. An inverted repeat folds the
        // sequence back on itself, producing a fold-back connection that looks
        // like an outgoing edge at the head (or incoming at the tail). A
        // directional criterion would walk past such an edge and disable it,
        // destroying real inverted-repeat structure. Stopping at any
        // inter-window edge preserves it; residual strand-mix at the ends is
        // pruned later by coverage-aware tip removal. Do not make this
        // directional.

        // Head trim: remove anchors before the first one with any
        // inter-window edge (incoming or outgoing).
        uint64_t headTrim = 0;
        for(uint64_t i = 0; i < positions.size(); i++) {
            const uint64_t aid = uint64_t(journey[positions[i]]);
            if(hasIncomingInterWindowEdge(aid) || hasOutgoingInterWindowEdge(aid)) break;
            ++headTrim;
        }
        // No anchor has any inter-window edge — window is isolated.
        // Don't trim (leave for the isolated-window filter).
        if(headTrim >= positions.size()) continue;

        // Tail trim: remove anchors after the last one with any
        // inter-window edge (incoming or outgoing).
        uint64_t tailTrim = 0;
        for(int64_t i = int64_t(positions.size()) - 1; i > int64_t(headTrim); i--) {
            const uint64_t aid = uint64_t(journey[positions[uint64_t(i)]]);
            if(hasIncomingInterWindowEdge(aid) || hasOutgoingInterWindowEdge(aid)) break;
            ++tailTrim;
        }

        if(headTrim == 0 && tailTrim == 0) continue;
        ++trimmedWindowCount;

        for(uint64_t i = 0; i < headTrim; i++) {
            const uint64_t aid = uint64_t(journey[positions[i]]);
            disableAllEdges(aid);
            disableAllEdges(aid ^ 1ULL);
            ++trimmedVertexCount;
        }

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
}


uint64_t Shasta2AnchorGraph::removeIsolatedWindows(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Journeys& journeys)
{
    Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = num_vertices(anchorGraph);

    auto normalizeW = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // Find which normalized windows have active inter-window edges.
    std::set<uint32_t> connectedWindows;
    BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraphBaseClass) {
        if(!anchorGraph[e].useForAssembly) continue;
        const uint64_t src = uint64_t(source(e, anchorGraph));
        const uint64_t dst = uint64_t(target(e, anchorGraph));
        if(src >= anchorCount || dst >= anchorCount) continue;
        const uint32_t srcW = anchorToWindow[src];
        const uint32_t dstW = anchorToWindow[dst];
        if(srcW == noWindow || dstW == noWindow) continue;
        const uint32_t srcNorm = normalizeW(srcW);
        const uint32_t dstNorm = normalizeW(dstW);
        if(srcNorm != dstNorm) {
            connectedWindows.insert(srcNorm);
            connectedWindows.insert(dstNorm);
        }
    }

    // Disable all edges of isolated windows.
    uint64_t removedCount = 0;
    for(uint32_t w = 0; w < windowCount; w++) {
        if(connectedWindows.count(w)) continue;

        const auto& window = anchorWindows[w];
        const auto backboneJourney = journeys[window.backboneOrientedReadId];

        auto disableAnchor = [&](uint32_t pos) {
            const uint64_t aid = uint64_t(backboneJourney[pos]);
            if(aid < anchorCount) {
                for(auto oe = boost::out_edges(aid, anchorGraph); oe.first != oe.second; ++oe.first) {
                    if(anchorGraph[*oe.first].useForAssembly) disableEdge(*oe.first);
                }
                const uint64_t rcAid = aid ^ 1ULL;
                if(rcAid < anchorCount) {
                    for(auto oe = boost::out_edges(rcAid, anchorGraph); oe.first != oe.second; ++oe.first) {
                        if(anchorGraph[*oe.first].useForAssembly) disableEdge(*oe.first);
                    }
                }
            }
        };

        if(!window.filteredBackbonePositions.empty()) {
            for(const uint32_t pos : window.filteredBackbonePositions) {
                disableAnchor(pos);
            }
        } else {
            for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                disableAnchor(pos);
            }
        }
        ++removedCount;
    }

    cout << "removeIsolatedWindows: " << removedCount << " isolated windows removed ("
         << connectedWindows.size() << " connected)." << endl;
    return removedCount;
}


uint64_t Shasta2AnchorGraph::removeTipWindows(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Journeys& journeys,
    uint32_t maxTipWindows)
{
    Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = num_vertices(anchorGraph);

    auto normalizeW = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // Build window-level directed adjacency from active inter-window edges.
    // Uses RAW window IDs (forward and RC windows are separate nodes).
    // This is essential: with normalized IDs, the RC mirror of every
    // inter-window edge creates a reverse edge between the same normalized
    // pair, making every window appear to have both predecessors and
    // successors, so no tips would ever be detected.
    std::map<uint32_t, std::set<uint32_t>> windowPredecessors;
    std::map<uint32_t, std::set<uint32_t>> windowSuccessors;

    BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraphBaseClass) {
        if(!anchorGraph[e].useForAssembly) continue;
        const uint64_t src = uint64_t(source(e, anchorGraph));
        const uint64_t dst = uint64_t(target(e, anchorGraph));
        if(src >= anchorCount || dst >= anchorCount) continue;
        const uint32_t srcW = anchorToWindow[src];
        const uint32_t dstW = anchorToWindow[dst];
        if(srcW == noWindow || dstW == noWindow) continue;
        if(srcW == dstW) continue;
        if(normalizeW(srcW) == normalizeW(dstW)) continue;
        windowSuccessors[srcW].insert(dstW);
        windowPredecessors[dstW].insert(srcW);
    }

    // Helper to look up a raw window's neighbor set.
    static const std::set<uint32_t> emptySet;
    auto getPredecessors = [&](uint32_t ww) -> const std::set<uint32_t>& {
        auto it = windowPredecessors.find(ww);
        return (it != windowPredecessors.end()) ? it->second : emptySet;
    };
    auto getSuccessors = [&](uint32_t ww) -> const std::set<uint32_t>& {
        auto it = windowSuccessors.find(ww);
        return (it != windowSuccessors.end()) ? it->second : emptySet;
    };

    // Walk a tip chain from a dead-end raw window.
    // Returns the chain (raw window IDs) if it's a valid linear tip,
    // empty vector otherwise.
    auto walkTipChain = [&](uint32_t w) -> std::vector<uint32_t> {
        const bool hasPred = !getPredecessors(w).empty();
        const bool hasSucc = !getSuccessors(w).empty();

        // Must be dangling (one side only).
        if(hasPred == hasSucc) return {};

        std::vector<uint32_t> chain;
        chain.push_back(w);
        uint32_t current = w;

        while(chain.size() <= maxTipWindows) {
            const auto& neighbors = hasPred
                ? getPredecessors(current)
                : getSuccessors(current);

            if(neighbors.size() != 1) break;

            const uint32_t next = *neighbors.begin();

            // Cycle check.
            bool inChain = false;
            for(const uint32_t c : chain) {
                if(c == next) { inChain = true; break; }
            }
            if(inChain) break;

            // Linearity: next must have exactly one neighbor facing back.
            const auto& backNeighbors = hasPred
                ? getSuccessors(next)
                : getPredecessors(next);
            if(backNeighbors.size() != 1) break;

            chain.push_back(next);
            current = next;
        }

        if(chain.size() <= maxTipWindows) return chain;
        return {};
    };

    // Phase 1: Find all candidate tips and measure their chain length.
    // Store as (chainLength, startRawWindow) sorted shortest-first.
    std::vector<std::pair<uint32_t, uint32_t>> candidates;

    // Collect all raw window IDs that have any inter-window edges.
    std::set<uint32_t> allRawWindows;
    for(const auto& [w, s] : windowPredecessors) allRawWindows.insert(w);
    for(const auto& [w, s] : windowSuccessors) allRawWindows.insert(w);

    for(const uint32_t w : allRawWindows) {
        auto chain = walkTipChain(w);
        if(!chain.empty()) {
            candidates.push_back({(uint32_t)chain.size(), w});
        }
    }

    // Sort shortest-first so removing short tips can expose longer ones.
    std::sort(candidates.begin(), candidates.end());

    // Phase 2: Process candidates shortest-first.
    // Re-walk each tip to account for topology changes from prior removals.
    // Track removed normalized windows to avoid double-processing.
    std::set<uint32_t> normalizedWindowsRemoved;

    // Helper to disable all edges of a normalized window and update
    // the adjacency maps so subsequent re-walks see the new topology.
    auto removeWindow = [&](uint32_t normW) {
        if(!normalizedWindowsRemoved.insert(normW).second) return;

        const auto& window = anchorWindows[normW];
        const auto backboneJourney = journeys[window.backboneOrientedReadId];

        // Collect all anchor IDs in this window.
        std::vector<uint64_t> anchorIds;
        if(!window.filteredBackbonePositions.empty()) {
            for(const uint32_t pos : window.filteredBackbonePositions) {
                anchorIds.push_back(uint64_t(backboneJourney[pos]));
            }
        } else {
            for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                anchorIds.push_back(uint64_t(backboneJourney[pos]));
            }
        }

        // Disable all edges of each anchor (disableEdge handles RC mirrors).
        for(const uint64_t aid : anchorIds) {
            if(aid >= anchorCount) continue;
            for(auto oe = boost::out_edges(aid, anchorGraph); oe.first != oe.second; ++oe.first) {
                if(anchorGraph[*oe.first].useForAssembly)
                    disableEdge(*oe.first);
            }
            for(auto ie = boost::in_edges(aid, anchorGraph); ie.first != ie.second; ++ie.first) {
                if(anchorGraph[*ie.first].useForAssembly)
                    disableEdge(*ie.first);
            }
        }

        // Update adjacency maps: remove this window (both raw IDs)
        // from all neighbor sets.
        const uint32_t fwdW = normW;
        const uint32_t rcW = normW + windowCount;
        for(const uint32_t rawW : {fwdW, rcW}) {
            // Remove rawW from successors of its predecessors.
            for(const uint32_t pred : getPredecessors(rawW)) {
                auto it = windowSuccessors.find(pred);
                if(it != windowSuccessors.end()) it->second.erase(rawW);
            }
            // Remove rawW from predecessors of its successors.
            for(const uint32_t succ : getSuccessors(rawW)) {
                auto it = windowPredecessors.find(succ);
                if(it != windowPredecessors.end()) it->second.erase(rawW);
            }
            windowPredecessors.erase(rawW);
            windowSuccessors.erase(rawW);
        }
    };

    for(const auto& [chainLen, startW] : candidates) {
        const uint32_t startNorm = normalizeW(startW);

        // Skip if already removed.
        if(normalizedWindowsRemoved.count(startNorm)) continue;

        // Re-walk: topology may have changed from prior removals.
        auto chain = walkTipChain(startW);
        if(chain.empty()) continue;

        // Remove all windows in the chain.
        for(const uint32_t cw : chain) {
            removeWindow(normalizeW(cw));
        }
    }

    cout << "removeTipWindows (maxTipWindows=" << maxTipWindows << "): "
         << normalizedWindowsRemoved.size() << " windows removed." << endl;
    return normalizedWindowsRemoved.size();
}



uint64_t Shasta2AnchorGraph::popSuperbubbles(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Journeys& journeys,
    uint64_t maxSize)
{
    Shasta2AnchorGraph& anchorGraph = *this;
    const uint64_t anchorCount = num_vertices(anchorGraph);

    auto normalizeW = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // ---- Step 1: Build a window-level directed graph. ----
    // Nodes are normalized window IDs.
    // Only include edges from FORWARD raw windows (srcW < windowCount)
    // to avoid RC mirror edges creating reverse edges that form
    // universal 2-cycles and prevent superbubble detection.
    // The forward edges define the graph topology; the RC edges
    // are the same topology in reverse and are handled implicitly
    // by disableEdge's RC mirror support.
    std::map<uint32_t, std::set<uint32_t>> windowSuccessorsMap;

    BGL_FORALL_EDGES(e, anchorGraph, Shasta2AnchorGraphBaseClass) {
        if(!anchorGraph[e].useForAssembly) continue;
        const uint64_t src = uint64_t(source(e, anchorGraph));
        const uint64_t dst = uint64_t(target(e, anchorGraph));
        if(src >= anchorCount || dst >= anchorCount) continue;
        const uint32_t srcW = anchorToWindow[src];
        const uint32_t dstW = anchorToWindow[dst];
        if(srcW == noWindow || dstW == noWindow) continue;
        // Only forward-strand edges.
        if(srcW >= windowCount) continue;
        const uint32_t srcNorm = normalizeW(srcW);
        const uint32_t dstNorm = normalizeW(dstW);
        if(srcNorm == dstNorm) continue;
        windowSuccessorsMap[srcNorm].insert(dstNorm);
    }

    // Collect all normalized window IDs.
    std::set<uint32_t> allNormWindows;
    for(const auto& [w, succs] : windowSuccessorsMap) {
        allNormWindows.insert(w);
        for(const uint32_t s : succs) allNormWindows.insert(s);
    }

    // Map normalized window IDs to dense indices.
    std::vector<uint32_t> indexToNormWindow(allNormWindows.begin(), allNormWindows.end());
    std::map<uint32_t, uint32_t> normWindowToIndex;
    for(uint32_t i = 0; i < indexToNormWindow.size(); i++) {
        normWindowToIndex[indexToNormWindow[i]] = i;
    }

    // Build the window-level Boost graph.
    using WindowGraph = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::bidirectionalS>;
    const uint32_t nWindowNodes = uint32_t(indexToNormWindow.size());
    WindowGraph windowGraph(nWindowNodes);

    for(const auto& [srcW, succs] : windowSuccessorsMap) {
        const uint32_t srcIdx = normWindowToIndex[srcW];
        for(const uint32_t dstW : succs) {
            const uint32_t dstIdx = normWindowToIndex[dstW];
            boost::add_edge(srcIdx, dstIdx, windowGraph);
        }
    }

    cout << "popSuperbubbles: window graph has " << nWindowNodes
         << " nodes, " << num_edges(windowGraph) << " edges." << endl;

    // ---- Step 2: Find and pop superbubbles in the window graph. ----
    std::set<uint32_t> poppedWindowIndices;
    uint64_t poppedCount = 0;
    uint64_t candidateCount = 0;

    for(uint32_t vi = 0; vi < nWindowNodes; vi++) {
        if(out_degree(vi, windowGraph) < 2) continue;
        if(poppedWindowIndices.count(vi)) continue;
        ++candidateCount;

        // Diagnostic: log candidate details.
        {
            const uint32_t normW = indexToNormWindow[vi];
            cout << "  Candidate W" << normW
                 << " out-degree=" << out_degree(vi, windowGraph)
                 << " in-degree=" << in_degree(vi, windowGraph)
                 << " succs=[";
            bool first = true;
            BGL_FORALL_OUTEDGES_T(vi, e, windowGraph, WindowGraph) {
                if(!first) cout << ",";
                cout << "W" << indexToNormWindow[target(e, windowGraph)];
                first = false;
            }
            cout << "]" << endl;
        }

        const auto viExit = findSuperbubbleOnodera(windowGraph, vi, maxSize);
        if(viExit == WindowGraph::null_vertex()) {
            cout << "    -> no superbubble found." << endl;
            continue;
        }

        // Collect interior vertices via BFS.
        std::set<uint32_t> interiorIndices;
        {
            std::vector<uint32_t> queue;
            queue.push_back(vi);
            interiorIndices.insert(vi);
            while(!queue.empty()) {
                const uint32_t u = queue.back();
                queue.pop_back();
                if(u == viExit) continue;
                BGL_FORALL_OUTEDGES_T(u, e, windowGraph, WindowGraph) {
                    const uint32_t w = target(e, windowGraph);
                    if(interiorIndices.insert(w).second) {
                        queue.push_back(w);
                    }
                }
            }
        }

        // Skip if overlapping with a previously popped superbubble.
        bool overlap = false;
        for(const uint32_t u : interiorIndices) {
            if(poppedWindowIndices.count(u)) { overlap = true; break; }
        }
        if(overlap) continue;

        // Enumerate all paths from source to sink (as sequences of window indices).
        std::vector<std::vector<uint32_t>> allPaths;
        {
            struct DfsFrame {
                uint32_t vertex;
                std::vector<uint32_t> path;
            };
            std::vector<DfsFrame> stack;
            stack.push_back({vi, {vi}});

            while(!stack.empty()) {
                auto [current, path] = std::move(stack.back());
                stack.pop_back();

                if(current == viExit) {
                    allPaths.push_back(std::move(path));
                    continue;
                }

                BGL_FORALL_OUTEDGES_T(current, e, windowGraph, WindowGraph) {
                    const uint32_t w = target(e, windowGraph);
                    if(!interiorIndices.count(w)) continue;
                    auto newPath = path;
                    newPath.push_back(w);
                    stack.push_back({w, std::move(newPath)});
                }
            }
        }

        if(allPaths.size() < 2) continue;

        // Score each path. Prefer the path that follows the source
        // window's backbone chain (backboneNextWindow).
        // Walk the backbone chain from the source to build the set of
        // windows the backbone passes through.
        const uint32_t sourceNormW = indexToNormWindow[vi];
        const uint32_t sinkNormW = indexToNormWindow[viExit];

        std::set<uint32_t> backboneChainNormWindows;
        {
            uint32_t current = sourceNormW;
            backboneChainNormWindows.insert(current);
            // Walk backboneNextWindow until we reach the sink or a dead end.
            for(uint32_t step = 0; step < maxSize + 2; step++) {
                if(current >= anchorWindows.size()) break;
                const uint32_t next = anchorWindows[current].backboneNextWindow;
                if(next == AnchorWindowReadInterval::noWindow) break;
                const uint32_t nextNorm = normalizeW(next);
                backboneChainNormWindows.insert(nextNorm);
                if(nextNorm == sinkNormW) break;
                current = nextNorm;
            }
        }

        uint64_t bestPathIndex = 0;
        uint64_t bestScore = 0;

        for(uint64_t pi = 0; pi < allPaths.size(); pi++) {
            const auto& path = allPaths[pi];
            uint64_t score = 0;

            for(const uint32_t wIdx : path) {
                const uint32_t normW = indexToNormWindow[wIdx];

                // High score for windows on the backbone chain.
                if(backboneChainNormWindows.count(normW)) {
                    score += 1000;
                }

                // Tiebreak: total transition read support.
                if(normW < anchorWindows.size()) {
                    const auto& window = anchorWindows[normW];
                    for(const auto& [key, reads] : window.transitionReads) {
                        score += reads.size();
                    }
                }
            }

            if(score > bestScore) {
                bestScore = score;
                bestPathIndex = pi;
            }
        }

        // Collect the set of normalized window IDs on the best path.
        std::set<uint32_t> bestPathNormWindows;
        for(const uint32_t wIdx : allPaths[bestPathIndex]) {
            bestPathNormWindows.insert(indexToNormWindow[wIdx]);
        }

        // Collect normalized window IDs on ALL paths.
        std::set<uint32_t> allPathNormWindows;
        for(const auto& path : allPaths) {
            for(const uint32_t wIdx : path) {
                allPathNormWindows.insert(indexToNormWindow[wIdx]);
            }
        }

        // Windows to disable: in the superbubble but NOT on the best path.
        // Don't disable the source or sink — they're shared by all paths.
        std::set<uint32_t> windowsToDisable;
        for(const uint32_t normW : allPathNormWindows) {
            if(bestPathNormWindows.count(normW)) continue;
            if(normW == sourceNormW) continue;
            if(normW == sinkNormW) continue;
            windowsToDisable.insert(normW);
        }

        if(windowsToDisable.empty()) continue;

        // Disable all anchor edges of the non-chosen windows.
        // This is safe because the superbubble property guarantees
        // interior vertices have no connections outside the bubble
        // except through source and sink.
        uint64_t disabledEdgeCount = 0;
        for(const uint32_t normW : windowsToDisable) {
            if(normW >= anchorWindows.size()) continue;
            const auto& window = anchorWindows[normW];
            const auto backboneJourney = journeys[window.backboneOrientedReadId];

            vector<uint32_t> positions;
            if(!window.filteredBackbonePositions.empty()) {
                positions = window.filteredBackbonePositions;
            } else {
                for(uint32_t p = window.backboneBegin; p < window.backboneEnd; p++) {
                    positions.push_back(p);
                }
            }

            for(const uint32_t pos : positions) {
                const uint64_t aid = uint64_t(backboneJourney[pos]);
                if(aid >= anchorCount) continue;

                for(auto oe = boost::out_edges(aid, anchorGraph); oe.first != oe.second; ++oe.first) {
                    if(anchorGraph[*oe.first].useForAssembly) {
                        disableEdge(*oe.first);
                        ++disabledEdgeCount;
                    }
                }
                for(auto ie = boost::in_edges(aid, anchorGraph); ie.first != ie.second; ++ie.first) {
                    if(anchorGraph[*ie.first].useForAssembly) {
                        disableEdge(*ie.first);
                        ++disabledEdgeCount;
                    }
                }
            }
        }

        // Remove disabled windows from the window graph so subsequent
        // superbubble detections see the correct topology.
        for(const uint32_t normW : windowsToDisable) {
            auto it = normWindowToIndex.find(normW);
            if(it == normWindowToIndex.end()) continue;
            boost::clear_vertex(it->second, windowGraph);
        }

        for(const uint32_t u : interiorIndices) {
            poppedWindowIndices.insert(u);
        }
        ++poppedCount;

        cout << "  Superbubble popped: source=W" << sourceNormW
             << " sink=W" << sinkNormW
             << " paths=" << allPaths.size()
             << " best=[";
        for(uint64_t i = 0; i < allPaths[bestPathIndex].size(); i++) {
            if(i > 0) cout << ",";
            cout << "W" << indexToNormWindow[allPaths[bestPathIndex][i]];
        }
        cout << "] disabled=" << windowsToDisable.size()
             << " windows, " << disabledEdgeCount << " edges." << endl;
    }

    cout << "popSuperbubbles: " << candidateCount << " candidates (out-degree >= 2), "
         << poppedCount << " superbubbles popped." << endl;
    return poppedCount;
}

