// Shasta.
#include "Shasta2AnchorGraph.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2AnchorPair.hpp"
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
            DINARA_ASSERT(anchorPair.size() == counts[i]);
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
    const Reads* reads) :
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
    const uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    const uint32_t windowCount = uint32_t(anchorWindows.size());
    vector<uint32_t> anchorToWindow(anchorCount, noWindow);
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

        // Use filtered backbone positions if available.
        const auto& positions = window.filteredBackbonePositions;
        if(!positions.empty()) {
            for(uint64_t i = 0; i + 1 < positions.size(); i++) {
                const Shasta2AnchorId anchorIdA = backboneJourney[positions[i]];
                const Shasta2AnchorId anchorIdB = backboneJourney[positions[i + 1]];
                addEdgeIfValid(anchorIdA, anchorIdB);
                // RC mirror edge.
                const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(anchorIdA) ^ 1ULL);
                const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(anchorIdB) ^ 1ULL);
                if(uint64_t(rcA) < anchorCount && uint64_t(rcB) < anchorCount) {
                    addEdgeIfValid(rcB, rcA);
                }
            }
        } else {
            // Fallback: all consecutive positions.
            for(uint32_t pos = window.backboneBegin; pos + 1 < window.backboneEnd; pos++) {
                const Shasta2AnchorId anchorIdA = backboneJourney[pos];
                const Shasta2AnchorId anchorIdB = backboneJourney[pos + 1];
                addEdgeIfValid(anchorIdA, anchorIdB);
                const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(anchorIdA) ^ 1ULL);
                const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(anchorIdB) ^ 1ULL);
                if(uint64_t(rcA) < anchorCount && uint64_t(rcB) < anchorCount) {
                    addEdgeIfValid(rcB, rcA);
                }
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

    const uint64_t intraEdgeCount = num_edges(anchorGraph) - alternatePathEdgeCount;

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

                // Track backbone read's transition.
                if(oid == w.backboneOrientedReadId) {
                    w.backbonePreviousWindow = prev;
                    w.backboneNextWindow = next;
                }
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

    // For each window pair, pick the candidate with the most shared reads.
    uint64_t interWindowZeroPairs = 0;
    uint64_t interWindowBelowCoverage = 0;
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
        Shasta2AnchorPair bestPair;
        uint64_t bestSize = 0;
        for(const auto& [apk, count] : candidates) {
            Shasta2AnchorPair anchorPair(anchors, apk.anchorIdA, apk.anchorIdB, false);
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
         << interWindowZeroPairs << " rejected (zero forward-flow reads), "
         << interWindowBelowCoverage << " rejected (below minInterWindowCoverage="
         << minInterWindowCoverage << ")." << endl;

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
                 << " (backbone: ";
            if(window.backbonePreviousWindow == noW) cout << "-"; else cout << window.backbonePreviousWindow;
            cout << "→" << w << "→";
            if(window.backboneNextWindow == noW) cout << "-"; else cout << window.backboneNextWindow;
            cout << ", " << totalTransitionReads << " reads): incoming=[";
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

    // Trim dangling backbone tails at unambiguous 1-to-1 window connections.
    // If window A connects to only one next window B, and window B has only
    // one incoming window A, then the backbone chain of A past the inter-window
    // edge anchor and the backbone chain of B before the inter-window edge
    // anchor are dangling tails that should be removed.
    {
        auto normalize = [&](uint32_t w) -> uint32_t {
            return (w >= windowCount) ? (w - windowCount) : w;
        };
        const uint32_t noW = AnchorWindow::noWindow;

        // Build per-window outgoing/incoming sets from transitionReads (through-flows only).
        std::map<uint32_t, std::set<uint32_t>> outgoing, incoming;
        for(uint32_t w = 0; w < windowCount; w++) {
            const auto& window = anchorWindows[w];
            for(const auto& [key, reads] : window.transitionReads) {
                if(key.second != noW) {
                    outgoing[w].insert(key.second);
                }
                if(key.first != noW) {
                    incoming[w].insert(key.first);
                }
            }
        }

        uint64_t trimmedEdgeCount = 0;
        for(const auto& edgeInfo : createdEdges) {
            const uint32_t srcW = normalize(edgeInfo.windowPair.first);
            const uint32_t dstW = normalize(edgeInfo.windowPair.second);

            // Check 1-to-1: source has exactly 1 outgoing, destination has exactly 1 incoming.
            if(outgoing[srcW].size() != 1 || incoming[dstW].size() != 1) continue;
            if(*outgoing[srcW].begin() != dstW || *incoming[dstW].begin() != srcW) continue;

            // Find the inter-window edge anchors in the backbone chains.
            const auto& srcWindow = anchorWindows[srcW];
            const auto& dstWindow = anchorWindows[dstW];
            const auto srcJourney = journeys[srcWindow.backboneOrientedReadId];
            const auto dstJourney = journeys[dstWindow.backboneOrientedReadId];

            const auto& srcPositions = srcWindow.filteredBackbonePositions;
            const auto& dstPositions = dstWindow.filteredBackbonePositions;

            // Find position of anchorIdA in source backbone.
            // Trim all edges after anchorIdA (the tail past the connection).
            int64_t srcIdx = -1;
            if(!srcPositions.empty()) {
                for(uint64_t i = 0; i < srcPositions.size(); i++) {
                    if(srcJourney[srcPositions[i]] == edgeInfo.anchorIdA) {
                        srcIdx = int64_t(i);
                        break;
                    }
                }
                if(srcIdx >= 0) {
                    for(uint64_t i = uint64_t(srcIdx) + 1; i < srcPositions.size(); i++) {
                        const Shasta2AnchorId aid = srcJourney[srcPositions[i]];
                        // Remove all out-edges and in-edges of this vertex
                        // that are intra-window (not inter-window).
                        boost::clear_vertex(uint64_t(aid), anchorGraph);
                        // Also clear the RC mirror vertex.
                        const uint64_t rcAid = uint64_t(aid) ^ 1ULL;
                        if(rcAid < anchorCount) {
                            boost::clear_vertex(rcAid, anchorGraph);
                        }
                        ++trimmedEdgeCount;
                    }
                }
            }

            // Find position of anchorIdB in destination backbone.
            // Trim all edges before anchorIdB (the head before the connection).
            int64_t dstIdx = -1;
            if(!dstPositions.empty()) {
                for(uint64_t i = 0; i < dstPositions.size(); i++) {
                    if(dstJourney[dstPositions[i]] == edgeInfo.anchorIdB) {
                        dstIdx = int64_t(i);
                        break;
                    }
                }
                if(dstIdx > 0) {
                    for(int64_t i = 0; i < dstIdx; i++) {
                        const Shasta2AnchorId aid = dstJourney[dstPositions[i]];
                        boost::clear_vertex(uint64_t(aid), anchorGraph);
                        const uint64_t rcAid = uint64_t(aid) ^ 1ULL;
                        if(rcAid < anchorCount) {
                            boost::clear_vertex(rcAid, anchorGraph);
                        }
                        ++trimmedEdgeCount;
                    }
                }
            }
        }

        if(trimmedEdgeCount > 0) {
            cout << "Trimmed " << trimmedEdgeCount
                 << " dangling backbone vertices at 1-to-1 window connections." << endl;
        }
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

    const uint64_t interEdgeCount =
        num_edges(anchorGraph) - intraEdgeCount - alternatePathEdgeCount;

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
         << " (" << intraEdgeCount << " intra-window, "
         << alternatePathEdgeCount << " alternate-path, "
         << interEdgeCount << " inter-window)." << endl;
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




