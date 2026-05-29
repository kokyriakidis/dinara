// Detangling: find all bypass candidates across all tangled windows,
// sort by common read count (strongest first), process in order.
//
// For each candidate (prevW -> X -> nextW):
// 1. Compute common reads between prevW's exit anchor and nextW's entry anchor
// 2. Create a bypass edge connecting them directly (+ RC mirror)
// 3. Remove the common reads from every backbone anchor of window X
//
// Processing strongest candidates first ensures the best connections
// are made with the most reads available.

#include "DinaraDetangle.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <set>

using namespace dinara;
using std::cout;
using std::endl;
using std::map;
using std::pair;
using std::set;
using std::vector;



uint64_t dinara::detangleWindows(
    Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const vector<AnchorWindow>& anchorWindows,
    uint64_t minFlowCoverage,
    vector<DetangleBypassEdge>& bypassEdges)
{
    cout << timestamp << "Detangling begins." << endl;

    const uint32_t noW = AnchorWindow::noWindow;
    bypassEdges.clear();

    // Copy all anchor data to a mutable vector so we can modify in place.
    const uint64_t anchorCount = anchors.anchorMarkerInfos.size();
    vector<vector<Shasta2AnchorMarkerInfo>> mutableAnchors(anchorCount);
    for(uint64_t i = 0; i < anchorCount; i++) {
        const auto span = anchors.anchorMarkerInfos[i];
        mutableAnchors[i].assign(span.begin(), span.end());
    }

    // Collect all bypass candidates from all tangled windows.
    // For each tangled window B with triplet (A -> B -> C), the candidate's
    // read count is the number of common reads between windows A and C
    // (from transitionReads).
    struct BypassCandidate {
        uint32_t windowId;       // The tangled window being bypassed (B).
        uint32_t prevW;          // Previous window (A).
        uint32_t nextW;          // Next window (C).
        Shasta2AnchorId prevExitAnchor;   // Last anchor in prev window.
        Shasta2AnchorId nextEntryAnchor;  // First anchor in next window.
        set<uint32_t> flowReads; // Common reads between A and C through B.
    };
    vector<BypassCandidate> candidates;

    for(uint32_t wIdx = 0; wIdx < anchorWindows.size(); wIdx++) {
        const AnchorWindow& window = anchorWindows[wIdx];

        // Find flows: (prev, next) pairs from transitionReads.
        // Merge directional pairs (A->B) and (B->A).
        using FlowKey = pair<uint32_t, uint32_t>;
        map<FlowKey, set<uint32_t>> mergedFlows; // canonical key -> read set

        for(const auto& [key, reads] : window.transitionReads) {
            if(key.first == noW || key.second == noW) continue;
            FlowKey canonical = {std::min(key.first, key.second),
                                 std::max(key.first, key.second)};
            for(const OrientedReadId& oid : reads) {
                mergedFlows[canonical].insert(oid.getValue());
            }
        }

        // Only tangled windows (≥2 distinct flows).
        if(mergedFlows.size() < 2) continue;

        // For each directional key in transitionReads, create a candidate.
        for(const auto& [key, reads] : window.transitionReads) {
            if(key.first == noW || key.second == noW) continue;

            Shasta2AnchorId prevExitAnchor = Shasta2AnchorId(0);
            bool foundPrev = false;
            for(const auto& ie : window.inEdges) {
                if(ie.otherWindow == key.first) {
                    prevExitAnchor = ie.anchorIdA;
                    foundPrev = true;
                    break;
                }
            }

            Shasta2AnchorId nextEntryAnchor = Shasta2AnchorId(0);
            bool foundNext = false;
            for(const auto& oe : window.outEdges) {
                if(oe.otherWindow == key.second) {
                    nextEntryAnchor = oe.anchorIdB;
                    foundNext = true;
                    break;
                }
            }

            if(foundPrev && foundNext) {
                // The flow reads are the merged reads for this canonical pair.
                FlowKey canonical = {std::min(key.first, key.second),
                                     std::max(key.first, key.second)};
                candidates.push_back({wIdx, key.first, key.second,
                                      prevExitAnchor, nextEntryAnchor,
                                      mergedFlows[canonical]});
            }
        }
    }

    // Sort candidates by flow read count (descending) — strongest first.
    std::sort(candidates.begin(), candidates.end(),
        [](const BypassCandidate& a, const BypassCandidate& b) {
            return a.flowReads.size() > b.flowReads.size();
        });

    // Process candidates in order of decreasing common reads.
    set<uint32_t> detangledWindows;
    uint64_t processedCount = 0;

    for(const auto& cand : candidates) {
        // Recompute common reads between the connection anchors
        // using the mutable copy (reflects prior removals).
        // Only keep reads that are still on both connection anchors
        // AND are in the original flow read set.
        set<uint32_t> prevReads, nextReads;
        for(const auto& mi : mutableAnchors[uint64_t(cand.prevExitAnchor)]) {
            prevReads.insert(mi.orientedReadId.getValue());
        }
        for(const auto& mi : mutableAnchors[uint64_t(cand.nextEntryAnchor)]) {
            nextReads.insert(mi.orientedReadId.getValue());
        }
        set<uint32_t> commonReads;
        for(const uint32_t r : cand.flowReads) {
            if(prevReads.count(r) && nextReads.count(r)) {
                commonReads.insert(r);
            }
        }

        if(commonReads.empty()) {
            continue;
        }

        // Create forward bypass edge.
        bypassEdges.push_back({cand.prevExitAnchor, cand.nextEntryAnchor});

        // Create RC bypass edge.
        const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(cand.prevExitAnchor) ^ 1ULL);
        const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(cand.nextEntryAnchor) ^ 1ULL);
        bypassEdges.push_back({rcB, rcA});

        cout << "  Bypass: window " << cand.windowId
             << " " << cand.prevExitAnchor << " -> " << cand.nextEntryAnchor
             << " (+ RC " << rcB << " -> " << rcA << ")"
             << " prev=" << cand.prevW << " next=" << cand.nextW
             << " commonReads=" << commonReads.size() << endl;

        // Remove common reads from every backbone anchor of the tangled window.
        const AnchorWindow& window = anchorWindows[cand.windowId];
        const auto backboneJourney = journeys[window.backboneOrientedReadId];
        const auto& positions = window.filteredBackbonePositions;

        set<Shasta2AnchorId> processedAnchors;
        for(const uint32_t pos : positions) {
            const Shasta2AnchorId anchorId = backboneJourney[pos];
            const Shasta2AnchorId canonicalId = anchorId & ~Shasta2AnchorId(1);
            const Shasta2AnchorId rcId = canonicalId | Shasta2AnchorId(1);

            if(processedAnchors.count(canonicalId)) continue;
            processedAnchors.insert(canonicalId);

            // Remove from canonical anchor.
            auto& canonicalData = mutableAnchors[uint64_t(canonicalId)];
            canonicalData.erase(
                std::remove_if(canonicalData.begin(), canonicalData.end(),
                    [&](const Shasta2AnchorMarkerInfo& mi) {
                        return commonReads.count(mi.orientedReadId.getValue());
                    }),
                canonicalData.end());

            // Remove from RC anchor (flip read IDs).
            auto& rcData = mutableAnchors[uint64_t(rcId)];
            rcData.erase(
                std::remove_if(rcData.begin(), rcData.end(),
                    [&](const Shasta2AnchorMarkerInfo& mi) {
                        return commonReads.count(mi.orientedReadId.getValue() ^ 1);
                    }),
                rcData.end());
        }

        detangledWindows.insert(cand.windowId);
        ++processedCount;
    }

    // Rebuild anchorMarkerInfos from the modified mutable copy.
    anchors.anchorMarkerInfos.clear();
    for(const auto& anchorData : mutableAnchors) {
        anchors.anchorMarkerInfos.appendVector(anchorData);
    }

    cout << timestamp << "Detangling: " << detangledWindows.size()
         << " windows detangled, "
         << processedCount << " bypasses processed, "
         << bypassEdges.size() << " bypass edges created." << endl;

    return detangledWindows.size();
}
