// ============================================================================
// Detangle: window bypass
// ============================================================================
//
// For each window X with through-flows (A -> X -> B), creates bypass
// edges connecting A directly to B, skipping X.
//
// Algorithm:
//   1. Build a local anchorToWindow map from backbone positions.
//   2. For each window X, iterate transitionReads to find (A, B) pairs
//      where reads flow A -> X -> B (both A and B are real windows).
//   3. For each (A, B) pair, walk each flow read's journey to find the
//      last anchor in A and the first anchor in B. Tally (anchorA, anchorB)
//      pairs and pick the one with the most reads as the bypass edge.
//   4. Create bypass edge (forward + RC mirror).
//   5. Remove flow reads from X's backbone anchors.
//
// No minimum number of distinct flows is required — even a single
// A -> X -> B flow is bypassed if it meets minFlowCoverage.
//
// Processing order: candidates sorted by flow read count (descending)
// so stronger flows get their full read sets first.
//
// RC handling: for each forward bypass edge (anchorA -> anchorB),
// an RC mirror (anchorB' -> anchorA') is also created. Read removals
// update both canonical and RC anchors.
// ============================================================================

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
    cout << timestamp << "Detangle (window bypass) begins." << endl;

    const uint32_t noW = AnchorWindow::noWindow;
    const uint32_t windowCount = uint32_t(anchorWindows.size());
    bypassEdges.clear();

    // ---- Step 1: Build local anchorToWindow map. ----
    const uint64_t anchorCount = anchors.anchorMarkerInfos.size();
    vector<uint32_t> anchorToWindow(anchorCount, noW);
    for(uint32_t windowId = 0; windowId < windowCount; windowId++) {
        const AnchorWindow& window = anchorWindows[windowId];
        const auto backboneJourney = journeys[window.backboneOrientedReadId];

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
            const uint64_t aid = uint64_t(backboneJourney[pos]);
            if(aid < anchorCount) {
                anchorToWindow[aid] = windowId;
            }
            const uint64_t rcAid = aid ^ 1ULL;
            if(rcAid < anchorCount) {
                anchorToWindow[rcAid] = windowId + windowCount;
            }
        }
    }

    auto normalize = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // ---- Step 2: Copy anchor data to a mutable structure. ----
    vector<vector<Shasta2AnchorMarkerInfo>> mutableAnchors(anchorCount);
    for(uint64_t i = 0; i < anchorCount; i++) {
        const auto span = anchors.anchorMarkerInfos[i];
        mutableAnchors[i].assign(span.begin(), span.end());
    }

    // ---- Step 3: Collect bypass candidates. ----
    // For each window X and each (A, B) flow through it, walk the flow
    // reads' journeys to find the best anchor pair between A and B.
    struct BypassCandidate {
        uint32_t windowId;                // The window being bypassed (X).
        uint32_t prevW;                   // Window A.
        uint32_t nextW;                   // Window B.
        Shasta2AnchorId bestAnchorA;      // Last anchor in A (best pair).
        Shasta2AnchorId bestAnchorB;      // First anchor in B (best pair).
        vector<uint32_t> flowReadIds;     // Oriented read IDs for this flow.
    };
    vector<BypassCandidate> candidates;

    for(uint32_t wIdx = 0; wIdx < windowCount; wIdx++) {
        const AnchorWindow& window = anchorWindows[wIdx];

        for(const auto& [key, reads] : window.transitionReads) {
            if(key.first == noW || key.second == noW) continue;
            if(reads.size() < minFlowCoverage) continue;

            const uint32_t prevW = key.first;
            const uint32_t nextW = key.second;

            // Walk each flow read's journey to find the last anchor in
            // prevW and the first anchor in nextW. Tally anchor pairs.
            map<pair<uint64_t, uint64_t>, uint64_t> anchorPairCounts;
            map<pair<uint64_t, uint64_t>, vector<uint32_t>> anchorPairReads;

            for(const OrientedReadId& oid : reads) {
                const auto journey = journeys[oid];
                if(journey.empty()) continue;

                // Find the span of window X (wIdx) in the journey,
                // then pick lastInPrev before that span and firstInNext
                // after it. This handles prevW == nextW correctly.

                // Find the first and last journey positions belonging to
                // the bypassed window X.
                uint32_t xFirst = uint32_t(journey.size());
                uint32_t xLast = 0;
                for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
                    const Shasta2AnchorId anchorId = journey[pos];
                    if(uint64_t(anchorId) >= anchorCount) continue;
                    const uint32_t aw = anchorToWindow[uint64_t(anchorId)];
                    if(aw == noW) continue;
                    if(normalize(aw) == wIdx) {
                        if(xFirst == uint32_t(journey.size())) xFirst = pos;
                        xLast = pos;
                    }
                }
                // Read must actually touch window X.
                if(xFirst == uint32_t(journey.size())) continue;

                // lastInPrev: last anchor in prevW before the X span.
                Shasta2AnchorId lastInPrev = Shasta2AnchorId(0);
                bool foundPrev = false;
                for(uint32_t pos = 0; pos < xFirst; pos++) {
                    const Shasta2AnchorId anchorId = journey[pos];
                    if(uint64_t(anchorId) >= anchorCount) continue;
                    const uint32_t aw = anchorToWindow[uint64_t(anchorId)];
                    if(aw == noW) continue;
                    if(normalize(aw) == prevW) {
                        lastInPrev = anchorId;
                        foundPrev = true;
                    }
                }

                // firstInNext: first anchor in nextW after the X span.
                Shasta2AnchorId firstInNext = Shasta2AnchorId(0);
                bool foundNext = false;
                for(uint32_t pos = xLast + 1; pos < uint32_t(journey.size()); pos++) {
                    const Shasta2AnchorId anchorId = journey[pos];
                    if(uint64_t(anchorId) >= anchorCount) continue;
                    const uint32_t aw = anchorToWindow[uint64_t(anchorId)];
                    if(aw == noW) continue;
                    if(normalize(aw) == nextW) {
                        firstInNext = anchorId;
                        foundNext = true;
                        break;
                    }
                }

                if(foundPrev && foundNext) {
                    auto pairKey = pair<uint64_t, uint64_t>(
                        uint64_t(lastInPrev), uint64_t(firstInNext));
                    anchorPairCounts[pairKey]++;
                    anchorPairReads[pairKey].push_back(oid.getValue());
                }
            }

            if(anchorPairCounts.empty()) continue;

            // Pick the anchor pair with the most reads.
            auto bestIt = std::max_element(
                anchorPairCounts.begin(), anchorPairCounts.end(),
                [](const auto& a, const auto& b) {
                    return a.second < b.second;
                });

            if(bestIt->second < minFlowCoverage) continue;

            candidates.push_back({
                wIdx, prevW, nextW,
                Shasta2AnchorId(bestIt->first.first),
                Shasta2AnchorId(bestIt->first.second),
                anchorPairReads[bestIt->first]
            });
        }
    }

    // ---- Step 4: Sort candidates by flow read count (descending). ----
    std::sort(candidates.begin(), candidates.end(),
        [](const BypassCandidate& a, const BypassCandidate& b) {
            return a.flowReadIds.size() > b.flowReadIds.size();
        });

    // ---- Step 5: Process each candidate. ----
    set<uint32_t> detangledWindows;
    uint64_t processedCount = 0;

    for(const auto& cand : candidates) {

        // Recompute: only keep reads still present on both anchors.
        set<uint32_t> anchorAReads, anchorBReads;
        for(const auto& mi : mutableAnchors[uint64_t(cand.bestAnchorA)]) {
            anchorAReads.insert(mi.orientedReadId.getValue());
        }
        for(const auto& mi : mutableAnchors[uint64_t(cand.bestAnchorB)]) {
            anchorBReads.insert(mi.orientedReadId.getValue());
        }
        set<uint32_t> commonReads;
        for(const uint32_t r : cand.flowReadIds) {
            if(anchorAReads.count(r) && anchorBReads.count(r)) {
                commonReads.insert(r);
            }
        }

        if(commonReads.size() < minFlowCoverage) continue;

        // Create bypass edges (forward + RC mirror).
        bypassEdges.push_back({cand.bestAnchorA, cand.bestAnchorB});

        const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(cand.bestAnchorA) ^ 1ULL);
        const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(cand.bestAnchorB) ^ 1ULL);
        bypassEdges.push_back({rcB, rcA});

        cout << "  Bypass: window " << cand.windowId
             << " anchor " << cand.bestAnchorA << " -> " << cand.bestAnchorB
             << " (+ RC " << rcB << " -> " << rcA << ")"
             << " prev=" << cand.prevW << " next=" << cand.nextW
             << " reads=" << commonReads.size() << endl;

        // Remove flow reads from backbone anchors of the bypassed window.
        const AnchorWindow& window = anchorWindows[cand.windowId];
        const auto backboneJourney = journeys[window.backboneOrientedReadId];

        // Use filtered positions if available, otherwise all positions.
        vector<uint32_t> allPositions;
        if(window.filteredBackbonePositions.empty()) {
            for(uint32_t p = window.backboneBegin; p < window.backboneEnd; p++) {
                allPositions.push_back(p);
            }
        }
        const auto& positions = window.filteredBackbonePositions.empty()
            ? allPositions : window.filteredBackbonePositions;

        set<Shasta2AnchorId> processedAnchors;
        for(const uint32_t pos : positions) {
            const Shasta2AnchorId anchorId = backboneJourney[pos];
            const Shasta2AnchorId canonicalId = anchorId & ~Shasta2AnchorId(1);
            const Shasta2AnchorId rcId = canonicalId | Shasta2AnchorId(1);

            if(processedAnchors.count(canonicalId)) continue;
            processedAnchors.insert(canonicalId);

            // Remove flow reads from both the canonical and RC anchor.
            // commonReads contains oriented read IDs which may be on
            // either strand, so check both oid and oid^1.
            auto isFlowRead = [&](const Shasta2AnchorMarkerInfo& mi) {
                const uint32_t oid = mi.orientedReadId.getValue();
                return commonReads.count(oid) || commonReads.count(oid ^ 1);
            };

            auto& canonicalData = mutableAnchors[uint64_t(canonicalId)];
            canonicalData.erase(
                std::remove_if(canonicalData.begin(), canonicalData.end(), isFlowRead),
                canonicalData.end());

            auto& rcData = mutableAnchors[uint64_t(rcId)];
            rcData.erase(
                std::remove_if(rcData.begin(), rcData.end(), isFlowRead),
                rcData.end());
        }

        detangledWindows.insert(cand.windowId);
        ++processedCount;
    }

    // ---- Step 6: Rebuild anchorMarkerInfos from the modified copy. ----
    anchors.anchorMarkerInfos.clear();
    for(const auto& anchorData : mutableAnchors) {
        anchors.anchorMarkerInfos.appendVector(anchorData);
    }

    cout << timestamp << "Detangle: " << detangledWindows.size()
         << " windows bypassed, "
         << processedCount << " bypasses processed, "
         << bypassEdges.size() << " bypass edges created." << endl;

    return detangledWindows.size();
}
