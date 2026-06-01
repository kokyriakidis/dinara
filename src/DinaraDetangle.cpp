// ============================================================================
// Detangle: tangle-aware window bypass
// ============================================================================
//
// For each window X that is a collapsed repeat (multiple predecessors
// AND/OR multiple successors with clean diagonal pairing in the tangle
// matrix), creates bypass edges connecting predecessor windows directly
// to successor windows, skipping X.
//
// Linear windows (≤1 predecessor AND ≤1 successor) are never bypassed.
//
// Algorithm:
//   1. Build a local anchorToWindow map from backbone positions.
//   2. For each window X, collect predecessors/successors from
//      transitionReads. Skip linear windows.
//   3. Build N×M tangle matrix (rows=predecessors, cols=successors,
//      cells=read count). Check if cleanly separable.
//   4. For each dominant (prev, succ) pairing, walk flow reads'
//      journeys to find the best anchor pair and create a bypass edge.
//   5. Remove flow reads from X's backbone anchors.
//
// Tangle resolution uses Verkko-style triplet counting: a triplet
// (pred, X, succ) is "solid" if read count >= minEdgeSupport. Every
// significant edge (coverage >= minEdgeCoverage) must be covered by
// at least one solid triplet for the window to be resolvable.
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
    vector<DetangleBypassEdge>& bypassEdges)
{
    cout << timestamp << "Detangle (tangle-aware window bypass) begins." << endl;

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

    // Verkko-style triplet parameters (matching Verkko naming/defaults).
    const vector<uint64_t> resolveSteps = {20, 10, 5};  // Iterated min_edge_support values.
    const uint64_t minEdgeCoverage = 5;  // Min total reads on an edge for it to be "significant".

    // ---- Step 2b: Compute per-window coverage (average reads per anchor). ----
    // A window is "removable" if its coverage < minEdgeCoverage (matching Verkko).
    vector<double> windowCoverage(windowCount, 0.);
    for(uint32_t windowId = 0; windowId < windowCount; windowId++) {
        const AnchorWindow& window = anchorWindows[windowId];
        const auto backboneJourney = journeys[window.backboneOrientedReadId];

        vector<uint32_t> allPositions;
        if(window.filteredBackbonePositions.empty()) {
            for(uint32_t p = window.backboneBegin; p < window.backboneEnd; p++) {
                allPositions.push_back(p);
            }
        }
        const auto& positions = window.filteredBackbonePositions.empty()
            ? allPositions : window.filteredBackbonePositions;

        if(positions.empty()) continue;
        uint64_t totalReads = 0;
        for(const uint32_t pos : positions) {
            const uint64_t aid = uint64_t(backboneJourney[pos]);
            if(aid < anchorCount) {
                totalReads += mutableAnchors[aid].size();
            }
        }
        windowCoverage[windowId] = double(totalReads) / double(positions.size());
    }

    auto isRemovable = [&](uint32_t w) -> bool {
        const uint32_t normW = normalize(w);
        if(normW >= windowCount) return true;
        return windowCoverage[normW] < double(minEdgeCoverage);
    };

    // ---- Steps 3-5: Iterative triplet resolution (matching Verkko). ----
    // Resolve easy tangles first (high minEdgeSupport), then progressively
    // harder ones. Each round finds candidates, creates bypass edges, and
    // removes flow reads before the next round.
    struct BypassCandidate {
        uint32_t windowId;
        uint32_t prevW;
        uint32_t nextW;
        Shasta2AnchorId bestAnchorA;
        Shasta2AnchorId bestAnchorB;
        vector<uint32_t> flowReadIds;
    };

    set<uint32_t> detangledWindows;
    uint64_t totalProcessed = 0;

    for(const uint64_t minEdgeSupport : resolveSteps) {
    cout << timestamp << "Resolve step: minEdgeSupport=" << minEdgeSupport << endl;

    vector<BypassCandidate> candidates;
    uint64_t linearSkipped = 0;
    uint64_t tripletFailed = 0;

    // Helper: walk flow reads' journeys to find anchor pairs for a
    // specific (prevW, nextW) pairing through window wIdx.
    auto findAnchorPairs = [&](
        uint32_t wIdx, uint32_t prevW, uint32_t nextW,
        const vector<OrientedReadId>& reads)
        -> BypassCandidate
    {
        map<pair<uint64_t, uint64_t>, uint64_t> anchorPairCounts;
        map<pair<uint64_t, uint64_t>, vector<uint32_t>> anchorPairReads;

        for(const OrientedReadId& oid : reads) {
            const auto journey = journeys[oid];
            if(journey.empty()) continue;

            const uint32_t jLen = uint32_t(journey.size());
            vector<uint32_t> posWindow(jLen, noW);
            for(uint32_t pos = 0; pos < jLen; pos++) {
                const Shasta2AnchorId anchorId = journey[pos];
                if(uint64_t(anchorId) >= anchorCount) continue;
                const uint32_t aw = anchorToWindow[uint64_t(anchorId)];
                if(aw == noW) continue;
                posWindow[pos] = normalize(aw);
            }

            uint32_t pos = 0;
            while(pos < jLen) {
                if(posWindow[pos] != wIdx) { pos++; continue; }
                const uint32_t runStart = pos;
                while(pos < jLen && (posWindow[pos] == wIdx || posWindow[pos] == noW)) pos++;

                uint32_t immPrev = noW;
                for(uint32_t p = runStart; p > 0; ) {
                    p--;
                    if(posWindow[p] != noW) { immPrev = posWindow[p]; break; }
                }
                if(immPrev != prevW) continue;

                uint32_t immNext = noW;
                for(uint32_t p = pos; p < jLen; p++) {
                    if(posWindow[p] != noW) { immNext = posWindow[p]; break; }
                }
                if(immNext != nextW) continue;

                Shasta2AnchorId lastInPrev = Shasta2AnchorId(0);
                bool foundPrev = false;
                for(uint32_t p = 0; p < runStart; p++) {
                    if(posWindow[p] == prevW) {
                        lastInPrev = journey[p];
                        foundPrev = true;
                    }
                }

                Shasta2AnchorId firstInNext = Shasta2AnchorId(0);
                bool foundNext = false;
                for(uint32_t p = pos; p < jLen; p++) {
                    if(posWindow[p] == nextW) {
                        firstInNext = journey[p];
                        foundNext = true;
                        break;
                    }
                }

                if(foundPrev && foundNext) {
                    auto pairKey = pair<uint64_t, uint64_t>(
                        uint64_t(lastInPrev), uint64_t(firstInNext));
                    anchorPairCounts[pairKey]++;
                    anchorPairReads[pairKey].push_back(oid.getValue());
                    break;
                }
            }
        }

        BypassCandidate result{wIdx, prevW, nextW, Shasta2AnchorId(0), Shasta2AnchorId(0), {}};
        if(!anchorPairCounts.empty()) {
            auto bestIt = std::max_element(
                anchorPairCounts.begin(), anchorPairCounts.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });
            result.bestAnchorA = Shasta2AnchorId(bestIt->first.first);
            result.bestAnchorB = Shasta2AnchorId(bestIt->first.second);
            result.flowReadIds = anchorPairReads[bestIt->first];
        }
        return result;
    };

    for(uint32_t wIdx = 0; wIdx < windowCount; wIdx++) {
        const AnchorWindow& window = anchorWindows[wIdx];

        // Collect distinct predecessors and successors (excluding noWindow).
        set<uint32_t> predecessors, successors;
        for(const auto& [key, reads] : window.transitionReads) {
            if(key.first != noW) predecessors.insert(key.first);
            if(key.second != noW) successors.insert(key.second);
        }

        // Skip linear windows.
        if(predecessors.size() <= 1 && successors.size() <= 1) {
            ++linearSkipped;
            continue;
        }

        // Skip hairpin windows: any neighbor that normalizes to this
        // window means the window transitions to/from its own RC mirror.
        bool isHairpin = false;
        for(const uint32_t pred : predecessors) {
            if(normalize(pred) == wIdx) { isHairpin = true; break; }
        }
        if(!isHairpin) {
            for(const uint32_t succ : successors) {
                if(normalize(succ) == wIdx) { isHairpin = true; break; }
            }
        }
        if(isHairpin) {
            ++tripletFailed;
            continue;
        }

        // ---- Verkko-style triplet resolution. ----
        // A triplet (pred, X, succ) is "solid" if its read count >= minEdgeSupport.
        // The window is resolvable if every significant edge (coverage >= minEdgeCoverage)
        // is covered by at least one solid triplet.

        // Collect solid triplets and per-edge coverage.
        set<pair<uint32_t, uint32_t>> solidTriplets;
        set<uint32_t> tripletCoveredPreds, tripletCoveredSuccs;

        // Count per-edge coverage from ALL transitions (including edge-of-graph
        // entries where one side is noW). Matching Verkko, which counts
        // covered_in/out_neighbors from pairwise adjacency, not just full triplets.
        map<uint32_t, uint64_t> coveredInNeighbors, coveredOutNeighbors;
        for(const auto& [key, reads] : window.transitionReads) {
            if(key.first != noW) coveredInNeighbors[key.first] += reads.size();
            if(key.second != noW) coveredOutNeighbors[key.second] += reads.size();
            // Only full triplets (both sides known) can be solid.
            if(key.first != noW && key.second != noW && reads.size() >= minEdgeSupport) {
                solidTriplets.insert({key.first, key.second});
                tripletCoveredPreds.insert(key.first);
                tripletCoveredSuccs.insert(key.second);
            }
        }

        // Check: every significant predecessor edge must be covered.
        // An edge is skippable only if BOTH its coverage < minEdgeCoverage
        // AND the neighbor window is removable (low coverage).
        // Matching Verkko: weak edges to disposable neighbors are ignored.
        bool tripletClean = true;
        for(const uint32_t pred : predecessors) {
            auto it = coveredInNeighbors.find(pred);
            uint64_t cov = (it != coveredInNeighbors.end()) ? it->second : 0;
            if(cov < minEdgeCoverage && isRemovable(pred)) continue;
            if(tripletCoveredPreds.find(pred) == tripletCoveredPreds.end()) {
                const uint32_t normPred = normalize(pred);
                cout << "    Blocked: in-edge from window " << normPred
                     << " (cov=" << cov
                     << ", windowCov=" << (normPred < windowCount ? windowCoverage[normPred] : 0.)
                     << ") has no solid triplet" << endl;
                tripletClean = false;
                break;
            }
        }

        // Check: every significant successor edge must be covered.
        if(tripletClean) {
            for(const uint32_t succ : successors) {
                auto it = coveredOutNeighbors.find(succ);
                uint64_t cov = (it != coveredOutNeighbors.end()) ? it->second : 0;
                if(cov < minEdgeCoverage && isRemovable(succ)) continue;
                if(tripletCoveredSuccs.find(succ) == tripletCoveredSuccs.end()) {
                    const uint32_t normSucc = normalize(succ);
                    cout << "    Blocked: out-edge to window " << normSucc
                         << " (cov=" << cov
                         << ", windowCov=" << (normSucc < windowCount ? windowCoverage[normSucc] : 0.)
                         << ") has no solid triplet" << endl;
                    tripletClean = false;
                    break;
                }
            }
        }

        if(!tripletClean) {
            ++tripletFailed;
            continue;
        }

        // Log the tangle.
        cout << "  Tangle at window " << wIdx << ": "
             << predecessors.size() << " preds, "
             << successors.size() << " succs, "
             << solidTriplets.size() << " solid triplets" << endl;

        // Create bypass candidates for each solid triplet.
        for(const auto& [prevW, nextW] : solidTriplets) {
            auto trIt = window.transitionReads.find({prevW, nextW});
            if(trIt == window.transitionReads.end() || trIt->second.empty()) continue;

            auto cand = findAnchorPairs(wIdx, prevW, nextW, trIt->second);
            if(!cand.flowReadIds.empty()) {
                candidates.push_back(std::move(cand));
            }
        }
    }

    cout << "  Candidates: " << candidates.size()
         << " bypass pairs (" << linearSkipped << " linear, "
         << tripletFailed << " triplet failed)." << endl;

    // Sort candidates by flow read count (descending).
    std::sort(candidates.begin(), candidates.end(),
        [](const BypassCandidate& a, const BypassCandidate& b) {
            return a.flowReadIds.size() > b.flowReadIds.size();
        });

    // Process each candidate.
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

        if(commonReads.empty()) continue;

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

    totalProcessed += processedCount;
    cout << "  Resolved " << processedCount << " bypasses this step." << endl;

    } // end of resolveSteps iteration

    // ---- Step 6: Rebuild anchorMarkerInfos from the modified copy. ----
    anchors.anchorMarkerInfos.clear();
    for(const auto& anchorData : mutableAnchors) {
        anchors.anchorMarkerInfos.appendVector(anchorData);
    }

    cout << timestamp << "Detangle: " << detangledWindows.size()
         << " windows bypassed, "
         << totalProcessed << " bypasses processed, "
         << bypassEdges.size() << " bypass edges created ("
         << resolveSteps.size() << " resolve steps)." << endl;

    return detangledWindows.size();
}
