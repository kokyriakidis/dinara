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
// Anchor pair selection (step 3):
//   For each flow read, find the span of window X in the journey
//   (first and last positions belonging to X), then search for:
//     - lastInPrev:  last anchor in A *before* the X span
//     - firstInNext: first anchor in B *after* the X span
//   This X-span-relative search handles all topologies:
//
//   Case 1: A != B (distinct prev/next windows)
//     Journey: [a1, a2, x1, x2, x3, b1, b2]
//     X span: x1..x3.  lastInPrev=a2, firstInNext=b1.
//     Bypass edge: a2 -> b1.
//
//   Case 2: A == B (flow exits to the same window it entered from)
//     Journey: [a1, a2, x1, x2, x3, a3, a4]
//     X span: x1..x3.  lastInPrev=a2, firstInNext=a3.
//     Bypass edge: a2 -> a3.
//     Without the X-span approach, searching for "first nextW anchor"
//     would find a1 (since nextW == prevW), leaving no room for
//     lastInPrev before it.
//
//   Case 3: read revisits X (A -> X -> B -> X -> C)
//     Journey: [a1, x1, x2, b1, x3, x4, c1]
//     X runs: [x1,x2] and [x3,x4].
//     transitionReads records this read twice for window X:
//       {A,B} (first visit) and {B,C} (second visit).
//     For {A,B}: run [x1,x2] has A before and B after → match.
//       lastInPrev=a1, firstInNext=b1. Bypass: a1 → b1.
//     For {B,C}: run [x1,x2] has no B before → skip.
//       run [x3,x4] has B before and C after → match.
//       lastInPrev=b1, firstInNext=c1. Bypass: b1 → c1.
//
//   Case 4: other windows interleaved before X
//     Journey: [a1, c1, a2, x1, x2, b1]
//     X run: [x1,x2]. Immediate prev=A, immediate next=B.
//     lastInPrev=a2 (c1 skipped), firstInNext=b1.
//     Only anchors matching prevW/nextW are considered.
//
//   Case 5: same prevW, multiple X visits (A -> X -> A -> X -> B)
//     Journey: [a1, x1, a2, x2, x3, b1]
//     X runs: [x1] and [x2,x3].
//     For {A,B}: run [x1] has immediate prev=A, immediate next=A.
//       A != B → skip. Run [x2,x3] has immediate prev=A, immediate
//       next=B → match. lastInPrev=a2, firstInNext=b1.
//     The immediate-neighbor check prevents run [x1] from falsely
//     matching transition {A,B} just because B exists later.
//
//   Case 6: read starts in X (no prevW before X)
//     Journey: [x1, x2, b1, b2]
//     X run: [x1,x2]. Immediate prev=noW.
//     prevW != noW (filtered at entry) → no run matches.
//     Read is correctly skipped — it doesn't flow through X.
//
//   Case 7: read ends in X (no nextW after X)
//     Journey: [a1, a2, x1, x2]
//     X run: [x1,x2]. Immediate next=noW.
//     nextW != noW → no run matches. Read correctly skipped.
//
//   Case 8: unmapped anchors (noW) between prevW and X
//     Journey: [a1, ?1, ?2, x1, x2, b1]
//     posWindow: [A, noW, noW, X, X, B]
//     X run: [x1,x2]. Immediate prev scan skips noW positions,
//     finds A → match. lastInPrev=a1, firstInNext=b1.
//     Unmapped anchors don't break the immediate-neighbor check.
//
//   Case 9: noW gap inside an X visit
//     Journey: [a1, x1, ?1, x2, x3, b1]
//     posWindow: [A, X, noW, X, X, B]
//     The run extends through X and noW positions, stopping at B.
//     Single run: [x1, ?1, x2, x3] (pos 1..4). Immediate prev=A,
//     immediate next=B → match. lastInPrev=a1, firstInNext=b1.
//     noW gaps from unmapped anchors don't split the run.
//
// Processing order: candidates sorted by flow read count (descending)
// so stronger flows get their full read sets first. Weaker flows that
// share reads with stronger ones see reduced read sets after prior
// removals.
//
// RC handling: for each forward bypass edge (anchorA -> anchorB),
// an RC mirror (anchorB' -> anchorA') is also created. Read removals
// check both oid and oid^1 to handle flow reads on either strand.
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

                // Build per-position window labels for this journey.
                const uint32_t jLen = uint32_t(journey.size());
                vector<uint32_t> posWindow(jLen, noW);
                for(uint32_t pos = 0; pos < jLen; pos++) {
                    const Shasta2AnchorId anchorId = journey[pos];
                    if(uint64_t(anchorId) >= anchorCount) continue;
                    const uint32_t aw = anchorToWindow[uint64_t(anchorId)];
                    if(aw == noW) continue;
                    posWindow[pos] = normalize(aw);
                }

                // Scan for runs of window X (contiguous X positions,
                // possibly separated by noW gaps). For each run,
                // verify that the immediate neighboring windows match
                // (prevW, nextW) for this transition, then find the
                // best anchor pair across that boundary.
                //
                // "Run" = maximal stretch of positions that are either
                // X or noW, containing at least one X. A noW gap
                // inside an X visit (from unmapped anchors) doesn't
                // split the run. Only a position belonging to a
                // different window ends the run.
                //
                // "Immediate neighbor" = first mapped window (not noW)
                // before/after the run. This ensures each run matches
                // only its own transition.
                uint32_t pos = 0;
                while(pos < jLen) {
                    // Find start of an X run (first X position).
                    if(posWindow[pos] != wIdx) { pos++; continue; }
                    const uint32_t runStart = pos;
                    // Extend through X and noW positions. Stop at a
                    // position belonging to a different window.
                    // After this, pos points to the first position past
                    // the run (either a different window or jLen).
                    while(pos < jLen && (posWindow[pos] == wIdx || posWindow[pos] == noW)) pos++;

                    // Check immediate previous window: last mapped
                    // window before runStart must be prevW.
                    uint32_t immPrev = noW;
                    for(uint32_t p = runStart; p > 0; ) {
                        p--;
                        if(posWindow[p] != noW) { immPrev = posWindow[p]; break; }
                    }
                    if(immPrev != prevW) continue;

                    // Check immediate next window: first mapped
                    // window at or after pos must be nextW.
                    uint32_t immNext = noW;
                    for(uint32_t p = pos; p < jLen; p++) {
                        if(posWindow[p] != noW) { immNext = posWindow[p]; break; }
                    }
                    if(immNext != nextW) continue;

                    // lastInPrev: last anchor in prevW before runStart.
                    Shasta2AnchorId lastInPrev = Shasta2AnchorId(0);
                    bool foundPrev = false;
                    for(uint32_t p = 0; p < runStart; p++) {
                        if(posWindow[p] == prevW) {
                            lastInPrev = journey[p];
                            foundPrev = true;
                        }
                    }

                    // firstInNext: first anchor in nextW at or after pos.
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
                        break; // Use the first matching run for this read.
                    }
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
