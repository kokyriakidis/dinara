// ============================================================================
// Detangle Case 1: 1-to-1 bypass
// ============================================================================
//
// Resolves tangled windows where multiple genomic paths (flows) pass
// through the same window. A window is tangled when it has ≥2 distinct
// flows, where a flow is a group of reads that enter from one neighbor
// window and exit to another.
//
// For each tangled window B with a flow triplet (A -> B -> C):
//   - The flow reads are the reads common between windows A and C
//     (reads that traverse B from A to C), taken from transitionReads.
//   - A bypass edge is created connecting A's exit anchor directly to
//     C's entry anchor, skipping over B.
//   - The flow reads are removed from every backbone anchor of B
//     (and their RC counterparts), so B's backbone retains only reads
//     that don't participate in any bypass flow.
//
// Processing order:
//   All bypass candidates across all tangled windows are collected,
//   sorted by flow read count (descending), and processed strongest
//   first. This ensures the most confident flows get their full read
//   sets. Weaker flows that share reads with stronger ones see the
//   reduced read sets after prior removals.
//
// Anchor modification:
//   Since anchorMarkerInfos is a MemoryMapped::VectorOfVectors that
//   cannot be modified in place, all anchor data is copied to a mutable
//   std::vector<std::vector<>> at the start. Read removals are applied
//   immediately so subsequent candidates see the updated state. The
//   VectorOfVectors is rebuilt once at the end.
//
// RC handling:
//   For each forward bypass edge (A -> C), an RC mirror edge (C' -> A')
//   is also created. When removing reads from backbone anchors, both
//   the canonical and RC anchors are updated (with strand-flipped read
//   IDs for the RC anchor).
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
    cout << timestamp << "Detangle Case 1 (1-to-1 bypass) begins." << endl;

    const uint32_t noW = AnchorWindow::noWindow;
    bypassEdges.clear();

    // ---- Step 1: Copy anchor data to a mutable structure. ----
    // This allows in-place read removal so each processed candidate's
    // changes are immediately visible to subsequent candidates.
    const uint64_t anchorCount = anchors.anchorMarkerInfos.size();
    vector<vector<Shasta2AnchorMarkerInfo>> mutableAnchors(anchorCount);
    for(uint64_t i = 0; i < anchorCount; i++) {
        const auto span = anchors.anchorMarkerInfos[i];
        mutableAnchors[i].assign(span.begin(), span.end());
    }

    // ---- Step 2: Collect bypass candidates from all tangled windows. ----
    // A bypass candidate represents a single directional flow
    // (prevW -> tangledW -> nextW). The flow reads are the reads
    // common between prevW and nextW, obtained from transitionReads.
    struct BypassCandidate {
        uint32_t windowId;                // The tangled window being bypassed.
        uint32_t prevW;                   // Previous window.
        uint32_t nextW;                   // Next window.
        Shasta2AnchorId prevExitAnchor;   // Last anchor in prev window.
        Shasta2AnchorId nextEntryAnchor;  // First anchor in next window.
        set<uint32_t> flowReads;          // Reads common between prevW and nextW.
    };
    vector<BypassCandidate> candidates;

    for(uint32_t wIdx = 0; wIdx < anchorWindows.size(); wIdx++) {
        const AnchorWindow& window = anchorWindows[wIdx];

        // Identify flows from transitionReads. Each (prev, next) pair
        // where both are real windows defines a directional flow.
        // Merge directional pairs (A->B) and (B->A) into a single
        // canonical flow to count distinct genomic paths.
        using FlowKey = pair<uint32_t, uint32_t>;
        map<FlowKey, set<uint32_t>> mergedFlows;

        for(const auto& [key, reads] : window.transitionReads) {
            if(key.first == noW || key.second == noW) continue;
            FlowKey canonical = {std::min(key.first, key.second),
                                 std::max(key.first, key.second)};
            for(const OrientedReadId& oid : reads) {
                mergedFlows[canonical].insert(oid.getValue());
            }
        }

        // A window is tangled only if it has ≥2 distinct flows.
        if(mergedFlows.size() < 2) continue;

        // Create a candidate for each directional flow.
        // Each directional key (prev -> next) gets its own bypass edge.
        for(const auto& [key, reads] : window.transitionReads) {
            if(key.first == noW || key.second == noW) continue;

            // Find the connection anchors from the inter-window edges.
            // inEdge from prevW: anchorIdA = last anchor in prevW.
            Shasta2AnchorId prevExitAnchor = Shasta2AnchorId(0);
            bool foundPrev = false;
            for(const auto& ie : window.inEdges) {
                if(ie.otherWindow == key.first) {
                    prevExitAnchor = ie.anchorIdA;
                    foundPrev = true;
                    break;
                }
            }

            // outEdge to nextW: anchorIdB = first anchor in nextW.
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
                FlowKey canonical = {std::min(key.first, key.second),
                                     std::max(key.first, key.second)};
                candidates.push_back({wIdx, key.first, key.second,
                                      prevExitAnchor, nextEntryAnchor,
                                      mergedFlows[canonical]});
            }
        }
    }

    // ---- Step 3: Sort candidates by flow read count (descending). ----
    // Strongest flows are processed first to get their full read sets.
    std::sort(candidates.begin(), candidates.end(),
        [](const BypassCandidate& a, const BypassCandidate& b) {
            return a.flowReads.size() > b.flowReads.size();
        });

    // ---- Step 4: Process each candidate. ----
    // For each candidate:
    //   a) Recompute common reads (prior removals may have reduced them).
    //   b) Create bypass edge (forward + RC).
    //   c) Remove common reads from backbone anchors of the tangled window.
    set<uint32_t> detangledWindows;
    uint64_t processedCount = 0;

    for(const auto& cand : candidates) {

        // (a) Recompute common reads between connection anchors.
        // Only keep reads that are still present on both anchors
        // AND were in the original flow read set.
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

        // (b) Create bypass edges.
        bypassEdges.push_back({cand.prevExitAnchor, cand.nextEntryAnchor});

        // RC mirror: B' -> A'.
        const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(cand.prevExitAnchor) ^ 1ULL);
        const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(cand.nextEntryAnchor) ^ 1ULL);
        bypassEdges.push_back({rcB, rcA});

        cout << "  Bypass: window " << cand.windowId
             << " " << cand.prevExitAnchor << " -> " << cand.nextEntryAnchor
             << " (+ RC " << rcB << " -> " << rcA << ")"
             << " prev=" << cand.prevW << " next=" << cand.nextW
             << " commonReads=" << commonReads.size() << endl;

        // (c) Remove common reads from backbone anchors of the tangled window.
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

            // Canonical anchor: remove reads matching flow read IDs.
            auto& canonicalData = mutableAnchors[uint64_t(canonicalId)];
            canonicalData.erase(
                std::remove_if(canonicalData.begin(), canonicalData.end(),
                    [&](const Shasta2AnchorMarkerInfo& mi) {
                        return commonReads.count(mi.orientedReadId.getValue());
                    }),
                canonicalData.end());

            // RC anchor: flow reads are forward-strand values,
            // RC anchor reads are strand-flipped.
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

    // ---- Step 5: Rebuild anchorMarkerInfos from the modified copy. ----
    anchors.anchorMarkerInfos.clear();
    for(const auto& anchorData : mutableAnchors) {
        anchors.anchorMarkerInfos.appendVector(anchorData);
    }

    cout << timestamp << "Detangle Case 1: " << detangledWindows.size()
         << " windows detangled, "
         << processedCount << " bypasses processed, "
         << bypassEdges.size() << " bypass edges created." << endl;

    return detangledWindows.size();
}
