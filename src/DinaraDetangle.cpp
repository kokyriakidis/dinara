// Detangling: split backbone anchors of tangled windows so that
// each path through the window gets its own anchor copies containing
// only that path's reads.
//
// A flow is defined by the reads common between a previous window and
// a next window. These reads traverse the current window's backbone.
// For each backbone anchor, a split copy is created containing the
// subset of flow reads present in that anchor.

#include "DinaraDetangle.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <iostream>
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
    map<Shasta2AnchorId, vector<Shasta2AnchorId>>& anchorSplitMap)
{
    cout << timestamp << "Detangling begins." << endl;

    const uint32_t noW = AnchorWindow::noWindow;
    uint64_t detangledCount = 0;
    anchorSplitMap.clear();

    // Collect all new anchors to append at the end, to avoid invalidating
    // spans into anchorMarkerInfos during iteration.
    vector<vector<Shasta2AnchorMarkerInfo>> newAnchors;

    // Deferred split: originalId -> index into newAnchors per path.
    // invalidIndex means the flow has no reads at this anchor.
    static constexpr uint64_t invalidIndex = uint64_t(-1);
    struct DeferredSplit {
        Shasta2AnchorId originalId;
        vector<uint64_t> newAnchorIndices;
    };
    vector<DeferredSplit> deferredSplits;

    for(const AnchorWindow& window : anchorWindows) {

        // Build flows from transitionReads: each (prev, next) pair where
        // both are real windows defines a flow. Merge directional pairs
        // (A->B) and (B->A) since they represent the same genomic path.
        using FlowKey = pair<uint32_t, uint32_t>;
        struct MergedFlow {
            vector<FlowKey> flowKeys;
            set<uint32_t> readSet; // OrientedReadId values.
        };
        map<FlowKey, uint64_t> canonicalToIdx;
        vector<MergedFlow> mergedFlows;

        for(const auto& [key, reads] : window.transitionReads) {
            if(key.first == noW || key.second == noW) continue;
            FlowKey canonical = {std::min(key.first, key.second),
                                 std::max(key.first, key.second)};
            auto it = canonicalToIdx.find(canonical);
            if(it == canonicalToIdx.end()) {
                canonicalToIdx[canonical] = mergedFlows.size();
                MergedFlow mf;
                mf.flowKeys.push_back(key);
                for(const OrientedReadId& oid : reads) {
                    mf.readSet.insert(oid.getValue());
                }
                mergedFlows.push_back(std::move(mf));
            } else {
                auto& mf = mergedFlows[it->second];
                mf.flowKeys.push_back(key);
                for(const OrientedReadId& oid : reads) {
                    mf.readSet.insert(oid.getValue());
                }
            }
        }

        if(mergedFlows.size() < 2) {
            continue;
        }

        const uint64_t pathCount = mergedFlows.size();

        cout << "  Window " << window.windowId
             << " backbone=" << window.backboneOrientedReadId
             << ": " << pathCount << " flows:";
        for(const auto& f : mergedFlows) {
            cout << " {";
            for(uint64_t k = 0; k < f.flowKeys.size(); k++) {
                if(k > 0) cout << "+";
                cout << f.flowKeys[k].first << "->" << f.flowKeys[k].second;
            }
            cout << "}=" << f.readSet.size();
        }
        cout << endl;

        // Get backbone info.
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = journeys[backboneOid];
        const auto& positions = window.filteredBackbonePositions;
        if(positions.empty()) continue;

        ++detangledCount;

        // Collect unique backbone anchor pairs (canonical, RC) in order.
        struct BackboneAnchorPair {
            Shasta2AnchorId canonicalId;
            Shasta2AnchorId rcId;
        };
        vector<BackboneAnchorPair> backbonePairs;
        set<Shasta2AnchorId> seen;
        for(const uint32_t pos : positions) {
            const Shasta2AnchorId anchorId = backboneJourney[pos];
            const Shasta2AnchorId canonicalId = anchorId & ~Shasta2AnchorId(1);
            if(seen.count(canonicalId)) continue;
            seen.insert(canonicalId);
            backbonePairs.push_back({canonicalId, canonicalId | Shasta2AnchorId(1)});
        }

        const uint64_t nPositions = backbonePairs.size();

        // Phase 1: For each flow and each backbone position, compute the
        // read subsets (canonical and RC). Store as sets of OrientedReadId
        // values for connectivity checking.
        // readSubsets[pathIdx][posIdx] = set of canonical-strand read values.
        vector<vector<set<uint32_t>>> readSubsets(pathCount, vector<set<uint32_t>>(nPositions));

        for(uint64_t pathIdx = 0; pathIdx < pathCount; pathIdx++) {
            const set<uint32_t>& flowReads = mergedFlows[pathIdx].readSet;
            for(uint64_t posIdx = 0; posIdx < nPositions; posIdx++) {
                const auto canonicalSpan = anchors[backbonePairs[posIdx].canonicalId];
                for(const auto& mi : canonicalSpan) {
                    if(flowReads.count(mi.orientedReadId.getValue())) {
                        readSubsets[pathIdx][posIdx].insert(mi.orientedReadId.getValue());
                    }
                }
            }
        }

        // Phase 2: For each flow, thin the chain to keep only positions
        // where consecutive pairs share reads. Use greedy bridging:
        // walk forward, skip positions that don't share reads with the
        // last kept position.
        // keep[pathIdx][posIdx] = true if this position is kept.
        vector<vector<bool>> keep(pathCount, vector<bool>(nPositions, false));

        for(uint64_t pathIdx = 0; pathIdx < pathCount; pathIdx++) {
            // Find positions where this flow has reads.
            vector<uint64_t> flowPositions;
            for(uint64_t posIdx = 0; posIdx < nPositions; posIdx++) {
                if(!readSubsets[pathIdx][posIdx].empty()) {
                    flowPositions.push_back(posIdx);
                }
            }
            if(flowPositions.size() < 2) continue;

            // Greedy bridging: keep positions that share reads with predecessor.
            vector<uint64_t> kept;
            kept.push_back(flowPositions[0]);
            for(uint64_t i = 1; i < flowPositions.size(); i++) {
                const uint64_t prevPos = kept.back();
                const uint64_t curPos = flowPositions[i];
                // Check if they share reads.
                const auto& prevReads = readSubsets[pathIdx][prevPos];
                const auto& curReads = readSubsets[pathIdx][curPos];
                bool sharesReads = false;
                for(const uint32_t r : prevReads) {
                    if(curReads.count(r)) {
                        sharesReads = true;
                        break;
                    }
                }
                if(sharesReads) {
                    kept.push_back(curPos);
                }
                // If not, skip curPos (it will not be split).
            }
            // Only keep if chain has at least 2 positions.
            if(kept.size() >= 2) {
                for(const uint64_t posIdx : kept) {
                    keep[pathIdx][posIdx] = true;
                }
            }
        }

        // Phase 3: Create split copies only for kept positions.
        uint64_t splitPairCount = 0;
        for(uint64_t posIdx = 0; posIdx < nPositions; posIdx++) {
            // Check if any flow keeps this position.
            bool anyKept = false;
            for(uint64_t pathIdx = 0; pathIdx < pathCount; pathIdx++) {
                if(keep[pathIdx][posIdx]) {
                    anyKept = true;
                    break;
                }
            }
            if(!anyKept) continue;

            const Shasta2AnchorId canonicalId = backbonePairs[posIdx].canonicalId;
            const Shasta2AnchorId rcId = backbonePairs[posIdx].rcId;

            // Copy anchor data.
            const auto canonicalSpan = anchors[canonicalId];
            const auto rcSpan = anchors[rcId];
            vector<Shasta2AnchorMarkerInfo> canonicalData(canonicalSpan.begin(), canonicalSpan.end());
            vector<Shasta2AnchorMarkerInfo> rcData(rcSpan.begin(), rcSpan.end());

            DeferredSplit canonicalSplit{canonicalId, {}};
            DeferredSplit rcSplit{rcId, {}};

            for(uint64_t pathIdx = 0; pathIdx < pathCount; pathIdx++) {
                if(!keep[pathIdx][posIdx]) {
                    // This flow doesn't keep this position — use original.
                    canonicalSplit.newAnchorIndices.push_back(invalidIndex);
                    rcSplit.newAnchorIndices.push_back(invalidIndex);
                    continue;
                }

                const set<uint32_t>& flowReads = mergedFlows[pathIdx].readSet;

                // Canonical subset.
                vector<Shasta2AnchorMarkerInfo> canonicalSubset;
                for(const auto& mi : canonicalData) {
                    if(flowReads.count(mi.orientedReadId.getValue())) {
                        canonicalSubset.push_back(mi);
                    }
                }

                // RC subset.
                vector<Shasta2AnchorMarkerInfo> rcSubset;
                for(const auto& mi : rcData) {
                    const uint32_t flippedValue = mi.orientedReadId.getValue() ^ 1;
                    if(flowReads.count(flippedValue)) {
                        rcSubset.push_back(mi);
                    }
                }

                if(!canonicalSubset.empty()) {
                    canonicalSplit.newAnchorIndices.push_back(newAnchors.size());
                    newAnchors.push_back(std::move(canonicalSubset));
                } else {
                    canonicalSplit.newAnchorIndices.push_back(invalidIndex);
                }

                if(!rcSubset.empty()) {
                    rcSplit.newAnchorIndices.push_back(newAnchors.size());
                    newAnchors.push_back(std::move(rcSubset));
                } else {
                    rcSplit.newAnchorIndices.push_back(invalidIndex);
                }
            }

            deferredSplits.push_back(std::move(canonicalSplit));
            deferredSplits.push_back(std::move(rcSplit));
            ++splitPairCount;
        }

        cout << "    Split " << splitPairCount << " anchor pairs into "
             << pathCount << " paths each." << endl;
    }

    // Append all new anchors at once. Original anchors are kept intact —
    // they serve as inter-window connection points with their full read sets.
    // Split copies contain flow-specific subsets for the parallel chains.
    const Shasta2AnchorId firstNewId = Shasta2AnchorId(anchors.anchorMarkerInfos.size());
    for(const auto& anchorData : newAnchors) {
        anchors.anchorMarkerInfos.appendVector(anchorData);
    }

    // Build the split map. invalidIndex entries get the original anchor ID.
    for(const auto& ds : deferredSplits) {
        vector<Shasta2AnchorId> splitIds;
        for(const uint64_t idx : ds.newAnchorIndices) {
            if(idx == invalidIndex) {
                splitIds.push_back(ds.originalId);
            } else {
                splitIds.push_back(Shasta2AnchorId(uint64_t(firstNewId) + idx));
            }
        }
        anchorSplitMap[ds.originalId] = splitIds;
    }

    cout << timestamp << "Detangling: " << detangledCount
         << " windows detangled, "
         << newAnchors.size() << " new anchors created, "
         << anchors.anchorMarkerInfos.size() << " total anchors." << endl;

    return detangledCount;
}
