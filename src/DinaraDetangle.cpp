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

        // Filter by minimum coverage.
        vector<MergedFlow> flows;
        for(auto& mf : mergedFlows) {
            if(mf.readSet.size() >= minFlowCoverage) {
                flows.push_back(std::move(mf));
            }
        }

        if(flows.size() < 2) {
            continue;
        }

        const uint64_t pathCount = flows.size();

        cout << "  Window " << window.windowId
             << " backbone=" << window.backboneOrientedReadId
             << ": " << pathCount << " flows:";
        for(const auto& f : flows) {
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
        set<Shasta2AnchorId> alreadySplit;

        for(const uint32_t pos : positions) {
            const Shasta2AnchorId anchorId = backboneJourney[pos];
            const Shasta2AnchorId canonicalId = anchorId & ~Shasta2AnchorId(1);
            const Shasta2AnchorId rcId = canonicalId | Shasta2AnchorId(1);

            if(alreadySplit.count(canonicalId)) continue;
            alreadySplit.insert(canonicalId);

            // Copy anchor data before any appending.
            const auto canonicalSpan = anchors[canonicalId];
            const auto rcSpan = anchors[rcId];
            vector<Shasta2AnchorMarkerInfo> canonicalData(canonicalSpan.begin(), canonicalSpan.end());
            vector<Shasta2AnchorMarkerInfo> rcData(rcSpan.begin(), rcSpan.end());

            DeferredSplit canonicalSplit{canonicalId, {}};
            DeferredSplit rcSplit{rcId, {}};

            for(uint64_t pathIdx = 0; pathIdx < pathCount; pathIdx++) {
                const set<uint32_t>& flowReads = flows[pathIdx].readSet;

                // Canonical: keep reads that are in the flow's read set.
                vector<Shasta2AnchorMarkerInfo> canonicalSubset;
                for(const auto& mi : canonicalData) {
                    if(flowReads.count(mi.orientedReadId.getValue())) {
                        canonicalSubset.push_back(mi);
                    }
                }

                // RC: the flow reads are stored as forward-strand values,
                // but RC anchor reads are strand-flipped.
                vector<Shasta2AnchorMarkerInfo> rcSubset;
                for(const auto& mi : rcData) {
                    const uint32_t flippedValue = mi.orientedReadId.getValue() ^ 1;
                    if(flowReads.count(flippedValue)) {
                        rcSubset.push_back(mi);
                    }
                }

                // Only create a new anchor if this flow has reads here.
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
        }

        cout << "    Split " << alreadySplit.size() << " anchor pairs into "
             << pathCount << " paths each." << endl;
    }

    // Append all new anchors at once.
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
