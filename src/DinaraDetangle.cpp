// Detangling: split backbone anchors of tangled windows so that
// each path through the window gets its own anchor copies containing
// only that path's reads.

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

    // Deferred split map entries: originalId -> indices into newAnchors per path.
    struct DeferredSplit {
        Shasta2AnchorId originalId;
        vector<uint64_t> newAnchorIndices;
    };
    vector<DeferredSplit> deferredSplits;

    for(const AnchorWindow& window : anchorWindows) {

        // Collect through-flows: (prev, next) pairs where both are real windows,
        // with sufficient read support.
        using Flow = pair<pair<uint32_t, uint32_t>, vector<OrientedReadId>>;
        vector<Flow> throughFlows;
        for(const auto& [key, reads] : window.transitionReads) {
            if(key.first != noW && key.second != noW &&
               reads.size() >= minFlowCoverage) {
                throughFlows.push_back({key, reads});
            }
        }

        if(throughFlows.size() < 2) {
            continue;
        }

        ++detangledCount;
        const uint64_t pathCount = throughFlows.size();

        cout << "  Window " << window.windowId
             << " backbone=" << window.backboneOrientedReadId
             << ": " << pathCount << " through-flows:";
        for(const auto& [key, reads] : throughFlows) {
            cout << " (" << key.first << "->" << key.second
                 << ")=" << reads.size();
        }
        cout << endl;

        // Build per-path read sets.
        vector<set<uint32_t>> pathReadSets(pathCount);
        for(uint64_t pathIdx = 0; pathIdx < pathCount; pathIdx++) {
            for(const OrientedReadId& oid : throughFlows[pathIdx].second) {
                pathReadSets[pathIdx].insert(oid.getValue());
            }
        }

        // Get backbone anchor IDs from the backbone read's journey.
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = journeys[backboneOid];
        const auto& positions = window.filteredBackbonePositions;
        if(positions.empty()) continue;

        set<Shasta2AnchorId> alreadySplit;

        for(const uint32_t pos : positions) {
            const Shasta2AnchorId anchorId = backboneJourney[pos];
            const Shasta2AnchorId canonicalId = anchorId & ~Shasta2AnchorId(1);
            const Shasta2AnchorId rcId = canonicalId | Shasta2AnchorId(1);

            if(alreadySplit.count(canonicalId)) continue;
            alreadySplit.insert(canonicalId);

            // Copy anchor data to local vectors before any appending.
            const auto canonicalSpan = anchors[canonicalId];
            const auto rcSpan = anchors[rcId];
            vector<Shasta2AnchorMarkerInfo> canonicalData(canonicalSpan.begin(), canonicalSpan.end());
            vector<Shasta2AnchorMarkerInfo> rcData(rcSpan.begin(), rcSpan.end());

            DeferredSplit canonicalSplit{canonicalId, {}};
            DeferredSplit rcSplit{rcId, {}};

            for(uint64_t pathIdx = 0; pathIdx < pathCount; pathIdx++) {
                const set<uint32_t>& pathReads = pathReadSets[pathIdx];

                // Canonical subset: reads in this path.
                vector<Shasta2AnchorMarkerInfo> canonicalSubset;
                for(const auto& mi : canonicalData) {
                    if(pathReads.count(mi.orientedReadId.getValue())) {
                        canonicalSubset.push_back(mi);
                    }
                }

                // RC subset: RC anchor has strand-flipped reads.
                vector<Shasta2AnchorMarkerInfo> rcSubset;
                for(const auto& mi : rcData) {
                    const uint32_t flippedValue = mi.orientedReadId.getValue() ^ 1;
                    if(pathReads.count(flippedValue)) {
                        rcSubset.push_back(mi);
                    }
                }

                canonicalSplit.newAnchorIndices.push_back(newAnchors.size());
                newAnchors.push_back(std::move(canonicalSubset));

                rcSplit.newAnchorIndices.push_back(newAnchors.size());
                newAnchors.push_back(std::move(rcSubset));
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

    // Build the split map with actual anchor IDs.
    for(const auto& ds : deferredSplits) {
        vector<Shasta2AnchorId> splitIds;
        for(const uint64_t idx : ds.newAnchorIndices) {
            splitIds.push_back(Shasta2AnchorId(uint64_t(firstNewId) + idx));
        }
        anchorSplitMap[ds.originalId] = splitIds;
    }

    cout << timestamp << "Detangling: " << detangledCount
         << " windows detangled, "
         << newAnchors.size() << " new anchors created, "
         << anchors.anchorMarkerInfos.size() << " total anchors." << endl;

    return detangledCount;
}
