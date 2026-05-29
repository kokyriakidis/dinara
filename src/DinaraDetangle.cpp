// Detangling: for each tangled window, identify flows (prev->current->next
// triplets), create bypass edges connecting prev directly to next, and
// remove flow reads from the current window's backbone anchors.
//
// Algorithm:
// 1. Find all triplets (W_prev -> X -> W_next) from transitionReads
// 2. Merge directional pairs (A->B) and (B->A) into single flows
// 3. Sort flows by common read count (descending)
// 4. For each flow: create a bypass edge from prev's exit anchor to
//    next's entry anchor, record the common reads
// 5. Remove all flow reads from every backbone anchor of window X

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
    uint64_t detangledCount = 0;
    bypassEdges.clear();

    // Track which anchors need flow reads removed, and which reads to remove.
    // Key: anchor ID. Value: set of OrientedReadId values to remove.
    map<Shasta2AnchorId, set<uint32_t>> readsToRemove;

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

        // Sort flows by read count (descending) — process largest first.
        std::sort(mergedFlows.begin(), mergedFlows.end(),
            [](const MergedFlow& a, const MergedFlow& b) {
                return a.readSet.size() > b.readSet.size();
            });

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

        ++detangledCount;

        // Collect all flow reads for removal from backbone.
        set<uint32_t> allFlowReads;

        for(const auto& flow : mergedFlows) {
            // For each directional key (prev->next), find the connection
            // anchors and create a bypass edge.
            for(const auto& fk : flow.flowKeys) {
                const uint32_t prevW = fk.first;
                const uint32_t nextW = fk.second;

                // Find the incoming edge from prevW: anchorIdA is the last
                // anchor in prevW (the exit point of prevW).
                Shasta2AnchorId prevExitAnchor = Shasta2AnchorId(0);
                bool foundPrev = false;
                for(const auto& ie : window.inEdges) {
                    if(ie.otherWindow == prevW) {
                        prevExitAnchor = ie.anchorIdA;
                        foundPrev = true;
                        break;
                    }
                }

                // Find the outgoing edge to nextW: anchorIdB is the first
                // anchor in nextW (the entry point of nextW).
                Shasta2AnchorId nextEntryAnchor = Shasta2AnchorId(0);
                bool foundNext = false;
                for(const auto& oe : window.outEdges) {
                    if(oe.otherWindow == nextW) {
                        nextEntryAnchor = oe.anchorIdB;
                        foundNext = true;
                        break;
                    }
                }

                if(foundPrev && foundNext) {
                    // Forward bypass edge.
                    bypassEdges.push_back({prevExitAnchor, nextEntryAnchor});

                    // RC bypass edge: B' -> A' (reversed direction).
                    const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(prevExitAnchor) ^ 1ULL);
                    const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(nextEntryAnchor) ^ 1ULL);
                    bypassEdges.push_back({rcB, rcA});

                    cout << "    Bypass: " << prevExitAnchor << " -> "
                         << nextEntryAnchor << " (+ RC " << rcB << " -> " << rcA << ")"
                         << " (prev=" << prevW << " next=" << nextW
                         << " reads=" << flow.readSet.size() << ")" << endl;
                }
            }

            // Accumulate flow reads.
            for(const uint32_t r : flow.readSet) {
                allFlowReads.insert(r);
            }
        }

        // Mark flow reads for removal from each backbone anchor.
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = journeys[backboneOid];
        const auto& positions = window.filteredBackbonePositions;

        set<Shasta2AnchorId> processedAnchors;
        for(const uint32_t pos : positions) {
            const Shasta2AnchorId anchorId = backboneJourney[pos];
            const Shasta2AnchorId canonicalId = anchorId & ~Shasta2AnchorId(1);
            const Shasta2AnchorId rcId = canonicalId | Shasta2AnchorId(1);

            if(processedAnchors.count(canonicalId)) continue;
            processedAnchors.insert(canonicalId);

            // For canonical anchor: remove reads whose OrientedReadId value
            // is in allFlowReads.
            for(const uint32_t r : allFlowReads) {
                readsToRemove[canonicalId].insert(r);
            }

            // For RC anchor: flow reads are stored as forward-strand values,
            // but RC anchor reads are strand-flipped.
            for(const uint32_t r : allFlowReads) {
                readsToRemove[rcId].insert(r ^ 1);
            }
        }

        cout << "    Marked " << allFlowReads.size() << " flow reads for removal from "
             << processedAnchors.size() << " backbone anchor pairs." << endl;
    }

    // Rebuild anchorMarkerInfos: copy all anchors, removing flow reads
    // from backbone anchors of detangled windows.
    if(!readsToRemove.empty()) {
        const uint64_t anchorCount = anchors.anchorMarkerInfos.size();
        vector<vector<Shasta2AnchorMarkerInfo>> allAnchors(anchorCount);

        for(uint64_t i = 0; i < anchorCount; i++) {
            const auto span = anchors.anchorMarkerInfos[i];
            auto it = readsToRemove.find(Shasta2AnchorId(i));
            if(it != readsToRemove.end()) {
                const auto& toRemove = it->second;
                for(const auto& mi : span) {
                    if(!toRemove.count(mi.orientedReadId.getValue())) {
                        allAnchors[i].push_back(mi);
                    }
                }
            } else {
                allAnchors[i].assign(span.begin(), span.end());
            }
        }

        // Rebuild the VectorOfVectors.
        anchors.anchorMarkerInfos.clear();
        for(const auto& anchorData : allAnchors) {
            anchors.anchorMarkerInfos.appendVector(anchorData);
        }
    }

    cout << timestamp << "Detangling: " << detangledCount
         << " windows detangled, "
         << bypassEdges.size() << " bypass edges created." << endl;

    return detangledCount;
}
