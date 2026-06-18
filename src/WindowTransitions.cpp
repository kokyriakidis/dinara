// Compute per-window transition counts from journeys.
// Populates AnchorWindow::transitionReads, per-read previousWindow/nextWindow,
// and backbonePreviousWindow/backboneNextWindow.

#include "WindowTransitions.hpp"
#include "timestamp.hpp"

#include <iostream>
#include <map>
#include <vector>

using namespace dinara;
using std::cout;
using std::endl;
using std::map;
using std::pair;
using std::vector;



void dinara::computeWindowTransitions(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    vector<AnchorWindow>& anchorWindows,
    const vector<uint32_t>* anchorDovetailWindow)
{
    cout << timestamp << "computeWindowTransitions begins." << endl;

    const uint32_t windowCount = uint32_t(anchorWindows.size());
    const uint64_t anchorCount = anchors.anchorMarkerInfos.size();
    const uint32_t noW = AnchorWindowReadInterval::noWindow;

    // Build anchorToWindow map from backbone positions (same logic as
    // the anchor graph constructor). Forward windows get their windowId,
    // RC mirrors get windowId + windowCount.
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

    // Whole-journey dovetail membership (additive; matches the anchor graph
    // constructor). Backbone fill above is authoritative; dovetails only fill
    // anchors still unmapped. RC twin -> mirror window.
    if(anchorDovetailWindow != nullptr && !anchorDovetailWindow->empty() &&
       anchorDovetailWindow->size() == anchorCount) {
        for(uint64_t aid = 0; aid < anchorCount; aid++) {
            const uint32_t windowId = (*anchorDovetailWindow)[aid];
            if(windowId == noW) continue;
            if(anchorToWindow[aid] == noW) {
                anchorToWindow[aid] = windowId;
            }
            const uint64_t rcAid = aid ^ 1ULL;
            if(rcAid < anchorCount && anchorToWindow[rcAid] == noW) {
                anchorToWindow[rcAid] = windowId + windowCount;
            }
        }
    }

    auto normalize = [&](uint32_t w) -> uint32_t {
        return (w >= windowCount) ? (w - windowCount) : w;
    };

    // Clear any existing transition data.
    for(auto& window : anchorWindows) {
        window.transitionReads.clear();
        window.backbonePreviousWindow = noW;
        window.backboneNextWindow = noW;
        for(auto& ri : window.readIntervals) {
            ri.previousWindow = noW;
            ri.nextWindow = noW;
        }
    }

    // Build lookup: (windowId, orientedReadId) -> index in readIntervals.
    map<pair<uint32_t, uint32_t>, uint64_t> windowReadIndex;
    for(uint32_t wid = 0; wid < windowCount; wid++) {
        const auto& window = anchorWindows[wid];
        for(uint64_t ri = 0; ri < window.readIntervals.size(); ri++) {
            windowReadIndex[{wid, window.readIntervals[ri].orientedReadId.getValue()}] = ri;
        }
    }

    // Walk each read's journey and collect the normalized window sequence.
    const uint64_t journeyCount = journeys.size();
    for(uint64_t oidValue = 0; oidValue < journeyCount; oidValue++) {
        const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
        const auto journey = journeys[oid];
        if(journey.empty()) continue;

        // Collect the sequence of distinct normalized windows this read visits.
        vector<uint32_t> windowSequence;
        for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
            const Shasta2AnchorId anchorId = journey[pos];
            if(uint64_t(anchorId) >= anchorCount) continue;
            const uint32_t windowId = anchorToWindow[uint64_t(anchorId)];
            if(windowId == noW) continue;
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

            auto it = windowReadIndex.find({wid, uint32_t(oidValue)});
            if(it != windowReadIndex.end()) {
                auto& interval = anchorWindows[wid].readIntervals[it->second];
                interval.previousWindow = prev;
                interval.nextWindow = next;
            }

            // Always add to transitionReads (even if read isn't in readIntervals,
            // it still provides transition evidence).
            anchorWindows[wid].transitionReads[{prev, next}].push_back(oid);
        }
    }

    // Populate backbonePreviousWindow / backboneNextWindow for each window.
    for(uint32_t wid = 0; wid < windowCount; wid++) {
        auto& window = anchorWindows[wid];
        const uint32_t bbOid = window.backboneOrientedReadId.getValue();
        for(const auto& ri : window.readIntervals) {
            if(ri.orientedReadId.getValue() == bbOid) {
                window.backbonePreviousWindow = ri.previousWindow;
                window.backboneNextWindow = ri.nextWindow;
                break;
            }
        }
    }

    // Diagnostic summary.
    uint64_t totalTransitionEntries = 0;
    uint64_t totalTransitionReads = 0;
    uint64_t fullTripletEntries = 0;
    uint64_t fullTripletReads = 0;
    for(const auto& window : anchorWindows) {
        for(const auto& [key, reads] : window.transitionReads) {
            totalTransitionEntries++;
            totalTransitionReads += reads.size();
            if(key.first != noW && key.second != noW) {
                fullTripletEntries++;
                fullTripletReads += reads.size();
            }
        }
    }

    cout << timestamp << "computeWindowTransitions: "
         << windowCount << " windows, "
         << totalTransitionEntries << " transition entries ("
         << fullTripletEntries << " full triplets), "
         << totalTransitionReads << " total read-transitions ("
         << fullTripletReads << " in full triplets)." << endl;
}
