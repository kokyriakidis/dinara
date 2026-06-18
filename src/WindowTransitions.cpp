// Compute per-window transition counts from journeys.
// Populates AnchorWindow::transitionReads, per-read previousWindow/nextWindow,
// and backbonePreviousWindow/backboneNextWindow.

#include "WindowTransitions.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

using namespace dinara;
using std::cout;
using std::endl;
using std::map;
using std::ofstream;
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
    // Also record each anchor's local index within its window's backbone order
    // (anchorLocalPos) and each window's backbone anchor count (windowAnchorLen)
    // so the branch analyzer can tell where along a window its neighbors attach.
    vector<uint32_t> anchorToWindow(anchorCount, noW);
    vector<uint32_t> anchorLocalPos(anchorCount, noW);
    vector<uint32_t> windowAnchorLen(windowCount, 0);
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

        windowAnchorLen[windowId] = uint32_t(positions.size());
        uint32_t localIdx = 0;
        for(const uint32_t pos : positions) {
            const uint64_t aid = uint64_t(backboneJourney[pos]);
            if(aid < anchorCount) {
                anchorToWindow[aid] = windowId;
                anchorLocalPos[aid] = localIdx;
            }
            const uint64_t rcAid = aid ^ 1ULL;
            if(rcAid < anchorCount) {
                anchorToWindow[rcAid] = windowId + windowCount;
                // RC mirror: local index is reversed (last backbone anchor
                // becomes first in the mirror frame).
                anchorLocalPos[rcAid] = localIdx;  // reversed below once len known
            }
            localIdx++;
        }
    }
    // Fix up RC-mirror local positions to reversed frame: pos' = len-1-pos.
    for(uint64_t aid = 0; aid < anchorCount; aid++) {
        const uint32_t raw = anchorToWindow[aid];
        if(raw == noW || raw < windowCount) continue;  // only mirrors
        const uint32_t wid = raw - windowCount;
        const uint32_t len = windowAnchorLen[wid];
        if(len > 0 && anchorLocalPos[aid] != noW) {
            anchorLocalPos[aid] = len - 1u - anchorLocalPos[aid];
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

    // ========================================================================
    // Window reach diagnostic (read-only).
    //
    // For each window W, in W's backbone orientation, find how far reads reach
    // to its left and right. Window IDs are creation order (not positional) and
    // strand reverses journeys, so reach is measured by BASE OFFSET along the
    // reads, in W's backbone frame, not by window-ID magnitude or step count:
    //
    //   - A read traverses W either forward (it visits W's forward anchors) or
    //     reversed (it visits W's RC mirror). The mirror flag tells us which.
    //   - Each window visit gets a base position from the read's anchor marker
    //     (anchors.getPosition). A neighbor at a LOWER base position than W (in
    //     W's backbone frame) is "left"; HIGHER is "right". For a read
    //     traversing W reversed, left/right are swapped (and base axis flips).
    //
    // Aggregated per window (option B: immediate + furthest):
    //   - immediate left/right: nearest neighbor window most reads agree on
    //     (consensus over reads, matching backbonePrev/Next semantics)
    //   - furthest left/right: the extreme neighbor any read reaches on that
    //     side, with the reach distance reported in BASES
    //
    // Output: WindowReach.csv. No graph or struct mutation.
    // ========================================================================
    {
        constexpr bool windowReachDiagnostic = true;
        if(windowReachDiagnostic) {
            // Immediate-neighbor votes: window -> (neighborWindow -> read count).
            vector<map<uint32_t, uint64_t>> leftImmediateVotes(windowCount);
            vector<map<uint32_t, uint64_t>> rightImmediateVotes(windowCount);
            // Furthest reach: window -> (neighborWindow, baseDistance).
            struct Reach { uint32_t window = noW; uint64_t bases = 0; };
            vector<Reach> leftFurthest(windowCount);
            vector<Reach> rightFurthest(windowCount);

            // Re-walk journeys, keeping orientation and a base position per visit.
            // basePos is the read-base position of the visit's first anchor.
            struct Visit { uint32_t normW; bool isMirror; uint32_t basePos; };
            for(uint64_t oidValue = 0; oidValue < journeyCount; oidValue++) {
                const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
                const auto journey = journeys[oid];
                if(journey.empty()) continue;

                // Distinct window visits with orientation and base position.
                vector<Visit> visits;
                for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
                    const Shasta2AnchorId anchorId = journey[pos];
                    if(uint64_t(anchorId) >= anchorCount) continue;
                    const uint32_t rawW = anchorToWindow[uint64_t(anchorId)];
                    if(rawW == noW) continue;
                    const uint32_t normW = normalize(rawW);
                    const bool isMirror = (rawW >= windowCount);
                    if(visits.empty() ||
                       visits.back().normW != normW ||
                       visits.back().isMirror != isMirror) {
                        const uint32_t basePos = anchors.getPosition(anchorId, oid);
                        visits.push_back({normW, isMirror, basePos});
                    }
                }

                // For each visited window, look at its neighbors in this read's
                // journey, mapped into the window's backbone frame, and measure
                // base-offset reach.
                for(uint64_t i = 0; i < visits.size(); i++) {
                    const uint32_t wid = visits[i].normW;
                    if(wid >= windowCount) continue;
                    const uint32_t wBase = visits[i].basePos;
                    const bool mirror = visits[i].isMirror;

                    // Journey-order neighbors mapped to W's backbone frame: a
                    // read traversing W via its RC mirror runs opposite to W's
                    // backbone, so journey-before is W's right and journey-after
                    // is W's left. Otherwise as-is.
                    const uint32_t beforeW = (i > 0) ? visits[i - 1].normW : noW;
                    const uint32_t afterW =
                        (i + 1 < visits.size()) ? visits[i + 1].normW : noW;
                    const uint32_t leftW = mirror ? afterW : beforeW;
                    const uint32_t rightW = mirror ? beforeW : afterW;

                    // Base distance helper (clamped to >= 0). Invalid positions
                    // (invalid<uint32_t>) yield 0.
                    auto baseDist = [&](uint32_t a, uint32_t b) -> uint64_t {
                        if(a == invalid<uint32_t> || b == invalid<uint32_t>) return 0;
                        return (a >= b) ? uint64_t(a - b) : uint64_t(b - a);
                    };

                    if(leftW != noW) {
                        leftImmediateVotes[wid][leftW]++;
                        // Furthest-left journey end (extreme on W's left side).
                        const uint64_t endIdx = mirror ? (visits.size() - 1) : 0;
                        const uint32_t endW = visits[endIdx].normW;
                        const uint64_t bases = baseDist(visits[endIdx].basePos, wBase);
                        if(endW != noW && endW < windowCount &&
                           bases >= leftFurthest[wid].bases) {
                            leftFurthest[wid] = {endW, bases};
                        }
                    }
                    if(rightW != noW) {
                        rightImmediateVotes[wid][rightW]++;
                        const uint64_t endIdx = mirror ? 0 : (visits.size() - 1);
                        const uint32_t endW = visits[endIdx].normW;
                        const uint64_t bases = baseDist(visits[endIdx].basePos, wBase);
                        if(endW != noW && endW < windowCount &&
                           bases >= rightFurthest[wid].bases) {
                            rightFurthest[wid] = {endW, bases};
                        }
                    }
                }
            }

            // Consensus immediate neighbor = most-voted, ties broken by smaller id.
            auto consensus = [&](const map<uint32_t, uint64_t>& votes,
                                 uint64_t& bestCount) -> uint32_t {
                uint32_t best = noW;
                bestCount = 0;
                for(const auto& [w, c] : votes) {
                    if(c > bestCount) { bestCount = c; best = w; }
                }
                return best;
            };

            ofstream csv("WindowReach.csv");
            csv << "WindowId,BackboneOrientedReadId,"
                   "ImmediateLeft,ImmediateLeftReads,"
                   "ImmediateRight,ImmediateRightReads,"
                   "FurthestLeft,FurthestLeftBases,"
                   "FurthestRight,FurthestRightBases\n";
            auto wstr = [&](uint32_t w) -> std::string {
                return (w == noW) ? std::string("none") : std::to_string(w);
            };
            for(uint32_t wid = 0; wid < windowCount; wid++) {
                uint64_t lCount = 0, rCount = 0;
                const uint32_t il = consensus(leftImmediateVotes[wid], lCount);
                const uint32_t ir = consensus(rightImmediateVotes[wid], rCount);
                csv << wid << ','
                    << anchorWindows[wid].backboneOrientedReadId.getValue() << ','
                    << wstr(il) << ',' << lCount << ','
                    << wstr(ir) << ',' << rCount << ','
                    << wstr(leftFurthest[wid].window) << ','
                    << leftFurthest[wid].bases << ','
                    << wstr(rightFurthest[wid].window) << ','
                    << rightFurthest[wid].bases << '\n';
            }
            csv.close();
            cout << timestamp << "Window reach diagnostic written to "
                    "WindowReach.csv (" << windowCount << " windows)." << endl;
        }
    }

    // ========================================================================
    // Window (prev, next) class diagnostic (read-only).
    //
    // Quantifies how many path copies each window would spawn under the
    // window-decomposition rewrite (docs/WindowDecompositionPlan.md). For each
    // window W, W.transitionReads groups reads by their (previousWindow,
    // nextWindow) route; every such class becomes one prev-next path under the
    // "support >= 1" rule. Partial classes (prev==noWindow or next==noWindow)
    // are the candidates for folding into full paths or becoming terminal tips.
    //
    // Two CSVs:
    //   - WindowClasses.csv: per-window summary (class counts + read totals,
    //     split into full / left-partial / right-partial / isolated).
    //   - WindowClassDetail.csv: one row per (window, prev, next) class with its
    //     read count and a kind label. Use to inspect individual routes.
    //
    // No graph or struct mutation.
    // ========================================================================
    {
        constexpr bool windowClassDiagnostic = true;
        if(windowClassDiagnostic) {
            auto wstr = [&](uint32_t w) -> std::string {
                return (w == noW) ? std::string("none") : std::to_string(w);
            };

            ofstream summary("WindowClasses.csv");
            summary << "WindowId,BackboneOrientedReadId,"
                       "TotalClasses,TotalReads,"
                       "FullClasses,FullReads,"
                       "LeftPartialClasses,LeftPartialReads,"
                       "RightPartialClasses,RightPartialReads,"
                       "IsolatedClasses,IsolatedReads,"
                       "MaxFullClassReads\n";

            ofstream detail("WindowClassDetail.csv");
            detail << "WindowId,Prev,Next,Kind,Reads\n";

            // Aggregate totals across all windows.
            uint64_t totFull = 0, totLeftPartial = 0;
            uint64_t totRightPartial = 0, totIsolated = 0;

            for(uint32_t wid = 0; wid < windowCount; wid++) {
                const auto& window = anchorWindows[wid];

                uint64_t fullClasses = 0, fullReads = 0;
                uint64_t leftPartialClasses = 0, leftPartialReads = 0;
                uint64_t rightPartialClasses = 0, rightPartialReads = 0;
                uint64_t isolatedClasses = 0, isolatedReads = 0;
                uint64_t maxFullClassReads = 0;
                uint64_t totalClasses = 0, totalReads = 0;

                for(const auto& [key, reads] : window.transitionReads) {
                    const uint32_t prev = key.first;
                    const uint32_t next = key.second;
                    const uint64_t n = reads.size();
                    totalClasses++;
                    totalReads += n;

                    const bool prevNone = (prev == noW);
                    const bool nextNone = (next == noW);
                    std::string kind;
                    if(!prevNone && !nextNone) {
                        kind = "full";
                        fullClasses++; fullReads += n;
                        maxFullClassReads = std::max(maxFullClassReads, n);
                    } else if(prevNone && !nextNone) {
                        kind = "leftPartial";
                        leftPartialClasses++; leftPartialReads += n;
                    } else if(!prevNone && nextNone) {
                        kind = "rightPartial";
                        rightPartialClasses++; rightPartialReads += n;
                    } else {
                        kind = "isolated";
                        isolatedClasses++; isolatedReads += n;
                    }

                    detail << wid << ',' << wstr(prev) << ',' << wstr(next)
                           << ',' << kind << ',' << n << '\n';
                }

                totFull += fullClasses;
                totLeftPartial += leftPartialClasses;
                totRightPartial += rightPartialClasses;
                totIsolated += isolatedClasses;

                summary << wid << ','
                        << window.backboneOrientedReadId.getValue() << ','
                        << totalClasses << ',' << totalReads << ','
                        << fullClasses << ',' << fullReads << ','
                        << leftPartialClasses << ',' << leftPartialReads << ','
                        << rightPartialClasses << ',' << rightPartialReads << ','
                        << isolatedClasses << ',' << isolatedReads << ','
                        << maxFullClassReads << '\n';
            }
            summary.close();
            detail.close();

            cout << timestamp << "Window class diagnostic written to "
                    "WindowClasses.csv / WindowClassDetail.csv ("
                 << windowCount << " windows; "
                 << totFull << " full, "
                 << totLeftPartial << " left-partial, "
                 << totRightPartial << " right-partial, "
                 << totIsolated << " isolated classes)." << endl;
        }
    }

    // ========================================================================
    // Doubled-vertex window arc diagnostic (read-only).
    //
    // Models the windows as a hifiasm-style doubled-vertex string graph, where
    // strand consistency is an invariant of the representation rather than a
    // cleanup step (see docs/WindowDecompositionPlan.md).
    //
    //   - Each window W becomes two oriented vertices:
    //         v = (W << 1) | strand,   strand 0 = forward, 1 = RC mirror.
    //         twin(v) = v ^ 1.
    //     This is exactly the forward window / RC-mirror pair already tracked
    //     in anchorToWindow (windowId vs windowId + windowCount).
    //   - A read that visits oriented window (W, sW) then (X, sX) emits a
    //     directed arc  vertex(W,sW) -> vertex(X,sX)  and its reverse-complement
    //     twin  vertex(X,sX)^1 -> vertex(W,sW)^1. Both are stored, so the graph
    //     is symmetric by construction and a +/- "strand contact" is not
    //     representable: orientation lives in the vertex id, not on the edge.
    //
    // Linearization = unitig extension over this doubled graph (the analogue of
    // asg_extend / stringGraphExtend): follow a vertex forward while it has a
    // unique outgoing arc whose target has a unique incoming arc. RC-chosen
    // windows are emitted from their strand-1 vertex (anchorId ^ 1 anchors).
    //
    // Outputs:
    //   - WindowArcs.csv: one row per directed arc (From, To with +/- strand,
    //     read support, and whether/how strongly its twin is present).
    //   - WindowUnitigs.csv: one row per maximal unitig (linear run), listing
    //     its oriented window path.
    // No graph or struct mutation.
    // ========================================================================
    {
        constexpr bool windowArcDiagnostic = true;
        if(windowArcDiagnostic) {
            // Minimum read support for an arc to be used (toggle).
            constexpr uint64_t minArcReads = 1;

            const uint32_t vertexCount = windowCount * 2u;
            auto vid = [&](uint32_t w, bool mirror) -> uint32_t {
                return (w << 1u) | (mirror ? 1u : 0u);
            };
            auto vWindow = [&](uint32_t v) -> uint32_t { return v >> 1u; };
            auto vStrandChar = [&](uint32_t v) -> char {
                return (v & 1u) ? '-' : '+';
            };

            // Accumulate directed arc weights (read counts) over journeys.
            // Key: (fromVertex, toVertex). Orientation is encoded in the
            // vertex ids, so each key is already strand-resolved.
            map<pair<uint32_t, uint32_t>, uint64_t> arcWeight;
            // Per-arc supporting read set (by oriented read value), used by the
            // branch analyzer to measure how cleanly competing branches at a
            // branch window separate by read membership.
            map<pair<uint32_t, uint32_t>, vector<uint32_t>> arcReads;
            // Per-arc exit attachment: the local index (within the from-window's
            // backbone order) of the last anchor the read occupied before
            // crossing to the next window. Normalized later by windowAnchorLen.
            // Spread of these values across a branch window's outgoing arcs
            // distinguishes a too-big window (exits spread along its length =>
            // positional split) from a true fan (exits clustered at the end =>
            // read-set split).
            map<pair<uint32_t, uint32_t>, vector<uint32_t>> arcExitLocal;

            for(uint64_t oidValue = 0; oidValue < journeyCount; oidValue++) {
                const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
                const auto journey = journeys[oid];
                if(journey.empty()) continue;

                // Distinct consecutive oriented window visits, recording the
                // local index of the last anchor in each visit (the exit point).
                vector<uint32_t> visits;       // oriented vertex ids
                vector<uint32_t> exitLocal;    // last local index per visit
                for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
                    const Shasta2AnchorId anchorId = journey[pos];
                    if(uint64_t(anchorId) >= anchorCount) continue;
                    const uint32_t rawW = anchorToWindow[uint64_t(anchorId)];
                    if(rawW == noW) continue;
                    const uint32_t normW = normalize(rawW);
                    const bool isMirror = (rawW >= windowCount);
                    const uint32_t v = vid(normW, isMirror);
                    const uint32_t lp = anchorLocalPos[uint64_t(anchorId)];
                    if(visits.empty() || visits.back() != v) {
                        visits.push_back(v);
                        exitLocal.push_back(lp);
                    } else {
                        exitLocal.back() = lp;  // advance to latest anchor in v
                    }
                }

                // Emit an arc per consecutive pair; the twin is added so the
                // graph is symmetric by construction (hifiasm arcId ^ 1 rule).
                for(uint64_t i = 0; i + 1 < visits.size(); i++) {
                    const uint32_t a = visits[i];
                    const uint32_t b = visits[i + 1];
                    if(a == b) continue;
                    arcWeight[{a, b}]++;
                    arcWeight[{b ^ 1u, a ^ 1u}]++;   // reverse-complement twin
                    arcReads[{a, b}].push_back(uint32_t(oidValue));
                    arcReads[{b ^ 1u, a ^ 1u}].push_back(uint32_t(oidValue));
                    if(exitLocal[i] != noW) {
                        arcExitLocal[{a, b}].push_back(exitLocal[i]);
                    }
                }
            }

            // Build filtered adjacency (degree-counted) for unitig extension.
            vector<vector<pair<uint32_t, uint64_t>>> outAdj(vertexCount);
            vector<uint32_t> outDeg(vertexCount, 0);
            vector<uint32_t> inDeg(vertexCount, 0);
            for(const auto& [key, w] : arcWeight) {
                if(w < minArcReads) continue;
                const uint32_t a = key.first;
                const uint32_t b = key.second;
                outAdj[a].push_back({b, w});
                outDeg[a]++;
                inDeg[b]++;
            }

            // Write the arc list (with twin support) to CSV.
            {
                ofstream arcsCsv("WindowArcs.csv");
                arcsCsv << "FromWindow,FromStrand,ToWindow,ToStrand,Reads,"
                           "TwinReads,TwinPresent\n";
                for(const auto& [key, w] : arcWeight) {
                    if(w < minArcReads) continue;
                    const uint32_t a = key.first;
                    const uint32_t b = key.second;
                    const auto twinIt = arcWeight.find({b ^ 1u, a ^ 1u});
                    const uint64_t twinW =
                        (twinIt != arcWeight.end()) ? twinIt->second : 0;
                    arcsCsv << vWindow(a) << ',' << vStrandChar(a) << ','
                            << vWindow(b) << ',' << vStrandChar(b) << ','
                            << w << ',' << twinW << ','
                            << ((twinW >= minArcReads) ? 1 : 0) << '\n';
                }
                arcsCsv.close();
            }

            // Unitig extension: a step a -> b is "mergeable" iff a has exactly
            // one outgoing arc and b has exactly one incoming arc. A unitig
            // starts at any vertex that is not such an unambiguous continuation
            // of a predecessor (in-degree != 1, or its sole predecessor
            // branches), as long as it has outgoing reach.
            auto uniqueOut = [&](uint32_t v, uint32_t& to) -> bool {
                if(outDeg[v] != 1) return false;
                to = outAdj[v].front().first;
                return true;
            };

            // Identify the unique predecessor of a vertex (only meaningful when
            // inDeg == 1). Built lazily from outAdj.
            vector<int64_t> solePred(vertexCount, -1);
            for(uint32_t a = 0; a < vertexCount; a++) {
                for(const auto& [b, w] : outAdj[a]) {
                    (void)w;
                    solePred[b] = (solePred[b] == -1) ? int64_t(a) : -2;
                }
            }
            auto isUnitigInterior = [&](uint32_t v) -> bool {
                // Continues a unitig iff in-degree 1 and that predecessor has
                // out-degree 1 (so the join is unambiguous on both sides).
                if(inDeg[v] != 1) return false;
                if(solePred[v] < 0) return false;
                return outDeg[uint32_t(solePred[v])] == 1;
            };

            vector<uint8_t> visited(vertexCount, 0);
            uint64_t branchVertices = 0, tipVertices = 0;
            uint64_t longestUnitig = 0;

            // Collect all unitig paths (as oriented vertex sequences) first, so
            // we can twin-dedup before emitting canonical linear paths.
            vector<vector<uint32_t>> unitigPaths;

            for(uint32_t start = 0; start < vertexCount; start++) {
                if(outDeg[start] > 1) branchVertices++;
                if(outDeg[start] == 0 && inDeg[start] == 0) continue;
                if(inDeg[start] == 0) tipVertices++;
                if(visited[start]) continue;
                if(isUnitigInterior(start)) continue;  // not a start

                // Walk forward.
                vector<uint32_t> path;
                uint32_t v = start;
                while(true) {
                    if(visited[v]) break;
                    visited[v] = 1;
                    path.push_back(v);
                    uint32_t to;
                    if(!uniqueOut(v, to)) break;
                    if(inDeg[to] != 1) break;       // target join is ambiguous
                    if(visited[to]) break;
                    v = to;
                }
                if(path.empty()) continue;

                longestUnitig = std::max(longestUnitig, uint64_t(path.size()));
                unitigPaths.push_back(std::move(path));
            }
            const uint64_t unitigCount = unitigPaths.size();

            // Emit raw unitigs (both strands) for inspection.
            {
                ofstream unitigsCsv("WindowUnitigs.csv");
                unitigsCsv << "UnitigId,WindowCount,WindowPath\n";
                for(uint64_t u = 0; u < unitigPaths.size(); u++) {
                    const auto& path = unitigPaths[u];
                    std::string pathStr;
                    for(uint64_t i = 0; i < path.size(); i++) {
                        if(i) pathStr += ' ';
                        pathStr += std::to_string(vWindow(path[i]));
                        pathStr += vStrandChar(path[i]);
                    }
                    unitigsCsv << (u + 1) << ',' << path.size() << ','
                               << pathStr << '\n';
                }
                unitigsCsv.close();
            }

            // ----------------------------------------------------------------
            // Canonical linearization: every contig appears twice in the
            // doubled graph (a path and its reverse-complement twin). The RC of
            // path [v0, v1, ..., vk] is [vk^1, ..., v1^1, v0^1]. Keep one
            // representative per twin pair; the survivor fixes each window's
            // orientation (the +/- it carries) and drops its RC duplicate. This
            // is the fw/rc separation: after dedup, each window appears in
            // exactly one orientation across the canonical paths, connected
            // linearly with no strand contacts.
            // ----------------------------------------------------------------
            {
                auto rcPath = [&](const vector<uint32_t>& p) {
                    vector<uint32_t> r(p.size());
                    for(uint64_t i = 0; i < p.size(); i++) {
                        r[p.size() - 1 - i] = p[i] ^ 1u;
                    }
                    return r;
                };
                // Canonical key: lexicographically smaller of (path, rcPath).
                auto canonicalKey = [&](const vector<uint32_t>& p) {
                    const vector<uint32_t> r = rcPath(p);
                    return (p <= r) ? p : r;
                };

                std::set<vector<uint32_t>> seen;
                vector<vector<uint32_t>> canonical;
                for(const auto& p : unitigPaths) {
                    auto key = canonicalKey(p);
                    if(seen.insert(key).second) {
                        canonical.push_back(key);
                    }
                }

                // Orientation assignment per window from the canonical paths.
                // -1 = unset, 0 = forward, 1 = RC, 2 = conflict (window appears
                // in both orientations across distinct canonical paths).
                vector<int8_t> windowOrient(windowCount, -1);
                uint64_t orientConflicts = 0;
                for(const auto& p : canonical) {
                    for(const uint32_t vtx : p) {
                        const uint32_t wid = vWindow(vtx);
                        const int8_t s = int8_t(vtx & 1u);
                        if(windowOrient[wid] == -1) windowOrient[wid] = s;
                        else if(windowOrient[wid] != s && windowOrient[wid] != 2) {
                            windowOrient[wid] = 2;
                            orientConflicts++;
                        }
                    }
                }

                ofstream linCsv("WindowLinearPaths.csv");
                linCsv << "PathId,WindowCount,OrientedWindowPath\n";
                uint64_t longestCanonical = 0;
                for(uint64_t u = 0; u < canonical.size(); u++) {
                    const auto& p = canonical[u];
                    longestCanonical = std::max(longestCanonical,
                                                uint64_t(p.size()));
                    std::string pathStr;
                    for(uint64_t i = 0; i < p.size(); i++) {
                        if(i) pathStr += ' ';
                        pathStr += std::to_string(vWindow(p[i]));
                        pathStr += vStrandChar(p[i]);
                    }
                    linCsv << (u + 1) << ',' << p.size() << ',' << pathStr
                           << '\n';
                }
                linCsv.close();

                // Per-window chosen orientation.
                {
                    ofstream orientCsv("WindowOrientation.csv");
                    orientCsv << "WindowId,Orientation\n";
                    for(uint32_t wid = 0; wid < windowCount; wid++) {
                        const char* o =
                            (windowOrient[wid] == 0) ? "+" :
                            (windowOrient[wid] == 1) ? "-" :
                            (windowOrient[wid] == 2) ? "conflict" : "unplaced";
                        orientCsv << wid << ',' << o << '\n';
                    }
                    orientCsv.close();
                }

                cout << timestamp << "Canonical linearization written to "
                        "WindowLinearPaths.csv / WindowOrientation.csv ("
                     << canonical.size() << " canonical paths (from "
                     << unitigCount << " unitigs), longest "
                     << longestCanonical << " windows; "
                     << orientConflicts << " orientation conflicts)." << endl;
            }

            // Per-window node report. For each window, record both strands'
            // degrees and classify it. The forward (+) and RC (-) vertices of a
            // window have mirror-symmetric degree by construction:
            //   outDeg(W+) == inDeg(W-)  and  inDeg(W+) == outDeg(W-).
            // So a window's actionable shape is fully described by its forward
            // vertex: leftDeg = inDeg(W+) (connections on the - side after
            // mirroring), rightDeg = outDeg(W+). Classes:
            //   linear : leftDeg <= 1 and rightDeg <= 1 (no ambiguity)
            //   branch : leftDeg > 1 or rightDeg > 1 (needs intra-window split)
            //   tip    : a side with degree 0 (assembly end)
            //   selfRC : window connects to its own mirror (inverted-repeat
            //            boundary); flagged because it breaks tip twin-pairing.
            uint64_t branchWindows = 0, tipWindows = 0, selfRcWindows = 0;
            uint64_t linearWindows = 0;
            {
                // Detect self-RC adjacency: an arc whose endpoints are the same
                // window in opposite strands.
                vector<uint8_t> selfRc(windowCount, 0);
                for(const auto& [key, w] : arcWeight) {
                    if(w < minArcReads) continue;
                    if(vWindow(key.first) == vWindow(key.second) &&
                       (key.first ^ 1u) == key.second) {
                        selfRc[vWindow(key.first)] = 1;
                    }
                }

                ofstream nodesCsv("WindowNodes.csv");
                nodesCsv << "WindowId,BackboneOrientedReadId,"
                            "LeftDeg,RightDeg,SelfRC,Class\n";
                for(uint32_t wid = 0; wid < windowCount; wid++) {
                    const uint32_t fwd = vid(wid, false);
                    const uint32_t leftDeg = inDeg[fwd];
                    const uint32_t rightDeg = outDeg[fwd];
                    const bool isSelfRc = selfRc[wid] != 0;

                    std::string cls;
                    if(isSelfRc) {
                        cls = "selfRC";
                        selfRcWindows++;
                    } else if(leftDeg > 1 || rightDeg > 1) {
                        cls = "branch";
                        branchWindows++;
                    } else if(leftDeg == 0 || rightDeg == 0) {
                        cls = "tip";
                        tipWindows++;
                    } else {
                        cls = "linear";
                        linearWindows++;
                    }

                    nodesCsv << wid << ','
                             << anchorWindows[wid].backboneOrientedReadId.getValue()
                             << ',' << leftDeg << ',' << rightDeg << ','
                             << (isSelfRc ? 1 : 0) << ',' << cls << '\n';
                }
                nodesCsv.close();
            }

            // Branch analyzer. For each branch side (a vertex with >1 outgoing
            // arcs), measure how cleanly the competing branches separate by
            // supporting read set. In the current pipeline per-window
            // readClusters / alternatePaths are not populated when this runs, so
            // read-set separation is the available signal for whether a branch
            // is a resolvable haplotype/context split (disjoint reads) or a
            // genuine shared locus / repeat (heavily shared reads).
            //
            // For a vertex v with outgoing targets {b1, b2, ...}, we report the
            // maximum pairwise Jaccard overlap of their read sets. Low max
            // Jaccard => branches are read-disjoint => splittable. High =>
            // branches share reads => not separable by read partition alone.
            {
                ofstream branchCsv("WindowBranches.csv");
                branchCsv << "Window,Strand,Side,Degree,Targets,"
                             "MaxPairJaccard,MinPairShared,"
                             "ExitMinFrac,ExitMaxFrac,ExitSpreadFrac,"
                             "Verdict,SplitKind\n";

                auto readsForArc =
                    [&](uint32_t a, uint32_t b) -> const vector<uint32_t>& {
                        static const vector<uint32_t> empty;
                        auto it = arcReads.find({a, b});
                        return (it != arcReads.end()) ? it->second : empty;
                    };
                // Jaccard over two sorted-unique-able read lists.
                auto jaccard = [&](vector<uint32_t> x, vector<uint32_t> y,
                                   uint64_t& shared) -> double {
                    std::sort(x.begin(), x.end());
                    x.erase(std::unique(x.begin(), x.end()), x.end());
                    std::sort(y.begin(), y.end());
                    y.erase(std::unique(y.begin(), y.end()), y.end());
                    uint64_t inter = 0;
                    uint64_t i = 0, j = 0;
                    while(i < x.size() && j < y.size()) {
                        if(x[i] == y[j]) { inter++; i++; j++; }
                        else if(x[i] < y[j]) i++;
                        else j++;
                    }
                    shared = inter;
                    const uint64_t uni = x.size() + y.size() - inter;
                    return uni ? double(inter) / double(uni) : 0.0;
                };

                // Median exit fraction of an arc (where along the from-window
                // its reads leave, 0=window start .. 1=window end).
                auto arcExitMedianFrac =
                    [&](uint32_t a, uint32_t b) -> double {
                        auto it = arcExitLocal.find({a, b});
                        if(it == arcExitLocal.end() || it->second.empty())
                            return -1.0;
                        const uint32_t len = windowAnchorLen[vWindow(a)];
                        if(len <= 1) return 1.0;
                        vector<uint32_t> v = it->second;
                        std::sort(v.begin(), v.end());
                        const uint32_t med = v[v.size() / 2];
                        return double(med) / double(len - 1u);
                    };

                uint64_t splittable = 0, sharedLocus = 0;
                uint64_t positionalSplits = 0, readSetSplits = 0;
                // Iterate the "outgoing" side of every vertex (forward branch).
                for(uint32_t v = 0; v < vertexCount; v++) {
                    if(outDeg[v] <= 1) continue;
                    const auto& targets = outAdj[v];

                    double maxJ = 0.0;
                    uint64_t minShared = std::numeric_limits<uint64_t>::max();
                    for(uint64_t i = 0; i < targets.size(); i++) {
                        for(uint64_t j = i + 1; j < targets.size(); j++) {
                            uint64_t shared = 0;
                            const double jc = jaccard(
                                readsForArc(v, targets[i].first),
                                readsForArc(v, targets[j].first),
                                shared);
                            maxJ = std::max(maxJ, jc);
                            minShared = std::min(minShared, shared);
                        }
                    }
                    if(minShared == std::numeric_limits<uint64_t>::max())
                        minShared = 0;

                    // Exit-attachment spread across the branch's arcs.
                    double exitMin = 2.0, exitMax = -1.0;
                    for(const auto& [tgt, w] : targets) {
                        (void)w;
                        const double f = arcExitMedianFrac(v, tgt);
                        if(f < 0.0) continue;
                        exitMin = std::min(exitMin, f);
                        exitMax = std::max(exitMax, f);
                    }
                    if(exitMax < 0.0) { exitMin = 0.0; exitMax = 0.0; }
                    const double exitSpread = exitMax - exitMin;

                    // Verdict: read-disjoint branches (no shared reads, low
                    // Jaccard) are splittable; otherwise a shared locus.
                    const bool isSplittable = (minShared == 0 && maxJ < 0.1);
                    if(isSplittable) splittable++; else sharedLocus++;

                    // Split kind: exits spread along the window length => the
                    // window absorbed sub-loci (positional split). Exits
                    // clustered (small spread) => a true fan (read-set split).
                    const char* splitKind =
                        (exitSpread > 0.25) ? "positional" : "readSet";
                    if(isSplittable) {
                        if(exitSpread > 0.25) positionalSplits++;
                        else                  readSetSplits++;
                    }

                    std::string tgtStr;
                    for(uint64_t i = 0; i < targets.size(); i++) {
                        if(i) tgtStr += ' ';
                        tgtStr += std::to_string(vWindow(targets[i].first));
                        tgtStr += vStrandChar(targets[i].first);
                    }

                    branchCsv << vWindow(v) << ',' << vStrandChar(v) << ','
                              << "out," << outDeg[v] << ',' << tgtStr << ','
                              << maxJ << ',' << minShared << ','
                              << exitMin << ',' << exitMax << ',' << exitSpread
                              << ',' << (isSplittable ? "splittable" : "sharedLocus")
                              << ',' << splitKind << '\n';
                }
                branchCsv.close();

                cout << timestamp << "Window branches written to "
                        "WindowBranches.csv ("
                     << (splittable + sharedLocus) << " branch sides: "
                     << splittable << " read-disjoint (splittable), "
                     << sharedLocus << " read-sharing (shared locus); "
                     << positionalSplits << " positional, "
                     << readSetSplits << " read-set split kind)." << endl;
            }

            // Strand-contact / symmetry sanity: every used arc must have its
            // twin present. By construction it does; report any violation.
            uint64_t arcsUsed = 0, twinMissing = 0;
            for(const auto& [key, w] : arcWeight) {
                if(w < minArcReads) continue;
                arcsUsed++;
                const auto twinIt =
                    arcWeight.find({key.second ^ 1u, key.first ^ 1u});
                if(twinIt == arcWeight.end() || twinIt->second < minArcReads) {
                    twinMissing++;
                }
            }

            cout << timestamp << "Window arc diagnostic written to "
                    "WindowArcs.csv / WindowUnitigs.csv ("
                 << vertexCount << " oriented vertices, "
                 << arcsUsed << " arcs, "
                 << unitigCount << " unitigs, longest "
                 << longestUnitig << " windows; "
                 << branchVertices << " branch vertices, "
                 << tipVertices << " tips; "
                 << twinMissing << " twin-missing (strand-contact) arcs)."
                 << endl;
            cout << timestamp << "Window nodes written to WindowNodes.csv ("
                 << linearWindows << " linear, "
                 << branchWindows << " branch, "
                 << tipWindows << " tip, "
                 << selfRcWindows << " self-RC windows)." << endl;
        }
    }
}
