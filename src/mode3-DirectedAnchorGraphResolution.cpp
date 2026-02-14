// ============================================================================
// DirectedAnchorGraph — Resolution engine.
//
// Verkko-style two-node representation: each segment has oriented IDs
// fwd = segId*2, rev = segId*2+1.  rcNode(n) = n ^ 1.
//
// Implements:
//   - getValidTriplets()  — count and validate triplets for a segment
//   - resolveNodes()      — create edge nodes, rewire graph
//   - resolveHairpins()   — handle palindromic (self-RC) segments
//   - resolveRound()      — single round processing segments by length
//   - unitigifyOne()      — merge linear chains
//   - unitigifyAll()      — unitigify all segments
//   - extendForward()     — walk forward while in-degree=1, out-degree=1
//   - replaceUnitig()     — merge a chain into one segment
//   - replacePathNodes()  — rewrite paths after resolution
//   - splitPathsAtBreaks()— fix paths with broken edges
//   - runResolution()     — full pipeline
// ============================================================================

#include "mode3-DirectedAnchorGraph.hpp"
#include "DINARA_ASSERT.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <limits>
#include <optional>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

using namespace dinara;
using namespace mode3;
using namespace std;


// ============================================================================
// getValidTriplets — equivalent to Verkko's get_valid_triplets().
// ============================================================================

uint64_t DirectedAnchorGraph::pathGroupWeight(uint64_t pathIdx) const
{
    if(pathIdx >= paths.size() || pathIdx >= pathRemoved.size() || pathRemoved[pathIdx]) {
        return 0;
    }
    if(pathIdx < pathReadIds.size() && !pathReadIds[pathIdx].empty()) {
        return pathReadIds[pathIdx].size();
    }
    if(pathIdx < pathWeights.size()) {
        return pathWeights[pathIdx];
    }
    return 0;
}

uint64_t DirectedAnchorGraph::getCrossingCount(uint64_t segId) const
{
    uint64_t total = 0;
    const auto& crossing = pathsCrossing.getPathsCrossing(segId);
    for(const auto& occ : crossing) {
        // If a path crosses twice, we count it twice. Matches MBG behavior.
        total += pathGroupWeight(occ.pathIdx);
    }
    return total;
}

uint64_t DirectedAnchorGraph::getTrimAmountToCheck(
    DagNodeId from,
    DagNodeId to) const
{
    if(!nodeExists(segmentOf(from)) || !nodeExists(segmentOf(to))) {
        return 0;
    }
    const uint64_t overlap = getBpOverlap(from, to);
    const uint64_t toLen = nodeBpLength(segmentOf(to));
    if(overlap >= toLen) {
        return 0;
    }

    uint64_t minClip = nodeBpLength(segmentOf(from));
    bool foundBoundarySupport = false;
    const auto& crossing = pathsCrossing.getPathsCrossing(segmentOf(from));
    for(const auto& occ : crossing) {
        const uint64_t pathIdx = occ.pathIdx;
        const size_t i = occ.offset;

        if(pathIdx >= paths.size() || pathRemoved[pathIdx] ||
           pathIdx >= pathReadIds.size()) {
            continue;
        }
        const auto& path = paths[pathIdx];
        if(pathReadIds[pathIdx].empty() || path.size() < 2) continue;

        // Verify that this occurrence corresponds to the edge (from -> to).
        // path[i] is segmentOf(from).
        // Check if forward or reverse match.
        bool match = false;
        
        // Forward: path[i] == from
        if(path[i] == from) {
            if(i + 1 < path.size() && path[i+1] == to) {
                // Determine if this is a boundary read.
                if(isForward(from)) {
                    // from is forward, so we look at rightClip?
                    // Original code: if(isForward(from)) { if(i+1 != path.size()-1) continue; ... }
                    // Wait, original code:
                    /*
                    for(size_t i = 0; i + 1 < path.size(); i++) {
                        if(path[i] != from || path[i + 1] != to) continue;
                        if(isForward(from)) {
                            if(i + 1 != path.size() - 1) continue;
                    */
                   // So we only care if this edge is the LAST edge of the path (for forward 'from')?
                   if(i + 1 == path.size() - 1) {
                       match = true;
                       for(const auto& read : pathReadIds[pathIdx]) {
                           minClip = min(minClip, read.rightClip);
                           foundBoundarySupport = true;
                       }
                   }
                } else {
                    // from is reverse. Original code:
                    /*
                    } else {
                        if(i != 0) continue;
                    */
                   if(i == 0) {
                       match = true;
                       for(const auto& read : pathReadIds[pathIdx]) {
                           minClip = min(minClip, read.leftClip); // Wait, original used leftClip for !isForward(from)?
                           // If from is reverse, it's <A. And edge is <A -> B.
                           // This means path starts with <A.
                           // So it's a start tip?
                           foundBoundarySupport = true;
                       }
                   }
                }
            }
        }
        
        // Reverse: path[i] == rc(from)?
        // Original loop checked `path[i] != from`.
        // So it implicitly required `path[i] == from`.
        // It did NOT handle `rc(from)` case in original code?
        // Let's re-read original code.
        /*
        for(size_t i = 0; i + 1 < path.size(); i++) {
            if(path[i] != from || path[i + 1] != to) continue;
            if(isForward(from)) { ... } else { ... }
        }
        */
        // Yes, it only checked `path[i] == from`.
        // It did not check `rc(from)`.
        // Does this mean `getTrimAmountToCheck` relies on `from` being the oriented node?
        // Yes.
        // And path must contain that oriented node.
        // So if path contains `rc(from)`, it is ignored.
        //
        // So here, `occ` gives us `path[i]`. `path[i]` must be `from`.
        // But `occ` is for `segmentOf(from)`. So `path[i]` could be `from` or `rc(from)`.
        
        if (!match && path[i] == from && i + 1 < path.size() && path[i+1] == to) {
             if(isForward(from)) {
                if(i + 1 == path.size() - 1) {
                    for(const auto& read : pathReadIds[pathIdx]) {
                        minClip = min(minClip, read.rightClip);
                        foundBoundarySupport = true;
                    }
                }
            } else {
                if(i == 0) {
                    for(const auto& read : pathReadIds[pathIdx]) {
                        minClip = min(minClip, read.leftClip);
                        foundBoundarySupport = true;
                    }
                }
            }
        }
        // Note: original code `if(isForward(from))` checked `i+1 != path.size()-1`.
        // `else` checked `i != 0`.
    }
    if(!foundBoundarySupport) {
        return 0;
    }
    if(minClip == 1) {
        return toLen - overlap;
    }
    if(minClip == 0) {
        if(toLen > overlap + 1) {
            return toLen - overlap - 1;
        }
        return toLen - overlap;
    }
    return 0;
}

uint64_t DirectedAnchorGraph::getAnchorSize(
    DagNodeId from,
    DagNodeId to) const
{
    if(!nodeExists(segmentOf(from)) || !nodeExists(segmentOf(to))) {
        return 0;
    }
    const uint64_t toLen = nodeBpLength(segmentOf(to));
    uint64_t result = 0;
    const auto& crossing = pathsCrossing.getPathsCrossing(segmentOf(from));
    for(const auto& occ : crossing) {
        const uint64_t pathIdx = occ.pathIdx;
        const size_t i = occ.offset;

        if(pathIdx >= paths.size() || pathRemoved[pathIdx] ||
           pathIdx >= pathReadIds.size()) {
            continue;
        }
        const auto& path = paths[pathIdx];
        if(path.size() < 2 || pathReadIds[pathIdx].empty()) continue;

        // Check forward: from -> to.
        // If occ refers to 'from' node.
        if(path[i] == from) {
            if(i + 1 < path.size() && path[i + 1] == to) {
                if(i + 1 < path.size() - 1) {
                    return toLen;
                }
                for(const auto& read : pathReadIds[pathIdx]) {
                    if(toLen > read.rightClip) {
                        result = max(result, toLen - read.rightClip);
                    }
                }
            }
        }

        // Check reverse: rc(to) -> rc(from).
        // This corresponds to path[i-1] == rc(to) && path[i] == rc(from).
        if(path[i] == rcNode(from)) {
            if(i > 0 && path[i - 1] == rcNode(to)) {
                if(i > 1) {
                    return toLen;
                }
                for(const auto& read : pathReadIds[pathIdx]) {
                    if(toLen > read.leftClip) {
                        result = max(result, toLen - read.leftClip);
                    }
                }
            }
        }
    }
    return result;
}

uint64_t DirectedAnchorGraph::createFakeEdgeNode(
    DagNodeId from,
    DagNodeId to,
    const unordered_set<uint64_t>& resolvables,
    const unordered_set<uint64_t>& unresolvables)
{
    (void)resolvables;
    (void)unresolvables;
    const uint64_t fromSeg = segmentOf(from);
    DagNodeInfo info = nodes[fromSeg];
    info.lengthBp = nodeBpLength(fromSeg);
    info.coverage = getPathCoverage(fromSeg);
    info.removed = false;
    if(!isForward(from)) {
        vector<DagNodeId> oriented;
        oriented.reserve(info.anchorChain.size());
        for(auto it = info.anchorChain.rbegin(); it != info.anchorChain.rend(); ++it) {
            oriented.push_back(*it ^ 1);
        }
        info.anchorChain = std::move(oriented);
    }
    const uint64_t newSeg = addNode(info);
    addEdge(fwdNodeId(newSeg), to, getBpOverlap(from, to));
    return newSeg;
}

uint64_t DirectedAnchorGraph::createEdgeNode(
    DagNodeId from,
    DagNodeId to,
    const unordered_set<uint64_t>& resolvables,
    const unordered_set<uint64_t>& unresolvables)
{
    const uint64_t fromSeg = segmentOf(from);
    const uint64_t toSeg = segmentOf(to);
    DagNodeInfo info = nodes[fromSeg];
    info.lengthBp = nodeBpLength(fromSeg);
    info.coverage = getPathCoverage(fromSeg);
    info.removed = false;
    if(!isForward(from)) {
        vector<DagNodeId> oriented;
        oriented.reserve(info.anchorChain.size());
        for(auto it = info.anchorChain.rbegin(); it != info.anchorChain.rend(); ++it) {
            oriented.push_back(*it ^ 1);
        }
        info.anchorChain = std::move(oriented);
    }

    const uint64_t overlap = getBpOverlap(from, to);
    const bool toIsResolvable =
        resolvables.count(toSeg) == 1 && unresolvables.count(toSeg) == 0;

    if(toIsResolvable) {
        // Exact merge length: L1 + L2 - overlap.
        // We already have L1 (info.lengthBp). Add L2 - overlap.
        if (nodeBpLength(toSeg) > overlap) {
            info.lengthBp += (nodeBpLength(toSeg) - overlap);
        }

        if(isForward(to)) {
            info.anchorChain.insert(
                info.anchorChain.end(),
                nodes[toSeg].anchorChain.begin(),
                nodes[toSeg].anchorChain.end());
        } else {
            for(auto it = nodes[toSeg].anchorChain.rbegin();
                it != nodes[toSeg].anchorChain.rend(); ++it) {
                info.anchorChain.push_back(*it ^ 1);
            }
        }
    } else {
        // MBG logic: if destination is unresolvable, we extend by 1 unit (kmer/bp)
        // to represent the step into the node, but do not satisfy the full length.
        info.lengthBp += 1;
    }

    const uint64_t newSeg = addNode(info);
    if(!toIsResolvable) {
        addEdge(fwdNodeId(newSeg), to, overlap + 1);
    }
    return newSeg;
}

vector<DagTriplet> DirectedAnchorGraph::getRawTriplets(
    uint64_t segId,
    uint64_t minEdgeSupport,
    bool partTriplets) const
{
    if(!nodeExists(segId)) return {};

    // Hash for pair<DagNodeId, DagNodeId>.
    struct PairHash {
        size_t operator()(pair<DagNodeId, DagNodeId> p) const {
            return hash<uint64_t>()(p.first) ^ (hash<uint64_t>()(p.second) << 32);
        }
    };

    unordered_map<pair<DagNodeId, DagNodeId>, uint64_t, PairHash> tripletCounts;

    const auto& crossingPaths = pathsCrossing.getPathsCrossing(segId);
    for(const auto& occ : crossingPaths) {
        const uint64_t pathIdx = occ.pathIdx;
        const size_t i = occ.offset;

        if(pathIdx >= paths.size() || pathRemoved[pathIdx]) {
            continue;
        }
        const auto& path = paths[pathIdx];
        if(path.size() < 2) continue;
        const bool haveReads = pathIdx < pathReadIds.size() && !pathReadIds[pathIdx].empty();
        const uint64_t fallbackWeight = pathGroupWeight(pathIdx);
        if(!haveReads && fallbackWeight == 0) continue;

        // occ logic: path[i] is this segment.
        if(segmentOf(path[i]) != segId) continue; // Should be redundant but safe.

            DagNodeId from = dagTipSentinel;
            DagNodeId to = dagTipSentinel;
            uint64_t coverageHere = 0;

            // Start tip.
            if(i == 0 && outDegree(rcNode(path[0])) == 0) {
                if(path[0] == fwdNodeId(segId)) {
                    if(i + 1 < path.size()) {
                        to = path[1];
                    }
                } else {
                    if(i + 1 < path.size()) {
                        from = rcNode(path[1]);
                    }
                }
                uint64_t shouldCheckTrimAmount = 0;
                if(path.size() == 2) {
                    shouldCheckTrimAmount = getTrimAmountToCheck(path[0], path[1]);
                }
                if(!haveReads) {
                    coverageHere = fallbackWeight;
                } else if(shouldCheckTrimAmount == 0) {
                    coverageHere = pathReadIds[pathIdx].size();
                } else {
                    for(const auto& read : pathReadIds[pathIdx]) {
                        if(read.rightClip < shouldCheckTrimAmount) {
                            coverageHere++;
                        }
                    }
                }
            }

            // End tip.
            if(i + 1 == path.size() && outDegree(path.back()) == 0) {
                from = dagTipSentinel;
                to = dagTipSentinel;
                if(path.back() == fwdNodeId(segId)) {
                    if(path.size() >= 2) {
                        from = path[path.size() - 2];
                    }
                } else {
                    if(path.size() >= 2) {
                        to = rcNode(path[path.size() - 2]);
                    }
                }
                uint64_t shouldCheckTrimAmount = 0;
                if(path.size() == 2) {
                    shouldCheckTrimAmount =
                        getTrimAmountToCheck(rcNode(path[1]), rcNode(path[0]));
                }
                uint64_t cov = 0;
                if(!haveReads) {
                    cov = fallbackWeight;
                } else if(shouldCheckTrimAmount == 0) {
                    cov = pathReadIds[pathIdx].size();
                } else {
                    for(const auto& read : pathReadIds[pathIdx]) {
                        if(read.leftClip < shouldCheckTrimAmount) {
                            cov++;
                        }
                    }
                }
                if(cov > coverageHere) {
                    coverageHere = cov;
                }
            }

            // Internal occurrence.
            if(i > 0 && i + 1 < path.size()) {
                if(path[i] == fwdNodeId(segId)) {
                    from = path[i - 1];
                    to = path[i + 1];
                } else {
                    from = rcNode(path[i + 1]);
                    to = rcNode(path[i - 1]);
                }

                uint64_t shouldCheckLeftTrimAmount = 0;
                uint64_t shouldCheckRightTrimAmount = 0;
                if(i == 1) {
                    shouldCheckLeftTrimAmount =
                        getTrimAmountToCheck(rcNode(path[1]), rcNode(path[0]));
                }
                if(i + 2 == path.size()) {
                    shouldCheckRightTrimAmount =
                        getTrimAmountToCheck(path[i], path[i + 1]);
                }

                uint64_t cov = 0;
                if(!haveReads) {
                    cov = fallbackWeight;
                } else if(shouldCheckLeftTrimAmount == 0 &&
                          shouldCheckRightTrimAmount == 0) {
                    cov = pathReadIds[pathIdx].size();
                } else {
                    for(const auto& read : pathReadIds[pathIdx]) {
                        if((shouldCheckLeftTrimAmount == 0 ||
                            read.leftClip < shouldCheckLeftTrimAmount) &&
                           (shouldCheckRightTrimAmount == 0 ||
                            read.rightClip < shouldCheckRightTrimAmount)) {
                            cov++;
                        }
                    }
                }
                if(cov > coverageHere) {
                    coverageHere = cov;
                }
            }

            if(coverageHere > 0) {
                tripletCounts[{from, to}] += coverageHere;
            }
        }
    unordered_map<pair<DagNodeId, DagNodeId>, uint64_t, PairHash> partTripletCounts =
        tripletCounts;

    const bool canTryPartTriplets =
        partTriplets &&
        outDegree(fwdNodeId(segId)) >= 1 &&
        outDegree(revNodeId(segId)) >= 1;
    if(canTryPartTriplets) {
        // readNameIndex -> (pathIdx, (readPosZeroOffset, readPosBoundaryIndex)).
        using PathReadPos = pair<uint64_t, pair<uint64_t, uint64_t>>;
        unordered_map<uint64_t, PathReadPos> fwReads;
        unordered_map<uint64_t, PathReadPos> bwReads;
        const uint64_t invalidPath = numeric_limits<uint64_t>::max();

        for(const auto& occ : crossingPaths) {
            uint64_t pathIdx = occ.pathIdx;
            if(pathIdx >= paths.size() || pathRemoved[pathIdx] ||
               pathIdx >= pathReadIds.size()) {
                continue;
            }
            const auto& path = paths[pathIdx];
            if(path.size() < 2 || pathReadIds[pathIdx].empty()) continue;
            if(segmentOf(path.front()) == segId &&
               segmentOf(path.back()) == segId) {
                continue;
            }

            uint64_t pos = occ.offset;
            // Derived from occ.offset.
            if(pos == 0) {
                for(const auto& read : pathReadIds[pathIdx]) {
                    auto it = fwReads.find(read.readNameIndex);
                    if(it != fwReads.end()) {
                        it->second = {invalidPath, {0, 0}};
                    } else {
                        fwReads[read.readNameIndex] = {
                            pathIdx,
                            {read.readPosZeroOffset, read.readPosStartIndex}
                        };
                    }
                }
            }
            if(pos + 1 == path.size()) {
                for(const auto& read : pathReadIds[pathIdx]) {
                    auto it = bwReads.find(read.readNameIndex);
                    if(it != bwReads.end()) {
                        it->second = {invalidPath, {0, 0}};
                    } else {
                        bwReads[read.readNameIndex] = {
                            pathIdx,
                            {read.readPosZeroOffset, read.readPosEndIndex}
                        };
                    }
                }
            }
        }

        unordered_set<pair<uint64_t, uint64_t>, PairHash> pathPairs;
        for(const auto& entry : fwReads) {
            const uint64_t readName = entry.first;
            const auto& fw = entry.second;
            if(fw.first == invalidPath) continue;
            auto bwIt = bwReads.find(readName);
            if(bwIt == bwReads.end()) continue;
            const auto& bw = bwIt->second;
            if(bw.first == invalidPath) continue;
            if(fw.first == bw.first) continue;

            // MBG ordering check:
            // bw.zeroOffset <= fw.zeroOffset, and if equal then bw.endPos <= fw.startPos.
            if(bw.second.first > fw.second.first) continue;
            if(bw.second.first == fw.second.first &&
               bw.second.second > fw.second.second) {
                continue;
            }
            if(paths[fw.first].empty() || paths[bw.first].empty()) continue;
            if(paths[fw.first].front() != paths[bw.first].back()) continue;

            pathPairs.emplace(bw.first, fw.first);
        }

        for(const auto& pair : pathPairs) {
            const uint64_t bwPathIdx = pair.first;
            const uint64_t fwPathIdx = pair.second;
            if(bwPathIdx >= paths.size() || fwPathIdx >= paths.size()) continue;
            if(pathRemoved[bwPathIdx] || pathRemoved[fwPathIdx]) continue;
            if(bwPathIdx >= pathReadIds.size() || fwPathIdx >= pathReadIds.size()) continue;
            const auto& bwPath = paths[bwPathIdx];
            const auto& fwPath = paths[fwPathIdx];
            if(bwPath.size() < 2 || fwPath.size() < 2) continue;
            if(bwPath.back() != fwPath.front()) continue;

            const uint64_t nodeLen = nodeBpLength(segId);
            uint64_t minLeftClip = nodeLen;
            uint64_t minRightClip = nodeLen;
            for(const auto& read : pathReadIds[bwPathIdx]) {
                minLeftClip = min(minLeftClip, read.rightClip);
            }
            for(const auto& read : pathReadIds[fwPathIdx]) {
                minRightClip = min(minRightClip, read.leftClip);
            }

            unordered_set<uint64_t> fwReadNames;
            fwReadNames.reserve(pathReadIds[fwPathIdx].size() * 2 + 1);
            for(const auto& read : pathReadIds[fwPathIdx]) {
                fwReadNames.insert(read.readNameIndex);
            }

            uint64_t sharedCoverage = 0;
            for(const auto& read : pathReadIds[bwPathIdx]) {
                if(fwReadNames.count(read.readNameIndex) == 1) {
                    sharedCoverage++;
                }
            }
            if(sharedCoverage == 0) continue;

            if(minLeftClip + minRightClip < nodeLen) {
                if(isForward(bwPath.back())) {
                    partTripletCounts[{bwPath[bwPath.size() - 2], fwPath[1]}] += sharedCoverage;
                } else {
                    partTripletCounts[{
                        rcNode(fwPath[1]),
                        rcNode(bwPath[bwPath.size() - 2])
                    }] += sharedCoverage;
                }
            }
        }


    }

    vector<DagTriplet> supported;
    for(const auto& [key, count] : tripletCounts) {
        if(count >= minEdgeSupport) {
            DagTriplet t;
            t.from = key.first;
            t.through = fwdNodeId(segId);
            t.to = key.second;
            t.support = count;
            supported.push_back(t);
        }
    }

    if(canTryPartTriplets) {
        bool canAddPartTriplets = true;
        unordered_set<DagNodeId> uniqueLeft;
        unordered_set<DagNodeId> uniqueRight;

        for(const auto& [key, count] : partTripletCounts) {
            if(count < minEdgeSupport) continue;
            const DagNodeId left = key.first;
            const DagNodeId right = key.second;
            if(left == dagTipSentinel || right == dagTipSentinel) {
                canAddPartTriplets = false;
                break;
            }
            if(uniqueLeft.count(left) == 1 || uniqueRight.count(right) == 1) {
                canAddPartTriplets = false;
                break;
            }
            if(outDegree(left) != 1) {
                canAddPartTriplets = false;
                break;
            }
            if(outDegree(rcNode(right)) != 1) {
                canAddPartTriplets = false;
                break;
            }
            uniqueLeft.insert(left);
            uniqueRight.insert(right);
        }

        if(canAddPartTriplets) {
            supported.clear();
            for(const auto& [key, count] : partTripletCounts) {
                if(count >= minEdgeSupport) {
                    DagTriplet t;
                    t.from = key.first;
                    t.through = fwdNodeId(segId);
                    t.to = key.second;
                    t.support = count;
                    supported.push_back(t);
                }
            }
        }
    }
    return supported;
}

vector<DagTriplet> DirectedAnchorGraph::getReadSupportedTriplets(
    const unordered_set<uint64_t>& resolvables,
    uint64_t segId,
    uint64_t minEdgeSupport,
    bool unconditional,
    bool guesswork) const
{
    (void)resolvables;
    const vector<DagTriplet> coveredTriplets =
        getRawTriplets(segId, minEdgeSupport, guesswork);
    if(unconditional && coveredTriplets.size() >= 2) {
        return coveredTriplets;
    }
    unordered_set<DagNodeId> coveredInNeighbors;
    unordered_set<DagNodeId> coveredOutNeighbors;
    for(const auto& t : coveredTriplets) {
        if(!t.isFromTip()) {
            coveredInNeighbors.insert(rcNode(t.from));
        }
        if(!t.isToTip()) {
            coveredOutNeighbors.insert(t.to);
        }
    }
    if(coveredInNeighbors.size() < outDegree(revNodeId(segId))) {
        return {};
    }
    if(coveredOutNeighbors.size() < outDegree(fwdNodeId(segId))) {
        return {};
    }
    return coveredTriplets;
}

vector<DagTriplet> DirectedAnchorGraph::getGuessworkTriplets(
    const unordered_set<uint64_t>& resolvables,
    uint64_t segId,
    uint64_t minEdgeSupport,
    bool unconditional,
    double averageCoverage) const
{
    (void)unconditional;
    if(!nodeExists(segId) || averageCoverage <= 0.0) {
        return {};
    }
    const uint64_t fwDegree = outDegree(fwdNodeId(segId));
    const uint64_t bwDegree = outDegree(revNodeId(segId));
    if(fwDegree == 1 || bwDegree == 1) return {};

    const double nodeCov = getPathCoverage(segId);
    const uint64_t nodeCopyCount =
        uint64_t((nodeCov + averageCoverage / 2.0) / averageCoverage);
    if(nodeCopyCount < 2) return {};

    struct NeighborCopy {
        DagNodeId oriented;
        uint64_t copyCount;
    };
    uint64_t outneighborCopyCountSum = 0;
    uint64_t inneighborCopyCountSum = 0;
    vector<NeighborCopy> outneighborCopyCounts;
    vector<NeighborCopy> inneighborCopyCounts;

    for(DagNodeId edge : getOutEdges(fwdNodeId(segId))) {
        const double cov = getPathCoverage(segmentOf(edge));
        const uint64_t copyCount = uint64_t((cov + averageCoverage / 2.0) / averageCoverage);
        outneighborCopyCountSum += copyCount;
        outneighborCopyCounts.push_back({edge, copyCount});
        if(resolvables.count(segmentOf(edge)) == 1) return {};
    }
    for(DagNodeId edge : getOutEdges(revNodeId(segId))) {
        const double cov = getPathCoverage(segmentOf(edge));
        const uint64_t copyCount = uint64_t((cov + averageCoverage / 2.0) / averageCoverage);
        inneighborCopyCountSum += copyCount;
        inneighborCopyCounts.push_back({rcNode(edge), copyCount});
        if(resolvables.count(segmentOf(edge)) == 1) return {};
    }
    if(outneighborCopyCountSum != nodeCopyCount) return {};
    if(inneighborCopyCountSum != nodeCopyCount) return {};

    auto coveredTriplets = getRawTriplets(segId, minEdgeSupport, true);
    if(coveredTriplets.empty()) return {};
    if(coveredTriplets.size() == 1 && minEdgeSupport == 1) {
        const auto test = getRawTriplets(segId, 2, true);
        if(test.empty()) return {};
    }

    unordered_set<DagNodeId> coveredInNeighbors;
    unordered_set<DagNodeId> coveredOutNeighbors;
    for(const auto& t : coveredTriplets) {
        if(!t.isFromTip()) coveredInNeighbors.insert(rcNode(t.from));
        if(!t.isToTip()) coveredOutNeighbors.insert(t.to);
    }
    if(outDegree(revNodeId(segId)) >= coveredInNeighbors.size() + 2 ||
       outDegree(fwdNodeId(segId)) >= coveredOutNeighbors.size() + 2) {
        return {};
    }

    uint64_t uncoveredOutCopycounts = 0;
    uint64_t uncoveredInCopycounts = 0;
    for(const auto& p : outneighborCopyCounts) {
        if(coveredOutNeighbors.count(p.oriented) == 0) uncoveredOutCopycounts += p.copyCount;
    }
    for(const auto& p : inneighborCopyCounts) {
        if(coveredInNeighbors.count(p.oriented) == 0) uncoveredInCopycounts += p.copyCount;
    }
    if(uncoveredInCopycounts != uncoveredOutCopycounts) return {};
    if(uncoveredOutCopycounts == 0) return coveredTriplets;
    if(uncoveredOutCopycounts >= 4) return {};

    uint64_t shortestOutAnchor = nodeBpLength(segId);
    uint64_t shortestInAnchor = nodeBpLength(segId);
    for(const auto& outpair : outneighborCopyCounts) {
        if(coveredOutNeighbors.count(outpair.oriented) == 1) continue;
        if(outpair.copyCount == 0) continue;
        const uint64_t anchorSize = getAnchorSize(
            rcNode(outpair.oriented),
            revNodeId(segId));
        shortestOutAnchor = min(shortestOutAnchor, anchorSize);
    }
    for(const auto& inpair : inneighborCopyCounts) {
        if(coveredInNeighbors.count(inpair.oriented) == 1) continue;
        if(inpair.copyCount == 0) continue;
        const uint64_t anchorSize = getAnchorSize(
            inpair.oriented,
            fwdNodeId(segId));
        shortestInAnchor = min(shortestInAnchor, anchorSize);
    }
    if(shortestInAnchor + shortestOutAnchor <= nodeBpLength(segId)) {
        return {};
    }

    uint64_t addedGuesses = 0;
    for(const auto& inpair : inneighborCopyCounts) {
        if(inpair.copyCount == 0) continue;
        if(coveredInNeighbors.count(inpair.oriented) == 1) continue;
        for(const auto& outpair : outneighborCopyCounts) {
            if(outpair.copyCount == 0) continue;
            if(coveredOutNeighbors.count(outpair.oriented) == 1) continue;
            if(segmentOf(inpair.oriented) == segmentOf(outpair.oriented)) return {};
            DagTriplet t;
            t.from = inpair.oriented;
            t.through = fwdNodeId(segId);
            t.to = outpair.oriented;
            t.support = 0;
            coveredTriplets.push_back(t);
            addedGuesses++;
        }
    }
    if(addedGuesses == 0) {
        return {};
    }
    return coveredTriplets;
}

void DirectedAnchorGraph::filterCopyCountTriplets(
    uint64_t segId,
    vector<DagTriplet>& originals,
    double averageCoverage) const
{
    if(!nodeExists(segId) || originals.empty() || averageCoverage <= 0.0) return;
    const double coverage = getPathCoverage(segId);
    const uint64_t copyCount = uint64_t(coverage / averageCoverage + 0.5);
    if(originals.size() <= copyCount || copyCount == 0) return;

    sort(originals.begin(), originals.end(),
         [](const DagTriplet& left, const DagTriplet& right) {
             return left.support > right.support;
         });

    unordered_set<DagNodeId> inneighborNeedsCovering;
    unordered_set<DagNodeId> outneighborNeedsCovering;
    unordered_set<DagNodeId> inneighborDoesntNeedCovering;
    unordered_set<DagNodeId> outneighborDoesntNeedCovering;
    for(DagNodeId edge : getOutEdges(fwdNodeId(segId))) {
        if(getPathCoverage(segmentOf(edge)) > averageCoverage * 0.25) {
            outneighborNeedsCovering.insert(edge);
        } else {
            outneighborDoesntNeedCovering.insert(edge);
        }
    }
    for(DagNodeId edge : getOutEdges(revNodeId(segId))) {
        DagNodeId pred = rcNode(edge);
        if(getPathCoverage(segmentOf(edge)) > averageCoverage * 0.25) {
            inneighborNeedsCovering.insert(pred);
        } else {
            inneighborDoesntNeedCovering.insert(pred);
        }
    }

    const bool maybeDontRemoveLast =
        (inneighborDoesntNeedCovering.size() == 1) &&
        (outneighborDoesntNeedCovering.size() == 1);
    uint64_t threshold = 0;
    const uint64_t maxI = min<uint64_t>(copyCount, originals.size());
    for(uint64_t i = 0; i + 1 < originals.size() && i < maxI; i++) {
        if(!originals[i].isFromTip()) {
            inneighborNeedsCovering.erase(originals[i].from);
            inneighborDoesntNeedCovering.erase(originals[i].from);
        }
        if(!originals[i].isToTip()) {
            outneighborNeedsCovering.erase(originals[i].to);
            outneighborDoesntNeedCovering.erase(originals[i].to);
        }
        if(!inneighborNeedsCovering.empty()) continue;
        if(!outneighborNeedsCovering.empty()) continue;
        if(!(double(originals[i].support) >= double(originals[i + 1].support) * 4.0)) continue;
        if(double(originals[i + 1].support) >= averageCoverage * 0.25) continue;
        if(double(originals[i].support) < averageCoverage * 0.5) break;
        threshold = i + 1;
        break;
    }
    if(threshold == 0 || threshold >= originals.size()) return;
    if(maybeDontRemoveLast && threshold == originals.size() - 1 &&
       inneighborDoesntNeedCovering.size() == 1 &&
       outneighborDoesntNeedCovering.size() == 1) {
        if(!originals[threshold].isFromTip() &&
           !originals[threshold].isToTip() &&
           originals[threshold].from == *inneighborDoesntNeedCovering.begin() &&
           originals[threshold].to == *outneighborDoesntNeedCovering.begin()) {
            return;
        }
    }
    originals.erase(originals.begin() + threshold, originals.end());
}

vector<DagTriplet> DirectedAnchorGraph::getValidTripletsForResolve(
    const unordered_set<uint64_t>& resolvables,
    uint64_t segId,
    uint64_t minEdgeSupport,
    bool unconditional,
    bool guesswork,
    bool copycountFilterHeuristic,
    double averageCoverage) const
{
    auto triplets = getReadSupportedTriplets(
        resolvables, segId, minEdgeSupport, unconditional, guesswork);
    if(!triplets.empty() && guesswork && unconditional && copycountFilterHeuristic) {
        filterCopyCountTriplets(segId, triplets, averageCoverage);
    } else if(triplets.empty() && guesswork) {
        triplets = getGuessworkTriplets(
            resolvables, segId, minEdgeSupport, unconditional, averageCoverage);
    }
    sort(triplets.begin(), triplets.end(),
         [](const DagTriplet& left, const DagTriplet& right) {
             if(left.from < right.from) return true;
             if(left.from > right.from) return false;
             if(left.to < right.to) return true;
             if(left.to > right.to) return false;
             return false;
         });
    return triplets;
}

vector<DagTriplet> DirectedAnchorGraph::getValidTriplets(
    uint64_t segId,
    uint64_t minEdgeSupport) const
{
    return getRawTriplets(segId, minEdgeSupport, false);
}


// ============================================================================
// resolveHairpins — handle palindromic segments with self-RC edges.
// ============================================================================

unordered_set<uint64_t> DirectedAnchorGraph::resolveHairpins(
    const vector<uint64_t>& candidates)
{
    unordered_set<uint64_t> newSegs;

    for(uint64_t segId : candidates) {
        if(!nodeExists(segId)) continue;

        DagNodeId fwd = fwdNodeId(segId);
        DagNodeId rev = revNodeId(segId);

        // Check for hairpin: edge fwd→rev or rev→fwd.
        const auto& fwdEdges = getOutEdges(fwd);
        const auto& revEdges = getOutEdges(rev);
        bool fwdToRev = find(fwdEdges.begin(), fwdEdges.end(), rev) != fwdEdges.end();
        bool revToFwd = find(revEdges.begin(), revEdges.end(), fwd) != revEdges.end();

        if(!fwdToRev && !revToFwd) continue;
        if(fwdToRev && revToFwd) continue;  // double hairpin — skip

        // Create forward and backward copies.
        DagNodeInfo fwdInfo = nodes[segId];
        DagNodeInfo bwdInfo = fwdInfo;

        uint64_t fwdCopy = addNode(fwdInfo);
        uint64_t bwdCopy = addNode(bwdInfo);
        newSegs.insert(fwdCopy);
        newSegs.insert(bwdCopy);

        // Copy non-hairpin edges from original fwd to fwdCopy fwd.
        vector<DagNodeId> fwdCopy_edges = getOutEdges(fwd);
        for(DagNodeId to : fwdCopy_edges) {
            if(segmentOf(to) == segId) continue;
            addEdge(fwdNodeId(fwdCopy), to, getBpOverlap(fwd, to));
        }

        // Copy non-hairpin edges from original rev to bwdCopy fwd.
        vector<DagNodeId> bwdCopy_edges = getOutEdges(rev);
        for(DagNodeId to : bwdCopy_edges) {
            if(segmentOf(to) == segId) continue;
            addEdge(fwdNodeId(bwdCopy), to, getBpOverlap(rev, to));
        }

        // Add the hairpin edge between copies.
        if(fwdToRev) {
            const uint64_t ov = (nodeBpLength(segId) > 0) ?
                (nodeBpLength(segId) - 1) : 0;
            addEdge(fwdNodeId(fwdCopy), revNodeId(bwdCopy), ov);
        } else {
            const uint64_t ov = (nodeBpLength(segId) > 0) ?
                (nodeBpLength(segId) - 1) : 0;
            addEdge(fwdNodeId(bwdCopy), revNodeId(fwdCopy), ov);
        }

        // Rewrite paths: forward strand → fwdCopy, reverse strand → bwdCopy.
        const auto& crossingPaths = pathsCrossing.getPathsCrossing(segId);
        
        static thread_local vector<uint64_t> pathsToUpdate;
        pathsToUpdate.clear();
        pathsToUpdate.reserve(crossingPaths.size());
        for(const auto& occ : crossingPaths) {
            pathsToUpdate.push_back(occ.pathIdx);
        }
        sort(pathsToUpdate.begin(), pathsToUpdate.end());
        pathsToUpdate.erase(unique(pathsToUpdate.begin(), pathsToUpdate.end()), pathsToUpdate.end());
        
        for(uint64_t pathIdx : pathsToUpdate) {
            if(pathRemoved[pathIdx]) continue;
            auto& path = paths[pathIdx];
            pathsCrossing.removePath(pathIdx, path);
            for(auto& n : path) {
                if(segmentOf(n) == segId) {
                    if(isForward(n)) {
                        n = fwdNodeId(fwdCopy);
                    } else {
                        n = fwdNodeId(bwdCopy);
                    }
                }
            }
            pathsCrossing.addPath(pathIdx, path);
        }

        removeNode(segId);
    }

    return newSegs;
}


// ============================================================================
// resolveNodes — create edge nodes for each valid triplet.
// ============================================================================




DirectedAnchorGraph::ResolveStats DirectedAnchorGraph::resolveNodes(
    const vector<uint64_t>& candidates,
    uint64_t minEdgeSupport,
    bool unconditional,
    bool copycountFilterHeuristic,
    bool guesswork,
    double averageCoverage)
{
    ResolveStats result;
    if(averageCoverage <= 0.0) averageCoverage = 1.0;

    vector<uint64_t> orderedCandidates = candidates;
    sort(orderedCandidates.begin(), orderedCandidates.end());
    const unordered_set<uint64_t> resolvablesSet(
        orderedCandidates.begin(), orderedCandidates.end());
    unordered_set<uint64_t> unresolvables;

    auto isExactPalindromeFrozen = [&](uint64_t segId) {
        if(!nodeExists(segId)) return false;
        const DagNodeId fw = fwdNodeId(segId);
        const DagNodeId bw = revNodeId(segId);
        const auto& fwEdges = getOutEdges(fw);
        const auto& bwEdges = getOutEdges(bw);
        return find(fwEdges.begin(), fwEdges.end(), bw) != fwEdges.end() ||
            find(bwEdges.begin(), bwEdges.end(), fw) != bwEdges.end();
    };

    for(uint64_t segId : orderedCandidates) {
        if(!nodeExists(segId)) continue;
        if(isExactPalindromeFrozen(segId)) {
            unresolvables.insert(segId);
        }
    }

    for(uint64_t segId : orderedCandidates) {
        if(!nodeExists(segId) || unresolvables.count(segId) == 1) continue;
        auto triplets = getValidTripletsForResolve(
            resolvablesSet,
            segId,
            minEdgeSupport,
            unconditional,
            guesswork,
            copycountFilterHeuristic,
            averageCoverage);
        if(triplets.empty()) {
            unresolvables.insert(segId);
        }
    }

    vector<uint64_t> check = orderedCandidates;
    while(true) {
        unordered_set<uint64_t> newCheck;
        bool anyChanged = false;
        for(uint64_t segId : check) {
            if(!nodeExists(segId) || unresolvables.count(segId) == 1) continue;
            auto triplets = getValidTripletsForResolve(
                resolvablesSet,
                segId,
                minEdgeSupport,
                unconditional,
                guesswork,
                copycountFilterHeuristic,
                averageCoverage);

            bool unresolve = false;
            for(const auto& t : triplets) {
                bool bwFake = false;
                bool fwFake = false;

                if(!t.isFromTip()) {
                    const uint64_t leftSeg = segmentOf(t.from);
                    if((resolvablesSet.count(leftSeg) == 0 ||
                        unresolvables.count(leftSeg) == 1) &&
                       nodeBpLength(leftSeg) ==
                       getBpOverlap(revNodeId(segId), rcNode(t.from)) + 1) {
                        bwFake = true;
                    }
                }
                if(!t.isToTip()) {
                    const uint64_t rightSeg = segmentOf(t.to);
                    if((resolvablesSet.count(rightSeg) == 0 ||
                        unresolvables.count(rightSeg) == 1) &&
                       nodeBpLength(rightSeg) ==
                       getBpOverlap(fwdNodeId(segId), t.to) + 1) {
                        fwFake = true;
                    }
                }
                if(bwFake && !t.isFromTip() && outDegree(t.from) >= 2) {
                    unresolve = true;
                }
                if(fwFake && !t.isToTip() && outDegree(rcNode(t.to)) >= 2) {
                    unresolve = true;
                }
                if(bwFake && fwFake) {
                    unresolve = true;
                }
                if(unresolve) break;
            }

            if(unresolve) {
                unresolvables.insert(segId);
                anyChanged = true;
                for(DagNodeId edge : getOutEdges(fwdNodeId(segId))) {
                    const uint64_t nSeg = segmentOf(edge);
                    if(resolvablesSet.count(nSeg) == 1 &&
                       unresolvables.count(nSeg) == 0) {
                        newCheck.insert(nSeg);
                    }
                }
                for(DagNodeId edge : getOutEdges(revNodeId(segId))) {
                    const uint64_t nSeg = segmentOf(edge);
                    if(resolvablesSet.count(nSeg) == 1 &&
                       unresolvables.count(nSeg) == 0) {
                        newCheck.insert(nSeg);
                    }
                }
            }
        }
        if(!anyChanged || newCheck.empty()) break;
        check.assign(newCheck.begin(), newCheck.end());
        sort(check.begin(), check.end());
    }

    unordered_set<uint64_t> actuallyResolvables;
    for(uint64_t segId : orderedCandidates) {
        if(!nodeExists(segId)) continue;
        if(unresolvables.count(segId) == 0) {
            actuallyResolvables.insert(segId);
        }
    }
    if(actuallyResolvables.empty()) return result;

    unordered_map<uint64_t, vector<DagTriplet>> tripletsPerNode;
    for(uint64_t segId : orderedCandidates) {
        if(actuallyResolvables.count(segId) == 0) continue;
        tripletsPerNode[segId] = getValidTripletsForResolve(
            resolvablesSet,
            segId,
            minEdgeSupport,
            unconditional,
            guesswork,
            copycountFilterHeuristic,
            averageCoverage);
    }

    unordered_map<pair<DagNodeId, DagNodeId>, uint64_t, DagEdgePairHash> newEdgeNodes;
    for(uint64_t segId : orderedCandidates) {
        if(actuallyResolvables.count(segId) == 0) continue;
        const auto& triplets = tripletsPerNode[segId];
        if(triplets.empty()) continue;

        unordered_set<DagNodeId> fwCovered;
        unordered_set<DagNodeId> bwCovered;
        for(const auto& t : triplets) {
            if(!t.isToTip()) {
                fwCovered.insert(t.to);
            }
            if(!t.isFromTip()) {
                bwCovered.insert(rcNode(t.from));
            }
        }

        vector<DagNodeId> fwEdges = getOutEdges(fwdNodeId(segId));
        vector<DagNodeId> bwEdges = getOutEdges(revNodeId(segId));
        sort(fwEdges.begin(), fwEdges.end());
        sort(bwEdges.begin(), bwEdges.end());

        for(DagNodeId edge : fwEdges) {
            const auto key = make_pair(fwdNodeId(segId), edge);
            const auto rcKey = make_pair(rcNode(edge), rcNode(fwdNodeId(segId)));
            if(fwCovered.count(edge) == 0) {
                result.maybeUnitigifiable.insert(segmentOf(edge));
                result.maybeTrimmable[rcNode(edge)] = max(
                    result.maybeTrimmable[rcNode(edge)],
                    getBpOverlap(fwdNodeId(segId), edge));
                continue;
            }
            if(newEdgeNodes.count(rcKey) == 1) {
                continue;
            }
            const bool isFake =
                (resolvablesSet.count(segmentOf(edge)) == 0 ||
                 unresolvables.count(segmentOf(edge)) == 1) &&
                nodeBpLength(segmentOf(edge)) ==
                getBpOverlap(fwdNodeId(segId), edge) + 1;
            const uint64_t newSeg = isFake ?
                createFakeEdgeNode(fwdNodeId(segId), edge, resolvablesSet, unresolvables) :
                createEdgeNode(fwdNodeId(segId), edge, resolvablesSet, unresolvables);
            newEdgeNodes[key] = newSeg;
            result.newSegs.insert(newSeg);
            result.maybeUnitigifiable.insert(newSeg);
            if(bwEdges.empty()) {
                result.maybeTrimmable[revNodeId(newSeg)] = max(
                    result.maybeTrimmable[revNodeId(newSeg)],
                    nodeBpLength(segId));
            }
        }

        for(DagNodeId edge : bwEdges) {
            const auto key = make_pair(revNodeId(segId), edge);
            const auto rcKey = make_pair(rcNode(edge), rcNode(revNodeId(segId)));
            if(bwCovered.count(edge) == 0) {
                result.maybeUnitigifiable.insert(segmentOf(edge));
                result.maybeTrimmable[rcNode(edge)] = max(
                    result.maybeTrimmable[rcNode(edge)],
                    getBpOverlap(revNodeId(segId), edge));
                continue;
            }
            if(newEdgeNodes.count(rcKey) == 1) {
                continue;
            }
            const bool isFake =
                (resolvablesSet.count(segmentOf(edge)) == 0 ||
                 unresolvables.count(segmentOf(edge)) == 1) &&
                nodeBpLength(segmentOf(edge)) ==
                getBpOverlap(revNodeId(segId), edge) + 1;
            const uint64_t newSeg = isFake ?
                createFakeEdgeNode(revNodeId(segId), edge, resolvablesSet, unresolvables) :
                createEdgeNode(revNodeId(segId), edge, resolvablesSet, unresolvables);
            newEdgeNodes[key] = newSeg;
            result.newSegs.insert(newSeg);
            result.maybeUnitigifiable.insert(newSeg);
            if(fwEdges.empty()) {
                result.maybeTrimmable[revNodeId(newSeg)] = max(
                    result.maybeTrimmable[revNodeId(newSeg)],
                    nodeBpLength(segId));
            }
        }
    }

    if(unconditional) {
        for(const auto& entry : newEdgeNodes) {
            const DagNodeId fromNode = entry.first.first;
            const DagNodeId toNode = entry.first.second;
            const uint64_t edgeSeg = entry.second;
            const uint64_t fromSeg = segmentOf(fromNode);
            const uint64_t toSeg = segmentOf(toNode);
            if(actuallyResolvables.count(fromSeg) == 0 ||
               actuallyResolvables.count(toSeg) == 0) {
                continue;
            }

            bool fromHasTo = false;
            bool toHasFrom = false;
            for(const auto& t : tripletsPerNode[fromSeg]) {
                if(isForward(fromNode) && !t.isToTip() && t.to == toNode) {
                    fromHasTo = true;
                }
                if(!isForward(fromNode) && !t.isFromTip() && rcNode(t.from) == toNode) {
                    fromHasTo = true;
                }
            }
            for(const auto& t : tripletsPerNode[toSeg]) {
                if(isForward(toNode) && !t.isToTip() && t.to == rcNode(fromNode)) {
                    toHasFrom = true;
                }
                if(!isForward(toNode) && !t.isFromTip() &&
                   rcNode(t.from) == rcNode(fromNode)) {
                    toHasFrom = true;
                }
            }
            if(fromHasTo && !toHasFrom) {
                result.maybeTrimmable[fwdNodeId(edgeSeg)] = max(
                    result.maybeTrimmable[fwdNodeId(edgeSeg)],
                    nodeBpLength(toSeg));
            }
            if(toHasFrom && !fromHasTo) {
                result.maybeTrimmable[revNodeId(edgeSeg)] = max(
                    result.maybeTrimmable[revNodeId(edgeSeg)],
                    nodeBpLength(fromSeg));
            }
        }
    }

    for(uint64_t segId : orderedCandidates) {
        if(actuallyResolvables.count(segId) == 0) continue;
        const DagNodeId pos = fwdNodeId(segId);
        const auto& triplets = tripletsPerNode[segId];
        for(const auto& t : triplets) {
            const bool beforeTip = t.isFromTip();
            const bool afterTip = t.isToTip();
            if(beforeTip && afterTip) {
                continue;
            }

            if(beforeTip) {
                const DagNodeId after = t.to;
                auto itA = newEdgeNodes.find({pos, after});
                auto itB = newEdgeNodes.find({rcNode(after), rcNode(pos)});
                if(itA != newEdgeNodes.end()) {
                    result.maybeTrimmable[revNodeId(itA->second)] = max(
                        result.maybeTrimmable[revNodeId(itA->second)],
                        nodeBpLength(segId));
                } else if(itB != newEdgeNodes.end()) {
                    result.maybeTrimmable[fwdNodeId(itB->second)] = max(
                        result.maybeTrimmable[fwdNodeId(itB->second)],
                        nodeBpLength(segId));
                } else {
                    result.maybeTrimmable[rcNode(after)] = max(
                        result.maybeTrimmable[rcNode(after)],
                        nodeBpLength(segmentOf(after)));
                }
                continue;
            }

            if(afterTip) {
                const DagNodeId before = t.from;
                auto itA = newEdgeNodes.find({rcNode(pos), rcNode(before)});
                auto itB = newEdgeNodes.find({before, pos});
                if(itB != newEdgeNodes.end()) {
                    result.maybeTrimmable[fwdNodeId(itB->second)] = max(
                        result.maybeTrimmable[fwdNodeId(itB->second)],
                        nodeBpLength(segId));
                } else if(itA != newEdgeNodes.end()) {
                    result.maybeTrimmable[revNodeId(itA->second)] = max(
                        result.maybeTrimmable[revNodeId(itA->second)],
                        nodeBpLength(segId));
                } else {
                    result.maybeTrimmable[before] = max(
                        result.maybeTrimmable[before],
                        nodeBpLength(segmentOf(before)));
                }
                continue;
            }

            DagNodeId leftNode = t.from;
            DagNodeId rightNode = t.to;

            auto leftA = newEdgeNodes.find({rcNode(pos), rcNode(t.from)});
            auto leftB = newEdgeNodes.find({t.from, pos});
            if(leftA != newEdgeNodes.end()) {
                leftNode = revNodeId(leftA->second);
            } else if(leftB != newEdgeNodes.end()) {
                leftNode = fwdNodeId(leftB->second);
            } else {
                continue;
            }

            auto rightA = newEdgeNodes.find({pos, t.to});
            auto rightB = newEdgeNodes.find({rcNode(t.to), rcNode(pos)});
            if(rightA != newEdgeNodes.end()) {
                rightNode = fwdNodeId(rightA->second);
            } else if(rightB != newEdgeNodes.end()) {
                rightNode = revNodeId(rightB->second);
            } else {
                continue;
            }

            result.maybeUnitigifiable.insert(segmentOf(leftNode));
            result.maybeUnitigifiable.insert(segmentOf(rightNode));
            const auto& out = getOutEdges(leftNode);
            if(find(out.begin(), out.end(), rightNode) == out.end()) {
                addEdge(leftNode, rightNode, nodeBpLength(segId));
            }
        }
    }

    for(uint64_t segId : orderedCandidates) {
        if(actuallyResolvables.count(segId) == 1 && nodeExists(segId)) {
            removeNode(segId);
        }
    }

    replacePathsFromEdgeNodes(actuallyResolvables, newEdgeNodes);
    splitPathsAtBreaks();

    result.nodesResolved = actuallyResolvables.size();
    result.nodesAdded = newEdgeNodes.size();
    return result;
}


// ============================================================================
// replacePathNodes — rewrite paths after resolution.
// ============================================================================

void DirectedAnchorGraph::replacePathNodes(
    const unordered_map<uint64_t,
        unordered_map<DagNodeId,
            unordered_map<DagNodeId, uint64_t>>>& resolutionMap)
{
    for(uint64_t pathIdx = 0; pathIdx < paths.size(); ++pathIdx) {
        if(pathRemoved[pathIdx]) continue;

        auto& oldPath = paths[pathIdx];
        if(oldPath.empty()) continue;
        const vector<DagPathReadSupport> originalSupport =
            (pathIdx < pathReadIds.size()) ? pathReadIds[pathIdx] :
            vector<DagPathReadSupport>{};

        vector<uint64_t> oldStarts(oldPath.size(), 0);
        vector<uint64_t> oldEnds(oldPath.size(), 0);
        for(size_t i = 0; i < oldPath.size(); i++) {
            if(i == 0) {
                oldStarts[i] = 0;
            } else {
                const uint64_t ov = getBpOverlap(oldPath[i - 1], oldPath[i]);
                oldStarts[i] = (oldEnds[i - 1] >= ov) ? (oldEnds[i - 1] - ov) : 0;
            }
            oldEnds[i] = oldStarts[i] + nodeBpLength(segmentOf(oldPath[i]));
        }
        vector<uint64_t> oldBoundaries;
        oldBoundaries.reserve(oldPath.size() + 1);
        oldBoundaries.push_back(0);
        for(uint64_t e : oldEnds) {
            oldBoundaries.push_back(e);
        }

        vector<DagNodeId> newPath;
        vector<uint64_t> newStarts;
        vector<uint64_t> newEnds;
        vector<uint64_t> forcedBreaks;
        bool modified = false;

        for(uint64_t i = 0; i < oldPath.size(); ++i) {
            const uint64_t segId = segmentOf(oldPath[i]);
            auto segIt = resolutionMap.find(segId);
            if(segIt == resolutionMap.end()) {
                newPath.push_back(oldPath[i]);
                newStarts.push_back(oldStarts[i]);
                newEnds.push_back(oldEnds[i]);
                continue;
            }

            modified = true;
            DagNodeId from = (i > 0) ? oldPath[i - 1] : dagTipSentinel;
            DagNodeId to = (i + 1 < oldPath.size()) ? oldPath[i + 1] : dagTipSentinel;

            auto fromIt = segIt->second.find(from);
            if(fromIt == segIt->second.end()) {
                if(!newPath.empty() &&
                   (forcedBreaks.empty() || forcedBreaks.back() != newPath.size())) {
                    forcedBreaks.push_back(newPath.size());
                }
                continue;
            }
            auto toIt = fromIt->second.find(to);
            if(toIt == fromIt->second.end()) {
                if(!newPath.empty() &&
                   (forcedBreaks.empty() || forcedBreaks.back() != newPath.size())) {
                    forcedBreaks.push_back(newPath.size());
                }
                continue;
            }

            newPath.push_back(fwdNodeId(toIt->second));
            newStarts.push_back(oldStarts[i]);
            newEnds.push_back(oldEnds[i]);
        }

        if(!modified) {
            continue;
        }

        pathsCrossing.removePath(pathIdx, oldPath);
        if(newPath.empty()) {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
            continue;
        }

        vector<uint64_t> breakPositions = forcedBreaks;
        for(uint64_t j = 1; j < newPath.size(); j++) {
            const auto& out = getOutEdges(newPath[j - 1]);
            if(find(out.begin(), out.end(), newPath[j]) == out.end()) {
                breakPositions.push_back(j);
            }
        }
        sort(breakPositions.begin(), breakPositions.end());
        breakPositions.erase(
            unique(breakPositions.begin(), breakPositions.end()),
            breakPositions.end());

        auto buildReadsForSegment = [&](
            uint64_t beginIdx,
            uint64_t endIdx) {
            vector<DagPathReadSupport> result;
            if(beginIdx >= endIdx || endIdx > newPath.size()) {
                return result;
            }
            vector<uint64_t> projectedBoundaries;
            projectedBoundaries.reserve((endIdx - beginIdx) + 1);
            projectedBoundaries.push_back(newStarts[beginIdx]);
            for(uint64_t i = beginIdx; i < endIdx; i++) {
                projectedBoundaries.push_back(newEnds[i]);
            }
            for(const auto& read : originalSupport) {
                auto projected = projectReadSupportToBoundaries(
                    read,
                    oldBoundaries,
                    projectedBoundaries);
                if(projected.has_value()) {
                    result.push_back(std::move(*projected));
                }
            }
            return result;
        };

        vector<pair<uint64_t, uint64_t>> segments;
        uint64_t segStartIdx = 0;
        for(uint64_t bp : breakPositions) {
            if(bp > segStartIdx) {
                segments.push_back({segStartIdx, bp});
            }
            segStartIdx = bp;
        }
        if(segStartIdx < newPath.size()) {
            segments.push_back({segStartIdx, uint64_t(newPath.size())});
        }

        if(segments.empty()) {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
            continue;
        }

        {
            const auto [beginIdx, endIdx] = segments.front();
            vector<DagNodeId> firstPath(
                newPath.begin() + beginIdx, newPath.begin() + endIdx);
            vector<DagPathReadSupport> segReads =
                buildReadsForSegment(beginIdx, endIdx);

            paths[pathIdx] = std::move(firstPath);
            pathRemoved[pathIdx] = false;
            if(pathIdx < pathWeights.size()) {
                pathWeights[pathIdx] = max<uint64_t>(1, segReads.size());
            }
            if(pathIdx < pathReadIds.size()) {
                pathReadIds[pathIdx] = std::move(segReads);
            }
            pathsCrossing.addPath(pathIdx, paths[pathIdx]);
        }

        for(size_t s = 1; s < segments.size(); s++) {
            const auto [beginIdx, endIdx] = segments[s];
            vector<DagNodeId> segPath(
                newPath.begin() + beginIdx, newPath.begin() + endIdx);
            vector<DagPathReadSupport> segReads =
                buildReadsForSegment(beginIdx, endIdx);
            addPathButFirstMaybeTrim(std::move(segPath), std::move(segReads));
        }
    }
}

void DirectedAnchorGraph::replacePathsFromEdgeNodes(
    const unordered_set<uint64_t>& actuallyResolvables,
    const unordered_map<pair<DagNodeId, DagNodeId>, uint64_t, DagEdgePairHash>& newEdgeNodes)
{
    for(uint64_t pathIdx = 0; pathIdx < paths.size(); ++pathIdx) {
        if(pathRemoved[pathIdx]) continue;
        const auto& oldPathRef = paths[pathIdx];
        if(oldPathRef.empty()) continue;

        bool relevant = false;
        for(DagNodeId n : oldPathRef) {
            if(actuallyResolvables.count(segmentOf(n)) == 1) {
                relevant = true;
                break;
            }
        }
        if(!relevant) continue;

        const vector<DagNodeId> oldPath = oldPathRef;
        const vector<DagPathReadSupport> originalSupport =
            (pathIdx < pathReadIds.size()) ? pathReadIds[pathIdx] :
            vector<DagPathReadSupport>{};
        const vector<uint64_t> oldBoundaries = computePathNodeBoundaries(oldPath);
        if(oldBoundaries.size() != oldPath.size() + 1) {
            continue;
        }

        vector<DagNodeId> newPath;
        vector<uint64_t> nodePosStarts;
        vector<uint64_t> nodePosEnds;
        vector<uint64_t> breakFromInvalidEdge;
        uint64_t runningKmerStartPos = 0;
        uint64_t runningKmerEndPos = 0;
        newPath.reserve(oldPath.size() * 2);
        nodePosStarts.reserve(oldPath.size() * 2);
        nodePosEnds.reserve(oldPath.size() * 2);

        auto maybePushBreak = [&]() {
            if(newPath.empty()) return;
            if(breakFromInvalidEdge.empty() ||
               breakFromInvalidEdge.back() != newPath.size()) {
                breakFromInvalidEdge.push_back(newPath.size());
            }
        };

        for(size_t j = 0; j < oldPath.size(); j++) {
            const DagNodeId cur = oldPath[j];
            runningKmerStartPos = runningKmerEndPos;
            runningKmerEndPos += nodeBpLength(segmentOf(cur));
            uint64_t overlapPrev = 0;
            if(j > 0) {
                overlapPrev = getBpOverlap(oldPath[j - 1], oldPath[j]);
                runningKmerEndPos =
                    (runningKmerEndPos >= overlapPrev) ?
                    (runningKmerEndPos - overlapPrev) : 0;
            }

            const bool curResolvable =
                actuallyResolvables.count(segmentOf(cur)) == 1;
            if(!curResolvable) {
                uint64_t start = 0;
                if(j > 0) {
                    start = runningKmerStartPos;
                    start = (start >= overlapPrev) ? (start - overlapPrev) : 0;
                }
                newPath.push_back(cur);
                nodePosStarts.push_back(start);
                nodePosEnds.push_back(runningKmerEndPos);
                continue;
            }

            const bool prevResolvable =
                (j > 0 && actuallyResolvables.count(segmentOf(oldPath[j - 1])) == 1);
            const bool nextResolvable =
                (j + 1 < oldPath.size() &&
                 actuallyResolvables.count(segmentOf(oldPath[j + 1])) == 1);

            bool incomingExistsRc = false;
            bool incomingExistsFw = false;
            bool outgoingExistsFw = false;
            bool outgoingExistsRc = false;
            uint64_t incomingSeg = 0;
            uint64_t outgoingSeg = 0;

            if(j > 0) {
                auto it = newEdgeNodes.find({rcNode(cur), rcNode(oldPath[j - 1])});
                if(it != newEdgeNodes.end()) {
                    incomingExistsRc = true;
                    incomingSeg = it->second;
                }
                incomingExistsFw =
                    newEdgeNodes.count({oldPath[j - 1], cur}) == 1;
            }
            if(j + 1 < oldPath.size()) {
                auto it = newEdgeNodes.find({cur, oldPath[j + 1]});
                if(it != newEdgeNodes.end()) {
                    outgoingExistsFw = true;
                    outgoingSeg = it->second;
                }
                outgoingExistsRc =
                    newEdgeNodes.count({rcNode(oldPath[j + 1]), rcNode(cur)}) == 1;
            }

            if(j > 0 && j + 1 < oldPath.size() &&
               !prevResolvable && !nextResolvable &&
               !incomingExistsRc && !incomingExistsFw &&
               !outgoingExistsFw && !outgoingExistsRc) {
                maybePushBreak();
                continue;
            }
            if(j > 0 && prevResolvable &&
               !incomingExistsRc && !incomingExistsFw) {
                maybePushBreak();
                continue;
            }

            if(j > 0 && incomingExistsRc) {
                uint64_t start = runningKmerStartPos;
                if(!prevResolvable) {
                    start = (start >= overlapPrev) ? (start - overlapPrev) : 0;
                } else {
                    const uint64_t prevLen = nodeBpLength(segmentOf(oldPath[j - 1]));
                    start = (start >= prevLen) ? (start - prevLen) : 0;
                }
                newPath.push_back(revNodeId(incomingSeg));
                nodePosStarts.push_back(start);
                nodePosEnds.push_back(runningKmerEndPos);
            }

            if(j + 1 < oldPath.size() && outgoingExistsFw) {
                uint64_t start = 0;
                if(j > 0) {
                    start = runningKmerStartPos;
                    start = (start >= overlapPrev) ? (start - overlapPrev) : 0;
                }
                uint64_t end = runningKmerEndPos;
                if(j + 1 < oldPath.size() && nextResolvable) {
                    const uint64_t ovNext = getBpOverlap(cur, oldPath[j + 1]);
                    end += nodeBpLength(segmentOf(oldPath[j + 1]));
                    end = (end >= ovNext) ? (end - ovNext) : 0;
                }
                newPath.push_back(fwdNodeId(outgoingSeg));
                nodePosStarts.push_back(start);
                nodePosEnds.push_back(end);
            }
        }

        pathsCrossing.removePath(pathIdx, paths[pathIdx]);
        if(newPath.empty()) {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
            continue;
        }

        reverse(breakFromInvalidEdge.begin(), breakFromInvalidEdge.end());

        if(nodePosStarts.empty() || nodePosEnds.empty()) {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
            continue;
        }

        const uint64_t kmerPathLength = runningKmerEndPos;
        const uint64_t startRemove = nodePosStarts.front();
        const uint64_t keepEnd = nodePosEnds.back();
        if(keepEnd <= startRemove) {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
            continue;
        }

        for(size_t j = 0; j < nodePosStarts.size(); j++) {
            nodePosStarts[j] = (nodePosStarts[j] >= startRemove) ?
                (nodePosStarts[j] - startRemove) : 0;
            nodePosEnds[j] = (nodePosEnds[j] >= startRemove) ?
                (nodePosEnds[j] - startRemove) : 0;
        }

        vector<uint64_t> fullPathBoundaries;
        fullPathBoundaries.reserve(nodePosEnds.size() + 1);
        fullPathBoundaries.push_back(0);
        for(uint64_t e : nodePosEnds) {
            fullPathBoundaries.push_back(e);
        }
        const uint64_t fullPathLen = fullPathBoundaries.back();
        if(fullPathLen == 0) {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
            continue;
        }

        vector<DagPathReadSupport> rewrittenReads;
        rewrittenReads.reserve(originalSupport.size());
        for(const auto& read : originalSupport) {
            DagPathReadSupport rewritten = read;
            if(startRemove > rewritten.leftClip) {
                rewritten.readPosStartIndex += (startRemove - rewritten.leftClip);
                rewritten.leftClip = 0;
            } else {
                rewritten.leftClip -= startRemove;
            }

            uint64_t coveredRightBoundary = 0;
            if(kmerPathLength > rewritten.rightClip) {
                coveredRightBoundary = kmerPathLength - rewritten.rightClip;
            }
            if(keepEnd < coveredRightBoundary) {
                const uint64_t extraClip = coveredRightBoundary - keepEnd;
                if(extraClip >=
                   (rewritten.readPosEndIndex - rewritten.readPosStartIndex)) {
                    continue;
                }
                rewritten.readPosEndIndex -= extraClip;
                rewritten.rightClip = 0;
            } else {
                const uint64_t rightRemove = kmerPathLength - keepEnd;
                if(rewritten.rightClip < rightRemove) {
                    continue;
                }
                rewritten.rightClip -= rightRemove;
            }

            if(!reconcileReadSupportToPath(rewritten, fullPathLen)) {
                continue;
            }

            rewritten.readPoses.resize(fullPathBoundaries.size());
            for(size_t i = 0; i < fullPathBoundaries.size(); i++) {
                const uint64_t oldPos = fullPathBoundaries[i] + startRemove;
                uint64_t mapped =
                    mapPathPositionFromSupport(read, oldBoundaries, oldPos);
                mapped = min(max(mapped, rewritten.readPosStartIndex),
                             rewritten.readPosEndIndex);
                rewritten.readPoses[i] = mapped;
            }
            if(!rewritten.readPoses.empty()) {
                rewritten.readPoses.front() = rewritten.readPosStartIndex;
                rewritten.readPoses.back() = rewritten.readPosEndIndex;
                for(size_t i = 1; i < rewritten.readPoses.size(); i++) {
                    if(rewritten.readPoses[i] < rewritten.readPoses[i - 1]) {
                        rewritten.readPoses[i] = rewritten.readPoses[i - 1];
                    }
                }
            }
            rewrittenReads.push_back(std::move(rewritten));
        }

        vector<pair<uint64_t, uint64_t>> ranges;
        uint64_t lastStart = 0;
        for(uint64_t j = 1; j < newPath.size(); j++) {
            bool forcedBreak =
                (!breakFromInvalidEdge.empty() &&
                 breakFromInvalidEdge.back() == j);
            const auto& out = getOutEdges(newPath[j - 1]);
            const bool hasEdge =
                find(out.begin(), out.end(), newPath[j]) != out.end();
            if(!hasEdge || forcedBreak) {
                if(forcedBreak) {
                    breakFromInvalidEdge.pop_back();
                }
                if(j > lastStart) {
                    ranges.push_back({lastStart, j});
                }
                lastStart = j;
            }
        }
        if(lastStart < newPath.size()) {
            ranges.push_back({lastStart, uint64_t(newPath.size())});
        }

        struct Segment {
            vector<DagNodeId> nodes;
            vector<DagPathReadSupport> reads;
        };
        vector<Segment> segments;
        for(const auto& [beginIdx, endIdx] : ranges) {
            if(beginIdx >= endIdx || endIdx > newPath.size()) continue;
            Segment seg;
            seg.nodes.insert(seg.nodes.end(),
                newPath.begin() + beginIdx,
                newPath.begin() + endIdx);
            vector<uint64_t> segBoundaries;
            segBoundaries.reserve(endIdx - beginIdx + 1);
            segBoundaries.push_back(nodePosStarts[beginIdx]);
            for(uint64_t i = beginIdx; i < endIdx; i++) {
                segBoundaries.push_back(nodePosEnds[i]);
            }
            const uint64_t segLen = segBoundaries.back() - segBoundaries.front();

            seg.reads.reserve(rewrittenReads.size());
            for(const auto& read : rewrittenReads) {
                uint64_t posesStart = nodePosStarts[beginIdx];
                uint64_t posesEnd = nodePosEnds[endIdx - 1];
                uint64_t leftClipRemove = 0;
                uint64_t rightClipRemove = 0;
                const uint64_t span = read.readPosEndIndex - read.readPosStartIndex;
                const uint64_t total = read.leftClip + read.rightClip + span;

                if(posesEnd < read.leftClip + span) {
                    rightClipRemove = read.rightClip;
                    if(posesEnd <= read.leftClip) {
                        continue;
                    }
                    posesEnd -= read.leftClip;
                } else if(posesEnd < total) {
                    rightClipRemove = total - posesEnd;
                    posesEnd = span;
                }

                if(posesStart > read.leftClip) {
                    leftClipRemove = read.leftClip;
                    posesStart -= read.leftClip;
                } else if(posesStart > 0) {
                    leftClipRemove = posesStart;
                    posesStart = 0;
                }

                if(posesEnd > span + posesStart) {
                    posesEnd = span + posesStart;
                }
                if(posesStart >= posesEnd) {
                    continue;
                }

                DagPathReadSupport out = read;
                out.readPosStartIndex = read.readPosStartIndex + posesStart;
                out.readPosEndIndex = read.readPosStartIndex + posesEnd;
                if(out.leftClip < leftClipRemove ||
                   out.rightClip < rightClipRemove) {
                    continue;
                }
                out.leftClip -= leftClipRemove;
                out.rightClip -= rightClipRemove;
                if(!reconcileReadSupportToPath(out, segLen)) {
                    continue;
                }

                out.readPoses.resize(segBoundaries.size());
                for(size_t i = 0; i < segBoundaries.size(); i++) {
                    uint64_t mapped =
                        mapPathPositionFromSupport(read, fullPathBoundaries, segBoundaries[i]);
                    mapped = min(max(mapped, out.readPosStartIndex),
                                 out.readPosEndIndex);
                    out.readPoses[i] = mapped;
                }
                if(!out.readPoses.empty()) {
                    out.readPoses.front() = out.readPosStartIndex;
                    out.readPoses.back() = out.readPosEndIndex;
                    for(size_t i = 1; i < out.readPoses.size(); i++) {
                        if(out.readPoses[i] < out.readPoses[i - 1]) {
                            out.readPoses[i] = out.readPoses[i - 1];
                        }
                    }
                }

                seg.reads.push_back(std::move(out));
            }
            if(!seg.nodes.empty()) {
                segments.push_back(std::move(seg));
            }
        }

        if(segments.empty()) {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
            continue;
        }

        paths[pathIdx] = std::move(segments[0].nodes);
        pathRemoved[pathIdx] = false;
        if(pathIdx < pathReadIds.size()) {
            pathReadIds[pathIdx] = std::move(segments[0].reads);
        }
        if(pathIdx < pathWeights.size()) {
            pathWeights[pathIdx] = max<uint64_t>(
                1,
                (pathIdx < pathReadIds.size()) ? pathReadIds[pathIdx].size() : 0);
        }
        pathsCrossing.addPath(pathIdx, paths[pathIdx]);

        for(size_t i = 1; i < segments.size(); i++) {
            addPathButFirstMaybeTrim(
                std::move(segments[i].nodes),
                std::move(segments[i].reads));
        }
    }
}


// ============================================================================
// splitPathsAtBreaks
// ============================================================================

void DirectedAnchorGraph::splitPathsAtBreaks()
{
    uint64_t splitCount = 0;
    uint64_t pathsToProcess = paths.size();

    for(uint64_t pathIdx = 0; pathIdx < pathsToProcess; ++pathIdx) {
        if(pathRemoved[pathIdx]) continue;
        const vector<DagPathReadSupport> originalSupport =
            (pathIdx < pathReadIds.size()) ? pathReadIds[pathIdx] :
            vector<DagPathReadSupport>{};
        const auto& path = paths[pathIdx];
        if(path.size() < 2) continue;

        vector<uint64_t> nodeStarts(path.size(), 0);
        vector<uint64_t> nodeEnds(path.size(), 0);
        for(size_t i = 0; i < path.size(); i++) {
            if(i == 0) {
                nodeStarts[i] = 0;
            } else {
                const uint64_t ov = getBpOverlap(path[i - 1], path[i]);
                nodeStarts[i] = (nodeEnds[i - 1] >= ov) ? (nodeEnds[i - 1] - ov) : 0;
            }
            nodeEnds[i] = nodeStarts[i] + nodeBpLength(segmentOf(path[i]));
        }
        vector<uint64_t> oldBoundaries;
        oldBoundaries.reserve(path.size() + 1);
        oldBoundaries.push_back(0);
        for(uint64_t e : nodeEnds) {
            oldBoundaries.push_back(e);
        }

        // Find break positions.
        vector<uint64_t> breakPositions;
        for(uint64_t i = 0; i + 1 < path.size(); ++i) {
            if(!nodeExists(segmentOf(path[i])) ||
               !nodeExists(segmentOf(path[i + 1]))) {
                breakPositions.push_back(i);
            } else {
                const auto& outNbrs = getOutEdges(path[i]);
                if(find(outNbrs.begin(), outNbrs.end(), path[i + 1]) == outNbrs.end()) {
                    breakPositions.push_back(i);
                }
            }
        }

        if(breakPositions.empty()) continue;

        pathsCrossing.removePath(pathIdx, path);

        vector<pair<uint64_t, uint64_t>> segmentRanges;
        uint64_t segStart = 0;
        for(uint64_t bp : breakPositions) {
            if(bp + 1 > segStart) {
                segmentRanges.push_back({segStart, bp + 1});
            }
            segStart = bp + 1;
        }
        if(segStart < path.size()) {
            segmentRanges.push_back({segStart, uint64_t(path.size())});
        }

        struct Segment {
            vector<DagNodeId> nodes;
            vector<uint64_t> boundaries;
        };
        vector<Segment> segments;
        for(const auto& [beginIdx, endIdx] : segmentRanges) {
            Segment seg;
            seg.nodes.reserve(endIdx - beginIdx);
            seg.boundaries.reserve((endIdx - beginIdx) + 1);
            for(uint64_t i = beginIdx; i < endIdx; i++) {
                if(!nodeExists(segmentOf(path[i]))) {
                    continue;
                }
                if(seg.nodes.empty()) {
                    seg.boundaries.push_back(nodeStarts[i]);
                }
                seg.nodes.push_back(path[i]);
                seg.boundaries.push_back(nodeEnds[i]);
            }
            if(!seg.nodes.empty()) {
                segments.push_back(std::move(seg));
            }
        }
        if(segments.empty()) {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
            continue;
        }

        auto supportForSegment = [&](const Segment& segment) {
            vector<DagPathReadSupport> result;
            for(const auto& read : originalSupport) {
                auto projected = projectReadSupportToBoundaries(
                    read,
                    oldBoundaries,
                    segment.boundaries);
                if(projected.has_value()) {
                    result.push_back(std::move(*projected));
                }
            }
            return result;
        };

        // First segment replaces original path; rest become new paths.
        vector<DagPathReadSupport> firstSupport = supportForSegment(segments[0]);
        if(!segments.empty() && segments[0].nodes.size() >= 2) {
            paths[pathIdx] = std::move(segments[0].nodes);
            pathRemoved[pathIdx] = false;
            if(pathIdx < pathWeights.size()) {
                pathWeights[pathIdx] = max<uint64_t>(1, firstSupport.size());
            }
            if(pathIdx < pathReadIds.size()) {
                pathReadIds[pathIdx] = std::move(firstSupport);
            }
            pathsCrossing.addPath(pathIdx, paths[pathIdx]);
        } else {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
        }

        for(uint64_t s = 1; s < segments.size(); ++s) {
            if(segments[s].nodes.size() >= 2) {
                vector<DagPathReadSupport> segSupport = supportForSegment(segments[s]);
                addPathButFirstMaybeTrim(
                    std::move(segments[s].nodes),
                    std::move(segSupport));
                splitCount++;
            }
        }
    }

    if(splitCount > 0) {
        performanceLog << timestamp
            << "splitPathsAtBreaks: created " << splitCount
            << " new path segments." << endl;
    }
}

void DirectedAnchorGraph::compactGraph()
{
    vector<uint64_t> remap(nodes.size(), numeric_limits<uint64_t>::max());
    vector<DagNodeInfo> compactNodes;
    compactNodes.reserve(nodeCount());
    for(uint64_t oldSeg = 0; oldSeg < nodes.size(); oldSeg++) {
        if(!nodeExists(oldSeg)) continue;
        remap[oldSeg] = compactNodes.size();
        compactNodes.push_back(nodes[oldSeg]);
    }

    auto remapOriented = [&](DagNodeId n) {
        const uint64_t newSeg = remap[segmentOf(n)];
        return isForward(n) ? fwdNodeId(newSeg) : revNodeId(newSeg);
    };

    vector<vector<DagNodeId>> compactEdges(compactNodes.size() * 2);
    unordered_map<pair<DagNodeId, DagNodeId>, uint64_t, DagEdgePairHash> compactOverlaps;
    for(uint64_t oldSeg = 0; oldSeg < nodes.size(); oldSeg++) {
        if(remap[oldSeg] == numeric_limits<uint64_t>::max()) continue;
        for(DagNodeId side : {fwdNodeId(oldSeg), revNodeId(oldSeg)}) {
            if(side >= edges.size()) continue;
            const DagNodeId newFrom = remapOriented(side);
            for(DagNodeId to : edges[side]) {
                if(segmentOf(to) >= remap.size() ||
                   remap[segmentOf(to)] == numeric_limits<uint64_t>::max()) {
                    continue;
                }
                const DagNodeId newTo = remapOriented(to);
                auto& out = compactEdges[newFrom];
                if(find(out.begin(), out.end(), newTo) == out.end()) {
                    out.push_back(newTo);
                }
                const auto oldCanon = canonEdge(side, to);
                const auto newCanon = canonEdge(newFrom, newTo);
                auto it = edgeOverlaps.find(oldCanon);
                if(it != edgeOverlaps.end()) {
                    auto jt = compactOverlaps.find(newCanon);
                    if(jt == compactOverlaps.end()) {
                        compactOverlaps[newCanon] = it->second;
                    } else {
                        jt->second = max(jt->second, it->second);
                    }
                }
            }
        }
    }
    for(auto& out : compactEdges) {
        sort(out.begin(), out.end());
    }

    for(uint64_t pathIdx = 0; pathIdx < paths.size(); pathIdx++) {
        if(pathRemoved[pathIdx]) continue;
        auto& path = paths[pathIdx];
        vector<DagNodeId> remappedPath;
        remappedPath.reserve(path.size());
        for(DagNodeId n : path) {
            if(segmentOf(n) >= remap.size()) continue;
            const uint64_t newSeg = remap[segmentOf(n)];
            if(newSeg == numeric_limits<uint64_t>::max()) continue;
            remappedPath.push_back(isForward(n) ? fwdNodeId(newSeg) : revNodeId(newSeg));
        }
        if(remappedPath.empty()) {
            pathRemoved[pathIdx] = true;
            path.clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
        } else {
            path = std::move(remappedPath);
        }
    }

    vector<uint64_t> remappedEverTippable;
    remappedEverTippable.reserve(everTippable.size());
    for(uint64_t segId : everTippable) {
        if(segId < remap.size() &&
           remap[segId] != numeric_limits<uint64_t>::max()) {
            remappedEverTippable.push_back(remap[segId]);
        }
    }
    sort(remappedEverTippable.begin(), remappedEverTippable.end());
    remappedEverTippable.erase(
        unique(remappedEverTippable.begin(), remappedEverTippable.end()),
        remappedEverTippable.end());
    everTippable = std::move(remappedEverTippable);
    lastTippableChecked = compactNodes.size();

    nodes = std::move(compactNodes);
    edges = std::move(compactEdges);
    edgeOverlaps = std::move(compactOverlaps);
    pathsCrossing.rebuild(paths, pathRemoved);
}


// ============================================================================
// extendForward — MBG's extend(): walk forward while out-degree=1 at
// current and in-degree=1 at successor.
// ============================================================================

vector<DagNodeId> DirectedAnchorGraph::extendForward(DagNodeId start) const
{
    vector<DagNodeId> chain;
    chain.push_back(start);

    DagNodeId current = start;
    while(true) {
        if(outDegree(current) != 1) break;

        const DagNodeId next = getOutEdges(current).front();
        if(inDegree(next) != 1) break;

        // MBG circular detection: if we loop back to start, append it and stop.
        if(next == start) {
            chain.push_back(next);
            break;
        }

        // Self-loop guard (same segment, different orientation).
        if(segmentOf(next) == segmentOf(current)) break;

        current = next;
        chain.push_back(current);
    }

    return chain;
}


// ============================================================================
// replaceUnitig — merge a linear chain into a single segment.
//
// Chain traversal direction detection is simple with flat IDs:
//   Forward:  path[i] == chain[idx]          → use fwdNodeId(newSeg)
//   Reverse:  path[i] == rcNode(chain[idx])  → use revNodeId(newSeg)
// ============================================================================

uint64_t DirectedAnchorGraph::replaceUnitig(const vector<DagNodeId>& chain)
{
    DINARA_ASSERT(chain.size() >= 2);

    // Build merged node info.
    DagNodeInfo mergedInfo;
    mergedInfo.lengthBp = 0;
    double totalCoverage = 0.0;

    for(size_t i = 0; i < chain.size(); i++) {
        const DagNodeId nodeId = chain[i];
        uint64_t segId = segmentOf(nodeId);
        const DagNodeInfo& info = nodes[segId];

        if(isForward(nodeId)) {
            for(DagNodeId aid : info.anchorChain) {
                mergedInfo.anchorChain.push_back(aid);
            }
        } else {
            for(auto it = info.anchorChain.rbegin();
                it != info.anchorChain.rend(); ++it) {
                mergedInfo.anchorChain.push_back(*it ^ 1);
            }
        }

        const uint64_t curLen = nodeBpLength(segId);
        if(i == 0) {
            mergedInfo.lengthBp += curLen;
        } else {
            const uint64_t ov = getBpOverlap(chain[i - 1], chain[i]);
            mergedInfo.lengthBp += (curLen > ov) ? (curLen - ov) : 1;
        }
        totalCoverage += info.coverage;
    }
    mergedInfo.coverage = totalCoverage / double(chain.size());

    uint64_t newSeg = addNode(mergedInfo);
    const DagNodeId newFwd = fwdNodeId(newSeg);
    const DagNodeId newRev = revNodeId(newSeg);

    // Sorted vectors for O(log n) lookups — faster than hash containers
    // for typical chain sizes (2–50) due to cache locality.
    vector<uint64_t> nodesInUnitig;
    nodesInUnitig.reserve(chain.size());
    vector<pair<uint64_t, uint64_t>> chainPosVec;
    chainPosVec.reserve(chain.size());
    for(size_t i = 0; i < chain.size(); i++) {
        const uint64_t segId = segmentOf(chain[i]);
        nodesInUnitig.push_back(segId);
        chainPosVec.emplace_back(segId, i);
    }
    sort(nodesInUnitig.begin(), nodesInUnitig.end());
    sort(chainPosVec.begin(), chainPosVec.end());

    auto inUnitig = [&](uint64_t seg) -> bool {
        return binary_search(nodesInUnitig.begin(), nodesInUnitig.end(), seg);
    };
    auto findChainPos = [&](uint64_t seg) -> pair<bool, uint64_t> {
        auto it = lower_bound(chainPosVec.begin(), chainPosVec.end(),
                              pair<uint64_t, uint64_t>(seg, 0));
        if(it != chainPosVec.end() && it->first == seg) {
            return {true, it->second};
        }
        return {false, 0};
    };

    const DagNodeId first = chain.front();
    const DagNodeId last = chain.back();

    // MBG-style endpoint rewiring, including corner-case self loops.
    for(DagNodeId edge : getOutEdges(rcNode(first))) {
        const uint64_t ov = getBpOverlap(rcNode(first), edge);
        if(edge == rcNode(last)) {
            addEdge(newRev, newRev, ov);
        } else if(edge == first) {
            addEdge(newRev, newFwd, ov);
        } else if(!inUnitig(segmentOf(edge))) {
            addEdge(newRev, edge, ov);
        }
    }
    for(DagNodeId edge : getOutEdges(last)) {
        const uint64_t ov = getBpOverlap(last, edge);
        if(edge == rcNode(last)) {
            addEdge(newFwd, newRev, ov);
        } else if(!inUnitig(segmentOf(edge))) {
            addEdge(newFwd, edge, ov);
        }
    }

    vector<uint64_t> leftClip(chain.size(), 0);
    vector<uint64_t> rightClip(chain.size(), 0);
    const uint64_t firstLen = nodeBpLength(segmentOf(chain.front()));
    if(mergedInfo.lengthBp > firstLen) {
        rightClip[0] = mergedInfo.lengthBp - firstLen;
    }
    for(size_t i = 1; i < chain.size(); i++) {
        const uint64_t prevLen = nodeBpLength(segmentOf(chain[i - 1]));
        const uint64_t curLen = nodeBpLength(segmentOf(chain[i]));
        const uint64_t ov = getBpOverlap(chain[i - 1], chain[i]);
        const uint64_t addPrev = (prevLen > ov) ? (prevLen - ov) : 1;
        const uint64_t addCur = (curLen > ov) ? (curLen - ov) : 1;
        leftClip[i] = leftClip[i - 1] + addPrev;
        rightClip[i] =
            (rightClip[i - 1] >= addCur) ? (rightClip[i - 1] - addCur) : 0;
    }

    // Collect affected paths using sorted vector + dedup.
    vector<uint64_t> affectedPaths;
    for(DagNodeId nodeId : chain) {
        const auto& crossing = pathsCrossing.getPathsCrossing(segmentOf(nodeId));
        for(const auto& occ : crossing) {
            affectedPaths.push_back(occ.pathIdx);
        }
    }
    sort(affectedPaths.begin(), affectedPaths.end());
    affectedPaths.erase(
        unique(affectedPaths.begin(), affectedPaths.end()),
        affectedPaths.end());

    for(uint64_t pathIdx : affectedPaths) {
        if(pathRemoved[pathIdx]) continue;
        // Move instead of copy — avoid deep-copying path + read support vectors.
        vector<DagNodeId> oldPath = std::move(paths[pathIdx]);
        if(oldPath.empty()) continue;
        const vector<uint64_t> oldBoundaries = computePathNodeBoundaries(oldPath);
        if(oldBoundaries.size() != oldPath.size() + 1) {
            paths[pathIdx] = std::move(oldPath);  // restore on failure
            continue;
        }
        vector<DagPathReadSupport> oldSupport =
            (pathIdx < pathReadIds.size()) ? std::move(pathReadIds[pathIdx]) :
            vector<DagPathReadSupport>{};
        pathsCrossing.removePath(pathIdx, oldPath);

        vector<DagNodeId> newPath;
        vector<uint64_t> projectedBoundaries;
        newPath.reserve(oldPath.size());
        projectedBoundaries.reserve(oldPath.size() + 1);
        projectedBoundaries.push_back(0);
        size_t j = 0;
        while(j < oldPath.size()) {
            auto [found, posVal] = findChainPos(segmentOf(oldPath[j]));
            if(!found) {
                newPath.push_back(oldPath[j]);
                projectedBoundaries.push_back(oldBoundaries[j + 1]);
                j++;
                continue;
            }

            const bool isEntry =
                (j == 0) ||
                (oldPath[j] == chain.front()) ||
                (oldPath[j] == rcNode(chain.back()));
            if(!isEntry) {
                j++;
                continue;
            }

            const size_t posIdx = posVal;
            DagNodeId replacement = newFwd;
            if(j == 0) {
                bool fw = isForward(chain[posIdx]);
                if(!isForward(oldPath[j])) {
                    fw = !fw;
                }
                replacement = fw ? newFwd : newRev;
            } else if(oldPath[j] == chain.front()) {
                replacement = newFwd;
            } else if(oldPath[j] == rcNode(chain.back())) {
                replacement = newRev;
            }

            size_t k = j + 1;
            while(k < oldPath.size()) {
                if(!inUnitig(segmentOf(oldPath[k]))) break;
                if(oldPath[k] == chain.front() ||
                   oldPath[k] == rcNode(chain.back())) {
                    break;
                }
                k++;
            }

            newPath.push_back(replacement);
            projectedBoundaries.push_back(oldBoundaries[k]);
            j = (k > j) ? k : (j + 1);
        }

        if(newPath.empty() ||
           projectedBoundaries.size() != newPath.size() + 1) {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
            continue;
        }

        uint64_t extraLeftClip = 0;
        uint64_t extraRightClip = 0;
        auto [fFound, fPos] = findChainPos(segmentOf(oldPath.front()));
        if(fFound) {
            const size_t posIdx = fPos;
            bool fw = isForward(chain[posIdx]);
            if(!isForward(oldPath.front())) {
                fw = !fw;
            }
            extraLeftClip += fw ? leftClip[posIdx] : rightClip[posIdx];
        }
        auto [lFound, lPos] = findChainPos(segmentOf(oldPath.back()));
        if(lFound) {
            const size_t posIdx = lPos;
            bool fw = isForward(chain[posIdx]);
            if(!isForward(oldPath.back())) {
                fw = !fw;
            }
            extraRightClip += fw ? rightClip[posIdx] : leftClip[posIdx];
        }

        paths[pathIdx] = std::move(newPath);
        pathRemoved[pathIdx] = false;
        const uint64_t bpLen = pathBpLength(paths[pathIdx]);
        if(pathIdx < pathReadIds.size()) {
            vector<DagPathReadSupport> remappedReads;
            remappedReads.reserve(oldSupport.size());
            for(const auto& read : oldSupport) {
                auto mapped =
                    projectReadSupportToBoundaries(read, oldBoundaries, projectedBoundaries);
                if(!mapped.has_value()) {
                    continue;
                }
                DagPathReadSupport out = std::move(*mapped);
                out.leftClip += extraLeftClip;
                out.rightClip += extraRightClip;
                if(!reconcileReadSupportToPath(out, bpLen)) {
                    continue;
                }
                if(!out.readPoses.empty()) {
                    out.readPoses.front() = out.readPosStartIndex;
                    out.readPoses.back() = out.readPosEndIndex;
                    for(size_t r = 1; r < out.readPoses.size(); r++) {
                        if(out.readPoses[r] < out.readPoses[r - 1]) {
                            out.readPoses[r] = out.readPoses[r - 1];
                        }
                    }
                }
                remappedReads.push_back(std::move(out));
            }
            pathReadIds[pathIdx] = std::move(remappedReads);
        }
        if(pathIdx < pathWeights.size()) {
            const uint64_t supportCount =
                (pathIdx < pathReadIds.size()) ? pathReadIds[pathIdx].size() : 0;
            pathWeights[pathIdx] = max<uint64_t>(1, supportCount);
        }
        pathsCrossing.addPath(pathIdx, paths[pathIdx]);
    }

    // Remove old chain segments.
    for(DagNodeId nodeId : chain) {
        removeNode(segmentOf(nodeId));
    }

    return newSeg;
}


// ============================================================================
// unitigifyOne — MBG's getUnitigPath() + unitigifyOne().
// ============================================================================

uint64_t DirectedAnchorGraph::unitigifyOne(uint64_t segId)
{
    if(!nodeExists(segId)) return segId;

    const DagNodeId fwd = fwdNodeId(segId);
    const DagNodeId rev = revNodeId(segId);

    // Quick pre-check: can this node extend in at least one direction?
    // Avoids 3 vector allocations for non-mergeable nodes (~80-90% of nodes).
    const bool canFwd = (outDegree(fwd) == 1) && (inDegree(getOutEdges(fwd).front()) == 1);
    const bool canBwd = (outDegree(rev) == 1) && (inDegree(getOutEdges(rev).front()) == 1);
    if(!canFwd && !canBwd) return segId;

    // MBG: extend forward from (unitig, true) and (unitig, false).
    const vector<DagNodeId> fwdExtension = extendForward(fwd);
    const vector<DagNodeId> bwdExtension = extendForward(rev);

    // MBG: reverse the backward extension → newUnitig.
    vector<DagNodeId> newUnitig = rcPath(bwdExtension);

    if(newUnitig.size() >= 2 && newUnitig.front() == newUnitig.back()) {
        // Circular: MBG verifies fwdExtension matches, then drops duplicate.
        newUnitig.pop_back();

        // MBG: rotate so breakpoint is at the lowest segment ID.
        if(!newUnitig.empty()) {
            size_t rotationBreakpoint = 0;
            for(size_t i = 1; i < newUnitig.size(); i++) {
                if(segmentOf(newUnitig[i]) < segmentOf(newUnitig[rotationBreakpoint])) {
                    rotationBreakpoint = i;
                }
            }
            if(rotationBreakpoint != 0) {
                rotate(
                    newUnitig.begin(),
                    newUnitig.begin() + rotationBreakpoint,
                    newUnitig.end());
            }
        }
    } else {
        // Linear: MBG asserts newUnitig.back() == (unitig, true) == fwd,
        // then appends fwdExtension[1:].
        DINARA_ASSERT(newUnitig.back() == fwd);
        newUnitig.insert(
            newUnitig.end(),
            fwdExtension.begin() + 1,
            fwdExtension.end());
    }

    if(newUnitig.size() < 2) return segId;
    return replaceUnitig(newUnitig);
}


// ============================================================================
// unitigifyAll — Batch unitigification.
//
// Instead of processing one chain at a time (which rewrites the same
// paths repeatedly), this collects ALL maximal chains first, applies
// topology changes, then rewrites all paths in a single pass.
//
// Phase 1: Collect all maximal chains.
// Phase 2: Apply topology changes (merged nodes, edge rewiring, clips).
// Phase 3: Build dense lookup tables for O(1) segment → merge mapping.
// Phase 4: Single-pass path rewriting across all merges.
// Phase 5: Rebuild pathsCrossing + remove old nodes.
// ============================================================================

void DirectedAnchorGraph::unitigifyAll()
{
    performanceLog << timestamp
        << "DirectedAnchorGraph::unitigifyAll begins. "
        << nodes.size() << " total nodes." << endl;

    const uint64_t origNodeCount = nodes.size();

    // ----------------------------------------------------------------
    // Phase 1: Collect all maximal chains.
    //
    // A maximal chain is a non-overlapping sequence of nodes where each 
    // internal edge is the ONLY out-edge of its source and the ONLY 
    // in-edge of its target. We use the degree-1 pre-check to avoid 
    // costly vector allocations for nodes that clearly cannot merge.
    // ----------------------------------------------------------------

    struct ChainInfo {
        vector<DagNodeId> chain;
    };
    vector<ChainInfo> allChains;
    vector<uint8_t> claimed(origNodeCount, 0);

    for(uint64_t i = 0; i < origNodeCount; i++) {
        if(nodes[i].removed || claimed[i]) continue;

        const DagNodeId fwd = fwdNodeId(i);
        const DagNodeId rev = revNodeId(i);
        
        // Basic degree-1 check for potential chain start/internal.
        const bool canFwd = (outDegree(fwd) == 1) &&
            (inDegree(getOutEdges(fwd).front()) == 1);
        const bool canBwd = (outDegree(rev) == 1) &&
            (inDegree(getOutEdges(rev).front()) == 1);
        if(!canFwd && !canBwd) continue;

        // Extend in both directions and assemble the maximal chain.
        const vector<DagNodeId> fwdExt = extendForward(fwd);
        const vector<DagNodeId> bwdExt = extendForward(rev);
        vector<DagNodeId> chain = rcPath(bwdExt);

        if(chain.size() >= 2 && chain.front() == chain.back()) {
            // Circular unitig: simplify and rotate to a canonical starting node.
            chain.pop_back();
            if(!chain.empty()) {
                size_t bp = 0;
                for(size_t k = 1; k < chain.size(); k++) {
                    if(segmentOf(chain[k]) < segmentOf(chain[bp])) bp = k;
                }
                if(bp != 0) {
                    rotate(chain.begin(), chain.begin() + bp, chain.end());
                }
            }
        } else {
            // Linear unitig: combine the backward and forward extensions.
            DINARA_ASSERT(chain.back() == fwd);
            chain.insert(chain.end(),
                fwdExt.begin() + 1, fwdExt.end());
        }

        if(chain.size() < 2) continue;

        // Mark segments as claimed so they aren't processed twice.
        for(const DagNodeId n : chain) {
            const uint64_t seg = segmentOf(n);
            if(seg < origNodeCount) claimed[seg] = 1;
        }
        allChains.push_back({std::move(chain)});
    }

    if(allChains.empty()) {
        performanceLog << timestamp
            << "DirectedAnchorGraph::unitigifyAll ends (no chains). "
            << nodes.size() << " total nodes." << endl;
        return;
    }

    // ----------------------------------------------------------------
    // Phase 2a: Create merged nodes and compute clips.
    //
    // We create ALL merged nodes first so that we can correctly handle 
    // cross-chain edges in Phase 2c. If chain A connects to chain B, 
    // the merged node for A will need to rewire its edge to the 
    // merged node for B.
    // ----------------------------------------------------------------

    struct MergeInfo {
        vector<DagNodeId> chain;
        uint64_t newSeg;
        DagNodeId newFwd, newRev;
        vector<uint64_t> leftClip, rightClip;
    };
    vector<MergeInfo> merges;
    merges.reserve(allChains.size());

    for(auto& ci : allChains) {
        MergeInfo mi;
        mi.chain = std::move(ci.chain);

        // Build merged node info (lengths, coverage, anchor chain).
        DagNodeInfo mergedInfo;
        mergedInfo.lengthBp = 0;
        double totalCoverage = 0.0;

        for(size_t j = 0; j < mi.chain.size(); j++) {
            const DagNodeId nodeId = mi.chain[j];
            const uint64_t segId = segmentOf(nodeId);
            const DagNodeInfo& info = nodes[segId];

            if(isForward(nodeId)) {
                for(DagNodeId aid : info.anchorChain) {
                    mergedInfo.anchorChain.push_back(aid);
                }
            } else {
                for(auto it = info.anchorChain.rbegin();
                    it != info.anchorChain.rend(); ++it) {
                    mergedInfo.anchorChain.push_back(*it ^ 1);
                }
            }

            const uint64_t curLen = nodeBpLength(segId);
            if(j == 0) {
                mergedInfo.lengthBp += curLen;
            } else {
                const uint64_t ov = getBpOverlap(mi.chain[j-1], mi.chain[j]);
                mergedInfo.lengthBp += (curLen > ov) ? (curLen - ov) : 1;
            }
            totalCoverage += info.coverage;
        }
        mergedInfo.coverage = totalCoverage / double(mi.chain.size());

        mi.newSeg = addNode(mergedInfo);
        mi.newFwd = fwdNodeId(mi.newSeg);
        mi.newRev = revNodeId(mi.newSeg);

        // Compute clip offsets for each member of the chain. These clips 
        // allow us to project read positions from old nodes to the merged node.
        mi.leftClip.resize(mi.chain.size(), 0);
        mi.rightClip.resize(mi.chain.size(), 0);
        const uint64_t firstLen = nodeBpLength(segmentOf(mi.chain.front()));
        if(mergedInfo.lengthBp > firstLen) {
            mi.rightClip[0] = mergedInfo.lengthBp - firstLen;
        }
        for(size_t j = 1; j < mi.chain.size(); j++) {
            const uint64_t prevLen = nodeBpLength(segmentOf(mi.chain[j-1]));
            const uint64_t curLen  = nodeBpLength(segmentOf(mi.chain[j]));
            const uint64_t ov = getBpOverlap(mi.chain[j-1], mi.chain[j]);
            const uint64_t addPrev = (prevLen > ov) ? (prevLen - ov) : 1;
            const uint64_t addCur  = (curLen > ov)  ? (curLen - ov)  : 1;
            mi.leftClip[j]  = mi.leftClip[j-1] + addPrev;
            mi.rightClip[j] =
                (mi.rightClip[j-1] >= addCur)
                    ? (mi.rightClip[j-1] - addCur) : 0;
        }

        merges.push_back(std::move(mi));
    }

    // ----------------------------------------------------------------
    // Phase 2b: Build dense lookup tables for O(1) merge mapping.
    // ----------------------------------------------------------------

    vector<uint64_t> segMergeIdx(origNodeCount, UINT64_MAX);
    vector<uint64_t> segChainPos(origNodeCount, 0);
    for(uint64_t m = 0; m < merges.size(); m++) {
        for(uint64_t p = 0; p < merges[m].chain.size(); p++) {
            const uint64_t seg = segmentOf(merges[m].chain[p]);
            if(seg < origNodeCount) {
                segMergeIdx[seg] = m;
                segChainPos[seg] = p;
            }
        }
    }

    // Helper: map a target node to its correct merged node ID.
    // Handles cross-chain edges by resolving to the new merge.
    auto mapTarget = [&](DagNodeId edge) -> DagNodeId {
        const uint64_t seg = segmentOf(edge);
        if(seg >= origNodeCount || segMergeIdx[seg] == UINT64_MAX) {
            return edge;  // not in any chain
        }
        const auto& otherMi = merges[segMergeIdx[seg]];
        if(edge == otherMi.chain.front()) {
            return otherMi.newFwd;
        }
        if(edge == rcNode(otherMi.chain.back())) {
            return otherMi.newRev;
        }
        // Internal chain nodes cannot have external edges.
        return DagNodeId(UINT64_MAX);
    };

    // ----------------------------------------------------------------
    // Phase 2c: Edge rewiring with cross-chain mapping.
    // ----------------------------------------------------------------

    for(auto& mi : merges) {
        const DagNodeId first = mi.chain.front();
        const DagNodeId last  = mi.chain.back();

        // Rewire incoming edges from 'first' to 'newRev'.
        for(DagNodeId edge : getOutEdges(rcNode(first))) {
            const uint64_t ov = getBpOverlap(rcNode(first), edge);
            if(edge == rcNode(last)) {
                addEdge(mi.newRev, mi.newRev, ov);    // self-loop
            } else if(edge == first) {
                addEdge(mi.newRev, mi.newFwd, ov);    // circular
            } else {
                const uint64_t seg = segmentOf(edge);
                if(seg < origNodeCount && claimed[seg]) {
                    DagNodeId mapped = mapTarget(edge);
                    if(mapped != DagNodeId(UINT64_MAX)) {
                        addEdge(mi.newRev, mapped, ov);
                    }
                } else {
                    addEdge(mi.newRev, edge, ov);
                }
            }
        }

        // Rewire outgoing edges from 'last' to 'newFwd'.
        for(DagNodeId edge : getOutEdges(last)) {
            const uint64_t ov = getBpOverlap(last, edge);
            if(edge == rcNode(last)) {
                addEdge(mi.newFwd, mi.newRev, ov);    // self-loop
            } else {
                const uint64_t seg = segmentOf(edge);
                if(seg < origNodeCount && claimed[seg]) {
                    DagNodeId mapped = mapTarget(edge);
                    if(mapped != DagNodeId(UINT64_MAX)) {
                        addEdge(mi.newFwd, mapped, ov);
                    }
                } else {
                    addEdge(mi.newFwd, edge, ov);
                }
            }
        }
    }

    // ----------------------------------------------------------------
    // Phase 3: Single-pass path rewriting.
    //
    // Instead of per-merge rewriting (O(merges * paths)), we process 
    // each path exactly once. For each path, we perform all required 
    // chain substitutions in a single linear traversal.
    // ----------------------------------------------------------------

    for(uint64_t pathIdx = 0; pathIdx < paths.size(); pathIdx++) {
        if(pathRemoved[pathIdx]) continue;
        if(paths[pathIdx].empty()) continue;

        // Optimized skip: check if the path crosses ANY merge before allocating.
        bool hasAnyMerge = false;
        for(const DagNodeId n : paths[pathIdx]) {
            const uint64_t seg = segmentOf(n);
            if(seg < origNodeCount && segMergeIdx[seg] != UINT64_MAX) {
                hasAnyMerge = true;
                break;
            }
        }
        if(!hasAnyMerge) continue;

        // Move the old path and support out for reference.
        vector<DagNodeId> oldPath = std::move(paths[pathIdx]);
        const vector<uint64_t> oldBoundaries =
            computePathNodeBoundaries(oldPath);
        if(oldBoundaries.size() != oldPath.size() + 1) {
            paths[pathIdx] = std::move(oldPath); // abort and restore
            continue;
        }
        vector<DagPathReadSupport> oldSupport =
            (pathIdx < pathReadIds.size())
                ? std::move(pathReadIds[pathIdx])
                : vector<DagPathReadSupport>{};

        // Build the new path by collapsing chain segments into merged nodes.
        vector<DagNodeId> newPath;
        vector<uint64_t> projectedBoundaries;
        newPath.reserve(oldPath.size());
        projectedBoundaries.reserve(oldPath.size() + 1);
        projectedBoundaries.push_back(0);

        size_t j = 0;
        while(j < oldPath.size()) {
            const uint64_t seg = segmentOf(oldPath[j]);
            const uint64_t mIdx =
                (seg < origNodeCount) ? segMergeIdx[seg] : UINT64_MAX;

            if(mIdx == UINT64_MAX) {
                // Node isn't part of any merged chain; copy as-is.
                newPath.push_back(oldPath[j]);
                projectedBoundaries.push_back(oldBoundaries[j + 1]);
                j++;
                continue;
            }

            const MergeInfo& mi = merges[mIdx];

            // A chain can only be substituted if we enter it at an endpoint 
            // OR if it's the very start of our path.
            const bool isEntry =
                (j == 0) ||
                (oldPath[j] == mi.chain.front()) ||
                (oldPath[j] == rcNode(mi.chain.back()));
            if(!isEntry) {
                j++;
                continue;
            }

            const size_t posIdx = segChainPos[seg];
            DagNodeId replacement = mi.newFwd;
            if(j == 0) {
                // If path starts mid-chain, infer orientation.
                bool fw = isForward(mi.chain[posIdx]);
                if(!isForward(oldPath[j])) fw = !fw;
                replacement = fw ? mi.newFwd : mi.newRev;
            } else if(oldPath[j] == mi.chain.front()) {
                replacement = mi.newFwd;
            } else if(oldPath[j] == rcNode(mi.chain.back())) {
                replacement = mi.newRev;
            }

            // Skip over the consecutive sequence of nodes belonging to this chain.
            size_t k = j + 1;
            while(k < oldPath.size()) {
                const uint64_t kseg = segmentOf(oldPath[k]);
                const uint64_t km =
                    (kseg < origNodeCount) ? segMergeIdx[kseg] : UINT64_MAX;
                if(km != mIdx) break;
                // Stop if we hit a node that could also serve as a chain entry.
                if(oldPath[k] == mi.chain.front() ||
                   oldPath[k] == rcNode(mi.chain.back())) {
                    break;
                }
                k++;
            }

            newPath.push_back(replacement);
            projectedBoundaries.push_back(oldBoundaries[k]);
            j = (k > j) ? k : (j + 1);
        }

        // Validity check for the rewritten path.
        if(newPath.empty() ||
           projectedBoundaries.size() != newPath.size() + 1) {
            pathRemoved[pathIdx] = true;
            paths[pathIdx].clear();
            if(pathIdx < pathWeights.size()) pathWeights[pathIdx] = 0;
            if(pathIdx < pathReadIds.size()) pathReadIds[pathIdx].clear();
            continue;
        }

        // Project path endpoint clips to the new merged nodes.
        uint64_t extraLeftClip = 0;
        uint64_t extraRightClip = 0;
        {
            const uint64_t fs = segmentOf(oldPath.front());
            if(fs < origNodeCount && segMergeIdx[fs] != UINT64_MAX) {
                const auto& mii = merges[segMergeIdx[fs]];
                const size_t pi = segChainPos[fs];
                bool fw = isForward(mii.chain[pi]);
                if(!isForward(oldPath.front())) fw = !fw;
                extraLeftClip += fw ? mii.leftClip[pi] : mii.rightClip[pi];
            }
        }
        {
            const uint64_t bs = segmentOf(oldPath.back());
            if(bs < origNodeCount && segMergeIdx[bs] != UINT64_MAX) {
                const auto& mii = merges[segMergeIdx[bs]];
                const size_t pi = segChainPos[bs];
                bool fw = isForward(mii.chain[pi]);
                if(!isForward(oldPath.back())) fw = !fw;
                extraRightClip += fw ? mii.rightClip[pi] : mii.leftClip[pi];
            }
        }

        // Apply path and support remapping.
        paths[pathIdx] = std::move(newPath);
        pathRemoved[pathIdx] = false;
        const uint64_t bpLen = pathBpLength(paths[pathIdx]);
        if(pathIdx < pathReadIds.size()) {
            vector<DagPathReadSupport> remappedReads;
            remappedReads.reserve(oldSupport.size());
            for(const auto& read : oldSupport) {
                auto mapped = projectReadSupportToBoundaries(
                    read, oldBoundaries, projectedBoundaries);
                if(!mapped.has_value()) continue;
                DagPathReadSupport out = std::move(*mapped);
                out.leftClip += extraLeftClip;
                out.rightClip += extraRightClip;
                if(!reconcileReadSupportToPath(out, bpLen)) continue;
                // Correct read position bounds.
                if(!out.readPoses.empty()) {
                    out.readPoses.front() = out.readPosStartIndex;
                    out.readPoses.back()  = out.readPosEndIndex;
                    for(size_t r = 1; r < out.readPoses.size(); r++) {
                        if(out.readPoses[r] < out.readPoses[r-1]) {
                            out.readPoses[r] = out.readPoses[r-1];
                        }
                    }
                }
                remappedReads.push_back(std::move(out));
            }
            pathReadIds[pathIdx] = std::move(remappedReads);
        }
        if(pathIdx < pathWeights.size()) {
            const uint64_t supportCount =
                (pathIdx < pathReadIds.size())
                    ? pathReadIds[pathIdx].size() : 0;
            pathWeights[pathIdx] = max<uint64_t>(1, supportCount);
        }
    }

    // ----------------------------------------------------------------
    // Phase 4: Rebuild pathsCrossing + remove old segments.
    // ----------------------------------------------------------------

    // Full rebuild is much faster than thousands of incremental updates.
    pathsCrossing.rebuild(paths, pathRemoved);

    for(const auto& mi : merges) {
        for(const DagNodeId n : mi.chain) {
            removeNode(segmentOf(n));
        }
    }

    // ----------------------------------------------------------------
    // Phase 5: Sequential second pass for secondary merges.
    //
    // Newly created merged nodes may now form new chains with each 
    // other or with external nodes. We scan only indices ≥ origNodeCount 
    // to handle these rare (~20) secondary merges.
    // ----------------------------------------------------------------

    for(uint64_t i = origNodeCount; i < nodes.size(); i++) {
        if(nodes[i].removed) continue;
        unitigifyOne(i);
    }

    performanceLog << timestamp
        << "DirectedAnchorGraph::unitigifyAll ends. "
        << "Batch merged " << merges.size() << " chains, "
        << nodes.size() << " total nodes." << endl;
}


// ============================================================================
// Graph Cleaning Operations (MBG-style)
// ============================================================================

// ============================================================================
// MBG-style path-derived coverage helpers.
// ============================================================================

double DirectedAnchorGraph::getPathCoverage(uint64_t segId) const
{
    if(!nodeExists(segId)) return 0.0;
    double coverage = 0.0;
    const auto& crossing = pathsCrossing.getPathsCrossing(segId);
    for(const auto& occ : crossing) {
        const uint64_t weight = pathGroupWeight(occ.pathIdx);
        if(weight > 0) {
            coverage += double(weight);
        }
    }
    return coverage;
}

uint64_t DirectedAnchorGraph::getEdgePathCoverage(
    DagNodeId from,
    DagNodeId to,
    uint64_t maxCount) const
{
    if(!nodeExists(segmentOf(from)) || !nodeExists(segmentOf(to))) {
        return 0;
    }

    if(getCrossingCount(segmentOf(to)) < getCrossingCount(segmentOf(from))) {
        std::swap(from, to);
        from = rcNode(from);
        to = rcNode(to);
    }

    uint64_t result = 0;
    const auto& crossing = pathsCrossing.getPathsCrossing(segmentOf(from));
    
    // We want to avoid duplicates if multiple occurrences map to the same path index,
    // although getEdgePathCoverage logic sums "occurrences that match the edge".
    // If a path crosses 'from' multiple times, 'crossing' has multiple entries.
    // For each entry, we check if it forms the edge (from, to).
    // This naturally handles looping edges correctly (if path goes A->B -> ... -> A->B).
    
    for(const auto& occ : crossing) {
        const uint64_t pathIdx = occ.pathIdx;
        const uint32_t j = occ.offset;
        
        if(pathIdx >= paths.size() || pathRemoved[pathIdx]) continue;
        const uint64_t weight = pathGroupWeight(pathIdx);
        if(weight == 0) continue;
        
        const auto& path = paths[pathIdx];
        
        // Check forward: from -> to
        // If occ refers to 'from' node.
        if(j + 1 < path.size() && path[j] == from && path[j + 1] == to) {
            result += weight;
            if(result > maxCount) return result;
        }
        
        // Check reverse: rc(to) -> rc(from)
        // This corresponds to path[j-1] == rc(to) && path[j] == rc(from)
        // Here 'path[j]' matches 'from' if path[j] == rcNode(rcNode(from)) = from.
        // Wait, 'from' passed to this function is an oriented ID.
        // And 'path[j]' matches 'segmentOf(from)'.
        // So path[j] is either 'from' or 'rcNode(from)'.
        //
        // If path[j] == from, we check path[j+1] == to.
        // If path[j] == rcNode(from), we check path[j-1] == rcNode(to).
        //
        // NOTE: The previous loop iterated j over the whole path and checked both conditions:
        // 1. path[j] == from && path[j+1] == to
        // 2. path[j] == rcNode(from) && path[j-1] == rcNode(to)
        //
        // Now 'j' points to segmentOf(from).
        // So path[j] is either from or rcNode(from).
        
        if (path[j] == from) {
            if (j + 1 < path.size() && path[j+1] == to) {
                // Already handled above? Yes.
                // The block above `if(j+1 < ...)` covers this.
            }
        }
        
        if (path[j] == rcNode(from)) {
            if (j > 0 && path[j-1] == rcNode(to)) {
                 result += weight;
                 if(result > maxCount) return result;
            }
        }
    }
    return result;
}


// ============================================================================
// removeLowCoverageTips — MBG-like tip removal with safe-edge checks.
// ============================================================================

DirectedAnchorGraph::CleaningStats DirectedAnchorGraph::removeLowCoverageTips(
    double maxRemovableCoverage,
    double minSafeCoverage,
    uint64_t maxRemovableLength,
    const unordered_set<uint64_t>& maybeUntippable)
{
    performanceLog << timestamp
        << "removeLowCoverageTips begins (maxRemCov=" << maxRemovableCoverage
        << ", minSafeCov=" << minSafeCoverage
        << ", maxRemLen=" << maxRemovableLength << ")." << endl;

    CleaningStats stats;
    auto isPotentialTip = [&](uint64_t segId) -> bool {
        if(!nodeExists(segId)) return false;
        const DagNodeId fwd = fwdNodeId(segId);
        const DagNodeId rev = revNodeId(segId);
        if(outDegree(fwd) >= 2) return false;
        if(outDegree(rev) >= 2) return false;
        if(outDegree(fwd) == 0 && outDegree(rev) == 0) return false;
        // Use the passed thresholds for candidate identification.
        if(nodes[segId].lengthBp > maxRemovableLength) return false;
        if(getPathCoverage(segId) > maxRemovableCoverage) return false;
        return true;
    };

    auto appendTipCandidate = [&](uint64_t segId) {
        if(find(everTippable.begin(), everTippable.end(), segId) == everTippable.end()) {
            everTippable.push_back(segId);
        }
    };

    // MBG-style incremental candidate list.
    for(uint64_t segId = lastTippableChecked; segId < nodes.size(); segId++) {
        if(isPotentialTip(segId)) {
            appendTipCandidate(segId);
        }
    }
    vector<uint64_t> deterministicMaybeUntippable(
        maybeUntippable.begin(),
        maybeUntippable.end());
    sort(deterministicMaybeUntippable.begin(), deterministicMaybeUntippable.end());
    for(uint64_t segId : deterministicMaybeUntippable) {
        if(segId >= lastTippableChecked) continue;
        if(isPotentialTip(segId)) {
            appendTipCandidate(segId);
        }
    }
    lastTippableChecked = nodes.size();

    auto tryRemoveTip = [&](uint64_t segId) {
        if(!nodeExists(segId)) return;

        const DagNodeId fwd = fwdNodeId(segId);
        const DagNodeId rev = revNodeId(segId);

        bool fwdHasSafeEdge = false;
        const auto fwdEdges = getOutEdges(fwd);
        for(DagNodeId edge : fwdEdges) {
            if(getPathCoverage(segmentOf(edge)) < minSafeCoverage) return;
            if(getEdgePathCoverage(fwd, edge, uint64_t(maxRemovableCoverage)) >
               uint64_t(maxRemovableCoverage)) {
                return;
            }
            for(DagNodeId edge2 : getOutEdges(rcNode(edge))) {
                if(segmentOf(edge2) == segId) continue;
                if(getEdgePathCoverage(rcNode(edge), edge2, uint64_t(minSafeCoverage)) >=
                   uint64_t(minSafeCoverage)) {
                    fwdHasSafeEdge = true;
                    break;
                }
            }
        }
        if(!fwdEdges.empty() && !fwdHasSafeEdge) return;

        bool revHasSafeEdge = false;
        const auto revEdges = getOutEdges(rev);
        for(DagNodeId edge : revEdges) {
            if(getPathCoverage(segmentOf(edge)) < minSafeCoverage) return;
            if(getEdgePathCoverage(rev, edge, uint64_t(maxRemovableCoverage)) >
               uint64_t(maxRemovableCoverage)) {
                return;
            }
            for(DagNodeId edge2 : getOutEdges(rcNode(edge))) {
                if(segmentOf(edge2) == segId) continue;
                if(getEdgePathCoverage(rcNode(edge), edge2, uint64_t(minSafeCoverage)) >=
                   uint64_t(minSafeCoverage)) {
                    revHasSafeEdge = true;
                    break;
                }
            }
        }
        if(!revEdges.empty() && !revHasSafeEdge) return;

        for(DagNodeId edge : fwdEdges) {
            stats.maybeUnitigifiable.insert(segmentOf(edge));
        }
        for(DagNodeId edge : revEdges) {
            stats.maybeUnitigifiable.insert(segmentOf(edge));
        }

        const uint64_t edgesBefore = outDegree(fwd) + outDegree(rev);
        removeNode(segId);
        stats.nodesRemoved++;
        stats.edgesRemoved += edgesBefore;
    };

    for(size_t index = everTippable.size(); index > 0; index--) {
        const size_t idx = index - 1;
        const uint64_t segId = everTippable[idx];
        if(!nodeExists(segId)) {
            std::swap(everTippable[idx], everTippable.back());
            everTippable.pop_back();
            continue;
        }
        tryRemoveTip(segId);
        if(!nodeExists(segId)) {
            std::swap(everTippable[idx], everTippable.back());
            everTippable.pop_back();
        }
    }

    if(stats.nodesRemoved > 0) {
        splitPathsAtBreaks();
    }

    performanceLog << timestamp
        << "removeLowCoverageTips ends. Removed " << stats.nodesRemoved
        << " nodes, " << stats.edgesRemoved << " edges." << endl;

    return stats;
}


// ============================================================================
// removeLowCoverageCrosslinks
//
// MBG-style crosslink filtering. Removes weak edges that connect two
// high-coverage unitigs. This is essential for untangling the graph
// before resolution.
// ============================================================================
DirectedAnchorGraph::CleaningStats DirectedAnchorGraph::removeLowCoverageCrosslinks(
    double maxRemovableCoverage,
    uint64_t minSafeCoverageRaw)
{
    performanceLog << timestamp
        << "removeLowCoverageCrosslinks begins (maxRemCov="
        << maxRemovableCoverage << ", minSafeCov="
        << minSafeCoverageRaw << ")." << endl;

    CleaningStats stats;
    const double minSafeCoverage = double(minSafeCoverageRaw);

    struct PairHash {
        size_t operator()(const pair<DagNodeId, DagNodeId>& p) const
        {
            return hash<uint64_t>()(p.first) ^ (hash<uint64_t>()(p.second) << 1);
        }
    };
    unordered_set<pair<DagNodeId, DagNodeId>, PairHash> edgesToRemove;

    auto edgeExists = [&](DagNodeId from, DagNodeId to) {
        const auto& out = getOutEdges(from);
        return find(out.begin(), out.end(), to) != out.end();
    };

    auto tryRemoveCrosslinks = [&](DagNodeId start) {
        if(!nodeExists(segmentOf(start))) return;
        if(outDegree(start) < 2) return;

        bool possibleToRemove = false;
        for(DagNodeId edge : getOutEdges(start)) {
            if(outDegree(rcNode(edge)) >= 2) {
                possibleToRemove = true;
                break;
            }
        }
        if(!possibleToRemove) return;

        bool hasSafe = false;
        vector<DagNodeId> checkThese;
        for(DagNodeId edge : getOutEdges(start)) {
            const uint64_t coverage = getEdgePathCoverage(start, edge);
            if(double(coverage) >= minSafeCoverage) {
                hasSafe = true;
            }
            if(double(coverage) <= maxRemovableCoverage &&
               outDegree(rcNode(edge)) >= 2) {
                checkThese.push_back(edge);
            }
        }
        if(!hasSafe) return;

        for(DagNodeId edge : checkThese) {
            bool otherHasSafe = false;
            for(DagNodeId edge2 : getOutEdges(rcNode(edge))) {
                if(edge2 == rcNode(start)) continue;
                if(double(getEdgePathCoverage(rcNode(edge), edge2)) >= minSafeCoverage) {
                    otherHasSafe = true;
                    break;
                }
            }
            if(otherHasSafe) {
                edgesToRemove.insert(canonEdge(start, edge));
                stats.maybeUnitigifiable.insert(segmentOf(edge));
                stats.maybeUnitigifiable.insert(segmentOf(start));
            }
        }
    };

    for(uint64_t segId = 0; segId < nodes.size(); segId++) {
        if(!nodeExists(segId)) continue;
        tryRemoveCrosslinks(fwdNodeId(segId));
        tryRemoveCrosslinks(revNodeId(segId));
    }

    for(const auto& edge : edgesToRemove) {
        if(edgeExists(edge.first, edge.second)) {
            removeEdge(edge.first, edge.second);
            stats.edgesRemoved++;
        }
    }

    if(stats.edgesRemoved > 0) {
        splitPathsAtBreaks();
    }

    performanceLog << timestamp
        << "removeLowCoverageCrosslinks ends. Removed " << stats.edgesRemoved
        << " edges." << endl;

    return stats;
}


// ============================================================================
// cleanComponentsByCopynumber — MBG-like component copycount validation.
// ============================================================================

DirectedAnchorGraph::CleaningStats DirectedAnchorGraph::cleanComponentsByCopynumber(
    double averageCoverage)
{
    return cleanComponentsByCopynumber(
        averageCoverage,
        50000,
        0,
        0,
        {},
        0);
}


DirectedAnchorGraph::CleaningStats DirectedAnchorGraph::cleanComponentsByCopynumber(
    double averageCoverage,
    uint64_t minLongLength,
    uint64_t minUnresolvableLength,
    uint64_t maxUnresolvableLength,
    const unordered_set<uint64_t>& checkThese,
    uint64_t minNew)
{
    performanceLog << timestamp
        << "cleanComponentsByCopynumber begins (avgCov=" << averageCoverage
        << ", minLongLength=" << minLongLength
        << ", minUnresolvableLength=" << minUnresolvableLength
        << ", maxUnresolvableLength=" << maxUnresolvableLength
        << ", checkThese=" << checkThese.size()
        << ", minNew=" << minNew
        << ")." << endl;

    CleaningStats stats;
    if(averageCoverage <= 0.0) {
        averageCoverage = 1.0;
    }

    auto estimateCopyCount = [&](double coverage) -> uint64_t {
        return uint64_t((coverage + averageCoverage / 2.0) / averageCoverage);
    };

    struct ComponentResult {
        uint64_t nodesRemoved = 0;
        uint64_t edgesRemoved = 0;
        unordered_set<uint64_t> maybeUnitigifiable;
        unordered_set<DagNodeId> checked;
    };

    auto cleanComponent = [&](DagNodeId start) -> ComponentResult {
        ComponentResult result;

        vector<DagNodeId> stack;
        stack.push_back(start);
        vector<DagNodeId> componentSides;
        bool allValid = true;

        while(!stack.empty()) {
            DagNodeId top = stack.back();
            stack.pop_back();
            if(result.checked.count(top) == 1) continue;
            if(!nodeExists(segmentOf(top))) continue;
            result.checked.insert(top);

            const uint64_t segId = segmentOf(top);
            const uint64_t length = nodes[segId].lengthBp;
            if(length > minUnresolvableLength && length < maxUnresolvableLength) {
                allValid = false;
            }

            componentSides.push_back(top);
            const double coverage = getPathCoverage(segId);
            const uint64_t estimatedCopyCount = estimateCopyCount(coverage);

            if(length < minLongLength || estimatedCopyCount == 0) {
                stack.push_back(rcNode(top));
            }
            for(DagNodeId edge : getOutEdges(top)) {
                stack.push_back(rcNode(edge));
            }
        }

        if(!allValid) {
            return result;
        }

        for(DagNodeId side : componentSides) {
            const uint64_t segId = segmentOf(side);
            const double coverage = getPathCoverage(segId);
            const uint64_t estimatedCopyCount = estimateCopyCount(coverage);

            if(estimatedCopyCount <= 1 &&
               coverage < averageCoverage * (double(estimatedCopyCount) - 0.4)) {
                allValid = false;
                break;
            }
            if(estimatedCopyCount <= 1 &&
               coverage > averageCoverage * (double(estimatedCopyCount) + 0.4)) {
                allValid = false;
                break;
            }

            uint64_t edgeCopyCountSum = 0;
            for(DagNodeId edge : getOutEdges(side)) {
                const double edgeCoverage = double(getEdgePathCoverage(side, edge));
                const uint64_t edgeCopyCount = estimateCopyCount(edgeCoverage);

                if(estimatedCopyCount <= 1 &&
                   edgeCoverage < averageCoverage * (double(edgeCopyCount) - 0.4)) {
                    allValid = false;
                    break;
                }
                if(estimatedCopyCount <= 1 &&
                   edgeCoverage > averageCoverage * (double(edgeCopyCount) + 0.4)) {
                    allValid = false;
                    break;
                }
                edgeCopyCountSum += edgeCopyCount;
            }
            if(!allValid) break;

            if(edgeCopyCountSum != estimatedCopyCount) {
                allValid = false;
                break;
            }
        }

        if(!allValid) {
            return result;
        }

        struct PairHash {
            size_t operator()(const pair<DagNodeId, DagNodeId>& p) const
            {
                return hash<uint64_t>()(p.first) ^ (hash<uint64_t>()(p.second) << 1);
            }
        };
        unordered_set<pair<DagNodeId, DagNodeId>, PairHash> removedEdges;
        unordered_set<uint64_t> removedNodes;

        for(DagNodeId side : componentSides) {
            const uint64_t segId = segmentOf(side);
            for(DagNodeId edge : getOutEdges(side)) {
                const uint64_t edgeCoverage = getEdgePathCoverage(side, edge);
                const uint64_t edgeCopyCount = estimateCopyCount(double(edgeCoverage));
                if(edgeCopyCount == 0) {
                    removedEdges.insert(canonEdge(side, edge));
                }
            }
            const uint64_t estimatedCopyCount =
                estimateCopyCount(getPathCoverage(segId));
            if(estimatedCopyCount == 0) {
                removedNodes.insert(segId);
            }
        }

        auto edgeExists = [&](DagNodeId from, DagNodeId to) {
            const auto& out = getOutEdges(from);
            return find(out.begin(), out.end(), to) != out.end();
        };

        for(const auto& edge : removedEdges) {
            if(removedNodes.count(segmentOf(edge.first)) == 0) {
                result.maybeUnitigifiable.insert(segmentOf(edge.first));
            }
            if(removedNodes.count(segmentOf(edge.second)) == 0) {
                result.maybeUnitigifiable.insert(segmentOf(edge.second));
            }
        }

        for(const auto& edge : removedEdges) {
            if(edgeExists(edge.first, edge.second)) {
                removeEdge(edge.first, edge.second);
                result.edgesRemoved++;
            }
        }
        for(uint64_t segId : removedNodes) {
            if(nodeExists(segId)) {
                removeNode(segId);
                result.nodesRemoved++;
            }
        }

        return result;
    };

    auto mergePart = [&](const ComponentResult& part) {
        stats.nodesRemoved += part.nodesRemoved;
        stats.edgesRemoved += part.edgesRemoved;
        stats.maybeUnitigifiable.insert(
            part.maybeUnitigifiable.begin(),
            part.maybeUnitigifiable.end());
    };

    if(checkThese.empty() && minNew == 0) {
        vector<bool> fwChecked(nodes.size(), false);
        vector<bool> bwChecked(nodes.size(), false);

        for(uint64_t segId = 0; segId < nodes.size(); segId++) {
            if(!nodeExists(segId)) continue;

            if(!fwChecked[segId]) {
                auto part = cleanComponent(fwdNodeId(segId));
                for(DagNodeId side : part.checked) {
                    const uint64_t checkedSeg = segmentOf(side);
                    if(checkedSeg >= fwChecked.size()) continue;
                    if(isForward(side)) fwChecked[checkedSeg] = true;
                    else bwChecked[checkedSeg] = true;
                }
                mergePart(part);
            }

            if(!nodeExists(segId)) continue;

            if(!bwChecked[segId]) {
                auto part = cleanComponent(revNodeId(segId));
                for(DagNodeId side : part.checked) {
                    const uint64_t checkedSeg = segmentOf(side);
                    if(checkedSeg >= fwChecked.size()) continue;
                    if(isForward(side)) fwChecked[checkedSeg] = true;
                    else bwChecked[checkedSeg] = true;
                }
                mergePart(part);
            }
        }
    } else {
        unordered_set<DagNodeId> checked;
        vector<uint64_t> deterministicCheckThese(checkThese.begin(), checkThese.end());
        sort(deterministicCheckThese.begin(), deterministicCheckThese.end());

        auto runTargeted = [&](uint64_t segId, bool forward) {
            if(!nodeExists(segId)) return;
            const DagNodeId side = forward ? fwdNodeId(segId) : revNodeId(segId);
            if(checked.count(side) == 1) return;
            auto part = cleanComponent(side);
            for(DagNodeId nodeSide : part.checked) {
                const uint64_t nodeSeg = segmentOf(nodeSide);
                if(checkThese.count(nodeSeg) == 1 || nodeSeg >= minNew) {
                    checked.insert(nodeSide);
                }
            }
            mergePart(part);
        };

        for(uint64_t segId : deterministicCheckThese) {
            runTargeted(segId, true);
            runTargeted(segId, false);
        }
        for(uint64_t segId = minNew; segId < nodes.size(); segId++) {
            runTargeted(segId, true);
            runTargeted(segId, false);
        }
    }

    if(stats.nodesRemoved > 0 || stats.edgesRemoved > 0) {
        splitPathsAtBreaks();
    }

    performanceLog << timestamp
        << "cleanComponentsByCopynumber ends. "
        << "Removed " << stats.nodesRemoved << " nodes, "
        << stats.edgesRemoved << " edges." << endl;

    return stats;
}


// ============================================================================
// MBG-style parity helpers used by resolveRound.
// ============================================================================

bool DirectedAnchorGraph::nodeIsPalindrome(uint64_t segId) const
{
    if(!nodeExists(segId)) return false;
    unordered_set<uint64_t> readsWhichCoverFw;
    unordered_set<uint64_t> readsWhichCoverBw;
    const auto& crossing = pathsCrossing.getPathsCrossing(segId);
    
    for(const auto& occ : crossing) {
        const uint64_t pathIdx = occ.pathIdx;
        const size_t pos = occ.offset;
        
        if(pathIdx >= paths.size() || pathRemoved[pathIdx] ||
           pathIdx >= pathReadIds.size()) {
            continue;
        }
        const auto& path = paths[pathIdx];
        if(path.size() < 2 || pathReadIds[pathIdx].empty()) continue;
        
        // Use occ.offset logic
        if(pos > 0 && pos + 1 < path.size()) {
            // Verify? path[pos] should be segmentOf segId.
            // if(segmentOf(path[pos]) != segId) continue; // Should be guaranteed by index.
            
            for(const auto& read : pathReadIds[pathIdx]) {
                const uint64_t readName = read.readNameIndex;
                if(isForward(path[pos])) {
                    if(readsWhichCoverBw.count(readName) == 1) return true;
                    readsWhichCoverFw.insert(readName);
                } else {
                    if(readsWhichCoverFw.count(readName) == 1) return true;
                    readsWhichCoverBw.insert(readName);
                }
            }
        }
    }
    return false;
}

bool DirectedAnchorGraph::isLocallyRepetitive(uint64_t segId) const
{
    if(!nodeExists(segId)) return false;
    const auto& crossing = pathsCrossing.getPathsCrossing(segId);

    // MBG: direct repeated occurrences on same path.
    // Since crossing now contains all occurrences, if pathIdx appears >= 2 times,
    // it means the path crosses the segment multiple times.
    // Optimization: Check for duplicate pathIdx in crossing.
    {
        // Vector to sort and check duplicates.
        // We only care about active paths.
        static thread_local vector<uint64_t> activePaths;
        activePaths.clear();
        activePaths.reserve(crossing.size());
        
        for(const auto& occ : crossing) {
             if(occ.pathIdx < paths.size() && !pathRemoved[occ.pathIdx]) {
                 activePaths.push_back(occ.pathIdx);
             }
        }
        sort(activePaths.begin(), activePaths.end());
        for(size_t i = 1; i < activePaths.size(); ++i) {
            if(activePaths[i] == activePaths[i-1]) return true;
        }
    }

    // MBG: repeated read support over internal occurrences.
    unordered_set<uint64_t> readsWhichCover;
    for(const auto& occ : crossing) {
        const uint64_t pathIdx = occ.pathIdx;
        const size_t pos = occ.offset;

        if(pathIdx >= paths.size() || pathRemoved[pathIdx] ||
           pathIdx >= pathReadIds.size()) {
            continue;
        }
        const auto& path = paths[pathIdx];
        if(path.size() < 2 || pathReadIds[pathIdx].empty()) continue;
        
        if(pos > 0 && pos + 1 < path.size()) {
             for(const auto& read : pathReadIds[pathIdx]) {
                const uint64_t readName = read.readNameIndex;
                if(readsWhichCover.count(readName) == 1) return true;
                readsWhichCover.insert(readName);
            }
        }
    }
    return false;
}

unordered_set<uint64_t> DirectedAnchorGraph::filterToOnlyLocallyRepetitives(
    const unordered_set<uint64_t>& unfilteredResolvables,
    uint64_t maxDist,
    bool resolvePalindromesGlobal) const
{
    unordered_set<uint64_t> result;
    for(uint64_t segId : unfilteredResolvables) {
        if(!nodeExists(segId)) continue;
        if(resolvePalindromesGlobal && nodeIsPalindrome(segId)) {
            result.insert(segId);
            continue;
        }
        if(maxDist == numeric_limits<uint64_t>::max()) {
            if(isLocallyRepetitive(segId)) {
                result.insert(segId);
            }
            continue;
        }

        bool found = false;
        const auto& crossing = pathsCrossing.getPathsCrossing(segId);

        // Filter to unique paths to avoid O(N^2) behavior on paths with many repeats.
        static thread_local vector<uint64_t> uniquePaths;
        uniquePaths.clear();
        uniquePaths.reserve(crossing.size());
        for(const auto& occ : crossing) {
            if(occ.pathIdx < paths.size() && !pathRemoved[occ.pathIdx]) {
                uniquePaths.push_back(occ.pathIdx);
            }
        }
        sort(uniquePaths.begin(), uniquePaths.end());
        uniquePaths.erase(unique(uniquePaths.begin(), uniquePaths.end()), uniquePaths.end());

        for(uint64_t pathIdx : uniquePaths) {
            const auto& path = paths[pathIdx];
            if(path.size() < 2) continue;
            uint64_t last = numeric_limits<uint64_t>::max();
            for(uint64_t i = 0; i < path.size(); i++) {
                if(segmentOf(path[i]) != segId) continue;
                if(last == numeric_limits<uint64_t>::max()) {
                    last = i;
                    continue;
                }
                uint64_t dist = 0;
                for(uint64_t j = last; j <= i; j++) {
                    if(j < i) {
                        dist += nodeBpLength(segmentOf(path[j]));
                    }
                    if(j > last) {
                        const uint64_t ov = getBpOverlap(path[j - 1], path[j]);
                        if(dist >= ov) dist -= ov;
                        else dist = 0;
                    }
                    if(dist >= maxDist) break;
                }
                if(dist < maxDist) {
                    result.insert(segId);
                    found = true;
                    break;
                }
                last = i;
            }
            if(found) break;
        }
    }
    return result;
}

unordered_set<uint64_t> DirectedAnchorGraph::trimNodes(
    const unordered_map<DagNodeId, uint64_t>& maybeTrimmable)
{
    unordered_set<uint64_t> maybeUnitigifiable;

    auto maxReadTrimAtBoundary = [&](DagNodeId pos) -> uint64_t {
        const uint64_t segId = segmentOf(pos);
        if(!nodeExists(segId)) return 0;
        uint64_t maxReadTrim = nodeBpLength(segId);
        bool hasBoundary = false;
        const auto crossing = pathsCrossing.getPathsCrossing(segId);
        for(const auto& occ : crossing) {
            const uint64_t pathIdx = occ.pathIdx;
            const size_t j = occ.offset;

            if(pathIdx >= paths.size() || pathRemoved[pathIdx]) continue;
            if(pathIdx >= pathReadIds.size()) continue;
            const auto& path = paths[pathIdx];
            
            // Check boundary condition using offset j.
            if(j == 0 && path[j] == rcNode(pos)) {
                hasBoundary = true;
                for(const auto& read : pathReadIds[pathIdx]) {
                    maxReadTrim = min(maxReadTrim, read.leftClip);
                }
                if(maxReadTrim == 0) return 0;
            } else if(j + 1 == path.size() && path[j] == pos) {
                hasBoundary = true;
                for(const auto& read : pathReadIds[pathIdx]) {
                    maxReadTrim = min(maxReadTrim, read.rightClip);
                }
                if(maxReadTrim == 0) return 0;
            }
        }
        return hasBoundary ? maxReadTrim : 0;
    };

    function<bool(DagNodeId, uint64_t, unordered_set<DagNodeId>&)> canTrimRecursive;
    canTrimRecursive = [&](DagNodeId pos,
                           uint64_t trimAmount,
                           unordered_set<DagNodeId>& visiting) -> bool {
        const uint64_t segId = segmentOf(pos);
        if(!nodeExists(segId)) return false;
        const uint64_t nodeLen = nodeBpLength(segId);
        if(trimAmount == 0 || trimAmount >= nodeLen) return false;
        if(visiting.count(pos) == 1) return false;
        visiting.insert(pos);

        for(DagNodeId edge : getOutEdges(pos)) {
            if(getBpOverlap(pos, edge) < trimAmount) {
                visiting.erase(pos);
                return false;
            }
        }
        for(DagNodeId edge : getOutEdges(rcNode(pos))) {
            const uint64_t overlap = getBpOverlap(rcNode(pos), edge);
            if(overlap >= nodeLen - trimAmount) {
                if(overlap >= nodeLen) {
                    visiting.erase(pos);
                    return false;
                }
                const uint64_t trimThere = overlap - (nodeLen - trimAmount) + 1;
                if(trimThere == 0 ||
                   trimThere >= nodeBpLength(segmentOf(edge))) {
                    visiting.erase(pos);
                    return false;
                }
                if(!canTrimRecursive(rcNode(edge), trimThere, visiting)) {
                    visiting.erase(pos);
                    return false;
                }
            }
        }

        const uint64_t maxReadTrim = maxReadTrimAtBoundary(pos);
        visiting.erase(pos);
        return maxReadTrim >= trimAmount;
    };

    function<void(DagNodeId, uint64_t, unordered_set<DagNodeId>&)> trimEndRecursive;
    trimEndRecursive = [&](DagNodeId pos,
                           uint64_t trimAmount,
                           unordered_set<DagNodeId>& visiting) {
        const uint64_t segId = segmentOf(pos);
        if(!nodeExists(segId)) return;
        if(visiting.count(pos) == 1) return;
        visiting.insert(pos);

        const uint64_t nodeLen = nodeBpLength(segId);
        vector<DagNodeId> outPos = getOutEdges(pos);
        for(DagNodeId edge : outPos) {
            const uint64_t ov = getBpOverlap(pos, edge);
            setBpOverlap(pos, edge, (ov >= trimAmount) ? (ov - trimAmount) : 0);
        }

        vector<pair<DagNodeId, uint64_t>> recurse;
        vector<DagNodeId> outRev = getOutEdges(rcNode(pos));
        for(DagNodeId edge : outRev) {
            const uint64_t overlap = getBpOverlap(rcNode(pos), edge);
            if(overlap >= nodeLen - trimAmount && overlap < nodeLen) {
                const uint64_t trimThere = overlap - (nodeLen - trimAmount) + 1;
                if(trimThere > 0 &&
                   trimThere < nodeBpLength(segmentOf(edge))) {
                    recurse.push_back({rcNode(edge), trimThere});
                }
            }
        }

        if(nodes[segId].lengthBp > trimAmount) {
            nodes[segId].lengthBp -= trimAmount;
            maybeUnitigifiable.insert(segId);
        } else {
            removeNode(segId);
            visiting.erase(pos);
            return;
        }

        const auto crossing = pathsCrossing.getPathsCrossing(segId);
        for(const auto& occ : crossing) {
            const uint64_t pathIdx = occ.pathIdx;
            const size_t j = occ.offset;

            if(pathIdx >= paths.size() || pathRemoved[pathIdx]) continue;
            if(pathIdx >= pathReadIds.size()) continue;
            const auto& path = paths[pathIdx];
            
            if(j == 0 && path[j] == rcNode(pos)) {
                for(auto& read : pathReadIds[pathIdx]) {
                    if(read.leftClip < trimAmount) {
                        read.leftClip = 0;
                    } else {
                        read.leftClip -= trimAmount;
                    }
                }
            } else if(j + 1 == path.size() && path[j] == pos) {
                for(auto& read : pathReadIds[pathIdx]) {
                    if(read.rightClip < trimAmount) {
                        read.rightClip = 0;
                    } else {
                        read.rightClip -= trimAmount;
                    }
                }
            }
        }

        for(const auto& p : recurse) {
            trimEndRecursive(p.first, p.second, visiting);
        }
        visiting.erase(pos);
    };

    function<void(DagNodeId, uint64_t)> maybeTrim;
    maybeTrim = [&](DagNodeId pos, uint64_t maxTrim) {
        const uint64_t segId = segmentOf(pos);
        if(!nodeExists(segId)) return;
        if(maxTrim == 0) return;
        if(outDegree(pos) > 0) return;

        uint64_t maxReadTrim = maxReadTrimAtBoundary(pos);
        if(maxReadTrim == 0 || maxReadTrim > maxTrim) return;
        if(maxReadTrim >= nodeBpLength(segId)) {
            return;
        }

        unordered_set<DagNodeId> visiting;
        if(!canTrimRecursive(pos, maxReadTrim, visiting)) {
            vector<pair<DagNodeId, uint64_t>> alsoTrim;
            for(DagNodeId edge : getOutEdges(rcNode(pos))) {
                maybeUnitigifiable.insert(segmentOf(edge));
                const uint64_t overlap = getBpOverlap(rcNode(pos), edge);
                if(overlap > 0) {
                    alsoTrim.push_back({rcNode(edge), overlap});
                }
            }
            removeNode(segId);
            for(const auto& p : alsoTrim) {
                if(nodeExists(segmentOf(p.first))) {
                    maybeTrim(p.first, p.second);
                }
            }
            return;
        }

        unordered_set<DagNodeId> trimVisiting;
        trimEndRecursive(pos, maxReadTrim, trimVisiting);
    };

    vector<pair<DagNodeId, uint64_t>> order;
    order.reserve(maybeTrimmable.size());
    for(const auto& entry : maybeTrimmable) {
        order.push_back(entry);
    }
    sort(order.begin(), order.end(),
         [](const auto& left, const auto& right) {
             if(left.first < right.first) return true;
             if(left.first > right.first) return false;
             return left.second > right.second;
         });

    for(const auto& p : order) {
        if(!nodeExists(segmentOf(p.first))) continue;
        maybeTrim(p.first, p.second);
    }
    return maybeUnitigifiable;
}

void DirectedAnchorGraph::addPathButFirstMaybeTrim(
    vector<DagNodeId>&& path,
    vector<DagPathReadSupport>&& reads)
{
    if(path.empty() || reads.empty()) return;
    if(path.size() == 1) {
        addPath(std::move(path), reads.size(), std::move(reads));
        return;
    }

    const vector<uint64_t> oldBoundaries = computePathNodeBoundaries(path);
    if(oldBoundaries.size() < 3) {
        addPath(std::move(path), reads.size(), std::move(reads));
        return;
    }
    const uint64_t leftCompare = oldBoundaries[1] - oldBoundaries[0];
    const uint64_t rightCompare =
        oldBoundaries.back() - oldBoundaries[oldBoundaries.size() - 2];

    vector<DagPathReadSupport> baseReads;
    vector<DagPathReadSupport> leftTrimReads;
    vector<DagPathReadSupport> rightTrimReads;
    vector<DagPathReadSupport> doubleTrimReads;
    baseReads.reserve(reads.size());
    for(const auto& read : reads) {
        if(path.size() >= 3 &&
           read.leftClip >= leftCompare &&
           read.rightClip >= rightCompare) {
            doubleTrimReads.push_back(read);
        } else if(read.leftClip >= leftCompare) {
            leftTrimReads.push_back(read);
        } else if(read.rightClip >= rightCompare) {
            rightTrimReads.push_back(read);
        } else {
            baseReads.push_back(read);
        }
    }

    auto addProjectedPath = [&](
        size_t beginNode,
        size_t endNode,
        vector<DagPathReadSupport>& groupReads) {
        if(groupReads.empty() || beginNode >= endNode || endNode > path.size()) return;
        vector<DagNodeId> subPath(path.begin() + beginNode, path.begin() + endNode);
        vector<uint64_t> subBoundaries(
            oldBoundaries.begin() + beginNode,
            oldBoundaries.begin() + endNode + 1);
        vector<DagPathReadSupport> projected;
        projected.reserve(groupReads.size());
        for(const auto& read : groupReads) {
            auto mapped = projectReadSupportToBoundaries(read, oldBoundaries, subBoundaries);
            if(mapped.has_value()) {
                projected.push_back(std::move(*mapped));
            }
        }
        if(!projected.empty()) {
            addPath(std::move(subPath), projected.size(), std::move(projected));
        }
    };

    addProjectedPath(0, path.size(), baseReads);
    addProjectedPath(1, path.size(), leftTrimReads);
    addProjectedPath(0, path.size() - 1, rightTrimReads);
    if(path.size() >= 3) {
        addProjectedPath(1, path.size() - 1, doubleTrimReads);
    }
}

void DirectedAnchorGraph::checkValidity() const
{
    return;
    for(uint64_t segId = 0; segId < nodes.size(); segId++) {
        if(nodes[segId].removed) {
            if(fwdNodeId(segId) < edges.size() && !edges[fwdNodeId(segId)].empty()) {
                throw runtime_error("Validity error: removed node has forward edges.");
            }
            if(revNodeId(segId) < edges.size() && !edges[revNodeId(segId)].empty()) {
                throw runtime_error("Validity error: removed node has reverse edges.");
            }
            continue;
        }
        for(DagNodeId side : {fwdNodeId(segId), revNodeId(segId)}) {
            if(side >= edges.size()) continue;
            for(DagNodeId to : edges[side]) {
                if(!nodeExists(segmentOf(to))) {
                    throw runtime_error("Validity error: edge points to removed node.");
                }
                const auto& rcOut = getOutEdges(rcNode(to));
                if(find(rcOut.begin(), rcOut.end(), rcNode(side)) == rcOut.end()) {
                    throw runtime_error("Validity error: missing reverse-complement edge.");
                }
                const uint64_t ov = getBpOverlap(side, to);
                if(ov >= nodeBpLength(segmentOf(side)) ||
                   ov >= nodeBpLength(segmentOf(to))) {
                    throw runtime_error("Validity error: overlap >= node length.");
                }
            }
        }
    }

    for(uint64_t pathIdx = 0; pathIdx < paths.size(); pathIdx++) {
        if(pathRemoved[pathIdx]) continue;
        if(paths[pathIdx].empty()) {
            throw runtime_error("Validity error: active path is empty.");
        }
        for(size_t i = 1; i < paths[pathIdx].size(); i++) {
            const auto& out = getOutEdges(paths[pathIdx][i - 1]);
            if(find(out.begin(), out.end(), paths[pathIdx][i]) == out.end()) {
                throw runtime_error("Validity error: path uses missing edge.");
            }
        }
        if(pathIdx >= pathReadIds.size()) {
            throw runtime_error("Validity error: missing path read supports.");
        }
        const vector<uint64_t> boundaries = computePathNodeBoundaries(paths[pathIdx]);
        const uint64_t bpLen = boundaries.empty() ? 0 : boundaries.back();
        for(const auto& read : pathReadIds[pathIdx]) {
            if(!isReadSupportStrictForPath(read, boundaries, bpLen)) {
                throw runtime_error("Validity error: inconsistent path read support.");
            }
        }
    }
}


// ============================================================================
// resolveRound
// ============================================================================

void DirectedAnchorGraph::resolveRound(
    uint64_t minEdgeSupport,
    uint64_t maxResolveLength,
    bool doCleaning,
    bool doGuesswork,
    uint64_t maxUnconditionalResolveLength,
    bool copycountFilterHeuristic,
    uint64_t maxLocalResolve,
    bool resolvePalindromesGlobal)
{
    performanceLog << timestamp
        << "DirectedAnchorGraph::resolveRound begins (minEdgeSupport="
        << minEdgeSupport << ", maxResolveLength=" << maxResolveLength
        << ", maxUnconditionalResolveLength=" << maxUnconditionalResolveLength
        << ", maxLocalResolve=" << maxLocalResolve
        << ", resolvePalindromesGlobal=" << resolvePalindromesGlobal
        << ", copycountFilterHeuristic=" << copycountFilterHeuristic
        << ")." << endl;
    cout << timestamp << "  resolveRound(minEdgeSupport=" << minEdgeSupport << ") begins." << endl;
    checkValidity();

    // Min-heap of (lengthBp, segId).
    using HeapEntry = pair<uint64_t, uint64_t>;
    priority_queue<HeapEntry, vector<HeapEntry>, greater<HeapEntry>> heap;

    for(uint64_t segId = 0; segId < nodes.size(); ++segId) {
        const auto& info = nodes[segId];
        if(info.removed) continue;
        heap.push({info.lengthBp, segId});
    }

    uint64_t resolvedNodes = 0;
    uint64_t removedSinceCompact = 0;

    auto addPlusOneComponent = [&](
        unordered_set<uint64_t>& resolvables,
        uint64_t startSeg,
        uint64_t maxLength) {
        if(!nodeExists(startSeg)) return;
        unordered_set<uint64_t> checked;
        vector<uint64_t> stack;
        stack.push_back(startSeg);

        while(!stack.empty()) {
            const uint64_t seg = stack.back();
            stack.pop_back();
            if(checked.count(seg) == 1) continue;
            if(!nodeExists(seg)) continue;
            checked.insert(seg);

            auto visitSide = [&](DagNodeId side) -> bool {
                for(DagNodeId edge : getOutEdges(side)) {
                    const uint64_t nextSeg = segmentOf(edge);
                    if(!nodeExists(nextSeg)) continue;
                    const uint64_t nextLen = nodeBpLength(nextSeg);
                    if(nextLen == getBpOverlap(side, edge) + 1) {
                        if(nextLen > maxLength) {
                            checked.clear();
                            return false;
                        }
                        stack.push_back(nextSeg);
                    }
                }
                return true;
            };
            if(!visitSide(fwdNodeId(seg))) {
                return;
            }
            if(!visitSide(revNodeId(seg))) {
                return;
            }
        }
        if(checked.size() > 1) {
            resolvables.insert(checked.begin(), checked.end());
        }
    };

    auto shouldAddResolvable = [&](uint64_t segId, uint64_t topSize) -> bool {
        if(!nodeExists(segId)) return false;
        const DagNodeId fwd = fwdNodeId(segId);
        const DagNodeId rev = revNodeId(segId);
        const uint64_t fwdOut = outDegree(fwd);
        const uint64_t revOut = outDegree(rev);

        if(fwdOut >= 2 || revOut >= 2) {
            return true;
        }
        if(fwdOut >= 1 && revOut >= 1) {
            if(fwdOut == 1 && revOut == 1) {
                const DagNodeId fwdNbr = getOutEdges(fwd).front();
                const DagNodeId revNbr = getOutEdges(rev).front();
                if(segmentOf(fwdNbr) != segId &&
                   (nodes[segmentOf(fwdNbr)].lengthBp == topSize ||
                    nodes[segmentOf(revNbr)].lengthBp == topSize)) {
                    return true;
                }
            }
        }
        return false;
    };

    auto computeAveragePathCoverage = [&]() -> double {
        double totalCov = 0.0;
        uint64_t count = 0;
        for(uint64_t segId = 0; segId < nodes.size(); ++segId) {
            if(!nodeExists(segId)) continue;
            totalCov += getPathCoverage(segId);
            count++;
        }
        return (count > 0) ? (totalCov / double(count)) : 1.0;
    };

    auto unitigifyCandidates = [&](const unordered_set<uint64_t>& maybeUnitigifiable) {
        vector<uint64_t> ordered(maybeUnitigifiable.begin(), maybeUnitigifiable.end());
        sort(ordered.begin(), ordered.end());
        for(uint64_t segId : ordered) {
            if(!nodeExists(segId)) continue;
            const uint64_t merged = unitigifyOne(segId);
            if(nodeExists(merged)) {
                heap.push({nodes[merged].lengthBp, merged});
            }
        }
    };

    while(!heap.empty()) {
        const uint64_t currentLength = heap.top().first;
        if(currentLength >= maxResolveLength &&
           currentLength >= maxLocalResolve) {
            break;
        }
        unordered_set<uint64_t> resolvables;
        unordered_set<uint64_t> thisLengthNodes;

        while(!heap.empty() && heap.top().first == currentLength) {
            const uint64_t segId = heap.top().second;
            heap.pop();
            if(nodeExists(segId)) {
                thisLengthNodes.insert(segId);
                addPlusOneComponent(resolvables, segId, currentLength);
                if(shouldAddResolvable(segId, currentLength)) {
                    resolvables.insert(segId);
                }
            }
        }

        if(resolvables.empty()) continue;

        if(maxLocalResolve > 0) {
            resolvables = filterToOnlyLocallyRepetitives(
                resolvables,
                maxLocalResolve,
                resolvePalindromesGlobal);
            unordered_set<uint64_t> oldResolvables = resolvables;
            for(uint64_t segId : oldResolvables) {
                addPlusOneComponent(resolvables, segId, currentLength);
            }
        }
        if(resolvables.empty()) continue;

        const uint64_t activeBeforeBatch = nodeCount();
        const uint64_t oldSize = nodes.size();

        vector<uint64_t> batch(resolvables.begin(), resolvables.end());
        sort(batch.begin(), batch.end());

        bool changedThisBatch = false;

        vector<uint64_t> remaining;
        remaining.reserve(batch.size());
        for(uint64_t segId : batch) {
            if(nodeExists(segId)) {
                remaining.push_back(segId);
            }
        }
        if(remaining.empty()) continue;

        const bool unconditional =
            (maxUnconditionalResolveLength > 0 &&
             currentLength < maxUnconditionalResolveLength);
        const double avgPathCov =
            (copycountFilterHeuristic && doGuesswork) ?
            computeAveragePathCoverage() : 1.0;

        auto resolveStats = resolveNodes(
            remaining,
            minEdgeSupport,
            unconditional,
            copycountFilterHeuristic,
            doGuesswork,
            avgPathCov);
        if(resolveStats.nodesResolved > 0) {
            changedThisBatch = true;
            resolvedNodes += resolveStats.nodesResolved;
        }

        if(!changedThisBatch) continue;

        unordered_set<uint64_t> maybeUnitigifiable = resolveStats.maybeUnitigifiable;

        if(!resolveStats.maybeTrimmable.empty()) {
            auto trimmed = trimNodes(resolveStats.maybeTrimmable);
            maybeUnitigifiable.insert(trimmed.begin(), trimmed.end());
        }
        unitigifyCandidates(maybeUnitigifiable);

        if(doCleaning) {
            if(doGuesswork) {
                unordered_set<uint64_t> cleanables = thisLengthNodes;
                cleanables.insert(
                    resolveStats.maybeUnitigifiable.begin(),
                    resolveStats.maybeUnitigifiable.end());
                const double avgCov = computeAveragePathCoverage();
                auto copyStats = cleanComponentsByCopynumber(
                    avgCov,
                    50000,
                    currentLength,
                    0,
                    cleanables,
                    oldSize);
                if(copyStats.nodesRemoved > 0 || copyStats.edgesRemoved > 0) {
                    unitigifyCandidates(copyStats.maybeUnitigifiable);
                }
            }

            auto tipStats1 = removeLowCoverageTips(3.0, 10.0, 10000, resolveStats.maybeUnitigifiable);
            resolveStats.maybeUnitigifiable.insert(
                tipStats1.maybeUnitigifiable.begin(),
                tipStats1.maybeUnitigifiable.end());
            auto tipStats2 = removeLowCoverageTips(
                2.0, 5.0, 10000, resolveStats.maybeUnitigifiable);
            auto crossStats1 = removeLowCoverageCrosslinks(1.0, 5);
            auto crossStats2 = removeLowCoverageCrosslinks(2.0, 10);

            if(tipStats1.nodesRemoved > 0 || tipStats2.nodesRemoved > 0 ||
               tipStats1.edgesRemoved > 0 || tipStats2.edgesRemoved > 0 ||
               crossStats1.edgesRemoved > 0 || crossStats2.edgesRemoved > 0) {
                unordered_set<uint64_t> cleanMaybe = tipStats1.maybeUnitigifiable;
                cleanMaybe.insert(
                    tipStats2.maybeUnitigifiable.begin(),
                    tipStats2.maybeUnitigifiable.end());
                cleanMaybe.insert(
                    crossStats1.maybeUnitigifiable.begin(),
                    crossStats1.maybeUnitigifiable.end());
                cleanMaybe.insert(
                    crossStats2.maybeUnitigifiable.begin(),
                    crossStats2.maybeUnitigifiable.end());
                unitigifyCandidates(cleanMaybe);
            }
        }

        const uint64_t activeAfterBatch = nodeCount();
        if(activeBeforeBatch > activeAfterBatch) {
            removedSinceCompact += (activeBeforeBatch - activeAfterBatch);
        }

        bool didCompact = false;
        if(totalNodeCount() > 0 &&
           removedSinceCompact > totalNodeCount() / 2 + 50) {
            compactGraph();
            removedSinceCompact = 0;
            didCompact = true;
        }

        if(didCompact) {
            while(!heap.empty()) {
                heap.pop();
            }
            for(uint64_t segId = 0; segId < nodes.size(); segId++) {
                if(!nodeExists(segId)) continue;
                heap.push({nodes[segId].lengthBp, segId});
            }
            continue;
        }

        for(uint64_t segId = oldSize; segId < nodes.size(); segId++) {
            if(!nodeExists(segId)) continue;
            heap.push({nodes[segId].lengthBp, segId});
        }
    }

    performanceLog << timestamp
        << "DirectedAnchorGraph::resolveRound ends. "
        << "Resolved " << resolvedNodes << " nodes, "
        << nodeCount() << " nodes remaining." << endl;
}


// ============================================================================
// runResolution — full pipeline.
// ============================================================================

void DirectedAnchorGraph::runResolution(const DirectedAnchorGraphConfig& config)
{
    performanceLog << timestamp
        << "DirectedAnchorGraph::runResolution begins." << endl;

    // Step 1: Initial unitigification
    unitigifyAll();
    writeSummary(cout);
    writeGfa("DirectedAnchorGraph-initial.gfa");

    // Step 2: MBG-style cleaning phase (if enabled)
    if(config.doCleaning) {
        cout << timestamp << "Starting MBG-style cleaning phase..." << endl;

        // 2a: Remove low-coverage tips (MBG pre-resolve pass)
        auto tipStats = removeLowCoverageTips(
            config.tipMinCoverage,
            double(config.tipMinLength),
            config.tipMaxLength);
        if(tipStats.nodesRemoved > 0) {
            cout << "  Removed " << tipStats.nodesRemoved << " tip nodes." << endl;
            unitigifyAll();
            cout << "  After tip removal: "
                 << nodeCount() << " nodes, "
                 << edgeCount() << " edges." << endl;
        }

        // 2b: Remove low-coverage crosslinks
        auto crosslinkStats = removeLowCoverageCrosslinks(
            config.crosslinkMinCoverage,
            config.crosslinkMinLength);
        if(crosslinkStats.edgesRemoved > 0) {
            cout << "  Removed " << crosslinkStats.edgesRemoved
                 << " crosslink edges." << endl;
            unitigifyAll();
            cout << "  After crosslink removal: "
                 << nodeCount() << " nodes, "
                 << edgeCount() << " edges." << endl;
        }

        // 2c: Copy-number based cleaning (if guesswork enabled)
        if(config.doGuesswork) {
            const double avgCov = [&]() {
                double totalCov = 0.0;
                uint64_t count = 0;
                for(uint64_t segId = 0; segId < nodes.size(); ++segId) {
                    if(nodeExists(segId)) {
                        totalCov += getPathCoverage(segId);
                        count++;
                    }
                }
                return count > 0 ? totalCov / double(count) : 1.0;
            }();

            auto copyStats = cleanComponentsByCopynumber(
                avgCov,
                50000,
                0,
                max(config.maxResolveLength, config.maxLocalResolve),
                {},
                0);
            if(copyStats.nodesRemoved > 0 || copyStats.edgesRemoved > 0) {
                cout << "  Copy-number cleaning removed "
                     << copyStats.nodesRemoved << " nodes, "
                     << copyStats.edgesRemoved << " edges." << endl;
                unitigifyAll();
                cout << "  After copy-number cleaning: "
                     << nodeCount() << " nodes, "
                     << edgeCount() << " edges." << endl;
            }
        }

        cout << timestamp << "Cleaning phase complete." << endl;
        writeSummary(cout);
    }

    // Step 3: MBG-style first resolve pass at configured min support.
    uint64_t initialStep = config.minEdgeSupport;
    if(!config.resolveSteps.empty()) {
        initialStep = config.resolveSteps.front();
    }
    cout << timestamp
         << "Resolution step with minEdgeSupport=" << initialStep << endl;
    resolveRound(
        initialStep,
        config.maxResolveLength,
        config.doCleaning,
        config.doGuesswork,
        config.maxUnconditionalResolveLength,
        config.copycountFilterHeuristic,
        config.maxLocalResolve,
        config.resolvePalindromesGlobal);
    unitigifyAll();
    cout << "  After resolution step " << initialStep << ": ";
    cout << nodeCount() << " nodes, "
         << edgeCount() << " edges, "
         << pathCount() << " paths." << endl;

    if(config.doFinalLowSupportResolvePass) {
        cout << timestamp
             << "Resolution final low-support pass with minEdgeSupport=1" << endl;
        resolveRound(
            1,
            config.maxResolveLength,
            config.doCleaning,
            config.doGuesswork,
            config.maxUnconditionalResolveLength,
            config.copycountFilterHeuristic,
            config.maxLocalResolve,
            config.resolvePalindromesGlobal);
        unitigifyAll();
        cout << "  After final low-support pass: "
             << nodeCount() << " nodes, "
             << edgeCount() << " edges, "
             << pathCount() << " paths." << endl;
    }

    writeSummary(cout);
    writeGfa(config.outputGfa);
    writePaths(config.outputPaths);

    performanceLog << timestamp
        << "DirectedAnchorGraph::runResolution ends." << endl;
}


void DirectedAnchorGraph::verifyEdgeConsistency() const
{
    performanceLog << timestamp << "DirectedAnchorGraph::verifyEdgeConsistency begins." << endl;
    
    uint64_t totalVerified = 0;
    uint64_t totalMismatches = 0;

    for(const auto& [edge, overlap] : edgeOverlaps) {
        const DagNodeId from = edge.first;
        const DagNodeId to = edge.second;
        
        vector<Base> seqFrom = getSegmentSequence(segmentOf(from));
        vector<Base> seqTo = getSegmentSequence(segmentOf(to));
        
        if(seqFrom.empty() || seqTo.empty()) continue;
        
        if(!isForward(from)) {
            reverseComplement(seqFrom);
        }
        if(!isForward(to)) {
            reverseComplement(seqTo);
        }

        if(seqFrom.size() < overlap || seqTo.size() < overlap) {
            totalMismatches++;
            continue;
        }

        const uint64_t fromStart = seqFrom.size() - overlap;
        for(uint64_t i = 0; i < overlap; i++) {
            if(seqFrom[fromStart + i] != seqTo[i]) {
                totalMismatches++;
                break;
            }
        }
        totalVerified++;
    }

    performanceLog << timestamp << "Verified " << totalVerified << " edges. "
                   << totalMismatches << " mismatches found." << endl;
    
    if(totalMismatches > 0) {
        cout << "WARNING: verifyEdgeConsistency found " << totalMismatches << " mismatches!" << endl;
    }
}
