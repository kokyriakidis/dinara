// ============================================================================
// DirectedAnchorGraph — Core graph operations, build from BRG, GFA I/O.
//
// Verkko-style two-node representation:
//   Each segment has two oriented DagNodeIds (even=fwd, odd=rev).
//   rcNode(n) = n ^ 1.  No OrientedDagNode struct needed.
// ============================================================================

#include "mode3-DirectedAnchorGraph.hpp"
#include "mode3-Anchor.hpp"
#include "mode3-DirectedAnchors.hpp"
#include "DINARA_ASSERT.hpp"
#include "Marker.hpp"
#include "MultithreadedObject.tpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace dinara;
using namespace mode3;
using namespace std;

namespace {

uint64_t medianValue(vector<uint64_t>& values)
{
    if(values.empty()) {
        return 0;
    }
    sort(values.begin(), values.end());
    const size_t n = values.size();
    if((n & 1ULL) != 0) {
        return values[n / 2];
    }
    return (values[n / 2 - 1] + values[n / 2]) / 2;
}

bool getAnchorIntervalOnOrientedRead(
    const Anchors& anchors,
    DagNodeId anchorId,
    OrientedReadId orientedReadId,
    uint64_t& start,
    uint64_t& end)
{
    const uint32_t ordinal0 = anchors.getFirstOrdinal(anchorId, orientedReadId);
    const uint32_t ordinal1 = ordinal0 + anchors.ordinalOffset(anchorId);
    const auto orientedMarkers = anchors.markers[orientedReadId.getValue()];
    if(ordinal0 >= orientedMarkers.size() || ordinal1 >= orientedMarkers.size()) {
        return false;
    }

    // Anchors are represented by marker midpoints.
    const uint64_t p0 = uint64_t(orientedMarkers[ordinal0].position) + uint64_t(anchors.kHalf);
    const uint64_t p1 = uint64_t(orientedMarkers[ordinal1].position) + uint64_t(anchors.kHalf);
    if(p1 < p0) {
        return false;
    }
    start = p0;
    end = p1;
    return true;
}

uint64_t intervalOverlap(uint64_t aStart, uint64_t aEnd, uint64_t bStart, uint64_t bEnd)
{
    const uint64_t left = max(aStart, bStart);
    const uint64_t right = min(aEnd, bEnd);
    return (right > left) ? (right - left) : 0;
}

uint64_t mapPathPositionToRead(
    const DagPathReadSupport& read,
    uint64_t pathPos,
    uint64_t pathLen)
{
    const uint64_t span =
        (read.readPosEndIndex > read.readPosStartIndex) ?
        (read.readPosEndIndex - read.readPosStartIndex) : 0;
    const uint64_t coveredStart = min(read.leftClip, pathLen);
    uint64_t coveredEnd = 0;
    if(pathLen > read.rightClip) {
        coveredEnd = pathLen - read.rightClip;
    }
    if(coveredEnd < coveredStart) {
        coveredEnd = coveredStart;
    }
    const uint64_t coveredSpan = coveredEnd - coveredStart;

    if(pathPos <= coveredStart || coveredSpan == 0 || span == 0) {
        return read.readPosStartIndex;
    }
    if(pathPos >= coveredEnd) {
        return read.readPosEndIndex;
    }
    const uint64_t delta = pathPos - coveredStart;
    const uint64_t clamped = min(delta, span);
    return read.readPosStartIndex + clamped;
}

} // namespace


// ============================================================================
// Static members and template instantiation.
// ============================================================================

template class dinara::MultithreadedObject<mode3::DirectedAnchorGraph>;

const vector<DagNodeId> DirectedAnchorGraph::emptyEdgeVec;
const unordered_set<uint64_t> DagPathsCrossingIndex::emptySet;


// ============================================================================
// DagPathsCrossingIndex — indexes by segment ID (= nodeId >> 1).
// ============================================================================

void DagPathsCrossingIndex::addPath(
    uint64_t pathIdx, const vector<DagNodeId>& path)
{
    for(DagNodeId n : path) {
        index[segmentOf(n)].insert(pathIdx);
    }
}

void DagPathsCrossingIndex::removePath(
    uint64_t pathIdx, const vector<DagNodeId>& path)
{
    for(DagNodeId n : path) {
        uint64_t seg = segmentOf(n);
        auto it = index.find(seg);
        if(it != index.end()) {
            it->second.erase(pathIdx);
            if(it->second.empty()) {
                index.erase(it);
            }
        }
    }
}

void DagPathsCrossingIndex::rebuild(
    const vector<vector<DagNodeId>>& paths,
    const vector<bool>& pathRemoved)
{
    index.clear();
    for(uint64_t i = 0; i < paths.size(); ++i) {
        if(!pathRemoved[i]) {
            addPath(i, paths[i]);
        }
    }
}

const unordered_set<uint64_t>& DagPathsCrossingIndex::getPathsCrossing(
    uint64_t segId) const
{
    auto it = index.find(segId);
    if(it != index.end()) {
        return it->second;
    }
    return emptySet;
}

bool DagPathsCrossingIndex::hasPathsCrossing(uint64_t segId) const
{
    auto it = index.find(segId);
    return it != index.end() && !it->second.empty();
}


// ============================================================================
// Graph accessors.
// ============================================================================

uint64_t DirectedAnchorGraph::nodeCount() const
{
    uint64_t count = 0;
    for(const auto& info : nodes) {
        if(!info.removed) ++count;
    }
    return count;
}

uint64_t DirectedAnchorGraph::totalNodeCount() const
{
    return nodes.size();
}

uint64_t DirectedAnchorGraph::edgeCount() const
{
    uint64_t count = 0;
    for(const auto& adjList : edges) {
        count += adjList.size();
    }
    return count;
}

uint64_t DirectedAnchorGraph::pathCount() const
{
    uint64_t count = 0;
    for(uint64_t i = 0; i < paths.size(); ++i) {
        if(!pathRemoved[i]) ++count;
    }
    return count;
}

bool DirectedAnchorGraph::nodeExists(uint64_t segId) const
{
    return segId < nodes.size() && !nodes[segId].removed;
}

const vector<DagNodeId>& DirectedAnchorGraph::getOutEdges(DagNodeId from) const
{
    if(from < edges.size() && !edges[from].empty()) {
        return edges[from];
    }
    return emptyEdgeVec;
}

uint64_t DirectedAnchorGraph::outDegree(DagNodeId from) const
{
    return getOutEdges(from).size();
}

uint64_t DirectedAnchorGraph::inDegree(DagNodeId from) const
{
    // In-degree of n = out-degree of n^1 (by RC invariant).
    return outDegree(from ^ 1);
}


// ============================================================================
// Edge manipulation — maintains the RC invariant via XOR.
// ============================================================================

void DirectedAnchorGraph::addEdge(
    DagNodeId from,
    DagNodeId to,
    uint64_t overlapBp)
{
    // Ensure edges vector is large enough.
    uint64_t maxId = max(from, to ^ 1) + 1;
    if(edges.size() < maxId) {
        edges.resize(maxId);
    }

    const auto canonical = canonEdge(from, to);
    const uint64_t fromLen = nodeBpLength(segmentOf(from));
    const uint64_t toLen = nodeBpLength(segmentOf(to));
    const uint64_t maxOverlap = (fromLen > 0 && toLen > 0) ?
        min(fromLen, toLen) - 1 : 0;
    const uint64_t finalOverlap = min<uint64_t>(overlapBp, maxOverlap);

    // Add forward edge (check for duplicates).
    auto& fwdList = edges[from];
    if(find(fwdList.begin(), fwdList.end(), to) == fwdList.end()) {
        fwdList.push_back(to);
    }

    // Add RC edge (check for duplicates).
    auto& rcList = edges[to ^ 1];
    DagNodeId rcFrom = from ^ 1;
    if(find(rcList.begin(), rcList.end(), rcFrom) == rcList.end()) {
        rcList.push_back(rcFrom);
    }

    auto it = edgeOverlaps.find(canonical);
    if(it == edgeOverlaps.end()) {
        edgeOverlaps.emplace(canonical, finalOverlap);
    } else {
        it->second = max(it->second, finalOverlap);
    }
}

void DirectedAnchorGraph::removeEdge(DagNodeId from, DagNodeId to)
{
    // Remove forward edge (swap-and-pop for O(1) instead of O(degree) erase).
    if(from < edges.size()) {
        auto& fwdList = edges[from];
        auto it = find(fwdList.begin(), fwdList.end(), to);
        if(it != fwdList.end()) {
            *it = fwdList.back();
            fwdList.pop_back();
        }
    }

    // Remove RC edge (swap-and-pop).
    DagNodeId rcTo = to ^ 1;
    if(rcTo < edges.size()) {
        auto& rcList = edges[rcTo];
        DagNodeId rcFrom = from ^ 1;
        auto it = find(rcList.begin(), rcList.end(), rcFrom);
        if(it != rcList.end()) {
            *it = rcList.back();
            rcList.pop_back();
        }
    }
    edgeOverlaps.erase(canonEdge(from, to));
}

void DirectedAnchorGraph::removeAllEdges(uint64_t segId)
{
    DagNodeId fwd = fwdNodeId(segId);
    DagNodeId rev = revNodeId(segId);

    // Pop edges from the back to avoid copying the neighbor list.
    while(fwd < edges.size() && !edges[fwd].empty()) {
        removeEdge(fwd, edges[fwd].back());
    }
    while(rev < edges.size() && !edges[rev].empty()) {
        removeEdge(rev, edges[rev].back());
    }
}

uint64_t DirectedAnchorGraph::getBpOverlap(DagNodeId from, DagNodeId to) const
{
    auto it = edgeOverlaps.find(canonEdge(from, to));
    if(it == edgeOverlaps.end()) return 0;
    return it->second;
}

void DirectedAnchorGraph::setBpOverlap(DagNodeId from, DagNodeId to, uint64_t overlapBp)
{
    const uint64_t fromLen = nodeBpLength(segmentOf(from));
    const uint64_t toLen = nodeBpLength(segmentOf(to));
    const uint64_t maxOverlap = (fromLen > 0 && toLen > 0) ?
        min(fromLen, toLen) - 1 : 0;
    edgeOverlaps[canonEdge(from, to)] = min<uint64_t>(overlapBp, maxOverlap);
}

uint64_t DirectedAnchorGraph::nodeBpLength(uint64_t segId) const
{
    if(segId >= nodes.size()) return 1;
    return max<uint64_t>(1, nodes[segId].lengthBp);
}

uint64_t DirectedAnchorGraph::pathBpLength(const vector<DagNodeId>& path) const
{
    if(path.empty()) return 0;
    uint64_t total = 0;
    for(size_t i = 0; i < path.size(); i++) {
        const uint64_t segId = segmentOf(path[i]);
        total += nodeBpLength(segId);
        if(i > 0) {
            const uint64_t ov = getBpOverlap(path[i - 1], path[i]);
            if(total >= ov) total -= ov;
            else total = 0;
        }
    }
    return total;
}

vector<uint64_t> DirectedAnchorGraph::computePathNodeBoundaries(
    const vector<DagNodeId>& path) const
{
    vector<uint64_t> boundaries;
    if(path.empty()) return boundaries;
    boundaries.reserve(path.size() + 1);
    uint64_t running = 0;
    boundaries.push_back(0);
    for(size_t i = 0; i < path.size(); i++) {
        running += nodeBpLength(segmentOf(path[i]));
        if(i > 0) {
            const uint64_t ov = getBpOverlap(path[i - 1], path[i]);
            if(running >= ov) {
                running -= ov;
            } else {
                running = 0;
            }
        }
        boundaries.push_back(running);
    }
    return boundaries;
}

bool DirectedAnchorGraph::reconcileReadSupportToPath(
    DagPathReadSupport& read,
    uint64_t bpLen) const
{
    if(read.readPosEndIndex < read.readPosStartIndex) {
        read.readPosEndIndex = read.readPosStartIndex;
    }

    uint64_t span = read.readPosEndIndex - read.readPosStartIndex;
    uint64_t total = span + read.leftClip + read.rightClip;

    if(total < bpLen) {
        read.rightClip += (bpLen - total);
    } else if(total > bpLen) {
        uint64_t extra = total - bpLen;
        const uint64_t cutRight = min<uint64_t>(extra, read.rightClip);
        read.rightClip -= cutRight;
        extra -= cutRight;

        const uint64_t cutLeft = min<uint64_t>(extra, read.leftClip);
        read.leftClip -= cutLeft;
        extra -= cutLeft;

        if(extra > 0) {
            span = read.readPosEndIndex - read.readPosStartIndex;
            if(span <= extra) {
                return false;
            }
            read.readPosEndIndex -= extra;
        }
    }

    span = read.readPosEndIndex - read.readPosStartIndex;
    total = span + read.leftClip + read.rightClip;
    if(total < bpLen) {
        read.rightClip += (bpLen - total);
    } else if(total > bpLen) {
        const uint64_t extra = total - bpLen;
        if(read.rightClip < extra) {
            return false;
        }
        read.rightClip -= extra;
    }
    return read.readPosEndIndex > read.readPosStartIndex;
}

bool DirectedAnchorGraph::hasUsableReadPoses(
    const DagPathReadSupport& read,
    const vector<uint64_t>& pathBoundaries) const
{
    if(read.readPoses.size() != pathBoundaries.size() || pathBoundaries.empty()) {
        return false;
    }
    for(size_t i = 1; i < read.readPoses.size(); i++) {
        if(read.readPoses[i] < read.readPoses[i - 1]) {
            return false;
        }
    }
    return true;
}

uint64_t DirectedAnchorGraph::mapPathPositionFromSupport(
    const DagPathReadSupport& read,
    const vector<uint64_t>& pathBoundaries,
    uint64_t pathPos) const
{
    if(hasUsableReadPoses(read, pathBoundaries)) {
        if(pathPos <= pathBoundaries.front()) {
            return read.readPoses.front();
        }
        if(pathPos >= pathBoundaries.back()) {
            return read.readPoses.back();
        }
        auto exactIt = lower_bound(pathBoundaries.begin(), pathBoundaries.end(), pathPos);
        if(exactIt != pathBoundaries.end() && *exactIt == pathPos) {
            return read.readPoses[size_t(exactIt - pathBoundaries.begin())];
        }
        auto upperIt = upper_bound(pathBoundaries.begin(), pathBoundaries.end(), pathPos);
        if(upperIt == pathBoundaries.begin() || upperIt == pathBoundaries.end()) {
            return read.readPoses.back();
        }
        const size_t idx = size_t(upperIt - pathBoundaries.begin() - 1);
        const uint64_t p0 = pathBoundaries[idx];
        const uint64_t p1 = pathBoundaries[idx + 1];
        const uint64_t r0 = read.readPoses[idx];
        const uint64_t r1 = read.readPoses[idx + 1];
        if(p1 <= p0 || r1 <= r0) {
            return r0;
        }
        const uint64_t num = (pathPos - p0) * (r1 - r0);
        return r0 + (num / (p1 - p0));
    }

    const uint64_t pathLen = pathBoundaries.empty() ? 0 : pathBoundaries.back();
    return mapPathPositionToRead(read, pathPos, pathLen);
}

optional<DagPathReadSupport> DirectedAnchorGraph::projectReadSupportToBoundaries(
    const DagPathReadSupport& read,
    const vector<uint64_t>& oldPathBoundaries,
    const vector<uint64_t>& projectedBoundaries) const
{
    if(oldPathBoundaries.empty() || projectedBoundaries.size() < 2) {
        return nullopt;
    }
    for(size_t i = 1; i < projectedBoundaries.size(); i++) {
        if(projectedBoundaries[i] < projectedBoundaries[i - 1]) {
            return nullopt;
        }
    }

    const uint64_t oldPathLen = oldPathBoundaries.back();
    const uint64_t rawStart = min(projectedBoundaries.front(), oldPathLen);
    const uint64_t rawEnd = min(projectedBoundaries.back(), oldPathLen);
    if(rawEnd <= rawStart) {
        return nullopt;
    }

    const uint64_t coveredStart = min(read.leftClip, oldPathLen);
    uint64_t coveredEnd = 0;
    if(oldPathLen > read.rightClip) {
        coveredEnd = oldPathLen - read.rightClip;
    }
    if(coveredEnd < coveredStart) {
        coveredEnd = coveredStart;
    }

    const uint64_t overlapStart = max(rawStart, coveredStart);
    const uint64_t overlapEnd = min(rawEnd, coveredEnd);
    if(overlapEnd <= overlapStart) {
        return nullopt;
    }

    DagPathReadSupport out = read;
    out.leftClip = overlapStart - rawStart;
    out.rightClip = rawEnd - overlapEnd;
    out.readPosStartIndex = mapPathPositionFromSupport(read, oldPathBoundaries, overlapStart);
    out.readPosEndIndex = mapPathPositionFromSupport(read, oldPathBoundaries, overlapEnd);
    if(out.readPosEndIndex <= out.readPosStartIndex) {
        return nullopt;
    }

    out.readPoses.resize(projectedBoundaries.size());
    for(size_t i = 0; i < projectedBoundaries.size(); i++) {
        const uint64_t p = projectedBoundaries[i];
        const uint64_t clamped = min(max(p, overlapStart), overlapEnd);
        out.readPoses[i] =
            mapPathPositionFromSupport(read, oldPathBoundaries, clamped);
    }
    return out;
}

bool DirectedAnchorGraph::isReadSupportStrictForPath(
    const DagPathReadSupport& read,
    const vector<uint64_t>& pathBoundaries,
    uint64_t pathBpLength) const
{
    DagPathReadSupport reconciled = read;
    if(!reconcileReadSupportToPath(reconciled, pathBpLength)) {
        return false;
    }
    if(reconciled.leftClip != read.leftClip ||
       reconciled.rightClip != read.rightClip ||
       reconciled.readPosStartIndex != read.readPosStartIndex ||
       reconciled.readPosEndIndex != read.readPosEndIndex) {
        return false;
    }
    if(!hasUsableReadPoses(read, pathBoundaries)) {
        return false;
    }
    if(read.readPoses.front() != read.readPosStartIndex ||
       read.readPoses.back() != read.readPosEndIndex) {
        return false;
    }
    for(uint64_t p : read.readPoses) {
        if(p < read.readPosStartIndex || p > read.readPosEndIndex) {
            return false;
        }
    }
    return true;
}

void DirectedAnchorGraph::refreshReadPosesForPath(uint64_t pathIdx)
{
    if(pathIdx >= paths.size() || pathIdx >= pathReadIds.size()) {
        return;
    }
    if(pathRemoved[pathIdx] || paths[pathIdx].empty()) {
        pathReadIds[pathIdx].clear();
        if(pathIdx < pathWeights.size()) {
            pathWeights[pathIdx] = 0;
        }
        return;
    }

    const vector<uint64_t> boundaries = computePathNodeBoundaries(paths[pathIdx]);
    const uint64_t bpLen = boundaries.empty() ? 0 : boundaries.back();

    vector<DagPathReadSupport> keptReads;
    keptReads.reserve(pathReadIds[pathIdx].size());
    for(auto read : pathReadIds[pathIdx]) {
        if(!reconcileReadSupportToPath(read, bpLen)) {
            continue;
        }
        if(!hasUsableReadPoses(read, boundaries) ||
           read.readPoses.front() != read.readPosStartIndex ||
           read.readPoses.back() != read.readPosEndIndex) {
            read.readPoses.resize(boundaries.size());
            for(size_t i = 0; i < boundaries.size(); i++) {
                read.readPoses[i] = mapPathPositionToRead(read, boundaries[i], bpLen);
            }
        }
        keptReads.push_back(std::move(read));
    }
    pathReadIds[pathIdx] = std::move(keptReads);
    if(pathIdx < pathWeights.size()) {
        pathWeights[pathIdx] = max<uint64_t>(1, pathReadIds[pathIdx].size());
    }
}


// ============================================================================
// Segment manipulation.
// ============================================================================

uint64_t DirectedAnchorGraph::addNode(const DagNodeInfo& info)
{
    uint64_t segId = nodes.size();
    nodes.push_back(info);

    // Ensure edges vector has space for both orientations.
    uint64_t requiredEdgeSize = (segId + 1) * 2;
    if(edges.size() < requiredEdgeSize) {
        edges.resize(requiredEdgeSize);
    }

    return segId;
}

void DirectedAnchorGraph::removeNode(uint64_t segId)
{
    if(segId >= nodes.size()) return;

    nodes[segId].removed = true;
    removeAllEdges(segId);
}


// ============================================================================
// Path manipulation.
// ============================================================================

uint64_t DirectedAnchorGraph::addPath(
    vector<DagNodeId> path,
    uint64_t weight,
    vector<DagPathReadSupport> readSupport)
{
    uint64_t idx = paths.size();
    pathsCrossing.addPath(idx, path);
    paths.push_back(std::move(path));
    pathRemoved.push_back(false);
    pathWeights.push_back(weight);
    const vector<uint64_t> boundaries = computePathNodeBoundaries(paths[idx]);
    const uint64_t bpLen = boundaries.empty() ? 0 : boundaries.back();
    if(readSupport.empty()) {
        readSupport.resize(weight);
        for(uint64_t k = 0; k < weight; k++) {
            readSupport[k].readId = (idx << 32) ^ k;
            readSupport[k].readNameIndex = readSupport[k].readId;
            readSupport[k].readInfoIndex = readSupport[k].readId;
            readSupport[k].readPosStartIndex = 0;
            readSupport[k].readPosEndIndex = bpLen;
            readSupport[k].readPoses = boundaries;
        }
    } else {
        for(size_t i = 0; i < readSupport.size(); i++) {
            auto& read = readSupport[i];
            read.readNameIndex = read.readId;
            read.readInfoIndex = read.readId;
            if(!isReadSupportStrictForPath(read, boundaries, bpLen)) {
                // Try to reconcile the read support instead of failing.
                if(!reconcileReadSupportToPath(read, bpLen)) {
                    throw runtime_error(
                        "DirectedAnchorGraph::addPath received totally inconsistent read support " +
                        to_string(i) + " for path " + to_string(idx));
                }
                // Ensure readPoses are consistent with the reconciled indices.
                if(!hasUsableReadPoses(read, boundaries) ||
                   read.readPoses.front() != read.readPosStartIndex ||
                   read.readPoses.back() != read.readPosEndIndex) {
                    read.readPoses.resize(boundaries.size());
                    for(size_t j = 0; j < boundaries.size(); j++) {
                        read.readPoses[j] = mapPathPositionToRead(read, boundaries[j], bpLen);
                    }
                }
            }
        }
        pathWeights[idx] = readSupport.size();
    }
    pathReadIds.push_back(std::move(readSupport));
    return idx;
}

void DirectedAnchorGraph::removePath(uint64_t pathIdx)
{
    if(pathIdx >= paths.size() || pathRemoved[pathIdx]) return;
    pathsCrossing.removePath(pathIdx, paths[pathIdx]);
    pathRemoved[pathIdx] = true;
    paths[pathIdx].clear();
}

const vector<DagNodeId>& DirectedAnchorGraph::getPath(uint64_t pathIdx) const
{
    return paths[pathIdx];
}

void DirectedAnchorGraph::setPath(uint64_t pathIdx, vector<DagNodeId>&& newPath)
{
    DINARA_ASSERT(pathIdx < paths.size());
    DINARA_ASSERT(!pathRemoved[pathIdx]);
    DINARA_ASSERT(pathIdx < pathWeights.size());
    DINARA_ASSERT(pathIdx < pathReadIds.size());

    const vector<uint64_t> oldBoundaries = computePathNodeBoundaries(paths[pathIdx]);
    const vector<uint64_t> newBoundaries = computePathNodeBoundaries(newPath);
    if(oldBoundaries != newBoundaries) {
        throw runtime_error(
            "DirectedAnchorGraph::setPath changed path boundaries for path " +
            to_string(pathIdx) +
            ". Use resolution remap functions instead of setPath for support-carrying paths.");
    }

    pathsCrossing.removePath(pathIdx, paths[pathIdx]);
    paths[pathIdx] = std::move(newPath);
    const uint64_t bpLen = newBoundaries.empty() ? 0 : newBoundaries.back();
    for(size_t i = 0; i < pathReadIds[pathIdx].size(); i++) {
        if(!isReadSupportStrictForPath(pathReadIds[pathIdx][i], newBoundaries, bpLen)) {
            throw runtime_error(
                "DirectedAnchorGraph::setPath found inconsistent existing read support " +
                to_string(i) + " for path " + to_string(pathIdx));
        }
    }
    pathsCrossing.addPath(pathIdx, paths[pathIdx]);
}


// ============================================================================
// buildFromAnchors — multithreaded build from mode3::Anchors.
//
// Phase 1 (sequential): Create segment metadata (nodes map).
// Phase 2 (parallel):   Each thread processes a batch of strand-0 journeys,
//                        collecting paths and edge pairs in per-thread vectors.
// Phase 3 (sequential): Merge per-thread paths into global paths vector,
//                        insert all collected edge pairs.
// Phase 4 (sequential): Build paths-crossing index.
// ============================================================================

void DirectedAnchorGraph::buildFromAnchors(
    const Anchors& anchors,
    uint64_t threadCount,
    const vector<DagReadInterval>* readIntervals)
{
    performanceLog << timestamp
        << "DirectedAnchorGraph::buildFromAnchors begins with "
        << threadCount << " threads." << endl;

    nodes.clear();
    edges.clear();
    edgeOverlaps.clear();
    paths.clear();
    pathRemoved.clear();
    pathWeights.clear();
    pathReadIds.clear();
    pathsCrossing.clear();
    everTippable.clear();
    lastTippableChecked = 0;
    anchorsPtr = &anchors;

    const uint64_t totalAnchors = anchors.size();
    const uint64_t segmentCount = totalAnchors / 2;

    // Phase 1: One segment per anchor pair (sequential).
    nodes.clear();
    nodes.reserve(segmentCount);
    for(uint64_t seg = 0; seg < segmentCount; ++seg) {
        DagNodeInfo info;
        DagNodeId fwdAnchor = fwdNodeId(seg);
        info.anchorChain.push_back(fwdAnchor);
        info.coverage = double(anchors[fwdAnchor].coverage());
        info.removed = false;

        // Calibrate node length from anchor marker-interval geometry.
        vector<uint64_t> observedLengths;
        const auto markerIntervals = anchors[fwdAnchor];
        observedLengths.reserve(markerIntervals.size());
        for(const auto& markerInterval : markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            const auto orientedMarkers = anchors.markers[orientedReadId.getValue()];
            const uint32_t ordinal0 = markerInterval.ordinal0;
            const uint32_t ordinal1 = ordinal0 + anchors.ordinalOffset(fwdAnchor);
            if(ordinal0 >= orientedMarkers.size() || ordinal1 >= orientedMarkers.size()) {
                continue;
            }
            const uint64_t p0 =
                uint64_t(orientedMarkers[ordinal0].position) + uint64_t(anchors.kHalf);
            const uint64_t p1 =
                uint64_t(orientedMarkers[ordinal1].position) + uint64_t(anchors.kHalf);
            if(p1 >= p0) {
                observedLengths.push_back(p1 - p0);
            }
        }
        info.lengthBp = max<uint64_t>(1, medianValue(observedLengths));
        nodes.push_back(std::move(info));
    }

    // Pre-size edges vector for all oriented nodes (2 per segment).
    edges.clear();
    edges.resize(segmentCount * 2);

    performanceLog << timestamp
        << "Phase 1 done: created " << segmentCount << " segments." << endl;

    // Build the list of strand-0 journey indices with >= 2 steps.
    strand0Indices.clear();
    strand0Indices.reserve(anchors.journeys.size() / 2);
    for(uint64_t oridValue = 0; oridValue < anchors.journeys.size(); oridValue += 2) {
        if(anchors.journeys[oridValue].size() >= 2) {
            strand0Indices.push_back(oridValue);
        }
    }

    performanceLog << timestamp
        << "Found " << strand0Indices.size()
        << " strand-0 journeys with >= 2 steps." << endl;

    // Phase 2: Parallel journey processing.
    buildAnchorsPtr = &anchors;
    buildReadIntervalsPtr = readIntervals;
    threadPaths.resize(threadCount);
    threadPathReadSupport.resize(threadCount);
    threadEdgePairs.resize(threadCount);
    threadEdgeOverlapSamples.resize(threadCount);
    for(uint64_t t = 0; t < threadCount; ++t) {
        threadPaths[t].clear();
        threadPathReadSupport[t].clear();
        threadEdgePairs[t].clear();
        threadEdgeOverlapSamples[t].clear();
    }

    const uint64_t batchSize = max(uint64_t(1), uint64_t(strand0Indices.size() / (threadCount * 100)));
    setupLoadBalancing(strand0Indices.size(), batchSize);
    runThreads(&DirectedAnchorGraph::buildFromAnchorsThreadFunction, threadCount);

    performanceLog << timestamp
        << "Phase 2 done: parallel journey processing complete." << endl;

    // Phase 3: Merge per-thread results (sequential).

    // 3a: Count total paths for reservation.
    uint64_t totalPaths = 0;
    for(uint64_t t = 0; t < threadCount; ++t) {
        totalPaths += threadPaths[t].size();
    }
    paths.clear();
    pathRemoved.clear();
    pathWeights.clear();
    pathReadIds.clear();

    // MBG PathGroup-like grouping: collapse identical paths and store multiplicity.
    map<vector<DagNodeId>, vector<DagPathReadSupport>> groupedPaths;
    for(uint64_t t = 0; t < threadCount; ++t) {
        DINARA_ASSERT(threadPaths[t].size() == threadPathReadSupport[t].size());
        for(size_t i = 0; i < threadPaths[t].size(); i++) {
            const auto& path = threadPaths[t][i];
            groupedPaths[path].push_back(threadPathReadSupport[t][i]);
        }
        threadPaths[t].clear();
        threadPaths[t].shrink_to_fit();
        threadPathReadSupport[t].clear();
        threadPathReadSupport[t].shrink_to_fit();
    }

    paths.reserve(groupedPaths.size());
    pathRemoved.reserve(groupedPaths.size());
    pathWeights.reserve(groupedPaths.size());
    pathReadIds.reserve(groupedPaths.size());
    for(const auto& [path, support] : groupedPaths) {
        paths.push_back(path);
        pathRemoved.push_back(false);
        pathWeights.push_back(support.size());
        pathReadIds.push_back(support);
    }
    for(uint64_t pathIdx = 0; pathIdx < paths.size(); pathIdx++) {
        refreshReadPosesForPath(pathIdx);
    }

    // 3b: Insert edges using temporary sets for deduplication, then convert to sorted vectors.
    // Use unordered_set during build to handle duplicates efficiently.
    vector<unordered_set<DagNodeId>> tempEdgeSets(edges.size());
    unordered_map<pair<DagNodeId, DagNodeId>, vector<uint64_t>, DagEdgePairHash> edgeOverlapSamples;

    uint64_t totalEdgePairs = 0;
    for(uint64_t t = 0; t < threadCount; ++t) {
        totalEdgePairs += threadEdgePairs[t].size();
        for(const auto& [from, to] : threadEdgePairs[t]) {
            // Bounds check to prevent crashes
            if(from >= tempEdgeSets.size() || to >= tempEdgeSets.size()) {
                throw runtime_error("Edge ID out of range: from=" + std::to_string(from) +
                                    " to=" + std::to_string(to) +
                                    " size=" + std::to_string(tempEdgeSets.size()));
            }
            tempEdgeSets[from].insert(to);
            tempEdgeSets[to ^ 1].insert(from ^ 1);
        }
        for(const auto& sample : threadEdgeOverlapSamples[t]) {
            const DagNodeId from = get<0>(sample);
            const DagNodeId to = get<1>(sample);
            const uint64_t overlapBp = get<2>(sample);
            edgeOverlapSamples[canonEdge(from, to)].push_back(overlapBp);
        }
        threadEdgePairs[t].clear();
        threadEdgePairs[t].shrink_to_fit();
        threadEdgeOverlapSamples[t].clear();
        threadEdgeOverlapSamples[t].shrink_to_fit();
    }

    // Convert temporary sets to sorted vectors.
    for(uint64_t i = 0; i < tempEdgeSets.size(); ++i) {
        if(!tempEdgeSets[i].empty()) {
            edges[i].assign(tempEdgeSets[i].begin(), tempEdgeSets[i].end());
            sort(edges[i].begin(), edges[i].end());
            for(DagNodeId to : edges[i]) {
                const auto canonical = canonEdge(i, to);
                if(edgeOverlaps.find(canonical) != edgeOverlaps.end()) {
                    continue;
                }
                uint64_t overlapBp = 0;
                auto it = edgeOverlapSamples.find(canonical);
                if(it != edgeOverlapSamples.end() && !it->second.empty()) {
                    overlapBp = medianValue(it->second);
                }
                const uint64_t fromLen = nodeBpLength(segmentOf(canonical.first));
                const uint64_t toLen = nodeBpLength(segmentOf(canonical.second));
                const uint64_t maxOverlap =
                    (fromLen > 0 && toLen > 0) ? (min(fromLen, toLen) - 1) : 0;
                edgeOverlaps.emplace(canonical, min(overlapBp, maxOverlap));
            }
        }
    }

    // Clean up thread-local storage.
    threadPaths.clear();
    threadPathReadSupport.clear();
    threadEdgePairs.clear();
    threadEdgeOverlapSamples.clear();
    strand0Indices.clear();
    buildAnchorsPtr = nullptr;
    buildReadIntervalsPtr = nullptr;

    performanceLog << timestamp
        << "Phase 3 done: merged " << totalPaths << " paths, "
        << totalEdgePairs << " edge pairs." << endl;

    // Phase 4: Build paths-crossing index.
    pathsCrossing.rebuild(paths, pathRemoved);

    performanceLog << timestamp
        << "DirectedAnchorGraph::buildFromAnchors ends. "
        << nodeCount() << " nodes, "
        << edgeCount() << " directed edges, "
        << pathCount() << " path groups (" << totalPaths << " raw paths)." << endl;
}


// Thread function for Phase 2: process journey batches.
void DirectedAnchorGraph::buildFromAnchorsThreadFunction(size_t threadId)
{
    const Anchors& anchors = *buildAnchorsPtr;
    const vector<DagReadInterval>* readIntervals = buildReadIntervalsPtr;
    auto& myPaths = threadPaths[threadId];
    auto& myPathSupport = threadPathReadSupport[threadId];
    auto& myEdges = threadEdgePairs[threadId];
    auto& myOverlapSamples = threadEdgeOverlapSamples[threadId];

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t idx = begin; idx < end; ++idx) {
            const uint64_t oridValue = strand0Indices[idx];
            const auto journey = anchors.journeys[oridValue];

            vector<DagNodeId> path;
            path.reserve(journey.size());
            for(uint64_t i = 0; i < journey.size(); ++i) {
                path.push_back(DagNodeId(journey[i]));
            }

            const OrientedReadId orientedReadId =
                OrientedReadId::fromValue(oridValue);

            // Collect edge pairs (forward direction only — merge does RC)
            // and per-read overlap observations for each edge.
            for(uint64_t i = 0; i + 1 < path.size(); ++i) {
                myEdges.push_back({path[i], path[i + 1]});
                uint64_t aStart = 0;
                uint64_t aEnd = 0;
                uint64_t bStart = 0;
                uint64_t bEnd = 0;
                if(getAnchorIntervalOnOrientedRead(
                        anchors,
                        journey[i],
                        orientedReadId,
                        aStart,
                        aEnd) &&
                   getAnchorIntervalOnOrientedRead(
                        anchors,
                        journey[i + 1],
                        orientedReadId,
                        bStart,
                        bEnd)) {
                    myOverlapSamples.emplace_back(
                        path[i],
                        path[i + 1],
                        intervalOverlap(aStart, aEnd, bStart, bEnd));
                }
            }

            myPaths.push_back(std::move(path));
            DagPathReadSupport support;
            support.readId = oridValue / 2; // strand-0 journey index -> read index.
            support.readNameIndex = support.readId;
            support.readInfoIndex = support.readId;
            support.readPosStartIndex = 0;

            vector<uint64_t> nodeStartsAbs;
            vector<uint64_t> nodeEndsAbs;
            nodeStartsAbs.resize(journey.size(), 0);
            nodeEndsAbs.resize(journey.size(), 0);
            bool haveAllIntervals = !journey.empty();
            for(size_t i = 0; i < journey.size(); i++) {
                if(!getAnchorIntervalOnOrientedRead(
                        anchors,
                        journey[i],
                        orientedReadId,
                        nodeStartsAbs[i],
                        nodeEndsAbs[i])) {
                    haveAllIntervals = false;
                    break;
                }
            }

            uint64_t pathBp = myPaths.back().size();
            uint64_t pathStartAbs = 0;
            if(haveAllIntervals &&
                nodeEndsAbs.back() > nodeStartsAbs.front()) {
                pathStartAbs = nodeStartsAbs.front();
                pathBp = nodeEndsAbs.back() - pathStartAbs;
                support.readPosZeroOffset = pathStartAbs;
                support.readPosStartIndex = 0;
                support.readPosEndIndex = pathBp;
                support.readPoses.resize(journey.size() + 1);
                for(size_t i = 0; i < journey.size(); i++) {
                    support.readPoses[i] = nodeStartsAbs[i] - pathStartAbs;
                }
                support.readPoses.back() = nodeEndsAbs.back() - pathStartAbs;
            } else {
                support.readPosEndIndex = pathBp;
                support.readPosZeroOffset = 0;
            }

            // If alignment envelopes are available, anchor this support on them.
            if(readIntervals && support.readId < readIntervals->size()) {
                const DagReadInterval& interval = (*readIntervals)[support.readId];
                if(interval.valid) {
                    if(haveAllIntervals && pathStartAbs >= interval.start) {
                        support.readPosZeroOffset = interval.start;
                        support.readPosStartIndex = pathStartAbs - interval.start;
                        support.readPosEndIndex = support.readPosStartIndex + pathBp;
                        if(!support.readPoses.empty()) {
                            const uint64_t shift = support.readPosStartIndex;
                            for(auto& p : support.readPoses) {
                                p += shift;
                            }
                        }
                    } else {
                        support.readPosZeroOffset = haveAllIntervals ? pathStartAbs : 0;
                        support.readPosStartIndex = 0;
                        support.readPosEndIndex = pathBp;
                    }
                }
            }
            myPathSupport.push_back(support);
        }
    }
}


DirectedAnchorGraph::CleaningStats DirectedAnchorGraph::removeLowCoverageSegments(double minCoverage)
{
    CleaningStats stats;
    for(uint64_t segId = 0; segId < nodes.size(); ++segId) {
        if(nodes[segId].removed) continue;
        if(getPathCoverage(segId) < minCoverage) {
            removeNode(segId);
            stats.nodesRemoved++;
        }
    }
    if(stats.nodesRemoved > 0) {
        unitigifyAll();
    }
    return stats;
}


vector<Base> DirectedAnchorGraph::getSegmentSequence(uint64_t segId) const
{
    if(!anchorsPtr || segId >= nodes.size()) {
        return {};
    }
    const auto& chain = nodes[segId].anchorChain;
    if(chain.empty()) return {};

    vector<Base> sequence;
    for(size_t i = 0; i < chain.size(); i++) {
        const DagNodeId aid = chain[i];
        const span<const Base> anchorSeq = anchorsPtr->anchorSequence(aid);
        
        // Initial anchor: get the whole sequence.
        if(i == 0) {
            sequence.insert(sequence.end(), anchorSeq.begin(), anchorSeq.end());
        } else {
            // Subsequent anchors: the overlap with the previous anchor is already accounted for
            // if we use anchorSequence() which represents the span between midpoints.
            // Wait, if anchorSequence(aid) is the sequence between midpoint(M0) and midpoint(M1)...
            // and aid = (M0, M1).
            // If they are adjacent: chain = [(M0, M1), (M1, M2), (M2, M3)]
            // sequence = seq(M0, M1) + seq(M1, M2) + seq(M2, M3).
            // This is exactly what we want.
            sequence.insert(sequence.end(), anchorSeq.begin(), anchorSeq.end());
        }
    }
    return sequence;
}


// ============================================================================
// GFA I/O.
// ============================================================================

void DirectedAnchorGraph::writeGfa(const string& filename, bool omitSequences) const
{
    ofstream gfa(filename);
    if(!gfa) {
        throw runtime_error("Cannot open " + filename + " for writing.");
    }

    gfa << "H\tVN:Z:1.0\n";

    // S lines — one per segment.
    for(uint64_t segId = 0; segId < nodes.size(); ++segId) {
        const auto& info = nodes[segId];
        if(info.removed) continue;
        
        string seqStr = "*";
        uint64_t len = info.lengthBp;

        if(!omitSequences) {
            const vector<Base> seq = getSegmentSequence(segId);
            if(!seq.empty()) {
                seqStr = string(reinterpret_cast<const char*>(seq.data()), seq.size());
                len = seq.size();
            }
        }

        gfa << "S\t" << segId
            << "\t" << seqStr
            << "\tLN:i:" << max(len, uint64_t(1))
            << "\tRC:i:" << uint64_t(info.coverage)
            << "\n";
    }

    // L lines — one per canonical edge.
    set<pair<DagNodeId, DagNodeId>> writtenEdges;
    for(DagNodeId from = 0; from < edges.size(); ++from) {
        if(!nodeExists(segmentOf(from))) continue;
        const auto& toList = edges[from];
        for(DagNodeId to : toList) {
            if(!nodeExists(segmentOf(to))) continue;
            auto ce = canonEdge(from, to);
            if(writtenEdges.count(ce)) continue;
            writtenEdges.insert(ce);

            gfa << "L\t"
                << segmentOf(ce.first) << "\t"
                << (isForward(ce.first) ? "+" : "-") << "\t"
                << segmentOf(ce.second) << "\t"
                << (isForward(ce.second) ? "+" : "-") << "\t"
                << getBpOverlap(ce.first, ce.second) << "M\n";
        }
    }

    performanceLog << timestamp
        << "Wrote GFA to " << filename
        << " (" << nodeCount() << " segments, "
        << writtenEdges.size() << " links)." << endl;
}

void DirectedAnchorGraph::writePaths(const string& filename) const
{
    ofstream gaf(filename);
    if(!gaf) {
        throw runtime_error("Cannot open " + filename + " for writing.");
    }

    for(uint64_t i = 0; i < paths.size(); ++i) {
        if(pathRemoved[i]) continue;
        const auto& path = paths[i];
        if(path.empty()) continue;

        gaf << i << "\t";
        for(DagNodeId n : path) {
            gaf << (isForward(n) ? ">" : "<") << segmentOf(n);
        }
        gaf << "\n";
    }
}


// ============================================================================
// Diagnostics.
// ============================================================================

void DirectedAnchorGraph::writeSummary(ostream& out) const
{
    out << "DirectedAnchorGraph summary:\n"
        << "  Nodes: " << nodeCount() << " (total with removed: " << totalNodeCount() << ")\n"
        << "  Directed edges: " << edgeCount() << "\n"
        << "  Paths: " << pathCount() << "\n";

    map<uint64_t, uint64_t> outDegreeDist;
    map<uint64_t, uint64_t> inDegreeDist;
    for(uint64_t segId = 0; segId < nodes.size(); ++segId) {
        const auto& info = nodes[segId];
        if(info.removed) continue;
        DagNodeId fwd = fwdNodeId(segId);
        DagNodeId rev = revNodeId(segId);
        outDegreeDist[outDegree(fwd)]++;
        outDegreeDist[outDegree(rev)]++;
        inDegreeDist[inDegree(fwd)]++;
        inDegreeDist[inDegree(rev)]++;
    }
    out << "  Out-degree distribution:\n";
    for(const auto& [deg, cnt] : outDegreeDist) {
        out << "    " << deg << ": " << cnt << "\n";
    }
    out << "  In-degree distribution:\n";
    for(const auto& [deg, cnt] : inDegreeDist) {
        out << "    " << deg << ": " << cnt << "\n";
    }
}


vector<uint64_t> DirectedAnchorGraph::getActiveNodeIds() const
{
    vector<uint64_t> result;
    result.reserve(nodes.size());
    for(uint64_t segId = 0; segId < nodes.size(); ++segId) {
        if(!nodes[segId].removed) {
            result.push_back(segId);
        }
    }
    // Already sorted since we iterate in order.
    return result;
}

// ============================================================================
// buildFromAnchors (Optimized version for DirectedAnchors)
//
// This is the primary entry point for MBG-style graph construction.
// It initializes the graph structure and spawns threads to process journeys.
// ============================================================================
void DirectedAnchorGraph::buildFromAnchors(
    const DirectedAnchors& anchors,
    uint64_t threadCount)
{
    performanceLog << timestamp << "DirectedAnchorGraph::buildFromAnchors (DirectedAnchors) begins." << endl;
    
    // Clear existing graph state and set the new anchor pointer.
    nodes.clear();
    edges.clear();
    edgeOverlaps.clear();
    paths.clear();
    pathRemoved.clear();
    pathWeights.clear();
    pathReadIds.clear();
    pathsCrossing.clear();
    anchorsPtr = &anchors;
    buildDirectedAnchorsPtr = &anchors;
    
    // The rest of the initialization follows.
    everTippable.clear();
    lastTippableChecked = 0;

    const uint64_t segmentCount = anchors.size() / 2;

    nodes.clear();
    nodes.reserve(segmentCount);
    for(uint64_t seg = 0; seg < segmentCount; ++seg) {
        DagNodeInfo info;
        DagNodeId fwdAnchor = fwdNodeId(seg);
        info.anchorChain.push_back(fwdAnchor);
        info.coverage = double(anchors[fwdAnchor].coverage());
        info.removed = false;

        vector<uint64_t> observedLengths;
        const auto markerIntervals = anchors[fwdAnchor];
        const uint32_t offset = anchors.ordinalOffset(fwdAnchor);
        for(const auto& mi : markerIntervals) {
            const auto& markers = anchors.markers[mi.orientedReadId.getValue()];
            if(mi.ordinal0 + offset < markers.size()) {
                observedLengths.push_back(markers[mi.ordinal0 + offset].position - markers[mi.ordinal0].position);
            }
        }
        info.lengthBp = max<uint64_t>(1, medianValue(observedLengths));
        nodes.push_back(std::move(info));
    }

    edges.clear();
    edges.resize(segmentCount * 2);

    strand0Indices.clear();
    strand0Indices.reserve(anchors.journeysWithPositions.size() / 2);
    for(uint64_t oridValue = 0; oridValue < anchors.journeysWithPositions.size(); oridValue += 2) {
        if(anchors.journeysWithPositions.size(oridValue) >= 2) {
            strand0Indices.push_back(oridValue);
        }
    }

    buildDirectedAnchorsPtr = &anchors;
    buildReadIntervalsPtr = nullptr;
    threadPaths.resize(threadCount);
    threadPathReadSupport.resize(threadCount);
    threadEdgePairs.resize(threadCount);
    threadEdgeOverlapSamples.resize(threadCount);
    for(uint64_t t = 0; t < threadCount; ++t) {
        threadPaths[t].clear();
        threadPathReadSupport[t].clear();
        threadEdgePairs[t].clear();
        threadEdgeOverlapSamples[t].clear();
    }

    setupLoadBalancing(strand0Indices.size(), 1000);
    runThreads(&DirectedAnchorGraph::buildFromDirectedAnchorsThreadFunction, threadCount);

    // Merge logic is identical to buildFromAnchors Phase 3.
    // Call a helper or just reimplement? Let's just finish Phase 4 logic here.
    // (Actually many parts are shared, but for a one-off refactor this works).
    
    // Phase 3: Merge per-thread results.
    map<vector<DagNodeId>, vector<DagPathReadSupport>> groupedPaths;
    for(uint64_t t = 0; t < threadCount; ++t) {
        for(size_t i = 0; i < threadPaths[t].size(); i++) {
            groupedPaths[threadPaths[t][i]].push_back(threadPathReadSupport[t][i]);
        }
        threadPaths[t].clear();
        threadPathReadSupport[t].clear();
    }

    paths.clear();
    pathRemoved.clear();
    pathWeights.clear();
    pathReadIds.clear();
    for(const auto& [path, support] : groupedPaths) {
        paths.push_back(path);
        pathRemoved.push_back(false);
        pathWeights.push_back(support.size());
        pathReadIds.push_back(support);
    }
    for(uint64_t pathIdx = 0; pathIdx < paths.size(); pathIdx++) {
        refreshReadPosesForPath(pathIdx);
    }

    vector<unordered_set<DagNodeId>> tempEdgeSets(edges.size());
    unordered_map<pair<DagNodeId, DagNodeId>, vector<uint64_t>, DagEdgePairHash> edgeOverlapSamples;

    for(uint64_t t = 0; t < threadCount; ++t) {
        for(const auto& [from, to] : threadEdgePairs[t]) {
            tempEdgeSets[from].insert(to);
            tempEdgeSets[to ^ 1].insert(from ^ 1);
        }
        for(const auto& sample : threadEdgeOverlapSamples[t]) {
            edgeOverlapSamples[canonEdge(get<0>(sample), get<1>(sample))].push_back(get<2>(sample));
        }
        threadEdgePairs[t].clear();
        threadEdgeOverlapSamples[t].clear();
    }

    for(uint64_t i = 0; i < tempEdgeSets.size(); ++i) {
        if(!tempEdgeSets[i].empty()) {
            edges[i].assign(tempEdgeSets[i].begin(), tempEdgeSets[i].end());
            sort(edges[i].begin(), edges[i].end());
            for(DagNodeId to : edges[i]) {
                const auto canonical = canonEdge(i, to);
                if(edgeOverlaps.find(canonical) != edgeOverlaps.end()) continue;
                uint64_t overlapBp = 0;
                auto it = edgeOverlapSamples.find(canonical);
                if(it != edgeOverlapSamples.end() && !it->second.empty()) {
                    overlapBp = medianValue(it->second);
                }
                const uint64_t fromLen = nodeBpLength(segmentOf(canonical.first));
                const uint64_t toLen = nodeBpLength(segmentOf(canonical.second));
                const uint64_t maxOverlap = (fromLen > 0 && toLen > 0) ? (min(fromLen, toLen) - 1) : 0;
                edgeOverlaps.emplace(canonical, min(overlapBp, maxOverlap));
            }
        }
    }

    pathsCrossing.rebuild(paths, pathRemoved);
    buildDirectedAnchorsPtr = nullptr;
    performanceLog << timestamp << "DirectedAnchorGraph::buildFromAnchors (DirectedAnchors) ends." << endl;
}

// ============================================================================
// buildFromDirectedAnchorsThreadFunction
//
// Processes oriented read journeys in parallel. Extracts edges and path support.
// Overlaps are calculated directly from base-pair coordinates:
// Overlap = (Node1End - Node2Start) if Node2Start < Node1End, else 0.
// ============================================================================
void DirectedAnchorGraph::buildFromDirectedAnchorsThreadFunction(size_t threadId)
{
    const DirectedAnchors& anchors = *buildDirectedAnchorsPtr;
    auto& myPaths = threadPaths[threadId];
    auto& myReadSupport = threadPathReadSupport[threadId];
    auto& myEdges = threadEdgePairs[threadId];
    auto& myOverlapSamples = threadEdgeOverlapSamples[threadId];

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t idx = begin; idx < end; ++idx) {
            const OrientedReadId orid0 = strand0Indices[idx];
            const OrientedReadId orid1 = orid0.getOppositeStrand();
            
            auto processOneStrand = [&](const OrientedReadId orid) {
                const auto journey = anchors.journeysWithPositions[orid.getValue()];
                if(journey.empty()) return;

                vector<DagNodeId> path;
                path.reserve(journey.size());
                for(size_t i = 0; i < journey.size(); i++) {
                    path.push_back(DagNodeId(journey[i].anchorId));
                }

                for(size_t i = 0; i + 1 < journey.size(); i++) {
                    myEdges.push_back({path[i], path[i+1]});
                    const uint32_t end0 = journey[i].end;
                    const uint32_t start1 = journey[i+1].start;
                    uint64_t overlap = (start1 < end0) ? (end0 - start1) : 0;
                    myOverlapSamples.emplace_back(path[i], path[i+1], overlap);
                }

                myPaths.push_back(std::move(path));
                DagPathReadSupport support;
                support.readId = orid.getValue() / 2;
                support.readNameIndex = support.readId;
                support.readInfoIndex = support.readId;
                support.readPosZeroOffset = journey.front().start;
                const uint64_t pathBp = journey.back().end - journey.front().start;
                support.readPosStartIndex = 0;
                support.readPosEndIndex = pathBp;
                support.readPoses.resize(journey.size() + 1);
                for(size_t i = 0; i < journey.size(); i++) {
                    support.readPoses[i] = journey[i].start - support.readPosZeroOffset;
                }
                support.readPoses.back() = journey.back().end - support.readPosZeroOffset;
                myReadSupport.push_back(std::move(support));
            };

            processOneStrand(orid0);
            processOneStrand(orid1);
        }
    }
}
