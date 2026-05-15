// Shasta2AnchorsFromSplitVertices.cpp
//
// Creates Shasta2Anchors by splitting marker graph vertices whose reads
// were merged by transitive closure but lack direct pairwise overlaps.
// Each vertex is partitioned into overlap-connected components using the
// read graph, producing anchors that only contain mutually overlapping reads.
//
// Two-pass multithreaded approach:
//   Pass 1: Split each canonical vertex, record per-anchor read counts.
//           Prefix sum assigns deterministic anchor IDs per vertex.
//   Pass 2: Re-split each canonical vertex, fill anchorMarkerInfos.
//           Each thread writes to non-overlapping anchor ID ranges.

#include "Shasta2AnchorsFromSplitVertices.hpp"
#include "DINARA_ASSERT.hpp"
#include "MarkerGraph.hpp"
#include "ReadGraph.hpp"
#include "Reads.hpp"
#include "findMarkerId.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <atomic>
#include <iostream>
#include <thread>
#include <vector>

using namespace dinara;
using namespace std;



// Split a single marker graph vertex into overlap-connected components.
// If the vertex is coherent (all reads overlap-connected), returns a single group.
static void splitVertex(
    const MarkerGraph& markerGraph,
    MarkerGraphVertexId vertexId,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const ReadGraph& readGraph,
    const Reads& reads,
    vector<vector<Shasta2AnchorMarkerInfo>>& groups)
{
    groups.clear();

    struct ReadEntry {
        OrientedReadId orientedReadId;
        uint32_t ordinal;
        uint64_t rawLength;
    };
    vector<ReadEntry> entries;
    for(const MarkerId markerId : markerGraph.getVertexMarkerIds(vertexId)) {
        OrientedReadId orientedReadId;
        uint32_t ordinal;
        tie(orientedReadId, ordinal) = findMarkerId(markerId, markers);
        entries.push_back({orientedReadId, ordinal,
            reads.getReadRawSequenceLength(orientedReadId.getReadId())});
    }

    // Sort by OrientedReadId for deduplication.
    sort(entries.begin(), entries.end(),
        [](const ReadEntry& a, const ReadEntry& b) {
            return a.orientedReadId < b.orientedReadId;
        });
    {
        auto last = unique(entries.begin(), entries.end(),
            [](const ReadEntry& a, const ReadEntry& b) {
                return a.orientedReadId == b.orientedReadId;
            });
        entries.erase(last, entries.end());
    }

    // Sort by raw read length descending, then by OrientedReadId for stability.
    sort(entries.begin(), entries.end(),
        [](const ReadEntry& a, const ReadEntry& b) {
            if(a.rawLength != b.rawLength) return a.rawLength > b.rawLength;
            return a.orientedReadId < b.orientedReadId;
        });

    if(entries.size() <= 1) {
        if(!entries.empty()) {
            groups.resize(1);
            groups[0].emplace_back(entries[0].orientedReadId, entries[0].ordinal);
        }
        return;
    }

    const uint64_t n = entries.size();
    vector<uint32_t> vertexReadsSorted(n);
    for(uint64_t i = 0; i < n; i++) {
        vertexReadsSorted[i] = entries[i].orientedReadId.getValue();
    }
    sort(vertexReadsSorted.begin(), vertexReadsSorted.end());

    const bool readGraphIsOpen = readGraph.connectivity.isOpen();
    vector<vector<uint32_t>> inVertexNeighbors(n);
    for(uint64_t i = 0; i < n; i++) {
        if(!readGraphIsOpen) continue;
        const auto orientedReadId = entries[i].orientedReadId;
        const auto adjacentEdges = readGraph.connectivity[orientedReadId.getValue()];
        for(const uint32_t edgeId : adjacentEdges) {
            const ReadGraphEdge& edge = readGraph.edges[edgeId];
            if(edge.crossesStrands) continue;
            const OrientedReadId neighbor = edge.getOther(orientedReadId);
            const uint32_t neighborVal = neighbor.getValue();
            if(binary_search(vertexReadsSorted.begin(), vertexReadsSorted.end(), neighborVal)) {
                inVertexNeighbors[i].push_back(neighborVal);
            }
        }
        sort(inVertexNeighbors[i].begin(), inVertexNeighbors[i].end());
    }

    // Find roots.
    vector<int64_t> groupAssignment(n, -1);
    vector<uint64_t> rootIndices;
    rootIndices.push_back(0);
    groupAssignment[0] = 0;

    for(uint64_t i = 1; i < n; i++) {
        bool overlapsAnyRoot = false;
        const auto& neighbors = inVertexNeighbors[i];
        for(const uint64_t rootIdx : rootIndices) {
            const uint32_t rootVal = entries[rootIdx].orientedReadId.getValue();
            if(binary_search(neighbors.begin(), neighbors.end(), rootVal)) {
                overlapsAnyRoot = true;
                break;
            }
        }
        if(!overlapsAnyRoot) {
            const int64_t groupIdx = static_cast<int64_t>(rootIndices.size());
            rootIndices.push_back(i);
            groupAssignment[i] = groupIdx;
        }
    }

    if(rootIndices.size() <= 1) {
        groups.resize(1);
        groups[0].reserve(n);
        for(const auto& entry : entries) {
            groups[0].emplace_back(entry.orientedReadId, entry.ordinal);
        }
        sort(groups[0].begin(), groups[0].end());
        return;
    }

    // Assign remaining reads to roots.
    for(uint64_t i = 0; i < n; i++) {
        if(groupAssignment[i] >= 0) continue;
        const auto& neighbors = inVertexNeighbors[i];
        int64_t bestGroup = -1;
        uint64_t bestLength = 0;
        for(uint64_t r = 0; r < rootIndices.size(); r++) {
            const uint64_t rootIdx = rootIndices[r];
            const uint32_t rootVal = entries[rootIdx].orientedReadId.getValue();
            if(binary_search(neighbors.begin(), neighbors.end(), rootVal)) {
                const uint64_t rootLen = entries[rootIdx].rawLength;
                if(rootLen > bestLength || (rootLen == bestLength && static_cast<int64_t>(r) < bestGroup)) {
                    bestLength = rootLen;
                    bestGroup = static_cast<int64_t>(r);
                }
            }
        }
        DINARA_ASSERT(bestGroup >= 0);
        groupAssignment[i] = bestGroup;
    }

    // Build output groups.
    const uint64_t numGroups = rootIndices.size();
    groups.resize(numGroups);
    for(uint64_t i = 0; i < n; i++) {
        const int64_t gIdx = groupAssignment[i];
        DINARA_ASSERT(gIdx >= 0);
        groups[gIdx].emplace_back(entries[i].orientedReadId, entries[i].ordinal);
    }
    for(auto& group : groups) {
        sort(group.begin(), group.end());
    }
    groups.erase(
        remove_if(groups.begin(), groups.end(),
            [](const vector<Shasta2AnchorMarkerInfo>& g) { return g.empty(); }),
        groups.end());
}



shared_ptr<Shasta2Anchors> dinara::createShasta2AnchorsFromSplitVertices(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    const ReadGraph& readGraph,
    uint64_t threadCount,
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage)
{
    cout << timestamp << "Creating Shasta2Anchors with vertex splitting." << endl;

    auto anchors = make_shared<Shasta2Anchors>(
        mappedMemoryOwner,
        reads,
        k,
        markers,
        markerGraph,
        threadCount);

    cout << timestamp << "Splitting anchors by overlap components using "
         << threadCount << " threads." << endl;

    const uint64_t vertexCount = markerGraph.vertexCount();

    // Identify canonical vertices.
    vector<MarkerGraphVertexId> canonicalVertices;
    canonicalVertices.reserve(vertexCount / 2);
    for(MarkerGraphVertexId v = 0; v < vertexCount; ++v) {
        if(v <= markerGraph.reverseComplementVertex[v]) {
            canonicalVertices.push_back(v);
        }
    }
    const uint64_t canonicalCount = canonicalVertices.size();

    // ---------------------------------------------------------------
    // Pass 1: Split each canonical vertex, record per-anchor sizes.
    // Each canonical vertex produces a contiguous block of anchors:
    //   [canonical_group_0, rc_group_0, canonical_group_1, rc_group_1, ...]
    // For self-complement vertices, no RC anchors are emitted.
    // ---------------------------------------------------------------

    // Per canonical vertex: number of anchors it produces.
    vector<uint64_t> anchorsPerVertex(canonicalCount, 0);
    // Per canonical vertex: read counts for each anchor it produces.
    // Stored as a flat vector of vectors for cache locality.
    vector<vector<uint64_t>> anchorSizesPerVertex(canonicalCount);

    struct ThreadCounters {
        uint64_t verticesSplit = 0;
        uint64_t anchorsUnsplit = 0;
        uint64_t anchorsFromSplits = 0;
        uint64_t anchorsDiscarded = 0;
    };
    vector<ThreadCounters> threadCounters(threadCount);

    atomic<uint64_t> nextBatch(0);
    const uint64_t batchSize = 1000;

    auto pass1Function = [&](uint64_t threadId) {
        auto& tc = threadCounters[threadId];
        vector<vector<Shasta2AnchorMarkerInfo>> groups;

        while(true) {
            const uint64_t batchBegin = nextBatch.fetch_add(batchSize);
            if(batchBegin >= canonicalCount) break;
            const uint64_t batchEnd = min(batchBegin + batchSize, canonicalCount);

            for(uint64_t ci = batchBegin; ci < batchEnd; ci++) {
                const MarkerGraphVertexId vertexId = canonicalVertices[ci];
                const MarkerGraphVertexId rcVertexId =
                    markerGraph.reverseComplementVertex[vertexId];
                const bool isSelfComplement = (vertexId == rcVertexId);

                splitVertex(markerGraph, vertexId, markers, readGraph, reads, groups);

                if(groups.empty()) continue;

                const bool wasSplit = (groups.size() > 1);
                if(wasSplit) ++tc.verticesSplit;

                auto& sizes = anchorSizesPerVertex[ci];
                uint64_t count = 0;

                for(const auto& group : groups) {
                    const uint64_t groupSize = group.size();
                    if(groupSize < minAnchorCoverage || groupSize > maxAnchorCoverage) {
                        ++tc.anchorsDiscarded;
                        if(!isSelfComplement) ++tc.anchorsDiscarded;
                        continue;
                    }
                    // Canonical anchor.
                    sizes.push_back(groupSize);
                    ++count;
                    if(wasSplit) ++tc.anchorsFromSplits;
                    else ++tc.anchorsUnsplit;

                    // RC anchor (same read count).
                    if(!isSelfComplement) {
                        sizes.push_back(groupSize);
                        ++count;
                        if(wasSplit) ++tc.anchorsFromSplits;
                        else ++tc.anchorsUnsplit;
                    }
                }
                anchorsPerVertex[ci] = count;
            }
        }
    };

    {
        vector<thread> threads;
        threads.reserve(threadCount);
        for(uint64_t t = 0; t < threadCount; t++) {
            threads.emplace_back(pass1Function, t);
        }
        for(auto& t : threads) t.join();
    }

    // Prefix sum to get anchor ID offsets per canonical vertex.
    vector<uint64_t> anchorIdOffset(canonicalCount + 1, 0);
    for(uint64_t ci = 0; ci < canonicalCount; ci++) {
        anchorIdOffset[ci + 1] = anchorIdOffset[ci] + anchorsPerVertex[ci];
    }
    const uint64_t totalAnchorCount = anchorIdOffset[canonicalCount];

    // Merge counters and report.
    uint64_t totalVerticesSplit = 0;
    uint64_t totalAnchorsUnsplit = 0;
    uint64_t totalAnchorsFromSplits = 0;
    uint64_t totalAnchorsDiscarded = 0;
    for(const auto& tc : threadCounters) {
        totalVerticesSplit += tc.verticesSplit;
        totalAnchorsUnsplit += tc.anchorsUnsplit;
        totalAnchorsFromSplits += tc.anchorsFromSplits;
        totalAnchorsDiscarded += tc.anchorsDiscarded;
    }

    cout << timestamp << "Processed " << canonicalCount << " canonical vertices." << endl;
    cout << timestamp << "Unsplit anchors: " << totalAnchorsUnsplit << endl;
    cout << timestamp << "Split " << totalVerticesSplit << " vertices into "
         << totalAnchorsFromSplits << " anchors." << endl;
    cout << timestamp << "Discarded " << totalAnchorsDiscarded
         << " groups outside coverage range ["
         << minAnchorCoverage << ", " << maxAnchorCoverage << "]." << endl;
    cout << timestamp << "Total anchors: " << totalAnchorCount << endl;

    // Build VectorOfVectors structure using sizes from Pass 1.
    anchors->anchorMarkerInfos.createNew(
        anchors->largeDataName("Shasta2Anchors-AnchorMarkerInfos"),
        anchors->largeDataPageSize);

    anchors->anchorMarkerInfos.beginPass1(totalAnchorCount);
    for(uint64_t ci = 0; ci < canonicalCount; ci++) {
        const auto& sizes = anchorSizesPerVertex[ci];
        uint64_t anchorId = anchorIdOffset[ci];
        for(const uint64_t sz : sizes) {
            anchors->anchorMarkerInfos.incrementCount(anchorId, sz);
            ++anchorId;
        }
    }
    anchors->anchorMarkerInfos.beginPass2();

    // Free pass 1 sizing data.
    anchorSizesPerVertex.clear();
    anchorSizesPerVertex.shrink_to_fit();
    anchorsPerVertex.clear();
    anchorsPerVertex.shrink_to_fit();

    // ---------------------------------------------------------------
    // Pass 2: Re-split each canonical vertex, fill anchorMarkerInfos.
    // Each thread writes to non-overlapping anchor ID ranges.
    // ---------------------------------------------------------------

    nextBatch.store(0);

    auto pass2Function = [&](uint64_t /*threadId*/) {
        vector<vector<Shasta2AnchorMarkerInfo>> groups;
        vector<Shasta2AnchorMarkerInfo> rcGroup;

        while(true) {
            const uint64_t batchBegin = nextBatch.fetch_add(batchSize);
            if(batchBegin >= canonicalCount) break;
            const uint64_t batchEnd = min(batchBegin + batchSize, canonicalCount);

            for(uint64_t ci = batchBegin; ci < batchEnd; ci++) {
                if(anchorIdOffset[ci] == anchorIdOffset[ci + 1]) continue;

                const MarkerGraphVertexId vertexId = canonicalVertices[ci];
                const MarkerGraphVertexId rcVertexId =
                    markerGraph.reverseComplementVertex[vertexId];
                const bool isSelfComplement = (vertexId == rcVertexId);

                splitVertex(markerGraph, vertexId, markers, readGraph, reads, groups);

                uint64_t anchorId = anchorIdOffset[ci];
                for(const auto& group : groups) {
                    if(group.size() < minAnchorCoverage || group.size() > maxAnchorCoverage) {
                        continue;
                    }

                    // Store canonical anchor (reverse order for store()).
                    for(auto it = group.rbegin(); it != group.rend(); ++it) {
                        anchors->anchorMarkerInfos.store(anchorId, *it);
                    }
                    ++anchorId;

                    // Store RC anchor.
                    if(!isSelfComplement) {
                        rcGroup.clear();
                        rcGroup.reserve(group.size());
                        for(const auto& info : group) {
                            OrientedReadId rcOrientedReadId = info.orientedReadId;
                            rcOrientedReadId.flipStrand();
                            const uint32_t rcOrdinal = uint32_t(
                                markers[info.orientedReadId.getValue()].size() - 1 - info.ordinal);
                            rcGroup.emplace_back(rcOrientedReadId, rcOrdinal);
                        }
                        sort(rcGroup.begin(), rcGroup.end());
                        for(auto it = rcGroup.rbegin(); it != rcGroup.rend(); ++it) {
                            anchors->anchorMarkerInfos.store(anchorId, *it);
                        }
                        ++anchorId;
                    }
                }
            }
        }
    };

    {
        vector<thread> threads;
        threads.reserve(threadCount);
        for(uint64_t t = 0; t < threadCount; t++) {
            threads.emplace_back(pass2Function, t);
        }
        for(auto& t : threads) t.join();
    }

    anchors->anchorMarkerInfos.endPass2();

    // Free offset data.
    anchorIdOffset.clear();
    anchorIdOffset.shrink_to_fit();
    canonicalVertices.clear();
    canonicalVertices.shrink_to_fit();

    cout << timestamp << "Shasta2Anchors with vertex splitting completed. "
         << totalAnchorCount << " anchors." << endl;

    return anchors;
}
