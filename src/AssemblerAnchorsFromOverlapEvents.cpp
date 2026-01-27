// Dinara.
#include "Assembler.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
#include "mode3-Anchor.hpp"

// Standard library.
#include "algorithm.hpp"
#include <thread>
#include <unordered_map>

using namespace dinara;
using namespace std;

namespace {
    inline MarkerGraph::VertexId asVertexId(MarkerGraph::CompressedVertexId x)
    {
        return MarkerGraph::VertexId(uint64_t(x));
    }
}



shared_ptr<mode3::Anchors> Assembler::createAnchorsFromMarkerGraphVerticesAtOverlapEvents(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount)
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkReadGraphIsOpen();
    if(!alignmentData.isOpen) {
        throw runtime_error("Alignment data are not accessible.");
    }
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    const auto& mgVertices = markerGraph.vertices();

    vector<vector<MarkerGraphVertexId>> threadSelected(threadCount);

    uint64_t chunk = orientedReadCount / threadCount;
    if(chunk == 0) chunk = 1;

    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = (t == threadCount - 1) ? orientedReadCount : min(orientedReadCount, (t+1) * chunk);
            auto& out = threadSelected[t];
            out.reserve((end - begin) * 4);

            vector<uint32_t> eventOrdinals;
            eventOrdinals.reserve(256);

            for(uint64_t v=begin; v<end; ++v) {
                const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
                const uint32_t markerCount = uint32_t(markers->size(orientedReadId.getValue()));
                if(markerCount == 0) {
                    continue;
                }

                eventOrdinals.clear();

                // Collect start/end events from read-graph overlaps incident to this oriented read.
                for(const uint32_t edgeId : readGraph.connectivity[orientedReadId.getValue()]) {
                    const ReadGraphEdge& edge = readGraph.edges[edgeId];
                    if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                        continue;
                    }
                    const uint64_t alignmentId = edge.alignmentId;
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if(ad.isDeleted()) {
                        continue;
                    }
                    if(!ad.info.isInReadGraph) {
                        continue;
                    }
                    const OrientedReadId other = edge.getOther(orientedReadId);
                    const AlignmentInfo info = ad.orient(orientedReadId, other);

                    const uint32_t first = info.data[0].firstOrdinal;
                    const uint32_t last = info.data[0].lastOrdinal;
                    if(first < markerCount) {
                        eventOrdinals.push_back(first);
                    }
                    const uint32_t afterLast = last + 1;
                    if(afterLast < markerCount) {
                        eventOrdinals.push_back(afterLast);
                    }
                }

                if(eventOrdinals.empty()) {
                    continue;
                }

                std::sort(eventOrdinals.begin(), eventOrdinals.end());
                eventOrdinals.erase(std::unique(eventOrdinals.begin(), eventOrdinals.end()), eventOrdinals.end());

                MarkerGraphVertexId lastSelected = MarkerGraph::invalidVertexId;
                for(const uint32_t ordinal : eventOrdinals) {
                    if(ordinal >= markerCount) {
                        continue;
                    }
                    const MarkerId markerId = getMarkerId(orientedReadId, ordinal);
                    const auto compressedVertex = markerGraph.vertexTable[markerId];
                    if(compressedVertex == MarkerGraph::invalidCompressedVertexId) {
                        continue;
                    }
                    const MarkerGraphVertexId vertexId = asVertexId(compressedVertex);
                    if(vertexId >= mgVertices.size()) {
                        continue;
                    }
                    const uint64_t cov = mgVertices.size(vertexId);
                    if(cov < minAnchorCoverage || cov > maxAnchorCoverage) {
                        continue;
                    }

                    const MarkerGraphVertexId rc = markerGraph.reverseComplementVertex[vertexId];
                    const MarkerGraphVertexId canonical = min(vertexId, rc);
                    if(canonical == lastSelected) {
                        continue;
                    }
                    out.push_back(canonical);
                    lastSelected = canonical;
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    // Merge and deduplicate selected canonical vertices.
    vector<MarkerGraphVertexId> selected;
    {
        size_t total = 0;
        for(const auto& v : threadSelected) total += v.size();
        selected.reserve(total);
        for(auto& v : threadSelected) {
            selected.insert(selected.end(), v.begin(), v.end());
        }
    }
    std::sort(selected.begin(), selected.end());
    selected.erase(std::unique(selected.begin(), selected.end()), selected.end());

    cout << timestamp << "Selected " << selected.size()
         << " marker graph vertices from overlap events for anchors." << endl;

    return make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        markerGraph,
        selected,
        minAnchorCoverage,
        maxAnchorCoverage,
        threadCount);
}



shared_ptr<mode3::Anchors> Assembler::createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount)
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkReadGraphIsOpen();
    if(!alignmentData.isOpen) {
        throw runtime_error("Alignment data are not accessible.");
    }
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    const auto& mgVertices = markerGraph.vertices();

    struct Event {
        uint32_t ordinal = 0;
        int32_t delta = 0; // +1 start, -1 end (end = lastOrdinal+1)
    };

    vector<vector<MarkerGraphVertexId>> threadSelected(threadCount);

    uint64_t chunk = orientedReadCount / threadCount;
    if(chunk == 0) chunk = 1;

    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = (t == threadCount - 1) ? orientedReadCount : min(orientedReadCount, (t+1) * chunk);
            auto& out = threadSelected[t];
            out.reserve((end - begin) * 2);

            vector<Event> events;
            events.reserve(512);

            vector<Event> merged;
            merged.reserve(512);

            std::unordered_map<MarkerGraphVertexId, bool> duplicateReadIdCache;
            duplicateReadIdCache.reserve(4096);

            for(uint64_t v=begin; v<end; ++v) {
                const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
                const uint32_t markerCount = uint32_t(markers->size(orientedReadId.getValue()));
                if(markerCount == 0) {
                    continue;
                }

                events.clear();

                // Collect start/end events from read-graph overlaps incident to this oriented read.
                for(const uint32_t edgeId : readGraph.connectivity[orientedReadId.getValue()]) {
                    const ReadGraphEdge& edge = readGraph.edges[edgeId];
                    if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                        continue;
                    }
                    const uint64_t alignmentId = edge.alignmentId;
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if(ad.isDeleted()) {
                        continue;
                    }
                    if(!ad.info.isInReadGraph) {
                        continue;
                    }

                    const OrientedReadId other = edge.getOther(orientedReadId);
                    const AlignmentInfo info = ad.orient(orientedReadId, other);

                    const uint32_t first = info.data[0].firstOrdinal;
                    const uint32_t last = info.data[0].lastOrdinal;
                    if(first < markerCount) {
                        events.push_back({first, +1});
                    }
                    const uint32_t afterLast = last + 1;
                    if(afterLast < markerCount) {
                        events.push_back({afterLast, -1});
                    }
                }

                if(events.empty()) {
                    continue;
                }

                std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
                    return a.ordinal < b.ordinal;
                });

                merged.clear();
                for(const auto& e : events) {
                    if(merged.empty() || merged.back().ordinal != e.ordinal) {
                        merged.push_back(e);
                    } else {
                        merged.back().delta += e.delta;
                    }
                }

                int32_t active = 0;
                MarkerGraphVertexId lastSelected = MarkerGraph::invalidVertexId;

                for(size_t i=0; i<merged.size(); ++i) {
                    active += merged[i].delta;

                    const uint32_t segmentStart = merged[i].ordinal;
                    const uint32_t segmentEnd = (i + 1 < merged.size()) ? merged[i + 1].ordinal : markerCount;
                    if(active <= 0) {
                        continue;
                    }
                    if(segmentEnd <= segmentStart) {
                        continue;
                    }

                    // Scan all ordinals in this interval to pick the vertex with maximum coverage.
                    // Canonicalize by RC to avoid duplicates.
                    const uint32_t mid = segmentStart + (segmentEnd - segmentStart) / 2;
                    // Prefer a "clean" vertex with no duplicate ReadIds (does not contain both strands of any read).
                    // If none exist in this interval, fall back to the best vertex regardless.
                    uint64_t bestAnyCov = 0;
                    MarkerGraphVertexId bestAny = MarkerGraph::invalidVertexId;
                    uint32_t bestAnyDistanceToMid = std::numeric_limits<uint32_t>::max();

                    uint64_t bestCleanCov = 0;
                    MarkerGraphVertexId bestClean = MarkerGraph::invalidVertexId;
                    uint32_t bestCleanDistanceToMid = std::numeric_limits<uint32_t>::max();

                    MarkerGraphVertexId previous = MarkerGraph::invalidVertexId;
                    for(uint32_t ordinal = segmentStart; ordinal < segmentEnd; ++ordinal) {
                        const MarkerId markerId = getMarkerId(orientedReadId, ordinal);
                        const auto compressedVertex = markerGraph.vertexTable[markerId];
                        if(compressedVertex == MarkerGraph::invalidCompressedVertexId) {
                            continue;
                        }
                        const MarkerGraphVertexId vertexId = asVertexId(compressedVertex);
                        if(vertexId >= mgVertices.size()) {
                            continue;
                        }
                        const uint64_t cov = mgVertices.size(vertexId);
                        if(cov < minAnchorCoverage || cov > maxAnchorCoverage) {
                            continue;
                        }

                        const MarkerGraphVertexId rc = markerGraph.reverseComplementVertex[vertexId];
                        const MarkerGraphVertexId canonical = min(vertexId, rc);
                        if(canonical == previous) {
                            continue;
                        }
                        previous = canonical;

                        const uint32_t distanceToMid = (ordinal > mid) ? (ordinal - mid) : (mid - ordinal);
                        // Update best-any (no cleanliness constraint).
                        if(cov > bestAnyCov ||
                           (cov == bestAnyCov && distanceToMid < bestAnyDistanceToMid) ||
                           (cov == bestAnyCov && distanceToMid == bestAnyDistanceToMid && canonical < bestAny)) {
                            bestAnyCov = cov;
                            bestAny = canonical;
                            bestAnyDistanceToMid = distanceToMid;
                        }

                        // Update best-clean only if this vertex has no duplicate ReadIds.
                        if(cov > bestCleanCov ||
                           (cov == bestCleanCov && distanceToMid < bestCleanDistanceToMid) ||
                           (cov == bestCleanCov && distanceToMid == bestCleanDistanceToMid && canonical < bestClean)) {
                            auto dupIt = duplicateReadIdCache.find(canonical);
                            bool hasDuplicateReadIds = false;
                            if(dupIt != duplicateReadIdCache.end()) {
                                hasDuplicateReadIds = dupIt->second;
                            } else {
                                hasDuplicateReadIds = markerGraph.vertexHasDuplicateReadIds(canonical, *markers);
                                duplicateReadIdCache.insert({canonical, hasDuplicateReadIds});
                            }
                            if(!hasDuplicateReadIds) {
                                bestCleanCov = cov;
                                bestClean = canonical;
                                bestCleanDistanceToMid = distanceToMid;
                                // Only safe to early-exit if we found a clean vertex at maximum allowed coverage.
                                if(bestCleanCov == maxAnchorCoverage) {
                                    break;
                                }
                            }
                        }
                    }

                    const MarkerGraphVertexId chosen =
                        (bestClean != MarkerGraph::invalidVertexId) ? bestClean : bestAny;
                    if(chosen == MarkerGraph::invalidVertexId) {
                        continue;
                    }
                    if(chosen == lastSelected) {
                        continue;
                    }
                    out.push_back(chosen);
                    lastSelected = chosen;
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    // Merge and deduplicate selected canonical vertices.
    vector<MarkerGraphVertexId> selected;
    {
        size_t total = 0;
        for(const auto& v : threadSelected) total += v.size();
        selected.reserve(total);
        for(auto& v : threadSelected) {
            selected.insert(selected.end(), v.begin(), v.end());
        }
    }
    std::sort(selected.begin(), selected.end());
    selected.erase(std::unique(selected.begin(), selected.end()), selected.end());

    cout << timestamp << "Selected " << selected.size()
         << " marker graph vertices (best per overlap interval) for anchors." << endl;

    return make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        markerGraph,
        selected,
        minAnchorCoverage,
        maxAnchorCoverage,
        threadCount);
}
