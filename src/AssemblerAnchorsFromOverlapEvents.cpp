// Dinara.
#include "Assembler.hpp"
#include "Kmer.hpp"
#include "KmerCounter.hpp"
#include "Reads.hpp"
#include "extractKmer.hpp"
#include "findMarkerId.hpp"
#include "timestamp.hpp"
#include "mode3-Anchor.hpp"

// Standard library.
#include "algorithm.hpp"
#include <functional>
#include <numeric>
#include <optional>
#include <thread>
#include <unordered_map>

using namespace dinara;
using namespace std;

namespace {
    inline MarkerGraph::VertexId asVertexId(MarkerGraph::CompressedVertexId x)
    {
        return MarkerGraph::VertexId(uint64_t(x));
    }

    using Interval = dinara::mode3::AnchorMarkerInterval;

    inline Interval reverseComplementInterval(
        Interval interval,
        const dinara::MemoryMapped::VectorOfVectors<dinara::CompressedMarker, uint64_t>& markers)
    {
        interval.orientedReadId.flipStrand();
        const uint64_t markerCount = markers.size(interval.orientedReadId.getValue());
        interval.ordinal0 = uint32_t(markerCount) - 1 - interval.ordinal0;
        return interval;
    }

    inline vector<Interval> reverseComplementAnchor(
        const vector<Interval>& anchor,
        const dinara::MemoryMapped::VectorOfVectors<dinara::CompressedMarker, uint64_t>& markers)
    {
        vector<Interval> rc = anchor;
        for(auto& interval : rc) {
            interval = reverseComplementInterval(interval, markers);
        }
        std::sort(rc.begin(), rc.end(), [](const Interval& a, const Interval& b) {
            return a.orientedReadId < b.orientedReadId;
        });
        return rc;
    }

    Kmer getMarkerKmer(
        const dinara::Reads& reads,
        const dinara::MemoryMapped::VectorOfVectors<dinara::CompressedMarker, uint64_t>& markers,
        uint64_t k,
        dinara::OrientedReadId orientedReadId,
        uint32_t ordinal)
    {
        const dinara::ReadId readId = orientedReadId.getReadId();
        const dinara::Strand strand = orientedReadId.getStrand();
        const dinara::LongBaseSequenceView readSequence = reads.getRead(readId);

        if(strand == 0) {
            const auto orientedMarkers = markers[orientedReadId.getValue()];
            const uint32_t position = orientedMarkers[ordinal].position;
            Kmer kmer;
            extractKmer(readSequence, position, k, kmer);
            return kmer;
        }

        // Strand 1: use the corresponding marker on strand 0 and reverse-complement.
        const dinara::OrientedReadId orientedReadId0(readId, 0);
        const auto orientedMarkers0 = markers[orientedReadId0.getValue()];
        const uint32_t markerCount0 = uint32_t(orientedMarkers0.size());
        const uint32_t ordinal0 = markerCount0 - 1 - ordinal;
        const uint32_t position0 = orientedMarkers0[ordinal0].position;
        Kmer kmer0;
        extractKmer(readSequence, position0, k, kmer0);
        return kmer0.reverseComplement(k);
    }

    bool mapMarkerOrdinalByOffsetAndKmer(
        const dinara::Reads& reads,
        const dinara::MemoryMapped::VectorOfVectors<dinara::CompressedMarker, uint64_t>& markers,
        uint64_t k,
        const dinara::AlignmentInfo& info,
        dinara::OrientedReadId orientedReadId0,
        uint32_t ordinal0,
        dinara::OrientedReadId orientedReadId1,
        uint32_t& ordinal1Out,
        const Kmer& seedKmer,
        uint32_t maxSearchRadius,
        uint32_t maxOffsetRange)
    {
        (void)orientedReadId0;

        if(ordinal0 < info.data[0].firstOrdinal || ordinal0 > info.data[0].lastOrdinal) {
            return false;
        }

        // ord0 - ord1 = offset  =>  ord1 = ord0 - offset
        int64_t lo = int64_t(ordinal0) - int64_t(info.maxOrdinalOffset);
        int64_t hi = int64_t(ordinal0) - int64_t(info.minOrdinalOffset);
        if(lo > hi) {
            std::swap(lo, hi);
        }

        // Clamp to aligned range on read1.
        lo = std::max<int64_t>(lo, int64_t(info.data[1].firstOrdinal));
        hi = std::min<int64_t>(hi, int64_t(info.data[1].lastOrdinal));
        if(lo > hi) {
            return false;
        }

        if(uint64_t(hi - lo) > maxOffsetRange) {
            return false;
        }

        const int64_t center = int64_t(ordinal0) - int64_t(info.averageOrdinalOffset);

        auto tryOrdinal1 = [&](int64_t candidate) -> bool {
            if(candidate < lo || candidate > hi) {
                return false;
            }
            const uint32_t o1 = uint32_t(candidate);
            const Kmer kmer1 = getMarkerKmer(reads, markers, k, orientedReadId1, o1);
            if(kmer1 == seedKmer) {
                ordinal1Out = o1;
                return true;
            }
            return false;
        };

        // Try center, then expand out.
        if(tryOrdinal1(center)) {
            return true;
        }
        for(uint32_t d=1; d<=maxSearchRadius; ++d) {
            if(tryOrdinal1(center - int64_t(d))) {
                return true;
            }
            if(tryOrdinal1(center + int64_t(d))) {
                return true;
            }
        }
        return false;
    }

    struct OverlapEvent {
        uint32_t ordinal = 0;
        int8_t delta = 0; // +1 start, -1 end (end = lastOrdinal+1)
        uint32_t edgeId = 0;
    };

    struct OverlapInterval {
        uint32_t start = 0;
        uint32_t end = 0; // one past last
    };

    inline bool findOverlapIntervalIndex(
        const vector<OverlapInterval>& intervals,
        uint32_t ordinal,
        uint32_t& indexOut)
    {
        if(intervals.empty()) {
            return false;
        }
        auto it = std::upper_bound(
            intervals.begin(),
            intervals.end(),
            ordinal,
            [](uint32_t value, const OverlapInterval& x) {
                return value < x.start;
            });
        if(it == intervals.begin()) {
            return false;
        }
        --it;
        if(ordinal >= it->end) {
            return false;
        }
        indexOut = uint32_t(it - intervals.begin());
        return true;
    }

    struct ArticulationResult {
        vector<uint8_t> isArticulation;
    };

    ArticulationResult findArticulationPoints(const vector<vector<uint32_t>>& adj)
    {
        const uint32_t n = uint32_t(adj.size());
        ArticulationResult result;
        result.isArticulation.assign(n, 0);
        if(n == 0) {
            return result;
        }

        vector<int32_t> disc(n, -1);
        vector<int32_t> low(n, -1);
        vector<int32_t> parent(n, -1);
        int32_t time = 0;

        std::function<void(uint32_t)> dfs = [&](uint32_t u) {
            disc[u] = low[u] = time++;
            uint32_t childCount = 0;
            for(const uint32_t v : adj[u]) {
                if(disc[v] == -1) {
                    parent[v] = int32_t(u);
                    ++childCount;
                    dfs(v);
                    low[u] = std::min(low[u], low[v]);

                    if(parent[u] == -1 && childCount > 1) {
                        result.isArticulation[u] = 1;
                    }
                    if(parent[u] != -1 && low[v] >= disc[u]) {
                        result.isArticulation[u] = 1;
                    }
                } else if(int32_t(v) != parent[u]) {
                    low[u] = std::min(low[u], disc[v]);
                }
            }
        };

        for(uint32_t i=0; i<n; ++i) {
            if(disc[i] == -1) {
                dfs(i);
            }
        }
        return result;
    }

    vector<vector<uint32_t>> connectedComponents(
        const vector<vector<uint32_t>>& adj,
        const vector<uint8_t>& removed)
    {
        const uint32_t n = uint32_t(adj.size());
        vector<uint8_t> visited(n, 0);
        vector<vector<uint32_t>> comps;
        comps.reserve(n);

        vector<uint32_t> stack;
        stack.reserve(n);

        for(uint32_t i=0; i<n; ++i) {
            if(removed[i] || visited[i]) {
                continue;
            }
            visited[i] = 1;
            stack.clear();
            stack.push_back(i);
            comps.emplace_back();
            auto& comp = comps.back();

            while(!stack.empty()) {
                const uint32_t u = stack.back();
                stack.pop_back();
                comp.push_back(u);
                for(const uint32_t v : adj[u]) {
                    if(removed[v] || visited[v]) {
                        continue;
                    }
                    visited[v] = 1;
                    stack.push_back(v);
                }
            }
        }
        return comps;
    }

    uint32_t countCommonSortedNeighbors(const vector<uint32_t>& a, const vector<uint32_t>& b)
    {
        uint32_t common = 0;
        size_t i = 0;
        size_t j = 0;
        while(i < a.size() && j < b.size()) {
            const uint32_t x = a[i];
            const uint32_t y = b[j];
            if(x == y) {
                ++common;
                ++i;
                ++j;
            } else if(x < y) {
                ++i;
            } else {
                ++j;
            }
        }
        return common;
    }

    vector<vector<uint32_t>> splitVertexByOverlapSupport(
        const vector<vector<uint32_t>>& adj,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        const uint32_t n = uint32_t(adj.size());
        if(n == 0) {
            return {};
        }
        if(n == 1) {
            if(minAnchorCoverage <= 1 && maxAnchorCoverage >= 1) {
                return { {0} };
            } else {
                return {};
            }
        }

        // Drop isolated nodes (no overlap support) if doing so keeps a valid anchor.
        vector<uint8_t> removedIsolated(n, 0);
        uint32_t keptAfterIsolated = 0;
        for(uint32_t i=0; i<n; ++i) {
            if(adj[i].empty()) {
                removedIsolated[i] = 1;
            } else {
                ++keptAfterIsolated;
            }
        }
        if(keptAfterIsolated < n &&
           keptAfterIsolated >= minAnchorCoverage && keptAfterIsolated <= maxAnchorCoverage) {
            return { [&]{
                vector<uint32_t> all;
                all.reserve(keptAfterIsolated);
                for(uint32_t i=0; i<n; ++i) {
                    if(!removedIsolated[i]) all.push_back(i);
                }
                return all;
            }() };
        }

        // If already disconnected, keep components.
        {
            vector<uint8_t> none(n, 0);
            const auto comps = connectedComponents(adj, none);
            if(comps.size() > 1) {
                vector<vector<uint32_t>> kept;
                for(const auto& comp : comps) {
                    if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                        kept.push_back(comp);
                    }
                }
                if(!kept.empty()) {
                    return kept;
                }
            }
        }

        // Remove articulation points to avoid "single-bridge" collapses.
        const auto art = findArticulationPoints(adj);
        vector<uint8_t> removed = art.isArticulation;
        const auto compsNoArt = connectedComponents(adj, removed);
        vector<vector<uint32_t>> keptNoArt;
        for(const auto& comp : compsNoArt) {
            if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                keptNoArt.push_back(comp);
            }
        }
        if(keptNoArt.size() >= 2) {
            return keptNoArt;
        }

        // If the vertex is held together by "weak" cross-edges with divergent neighbor sets,
        // split using only edges that have at least one common neighbor (triangle support).
        // This catches multi-bridge cases that are not single articulation points.
        {
            vector<vector<uint32_t>> strongAdj(n);
            for(uint32_t u=0; u<n; ++u) {
                strongAdj[u].reserve(adj[u].size());
            }
            for(uint32_t u=0; u<n; ++u) {
                for(const uint32_t v : adj[u]) {
                    if(v <= u) {
                        continue;
                    }
                    const uint32_t common = countCommonSortedNeighbors(adj[u], adj[v]);
                    if(common >= 1) {
                        strongAdj[u].push_back(v);
                        strongAdj[v].push_back(u);
                    }
                }
            }
            for(uint32_t u=0; u<n; ++u) {
                auto& nbr = strongAdj[u];
                std::sort(nbr.begin(), nbr.end());
                nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
            }

            vector<uint8_t> none(n, 0);
            const auto compsStrong = connectedComponents(strongAdj, none);
            vector<vector<uint32_t>> keptStrong;
            for(const auto& comp : compsStrong) {
                if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                    keptStrong.push_back(comp);
                }
            }
            if(keptStrong.size() >= 2) {
                return keptStrong;
            }
        }

        // Fall back to keeping all nodes if it yields a valid anchor.
        if(n >= minAnchorCoverage && n <= maxAnchorCoverage) {
            vector<uint32_t> all;
            all.reserve(n);
            for(uint32_t i=0; i<n; ++i) {
                all.push_back(i);
            }
            return {all};
        }

        // Otherwise, give up.
        return {};
    }

    // Like splitVertexByOverlapSupport, but uses phasing information when available:
    // - If there are at least two non-trivial connected components in the CIS-only overlap graph,
    //   we treat them as the core clusters.
    // - Reads without CIS edges (typically no het sites / unphased) are attached to exactly one
    //   core cluster based on the number of CIS edges they have to that cluster.
    // - This prevents unphased reads from "gluing" two haplotype clusters together via triangles.
    vector<vector<uint32_t>> splitVertexByOverlapSupportWithPhasing(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        if(!hasAnyCisEdge) {
            return splitVertexByOverlapSupport(adjAll, minAnchorCoverage, maxAnchorCoverage);
        }

        const uint32_t n = uint32_t(adjAll.size());
        vector<uint8_t> none(n, 0);
        const auto cisComps = connectedComponents(adjCis, none);

        vector<vector<uint32_t>> core;
        core.reserve(cisComps.size());
        for(const auto& comp : cisComps) {
            if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                core.push_back(comp);
            }
        }

        // If we don't have enough CIS structure to split, fall back to topology-only splitting.
        if(core.size() < 2) {
            return splitVertexByOverlapSupport(adjAll, minAnchorCoverage, maxAnchorCoverage);
        }

        vector<int32_t> assigned(n, -1);
        for(uint32_t ci=0; ci<uint32_t(core.size()); ++ci) {
            for(const uint32_t u : core[ci]) {
                assigned[u] = int32_t(ci);
            }
        }

        // Attach remaining nodes to exactly one core component using CIS support counts.
        vector<uint32_t> counts(core.size());
        for(uint32_t u=0; u<n; ++u) {
            if(assigned[u] != -1) {
                continue;
            }
            std::fill(counts.begin(), counts.end(), 0);
            for(const uint32_t v : adjCis[u]) {
                const int32_t ci = assigned[v];
                if(ci >= 0) {
                    ++counts[uint32_t(ci)];
                }
            }
            uint32_t bestCi = std::numeric_limits<uint32_t>::max();
            uint32_t bestCount = 0;
            uint32_t secondCount = 0;
            for(uint32_t ci=0; ci<counts.size(); ++ci) {
                const uint32_t c = counts[ci];
                if(c > bestCount) {
                    secondCount = bestCount;
                    bestCount = c;
                    bestCi = ci;
                } else if(c == bestCount && c != 0) {
                    // Tie.
                    secondCount = bestCount;
                } else if(c > secondCount) {
                    secondCount = c;
                }
            }
            if(bestCount == 0) {
                continue; // no CIS support -> leave unassigned
            }
            if(secondCount == bestCount) {
                continue; // ambiguous -> leave unassigned
            }
            if(core[bestCi].size() >= maxAnchorCoverage) {
                continue;
            }
            assigned[u] = int32_t(bestCi);
            core[bestCi].push_back(u);
        }

        // Return the split only if we still have at least two valid groups.
        vector<vector<uint32_t>> kept;
        kept.reserve(core.size());
        for(auto& comp : core) {
            if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                kept.push_back(comp);
            }
        }
        if(kept.size() >= 2) {
            return kept;
        }

        return splitVertexByOverlapSupport(adjAll, minAnchorCoverage, maxAnchorCoverage);
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



shared_ptr<mode3::Anchors> Assembler::createAnchorsFromMarkerGraphVerticesBestPerOverlapIntervalDecomposed(
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

    // Phase 1: select canonical marker graph vertices (same logic as best-per-interval).
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
                for(const uint32_t edgeId : readGraph.connectivity[orientedReadId.getValue()]) {
                    const ReadGraphEdge& edge = readGraph.edges[edgeId];
                    if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                        continue;
                    }
                    const uint64_t alignmentId = edge.alignmentId;
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if(!ad.keptByBothSides()) {
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

                    const uint32_t mid = segmentStart + (segmentEnd - segmentStart) / 2;
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
                        if(cov > bestAnyCov ||
                           (cov == bestAnyCov && distanceToMid < bestAnyDistanceToMid) ||
                           (cov == bestAnyCov && distanceToMid == bestAnyDistanceToMid && canonical < bestAny)) {
                            bestAnyCov = cov;
                            bestAny = canonical;
                            bestAnyDistanceToMid = distanceToMid;
                        }

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
         << " marker graph vertices (best per overlap interval) for anchor decomposition." << endl;

    // Phase 2: decompose each selected vertex using overlap support among its oriented reads.
    vector<vector<vector<Interval>>> threadAnchors(threadCount);
    uint64_t vChunk = selected.size() / threadCount;
    if(vChunk == 0) vChunk = 1;

    threads.clear();
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * vChunk;
            const uint64_t end = (t == threadCount - 1) ? selected.size() : min(selected.size(), (t + 1) * vChunk);
            auto& outAnchors = threadAnchors[t];
            outAnchors.reserve((end - begin) * 2);

            std::unordered_map<OrientedReadId::Int, uint32_t> index;
            index.reserve(256);

            vector<vector<uint32_t>> adjAll;
            adjAll.reserve(256);
            vector<vector<uint32_t>> adjCis;
            adjCis.reserve(256);

            vector<Interval> vertexIntervals;
            vertexIntervals.reserve(256);

            for(uint64_t i=begin; i<end; ++i) {
                const MarkerGraphVertexId vertexId = selected[i];
                const MarkerGraphVertexId rcVertexId = markerGraph.reverseComplementVertex[vertexId];
                if(vertexId > rcVertexId) {
                    continue; // not canonical
                }

                const auto vertexMarkerIds = mgVertices[vertexId];
                if(vertexMarkerIds.size() < minAnchorCoverage) {
                    continue;
                }

                vertexIntervals.clear();
                // Keep at most one marker per oriented read.
                OrientedReadId::Int lastOrientedReadValue = std::numeric_limits<OrientedReadId::Int>::max();
                for(const MarkerId markerId : vertexMarkerIds) {
                    OrientedReadId orientedReadId;
                    uint32_t ordinal0;
                    tie(orientedReadId, ordinal0) = dinara::findMarkerId(markerId, *markers);
                    const OrientedReadId::Int orientedReadValue = orientedReadId.getValue();
                    if(orientedReadValue == lastOrientedReadValue) {
                        continue;
                    }
                    lastOrientedReadValue = orientedReadValue;
                    vertexIntervals.emplace_back(orientedReadId, ordinal0);
                }

                if(vertexIntervals.size() < minAnchorCoverage) {
                    continue;
                }

                // Build an overlap-support graph among the oriented reads in this vertex.
                index.clear();
                for(uint32_t u=0; u<uint32_t(vertexIntervals.size()); ++u) {
                    index.emplace(vertexIntervals[u].orientedReadId.getValue(), u);
                }

                const uint32_t n = uint32_t(vertexIntervals.size());
                adjAll.assign(n, {});
                adjCis.assign(n, {});
                bool hasAnyCisEdge = false;
                for(uint32_t u=0; u<n; ++u) {
                    const OrientedReadId orientedReadId = vertexIntervals[u].orientedReadId;
                    for(const uint32_t edgeId : readGraph.connectivity[orientedReadId.getValue()]) {
                        const ReadGraphEdge& edge = readGraph.edges[edgeId];
                        if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                            continue;
                        }
                        const uint64_t alignmentId = edge.alignmentId;
                        const AlignmentData& ad = alignmentData[alignmentId];
                        if(!ad.keptByBothSides()) {
                            continue;
                        }
                        if(!ad.info.isInReadGraph) {
                            continue;
                        }
                        if(ad.cisTransStatus == CisTransStatus::Trans) {
                            continue;
                        }
                        const OrientedReadId other = edge.getOther(orientedReadId);
                        const auto it = index.find(other.getValue());
                        if(it == index.end()) {
                            continue;
                        }
                        const uint32_t vIdx = it->second;
                        if(vIdx == u) {
                            continue;
                        }
                        const uint32_t ordinalU = vertexIntervals[u].ordinal0;
                        const uint32_t ordinalV = vertexIntervals[vIdx].ordinal0;
                        const AlignmentInfo info = ad.orient(orientedReadId, other);
                        if(ordinalU < info.data[0].firstOrdinal || ordinalU > info.data[0].lastOrdinal) {
                            continue;
                        }
                        if(ordinalV < info.data[1].firstOrdinal || ordinalV > info.data[1].lastOrdinal) {
                            continue;
                        }
                        adjAll[u].push_back(vIdx);
                        if(ad.cisTransStatus == CisTransStatus::Cis) {
                            adjCis[u].push_back(vIdx);
                            hasAnyCisEdge = true;
                        }
                    }
                }
                for(uint32_t u=0; u<n; ++u) {
                    auto& nbr = adjAll[u];
                    std::sort(nbr.begin(), nbr.end());
                    nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
                }
                for(uint32_t u=0; u<n; ++u) {
                    auto& nbr = adjCis[u];
                    std::sort(nbr.begin(), nbr.end());
                    nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
                }

                const auto groups = splitVertexByOverlapSupportWithPhasing(
                    adjAll,
                    adjCis,
                    hasAnyCisEdge,
                    minAnchorCoverage,
                    maxAnchorCoverage);
                for(const auto& group : groups) {
                    vector<Interval> anchor;
                    anchor.reserve(group.size());
                    for(const uint32_t idxInVertex : group) {
                        anchor.push_back(vertexIntervals[idxInVertex]);
                    }
                    std::sort(anchor.begin(), anchor.end(), [](const Interval& a, const Interval& b) {
                        return a.orientedReadId < b.orientedReadId;
                    });
                    if(anchor.size() < minAnchorCoverage || anchor.size() > maxAnchorCoverage) {
                        continue;
                    }
                    outAnchors.push_back(anchor);
                    outAnchors.push_back(reverseComplementAnchor(anchor, *markers));
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    vector<vector<Interval>> anchorsExplicit;
    {
        size_t total = 0;
        for(const auto& v : threadAnchors) total += v.size();
        anchorsExplicit.reserve(total);
        for(auto& v : threadAnchors) {
            anchorsExplicit.insert(anchorsExplicit.end(), v.begin(), v.end());
        }
    }

    cout << timestamp << "Constructed " << anchorsExplicit.size()
         << " explicit anchors (including reverse complements) after decomposition." << endl;

    return make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        anchorsExplicit,
        /*ordinalOffset*/ 0,
        threadCount);
}



shared_ptr<mode3::Anchors> Assembler::createAnchorsFromOverlapsBestPerOverlapInterval(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount)
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    computeMarkerKmerIds(threadCount);
    checkReadGraphIsOpen();
    if(!alignmentData.isOpen) {
        throw runtime_error("Alignment data are not accessible.");
    }

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    // Only generate anchors from strand-0 oriented reads (read+).
    // We always emit the reverse-complement anchor explicitly, so processing strand-1 reads would
    // create duplicate anchors and slow down anchor generation.
    const uint64_t orientedReadCount = readCount;
    const uint64_t k = assemblerInfo->k;
    // Same defaults as Shasta2 (--max-anchor-repeat-length).
    const vector<uint64_t> maxAnchorRepeatLength = {8, 3, 3, 3, 3};

    // Experimental: If we have informative het coverage (set by AssemblerHifiasmEC),
    // only generate anchors from reads that cover at least one informative het site,
    // and only using their cis overlaps (trans overlaps are excluded by readGraph filtering,
    // but we also explicitly skip those if cisTransStatus is populated).
    bool restrictToInformativeHetReads = false;
    vector<uint8_t> readHasInformativeHet(readCount, 0);
    {
        // We treat the dataset as "informative-het mode" if ANY alignment is flagged as
        // covering an informative het/SV site. This is independent of read graph construction
        // and avoids missing directional phasing cases that were filtered out by conservative AND.
        for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); ++alignmentId) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if(!ad.coversHetSite) {
                continue;
            }
            restrictToInformativeHetReads = true;
            if(ad.readIds[0] < readCount) {
                readHasInformativeHet[ad.readIds[0]] = 1;
            }
            if(ad.readIds[1] < readCount) {
                readHasInformativeHet[ad.readIds[1]] = 1;
            }
        }
    }
    if(restrictToInformativeHetReads) {
        cout << timestamp << "Using informative-het-only mode for overlap-only anchors: "
             << "only reads covering an informative het site will contribute anchors, using cis overlaps only." << endl;
    } else {
        cout << timestamp << "No informative het coverage detected in the filtered read graph; "
             << "generating overlap-only anchors from all reads." << endl;
    }

    if((not markerKmerIds) or (not markerKmerIds->isOpen())) {
        throw runtime_error("MarkerKmerIds are required for FromOverlapsBestPerOverlapInterval (duplicate ReadId filter).");
    }

    // Precompute the set of canonical marker k-mers that have duplicate ReadIds in MarkerKmerIds
    // (meaning the same canonical k-mer appears more than once in at least one read).
    // These are problematic as anchor seeds because the seed k-mer is not unique within some read.
    vector<vector<KmerId>> duplicateCanonicalKmerIdsByThread(threadCount);
    {
        vector<thread> dupThreads;
        dupThreads.reserve(threadCount);

        uint64_t chunk2 = readCount / threadCount;
        if(chunk2 == 0) chunk2 = 1;

        for(uint64_t t=0; t<threadCount; t++) {
            dupThreads.emplace_back([&, t]() {
                const uint64_t begin = t * chunk2;
                const uint64_t end = (t == threadCount - 1) ? readCount : min(readCount, (t+1) * chunk2);

                auto& duplicates = duplicateCanonicalKmerIdsByThread[t];
                duplicates.clear();
                duplicates.reserve(1024);

                unordered_set<KmerId, KmerIdHasher> seen;
                seen.reserve(4096);

                for(uint64_t r=begin; r<end; ++r) {
                    const ReadId readId = ReadId(r);
                    const OrientedReadId or0(readId, 0);
                    const OrientedReadId or1(readId, 1);
                    const auto kmerIds0 = (*markerKmerIds)[or0.getValue()];
                    const auto kmerIds1 = (*markerKmerIds)[or1.getValue()];
                    const uint32_t markerCount = uint32_t(kmerIds0.size());
                    if(markerCount == 0) {
                        continue;
                    }
                    DINARA_ASSERT(kmerIds1.size() == markerCount);

                    seen.clear();
                    if(seen.bucket_count() < size_t(markerCount) * 2ULL) {
                        seen.reserve(size_t(markerCount) * 2ULL);
                    }

                    for(uint32_t ordinal0=0; ordinal0<markerCount; ++ordinal0) {
                        const uint32_t ordinal1 = markerCount - 1U - ordinal0;
                        const KmerId id0 = kmerIds0[ordinal0];
                        const KmerId id1 = kmerIds1[ordinal1];
                        const KmerId canonical = (id0 <= id1) ? id0 : id1;
                        const auto inserted = seen.insert(canonical).second;
                        if(not inserted) {
                            duplicates.push_back(canonical);
                        }
                    }
                }
            });
        }
        for(auto& th : dupThreads) {
            th.join();
        }
    }

    vector<KmerId> duplicateCanonicalKmerIds;
    {
        size_t total = 0;
        for(const auto& v : duplicateCanonicalKmerIdsByThread) total += v.size();
        duplicateCanonicalKmerIds.reserve(total);
        for(auto& v : duplicateCanonicalKmerIdsByThread) {
            duplicateCanonicalKmerIds.insert(duplicateCanonicalKmerIds.end(), v.begin(), v.end());
        }
        sort(duplicateCanonicalKmerIds.begin(), duplicateCanonicalKmerIds.end());
        duplicateCanonicalKmerIds.erase(
            unique(duplicateCanonicalKmerIds.begin(), duplicateCanonicalKmerIds.end()),
            duplicateCanonicalKmerIds.end());
    }

    unordered_set<KmerId, KmerIdHasher> canonicalKmerIdsWithDuplicateReadIds;
    canonicalKmerIdsWithDuplicateReadIds.reserve(duplicateCanonicalKmerIds.size() * 2ULL + 1ULL);
    for(const KmerId& id : duplicateCanonicalKmerIds) {
        canonicalKmerIdsWithDuplicateReadIds.insert(id);
    }
    cout << timestamp << "Identified " << canonicalKmerIdsWithDuplicateReadIds.size()
         << " canonical marker k-mers with duplicate ReadIds (will not be used as anchor seeds)." << endl;

    struct CandidateAnchor {
        ReadId seedReadId;
        uint32_t overlapIntervalIndex = 0;
        uint32_t intervalStart = 0;
        uint32_t intervalEnd = 0; // one past last
        uint32_t seedOrdinal = 0;
        uint32_t support = 0;
        vector<Interval> anchor;
    };

    // For each read (strand 0), store the overlap-event intervals (where the active overlap set is constant and non-empty).
    vector<vector<OverlapInterval>> overlapIntervalsPerRead(readCount);

    vector<vector<CandidateAnchor>> threadCandidates(threadCount);

    uint64_t chunk = orientedReadCount / threadCount;
    if(chunk == 0) chunk = 1;

    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = (t == threadCount - 1) ? orientedReadCount : min(orientedReadCount, (t+1) * chunk);
            auto& outCandidates = threadCandidates[t];
            outCandidates.reserve(end - begin);

            vector<OverlapEvent> events;
            events.reserve(512);
            vector<uint32_t> activeEdgeIds;
            activeEdgeIds.reserve(256);

            vector<uint32_t> segmentEdgeIds;
            segmentEdgeIds.reserve(256);

            for(uint64_t v=begin; v<end; ++v) {
                const OrientedReadId orientedReadId0(ReadId(v), 0);
                const uint32_t markerCount0 = uint32_t(markers->size(orientedReadId0.getValue()));
                if(markerCount0 == 0) {
                    continue;
                }
                const OrientedReadId orientedReadId0rc(ReadId(v), 1);
                const auto orientedReadKmerIds0 = (*markerKmerIds)[orientedReadId0.getValue()];
                const auto orientedReadKmerIds1 = (*markerKmerIds)[orientedReadId0rc.getValue()];
                DINARA_ASSERT(orientedReadKmerIds0.size() == markerCount0);
                DINARA_ASSERT(orientedReadKmerIds1.size() == markerCount0);

                auto& outIntervals = overlapIntervalsPerRead[v];
                outIntervals.clear();

                events.clear();
                if(restrictToInformativeHetReads) {
                    if(v >= readHasInformativeHet.size() || !readHasInformativeHet[v]) {
                        continue;
                    }
                }

                // Collect start/end events from filtered readGraph overlaps.
                for(const uint32_t edgeId : readGraph.connectivity[orientedReadId0.getValue()]) {
                    const ReadGraphEdge& edge = readGraph.edges[edgeId];
                    if(edge.crossesStrands) {
                        continue;
                    }
                    if(edge.hasInconsistentAlignment) {
                        continue;
                    }
                    const uint64_t alignmentId = edge.alignmentId;
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if(!ad.keptByBothSides()) {
                        continue;
                    }
                    if(!ad.info.isInReadGraph) {
                        continue;
                    }
                    if(restrictToInformativeHetReads) {
                        // Use directional phasing decisions:
                        // if this read marked the overlap as trans, it will have DeleteReasonPhase set
                        // from this read's perspective (even if the other read kept it).
                        const ReadId seedReadId = orientedReadId0.getReadId();
                        const AlignmentData::DeleteReasonMask reasons =
                            (ad.readIds[0] == seedReadId) ? ad.deleteReasons0 : ad.deleteReasons1;
                        if(reasons & AlignmentData::DeleteReasonPhase) {
                            continue;
                        }
                    }

                    const OrientedReadId orientedReadId1 = edge.getOther(orientedReadId0);
                    const AlignmentInfo info = ad.orient(orientedReadId0, orientedReadId1);
                    const uint32_t first = info.data[0].firstOrdinal;
                    const uint32_t last = info.data[0].lastOrdinal;
                    if(first < markerCount0) {
                        events.push_back({first, +1, edgeId});
                    }
                    const uint32_t afterLast = last + 1;
                    if(afterLast < markerCount0) {
                        events.push_back({afterLast, -1, edgeId});
                    }
                }

                if(events.empty()) {
                    continue;
                }

                std::sort(events.begin(), events.end(), [](const OverlapEvent& a, const OverlapEvent& b) {
                    if(a.ordinal != b.ordinal) {
                        return a.ordinal < b.ordinal;
                    }
                    // Deterministic: process end events before start events at the same ordinal.
                    return a.delta < b.delta;
                });

                auto deactivate = [&](uint32_t edgeIdToRemove) {
                    const auto it = std::find(activeEdgeIds.begin(), activeEdgeIds.end(), edgeIdToRemove);
                    if(it != activeEdgeIds.end()) {
                        *it = activeEdgeIds.back();
                        activeEdgeIds.pop_back();
                    }
                };

                auto shouldSkipKmerDueToRepeats = [&](const Kmer& kmer0) -> bool {
                    for(uint64_t i=0; i<maxAnchorRepeatLength.size(); i++) {
                        const uint64_t period = i + 1;
                        const uint64_t maxAllowedCopyNumber = maxAnchorRepeatLength[i];
                        uint64_t copies = 0;
                        switch(period) {
                        case 1: copies = kmer0.countExactRepeatCopies<1>(k); break;
                        case 2: copies = kmer0.countExactRepeatCopies<2>(k); break;
                        case 3: copies = kmer0.countExactRepeatCopies<3>(k); break;
                        case 4: copies = kmer0.countExactRepeatCopies<4>(k); break;
                        case 5: copies = kmer0.countExactRepeatCopies<5>(k); break;
                        case 6: copies = kmer0.countExactRepeatCopies<6>(k); break;
                        default:
                            // Shasta2 only supports up to period 6.
                            copies = 0;
                            break;
                        }
                        if(copies > maxAllowedCopyNumber) {
                            return true;
                        }
                    }
                    return false;
                };

                auto buildAnchorAtSeed = [&](uint32_t seedOrdinal, const Kmer& seedKmer, const vector<uint32_t>& edgeIds, vector<Interval>& anchorOut) -> void {
                    anchorOut.clear();
                    anchorOut.reserve(64);
                    anchorOut.emplace_back(orientedReadId0, seedOrdinal);

                    std::unordered_set<ReadId> usedReadIds;
                    usedReadIds.reserve(128);
                    usedReadIds.insert(orientedReadId0.getReadId());

                    for(const uint32_t activeEdgeId : edgeIds) {
                        const ReadGraphEdge& edge = readGraph.edges[activeEdgeId];
                        if(edge.crossesStrands) {
                            continue;
                        }
                        if(edge.hasInconsistentAlignment) {
                            continue;
                        }
                        const uint64_t alignmentId = edge.alignmentId;
                        const AlignmentData& ad = alignmentData[alignmentId];
                        if(!ad.keptByBothSides()) {
                            continue;
                        }
                        if(!ad.info.isInReadGraph) {
                            continue;
                        }

                        const OrientedReadId orientedReadId1 = edge.getOther(orientedReadId0);
                        const ReadId readId1 = orientedReadId1.getReadId();
                        if(usedReadIds.contains(readId1)) {
                            continue;
                        }

                        const AlignmentInfo info = ad.orient(orientedReadId0, orientedReadId1);
                        uint32_t ordinal1 = 0;
                        if(!mapMarkerOrdinalByOffsetAndKmer(
                            *reads,
                            *markers,
                            k,
                            info,
                            orientedReadId0,
                            seedOrdinal,
                            orientedReadId1,
                            ordinal1,
                            seedKmer,
                            /*maxSearchRadius*/ 8,
                            /*maxOffsetRange*/ 32)) {
                            continue;
                        }

                        anchorOut.emplace_back(orientedReadId1, ordinal1);
                        usedReadIds.insert(readId1);
                        if(anchorOut.size() >= maxAnchorCoverage) {
                            break;
                        }
                    }
                };

                vector<Interval> tmpAnchor;
                vector<Interval> bestAnchor;

                activeEdgeIds.clear();
                size_t ei = 0;
                while(ei < events.size()) {
                    const uint32_t ordinal = events[ei].ordinal;

                    // Apply all events at this ordinal.
                    while(ei < events.size() && events[ei].ordinal == ordinal) {
                        if(events[ei].delta > 0) {
                            activeEdgeIds.push_back(events[ei].edgeId);
                        } else {
                            deactivate(events[ei].edgeId);
                        }
                        ++ei;
                    }

                    const uint32_t nextOrdinal = (ei < events.size()) ? events[ei].ordinal : markerCount0;
                    if(activeEdgeIds.empty()) {
                        continue;
                    }
                    if(nextOrdinal <= ordinal) {
                        continue;
                    }

                    const uint32_t segmentStart = ordinal;
                    const uint32_t segmentEnd = nextOrdinal;

                    // This overlap-event segment has a constant active overlap set.
                    // We further split long segments into shorter anchor intervals to avoid large gaps.
                    constexpr uint32_t maxIntervalMarkers = 200;

                    // Early exit: even with perfect k-mer mapping, we cannot exceed 1 + (# unique reads in active edges).
                    segmentEdgeIds = activeEdgeIds;
                    std::sort(segmentEdgeIds.begin(), segmentEdgeIds.end(), [&](uint32_t a, uint32_t b) {
                        const ReadId ra = readGraph.edges[a].getOther(orientedReadId0).getReadId();
                        const ReadId rb = readGraph.edges[b].getOther(orientedReadId0).getReadId();
                        if(ra != rb) {
                            return ra < rb;
                        }
                        return a < b;
                    });
                    {
                        ReadId prev = invalid<ReadId>;
                        uint32_t uniqueOtherReads = 0;
                        for(const uint32_t edgeId : segmentEdgeIds) {
                            const ReadId r = readGraph.edges[edgeId].getOther(orientedReadId0).getReadId();
                            if(r != prev) {
                                ++uniqueOtherReads;
                                prev = r;
                            }
                        }
                        const uint32_t maxPossible = std::min<uint32_t>(uint32_t(maxAnchorCoverage), 1 + uniqueOtherReads);
                        if(maxPossible < minAnchorCoverage) {
                            continue;
                        }
                    }

                    // Select one anchor per (possibly split) anchor interval.
                    for(uint32_t intervalStart = segmentStart; intervalStart < segmentEnd; intervalStart += maxIntervalMarkers) {
                        const uint32_t intervalEnd = std::min(segmentEnd, intervalStart + maxIntervalMarkers);
                        const uint32_t intervalLen = intervalEnd - intervalStart;
                        if(intervalLen == 0) {
                            continue;
                        }

                        outIntervals.push_back({intervalStart, intervalEnd});
                        const uint32_t overlapIntervalIndex = uint32_t(outIntervals.size() - 1);

                        // Pick the seed marker with maximum support by scanning all marker ordinals in the interval.
                        vector<uint32_t> seedCandidates;
                        seedCandidates.reserve(intervalLen);
                        for(uint32_t o = intervalStart; o < intervalEnd; ++o) {
                            seedCandidates.push_back(o);
                        }

                        uint32_t bestSeed = invalid<uint32_t>;
                        uint32_t bestSupport = 0;
                        bestAnchor.clear();

                        auto seedBetter = [&](uint32_t a, uint32_t b) -> bool {
                            // Prefer the seed closer to the interval center (more robust), then lower ordinal for determinism.
                            const int64_t center = int64_t(intervalStart) + int64_t(intervalLen) / 2;
                            const int64_t da = (int64_t(a) >= center) ? (int64_t(a) - center) : (center - int64_t(a));
                            const int64_t db = (int64_t(b) >= center) ? (int64_t(b) - center) : (center - int64_t(b));
                            if(da != db) {
                                return da < db;
                            }
                            return a < b;
                        };

                        for(const uint32_t seedOrdinal : seedCandidates) {
                            if(seedOrdinal < intervalStart || seedOrdinal >= intervalEnd) {
                                continue;
                            }

                            // Skip this seed if its canonical marker k-mer has duplicate ReadIds globally
                            // (it appears more than once in at least one read).
                            const uint32_t seedOrdinalRc = markerCount0 - 1U - seedOrdinal;
                            const KmerId id0 = orientedReadKmerIds0[seedOrdinal];
                            const KmerId id1 = orientedReadKmerIds1[seedOrdinalRc];
                            const KmerId canonicalId = (id0 <= id1) ? id0 : id1;
                            if(canonicalKmerIdsWithDuplicateReadIds.contains(canonicalId)) {
                                continue;
                            }

                            const Kmer seedKmer = getMarkerKmer(*reads, *markers, k, orientedReadId0, seedOrdinal);
                            if(shouldSkipKmerDueToRepeats(seedKmer)) {
                                continue;
                            }
                            buildAnchorAtSeed(seedOrdinal, seedKmer, segmentEdgeIds, tmpAnchor);
                            const uint32_t support = uint32_t(tmpAnchor.size());
                            if(support > bestSupport || (support == bestSupport && bestSeed != invalid<uint32_t> && seedBetter(seedOrdinal, bestSeed))) {
                                bestSupport = support;
                                bestSeed = seedOrdinal;
                                bestAnchor = tmpAnchor;
                                if(bestSupport >= maxAnchorCoverage) {
                                    break;
                                }
                            }
                        }

                        if(bestSupport < minAnchorCoverage || bestSupport > maxAnchorCoverage) {
                            continue;
                        }

                        std::sort(bestAnchor.begin(), bestAnchor.end(), [](const Interval& a, const Interval& b) {
                            return a.orientedReadId < b.orientedReadId;
                        });

                        CandidateAnchor candidate;
                        candidate.seedReadId = ReadId(v);
                        candidate.overlapIntervalIndex = overlapIntervalIndex;
                        candidate.intervalStart = intervalStart;
                        candidate.intervalEnd = intervalEnd;
                        candidate.seedOrdinal = bestSeed;
                        candidate.support = bestSupport;
                        candidate.anchor = std::move(bestAnchor);
                        outCandidates.push_back(std::move(candidate));
                    }
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    vector<CandidateAnchor> candidates;
    {
        size_t total = 0;
        for(const auto& v : threadCandidates) total += v.size();
        candidates.reserve(total);
        for(auto& v : threadCandidates) {
            candidates.insert(candidates.end(), v.begin(), v.end());
        }
    }

    cout << timestamp << "Constructed " << candidates.size()
         << " candidate anchors from overlaps (one per overlap-event interval, with long-interval splitting)." << endl;

    // Select anchors with cross-read overlap-interval claiming:
    // if a read already has an anchor in a given overlap-event interval, skip selecting another anchor for that interval.
    // This reduces duplicates and spreads anchors across overlap-covered regions.
    vector<vector<uint8_t>> intervalClaimed(readCount);
    for(uint64_t v=0; v<readCount; ++v) {
        intervalClaimed[v].assign(overlapIntervalsPerRead[v].size(), 0);
    }

    std::sort(candidates.begin(), candidates.end(), [](const CandidateAnchor& a, const CandidateAnchor& b) {
        if(a.support != b.support) {
            return a.support > b.support; // stronger anchors first
        }
        if(a.seedReadId != b.seedReadId) {
            return a.seedReadId < b.seedReadId;
        }
        if(a.intervalStart != b.intervalStart) {
            return a.intervalStart < b.intervalStart;
        }
        if(a.intervalEnd != b.intervalEnd) {
            return a.intervalEnd < b.intervalEnd;
        }
        return a.seedOrdinal < b.seedOrdinal;
    });

    vector<vector<Interval>> selected;
    selected.reserve(candidates.size());

    auto claimIntervalForMarker = [&](const Interval& interval) -> std::optional<pair<ReadId, uint32_t>> {
        const ReadId readId = interval.orientedReadId.getReadId();
        if(uint64_t(readId) >= readCount) {
            return std::nullopt;
        }
        uint32_t idx = 0;
        if(!findOverlapIntervalIndex(overlapIntervalsPerRead[readId], interval.ordinal0, idx)) {
            return std::nullopt;
        }
        return pair<ReadId, uint32_t>(readId, idx);
    };

    auto selectCandidate = [&](const CandidateAnchor& candidate, bool allowClaimedInRepair) -> bool {
        const uint64_t v = uint64_t(candidate.seedReadId);
        if(v >= readCount) {
            return false;
        }
        if(candidate.overlapIntervalIndex >= overlapIntervalsPerRead[v].size()) {
            return false;
        }
        if(intervalClaimed[v][candidate.overlapIntervalIndex]) {
            return false;
        }

        vector<Interval> anchor;
        anchor.reserve(candidate.anchor.size());

        vector<pair<ReadId, uint32_t>> toClaim;
        toClaim.reserve(candidate.anchor.size());

        // Always include the seed marker, so we anchor this interval on the seed read.
        anchor.emplace_back(OrientedReadId(candidate.seedReadId, 0), candidate.seedOrdinal);
        toClaim.emplace_back(candidate.seedReadId, candidate.overlapIntervalIndex);

        // Add additional read intervals from this candidate, preferring unclaimed overlap-intervals.
        // If allowClaimedInRepair is true, we can fall back to using already-claimed intervals to reach min coverage.
        vector<pair<Interval, pair<ReadId, uint32_t>>> claimedDeferred;
        claimedDeferred.reserve(candidate.anchor.size());

        for(const Interval& interval : candidate.anchor) {
            const ReadId readId = interval.orientedReadId.getReadId();
            if(readId == candidate.seedReadId) {
                continue;
            }

            const auto keyOpt = claimIntervalForMarker(interval);
            if(!keyOpt) {
                continue;
            }
            const auto [rid, idx] = *keyOpt;
            if(intervalClaimed[rid][idx]) {
                if(allowClaimedInRepair) {
                    claimedDeferred.push_back({interval, {rid, idx}});
                }
                continue;
            }
            anchor.push_back(interval);
            toClaim.push_back({rid, idx});
            if(anchor.size() >= maxAnchorCoverage) {
                break;
            }
        }

        if(allowClaimedInRepair && anchor.size() < minAnchorCoverage) {
            for(const auto& x : claimedDeferred) {
                anchor.push_back(x.first);
                if(anchor.size() >= minAnchorCoverage || anchor.size() >= maxAnchorCoverage) {
                    break;
                }
            }
        }

        if(anchor.size() < minAnchorCoverage || anchor.size() > maxAnchorCoverage) {
            return false;
        }

        std::sort(anchor.begin(), anchor.end(), [](const Interval& a, const Interval& b) {
            return a.orientedReadId < b.orientedReadId;
        });

        selected.push_back(anchor);

        // Claim all unclaimed overlap-event intervals used by this anchor.
        for(const auto& [rid, idx] : toClaim) {
            intervalClaimed[rid][idx] = 1;
        }
        return true;
    };

    // Pass 1: strict (only uses unclaimed intervals on all reads).
    uint64_t selectedStrict = 0;
    for(const CandidateAnchor& candidate : candidates) {
        if(selectCandidate(candidate, /*allowClaimedInRepair*/ false)) {
            ++selectedStrict;
        }
    }

    // Pass 2 (repair) disabled for now.
    // It can increase coverage by allowing already-claimed intervals on other reads,
    // but it changes the strict "one anchor claims one interval" behavior.
    uint64_t selectedRepair = 0;

    cout << timestamp << "Selected " << selected.size()
         << " anchors (strict=" << selectedStrict << ", repair=" << selectedRepair << ")." << endl;

    uint64_t totalIntervals = 0;
    uint64_t claimedIntervals = 0;
    for(uint64_t v=0; v<readCount; ++v) {
        totalIntervals += intervalClaimed[v].size();
        for(const uint8_t x : intervalClaimed[v]) {
            claimedIntervals += (x != 0);
        }
    }
    if(totalIntervals) {
        cout << timestamp << "Claimed " << claimedIntervals << " / " << totalIntervals
             << " overlap intervals (" << double(claimedIntervals) / double(totalIntervals) << ")." << endl;
    }

    // Emit explicit reverse complements.
    vector<vector<Interval>> anchorsExplicit;
    anchorsExplicit.reserve(selected.size() * 2);
    for(const auto& anchor : selected) {
        anchorsExplicit.push_back(anchor);
        anchorsExplicit.push_back(reverseComplementAnchor(anchor, *markers));
    }

    // Make output deterministic.
    auto anchorLessLex = [](const vector<Interval>& a, const vector<Interval>& b) -> bool {
        const size_t n = std::min(a.size(), b.size());
        for(size_t i=0; i<n; ++i) {
            if(a[i].orientedReadId != b[i].orientedReadId) {
                return a[i].orientedReadId < b[i].orientedReadId;
            }
            if(a[i].ordinal0 != b[i].ordinal0) {
                return a[i].ordinal0 < b[i].ordinal0;
            }
        }
        return a.size() < b.size();
    };
    std::sort(anchorsExplicit.begin(), anchorsExplicit.end(), [&](const vector<Interval>& a, const vector<Interval>& b) {
        if(a.size() != b.size()) {
            return a.size() > b.size();
        }
        return anchorLessLex(a, b);
    });
    anchorsExplicit.erase(std::unique(anchorsExplicit.begin(), anchorsExplicit.end(),
        [&](const vector<Interval>& a, const vector<Interval>& b) {
            return a.size() == b.size() && !anchorLessLex(a, b) && !anchorLessLex(b, a);
        }), anchorsExplicit.end());

    cout << timestamp << "Selected " << anchorsExplicit.size()
         << " anchors (including reverse complements) after deduplication." << endl;

    auto anchors = make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        anchorsExplicit,
        /*ordinalOffset*/ 0,
        threadCount);

    // These anchors are used for mode3-style exploration (and sometimes for downstream processing).
    // The local anchor graph and other navigation features require journeys.
    anchors->computeJourneys(threadCount);
    return anchors;
}
