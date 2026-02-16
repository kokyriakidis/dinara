// ============================================================================
// BRG-native anchor creation, journey computation, and edge table.
//
// Creates one BrgAnchor per marker graph vertex (no RC doubling).
// Self-RC vertices naturally contain both-strand markers.
// Non-self-RC vertices are handled but expected to be rare after merging.
//
// Journeys are per-ReadId: the ordered sequence of (BidirectedAnchorId, strand)
// entries along each physical read.
//
// The edge table precomputes directed edges from journey adjacencies,
// with the bidirectional symmetry invariant enforced:
//   if edge (A,sA)→(B,sB) exists, then (reverse(B),reverse(A)) also exists.
// ============================================================================

#include "mode3-BidirectedAnchor.hpp"
#include "DINARA_ASSERT.hpp"
#include "findMarkerId.hpp"
#include "invalid.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"

// For MultithreadedObject explicit instantiation.
#include "MultithreadedObject.tpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>

using namespace dinara;
using namespace mode3;
using namespace std;

namespace dinara {
    template class MultithreadedObject<BidirectedAnchors>;
}


// ============================================================================
// Construction from marker graph vertices.
// ============================================================================

BidirectedAnchors::BidirectedAnchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    uint64_t minCoverage,
    uint64_t maxCoverage,
    uint64_t threadCount) :
    MultithreadedObject<BidirectedAnchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    performanceLog << timestamp << "BidirectedAnchors creation begins." << endl;

    DINARA_ASSERT(markers.isOpen());
    DINARA_ASSERT(markerGraph.vertices().isOpen());
    DINARA_ASSERT(markerGraph.vertexTable.isOpen);
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);

    using VertexId = MarkerGraph::VertexId;
    const auto& mgVertices = markerGraph.vertices();
    const VertexId vertexCount = markerGraph.vertexCount();

    // Collect BidirectedAnchors into a temporary vector, then populate the
    // memory-mapped VectorOfVectors.
    vector<vector<BidirectedAnchorMarkerInterval>> anchorsTemp;
    anchorsTemp.reserve(vertexCount / 2 + 1);  // rough upper bound

    uint64_t selfRcVertices = 0;
    uint64_t nonSelfRcVertices = 0;
    uint64_t skippedCoverage = 0;

    for(VertexId vertexId = 0; vertexId < vertexCount; ++vertexId) {

        const VertexId rcVertexId = markerGraph.reverseComplementVertex[vertexId];

        // Only process canonical vertices.
        // For self-RC: vertexId == rcVertexId, always processed.
        // For non-self-RC pairs: process only when vertexId < rcVertexId.
        if(vertexId > rcVertexId) {
            continue;
        }

        const auto vertexMarkerIds = mgVertices[vertexId];
        const bool isSelfRc = (vertexId == rcVertexId);

        if(isSelfRc) {
            // Self-RC vertex: contains markers from both strands.
            // Extract one interval per physical read.
            // For each ReadId, take the strand-0 marker (the strand indicates
            // direction of traversal).

            vector<BidirectedAnchorMarkerInterval> anchor;
            anchor.reserve(vertexMarkerIds.size() / 2 + 1);

            for(const MarkerId markerId : vertexMarkerIds) {
                OrientedReadId orientedReadId;
                uint32_t ordinal;
                tie(orientedReadId, ordinal) = findMarkerId(markerId, markers);

                const ReadId readId = orientedReadId.getReadId();
                const Strand strand = orientedReadId.getStrand();

                // Store ordinal on strand 0 for consistent coordinates.
                uint32_t ordinalOnStrand0;
                if(strand == 0) {
                    ordinalOnStrand0 = ordinal;
                } else {
                    // Convert strand-1 ordinal to strand-0 ordinal.
                    const uint64_t markerCount = markers.size(
                        OrientedReadId(readId, 0).getValue());
                    ordinalOnStrand0 = uint32_t(markerCount) - 1 - ordinal;
                }

                anchor.push_back(BidirectedAnchorMarkerInterval(readId, ordinalOnStrand0, strand));
            }

            // Sort by ReadId for deduplication.
            sort(anchor.begin(), anchor.end(),
                [](const BidirectedAnchorMarkerInterval& a, const BidirectedAnchorMarkerInterval& b) {
                    return a.readId < b.readId;
                });

            // Deduplicate by ReadId: keep the strand-0 entry if both strands present.
            // (Both strands of the same read will be adjacent after sorting.)
            {
                vector<BidirectedAnchorMarkerInterval> deduped;
                deduped.reserve(anchor.size());
                for(size_t i = 0; i < anchor.size(); ) {
                    // Find run of same ReadId.
                    size_t j = i + 1;
                    while(j < anchor.size() && anchor[j].readId == anchor[i].readId) {
                        ++j;
                    }
                    // Pick strand-0 if available, else strand-1.
                    const BidirectedAnchorMarkerInterval* best = &anchor[i];
                    for(size_t ii = i; ii < j; ++ii) {
                        if(anchor[ii].strand == 0) {
                            best = &anchor[ii];
                            break;
                        }
                    }
                    deduped.push_back(*best);
                    i = j;
                }
                anchor = std::move(deduped);
            }

            // Coverage check.
            if(anchor.size() < minCoverage || anchor.size() > maxCoverage) {
                ++skippedCoverage;
                continue;
            }

            anchorsTemp.push_back(std::move(anchor));
            ++selfRcVertices;

        } else {
            // Non-self-RC vertex pair: use the canonical vertex (vertexId),
            // which contains markers from one strand only.
            // We still create a single BrgAnchor from it.

            vector<BidirectedAnchorMarkerInterval> anchor;
            anchor.reserve(vertexMarkerIds.size());

            for(const MarkerId markerId : vertexMarkerIds) {
                OrientedReadId orientedReadId;
                uint32_t ordinal;
                tie(orientedReadId, ordinal) = findMarkerId(markerId, markers);

                const ReadId readId = orientedReadId.getReadId();
                const Strand strand = orientedReadId.getStrand();

                uint32_t ordinalOnStrand0;
                if(strand == 0) {
                    ordinalOnStrand0 = ordinal;
                } else {
                    const uint64_t markerCount = markers.size(
                        OrientedReadId(readId, 0).getValue());
                    ordinalOnStrand0 = uint32_t(markerCount) - 1 - ordinal;
                }

                anchor.push_back(BidirectedAnchorMarkerInterval(readId, ordinalOnStrand0, strand));
            }

            // Sort and deduplicate by ReadId.
            sort(anchor.begin(), anchor.end(),
                [](const BidirectedAnchorMarkerInterval& a, const BidirectedAnchorMarkerInterval& b) {
                    return a.readId < b.readId;
                });
            anchor.erase(
                unique(anchor.begin(), anchor.end(),
                    [](const BidirectedAnchorMarkerInterval& a, const BidirectedAnchorMarkerInterval& b) {
                        return a.readId == b.readId;
                    }),
                anchor.end());

            // Check for true ReadId duplicates in the original vertex.
            if(markerGraph.vertexHasDuplicateReadIds(vertexId, markers)) {
                continue;
            }

            if(anchor.size() < minCoverage || anchor.size() > maxCoverage) {
                ++skippedCoverage;
                continue;
            }

            anchorsTemp.push_back(std::move(anchor));
            ++nonSelfRcVertices;
        }
    }

    anchorCount = anchorsTemp.size();

    cout << timestamp << "BidirectedAnchors: " << anchorCount << " anchors created ("
         << selfRcVertices << " from self-RC vertices, "
         << nonSelfRcVertices << " from non-self-RC vertex pairs, "
         << skippedCoverage << " skipped for coverage)." << endl;

    // Populate memory-mapped storage.
    anchorMarkerIntervals.createNew(
        largeDataName("BidirectedAnchorMarkerIntervals"), largeDataPageSize);
    for(const auto& anchor : anchorsTemp) {
        anchorMarkerIntervals.appendVector();
        for(const auto& interval : anchor) {
            anchorMarkerIntervals.append(interval);
        }
    }

    anchorInfos.createNew(
        largeDataName("BidirectedAnchorInfos"), largeDataPageSize);
    anchorInfos.resize(anchorCount);
    for(uint64_t i = 0; i < anchorCount; ++i) {
        anchorInfos[i].componentId = invalid<uint32_t>;
        anchorInfos[i].localAnchorIdInComponent = invalid<uint64_t>;
    }

    performanceLog << timestamp << "BidirectedAnchors creation ends." << endl;
}



// ============================================================================
// Access existing BidirectedAnchors from binary data (for HTTP server).
// ============================================================================

BidirectedAnchors::BidirectedAnchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) :
    MultithreadedObject<BidirectedAnchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    anchorMarkerIntervals.accessExistingReadOnly(
        largeDataName("BidirectedAnchorMarkerIntervals"));
    anchorInfos.accessExistingReadOnly(
        largeDataName("BidirectedAnchorInfos"));
    journeys.accessExistingReadOnly(
        largeDataName("BidirectedJourneys"));

    anchorCount = anchorMarkerIntervals.size();

    cout << timestamp << "Accessed existing BidirectedAnchors: "
         << anchorCount << " anchors, "
         << journeys.size() << " journey entries." << endl;

    // Recompute in-memory edge table from journeys.
    computeEdges(1);
}



// ============================================================================
// Accessors.
// ============================================================================

span<const BidirectedAnchorMarkerInterval> BidirectedAnchors::operator[](BidirectedAnchorId anchorId) const
{
    return anchorMarkerIntervals[anchorId];
}


uint64_t BidirectedAnchors::coverage(BidirectedAnchorId anchorId) const
{
    return anchorMarkerIntervals.size(anchorId);
}


// ============================================================================
// Journey computation.
//
// For each physical read, compute the ordered sequence of BidirectedAnchors it
// visits, along with the strand (direction) at each anchor.
//
// Algorithm:
// Pass 1-2: For each anchor, for each marker interval, record
//           (readId → {anchorId, strand, ordinal}).
// Pass 3:   Sort each read's entries by ordinal.
// Pass 4:   Copy sorted entries to journeys.
// ============================================================================

void BidirectedAnchors::computeJourneys(uint64_t threadCount)
{
    performanceLog << timestamp << "BidirectedAnchors::computeJourneys begins." << endl;

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads.readCount();

    // Pass 1: count entries per ReadId.
    journeysWithOrdinals.createNew(
        largeDataName("tmp-BidirectedJourneysWithOrdinals"), largeDataPageSize);
    journeysWithOrdinals.beginPass1(readCount);
    const uint64_t batchSize = 1000;
    setupLoadBalancing(anchorCount, batchSize);
    runThreads(&BidirectedAnchors::computeJourneysThreadFunction1, threadCount);

    // Pass 2: store entries.
    journeysWithOrdinals.beginPass2();
    setupLoadBalancing(anchorCount, batchSize);
    runThreads(&BidirectedAnchors::computeJourneysThreadFunction2, threadCount);
    journeysWithOrdinals.endPass2();

    // Pass 3: sort each read's entries by ordinal and allocate journeys.
    journeys.createNew(largeDataName("BidirectedJourneys"), largeDataPageSize);
    journeys.beginPass1(readCount);
    const uint64_t readBatchSize = 1000;
    setupLoadBalancing(readCount, readBatchSize);
    runThreads(&BidirectedAnchors::computeJourneysThreadFunction3, threadCount);

    // Pass 4: copy sorted entries into journeys.
    journeys.beginPass2();
    setupLoadBalancing(readCount, readBatchSize);
    runThreads(&BidirectedAnchors::computeJourneysThreadFunction4, threadCount);
    journeys.endPass2(false, true);

    journeysWithOrdinals.remove();

    performanceLog << timestamp << "BidirectedAnchors::computeJourneys ends." << endl;

    writeJourneys();
}


void BidirectedAnchors::computeJourneysThreadFunction1(uint64_t /* threadId */)
{
    computeJourneysThreadFunction12(1);
}


void BidirectedAnchors::computeJourneysThreadFunction2(uint64_t /* threadId */)
{
    computeJourneysThreadFunction12(2);
}


void BidirectedAnchors::computeJourneysThreadFunction12(uint64_t pass)
{
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(BidirectedAnchorId anchorId = begin; anchorId < end; ++anchorId) {
            const auto intervals = anchorMarkerIntervals[anchorId];
            for(const auto& interval : intervals) {
                const uint64_t readIdValue = uint64_t(interval.readId);
                if(pass == 1) {
                    journeysWithOrdinals.incrementCountMultithreaded(readIdValue);
                } else {
                    JourneyWithOrdinal entry;
                    entry.anchorId = anchorId;
                    entry.strand = interval.strand;
                    entry.ordinal = interval.ordinal;
                    journeysWithOrdinals.storeMultithreaded(readIdValue, entry);
                }
            }
        }
    }
}


void BidirectedAnchors::computeJourneysThreadFunction3(uint64_t /* threadId */)
{
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t readIdValue = begin; readIdValue < end; ++readIdValue) {
            auto v = journeysWithOrdinals[readIdValue];
            // Sort by ordinal so the journey follows the read's marker order.
            sort(v.begin(), v.end(),
                [](const JourneyWithOrdinal& a, const JourneyWithOrdinal& b) {
                    return a.ordinal < b.ordinal;
                });
            journeys.incrementCountMultithreaded(readIdValue, v.size());
        }
    }
}


void BidirectedAnchors::computeJourneysThreadFunction4(uint64_t /* threadId */)
{
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t readIdValue = begin; readIdValue < end; ++readIdValue) {
            const auto v = journeysWithOrdinals[readIdValue];
            const auto jrn = journeys[readIdValue];
            DINARA_ASSERT(jrn.size() == v.size());
            for(uint64_t i = 0; i < v.size(); ++i) {
                jrn[i] = BidirectedJourneyEntry(v[i].anchorId, v[i].strand);
            }
        }
    }
}


span<const BidirectedJourneyEntry> BidirectedAnchors::journey(ReadId readId) const
{
    return journeys[uint64_t(readId)];
}


void BidirectedAnchors::writeJourneys() const
{
    const ReadId readCount = reads.readCount();

    // Find the maximum journey length to determine the number of columns.
    uint64_t maxJourneyLength = 0;
    for(ReadId r = 0; r < readCount; ++r) {
        const auto jrn = journeys[uint64_t(r)];
        if(jrn.size() > maxJourneyLength) {
            maxJourneyLength = jrn.size();
        }
    }

    ofstream csv("BidirectedJourneys.csv");

    // Header: ReadId, then one column per step.
    // Values use Verkko notation: >anchorId (strand 0) or <anchorId (strand 1).
    csv << "ReadId";
    for(uint64_t i = 0; i < maxJourneyLength; i++) {
        csv << ",Step_" << i;
    }
    csv << "\n";

    // One row per read.
    for(ReadId r = 0; r < readCount; ++r) {
        const auto jrn = journeys[uint64_t(r)];
        if(jrn.empty()) continue;
        csv << r;
        for(uint64_t i = 0; i < jrn.size(); i++) {
            csv << "," << (jrn[i].strand == 0 ? ">" : "<") << jrn[i].anchorId;
        }
        csv << "\n";
    }
}


// ============================================================================
// Edge table computation.
//
// For every consecutive pair of journey entries (A,sA) → (B,sB), record
// a directed edge.  Enforce the bidirected symmetry invariant:
//   edge (A,sA)→(B,sB)  implies  edge reverse(B,sB)→reverse(A,sA)
// with the same coverage.
//
// Uses canonical edge representation to accumulate coverage, then
// stores edges in per-directed-node adjacency lists for O(1) lookup.
// ============================================================================

void BidirectedAnchors::computeEdges(uint64_t threadCount, uint64_t minEdgeCoverage)
{
    performanceLog << timestamp << "BidirectedAnchors::computeEdges begins."
        << " minEdgeCoverage=" << minEdgeCoverage << endl;

    (void)threadCount;  // TODO: parallelize if needed.

    const uint64_t readCount = reads.readCount();

    // Phase 1: Count edges using canonical representation.
    // canonical edge → total coverage.
    map<pair<OrientedBidirectedAnchor, OrientedBidirectedAnchor>, uint64_t> edgeCoverage;

    for(uint64_t readIdValue = 0; readIdValue < readCount; ++readIdValue) {
        const auto jrn = journeys[readIdValue];
        if(jrn.size() < 2) continue;

        for(uint64_t i = 0; i + 1 < jrn.size(); ++i) {
            OrientedBidirectedAnchor from(jrn[i].anchorId, jrn[i].strand);
            OrientedBidirectedAnchor to(jrn[i + 1].anchorId, jrn[i + 1].strand);

            // Always store in canonical form.
            auto canonEdge = canon(from, to);
            ++edgeCoverage[canonEdge];
        }
    }

    // Phase 2: Expand canonical edges into both directions and populate
    // the adjacency lists.  Edges below minEdgeCoverage are discarded.
    outEdges.resize(anchorCount);
    inEdges.resize(anchorCount);

    uint64_t removedEdgeCount = 0;

    for(const auto& [canonPair, cov] : edgeCoverage) {

        // Filter out low-coverage edges.
        if(cov < minEdgeCoverage) {
            ++removedEdgeCount;
            continue;
        }

        const OrientedBidirectedAnchor& from = canonPair.first;
        const OrientedBidirectedAnchor& to = canonPair.second;

        // Forward direction: from → to.
        outEdges[from].push_back(BidirectedEdge(from, to, cov));
        inEdges[to].push_back(BidirectedEdge(from, to, cov));

        // Reverse direction: reverse(to) → reverse(from).
        OrientedBidirectedAnchor revTo = reverse(to);
        OrientedBidirectedAnchor revFrom = reverse(from);
        if(revTo != from || revFrom != to) {
            // Only add if it's a different edge (avoid double-counting palindromic edges).
            outEdges[revTo].push_back(BidirectedEdge(revTo, revFrom, cov));
            inEdges[revFrom].push_back(BidirectedEdge(revTo, revFrom, cov));
        }
    }

    // Sort each adjacency list for deterministic output.
    for(uint64_t i = 0; i < anchorCount; ++i) {
        for(Strand s = 0; s <= 1; ++s) {
            OrientedBidirectedAnchor oa(i, s);
            auto& out = outEdges[oa];
            sort(out.begin(), out.end());
            auto& in = inEdges[oa];
            sort(in.begin(), in.end());
        }
    }

    edgesComputed = true;

    // Count total edges.
    uint64_t totalEdges = 0;
    for(uint64_t i = 0; i < anchorCount; ++i) {
        totalEdges += outEdges[OrientedBidirectedAnchor(i, 0)].size();
        totalEdges += outEdges[OrientedBidirectedAnchor(i, 1)].size();
    }

    performanceLog << timestamp << "BidirectedAnchors::computeEdges ends. "
        << edgeCoverage.size() << " canonical edges, "
        << removedEdgeCount << " removed (coverage < " << minEdgeCoverage << "), "
        << totalEdges << " directed edges." << endl;

    writeEdges();
}


span<const BidirectedEdge> BidirectedAnchors::getOutEdges(OrientedBidirectedAnchor oa) const
{
    DINARA_ASSERT(edgesComputed);
    const auto& v = outEdges[oa];
    return span<const BidirectedEdge>(v.data(), v.size());
}


span<const BidirectedEdge> BidirectedAnchors::getInEdges(OrientedBidirectedAnchor oa) const
{
    DINARA_ASSERT(edgesComputed);
    const auto& v = inEdges[oa];
    return span<const BidirectedEdge>(v.data(), v.size());
}


void BidirectedAnchors::writeEdges() const
{
    ofstream csv("BidirectedEdges.csv");
    csv << "FromAnchor,FromStrand,ToAnchor,ToStrand,Coverage\n";
    for(uint64_t i = 0; i < anchorCount; ++i) {
        for(Strand s = 0; s <= 1; ++s) {
            const auto edges = getOutEdges(OrientedBidirectedAnchor(i, s));
            for(const auto& e : edges) {
                csv << e.from.anchorId << "," << uint32_t(e.from.strand)
                    << "," << e.to.anchorId << "," << uint32_t(e.to.strand)
                    << "," << e.coverage << "\n";
            }
        }
    }
}



// ============================================================================
// Transitive reduction.
//
// Performance-critical: uses flat vectors for BFS visited state and
// removed-edge tracking instead of std::map/std::set, giving O(1) lookups
// and avoiding per-BFS heap allocation.
// ============================================================================

bool BidirectedAnchors::transitiveReductionCanRemove(
    OrientedBidirectedAnchor from,
    OrientedBidirectedAnchor to,
    uint64_t edgeCoverage,
    uint64_t maxDistance,
    const set<pair<OrientedBidirectedAnchor, OrientedBidirectedAnchor>>& /*removedSet*/) const
{
    // Unused — kept for interface compatibility.
    // The optimized path uses the flat-vector overload below.
    (void)from; (void)to; (void)edgeCoverage; (void)maxDistance;
    return false;
}


void BidirectedAnchors::transitiveReduction(
    uint64_t maxEdgeCoverage,
    uint64_t maxDistance)
{
    DINARA_ASSERT(edgesComputed);

    performanceLog << timestamp << "BidirectedAnchors::transitiveReduction begins."
        << " maxEdgeCoverage=" << maxEdgeCoverage
        << " maxDistance=" << maxDistance << endl;

    const uint64_t orientedCount = 2 * anchorCount;

    // Helper: oriented anchor → flat index.
    auto oaIndex = [](OrientedBidirectedAnchor oa) -> uint64_t {
        return 2 * oa.anchorId + uint64_t(oa.strand);
    };

    // ----------------------------------------------------------------
    // Step 1: Assign each canonical edge a numeric ID and collect them
    //         sorted by coverage for efficient per-level iteration.
    // ----------------------------------------------------------------

    struct CanonicalEdgeInfo {
        OrientedBidirectedAnchor from;
        OrientedBidirectedAnchor to;
        uint64_t coverage;
        uint64_t canonId;
    };

    // Map canonical pair → ID.
    map<pair<OrientedBidirectedAnchor, OrientedBidirectedAnchor>, uint64_t> canonToId;
    vector<CanonicalEdgeInfo> canonEdges;

    for(BidirectedAnchorId i = 0; i < anchorCount; ++i) {
        for(Strand s = 0; s <= 1; ++s) {
            const OrientedBidirectedAnchor oa(i, s);
            for(const auto& e : getOutEdges(oa)) {
                auto ce = canon(e.from, e.to);
                if(!canonToId.contains(ce)) {
                    const uint64_t id = canonEdges.size();
                    canonToId[ce] = id;
                    canonEdges.push_back({ce.first, ce.second, e.coverage, id});
                }
            }
        }
    }

    const uint64_t totalCanonicalBefore = canonEdges.size();

    // Pre-map each directed edge to its canonical ID for O(1) lookup during BFS.
    // directedEdgeCanonId[oaIndex(from)][j] = canonId of the j-th out-edge from `from`.
    vector<vector<uint64_t>> directedEdgeCanonId(orientedCount);
    for(BidirectedAnchorId i = 0; i < anchorCount; ++i) {
        for(Strand s = 0; s <= 1; ++s) {
            const OrientedBidirectedAnchor oa(i, s);
            const auto edges = getOutEdges(oa);
            auto& ids = directedEdgeCanonId[oaIndex(oa)];
            ids.resize(edges.size());
            for(uint64_t j = 0; j < edges.size(); ++j) {
                ids[j] = canonToId[canon(edges[j].from, edges[j].to)];
            }
        }
    }

    // Sort canonical edges by coverage for per-level grouping.
    sort(canonEdges.begin(), canonEdges.end(),
        [](const CanonicalEdgeInfo& a, const CanonicalEdgeInfo& b) {
            return a.coverage < b.coverage;
        });

    // ----------------------------------------------------------------
    // Step 2: Flat removed-edge vector (O(1) lookup by canonId).
    // ----------------------------------------------------------------

    vector<bool> isRemoved(totalCanonicalBefore, false);

    // ----------------------------------------------------------------
    // Step 3: Reusable flat BFS state (avoids per-BFS allocation).
    //         bfsDist[oaIndex] stores the BFS distance, or INVALID if unvisited.
    //         bfsTouched collects indices to reset after each BFS.
    // ----------------------------------------------------------------

    const uint64_t UNVISITED = invalid<uint64_t>;
    vector<uint64_t> bfsDist(orientedCount, UNVISITED);
    vector<uint64_t> bfsTouched;
    bfsTouched.reserve(1024);

    // BFS queue — reused across calls.
    queue<OrientedBidirectedAnchor> q;

    // ----------------------------------------------------------------
    // Step 4: Process by coverage level.
    // ----------------------------------------------------------------

    uint64_t cursor = 0;  // into sorted canonEdges

    for(uint64_t cov = 1; cov <= maxEdgeCoverage; cov++) {

        // Advance cursor to find range [cursor, rangeEnd) with this coverage.
        const uint64_t rangeStart = cursor;
        while(cursor < totalCanonicalBefore && canonEdges[cursor].coverage == cov) {
            ++cursor;
        }
        const uint64_t rangeEnd = cursor;

        if(rangeStart == rangeEnd) {
            continue;
        }

        uint64_t removedAtThisCoverage = 0;

        for(uint64_t idx = rangeStart; idx < rangeEnd; ++idx) {
            const auto& ce = canonEdges[idx];
            const OrientedBidirectedAnchor from = ce.from;
            const OrientedBidirectedAnchor to = ce.to;

            // ---- BFS from `from` looking for `to` ----

            bfsDist[oaIndex(from)] = 0;
            bfsTouched.push_back(oaIndex(from));
            q.push(from);

            bool found = false;

            while(!q.empty()) {
                const OrientedBidirectedAnchor cur = q.front();
                q.pop();

                const uint64_t curDist = bfsDist[oaIndex(cur)];
                if(curDist >= maxDistance) {
                    continue;
                }

                const auto edges = getOutEdges(cur);
                const auto& edgeCanonIds = directedEdgeCanonId[oaIndex(cur)];

                for(uint64_t j = 0; j < edges.size(); ++j) {
                    const auto& e = edges[j];

                    // Skip the candidate edge itself.
                    if(cur == from && e.to == to) {
                        continue;
                    }

                    // Only use edges with strictly higher coverage.
                    if(e.coverage <= cov) {
                        continue;
                    }

                    // Skip already-removed edges (O(1) lookup).
                    if(isRemoved[edgeCanonIds[j]]) {
                        continue;
                    }

                    const OrientedBidirectedAnchor next = e.to;
                    const uint64_t nextIdx = oaIndex(next);

                    // Found the target via an alternate path.
                    if(next == to) {
                        found = true;
                        goto bfs_done;
                    }

                    // Already visited.
                    if(bfsDist[nextIdx] != UNVISITED) {
                        continue;
                    }

                    bfsDist[nextIdx] = curDist + 1;
                    bfsTouched.push_back(nextIdx);
                    q.push(next);
                }
            }

            bfs_done:

            // Clean up BFS state for reuse.
            for(const uint64_t ti : bfsTouched) {
                bfsDist[ti] = UNVISITED;
            }
            bfsTouched.clear();

            // Drain the queue (may have leftover entries after goto).
            while(!q.empty()) {
                q.pop();
            }

            if(found) {
                isRemoved[ce.canonId] = true;
                ++removedAtThisCoverage;
            }
        }

        performanceLog << "  Coverage " << cov
            << ": " << (rangeEnd - rangeStart) << " canonical edges, "
            << removedAtThisCoverage << " removed." << endl;
    }

    // ----------------------------------------------------------------
    // Step 5: Rebuild adjacency lists excluding removed edges.
    // ----------------------------------------------------------------

    uint64_t directedBefore = 0;
    uint64_t directedAfter = 0;

    for(BidirectedAnchorId i = 0; i < anchorCount; ++i) {
        for(Strand s = 0; s <= 1; ++s) {
            const OrientedBidirectedAnchor oa(i, s);
            const uint64_t oi = oaIndex(oa);

            auto& out = outEdges[oa];
            const auto& outCanonIds = directedEdgeCanonId[oi];
            directedBefore += out.size();

            // In-place removal using pre-computed canon IDs.
            uint64_t writePos = 0;
            for(uint64_t j = 0; j < out.size(); ++j) {
                if(!isRemoved[outCanonIds[j]]) {
                    if(writePos != j) {
                        out[writePos] = out[j];
                    }
                    ++writePos;
                }
            }
            out.resize(writePos);
            directedAfter += writePos;

            // Rebuild inEdges similarly — but inEdges don't have pre-computed
            // canon IDs, so compute on the fly (done once, not in BFS).
            auto& in = inEdges[oa];
            uint64_t inWritePos = 0;
            for(uint64_t j = 0; j < in.size(); ++j) {
                auto ce = canon(in[j].from, in[j].to);
                auto it = canonToId.find(ce);
                if(it == canonToId.end() || !isRemoved[it->second]) {
                    if(inWritePos != j) {
                        in[inWritePos] = in[j];
                    }
                    ++inWritePos;
                }
            }
            in.resize(inWritePos);
        }
    }

    uint64_t removedCount = 0;
    for(uint64_t i = 0; i < totalCanonicalBefore; ++i) {
        if(isRemoved[i]) ++removedCount;
    }

    performanceLog << timestamp << "BidirectedAnchors::transitiveReduction ends. "
        << removedCount << " of " << totalCanonicalBefore
        << " canonical edges removed. "
        << directedBefore << " directed edges before, "
        << directedAfter << " after." << endl;
}



void BidirectedAnchors::writeGfa(const string& fileName, bool includePaths) const
{
    if(!edgesComputed) {
        throw runtime_error(
            "BidirectedAnchors::writeGfa called before edges are computed.");
    }

    ofstream gfa(fileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }

    gfa << "H\tVN:Z:1.0\n";

    // S lines: one segment per anchor.
    for(BidirectedAnchorId anchorId = 0; anchorId < anchorCount; ++anchorId) {
        gfa << "S\t" << anchorId
            << "\t*"
            << "\tLN:i:" << max<uint64_t>(k, 1)
            << "\tRC:i:" << coverage(anchorId)
            << "\n";
    }

    // L lines: one per canonical oriented edge.
    map<pair<OrientedBidirectedAnchor, OrientedBidirectedAnchor>, uint64_t> canonicalEdges;
    for(BidirectedAnchorId anchorId = 0; anchorId < anchorCount; ++anchorId) {
        for(Strand strand = 0; strand <= 1; ++strand) {
            const OrientedBidirectedAnchor oa(anchorId, strand);
            for(const BidirectedEdge& edge : getOutEdges(oa)) {
                const auto canonicalEdge = canon(edge.from, edge.to);
                auto it = canonicalEdges.find(canonicalEdge);
                if(it == canonicalEdges.end()) {
                    canonicalEdges.insert({canonicalEdge, edge.coverage});
                } else {
                    // Symmetric entries should have equal coverage.
                    // Keep max for robustness.
                    it->second = max(it->second, edge.coverage);
                }
            }
        }
    }

    for(const auto& [edge, cov] : canonicalEdges) {
        const OrientedBidirectedAnchor from = edge.first;
        const OrientedBidirectedAnchor to = edge.second;
        gfa << "L\t"
            << from.anchorId << "\t" << (from.strand == 0 ? "+" : "-") << "\t"
            << to.anchorId << "\t" << (to.strand == 0 ? "+" : "-") << "\t"
            << "0M"
            << "\tRC:i:" << cov
            << "\n";
    }

    // Optional P lines: one path per physical read journey.
    if(includePaths) {
        const ReadId readCount = reads.readCount();
        for(ReadId readId = 0; readId < readCount; ++readId) {
            const auto jrn = journey(readId);
            if(jrn.empty()) {
                continue;
            }
            gfa << "P\tread-" << readId << "\t";
            for(uint64_t i = 0; i < jrn.size(); ++i) {
                if(i > 0) {
                    gfa << ",";
                }
                gfa << jrn[i].anchorId << (jrn[i].strand == 0 ? "+" : "-");
            }
            gfa << "\t*\n";
        }
    }

    performanceLog << timestamp
        << "BidirectedAnchors::writeGfa wrote " << fileName
        << " (" << anchorCount << " segments, "
        << canonicalEdges.size() << " links, "
        << (includePaths ? "with" : "without") << " paths)." << endl;
}


// ============================================================================
// Neighbor finding — uses the precomputed edge table.
// ============================================================================

void BidirectedAnchors::findForwardNeighbors(
    OrientedBidirectedAnchor oa,
    vector<Neighbor>& neighbors,
    uint64_t minCount) const
{
    neighbors.clear();

    if(edgesComputed) {
        // O(1) lookup from edge table.
        const auto edges = getOutEdges(oa);
        for(const auto& e : edges) {
            if(e.coverage >= minCount) {
                neighbors.push_back(Neighbor{e.to.anchorId, e.to.strand, e.coverage});
            }
        }
    } else {
        // Fallback: scan journeys (used before computeEdges is called).
        const auto intervals = anchorMarkerIntervals[oa.anchorId];
        for(const auto& interval : intervals) {
            if(interval.strand != oa.strand) continue;
            const auto jrn = journeys[uint64_t(interval.readId)];
            for(uint64_t pos = 0; pos < jrn.size(); ++pos) {
                if(jrn[pos].anchorId == oa.anchorId) {
                    if(pos + 1 < jrn.size()) {
                        neighbors.push_back(Neighbor{
                            jrn[pos + 1].anchorId,
                            jrn[pos + 1].strand, 1});
                    }
                    break;
                }
            }
        }
        // Deduplicate and count.
        sort(neighbors.begin(), neighbors.end(),
            [](const Neighbor& a, const Neighbor& b) {
                if(a.anchorId != b.anchorId) return a.anchorId < b.anchorId;
                return a.strand < b.strand;
            });
        vector<Neighbor> merged;
        for(size_t i = 0; i < neighbors.size(); ) {
            size_t j = i + 1;
            while(j < neighbors.size() &&
                  neighbors[j].anchorId == neighbors[i].anchorId &&
                  neighbors[j].strand == neighbors[i].strand) {
                ++j;
            }
            const uint64_t cnt = j - i;
            if(cnt >= minCount) {
                merged.push_back(Neighbor{neighbors[i].anchorId, neighbors[i].strand, cnt});
            }
            i = j;
        }
        neighbors = std::move(merged);
    }
}


void BidirectedAnchors::findBackwardNeighbors(
    OrientedBidirectedAnchor oa,
    vector<Neighbor>& neighbors,
    uint64_t minCount) const
{
    neighbors.clear();

    if(edgesComputed) {
        // O(1) lookup from edge table.
        const auto edges = getInEdges(oa);
        for(const auto& e : edges) {
            if(e.coverage >= minCount) {
                neighbors.push_back(Neighbor{e.from.anchorId, e.from.strand, e.coverage});
            }
        }
    } else {
        // Fallback: scan journeys.
        const auto intervals = anchorMarkerIntervals[oa.anchorId];
        for(const auto& interval : intervals) {
            if(interval.strand != oa.strand) continue;
            const auto jrn = journeys[uint64_t(interval.readId)];
            for(uint64_t pos = 0; pos < jrn.size(); ++pos) {
                if(jrn[pos].anchorId == oa.anchorId) {
                    if(pos > 0) {
                        neighbors.push_back(Neighbor{
                            jrn[pos - 1].anchorId,
                            jrn[pos - 1].strand, 1});
                    }
                    break;
                }
            }
        }
        sort(neighbors.begin(), neighbors.end(),
            [](const Neighbor& a, const Neighbor& b) {
                if(a.anchorId != b.anchorId) return a.anchorId < b.anchorId;
                return a.strand < b.strand;
            });
        vector<Neighbor> merged;
        for(size_t i = 0; i < neighbors.size(); ) {
            size_t j = i + 1;
            while(j < neighbors.size() &&
                  neighbors[j].anchorId == neighbors[i].anchorId &&
                  neighbors[j].strand == neighbors[i].strand) {
                ++j;
            }
            const uint64_t cnt = j - i;
            if(cnt >= minCount) {
                merged.push_back(Neighbor{neighbors[i].anchorId, neighbors[i].strand, cnt});
            }
            i = j;
        }
        neighbors = std::move(merged);
    }
}


// ============================================================================
// Common read counting.
// ============================================================================

uint64_t BidirectedAnchors::countCommon(BidirectedAnchorId a, BidirectedAnchorId b) const
{
    const auto intervalsA = anchorMarkerIntervals[a];
    const auto intervalsB = anchorMarkerIntervals[b];

    // Both are sorted by ReadId.
    uint64_t count = 0;
    size_t ia = 0, ib = 0;
    while(ia < intervalsA.size() && ib < intervalsB.size()) {
        if(intervalsA[ia].readId < intervalsB[ib].readId) {
            ++ia;
        } else if(intervalsA[ia].readId > intervalsB[ib].readId) {
            ++ib;
        } else {
            ++count;
            ++ia;
            ++ib;
        }
    }
    return count;
}



// ============================================================================
// Unitigification.
//
// Merge maximal linear chains into single unitig anchors.
// A maximal chain is a path of oriented anchors where every internal
// transition is the only out-edge of its source and the only in-edge
// of its target (out-degree 1, in-degree 1).
//
// After unitigification the graph is smaller: the number of anchors
// equals the number of unitigs.  Journeys, marker intervals, and the
// edge table are rebuilt to use unitig IDs.
// ============================================================================

void BidirectedAnchors::unitigifyAll()
{
    DINARA_ASSERT(edgesComputed);

    performanceLog << timestamp << "BidirectedAnchors::unitigifyAll begins." << endl;

    const uint64_t oldAnchorCount = anchorCount;
    const uint64_t readCount = reads.readCount();

    // ----------------------------------------------------------------
    // Step 1: Find maximal unitig chains.
    // ----------------------------------------------------------------

    vector<bool> visited(oldAnchorCount, false);
    vector<vector<OrientedBidirectedAnchor>> chains;
    chains.reserve(oldAnchorCount);  // upper bound; shrinks after

    for(BidirectedAnchorId anchorId = 0; anchorId < oldAnchorCount; ++anchorId) {
        if(visited[anchorId]) continue;
        visited[anchorId] = true;

        OrientedBidirectedAnchor start(anchorId, 0);

        // Extend forward.
        vector<OrientedBidirectedAnchor> fwd;
        fwd.push_back(start);
        {
            OrientedBidirectedAnchor cur = start;
            while(true) {
                const auto outgoing = getOutEdges(cur);
                if(outgoing.size() != 1) break;
                const OrientedBidirectedAnchor next = outgoing[0].to;
                if(visited[next.anchorId]) break;
                const auto incoming = getInEdges(next);
                if(incoming.size() != 1) break;
                fwd.push_back(next);
                visited[next.anchorId] = true;
                cur = next;
            }
        }

        // Extend backward.
        vector<OrientedBidirectedAnchor> bwd;
        {
            OrientedBidirectedAnchor cur = start;
            while(true) {
                const auto incoming = getInEdges(cur);
                if(incoming.size() != 1) break;
                const OrientedBidirectedAnchor prev = incoming[0].from;
                if(visited[prev.anchorId]) break;
                const auto outgoing = getOutEdges(prev);
                if(outgoing.size() != 1) break;
                bwd.push_back(prev);
                visited[prev.anchorId] = true;
                cur = prev;
            }
        }

        // Combine: reverse(backward) + forward.
        vector<OrientedBidirectedAnchor> chain;
        chain.reserve(bwd.size() + fwd.size());
        for(auto it = bwd.rbegin(); it != bwd.rend(); ++it) {
            chain.push_back(*it);
        }
        for(const auto& oa : fwd) {
            chain.push_back(oa);
        }

        chains.push_back(std::move(chain));
    }
    chains.shrink_to_fit();

    const uint64_t newAnchorCount = chains.size();

    // Statistics.
    uint64_t singletonCount = 0;
    uint64_t maxChainLen = 0;
    uint64_t totalChainAnchors = 0;
    for(const auto& chain : chains) {
        if(chain.size() == 1) ++singletonCount;
        maxChainLen = max(maxChainLen, chain.size());
        totalChainAnchors += chain.size();
    }

    cout << timestamp << "Unitigification: " << oldAnchorCount << " anchors -> "
         << newAnchorCount << " unitigs ("
         << singletonCount << " singletons, longest chain " << maxChainLen
         << ", avg " << (double(totalChainAnchors) / max(newAnchorCount, uint64_t(1)))
         << ")." << endl;

    if(newAnchorCount == oldAnchorCount) {
        cout << timestamp << "No chains to merge; graph is already unitigified." << endl;
        unitigChains.resize(newAnchorCount);
        for(BidirectedAnchorId i = 0; i < newAnchorCount; ++i) {
            unitigChains[i] = { OrientedBidirectedAnchor(i, 0) };
        }
        performanceLog << timestamp << "BidirectedAnchors::unitigifyAll ends (no-op)." << endl;
        return;
    }

    // ----------------------------------------------------------------
    // Step 2: Build old-anchor -> unitig mapping.
    // ----------------------------------------------------------------

    struct AnchorMapping {
        uint64_t unitigId;
        Strand strandInChain;  // strand of this anchor within its chain
    };
    vector<AnchorMapping> anchorToUnitig(oldAnchorCount);

    for(uint64_t unitigId = 0; unitigId < newAnchorCount; ++unitigId) {
        for(const auto& oa : chains[unitigId]) {
            anchorToUnitig[oa.anchorId] = {unitigId, oa.strand};
        }
    }

    // ----------------------------------------------------------------
    // Step 3: Compute unitig-level edges directly from the old edge table.
    // Only boundary edges (between different unitigs) survive.
    // This avoids rescanning all journey entries via computeEdges().
    //
    // Collect into a flat vector, sort, and merge -- much faster than
    // std::map for cache locality on large edge sets.
    // ----------------------------------------------------------------

    struct RawUnitigEdge {
        uint64_t fromEncoded;  // (unitigId << 1) | strand
        uint64_t toEncoded;
        uint64_t coverage;
    };
    vector<RawUnitigEdge> rawEdges;
    rawEdges.reserve(oldAnchorCount);  // rough estimate

    for(uint64_t oldAnchorId = 0; oldAnchorId < oldAnchorCount; ++oldAnchorId) {
        for(Strand s = 0; s <= 1; ++s) {
            OrientedBidirectedAnchor oldOa(oldAnchorId, s);
            for(const auto& edge : getOutEdges(oldOa)) {
                const auto& fromMap = anchorToUnitig[edge.from.anchorId];
                const auto& toMap = anchorToUnitig[edge.to.anchorId];

                // Skip internal edges (both endpoints in the same unitig).
                if(fromMap.unitigId == toMap.unitigId) continue;

                // Compute unitig-level orientations.
                Strand fromStrand = (edge.from.strand == fromMap.strandInChain) ? 0 : 1;
                Strand toStrand = (edge.to.strand == toMap.strandInChain) ? 0 : 1;

                // Canonicalize for deduplication.
                OrientedBidirectedAnchor uf(fromMap.unitigId, fromStrand);
                OrientedBidirectedAnchor ut(toMap.unitigId, toStrand);
                auto [cf, ct] = canon(uf, ut);

                rawEdges.push_back({
                    (cf.anchorId << 1) | cf.strand,
                    (ct.anchorId << 1) | ct.strand,
                    edge.coverage});
            }
        }
    }

    // Sort by (from, to) so duplicates are adjacent, then merge.
    sort(rawEdges.begin(), rawEdges.end(),
        [](const RawUnitigEdge& a, const RawUnitigEdge& b) {
            if(a.fromEncoded != b.fromEncoded) return a.fromEncoded < b.fromEncoded;
            return a.toEncoded < b.toEncoded;
        });

    // Deduplicate: keep max coverage for each canonical edge.
    vector<RawUnitigEdge> mergedEdges;
    mergedEdges.reserve(rawEdges.size());
    for(size_t i = 0; i < rawEdges.size(); ) {
        size_t j = i + 1;
        uint64_t maxCov = rawEdges[i].coverage;
        while(j < rawEdges.size() &&
              rawEdges[j].fromEncoded == rawEdges[i].fromEncoded &&
              rawEdges[j].toEncoded == rawEdges[i].toEncoded) {
            maxCov = max(maxCov, rawEdges[j].coverage);
            ++j;
        }
        mergedEdges.push_back({rawEdges[i].fromEncoded, rawEdges[i].toEncoded, maxCov});
        i = j;
    }
    rawEdges.clear();
    rawEdges.shrink_to_fit();

    // ----------------------------------------------------------------
    // Step 4: Build new marker intervals.
    // Singleton unitigs (majority): just copy with strand adjustment.
    // Multi-anchor chains: collect, sort by ReadId, deduplicate.
    // ----------------------------------------------------------------

    vector<vector<BidirectedAnchorMarkerInterval>> newIntervals(newAnchorCount);

    for(uint64_t unitigId = 0; unitigId < newAnchorCount; ++unitigId) {
        const auto& chainAnchors = chains[unitigId];

        if(chainAnchors.size() == 1) {
            // Fast path for singletons (typically 60-80% of all unitigs).
            const auto& oa = chainAnchors[0];
            const auto intervals = anchorMarkerIntervals[oa.anchorId];
            newIntervals[unitigId].reserve(intervals.size());
            if(oa.strand == 0) {
                // No strand flip needed — direct copy.
                for(const auto& interval : intervals) {
                    newIntervals[unitigId].push_back(interval);
                }
            } else {
                for(const auto& interval : intervals) {
                    newIntervals[unitigId].push_back(BidirectedAnchorMarkerInterval(
                        interval.readId, interval.ordinal, Strand(1 - interval.strand)));
                }
            }
        } else {
            // Multi-anchor chain: collect all, sort, deduplicate.
            vector<BidirectedAnchorMarkerInterval> all;
            for(const auto& oa : chainAnchors) {
                const auto intervals = anchorMarkerIntervals[oa.anchorId];
                for(const auto& interval : intervals) {
                    Strand adjustedStrand = (oa.strand == 0)
                        ? interval.strand
                        : Strand(1 - interval.strand);
                    all.push_back(BidirectedAnchorMarkerInterval(
                        interval.readId, interval.ordinal, adjustedStrand));
                }
            }

            sort(all.begin(), all.end(),
                [](const BidirectedAnchorMarkerInterval& a, const BidirectedAnchorMarkerInterval& b) {
                    return a.readId < b.readId;
                });

            newIntervals[unitigId].reserve(all.size());
            for(size_t i = 0; i < all.size(); ) {
                newIntervals[unitigId].push_back(all[i]);
                const ReadId r = all[i].readId;
                while(i < all.size() && all[i].readId == r) ++i;
            }
        }
    }

    // ----------------------------------------------------------------
    // Step 5: Build new journeys.
    // Two-pass: count entries first, then fill.
    // Avoids allocating millions of small std::vector objects.
    // ----------------------------------------------------------------

    // Pass 1: count collapsed journey entries per read.
    vector<uint64_t> journeySizes(readCount, 0);

    for(uint64_t readIdValue = 0; readIdValue < readCount; ++readIdValue) {
        const auto jrn = journeys[readIdValue];
        if(jrn.empty()) continue;

        uint64_t lastUnitigId = invalid<uint64_t>;
        Strand lastStrand = 0;

        for(uint64_t i = 0; i < jrn.size(); ++i) {
            const auto& mapping = anchorToUnitig[jrn[i].anchorId];
            const uint64_t unitigId = mapping.unitigId;
            const Strand unitigStrand =
                (jrn[i].strand == mapping.strandInChain) ? Strand(0) : Strand(1);

            if(unitigId != lastUnitigId || unitigStrand != lastStrand) {
                ++journeySizes[readIdValue];
                lastUnitigId = unitigId;
                lastStrand = unitigStrand;
            }
        }
    }

    // Pass 2: fill a flat buffer.
    // Compute offsets from sizes.
    vector<uint64_t> journeyOffsets(readCount + 1, 0);
    for(uint64_t i = 0; i < readCount; ++i) {
        journeyOffsets[i + 1] = journeyOffsets[i] + journeySizes[i];
    }
    const uint64_t totalJourneyEntries = journeyOffsets[readCount];
    vector<BidirectedJourneyEntry> flatJourneys(totalJourneyEntries);

    for(uint64_t readIdValue = 0; readIdValue < readCount; ++readIdValue) {
        const auto jrn = journeys[readIdValue];
        if(jrn.empty()) continue;

        uint64_t pos = journeyOffsets[readIdValue];
        uint64_t lastUnitigId = invalid<uint64_t>;
        Strand lastStrand = 0;

        for(uint64_t i = 0; i < jrn.size(); ++i) {
            const auto& mapping = anchorToUnitig[jrn[i].anchorId];
            const uint64_t unitigId = mapping.unitigId;
            const Strand unitigStrand =
                (jrn[i].strand == mapping.strandInChain) ? Strand(0) : Strand(1);

            if(unitigId != lastUnitigId || unitigStrand != lastStrand) {
                flatJourneys[pos++] = BidirectedJourneyEntry(unitigId, unitigStrand);
                lastUnitigId = unitigId;
                lastStrand = unitigStrand;
            }
        }
    }

    // ----------------------------------------------------------------
    // Step 6: Replace internal data structures.
    // ----------------------------------------------------------------

    // Store chain composition for traceability.
    unitigChains = std::move(chains);
    anchorCount = newAnchorCount;

    // Rebuild marker intervals.
    anchorMarkerIntervals.remove();
    anchorMarkerIntervals.createNew(
        largeDataName("BidirectedAnchorMarkerIntervals"), largeDataPageSize);
    for(const auto& intervals : newIntervals) {
        anchorMarkerIntervals.appendVector();
        for(const auto& interval : intervals) {
            anchorMarkerIntervals.append(interval);
        }
    }
    newIntervals.clear();
    newIntervals.shrink_to_fit();

    // Rebuild anchor infos.
    anchorInfos.remove();
    anchorInfos.createNew(
        largeDataName("BidirectedAnchorInfos"), largeDataPageSize);
    anchorInfos.resize(newAnchorCount);
    for(uint64_t i = 0; i < newAnchorCount; ++i) {
        anchorInfos[i].componentId = invalid<uint32_t>;
        anchorInfos[i].localAnchorIdInComponent = invalid<uint64_t>;
    }

    // Rebuild journeys from flat buffer.
    journeys.remove();
    journeys.createNew(
        largeDataName("BidirectedJourneys"), largeDataPageSize);
    for(uint64_t readIdValue = 0; readIdValue < readCount; ++readIdValue) {
        journeys.appendVector();
        const uint64_t begin = journeyOffsets[readIdValue];
        const uint64_t end = journeyOffsets[readIdValue + 1];
        for(uint64_t j = begin; j < end; ++j) {
            journeys.append(flatJourneys[j]);
        }
    }
    flatJourneys.clear();
    flatJourneys.shrink_to_fit();
    journeyOffsets.clear();
    journeyOffsets.shrink_to_fit();
    journeySizes.clear();
    journeySizes.shrink_to_fit();

    // ----------------------------------------------------------------
    // Step 7: Build unitig edge table directly from precomputed edges.
    // No need to rescan journeys — we computed boundary edges in Step 3.
    // ----------------------------------------------------------------

    outEdges.clear();
    inEdges.clear();
    outEdges.resize(newAnchorCount);
    inEdges.resize(newAnchorCount);

    for(const auto& me : mergedEdges) {
        OrientedBidirectedAnchor from(me.fromEncoded >> 1, Strand(me.fromEncoded & 1));
        OrientedBidirectedAnchor to(me.toEncoded >> 1, Strand(me.toEncoded & 1));

        // Forward direction.
        outEdges[from].push_back(BidirectedEdge(from, to, me.coverage));
        inEdges[to].push_back(BidirectedEdge(from, to, me.coverage));

        // Reverse direction.
        OrientedBidirectedAnchor revTo = reverse(to);
        OrientedBidirectedAnchor revFrom = reverse(from);
        if(revTo != from || revFrom != to) {
            outEdges[revTo].push_back(BidirectedEdge(revTo, revFrom, me.coverage));
            inEdges[revFrom].push_back(BidirectedEdge(revTo, revFrom, me.coverage));
        }
    }
    mergedEdges.clear();
    mergedEdges.shrink_to_fit();

    // Sort adjacency lists for deterministic output.
    for(uint64_t i = 0; i < newAnchorCount; ++i) {
        for(Strand s = 0; s <= 1; ++s) {
            OrientedBidirectedAnchor oa(i, s);
            sort(outEdges[oa].begin(), outEdges[oa].end());
            sort(inEdges[oa].begin(), inEdges[oa].end());
        }
    }
    edgesComputed = true;

    // Count total edges for logging.
    uint64_t totalEdges = 0;
    for(uint64_t i = 0; i < newAnchorCount; ++i) {
        totalEdges += outEdges[OrientedBidirectedAnchor(i, 0)].size();
        totalEdges += outEdges[OrientedBidirectedAnchor(i, 1)].size();
    }

    performanceLog << timestamp << "BidirectedAnchors::unitigifyAll ends. "
        << newAnchorCount << " unitigs, "
        << totalEdges << " directed edges." << endl;
}
