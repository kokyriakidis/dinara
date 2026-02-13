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

void BidirectedAnchors::computeEdges(uint64_t threadCount)
{
    performanceLog << timestamp << "BidirectedAnchors::computeEdges begins." << endl;

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
    // the adjacency lists.
    outEdges.resize(anchorCount);
    inEdges.resize(anchorCount);

    for(const auto& [canonPair, cov] : edgeCoverage) {
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
