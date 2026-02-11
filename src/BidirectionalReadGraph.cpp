// Dinara.
#include "BidirectionalReadGraph.hpp"
using namespace dinara;



const uint32_t BidirectionalReadGraph::infiniteDistance =
    std::numeric_limits<uint32_t>::max();



// ============================================================================
// buildConnectivity
// ============================================================================
// After all edges have been pushed into edges[], this builds the adjacency
// index.  Each ReadId gets a list of incident edge indices.
//
// Algorithm:
//   1. Count the degree of each read (two increments per edge: one per endpoint).
//   2. Use VectorOfVectors::beginPass / incrementCount / endPass to build
//      the compressed sparse structure, then scatter edge indices.
//
// Deleted edges are included in connectivity — callers skip them during traversal.
// This avoids rebuilding the index after soft-deletes.
void BidirectionalReadGraph::buildConnectivity(uint64_t nReads)
{
    // Pass 1: count degrees.
    connectivity.beginPass1(nReads);
    for(uint64_t edgeId = 0; edgeId < edges.size(); ++edgeId) {
        const auto& edge = edges[edgeId];
        connectivity.incrementCount(edge.readIds[0]);
        connectivity.incrementCount(edge.readIds[1]);
    }
    connectivity.beginPass2();

    // Pass 2: scatter edge indices.
    for(uint64_t edgeId = 0; edgeId < edges.size(); ++edgeId) {
        const auto& edge = edges[edgeId];
        connectivity.store(edge.readIds[0], uint32_t(edgeId));
        connectivity.store(edge.readIds[1], uint32_t(edgeId));
    }
    connectivity.endPass2();
}



// ============================================================================
// findNeighbors (distance 1)
// ============================================================================
void BidirectionalReadGraph::findNeighbors(
    ReadId readId,
    vector<ReadId>& neighbors) const
{
    neighbors.clear();
    for(const uint32_t edgeId : connectivity[readId]) {
        const auto& edge = edges[edgeId];
        if(edge.isDeleted) {
            continue;
        }
        neighbors.push_back(edge.getOther(readId));
    }
    sort(neighbors.begin(), neighbors.end());
}



// ============================================================================
// findNeighbors (multi-hop BFS)
// ============================================================================
void BidirectionalReadGraph::findNeighbors(
    ReadId readId,
    uint64_t maxDistance,
    vector<ReadId>& neighbors) const
{
    std::queue<ReadId> q;
    std::map<ReadId, uint64_t> distanceMap;
    q.push(readId);
    distanceMap[readId] = 0;

    neighbors.clear();
    while(!q.empty()) {
        const ReadId r0 = q.front();
        q.pop();
        const uint64_t d0 = distanceMap[r0];
        const uint64_t d1 = d0 + 1;
        DINARA_ASSERT(d1 <= maxDistance);

        for(const uint32_t edgeId : connectivity[r0]) {
            const auto& edge = edges[edgeId];
            if(edge.isDeleted) {
                continue;
            }
            const ReadId r1 = edge.getOther(r0);
            if(distanceMap.find(r1) != distanceMap.end()) {
                continue;
            }
            neighbors.push_back(r1);
            distanceMap[r1] = d1;
            if(d1 < maxDistance) {
                q.push(r1);
            }
        }
    }
    sort(neighbors.begin(), neighbors.end());
}



// ============================================================================
// computeShortPath
// ============================================================================
// BFS shortest path between two reads, skipping deleted edges.
// Work arrays must be pre-sized to readCount() and initialised to
// infiniteDistance / invalidReadId by the caller.
void BidirectionalReadGraph::computeShortPath(
    ReadId readId0,
    ReadId readId1,
    size_t maxDistance,
    vector<uint32_t>& path,
    vector<uint32_t>& distance,
    vector<ReadId>& reachedVertices,
    vector<uint32_t>& parentEdges) const
{
    path.clear();
    std::queue<ReadId> q;
    q.push(readId0);
    distance[readId0] = 0;
    reachedVertices.clear();
    reachedVertices.push_back(readId0);

    while(!q.empty()) {
        const ReadId v0 = q.front();
        q.pop();
        const uint32_t d0 = distance[v0];
        const uint32_t d1 = d0 + 1;

        bool found = false;
        for(const uint32_t edgeId : connectivity[v0]) {
            const auto& edge = edges[edgeId];
            if(edge.isDeleted) {
                continue;
            }
            const ReadId v1 = edge.getOther(v0);

            if(distance[v1] == infiniteDistance) {
                distance[v1] = d1;
                reachedVertices.push_back(v1);
                parentEdges[v1] = edgeId;
                if(d1 < maxDistance) {
                    q.push(v1);
                }
            }

            if(v1 == readId1) {
                found = true;
                ReadId v = v1;
                while(v != readId0) {
                    const uint32_t eid = parentEdges[v];
                    path.push_back(eid);
                    v = edges[eid].getOther(v);
                }
                std::reverse(path.begin(), path.end());
                break;
            }
        }
        if(found) {
            break;
        }
    }

    // Clean up work areas.
    for(const ReadId v : reachedVertices) {
        distance[v] = infiniteDistance;
    }
    reachedVertices.clear();
}



// ============================================================================
// propagateStrands  —  orientation-aware BFS
// ============================================================================
// Starting from (seedReadId, seedStrand), performs BFS over non-deleted edges.
// Each edge propagates derived strand:
//    toStrand = fromStrand                    for same-strand edges
//    toStrand = fromStrand ^ 1                for cross-strand edges
//
// If a read is reached a second time with a DIFFERENT implied strand, that
// constitutes a "conflict" — the read sits at an inversion or palindrome
// boundary.  The first-discovered strand wins (consistent with legacy
// computeReadGraphStrandsFromSeed behaviour).  The total conflict count is
// returned so callers can detect problematic regions.
BidirectionalReadGraph::StrandResult
BidirectionalReadGraph::propagateStrands(
    ReadId seedReadId,
    Strand seedStrand) const
{
    const uint64_t n = readCount();
    StrandResult result;
    result.strandByRead.assign(n, int8_t(-1));
    result.conflictCount = 0;

    result.strandByRead[seedReadId] = int8_t(seedStrand);

    std::queue<ReadId> q;
    q.push(seedReadId);

    while(!q.empty()) {
        const ReadId r0 = q.front();
        q.pop();
        const Strand s0 = Strand(result.strandByRead[r0]);

        for(const uint32_t edgeId : connectivity[r0]) {
            const auto& edge = edges[edgeId];
            if(edge.isDeleted) {
                continue;
            }
            const auto [r1, s1] = edge.traverse(r0, s0);

            if(result.strandByRead[r1] == int8_t(-1)) {
                // First visit — assign strand and enqueue.
                result.strandByRead[r1] = int8_t(s1);
                q.push(r1);
            } else if(result.strandByRead[r1] != int8_t(s1)) {
                // Already visited with a different strand — conflict.
                ++result.conflictCount;
            }
        }
    }
    return result;
}



// ============================================================================
// Lifetime helpers
// ============================================================================
void BidirectionalReadGraph::unreserve()
{
    if(edges.isOpenWithWriteAccess) {
        edges.unreserve();
    }
    if(connectivity.isOpenWithWriteAccess()) {
        connectivity.unreserve();
    }
}

void BidirectionalReadGraph::remove()
{
    if(edges.isOpen) {
        edges.remove();
    }
    if(connectivity.isOpen()) {
        connectivity.remove();
    }
}
