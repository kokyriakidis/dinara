// Dinara.
#include "Assembler.hpp"
#include "Kmer.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"

// Standard library.
#include "algorithm.hpp"
#include <thread>
#include <vector>

using namespace dinara;
using namespace std;



void Assembler::filterMarkerGraphVerticesByDistinctSubkmerCount(uint64_t threadCount)
{
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    reads->checkReadFlagsAreOpen();

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    const uint64_t k = assemblerInfo->k;

    // Same defaults as Shasta2 (--min-anchor-distinct-subkmer-count).
    const vector<uint64_t> minAnchorDistinctSubkmerCount = {4, 12, 24};

    auto shouldSkipKmerDueToLowComplexity = [&](const Kmer& kmer) -> bool {
        for(uint64_t i=0; i<minAnchorDistinctSubkmerCount.size(); i++) {
            const uint64_t subKmerLength = i + 1;
            const uint64_t minAllowedCount = minAnchorDistinctSubkmerCount[i];
            if(kmer.count(subKmerLength, k) < minAllowedCount) {
                return true;
            }
        }
        return false;
    };

    const MarkerGraph::VertexId vertexCount = markerGraph.vertexCount();
    if(vertexCount == 0) {
        return;
    }

    // Parallel scan of vertices to identify ones to be removed.
    vector<vector<MarkerGraph::VertexId>> toRemoveByThread(threadCount);
    uint64_t chunk = vertexCount / threadCount;
    if(chunk == 0) {
        chunk = 1;
    }

    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; ++t) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = (t == threadCount - 1) ? vertexCount : min<uint64_t>(vertexCount, (t + 1) * chunk);
            auto& out = toRemoveByThread[t];
            out.reserve((end - begin) / 4 + 1);

            for(MarkerGraph::VertexId vertexId=begin; vertexId<end; ++vertexId) {
                const span<const MarkerId> markerIds = markerGraph.getVertexMarkerIds(vertexId);
                if(markerIds.empty()) {
                    continue;
                }
                const MarkerId markerId = markerIds[0];
                OrientedReadId orientedReadId;
                uint32_t ordinal = 0;
                tie(orientedReadId, ordinal) = findMarkerId(markerId);
                const Kmer kmer = getOrientedReadMarkerKmer(orientedReadId, ordinal);
                if(!shouldSkipKmerDueToLowComplexity(kmer)) {
                    continue;
                }

                out.push_back(vertexId);

                // Also remove its reverse complement vertex, if present.
                const MarkerId markerIdRc = findReverseComplement(markerId);
                const MarkerGraph::CompressedVertexId vRcCompressed = markerGraph.vertexTable[markerIdRc];
                if(vRcCompressed != MarkerGraph::invalidCompressedVertexId) {
                    out.push_back(MarkerGraph::VertexId(uint64_t(vRcCompressed)));
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    vector<uint8_t> remove(vertexCount, 0);
    for(const auto& v : toRemoveByThread) {
        for(const MarkerGraph::VertexId vertexId : v) {
            if(vertexId < vertexCount) {
                remove[vertexId] = 1;
            }
        }
    }

    vector<MarkerGraph::VertexId> kept;
    kept.reserve(vertexCount);
    for(MarkerGraph::VertexId vertexId=0; vertexId<vertexCount; ++vertexId) {
        if(!remove[vertexId]) {
            kept.push_back(vertexId);
        }
    }

    const uint64_t removedCount = vertexCount - kept.size();
    cout << timestamp << "Distinct-subkmer-based marker-graph vertex filter removed "
         << removedCount << " / " << vertexCount << " vertices." << endl;
    if(removedCount == 0) {
        return;
    }

    // Keep only the selected vertices.
    MemoryMapped::Vector<MarkerGraph::VertexId> verticesToBeKept;
    verticesToBeKept.createNew(
        largeDataName("tmp-MarkerGraphVerticesToBeKeptSubkmerFilter"),
        largeDataPageSize);
    verticesToBeKept.reserveAndResize(kept.size());
    copy(kept.begin(), kept.end(), verticesToBeKept.begin());

    markerGraph.removeVertices(verticesToBeKept, largeDataPageSize, threadCount);
    verticesToBeKept.remove();
}
