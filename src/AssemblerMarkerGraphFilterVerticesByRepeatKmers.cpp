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



void Assembler::filterMarkerGraphVerticesByRepeatKmers(uint64_t threadCount)
{
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    reads->checkReadFlagsAreOpen();

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    const uint64_t k = assemblerInfo->k;

    // Same defaults as Shasta2 (--max-anchor-repeat-length).
    const vector<uint64_t> maxAnchorRepeatLength = {8, 3, 3, 3, 3};

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
                if(!shouldSkipKmerDueToRepeats(kmer)) {
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
    cout << timestamp << "Repeat-based marker-graph vertex filter removed "
         << removedCount << " / " << vertexCount << " vertices." << endl;
    if(removedCount == 0) {
        return;
    }

    // Keep only the selected vertices.
    MemoryMapped::Vector<MarkerGraph::VertexId> verticesToBeKept;
    verticesToBeKept.createNew(
        largeDataName("tmp-MarkerGraphVerticesToBeKeptRepeatFilter"),
        largeDataPageSize);
    verticesToBeKept.reserveAndResize(kept.size());
    copy(kept.begin(), kept.end(), verticesToBeKept.begin());

    markerGraph.removeVertices(verticesToBeKept, largeDataPageSize, threadCount);
    verticesToBeKept.remove();
}
