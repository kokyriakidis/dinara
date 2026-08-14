// Filter marker graph vertices where reads were grouped by transitive collapse
// at k-mer positions outside their chaining range.
//
// For each vertex, check all unique pairs of reads. If any pair has an
// alignment where the vertex ordinals fall outside the chaining range,
// the entire vertex is removed.

#include "Assembler.hpp"
#include "findMarkerId.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <fstream>
#include <thread>
#include <vector>

using namespace dinara;
using namespace std;



void Assembler::filterMarkerGraphVerticesByChainConsistency(uint64_t threadCount)
{
    cout << timestamp << "Filtering marker graph vertices by chain consistency." << endl;

    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    checkAlignmentDataAreOpen();
    DINARA_ASSERT(alignmentTable.isOpen());

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const MarkerGraph::VertexId vertexCount = markerGraph.vertexCount();
    if(vertexCount == 0) {
        return;
    }

    // Parallel scan of vertices.
    vector<vector<MarkerGraph::VertexId>> toRemoveByThread(threadCount);
    const uint64_t chunk = max<uint64_t>(1, (vertexCount + threadCount - 1) / threadCount);

    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t = 0; t < threadCount; ++t) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = min<uint64_t>(vertexCount, (t + 1) * chunk);
            auto& out = toRemoveByThread[t];

            // Thread-local reusable buffers.
            vector<pair<OrientedReadId, uint32_t>> readOrdinals;

            for(MarkerGraph::VertexId vertexId = begin; vertexId < end; ++vertexId) {
                const span<const MarkerId> markerIds = markerGraph.getVertexMarkerIds(vertexId);
                const uint64_t n = markerIds.size();
                if(n < 2) {
                    continue;
                }

                // Get (orientedReadId, ordinal) for each marker in this vertex.
                readOrdinals.resize(n);
                for(uint64_t i = 0; i < n; ++i) {
                    readOrdinals[i] = dinara::findMarkerId(markerIds[i], *markers);
                }

                // Check all unique pairs (i, j) with i < j.
                // For each i, scan its candidate table once and check
                // all j's that appear in it.
                bool foundInconsistent = false;

                for(uint64_t i = 0; i < n && !foundInconsistent; ++i) {
                    const OrientedReadId readIdI = readOrdinals[i].first;
                    const uint32_t ordinalI = readOrdinals[i].second;

                    // Chain endpoints per read pair come from the STORED alignment
                    // (alignmentData / alignmentTable), not from the precomputed
                    // chaining store: in direct-chaining mode the latter is empty,
                    // and AlignmentInfo::data[slot].firstOrdinal/lastOrdinal capture
                    // the same chain start/end ordinals in canonical orientation in
                    // both modes.
                    const auto alignmentIndices = alignmentTable[readIdI.getValue()];

                    // For each stored alignment of readIdI, check if the other
                    // read is one of the reads j > i in this vertex.
                    for(const uint32_t alignmentIndex : alignmentIndices) {
                        if(foundInconsistent) break;

                        const AlignmentData& ad = alignmentData[alignmentIndex];
                        const OrientedReadId canonicalRead0(ad.readIds[0], 0);
                        const OrientedReadId canonicalRead1(ad.readIds[1],
                            ad.isSameStrand ? 0 : 1);
                        // The partner of readIdI in this alignment (canonical frame).
                        const OrientedReadId otherReadId =
                            (readIdI.getReadId() == ad.readIds[0]) ?
                                canonicalRead1 : canonicalRead0;

                        // Find otherReadId among reads j > i in this vertex.
                        for(uint64_t j = i + 1; j < n; ++j) {
                            if(readOrdinals[j].first != otherReadId) {
                                continue;
                            }

                            // Found a pair (i, j) with an alignment. Check consistency.
                            const uint32_t ordinalJ = readOrdinals[j].second;

                            if(ad.info.markerCount == 0) {
                                foundInconsistent = true;
                                break;
                            }

                            // Determine which alignment slot (0 or 1) corresponds to readIdI and readIdJ.
                            int slotI, slotJ;
                            if(readIdI.getReadId() == ad.readIds[0]) {
                                slotI = 0;
                                slotJ = 1;
                            } else {
                                slotI = 1;
                                slotJ = 0;
                            }

                            // Convert vertex ordinals to canonical ordinal space if needed.
                            const OrientedReadId canonI = (slotI == 0) ? canonicalRead0 : canonicalRead1;
                            const OrientedReadId canonJ = (slotJ == 0) ? canonicalRead0 : canonicalRead1;

                            uint32_t canonOrdinalI = ordinalI;
                            if(readIdI != canonI) {
                                // readIdI is on the opposite strand from canonical.
                                const uint32_t mc = uint32_t(markers->size(readIdI.getValue()));
                                canonOrdinalI = mc - 1 - ordinalI;
                            }

                            uint32_t canonOrdinalJ = ordinalJ;
                            if(readOrdinals[j].first != canonJ) {
                                const uint32_t mc = uint32_t(markers->size(readOrdinals[j].first.getValue()));
                                canonOrdinalJ = mc - 1 - ordinalJ;
                            }

                            const uint32_t chainStartI = ad.info.data[slotI].firstOrdinal;
                            const uint32_t chainEndI   = ad.info.data[slotI].lastOrdinal;
                            const uint32_t chainStartJ = ad.info.data[slotJ].firstOrdinal;
                            const uint32_t chainEndJ   = ad.info.data[slotJ].lastOrdinal;

                            if(canonOrdinalI < chainStartI || canonOrdinalI > chainEndI ||
                               canonOrdinalJ < chainStartJ || canonOrdinalJ > chainEndJ) {
                                foundInconsistent = true;
                            }

                            break; // Only one j can match this otherReadId.
                        }
                    }
                }

                if(foundInconsistent) {
                    out.push_back(vertexId);

                    // Also remove its reverse complement vertex, if present.
                    const MarkerId markerIdRc = findReverseComplement(markerIds[0]);
                    const MarkerGraph::CompressedVertexId vRcCompressed = markerGraph.vertexTable[markerIdRc];
                    if(vRcCompressed != MarkerGraph::invalidCompressedVertexId) {
                        out.push_back(MarkerGraph::VertexId(uint64_t(vRcCompressed)));
                    }
                }
            }
        });
    }

    for(auto& th : threads) {
        th.join();
    }

    // Build removal mask.
    vector<uint8_t> remove(vertexCount, 0);
    for(const auto& v : toRemoveByThread) {
        for(const MarkerGraph::VertexId vertexId : v) {
            if(vertexId < vertexCount) {
                remove[vertexId] = 1;
            }
        }
    }

    // Collect kept vertices.
    vector<MarkerGraph::VertexId> kept;
    kept.reserve(vertexCount);
    for(MarkerGraph::VertexId vertexId = 0; vertexId < vertexCount; ++vertexId) {
        if(!remove[vertexId]) {
            kept.push_back(vertexId);
        }
    }

    const uint64_t removedCount = vertexCount - kept.size();
    cout << timestamp << "Chain-consistency filter: removed "
         << removedCount << " / " << vertexCount << " marker graph vertices." << endl;

    // Diagnostic (env DINARA_CHAIN_CONSISTENCY_DUMP=1): dump the flagged vertex
    // ids, one per line, for comparison against other candidate filters that
    // target the same class of transitive-collapse errors.
    if(getenv("DINARA_CHAIN_CONSISTENCY_DUMP") != nullptr) {
        ofstream out("ChainConsistencyFlaggedVertices.csv");
        for(MarkerGraph::VertexId vertexId = 0; vertexId < vertexCount; ++vertexId) {
            if(remove[vertexId]) {
                out << vertexId << "\n";
            }
        }
    }

    if(removedCount == 0) {
        return;
    }

    MemoryMapped::Vector<MarkerGraph::VertexId> verticesToBeKept;
    verticesToBeKept.createNew(
        largeDataName("tmp-MarkerGraphVerticesToBeKeptChainConsistencyFilter"),
        largeDataPageSize);
    verticesToBeKept.reserveAndResize(kept.size());
    copy(kept.begin(), kept.end(), verticesToBeKept.begin());

    markerGraph.removeVertices(verticesToBeKept, largeDataPageSize, threadCount);
    verticesToBeKept.remove();
}
