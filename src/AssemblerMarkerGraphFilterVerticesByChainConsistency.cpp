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
#include <thread>
#include <vector>

using namespace dinara;
using namespace std;



void Assembler::filterMarkerGraphVerticesByChainConsistency(uint64_t threadCount)
{
    cout << timestamp << "Filtering marker graph vertices by chain consistency." << endl;

    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(alignmentCandidates.candidateTable.isOpen());

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

                    const auto candidateIndices =
                        alignmentCandidates.candidateTable[readIdI.getValue()];

                    // For each candidate of readIdI, check if the other read
                    // is one of the reads j > i in this vertex.
                    for(const uint64_t candidateIndex : candidateIndices) {
                        if(foundInconsistent) break;

                        const OrientedReadPair& pair =
                            alignmentCandidates.candidates[candidateIndex];
                        const OrientedReadId otherReadId = pair.getOther(readIdI);

                        // Find otherReadId among reads j > i in this vertex.
                        for(uint64_t j = i + 1; j < n; ++j) {
                            if(readOrdinals[j].first != otherReadId) {
                                continue;
                            }

                            // Found a pair (i, j) with an alignment. Check consistency.
                            const uint32_t ordinalJ = readOrdinals[j].second;
                            const Alignment& alignment =
                                alignmentCandidatesAlignmentsData.alignments[candidateIndex];

                            if(alignment.ordinals.empty()) {
                                foundInconsistent = true;
                                break;
                            }

                            // The alignment ordinals are stored for the canonical orientation:
                            //   [0] -> OrientedReadId(readIds[0], 0)
                            //   [1] -> OrientedReadId(readIds[1], isSameStrand ? 0 : 1)
                            // If readIdI or readIdJ is on a different strand, we must
                            // convert the vertex ordinal to the canonical ordinal space.
                            const OrientedReadId canonicalRead0(pair.readIds[0], 0);
                            const OrientedReadId canonicalRead1(pair.readIds[1],
                                pair.isSameStrand ? 0 : 1);

                            // Determine which alignment slot (0 or 1) corresponds to readIdI and readIdJ.
                            int slotI, slotJ;
                            if(readIdI.getReadId() == pair.readIds[0]) {
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

                            const uint32_t chainStartI = alignment.ordinals.front()[slotI];
                            const uint32_t chainEndI   = alignment.ordinals.back()[slotI];
                            const uint32_t chainStartJ = alignment.ordinals.front()[slotJ];
                            const uint32_t chainEndJ   = alignment.ordinals.back()[slotJ];

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
