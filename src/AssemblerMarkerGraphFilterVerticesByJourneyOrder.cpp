// Remove marker graph vertices created by transitive-collapse false merges,
// detected via a much cheaper signal than
// filterMarkerGraphVerticesByChainConsistency's per-pair alignment-ordinal-range
// check.
//
// Idea: walk each (strand 0) oriented read's own markers in ordinal order and
// map each to its marker graph vertex, skipping markers not in any vertex. This
// gives, for that read, the sequence of vertices it visits IN READ ORDER -- for
// free, no alignment lookups needed, since ordinal order on a single read is
// trivially position order on that read. For every CONSECUTIVE pair (A, B) in
// that sequence, the read is asserting "A comes before B". If some OTHER read
// asserts the opposite ("B comes before A") for the same pair, that is a direct,
// cheap witness that A and/or B cannot both be single, correctly-placed
// genomic loci -- one of them is a transitive-collapse false merge. This is the
// same class of error filterMarkerGraphVerticesByChainConsistency targets, but
// detected via cross-read order agreement on consecutive vertices instead of
// reconstructing alignment chain ranges. A windowed sweep (checking gaps up to
// 50 journey positions apart, not just immediate neighbors) showed recall
// against the chain-consistency filter's flagged set grows roughly linearly
// with window size and stays far below it even at window 50, confirming
// consecutive-only is the intended scope (not a stepping stone to something
// wider): this is a narrower, cheaper signal, not an attempt to replicate that
// filter's coverage.
//
// This filter is complementary to, not a subset of, filterMarkerGraphVerticesByChainConsistency:
// that filter only checks read pairs that are direct alignment CANDIDATES; this
// one doesn't care whether the two witnessing reads ever aligned to each other,
// only that each independently visited the same two vertices. On one test
// region every flagged vertex here was also independently flagged by that
// filter; on another, ~30% were NOT (the chain-consistency filter missed them
// because the witnessing read pair never became an alignment candidate) --
// i.e. this filter can catch real errors that one misses, not just a cheaper
// echo of it. Runs before that filter so it sees a slightly smaller,
// cheaply-pre-cleaned vertex set.
//
// Implementation: emit each read's directed (from, to) adjacent-vertex pairs in
// parallel (embarrassingly parallel -- reads are independent), merge into one
// sorted, deduplicated list, then scan it in parallel for pairs whose mirror
// (to, from) is also present. Flagged vertices (plus their reverse-complement
// mirror, since the marker graph is strand-symmetric by construction) are
// removed.

#include "Assembler.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <fstream>
#include <thread>
#include <vector>

using namespace dinara;
using namespace std;



void Assembler::filterMarkerGraphVerticesByJourneyOrderConsistency(uint64_t threadCount)
{
    cout << timestamp << "Journey-order-consistency check begins." << endl;
    const auto t0 = std::chrono::steady_clock::now();
    auto secondsSince = [](const std::chrono::steady_clock::time_point& t) {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
    };

    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t vertexCount = markerGraph.vertexCount();
    if(vertexCount == 0) {
        cout << timestamp << "Journey-order-consistency check: no vertices, skipping." << endl;
        return;
    }
    // Packing (from, to) into one uint64_t (32 bits each) requires vertex ids
    // to fit in 32 bits. Comfortably true for any real genome; guard it anyway.
    DINARA_ASSERT(vertexCount < (uint64_t(1) << 32));

    const uint64_t readCount = reads->readCount();

    // Phase 1 (parallel over reads): emit each read's adjacent-vertex pairs.
    const auto tPhase1 = std::chrono::steady_clock::now();
    vector<vector<uint64_t>> pairsByThread(threadCount);
    {
        const uint64_t chunk = max<uint64_t>(1, (readCount + threadCount - 1) / threadCount);
        vector<thread> threads;
        threads.reserve(threadCount);
        for(uint64_t t = 0; t < threadCount; ++t) {
            threads.emplace_back([&, t]() {
                const uint64_t begin = t * chunk;
                const uint64_t end = min<uint64_t>(readCount, (t + 1) * chunk);
                auto& out = pairsByThread[t];

                vector<uint64_t> seq;
                for(uint64_t readId = begin; readId < end; ++readId) {
                    const OrientedReadId oid(ReadId(readId), 0);
                    const uint64_t markerCount = markers->size(oid.getValue());

                    seq.clear();
                    for(uint32_t ordinal = 0; ordinal < markerCount; ++ordinal) {
                        const MarkerId markerId = getMarkerId(oid, ordinal);
                        const MarkerGraph::CompressedVertexId vCompressed = markerGraph.vertexTable[markerId];
                        if(vCompressed == MarkerGraph::invalidCompressedVertexId) {
                            continue;
                        }
                        seq.push_back(uint64_t(vCompressed));
                    }

                    // Only strictly consecutive vertices in this read's own
                    // sequence: a windowed sweep (checking gaps up to 50) showed
                    // recall against filterMarkerGraphVerticesByChainConsistency's
                    // flagged set grows roughly linearly with window size and
                    // stays far below it even at window 50 -- that filter is
                    // catching a different class of thing, not just the same
                    // violations further apart, so widening the window here
                    // doesn't buy anything and only adds cost.
                    for(uint64_t i = 1; i < seq.size(); ++i) {
                        const uint64_t from = seq[i - 1];
                        const uint64_t to = seq[i];
                        if(from == to) {
                            continue;
                        }
                        out.push_back((from << 32) | to);
                    }
                }
            });
        }
        for(auto& th : threads) {
            th.join();
        }
    }
    const double phase1Seconds = secondsSince(tPhase1);

    // Phase 2 (serial): merge, sort, deduplicate.
    const auto tPhase2 = std::chrono::steady_clock::now();
    uint64_t totalPairs = 0;
    for(const auto& v : pairsByThread) {
        totalPairs += v.size();
    }
    vector<uint64_t> allPairs;
    allPairs.reserve(totalPairs);
    for(auto& v : pairsByThread) {
        allPairs.insert(allPairs.end(), v.begin(), v.end());
        vector<uint64_t>().swap(v);
    }
    std::sort(allPairs.begin(), allPairs.end());
    allPairs.erase(std::unique(allPairs.begin(), allPairs.end()), allPairs.end());
    const double phase2Seconds = secondsSince(tPhase2);

    // Phase 3 (parallel over the sorted unique pairs): flag any pair whose
    // mirror is also present.
    const auto tPhase3 = std::chrono::steady_clock::now();
    const uint64_t n = allPairs.size();
    vector<vector<uint64_t>> conflictsByThread(threadCount);
    {
        const uint64_t chunk = max<uint64_t>(1, (n + threadCount - 1) / threadCount);
        vector<thread> threads;
        threads.reserve(threadCount);
        for(uint64_t t = 0; t < threadCount; ++t) {
            threads.emplace_back([&, t]() {
                const uint64_t begin = t * chunk;
                const uint64_t end = min<uint64_t>(n, (t + 1) * chunk);
                auto& out = conflictsByThread[t];
                for(uint64_t i = begin; i < end; ++i) {
                    const uint64_t key = allPairs[i];
                    const uint64_t from = key >> 32;
                    const uint64_t to = key & 0xFFFFFFFFULL;
                    const uint64_t mirror = (to << 32) | from;
                    if(std::binary_search(allPairs.begin(), allPairs.end(), mirror)) {
                        out.push_back(from);
                        out.push_back(to);
                    }
                }
            });
        }
        for(auto& th : threads) {
            th.join();
        }
    }
    const double phase3Seconds = secondsSince(tPhase3);

    // Collect flagged vertices, then add their reverse-complement mirrors
    // (the marker graph is strand-symmetric by construction, matching how
    // filterMarkerGraphVerticesByChainConsistency also removes the RC vertex).
    const auto tPhase4 = std::chrono::steady_clock::now();
    vector<uint8_t> flagged(vertexCount, 0);
    vector<uint64_t> flaggedList;
    for(const auto& v : conflictsByThread) {
        for(const uint64_t vid : v) {
            if(!flagged[vid]) {
                flagged[vid] = 1;
                flaggedList.push_back(vid);
            }
        }
    }
    const uint64_t directCount = flaggedList.size();
    for(const uint64_t vid : flaggedList) {
        const MarkerId markerId0 = markerGraph.getVertexMarkerIds(MarkerGraph::VertexId(vid))[0];
        const MarkerId markerIdRc = findReverseComplement(markerId0);
        const MarkerGraph::CompressedVertexId vRc = markerGraph.vertexTable[markerIdRc];
        if(vRc != MarkerGraph::invalidCompressedVertexId) {
            flagged[uint64_t(vRc)] = 1;
        }
    }
    const uint64_t flaggedCount = std::count(flagged.begin(), flagged.end(), uint8_t(1));
    const double phase4Seconds = secondsSince(tPhase4);

    // Dump for cross-checking against filterMarkerGraphVerticesByChainConsistency's
    // own flagged-vertex dump (DINARA_CHAIN_CONSISTENCY_DUMP=1): same file
    // format, one vertex id per line.
    {
        ofstream out("JourneyOrderFlaggedVertices.csv");
        for(uint64_t vertexId = 0; vertexId < vertexCount; ++vertexId) {
            if(flagged[vertexId]) {
                out << vertexId << "\n";
            }
        }
    }

    cout << timestamp << "Journey-order-consistency check found"
         << " directPairs=" << totalPairs
         << " uniquePairs=" << n
         << " directConflictVertices=" << directCount
         << " conflictVerticesWithRcMirror=" << flaggedCount
         << " / " << vertexCount << " total vertices." << endl;
    cout << timestamp << "Journey-order-consistency timing:"
         << " phase1(emit)=" << phase1Seconds << "s"
         << " phase2(sort+dedup)=" << phase2Seconds << "s"
         << " phase3(conflict scan)=" << phase3Seconds << "s"
         << " phase4(rc mirror)=" << phase4Seconds << "s" << endl;

    // Diagnostic escape hatch (env DINARA_JOURNEY_ORDER_DRYRUN=1): report and
    // dump the flagged set as above, but skip the actual removal below. Used
    // to compare this filter's flagged set against
    // filterMarkerGraphVerticesByChainConsistency's on the SAME unmodified
    // graph (that filter runs immediately afterward and would otherwise see a
    // graph this filter already pruned).
    if(getenv("DINARA_JOURNEY_ORDER_DRYRUN") != nullptr) {
        cout << timestamp << "Journey-order-consistency check: dry run, not removing." << endl;
        return;
    }

    // Remove the flagged vertices, same mechanism filterMarkerGraphVerticesByChainConsistency uses.
    const auto tPhase5 = std::chrono::steady_clock::now();
    if(flaggedCount == 0) {
        cout << timestamp << "Journey-order-consistency check: nothing to remove." << endl;
        return;
    }

    vector<MarkerGraph::VertexId> kept;
    kept.reserve(vertexCount - flaggedCount);
    for(uint64_t vertexId = 0; vertexId < vertexCount; ++vertexId) {
        if(!flagged[vertexId]) {
            kept.push_back(MarkerGraph::VertexId(vertexId));
        }
    }

    MemoryMapped::Vector<MarkerGraph::VertexId> verticesToBeKept;
    verticesToBeKept.createNew(
        largeDataName("tmp-MarkerGraphVerticesToBeKeptJourneyOrderFilter"),
        largeDataPageSize);
    verticesToBeKept.reserveAndResize(kept.size());
    copy(kept.begin(), kept.end(), verticesToBeKept.begin());

    markerGraph.removeVertices(verticesToBeKept, largeDataPageSize, threadCount);
    verticesToBeKept.remove();
    const double phase5Seconds = secondsSince(tPhase5);

    const double totalSeconds = secondsSince(t0);
    cout << timestamp << "Journey-order-consistency filter: removed "
         << flaggedCount << " / " << vertexCount << " marker graph vertices"
         << " (phase5(remove)=" << phase5Seconds << "s"
         << " total=" << totalSeconds << "s)." << endl;
}
