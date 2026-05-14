#pragma once

// Count-then-scatter inverted index construction (hifiasm-style).
//
// Replaces the old sort-based approach that materialized a 24-byte-per-occurrence
// intermediate array and radix-sorted it. This version:
//   1. Counts k-mer occurrences via per-thread hash tables (no positions stored).
//   2. Merges thread tables, computes prefix sums to assign each k-mer a range.
//   3. Scatters (readId, position) directly into the 8-byte CompactOccurrence array.
//
// Peak memory: ~8 bytes/occurrence (final array) + hash table overhead.
// The old approach peaked at ~48 bytes/occurrence (24-byte sort src + dst buffers).

#include "Assembler.hpp"
#include "bitReversal.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <atomic>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace dinara {

/// Compute the reverse complement of a k-mer stored in a KmerId.
/// KmerId packs a k-mer as two k-bit halves (LSB and MSB). The reverse
/// complement flips and reverses both halves, then swaps them.
inline KmerId getRcKmerId(KmerId id, uint64_t k) {
    const KmerId mask = (KmerId(1) << k) - 1;
    const KmerId lsb = id & mask;
    const KmerId msb = (id >> k) & mask;

    auto reverseBits = [&](KmerId x) -> KmerId {
        if constexpr (sizeof(KmerId) <= 8)
            return KmerId(bitReversal(uint64_t(x)) >> (64 - k));
        else
            return KmerId(bitReversal((__uint128_t)x) >> (128 - k));
    };

    KmerId rc_lsb = (~reverseBits(lsb)) & mask;
    KmerId rc_msb = (~reverseBits(msb)) & mask;
    return (rc_msb << k) | rc_lsb;
}

/// Hash function for K-mer IDs.
/// Uses Boost-style hash combine to fold 128-bit KmerId (or wider) into
/// a 64-bit hash. Used for hash-table construction in the inverted index.
inline uint64_t hashKmer(KmerId k) {
    const uint64_t* p = reinterpret_cast<const uint64_t*>(&k);
    uint64_t k1 = p[0];
    uint64_t k2 = sizeof(KmerId) > 8 ? p[1] : 0;
    return k1 ^ (k2 + 0x9e3779b9 + (k1 << 6) + (k1 >> 2));
}

namespace inverted_index_builder {

// ---------------------------------------------------------------------------
// Counting hash table used only during index construction.
// ---------------------------------------------------------------------------

struct CountEntry {
    KmerId key{};
    uint64_t start = 0;     // Prefix-sum offset into compactOccurrences.
    uint32_t count = 0;     // Number of occurrences of this k-mer.
    bool empty = true;
};

// Atomic write cursor, stored separately from CountEntry to avoid
// polluting the merge phase with atomic overhead.
struct ScatterState {
    std::atomic<uint64_t> writeCursor{0};
};

inline CountEntry* findOrInsert(
    std::vector<CountEntry>& table,
    const uint64_t mask,
    const KmerId key)
{
    uint64_t slot = hashKmer(key) & mask;
    while (true) {
        if (table[slot].empty) {
            table[slot].key = key;
            table[slot].empty = false;
            table[slot].count = 0;
            return &table[slot];
        }
        if (table[slot].key == key) {
            return &table[slot];
        }
        slot = (slot + 1) & mask;
    }
}

// ---------------------------------------------------------------------------
// Main builder function.
// ---------------------------------------------------------------------------

/// Build the inverted index using a count-then-scatter approach.
///
/// Populates:
///   - data.compactOccurrences  (8 bytes per occurrence)
///   - data.hashTable           (query-phase hash table)
///   - data.strand0CanonicalKmerIds / IsRc / Offsets
///   - data.k
inline void build(
    Assembler::AlignmentCandidatesInvertedIndexData& data,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MemoryMapped::VectorOfVectors<KmerId, uint64_t>& markerKmerIds,
    const uint64_t k,
    const uint64_t threadCount)
{
    using std::vector;
    using std::thread;
    using std::cout;
    using std::endl;

    data.k = k;
    const ReadId readCount = ReadId(markers.size() / 2);

    // =====================================================================
    // Phase 1: Compute per-read marker counts + canonical cache
    // =====================================================================

    vector<uint64_t> readMarkerCounts(size_t(readCount), 0);
    for (ReadId rId = 0; rId < readCount; ++rId) {
        const auto& rMarkers = markers[size_t(rId) << 1];
        const auto& rKmerIds = markerKmerIds[size_t(rId) << 1];
        readMarkerCounts[size_t(rId)] = std::min<size_t>(rMarkers.size(), rKmerIds.size());
    }

    data.strand0CanonicalOffsets.resize(size_t(readCount) + 1, 0);
    for (size_t r = 0; r < size_t(readCount); ++r) {
        data.strand0CanonicalOffsets[r + 1] =
            data.strand0CanonicalOffsets[r] + readMarkerCounts[r];
    }
    const uint64_t totalMarkersFound = data.strand0CanonicalOffsets.back();

    data.strand0CanonicalKmerIds.resize(totalMarkersFound);
    data.strand0CanonicalIsRc.resize(totalMarkersFound);

    cout << "Building Inverted Index for " << totalMarkersFound
         << " markers (" << readCount << " reads)." << endl;

    // Parallel canonicalization.
    {
        vector<thread> threads;
        threads.reserve(threadCount);
        for (size_t tid = 0; tid < threadCount; ++tid) {
            threads.emplace_back([&, tid]() {
                const ReadId startRead = ReadId((uint64_t(readCount) * tid) / threadCount);
                const ReadId endRead = ReadId((uint64_t(readCount) * (tid + 1)) / threadCount);

                for (ReadId rId = startRead; rId < endRead; ++rId) {
                    const auto& rKmerIds = markerKmerIds[size_t(rId) << 1];
                    const size_t n = readMarkerCounts[size_t(rId)];
                    size_t writeOff = size_t(data.strand0CanonicalOffsets[size_t(rId)]);

                    for (size_t i = 0; i < n; ++i) {
                        KmerId kId = rKmerIds[i];
                        KmerId rcKId = getRcKmerId(kId, k);
                        KmerId canonicalKId = (kId < rcKId) ? kId : rcKId;
                        uint8_t isRc = uint8_t(kId > rcKId);

                        data.strand0CanonicalKmerIds[writeOff] = canonicalKId;
                        data.strand0CanonicalIsRc[writeOff] = isRc;
                        ++writeOff;
                    }
                }
            });
        }
        for (auto& t : threads) t.join();
    }

    // =====================================================================
    // Phase 2: Count occurrences per canonical k-mer (per-thread tables)
    // =====================================================================

    struct ThreadCountTable {
        vector<CountEntry> table;
        uint64_t mask = 0;
        uint64_t distinctCount = 0;
    };
    vector<ThreadCountTable> threadTables(threadCount);

    {
        vector<thread> threads;
        threads.reserve(threadCount);
        for (size_t tid = 0; tid < threadCount; ++tid) {
            threads.emplace_back([&, tid]() {
                const ReadId startRead = ReadId((uint64_t(readCount) * tid) / threadCount);
                const ReadId endRead = ReadId((uint64_t(readCount) * (tid + 1)) / threadCount);

                uint64_t threadMarkers = 0;
                for (ReadId rId = startRead; rId < endRead; ++rId) {
                    threadMarkers += readMarkerCounts[size_t(rId)];
                }

                // Size local table at ~4x marker count for low collision rate.
                uint64_t tableSize = 1024;
                while (tableSize < threadMarkers * 4) tableSize *= 2;

                auto& tt = threadTables[tid];
                tt.table.resize(tableSize);
                tt.mask = tableSize - 1;
                tt.distinctCount = 0;

                for (ReadId rId = startRead; rId < endRead; ++rId) {
                    const size_t n = readMarkerCounts[size_t(rId)];
                    const size_t base = size_t(data.strand0CanonicalOffsets[size_t(rId)]);

                    for (size_t i = 0; i < n; ++i) {
                        const KmerId canonicalKId = data.strand0CanonicalKmerIds[base + i];
                        CountEntry* entry = findOrInsert(tt.table, tt.mask, canonicalKId);
                        if (entry->count == 0) {
                            ++tt.distinctCount;
                        }
                        ++entry->count;
                    }
                }
            });
        }
        for (auto& t : threads) t.join();
    }

    // =====================================================================
    // Phase 3: Merge per-thread tables into global counting table
    // =====================================================================

    uint64_t totalDistinctUpperBound = 0;
    for (const auto& tt : threadTables) {
        totalDistinctUpperBound += tt.distinctCount;
    }

    uint64_t globalTableSize = 1024;
    while (globalTableSize < totalDistinctUpperBound * 4) globalTableSize *= 2;

    vector<CountEntry> globalTable(globalTableSize);
    const uint64_t globalMask = globalTableSize - 1;
    uint64_t numDistinctKmers = 0;

    for (auto& tt : threadTables) {
        for (auto& entry : tt.table) {
            if (!entry.empty) {
                CountEntry* ge = findOrInsert(globalTable, globalMask, entry.key);
                if (ge->count == 0) {
                    ++numDistinctKmers;
                }
                ge->count += entry.count;
            }
        }
        // Free thread-local table immediately.
        { vector<CountEntry>().swap(tt.table); }
    }
    { vector<ThreadCountTable>().swap(threadTables); }

    cout << "Distinct K-mers found: " << numDistinctKmers << endl;

    // =====================================================================
    // Phase 4: Prefix-sum → assign each k-mer its range
    // =====================================================================

    // Build a parallel array of atomic write cursors for the scatter phase.
    vector<ScatterState> scatterStates(globalTableSize);
    {
        uint64_t offset = 0;
        for (size_t i = 0; i < globalTableSize; ++i) {
            if (!globalTable[i].empty) {
                globalTable[i].start = offset;
                scatterStates[i].writeCursor.store(offset, std::memory_order_relaxed);
                offset += globalTable[i].count;
            }
        }
        if (offset != totalMarkersFound) {
            throw std::runtime_error(
                "buildInvertedIndex count-scatter: prefix sum mismatch: "
                + std::to_string(offset) + " vs " + std::to_string(totalMarkersFound));
        }
    }

    // =====================================================================
    // Phase 5: Allocate compactOccurrences and scatter positions
    // =====================================================================

    data.compactOccurrences.resize(totalMarkersFound);

    {
        vector<thread> threads;
        threads.reserve(threadCount);
        for (size_t tid = 0; tid < threadCount; ++tid) {
            threads.emplace_back([&, tid]() {
                const ReadId startRead = ReadId((uint64_t(readCount) * tid) / threadCount);
                const ReadId endRead = ReadId((uint64_t(readCount) * (tid + 1)) / threadCount);

                for (ReadId rId = startRead; rId < endRead; ++rId) {
                    const auto& rMarkers = markers[size_t(rId) << 1];
                    const size_t n = readMarkerCounts[size_t(rId)];
                    const size_t base = size_t(data.strand0CanonicalOffsets[size_t(rId)]);

                    for (size_t i = 0; i < n; ++i) {
                        const KmerId canonicalKId = data.strand0CanonicalKmerIds[base + i];
                        const uint8_t isRc = data.strand0CanonicalIsRc[base + i];
                        const uint32_t position = rMarkers[i].position;
                        const uint32_t encodedPosition = position | (uint32_t(isRc) << 31);

                        // Look up the k-mer's slot in the global table.
                        uint64_t slot = hashKmer(canonicalKId) & globalMask;
                        while (globalTable[slot].key != canonicalKId) {
                            slot = (slot + 1) & globalMask;
                        }

                        // Atomically claim a write position.
                        const uint64_t writeIdx =
                            scatterStates[slot].writeCursor.fetch_add(1, std::memory_order_relaxed);

                        data.compactOccurrences[writeIdx] = {rId, encodedPosition};
                    }
                }
            });
        }
        for (auto& t : threads) t.join();
    }

    // Free scatter state.
    { vector<ScatterState>().swap(scatterStates); }

    // =====================================================================
    // Phase 6: Build the final query hash table
    // =====================================================================

    {
        uint64_t queryTableSize = 1;
        while (queryTableSize < numDistinctKmers * 2) queryTableSize *= 2;

        data.hashTable.resize(queryTableSize);
        const uint64_t queryMask = queryTableSize - 1;

        for (const auto& entry : globalTable) {
            if (!entry.empty) {
                uint64_t slot = hashKmer(entry.key) & queryMask;
                while (!data.hashTable[slot].empty) {
                    slot = (slot + 1) & queryMask;
                }
                data.hashTable[slot] = {entry.key, entry.start, entry.count, false};
            }
        }
    }

    // Free the global counting table.
    { vector<CountEntry>().swap(globalTable); }

    cout << "Index construction complete (count-scatter)." << endl;
}

} // namespace inverted_index_builder
} // namespace dinara
