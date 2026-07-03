#pragma once

// Count-then-scatter inverted index construction (hifiasm-style).
//
// This version:
//   1. Counts k-mer occurrences via per-thread hash tables (no positions stored).
//   2. Merges thread tables into a global counting table, then builds the query
//      hash table with prefix-sum ranges (each k-mer's slice of the output).
//   3. Scatters (readId, position) directly into the 8-byte CompactOccurrence
//      array, using a write cursor parallel to the query hash table.
//
// Memory model (N = total occurrences, D = distinct k-mers, T = threads,
// s = sizeof(KmerId), typically 16):
//   - Counting/merge tables are sized by DISTINCT k-mers (load factor ~0.7),
//     grown on demand, not sized by marker count. Footprint is O(D) per thread
//     and O(D) global, using parallel key/count arrays (no 32-byte entry, no
//     unused "start" field during counting).
//   - The scatter cursor array is parallel to the query hash table, i.e. O(D),
//     replacing the earlier design's cursor array sized to the (much larger)
//     counting table.
//   - The per-occurrence canonical k-mer cache (s + 1 bytes each) is optional
//     (buildCanonicalCache). Counting and scatter recompute canonicalization
//     inline (cheap ALU, no memory), so the cache only exists to speed the
//     query phase and can be disabled to save ~17N bytes on low-RAM machines.
//
// Peak heap is therefore dominated by the 8N output array plus O(D) tables.
// Measured on synthetic data: ~25 bytes/occurrence with the canonical cache on,
// ~8 bytes/occurrence with it off, versus the previous design's ~190-245
// bytes/occurrence (per-thread tables at 4-8x marker count with 32-byte entries
// plus a counting-table-sized atomic cursor array).

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
// Compact counting hash table used only during index construction.
//
// Open-addressing, linear probing, parallel arrays:
//   keys[slot]   : the canonical k-mer for an occupied slot.
//   counts[slot] : occurrence count; 0 means the slot is empty (sentinel).
//
// Sizing is driven by the number of DISTINCT k-mers via a load-factor trigger
// (grow at ~0.7), so the footprint tracks the working set rather than the
// marker count. Storing count directly (no 32-byte struct, no unused "start"
// field during counting) minimizes bytes touched per probe.
// ---------------------------------------------------------------------------

struct CompactCountTable {
    std::vector<KmerId>   keys;
    std::vector<uint32_t> counts;   // 0 => empty
    uint64_t mask = 0;
    uint64_t size = 0;              // capacity (power of two)
    uint64_t occupied = 0;         // number of distinct keys stored
    uint64_t growThreshold = 0;    // occupied count that triggers a grow

    void init(uint64_t initialCapacityPow2) {
        size = initialCapacityPow2;
        mask = size - 1;
        keys.assign(size, KmerId(0));
        counts.assign(size, 0u);
        occupied = 0;
        growThreshold = (size * 7) / 10;   // load factor 0.7
    }

    void grow() {
        std::vector<KmerId>   oldKeys;
        std::vector<uint32_t> oldCounts;
        oldKeys.swap(keys);
        oldCounts.swap(counts);
        const uint64_t oldSize = size;

        size <<= 1;
        mask = size - 1;
        keys.assign(size, KmerId(0));
        counts.assign(size, 0u);
        growThreshold = (size * 7) / 10;

        for (uint64_t i = 0; i < oldSize; ++i) {
            if (oldCounts[i] != 0) {
                uint64_t slot = hashKmer(oldKeys[i]) & mask;
                while (counts[slot] != 0) slot = (slot + 1) & mask;
                keys[slot] = oldKeys[i];
                counts[slot] = oldCounts[i];
            }
        }
    }

    // Increment the count for key, inserting if absent. Grows on demand.
    inline void add(const KmerId key) {
        if (occupied >= growThreshold) grow();
        uint64_t slot = hashKmer(key) & mask;
        while (true) {
            if (counts[slot] == 0) {
                keys[slot] = key;
                counts[slot] = 1;
                ++occupied;
                return;
            }
            if (keys[slot] == key) { ++counts[slot]; return; }
            slot = (slot + 1) & mask;
        }
    }

    // Merge another key with an explicit count (used when merging thread tables).
    inline void addCount(const KmerId key, const uint32_t c) {
        if (occupied >= growThreshold) grow();
        uint64_t slot = hashKmer(key) & mask;
        while (true) {
            if (counts[slot] == 0) {
                keys[slot] = key;
                counts[slot] = c;
                ++occupied;
                return;
            }
            if (keys[slot] == key) { counts[slot] += c; return; }
            slot = (slot + 1) & mask;
        }
    }
};

// ---------------------------------------------------------------------------
// Main builder function.
// ---------------------------------------------------------------------------

/// Build the inverted index using a count-then-scatter approach.
///
/// Populates:
///   - data.compactOccurrences  (8 bytes per occurrence)
///   - data.hashTable           (query-phase hash table)
///   - data.k
/// and, when buildCanonicalCache is true, the optional per-marker canonical
/// cache (data.strand0CanonicalKmerIds / IsRc / Offsets) that lets the query
/// phase skip recomputing reverse complements. When false, those vectors are
/// left empty and the query phase recomputes canonicalization on the fly.
///
/// Internal counting/scatter passes recompute canonicalization inline from
/// markerKmerIds, so construction never depends on the canonical cache.
inline void build(
    Assembler::AlignmentCandidatesInvertedIndexData& data,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MemoryMapped::VectorOfVectors<KmerId, uint64_t>& markerKmerIds,
    const uint64_t k,
    const uint64_t threadCount,
    const bool buildCanonicalCache = true)
{
    using std::vector;
    using std::thread;
    using std::cout;
    using std::endl;

    data.k = k;
    const ReadId readCount = ReadId(markers.size() / 2);

    // =====================================================================
    // Phase 1: Per-read marker counts (and optional canonical cache)
    // =====================================================================

    vector<uint64_t> readMarkerCounts(size_t(readCount), 0);
    for (ReadId rId = 0; rId < readCount; ++rId) {
        const auto& rMarkers = markers[size_t(rId) << 1];
        const auto& rKmerIds = markerKmerIds[size_t(rId) << 1];
        readMarkerCounts[size_t(rId)] = std::min<size_t>(rMarkers.size(), rKmerIds.size());
    }

    // Per-read base offsets into the flat per-marker arrays. Always needed for
    // the optional cache; also serves as the strand0CanonicalOffsets output.
    vector<uint64_t> readBase(size_t(readCount) + 1, 0);
    for (size_t r = 0; r < size_t(readCount); ++r) {
        readBase[r + 1] = readBase[r] + readMarkerCounts[r];
    }
    const uint64_t totalMarkersFound = readBase.back();

    cout << "Building Inverted Index for " << totalMarkersFound
         << " markers (" << readCount << " reads)." << endl;

    if (buildCanonicalCache) {
        data.strand0CanonicalOffsets = readBase;   // copy: also used by query phase
        data.strand0CanonicalKmerIds.resize(totalMarkersFound);
        data.strand0CanonicalIsRc.resize(totalMarkersFound);

        vector<thread> threads;
        threads.reserve(threadCount);
        for (size_t tid = 0; tid < threadCount; ++tid) {
            threads.emplace_back([&, tid]() {
                const ReadId startRead = ReadId((uint64_t(readCount) * tid) / threadCount);
                const ReadId endRead = ReadId((uint64_t(readCount) * (tid + 1)) / threadCount);
                for (ReadId rId = startRead; rId < endRead; ++rId) {
                    const auto& rKmerIds = markerKmerIds[size_t(rId) << 1];
                    const size_t n = readMarkerCounts[size_t(rId)];
                    size_t writeOff = size_t(readBase[size_t(rId)]);
                    for (size_t i = 0; i < n; ++i) {
                        const KmerId kId = rKmerIds[i];
                        const KmerId rcKId = getRcKmerId(kId, k);
                        data.strand0CanonicalKmerIds[writeOff] = (kId < rcKId) ? kId : rcKId;
                        data.strand0CanonicalIsRc[writeOff] = uint8_t(kId > rcKId);
                        ++writeOff;
                    }
                }
            });
        }
        for (auto& t : threads) t.join();
    }

    // Helper: canonical k-mer for marker i of a read (cache or recompute).
    // Generic over the row type so it works with the memory-mapped span the
    // real VectorOfVectors returns.
    const bool useCache = buildCanonicalCache;
    auto canonicalAt = [&](size_t i, const auto& rKmerIds, uint64_t base) -> KmerId {
        if (useCache) {
            return data.strand0CanonicalKmerIds[base + i];
        }
        const KmerId kId = rKmerIds[i];
        const KmerId rcKId = getRcKmerId(kId, k);
        return (kId < rcKId) ? kId : rcKId;
    };

    // Reasonable initial per-thread table capacity: enough to hold a small
    // fraction of this thread's markers before the first grow, capped so we
    // never pre-allocate anywhere near marker-count territory.
    auto initialCapacityFor = [](uint64_t threadMarkers) -> uint64_t {
        uint64_t cap = 1024;
        const uint64_t target = std::max<uint64_t>(1024, threadMarkers / 8);
        while (cap < target) cap <<= 1;
        return cap;
    };

    // =====================================================================
    // Phase 2: Count occurrences per canonical k-mer (per-thread tables)
    // =====================================================================

    vector<CompactCountTable> threadTables(threadCount);
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

                auto& tt = threadTables[tid];
                tt.init(initialCapacityFor(threadMarkers));

                for (ReadId rId = startRead; rId < endRead; ++rId) {
                    const auto& rKmerIds = markerKmerIds[size_t(rId) << 1];
                    const size_t n = readMarkerCounts[size_t(rId)];
                    const uint64_t base = readBase[size_t(rId)];
                    for (size_t i = 0; i < n; ++i) {
                        tt.add(canonicalAt(i, rKmerIds, base));
                    }
                }
            });
        }
        for (auto& t : threads) t.join();
    }

    // =====================================================================
    // Phase 3: Merge per-thread tables into a global counting table
    // =====================================================================

    uint64_t distinctUpperBound = 0;
    for (const auto& tt : threadTables) distinctUpperBound += tt.occupied;

    CompactCountTable globalTable;
    {
        uint64_t cap = 1024;
        const uint64_t target = std::max<uint64_t>(1024, (distinctUpperBound * 10) / 7 + 1);
        while (cap < target) cap <<= 1;
        globalTable.init(cap);
    }

    for (auto& tt : threadTables) {
        for (uint64_t s = 0; s < tt.size; ++s) {
            if (tt.counts[s] != 0) {
                globalTable.addCount(tt.keys[s], tt.counts[s]);
            }
        }
        // Free thread-local table immediately.
        { vector<KmerId>().swap(tt.keys); }
        { vector<uint32_t>().swap(tt.counts); }
    }
    { vector<CompactCountTable>().swap(threadTables); }

    const uint64_t numDistinctKmers = globalTable.occupied;
    cout << "Distinct K-mers found: " << numDistinctKmers << endl;

    // =====================================================================
    // Phase 4: Build query hash table with prefix-sum ranges; init cursors
    // =====================================================================
    //
    // The query hash table (indexed by hash of the canonical k-mer) doubles as
    // the scatter lookup structure. A cursor array parallel to it (indexed by
    // query slot) holds each k-mer's next free write position. This keeps the
    // cursor footprint O(distinct), not O(counting-table-size).

    uint64_t queryTableSize = 1;
    while (queryTableSize < numDistinctKmers * 2) queryTableSize <<= 1;
    const uint64_t queryMask = queryTableSize - 1;

    data.hashTable.assign(queryTableSize, {});
    vector<std::atomic<uint64_t>> cursors(queryTableSize);

    {
        uint64_t offset = 0;
        for (uint64_t s = 0; s < globalTable.size; ++s) {
            if (globalTable.counts[s] == 0) continue;
            const KmerId key = globalTable.keys[s];
            const uint32_t count = globalTable.counts[s];

            uint64_t slot = hashKmer(key) & queryMask;
            while (!data.hashTable[slot].empty) slot = (slot + 1) & queryMask;
            data.hashTable[slot] = {key, offset, count, false};
            cursors[slot].store(offset, std::memory_order_relaxed);
            offset += count;
        }
        if (offset != totalMarkersFound) {
            throw std::runtime_error(
                "buildInvertedIndex count-scatter: prefix sum mismatch: "
                + std::to_string(offset) + " vs " + std::to_string(totalMarkersFound));
        }
    }

    // Free the global counting table before allocating the output array.
    { vector<KmerId>().swap(globalTable.keys); }
    { vector<uint32_t>().swap(globalTable.counts); }

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
                    const auto& rKmerIds = markerKmerIds[size_t(rId) << 1];
                    const size_t n = readMarkerCounts[size_t(rId)];
                    const uint64_t base = readBase[size_t(rId)];

                    for (size_t i = 0; i < n; ++i) {
                        KmerId canonicalKId;
                        uint8_t isRc;
                        if (useCache) {
                            canonicalKId = data.strand0CanonicalKmerIds[base + i];
                            isRc = data.strand0CanonicalIsRc[base + i];
                        } else {
                            const KmerId kId = rKmerIds[i];
                            const KmerId rcKId = getRcKmerId(kId, k);
                            canonicalKId = (kId < rcKId) ? kId : rcKId;
                            isRc = uint8_t(kId > rcKId);
                        }
                        const uint32_t position = rMarkers[i].position;
                        const uint32_t encodedPosition = position | (uint32_t(isRc) << 31);

                        uint64_t slot = hashKmer(canonicalKId) & queryMask;
                        while (data.hashTable[slot].key != canonicalKId) {
                            slot = (slot + 1) & queryMask;
                        }
                        const uint64_t writeIdx =
                            cursors[slot].fetch_add(1, std::memory_order_relaxed);
                        data.compactOccurrences[writeIdx] = {rId, encodedPosition};
                    }
                }
            });
        }
        for (auto& t : threads) t.join();
    }

    // Free scatter cursors.
    { vector<std::atomic<uint64_t>>().swap(cursors); }

    cout << "Index construction complete (count-scatter)." << endl;
}

} // namespace inverted_index_builder
} // namespace dinara
