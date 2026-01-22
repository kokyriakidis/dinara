#pragma once

// KmerCounter is a class that holds a hash table capable of returning
// the number of time a marker KmerId appears in
// the markers of this assembly.

// Dinara.
#include "Kmer.hpp"
#include "MappedMemoryOwner.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultithreadedObject.hpp"
#include "dinaraTypes.hpp"

#include <unordered_map>



namespace dinara {
    class KmerCounter;

    class CompressedMarker;
    class KmerDistributionInfo;
    class Reads;
    
    // Custom hasher for KmerId (which may be __uint128_t or larger).
    // std::hash does not provide a specialization for 128-bit types.
    struct KmerIdHasher {
        size_t operator()(const KmerId& k) const noexcept {
            // XOR the high and low 64-bit halves of the KmerId.
            const uint64_t* p = reinterpret_cast<const uint64_t*>(&k);
            uint64_t h = p[0];
            if constexpr (sizeof(KmerId) > 8) {
                h ^= p[1] + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return static_cast<size_t>(h);
        }
    };
}


class dinara::KmerCounter :
    public MultithreadedObject<KmerCounter>,
    public MappedMemoryOwner {
public:

    // This constructor creates the KmerIdFrequencies hash table.
    KmerCounter(
        uint64_t k,
        const Reads&,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MappedMemoryOwner& mappedMemoryOwner,
        uint64_t threadCount
        );

    // This constructor creates the KmerIdFrequencies hash table
    // from pre-calculated KmerIds (markerKmerIds).
    KmerCounter(
        uint64_t k,
        const MemoryMapped::VectorOfVectors<KmerId, uint64_t>& markerKmerIds,
        const MappedMemoryOwner& mappedMemoryOwner,
        uint64_t threadCount
        );

    // This constructor accesses an existing KmerIdFrequencies hash table.
    KmerCounter(
        uint64_t k,
        const MappedMemoryOwner& mappedMemoryOwner
        );

    bool isAvailable() const
    {
        return kmerIdFrequencies.isOpen();
    }

    uint64_t getFrequency(KmerId) const;
    uint64_t getFrequency(const Kmer&) const;
    
    // Optimized O(1) lookup for pre-canonicalized KmerIds.
    // Caller is responsible for computing the canonical form before calling.
    uint64_t getFrequencyFast(KmerId canonicalKmerId) const;
    
    // Build the frequency LUT for O(1) lookups. Call after createHistogram().
    void buildFrequencyLUT();

    // Override the frequencies stored in this KmerCounter
    // with the ones obtained from another KmerCounter.
    void overrideFrequencies(const KmerCounter&);

private:

    // Data passed in to the constructor.
    uint64_t k;
    Reads const* readsPointer = 0;
    MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t> const* markersPointer = 0;
    MemoryMapped::VectorOfVectors<KmerId, uint64_t> const* markerKmerIdsPointer = 0; // New member

    // Hashing.
    const uint32_t hashSeed = 12771;
    uint64_t hashMask;

    // A temporary hash table that stores KmerIds of each marker k-mer we encounter.
    // KmerIds are stored in canonical form. Canonical KmerId is the
    // lowest of the KmerId and the KmerId of its reverse complement.
    MemoryMapped::VectorOfVectors<KmerId, uint64_t> kmerIds;

    // This is the persistent hash table that contains in each bucket
    // pairs(KmerId, frequency).
public:
    MemoryMapped::VectorOfVectors<pair<KmerId, uint64_t>, uint64_t> kmerIdFrequencies;
    
    // Optimized LUT for O(1) frequency lookups. Built by buildFrequencyLUT().
    std::unordered_map<KmerId, uint64_t, KmerIdHasher> frequencyLUT;
private:


    // Passes 1 and 2 gather marker KmerIds in the kmerIds hash table.
    void threadFunction1(uint64_t threadId);
    void threadFunction2(uint64_t threadId);
    void threadFunction12(uint64_t pass);

    // Optimized threading functions that use markerKmerIds.
    void threadFunction1FromIds(uint64_t threadId);
    void threadFunction2FromIds(uint64_t threadId);
    void threadFunction12FromIds(uint64_t pass);

    // Passes 3 and 4 gather fill in the KmerIdFrequencies hash table.
    void threadFunction3(uint64_t threadId);
    void threadFunction4(uint64_t threadId);
    void threadFunction34(uint64_t pass);

public:

    // The frequency histogram.
    // It consists of pairs (coverage, frequency).
    void createHistogram();
    void writeHistogram(ostream&) const;
    MemoryMapped::Vector< pair<uint64_t, uint64_t> > histogram;

    void getHistogramInfo(KmerDistributionInfo&) const;
};
