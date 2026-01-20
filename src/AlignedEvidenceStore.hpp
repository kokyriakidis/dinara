#ifndef DINARA_ALIGNED_EVIDENCE_STORE_HPP
#define DINARA_ALIGNED_EVIDENCE_STORE_HPP

#include "cstdint.hpp"
#include "vector.hpp"
#include "span.hpp"
#include "DINARA_ASSERT.hpp"

namespace dinara {

    // APES / TASSD Architecture
    // Target-Aligned Separated Sparse Diffs
    
    // Stream 1: SNP Evidence (16-bit)
    // [ Delta (14 bits) | Base (2 bits) ]
    // Delta is relative to the previous SNP position on the TARGET read.
    struct SnpEvidence {
        uint16_t data;

        static const uint16_t DELTA_MASK = 0xFFFC;
        static const uint16_t BASE_MASK = 0x0003;
        static const int BASE_SHIFT = 0;
        static const int DELTA_SHIFT = 2;

        SnpEvidence() : data(0) {}
        SnpEvidence(uint16_t delta, uint8_t base) {
            DINARA_ASSERT(delta < (1<<14));
            DINARA_ASSERT(base < 4);
            data = (delta << DELTA_SHIFT) | (base & BASE_MASK);
        }

        uint16_t delta() const { return data >> DELTA_SHIFT; }
        uint8_t base() const { return data & BASE_MASK; }
    };

    // Stream 2: Indel Evidence (64-bit)
    // [ Position (32 bits) | Length (24 bits) | Type (8 bits) ]
    // Position: Absolute on Target Read. Supports 4GB (covers >1MB reads).
    // Length: Supports ~16MB Indels/SVs.
    struct IndelEvidence {
        uint64_t data;

        static const uint64_t TYPE_MASK = 0xFF;
        static const uint64_t LEN_MASK = 0xFFFFFF;
        static const uint64_t POS_MASK = 0xFFFFFFFF;
        
        static const int TYPE_SHIFT = 0;
        static const int LEN_SHIFT = 8;
        static const int POS_SHIFT = 32;

        IndelEvidence() : data(0) {}
        IndelEvidence(uint32_t targetPos, uint32_t length, uint8_t type) {
            // DINARA_ASSERT(targetPos < (1ULL<<32)); Implicit by type
            DINARA_ASSERT(length < (1<<24));    // ~16MB limit
            // DINARA_ASSERT(type < 256); Implicit by type
            
            data = ((uint64_t)targetPos << POS_SHIFT) | ((uint64_t)length << LEN_SHIFT) | ((uint64_t)type & TYPE_MASK);
        }

        uint32_t pos() const { return (uint32_t)((data >> POS_SHIFT) & POS_MASK); }
        uint32_t len() const { return (uint32_t)((data >> LEN_SHIFT) & LEN_MASK); }
        uint8_t type() const { return (uint8_t)((data >> TYPE_SHIFT) & TYPE_MASK); }
        
        bool isInsertion() const { return type() == 0; }
        bool isDeletion() const { return type() == 1; }
    };

    class AlignedEvidenceStore {
    public:
        // The Global Arenas
        // Stream 0: Projected to Read 1 (Target)
        vector<SnpEvidence> snpStream0;
        vector<IndelEvidence> indelStream0;

        // Stream 1: Projected to Read 0 (Query)
        vector<SnpEvidence> snpStream1;
        vector<IndelEvidence> indelStream1;

        // The Index
        // Maps AlignmentID -> Range in streams
        struct IndexEntry {
            uint64_t snpOffset0;
            uint64_t indelOffset0;
            uint32_t snpCount0;
            uint32_t indelCount0;
            
            uint64_t snpOffset1;
            uint64_t indelOffset1;
            uint32_t snpCount1;
            uint32_t indelCount1;
        };
        vector<IndexEntry> index;

        void clear() {
            snpStream0.clear();
            indelStream0.clear();
            snpStream1.clear();
            indelStream1.clear();
            index.clear();
        }

        // Reserve space to avoid reallocations
        void reserve(size_t indexSize, size_t snpSize, size_t indelSize) {
            index.reserve(indexSize);
            snpStream0.reserve(snpSize);
            indelStream0.reserve(indelSize);
            snpStream1.reserve(snpSize);
            indelStream1.reserve(indelSize);
        }
        
        // Start a new alignment entry
        // Returns the alignment ID
        uint32_t beginAlignment() {
            IndexEntry entry;
            entry.snpOffset0 = snpStream0.size();
            entry.indelOffset0 = indelStream0.size();
            entry.snpCount0 = 0;
            entry.indelCount0 = 0;
            
            entry.snpOffset1 = snpStream1.size();
            entry.indelOffset1 = indelStream1.size();
            entry.snpCount1 = 0;
            entry.indelCount1 = 0;
            
            index.push_back(entry);
            return (uint32_t)(index.size() - 1);
        }

        void addSnp0(uint16_t delta, uint8_t base) {
             snpStream0.emplace_back(delta, base);
             index.back().snpCount0++;
        }

        void addIndel0(uint32_t pos, uint32_t len, uint8_t type) {
            indelStream0.emplace_back(pos, len, type);
            index.back().indelCount0++;
        }

        void addSnp1(uint16_t delta, uint8_t base) {
             snpStream1.emplace_back(delta, base);
             index.back().snpCount1++;
        }

        void addIndel1(uint32_t pos, uint32_t len, uint8_t type) {
            indelStream1.emplace_back(pos, len, type);
            index.back().indelCount1++;
        }

        // Accessors
        span<const SnpEvidence> getSnps0(uint32_t alignmentId) const {
            if (alignmentId >= index.size()) return {};
            const auto& entry = index[alignmentId];
            return { snpStream0.data() + entry.snpOffset0, entry.snpCount0 };
        }

        span<const IndelEvidence> getIndels0(uint32_t alignmentId) const {
             if (alignmentId >= index.size()) return {};
             const auto& entry = index[alignmentId];
             return { indelStream0.data() + entry.indelOffset0, entry.indelCount0 };
        }

        span<const SnpEvidence> getSnps1(uint32_t alignmentId) const {
            if (alignmentId >= index.size()) return {};
            const auto& entry = index[alignmentId];
            return { snpStream1.data() + entry.snpOffset1, entry.snpCount1 };
        }

        span<const IndelEvidence> getIndels1(uint32_t alignmentId) const {
             if (alignmentId >= index.size()) return {};
             const auto& entry = index[alignmentId];
             return { indelStream1.data() + entry.indelOffset1, entry.indelCount1 };
        }

    };

}

#endif
