#ifndef DINARA_ALIGNED_EVIDENCE_STORE_HPP
#define DINARA_ALIGNED_EVIDENCE_STORE_HPP

#include "cstdint.hpp"
#include "vector.hpp"
#include "span.hpp"
#include "DINARA_ASSERT.hpp"

namespace dinara {

    // Optional: enable per-alignment monotonicity checks for appended evidence.
    // This is useful when changing evidence generation code, but can be expensive
    // in hot loops. Enable with -DDINARA_ALIGNED_EVIDENCE_STORE_ORDER_CHECK=1.
    #ifndef DINARA_ALIGNED_EVIDENCE_STORE_ORDER_CHECK
    #define DINARA_ALIGNED_EVIDENCE_STORE_ORDER_CHECK 0
    #endif

    // APES / TASSD Architecture
    // Target-Aligned Separated Sparse Diffs
    
    // Stream 1: SNP Evidence (16-bit, delta-encoded)
    // [ Delta (14 bits) | Base (2 bits) ]
    // Delta is relative to the current position while decoding.
    // A token with delta()==MAX_DELTA is a "hop" (no SNP) used to bridge large gaps.
    // To represent an exact delta==MAX_DELTA mismatch, we emit a hop token followed
    // by a SNP token with delta==0.
    struct SnpEvidence {
        uint16_t data;

        static constexpr uint16_t BASE_MASK = 0x3;
        static constexpr int BASE_SHIFT = 0;
        static constexpr int DELTA_SHIFT = 2;
        static constexpr uint16_t DELTA_MASK = 0x3FFF; // 14 bits.

        static constexpr uint16_t MAX_DELTA = DELTA_MASK; // Reserved "hop".

        // Kept for legacy bounds checks in alignment code.
        static constexpr uint32_t POS_MASK = 0x3FFFFFFF; // 30 bits.

        SnpEvidence() : data(0) {}
        SnpEvidence(uint16_t delta, uint8_t base) {
            DINARA_ASSERT(delta <= DELTA_MASK);
            DINARA_ASSERT(base < 4);
            data = (uint16_t(delta << DELTA_SHIFT) | (base & BASE_MASK));
        }

        uint16_t delta() const { return uint16_t(data >> DELTA_SHIFT); }
        uint8_t base() const { return uint8_t(data & BASE_MASK); }
        bool isHop() const { return delta() == MAX_DELTA; }
    };

    struct SnpCheckpoint {
        // Token index within the per-alignment SNP token slice.
        uint32_t tokenIndex;
        // Absolute position at the start of tokenIndex (before applying tokenIndex delta).
        uint32_t pos;
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
        vector<SnpCheckpoint> snpCheckpoints0;
        vector<IndelEvidence> indelStream0;

        // Stream 1: Projected to Read 0 (Query)
        vector<SnpEvidence> snpStream1;
        vector<SnpCheckpoint> snpCheckpoints1;
        vector<IndelEvidence> indelStream1;

        // The Index
        // Maps AlignmentID -> Range in streams
        struct IndexEntry {
            uint64_t snpOffset0;
            uint64_t indelOffset0;
            uint32_t snpCount0;
            uint32_t indelCount0;
            uint64_t snpCheckpointOffset0;
            uint32_t snpCheckpointCount0;
            
            uint64_t snpOffset1;
            uint64_t indelOffset1;
            uint32_t snpCount1;
            uint32_t indelCount1;
            uint64_t snpCheckpointOffset1;
            uint32_t snpCheckpointCount1;
        };
        vector<IndexEntry> index;

        // Checkpoint stride in tokens (includes hops). Smaller -> faster seeks, more memory.
        static constexpr uint32_t SnpCheckpointStride = 64;

        void clear() {
            snpStream0.clear();
            snpCheckpoints0.clear();
            indelStream0.clear();
            snpStream1.clear();
            snpCheckpoints1.clear();
            indelStream1.clear();
            index.clear();
            currentSnpPos0 = 0;
            currentSnpPos1 = 0;
        }

        // Reserve space to avoid reallocations.
        void reserve(
            size_t indexSize,
            size_t snpSize0,
            size_t snpCheckpointSize0,
            size_t indelSize0,
            size_t snpSize1,
            size_t snpCheckpointSize1,
            size_t indelSize1
        ) {
            index.reserve(indexSize);
            snpStream0.reserve(snpSize0);
            snpCheckpoints0.reserve(snpCheckpointSize0);
            indelStream0.reserve(indelSize0);
            snpStream1.reserve(snpSize1);
            snpCheckpoints1.reserve(snpCheckpointSize1);
            indelStream1.reserve(indelSize1);
        }
        void reserve(size_t indexSize, size_t snpSize, size_t indelSize) {
            const size_t ckptSize = (snpSize / SnpCheckpointStride) + 2;
            reserve(indexSize, snpSize, ckptSize, indelSize, snpSize, ckptSize, indelSize);
        }
        
        // Start a new alignment entry
        // Returns the alignment ID
        uint32_t beginAlignment() {
            IndexEntry entry;
            entry.snpOffset0 = snpStream0.size();
            entry.indelOffset0 = indelStream0.size();
            entry.snpCount0 = 0;
            entry.indelCount0 = 0;
            entry.snpCheckpointOffset0 = snpCheckpoints0.size();
            entry.snpCheckpointCount0 = 0;
            
            entry.snpOffset1 = snpStream1.size();
            entry.indelOffset1 = indelStream1.size();
            entry.snpCount1 = 0;
            entry.indelCount1 = 0;
            entry.snpCheckpointOffset1 = snpCheckpoints1.size();
            entry.snpCheckpointCount1 = 0;
            
            index.push_back(entry);
            currentSnpPos0 = 0;
            currentSnpPos1 = 0;
            return (uint32_t)(index.size() - 1);
        }

        void addSnp0(uint32_t position, uint8_t base) {
            DINARA_ASSERT(position <= SnpEvidence::POS_MASK);
            DINARA_ASSERT(base < 4);
            auto& entry = index.back();
#if DINARA_ALIGNED_EVIDENCE_STORE_ORDER_CHECK
            DINARA_ASSERT(position >= currentSnpPos0);
#endif
            uint32_t remaining = position - currentSnpPos0;
            while (remaining >= SnpEvidence::MAX_DELTA) {
                maybeAddSnpCheckpoint0(entry);
                snpStream0.emplace_back(SnpEvidence::MAX_DELTA, 0);
                entry.snpCount0++;
                currentSnpPos0 += SnpEvidence::MAX_DELTA;
                remaining -= SnpEvidence::MAX_DELTA;
            }
            maybeAddSnpCheckpoint0(entry);
            snpStream0.emplace_back(uint16_t(remaining), base);
            entry.snpCount0++;
            currentSnpPos0 += remaining;
        }

        void addIndel0(uint32_t pos, uint32_t len, uint8_t type) {
            auto& entry = index.back();
#if DINARA_ALIGNED_EVIDENCE_STORE_ORDER_CHECK
            if (entry.indelCount0 > 0) {
                const uint32_t lastPos = indelStream0[entry.indelOffset0 + entry.indelCount0 - 1].pos();
                DINARA_ASSERT(pos >= lastPos);
            }
#endif
            indelStream0.emplace_back(pos, len, type);
            entry.indelCount0++;
        }

        void addSnp1(uint32_t position, uint8_t base) {
            DINARA_ASSERT(position <= SnpEvidence::POS_MASK);
            DINARA_ASSERT(base < 4);
            auto& entry = index.back();
#if DINARA_ALIGNED_EVIDENCE_STORE_ORDER_CHECK
            DINARA_ASSERT(position >= currentSnpPos1);
#endif
            uint32_t remaining = position - currentSnpPos1;
            while (remaining >= SnpEvidence::MAX_DELTA) {
                maybeAddSnpCheckpoint1(entry);
                snpStream1.emplace_back(SnpEvidence::MAX_DELTA, 0);
                entry.snpCount1++;
                currentSnpPos1 += SnpEvidence::MAX_DELTA;
                remaining -= SnpEvidence::MAX_DELTA;
            }
            maybeAddSnpCheckpoint1(entry);
            snpStream1.emplace_back(uint16_t(remaining), base);
            entry.snpCount1++;
            currentSnpPos1 += remaining;
        }

        void addIndel1(uint32_t pos, uint32_t len, uint8_t type) {
            auto& entry = index.back();
#if DINARA_ALIGNED_EVIDENCE_STORE_ORDER_CHECK
            if (entry.indelCount1 > 0) {
                const uint32_t lastPos = indelStream1[entry.indelOffset1 + entry.indelCount1 - 1].pos();
                DINARA_ASSERT(pos >= lastPos);
            }
#endif
            indelStream1.emplace_back(pos, len, type);
            entry.indelCount1++;
        }

        // Iterate SNP evidence (decoded absolute positions) in [beginPos, endPos).
        // The callback is invoked with (pos, base) for true SNP tokens; hop tokens
        // are not reported.
        // Raw accessors return the delta-encoded token slice (including hops).
        span<const SnpEvidence> getSnps0(uint32_t alignmentId) const {
            if (alignmentId >= index.size()) return {};
            const auto& entry = index[alignmentId];
            return { snpStream0.data() + entry.snpOffset0, entry.snpCount0 };
        }
        template<class F>
        void forEachSnp0InRange(uint32_t alignmentId, uint32_t beginPos, uint32_t endPos, F&& f) const {
            if (alignmentId >= index.size() || beginPos >= endPos) return;
            const auto& entry = index[alignmentId];
            if (entry.snpCount0 == 0) return;
            const SnpEvidence* tokens = snpStream0.data() + entry.snpOffset0;
            const SnpCheckpoint* checkpoints = snpCheckpoints0.data() + entry.snpCheckpointOffset0;
            const uint32_t checkpointCount = entry.snpCheckpointCount0;

            uint32_t tokenIndex = 0;
            uint32_t pos = 0;
            if (checkpointCount) {
                auto it = std::upper_bound(
                    checkpoints, checkpoints + checkpointCount, beginPos,
                    [](uint32_t value, const SnpCheckpoint& cp) { return value < cp.pos; }
                );
                if (it != checkpoints) --it;
                tokenIndex = it->tokenIndex;
                pos = it->pos;
            }

            for (uint32_t i = tokenIndex; i < entry.snpCount0; ++i) {
                const SnpEvidence ev = tokens[i];
                pos += ev.delta();
                if (ev.isHop()) continue;
                if (pos < beginPos) continue;
                if (pos >= endPos) break;
                f(pos, ev.base());
            }
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
        template<class F>
        void forEachSnp1InRange(uint32_t alignmentId, uint32_t beginPos, uint32_t endPos, F&& f) const {
            if (alignmentId >= index.size() || beginPos >= endPos) return;
            const auto& entry = index[alignmentId];
            if (entry.snpCount1 == 0) return;
            const SnpEvidence* tokens = snpStream1.data() + entry.snpOffset1;
            const SnpCheckpoint* checkpoints = snpCheckpoints1.data() + entry.snpCheckpointOffset1;
            const uint32_t checkpointCount = entry.snpCheckpointCount1;

            uint32_t tokenIndex = 0;
            uint32_t pos = 0;
            if (checkpointCount) {
                auto it = std::upper_bound(
                    checkpoints, checkpoints + checkpointCount, beginPos,
                    [](uint32_t value, const SnpCheckpoint& cp) { return value < cp.pos; }
                );
                if (it != checkpoints) --it;
                tokenIndex = it->tokenIndex;
                pos = it->pos;
            }

            for (uint32_t i = tokenIndex; i < entry.snpCount1; ++i) {
                const SnpEvidence ev = tokens[i];
                pos += ev.delta();
                if (ev.isHop()) continue;
                if (pos < beginPos) continue;
                if (pos >= endPos) break;
                f(pos, ev.base());
            }
        }

        span<const IndelEvidence> getIndels1(uint32_t alignmentId) const {
             if (alignmentId >= index.size()) return {};
             const auto& entry = index[alignmentId];
             return { indelStream1.data() + entry.indelOffset1, entry.indelCount1 };
        }

    private:
        uint32_t currentSnpPos0 = 0;
        uint32_t currentSnpPos1 = 0;

        void maybeAddSnpCheckpoint0(IndexEntry& entry) {
            if ((entry.snpCount0 % SnpCheckpointStride) == 0) {
                snpCheckpoints0.push_back({entry.snpCount0, currentSnpPos0});
                entry.snpCheckpointCount0++;
            }
        }
        void maybeAddSnpCheckpoint1(IndexEntry& entry) {
            if ((entry.snpCount1 % SnpCheckpointStride) == 0) {
                snpCheckpoints1.push_back({entry.snpCount1, currentSnpPos1});
                entry.snpCheckpointCount1++;
            }
        }

    };

}

#endif
