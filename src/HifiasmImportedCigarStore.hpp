#ifndef DINARA_HIFIASM_IMPORTED_CIGAR_STORE_HPP
#define DINARA_HIFIASM_IMPORTED_CIGAR_STORE_HPP

// Holds the NATIVE (alignment-frame) hifiasm CIGAR for each imported overlap,
// keyed by the canonical pair key and strand used by the alignment-candidate
// pipeline (see makePafEntry / OrientedReadPair: read0 = min(ReadId)).
//
// importAlignmentCandidatesFromMemory copies, for each overlap that survives
// dedup, hifiasm's packed uint16_t CIGAR tokens plus the metadata needed to
// reframe them later (the query/target ReadIds, forward-strand overlap spans,
// and strand). computeBaseAlignmentsAndStore looks the CIGAR up per candidate
// and, via normalizeHifiasmCigar, reframes + trims it to the marker interval
// instead of recomputing the base alignment with A*PA2.
//
// The tokens are stored as-is (query forward, target in alignment orientation);
// reframing into dinara's read0/read1 canonical frame happens at consumption
// time, where the read lengths are readily available.

#include "OverlapCigarStore.hpp"
#include "cstdint.hpp"
#include "vector.hpp"
#include "span.hpp"

#include <unordered_map>

namespace dinara {

    class HifiasmImportedCigarStore {
    public:
        // Per-overlap CIGAR record in hifiasm's native alignment frame.
        struct Record {
            uint64_t cigarOffset = 0;      // start in the flat token arena
            uint32_t cigarTokenCount = 0;  // token count; 0 => no CIGAR available
            // hifiasm query/target read ids (dinara ReadIds) and forward-strand
            // overlap spans, needed to reframe into the canonical read0/read1.
            uint32_t readIdQ = 0;
            uint32_t readIdT = 0;
            uint32_t qStart = 0;
            uint32_t qEnd = 0;
            uint32_t tStart = 0;
            uint32_t tEnd = 0;
            bool isSameStrand = true;
            // Native dense chain anchors for this overlap (hifiasm's colinear-DP
            // seeds), each packed (q_start<<32)|t_start in the query-forward /
            // target-alignment-orientation frame. Slice into chainArena.
            // chainCount == 0 when no native chain was imported for this pair.
            uint64_t chainOffset = 0;
            uint32_t chainCount = 0;
        };

        void clear() {
            arena.clear();
            chainArena.clear();
            sameStrand.clear();
            reverseStrand.clear();
        }

        void reserve(size_t overlapCount, size_t arenaTokens) {
            sameStrand.reserve(overlapCount);
            reverseStrand.reserve(overlapCount);
            arena.reserve(arenaTokens);
        }

        void reserveChain(size_t chainAnchors) {
            chainArena.reserve(chainAnchors);
        }

        // Append one overlap's CIGAR (raw hifiasm tokens) and its metadata.
        // Not thread-safe: call from a single thread after dedup.
        //
        // hifiasm's exported op2/op3 are the transpose of dinara's convention:
        // in hifiasm's bit_extz_t frame op2 consumes the TARGET and op3 consumes
        // the QUERY, whereas dinara's OverlapCigarStore defines op2 (CigarOpIns)
        // as query-consuming and op3 (CigarOpDel) as target-consuming. Verified
        // by base content: walking the raw tokens with op2=target/op3=query makes
        // every op0 (match) column pair identical bases, while the opposite
        // interpretation mismatches ~70% of them. Transpose op2<->op3 here, at
        // the single ingest boundary, so every downstream consumer (and the
        // recorded qStart/qEnd/tStart/tEnd spans) share dinara's convention.
        void add(uint64_t pairKey, bool isSameStrand,
                 span<const uint16_t> tokens,
                 uint32_t readIdQ, uint32_t readIdT,
                 uint32_t qStart, uint32_t qEnd,
                 uint32_t tStart, uint32_t tEnd) {
            Record rec;
            rec.cigarOffset = arena.size();
            rec.cigarTokenCount = uint32_t(tokens.size());
            rec.readIdQ = readIdQ;
            rec.readIdT = readIdT;
            rec.qStart = qStart;
            rec.qEnd = qEnd;
            rec.tStart = tStart;
            rec.tEnd = tEnd;
            rec.isSameStrand = isSameStrand;
            for(size_t i = 0; i < tokens.size(); ++i) {
                const CigarToken raw(tokens[i]);
                const uint8_t op = raw.op();
                const uint8_t dinaraOp =
                    (op == CigarOpIns) ? CigarOpDel :
                    (op == CigarOpDel) ? CigarOpIns : op;
                arena.emplace_back(CigarToken(dinaraOp, raw.len()));
            }
            (isSameStrand ? sameStrand : reverseStrand)[pairKey] = rec;
        }

        // Attach the native chain anchors to an already-added record (same
        // pairKey/strand). Copies the packed (q_start<<32)|t_start anchors into
        // the chain arena. Not thread-safe; call from the single ingest thread.
        void addChain(uint64_t pairKey, bool isSameStrand,
                      span<const uint64_t> anchors) {
            auto& m = isSameStrand ? sameStrand : reverseStrand;
            auto it = m.find(pairKey);
            if(it == m.end()) return;
            it->second.chainOffset = chainArena.size();
            it->second.chainCount = uint32_t(anchors.size());
            for(size_t i = 0; i < anchors.size(); ++i)
                chainArena.push_back(anchors[i]);
        }

        // Look up the CIGAR record for a candidate. Returns nullptr if absent.
        const Record* find(uint64_t pairKey, bool isSameStrand) const {
            const auto& m = isSameStrand ? sameStrand : reverseStrand;
            auto it = m.find(pairKey);
            return it == m.end() ? nullptr : &it->second;
        }

        // Token slice for a record.
        span<const CigarToken> tokensOf(const Record& rec) const {
            return { arena.data() + rec.cigarOffset, rec.cigarTokenCount };
        }

        // Native chain-anchor slice for a record (empty if none imported).
        span<const uint64_t> chainOf(const Record& rec) const {
            return { chainArena.data() + rec.chainOffset, rec.chainCount };
        }

        bool empty() const { return sameStrand.empty() && reverseStrand.empty(); }

    private:
        // Flat token arena (native hifiasm frame).
        std::vector<CigarToken> arena;
        // Flat native-chain anchor arena (packed (q_start<<32)|t_start).
        std::vector<uint64_t> chainArena;
        std::unordered_map<uint64_t, Record> sameStrand;
        std::unordered_map<uint64_t, Record> reverseStrand;
    };

} // namespace dinara

#endif
