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
        };

        void clear() {
            arena.clear();
            sameStrand.clear();
            reverseStrand.clear();
        }

        void reserve(size_t overlapCount, size_t arenaTokens) {
            sameStrand.reserve(overlapCount);
            reverseStrand.reserve(overlapCount);
            arena.reserve(arenaTokens);
        }

        // Append one overlap's CIGAR (raw hifiasm tokens) and its metadata.
        // Not thread-safe: call from a single thread after dedup.
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
                arena.emplace_back(CigarToken(tokens[i]));
            }
            (isSameStrand ? sameStrand : reverseStrand)[pairKey] = rec;
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

        bool empty() const { return sameStrand.empty() && reverseStrand.empty(); }

    private:
        // Flat token arena (native hifiasm frame).
        std::vector<CigarToken> arena;
        std::unordered_map<uint64_t, Record> sameStrand;
        std::unordered_map<uint64_t, Record> reverseStrand;
    };

} // namespace dinara

#endif
