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

#include <cstdlib>
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
            // What the CIGAR itself accounts for: the summed query-consuming and
            // target-consuming op lengths. These are NOT always qEnd-qStart /
            // tEnd-tStart. hifiasm's in-memory record pairs the overlap BOX
            // (x_pos_s..x_pos_e) with the concatenation of only the windows that
            // aligned -- ecovlp.cpp skips is_ualn_win when building the token
            // stream -- so a window that failed to align leaves the box and the
            // tokens describing different spans, with nothing marking where the
            // hole was. (hifiasm's file PAF cannot show this: it emits one record
            // per window, each carrying its own window's bounds, so coordinates
            // and CIGAR always agree there.) Measured: 36/29851 records (0.12%)
            // on the small GIAB fixture, 1201/251791 (0.48%, mean ~1068bp) at
            // 251k overlaps -- it scales the wrong way.
            //
            // So before walking the tokens, compare these against the box. When
            // they disagree the token stream cannot be walked from qStart/tStart
            // without desynchronising, and the record carries no information
            // about whether the missing span is leading, trailing or internal.
            uint32_t cigarQuerySpan = 0;
            uint32_t cigarTargetSpan = 0;
            // Native dense chain anchors for this overlap (hifiasm's colinear-DP
            // seeds), each packed (q_start<<32)|t_start in the query-forward /
            // target-alignment-orientation frame. Slice into chainArena.
            // chainCount == 0 when no native chain was imported for this pair.
            uint64_t chainOffset = 0;
            uint32_t chainCount = 0;
        };

        ~HifiasmImportedCigarStore() { releaseChainArena(); }

        // The store OWNS the adopted chain arena (a plain C allocation), so
        // copying it would double-free. Move-only.
        HifiasmImportedCigarStore() = default;
        HifiasmImportedCigarStore(const HifiasmImportedCigarStore&) = delete;
        HifiasmImportedCigarStore& operator=(const HifiasmImportedCigarStore&) = delete;

        void clear() {
            arena.clear();
            releaseChainArena();
            sameStrand.clear();
            reverseStrand.clear();
        }

        void reserve(size_t overlapCount, size_t arenaTokens) {
            sameStrand.reserve(overlapCount);
            reverseStrand.reserve(overlapCount);
            arena.reserve(arenaTokens);
        }

        // Take ownership of hifiasm's own native-chain arena instead of copying
        // it. hifiasm hands back one flat array of packed anchors plus, per
        // overlap, the offset/length of that overlap's slice within it -- which
        // is exactly the layout this store needs. Copying it bought nothing but
        // a lifetime: the caller used to free the arena right after import,
        // while chainOf() is not consumed until computeBaseAlignmentsAndStore.
        // Adopting it removes the copy AND the duplicate residency.
        //
        // Measured on E821 chr12:11-17Mb (161.6M anchors, 1.20 GiB): the copy
        // cost 715 ms, of which only 173 ms was the memcpy itself (7.5 GB/s) --
        // the other 542 ms was 315,663 minor page faults at ~1.7 us each,
        // zeroing pages that were about to be overwritten. Pre-faulting the
        // reservation confirmed it (715 ms -> 173 ms with the cost merely moved
        // into the reserve). Not copying at all removes both halves.
        //
        // `arena` MUST be a malloc/realloc allocation (it is: hifiasm's sink
        // grows it with realloc and the caller receives it via out_chain). The
        // store frees it in clear() and in its destructor; the caller must not.
        void adoptChainArena(uint64_t* arena, uint64_t anchorCount) {
            releaseChainArena();
            chainArenaData = arena;
            chainArenaCount = anchorCount;
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
            // The transpose already visits every token, so accumulate what the
            // CIGAR spans here rather than in a second pass. Lengths are summed
            // in DINARA's convention (post-transpose), matching qStart/qEnd and
            // tStart/tEnd.
            uint64_t querySpan = 0;
            uint64_t targetSpan = 0;
            for(size_t i = 0; i < tokens.size(); ++i) {
                const CigarToken raw(tokens[i]);
                const uint8_t op = raw.op();
                const uint8_t dinaraOp =
                    (op == CigarOpIns) ? uint8_t(CigarOpDel) :
                    (op == CigarOpDel) ? uint8_t(CigarOpIns) : op;
                arena.emplace_back(CigarToken(dinaraOp, raw.len()));
                if(opConsumesQuery(dinaraOp))  querySpan  += raw.len();
                if(opConsumesTarget(dinaraOp)) targetSpan += raw.len();
            }
            rec.cigarQuerySpan  = uint32_t(querySpan);
            rec.cigarTargetSpan = uint32_t(targetSpan);
            (isSameStrand ? sameStrand : reverseStrand)[pairKey] = rec;
        }

        // Point an already-added record at its slice of the adopted chain arena
        // (see adoptChainArena). Records the offset/count hifiasm already
        // computed -- no anchor is copied or touched. Out-of-range slices are
        // ignored rather than trusted, so a record can never hand out a span
        // past the arena. Not thread-safe; call from the single ingest thread.
        void setChain(uint64_t pairKey, bool isSameStrand,
                      uint64_t anchorOffset, uint32_t anchorCount) {
            auto& m = isSameStrand ? sameStrand : reverseStrand;
            auto it = m.find(pairKey);
            if(it == m.end()) return;
            if(chainArenaData == nullptr) return;
            if(anchorOffset + anchorCount > chainArenaCount) return;
            it->second.chainOffset = anchorOffset;
            it->second.chainCount = anchorCount;
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
            if(chainArenaData == nullptr or rec.chainCount == 0) {
                return {};
            }
            return { chainArenaData + rec.chainOffset, rec.chainCount };
        }

        bool empty() const { return sameStrand.empty() && reverseStrand.empty(); }

        // Visit every stored record. The store is keyed by (pairKey, strand)
        // for per-pair lookup, so a consumer that needs "all overlaps involving
        // read r" has to build that index itself; this is the way in.
        template<class F> void forEachRecord(F&& f) const {
            for(const auto& [key, rec]: sameStrand)    { (void)key; f(rec); }
            for(const auto& [key, rec]: reverseStrand) { (void)key; f(rec); }
        }

        uint64_t recordCount() const {
            return sameStrand.size() + reverseStrand.size();
        }

    private:
        void releaseChainArena() {
            std::free(chainArenaData);
            chainArenaData = nullptr;
            chainArenaCount = 0;
        }

        // Flat token arena (native hifiasm frame).
        std::vector<CigarToken> arena;
        // hifiasm's own flat native-chain anchor arena (packed
        // (q_start<<32)|t_start), ADOPTED not copied -- this store frees it.
        // Record::chainOffset/chainCount index into it directly, using the
        // offsets hifiasm itself assigned.
        uint64_t* chainArenaData = nullptr;
        uint64_t chainArenaCount = 0;
        std::unordered_map<uint64_t, Record> sameStrand;
        std::unordered_map<uint64_t, Record> reverseStrand;
    };

} // namespace dinara

#endif
