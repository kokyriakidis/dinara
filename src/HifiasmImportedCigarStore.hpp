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
#include <cstring>
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

        ~HifiasmImportedCigarStore() { releaseCigarArena(); releaseChainArena(); }

        // The store OWNS the adopted CIGAR and chain arenas (plain C
        // allocations), so copying it would double-free. Move-only.
        HifiasmImportedCigarStore() = default;
        HifiasmImportedCigarStore(const HifiasmImportedCigarStore&) = delete;
        HifiasmImportedCigarStore& operator=(const HifiasmImportedCigarStore&) = delete;

        void clear() {
            releaseCigarArena();
            releaseChainArena();
            sameStrand.clear();
            reverseStrand.clear();
        }

        void reserve(size_t overlapCount) {
            sameStrand.reserve(overlapCount);
            reverseStrand.reserve(overlapCount);
        }

        // Take ownership of hifiasm's own CIGAR token arena, exactly as
        // adoptChainArena does for the chain anchors. hifiasm's uint16_t token
        // and dinara's CigarToken are the SAME 16 bits (op in [15:14], length
        // in [13:0]) -- CigarToken is a one-field struct over that word and its
        // uint16_t constructor takes a raw hifiasm token -- so the arena needs
        // no re-encoding to be readable as CigarToken, only the op2/op3
        // transpose, which transposeInPlace does WITHOUT copying.
        //
        // Must be a malloc/realloc allocation (out_cigar is). The store frees
        // it; the caller must not, and must not pass it to
        // hifiasm_overlaps_mem_free.
        void adoptCigarArena(uint16_t* arena, uint64_t tokenCount) {
            releaseCigarArena();
            cigarArenaData = arena;
            cigarArenaCount = tokenCount;
        }

        // Register one overlap's record WITHOUT touching a single token: the
        // tokens already sit in the adopted arena at hifiasm's own offset.
        // cigarQuerySpan/cigarTargetSpan are left zero for transposeInPlace +
        // setCigarSpans to fill.
        void addRecord(uint64_t pairKey, bool isSameStrand,
                       uint64_t cigarOffset, uint32_t cigarTokenCount,
                       uint32_t readIdQ, uint32_t readIdT,
                       uint32_t qStart, uint32_t qEnd,
                       uint32_t tStart, uint32_t tEnd) {
            Record rec;
            rec.cigarOffset = cigarOffset;
            rec.cigarTokenCount = cigarTokenCount;
            rec.readIdQ = readIdQ;
            rec.readIdT = readIdT;
            rec.qStart = qStart;
            rec.qEnd = qEnd;
            rec.tStart = tStart;
            rec.tEnd = tEnd;
            rec.isSameStrand = isSameStrand;
            (isSameStrand ? sameStrand : reverseStrand)[pairKey] = rec;
        }

        // Transpose one overlap's tokens IN PLACE in the adopted arena, from
        // hifiasm's op convention to dinara's, and report what the tokens span.
        //
        // hifiasm's op2 consumes the TARGET and op3 the QUERY; dinara's op2
        // (CigarOpIns) consumes the query and op3 (CigarOpDel) the target. So
        // the fix is to swap 2 <-> 3 and leave 0/1 alone. In the packed word the
        // op is bits [15:14], so op2 = 0b10 and op3 = 0b11 differ only in bit
        // 14, and ops 0/1 both have bit 15 clear: flipping bit 14 exactly when
        // bit 15 is set does the swap with no branch and cannot disturb the
        // 14-bit length.
        //
        // Touches only [offset, offset+count), which hifiasm assigns disjointly
        // per overlap, so callers may run this over different records
        // concurrently. It is NOT idempotent -- swapping twice restores
        // hifiasm's convention -- so call it exactly once per record.
        void transposeInPlace(uint64_t cigarOffset, uint32_t cigarTokenCount,
                              uint32_t& querySpan, uint32_t& targetSpan) const {
            uint64_t q = 0, t = 0;
            uint16_t* const first = cigarArenaData + cigarOffset;
            for(uint32_t i = 0; i < cigarTokenCount; i++) {
                uint16_t data = first[i];
                data = uint16_t(data ^ ((data >> 1) & 0x4000u));
                first[i] = data;
                const CigarToken token(data);
                const uint8_t op = token.op();
                if(opConsumesQuery(op))  q += token.len();
                if(opConsumesTarget(op)) t += token.len();
            }
            querySpan = uint32_t(q);
            targetSpan = uint32_t(t);
        }

        // Convenience for a producer that has no arena to adopt and must build
        // one overlap at a time: grow the owned arena by these raw hifiasm
        // tokens, register the record, transpose, record the spans. Exactly
        // what importAlignmentCandidatesFromMemory does across its three
        // passes, for one overlap -- and the ONLY entry point that copies a
        // token, which is why the import path does not use it. There is still
        // just one arena, and the store still owns and frees it.
        void addCopyingTokens(uint64_t pairKey, bool isSameStrand,
                              span<const uint16_t> tokens,
                              uint32_t readIdQ, uint32_t readIdT,
                              uint32_t qStart, uint32_t qEnd,
                              uint32_t tStart, uint32_t tEnd) {
            const uint64_t offset = cigarArenaCount;
            if(not tokens.empty()) {
                // realloc, so the arena stays the plain C allocation that
                // releaseCigarArena() frees -- same shape as an adopted one.
                uint16_t* const grown = static_cast<uint16_t*>(std::realloc(
                    cigarArenaData,
                    size_t(cigarArenaCount + tokens.size()) * sizeof(uint16_t)));
                DINARA_ASSERT(grown != nullptr);
                cigarArenaData = grown;
                std::memcpy(cigarArenaData + cigarArenaCount, &*tokens.begin(),
                    tokens.size() * sizeof(uint16_t));
                cigarArenaCount += tokens.size();
            }
            addRecord(pairKey, isSameStrand, offset, uint32_t(tokens.size()),
                readIdQ, readIdT, qStart, qEnd, tStart, tEnd);
            uint32_t querySpan = 0, targetSpan = 0;
            transposeInPlace(offset, uint32_t(tokens.size()), querySpan, targetSpan);
            setCigarSpans(pairKey, isSameStrand, querySpan, targetSpan);
        }

        void setCigarSpans(uint64_t pairKey, bool isSameStrand,
                           uint32_t querySpan, uint32_t targetSpan) {
            auto& m = isSameStrand ? sameStrand : reverseStrand;
            auto it = m.find(pairKey);
            if(it == m.end()) return;
            it->second.cigarQuerySpan = querySpan;
            it->second.cigarTargetSpan = targetSpan;
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

        // Token slice for a record, viewed over the adopted arena. Safe to
        // reinterpret: CigarToken is a single uint16_t with hifiasm's own bit
        // layout (static_assert(sizeof(CigarToken) == 2) in OverlapCigarStore).
        span<const CigarToken> tokensOf(const Record& rec) const {
            if(cigarArenaData == nullptr or rec.cigarTokenCount == 0) {
                return {};
            }
            return {
                reinterpret_cast<const CigarToken*>(cigarArenaData + rec.cigarOffset),
                rec.cigarTokenCount };
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
        void releaseCigarArena() {
            std::free(cigarArenaData);
            cigarArenaData = nullptr;
            cigarArenaCount = 0;
        }

        void releaseChainArena() {
            std::free(chainArenaData);
            chainArenaData = nullptr;
            chainArenaCount = 0;
        }

        // hifiasm's own flat CIGAR token arena, ADOPTED not copied -- this store
        // frees it. Record::cigarOffset/cigarTokenCount index into it using the
        // offsets hifiasm itself assigned; transposeInPlace rewrites the op bits
        // there rather than into a second arena.
        uint16_t* cigarArenaData = nullptr;
        uint64_t cigarArenaCount = 0;
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
