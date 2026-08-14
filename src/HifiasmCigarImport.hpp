#ifndef DINARA_HIFIASM_CIGAR_IMPORT_HPP
#define DINARA_HIFIASM_CIGAR_IMPORT_HPP

// Normalize a hifiasm per-overlap CIGAR into dinara's canonical read0/read1
// frame.
//
// This operates on tokens already in dinara's OverlapCigarStore convention
// (op 0='=', 1='X', 2='I' consumes read0/query, 3='D' consumes read1/target).
// hifiasm's raw export uses the TRANSPOSE of that for the indel ops, but the
// swap is applied once at ingest in HifiasmImportedCigarStore::add(), so by the
// time a CIGAR reaches this function op I already consumes the query and op D
// the target.
//
// In that convention a hifiasm overlap's CIGAR is in the ALIGNMENT frame:
//   - the query read (q_id) runs forward,
//   - the target read (t_id) runs in the alignment orientation (forward for a
//     same-strand overlap, reverse-complemented for a reverse overlap).
//
// dinara canonicalizes every read pair to read0 = min(ReadId), read1 =
// max(ReadId) (see makePafEntry). read0 always runs
// forward; read1 runs in the alignment orientation (reverse-complemented for a
// reverse overlap), exactly as consumers walk it (e.g. read1Start = rlen - te
// in AssemblerWindowProjectedAlignmentLeafSnarls.cpp).
//
// So the hifiasm CIGAR needs no op remapping, but it may need reframing when
// dinara's read0 is hifiasm's TARGET (i.e. t_id has the smaller ReadId):
//
//   readIdQ < readIdT  (no swap): read0 = query, read1 = target.
//       Tokens unchanged.
//   readIdQ > readIdT  (swap), same strand: read0 = target, read1 = query.
//       Transpose the alignment (swap I<->D; column order preserved).
//   readIdQ > readIdT  (swap), reverse strand: read0 = target, read1 = query.
//       Reverse-complement the whole alignment (reverse column order) AND
//       transpose (swap I<->D). Match/mismatch ops are unaffected by either
//       transform.
//
// The output anchors (read0Start/End, read1Start/End) are half-open spans:
// read0 in forward coordinates of dinara's read0; read1 in the alignment
// orientation coordinates of dinara's read1 (for a reverse overlap that is the
// reverse-complemented frame, i.e. readLen - end).

#include "OverlapCigarStore.hpp"
#include "cstdint.hpp"
#include "vector.hpp"
#include "span.hpp"

namespace dinara {

    struct NormalizedHifiasmCigar {
        std::vector<CigarToken> tokens;   // canonical read0/read1 frame
        uint32_t read0Start = 0;          // forward coords on dinara read0
        uint32_t read0End   = 0;
        uint32_t read1Start = 0;          // alignment-orientation coords on read1
        uint32_t read1End   = 0;
    };

    // Merge a raw op onto the end of a token vector, coalescing with the
    // previous token when the op matches (keeps runs contiguous; the store's
    // MAX_LEN splitting still applies on final push).
    inline void appendCoalesced(std::vector<CigarToken>& out, uint8_t op, uint32_t len)
    {
        if(len == 0) return;
        while(len > CigarToken::MAX_LEN) {
            out.emplace_back(op, CigarToken::MAX_LEN);
            len -= CigarToken::MAX_LEN;
        }
        out.emplace_back(op, uint16_t(len));
    }

    // Transpose a single op: an op that consumed only the query becomes one that
    // consumes only the target, and vice versa. Match/mismatch are symmetric.
    inline uint8_t transposeOp(uint8_t op)
    {
        if(op == CigarOpIns) return CigarOpDel;
        if(op == CigarOpDel) return CigarOpIns;
        return op; // match / mismatch
    }

    // Normalize the hifiasm CIGAR for one overlap into dinara's canonical frame.
    //
    //   rawTokens         : hifiasm CIGAR tokens (alignment frame, see above)
    //   readIdQ, readIdT  : dinara ReadIds of hifiasm's query and target reads
    //   qStart,qEnd       : query forward-strand overlap span (half-open)
    //   tStart,tEnd       : target forward-strand overlap span (half-open)
    //   qLen, tLen        : full base lengths of the query and target reads
    //   isSameStrand      : true for a same-strand overlap, false for reverse
    inline NormalizedHifiasmCigar normalizeHifiasmCigar(
        span<const CigarToken> rawTokens,
        uint64_t readIdQ, uint64_t readIdT,
        uint32_t qStart, uint32_t qEnd,
        uint32_t tStart, uint32_t tEnd,
        uint32_t qLen, uint32_t tLen,
        bool isSameStrand)
    {
        NormalizedHifiasmCigar r;
        const bool swap = (readIdQ > readIdT); // dinara read0 = hifiasm target

        if(!swap) {
            // read0 = query (forward), read1 = target (alignment orientation).
            r.tokens.assign(rawTokens.begin(), rawTokens.end());
            r.read0Start = qStart;
            r.read0End   = qEnd;
            if(isSameStrand) {
                r.read1Start = tStart;
                r.read1End   = tEnd;
            } else {
                r.read1Start = tLen - tEnd;   // == bridge cigar_t_start
                r.read1End   = tLen - tStart;
            }
            return r;
        }

        // swap: read0 = target (forward), read1 = query (alignment orientation).
        r.read0Start = tStart;
        r.read0End   = tEnd;
        if(isSameStrand) {
            r.read1Start = qStart;
            r.read1End   = qEnd;
        } else {
            r.read1Start = qLen - qEnd;
            r.read1End   = qLen - qStart;
        }

        if(isSameStrand) {
            // Transpose only: swap I<->D, preserve column order.
            r.tokens.reserve(rawTokens.size());
            for(size_t i = 0; i < rawTokens.size(); ++i) {
                appendCoalesced(r.tokens, transposeOp(rawTokens[i].op()), rawTokens[i].len());
            }
        } else {
            // Reverse-complement the whole alignment (reverse column order) and
            // transpose (swap I<->D).
            r.tokens.reserve(rawTokens.size());
            for(size_t i = rawTokens.size(); i-- > 0; ) {
                appendCoalesced(r.tokens, transposeOp(rawTokens[i].op()), rawTokens[i].len());
            }
        }
        return r;
    }

} // namespace dinara

#endif
