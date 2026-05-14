#ifndef DINARA_OVERLAP_CIGAR_STORE_HPP
#define DINARA_OVERLAP_CIGAR_STORE_HPP

// Per-overlap CIGAR storage using hifiasm-style packed uint16_t tokens.
//
// Each token is a uint16_t encoding (op << 14 | length):
//   op 0 = match      (both sequences identical at these positions)
//   op 1 = mismatch   (substitution: both consume 1 base each per unit)
//   op 2 = insertion   (bases in read1/target not in read0/query)
//   op 3 = deletion    (bases in read0/query not in read1/target)
//
// Length field is 14 bits (max 16383). Runs longer than that are split
// across consecutive tokens with the same op; pop_token() coalesces them.
//
// All tokens for all alignments live in a single flat arena.
// An index maps alignmentId -> (offset, count) into the arena.

#include "cstdint.hpp"
#include "vector.hpp"
#include "span.hpp"
#include "DINARA_ASSERT.hpp"

namespace dinara {

    // Packed CIGAR token: 2-bit op + 14-bit length.
    struct CigarToken {
        uint16_t data;

        static constexpr int OP_SHIFT = 14;
        static constexpr uint16_t LEN_MASK = 0x3FFF; // 14 bits
        static constexpr uint32_t MAX_LEN = LEN_MASK;

        CigarToken() : data(0) {}
        CigarToken(uint16_t raw) : data(raw) {}
        CigarToken(uint8_t op, uint16_t len) {
            DINARA_ASSERT(op < 4);
            DINARA_ASSERT(len <= MAX_LEN);
            data = (uint16_t(op) << OP_SHIFT) | len;
        }

        uint8_t op() const { return uint8_t(data >> OP_SHIFT); }
        uint16_t len() const { return data & LEN_MASK; }
    };
    static_assert(sizeof(CigarToken) == 2);

    class OverlapCigarStore {
    public:
        // The flat token arena shared across all alignments.
        vector<CigarToken> arena;

        // Per-alignment index into the arena.
        struct IndexEntry {
            uint64_t offset; // Start position in arena.
            uint32_t count;  // Number of tokens.
        };
        vector<IndexEntry> index;

        void clear() {
            arena.clear();
            index.clear();
        }

        void reserve(size_t indexSize, size_t arenaSize) {
            index.reserve(indexSize);
            arena.reserve(arenaSize);
        }

        // Begin a new alignment's CIGAR. Returns the alignment's cigar ID
        // (index into this store's index, NOT the global alignmentId).
        uint32_t beginAlignment() {
            IndexEntry entry;
            entry.offset = arena.size();
            entry.count = 0;
            index.push_back(entry);
            return uint32_t(index.size() - 1);
        }

        // Push a single CIGAR operation. Runs exceeding MAX_LEN are
        // automatically split into multiple tokens.
        void pushOp(uint8_t op, uint32_t length) {
            DINARA_ASSERT(op < 4);
            DINARA_ASSERT(!index.empty());
            auto& entry = index.back();
            while (length > CigarToken::MAX_LEN) {
                arena.emplace_back(op, CigarToken::MAX_LEN);
                entry.count++;
                length -= CigarToken::MAX_LEN;
            }
            if (length > 0) {
                arena.emplace_back(op, uint16_t(length));
                entry.count++;
            }
        }

        // Convenience: push a match run.
        void pushMatch(uint32_t length) { pushOp(0, length); }

        // Convenience: push a mismatch run.
        void pushMismatch(uint32_t length) { pushOp(1, length); }

        // Convenience: push an insertion (extra bases in read1/target).
        void pushInsertion(uint32_t length) { pushOp(2, length); }

        // Convenience: push a deletion (extra bases in read0/query).
        void pushDeletion(uint32_t length) { pushOp(3, length); }

        // Retrieve the token slice for a given cigar ID.
        span<const CigarToken> getTokens(uint32_t cigarId) const {
            if (cigarId >= index.size()) return {};
            const auto& entry = index[cigarId];
            return { arena.data() + entry.offset, entry.count };
        }

        // Resumable cursor for walking a CIGAR across successive query
        // ranges. Saves (tokenIndex, xk, yk) so that the next
        // walkRangeWithCursor call can resume from the saved position
        // instead of scanning from the start. This matches hifiasm's
        // ovlp_cur_xoff/yoff/coff pattern used in the sliding-window
        // phasing loop.
        //
        // Usage:
        //   OverlapCigarStore::Cursor cursor;
        //   cursor.reset(cigarId, read0Start, read1Start, store);
        //   for each window [s, e):
        //       store.walkRangeWithCursor(cursor, s, e, callback);
        //
        // The cursor is valid as long as the underlying store is not modified.
        struct Cursor {
            uint32_t cigarId = uint32_t(-1);
            uint32_t tokenIndex = 0;  // Current position in the token slice.
            uint64_t xk = 0;         // Read0 (query) position at tokenIndex.
            uint64_t yk = 0;         // Read1 (target) position at tokenIndex.

            // Token slice bounds (cached from the store's index).
            const CigarToken* tokens = nullptr;
            uint32_t tokenCount = 0;

            void reset(uint32_t id, uint64_t read0Start, uint64_t read1Start,
                       const OverlapCigarStore& store) {
                cigarId = id;
                auto span = store.getTokens(id);
                tokens = span.data();
                tokenCount = uint32_t(span.size());
                tokenIndex = 0;
                xk = read0Start;
                yk = read1Start;
            }

            bool valid() const { return cigarId != uint32_t(-1) && tokens != nullptr; }
        };

        // Iterate the CIGAR for a given cigar ID, coalescing consecutive
        // tokens with the same op (as hifiasm's pop_trace does).
        // Callback signature: void(uint8_t op, uint32_t length)
        template<class F>
        void forEachOp(uint32_t cigarId, F&& f) const {
            const auto tokens = getTokens(cigarId);
            uint32_t i = 0;
            while (i < tokens.size()) {
                const uint8_t op = tokens[i].op();
                uint32_t totalLen = tokens[i].len();
                i++;
                // Coalesce consecutive tokens with the same op.
                while (i < tokens.size() && tokens[i].op() == op) {
                    totalLen += tokens[i].len();
                    i++;
                }
                f(op, totalLen);
            }
        }

        // Walk the CIGAR with full position tracking on both reads.
        // Like hifiasm's extract_sub_cigar_hc, this maintains xk (read0/query
        // position) and yk (read1/target position) as it walks each op.
        //
        // Callback signature:
        //   void(uint8_t op, uint32_t len, uint64_t read0Pos, uint64_t read1Pos)
        //
        // read0Pos/read1Pos are the starting positions of the op on each read.
        // For match/mismatch: both advance by len.
        // For insertion (op 2): only read1Pos advances by len.
        // For deletion  (op 3): only read0Pos advances by len.
        //
        // read0Start/read1Start are the base positions of the first aligned
        // base on each read (i.e. the midpoint of the first aligned marker).
        template<class F>
        void forEachOpWithPositions(
            uint32_t cigarId,
            uint64_t read0Start,
            uint64_t read1Start,
            F&& f) const
        {
            const auto tokens = getTokens(cigarId);
            uint64_t xk = read0Start;
            uint64_t yk = read1Start;
            uint32_t i = 0;
            while (i < tokens.size()) {
                const uint8_t op = tokens[i].op();
                uint32_t totalLen = tokens[i].len();
                i++;
                while (i < tokens.size() && tokens[i].op() == op) {
                    totalLen += tokens[i].len();
                    i++;
                }
                f(op, totalLen, xk, yk);
                switch(op) {
                    case 0: case 1: xk += totalLen; yk += totalLen; break;
                    case 2: yk += totalLen; break;
                    case 3: xk += totalLen; break;
                }
            }
        }

        // Walk only the sub-range [queryStart, queryEnd) on read0, calling
        // the callback for each op that overlaps that range. Ops are clipped
        // to the range boundaries.
        //
        // This is the equivalent of hifiasm's extract_sub_cigar_hc seeking
        // to position s and walking until position e.
        //
        // Callback signature:
        //   void(uint8_t op, uint32_t len, uint64_t read0Pos, uint64_t read1Pos)
        //
        // For ops partially overlapping the range, len is the clipped length
        // and read0Pos/read1Pos are adjusted to the clipped start.
        template<class F>
        void walkRange(
            uint32_t cigarId,
            uint64_t read0Start,
            uint64_t read1Start,
            uint64_t queryStart,
            uint64_t queryEnd,
            F&& f) const
        {
            const auto tokens = getTokens(cigarId);
            uint64_t xk = read0Start;
            uint64_t yk = read1Start;
            uint32_t i = 0;
            while (i < tokens.size()) {
                const uint8_t op = tokens[i].op();
                uint32_t totalLen = tokens[i].len();
                i++;
                while (i < tokens.size() && tokens[i].op() == op) {
                    totalLen += tokens[i].len();
                    i++;
                }

                const uint64_t xkEnd = xk + ((op != 2) ? totalLen : 0);
                const uint64_t ykEnd = yk + ((op != 3) ? totalLen : 0);

                if(op == 2) {
                    // Insertion: doesn't consume read0. Report it if we're
                    // within the query range (insertion is "between" two
                    // read0 positions, anchored at xk).
                    if(xk >= queryStart && xk < queryEnd) {
                        f(op, totalLen, xk, yk);
                    }
                } else {
                    // Match, mismatch, or deletion: consumes read0.
                    // Clip to [queryStart, queryEnd).
                    if(xkEnd > queryStart && xk < queryEnd) {
                        const uint64_t clipStart = (xk < queryStart) ? queryStart : xk;
                        const uint64_t clipEnd = (xkEnd > queryEnd) ? queryEnd : xkEnd;
                        const uint32_t clipLen = uint32_t(clipEnd - clipStart);
                        const uint64_t skipBases = clipStart - xk;

                        uint64_t adjYk = yk;
                        if(op == 0 || op == 1) {
                            adjYk += skipBases; // match/mismatch: read1 advances in lockstep
                        }
                        // deletion: read1 doesn't advance, so adjYk stays at yk

                        f(op, clipLen, clipStart, adjYk);
                    }
                }

                // Past the range — stop early.
                if(xkEnd >= queryEnd && op != 2) break;

                xk = xkEnd;
                yk = ykEnd;
            }
        }

        // Walk a query range [queryStart, queryEnd) using a resumable cursor.
        //
        // Like walkRange, but instead of scanning from the start of the
        // CIGAR, resumes from the cursor's saved position. After the call,
        // the cursor is positioned at the end of the walked range, ready
        // for the next call.
        //
        // This matches hifiasm's pattern in extract_sub_cigar_hc where
        // (ck, xk, yk) are saved between successive calls for the same
        // overlap across sliding windows.
        //
        // If queryStart < cursor.xk, the cursor seeks backward (like
        // hifiasm's `while (ck > 0 && xk > s)` loop). For the common
        // forward-sliding-window case, no backward seek is needed.
        //
        // Callback signature:
        //   void(uint8_t op, uint32_t len, uint64_t read0Pos, uint64_t read1Pos)
        template<class F>
        void walkRangeWithCursor(
            Cursor& cur,
            uint64_t queryStart,
            uint64_t queryEnd,
            F&& f) const
        {
            DINARA_ASSERT(cur.valid());

            // Seek backward if the cursor is past the start of the range.
            // Walk tokens in reverse, undoing their position effects.
            while (cur.tokenIndex > 0 && cur.xk > queryStart) {
                cur.tokenIndex--;
                const uint8_t op = cur.tokens[cur.tokenIndex].op();
                const uint32_t len = cur.tokens[cur.tokenIndex].len();
                if(op != 2) cur.xk -= len; // not insertion → undo read0
                if(op != 3) cur.yk -= len; // not deletion → undo read1
            }

            // Walk forward through the range [queryStart, queryEnd).
            while (cur.tokenIndex < cur.tokenCount) {
                // Coalesce consecutive tokens with the same op.
                const uint8_t op = cur.tokens[cur.tokenIndex].op();
                uint32_t totalLen = cur.tokens[cur.tokenIndex].len();
                uint32_t peekIdx = cur.tokenIndex + 1;
                while (peekIdx < cur.tokenCount && cur.tokens[peekIdx].op() == op) {
                    totalLen += cur.tokens[peekIdx].len();
                    peekIdx++;
                }

                const uint64_t xkEnd = cur.xk + ((op != 2) ? totalLen : 0);
                const uint64_t ykEnd = cur.yk + ((op != 3) ? totalLen : 0);

                if(op == 2) {
                    // Insertion: report if anchored within the query range.
                    if(cur.xk >= queryStart && cur.xk < queryEnd) {
                        f(op, totalLen, cur.xk, cur.yk);
                    }
                } else {
                    // Match, mismatch, or deletion: consumes read0.
                    if(xkEnd > queryStart && cur.xk < queryEnd) {
                        const uint64_t clipStart = (cur.xk < queryStart) ? queryStart : cur.xk;
                        const uint64_t clipEnd = (xkEnd > queryEnd) ? queryEnd : xkEnd;
                        const uint32_t clipLen = uint32_t(clipEnd - clipStart);
                        const uint64_t skipBases = clipStart - cur.xk;

                        uint64_t adjYk = cur.yk;
                        if(op == 0 || op == 1) {
                            adjYk += skipBases;
                        }

                        f(op, clipLen, clipStart, adjYk);
                    }
                }

                // Past the range — stop. Don't advance the cursor past
                // the current op so the next call can resume here.
                if(xkEnd >= queryEnd && op != 2) break;

                // Advance cursor past the coalesced tokens.
                cur.tokenIndex = peekIdx;
                cur.xk = xkEnd;
                cur.yk = ykEnd;
            }
        }

        // Map a single read0 (query) position to the corresponding read1
        // (target) position. Returns the read1 position, or uint64_t(-1)
        // if the query position falls inside a deletion (no corresponding
        // target base).
        //
        // This is the equivalent of hifiasm's xk/yk coordinate tracking.
        uint64_t queryToTarget(
            uint32_t cigarId,
            uint64_t read0Start,
            uint64_t read1Start,
            uint64_t queryPos) const
        {
            const auto tokens = getTokens(cigarId);
            uint64_t xk = read0Start;
            uint64_t yk = read1Start;
            for (uint32_t i = 0; i < tokens.size(); ) {
                const uint8_t op = tokens[i].op();
                uint32_t totalLen = tokens[i].len();
                i++;
                while (i < tokens.size() && tokens[i].op() == op) {
                    totalLen += tokens[i].len();
                    i++;
                }
                switch(op) {
                    case 0: case 1: // match/mismatch
                        if(queryPos >= xk && queryPos < xk + totalLen) {
                            return yk + (queryPos - xk);
                        }
                        xk += totalLen; yk += totalLen;
                        break;
                    case 2: // insertion — read0 doesn't advance
                        yk += totalLen;
                        break;
                    case 3: // deletion — read1 doesn't advance
                        if(queryPos >= xk && queryPos < xk + totalLen) {
                            return uint64_t(-1); // inside a deletion
                        }
                        xk += totalLen;
                        break;
                }
            }
            return uint64_t(-1); // out of range
        }

        // Map a single read1 (target) position to the corresponding read0
        // (query) position. Returns uint64_t(-1) if the target position
        // falls inside an insertion.
        uint64_t targetToQuery(
            uint32_t cigarId,
            uint64_t read0Start,
            uint64_t read1Start,
            uint64_t targetPos) const
        {
            const auto tokens = getTokens(cigarId);
            uint64_t xk = read0Start;
            uint64_t yk = read1Start;
            for (uint32_t i = 0; i < tokens.size(); ) {
                const uint8_t op = tokens[i].op();
                uint32_t totalLen = tokens[i].len();
                i++;
                while (i < tokens.size() && tokens[i].op() == op) {
                    totalLen += tokens[i].len();
                    i++;
                }
                switch(op) {
                    case 0: case 1:
                        if(targetPos >= yk && targetPos < yk + totalLen) {
                            return xk + (targetPos - yk);
                        }
                        xk += totalLen; yk += totalLen;
                        break;
                    case 2: // insertion — read0 doesn't advance
                        if(targetPos >= yk && targetPos < yk + totalLen) {
                            return uint64_t(-1); // inside an insertion
                        }
                        yk += totalLen;
                        break;
                    case 3: // deletion — read1 doesn't advance
                        xk += totalLen;
                        break;
                }
            }
            return uint64_t(-1);
        }

        // Merge another store into this one. Token data is appended and
        // index entries are offset-adjusted. Returns the starting cigar ID
        // in this store for the merged entries.
        uint32_t merge(const OverlapCigarStore& other) {
            const uint32_t baseId = uint32_t(index.size());
            const uint64_t arenaOffset = arena.size();

            arena.insert(arena.end(), other.arena.begin(), other.arena.end());
            index.reserve(index.size() + other.index.size());
            for (const auto& entry : other.index) {
                index.push_back({ entry.offset + arenaOffset, entry.count });
            }
            return baseId;
        }

        size_t tokenCount() const { return arena.size(); }
        size_t alignmentCount() const { return index.size(); }

        // Memory usage in bytes (approximate).
        size_t memoryUsage() const {
            return arena.capacity() * sizeof(CigarToken)
                 + index.capacity() * sizeof(IndexEntry);
        }
    };

} // namespace dinara

#endif
