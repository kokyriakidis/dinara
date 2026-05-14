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
// across consecutive tokens with the same op.
//
// All tokens for all alignments live in a single flat arena.
// Each overlap stores its own (cigarOffset, cigarTokenCount) into the
// arena directly in AlignmentInfo — no separate index structure.

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

        void clear() { arena.clear(); }

        void reserve(size_t arenaSize) { arena.reserve(arenaSize); }

        // Begin a new alignment's CIGAR. Returns the arena offset where
        // this alignment's tokens will start. Store this in
        // AlignmentInfo::cigarOffset.
        uint32_t beginAlignment() {
            return uint32_t(arena.size());
        }

        // Return the number of tokens pushed since the given offset.
        // Call after all pushOp calls for an alignment; store the result
        // in AlignmentInfo::cigarTokenCount.
        uint32_t tokensSince(uint32_t offset) const {
            DINARA_ASSERT(arena.size() >= offset);
            return uint32_t(arena.size()) - offset;
        }

        // Push a single CIGAR operation. Runs exceeding MAX_LEN are
        // automatically split into multiple tokens.
        void pushOp(uint8_t op, uint32_t length) {
            DINARA_ASSERT(op < 4);
            while (length > CigarToken::MAX_LEN) {
                arena.emplace_back(op, CigarToken::MAX_LEN);
                length -= CigarToken::MAX_LEN;
            }
            if (length > 0) {
                arena.emplace_back(op, uint16_t(length));
            }
        }

        void pushMatch(uint32_t length) { pushOp(0, length); }
        void pushMismatch(uint32_t length) { pushOp(1, length); }
        void pushInsertion(uint32_t length) { pushOp(2, length); }
        void pushDeletion(uint32_t length) { pushOp(3, length); }

        // Retrieve the token slice for a given (offset, count).
        span<const CigarToken> getTokens(uint32_t offset, uint32_t count) const {
            if (offset == uint32_t(-1) || count == uint32_t(-1)) return {};
            DINARA_ASSERT(uint64_t(offset) + count <= arena.size());
            return { arena.data() + offset, count };
        }

        // Resumable cursor for walking a CIGAR across successive query
        // ranges. Matches hifiasm's ovlp_cur_xoff/yoff/coff pattern.
        struct Cursor {
            uint32_t tokenIndex = 0;
            uint64_t xk = 0;         // Read0 (query) position.
            uint64_t yk = 0;         // Read1 (target) position.

            const CigarToken* tokens = nullptr;
            uint32_t tokenCount = 0;

            void reset(uint32_t offset, uint32_t count,
                       uint64_t read0Start, uint64_t read1Start,
                       const OverlapCigarStore& store) {
                auto s = store.getTokens(offset, count);
                tokens = s.data();
                tokenCount = uint32_t(s.size());
                tokenIndex = 0;
                xk = read0Start;
                yk = read1Start;
            }

            bool valid() const { return tokens != nullptr && tokenCount > 0; }
        };

        // Iterate the CIGAR, coalescing consecutive tokens with the same op.
        // Callback: void(uint8_t op, uint32_t length)
        template<class F>
        static void forEachOp(span<const CigarToken> tokens, F&& f) {
            uint32_t i = 0;
            while (i < tokens.size()) {
                const uint8_t op = tokens[i].op();
                uint32_t totalLen = tokens[i].len();
                i++;
                while (i < tokens.size() && tokens[i].op() == op) {
                    totalLen += tokens[i].len();
                    i++;
                }
                f(op, totalLen);
            }
        }

        // Convenience overload taking (offset, count).
        template<class F>
        void forEachOp(uint32_t offset, uint32_t count, F&& f) const {
            forEachOp(getTokens(offset, count), std::forward<F>(f));
        }

        // Walk with full position tracking on both reads.
        // Callback: void(uint8_t op, uint32_t len, uint64_t read0Pos, uint64_t read1Pos)
        template<class F>
        void forEachOpWithPositions(
            uint32_t offset, uint32_t count,
            uint64_t read0Start, uint64_t read1Start,
            F&& f) const
        {
            const auto tokens = getTokens(offset, count);
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

        // Walk only the sub-range [queryStart, queryEnd) on read0.
        // Ops are clipped to the range boundaries.
        // Callback: void(uint8_t op, uint32_t len, uint64_t read0Pos, uint64_t read1Pos)
        template<class F>
        void walkRange(
            uint32_t offset, uint32_t count,
            uint64_t read0Start, uint64_t read1Start,
            uint64_t queryStart, uint64_t queryEnd,
            F&& f) const
        {
            const auto tokens = getTokens(offset, count);
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
                    if(xk >= queryStart && xk < queryEnd) {
                        f(op, totalLen, xk, yk);
                    }
                } else {
                    if(xkEnd > queryStart && xk < queryEnd) {
                        const uint64_t clipStart = (xk < queryStart) ? queryStart : xk;
                        const uint64_t clipEnd = (xkEnd > queryEnd) ? queryEnd : xkEnd;
                        const uint32_t clipLen = uint32_t(clipEnd - clipStart);
                        const uint64_t skipBases = clipStart - xk;

                        uint64_t adjYk = yk;
                        if(op == 0 || op == 1) adjYk += skipBases;

                        f(op, clipLen, clipStart, adjYk);
                    }
                }

                if(xkEnd >= queryEnd && op != 2) break;

                xk = xkEnd;
                yk = ykEnd;
            }
        }

        // Walk a query range using a resumable cursor.
        // Callback: void(uint8_t op, uint32_t len, uint64_t read0Pos, uint64_t read1Pos)
        template<class F>
        void walkRangeWithCursor(
            Cursor& cur,
            uint64_t queryStart,
            uint64_t queryEnd,
            F&& f) const
        {
            DINARA_ASSERT(cur.valid());

            while (cur.tokenIndex > 0 && cur.xk > queryStart) {
                cur.tokenIndex--;
                const uint8_t op = cur.tokens[cur.tokenIndex].op();
                const uint32_t len = cur.tokens[cur.tokenIndex].len();
                if(op != 2) cur.xk -= len;
                if(op != 3) cur.yk -= len;
            }

            while (cur.tokenIndex < cur.tokenCount) {
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
                    if(cur.xk >= queryStart && cur.xk < queryEnd) {
                        f(op, totalLen, cur.xk, cur.yk);
                    }
                } else {
                    if(xkEnd > queryStart && cur.xk < queryEnd) {
                        const uint64_t clipStart = (cur.xk < queryStart) ? queryStart : cur.xk;
                        const uint64_t clipEnd = (xkEnd > queryEnd) ? queryEnd : xkEnd;
                        const uint32_t clipLen = uint32_t(clipEnd - clipStart);
                        const uint64_t skipBases = clipStart - cur.xk;

                        uint64_t adjYk = cur.yk;
                        if(op == 0 || op == 1) adjYk += skipBases;

                        f(op, clipLen, clipStart, adjYk);
                    }
                }

                if(xkEnd >= queryEnd && op != 2) break;

                cur.tokenIndex = peekIdx;
                cur.xk = xkEnd;
                cur.yk = ykEnd;
            }
        }

        // Map read0 position to read1 position.
        // Returns uint64_t(-1) if inside a deletion.
        uint64_t queryToTarget(
            uint32_t offset, uint32_t count,
            uint64_t read0Start, uint64_t read1Start,
            uint64_t queryPos) const
        {
            const auto tokens = getTokens(offset, count);
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
                        if(queryPos >= xk && queryPos < xk + totalLen)
                            return yk + (queryPos - xk);
                        xk += totalLen; yk += totalLen;
                        break;
                    case 2:
                        yk += totalLen;
                        break;
                    case 3:
                        if(queryPos >= xk && queryPos < xk + totalLen)
                            return uint64_t(-1);
                        xk += totalLen;
                        break;
                }
            }
            return uint64_t(-1);
        }

        // Map read1 position to read0 position.
        // Returns uint64_t(-1) if inside an insertion.
        uint64_t targetToQuery(
            uint32_t offset, uint32_t count,
            uint64_t read0Start, uint64_t read1Start,
            uint64_t targetPos) const
        {
            const auto tokens = getTokens(offset, count);
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
                        if(targetPos >= yk && targetPos < yk + totalLen)
                            return xk + (targetPos - yk);
                        xk += totalLen; yk += totalLen;
                        break;
                    case 2:
                        if(targetPos >= yk && targetPos < yk + totalLen)
                            return uint64_t(-1);
                        yk += totalLen;
                        break;
                    case 3:
                        xk += totalLen;
                        break;
                }
            }
            return uint64_t(-1);
        }

        // Merge another store's arena into this one. Returns the offset
        // that must be added to all cigarOffset values from the other store.
        uint32_t merge(const OverlapCigarStore& other) {
            const uint32_t arenaOffset = uint32_t(arena.size());
            arena.insert(arena.end(), other.arena.begin(), other.arena.end());
            return arenaOffset;
        }

        size_t tokenCount() const { return arena.size(); }

        // Memory usage in bytes.
        size_t memoryUsage() const {
            return arena.capacity() * sizeof(CigarToken);
        }
    };

} // namespace dinara

#endif
