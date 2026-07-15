#pragma once

// A minimal POD (trivially copyable) 256-bit unsigned integer.
//
// Purpose: serve as the packed KmerId type when DINARA_LONG_MARKERS is
// enabled (Kmer capacity 128 -> KmerId needs 2*128 = 256 bits).
//
// Rationale: dinara packs a k-mer into a single integer via
// ShortBaseSequence::id(k), and that packed value (KmerId) is:
//   - stored in memory-mapped containers (MemoryMapped::VectorOfVectors<KmerId>)
//   - hashed by raw bytes (MurmurHash2(&kmerId, sizeof(kmerId), ...))
//   - used as a hash-table / std::map key
// boost::multiprecision::uint256_t is NOT trivially copyable and has no
// guaranteed flat byte layout, so it is unsound for all three uses above.
// This struct is a flat little-endian array of four uint64_t words, so it is
// trivially copyable, memory-mappable, and safe to hash by raw bytes.
//
// Only the operations actually used by dinara's KmerId consumers are provided.
// Word order: w[0] is the least significant 64 bits.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <ostream>

namespace dinara {

    class Uint256 {
    public:
        uint64_t w[4];

        // Zero.
        Uint256() : w{0, 0, 0, 0} {}

        // From small unsigned/signed integers.
        Uint256(uint64_t x) : w{x, 0, 0, 0} {}
        Uint256(int x) : w{uint64_t(x), 0, 0, 0} {}
        Uint256(unsigned int x) : w{uint64_t(x), 0, 0, 0} {}

        // From a 128-bit integer (the Kmer word type when capacity == 128).
        Uint256(__uint128_t x) :
            w{uint64_t(x), uint64_t(x >> 64), 0, 0} {}

        // Truncating conversion to __uint128_t. Used by ShortBaseSequence's
        // id-constructor: data[i] = Int((id & mask) << shift), Int = __uint128_t.
        explicit operator __uint128_t() const {
            return (__uint128_t(w[1]) << 64) | __uint128_t(w[0]);
        }
        explicit operator uint64_t() const { return w[0]; }
        explicit operator bool() const {
            return (w[0] | w[1] | w[2] | w[3]) != 0;
        }

        // Bitwise.
        Uint256 operator~() const {
            Uint256 r;
            r.w[0] = ~w[0]; r.w[1] = ~w[1]; r.w[2] = ~w[2]; r.w[3] = ~w[3];
            return r;
        }
        Uint256 operator&(const Uint256& o) const {
            Uint256 r;
            r.w[0] = w[0] & o.w[0]; r.w[1] = w[1] & o.w[1];
            r.w[2] = w[2] & o.w[2]; r.w[3] = w[3] & o.w[3];
            return r;
        }
        Uint256 operator|(const Uint256& o) const {
            Uint256 r;
            r.w[0] = w[0] | o.w[0]; r.w[1] = w[1] | o.w[1];
            r.w[2] = w[2] | o.w[2]; r.w[3] = w[3] | o.w[3];
            return r;
        }
        Uint256 operator^(const Uint256& o) const {
            Uint256 r;
            r.w[0] = w[0] ^ o.w[0]; r.w[1] = w[1] ^ o.w[1];
            r.w[2] = w[2] ^ o.w[2]; r.w[3] = w[3] ^ o.w[3];
            return r;
        }

        // Left shift by s bits (0 <= s < 256; s >= 256 yields 0).
        Uint256 operator<<(uint64_t s) const {
            Uint256 r;
            if(s >= 256) {
                return r; // zero
            }
            const uint64_t wordShift = s >> 6;      // s / 64
            const uint64_t bitShift = s & 63;       // s % 64
            if(bitShift == 0) {
                for(int i = 3; i >= int(wordShift); --i) {
                    r.w[i] = w[i - wordShift];
                }
            } else {
                for(int i = 3; i >= 0; --i) {
                    const int src = i - int(wordShift);
                    if(src < 0) {
                        continue;
                    }
                    uint64_t v = w[src] << bitShift;
                    if(src - 1 >= 0) {
                        v |= w[src - 1] >> (64 - bitShift);
                    }
                    r.w[i] = v;
                }
            }
            return r;
        }

        // Right shift by s bits (0 <= s < 256; s >= 256 yields 0).
        Uint256 operator>>(uint64_t s) const {
            Uint256 r;
            if(s >= 256) {
                return r; // zero
            }
            const uint64_t wordShift = s >> 6;
            const uint64_t bitShift = s & 63;
            if(bitShift == 0) {
                for(int i = 0; i + int(wordShift) < 4; ++i) {
                    r.w[i] = w[i + wordShift];
                }
            } else {
                for(int i = 0; i < 4; ++i) {
                    const int src = i + int(wordShift);
                    if(src > 3) {
                        continue;
                    }
                    uint64_t v = w[src] >> bitShift;
                    if(src + 1 <= 3) {
                        v |= w[src + 1] << (64 - bitShift);
                    }
                    r.w[i] = v;
                }
            }
            return r;
        }

        // Addition (wraps mod 2^256), used for mask arithmetic and small sums.
        Uint256 operator+(const Uint256& o) const {
            Uint256 r;
            unsigned __int128 carry = 0;
            for(int i = 0; i < 4; ++i) {
                const unsigned __int128 s =
                    (unsigned __int128)w[i] + o.w[i] + carry;
                r.w[i] = uint64_t(s);
                carry = s >> 64;
            }
            return r;
        }

        // Subtraction (wraps mod 2^256). Used for mask = (1<<n) - 1.
        Uint256 operator-(const Uint256& o) const {
            Uint256 r;
            unsigned __int128 borrow = 0;
            for(int i = 0; i < 4; ++i) {
                const unsigned __int128 d =
                    (unsigned __int128)w[i] - o.w[i] - borrow;
                r.w[i] = uint64_t(d);
                borrow = (d >> 64) & 1; // 1 if underflow
            }
            return r;
        }

        // Comparisons.
        bool operator==(const Uint256& o) const {
            return w[0] == o.w[0] && w[1] == o.w[1] &&
                   w[2] == o.w[2] && w[3] == o.w[3];
        }
        bool operator!=(const Uint256& o) const { return !(*this == o); }
        bool operator<(const Uint256& o) const {
            for(int i = 3; i >= 0; --i) {
                if(w[i] != o.w[i]) {
                    return w[i] < o.w[i];
                }
            }
            return false;
        }
        bool operator>(const Uint256& o) const { return o < *this; }
        bool operator<=(const Uint256& o) const { return !(o < *this); }
        bool operator>=(const Uint256& o) const { return !(*this < o); }
    };

    static_assert(sizeof(Uint256) == 32, "Uint256 must be 32 bytes.");

    // Stream as hexadecimal (most significant word first), matching the
    // debug-output expectations of the KmerId consumers.
    inline std::ostream& operator<<(std::ostream& s, const Uint256& x) {
        const std::ios_base::fmtflags f = s.flags();
        s << std::hex;
        bool started = false;
        for(int i = 3; i >= 0; --i) {
            if(started) {
                // Zero-pad continuation words to 16 hex digits.
                const char fill = s.fill('0');
                s.width(16);
                s << x.w[i];
                s.fill(fill);
            } else if(x.w[i] != 0 || i == 0) {
                s << x.w[i];
                started = true;
            }
        }
        s.flags(f);
        return s;
    }
}

// Hash support so Uint256 can be used as a std::unordered_* key and by
// generic code that expects std::hash. dinara also has an explicit
// KmerIdHasher (MurmurHash over raw bytes); this mirrors that intent.
namespace std {
    template<> struct hash<dinara::Uint256> {
        size_t operator()(const dinara::Uint256& x) const {
            // FNV-1a over the four words. Deterministic and byte-layout stable.
            uint64_t h = 1469598103934665603ULL;
            for(int i = 0; i < 4; ++i) {
                h ^= x.w[i];
                h *= 1099511628211ULL;
            }
            return size_t(h);
        }
    };
}
