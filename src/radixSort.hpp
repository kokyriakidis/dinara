#ifndef DINARA_RADIX_SORT_HPP
#define DINARA_RADIX_SORT_HPP

/// @file radixSort.hpp
/// @brief In-place MSD radix sort with insertion sort fallback.
///
/// Port of hifiasm's KRADIX_SORT_INIT (ksort.h) as a C++ template.
/// Sorts by an unsigned integer key extracted via a user-supplied functor.
///
/// Properties:
///   - In-place (no temporary buffer).
///   - MSD (most-significant-digit first), 8 bits per pass.
///   - Insertion sort for partitions <= 64 elements.
///   - Stack-allocated 256-bucket array per recursion level.
///
/// Usage:
///   dinara::radixSort(begin, end, [](const T& x) { return x.key; });

#include <cstddef>
#include <cstdint>

namespace dinara {

namespace detail {

static constexpr int RS_MIN_SIZE = 64;

/// Insertion sort fallback for small ranges.
/// Exact port of hifiasm's rs_insertsort (ksort.h).
template<typename T, typename KeyFn>
inline void rsInsertionSort(T* beg, T* end, KeyFn key) {
    for (T* i = beg + 1; i < end; ++i) {
        if (key(*i) < key(*(i - 1))) {
            T tmp = *i;
            T* j = i;
            for (; j > beg && key(tmp) < key(*(j - 1)); --j)
                *j = *(j - 1);
            *j = tmp;
        }
    }
}

/// MSD radix sort, 8 bits per pass. Recurses on each bucket.
/// Exact port of hifiasm's rs_sort (ksort.h), using pointer-based
/// bucket iteration to match the original's access pattern.
template<typename T, typename KeyFn>
void rsSort(T* beg, T* end, int n_bits, int s, KeyFn key) {
    const int size = 1 << n_bits;
    const int m = size - 1;

    struct Bucket { T *b, *e; };
    Bucket b[256]; // n_bits is always 8
    Bucket* be = b + size;

    for (Bucket* k = b; k != be; ++k) k->b = k->e = beg;
    for (T* i = beg; i != end; ++i) ++b[(key(*i) >> s) & m].e;
    for (Bucket* k = b + 1; k != be; ++k) {
        k->e += (k - 1)->e - beg;
        k->b = (k - 1)->e;
    }

    // In-place permutation via cycle-leader algorithm.
    for (Bucket* k = b; k != be;) {
        if (k->b != k->e) {
            Bucket* l;
            if ((l = b + ((key(*k->b) >> s) & m)) != k) {
                T tmp = *k->b, swap;
                do {
                    swap = tmp; tmp = *l->b; *l->b++ = swap;
                    l = b + ((key(tmp) >> s) & m);
                } while (l != k);
                *k->b++ = tmp;
            } else {
                ++k->b;
            }
        } else {
            ++k;
        }
    }

    // Restore b pointers for recursion.
    b[0].b = beg;
    for (Bucket* k = b + 1; k != be; ++k) k->b = (k - 1)->e;

    // Recurse on each bucket (next byte down).
    if (s) {
        const int nextS = s > n_bits ? s - n_bits : 0;
        for (Bucket* k = b; k != be; ++k) {
            if (k->e - k->b > RS_MIN_SIZE)
                rsSort(k->b, k->e, n_bits, nextS, key);
            else if (k->e - k->b > 1)
                rsInsertionSort(k->b, k->e, key);
        }
    }
}

} // namespace detail

/// In-place MSD radix sort.
///
/// @tparam T       Element type.
/// @tparam KeyFn   Functor returning an unsigned integer key for each element.
///                 The key size is deduced automatically (1, 2, 4, or 8 bytes).
/// @param beg      Pointer to first element.
/// @param end      Pointer past last element.
/// @param key      Key extraction functor.
///
/// Example:
///   radixSort(vec.data(), vec.data() + vec.size(),
///             [](const Ev& e) { return e.site; });
template<typename T, typename KeyFn>
inline void radixSort(T* beg, T* end, KeyFn key) {
    if (end - beg <= detail::RS_MIN_SIZE) {
        if (end - beg > 1)
            detail::rsInsertionSort(beg, end, key);
        return;
    }
    constexpr int keySize = sizeof(decltype(key(*beg)));
    detail::rsSort(beg, end, 8, keySize * 8 - 8, key);
}

} // namespace dinara

#endif // DINARA_RADIX_SORT_HPP
