#ifndef DINARA_PAF_IMPORT_HPP
#define DINARA_PAF_IMPORT_HPP

// Pure, allocation-free helpers for canonicalizing and deduplicating overlap
// entries collected from hifiasm's in-memory overlaps.
//
// These functions are deliberately independent of the Assembler so they can be
// unit tested in isolation. The Assembler side (importAlignmentCandidatesFromMemory)
// resolves hifiasm read ids to ReadIds, builds PafEntry records via makePafEntry,
// and merges them using the helpers below.
//
// Design notes:
//  - Canonical keying: q* always refers to min(readId), t* to max(readId).
//  - Deduplication is deterministic and independent of thread count/scheduling.

#include "ReadId.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace dinara {

// The overlap interval a PAF record agreed on, stored per canonical read pair.
// Coordinates are half-open base positions. q* refers to the smaller read id,
// t* to the larger read id (forward strand, matching Alignment::ts/te).
struct PafCandidateInterval {
    uint32_t qStart = 0;
    uint32_t qEnd = 0;
    uint32_t tStart = 0;
    uint32_t tEnd = 0;
    uint32_t blockLen = 0;   // PAF alignment block length (used as chain score).
    bool isSameStrand = true;
};

// A merged candidate entry: the canonical key plus its interval.
// key packs (readId0 << 32) | readId1 with readId0 < readId1.
//
// sourceIndex optionally identifies the record this entry came from, so a
// caller can recover per-record data (e.g. the hifiasm CIGAR) for the entry
// that survives dedup. It is not used for keying, ordering, or dedup and is
// left as uint64_t(-1) by the PAF-file path, which has no such side data.
struct PafEntry {
    uint64_t key = 0;
    PafCandidateInterval iv;
    uint64_t sourceIndex = uint64_t(-1);
};

// Canonicalize a resolved overlap into a PafEntry. readId0/readId1 are the
// resolved query/target ids (already known distinct). The interval is stored so
// q* always refers to min(id) and t* to max(id), matching the forward-coordinate
// convention. blockLen and strand are copied through.
inline PafEntry makePafEntry(
    ReadId readId0, ReadId readId1,
    uint32_t qStart, uint32_t qEnd,
    uint32_t tStart, uint32_t tEnd,
    uint32_t blockLen, bool isSameStrand)
{
    PafEntry entry;
    entry.iv.isSameStrand = isSameStrand;
    entry.iv.blockLen = blockLen;
    if(readId0 < readId1) {
        entry.iv.qStart = qStart;
        entry.iv.qEnd   = qEnd;
        entry.iv.tStart = tStart;
        entry.iv.tEnd   = tEnd;
        entry.key = (uint64_t(readId0) << 32) | uint64_t(readId1);
    } else {
        entry.iv.qStart = tStart;
        entry.iv.qEnd   = tEnd;
        entry.iv.tStart = qStart;
        entry.iv.tEnd   = qEnd;
        entry.key = (uint64_t(readId1) << 32) | uint64_t(readId0);
    }
    return entry;
}


// Strict ordering used to make deduplication deterministic:
//   key ascending, then strand (same before diff), then blockLen DESCENDING (so
//   the longest overlap sorts first for each (key, strand)), then remaining
//   fields to break ties reproducibly.
inline bool pafEntryLess(const PafEntry& a, const PafEntry& b)
{
    if(a.key != b.key) return a.key < b.key;
    // same-strand (true) sorts before reverse (false).
    if(a.iv.isSameStrand != b.iv.isSameStrand) return int(a.iv.isSameStrand) > int(b.iv.isSameStrand);
    if(a.iv.blockLen != b.iv.blockLen) return a.iv.blockLen > b.iv.blockLen;
    if(a.iv.qStart != b.iv.qStart) return a.iv.qStart < b.iv.qStart;
    if(a.iv.qEnd != b.iv.qEnd) return a.iv.qEnd < b.iv.qEnd;
    if(a.iv.tStart != b.iv.tStart) return a.iv.tStart < b.iv.tStart;
    return a.iv.tEnd < b.iv.tEnd;
}


// Sort `entries` in place and collapse duplicates by (key, strand), keeping the
// entry with the largest blockLen (the longest overlap) for each (key, strand).
// A read pair that overlaps in both orientations therefore keeps up to two
// entries (one +, one -). After this call `entries` is in ascending key order,
// same-strand before reverse within a key (deterministic).
inline void dedupPafEntriesKeepLongest(std::vector<PafEntry>& entries)
{
    std::sort(entries.begin(), entries.end(), pafEntryLess);
    size_t w = 0;
    for(size_t r = 0; r < entries.size(); ++r) {
        const bool distinct = (w == 0) ||
            (entries[r].key != entries[w - 1].key) ||
            (entries[r].iv.isSameStrand != entries[w - 1].iv.isSameStrand);
        if(distinct) {
            entries[w++] = entries[r];   // First per (key, strand) = longest.
        }
    }
    entries.resize(w);
}

} // namespace dinara

#endif
