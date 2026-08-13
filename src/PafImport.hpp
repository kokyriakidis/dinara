#ifndef DINARA_PAF_IMPORT_HPP
#define DINARA_PAF_IMPORT_HPP

// Pure, allocation-free helpers for parsing a PAF overlap file in parallel.
//
// These functions are deliberately independent of the Assembler so they can be
// unit tested in isolation. The Assembler side (importAlignmentCandidatesFromPaf)
// only mmaps the file, splits it into line-aligned chunks, resolves read names to
// ReadIds, and merges the per-thread results using the helpers below.
//
// Design notes:
//  - Parsing is tab/whitespace tolerant to match the previous stringstream-based
//    behavior (operator>> splits on any whitespace run).
//  - Integer fields use std::from_chars: no locale, no allocation, no exceptions.
//  - Names are returned as span<const char> into the mmap; no copies are made.
//  - Chunk ranges always start and end on a line boundary, so every line is
//    processed exactly once regardless of how many chunks are used.
//  - Deduplication is deterministic and independent of thread count/scheduling.

#include "ReadId.hpp"
#include "span.hpp"

#include <charconv>
#include <cstdint>
#include <utility>
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

// A parsed PAF line, before read-name resolution.
struct PafRecord {
    span<const char> qName;
    span<const char> tName;
    uint64_t qLen = 0;
    uint64_t qStart = 0;
    uint64_t qEnd = 0;
    uint64_t tLen = 0;
    uint64_t tStart = 0;
    uint64_t tEnd = 0;
    uint64_t mapQ = 0;
    uint64_t alignLen = 0;
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

// Both orientations a read pair may overlap in. A pair can legitimately appear
// as both a same-strand (+) and a reverse (-) overlap (e.g. inverted repeats),
// so we keep them separately instead of collapsing to one. Each is optional;
// the longest overlap is kept per orientation.
struct PafPairIntervals {
    bool haveSame = false;
    bool haveDiff = false;
    PafCandidateInterval same;   // valid iff haveSame
    PafCandidateInterval diff;   // valid iff haveDiff

    // Return the interval for the requested orientation, or nullptr if absent.
    const PafCandidateInterval* get(bool isSameStrand) const {
        if(isSameStrand) return haveSame ? &same : nullptr;
        return haveDiff ? &diff : nullptr;
    }
};


// Return true if c is a field delimiter (space, tab, CR, LF).
inline bool pafIsDelim(char c)
{
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}


// Parse a single PAF line given as [lineBegin, lineEnd).
// Returns true and fills `out` if at least the 11 mandatory columns are present
// and the integer columns parse. Columns 12+ (e.g. cg:Z: CIGAR) are ignored.
// Makes no allocations; names point into the input buffer.
inline bool parsePafLine(const char* lineBegin, const char* lineEnd, PafRecord& out)
{
    const char* p = lineBegin;

    // Advance p past any leading delimiters, then return [tokBegin, tokEnd)
    // as the next whitespace-delimited token. Returns false if none remains.
    auto nextToken = [&](const char*& tokBegin, const char*& tokEnd) -> bool {
        while(p < lineEnd && pafIsDelim(*p)) ++p;
        if(p >= lineEnd) return false;
        tokBegin = p;
        while(p < lineEnd && !pafIsDelim(*p)) ++p;
        tokEnd = p;
        return true;
    };

    auto parseUint = [](const char* b, const char* e, uint64_t& v) -> bool {
        if(b >= e) return false;
        const auto res = std::from_chars(b, e, v);
        return res.ec == std::errc() && res.ptr == e;
    };

    const char* b;
    const char* e;

    // Col 1: query name.
    if(!nextToken(b, e)) return false;
    out.qName = span<const char>(b, size_t(e - b));

    // Col 2: query length.
    if(!nextToken(b, e) || !parseUint(b, e, out.qLen)) return false;
    // Col 3: query start.
    if(!nextToken(b, e) || !parseUint(b, e, out.qStart)) return false;
    // Col 4: query end.
    if(!nextToken(b, e) || !parseUint(b, e, out.qEnd)) return false;

    // Col 5: strand.
    if(!nextToken(b, e)) return false;
    if(e - b != 1 || (*b != '+' && *b != '-')) return false;
    out.isSameStrand = (*b == '+');

    // Col 6: target name.
    if(!nextToken(b, e)) return false;
    out.tName = span<const char>(b, size_t(e - b));

    // Col 7: target length.
    if(!nextToken(b, e) || !parseUint(b, e, out.tLen)) return false;
    // Col 8: target start.
    if(!nextToken(b, e) || !parseUint(b, e, out.tStart)) return false;
    // Col 9: target end.
    if(!nextToken(b, e) || !parseUint(b, e, out.tEnd)) return false;
    // Col 10: residue matches (mapQ position in minimap PAF is col 12; hifiasm
    // uses col 10 here as a numeric field, kept for interface compatibility).
    if(!nextToken(b, e) || !parseUint(b, e, out.mapQ)) return false;
    // Col 11: alignment block length.
    if(!nextToken(b, e) || !parseUint(b, e, out.alignLen)) return false;

    return true;
}


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


// Split a buffer of `size` bytes into `nChunks` byte ranges [start, end) that
// each begin and end on a line boundary, so concatenating the lines of all
// chunks reproduces every line of the file exactly once. A chunk may be empty
// (start == end) when lines are long relative to the chunk count.
//
// Rule: chunk i tentatively covers [i*size/n, (i+1)*size/n). The real start is
// advanced to just past the previous newline (except chunk 0, which starts at 0);
// the real end is the real start of the next chunk. The final chunk always ends
// at `size`.
inline std::vector<std::pair<size_t, size_t>> computePafChunkRanges(
    const char* data, size_t size, size_t nChunks)
{
    std::vector<std::pair<size_t, size_t>> ranges;
    if(nChunks == 0) nChunks = 1;
    if(size == 0) {
        return ranges;   // No lines.
    }

    // Compute a line-aligned start offset for each raw boundary.
    std::vector<size_t> starts(nChunks + 1, 0);
    starts[0] = 0;
    for(size_t i = 1; i < nChunks; ++i) {
        size_t pos = (size * i) / nChunks;
        if(pos >= size) {
            pos = size;
        } else {
            // Advance to the byte just after the next newline at or after pos.
            while(pos < size && data[pos] != '\n') ++pos;
            if(pos < size) ++pos;   // Step over the '\n'.
        }
        starts[i] = pos;
    }
    starts[nChunks] = size;

    // Make starts monotonically non-decreasing (they already are by construction)
    // and emit ranges.
    ranges.reserve(nChunks);
    for(size_t i = 0; i < nChunks; ++i) {
        size_t s = starts[i];
        size_t e = starts[i + 1];
        if(e < s) e = s;
        ranges.emplace_back(s, e);
    }
    return ranges;
}

} // namespace dinara

#endif
