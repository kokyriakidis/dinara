/**
 * @file test_paf_import.cpp
 * @brief Unit tests for the parallel PAF import helpers in src/PafImport.hpp.
 *
 * These tests exercise the pure, allocation-free building blocks used by
 * Assembler::importAlignmentCandidatesFromPaf:
 *   - parsePafLine: tokenization + from_chars integer parsing, malformed input.
 *   - computePafChunkRanges: line-aligned splitting for every chunk count, so
 *     each line is processed exactly once regardless of thread count.
 *   - makePafEntry / pafEntryLess / dedupPafEntriesKeepLongest: canonical keying
 *     and longest-overlap-wins deduplication, independent of input order.
 *
 * A key property verified here is chunk-count invariance: parsing the whole
 * buffer as N line-aligned chunks (for every N from 1..lines+2) and merging must
 * always yield the identical deduplicated result.
 */

#include "catch.hpp"
#include "PafImport.hpp"

#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

using namespace dinara;
using std::string;
using std::vector;

namespace {

// A tiny in-memory "read name -> id" table for tests, mirroring what
// Reads::getReadId does for the real importer.
struct NameTable {
    std::unordered_map<string, ReadId> ids;
    ReadId get(span<const char> name) const {
        string key(name.data(), name.size());
        auto it = ids.find(key);
        return it == ids.end() ? invalidReadId : it->second;
    }
};

// Parse an entire buffer as `nChunks` line-aligned chunks, resolve names, and
// return deduplicated entries. This mirrors the Assembler's parallel path but
// runs the chunks serially so tests are deterministic and easy to reason about.
vector<PafEntry> importBuffer(
    const string& buffer, const NameTable& names, size_t nChunks,
    uint64_t minBlock = 200)
{
    const char* data = buffer.data();
    const size_t size = buffer.size();
    const auto ranges = computePafChunkRanges(data, size, nChunks);

    vector<PafEntry> entries;
    for (const auto& [rBegin, rEnd] : ranges) {
        const char* p = data + rBegin;
        const char* const chunkEnd = data + rEnd;
        while (p < chunkEnd) {
            const char* lineEnd = p;
            while (lineEnd < chunkEnd && *lineEnd != '\n') ++lineEnd;
            PafRecord rec;
            if (parsePafLine(p, lineEnd, rec) && rec.alignLen >= minBlock) {
                const ReadId r0 = names.get(rec.qName);
                const ReadId r1 = names.get(rec.tName);
                if (r0 != invalidReadId && r1 != invalidReadId && r0 != r1) {
                    entries.push_back(makePafEntry(
                        r0, r1,
                        uint32_t(rec.qStart), uint32_t(rec.qEnd),
                        uint32_t(rec.tStart), uint32_t(rec.tEnd),
                        uint32_t(rec.alignLen), rec.isSameStrand));
                }
            }
            p = (lineEnd < chunkEnd) ? lineEnd + 1 : chunkEnd;
        }
    }
    dedupPafEntriesKeepLongest(entries);
    return entries;
}

// Count total lines (including the last unterminated one) in a buffer.
size_t countLines(const string& s) {
    if (s.empty()) return 0;
    size_t n = 0;
    for (char c : s) if (c == '\n') ++n;
    if (s.back() != '\n') ++n;
    return n;
}

} // namespace


TEST_CASE("parsePafLine: valid 11-column line", "[paf]") {
    const string line = "readA\t1000\t100\t900\t+\treadB\t1200\t50\t850\t700\t800";
    PafRecord rec;
    REQUIRE(parsePafLine(line.data(), line.data() + line.size(), rec));
    REQUIRE(string(rec.qName.data(), rec.qName.size()) == "readA");
    REQUIRE(rec.qLen == 1000);
    REQUIRE(rec.qStart == 100);
    REQUIRE(rec.qEnd == 900);
    REQUIRE(rec.isSameStrand == true);
    REQUIRE(string(rec.tName.data(), rec.tName.size()) == "readB");
    REQUIRE(rec.tLen == 1200);
    REQUIRE(rec.tStart == 50);
    REQUIRE(rec.tEnd == 850);
    REQUIRE(rec.mapQ == 700);
    REQUIRE(rec.alignLen == 800);
}

TEST_CASE("parsePafLine: extra columns (cg:Z: CIGAR) are ignored", "[paf]") {
    const string line =
        "readA\t1000\t100\t900\t-\treadB\t1200\t50\t850\t700\t800\t255\tcg:Z:800M";
    PafRecord rec;
    REQUIRE(parsePafLine(line.data(), line.data() + line.size(), rec));
    REQUIRE(rec.isSameStrand == false);   // '-' strand
    REQUIRE(rec.alignLen == 800);
}

TEST_CASE("parsePafLine: space-delimited also works", "[paf]") {
    const string line = "readA 1000 100 900 + readB 1200 50 850 700 800";
    PafRecord rec;
    REQUIRE(parsePafLine(line.data(), line.data() + line.size(), rec));
    REQUIRE(rec.qEnd == 900);
    REQUIRE(rec.alignLen == 800);
}

TEST_CASE("parsePafLine: trailing carriage return handled", "[paf]") {
    const string line = "readA\t1000\t100\t900\t+\treadB\t1200\t50\t850\t700\t800\r";
    PafRecord rec;
    REQUIRE(parsePafLine(line.data(), line.data() + line.size(), rec));
    REQUIRE(rec.alignLen == 800);
}

TEST_CASE("parsePafLine: malformed lines are rejected", "[paf]") {
    auto bad = [](const string& s) {
        PafRecord rec;
        return parsePafLine(s.data(), s.data() + s.size(), rec);
    };
    REQUIRE_FALSE(bad(""));                                   // empty
    REQUIRE_FALSE(bad("readA\t1000\t100"));                   // too few columns
    REQUIRE_FALSE(bad("readA\tXXX\t100\t900\t+\treadB\t1200\t50\t850\t700\t800")); // non-numeric len
    REQUIRE_FALSE(bad("readA\t1000\t100\t900\t*\treadB\t1200\t50\t850\t700\t800")); // bad strand
    REQUIRE_FALSE(bad("readA\t1000\t100\t900\t+\treadB\t1200\t50\t850\t700\t80x")); // trailing junk in int
}

TEST_CASE("makePafEntry: canonicalizes read pair and swaps intervals", "[paf]") {
    // readId0 > readId1 -> must swap so key uses (min<<32|max) and q* = min id.
    PafEntry e = makePafEntry(
        /*readId0=*/5, /*readId1=*/2,
        /*qStart=*/100, /*qEnd=*/900,   // these are read 5's coords
        /*tStart=*/50,  /*tEnd=*/850,   // these are read 2's coords
        /*blockLen=*/700, /*isSameStrand=*/true);
    REQUIRE(e.key == ((uint64_t(2) << 32) | 5));
    // After swap: q* refers to read 2 (=old target), t* to read 5 (=old query).
    REQUIRE(e.iv.qStart == 50);
    REQUIRE(e.iv.qEnd == 850);
    REQUIRE(e.iv.tStart == 100);
    REQUIRE(e.iv.tEnd == 900);
    REQUIRE(e.iv.blockLen == 700);
}

TEST_CASE("makePafEntry: reciprocal A->B and B->A produce identical entry", "[paf]") {
    // A=5, B=2. Record A->B: q=A coords, t=B coords.
    PafEntry ab = makePafEntry(5, 2, 100, 900, 50, 850, 700, true);
    // Record B->A: q=B coords, t=A coords.
    PafEntry ba = makePafEntry(2, 5, 50, 850, 100, 900, 700, true);
    REQUIRE(ab.key == ba.key);
    REQUIRE(ab.iv.qStart == ba.iv.qStart);
    REQUIRE(ab.iv.qEnd == ba.iv.qEnd);
    REQUIRE(ab.iv.tStart == ba.iv.tStart);
    REQUIRE(ab.iv.tEnd == ba.iv.tEnd);
}

TEST_CASE("dedupPafEntriesKeepLongest: keeps the longest block per pair", "[paf]") {
    vector<PafEntry> entries = {
        makePafEntry(1, 2, 0, 300, 0, 300, 300, true),
        makePafEntry(1, 2, 0, 700, 0, 700, 700, true),   // longest for (1,2)
        makePafEntry(1, 2, 0, 150, 0, 150, 150, true),
        makePafEntry(3, 4, 0, 500, 0, 500, 500, false),
    };
    dedupPafEntriesKeepLongest(entries);
    REQUIRE(entries.size() == 2);
    // Sorted by key ascending: (1,2) then (3,4).
    REQUIRE(entries[0].key == ((uint64_t(1) << 32) | 2));
    REQUIRE(entries[0].iv.blockLen == 700);
    REQUIRE(entries[0].iv.qEnd == 700);
    REQUIRE(entries[1].key == ((uint64_t(3) << 32) | 4));
    REQUIRE(entries[1].iv.blockLen == 500);
}

TEST_CASE("dedupPafEntriesKeepLongest: result is order-independent", "[paf]") {
    auto build = [](bool reversed) {
        vector<PafEntry> v = {
            makePafEntry(1, 2, 0, 300, 0, 300, 300, true),
            makePafEntry(1, 2, 0, 700, 0, 700, 700, true),
            makePafEntry(2, 5, 0, 400, 0, 400, 400, false),
        };
        if (reversed) std::reverse(v.begin(), v.end());
        dedupPafEntriesKeepLongest(v);
        return v;
    };
    const auto a = build(false);
    const auto b = build(true);
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        REQUIRE(a[i].key == b[i].key);
        REQUIRE(a[i].iv.blockLen == b[i].iv.blockLen);
        REQUIRE(a[i].iv.qStart == b[i].iv.qStart);
        REQUIRE(a[i].iv.qEnd == b[i].iv.qEnd);
    }
}

TEST_CASE("computePafChunkRanges: empty buffer yields no ranges", "[paf]") {
    REQUIRE(computePafChunkRanges(nullptr, 0, 4).empty());
}

TEST_CASE("computePafChunkRanges: ranges are contiguous and line-aligned", "[paf]") {
    const string buf =
        "l1\n"
        "line2\n"
        "third-line\n"
        "l4\n"
        "last";   // no trailing newline
    for (size_t n = 1; n <= 8; ++n) {
        const auto ranges = computePafChunkRanges(buf.data(), buf.size(), n);
        REQUIRE_FALSE(ranges.empty());
        // Contiguity: first starts at 0, last ends at size, each end == next start.
        REQUIRE(ranges.front().first == 0);
        REQUIRE(ranges.back().second == buf.size());
        for (size_t i = 1; i < ranges.size(); ++i) {
            REQUIRE(ranges[i - 1].second == ranges[i].first);
        }
        // Every non-empty range must start at buffer start or just after a '\n',
        // and end at buffer end or just after a '\n'.
        for (const auto& [s, e] : ranges) {
            if (s != e) {
                REQUIRE((s == 0 || buf[s - 1] == '\n'));
                REQUIRE((e == buf.size() || buf[e - 1] == '\n'));
            }
        }
    }
}

TEST_CASE("import: chunk-count invariance over a realistic PAF", "[paf]") {
    NameTable names;
    names.ids = {{"r0", 0}, {"r1", 1}, {"r2", 2}, {"r3", 3}, {"r4", 4}};

    // Mixed strands, both orderings, duplicates, a short overlap (<200, dropped),
    // a self-overlap (dropped), an unknown read (dropped), a blank line, and a
    // malformed line.
    const string buf =
        "r0\t1000\t100\t900\t+\tr1\t1000\t50\t850\t700\t800\n"      // (0,1) block 800
        "r1\t1000\t60\t860\t+\tr0\t1000\t110\t910\t690\t790\n"      // (0,1) reciprocal block 790
        "r0\t1000\t0\t500\t-\tr2\t1000\t20\t480\t400\t490\n"        // (0,2) block 490
        "r2\t1000\t30\t520\t-\tr0\t1000\t40\t600\t500\t560\n"       // (0,2) dup, block 560 (longer)
        "r3\t1000\t10\t250\t+\tr4\t1000\t10\t250\t210\t150\n"       // (3,4) alignLen 150 < 200 -> drop
        "r1\t1000\t0\t400\t+\tr1\t1000\t0\t400\t380\t400\n"         // self overlap -> drop
        "rX\t1000\t0\t400\t+\tr2\t1000\t0\t400\t380\t400\n"         // unknown read -> drop
        "\n"                                                          // blank line -> ignored
        "garbage line that does not parse\n"                        // malformed -> ignored
        "r2\t1000\t100\t800\t+\tr3\t1000\t100\t800\t650\t700\n"     // (2,3) block 700
        "r3\t1000\t200\t500\t+\tr4\t1000\t200\t500\t280\t300";      // (3,4) block 300, no trailing NL

    // Expected canonical result, computed independent of chunking:
    //   (0,1): longest block 800  -> from first line, q*=r0[100,900], t*=r1[50,850]
    //   (0,2): longest block 560  -> from dup line "r2..r0": swapped q*=r0[40,600], t*=r2[30,520]
    //   (2,3): block 700          -> q*=r2[100,800], t*=r3[100,800]
    //   (3,4): block 300          -> q*=r3[200,500], t*=r4[200,500]
    const size_t lines = countLines(buf);

    vector<PafEntry> reference;
    for (size_t n = 1; n <= lines + 2; ++n) {
        auto got = importBuffer(buf, names, n);
        if (n == 1) {
            reference = got;
            REQUIRE(reference.size() == 4);
        }
        REQUIRE(got.size() == reference.size());
        for (size_t i = 0; i < got.size(); ++i) {
            REQUIRE(got[i].key == reference[i].key);
            REQUIRE(got[i].iv.blockLen == reference[i].iv.blockLen);
            REQUIRE(got[i].iv.qStart == reference[i].iv.qStart);
            REQUIRE(got[i].iv.qEnd == reference[i].iv.qEnd);
            REQUIRE(got[i].iv.tStart == reference[i].iv.tStart);
            REQUIRE(got[i].iv.tEnd == reference[i].iv.tEnd);
            REQUIRE(got[i].iv.isSameStrand == reference[i].iv.isSameStrand);
        }
    }

    // Spot-check the actual values on the reference (single-chunk) result.
    auto find = [&](ReadId a, ReadId b) -> const PafEntry* {
        uint64_t key = (uint64_t(a) << 32) | b;
        for (const auto& e : reference) if (e.key == key) return &e;
        return nullptr;
    };
    const PafEntry* e01 = find(0, 1);
    REQUIRE(e01 != nullptr);
    REQUIRE(e01->iv.blockLen == 800);
    REQUIRE(e01->iv.qStart == 100); REQUIRE(e01->iv.qEnd == 900);
    REQUIRE(e01->iv.tStart == 50);  REQUIRE(e01->iv.tEnd == 850);
    REQUIRE(e01->iv.isSameStrand == true);

    const PafEntry* e02 = find(0, 2);
    REQUIRE(e02 != nullptr);
    REQUIRE(e02->iv.blockLen == 560);
    REQUIRE(e02->iv.qStart == 40);  REQUIRE(e02->iv.qEnd == 600);   // r0 coords from dup line
    REQUIRE(e02->iv.tStart == 30);  REQUIRE(e02->iv.tEnd == 520);   // r2 coords
    REQUIRE(e02->iv.isSameStrand == false);

    const PafEntry* e23 = find(2, 3);
    REQUIRE(e23 != nullptr);
    REQUIRE(e23->iv.blockLen == 700);

    const PafEntry* e34 = find(3, 4);
    REQUIRE(e34 != nullptr);
    REQUIRE(e34->iv.blockLen == 300);   // the 150 line was dropped (<200)
}

TEST_CASE("import: CRLF line endings across chunk boundaries", "[paf]") {
    NameTable names;
    names.ids = {{"r0", 0}, {"r1", 1}, {"r2", 2}};
    const string buf =
        "r0\t1000\t100\t900\t+\tr1\t1000\t50\t850\t700\t800\r\n"
        "r1\t1000\t0\t400\t-\tr2\t1000\t0\t400\t380\t400\r\n";
    const size_t lines = countLines(buf);
    vector<PafEntry> reference = importBuffer(buf, names, 1);
    REQUIRE(reference.size() == 2);
    for (size_t n = 1; n <= lines + 2; ++n) {
        auto got = importBuffer(buf, names, n);
        REQUIRE(got.size() == reference.size());
        for (size_t i = 0; i < got.size(); ++i) {
            REQUIRE(got[i].key == reference[i].key);
            REQUIRE(got[i].iv.blockLen == reference[i].iv.blockLen);
        }
    }
}
