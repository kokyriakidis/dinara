# OverlapCigarStore: Per-Overlap Packed CIGAR Storage

## 1. Motivation

During overlap computation, dinara aligns pairs of reads and produces base-level
alignments between them. Previously, these alignments were immediately decomposed
into sparse SNP and indel evidence streams (`AlignedEvidenceStore`), discarding
the full alignment topology. While memory-efficient, this approach loses the
ability to:

1. **Map coordinates** between reads at arbitrary positions within an overlap.
2. **Extract subsequences** from overlapping reads around a region of interest,
   which is required for multiple sequence alignment (MSA).
3. **Walk the alignment** in a sliding-window fashion for haplotype phasing,
   where the phasing algorithm needs to identify mismatch positions, map them
   to the other read, and extract the corresponding base.

The `OverlapCigarStore` retains the full per-overlap CIGAR using a compact
encoding inspired by hifiasm's `asg16_v` data structure. It stores the complete
alignment topology at roughly half the memory cost of the existing evidence
store, while providing the query primitives needed for phasing and MSA.


## 2. Token Encoding

Each CIGAR operation is encoded as a single `uint16_t` (2 bytes):

```
  Bits 15-14:  Operation code (2 bits)
  Bits 13-0:   Run length    (14 bits, max 16383)
```

| Op code | Name       | Read0 (query) | Read1 (target) | Description                        |
|---------|------------|---------------|----------------|------------------------------------|
| 0       | Match      | consumes      | consumes       | Bases are identical                |
| 1       | Mismatch   | consumes      | consumes       | Substitution (both advance by 1)   |
| 2       | Insertion  | —             | consumes       | Extra bases in read1, not in read0 |
| 3       | Deletion   | consumes      | —              | Extra bases in read0, not in read1 |

Runs exceeding 16383 bases are automatically split across consecutive tokens
with the same operation code. All read operations coalesce consecutive
same-op tokens transparently, so callers never observe the split.

This encoding is identical to hifiasm's `push_trace` / `pop_trace` format
(`Levenshtein_distance.h`), where `c << 14 | len` packs the operation and
length into a single `uint16_t`.


## 3. Data Structure

### 3.1 Arena and Index

All tokens for all alignments are stored in a single flat array (the *arena*).
A separate index maps each alignment to its slice of the arena:

```
Arena:   [tok₀ tok₁ tok₂ ... tokₙ]
          ↑         ↑              ↑
Index:   [offset=0, count=3]  [offset=3, count=k]  ...
          alignment 0          alignment 1
```

```cpp
struct IndexEntry {
    uint64_t offset;  // Start position in the arena
    uint32_t count;   // Number of tokens
};

vector<CigarToken> arena;     // Flat token storage
vector<IndexEntry> index;     // Per-alignment index
```

This is the same architecture as hifiasm's `window_list_alloc`, where
`kvec_t(uint16_t) c` is the shared arena and each `window_list` entry has
`cidx` (offset) and `clen` (count).

### 3.2 Memory Layout

For a dataset of *N* alignments with *T* total tokens:

| Component | Size | Formula |
|-----------|------|---------|
| Arena     | 2*T* bytes | *T* × `sizeof(uint16_t)` |
| Index     | 16*N* bytes | *N* × `sizeof(IndexEntry)` |
| **Total** | **2*T* + 16*N*** | |

### 3.3 Empirical Measurements

Using the GIAB HG002 test dataset (1057 reads, chr1:15–15.4 Mb):

| Metric | Value |
|--------|-------|
| Alignments stored | 33,550 |
| Total tokens | 16,279,797 |
| Average tokens per alignment | 485 |
| Arena memory | 31.1 MB |
| Index memory | 0.5 MB |
| **Total memory** | **31.6 MB** |

For comparison, the `AlignedEvidenceStore` for the same dataset uses 60.9 MB,
making the CIGAR store **0.52× the memory** while retaining strictly more
information (full alignment topology vs. sparse diffs only).

| Store | Memory | Information retained |
|-------|--------|---------------------|
| AlignedEvidenceStore | 60.9 MB | Sparse SNPs + indels (two projections) |
| OverlapCigarStore | 31.6 MB | Full alignment topology (single bidirectional trace) |

The evidence store's memory is dominated by indel streams (50.3 MB), which
store two separate projections (query-view and target-view) with absolute
32-bit positions per indel. The CIGAR store encodes the same information as
2-byte run-length tokens in a single bidirectional trace.


## 4. CIGAR Construction

### 4.1 Projected Alignment Segments

Dinara's overlap alignments are computed in marker space: pairs of aligned
k-mer markers define anchor points, and the regions between consecutive
marker midpoints are aligned at the base level using A\*PA2 (an optimal
edit-distance aligner).

For an overlap with *N* aligned markers, there are *N*−1 segments. Each
segment spans from the midpoint of marker *i* to the midpoint of marker
*i*+1 on both reads:

```
Read0:  ──[midₐ]════════[mid_b]════════[mid_c]──
Read1:  ──[midₐ]════════[mid_b]════════[mid_c]──
           ↑     segment 0     ↑     segment 1    ↑
```

Segments tile contiguously at marker midpoints — no gaps, no overlaps.
Each marker's k-mer bases are split across two adjacent segments at the
midpoint. The CIGAR covers the span from the first marker midpoint to the
last marker midpoint on both reads.

### 4.2 Per-Segment CIGAR Generation

For each segment, one of the following cases applies:

| Case | Condition | CIGAR emitted |
|------|-----------|---------------|
| 1 | `len0 = 0, len1 = 0` | Nothing (empty segment) |
| 2 | `len0 > 0, len1 = 0` | Deletion of `len0` bases |
| 3 | `len0 = 0, len1 > 0` | Insertion of `len1` bases |
| 4 | Sequences identical | Match run of `len0` bases |
| 5 | Same length, differ | A\*PA2 CIGAR, split into match/mismatch runs |
| 6 | Different lengths | A\*PA2 CIGAR with match/mismatch/indel ops |

Case 4 is an optimization: when the byte-level comparison
`asciiSequence0 == asciiSequence1` succeeds, the segment is a perfect match
and no alignment is needed. The `vector::operator==` checks both size and
content, so this only triggers when every base matches.

For cases 5 and 6, the A\*PA2 CIGAR string (e.g., `50M2I48M`) is parsed.
`M`/`=`/`X` blocks are further split into contiguous match (op 0) and
mismatch (op 1) runs by comparing the underlying bases at each position.
This per-base distinction matches hifiasm's encoding and is necessary for
phasing, where mismatch positions are the candidate SNP sites.

### 4.3 Segment Stitching

The per-segment CIGARs are concatenated directly into the arena. Since
segments tile contiguously, the concatenation produces a single CIGAR
covering the full overlap. Consecutive tokens with the same operation code
(e.g., a match at the end of segment *i* followed by a match at the start
of segment *i*+1) are coalesced transparently during read operations.

### 4.4 Runtime Verification

After constructing each CIGAR, a verification pass walks the complete token
sequence and asserts that:

- Total bases consumed on read0 (match + mismatch + deletion) equals `totalLength[0]`
- Total bases consumed on read1 (match + mismatch + insertion) equals `totalLength[1]`

This assertion runs for every alignment during `computeAlignmentsWithEvidence`
and catches any inconsistency between the CIGAR and the projected alignment.


## 5. Thread-Parallel Construction

### 5.1 Thread-Local Stores

Each thread maintains its own `OverlapCigarStore` instance, avoiding
synchronization during the alignment phase. The store is passed to
`ProjectedAlignment` via a pointer parameter; when non-null,
`constructQuickRawSparse()` packs tokens into it alongside the existing
sparse evidence extraction. No additional alignment work is performed —
the CIGAR packing piggybacks on the existing A\*PA2 CIGAR parse loop.

### 5.2 Global Merge

After all threads complete, the thread-local stores are merged into the
global `Assembler::overlapCigarStore`:

1. The global store pre-allocates arena and index capacity based on the
   sum of all thread-local sizes.
2. Each thread's arena is appended to the global arena.
3. Each thread's index entries are appended with their offsets adjusted
   by the current global arena size.
4. The `cigarId` field in each `AlignmentInfo` is adjusted by the
   cumulative global index count.

This produces a single global store where `alignmentData[i].info.cigarId`
directly indexes into `overlapCigarStore`.

### 5.3 Performance Overhead

Measured on the GIAB HG002 test dataset (4 threads):

| Configuration | Mean time | Overhead |
|---------------|-----------|----------|
| Without CIGAR store | 3.029 s | — |
| With CIGAR store | 3.113 s | +84 ms (+2.8%) |

The overhead is negligible because the CIGAR packing consists of integer
shifts and vector appends, while the A\*PA2 alignment dominates the cost.


## 6. Query API

### 6.1 Basic Iteration

**`forEachOp(cigarId, callback)`** — Walks the CIGAR, coalescing consecutive
same-op tokens. The callback receives `(op, length)` for each coalesced
operation.

```cpp
store.forEachOp(cigarId, [](uint8_t op, uint32_t len) {
    // op: 0=match, 1=mismatch, 2=insertion, 3=deletion
});
```

### 6.2 Position-Tracking Iteration

**`forEachOpWithPositions(cigarId, read0Start, read1Start, callback)`** —
Like `forEachOp`, but maintains read0 (query) and read1 (target) positions
as the CIGAR is walked. The callback receives `(op, length, read0Pos, read1Pos)`
where the positions are the starting coordinates of each operation.

This corresponds to hifiasm's `xk`/`yk` tracking in `extract_sub_cigar_hc`.

```cpp
store.forEachOpWithPositions(cigarId, r0Start, r1Start,
    [](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        // r0, r1 are the positions at the start of this op
    });
```

### 6.3 Sub-Range Query

**`walkRange(cigarId, read0Start, read1Start, queryStart, queryEnd, callback)`** —
Walks only the sub-range `[queryStart, queryEnd)` on read0. Operations are
clipped at the range boundaries. Insertions anchored within the range are
reported with their full length.

This is the equivalent of hifiasm's `extract_sub_cigar_hc` seeking to
position `s` and walking until position `e`.

```cpp
store.walkRange(cigarId, r0Start, r1Start, windowStart, windowEnd,
    [](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        // Only ops overlapping [windowStart, windowEnd) are reported.
        // len is clipped to the range boundaries.
    });
```

### 6.4 Coordinate Mapping

**`queryToTarget(cigarId, read0Start, read1Start, queryPos)`** — Maps a
single read0 position to the corresponding read1 position by walking the
CIGAR. Returns `uint64_t(-1)` if the position falls inside a deletion
(no corresponding target base).

**`targetToQuery(cigarId, read0Start, read1Start, targetPos)`** — The
inverse mapping. Returns `uint64_t(-1)` if the position falls inside an
insertion.

These correspond to hifiasm's inline `t - xk + yk` arithmetic used in
`extract_sub_cigar_hc` to map query mismatch positions to target positions
for base extraction.

### 6.5 Resumable Cursor

**`Cursor`** — A struct that saves `(tokenIndex, xk, yk)` between successive
`walkRangeWithCursor` calls. This avoids re-scanning the CIGAR from the
beginning for each window in a sliding-window traversal.

**`walkRangeWithCursor(cursor, queryStart, queryEnd, callback)`** — Like
`walkRange`, but resumes from the cursor's saved position. If the new range
starts before the cursor, it seeks backward (like hifiasm's
`while (ck > 0 && xk > s)` loop). For the common forward-sliding-window
case, no backward seek is needed.

This matches hifiasm's `ovlp_cur_xoff`/`ovlp_cur_yoff`/`ovlp_cur_coff`
pattern used in `rphase_hc`, where cursor state is saved per overlap and
reused across successive phasing windows.

```cpp
OverlapCigarStore::Cursor cursor;
cursor.reset(cigarId, r0Start, r1Start, store);

for (uint64_t s = 0; s < readLength; s += windowStep) {
    uint64_t e = std::min(s + windowSize, readLength);
    store.walkRangeWithCursor(cursor, s, e,
        [](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
            // Process ops in [s, e)
        });
}
```


## 7. Comparison with hifiasm

### 7.1 Encoding

| Aspect | hifiasm | dinara |
|--------|---------|-------|
| Token type | `uint16_t` | `uint16_t` (`CigarToken`) |
| Op field | bits 15–14 (2 bits) | bits 15–14 (2 bits) |
| Length field | bits 13–0 (14 bits) | bits 13–0 (14 bits) |
| Max run per token | 16383 | 16383 |
| Long run handling | Split across tokens (`push_trace`) | Split across tokens (`pushOp`) |
| Coalescing on read | `pop_trace` | `forEachOp` / all walkers |

The encoding is identical.

### 7.2 Storage Architecture

| Aspect | hifiasm | dinara |
|--------|---------|-------|
| Arena | `kvec_t(uint16_t) c` in `window_list_alloc` | `vector<CigarToken> arena` |
| Per-alignment index | `cidx`/`clen` in `window_list` | `IndexEntry` (offset + count) |
| Scope | Per-overlap (one arena per `overlap_region`) | Global (one arena for all overlaps) |
| Thread safety | Per-thread overlap allocation | Per-thread store, merged post-alignment |

The main architectural difference is scope: hifiasm stores CIGARs per-overlap
with a per-window sub-index, while dinara uses a single global arena with a
flat per-alignment index. This is because dinara's segments (marker-midpoint
to marker-midpoint) are concatenated into a single CIGAR per overlap, whereas
hifiasm maintains separate CIGARs per window within each overlap.

### 7.3 Query API Parity

| hifiasm operation | hifiasm function | dinara API |
|-------------------|-----------------|------------|
| Walk CIGAR in `[s, e)` | `extract_sub_cigar_hc` | `walkRange` / `walkRangeWithCursor` |
| Track xk/yk positions | Inline in `extract_sub_cigar_hc` | `forEachOpWithPositions` |
| Map query → target | `t - xk + yk` (inline) | `queryToTarget` |
| Map target → query | — | `targetToQuery` |
| Resumable cursor | `ovlp_cur_xoff/yoff/coff` | `Cursor` struct |
| Backward seek | `while (ck > 0 && xk > s)` | `walkRangeWithCursor` backward seek |
| Identify mismatches | `if(op == 1) f[t-s]++` | Callback checks `op == 1` |
| Identify indels | `if(op != 0)` in `iter_sub_cigar_sv` | Callback checks `op == 2` or `op == 3` |
| Count errors in range | `extract_sub_err` | `walkRange` with accumulator |
| Base extraction | `qstr[t]`, `tu->seq[yk]` | Direct sequence access (not store's job) |

dinara provides `targetToQuery` which hifiasm does not have as a standalone
operation (hifiasm only maps query→target).

### 7.4 Optimizations

| Optimization | hifiasm | dinara |
|--------------|---------|-------|
| 2-byte tokens | ✅ | ✅ |
| Flat shared arena | ✅ | ✅ |
| Run-length splitting at 16383 | ✅ | ✅ |
| Coalescing on read | ✅ | ✅ |
| Resumable cursor | ✅ | ✅ |
| Backward seek from cursor | ✅ | ✅ |
| Early termination at range end | ✅ | ✅ |
| Pre-allocated arena (known size) | ✅ | ✅ |
| Base-packed tokens | ✅ (error correction only) | Not implemented (not needed for phasing) |

The base-packed token variant (`push_trace_bp_f`) embeds the actual
mismatched bases into the token. hifiasm uses this during error correction
but not during phasing — the phasing code reads bases from the sequence
store directly. dinara follows the same approach.


## 8. Phasing Workflow

The `OverlapCigarStore` provides all primitives needed to implement
hifiasm-style haplotype phasing. The phasing workflow proceeds as follows:

### 8.1 Pass 1: SNP Site Discovery

For each query read, slide a window `[s, e)` across its length. For each
window, iterate over all overlapping alignments and walk their CIGARs in
the window range:

```cpp
Cursor cursor;
cursor.reset(cigarId, r0Start, r1Start, store);

// For each window [s, e):
store.walkRangeWithCursor(cursor, s, e,
    [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        if (op == 1) {  // mismatch
            for (uint32_t i = 0; i < len; i++) {
                mismatchCount[r0 + i - s]++;
            }
        }
    });
```

Positions where the mismatch count exceeds a threshold are candidate SNP
sites.

### 8.2 Pass 2: Allele Evidence Collection

Walk the same CIGARs again. At candidate SNP sites:

- **Match (op 0)**: The overlapping read agrees with the query. Record the
  query base as evidence for the reference allele.
- **Mismatch (op 1)**: The overlapping read disagrees. Use `queryToTarget`
  (or the callback's `r1Pos`) to find the corresponding position on the
  target read, extract the target base from the sequence store, and record
  it as evidence for the alternate allele.

### 8.3 Indel Detection

Walk the CIGAR and collect all non-match operations:

```cpp
store.forEachOpWithPositions(cigarId, r0Start, r1Start,
    [&](uint8_t op, uint32_t len, uint64_t r0, uint64_t r1) {
        if (op == 2) {
            // Insertion: len bases in read1 at position r1,
            // anchored at read0 position r0.
            recordIndel(r0, r1, len, INSERTION);
        } else if (op == 3) {
            // Deletion: len bases in read0 at position r0,
            // no corresponding bases in read1 (r1 is the
            // position where the deletion is anchored).
            recordIndel(r0, r1, len, DELETION);
        }
    });
```

### 8.4 Region Extraction for MSA

To extract the subsequence of an overlapping read corresponding to a query
region `[qs, qe)`:

1. Map the query boundaries to target coordinates:
   ```cpp
   uint64_t ts = store.queryToTarget(cigarId, r0Start, r1Start, qs);
   uint64_t te = store.queryToTarget(cigarId, r0Start, r1Start, qe - 1) + 1;
   ```
2. Extract `targetRead[ts..te)` from the sequence store.
3. Feed the extracted sequences into an MSA algorithm.


## 9. Files

| File | Description |
|------|-------------|
| `src/OverlapCigarStore.hpp` | Store class, token encoding, cursor, all query methods |
| `src/Alignment.hpp` | `AlignmentInfo::cigarId` field linking alignments to CIGARs |
| `src/ProjectedAlignment.hpp` | `OverlapCigarStore*` parameter and `cigarId` member |
| `src/ProjectedAlignment.cpp` | CIGAR packing in `constructQuickRawSparse()` |
| `src/Assembler.hpp` | Global `overlapCigarStore` and thread-local stores |
| `src/AssemblerComputeAlignments.cpp` | Thread init, store passing, merge, logging |
| `tests/test_overlap_cigar_store.cpp` | 43 unit tests (861 assertions) |


## 10. Testing

The implementation is verified at two levels:

### 10.1 Unit Tests (43 test cases, 861 assertions)

| Category | Tests | What is verified |
|----------|-------|-----------------|
| Token encoding | 1 | Round-trip op/len through `CigarToken` |
| Segment cases 1–6 | 6 | All segment types produce correct tokens and consumed-base counts |
| Multi-segment stitching | 2 | Mixed identical + aligned segments, deletion between matches |
| Long run splitting | 2 | Runs of 20K and 50K bases split and coalesce correctly |
| Thread merge | 2 | Offset adjustment across single and multi-alignment merges |
| Coalescing | 3 | Same-op coalescing, no cross-op coalescing |
| Edge cases | 4 | Single-base mismatch, alternating ops, empty alignment, invalid ID |
| Realistic overlap | 1 | 10 segments with mixed cases, total consumed bases verified |
| Position tracking | 3 | `forEachOpWithPositions` with match, insertion, deletion |
| Coordinate mapping | 7 | `queryToTarget`, `targetToQuery`, inverse property |
| Sub-range query | 5 | `walkRange` clipping, insertions, deletions, sub-op ranges |
| Resumable cursor | 7 | Equivalence with `walkRange`, forward sliding, backward seek, insertions, deletions, phasing-like pattern |

### 10.2 Runtime Assertions

Every alignment produced during `computeAlignmentsWithEvidence` is verified
by walking its CIGAR and asserting that consumed bases match `totalLength`
on both reads. This runs on all 33,550 alignments in the test dataset.
