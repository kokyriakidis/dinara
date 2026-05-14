# Overlap Candidate Finding in Dinara

This document describes how Dinara discovers overlap candidates between reads.

## 1. Overview

The candidate finding pipeline identifies pairs of reads that likely share a
true genomic overlap. It produces scored overlap candidates that are passed to
downstream DP chaining and base-level alignment. The pipeline has six stages:

1. **Marker extraction** — select informative k-mer positions from each read
   using closed syncmers in run-length encoded (RLE) space.
2. **Frequency filtering** — remove k-mers whose genome-wide frequency falls
   outside `[minFreq, maxFreq]`, eliminating sequencing errors (singletons)
   and extremely repetitive k-mers (count > 2000).
3. **Inverted index construction** — build a count-then-scatter lookup table
   mapping each canonical k-mer to all (readId, position) pairs where it occurs.
4. **Hit collection** — for each query read, look up every marker in the index,
   apply frequency-based weighting, and downsample consecutive high-frequency
   marker streaks.
5. **Hit grouping** — radix-sort hits by (partnerReadId, posA), then dispatch
   each per-partner group to the DP chaining pipeline.
6. **Post-filtering** — apply per-overlap-type caps (max_n_chain), coverage
   window rescue (COV_W), and weak-chain suppression (R485).

The design follows hifiasm's overlap discovery approach (ecovlp.cpp, anchor.cpp,
Hash_Table.cpp, sketch.cpp) but replaces hifiasm's minimizer-based position
table with a closed-syncmer inverted index built using a count-then-scatter
algorithm that reduces peak memory during construction.

### Pipeline Diagram

```
  Reads (FASTQ)
       │
       ▼
  ┌─────────────────────────────┐
  │  RLE compression            │  Homopolymer runs → single bases
  │  Closed syncmer extraction  │  k=50, s=11, SIMD-accelerated
  └─────────────┬───────────────┘
                │
                ▼
  ┌─────────────────────────────┐
  │  Frequency filtering        │  Keep k-mers with count ∈ [2, 2000]
  │  Coverage peak detection    │  Find modal k-mer frequency
  └─────────────┬───────────────┘
                │
                ▼
  ┌─────────────────────────────┐
  │  Count-then-scatter index   │  Phase 1: parallel counting
  │  construction               │  Phase 2: merge + prefix sum
  │                             │  Phase 3: parallel scatter
  │                             │  Phase 4: query hash table
  └─────────────┬───────────────┘
                │
                ▼
  ┌─────────────────────────────┐
  │  Per-read hit collection    │  Hash table lookup per marker
  │  + frequency weighting      │  3-tier weight: rare/normal/repetitive
  │  + high-freq downsampling   │  Streak-based, keep ≤ 16 per streak
  └─────────────┬───────────────┘
                │
                ▼
  ┌─────────────────────────────┐
  │  Radix sort by              │  LSD radix sort, 8 passes
  │  (partnerReadId, posA)      │  O(N) time, O(N) space
  └─────────────┬───────────────┘
                │
                ▼
  ┌─────────────────────────────┐
  │  Per-partner DP chaining    │  → see chaining.md
  │  + mcopy extraction         │
  │  + post-filtering           │
  └─────────────┬───────────────┘
                │
                ▼
  Overlap candidates (scored, filtered)
```

## 2. Marker Extraction

Dinara operates in **run-length encoded (RLE) space**, where homopolymer runs
are collapsed to single bases before any k-mer analysis. This is equivalent to
hifiasm's homopolymer compression (HPC) and reduces noise from homopolymer
length errors, which are the dominant error mode in Oxford Nanopore sequencing.

### 2.1 Closed Syncmers

When `Kmers.useSimdClosedSyncmers = true` (the recommended and default
setting), Dinara extracts markers using **closed syncmers** with parameters:

- **k = 50** — k-mer length in RLE bases
- **s = 11** — sub-kmer (s-mer) size for syncmer selection

**Definition.** A k-mer is a *closed (s, k)-syncmer* if the minimum-hash s-mer
within the k-mer occurs at the first or last position. Formally, given a k-mer
`K = K[0..k-1]` and a hash function `h`, let `p = argmin_{i ∈ [0, k-s]}
h(K[i..i+s-1])`. Then `K` is a closed syncmer if and only if `p = 0` or
`p = k - s`.

**Properties:**

- **Deterministic**: The same k-mer always produces the same selection decision,
  regardless of context.
- **Strand-symmetric**: A k-mer and its reverse complement are both selected or
  both rejected, because the minimum s-mer position is preserved under reversal.
- **Approximately uniform**: Closed syncmers sample roughly `2/(k-s+1)` of all
  k-mer positions. With k=50 and s=11, the sampling rate is approximately
  `2/40 = 5%` of positions, comparable to hifiasm's minimizer density with
  k=51, w=51.
- **Context-free**: Unlike minimizers, syncmer selection depends only on the
  k-mer itself, not on neighboring k-mers. This means insertions or deletions
  in one part of a read do not affect syncmer selection in distant parts.

The extraction is SIMD-accelerated (`findMarkersSimdClosedSyncmers` in
`AssemblerMarkers.cpp`) and processes reads in parallel across threads using a
two-pass approach: pass 1 counts markers per read to allocate storage, pass 2
extracts and stores them.

### 2.2 Canonicalization

Each extracted k-mer is stored as a `KmerId` — a 128-bit integer encoding the
two-bit-per-base representation (A=00, C=01, G=10, T=11). The k-mer is split
into two k-bit halves (LSB and MSB) packed into the 128-bit value.

For each k-mer, we compute its reverse complement using bit-reversal
(`getRcKmerId` in `InvertedIndexBuilder.hpp`):

```
RC(K) = complement(reverse(K))
```

The **canonical form** is the lexicographically smaller of the k-mer and its
reverse complement:

```
canonical(K) = min(K, RC(K))
```

This ensures that overlapping reads match regardless of which strand they
originate from. A boolean flag `isRc` records whether the observed k-mer was
the reverse complement of the canonical form.

### 2.3 Canonical Cache

To avoid redundant reverse-complement computation during hit collection, the
builder pre-computes and caches canonical k-mer IDs and RC flags for all
strand-0 markers in flat arrays:

| Array | Type | Contents |
|-------|------|----------|
| `strand0CanonicalKmerIds` | `KmerId[]` | Canonical k-mer ID for each marker |
| `strand0CanonicalIsRc` | `uint8_t[]` | 1 if observed k-mer is RC of canonical |
| `strand0CanonicalOffsets` | `uint64_t[]` | Per-read offset into the above arrays |

For read `r`, its markers occupy positions
`[strand0CanonicalOffsets[r], strand0CanonicalOffsets[r+1])` in the flat arrays.
This cache is populated in parallel during Phase 1 of index construction.

### 2.4 Frequency Filtering

Before index construction, `applyKmerCountFilter` removes markers whose
genome-wide frequency falls outside `[minFreq, maxFreq]`:

- **minFreq = 2**: Removes singleton k-mers, which are almost certainly
  sequencing errors. A k-mer that appears only once in the entire dataset
  cannot produce any useful hit.
- **maxFreq = invertedIndexMaxKmerCount = 2000**: Removes extremely
  high-frequency k-mers (centromeric repeats, telomeric sequences, satellite
  DNA). This matches hifiasm's `max_kmer_cnt = 2000` hard cutoff
  (`htab.cpp:gen_hh`).

**Relationship to hifiasm's three-tier k-mer handling:**

Hifiasm uses three data structures with different frequency thresholds:

1. **Position index** (`ha_pt_gen`): Keeps all k-mers with count ∈ [2, ∞)
   when the filter table exists. The position index itself does not apply an
   upper frequency bound.
2. **Filter table** (`ha_flt_tab`, built by `gen_hh`): Records k-mers with
   count ≥ `5 × hom_cov`. K-mers with count > `max_kmer_cnt = 2000` are
   marked with `INT16_MAX` (effectively excluded).
3. **Sketching** (`minimizers_qgen0`): During marker extraction, k-mers with
   `ha_ft_cnt ≥ 1<<28` (i.e., count > 2000) are excluded entirely. K-mers
   with count ∈ `[5 × hom_cov, 2000]` are kept but annotated with their
   count in the `rid` field for downstream downsampling.

Dinara replicates this behavior by:
- Hard-filtering at `maxFreq = 2000` (equivalent to hifiasm's `max_kmer_cnt`)
- Keeping k-mers with count ∈ `[5 × coveragePeak, 2000]` in the marker set
  and handling them via streak-based downsampling during hit collection
  (Section 5)

### 2.5 Coverage Peak Detection

The k-mer frequency distribution is analyzed to find the **coverage peak**
(`coveragePeak`), which represents the expected k-mer frequency for
single-copy genomic sequence. For a genome sequenced at 30× coverage, the
peak is typically around 30.

The coverage peak is used downstream to set coverage-relative thresholds:
- `lowFreqThreshold = coveragePeak × 0.333` (rare k-mers, weight = 2)
- `highFreqThreshold = coveragePeak × 1.667` (repetitive k-mers, weight > 1)
- `max_n_chain = max(coveragePeak × 5, 100)` (per-type overlap cap)

## 3. Inverted Index Construction

The inverted index maps each canonical k-mer to the list of (readId, position)
pairs where it occurs. Only strand-0 markers are indexed because canonical
k-mers are strand-symmetric — reverse-complement matches are detected during
hit collection by comparing the RC flags of the query and target observations.

### 3.1 Count-Then-Scatter Algorithm

The index is built using a **count-then-scatter** approach
(`InvertedIndexBuilder.hpp`) that avoids materializing a large intermediate
array and sorting it. The algorithm has six phases:

#### Phase 1: Parallel Canonicalization

Each thread processes a range of reads and computes the canonical k-mer ID and
RC flag for every marker. Results are written to the pre-allocated flat arrays
(`strand0CanonicalKmerIds`, `strand0CanonicalIsRc`) using per-read offsets
computed from marker counts. No synchronization is needed because each thread
writes to a disjoint range.

#### Phase 2: Parallel Counting

Each thread processes a range of reads and builds a **thread-local hash table**
mapping canonical k-mer → occurrence count. No positions are stored in this
phase — only the count is incremented.

The thread-local tables use open-addressing with linear probing, sized at 4×
the expected distinct k-mer count to keep the load factor below 0.25 and
collision rates low:

```
tableSize = smallest power of 2 ≥ 4 × threadMarkerCount
```

The `findOrInsert` function probes linearly from `hashKmer(key) & mask`:

```cpp
uint64_t slot = hashKmer(key) & mask;
while (true) {
    if (table[slot].empty)  → insert new entry, return
    if (table[slot].key == key) → return existing entry
    slot = (slot + 1) & mask;   → linear probe
}
```

#### Phase 3: Merge

Thread-local tables are merged into a single global counting table. This
touches only distinct k-mers (not all occurrences), so it is fast even
single-threaded. Each thread-local table is freed immediately after merging
to reduce peak memory.

#### Phase 4: Prefix Sum

A sequential scan over the global table assigns each k-mer a contiguous range
in the final compact array. The starting offset for k-mer *i* is the
cumulative sum of counts for all k-mers before it:

```
globalTable[i].start = Σ_{j < i} globalTable[j].count
```

A parallel array of `ScatterState` structs is initialized with atomic write
cursors set to each k-mer's start offset.

#### Phase 5: Parallel Scatter

Each thread re-iterates its reads. For each marker, it looks up the k-mer in
the global table and atomically claims a slot in the pre-allocated compact
array using `fetch_add`:

```cpp
uint64_t writeIdx = scatterStates[slot].writeCursor.fetch_add(1);
compactOccurrences[writeIdx] = {readId, encodedPosition};
```

The compact array stores 8-byte entries:

```cpp
struct CompactOccurrence {
    ReadId   readId;    // 4 bytes — which read contains this k-mer
    uint32_t position;  // 4 bytes — base position on the read
                        //   bits [30:0] = position
                        //   bit  [31]   = isRc flag (1 if observed k-mer
                        //                 is RC of canonical form)
};
```

The top bit of `position` encodes the RC flag, allowing the hit collection
phase to determine strand orientation without a separate lookup. During hit
collection, the position is extracted as `position & 0x7FFFFFFF` and the RC
flag as `position >> 31`.

#### Phase 6: Query Hash Table

A separate open-addressing hash table is built for O(1) k-mer lookups during
hit collection. This table is sized at 2× the distinct k-mer count (load
factor ≤ 0.5) and uses linear probing:

```cpp
struct HashEntry {
    KmerId   key;    // 16 bytes — canonical k-mer
    uint64_t start;  // 8 bytes  — offset into compactOccurrences
    uint32_t count;  // 4 bytes  — number of occurrences
    bool     empty;  // 1 byte   — slot occupancy flag
};
```

The hash function (`hashKmer`) uses Boost-style hash combining to fold the
128-bit KmerId into a 64-bit hash:

```cpp
uint64_t hashKmer(KmerId k) {
    uint64_t k1 = lower_64_bits(k);
    uint64_t k2 = upper_64_bits(k);
    return k1 ^ (k2 + 0x9e3779b9 + (k1 << 6) + (k1 >> 2));
}
```

The constant `0x9e3779b9` is derived from the golden ratio and provides good
bit mixing.

### 3.2 Memory Analysis

| Metric | Old (sort-based) | New (count-scatter) |
|--------|-----------------|---------------------|
| Intermediate per occurrence | 24 bytes (KmerId + ReadId + pos) | 0 bytes |
| Sort buffer | 24 bytes (radix sort dst) | 0 bytes |
| Peak transient | ~48 bytes/occurrence | ~8 bytes/occurrence |
| Final compact array | 8 bytes/occurrence | 8 bytes/occurrence |
| Radix sort passes | ⌈128/8⌉ = 16 for 128-bit keys | 0 |
| Thread-local tables | 0 | ~32 bytes/distinct k-mer |

For a dataset with 100M marker occurrences, the old approach required ~4.5 GB
of transient memory for sorting alone. The new approach allocates only the
final 800 MB compact array plus modest hash table overhead.

The thread-local counting tables are the main transient cost. With T threads
and D distinct k-mers per thread, the overhead is approximately
`T × 4 × D × sizeof(CountEntry)`. Since `CountEntry` is ~30 bytes and each
thread sees roughly `D_total / T` distinct k-mers, the total is approximately
`4 × D_total × 30 ≈ 120 × D_total` bytes, which is typically much smaller
than the compact occurrence array.

## 4. Hit Collection

For each query read A, the hit collection phase finds all reads that share at
least one k-mer with A. This is the most compute-intensive part of candidate
finding.

### 4.1 Per-Marker Lookup

For each marker position `i` in read A:

1. **Retrieve canonical k-mer**: Either from the pre-computed cache
   (`canonicalIdsA[i]`, `canonicalIsRcA[i]`) or by computing it on the fly
   via `getRcKmerId` and taking the lexicographic minimum.

2. **Hash table probe**: Compute `slot = hashKmer(canonicalKId) & hashMask`
   and probe linearly until finding the matching key or an empty slot:

   ```cpp
   while (!hashTable[slot].empty) {
       if (hashTable[slot].key == canonicalKId) {
           startIdx = hashTable[slot].start;
           count = hashTable[slot].count;
           found = true;
           break;
       }
       slot = (slot + 1) & hashMask;
   }
   ```

3. **Decision branch**: Based on the k-mer's frequency and the downsampling
   configuration:
   - **Not found**: The k-mer was filtered out during frequency filtering.
     If downsampling is active, flush any pending high-frequency streak and
     update the boundary position.
   - **High-frequency** (`count > highFreqDownsampleThreshold`): Buffer the
     marker into the current high-frequency streak for later downsampling
     (Section 5). Do not emit hits yet.
   - **Normal/rare frequency**: Flush any pending streak, then emit hits
     for all occurrences of this k-mer.

4. **Hit emission**: For each occurrence in `compactOccurrences[startIdx ..
   startIdx + count - 1]`, skip self-hits (`readId == readIdA`) and emit a
   flat hit record:

   ```cpp
   struct InvertedIndexTempHit {
       ReadId   partnerReadId;  // partner read (Read B)
       uint32_t posA;           // base position in Read A
       uint32_t posB;           // base position in Read B (bits [30:0])
       uint32_t ordinalA;       // marker ordinal in Read A
       uint32_t weight;         // frequency-based weight
       uint8_t  isRcA;          // RC flag for Read A's k-mer
       uint8_t  isRcB;          // RC flag for Read B's k-mer
   };
   ```

### 4.2 Frequency-Based Weighting

Each hit is assigned a weight based on the genome-wide frequency of its k-mer.
The weight normalizes the contribution of each anchor during DP chaining,
following hifiasm's `normal_w` convention (Hash_Table.cpp:20).

The weight computation uses three tiers:

| Tier | Condition | Weight | Rationale |
|------|-----------|--------|-----------|
| Rare | `count ≤ lowFreqThreshold` | 2 | More informative — likely unique genomic position |
| Normal | `lowFreq < count < highFreq` | 1 | Baseline — typical single-copy k-mer |
| Repetitive | `count ≥ highFreqThreshold` | `⌊(1 + ⌈count / (2 × highFreq)⌉)^1.1⌋` | Downweight repeats |

Where:
- `lowFreqThreshold = max(2, ⌊coveragePeak × 0.333⌋)` — hifiasm's
  `HA_KMER_GOOD_RATIO = 0.333`
- `highFreqThreshold = max(1, ⌊coveragePeak × 1.667⌋)` — hifiasm's
  `2.0 - HA_KMER_GOOD_RATIO = 1.667`
- `highFreqWeightUnit = max(1, highFreqThreshold × 2)`

**Weight lookup table (LUT):** To avoid expensive `pow()` calls in the hot
path, a 512-entry lookup table is pre-computed at initialization:

```cpp
weightLut[i] = ⌊i^1.1⌋   for i ∈ [0, 511]
```

For normalized counts beyond 511, the weight falls back to direct `pow()`
computation. In practice, normalized counts rarely exceed 512.

**How weight is used in DP scoring:** During chaining, the weight is packed
into the upper 24 bits of the `cnt` field of `HifiasmKmerHit`:

```
cnt = (weight << 8) | span
```

The scoring function `comput_sc_ch_ec` divides the base score by the weight
via `normal_w(score, weight)`, so high-frequency k-mers contribute less per
anchor. This prevents repetitive regions from dominating chain scores.

### 4.3 Strand Determination

The strand of a hit is determined by comparing the RC flags of the k-mer
observation in reads A and B:

```
strand = isRcA ⊕ isRcB
```

- `strand = 0` (same-strand): Both reads observed the same canonical
  orientation (both forward or both RC). In a true overlap, target positions
  increase monotonically with query positions.
- `strand = 1` (opposite-strand): Reads observed opposite orientations. In a
  true overlap, target positions decrease as query positions increase.

This is equivalent to hifiasm's `z->rev == y->rev` check in anchor.cpp.

## 5. High-Frequency Marker Downsampling

Consecutive high-frequency markers (repetitive regions like tandem repeats,
segmental duplications, or satellite DNA) can produce an explosion of hits
that overwhelms memory and compute. To control this, Dinara implements
hifiasm's streak-based downsampling algorithm (sketch.cpp: `select_mz_h` +
`hf_select`).

### 5.1 Streak Detection

As markers are processed left-to-right along the query read, consecutive
markers whose frequency exceeds `highFreqDownsampleThreshold` are buffered
into a **streak**. The streak ends when either:

- A non-high-frequency marker is encountered, or
- A marker is not found in the index (filtered out), or
- The end of the read is reached.

Each buffered marker is stored as a `PendingHighFrequencyMarker`:

```cpp
struct PendingHighFrequencyMarker {
    uint64_t startIdx;   // offset into compactOccurrences
    uint32_t count;      // genome-wide occurrence count
    uint64_t hashKey;    // yak_hash64_64(fold(canonicalKmerId))
    uint32_t posA;       // base position on query read
    uint32_t ordinalA;   // marker ordinal on query read
    uint32_t weight;     // pre-computed frequency weight
    uint8_t  isRcA;      // RC flag
};
```

The `highFreqDownsampleThreshold` is set to
`max(3, ⌊coveragePeak × highFreqMultiplier⌋)`. The minimum of 3 prevents
treating nearly all k-mers as high-frequency in low-coverage datasets.

### 5.2 Streak Flushing Algorithm

When a streak ends, the `flushHighFrequencyStreak` function decides how many
markers to retain:

**Step 1 — Compute span and keep count:**

```
leftPos  = position of last non-high-freq marker (or 0 if none)
rightPos = position of next non-high-freq marker (or read length)
span     = rightPos - leftPos
keep     = round(span / sampleDistance)     // sampleDistance = 500
if keep > maxPerStreak: keep = maxPerStreak  // maxPerStreak = 16
```

**Step 2 — Handle edge cases:**

- If `keep == 0` (streak shorter than ~250 bp): **discard all markers** in the
  streak. This matches hifiasm's behavior where `select_mz_h` sets
  `max_high_occ = 0`, causing `hf_select` to never be called, so all markers
  in the streak retain their non-zero `rid` and are squeezed out.

- If `keep ≥ streak.size()`: **emit all markers** (no downsampling needed).

**Step 3 — Select top-k markers:**

Use `std::nth_element` (O(n) average) to select the `keep` markers with the
smallest `(count, hashKey)` composite key. This prefers:
1. Markers with lower occurrence count (more informative)
2. Among equal counts, markers with smaller hash (deterministic tie-breaking)

The hash key is computed as `yak_hash64_64(fold(canonicalKmerId))`, where:
- `fold` XOR-folds the 128-bit KmerId into 64 bits
- `yak_hash64_64` is hifiasm's Thomas Wang-style bijective integer hash

**Step 4 — Sort selected markers by position:**

The selected markers are sorted by `posA` before emission so that downstream
radix sort sees monotonically increasing query positions within each streak.

**Step 5 — Emit with span guard:**

Each selected marker is emitted through a callback that applies one final
filter: if `count ≥ span`, the marker is discarded. This matches hifiasm's
`rid < pe - ps` guard — a marker whose occurrence count exceeds the local
genomic span is too repetitive even relative to the local context.

### 5.3 Downsampling Example

Consider a 5 kb tandem repeat region with 100 high-frequency markers:

```
span = 5000 bp
keep = round(5000 / 500) = 10
maxPerStreak = 16 → keep = min(10, 16) = 10
```

The algorithm selects the 10 markers with the lowest occurrence counts,
reducing the hit explosion from `100 × avg_count` hits to `10 × avg_count`
hits while preserving positional signal across the repeat.

### 5.4 Shared Implementation

The `flushHighFrequencyStreak` function is a C++ template parameterized on
the emit callback:

```cpp
template<typename EmitFn>
static inline void flushHighFrequencyStreak(
    vector<PendingHighFrequencyMarker>& streak,
    vector<size_t>& workspace,
    int64_t lastNonHighPos,
    uint32_t rightBoundaryPos,
    uint32_t sampleDistance,
    uint32_t maxPerStreak,
    EmitFn&& emitFn);
```

Both the discovery path and the PAF path use this same function with different
emit callbacks:
- **Discovery path**: Emits hits for all partner reads (`readId != readIdA`)
- **PAF path**: Emits hits only for the specific partner read (`readId == readIdB`)

## 6. Hit Grouping and Per-Partner Dispatch

After collecting all hits for a query read, they must be grouped by partner
read and sorted for the DP chaining algorithm.

### 6.1 LSD Radix Sort

Hits are sorted by a packed 64-bit key:

```
key = (partnerReadId << 32) | posA
```

This groups hits by partner read (high 32 bits) and within each partner,
orders them by query position (low 32 bits).

The sort uses **LSD (Least Significant Digit) radix sort** with 8 passes,
one per byte of the 64-bit key, starting from the least significant byte.
Each pass is a stable counting sort:

```
for pass = 0 to 7:
    shift = pass × 8
    1. Histogram: count occurrences of each byte value
    2. Prefix sum: convert histogram to output offsets
    3. Scatter: place each hit at its sorted position in tmp
    4. Swap hits ↔ tmp
```

After 8 passes (even number), the sorted result is in the original `hits`
vector. The sort is O(N) time and O(N) space.

**Why radix sort over `std::sort`?** With millions of hits per read, the O(N)
complexity of radix sort provides a measurable speedup over O(N log N)
comparison-based sorts. Hifiasm uses radix sorts in the same hot path for the
same reason.

### 6.2 Per-Partner Iteration

After sorting, hits for the same partner read are contiguous. The code
iterates through the sorted array, identifying partner-read boundaries:

```cpp
size_t hitIter = 0;
while (hitIter < flatHits.size()) {
    ReadId readIdB = flatHits[hitIter].partnerReadId;

    // Skip mirrored pairs: only process readIdB > readIdA
    if (readIdB <= readIdA) { skip group; continue; }

    // Collect all hits for this partner
    size_t startInFlat = hitIter;
    while (hitIter < flatHits.size() &&
           flatHits[hitIter].partnerReadId == readIdB)
        ++hitIter;
    size_t numHits = hitIter - startInFlat;

    // Dispatch to DP chaining...
}
```

The `readIdB <= readIdA` check ensures each pair is processed only once
(the pair (A, B) with A < B is canonical).

### 6.3 Per-Strand Minimum Hit Count

Before running the DP, hits are counted per strand orientation:

```cpp
size_t revCount[2] = {0, 0};
for (each hit h) {
    uint8_t rev = h.isRcA ^ h.isRcB;
    ++revCount[rev];
}
bool keepRev0 = (revCount[0] >= chain_cutoff);  // chain_cutoff = 2
bool keepRev1 = (revCount[1] >= chain_cutoff);
```

If neither strand has at least `chain_cutoff = 2` hits, the entire partner
group is skipped. This matches hifiasm's `ecovlp.cpp:3274` where
`chain_cutoff = 2` is passed to `lchain_qgen_mcopy_fast`.

### 6.4 Coordinate Transform to Hifiasm Format

Hits that pass the strand filter are converted to `HifiasmKmerHit` anchors
for the DP chaining algorithm. The key transformation is converting Dinara's
start-position coordinates to hifiasm's end-position coordinates:

```cpp
// Dinara stores start positions; hifiasm uses end positions
uint32_t seedSpan = min(kmerLen, 255);
uint32_t selfOff  = posA + (seedSpan - 1);        // query end position
uint32_t offSame  = posB + (seedSpan - 1);         // target end (same strand)
uint32_t offDiff  = readLenB - 1 - posB;           // target end (diff strand)

HifiasmKmerHit kh;
kh.self_offset = selfOff;
kh.offset = (strand == 0) ? offSame : offDiff;
kh.cnt = (weight << 8) | min(seedSpan, 255);
kh.strand = strand;
```

For opposite-strand hits, the target position is reflected:
`offset = readLenB - 1 - posB`. This converts from forward-strand coordinates
to the reverse-complement coordinate system that hifiasm's DP expects.

### 6.5 Anchor Sorting

Anchors are split by strand and sorted by `(self_offset, offset)` within each
strand group. Since hits are already sorted by `posA` (from the radix sort),
`self_offset` is already non-decreasing. Only ties on `self_offset` need
sub-sorting by `offset`, which is done efficiently by
`sortHifiasmHitsBySelfOffsetThenOffsetRuns`:

```
1. Find runs of anchors with identical self_offset
2. Within each run, sort by offset
3. Move to next run
```

This is O(n) in the common case (no ties) and O(n log k) worst case where k
is the maximum run length.

The two strand groups are then concatenated (strand 0 first, strand 1 second)
into a single array for the DP, matching hifiasm's convention.

### 6.6 Marker Ordinal Resolution

The DP chaining produces chains of anchors with base positions. To construct
Dinara `Alignment` objects (which use marker ordinals), we need to map each
hit's `posB` to its marker ordinal in read B.

This is done by `mapHitPositionsToMarkerOrdinals` using a **tandem scan**:

1. Sort hit indices by `posB`
2. Walk markers of read B in position order
3. For each hit (in posB order), advance the marker pointer until it matches

This is O(n log n + m) where n = number of hits and m = number of markers in
read B, avoiding O(n log m) repeated binary searches.

## 7. Two Entry Points: Discovery and PAF

The candidate finding pipeline has two entry points that share the same
index, hit collection, weighting, and downsampling code but differ in how
read pairs are selected.

### 7.1 Discovery Path (`findAlignmentCandidatesInvertedIndex`)

The discovery path performs **all-vs-all** overlap detection:

1. `buildInvertedIndex()` — constructs the count-then-scatter index
2. `chainAlignmentCandidates()` — for each read A, queries every marker
   against the index, collects hits for all partner reads, and runs DP
   chaining per partner

This is the primary path used during assembly. Each query read produces hits
against every other read that shares at least one k-mer, and the DP chaining
+ post-filtering pipeline selects the best overlaps.

### 7.2 PAF Path (`chainPafCandidates`)

The PAF path receives **pre-determined read pairs** from an external PAF file
(e.g., produced by minimap2) and re-chains them using Dinara's hifiasm-
compatible DP scoring:

1. `buildInvertedIndex()` — same index construction
2. `importAlignmentCandidatesFromPaf()` — reads PAF records into candidate
   pairs (external to this file)
3. `chainPafCandidates()` — for each imported pair (A, B), collects only
   hits between A and B, runs DP chaining for the PAF-specified orientation

Key differences from the discovery path:

| Aspect | Discovery | PAF |
|--------|-----------|-----|
| Partner selection | All reads sharing a k-mer | Fixed pair from PAF |
| Hit emission | `readId != readIdA` | `readId == readIdB` |
| Strand | Both orientations | PAF-specified only |
| Post-filter | Per-read max_n_chain | Per-read max_n_chain + dedup |
| Deduplication | Implicit (one chain per pair) | Explicit best-per-partner |

The PAF path also applies an additional deduplication step: when multiple PAF
records produce overlaps between the same (query, partner) pair, only the
highest-scoring overlap is retained. This deduplication happens before
max_n_chain filtering so the per-read cap is not wasted on duplicate chains.

## 8. Threading Model

### 8.1 Discovery Path Threading

The discovery path uses the `InvertedIndexFinder` class, which extends
`MultithreadedObject` for dynamic batch scheduling:

1. **Batch partitioning**: Reads are divided into batches of ~100 reads each.
   Each thread picks batches from a shared atomic counter (dynamic scheduling),
   which provides natural load balancing — threads that finish fast pick up
   more work.

2. **Thread-local accumulation**: Each thread has its own:
   - `ThreadScratchpad` — reusable DP arrays, hit buffers, chain workspace
   - `vector<OrientedReadPair>` — accumulated candidate pairs
   - `vector<Alignment>` — accumulated alignments
   - `vector<int32_t>` — accumulated chain scores

   The scratchpad's `clear()` method resets all vectors without releasing
   memory (`vector::clear()` preserves capacity), so allocation cost is
   amortized across reads.

3. **Single merge**: After all threads finish, results are merged into the
   global output vectors in a single sequential pass. The global vectors are
   pre-sized to the total count, and each thread's results are copied
   contiguously.

### 8.2 PAF Path Threading

The PAF path uses raw `std::thread` with the same batch-scheduling mechanism.
Each thread processes a range of imported PAF pairs independently. After all
threads join, results are merged, deduplicated, and filtered sequentially.

### 8.3 Index Construction Threading

Index construction (Phase 2 counting, Phase 1 canonicalization, Phase 5
scatter) uses `std::thread` with static range partitioning — each thread
processes `readCount / threadCount` reads. The merge phase (Phase 3) is
single-threaded because it touches only distinct k-mers, not all occurrences.

## 9. Parameters

### 9.1 Marker Extraction Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Kmers.k` | 50 | K-mer length in RLE space |
| `Kmers.syncmerS` | 11 | Sub-kmer size for closed syncmer selection |
| `Kmers.useSimdClosedSyncmers` | true | Enable SIMD-accelerated syncmer extraction |

### 9.2 Frequency Filtering Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `minFreq` | 2 | Minimum k-mer frequency (removes singletons) |
| `invertedIndexMaxKmerCount` | 2000 | Maximum k-mer frequency (hifiasm `max_kmer_cnt`) |

### 9.3 Hit Weighting Parameters

| Parameter | Default | Hifiasm Equivalent | Description |
|-----------|---------|-------------------|-------------|
| `invertedIndexWeightExponent` | 1.1 | — | Exponent for repetitive k-mer weight LUT |
| `invertedIndexLowFreqMultiplier` | 0.333 | `HA_KMER_GOOD_RATIO` | `lowFreq = coveragePeak × 0.333` |
| `invertedIndexHighFreqMultiplier` | 1.667 | `2.0 - HA_KMER_GOOD_RATIO` | `highFreq = coveragePeak × 1.667` |
| `invertedIndexRareKmerWeight` | 2 | — | Weight for rare k-mers (count ≤ lowFreq) |

### 9.4 High-Frequency Downsampling Parameters

| Parameter | Default | Hifiasm Equivalent | Description |
|-----------|---------|-------------------|-------------|
| `invertedIndexDownsampleHighFrequencyMarkers` | true | — | Enable streak-based downsampling |
| `invertedIndexHighFrequencySampleDistance` | 500 | `sample_dist` | Spacing between retained markers |
| `invertedIndexMaxHighFrequencyPerStreak` | 16 | `MAX_MAX_HIGH_OCC` | Hard cap on markers per streak |

### 9.5 Chaining and Post-Filter Parameters

| Parameter | Default | Hifiasm Equivalent | Description |
|-----------|---------|-------------------|-------------|
| `driftRateTolerance` | 0.05 | `bw_rate` | Bandwidth rate for ONT reads |
| `invertedIndexLchainIsAccurate` | true | `is_accurate` | Use ONT-optimized DP parameters |
| `invertedIndexUseEcScoring` | true | — | Use `comput_sc_ch_ec` scoring |
| `invertedIndexHighFactor` | 5.0 | `high_factor` | `max_n_chain = max(coveragePeak × 5, 100)` |
| `invertedIndexMinNChain` | 100 | `MIN_N_CHAIN` | Minimum max_n_chain |
| `invertedIndexEnableMcopyFast` | true | — | Enable multi-copy chain extraction |
| `invertedIndexMcopyNum` | 3 | `mcopy_num` | Max chains per read pair |
| `invertedIndexMcopyRate` | 0.70 | `mcopy_rate` | Min score ratio for secondary chains |
| `invertedIndexMcopyKhitCutoff` | 32 | `mcopy_khit_cut` | Min anchors to enable mcopy |
| `invertedIndexMcopyOcvWindow` | 3072 | `COV_W` | Window size for coverage rescue |

## 10. Complexity Analysis

### 10.1 Index Construction

| Phase | Time | Space |
|-------|------|-------|
| Canonicalization | O(N) parallel | O(N) for cache arrays |
| Counting | O(N) parallel | O(D × T) for thread tables |
| Merge | O(D) sequential | O(D) for global table |
| Prefix sum | O(D) sequential | O(D) for scatter states |
| Scatter | O(N) parallel | O(N) for compact array |
| Query table | O(D) sequential | O(D) for hash table |

Where N = total marker occurrences, D = distinct k-mers, T = thread count.

### 10.2 Hit Collection (per read)

| Operation | Time | Space |
|-----------|------|-------|
| Hash table lookups | O(M) expected | O(1) per lookup |
| Hit emission | O(M × C̄) | O(H) for flat hits |
| Streak downsampling | O(S log S) per streak | O(S) per streak |
| Radix sort | O(H) | O(H) for tmp buffer |

Where M = markers in query read, C̄ = average k-mer count, H = total hits,
S = max streak size.

### 10.3 Overall

For a dataset with R reads, M markers per read, and C̄ average k-mer count:

- **Index construction**: O(R × M) work, parallelized across T threads
- **Hit collection**: O(R × M × C̄) work, parallelized across T threads
- **Memory**: O(R × M × 8) bytes for compact array + O(D × 30) for hash table

The hit collection phase dominates runtime. With typical parameters (M ≈ 500
markers per read, C̄ ≈ 30), each read generates ~15,000 hits before
downsampling.

## 11. Implementation Files

| File | Contents |
|------|----------|
| `AssemblerMarkers.cpp` | Marker extraction (syncmers, frequency filtering, coverage peak) |
| `InvertedIndexBuilder.hpp` | Count-then-scatter index construction (all 6 phases) |
| `AssemblerInvertedIndex.cpp` | Hit collection, weighting, downsampling, grouping, chaining orchestration |
| `AlignmentCanonicalization.hpp` | Scoring functions (`comput_sc_ch`, `comput_sc_ch_ec`, `cal_bw`, `NORMAL_W`) |
| `hifiasmCoordinateTransforms.hpp` | RC interval coordinate conversion |
| `Assembler.hpp` | Data structure definitions (CompactOccurrence, HashEntry, InvertedIndexData) |
| `AssemblerOptions.hpp` | Parameter definitions and defaults |
| `AssemblerOptions.cpp` | Parameter registration and display |

## 12. Relationship to Hifiasm

The candidate finding pipeline is designed for behavioral parity with
hifiasm v0.25.0-r726 (commit ec9a8b2). The following table maps Dinara
components to their hifiasm counterparts:

### 12.1 K-mer Frequency Handling

| Dinara | Hifiasm | Description |
|--------|---------|-------------|
| `minFreq = 2` | `ha_ct_shrink(ct, 2, ...)` | Remove singleton k-mers |
| `maxFreq = 2000` | `max_kmer_cnt = 2000` in `gen_hh` | Hard cutoff for extremely repetitive k-mers |
| `highFreqThreshold` | `5 × hom_cov` in `ha_flt_tab` | Threshold for downsampling annotation |
| Weight LUT | `normal_w` macro | Frequency-based score normalization |

### 12.2 Downsampling

| Dinara | Hifiasm | Description |
|--------|---------|-------------|
| `flushHighFrequencyStreak` | `hf_select` + `select_mz_h` | Per-streak marker selection |
| `sampleDistance = 500` | `sample_dist` | Spacing between retained markers |
| `maxPerStreak = 16` | `MAX_MAX_HIGH_OCC` | Hard cap per streak |
| `keep == 0` → discard all | `max_high_occ == 0` → skip `hf_select` | Short streak handling |
| `count >= span` guard | `rid < pe - ps` guard | Local repetitiveness filter |

### 12.3 Index Construction

| Dinara | Hifiasm | Description |
|--------|---------|-------------|
| Count-then-scatter | `ha_pt_gen` | Position table construction |
| `CompactOccurrence` (8 bytes) | Position table entry | Per-occurrence storage |
| Open-addressing hash table | `ha_pt` hash table | O(1) k-mer lookup |

### 12.4 Hit Collection

| Dinara | Hifiasm | Description |
|--------|---------|-------------|
| `InvertedIndexTempHit` | `k_mer_hit` | Per-hit record |
| `isRcA ⊕ isRcB` | `z->rev ^ y->rev` | Strand determination |
| `chain_cutoff = 2` | `ecovlp.cpp:3274` | Per-strand minimum hit count |
| Radix sort by `(readIdB, posA)` | Radix sort in anchor.cpp | Hit grouping |
