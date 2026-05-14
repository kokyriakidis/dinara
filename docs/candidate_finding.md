# Overlap Candidate Finding in Dinara

This document describes how Dinara discovers overlap candidates between reads.
The pipeline transforms raw sequencing reads into a set of read pairs that
likely share a true genomic overlap, suitable for downstream DP chaining and
alignment.

## Overview

The candidate finding pipeline has four stages:

1. **Marker extraction** — select informative k-mer positions from each read
2. **Inverted index construction** — build a lookup table mapping each k-mer to
   all reads that contain it
3. **Hit collection** — for each read, query the index to find shared k-mers
   with other reads
4. **Hit grouping** — group hits by partner read and strand, producing
   candidate pairs for chaining

The design follows hifiasm's overlap discovery approach (anchor.cpp,
Hash_Table.cpp) but replaces the position table with a count-then-scatter
inverted index that reduces peak memory during construction.

## Stage 1: Marker Extraction

Dinara operates in **run-length encoded (RLE) space**, where homopolymer runs
are collapsed to single bases. This is equivalent to hifiasm's homopolymer
compression (HPC) and reduces noise from homopolymer length errors common in
nanopore sequencing.

### Closed Syncmers

When `Kmers.useSimdClosedSyncmers = true` (the recommended setting), Dinara
extracts markers using **closed syncmers** with parameters:

- **k = 50** — k-mer length in RLE bases
- **s = 11** — sub-kmer (s-mer) size for syncmer selection

A closed syncmer is a k-mer where the minimum s-mer hash occurs at the first
or last position within the k-mer. This selection criterion is:

- **Deterministic** — the same k-mer always produces the same selection decision
- **Strand-symmetric** — a k-mer and its reverse complement are both selected
  or both rejected
- **Approximately uniform** — syncmers sample roughly 2/(k-s+1) of all k-mer
  positions, providing even coverage

With k=50 and s=11, the sampling rate is approximately 2/40 = 5% of positions,
which is comparable to hifiasm's minimizer density with k=51, w=51.

The extraction is SIMD-accelerated (`findMarkersSimdClosedSyncmers`) and
processes reads in parallel across threads.

### Canonicalization

Each extracted k-mer is stored as a `KmerId` — a 128-bit integer encoding the
two-bit-per-base representation. For each k-mer, we compute its reverse
complement and take the lexicographically smaller of the two as the
**canonical form**. This ensures that overlapping reads match regardless of
which strand they originate from.

The canonical k-mer ID and a flag indicating whether the observed k-mer was
the reverse complement are cached in flat arrays
(`strand0CanonicalKmerIds`, `strand0CanonicalIsRc`) indexed by a per-read
offset table (`strand0CanonicalOffsets`). This cache eliminates redundant
reverse-complement computation during hit collection.

### Frequency Filtering

Before index construction, `applyKmerCountFilter` removes markers whose
genome-wide frequency falls outside `[minFreq, maxFreq]`. This eliminates:

- **Singleton k-mers** (sequencing errors) that would never produce useful hits
- **Extremely high-frequency k-mers** (centromeric/telomeric repeats) that
  would dominate memory and compute

The frequency distribution is analyzed to find the coverage peak
(`coveragePeak`), which represents the expected k-mer frequency for
single-copy genomic sequence. The peak is used downstream to set
coverage-relative thresholds for hit weighting.

## Stage 2: Inverted Index Construction

The inverted index maps each canonical k-mer to the list of (readId, position)
pairs where it occurs. Only strand-0 markers are indexed because canonical
k-mers are strand-symmetric — reverse-complement matches are detected during
chaining by observing the RC flags.

### Count-Then-Scatter Algorithm

The index is built using a **count-then-scatter** approach inspired by
hifiasm's position table construction. This avoids materializing a large
intermediate array and sorting it.

**Phase 1 — Parallel counting:**
Each thread processes a range of reads and builds a thread-local hash table
mapping canonical k-mer → occurrence count. No positions are stored in this
phase. The thread-local tables are sized at 4× the expected distinct k-mer
count to keep collision rates low.

**Phase 2 — Merge:**
Thread-local tables are merged into a single global counting table. This
touches only distinct k-mers (not all occurrences), so it is fast even
single-threaded.

**Phase 3 — Prefix sum:**
A sequential scan over the global table assigns each k-mer a contiguous range
in the final compact array. The starting offset for k-mer *i* is the sum of
counts for all k-mers before it.

**Phase 4 — Parallel scatter:**
Each thread re-iterates its reads. For each marker, it looks up the k-mer in
the global table and atomically claims a slot in the pre-allocated compact
array using `fetch_add`. The compact array stores 8-byte entries:

```
struct CompactOccurrence {
    ReadId   readId;    // 4 bytes
    uint32_t position;  // 4 bytes (top bit encodes RC flag)
};
```

**Phase 5 — Query hash table:**
A separate open-addressing hash table is built for O(1) k-mer lookups during
hit collection. Each entry stores:

```
struct HashEntry {
    KmerId   key;    // 16 bytes (canonical k-mer)
    uint64_t start;  // offset into CompactOccurrence array
    uint32_t count;  // number of occurrences
    bool     empty;  // slot occupancy flag
};
```

### Memory Comparison

| Metric | Old (sort-based) | New (count-scatter) |
|--------|-----------------|---------------------|
| Intermediate per occurrence | 24 bytes (KmerId + ReadId + position) | 0 bytes (no intermediate) |
| Sort buffer | 24 bytes (second array for radix sort) | 0 bytes (no sort) |
| Peak transient | ~48 bytes/occurrence | ~8 bytes/occurrence |
| Final (query phase) | 8 bytes/occurrence | 8 bytes/occurrence |
| Radix sort passes | ceil(2k/8) = 13 for k=50 | 0 |

For a dataset with 100M marker occurrences, the old approach required ~4.5 GB
of transient memory for sorting alone. The new approach allocates only the
final 800 MB compact array.

## Stage 3: Hit Collection

For each query read A, the hit collection phase finds all reads that share at
least one k-mer with A. This is the most compute-intensive part of candidate
finding.

### Per-Read Processing

For each marker in read A:

1. Look up the canonical k-mer in the hash table → get (start, count) into the
   compact occurrence array
2. Iterate all occurrences, skipping self-hits (readId == readIdA)
3. For each hit, emit a flat hit record:

```
struct FlatHit {
    ReadId   readIdB;   // partner read
    uint32_t posA;      // position in read A
    uint32_t posB;      // position in read B
    uint32_t ordinalA;  // marker ordinal in read A
    uint32_t weight;    // frequency-based weight
    uint8_t  isRcA;     // RC flag for read A's k-mer
    uint8_t  isRcB;     // RC flag for read B's k-mer
};
```

### Frequency-Based Weighting

Each hit is assigned a weight based on the genome-wide frequency of its k-mer.
The weight normalizes the contribution of each anchor during DP chaining,
following hifiasm's `normal_w` convention:

| Frequency tier | Condition | Weight |
|---------------|-----------|--------|
| Rare | count ≤ lowFreqThreshold | 2 (more informative) |
| Normal | lowFreq < count < highFreq | 1 (baseline) |
| Repetitive | count ≥ highFreqThreshold | pow(normalized_count, 1.1) |

The thresholds are derived from the coverage peak:
- `lowFreqThreshold = coveragePeak * lowFreqMultiplier`
- `highFreqThreshold = coveragePeak * highFreqMultiplier`

Rare k-mers receive double weight because they are more likely to represent
unique genomic positions. Repetitive k-mers receive increasing weight to
account for the fact that their score contribution is divided by their
occurrence count during chaining (via `normal_w(score, count)`).

### High-Frequency Marker Downsampling

Consecutive high-frequency markers (repetitive regions like tandem repeats)
can produce an explosion of hits. To control memory and compute, Dinara
implements hifiasm's streak-based downsampling:

1. **Buffer** consecutive high-frequency markers into a streak
2. When the streak ends, compute:
   - `span = rightBoundary - leftBoundary` (genomic span of the streak)
   - `keep = round(span / sampleDistance)` where `sampleDistance = 500`
   - Cap at `maxPerStreak = 16`
3. **Select** the `keep` markers with the smallest (count, hash) keys using
   `nth_element` — this prefers the most informative markers within the repeat
4. **Filter**: discard any selected marker whose count ≥ span (too repetitive
   even relative to the local context)
5. **Emit** selected markers in position order

This reduces a 5 kb tandem repeat with 100 high-frequency markers to ~10
evenly-spaced markers, preserving positional signal while controlling memory.

## Stage 4: Hit Grouping

After collecting all hits for read A, they are sorted by a packed key
`(readIdB << 32 | posA)` using LSD radix sort. This groups hits by partner
read, and within each partner read, orders them by query position.

Within each partner-read group, hits are split by **strand orientation**:

- **Same-strand** (isRcA == isRcB): target positions increase with query
  positions
- **Opposite-strand** (isRcA != isRcB): target positions decrease with query
  positions

Each strand bucket is then converted to hifiasm-format anchor arrays and
passed to the DP chaining phase (described in the chaining document).

### Strand Determination

The strand of a hit is determined by comparing the RC flags of the k-mer
observation in reads A and B:

- If both reads observed the same canonical orientation (both forward or both
  RC), the reads overlap on the same strand
- If they observed opposite orientations, the reads overlap on opposite
  strands

This is equivalent to hifiasm's `z->rev == y->rev` check in anchor.cpp.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Kmers.k` | 50 | K-mer length in RLE space |
| `Kmers.syncmerS` | 11 | Sub-kmer size for closed syncmer selection |
| `Kmers.useSimdClosedSyncmers` | true | Enable SIMD-accelerated syncmer extraction |
| `OverlapCandidates.method` | InvertedIndex | Candidate finding method |
| `OverlapCandidates.driftRateTolerance` | 0.05 | Bandwidth rate for chaining (ONT) |

## Implementation Files

| File | Contents |
|------|----------|
| `AssemblerMarkers.cpp` | Marker extraction (syncmers, minimizers, frequency filtering) |
| `InvertedIndexBuilder.hpp` | Count-then-scatter index construction |
| `AssemblerInvertedIndex.cpp` | Hit collection, grouping, and chaining orchestration |
| `Assembler.hpp` | Data structure definitions (CompactOccurrence, HashEntry, etc.) |
