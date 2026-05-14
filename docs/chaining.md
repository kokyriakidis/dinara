# Overlap Chaining in Dinara

This document describes how Dinara chains k-mer anchors into overlap
candidates.

## 1. Overview

After candidate finding produces a set of k-mer hits between read pairs, the
chaining phase determines which hits form coherent, collinear chains that
represent true genomic overlaps. The chaining algorithm is a faithful port of
hifiasm's ONT error correction path (ecovlp.cpp, anchor.cpp, Hash_Table.cpp)
adapted to Dinara's inverted index data structures.

The pipeline has four stages:

1. **Anchor preparation** — convert hits to hifiasm-format anchors with
   end-position coordinates, split by strand, sort by `(self_offset, offset)`.
2. **DP chaining** — O(N × max_iter) dynamic programming with three pruning
   strategies (lookback window, max-skip, sliding distance window) to find
   optimal anchor chains.
3. **Multi-copy extraction** — extract up to `mcopy_num = 3` high-scoring
   independent chains per read pair for repetitive regions.
4. **Post-filtering** — apply per-overlap-type caps (max_n_chain), coverage
   window rescue (COV_W), and weak-chain suppression (R485).

### Pipeline Diagram

```
  K-mer hits (from candidate finding)
       │
       ▼
  ┌─────────────────────────────────┐
  │  Convert to HifiasmKmerHit      │  Start → end positions
  │  Split by strand (0/1)          │  RC coordinate transform
  │  Sort by (self_offset, offset)  │  Run-based tie-breaking
  └─────────────┬───────────────────┘
                │
                ▼
  ┌─────────────────────────────────┐
  │  Quick check (O(N))             │  Greedy collinear scan
  │  → narrows DP range [si, ei)   │  Skips full DP if monotonic
  └─────────────┬───────────────────┘
                │
                ▼
  ┌─────────────────────────────────┐
  │  Main DP loop (O(N × max_iter)) │  f[i] = max(f[j] + score(i,j))
  │  + max_skip pruning             │  Backtrack via p[i] pointers
  │  + sliding distance window       │  max_ii rescue for skipped anchors
  └─────────────┬───────────────────┘
                │
                ▼
  ┌─────────────────────────────────┐
  │  Backtrack best chain           │  Follow p[i] from argmax(f[i])
  │  Multi-copy extraction          │  Up to 3 independent chains
  └─────────────┬───────────────────┘
                │
                ▼
  ┌─────────────────────────────────┐
  │  Post-filter                    │
  │  Stage 1: max_n_chain per type  │  Per-overlap-type score cap
  │  Stage 2: COV_W window rescue   │  Coverage-aware type-3 rescue
  │  Stage 3: R485 weak suppression │  Remove dominated short chains
  └─────────────┬───────────────────┘
                │
                ▼
  Overlap candidates (HifiasmOverlapRegion → Dinara Alignment)
```

## 2. Anchor Preparation

Each k-mer hit between reads A and B is converted to a `HifiasmKmerHit`
anchor, which is the fundamental unit of the DP chaining algorithm.

### 2.1 HifiasmKmerHit Structure

```cpp
struct HifiasmKmerHit {
    uint32_t readID;       // Target read ID (y_id in hifiasm)
    uint32_t self_offset;  // Query end position (x coordinate)
    uint32_t offset;       // Target end position (y coordinate)
    uint32_t cnt;          // Bit-packed: (weight << 8) | span
    uint8_t  strand;       // 0 = same-strand, 1 = opposite-strand

    // Dinara extensions (not used by DP):
    uint32_t ordinalA;     // Query marker ordinal
    uint32_t ordinalB;     // Target marker ordinal
    uint32_t globalIndex;  // Index into global hit array
};
```

**Bit-packing scheme for `cnt`:**

```
cnt = (occurrence_weight << 8) | k-mer_span

Bits [7:0]   = k-mer span (min(kmerLen, 255))
Bits [31:8]  = occurrence weight (from frequency-based weighting)
```

The span is capped at 255 because it occupies only 8 bits. The weight
occupies the upper 24 bits, allowing values up to 16,777,215.

### 2.2 Coordinate System

Hifiasm uses **end positions** (inclusive) for anchor coordinates, while
Dinara's markers store **start positions**. The conversion is:

```
self_offset = posA + (seedSpan - 1)     // query end position
```

For same-strand hits:
```
offset = posB + (seedSpan - 1)          // target end position
```

For opposite-strand hits, the target position is reflected into
reverse-complement coordinates:
```
offset = readLenB - 1 - posB            // RC-transformed end position
```

This reflection ensures that in a true opposite-strand overlap, both
`self_offset` and `offset` increase monotonically along the chain, allowing
the same DP recurrence to handle both orientations.

### 2.3 Strand Splitting and Sorting

Anchors are split into two groups by strand orientation:

- **Strand 0 (same-strand)**: `isRcA ⊕ isRcB = 0`. Both reads observed the
  k-mer in the same canonical orientation.
- **Strand 1 (opposite-strand)**: `isRcA ⊕ isRcB = 1`. Reads observed
  opposite orientations.

Within each strand group, anchors are sorted by `(self_offset, offset)` —
query position first, then target position for ties. Since hits arrive
pre-sorted by `posA` from the radix sort, `self_offset` is already
non-decreasing. Only ties need sub-sorting, which is done by
`sortHifiasmHitsBySelfOffsetThenOffsetRuns` in O(n) expected time.

The two strand groups are concatenated (strand 0 first, strand 1 second)
into a single array. The DP loop uses the `strand` field to skip anchors
on the wrong strand during predecessor search.

## 3. DP Chaining

The core chaining algorithm finds the highest-scoring chain of collinear
anchors. It is a direct port of hifiasm's `lchain_qdp_mcopy_fast`
(Hash_Table.cpp:2097).

### 3.1 DP Parameters

The DP parameters are configured by `getHifiasmLchainDpOptions`, which
implements hifiasm's `set_lchain_dp_op` (anchor.cpp:2272). For ONT reads
(`isAccurate = true`):

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `max_skip` | 25 | Max consecutive predecessors to skip before pruning |
| `max_iter` | 5000 | Lookback window size |
| `max_dis` | 5000 | Max gap distance for sliding window |
| `quick_check` | true | Enable O(N) fast path |
| `chn_pen_gap` | `0.5 × exp(-0.01 × k)` | Gap penalty coefficient |
| `chn_pen_skip` | `0.0005 × exp(-0.01 × k)` | Skip penalty coefficient |

**Penalty formula:**

```
penalty_coefficient = base × exp(-div × k)
```

Where `div = 0.01` for ONT (slower decay, higher penalties for noisy reads)
and `div = 0.1` for HiFi (faster decay, lower penalties for accurate reads).

The exponential decay with k-mer length reflects the fact that longer k-mers
produce fewer, more widely-spaced anchors, so the gap penalty should be
reduced to avoid over-penalizing natural spacing.

**Example values for k=50:**
- `chn_pen_gap = 0.5 × exp(-0.01 × 50) = 0.5 × 0.6065 ≈ 0.303`
- `chn_pen_skip = 0.0005 × exp(-0.01 × 50) ≈ 0.000303`

### 3.2 Scoring Function (`comput_sc_ch_ec`)

The scoring function evaluates the quality of connecting anchor *i* to
predecessor anchor *j*. It returns `INT32_MIN` for invalid connections or a
positive score for valid ones.

**Reference:** hifiasm Hash_Table.cpp:1515-1541

#### Step 1 — Monotonicity Check

```
dq = a[i].self_offset - a[j].self_offset    (query gap)
dr = a[i].offset - a[j].offset              (target gap)
```

Both `dq` and `dr` must be strictly positive. If either is ≤ 0, the anchors
are not collinear → return `INT32_MIN`.

#### Step 2 — Bandwidth Check

```
dd = |dq - dr|    (diagonal deviation)
```

If `dd > 16` and `dd > cal_bw(a[i], a[j], bw_rate, readLenA, readLenB)`,
the connection deviates too far from the expected diagonal → return
`INT32_MIN`.

The bandwidth `cal_bw` computes the maximum allowed diagonal deviation based
on the alignable span of the overlap region:

```
sf_s = a[j].self_offset          (query start)
sf_e = a[i].self_offset + 1      (query end)
ot_s = a[j].offset               (target start)
ot_e = a[i].offset + 1           (target end)

// Left-normalize: shift query start to 0 or offset by target start
if sf_s ≤ ot_s: sf_s = 0
else:           sf_s -= ot_s

// Right-extend: extend query end to sequence boundary
sf_r = readLenA - sf_e
ot_r = readLenB - ot_e
if sf_r ≤ ot_r: sf_e = readLenA
else:           sf_e += ot_r

bandwidth = (sf_e - sf_s) × bw_rate
```

With `bw_rate = 0.05` (ONT), a 10 kb overlap allows up to 500 bp of diagonal
drift. Near sequence ends, the bandwidth is larger because the alignable span
extends to the read boundaries.

The `dd > 16` short-circuit avoids the `cal_bw` computation for small
deviations, which are always accepted.

#### Step 3 — Base Score

```
dg = min(dq, dr)                    (effective gap length)
q_span = a[i].cnt & 0xFF           (k-mer span, typically k)
score = min(q_span, dg)            (base score)
```

The base score is the minimum of the k-mer span and the gap length. This
ensures that overlapping anchors (dg < q_span) don't double-count bases.

#### Step 4 — Occurrence Normalization

```
weight = a[i].cnt >> 8             (occurrence weight)
score = normal_w(score, weight)    = max(1, score / weight)
```

High-frequency k-mers contribute less per anchor. A k-mer appearing 100 times
in the genome has weight ≈ `⌊(1 + ⌈100/(2×highFreq)⌉)^1.1⌋`, so its score
contribution is divided by this factor.

#### Step 5 — Gap Penalty (Split Strategy)

The penalty has two components:

```
lin_pen = chn_pen_gap × dd         (linear penalty, proportional to deviation)
a_pen = score × (dd / dg) / bw_rate  (adaptive penalty, proportional to drift ratio)
```

The split strategy selects between them based on gap size:

```
if dd < 4:   penalty = min(lin_pen, a_pen)   // small gaps: be lenient
if dd ≥ 4:   penalty = max(lin_pen, a_pen)   // large gaps: be strict
```

This is the key difference between `comput_sc_ch_ec` (ONT) and `comput_sc_ch`
(HiFi). The HiFi version always uses `min(lin_pen, a_pen)` regardless of gap
size, because HiFi reads have fewer indels and don't need the stricter
large-gap penalty.

A skip penalty is added proportional to the gap length:

```
penalty += chn_pen_skip × dg
```

This discourages chains that skip many intermediate anchors, preferring dense
chains.

#### Step 6 — Final Score

```
chain_score = score - penalty
```

### 3.3 Quick Check (O(N) Fast Path)

Before running the full O(N × max_iter) DP, a linear-time quick check
(`hifiasm_quick_ck_lchain`) attempts to find a high-quality chain by scanning
anchors left-to-right and greedily extending collinear runs.

**Reference:** hifiasm Hash_Table.cpp:2007-2095

**Algorithm:**

1. Partition the anchor array into **strand-homogeneous segments** — maximal
   contiguous runs where all anchors have the same strand value.

2. For each segment, check if it is **monotonic**: both `self_offset` and
   `offset` strictly increase from one anchor to the next.

3. If a segment is monotonic, attempt **greedy chaining**: link each anchor
   to its immediate predecessor using the same `comput_sc_ch_ec` scoring.
   If any connection scores below the anchor's own span (starting fresh is
   better), stop.

4. If the entire segment is successfully chained and the best anchor is the
   last one, validate the cumulative diagonal drift against the bandwidth.
   If valid, the segment is **solved** — the full DP can skip it.

5. The quick check narrows the unsolved DP range `[si, ei)`:
   - If the solved segment is at the beginning: `si = segment_end`
   - If the solved segment is at the end: `ei = segment_start`

**When the quick check helps:** For high-quality reads with dense, collinear
anchors, the quick check often solves the entire anchor array in O(N) time,
completely skipping the O(N × max_iter) DP. This provides a significant
speedup for the common case.

**When it doesn't help:** If anchors are not monotonic (e.g., due to
structural variants, chimeric reads, or high error rates), the quick check
falls through to the full DP with `si = 0, ei = N`.

### 3.4 Main DP Loop

The main DP loop processes anchors in query-position order over the unsolved
range `[si, ei)`:

**Reference:** hifiasm Hash_Table.cpp:2097-2170

```
for i = si to ei - 1:
    // Restrict lookback window
    st = max(si, i - max_iter)
    while a[st].strand ≠ a[i].strand: st++

    max_f = a[i].cnt & 0xFF    // base score (anchor's own span)
    max_j = -1                  // best predecessor
    n_skip = 0                  // consecutive skips

    // Inner loop: find best predecessor
    for j = i-1 down to st:
        if a[j].strand ≠ a[i].strand: continue

        sc = comput_sc_ch_ec(a[i], a[j], ...)
        if sc == INT32_MIN: continue
        sc += f[j]

        if sc > max_f:
            max_f = sc
            max_j = j
            if n_skip > 0: n_skip--
        else if t[j] == i:        // j's predecessor was visited
            if ++n_skip > max_skip: break

        if p[j] ≥ 0: t[p[j]] = i  // mark j's predecessor as visited

    f[i] = max_f
    p[i] = max_j
```

**Three pruning strategies:**

1. **Lookback window** (`max_iter = 5000`): The inner loop only considers
   the most recent 5000 anchors as potential predecessors. This limits
   worst-case complexity to O(N × 5000) instead of O(N²).

2. **Max-skip** (`max_skip = 25`): Once 25 consecutive predecessors fail to
   improve the score, the inner loop terminates. The "consecutive" count is
   tracked via the visit array `t[]`: a predecessor `j` counts as a "skip"
   only if `t[j] == i`, meaning `j`'s own predecessor was already considered
   for anchor `i`. This prevents O(N²) behavior in practice — most anchors
   find their best predecessor within a few iterations.

3. **Sliding distance window** (`max_dis = 5000`): A cached best-score
   anchor `max_ii` within `max_dis` base pairs is maintained. If the
   max-skip pruning terminates the inner loop before reaching `max_ii`,
   a rescue attempt tries `max_ii` as a predecessor. This prevents the
   max-skip heuristic from missing good predecessors that are far back in
   the array but close in genomic distance.

**Tie-breaking:** When multiple anchors have the same DP score, the algorithm
prefers the one with the shorter effective overlap length (computed by
`get_chainLen`). This favors tighter, more compact overlaps.

### 3.5 Backtracking

After the DP loop, the best chain is extracted by following predecessor
pointers from the anchor with the highest score:

```
msc_i = argmax(f[i])
msc = f[msc_i]

// Backtrack
k = msc_i
cL = 0
while k ≥ 0:
    ii[k] = 1          // mark as in-chain
    t[cL++] = k        // store in reverse order
    k = p[k]

// t[0..cL-1] contains chain indices in reverse order
// (t[0] = last anchor, t[cL-1] = first anchor)
```

The chain membership array `ii[]` is used by the multi-copy extraction phase
to identify which anchors are already claimed.

## 4. Multi-Copy Extraction

For read pairs with multiple valid overlaps (e.g., reads spanning a tandem
duplication or segmental duplication), the algorithm extracts up to
`mcopy_num = 3` independent chains.

**Reference:** hifiasm Hash_Table.cpp:2170-2284

### 4.1 Discovery Path (Hifiasm-Native mcopy)

In the discovery path, multi-copy extraction is performed inside
`hifiasm_lchain_qdp_mcopy_fast` using hifiasm's native algorithm:

**Precondition:** The best chain must have at least `mcopy_khit_cutoff = 32`
anchors. If not, only the best chain is emitted.

**Algorithm:**

1. **Normalize scores**: Subtract the minimum DP score (`plus`) from all
   `f[i]` values to make them non-negative. Compute the minimum acceptable
   score: `min_sc = msc × mcopy_rate` (default: 0.70 × best score).

2. **Collect candidates**: Pack `(score << 32 | index << 1)` into the `t[]`
   array for all anchors with `f[i] ≥ min_sc` that are not in the best chain.
   Sort by packed value (ascending, so highest scores are at the end).

3. **Extract chains**: Iterate from highest to lowest score. For each
   candidate:
   - Follow predecessor pointers until hitting a used anchor (marked by
     `t[i] |= 1`)
   - Compute the independent score: `sc = candidate_score - f[used_ancestor]`
   - If `sc ≥ min_sc` and the chain has more than 1 anchor, accept it
   - Mark all anchors in the accepted chain as used
   - Stop after `mcopy_num` chains

4. **Emit chains**: Each accepted chain is converted to a
   `HifiasmOverlapRegion` via `push_ovlp_chain_qgen`.

### 4.2 PAF Path (applyMcopyFastSelection)

The PAF path uses a separate `applyMcopyFastSelection` function that operates
on the `ChainCandidate` array produced by the PAF path's DP:

1. Sort candidates by `(score descending, chainLen ascending)`
2. Mark all anchors in the best chain as used
3. For each subsequent candidate with `score ≥ best × mcopyRate`:
   - Trace backward until hitting a used anchor
   - Compute adjusted score (subtract the used ancestor's score)
   - Accept if adjusted score ≥ threshold and chain has > 1 anchor
4. Replace the candidate array with the selected chains

### 4.3 Chain-to-Overlap Conversion

Each extracted chain is converted to an overlap region via
`hifiasm_push_ovlp_chain_qgen` (Hash_Table.cpp:1752), which projects the
chain's anchor coordinates to read-level overlap boundaries:

**Left normalization** (one start becomes 0):

```
if x_pos_s ≤ y_pos_s:
    y_pos_s -= x_pos_s
    x_pos_s = 0
else:
    x_pos_s -= y_pos_s
    y_pos_s = 0
```

**Right extension** (one end reaches its sequence boundary):

```
xr = readLenA - x_pos_e - 1    (remaining bases on query)
yr = readLenB - y_pos_e - 1    (remaining bases on target)

if xr ≤ yr:
    x_pos_e = readLenA - 1     (query reaches its end)
    y_pos_e += xr
else:
    y_pos_e = readLenB - 1     (target reaches its end)
    x_pos_e += yr
```

This projection estimates the full overlap extent from the chain endpoints,
assuming the overlap extends to whichever read boundary is closer.

### 4.4 Effective Overlap Length (get_chainLen)

The effective overlap length is used for tie-breaking when chains have
identical scores. It is computed by `hifiasm_get_chainLen`
(Hash_Table.cpp:779):

```
// Apply the same left-normalization and right-extension
// as push_ovlp_chain_qgen, then return x_end - x_begin + 1
chainLen = x_pos_e - x_pos_s + 1
```

Among equal-score chains, the one with the shorter effective length is
preferred (tighter, more compact overlap).

## 5. Post-Filtering

The post-filter (`hifiasm_lchain_qgen_mcopy_fast_postfilter`) applies three
stages of filtering to control the number of overlap candidates per read.
This is a direct port of hifiasm's postfilter in `lchain_qgen_mcopy_fast`
(anchor.cpp:1953-2100).

**Reference:** hifiasm anchor.cpp:1920-2100

### 5.1 Stage 1: max_n_chain Per-Type Cap

Overlaps are classified into four types based on how the query read's
coordinates relate to the overlap boundaries:

| Type | Condition | Meaning | Typical scenario |
|------|-----------|---------|-----------------|
| 0 | `x_pos_s = 0, x_pos_e < len-1` | Prefix overlap | Read A's left end overlaps |
| 1 | `x_pos_s > 0, x_pos_e = len-1` | Suffix overlap | Read A's right end overlaps |
| 2 | `x_pos_s = 0, x_pos_e = len-1` | Contained | Read A is fully contained in B |
| 3 | `x_pos_s > 0, x_pos_e < len-1` | Internal/Containing | Read A contains Read B |

The classification uses `hifiasm_ha_ov_type` (anchor.cpp:86), which checks
the query coordinates against 0 and `readLen - 1`.

If the total number of overlaps exceeds `max_n_chain` (default:
`max(coveragePeak × 5, 100)`):

1. Sort overlaps by score (descending)
2. Count overlaps per type
3. For each type, record the score threshold when `max_n_chain` overlaps of
   that type have been seen
4. Keep an overlap only if its score ≥ the threshold for its type

This ensures that each overlap type retains at most `max_n_chain` candidates,
preventing a single type (typically type 3 in repetitive regions) from
monopolizing the output.

### 5.2 Stage 2: COV_W Window Rescue

If type-3 (internal/containing) overlaps exceed `max_n_chain`, an additional
window-based filter rescues below-threshold overlaps that cover under-served
regions of the query read.

**Window setup:**

```
window_count = ⌈readLen / COV_W⌉        (COV_W = 3072 bp)
```

Each window tracks capacity and usage in a packed 64-bit value:

```
cc[w] = (capacity << 32) | used

capacity = min(window_span × (max_n_chain >> 1), UINT32_MAX)
```

**Rescue criterion:** A below-threshold type-3 overlap is rescued if ≥ 70%
of its query span falls in windows that are not yet at capacity:

```
for each window w overlapping the overlap's query span:
    overlap_in_window = intersection(overlap_span, window_span)
    if used[w] + overlap_in_window ≥ capacity[w]:
        bad += overlap_in_window
    else:
        good += overlap_in_window

rescue = (good ≥ 0.7 × (good + bad))
```

This prevents regions with extreme coverage from being over-represented while
ensuring that low-coverage windows retain their overlaps. The 70% threshold
means an overlap is rescued only if most of its span covers under-served
regions.

**Window update:** When an overlap is kept (either above threshold or
rescued), its contribution is added to the usage counters of all windows it
spans.

### 5.3 Stage 3: R485 Weak-Chain Suppression

The final filter removes weak chains that are likely false positives. A chain
is considered **weak** if its anchor count (`align_length`) is below
`chain_cutoff = 2`.

**Precondition:** The `lch` flag must be set, indicating that at least one
weak chain exists among the surviving overlaps. If no weak chains exist, this
stage is skipped entirely.

**Algorithm:** For each weak chain *w*:

1. Compute the weak chain's query span `[zs, ze)` and overlap threshold:
   ```
   ob = max(16, ⌊(ze - zs) × 0.95⌋)    (95% overlap required)
   ```

2. Compute dominance thresholds:
   ```
   osc = w.shared_seed × 16             (strong chain needs 16× score)
   ocn = w.align_length << 4            (strong chain needs 16× anchors)
   ```

3. Search for a **strong chain** (anchor count ≥ `chain_cutoff`) that:
   - Overlaps the weak chain's query span by at least `ob` bases
   - Has `align_length ≥ ocn` (enough anchors)
   - Has `shared_seed ≥ osc` (much better score)
   - Has at least `ocn` anchors within the overlap zone (verified by
     scanning the strong chain's anchor positions)

4. If such a strong chain exists, **suppress** the weak chain.

**Anchor position verification:** The strong chain's anchors are stored in
`chainHitIndexFlat` (discovery path) or reconstructed from backtrack pointers
(PAF path). The algorithm counts how many of the strong chain's anchors fall
within the overlap zone `[os, oe)` and requires at least `ocn` of them.

This removes spurious short chains in regions where a much better chain
already exists, reducing false positive overlaps.

### 5.4 Constants

| Constant | Value | Hifiasm Name | Purpose |
|----------|-------|-------------|---------|
| `OFL` | 0.95 | `OFL` | Minimum overlap fraction for suppression |
| `CH_OCC` | 4 | `CH_OCC` | Bit shift for anchor count comparison (`<< 4`) |
| `CH_SC` | 16 | `CH_SC` | Score multiplier for dominance check |
| `chain_cutoff` | 2 | `chain_cutoff` | Weak chain threshold |
| `COV_W` | 3072 | `COV_W` | Window size for coverage rescue |

## 6. Overlap Region Output

Each surviving chain produces a `HifiasmOverlapRegion` that is then converted
to a Dinara `Alignment` object.

### 6.1 HifiasmOverlapRegion Structure

```cpp
struct HifiasmOverlapRegion {
    uint32_t x_id;            // Query read ID
    uint32_t y_id;            // Target read ID
    uint8_t  x_pos_strand;    // Always 0 (query is forward)
    uint8_t  y_pos_strand;    // 0 = same-strand, 1 = reverse-complement

    uint32_t x_pos_s;         // Query start (left-normalized, inclusive)
    uint32_t x_pos_e;         // Query end (right-extended, inclusive)
    uint32_t y_pos_s;         // Target start (inclusive)
    uint32_t y_pos_e;         // Target end (inclusive)

    int32_t  shared_seed;     // Chain DP score
    uint32_t align_length;    // Number of anchors in chain
    uint32_t non_homopolymer_errors;  // Offset into chainHitIndexFlat
};
```

The `non_homopolymer_errors` field is repurposed in Dinara to store the
offset into the `chainHitIndexFlat` array, which records the global indices
of all anchors in the chain. Combined with `align_length`, this defines the
chain's anchor subarray: `chainHitIndexFlat[offset .. offset + align_length)`.

### 6.2 Conversion to Dinara Alignment

The conversion from `HifiasmOverlapRegion` to Dinara `Alignment` involves:

1. **Coordinate system**: Hifiasm uses inclusive endpoints `[x_pos_s, x_pos_e]`.
   Dinara uses half-open intervals `[qs, qe)`. The conversion is:
   ```
   qs = x_pos_s
   qe = x_pos_e + 1
   ```

2. **RC interval transform**: For opposite-strand overlaps, target coordinates
   are in RC space. They are converted to forward-strand coordinates using
   `rcIntervalToForward`:
   ```
   [s, e) → [readLen - e, readLen - s)
   ```

3. **Marker ordinals**: The chain's anchor indices (from `chainHitIndexFlat`)
   are used to look up the `ordinalA` and `ordinalB` fields of each
   `HifiasmKmerHit`. For opposite-strand overlaps, `ordinalB` is reflected:
   ```
   ordinalB_rc = markerCountB - 1 - ordinalB
   ```

4. **Canonicalization**: The candidate pair is canonicalized so that
   `readIds[0] < readIds[1]`. If the pair is swapped, the alignment
   coordinates are also swapped (query ↔ target) and ordinals are adjusted.

### 6.3 Optional Filters

Two optional filters are applied before emitting the final candidate:

**Minimum overlap length** (`minOverlapLength`): Rejects candidates where
`min(querySpan, targetSpan) < minOverlapLength`. This removes very short
overlaps that are unlikely to be useful for assembly.

**Maximum end fuzz** (`maxEndFuzz`): Rejects candidates with large softclips
on both sides. The joint check ensures that for each side (left and right),
at least one of the two reads reaches near the overlap boundary:
```
leftNeed  = min(qs, ts)
rightNeed = min(readLenA - qe, readLenB - te)
reject if leftNeed > maxEndFuzz or rightNeed > maxEndFuzz
```

### 6.4 Output Fields

Each surviving candidate produces:

| Field | Description |
|-------|-------------|
| `OrientedReadPair` | Canonicalized (readId0, readId1, isSameStrand) |
| `Alignment` | Marker ordinal pairs + coordinate intervals |
| `shared_seed` | Chain DP score (used for downstream prioritization) |

## 7. PAF Path Differences

The PAF path (`chainPafCandidates`) uses the same DP scoring and post-filtering
algorithms but differs in several structural ways from the discovery path.

### 7.1 Single-Orientation DP

The PAF record specifies whether the overlap is same-strand or opposite-strand.
The PAF path runs the DP only for the requested orientation, skipping the other.
This halves the DP work per pair compared to the discovery path, which processes
both orientations.

### 7.2 Inline DP (No hifiasm_lchain_qdp_mcopy_fast)

The PAF path does not call `hifiasm_lchain_qdp_mcopy_fast`. Instead, it
implements the DP loop inline using the `ThreadScratchpad` SoA arrays. This
allows it to:

- Run same-strand and diff-strand DP in a single pass over the anchor array
  (though only one orientation is active per pair)
- Use the `runQuickLinearChainPrefix` function for the quick check, which
  operates on SoA arrays rather than the `HifiasmKmerHit` AoS format
- Maintain separate DP state arrays for each strand (`dpSame`, `dpDiff`,
  `parentSame`, `parentDiff`)

### 7.3 Diff-Strand Coordinate Handling

For opposite-strand overlaps in the PAF path, the DP loop uses original
strand-0 coordinates for `posB` (not RC-transformed). The monotonicity check
is reversed: `posB[j] > posB[i]` (target positions decrease as query positions
increase). The bandwidth check converts to forward coordinates on the fly:

```cpp
uint32_t posBi_fwd = readLenB - 1 - posBi;
uint32_t posBj_fwd = readLenB - 1 - posBj;
int32_t bw = hifiasm_cal_bw(posAi, posBi_fwd, posAj, posBj_fwd, ...);
```

### 7.4 Post-Processing Pipeline

After DP chaining, the PAF path applies:

1. **Multi-copy extraction** via `applyMcopyFastSelection`
2. **Weak-chain suppression** (R485) using the same OFL/CH_OCC/CH_SC constants
3. **Best-per-partner deduplication**: Sort by (queryReadId, partnerReadId,
   score desc), keep only the first entry per (query, partner) group
4. **Per-read max_n_chain filtering**: Same per-type cap + COV_W rescue as
   the discovery path

The deduplication step is unique to the PAF path because multiple PAF records
may specify the same read pair with different coordinates.

## 8. Parameters

### 8.1 DP Parameters

| Parameter | ONT Value | HiFi Value | Hifiasm Reference |
|-----------|-----------|------------|-------------------|
| `max_skip` | 25 | 25 | `set_lchain_dp_op` |
| `max_iter` | 5000 | 5000 | `set_lchain_dp_op` |
| `max_dis` | 5000 | 5000 | `set_lchain_dp_op` |
| `quick_check` | true | false | `set_lchain_dp_op` |
| `chn_pen_gap` | `0.5 × e^(-0.01k)` | `0.5 × e^(-0.1k)` | `set_lchain_dp_op` |
| `chn_pen_skip` | `0.0005 × e^(-0.01k)` | `0.0005 × e^(-0.1k)` | `set_lchain_dp_op` |
| `bw_rate` | 0.05 | 0.05 | `driftRateTolerance` |

### 8.2 Multi-Copy Parameters

| Parameter | Default | Hifiasm Reference | Description |
|-----------|---------|-------------------|-------------|
| `mcopyNum` | 3 | `ecovlp.cpp:3274` | Max chains per read pair |
| `mcopyRate` | 0.70 | `ecovlp.cpp:3274` | Min score ratio for secondary chains |
| `mcopyKhitCutoff` | 32 | `ecovlp.cpp:3274` | Min anchors to enable mcopy |

### 8.3 Post-Filter Parameters

| Parameter | Default | Hifiasm Reference | Description |
|-----------|---------|-------------------|-------------|
| `max_n_chain` | `max(coveragePeak × 5, 100)` | `CommandLines.cpp:271` | Per-type overlap cap |
| `chain_cutoff` | 2 | `ecovlp.cpp:3274` | Weak chain threshold |
| `COV_W` | 3072 | `anchor.cpp` | Window size for coverage rescue |
| `OFL` | 0.95 | `anchor.cpp` | Min overlap fraction for suppression |
| `CH_OCC` | 4 | `anchor.cpp` | Anchor count bit shift |
| `CH_SC` | 16 | `anchor.cpp` | Score multiplier for dominance |

## 9. Complexity Analysis

### 9.1 Per Read-Pair Complexity

| Phase | Time | Space |
|-------|------|-------|
| Anchor sort | O(N) expected, O(N log k) worst | O(1) in-place |
| Quick check | O(N) | O(N) for DP arrays |
| Main DP | O(N × max_iter) worst, O(N × max_skip) typical | O(N) |
| Backtrack | O(C) where C = chain length | O(C) |
| Mcopy extraction | O(mcopy_num × N) | O(N) |

Where N = number of anchors for this read pair.

### 9.2 Effective Complexity

In practice, the DP loop is the dominant cost. With `max_iter = 5000` and
`max_skip = 25`, most anchor pairs are resolved within 25 iterations of the
inner loop, making the effective complexity closer to O(N × 25) per read pair.

The quick check further reduces this: for clean read pairs with monotonic
anchors, the entire chaining is O(N).

### 9.3 Per-Read Complexity

| Phase | Time | Space |
|-------|------|-------|
| Post-filter (max_n_chain) | O(M log M) | O(M) |
| Post-filter (COV_W) | O(M × W) | O(W) |
| Post-filter (R485) | O(M² × C̄) worst case | O(M) |

Where M = total overlaps per read, W = number of windows, C̄ = average chain
length.

### 9.4 Memory Layout

The `ThreadScratchpad` uses Structure-of-Arrays (SoA) layout for the DP
arrays. This is important for cache efficiency: the inner DP loop accesses
`hitPosA`, `hitPosB`, `dpSame/dpDiff`, and `parentSame/parentDiff`
sequentially. With SoA layout, each array occupies contiguous memory, so
hardware prefetching loads useful data into cache lines.

The scratchpad is recycled across reads via `clear()` (which preserves
vector capacity), amortizing allocation cost to O(1) per read after the
first few reads establish the high-water mark.

## 10. Implementation Files

| File | Contents |
|------|----------|
| `AssemblerInvertedIndex.cpp` | All chaining functions: DP loop, quick check, mcopy, postfilter, hit-to-anchor conversion, discovery and PAF entry points |
| `AlignmentCanonicalization.hpp` | Scoring functions: `comput_sc_ch`, `comput_sc_ch_ec`, `cal_bw`, `NORMAL_W` macro |
| `hifiasmCoordinateTransforms.hpp` | `rcIntervalToForward` for RC coordinate conversion |
| `AssemblerOptions.hpp` | Parameter definitions and defaults |
| `AssemblerOptions.cpp` | Parameter registration with boost::program_options |
| `Assembler.hpp` | Data structure definitions |

### Key Functions in AssemblerInvertedIndex.cpp

| Function | Lines | Purpose |
|----------|-------|---------|
| `getHifiasmLchainDpOptions` | ~230 | Configure DP parameters for ONT/HiFi |
| `hifiasm_quick_ck_lchain` | ~1200 | O(N) quick check for collinear segments |
| `run_main_dp_loop` | ~1310 | Main O(N × max_iter) DP loop |
| `backtrack_best_chain` | ~1400 | Follow predecessor pointers |
| `hifiasm_lchain_qdp_mcopy_fast` | ~1440 | Full DP pipeline with mcopy |
| `hifiasm_lchain_qgen_mcopy_fast_postfilter` | ~1560 | Three-stage post-filter |
| `applyMcopyFastSelection` | ~750 | PAF path mcopy extraction |
| `runQuickLinearChainPrefix` | ~610 | PAF path quick check (SoA format) |

## 11. Relationship to Hifiasm

The chaining algorithm is a line-by-line port of hifiasm's ONT error
correction path. The following table maps Dinara functions to their hifiasm
counterparts:

### 11.1 Core DP Functions

| Dinara | Hifiasm | File:Line |
|--------|---------|-----------|
| `getHifiasmLchainDpOptions` | `set_lchain_dp_op` | anchor.cpp:2272 |
| `hifiasm_cal_bw` | `cal_bw` | Hash_Table.cpp:1475 |
| `hifiasm_comput_sc_ch_ec` | `comput_sc_ch_ec` | Hash_Table.cpp:1515 |
| `hifiasm_comput_sc_ch` | `comput_sc_ch` | Hash_Table.cpp:1490 |
| `hifiasm_normal_w` / `HIFIASM_NORMAL_W` | `normal_w` | Hash_Table.cpp:20 |
| `hifiasm_get_chainLen` | `get_chainLen` | Hash_Table.cpp:779 |
| `hifiasm_push_ovlp_chain_qgen` | `push_ovlp_chain_qgen` | Hash_Table.cpp:1752 |
| `hifiasm_quick_ck_lchain` | `quick_ck_lchain` | Hash_Table.cpp:2007 |
| `hifiasm_lchain_qdp_mcopy_fast` | `lchain_qdp_mcopy_fast` | Hash_Table.cpp:2097 |

### 11.2 Post-Filter Functions

| Dinara | Hifiasm | File:Line |
|--------|---------|-----------|
| `hifiasm_ha_ov_type` | `ha_ov_type` | anchor.cpp:86 |
| `hifiasm_lchain_qgen_mcopy_fast_postfilter` | postfilter in `lchain_qgen_mcopy_fast` | anchor.cpp:1953 |

### 11.3 Scoring Variants

Dinara implements two scoring functions:

- **`comput_sc_ch_ec`** (ONT): Uses split penalty strategy — `min(lin, adaptive)`
  for small gaps (dd < 4), `max(lin, adaptive)` for large gaps (dd ≥ 4). This
  is more conservative, penalizing large diagonal deviations more heavily.

- **`comput_sc_ch`** (HiFi): Always uses `min(lin, adaptive)` regardless of
  gap size. This is more lenient, appropriate for high-accuracy reads where
  large gaps are rare.

The `invertedIndexUseEcScoring` parameter (default: true) selects between them.
When true, `comput_sc_ch_ec` is used; when false, `comput_sc_ch` is used.

### 11.4 Behavioral Parity Notes

The following behaviors were verified for exact parity with hifiasm
v0.25.0-r726 (commit ec9a8b2):

1. **`keep == 0` in downsampling**: When a high-frequency streak is too short
   for downsampling (`round(span/sampleDistance) == 0`), all markers in the
   streak are discarded. Hifiasm's `select_mz_h` sets `max_high_occ = 0`,
   which means `hf_select` is never called.

2. **`maxFreq = 2000`**: The hard cutoff for k-mer frequency is fixed at 2000
   (`invertedIndexMaxKmerCount`), matching hifiasm's `max_kmer_cnt = 2000`.
   This is not scaled by coverage — even at 500× coverage, k-mers with
   count > 2000 are excluded.

3. **`chain_cutoff = 2`**: The per-strand minimum hit count and the R485 weak
   chain threshold are both 2, matching hifiasm's `ecovlp.cpp:3274`.

4. **Score normalization in mcopy**: DP scores are shifted by `plus` (minimum
   score) before threshold comparison, ensuring non-negative values for the
   packed `(score << 32 | index << 1)` representation.

5. **Tie-breaking**: Among equal-score chains, the one with shorter effective
   overlap length (`get_chainLen`) is preferred, matching hifiasm's comparison
   in the DP loop and chain candidate sorting.
