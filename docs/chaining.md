# Overlap Chaining in Dinara

This document describes how Dinara chains k-mer anchors into overlap
candidates. The chaining algorithm is a faithful port of hifiasm's ONT error
correction path (ecovlp.cpp, anchor.cpp, Hash_Table.cpp) adapted to Dinara's
inverted index data structures.

## Overview

After candidate finding produces a set of k-mer hits between read pairs, the
chaining phase determines which hits form coherent, collinear chains that
represent true genomic overlaps. The pipeline has four stages:

1. **Anchor preparation** — convert hits to hifiasm-format anchors, sorted by
   query position
2. **DP chaining** — O(N × max_iter) dynamic programming to find optimal
   anchor chains
3. **Multi-copy extraction** — extract up to `mcopy_num` high-scoring chains
   per read pair
4. **Post-filtering** — apply coverage-aware filters to control output volume

## Stage 1: Anchor Preparation

Each k-mer hit between reads A and B is converted to a **HifiasmKmerHit**
anchor:

```
struct HifiasmKmerHit {
    uint32_t self_offset;  // position in query read (read A)
    uint32_t offset;       // position in target read (read B)
    uint8_t  strand;       // 0 = same-strand, 1 = opposite-strand
    uint32_t cnt;          // packed: lower 8 bits = weight, upper 24 bits = occurrence count
};
```

Anchors are split into two groups by strand orientation:

- **Same-strand anchors**: both reads observed the k-mer in the same canonical
  orientation. In a true overlap, target positions increase monotonically with
  query positions.
- **Opposite-strand anchors**: reads observed opposite orientations. In a true
  overlap, target positions decrease as query positions increase.

Within each strand group, anchors are sorted by `(self_offset, offset)` — query
position first, then target position for ties. This sort order is required by
the DP algorithm, which assumes anchors are ordered by query coordinate.

The two strand groups are concatenated into a single array with a strand
boundary marker, matching hifiasm's convention of processing both orientations
in one DP pass.

## Stage 2: DP Chaining

The core chaining algorithm finds the highest-scoring chain of collinear
anchors. It is a direct port of hifiasm's `lchain_qdp_mcopy_fast`
(Hash_Table.cpp:2097).

### DP Parameters

The DP parameters are set by `getHifiasmLchainDpOptions`, which implements
hifiasm's `set_lchain_dp_op` (anchor.cpp:2272). For ONT reads
(`isAccurate = true`):

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `max_skip` | 25 | Max consecutive predecessors to skip before pruning |
| `max_iter` | 5000 | Lookback window size (max anchors to consider as predecessors) |
| `max_dis` | 5000 | Max gap distance between consecutive anchors |
| `quick_check` | 1 | Enable O(N) fast path for collinear hits |
| `chn_pen_gap` | 0.5 × exp(-0.01 × k) | Gap penalty coefficient |
| `chn_pen_skip` | 0.0005 × exp(-0.01 × k) | Skip penalty coefficient |

The penalties decay exponentially with k-mer length. Longer k-mers produce
fewer, more widely-spaced anchors, so the gap penalty is reduced to avoid
over-penalizing natural spacing.

### Quick Check (O(N) Fast Path)

Before running the full O(N × max_iter) DP, a linear-time quick check
attempts to find a high-quality chain by scanning anchors left-to-right and
greedily extending collinear runs.

The quick check (`hifiasm_quick_ck_lchain`) works as follows:

1. Initialize with the first anchor as a chain of length 1
2. For each subsequent anchor *i*, check if it extends the current chain:
   - Query position must strictly increase (`dq > 0`)
   - Target position must strictly increase for same-strand, decrease for
     opposite-strand (`dr > 0` after strand adjustment)
   - The diagonal drift `|dq - dr|` must be within bandwidth
3. If the anchor extends the chain, compute its score using `comput_sc_ch_ec`
   and add it
4. If not, start a new chain segment
5. Track the best chain seen so far

If the quick check finds a chain that spans most of the overlap region, the
full DP is skipped. This provides a significant speedup for high-quality reads
with dense, collinear anchors.

### Scoring Function

The scoring function `hifiasm_comput_sc_ch_ec` (porting hifiasm's
`comput_sc_ch_ec`, Hash_Table.cpp:1515) evaluates the quality of connecting
anchor *i* to predecessor anchor *j*:

**Step 1 — Monotonicity check:**
```
dq = ai.self_offset - aj.self_offset    // query gap
dr = ai.offset - aj.offset              // target gap (strand-adjusted)
```
Both `dq` and `dr` must be positive (anchors must advance on both axes).
If either is non-positive, return `INT32_MIN` (invalid connection).

**Step 2 — Bandwidth check:**
```
dd = |dq - dr|                          // diagonal deviation
```
If `dd > 16` and `dd > cal_bw(ai, aj, bw_rate)`, the connection deviates too
far from the expected diagonal — return `INT32_MIN`.

The bandwidth `cal_bw` is computed as:
```
alignable_span = min(query_span, target_span)
bandwidth = alignable_span × bw_rate
```
where the spans are extended to the read boundaries to account for the full
potential overlap region. With `bw_rate = 0.05` (ONT), a 10 kb overlap allows
up to 500 bp of diagonal drift.

**Step 3 — Base score:**
```
dg = min(dq, dr)                        // gap length
score = normal_w(dg, 1) = dg            // base score = gap length
```
The base score equals the number of bases covered by the gap between anchors.
The `normal_w(x, y)` function returns `x/y` if `x >= y`, else `1`.

**Step 4 — Occurrence normalization:**
```
score = normal_w(score, occurrence_weight)
```
The score is divided by the k-mer occurrence weight. High-frequency k-mers
contribute less per anchor, preventing repetitive regions from dominating
chain scores.

**Step 5 — Gap penalty:**
```
log_gap = log2(dg)
penalty = (log_gap × chn_pen_gap) + (log_gap × chn_pen_skip × n_skipped)
```
The penalty grows logarithmically with gap size and linearly with the number
of skipped anchors. The `chn_pen_skip` term discourages chains that skip many
intermediate anchors, preferring dense chains.

**Step 6 — Adaptive penalty (bandwidth-aware):**
```
if dd > 0:
    adaptive_penalty = score × (dd / dg) / bw_rate
    penalty = max(penalty, adaptive_penalty)
```
If the connection drifts off-diagonal, an additional penalty proportional to
the drift ratio is applied. This penalizes chains that approach the bandwidth
limit.

**Final score:**
```
chain_score = score - penalty
```

### DP Loop

The main DP loop processes anchors in query-position order:

```
for i = 0 to N-1:
    // Restrict lookback window
    st = max(si, i - max_iter)

    // Skip anchors on different strand
    while a[st].strand != a[i].strand: st++

    max_j = -1          // best predecessor
    max_f = a[i].cnt    // base score (anchor's own weight)
    n_skip = 0          // consecutive skips without improvement

    for j = i-1 down to st:
        // Skip different-strand anchors
        if a[j].strand != a[i].strand: continue

        sc = comput_sc_ch_ec(a[i], a[j], ...)
        if sc == INT32_MIN: continue

        // Candidate score = predecessor's chain score + connection score
        candidate = f[j] + sc
        if candidate > max_f:
            max_f = candidate
            max_j = j
            n_skip = 0
        elif p[j] >= 0:
            // j's predecessor was already used — count as skip
            n_skip++
            if n_skip > max_skip: break
    
    f[i] = max_f        // best chain score ending at i
    p[i] = max_j        // predecessor pointer
```

The `max_skip` pruning is the key optimization: once 25 consecutive
predecessors fail to improve the score, the inner loop terminates. This
prevents O(N²) behavior in practice — most anchors find their best
predecessor within a few iterations.

### Backtracking

After the DP loop, the best chain is extracted by following predecessor
pointers from the anchor with the highest score:

```
// Find global best
msc_i = argmax(f[i])
msc = f[msc_i]

// Backtrack
k = msc_i
while k >= 0:
    t[cL++] = k
    k = p[k]

// Reverse to get chain in query order
reverse(t[0..cL-1])
```

## Stage 3: Multi-Copy Extraction

For read pairs with multiple valid overlaps (e.g., reads spanning a tandem
duplication), the algorithm extracts up to `mcopy_num = 3` chains.

After extracting the best chain, the algorithm:

1. **Marks** all anchors in the best chain as used (via the `ii[]` visit array)
2. **Re-runs** the DP on the remaining unmarked anchors
3. **Accepts** the new chain if:
   - Its score ≥ `mcopy_rate × best_score` (default: 0.70 × best)
   - The total anchor count ≥ `mcopy_khit_cut` (default: 32)
4. **Repeats** until `mcopy_num` chains are extracted or no qualifying chain
   remains

Each extracted chain is converted to an overlap region via
`hifiasm_push_ovlp_chain_qgen`, which projects the chain's anchor coordinates
to read-level overlap boundaries:

```
// Left boundary: extend to read start
if query_start <= target_start:
    x_pos_s = 0
    y_pos_s = target_start - query_start
else:
    y_pos_s = 0
    x_pos_s = query_start - target_start

// Right boundary: extend to read end
query_remaining = query_length - query_end - 1
target_remaining = target_length - target_end - 1
if query_remaining <= target_remaining:
    x_pos_e = query_length - 1
    y_pos_e = target_end + query_remaining
else:
    y_pos_e = target_length - 1
    x_pos_e = query_end + target_remaining
```

This projection estimates the full overlap extent from the chain endpoints,
assuming the overlap extends to whichever read boundary is closer.

## Stage 4: Post-Filtering

The post-filter (`hifiasm_lchain_qgen_mcopy_fast_postfilter`) applies three
stages of filtering to control the number of overlap candidates per read.
This is a direct port of hifiasm's postfilter in `lchain_qgen_mcopy_fast`
(anchor.cpp:1953-2100).

### Stage 1: max_n_chain Per-Type Cap

Overlaps are classified into four types based on how the query read's
coordinates relate to the overlap boundaries:

| Type | Condition | Meaning |
|------|-----------|---------|
| 0 | Query start > 0, query end < len-1 | Query is a suffix overlap |
| 1 | Query start = 0, query end < len-1 | Query is a prefix overlap |
| 2 | Query start = 0, query end = len-1 | Query is contained |
| 3 | Query start > 0, query end = len-1 | Query is containing |

If the total number of overlaps exceeds `max_n_chain` (default:
`max(coveragePeak × 5, 100)`):

1. Sort overlaps by score (descending)
2. Count overlaps per type
3. For each type, record the score threshold when `max_n_chain` overlaps of
   that type have been seen
4. Keep an overlap only if its score exceeds the threshold for its type

This ensures that each overlap type retains at most `max_n_chain` candidates,
preventing a single type from monopolizing the output.

### Stage 2: COV_W Window Saturation

If type-3 (containing) overlaps exceed `max_n_chain`, an additional
window-based filter is applied:

1. Divide the query read into windows of size `COV_W = 3072` bases
2. For each window, count how many type-3 overlaps cover it
3. An overlap is **rescued** if it covers at least one window with fewer than
   `max_n_chain` overlaps

This prevents regions with extreme coverage from being over-represented while
ensuring that low-coverage windows retain their overlaps.

### Stage 3: R485 Weak-Chain Suppression

The final filter removes weak chains that are likely false positives. A chain
is considered **weak** if its anchor count (`align_length`) is below
`chain_cutoff = 2`.

For each weak chain *w*:

1. Compute the weak chain's query span
2. Search for **strong chains** (anchor count ≥ `chain_cutoff`) that overlap
   the weak chain's query span by at least `OFL = 95%`
3. If a strong chain exists with:
   - Occurrence count ≥ `CH_OCC = 4` (the strong chain has enough anchors)
   - Score ≥ `CH_SC = 16` × weak chain's score (the strong chain is much
     better)
4. Then **suppress** the weak chain

This removes spurious short chains in regions where a much better chain
already exists, reducing false positive overlaps.

## Overlap Region Output

Each surviving chain produces an overlap candidate with:

| Field | Description |
|-------|-------------|
| `readIdA`, `readIdB` | The two reads |
| `isSameStrand` | Whether the overlap is same-strand or opposite-strand |
| `x_pos_s`, `x_pos_e` | Overlap boundaries on read A |
| `y_pos_s`, `y_pos_e` | Overlap boundaries on read B |
| `shared_seed` | Chain score (DP score) |
| `align_length` | Number of anchors in the chain |

These candidates are passed to the alignment computation phase, which uses
A*PA2 (a SIMD-accelerated optimal edit distance aligner) to compute the
actual base-level alignment and CIGAR string.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `driftRateTolerance` | 0.05 | Bandwidth rate (`bw_rate`) for ONT reads |
| `invertedIndexLchainIsAccurate` | true | Use ONT-optimized DP parameters |
| `invertedIndexMcopyNum` | 3 | Max chains per read pair |
| `invertedIndexMcopyRate` | 0.70 | Min score ratio for secondary chains |
| `invertedIndexMcopyKhitCutoff` | 32 | Min anchors to enable multi-copy |
| `invertedIndexMcopyOcvWindow` | 3072 | COV_W window size for rescue |
| `invertedIndexHighFactor` | 5.0 | max_n_chain = max(coveragePeak × 5, 100) |
| `invertedIndexMinNChain` | 100 | Minimum max_n_chain |

## Complexity

| Phase | Time | Space |
|-------|------|-------|
| Anchor sort | O(N log N) per read pair | O(N) |
| Quick check | O(N) | O(N) |
| DP chaining | O(N × max_iter) worst case, O(N × max_skip) typical | O(N) |
| Mcopy extraction | O(mcopy_num × N) | O(N) |
| Post-filtering | O(M log M) where M = total overlaps per read | O(M) |

In practice, the DP loop is the dominant cost. With `max_iter = 5000` and
`max_skip = 25`, most anchor pairs are resolved within 25 iterations of the
inner loop, making the effective complexity closer to O(N × 25) per read pair.

## Implementation Files

| File | Contents |
|------|----------|
| `AssemblerInvertedIndex.cpp` | All chaining functions, postfilter, hit-to-anchor conversion |
| `AssemblerOptions.hpp` | Parameter definitions and defaults |
| `Assembler.hpp` | Data structure definitions |

## Relationship to Hifiasm

The chaining algorithm is a line-by-line port of hifiasm's ONT error
correction path. The following table maps Dinara functions to their hifiasm
counterparts:

| Dinara | Hifiasm | File |
|--------|---------|------|
| `getHifiasmLchainDpOptions` | `set_lchain_dp_op` | anchor.cpp:2272 |
| `hifiasm_cal_bw` | `cal_bw` | Hash_Table.cpp:1475 |
| `hifiasm_comput_sc_ch_ec` | `comput_sc_ch_ec` | Hash_Table.cpp:1515 |
| `hifiasm_normal_w` | `normal_w` | Hash_Table.cpp:20 |
| `hifiasm_get_chainLen` | `get_chainLen` | Hash_Table.cpp:779 |
| `hifiasm_push_ovlp_chain_qgen` | `push_ovlp_chain_qgen` | Hash_Table.cpp:1752 |
| `hifiasm_quick_ck_lchain` | `quick_ck_lchain` | Hash_Table.cpp:2007 |
| `hifiasm_lchain_qdp_mcopy_fast` | `lchain_qdp_mcopy_fast` | Hash_Table.cpp:2097 |
| `hifiasm_ha_ov_type` | `ha_ov_type` | anchor.cpp:86 |
| `hifiasm_lchain_qgen_mcopy_fast_postfilter` | postfilter in `lchain_qgen_mcopy_fast` | anchor.cpp:1953 |
