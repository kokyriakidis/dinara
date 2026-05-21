# Noisy-Region MSA using abPOA — Implementation Plan

Port of pgphase step 4 (`collect_noisy_vars_step4`) into Dinara's k-means
phasing pipeline.

## Overview

During initial k-means phasing, variant sites are discovered from pairwise
CIGAR alignments between the backbone read and each overlapping read.  Some
regions of the backbone are too noisy for pairwise alignment to reliably call
variants — these are the "noisy regions" detected by the sliding-window
filter in `kmParseCigars`.

This step recovers hidden variants in noisy regions by running a multiple
sequence alignment (MSA) of all reads spanning each region.  The MSA reveals
variant positions that pairwise alignment missed, and the reads are scored at
those positions.  The new variant profiles are merged into the candidate table
and the k-means clustering is re-run with the expanded variant set.

### Why abPOA

pgphase uses abPOA with `max_n_cons=2` for de-novo haplotype clustering when
no prior haplotype labels exist.  Theseus does not have built-in multi-
consensus clustering.  Using abPOA for both paths (per-haplotype and combined)
keeps the implementation uniform and faithful to the pgphase port.

### Why chained ordinals instead of edlib

pgphase uses edlib to find the overlap region between read 0 (seed) and each
subsequent read in the POA.  In Dinara, the chained alignment ordinals from
the marker-based aligner already provide exact anchor correspondences between
backbone and read.  These are more precise than edlib's approximate endpoint
detection and are already computed.

---

## Files

| File | Role |
|---|---|
| `PhasingKmeansAlign.hpp` | Types and function declarations |
| `PhasingKmeansAlign.cpp` | Read collection, abPOA seeding, subgraph alignment |
| `PhasingKmeansNoisyMsa.cpp` | Variant extraction, read scoring, cross-haplotype scoring, profile merge, outer loop |

---

## Data flow

```
kmNoisyMsaStep4 (outer loop over noisy regions)
  └─ collectNoisyRegMsa (per region)
       ├─ collectNoisyReadInfo
       │    ├─ Walk chained ordinals to find flanking markers
       │    ├─ Extract backbone seed sequence (widest envelope)
       │    └─ Extract per-read subsequences + backbone-relative positions
       │
       ├─ [has both haplotypes] Per-haplotype path
       │    ├─ abpoaMsaRun(hap1 reads, max_n_cons=1, includeBb=true)
       │    └─ abpoaMsaRun(hap2 reads, max_n_cons=1, includeBb=false)
       │
       └─ [else] Combined fallback
            └─ abpoaMsaRun(all reads, max_n_cons=2, includeBb=true)

  └─ [nCons == 2] processTwoHapResults
       ├─ extractVariantsFromMsa (hap1 MSA)
       ├─ extractVariantsFromMsa (hap2 MSA)
       ├─ Merge variant lists (union, dedup by KmVarKey)
       ├─ crossScoreReads (hap1 reads at all variants)
       ├─ crossScoreReads (hap2 reads at all variants)
       └─ Merge into scratch.candidates + overlapProfiles

  └─ [nCons == 1] processOneMsaResult
       ├─ extractVariantsFromMsa
       ├─ scoreReadsAtVariants
       └─ mergeNoisyVariants into scratch
```

---

## Step 1 — `collectNoisyReadInfo`

**Goal:** Identify all reads spanning the noisy region, extract their
sequences, and determine the backbone seed.

1. **Determine the backbone read** — always strand 0 of the query read.

2. **For each overlap, walk its chained ordinal pairs** `(bbOrd, rdOrd)`:
   - Each pair represents a shared marker — a position where we know the
     exact correspondence between backbone and read.
   - Find the **left anchor**: last pair where the backbone marker position
     ≤ `regStart`.
   - Find the **right anchor**: first pair where the backbone marker position
     ≥ `regEnd`.
   - These two pairs define the read's entry and exit points around the noisy
     region.
   - When `bbNeedsFlip` is true (backbone in slot 1, opposite strand), the
     flipped backbone ordinals decrease as j increases.  We sort by backbone
     position before searching to handle this correctly.

3. **Track the widest envelope** across all reads:
   `[globalLeftOrd, globalRightOrd]` — the leftmost left-anchor and rightmost
   right-anchor backbone ordinals.  This defines the seed region.

4. **Extract the backbone seed sequence** between `globalLeftOrd` and
   `globalRightOrd` (marker midpoints), encoded as 0123 (ACGT) for abPOA.

5. **For each read:**
   - Extract its subsequence between its own flanking marker ordinals.
   - Record the backbone-relative base positions of those markers (`bbLeftPos`,
     `bbRightPos`).  These are 0-based offsets into the backbone seed sequence.
   - Map the overlap index from `scratch.overlaps` into `info.overlapIndices`
     for downstream profile updates.

### Ordinal orientation

Alignment ordinals are stored for the canonical orientation:
- Slot `[0]` → `OrientedReadId(readIds[0], 0)` (always strand 0)
- Slot `[1]` → `OrientedReadId(readIds[1], isSameStrand ? 0 : 1)`

The read's `OrientedReadId` is constructed to match its ordinal slot:
- `queryIsRead0 == true`:  read is `readIds[1]` on strand
  `isSameStrand ? 0 : 1`
- `queryIsRead0 == false`: read is `readIds[0]` on strand 0

When the read is on the RC strand, its ordinals decrease as backbone ordinals
increase.  The `leftRdOrd > rightRdOrd` swap handles this.

---

## Step 2 — `abpoaMsaRun`

**Goal:** Build the POA graph seeded with the backbone, then align each read
to the subgraph between its shared anchor points.

### Parameters

- `info`: read collection from step 1
- `readIndices`: which reads from `info` to include (info-level indices)
- `maxNCons`: 1 for per-haplotype, 2 for de-novo clustering
- `includeBackboneInReads`: whether to include the backbone MSA row in the
  result's `readMsaRows` (true for hap1, false for hap2)

### Seeding

On an empty graph (`node_n == 2`, only SRC and SINK),
`abpoa_align_sequence_to_subgraph` returns -1 and `res.n_cigar` stays 0.
`abpoa_add_subgraph_alignment` detects the empty graph internally and calls
`abpoa_add_graph_sequence`, creating a linear chain:
`SRC(0) → node2 → node3 → ... → node(L+1) → SINK(1)`.

After seeding, backbone base position `p` (0-based) → node ID `p + 2`.
This mapping is **stable** — subsequent reads add new nodes with higher IDs
(for insertions and mismatches), but the backbone nodes keep their original
IDs.

### Subgraph alignment

For each read:
- `incBeg = bbLeftPos + 2` (inclusive start node)
- `incEnd = bbRightPos + 2 - 1` (inclusive end node; `bbRightPos` is
  exclusive, so the last included base is at `bbRightPos - 1`)
- `abpoa_subgraph_nodes` converts inclusive node IDs to exclusive boundaries
  encompassing the full subgraph (including branches from previous reads)
- `abpoa_align_sequence_to_subgraph` aligns the read to that subgraph
- `abpoa_add_subgraph_alignment` adds the alignment to the graph

### Output

`abpoa_output(ab, abpt, nullptr)` generates consensus + MSA.

MSA row layout in `abc->msa_base`:
- `[0]` = backbone (read 0)
- `[1..nReads]` = reads (in input order)
- `[n_seq + 0]` = consensus 0
- `[n_seq + 1]` = consensus 1 (if `n_cons == 2`)

MSA values: 0=A, 1=C, 2=G, 3=T, 4=gap.

### Result construction

The `KmAbpoaMsaResult` struct stores:
- `readIndices`: **overlap indices** into `scratch.overlaps` (mapped from
  info indices via `info.overlapIndices[infoIdx]`).  `-1` = backbone.
- `readMsaRows`: MSA rows for each entry in `readIndices`.
- `consensusMsaRow`: MSA row for the consensus.
- `backboneMsaRow`: MSA row for the backbone.
- `backboneStartPos`: absolute backbone base position of the MSA's left edge.

When `includeBackboneInReads` is true, the backbone is added as the first
entry in `readMsaRows` with `readIndices = -1`.

When `max_n_cons == 2` (combined fallback), the backbone's cluster assignment
determines hap orientation: the cluster containing the backbone becomes
`results[0]` (hap1), the other becomes `results[1]` (hap2).  The backbone is
only included in the hap1 cluster's `readMsaRows`.

---

## Step 3 — `collectNoisyRegMsa`

**Goal:** Decide whether to run per-haplotype or combined POA.

### Per-haplotype path (has both haplotypes)

Requires ≥ `minHapReads` reads from each haplotype.

- Run `abpoaMsaRun` for hap1 reads with `max_n_cons=1`,
  `includeBackboneInReads=true`.  The backbone is a hap1 read.
- Run `abpoaMsaRun` for hap2 reads with `max_n_cons=1`,
  `includeBackboneInReads=false`.  The backbone seeds the graph but is not
  included in hap2's read list.
- Both must succeed; if either fails, return 0 (matches pgphase
  `with_ps_hap`).
- Returns 2.

### Combined fallback (insufficient haplotype labels)

Used when too few reads from one haplotype, or all reads are `hap=0`
(unclassified from initial k-means).

- Run `abpoaMsaRun` with all reads and `max_n_cons=2`.
- abPOA's internal clustering separates reads into two groups de-novo.
- Backbone's cluster → hap1, other → hap2.
- Returns 1 or 2 depending on what abPOA produces.

This matches pgphase's `no_ps_hap` path, which handles the case where the
initial clean variant sites didn't produce any phasing — the noisy region MSA
itself discovers the haplotype separation.

---

## Step 4 — Variant extraction and scoring

### `extractVariantsFromMsa`

Walks the MSA column by column, comparing the consensus row against the
backbone row.  Discovers:
- **SNPs**: backbone has base X, consensus has base Y (both non-gap)
- **Insertions**: backbone is gap, consensus has bases
- **Deletions**: backbone has bases, consensus is gap

Each variant records:
- `backbonePos`: absolute backbone base position
- `msaColStart`, `msaColEnd`: MSA column range (for same-MSA scoring)
- `refBase`/`altBase` (SNPs) or `refSeq`/`altSeq` (indels)

### `scoreReadsAtVariants` (same-MSA scoring)

Used when all reads are in the same MSA (combined fallback, `nCons == 1`).
Scores each read at each variant's MSA column range:
- SNP: read has ref base → 0, alt base → 1, gap/other → -1
- Insertion: all gaps → 0 (ref), all bases → 1 (alt), mixed → -1
- Deletion: all bases → 0 (ref), all gaps → 1 (alt), mixed → -1

### `crossScoreReads` (cross-MSA scoring)

Used when reads are in different MSAs (per-haplotype path, `nCons == 2`).
The two per-haplotype MSAs have different column alignments, so we cannot use
`msaColStart`/`msaColEnd` across MSAs.

Instead, we build a backbone-position-to-MSA-column map (`buildBbPosToMsaCol`)
by walking the backbone MSA row and recording which column each non-gap base
corresponds to.  Then for each variant, we find the corresponding column(s)
in this MSA by backbone position and score the reads there.

This ensures hap1 reads are scored at variants discovered from hap2's MSA
and vice versa.

### `processTwoHapResults` (per-haplotype path)

1. Extract variants from both haplotypes' MSAs.
2. Merge variant lists (union, deduplicated by `KmVarKey::operator==`).
3. Determine category: variant found in only one haplotype → `NoisyCandHet`,
   found in both → `NoisyCandHom`.
4. Cross-score all reads from both MSAs at all variant positions.
5. Combine scores and readIndices from both haplotypes.
6. For each variant: create or find existing `KmCandidate`, update overlap
   profiles.

### `processOneMsaResult` (combined fallback, `nCons == 1`)

Simpler path: all reads are in the same MSA.
1. Extract variants.
2. Score reads using `scoreReadsAtVariants` (same-MSA, uses `msaColStart`).
3. Merge via `mergeNoisyVariants`.

---

## Backbone handling

The backbone read is always strand 0 of the query read.  It is a **hap1
read** — it represents the reference haplotype.

### In the POA graph

The backbone always seeds the abPOA graph as read 0.  This provides:
- A stable coordinate system (position `p` → node `p + 2`)
- A reference for subgraph alignment endpoints
- A reference row in the MSA for variant calling

### In the MSA results

The backbone is included in `readMsaRows` **only for hap1** (controlled by
`includeBackboneInReads`):
- Per-haplotype path: included in hap1's result, excluded from hap2's
- Combined fallback: included in the cluster mapped to hap1

The backbone's `readIndices` entry is `-1` (sentinel) since it has no overlap
in `scratch.overlaps` to update.

### In allele counting

The backbone's allele score counts toward `variantToCandidate` coverage
(it's a real read).  At a het site, the backbone shows the hap1 (ref) allele.

### In profile updates

Code that updates `scratch.overlapProfiles` skips entries where
`readIndices == -1`.

---

## Profile update mechanics

When a new variant is added to `scratch.candidates`, each read's overlap
profile must be extended to cover it.

- `candIdx`: index of the new candidate in `scratch.candidates`
- `off = candIdx - prof.startVarIdx`: offset into `prof.alleles`
- If `prof.startVarIdx < 0` (uninitialized), set it to `candIdx`
- Extend `prof.alleles` with `-1` (no observation) up to `off`
- Set `prof.alleles[off]` to the read's score (0=ref, 1=alt, -1=unknown)
- Update `prof.endVarIdx`

Overlaps not present in any noisy-region MSA are not extended.  The k-means
code handles this via bounds checking: `if(off >= prof.alleles.size())
return -1`.

---

## Verification checklist

All items verified against abPOA source code, pgphase reference, and
end-to-end traces with concrete examples.

### abPOA API

- [x] `abpoa_init_graph` sets `node_n = 2` (SRC=0, SINK=1)
- [x] `abpoa_add_graph_node` returns sequential IDs starting from 2
- [x] `abpoa_add_subgraph_alignment` handles empty graph (`node_n == 2`)
- [x] `abpoa_subgraph_nodes` takes inclusive node IDs as input
- [x] MSA layout: `msa_base[0..n_seq-1]` = seqs,
      `[n_seq..n_seq+n_cons-1]` = consensus
- [x] `clu_read_ids[ci][k]` = original input sequence indices (0-based)
- [x] All data copied from `abc` before `abpoa_free`
- [x] `abpoa_output(ab, abpt, nullptr)` generates data without printing

### Ordinal handling

- [x] Ordinal slot assignment matches `readOid` for all 4 combinations of
      `queryIsRead0 × isSameStrand`
- [x] `bbNeedsFlip` correct: only true when `bbSlot == 1 && !isSameStrand`
- [x] Anchor search sorts by bbPos to handle reversed ordinal order
- [x] Read ordinal swap handles RC strand correctly
- [x] `extractSegmentSeq0123` uses `getOrientedReadBase` with correct `oid`

### Coordinate mapping

- [x] Backbone and read extraction both use `position + kHalf` midpoint
- [x] Both extract `[start, end)` — exclusive end
- [x] `relLeft ≥ 0`, `relRight ≤ bbLen` guaranteed by envelope construction
- [x] Node ID: `incBeg = bbLeftPos + 2`, `incEnd = bbRightPos + 2 - 1`
- [x] Edge case: full-span read produces correct node range

### Backbone handling

- [x] Backbone included in hap1's `readMsaRows` only, not hap2's
- [x] Backbone cluster assignment orients hap1/hap2 in combined path
- [x] Backbone allele counts in `variantToCandidate`
- [x] Backbone skipped in profile updates (`readIndices == -1`)

### Variant extraction and scoring

- [x] `extractVariantsFromMsa` tracks backbone position correctly
- [x] SNP, insertion, deletion all handled
- [x] `scoreReadsAtVariants` scores all reads at MSA column ranges
- [x] `crossScoreReads` uses backbone-position-to-MSA-column mapping
- [x] `crossScoreReads` always appends exactly `nReads` entries per variant
- [x] `processTwoHapResults` merges variants from both haplotypes
- [x] Category: het if in one haplotype, hom if in both
- [x] `varToKey` sets all `KmVarKey` fields correctly (refLen, altLen)
- [x] Duplicate detection uses `KmVarKey::operator==` (all 5 fields)
- [x] Profile extension relative to `startVarIdx`, not absolute

### End-to-end trace

- [x] Per-haplotype path: het SNP discovered from hap2's MSA, hap1 reads
      cross-scored as ref, hap2 reads scored as alt
- [x] Combined fallback: backbone clusters with hap1, correct orientation
- [x] All 4 `queryIsRead0 × isSameStrand` ordinal cases traced

---

## Bugs found and fixed during development

| # | Bug | Fix |
|---|---|---|
| 1 | `bbNeedsFlip` caused incorrect anchor search | Sort ordinal pairs by backbone position |
| 2 | `backboneStartPos` passed as `regStart` | Added `backboneStartPos` to `KmAbpoaMsaResult` |
| 3 | `readIndices` were info indices, not overlap indices | Map through `info.overlapIndices` |
| 4 | Profile extension absolute, not relative to `startVarIdx` | Use `off = candIdx - startVarIdx` |
| 5 | Backbone included in hap2's read list | Added `includeBackboneInReads` parameter |
| 6 | Backbone excluded from allele counting | Include in `variantToCandidate` |
| 7 | Hap1 reads not scored at hap2's variants | Added `crossScoreReads` + `processTwoHapResults` |
| 8 | `crossScoreReads` could skip variants | Always append `nReads` entries |
| 9 | `KmVarKey.refLen` never set | Added `varToKey` helper |
| 10 | Duplicate detection incomplete | Use `KmVarKey::operator==` |
| 11 | Hap2 reads not profiled for shared variants | Always update profiles, even for existing candidates |
