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
no prior haplotype labels exist.  Using abPOA for both paths (per-haplotype
and combined) keeps the implementation uniform and faithful to the pgphase
port.

### Why chained ordinals instead of edlib

pgphase uses edlib to find the overlap region between read 0 (seed) and each
subsequent read in the POA.  In Dinara, the chained alignment ordinals from
the marker-based aligner already provide exact anchor correspondences between
backbone and read.  These are more precise than edlib's approximate endpoint
detection and are already computed.

### Why MSA-column scoring instead of AlnStr

pgphase builds intermediate `AlnStr` (pairwise alignment string) objects
between consensus and each read, then walks those to score variants.  This
requires tracking a cumulative `delta_ref_alt` offset to map between ref and
consensus coordinate systems as upstream indels shift positions.

Dinara scores directly from the MSA matrix.  Each variant has exact
`msaColStart`/`msaColEnd` column indices, and each read's MSA row can be
indexed directly.  This eliminates the fragile cumulative offset tracking
and the intermediate AlnStr construction, while producing identical results.

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
kmNoisyMsaStep4 (iterative outer loop)
  ├─ sortNoisyRegs (by label asc, length asc)
  └─ while(true)
       ├─ for each undone region (in sorted order):
       │    └─ collectNoisyVars1
       │         ├─ collectNoisyRegMsa (per region)
       │         │    ├─ collectNoisyReadInfo
       │         │    │    ├─ Walk chained ordinals to find flanking markers
       │         │    │    ├─ Extract backbone seed sequence (widest envelope)
       │         │    │    └─ Extract per-read subsequences + backbone-relative positions
       │         │    │
       │         │    ├─ [has both haplotypes] Per-haplotype path
       │         │    │    ├─ abpoaMsaRun(hap1 reads, max_n_cons=1, includeBb=true)
       │         │    │    └─ abpoaMsaRun(hap2 reads, max_n_cons=1, includeBb=false)
       │         │    │
       │         │    └─ [else] Combined fallback
       │         │         └─ abpoaMsaRun(all reads, max_n_cons=2, includeBb=true)
       │         │
       │         ├─ [nCons == 2] processTwoHapResults
       │         │    ├─ extractVariantsFromMsa (hap0 MSA)
       │         │    ├─ extractVariantsFromMsa (hap1 MSA)
       │         │    ├─ Sorted merge → unified variant list with varFromCons bitmask
       │         │    ├─ scoreReadsAtVariants (same-hap, full scoring)
       │         │    ├─ crossCoverageCheck (cross-hap, coverage-only ref)
       │         │    └─ mergeNoisyVariants into scratch
       │         │
       │         └─ [nCons == 1] processOneMsaResult
       │              ├─ extractVariantsFromMsa
       │              ├─ scoreReadsAtVariants
       │              └─ mergeNoisyVariants into scratch
       │
       ├─ if any new variants: kmRunKmeans (re-phase with expanded variant set)
       └─ if no region made progress: break
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
- **SNPs**: backbone has base X, consensus has base Y (both non-gap).
  Suppressed if the next column has a gap in either row (complex indel
  boundary — pgphase SNP adjacency filter).
- **Insertions**: backbone is gap, consensus has bases.  Consecutive
  insertion columns are grouped into a single variant.
- **Deletions**: backbone has bases, consensus is gap.  Consecutive
  deletion columns are grouped into a single variant.

Each variant records:
- `backbonePos`: absolute backbone base position
- `msaColStart`, `msaColEnd`: MSA column range (for same-MSA scoring)
- `refBase`/`altBase` (SNPs) or `refSeq`/`altSeq` (indels)

### `readFullyCoversVariant` (partial-cover filter)

Before scoring a read at an indel variant, checks that the read spans the
variant region by finding the nearest backbone-base columns flanking the
variant's MSA column range and verifying the read has non-gap bases there.
Mirrors pgphase's `full_cover` (`cover_start && cover_end`) check.

Reads that don't fully cover the variant get score -1 (unknown).

### `scoreReadsAtVariants` (same-MSA scoring)

Scores each read at each variant's MSA column range.  Used for same-haplotype
scoring (reads in the same MSA as the consensus that discovered the variant).

**SNPs:**
- Read has alt base → 1 (alt)
- Read has ref base → 0 (ref)
- Read has gap or other → -1 (unknown)

**Insertions:**
- Counts `nMatch` (read base == consensus base), `nMismatch` (read has
  different base), `nGaps` across the insertion columns.
- All gaps → 0 (ref: read doesn't have the insertion)
- No bases at all → -1 (unknown)
- Long insertions (≥ `kLongIndelLen` = 10bp): `nMatch ≥ 90%` of span → 1
  (alt); else if fully covered → 0 (ref); else -1
- Short insertions: exact match (`nMatch == span`, `nMismatch == 0`) → 1
  (alt); else if fully covered → 0 (ref); else -1
- The 90% similarity threshold (`kConsSimilarityThreshold`) matches pgphase's
  `cons_sim_thres`.

**Deletions:**
- All bases → 0 (ref: read has the reference sequence)
- All gaps → 1 (alt: read has the deletion)
- Mixed → -1 (unknown)

### `crossCoverageCheck` (cross-haplotype scoring)

Used in the two-haplotype path for variants from the OTHER consensus.
pgphase only checks coverage and assigns ref (0) for cross-haplotype
variants — it never does full allele scoring across haplotypes.

Uses `buildBbPosToMsaCol` to map the variant's backbone position to the
corresponding MSA column in this haplotype's MSA.

- **Deletions**: finds the backbone-base columns spanning the deletion
  length, then checks `readFullyCoversVariant`.  If covered → 0 (ref).
- **SNPs/Insertions**: checks if the read has a non-gap base at the anchor
  column.  If so → 0 (ref).
- If not covered → -1 (unknown).

### `processTwoHapResults` (two-haplotype path)

Port of pgphase `update_cand_var_profile_from_cons_aln_str2` +
`update_cand_var_profile_from_cons_aln_str21`.

1. Extract variants from both haplotypes' MSAs independently.
2. Sorted merge into a unified variant list with `varFromCons` bitmask:
   - 1 = from hap0 consensus only → `NoisyCandHet`
   - 2 = from hap1 consensus only → `NoisyCandHet`
   - 3 = from both consensuses → `NoisyCandHom`
3. Full scoring (`scoreReadsAtVariants`) for same-haplotype reads.
4. Coverage-only ref (`crossCoverageCheck`) for cross-haplotype reads.
5. Skip variants with zero alt-allele reads.
6. Merge via `mergeNoisyVariants`.

### `processOneMsaResult` (single-consensus path)

Port of pgphase `update_cand_var_profile_from_cons_aln_str1`.

All reads are in the same MSA.  All variants are `NoisyCandHom`.
1. Extract variants.
2. Score reads using `scoreReadsAtVariants`.
3. Skip variants with zero alt-allele reads.
4. Merge via `mergeNoisyVariants`.

---

## Sorted merge — `mergeNoisyVariants`

Port of pgphase `merge_var_profile` + `merge_read_var_profile_entries`.

Merges new variants into `scratch.candidates` and rebuilds all overlap
profiles with remapped indices.

### Algorithm

1. Sort new variants by `KmVarKey`.
2. Two-pointer merge with existing `scratch.candidates` (already sorted).
3. Build `oldToMerged[i]` and `newToMerged[i]` index maps.
4. **Collision** (same `KmVarKey`): keep old candidate, `newToMerged[i]`
   stays -1 (new variant's scores are discarded).
5. Rebuild overlap profiles:
   - Transfer old profile entries using `oldToMerged`.
   - Apply new variant scores using `newToMerged` (skipped for collisions).
6. Remap `scratch.validVarIdx` through `oldToMerged`.

---

## Iterative outer loop — `kmNoisyMsaStep4`

Port of pgphase `collect_noisy_vars_step4`.

### Region sorting

`sortNoisyRegs` sorts regions by (label ascending, length ascending) using
pgphase's bubble sort.  Label represents the number of existing variant sites
in the region — regions with fewer known variants are processed first.

### Retry loop

```
while(true):
    for each undone region (in sorted order):
        ret = collectNoisyVars1(region)
        if ret >= 0:  mark done, track if new vars found
        if ret < 0:   leave undone (MSA failed, may retry)
    if any new vars: kmRunKmeans (re-phase)
    if no region made progress: break
```

Regions that fail MSA (`ret == -1`) are left undone and retried on the next
pass.  The rationale: another region's new variants may change the k-means
phasing, which changes the haplotype labels, which may allow the failed
region to succeed on retry (e.g., by providing enough reads per haplotype
for the per-haplotype path).

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

## Constants

| Constant | Value | Purpose |
|---|---|---|
| `kGap` | 4 | abPOA gap value in MSA rows |
| `kLongIndelLen` | 10 | Threshold for long vs short indel scoring |
| `kConsSimilarityThreshold` | 0.9 | Minimum consensus-matching fraction for long insertions |

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
- [x] SNP adjacency filter suppresses SNPs at complex indel boundaries
- [x] SNP, insertion, deletion all handled
- [x] `readFullyCoversVariant` checks flanking backbone-base columns
- [x] Insertion scoring uses `nMatch` (consensus-matching), not any-base
- [x] Long insertion similarity threshold (90%) matches pgphase
- [x] Short insertion requires exact match
- [x] `scoreReadsAtVariants` scores all reads at MSA column ranges
- [x] `crossCoverageCheck` is coverage-only ref (matches pgphase cross-hap)
- [x] `processTwoHapResults` merges variants with `varFromCons` bitmask
- [x] Category: het if in one haplotype, hom if in both
- [x] `varToKey` sets all `KmVarKey` fields correctly (refLen, altLen)
- [x] Duplicate detection uses `KmVarKey::operator==` (all 5 fields)
- [x] Profile extension relative to `startVarIdx`, not absolute

### Sorted merge

- [x] `mergeNoisyVariants` uses two-pointer sorted merge
- [x] Collision keeps old candidate, discards new scores
- [x] `oldToMerged`/`newToMerged` index maps rebuild profiles correctly
- [x] `validVarIdx` remapped through `oldToMerged`

### Outer loop

- [x] `sortNoisyRegs` sorts by (label asc, length asc) — bubble sort
- [x] Iterative retry: failed regions (`ret == -1`) retried on next pass
- [x] `kmRunKmeans` called after each pass with new variants
- [x] Loop terminates when no region makes progress

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
| 7 | Hap1 reads not scored at hap2's variants | Added cross-haplotype scoring |
| 8 | `crossScoreReads` could skip variants | Always append `nReads` entries |
| 9 | `KmVarKey.refLen` never set | Added `varToKey` helper |
| 10 | Duplicate detection incomplete | Use `KmVarKey::operator==` |
| 11 | Dead code: `scoreReadsAtVariants` called then discarded | Removed |
| 12 | Missing `<array>` include | Added |
| 13 | Unused `bbRow` variable | Removed |
| 14 | Stale comment on `readIndices` | Updated |
| 15 | Single-pass outer loop (no retry) | Added iterative `while(true)` with done tracking |
| 16 | No region sorting | Added `sortNoisyRegs` by (label, length) |
| 17 | No inter-region k-means re-run | Added `kmRunKmeans` call after each pass |
| 18 | Candidates not sorted after merge | Rewrote as sorted merge with index maps |
| 19 | SNP adjacency filter missing | Added next-column gap check |
| 20 | No similarity threshold for long indels | Added 90% threshold for insertions ≥10bp |
| 21 | No partial-cover filtering | Added `readFullyCoversVariant` |
| 22 | Cross-haplotype scoring too aggressive | Replaced with `crossCoverageCheck` (coverage-only ref) |
| 23 | Collision merge applied new scores | `newToMerged` stays -1 on collision |
| 24 | Insertion scoring counted any base, not consensus-matching | Changed to `nMatch` (read == consensus) |
| 25 | Insertion mixed base/gap scored as ref instead of unknown | Added `nBases + nGaps == span` check |
