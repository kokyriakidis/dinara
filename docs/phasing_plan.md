# ONT Overlap Phasing

Port of hifiasm v0.25.0's `rphase_hc` pipeline (Correct.cpp:20191) to Dinara.
Reads directly from `OverlapCigarStore` — no `AlignedEvidenceStore` needed.

## Architecture

```
computeBaseAlignmentsAndStore()     ← populates OverlapCigarStore
        │
        ▼
phaseOverlaps(threadCount)          ← entry point (AssemblerPhasing.cpp)
        │
        │  Per query read (thread-local PhasingScratchpad):
        │
        ├── 1. gatherOverlaps         — collect overlaps, init isMatch=1 (cis)
        ├── 2. unpackQuerySequence     — 2-bit → uint8_t for base lookups
        ├── 3. detectSnpSites          — sliding-window SNP detection (428 bp)
        │       ├── Pre-walk: single CIGAR walk per overlap → CigarEvent array
        │       ├── Pass 1: count mismatch votes per position
        │       │    └── periodic repeat mask (period 1–4, hifiasm hpc_mask_ff)
        │       └── Pass 2: emit PhasingEvidence at candidate positions
        ├── 4. buildSnpMatrix          — group evidence by site, filter, create PhasingSite
        │       └── multi-alt: separate PhasingSite per qualifying alt base
        ├── 4b. filterAdjacentSites    — remove sites at positions p and p+1
        ├── 5. runDpPhasing            — O(n²) longest compatible chain
        ├── 6. labelCisTrans           — greedy cis/trans labeling
        │       ├── Phase 1: count confirmed mismatches per overlap
        │       ├── Phase 2: greedy labeling (most mismatches → trans first)
        │       │    └── strand bias filter (skip biased sites)
        │       ├── Phase 3: consistency check (flip cis at trans-confirmed sites)
        │       └── Phase 4: multi_check (weak site promotion)
        ├── 7. phaseLargeIndels        — ≥16 bp SV detection + BFS clustering
        ├── 8. dedupChains             — best overlap per target read
        └── 9. writeResults            — per-read-perspective hifiasmEcMatchState
                │
                ▼
        hifiasmEcMatchState0/1 on AlignmentData: 1=cis, 2=trans
```

## Files

| File | Role |
|------|------|
| `src/AssemblerPhasing.cpp` | Full pipeline implementation |
| `src/PhasingTypes.hpp` | Data structures and constants |
| `src/Alignment.hpp` | `AlignmentData` with `hifiasmEcMatchState0/1` output fields |
| `src/Assembler.hpp` | `phaseOverlaps(threadCount)` declaration |

## Data Flow

### Input

- **`OverlapCigarStore`** — packed uint16_t CIGAR tokens per overlap
  (op: 0=match, 1=mismatch, 2=insertion, 3=deletion).
- **`AlignmentData`** — per-overlap metadata: `readIds[0/1]`, `isSameStrand`,
  `qs/qe/ts/te` (marker-based forward-strand coordinates),
  `cigarOffset`, `cigarTokenCount`.
- **`alignmentTable`** — `VectorOfVectors<uint32_t>` indexed by
  `OrientedReadId::getValue()`, mapping reads to alignment IDs.
- **`Reads`** — read sequences for base lookup at mismatch positions.

### Output

- **`hifiasmEcMatchState0`** — phasing label from read0's perspective (1=cis, 2=trans).
- **`hifiasmEcMatchState1`** — phasing label from read1's perspective.
- Default is 1 (cis). Each thread writes only its query read's perspective,
  so no locking is needed.
- Downstream consumers call `getCisTransStatusFromReadPerspective(queryReadId)`
  which maps 1→Cis, 2→Trans.

## Coordinate Handling

`AlignmentData.qs/ts/qe/te` store marker-based forward-strand coordinates.
The CIGAR's `yk` (read1 axis) is in oriented-read coordinates — for RC
overlaps, `oriented_start = targetLen - ad.te`. The helper `cigarRead1Start()`
handles this conversion.

When `queryIsRead0 == 0` (query is read1) and the overlap is RC, CIGAR `yk`
values are in RC coordinates. `detectSnpSites` converts these to forward:
`qpos_fwd = queryLen - 1 - qpos_rc`, then sorts the per-overlap events.

## Pipeline Steps

### 1. gatherOverlaps

Collects all non-deleted overlaps with CIGAR data for the query read.
Initializes `isMatch = 1` (cis default, matching hifiasm's post-alignment state).
Sets `queryIsRead0` based on which read in the canonical pair is the query,
and swaps qs/qe ↔ ts/te accordingly.

### 2. unpackQuerySequence

Converts the query read's 2-bit packed bases to a `uint8_t` array for
O(1) base lookups during SNP detection.

### 3. detectSnpSites

Port of hifiasm's `hc_phase_robust_rr` (Correct.cpp:10200).

**Pre-walk**: Each overlap's CIGAR is walked once, collecting match/mismatch
events as `CigarEvent` structs (qpos, tpos, op) in forward-strand query
coordinates. Events are stored contiguously with per-overlap index ranges.

**Pass 1** (per 428 bp window): Count mismatch votes per position. Positions
with > 1 vote are candidates. Candidates in periodic repeat regions
(period 1–4, span ≥ period × 2) are masked (flag=3).

**Pass 2**: At candidate positions (flag=1), emit `PhasingEvidence` entries.
Match ops record the query base; mismatch ops look up the target base from
the read sequence (with RC complement when needed).

### 4. buildSnpMatrix

Port of hifiasm's `SetSnpMatrix` + `push_info` (Correct.cpp:10511).

1. Sort evidence by (site, overlapIdx). Dedup within each (site, overlap) group.
2. Group by site. Count matchCount and per-base altCount[4].
3. Add +1 to matchCount for the query read itself.
4. For each alt base with altCount ≥ 2: require matchCount ≥ 3 (`S_HAP_COV`)
   and altCount ≥ 3 (`INFOR_COV`). Create a separate `PhasingSite` per
   qualifying alt base (multi-alt support).
5. Track `fwdStrandCount` for strand bias filtering in labelCisTrans.

### 4b. filterAdjacentSites

Port of hifiasm's adjacent-site filter (Correct.cpp:8855).
If sites exist at positions p and p+1, both are removed. Prevents alignment
artifacts (indels manifesting as adjacent mismatches) from being treated as
het sites.

### 5. runDpPhasing

Port of hifiasm's `gen_rphase_dp0_single_path` (Correct.cpp:9648).

O(n²) longest-increasing-subsequence DP over SNP sites. Two sites are
compatible (`checkCompatibility`) if overlaps covering both show consistent
allele assignments — both match or both mismatch. Mixed assignments are
skipped (not counted as incompatible).

Chains are extracted from longest to shortest (up to 15 chains). Sites in
chains of length ≥ 2 are confirmed. Single-site chains require stricter
thresholds: matchCount ≥ max(6, hetCov × 0.7).

### 6. labelCisTrans

Port of hifiasm's `generate_haplotypes_naive_HiFi` (Correct.cpp:8845).

**Phase 1**: Count confirmed mismatches per overlap at DP-confirmed,
non-strand-biased sites.

**Phase 2**: Sort overlaps by mismatch count descending. Overlaps with
confirmed mismatches → trans. Others → cis.

**Phase 3** (consistency check): Build set of "trans-confirmed" sites
(confirmed sites where a trans overlap has a mismatch). Flip any cis overlap
that mismatches at a trans-confirmed site.

**Phase 4** (multi_check): Identify weak sites (occ_0 ≥ 2, occ_1 ≥ 2, not
DP-confirmed). For each cis overlap, collect weak mismatch positions. Apply
density filter (≥ 4% of alignment length) and 32bp proximity filter (≥ 2
surviving positions). Positions appearing in ≥ 2 overlaps get promoted.
Cis overlaps with mismatches at promoted positions → trans.

#### Strand Bias Filter

Port of hifiasm's `is_st_bs` macro. A site is strand-biased if nearly all
ref-matching overlaps come from one strand. Parameters: `ST_RATE = 0.05`,
`ST_MAX = 2`. Biased sites are excluded from phase 1 counting.

### 7. phaseLargeIndels

Port of hifiasm's `rphase_lidel` (Correct.cpp:20155).

**Step 1**: Walk each cis overlap's CIGAR. Contiguous indel runs (ins/del ops)
whose span or base count ≥ 16 bp are recorded as `PhasingSvEvent`. RC
coordinates are converted to forward afterward.

**Step 2**: Sort SV events by query position.

**Step 3**: BFS connected-component clustering. Events are merged when they
overlap by ≥ 50% of the smaller span. Within each cluster, deduplicate by
target read (keep first per target).

**Step 4**: For each cluster with ≥ 3 unique targets and ≥ 3 cis overlaps
spanning the consensus region, label SV-carrying overlaps as trans.

### 8. dedupChains

Port of hifiasm's `dedup_chains` (ecovlp.cpp:2984).

Reduce to one overlap per target read. Sort by:
1. `targetReadId`
2. `isMatch` ascending (1=cis preferred over 2=trans)
3. Quality score descending: `score = (qe - qs) - 12 × errorCount`
4. Span descending (tiebreaker)

Keep first per target.

### 9. writeResults

Write `isMatch` back to `AlignmentData.hifiasmEcMatchState{0,1}` via
`setHifiasmEcMatchStateFromReadPerspective(queryReadId, isMatch)`.
Each thread writes only its query read's field, so no synchronization needed.

## Threading

Static block scheduling: reads are divided into contiguous chunks across
threads. Each thread owns a `PhasingScratchpad` (cleared between reads,
memory reused). Bounded by `threadCount > readCount` check.

## Constants

| Constant | Value | Hifiasm equivalent | Purpose |
|----------|-------|--------------------|---------|
| `PHASING_WINDOW_SIZE` | 428 | `WINDOW_MAX_SIZE` (375 + 50 + 3) | Sliding window size |
| `PHASING_OCC_THRES` | 1 | `occ_thres` | Min mismatch votes − 1 to flag position |
| `PHASING_S_HAP_COV` | 3 | `s_hap_cov` | Min ref allele count for confirmed site |
| `PHASING_INFOR_COV` | 3 | `infor_cov` | Min alt allele count for confirmed site |
| `PHASING_SV_MIN_LEN` | 16 | `min_err` in `extract_sub_cigar_sv` | Min indel size for SV phasing |
| `PHASING_MAX_DP_CHAINS` | 15 | hardcoded in `gen_rphase_dp0_single_path` | Max DP chains to extract |
| `HPC_RR` | 4 | `hpc_rr` in `hpc_mask_ff` | Max repeat period for masking |
| `HPC_CC` | 2 | `hpc_cc` in `hpc_mask_ff` | Repeat span cutoff multiplier |
| `PHASING_ST_RATE` | 0.05 | `st_rate` | Strand bias rate threshold |
| `PHASING_ST_MAX` | 2 | `st_max` | Strand bias max reverse-strand count |

## Hifiasm Function Mapping

| Dinara | Hifiasm | Location |
|--------|---------|----------|
| `phaseOverlaps` | `rphase_hc` | Correct.cpp:20191 |
| `detectSnpSites` | `hc_phase_robust_rr` | Correct.cpp:10200 |
| `isPeriodicRepeat` | `hpc_mask_ff` | Correct.cpp:10100 |
| `isStrandBiased` | `is_st_bs` macro | Correct.cpp:8800 |
| `buildSnpMatrix` | `SetSnpMatrix` + `push_info` | Correct.cpp:10511 |
| `filterAdjacentSites` | adjacent-site filter | Correct.cpp:8855 |
| `checkCompatibility` | `comput_sc_rphase` | Correct.cpp:9600 |
| `runDpPhasing` | `gen_rphase_dp0_single_path` | Correct.cpp:9648 |
| `labelCisTrans` | `generate_haplotypes_naive_HiFi` | Correct.cpp:8845 |
| `phaseLargeIndels` | `rphase_lidel` + `rphase_lidel_cc` | Correct.cpp:20155 |
| `dedupChains` | `dedup_chains` | ecovlp.cpp:2984 |

## Data Structures

Defined in `src/PhasingTypes.hpp`:

| Struct | Purpose |
|--------|---------|
| `PhasingOverlap` | Per-overlap state: coordinates, CIGAR ref, isMatch label, error count |
| `PhasingEvidence` | One observation at a candidate SNP site from one overlap |
| `PhasingSite` | Confirmed het site: position, alleles, counts, DP chain assignment |
| `PhasingSvEvent` | Contiguous indel region ≥ 16 bp from one overlap |
| `PhasingSvCluster` | Cluster of SV events at similar positions |
| `CigarEvent` | Single match/mismatch event from CIGAR pre-walk |
| `OverlapEventRange` | Index range into CigarEvent array for one overlap |
| `PhasingScratchpad` | Thread-local workspace, cleared between reads |

## Differences from AssemblerHifiasmEC.cpp

`AssemblerHifiasmEC.cpp` reads from `AlignedEvidenceStore` (pre-extracted SNPs
and indels as separate streams). `AssemblerPhasing.cpp` reads directly from
`OverlapCigarStore`, which matches hifiasm's data flow and avoids the
intermediate evidence store.
