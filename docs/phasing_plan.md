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
        ├── 3. detectSnpSites          — single-pass SNP detection
        │       ├── Walk each overlap's CIGAR once: vote counting +
        │       │    mismatch recording + indel range tracking
        │       ├── Classify: positions with > OCC_THRES votes → candidates
        │       │    └── periodic repeat mask (period 1–4, ±12 bp window)
        │       └── Emit evidence from recorded data (no second CIGAR walk)
        ├── 4. buildSnpMatrix          — group evidence by site, filter, create PhasingSite
        │       ├── multi-alt: separate PhasingSite per qualifying alt base
        │       └── reject position if all ref matches are forward-strand
        ├── 5. runDpPhasing            — O(n²) longest compatible chain
        │       ├── pre-DP: exclude strand-biased sites (is_st_bs)
        │       ├── fill_incom: add synthetic match evidence at gap positions
        │       ├── checkCompatibility: fi/fj classification via ev.siteIdx
        │       ├── isHpcDependent: reject single-site chains needing HPC evidence
        │       ├── per-site matchCount >= cc check within multi-site chains
        │       └── chain sort tie-breaking: bufStart as secondary key
        ├── 5b. filterAdjacentSites    — remove sites at positions p and p+1
        │       (runs after DP so dpChainId is available for labelCisTrans)
        ├── 6. labelCisTrans           — greedy cis/trans labeling (ONT path)
        │       ├── Step A: count confirmed mismatches per overlap
        │       │    └── dpChainId guard: only count sites on DP chain
        │       ├── Step B: greedy labeling (sorted by count desc)
        │       │    ├── dpChainId guard on recount
        │       │    └── decrement labelMatchCount/labelFwdStrandCount for matches
        │       ├── Step C: consistency check (flip cis at transConfirmed sites)
        │       ├── Step D: set cisReset for mismatch sites at cis overlaps
        │       │    └── dpChainId guard: only process sites on DP chain
        │       ├── Step E: disabled on ONT path (multi_check not used)
        │       └── Step F: final loop (set strong, require score==1)
        │            └── dpChainId guard: only process sites on DP chain
        ├── 7. phaseLargeIndels        — ≥16 bp SV detection + DFS clustering
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
| `src/radixSort.hpp` | In-place MSD radix sort (port of hifiasm's ksort.h `KRADIX_SORT_INIT`) |
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

Port of hifiasm's `extract_sub_cigar_hc` / `hc_phase_robust_rr`
(Correct.cpp:18541, 19065).

Fused single-pass approach: each overlap's CIGAR is walked exactly once.
During the walk, mismatch positions are vote-counted AND recorded as
`(qpos, tpos)` pairs. Indel ranges are also tracked per overlap. After
the walk, candidate positions are classified and evidence is emitted
from the recorded data — no second CIGAR walk is needed.

**Single-pass CIGAR walk** — vote counting + mismatch recording:

```
flag[0..queryLen) = 0                   // full-length array
mismatchRecords = []                    // flat (qpos, tpos) pairs
indelRanges = []                        // per-overlap [start, end) ranges

for each overlap oi:
    recordStart = len(mismatchRecords)
    walk CIGAR via walkRangeWithCursor:
        if op == 1 (mismatch):
            for each base b in mismatch run:
                qpos = query-forward position
                if flag[qpos] < 255: flag[qpos]++
                record (qpos, tpos) in mismatchRecords
        if op == 2 or op == 3 (indel):
            record indel range [qpos, qpos+len) in indelRanges
    recordEnd = len(mismatchRecords)
    store (recordStart, recordEnd) for overlap oi

// Classify
for each position p in [0, queryLen):
    if flag[p] > OCC_THRES:             // ≥2 mismatch votes
        if isPeriodicRepeat(querySeq, p):
            flag[p] = 3                 // HPC-masked candidate
        else:
            flag[p] = 1                 // candidate SNP
    else:
        flag[p] = 0                     // not a candidate
```

**Evidence emission** — from recorded data (no second CIGAR walk):

```
// Alt evidence: from recorded mismatches at candidate positions
for each overlap oi:
    for each (qpos, tpos) in mismatchRecords[recordStart..recordEnd):
        if flag[qpos] == 0: skip
        emit PhasingEvidence:
            site = qpos, overlapIdx = oi, isAlt = 1
            base = target read base at tpos
            isHpc = (flag[qpos] == 3) ? 1 : 0

// Match evidence: overlaps whose query range covers a candidate
// position and has no mismatch or indel there
for each overlap oi:
    for each candidate position p in [ov.qs, ov.qe):
        if flag[p] == 0: skip
        if p falls inside an indel range for oi: skip
        if p has a recorded mismatch for oi: skip
        emit PhasingEvidence:
            site = p, overlapIdx = oi, isAlt = 0
            base = queryBases[p]
            isHpc = (flag[p] == 3) ? 1 : 0
```

Flag==3 positions still have evidence collected (with `isHpc=1`). The
`isHpc` flag is used later by `isHpcDependent` to test whether a single-site
DP chain depends on periodic-repeat evidence to pass thresholds.

#### isPeriodicRepeat

Port of hifiasm's `hpc_mask_ff` (f=NULL path, Correct.cpp:10399).

Checks if position p is in a periodic repeat by scanning the flanking
window `[p - HPC_PL, p + HPC_PL)` where `HPC_PL = 12`. Four scans check
tandem repeats of period 1–4:

```
for period in {1, 2, 3, 4}:
    // Scan 0: forward from p, check seq[j] == seq[j - period]
    extend forward while matching → ze
    extend backward while matching → zs
    clamp ze, zs to [p - 12, p + 12)
    if (ze - zs) >= period * 2: return true

    // Scan 1: backward from p
    extend backward while matching → zs
    extend forward while matching → ze
    clamp to window
    if (ze - zs) >= period * 2: return true

    // Scan 2: forward from p, check seq[j] == seq[j + period]
    (same structure, different comparison direction)

    // Scan 3: backward from p, check seq[j] == seq[j + period]
    (same structure)

return false
```

All scan boundaries are clamped to `[p - HPC_PL, p + HPC_PL)`.

### 4. buildSnpMatrix

Port of hifiasm's `radix_sort_haplotype_evdience_srt` (by site) followed
by per-site `push_info` (Correct.cpp:10511) which does radix sort by
overlapID + fused dedup/count/filter/emit.

Two-level sort matching hifiasm's pattern:

1. **Counting sort by site position** — O(n + maxSite). Site values are
   bounded by query read length (~20k), so the histogram is small and
   a single-pass counting sort outperforms MSD radix sort (which would
   do 4 passes for a uint32_t key).
   Hifiasm: `radix_sort_haplotype_evdience_srt`.

2. **Per-site `radixSort` by overlapIdx** — uses `radixSort.hpp`, a
   port of hifiasm's `KRADIX_SORT_INIT` (ksort.h). In-place MSD radix
   sort with insertion sort fallback for groups ≤ 64 (`RS_MIN_SIZE`).
   In practice all per-site groups are ≤ 64 (avg ~35), so the insertion
   sort path is always taken.
   Hifiasm: `radix_sort_haplotype_evdience_id_srt` inside `push_info`.

3. **Fused dedup + count** — the dedup loop simultaneously compacts
   duplicates and accumulates match/alt/strand counts, eliminating a
   separate counting pass. Matches hifiasm's `push_info` lines 10517-10533.

```
// Step 1: counting sort all evidence by site — O(n + maxSite)
histogram countBuf[0..maxSite+1] from evidence[].site
prefix sum → scatter into evidenceTmp[]
swap evidence ↔ evidenceTmp

// Step 2: per-site fused sort/dedup/count/filter/emit
//         (mirrors hifiasm's push_info called per site group)
for each group of evidence at the same site position:

    // (a) radixSort by overlapIdx (insertion sort for ≤ 64)
    radixSort(group.begin, group.end, [](ev) { return ev.overlapIdx; })

    // (b) fused dedup + count (hifiasm push_info lines 10517-10533)
    for each overlapIdx run (sorted, adjacent duplicates):
        keep first entry, accumulate:
            matchCount, fwdStrandRefCount, altCount[4], isHpc

    matchCount += 1                     // +1 for query read itself

    // (c) filter
    if fwdStrandRefCount == matchCount - 1: skip  // strand bias
    for each alt base with altCount[b] >= INFOR_COV and matchCount >= S_HAP_COV:
        create PhasingSite (same as before)

    // (d) emit qualifying evidence with siteIdx assigned
    //     match → lastSiteIdx, mismatch → baseSiteIdx[base]
    append to evidenceTmp (output buffer), set evidenceBegin/End

// Step 3: swap evidenceTmp → evidence (output)
```

**Performance**: 4.7× speedup over the original `std::sort` — from
~1,870ms to ~408ms (sum over 4 threads on the 989-read test dataset).
Phasing output is identical.

The `siteIdx` assignment is the key multi-alt mechanism. It allows
downstream functions to distinguish "mismatch for alt C" from "mismatch
for alt G" at the same position, matching hifiasm's `overlapSite` field.

### 5. runDpPhasing

Port of hifiasm's `gen_rphase_dp` + `gen_rphase_dp0_single_path`
(Correct.cpp:9648, 9428).

#### 5a. Pre-DP strand bias filter

Hifiasm's `gen_rphase_dp` (Correct.cpp:9665) physically removes
strand-biased sites before running the DP. We mark them as excluded
instead and skip them during scoring and chain extraction.

```
// isStrandBiased: port of hifiasm's is_st_bs macro
// st_rate = 0.05, st_max = 2
// True when: reverse-strand count <= 2 AND forward fraction >= 95%
for each site i:
    if isStrandBiased(site[i]):
        dpExcluded[i] = true
```

Excluded sites do not participate in the DP at all — they cannot be
chain endpoints, chain members, or affect `checkCompatibility` calls.

#### 5b. Incomplete evidence filling (fill_incom)

Port of hifiasm's `fill_incom` (Correct.cpp:9597, called at 9749).

If an overlap has evidence at site positions A and C but not at
intermediate position B, it likely matches the query at B. This step
adds synthetic match entries for those gaps, then recounts
`matchCount`, `fwdStrandCount`, and `altCount` from the expanded
evidence.

```
// Collect positions of non-excluded sites (deduplicated, sorted)
activeSitePositions = unique sorted positions of non-excluded sites

// Per-overlap: track last seen position index
lastPosIdx[oi] = UINT32_MAX for all overlaps

// Walk evidence (sorted by site position)
for each evidence entry ev:
    if ev.site not in activeSitePositions: skip
    curPosIdx = index of ev.site in activeSitePositions
    oi = ev.overlapIdx

    if lastPosIdx[oi] == UINT32_MAX:
        lastPosIdx[oi] = curPosIdx          // first time
    else if curPosIdx > lastPosIdx[oi] + 1:
        // Gap: add match entries for intermediate positions
        for pi = lastPosIdx[oi]+1 to curPosIdx-1:
            add synthetic match entry (oi, activeSitePositions[pi])
        lastPosIdx[oi] = curPosIdx
    else:
        lastPosIdx[oi] = curPosIdx

// Append synthetic entries, re-sort evidence, rebuild evidenceBegin/End
// Assign siteIdx = lastSiteIdx at each position (same as buildSnpMatrix)

// Recount from expanded evidence (hifiasm recounts occ_0/occ_1/overlap_num)
for each site i:
    matchCount = (count of isAlt==0 in evidence range) + 1
    fwdStrandCount = (count of forward-strand matches) + 1
    altCount = (count of isAlt==1 where siteIdx==i)
    labelMatchCount = matchCount
    labelFwdStrandCount = fwdStrandCount

// Re-evaluate strand bias with updated counts
for each site i:
    dpExcluded[i] = isStrandBiased(site[i])
```

#### 5c. checkCompatibility (comput_sc_rphase)

Tests whether two sites i and j are on the same haplotype block by
checking if overlaps covering both show consistent allele assignments.

```
checkCompatibility(siteI, siteJ):
    if sites[siteI].site == sites[siteJ].site: return false

    // Merge-join evidence ranges (both sorted by overlapIdx)
    for each overlap covering both sites:
        // Classify evidence at site I
        fi = 2                          // default: mismatch for different site
        if ev_i.isAlt == 0: fi = 0      // match
        else if ev_i.siteIdx == siteI: fi = 1  // mismatch for THIS site

        // Classify evidence at site J (same logic)
        fj = 2
        if ev_j.isAlt == 0: fj = 0
        else if ev_j.siteIdx == siteJ: fj = 1

        // Both fi=2: both are mismatches for OTHER sites at same position.
        // If both have valid siteIdx, treat as both matching (rare case).
        if fi == 2 and fj == 2 and both siteIdx != UINT32_MAX:
            fi = fj = 0

        if fi == 2 or fj == 2: return false   // one unknown → incompatible
        if fi != fj: return false              // mixed → incompatible
        nn[fi]++                               // fi==fj: count agreement

    return nn[0] > 0 and nn[1] > 0     // need evidence from both haplotypes
```

#### 5d. DP and chain extraction

```
cc = max(6, coveragePeak / 2 * 0.7)    // hifiasm: max(cut_bd, hom_cov/n_hap * cut_rate)

// DP: iterate j from i-1 DOWN to 0 (hifiasm tie-breaking)
// Excluded (strand-biased) sites are skipped entirely.
for i = 1 to n-1:
    if dpExcluded[i]: continue
    maxF = 1, maxJ = -1
    for j = i-1 down to 0:
        if dpExcluded[j]: continue
        if not checkCompatibility(i, j): continue
        sc = dpScore[j] + 1
        if sc > maxF: maxF = sc, maxJ = j
    dpScore[i] = maxF, dpParent[i] = maxJ

// Sort non-excluded sites by score descending (not just endpoints)
sort non-excluded sites by dpScore desc, index asc

// Greedy chain extraction
for each site in sorted order:
    if already claimed: skip
    follow parent pointers, collecting chain, marking claimed
    record (chainLength, chainBuffer)

// Sort chains by length descending, then by bufStart ascending
// (hifiasm tie-breaking: gen_rphase_dp0_single_path, Correct.cpp:9429)
sort chains by (length desc, bufStart asc)
```

#### 5e. Chain confirmation

```
for each chain (longest first):
    if chainLength > 1:
        plus = 1                        // multi-site → confirmed
    else:
        // Single-site: reject if HPC-dependent or low matchCount
        if isHpcDependent(site) or matchCount < cc:
            plus = -1                   // rejected
        else:
            plus = 1

    // Per-site check within confirmed chains
    for each site in chain:
        if matchCount >= cc and plus > 0:
            dpChainId = chainId         // confirmed
        else:
            dpChainId stays -1          // rejected even in confirmed chain
```

#### 5f. isHpcDependent (is_hpc_vec)

```
isHpcDependent(siteIdx):
    mc = matchCount, ac = altCount      // local copies
    for each evidence entry at this site:
        if not isHpc: skip              // only subtract HPC evidence
        if isAlt == 0: mc--             // HPC match
        else if siteIdx == siteIdx: ac-- // HPC mismatch for THIS site
    return mc < 2 or ac < 2 or mc < S_HAP_COV or ac < INFOR_COV
```

### 5g. filterAdjacentSites

Port of hifiasm's adjacent-site filter (Correct.cpp:8855). Runs after
the DP (not before) so that `dpChainId` is available for `labelCisTrans`.

Adjacent mismatches (positions p and p+1) are typically alignment artifacts
— an indel that the aligner represents as two consecutive substitutions —
rather than real heterozygous variants. When a position is removed, ALL its
PhasingSite entries are removed (multi-alt: if position p has sites for
alt C and alt G, both are removed).

Three arrays must stay consistent after removal: `sites[]` (compacted,
indices shift down), `evidence[]` (entries for removed sites dropped),
and `ev.siteIdx` (remapped from old to new site indices).

```
// 1. Identify adjacent positions
group sites by position (multi-alt entries at same position are contiguous)
for each pair of consecutive position-groups:
    if group[i+1].position == group[i].position + 1:
        mark ALL sites in both groups for removal

// 2. Build old-to-new site index mapping
for each site i:
    if not removed: oldToNew[i] = next new index
    else:           oldToNew[i] = UINT32_MAX

// 3. Compact sites array
keep only non-removed sites, shifting down to fill gaps

// 4. Compact evidence array and remap siteIdx
for each evidence entry:
    if oldToNew[ev.siteIdx] == UINT32_MAX:
        drop entry (its site was removed)
    else:
        ev.siteIdx = oldToNew[ev.siteIdx]   // remap to new index
        keep entry

// 5. Rebuild evidenceBegin/evidenceEnd on surviving sites
walk sites and evidence in tandem (both sorted by position)
for each position with surviving sites:
    find evidence range [begin, end) for this position
    set evidenceBegin/evidenceEnd on ALL sites at this position
    (multi-alt sites at same position share the same evidence range)
```

Hifiasm (Correct.cpp:8855) does the same compaction: it compacts
`snp_stat[]` and the evidence `list[]`, rewriting
`overlapSite -= m_off` where `m_off` is the cumulative count of removed
entries before each surviving group. Our approach uses an explicit
`oldToNew[]` mapping instead, which is equivalent.

### 6. labelCisTrans

Port of hifiasm's `generate_haplotypes_naive_HiFi` (Correct.cpp:8845),
ONT path (`!qual_a`).

Builds a per-overlap evidence index (`overlapEvIdx`, `overlapEvBegin/End`)
to support hifiasm's per-overlap iteration pattern while keeping the
per-site evidence layout from `buildSnpMatrix`.

All steps use a `dpChainId` guard: only sites with `dpChainId >= 0`
(confirmed by the DP in step 5) participate in mismatch counting and
labeling decisions. This matches hifiasm's ONT path where
`filterAdjacentSites` runs after the DP.

#### Threshold helper (used in Steps A, B, C, E, F)

```
passesThresholds(site):
    if site.labelMatchCount < 2 or site.altCount < 2: return false
    if isLabelStrandBiased(site): return false
    return true

passesFullThresholds(site):
    if not passesThresholds(site): return false
    if site.labelMatchCount < S_HAP_COV or site.altCount < INFOR_COV: return false
    return true
```

#### Step A — Count confirmed mismatches per overlap (Correct.cpp:8889)

```
// Build per-overlap evidence index (hifiasm sorts evidence by overlapID)
sort evidence indices by overlapIdx → overlapEvIdx[]
build overlapEvBegin[oi], overlapEvEnd[oi] ranges

// Count mismatches at DP-confirmed sites for each overlap
for each overlap oi:
    o = 0
    for each evidence entry at oi (via overlapEvIdx):
        if isAlt == 0: skip                     // only count mismatches
        site = sites[ev.siteIdx]
        if site.dpChainId < 0: skip             // dpChainId guard
        if not passesFullThresholds(site): skip
        o++
    confirmedMismatchCount[oi] = o

// Build sorted list: overlaps with o > 0, sorted by count descending
sortedOverlapIndices = [oi for oi where confirmedMismatchCount[oi] > 0]
sort sortedOverlapIndices by confirmedMismatchCount desc
```

#### Step B — Greedy labeling (Correct.cpp:8949)

This is the core greedy mechanism. Overlaps are processed from most
mismatches to fewest. As each trans overlap is processed, match counts
on its sites are decremented, potentially causing sites to fail threshold
checks for later overlaps.

On the ONT path, hifiasm does NOT do a forward walk to find additional
sites — only the sites already in the overlap's evidence are processed.

```
for each oi in sortedOverlapIndices:
    // RE-COUNT with current label counts (may have changed from decrements)
    o = 0
    for each evidence entry at oi:
        if isAlt == 0: skip
        site = sites[ev.siteIdx]
        if site.dpChainId < 0: skip             // dpChainId guard
        if not passesFullThresholds(site): skip
        o++
    if o == 0: skip                             // thresholds no longer met

    if overlaps[oi].isMatch == 1:
        overlaps[oi].isMatch = 2                // mark trans

    // Process this overlap's evidence
    for each evidence entry at oi:
        if isAlt == 1:
            sites[ev.siteIdx].transConfirmed = 1    // mark site

        else if isAlt == 0:                     // match evidence
            pos = sites[ev.siteIdx].site
            // Decrement labelMatchCount on ALL sites at this position
            // (walk backwards from ev.siteIdx — it's the last site at pos)
            for si = ev.siteIdx down to 0:
                if sites[si].site != pos: break
                sites[si].labelMatchCount--     // clamped to ≥1
                if overlap is forward-strand:
                    sites[si].labelFwdStrandCount--  // clamped to ≥1
```

#### Step C — Consistency check (Correct.cpp:8998)

Second pass over the same sorted overlap list. Catches cis overlaps that
have mismatches at sites already marked `transConfirmed` by Step B.

```
for each oi in sortedOverlapIndices:
    o = 0
    for each evidence entry at oi:
        if isAlt == 0: skip
        site = sites[ev.siteIdx]
        if site.dpChainId < 0: skip
        if not passesThresholds(site): skip     // uses decremented counts
        if site.transConfirmed: o++
    if overlaps[oi].isMatch == 1 and o > 0:
        overlaps[oi].isMatch = 2                // flip cis → trans
```

Note: `cisReset` has not been set yet at this point, so checking
`transConfirmed` alone is correct (matches hifiasm's `score == 1`).

#### Step D — Reset for cis overlaps (Correct.cpp:9015)

Marks sites that have mismatches at cis overlaps. This distinguishes
sites with mismatches exclusively at trans overlaps (remain confirmed)
from sites with mismatches at both trans and cis overlaps (reset).

```
for each overlap oi (ALL overlaps, not just sorted list):
    if overlaps[oi].isMatch != 1: skip          // only cis overlaps
    for each evidence entry at oi:
        if isAlt == 1:
            site = sites[ev.siteIdx]
            if site.dpChainId < 0: skip         // dpChainId guard
            sites[ev.siteIdx].cisReset = 1
```

After Step D, `isLabelConfirmed()` = `(transConfirmed && !cisReset)`.
Only sites with mismatches exclusively at trans overlaps remain confirmed.

#### Step E — multi_check: disabled on ONT path

Hifiasm's `multi_check` (Correct.cpp:9035) is a weak-site promotion
mechanism that catches overlaps the DP missed. On the ONT path
(`multi_check=0`), this step is skipped entirely.

`isLabelConfirmed()` = `(transConfirmed && !cisReset)` (no `promoted`
flag is set).

#### Step F — Final loop (Correct.cpp:9085)

Sets the `strong` flag and applies final cis→trans flips. On the ONT
path (no Step E promotions), `isLabelConfirmed()` is equivalent to
hifiasm's `score == 1` check: `transConfirmed && !cisReset`.

```
for each overlap oi:
    if isMatch == 2:                            // already trans
        strong = 1

    else if isMatch == 1:                       // cis
        for each evidence entry at oi:
            site = sites[ev.siteIdx]
            if site.dpChainId < 0: skip         // dpChainId guard
            if not passesThresholds(site): skip
            // Explicit score==1 check (Correct.cpp:9098):
            // transConfirmed && !cisReset
            if not site.isLabelConfirmed(): skip
            // Evidence at a confirmed site (match or mismatch)
            strong = 1
            if isAlt == 1:                      // mismatch at confirmed site
                isMatch = 2                     // flip cis → trans
                break
```

#### Mutable vs Immutable Fields

Hifiasm reuses `occ_0`, `overlap_num`, and `score` fields with different
semantics across pipeline stages. We use separate fields for clarity:

| Immutable (set by buildSnpMatrix) | Mutable (labelCisTrans) | Hifiasm equivalent |
|-----------------------------------|-------------------------|--------------------|
| `matchCount` | `labelMatchCount` | `occ_0` (decremented in Step B) |
| `fwdStrandCount` | `labelFwdStrandCount` | `overlap_num` (decremented in Step B) |

Hifiasm's `score` field is reused across steps with three different meanings.
We split it into three separate flags:

| Flag | Set in | Meaning | Hifiasm `score` equivalent |
|------|--------|---------|---------------------------|
| `transConfirmed` | Step B | Mismatch at this site in a trans overlap | `score = 1` (Step B) |
| `cisReset` | Step D | Mismatch at this site in a cis overlap | `score = -1` (Step D) |
| `promoted` | Step E | multi_check promoted this weak site (disabled on ONT) | `score = 1` (Step E) |

On the ONT path (Step E disabled), `isLabelConfirmed()` returns
`(transConfirmed && !cisReset)`, equivalent to hifiasm's `score == 1`
at the point Step F runs.

#### Strand Bias Filter

Port of hifiasm's `is_st_bs` macro (Correct.cpp:8800).

```
isStrandBiased(matchCount, fwdStrandCount):
    // Biased if reverse-strand count ≤ ST_MAX AND forward fraction ≥ (1 - ST_RATE)
    if fwdStrandCount + ST_MAX >= matchCount
       and matchCount * ST_RATE + fwdStrandCount >= matchCount:
        return true
    return false
```

Parameters: `ST_RATE = 0.05`, `ST_MAX = 2`. Two variants:
- `isStrandBiased(site)` — uses immutable `matchCount`/`fwdStrandCount`
  (pre-DP exclusion in §5a, and Steps A/B initial filtering).
- `isLabelStrandBiased(site)` — uses mutable `labelMatchCount`/`labelFwdStrandCount`
  (greedy labeling, where counts change as trans overlaps are processed).

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
| `HPC_PL` | 12 | `HPC_PL` / `hpc_flk` in `hpc_mask_ff` | Flanking window for repeat detection |
| `HPC_RR` | 4 | `hpc_rr` in `hpc_mask_ff` | Max repeat period for masking |
| `HPC_CC` | 2 | `hpc_cc` in `hpc_mask_ff` | Repeat span cutoff multiplier |
| `PHASING_ST_RATE` | 0.05 | `st_rate` | Strand bias rate threshold |
| `PHASING_ST_MAX` | 2 | `st_max` | Strand bias max reverse-strand count |

## Hifiasm Function Mapping

| Dinara | Hifiasm | Location |
|--------|---------|----------|
| `phaseOverlaps` | `rphase_hc` | Correct.cpp:20191 |
| `detectSnpSites` | `hc_phase_robust_rr` | Correct.cpp:10200 |
| `isPeriodicRepeat` | `hpc_mask_ff` (f=NULL path) | Correct.cpp:10399 |
| `isStrandBiased` / `isLabelStrandBiased` | `is_st_bs` macro | Correct.cpp:8800 |
| `buildSnpMatrix` | `radix_sort_haplotype_evdience_srt` + per-site `push_info` (with `radix_sort_haplotype_evdience_id_srt`) | Correct.cpp:10511, 19706 |
| `filterAdjacentSites` | adjacent-site filter | Correct.cpp:8855 |
| `checkCompatibility` | `comput_sc_rphase` | Correct.cpp:9244 |
| `isHpcDependent` | `is_hpc_vec` | Correct.cpp:9383 |
| `runDpPhasing` | `gen_rphase_dp` + `fill_incom` + `gen_rphase_dp0_single_path` | Correct.cpp:9648, 9597, 9428 |
| `labelCisTrans` | `generate_haplotypes_naive_HiFi` | Correct.cpp:8845 |
| `phaseLargeIndels` | `rphase_lidel` + `rphase_lidel_cc` | Correct.cpp:20155 |
| `dedupChains` | `dedup_chains` | ecovlp.cpp:2984 |

## Data Structures

Defined in `src/PhasingTypes.hpp`:

| Struct | Purpose |
|--------|---------|
| `PhasingOverlap` | Per-overlap state: coordinates, CIGAR ref, isMatch label, strong flag, error count |
| `PhasingEvidence` | One observation at a candidate SNP site from one overlap. Includes `siteIdx` (index into PhasingSite array, hifiasm's `overlapSite`) for multi-alt site tracking |
| `PhasingSite` | Confirmed het site: position, alleles, immutable counts (`matchCount`, `altCount`, `fwdStrandCount`), mutable label counts (`labelMatchCount`, `labelFwdStrandCount`), per-step flags (`transConfirmed`, `cisReset`, `promoted`), DP chain assignment |
| `PhasingSvEvent` | Contiguous indel region ≥ 16 bp from one overlap |
| `PhasingSvCluster` | Cluster of SV events at similar positions |
| `PhasingScratchpad` | Thread-local workspace, cleared between reads. Includes `evidenceTmp` and `countBuf` for global counting sort by site |

## Differences from AssemblerHifiasmEC.cpp

### Evidence source

`AssemblerHifiasmEC.cpp` reads from `AlignedEvidenceStore` (pre-extracted SNPs
and indels as separate streams). `AssemblerPhasing.cpp` reads directly from
`OverlapCigarStore`, which matches hifiasm's data flow and avoids the
intermediate evidence store.

### Alt-strand-bias filter

`performHifiasmECParity` (in `detectHetSites`, AssemblerHifiasmEC.cpp) applies
an alt-strand-bias filter during SNP detection that rejects candidate sites
where ≥95% of alt-supporting reads are on the same strand:

```cpp
if (misCount > 2) {
    if ((altReadsFwd + 2 >= misCount) && (altReadsFwd * 100 >= misCount * 95)) continue;
    if ((altReadsFwd <= 2) && (altReadsFwd * 100 <= misCount * 5)) continue;
}
```

**This filter is not present in hifiasm.** Hifiasm's `push_info`
(Correct.cpp:10511) only applies a ref-strand-bias filter
(`rev_n == occ_0`, rejecting sites where all ref overlaps are on the same
strand). The `is_st_bs` macro applied later in the DP checks ref-strand-bias
only.

`phaseOverlaps` matches hifiasm's behavior: it applies only a ref-strand-bias
filter in `buildSnpMatrix` (reject if all ref overlaps are forward-strand)
and the `is_st_bs`-equivalent check in `runDpPhasing`. No alt-strand-bias
filter is applied during detection.

**Impact**: The extra alt-strand-bias filter in `performHifiasmECParity`
causes it to reject SNP sites that hifiasm and `phaseOverlaps` would keep.
On the test dataset (read 100), this removes 11 confirmed sites, reducing
trans calls by ~43 overlaps from that read's perspective. Overall,
`performHifiasmECParity` produces ~3,200 fewer trans labels than
`phaseOverlaps` (15,345 vs 18,602 out of 38,542 alignments).


---

# K-means Overlap Phasing (`phaseOverlapsKmeans`)

Port of pgphase/longcallD's k-means phasing pipeline. Replaces the DP-based
`phaseOverlaps` above. Reads CIGARs from `OverlapCigarStore`.

## Architecture

```
computeBaseAlignmentsAndStore()     ← populates OverlapCigarStore
        │
        ▼
phaseOverlapsKmeans(threadCount)    ← entry point (AssemblerPhasingKmeans.cpp)
        │
        │  Per backbone read (thread-local KmScratchpad):
        │
        ├── 1. kmGatherOverlaps       — collect overlaps from alignmentTable
        ├── 2. kmParseCigars          — walk CIGARs, emit KmDigarOps + noisy regions
        │       ├── SNPs: vote-counted mismatches with actual alt base
        │       ├── Deletions: backbone-consumed ranges
        │       ├── Insertions: non-backbone bases with altSeq extraction
        │       └── Noisy region sliding window (pgphase XidQueue)
        ├── 3a. kmCollectCandidates   — sort + fuzzy dedup → KmCandidate list
        ├── 3b. kmAlleleCounting      — merge-walk digars against candidates
        ├── 4. kmClassifyVariants     — per-site classification + noisy region handling
        │       ├── CleanHetSnp / CleanHetIndel → used for k-means
        │       ├── RepeatHetIndel → excluded from k-means, seeds noisy regions
        │       ├── CleanHom → used for k-means
        │       └── NoisyCandHet/Hom → excluded
        ├── 5. kmRunKmeans            — 2-means clustering on overlap profiles
        └── 6. writeResults           — set matchState on AlignmentData
                                        1=cis, 2=trans (unclassified → cis)
```

## Files

| File | Role |
|------|------|
| `src/AssemblerPhasingKmeans.cpp` | Full pipeline implementation |
| `src/PhasingKmeansTypes.hpp` | Data structures, constants, thresholds |
| `src/PhasingKmeansNoisyMsa.cpp` | Noisy-region MSA step (currently disabled) |
| `src/PhasingKmeansNoisyMsa.md` | Design doc for noisy MSA |

## Variant Classification Thresholds (`KmPhasingOptions`)

| Parameter | Default | Effect |
|-----------|---------|--------|
| `minDepth` | 6 | Total coverage at site must be ≥ 6 |
| `minAltDepth` | 3 | Alt allele coverage must be ≥ 3 |
| `minAf` | 0.12 | Allele fraction must be ≥ 12% |
| `maxAf` | 1.0 | Disabled (accepts any AF above minAf) |
| `strandBiasPval` | 0.01 | (ONT only) Fisher exact test p-value threshold |
| `noisyRegMaxXgaps` | 5 | Max indel span for repeat classification |

A site passes as **CleanHet** when: `totalCov ≥ minDepth`, `altCov ≥ minAltDepth`,
`alleleFraction ≥ minAf`, no strand bias (ONT), and for indels, not in a
homopolymer/repeat region.

## Repeat Classification (indels only)

Two checks determine if a het indel is in a repetitive context, making it
unreliable for phasing. Both are gated on indel span ≤ `noisyRegMaxXgaps` (5bp).

### `kmIsHomopolymer` (pgphase `var_is_homopolymer_pg`)

Checks if the indel sits adjacent to a tandem repeat of unit length 1–6bp.
Scans both forward and backward from the indel boundary for ≥ 3 consecutive
copies of a repeat unit. Examples: 1bp deletion next to `AAAAAAA`, 2bp
insertion next to `ATATAT`.

### `kmIsRepeatRegion` (pgphase `var_is_repeat_region_pg`)

**Deletions**: compares the deleted bases `[pos, pos+delLen)` against the
next `3 × delLen` reference bases. If they match, the deletion is in a
tandem repeat.

**Insertions**: builds an expected sequence by copying the reference at
`[pos, pos+3*insLen)`, propagating a tandem repeat from the reference prefix,
then overwriting `[0, insLen)` with the actual inserted bases (from
`KmVarKey.altSeq`). If this matches the reference, the insertion is a tandem
copy of the flanking sequence.

The insertion case requires the alt sequence, which is extracted during CIGAR
parsing and stored in `KmDigarOp.altSeq` / `KmVarKey.altSeq`. For RC overlaps,
the inserted bases are complemented and reversed to put them in the backbone's
forward frame.

## Insertion Alt Sequence Extraction

During `kmParseCigars`, insertion events extract the actual inserted bases
from the non-backbone read:

- **`qIsR0 = true`** (backbone is query/read0): inserted bases come from
  read1 at positions `yk..yk+len-1`. `kmGetBase` handles RC conversion.
- **`qIsR0 = false`** (backbone is target/read1): inserted bases come from
  read0 at positions `xk..xk+len-1`. Bases are complemented if RC, then
  reversed to match the backbone's forward frame.

Stored as ACGT characters in `KmDigarOp.altSeq`, propagated to
`KmVarKey.altSeq` during candidate construction.

## Unclassified Overlaps

Overlaps not assigned to either haplotype by k-means are labeled **cis**
(`matchState = 1`), not trans. This avoids incorrectly breaking overlaps
that simply lack phasing signal.

## Noisy-Region MSA (Step 4)

Currently **disabled** (commented out in `phaseOverlapsKmeans`). The step
recovers variants in noisy regions via abPOA multiple sequence alignment.
See `src/PhasingKmeansNoisyMsa.md` for the design.
