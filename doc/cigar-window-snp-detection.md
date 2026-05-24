# CIGAR-Based Per-Window SNP Detection

## Overview

`cigarDetectSnpsInWindow` detects heterozygous SNPs within an anchor window
by building per-read variant profiles from pairwise CIGARs against the
backbone read. It returns the count of clean het SNPs that pass all quality
filters. This count is stored in `AnchorWindow::cleanHetSnpCount` and used
to gate alternate-path output in the GFA.

**Status**: SNP counting is complete and integrated into the main pipeline.
Phasing (read clustering into haplotype groups, cis/trans classification)
is not yet implemented — the per-read profiles are built but discarded
after counting.

## Files

| File | Role |
|------|------|
| `src/AssemblerWindowCigarMSA.cpp` | All SNP detection logic (927 lines) |
| `src/AssemblerTestAnchorWindowsClean.cpp` | `writeAnchorWindowsCleanGfa` + legacy test function |
| `src/AnchorWindows.hpp` | `AnchorWindow` struct with `cleanHetSnpCount` field |
| `src/Assembler.hpp` | Declarations for `cigarDetectSnpsInWindow`, `writeAnchorWindowsCleanGfa` |
| `src/PhasingKmeansTypes.hpp` | `KmPhasingOptions`, `kmIsHomopolymer`, `kmFisherExactTwoTail` |
| `srcMain/main.cpp` | Pipeline integration (single window computation + SNP detection) |

## Pipeline Integration (main.cpp)

```
computeAnchorWindowsClean(...)        // Window planning, anchor claiming
    ↓
for each window:
    cigarDetectSnpsInWindow(...)       // Sets window.cleanHetSnpCount
    ↓
writeAnchorWindowsCleanGfa(...)       // GFA with het/hom-gated alternate paths
    ↓
Shasta2AnchorGraph(...)               // Does NOT use cleanHetSnpCount yet
```

`Shasta2AnchorGraph` and everything downstream (assembly graph, transitive
reduction, consensus) currently ignores `cleanHetSnpCount`. The het/hom
classification only affects `AnchorWindowsClean.gfa`.

## Data Structures

### CwVariant
A single variant event extracted from a CIGAR against the backbone.

```cpp
struct CwVariant {
    uint32_t bbPos;       // backbone position (oriented frame)
    KmVarType type;       // Snp, Insertion, Deletion
    uint8_t altBase;      // for SNPs: 0=A, 1=C, 2=G, 3=T
    uint16_t len;         // refLen for DEL, altLen for INS, 1 for SNP
    string insSeq;        // inserted bases (INS only)
};
```

### CwReadProfile
Per-read variant profile within a window. Contains all variants the read
has relative to the backbone, plus the backbone coordinate range the read
covers.

```cpp
struct CwReadProfile {
    OrientedReadId oid;
    vector<CwVariant> variants;   // sorted by bbPos
    bool isDirect;                // true = direct CIGAR, false = transitive/fallback
    uint32_t bbCovBegin = 0;      // backbone coverage range [begin, end)
    uint32_t bbCovEnd = 0;        // in oriented backbone coordinates
};
```

The coverage range `[bbCovBegin, bbCovEnd)` tracks which backbone positions
this read actually spans within the window. Used by the sweep-line algorithm
to compute per-position spanning counts. A read with a deletion at a position
still spans it — the coverage range is the alignment's overall footprint,
not individual base matches.

### CwPosMap
Sparse position map from intermediary read positions to backbone positions.
Built from the intermediary's CIGAR against the backbone. Used for transitive
read projection.

```cpp
struct CwPosMap {
    unordered_map<uint32_t, uint32_t> toBb;  // intermediary pos → backbone pos
};
```

Keys are in the intermediary's forward frame. Values are in the oriented
backbone frame, clamped to `[windowBbBegin, windowBbEnd)`.

## Function Reference

### Helper Functions (static, file-scope)

#### `cwCigarRead1Start(assembler, ad) → uint64_t`
Computes the read1 start position for CIGAR walking. When `isSameStrand`,
returns `ad.ts`. When different strands, returns `targetLen - ad.te` so the
CIGAR walks read1 in its RC frame.

#### `cwGetBaseOriented(assembler, readId, position, isReverseComplement) → uint8_t`
Returns a base from a read, handling RC orientation. When `isReverseComplement`,
reads from `seqLen - 1 - position` and complements.

#### `cwBuildPosMap(assembler, backboneReadId, backboneLen, backboneIsStrand1, alignmentId, windowBbBegin, windowBbEnd) → CwPosMap`
Builds a position map from an intermediary read B's CIGAR against the backbone.
Walks match/mismatch ops, mapping B's positions to backbone positions.

**Coordinate handling:**
- `bbIsRead0`: backbone is read0 → `bbPos = xk`, partner pos = `yk`
- `!bbIsRead0`: backbone is read1 → `bbPos = yk`, partner pos = `xk`
- `needsRc = !bbIsRead0 && isRev`: mirrors `yk` to forward frame
- `mirrorPartner = bbIsRead0 && isRev`: mirrors `yk` (partner) to forward frame
- `mirrorBb = backboneIsStrand1`: converts backbone pos to oriented frame
- Only positions within `[windowBbBegin, windowBbEnd)` are stored

#### `cwParseTransitiveCigar(assembler, intermediaryReadId, intermediaryLen, alignmentId, posMap, variants, covBeginOut, covEndOut)`
Parses a transitive read C's CIGAR against intermediary B, projecting
mismatches onto backbone coordinates via the posMap.

- Only extracts SNPs (op 1 = mismatch). Indels through transitive
  composition are unreliable and skipped.
- Coverage range: min/max backbone positions from posMap lookups across
  all match/mismatch ops. This iterates per-base (not just endpoints)
  because posMap is sparse.
- Output: variants appended and sorted by bbPos; covBeginOut/covEndOut set.

#### `cwParseCigarForOverlap(assembler, backboneReadId, backboneLen, backboneIsStrand1, alignmentId, windowBbBegin, windowBbEnd, variants, covBeginOut, covEndOut)`
Parses a direct overlap CIGAR against the backbone, extracting all variant
types (SNPs, insertions, deletions) within the window range.

**Coverage range tracking:**
Uses `updateCovRange` lambda that takes raw backbone start + length and
applies `needsRc` / `mirrorBb` transforms. Called for backbone-consuming ops:
- Match/mismatch (op 0, 1): always consume both reads
- Deletion (op 3 when `bbIsRead0`, op 2 when `!bbIsRead0`): consume backbone only

Coverage is clamped to `[windowBbBegin, windowBbEnd)` on output.

**Variant extraction:**
- Op 1 (mismatch): one SNP per base, filtered to window range
- Op 2/3 (ins/del): `bbConsumed = (op==3 && bbIsRead0) || (op==2 && !bbIsRead0)`
  - Backbone-consumed → deletion variant
  - Not backbone-consumed → insertion variant at anchor position

### Main Function

#### `Assembler::cigarDetectSnpsInWindow(window, anchors, journeys) → uint32_t`
Entry point. Returns the number of clean het SNPs in the window.

**Inputs:**
- `window`: `AnchorWindow` with backbone read, journey range, read intervals
- `anchors`: `Shasta2Anchors` for ordinal/position lookups
- `journeys`: `Shasta2Journeys` for journey traversal

**Early exits:**
- Window has ≤ 1 anchor: return 0
- Invalid ordinals: return 0
- `profiles.size() + 1 < cwMinReadCoverage` (6): return 0

## Three-Tier Read Recruitment

All reads come from `window.readIntervals` (index 0 is the backbone, skipped).
Each tier only processes reads not already handled by a previous tier.

### Tier 1: Direct CIGAR Overlaps (~lines 500–520)

For each read in the window, look up its `ReadId` in `partnerReadToAlnIds`
(built from the backbone's alignment table entry). If found, call
`cwParseCigarForOverlap` to extract variants and coverage range.

**Alignment selection:** `findBestAlnForWindow` picks the alignment whose
backbone coverage overlaps the window region. Falls back to the first
available chain. Keeps all chains (including secondary/deduped) via
`partnerReadToAlnIds` multi-value map.

**`partnerReadToAlnIds` construction:**
- Iterates `alnTable[bbOid.getValue()]`
- Skips alignments where both sides are deleted
- Skips alignments with no CIGAR (`cigarOffset == -1` or `cigarTokenCount == 0`)
- Keys by partner `ReadId` (not `OrientedReadId`)

### Tier 2: Transitive Projection (~lines 525–600)

For reads with no direct backbone CIGAR, find an intermediary read B that:
1. Is already in `profiles` (has a direct backbone overlap)
2. Has a pairwise CIGAR with the transitive read C (`keptByBothSides()`)

**Process:**
1. Look up C's alignment table (`alnTable[cOid.getValue()]`)
2. For each alignment partner B, check if B is in `readToProfile` (either strand)
3. Find B's backbone alignment via `findBestAlnForWindow(bReadId)`
4. Build or retrieve cached `CwPosMap` for B's backbone alignment
5. Call `cwParseTransitiveCigar` to project C's mismatches through the posMap

**PosMap caching:** keyed by `bBbAlnId` (B's backbone alignment ID). Avoids
rebuilding the same map for multiple transitive reads sharing the same
intermediary.

### Tier 3: Anchor-Based Segment Comparison (~lines 610–760)

For reads with no usable CIGAR (direct or transitive). Compares raw sequence
between consecutive shared anchors.

**Process:**
1. Build `bbAnchorToJP`: backbone anchor IDs → journey positions in window
2. For each unprofiled read C, find anchors shared with backbone (only backbone
   anchors, not intermediate anchors from other reads)
3. Require ≥ 2 shared anchors (reads with < 2 are skipped entirely)
4. Sort shared anchors by backbone journey position
5. Between consecutive anchors, compare sequence:
   - Same length segments: base-by-base SNP detection
   - Different length segments (indel): compare prefix (half of shorter
     segment from left anchor) and suffix (half from right anchor) to
     detect SNPs in non-indel flanks

**Coverage range:** `[first_anchor_pos, last_anchor_pos + k)` clamped to window.

### Why Reads Reach Tier 3

A read in the window shares anchors with the backbone (it was claimed during
window planning). But it may lack a usable CIGAR because:
- The CIGAR was never computed for that alignment pair (`cigarOffset == -1`)
- The CIGAR exists but was filtered (`!keptByBothSides()`)
- No intermediary read B has both a direct backbone CIGAR and a CIGAR with C

This is a CIGAR availability gap, not an alignment gap.

## SNP Aggregation

After building all profiles, SNPs are aggregated into a map:

```
Key:   (bbPos << 8) | altBase    (uint64_t)
Value: SnpAccum { fwd, rev, total }
```

- `fwd`: count of reads on strand 0 with this SNP
- `rev`: count of reads on strand 1 with this SNP
- `total`: fwd + rev

Only `KmVarType::Snp` variants are aggregated. Insertions and deletions
are extracted but not used for het/hom classification.

## Sweep-Line Per-Position Coverage

### Problem
Dovetail reads only cover part of the window. Using `totalReads` as the
denominator for allele frequency inflates ref support at window edges
where reads drop off, suppressing real het SNPs.

### Solution
Each `CwReadProfile` carries `[bbCovBegin, bbCovEnd)` — the backbone
positions the read actually covers within the window.

**Event construction:**
```
For each profile with bbCovBegin < bbCovEnd:
    emit (bbCovBegin, +1)
    emit (bbCovEnd,   -1)
```

**Sort order:** by position, then `+1` before `-1` at the same position.
This ensures a read starting at position p is counted as spanning p,
while a read ending at position p (half-open) is not.

**Sweep:** process events in order, maintaining a running count. For each
SNP position (sorted), advance the event pointer to include all events
at positions ≤ the SNP position. The running count is the number of reads
spanning that position.

**Complexity:** O(R log R + S log S) for R reads and S SNP positions.

### Correctness for Deletions
A read with a deletion at position p still spans p. The deletion is a
variant event (the read has information about that position — it disagrees
about length). The coverage range tracks the alignment's overall backbone
footprint, so deletions are included in the spanning count.

## SNP Classification Filters

For each SNP position × alt allele, applied in order:

| # | Filter | Threshold | Source |
|---|--------|-----------|--------|
| 1 | Min alt support | `acc.total >= 3` | `cwMinSnpAltSupport` |
| 2 | Min ref support | `refCov >= 3` (per-position) | `cwMinSnpRefSupport` |
| 3 | Allele frequency | `AF ∈ [0.20, 0.80]` (per-position) | `KmPhasingOptions::minAf/maxAf` |
| 4 | Strand bias | Fisher exact two-tail `p >= 0.01` | `KmPhasingOptions::strandBiasPval` |
| 5 | Homopolymer context | `kmIsHomopolymer(bbSeq, bbLen, vkey, 0)` | `PhasingKmeansTypes.hpp` |
| 6 | Repeat region | `kmIsRepeatRegion(bbSeq, bbLen, vkey, 0)` | No-op for SNPs (returns false) |

**Per-position calculations:**
- `spanning` = sweep-line count of reads covering this position
- `refCov = spanning - acc.total` (clamped to 0)
- `AF = acc.total / spanning`

**Strand bias test:**
- Expected: `acc.total / 2` on each strand
- Fisher exact two-tail test comparing observed (fwd, rev) vs expected
- SNP rejected if `p < 0.01`

**Homopolymer check:**
- Checks if the SNP is adjacent to a tandem repeat of unit length 1–6
  with ≥ 3 copies, looking both forward and backward from the SNP position
- Uses the full oriented backbone sequence with bounds-safe access

**Repeat region check:**
- Only handles Deletion and Insertion types
- Returns `false` for SNPs — effectively a no-op

## Coordinate Frames

All positions in the pipeline use the **oriented backbone frame** — the
backbone read as seen on its assigned strand (strand 0 or strand 1).

**CIGAR positions** are in the read0-forward frame (read0 is always strand 0
in `AlignmentData`). Two transforms convert to oriented backbone frame:

1. `needsRc`: when backbone is read1 and strands differ, `yk` is in RC frame.
   Mirror: `bbPos = backboneLen - 1 - yk`

2. `mirrorBb`: when backbone oriented read is on strand 1, the read0-forward
   frame doesn't match the oriented frame. Mirror: `bbPos = backboneLen - 1 - bbPos`

Window coordinates (`windowBbBegin`, `windowBbEnd`) are in oriented backbone
frame, derived from marker positions of the first/last anchors in the window.

## GFA Output (writeAnchorWindowsCleanGfa)

Writes `AnchorWindowsClean.gfa` and `AnchorWindowsClean.csv`.

**Vertices:** all backbone anchors from all windows (deduplicated).

**Intra-window edges:** consecutive backbone anchor pairs within each window.
Always emitted regardless of het/hom status.

**Alternate-path edges:** intermediate anchors from `window.alternatePaths`,
forming parallel chains between backbone pillar anchors. Only emitted for
windows where `cleanHetSnpCount > 0` (het windows).

**Inter-window edges:** connecting edges between windows, found by walking
all read journeys and tracking window transitions. For each ordered window
pair (A, B), keeps the latest backbone anchor in A and earliest in B.

**Edge weights:** `RC:i:N` tag with the count of common oriented reads
between the two anchor endpoints (two-pointer merge on sorted marker info).

**CSV:** Bandage coloring file. Each backbone anchor gets an HSL-derived
color based on its window ID.

## Constants

| Name | Value | Purpose |
|------|-------|---------|
| `cwMinReadCoverage` | 6 | Min profiled reads + backbone to attempt SNP detection |
| `cwMinSnpAltSupport` | 3 | Min reads supporting alt allele |
| `cwMinSnpRefSupport` | 3 | Min reads supporting ref allele (per-position) |
| `KmPhasingOptions::minAf` | 0.20 | Min allele frequency for het SNP |
| `KmPhasingOptions::maxAf` | 0.80 | Max allele frequency for het SNP |
| `KmPhasingOptions::strandBiasPval` | 0.01 | Strand bias p-value threshold |
| k | 14 (typical) | Anchor k-mer length, from `assemblerInfo->k` |

## Known Limitations

1. **No phasing.** Per-read profiles are built but discarded after counting.
   The read × SNP matrix needed for k-means clustering exists transiently
   inside `cigarDetectSnpsInWindow` but is not returned or stored.

2. **Shasta2AnchorGraph ignores het/hom.** The `cleanHetSnpCount` field is
   set but not consumed by the downstream anchor graph, assembly graph,
   or consensus pipeline.

3. **Only SNPs classified.** Insertions and deletions are extracted in
   variant profiles but not used for het/hom classification. The
   `kmIsRepeatRegion` filter is a no-op for SNPs.

4. **Transitive coverage tracking is per-base.** The posMap is sparse
   (hash map), so coverage range computation iterates every base of every
   match/mismatch op to do posMap lookups. For a 10kb alignment at 95%
   identity, this is ~9.5k hash lookups per transitive read. Could be
   optimized by tracking posMap min/max during construction.

5. **Anchor fallback only uses backbone anchors.** Reads in tier 3 are
   matched against backbone anchors only, not intermediate anchors from
   other reads' journeys. Reads with < 2 shared backbone anchors in the
   window are dropped entirely.

## Next Steps

1. **Return profiles from `cigarDetectSnpsInWindow`** — either as a return
   value or by storing them in the `AnchorWindow` struct.

2. **Build read × SNP matrix** — rows = reads, columns = het SNP positions
   (those passing all filters), values = ref (0) / alt (1) / missing (-1).

3. **K-means clustering** — split reads into 2 haplotype groups. The
   infrastructure exists in `AssemblerPhasingKmeans.cpp` (`kmRunKmeans`).

4. **Cis/trans classification** — two reads in the same cluster = cis,
   different clusters = trans. Write results to `AlignmentData` via
   `kmWriteResults`.

5. **Integrate into Shasta2AnchorGraph** — use het/hom classification
   and phasing results to control graph topology (split anchors at het
   sites, keep single anchors at hom sites).
