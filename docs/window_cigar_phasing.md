# Window CIGAR-Based Phasing

Per-window haplotype phasing using CIGAR-derived variant profiles.
Implemented in `src/AssemblerWindowCigarMSA.cpp`.

## Overview

Each anchor window has a backbone read and a set of overlapping reads with
CIGAR-based variant profiles (SNPs detected from pairwise alignments against
the backbone). The phasing pipeline finds heterozygous SNP sites, chains
them via DP, and classifies reads into haplotype clusters using transitive
closure.

The pipeline runs in two rounds. Round 0 splits reads relative to the
original backbone. Round 1 re-splits the trans group relative to the
longest-spanning trans read, using new het sites.

## Architecture

```
cigarDetectSnpsInWindow(window, anchors, journeys)
    │
    ├── Build variant profiles from OverlapCigarStore CIGARs
    │     ├── Direct overlaps: parse CIGAR vs backbone
    │     ├── Transitive overlaps: project through shared anchors
    │     └── Anchor-based fallback: recruit via anchor sharing
    │
    ├── Detect het SNPs (outer function)
    │     ├── Sweep-line per-position coverage (deletion-aware)
    │     ├── AF filter, strand bias (Fisher exact), homopolymer/repeat
    │     └── Track strand bias / homopolymer rejected positions
    │
    ├── Build per-read allele evidence at passing sites
    │     ├── allele=0 (ref), allele=1 (alt), allele=-1 (missing)
    │     ├── allele=-2 (different alt — multi-allelic)
    │     └── Backbone added as virtual evidence (allele=0 at all sites)
    │
    ├── cwPhaseSplit (reusable lambda)
    │     ├── Step 1: Re-detect het SNPs within subset
    │     ├── Step 2: Build allele matrix
    │     ├── Step 3: DP chaining (gen_rphase_dp0_single_path port)
    │     ├── Step 4: Path extraction + ONT branch scoring
    │     ├── Step 5: Adjacent-site filter
    │     └── Step 6: Transitive closure (generate_haplotypes_naive_HiFi port)
    │
    ├── Round 0: cwPhaseSplit on all reads
    │     ├── Backbone as virtual evidence
    │     ├── Full filters (AF, coverage, strand bias, homopolymer)
    │     └── Result: cis + trans + unclassified
    │
    └── Round 1: cwPhaseSplit on trans + unclassified
          ├── Longest-spanning trans read as virtual backbone
          ├── Exclude round 0 valid site positions
          ├── Re-evaluate strand bias / homopolymer positions
          └── Result: cis + trans + unclassified
```

## Coverage Parameters

```
coverageHet  = kmerDistributionInfo.coverageHet   (haploid k-mer peak)
coverageHom  = kmerDistributionInfo.coverageHom   (homozygous k-mer coverage)
hetCov       = coverageHet                        (expected haplotype read count)
cc           = hetCov * 0.7, min 6                (min ref support for DP sites)
```

## Split Guard

A group is only split if it has at least diploid-level coverage:

```
nReads >= coverageHom / 1.5
```

Below this threshold, the group is already at or near single-haplotype
coverage and splitting would be noise.

## cwPhaseSplit

Reusable lambda that runs the full hifiasm-style phasing on a read subset.

### Parameters

| Parameter | Description |
|-----------|-------------|
| `subsetIdx` | Indices into `profiles[]` |
| `bbCovBegin/End` | Backbone coverage range for virtual evidence (0,0 to disable) |
| `excludePos` | Backbone positions to skip (used in prior rounds) |
| `result` | Output: per-read classification (1=cis, 2=trans, 0=unclassified) |
| `usedPositions` | Output: backbone positions of valid sites used in this split |

### Step 1: Het SNP Detection

Re-detects het SNPs within the subset. Builds a sweep-line coverage map
from the subset's `bbCovBegin/bbCovEnd` ranges, then filters SNPs by:

- Min alt support (`cwMinSnpAltSupport`)
- Min ref support (`cwMinSnpRefSupport`)
- AF range (`opts.minAf` to `opts.maxAf`)
- Strand bias (Fisher exact test, `opts.strandBiasPval`)
- Homopolymer / repeat context

Positions in `excludePos` are skipped entirely.

### Step 2: Allele Matrix

For each passing site, each read gets an allele value:

| Value | Meaning |
|-------|---------|
| 0 | Ref (no variant at this position) |
| 1 | Alt (matches the passing SNP's alt base) |
| -1 | Missing (read doesn't span this position) |
| -2 | Different alt (multi-allelic — excluded from ref and alt counts) |

If `bbCovBegin > 0`, a virtual backbone entry is added at index `nReads`
with allele=0 at all sites within its coverage range.

### Step 3: DP Chaining

Port of hifiasm's `gen_rphase_dp0_single_path`. Finds the longest chain
of mutually consistent het sites.

Two sites are consistent (`comput_sc_rphase`) if shared reads agree:
need ≥1 shared read with ref at both sites AND ≥1 shared read with alt
at both sites. Any disagreement (ref at one, alt at other) → score
`INT64_MIN`. Both-ambiguous with valid overlap → treated as ref.

### Step 4: Path Extraction + ONT Branch Scoring

Extracts the DP chain path and validates each site:

- **Dense check**: if all sites within 8bp of each other, apply strict
  consecutive-position run rejection (runs of ≥3 consecutive positions
  are rejected).
- **Per-site scoring**: each site needs `siteRefCov >= cc` to be retained.
  Sites failing this check get `siteScore = -1`.

### Step 5: Adjacent-Site Filter

Sites at distance 1 from each other are dropped. Applied after DP
(in the transitive closure phase, not before).

### Step 6: Transitive Closure

Port of hifiasm's `generate_haplotypes_naive_HiFi`. Seed-and-propagate
classification:

1. **Seed trans**: reads sorted by nAlt descending. Reads with alt at
   valid sites where `siteRefCov >= 2` are seeded as trans.
2. **De-promote**: non-trans reads with alt at valid sites cause those
   sites to get `siteScore = -1`.
3. **Propagate**: remaining reads with alt at still-valid sites → trans.
4. **Remaining**: reads with `nObs > 0` that aren't trans → cis.
5. **Unclassified**: reads with `nObs == 0` stay at result=0.

## Two-Round Splitting

### Round 0

All reads, backbone as virtual evidence, full filters.

- **Cis** → final cluster (same haplotype as backbone)
- **Trans** → passed to round 1
- **Unclassified** → passed to round 1

### Round 1

Trans + unclassified from round 0. Uses the longest-spanning trans read's
coverage range as virtual backbone evidence.

**Site exclusion**: round 0's valid site positions (`r0usedPos`) are
excluded. These sites already separated cis from trans — re-using them
would cause the same split again.

**Strand bias / homopolymer re-evaluation**: positions that failed strand
bias or homopolymer filters in the outer function are NOT excluded. They
are re-evaluated with full filters on the trans subset, since the subset
has different read composition and strand distribution.

- **Cis** → final cluster
- **Trans** → final cluster
- **Unclassified** → separate cluster (marked `clusterIsUnclassified`)

### Cluster Layout

| Scenario | Clusters |
|----------|----------|
| No split | `[all]` |
| Round 0 only | `[r0 cis] [r0 trans] [r0 unc]` |
| Round 0 + round 1 | `[r0 cis] [r1 cis] [r1 trans] [r1 unc]` |

Unclassified clusters contain reads that don't span any het site across
both rounds. They are compatible with all haplotype clusters.

## Multi-Allelic Handling

Multi-allelic positions are decomposed: each qualifying alt allele at the
same backbone position becomes a separate het site with its own
`PassingSnp` entry. For example, position 100 with both C→G and C→T
passing filters produces two independent sites.

For each site, reads are classified relative to that site's specific alt:

```cpp
if (hasThisAlt) allele = 1;       // matches this site's alt
else if (hasOtherAlt) allele = -2; // carries a different alt — excluded
else allele = 0;                   // ref
```

Reads with `allele = -2` are excluded from both ref and alt counts for
that site. They are not missing — they have a base at the position, but
it's a different variant.

## Output

### AnchorWindow Fields

| Field | Description |
|-------|-------------|
| `readClusters` | Vector of clusters, each a vector of `OrientedReadId` |
| `clusterIsUnclassified` | `true` if the cluster contains unclassified reads |
| `readHaplotypes` | Per-read cluster index assignment |

### WindowClusters.csv

Written by `main.cpp`:

```
Cluster,OrientedReadId
0,read1-0
0,read2-1
1,read3-0
...
```

## Source Files

| File | Contents |
|------|----------|
| `src/AssemblerWindowCigarMSA.cpp` | `cigarDetectSnpsInWindow`, `cwPhaseSplit`, two-round caller |
| `src/AnchorWindows.hpp` | `readClusters`, `clusterIsUnclassified`, `readHaplotypes` |
| `src/KmerDistributionInfo.hpp` | `coverageHet`, `coverageHom` definitions |
