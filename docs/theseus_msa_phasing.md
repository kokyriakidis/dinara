# Theseus MSA-Based Phasing

Design for `phaseOverlapsMSA`, a standalone phasing function that uses
POA-based multiple sequence alignment to detect variant sites and classify
overlaps as cis/trans. Coexists with `phaseOverlapsKmeans` — both write
the same `hifiasmEcMatchState` fields so downstream code is agnostic to
which phaser ran.

## Motivation

`phaseOverlapsKmeans` discovers variants from pairwise CIGARs: each
backbone read parses its overlaps independently, counts alleles at each
position, classifies candidates, builds per-overlap allele profiles, and
runs k-means. This works well but has two limitations:

1. **Pairwise-only view.** Variant sites are detected from individual
   CIGARs against the backbone. Systematic errors shared by multiple
   reads can masquerade as het sites, and true low-frequency variants
   can be missed when a single overlap has ambiguous signal.

2. **No multi-read consensus.** Each overlap's allele call is binary
   (ref/alt) from its own CIGAR. There is no joint alignment to
   arbitrate between reads at the same locus.

A POA-based MSA aligns all reads in a local region simultaneously.
Variant sites emerge from graph topology (columns where reads disagree),
giving a multi-read consensus view that is more robust to systematic
errors and more sensitive to true variants.

## Pipeline Position

`phaseOverlapsKmeans` runs early in the pipeline (line ~1044 in
`main.cpp`), before marker graph vertices, anchors, or journeys exist.
It only needs pairwise CIGARs from `OverlapCigarStore`.

`phaseOverlapsMSA` requires anchors and journeys, so it runs later —
after `Shasta2Anchors` and `Shasta2Journeys` are constructed (line
~1181). This means it runs in a different part of the pipeline:

```
main.cpp pipeline order:
    ...
    computeBaseAlignmentsAndStore()
    computeCandidateTable()
    phaseOverlapsKmeans(threadCount)        ← CIGAR-based, runs here
    dedupChainsPrePhasing(threadCount)
    ...
    createMarkerGraphVertices(...)
    filterMarkerGraphVertices*(...)
    Shasta2Anchors(...)
    Shasta2Journeys(...)                    ← journeys now available
    phaseOverlapsMSA(threadCount)           ← MSA-based, runs here
    ...
    computeAnchorWindowsClean(...)          ← uses phasing results
    Shasta2AnchorGraph(...)
```

When testing `phaseOverlapsMSA`, `phaseOverlapsKmeans` should still
run first. The MSA phaser overwrites `hifiasmEcMatchState` for any
overlap it can classify, leaving the k-means result for overlaps not
covered by any anchor window. Alternatively, `phaseOverlapsKmeans`
can be skipped entirely if the MSA phaser proves sufficient.

The `dedupChainsPrePhasing` step currently runs right after
`phaseOverlapsKmeans`. When using the MSA phaser, it should be moved
to run after `phaseOverlapsMSA` instead, since it depends on phasing
results.

## Existing Infrastructure

### Anchor windows (`computeAnchorWindowsClean`)

Partition the genome into disjoint windows, each owned by a backbone
read. Each window covers a contiguous interval of the backbone's anchor
journey and tracks which other reads overlap it. This is the natural
unit of work for MSA-based phasing: one MSA per window (or per
consecutive anchor pair within a window).

### Theseus MSA prototype (`AssemblerTheseusReadWindowMSA.cpp`)

Already implements:

- **Anchor-interval window planning** — greedy longest-first claiming of
  anchor intervals, identical in spirit to `computeAnchorWindowsClean`.
- **Backbone-pair MSA jobs** — for each consecutive anchor pair on the
  backbone, extract the inter-anchor segment from the backbone and all
  overlapping reads, build a POA graph with Theseus, align reads.
- **Variant site detection** (`printReadWindowVariationSitesFromMsa`) —
  parse the MSA output, find "dirty" columns (mismatches/indels vs the
  backbone), group them into variant sites, count ref/alt alleles per
  read, filter by support thresholds and homopolymer context.
- **Read classification** — B-reads (span both anchors), L-reads (left
  anchor only, ends-free right), R-reads (right anchor only, ends-free
  left).
- **Parallel execution** — atomic work-queue, per-thread Theseus
  instances, no shared mutable state.

### Multi-segment MSA (`AssemblerMultiSegmentMSA.cpp`)

An alternative approach that builds a single Theseus graph spanning all
segments in a window (not just one anchor pair). Reads are aligned via
`align_from` starting at known anchor nodes. This gives a window-wide
MSA in one pass but is more expensive per window.

## Design

### Entry point

```cpp
void Assembler::phaseOverlapsMSA(uint64_t threadCount);
```

New file: `src/AssemblerPhasingMSA.cpp`.

Called from `main.cpp` after journeys are constructed:

```cpp
// After Shasta2Journeys construction:
assembler.phaseOverlapsMSA(threadCount);
```

Requires `assembler.shasta2Anchors` and `assembler.shasta2Journeys`
to already be populated. Does not construct them internally — they are
shared with `computeAnchorWindowsClean` and the anchor graph.

Both `phaseOverlapsKmeans` and `phaseOverlapsMSA` write the same
output (`hifiasmEcMatchState0/1` on `AlignmentData`). When both run,
the MSA phaser's results take precedence for overlaps it classifies.

### Pipeline overview

```
phaseOverlapsMSA(threadCount)
    │
    │   Preconditions: shasta2Anchors and shasta2Journeys already built
    │
    ├── 1. Plan anchor windows (reuse computeAnchorWindowsClean)
    │
    ├── 2. Per-window MSA and variant detection (parallel)
    │      For each window, for each consecutive anchor pair:
    │        - Extract inter-anchor segments from backbone + reads
    │        - Build Theseus POA, align all reads
    │        - Detect variant sites from MSA columns
    │        - Record per-read allele at each site
    │      Merge sites across pairs → window-wide variant table
    │
    ├── 3. Per-window phasing (parallel, same threads)
    │        - Filter variant sites (AF, coverage, homopolymer)
    │        - Build overlap allele profiles
    │        - K-means → hap1/hap2 per read
    │        - Classify overlaps: cis / trans
    │
    ├── 4. Write hifiasmEcMatchState to AlignmentData
    │
    └── 5. Optional: cis refinement (cisDifferentCopy detection)
```

### Precondition: Anchors and journeys

`phaseOverlapsMSA` expects `assembler.shasta2Anchors` and
`assembler.shasta2Journeys` to already be populated by the time it
runs. These are constructed earlier in `main.cpp`:

```cpp
assembler.shasta2Anchors = make_shared<Shasta2Anchors>(
    shasta2Owner, assembler.getReads(), assembler.assemblerInfo->k,
    *assembler.markers, assembler.markerGraph,
    threadCount, minAnchorCoverage, maxAnchorCoverage);

assembler.shasta2Journeys = make_shared<Shasta2Journeys>(
    2 * assembler.getReads().readCount(),
    shasta2Anchors, threadCount, shasta2Owner);
```

The function also requires `alignmentData` (for writing results) and
`markers` (for segment extraction).

### Step 1: Window planning

Use `computeAnchorWindowsClean` to produce `vector<AnchorWindow>`.
The `AnchorWindow` struct already contains:

- `backboneOrientedReadId`, `backboneBegin`, `backboneEnd`
- `readIntervals` — which reads overlap and their journey ranges
- `alternatePaths` — parallel chains at het sites (from the anchor
  window algorithm)

The alternate paths are not directly used for phasing but provide a
cross-check: windows with alternate paths are expected to contain
het variant sites.

### Step 2: Per-window MSA and variant detection

This is the core new work. For each window, iterate over consecutive
backbone anchor pairs and run a Theseus MSA.

#### 2a. Segment extraction

For each anchor pair `(journey[j], journey[j+1])` on the backbone:

1. Extract the backbone segment (midpoint of left anchor marker to
   midpoint of right anchor marker) using `extractMsaSegmentFromOrdinals`.
2. For each read in the window's `readIntervals`, check if it contains
   both anchors (B-read), only the left (L-read), or only the right
   (R-read). Extract the corresponding segment.
3. Cap reads per anchor pair at `maxReadsPerWindowAnchorPair` (default
   200, matching the prototype).

#### 2b. POA alignment

```cpp
theseus::Penalties penalties(0, 2, 3, 1);
theseus::Heuristics heuristics(false, false);
theseus::TheseusMSA aligner(penalties, heuristics, backboneSegment, 1, false);

for (size_t i = 1; i < sequences.size(); i++) {
    auto [leftFree, rightFree] = theseusAlignEndsFreeFlags(info[i]);
    aligner.align(sequences[i], 1, leftFree, rightFree);
}
```

B-reads are aligned first (end-to-end), then L-reads (ends-free right),
then R-reads (ends-free left). This matches the prototype's ordering
constraint.

#### 2c. Variant site detection

Reuse the logic from `printReadWindowVariationSitesFromMsa`, refactored
to return structured data instead of printing to stdout:

```cpp
struct MsaVariantSite {
    uint32_t backbonePosition;    // position on backbone read
    uint32_t anchorPairIndex;     // which anchor pair this came from
    string refAllele;
    vector<MsaAltAllele> altAlleles;  // each with sequence + read list
    vector<uint64_t> refReads;        // MSA row indices
};

struct MsaAltAllele {
    string type;       // "SNP", "INS", "DEL", "MNP"
    string sequence;
    vector<uint64_t> readRows;
};
```

Filtering (same as prototype):
- Skip sites where ref support < `rwMinSnpRefSupport` (3)
- Skip sites where max alt support < `rwMinSnpAltSupport` (3)
- Skip SNPs touching homopolymer runs >= `rwMinFilteredHomopolymerRunLength` (5)
- Skip non-SNP alleles shorter than `rwMinReportedAltLength` (2)

#### 2d. Per-read allele tracking

For each variant site, record each read's allele:

```cpp
struct ReadAlleleObservation {
    OrientedReadId readId;
    int allele;  // 0 = ref, 1 = primary alt, -1 = no observation
};
```

Across all anchor pairs in a window, accumulate a
`vector<MsaVariantSite>` and a per-read allele matrix.

### Step 3: Per-window phasing

#### 3a. AF filtering

Apply allele-fraction filtering identical to `phaseOverlapsKmeans`:
- `minAf = 0.20`, `maxAf = 0.80`
- Sites outside this range are likely homozygous or artifacts

#### 3b. Build overlap profiles

For each read that appears in the window, build an allele profile
vector across all filtered het sites. This is analogous to
`kmBuildOverlapProfiles` but uses MSA-derived alleles instead of
CIGAR-derived ones.

```cpp
struct MsaOverlapProfile {
    OrientedReadId readId;
    vector<int> alleles;  // one per filtered het site, 0/1/-1
    int firstSiteIdx;
    int lastSiteIdx;
};
```

#### 3c. K-means clustering

Reuse the k-means logic from `phaseOverlapsKmeans` (`kmRunKmeans`).
The input is the same shape: a set of overlap profiles with allele
vectors. The output is a haplotype assignment (hap1 or hap2) for each
overlap.

The k-means implementation operates on allele profiles and is
independent of how those profiles were built (CIGAR vs MSA). It can
be called directly or factored into a shared utility.

#### 3d. Overlap classification

Map haplotype assignments back to `AlignmentData`:

- hap1 (same as backbone) → `hifiasmEcMatchState = 1` (cis)
- hap2 (different from backbone) → `hifiasmEcMatchState = 2` (trans)
- unassigned → `hifiasmEcMatchState = 0` (unclassified)

### Step 4: Write results

For each overlap in the window, find its `AlignmentData` entry and
call `setHifiasmEcMatchState`. This is identical to `kmWriteResults`.

A read may appear in multiple windows. Conflict resolution:
- If a read-pair gets different classifications from different windows,
  the window where the read has more variant site observations wins.
- Ties: prefer cis over unclassified, trans over unclassified.

### Step 5: Cis refinement (optional)

After the main phasing pass, run the same cis-refinement logic as
`kmRefineCis`: within the cis set, re-run k-means to detect
`cisDifferentCopy` (state 3) overlaps from paralogous regions.

## Data flow

```
AnchorWindow
    │
    ├── backbone journey positions [begin, end)
    │
    ├── for each anchor pair (j, j+1):
    │       │
    │       ├── extract segments from backbone + reads
    │       ├── Theseus POA → MSA
    │       └── detect variant sites → MsaVariantSite[]
    │
    ├── merge sites across pairs → window variant table
    │
    ├── build per-read allele profiles
    │
    ├── AF filter → het sites only
    │
    ├── k-means → hap1/hap2 per read
    │
    └── write hifiasmEcMatchState
```

## Threading model

Same as the prototype: atomic work-queue over windows. Each thread
owns its own Theseus instance and scratchpad. No shared mutable state
except the final `AlignmentData` writes, which use per-alignment
atomics or are serialized by read ownership (each backbone read is
processed by exactly one thread).

```cpp
atomic<uint64_t> nextWindow(0);
vector<thread> threads(threadCount);
for (uint64_t t = 0; t < threadCount; t++) {
    threads[t] = thread([&]() {
        while (true) {
            uint64_t wi = nextWindow.fetch_add(1);
            if (wi >= windows.size()) break;
            processWindow(windows[wi], ...);
        }
    });
}
```

## Differences from phaseOverlapsKmeans

| Aspect | phaseOverlapsKmeans | phaseOverlapsMSA |
|--------|-------------------|-----------------|
| Pipeline position | Early (before marker graph) | Late (after journeys) |
| Prerequisites | Pairwise CIGARs only | Anchors + journeys |
| Variant source | Pairwise CIGARs | Multi-read POA MSA |
| Unit of work | One backbone read + its overlaps | One anchor window |
| Allele calls | Binary from individual CIGAR | Consensus from MSA columns |
| Sensitivity | Limited by pairwise noise | Multi-read arbitration |
| Cost | O(overlaps × CIGAR length) | O(reads × segment length) for POA |
| Noisy regions | Separate MSA step (disabled) | Naturally handled by POA |

## Reusable components from phaseOverlapsKmeans

- `kmRunKmeans` — k-means clustering on allele profiles
- `kmRefineCis` — cis-refinement for cisDifferentCopy detection
- `kmWriteResults` — writing hifiasmEcMatchState to AlignmentData
- `KmPhasingOptions` — AF thresholds, filtering parameters
- `KmCandidate` / `KmOverlapProfile` — data structures (or MSA-specific
  equivalents with the same shape)

These should be factored into shared utilities callable by both
phasing functions.

## Reusable components from the Theseus MSA prototype

- `extractMsaSegmentFromOrdinals` / `extractMsaSegmentFromBases` —
  segment extraction
- `theseusAlignEndsFreeFlags` — ends-free flag logic for B/L/R reads
- `prepareBackbonePairMsaJob` — read gathering and ordinal lookup
- `printReadWindowVariationSitesFromMsa` — variant site detection
  (refactored to return data instead of printing)
- `rwParseMsaFasta` — MSA output parsing
- Homopolymer and allele-type filtering helpers

## Output contract

Both `phaseOverlapsKmeans` and `phaseOverlapsMSA` write the same
fields on `AlignmentData`:

```cpp
alignmentData[id].hifiasmEcMatchState0  // from read0's perspective
alignmentData[id].hifiasmEcMatchState1  // from read1's perspective
```

Values: 0=unclassified, 1=cis, 2=trans, 3=cisDifferentCopy.

All downstream consumers (`removeContainedReads`, `dedupChainsPrePhasing`,
`computeAnchorWindowsClean` non-cis filter, read graph construction)
read these fields and are agnostic to which phaser produced them.

## Implementation plan

1. **Refactor shared k-means utilities** — extract `kmRunKmeans`,
   `kmRefineCis`, `kmWriteResults`, and profile-building helpers into
   a shared header/source so both phasers can use them.

2. **Refactor variant site detection** — change
   `printReadWindowVariationSitesFromMsa` to return
   `vector<MsaVariantSite>` instead of printing. Keep the print
   version as a wrapper for debugging.

3. **Implement `phaseOverlapsMSA`** — new file
   `src/AssemblerPhasingMSA.cpp`. Takes `shasta2Anchors` and
   `shasta2Journeys` as preconditions. Wire up window planning →
   MSA → variant detection → profile building → k-means → result
   writing.

4. **Wire into main.cpp** — call `phaseOverlapsMSA` after
   `Shasta2Journeys` construction. When both phasers run,
   `phaseOverlapsKmeans` runs first (early pipeline), then
   `phaseOverlapsMSA` overwrites results for overlaps it covers.
   Move `dedupChainsPrePhasing` to run after whichever phaser
   runs last.

5. **Add command-line selection** — config flag to choose between
   running only `phaseOverlapsKmeans`, only `phaseOverlapsMSA`,
   or both (k-means first, MSA overwrites).

6. **Test on chr1 15-15.4 Mbp region** — compare cis/trans/unclassified
   counts and downstream assembly quality between the two phasers.
