# Anchor Window Partitioning — Clean Algorithm

## Overview

`computeAnchorWindowsClean` partitions the genome's anchor set into disjoint
**windows**. Each window is defined by a contiguous interval on a backbone
read's anchor journey. The algorithm is greedy: it processes the longest
available interval first, claims anchors for that window, and re-pushes any
remaining unclaimed intervals for affected reads.

Compared to the original `computeAnchorWindows`, the clean version adds two
features:

1. **LIS-based backbone ordering.** For each read that shares anchors with the
   backbone, a Longest Increasing Subsequence (LIS) of backbone positions
   enforces consistent ordering. Only LIS-selected anchors count as shared.

2. **Alternate paths from non-direct overlaps.** Reads that are *not* direct
   cis overlaps of the backbone (i.e., not in the read graph) contribute
   alternate paths between consecutive LIS pillars. These paths form bubbles
   at heterozygous sites.

The algorithm is implemented in `Assembler::computeAnchorWindowsClean`
(`src/AssemblerAnchorWindowsClean.cpp`). Data structures are in
`src/AnchorWindows.hpp`.

## Prerequisites

- **Reads** — loaded and accessible via `reads->`.
- **Markers** — loaded (`markers`), used for base-span computation.
- **Shasta2Anchors** — the anchor set.
- **Shasta2Journeys** — per-oriented-read ordered anchor sequences.
- **Read graph** — must contain only cis overlaps (from
  `createReadGraphFromPhasingCisOverlaps`). Used to distinguish direct from
  non-direct overlaps.
- **readIdsSortedByLength** — all ReadIds sorted by raw sequence length,
  longest first.

## Data Structures

### AnchorWindowReadInterval

```cpp
struct AnchorWindowReadInterval {
    OrientedReadId orientedReadId;
    uint32_t begin;              // Inclusive journey position.
    uint32_t end;                // Exclusive journey position.
    uint32_t touchedAnchorCount; // Anchors shared with backbone (LIS count).
};
```

One read's participation in a window. `touchedAnchorCount` is the number of
anchors that survived LIS filtering — it may be less than `end - begin`.

### AnchorWindowAlternatePath

```cpp
struct AnchorWindowAlternatePath {
    Shasta2AnchorId anchorIdA;   // LIS pillar (backbone anchor) at start.
    Shasta2AnchorId anchorIdB;   // LIS pillar (backbone anchor) at end.
    vector<Shasta2AnchorId> intermediateAnchorIds;
};
```

A chain of non-backbone anchors between two consecutive LIS pillars, extracted
from a non-direct overlap read. The full path is
`anchorIdA → intermediateAnchorIds[0] → ... → anchorIdB`.

### AnchorWindow

```cpp
struct AnchorWindow {
    uint32_t windowId;
    OrientedReadId backboneOrientedReadId;
    uint32_t backboneBegin;      // Inclusive journey position on backbone.
    uint32_t backboneEnd;        // Exclusive journey position on backbone.
    uint32_t claimedAnchorCount;
    vector<AnchorWindowReadInterval> readIntervals;
    vector<AnchorWindowAlternatePath> alternatePaths;
};
```

## Algorithm

### Priority Queue and Stale-Candidate Detection

A max-heap orders candidates by:
1. **Base span** (descending) — genomic distance between first and last anchor.
2. **Read length** (descending) — raw sequence length.
3. **OrientedReadId** (ascending) — deterministic tie-breaking.

Each read has a **generation counter**. Every candidate records the generation
at push time. When a candidate is popped, if its generation doesn't match the
read's current generation, it is stale and discarded. This avoids processing
intervals that were invalidated by a later claiming pass.

### Initialization

For each read in `readIdsSortedByLength` order, push its full strand-0 journey
`[0, journeySize)` as a candidate. Only intervals with ≥ 2 anchors
(`minBackboneWindowAnchors`) are pushed.

### Main Loop

Repeat until the heap is empty:

1. **Pop** the top candidate.
2. **Generation check.** Discard if stale.
3. **Unclaimed check.** Scan every anchor in the interval. If any has been
   claimed since this candidate was pushed:
   - Increment the read's generation (invalidating older entries).
   - Scan the read's full journey for contiguous unclaimed runs.
   - Push each qualifying run as a new candidate.
   - Continue to the next candidate.
4. **Create window** (see below).

### Window Creation

Given a validated backbone interval `[seedBegin, seedEnd)` on `backboneOid`:

#### 1. Build backbone lookup

Build a hash map from anchor ID to backbone journey position for the interval
`[seedBegin, seedEnd)`. Anchors whose k-mer appears more than once in the
backbone interval are excluded from the lookup (except the first and last
anchors, which are preserved to avoid disconnecting the chain). This prevents
repeat-induced loops in the backbone chain.

#### 2. Record backbone read interval

The backbone itself is always the first entry in `readIntervals`.

#### 3. Collect deferred claim spans

All claiming is deferred to the end of window creation. A list of
`(orientedReadId, begin, end)` spans accumulates throughout the process. The
backbone span is added first.

#### 4. Find touching reads (epoch-based)

Iterate over every anchor in the backbone interval. For each anchor, look up
all oriented reads that pass through it. Track which reads have been seen using
an epoch counter (avoids clearing per-read arrays between windows).

#### 5. Build direct-overlap set

Query the read graph for all edges incident to `backboneOid`. The partners are
the **direct cis overlap** reads. This set determines whether a touching read
contributes alternate paths.

#### 6. Process each touching read

For each touching read (excluding the backbone):

**a. Find shared anchors.** First, identify anchors with duplicate k-mers in
the read's journey (excluding edge anchors). Then walk the journey: for each
anchor that is not a read duplicate, appears in the backbone lookup, and is
still unclaimed, record the read position and the corresponding backbone
position.

**b. Compute LIS.** Run patience sorting (O(n log n)) on the backbone positions
to find the longest increasing subsequence. This enforces backbone order —
anchors that would require going backward on the backbone are discarded.

**c. Record span.** The read's span is `[readPos of first LIS anchor,
readPos of last LIS anchor + 1)`. This span is added to the deferred claim
list. A `readInterval` is also recorded with `touchedAnchorCount` set to the
LIS length.

**d. Extract alternate paths (non-direct overlaps only).** If the read is
*not* a direct cis overlap and the LIS has ≥ 2 anchors: for each pair of
consecutive LIS pillars, collect the read's intermediate anchors (those between
the two pillar positions in the read's journey).

Each intermediate anchor is filtered in two ways:

1. **Duplicate k-mer filter.** Anchors whose k-mer appears more than once in
   the read's journey are skipped (same duplicate detection as step 6a).

2. **Non-cis read filter.** Only anchors that have at least one oriented read
   that is **not** the backbone and **not** a direct cis overlap of the
   backbone are kept. This removes intermediates that merely echo the backbone
   haplotype — anchors whose reads are all cis with the backbone carry the
   same haplotype and add no structural information. The surviving
   intermediates are anchors visited by reads from a different haplotype or
   structural path.

If any intermediates survive both filters, an `AnchorWindowAlternatePath` is
created and added to the window.

```
Backbone:  ... A --------- B ...     (LIS pillars)
Read:      ... A  x  y  z  B ...     (x, y, z are intermediates)

If x has a duplicate k-mer in the read → filtered out.
If y has only cis reads → filtered out.
If z has a non-cis read → kept.
                                      → alternate path: A → z → B
```

Direct cis overlaps skip this step entirely because they are expected to agree
with the backbone — their anchors between LIS pillars are redundant.

#### 6e. Deduplicate shared intermediates

After all touching reads have been processed, the same intermediate anchor may
appear in multiple alternate paths (from different reads with different LIS
pillar pairs). Each intermediate must belong to exactly one path to avoid
creating edges to multiple backbone anchors.

Paths are sorted by **pillar B backbone position**, furthest first. Each
intermediate is assigned to the first path that contains it — the one whose
pillar B reaches furthest forward on the backbone. This maximizes forward
connectivity for long-range phasing decisions. Paths that lose all their
intermediates after deduplication are discarded.

```
Path 1:  A → X → B    (B at backbone pos 8)
Path 2:  C → X → D    (D at backbone pos 15)

X is kept only in Path 2 (furthest pillar B).
Path 1 loses X and is discarded if no other intermediates remain.
```

#### 7. Claim anchors (deferred)

After all touching reads have been processed, iterate over all accumulated
spans. For each anchor in each span that is still unclaimed, claim it for this
window. Claiming also claims the reverse-complement anchor (`anchorId ^ 1`).

#### 8. Bump generations and re-push

For every touching read, increment its generation counter (invalidating any
heap entries pushed before this window) and scan its journey for remaining
contiguous unclaimed intervals. Push each qualifying interval as a new
candidate.

### Why Claiming is Deferred

If anchors were claimed during the per-read loop, a read processed early could
claim anchors that a later read also shares with the backbone. The later read
would then see those anchors as already claimed and exclude them from its LIS
computation, potentially missing valid shared anchors. Deferring all claims to
the end ensures every touching read sees the same unclaimed state during
discovery.

## Graph Construction from Anchor Windows

The `Shasta2AnchorGraph` constructor in `src/Shasta2AnchorGraph.cpp` builds the
anchor graph from the computed windows. It creates three types of edges:

### Intra-window edges

For each window, consecutive backbone anchors `(journey[pos], journey[pos+1])`
become edges. These form a linear chain — the backbone of the window.

### Alternate-path edges

For each `AnchorWindowAlternatePath`, a chain of edges is created:
`anchorIdA → intermediate[0] → ... → intermediate[N-1] → anchorIdB`. These
run parallel to the backbone chain between the same two LIS pillars, forming
bubbles at heterozygous sites.

### Inter-window edges

All read journeys are walked. For each read, the algorithm tracks which window
each anchor belongs to (via `anchorToWindow` lookup). When the window changes,
an edge is added from the last backbone anchor seen in the previous window to
the first backbone anchor in the new window.

Only **backbone anchors** are registered in `anchorToWindow`. Alternate-path
intermediate anchors are not registered and are invisible to the inter-window
walk. This is correct because any read traversing an alternate path must also
pass through the LIS pillar anchors on both sides, which are backbone anchors
with valid window assignments. The inter-window transitions are detected at
those pillars.

Within a window, the tracker only moves **forward** in backbone order. If a
read visits backbone anchors out of order (due to alignment inconsistencies),
backward positions are ignored. This prevents spurious backward edges.

Inter-window edges are deduplicated by `(anchorIdA, anchorIdB)` pair.

## Properties

- **Disjoint anchor ownership.** Each anchor belongs to at most one window.
  Reverse-complement pairs are co-claimed.

- **Greedy longest-first.** The largest base-span interval becomes the next
  backbone. This produces windows covering large genomic regions.

- **LIS consistency.** Shared anchors are filtered to respect backbone order.
  Anchors that would require going backward on the backbone are excluded.

- **Bubble formation.** Non-direct overlap reads contribute alternate paths
  between LIS pillars. Intermediate anchors are filtered to keep only those
  with non-cis reads and non-duplicate k-mers, and deduplicated across paths
  to prevent backward connections. The result is clean bubbles representing
  genuinely divergent haplotypes.

- **Duplicate k-mer exclusion.** Anchors whose k-mer appears more than once
  in a journey are excluded from the backbone lookup, shared anchor detection,
  and alternate path intermediates. Edge anchors (first/last in interval) are
  preserved to avoid disconnecting chains.

- **Reads span multiple windows.** A single read may contribute intervals to
  several windows.

- **Backbone reads can seed multiple windows.** If a backbone candidate is
  partially claimed, its remaining unclaimed sub-intervals are re-pushed.

## Test Function

`testAnchorWindowsCleanLongestRead` (`src/AssemblerTestAnchorWindowsClean.cpp`)
runs the algorithm on all reads and writes:

- **AnchorWindowsClean.gfa** — backbone chains, alternate paths, and
  inter-window connecting edges. Vertices are deduplicated. All edges carry
  `RC:i:N` tags indicating the number of common oriented reads between the
  two anchors, computed by two-pointer merge on sorted marker info spans.
- **AnchorWindowsClean.csv** — per-window HSL-based colors for Bandage
  visualization.

## Pipeline Integration

In `main.cpp`, the clean algorithm replaces the original anchor window
computation:

```cpp
vector<AnchorWindow> anchorWindows;
assembler.computeAnchorWindowsClean(
    assembler.shasta2Anchors, assembler.shasta2Journeys,
    readIdsSortedByLength, anchorWindows, threadCount);
assembler.shasta2AnchorGraph = make_shared<Shasta2AnchorGraph>(
    *shasta2Anchors, *shasta2Journeys, anchorWindows, threadCount);
```

The resulting `Shasta2AnchorGraph` feeds into transitive reduction, assembly
graph construction, and superbubble detection.

## Source Files

| File | Contents |
|------|----------|
| `src/AnchorWindows.hpp` | `AnchorWindow`, `AnchorWindowReadInterval`, `AnchorWindowAlternatePath` |
| `src/AssemblerAnchorWindowsClean.cpp` | `computeAnchorWindowsClean` implementation |
| `src/AssemblerTestAnchorWindowsClean.cpp` | Test function with GFA/CSV output |
| `src/Shasta2AnchorGraph.cpp` | Graph constructor from anchor windows |
