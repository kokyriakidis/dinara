# Anchor Window Partitioning

## Overview

Anchor window partitioning divides the genome into disjoint regions called
**windows**, each defined by a contiguous interval on a backbone read's anchor
journey. Every anchor in the genome is assigned to exactly one window (or
remains unclaimed). The windows, together with the list of reads that
participate in each one, provide the substrate for per-window phasing.

The algorithm is implemented in `Assembler::computeAnchorWindows`
(`src/AssemblerAnchorWindows.cpp`). The output data structures are defined in
`src/AnchorWindows.hpp`.

## Prerequisites

The following must be available before calling `computeAnchorWindows`:

- **Reads** — loaded and accessible via `reads->`.
- **Markers** — loaded (`markers`), used to convert journey positions to base
  positions for base-span computation.
- **Shasta2Anchors** — the anchor set. Each anchor knows which oriented reads
  pass through it and at which journey position.
- **Shasta2Journeys** — for each oriented read, the ordered sequence of anchor
  IDs it passes through.
- **readIdsSortedByLength** — a vector of all ReadIds sorted by raw sequence
  length (longest first). This controls the initial heap seeding order.

## Data Structures

### AnchorWindowReadInterval

```cpp
struct AnchorWindowReadInterval {
    OrientedReadId orientedReadId;
    uint32_t begin;             // Inclusive journey position.
    uint32_t end;               // Exclusive journey position.
    uint32_t touchedAnchorCount; // Anchors shared with the backbone.
};
```

Represents one read's participation in a window. The interval `[begin, end)` is
in journey coordinates (indices into the read's anchor journey). The
`touchedAnchorCount` records how many of the backbone's anchors this read
actually shares — it may be less than `end - begin` if the read has additional
anchors between the shared ones.

### AnchorWindow

```cpp
struct AnchorWindow {
    uint32_t windowId;
    OrientedReadId backboneOrientedReadId;
    uint32_t backboneBegin;     // Inclusive journey position on backbone.
    uint32_t backboneEnd;       // Exclusive journey position on backbone.
    uint32_t claimedAnchorCount;
    vector<AnchorWindowReadInterval> readIntervals;
};
```

Represents one window. The backbone interval `[backboneBegin, backboneEnd)` is
the seed that defined this window. `readIntervals` contains entries for all
reads that participate, including the backbone itself (always the first entry).
`claimedAnchorCount` is the total number of anchors owned by this window across
all reads.

## Algorithm

### Step 1: Initialization

A max-priority-queue (heap) is seeded with one candidate per read. For each
read in `readIdsSortedByLength` order, we create a candidate representing its
full strand-0 journey interval `[0, journeySize)` and push it onto the heap.

The heap is ordered by:
1. **Base span** (primary, descending) — the genomic distance in bases between
   the first and last anchor in the interval, computed from marker positions.
2. **Read length** (secondary, descending) — raw sequence length of the read.
3. **OrientedReadId** (tertiary, ascending) — deterministic tie-breaking.

This ensures the largest genomic interval is always processed first.

Each read also has a **generation counter**, initialized to 0. Every candidate
records the generation at the time it was pushed. This is used to detect stale
candidates (see Step 3).

### Step 2: Anchor ownership tracking

A vector `anchorOwner[anchorId]` tracks which window owns each anchor. All
entries start as `anchorUnclaimed` (= `uint32_t::max`). Once an anchor is
claimed by a window, its entry is set to that window's ID and it is never
reassigned.

### Step 3: Main loop — pop, validate, create or re-push

Repeat until the heap is empty:

1. **Pop** the top candidate. It specifies a backbone oriented read and a
   journey interval `[begin, end)`.

2. **Generation check.** If the candidate's generation does not match the
   read's current generation counter, discard it — it is stale (the read was
   re-pushed with updated intervals after a partial claim). This avoids
   processing outdated intervals.

3. **Unclaimed check.** Scan every anchor in the candidate's interval. If all
   are still unclaimed, proceed to Step 4. If any anchor has been claimed since
   this candidate was pushed:
   - Increment the read's generation counter (invalidating all its older heap
     entries).
   - Scan the read's full journey for remaining contiguous unclaimed runs.
   - Push each run that has ≥ `minBackboneWindowAnchors` (= 2) anchors as a
     new candidate with the updated generation.
   - Go back to popping the next candidate.

4. **Create window.** See Step 4 below.

### Step 4: Window creation

Given a validated backbone interval `[seedBegin, seedEnd)` on oriented read
`backboneOid`:

#### 4a. Claim backbone anchors

For each journey position in `[seedBegin, seedEnd)`, look up the anchor ID
from the backbone's journey. If unclaimed, set `anchorOwner[anchorId] =
windowId`.

#### 4b. Find touching reads

Iterate over every anchor in the backbone interval. For each anchor, look up
all oriented reads that pass through it (via `Shasta2AnchorMarkerInfo`). For
each such read (excluding the backbone itself and entries with invalid journey
positions):

- Track the **minimum** and **maximum** journey positions where this read
  touches the backbone anchors.
- Track the **count** of backbone anchors this read shares.

This uses an epoch-based visited mechanism to avoid re-initializing per-read
arrays: a global `epoch` counter is incremented per window, and per-read
`touchedEpoch` values are compared against it.

#### 4c. Claim touching-read anchors

For each touching read, consider its journey range `[touchedMin, touchedMax+1)`.
Within that range, find contiguous runs of unclaimed anchors and claim them for
this window. Each such run becomes an `AnchorWindowReadInterval` added to the
window's `readIntervals`.

This means the window claims **all** unclaimed anchors in the touching read's
range, not just the ones shared with the backbone. If a touching read's journey
has anchors A-B-C-D-E and only B and D are shared with the backbone, anchors A,
C, and E are also claimed (if unclaimed) because they fall within the touched
range `[A, E]`. However, read intervals are only created for the **unclaimed
runs** — anchors already claimed by the backbone (B and D) are skipped. So the
touching read would get up to three read intervals: `[A, A]`, `[C, C]`, and
`[E, E]` (the gaps between the already-claimed backbone anchors).

#### 4d. Record the window

The completed `AnchorWindow` (backbone info, claimed anchor count, all read
intervals) is appended to the output vector.

### Step 5: Completion

After the heap is empty, every anchor is either owned by exactly one window or
unclaimed. The function logs summary statistics: total windows created, anchors
claimed, anchors unclaimed, total read intervals, and elapsed time.

## Properties

- **Disjoint anchor ownership.** Each anchor belongs to at most one window.
  Once claimed, an anchor is never reassigned.

- **Greedy longest-first.** The largest available base-span interval always
  becomes the next backbone. This tends to produce windows that cover large
  genomic regions with high-quality backbone reads.

- **Reads can span multiple windows.** A single read may contribute intervals
  to several windows if different parts of its journey are claimed by different
  windows.

- **Backbone reads can seed multiple windows.** If a backbone candidate is
  stale (partially claimed), its remaining unclaimed sub-intervals are
  re-pushed and may each become backbones for separate windows.

- **Minimum window size.** An interval must have ≥ 2 anchors
  (`minBackboneWindowAnchors`) to be pushed as a candidate. Single-anchor
  intervals are discarded.

## Example

Consider 3 reads with journeys (each letter is an anchor):

```
Read 1 (longest):  A B C D E F G H
Read 2:                C D E F
Read 3:                    E F G H I J
```

1. Read 1's full interval `[A, H]` has the largest base span. It becomes the
   backbone of Window 0. Anchors A–H are claimed.

2. Reads touching Window 0's anchors:
   - Read 2 touches C, D, E, F. Its touched range is `[C, F]`. All anchors in
     that range are already claimed by the backbone — no unclaimed runs exist.
     Read 2 gets **no** read interval and does not appear in this window.
   - Read 3 touches E, F, G, H. Its touched range is `[E, J]`. Within that
     range, E–H are already claimed (skipped), but I and J are unclaimed →
     claimed for Window 0. Read 3 gets a read interval `[I, J]` (only the
     unclaimed run).

3. All anchors are now claimed. The heap has no more valid candidates. Done.

Result: one window with backbone Read 1 `[A, H]` and Read 3 `[I, J]`. Read 2
does not appear because all its anchors were already claimed by the backbone.

**Note:** This means touching reads only contribute read intervals for anchors
that are **not** on the backbone. Reads fully contained within the backbone's
anchor range (like Read 2) will not appear in the window's `readIntervals` at
all — only the backbone covers those anchors.

## Pipeline Integration

In `main.cpp`, `computeAnchorWindows` is called after journey creation and
before anchor graph construction:

```cpp
vector<AnchorWindow> anchorWindows;
assembler.computeAnchorWindows(
    shasta2Anchors,
    shasta2Journeys,
    readIdsSortedByLength,
    anchorWindows,
    threadCount);
```

The resulting `anchorWindows` vector is available for downstream per-window
phasing.
