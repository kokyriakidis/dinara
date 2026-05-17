# Multi-Segment MSA

## Overview

`AssemblerMultiSegmentMSA.cpp` implements per-anchor-window multiple sequence
alignment using the multi-segment TheseusMSA API. Given an anchor window
(a contiguous interval on a backbone read's journey plus overlapping reads),
it builds a POA graph segmented at anchor boundaries and progressively aligns
each overlapping read into the graph.

## Graph Construction

An anchor window has a backbone read whose journey defines a sequence of
anchors A0, A1, ..., An. Between each pair of consecutive anchors there is a
**segment** — the substring of the backbone read from the midpoint of one
anchor's k-mer to the midpoint of the next:

```
Backbone read:  ──────[A0]──────[A1]──────[A2]──────[A3]──────
                       │  seg 0  │  seg 1  │  seg 2  │
```

The midpoint (`marker.position + k/2`) is used as the cut point so that
adjacent segments don't overlap — each anchor k-mer is split in half between
the preceding and following segments.

These segments are passed to the `TheseusMSA` multi-segment constructor, which
builds a compact POA graph:

```
source → [seg 0 bases] → [seg 1 bases] → [seg 2 bases] → sink
          node 1           node 2           node 3
```

Each node contains the full base sequence of one inter-anchor segment. The
constructor returns the `NodeId` for each segment node, which is used to
target `align_from` calls.

## Read Alignment

Each non-backbone read in the window shares some subset of the backbone's
anchors. The function identifies the read's **entry boundary** (first shared
anchor) and **exit boundary** (last shared anchor), then:

1. Extracts the read's base sequence between its own marker ordinals for the
   entry and exit anchors.
2. Calls `align_from(seq, nodeIds[entryBoundary])` to align the read starting
   at the correct segment node.

The POA aligner handles the fact that the read's segment lengths may differ
from the backbone's due to insertions and deletions — the alignment naturally
crosses node boundaries at whatever position produces the best score.

Each `align_from` call updates the POA graph: matching bases increment node
weights, while divergent bases create new branches. After all reads are
aligned, the graph encodes all variation (SNPs, indels) across reads within
the window.

## Partial Coverage

Reads typically don't span the entire window. A read entering at anchor A1
and exiting at A2 only covers segments 1. The `align_from` call starts at
segment 1's node with `is_ends_free=true`, so the aligner doesn't penalize
the read for not reaching the sink. The partial backtrace mechanism handles
this case.

## Output

The function produces two output files per window:

- `testMultiSegmentMSA_window{id}.fasta` — Column-aligned MSA of all
  sequences (backbone + reads) in FASTA format.
- `testMultiSegmentMSA_window{id}.gfa` — The final POA graph in GFA format,
  showing all branches created by read variation.

## Source Files

| File | Role |
|------|------|
| `src/AssemblerMultiSegmentMSA.cpp` | Implementation of `testMultiSegmentMSA` |
| `src/AnchorWindows.hpp` | `AnchorWindow` and `AnchorWindowReadInterval` structs |
| `src/AssemblerAnchorWindows.cpp` | `computeAnchorWindows` — window partitioning |

## Dependencies

Requires the `theseus-lib` pericles branch with multi-segment extensions
([kokyriakidis/theseus-lib-multi-segment](https://github.com/kokyriakidis/theseus-lib-multi-segment),
branch `pericles`), which provides:

- `TheseusMSA(penalties, heuristics, segments, node_ids, weight)` —
  multi-segment constructor
- `align_from(seq, start_node, weight, is_ends_free, start_offset)` —
  alignment starting at an interior graph node
