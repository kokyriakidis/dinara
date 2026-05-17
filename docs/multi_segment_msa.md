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

## Read Discovery

Overlapping reads are discovered directly from the anchor marker intervals
rather than from the window's `readIntervals`. For each backbone boundary
anchor, the function queries `anchorMarkerInfos[anchorId]` to get all
oriented reads containing that anchor, along with their marker ordinals.
This builds a map from each read to its sorted list of **boundary hits**
(backbone boundary index + marker ordinal) in a single pass over the
backbone anchors.

This approach avoids:
- Deduplicating reads from `readIntervals`
- Walking each read's full journey to find shared backbone anchors
- Calling `getOrdinal()` per shared anchor (ordinals come directly from
  the marker info struct)

Reads with fewer than 2 boundary hits are discarded — at least two shared
backbone anchors are needed to define a segment to align.

## Read Alignment

For each read with ≥2 shared backbone anchors, the function aligns one
segment per pair of consecutive boundary hits. Given consecutive hits at
backbone boundaries `i` and `j` with marker ordinals `o_i` and `o_j`:

1. Extracts the read's base sequence between ordinals `o_i` and `o_j`.
2. Calls `align_from(seq, nodeIds[i])` to align the read starting at
   boundary `i`'s segment node.

When `j - i = 1` the read segment corresponds to exactly one backbone
segment. When `j - i > 1` the read skips intermediate backbone anchors
and the alignment traverses multiple backbone segments in the graph.

Each `align_from` call updates the POA graph: matching bases increment node
weights, while divergent bases create new branches. After all reads are
aligned, the graph encodes all variation (SNPs, indels) across reads within
the window.

## Partial Coverage

Reads typically don't span the entire window. A read sharing anchors A5
and A8 only contributes segments between those boundaries. The `align_from`
call uses `is_ends_free=true`, so the aligner doesn't penalize the read
for not reaching the sink.

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
