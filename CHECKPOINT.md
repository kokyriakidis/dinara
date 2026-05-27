# Dinara Anchor Window Pipeline — Checkpoint

## Overview

Dinara is a genome assembler. The current work focuses on building an **anchor graph from anchor windows** that is compatible with shasta2. Each window represents a contiguous region on a backbone read's journey, and the anchor graph connects these windows via inter-window edges discovered from read journeys.

## Architecture

### Pipeline Flow (in `srcMain/main.cpp`)

1. **Create anchors and journeys** — `Shasta2Anchors`, `Shasta2Journeys`
2. **Flag contained reads** — `flagContainedReads(1000, 0.8, 0, threadCount)` before window creation
3. **Compute anchor windows** — `computeAnchorWindowsClean()` with backbone filtering parameters `minCommonForBackbone` (default 2) and `maxSkipForBackbone` (default 10)
4. **Phasing / clustering / alternate paths** — per-window het SNP detection, read clustering, alternate path creation
5. **Build anchor graph** — `Shasta2AnchorGraph` constructor from anchor windows
6. **Export for shasta2** — external anchors (`Shasta2ExternalAnchors`) and anchor graph (`Shasta2ExternalAnchorGraph`) in shasta2-native binary format
7. **Early return** — post-graph steps (transitive reduction, assembly graph) are currently disabled via `return;`

### Key Data Structures

#### `AnchorWindow` (`src/AnchorWindows.hpp`)

Central data structure for a window. Key fields:

- `windowId`, `backboneOrientedReadId`, `backboneBegin`, `backboneEnd` — identity and span
- `filteredBackbonePositions` — DP-filtered backbone positions where consecutive pairs have ≥ `minCommonForBackbone` common reads. Intra-window edges use these positions.
- `readIntervals` — all reads overlapping this window, with `previousWindow`/`nextWindow` tracking
- `alternatePaths` — parallel chains between backbone anchors at het sites
- `transitionReads` — `map<pair<uint32_t, uint32_t>, vector<OrientedReadId>>` grouping reads by `(previousWindow, nextWindow)` transition. Used for detangling.
- `outEdges` / `inEdges` — `vector<InterWindowEdge>` tracking inter-window edges incident to this window. Each stores `otherWindow`, `anchorIdA`, `anchorIdB`, `readCount`.
- `readClusters` — per-window phasing clusters
- `computeClusterTransitions()` — aggregates transition counts per cluster

#### `Shasta2AnchorGraph` (`src/Shasta2AnchorGraph.hpp`)

Boost adjacency_list (bidirectional, vecS) with `Shasta2AnchorGraphEdge` properties. Key members:

- `anchorToWindow` — `vector<uint32_t>` mapping each anchor ID to its window ID (including RC mirrors at `windowId + windowCount`). Stored as member for use by `writeCsv`.
- `windowCount` — number of forward windows
- `noWindow` — sentinel value (`uint32_t max`)

### Anchor ID Scheme

Paired anchor IDs: `2*i` = canonical (strand 0), `2*i+1` = RC (strand 1). The RC of anchor `a` is `a ^ 1`.

### RC Mirror Windows

For each forward window `W` (ID < windowCount), a virtual RC mirror window exists at ID `W + windowCount`. Its backbone anchors are `anchorId ^ 1` of the forward window's anchors. This lets strand-1 reads discover inter-window transitions. RC mirror windows don't have `AnchorWindow` entries — they exist only in `anchorToWindow`.

## Anchor Graph Construction (`src/Shasta2AnchorGraph.cpp`)

The anchor-window constructor builds the graph in these steps:

### 1. Vertex Creation
One vertex per anchor (all anchors, including unmapped ones).

### 2. `anchorToWindow` Map
Maps each backbone anchor to its window ID, and each RC anchor to `windowId + windowCount`.

### 3. Intra-window Edges
For each window, creates edges between consecutive `filteredBackbonePositions` entries. Also creates RC mirror edges (reversed direction, flipped anchor IDs).

### 4. Alternate Path Edges
For het windows, creates edges along alternate path chains (forward + RC mirror).

### 5. Per-read Window Transition Tracking
Walks all read journeys, builds `windowSequence` (normalized window IDs per read), populates:
- `readIntervals[].previousWindow / nextWindow`
- `transitionReads` map

### 6. Inter-window Edge Discovery
Walks all read journeys (skipping contained reads), collects candidate anchor pairs per window pair, picks the candidate with the highest `Shasta2AnchorPair.size()`. Rejects pairs with zero shared reads or below `minInterWindowCoverage`.

### 7. `removeNegativeOffsets`
Called on every `Shasta2AnchorPair` at all 3 edge construction sites (journey-based, intra-window via `addEdgeIfValid`, inter-window). Removes oriented reads where the base position on anchor A exceeds the base position on anchor B, ensuring all edges have forward offsets.

### 8. Per-window `outEdges` / `inEdges` Population
Populates from `createdEdges`, only for forward windows (< windowCount).

### 9. Per-window Connectivity Diagnostic
Prints per-window incoming/outgoing counts and transition flows.

### 10. Rule 1: Trim Backbone Outside Bounding Inter-window Span
For each forward window with inter-window edges:
- **Start bound:** min journey position of all incoming connection anchors (`anchorIdB`). If no incoming edges, defaults to backbone start (no head trimming).
- **End bound:** max journey position of all outgoing connection anchors (`anchorIdA`). If no outgoing edges, defaults to backbone end (no tail trimming).
- Trims backbone vertices outside `[start, end]` by removing only **intra-window** edges (edges where both endpoints belong to the same window). Inter-window edges are preserved.
- RC mirror vertices are also trimmed.

### 11. Edge Verification
After construction, verifies all edges: for each oriented read on an edge, checks that its base position on anchor A < position on anchor B (forward offset). Reports backward, onlyA, onlyB, and neither counts.

### 12. Edge Count and Validation
Recounts edge types (intra/inter/alt-path) after trimming. Warns about edges with no shared reads.

## Shasta2 Export (`src/Shasta2AnchorGraphExport.cpp`)

### External Anchors
`Shasta2Anchors::writeExternalAnchors()` exports only canonical (even) anchors. Each anchor's markers are written with raw positions (midpoint - kHalf). Shasta2 loads these and creates paired anchors (canonical + RC) from each.

### Anchor Graph
`Shasta2AnchorGraph::saveForShasta2()` serializes the graph in shasta2's `MemoryMapped::Vector<char>` format with boost binary archive. The graph uses shasta2's `AnchorGraphBaseClass` (boost adjacency_list with vecS). Anchor IDs in the exported graph directly match shasta2's scheme: external anchor `n` (dinara anchor `2*n`) → shasta2 anchors `2*n` and `2*n+1`.

### Round-trip Verification
After serialization, the export deserializes the binary data back and compares every edge's `anchorIdA`, `anchorIdB`, and `orientedReadIds` against the original graph.

## Shasta2 Integration

### Position Convention
- Dinara stores anchor positions as midpoints: `rawMarkerPosition + k/2`
- External anchors store raw positions: `midpoint - k/2`
- Shasta2 loads raw positions and adds `k/2` to get midpoints
- RC anchor positions in shasta2: `readLength - canonicalMidpoint`

### Known Issue (Fixed)
When shasta2 loads external anchors, `assemblerInfo->k` must be set before calling `readExternalAnchors()`. Without this, `kHalf=0` and all positions are off by `k/2`, causing backward offsets on edges between `+` and `-` anchors. Fix: add `assemblerInfo->k = options.k;` in the external anchors branch of `Assembler::assemble()` (committed to `kokyriakidis/shasta2`).

### Running with Shasta2
```bash
# 1. Run dinara
dinara --input reads.fastq --threads N

# 2. Run shasta2 with dinara's output
shasta2 --input reads.fastq \
  --external-anchors-name /path/to/DinaraRun/Shasta2ExternalAnchors \
  --external-anchor-graph-name /path/to/DinaraRun/Shasta2ExternalAnchorGraph \
  --min-read-length 1000 --k 50 --threads N
```

The `--k` value must match dinara's k-mer length (default 50). The `--min-read-length` must match dinara's read length cutoff so that read IDs are consistent.

## Output Files

- `Shasta2AnchorGraph.gfa` — GFA with non-isolated vertices only, all edges `+/+` orientation, `RC:i:` coverage tag
- `Shasta2AnchorGraph.csv` — Bandage color CSV, HSL color per window, RC mirrors get same color as forward window
- `Shasta2ExternalAnchors` — binary external anchors for shasta2
- `Shasta2ExternalAnchorGraph` — binary anchor graph for shasta2

## Key Files

| File | Purpose |
|---|---|
| `src/AnchorWindows.hpp` | `AnchorWindow`, `AnchorWindowReadInterval`, `InterWindowEdge` structs |
| `src/Shasta2AnchorGraph.hpp` | Graph class declaration, `anchorToWindow`, `windowCount` |
| `src/Shasta2AnchorGraph.cpp` | Anchor-window constructor: edge creation, transitions, Rule 1 |
| `src/Shasta2AnchorGraphExport.cpp` | `saveForShasta2()` — binary export for shasta2 |
| `src/Shasta2AnchorGraphGfa.cpp` | `writeGfa()` and `writeCsv()` implementations |
| `src/Shasta2AnchorPair.hpp/cpp` | `Shasta2AnchorPair` with `removeNegativeOffsets()` |
| `src/Shasta2Anchors.hpp/cpp` | Anchor construction, `writeExternalAnchors()` |
| `src/Shasta2Journeys.hpp/cpp` | Journey construction from anchors |
| `src/AssemblerAnchorWindowsClean.cpp` | Window creation, backbone DP filtering |
| `srcMain/main.cpp` | Pipeline orchestration, early return after graph save |
| `src/AssemblerOptions.hpp/cpp` | CLI options including `minCommonForBackbone`, `maxSkipForBackbone` |
| `src/Assembler.hpp` | `computeAnchorWindowsClean` signature |

## CLI Options (Assembly.mode3)

- `--Assembly.mode3.minCommonForBackbone` (default: 2) — min common reads for consecutive backbone pairs in DP filtering
- `--Assembly.mode3.maxSkipForBackbone` (default: 10) — max forward look in backbone DP
- `--Assembly.mode3.minInterWindowCoverage` — min shared reads for inter-window edges

## Build

```bash
cd /workspaces/dinara/build && make -j$(nproc)
```

Binary: `build/Executable/dinara`

Test data: `tests/GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq`

## Shasta2 Fork

Repository: `kokyriakidis/shasta2` (fork of `paoloshasta/shasta2`)

Build:
```bash
cd /tmp/shasta2_work/build && cmake .. -DBUILD_ID="dinara-debug" && make -j$(nproc)
```

Key modifications in the fork:
- `src/Assembler.cpp` — Set `assemblerInfo->k = options.k` before `readExternalAnchors()` (the kHalf fix)
- `src/AnchorGraph.hpp/cpp` — Added `verifyEdgeOffsets()` for edge offset validation
- `src/AssemblyGraph.cpp` — Added pre-assembly diagnostic for failing edges

## Reverted/Abandoned Approaches

- **Contradictory pair removal** (bidirectional window connections W1→W2 and W2→W1) — caused full graph fragmentation
- **Backbone-only inter-window edge discovery** — caused breaks in connectivity
- **`clear_vertex` in Rule 1** — removed ALL edges including inter-window, fragmenting the graph. Fixed to only remove intra-window edges.
- **`transitionReads`-based 1-to-1 check** — stray reads added extra neighbors, causing zero 1-to-1 matches. Replaced with direct `outEdges`/`inEdges` counting, then generalized to bounding span approach.
- **`filterByShasta2HashedKmerChecker`** — removed entirely; was filtering anchors by shasta2's hashed k-mer checker but caused issues and was unnecessary.

## Current State

- Rule 1 (bounding span trimming) is active and working
- `removeNegativeOffsets` is called at all 3 edge construction sites
- Edge verification runs after graph construction (0 backward edges expected)
- Export to shasta2 format is working with round-trip verification
- Shasta2 assembly completes successfully with dinara's exported anchors and graph
- Post-graph steps (transitive reduction, assembly graph, etc.) are disabled via `return;` in main.cpp
