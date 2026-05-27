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
6. **Output** — `Shasta2AnchorGraph.gfa`, `Shasta2AnchorGraph.csv` (Bandage colors), binary save
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

### 7. Per-window `outEdges` / `inEdges` Population
Populates from `createdEdges`, only for forward windows (< windowCount).

### 8. Per-window Connectivity Diagnostic
Prints per-window incoming/outgoing counts and transition flows.

### 9. Rule 1: Trim Backbone Outside Bounding Inter-window Span
For each forward window with inter-window edges:
- **Start bound:** min journey position of all incoming connection anchors (`anchorIdB`). If no incoming edges, defaults to backbone start (no head trimming).
- **End bound:** max journey position of all outgoing connection anchors (`anchorIdA`). If no outgoing edges, defaults to backbone end (no tail trimming).
- Trims backbone vertices outside `[start, end]` by removing only **intra-window** edges (edges where both endpoints belong to the same window). Inter-window edges are preserved.
- RC mirror vertices are also trimmed.

### 10. Edge Count and Validation
Recounts edge types (intra/inter/alt-path) after trimming. Warns about edges with no shared reads.

## Output Files

- `Shasta2AnchorGraph.gfa` — GFA with non-isolated vertices only, all edges `+/+` orientation, `RC:i:` coverage tag
- `Shasta2AnchorGraph.csv` — Bandage color CSV, HSL color per window, RC mirrors get same color as forward window
- `Shasta2AnchorGraph` (binary) — serialized anchor graph for shasta2

## Key Files

| File | Purpose |
|---|---|
| `src/AnchorWindows.hpp` | `AnchorWindow`, `AnchorWindowReadInterval`, `InterWindowEdge` structs |
| `src/Shasta2AnchorGraph.hpp` | Graph class declaration, `anchorToWindow`, `windowCount` |
| `src/Shasta2AnchorGraph.cpp` | Anchor-window constructor: edge creation, transitions, Rule 1 |
| `src/Shasta2AnchorGraphGfa.cpp` | `writeGfa()` and `writeCsv()` implementations |
| `src/AssemblerAnchorWindowsClean.cpp` | Window creation, backbone DP filtering |
| `src/AssemblerTestAnchorWindowsClean.cpp` | `writeAnchorWindowsCleanGfa()` (still exists but no longer called) |
| `srcMain/main.cpp` | Pipeline orchestration, early return after graph save |
| `src/AssemblerOptions.hpp/cpp` | CLI options including `minCommonForBackbone`, `maxSkipForBackbone` |
| `src/Assembler.hpp` | `computeAnchorWindowsClean` signature |
| `srcMain/AssemblerShasta2Anchors.cpp` | Shasta2 `Anchors` constructor call |

## CLI Options (Assembly.mode3)

- `--Assembly.mode3.minCommonForBackbone` (default: 2) — min common reads for consecutive backbone pairs in DP filtering
- `--Assembly.mode3.maxSkipForBackbone` (default: 10) — max forward look in backbone DP
- `--Assembly.mode3.minInterWindowCoverage` — min shared reads for inter-window edges

## Build

```bash
cd /workspaces/dinara/build && make -j$(nproc)
```

Binary: `build/Executable/dinara`

## Reverted/Abandoned Approaches

- **Contradictory pair removal** (bidirectional window connections W1→W2 and W2→W1) — caused full graph fragmentation
- **Backbone-only inter-window edge discovery** — caused breaks in connectivity
- **`clear_vertex` in Rule 1** — removed ALL edges including inter-window, fragmenting the graph. Fixed to only remove intra-window edges.
- **`transitionReads`-based 1-to-1 check** — stray reads added extra neighbors, causing zero 1-to-1 matches. Replaced with direct `outEdges`/`inEdges` counting, then generalized to bounding span approach.

## Current State

- Rule 1 (bounding span trimming) is active and working
- Post-graph steps (transitive reduction, assembly graph, etc.) are disabled via `return;` in main.cpp
- `writeAnchorWindowsCleanGfa` still exists in code but is no longer called — replaced by `Shasta2AnchorGraph.gfa` + `.csv`
- The graph is compatible with shasta2's `readExternalAnchors` via paired anchor IDs
