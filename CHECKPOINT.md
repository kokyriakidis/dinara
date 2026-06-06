# Dinara — Checkpoint

## Overview

Dinara is a genome assembler with two active workstreams:

1. **Anchor window pipeline** — builds an anchor graph compatible with shasta2 from long reads
2. **SV detection from short reads** — the `--command svanchors` subcommand maps short reads to a reference and detects structural variants

## Architecture

### Pipeline Flow (in `srcMain/main.cpp`)

1. **Create anchors and journeys** — `Shasta2Anchors`, `Shasta2Journeys`
2. **Flag contained reads** — `flagContainedReads(1000, 0.8, 0, threadCount)` before window creation
3. **Compute anchor windows** — `computeAnchorWindowsClean()` with backbone filtering parameters `minCommonForBackbone` (default 2) and `maxSkipForBackbone` (default 10)
4. **Phasing / clustering / alternate paths** — per-window het SNP detection, read clustering, alternate path creation
5. **Compute window transitions** — `computeWindowTransitions()` walks all journeys and populates `transitionReads`, per-read `previousWindow`/`nextWindow`, and `backbonePreviousWindow`/`backboneNextWindow` on each `AnchorWindow`. Runs before graph construction so detangling has transition data available.
6. **Build anchor graph** — `Shasta2AnchorGraph` constructor from anchor windows. Recomputes transitions from its own `anchorToWindow` (clears and repopulates).
7. **Detangle** — `detangleWindows()` splits backbone anchors of tangled windows, then rebuilds the graph with parallel chains per path
8. **Export for shasta2** — external anchors (`Shasta2ExternalAnchors`) and anchor graph (`Shasta2ExternalAnchorGraph`) in shasta2-native binary format
9. **Early return** — post-graph steps (transitive reduction, assembly graph) are currently disabled via `return;`

### Key Data Structures

#### `AnchorWindow` (`src/AnchorWindows.hpp`)

Central data structure for a window. Key fields:

- `windowId`, `backboneOrientedReadId`, `backboneBegin`, `backboneEnd` — identity and span
- `baseSpan` — base distance from first to last backbone anchor on the backbone read
- `backbonePreviousWindow` / `backboneNextWindow` — normalized window IDs the backbone read transitions to/from. Recomputed after each filter step via `recomputeBackboneEndpoints()`. `noWindow` if the backbone starts/ends here.
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

### 10. Inter-window Edge Filter Pipeline

All edge disabling goes through `disableEdge()`, which sets `useForAssembly=false` on the edge and its RC mirror (`dst^1 → src^1`). Edges remain in the graph but are excluded from assembly and GFA/CSV output.

The pipeline runs `trimBackbones()` between each filter and `recomputeBackboneEndpoints()` before filters that use backbone endpoint fields.

```
trimBackbones()
recomputeBackboneEndpoints()
removeDeadEndSpurs()            // Remove single-neighbor internal spurs
trimBackbones()
recomputeBackboneEndpoints()
runSingleEdgeFilter()           // Case 2: remove single-point connections
trimBackbones()
recomputeBackboneEndpoints()
runBypassDetourFilter()         // Case 1: bypass detours through other windows
trimBackbones()
recomputeBackboneEndpoints()
runBubblePopFilter()            // Case 3: pop bubbles returning to same window
trimBackbones()
recomputeBackboneEndpoints()
runShortcutFilter()
trimBackbones()
recomputeBackboneEndpoints()
runCrossWindowFilter()
trimBackbones()
removeIsolatedWindows()
trimBackbones()
removeSmallWindows(2)
trimBackbones()
removeDanglingWindowsIterative("post-filter")
```

#### Graph Surgery Filters (Cases 1–3)

Three filters that restructure inter-window edges to linearize the graph:

#### `runSingleEdgeFilter()` (Case 2)
If two windows are connected by exactly one inter-window edge, it's a spurious single-point connection. Delete it and its RC mirror.

#### `runBypassDetourFilter()` (Case 1)
Walking window `w`'s backbone, if window X has an incoming edge at backbone anchor `a_i` and an outgoing edge at a later anchor `a_j`, then X's path detours through `w`. Create a bypass edge in X (connecting X's anchors on either side of the detour), then remove the inter-window edges. Window `w`'s backbone stays intact.

```
Before: Window X: ... → x1          x2 → ...
                         |           ↑
        Window w: a0 → a_i → ... → a_j → a4

After:  Window X: ... → x1 ------→ x2 → ...  (bypass edge created)
        Window w: a0 → a_i → ... → a_j → a4  (backbone intact)
```

#### `runBubblePopFilter()` (Case 3)
For each backbone anchor with inter-window edges, BFS through the anchor graph (limited by windows traversed, max 3). If the BFS reaches a later backbone anchor of the same window, a bubble exists. Keep the inter-window path, disable intra-window edges of intermediate backbone anchors.

```
Before: a0 → a1 → a2 → a3 → a4  (backbone of w)
              |              ↑
              → [X] ------→ |

After:  a0 → a1             a3 → a4  (linearized)
              |              ↑
              → [X] ------→ |
        a2 disconnected from backbone chain
```

#### `trimBackbones()`
Directional head/tail trim respecting backbone orientation. The backbone chain always goes forward (position 0 → position N-1), so incoming inter-window edges arrive at the head and outgoing inter-window edges depart from the tail.

- **Head trim**: walks from position 0 inward, trims anchors before the first one with an **incoming** inter-window edge. If no anchor has an incoming inter-window edge, head trim is skipped (resets to 0).
- **Tail trim**: walks from the last position inward, trims anchors after the last one with an **outgoing** inter-window edge. If no anchor has an outgoing inter-window edge, tail trim is skipped.
- **RC direction flip**: when checking the RC mirror anchor `aid^1`, edge direction is inverted — an out-edge on `aid^1` counts as an in-edge on `aid`, and vice versa, because the RC of edge `A → B` is `B^1 → A^1`.
- **Deletion**: disables all edges of trimmed anchors and their RC mirrors via `disableAllEdges(aid)` + `disableAllEdges(aid^1)`. Trimmed anchors are guaranteed to have no inter-window edges (otherwise the trim would have stopped), so only intra-window edges are destroyed. `disableEdge` ensures each disabled edge also disables its RC mirror.

#### `recomputeBackboneEndpoints()`
For each window, walks the backbone read's full journey, builds a window sequence from anchors that still have active edges, and updates `backbonePreviousWindow`/`backboneNextWindow`. Must run before any filter that uses these fields.

#### `runShortcutFilter()`
A window is a **shortcut** if `prevW ≠ nextW`, both are not `noWindow`, and `prevW` and `nextW` are directly connected (have an anchor-level edge between them). The window is a redundant bypass. Removes only edges between `w↔prevW` and `w↔nextW`; edges to other windows are preserved.

#### `runCrossWindowFilter()`
A window is a **cross-window** if `prevW ≠ nextW`, both are not `noWindow`, and `prevW` and `nextW` are NOT connected. The backbone bridges unrelated regions. Removes all inter-window edges of `w`, making it isolated for cleanup by `removeIsolatedWindows()`.

#### `removeIsolatedWindows()`
Finds windows with no remaining active inter-window edges and disables all their backbone edges (forward + RC).

#### `removeSmallWindows(maxSize)`
Removes windows with ≤ `maxSize` backbone anchors. Safety check: neighbor must have other connections so removal doesn't disconnect the graph.

#### `removeDanglingWindowsIterative(label)`
Iteratively removes windows with edges on only one side. Safety check: neighbor must have other connections on the same side from a non-dangling window.

### 11. Edge Verification
After construction, verifies all edges: for each oriented read on an edge, checks that its base position on anchor A < position on anchor B (forward offset). Reports backward, onlyA, onlyB, and neither counts.

### 12. Edge Count and Validation
Recounts edge types (intra/inter/alt-path) after trimming. Warns about edges with no shared reads.

## Detangling (`src/DinaraDetangle.hpp/cpp`)

After the initial anchor graph is built (including `transitionReads` population), detangling creates bypass edges around tangled windows using Verkko-style triplet resolution.

### When a Window Is Tangled

A window is a detangling candidate if it has ≥ 2 distinct predecessors OR ≥ 2 distinct successors (from `transitionReads` entries where `key != noWindow`). Linear windows (≤1 pred AND ≤1 succ) are skipped. Hairpin windows (a neighbor normalizes to the window itself, i.e., transitions to/from its own RC mirror) are also skipped.

### Verkko-Style Triplet Resolution

The algorithm uses iterative resolve steps with decreasing `minEdgeSupport` thresholds (default: {20, 10, 5}), resolving easy tangles first.

For each non-linear, non-hairpin window X:

1. **Count edge coverage**: Only full triplets `(pred, succ)` where both sides are known contribute. Entries with `noW` on either side (reads starting/ending inside the window) are excluded — they cannot form triplets and would inflate edge coverage, blocking resolution.

2. **Identify solid triplets**: A triplet `(pred, X, succ)` is "solid" if its read count ≥ `minEdgeSupport`.

3. **Check resolvability**: Every significant edge (coverage ≥ `minEdgeCoverage`, default 5) must be covered by at least one solid triplet. An edge is skippable if its coverage < `minEdgeCoverage` AND the neighbor window is "removable" (average anchor coverage < `minEdgeCoverage`).

4. **Find anchor pairs**: For each solid triplet, walk flow reads' journeys to find the best `(lastAnchorInPred, firstAnchorInSucc)` pair. The most popular pair across reads is selected.

5. **Create bypass edges**: Forward bypass `A→B` plus RC mirror `rc(B)→rc(A)`.

6. **Remove flow reads**: Flow reads are removed from the bypassed window's backbone anchors (both canonical and RC) in a mutable copy of `anchorMarkerInfos`. A `commonReads` recomputation ensures only reads still present on both bypass anchors are removed.

7. **Rebuild**: After all resolve steps, `anchorMarkerInfos` is rebuilt from the modified copy.

### Key Design Decisions

- **Edge coverage excludes `noW` entries**: Reads that enter a window but don't exit (or vice versa) cannot form triplets. Including them inflated edge coverage, making edges look significant but uncoverable, which blocked nearly all resolutions.
- **Significance check iterates `coveredInNeighbors`/`coveredOutNeighbors`** (maps built from full triplets only), not the broader `predecessors`/`successors` sets. Predecessors that only appear in `(pred, noW)` entries are invisible to the check.
- **No "borders unresolvable" check**: Verkko prevents resolving adjacent tangles simultaneously. We don't implement this — the `commonReads` recomputation provides some protection, and our windows are larger than Verkko's k-mer nodes.
- **No hairpin resolution**: Verkko resolves some hairpins by unfolding them into paired fw/bw nodes. We skip hairpins entirely (conservative).
- **Iterative steps allow progressive resolution**: A window partially resolved at `minEdgeSupport=20` can have additional triplets resolved at lower thresholds.

### 0-Full-Triplet Windows

Many windows show "2 preds, 2 succs" but have 0 full-triplet predecessors/successors. This means no reads span the entire window from predecessor to successor — the "2 preds" come from `(pred, noW)` entries (reads entering but ending inside) and the "2 succs" from `(noW, succ)` entries. These aren't real tangles — the windows are too long for reads to span. They pass the significance check harmlessly (no full-triplet edges to block on) and produce no candidates.

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

- `Shasta2AnchorGraph.gfa` — GFA with non-isolated vertices only, all edges `+/+` orientation, tags: `RC:i:` (coverage), `wn:i:` (window ID, `noWindow` if unassigned), `ws:Z:` (strand: `fw`/`rc`/`none`)
- `Shasta2AnchorGraph.csv` — Bandage color CSV, HSL color per window, RC mirrors get same color as forward window
- `Shasta2AnchorGraph-bidirected-pre-bypass.gfa` — bidirected graph before bypass detour filter
- `Shasta2AnchorGraph-bidirected.gfa` — bidirected graph after bypass detour filter
- `Shasta2AnchorGraph-bidirected.csv` — Bandage color CSV for bidirected graph (by window)
- `Shasta2AnchorGraph-bidirected-byread.csv` — Bandage color CSV for bidirected graph (by backbone read)
- `Shasta2AnchorGraph-unitigs.gfa` — unitigified bidirected graph with `LN:i:`, `RC:i:`, `wn:Z:` tags
- `Shasta2AnchorGraph-unitigs.csv` — Bandage color CSV for unitig graph (by window)
- `Shasta2AnchorGraph-unitigs-byread.csv` — Bandage color CSV for unitig graph (by backbone read)
- `Shasta2ExternalAnchors` — binary external anchors for shasta2
- `Shasta2ExternalAnchorGraph` — binary anchor graph for shasta2

## Key Files

| File | Purpose |
|---|---|
| `src/AnchorWindows.hpp` | `AnchorWindow`, `AnchorWindowReadInterval`, `InterWindowEdge` structs |
| `src/Shasta2AnchorGraph.hpp` | Graph class declaration, `anchorToWindow`, `windowCount` |
| `src/Shasta2AnchorGraph.cpp` | Anchor-window constructor: edge creation, transitions, Rule 1 |
| `src/Shasta2AnchorGraphExport.cpp` | `saveForShasta2()` — binary export for shasta2 |
| `src/Shasta2AnchorGraphGfa.cpp` | `writeGfa()`, `writeCsv()`, and `toBidirected()` implementations |
| `src/BidirectedAnchorGraph.hpp` | Bidirected graph class, `Unitig` struct, `OrientedAnchor`, `BidirectedAnchorId` |
| `src/BidirectedAnchorGraph.cpp` | `unitigify()`, `writeUnitigGfa()`, `bypassDetourFilter()`, GFA/CSV writers |
| `src/WindowTransitions.hpp/cpp` | `computeWindowTransitions()` — standalone transition counting from journeys, decoupled from graph |
| `src/DinaraDetangle.hpp/cpp` | `detangleWindows()` — split backbone anchors of tangled windows |
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

## Bidirected Anchor Graph

After the directed `Shasta2AnchorGraph` is built and filtered, it is converted to a bidirected graph (`BidirectedAnchorGraph`) matching MBG/Verkko's design. This collapses each anchor and its RC mirror into a single bidirected node.

### Core Types

- **`BidirectedAnchorId`**: `directedId / 2`. Conversions: `fromDirected(d) = d/2`, `toDirected(b) = 2*b`, `forwardDirected(b) = 2*b`, `rcDirected(b) = 2*b+1`.
- **`OrientedAnchor`**: `pair<BidirectedAnchorId, bool>` — node + traversal direction (true = forward, false = reverse).
- **`reverseAnchor({v, s})`**: returns `{v, !s}`.

### Edge Semantics

An edge `(u, s) → (v, t)` means: leave node `u` on side `s`, arrive at node `v` on side `!t`, traverse `v` in direction `t`. The RC mirror of this edge is `(v, !t) → (u, !s)`.

- **`getNeighbors({v, s})`**: returns all `OrientedAnchor` targets reachable by edges leaving `v` on side `s`.
- **Entry-side neighbors**: `getNeighbors({v, !orient})` returns RC-mirror targets; the real source is `reverseAnchor(neighbor)`.

### Node and Edge Properties

- `NodeProperties`: `windowId`, `backboneReadId`.
- `EdgeProperties`: `useForAssembly`, `offset` (base-pair distance between anchors).

### Conversion from Directed Graph

`Shasta2AnchorGraph::toBidirected()` in `src/Shasta2AnchorGraphGfa.cpp`:
1. Creates one bidirected node per anchor pair, transfers `windowId` and `backboneReadId`.
2. Converts each directed edge `A → B` to bidirected edge with proper orientation, transfers `offset`.
3. Only edges with `useForAssembly=true` are converted.

### Bypass Detour Filter (Bidirected)

`BidirectedAnchorGraph::bypassDetourFilter()` runs on the bidirected graph before unitigification. Same concept as the directed-graph `runBypassDetourFilter()` but operates on bidirected edges.

For each window `w`, walks the backbone anchors and collects inter-window edges grouped by neighbor window. For each neighbor window X that enters `w` at backbone position `i` and exits at position `j > i`:
1. Finds the closest exit after each entry (may be multiple exits at the same backbone index).
2. Creates a bypass edge `reverseAnchor(inE.xNeighbor) → outE.xNeighbor` in X.
3. Removes the entry and exit inter-window edges.

Output: `Shasta2AnchorGraph-bidirected-pre-bypass.gfa` (before) and `Shasta2AnchorGraph-bidirected.gfa` (after).

### Unitigification

`BidirectedAnchorGraph::unitigify()` collapses linear chains into unitig segments:

1. **Internal node detection**: A node is internal if it has degree 1 on both sides and its sole neighbor on each side also has degree 1 on the connecting side.
2. **Chain walking**: Starting from non-internal nodes, walks forward through internal nodes to build chains.
3. **RC-mirror deduplication**: Each chain and its RC mirror represent the same unitig. Chains are canonicalized by lexicographic comparison; duplicates are skipped.

Each `Unitig` stores:
- `chain` — vector of `OrientedAnchor` forming the unitig path.
- `windowSequence` — ordered list of distinct window IDs traversed.
- `averageCoverage` — mean anchor coverage across the chain.
- `totalOffset` — sum of edge offsets along the chain (base-pair length).

### Unitig GFA Output

`writeUnitigGfa()` writes segments and links:
- **Segments**: one S-line per unitig with `LN:i:` (total offset), `RC:i:` (average coverage), `wn:Z:` (comma-separated window sequence).
- **Links**: discovered by walking exit-side neighbors of each unitig's endpoints and looking up the neighbor in an entry map.

The entry map registers four entry points per unitig to account for RC-mirror deduplication:
```
front                → i+  (enter unitig forward from left)
reverseAnchor(back)  → i-  (enter unitig reversed from right)
reverseAnchor(front) → i-  (RC mirror's back entry)
back                 → i+  (RC mirror's front entry)
```

Links are deduplicated: each link and its RC mirror are the same, so only the canonical form is emitted.

### CSV Outputs

- `writeCsv()` / `writeCsvByRead()` — per-node Bandage color CSVs for the raw bidirected graph. `writeCsvByRead` colors nodes by backbone read using golden-angle hue spacing.
- `writeUnitigCsv()` / `writeUnitigCsvByRead()` — same for unitig segments.

### Pipeline Order (in `main.cpp`)

```cpp
auto bidirectedGraph = shasta2AnchorGraph->toBidirected(anchorWindows, *shasta2Journeys);
bidirectedGraph.writeGfa("Shasta2AnchorGraph-bidirected-pre-bypass.gfa");
bidirectedGraph.bypassDetourFilter(anchorWindows, *shasta2Journeys);
bidirectedGraph.writeGfa("Shasta2AnchorGraph-bidirected.gfa");
bidirectedGraph.writeCsv(...);
bidirectedGraph.writeCsvByRead(...);
const auto unitigs = bidirectedGraph.unitigify();
bidirectedGraph.writeUnitigGfa("Shasta2AnchorGraph-unitigs.gfa", unitigs, ...);
bidirectedGraph.writeUnitigCsv(...);
bidirectedGraph.writeUnitigCsvByRead(...);
```

### Key Files

| File | Purpose |
|---|---|
| `src/BidirectedAnchorGraph.hpp` | Class declaration, `Unitig` struct, `OrientedAnchor`, `BidirectedAnchorId` |
| `src/BidirectedAnchorGraph.cpp` | `unitigify()`, `writeUnitigGfa()`, `bypassDetourFilter()`, GFA/CSV writers |
| `src/Shasta2AnchorGraphGfa.cpp` | `toBidirected()` conversion, directed-graph `writeGfa()`/`writeCsv()` |

## Current State (Anchor Windows)

- Inter-window edge filter pipeline is active: shortcut, parallel, cross-window filters using backbone endpoints, with trimBackbones and recomputeBackboneEndpoints between each step
- All edge disabling uses `disableEdge()` member function for RC mirror symmetry
- Bypass edges and Case 2 detangle are disabled (`#if 0`)
- Detangling splits tangled windows by through-flow paths (through-flows only; start/end reads not yet assigned)
- `removeNegativeOffsets` is called at all 3 edge construction sites
- Edge verification runs after graph construction (0 backward edges expected)
- Export to shasta2 format is working with round-trip verification
- Shasta2 assembly completes successfully with dinara's exported anchors and graph
- Bidirected graph conversion, bypass detour filter, unitigification, and RC-mirror deduplication are active
- Post-graph steps (transitive reduction, assembly graph, etc.) are disabled via `return;` in main.cpp

---

## SV Detection Pipeline (`--command svanchors`)

### Overview

Detects structural variants (insertions and deletions) from short reads aligned to a local reference. The pipeline uses k-mer-based chaining (hifiasm-derived DP) with minimap2-sr scoring to align reads, then applies multiple detection layers.

### Running

```bash
dinara --command svanchors \
  --reference reference.fa \
  --input reads.fa \
  --bam region.bam \
  --assemblyDirectory outdir \
  --Kmers.k 10 --Kmers.minimizerW 6
```

| Flag | Required | Description |
|------|----------|-------------|
| `--command svanchors` | yes | SV calling mode (not the default assembly mode) |
| `--reference` | yes | Local reference FASTA for the region |
| `--input` | yes | Reads FASTA extracted for the region |
| `--bam` | yes | Region BAM — provides CIGAR indels, depth, and SA-tag evidence. Without this, most detection sources (early-CIGAR, covdrop, depth-scan) are missing. |
| `--assemblyDirectory` | yes | Output directory. Must be unique per concurrent run to avoid `DinaraRun` conflicts. |
| `--Kmers.k 10` | yes | K-mer size for minimizer markers |
| `--Kmers.minimizerW 6` | yes | Minimizer window size |

Output is written to stdout/stderr. SV calls appear as lines matching `>>> {INS,DEL} CALL`. Parse with:
```bash
grep -E ">>> .* CALL" output.log
```

Test cases: `test_cases/` in the repo root, and on the cluster under `.../full_ins_eval/cases/` and `.../full_del_eval/cases/`.

### Pipeline Steps (in `srcMain/main.cpp`, svanchors section)

1. **Load reference** (read 0) and short reads
2. **Find minimizer markers** (k=10, w=6)
3. **Remove non-unique reference k-mers** — blacklist k-mers appearing >1 time in the reference, then rescue freq-2 k-mers in large VNTR-depleted regions (≥5 consecutive depleted 50bp windows)
4. **Build inverted index** — hash table mapping canonical k-mers to (readId, position) occurrences
5. **K-mer hit depth profiling** — per-window hit depth along the reference for breakpoint detection
6. **DP chaining** — minimap2-sr scoring mode (chainingMode=1), bw=100, maxGap=5000, multi-chain extraction
7. **Replace extended coordinates** — use raw anchor positions instead of hifiasm-extended qs/qe/ts/te
8. **Classify chains** — primary/supplementary/secondary based on query overlap (minimap2 mm_set_parent)
9. **Build read graph** — connect reads sharing chains for path-based insertion sizing
10. **Build SV MSA** (`buildSvMSA` in `src/AssemblerSvAnchors.cpp`) — the main detection engine

### Chaining Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| chainingMode | 1 | minimap2-sr scoring |
| minimap2Bw | 100 | Fixed bandwidth; limits single-chain deletion detection to ~100bp |
| minimap2MaxGap | 5000 | Large gaps for SV detection |
| minimap2MinChainScore | 25 | minimap2 -m25 |
| mcopyNum | 20 | minimap2 -N20 (max secondary chains) |
| mcopyRate | 0.50 | minimap2 -p.5 (secondary ratio) |
| downsampleHighFrequencyMarkers | false | Disabled; coverageHet not computed in svanchors |
| referenceReadCount | 0 | Chain all-vs-all (needed for read graph) |

### Detection Layers (in `buildSvMSA`)

#### Per-read diagonal analysis
For each read's chain anchors, compute diagonal (refPos - readPos). A large rise/drop in diagonal indicates a deletion/insertion.

#### Coverage-drop breakpoint detection
Per-50bp-window analysis of chain endpoints and spanning counts. Left breakpoints (chain ends spike) and right breakpoints (chain starts spike) are paired.

#### K-mer hit depth profiling
Windows where hit depth drops below 50% of median indicate marker-depleted regions (VNTRs, repeats). Used for VNTR gap detection and marker-density gating.

#### VNTR gap detection
When left/right breakpoints are >500bp apart and >50% of intervening windows have low hit depth, the gap is classified as a VNTR. Path-based sizing is skipped for VNTRs.

#### Marker-density gating (commit `b2677da`)
Even when breakpoints are close (<500bp), if >60% of windows in a ±300bp radius have zero reference markers, the region is treated as VNTR-like and path-based sizing is skipped. Prevents false insertion calls from read-to-read chain shortcuts in repetitive regions.

#### Path-based insertion sizing
For non-VNTR breakpoint pairs, search the read graph for paths from left-flank reads to right-flank reads (1-3 hops). Insertion size = left overhang + intermediate extensions + right overhang - overlap spans. Gap-capped overlap prevents breakpoint-spanning gaps from inflating the overlap estimate.

#### Adaptive multi-k gap filling
For deletion sizing, progressively try k=60 down to k=10 to find skeleton anchors in marker-depleted regions. Larger k provides more specific anchors; smaller k fills gaps.

#### Split-read deletion detection (commit `e1284fc`)
Groups chains by read ID. For reads with 2+ same-strand chains, compares median diagonals between consecutive non-overlapping chain pairs. Diagonal differences >50bp and <1000bp with consistent reference gap indicate deletions. Requires ≥2 supporting reads.

#### Breakpoint pairing with fallback candidates
Breakpoint pairs are formed nearest-first. When the nearest right BP is out of range, remaining candidates are tried in strength-order (fold enrichment × ovhReadCount²). This allows strong distant BPs to be reached when the nearest is invalid.

#### Deletion-like pair suppression of single-breakpoint insertion
When a strong left-right BP pair exists (both fold ≥ 1.5, both ovhReads ≥ 2, refGap ≥ 100bp) but no path is found, the pattern is more consistent with a deletion than an insertion. The single-breakpoint insertion detector is suppressed so the downstream coverage-drop deletion detector can fire instead. Fixes cases like DEL137 that were previously misclassified as insertions.

#### Hit-depth-only breakpoint generation
When chain-endpoint breakpoints are missing on one or both sides, breakpoints are generated from hit-depth cluster edges. The left edge of a contiguous hit-depth drop cluster (minRatio < 0.3) becomes a left BP and the right edge becomes a right BP. Overhang reads are collected from reads spanning across or starting just after the drop zone.

#### Alternative-pairing insertion refinement
After an initial insertion call from the nearest BP pairing, alternative right BPs further away are tried. If a further BP produces a larger valid path distance, the larger insertion call is emitted instead. This handles insertions larger than the read length where the nearest BP captures only a partial insertion.

#### Tandem repeat weighted median cluster selection (commit `e1284fc`)
In marker-depleted regions (markerDepleted flag), uses weighted median of qualifying deletion clusters (≥100bp) instead of highest-support cluster. Improves sizing accuracy for tandem repeat deletions.

#### Per-anchor pairwise diagonal difference analysis (commit `69e589e`)
For each consecutive anchor pair in a chain, computes the diagonal difference. Collects all pairwise differences and clusters them for deletion sizing.

#### Inversion filter (commit `69e589e`)
Requires ≥15 anchors and ≥80bp span for inversion calls. Prevents small spurious inversion clusters.

#### DUST-gated STR detection (commit `4f3b22f`)
Runs SDUST (symmetric DUST algorithm, `src/Sdust.hpp`) on the reference to find low-complexity intervals ≥200bp. For each, finds reads with chain anchors on both flanks and compares read-space gap to reference-space gap. The difference gives insertion/deletion size.

#### VNTR depth estimation
For VNTR gaps where path-based sizing fails, estimates insertion size from total read bases in the VNTR vs expected bases at flanking coverage. Currently limited by BAM extraction depleting VNTR reads (short-read aligners assign low MAPQ in tandem repeats).

#### Breakpoint pairing with fallback candidates
Breakpoint pairs are formed nearest-first. When the nearest right BP is out of range, remaining candidates are tried in strength-order (fold enrichment × ovhReadCount²). This allows strong distant BPs to be reached when the nearest is invalid.

#### Deletion-like pair suppression of single-breakpoint insertion
When a strong left-right BP pair exists (both fold ≥ 1.5, both ovhReads ≥ 2, 100bp ≤ refGap ≤ 500bp) near the strongest BP, the single-breakpoint insertion detector is suppressed so the downstream coverage-drop deletion detector can fire instead. Fixes cases like DEL137 that were previously misclassified as insertions.

#### Hit-depth-only breakpoint generation
When chain-endpoint breakpoints are missing on one or both sides, breakpoints are generated from hit-depth cluster edges. The left edge of a contiguous hit-depth drop cluster (minRatio < 0.3) becomes a left BP and the right edge becomes a right BP. Overhang reads are collected from reads spanning across or starting just after the drop zone.

#### Alternative-pairing insertion refinement
After an initial insertion call from the nearest BP pairing, alternative right BPs further away are tried. If a further BP produces a larger valid path distance, the larger insertion call is emitted instead. Handles insertions larger than the read length where the nearest BP captures only a partial insertion.

#### Deletion fallback from per-read classifications
When a breakpoint pair has refGap ≥ 50bp, no path found, and no spanning chains for diagonal-shift detection (deletion larger than read length), falls back to per-read DEL classifications. Reads classified as DEL with breakpoints near the BP pair zone provide the deletion size estimate. Fixes DEL324 which had no primary detection.

#### Marker-depleted region detection
A coverage-drop region is marker-depleted if either all windows have zero reference markers or >50% of windows with markers have low hit depth. Previously, regions with zero reference markers were not flagged.

#### VNTR-like coverage-drop suppression
Large coverage-drop regions (>500bp) with minRatio ≈ 0 and strong edge BPs with low spanning counts (<10) are suppressed as VNTR-like. Prevents false deletion calls from VNTR regions where chains don't span due to repeat structure rather than a real deletion.

#### SA-tag DEL suppression in VNTR regions
SA-tag DEL calls are suppressed when the region is marker-depleted (≥10 hit-depth BPs with ≥30 unanchored reads) or when a VNTR gap was detected during breakpoint pairing. In VNTRs, the aligner maps split reads to different repeat copies, producing false DEL calls. Fixes INS235 and INS65 which had false SA-tag DEL calls.

#### Hit-depth drop zone span for large insertions
When a path-based insertion call is small (<200bp) but a further right BP exists with a hit-depth drop zone between L and R, the insertion size is re-estimated from the drop zone span. The contiguous region of low hit-depth (below 50% of median) between the breakpoints approximates the insertion size. Tries all valid candidate right BPs and picks the one with the largest drop span. Fixes INS454 (156bp → 400bp) and INS277 (56bp → 250bp).

#### Adaptive marker rescue in VNTR-depleted regions
After blacklisting non-unique reference k-mers, some reference windows have zero remaining markers (common in VNTRs where all k-mers are repetitive). Phase 2b of `removeNonUniqueReferenceMarkers` identifies contiguous depleted regions (≥5 consecutive 50bp windows with zero markers, i.e. ≥250bp) and rescues k-mers with reference frequency exactly 2 by removing them from the blacklist. Only k-mers that have at least one occurrence in a VNTR-depleted window are rescued. This provides some anchoring in otherwise marker-free VNTR regions while avoiding false chains in non-VNTR regions. The conservative thresholds (freq=2 only, ≥5 consecutive windows) prevent regressions in cases where rescued k-mers would create ambiguous chains.

#### SA-tag size refinement for coverage-drop DEL calls
When an adaptive-bimodal or flank-gap DEL call is made, the size may be inflated by repeat-unit slippage in the diagonal analysis. If an SA-tag DEL call exists nearby (within 500bp) with ≥2 supporting reads and a compatible size, the SA-tag size replaces the diagonal-based size. SA-tag uses BAM aligner coordinates which handle repeats better than k-mer chaining. Size range: 0.3×–1.5× for weak SA-tag (2-4 reads), 0.3×–5× for strong SA-tag (≥5 reads). Fixes DEL182 (240bp → 183bp), DEL910 (252bp → 859bp).

#### SA-tag out-of-region filtering
SA-tag calls with breakpoints outside the local reference region are discarded. These arise when supplementary alignments map to distant genomic locations, producing bogus calls with absolute genomic coordinates instead of local coordinates.

#### Single-cluster promotion
Per-read DEL/INS clusters with ≥15 reads are promoted to calls even without merging with other clusters. Large SVs (>1000bp) often produce a single cluster where all reads see the same breakpoint.

#### Soft-clip breakpoint evidence
Reads with soft clips ≥20bp at consistent positions (≥3 reads within 5bp) are clustered into breakpoint evidence. Soft-clip sequences are stored for downstream assembly. Left clips indicate right breakpoints; right clips indicate left breakpoints.

#### CIGAR indel calls
Large CIGAR I/D operations (≥30bp) from BAM reads are clustered by position (within 20bp) and type. Clusters with ≥3 reads and size ≥50bp are emitted as calls. This directly detects small-medium SVs that are captured within a single read alignment, bypassing the k-mer chaining pipeline entirely.

#### Soft-clip de Bruijn assembly for INS sizing
At paired soft-clip breakpoints (right-clip and left-clip clusters within 50bp), the clipped sequences are assembled using a greedy de Bruijn graph (k=21). Each side's contig length estimates how far into the insertion the reads extend. The insertion size is estimated from the sum of contig lengths minus overlap. This handles insertions larger than read length where path-based sizing fails.

#### Repetitive region INS call suppression
In highly repetitive regions (≥5 left BPs AND ≥5 right BPs), path-based insertion calls are suppressed unless both breakpoints have very strong overhang support (≥20 reads each). Many breakpoints on both sides indicate a deeply repetitive region where the read graph has many false connections through rescued k-mers, producing artifact insertion calls. Fixes DEL379 which had three false INS calls (264bp, 168bp, 238bp) from rescued k-mer chains.

### Key Files (SV Detection)

| File | Purpose |
|---|---|
| `src/AssemblerSvAnchors.cpp` | `buildSvMSA` — main SV detection engine |
| `src/Sdust.hpp` | Standalone SDUST low-complexity filter (Heng Li's algorithm) |
| `src/Assembler.hpp` | `AlignmentCandidatesInvertedIndexData` with chaining parameters |
| `src/AssemblerInvertedIndex.cpp` | Inverted index build, chaining DP, hit collection |
| `src/InvertedIndexBuilder.hpp` | Count-then-scatter index construction |
| `src/AssemblerMarkers.cpp` | `removeNonUniqueReferenceMarkers` — k-mer blacklisting + adaptive VNTR rescue |
| `srcMain/main.cpp` | svanchors pipeline orchestration (line ~3830+) |

### Test Results (18 cases)

| Case | True Size | Detected | Method |
|------|-----------|----------|--------|
| INS268 | 268bp | 278bp | path-based (1 hop) |
| INS254 | 254bp | 263bp | path-based (2 hops) |
| INS148 | 148bp | 149bp | path-based (1 hop) |
| DEL277 | 277bp | 277bp | split-read (14 reads) |
| DEL347 | 347bp | 328bp | split-read (9 reads) + DUST-STR |
| DEL324 | 324bp | 345bp | diagonal-shift from per-read DEL classifications (was SA-tag only) |
| DEL160 | 160bp | 160bp | split-read (3 reads) |
| DEL137 | 137bp | 96bp | flank-gap DEL (marker rescue improved anchoring) |
| DEL182 | 182bp | 183bp | adaptive-bimodal + SA-tag refinement |
| DEL119a | 119bp | 121bp | adaptive-bimodal (marker rescue resolved VNTR chains) |
| DEL119b | 119bp | 121bp | adaptive-bimodal (marker rescue resolved VNTR chains) |
| DEL147 | 234bp | 229bp | split-read (2 reads) |
| INS454 | 454bp | 400bp | hit-depth drop zone span (was 156bp path-based) |
| INS277 | 277bp | 250bp | hit-depth drop zone span (was 56bp path-based) |

**Undetected (known limitations):**

| Case | Root Cause |
|------|-----------|
| INS57 (57bp) | VNTR, detected as INS 342bp (wrong size, repeat-unit inflation) |
| INS62 (62bp) | 1685bp VNTR (75× core motif), VNTR gap detected but no sizing |
| INS65 (65bp) | 1700bp VNTR, false SA-tag DEL suppressed, no positive call |
| INS235 (235bp) | VNTR, false SA-tag DEL suppressed, no positive call |
| INS61 (61bp) | 2256bp AT-repeat microsatellite |
| INS108 (108bp) | Small diagonal shifts below detection threshold |
| DEL379 (379bp) | 98% repetitive, DEL 239bp (wrong size), false INS calls suppressed |

### Full SBX-D TP Benchmark (80 cases)

80 TP cases from SBX-D.30X.bam (Roche 2×250bp, GRCh38, ~30× coverage), 10 per size bin.
Test cases in `test_cases/sbxd_tp/`. Runner: `test_cases/run_sbxd_tp_tests.sh`.

Scoring: ✅ = correct type, size within 30%. ⚠️ = correct type, size off >30%. ❌ = wrong type or no call.

| Bin | Total | ✅ | ⚠️ | ❌ | ✅% |
|-----|-------|---|---|---|-----|
| DEL <100bp | 10 | 9 | 1 | 0 | 90% |
| DEL 100-500bp | 10 | 9 | 1 | 0 | 90% |
| DEL 500-1000bp | 10 | 7 | 3 | 0 | 70% |
| DEL >1000bp | 10 | 8 | 1 | 1 | 80% |
| INS <100bp | 10 | 9 | 1 | 0 | 90% |
| INS 100-500bp | 10 | 7 | 3 | 0 | 70% |
| INS 500-1000bp | 10 | 7 | 2 | 1 | 70% |
| INS >1000bp | 10 | 0 | 9 | 1 | 0% |
| **TOTAL** | **80** | **56** | **21** | **3** | **70%** |

**Correct type (✅+⚠️): 77/80 = 96%**

Key observations:
- **DEL <100bp at 90%**: compound CIGAR merging (35D+55D→90D) and CIGAR-corroborated k-mer clusters fix tandem repeat and STR cases
- **DEL 100-500bp at 90%**: size-gated CIGAR clustering separates mixed-size deletion clusters (e.g. 57D vs 114D)
- **DEL 500-1000bp at 70%**: SA-tag VNTR suppression relaxed for strong calls (≥10 reads); coverage-drop corroborated k-mer clusters for marker-depleted regions; remaining 3 ⚠️ are tandem repeat copy-number ambiguities (DEL662, DEL605, DEL500)
- **DEL >1000bp at 80%**: 1 ❌ (DEL1432 — no call), 1 ⚠️ (DEL6191 — partial detection)
- **INS <100bp at 90%**: no-covdrop-flip detects INS58 (diagonal-shift DEL with no coverage-drop support → INS)
- **INS 100-500bp at 70%**: reversed-BP INS (INS330), no-covdrop-flip (INS313), flank-gap covDrop sizing (INS310), CIGAR-covdrop corroboration (INS266), SA-DEL-flip (INS306), flank-gap-rounded (INS252)
- **INS 500-1000bp at 70%**: no-covdrop-flip detects INS509 (split-read DEL with no coverage-drop → INS)
- **INS >1000bp** has 0% exact: assembly contigs max out at ~600bp with 250bp reads; flank-gap-alt provides ⚠️ for INS2399

### Q100 DEL 5011-Case Benchmark

5011 DEL cases (≥50bp) from HG002 GRCh38 v5.0q stvar truthset, extracted from SBX-D.30X.bam.
Cases in `/tmp/truthset/del_5000/`. Runner: `/tmp/truthset/run_del5000.sh`.

| Bin | Total | ✅ | ⚠️ | ❌ | ✅% |
|---|---|---|---|---|---|
| DEL <100bp | 1833 | 1556 | 230 | 47 | 84% |
| DEL 100-500bp | 1842 | 1505 | 293 | 44 | 81% |
| DEL 500-1000bp | 508 | 396 | 105 | 7 | 77% |
| DEL >1000bp | 828 | 764 | 44 | 20 | 92% |
| **TOTAL** | **5011** | **4221** | **672** | **118** | **84%** |

**Correct type (✅+⚠️): 4893/5011 = 97%**

Best call sources (✅): CIGAR 1705, cluster 770, diagonal 597, SA-tag 394, adaptive-bimodal 253, flank-gap 178, split-read 131, SDUST-STR 56, path-based 42, SDUST-VNTR 30, VNTR-depth 25, coverage 16, merged-clusters 14, SA-DEL-flip 5, no-covdrop-flip 5.

Remaining ❌ (118 cases): INS-only calls in tandem repeats (66), no calls at all (52 — zero spanning coverage, homozygous deletions, unmappable regions).

Remaining ⚠️ (672 cases): tandem repeat copy-number ambiguity across all sources (CIGAR sees partial repeat units, SA-tag/adaptive-bimodal overshoot by 1-2 repeat units, flank-gap measures one repeat unit shift).

### Q100 DEL 10,160-Case Exhaustive Benchmark

Full 10,160 DEL cases (≥50bp) from HG002 GRCh38 v5.0q stvar truthset.
Cases in `/tmp/truthset/del_5000/` (5011), `/tmp/truthset/del_new_2000/` (1939), `/tmp/truthset/del_remaining/` (3210).
Runner: `/tmp/truthset/run_full_v7.sh`. Per-case output: `/tmp/truthset/v7_per_case/`.

| Metric | Count | % |
|---|---|---|
| ✅ (ratio ≥ 0.7) | 8,883 | 87% |
| ⚠️ (correct type, wrong size) | 1,187 | 12% |
| ❌ (wrong type or no call) | 88 | 1% |
| **Total scored** | **10,158** | |

**⚠️ by best-call source (1,187 total):**
- CIGAR: 369, cluster: 230, bp-pair: 121, SA-tag: 107, adaptive-bimodal: 86
- flank-gap: 72, path-mirror: 54, diagonal: 47, SDUST-STR: 33, split-read: 25
- SDUST-VNTR: 17, VNTR-depth: 10, coverage: 7, merged-clusters: 5, SA-DEL-flip: 2, no-covdrop-flip: 2

**bp-pair/path-mirror ⚠️ analysis (175 cases):**
- All sizes quantized to 50bp multiples (window center distance)
- ~40 in tandem repeats (inherent copy-number ambiguity)
- ~135 in non-TR regions (potentially fixable)
- 58/60 close cases (ratio 0.6-0.7) have breakpoints in wrong windows
- Multiple refinement approaches tried and failed (see Known DEL Failure Patterns #4)

### Known Limitations

1. **bw=100** prevents single-chain deletion detection >100bp; split-read detection handles larger deletions
2. **VNTR insertions** are fundamentally limited: short reads can't span VNTRs, k-mers are non-unique, BAM extraction depletes VNTR reads
3. **Microsatellite insertions** (e.g., AT-repeats): reads are too short to span from unique flanking sequence past the breakpoint
4. **Highly repetitive regions** (>90% non-unique 10-mers): no reliable anchoring possible with k=10

### DEL 500-1000bp Remaining ⚠️ Cases

Three cases remain ⚠️ due to tandem repeat copy-number ambiguity:

- **DEL662**: Coverage-drop 1300bp (marker-depleted), SA-tag 396bp (10 reads). SA-tag undersizes by ~5 repeat units (unit≈53bp). Coverage-drop oversizes. A false BP-pair insertion call at the deletion boundary was suppressing the coverage-drop — fixed by skipping small insertions (<1/3 of coverage-drop size) in the `overlapsInsertion` check. SA-tag refinement now works via widened proximity (coverage-drop region boundaries instead of fixed 500bp) and widened size ratio for marker-depleted regions.

- **DEL605**: medianSpanning=0, no SA-tag (supplementary alignments fall outside extracted region), no coverage-drop detected (requires medianSpanning>5). Only evidence is scattered CIGAR DELs (33-83bp) and HitDepth BPs. Fundamentally limited by lack of spanning reads and SA-tag evidence.

- **DEL500**: Coverage-drop 1650bp, SA-tag 675bp (13 reads, 35% over truth). Adaptive-bimodal picks 121bp (one repeat unit artifact). SA-tag refinement now works via widened proximity and ratio, emitting 675bp. Still ⚠️ because 675/500=1.35 exceeds 30% threshold by 5%.

Code changes (not yet committed):
1. `overlapsInsertion` skip for small insertions relative to coverage-drop size
2. SA-tag proximity uses coverage-drop region boundaries (±300bp) instead of fixed 500bp from center
3. SA-tag size ratio widened to `delSize/bestShift` for strong SA-tag support (≥5 reads) in both flank-gap and adaptive-bimodal paths

### Investigated but Not Viable

- **Relaxed k-mer blacklist** (maxRefKmerCount=50): chains still don't form in VNTRs because k-mers match too many positions; signal-to-noise ratio is too low for the chaining DP
- **Inverted index frequency filters**: downsampling is already disabled for svanchors; the bottleneck is in the chaining DP, not candidate generation
- **VNTR depth estimation from read bases**: fails because BAM extraction depletes VNTR reads (mappers assign low MAPQ)

### VNTR Alignment Strategy Investigation

Extensive investigation into alternative alignment strategies for detecting SVs in tandem repeat regions. All failing cases are VNTR insertions (INS57, INS62, INS65, INS235, INS61) or a highly repetitive deletion (DEL379).

#### Repeat Structure Characterization

| Case | Period | Autocorr | VNTR Length | k=50 Unique |
|------|--------|----------|-------------|-------------|
| INS57 | 57bp | 0.93 | 1950bp | 56% |
| INS62 | 21bp | 0.65 | 1685bp | 54% |
| INS65 | 62bp | 0.98 | 1136bp | 0% (k=80) |
| INS235 | 60bp | 0.96 | 400bp | 8% (k=100) |
| INS61 | 70bp | 0.97 | 218bp | 54% |
| DEL379 | 39bp | 0.67 | 1298bp | 24% |

#### Approaches Tested

1. **k=50 read-to-reference chaining**: 94-96% of reads are chainable with k=50 unique k-mers (vs 7% with k=10). However, the chains don't produce clean SV signals because inserted repeat copies have k-mers that match existing reference positions.

2. **k=50 split-read detection**: Diagonal clustering with k=50 unique k-mers produces 155 split reads for INS57, but the diagonal difference distribution is nearly uniform (no peak at the true 57bp insertion size). Only 3 reads show differences near 57bp. The repeat structure creates false diagonal clusters that dominate the signal.

3. **Reference position gap analysis**: Sorting unique k-mer matches by reference position reveals gaps of ~51bp (close to the 57bp repeat period). However, this gap is a property of the reference's unique k-mer distribution, not the insertion. Both alleles show the same gap pattern.

4. **Read-gap vs reference-gap comparison**: For each read, comparing read-space gaps to reference-space gaps between consecutive unique k-mer matches. The LIS chain absorbs the insertion signal because inserted k-mers match at shifted reference positions.

5. **k=50 read-to-read overlap graph**: True overlaps show 47-196 shared k=50 k-mers vs 7 for non-overlapping reads (20:1 ratio). However, with ~4% sequencing error rate, the expected shared k=50 k-mers between two overlapping reads is only 0.8 (because (1-0.04)^100 ≈ 0.016). The overlap graph has zero paths from left-flank to right-flank reads.

6. **Repeat copy counting**: Per-read repeat unit counts show 0.065 copies/read difference between alleles — too weak to detect with 292 reads.

7. **Full pipeline with k=50**: Running dinara with k=50 produces worse results than k=10 for non-VNTR cases (INS268: 355bp vs 278bp true 268bp) and noisy results for VNTR cases (INS57: 432bp vs true 57bp).

#### Root Cause Analysis

The fundamental limitation is that **VNTR insertions add copies of existing repeat units**. The inserted sequence's k-mers match the reference at positions of existing copies, regardless of k. This makes the insertion invisible to any k-mer-based alignment approach.

The only viable approaches for VNTR insertion detection with short reads require:
- **Paired-end insert size analysis** (not available in current pipeline)
- **Local assembly** with error-tolerant overlap detection (k ≤ 20 for read-to-read matching, but k=20 gives only 19% unique k-mers in the VNTR)
- **Long reads** that can span the entire VNTR

#### Error Rate vs k Trade-off

With ~4% per-base error rate (typical for Illumina short reads):

| k | VNTR Unique | Read-to-Ref Match Rate | Read-to-Read Match Rate | Shared k-mers (100bp overlap) |
|---|-------------|----------------------|------------------------|-------------------------------|
| 10 | 7% | 66% | 44% | 22 |
| 20 | 19% | 44% | 19% | 16 |
| 30 | 32% | 29% | 8% | 6 |
| 50 | 56% | 13% | 2% | 0.8 |

The sweet spot for read-to-reference matching is k=50 (56% unique), but for read-to-read matching it's k=10-20. No single k value works for both.

---

## SV Detection — Latest State (May 2026)

### Latest Commit

`8c9ccd5` on `main` — depth-deficit DEL source (V36k)

V36l changes:
- `src/Assembler.hpp` — `DepthScanDelCall` struct + `depthScanDelCalls()` method declaration
- `src/AssemblerSvAnchors.cpp` — `depthScanDelCalls()` implementation (~250 lines, sweep line depth, minCallSize scaling) + cigar-covdrop source (~40 lines in HitDepth cluster handler) + source priority table (cigar-covdrop=2, early-CIGAR=3, depth-scan-hom/het/sub at 14-16)
- `srcMain/main.cpp` — depth-scan call site after early-CIGAR block, before DP chaining

### Git Info

- **Commit authorship**: `kokyriakidis <kokyriakidis@gmail.com>` — no Co-authored-by trailer
- **Branch**: `main`

### Parameters

- k=10, w=6 (minimizer window; w=6 works best, w=10 causes regressions in INS_130bp and DEL_164bp)
- SDUST: T=20, W=64
- maxK=62 for adaptive multi-k gap filling
- minimap2Bw=100, minimap2MaxGap=5000, chainingMode=1

### SV Detection Pipeline Layers (complete, in order)

**Pre-chaining (in `main.cpp`, run before DP chaining):**

0a. Early-CIGAR DEL detection from BAM CIGAR strings
0b. **Depth-scan DEL detection** — windowed BAM depth scanning (V36l)

**Post-chaining (in `main.cpp` and `AssemblerSvAnchors.cpp`):**

1. Per-read DP chaining (minimap2-sr scoring) → primary + supplementary chains
2. Split-read classification → sv_split_reads.tsv
3. TheseusMSA construction from chain anchors
4. Per-read diagonal analysis → DEL/INS clusters
5. Indirect read alignment via read graph BFS
6. Coverage analysis → breakpoints (chain-endpoint, k-mer hit depth)
7. Breakpoint pair analysis → path-based INS/DEL calls
8. Large-insertion detection via relaxed breakpoints
9. Coverage-drop detection → candidate deletion regions
10. VNTR gap detection + SDUST-gated STR detection (Strategy 1: spanning, Strategy 2: base-count)
11. Adaptive multi-k gap filling (k=62 down to k=10)
12. Pairwise diagonal difference analysis (for non-marker-depleted regions)
13. Flank-gap analysis (for marker-depleted coverage-drop regions)
14. VNTR-depth deletion detection (negative depth → DEL call)
15. Covdrop-indirect insertion detection
16. Split-read deletion detection
17. Merged-cluster calls
18. SA tag parsing from BAM → clustered SV calls
19. Depth-deficit DEL detection (integrated deficit across region)

### Test Case Datasets

#### Dataset 1: Original hand-curated cases (38 cases)
- **Location**: `/tmp/tp_cases/` (10), `/tmp/tp_cases2/` (10), `/tmp/sv_cases/chr*/` (8 INS), `/tmp/sv_cases/DEL_medium_100_500bp/*/` (10 DEL)
- **Also in repo**: `/workspaces/dinara/test_cases/` (6 FN cases)
- **Read length**: ~150bp Illumina
- **Reference**: GRCh38 with `chr` prefix
- **Results**: 22/38 pass (57%), 16 failing
- **No BAM files** for most of these (no SA-tag support)

#### Dataset 2: GIAB HG002 TP cases (30 cases)
- **Location**: `/tmp/giab_cases/` — 30 directories
- **Read length**: 2x250bp Illumina (GIAB HG002)
- **Reference**: GRCh37 (no `chr` prefix, chromosomes named "1", "2")
- **Truth**: GIAB Tier 1 v0.6 SV benchmark
- **Each dir contains**: `reference.fa`, `reads.fa`, `region.bam`, `region.bam.bai`, `info.txt`
- **Results**: 15/30 pass (50%), DEL 12/16 pass, INS 3/14 pass
- **GIAB source data**:
  - VCF: `/tmp/giab/HG002_SVs_Tier1_v0.6.vcf.gz`
  - Reference: `/tmp/giab/ref_chr1.fa`, `/tmp/giab/ref_chr2.fa`
  - BAM: remote `https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.hs37d5.2x250.bam`

#### Dataset 3: SBX-D Roche cases (129 cases) ← PRIMARY EVALUATION SET
- **Location**: `/tmp/sbxd_cases/tp/` (80 TP) and `/tmp/sbxd_cases/fn/` (49 FN)
- **Read length**: ~350-500bp (variable, Roche sequencing)
- **Reference**: GRCh38 with `chr` prefix
- **Truth**: T2T Q100 v1.1 stvar benchmark, benchmarked via truvari
- **Each dir contains**: `reference.fa`, `reads.fa`, `region.bam`, `region.bam.bai`, `info.txt`
- **Directory structure**: `{tp,fn}/{SV_TYPE}_{SIZE_BIN}/{chrN_POS_TYPELEN}/`
  - Size bins: `small_lt100bp`, `medium_100_500bp`, `large_500_1000bp`, `xlarge_gt1000bp`
- **TP breakdown**: 10 per bin × 8 bins (DEL/INS × 4 size categories) = 80
- **FN breakdown**: 49 total (11 DEL_small, 10 DEL_medium, 3 DEL_xlarge, 12 INS_small, 10 INS_medium, 2 INS_large, 1 INS_xlarge)
- **Raw results file**: `/tmp/sbxd_all_results2.txt`
- **Evaluation script**: `/tmp/eval_sbxd2.py` (run: `python3 /tmp/eval_sbxd2.py < /tmp/sbxd_all_results2.txt`)

### SBX-D Results Summary (129 cases)

|  | PASS | WRONG_SIZE | WRONG_TYPE | NO_CALL | Total | Pass% | Correct-type% |
|---|---|---|---|---|---|---|---|
| **TP DEL** | 23 | 12 | 5 | 0 | 40 | 57% | 87% |
| **TP INS** | 8 | 19 | 9 | 4 | 40 | 20% | 67% |
| **FN DEL** | 13 | 2 | 5 | 4 | 24 | 54% | 62% |
| **FN INS** | 3 | 7 | 11 | 4 | 25 | 12% | 40% |
| **ALL** | **47** | **40** | **30** | **12** | **129** | **36%** | **67%** |

Overall by SV type: DEL 36/64 pass (56%), INS 11/65 pass (16%)

### Remote Machine Access

- **SSH**: `ssh kyriakik@ec-hub.sc1.science.roche.com`
- **BAM**: `/sc1/groups/sbx/workspace/kyriakik/data/bams/SBX-D.30X.bam` (23GB, GRCh38, ~350-500bp reads)
- **Reference**: `/sc1/groups/sbx/workspace/kyriakik/data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`
- **Truth VCF**: `/sc1/groups/sbx/workspace/kyriakik/data/truth/GRCh38_HG2-T2TQ100-V1.1_stvar.filt.vcf.gz`
- **Benchmark bed**: `/sc1/groups/sbx/workspace/kyriakik/data/truth/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed`
- **Truvari results**: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/test_runs/bench_20260519_155309/SBX-D/truvari/`
  - `summary.json`: TP=26763, FP=9450, FN=1425, P=0.899, R=0.949, F1=0.923
  - `tp-base.vcf.gz`, `fn.vcf.gz`, `fp.vcf.gz`
- **Extracted cases on remote**: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/cases/SBX-D/{tp,fn}/`
- **Assembly contigs**: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/test_runs/bench_20260519_155309/SBX-D/final_assemblies.fa`
- **SBX tools**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/{assemble,contig_aligner,contig-variant-caller}`
- **Modules**: `SAMtools/1.21-GCC-13.3.0`, `BCFtools/1.21-GCC-13.3.0`, `Micromamba/2.0.7-0`
- **SLURM**: Use `--qos=3h` (not `1h`)

### How to Run Full SBX-D Evaluation

```bash
cd /workspaces/dinara
DINARA=build/Executable/dinara

{
echo "========== SBX-D TP CASES =========="
for d in /tmp/sbxd_cases/tp/*/*/; do
  name=$(basename "$d")
  svtype=$(grep "SV Type:" "$d/info.txt" | awk '{print $NF}')
  svsize=$(grep "SV Length:" "$d/info.txt" | sed 's/[^0-9]//g')
  rm -rf "$d/outdir_eval"
  calls=$($DINARA --command svanchors --reference "$d/reference.fa" --input "$d/reads.fa" --bam "$d/region.bam" --assemblyDirectory "$d/outdir_eval" --Kmers.k 10 --Kmers.minimizerW 6 2>&1 | grep -E ">>> .* CALL" | tail -5)
  echo "--- $name (truth: $svtype ${svsize}bp) ---"
  [ -z "$calls" ] && echo "  (no call)" || echo "$calls"
done
echo "========== SBX-D FN CASES =========="
for d in /tmp/sbxd_cases/fn/*/*/; do
  name=$(basename "$d")
  svtype=$(grep "SV Type:" "$d/info.txt" | awk '{print $NF}')
  svsize=$(grep "SV Length:" "$d/info.txt" | sed 's/[^0-9]//g')
  rm -rf "$d/outdir_eval"
  calls=$($DINARA --command svanchors --reference "$d/reference.fa" --input "$d/reads.fa" --bam "$d/region.bam" --assemblyDirectory "$d/outdir_eval" --Kmers.k 10 --Kmers.minimizerW 6 2>&1 | grep -E ">>> .* CALL" | tail -5)
  echo "--- $name (truth: $svtype ${svsize}bp) ---"
  [ -z "$calls" ] && echo "  (no call)" || echo "$calls"
done
} > /tmp/sbxd_all_results.txt 2>&1

python3 /tmp/eval_sbxd2.py < /tmp/sbxd_all_results.txt
```

### Known INS Failure Patterns (main area for improvement)

1. **Size overestimation for small INS (50-100bp)**: Path traversal overshoots through repetitive sequence → calls 200-700bp instead of 50-100bp
2. **Wrong type for medium INS (200-400bp)**: Coverage-drop from insertion misinterpreted as deletion
3. **Size underestimation for large INS (>500bp)**: Reads can't span full insertion → partial size
4. **NO_CALL for very large INS (>1000bp)**: No breakpoints detected at all
5. **chr1:155188986-155189398 cluster**: 6 FN INS cases at same locus all get DEL from SA-tag — complex region

### Known DEL Failure Patterns

1. **Some small DELs (50-90bp) called as INS**: Breakpoint pair with pathDist > refGap looks like insertion (partially fixed by deletion-like pair suppression)
2. **Large DEL size errors**: SA-tag sometimes gives wrong size for complex regions
3. **XLarge DELs (>5000bp)**: Often only partial detection or no call
4. **bp-pair/path-mirror 50bp quantization**: All bp-pair and path-mirror sizes are multiples of 50bp (window center distance). 175 ⚠️ cases are bp-pair/path-mirror. Of these, ~135 are non-TR. Analysis shows 58/60 close cases (ratio 0.6-0.7) have breakpoints in the **wrong windows**, not just quantization noise. Root causes:
   - **Small DELs (<100bp)**: Chain absorbs deletion within bw=100. No diagonal shifts detected. Per-read analysis shows max drops of ~18bp for 56bp deletions. The chainer treats the deletion as normal gap.
   - **Larger DELs (100-500bp)**: Breakpoints placed in wrong windows. Called size ≠ nearest 50bp multiple to truth.
   - **Approaches tried and failed**: (a) Raw chain endpoint positions — overestimates by variable amount. (b) K-mer diagonal sizing — unreliable, net -5 on full benchmark. (c) Shared-read diagonal differences — chain endpoints have near-zero diagonal shift because chainer absorbs deletion. (d) astarpa2 global alignment — prefers many small edits over one big DEL in repeats.
   - **Fundamental blocker**: bp-pair fires when <2 spanning reads have diagonal shifts. Without spanning reads, no direct measurement of deletion size exists. The window-center distance is the only available estimate.

---

## Full DEL Evaluation — HG002 Q100 v5.0q (May 2026)

### Dataset

- **Truth**: HG002 Q100 v5.0q stvar truthset, 10,173 DEL entries ≥50bp
- **Reads**: Roche 2×250bp variable-length consensus reads, GRCh38, ~30x coverage
- **Matching criterion**: `min(call_size, truth_size) / max(call_size, truth_size) >= 0.7`
- **Cases directory**: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval/cases/` (10,173 directories, each with `reference.fa`, `reads.fa`, `region.bam`, `info.txt`)
- **Eval script**: `analyze.py cases results_<version>`

### Version History

| Version | Key Change | Recall | Avg calls/case |
|---|---|---|---|
| V35 | Baseline (no journey) | 98.7% | 11.0 |
| V36d | K-mer journey: two-diagonal only | 98.9% | 11.4 |
| V36e | K-mer journey: all diagonal pairs | 98.9% | 11.8 |
| V36g | K-mer journey: two-diagonal + triplet | 98.9% | 11.5 |
| V36h | + removed MAPQ0 filter | 99.1% | 11.6 |
| V36i | + removed depth-fc filter | 99.2% | 11.9 |
| V36j | + source-priority dedup + SA-DEL removal | 99.2% | 11.8 |
| V36k | + depth-deficit DEL source | 99.4% | 13.9 |
| V36l | + depth-scan DEL source (windowed BAM depth) | 99.4%* | ~15 |
| V36m | + cigar-covdrop + sweep line + minCallSize + skip MSA + ref-only chaining | 99.4%* | ~15 |

*V36l/V36m recall on HG002 is unchanged (depth-scan adds calls but doesn't recover new TPs for germline). The depth-scan and cigar-covdrop sources are primarily impactful for somatic analysis (HG008-T: 81.7% → 100%).

### Current State (V36m)

- **Germline recall (HG002)**: 99.4% (10,112 / 10,173)
- **Somatic recall (HG008-T)**: 100% (60 / 60) at pctsize=0.0, 100% (60 / 60) at pctsize=0.7
- **Somatic precision (HG008-T)**: 85.7% at pctsize=0.0
- **61 remaining germline misses**: Small DELs in large tandem repeat arrays where depth deficit is dominated by mappability loss. 0 cases with no calls at all. 100% detection sensitivity.
- **0 remaining somatic misses**: All 60 PASS DELs detected.
- **Avg calls/case**: ~15 (germline), ~20 (somatic, higher coverage)

### Performance (V36m)

Three optimizations reduced per-case runtime by 4-25×:

1. **Skip Theseus MSA** — not used for DEL calling (saves 19-146s per case). Commented out, re-enable for INS calling.
2. **Read-vs-reference chaining only** — `referenceReadCount` set to actual count instead of 0 (all-vs-all). Skips O(n²) read-vs-read pairs that are only needed for Theseus MSA. Saves 40-75% of chaining time.
3. **Sweep line depth computation** — O(n_reads + region_length) instead of O(n_reads × read_length). Critical for the chr13 6.3Mb case (7.3M reads).

| Case | Before | After | Speedup |
|------|--------|-------|---------|
| chr16 (168bp DEL, 12K reads) | 296s | 12s | 25× |
| chr20 (44kb DEL, 63K reads) | 556s | 137s | 4× |
| chr13 (6.3Mb DEL, 7.3M reads) | timeout | ~40min | now completes |

**Profiled breakdown (chr20, 63K reads):**

| Step | Time |
|------|------|
| Read loading | 3.4s |
| Markers + index | 4.1s |
| parseBamEvidence | 1.5s |
| depthScanDelCalls | 0.03s |
| DP chaining | 131s |
| Split-read classification | 0.2s |
| **Total** | **137s** |

**Remaining bottleneck**: DP chaining (131s for chr20). This is inherent — all 63K reads share markers with the same reference, so every read is a candidate pair. Running 60 cases in parallel, wall time is ~2-3 minutes.

### Active DEL Sources

**V36k numbers (with source-priority dedup + depth-deficit):**

| Source | Calls | TPs | Precision |
|---|---|---|---|
| merged-clusters | 3,519 | 3,128 | 88.9% |
| diagonal | 2,110 | 1,745 | 82.7% |
| early-CIGAR | 12,499 | 6,070 | 48.6% |
| SA-tag | 3,147 | 1,396 | 44.4% |
| split-read | 642 | 343 | 53.4% |
| cluster | 5,340 | 1,299 | 24.3% |
| flank-gap | 2,548 | 797 | 31.3% |
| kmer-journey | 10,882 | 1,965 | 18.1% |
| per-read-DEL | 1,614 | 548 | 34.0% |
| path-based | 664 | 146 | 22.0% |
| INV-cluster | 5,210 | 862 | 16.5% |
| depth-deficit-hom | 17,072 | 2,370 | 13.9% |
| depth-deficit-het | 16,701 | 1,286 | 7.7% |
| multi-k | 59,285 | 6,550 | 11.0% |

**Source accuracy when disagreeing with dedup winner (V36i analysis, pre-source-priority):**

| Source | Right | Wrong | Accuracy |
|---|---|---|---|
| merged-clusters | 916 | 10 | 98.9% |
| diagonal | 506 | 27 | 94.9% |
| early-CIGAR | 3,605 | 495 | 87.9% |
| SA-tag | 830 | 115 | 87.8% |
| split-read | 16 | 3 | 84.2% |
| cluster | 693 | 323 | 68.2% |
| flank-gap | 396 | 183 | 68.4% |
| kmer-journey | 1,049 | 566 | 65.0% |
| per-read-DEL | 427 | 274 | 60.9% |
| multi-k | 99 | 908 | 9.8% |

This accuracy metric drove the source-priority ranking in the dedup step.

### Removed Sources/Filters

- **SA-DEL** (commit `7e340a3`): 18 calls, 2 TPs, 0 unique TPs. Generated in two early-exit paths (lines 1224, 1293) via `"SA-" + sc.svType` when region lacked sufficient reference markers for MSA. Removed — all TPs covered by other sources.
- **MAPQ0 filter** (V36h): Removed. Was filtering calls where >80% of reads had MAPQ=0. Removal gained +25 TPs, 0 regressions.
- **Depth fold-change filter** (V36i): Removed. Was filtering multi-k calls ≥300bp with DHFFC ≥0.95. Removal gained +5 TPs (het DELs where other haplotype fills depth), 0 regressions.

### Depth-Deficit DEL Source (V36k)

Infers deletion size from the integrated depth deficit across the entire region, without needing to identify repeat units or localize breakpoints.

**Algorithm:**
1. Compute per-base depth from BAM across the full region
2. Measure flanking depth (median of left/right edges)
3. `deficit = flank_depth × region_length - sum(all_depths)`
4. `hom_del_size = deficit / flank_depth`
5. `het_del_size = deficit / (flank_depth / 2)`
6. Emit both as separate calls with source "depth-deficit-hom" / "depth-deficit-het"

**Why it works:** In tandem repeat regions, k-mer-based sources produce wrong sizes because anchors land on wrong repeat copies. But the total integrated depth across the region still reflects the missing bases — a deletion means fewer total bases mapped, regardless of where individual reads land.

**Tries 4 flank sizes** (200, 300, 500, 800bp) to handle varying region sizes. Each flank size produces independent hom/het calls.

**Impact:** +23 TPs (99.2% → 99.4%). Largest gain in 1000+bp bin (98.6% → 99.9%). Adds ~33K calls with low precision (7.7-13.9%), acceptable since downstream ML model handles classification.

**Limitations:** Fails for small DELs in large repeat arrays where mappability loss across the entire repeat dominates the deficit signal, swamping the small deletion contribution.

### Source-Priority Dedup (V36j, commit `c9d5f94`)

**Problem**: The old dedup sorted calls by read count descending. multi-k always won (most reads) but was correct only 9.8% of the time when disagreeing with other sources. In 5,791 of 10,173 loci, multi-k was the wrong-size winner with a correct-size call from another source available.

**Solution**: Sort by source accuracy instead of read count. Priority order: merged-clusters > diagonal > cigar-covdrop > early-CIGAR > SA-tag > split-read > cluster > flank-gap > kmer-journey > per-read-DEL > INV-cluster > path-based > depth-deficit-hom > depth-deficit-het > depth-scan-hom > depth-scan-het > depth-scan-sub > multi-k. Read count is tiebreaker within same source.

**Impact on recall**: None — recall is 99.2% for both V36i and V36j because the eval counts ALL emitted calls, not just the top one.

**Impact on top-N call selection** (relevant for downstream ML model):

| Top-N | V36i (read-count) | V36j (source-priority) |
|---|---|---|
| Top-1 | 51.1% | 64.2% |
| Top-2 | 63.7% | 77.7% |
| Top-3 | 71.5% | 83.7% |
| Top-5 | 82.6% | 90.1% |
| All | 99.2% | 99.2% |

### Multi-Source Overlap Analysis

**96.6% of loci** (9,826/10,173) have calls from multiple sources. Source count distribution:

| Sources | Loci |
|---|---|
| 1 | 347 |
| 2 | 1,829 |
| 3 | 3,481 |
| 4 | 2,736 |
| 5 | 1,299 |
| 6 | 408 |
| 7 | 66 |
| 8 | 7 |

### Call Selection Strategy Analysis

Tested different strategies for ranking the ~12 calls per locus to surface the correct one. All strategies operate on the same set of emitted calls (99.2% recall).

| Strategy | Top-1 | Top-3 | Top-5 |
|---|---|---|---|
| by-reads (old V36i) | 22.0% | 52.3% | 73.3% |
| source-only (V36j) | 61.9% | 79.3% | 87.2% |
| tier+dhffc | 63.0% | 80.2% | 88.4% |
| **dedup(tier+dhffc)+rank** | **63.0%** | **82.1%** | **90.2%** |
| dhffc-only | 27.7% | 68.4% | 87.1% |
| oracle | 99.2% | 99.2% | 99.2% |

**Best strategy: `dedup(tier+dhffc)+rank`** — dedup using source tiers + dhffc, then rank survivors the same way.

**Source tiers** (grouping by accuracy class instead of 12 strict ranks):

| Tier | Sources | Accuracy |
|---|---|---|
| 0 | merged-clusters, diagonal | >90% |
| 1 | early-CIGAR, SA-tag, split-read | 45-88% |
| 2 | cluster, flank-gap, kmer-journey, per-read-DEL, INV-cluster, path-based | 25-68% |
| 3 | multi-k | <15% |

Within same tier, DHFFC (depth fold-change) decides — lower dhffc = stronger depth evidence = better.

**Per size bin (best strategy vs source-only):**

| Size Bin | source-only Top-1 | dedup(tier+dhffc) Top-1 | source-only Top-5 | dedup(tier+dhffc) Top-5 |
|---|---|---|---|---|
| 50-100bp | 70.6% | 69.4% | 89.3% | 94.9% |
| 100-500bp | 55.6% | 58.5% | 85.8% | 89.3% |
| 500-1000bp | 34.7% | 44.3% | 75.3% | 70.0% |
| 1000+bp | 70.1% | 67.6% | 91.4% | 84.0% |

### Gap Analysis (37% where correct call isn't top-1)

3,790 cases where source-priority picks the wrong call as #1 but a correct call exists:

- **early-CIGAR is the main blocker** (63.6% of gap cases) — high priority but often wrong size
- **multi-k has the correct call** in 49% of gap cases — produces many calls, one often matches truth
- **DHFFC distinguishes correct from wrong** in 63.6% of gap cases (correct call has lower dhffc). In 54.3%, the dhffc gap is >0.2 (strong signal).
- **Read count is uninformative**: correct call has more reads only 49.3% of the time (coin flip)
- **Correct call rank distribution**: 30% at rank 2, 47% within top-3, 68% within top-5, 94% within top-10

### Downstream ML Model Consideration

The downstream ML model classifies calls. For maximum sensitivity, report ALL calls (99.2% recall) and let the ML model select. The source-priority dedup and tier+dhffc ranking are useful if the ML model uses emission order or rank as a feature, but don't affect the set of calls available.

### K-mer Journey Implementation Details

Located in `src/AssemblerSvAnchors.cpp` (~lines 6842-7150).

**Approach**: Multi-resolution k-mer scanning (k=60,50,40,30,20,14) builds per-read anchor positions against the reference. Two methods extract deletion signals:

1. **Two-diagonal**: Cluster anchor diagonals (refPos-readPos), find two most populated clusters, difference = DEL size. Best precision (30%).
2. **Triplet**: For anchor pairs A,C with transition-zone anchor B between them, compute refGap-readGap across A-C. Adds ~255 TPs but also ~1,027 FPs.

**Key implementation detail**: `covered[]` array marks only the start position of each k-mer (not the full k-mer span), allowing denser anchors at multiple k values.

**Votes are clustered by median** for robustness to outlier votes.

### Tandem Repeat Limitation

The 61 remaining misses (down from 84 after depth-deficit source) are small DELs in large tandem repeat arrays where:
- The repeat unit is shorter than the read length
- Every read fits entirely within the repeat — no unique-sequence flanking anchors exist
- K-mer anchors match repeat copies at wrong reference positions
- All sequence-based sources produce wrong sizes (multi-k harmonics)
- The depth-deficit approach also fails because the deletion is tiny relative to the repeat-induced mappability loss across the whole region (e.g., 67bp DEL in a 4kb repeat array)
- 10 of 13 analyzed strong cases had NO reads with the correct diagonal pair
- This is a fundamental limitation of short-read sequencing in tandem repeats

**Size distribution of remaining 61 misses:** 43 are 50-100bp, 11 are 100-500bp, 7 are 1000+bp.

**Depth-deficit recovered 23 of the original 84:** Mostly larger DELs (1000+bp) and cases where the repeat array is small enough that the deletion contributes a measurable fraction of the total depth deficit.

### Depth-Scan DEL Source (V36l, improved in V36m)

Windowed BAM depth scanning to detect deletions from depth drops. Unlike depth-deficit (which integrates total deficit across the whole region), depth-scan identifies the specific contiguous region of low depth and estimates DEL size from its width.

**Algorithm (two-pass detection):**
1. Open BAM, compute per-base depth via **sweep line** (O(n_reads + region_length)) — record +1 at each aligned segment start and -1 at end, then prefix sum
2. Aggregate into adaptive windows (50bp for regions <3kb, 100bp for 3-10kb, 200bp for 10-50kb, 500bp for >50kb)
3. **Pass 1 — Sliding flanking:** For each position, compute left-flanking median depth (up to 10 windows back). Detect contiguous runs of windows where `depth / flank_depth ≤ threshold`. Good for small/medium DELs where local context is reliable.
4. **Pass 2 — Edge-based flanking:** Use median depth of first/last N windows of the region as a stable baseline. Detect drops relative to this global baseline with **gap tolerance** (up to `max(3, nWindows/50)` consecutive windows above threshold allowed within a drop). Handles large DELs where the sliding flanking gets contaminated by the deletion itself.
5. Both passes scan at 3 sensitivity levels (thresholds: 0.65, 0.75, 0.85)
6. **Minimum call size** scales with region length: `max(50bp, regionLength/500)`. Suppresses noise from natural depth variation in large regions (e.g., 100+ false 1kb calls in a 6.3Mb region).
7. Deduplication: a new call overlapping an existing call is kept if it is >1.5× larger (better estimate supersedes partial detection)
8. Classify by depth ratio: `<0.05` → hom, `<0.55` → het, else → sub (subclonal)
9. Emit calls as `depth-scan-hom` / `depth-scan-het` / `depth-scan-sub`

**Key design decision:** Runs in `main.cpp` **before DP chaining** (right after early-CIGAR), so calls are emitted even when assembly times out on high-read-count regions. This is the primary motivation — large DELs (5kb+) often have too many reads for the DP chaining to complete within wall time.

**Implementation:** `Assembler::depthScanDelCalls()` in `src/AssemblerSvAnchors.cpp` (~250 lines). Uses htslib to read BAM directly. Called from `srcMain/main.cpp` after the early-CIGAR block.

**V36m improvements:**
- **Sweep line depth computation** replaces per-base CIGAR walk. The old approach was O(n_reads × read_length) — for the chr13 6.3Mb case (7.3M reads × 250bp), that's ~1.8 billion individual depth increments. The sweep line is O(n_reads + region_length).
- **Minimum call size scaling** (`regionLength/500`) prevents false positives from natural depth variation in large regions.

**Why two passes:** The sliding-flanking pass (Pass 1) works well for DELs up to ~20kb but fails for larger ones. For example, a 64kb het DEL at 34% VAF has ~130 windows of low depth. After scanning ~50 windows into the deletion, the sliding flanking estimate (10 windows back) is entirely within the deletion, so it stops detecting a drop. The edge-based pass (Pass 2) uses the region edges as baseline, which stays stable regardless of deletion size. Gap tolerance handles local depth spikes inside large deletions (e.g., 4 windows spiked above threshold inside a 64kb deletion due to mappability artifacts).

**Somatic impact:** Recovers large DELs (5kb-6.3Mb) that assembly cannot handle within wall time. See HG008-T evaluation below.

### CIGAR-Guided Covdrop Source (V36m)

When a HitDepth cluster (k-mer depth drop) has no spanning chains — indicating a marker-depleted tandem repeat where chains can't cross the deletion — the code searches `cigarIndels` for a DEL call within 100-1500bp of the HitDepth position. If found, it emits a call at the HitDepth position (correct) using the CIGAR size (correct).

**Motivation:** In tandem repeats, the aligner can place CIGAR D operations at a different repeat copy than where the deletion actually occurs. The early-CIGAR source finds the right size but wrong position. The HitDepth cluster finds the right position but has no size estimate (no spanning chains). Combining them produces a correct call.

**Example:** chr3:169770590 (108bp DEL, VAF=10%) — the aligner placed the CIGAR D 672bp downstream. The HitDepth cluster at offset 2053 was 53bp from truth. The cigar-covdrop call (89bp at offset 2053) matched truvari with refdist=500.

**Priority:** cigar-covdrop gets priority 2 (above early-CIGAR at 3) so the position-corrected call wins dedup over the misplaced original.

**Implementation:** ~40 lines in the HitDepth cluster handler in `src/AssemblerSvAnchors.cpp`.

### Dinara V36k vs sbx-assemble — DEL Comparison (HG002)

Direct comparison using identical truvari parameters on HG002 Q100 v5.0q stvar truthset (10,175 DEL entries ≥50bp in benchmark BED regions).

**Truvari parameters (identical for both tools):**
```
refdist=2000, pctseq=0.0, pctovl=0.0, sizemin=50, sizefilt=30, sizemax=50000
passonly, pick=multi, dup-to-ins, no-ref=a, chunksize=5000
```

**With `pctsize=0.0` (no size matching — sbx-assemble's default params):**

| Metric | Dinara V36k | sbx-assemble | Δ |
|---|---|---|---|
| TP-base | **10,175** | 9,750 | **+425** |
| FN | **0** | 425 | **−425** |
| Recall | **100.0%** | 95.8% | **+4.2pp** |
| Total DEL calls | 128,484 | 76,339 | — |

**With `pctsize=0.7` (size accuracy required — ≥70% reciprocal size match):**

| Metric | Dinara V36k | sbx-assemble | Δ |
|---|---|---|---|
| TP-base | **10,121** | 9,163 | **+958** |
| FN | **54** | 1,012 | **−958** |
| Recall | **99.5%** | 90.1% | **+9.4pp** |

**Notes:**
- sbx-assemble's ML filtering does not remove any DELs — the unfiltered and filtered DEL counts are identical (76,339). The comparison is already apples-to-apples.
- Dinara emits all candidate calls from all sources (128K), not a single filtered call per locus. Precision comparison is not meaningful until final call selection is applied.
- The truvari pctsize=0.0 result (100%) means dinara produces at least one call near every truth DEL, but not always the right size. The internal eval with 70% reciprocal overlap shows 99.4%.

**Cluster paths:**
- Dinara VCF: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval/dinara_v36k_dels.vcf.gz`
- Dinara truvari (pctsize=0.0): `.../full_del_eval/truvari_v36k/`
- Dinara truvari (pctsize=0.7): `.../full_del_eval/truvari_v36k_pctsize07/`
- sbx-assemble truvari (pctsize=0.0): `.../test_runs/bench_20260519_155309/SBX-D/truvari/`
- sbx-assemble truvari (pctsize=0.7): `.../test_runs/bench_20260519_155309/SBX-D/truvari_pctsize07/`
- VCF generation script: `.../full_del_eval/logs_to_vcf.py`

---

## Somatic DEL Evaluation — HG008-T V0.5 (June 2026)

### Dataset

- **Truth**: GIAB HG008-T V0.5 draft somatic SV benchmark, 60 PASS DELs ≥50bp
- **BED regions**: `GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed`
- **BAM**: `/sc1/groups/sbx/workspace/eizengaj/structural_variants/somatic_mapping/dedupped_sorted/HG008_T.intra-consensus.personalized.dedup.sorted.bam`
- **Reference**: GRCh38 (`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`)
- **GIAB reference** (for truvari): `GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta`
- **Coverage**: ~50-80x
- **VAF range**: 5% to 100% (somatic, includes subclonal)
- **DEL size range**: 54bp to 6.3Mb

### Truth Set Details

The ALL VCF has 72 DELs total:
- 60 PASS (benchmark set)
- 8 VAFlt5percent (VAF <5%, below detection threshold)
- 4 LT50 (below 50bp size threshold)

### Results

| Version | Binary | Recall | TP | FN |
|---|---|---|---|---|
| V36k (assembly only) | `dinara_v36k_depthdeficit` | 81.7% | 49 | 11 |
| V36l v1 (+ depth-scan, sliding only) | `dinara_v36k_depthscan` | 95.0% | 57 | 3 |
| V36l v2 (+ edge-based + gap tolerance) | `dinara_v36l_depthscan2` | 96.7% | 58 | 2 |
| V36m (+ cigar-covdrop + sweep line + minCallSize) | `dinara_v36m_cigarcovdrop` | **100%** | **60** | **0** |

### Failure Analysis (resolved in V36m)

**7 cases timed out** during DP chaining (assembly never completed):
- All had large regions with many reads (5.8kb-64kb DELs, 50K-262K reads)
- All had clear depth signals (depth ratio 0.00-0.70)
- The depth-scan source recovers all 7 because it runs before DP chaining

**3 cases had wrong-size calls** from assembly:
- chr17:31977413 (27kb, VAF=40%) — assembly found 31bp, depth-scan found 27kb ✅
- chr7:6404812 (1.2kb, VAF=28%) — assembly found 43bp, depth-scan found 1.2kb ✅
- chr7:6408188 (62bp, VAF=26%) — assembly found 35+36bp split calls. Depth-scan can't resolve 62bp (below window granularity). Truvari matched via other means ✅

**chr11:58991970 (64kb, VAF=34%)** — Initially missed by depth-scan v1 because the sliding flanking got contaminated 25kb into the deletion. Fixed in v2 with edge-based flanking + gap tolerance (4 windows spiked above threshold inside the deletion). Now detected as 65kb call (ratio=0.99). ✅

**chr3:169770590 (108bp, VAF=10%)** — In a tandem repeat region, the aligner placed the CIGAR D 672bp downstream at a different repeat copy. The HitDepth cluster correctly identified the depth drop at the true position but had no spanning chains (marker-depleted repeat). Fixed in V36m with cigar-covdrop: combines the CIGAR size (89bp) with the HitDepth position (53bp from truth). ✅

**chr13:66842467 (6.3Mb)** — Region has 7.3M reads. Previously the depth-scan BAM iteration used per-base CIGAR walk (O(n_reads × read_length) = ~1.8 billion increments), causing timeout. Fixed in V36m with sweep line depth computation (O(n_reads + region_length)). Also added minCallSize scaling (regionLength/500) to suppress ~100 false 1kb calls from natural depth variation in the 6.3Mb region. Now detected as 6.3Mb call (dhffc=0.47). ✅

### Somatic-Specific Observations

- Dinara detects somatic DELs down to **VAF ~7%** (chr6:72397103 at 6.6%, chr10:132826517 at 8.8%)
- The depth-scan source is essential for somatic analysis — large somatic DELs (5kb+) consistently time out during assembly
- The depth-scan `sub` (subclonal) threshold (ratio <0.85) catches het-like and subclonal depth drops that the existing depth-deficit source (which uses fixed flank sizes) may miss
- The cigar-covdrop source handles tandem repeat position misplacement by combining CIGAR size with HitDepth position

### Cluster Paths

- **Cases**: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/somatic_hg008/cases/`
- **Results (assembly only)**: `.../somatic_hg008/results/`
- **Results (depth-scan v1)**: `.../somatic_hg008/results_dscan/`
- **Results (depth-scan v2)**: `.../somatic_hg008/results_dscan2/`
- **Results (V36m)**: `.../somatic_hg008/results_cigarcovdrop/`
- **Truth set**: `.../somatic_hg008/truth/`
- **Manifest**: `.../somatic_hg008/manifest.tsv`
- **VCFs**: `.../somatic_hg008/dinara_somatic_dels.vcf.gz` (assembly only), `.../somatic_hg008/dinara_dscan2_dels.vcf.gz` (depth-scan v2), `.../somatic_hg008/dinara_cigarcovdrop_dels.vcf.gz` (V36m)
- **Truvari**: `.../somatic_hg008/truvari_all/` (assembly only), `.../somatic_hg008/truvari_dscan2/` (depth-scan v2)
- **Binary (v1)**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36k_depthscan`
- **Binary (v2)**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36l_depthscan2`
- **Binary (V36m)**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36m_cigarcovdrop`
- **Prepare script**: `.../somatic_hg008/prepare_cases_v2.sh`
- **Run script**: `.../somatic_hg008/run_depthscan2.sh`

### Truvari Command (GIAB-recommended for HG008-T)

Two configurations recommended by the GIAB README (V0.5):

**Position-only matching (pctsize=0):**
```bash
truvari bench \
  -b GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz \
  -c {calls}.vcf.gz \
  --reference GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
  --includebed GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed \
  --refdist 1000 --pctseq 0 --pctsize 0 --pctovl 0 \
  --sizemax -1 --passonly --pick multi --typeignore \
  -o {output_ps0}
```

**Size-aware matching (pctsize=0.7):**
```bash
truvari bench \
  -b GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz \
  -c {calls}.vcf.gz \
  --reference GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
  --includebed GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed \
  --refdist 1000 --pctseq 0 --pctsize 0.7 --pctovl 0 \
  --sizemax -1 --passonly --pick multi --typeignore \
  -o {output_ps07}
```

Source: GIAB HG008-T V0.5 README — `--refdist 1000`, `--typeignore`, `--pick multi`, `--sizemax -1` (no upper limit), `--pctseq 0` (no sequence similarity). The README recommends two runs: pctsize=0 (position-only) and pctsize=0.7 (standard size matching).

---

### Cluster Paths and Binaries

- **Binaries**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36{d..k}_*`, `dinara_v36k_depthscan`, `dinara_v36l_depthscan2`
- **Results**: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval/results_v36{d..k}/`
- **Somatic results**: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/somatic_hg008/`
- **Analysis scripts**: `analyze.py`, `source_essentiality.py`, `analyze_multisource_v2.py`, `analyze_gap.py`, `analyze_adaptive.py` in `/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval/`
- **SLURM template**: `run_eval_v36k.sh` — array job, 100 cases per task, 102 tasks

### Build and Deploy

```bash
# Build locally
cd /workspaces/dinara/build_release && make -j$(nproc) dinaraExecutable

# Deploy to cluster
scp build_release/Executable/dinara kyriakik@ec-hub.sc1.science.roche.com:/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_<version>

# Run eval (germline HG002 DEL)
ssh kyriakik@ec-hub.sc1.science.roche.com
cd /sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval
sbatch run_eval_<version>.sh

# Run eval (germline HG002 INS)
cd /sc1/groups/sbx/workspace/kyriakik/structural_variants/full_ins_eval
sbatch run_eval_<version>.sh

# Run eval (somatic HG008-T)
cd /sc1/groups/sbx/workspace/kyriakik/structural_variants/somatic_hg008
sbatch run_depthscan.sh

# Analyze
python3 analyze.py cases results_<version>
```

## Full INS Evaluation — HG002 Q100 v5.0q (June 2026)

### Dataset

- **Truth**: HG002 Q100 v5.0q stvar truthset, INS ≥50bp (28,942 entries, 26,101 unique cases after dedup by chr_pos_size)
- **BAM**: SBX-D.30X.bam (Roche 2×250bp variable-length consensus reads, GRCh38)
- **Benchmark BED**: GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
- **Region extraction**: ±2000bp around each truth position → reference.fa, reads.fa, region.bam per case

### Truvari Parameters

See [Standardized Truvari Benchmark Parameters](#standardized-truvari-benchmark-parameters) for the GIAB-recommended commands.

Parameters used for INS evaluation: `refdist=2000, chunksize=5000, sizemax=50000, passonly, pick=multi, pctseq=0, pctsize={0, 0.7}`.

Note: The GIAB v5.0q README recommends `--pick ac` and `--refine`. Our runs use `--pick multi` and no `--refine` for simpler FN analysis. Results are comparable for recall measurement.

### Version History

| Version | Commit | INS recall (p=0) | INS FN | Key changes |
|---------|--------|-----------------|--------|-------------|
| V36o | — | ~52% (pilot) | — | Initial INS calling, read-graph + soft-clip + reversed-BP + large-ins + hit-depth + flank-gap + covdrop-indirect + CIGAR-covdrop |
| V36p | `9dc0282` | ~78% (pilot) | — | Added early-CIGAR INS calls (≥3 reads, ≥30bp) and indirect-covdrop fallback |
| V36q | `91e2134` | 96.6% | 606 | Widened soft-clip pairing gap [-200,200], improved indirect-covdrop position (soft-clip midpoint), added softclip-unpaired fallback |
| V36r | `66ce950` | 98.33% | 353 | Relaxed indirect-covdrop from bestStrong==nullptr to insertionCallRegions.empty(), added zero-anchor soft-clip INS detection |
| V36s | `4b396d0` | **98.64%** | **287** | Fixed bad-segment soft-clip INS detection, fixed parseBamEvidence genomic coordinate leak |

### Current State (V36s)

**Truvari results (stvar v5.0q truth, all cases extracted):**

Truth set: HG002 Q100 v5.0q stvar — 28,942 INS ≥50bp (26,101 unique positions), 17,562 DEL ≥50bp (15,472 unique positions).
Truvari params: `refdist=2000, chunksize=5000, sizemax=50000, passonly, pick=multi, pctseq=0`.

| | pctsize=0 | pctsize=0.7 |
|---|---|---|
| INS TP-base | 20,793 | 10,668 |
| INS FN | 287 | 10,412 |
| INS recall | **98.64%** | 50.61% |
| INS precision | 54.54% | 16.34% |
| DEL TP-base | 11,936 | 11,811 |
| DEL FN | 2 | 127 |
| DEL recall | **99.98%** | **98.94%** |
| DEL precision | 79.77% | 21.89% |

Note: INS truth has 21,080 entries within truvari's matching scope at pctsize=0 (the rest are outside sizemax or unmatched). DEL truth has 11,938 entries in scope.

### Active INS Sources (in order of call volume)

| Source | Calls | Description |
|--------|-------|-------------|
| soft-clip | 912,000+ | Paired soft-clip de Bruijn assembly; widened gap [-200,200] |
| early-CIGAR | 15,800+ | CIGAR-based INS ≥30bp with ≥3 reads |
| indirect-covdrop | 15,200+ | Fallback: indirect read bases / coverage when no INS call exists |
| read-graph | 10,800+ | BFS path measurement through read graph |
| reversed-BP | 2,900+ | Reversed breakpoint orientation detection |
| CIGAR-covdrop | 1,800+ | CIGAR INS corroborated by coverage drop |
| large-ins | 1,000+ | Strong breakpoint with partner, path-based sizing |
| large-ins-single | 1,000+ | Single strong breakpoint, no partner |
| large-ins-het | 850+ | Heterozygous large insertion |
| large-ins-single-het | 900+ | Heterozygous single-BP large insertion |
| softclip-unpaired | 580+ | Unpaired soft-clip fallback when no other INS source fires |
| flank-gap | 400 | Flank gap measurement |
| flank-gap-rounded | 390 | Rounded flank gap |
| hit-depth | 110 | Hit-depth-only breakpoint |
| covdrop-indirect | 9 | Coverage drop with indirect reads |

### INS Detection Architecture

1. **parseBamEvidence** — extracts soft-clip breakpoints and CIGAR indels from BAM reads
2. **early-CIGAR** — emits INS calls from CIGAR clusters (≥3 reads, ≥30bp) before MSA
3. **Zero-anchor early exit** — when allRefOrdinals < 2, checks paired/unpaired soft-clips before skipping MSA
4. **MSA alignment + read classification** — classifies reads as REF/INS/DEL/UNK via diagonal analysis
5. **Read-graph BFS** — finds insertion-internal reads via graph traversal
6. **Breakpoint analysis** — coverage-drop, hit-depth, VNTR detection
7. **Path-based INS sizing** — BFS paths through read graph for size estimation
8. **indirect-covdrop fallback** — when no INS call exists but ≥10 indirect reads, estimates size from indirect bases / coverage; position from soft-clip midpoint > HitDepth drop > region center
9. **softclip-unpaired fallback** — when allInsCalls is empty and strong soft-clip clusters exist

### Key Design Decisions

- **All-vs-all chaining** (`referenceReadCount = 0`): Required for INS read graph. Read-vs-ref-only chaining (used for DEL-only V36n) misses read-to-read connections needed for BFS path measurement.
- **Soft-clip pairing gap [-200, 200]**: Widened from [-10, 50] to catch cases where L/R clips are further apart or in unexpected order due to insertion complexity.
- **indirect-covdrop position hierarchy**: soft-clip midpoint > HitDepth drop > region center. Soft-clip positions are base-precise; HitDepth drops can be hundreds of bp off.
- **insertionCallRegions.empty() guard**: The indirect-covdrop fallback fires whenever no INS call was emitted, regardless of whether a strong breakpoint was found. Previously, a false-positive strong BP at the region edge would block the fallback.

### Remaining 183 FN (pctsize=0)

Primarily cases with no detectable evidence in the reads — no soft-clips, no CIGAR indels, no indirect reads. The insertion sequence doesn't share k-mers with the reference, and reads spanning the insertion don't produce alignment artifacts detectable by the current pipeline.

### Size Estimation Gap

Recall drops from 99.0% → 52.0% when requiring size match within 30% (pctsize=0.7). Root causes:
- **soft-clip**: clip length underestimates true insertion size (only captures the portion of the insertion that extends beyond the read alignment)
- **indirect-covdrop**: size = indirect_bases / coverage is a rough estimate
- **read-graph BFS**: limited by 250bp read length; can't span insertions >~500bp

Improving size accuracy would require assembled contigs spanning the insertion, which is fundamentally limited by 2×250bp read lengths for large insertions.

### HG002 Germline SV Benchmark

**Truth set**: HG002 Q100 v5.0q stvar, GRCh38.
- INS ≥50bp: 28,942 entries (23,086 unique positions, 26,101 unique case names)
- DEL ≥50bp: 17,562 entries (15,319 unique positions, 15,472 unique case names)
- Truvari in-scope (after sizemax=50000 and pick=single dedup): 21,080 INS, 11,938 DEL
- Duplicate case names (same chrom+pos+size, different het alleles): 1,400 INS, 2,090 DEL — these share the same extracted case directory since both alleles are at the same position.

**V36u results** (old truvari params: refdist=500, pctseq=0, no includebed):

| Metric | pctsize=0 | pctsize=0.7 |
|---|---|---|
| **INS recall** | **88.29%** (TP=18,611, FN=2,469, total=21,080) | 40.63% (TP=8,564, FN=12,516) |
| **DEL recall** | **99.97%** (TP=11,934, FN=4, total=11,938) | 97.94% (TP=11,692, FN=246) |

## Dinara V36v — INS Size Estimation Fixes (June 2026)

### Changes

Three fixes to INS size estimation in `src/AssemblerSvAnchors.cpp`:

1. **early-CIGAR: max instead of mean** — INS CIGAR clusters now report the maximum observation size instead of the mean. Partial alignments under-report insertion size; the largest CIGAR observation is closest to truth. DEL clusters still use the mean (deletions are consistently sized across reads).

2. **indirect-covdrop: bounded by CIGAR/soft-clip evidence** — When CIGAR or soft-clip calls exist nearby, the `indirectBases/coverage` estimate is floored at the max observed size and capped at 3× to prevent gross over-estimation for small INS.

3. **soft-clip: actual sequence overlap** — Replaced the hardcoded `k-1=20` overlap subtraction with actual suffix-prefix sequence overlap detection between right-clip and left-clip contigs. Added `assembleClipContig()` (returns the contig sequence) and `suffixPrefixOverlap()`. Applied to all 3 soft-clip sizing paths (main, zero-anchor, bad-segment).

### V36v vs V36u vs sbx-assemble — HG002 Q100 v5.0q

Truvari params: `refdist=2000, chunksize=5000, pctseq=0, pick=multi, passonly, includebed`.

**pctsize=0 (position-only — "did we find it?"):**

| Metric | V36v | V36u | sbx-assemble |
|---|---|---|---|
| **ALL recall** | **99.66%** | 99.60% | 94.94% |
| INS recall | 99.47% | 99.37% | 94.45% |
| DEL recall | **100.00%** | 100.00% | 95.82% |
| ALL FN | 96 | 114 | 1,425 |

**pctsize=0.7 (size-aware — "did we get the right size?"):**

| Metric | V36v | V36u | Δ (v→u) | sbx-assemble |
|---|---|---|---|---|
| **ALL recall** | **70.68%** | 69.13% | **+1.6pp** | 86.88% |
| **INS recall** | **54.41%** | 51.98% | **+2.4pp** | 85.09% |
| DEL recall | 99.49% | 99.49% | 0 | 90.05% |
| ALL FN | 8,268 | 8,706 | −438 | 3,698 |
| INS FN | 8,208 | 8,650 | −442 | 2,688 |

### Analysis of remaining INS sizing gap

8,208 INS cases where V36v finds the right position but wrong size (ratio < 0.7). The gap to sbx-assemble (−30.7pp) is fundamental: 250bp reads can't span insertions >500bp, while assembly-based approaches produce full-length contigs.

**Breakdown by best-matching source:**

| Source | FN cases | Pattern |
|---|---|---|
| indirect-covdrop | ~47% | Small: over-estimates. Large: under-estimates |
| read-graph | ~18% | BFS paths can't span large insertions |
| soft-clip | ~14% | Max contig ~1kb (read length ceiling) |
| early-CIGAR | ~11% | Partial alignments, improved by max fix |
| large-ins-* | ~7% | Het calls under-estimate |

**By truth size:**

| Size range | FN cases | Issue |
|---|---|---|
| 50-99bp | ~25% | Over-estimation by indirect-covdrop |
| 100-199bp | ~21% | Mixed over/under |
| 200-499bp | ~18% | Transition zone |
| 500-999bp | ~13% | Under-estimation begins |
| 1000-4999bp | ~20% | All sources under-estimate |
| 5000+bp | ~4% | Unfixable with short reads |

### Cluster execution improvements

- **Shuffled case lists** prevent slow cases (150s each in complex repeat regions) from clustering in one chunk
- **Separate INS/DEL jobs** run fully in parallel instead of sequentially
- **20 INS + 10 DEL jobs** (16 CPUs, 14 parallel each) complete the full 36K-case benchmark in ~5 minutes
- Shuffled lists: `case_list_shuffled.txt`, `all_cases_shuffled.txt`
- Use `mktemp -u` for `--assemblyDirectory` (generates unique name without creating the directory — dinara creates it internally)

### Cluster Paths and Binaries

- **Binaries**: `dinara_v36{o,p,q,r,s,t,u,v}_ins` in `/sc1/groups/sbx/workspace/kyriakik/data/tools/`
- **samtools**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/samtools` (v1.21, compiled from source)

#### HG002 Case Extraction

Each case is a ±2000bp region around a truth position, containing:
- `reads.fa` — consensus reads overlapping the region
- `reference.fa` — reference sequence for the region
- `region.bam` + `region.bam.bai` — original BAM alignments for the region

| Dataset | Cases | Path |
|---|---|---|
| **HG002 Full INS** | 26,101 | `/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_ins_eval/cases/` |
| **HG002 Full DEL** | 15,472 | `/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval/cases/` |
| **HG002 INS pilot** | 100 | `/sc1/groups/sbx/workspace/kyriakik/structural_variants/ins_pilot/cases/` |

Case naming: `{chrom}_{pos}_{type}{size}` (e.g., `chr1_100025966_INS51`, `chr10_100092115_DEL68`).

The 28,942 INS truth entries map to 26,101 unique case names (1,400 positions have two het alleles with identical size → same case name). Similarly, 17,562 DEL entries map to 15,472 case names (2,090 duplicates). Both alleles share the same extracted region.

### HG008-T Somatic SV Benchmark

**Truth set**: HG008-T V0.5 Draft Somatic SV Benchmark, GRCh38.
- Source: [GIAB FTP](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/NIST_HG008-T_somatic-stvar-CNV_DraftBenchmark_V0.5-20260318/)
- PASS variants: 210 total (22 INS, 60 DEL, 54 DUP, 74 BND)
- All INS and DEL are ≥50bp (INS range: 60-4,775bp, DEL range: 54bp-6.3Mbp)
- BAM: `/sc1/groups/sbx/workspace/eizengaj/structural_variants/somatic_mapping/dedupped_sorted/HG008_T.intra-consensus.personalized.dedup.sorted.bam`

**Truth files on cluster**: `/sc1/groups/sbx/workspace/kyriakik/data/truth/HG008T_somatic_V0.5/`
- `GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz` — PASS variants (use for benchmarking)
- `GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf.gz` — all variants including filtered
- `GRCh38_HG008-T-V0.5_somatic-stvar-clonal.draftbenchmark.bed` — clonal benchmark regions
- `GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed` — clonal + subclonal regions
- `README.md` — benchmark description and truvari usage instructions

#### HG008-T Case Extraction

| Dataset | Cases | Path |
|---|---|---|
| **HG008-T INS** | 22 | `/sc1/groups/sbx/workspace/kyriakik/structural_variants/hg008t_ins_eval/cases/` |
| **HG008-T DEL** | 60 | `/sc1/groups/sbx/workspace/kyriakik/structural_variants/hg008t_del_eval/cases/` |

Case list files: `ins_cases.txt` in `.../hg008t_ins_eval/`, `del_cases.txt` in `.../hg008t_del_eval/`.
Extraction script: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/extract_hg008t_cases.py`.

#### Results

Results are per-version log files in `results_v36{version}/` under each eval directory.

| Version | INS results | DEL results |
|---|---|---|
| V36s | `.../full_ins_eval/results_v36s/` | `.../full_del_eval/results_v36s/` |
| V36t | `.../full_ins_eval/results_v36t/` | `.../full_del_eval/results_v36t/` |
| V36u | `.../full_ins_eval/results_v36u/` | `.../full_del_eval/results_v36u/` |
| V36v | `.../full_ins_eval/results_v36v/` | `.../full_del_eval/results_v36v/` |

#### Supporting Files

- **Case lists** (for Slurm array jobs): `ins_cases.txt` in `.../full_ins_eval/`, `del_cases.txt` in `.../full_del_eval/`
- **Truth positions**: `all_ins_positions.tsv` (28,942 entries) in `.../full_ins_eval/`
- **VCF generation scripts**: `ins_logs_to_vcf.py` in `.../full_ins_eval/`, `logs_to_vcf.py` in `.../full_del_eval/`
- **Analysis scripts**: `analyze_fast.py`, `analyze_fn.py`, `analyze_fn_deep.py` in `.../full_ins_eval/`
- **Truth VCFs**: `/sc1/groups/sbx/workspace/kyriakik/data/truth/GRCh38_HG2-T2TQ100-V1.1_stvar.filt.vcf.gz`
- **Filtered truth**: `stvar_INS50.vcf.gz`, `stvar_DEL50.vcf.gz` in `.../structural_variants/`

## Dinara V36z — Geometric Mean INS Sizing (June 2026)

### Problem

indirect-covdrop INS sizing has a systematic size-dependent bias:
- Small INS (50-100bp): over-estimates by 3-4× (read length counted, not insertion length)
- Medium INS (200-300bp): approximately correct
- Large INS (1000+bp): under-estimates by 4-10× (bounded by read length)

CIGAR/soft-clip sizes have the opposite bias: under-estimate for tandem repeats (see one repeat unit, truth is N units).

### Approach

When both signals exist at the same locus, `sqrt(CIGAR_size × indirect_size)` lands near truth. Added `geomeanRefineInsCalls` lambda: for each `indirect-covdrop` call where the best CIGAR INS (or soft-clip avgClipLen fallback) diverges by >40%, emit an additional `+geomean` call. Purely additive — original calls preserved.

The soft-clip fallback (when CIGAR INS < 30bp) is why this works far better than predicted: soft-clip evidence exists at thousands of loci where CIGAR INS operations don't.

### V36z Benchmark Results — HG002 Q100 v5.0q

| Metric | V36w (baseline) | V36y | V36z |
|--------|---:|---:|---:|
| INS recall ps=0 | 99.45% (TP=17882, FN=98) | 99.46% (TP=17916, FN=97) | 99.46% (TP=17915, FN=98) |
| INS recall ps=0.7 | 54.49% (TP=9797, FN=8183) | 54.89% (TP=9887, FN=8126) | **60.81% (TP=10953, FN=7060)** |
| DEL recall ps=0 | 100.00% (TP=10143, FN=0) | 100.00% (TP=10175, FN=0) | 100.00% (TP=10175, FN=0) |
| DEL recall ps=0.7 | 99.54% (TP=10096, FN=47) | 99.55% (TP=10129, FN=46) | 99.55% (TP=10129, FN=46) |

**V36z vs V36y: +1066 INS TP at ps=0.7 (+5.92% recall)**. No DEL regression.

Breakdown: 1196 gained, 150 lost (truvari matching artifacts). 806 gained directly from `+geomean` calls, 390 from improved matching landscape. Geomean gains by truth size: 187 (50-100bp), 394 (100-200bp), 201 (200-500bp), 24 (500bp+). Average pss of geomean gains: 0.86.

### Remaining FN Analysis (V36y baseline)

**DEL ps=0.7 (46 FN):** Dead end. 38/46 are 50-100bp tandem repeat DELs. 0/46 have CIGAR D ops within 30% of truth. No actionable signal.

**INS ps=0 (97 FN):** 72 emit only DEL calls (no INS call). 64 of those have spanning reads but zero diagonal shift — k-mers from inserted tandem repeat match reference. 46/97 are in multi-truth clusters. ~2 fixable (position near BED boundary).

**INS ps=0.7 (sizing errors):** Dominated by indirect-covdrop (2429 cases). Systematic bias: over for small INS, under for large INS, bounded by read length (~250bp). Simple correction factors (×0.5, ×0.7, estSize²/250) all make things worse — they help one direction but hurt the other. The geometric mean approach is the only method that works because it combines two signals with opposite biases.

### Cluster Paths

- **Binary**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36z`
- **INS results**: `.../full_ins_eval/results_v36z/` (26,101 logs)
- **DEL results**: `.../full_del_eval/results_v36z/` (10,173 logs)
- **INS VCF**: `.../full_ins_eval/dinara_v36z_ins_sorted.vcf.gz` (247,609 calls, 8,014 +geomean)
- **DEL VCF**: `.../full_del_eval/dinara_v36z_dels_sorted.vcf.gz`

### V36z DEL Recall Ceiling

DEL detection is at 100.00% recall at ps=0 (10,175 / 10,175, 0 FN) and 99.55% at ps=0.7 (10,129 / 10,175). The 46 FN at ps=0.7 are 50-100bp tandem repeat DELs with no actionable CIGAR signal — this is the practical recall ceiling for DEL.

---

## Dinara V37a — Source Elimination (June 2026)

### Motivation

Greedy elimination analysis on V36z identified sources with few or zero sole contributions to recall. Removing low-value sources reduces false positives and call volume without meaningful recall loss.

### Eliminated Sources

**INS (11 of 22 sources removed):** `large-ins-het+CIGAR`, `reversed-BP+CIGAR`, `large-ins-single-het+CIGAR`, `CIGAR-covdrop+CIGAR`, `large-ins+CIGAR`, `softclip-unpaired`, `large-ins-single+CIGAR`, `hit-depth`, `soft-clip+CIGAR`, `reversed-BP`, `large-ins-het`.

**DEL (9 of 18 sources removed):** `depth-scan-sub`, `depth-scan-hom`, `depth-scan-het`, `path-based`, `split-read`, `cigar-covdrop`, `INV-cluster`, `depth-deficit-het`, `per-read-DEL`.

Implementation: source-filtering lambdas in `AssemblerSvAnchors.cpp` erase disabled sources before output. Depth-scan call block in `main.cpp` disabled entirely (all 3 depth-scan sources removed).

### V37a Benchmark Results — HG002 Q100 v5.0q

| Metric | V36z (baseline) | V37a |
|--------|---:|---:|
| INS recall ps=0 | 99.46% (TP=17915, FN=98) | 99.27% (TP=17882, FN=131) |
| INS recall ps=0.7 | 60.81% (TP=10953, FN=7060) | 59.62% (TP=10739, FN=7274) |
| DEL recall ps=0 | 100.00% (TP=10175, FN=0) | 100.00% (TP=10175, FN=0) |
| DEL recall ps=0.7 | 99.55% (TP=10129, FN=46) | 99.41% (TP=10115, FN=60) |

INS ps=0.7: −214 TP (−1.19%), worse than predicted −0.91%. Source interactions during dedup cause non-additive losses. INS ps=0: −33 TP (−0.18%), not predicted by the ps=0.7 analysis. DEL ps=0: no regression. DEL ps=0.7: −14 TP (−0.14%).

FP reduction: INS ps=0.7 −1,166 FP, DEL ps=0.7 −24,740 FP. Total call volume: INS 247K → 231K, DEL 173K → 143K.

### Cluster Paths

- **Binary**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v37a`
- **INS results**: `.../full_ins_eval/results_v37a/` (26,101 logs)
- **DEL results**: `.../full_del_eval/results_v37a/` (10,173 logs)
- **INS VCF**: `.../full_ins_eval/dinara_v37a_ins_sorted.vcf.gz` (230,911 calls)
- **DEL VCF**: `.../full_del_eval/dinara_v37a_dels_sorted.vcf.gz` (144,700 calls)

---

## Dinara V37b — Bidirectional Geomean + Emission Bug Fix (June 2026)

### Changes

1. **Bidirectional geomean**: `geomeanRefineInsCalls` now handles both directions — indirect-covdrop over-estimates (IC > CIGAR) and under-estimates (IC < CIGAR/soft-clip). Previously only the over-estimate direction was covered. The ratio and guard checks now use `min/max` to be symmetric.

2. **Emission bug fix**: `cigar-covdrop` DEL calls had an inline `cout` before `filterDisabledDelSources()` ran, so they appeared in VCF output despite being disabled. Wrapped in `disabledDelSources.count()` check.

### V37b Benchmark Results — HG002 Q100 v5.0q

| Metric | V36z (baseline) | V37a | V37b |
|--------|---:|---:|---:|
| INS recall ps=0 | 99.46% (TP=17915, FN=98) | 99.27% (TP=17882, FN=131) | 99.28% (TP=17884, FN=129) |
| INS recall ps=0.7 | 60.81% (TP=10953, FN=7060) | 59.62% (TP=10739, FN=7274) | 59.66% (TP=10746, FN=7267) |
| DEL recall ps=0 | 100.00% (TP=10175, FN=0) | 100.00% (TP=10175, FN=0) | 100.00% (TP=10175, FN=0) |
| DEL recall ps=0.7 | 99.55% (TP=10129, FN=46) | 99.41% (TP=10115, FN=60) | 99.41% (TP=10115, FN=60) |

V37b vs V37a: +7 INS TP at ps=0.7, +2 at ps=0. The bidirectional geomean gained less than the predicted 61 because the analysis used soft-clip VCF call sizes, but the geomean lambda uses `softClipBPs.avgClipLen` from BAM parsing — these differ.

### Source Analysis Summary

INS FN at ps=0.7 (7,060 total): 6,962 are sizing misses (found at ps=0), only 98 are position misses. The dominant problem is size estimation, not detection.

- 1,824 sizing misses have only `indirect-covdrop` — no CIGAR, no soft-clip, no second signal. No current path to fix these.
- `indirect-covdrop` sizing bias: 3.5× over-estimate for 50-100bp truth, 0.19× under-estimate for 1000-5000bp truth.
- Cross-source geomean could theoretically fix ~453 sizing misses, but most overlap with existing geomean and the net gain is small.

DEL is at the practical recall ceiling (100% ps=0, 99.41% ps=0.7). The 60 FN are dominated by 50-100bp tandem repeat DELs.

### Cluster Paths

- **Binary**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v37b`
- **INS results**: `.../full_ins_eval/results_v37b/` (26,101 logs)
- **DEL results**: `.../full_del_eval/results_v37b/` (10,173 logs)
- **INS VCF**: `.../full_ins_eval/dinara_v37b_ins_sorted.vcf.gz` (230,970 calls)
- **DEL VCF**: `.../full_del_eval/dinara_v37b_dels_sorted.vcf.gz` (140,425 calls)

---

## Dinara V36y — CIGAR-Guided INS Refinement (June 2026)

### Changes

Three changes over V36w baseline:

1. **`cigarRefineInsCalls` lambda (additive)**: For each INS call whose size significantly exceeds the best CIGAR INS size (ratio < 0.7), emits an additional `+CIGAR` call alongside the original. Truvari `--pick multi` selects whichever matches truth better. Called at all 3 INS emission sites.

2. **Assembly overlap improvement**: `assembleClipSequences` refactored to `assembleClipContig` (returns full contig string). New `suffixPrefixOverlap` function replaces hardcoded `k-1=20` overlap estimate with actual suffix-prefix overlap detection between assembled soft-clip contigs.

3. **Removed destructive CIGAR refinement** (was in V36x, caused regression): Inline CIGAR refinement in both indirect-covdrop paths replaced `estSize` with `maxObservedSize`. For tandem repeats, CIGAR sees one repeat unit while truth is N units — this destroyed accurate indirect estimates. Removed entirely in V36y.

### V36y Benchmark Results — HG002 Q100 v5.0q

| Metric | V36w (baseline) | V36x (regression) | V36y |
|--------|---:|---:|---:|
| INS recall ps=0 | 99.45% (TP=17882, FN=98) | 99.45% (TP=17881, FN=99) | **99.46% (TP=17916, FN=97)** |
| INS recall ps=0.7 | 54.49% (TP=9797, FN=8183) | 52.44% (TP=9429, FN=8551) | **54.89% (TP=9887, FN=8126)** |
| DEL recall ps=0 | 100.00% (TP=10143, FN=0) | — | **100.00% (TP=10175, FN=0)** |
| DEL recall ps=0.7 | 99.54% (TP=10096, FN=47) | — | **99.55% (TP=10129, FN=46)** |

Net gain over V36w: **+90 INS TP at ps=0.7** (+0.40%), no regressions.

### V36x Regression Analysis

V36x applied destructive CIGAR refinement in both indirect-covdrop emission paths: when `maxObservedSize ∈ [50, 500)` and `estSize > maxObservedSize`, it replaced `estSize` with `maxObservedSize`. This caused:

- **1260 lost** truth variants (pss dropped below 0.7)
- **903 gained** truth variants (pss rose above 0.7)
- **Net: -357 TP**

Root cause: for tandem repeat insertions, CIGAR sees one repeat unit (e.g., 99bp) while truth is N units (e.g., 297bp). The indirect-covdrop estimate (`indirectBases/coverage`) was already close to truth but got replaced with the single-unit CIGAR size. The `estSize/maxObservedSize` ratio distributions for lost and gained cases overlap completely (both peak at 2-4×), so no threshold can separate them.

**Lesson**: Never replace an existing estimate destructively. Always emit both original and refined as separate calls, letting truvari `--pick multi` choose the best match. See `AGENTS.md` for coding guidelines.

### Cluster Paths

- **Binary**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36y`
- **INS results**: `.../full_ins_eval/results_v36y/` (26,101 logs)
- **DEL results**: `.../full_del_eval/results_v36y/` (10,173 logs)
- **INS VCF**: `.../full_ins_eval/dinara_v36y_ins_sorted.vcf.gz`
- **DEL VCF**: `.../full_del_eval/dinara_v36y_dels_sorted.vcf.gz`

### Standardized Truvari Benchmark Parameters

#### HG002 Germline — Q100 v5.0q (GIAB-recommended)

Source: `NIST_HG002_v5.0q_variant-benchmarksets_README.md` on GIAB FTP.

**GIAB-recommended command (with refine):**
```bash
truvari bench \
  -b GRCh38_HG2-T2TQ100-V1.1_stvar.filt.vcf.gz \
  -c {calls}.vcf.gz \
  -f {reference.fasta} \
  --includebed GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed \
  --pick ac --passonly \
  -r 2000 -C 5000 \
  --refine \
  -o {output}
```

**Our evaluation command (without refine, for FN analysis):**
```bash
truvari bench \
  -b /tmp/stvar_{INS,DEL}50.vcf.gz \
  -c {calls}.vcf.gz \
  --includebed /tmp/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed \
  --passonly --pick multi --pctseq 0 --pctsize {0,0.7} \
  -r 2000 -C 5000 \
  -o {output}
```

| Parameter | GIAB recommended | Our runs | Notes |
|-----------|-----------------|----------|-------|
| `--pick` | `ac` | `multi` | `ac` matches by allele count (diploid-aware); `multi` allows multiple matches per call. Both are valid for recall. |
| `--refdist` | 2000 | 2000 | ✅ Matches |
| `--chunksize` | 5000 | 5000 | ✅ Matches |
| `--pctseq` | 0.7 (default) | 0 | We disable sequence similarity — our calls lack ALT sequences |
| `--pctsize` | 0.7 (default) | 0 and 0.7 | We run both: pctsize=0 for position-only recall, pctsize=0.7 for size accuracy |
| `--refine` | yes | no | Refine resolves complex representations but makes FN analysis harder |
| `--includebed` | yes | yes | ✅ Matches |
| `--passonly` | yes | yes | ✅ Matches |
| `ALT=*` filter | yes (pre-filter VCF) | yes | `stvar.filt.vcf.gz` has `ALT="."` removed via `bcftools view -e 'ALT="."'` |

**Pre-filtering note:** The v5.0q README warns that `ALT=*` variants cause incorrect truvari categorization. The truth VCF on the cluster (`GRCh38_HG2-T2TQ100-V1.1_stvar.filt.vcf.gz`) already has these filtered out. The per-type VCFs (`stvar_INS50.vcf.gz`, `stvar_DEL50.vcf.gz`) are extracted from the filtered file.

#### HG008-T Somatic — V0.5 (GIAB-recommended)

Source: `README.md` in HG008-T V0.5 truth set directory.

**GIAB-recommended command (two pctsize configurations):**
```bash
# Position-only (pctsize=0)
truvari bench \
  -b GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz \
  -c {calls}.vcf.gz \
  --reference GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
  --includebed GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed \
  --refdist 1000 --pctseq 0 --pctsize 0 --pctovl 0 \
  --sizemax -1 --passonly --pick multi --typeignore \
  -o {output_ps0}

# Size-aware (pctsize=0.7)
truvari bench \
  -b GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz \
  -c {calls}.vcf.gz \
  --reference GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
  --includebed GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed \
  --refdist 1000 --pctseq 0 --pctsize 0.7 --pctovl 0 \
  --sizemax -1 --passonly --pick multi --typeignore \
  -o {output_ps07}
```

| Parameter | Value | Notes |
|-----------|-------|-------|
| `--refdist` | 1000 | Somatic SVs have less precise breakpoints |
| `--pctseq` | 0 | No sequence similarity |
| `--pctsize` | 0 and 0.7 | Two runs recommended |
| `--pctovl` | 0 | No reciprocal overlap required |
| `--sizemax` | -1 | No upper size limit (somatic SVs can be very large) |
| `--passonly` | yes | Only PASS variants in truth |
| `--pick` | multi | Multiple matches allowed |
| `--typeignore` | yes | Don't enforce type matching (DUP/INS/DEL interchangeable) |
| `--reference` | GIAB GRCh38 | Required for BND comparison and symbolic variant resolution |
| `--includebed` | clonal+subclonal BED | Restricts evaluation to benchmark regions |

**Key differences from HG002 germline:**
- `--typeignore` — somatic truth set has DUP/BND types that callers may represent differently
- `--sizemax -1` — no upper limit (germline uses 50000)
- `--refdist 1000` — tighter than germline's 2000 (somatic truth set is smaller, positions are well-characterized)
- `--reference` required — needed for BND-to-BND comparison in truvari v5+
- No `--refine` — not recommended for somatic benchmark

---

### How to Run a Full Benchmark

The benchmark runs ~36,000 independent cases across 4 datasets. Each case takes <1s (median ~0.7s). With 500 concurrent Slurm array tasks, the full benchmark completes in ~15 minutes wall time.

#### Cluster resources

- **Partition**: `batch_cpu` — 137 nodes, 64 CPUs / 348 GB RAM each
- **QoS**: `3h` (max wall time 3 hours, more than enough)
- **Concurrency**: no cap — the Slurm scheduler manages task scheduling across available nodes

#### Case counts

| Dataset | Cases | Case list file | Cases dir |
|---------|-------|---------------|-----------|
| HG002 INS | 26,101 | `full_ins_eval/case_list.txt` | `full_ins_eval/cases/` |
| HG002 DEL | 10,173 | `full_del_eval/all_cases.txt` | `full_del_eval/cases/` |
| HG008-T INS | 22 | `hg008t_ins_eval/ins_cases.txt` | `hg008t_ins_eval/cases/` |
| HG008-T DEL | 60 | `hg008t_del_eval/del_cases.txt` | `hg008t_del_eval/cases/` |

All paths relative to `/sc1/groups/sbx/workspace/kyriakik/structural_variants/`.

#### Step 1: Build and deploy binary

```bash
# Build locally (in dev container)
cd /workspaces/dinara/build_release && make -j$(nproc) -C Executable

# Copy to cluster
scp build_release/Executable/dinara ec-hub:/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36X
```

#### Step 2: Create parallel Slurm scripts

Each script processes a chunk of the shuffled case list using `bash` background jobs (14 parallel per 16-CPU node). Use **shuffled** case lists to distribute slow cases (complex repeat regions, ~150s each) evenly across chunks. Submit INS and DEL as **separate** jobs so they run on different nodes in parallel.

```bash
#!/bin/bash
#SBATCH --job-name=v36X-ins
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --qos=3h
#SBATCH --partition=batch_cpu

set -uo pipefail
DINARA=/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36X
BASE=/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_ins_eval
RESULTS=$BASE/results_v36X
CASE_LIST=$BASE/case_list_shuffled.txt
NCHUNKS=20
NCPU=14
CHUNK=$SLURM_ARRAY_TASK_ID

total=$(wc -l < "$CASE_LIST")
per_chunk=$(( (total + NCHUNKS - 1) / NCHUNKS ))
start=$(( CHUNK * per_chunk + 1 ))
end=$(( start + per_chunk - 1 ))
[ $end -gt $total ] && end=$total
[ $start -gt $total ] && exit 0

mkdir -p "$RESULTS"
running=0
while IFS= read -r name; do
    outfile="$RESULTS/${name}.log"
    casedir="$BASE/cases/$name"
    [ -f "$outfile" ] && [ -s "$outfile" ] && grep -q "buildSvMSA completed" "$outfile" 2>/dev/null && continue
    [ ! -f "$casedir/reads.fa" ] && continue
    (
        tmpdir=$(mktemp -u /tmp/dinara_v36X_XXXXXX)
        "$DINARA" --command svanchors \
            --reference "$casedir/reference.fa" \
            --input "$casedir/reads.fa" \
            --bam "$casedir/region.bam" \
            --assemblyDirectory "$tmpdir" \
            --Kmers.k 10 --Kmers.minimizerW 6 \
            > "$outfile" 2>&1 || true
        rm -rf "$tmpdir"
    ) &
    running=$((running + 1))
    if [ $running -ge $NCPU ]; then
        wait -n 2>/dev/null || true
        running=$((running - 1))
    fi
done < <(sed -n "${start},${end}p" "$CASE_LIST")
wait
```

Create variants for each dataset:
- **HG002 INS**: `NCHUNKS=20`, `CASE_LIST=case_list_shuffled.txt`, `--array=0-19`
- **HG002 DEL**: `NCHUNKS=10`, `CASE_LIST=all_cases_shuffled.txt`, `--array=0-9`
- **HG008-T INS/DEL**: single job each (22/60 cases, no chunking needed)

Key design choices:
- **`mktemp -u`** for `--assemblyDirectory`: generates a unique name **without** creating the directory. Dinara creates it internally and fails if it already exists.
- **Shuffled case lists**: prevents slow cases from clustering in one chunk. Create with `shuf case_list.txt > case_list_shuffled.txt`.
- **Skip-if-done**: checks for `buildSvMSA completed` in existing log — safe to resubmit after partial failures.
- **Separate INS/DEL jobs**: run on different nodes in parallel instead of sequentially.

#### Step 3: Submit all jobs

```bash
ssh ec-hub 'cd /sc1/groups/sbx/workspace/kyriakik/structural_variants && \
  mkdir -p full_ins_eval/results_v36X full_del_eval/results_v36X \
           hg008t_ins_eval/results_v36X hg008t_del_eval/results_v36X && \
  sbatch --array=0-19 full_ins_eval/run_v36X_ins.sh && \
  sbatch --array=0-9  full_del_eval/run_v36X_del.sh'
```

20 INS + 10 DEL jobs × 14 parallel cases each = ~420 concurrent cases. Full 36K-case benchmark completes in ~5 minutes.

#### Step 4: Monitor progress

```bash
# Job status
ssh ec-hub 'squeue -u kyriakik -o "%.10i %.12j %.8T %.10M %.4C" | head -20'

# Count completed files
ssh ec-hub 'echo "INS: $(ls .../full_ins_eval/results_v36X/ | wc -l) / 26101"'
ssh ec-hub 'echo "DEL: $(ls .../full_del_eval/results_v36X/ | wc -l) / 10173"'

# Resubmit (safe — skip-if-done logic)
ssh ec-hub 'sbatch --array=0-19 .../full_ins_eval/run_v36X_ins.sh'
```

#### Step 5: Generate VCFs, sort, compress, index

Run on ec-hub (or submit as a Slurm job for large datasets):

```bash
ssh ec-hub 'cd /sc1/groups/sbx/workspace/kyriakik/structural_variants && \
  module load BCFtools/1.21-GCC-13.3.0 && \
  module load HTSlib/1.21-GCC-13.3.0 && \
  cd full_ins_eval && \
  python3 ins_logs_to_vcf.py all_ins_positions.tsv results_v36X dinara_v36X_ins.vcf && \
  bcftools sort dinara_v36X_ins.vcf -o dinara_v36X_ins_sorted.vcf && \
  bgzip -c dinara_v36X_ins_sorted.vcf > dinara_v36X_ins_sorted.vcf.gz && \
  tabix -p vcf dinara_v36X_ins_sorted.vcf.gz && \
  cd ../full_del_eval && \
  python3 logs_to_vcf.py cases results_v36X dinara_v36X_dels.vcf && \
  bcftools sort dinara_v36X_dels.vcf -o dinara_v36X_dels_sorted.vcf && \
  bgzip -c dinara_v36X_dels_sorted.vcf > dinara_v36X_dels_sorted.vcf.gz && \
  tabix -p vcf dinara_v36X_dels_sorted.vcf.gz'
```

#### Step 6: Copy VCFs locally and run truvari

Truvari runs locally in the dev container (`/home/vscode/.local/bin/truvari`).

```bash
# Copy result VCFs from cluster
scp ec-hub:/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_ins_eval/dinara_v36X_ins_sorted.vcf.gz /tmp/
scp ec-hub:/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_ins_eval/dinara_v36X_ins_sorted.vcf.gz.tbi /tmp/
scp ec-hub:/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval/dinara_v36X_dels_sorted.vcf.gz /tmp/
scp ec-hub:/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval/dinara_v36X_dels_sorted.vcf.gz.tbi /tmp/

# Copy truth/benchmark files (if not already local)
scp ec-hub:/sc1/groups/sbx/workspace/kyriakik/data/truth/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed /tmp/

# HG002 germline INS benchmark (see Standardized Truvari Parameters for rationale)
truvari bench -b /tmp/stvar_INS50.vcf.gz -c /tmp/dinara_v36X_ins_sorted.vcf.gz \
  --includebed /tmp/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed \
  --passonly --pick multi --pctseq 0 --pctsize 0 \
  -r 2000 -C 5000 \
  -o /tmp/truvari_v36X/ins_ps0

truvari bench -b /tmp/stvar_INS50.vcf.gz -c /tmp/dinara_v36X_ins_sorted.vcf.gz \
  --includebed /tmp/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed \
  --passonly --pick multi --pctseq 0 --pctsize 0.7 \
  -r 2000 -C 5000 \
  -o /tmp/truvari_v36X/ins_ps07

# HG002 germline DEL benchmark
truvari bench -b /tmp/stvar_DEL50.vcf.gz -c /tmp/dinara_v36X_dels_sorted.vcf.gz \
  --includebed /tmp/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed \
  --passonly --pick multi --pctseq 0 --pctsize 0 \
  -r 2000 -C 5000 \
  -o /tmp/truvari_v36X/del_ps0

truvari bench -b /tmp/stvar_DEL50.vcf.gz -c /tmp/dinara_v36X_dels_sorted.vcf.gz \
  --includebed /tmp/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed \
  --passonly --pick multi --pctseq 0 --pctsize 0.7 \
  -r 2000 -C 5000 \
  -o /tmp/truvari_v36X/del_ps07

# HG008-T somatic INS benchmark
truvari bench -b /tmp/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz \
  -c /tmp/dinara_v36X_hg008t_ins.vcf.gz \
  --reference /tmp/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
  --includebed /tmp/GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed \
  --refdist 1000 --pctseq 0 --pctsize 0 --pctovl 0 \
  --sizemax -1 --passonly --pick multi --typeignore \
  -o /tmp/truvari_v36X/hg008t_ins_ps0

truvari bench -b /tmp/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz \
  -c /tmp/dinara_v36X_hg008t_ins.vcf.gz \
  --reference /tmp/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
  --includebed /tmp/GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed \
  --refdist 1000 --pctseq 0 --pctsize 0.7 --pctovl 0 \
  --sizemax -1 --passonly --pick multi --typeignore \
  -o /tmp/truvari_v36X/hg008t_ins_ps07

# HG008-T somatic DEL benchmark
truvari bench -b /tmp/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz \
  -c /tmp/dinara_v36X_hg008t_dels.vcf.gz \
  --reference /tmp/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
  --includebed /tmp/GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed \
  --refdist 1000 --pctseq 0 --pctsize 0 --pctovl 0 \
  --sizemax -1 --passonly --pick multi --typeignore \
  -o /tmp/truvari_v36X/hg008t_del_ps0

truvari bench -b /tmp/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz \
  -c /tmp/dinara_v36X_hg008t_dels.vcf.gz \
  --reference /tmp/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
  --includebed /tmp/GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed \
  --refdist 1000 --pctseq 0 --pctsize 0.7 --pctovl 0 \
  --sizemax -1 --passonly --pick multi --typeignore \
  -o /tmp/truvari_v36X/hg008t_del_ps07

# Read all results
for d in ins_ps0 ins_ps07 del_ps0 del_ps07 hg008t_ins_ps0 hg008t_ins_ps07 hg008t_del_ps0 hg008t_del_ps07; do
  [ -f "/tmp/truvari_v36X/$d/summary.json" ] && \
  python3 -c "import json; d=json.load(open('/tmp/truvari_v36X/$d/summary.json')); print(f'$d: TP={d[\"TP-base\"]}, FN={d[\"FN\"]}, recall={d[\"recall\"]:.4f}')"
done
```

#### Timeline summary

| Phase | Wall time | Notes |
|-------|-----------|-------|
| Build + deploy | ~2 min | `make -j$(nproc)` + scp |
| Create scripts | ~2 min | One-time per version |
| Slurm execution | ~5 min | 30 jobs × 14 parallel, shuffled case lists |
| VCF generation | ~5 min | Python log parsing + bcftools sort/bgzip |
| Copy + truvari | ~5 min | 8 truvari runs, each <30s |
| **Total** | **~20 min** | End-to-end from code change to results |

### Running sbx-assemble (Competitor Baseline)

sbx-assemble is the competing assembly-based SV caller. It runs as a single whole-genome job (not per-case like dinara) with a 5-step pipeline: assemble → align contigs → call SVs (SVIM-asm) → truvari bench → truvari refine.

#### Prerequisites

- **Conda env**: `sv_bench` on the cluster (has truvari, svim-asm, pysam==0.22.1)
- **Binaries**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/assemble` and `contig_aligner`
- **Pipeline script**: `/sc1/groups/sbx/workspace/kyriakik/sbx-assemble/benchmark_full.sh`

#### HG002 Germline Benchmark

Single Slurm job, ~1-2 hours, 32 CPUs:

```bash
# Run on SBX-D (primary evaluation BAM)
ssh ec-hub 'cd /sc1/groups/sbx/workspace/kyriakik/sbx-assemble && \
  sbatch benchmark_full.sh SBX-D'

# Or run both BAMs in parallel
ssh ec-hub 'cd /sc1/groups/sbx/workspace/kyriakik/sbx-assemble && \
  sbatch benchmark_full.sh SBX-S && \
  sbatch benchmark_full.sh SBX-D'
```

Output goes to `/sc1/groups/sbx/workspace/kyriakik/structural_variants/test_runs/bench_<timestamp>/SBX-D/`.

The script runs truvari with: `--pick multi -p 0 -P 0 --passonly -r 2000 -C 5000 --includebed`. Then runs `truvari refine --use-original-vcfs --align=mafft`.

#### HG008-T Somatic Benchmark

sbx-assemble doesn't have a somatic benchmark script. To run it on HG008-T, modify the pipeline manually:

```bash
# On ec-hub: create a somatic benchmark script
cat > /sc1/groups/sbx/workspace/kyriakik/sbx-assemble/benchmark_hg008t.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name="sbx-hg008t"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --time=03:00:00
#SBATCH --qos=3h

set -euo pipefail

DATA=/sc1/groups/sbx/workspace/kyriakik/data
OUTDIR=/sc1/groups/sbx/workspace/kyriakik/structural_variants/test_runs/sbx_hg008t_$(date +%Y%m%d_%H%M%S)
mkdir -p "$OUTDIR"

EXE=$DATA/tools/assemble
ALIGNER=$DATA/tools/contig_aligner
REF=$DATA/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
BAM=/sc1/groups/sbx/workspace/eizengaj/structural_variants/somatic_mapping/dedupped_sorted/HG008_T.intra-consensus.personalized.dedup.sorted.bam
TRUTH=$DATA/truth/HG008T_somatic_V0.5/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz
TRUTH_BED=$DATA/truth/HG008T_somatic_V0.5/GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed
GIAB_REF=$DATA/reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
THREADS=32
NAME="final_assemblies"

module load GCC/13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load Micromamba/2.0.7-0
eval "$(micromamba shell hook --shell bash)"
micromamba activate sv_bench

exec > >(tee "$OUTDIR/benchmark.log") 2>&1

# Step 1: Assembly
$EXE -x "$REF" -i "$BAM" -t $THREADS -o "$OUTDIR"

# Step 2: Align contigs
"$ALIGNER" -t $THREADS "$REF" "$OUTDIR/$NAME.fa" \
    | samtools sort -@ 4 -O BAM > "$OUTDIR/$NAME.GRCh38.bam"
samtools index -@ $THREADS "$OUTDIR/$NAME.GRCh38.bam"

# Step 3: Call SVs
svim-asm haploid --min_mapq 0 --min_sv_size 30 \
    --query_gap_tolerance 50000 --reference_gap_tolerance 50000 \
    --query_overlap_tolerance 50000 --reference_overlap_tolerance 50000 \
    "$OUTDIR" "$OUTDIR/$NAME.GRCh38.bam" "$REF" > /dev/null

bcftools view -e 'ALT="." || TYPE="bnd"' \
    -Oz -o "$OUTDIR/variants.filtered.vcf.gz" "$OUTDIR/variants.vcf"
tabix -p vcf "$OUTDIR/variants.filtered.vcf.gz"

# Step 4: Truvari bench (pctsize=0)
truvari bench \
    --base="$TRUTH" --comp="$OUTDIR/variants.filtered.vcf.gz" \
    --reference "$GIAB_REF" \
    --includebed "$TRUTH_BED" \
    --refdist 1000 --pctseq 0 --pctsize 0 --pctovl 0 \
    --sizemax -1 --passonly --pick multi --typeignore \
    --output="$OUTDIR/truvari_ps0"

# Step 5: Truvari bench (pctsize=0.7)
truvari bench \
    --base="$TRUTH" --comp="$OUTDIR/variants.filtered.vcf.gz" \
    --reference "$GIAB_REF" \
    --includebed "$TRUTH_BED" \
    --refdist 1000 --pctseq 0 --pctsize 0.7 --pctovl 0 \
    --sizemax -1 --passonly --pick multi --typeignore \
    --output="$OUTDIR/truvari_ps07"

# Results
for d in truvari_ps0 truvari_ps07; do
  echo "=== $d ==="
  python3 -c "import json; s=json.load(open('$OUTDIR/$d/summary.json')); print(f'TP={s[\"TP-base\"]} FP={s[\"FP\"]} FN={s[\"FN\"]} P={s[\"precision\"]:.4f} R={s[\"recall\"]:.4f}')"
done
EOF
```

Submit with:
```bash
ssh ec-hub 'sbatch /sc1/groups/sbx/workspace/kyriakik/sbx-assemble/benchmark_hg008t.sh'
```

#### Reading sbx-assemble results

```bash
# HG002 germline (latest run)
ssh ec-hub 'cat /sc1/groups/sbx/workspace/kyriakik/structural_variants/test_runs/bench_20260519_155309/SBX-D/truvari/summary.json' | python3 -c "
import json, sys; s=json.load(sys.stdin)
print(f'TP={s[\"TP-base\"]} FP={s[\"FP\"]} FN={s[\"FN\"]} P={s[\"precision\"]:.4f} R={s[\"recall\"]:.4f} F1={s[\"f1\"]:.4f}')
"

# HG008-T somatic (after running benchmark_hg008t.sh)
ssh ec-hub 'cat /sc1/.../sbx_hg008t_*/truvari_ps0/summary.json' | python3 -c "
import json, sys; s=json.load(sys.stdin)
print(f'TP={s[\"TP-base\"]} FP={s[\"FP\"]} FN={s[\"FN\"]} P={s[\"precision\"]:.4f} R={s[\"recall\"]:.4f}')
"
```

#### Key differences: dinara vs sbx-assemble

| Aspect | dinara | sbx-assemble |
|--------|--------|-------------|
| Architecture | Per-case (36K independent jobs) | Whole-genome (single job) |
| Wall time | ~5-15 min (massively parallel) | ~1-2 hours (single node) |
| SV calling | Direct from reads (k-mer chaining + multi-source detection) | Assembly → contig alignment → SVIM-asm |
| Dependencies | None (static binary) | conda env (truvari, svim-asm, pysam, mafft) |
| Truvari refine | Not used | Used (resolves complex representations) |
| GBZ database | Not needed | Required per sample type |

---

### Connecting to the Cluster

#### SSH Configuration

The SSH config is at `~/.ssh/config` in the dev container. Key: `~/.ssh/id_ed25519` (label: `ona-dinara-agent`).

```
Host ec-hub
  Hostname ec-hub.sc1.science.roche.com
  User kyriakik
  IdentityFile ~/.ssh/id_ed25519

Host sc1
  Hostname lb022dev.eth.rsshpc1.sc1.science.roche.com
  ProxyJump ec-hub
  User kyriakik
  IdentityFile ~/.ssh/id_ed25519

Host *
  StrictHostKeyChecking no
  UserKnownHostsFile /dev/null
  ForwardAgent yes
```

#### Which host to use

- **ec-hub** — jump host with access to the shared filesystem (`/sc1/...`) and Slurm (`sbatch`, `squeue`, `scancel`). Use for all operations: running commands, submitting jobs, copying files.
- **sc1** (lb022dev) — direct login node via ProxyJump through ec-hub. May be unreachable at times — prefer ec-hub.

#### Common operations

```bash
# Interactive shell
ssh ec-hub

# Run a remote command
ssh ec-hub 'ls /sc1/groups/sbx/workspace/kyriakik/data/tools/'

# Copy files to cluster
scp build_release/Executable/dinara ec-hub:/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36X

# Copy files from cluster
scp ec-hub:/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_ins_eval/dinara_v36X_ins_sorted.vcf.gz /tmp/

# Submit a Slurm job
ssh ec-hub 'cd /sc1/groups/sbx/workspace/kyriakik/structural_variants/full_ins_eval && sbatch --array=1-26101 run_eval_v36X.sh'

# Check job status
ssh ec-hub 'squeue -u kyriakik'

# Cancel jobs
ssh ec-hub 'scancel -u kyriakik'
```

#### Available modules on ec-hub

Load bioinformatics tools with `module load`:

```bash
module load BCFtools/1.21-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0    # provides bgzip, tabix
module load SAMtools/1.21-GCC-13.3.0
```

These must also be loaded inside Slurm job scripts if the job uses `bcftools`, `bgzip`, etc.

#### Filesystem layout

All data lives on the shared filesystem under `/sc1/groups/sbx/workspace/kyriakik/`:

```
structural_variants/
├── full_del_eval/          # HG002 germline DEL benchmark
│   ├── cases/              # Extracted per-case dirs (reference.fa, reads.fa, region.bam)
│   ├── results_v36{X}/     # Per-version log files
│   ├── del_cases.txt       # Case list for Slurm array
│   └── logs_to_vcf.py      # Log → VCF conversion
├── full_ins_eval/          # HG002 germline INS benchmark
│   ├── cases/
│   ├── results_v36{X}/
│   ├── ins_cases.txt
│   ├── all_ins_positions.tsv
│   └── ins_logs_to_vcf.py
├── hg008t_del_eval/        # HG008-T somatic DEL benchmark
│   ├── cases/
│   └── del_cases.txt
├── hg008t_ins_eval/        # HG008-T somatic INS benchmark
│   ├── cases/
│   └── ins_cases.txt
├── somatic_hg008/          # Legacy somatic eval (pre-extraction)
└── extract_hg008t_cases.py # Case extraction script
data/
├── tools/                  # Dinara binaries (dinara_v36{d..u}_*)
└── truth/                  # Truth sets
    ├── GRCh38_HG2-T2TQ100-V1.1_stvar.filt.vcf.gz      # HG002 germline (ALT=* filtered)
    ├── GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed     # HG002 benchmark regions
    └── HG008T_somatic_V0.5/                             # HG008-T somatic truth set
        ├── GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz
        ├── GRCh38_HG008-T-V0.5_somatic-stvar_ALL.draftbenchmark.vcf.gz
        ├── GRCh38_HG008-T-V0.5_somatic-stvar-clonal.draftbenchmark.bed
        ├── GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed
        └── README.md
```

#### Troubleshooting

- **"Permission denied (publickey)"**: The SSH key (`~/.ssh/id_ed25519`) must be registered on the cluster. Check with `ssh-keygen -l -f ~/.ssh/id_ed25519`.
- **Hanging SSH**: ec-hub may be slow on first connect (host key verification). The config sets `StrictHostKeyChecking no` to avoid interactive prompts.
- **"module not found"**: Use exact module names with version (e.g., `BCFtools/1.21-GCC-13.3.0`, not `bcftools`).
- **Slurm job failures**: Check logs with `ssh ec-hub 'cat /sc1/.../results_v36X/{case}.log'`. Common issues: missing `--assemblyDirectory` (causes `DinaraRun` conflicts), `/tmp` not shared across nodes.
- **scp slow**: Large VCF files (>100MB) may take minutes. Use `rsync` for resumable transfers: `rsync -avP ec-hub:/sc1/.../ /tmp/`.

# Anchor Window Pipeline — Endpoint Anchor Design

## Two-Pass Inter-Window Edge Creation

The anchor graph constructor uses a two-pass approach for inter-window edges:

- **Pass 1 (endpoint edges):** For each backbone transition (`backbonePreviousWindow` / `backboneNextWindow`), create the single best-sharing edge between the two windows. Reserve the anchors used (+ RC mirrors) so pass 2 cannot reuse them.
- **Early trim:** After pass 1, disable backbone anchors beyond the endpoint positions. This constrains the graph to the region between endpoints.
- **Pass 2 (internal edges):** Create edges for all remaining inter-window anchor pairs, skipping reserved anchors.

## Endpoint Anchor Invariants

Each window has **at most 2 endpoint anchors** — one at the head (connecting to `backbonePreviousWindow`) and one at the tail (connecting to `backboneNextWindow`). Fewer if the window is at the start/end of a backbone chain.

Key properties:
1. **Never filtered.** All filters (singleEdge, bypass/detour, bubble, spur) check `isEndpointAnchorPrev || isEndpointAnchorNext` and skip those edges.
2. **Single edge per side.** Each endpoint anchor has exactly one endpoint edge, connecting to exactly one adjacent window.
3. **Highest sharing.** Pass 1 selects the anchor pair with the most shared reads between the two backbone reads of the adjacent windows.
4. **Internal edges are bounded.** After early trim, internal edges exist strictly between the two endpoint anchors of each window.

## Per-Anchor Endpoint Flags

Each edge carries two flags (not serialized):
- `isEndpointAnchorPrev` — source anchor is in the `endpointAnchors` set (reserved during pass 1).
- `isEndpointAnchorNext` — target anchor is in the `endpointAnchors` set.

These flags are set directly from the `endpointAnchors` set populated during pass 1 — no re-derivation needed. The flags drive both GFA tag output and filter protection.

## GFA Tags

Each link line gets `pw:Z:` and `nw:Z:` tags:
- `Endpoint` — anchor is in the `endpointAnchors` set
- `internal` — inter-window but not an endpoint anchor
- `intra` — both anchors in the same window

## Key Data Structures

- `endpointWindowPairs` — set of `{min, max}` normalized window pairs that correspond to backbone transitions. Used to identify which window pairs get endpoint edges in pass 1.
- `endpointAnchors` — set of anchor IDs (+ RC mirrors) reserved by pass 1 endpoint edges. Used for per-anchor flags and to prevent pass 2 from reusing these anchors.
- `filteredBackbonePositions` — per-window backbone chain: the longest subsequence of the backbone read's journey positions where every consecutive pair has sufficient read support. Defines the spine of each window's subgraph.

## Removed

- `recomputeBackboneEndpoints` — removed. `backbonePreviousWindow` / `backboneNextWindow` are set once during construction and do not change.
- Position-based endpoint detection — removed. Previously walked `filteredBackbonePositions` to find first/last active backbone anchor per window. Replaced by direct `endpointAnchors` lookup, which is the ground truth from pass 1.

## Files Changed

- `src/Shasta2AnchorGraph.hpp` — added `isEndpointAnchorPrev`, `isEndpointAnchorNext` to edge; added `endpointWindowPairs`, `endpointAnchors` to graph.
- `src/Shasta2AnchorGraph.cpp` — two-pass construction, early trim, `endpointAnchors`-based endpoint flags, updated filters, removed `recomputeBackboneEndpoints`.
- `src/Shasta2AnchorGraphGfa.cpp` — `writeGfa` accepts `anchorWindows`, emits `pw:Z:`/`nw:Z:` tags using per-anchor flags.
- `srcMain/main.cpp` — updated `writeGfa` call sites to pass `&anchorWindows`.
