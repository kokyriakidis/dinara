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

## Current State (Anchor Windows)

- Inter-window edge filter pipeline is active: shortcut, parallel, cross-window filters using backbone endpoints, with trimBackbones and recomputeBackboneEndpoints between each step
- All edge disabling uses `disableEdge()` member function for RC mirror symmetry
- Bypass edges and Case 2 detangle are disabled (`#if 0`)
- Detangling splits tangled windows by through-flow paths (through-flows only; start/end reads not yet assigned)
- `removeNegativeOffsets` is called at all 3 edge construction sites
- Edge verification runs after graph construction (0 backward edges expected)
- Export to shasta2 format is working with round-trip verification
- Shasta2 assembly completes successfully with dinara's exported anchors and graph
- Post-graph steps (transitive reduction, assembly graph, etc.) are disabled via `return;` in main.cpp

---

## SV Detection Pipeline (`--command svanchors`)

### Overview

Detects structural variants (insertions and deletions) from short reads aligned to a local reference. The pipeline uses k-mer-based chaining (hifiasm-derived DP) with minimap2-sr scoring to align reads, then applies multiple detection layers.

### Running

```bash
dinara --command svanchors \
  --reference reference.fa --input reads.fa \
  --assemblyDirectory outdir \
  --Reads.minReadLength 50 --Kmers.k 10 --Kmers.minimizerW 6
```

Test cases: `/tmp/sv_cases/` with insertion directories at top level and deletion directories under `DEL_medium_100_500bp/`.

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

Subsequent uncommitted changes (V36l):
- `src/Assembler.hpp` — `DepthScanDelCall` struct + `depthScanDelCalls()` method declaration
- `src/AssemblerSvAnchors.cpp` — `depthScanDelCalls()` implementation (~200 lines) + source priority table update (depth-scan-hom/het/sub at priorities 13-15)
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

*V36l recall on HG002 is unchanged (depth-scan adds calls but doesn't recover new TPs for germline). The depth-scan source is primarily impactful for somatic analysis (HG008-T: 81.7% → 95.0%).

### Current State (V36l)

- **Germline recall (HG002)**: 99.4% (10,112 / 10,173)
- **Somatic recall (HG008-T)**: 95.0% (57 / 60)
- **61 remaining germline misses**: Small DELs in large tandem repeat arrays where depth deficit is dominated by mappability loss. 0 cases with no calls at all. 100% detection sensitivity.
- **3 remaining somatic misses**: 1 too large (6.3Mb), 1 partial depth-scan (64kb), 1 position accuracy at 10% VAF
- **Avg calls/case**: ~15 (germline), ~20 (somatic, higher coverage)

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

**Solution**: Sort by source accuracy instead of read count. Priority order: merged-clusters > diagonal > early-CIGAR > SA-tag > split-read > cluster > flank-gap > kmer-journey > per-read-DEL > INV-cluster > path-based > depth-deficit-hom > depth-deficit-het > depth-scan-hom > depth-scan-het > depth-scan-sub > multi-k. Read count is tiebreaker within same source.

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

### Depth-Scan DEL Source (V36l)

Windowed BAM depth scanning to detect deletions from depth drops. Unlike depth-deficit (which integrates total deficit across the whole region), depth-scan identifies the specific contiguous region of low depth and estimates DEL size from its width.

**Algorithm:**
1. Open BAM, compute per-base depth across the full region in a single pass
2. Aggregate into adaptive windows (100bp for regions <10kb, 200bp for 10-50kb, 500bp for >50kb)
3. For each position, compute left-flanking median depth (up to 10 windows back)
4. Detect contiguous runs of windows where `depth / flank_depth ≤ threshold`
5. Scan at 3 sensitivity levels (thresholds: 0.65, 0.75, 0.85) with deduplication
6. Classify by depth ratio: `<0.05` → hom, `<0.55` → het, else → sub (subclonal)
7. Emit calls as `depth-scan-hom` / `depth-scan-het` / `depth-scan-sub`

**Key design decision:** Runs in `main.cpp` **before DP chaining** (right after early-CIGAR), so calls are emitted even when assembly times out on high-read-count regions. This is the primary motivation — large DELs (5kb+) often have too many reads for the DP chaining to complete within wall time.

**Implementation:** `Assembler::depthScanDelCalls()` in `src/AssemblerSvAnchors.cpp` (~200 lines). Uses htslib to read BAM directly. Called from `srcMain/main.cpp` after the early-CIGAR block.

**Limitations:**
- Window granularity limits size accuracy (±window_size). A 5.8kb DEL may be called as 5.6kb or 6.0kb.
- Cannot detect small DELs (<~400bp) because they don't span enough consecutive windows.
- Subclonal DELs at very low VAF (<10%) may not produce enough depth drop to trigger detection.
- Does not resolve breakpoints precisely — position accuracy is ±window_size.

**Somatic impact:** Recovers large DELs (5kb-64kb) that assembly cannot handle within wall time. See HG008-T evaluation below.

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
| V36l (+ depth-scan) | `dinara_v36k_depthscan` | **95.0%** | **57** | **3** |

### Failure Analysis

**7 cases timed out** during DP chaining (assembly never completed):
- All had large regions with many reads (5.8kb-64kb DELs, 50K-262K reads)
- All had clear depth signals (depth ratio 0.00-0.70)
- The depth-scan source recovers all 7 because it runs before DP chaining

**3 cases had wrong-size calls** from assembly:
- chr17:31977413 (27kb, VAF=40%) — assembly found 31bp, depth-scan found 27kb ✅
- chr7:6404812 (1.2kb, VAF=28%) — assembly found 43bp, depth-scan found 1.2kb ✅
- chr7:6408188 (62bp, VAF=26%) — assembly found 35+36bp split calls. Depth-scan can't resolve 62bp (below window granularity). Truvari matched via other means ✅

### Remaining 3 FNs

| Locus | Size | VAF | Issue |
|---|---|---|---|
| chr13:66842467 | 6.3Mb | NA | Too large for current pipeline (ignored) |
| chr11:58991970 | 64kb | 34% | Depth-scan found 25kb partial call; full 64kb not resolved due to internal depth variation |
| chr3:169770590 | 108bp | 10% | Position accuracy — closest call (89bp early-CIGAR) lands 672bp downstream. Very low VAF. |

### Somatic-Specific Observations

- Dinara detects somatic DELs down to **VAF ~7%** (chr6:72397103 at 6.6%, chr10:132826517 at 8.8%)
- The depth-scan source is essential for somatic analysis — large somatic DELs (5kb+) consistently time out during assembly
- The depth-scan `sub` (subclonal) threshold (ratio <0.85) catches het-like and subclonal depth drops that the existing depth-deficit source (which uses fixed flank sizes) may miss

### Cluster Paths

- **Cases**: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/somatic_hg008/cases/`
- **Results (assembly only)**: `.../somatic_hg008/results/`
- **Results (+ depth-scan)**: `.../somatic_hg008/results_dscan/`
- **Truth set**: `.../somatic_hg008/truth/`
- **Manifest**: `.../somatic_hg008/manifest.tsv`
- **VCFs**: `.../somatic_hg008/dinara_somatic_dels.vcf.gz` (assembly only), `.../somatic_hg008/dinara_dscan_dels.vcf.gz` (+ depth-scan)
- **Truvari**: `.../somatic_hg008/truvari_all/` (assembly only), `.../somatic_hg008/truvari_dscan/` (+ depth-scan)
- **Binary**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36k_depthscan`
- **Prepare script**: `.../somatic_hg008/prepare_cases_v2.sh`
- **Run script**: `.../somatic_hg008/run_depthscan.sh`

### Truvari Command (GIAB-recommended for HG008-T)

```bash
truvari bench \
  -b GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz \
  -c {calls}.vcf.gz \
  --reference GRCh38_GIABv3_...fasta \
  --includebed GRCh38_HG008-T-V0.5_somatic-stvar-clonal_and_subclonal.draftbenchmark.bed \
  --sizemax -1 --passonly --pick multi \
  -o {output}
```

---

### Cluster Paths and Binaries

- **Binaries**: `/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_v36{d..k}_*`, `dinara_v36k_depthscan`
- **Results**: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval/results_v36{d..k}/`
- **Somatic results**: `/sc1/groups/sbx/workspace/kyriakik/structural_variants/somatic_hg008/`
- **Analysis scripts**: `analyze.py`, `source_essentiality.py`, `analyze_multisource_v2.py`, `analyze_gap.py`, `analyze_adaptive.py` in `/sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval/`
- **SLURM template**: `run_eval_v36k.sh` — array job, 100 cases per task, 102 tasks

### Build and Deploy

```bash
# Build locally
cd /workspaces/dinara/build && make -j$(nproc)

# Deploy to cluster
scp build/Executable/dinara kyriakik@ec-hub.sc1.science.roche.com:/sc1/groups/sbx/workspace/kyriakik/data/tools/dinara_<version>

# Run eval (germline HG002)
ssh kyriakik@ec-hub.sc1.science.roche.com
cd /sc1/groups/sbx/workspace/kyriakik/structural_variants/full_del_eval
sbatch run_eval_<version>.sh

# Run eval (somatic HG008-T)
cd /sc1/groups/sbx/workspace/kyriakik/structural_variants/somatic_hg008
sbatch run_depthscan.sh

# Analyze
python3 analyze.py cases results_<version>
```

