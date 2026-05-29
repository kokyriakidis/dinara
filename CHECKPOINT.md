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
5. **Build anchor graph** — `Shasta2AnchorGraph` constructor from anchor windows
6. **Detangle** — `detangleWindows()` splits backbone anchors of tangled windows, then rebuilds the graph with parallel chains per path
7. **Export for shasta2** — external anchors (`Shasta2ExternalAnchors`) and anchor graph (`Shasta2ExternalAnchorGraph`) in shasta2-native binary format
8. **Early return** — post-graph steps (transitive reduction, assembly graph) are currently disabled via `return;`

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

## Detangling (`src/DinaraDetangle.hpp/cpp`)

After the initial anchor graph is built (including Rule 1 and `transitionReads` population), detangling splits backbone anchors of tangled windows so that each path through the window gets its own anchor copies.

### When a Window Is Tangled

A window is a detangling candidate if it has ≥ 2 distinct **through-flows** — `(prev, next)` pairs in `transitionReads` where both `prev ≠ noWindow` and `next ≠ noWindow`, each with ≥ `minInterWindowCoverage` reads.

### Algorithm

1. **Identify candidates**: Scan each window's `transitionReads` for through-flows.
2. **Partition reads by path**: Each through-flow `(prev, next)` defines a path. Reads are assigned to their path based on their `(prev, next)` pair.
3. **Split backbone anchors**: For each backbone anchor pair (canonical + RC) in a tangled window, create new anchor copies — one per path — containing only that path's reads. RC subsets are built by matching strand-flipped read IDs.
4. **Deferred appending**: All new anchors are collected first, then appended to `anchorMarkerInfos` in bulk to avoid invalidating spans during iteration.
5. **Build split map**: `anchorSplitMap[originalId] = [newId_path0, newId_path1, ...]`.
6. **Rebuild the graph**: The `Shasta2AnchorGraph` constructor accepts the split map. For windows whose backbone anchors appear in the map, parallel intra-window chains are created — one per path — using the new anchor IDs.

### Current Limitations

- Only handles through-flow reads (both `prev` and `next` are real windows). Reads that start or end at a tangled window are excluded from all split anchor copies.
- Runs once (no iterative detangling for cascading tangles).
- Uses `minInterWindowCoverage` as the flow threshold.

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

- Rule 1 (bounding span trimming) is active and working
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
| DEL 100-500bp | 10 | 8 | 2 | 0 | 80% |
| DEL 500-1000bp | 10 | 7 | 3 | 0 | 70% |
| DEL >1000bp | 10 | 6 | 2 | 2 | 60% |
| INS <100bp | 10 | 7 | 2 | 1 | 70% |
| INS 100-500bp | 10 | 1 | 6 | 3 | 10% |
| INS 500-1000bp | 10 | 5 | 3 | 2 | 50% |
| INS >1000bp | 10 | 0 | 5 | 5 | 0% |
| **TOTAL** | **80** | **43** | **24** | **13** | **54%** |

**Correct type (✅+⚠️): 67/80 = 84%**

Key observations:
- **DEL <100bp at 90%**: compound CIGAR merging (35D+55D→90D) and CIGAR-corroborated k-mer clusters fix tandem repeat and STR cases
- **DEL 100-500bp at 80%**: size-gated CIGAR clustering separates mixed-size deletion clusters (e.g. 57D vs 114D)
- **DEL 500-1000bp at 70%**: SA-tag VNTR suppression relaxed for strong calls (≥10 reads); coverage-drop corroborated k-mer clusters for marker-depleted regions
- **INS detection limited**: 13/40 exact (33%), 27/40 correct type (68%)
- **INS 100-500bp** remains weak (10% exact): insertions exceed read length and soft-clip assembly can only partially span them
- **INS >1000bp** has 0% exact: assembly contigs max out at ~600bp with 250bp reads
- **Type confusion**: ~8 INS cases still detected as DEL in tandem repeats

### Known Limitations

1. **bw=100** prevents single-chain deletion detection >100bp; split-read detection handles larger deletions
2. **VNTR insertions** are fundamentally limited: short reads can't span VNTRs, k-mers are non-unique, BAM extraction depletes VNTR reads
3. **Microsatellite insertions** (e.g., AT-repeats): reads are too short to span from unique flanking sequence past the breakpoint
4. **Highly repetitive regions** (>90% non-unique 10-mers): no reliable anchoring possible with k=10

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

`dedca41` on `main` — "Add SA tag parsing, flank-gap analysis, and VNTR-depth deletion detection"

**Changes in this commit (545 insertions, 5 deletions):**
- Flank-gap analysis for marker-depleted coverage-drop regions
- VNTR-depth deletion detection (negative depth → DEL call, guarded: ≤30% of VNTR length, ≥50bp)
- Covdrop-indirect insertion detection
- SA tag parsing from BAM via htslib (`parseSaTagSvCalls()`)
- `--bam` CLI option for SA tag input
- maxK increased from 60 to 62
- htslib link dependencies added to CMakeLists.txt

### Git Info

- **Commit authorship**: `kokyriakidis <kokyriakidis@gmail.com>` — no Co-authored-by trailer
- **Branch**: `main`

### Parameters

- k=10, w=6 (minimizer window; w=6 works best, w=10 causes regressions in INS_130bp and DEL_164bp)
- SDUST: T=20, W=64
- maxK=62 for adaptive multi-k gap filling
- minimap2Bw=100, minimap2MaxGap=5000, chainingMode=1

### SV Detection Pipeline Layers (complete, in order)

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
