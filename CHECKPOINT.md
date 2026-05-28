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

## Current State (Anchor Windows)

- Rule 1 (bounding span trimming) is active and working
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
3. **Remove non-unique reference k-mers** — blacklist k-mers appearing >1 time in the reference
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

### Key Files (SV Detection)

| File | Purpose |
|---|---|
| `src/AssemblerSvAnchors.cpp` | `buildSvMSA` — main SV detection engine |
| `src/Sdust.hpp` | Standalone SDUST low-complexity filter (Heng Li's algorithm) |
| `src/Assembler.hpp` | `AlignmentCandidatesInvertedIndexData` with chaining parameters |
| `src/AssemblerInvertedIndex.cpp` | Inverted index build, chaining DP, hit collection |
| `src/InvertedIndexBuilder.hpp` | Count-then-scatter index construction |
| `src/AssemblerMarkers.cpp` | `removeNonUniqueReferenceMarkers` — k-mer blacklisting |
| `srcMain/main.cpp` | svanchors pipeline orchestration (line ~3830+) |

### Test Results (18 cases)

| Case | True Size | Detected | Method |
|------|-----------|----------|--------|
| INS268 | 268bp | 278bp | path-based (1 hop) |
| INS254 | 254bp | 263bp | path-based (2 hops) |
| INS148 | 148bp | 149bp | path-based (1 hop) |
| DEL277 | 277bp | 277bp | split-read (14 reads) |
| DEL347 | 347bp | 328bp | split-read (9 reads) + DUST-STR |
| DEL324 | 324bp | 344bp | DEL CLUSTER |
| DEL160 | 160bp | 160bp | split-read (3 reads) |
| DEL137 | 137bp | 136bp | adaptive-bimodal |
| DEL182 | 182bp | 195bp | adaptive-bimodal |
| DEL119a | 119bp | 121bp | adaptive-bimodal |
| DEL119b | 119bp | 121bp | adaptive-bimodal |
| DEL147 | 234bp | 229bp | split-read (2 reads) |

**Undetected (known limitations):**

| Case | Root Cause |
|------|-----------|
| INS57 (57bp) | VNTR, 10% unique 10-mers at breakpoint |
| INS62 (62bp) | 1685bp VNTR (75× core motif), BAM-depleted |
| INS65 (65bp) | 1700bp VNTR, BAM-depleted |
| INS235 (235bp) | VNTR, 16.5x coverage (heavily depleted) |
| INS61 (61bp) | 2256bp AT-repeat microsatellite |
| DEL379 (379bp) | 98% repetitive 10-mers |

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
