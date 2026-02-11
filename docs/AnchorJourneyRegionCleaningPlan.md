# Anchor-Journey Region Seeding, Graph Cleaning, Bubble Detection, Diploid Phasing (ONT)

## Scope and Goal
We already compute `journeys`: for each oriented read, an ordered list of anchors it visits.
The local anchor graph (derived from adjacency in journeys) is currently messy (repeats, spurious edges, tips, wrong joins).

This document proposes an implementation plan that uses the fact that each read is a linear witness to:
1. Define *candidate continuous genomic regions* (in anchor space) from read-journey terminals.
2. Induce a region-restricted subgraph from supporting reads.
3. Clean that subgraph (remove wrong branches, tips, transitive edges).
4. Detect bubbles and phase them assuming *diploid* structure (two branches).

We target ONT data; thresholds should be robust to higher error rates and occasional chimeras.


## Backbone-First: Use The Longest Read To Define A Global Region
This plan supports a **backbone-first** workflow:
1. Pick a (near-)longest read as a backbone.
2. Treat its journey as an ordered anchor scaffold for a "global genomic region" within a connected component.
3. Recruit other reads that share anchors with the backbone and use them to clean the local anchor graph and determine:
   - **homozygous** region: one dominant linear path, no stable two-branch structure.
   - **heterozygous** region: stable 2-branch bubbles with consistent read partition.
   - **misjoin / wrong connection**: junctions on the backbone not supported by enough bridging reads, with inconsistent alternatives.

Why this works:
- A read journey is a linear witness; the longest reads can span far enough to define long-range context.

Risks (ONT-specific) and how we control them:
- Long reads can be chimeric or repeat-traversing, so **the backbone defines a candidate region**, but every backbone adjacency must be validated by *multiple* supporting reads.

Backbone policy:
- Use the longest read **per connected component** (or top `B` longest) to propose a backbone.
- Break the backbone into sub-regions at "low confidence junctions" (defined below), rather than forcing one read to define a single uninterrupted region.


## Definitions
- `journey(r) = [A0, A1, ..., Ak]` for oriented read `r`.
- `terminal window size n`: number of anchors taken from the start and end of a journey.
- `seed`: a (canonical) pair of anchors `(p, s)` where `p` appears in the first `n` anchors and `s` appears in the last `n` anchors of many reads.
- `region`: the anchor set and adjacency edges induced by reads supporting a seed, restricted to the in-read interval between `p` and `s`.
- `region graph`: a directed adjacency graph on anchors inside the region.
- `backbone read`: a chosen long read whose journey defines the scaffold order for a candidate region.
- `backbone path`: the ordered anchor list from the backbone journey (possibly split into sub-regions).

Canonicalization:
- Anchors have reverse-complement ids; treat `(p, s)` and `(rc(s), rc(p))` as the same physical seed.


## High-Level Architecture
Two compatible entry points:
1. **Backbone-first (primary for this doc):** longest read defines a candidate region; then we induce/clean/phase using recruited reads.
2. **Seed-first (optional):** terminal co-occurrence seeds define regions without choosing a single backbone read first.

Both build the same evidence layers:
- **Adjacency evidence (local):** edges from consecutive anchors in journeys.
- **Terminal/long-range evidence:** links from start-window anchors to end-window anchors per read.

Cleaning and phasing are run on **region-induced subgraphs** rather than the entire messy graph.


## Phase 0: Instrumentation and Minimal Data Products
Add (or reuse existing) CSV/TSV outputs to support iteration:
- `AdjacencyEdges.csv`: `A,B,support,medianOffset,madOffset,rcSupport,qualityScore`.
- `TerminalPairs.csv`: `p,s,support,uniquePartnerCountP,uniquePartnerCountS,entropyP,entropyS`.
- `SeedReads.csv`: per seed `(p,s)`, list of supporting reads (or counts per component).
- `RegionSummary.csv`: per region, `seed,supportingReads,anchorCount,edgeCount,branchNodesCount`.
- `BubbleReport.csv`: per bubble, `entrance,exit,branch1Support,branch2Support,classification`.
- `PhaseBlocks.csv`: per region, phase block boundaries and read-cluster sizes.

This phase should not change assembly behavior; it is for observability.


## Phase 1: Backbone Selection and Backbone Confidence Profiling
Inputs: `journeys`, read lengths, and (if available) read-level chimera/contained flags.

1. Partition work by connected component (anchors x reads bipartite component).
2. Choose backbone candidates:
   - Start with the longest read in the component (strand 0), optionally keep top `B` for fallback.
3. Backbone journey preprocessing:
   - Collapse immediate duplicates if they occur (repeat noise).
   - Annotate each backbone anchor with:
     - coverage (anchor coverage),
     - degree in the raw adjacency graph (repeat-ness proxy),
     - partner entropy (optional; see terminal co-occurrence phase).
4. Define **backbone junctions**: each consecutive pair `(Ai -> Ai+1)` on the backbone.
5. For each junction, compute a confidence score from initial evidence:
   - `support`: number of other reads that contain `Ai` immediately followed by `Ai+1` (same direction) anywhere in their journey.
   - `rcSupport`: number of reads supporting `rc(Ai+1) -> rc(Ai)`.
   - `altOut(Ai)`: strength of best alternative outgoing edge from `Ai` (not `Ai+1`).
   - `altIn(Ai+1)`: strength of best alternative incoming edge to `Ai+1` (not `Ai`).

Low confidence junction heuristic (tunable):
- `support` is low relative to anchor coverage, OR
- strong competing alternatives exist (`altOut` close to `support`), OR
- poor RC symmetry.

Outcome:
- Split the backbone path into **sub-regions** separated at low-confidence junctions.
- Each sub-region is a candidate "global region" to be cleaned/phased.


## Phase 2: Read Recruitment Around A Backbone Sub-Region
Goal: collect reads that truly span the sub-region so they can validate edges and reveal diploid structure.

For a backbone sub-region defined by anchors `[Astart ... Aend]`:
1. Recruit reads with enough overlap to be informative:
   - They share at least `mShared` anchors with the backbone sub-region, AND
   - Their shared anchors are spread across the region (not all clustered in a small area).
2. Prefer "bridging" reads:
   - Reads that contain anchors near both ends of the sub-region (within a window).
3. For each recruited read, define its **projected interval** on the backbone:
   - Identify shared anchors and their positions in both journeys.
   - Ensure order-consistency (shared anchors appear in mostly the same order).
   - If order is inconsistent, classify the read as repeat-traversing/noisy for this region (down-weight or drop).

Outcome:
- A set of reads with a backbone-projected span and an extracted journey interval relevant to the region.


## Phase 3: Induce The Region Subgraph (Backbone-Conditioned)
From recruited reads:
1. Extract each read’s journey slice corresponding to its projected interval on the backbone.
2. Add adjacency edges `(X -> Y)` from consecutive anchors in these slices.
3. Store per-edge:
   - support,
   - offset samples (markers/bases) if available,
   - RC symmetry info.

This produces a region-induced adjacency multigraph strongly conditioned on the backbone region.


## Phase 4: Clean The Region Graph (ONT, Diploid-Aware)
Use the pruning rules from the original plan (support, coherence, RC reciprocity, competition, tips, transitive reduction),
but compute all statistics **within the recruited-read set**.

Key backbone-specific check:
- For each backbone junction `(Ai -> Ai+1)`, require a minimum number of bridging reads that support it.
- If a junction fails, treat it as a **candidate misjoin** and split the region there.


## Phase 5: Decide Homozygous vs Heterozygous vs Misjoin
We make this decision on cleaned graphs and read traversal patterns.

### 5.1 Homozygous Region Signature
- Cleaned region graph is mostly linear: most nodes have `inDegree<=1` and `outDegree<=1`.
- Branch competition prunes alternatives strongly (dominant single edge everywhere).
- Recruited reads mostly agree on a single traversal.

### 5.2 Heterozygous Region Signature (Diploid)
- Multiple stable 2-branch bubbles exist:
  - entrance has 2 strong outgoing edges,
  - branches reconverge to a common exit,
  - both branches have coherent offsets and sufficient support.
- Reads partition consistently by bubble traversal (two read clusters).
- Along the backbone order, bubble choices are consistent within each cluster (phase blocks).

### 5.3 Misjoin / Wrong Connection Signature
Backbone contains one or more junctions where:
- bridging support is insufficient (few reads span across the junction),
- strong inconsistent alternatives exist (high `altOut` or `altIn`),
- recruited reads that span the junction disagree on ordering or offsets,
- long-range/terminal evidence contradicts the backbone continuation.

Action:
- Cut the region at misjoin junctions.
- Optionally re-run backbone selection within the resulting subcomponents.


## Phase 6: Diploid Phasing Along The Backbone
1. Detect bubbles in the cleaned region graph.
2. For each recruited read, record bubble-branch choices as a sparse binary signature.
3. Cluster reads into 2 haplotypes by agreement on signatures.
4. Assign bubble branches to haplotypes.
5. Emit phase blocks along backbone order; break blocks at:
   - complex regions (>2-branch),
   - low support,
   - inconsistent clustering.


## Phase 7: Optional Terminal Co-Occurrence (Seed-First) To Stabilize/Extend
Terminal co-occurrence remains useful to:
- propose alternate backbones if the longest read is chimeric,
- connect adjacent backbone sub-regions,
- stabilize across sparse anchor stretches.


## Phase X: Build Terminal Co-Occurrence Table (Optional Support Layer)
Inputs: `journeys`, read lengths, RC mapping.

For each physical read (use strand 0 to avoid duplicates):
1. If `|journey| < minJourneyAnchors`, skip for seeding (still allow it for adjacency support).
2. Let `P = first n anchors`, `S = last n anchors`.
3. For each `(p in P, s in S)` increment `count(p,s)`:
   - Apply rank weights: `w = wStart(rank(p)) * wEnd(rank(s)) * wLen(readLength)`.
   - Suggested weights:
     - `wStart(i) = exp(-i / tau)` with small `tau` (e.g. 1.0)
     - `wEnd(i) = exp(-i / tau)`
     - `wLen(L) = clamp(L / L0, 0.5, 2.0)` (avoid extreme domination by ultra-long reads)
4. Canonicalize `(p,s)` with RC: map to `min((p,s), (rc(s), rc(p)))`.

Post-processing on `(p,s)` table:
- For each anchor `p`, compute:
  - number of distinct `s` partners above a tiny support threshold (repeat-ness proxy),
  - partner entropy `H(p) = -sum q(s|p) log q(s|p)` where `q` is normalized counts.
- Same for `s`.

Seed candidates are `(p,s)` pairs with:
- high support,
- low partner entropy for `p` and `s` (anchors behave uniquely at terminals),
- acceptable component consistency (most supporting reads in the same connected component).


## Phase Y: Choose Seeds and Induce Regions (Optional Alternative Entry)
Seed selection strategy:
- Choose top `K` seeds by support per connected component.
- Avoid near-duplicates: if two seeds share most supporting reads and induced anchor sets, keep one.

For a seed `(p,s)`:
1. Gather supporting reads: reads where `p` appears in first `n` and `s` appears in last `n` (same orientation).
2. For each supporting read:
   - choose an occurrence of `p` constrained to the first `n` positions of the journey,
   - choose an occurrence of `s` constrained to the last `n` positions of the journey,
   - ensure the chosen `p` position < chosen `s` position; otherwise discard or flip via RC canonicalization.
3. Extract the journey interval `[p..s]` (inclusive) as that read’s region witness.
4. Accumulate:
   - anchor set `R`,
   - adjacency edges `(Ai -> Ai+1)` within this interval.

Quality checks while inducing:
- Discard reads with:
  - huge internal gaps (anchor-to-anchor base gap above a percentile threshold),
  - inconsistent directionality (e.g. multiple inversions implied by ordinal offsets),
  - very high repeated anchor occurrences in the interval (repeat traversal).

Output: region-induced multigraph with per-edge support and offset samples.


## Phase 3: Clean the Region Graph (ONT, Diploid-Aware)
We want to remove:
- tips (short dead ends),
- spurious edges caused by repeats,
- transitive edges,
- wrong branches (non-diploid complexity).

### 3.1 Edge Scoring
For an edge `A->B` in the region:
- `support(A->B)` = number of region-supporting reads that show `A,B` consecutively.
- `offsetSamples(A->B)` = base or marker offsets observed on those reads.
- `offsetMedian`, `offsetMAD`.
- `rcSupportRatio` = support of `rc(B)->rc(A)` divided by support of `A->B` (bounded).
- `compositionScore` = corrected Jaccard (or similar) from existing anchor-pair similarity utilities, computed on-demand for high-support edges.

Define a simple score:
`score = log(1+support) + alpha*compositionScore - beta*log(1+offsetMAD) - gamma*rcPenalty`

### 3.2 Pruning Rules
Apply in this order:
1. **Support filter:** drop edges with `support < Smin(region)` (coverage-aware).
2. **Offset coherence:** drop edges with `offsetMAD > MADmax` or with obvious multimodality.
3. **RC reciprocity:** drop edges with `rcSupportRatio` outside `[rMin, rMax]` unless support is extremely high.
4. **Branch competition:**
   - For each node `A`, rank outgoing edges by score.
   - Keep best edge.
   - Keep second edge only if it passes a diploid allowance rule:
     - `score2 >= diploidFrac * score1` and `support2 >= supportFloor`.
   - Mark additional edges as candidates for removal (or for “complex region” labeling).
5. **Tip removal:** iteratively remove nodes/edges on short dead-end paths below length/support thresholds.
6. **Transitive reduction:** if `A->C` is explained by `A->B->C` on most supporting reads with coherent offsets, remove `A->C`.

The result should be a mostly-linear graph with occasional 2-branch bubbles.


## Phase 4: Bubble Detection (Superbubbles) and Classification
Within each cleaned region graph:
- Identify candidate entrances: nodes with `outDegree == 2`.
- Find exits where paths reconverge (classic superbubble detection).

For each candidate bubble `(entrance E, exit X)` with two branches:
- Compute branch support:
  - number of reads whose region interval traverses branch 1 vs branch 2.
- Check branch coherence:
  - offsets coherent within branch edges,
  - no extra branching inside branch (or limited).

Classification:
- **Likely diploid bubble** if both branches have:
  - sufficient support (not necessarily balanced),
  - coherent offsets,
  - and both branches reconnect to the same exit.
- **Likely wrong branch** if one branch:
  - is a tip (no reconvergence),
  - has low support or high offset dispersion,
  - or conflicts with terminal link evidence (next phase).


## Phase 5: Diploid Phasing Inside Regions
We use bubble traversals as a sparse binary signature.

For each read supporting a region:
- For each detected bubble it touches, record which branch it took (0/1).

Clustering:
- Cluster reads into two groups by maximizing agreement across bubble signatures (greedy + refinement is fine).
- Assign each branch to haplotype based on which cluster supports it.

Phase blocks:
- A phase block is a contiguous stretch of the region where bubble assignments are consistent.
- Break a block where:
  - reads frequently switch clusters,
  - bubble structure becomes complex (>2 branches),
  - or support drops sharply.

Use terminal co-occurrence as stabilizer:
- For a seed `(p,s)`, the region-supporting reads are already terminal-consistent.
- If a bubble branch is supported by reads that are not terminal-consistent with the seed, penalize that branch (possible error).


## Phase 6: Iterate Seeds, Merge Regions, and Expand
After cleaning/phasing a region:
- Expand seeds:
  - derive new terminal-like anchors from the region ends (in the cleaned graph),
  - look for neighboring seeds that overlap strongly, merge regions.
- Repeat per component until coverage saturates or complexity increases.

Stop conditions per component:
- no new high-confidence seeds,
- remaining graph labeled complex/repeat-rich,
- marginal gain in cleaned anchors per iteration is small.


## Phase 7: Success Criteria and Debugging Loops
Success criteria:
- In cleaned region graphs, most nodes have `inDegree<=1` and `outDegree<=1` except at bubbles.
- Bubble count is plausible and bubbles are mostly 2-branch with reconvergence.
- Phase blocks are long and have strong read support.

Debugging loop:
- If too many branches remain: tighten support/entropy filters for seeds and branch competition thresholds.
- If bubbles collapse incorrectly: relax diploid allowance or use more reads (increase seed support threshold).
- If regions are too small: increase `n`, allow more seeds per component, or down-weight entropy less aggressively.


## Initial Parameter Suggestions (ONT)
These are starting points only:
- `n` (terminal anchors per end): 3 (try 2..5).
- `minJourneyAnchors` for seeding: 8.
- seed support threshold: adaptive; start with `>= 10` reads (or weighted equivalent) per component.
- diploid allowance: keep 2nd outgoing edge if `support2 >= 0.25 * support1` and both pass coherence checks.
- tip removal: remove dead ends shorter than 3 anchors unless support is very high.

Backbone-specific (starting points):
- `B` backbone candidates per component: 1..3
- `mShared` anchors to recruit a read: 5 (or 3 if anchors are sparse)
- bridging requirement: read must share anchors within `w` of both ends of the backbone sub-region (try `w=10` anchors)
- misjoin junction: backbone edge support < `max(3, 0.05 * coverage(Ai))` OR strong competing alternative within 0.5x support


## Next Iteration Questions
1. What is typical ONT coverage for your runs (30x, 60x, 100x)?
2. Are chimeric reads already flagged/removed before journeys are computed?
3. Do you want to run this region cleaning before mode3 assembly, or as an analysis/diagnostics path first?
