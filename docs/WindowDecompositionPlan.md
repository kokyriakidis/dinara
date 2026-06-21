# Window Layout as a Doubled-Vertex String Graph

## What changed and why

An earlier version of this plan proposed splitting each window by its
`(prev, next)` traversal context ("prev-next paths"). Diagnostics on the test
assembly showed that premise does not hold for the current window model:

- **Full triplets are nearly absent.** In `WindowClasses.csv`, most windows have
  `FullClasses = 0`; the best window reaches `MaxFullClassReads = 5`. Reads do
  not pass *through* windows often enough to define `(prev, next)` contexts.
- **Isolated reads dominate.** The single `(noWindow, noWindow)` class holds the
  majority of reads in nearly every window (e.g. 65-80%).
- **Connectivity is pairwise, not path-shaped.** What reads *do* provide is
  window-to-window adjacency: "W is next to X." Left- and right-partial classes
  are symmetric per window (`LeftPartial == RightPartial`), which is the
  forward/RC mirror of the *same* physical connection.

The reason is by design: `computeAnchorWindowsClean` is greedy on base span, so a
large window claims a long backbone interval and **absorbs the very bridge
spans** that smaller windows would have used to pass reads through. The merged
genomic contexts the old plan wanted to separate are therefore *inside* each
window already - captured by per-window `readClusters` (MSA phasing) and
`alternatePaths` - not strung across windows as triplets.

So the open problem is not intra-window decomposition. It is **layout**: connect
the windows into a linear order, mirror the reverse-complement ones, and never
create a strand contact (an edge joining the `+` end of one locus to the `-` end
of another).

## Findings from the doubled-graph diagnostics

Running the diagnostics (below) on the test assembly (69 windows) established:

- **Layout already linearizes cleanly.** Twin-dedup of the doubled-vertex
  unitigs yields **11 canonical oriented paths from 22 unitigs, 0 orientation
  conflicts, 0 strand contacts** (`WindowLinearPaths.csv`,
  `WindowOrientation.csv`). The fw/rc separation works: each window resolves to a
  single orientation. This is the answer to "connect windows linearly, mirror
  RC, no contacts."
- **57 / 69 windows are linear**; the residual fragmentation into 11 paths is
  caused entirely by **8 branch windows** (+ 4 genuine tips: windows 28, 50, 67
  and one more).
- **Read-threading does not work and is not needed.** Almost no read spans a
  window (14 full triplets total), so matching a window's in-arc to its out-arc
  by shared reads fails on nearly every window, *including trivially linear
  ones* (`WindowThreading.csv`: only 4/68 thread). Linearization succeeds anyway
  because it relies on **adjacency uniqueness** (in-deg/out-deg == 1), not on
  any single read spanning the window. Threading is a dead end; do not use it.

### The branch windows are one root cause, not two

`WindowBranches.csv` classifies each branch by where its exits attach along the
window (`ExitSpreadFrac`):

- **Positional branches** (spread large): windows **0** (sides 0.92 / 0.97),
  **19+** (0.99), **5+** (0.66), **13+** (0.44). Their exits are spread along the
  whole window length - the window **absorbed multiple sub-loci** (greedy
  over-claiming).
- **read-set branches** (spread ~0): windows 9+, 19-, 31-, 37+, 41+. Every one of
  these attaches at fraction ~1.0 **into window 0 or window 19** (targets `0-`,
  `19-`).

These are not independent. The read-set "forks" are other windows connecting to
**different internal positions of the over-sized windows 0 and 19**. From a small
window's view, "connects to 0- or 67+" looks like a fork at its end, but it is
really that small window attaching to one sub-locus of window 0 while another
window continues from a different sub-locus.

**Consequence:** splitting the over-claimed windows (0, 19, and likely 5, 13)
*positionally* at their exit-attachment points resolves the read-set branches
too - the single hub becomes a chain `0a -> 0b -> 0c`, and each former fork now
targets a distinct sub-window. One fix, not two.

### Fix location

The true root cause is upstream: `computeAnchorWindowsClean` over-claims base
span when building the largest windows. Two options:

1. **Upstream** - cap window span / break at coverage discontinuities in the
   partitioner. Cleaner but riskier (the partitioner feeds everything).
2. **Post-hoc positional split** - split the few positional windows at their
   attachment fractions before linearization. Localized, reversible, behind a
   toggle.

Recommendation: prototype with **(2)** to confirm positional splitting collapses
the 11 paths toward the expected ~2-4 contigs, then fold the logic upstream into
**(1)**.

## The model: borrow hifiasm's invariant

hifiasm does not *resolve* strand contacts on its string graph - it makes them
**unrepresentable**. Strand consistency is a property of the data structure, not
a cleanup pass. Dinara already mirrors this for reads in `StringGraph`
(`src/StringGraph.hpp`, `src/AssemblerStringGraph.cpp`). We apply the same model
to windows.

### Doubled vertices

Every window `W` becomes two oriented vertices:

```
v         = (W << 1) | strand       strand 0 = forward, 1 = RC mirror
window(v) = v >> 1
twin(v)   = v ^ 1                    same window, opposite orientation
```

This is exactly the forward / RC-mirror pair already tracked in
`anchorToWindow`: a forward window is `windowId`, its RC mirror is
`windowId + windowCount`. The strand bit *is* the orientation; there is no
separate edge attribute to get wrong.

### Twin arcs

A read that visits oriented window `(W, sW)` then `(X, sX)` emits **two** arcs:

```
arc:   vertex(W,sW) -> vertex(X,sX)
twin:  vertex(X,sX)^1 -> vertex(W,sW)^1
```

The twin is the same connection read in the opposite direction. Storing both
makes the graph symmetric by construction and accounts for the observed
left/right partial symmetry directly. A `+`-to-`-` strand contact cannot be
written down: orientation is carried by the vertex id, and the RC of any walk
`v0 -> ... -> vk` is the guaranteed-present twin walk `vk^1 -> ... -> v0^1`.

### Linearization = unitig extension

With the doubled graph in place, layout is unitig construction - the analogue of
hifiasm's `asg_extend` (Dinara: `stringGraphExtend`,
`src/AssemblerStringGraphClean.cpp:79`). A step `a -> b` is *mergeable* iff `a`
has a unique outgoing arc and `b` has a unique incoming arc. Follow mergeable
steps to grow a maximal linear run. A window chosen RC is emitted from its
strand-1 vertex, i.e. with `anchorId ^ 1` anchors in reversed backbone order -
which the existing `anchorToWindow` mirror convention already encodes.

No spanning-tree + BFS orientation pass is needed. That technique is for an
undirected *signed* graph; the doubled representation removes the sign entirely.

## Existing substrate (no new data collection)

- **`anchorToWindow`** (`WindowTransitions.cpp`) - forward windows as
  `windowId`, RC mirrors as `windowId + windowCount`. This is the doubled
  vertex set.
- **`Shasta2Journeys`** - per-oriented-read ordered anchor sequences. Walking a
  journey and reading each anchor's `anchorToWindow` entry yields the oriented
  window visit sequence, hence the arcs.
- **`AnchorWindow::backboneOrientedReadId` / `filteredBackbonePositions`** - the
  per-window anchor chain to emit once orientation is chosen.
- **Per-window `readClusters` / `alternatePaths`** - already hold the
  intra-window haplotype structure; layout is orthogonal to them.

## Diagnostic (implemented, read-only)

`computeWindowTransitions` (`src/WindowTransitions.cpp`) now emits, behind
`windowArcDiagnostic`:

- **`WindowArcs.csv`** - one row per directed arc:
  `FromWindow, FromStrand(+/-), ToWindow, ToStrand(+/-), Reads, TwinReads,
  TwinPresent`. Lets you confirm twin symmetry and inspect edge support.
- **`WindowUnitigs.csv`** - one row per maximal unitig:
  `UnitigId, WindowCount, WindowPath` (oriented windows, e.g. `12+ 7- 31+`).
- **Console summary** - oriented vertex count, used arcs, unitig count, longest
  unitig (in windows), branch vertices, tips, and `twin-missing` arcs. The
  `twin-missing` count is the strand-contact check: it must be 0 by
  construction.

Read these to answer, before any graph rewrite:

1. **How linear is the layout?** longest unitig / total windows, and the number
   of branch vertices. Few branches and long unitigs mean the windows already
   lay out as near-linear contigs.
2. **Is the representation sound?** `twin-missing` must be 0.
3. **Where are the tangles?** branch vertices are the windows that need
   intra-window separation (`readClusters` / `alternatePaths`) before they can
   be laid out unambiguously.

## Implementation plan

1. **Diagnostic (done).** Doubled-vertex arc list + unitig extension, read-only,
   behind a toggle in `WindowTransitions.cpp`.
2. **Construction (behind a toggle).** Build the doubled window string graph as a
   first-class structure and run hifiasm-style cleanup on it: symmetric arc
   deletion, tip cutting, and unitig assembly - reusing the existing
   `StringGraph` machinery with windows as vertices instead of reads.
3. **Emit oriented window paths.** For each unitig, concatenate each window's
   backbone chain in the chosen orientation (forward anchors, or `anchorId ^ 1`
   reversed for RC), stitched at the boundaries.
4. **Verify** on the test assembly: `twin-missing == 0`, total length preserved,
   no self-RC crossings, unitig count and N50 vs. the current graph.

## Risks / open questions

- **Branch resolution.** Where a window has out-degree > 1, layout is ambiguous.
  This is exactly where intra-window structure (`readClusters`,
  `alternatePaths`) must split the window into per-haplotype copies *before*
  layout. Sequencing: separate inside, then lay out.
- **Arc support threshold.** `minArcReads` (currently 1) trades connectivity for
  noise. The diagnostic reports arc support so the threshold can be chosen from
  data.
- **Cycles.** Repeats can make the window graph non-acyclic; unitig extension
  stops at revisits (handled), but contig extraction across cycles needs the
  same tip/bubble handling hifiasm applies to reads.
- **Coordinate with detangling.** `findDetourWindowPairs`
  (`src/Shasta2AnchorGraph.cpp`) already reasons about cross-window detours;
  the doubled-vertex layout should subsume or feed it, not duplicate it.

## Files touched

- `src/WindowTransitions.cpp` - diagnostic (done); doubled-graph builder (step 2).
- `src/Shasta2AnchorGraph.cpp` - consume oriented unitig paths in construction
  (steps 2/3), behind toggle; reconcile with `findDetourWindowPairs`.
- Reuse `src/StringGraph.*` / `src/AssemblerStringGraphClean.cpp` cleanup
  helpers with windows as vertices where practical.
