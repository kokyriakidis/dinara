# Window Decomposition into Prev-Next Paths

## Motivation

In the disjoint-core model each window `W` is a single backbone anchor chain,
shared by every read that passes through it. But reads enter `W` from different
previous windows and leave to different next windows. Collapsing all of them
into one chain merges distinct genomic contexts (haplotypes, repeat copies,
strands) into a single shared segment, which is exactly where the graph tangles.

**Decomposition** replaces each window with one anchor path *per traversal
context*. A read group that enters from `A` and exits to `B` gets its own copy
of the window's interior, supported by that group's reads. Read-coherent threads
are kept separate instead of being forced through a shared bottleneck.

This subsumes strand separation: because reads are oriented, forward and RC
traffic fall into different `(prev, next)` classes automatically, so the two
strands separate as a *consequence* of following reads — no special-case guard.

## The Existing Substrate

All inputs are already computed, independently of inter-window edges:

- **`AnchorWindow::transitionReads`** — `map<(prev,next), vector<OrientedReadId>>`
  (populated in `WindowTransitions.cpp:141`). For window `W`, every observed
  `(previousWindow, nextWindow)` route and the reads that follow it. `noWindow`
  marks a read that starts/ends inside `W`.
- **`Shasta2AnchorGraph::readWindows`** — `map<read, vector<windowId>>`, each
  read's ordered window path (the trace to follow).
- **`AnchorWindow::readIntervals`** + the `windowReadIndex` lookup
  (`WindowTransitions.cpp:97`) — for a `(window, read)` pair, the read's
  `[begin, end)` journey interval inside the window.
- **`AnchorWindow::backboneOrientedReadId` / `filteredBackbonePositions`** — the
  window's canonical anchor chain and ordering frame.

Decomposition is therefore a **rewrite of window→graph construction**, not new
data collection.

## Core Model

For each window `W`, for each `(prev, next)` class in `W.transitionReads`:

1. Gather the reads in that class.
2. Within `W`, reconstruct the ordered anchor sub-path those reads traverse
   (their agreed sequence of `W`'s anchors, in backbone order).
3. Emit that sub-path as a distinct chain — a **prev-next path** `P(W, prev, next)`.
4. Stitch paths across windows: the tail of `P(W, prev, next)` connects to the
   head of `P(next, W, nextnext)` for reads that continue, and the head of
   `P(W, prev, next)` connects to the tail of `P(prev, prevprev, W)`.

The result is a graph woven from read-coherent threads. A locus that is
traversed in multiple contexts gets multiple copies (one per context) rather
than one collapsed window.

## Decisions (resolved)

- **Support threshold: any support >= 1.** Every observed `(prev, next)` class
  becomes a path; a single full-triplet read is enough to instantiate one.
- **Partial reads (`prev==noWindow` or `next==noWindow`): fold into matching
  full paths, fall back to terminal paths.**
  - A `(noWindow, B)` read started inside `W` only because of where the read
    begins, not biology. It is real evidence for the `W→B` exit and is
    compatible with any full `(X, B)` path that shares its sub-path through `W`.
  - **Rule:** a partial folds into (contributes coverage to) any anchor-compatible
    full path. Only if *no* full path matches does it become its own terminal
    path — which is then a genuine contig end (true assembly tip), not a
    per-read artifact.
  - Discriminator: "does any read carry a full triplet through this side of
    `W`?" Yes → absorb partials; No → the partial is the boundary.

## Algorithm Sketch

```
for each window W:
    classes = W.transitionReads                      # (prev,next) -> reads
    fullClasses = { (p,n): reads | p != noW and n != noW }
    leftPartials  = { (noW,n): reads }
    rightPartials = { (p,noW): reads }

    # 1. Build a path per full class.
    paths = []
    for (prev, next), reads in fullClasses:
        subpath = agreedAnchorSubpath(W, reads)      # ordered W anchors
        paths.append(Path(W, prev, next, subpath, reads))

    # 2. Fold partials into compatible full paths; else make terminals.
    for (noW, n), reads in leftPartials:
        matches = [p for p in paths if p.next == n and compatible(p, reads)]
        if matches: addCoverage(best(matches), reads)
        else:       paths.append(Path(W, noW, n, agreedAnchorSubpath(W,reads), reads))
    # symmetric for rightPartials

    emit(paths)

# 3. Stitch across windows by shared (window, boundary-window) endpoints.
for each read r:
    seq = readWindows[r]                              # ordered windows
    for consecutive (W, X) in seq:
        connect tail of P(W, ?, X)  ->  head of P(X, W, ?)
        choosing the path copies whose prev/next match r's actual neighbors
```

`agreedAnchorSubpath(W, reads)` reuses the within-window ordering already used
for the backbone chain (`filteredBackbonePositions` / backbone-position order),
intersected with the anchors the class's reads actually carry.

## Anchor Ownership Question (open)

Today an anchor belongs to exactly one window (`anchorToWindow`), and the graph
is built over global `anchorId` vertices. Decomposition wants the *same* locus
to appear in multiple path copies. Two options:

- **(D1) Shared anchor vertices, multiple edges.** Keep one vertex per anchor;
  each path copy is a distinct *edge route* over shared vertices. Simpler graph,
  but copies of a locus still share a vertex — partial separation only.
- **(D2) Path-local anchor instances.** Give each path copy its own vertices
  (anchor instance = `(anchorId, pathId)`). True separation; larger graph; needs
  an instance→anchor map for sequence/consensus and for `check()`.

Recommendation: prototype with **D1** (cheap, reuses `addEdgeIfValid` and the
existing bidirected fold), measure how much tangling remains, then move to **D2**
only if shared vertices still merge contexts.

## Implementation Plan

1. **Diagnostic (read-only, behind a toggle).** Per window, report: number of
   `(prev,next)` classes, reads per class, full vs partial split, and whether
   each class's reads agree on a sub-path. Validates the premise (do classes
   actually carve clean threads?) before any graph change. Output to CSV.
2. **D1 construction (behind a toggle).** Replace single intra-window chain +
   inter-window edges with per-class path emission and cross-window stitching,
   over shared anchor vertices. Keep the current path available via toggle.
3. **Verify** on the test assembly: total length preserved, strands separated
   (no self-RC crossings), tangles reduced vs. the shared-window graph.
4. **D2 (only if needed).** Path-local anchor instances for full separation.

## Risks / Open Questions

- **Sub-path agreement.** Reads in a class may not perfectly agree on `W`'s
  interior (indels, missed anchors). Need a consensus/LIS step like the existing
  backbone filtering; "agreed sub-path" must be defined precisely.
- **Stitching ambiguity.** When `W` has multiple copies for boundary `X`,
  choosing which copy a read's `W→X` transition connects to requires matching
  the read's *full* neighbor context, not just the adjacent window.
- **Coverage fragmentation at support 1.** Support >= 1 maximizes separation but
  admits noise; the diagnostic (step 1) should quantify how many singleton
  classes exist before committing.
- **Downstream contracts.** `addConfidentBridges`, `toBidirected`, and `check()`
  assume one-vertex-per-anchor; D2 would require touching those.

## Files Touched (anticipated)

- `src/WindowTransitions.cpp` / `.hpp` — already produces `transitionReads`;
  may add the diagnostic here.
- `src/Shasta2AnchorGraph.cpp` — construction rewrite (steps 2/4), behind toggle.
- `src/AnchorWindows.hpp` — possible per-path/per-instance structures (D2).
- `srcMain/main.cpp` — wire the diagnostic and the toggle.
