# Gate 2 Design — Strand-Separated, Read-Corroborated Window Backbone

## Goal

Build the inter-window backbone so that a single bad/chimeric alignment cannot
collapse a locus with its own reverse complement, and ambiguous loci are
identified and deferred rather than guessed. Construction is **incremental and
evidence-gated**: contract what is unambiguous first, corroborate the rest, defer
the hard cases.

Phase 1 builds the clean backbone and **collects** (does not resolve) the
difficult places for Phase 2.

## Founding principles (locked with user)

1. A window is the biggest unfolded region supported by one backbone read; it is
   strand-pure internally. Strand contamination enters only through inter-window
   edges and through over-merged windows.
2. The two strands of a locus (window `W_fw` and its mirror `W_rc`) must remain
   two separate paths. An edge forcing a window's forward strand into the same
   component as its own reverse strand is the artifact.
3. We cannot pre-screen chimeras (only detectable by contradiction). Therefore no
   single read may create a backbone link: every union must be **corroborated**.
4. Start from the unambiguous interiors (one-in / one-out windows) — these are
   safe by definition. Whatever is left defines the difficult set.

## Confidence / corroboration gate (locked)

A window→window connection may be unioned into the backbone only if BOTH:

- **count** ≥ `minConnectionReads` (distinct physical reads, N) — chimera defense:
  a lone false join has count 1 and is rejected.
- **span product** ≥ `minConnectionSpanProduct` (T) — depth: cumulative
  `supportingSpanProduct` over the connection's transitions.

(N, T are tunable; start conservative, e.g. N=2, T from a span-product histogram.)
This gate admits **1b junction unions only**. It is NOT used to compute window
degree — 1a degree is pure topology over all distinct neighbor windows.

## Phase 1a — Contract the unambiguous interiors (safe, pure topology)

Degree counts **all distinct neighbor windows** (normalized strand), regardless of
how many anchor-edges or reads back each connection. A window is contracted only
if it is topologically clean — exactly one predecessor window and one successor
window. Any extra connection (even a single weak/chimeric edge) disqualifies the
window from the safe core and sends it to 1b / the difficult set. This is
deliberately conservative: under-contract rather than absorb topology that looks
clean but hides a stray edge.

(The count-and-span corroboration gate is NOT used here; it governs 1b unions
only. 1a is pure topology + guards.)

```
for each window W (forward strand units):
    inDeg(W)  = # distinct neighbor windows with an edge -> W   (all, normalized)
    outDeg(W) = # distinct neighbor windows with an edge W ->   (all, normalized)

for each window W with inDeg(W)==1 and outDeg(W)==1:
    if W is a GUARD violation: skip  (see guards below; -> difficult set)
    else:
        union(pred(W), W, succ(W))         # forced linear path — safe
        mirror-union the RC complements    # keep fw/rc exact reverse complements
```

Result: maximal linear components of strictly-(1,1) windows = the safe backbone
core. Equivalent to window-level `findLinearChains` with the strand-mirror
constraint, gated by the two guards.

### Guards applied in 1a (locked)

A (1,1) window is NOT contracted (diverted to the difficult set) if either:

- **strand conflict**: contracting it would put a window's fw strand in the same
  component as its own rc (`find(X) == find(complement(neighbor))`).
- **fw/rc same read**: the window contains both `oid` and `oid^1` of the SAME
  physical read among its supporting reads (self-fold / palindrome / collapse
  signature). Detected per window before contraction.

## Phase 1b — Join across junctions (corroborated, strand-guarded)

Reads are sorted **longest first** (exploration order = maximal total span; this
is the span metric realized at journey granularity). No chimera pre-screen —
corroboration does the screening.

```
# Accumulate corroboration over ALL reads first.
for each read R:
    project journey -> window chain w0..wk
    for consecutive (wi, wi+1):
        support[(wi,wi+1)].count += 1
        support[(wi,wi+1)].spanProduct += this transition's supportingSpanProduct

# Lay junctions, longest read first, union only corroborated + strand-safe steps.
for each read R (longest first):
    for consecutive (wi, wi+1):
        s = support[(wi,wi+1)]
        if s.count < N or s.spanProduct < T:          -> uncorroborated -> difficult
        elif find(wi) == find(complement(wi+1)):       -> strand conflict -> difficult
        elif wi or wi+1 is fw/rc-same-read window:      -> difficult
        else: union(wi, wi+1); mirror-union complements
```

A chimeric read contributes to `support` for its false step, but that step stays
at count 1 (no other read corroborates) -> never unioned. The chimera screens
itself out.

## Difficult set (Phase 2 input — preserved, never discarded)

Tagged by reason:
- `strand-conflict` (parity violation against established backbone)
- `fw-rc-same-read` (self-fold / palindrome / collapse)
- `uncorroborated` (count<N or span<T; includes chimeric joins)
- `branch` (legitimate degree>1 biology — distinguished in Phase 2)

Phase 2 (later) reclassifies high-confidence flank-coherent fold-backs as real
inverted repeats (explicit cross-strand links that do NOT merge components), and
discards low-confidence artifacts. Phase 1 must keep these intact (with their
transitions).

## Strand mechanism

Units are `(window, strand)` = raw window IDs `[0, 2*windowCount)`, IDs ≥
`windowCount` = rc. `complement(w) = (w>=windowCount)? w-windowCount : w+windowCount`
(existing `rcWindow` lambda). DSU over 2N units + mirror-union keeps fw and rc
graphs as exact reverse complements; "parity conflict" surfaces as
`find(X) == find(complement(Y))`. Reuses `dset64`.

## Data available (verified)

- `windowPairTransitions: map<pair<uint32_t,uint32_t>, vector<ReadTransition>>`,
  window IDs encode strand.
- `ReadTransition.supportingSpanProduct`, `supportingSpanA/B`, `oidValue`.
- fw-fw/rc-rc/fw-rc/rc-fw classification + RC-mirror canonicalization assert.
- `dset64` / `Shasta2DisjointSets`.

## Where it plugs in

`Shasta2AnchorGraph` constructor inter-window edge creation
(`src/Shasta2AnchorGraph.cpp` ~626–672). Replaces admission: build candidates,
run 1a then 1b, create edges only for unioned connections (reuse existing
best-pair + `addEdgeIfValid`). Intra-window backbone untouched. Directed pipeline.

## Verification (structural, no truth set)

`--command assemble` -> `scripts/Gate1AnchorGraphCollapseStats.pl` on
`Shasta2AnchorGraph.gfa`, before vs after:
- strand-strand (fw↔rc) edges -> drop toward ≈ #real IRs
- hairpin windows -> drop
- collapsed-window proxy (>3× median verts) -> drop
- connected components / contig count -> must NOT fragment the genome
New counters to emit: candidates, 1a-contracted, 1b-unioned, difficult by reason.
Success = strand-strand/hairpin/collapse fall sharply with no loss of legitimate
backbone connectivity.

## Scope boundaries

- Phase 1 (1a + 1b) only. No IR resolution, no interior-attachment bypass/split,
  no `detangleWindows`/`GTest` deletion (separate step, only after Phase 1 is
  verified to subsume their strand role).
- I cannot build/run here (no spoa/seqan/abpoa/rust). User builds, runs
  `--command assemble`, runs Gate-1 script before/after. Mandatory per V36x lesson.
