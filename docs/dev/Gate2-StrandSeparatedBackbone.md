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

## Phase 1b — Resolve junctions (corroborated, strand-guarded, shortest-first)

Adopts Verkko's control flow: NO single read seeds the backbone (a long chimeric
read must never be a seed). Support is accumulated from ALL read paths
(order-independent), then nodes are resolved **shortest-first**, length-stepped,
re-contracting safe unitigs after each round and iterating at increasing length.
This resolves the easiest (short, well-spanned) ambiguities first and simplifies
larger structure before attempting it.

**Triplet** = `(prevWindow, window, nextWindow)` from a read path crossing the
window. A triplet is **solid** iff count >= N AND span product >= T. A junction is
**resolvable** iff EVERY significant incident edge is covered by a solid triplet
(all-or-nothing, per Verkko); otherwise the whole node is deferred.

```
# 1. Accumulate triplet support from ALL read paths (order-independent).
for each read R:
    project journey -> window chain w0..wk
    for each interior wi (with prev wi-1, next wi+1):
        triplet[(wi-1, wi, wi+1)].count += 1
        triplet[(wi-1, wi, wi+1)].spanProduct += transition supportingSpanProduct

# 2. Resolve shortest-first, length-stepped, iterate.
heap = windows not in safe unitigs, keyed by window length (SHORTEST first)
while heap not empty and currentLength <= maxResolveLength:
    for each window W at this length:
        if W is hairpin / fw-rc-same-read: defer (Phase 2); continue
        incident = significant incident edges of W
        if every edge in incident is covered by a solid triplet
           AND all implied unions are strand-safe
           (find(side) != find(complement(otherSide))):
              for each solid triplet through W:
                  union(prev, W); union(W, next); mirror-union complements
        else:
              defer W (difficult, reason = uncorroborated/strand/partial)
    re-contract safe (1,1) unitigs; re-push affected windows
    advance currentLength

A chimeric read contributes count 1 to its false triplet; no other read
corroborates -> never solid -> never unioned. The chimera screens itself out, and
because no read seeds the backbone, a long chimera cannot anchor a false join.
```

### Safe-unitig immunity (per Verkko)

Windows contracted into safe (1,1) unitigs in 1a are **protected**: they are not
subject to the corroboration gate or deferral in 1b. Only non-safe (junction)
windows are resolved or deferred.

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

### Bidirected vs explicit-mirror: where the parity guard fires

Verkko uses a bidirected graph (one node per locus, strand is an edge property),
so a strand-strand contact is a hairpin `>node -> <node` on a single node id, and
its resolution is to **split** that node into fw/bw copies. We use a directed
graph with **explicit RC mirror windows** — strand is a node identity, two nodes
(`w`, `complement(w)`) per locus. The mirror split Verkko performs at resolution
time is therefore **already done at construction** for us; Verkko's "split" maps
to a **no-op** here, and its read-support gate ports unchanged (it is the real
engine in both systems).

Consequence — **the parity guard must fire on UNITE/contraction, never on edge
admission.** Because strands are separate nodes, a *real inverted-repeat
fold-back is an edge from a fw-window to an rc-window* (`A -> A'`), and for any
fold-back `complement(A') == A`, so `find(A) == find(complement(A'))` fires
trivially. Applying the parity test at edge admission would delete every genuine
IR. The correct rule:

- **Admit** a corroborated fw->rc edge — it is a legitimate edge between two
  distinct nodes (component-of-`A` and component-of-`A'`).
- **Forbid only the UNITE** that would merge a component with its own mirror
  (`union` whose two sides satisfy `find(X) == find(complement(Y))`). Fusing a
  component into its mirror is the marker-graph collapse signature; edges
  *between* mirror components are fine, *fusing* them is the violation.

So the discriminator between "real IR" and "collapse artifact" is the same as
Verkko's — read support — and the only model-dependent difference is the verdict:
Verkko splits; we simply decline to fuse (the split already exists). Both 1a
mirror-union and 1b resolution unions are subject to this guard; edge admission
is not.

### Decision: strand work stays in directed space; fold last

We keep all Phase 1 strand work in **directed space** (fw and rc as separate
nodes) and let the existing `Shasta2AssemblyGraph::toBidirected` perform the
natural fold **at the end of the pipeline**, after read support has decided what
is a real IR and what is a collapse artifact. Rationale:

- The graph model does NOT discriminate IR vs collapse — both are strand-strand
  contacts. **Read support** discriminates, and counting it requires fw and rc to
  be separate, countable nodes. Directed space preserves that; folding first
  destroys it.
- Folding (bidirected) merges `A` and `A'` into one node, collapsing a real IR
  and a collapse artifact into the *same* hairpin *before* the read-support gate
  runs. Verkko's `resolve_hairpins` exists precisely to un-fold (split back) what
  the bidirected model folded too eagerly — evidence that the fold belongs at the
  end, not the start.
- A bidirected migration is a large, risky rewrite of working machinery
  (`windowPairTransitions`, `readWindows`, mirror window IDs) that does NOT solve
  the core problem. The legitimate fold we want already exists at the right place
  (`toBidirected`, end of pipeline).

Ordering: **discriminate in directed space (Phase 1), then fold at
`toBidirected`** — never fold-then-guess.

## Data available (verified)

- `windowPairTransitions: map<pair<uint32_t,uint32_t>, vector<ReadTransition>>`,
  window IDs encode strand.
- `ReadTransition.supportingSpanProduct`, `supportingSpanA/B`, `oidValue`.
- fw-fw/rc-rc/fw-rc/rc-fw classification + RC-mirror canonicalization assert.
- `dset64` / `Shasta2DisjointSets`.

## Where it plugs in

**Separate post-pass, NOT constructor surgery.** The working constructor
(`src/Shasta2AnchorGraph.cpp`) builds the full directed graph (including the
inter-window edges at ~626–672) unchanged. Phase 1 runs as a new method invoked
*after* construction, consuming `windowPairTransitions` and the built edge set,
and rewriting edges in place: it contracts (1a), resolves (1b), and then keeps
only edges backed by a union; un-unioned edges are disabled via the existing
`disableEdge` (RC-mirror aware). This isolates the new logic from the proven
build path — if the post-pass is a no-op, the graph is identical to today's.
Directed pipeline; intra-window backbone untouched.

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

## Cross-check against Verkko (resolve_triplets_kmerify.py)

Verkko's ONT resolution stage independently arrives at the same structure,
validating this design and suggesting refinements:

- `add_safe_unitig` / `get_safe_unitigs_and_edges`: walks one-in/one-out linear
  runs and marks them **safe** — this is our 1a. Refinement: safe unitigs are
  **protected** from the low-coverage pruning that governs junctions. Adopt:
  the 1a backbone should be immune to the 1b corroboration gate.
- `get_valid_triplets`: collects `(prev, node, next)` triplets from crossing read
  paths; a triplet is **solid** iff path count >= `min_edge_support` (our count N).
  A node is resolvable iff **every** significant incident edge is covered by a
  solid triplet, else it is NOT resolved (-> difficult). Refinement: adopt
  **all-or-nothing resolvability** — do not partially resolve a junction.
  Difference: Verkko uses path **count** only; we additionally require span
  product T (more signal). Keep T as our addition.
- `resolve_hairpins` / `is_hairpin`: Verkko's term for a strand-strand contact is
  a **hairpin** = a node whose forward end connects to its own reverse end
  (`>node -> <node`). Handling (Phase 2 should mirror this exactly):
  - **Detect** three categories: clean hairpin (`>node` has exactly one edge, to
    `<node`); `unresolvable_hairpin` (fold-back PLUS other branches,
    `len(edges)>=2 and <node in edges`); double hairpin (`revnode(node)` also a
    hairpin, folds at both ends).
  - **Resolve only if read-supported, all-or-nothing**: gather read paths that
    physically execute the fold-back (`...prev, node, revnode(node), next...`),
    count per `canon(prev,next)` resolution; a resolution is solid iff
    count >= min_edge_support; resolve ONLY if every incident edge is covered by a
    solid resolution, else leave in place.
  - **Resolution = NODE SPLIT, not delete**: split the node into separate fw/bw
    copies (`node+hairpinN+fw/bw`), route each read's supported fold-back through
    its own copy, remove the original. This physically separates the two strands
    into two nodes connected only via the supported fold path. Verkko does NOT
    delete strand-strand contacts — it splits the read-supported ones and leaves
    the rest. Directly supports our inverted-repeat requirement: a real IR
    fold-back is read-supported -> split into clean separate strands, never
    destroyed.
  - **Quarantine rules (adopt)**: skip double hairpins ("too hard"); skip
    disconnected hairpins (`len(edges[revnode])==0`, "spurious"); skip
    unresolvable hairpins; and do NOT resolve normal nodes ADJACENT to an
    unresolved hairpin (Verkko removes hairpin-bordering nodes from
    maybe_resolvable).
  Phase 1 defers all fw/rc-same-read & self-fold windows to Phase 2, which applies
  the split + quarantine rules above.
- `resolve` driver: heap of nodes by unitig length, **shortest first**,
  length-stepped with a `max_resolve_length` cap, re-unitigify and re-push after
  each resolution; outer loop over increasing `resolve_steps`. ADOPTED: 1b uses
  shortest-first, length-stepped, iterate (replacing the earlier
  longest-read-first seeding, which is unnecessary and risks a long chimera
  seeding a false join). Support is accumulated from all read paths
  order-independently; no single read seeds the backbone.

## Scope boundaries

- Phase 1 (1a + 1b) only. No IR resolution, no interior-attachment bypass/split,
  no `detangleWindows`/`GTest` deletion (separate step, only after Phase 1 is
  verified to subsume their strand role).
- I cannot build/run here (no spoa/seqan/abpoa/rust). User builds, runs
  `--command assemble`, runs Gate-1 script before/after. Mandatory per V36x lesson.
