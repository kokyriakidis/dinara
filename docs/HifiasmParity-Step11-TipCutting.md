# Step 11 Parity Check: Tip Cutting (hifiasm `asg_arc_cut_tips`)

Target upstream: hifiasm `0.25.0-r726` commit `ec9a8b222d149d25b7355e83765698640d59b189` in `/home/kokyriakidis/Downloads/hifiasm`.

This document verifies Dinara’s Step 11 (“cut short tips”) against hifiasm’s implementation and records any intentional deviations.

## Where Step 11 Happens

### hifiasm call site (pipeline order)

- `gfa_ut.cpp:3027` `ul_clean_gfa(...)`
- `gfa_ut.cpp:3059` first operation after (optional) `asg_arc_cut_weak`:
  - `asg_arc_cut_tips(sg, max_tip, &bu, is_ou, is_ou?rI:NULL, uopt->te);`

Notes:
- `ul_clean_gfa` calls `asg_arc_cut_tips` multiple times (inside and after the main cleaning loop), but Step 11 parity here refers to the *initial* call before the loop.

### Dinara call site (pipeline order)

- `src/AssemblerReadGraph6.cpp:330` Step 11:
  - `cutStringGraphTips(/*maxShortTipReads*/3);`

## hifiasm Call Graph (Step 11)

Entry:
- `gfa_ut.cpp:554` `asg_arc_cut_tips(asg_t *g, uint32_t max_ext, asg64_v *in, uint32_t is_ou, R_to_U *ru, telo_end_t *te)`

Helpers used by `asg_arc_cut_tips`:
- `gfa_ut.cpp:527` `asg_end(const asg_t *g, uint32_t v, uint64_t *lw, uint32_t *ou)`
  - This is the “unitig extension step” predicate used by `asg_arc_cut_tips`.
  - It is structurally equivalent to hifiasm’s classic `asg_is_utg_end` but with an optional `ou` output.
- `Overlaps.cpp:2647` `asg_is_utg_end(const asg_t *g, uint32_t v, uint64_t *lw)` (reference implementation)
- `Overlaps.cpp:2769` `asg_extend(const asg_t *g, uint32_t v, int max_ext, asg64_v *a)` (reference implementation)
  - `asg_arc_cut_tips` does *not* call `asg_extend` directly; it inlines a similar extension loop using `asg_end`.

Post-processing used by `asg_arc_cut_tips`:
- `gfa_ut.cpp:680` `asg_cleanup(g);` (only if it deleted at least one tip-unitig)

## hifiasm Behavior (Line-by-Line Summary)

Implementation: `gfa_ut.cpp:554..683`

Phase 1: gather candidate tips (short tip-unitigs)
- Iterate all oriented vertices `v` (`n_vtx = g->n_seq << 1`).
- Skip if read deleted: `g->seq[v>>1].del`.
- Skip if telomere-protected (if enabled): `te && te->hh[v>>1]`.
- Require “tip start” condition: no non-deleted arcs out of `v^1` (`asg_arc_a(g, v^1)` has zero non-deleted entries).
- Extend forward along a unique unitig path for up to `max_ext` steps using:
  - `asg_end(g, w^1, &lw, is_ou?&ou:NULL) == ASG_ET_MERGEABLE`
  - Update `w = (uint32_t)lw` each step.
- If extension stops early (`i < max_ext`) and the traversed path did not encounter a telomere, record candidate:
  - Key is `(((uint64_t)kv) << 32) | v`, where `kv` is the number of vertices in the path (used for sorting shortest-first).

Phase 2: apply deletions (shortest-first, re-checking)
- Sort candidates (`radix_sort_srt64` on the 64-bit key).
- For each candidate `v`, re-check:
  - read not deleted
  - still a tip start (`outdegree(v^1) == 0`)
- Re-traverse the same unitig path (up to `max_ext`).
- If still short (`i < max_ext`) and not telomere, delete the entire tip-unitig by marking reads deleted:
  - `asg_seq_del(g, ((uint32_t)pathVertex)>>1)`
- After processing all candidates, call `asg_cleanup(g)` if anything was deleted.

Optional Phase 3 (only when `ru && is_ou`)
- Additional “OU” specific deletion logic that consults `R_to_U` and re-traverses from tips while all traversed vertices remain non-unitig (`is_u == 0`) and mapped (`rr != -1`).

## Dinara Implementation (Step 11)

Entry:
- `src/AssemblerStringGraphClean.cpp:526` `Assembler::cutStringGraphTips(uint32_t maxShortTipReads)`
  - Calls static helper `stringGraphCutTips(*this, maxShortTipReads, vertexCount)`.

Helper call graph:
- `src/AssemblerStringGraphClean.cpp:337` `stringGraphCutTips(...)` (tip cutting)
- `src/AssemblerStringGraphClean.cpp:52` `stringGraphIsUtgEnd(...)` (parity helper for `asg_is_utg_end` / `asg_end`)
- `src/AssemblerStringGraphClean.cpp:20` `symmetrizeArcDeletion(...)` (enforces Dinara’s `arcId^1` invariant)
- `src/AssemblerStringGraphClean.cpp:115` `rebuildStringGraphAdjacency(...)` (Dinara equivalent of `asg_cleanup` + maintaining adjacency invariants)

Core behavior (non-OU / non-telomere parity mode):
- Phase 1: collect candidate tip-start vertices and compute their unitig length `kv` by extending up to `maxShortTipReads`.
- Sort candidates by `(kv, v)` ascending (equivalent to hifiasm’s key sort).
- Phase 2: re-check each candidate and delete the entire tip-unitig by:
  - marking `stringGraph.readDeleted[readId] = 1`
  - marking all outgoing arcs from both orientations (`v0`, `v1`) as deleted
  - calling `symmetrizeArcDeletion` so `arcId` and `arcId^1` deletions stay consistent
  - rebuilding adjacency to exclude deleted arcs

## Verified Parity and Known Gaps

Verified parity (for the common case in Dinara’s pipeline):
- “Tip start” definition: outdegree of `v^1` equals zero (matches hifiasm’s `asg_arc_a(g, v^1)` scan).
- “Short tip-unitig” criterion: extend up to `maxShortTipReads` and cut only if extension stops early (`i < max_ext`).
- Two-phase algorithm: candidates gathered first, sorted shortest-first, then re-validated before deletion (matches hifiasm’s structure).
- Cleanup effect: Dinara rebuilds adjacency lists and enforces RC-arc pairing; hifiasm calls `asg_cleanup(g)` (Dinara’s extra symmetrization is required by its `arcId^1` invariant).

Intentional gaps (not representable in current Dinara `StringGraphArc`):
- Telomere protection (`telo_end_t *te`, `te->hh[...]`) is not modeled in Dinara’s string graph, so Dinara cannot skip telomere-labeled nodes the way hifiasm can.
- “OU” mode fields and logic (`is_ou`, `R_to_U *ru`, arc `.ou`, `mm_ou` adjustments, and the `ru && is_ou` extra deletion phase) are not currently implemented.

If full parity for ONT UL cleaning is required (including `is_ou` + `ru` + `te`), Dinara will need:
- An `ou`-like per-arc metric in `StringGraphArc` (or a parallel structure),
- A representation for telomere-labeled reads,
- A translation of `R_to_U` semantics into Dinara’s data model and pipeline.

## Related ONT Step In hifiasm: `asg_arc_cut_weak`

In `ul_clean_gfa`, hifiasm performs an ONT-only weak-arc pruning step immediately before the initial tip cut:
- `gfa_ut.cpp:3057`: `if(asm_opt.is_ont) asg_arc_cut_weak(...)`

Dinara implements a best-effort approximation of that ONT pre-tip-cut stage as:
- `src/AssemblerReadGraph6.cpp`: Step 10b calls `cutStringGraphWeakArcsOntHifiasm(3, 0.975, 16)`
- `src/AssemblerStringGraphClean.cpp`: `Assembler::cutStringGraphWeakArcsOntHifiasm`

Current limitations (parity gaps):
- hifiasm’s `asg_arc_cut_weak` uses per-overlap `strong` (`ma_hit_t.ml`) and optional `ou`/`R_to_U`/trio logic.
- Dinara does not currently store the `ml`/`ou` fields on `StringGraphArc`, so the implementation uses overlap-length heuristics and checks for the existence of overlaps between candidate targets via the `alignmentTable`.
