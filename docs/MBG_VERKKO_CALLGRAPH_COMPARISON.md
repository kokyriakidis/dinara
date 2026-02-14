# DirectedAnchorGraph vs MBG/Verkko: Call Graph + Step Check

## Scope checked
- Local code path: `srcMain/main.cpp:1225` -> `buildFromAnchors()` -> `unitigifyAll()`.
- MBG version: `2e1d44b7a001e35d2c40c31520d11f7c621a7506`.
- Verkko repo HEAD checked: `e479f0e` (2026-02-06), with submodule pointer:
  `src/MBG -> 2e1d44b7a001e35d2c40c31520d11f7c621a7506`.

## Implementation Update
- The local DAG pipeline has now been extended to include:
  - in-round MBG-style cleaning passes during `resolveRound` (tips 3/10, tips 2/5, crosslinks 1/5 and 2/10),
  - a final low-support resolve pass (`minEdgeSupport=1`),
  - a non-stub `cleanComponentsByCopynumber` implementation.
- Main wiring now calls the upgraded `resolveRound(..., doCleaning, doGuesswork)`.

## 1) Local call graph (your current snippet)

```text
main.cpp
  -> DirectedAnchorGraph::buildFromAnchors(anchors, threadCount)
     -> setupLoadBalancing(...)
     -> runThreads(buildFromAnchorsThreadFunction, threadCount)
        -> buildFromAnchorsThreadFunction(threadId)
           -> build path from anchors.journeys[orid]
           -> collect threadEdgePairs (path[i] -> path[i+1])
     -> merge thread paths
     -> deduplicate + insert edges (+RC edges)
     -> pathsCrossing.rebuild(paths, pathRemoved)
  -> DirectedAnchorGraph::unitigifyAll()
     -> extendForward(fwd), extendForward(rev) for each candidate segment
     -> batch-create merged nodes and substitution map
     -> one global path-rewrite pass
     -> removeNode(old chain segments)
```

References:
- `src/mode3-DirectedAnchorGraph.cpp:300`
- `src/mode3-DirectedAnchorGraph.cpp:434`
- `src/mode3-DirectedAnchorGraphResolution.cpp:633`
- `src/mode3-DirectedAnchorGraphResolution.cpp:423`
- `src/mode3-DirectedAnchorGraphResolution.cpp:458`

## 2) MBG call graph (as used by Verkko)

```text
MBG main()
  -> runMBG(...)
     -> loadReadsAsHashesMultithread(...)
     -> getUnitigGraph(...)                  # builds + unitigifies k-mer graph
     -> (optional) filterUnitigsByCoverage + getUnitigs(...)
     -> getReadPaths(...)
     -> if maxResolveLength > 0:
          resolveUnitigs(...)
            -> getUnitigs(initial, ...)
            -> unitigifyAll(resolvableGraph, readPaths)
               -> for each node: unitigifyOne(...)
                  -> getUnitigPath(...)
                  -> replaceUnitigPath(...)
                  -> replacePathNodes(...)
            -> resolveRound(...)
            -> resolveRound(...)
```

References:
- `/tmp/MBG/src/main.cpp:235`
- `/tmp/MBG/src/MBG.cpp:1546`
- `/tmp/MBG/src/MBG.cpp:1583`
- `/tmp/MBG/src/MBG.cpp:1602`
- `/tmp/MBG/src/UnitigResolver.cpp:4271`
- `/tmp/MBG/src/UnitigResolver.cpp:1525`
- `/tmp/MBG/src/UnitigResolver.cpp:1514`

## 3) Verkko pipeline call graph around this stage

```text
Snakefile rule buildGraph
  -> MBG binary (same MBG commit above)

Snakefile rule processGraph
  -> scripts/unitigify.py                   # extra unitigify on transformed graph

Snakefile rule processONT
  -> scripts/resolve_triplets_kmerify.py
     -> unitigify_all()
        -> unitigify_one()
           -> extend_forward()
           -> replace_unitig() + path rewrite
     -> for coverage in resolve_steps: resolve(...)
  -> scripts/unitigify.py
```

References:
- `/tmp/verkko/src/Snakefiles/1-buildGraph.sm:88`
- `/tmp/verkko/src/Snakefiles/2-processGraph.sm:106`
- `/tmp/verkko/src/Snakefiles/4-processONT.sm:323`
- `/tmp/verkko/src/scripts/resolve_triplets_kmerify.py:987`
- `/tmp/verkko/src/scripts/unitigify.py:11`

## 4) Step-by-step equivalence check

1. `buildFromAnchors()` vs MBG graph construction:
   - Local: node = anchor pair, edges/paths from anchor journeys.
   - MBG: node = selected minimizer k-mer, edges from sequence overlaps, then unitig graph.
   - Result: **not exact**.

2. `unitigifyAll()` placement:
   - Local: explicit call right after graph build.
   - MBG: initial unitigification already inside `getUnitigGraph`; later resolver does another `unitigifyAll`.
   - Result: **not exact sequencing**.

3. Unitigification execution model:
   - Local: batched all-chains, one global path rewrite pass.
   - MBG resolver: per-node `unitigifyOne` loop, repeated path rewrites.
   - Verkko `resolve_triplets_kmerify.py`: also per-node style.
   - Result: **not exact implementation**.

4. Chain-growth rule (`extendForward` / degree-1 walking):
   - Local and MBG/Verkko use the same core idea (linear extension using 1-in/1-out constraints and loop guards).
   - Result: **conceptually similar**.

5. Path handling:
   - Local: `pathsCrossing` rebuilt once after build, then rewritten in a single substitution pass in `unitigifyAll`.
   - MBG resolver: `replacePathNodes` with additional clip/trim logic and incremental updates.
   - Result: **not exact**.

6. Current runtime pipeline difference:
   - In `srcMain/main.cpp:1247` there is an early `return;` immediately after initial `unitigifyAll`, so your resolve rounds are currently disabled.
   - MBG/Verkko continue with additional resolve rounds.

## Bottom line
- Verkko is using the exact MBG commit you linked.
- Your `buildFromAnchors -> unitigifyAll` is **not exactly the same** as MBG/Verkko function-by-function.
- It is **algorithmically inspired/similar** in some rules (degree-1 unitig growth, RC-aware edges), but differs in data model, call order, and path-rewrite strategy.
