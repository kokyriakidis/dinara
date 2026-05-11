# Backbone-Disjoint Read-Window Theseus Prototype

## Goal
Replace repeated marker-pair MSAs with a simpler global read-window prototype:

- Select backbone reads greedily, longest first.
- A read can become a backbone only if it has not already been claimed by an earlier window.
- For each selected backbone, use all valid direct read-graph overlaps as window evidence.
- Claim only previously unclaimed direct overlaps so they cannot become future backbones.
- Allow already-claimed reads to appear again as borrowed/context evidence for a later window MSA.
- Keep enough metadata to later stitch windows using shared evidence reads and cross-window read graph edges.

## Core Algorithm
1. Build/access `alignmentData`.
2. Build/access the existing `readGraph`.
3. Sort physical reads by length descending.
4. Maintain read ownership/backbone eligibility:
   - `readOwner[readId] = invalid` means unclaimed.
   - `readOwner[readId] = windowId` means the physical read was claimed by that window and cannot become a future backbone.
5. For each read in length order:
   - If already claimed, skip it.
   - Assign a new stable `windowId`.
   - Use `OrientedReadId(readId, 0)` as the backbone seed.
   - Claim the seed physical read.
   - Scan the seed's strand-0 read graph adjacency list.
   - Skip read graph edges flagged `crossesStrands` or `hasInconsistentAlignment`.
   - Store every valid direct read-graph edge as window evidence, even if the other read was already claimed by an earlier window.
   - Do not deduplicate evidence by physical or oriented read id: multiple good chains/placements of the same read are meaningful in repetitive or inverted regions.
   - If the other physical read is still unclaimed, claim it for this window.
   - If the other physical read was already claimed by another window, record/count it as borrowed evidence and a cross-window connection.
6. Store a `WindowTask` for each selected seed.

This is a one-hop evidence window: seed plus direct overlaps only. Reads that overlap recruited reads but not the seed are not pulled in transitively. A read can appear as evidence in more than one window, but it can only be claimed once and therefore can only lose backbone eligibility once.

## Read Graph Neighbor Access
Do not build a separate CSR graph from `alignmentData`. The existing `readGraph` already stores the selected good alignments and their adjacency lists:

```cpp
for(const uint32_t edgeId: readGraph.connectivity[seedOid.getValue()]) {
    const ReadGraphEdge& edge = readGraph.edges[edgeId];
}
```

Neighbor rules:

- Use only the strand-0 seed vertex: `OrientedReadId(seedReadId, 0)`.
- Skip `edge.crossesStrands`.
- Skip `edge.hasInconsistentAlignment`.
- Use `edge.getOther(seedOid)` to get the neighbor oriented read for the MSA row orientation.
- Store `edge.alignmentId` so the window can later recover overlap details.

This avoids duplicating graph construction, stays consistent with `createReadGraphAllAlignments()`, and uses the same filtering flags as the rest of the read graph code.

## Window Identity
Each window gets a stable deterministic id:

```cpp
uint32_t windowId = uint32_t(windows.size());
```

Use `windowId` in:

- `readOwner[readId]`
- `WindowTask::windowId`
- future cross-window stitching edges

## WindowTask Shape
Define a local task struct in `src/AssemblerTheseusReadWindowMSA.cpp`:

```cpp
struct ReadWindowTask {
    uint32_t windowId = 0;
    ReadId backboneReadId;
    vector<OrientedReadId> orientedReads; // backbone plus one row per valid 1-hop edge
    vector<ReadId> claimedReads;          // backbone plus newly claimed reads only
    vector<uint32_t> alignmentIds;
};
```

`orientedReads[0]` should be `OrientedReadId(backboneReadId, 0)`.

For evidence reads, store the oriented read id returned by the read graph edge relative to the seed. This includes both newly claimed reads and already-claimed borrowed/context reads. The same physical or oriented read can appear multiple times if there are multiple valid read graph edges/alignment placements.

`claimedReads` is not the same as the eventual MSA rows. It is the ownership/backbone-blocking set.

## Sequence Policy
For this first prototype, keep it simple:

- The whole backbone read belongs to the window.
- Every evidence read in `orientedReads` contributes its whole oriented read sequence.
- Do not extract only overlap-coordinate segments yet.
- Add the backbone as the first Theseus sequence.
- Add evidence reads as ends-free initially.

This may create noisy overhangs, but it keeps the partitioning prototype easy to validate.

## Cross-Window Metadata
During planning, when scanning a seed neighbor:

- Always add the other oriented read as evidence for this window.
- If the other physical read is unclaimed, claim it into this window and block it from future backbone selection.
- If the other physical read is already owned by a different window, count or store a cross-window edge and treat it as borrowed evidence:

```cpp
struct ReadWindowEdge {
    uint32_t windowId0 = 0;
    uint32_t windowId1 = 0;
    uint32_t alignmentId = 0;
};
```

These edges are not used in the first verification prototype. They are kept for later stitching. Shared evidence reads are also stitching evidence.

## Parallel Execution
Planning is deterministic and single-threaded over the read graph. This keeps `readOwner` race-free and makes window ids stable.

MSA execution is parallel:

1. Create all `ReadWindowTask`s.
2. Use an atomic task index or fixed chunks across `threadCount` workers.
3. Each worker creates its own `theseus::TheseusMSA`.
4. Each worker builds sequences only for the task it is processing.
5. Each worker accumulates thread-local counters.
6. Merge counters after all workers finish.

No shared write state is needed during MSA execution.

For the current verification step, the Theseus MSA execution block is disabled. The prototype only creates windows and prints planning/integrity statistics.

## Output
Keep diagnostic output minimal:

- No per-site logs by default.
- No MSA file output by default.
- Print final summary:
  - read count
  - read graph edge count
  - scanned read graph edges
  - planned windows
  - claimed reads
  - unclaimed reads
  - singleton windows
  - max/average claimed reads per window
  - total borrowed evidence reads
  - max/average evidence reads per window
  - cross-window overlap count
  - backbone conflict edge count
  - owner mismatch count
  - planning seconds
  - total seconds

Expected verification invariants:

- `unclaimedReads == 0`
- `ownerMismatches == 0`
- `backboneConflictEdges == 0`
- `runMsa == 0` while the MSA block is disabled

## Integration Points
- Add a new method in `src/Assembler.hpp`:

```cpp
void computeTheseusReadWindowMSAPrototype(
    uint64_t threadCount);
```

- Implement it in `src/AssemblerTheseusReadWindowMSA.cpp`.
- Update the diagnostic path in `srcMain/main.cpp` to call the new prototype.
- Keep the current marker-pair and long-backbone prototypes available for comparison.

## Efficiency Notes
- No dense overlap matrix is needed.
- No duplicate CSR construction from `alignmentData` is needed.
- `readGraph.connectivity` gives direct neighbor iteration.
- Claim lookup is O(1) via `readOwner`.
- Planning cost is roughly proportional to the neighbor entries scanned for selected seed reads.
- Expensive work is the per-window Theseus MSA phase, which parallelizes cleanly.

## Future Extensions
- Use cross-window edges to build a window graph for stitching.
- Add optional max window size if whole-read MSAs become too large.
- Add overlap-coordinate segment extraction if full-read overhangs are too noisy.
- Add limited two-hop expansion only if one-hop windows fragment too much.
