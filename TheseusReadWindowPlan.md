# Anchor-Interval Theseus Window Prototype

## Goal
Replace repeated marker-pair MSAs with windows defined directly on anchor journeys:

- Build **Shasta2** anchors from marker-graph vertices, then **Shasta2Journeys** (same journey algorithm as upstream [shasta2](https://github.com/paoloshasta/shasta2); Dinara uses a marker-graph-based anchor constructor, not upstream’s default MarkerKmers path—see **Anchor source** below).
- Run planning in `computeTheseusReadWindowMSAPrototype`; **creation of anchors and journeys stays outside** that call (typically `srcMain/main.cpp`) so the same `assembler.shasta2Anchors` / `assembler.shasta2Journeys` can feed later Shasta2 assembly steps.
- Use only strand-0 as a canonical backbone orientation.
- Select backbone intervals from a dynamic max-heap ordered by current unclaimed base span.
- Claim anchors, not whole reads.
- Expand each seed interval through reads present in the seed anchors, claiming contiguous unclaimed intervals on those reads.
- Later, run one MSA per anchor-interval window and detect het sites with direct mapping back to read positions and anchor intervals.

## Core Algorithm
1. Build **Shasta2Anchors** from the marker graph (`minAnchorCoverage` / `maxAnchorCoverage` on vertices), then build **Shasta2Journeys** (fills `positionInJourney` on anchor marker infos). Optionally call `computeTheseusReadWindowMSAPrototype(anchors, journeys, threadCount)` or use the one-arg overload that constructs Shasta2 objects locally.
2. Sort physical reads by length descending for deterministic initialization/tie-breaking.
3. Maintain anchor ownership:
   - `anchorOwner[anchorId] = invalid` means unclaimed.
   - `anchorOwner[anchorId] = windowId` means that anchor belongs to a window.
4. Initialize a max-heap with one full-journey candidate per strand-0 read:
   - Use `OrientedReadId(readId, 0)` as the canonical backbone.
   - Candidate priority is base span descending, then anchor count, then read length.
5. Repeatedly pop the best candidate:
   - Discard old-generation candidates.
   - Revalidate that every anchor in the candidate interval is still unclaimed.
   - If stale, rescan that read's journey, push its current maximal unclaimed intervals, and continue.
   - If valid and long enough, create a window from that interval.
6. To create a window:
   - Claim the seed backbone anchor interval.
   - Inspect all anchor members in the seed interval.
   - For every oriented read present in those seed anchors, find the leftmost and rightmost positions where it touches the seed interval.
   - In that read's journey, scan between those bounds and claim every maximal contiguous unclaimed anchor run.
   - Store each claimed run as an `AnchorWindowReadInterval`.
7. Store an `AnchorWindowTask` for each selected seed interval.

This makes the partition unit an anchor-journey interval rather than a whole physical read. A physical read can participate in multiple windows through disjoint journey intervals. Backbones can overlap globally, but they cannot seed on already claimed anchor intervals. Interval splitting is exact by journey index, while candidate priority uses the estimated base span between interval endpoints.

## Candidate Priority
The heap is a `std::priority_queue`, so it is not fully sorted. It only guarantees that `top()` is the currently best known candidate. Candidates are ordered by:

1. Larger estimated base span.
2. Larger anchor count.
3. Longer physical read length.
4. Smaller physical read id.

Base span is a score in bases, not in anchor count. To compute it, the planner:

1. Takes the left and right endpoint **Shasta2** anchors of the candidate interval.
2. Uses `shasta2Anchors->getOrdinal(anchorId, orientedReadId)` for each endpoint (Shasta2 anchors are one marker-graph vertex per anchor; there is no separate `ordinalOffset` step like mode3 two-marker anchors).
3. Looks up the corresponding marker base positions in `markers`.
4. Uses the difference between those base positions as the candidate span (with the same end-marker position semantics as before on the backbone).

The interval is still represented by journey indexes `[begin, end)`. The base span is only a priority score for choosing which interval to process next. All claiming and splitting still happens by anchor ownership in the journey.

## Lazy Revalidation
Window creation changes the ownership state of many anchors. Updating every candidate interval for every affected read immediately would be complicated and unnecessary, so the planner uses lazy revalidation:

1. A candidate is popped from the heap.
2. Its generation is checked. If a newer rescan of the same read already happened, this old candidate is discarded.
3. If the generation is current, the candidate interval is scanned in `anchorOwner`.
4. If every anchor is still unclaimed, the candidate is valid and becomes a window.
5. If any anchor is already claimed, the candidate is stale:
   - increment that read's generation,
   - rescan the read's whole journey,
   - push all current maximal unclaimed intervals for that read,
   - discard the stale candidate.

`staleCandidateIntervals` counts candidates that reached the top of the heap but were partially claimed. `discardedOldGenerationCandidates` counts obsolete candidates skipped cheaply because a newer rescan already replaced that read's interval candidates.

This keeps memory and update cost practical while preserving the correct behavior: every accepted window is created only from a currently valid unclaimed interval.

## Window Expansion
After a candidate interval is accepted, the planner creates an anchor window:

1. Claim all anchors in the backbone interval.
2. Visit all anchor members in that backbone interval.
3. For each touched oriented read, record the leftmost and rightmost `positionInJourney` where it touches the seed interval.
4. For that touched read, scan the bounded range between those positions.
5. Claim all currently unclaimed anchors in that bounded range.
6. Store the claimed anchors as one or more maximal contiguous `AnchorWindowReadInterval`s.

In other words, for a touched read we claim all unclaimed anchors between the left/right bounds. If claimed anchors already exist inside the bounded range, the stored result is split into separate contiguous intervals. This is how reads are progressively split across windows without maintaining an eager mutable interval list per read.

## Invariants
The planner is intended to maintain these invariants:

- Each `anchorId` is claimed by at most one window.
- Each `AnchorWindowReadInterval` is contiguous in its oriented read journey.
- A physical read can appear in multiple windows, but only through disjoint claimed anchor intervals.
- Backbone reads are always strand 0, while non-backbone intervals can be either strand depending on anchor membership.
- The queue can contain stale candidates, but every candidate is validated before it can create a window.

## Anchor source
- **Dinara Shasta2 path:** `Shasta2Anchors` selects marker-graph vertices by coverage and stores `Shasta2AnchorMarkerInfo` per oriented read at that vertex. `Shasta2Journeys` is then built from that table (ported from upstream `Journeys.cpp`, with mmap names prefixed, e.g. `Shasta2Journeys`, under `shasta2MappedMemoryOwner()`).
- **Upstream Shasta2 default:** anchors are normally created from **MarkerKmers** plus repeat / distinct-subkmer filters. The journey **algorithm** matches; the **set of anchors** will match upstream only if you use the same definition (e.g. external anchors export/import).

## Anchor Access

```cpp
// Typical wiring from main (after marker graph is ready):
const MappedMemoryOwner shasta2Owner = assembler.shasta2MappedMemoryOwner();
assembler.shasta2Anchors = make_shared<Shasta2Anchors>(
    shasta2Owner, assembler.getReads(), assembler.assemblerInfo->k,
    *assembler.markers, assembler.markerGraph, threadCount,
    minPrimaryCoverage, maxPrimaryCoverage);
assembler.shasta2Journeys = make_shared<Shasta2Journeys>(
    2 * assembler.getReads().readCount(),
    assembler.shasta2Anchors, threadCount, shasta2Owner);

assembler.computeTheseusReadWindowMSAPrototype(
    assembler.shasta2Anchors, assembler.shasta2Journeys, threadCount);
```

Important data:

- `(*shasta2Journeys)[orientedReadId]`: ordered `Shasta2AnchorId`s on that oriented read (`span<const Shasta2AnchorId>`).
- `(*shasta2Anchors)[anchorId]`: anchor as `span<const Shasta2AnchorMarkerInfo>`.
- `Shasta2AnchorMarkerInfo::orientedReadId`, `::ordinal`, `::positionInJourney` (filled by journeys).

## Window Identity
Each window gets a stable deterministic id:

```cpp
uint32_t windowId = uint32_t(anchorWindows.size());
```

Use `windowId` in:

- `anchorOwner[anchorId]`
- `AnchorWindowTask::windowId`
- future window graph/stitching metadata

## Window Shape
Define a local task struct in `src/AssemblerTheseusReadWindowMSA.cpp`:

```cpp
struct AnchorWindowReadInterval {
    OrientedReadId orientedReadId;
    uint32_t begin = 0; // inclusive journey position
    uint32_t end = 0;   // exclusive journey position
    uint32_t sharedBackboneAnchors = 0;
};

struct AnchorWindowTask {
    uint32_t windowId = 0;
    OrientedReadId backboneOrientedReadId;
    uint32_t backboneBegin = 0;
    uint32_t backboneEnd = 0;
    uint64_t claimedAnchorCount = 0;
    vector<AnchorWindowReadInterval> readIntervals;
};
```

Each `AnchorWindowReadInterval` is a contiguous interval on one oriented read journey. The first interval is the backbone seed interval. Other intervals are bounded by the read's participation in the seed anchors and split into maximal unclaimed runs before claiming.

## Sequence Policy
For MSA, the intended sequence policy is:

- Extract the sequence segment corresponding to each `AnchorWindowReadInterval`.
- Use the backbone interval as the MSA backbone row.
- Add other read intervals as ends-free rows.
- Detect het sites in the MSA and map them back to read positions through the interval metadata.

The current implementation is planning-only and reports window/interval statistics before enabling MSA.

## Parallel Execution
Planning is deterministic and single-threaded over anchor journeys. This keeps `anchorOwner` race-free and makes window ids stable.

MSA execution is parallel:

1. Create all `AnchorWindowTask`s.
2. Use an atomic task index or fixed chunks across `threadCount` workers.
3. Each worker creates its own `theseus::TheseusMSA`.
4. Each worker extracts interval sequences only for the task it is processing.
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
  - anchor count
  - planned windows
  - claimed anchors
  - unclaimed anchors
  - backbone intervals
  - candidate intervals pushed/popped
  - stale candidate intervals
  - old-generation candidates discarded
  - read intervals
  - split read intervals
  - touched reads
  - max/average read intervals per window
  - histograms: `readIntervalsPerWindowHistogram`, `claimedAnchorsPerWindowHistogram`, `backboneAnchorSpanHistogram` (see **Interpreting `readIntervalsPerWindowHistogram`**)
  - backbone vs expansion split: `backboneClaimedAnchors`, `nonBackboneClaimedAnchors`
  - skipped no-journey reads
  - skipped short runs
  - anchor-window planning seconds (`anchorWindowSeconds`)
  - total seconds (`totalSeconds` — see **Timing** below)

Expected verification checks:

- `claimedAnchors + unclaimedAnchors == anchors`
- Every stored `AnchorWindowReadInterval` is contiguous in its oriented read journey.
- The MSA block remains disabled until interval extraction and het-site mapping are validated.

## Timing
When **Shasta2Anchors** and **Shasta2Journeys** are built in `main` before the prototype, the log line **`totalSeconds`** inside `[TheseusReadWindowMSA]` measures **only** the anchor-window planner (heap + claiming + stats)—typically a few hundred milliseconds for ~1M anchors. The **several seconds** of anchor + journey construction appear **earlier** in the log (Shasta2 constructor / `Journeys creation`). For an end-to-end number, time from before `make_shared<Shasta2Anchors>` through the prototype summary (or add an explicit timestamp pair in `main`).

Older runs that built **mode3** anchors **inside** the prototype after “Prototype begins” included that work in the same `totalSeconds`, which made planning look much slower than it is.

## Interpreting `readIntervalsPerWindowHistogram`
The histogram buckets **`readIntervals.size()`** per window—the number of stored **`AnchorWindowReadInterval`** rows—not the backbone anchor count.

- Every window starts with **at least one** row: the **backbone** seed `[backboneBegin, backboneEnd)`.
- Additional rows appear only when expansion finds **non-empty runs of still-unclaimed** anchors on **other** reads between their min/max journey touch with the seed.

A large count in bucket **`1`** (e.g. `1:4955`) means: for those windows, **no other read** contributed a recorded interval—usually because **earlier** (larger-span) windows already claimed every overlapping anchor on neighbors, so only the backbone strip remains unclaimed for this window. That is expected under global anchor ownership and largest-first scheduling.

Use **`backboneAnchorSpanHistogram`** for how many **anchor positions** the backbone interval spans; a window can have `readIntervals == 1` and still have a long backbone (many anchors).

## Current Observations
Representative run (~25k reads, ~1.06M Shasta2 anchors):

- Top reads by length show the intended contrast between **journey anchor count** and **`journeyBaseSpan`** (e.g. long read with moderate anchor count vs shorter read with many dense anchors)—base span remains the right primary heap key.
- Coverage is nearly complete: on the order of **~1.06M** anchors claimed, **thousands** unclaimed, **~5.7k** windows (exact numbers depend on marker graph / filters).
- **`backboneClaimedAnchors`** dominates **`nonBackboneClaimedAnchors`**: expansion adds material support but backbone intervals drive most claims.
- Planning alone stays **sub-second** when timed inside the prototype; full Shasta2 anchor + journey build is separate.
- **`readIntervalsPerWindowHistogram`**: many windows have only **one** read-interval row (backbone-only expansion); tails extend to hundreds of intervals for a few windows (manageable for future MSA batching).
- Heap churn remains visible (`staleCandidateIntervals`, `discardedOldGenerationCandidates`) but is an intentional tradeoff for lazy revalidation.

## Integration Points
- `src/Assembler.hpp`:
  - `void computeTheseusReadWindowMSAPrototype(uint64_t threadCount);` — builds local Shasta2 anchors + journeys, then calls the overload below.
  - `void computeTheseusReadWindowMSAPrototype(shared_ptr<Shasta2Anchors>, shared_ptr<Shasta2Journeys>, uint64_t threadCount);`
- `src/AssemblerTheseusReadWindowMSA.cpp`: planner implementation.
- `srcMain/main.cpp`: diagnostic path creates `assembler.shasta2Anchors` / `shasta2Journeys`, then calls the two-argument overload (nested scope avoids duplicate coverage constant names with the full assembly path).
- Marker-pair and target-backbone Theseus prototypes remain available for comparison.

## Efficiency Notes
- No dense overlap matrix is needed.
- No duplicate CSR construction from `alignmentData` is needed.
- Claim lookup is O(1) via `anchorOwner[anchorId]`.
- Per-window temporary state is stored in reusable vectors indexed by `OrientedReadId`.
- Planning cost is roughly proportional to the anchor members scanned for selected backbone intervals.
- Expensive work is the per-window Theseus MSA phase, which parallelizes cleanly.

## Future Extensions
- Extract sequence intervals from anchor bounds.
- Run Theseus MSA per `AnchorWindowTask`.
- Detect het sites and store per-oriented-read positions.
- Use het sites to filter/split anchors inside a window.
- Build a window graph for stitching/phasing from shared reads and neighboring anchor intervals.
