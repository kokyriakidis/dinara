# Anchor-Interval Theseus Window Prototype

## Goal
Replace repeated marker-pair MSAs with windows defined directly on anchor journeys:

- Build mode3 anchors/journeys from marker-graph vertices.
- Use only strand-0 as a canonical backbone orientation.
- Select backbone intervals from a dynamic max-heap ordered by current unclaimed base span.
- Claim anchors, not whole reads.
- Expand each seed interval through reads present in the seed anchors, claiming contiguous unclaimed intervals on those reads.
- Later, run one MSA per anchor-interval window and detect het sites with direct mapping back to read positions and anchor intervals.

## Core Algorithm
1. Build mode3 anchors from marker-graph vertices and compute journeys.
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
   - Store each claimed run as a `WindowReadInterval`.
7. Store a `WindowAnchorInterval`/`WindowTask` for each selected seed interval.

This makes the partition unit an anchor-journey interval rather than a whole physical read. A physical read can participate in multiple windows through disjoint journey intervals. Backbones can overlap globally, but they cannot seed on already claimed anchor intervals. Interval splitting is exact by journey index, while candidate priority uses the estimated base span between interval endpoints.

## Candidate Priority
The heap is a `std::priority_queue`, so it is not fully sorted. It only guarantees that `top()` is the currently best known candidate. Candidates are ordered by:

1. Larger estimated base span.
2. Larger anchor count.
3. Longer physical read length.
4. Smaller physical read id.

Base span is a score in bases, not in anchor count. To compute it, the planner:

1. Takes the left and right endpoint anchors of the candidate interval.
2. Uses `anchors->getFirstOrdinal(anchorId, orientedReadId)` to find the marker ordinal of each endpoint anchor on the candidate backbone read.
3. Looks up the corresponding marker base positions in `markers`.
4. Uses the difference between those base positions as the candidate span.

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

## Anchor Access
Use the same mode3 marker-graph-vertex anchors/journeys as the existing Theseus marker-graph prototype:

```cpp
auto anchors = make_shared<mode3::Anchors>(
    MappedMemoryOwner(*this),
    getReads(),
    assemblerInfo->k,
    *markers,
    markerGraph,
    2,
    numeric_limits<uint64_t>::max(),
    threadCount,
    true);
anchors->computeJourneys(threadCount);
```

Important data:

- `anchors->journeys[orientedReadId.getValue()]`: ordered anchor ids on an oriented read.
- `(*anchors)[anchorId]`: anchor members as `AnchorMarkerInterval`s.
- `AnchorMarkerInterval::orientedReadId`: oriented read containing the anchor.
- `AnchorMarkerInterval::positionInJourney`: position of this anchor in that oriented read's journey.

## Window Identity
Each window gets a stable deterministic id:

```cpp
uint32_t windowId = uint32_t(windows.size());
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
  - max/average claimed anchors per window
  - skipped no-journey reads
  - skipped short runs
  - anchor-window planning seconds
  - total seconds

Expected verification checks:

- `claimedAnchors + unclaimedAnchors == anchors`
- Every stored `AnchorWindowReadInterval` is contiguous in its oriented read journey.
- The MSA block remains disabled until interval extraction and het-site mapping are validated.

## Current Observations
On the current test dataset:

- The heap/base-span priority is behaving as intended. A read can have fewer anchors but a much larger sequence span than another read with many dense anchors. For example, one run showed:
  - `readId=3729`: `journeyAnchors=493`, `journeyBaseSpan=298921`.
  - `readId=16296`: `journeyAnchors=1368`, `journeyBaseSpan=47192`.
  This supports using base span as the main priority score and anchor count only as a tie-breaker.
- Anchor coverage is essentially complete:
  - `claimedAnchors=782717`
  - `unclaimedAnchors=85`
  - `anchors=782802`
- Most claimed anchors come from accepted backbone intervals, with bounded expansion adding supporting context:
  - `backboneClaimedAnchors=621898`
  - `nonBackboneClaimedAnchors=160819`
  This suggests windows are driven primarily by coherent backbone intervals rather than by uncontrolled expansion.
- The heap/lazy-revalidation churn is visible but acceptable:
  - many candidates can become stale after windows claim anchors,
  - old-generation candidates are discarded cheaply,
  - planning still takes about a quarter second on this dataset.
- Window sizes look manageable for an MSA dry run:
  - most windows have at most 128 read intervals,
  - only a small number reach 257-512 read intervals,
  - no current window exceeds 512 read intervals.

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
