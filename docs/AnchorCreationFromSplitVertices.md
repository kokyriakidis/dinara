# Anchor Creation from Split Vertices

## Problem Statement

The marker graph vertex construction in `createMarkerGraphVertices` uses a
disjoint-set (union-find) data structure to merge markers across alignments.
The union operation is **transitive**: if read A aligns with read B at some
marker position, and read B aligns with read C at the same marker position,
then A, B, and C all end up in the same vertex — even if A and C share **no
alignment** between them.

This means a single marker graph vertex can contain reads that do not actually
overlap at that genomic position. When such a vertex becomes a Shasta2 anchor,
the anchor's read set is incoherent: it claims reads pass through the same
locus when some of them have no evidence of doing so. This causes downstream
problems in journey construction, anchor graph topology, and assembly.

### Concrete example

```
Read A  ──────────────────
Read B       ──────────────────
Read C            ──────────────────

Alignments: A↔B (at marker m), B↔C (at marker m)
No alignment: A↔C
```

The disjoint-set union merges markers from A, B, and C into one vertex because
A↔B and B↔C are both present. But A and C never aligned — B is the "bridge"
read that caused the transitive merge. The resulting vertex contains reads
spanning a wider genomic region than any single overlap supports.

## Solution: Overlap-Connected Component Splitting

We split each marker graph vertex into **overlap-connected components** using
the read graph. Each component becomes a separate anchor containing only reads
that are mutually connected through direct overlaps.

### Data available at splitting time

| Data structure | Description |
|---|---|
| `markerGraph.vertices` | For each vertex: list of MarkerId values |
| `markerGraph.reverseComplementVertex` | Maps each vertex to its RC partner |
| `findMarkerId(markerId, markers)` | Converts MarkerId → (OrientedReadId, ordinal) |
| `readGraph.connectivity[orientedReadId]` | Read graph adjacency: edge IDs for each oriented read |
| `readGraph.edges[edgeId]` | Edge details: the two OrientedReadIds, crossesStrands flag |
| `reads.getReadRawSequenceLength(readId)` | Raw sequence length of a read |

The read graph at this point contains only **cis overlaps**. The phasing
pipeline classifies every overlap as either cis (state 1) or trans (state 2)
— it never leaves overlaps unlabeled (state 0). Then
`createReadGraphFromPhasingCisOverlaps` removes any edge where either side
is trans. The result is that every edge in the read graph is confirmed cis
from both read perspectives. No additional cis/trans filtering is needed
during splitting.

## Algorithm

### Overview

For each **canonical vertex** (one per reverse-complement pair), we:

1. Collect the vertex's oriented reads
2. Sort them by read length (longest first)
3. Find "root" reads that are not connected to each other (Pass 1)
4. Assign remaining reads to root groups (Pass 2)
5. Derive the reverse-complement vertex's split from the canonical split
6. Apply coverage filters and emit anchors

If a vertex has only one root, it is coherent and passes through unchanged.
This is the common case for the vast majority of vertices.

### Step 0: Collect oriented reads

Extract the set of oriented reads from the vertex's marker IDs using
`findMarkerId`. Deduplicate by OrientedReadId (the marker graph can have
multiple markers from the same oriented read in a vertex; we keep one).

### Step 1: Sort by raw read length (descending)

Sort oriented reads by `reads.getReadRawSequenceLength(readId)`, longest
first. Ties are broken by OrientedReadId value for determinism.

**Rationale:** Longer reads are more likely to be genuine genomic reads and
provide better "root" anchors. Starting from the longest ensures that the
root reads represent the most informative overlap clusters.

### Step 2: Precompute in-vertex neighbors

For each read in the vertex, compute its **in-vertex neighbors**: the
intersection of its read-graph neighbors with the set of reads in this vertex.
This is done by iterating over `readGraph.connectivity[orientedReadId]` and
checking membership in a sorted vector of the vertex's reads via binary
search. The resulting neighbor lists are also stored as sorted vectors for
cache-friendly binary search in subsequent passes.

This precomputation avoids redundant read-graph lookups in the subsequent
passes.

### Step 3: Find roots (Pass 1)

Process reads in length-descending order. Maintain a list of roots.

The **first read (longest)** is always added as a root immediately — there
are no existing roots to check against, so the overlap check is trivially
false. This avoids a redundant empty-loop iteration.

For each subsequent read R:
- Check if R has a read-graph overlap with **any existing root** (using the
  precomputed in-vertex neighbors)
- If R overlaps no existing root: **R becomes a new root**
- If R overlaps one or more roots: skip R (handled in Pass 2)

After this pass, the roots are the "seed" reads that represent independent
overlap clusters within the vertex. By construction, no two roots overlap
with each other.

**If only one root is found:** The vertex is coherent — all reads are
reachable from the single root through overlaps. No splitting needed. This
is the fast path for the vast majority of vertices.

### Step 4: Assign remaining reads (Pass 2)

Process non-root reads in the same length-descending order. For each
unassigned read R:

1. Check which **roots** R has a **direct** read-graph overlap with. Only
   direct edges between R and a root count — transitive connections through
   other already-assigned reads do not.

2. Among the directly overlapping roots, assign R to the one with the
   **longest raw sequence length**. Ties are broken by group index.

Every read in the vertex must overlap at least one root. The vertex was
built using the cis-only read graph (constructed after phasing), so every
read was merged into the vertex through cis overlaps that still exist.
An assert verifies this invariant at runtime.

**Why direct overlap only?** Using transitive assignment (where overlap
with any already-assigned read counts) would risk re-merging groups that
should stay separate — exactly the transitive-closure problem we are
trying to fix. By requiring a direct edge to a root, each non-root read
is anchored to a specific overlap cluster based on its own alignment
evidence.

**Why single-assignment?** A read that overlaps multiple roots is the
"bridge" read that caused the transitive merge in the first place. Assigning
it to only one group is correct because the whole point of splitting is to
break these transitive chains. The read goes to the root with the longest
raw sequence length because that root is most likely to represent the
primary genomic locus.

### Step 5: Reverse complement consistency

The marker graph maintains reverse complement symmetry: for every vertex V,
there exists a vertex V' = `reverseComplementVertex[V]` containing the RC
of every read in V. When we split V into groups {G0, G1, ...}, we must
split V' identically.

We process only **canonical vertices** (where vertexId ≤ rcVertexId),
using the precomputed `markerGraph.reverseComplementVertex` mapping (a
single array lookup per vertex).

For each canonical vertex's split, the RC groups are derived by flipping
each read's strand and computing the RC ordinal directly:

```
rcOrientedReadId = orientedReadId with flipped strand
rcOrdinal = markers[orientedReadId].size() - 1 - ordinal
```

This avoids iterating the RC vertex's markers entirely. The RC ordinal
formula follows from how reverse complement markers are defined in the
marker graph (see `Assembler::findReverseComplement`).

Self-complementary vertices (V = V') are handled as a special case: only
one set of groups is emitted.

### Step 6: Coverage filter and anchor emission

Each group becomes a candidate anchor. Groups with fewer than
`minAnchorCoverage` or more than `maxAnchorCoverage` unique oriented reads
are discarded.

## Complexity

### Per vertex

- Collecting reads: O(m) where m = number of markers in the vertex
- Sorting: O(n log n) where n = number of unique oriented reads
- Precomputing in-vertex neighbors: O(n × d) where d = average read-graph
  degree
- Root finding: O(n × r) where r = number of roots (typically 1-3)
- Assignment: O(n × r) overlap lookups via the precomputed neighbor sets

Total per vertex: **O(n × d)** dominated by the neighbor precomputation.

### Total across all vertices

O(V × n × d) where V = number of vertices. However:

- Most vertices have r = 1 (coherent), so the splitting logic is a no-op
  after root finding
- n is bounded by `maxAnchorCoverage` (typically 5 × coveragePeak ≈ 150)
- Only canonical vertices are processed (V/2)

The algorithm is **embarrassingly parallel** across vertices since each
vertex is processed independently. The implementation uses a two-pass
multithreaded approach with dynamic load balancing (atomic batch counter,
batch size 1000):

- **Pass 1:** Each thread splits vertices and records per-anchor read
  counts. A prefix sum over per-vertex anchor counts assigns deterministic
  anchor IDs.
- **Pass 2:** Each thread re-splits vertices and fills `anchorMarkerInfos`
  directly. Each thread writes to non-overlapping anchor ID ranges, so no
  synchronization is needed.

`splitVertex` is called exactly twice per canonical vertex (once per pass).
No intermediate anchor storage is needed — Pass 1 stores only `uint64_t`
sizes, and Pass 2 writes directly into the memory-mapped VectorOfVectors.

## Integration

### Pipeline position

The splitting runs **after** marker graph vertex construction and RC vertex
finding, and **replaces** the standard `Shasta2Anchors` constructor call:

```
createMarkerGraphVertices(...)
findMarkerGraphReverseComplementVertices(...)

// OLD: direct anchor construction
// assembler.shasta2Anchors = make_shared<Shasta2Anchors>(...);

// NEW: anchor construction with vertex splitting
assembler.shasta2Anchors = createShasta2AnchorsFromSplitVertices(
    shasta2Owner,
    assembler.getReads(),
    assembler.assemblerInfo->k,
    *assembler.markers,
    assembler.markerGraph,
    assembler.readGraph,    // <-- additional input
    threadCount,
    minAnchorCoverage,
    maxAnchorCoverage);
```

### Implementation files

| File | Role |
|---|---|
| `src/Shasta2AnchorsFromSplitVertices.hpp` | Declares `createShasta2AnchorsFromSplitVertices()` |
| `src/Shasta2AnchorsFromSplitVertices.cpp` | Implements the splitting algorithm |
| `src/Shasta2Anchors.hpp` | Added minimal constructor declaration |
| `src/Shasta2Anchors.cpp` | Added minimal constructor implementation |
| `src/AssemblerOptions.cpp` | Changed default `Reads.representation` from 1 (RLE) to 0 (raw) |

The function uses the minimal `Shasta2Anchors` constructor, which only
initializes member variables and builds `MarkerKmers`. It does not build
`anchorMarkerInfos` — the splitting function populates that directly from
the split groups, avoiding any redundant vertex iteration.

### Output

The function returns a `shared_ptr<Shasta2Anchors>` that is a drop-in
replacement for the standard constructor's output. Downstream code
(journey construction, anchor graph, assembly graph) works unchanged.

**Note:** The `anchorVertexIds` vector (mapping anchor index → original
marker graph vertex ID) is not populated since the minimal constructor is
used and the 1:1 vertex-to-anchor correspondence no longer holds after
splitting. If any downstream code depends on this mapping, it would need
updating.

## Diagnostics

The splitting logs the following statistics:

- **Canonical vertices processed** — number of vertices examined (half of total, one per RC pair)
- **Unsplit anchors** — anchors from coherent vertices (including their RC counterparts)
- **Vertices split** — number of vertices that had > 1 root
- **Anchors from splits** — total anchors produced from split vertices (including their RC counterparts)
- **Discarded** — groups that fell outside `[minAnchorCoverage, maxAnchorCoverage]` (not included in any other count)
- **Total anchors** — final anchor count (`unsplit + from splits`)

Example output:
```
Processed 5000000 canonical vertices.
Unsplit anchors: 9990000
Split 5000 vertices into 12000 anchors.
Discarded 300 groups outside coverage range [8, 150].
Total anchors: 10002000
```

## Future refinements

### Bridge read handling

Currently, a bridge read (one that overlaps multiple roots) is assigned to
the group with the longest root. An alternative strategy would be to
duplicate the bridge read into all overlapping groups, since it genuinely
spans the boundary. This would increase anchor coverage but could
reintroduce some of the transitive-merge artifacts we are trying to
eliminate.

### Minimum root separation

An optional refinement: require that roots be separated by a minimum
genomic distance (estimated from marker positions) to avoid splitting
vertices where the roots are very close together and the lack of a direct
alignment is likely due to alignment sensitivity rather than genuine
non-overlap.

### Iterative splitting

After the initial split, some groups may still contain transitively-merged
reads (if the bridge read was assigned to a group and then serves as a new
bridge within that group). An iterative application of the splitting
algorithm could catch these cases, though they should be rare in practice.

### Integration with phasing

The splitting relies on the read graph containing only confirmed-cis edges.
The phasing pipeline guarantees this: every overlap is classified as cis or
trans (never unlabeled), and `createReadGraphFromPhasingCisOverlaps` removes
any edge where either side is trans. If phasing improves (e.g., through
het-site detection), the read graph will be cleaner, and the splitting will
produce better results automatically.
