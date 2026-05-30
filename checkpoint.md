# Anchor Window Pipeline — Endpoint Anchor Design

## Two-Pass Inter-Window Edge Creation

The anchor graph constructor uses a two-pass approach for inter-window edges:

- **Pass 1 (endpoint edges):** For each backbone transition (`backbonePreviousWindow` / `backboneNextWindow`), create the single best-sharing edge between the two windows. Reserve the anchors used (+ RC mirrors) so pass 2 cannot reuse them.
- **Early trim:** After pass 1, disable backbone anchors beyond the endpoint positions. This constrains the graph to the region between endpoints.
- **Pass 2 (internal edges):** Create edges for all remaining inter-window anchor pairs, skipping reserved anchors.

## Endpoint Anchor Invariants

Each window has **at most 2 endpoint anchors** — one at the head (connecting to `backbonePreviousWindow`) and one at the tail (connecting to `backboneNextWindow`). Fewer if the window is at the start/end of a backbone chain.

Key properties:
1. **Never filtered.** All filters (singleEdge, bypass/detour, bubble, spur) check `isEndpointAnchorPrev || isEndpointAnchorNext` and skip those edges.
2. **Single edge per side.** Each endpoint anchor has exactly one endpoint edge, connecting to exactly one adjacent window.
3. **Highest sharing.** Pass 1 selects the anchor pair with the most shared reads for each endpoint window pair.
4. **Internal edges are bounded.** After early trim, internal edges exist strictly between the two endpoint anchors of each window.

## Per-Anchor Endpoint Flags

Each edge carries two flags (not serialized):
- `isEndpointAnchorPrev` — source anchor is at the first or last active backbone position of its window.
- `isEndpointAnchorNext` — target anchor is at the first or last active backbone position of its window.

These are set after all edges are created, using position-based detection (first/last active backbone anchor per window). The flags drive both GFA tag output and filter protection.

## GFA Tags

Each link line gets `pw:Z:` and `nw:Z:` tags:
- `Endpoint` — anchor is an endpoint anchor
- `internal` — inter-window but not an endpoint anchor
- `intra` — both anchors in the same window

## Removed: `recomputeBackboneEndpoints`

The `recomputeBackboneEndpoints` lambda was removed. `backbonePreviousWindow` / `backboneNextWindow` are set once during construction and do not change — the backbone read's path through windows is a property of the read, not the current edge state.

## Files Changed

- `src/Shasta2AnchorGraph.hpp` — added `isEndpointAnchorPrev`, `isEndpointAnchorNext` to edge; added `endpointWindowPairs`, `endpointAnchors` to graph.
- `src/Shasta2AnchorGraph.cpp` — two-pass construction, early trim, position-based endpoint flags, updated filters, removed `recomputeBackboneEndpoints`.
- `src/Shasta2AnchorGraphGfa.cpp` — `writeGfa` accepts `anchorWindows`, emits `pw:Z:`/`nw:Z:` tags using per-anchor flags.
- `srcMain/main.cpp` — updated `writeGfa` call sites to pass `&anchorWindows`.
