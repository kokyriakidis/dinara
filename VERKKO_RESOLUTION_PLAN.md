# Plan: Verkko-style Directed Graph Resolution for Dinara

## Executive Summary

This document lays out a comprehensive plan to add a **new resolution mode** to Dinara that follows Verkko/MBG's approach: build a directed graph of anchors (analogous to MBG unitigs), thread read paths through it, then **iteratively resolve** tangles using triplet evidence from the paths. The existing BRG anchor system remains untouched for backward compatibility.

---

## 1. Verkko's Full Pipeline (Reference)

Verkko's assembly pipeline, as defined in its Snakefile, proceeds through these major stages:

| Stage | Snakefile | What it does |
|-------|-----------|-------------|
| **c1-c4** | `c1-buildStore`, `c2-findOverlaps`, `c4-findErrors` | Canu-based HiFi error correction |
| **1-buildGraph** | `1-buildGraph.sm` | MBG builds a multiplex de Bruijn graph → GFA + read paths (GAF) |
| **2-processGraph** | `2-processGraph.sm` | HiFi-path-based resolution: `resolve_triplets_kmerify.py` → `pop_bubbles` → `estimate_unique` → `connect_uniques` → `forbid_unbridged` → `merge_unresolved` → `unitigify` |
| **3-align** | `3-splitONT`, `3-alignONT`, `3-combineONT`, `3-alignTips` | Align ONT reads to the HiFi-resolved graph (GraphAligner) |
| **4-processONT** | `4-processONT.sm` | ONT-path-based resolution: same scripts as stage 2, with ONT paths and adjusted coverage thresholds |
| **5-untip** | `5-untip.sm` | Final tip clipping: `untip_relative.py`, `chop_misassemblies.py`, coverage injection |
| **6-rukki** | `6-rukki.sm` | Haplotype-aware path walking (trio/Hi-C/Pore-C) |
| **6-layoutContigs** | `6-layoutContigs.sm` | Contig layout from resolved graph |
| **7-consensus** | `7-*` | Extract reads, build packages, generate consensus |

The **graph flows** through these named files (tracing the data lineage):
```
1-buildGraph/              → unitig-unrolled.gfa + .gaf (MBG output)
2-processGraph/            → unitig-unrolled-hifi-resolved.gfa
                           → unitig-unrolled-popped.gfa  (after bubble popping)
                           → unitig-unrolled-popped-unitig-normal-connected.gfa
4-processONT/              → unitig-unrolled-ont-resolved.gfa
                           → unitig-unrolled-unitig-unrolled-popped...gfa
5-untip/                   → unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa
```

---

## APPENDIX A: Complete Verkko Function Call Graph

### A.1 Top-Level Pipeline Flow (resolve_triplets_kmerify.py main)

```
main():
│
├── read_graph(input_gfa)                           # L19-37: Parse GFA → (base_seqs, node_seqs, edges, edge_overlaps)
│   ├── For each S line: base_seqs[name]=seq, node_seqs[name]=([">name"], 0, 0)
│   └── For each L line: edges[from].add(to) + edges[revnode(to)].add(revnode(from))
│                         edge_overlaps[canon(from,to)] = overlap_bp
│
├── Parse paths from stdin → add_path(paths_crossing, path)   # L1010-1031
│   └── add_path()                                  # L822-827: For each node in path, paths_crossing[node][id(path)] = path
│
├── read_node_coverages(file)                       # L717-727: Parse TSV → {node: float_coverage}
│
├── Mark removable_nodes (coverage < min_allowed)   # L1035-1037
│
├── remove_and_split_low_coverage(...)              # L653-716
│   ├── Compute edge_coverage from initial_paths
│   ├── get_safe_unitigs_and_edges(...)             # L644-651: Walk unitigs from high-cov nodes
│   │   └── add_safe_unitig(node, edges, ...)       # L632-641: Forward walk while deg=1
│   ├── Remove unsafe low-cov edges
│   ├── Mark removable nodes (not safe, low-cov)
│   ├── Split paths at removed nodes/edges
│   └── remove_graph_node(node, node_seqs, edges)   # L123-144: Delete node, clean bidirectional edges
│
├── unitigify_all(...)                              # L986-991: Initial compaction
│   └── For each node: unitigify_one(...)           # L972-984
│       ├── extend_forward(">" + node, edges)       # L954-970: Walk fwd while in/out-deg=1
│       ├── extend_forward("<" + node, edges)       # Walk backward similarly
│       └── replace_unitig(unitig_chain)            # L829-953: Merge chain into single node
│
├── FOR EACH coverage IN resolve_steps:             # L1041-1044 (e.g., [20, 10, 5])
│   └── resolve(node_lens, edge_overlaps, node_seqs, edges, paths_crossing, coverage, min_cov, removable)
│       └── [See A.2 below]
│
├── Assign unitig names                             # L1046-1058
├── Write paths to out_path_file                    # L1060-1067
├── Write edges (GFA L lines) to stdout             # L1075-1077
├── read_graph_only_bases(input_gfa)                # Re-read base sequences
└── Write sequences (GFA S lines) to stdout         # L1079-end: get_seq() for composite nodes
```

### A.2 resolve() — The Main Resolution Loop (L729-772)

```
resolve(node_lens, edge_overlaps, node_seqs, edges, paths_crossing, min_edge_support, min_coverage, removable_nodes):
│
├── maybe_resolvable = set(node_seqs.keys())
├── nodes_by_len = min-heap [(get_unitig_len(n), n) for n in maybe_resolvable]
│   └── get_unitig_len(node_lens, edge_overlaps, node_seqs, node)  # L774-777
│       └── get_unitig_len_path(node_lens, edge_overlaps, node_seqs[node][0], ...)  # L779-805
│
├── WHILE heap not empty AND current_length <= max_resolve_length:
│   │
│   ├── Pop all nodes with same length → current_nodes set
│   │
│   ├── TRY HAIRPINS FIRST:
│   │   (new_nodes, resolved) = resolve_hairpins(current_length, current_nodes, ...)
│   │   └── [See A.3 below]
│   │
│   │   IF hairpins resolved:
│   │       FOR each new_node in new_nodes:
│   │           new_unitig = unitigify_one(new_node)
│   │           heappush(nodes_by_len, (get_unitig_len(new_unitig), new_unitig))
│   │           maybe_resolvable.add(new_unitig)
│   │       CONTINUE  # skip resolve_nodes for this length level
│   │
│   ├── TRY NODE RESOLUTION:
│   │   (new_nodes, resolved) = resolve_nodes(current_length, current_nodes, ...)
│   │   └── [See A.4 below]
│   │
│   │   IF nodes resolved:
│   │       FOR each new_node in new_nodes:
│   │           new_unitig = unitigify_one(new_node)
│   │           heappush(nodes_by_len, (get_unitig_len(new_unitig), new_unitig))
│   │           maybe_resolvable.add(new_unitig)
│   │
│   └── CONTINUE to next length level
```

### A.3 resolve_hairpins() — Palindromic Node Resolution (L229-333)

```
resolve_hairpins(nodelength, nodes, paths_crossing, node_seqs, node_lens, edges, maybe_resolvable, ...):
│
├── DETECT HAIRPINS:
│   FOR each node in nodes:
│       IF edges[">node"] == {["<node"]} → hairpin (node's only out-edge goes to its own RC)
│       IF edges["<node"] == {[">node"]} → hairpin (reverse direction)
│       IF edges[">node"] has "<node" AND len >= 2 → unresolvable_hairpin
│
├── FOR each hairpin node (skip double-hairpins, disconnected, unresolvable):
│   │
│   ├── Count resolutions from paths:
│   │   FOR each path crossing node (len >= 4):
│   │       IF path[i] == node AND path[i+1] == revnode(node):  # hairpin traversal
│   │           resolution_key = canon(path[i-1], path[i+2])
│   │           resolutions[key] += 1
│   │
│   ├── Filter: require min_edge_support per resolution
│   ├── Verify: all edges covered by solid resolutions
│   │
│   ├── CREATE fw/bw COPIES for each solid resolution:
│   │   fwname = "node_hairpinN_fw"
│   │   bwname = "node_hairpinN_bw"
│   │   node_seqs[fwname] = node_seqs[node] (or RC if "<")
│   │   node_seqs[bwname] = RC of fwname
│   │
│   ├── WIRE EDGES:
│   │   edges[">fwname"] → {">bwname"}       # fw→bw
│   │   edges["<bwname"] → {"<fwname"}       # RC of above
│   │   edges["<fwname"] → {revnode(key[0])} # predecessor
│   │   edges[">bwname"] → {key[1]}          # successor
│   │   + RC counterparts
│   │   edge_overlaps for fw→bw, pred→fw, bw→succ
│   │
│   ├── remove_graph_node(node) from graph
│   │
│   └── REWRITE PATHS:
│       FOR each path crossing node:
│           Replace "node, revnode(node)" pairs with "fwname, bwname" (or RC)
│           Handle first/last positions specially
│           Split paths if resolution_number not found
│
└── RETURN (new_node_names, resolved_nodes)
```

### A.4 resolve_nodes() — The Main Triplet Resolution Engine (L342-631)

```
resolve_nodes(nodelength, nodes, paths_crossing, node_seqs, node_lens, edges, maybe_resolvable, ...):
│
├── PHASE 1: COLLECT VALID TRIPLETS (L342-352)
│   FOR each node in nodes:
│       triplets_here = get_valid_triplets(node, edges, paths_crossing, ...)
│       └── [See A.5 below]
│       IF empty → remove from maybe_resolvable
│       ELSE → add to resolvable set + collect triplets
│
├── PHASE 2: FILTER UNRESOLVABLE (L353-377)
│   REPEAT until no changes:
│   │   FOR each triplet:
│   │       IF triplet neighbors have effective_length == 1 after overlap subtraction:
│   │           Mark middle node unresolvable  # "borders unresolvable"
│   │       IF neighbor is_hairpin():
│   │           Mark middle node unresolvable
│   │   └── is_hairpin(node_name, edges)  # L335-340: Check if node has self-RC edge
│   │
│   RE-COLLECT triplets only for still-resolvable nodes
│
├── PHASE 3: COMPUTE EXTENSION LENGTHS (L378-396)
│   FOR each resolvable node:
│       longest_extension_per_node[">node"] = min(1000, neighbor_lengths)
│       longest_extension_per_node["<node"] = min(1000, neighbor_lengths)
│       # Extension = how much to extend edge nodes into neighbors
│       # Capped at 1000bp to avoid strange multi-megabase bubbles
│
├── PHASE 4: CREATE EDGE NODES (L397-533)
│   FOR each triplet (from, through, to):
│   │   extend_amount = longest_extension_per_node["<" + through]
│   │   IF from is not None:
│   │       key = (revnode(through), revnode(from))
│   │       nodename = "edge_KEY0_KEY1"
│   │       new_edgenodes[key] = nodename
│   │       node_seqs[nodename] = EXTENDED version of from's sequence
│   │       # Complex sequence merging with overlaps and clips (L418-460)
│   │
│   │   extend_amount = longest_extension_per_node[">" + through]
│   │   IF to is not None:
│   │       key = (through, to)
│   │       nodename = "edge_KEY0_KEY1"
│   │       new_edgenodes[key] = nodename
│   │       node_seqs[nodename] = EXTENDED version of to's sequence
│   │
│   FOR each new edge node:
│       edges[">" + nodename] = set()
│       edges["<" + nodename] = set()
│
├── PHASE 5: WIRE EDGES BETWEEN EDGE NODES (L534-596)
│   FOR each triplet (from, through, to):
│       # Determine left_key and right_key:
│       # If 'from' is a tiny unresolvable node → left_key = from directly
│       # Else → left_key = "<" + new_edgenodes[(rev_through, rev_from)]
│       # Similarly for right_key
│       edges[left_key].add(right_key)
│       edges[revnode(right_key)].add(revnode(left_key))
│       edge_overlaps[canon(left_key, right_key)] = overlap
│       # Also wire to outer neighbors (lefter_key, righter_key)
│
├── PHASE 6: UPDATE PATHS (L596-600)
│   replace_path_nodes(resolvable, paths_crossing, new_edgenodes, ...)
│   └── [See A.6 below]
│
├── PHASE 7: CLEAN UP (L601-631)
│   Collect new_paths from edge nodes
│   FOR each resolvable node:
│       remove_graph_node(node, node_seqs, edges)
│   Split any paths with broken edges
│
└── RETURN (new_node_names, resolvable)
```

### A.5 get_valid_triplets() — Triplet Counting and Validation (L150-227)

```
get_valid_triplets(node, edges, paths_crossing, min_edge_support, min_coverage, removable_nodes, node_seqs):
│
├── EARLY EXITS:
│   IF no edges on both sides → return []
│   IF out-degree <= 1 AND in-degree <= 1 → return []
│   IF only edge is self-RC → return []
│   IF node not in paths_crossing → return []
│
├── COUNT TRIPLETS from paths:
│   FOR each path crossing node:
│       FOR each position in path:
│           IF path[i] == ">node":
│               covered_in_neighbors[path[i-1]] += 1
│               covered_out_neighbors[path[i+1]] += 1
│               triplet = (path[i-1], ">node", path[i+1])
│           IF path[i] == "<node":
│               # Reverse: covered_out_neighbors gets revnode(path[i-1])
│               #          covered_in_neighbors gets revnode(path[i+1])
│               triplet = (revnode(path[i+1]), ">node", revnode(path[i-1]))
│           # Handle path endpoints: None for missing predecessor/successor
│           triplets[triplet] += 1
│
├── FILTER by min_edge_support:
│   solid_triplets = {t for t in triplets if triplets[t] >= min_edge_support}
│
├── VALIDATE COVERAGE of neighbors:
│   FOR each out-edge of ">node":
│       IF edge NOT covered by any solid triplet:
│           # Check if removable (all parts are low-coverage)
│           IF not removable → return []  # Can't resolve without covering all edges
│   FOR each in-edge of "<node":
│       Same check for in-direction
│
└── RETURN sorted list of solid triplets
```

### A.6 replace_path_nodes() — Path Rewriting After Resolution (L50-86)

```
replace_path_nodes(resolvable, paths_crossing, new_edgenodes, maybe_resolvable, edges):
│
├── Collect all paths crossing any resolvable node (deduplicated by id())
│
├── FOR each affected path:
│   Build new_path:
│       FOR each position j in path:
│           IF path[j] not in resolvable → keep it
│           IF (revnode(path[j]), revnode(path[j-1])) in new_edgenodes:
│               → append "<" + new_edgenodes[...]  # backward edge node
│           IF (path[j], path[j+1]) in new_edgenodes:
│               → append ">" + new_edgenodes[...]  # forward edge node
│
├── Remove old paths from paths_crossing
│
├── FOR each new path:
│   Split at broken edges (where consecutive nodes lack an edge)
│   Add valid segments via add_path()
│
└── Assert: no paths remain crossing resolved nodes
```

### A.7 replace_path_node() — Single-Node Hairpin Path Rewriting (L87-121)

```
replace_path_node(paths_crossing, node, left_added, right_added):
│
├── FOR each path crossing node:
│   Build new_path:
│       FOR each position j:
│           IF path[j] != oriented_node → keep it
│           IF path[j] == ">node":
│               IF j > 0: new_path.append(">" + left_added[path[j-1]])
│               IF j < end: new_path.append(">" + right_added[path[j+1]])
│           IF path[j] == "<node":
│               IF j > 0: new_path.append("<" + right_added[revnode(path[j-1])])
│               IF j < end: new_path.append("<" + left_added[revnode(path[j+1])])
│
├── Remove old paths, add new paths (with splitting at broken edges)
```

### A.8 Unitigification Functions (L829-991)

```
unitigify_all(node_seqs, node_lens, edges, paths_crossing):    # L986-991
│   FOR each node (deterministic order):
│       IF node still exists → unitigify_one(node)

unitigify_one(node_seqs, node_lens, edges, paths_crossing, node):  # L972-984
│   forward_ext = extend_forward(">" + node, edges)
│   IF not circular → backward_ext = extend_forward("<" + node, edges)
│   IF total chain length == 1 → return node (nothing to merge)
│   unitig = reversed(backward_ext) + forward_ext  (deduplicated)
│   IF circular: remove duplicate endpoint
│   return replace_unitig(unitig_chain)

extend_forward(node, edges):                        # L954-970
│   result = [node]
│   WHILE True:
│       IF out-degree != 1 → break
│       IF palindrome hairpin (self-RC) → break
│       next = getone(edges[pos])
│       IF in-degree of next != 1 → break
│       IF circular (next == start) → append start, break
│       IF circular RC → break
│       result.append(next)
│   RETURN result

replace_unitig(node_seqs, node_lens, edges, paths_crossing, unitig):  # L829-953
│   new_node = "unitig_N" (global counter)
│   MERGE SEQUENCES:
│       For each node in chain:
│           Handle strand (forward/reverse of node_seqs)
│           Find overlap between consecutive elements
│           Concatenate non-overlapping suffix
│           Track start_clip and end_clip
│   node_seqs[new_node] = (merged_seq, start_clip, last_end_clip)
│
│   REWIRE EDGES:
│       edges[">" + new_node] = edges[unitig[-1]]  (last node's out-edges)
│       edges["<" + new_node] = edges[revnode(unitig[0])]  (first node's in-edges)
│       Update edge_overlaps for new edges
│
│   REMOVE old nodes: remove_graph_node() for each in chain
│
│   REWRITE PATHS:
│       Collect all paths touching any chain node
│       Build is_forward map for strand tracking
│       Replace contiguous runs of chain nodes with single new_node
│       Handle circular unitigs at boundaries
│       Remove old paths, add new paths via add_path()
│
│   RETURN new_node
```

### A.9 Low-Coverage Removal (L632-716)

```
remove_and_split_low_coverage(node_seqs, edges, initial_paths, paths_crossing, min_coverage, node_coverage):
│
├── Compute edge_coverage from initial_paths
│
├── get_safe_unitigs_and_edges(node_seqs, edges, min_coverage, node_coverage, edge_coverage):
│   └── add_safe_unitig(node, edges, safe_nodes, safe_edges):
│       Walk forward while out-deg=1, in-deg=1 (of next), mark as safe
│
├── Remove low-cov EDGES (not safe, cov < min_coverage)
│   edges[from].remove(to) + edges[revnode(to)].remove(revnode(from))
│
├── Mark removable NODES (not safe, cov < min_coverage)
│
├── SPLIT PATHS at removed nodes/edges:
│   FOR each initial_path:
│       Walk through, split at any removed node or broken edge
│       Add valid segments via add_path()
│
└── remove_graph_node() for each removable node
```

### A.10 graph_functions.py — Complete Function Index

```
GRAPH PRIMITIVES:
  revnode(n)                          # L173-176: ">" ↔ "<" prefix swap
  canon(left, right)                  # L189-193: Lexicographic canonical edge pair
  canontip(left, right)               # L195-198: Canonical pair without revnode
  getone(s)                           # L200-202: Get any element from a set
  pathstr(p)                          # L204-205: Join path to string
  iterate_deterministic(l, end="")    # L207-217: Sort + yield for deterministic iteration
  find(parent, key)                   # L178-181: Union-Find with path compression
  merge(parent, left, right)          # L183-186: Union-Find merge
  nor_node(node)                      # L341: Strip +-
  nor_path_id(path_id)               # L344: Strip +-
  revcomp(s)                          # L167-168: Reverse complement DNA sequence
  str2bool(v)                         # L170-171: String to boolean
  rc_seq(seq)                         # L397-398: Reverse complement string
  rc_orientation(c)                   # L400-404: + ↔ -
  rc_path_id(path_id)                # L406-407: Reverse path ID orientation
  rc_path(path)                       # L409-423: Reverse entire path with strand flip

GRAPH I/O:
  load_direct_graph(gfa_file, G)      # ~L100-130: Load GFA into NetworkX DiGraph (nodes with length+coverage)
  load_indirect_graph(gfa_file, G)    # ~L80-100: Load GFA into NetworkX undirected Graph
  get_lengths(fasta_file)             # ~L132-145: Parse FASTA → {name: length}

GRAPH ANALYSIS:
  get_component_length(G, component)  # L4-8: Sum node lengths in component
  nodes_in_tangles(G, MAX_LEN, MIN)   # L57-68: Find short-node clusters (tangles)
  remove_large_tangles(G, ...)        # L14-55: Remove large short-node components (rDNA)
  topological_sort(nodelens, edges)   # L280-308: Tarjan's SCC → topological order
  strong_connect_iterative(...)       # L219-278: Iterative Tarjan's SCC algorithm
  getComponentColors(G)               # L347-354: Assign color ID per connected component
  get_telomeric_nodes(telo_file, G)   # L359-395: Add telomere pseudo-nodes to graph

SCAFFOLDING:
  loadHiCGraph(hic_byread_file)       # ~L310-339: Load Hi-C data (deprecated?)
  tsv2gaf(tsv_path)                   # L70-78: Convert TSV path format to GAF
```

### A.11 Supporting Scripts — Function Index

#### pop_bubbles_coverage_based.py (Bubble Popping)
```
FUNCTIONS:
  find_component_coverage(key, ...)      # Component average coverage via union-find
  merge(parent, cov_sum, len_sum, l, r)  # Union-find merge with coverage tracking
  remove_graph_node(node, edges)         # Remove node + clean edges
  find_bubble(s, edges, nodelens, max)   # Superbubble detection (Onodera 2013, fig. 5)

ALGORITHM:
  1. Load graph, node coverages
  2. Compute per-component coverage (NetworkX connected components)
  3. Find superbubbles: walk from each ">node", track visited/seen
  4. Detect tips (dead-end nodes not in bubbles)
  5. Union-find: merge nodes across bubbles/tips into chains
  6. Classify chains: unique (0.5x-2.5x component coverage), tip (≤2.5x)
  7. Output: nodes in unique_chains

PARAMS: max_bubble_pop_size=10, max_poppable_node_size=200000, max_coverage_delta=1.5
```

#### estimate_unique_local.py (Unique Node Estimation)
```
FUNCTIONS:
  merge(parent, rank, left, right)       # Union-find merge
  find_bubble_end(edges, s)              # Superbubble detection (same Onodera algorithm)
  is_tip_unique_candidate(edges, s)      # Check if tip's neighbors are also tips

ALGORITHM:
  1. Load graph, coverages, alignments
  2. Identify long nodes (> threshold), tips, bubbles
  3. Group nodes by bubble/tip relationships via union-find
  4. Compute roughly_average_coverage_nodes: coverage within 0.5x-2.5x global average
  5. Path consistency analysis:
     - For each node, check fw/bw paths from alignment file
     - Count consistent paths (paths that agree with graph topology)
     - Require path_consistency_threshold fraction consistent
  6. Final unique = path_consistent ∩ roughly_average_coverage

PARAMS: long_node_threshold, solid_edge_threshold, path_consistency_threshold
```

#### connect_uniques.py (Tangle Bridging)
```
FUNCTIONS:
  merge(parent, left, right)             # Union-find for tangle grouping

ALGORITHM:
  1. Load graph, unique nodes, forbidden ends, resolving paths
  2. Union-find: group connected non-unique nodes into tangles
  3. Identify resolvable_ends: path endpoints at unique nodes
  4. Remove forbidden ends
  5. For each resolvable tangle:
     - For each path connecting unique ends through the tangle:
       Create node copies, wire-in copies, update paths
  6. Output new graph with resolved tangles

KEY DETAIL: Uses gf.find() and gf.canontip() for tangle grouping
```

#### forbid_unbridged_tangles.py (Prevent False Resolution)
```
FUNCTIONS:
  union(s1, s2, parent, rank)            # Union-find

ALGORITHM:
  1. Load unique nodes, graph, connections, paths, ONT/HiFi coverages
  2. Compute edge_coverage from paths
  3. Identify solid nodes (length > 200Kbp)
  4. Union-find: group nodes into tangles
  5. For each tangle connecting unique nodes:
     - Check if bridged (sufficient path support)
     - Check coverage (ONT >= min_ont_solid, HiFi >= min_hifi_solid)
     - If unbridged → add to forbidden connections
  6. Output: forbidden connection pairs

PARAMS: length_solid_node_threshold=200000, min_ont_solid_coverage, min_hifi_solid_coverage
```

#### merge_unresolved_dbg_nodes.py (Undo Failed Resolutions)
```
FUNCTIONS:
  remove_graph_node(node, nodeseqs, edges)
  get_base_name(name)                    # Extract original node name from resolved name

ALGORITHM:
  1. Load graph from stdin
  2. Group nodes by base name (before resolution suffixes)
  3. For each group with multiple copies:
     - Check if copies share same in-neighbor set
     - If yes → merge back into single node (undo resolution)
  4. Rename merged nodes, output GFA

KEY DETAIL: Uses gf.canontip() for edge overlap mapping
```

#### untip_relative.py (Relative Tip Removal)
```
FUNCTIONS:
  get_node_depths(order, belongs_to_component, nodelens, edges)
  remove_rec(kept, node, edges)          # Recursive removal from "kept" set
  get_keepers(order, belongs_to_component, depths, nodelens, edges)

ALGORITHM:
  1. Topological sort via gf.topological_sort()
  2. Compute depths: each node's max-depth descendant length
  3. Walk topological order:
     - For each node, compare its depth to neighbor depths
     - If depth < fraction * neighbor_depth AND length < max_removable_len:
       → Remove (it's a tip relative to its deeper neighbor)
  4. Output kept nodes as GFA

PARAMS: max_removable_len, min_safe_len, fraction
```

#### chop_misassemblies.py (Break at Misassembly Points)
```
ALGORITHM:
  1. Analyze read alignment coverage across nodes
  2. Find coverage drops / clustering breakpoints
  3. Split nodes at breakpoints
  4. Output modified GFA
```

#### unroll_tip_loops.py (Unroll Loops at Tips)
```
ALGORITHM:
  1. Load graph, coverages, paths
  2. Compute average coverage (weighted by node length)
  3. Identify tip loops: dead-end nodes with self-edges
  4. If loop coverage suggests > 1 copy: unroll
  5. Max unroll length = 200000bp
```

---

## 2. Key Algorithms (Detailed)

### 2.1 Data Structures (from `graph_functions.py` and `resolve_triplets_kmerify.py`)

**Node orientation convention:**
- `">nodeName"` = node in forward orientation
- `"<nodeName"` = node in reverse orientation
- `revnode(">X")` = `"<X"` and vice versa
- `canon(from, to)` = lexicographic canonical edge pair for deduplication

**Core mutable state:**
```python
node_seqs[name]     = (list_of_oriented_kmer_nodes, left_clip, right_clip)
edges[">name"]      = set(["<Y", ">Z", ...])  # directed adjacency
edge_overlaps[key]  = int  # bp overlap, keyed by canon(from, to)
paths_crossing[name]= {id(path): path_list}  # all paths that cross node "name"
node_lens[name]     = int  # bp length
```

**Invariant:**  edges are always stored bidirectionally:
```python
edges[from].add(to)
edges[revnode(to)].add(revnode(from))
```

### 2.2 resolve_triplets_kmerify.py — The Core Resolution Algorithm

**Main flow** (lines 1023–1044):
1. `read_graph(input_gfa)` → parse GFA into `(base_seqs, node_seqs, edges, edge_overlaps)`
2. Read paths from stdin (GAF/path format) → `add_path(paths_crossing, path)`
3. `read_node_coverages(file)` → per-node coverage
4. Mark `removable_nodes` = nodes with coverage < `min_allowed_coverage`
5. `remove_and_split_low_coverage()` — remove low-cov nodes, split paths at breakpoints
6. `unitigify_all()` — initial unitigification (merge linear chains)
7. **For each coverage threshold** in `resolve_steps` (e.g., `[20, 10, 5]`):
   - `resolve(node_lens, edge_overlaps, node_seqs, edges, paths_crossing, coverage, min_allowed_cov, removable_nodes)`
8. Assign unitig names, write paths, write GFA output

**`resolve()` function** (lines 729–772):
- Min-heap by unitig length. For each length level:
  1. Try `resolve_hairpins()` first — handle palindromic nodes (self-RC edges)
  2. If hairpins resolved → `unitigify_one` new nodes → push back to heap
  3. Try `resolve_nodes()` — the main triplet resolution
  4. If nodes resolved → `unitigify_one` new nodes → push back to heap
  5. Continue until `max_resolve_length` reached

**`get_valid_triplets(node)` (lines 150–227):**
- For a node, look at all paths crossing it via `paths_crossing[node]`
- Count (predecessor, node, successor) triplets
- Filter: require `min_edge_support` occurrences for each triplet
- Verify all in/out neighbors are covered by at least one triplet
- Check coverage thresholds, skip removable nodes
- Return list of `(left, through, right)` triplets (left/right can be `None` for tips)

**`resolve_nodes()` (lines 342–631):**
1. For each candidate node, call `get_valid_triplets()`
2. Filter out unresolvable: nodes with hairpin neighbors at borders
3. For each valid `(from→through→to)` triplet, create an "edge node" named `"edge_FROM_TO"`
4. Set up edge node sequences by extending the through-node's kmer sequence
5. Wire edges: `left_key → right_key` + reverse complement edges
6. `replace_path_nodes()` — update all paths that cross resolved nodes
7. Remove resolved graph nodes
8. Split paths at broken edges

**`unitigify_one(node)` (lines 972–984):**
- `extend_forward(">" + node)` — walk forward while in-degree=1 and out-degree=1
- `extend_forward("<" + node)` — walk backward similarly
- Merge the linear chain via `replace_unitig()`

**`replace_unitig(unitig_chain)` (lines 829–953):**
- Create a new node merging the chain of oriented nodes
- Handle sequence overlaps/clips, rewire edges, update all crossing paths

### 2.3 Supporting Scripts

| Script | Purpose | Key Algorithm |
|--------|---------|--------------|
| `pop_bubbles_coverage_based.py` | Remove bubble/tip alleles | Superbubble detection (Onodera 2013), union-find chain merging, coverage ratio filtering (0.5x–2.5x) |
| `estimate_unique_local.py` | Identify "unique" (non-repetitive) nodes | Superbubble + tip analysis, union-find grouping, uses `long_node_threshold` and `solid_edge_threshold` |
| `connect_uniques.py` | Bridge tangles using paths that connect unique nodes | Union-find on graph edges, create node copies for each path insertion through tangles |
| `forbid_unbridged_tangles.py` | Prevent false tangle resolution | Edge coverage from paths, solid nodes (>200Kbp), ONT/HiFi coverage thresholds |
| `merge_unresolved_dbg_nodes.py` | Merge failed resolution copies | Group resolved copies of same base node, merge back if same in-neighbor set |
| `unitigify.py` | Standalone unitigification | `start_unitig()` / `start_circular_unitig()`, outputs GFA + mapping |
| `untip_relative.py` | Remove short tips relative to neighbors | Topological sort, depth-based tip identification |
| `chop_misassemblies.py` | Break nodes at misassembly points | Read alignment clustering, cut at coverage-supported breakpoints |
| `unroll_tip_loops.py` | Unroll loops at tips | Coverage analysis, loop detection at tip nodes |

---

## 3. What Dinara Already Has (Existing BRG Infrastructure)

From `mode3-BrgAnchor.hpp` / `.cpp`:
- ✅ `OrientedBrgAnchor` = `(BrgAnchorId, Strand)` — equivalent to Verkko's `">nodeId"` / `"<nodeId"`
- ✅ `reverse(o)` — flip strand (equivalent to `revnode()`)
- ✅ `canon(from, to)` — canonical edge pair (equivalent to `gf.canon()`)
- ✅ `revCompPath(path)` — reverse path + flip strands (equivalent to `gf.rc_path()`)
- ✅ `BrgVectorWithDirection<T>` — container indexed by `(anchorId, strand)` (equivalent to `edges[">X"]` / `edges["<X"]`)
- ✅ `BrgEdge` — `(from, to OrientedBrgAnchor + coverage)`
- ✅ Precomputed bidirectional edge table with canonical accumulation
- ✅ Per-ReadId journeys with `BrgJourneyEntry` = (BrgAnchorId, strand)
- ✅ HTTP handlers for visualization

**What Dinara does NOT have yet:**
- ❌ Triplet counting / validation (`get_valid_triplets`)
- ❌ Node resolution (creating edge nodes, rewiring)
- ❌ Unitigification (merging linear chains)
- ❌ Path update after resolution (the `paths_crossing` index + `replace_path_nodes`)
- ❌ Iterative multi-round resolution
- ❌ Bubble popping / tip removal / unique estimation
- ❌ GFA I/O for the resolved graph

---

## 4. Implementation Plan: New Resolution Mode

### Design Principle
- **New mode**, activated separately from existing BRG anchors
- Reuse existing `OrientedBrgAnchor`, `reverse()`, `canon()`, `revCompPath()`, and `BrgVectorWithDirection<T>`
- All new code uses the naming convention **without** the "Brg" prefix (as agreed)
- Backward compatibility: existing BRG code untouched

### Phase 1: Directed Anchor Graph

**New file:** `src/mode3-AnchorGraph.hpp` + `.cpp`

```
class AnchorGraph {
    // ===== Core data structures (mirroring Verkko's) =====

    // Node sequences — each node is an oriented chain of base anchors
    // Equivalent to Verkko's node_seqs[name] = (list_of_oriented_nodes, left_clip, right_clip)
    struct NodeSequence {
        vector<OrientedBrgAnchor> anchors;  // the oriented anchor chain
        uint64_t leftClip = 0;              // bp trimmed from left
        uint64_t rightClip = 0;             // bp trimmed from right
    };
    map<NodeId, NodeSequence> nodeSeqs;

    // Directed adjacency — equivalent to Verkko's edges[">name"] = set()
    // Uses OrientedNodeId = (NodeId, strand) to represent ">X" / "<X"
    using OrientedNodeId = pair<NodeId, Strand>;
    map<OrientedNodeId, set<OrientedNodeId>> edges;

    // Edge overlaps — keyed by canonical pair
    map<pair<OrientedNodeId, OrientedNodeId>, uint64_t> edgeOverlaps;

    // Node lengths in bp
    map<NodeId, uint64_t> nodeLens;

    // Node coverages
    map<NodeId, double> nodeCoverages;

    // ===== Path index (equivalent to Verkko's paths_crossing) =====
    // For each node name, the set of paths that cross it
    // paths_crossing[nodeName] = {pathId -> path}
    // In C++: for each node, store indices into the path vector
    vector<vector<OrientedNodeId>> paths;  // all paths
    map<NodeId, vector<uint64_t>> pathsCrossing;  // node → path indices

    // ===== Resolution state =====
    uint64_t nextNodeId = 0;
};
```

**Key methods:**
```
// Graph I/O
void buildFromBrgAnchors(const BrgAnchors&);  // Convert BRG anchors to directed graph
void writeGfa(const string& filename) const;
void readGfa(const string& filename);

// Fundamental operations (direct translations of Verkko functions)
static OrientedNodeId revnode(OrientedNodeId);
static pair<OrientedNodeId, OrientedNodeId> canonEdge(OrientedNodeId, OrientedNodeId);
void addEdge(OrientedNodeId from, OrientedNodeId to, uint64_t overlap);
void removeNode(NodeId);

// Path management
void addPath(const vector<OrientedNodeId>& path);
void removePath(uint64_t pathIndex);
void rebuildPathIndex();
```

### Phase 2: Paths Crossing Index

The `paths_crossing` structure is **central** to Verkko's resolution. It allows O(1) lookup of all paths that traverse a given node — essential for counting triplets.

**In C++ this becomes:**

```cpp
// For each node, store the set of path indices whose path list mentions this node.
// When resolving, iterate these to count (predecessor, node, successor) triplets.
class PathsCrossingIndex {
    // node_name → set<path_index>
    unordered_map<NodeId, unordered_set<uint64_t>> index;

    void addPath(uint64_t pathIdx, const vector<OrientedNodeId>& path);
    void removePath(uint64_t pathIdx, const vector<OrientedNodeId>& path);

    // Iterate all paths crossing a node
    const unordered_set<uint64_t>& getPathsCrossing(NodeId) const;
};
```

### Phase 3: Resolution Engine (Core — Mirrors `resolve_triplets_kmerify.py`)

**New file:** `src/mode3-AnchorGraphResolution.cpp`

#### 3a. `getValidTriplets(NodeId node)`
- For each path crossing `node`, extract the (predecessor, node, successor) triplet
- Count occurrences of each unique triplet
- Filter by `minEdgeSupport`
- Verify all in/out neighbors are covered
- Return list of validated `(left, through, right)` triplets

#### 3b. `resolveHairpins(lengthLevel, candidateNodes)`
- Detect palindromic nodes: `">X"` has edge to `"<X"`
- For each hairpin, create fw/bw copies, rewire edges, rewrite paths
- Return set of new node names and resolved nodes

#### 3c. `resolveNodes(lengthLevel, candidateNodes)`
- For each candidate, call `getValidTriplets()`
- Filter unresolvable: hairpin neighbors at borders
- For each valid `(from→through→to)` triplet:
  - Create `edge_node` = new node representing the resolution
  - Set node sequence by extending through-node's anchor chain
  - Wire edges: `leftKey → rightKey` + RC complement
- `replacePathNodes()` — update all affected paths
- Remove resolved nodes
- Split paths at broken edges

#### 3d. `resolve(minEdgeSupport, minCoverage)`
- Min-heap of nodes by length
- Pop shortest node group
- Try `resolveHairpins()` first, then `resolveNodes()`
- After each successful resolution: `unitigifyOne()` on new nodes → push back
- Stop at `maxResolveLength`

#### 3e. `unitigifyOne(NodeId)` / `unitigifyAll()`
- `extendForward(">node")` — walk while in-degree=1 and out-degree=1
- `extendForward("<node")` — backward walk
- Merge chain via `replaceUnitig()`

#### 3f. `replaceUnitig(chain)`
- Create single new node from chain
- Merge anchor sequences with overlaps
- Rewire edges from chain endpoints to new node
- Update all crossing paths

### Phase 4: Graph Simplification

These are secondary operations that Verkko applies between resolution rounds:

#### 4a. Bubble Popping (`popBubbles()`)
- Superbubble detection (Onodera 2013 algorithm)
- Coverage-based chain merging via union-find
- Remove allele branches with coverage ratio outside 0.5x–2.5x of component average
- Max bubble pop size = 10 nodes, max node size = 200Kbp

#### 4b. Unique Node Estimation (`estimateUnique()`)
- Local superbubble analysis
- Tip analysis (`isTipUniqueCandidate`)
- Union-find grouping
- Output set of "unique" (non-repetitive) node names

#### 4c. Connect Uniques (`connectUniques()`)
- Union-find on graph edges
- Identify resolvable tangle ends from path evidence
- For each path bridging between unique nodes through a tangle: create node copies
- Handle forbidden connections

#### 4d. Forbid Unbridged Tangles (`forbidUnbridgedTangles()`)
- Compute edge coverage from paths
- Identify solid nodes (>200Kbp)
- Use union-find for tangle detection
- Check ONT and HiFi coverage thresholds
- Output forbidden connections to prevent false resolution

#### 4e. Merge Unresolved (`mergeUnresolved()`)
- Group resolved copies of same base node
- If copies share the same in-neighbor set → merge back (undo failed resolution)

#### 4f. Tip Removal (`untipRelative()`)
- Topological sort of strongly-connected components
- Compute node depths
- Remove tips shorter than threshold, relative to neighboring depth

### Phase 5: Orchestration — Full Pipeline

**New file:** `src/mode3-AnchorGraphPipeline.cpp`

```cpp
void AnchorGraph::runResolutionPipeline(const ResolutionConfig& config) {
    // ===== Stage 1: Build from BRG =====
    buildFromBrgAnchors(brgAnchors);       // Convert existing BRG to directed graph
    buildPathsFromJourneys(brgAnchors);    // Convert BRG journeys to directed paths
    computeNodeCoverages();
    unitigifyAll();                         // Initial compaction

    // ===== Stage 2: HiFi-path resolution (equivalent to 2-processGraph) =====
    for (uint64_t minSupport : config.resolveSteps) {   // e.g., {20, 10, 5}
        resolve(minSupport, config.minAllowedCoverage);
    }
    unitigifyAll();

    // ===== Stage 3: Bubble popping =====
    popBubbles(config.haploid);
    unitigifyAll();

    // ===== Stage 4: Unique estimation + tangle connection =====
    auto uniqueNodes = estimateUnique(config.longNodeThreshold, config.solidEdgeThreshold);
    auto forbiddenConnections = forbidUnbridgedTangles(uniqueNodes);
    connectUniques(uniqueNodes, forbiddenConnections);

    // ===== Stage 5: Merge back unresolvable =====
    mergeUnresolved();
    unitigifyAll();

    // ===== Stage 6: Tip removal =====
    untipRelative(config.maxRemovableLen, config.minSafeLen, config.tipFraction);

    // ===== Stage 7: Output =====
    writeGfa(config.outputGfa);
    writePaths(config.outputPaths);
    writeCoverages(config.outputCoverage);
}
```

**Configuration struct:**
```cpp
struct ResolutionConfig {
    vector<uint64_t> resolveSteps = {20, 10, 5};  // from Verkko's pop_resolve_steps
    double minAllowedCoverage = 5.0;               // from pop_min_allowed_cov
    uint64_t maxResolveLength = 500000;
    bool haploid = false;

    // Unique estimation
    uint64_t longNodeThreshold = 200000;
    uint64_t solidEdgeThreshold = 3;
    double pathConsistencyThreshold = 0.5;

    // Tip removal
    uint64_t maxRemovableLen = 200000;
    uint64_t minSafeLen = 500000;
    double tipFraction = 0.1;

    // Output
    string outputGfa;
    string outputPaths;
    string outputCoverage;
};
```

---

## 5. Implementation Order

| Priority | Component | Files | Depends On | Complexity |
|----------|-----------|-------|-----------|------------|
| **P0** | `AnchorGraph` data structures | `mode3-AnchorGraph.hpp` | Existing BRG types | Medium |
| **P0** | Build from BRG + path threading | `mode3-AnchorGraph.cpp` | P0 data structures | Medium |
| **P0** | GFA I/O | `mode3-AnchorGraph.cpp` | P0 data structures | Low |
| **P1** | `PathsCrossingIndex` | `mode3-AnchorGraph.hpp` | P0 | Low |
| **P1** | `getValidTriplets()` | `mode3-AnchorGraphResolution.cpp` | P1 | Medium |
| **P1** | `resolveNodes()` + `resolveHairpins()` | `mode3-AnchorGraphResolution.cpp` | P1 | High |
| **P1** | `unitigifyOne()` + `unitigifyAll()` | `mode3-AnchorGraphResolution.cpp` | P0 | Medium |
| **P1** | `resolve()` (main loop) | `mode3-AnchorGraphResolution.cpp` | P1 | Medium |
| **P2** | `popBubbles()` | `mode3-AnchorGraphSimplify.cpp` | P0 | High |
| **P2** | `estimateUnique()` | `mode3-AnchorGraphSimplify.cpp` | P0 | Medium |
| **P2** | `connectUniques()` | `mode3-AnchorGraphSimplify.cpp` | P2 unique | High |
| **P2** | `forbidUnbridgedTangles()` | `mode3-AnchorGraphSimplify.cpp` | P2 unique | Medium |
| **P2** | `mergeUnresolved()` | `mode3-AnchorGraphSimplify.cpp` | P1 resolve | Medium |
| **P3** | `untipRelative()` | `mode3-AnchorGraphSimplify.cpp` | P0 | Medium |
| **P3** | Full pipeline orchestration | `mode3-AnchorGraphPipeline.cpp` | All above | Low |
| **P3** | HTTP handlers for new mode | `AssemblerHttpServer-AnchorGraph.cpp` | P0 | Low |
| **P3** | Integration into main.cpp | `srcMain/main.cpp` | All above | Low |

---

## 6. Mapping: Verkko Concepts → Dinara Types

| Verkko (Python) | MBG (C++) | Dinara (C++) |
|----------------|-----------|-------------|
| `">nodeName"` / `"<nodeName"` | `pair<size_t, bool>` | `OrientedBrgAnchor` (anchorId + strand) |
| `gf.revnode(x)` | – | `reverse(OrientedBrgAnchor)` |
| `gf.canon(from, to)` | – | `canon(OrientedBrgAnchor, OrientedBrgAnchor)` |
| `gf.rc_path(path)` | – | `revCompPath(vector<OrientedBrgAnchor>)` |
| `edges[">X"]` = `set(...)` | `VectorWithDirection<T>` | `BrgVectorWithDirection<T>` or `map<OrientedNodeId, set<OrientedNodeId>>` |
| `paths_crossing[name]` | `readsCrossingNode` | `PathsCrossingIndex` |
| `node_seqs[name]` | `UnitigGraph::unitigs` | `AnchorGraph::NodeSequence` |
| `edge_overlaps[canon(a,b)]` | `overlaps` map | `AnchorGraph::edgeOverlaps` |
| `resolve_triplets_kmerify.py` | `UnitigResolver::resolve()` | `AnchorGraph::resolve()` |
| `unitigify_one()` | `unitigifyOne()` | `AnchorGraph::unitigifyOne()` |
| `get_valid_triplets()` | `getValidTriplets()` | `AnchorGraph::getValidTriplets()` |
| `pop_bubbles_coverage_based.py` | – | `AnchorGraph::popBubbles()` |
| GAF path format | `ReadPath` | `vector<OrientedNodeId>` (per read) |

---

## 7. Key Design Decisions

### 7.1 Node Identity
Verkko uses string names (`"utig4-123"`). MBG uses integer node IDs. We should use **integer NodeId** (uint64_t) for performance, with a name lookup table for GFA output.

### 7.2 Edge Storage
Verkko uses `dict → set`. For C++ we have two options:
- **Option A:** `unordered_map<OrientedNodeId, unordered_set<OrientedNodeId>>` — most faithful to Verkko, easy to add/remove
- **Option B:** `BrgVectorWithDirection<vector<NodeId>>` — denser, but harder to do dynamic add/remove during resolution

**Recommendation:** Option A for the resolution phase (dynamic). Convert to Option B for read-only query after resolution stabilizes.

### 7.3 Path Storage
Verkko stores raw Python lists and uses `id(path)` as identity. In C++:
- Store paths in a `vector<vector<OrientedNodeId>>` with a stable index
- The `PathsCrossingIndex` maps `NodeId → set<pathIndex>`
- Mark deleted paths with a tombstone rather than removing (to keep indices stable during resolution)

### 7.4 Resolution Granularity
Verkko resolves at the MBG k-mer level (small nodes). Dinara anchors are larger. The resolution loop works the same way — process nodes from smallest to largest — but the effective resolution will be coarser. This is acceptable for a first implementation; we can add sub-anchor resolution later.

### 7.5 Bidirectional Invariant
Every operation that adds an edge `A → B` must also add `revnode(B) → revnode(A)`. This is Verkko's fundamental invariant. Our `addEdge()` method must enforce this atomically.

---

## 8. Testing Strategy

1. **Unit tests** for graph primitives: `revnode`, `canon`, `addEdge`/`removeNode` bidir invariant
2. **Triplet counting** on hand-constructed small graphs with known paths
3. **Unitigification** on linear chain graphs — verify merging
4. **Full resolution** on small synthetic graphs with known ground truth
5. **Regression** against existing BRG assembly — ensure new mode doesn't break old mode
6. **End-to-end** on real data — compare resolution before/after

---

## 9. Risks and Mitigations

| Risk | Mitigation |
|------|-----------|
| Anchor granularity too coarse for effective resolution | Can add sub-anchor nodes by splitting at divergence points |
| Path evidence too sparse for triplet counting | Lower `minEdgeSupport` thresholds; combine with coverage signals |
| Performance of `map`/`set` during resolution | Profile first; switch to sorted vectors or flat hash maps if needed |
| Bubble popping too aggressive on diploid regions | Copy Verkko's coverage ratio bounds (0.5x–2.5x) directly |
| Existing BRG code interference | Completely separate class hierarchy; shared only via `OrientedBrgAnchor` type |

---

## 10. Summary

The plan implements Verkko's resolution procedure in Dinara through 5 phases:

1. **P0 (Foundation):** `AnchorGraph` data structures + build from BRG + GFA I/O
2. **P1 (Core Resolution):** Triplet counting, node resolution, hairpin resolution, unitigification, iterative resolve loop
3. **P2 (Simplification):** Bubble popping, unique estimation, tangle connection, unbridged tangle protection, merge unresolved
4. **P3 (Polish):** Tip removal, full pipeline orchestration, HTTP handlers, integration

The approach preserves full backward compatibility with existing BRG anchors while adding a principled, Verkko-proven resolution pipeline.

---

## APPENDIX B: Gaps Detected — Verkko Features Missing From Current Implementation

After deep re-analysis of Verkko's codebase, the following features/details are missing or incomplete in our current `DirectedAnchorGraph` implementation.

### B.1 Missing: Edge Overlaps

**Verkko:** Tracks `edge_overlaps[canon(from, to)] = int` for every edge. This overlap is used throughout:
- `get_unitig_len_path()` subtracts overlaps when computing total length
- `replace_unitig()` uses overlaps to correctly merge sequences, determining how much to clip from consecutive node sequences
- `resolve_nodes()` copies overlaps to newly wired edges
- `resolve_hairpins()` copies overlaps to hairpin copy edges

**Dinara (current):** `DagEdgeSet` / `DagEdgeMap` have no overlap field. Edges are just presence/absence. Node lengths are stored in `DagNodeInfo::lengthBp` but edge overlaps are not tracked.

**Impact:** Without edge overlaps, unitig length calculations are overestimates, and sequence merging during unitigification will produce incorrect assemblies.

**Fix:** Add `std::unordered_map<CanonicalEdge, uint64_t> edgeOverlaps;` to `DirectedAnchorGraph`.

### B.2 Missing: Node Sequence Hierarchy (node_seqs with clips)

**Verkko:** `node_seqs[name] = (sequence_of_oriented_nodes, left_clip, right_clip)` — each node tracks:
1. An ordered list of sub-nodes (kmer nodes from MBG)
2. Integer left_clip (how many bp clipped from first sub-node)
3. Integer right_clip (how many bp clipped from last sub-node)

This is critical for `replace_unitig()` where merged chains must correctly handle partial overlaps.

**Dinara (current):** `DagNodeInfo::anchorChain` is `vector<OrientedBrgAnchor>` but has no clip tracking (`leftClip`, `rightClip` fields).

**Fix:** Add `uint64_t leftClip = 0; uint64_t rightClip = 0;` to `DagNodeInfo`.

### B.3 Missing: remove_and_split_low_coverage() Pre-Resolution Step

**Verkko:** Before resolution, runs a cleanup step that:
1. Computes edge coverage from paths
2. Finds "safe" unitigs (high-coverage, unambiguous) via `get_safe_unitigs_and_edges()`
3. Removes low-coverage edges not in safe unitigs
4. Removes low-coverage nodes not in safe unitigs
5. Splits paths at removed nodes/edges

**Dinara (current):** `runResolution()` goes straight to `unitigifyAll()` → `resolveRound()` with no pre-cleanup.

**Impact:** Resolution will attempt to resolve using noisy, low-coverage edges that Verkko would have already removed.

**Fix:** Implement `removeAndSplitLowCoverage()` method, call it before `unitigifyAll()` in `runResolution()`.

### B.4 Missing: Removable Nodes Set

**Verkko:** Marks nodes with coverage < `min_allowed_coverage` as "removable". In `get_valid_triplets()`, when checking neighbor coverage, removable nodes can be skipped — they don't block resolution even if no triplet covers them.

**Dinara (current):** No concept of removable nodes. All neighbors must be covered by triplets or resolution fails.

**Impact:** Resolution will be too conservative. Many nodes will fail to resolve because a low-coverage garbage neighbor isn't covered by any triplet.

**Fix:** Add `std::unordered_set<DagNodeId> removableNodes;` and modify `getValidTriplets()` to use it.

### B.5 Missing: longest_extension_per_node Computation

**Verkko:** In `resolve_nodes()`, computes how far edge nodes should extend into their neighbors:
```python
longest_extension_per_node[node] = min(1000, neighbor_lengths)
```
This is capped at 1000bp to prevent strange multi-megabase bubble artifacts. Edge nodes get EXTENDED sequences from their neighbors, not just the through-node's sequence.

**Dinara (current):** Edge nodes in `resolveNodes()` inherit the through-node's anchor chain but don't extend into neighbors. No extension length computation.

**Impact:** Edge nodes won't carry enough sequence context, potentially causing issues in subsequent resolution rounds.

**Fix:** Add extension computation and sequence extension logic to `resolveNodes()`.

### B.6 Missing: Two-Tier Filtering (maybe_resolvable → resolvable)

**Verkko:** Uses two sets:
- `maybe_resolvable`: All nodes that _could_ be resolved
- `resolvable`: Subset that passes iterative filtering (no single-bp neighbors, no hairpin neighbors)

The filtering loop iterates until convergence, removing nodes whose resolution would create problems.

**Dinara (current):** `resolveNodes()` has a single-pass filter but no iterative convergence loop.

**Fix:** Add iterative filtering in `resolveNodes()` that removes nodes bordering single-bp or hairpin neighbors.

### B.7 Missing: Proper Edge Node Naming Convention

**Verkko:** Edge nodes are named `"edge_FROM_TO"` using oriented node names (e.g., `"edge_>node1_>node2"`). This naming enables debugging and tracing.

**Dinara (current):** Edge nodes get generic `nextNodeId++` identifiers with no semantic naming.

**Impact:** Debugging difficulty only; no algorithmic impact.

### B.8 Missing: Path Splitting at Broken Edges

**Verkko:** After any path rewriting (replace_path_nodes, replace_unitig, hairpin resolution), explicitly checks for "broken" edges — consecutive path nodes that lack an edge — and splits paths at those points.

**Dinara (current):** `splitPathsAtBreaks()` exists but may not be called consistently after every operation.

**Fix:** Audit all path-modifying operations to ensure `splitPathsAtBreaks()` is called after each.

### B.9 Missing: Handling of Circular Unitigs

**Verkko:** `unitigify_one()` explicitly handles circular paths (when `extend_forward` returns to the start node). Removes duplicate endpoint before calling `replace_unitig()`.

**Dinara (current):** `unitigifyOne()` doesn't check for circular chains.

**Fix:** Add circular detection in `extendForward()` and `unitigifyOne()`.

### B.10 Missing: Safe Unitig Protection During Low-Coverage Removal

**Verkko:** `add_safe_unitig()` walks from high-coverage nodes along unambiguous chains, marking all nodes and edges in the chain as "safe". These are never removed even if their neighbors have low coverage.

**Dinara (current):** No concept of safe unitigs.

**Fix:** Implement as part of `removeAndSplitLowCoverage()`.

### B.11 Missing: Phase 2 Supporting Scripts (Not Yet Implemented)

These Verkko scripts have NO equivalent in Dinara yet:

| Script | Function | Priority |
|--------|----------|----------|
| `pop_bubbles_coverage_based.py` | Superbubble detection + coverage-based bubble popping | High |
| `estimate_unique_local.py` | Identify "unique" (non-repeat) nodes via coverage+paths | High |
| `connect_uniques.py` | Bridge tangles between unique nodes using read paths | High |
| `forbid_unbridged_tangles.py` | Prevent resolution of tangles without bridge evidence | Medium |
| `merge_unresolved_dbg_nodes.py` | Undo failed resolutions (merge back node copies) | Medium |
| `untip_relative.py` | Relative tip removal using topological sort + depth | Medium |
| `unroll_tip_loops.py` | Unroll loops at graph tips | Low |
| `chop_misassemblies.py` | Break nodes at coverage drops | Low |

---

## APPENDIX C: Complete Verkko 2-processGraph Snakefile Pipeline

The exact sequence of operations in Verkko's `2-processGraph.sm` (HiFi resolution):

```
Step 1: resolve_triplets_kmerify.py
  Input:  unitig-unrolled.gfa + .gaf + .csv (coverages)
  Params: resolve_steps="20 10 5", min_allowed_coverage=auto
  Output: unitig-unrolled-hifi-resolved.gfa + .gaf + .noseq.gfa

Step 2: pop_bubbles_coverage_based.py
  Input:  hifi-resolved.gfa + .noseq.gfa
  Output: unitig-unrolled-popped.gfa (nodes in unique chains)

Step 3: estimate_unique_local.py
  Input:  popped.gfa + alignments + nodecov.csv
  Output: unique-nodes-list.txt

Step 4: connect_uniques.py
  Input:  popped.gfa + unique-nodes + paths + forbidden-connections
  Output: connected.gfa + .noseq.gfa

Step 5: forbid_unbridged_tangles.py
  Input:  connected.gfa + unique-nodes + paths + ONT alns + HiFi alns
  Output: forbidden-connections.txt (fed back to step 4 iteratively)

Step 6: merge_unresolved_dbg_nodes.py
  Input:  connected.gfa
  Output: merged.gfa (undo unnecessary resolution copies)

Step 7: unitigify.py
  Input:  merged.gfa
  Output: unitig-normal-connected.gfa (final compaction)
```

The `4-processONT.sm` runs the SAME scripts but with ONT paths and different coverage thresholds.

The `5-untip.sm` adds:
```
Step 8: untip_relative.py (multiple rounds with different params)
Step 9: chop_misassemblies.py
Step 10: unroll_tip_loops.py
Step 11: pop_bubbles_coverage_based.py (second round)
Step 12: Final unitigify
```