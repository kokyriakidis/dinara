# Snarl Finder

Vendored from [vg](https://github.com/vgteam/vg) with protobuf, SnarlManager,
bdsg overlay, and sonLib dependencies removed. Only `traverse_decomposition`
is exposed.

## Setup

Before building, clone the two CMake dependencies into this directory:

```sh
cd external/snarls
git clone --depth 1 https://github.com/vgteam/libhandlegraph.git
git clone --depth 1 https://github.com/vgteam/structures.git
```

These are gitignored and not committed.

## What's vendored

- `three_edge_connected_components.hpp/.cpp` — Tsin's algorithm (pure C++, sonLib `_cactus` variant removed)
- `integrated_snarl_finder.hpp/.cpp` — cactus graph + snarl decomposition (takes `RankedHandleGraph` directly)
- `snarl_finder_base.hpp` — minimal base class replacing vg's `snarls.hpp`
- `dinara_handle_graph.hpp` — adapter wrapping Dinara graphs as `RankedHandleGraph`
