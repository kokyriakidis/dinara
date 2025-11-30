# Tutorial: Integrating A*PA2 into Shasta

This document details the steps taken to integrate the A*PA2 (A-star Pairwise Aligner 2) library into the Shasta assembler. It covers building the Rust library with C bindings, linking it via CMake, and modifying the C++ source code to use it.

## Prerequisites

*   **Rust Toolchain**: You need `cargo` and `rustc` installed. We used the `nightly` toolchain, but `stable` may work depending on the A*PA2 version.
*   **Shasta Build Environment**: Standard dependencies for building Shasta (CMake, GCC/Clang, etc.).

---

## Step 1: Setting up the External Library

We use the `astar-pairwise-aligner` library, which is written in Rust. To use it in C++, we need to compile its C bindings (`astarpa-c`) into a static library.

### 1.1 Clone the Repository
We cloned the repository into `shasta/external/`:
```bash
cd shasta/external
git clone https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner
```

### 1.2 Modify C Bindings for Safety (Crucial Step)
The `astarpa-c` crate exposes Rust functions to C. However, modern Rust compilers (especially nightly) are very strict about `unsafe` code. We encountered errors because `unsafe` operations were performed inside `unsafe` functions without explicit `unsafe` blocks (deprecated behavior), and `#[no_mangle]` is now considered an unsafe attribute.

We had to modify `shasta/external/astar-pairwise-aligner/astarpa-c/src/lib.rs`:

1.  **Wrap unsafe operations**: Any raw pointer dereference or call to `from_raw` must be inside an `unsafe { ... }` block, even if the function itself is marked `unsafe`.
2.  **Mark `no_mangle` as unsafe**: Change `#[no_mangle]` to `#[unsafe(no_mangle)]`.

**Example Change:**
```rust
// Before
#[no_mangle]
pub unsafe extern "C" fn astarpa_free_cigar(cigar: *mut u8) {
    drop(CString::from_raw(cigar as *mut i8))
}

// After
#[unsafe(no_mangle)]
pub unsafe extern "C" fn astarpa_free_cigar(cigar: *mut u8) {
    unsafe { drop(CString::from_raw(cigar as *mut i8)) }
}
```

### 1.3 Build the Static Library
We need to build the library so it can be linked into Shasta's shared library (`shasta.so`). This requires Position Independent Code (PIC).

Run the build command from `shasta/external/astar-pairwise-aligner/astarpa-c`:
```bash
# RUSTFLAGS="-C relocation-model=pic" is REQUIRED for linking into shared libraries
RUSTFLAGS="-C relocation-model=pic" cargo build --release
```

This produces the static library at:
`../target/release/libastarpa_c.a` (relative to the workspace root).

---

## Step 2: Integrating into CMake

We need to tell Shasta's build system where to find the header files and the library.

### 2.1 Update `CMakeLists.txt`
We modified `staticLibrary/CMakeLists.txt`, `dynamicLibrary/CMakeLists.txt`, and `staticExecutable/CMakeLists.txt`.

**1. Include Directories:**
Add the path to the C header (`astarpa.h`):
```cmake
include_directories(${CMAKE_SOURCE_DIR}/external/astar-pairwise-aligner/astarpa-c)
```

**2. Link Libraries:**
Link the static library we built. We used `${CMAKE_SOURCE_DIR}` to ensure the path is absolute and robust:
```cmake
target_link_libraries(
    ...
    ${CMAKE_SOURCE_DIR}/external/astar-pairwise-aligner/target/release/libastarpa_c.a
)
```

---

## Step 3: C++ Code Integration

Now we modify Shasta to use the library.

### 3.1 Add Configuration Option
In `src/AssemblerOptions.hpp` and `src/AssemblerOptions.cpp`:
*   Added `bool useAstarPa;` to the `AlignOptions` struct.
*   Registered it as a command-line option `--Align.useAstarPa`.

### 3.2 Update `ProjectedAlignment` Class
This is where the alignment actually happens.

**1. Include Header:**
In `src/ProjectedAlignment.cpp`, we include the C header inside an `extern "C"` block to prevent C++ name mangling issues:
```cpp
extern "C" {
#include "astarpa.h"
}
```

**2. Pass the Option:**
We updated the `ProjectedAlignment` constructor to accept `bool useAstarPa` and store it.

**3. Implement the Logic:**
In `ProjectedAlignmentSegment::computeAlignment`, we added a branch:

```cpp
if (useAstarPa) {
    // 1. Call A*PA2
    char* cigar = nullptr;
    size_t cigarLen = 0;
    int64_t cost = astarpa2_simple(
        sequence0.data(), sequence0.size(),
        sequence1.data(), sequence1.size(),
        (unsigned char**)&cigar, &cigarLen);

    // 2. Parse CIGAR string to Shasta format
    // A*PA2 returns a CIGAR string (e.g., "10M2I5M").
    // Shasta uses a vector<pair<bool, bool>> where:
    // (true, true) = Match/Mismatch
    // (false, true) = Insertion (gap in seq0)
    // (true, false) = Deletion (gap in seq1)
    
    alignment.clear();
    string cigarString(cigar, cigarLen);
    astarpa_free_cigar((unsigned char*)cigar); // Don't forget to free!

    // ... parsing logic ...

} else {
    // Original SeqAn logic
    seqanAlign(...);
}
```

### 3.3 Wire it Up
Finally, in `src/AssemblerAlign.cpp`, we pass the option from `AlignOptions` to the `ProjectedAlignment` constructor:

```cpp
const ProjectedAlignment projectedAlignment(
    ...
    data.alignOptions->useAstarPa // Pass the flag
);
```

---

## Summary of "Why"

1.  **Why `unsafe` changes?** Rust is becoming safer. Calling C functions is inherently unsafe, and the compiler now demands you explicitly acknowledge this with `unsafe` blocks, even inside functions already marked `unsafe`.
2.  **Why `-fPIC`?** Shasta builds a shared library (`shasta.so`) for its Python API. You cannot link a static library built without `-fPIC` (position-dependent) into a shared library.
3.  **Why `extern "C"`?** C++ compilers "mangle" function names (e.g., `_Z8functionv`). C compilers don't. To link C code (or Rust exposing C ABI) into C++, you must tell the C++ compiler to expect unmangled C names.
4.  **Why CIGAR parsing?** A*PA2 outputs a standard CIGAR string. Shasta uses a custom internal representation for alignments. We had to write a translator to bridge the two.
