# Tutorial: Integrating WFA2-lib into Shasta

This document details the steps taken to integrate the WFA2-lib (Wavefront Alignment Algorithm) into the Shasta assembler. It covers adding the library as a submodule, configuring CMake, and modifying the C++ source code to use it as an alternative aligner.

## Prerequisites

*   **Shasta Build Environment**: Standard dependencies for building Shasta (CMake, GCC/Clang, etc.).
*   **Git**: To clone the WFA2-lib repository.

---

## Step 1: Setting up the External Library

We use the `WFA2-lib` library, which provides a fast and memory-efficient implementation of the wavefront alignment algorithm.

### 1.1 Clone the Repository
We cloned the repository into `shasta/external/`:
```bash
cd shasta/external
git clone https://github.com/smarco/WFA2-lib
```

Unlike A*PA2, we do not need to manually build a static library beforehand. We can integrate it directly into the CMake build system, which will handle the compilation.

---

## Step 2: Integrating into CMake

We need to tell Shasta's build system how to build WFA2-lib and where to find its headers and libraries.

### 2.1 Update Root `CMakeLists.txt`
We added WFA2-lib as a subdirectory so CMake knows to build it:
```cmake
add_subdirectory(external/WFA2-lib wfa2lib)
```

### 2.2 Update Component `CMakeLists.txt`
We modified `staticLibrary/CMakeLists.txt`, `dynamicLibrary/CMakeLists.txt`, and `staticExecutable/CMakeLists.txt`.

**1. Include Directories:**
Add the path to the WFA2-lib headers and its C++ bindings:
```cmake
include_directories(../external/WFA2-lib)
include_directories(../external/WFA2-lib/bindings/cpp)
```

**2. Link Libraries:**
Link the static libraries provided by WFA2-lib (`wfa2cpp_static` for C++ bindings and `wfa2_static` for the core library):
```cmake
target_link_libraries(
    ...
    wfa2cpp_static
    wfa2_static
)
```

---

## Step 3: C++ Code Integration

Now we modify Shasta to use the library.

### 3.1 Add Configuration Option
In `src/AssemblerOptions.hpp` and `src/AssemblerOptions.cpp`:
*   Added `bool useWfa;` to the `AlignOptions` struct.
*   Registered it as a command-line option `--Align.useWfa`.

### 3.2 Update `ProjectedAlignment` Class
This is where the alignment actually happens.

**1. Include Header:**
In `src/ProjectedAlignment.cpp`, we include the WFA2 C++ bindings header:
```cpp
#include "WFAligner.hpp"
```
*Note: We use the C++ bindings for easier integration.*

**2. Pass the Option:**
We updated the `ProjectedAlignment` constructor to accept `bool useWfa` and store it.

**3. Implement the Logic:**
In `ProjectedAlignmentSegment::computeAlignment`, we added a branch to use WFA2 when the flag is set:

```cpp
if (useWfa) {
    // 1. Initialize WFA2 Aligner
    // We use WFAlignerGapAffine for affine gap penalties
    wfa::WFAlignerGapAffine aligner(
        mismatchScore,
        gapScore,
        gapScore,
        wfa::WFAligner::Alignment,
        wfa::WFAligner::MemoryHigh);

    // 2. Perform Alignment
    aligner.alignEnd2End(
        (const char*)sequence0.data(), sequence0.size(),
        (const char*)sequence1.data(), sequence1.size());

    editDistance = aligner.getAlignmentScore();

    // 3. Convert CIGAR to Shasta Format
    // WFA2 returns a CIGAR string (e.g., "10M2I5M").
    string cigar = aligner.getCIGAR(true); // true = show mismatches (X/= instead of M)
    
    alignment.clear();
    // ... parsing logic to convert CIGAR to vector<pair<bool, bool>> ...

} else if (useAstarPa) {
    // A*PA2 logic
} else {
    // Original SeqAn logic
}
```

### 3.3 Wire it Up
Finally, in `src/AssemblerVariantClustering.cpp` (and other call sites), we pass the option from `AlignOptions` to the `ProjectedAlignment` constructor:

```cpp
const ProjectedAlignment projectedAlignment(
    ...
    data.alignOptions->useAstarPa,
    data.alignOptions->useWfa // Pass the new flag
);
```

---

## Summary

1.  **Why WFA2-lib?** It offers a modern, efficient alignment algorithm (Wavefront) that can be faster than traditional dynamic programming methods like Needleman-Wunsch, especially for similar sequences.
2.  **Why C++ Bindings?** WFA2-lib provides convenient C++ wrappers (`WFAligner`) that handle memory management and object-oriented usage, simplifying the integration compared to raw C calls.
3.  **Why CMake `add_subdirectory`?** This allows WFA2-lib to be built as part of the Shasta build process, ensuring compiler compatibility and simplifying dependency management.
