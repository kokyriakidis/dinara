# Dinara long read assembler

De novo assembler for long reads, optimized for Oxford Nanopore (ONT) reads.

For quick start information see [here](https://kokyriakidis.github.io/dinara/QuickStart.html).

The complete user documentation is available [here](https://kokyriakidis.github.io/dinara/).

Please file GitHub issues to report problems, request help, or ask questions.
Please keep each issue on a single topic when possible. 


## Dependencies

Dinara requires the following dependencies to be installed on your system:

*   **astar-pairwise-aligner (C bindings)**: You must install the `astar-pairwise-aligner` library (specifically `libastarpa_c` and `astarpa.h`).
    *   The C bindings are provided by the `astarpa-c` crate in the [astar-pairwise-aligner repository](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner).
    *   The header `astarpa.h` must be in a system include path (e.g., `/usr/include`, `/usr/local/include`).
    *   The library `libastarpa_c` must be in a system library path (e.g., `/usr/lib`, `/usr/local/lib`).
    *   **Note**: If the library is in a non-standard location, you may need to set `LD_LIBRARY_PATH` (e.g., `export LD_LIBRARY_PATH=/path/to/astarpa/target/release`) for the application to run.
