#pragma once

#include <cstdint>

namespace dinara {
    class Assembler;
    class AssemblerOptions;
}

namespace dinara {
    // Helper function to create shasta2 anchors efficiently (multithreaded)
    void createShasta2Anchors(
        dinara::Assembler& assembler,
        const dinara::AssemblerOptions& options,
        uint64_t threadCount
    );
}
