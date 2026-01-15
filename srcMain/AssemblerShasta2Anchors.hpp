#pragma once

#include <cstdint>

namespace dinara {
    class Assembler;
}

namespace dinara {
    // Helper function to create shasta2 anchors efficiently (multithreaded)
    void createShasta2Anchors(
        dinara::Assembler& assembler,
        uint64_t threadCount
    );
}
