#pragma once

#include <cstdint>
#include <memory>

namespace dinara {
    class Assembler;
    class AssemblerOptions;
    namespace mode3 {
        class Anchors;
    }
}

namespace dinara {
    // Helper function to create shasta2 anchors efficiently (multithreaded)
    void createShasta2Anchors(
        dinara::Assembler& assembler,
        const dinara::AssemblerOptions& options,
        uint64_t threadCount
    );

    // Same as above, but uses an already built Dinara anchor set verbatim for conversion.
    void createShasta2Anchors(
        dinara::Assembler& assembler,
        const dinara::AssemblerOptions& options,
        uint64_t threadCount,
        const std::shared_ptr<const dinara::mode3::Anchors>& precomputedDinaraAnchors
    );
}
