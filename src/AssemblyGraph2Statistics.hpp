#ifndef DINARA_ASSEMBLY_GRAPH2_STATISTICS_HPP
#define DINARA_ASSEMBLY_GRAPH2_STATISTICS_HPP

#include "cstdint.hpp"

namespace dinara {
    class AssemblyGraph2Statistics;
}

class dinara::AssemblyGraph2Statistics {
public:
    uint64_t totalBubbleChainLength = 0;
    uint64_t bubbleChainN50 = 0;
    uint64_t totalDiploidLengthBothHaplotypes = 0;
    uint64_t diploidN50 = 0;
    uint64_t totalHaploidLength = 0;
    uint64_t haploidN50 = 0;
    uint64_t outsideBubbleChainsLength = 0;
    uint64_t simpleSnpBubbleTransitionCount = 0;
    uint64_t simpleSnpBubbleTransversionCount = 0;
    uint64_t nonSimpleSnpBubbleCount = 0;
};

#endif

