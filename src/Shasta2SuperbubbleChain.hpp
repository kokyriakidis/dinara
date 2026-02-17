#pragma once

#include "Shasta2Superbubble.hpp"

namespace dinara {
    class Shasta2SuperbubbleChain;
}

class dinara::Shasta2SuperbubbleChain : public vector<Shasta2Superbubble> {
public:
    uint64_t phase1(
        Shasta2AssemblyGraph& assemblyGraph,
        uint64_t superbubbleChainId) const;
};
