#pragma once

#include "cstdint.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace dinara {
    class Base;
    class AlignedBase;

    void shasta2Abpoa(
        const vector< vector<Base> >& sequences,
        vector< pair<Base, uint64_t> >& consensus,
        vector< vector<AlignedBase> >& alignment,
        vector<AlignedBase>& alignedConsensus,
        bool computeAlignment);

    void shasta2Abpoa(
        const vector< pair<vector<Base>, uint64_t> >& sequencesWithWeights,
        vector< pair<Base, uint64_t> >& consensus,
        vector< vector<AlignedBase> >& alignment,
        vector<AlignedBase>& alignedConsensus);
}
