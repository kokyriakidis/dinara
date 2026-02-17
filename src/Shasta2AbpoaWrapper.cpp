#include "Shasta2AbpoaWrapper.hpp"

#include "Base.hpp"
#include "globalMsa.hpp"

#include <algorithm>
#include <array>

using namespace dinara;
using namespace std;

namespace {

void computeConsensusFromAlignment(
    const vector<uint64_t>& rowWeights,
    const vector< vector<AlignedBase> >& alignment,
    vector< pair<Base, uint64_t> >& consensus,
    vector<AlignedBase>& alignedConsensus)
{
    const uint64_t n = alignment.size();
    if(n == 0) {
        consensus.clear();
        alignedConsensus.clear();
        return;
    }

    const uint64_t alignmentLength = alignment.front().size();
    for(const auto& row: alignment) {
        if(row.size() != alignmentLength) {
            consensus.clear();
            alignedConsensus.clear();
            return;
        }
    }

    consensus.clear();
    alignedConsensus.resize(alignmentLength);

    for(uint64_t i=0; i<alignmentLength; i++) {
        array<uint64_t, 5> baseCoverage;
        fill(baseCoverage.begin(), baseCoverage.end(), 0);

        for(uint64_t j=0; j<n; j++) {
            const AlignedBase alignedBase = alignment[j][i];
            const uint64_t weight = (j < rowWeights.size()) ? rowWeights[j] : 1;
            baseCoverage[alignedBase.value] += weight;
        }

        const auto it = max_element(baseCoverage.begin(), baseCoverage.end());
        const AlignedBase consensusBase = AlignedBase::fromInteger(uint64_t(it - baseCoverage.begin()));
        const uint64_t coverage = *it;

        alignedConsensus[i] = consensusBase;
        if(!consensusBase.isGap()) {
            consensus.push_back(make_pair(Base(consensusBase), coverage));
        }
    }
}

} // namespace

void dinara::shasta2Abpoa(
    const vector< vector<Base> >& sequences,
    vector< pair<Base, uint64_t> >& consensus,
    vector< vector<AlignedBase> >& alignment,
    vector<AlignedBase>& alignedConsensus,
    const bool computeAlignment)
{
    vector< pair<vector<Base>, uint64_t> > sequencesWithWeights;
    sequencesWithWeights.reserve(sequences.size());
    for(const vector<Base>& sequence: sequences) {
        sequencesWithWeights.push_back(make_pair(sequence, 1));
    }

    globalMsaSpoa(sequencesWithWeights, alignment);

    vector<uint64_t> rowWeights(sequences.size(), 1);
    computeConsensusFromAlignment(rowWeights, alignment, consensus, alignedConsensus);

    if(!computeAlignment) {
        alignment.clear();
        alignedConsensus.clear();
    }
}

void dinara::shasta2Abpoa(
    const vector< pair<vector<Base>, uint64_t> >& sequencesWithWeights,
    vector< pair<Base, uint64_t> >& consensus,
    vector< vector<AlignedBase> >& alignment,
    vector<AlignedBase>& alignedConsensus)
{
    globalMsaSpoa(sequencesWithWeights, alignment);

    vector<uint64_t> rowWeights;
    rowWeights.reserve(sequencesWithWeights.size());
    for(const auto& p: sequencesWithWeights) {
        rowWeights.push_back(p.second);
    }

    computeConsensusFromAlignment(rowWeights, alignment, consensus, alignedConsensus);
}
