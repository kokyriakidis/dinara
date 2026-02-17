#include "Shasta2PoastaWrapper.hpp"

#include "Base.hpp"
#include "DINARA_ASSERT.hpp"

#include "poasta/poasta.h"

#include <algorithm>
#include <array>

using namespace dinara;

void dinara::shasta2Poasta(
    const vector< vector<Base> >& sequences,
    vector< pair<Base, uint64_t> >& consensus,
    vector< vector<AlignedBase> >& alignment,
    vector<AlignedBase>& alignedConsensus)
{
    const int mismatchScore = 4;
    const int gapExtendScore = 2;
    const int gapOpenScore = 6;

    PoastaGraph* graph = poasta_create_graph();

    string sequenceString;
    for(const vector<Base>& sequence: sequences) {
        sequenceString.clear();
        for(const Base base: sequence) {
            sequenceString.push_back(base.character());
        }
        poasta_add_sequence(graph, sequenceString.data(), sequence.size(),
            mismatchScore, gapExtendScore, gapOpenScore);
    }

    const PoastaMsa msa = poasta_get_msa(graph);
    DINARA_ASSERT(msa.num_sequences == sequences.size());

    const uint64_t n = sequences.size();
    alignment.clear();
    for(uint64_t i=0; i<n; i++) {
        const char* sequenceCharacters = msa.sequences[i];
        vector<AlignedBase>& alignmentRow = alignment.emplace_back();
        for(uint64_t j=0; ; j++) {
            const char c = sequenceCharacters[j];
            if(c == 0) {
                break;
            }
            alignmentRow.push_back(AlignedBase::fromCharacter(c));
        }
    }

    const uint64_t alignmentLength = alignment.front().size();
    for(uint64_t i=1; i<alignment.size(); i++) {
        DINARA_ASSERT(alignment[i].size() == alignmentLength);
    }

    consensus.clear();
    alignedConsensus.resize(alignmentLength);
    for(uint64_t i=0; i<alignmentLength; i++) {
        std::array<uint64_t, 5> baseCount;
        std::fill(baseCount.begin(), baseCount.end(), 0);
        for(uint64_t j=0; j<n; j++) {
            const AlignedBase alignedBase = alignment[j][i];
            ++baseCount[alignedBase.value];
        }

        const auto it = std::max_element(baseCount.begin(), baseCount.end());
        const AlignedBase consensusBase = AlignedBase::fromInteger(uint64_t(it - baseCount.begin()));
        const uint64_t coverage = *it;

        alignedConsensus[i] = consensusBase;
        if(not consensusBase.isGap()) {
            consensus.push_back(make_pair(Base(consensusBase), coverage));
        }
    }

    poasta_free_msa(msa);
    poasta_free_graph(graph);
}

void dinara::shasta2Poasta(
    const vector< pair<vector<Base>, uint64_t> >& sequencesWithWeights,
    vector< pair<Base, uint64_t> >& consensus,
    vector< vector<AlignedBase> >& alignment,
    vector<AlignedBase>& alignedConsensus)
{
    const int mismatchScore = 4;
    const int gapExtendScore = 2;
    const int gapOpenScore = 6;

    PoastaGraph* graph = poasta_create_graph();

    string sequenceString;
    for(const auto& [sequence, weight]: sequencesWithWeights) {
        sequenceString.clear();
        for(const Base base: sequence) {
            sequenceString.push_back(base.character());
        }
        poasta_add_sequence_with_weight(graph, sequenceString.data(), sequence.size(), int(weight),
            mismatchScore, gapExtendScore, gapOpenScore);
    }

    const PoastaMsa msa = poasta_get_msa(graph);
    DINARA_ASSERT(msa.num_sequences == sequencesWithWeights.size());

    const uint64_t n = sequencesWithWeights.size();
    alignment.clear();
    for(uint64_t i=0; i<n; i++) {
        const char* sequenceCharacters = msa.sequences[i];
        vector<AlignedBase>& alignmentRow = alignment.emplace_back();
        for(uint64_t j=0; ; j++) {
            const char c = sequenceCharacters[j];
            if(c == 0) {
                break;
            }
            alignmentRow.push_back(AlignedBase::fromCharacter(c));
        }
    }

    const uint64_t alignmentLength = alignment.front().size();
    for(uint64_t i=1; i<alignment.size(); i++) {
        DINARA_ASSERT(alignment[i].size() == alignmentLength);
    }

    consensus.clear();
    alignedConsensus.resize(alignmentLength);
    for(uint64_t i=0; i<alignmentLength; i++) {
        std::array<uint64_t, 5> baseCoverage;
        std::fill(baseCoverage.begin(), baseCoverage.end(), 0);
        for(uint64_t j=0; j<n; j++) {
            const AlignedBase alignedBase = alignment[j][i];
            baseCoverage[alignedBase.value] += sequencesWithWeights[j].second;
        }

        const auto it = std::max_element(baseCoverage.begin(), baseCoverage.end());
        const AlignedBase consensusBase = AlignedBase::fromInteger(uint64_t(it - baseCoverage.begin()));
        const uint64_t coverage = *it;

        alignedConsensus[i] = consensusBase;
        if(not consensusBase.isGap()) {
            consensus.push_back(make_pair(Base(consensusBase), coverage));
        }
    }

    poasta_free_msa(msa);
    poasta_free_graph(graph);
}
