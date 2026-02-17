#include "Shasta2MarkerKmerPair.hpp"
#include "Shasta2AbpoaWrapper.hpp"
#include "Marker.hpp"
#include "Reads.hpp"
#include "orderPairs.hpp"
#include "DINARA_ASSERT.hpp"

#include <algorithm>

using namespace dinara;

Shasta2MarkerKmerPair::Shasta2MarkerKmerPair(
    const MarkerKmers& markerKmers,
    const Kmer& kmer0,
    const Kmer& kmer1,
    uint64_t maxPositionOffset) :
    kmer0(kmer0),
    kmer1(kmer1)
{
    getMarkerInfos(markerKmers);
    gatherCommonOrientedReads(markerKmers, maxPositionOffset);
    gatherSequences(markerKmers.getReads());
    rankSequences();
    align();
}

void Shasta2MarkerKmerPair::getMarkerInfos(const MarkerKmers& markerKmers)
{
    markerKmers.get(kmer0, markerInfos0);
    markerKmers.get(kmer1, markerInfos1);
}

void Shasta2MarkerKmerPair::gatherCommonOrientedReads(
    const MarkerKmers& markerKmers,
    uint64_t maxPositionOffset)
{
    const auto& markers = markerKmers.getMarkers();
    const uint32_t kHalf = uint32_t(markerKmers.getK() / 2);

    auto it0 = markerInfos0.begin();
    const auto end0 = markerInfos0.end();
    auto it1 = markerInfos1.begin();
    const auto end1 = markerInfos1.end();
    while((it0 != end0) && (it1 != end1)) {
        if(it0->orientedReadId < it1->orientedReadId) {
            ++it0;
            continue;
        }
        if(it1->orientedReadId < it0->orientedReadId) {
            ++it1;
            continue;
        }

        const OrientedReadId orientedReadId = it0->orientedReadId;
        DINARA_ASSERT(orientedReadId == it1->orientedReadId);

        const auto streakBegin0 = it0;
        const auto streakBegin1 = it1;
        auto streakEnd0 = streakBegin0;
        auto streakEnd1 = streakBegin1;
        while((streakEnd0 != end0) && (streakEnd0->orientedReadId == orientedReadId)) {
            ++streakEnd0;
        }
        while((streakEnd1 != end1) && (streakEnd1->orientedReadId == orientedReadId)) {
            ++streakEnd1;
        }
        const uint64_t streakLength0 = streakEnd0 - streakBegin0;
        const uint64_t streakLength1 = streakEnd1 - streakBegin1;

        if(streakLength0 == 1 && streakLength1 == 1) {
            if(it0->ordinal < it1->ordinal) {
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];
                const uint32_t position0 = orientedReadMarkers[it0->ordinal].position + kHalf;
                const uint32_t position1 = orientedReadMarkers[it1->ordinal].position + kHalf;
                if(position1 - position0 <= maxPositionOffset) {
                    commonOrientedReads.emplace_back(orientedReadId, it0->ordinal, it1->ordinal, position0, position1);
                }
            }
        } else if(streakLength0 == 1 && streakLength1 > 1) {
            const uint32_t ordinal0 = it0->ordinal;
            auto kt1 = streakBegin1;
            for(; kt1 != streakEnd1; ++kt1) {
                if(kt1->ordinal > ordinal0) {
                    break;
                }
            }
            if(kt1 != streakEnd1) {
                const uint32_t ordinal1 = kt1->ordinal;
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];
                const uint32_t position0 = orientedReadMarkers[ordinal0].position + kHalf;
                const uint32_t position1 = orientedReadMarkers[ordinal1].position + kHalf;
                if(position1 - position0 <= maxPositionOffset) {
                    commonOrientedReads.emplace_back(orientedReadId, ordinal0, ordinal1, position0, position1);
                }
            }
        } else if(streakLength1 == 1 && streakLength0 > 1) {
            const uint32_t ordinal1 = it1->ordinal;
            const auto streakReverseBegin0 = streakEnd0 - 1;
            const auto streakReverseEnd0 = streakBegin0 - 1;
            auto kt0 = streakReverseBegin0;
            for(; kt0 != streakReverseEnd0; --kt0) {
                if(kt0->ordinal < ordinal1) {
                    break;
                }
            }
            if(kt0 != streakReverseEnd0) {
                const uint32_t ordinal0 = kt0->ordinal;
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];
                const uint32_t position0 = orientedReadMarkers[ordinal0].position + kHalf;
                const uint32_t position1 = orientedReadMarkers[ordinal1].position + kHalf;
                if(position1 - position0 <= maxPositionOffset) {
                    commonOrientedReads.emplace_back(orientedReadId, ordinal0, ordinal1, position0, position1);
                }
            }
        }

        it0 = streakEnd0;
        it1 = streakEnd1;
    }
}

void Shasta2MarkerKmerPair::gatherSequences(const Reads& reads)
{
    vector<Base> sequence;
    for(uint64_t i=0; i<commonOrientedReads.size(); i++) {
        CommonOrientedRead& commonOrientedRead = commonOrientedReads[i];
        commonOrientedRead.getSequence(reads, sequence);

        auto it = sequenceMap.find(sequence);
        if(it == sequenceMap.end()) {
            tie(it, ignore) = sequenceMap.insert(make_pair(sequence, SequenceInfo()));
        }
        it->second.orientedReadIndexes.push_back(i);
        commonOrientedRead.sequenceMapIterator = it;
    }
}

void Shasta2MarkerKmerPair::CommonOrientedRead::getSequence(
    const Reads& reads,
    vector<Base>& sequence) const
{
    sequence.clear();
    for(uint32_t position=position0; position!=position1; position++) {
        sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
    }
}

void Shasta2MarkerKmerPair::rankSequences()
{
    vector< pair<SequenceMap::iterator, uint64_t> > sequenceTable;
    for(auto it=sequenceMap.begin(); it!=sequenceMap.end(); ++it) {
        sequenceTable.push_back(make_pair(it, it->second.coverage()));
    }

    std::ranges::sort(sequenceTable, OrderPairsBySecondOnlyGreater<SequenceMap::iterator, uint64_t>());

    for(uint64_t rank=0; rank<sequenceTable.size(); rank++) {
        const auto it = sequenceTable[rank].first;
        it->second.rank = rank;
        sequencesByRank.push_back(it);
    }
}

void Shasta2MarkerKmerPair::align()
{
    vector< pair<vector<Base>, uint64_t> > sequencesWithCoverage;
    for(const auto it: sequencesByRank) {
        sequencesWithCoverage.push_back(make_pair(it->first, it->second.coverage()));
    }
    shasta2Abpoa(sequencesWithCoverage, consensus, alignment, alignedConsensus);
}

uint64_t Shasta2MarkerKmerPair::editDistance(uint64_t rank0, uint64_t rank1) const
{
    DINARA_ASSERT(rank0 < sequencesByRank.size());
    DINARA_ASSERT(rank1 < sequencesByRank.size());
    DINARA_ASSERT(rank0 < alignment.size());
    DINARA_ASSERT(rank1 < alignment.size());
    const vector<AlignedBase>& alignedSequence0 = alignment[rank0];
    const vector<AlignedBase>& alignedSequence1 = alignment[rank1];
    DINARA_ASSERT(alignedSequence0.size() == alignedSequence1.size());

    uint64_t editDistance = 0;
    for(uint64_t i=0; i<alignedSequence0.size(); i++) {
        const AlignedBase b0 = alignedSequence0[i];
        const AlignedBase b1 = alignedSequence1[i];
        if((not b0.isGap()) && (not b1.isGap()) && (b0 != b1)) {
            ++editDistance;
        }
    }
    return editDistance;
}
