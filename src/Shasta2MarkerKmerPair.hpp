#pragma once

#include "Base.hpp"
#include "Kmer.hpp"
#include "MarkerKmers.hpp"
#include "ReadId.hpp"
#include "invalid.hpp"

#include <map>
#include "vector.hpp"

namespace dinara {
    class Shasta2MarkerKmerPair;
}

class dinara::Shasta2MarkerKmerPair {
public:
    Shasta2MarkerKmerPair(
        const MarkerKmers&,
        const Kmer& kmer0,
        const Kmer& kmer1,
        uint64_t maxPositionOffset);

    Kmer kmer0;
    Kmer kmer1;

    using MarkerInfo = MarkerKmers::MarkerInfo;
    vector<MarkerInfo> markerInfos0;
    vector<MarkerInfo> markerInfos1;
    void getMarkerInfos(const MarkerKmers&);

    class SequenceInfo {
    public:
        vector<uint64_t> orientedReadIndexes;
        uint64_t coverage() const
        {
            return orientedReadIndexes.size();
        }
        uint64_t rank = invalid<uint64_t>;
    };
    using SequenceMap = std::map<vector<Base>, SequenceInfo>;
    SequenceMap sequenceMap;
    vector<SequenceMap::iterator> sequencesByRank;
    void gatherSequences(const Reads&);
    void rankSequences();

    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;
    void align();

    uint64_t editDistance(uint64_t rank0, uint64_t rank1) const;

    class CommonOrientedRead {
    public:
        OrientedReadId orientedReadId;
        uint32_t ordinal0;
        uint32_t ordinal1;
        uint32_t ordinalOffset() const
        {
            return ordinal1 - ordinal0;
        }
        uint32_t position0;
        uint32_t position1;
        uint32_t positionOffset() const
        {
            return position1 - position0;
        }
        void getSequence(const Reads&, vector<Base>&) const;
        SequenceMap::const_iterator sequenceMapIterator;

        CommonOrientedRead() {}
        CommonOrientedRead(
            OrientedReadId orientedReadId,
            uint32_t ordinal0,
            uint32_t ordinal1,
            uint32_t position0,
            uint32_t position1) :
            orientedReadId(orientedReadId),
            ordinal0(ordinal0),
            ordinal1(ordinal1),
            position0(position0),
            position1(position1)
        {}

        bool operator<(const CommonOrientedRead& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };
    vector<CommonOrientedRead> commonOrientedReads;
    void gatherCommonOrientedReads(const MarkerKmers&, uint64_t maxPositionOffset);
};
