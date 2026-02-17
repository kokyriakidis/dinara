#pragma once

#include "Shasta2AssemblyGraph.hpp"

#include <limits>

namespace dinara {
    class Shasta2TangleMatrix1;
}

class dinara::Shasta2TangleMatrix1 {
public:
    using edge_descriptor = Shasta2AssemblyGraph::edge_descriptor;

    Shasta2TangleMatrix1(
        const Shasta2AssemblyGraph& assemblyGraph,
        vector<edge_descriptor> entrances,
        vector<edge_descriptor> exits,
        ostream& html);

    const Shasta2AssemblyGraph& assemblyGraph;
    vector<edge_descriptor> entrances;
    vector<edge_descriptor> exits;

    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;
        uint64_t stepCount = 0;
        uint32_t positionInJourney = 0;
        bool operator<(const OrientedReadInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };
    vector< vector<OrientedReadInfo> > entranceOrientedReadInfos;
    vector< vector<OrientedReadInfo> > exitOrientedReadInfos;

    class CommonOrientedReadInfo {
    public:
        OrientedReadId orientedReadId;
        vector<uint64_t> entranceStepCount;
        vector<uint64_t> exitStepCount;
        uint32_t maxPositionInJourneyOnEntrances = 0;
        uint32_t minPositionInJourneyOnExits = std::numeric_limits<uint32_t>::max();
        vector< vector<double> > tangleMatrix;

        CommonOrientedReadInfo(
            OrientedReadId orientedReadId,
            uint64_t entranceCount = 0,
            uint64_t exitCount = 0);
        bool operator<(const CommonOrientedReadInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
        bool goesBackward() const
        {
            return maxPositionInJourneyOnEntrances > minPositionInJourneyOnExits;
        }
        void computeTangleMatrix();
    };
    vector<CommonOrientedReadInfo> commonOrientedReadInfos;

    vector< vector<double> > tangleMatrix;

    void gatherOrientedReads(uint64_t representativeRegionStepCount);
    void gatherEntranceOrientedReads(uint64_t iEntrance, uint64_t representativeRegionStepCount);
    void gatherExitOrientedReads(uint64_t iExit, uint64_t representativeRegionStepCount);
    void gatherCommonOrientedReads();
    void computeTotalTangleMatrix();

    bool goesBackward(OrientedReadId orientedReadId) const;
    uint64_t getCommonOrientedReadIdIndex(OrientedReadId orientedReadId) const;
};
