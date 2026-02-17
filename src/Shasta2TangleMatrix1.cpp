#include "Shasta2TangleMatrix1.hpp"

#include "deduplicate.hpp"

#include <algorithm>
#include <numeric>

using namespace dinara;
using namespace std;

Shasta2TangleMatrix1::Shasta2TangleMatrix1(
    const Shasta2AssemblyGraph& assemblyGraph,
    vector<edge_descriptor> entrances,
    vector<edge_descriptor> exits,
    ostream& html) :
    assemblyGraph(assemblyGraph),
    entrances(entrances),
    exits(exits)
{
    sort(
        this->entrances.begin(),
        this->entrances.end(),
        [&assemblyGraph](const edge_descriptor e0, const edge_descriptor e1) {
            return assemblyGraph[e0].id < assemblyGraph[e1].id;
        });
    sort(
        this->exits.begin(),
        this->exits.end(),
        [&assemblyGraph](const edge_descriptor e0, const edge_descriptor e1) {
            return assemblyGraph[e0].id < assemblyGraph[e1].id;
        });

    gatherOrientedReads(uint32_t(assemblyGraph.getOptions().representativeRegionStepCount));
    gatherCommonOrientedReads();
    computeTotalTangleMatrix();
    (void)html;
}

void Shasta2TangleMatrix1::gatherOrientedReads(const uint64_t representativeRegionStepCount)
{
    entranceOrientedReadInfos.resize(entrances.size());
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        gatherEntranceOrientedReads(iEntrance, representativeRegionStepCount);
    }

    exitOrientedReadInfos.resize(exits.size());
    for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
        gatherExitOrientedReads(iExit, representativeRegionStepCount);
    }
}

void Shasta2TangleMatrix1::gatherEntranceOrientedReads(
    const uint64_t iEntrance,
    const uint64_t representativeRegionStepCount)
{
    DINARA_ASSERT(assemblyGraph.getAnchorsPointer() != nullptr);
    const Shasta2Anchors& anchors = *assemblyGraph.getAnchorsPointer();

    const edge_descriptor e = entrances[iEntrance];
    const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
    const uint64_t stepCount = edge.size();
    const uint64_t begin = (representativeRegionStepCount >= stepCount) ? 0 : (stepCount - representativeRegionStepCount);

    struct WorkInfo {
        OrientedReadId orientedReadId;
        uint32_t positionInJourney;
        bool operator<(const WorkInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };
    vector<WorkInfo> workInfos;

    for(uint64_t stepId=begin; stepId<stepCount; stepId++) {
        const Shasta2AssemblyGraphEdgeStep& step = edge[stepId];
        const Shasta2AnchorId anchorIdB = step.anchorPair.anchorIdB;
        for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
            const auto& infoB = anchors.getAnchorMarkerInfo(anchorIdB, orientedReadId);
            workInfos.push_back({orientedReadId, infoB.positionInJourney});
        }
    }

    sort(workInfos.begin(), workInfos.end());
    auto& orientedReadInfos = entranceOrientedReadInfos[iEntrance];
    for(auto streakBegin=workInfos.begin(); streakBegin!=workInfos.end();) {
        const OrientedReadId orientedReadId = streakBegin->orientedReadId;
        auto streakEnd = streakBegin;
        while(streakEnd != workInfos.end() && streakEnd->orientedReadId == orientedReadId) {
            ++streakEnd;
        }

        OrientedReadInfo info;
        info.orientedReadId = orientedReadId;
        info.stepCount = streakEnd - streakBegin;
        info.positionInJourney = 0;
        for(auto it=streakBegin; it!=streakEnd; ++it) {
            info.positionInJourney = max(info.positionInJourney, it->positionInJourney);
        }
        orientedReadInfos.push_back(info);
        streakBegin = streakEnd;
    }
}

void Shasta2TangleMatrix1::gatherExitOrientedReads(
    const uint64_t iExit,
    const uint64_t representativeRegionStepCount)
{
    DINARA_ASSERT(assemblyGraph.getAnchorsPointer() != nullptr);
    const Shasta2Anchors& anchors = *assemblyGraph.getAnchorsPointer();

    const edge_descriptor e = exits[iExit];
    const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
    const uint64_t stepCount = edge.size();
    const uint64_t end = min(stepCount, representativeRegionStepCount);

    struct WorkInfo {
        OrientedReadId orientedReadId;
        uint32_t positionInJourney;
        bool operator<(const WorkInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };
    vector<WorkInfo> workInfos;

    for(uint64_t stepId=0; stepId<end; stepId++) {
        const Shasta2AssemblyGraphEdgeStep& step = edge[stepId];
        const Shasta2AnchorId anchorIdA = step.anchorPair.anchorIdA;
        for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
            const auto& infoA = anchors.getAnchorMarkerInfo(anchorIdA, orientedReadId);
            workInfos.push_back({orientedReadId, infoA.positionInJourney});
        }
    }

    sort(workInfos.begin(), workInfos.end());
    auto& orientedReadInfos = exitOrientedReadInfos[iExit];
    for(auto streakBegin=workInfos.begin(); streakBegin!=workInfos.end();) {
        const OrientedReadId orientedReadId = streakBegin->orientedReadId;
        auto streakEnd = streakBegin;
        while(streakEnd != workInfos.end() && streakEnd->orientedReadId == orientedReadId) {
            ++streakEnd;
        }

        OrientedReadInfo info;
        info.orientedReadId = orientedReadId;
        info.stepCount = streakEnd - streakBegin;
        info.positionInJourney = numeric_limits<uint32_t>::max();
        for(auto it=streakBegin; it!=streakEnd; ++it) {
            info.positionInJourney = min(info.positionInJourney, it->positionInJourney);
        }
        orientedReadInfos.push_back(info);
        streakBegin = streakEnd;
    }
}

Shasta2TangleMatrix1::CommonOrientedReadInfo::CommonOrientedReadInfo(
    const OrientedReadId orientedReadId,
    const uint64_t entranceCount,
    const uint64_t exitCount) :
    orientedReadId(orientedReadId),
    entranceStepCount(entranceCount, 0),
    exitStepCount(exitCount, 0)
{}

void Shasta2TangleMatrix1::CommonOrientedReadInfo::computeTangleMatrix()
{
    const double entranceSum = double(accumulate(entranceStepCount.begin(), entranceStepCount.end(), uint64_t(0)));
    const double exitSum = double(accumulate(exitStepCount.begin(), exitStepCount.end(), uint64_t(0)));
    tangleMatrix.assign(entranceStepCount.size(), vector<double>(exitStepCount.size(), 0.));

    if(entranceSum == 0. || exitSum == 0.) {
        return;
    }
    for(uint64_t i=0; i<entranceStepCount.size(); i++) {
        for(uint64_t j=0; j<exitStepCount.size(); j++) {
            tangleMatrix[i][j] =
                double(entranceStepCount[i] * exitStepCount[j]) / (entranceSum * exitSum);
        }
    }
}

void Shasta2TangleMatrix1::gatherCommonOrientedReads()
{
    vector<OrientedReadId> entranceOrientedReadIds;
    for(const auto& orientedReadInfos: entranceOrientedReadInfos) {
        for(const auto& orientedReadInfo: orientedReadInfos) {
            entranceOrientedReadIds.push_back(orientedReadInfo.orientedReadId);
        }
    }
    deduplicate(entranceOrientedReadIds);

    vector<OrientedReadId> exitOrientedReadIds;
    for(const auto& orientedReadInfos: exitOrientedReadInfos) {
        for(const auto& orientedReadInfo: orientedReadInfos) {
            exitOrientedReadIds.push_back(orientedReadInfo.orientedReadId);
        }
    }
    deduplicate(exitOrientedReadIds);

    vector<OrientedReadId> commonOrientedReadIds;
    set_intersection(
        entranceOrientedReadIds.begin(),
        entranceOrientedReadIds.end(),
        exitOrientedReadIds.begin(),
        exitOrientedReadIds.end(),
        back_inserter(commonOrientedReadIds));

    commonOrientedReadInfos.clear();
    vector<vector<OrientedReadInfo>::const_iterator> entranceIt(entranceOrientedReadInfos.size());
    for(uint64_t i=0; i<entranceOrientedReadInfos.size(); i++) {
        entranceIt[i] = entranceOrientedReadInfos[i].begin();
    }
    vector<vector<OrientedReadInfo>::const_iterator> exitIt(exitOrientedReadInfos.size());
    for(uint64_t i=0; i<exitOrientedReadInfos.size(); i++) {
        exitIt[i] = exitOrientedReadInfos[i].begin();
    }

    for(const OrientedReadId orientedReadId: commonOrientedReadIds) {
        commonOrientedReadInfos.emplace_back(
            orientedReadId,
            entranceOrientedReadInfos.size(),
            exitOrientedReadInfos.size());
        CommonOrientedReadInfo& commonInfo = commonOrientedReadInfos.back();

        for(uint64_t i=0; i<entranceOrientedReadInfos.size(); i++) {
            auto& it = entranceIt[i];
            while((it != entranceOrientedReadInfos[i].end()) && (it->orientedReadId < orientedReadId)) {
                ++it;
            }
            if((it != entranceOrientedReadInfos[i].end()) && (it->orientedReadId == orientedReadId)) {
                commonInfo.entranceStepCount[i] = it->stepCount;
                commonInfo.maxPositionInJourneyOnEntrances =
                    max(commonInfo.maxPositionInJourneyOnEntrances, it->positionInJourney);
            }
        }

        for(uint64_t i=0; i<exitOrientedReadInfos.size(); i++) {
            auto& it = exitIt[i];
            while((it != exitOrientedReadInfos[i].end()) && (it->orientedReadId < orientedReadId)) {
                ++it;
            }
            if((it != exitOrientedReadInfos[i].end()) && (it->orientedReadId == orientedReadId)) {
                commonInfo.exitStepCount[i] = it->stepCount;
                commonInfo.minPositionInJourneyOnExits =
                    min(commonInfo.minPositionInJourneyOnExits, it->positionInJourney);
            }
        }

        commonInfo.computeTangleMatrix();
    }
}

void Shasta2TangleMatrix1::computeTotalTangleMatrix()
{
    tangleMatrix.assign(entrances.size(), vector<double>(exits.size(), 0.));
    for(const auto& commonInfo: commonOrientedReadInfos) {
        if(commonInfo.goesBackward()) {
            continue;
        }
        for(uint64_t i=0; i<entrances.size(); i++) {
            for(uint64_t j=0; j<exits.size(); j++) {
                tangleMatrix[i][j] += commonInfo.tangleMatrix[i][j];
            }
        }
    }
}

bool Shasta2TangleMatrix1::goesBackward(const OrientedReadId orientedReadId) const
{
    const CommonOrientedReadInfo target(orientedReadId);
    const auto it = lower_bound(commonOrientedReadInfos.begin(), commonOrientedReadInfos.end(), target);
    if(it == commonOrientedReadInfos.end() || it->orientedReadId != orientedReadId) {
        return false;
    }
    return it->goesBackward();
}

uint64_t Shasta2TangleMatrix1::getCommonOrientedReadIdIndex(const OrientedReadId orientedReadId) const
{
    const CommonOrientedReadInfo target(orientedReadId);
    const auto it = lower_bound(commonOrientedReadInfos.begin(), commonOrientedReadInfos.end(), target);
    if(it == commonOrientedReadInfos.end() || it->orientedReadId != orientedReadId) {
        return invalid<uint64_t>;
    }
    return uint64_t(it - commonOrientedReadInfos.begin());
}
