#include "Shasta2SegmentStepSupport.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>

using namespace dinara;
using namespace std;

void Shasta2SegmentStepSupport::get(
    const Shasta2AssemblyGraph& assemblyGraph,
    const edge_descriptor e,
    const uint32_t stepBegin,
    const uint32_t stepEnd,
    vector<Shasta2SegmentStepSupport>& v)
{
    v.clear();
    for(uint32_t stepId=stepBegin; stepId<stepEnd; stepId++) {
        append(assemblyGraph, e, stepId, v);
    }
}

void Shasta2SegmentStepSupport::getInitial(
    const Shasta2AssemblyGraph& assemblyGraph,
    const edge_descriptor e,
    const uint32_t stepCount,
    vector<Shasta2SegmentStepSupport>& v)
{
    const uint32_t totalStepCount = uint32_t(assemblyGraph[e].size());
    get(assemblyGraph, e, 0, min(stepCount, totalStepCount), v);
}

void Shasta2SegmentStepSupport::getInitialFirst(
    const Shasta2AssemblyGraph& assemblyGraph,
    const edge_descriptor e,
    const uint32_t stepCount,
    vector<Shasta2SegmentStepSupport>& v)
{
    getInitial(assemblyGraph, e, stepCount, v);
    keepFirst(v);
}

void Shasta2SegmentStepSupport::getFinal(
    const Shasta2AssemblyGraph& assemblyGraph,
    const edge_descriptor e,
    const uint32_t stepCount,
    vector<Shasta2SegmentStepSupport>& v)
{
    const uint32_t totalStepCount = uint32_t(assemblyGraph[e].size());
    const uint32_t begin = (totalStepCount >= stepCount) ? (totalStepCount - stepCount) : 0;
    get(assemblyGraph, e, begin, totalStepCount, v);
}

void Shasta2SegmentStepSupport::getFinalLast(
    const Shasta2AssemblyGraph& assemblyGraph,
    const edge_descriptor e,
    const uint32_t stepCount,
    vector<Shasta2SegmentStepSupport>& v)
{
    getFinal(assemblyGraph, e, stepCount, v);
    keepLast(v);
}

void Shasta2SegmentStepSupport::append(
    const Shasta2AssemblyGraph& assemblyGraph,
    const edge_descriptor e,
    const uint32_t stepId,
    vector<Shasta2SegmentStepSupport>& v)
{
    DINARA_ASSERT(assemblyGraph.getAnchorsPointer() != nullptr);
    const Shasta2Anchors& anchors = *assemblyGraph.getAnchorsPointer();
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
    DINARA_ASSERT(stepId < uint32_t(edge.size()));
    const Shasta2AssemblyGraphEdgeStep& step = edge[stepId];
    const Shasta2AnchorPair& anchorPair = step.anchorPair;

    for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
        const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

        v.emplace_back();
        Shasta2SegmentStepSupport& stepSupport = v.back();
        stepSupport.e = e;
        stepSupport.stepId = stepId;
        stepSupport.orientedReadId = orientedReadId;

        const Shasta2AnchorId anchorIdA = anchorPair.anchorIdA;
        const Shasta2AnchorMarkerInfo& infoA = anchors.getAnchorMarkerInfo(anchorIdA, orientedReadId);
        stepSupport.positionInJourneyA = infoA.positionInJourney;
        stepSupport.ordinalA = infoA.ordinal;
        stepSupport.positionA = orientedReadMarkers[infoA.ordinal].position + kHalf;

        const Shasta2AnchorId anchorIdB = anchorPair.anchorIdB;
        const Shasta2AnchorMarkerInfo& infoB = anchors.getAnchorMarkerInfo(anchorIdB, orientedReadId);
        stepSupport.positionInJourneyB = infoB.positionInJourney;
        stepSupport.ordinalB = infoB.ordinal;
        stepSupport.positionB = orientedReadMarkers[infoB.ordinal].position + kHalf;
    }
}

void Shasta2SegmentStepSupport::keepLast(vector<Shasta2SegmentStepSupport>& v)
{
    ranges::sort(
        v,
        ranges::less(),
        [](const Shasta2SegmentStepSupport& s) {
            return tie(s.orientedReadId, s.stepId);
        });

    auto it = v.begin();
    const auto end = v.end();
    auto out = v.begin();
    while(it != end) {
        const OrientedReadId orientedReadId = it->orientedReadId;
        auto streakBegin = it;
        auto streakEnd = streakBegin + 1;
        while((streakEnd != end) && (streakEnd->orientedReadId == orientedReadId)) {
            ++streakEnd;
        }
        *out++ = *(streakEnd - 1);
        it = streakEnd;
    }
    v.resize(out - v.begin());
}

void Shasta2SegmentStepSupport::keepFirst(vector<Shasta2SegmentStepSupport>& v)
{
    ranges::sort(
        v,
        ranges::less(),
        [](const Shasta2SegmentStepSupport& s) {
            return tie(s.orientedReadId, s.stepId);
        });

    auto it = v.begin();
    const auto end = v.end();
    auto out = v.begin();
    while(it != end) {
        const OrientedReadId orientedReadId = it->orientedReadId;
        auto streakBegin = it;
        auto streakEnd = streakBegin + 1;
        while((streakEnd != end) && (streakEnd->orientedReadId == orientedReadId)) {
            ++streakEnd;
        }
        *out++ = *streakBegin;
        it = streakEnd;
    }
    v.resize(out - v.begin());
}

Shasta2SegmentPairInformation Shasta2SegmentStepSupport::analyzeSegmentPair(
    ostream& html,
    const Shasta2AssemblyGraph& assemblyGraph,
    const edge_descriptor e0,
    const edge_descriptor e1,
    const uint32_t representativeRegionStepCount)
{
    const Shasta2AssemblyGraphEdge& edge0 = assemblyGraph[e0];
    const Shasta2AssemblyGraphEdge& edge1 = assemblyGraph[e1];

    vector<Shasta2SegmentStepSupport> support0;
    getFinalLast(assemblyGraph, e0, representativeRegionStepCount, support0);
    vector<Shasta2SegmentStepSupport> support1;
    getInitialFirst(assemblyGraph, e1, representativeRegionStepCount, support1);

    auto it0 = support0.begin();
    auto it1 = support1.begin();
    const auto end0 = support0.end();
    const auto end1 = support1.end();

    uint32_t commonCount = 0;
    int64_t offsetSum = 0;
    while((it0 != end0) && (it1 != end1)) {
        if(it0->orientedReadId < it1->orientedReadId) {
            ++it0;
            continue;
        }
        if(it1->orientedReadId < it0->orientedReadId) {
            ++it1;
            continue;
        }
        ++commonCount;
        offsetSum += estimateOffset(assemblyGraph, *it0, *it1);
        ++it0;
        ++it1;
    }

    if(commonCount == 0) {
        return {};
    }

    const int32_t segmentOffset = int32_t(round(double(offsetSum) / double(commonCount)));
    const uint64_t intersectionCount = commonCount;
    const uint64_t unionCount = support0.size() + support1.size() - intersectionCount;
    const double jaccard = double(intersectionCount) / double(unionCount);

    uint32_t missing0 = 0;
    it0 = support0.begin();
    it1 = support1.begin();
    while(it0 != end0) {
        while((it1 != end1) && (it1->orientedReadId < it0->orientedReadId)) {
            ++it1;
        }
        if((it1 != end1) && (it1->orientedReadId == it0->orientedReadId)) {
            ++it0;
            continue;
        }

        int32_t position = it0->positionB;
        for(uint32_t stepId=it0->stepId + 1; stepId<edge0.size(); stepId++) {
            const uint64_t stepOffset = edge0.wasAssembled ? edge0[stepId].sequence.size() : edge0[stepId].offset;
            position += int32_t(stepOffset);
        }
        position += segmentOffset;
        const uint64_t step1Offset = edge1.wasAssembled ? edge1.front().sequence.size() : edge1.front().offset;
        position += int32_t(step1Offset);
        const int32_t readLength = int32_t(
            assemblyGraph.getAnchorsPointer()->reads.getRead(it0->orientedReadId.getReadId()).baseCount);
        const bool isShort = readLength < position;
        if(!isShort) {
            ++missing0;
        }
        ++it0;
    }

    uint32_t missing1 = 0;
    it0 = support0.begin();
    it1 = support1.begin();
    while(it1 != end1) {
        while((it0 != end0) && (it0->orientedReadId < it1->orientedReadId)) {
            ++it0;
        }
        if((it0 != end0) && (it0->orientedReadId == it1->orientedReadId)) {
            ++it1;
            continue;
        }

        int32_t position = it1->positionA;
        for(uint32_t stepId=0; stepId<it1->stepId; stepId++) {
            const uint64_t stepOffset = edge1.wasAssembled ? edge1[stepId].sequence.size() : edge1[stepId].offset;
            position -= int32_t(stepOffset);
        }
        position -= segmentOffset;
        const uint64_t step0Offset = edge0.wasAssembled ? edge0.back().sequence.size() : edge0.back().offset;
        position -= int32_t(step0Offset);
        const bool isShort = position < 0;
        if(!isShort) {
            ++missing1;
        }
        ++it1;
    }

    const uint64_t correctedUnionSize = missing0 + missing1 + commonCount;
    const double correctedJaccard = double(intersectionCount) / double(correctedUnionSize);

    if(html) {
        html << "common=" << commonCount << " jaccard=" << fixed << setprecision(3) << jaccard
            << " correctedJaccard=" << fixed << setprecision(3) << correctedJaccard;
    }

    Shasta2SegmentPairInformation info;
    info.commonCount = commonCount;
    info.missing0 = missing0;
    info.missing1 = missing1;
    info.segmentOffset = segmentOffset;
    info.jaccard = jaccard;
    info.correctedJaccard = correctedJaccard;
    return info;
}

int32_t Shasta2SegmentStepSupport::estimateOffset(
    const Shasta2AssemblyGraph& assemblyGraph,
    Shasta2SegmentStepSupport& s0,
    Shasta2SegmentStepSupport& s1)
{
    DINARA_ASSERT(s0.orientedReadId == s1.orientedReadId);

    const Shasta2AssemblyGraphEdge& edge0 = assemblyGraph[s0.e];
    const Shasta2AssemblyGraphEdge& edge1 = assemblyGraph[s1.e];

    const uint32_t position0 = s0.positionB;
    const uint32_t position1 = s1.positionA;

    uint64_t offset0 = 0;
    for(uint32_t stepId=s0.stepId + 1; stepId<edge0.size(); stepId++) {
        offset0 += edge0.wasAssembled ? edge0[stepId].sequence.size() : edge0[stepId].offset;
    }

    uint64_t offset1 = 0;
    for(uint32_t stepId=0; stepId<s1.stepId; stepId++) {
        offset1 += edge1.wasAssembled ? edge1[stepId].sequence.size() : edge1[stepId].offset;
    }

    return int32_t(position1) - int32_t(position0) - int32_t(offset0) - int32_t(offset1);
}
