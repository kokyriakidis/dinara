#pragma once

#include "Shasta2AssemblyGraph.hpp"

namespace dinara {
    class Shasta2SegmentStepSupport;
    class Shasta2SegmentPairInformation;
}

class dinara::Shasta2SegmentPairInformation {
public:
    uint64_t commonCount = 0;
    uint64_t missing0 = 0;
    uint64_t missing1 = 0;
    int32_t segmentOffset = invalid<int32_t>;
    double jaccard = 0.;
    double correctedJaccard = 0.;
};

class dinara::Shasta2SegmentStepSupport {
public:
    using edge_descriptor = Shasta2AssemblyGraph::edge_descriptor;

    edge_descriptor e;
    uint32_t stepId = 0;
    OrientedReadId orientedReadId;

    uint32_t positionInJourneyA = 0;
    uint32_t ordinalA = 0;
    uint32_t positionA = 0;

    uint32_t positionInJourneyB = 0;
    uint32_t ordinalB = 0;
    uint32_t positionB = 0;

    uint32_t positionInJourneyOffset() const
    {
        return positionInJourneyB - positionInJourneyA;
    }
    uint32_t ordinalOffset() const
    {
        return ordinalB - ordinalA;
    }
    uint32_t positionOffset() const
    {
        return positionB - positionA;
    }

    static void get(
        const Shasta2AssemblyGraph&,
        edge_descriptor,
        uint32_t stepBegin,
        uint32_t stepEnd,
        vector<Shasta2SegmentStepSupport>&);
    static void getInitial(
        const Shasta2AssemblyGraph&,
        edge_descriptor,
        uint32_t stepCount,
        vector<Shasta2SegmentStepSupport>&);
    static void getInitialFirst(
        const Shasta2AssemblyGraph&,
        edge_descriptor,
        uint32_t stepCount,
        vector<Shasta2SegmentStepSupport>&);
    static void getFinal(
        const Shasta2AssemblyGraph&,
        edge_descriptor,
        uint32_t stepCount,
        vector<Shasta2SegmentStepSupport>&);
    static void getFinalLast(
        const Shasta2AssemblyGraph&,
        edge_descriptor,
        uint32_t stepCount,
        vector<Shasta2SegmentStepSupport>&);
    static void append(
        const Shasta2AssemblyGraph&,
        edge_descriptor,
        uint32_t stepId,
        vector<Shasta2SegmentStepSupport>&);
    static void keepLast(vector<Shasta2SegmentStepSupport>&);
    static void keepFirst(vector<Shasta2SegmentStepSupport>&);

    // Output a vector of Shasta2SegmentStepSupport to a html table.
    static void writeHtml(
        ostream& html,
        const Shasta2AssemblyGraph&,
        const vector<Shasta2SegmentStepSupport>&);

    static Shasta2SegmentPairInformation analyzeSegmentPair(
        ostream& html,
        const Shasta2AssemblyGraph&,
        edge_descriptor e0,
        edge_descriptor e1,
        uint32_t representativeRegionStepCount);

    static int32_t estimateOffset(
        const Shasta2AssemblyGraph&,
        Shasta2SegmentStepSupport&,
        Shasta2SegmentStepSupport&);
};
