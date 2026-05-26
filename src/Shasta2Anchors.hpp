#pragma once

// Shasta2Anchors.hpp

#include "Reads.hpp"
#include "Marker.hpp"
#include "MarkerKmers.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
#include "MemoryMappedVectorOfVectors.hpp"

#include "MarkerGraph.hpp"
#include "span.hpp"
#include "invalid.hpp"

#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <limits>

namespace dinara {
    using Shasta2AnchorId = uint64_t;

    // Forward declarations.
    class Shasta2Anchors;
    class Shasta2Journeys;
    class Shasta2AnchorMarkerInfo;
    class Shasta2Anchor;
    class Shasta2AnchorPairInfo;

    using Shasta2AnchorBaseClass = span<const Shasta2AnchorMarkerInfo>;
    string shasta2AnchorIdToString(Shasta2AnchorId);
    Shasta2AnchorId shasta2AnchorIdFromString(const string&);
}


class dinara::Shasta2AnchorMarkerInfo {
public:
    OrientedReadId orientedReadId;
    uint32_t ordinal;
    uint32_t positionInJourney = invalid<uint32_t>;

    Shasta2AnchorMarkerInfo() {}

    Shasta2AnchorMarkerInfo(OrientedReadId orientedReadId, uint32_t ordinal) :
        orientedReadId(orientedReadId),
        ordinal(ordinal)
    {}

    bool operator<(const Shasta2AnchorMarkerInfo& that) const
    {
        return orientedReadId < that.orientedReadId;
    }
};


class dinara::Shasta2Anchor : public Shasta2AnchorBaseClass {
public:
    Shasta2Anchor(const Shasta2AnchorBaseClass& s) : Shasta2AnchorBaseClass(s) {}

    void check() const;

    uint64_t coverage() const
    {
        return size();
    }
};


class dinara::Shasta2AnchorPairInfo {
public:
    uint64_t totalA = 0;
    uint64_t totalB = 0;
    uint64_t common = 0;
    uint64_t onlyA = 0;
    uint64_t onlyB = 0;

    int64_t offsetInMarkers = invalid<int64_t>;
    int64_t offsetInBases = invalid<int64_t>;
    uint64_t onlyAShort = invalid<uint64_t>;
    uint64_t onlyBShort = invalid<uint64_t>;

    uint64_t intersectionCount() const { return common; }
    uint64_t unionCount() const { return totalA + totalB - common; }
    uint64_t correctedUnionCount() const { return unionCount() - onlyAShort - onlyBShort; }
    double jaccard() const { return double(intersectionCount()) / double(unionCount()); }
    double correctedJaccard() const { return double(intersectionCount()) / double(correctedUnionCount()); }

    void reverse()
    {
        std::swap(totalA, totalB);
        std::swap(onlyA, onlyB);
        std::swap(onlyAShort, onlyBShort);
        offsetInMarkers = - offsetInMarkers;
        offsetInBases = - offsetInBases;
    }
};


class dinara::Shasta2Anchors :
    public MultithreadedObject<Shasta2Anchors>,
    public MappedMemoryOwner {
public:


    // Constructor from MarkerGraph vertices.
    Shasta2Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph& markerGraph,
        uint64_t threadCount,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage);

    // Minimal constructor: initializes members and builds MarkerKmers only.
    // Does not build anchorMarkerInfos. Caller is responsible for populating them.
    Shasta2Anchors(
        const MappedMemoryOwner&,
        const Reads& reads,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph& markerGraph,
        uint64_t threadCount);

    // Access existing.
    Shasta2Anchors(
        const string& baseName,
        const MappedMemoryOwner&, // mappedMemoryOwner
        const Reads&,
        uint64_t k,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>&
        );

    Shasta2Anchor operator[](Shasta2AnchorId) const;
    uint64_t size() const;
    
    // Kmer helpers (mimicking Markers::getKmer).
    Kmer getKmer(OrientedReadId, uint32_t ordinal) const;
    // Helper to get sequence (for debug/display).
    vector<Base> anchorKmerSequence(Shasta2AnchorId) const;
    Kmer anchorKmer(Shasta2AnchorId) const;
    const MarkerKmers* getMarkerKmers() const
    {
        return markerKmers.get();
    }
    uint64_t filterByShasta2HashedKmerChecker(double markerDensity);
    uint64_t writeExternalAnchors(const string& name, bool canonicalOnly = true) const;


    // Read composition analysis.
    uint64_t countCommon(Shasta2AnchorId anchorId0, Shasta2AnchorId anchorId1) const;
    uint64_t countCommon(Shasta2AnchorId anchorId0, Shasta2AnchorId anchorId1, uint64_t& baseOffset) const;
    
    void analyzeAnchorPair(Shasta2AnchorId, Shasta2AnchorId, Shasta2AnchorPairInfo&) const;
    
    // HTML writing omitted for now unless requested, to keep port focused.

    MemoryMapped::VectorOfVectors<Shasta2AnchorMarkerInfo, uint64_t> anchorMarkerInfos;

    uint32_t getOrdinal(Shasta2AnchorId, OrientedReadId) const;
    uint32_t getPositionInJourney(Shasta2AnchorId, OrientedReadId) const;
    const Shasta2AnchorMarkerInfo& getAnchorMarkerInfo(Shasta2AnchorId, OrientedReadId) const;
    bool anchorContains(Shasta2AnchorId, OrientedReadId) const;

    // Journeys interaction.
    void findChildren(
        const Shasta2Journeys&,
        Shasta2AnchorId,
        vector<Shasta2AnchorId>&,
        vector<uint64_t>& count,
        uint64_t minCoverage = 0) const;

    void findParents(
        const Shasta2Journeys&,
        Shasta2AnchorId,
        vector<Shasta2AnchorId>&,
        vector<uint64_t>& count,
        uint64_t minCoverage = 0) const;

    const Reads& reads;
    const uint64_t k;
    const uint64_t kHalf;
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers;
    const MarkerGraph& markerGraph;

private:
    shared_ptr<MarkerKmers> markerKmers;
    vector<MarkerGraphVertexId> anchorVertexIds;
   
    // Data and functions used when constructing the Anchors.
    class ConstructData {
    public:
        uint64_t minAnchorCoverage;
        uint64_t maxAnchorCoverage;
        // vector<uint64_t> maxAnchorRepeatLength;
        MemoryMapped::Vector<uint64_t> coverage;
        vector<MarkerGraphVertexId> selectedVertexIds;
    };
    ConstructData constructData;
    void constructThreadFunctionPass1(uint64_t threadId);
    void constructThreadFunctionPass2(uint64_t threadId);
    
    // Helpers.
    Kmer getKmerStrand0(ReadId, uint32_t ordinal) const;
    Kmer getKmerStrand1(ReadId, uint32_t ordinal) const;
};
