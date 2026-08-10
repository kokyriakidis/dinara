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
#include <unordered_map>

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

    // The position of the middle of the marker relative to the beginning
    // of the oriented read. This equals the position of the first base
    // of the marker plus k/2. Matches shasta2's MarkerInfo::position.
    // This is the field that gets serialized for shasta2 compatibility.
    uint32_t position;

    uint32_t positionInJourney = invalid<uint32_t>;

    // The marker ordinal (index into the oriented read's markers array).
    // Not serialized for shasta2 — kept for internal use (alignment chaining,
    // k-mer lookups, HTTP server display).
    uint32_t ordinal = invalid<uint32_t>;

    Shasta2AnchorMarkerInfo() {}

    Shasta2AnchorMarkerInfo(OrientedReadId orientedReadId, uint32_t position, uint32_t ordinal) :
        orientedReadId(orientedReadId),
        position(position),
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
    
    // Kmer helpers.
    // getKmer takes an ordinal (index into the markers array).
    Kmer getKmer(OrientedReadId, uint32_t ordinal) const;
    // getKmerAtPosition takes a midpoint position (as stored in Shasta2AnchorMarkerInfo::position).
    Kmer getKmerAtPosition(OrientedReadId, uint32_t midpointPosition) const;
    // Helper to get sequence (for debug/display).
    vector<Base> anchorKmerSequence(Shasta2AnchorId) const;
    Kmer anchorKmer(Shasta2AnchorId) const;
    const MarkerKmers* getMarkerKmers() const
    {
        return markerKmers.get();
    }
    // Optional per-canonical-anchor member drop set: canonicalAnchorId ->
    // ReadIds to omit from that anchor at export. Used to resolve global
    // per-read journey position ties (two anchors at the same exported position
    // on one read), which shasta2 rejects with "Invalid Journey". Keyed on the
    // canonical (even) id and ReadId only: shasta2 loads canonicals and
    // regenerates each RC, so one canonical drop removes both the direct and the
    // RC-induced occurrence on that read. See resolveJourneyPositionTies in
    // main.cpp.
    using ExternalAnchorDropMap =
        std::unordered_map<Shasta2AnchorId, std::vector<ReadId>>;
    uint64_t writeExternalAnchors(
        const string& name,
        bool canonicalOnly = true,
        const ExternalAnchorDropMap* dropMap = nullptr) const;

    // Append a k=2 het anchor and its reverse complement to the store, returning
    // the canonical (even) anchor id; the RC is at id+1. Members are (read,
    // rawPosition) where rawPosition is the read base position of the first base
    // of the 2-base marker [predBase, alleleBase], in the oriented read's own
    // coordinate frame. The stored midpoint position follows the k-store
    // convention (rawPosition + k/2) so writeExternalAnchors recovers the true
    // raw position; shasta2 re-derives a 2-base k-mer from it. Must be called in
    // a serial pass AFTER all window processing (it grows the store and thus
    // invalidates any outstanding anchor spans). The first call records
    // hetAnchorFirstId so the export can skip the k=50 k-mer consistency check
    // for het anchors (their members agree on only 2 bases by design).
    Shasta2AnchorId appendHetAnchorPair(
        const vector<std::pair<OrientedReadId, uint32_t>>& members);

    // First anchor id that is a k=2 het anchor (invalid if none). All ids >=
    // this are het anchors, appended contiguously after the primary anchors.
    Shasta2AnchorId hetAnchorFirstId = invalid<Shasta2AnchorId>;


    // Read composition analysis.
    uint64_t countCommon(Shasta2AnchorId anchorId0, Shasta2AnchorId anchorId1) const;
    uint64_t countCommon(Shasta2AnchorId anchorId0, Shasta2AnchorId anchorId1, uint64_t& baseOffset) const;
    
    void analyzeAnchorPair(Shasta2AnchorId, Shasta2AnchorId, Shasta2AnchorPairInfo&) const;
    
    // HTML writing omitted for now unless requested, to keep port focused.

    MemoryMapped::VectorOfVectors<Shasta2AnchorMarkerInfo, uint64_t> anchorMarkerInfos;

    // Each anchor's members are stored sorted ascending by OrientedReadId (see
    // the anchor-construction code and appendHetAnchorPair), so the four
    // lookups below binary search -- O(log coverage) -- rather than scanning
    // the whole anchor.
    uint32_t getPosition(Shasta2AnchorId, OrientedReadId) const;
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
