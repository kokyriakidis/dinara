#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "deduplicate.hpp"
#include "findMarkerId.hpp"
#include "MultithreadedObject.tpp"
#include "Reads.hpp"

#include "Marker.hpp"
#include "MarkerGraph.hpp"

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

using namespace dinara;
using namespace std;

namespace {
    const MarkerGraph& shastaAnchorsDummyMarkerGraph()
    {
        static const MarkerGraph dummy;
        return dummy;
    }
}


Shasta2Anchors::Shasta2Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,

    const MarkerGraph& markerGraph,
    uint64_t threadCount) :
    MultithreadedObject<Shasta2Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    kHalf(k/2),
    markers(markers),
    markerGraph(markerGraph)
{
    cout << "Computing Shasta2 MarkerKmers..." << endl;
    markerKmers = make_shared<MarkerKmers>(
        k,
        *this,
        reads,
        markers,
        threadCount);

    // Pass 1: compute coverage for each anchor (VertexId).
    // The number of anchors is determined by the number of vertices in the markerGraph.
    const uint64_t vertexCount = markerGraph.vertexCount();

    anchorMarkerInfos.createNew(
        largeDataName("Shasta2Anchors-AnchorMarkerInfos"),
        largeDataPageSize);
    
    // Initialize Pass 1 for VectorOfVectors
    anchorMarkerInfos.beginPass1(vertexCount);
    
    const uint64_t batchSize = 1000;

    // Run Pass 1 in parallel to count entries.
    setupLoadBalancing(vertexCount, batchSize);
    runThreads(&Shasta2Anchors::constructThreadFunctionPass1, threadCount);
    
    // Initialize Pass 2
    anchorMarkerInfos.beginPass2();
    
    // Run Pass 2 in parallel to fill entries.
    setupLoadBalancing(vertexCount, batchSize);
    runThreads(&Shasta2Anchors::constructThreadFunctionPass2, threadCount);
    
    // Finalize VectorOfVectors
    anchorMarkerInfos.endPass2();
}

Shasta2Anchors::Shasta2Anchors(
    const string& baseName,
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) :
    MultithreadedObject<Shasta2Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    kHalf(k/2),
    markers(markers),
    markerGraph(shastaAnchorsDummyMarkerGraph())
{
    const string objectName = baseName.empty() ?
        "Shasta2Anchors-AnchorMarkerInfos" :
        (baseName + "-AnchorMarkerInfos");
    anchorMarkerInfos.accessExistingReadOnly(largeDataName(objectName));

    try {
        markerKmers = make_shared<MarkerKmers>(
            k,
            *this,
            reads,
            markers);
    } catch(...) {
        markerKmers.reset();
    }
}

void Shasta2Anchors::constructThreadFunctionPass1(uint64_t threadId)
{
    uint64_t begin, end;

    vector<OrientedReadId> orientedReadIds;

    while(getNextBatch(begin, end)) {
        for(uint64_t vertexId=begin; vertexId!=end; vertexId++) {
            
            orientedReadIds.clear();

            // This anchor corresponds exactly to the MarkerGraph vertex.
            // Do not mix in the reverse complement vertex.
            for(const MarkerId markerId : markerGraph.getVertexMarkerIds(vertexId)) {
                OrientedReadId orientedReadId;
                uint32_t ordinal;
                tie(orientedReadId, ordinal) = findMarkerId(markerId, markers);
                (void)ordinal;
                orientedReadIds.push_back(orientedReadId);
            }

            // Deduplicate orientedReadIds to get unique count.
            // Shasta2Anchors enforce one marker per oriented readID per anchor.
            std::sort(orientedReadIds.begin(), orientedReadIds.end());
            auto last = std::unique(orientedReadIds.begin(), orientedReadIds.end());
            uint64_t count = std::distance(orientedReadIds.begin(), last);
            
            anchorMarkerInfos.incrementCount(vertexId, count);
        }
    }
}

void Shasta2Anchors::constructThreadFunctionPass2(uint64_t threadId)
{
    uint64_t begin, end;
    
    vector<Shasta2AnchorMarkerInfo> buffer;

    while(getNextBatch(begin, end)) {
        for(uint64_t vertexId=begin; vertexId!=end; vertexId++) {
            
            buffer.clear();
            
            // Add markers from this vertex only.
            // Keep the OrientedReadId and ordinal exactly as stored in markers.
            for(const MarkerId markerId : markerGraph.getVertexMarkerIds(vertexId)) {
                OrientedReadId orientedReadId;
                uint32_t ordinal;
                tie(orientedReadId, ordinal) = findMarkerId(markerId, markers);
                buffer.emplace_back(orientedReadId, ordinal);
            }
            
            // Sort by OrientedReadId to ensure canonical order for the Anchor.
            std::sort(buffer.begin(), buffer.end());

            // Deduplicate: keep only unique OrientedReadIds.
            auto last = std::unique(buffer.begin(), buffer.end(),
                [](const Shasta2AnchorMarkerInfo& a, const Shasta2AnchorMarkerInfo& b){
                    return a.orientedReadId == b.orientedReadId;
                });
            buffer.erase(last, buffer.end());
            
            // Store in reverse because store() fills each vector from end to begin.
            // This preserves ascending OrientedReadId order in final storage.
            for(auto it = buffer.rbegin(); it != buffer.rend(); ++it) {
                anchorMarkerInfos.store(vertexId, *it);
            }
        }
    }
}




Shasta2Anchor Shasta2Anchors::operator[](Shasta2AnchorId anchorId) const {
    return Shasta2Anchor(anchorMarkerInfos[anchorId]);
}

uint64_t Shasta2Anchors::size() const {
    return anchorMarkerInfos.size();
}

// Helpers.
Kmer Shasta2Anchors::getKmer(OrientedReadId orientedReadId, uint32_t ordinal) const {
    if(orientedReadId.getStrand() == 0) {
        return getKmerStrand0(orientedReadId.getReadId(), ordinal);
    } else {
        return getKmerStrand1(orientedReadId.getReadId(), ordinal);
    }
}

Kmer Shasta2Anchors::getKmerStrand0(ReadId readId, uint32_t ordinal) const {
    const auto read = reads.getRead(readId);
    const auto orientedReadMarkers = markers[OrientedReadId(readId, 0).getValue()];
    uint64_t position = orientedReadMarkers[ordinal].position;
    
    Kmer kmer;
    for(uint64_t i=0; i<k; i++) {
        kmer.set(i, read[position + i]);
    }
    return kmer;
}

Kmer Shasta2Anchors::getKmerStrand1(ReadId readId, uint32_t ordinal1) const {
    const auto read = reads.getRead(readId);
    const auto orientedReadMarkers0 = markers[OrientedReadId(readId, 0).getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    const uint64_t ordinal0 = readMarkerCount - 1 - ordinal1;
    
    uint64_t position = orientedReadMarkers0[ordinal0].position;
    
    Kmer kmer0;
    for(uint64_t i=0; i<k; i++) {
        kmer0.set(i, read[position + i]);
    }
    return kmer0.reverseComplement(k);
}


vector<Base> Shasta2Anchors::anchorKmerSequence(Shasta2AnchorId anchorId) const
{
    vector<Base> sequence;
    const Shasta2Anchor anchor = (*this)[anchorId];
    if(anchor.empty()) {
        return sequence;
    }

    const Shasta2AnchorMarkerInfo& markerInfo = anchor.front();
    const Kmer kmerValue = getKmer(markerInfo.orientedReadId, markerInfo.ordinal);
    sequence.reserve(k);
    for(uint64_t i=0; i<k; i++) {
        sequence.push_back(kmerValue[i]);
    }
    return sequence;
}


Kmer Shasta2Anchors::anchorKmer(Shasta2AnchorId anchorId) const
{
    const Shasta2Anchor anchor = (*this)[anchorId];
    if(anchor.empty()) {
        return Kmer();
    }
    const Shasta2AnchorMarkerInfo& markerInfo = anchor.front();
    return getKmer(markerInfo.orientedReadId, markerInfo.ordinal);
}


// findChildren
void Shasta2Anchors::findChildren(
    const Shasta2Journeys& journeys,
    Shasta2AnchorId anchorId,
    vector<Shasta2AnchorId>& children,
    vector<uint64_t>& count,
    uint64_t minCoverage) const
{
    children.clear();
    for(const auto& markerInfo: (*this)[anchorId]) {
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const auto journey = journeys[orientedReadId];
        const uint64_t position = markerInfo.positionInJourney;
        const uint64_t nextPosition = position + 1;
        if(nextPosition < journey.size()) {
            const Shasta2AnchorId nextAnchorId = journey[nextPosition];
            children.push_back(nextAnchorId);
        }
    }
    deduplicateAndCountWithThreshold(children, count, minCoverage);
}

// findParents
void Shasta2Anchors::findParents(
    const Shasta2Journeys& journeys,
    Shasta2AnchorId anchorId,
    vector<Shasta2AnchorId>& parents,
    vector<uint64_t>& count,
    uint64_t minCoverage) const
{
    parents.clear();
    for(const auto& markerInfo: (*this)[anchorId]) {
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const auto journey = journeys[orientedReadId];
        const uint64_t position = markerInfo.positionInJourney;
        if(position > 0) {
            const uint64_t previousPosition = position - 1;
            const Shasta2AnchorId previousAnchorId = journey[previousPosition];
            parents.push_back(previousAnchorId);
        }
    }
    deduplicateAndCountWithThreshold(parents, count, minCoverage);
}


void Shasta2Anchors::analyzeAnchorPair(
    Shasta2AnchorId anchorId0, Shasta2AnchorId anchorId1, Shasta2AnchorPairInfo& info) const 
{
    info.totalA = (*this)[anchorId0].size();
    info.totalB = (*this)[anchorId1].size();
    info.common = countCommon(anchorId0, anchorId1);
}

uint64_t Shasta2Anchors::countCommon(Shasta2AnchorId anchorId0, Shasta2AnchorId anchorId1) const
{
    const Shasta2Anchor anchor0 = (*this)[anchorId0];
    const Shasta2Anchor anchor1 = (*this)[anchorId1];

    auto it0 = anchor0.begin();
    auto it1 = anchor1.begin();
    const auto end0 = anchor0.end();
    const auto end1 = anchor1.end();

    uint64_t count = 0;
    while((it0 != end0) and (it1 != end1)) {
        const OrientedReadId orientedReadId0 = it0->orientedReadId;
        const OrientedReadId orientedReadId1 = it1->orientedReadId;
        if(orientedReadId0 < orientedReadId1) {
            ++it0;
        } else if(orientedReadId1 < orientedReadId0) {
            ++it1;
        } else {
             // Same read.
             // Only count if ordinal on 1 > ordinal on 0 (forward flow).
            if(it0->ordinal < it1->ordinal) {
                ++count;
            }
            ++it0;
            ++it1;
        }
    }
    return count; 
}



uint32_t Shasta2Anchors::getOrdinal(Shasta2AnchorId anchorId, OrientedReadId orientedReadId) const {
    const Shasta2Anchor anchor = (*this)[anchorId];
    for(const auto& info: anchor) {
        if(info.orientedReadId == orientedReadId) {
            return info.ordinal;
        }
    }
    return invalid<uint32_t>;
}

uint32_t Shasta2Anchors::getPositionInJourney(Shasta2AnchorId anchorId, OrientedReadId orientedReadId) const {
    const Shasta2Anchor anchor = (*this)[anchorId];
    for(const auto& info: anchor) {
        if(info.orientedReadId == orientedReadId) {
            return info.positionInJourney;
        }
    }
    return invalid<uint32_t>;
}

const Shasta2AnchorMarkerInfo& Shasta2Anchors::getAnchorMarkerInfo(Shasta2AnchorId anchorId, OrientedReadId orientedReadId) const {
    const Shasta2Anchor anchor = (*this)[anchorId];
    for(const auto& info: anchor) {
        if(info.orientedReadId == orientedReadId) {
            return info;
        }
    }
    throw runtime_error("Anchor marker info not found.");
}

bool Shasta2Anchors::anchorContains(Shasta2AnchorId anchorId, OrientedReadId orientedReadId) const {
    const Shasta2Anchor anchor = (*this)[anchorId];
    for(const auto& info: anchor) {
        if(info.orientedReadId == orientedReadId) {
            return true;
        }
    }
    return false;
}
