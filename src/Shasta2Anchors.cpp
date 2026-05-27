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
#include <unordered_map>
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

    class ExternalAnchorOrientedRead {
    public:
        OrientedReadId orientedReadId;
        uint32_t position;

        ExternalAnchorOrientedRead() {}

        ExternalAnchorOrientedRead(
            OrientedReadId orientedReadId,
            uint32_t position) :
            orientedReadId(orientedReadId),
            position(position)
        {}
    };

}

string dinara::shasta2AnchorIdToString(Shasta2AnchorId anchorId)
{
    return std::to_string(anchorId);
}

dinara::Shasta2AnchorId dinara::shasta2AnchorIdFromString(const string& s)
{
    if(s.empty()) {
        return invalid<Shasta2AnchorId>;
    }

    string t = s;
    if(t.back() == '+' || t.back() == '-') {
        t.pop_back();
    }
    try {
        size_t used = 0;
        const uint64_t x = std::stoull(t, &used);
        if(used != t.size()) {
            return invalid<Shasta2AnchorId>;
        }
        return Shasta2AnchorId(x);
    } catch(...) {
        return invalid<Shasta2AnchorId>;
    }
}


// Minimal constructor: initializes members and builds MarkerKmers only.
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
}


Shasta2Anchors::Shasta2Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,

    const MarkerGraph& markerGraph,
    uint64_t threadCount,
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage) :
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

    // Build anchors in canonical/RC pairs so that anchor 2*i is the canonical
    // (forward) anchor and 2*i+1 is its reverse complement. This matches
    // shasta2's anchor ID scheme from readExternalAnchors.
    //
    // Step 1: Select canonical vertices that pass the coverage filter.
    // A vertex is canonical if vertexId <= reverseComplementVertex[vertexId].
    const uint64_t vertexCount = markerGraph.vertexCount();
    const bool hasRcVertex = markerGraph.reverseComplementVertex.isOpen;
    DINARA_ASSERT(hasRcVertex);

    auto& data = constructData;
    data.minAnchorCoverage = minAnchorCoverage;
    data.maxAnchorCoverage = maxAnchorCoverage;

    // Collect canonical vertex IDs that pass coverage filter.
    vector<MarkerGraphVertexId> canonicalVertexIds;
    canonicalVertexIds.reserve(vertexCount / 2);

    vector<OrientedReadId> orientedReadIds;
    for(MarkerGraphVertexId vertexId=0; vertexId<vertexCount; ++vertexId) {
        const MarkerGraphVertexId rcVertexId = markerGraph.reverseComplementVertex[vertexId];
        if(rcVertexId < vertexId) {
            continue; // Not canonical — skip.
        }

        orientedReadIds.clear();
        for(const MarkerId markerId : markerGraph.getVertexMarkerIds(vertexId)) {
            OrientedReadId orientedReadId;
            uint32_t ordinal;
            tie(orientedReadId, ordinal) = findMarkerId(markerId, markers);
            (void)ordinal;
            orientedReadIds.push_back(orientedReadId);
        }
        std::sort(orientedReadIds.begin(), orientedReadIds.end());
        auto last = std::unique(orientedReadIds.begin(), orientedReadIds.end());
        const uint64_t coverage = uint64_t(std::distance(orientedReadIds.begin(), last));
        if(coverage >= minAnchorCoverage && coverage <= maxAnchorCoverage) {
            canonicalVertexIds.push_back(vertexId);
        }
    }

    const uint64_t canonicalCount = canonicalVertexIds.size();
    cout << "Selected " << canonicalCount << " canonical Shasta2 anchors from "
         << vertexCount << " markerGraph vertices using coverage range ["
         << minAnchorCoverage << ", " << maxAnchorCoverage << "]." << endl;

    // Step 2: Build selectedVertexIds with paired layout:
    // anchor 2*i = canonicalVertexIds[i], anchor 2*i+1 = RC of canonicalVertexIds[i].
    data.selectedVertexIds.clear();
    data.selectedVertexIds.reserve(2 * canonicalCount);
    for(const MarkerGraphVertexId vertexId : canonicalVertexIds) {
        const MarkerGraphVertexId rcVertexId = markerGraph.reverseComplementVertex[vertexId];
        data.selectedVertexIds.push_back(vertexId);    // 2*i: canonical (forward)
        data.selectedVertexIds.push_back(rcVertexId);   // 2*i+1: RC
    }

    cout << "Created " << data.selectedVertexIds.size() << " Shasta2 anchors ("
         << canonicalCount << " canonical/RC pairs)." << endl;

    anchorMarkerInfos.createNew(
        largeDataName("Shasta2Anchors-AnchorMarkerInfos"),
        largeDataPageSize);
    
    // Initialize Pass 1 for VectorOfVectors
    anchorMarkerInfos.beginPass1(data.selectedVertexIds.size());
    
    const uint64_t batchSize = 1000;

    // Run Pass 1 in parallel to count entries.
    setupLoadBalancing(data.selectedVertexIds.size(), batchSize);
    runThreads(&Shasta2Anchors::constructThreadFunctionPass1, threadCount);
    
    // Initialize Pass 2
    anchorMarkerInfos.beginPass2();
    
    // Run Pass 2 in parallel to fill entries.
    setupLoadBalancing(data.selectedVertexIds.size(), batchSize);
    runThreads(&Shasta2Anchors::constructThreadFunctionPass2, threadCount);
    
    // Finalize VectorOfVectors
    anchorMarkerInfos.endPass2();
    anchorVertexIds = std::move(data.selectedVertexIds);
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
        for(uint64_t anchorId=begin; anchorId!=end; anchorId++) {
            const MarkerGraphVertexId vertexId = constructData.selectedVertexIds[anchorId];
            
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
            
            anchorMarkerInfos.incrementCount(anchorId, count);
        }
    }
}

void Shasta2Anchors::constructThreadFunctionPass2(uint64_t threadId)
{
    uint64_t begin, end;
    
    vector<Shasta2AnchorMarkerInfo> buffer;

    while(getNextBatch(begin, end)) {
        for(uint64_t anchorId=begin; anchorId!=end; anchorId++) {
            const MarkerGraphVertexId vertexId = constructData.selectedVertexIds[anchorId];
            
            buffer.clear();
            
            // Add markers from this vertex only.
            // Store both the midpoint position (marker position + k/2, matching
            // shasta2's MarkerInfo::position) and the ordinal (for internal use).
            for(const MarkerId markerId : markerGraph.getVertexMarkerIds(vertexId)) {
                OrientedReadId orientedReadId;
                uint32_t ordinal;
                tie(orientedReadId, ordinal) = findMarkerId(markerId, markers);
                const uint32_t position =
                    markers[orientedReadId.getValue()][ordinal].position + uint32_t(kHalf);
                buffer.emplace_back(orientedReadId, position, ordinal);
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
                anchorMarkerInfos.store(anchorId, *it);
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
Kmer Shasta2Anchors::getKmerAtPosition(OrientedReadId orientedReadId, uint32_t midpointPosition) const {
    // The midpoint position is marker position + k/2.
    // The raw position (first base of k-mer) is midpointPosition - k/2.
    const uint32_t kHalf = uint32_t(k / 2);
    const uint32_t rawPosition = midpointPosition - kHalf;
    const auto read = reads.getRead(orientedReadId.getReadId());

    if(orientedReadId.getStrand() == 0) {
        Kmer kmer;
        for(uint64_t i = 0; i < k; i++) {
            kmer.set(i, read[rawPosition + i]);
        }
        return kmer;
    } else {
        // For strand 1, convert to strand 0 position and RC.
        const uint64_t readLength = reads.getReadRawSequenceLength(orientedReadId.getReadId());
        const uint32_t rawPosition0 = uint32_t(readLength) - rawPosition - uint32_t(k);
        Kmer kmer0;
        for(uint64_t i = 0; i < k; i++) {
            kmer0.set(i, read[rawPosition0 + i]);
        }
        return kmer0.reverseComplement(k);
    }
}

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
    const Kmer kmerValue = getKmerAtPosition(markerInfo.orientedReadId, markerInfo.position);
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
    return getKmerAtPosition(markerInfo.orientedReadId, markerInfo.position);
}


uint64_t Shasta2Anchors::writeExternalAnchors(const string& name, bool canonicalOnly) const
{
    MemoryMapped::VectorOfVectors<ExternalAnchorOrientedRead, uint64_t> data;
    MemoryMapped::VectorOfVectors<char, uint64_t> names;
    data.createNew(name, 4096);
    names.createNew(name + "-Names", 4096);

    uint64_t exportedCount = 0;

    // With paired anchor IDs (2*i = canonical, 2*i+1 = RC),
    // canonical anchors are at even indices.
    for(Shasta2AnchorId anchorId=0; anchorId<size(); ++anchorId) {
        const Shasta2Anchor anchor = (*this)[anchorId];
        if(anchor.empty()) {
            continue;
        }

        // Skip RC anchors (odd indices) when canonicalOnly is set.
        if(canonicalOnly && (uint64_t(anchorId) % 2 != 0)) {
            continue;
        }

        const Kmer expectedKmer = getKmerAtPosition(anchor.front().orientedReadId, anchor.front().position);
        vector<ReadId> readIds;
        readIds.reserve(anchor.size());
        for(const Shasta2AnchorMarkerInfo& markerInfo : anchor) {
            const Kmer kmerValue = getKmerAtPosition(markerInfo.orientedReadId, markerInfo.position);
            if(kmerValue != expectedKmer) {
                throw runtime_error(
                    "Shasta2 external-anchor export failed: anchor " +
                    shasta2AnchorIdToString(anchorId) +
                    " contains inconsistent marker k-mers.");
            }

            const ReadId readId = markerInfo.orientedReadId.getReadId();
            if(std::find(readIds.begin(), readIds.end(), readId) != readIds.end()) {
                throw runtime_error(
                    "Shasta2 external-anchor export failed: anchor " +
                    shasta2AnchorIdToString(anchorId) +
                    " contains the same ReadId on both strands.");
            }
            readIds.push_back(readId);
        }

        const string anchorName = "anchor-" + shasta2AnchorIdToString(anchorId);
        data.appendVector();
        names.appendVector(anchorName.begin(), anchorName.end());
        for(const Shasta2AnchorMarkerInfo& markerInfo : anchor) {
            // External anchors store the raw position (first base of k-mer).
            const uint32_t rawPosition = markerInfo.position - uint32_t(k / 2);
            data.append(ExternalAnchorOrientedRead(markerInfo.orientedReadId, rawPosition));
        }
        ++exportedCount;
    }

    return exportedCount;
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
    const Shasta2Anchor anchor0 = (*this)[anchorId0];
    const Shasta2Anchor anchor1 = (*this)[anchorId1];

    const auto begin0 = anchor0.begin();
    const auto begin1 = anchor1.begin();
    const auto end0 = anchor0.end();
    const auto end1 = anchor1.end();

    info.totalA = end0 - begin0;
    info.totalB = end1 - begin1;

    info.common = 0;
    int64_t sumBaseOffsets = 0;

    auto it0 = begin0;
    auto it1 = begin1;
    while(it0 != end0 && it1 != end1) {
        if(it0->orientedReadId < it1->orientedReadId) {
            ++it0;
            continue;
        }
        if(it1->orientedReadId < it0->orientedReadId) {
            ++it1;
            continue;
        }

        ++info.common;

        const int64_t basePosition0 = int64_t(it0->position);
        const int64_t basePosition1 = int64_t(it1->position);
        sumBaseOffsets += basePosition1 - basePosition0;

        ++it0;
        ++it1;
    }

    info.onlyA = info.totalA - info.common;
    info.onlyB = info.totalB - info.common;

    if(info.common == 0) {
        info.offsetInMarkers = invalid<int64_t>;
        info.offsetInBases = invalid<int64_t>;
        info.onlyAShort = invalid<uint64_t>;
        info.onlyBShort = invalid<uint64_t>;
        return;
    }

    info.offsetInMarkers = invalid<int64_t>;  // Ordinals no longer stored.
    info.offsetInBases = int64_t(std::llround(double(sumBaseOffsets) / double(info.common)));

    it0 = begin0;
    it1 = begin1;
    uint64_t onlyACheck = 0;
    uint64_t onlyBCheck = 0;
    info.onlyAShort = 0;
    info.onlyBShort = 0;
    while(true) {
        if(it0 == end0 && it1 == end1) {
            break;
        }

        if(it1 == end1 || ((it0 != end0) && (it0->orientedReadId < it1->orientedReadId))) {
            ++onlyACheck;
            const OrientedReadId orientedReadId = it0->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadRawSequenceLength(orientedReadId.getReadId()));

            const int64_t basePosition0 = int64_t(it0->position);
            const int64_t hypotheticalPosition1 = basePosition0 + info.offsetInBases;
            if(hypotheticalPosition1 < 0 || hypotheticalPosition1 >= lengthInBases) {
                ++info.onlyAShort;
            }

            ++it0;
            continue;
        }

        if(it0 == end0 || ((it1 != end1) && (it1->orientedReadId < it0->orientedReadId))) {
            ++onlyBCheck;
            const OrientedReadId orientedReadId = it1->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadRawSequenceLength(orientedReadId.getReadId()));

            const int64_t basePosition1 = int64_t(it1->position);
            const int64_t hypotheticalPosition0 = basePosition1 - info.offsetInBases;
            if(hypotheticalPosition0 < 0 || hypotheticalPosition0 >= lengthInBases) {
                ++info.onlyBShort;
            }

            ++it1;
            continue;
        }

        ++it0;
        ++it1;
    }

    DINARA_ASSERT(onlyACheck == info.onlyA);
    DINARA_ASSERT(onlyBCheck == info.onlyB);
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
             // Only count if position on 1 > position on 0 (forward flow).
            if(it0->position < it1->position) {
                ++count;
            }
            ++it0;
            ++it1;
        }
    }
    return count; 
}



uint32_t Shasta2Anchors::getPosition(Shasta2AnchorId anchorId, OrientedReadId orientedReadId) const {
    const Shasta2Anchor anchor = (*this)[anchorId];
    for(const auto& info: anchor) {
        if(info.orientedReadId == orientedReadId) {
            return info.position;
        }
    }
    return invalid<uint32_t>;
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
