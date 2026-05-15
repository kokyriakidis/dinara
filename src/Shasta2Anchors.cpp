#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "MurmurHash2.hpp"
#include "deduplicate.hpp"
#include "findMarkerId.hpp"
#include "MultithreadedObject.tpp"
#include "Reads.hpp"
#include "shasta2/ShortBaseSequence.hpp"

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

    using Shasta2StyleKmer = shasta2::ShortBaseSequence128;

    uint32_t computeShasta2HashThreshold(double markerDensity)
    {
        if(markerDensity < 0. || markerDensity > 1.) {
            throw runtime_error("Invalid marker density " +
                to_string(markerDensity) + " requested.");
        }

        const double p = 1. - std::sqrt(1. - markerDensity);
        const double hashMax = std::numeric_limits<uint32_t>::max();
        return uint32_t(std::round(double(hashMax) * p));
    }

    bool isShasta2Marker(const Shasta2StyleKmer& kmer, uint64_t k, uint32_t hashThreshold)
    {
        if(MurmurHash2(&kmer, sizeof(Shasta2StyleKmer), 267457831) < hashThreshold) {
            return true;
        }

        const Shasta2StyleKmer kmerRc = kmer.reverseComplement(k);
        return MurmurHash2(&kmerRc, sizeof(Shasta2StyleKmer), 267457831) < hashThreshold;
    }
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

    // Pass 1: compute coverage for each anchor (VertexId).
    const uint64_t vertexCount = markerGraph.vertexCount();
    auto& data = constructData;
    data.minAnchorCoverage = minAnchorCoverage;
    data.maxAnchorCoverage = maxAnchorCoverage;
    data.selectedVertexIds.clear();
    data.selectedVertexIds.reserve(vertexCount);

    vector<OrientedReadId> orientedReadIds;
    for(MarkerGraphVertexId vertexId=0; vertexId<vertexCount; ++vertexId) {
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
            data.selectedVertexIds.push_back(vertexId);
        }
    }

    cout << "Selected " << data.selectedVertexIds.size() << " Shasta2 anchors from "
         << vertexCount << " markerGraph vertices using coverage range ["
         << minAnchorCoverage << ", " << maxAnchorCoverage << "]." << endl;

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


uint64_t Shasta2Anchors::filterByShasta2HashedKmerChecker(double markerDensity)
{
    const uint32_t hashThreshold = computeShasta2HashThreshold(markerDensity);
    const uint64_t anchorCount = size();
    vector<vector<Shasta2AnchorMarkerInfo> > keptAnchors;
    keptAnchors.reserve(anchorCount);

    uint64_t removedCount = 0;
    for(Shasta2AnchorId anchorId=0; anchorId<anchorCount; ++anchorId) {
        const Shasta2Anchor anchor = (*this)[anchorId];
        if(anchor.empty()) {
            ++removedCount;
            continue;
        }

        const uint32_t position0 =
            markers[anchor.front().orientedReadId.getValue()][anchor.front().ordinal].position;
        Shasta2StyleKmer shasta2Kmer;
        for(uint64_t i=0; i<k; i++) {
            shasta2Kmer.set(
                i,
                shasta2::Base::fromInteger(
                    reads.getOrientedReadBase(
                        anchor.front().orientedReadId,
                        uint32_t(position0 + i)).value));
        }
        if(!isShasta2Marker(shasta2Kmer, k, hashThreshold)) {
            ++removedCount;
            continue;
        }

        keptAnchors.emplace_back(anchor.begin(), anchor.end());
    }

    cout << "Shasta2 hashed k-mer checker retained " << keptAnchors.size() << " / "
         << anchorCount << " Shasta2 anchors." << endl;
    if(removedCount == 0) {
        return 0;
    }

    vector<MarkerGraphVertexId> keptAnchorVertexIds;
    if(!anchorVertexIds.empty()) {
        keptAnchorVertexIds.reserve(keptAnchors.size());
        for(Shasta2AnchorId anchorId=0; anchorId<anchorCount; ++anchorId) {
            const Shasta2Anchor anchor = (*this)[anchorId];
            if(anchor.empty()) {
                continue;
            }

            const uint32_t position0 =
                markers[anchor.front().orientedReadId.getValue()][anchor.front().ordinal].position;
            Shasta2StyleKmer shasta2Kmer;
            for(uint64_t i=0; i<k; i++) {
                shasta2Kmer.set(
                    i,
                    shasta2::Base::fromInteger(
                        reads.getOrientedReadBase(
                            anchor.front().orientedReadId,
                            uint32_t(position0 + i)).value));
            }
            if(isShasta2Marker(shasta2Kmer, k, hashThreshold)) {
                keptAnchorVertexIds.push_back(anchorVertexIds[anchorId]);
            }
        }
    }

    anchorMarkerInfos.clear();
    for(const auto& anchor : keptAnchors) {
        anchorMarkerInfos.appendVector(anchor);
    }
    if(!anchorVertexIds.empty()) {
        anchorVertexIds = std::move(keptAnchorVertexIds);
    }

    return removedCount;
}


uint64_t Shasta2Anchors::writeExternalAnchors(const string& name, bool canonicalOnly) const
{
    MemoryMapped::VectorOfVectors<ExternalAnchorOrientedRead, uint64_t> data;
    MemoryMapped::VectorOfVectors<char, uint64_t> names;
    data.createNew(name, 4096);
    names.createNew(name + "-Names", 4096);

    uint64_t exportedCount = 0;

    for(Shasta2AnchorId anchorId=0; anchorId<size(); ++anchorId) {
        const Shasta2Anchor anchor = (*this)[anchorId];
        if(anchor.empty()) {
            continue;
        }

        if(canonicalOnly) {
            bool skip = false;
            if(!anchorVertexIds.empty() && markerGraph.reverseComplementVertex.isOpen) {
                const MarkerGraphVertexId vertexId = anchorVertexIds[anchorId];
                const MarkerGraphVertexId vertexIdRc = markerGraph.reverseComplementVertex[vertexId];
                skip = vertexIdRc < vertexId;
            } else {
                const Kmer kmerValue = anchorKmer(anchorId);
                const Kmer kmerRc = kmerValue.reverseComplement(k);
                skip = kmerRc.id(k) < kmerValue.id(k);
            }
            if(skip) {
                continue;
            }
        }

        const Kmer expectedKmer = getKmer(anchor.front().orientedReadId, anchor.front().ordinal);
        vector<ReadId> readIds;
        readIds.reserve(anchor.size());
        for(const Shasta2AnchorMarkerInfo& markerInfo : anchor) {
            const Kmer kmerValue = getKmer(markerInfo.orientedReadId, markerInfo.ordinal);
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
            const uint32_t position = markers[markerInfo.orientedReadId.getValue()][markerInfo.ordinal].position;
            data.append(ExternalAnchorOrientedRead(markerInfo.orientedReadId, position));
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
    int64_t sumMarkerOffsets = 0;
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
        const OrientedReadId orientedReadId = it0->orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        const uint32_t ordinal0 = it0->ordinal;
        const uint32_t ordinal1 = it1->ordinal;
        sumMarkerOffsets += int64_t(ordinal1) - int64_t(ordinal0);

        const int64_t basePosition0 = int64_t(orientedReadMarkers[ordinal0].position) + int64_t(kHalf);
        const int64_t basePosition1 = int64_t(orientedReadMarkers[ordinal1].position) + int64_t(kHalf);
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

    info.offsetInMarkers = int64_t(std::llround(double(sumMarkerOffsets) / double(info.common)));
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

            const uint32_t ordinal0 = it0->ordinal;
            const int64_t basePosition0 = int64_t(orientedReadMarkers[ordinal0].position) + int64_t(kHalf);
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

            const uint32_t ordinal1 = it1->ordinal;
            const int64_t basePosition1 = int64_t(orientedReadMarkers[ordinal1].position) + int64_t(kHalf);
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
