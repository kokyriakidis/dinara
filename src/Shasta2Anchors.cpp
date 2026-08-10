#include "Shasta2Anchors.hpp"
#include "HetAnchorK.hpp"
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

    while(getNextBatch(begin, end)) {
        for(uint64_t anchorId=begin; anchorId!=end; anchorId++) {
            const MarkerGraphVertexId vertexId = constructData.selectedVertexIds[anchorId];
            
            const auto vertexMarkerIds = markerGraph.getVertexMarkerIds(vertexId);
            anchorMarkerInfos.incrementCount(anchorId, vertexMarkerIds.size());
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


uint64_t Shasta2Anchors::writeExternalAnchors(
    const string& name,
    bool canonicalOnly,
    const ExternalAnchorDropMap* dropMap) const
{
    MemoryMapped::VectorOfVectors<ExternalAnchorOrientedRead, uint64_t> data;
    MemoryMapped::VectorOfVectors<char, uint64_t> names;
    data.createNew(name, 4096);
    names.createNew(name + "-Names", 4096);

    uint64_t exportedCount = 0;
    uint64_t droppedMemberCount = 0;

    // Predicate: is this (canonical anchor, read) member dropped to resolve a
    // journey position tie? The drop set is keyed on the canonical (even) id.
    const auto isDropped = [&](Shasta2AnchorId anchorId, ReadId readId) -> bool {
        if(dropMap == nullptr) {
            return false;
        }
        const auto it = dropMap->find(anchorId);
        if(it == dropMap->end()) {
            return false;
        }
        return std::find(it->second.begin(), it->second.end(), readId) !=
            it->second.end();
    };

    // With paired anchor IDs (2*i = canonical, 2*i+1 = RC),
    // canonical anchors are at even indices.
    // We must not skip empty anchors: shasta2 assigns sequential IDs to
    // external anchors, so skipping would break the identity mapping
    // between dinara anchor IDs and shasta2 anchor IDs.
    for(Shasta2AnchorId anchorId=0; anchorId<size(); ++anchorId) {

        // Skip RC anchors (odd indices) when canonicalOnly is set.
        if(canonicalOnly && (uint64_t(anchorId) % 2 != 0)) {
            continue;
        }

        const Shasta2Anchor anchor = (*this)[anchorId];
        const string anchorName = "anchor-" + shasta2AnchorIdToString(anchorId);

        if(anchor.empty()) {
            // Export an empty anchor to preserve ID numbering.
            data.appendVector();
            names.appendVector(anchorName.begin(), anchorName.end());
            ++exportedCount;
            continue;
        }

        // Het anchors (k=2) are grouped by only 2 shared bases, so their
        // members' k=50 k-mers legitimately differ. Skip the k=50 consistency
        // check for them; shasta2 re-derives a 2-base k-mer from rawPosition.
        const bool isHetAnchor =
            (hetAnchorFirstId != invalid<Shasta2AnchorId>) &&
            (anchorId >= hetAnchorFirstId);

        // Only compute the k=50 reference k-mer for primary anchors. A het
        // anchor's stored position is a k=2 midpoint (rawPosition + 1), so
        // getKmerAtPosition -- which subtracts k/2 = 25 and reads k=50 bases --
        // would underflow the raw position and read out of bounds (segfault)
        // near a read end. Het anchors skip the k=50 consistency check anyway.
        const Kmer expectedKmer = isHetAnchor ? Kmer() :
            getKmerAtPosition(anchor.front().orientedReadId, anchor.front().position);
        vector<ReadId> readIds;
        readIds.reserve(anchor.size());
        uint64_t keptMembers = 0;
        for(const Shasta2AnchorMarkerInfo& markerInfo : anchor) {
            const ReadId readId = markerInfo.orientedReadId.getReadId();

            // Skip members dropped to resolve a journey position tie. Do this
            // before the k-mer / duplicate-ReadId checks: a dropped member is
            // not exported, so it need not satisfy them.
            if(isDropped(anchorId, readId)) {
                ++droppedMemberCount;
                continue;
            }

            if(!isHetAnchor) {
                const Kmer kmerValue = getKmerAtPosition(markerInfo.orientedReadId, markerInfo.position);
                if(kmerValue != expectedKmer) {
                    throw runtime_error(
                        "Shasta2 external-anchor export failed: anchor " +
                        shasta2AnchorIdToString(anchorId) +
                        " contains inconsistent marker k-mers.");
                }
            }

            if(std::find(readIds.begin(), readIds.end(), readId) != readIds.end()) {
                throw runtime_error(
                    "Shasta2 external-anchor export failed: anchor " +
                    shasta2AnchorIdToString(anchorId) +
                    " contains the same ReadId on both strands.");
            }
            readIds.push_back(readId);
            ++keptMembers;
        }

        // A drop can empty an anchor (or leave it a singleton). shasta2 accepts
        // low-coverage external anchors; still export the (possibly empty)
        // anchor to preserve the sequential id-to-id mapping shasta2 relies on.
        (void)keptMembers;

        data.appendVector();
        names.appendVector(anchorName.begin(), anchorName.end());
        // Every anchor's rawPosition is its stored midpoint minus the uniform
        // export shift, so shasta2 (which stores midpoint = rawPosition + k/2 on
        // load) recovers exactly the original midpoint for backbone, het, and
        // hom anchors alike. The recovered midpoint is the SAME in both modes;
        // only the exported raw number and the loader --k differ.
        //
        // At k=2 the exported 2-base k-mer differs by anchor class but is always
        // consistent across a given anchor's member reads:
        //  - Het/hom anchors have 2 shared bases by construction.
        //  - Primary anchors agree over the full k=50 window (verified above),
        //    so the centered 2-base subset at [position-1, position] is
        //    identical across all members too.
        // At k=0 there is no k-mer at all -- each anchor is a bare position
        // marker -- so no per-class k-mer consistency applies.
        // Export shift is UNIFORM across every anchor class and equals the k/2
        // that the shasta2 loader re-adds. shasta2 uses a single --k for the
        // whole external anchor set, so the shift cannot differ between primary
        // and het anchors without breaking their relative ordering. Default
        // k=2 -> subtract 1 (shasta2 --k 2 re-adds 1); experimental k=0 ->
        // subtract 0 (shasta2 --k 0 re-adds 0, positions land exactly on the
        // stored midpoints, which for het anchors are the exact SNP bases).
        // Primaries are k=2-clipped at k=2 and become point markers at their
        // k=50 center at k=0; either way the stored-midpoint ordering used by
        // the whole internal graph is preserved. The export k must match: see
        // the writeExternalAnchors caller, which passes --k hetAnchorK().
        const uint32_t exportShift = hetAnchorKHalf();
        for(const Shasta2AnchorMarkerInfo& markerInfo : anchor) {
            // Same drop filter as the validation loop above, so the exported
            // members exactly match the validated set.
            if(isDropped(anchorId, markerInfo.orientedReadId.getReadId())) {
                continue;
            }
            // External anchors store the raw position (first base of the k-mer).
            const uint32_t rawPosition = markerInfo.position - exportShift;
            data.append(ExternalAnchorOrientedRead(markerInfo.orientedReadId, rawPosition));
        }
        ++exportedCount;
    }

    if(droppedMemberCount > 0) {
        cout << "Dropped " << droppedMemberCount
             << " anchor member(s) to resolve journey position ties." << endl;
    }

    return exportedCount;
}


Shasta2AnchorId Shasta2Anchors::appendHetAnchorPair(
    const vector<std::pair<OrientedReadId, uint32_t>>& members)
{
    // Marker length of a het anchor. Default k=2: shasta2 re-derives [predBase,
    // alleleBase] from rawPosition; only 2 bases are meaningful. Experimental
    // k=0 (DINARA_HET_K=0): a zero-length position marker at the exact SNP base.
    // See HetAnchorK.hpp. The RC mirror formula below (readLen - rawPosition -
    // hetK, then + hetKHalf) keeps fwd+rc stored positions summing to readLen
    // for both k, so no per-k special case is needed.
    const uint32_t hetK = hetAnchorK();
    const uint32_t hetKHalf = hetAnchorKHalf();

    // Safety net: a het/hom anchor must have at least 2 members. A coverage-1
    // anchor (a single read, or only the backbone) is a spurious branch in the
    // assembly graph and must never be exported. The emit tail enforces this
    // upstream (per-arm minSupport gate; homs require >1) and the interval engine
    // now places one-sided reads by default so real allele support is not lost to
    // coverage holes. Assert it here so any future caller that regresses is
    // caught immediately rather than emitting a bad external anchor.
    DINARA_ASSERT(members.size() >= 2);

    // Record where het anchors begin so the export can bypass the k=50 k-mer
    // consistency check for them AND use hetK/2 (not k/2) when recovering the
    // raw position.
    if(hetAnchorFirstId == invalid<Shasta2AnchorId>) {
        hetAnchorFirstId = size();
    }

    // Build the canonical (forward) member list. Store the midpoint position
    // using the het marker's own half-length (hetK/2 = 1), NOT the store's k/2.
    // This keeps the het anchor's position properly ordered between its flanking
    // backbone anchors in the anchor graph (a k/2=25 offset would push it past
    // nearby anchors and create spurious backward edges). writeExternalAnchors
    // recovers rawPosition via (position - hetK/2) for het anchors.
    // ordinal/positionInJourney are not meaningful for a het anchor and are
    // left invalid.
    vector<Shasta2AnchorMarkerInfo> fwd;
    fwd.reserve(members.size());
    for(const auto& [orientedReadId, rawPosition] : members) {
        Shasta2AnchorMarkerInfo info;
        info.orientedReadId = orientedReadId;
        info.position = rawPosition + hetKHalf;
        fwd.push_back(info);
    }
    std::sort(fwd.begin(), fwd.end());

    // Build the reverse-complement member list: flip strand and mirror the raw
    // position to the opposite strand's coordinate frame, then re-apply the
    // store's midpoint convention. For a k-base marker at rawPosition on strand
    // s, the RC raw position on strand s^1 is readLen - rawPosition - hetK. This
    // is general in hetK: k=0 gives readLen - rawPosition (mirror of a boundary
    // position), keeping fwd+rc stored midpoints summing to readLen for both k.
    vector<Shasta2AnchorMarkerInfo> rc;
    rc.reserve(members.size());
    for(const auto& [orientedReadId, rawPosition] : members) {
        OrientedReadId rcOid = orientedReadId;
        rcOid.flipStrand();
        const uint64_t readLen =
            reads.getReadRawSequenceLength(orientedReadId.getReadId());
        const uint32_t rcRaw =
            uint32_t(readLen) - rawPosition - hetK;
        Shasta2AnchorMarkerInfo info;
        info.orientedReadId = rcOid;
        info.position = rcRaw + hetKHalf;
        rc.push_back(info);
    }
    std::sort(rc.begin(), rc.end());

    // Append canonical (even id) then RC (odd id), preserving the store's
    // 2i/2i+1 pairing. append() pushes to the end of the last vector, so
    // forward iteration keeps the already-sorted ascending OrientedReadId order.
    const Shasta2AnchorId canonicalId = size();
    anchorMarkerInfos.appendVector();
    for(const auto& info : fwd) {
        anchorMarkerInfos.append(info);
    }
    anchorMarkerInfos.appendVector();
    for(const auto& info : rc) {
        anchorMarkerInfos.append(info);
    }

    return canonicalId;
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
        // A read whose filtered journey no longer includes this anchor has
        // positionInJourney == invalid here; it does not transition out of this
        // anchor, so skip it (guards against invalid+1 == 0 aliasing to
        // journey[0], which would add a spurious, non-RC-symmetric child).
        if(position == invalid<uint32_t>) {
            continue;
        }
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
        // A read whose filtered journey no longer includes this anchor has
        // positionInJourney == invalid here; it does not transition into this
        // anchor, so skip it (guards against invalid-1 indexing journey out of
        // bounds and producing a non-RC-symmetric parent).
        if(position == invalid<uint32_t>) {
            continue;
        }
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



namespace {
    // Every anchor's member list is stored sorted ascending by OrientedReadId --
    // primary anchors are explicitly sorted before storage (see the "Sort by
    // OrientedReadId to ensure canonical order" block above), and appendHetAnchorPair
    // sorts both the forward and RC member lists the same way before appending.
    // So a lookup by OrientedReadId can binary search instead of scanning the
    // whole anchor, turning an O(coverage) scan into O(log coverage). This is the
    // single lookup primitive behind getPosition/getOrdinal/getPositionInJourney/
    // getAnchorMarkerInfo/anchorContains, all called in per-anchor hot loops
    // throughout window creation, window transitions, and anchor graph construction.
    const Shasta2AnchorMarkerInfo* findMarkerInfo(
        const Shasta2Anchor& anchor, OrientedReadId orientedReadId)
    {
        const auto it = std::lower_bound(anchor.begin(), anchor.end(), orientedReadId,
            [](const Shasta2AnchorMarkerInfo& info, OrientedReadId oid) {
                return info.orientedReadId < oid;
            });
        if(it != anchor.end() and it->orientedReadId == orientedReadId) {
            return &(*it);
        }
        return nullptr;
    }
}

uint32_t Shasta2Anchors::getPosition(Shasta2AnchorId anchorId, OrientedReadId orientedReadId) const {
    const Shasta2AnchorMarkerInfo* info = findMarkerInfo((*this)[anchorId], orientedReadId);
    return info ? info->position : invalid<uint32_t>;
}

uint32_t Shasta2Anchors::getOrdinal(Shasta2AnchorId anchorId, OrientedReadId orientedReadId) const {
    const Shasta2AnchorMarkerInfo* info = findMarkerInfo((*this)[anchorId], orientedReadId);
    return info ? info->ordinal : invalid<uint32_t>;
}

uint32_t Shasta2Anchors::getPositionInJourney(Shasta2AnchorId anchorId, OrientedReadId orientedReadId) const {
    const Shasta2AnchorMarkerInfo* info = findMarkerInfo((*this)[anchorId], orientedReadId);
    return info ? info->positionInJourney : invalid<uint32_t>;
}

const Shasta2AnchorMarkerInfo& Shasta2Anchors::getAnchorMarkerInfo(Shasta2AnchorId anchorId, OrientedReadId orientedReadId) const {
    const Shasta2AnchorMarkerInfo* info = findMarkerInfo((*this)[anchorId], orientedReadId);
    if(info == nullptr) {
        throw runtime_error("Anchor marker info not found.");
    }
    return *info;
}

bool Shasta2Anchors::anchorContains(Shasta2AnchorId anchorId, OrientedReadId orientedReadId) const {
    return findMarkerInfo((*this)[anchorId], orientedReadId) != nullptr;
}
