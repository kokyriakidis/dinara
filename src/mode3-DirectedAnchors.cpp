#include "mode3-DirectedAnchors.hpp"
#include "mode3-DirectedAnchorGraph.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"

namespace dinara {
namespace mode3 {

DirectedAnchors::DirectedAnchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    uint64_t minPrimaryCoverage,
    uint64_t maxPrimaryCoverage,
    uint64_t threadCount,
    bool createFromVertices) :
    Anchors(mappedMemoryOwner, reads, k, markers, markerGraph,
            minPrimaryCoverage, maxPrimaryCoverage, threadCount, createFromVertices)
{
    computeJourneysWithPositions(threadCount);
}

void DirectedAnchors::computeJourneysWithPositions(uint64_t threadCount)
{
    performanceLog << timestamp << "DirectedAnchors::computeJourneysWithPositions begins." << endl;

    const uint64_t orientedReadCount = 2 * reads.readCount();

    // Pass 1 & 2: Collect JourneyAnchor objects (ID, start, end) for each oriented read.
    MemoryMapped::VectorOfVectors<JourneyAnchor, uint64_t> tmpJourneys;
    tmpJourneys.createNew(largeDataName("tmp-JourneysWithPositions"), largeDataPageSize);

    tmpJourneys.beginPass1(orientedReadCount);
    for(AnchorId anchorId = 0; anchorId < size(); ++anchorId) {
        for(const auto& mi : (*this)[anchorId]) {
            tmpJourneys.incrementCountMultithreaded(mi.orientedReadId.getValue());
        }
    }

    tmpJourneys.beginPass2();
    const uint64_t kHalf = k / 2;
    for(AnchorId anchorId = 0; anchorId < size(); ++anchorId) {
        const uint32_t offset = ordinalOffset(anchorId);
        for(const auto& mi : (*this)[anchorId]) {
            const auto& orientedMarkers = markers[mi.orientedReadId.getValue()];
            const uint32_t p0 = orientedMarkers[mi.ordinal0].position + kHalf;
            const uint32_t p1 = orientedMarkers[mi.ordinal0 + offset].position + kHalf;
            tmpJourneys.storeMultithreaded(mi.orientedReadId.getValue(), {anchorId, p0, p1});
        }
    }
    tmpJourneys.endPass2();

    // Pass 3: Sort by position and store in final structures.
    journeysWithPositions.createNew(largeDataName("JourneysWithPositions"), largeDataPageSize);
    journeys.createNew(largeDataName("Journeys"), largeDataPageSize); // Maintain compatibility

    journeysWithPositions.beginPass1(orientedReadCount);
    journeys.beginPass1(orientedReadCount);

    for(uint64_t i = 0; i < orientedReadCount; ++i) {
        const uint64_t sz = tmpJourneys.size(i);
        journeysWithPositions.incrementCountMultithreaded(i, sz);
        journeys.incrementCountMultithreaded(i, sz);
    }

    journeysWithPositions.beginPass2();
    journeys.beginPass2();
    for(uint64_t i = 0; i < orientedReadCount; ++i) {
        auto journey = tmpJourneys[i];
        vector<JourneyAnchor> sortedJourney(journey.begin(), journey.end());
        std::sort(sortedJourney.begin(), sortedJourney.end());
        for(const auto& ja : sortedJourney) {
            journeysWithPositions.storeMultithreaded(i, ja);
            journeys.storeMultithreaded(i, ja.anchorId);
        }
    }
    journeysWithPositions.endPass2(false, true);
    journeys.endPass2(false, true);

    tmpJourneys.remove();

    performanceLog << timestamp << "DirectedAnchors::computeJourneysWithPositions ends. Loaded "
                   << journeysWithPositions.size() << " journeys." << endl;
}

} // namespace mode3
} // namespace dinara
