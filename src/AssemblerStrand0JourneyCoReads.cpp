// Per-strand-0-read journey co-read CSR: structural input before het detection / phasing.

#include "Assembler.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"

#include "algorithm.hpp"
#include "chrono.hpp"
#include "iostream.hpp"
#include <unordered_set>

#include "vector.hpp"

using namespace dinara;
using namespace std;

void Assembler::computeStrand0JourneyCoReadsTable()
{
    reads->checkReadsAreOpen();
    DINARA_ASSERT(shasta2Anchors);
    DINARA_ASSERT(shasta2Journeys);
    DINARA_ASSERT(shasta2Journeys->isOpen());

    if(strand0JourneyCoReads.isOpen()) {
        strand0JourneyCoReads.remove();
    }

    const uint64_t readCount = reads->readCount();
    cout << timestamp << "Computing strand0JourneyCoReads: reads=" << readCount << endl;
    const auto t0 = chrono::steady_clock::now();

    strand0JourneyCoReads.createNew(
        largeDataName("Strand0JourneyCoReads"),
        largeDataPageSize);
    strand0JourneyCoReads.beginPass1(readCount);

    unordered_set<uint32_t> unique;
    vector<uint32_t> sortedPartners;
    for(ReadId readId = 0; readId < readCount; readId++) {
        const OrientedReadId focal(readId, 0);
        const auto journey = (*shasta2Journeys)[focal];
        unique.clear();
        unique.reserve(min<size_t>(size_t(1) << 20, journey.size() * 8 + 16));
        for(const Shasta2AnchorId anchorId: journey) {
            const Shasta2Anchor anchor = (*shasta2Anchors)[anchorId];
            for(const Shasta2AnchorMarkerInfo& info: anchor) {
                if(info.orientedReadId != focal) {
                    unique.insert(uint32_t(info.orientedReadId.getValue()));
                }
            }
        }
        strand0JourneyCoReads.incrementCount(uint64_t(readId), uint64_t(unique.size()));
    }

    strand0JourneyCoReads.beginPass2();
    for(ReadId readId = 0; readId < readCount; readId++) {
        const OrientedReadId focal(readId, 0);
        const auto journey = (*shasta2Journeys)[focal];
        unique.clear();
        unique.reserve(min<size_t>(size_t(1) << 20, journey.size() * 8 + 16));
        for(const Shasta2AnchorId anchorId: journey) {
            const Shasta2Anchor anchor = (*shasta2Anchors)[anchorId];
            for(const Shasta2AnchorMarkerInfo& info: anchor) {
                if(info.orientedReadId != focal) {
                    unique.insert(uint32_t(info.orientedReadId.getValue()));
                }
            }
        }
        sortedPartners.assign(unique.begin(), unique.end());
        sort(sortedPartners.begin(), sortedPartners.end());
        for(const uint32_t v: sortedPartners) {
            strand0JourneyCoReads.store(uint64_t(readId), v);
        }
    }
    strand0JourneyCoReads.endPass2();
    strand0JourneyCoReads.unreserve();

    const auto t1 = chrono::steady_clock::now();
    const double seconds = chrono::duration<double>(t1 - t0).count();
    cout << timestamp << "strand0JourneyCoReads done in " << seconds << " s." << endl;
}
