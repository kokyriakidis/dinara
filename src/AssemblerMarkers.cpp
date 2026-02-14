// Dinara.
#include "Assembler.hpp"
#include "Align6Marker.hpp"
#include "extractKmer.hpp"
#include "findMarkerId.hpp"
#include "KmerCounter.hpp"
#include "KmerDistributionInfo.hpp"
#include "KmerChecker.hpp"
#include "MarkerFinder.hpp"
#include "MarkerKmers.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
using namespace dinara;

// SIMD minimizers library.
#include <simd-minimizers/simd_minimizers.h>

// Standard library.
#include "fstream.hpp"
#include <algorithm>
#include <chrono>
#include <sstream>
#include <vector>


void Assembler::findMarkers(uint64_t threadCount)
{
    reads->checkReadsAreOpen();
    DINARA_ASSERT(kmerChecker);

    markers->createNew(largeDataName("Markers"), largeDataPageSize);
    MarkerFinder markerFinder(
        assemblerInfo->k,
        *kmerChecker,
        getReads(),
        *markers,
        threadCount);

}



// Helper: compute deduplicated canonical closed syncmer markers for a read.
// The sketcher must be configured with (k_scan, s) where k_scan = k for odd k, or k+1 for even k.
// If validMarkers is provided, it is filled with (position, kmerId).
// If validMarkers is nullptr, it returns the count of markers without filling a buffer.
static size_t getSyncmerMarkersForRead(
    ReadId readId,
    const Reads& reads,
    int k,
    SimdSketcher* sketcher,
    const shared_ptr<KmerChecker>& kmerChecker,
    string& readSequence,
    std::vector<uint32_t>& positionBuffer,
    std::vector<std::pair<uint32_t, KmerId>>* validMarkers = nullptr
) {
    if(validMarkers) validMarkers->clear();
    const LongBaseSequenceView read = reads.getRead(readId);
    const uint64_t baseCount = read.baseCount;

    // k_scan is the k used by the sketcher (odd, as required by the library).
    const int k_scan = (k % 2 == 0) ? k + 1 : k;

    if(baseCount < uint64_t(k_scan)) {
        return 0;
    }

    // Convert read to string.
    readSequence.resize(baseCount);
    for(uint64_t i = 0; i < baseCount; i++) {
        readSequence[i] = read[i].character();
    }

    // Compute canonical closed syncmer positions and merge with forced read ends (0 and baseCount - k_scan).
    // Read end forcing ensures that the very first and last possible k-mers of every read
    // are included as markers, which is critical for graph connectivity at read boundaries.
    SyncmerList syncmerList = canonical_syncmer_positions(
        sketcher, readSequence.c_str(), readSequence.size());
    
    positionBuffer.assign(syncmerList.data, syncmerList.data + syncmerList.len);
    positionBuffer.push_back(0);
    positionBuffer.push_back(uint32_t(baseCount - k_scan));
    free_syncmer_list(syncmerList);

    // Sort and deduplicate to handle cases where read ends were already selected as syncmers.
    std::sort(positionBuffer.begin(), positionBuffer.end());
    positionBuffer.erase(std::unique(positionBuffer.begin(), positionBuffer.end()), positionBuffer.end());

    const size_t uniqueCandidateCount = positionBuffer.size();
    
    // Pass 1 optimization: if kmerChecker is null (initial discovery) and we only need the count,
    // we can return uniqueCandidateCount directly for odd k.
    if(!validMarkers && !kmerChecker && k_scan == k) {
        return uniqueCandidateCount;
    }

    if(validMarkers) validMarkers->reserve(uniqueCandidateCount);
    size_t validCount = 0;

    if(k_scan == k) {
        // Odd k: positions are directly usable as k-mer start positions.
        for(size_t i = 0; i < uniqueCandidateCount; i++) {
            const uint32_t position = positionBuffer[i];
            if(uint64_t(position) + uint64_t(k) > baseCount) continue;

            if(!kmerChecker) {
                // Initial discovery: everything is a marker.
                if(validMarkers) {
                    Kmer kmer;
                    extractKmer(read, uint64_t(position), uint64_t(k), kmer);
                    validMarkers->push_back({position, kmer.id(uint64_t(k))});
                }
                validCount++;
            } else {
                // Filtering mode: check if k-mer is a known marker.
                Kmer kmer;
                extractKmer(read, uint64_t(position), uint64_t(k), kmer);
                const KmerId kmerId = kmer.id(uint64_t(k));
                if(kmerChecker->isMarker(kmerId)) {
                    if(validMarkers) validMarkers->push_back({position, kmerId});
                    validCount++;
                }
            }
        }
    } else {
        // Even k: positions are (k+1)-mer positions. Adjust to k-mer positions
        // based on canonicality of each (k+1)-mer.
        // Since k_scan is odd, palindromes are impossible.
        uint32_t lastPosition = UINT32_MAX;
        for(size_t i = 0; i < uniqueCandidateCount; i++) {
            uint32_t position = positionBuffer[i];
            if(uint64_t(position) + uint64_t(k_scan) > baseCount) continue;

            // Efficient string-based canonicality check for (k+1)-mer.
            // If canonical (forward < RC): take the prefix (position P). 
            // If non-canonical: take the suffix (position P+1).
            bool isCanonical = true;
            for(int j=0; j < k_scan; ++j) {
                const char f = readSequence[position + j];
                const char r = readSequence[position + k_scan - 1 - j];
                const char rc = (r == 'A') ? 'T' : (r == 'T') ? 'A' : (r == 'C') ? 'G' : 'C';
                if(f < rc) { isCanonical = true; break; }
                if(f > rc) { isCanonical = false; break; }
            }
            if(!isCanonical) position++;

            // Skip duplicate positions from the adjustment.
            if(position == lastPosition) continue;
            lastPosition = position;

            if(uint64_t(position) + uint64_t(k) > baseCount) continue;

            if(!kmerChecker) {
                if(validMarkers) {
                    Kmer kmer;
                    extractKmer(read, uint64_t(position), uint64_t(k), kmer);
                    validMarkers->push_back({position, kmer.id(uint64_t(k))});
                }
                validCount++;
            } else {
                Kmer kmer;
                extractKmer(read, uint64_t(position), uint64_t(k), kmer);
                const KmerId kmerId = kmer.id(uint64_t(k));
                if(kmerChecker->isMarker(kmerId)) {
                    if(validMarkers) validMarkers->push_back({position, kmerId});
                    validCount++;
                }
            }
        }
    }
    return validCount;
}

void Assembler::findMarkersSimdClosedSyncmers(uint64_t threadCount, int k, int s)
{
    reads->checkReadsAreOpen();

    // The library requires odd k for canonical closed syncmers.
    // For even k, we scan with k+1 (odd) and then reduce to k-mers.
    const int k_scan = (k % 2 == 0) ? k + 1 : k;

    performanceLog << timestamp << "Finding markers using SIMD closed syncmers (k=" << k
        << ", k_scan=" << k_scan << ", s=" << s << ") in "
        << reads->readCount() << " reads." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // Store parameters.
    assemblerInfo->k = k;
    findMarkersSimdClosedSyncmersData.k = k_scan;
    findMarkersSimdClosedSyncmersData.w = s;

    // Create the markers and markerKmerIds data structures.
    markers->createNew(largeDataName("Markers"), largeDataPageSize);
    markerKmerIds->createNew(largeDataName("MarkerKmerIds"), largeDataPageSize);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t batchSize = 100;

    // Pass 1: Count markers for each oriented read.
    markers->beginPass1(2 * readCount);
    markerKmerIds->beginPass1(2 * readCount);
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::findMarkersSimdClosedSyncmersPass1, threadCount);

    // Pass 2: Store markers.
    markers->beginPass2();
    markerKmerIds->beginPass2();
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::findMarkersSimdClosedSyncmersPass2, threadCount);

    markers->endPass2(false);
    markerKmerIds->endPass2(false);

    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    performanceLog << timestamp << "Finding markers using SIMD closed syncmers completed in " << tTotal << " s." << endl;
    cout << "Created " << markers->totalSize() << " markers using SIMD closed syncmers." << endl;
}

void Assembler::findMarkersSimdClosedSyncmersPass1(size_t /* threadId */)
{
    const int k = assemblerInfo->k;
    const int k_scan = findMarkersSimdClosedSyncmersData.k;
    const int s = findMarkersSimdClosedSyncmersData.w;
    SimdSketcher* sketcher = simd_sketcher_new(
        static_cast<uint8_t>(k_scan), static_cast<uint8_t>(s));
    string readSequence;
    std::vector<uint32_t> positionBuffer;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            const size_t count = getSyncmerMarkersForRead(
                readId, *reads, k, sketcher, kmerChecker, readSequence, positionBuffer, nullptr);

            this->markers->incrementCount(OrientedReadId(readId, 0).getValue(), count);
            this->markers->incrementCount(OrientedReadId(readId, 1).getValue(), count);
            markerKmerIds->incrementCount(OrientedReadId(readId, 0).getValue(), count);
            markerKmerIds->incrementCount(OrientedReadId(readId, 1).getValue(), count);
        }
    }
    simd_sketcher_free(sketcher);
}

void Assembler::findMarkersSimdClosedSyncmersPass2(size_t /* threadId */)
{
    const int k = assemblerInfo->k;
    const int k_scan = findMarkersSimdClosedSyncmersData.k;
    const int s = findMarkersSimdClosedSyncmersData.w;
    SimdSketcher* sketcher = simd_sketcher_new(
        static_cast<uint8_t>(k_scan), static_cast<uint8_t>(s));
    string readSequence;
    std::vector<uint32_t> positionBuffer;
    std::vector<std::pair<uint32_t, KmerId>> markerBuffer;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            const LongBaseSequenceView read = reads->getRead(readId);
            getSyncmerMarkersForRead(
                readId, *reads, k, sketcher, kmerChecker, readSequence, positionBuffer, &markerBuffer);

            if(markerBuffer.empty()) continue;

            CompressedMarker* markerPointerStrand0 = this->markers->begin(OrientedReadId(readId, 0).getValue());
            CompressedMarker* markerPointerStrand1 = this->markers->end(OrientedReadId(readId, 1).getValue()) - 1;

            KmerId* kmerIdPointerStrand0 = markerKmerIds->begin(OrientedReadId(readId, 0).getValue());
            KmerId* kmerIdPointerStrand1 = markerKmerIds->end(OrientedReadId(readId, 1).getValue()) - 1;

            for(const auto& [position, kmerId] : markerBuffer) {
                // Strand 0.
                markerPointerStrand0->position = position;
                ++markerPointerStrand0;
                *kmerIdPointerStrand0 = kmerId;
                ++kmerIdPointerStrand0;

                // Strand 1: reverse complement.
                Kmer kmer(kmerId, k);
                markerPointerStrand1->position = static_cast<uint32_t>(read.baseCount - k - position);
                --markerPointerStrand1;
                *kmerIdPointerStrand1 = kmer.reverseComplement(k).id(k);
                --kmerIdPointerStrand1;
            }
        }
    }
    simd_sketcher_free(sketcher);
}



// ============================================================================
// SIMD MINIMIZERS IMPLEMENTATION
// ============================================================================

// Helper: compute deduplicated canonical minimizer positions and their KmerIds for a read.
// If validMarkers is provided, it is filled with (position, kmerId).
// If validMarkers is nullptr, it returns the count of markers without filling a buffer.
static size_t getMinimizerMarkersForRead(
    ReadId readId,
    const Reads& reads,
    int k,
    SimdSketcher* sketcher,
    const shared_ptr<KmerChecker>& kmerChecker,
    string& readSequence,
    std::vector<uint32_t>& positionBuffer,
    std::vector<std::pair<uint32_t, KmerId>>* validMarkers = nullptr
) {
    if(validMarkers) validMarkers->clear();
    const LongBaseSequenceView read = reads.getRead(readId);
    const uint64_t baseCount = read.baseCount;

    if(baseCount < uint64_t(k)) {
        return 0;
    }

    // Convert read to string for simd-minimizers.
    readSequence.resize(baseCount);
    for(uint64_t i = 0; i < baseCount; i++) {
        readSequence[i] = read[i].character();
    }

    // Compute canonical minimizer positions and merge with forced read ends (0 and baseCount - k).
    // Read end forcing ensures that the very first and last possible k-mers of every read
    // are included as markers, which is critical for graph connectivity at read boundaries.
    MinimizerList minimizerList = canonical_minimizer_positions(
        sketcher,
        readSequence.c_str(),
        readSequence.size());
    positionBuffer.assign(minimizerList.data, minimizerList.data + minimizerList.len);
    positionBuffer.push_back(0);
    positionBuffer.push_back(uint32_t(baseCount - k));
    free_minimizer_list(minimizerList);

    // Sort and deduplicate to handle cases where read ends were already selected as minimizers.
    std::sort(positionBuffer.begin(), positionBuffer.end());
    positionBuffer.erase(std::unique(positionBuffer.begin(), positionBuffer.end()), positionBuffer.end());

    const size_t uniqueCandidateCount = positionBuffer.size();

    // Pass 1 optimization: if kmerChecker is null (initial discovery) and we only need the count,
    // we can return uniqueCandidateCount directly.
    if(!validMarkers && !kmerChecker) {
        return uniqueCandidateCount;
    }

    if(validMarkers) validMarkers->reserve(uniqueCandidateCount);
    size_t validCount = 0;

    for(size_t i = 0; i < uniqueCandidateCount; i++) {
        const uint32_t position = positionBuffer[i];
        if(uint64_t(position) + uint64_t(k) > baseCount) continue;

        if(!kmerChecker) {
            if(validMarkers) {
                Kmer kmer;
                extractKmer(read, uint64_t(position), uint64_t(k), kmer);
                validMarkers->push_back({position, kmer.id(uint64_t(k))});
            }
            validCount++;
        } else {
            Kmer kmer;
            extractKmer(read, uint64_t(position), uint64_t(k), kmer);
            const KmerId kmerId = kmer.id(uint64_t(k));
            if(kmerChecker->isMarker(kmerId)) {
                if(validMarkers) validMarkers->push_back({position, kmerId});
                validCount++;
            }
        }
    }
    return validCount;
}

void Assembler::findMarkersSimdMinimizers(uint64_t threadCount, int k, int w)
{
    reads->checkReadsAreOpen();

    performanceLog << timestamp << "Finding markers using SIMD minimizers (k=" << k << ", w=" << w << ") in "
        << reads->readCount() << " reads." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // Store parameters.
    assemblerInfo->k = k;
    findMarkersSimdMinimizersData.k = k;
    findMarkersSimdMinimizersData.w = w;

    // Create the markers and markerKmerIds data structures.
    markers->createNew(largeDataName("Markers"), largeDataPageSize);
    markerKmerIds->createNew(largeDataName("MarkerKmerIds"), largeDataPageSize);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t batchSize = 100;

    // Pass 1: Count markers for each oriented read.
    markers->beginPass1(2 * readCount);
    markerKmerIds->beginPass1(2 * readCount);
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::findMarkersSimdMinimizersPass1, threadCount);

    // Pass 2: Store markers.
    markers->beginPass2();
    markerKmerIds->beginPass2();
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::findMarkersSimdMinimizersPass2, threadCount);

    markers->endPass2(false);
    markerKmerIds->endPass2(false);

    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    performanceLog << timestamp << "Finding markers using SIMD minimizers completed in " << tTotal << " s." << endl;
    cout << "Created " << markers->totalSize() << " markers using SIMD minimizers." << endl;
}

void Assembler::findMarkersSimdMinimizersPass1(size_t /* threadId */)
{
    const int k = findMarkersSimdMinimizersData.k;
    const int w = findMarkersSimdMinimizersData.w;
    SimdSketcher* sketcher = simd_sketcher_new(
        static_cast<uint8_t>(k), static_cast<uint8_t>(w));
    string readSequence;
    std::vector<uint32_t> positionBuffer;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            const size_t count = getMinimizerMarkersForRead(
                readId, *reads, k, sketcher, kmerChecker, readSequence, positionBuffer, nullptr);

            const uint64_t count64 = count;
            this->markers->incrementCount(OrientedReadId(readId, 0).getValue(), count64);
            this->markers->incrementCount(OrientedReadId(readId, 1).getValue(), count64);
            markerKmerIds->incrementCount(OrientedReadId(readId, 0).getValue(), count64);
            markerKmerIds->incrementCount(OrientedReadId(readId, 1).getValue(), count64);
        }
    }
    simd_sketcher_free(sketcher);
}

void Assembler::findMarkersSimdMinimizersPass2(size_t /* threadId */)
{
    const int k = findMarkersSimdMinimizersData.k;
    const int w = findMarkersSimdMinimizersData.w;
    SimdSketcher* sketcher = simd_sketcher_new(
        static_cast<uint8_t>(k), static_cast<uint8_t>(w));
    string readSequence;
    std::vector<uint32_t> positionBuffer;
    std::vector<std::pair<uint32_t, KmerId>> markerBuffer;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            const LongBaseSequenceView read = reads->getRead(readId);
            getMinimizerMarkersForRead(
                readId, *reads, k, sketcher, kmerChecker, readSequence, positionBuffer, &markerBuffer);

            if(markerBuffer.empty()) continue;

            CompressedMarker* markerPointerStrand0 = this->markers->begin(OrientedReadId(readId, 0).getValue());
            CompressedMarker* markerPointerStrand1 = this->markers->end(OrientedReadId(readId, 1).getValue()) - 1;

            KmerId* kmerIdPointerStrand0 = markerKmerIds->begin(OrientedReadId(readId, 0).getValue());
            KmerId* kmerIdPointerStrand1 = markerKmerIds->end(OrientedReadId(readId, 1).getValue()) - 1;

            for(const auto& [position, kmerId] : markerBuffer) {
                // Strand 0.
                markerPointerStrand0->position = position;
                ++markerPointerStrand0;
                *kmerIdPointerStrand0 = kmerId;
                ++kmerIdPointerStrand0;

                // Strand 1: reverse complement.
                Kmer kmer(kmerId, k);
                markerPointerStrand1->position = static_cast<uint32_t>(read.baseCount - k - position);
                --markerPointerStrand1;
                *kmerIdPointerStrand1 = kmer.reverseComplement(k).id(k);
                --kmerIdPointerStrand1;
            }
        }
    }
    simd_sketcher_free(sketcher);
}



void Assembler::accessMarkers()
{
    markers->accessExistingReadOnly(largeDataName("Markers"));
}

void Assembler::checkMarkersAreOpen() const
{
    if(!markers->isOpen()) {
        throw runtime_error("Markers are not accessible.");
    }
}


void Assembler::writeMarkers(ReadId readId, Strand strand, const string& fileName)
{
    // Check that we have what we need.
    DINARA_ASSERT(kmerChecker);
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    reads->checkReadId(readId);

    // Get the markers.
    const OrientedReadId orientedReadId(readId, strand);
    const auto orientedReadMarkers = (*markers)[orientedReadId.getValue()];

    // Write them out.
    ofstream csv(fileName);
    csv << "MarkerId,Ordinal,KmerId,Kmer,Position\n";
    for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const CompressedMarker& marker = orientedReadMarkers[ordinal];
        const MarkerId markerId = getMarkerId(orientedReadId, ordinal);
        const KmerId kmerId = getOrientedReadMarkerKmerId(orientedReadId, ordinal);
        const Kmer kmer(kmerId, assemblerInfo->k);
        csv << markerId << ",";
        csv << ordinal << ",";
        csv << kmerId << ",";
        csv << kmer << ",";
        csv << marker.position << "\n";
    }
}



// Get markers sorted by KmerId for a given OrientedReadId.
void Assembler::getMarkersSortedByKmerId(
    OrientedReadId orientedReadId,
    vector<MarkerWithOrdinal>& markersSortedByKmerId) const
{
    markersSortedByKmerId.resize(markers->size(orientedReadId.getValue()));
    getOrientedReadMarkers(orientedReadId, markersSortedByKmerId);
    sort(markersSortedByKmerId.begin(), markersSortedByKmerId.end());
}



// Given a marker by its OrientedReadId and ordinal,
// return the corresponding global marker id.
MarkerId Assembler::getMarkerId(
    OrientedReadId orientedReadId, uint32_t ordinal) const
{
    return
        (markers->begin(orientedReadId.getValue()) - markers->begin())
        + ordinal;
}

MarkerId Assembler::getReverseComplementMarkerId(
    OrientedReadId orientedReadId, uint32_t ordinal) const
{
    OrientedReadId orientedReadIdRc = orientedReadId;
    orientedReadIdRc.flipStrand();

    const uint32_t markerCount = uint32_t(markers->size(orientedReadId.getValue()));

    return getMarkerId(orientedReadIdRc, markerCount - 1 - ordinal);

}


// Inverse of the above: given a global marker id,
// return its OrientedReadId and ordinal.
// This requires a binary search in the markers toc.
// This could be avoided, at the cost of storing
// an additional 4 bytes per marker.
pair<OrientedReadId, uint32_t>
    Assembler::findMarkerId(MarkerId markerId) const
{
    return dinara::findMarkerId(markerId, *markers);
}



// Given a MarkerId, compute the MarkerId of the
// reverse complemented marker.
MarkerId Assembler::findReverseComplement(MarkerId markerId) const
{
	// Find the oriented read id and marker ordinal.
	OrientedReadId orientedReadId;
	uint32_t ordinal;
	tie(orientedReadId, ordinal) = findMarkerId(markerId);

	// Reverse complement.
	ordinal = uint32_t(markers->size(orientedReadId.getValue()) - 1 - ordinal);
	orientedReadId.flipStrand();

	// Return the corresponding Markerid.
	return getMarkerId(orientedReadId, ordinal);
}



void Assembler::computeMarkerKmerIds(uint64_t threadCount)
{
    performanceLog << timestamp << "Gathering marker KmerIds." << endl;

    // optimization: if we already have them (from findMarkersSimdClosedSyncmers), don't recompute.
    if(markerKmerIds->isOpen()) {
        performanceLog << timestamp << "Marker KmerIds are already present. Skipping computation." << endl;
        return;
    }

    // Check that we have what we need.
    checkMarkersAreOpen();
    const uint64_t readCount = reads->readCount();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Do it.
    // The layout is identical to that used by the markers.
    markerKmerIds->createNew(largeDataName("MarkerKmerIds"), largeDataPageSize);
    for(uint64_t readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(uint32_t(readId), 0);
        const OrientedReadId orientedReadId1(uint32_t(readId), 1);
        const uint64_t readMarkerCount = markers->size(orientedReadId0.getValue());
        DINARA_ASSERT(markers->size(orientedReadId1.getValue()) == readMarkerCount);
        for(uint64_t strand=0; strand<2; strand++) {
            markerKmerIds->appendVector(readMarkerCount);
        }
    }
    markerKmerIds->unreserve();
    const uint64_t batchSize = 100;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::computeMarkerKmerIdsThreadFunction, threadCount);



#if 0
    // Test the low level functions to extract Kmers/KmerIds.
    const uint64_t k = assemblerInfo->k;
    vector<Kmer> kmerVector;
    vector<KmerId> kmerIdVector;
    performanceLog << timestamp << "Testing." << endl;
    for(uint64_t readId=0; readId<readCount; readId++) {
        for(uint64_t strand=0; strand<2; strand++) {

            const OrientedReadId orientedReadId = OrientedReadId(ReadId(readId), Strand(strand));
            const auto orientedReadMarkers = (*markers)[orientedReadId.getValue()];
            const auto orientedReadMarkerKmerIds = (*markerKmerIds)[orientedReadId.getValue()];
            const uint64_t orientedReadMarkerCount = orientedReadMarkers.size();
            DINARA_ASSERT(orientedReadMarkerKmerIds.size() == orientedReadMarkerCount);

            kmerVector.resize(orientedReadMarkerCount);
            kmerIdVector.resize(orientedReadMarkerCount);
            const span<Kmer> kmerSpan(kmerVector);
            const span<KmerId> kmerIdSpan(kmerIdVector);

            getOrientedReadMarkerKmers(orientedReadId, kmerSpan);
            getOrientedReadMarkerKmerIds(orientedReadId, kmerIdSpan);

            for(uint64_t ordinal=0; ordinal<orientedReadMarkerCount; ordinal++) {
                DINARA_ASSERT(kmerVector[ordinal].id(k) == orientedReadMarkers[ordinal].kmerId);
                DINARA_ASSERT(kmerIdVector[ordinal] == orientedReadMarkers[ordinal].kmerId);

                DINARA_ASSERT(kmerVector[ordinal] == getOrientedReadMarkerKmer(orientedReadId, ordinal));
                DINARA_ASSERT(kmerIdVector[ordinal] == getOrientedReadMarkerKmerId(orientedReadId, ordinal));
            }
        }
    }
#endif

}



void Assembler::cleanupMarkerKmerIds()
{
    markerKmerIds->remove();
}



void Assembler::computeMarkerKmerIdsThreadFunction(uint64_t)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads in this batch.
        for(uint64_t readId=begin; readId!=end; ++readId) {

            const OrientedReadId orientedReadId0(uint32_t(readId), 0);
            const OrientedReadId orientedReadId1(uint32_t(readId), 1);

            getReadMarkerKmerIds(
                ReadId(readId),
                (*markerKmerIds)[orientedReadId0.getValue()],
                (*markerKmerIds)[orientedReadId1.getValue()]);
        }
    }

}



Kmer Assembler::getOrientedReadMarkerKmer(OrientedReadId orientedReadId, uint32_t ordinal) const
{
    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();

    if(strand == 0) {
        return getOrientedReadMarkerKmerStrand0(readId, ordinal);
    } else {
        return getOrientedReadMarkerKmerStrand1(readId, ordinal);
    }

}



Kmer Assembler::getOrientedReadMarkerKmerStrand0(ReadId readId, uint32_t ordinal0) const
{
    const uint64_t k = assemblerInfo->k;
    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];

    Kmer kmer0;
    extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);

    return kmer0;
}



Kmer Assembler::getOrientedReadMarkerKmerStrand1(ReadId readId, uint32_t ordinal1) const
{
    const uint64_t k = assemblerInfo->k;

    // We only have the read stored without reverse complement, so get it from there...
    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    const uint64_t ordinal0 = readMarkerCount - 1 - ordinal1;
    Kmer kmer0;
    extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);

    // ... then do the reverse complement.
    const Kmer kmer1 = kmer0.reverseComplement(k);
    return kmer1;
}



KmerId Assembler::getOrientedReadMarkerKmerId(OrientedReadId orientedReadId, uint32_t ordinal) const
{
    const Kmer kmer = getOrientedReadMarkerKmer(orientedReadId, ordinal);
    return KmerId(kmer.id(assemblerInfo->k));
}



// Get all marker Kmers for an oriented read.
void Assembler::getOrientedReadMarkerKmers(
    OrientedReadId orientedReadId,
    const span<Kmer>& kmers) const
{
    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();

    if(strand == 0) {
        getOrientedReadMarkerKmersStrand0(readId, kmers);
    } else {
        getOrientedReadMarkerKmersStrand1(readId, kmers);
    }
}



void Assembler::getOrientedReadMarkerKmersStrand0(ReadId readId, const span<Kmer>& kmers0) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    DINARA_ASSERT(kmers0.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        kmers0[ordinal0] = kmer0;
    }

}



void Assembler::getOrientedReadMarkerKmersStrand1(ReadId readId, const span<Kmer>& kmers1) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    DINARA_ASSERT(kmers1.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        const Kmer kmer1 = kmer0.reverseComplement(k);
        const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
        kmers1[ordinal1] = kmer1;
    }

}



// Get all marker KmerIds for an oriented read.
void Assembler::getOrientedReadMarkerKmerIds(
    OrientedReadId orientedReadId,
    const span<KmerId>& kmerIds) const
{
    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();

    if(strand == 0) {
        getOrientedReadMarkerKmerIdsStrand0(readId, kmerIds);
    } else {
        getOrientedReadMarkerKmerIdsStrand1(readId, kmerIds);
    }
}



void Assembler::getOrientedReadMarkerKmerIdsStrand0(ReadId readId, const span<KmerId>& kmerIds0) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    DINARA_ASSERT(kmerIds0.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        kmerIds0[ordinal0] = KmerId(kmer0.id(k));
    }

}



void Assembler::getOrientedReadMarkerKmerIdsStrand1(ReadId readId, const span<KmerId>& kmerIds1) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    DINARA_ASSERT(kmerIds1.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        const Kmer kmer1 = kmer0.reverseComplement(k);
        const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
        kmerIds1[ordinal1] = KmerId(kmer1.id(k));
    }

}



void Assembler::getOrientedReadMarkers(
    OrientedReadId orientedReadId,
    const span<MarkerWithOrdinal>& markers) const
{
    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();

    if(strand == 0) {
        getOrientedReadMarkersStrand0(readId, markers);
    } else {
        getOrientedReadMarkersStrand1(readId, markers);
    }

}



void Assembler::getOrientedReadMarkersStrand0(
    ReadId readId,
    const span<MarkerWithOrdinal>& markers0) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    DINARA_ASSERT(markers0.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        const CompressedMarker& compressedMarker0 = orientedReadMarkers0[ordinal0];
        const uint32_t position = compressedMarker0.position;
        Kmer kmer0;
        extractKmer(read, uint64_t(position), k, kmer0);
        markers0[ordinal0] = MarkerWithOrdinal(KmerId(kmer0.id(k)), position, uint32_t(ordinal0));
    }

}



void Assembler::getOrientedReadMarkersStrand1(
    ReadId readId,
    const span<MarkerWithOrdinal>& markers1) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const OrientedReadId orientedReadId1(readId, 1);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const auto orientedReadMarkers1 = (*markers)[orientedReadId1.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    DINARA_ASSERT(markers1.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
        const CompressedMarker& compressedMarker1 = orientedReadMarkers1[ordinal1];
        const uint32_t position1 = compressedMarker1.position;
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        const Kmer kmer1 = kmer0.reverseComplement(k);
        markers1[ordinal1] = MarkerWithOrdinal(KmerId(kmer1.id(k)), position1, uint32_t(ordinal1));
    }

}



void Assembler::getOrientedReadAlign6Markers(
    OrientedReadId orientedReadId,
    const span<Align6Marker>& align6Markers) const
{
    DINARA_ASSERT(kmerCounter and kmerCounter->isAvailable());

    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();

    if(strand == 0) {
        getOrientedReadAlign6MarkersStrand0(readId, align6Markers);
    } else {
        getOrientedReadAlign6MarkersStrand1(readId, align6Markers);
    }

    sort(align6Markers.begin(), align6Markers.end());
}



void Assembler::getOrientedReadAlign6MarkersStrand0(
    ReadId readId,
    const span<Align6Marker>& align6Markers) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    DINARA_ASSERT(align6Markers.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        const CompressedMarker& compressedMarker0 = orientedReadMarkers0[ordinal0];
        const uint32_t position = compressedMarker0.position;
        Kmer kmer0;
        extractKmer(read, uint64_t(position), k, kmer0);

        Align6Marker& align6Marker = align6Markers[ordinal0];
        align6Marker.kmerId = KmerId(kmer0.id(k));
        align6Marker.ordinal = uint32_t(ordinal0);
        align6Marker.setGlobalFrequency(kmerCounter->getFrequency(align6Marker.kmerId));
    }
}



void Assembler::getOrientedReadAlign6MarkersStrand1(
    ReadId readId,
    const span<Align6Marker>& align6Markers) const
{
    const uint64_t k = assemblerInfo->k;

    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(readId, 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    DINARA_ASSERT(align6Markers.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {
        const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        const Kmer kmer1 = kmer0.reverseComplement(k);

        Align6Marker& align6Marker = align6Markers[ordinal1];
        align6Marker.kmerId = KmerId(kmer1.id(k));
        align6Marker.ordinal = uint32_t(ordinal1);
        align6Marker.setGlobalFrequency(kmerCounter->getFrequency(align6Marker.kmerId));
    }
}



// Get all marker Kmers for a read in both orientations.
void Assembler::getReadMarkerKmers(
    ReadId readId,
    const span<Kmer>& kmers0,
    const span<Kmer>& kmers1) const
{
    const uint64_t k = assemblerInfo->k;

    // Access the information we need for this read.
    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(uint32_t(readId), 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    DINARA_ASSERT(kmers0.size() == readMarkerCount);
    DINARA_ASSERT(kmers1.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {

        // Strand 0.
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        kmers0[ordinal0] = kmer0;

        // Strand 1.
        const Kmer kmer1 = kmer0.reverseComplement(k);
        const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
        kmers1[ordinal1] = kmer1;
    }

}



// Get all marker KmerIds for a read in both orientations.
void Assembler::getReadMarkerKmerIds(
    ReadId readId,
    const span<KmerId>& kmerIds0,
    const span<KmerId>& kmerIds1) const
{
    // Get the marker length.
    const uint64_t k = assemblerInfo->k;

    // Access the information we need for this read.
    const auto read = reads->getRead(uint32_t(readId));
    const OrientedReadId orientedReadId0(uint32_t(readId), 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];
    const uint64_t readMarkerCount = orientedReadMarkers0.size();
    DINARA_ASSERT(kmerIds0.size() == readMarkerCount);
    DINARA_ASSERT(kmerIds1.size() == readMarkerCount);

    // Loop over all markers.
    for(uint64_t ordinal0=0; ordinal0<readMarkerCount; ordinal0++) {

        // Strand 0.
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        kmerIds0[ordinal0] = KmerId(kmer0.id(k));

        // Strand 1.
        const Kmer kmer1 = kmer0.reverseComplement(k);
        const uint64_t ordinal1 = readMarkerCount - 1 - ordinal0;
        kmerIds1[ordinal1] = KmerId(kmer1.id(k));
    }

}



// Get the Kmer for an oriented read at a given marker ordinal.
Kmer Assembler::getOrientedReadMarkerKmer(OrientedReadId orientedReadId, uint64_t ordinal) const
{
    const uint64_t k = assemblerInfo->k;

    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();
    const auto read = reads->getRead(readId);
    const OrientedReadId orientedReadId0(uint32_t(readId), 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];

    if(strand == 0) {

        const uint64_t ordinal0 = ordinal;
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        return kmer0;

    } else {

        const uint64_t ordinal0 = orientedReadMarkers0.size() - 1 - ordinal;
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        return kmer0.reverseComplement(k);

    }
}



// Get the KmerId for an oriented read at a given marker ordinal.
KmerId Assembler::getOrientedReadMarkerKmerId(OrientedReadId orientedReadId, uint64_t ordinal) const
{
    const uint64_t k = assemblerInfo->k;

    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();
    const auto read = reads->getRead(readId);
    const OrientedReadId orientedReadId0(uint32_t(readId), 0);
    const auto orientedReadMarkers0 = (*markers)[orientedReadId0.getValue()];

    if(strand == 0) {

        const uint64_t ordinal0 = ordinal;
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        return KmerId(kmer0.id(k));

    } else {

        const uint64_t ordinal0 = orientedReadMarkers0.size() - 1 - ordinal;
        Kmer kmer0;
        extractKmer(read, uint64_t(orientedReadMarkers0[ordinal0].position), k, kmer0);
        return KmerId(kmer0.reverseComplement(k).id(k));

    }
}



void Assembler::countKmers(
    uint64_t threadCount,
    const string& globalFrequencyOverrideDirectory)
{
    DINARA_ASSERT(markers->isOpen());
    kmerCounter = make_shared<KmerCounter>(
        assemblerInfo->k, getReads(), *markers, *this, threadCount);



    if(globalFrequencyOverrideDirectory.empty()) {
        // No override, just create the histogram.
        kmerCounter->createHistogram();
    } else {

        // The globalFrequencyOverrideDirectory must specify an absolute path.
        if(globalFrequencyOverrideDirectory[0] != '/') {
            throw runtime_error("Option --Kmers.globalFrequencyOverrideDirectory must specify an absolute path. "
                "A relative path is not accepted.");
        }

        // Override the frequencies.
        MappedMemoryOwner mappedMemoryOwner;
        mappedMemoryOwner.largeDataFileNamePrefix = globalFrequencyOverrideDirectory + "/";
        mappedMemoryOwner.largeDataPageSize = 4096;

        KmerCounter otherKmerCounter(assemblerInfo->k, mappedMemoryOwner);

        cout << timestamp << "Overriding k-mer global frequencies." << endl;
        kmerCounter->overrideFrequencies(otherKmerCounter);
        cout << timestamp << "Done overriding k-mer global frequencies." << endl;

        // Copy the histogram.
        kmerCounter->histogram.createNew(largeDataName("KmerCounterHistogram"), largeDataPageSize);
        kmerCounter->histogram.resize(otherKmerCounter.histogram.size());
        copy(otherKmerCounter.histogram.begin(), otherKmerCounter.histogram.end(), kmerCounter->histogram.begin());
    }


    ofstream csv("KmerFrequencyHistogram.csv");
    kmerCounter->writeHistogram(csv);
    kmerCounter->getHistogramInfo(assemblerInfo->kmerDistributionInfo);

    cout << "Marker k-mer coverage distribution:"
        " low "   << assemblerInfo->kmerDistributionInfo.coverageLow <<
        ", peak " << assemblerInfo->kmerDistributionInfo.coveragePeak <<
        ", high " << assemblerInfo->kmerDistributionInfo.coverageHigh << endl;
}



void Assembler::accessKmerCounts()
{
    DINARA_ASSERT(markers->isOpen());
    kmerCounter = make_shared<KmerCounter>(assemblerInfo->k, *this);
}



void Assembler::createMarkerKmers(uint64_t threadCount)
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        getReads(),
        *markers,
        threadCount);
}



void Assembler::accessMarkerKmers()
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        getReads(),
        *markers);
}


// Count k-mers from pre-calculated Marker KmerIds.
void Assembler::countKmersFromMarkerKmerIds(uint64_t threadCount)
{
    DINARA_ASSERT(markerKmerIds->isOpen());
    
    // Create KmerCounter from markerKmerIds.
    kmerCounter = make_shared<KmerCounter>(
        assemblerInfo->k,
        *markerKmerIds, 
        *this, 
        threadCount);
    
    kmerCounter->createHistogram(); 

    ofstream csv("KmerFrequencyHistogram.csv");
    kmerCounter->writeHistogram(csv);
    kmerCounter->getHistogramInfo(assemblerInfo->kmerDistributionInfo);
    
    // Build the frequency LUT for O(1) lookups in applyKmerCountFilter.
    kmerCounter->buildFrequencyLUT();

    cout << "Marker k-mer coverage distribution:"
        " low "   << assemblerInfo->kmerDistributionInfo.coverageLow <<
        ", peak " << assemblerInfo->kmerDistributionInfo.coveragePeak <<
        ", high " << assemblerInfo->kmerDistributionInfo.coverageHigh << endl;
}

/**
 * @brief High-performance marker frequency filtering using precomputed caches.
 *
 * This function filters markers based on k-mer frequency, keeping only those
 * within [minFreq, maxFreq]. It uses a two-pass algorithm with extensive caching
 * to minimize redundant computation:
 *
 * **Algorithm Overview:**
 *   1. **Parallel Initialization**: Allocate per-read caches for validity bits,
 *      read lengths, and reverse complement KmerIds.
 *   2. **Pass 1 (Count & Cache)**: For each marker, compute canonical KmerId,
 *      lookup frequency in O(1) LUT, and if valid:
 *      - Set bit in packed validity bitset
 *      - Cache the rcKmerId for Pass 2
 *      - Increment marker count for VectorOfVectors allocation
 *   3. **Pass 2 (Store)**: Using only cached data (no frequency lookups):
 *      - Check validity from bitset
 *      - Copy marker to new structure
 *      - Use cached rcKmerId for strand 1 (no reverseComplement needed)
 *
 * **Performance Optimizations:**
 *   - O(1) frequency lookup via pre-built unordered_map LUT
 *   - Packed validity bitset (8 markers per byte)
 *   - Cached read lengths (avoids Reads data structure access)
 *   - Cached rcKmerIds (eliminates reverseComplement in Pass 2)
 *   - Parallel initialization of per-read caches
 *   - Bitwise operations for bitset indexing (>> 3, & 7)
 *
 * @param minFreq Minimum k-mer frequency threshold (inclusive)
 * @param maxFreq Maximum k-mer frequency threshold (inclusive)
 * @param threadCount Number of worker threads (0 = auto-detect)
 */
void Assembler::applyKmerCountFilter(uint64_t minFreq, uint64_t maxFreq, uint64_t threadCount, bool filterPalindromes)
{
    performanceLog << timestamp << "Filtering markers by frequency [" 
                   << minFreq << ", " << maxFreq << "] and palindromes: " 
                   << (filterPalindromes ? "yes" : "no") << "." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // =========================================================================
    // Phase 0: Prerequisites and Data Structure Setup
    // =========================================================================
    checkMarkersAreOpen();
    DINARA_ASSERT(markerKmerIds->isOpen());
    if(!kmerCounter) {
        throw runtime_error("KmerCounter is required for marker filtering.");
    }

    // Swap current markers to "old" pointers for in-place filtering.
    const string markersName = markers->getName(); 
    const string kmerIdsName = markerKmerIds->getName();

    applyKmerCountFilterData.oldMarkers = markers;
    if(!markersName.empty()) {
        applyKmerCountFilterData.oldMarkers->rename(markersName + ".old");
    }
    markers = make_shared<MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>>();
    
    applyKmerCountFilterData.oldMarkerKmerIds = markerKmerIds;
    if(!kmerIdsName.empty()) {
        applyKmerCountFilterData.oldMarkerKmerIds->rename(kmerIdsName + ".old");
    }
    markerKmerIds = make_shared<MemoryMapped::VectorOfVectors<KmerId, uint64_t>>();

    applyKmerCountFilterData.minFreq = minFreq;
    applyKmerCountFilterData.maxFreq = maxFreq;
    applyKmerCountFilterData.filterPalindromes = filterPalindromes;

    markers->createNew(markersName, largeDataPageSize);
    markerKmerIds->createNew(kmerIdsName, largeDataPageSize);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t batchSize = 100;

    // =========================================================================
    // Phase 1: Parallel Cache Initialization
    // =========================================================================
    // Allocate per-read data structures in parallel to reduce init time.
    // Each read gets:
    //   - validityBits: packed bitset (1 bit per marker)
    //   - readLength: cached to avoid Reads access in Pass 2
    //   - rcKmerIds: vector of rcKmerIds for valid markers only
    applyKmerCountFilterData.markerValidity.resize(readCount);
    applyKmerCountFilterData.readLengths.resize(readCount);
    applyKmerCountFilterData.rcKmerIds.resize(readCount);
    
    {
        std::vector<std::thread> initThreads;
        initThreads.reserve(threadCount);
        
        auto initFunction = [this, readCount, threadCount](size_t tid) {
            const ReadId startRead = static_cast<ReadId>((readCount * tid) / threadCount);
            const ReadId endRead = static_cast<ReadId>((readCount * (tid + 1)) / threadCount);
            
            for(ReadId rId = startRead; rId < endRead; ++rId) {
                const size_t n = (*applyKmerCountFilterData.oldMarkers)[OrientedReadId(rId, 0).getValue()].size();
                applyKmerCountFilterData.markerValidity[rId].resize((n + 7) >> 3, 0);  // ceil(n/8)
                applyKmerCountFilterData.readLengths[rId] = reads->getReadRawSequenceLength(rId);
                applyKmerCountFilterData.rcKmerIds[rId].reserve(n >> 2);  // ~25% expected valid
            }
        };
        
        for(size_t tid = 0; tid < threadCount; ++tid) {
            initThreads.emplace_back(initFunction, tid);
        }
        for(auto& t : initThreads) {
            t.join();
        }
    }

    // =========================================================================
    // Phase 2: Pass 1 - Count Valid Markers & Populate Caches
    // =========================================================================
    markers->beginPass1(2 * readCount);
    markerKmerIds->beginPass1(2 * readCount);
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::applyKmerCountFilterThreadFunctionPass1, threadCount);

    // =========================================================================
    // Phase 3: Pass 2 - Store Valid Markers Using Cached Data
    // =========================================================================
    markers->beginPass2();
    markerKmerIds->beginPass2();
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::applyKmerCountFilterThreadFunctionPass2, threadCount);
    
    markers->endPass2(false);
    markerKmerIds->endPass2(false);

    // =========================================================================
    // Phase 4: Cleanup
    // =========================================================================
    const uint64_t oldTotalSize = applyKmerCountFilterData.oldMarkers->totalSize();

    applyKmerCountFilterData.oldMarkers->remove();
    applyKmerCountFilterData.oldMarkerKmerIds->remove();
    applyKmerCountFilterData.oldMarkers.reset();
    applyKmerCountFilterData.oldMarkerKmerIds.reset();
    
    // Release cache memory.
    applyKmerCountFilterData.markerValidity.clear();
    applyKmerCountFilterData.markerValidity.shrink_to_fit();
    applyKmerCountFilterData.readLengths.clear();
    applyKmerCountFilterData.readLengths.shrink_to_fit();
    applyKmerCountFilterData.rcKmerIds.clear();
    applyKmerCountFilterData.rcKmerIds.shrink_to_fit();

    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    performanceLog << timestamp << "Marker filtering completed in " << tTotal << " s." << endl;
    cout << "Filtered markers: " << markers->totalSize() << " / " << oldTotalSize 
         << " (" << (100.0 * markers->totalSize() / oldTotalSize) << "%)." << endl;
}

/**
 * @brief Pass 1: Count valid markers and populate caches.
 *
 * For each marker in each read:
 *   1. Compute canonical KmerId (min of forward and RC)
 *   2. Lookup frequency in O(1) LUT
 *   3. If within [minFreq, maxFreq]:
 *      - Set validity bit in packed bitset
 *      - Cache rcKmerId for use in Pass 2
 *      - Increment marker count
 *
 * This is the computationally expensive pass, but it only runs once.
 * All computed data is cached for Pass 2 to consume.
 */
void Assembler::applyKmerCountFilterThreadFunctionPass1(size_t /* threadId */)
{
    const uint64_t k = assemblerInfo->k;
    const uint64_t minF = applyKmerCountFilterData.minFreq;
    const uint64_t maxF = applyKmerCountFilterData.maxFreq;
    
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            
            const OrientedReadId orientedReadId(readId, 0);
            const auto oldMarkers = (*applyKmerCountFilterData.oldMarkers)[orientedReadId.getValue()];
            const auto oldKmerIds = (*applyKmerCountFilterData.oldMarkerKmerIds)[orientedReadId.getValue()];
            const size_t n = oldMarkers.size();
            
            if(n == 0) continue;
            
            auto& validityBits = applyKmerCountFilterData.markerValidity[readId];
            auto& rcCache = applyKmerCountFilterData.rcKmerIds[readId];
            
            uint64_t validCount = 0;
            
            for(size_t i = 0; i < n; ++i) {
                // Compute canonical form: min(kmerId, rcKmerId).
                const KmerId kmerId = oldKmerIds[i];
                const Kmer kmer(kmerId, k);
                const KmerId rcKmerId = kmer.reverseComplement(k).id(k);
                const KmerId canonical = std::min(kmerId, rcKmerId);

                // O(1) frequency lookup.
                const uint64_t freq = kmerCounter->getFrequencyFast(canonical);
                
                if(freq >= minF && freq <= maxF) {
                    
                    // Palindrome check: if requested, skip markers where kmerId == rcKmerId.
                    // Palindromic k-mers can cause ambiguity in directed graph construction 
                    // because their forward and reverse orientations are identical.
                    if(applyKmerCountFilterData.filterPalindromes && kmerId == rcKmerId) {
                        continue;
                    }

                    // Set validity bit using bitwise ops (>> 3 = /8, & 7 = %8).
                    validityBits[i >> 3] |= (uint8_t(1) << (i & 7));
                    // Cache rcKmerId for Pass 2 (avoids recomputing reverseComplement).
                    rcCache.push_back(rcKmerId);
                    ++validCount;
                }
            }

            // Update marker counts for both strands.
            const uint64_t orid0 = orientedReadId.getValue();
            const uint64_t orid1 = OrientedReadId(readId, 1).getValue();
            markers->incrementCount(orid0, validCount);
            markers->incrementCount(orid1, validCount);
            markerKmerIds->incrementCount(orid0, validCount);
            markerKmerIds->incrementCount(orid1, validCount);
        }
    }
}


/**
 * @brief Pass 2: Store valid markers using only cached data.
 *
 * This pass is highly optimized:
 *   - Validity is checked from packed bitset (no frequency lookup)
 *   - Read length is from cache (no Reads access)
 *   - rcKmerId is from cache (no reverseComplement computation)
 *
 * The only work per valid marker is:
 *   - Copy position and kmerId to strand 0
 *   - Compute strand 1 position and copy cached rcKmerId
 */
void Assembler::applyKmerCountFilterThreadFunctionPass2(size_t /* threadId */)
{
    const uint64_t k = assemblerInfo->k;
    
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            
            const OrientedReadId orientedReadId0(readId, 0);
            const auto oldMarkers = (*applyKmerCountFilterData.oldMarkers)[orientedReadId0.getValue()];
            const auto oldKmerIds = (*applyKmerCountFilterData.oldMarkerKmerIds)[orientedReadId0.getValue()];
            const size_t n = oldMarkers.size();
            
            if(n == 0) continue;

            // All data from caches - no external lookups.
            const uint64_t readLen = applyKmerCountFilterData.readLengths[readId];
            const auto& validityBits = applyKmerCountFilterData.markerValidity[readId];
            const auto& rcCache = applyKmerCountFilterData.rcKmerIds[readId];

            // Output pointers.
            CompressedMarker* outMarker0 = markers->begin(orientedReadId0.getValue());
            CompressedMarker* outMarker1 = markers->end(OrientedReadId(readId, 1).getValue()) - 1;
            KmerId* outKmerId0 = markerKmerIds->begin(orientedReadId0.getValue());
            KmerId* outKmerId1 = markerKmerIds->end(OrientedReadId(readId, 1).getValue()) - 1;

            size_t rcIdx = 0;  // Index into dense rcCache array.

            for(size_t i = 0; i < n; ++i) {
                // Check validity bit using bitwise ops.
                if((validityBits[i >> 3] >> (i & 7)) & 1) {
                    const uint32_t pos = oldMarkers[i].position;
                    const KmerId kmerId = oldKmerIds[i];
                    
                    // Strand 0: direct copy.
                    outMarker0->position = pos;
                    *outKmerId0 = kmerId;
                    ++outMarker0;
                    ++outKmerId0;

                    // Strand 1: reverse complement position, cached rcKmerId.
                    outMarker1->position = static_cast<uint32_t>(readLen - k - pos);
                    *outKmerId1 = rcCache[rcIdx++];
                    --outMarker1;
                    --outKmerId1;
                }
            }
        }
    }
}
