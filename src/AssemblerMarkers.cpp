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

    markers.createNew(largeDataName("Markers"), largeDataPageSize);
    MarkerFinder markerFinder(
        assemblerInfo->k,
        *kmerChecker,
        getReads(),
        markers,
        threadCount);

}



// Find markers using SIMD-accelerated minimizers.
// This uses the simd-minimizers-c library to compute canonical minimizers
// for each read, and stores the minimizer positions as markers.
// Helper to get sorted unique positions and KmerIds.
static std::vector<std::pair<uint32_t, KmerId>> getSortedUniquePositionsAndIds(
    ReadId readId,
    const Reads& reads,
    int k,
    SimdSketcher* sketcher,
    const shared_ptr<KmerChecker>& kmerChecker,
    string& readSequence // reusable buffer
) {
    const LongBaseSequenceView read = reads.getRead(readId);
    
    if(read.baseCount < uint64_t(k)) {
        return {};
    }

    // Convert read to string for simd-minimizers.
    readSequence.clear();
    for(uint64_t i = 0; i < read.baseCount; i++) {
        readSequence.push_back(read[i].character());
    }

    // Compute minimizer positions using simd-minimizers.
    MinimizerList minimizerPositions = canonical_minimizer_positions(
        sketcher,
        readSequence.c_str(),
        readSequence.size());

    // Copy to vector, sort, and remove duplicates based on position.
    std::vector<uint32_t> positions(minimizerPositions.data, 
                                     minimizerPositions.data + minimizerPositions.len);
    
    // Free the minimizer list immediately.
    free_minimizer_list(minimizerPositions);
    
    // Sort and uniq.
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
    
    std::vector<std::pair<uint32_t, KmerId>> validMarkers;
    validMarkers.reserve(positions.size());

    for(const uint32_t position : positions) {
        Kmer kmer;
        extractKmer(read, uint64_t(position), uint64_t(k), kmer);
        const KmerId kmerId = kmer.id(uint64_t(k));
        
        if(!kmerChecker || kmerChecker->isMarker(kmerId)) {
            validMarkers.push_back({position, kmerId});
        }
    }

    return validMarkers;
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

    // Create the markers data structure.
    markers.createNew(largeDataName("Markers"), largeDataPageSize);
    
    // Create the markerKmerIds data structure.
    markerKmerIds.createNew(largeDataName("MarkerKmerIds"), largeDataPageSize);

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Pass 1: Count markers for each oriented read.
    const uint64_t readCount = reads->readCount();
    const uint64_t batchSize = 100; // Adjust batch size as needed.

    markers.beginPass1(2 * readCount);
    markerKmerIds.beginPass1(2 * readCount);
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::findMarkersSimdMinimizersPass1, threadCount);

    // Pass 2: Store markers.
    markers.beginPass2();
    markerKmerIds.beginPass2();
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::findMarkersSimdMinimizersPass2, threadCount);
    
    markers.endPass2(false);
    markerKmerIds.endPass2(false);

    // Report.
    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    performanceLog << timestamp << "Finding markers using SIMD minimizers completed in " << tTotal << " s." << endl;
    cout << "Created " << markers.totalSize() << " markers using SIMD minimizers." << endl;
}

void Assembler::findMarkersSimdMinimizersPass1(size_t /* threadId */)
{
    // Initialize thread-local sketcher and buffer.
    SimdSketcher* sketcher = simd_sketcher_new(
        static_cast<uint8_t>(findMarkersSimdMinimizersData.k), 
        static_cast<uint8_t>(findMarkersSimdMinimizersData.w));
    string readSequence;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            const auto markers = getSortedUniquePositionsAndIds(
                readId, *reads, findMarkersSimdMinimizersData.k, sketcher, kmerChecker, readSequence);

            this->markers.incrementCount(OrientedReadId(readId, 0).getValue(), markers.size());
            this->markers.incrementCount(OrientedReadId(readId, 1).getValue(), markers.size());
            
            // MarkerKmerIds must match markers counts exactly.
            markerKmerIds.incrementCount(OrientedReadId(readId, 0).getValue(), markers.size());
            markerKmerIds.incrementCount(OrientedReadId(readId, 1).getValue(), markers.size());
        }
    }
    simd_sketcher_free(sketcher);
}

void Assembler::findMarkersSimdMinimizersPass2(size_t /* threadId */)
{
    // Initialize thread-local sketcher and buffer.
    SimdSketcher* sketcher = simd_sketcher_new(
        static_cast<uint8_t>(findMarkersSimdMinimizersData.k), 
        static_cast<uint8_t>(findMarkersSimdMinimizersData.w));
    string readSequence;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            const LongBaseSequenceView read = reads->getRead(readId); // Need read for baseCount
            const auto markers = getSortedUniquePositionsAndIds(
                readId, *reads, findMarkersSimdMinimizersData.k, sketcher, kmerChecker, readSequence);

            if(markers.empty()) continue;

            CompressedMarker* markerPointerStrand0 = this->markers.begin(OrientedReadId(readId, 0).getValue());
            CompressedMarker* markerPointerStrand1 = this->markers.end(OrientedReadId(readId, 1).getValue()) - 1;
            
            KmerId* kmerIdPointerStrand0 = markerKmerIds.begin(OrientedReadId(readId, 0).getValue());
            KmerId* kmerIdPointerStrand1 = markerKmerIds.end(OrientedReadId(readId, 1).getValue()) - 1;

            for(const auto& val : markers) {
                uint32_t position = val.first;
                KmerId kmerId = val.second;
                
                // Strand 0
                markerPointerStrand0->position = position;
                ++markerPointerStrand0;
                
                *kmerIdPointerStrand0 = kmerId;
                ++kmerIdPointerStrand0;

                // Strand 1: reverse complement position.
                // KmerId should be the RC of kmerId?
                // Wait, getSortedUniquePositionsAndIds returns Canonical KmerId if check was canonical?
                // MarkerKmers usually stores RC KmerId for Strand 1?
                // Or does it store the canonical ID?
                // Assembler::computeMarkerKmerIds stores the ID of the kmer *on that strand*.
                // So for Strand 1, we need RC of kmerId.
                
                // Let's compute RC KmerId.
                Kmer kmer(kmerId, assemblerInfo->k);
                Kmer rcKmer = kmer.reverseComplement(assemblerInfo->k);
                KmerId rcKmerId = rcKmer.id(assemblerInfo->k);

                markerPointerStrand1->position = static_cast<uint32_t>(read.baseCount - findMarkersSimdMinimizersData.k - position);
                --markerPointerStrand1;
                
                *kmerIdPointerStrand1 = rcKmerId;
                --kmerIdPointerStrand1;
            }
        }
    }
    simd_sketcher_free(sketcher);
}



void Assembler::accessMarkers()
{
    markers.accessExistingReadOnly(largeDataName("Markers"));
}

void Assembler::checkMarkersAreOpen() const
{
    if(!markers.isOpen()) {
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
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];

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
    markersSortedByKmerId.resize(markers.size(orientedReadId.getValue()));
    getOrientedReadMarkers(orientedReadId, markersSortedByKmerId);
    sort(markersSortedByKmerId.begin(), markersSortedByKmerId.end());
}



// Given a marker by its OrientedReadId and ordinal,
// return the corresponding global marker id.
MarkerId Assembler::getMarkerId(
    OrientedReadId orientedReadId, uint32_t ordinal) const
{
    return
        (markers.begin(orientedReadId.getValue()) - markers.begin())
        + ordinal;
}

MarkerId Assembler::getReverseComplementMarkerId(
    OrientedReadId orientedReadId, uint32_t ordinal) const
{
    OrientedReadId orientedReadIdRc = orientedReadId;
    orientedReadIdRc.flipStrand();

    const uint32_t markerCount = uint32_t(markers.size(orientedReadId.getValue()));

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
    return dinara::findMarkerId(markerId, markers);
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
	ordinal = uint32_t(markers.size(orientedReadId.getValue()) - 1 - ordinal);
	orientedReadId.flipStrand();

	// Return the corresponding Markerid.
	return getMarkerId(orientedReadId, ordinal);
}



void Assembler::computeMarkerKmerIds(uint64_t threadCount)
{
    performanceLog << timestamp << "Gathering marker KmerIds." << endl;

    // optimization: if we already have them (from findMarkersSimdMinimizers), don't recompute.
    if(markerKmerIds.isOpen()) {
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
    markerKmerIds.createNew(largeDataName("MarkerKmerIds"), largeDataPageSize);
    for(uint64_t readId=0; readId<readCount; readId++) {
        const OrientedReadId orientedReadId0(uint32_t(readId), 0);
        const OrientedReadId orientedReadId1(uint32_t(readId), 1);
        const uint64_t readMarkerCount = markers.size(orientedReadId0.getValue());
        DINARA_ASSERT(markers.size(orientedReadId1.getValue()) == readMarkerCount);
        for(uint64_t strand=0; strand<2; strand++) {
            markerKmerIds.appendVector(readMarkerCount);
        }
    }
    markerKmerIds.unreserve();
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
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const auto orientedReadMarkerKmerIds = markerKmerIds[orientedReadId.getValue()];
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
    markerKmerIds.remove();
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
                markerKmerIds[orientedReadId0.getValue()],
                markerKmerIds[orientedReadId1.getValue()]);
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];

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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
    const auto orientedReadMarkers1 = markers[orientedReadId1.getValue()];
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
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
    const OrientedReadId orientedReadId1(readId, 1);
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
    const auto orientedReadMarkers1 = markers[orientedReadId1.getValue()];
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];
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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];

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
    const auto orientedReadMarkers0 = markers[orientedReadId0.getValue()];

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
    DINARA_ASSERT(markers.isOpen());
    kmerCounter = make_shared<KmerCounter>(
        assemblerInfo->k, getReads(), markers, *this, threadCount);



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
    DINARA_ASSERT(markers.isOpen());
    kmerCounter = make_shared<KmerCounter>(assemblerInfo->k, *this);
}



void Assembler::createMarkerKmers(uint64_t threadCount)
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        getReads(),
        markers,
        threadCount);
}



void Assembler::accessMarkerKmers()
{
    const MappedMemoryOwner& mappedMemoryOwner = *this;

    markerKmers = make_shared<MarkerKmers>(
        assemblerInfo->k,
        mappedMemoryOwner,
        getReads(),
        markers);
}


// Count k-mers from pre-calculated Marker KmerIds.
void Assembler::countKmersFromMarkerKmerIds(uint64_t threadCount)
{
    DINARA_ASSERT(markerKmerIds.isOpen());
    
    // Create KmerCounter from markerKmerIds.
    kmerCounter = make_shared<KmerCounter>(
        assemblerInfo->k,
        markerKmerIds, 
        *this, 
        threadCount);
    
    kmerCounter->createHistogram(); 

    ofstream csv("KmerFrequencyHistogram.csv");
    kmerCounter->writeHistogram(csv);
    kmerCounter->getHistogramInfo(assemblerInfo->kmerDistributionInfo);

    cout << "Marker k-mer coverage distribution:"
        " low "   << assemblerInfo->kmerDistributionInfo.coverageLow <<
        ", peak " << assemblerInfo->kmerDistributionInfo.coveragePeak <<
        ", high " << assemblerInfo->kmerDistributionInfo.coverageHigh << endl;
}

// Prune existing markers based on KmerCounter frequencies using markerKmerIds.
void Assembler::applyKmerCountFilter(uint64_t minFreq, uint64_t maxFreq, uint64_t threadCount)
{
    performanceLog << timestamp << "Filtering markers based on KmerCounter frequency (Fast ID path)." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // Check prerequisites.
    checkMarkersAreOpen();
    DINARA_ASSERT(markerKmerIds.isOpen());
    if(!kmerCounter) {
        throw runtime_error("KmerCounter is required for marker filtering.");
    }

    // Move current markers/ids to oldMarkers by renaming files on disk.
    // VectorOfVectors does not support move assignment, so we use renaming.
    const string markersName = markers.getName(); 
    const string markersNameOld = markersName + "-Old";
    markers.rename(markersNameOld);
    applyKmerCountFilterData.oldMarkers.accessExistingReadOnly(markersNameOld);
    markers.close(); // Close the current handle so we can create a new one.

    const string kmerIdsName = markerKmerIds.getName();
    const string kmerIdsNameOld = kmerIdsName + "-Old";
    markerKmerIds.rename(kmerIdsNameOld);
    applyKmerCountFilterData.oldMarkerKmerIds.accessExistingReadOnly(kmerIdsNameOld);
    markerKmerIds.close(); // Close the current handle so we can create a new one.

    applyKmerCountFilterData.minFreq = minFreq;
    applyKmerCountFilterData.maxFreq = maxFreq;

    // Create new markers structure (overwriting/creating fresh files).
    markers.createNew(markersName, largeDataPageSize);
    markerKmerIds.createNew(kmerIdsName, largeDataPageSize);

    // Adjust threads.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t batchSize = 100;

    // Pass 1: Count valid markers.
    markers.beginPass1(2 * readCount);
    markerKmerIds.beginPass1(2 * readCount);
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::applyKmerCountFilterThreadFunctionPass1, threadCount);

    // Pass 2: Store valid markers.
    markers.beginPass2();
    markerKmerIds.beginPass2();
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::applyKmerCountFilterThreadFunctionPass2, threadCount);
    
    markers.endPass2(false);
    markerKmerIds.endPass2(false);

    // Capture old size before removal.
    const uint64_t oldTotalSize = applyKmerCountFilterData.oldMarkers.totalSize();

    // Clean up old markers.
    applyKmerCountFilterData.oldMarkers.remove();
    applyKmerCountFilterData.oldMarkerKmerIds.remove();

    // Report.
    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    performanceLog << timestamp << "Marker filtering completed in " << tTotal << " s." << endl;
    cout << "Filtered markers: kept " << markers.totalSize() << " out of " 
         << oldTotalSize << "." << endl;
}

void Assembler::applyKmerCountFilterThreadFunctionPass1(size_t /* threadId */)
{
    const uint64_t k = assemblerInfo->k;
    
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            
            // Get markers for this read from oldMarkers.
            const OrientedReadId orientedReadId(readId, 0);
            const auto oldReadMarkers = applyKmerCountFilterData.oldMarkers[orientedReadId.getValue()];
            const auto oldReadMarkerKmerIds = applyKmerCountFilterData.oldMarkerKmerIds[orientedReadId.getValue()];
            
            uint64_t validCount = 0;
            if(oldReadMarkers.size() > 0) {
                for(size_t i=0; i<oldReadMarkers.size(); i++) {
                    // Check frequency using ID (no Read access).
                    KmerId kmerId = oldReadMarkerKmerIds[i]; // Strand 0 ID.
                    
                    // Canonicalize for checking (since KmerCounter tracks canonical).
                    // Although KmerCounter could be built non-canonical, it's safer to query canonical.
                    Kmer kmer(kmerId, k);
                    KmerId rcKmerId = kmer.reverseComplement(k).id(k);
                    KmerId canonical = std::min(kmerId, rcKmerId);

                    const uint64_t freq = kmerCounter->getFrequency(canonical);
                    
                    if(freq >= applyKmerCountFilterData.minFreq && freq <= applyKmerCountFilterData.maxFreq) {
                        validCount++;
                    }
                }
            }

            markers.incrementCount(OrientedReadId(readId, 0).getValue(), validCount);
            markers.incrementCount(OrientedReadId(readId, 1).getValue(), validCount);
            
            markerKmerIds.incrementCount(OrientedReadId(readId, 0).getValue(), validCount);
            markerKmerIds.incrementCount(OrientedReadId(readId, 1).getValue(), validCount);
        }
    }
}

void Assembler::applyKmerCountFilterThreadFunctionPass2(size_t /* threadId */)
{
    const uint64_t k = assemblerInfo->k;
    
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            
            const OrientedReadId orientedReadId0(readId, 0);
            const auto oldReadMarkers = applyKmerCountFilterData.oldMarkers[orientedReadId0.getValue()];
            const auto oldReadMarkerKmerIds = applyKmerCountFilterData.oldMarkerKmerIds[orientedReadId0.getValue()];
            
            if(oldReadMarkers.size() == 0) continue;

            const LongBaseSequenceView read = reads->getRead(readId);

            // Get pointers for storing markers.
            CompressedMarker* markerPointerStrand0 = markers.begin(orientedReadId0.getValue());
            CompressedMarker* markerPointerStrand1 = markers.end(OrientedReadId(readId, 1).getValue()) - 1;
            
            KmerId* kmerIdPointerStrand0 = markerKmerIds.begin(orientedReadId0.getValue());
            KmerId* kmerIdPointerStrand1 = markerKmerIds.end(OrientedReadId(readId, 1).getValue()) - 1;

            for(size_t i=0; i<oldReadMarkers.size(); i++) {
                KmerId kmerId = oldReadMarkerKmerIds[i]; // Strand 0 ID
                
                Kmer kmer(kmerId, k);
                KmerId rcKmerId = kmer.reverseComplement(k).id(k);
                KmerId canonical = std::min(kmerId, rcKmerId);

                const uint64_t freq = kmerCounter->getFrequency(canonical);
                
                if(freq >= applyKmerCountFilterData.minFreq && freq <= applyKmerCountFilterData.maxFreq) {
                    const uint32_t position = oldReadMarkers[i].position;
                    
                    // Strand 0
                    markerPointerStrand0->position = position;
                    ++markerPointerStrand0;
                    
                    *kmerIdPointerStrand0 = kmerId;
                    ++kmerIdPointerStrand0;

                    // Strand 1 (Derived from Strand 0 info + read len)
                    markerPointerStrand1->position = static_cast<uint32_t>(read.baseCount - k - position);
                    --markerPointerStrand1;
                    
                    *kmerIdPointerStrand1 = rcKmerId; // Store correct Strand 1 ID (RC)
                    --kmerIdPointerStrand1;
                }
            }
        }
    }
}
