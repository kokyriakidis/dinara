// PngImage.hpp must be included first because of png issues on Ubuntu 16.04.
#include "PngImage.hpp"

// Dinara.
#include "Assembler.hpp"
#include "hifiasmCoordinateTransforms.hpp"
#include "Alignment.hpp"
#include "AlignmentGraph.hpp"
#include "Align4.hpp"
#include "Align6.hpp"
#include "AssemblerOptions.hpp"
#include "compressAlignment.hpp"
#include "performanceLog.hpp"
#include "ProjectedAlignment.hpp"
#include "Reads.hpp"
#include "span.hpp"
#include "timestamp.hpp"
using namespace dinara;

// Standard libraries.
#include "chrono.hpp"
#include <cmath>
#include "iterator.hpp"
#include "tuple.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>

// hifiasm in-memory overlap bridge (hifiasm_overlap_t).
#include "hifiasm_overlaps.h"



// Compute a marker alignment of two oriented reads.
void Assembler::alignOrientedReads(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    size_t maxSkip,     // Maximum ordinal skip allowed.
    size_t maxDrift,    // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency
)
{
    alignOrientedReads(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        maxSkip, maxDrift, maxMarkerFrequency
        );
}
void Assembler::alignOrientedReads(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    size_t maxSkip, // Maximum ordinal skip allowed.
    size_t maxDrift,    // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency
)
{
    reads->checkReadsAreOpen();
    reads->checkReadNamesAreOpen();
    checkMarkersAreOpen();


    // Get the markers sorted by kmerId.
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    vector<MarkerWithOrdinal> markers1SortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId[0]);
    getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId[1]);

    // Call the lower level function.
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    const bool debug = true;
    alignOrientedReads(
        markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);

    // Compute the AlignmentInfo.
    uint32_t leftTrim;
    uint32_t rightTrim;
    tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();
    cout << orientedReadId0 << " has " << reads->getRead(orientedReadId0.getReadId()).baseCount;
    cout << " bases and " << markersSortedByKmerId[0].size() << " markers." << endl;
    cout << orientedReadId1 << " has " << reads->getRead(orientedReadId1.getReadId()).baseCount;
    cout << " bases and " << markersSortedByKmerId[1].size() << " markers." << endl;
    cout << "The alignment has " << alignmentInfo.markerCount;
    cout << " markers. Left trim " << leftTrim;
    cout << " markers, right trim " << rightTrim << " markers." << endl;

#if 0
    // For convenience, also write the two oriented reads.
    ofstream fasta("AlignedOrientedReads.fasta");
    writeOrientedRead(orientedReadId0, fasta);
    writeOrientedRead(orientedReadId1, fasta);
#endif
}



// This lower level version takes as input vectors of
// markers already sorted by kmerId.
void Assembler::alignOrientedReads(
    const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
    size_t maxSkip,             // Maximum ordinal skip allowed.
    size_t maxDrift,            // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency
)
{
    // Compute the alignment.
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    const bool debug = true;
    alignOrientedReads(
        markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);
}



void Assembler::alignOrientedReads(
    const array<vector<MarkerWithOrdinal>, 2>& markersSortedByKmerId,
    size_t maxSkip,             // Maximum ordinal skip allowed.
    size_t maxDrift,            // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency,
    bool debug,
    AlignmentGraph& graph,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo
)
{
    align(markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);
}


bool Assembler::computeAlignmentParity(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    Alignment& alignment
)
{
    // Parity Alignment Reconstruction
    // Uses actual assembler parameters to reconstruct exact alignment path.
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();

    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId[0]);
    getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId[1]);

    AlignmentGraph graph;
    AlignmentInfo alignmentInfo;
    const bool debug = false;
    
    // Use stored alignment parameters
    size_t maxSkip = assemblerInfo->actualMaxSkip > 0 ? assemblerInfo->actualMaxSkip : 30; // Default fallback
    size_t maxDrift = assemblerInfo->actualMaxDrift > 0 ? assemblerInfo->actualMaxDrift : 15;
    uint32_t maxMarkerFrequency = 100; // Standard

    align(
        markersSortedByKmerId,
        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);
        
    return alignmentInfo.markerCount > 0;
}



// Compute marker alignments of an oriented read with all reads
// for which we have an alignment.
void Assembler::alignOverlappingOrientedReads(
    ReadId readId, Strand strand,
    size_t maxSkip,                 // Maximum ordinal skip allowed.
    size_t maxDrift,                // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency,
    size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
    size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
)
{
    alignOverlappingOrientedReads(
        OrientedReadId(readId, strand),
        maxSkip, maxDrift, maxMarkerFrequency, minAlignedMarkerCount, maxTrim);
}



void Assembler::alignOverlappingOrientedReads(
    OrientedReadId orientedReadId0,
    size_t maxSkip,                 // Maximum ordinal skip allowed.
    size_t maxDrift,                // Maximum ordinal drift allowed.
    uint32_t maxMarkerFrequency,
    size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
    size_t maxTrim                  // Maximum trim (number of markers) allowed in an alignment.
)
{
    // Check that we have what we need.
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Get the markers for orientedReadId0.
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId[0]);

    // Loop over all alignments involving this oriented read.
    size_t goodAlignmentCount = 0;
    for(const uint64_t i: alignmentTable[orientedReadId0.getValue()]) {
        const AlignmentData& ad = alignmentData[i];

        // Get the other oriented read involved in this alignment.
        const OrientedReadId orientedReadId1 = ad.getOther(orientedReadId0);

        // Get the markers for orientedReadId1.
        getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId[1]);

        // Compute the alignment.
        AlignmentGraph graph;
        Alignment alignment;
        AlignmentInfo alignmentInfo;
        const bool debug = false;
        alignOrientedReads(
            markersSortedByKmerId,
            maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);

        uint32_t leftTrim;
        uint32_t rightTrim;
        tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();

        cout << orientedReadId0 << " " << orientedReadId1 << " " << alignmentInfo.markerCount;
        if(alignmentInfo.markerCount) {
            cout << " " << leftTrim << " " << rightTrim;
            if(alignmentInfo.markerCount >= minAlignedMarkerCount && leftTrim<=maxTrim && rightTrim<=maxTrim) {
                cout << " good";
                ++goodAlignmentCount;
            }
        }
        cout << endl;

    }
    cout << "Found " << goodAlignmentCount << " alignments out of ";
    cout << alignmentTable[orientedReadId0.getValue()].size() << "." << endl;

}



// Minimum PAF alignment block length (column 11) for an overlap to be kept.
static constexpr uint64_t pafMinAlignmentBlockLength = 200;

// Deduplicate the collected overlap entries (one entry per (read pair, strand),
// keeping the overlap with the highest hifiasm chain DP score) and publish them
// into alignmentCandidates.candidates. Returns the number of duplicate entries
// that were merged away. The caller must have created
// alignmentCandidates.candidates.
uint64_t Assembler::publishPafEntries(vector<PafEntry>& entries)
{
    const size_t rawCount = entries.size();

    // One entry per (read pair, strand), keeping the overlap with the highest
    // sharedSeedScore (hifiasm's oreg_ss_lt selection). A pair overlapping in
    // both orientations keeps up to two entries.
    dedupPafEntriesKeepBestScore(entries);
    const uint64_t duplicateCount = uint64_t(rawCount - entries.size());

    // Publish (deterministic: entries are in ascending key order, same-strand
    // before reverse within a key). Each surviving entry becomes one oriented
    // candidate.
    for (const PafEntry& e : entries) {
        const ReadId readId0 = ReadId(e.key >> 32);
        const ReadId readId1 = ReadId(e.key & 0xffffffffULL);
        alignmentCandidates.candidates.push_back(
            OrientedReadPair(readId0, readId1, e.iv.isSameStrand));
    }

    alignmentCandidates.unreserve();
    return duplicateCount;
}


// Import alignment candidates directly from hifiasm's in-memory overlaps,
// avoiding the PAF file round-trip. `overlaps` and the name table come from
// hifiasm_detect_overlaps_mem(): each overlap's q_id/t_id index into the name
// table, and read ids are resolved BY NAME (via reads->getReadId) exactly as
// the PAF-file path does, so load-order differences do not matter.
void Assembler::importAlignmentCandidatesFromMemory(
    const hifiasm_overlap_t* overlaps,
    uint64_t overlapCount,
    const char* names,
    const uint64_t* nameOffsets,
    uint64_t readCountFromHifiasm,
    const uint16_t* cigar,
    uint64_t cigarLen,
    const uint64_t* chain,
    uint64_t chainLen,
    uint64_t threadCount,
    uint32_t minOverlapLength,
    uint32_t maxEndFuzz)
{
    cout << timestamp << "Importing " << overlapCount
         << " hifiasm overlaps from memory..." << endl;
    const auto tBegin = steady_clock::now();

    alignmentCandidates.candidates.createNew(largeDataName("AlignmentCandidates"), largeDataPageSize);
    hifiasmImportedCigarStore.clear();

    if(threadCount == 0) threadCount = std::thread::hardware_concurrency();
    if(threadCount == 0) threadCount = 1;

    // Resolve hifiasm read index -> dinara ReadId once, by name. Names in the
    // table are not individually NUL-terminated: read i occupies the byte range
    // [nameOffsets[i], nameOffsets[i+1]).
    vector<ReadId> hifiToDinara(readCountFromHifiasm, invalidReadId);
    for(uint64_t i = 0; i < readCountFromHifiasm; i++) {
        const uint64_t b = nameOffsets[i];
        const uint64_t e = nameOffsets[i + 1];
        hifiToDinara[i] = reads->getReadId(span<const char>(names + b, size_t(e - b)));
    }

    // Overlap length + dovetail gates applied to hifiasm's candidate overlaps
    // before they become dinara candidates. hifiasm overlap
    // coordinates (q_start/q_end, t_start/t_end) are RAW forward-strand read
    // coordinates (the PAF path reports them alongside the raw read length), so
    // these are base-count tests:
    //   - kMinOverlapLen (OverlapCandidates.minOverlapLength): keep only when
    //     BOTH the query and target spans reach this many bases.
    //   - kDovetailHang (OverlapCandidates.maxEndFuzz): keep only when the
    //     aligned interval reaches within this many bases of an END of EACH
    //     read (dovetail test), dropping internal / repeat matches interior to
    //     both reads.
    // 0 disables a gate. The dovetail test fails OPEN when a read length is
    // unknown, so a name-resolution miss never silently drops a real overlap.
    const uint32_t kMinOverlapLen = minOverlapLength;
    const uint32_t kDovetailHang  = maxEndFuzz;
    cout << timestamp << "Candidate filters: minOverlapLength="
         << kMinOverlapLen << " (min of query/target span), maxEndFuzz="
         << kDovetailHang << " (dovetail hang); 0 disables." << endl;
    auto nearEnd = [](uint32_t start, uint32_t end, uint32_t len,
                      uint32_t hang) -> bool {
        return start < hang || (len > 0 && end + hang > len);
    };

    // Build PafEntry records in parallel, mirroring the per-record filtering of
    // the PAF importer (block length floor, distinct, non-palindromic).
    const uint64_t batch = std::max<uint64_t>(1, overlapCount / threadCount);
    vector<vector<PafEntry>> threadEntries(threadCount);
    vector<uint64_t> threadKept(threadCount, 0);
    vector<std::thread> workers;
    workers.reserve(threadCount);
    for(uint64_t ti = 0; ti < threadCount; ti++) {
        workers.emplace_back([&, ti]() {
            auto& out = threadEntries[ti];
            uint64_t kept = 0;
            for(uint64_t i = ti * batch;
                i < overlapCount && (ti == threadCount - 1 || i < (ti + 1) * batch);
                i++) {
                const hifiasm_overlap_t& o = overlaps[i];
                if(o.block_len < pafMinAlignmentBlockLength) continue;
                if(o.q_id >= readCountFromHifiasm || o.t_id >= readCountFromHifiasm) continue;
                const ReadId readId0 = hifiToDinara[o.q_id];
                const ReadId readId1 = hifiToDinara[o.t_id];
                const bool valid =
                    readId0 != invalidReadId &&
                    readId1 != invalidReadId &&
                    readId0 != readId1 &&
                    !reads->getFlags(readId0).isPalindromic &&
                    !reads->getFlags(readId1).isPalindromic;
                if(!valid) continue;

                // Minimum overlap length: both spans must reach kMinOverlapLen.
                if(kMinOverlapLen > 0) {
                    const uint32_t qSpan = (o.q_end > o.q_start) ? (o.q_end - o.q_start) : 0;
                    const uint32_t tSpan = (o.t_end > o.t_start) ? (o.t_end - o.t_start) : 0;
                    if(qSpan < kMinOverlapLen || tSpan < kMinOverlapLen) continue;
                }

                // Dovetail gate: require end-proximity on BOTH reads. Fail open
                // when a raw length is unknown (0).
                if(kDovetailHang > 0) {
                    const uint32_t qLen = uint32_t(reads->getReadRawSequenceLength(readId0));
                    const uint32_t tLen = uint32_t(reads->getReadRawSequenceLength(readId1));
                    if(qLen > 0 && tLen > 0) {
                        const bool qNear = nearEnd(o.q_start, o.q_end, qLen, kDovetailHang);
                        const bool tNear = nearEnd(o.t_start, o.t_end, tLen, kDovetailHang);
                        if(!(qNear && tNear)) continue;
                    }
                }
                PafEntry entry = makePafEntry(
                    readId0, readId1,
                    o.q_start, o.q_end,
                    o.t_start, o.t_end,
                    o.block_len, o.shared_seed, o.is_same_strand != 0);
                // Remember which overlap this entry came from so the CIGAR of
                // the entry that survives dedup can be recovered below.
                entry.sourceIndex = i;
                out.push_back(entry);
                ++kept;
            }
            threadKept[ti] = kept;
        });
    }
    for(auto& t : workers) t.join();

    // Concatenate thread-local entries (order-independent: dedup sorts).
    size_t total = 0;
    uint64_t keptTotal = 0;
    for(uint64_t ti = 0; ti < threadCount; ti++) {
        total += threadEntries[ti].size();
        keptTotal += threadKept[ti];
    }
    vector<PafEntry> entries;
    entries.reserve(total);
    for(uint64_t ti = 0; ti < threadCount; ti++) {
        auto& te = threadEntries[ti];
        entries.insert(entries.end(), te.begin(), te.end());
        vector<PafEntry>().swap(te);
    }

    const uint64_t duplicateCount = publishPafEntries(entries);

    // publishPafEntries deduped `entries` in place, so it now holds exactly the
    // overlaps that became candidates (one per key+strand, highest chain DP
    // score). Copy
    // each survivor's hifiasm CIGAR into the imported-CIGAR store, keyed by
    // canonical read pair, so computeBaseAlignmentsAndStore can reuse hifiasm's
    // base alignment instead of recomputing it. Stored in the native
    // (query,target) alignment frame; reframing to read0/read1 happens at use.
    if(cigar != nullptr) {
        uint64_t cigarOverlaps = 0;
        hifiasmImportedCigarStore.reserve(entries.size(), cigarLen);
        if(chain != nullptr) hifiasmImportedCigarStore.reserveChain(chainLen);
        for(const PafEntry& e : entries) {
            if(e.sourceIndex == uint64_t(-1)) continue;
            const hifiasm_overlap_t& o = overlaps[e.sourceIndex];
            if(o.cigar_len == 0) continue;
            if(o.cigar_offset + o.cigar_len > cigarLen) continue; // defensive
            const ReadId readId0 = hifiToDinara[o.q_id];
            const ReadId readId1 = hifiToDinara[o.t_id];
            hifiasmImportedCigarStore.add(
                e.key, o.is_same_strand != 0,
                span<const uint16_t>(cigar + o.cigar_offset, size_t(o.cigar_len)),
                uint32_t(readId0), uint32_t(readId1),
                o.q_start, o.q_end, o.t_start, o.t_end);
            ++cigarOverlaps;
            if(chain != nullptr && o.chain_len > 0 &&
               o.chain_offset + o.chain_len <= chainLen) {
                hifiasmImportedCigarStore.addChain(
                    e.key, o.is_same_strand != 0,
                    span<const uint64_t>(chain + o.chain_offset, size_t(o.chain_len)));
            }
        }
        cout << timestamp << "Imported hifiasm CIGARs for " << cigarOverlaps
             << " overlaps (" << cigarLen << " tokens)." << endl;
    } else {
        // Interval-only import: hifiasm's base CIGAR is discarded (cigar ==
        // nullptr). Store each surviving overlap's interval with ZERO CIGAR
        // tokens plus its native dense chain, so computeBaseAlignmentsAndStore
        // can find the pair, map the chain anchors to marker ordinals
        // (mapNativeChainToOrdinals), and build the CIGAR per-segment with A*PA2
        // (constructQuickRawSparse). Without this record the pair has no chain
        // and the empty-ordinals guard downstream would drop every overlap.
        uint64_t intervalOverlaps = 0;
        uint64_t chainOverlaps = 0, chainAnchorsTotal = 0;
        hifiasmImportedCigarStore.reserve(entries.size(), 0);
        if(chain != nullptr) hifiasmImportedCigarStore.reserveChain(chainLen);
        for(const PafEntry& e : entries) {
            if(e.sourceIndex == uint64_t(-1)) continue;
            const hifiasm_overlap_t& o = overlaps[e.sourceIndex];
            const ReadId readId0 = hifiToDinara[o.q_id];
            const ReadId readId1 = hifiToDinara[o.t_id];
            hifiasmImportedCigarStore.add(
                e.key, o.is_same_strand != 0,
                span<const uint16_t>(),   // no tokens: interval-only record
                uint32_t(readId0), uint32_t(readId1),
                o.q_start, o.q_end, o.t_start, o.t_end);
            ++intervalOverlaps;
            // Attach hifiasm's native dense chain for this pair (query-forward /
            // target-alignment frame). Consumed by mapNativeChainToOrdinals.
            if(chain != nullptr && o.chain_len > 0 &&
               o.chain_offset + o.chain_len <= chainLen) {
                hifiasmImportedCigarStore.addChain(
                    e.key, o.is_same_strand != 0,
                    span<const uint64_t>(chain + o.chain_offset, size_t(o.chain_len)));
                ++chainOverlaps;
                chainAnchorsTotal += o.chain_len;
            }
        }
        cout << timestamp << "Imported interval-only overlaps (no CIGAR) for "
             << intervalOverlaps << " pairs." << endl;
        if(chain != nullptr) {
            cout << timestamp << "Imported native chain for " << chainOverlaps
                 << " pairs (" << chainAnchorsTotal << " anchors)." << endl;
        }
    }

    const double seconds = 1.e-9 * double(std::chrono::duration_cast<std::chrono::nanoseconds>(
        steady_clock::now() - tBegin).count());
    cout << timestamp << "Imported " << keptTotal << " hifiasm overlaps ("
         << duplicateCount << " duplicate pairs merged) in " << seconds
         << " s using " << threadCount << " threads." << endl;
    cout << timestamp << "Total candidates: " << alignmentCandidates.candidates.size() << endl;
}


// Compute alignments.
void Assembler::computeAlignments(

    const AlignOptions& alignOptions,

    // Number of threads. If zero, a number of threads equal to
    // the number of virtual processors is used.
    uint64_t threadCount
)
{

    const auto tBegin = steady_clock::now();


    performanceLog << timestamp << "Begin computing alignments for ";
    performanceLog << alignmentCandidates.candidates.size() << " alignment candidates." << endl;

    // Check that we have what we need.
    reads->checkReadsAreOpen();
    // Note: kmerChecker is not required when using SIMD minimizers for marker generation.
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Store parameters so they are accessible to the threads.
    auto& data = computeAlignmentsData;
    data.alignOptions = &alignOptions;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // For alignment method 4, compute sorted markers.
    if(alignOptions.alignMethod == 4) {
        cout << timestamp << "Computing sorted markers." << endl;
        computeSortedMarkers(threadCount);
    }

    // For alignment method 5, compute low frequency markers.
    if(alignOptions.alignMethod == 5) {
        cout << timestamp << "Computing unique markers." << endl;
        computeLowFrequencyMarkers(1, threadCount);
    }

    // For alignment method 6, compute Align6Markers.
    if(alignOptions.alignMethod == 6) {
        performanceLog << timestamp << "Computing Align6Markers." << endl;
        computeAlign6Markers(threadCount);
        performanceLog << timestamp << "Done computing Align6Markers." << endl;
    }

    // Pick the batch size for computing alignments.
    size_t batchSize = 10;
    if(batchSize > alignmentCandidates.candidates.size()/threadCount) {
        batchSize = alignmentCandidates.candidates.size()/threadCount;
    }
    if(batchSize == 0) {
        batchSize = 1;
    }


    // Compute the alignments.
    data.threadAlignmentData.resize(threadCount);
    data.threadCompressedAlignments.resize(threadCount);
    
    // Always resize these to avoid SIGSEGV during aggregation, even if not used.
    data.threadProjectedAlignmentTime.assign(threadCount, 0.0);
    data.threadCollectionTime.assign(threadCount, 0.0);
    data.threadFilteredByErrorRate.assign(threadCount, 0);
    data.threadFilteredByErrorRateGap.assign(threadCount, 0);
    data.threadFilteredByGapCount.assign(threadCount, 0);

    if (assemblerInfo->readGraphCreationMethod == 5) {
        data.threadVariantClusteringPositionPairs.resize(threadCount);
    }
    
    performanceLog << timestamp << "Alignment computation begins." << endl;
    cout << timestamp << "Alignment computation begins." << endl;
    setupLoadBalancing(alignmentCandidates.candidates.size(), batchSize);
    runThreads(&Assembler::computeAlignmentsThreadFunction, threadCount);
    performanceLog << timestamp << "Alignment computation completed." << endl;
    cout << timestamp << "Alignment computation completed." << endl;

    // Store the alignments found by each thread.
    performanceLog << timestamp << "Storing the alignment found by each thread." << endl;

    // Pre-size global containers to avoid repeated remaps/reallocations.
    uint64_t totalAlignments = 0;
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        totalAlignments += data.threadAlignmentData[threadId].size();
    }

    alignmentData.createNew(largeDataName("AlignmentData"), largeDataPageSize, 0, totalAlignments);
    compressedAlignments.createNew(largeDataName("CompressedAlignments"), largeDataPageSize);

    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
        const size_t idShift = alignmentData.size(); // Current global count serves as offset
        for(const AlignmentData& ad: threadAlignmentData) {
            alignmentData.push_back(ad);
            alignmentData.back().info.alignmentId += idShift;
        }

        const auto threadCompressedAlignments = data.threadCompressedAlignments[threadId];
        const auto size = threadCompressedAlignments->size();
        for(size_t i=0; i<size; i++) {
            compressedAlignments.appendVector(
                (*threadCompressedAlignments)[i].begin(),
                (*threadCompressedAlignments)[i].end()
            );
        }



    }

    // Release unused allocated memory.
    alignmentData.unreserve();
    compressedAlignments.unreserve();





    // Aggregate timing statistics from all threads
    double totalProjectedAlignmentTime = 0.0;
    double totalCollectionTime = 0.0;
    [[maybe_unused]] double estimatedProjectedWallTime = 0.0;
    [[maybe_unused]] double estimatedCollectionWallTime = 0.0;
    
    if (assemblerInfo->readGraphCreationMethod == 5) {
        for (size_t i = 0; i < threadCount; i++) {
            totalProjectedAlignmentTime += data.threadProjectedAlignmentTime[i];
            totalCollectionTime += data.threadCollectionTime[i];
        }
    
        // Estimate wall-clock time by dividing summed thread time by thread count
        estimatedProjectedWallTime = totalProjectedAlignmentTime / threadCount;
        estimatedCollectionWallTime = totalCollectionTime / threadCount;
    }

    if (assemblerInfo->readGraphCreationMethod == 5) {
        uint64_t totalFilteredByErrorRate = 0;
        uint64_t totalFilteredByErrorRateGap = 0;
        uint64_t totalFilteredByGapCount = 0;
        for (size_t i = 0; i < threadCount; i++) {
            totalFilteredByErrorRate += data.threadFilteredByErrorRate[i];
            totalFilteredByErrorRateGap += data.threadFilteredByErrorRateGap[i];
            totalFilteredByGapCount += data.threadFilteredByGapCount[i];
        }
        cout << "Time spent in projected alignment construction (all threads): " << totalProjectedAlignmentTime << " s" << endl;
        cout << "Time spent collecting variant clustering position pairs (all threads): " << totalCollectionTime << " s" << endl;
        cout << "Alignments filtered by error rate (> 0.07): " << totalFilteredByErrorRate << endl;
        cout << "Alignments filtered by gap error rate (> 0.006): " << totalFilteredByErrorRateGap << endl;
        cout << "Alignments filtered by gap count (> 64): " << totalFilteredByGapCount << endl;
    }
    if (assemblerInfo->readGraphCreationMethod == 5) {
#if DINARA_ENABLE_VARIANT_CLUSTERING
        variantClusteringProjectedAlignmentTime = estimatedProjectedWallTime;
        variantClusteringCollectionTime = estimatedCollectionWallTime;

        cout << "\nVariant clustering collection timing (estimated wall-clock):" << endl;
        cout << "  Projected alignment: ~" << estimatedProjectedWallTime << " s" << endl;
        cout << "  Position pair collection: ~" << estimatedCollectionWallTime << " s" << endl;
        performanceLog << timestamp << "Projected alignment (estimated wall-clock): " << estimatedProjectedWallTime << " s" << endl;
        performanceLog << timestamp << "Collection (estimated wall-clock): " << estimatedCollectionWallTime << " s" << endl;

        // Store position pairs collected by each thread
        performanceLog << timestamp << "Storing position pairs for variant clustering." << endl;
        cout << timestamp << "Storing position pairs for variant clustering." << endl;
        storeVariantClusteringPositionPairs(threadCount, data);
        performanceLog << timestamp << "Done storing position pairs for variant clustering." << endl;
        cout << timestamp << "Done storing position pairs for variant clustering." << endl;
#else
        throw runtime_error("readGraphCreationMethod=5 requires DINARA_ENABLE_VARIANT_CLUSTERING=ON.");
#endif
    }



    // Cleanup.
    if(alignOptions.alignMethod == 4) {
        sortedMarkers->remove();
    }
    if(alignOptions.alignMethod == 5) {
        lowFrequencyMarkers->remove();
    }
    if(alignOptions.alignMethod == 6) {
        align6Markers->remove();
    }

    cout << "Found and stored " << alignmentData.size() << " good alignments." << endl;
    performanceLog << timestamp << "Creating alignment table." << endl;
    computeAlignmentTable();

    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    performanceLog << timestamp << "Computation of alignments ";
    performanceLog << "completed in " << tTotal << " s." << endl;

    performanceLog << timestamp;
}



void Assembler::computeAlignmentsThreadFunction(size_t threadId)
{

    array<OrientedReadId, 2> orientedReadIds;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    string compressedAlignment;

    const bool debug = false;
    auto& data = computeAlignmentsData;
    const size_t alignmentMethod = data.alignOptions->alignMethod;
    const uint32_t maxMarkerFrequency = data.alignOptions->maxMarkerFrequency;
    const size_t maxSkip = data.alignOptions->maxSkip;
    const size_t maxDrift = data.alignOptions->maxDrift;
    const size_t minAlignedMarkerCount = data.alignOptions->minAlignedMarkerCount;
    const double minAlignedFraction = data.alignOptions->minAlignedFraction;
    const size_t maxTrim = data.alignOptions->maxTrim;
    const int matchScore = data.alignOptions->matchScore;
    const int mismatchScore = data.alignOptions->mismatchScore;
    const int gapScore = data.alignOptions->gapScore;
    const double downsamplingFactor = data.alignOptions->downsamplingFactor;
    const int bandExtend = data.alignOptions->bandExtend;
    const int maxBand = data.alignOptions->maxBand;
    const double align5DriftRateTolerance = data.alignOptions->align5DriftRateTolerance;
    const uint64_t align5MinBandExtend = data.alignOptions->align5MinBandExtend;


    // Align4-specific items.
    Align4::Options align4Options;
    MemoryMapped::ByteAllocator byteAllocator;
    if(alignmentMethod == 4) {
        align4Options.deltaX = data.alignOptions->align4DeltaX;
        align4Options.deltaY = data.alignOptions->align4DeltaY;
        align4Options.minEntryCountPerCell = data.alignOptions->align4MinEntryCountPerCell;
        align4Options.maxDistanceFromBoundary = data.alignOptions->align4MaxDistanceFromBoundary;
        align4Options.minAlignedMarkerCount = minAlignedMarkerCount;
        align4Options.minAlignedFraction = minAlignedFraction;
        align4Options.maxSkip = maxSkip;
        align4Options.maxDrift = maxDrift;
        align4Options.maxTrim = maxTrim;
        align4Options.maxBand = maxBand;
        align4Options.matchScore = matchScore;
        align4Options.mismatchScore = mismatchScore;
        align4Options.gapScore = gapScore;
        byteAllocator.createNew(
            largeDataName("tmp-ByteAllocator-" + to_string(threadId)),
            largeDataPageSize, 2ULL * 1024 * 1024 * 1024);
    }

    ofstream nullStream;
    Align6 align6(
        assemblerInfo->k,
        maxSkip,
        maxDrift,
        data.alignOptions->align6Options,
        assemblerInfo->kmerDistributionInfo,
        nullStream);
    if((threadId == 0) and (alignmentMethod == 6)) {
        std::lock_guard<std::mutex> lock(mutex);
        align6.writeGlobalFrequencyCriteria(cout);
    }

    vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
    
    shared_ptr< MemoryMapped::VectorOfVectors<char, uint64_t> > thisThreadCompressedAlignmentsPointer =
        make_shared< MemoryMapped::VectorOfVectors<char, uint64_t> >();
    data.threadCompressedAlignments[threadId] = thisThreadCompressedAlignmentsPointer;
    auto& thisThreadCompressedAlignments = *thisThreadCompressedAlignmentsPointer;
    thisThreadCompressedAlignments.createNew(
        largeDataName("tmp-ThreadGlobalCompressedAlignments-" + to_string(threadId)),
        largeDataPageSize);

    if (assemblerInfo->readGraphCreationMethod == 5) {
        // Create the vector to contain the position pairs collected by this thread.
        shared_ptr< MemoryMapped::Vector< pair<OrientedReadId, uint32_t> > > thisThreadVariantClusteringPositionPairsPointer =
            make_shared< MemoryMapped::Vector< pair<OrientedReadId, uint32_t> > >();
        data.threadVariantClusteringPositionPairs[threadId] = thisThreadVariantClusteringPositionPairsPointer;
        auto& thisThreadVariantClusteringPositionPairs = *thisThreadVariantClusteringPositionPairsPointer;
        thisThreadVariantClusteringPositionPairs.createNew(
            largeDataName("tmp-ThreadVariantClusteringPositionPairs-" + to_string(threadId)),
            largeDataPageSize);
    }



#if 0
    // A vector to store the time taken to compute each alignment.
    vector< pair<uint64_t, double> > elapsedTime;
#endif

    const uint64_t messageFrequency = max(1UL, min(1000000UL, alignmentCandidates.candidates.size()/20));

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        for(size_t i=begin; i!=end; i++) {
#if 0
            const auto steadyClock0 = std::chrono::steady_clock::now();
#endif

            const OrientedReadPair& candidate = alignmentCandidates.candidates[i];
            DINARA_ASSERT(candidate.readIds[0] < candidate.readIds[1]);

            // Get the oriented read ids, with the first one on strand 0.
            orientedReadIds[0] = OrientedReadId(candidate.readIds[0], 0);
            orientedReadIds[1] = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);

            if((i % messageFrequency) == 0){
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << "Working on alignment " << i;
                performanceLog << " of " << alignmentCandidates.candidates.size();
                // performanceLog << ": " << orientedReadIds[0] << " " << orientedReadIds[1];
                performanceLog << endl;
            }


            // Compute the alignment.
            const bool precomputedUsed = false;
            try {
                if(alignmentMethod == 0) {

                    // Get the markers for the two oriented reads in this candidate.
                    for(size_t j=0; j<2; j++) {
                        getMarkersSortedByKmerId(orientedReadIds[j], markersSortedByKmerId[j]);
                    }

                    // Compute the Alignment.
                    alignOrientedReads(
                        markersSortedByKmerId,
                        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);

                } else if(alignmentMethod == 1) {
                    alignOrientedReads1(orientedReadIds[0], orientedReadIds[1],
                        matchScore, mismatchScore, gapScore,
                        alignment, alignmentInfo);
                } else if(alignmentMethod == 3) {
                    alignOrientedReads3(orientedReadIds[0], orientedReadIds[1],
                        matchScore, mismatchScore, gapScore,
                        downsamplingFactor, bandExtend, maxBand,
                        alignment, alignmentInfo);
                } else if(alignmentMethod == 4) {
                    alignOrientedReads4(orientedReadIds[0], orientedReadIds[1],
                        align4Options,
                        byteAllocator,
                        alignment, alignmentInfo,
                        false);
                    DINARA_ASSERT(byteAllocator.isEmpty());
                } else if(alignmentMethod == 5) {
                    ofstream nullStream;
                    alignOrientedReads5(orientedReadIds[0], orientedReadIds[1],
                        matchScore, mismatchScore, gapScore,
                        align5DriftRateTolerance, align5MinBandExtend,
                        alignment, alignmentInfo,
                        nullStream);
                } else if(alignmentMethod == 6) {
                    ofstream nullStream;
                    alignOrientedReads6(orientedReadIds[0], orientedReadIds[1],
                        alignment, alignmentInfo, align6);
                } else {
                    DINARA_ASSERT(0);
                }
            } catch (std::exception& e) {
                std::lock_guard<std::mutex> lock(mutex);
                cout <<
                    "An error occurred while computing a marker alignment "
                    " of oriented reads " << orientedReadIds[0] << " and " << orientedReadIds[1] <<
                    ". This alignment candidate will be skipped. Error description is: " <<
                    e.what() << endl;
                continue;
            } catch(...) {
                std::lock_guard<std::mutex> lock(mutex);
                cout <<
                    "An error occurred while computing a marker alignment "
                    " of oriented reads " << orientedReadIds[0] << " and " << orientedReadIds[1] <<
                    ". This alignment candidate will be skipped. " << endl;
                continue;
            }

#if 0
            const auto steadyClock1 = std::chrono::steady_clock::now();
            const auto deltaT = 1.e-9 * double(std::chrono::duration_cast<std::chrono::nanoseconds>(steadyClock1 - steadyClock0).count());
            elapsedTime.push_back({i, deltaT});
#endif

            // If the alignment has too few markers, skip it.
            if(alignment.ordinals.size() < minAlignedMarkerCount) {
                continue;
            }

            // If the aligned fraction is too small, skip it.
            if(min(alignmentInfo.alignedFraction(0), alignmentInfo.alignedFraction(1)) < minAlignedFraction) {
                continue;
            }

            // If the alignment has too much trim, skip it.
            uint32_t leftTrim;
            uint32_t rightTrim;
            tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();
            if(leftTrim>maxTrim || rightTrim>maxTrim) {
                continue;
            }

            // For alignment methods other than method 0, we also need to check for
            // maxSip and maxDrift. Method 0 does that automatically.
            if(alignmentMethod != 0) {
                if(alignment.maxSkip() > maxSkip) {
                    continue;
                }
                if(alignment.maxDrift() > maxDrift) {
                    continue;
                }
            }

            // If getting here, this is a good alignment.

            // Compute projected alignment metrics.
            const auto tProjStart = steady_clock::now();
            const ProjectedAlignment projectedAlignment(
                *this,
                orientedReadIds,
                alignment,
                ProjectedAlignment::Method::QuickRawSparse,
                data.alignOptions->overlapDpMatchScore,
                data.alignOptions->overlapDpMismatchScore,
                data.alignOptions->overlapDpGapOpen1,
                data.alignOptions->overlapDpGapExtend1);
            const auto tProjEnd = steady_clock::now();
            
            alignmentInfo.errorRate = float(projectedAlignment.errorRate());
            alignmentInfo.mismatchCount = uint32_t(projectedAlignment.mismatchCount);
            alignmentInfo.errorRateGaps = float(projectedAlignment.errorRateGaps());
            alignmentInfo.gapCount = uint32_t(projectedAlignment.totalIndelBaseCount);
            alignmentInfo.gapEventCount = uint32_t(projectedAlignment.totalGapEventCount); // Transfer gap events
            alignmentInfo.dpScore = projectedAlignment.totalDpScore;
            

            data.threadProjectedAlignmentTime[threadId] += seconds(tProjEnd - tProjStart);

            const uint64_t ql = projectedAlignment.totalLength[0];
            const uint64_t tl = projectedAlignment.totalLength[1];
            
            const double errorRateThreshold = 0.07;
            const uint64_t totalErrors = uint64_t(projectedAlignment.mismatchCount) + uint64_t(projectedAlignment.totalIndelBaseCount);
            if ((totalErrors > (ql * errorRateThreshold)) || (totalErrors > (tl * errorRateThreshold))) {
                data.threadFilteredByErrorRate[threadId]++;
                continue;
            }

            // const double gapRateThreshold = 0.006;
            // const uint64_t totalGapCount = uint64_t(projectedAlignment.totalIndelBaseCount);
            // if ((totalGapCount > (ql * gapRateThreshold)) || (totalGapCount > (tl * gapRateThreshold))) {
            //     data.threadFilteredByErrorRateGap[threadId]++;
            //     continue;
            // }

            // // Collect position pairs for variant clustering
            // // Only collect those with SNP differences (No indels)
            // const auto tCollectStart = steady_clock::now();
            // collectVariantClusteringPositionPairs(
            //     projectedAlignment,
            //     orientedReadIds,
            //     *data.threadVariantClusteringPositionPairs[threadId]
            // );
            // const auto tCollectEnd = steady_clock::now();
            // data.threadCollectionTime[threadId] += seconds(tCollectEnd - tCollectStart);

            // cout << orientedReadIds[0] << " " << orientedReadIds[1] << " good." << endl;
            AlignmentData thisAlignmentData(candidate, alignmentInfo);
            
            // Calculate coordinates efficiently using direct marker access
            if (precomputedUsed) {
                // Option A: Use stored coordinates from Alignment (computed during InvertedIndex chaining)
                // This avoids re-fetching markers and re-computing extension logic.
                // Parity logic (including -1 for Diff) was already applied in InvertedIndexFinder.
                thisAlignmentData.qs = alignment.qs;
                thisAlignmentData.qe = alignment.qe;
                thisAlignmentData.ts = alignment.ts;
                thisAlignmentData.te = alignment.te;
            } else {
                // Fallback: Compute coordinates from markers (Legacy path for other aligners)
                const auto markers0 = (*markers)[orientedReadIds[0].getValue()];
                uint32_t qs_marker = markers0[alignmentInfo.data[0].firstOrdinal].position;
                uint32_t qe_marker = markers0[alignmentInfo.data[0].lastOrdinal].position + uint32_t(assemblerInfo->k);
    
                const auto markers1 = (*markers)[orientedReadIds[1].getValue()];
                uint32_t ts_marker = markers1[alignmentInfo.data[1].firstOrdinal].position;
                uint32_t te_marker = markers1[alignmentInfo.data[1].lastOrdinal].position + uint32_t(assemblerInfo->k);
    
                // Read lengths for coordinate extension
                const uint64_t len0 = reads->getRead(orientedReadIds[0].getReadId()).baseCount;
                const uint64_t len1 = reads->getRead(orientedReadIds[1].getReadId()).baseCount;
                
                // --- Hifiasm Parity: Extend coordinates to read boundaries ---
                // This matches append_inexact_overlap_region_alloc in Hash_Table.cpp:374-398
                // The logic extends alignment coordinates so they reach read tips symmetrically.
                
                uint32_t qs_ext = qs_marker;
                uint32_t qe_ext = qe_marker;
                uint32_t ts_ext = ts_marker;
                uint32_t te_ext = te_marker;
                
                // Extend start: whoever has shorter overhang extends to 0, other adjusts proportionally
                if (qs_ext <= ts_ext) {
                    ts_ext = ts_ext - qs_ext;  // Reduce ts by qs amount
                    qs_ext = 0;                 // qs goes to 0
                } else {
                    qs_ext = qs_ext - ts_ext;  // Reduce qs by ts amount
                    ts_ext = 0;                 // ts goes to 0
                }
                
                // Extend end: whoever has shorter remaining distance extends to read length
                // Hifiasm uses closed coordinates: x_right_length = xLen - x_pos_e - 1
                // Dinara uses half-open [start, end): equivalent is len0 - qe_ext
                // But we need to be careful about the semantics:
                // Hifiasm's x_pos_e is inclusive (last base), Dinara's qe_ext is exclusive (one past last)
                // So: Hifiasm's "remaining" = xLen - x_pos_e - 1 = positions after x_pos_e
                // Dinara's "remaining" = len0 - qe_ext = same thing (qe_ext points one past last aligned base)
                
                // For closed-to-half-open conversion: closed_end + 1 = half_open_end
                // So if Hifiasm has x_pos_e as last aligned base, Dinara has qe_ext = x_pos_e + 1
                // Thus: len0 - qe_ext = len0 - (x_pos_e + 1) = (len0 - 1) - x_pos_e = x_right_length ✓
                
                int64_t q_right = (int64_t)len0 - (int64_t)qe_ext;  // Safe signed arithmetic
                int64_t t_right = (int64_t)len1 - (int64_t)te_ext;
                
                if (q_right <= t_right) {
                    qe_ext = (uint32_t)len0;           // qe extends to read end (half-open)
                    te_ext = te_ext + (uint32_t)q_right;  // te extends by same amount
                } else {
                    te_ext = (uint32_t)len1;           // te extends to read end (half-open)
                    qe_ext = qe_ext + (uint32_t)t_right;  // qe extends by same amount
                }
                
                thisAlignmentData.qs = qs_ext;
                thisAlignmentData.qe = qe_ext;
                
                // --- Convert target coordinates to FORWARD STRAND (hifiasm convention) ---
                // When isSameStrand=false (reverse), the ts/te are on the reverse-complement strand.
                // We need to flip them to represent positions on the forward strand.
                if (!candidate.isSameStrand) {
                    // Convert reverse-complement half-open interval [ts_ext, te_ext) to forward [ts, te).
                    const auto p = dinara::rcIntervalToForward(uint32_t(len1), ts_ext, te_ext);
                    thisAlignmentData.ts = p.first;
                    thisAlignmentData.te = p.second;
                } else {
                    thisAlignmentData.ts = ts_ext;
                    thisAlignmentData.te = te_ext;
                }
            }

            // Check for large indels (>= 6 bases)
            thisAlignmentData.hasLargeIndel = projectedAlignment.hasLargeIndel;

            // Alignment covers an informative het site (counted per read perspective during EC parity)
            thisAlignmentData.informativeHetSiteCount0 = 0;
            thisAlignmentData.informativeHetSiteCount1 = 0;
            thisAlignmentData.informativeHetSiteScore = 0;
            
            // Alignment deletion reasons (none by default)
            thisAlignmentData.deleteReasons0 = AlignmentData::DeleteReasonNone;
            thisAlignmentData.deleteReasons1 = AlignmentData::DeleteReasonNone;

            // alignmentId is kept as the running per-thread alignment index
            // (later shifted to a global index during the merge above).
            thisAlignmentData.info.alignmentId = threadAlignmentData.size();

            // Store AlignmentData and the corresponding compressed alignment (same order).
            threadAlignmentData.push_back(thisAlignmentData);

            // Store the alignment in compressed form.
            dinara::compress(alignment, compressedAlignment);
            thisThreadCompressedAlignments.appendVector(
                compressedAlignment.c_str(),
                compressedAlignment.c_str() + compressedAlignment.size()
            );
            
        }
    }

    if(alignmentMethod == 4) {

        std::lock_guard<std::mutex> lock(mutex);
        cout << "Thread " << threadId << " byte allocator: " <<
            byteAllocator.getMaxAllocatedByteCount() << "/" <<
            2ULL * 1024 * 1024 * 1024 << endl;
    }

    thisThreadCompressedAlignments.unreserve();

#if 0
    // Write the elapsed time taken by each alignment.
    ofstream csv("AlignTime-Thread" + to_string(threadId) + ".csv");
    for(const auto& p: elapsedTime) {
        const uint64_t i = p.first;
        const OrientedReadPair& candidate = alignmentCandidates.candidates[i];
        DINARA_ASSERT(candidate.readIds[0] < candidate.readIds[1]);
        const OrientedReadId orientedReadId0 = OrientedReadId(candidate.readIds[0], 0);
        const OrientedReadId orientedReadId1 = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);
        const double t = p.second;
        csv << i << ",";
        csv << orientedReadId0 << ",";
        csv << orientedReadId1 << ",";
        csv << t << "\n";
    }
#endif
}



void Assembler::accessCompressedAlignments()
{
    compressedAlignments.accessExistingReadOnly(
        largeDataName("CompressedAlignments"));
}



// Compute alignmentTable from alignmentData.
// This could be made multithreaded if it becomes a bottleneck.
void Assembler::computeAlignmentTable()
{
    // Avoid rebuilding if already computed in this process.
    // Some pipelines call alignment computation routines more than once during experimentation;
    // rebuilding the alignment table is expensive and also fails if the memory-mapped object
    // is still open.
    if(alignmentTable.isOpen()) {
        performanceLog << timestamp << "Alignment table already exists - skipping recomputation." << endl;
        return;
    }

    alignmentTable.createNew(largeDataName("AlignmentTable"), largeDataPageSize);
    alignmentTable.beginPass1(ReadId(2 * reads->readCount()));
    for(const AlignmentData& ad: alignmentData) {
        const auto& readIds = ad.readIds;
        OrientedReadId orientedReadId0(readIds[0], 0);
        OrientedReadId orientedReadId1(readIds[1], ad.isSameStrand ? 0 : 1);
        alignmentTable.incrementCount(orientedReadId0.getValue());
        alignmentTable.incrementCount(orientedReadId1.getValue());
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        alignmentTable.incrementCount(orientedReadId0.getValue());
        alignmentTable.incrementCount(orientedReadId1.getValue());
    }
    alignmentTable.beginPass2();
    for(uint32_t i=0; i<alignmentData.size(); i++) {
        const AlignmentData& ad = alignmentData[i];
        const auto& readIds = ad.readIds;
        OrientedReadId orientedReadId0(readIds[0], 0);
        OrientedReadId orientedReadId1(readIds[1], ad.isSameStrand ? 0 : 1);
        alignmentTable.store(orientedReadId0.getValue(), i);
        alignmentTable.store(orientedReadId1.getValue(), i);
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        alignmentTable.store(orientedReadId0.getValue(), i);
        alignmentTable.store(orientedReadId1.getValue(), i);
    }
    alignmentTable.endPass2();



    // Sort each section of the alignment table by OrientedReadId.
    vector< pair<OrientedReadId, uint32_t> > v;
    for(ReadId readId0=0; readId0<reads->readCount(); readId0++) {
        for(Strand strand0=0; strand0<2; strand0++) {
            const OrientedReadId orientedReadId0(readId0, strand0);

            // Access the section of the alignment table for this oriented read.
            const span<uint32_t> alignmentTableSection =
                alignmentTable[orientedReadId0.getValue()];

            // Store pairs(OrientedReadId, alignmentIndex).
            v.clear();
            for(uint32_t alignmentIndex: alignmentTableSection) {
                const AlignmentData& alignment = alignmentData[alignmentIndex];
                const OrientedReadId orientedReadId1 = alignment.getOther(orientedReadId0);
                v.push_back(make_pair(orientedReadId1, alignmentIndex));
            }

            // Sort.
            sort(v.begin(), v.end());

            // Store the sorted alignmentIndex.
            for(size_t i=0; i<v.size(); i++) {
                alignmentTableSection[i] = v[i].second;
            }
        }
    }

    alignmentTable.unreserve();

}



void Assembler::accessAlignmentData()
{
    alignmentData.accessExistingReadOnly(largeDataName("AlignmentData"));
    alignmentTable.accessExistingReadOnly(largeDataName("AlignmentTable"));
}
void Assembler::accessAlignmentDataReadWrite()
{
    alignmentData.accessExistingReadWrite(largeDataName("AlignmentData"));
    alignmentTable.accessExistingReadWrite(largeDataName("AlignmentTable"));
}



void Assembler::checkAlignmentDataAreOpen() const
{
    if(!alignmentData.isOpen || !alignmentTable.isOpen()) {
        throw runtime_error("Alignment data are not accessible.");
    }
}





// Find in the alignment table the alignments involving
// a given oriented read, and return them with the correct
// orientation (this may involve a swap and/or reverse complement
// of the AlignmentInfo stored in the alignmentTable).
vector< pair<OrientedReadId, dinara::AlignmentInfo> >
    Assembler::findOrientedAlignments(
        OrientedReadId orientedReadId0Argument,
        bool inReadGraphOnly,
        vector<uint32_t>* alignmentIds) const
{
    const ReadId readId0 = orientedReadId0Argument.getReadId();
    const ReadId strand0 = orientedReadId0Argument.getStrand();

    vector< pair<OrientedReadId, AlignmentInfo> > result;

    // Loop over alignment involving this read, as stored in the
    // alignment table.
    const auto alignmentTable0 = alignmentTable[orientedReadId0Argument.getValue()];
    for(const auto i: alignmentTable0) {
        const AlignmentData& ad = alignmentData[i];

        // Get the oriented read ids that the AlignmentData refers to.
        OrientedReadId orientedReadId0(ad.readIds[0], 0);
        OrientedReadId orientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
        AlignmentInfo alignmentInfo = ad.info;

        // Skip it if it is not in the read graph and only alignments
        // in the read graph were requested.
        if(inReadGraphOnly and (not alignmentInfo.isInReadGraph)) {
            continue;
        }

        // Swap oriented reads, if necessary.
        if(orientedReadId0.getReadId() != readId0) {
            swap(orientedReadId0, orientedReadId1);
            alignmentInfo.swap();
        }
        DINARA_ASSERT(orientedReadId0.getReadId() == readId0);

        // Reverse complement, if necessary.
        if(orientedReadId0.getStrand() != strand0) {
            orientedReadId0.flipStrand();
            orientedReadId1.flipStrand();
            alignmentInfo.reverseComplement();
        }
        DINARA_ASSERT(orientedReadId0.getStrand() == strand0);
        DINARA_ASSERT(orientedReadId0 == orientedReadId0Argument);

        result.push_back(make_pair(orientedReadId1, alignmentInfo));
        if(alignmentIds) {
            alignmentIds->push_back(i);
        }
    }
    return result;
}



void Assembler::getAlignmentIdsSortedByInformativeSites(
    OrientedReadId orientedReadId0,
    vector<uint32_t>& alignmentIds,
    bool inReadGraphOnly) const
{
    checkAlignmentDataAreOpen();
    DINARA_ASSERT(alignmentTable.isOpen());

    const ReadId readId0 = orientedReadId0.getReadId();
    const span<const uint32_t> table = alignmentTable[orientedReadId0.getValue()];

    alignmentIds.clear();
    alignmentIds.reserve(table.size());
    for(const uint32_t alignmentId : table) {
        const AlignmentData& ad = alignmentData[alignmentId];
        if(inReadGraphOnly && (ad.info.isInReadGraph == 0)) {
            continue;
        }
        alignmentIds.push_back(alignmentId);
    }

    auto overlapLenOnRead = [&](const AlignmentData& ad) -> uint64_t {
        if(ad.readIds[0] == readId0) {
            return uint64_t(ad.qe) - uint64_t(ad.qs);
        } else {
            return uint64_t(ad.te) - uint64_t(ad.ts);
        }
    };

    sort(
        alignmentIds.begin(),
        alignmentIds.end(),
        [&](uint32_t aId, uint32_t bId) {
            const AlignmentData& a = alignmentData[aId];
            const AlignmentData& b = alignmentData[bId];
            const uint32_t aCount = a.getInformativeHetSiteCountFromReadPerspective(readId0);
            const uint32_t bCount = b.getInformativeHetSiteCountFromReadPerspective(readId0);
            if(aCount != bCount) {
                return aCount > bCount;
            }
            const uint64_t aLen = overlapLenOnRead(a);
            const uint64_t bLen = overlapLenOnRead(b);
            if(aLen != bLen) {
                return aLen > bLen;
            }
            return aId < bId;
        });
}



void Assembler::getCisAlignmentIdsSortedByInformativeSites(
    OrientedReadId orientedReadId0,
    vector<uint32_t>& alignmentIds,
    bool keptByBothSidesOnly) const
{
    checkAlignmentDataAreOpen();
    DINARA_ASSERT(alignmentTable.isOpen());

    const ReadId readId0 = orientedReadId0.getReadId();
    const span<const uint32_t> table = alignmentTable[orientedReadId0.getValue()];

    alignmentIds.clear();
    alignmentIds.reserve(table.size());
    for(const uint32_t alignmentId : table) {
        const AlignmentData& ad = alignmentData[alignmentId];
        if(!ad.isCisByBothSides()) {
            continue;
        }
        if(keptByBothSidesOnly && !ad.keptByBothSides()) {
            continue;
        }
        alignmentIds.push_back(alignmentId);
    }

    auto overlapLenOnRead = [&](const AlignmentData& ad) -> uint64_t {
        if(ad.readIds[0] == readId0) {
            return uint64_t(ad.qe) - uint64_t(ad.qs);
        } else {
            return uint64_t(ad.te) - uint64_t(ad.ts);
        }
    };

    sort(
        alignmentIds.begin(),
        alignmentIds.end(),
        [&](uint32_t aId, uint32_t bId) {
            const AlignmentData& a = alignmentData[aId];
            const AlignmentData& b = alignmentData[bId];
            if(a.informativeHetSiteScore != b.informativeHetSiteScore) {
                return a.informativeHetSiteScore > b.informativeHetSiteScore;
            }
            const uint64_t aLen = overlapLenOnRead(a);
            const uint64_t bLen = overlapLenOnRead(b);
            if(aLen != bLen) {
                return aLen > bLen;
            }
            return aId < bId;
        });
}

void Assembler::getAllCisAlignmentIdsSortedByInformativeSites(
    vector<uint32_t>& alignmentIds,
    bool keptByBothSidesOnly) const
{
    checkAlignmentDataAreOpen();

    const uint64_t alignmentCount = alignmentData.size();

    alignmentIds.clear();
    alignmentIds.reserve(alignmentCount);
    for(uint32_t alignmentId = 0; alignmentId < alignmentCount; ++alignmentId) {
        const AlignmentData& ad = alignmentData[alignmentId];
        if(!ad.isCisByBothSides()) {
            continue;
        }
        if(keptByBothSidesOnly && !ad.keptByBothSides()) {
            continue;
        }
        alignmentIds.push_back(alignmentId);
    }

    sort(
        alignmentIds.begin(),
        alignmentIds.end(),
        [&](uint32_t aId, uint32_t bId) {
            const AlignmentData& a = alignmentData[aId];
            const AlignmentData& b = alignmentData[bId];
            if(a.informativeHetSiteScore != b.informativeHetSiteScore) {
                return a.informativeHetSiteScore > b.informativeHetSiteScore;
            }
            return aId < bId;
        });
}






void Assembler::analyzeAlignmentMatrix(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1)
{
    // Get the oriented reads.
    const OrientedReadId orientedReadId0(readId0, strand0);
    const OrientedReadId orientedReadId1(readId1, strand1);

    // Get the markers sorted by kmerId.
    vector<MarkerWithOrdinal> markers0;
    vector<MarkerWithOrdinal> markers1;
    getMarkersSortedByKmerId(orientedReadId0, markers0);
    getMarkersSortedByKmerId(orientedReadId1, markers1);

    // Some iterators we will need.
    using MarkerIterator = vector<MarkerWithOrdinal>::const_iterator;
    const MarkerIterator begin0 = markers0.begin();
    const MarkerIterator end0   = markers0.end();
    const MarkerIterator begin1 = markers1.begin();
    const MarkerIterator end1   = markers1.end();

    // The number of markers in each oriented read.
    const int64_t n0 = end0 - begin0;
    const int64_t n1 = end1 - begin1;

    // We will use coordinates
    // x = ordinal0 + ordinal1
    // y = ordinal0 - ordinal1 (offset)
    // In these coordinates, diagonals in the alignment matrix
    // are lines of constant y and so they become horizontal.
    const int64_t xMin = 0;
    const int64_t xMax = n0 + n1 - 2;
    const int64_t yMin = -n1;
    const int64_t yMax = n0 - 1;
    const int64_t nx= xMax -xMin + 1;
    const int64_t ny= yMax -yMin + 1;

    // Create a histogram in cells of size (dx, dy).
    const int64_t dx = 100;
    const int64_t dy = 20;
    const int64_t nxCells = (nx-1)/dx + 1;
    const int64_t nyCells = (ny-1)/dy + 1;
    vector< vector<uint64_t> > histogram(nxCells, vector<uint64_t>(nyCells, 0));
    cout << "nxCells " << nxCells << endl;
    cout << "nyCells " << nyCells << endl;



    // Joint loop over the markers, looking for common k-mer ids.
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->kmerId < it1->kmerId) {
            ++it0;
        } else if(it1->kmerId < it0->kmerId) {
            ++it1;
        } else {

            // We found a common k-mer id.
            const KmerId kmerId = it0->kmerId;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in kmers0 and kmers1.
            MarkerIterator it0Begin = it0;
            MarkerIterator it1Begin = it1;
            MarkerIterator it0End = it0Begin;
            MarkerIterator it1End = it1Begin;
            while(it0End!=end0 && it0End->kmerId==kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->kmerId==kmerId) {
                ++it1End;
            }


            // Loop over pairs in the streaks.
            for(MarkerIterator jt0=it0Begin; jt0!=it0End; ++jt0) {
                const int64_t ordinal0 = jt0->ordinal;
                for(MarkerIterator jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const int64_t ordinal1 = int64_t(jt1->ordinal);

                    const int64_t x = ordinal0 + ordinal1;
                    const int64_t y = ordinal0 - ordinal1;
                    const int64_t ix = (x-xMin) / dx;
                    const int64_t iy = (y-yMin) / dy;
                    DINARA_ASSERT(ix >= 0);
                    DINARA_ASSERT(iy >= 0);
                    DINARA_ASSERT(ix < nxCells);
                    DINARA_ASSERT(iy < nyCells);

                    ++histogram[ix][iy];
                }
            }


            // Continue joint loop over k-mers.
            it0 = it0End;
            it1 = it1End;
        }
    }

    ofstream csv("Histogram.csv");
    PngImage image = PngImage(int(nxCells), int(nyCells));
    uint64_t activeCellCount = 0;
    for(int64_t iy=0; iy<nyCells; iy++) {
        for(int64_t ix=0; ix<nxCells; ix++) {
            const int64_t frequency = histogram[ix][iy];
            if(frequency > 0) {
                csv << frequency;
            }
            if(frequency >= 10) {
                ++activeCellCount;
            }
            // const int r = (frequency <= 10) ? 0 : min(int(255), int(10*frequency));
            const int r = (frequency>=10) ? 255 : 0;
            // const int r = min(int(255), int(10*frequency));
            const int g = r;
            const int b = r;
            image.setPixel(int(ix), int(iy), r, g, b);
            csv << ",";
        }
        csv << "\n";
    }
    image.write("Histogram.png");
    cout << activeCellCount << " active cells out of " << nxCells*nyCells << endl;


}


// Count the common marker near a given ordinal offset for
// two oriented reads. This can be used to check
// whether an alignmnent near the specified ordinal offset exists.
uint32_t Assembler::countCommonMarkersNearOffset(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int32_t offset,
    int32_t offsetTolerance
)
{
    const int32_t minOffset = offset - offsetTolerance;
    const int32_t maxOffset = offset + offsetTolerance;
    return countCommonMarkersWithOffsetIn(
        orientedReadId0, orientedReadId1,
        minOffset, maxOffset);
}
uint32_t Assembler::countCommonMarkersWithOffsetIn(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    int32_t minOffset,
    int32_t maxOffset
)
{
#if 1
    const bool debug = false;
#else
    const bool debug =
        orientedReadId0.getReadId() == 5 and
        orientedReadId0.getStrand() == 0 and
        orientedReadId1.getReadId() == 32 and
        orientedReadId1.getStrand() == 1;
#endif
    if(debug) {
        cout << "countCommonMarkersWithOffsetIn" << endl;
    }

    // Get the markers sorted by kmerId.
    checkMarkersAreOpen();
    vector<MarkerWithOrdinal> markers0;
    vector<MarkerWithOrdinal> markers1;
    getMarkersSortedByKmerId(orientedReadId0, markers0);
    getMarkersSortedByKmerId(orientedReadId1, markers1);

    // Some iterators we will need.
    using MarkerIterator = vector<MarkerWithOrdinal>::const_iterator;
    const MarkerIterator begin0 = markers0.begin();
    const MarkerIterator end0   = markers0.end();
    const MarkerIterator begin1 = markers1.begin();
    const MarkerIterator end1   = markers1.end();



    // Main loop, looking for common markers.
    uint32_t count = 0;
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->kmerId < it1->kmerId) {
            ++it0;
        } else if(it1->kmerId < it0->kmerId) {
            ++it1;
        } else {

            // We found a common marker.
            const KmerId kmerId = it0->kmerId;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in markers0 and markers1.
            MarkerIterator it0Begin = it0;
            MarkerIterator it1Begin = it1;
            MarkerIterator it0End = it0Begin;
            MarkerIterator it1End = it1Begin;
            while(it0End!=end0 && it0End->kmerId==kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->kmerId==kmerId) {
                ++it1End;
            }

            // Loop over pairs of markers in the streaks.
            for(MarkerIterator jt0=it0Begin; jt0!=it0End; ++jt0) {
                for(MarkerIterator jt1=it1Begin; jt1!=it1End; ++jt1) {
                    const int32_t ordinalOffset = int32_t(jt0->ordinal) - int32_t(jt1->ordinal);
                    if(debug) {
                        cout << jt0->ordinal << " " <<
                            jt1->ordinal << " " <<
                            ordinalOffset<< endl;
                    }
                    if(ordinalOffset >= minOffset and ordinalOffset <= maxOffset) {
                        ++count;
                    }

                }
            }
            // Continue main loop.
            it0 = it0End;
            it1 = it1End;
        }

    }

    return count;
}



// Check if an alignment between two reads should be suppressed,
// bases on the setting of command line option
// --Align.sameChannelReadAlignment.suppressDeltaThreshold.
bool Assembler::suppressAlignment(
    ReadId readId0,
    ReadId readId1,
    uint64_t delta)
{

    // If the ch meta data fields of the two reads are missing or different,
    // don't suppress the alignment.
    // Check the channel first for efficiency,
    // so we can return faster in most cases.
    const auto ch0 = reads->getMetaData(readId0, "ch");
    if(ch0.empty()) {
        return false;
    }
    const auto ch1 = reads->getMetaData(readId1, "ch");
    if(ch1.empty()) {
        return false;
    }
    if(ch0 not_eq ch1) {
        return false;
    }



    // If the sampleid meta data fields of the two reads are missing or different,
    // don't suppress the alignment.
    const auto sampleid0 = reads->getMetaData(readId0, "sampleid");
    if(sampleid0.empty()) {
        return false;
    }
    const auto sampleid1 = reads->getMetaData(readId1, "sampleid");
    if(sampleid1.empty()) {
        return false;
    }
    if(sampleid0 not_eq sampleid1) {
        return false;
    }



    // If the runid meta data fields of the two reads are missing or different,
    // don't suppress the alignment.
    const auto runid0 = reads->getMetaData(readId0, "runid");
    if(runid0.empty()) {
        return false;
    }
    const auto runid1 = reads->getMetaData(readId1, "runid");
    if(runid1.empty()) {
        return false;
    }
    if(runid0 not_eq runid1) {
        return false;
    }

    // cout << "Checking " << readId0 << " " << readId1 << endl;


    // If the read meta data fields of the two reads are missing,
    // don't suppress the alignment.
    const auto read0 = reads->getMetaData(readId0, "read");
    if(read0.empty()) {
        return false;
    }
    const auto read1 = reads->getMetaData(readId1, "read");
    if(read1.empty()) {
        return false;
    }
    // cout << read0 << " " << read1 << endl;




    // Convert the read meta data fields to integers.
    // Keep in mind the span<char> is not null-terminated.
    const int64_t r0 = int64_t(atoul(read0));
    const int64_t r1 = int64_t(atoul(read1));
    // cout << r0 << " " << r1 << endl;



    // Suppress the alignment if the absolute difference of the
    // read meta data fields is less than delta.
    return std::abs(r0 - r1) < int64_t(delta);

}



// Remove all alignment candidates for which suppressAlignment
// returns false.
void Assembler::suppressAlignmentCandidates(
    uint64_t delta,
    uint64_t threadCount)
{
    performanceLog << timestamp << "Suppressing alignment candidates." << endl;

    // Allocate memory for flags to keep track of which alignments
    // should be suppressed.
    suppressAlignmentCandidatesData.suppress.createNew(
        largeDataName("tmp-suppressAlignmentCandidates"), largeDataPageSize);
    const uint64_t candidateCount = alignmentCandidates.candidates.size();
    suppressAlignmentCandidatesData.suppress.resize(candidateCount);

    // Figure out which candidates should be suppressed.
    suppressAlignmentCandidatesData.delta = delta;
    const uint64_t batchSize = 10000;
    setupLoadBalancing(alignmentCandidates.candidates.size(), batchSize);
    runThreads(&Assembler::suppressAlignmentCandidatesThreadFunction, threadCount);

    ofstream csv("SuppressedAlignmentCandidates.csv");
    csv << "ReadId0,ReadId1,SameStrand,Name0,Name1,MetaData0,MetaData1" << endl;

    // Suppress the alignment candidates we flagged.
    cout << "Number of alignment candidates before suppression is " << candidateCount << endl;
    uint64_t j = 0;
    uint64_t suppressCount = 0;
    for(uint64_t i=0; i<candidateCount; i++) {
        if(suppressAlignmentCandidatesData.suppress[i]) {
            ++suppressCount;
            const ReadId readId0 = alignmentCandidates.candidates[i].readIds[0];
            const ReadId readId1 = alignmentCandidates.candidates[i].readIds[1];
            csv << readId0 << "," << readId1 << ","
                << (alignmentCandidates.candidates[i].isSameStrand ? "Yes" : "No") << ","
                << reads->getReadName(readId0) << "," << reads->getReadName(readId1) << ","
                << reads->getReadMetaData(readId0) << "," << reads->getReadMetaData(readId1) << endl;
        } else {
            alignmentCandidates.candidates[j++] =
                alignmentCandidates.candidates[i];
        }
    }
    DINARA_ASSERT(j + suppressCount == candidateCount);
    alignmentCandidates.candidates.resize(j);
    cout << "Suppressed " << suppressCount << " alignment candidates." << endl;
    cout << "Number of alignment candidates after suppression is " << j << endl;


    // Clean up.
    suppressAlignmentCandidatesData.suppress.remove();

    performanceLog << timestamp << "Done suppressing alignment candidates." << endl;
}



void Assembler::suppressAlignmentCandidatesThreadFunction(uint64_t)
{
    const uint64_t delta = suppressAlignmentCandidatesData.delta;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over candidate alignments in this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const OrientedReadPair& p = alignmentCandidates.candidates[i];
            // cout << "Checking " << p.readIds[0] << " " << p.readIds[1] <<  " " << int(p.isSameStrand) << endl;
            suppressAlignmentCandidatesData.suppress[i] =
                suppressAlignment(p.readIds[0], p.readIds[1], delta);
        }

    }
}
