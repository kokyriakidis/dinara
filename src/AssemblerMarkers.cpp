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

// hifiasm minimizer bridge (no-HPC position source, benchmarking alternative
// to simd-minimizers). Positions are still resolved to canonical KmerIds here.
#include "hifiasm_overlaps.h"

// myloasm marker bridge: open syncmers + SNPmers as a marker POSITION source.
// SNPmers are detected from the global k-mer spectrum, so this is a single
// whole-run call. Only positions are used; KmerIds are encoded by dinara.
#include "myloasm_ffi.h"

// Standard library.
#include "fstream.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
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





// ============================================================================
// SIMD MINIMIZERS IMPLEMENTATION
// ============================================================================

// Helper: compute deduplicated canonical minimizer positions and their KmerIds for a read.
// If validMarkers is provided, it is filled with (position, kmerId).
// If validMarkers is nullptr, it returns the count of markers without filling a buffer.
static size_t getMinimizerMarkersForRead(
    ReadId readId,
    const Reads& reads,
    int k,          // marker-encoding k-mer length (what dinara stores/indexes)
    int w,          // window for the simd backend
    int sketchK,    // hifiasm position-selection k (its tuned 51)
    int sketchW,    // hifiasm position-selection w (its tuned 51)
    bool useHifiasm,
    SimdSketcher* sketcher,
    hifiasm_sketch_ctx_t* hifiasmCtx,
    const hifiasm_filter_t* hifiasmFilter,
    int hifiasmSampleDist,
    const shared_ptr<KmerChecker>& kmerChecker,
    string& readSequence,
    std::vector<uint32_t>& positionBuffer,
    std::vector<std::pair<uint32_t, KmerId>>* validMarkers = nullptr
) {
    if(validMarkers) validMarkers->clear();
    const LongBaseSequenceView read = reads.getRead(readId);
    const uint64_t baseCount = read.baseCount;

    // A read must be long enough for BOTH the encoding k-mer and, on the
    // hifiasm path, the (longer) selection k-mer. Using the max avoids sketching
    // a read that can't host a full selection k-mer.
    const int minLen = useHifiasm ? std::max(k, sketchK) : k;
    if(baseCount < uint64_t(minLen)) {
        return 0;
    }

    // Convert read to string. Both minimizer backends read raw ASCII bases
    // (no HPC); positions are k-mer START offsets on the forward strand.
    readSequence.resize(baseCount);
    for(uint64_t i = 0; i < baseCount; i++) {
        readSequence[i] = read[i].character();
    }

    if(useHifiasm) {
        // hifiasm sketcher (no-HPC) via a reusable per-thread context, so
        // steady-state sketching does no per-read heap allocation. The context
        // returns START positions already sorted and deduplicated, so the
        // std::sort/std::unique below is skipped on this path.
        //
        // Position selection runs at hifiasm's tuned sketchK/sketchW (=51), the
        // same k/w it uses to seed overlaps, so marker seeds match overlap
        // seeds. dinara then encodes its own k-mer (length k, e.g. 50) at each
        // returned position in the shared tail below: a 51-mer at offset p has
        // dinara's 50-mer as its prefix, so the position is identical and no
        // realignment is needed.
        const hifiasm_minimizer_t* mz = nullptr;
        int n = 0;
        // When a filter is supplied, use the overlap-parity variant: it applies
        // hifiasm's high-occurrence k-mer filter (rarest-k-mer window selection)
        // and distance subsampling, matching the seeds hifiasm uses for overlap
        // detection. Otherwise fall back to the plain unfiltered sketch.
        const int rc = hifiasmFilter ?
            hifiasm_sketch_minimizers_ctx_filtered(
                hifiasmCtx, readSequence.c_str(), int(baseCount), sketchW, sketchK,
                /*is_hpc*/ 0, hifiasmFilter, hifiasmSampleDist, &mz, &n) :
            hifiasm_sketch_minimizers_ctx(
                hifiasmCtx, readSequence.c_str(), int(baseCount), sketchW, sketchK,
                /*is_hpc*/ 0, &mz, &n);
        if(rc != 0) {
            // Treat a sketch failure as "no markers" rather than aborting the
            // whole run; the read simply contributes nothing.
            positionBuffer.clear();
        } else {
            positionBuffer.resize(size_t(n));
            for(int i = 0; i < n; i++) {
                positionBuffer[size_t(i)] = mz[i].pos;
            }
        }
    } else {
        // simd-minimizers library. Compute canonical minimizer positions
        // (START offsets), then sort + deduplicate to guarantee monotonic
        // ordering (the library does not guarantee sorted, unique output).
        MinimizerList minimizerList = canonical_minimizer_positions(
            sketcher,
            readSequence.c_str(),
            readSequence.size());
        positionBuffer.assign(minimizerList.data, minimizerList.data + minimizerList.len);
        free_minimizer_list(minimizerList);

        std::sort(positionBuffer.begin(), positionBuffer.end());
        positionBuffer.erase(
            std::unique(positionBuffer.begin(), positionBuffer.end()),
            positionBuffer.end());
    }

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

void Assembler::findMarkersSimdMinimizers(uint64_t threadCount, int k, int w,
    bool useHifiasm, const void* hifiasmFilter, int hifiasmSampleDist,
    int sketchK, int sketchW)
{
    reads->checkReadsAreOpen();

    // The hifiasm overlap-path filter only applies on the hifiasm sketch path.
    const bool useFilter = useHifiasm && (hifiasmFilter != nullptr);

    // Effective position-selection k/w on the hifiasm path (fall back to k/w).
    const int effSketchK = (useHifiasm && sketchK > 0) ? sketchK : k;
    const int effSketchW = (useHifiasm && sketchW > 0) ? sketchW : w;

    performanceLog << timestamp << "Finding markers using "
        << (useHifiasm ? "hifiasm (no-HPC)" : "SIMD") << " minimizers (k=" << k
        << ", w=" << w << ")"
        << (useHifiasm && (effSketchK != k || effSketchW != w) ?
            " selecting positions at hifiasm k=w=" + std::to_string(effSketchK) : "")
        << (useFilter ? " with overlap-path filter (sampleDist=" : "")
        << (useFilter ? std::to_string(hifiasmSampleDist) : "")
        << (useFilter ? ")" : "")
        << " in " << reads->readCount() << " reads." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // Store parameters. assemblerInfo->k is the marker-encoding k (what the rest
    // of the pipeline indexes with); the selection k/w only affect which
    // positions the hifiasm sketch picks.
    assemblerInfo->k = k;
    findMarkersSimdMinimizersData.k = k;
    findMarkersSimdMinimizersData.w = w;
    findMarkersSimdMinimizersData.sketchK = effSketchK;
    findMarkersSimdMinimizersData.sketchW = effSketchW;
    findMarkersSimdMinimizersData.useHifiasm = useHifiasm;
    findMarkersSimdMinimizersData.hifiasmFilter = useFilter ? hifiasmFilter : nullptr;
    findMarkersSimdMinimizersData.hifiasmSampleDist = useFilter ? hifiasmSampleDist : 0;

    // Create the markers and markerKmerIds data structures.
    markers->createNew(largeDataName("Markers"), largeDataPageSize);
    markerKmerIds->createNew(largeDataName("MarkerKmerIds"), largeDataPageSize);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t batchSize = 100;

    // Sketching is the dominant cost here, so each read is sketched exactly once.
    // Pass 1 sketches and stages the forward-strand markers per read (in
    // stagedMarkers) while incrementing the CSR counts; pass 2 just fills the
    // memory-mapped vectors from the staged data (deriving strand 1) and frees
    // each read's stage as it goes.
    //
    // Memory tradeoff: staging holds the forward-strand markers for all reads in
    // RAM between the two passes (~half of what the final memory-mapped
    // markers/markerKmerIds hold, since those store both strands). This is the
    // cost of not sketching twice; the peak is transient and freed during pass 2.
    findMarkersSimdMinimizersData.stagedMarkers.clear();
    findMarkersSimdMinimizersData.stagedMarkers.resize(readCount);

    // Pass 1: sketch once, stage forward markers, count.
    markers->beginPass1(2 * readCount);
    markerKmerIds->beginPass1(2 * readCount);
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::findMarkersSimdMinimizersPass1, threadCount);

    // Pass 2: store markers from the staged data (no re-sketch).
    markers->beginPass2();
    markerKmerIds->beginPass2();
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::findMarkersSimdMinimizersPass2, threadCount);

    markers->endPass2(false);
    markerKmerIds->endPass2(false);

    // Release staging memory.
    findMarkersSimdMinimizersData.stagedMarkers.clear();
    findMarkersSimdMinimizersData.stagedMarkers.shrink_to_fit();

    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    const char* backend = findMarkersSimdMinimizersData.useHifiasm ?
        "hifiasm (no-HPC)" : "SIMD";
    performanceLog << timestamp << "Finding markers using " << backend
        << " minimizers completed in " << tTotal << " s." << endl;
    cout << "Created " << markers->totalSize() << " markers using " << backend
        << " minimizers." << endl;
}

void Assembler::findMarkersSimdMinimizersPass1(size_t /* threadId */)
{
    const int k = findMarkersSimdMinimizersData.k;
    const int w = findMarkersSimdMinimizersData.w;
    const bool useHifiasm = findMarkersSimdMinimizersData.useHifiasm;
    // Position-selection k/w for the hifiasm sketch (fall back to k/w).
    const int sketchK = findMarkersSimdMinimizersData.sketchK > 0 ?
        findMarkersSimdMinimizersData.sketchK : k;
    const int sketchW = findMarkersSimdMinimizersData.sketchW > 0 ?
        findMarkersSimdMinimizersData.sketchW : w;
    const hifiasm_filter_t* hifiasmFilter =
        static_cast<const hifiasm_filter_t*>(findMarkersSimdMinimizersData.hifiasmFilter);
    const int hifiasmSampleDist = findMarkersSimdMinimizersData.hifiasmSampleDist;
    // Allocate only the backend actually used on this path: the hifiasm sketch
    // context or the simd sketcher (each is per-thread; neither is thread-safe).
    SimdSketcher* sketcher = useHifiasm ? nullptr : simd_sketcher_new(
        static_cast<uint8_t>(k), static_cast<uint8_t>(w));
    hifiasm_sketch_ctx_t* hifiasmCtx = useHifiasm ? hifiasm_sketch_ctx_init() : nullptr;
    string readSequence;
    std::vector<uint32_t> positionBuffer;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            // Sketch once and stage the forward-strand markers for this read.
            // stagedMarkers is indexed by ReadId; each thread owns a disjoint
            // set of read IDs via the load balancer, so no locking is needed.
            std::vector<std::pair<uint32_t, KmerId>>& staged =
                findMarkersSimdMinimizersData.stagedMarkers[readId];
            getMinimizerMarkersForRead(
                readId, *reads, k, w, sketchK, sketchW, useHifiasm, sketcher,
                hifiasmCtx, hifiasmFilter, hifiasmSampleDist,
                kmerChecker, readSequence, positionBuffer, &staged);

            const uint64_t count64 = staged.size();
            this->markers->incrementCount(OrientedReadId(readId, 0).getValue(), count64);
            this->markers->incrementCount(OrientedReadId(readId, 1).getValue(), count64);
            markerKmerIds->incrementCount(OrientedReadId(readId, 0).getValue(), count64);
            markerKmerIds->incrementCount(OrientedReadId(readId, 1).getValue(), count64);
        }
    }
    if(sketcher) simd_sketcher_free(sketcher);
    if(hifiasmCtx) hifiasm_sketch_ctx_destroy(hifiasmCtx);
}

void Assembler::findMarkersSimdMinimizersPass2(size_t /* threadId */)
{
    // No sketching here: consume the forward-strand markers staged in pass 1.
    const int k = findMarkersSimdMinimizersData.k;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(ReadId readId = ReadId(begin); readId != ReadId(end); ++readId) {
            std::vector<std::pair<uint32_t, KmerId>>& staged =
                findMarkersSimdMinimizersData.stagedMarkers[readId];

            if(!staged.empty()) {
                const LongBaseSequenceView read = reads->getRead(readId);

                CompressedMarker* markerPointerStrand0 = this->markers->begin(OrientedReadId(readId, 0).getValue());
                CompressedMarker* markerPointerStrand1 = this->markers->end(OrientedReadId(readId, 1).getValue()) - 1;

                KmerId* kmerIdPointerStrand0 = markerKmerIds->begin(OrientedReadId(readId, 0).getValue());
                KmerId* kmerIdPointerStrand1 = markerKmerIds->end(OrientedReadId(readId, 1).getValue()) - 1;

                for(const auto& [position, kmerId] : staged) {
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

            // Free this read's staged markers now that they are stored.
            std::vector<std::pair<uint32_t, KmerId>>().swap(staged);
        }
    }
}


// ============================================================================
// MYLOASM MARKER IMPLEMENTATION (open syncmers + SNPmers as position source)
// ============================================================================
//
// Design: hifiasm still finds overlap pairs + intervals; this only changes what
// fills the marker table. myloasm supplies marker POSITIONS (open syncmers +
// SNPmers, the latter detected from the global k-mer spectrum across all reads),
// and dinara encodes its OWN canonical KmerId (length assemblerInfo->k) at each
// position -- exactly as the simd/hifiasm marker paths do. So markers,
// markerKmerIds, MarkerKmers and every downstream consumer keep the same shape
// and the same "KmerId == canonical k-mer re-extracted from the read" invariant.
// deriveChainFromInterval then chains syncmer+SNPmer anchors inside hifiasm's
// intervals with no change, and anchors/phasing inherit SNPmer positions.
//
// Because SNPmer detection is global, myloasm_index_reads is a single whole-run
// call (not per-read). We map each myloasm read back to a dinara ReadId by name
// (dinara stores the first header token; myloasm returns the same base id),
// merge that read's syncmer + SNPmer positions, sort + deduplicate by position,
// drop positions where pos + k would run off the read, encode dinara's KmerId at
// each surviving position, and stage the forward-strand (position, KmerId)
// pairs. The existing two-pass CSR store (findMarkersSimdMinimizersPass2)
// derives strand 1 by reverse-complement and fills the memory-mapped vectors.
void Assembler::findMarkersMyloasm(
    uint64_t threadCount,
    const vector<string>& inputFileNames)
{
    reads->checkReadsAreOpen();

    // myloasm places markers with a k=21 open-syncmer / SNPmer scheme and matches
    // them internally at k=20 (its 21-mer key clipped to 20 by masking the low
    // 40 bits, which keeps the middle base so a SNPmer's two alleles stay
    // distinct). We mirror that: dinara encodes its OWN canonical 20-mer at each
    // myloasm position. Encode k is FIXED at 20 here, independent of Kmers.k --
    // it is myloasm's internal marker k, not dinara's global k-mer length. This
    // is deliberately decoupled from:
    //   * hifiasm overlap detection, which runs at k=51 in HPC space on a
    //     separate path (hifiasm_detect_overlaps_mem) and never reads this k; and
    //   * Kmers.k, which stays at its own value for any non-myloasm use.
    // The het base of a SNPmer (myloasm index 10, i.e. the 11th base of the
    // 21-mer) lies inside [pos, pos+20), so the two alleles remain distinct in
    // dinara's encoding too. dinara's MarkerKmers and the anchor builders
    // re-extract and re-canonicalize k-mers from the read sequence, so using
    // dinara's own 20-mer (rather than storing myloasm's masked key) keeps
    // chaining and anchors self-consistent.
    constexpr int k = 20;

    if(inputFileNames.empty()) {
        throw runtime_error(
            "Kmers.useMyloasmMarkers requires at least one --input read file "
            "(myloasm re-reads the inputs to detect SNPmers).");
    }

    performanceLog << timestamp
        << "Finding markers using myloasm open syncmers + SNPmers (encode k="
        << k << ", myloasm marker k=21 clipped to 20) in "
        << reads->readCount() << " reads." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    // Marker encoding k is fixed at 20 on this path (see above); record it so
    // downstream consumers (MarkerKmers, anchors) encode/re-extract at the same
    // length. This overrides Kmers.k for the marker/anchor index only.
    assemblerInfo->k = k;

    // Create the markers and markerKmerIds data structures.
    markers->createNew(largeDataName("Markers"), largeDataPageSize);
    markerKmerIds->createNew(largeDataName("MarkerKmerIds"), largeDataPageSize);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();

    // --- Single whole-run FFI call: detect SNPmers + index every read. --------
    vector<const char*> readFiles;
    readFiles.reserve(inputFileNames.size());
    for(const string& f: inputFileNames) {
        readFiles.push_back(f.c_str());
    }

    MyloReadIndex index;
    performanceLog << timestamp
        << "Indexing reads with myloasm (syncmers + SNPmers) from "
        << inputFileNames.size() << " input file(s)." << endl;
    const int rc = myloasm_index_reads(
        readFiles.data(), readFiles.size(),
        // myloasm indexes at its native marker k=21 (its own default scheme);
        // dinara then encodes a 20-mer at each returned position (k above).
        /*kmer_size*/ 21, /*c*/ 0, int(threadCount), &index);
    if(rc != 0) {
        myloasm_read_index_free(&index);
        throw runtime_error(
            "myloasm_index_reads failed (code " + std::to_string(rc) + ").");
    }
    performanceLog << timestamp
        << "myloasm indexed " << index.n_reads << " reads (marker k="
        << index.k << ")." << endl;

    // --- Stage forward-strand (position, KmerId) per dinara ReadId. -----------
    // Encoding runs here (myloasm returns one big index; the per-position k-mer
    // extraction is cheap relative to the FFI's global counting).
    findMarkersSimdMinimizersData.k = k;
    findMarkersSimdMinimizersData.stagedMarkers.clear();
    findMarkersSimdMinimizersData.stagedMarkers.resize(readCount);

    uint64_t unmatchedReads = 0;
    uint64_t totalStaged = 0;
    std::vector<uint32_t> positions;
    for(size_t r = 0; r < index.n_reads; r++) {
        const MyloReadMarkers& mr = index.reads[r];

        // Map myloasm read (by name) to a dinara ReadId. The name is the first
        // header token, not NUL-terminated; use (ptr,len) directly.
        const span<const char> name(mr.name, mr.name + mr.name_len);
        const ReadId readId = reads->getReadId(name);
        if(readId == invalidReadId) {
            ++unmatchedReads;
            continue;
        }

        const LongBaseSequenceView read = reads->getRead(readId);
        const uint64_t baseCount = read.baseCount;

        // Merge syncmer + SNPmer positions (both are forward-strand starts).
        positions.clear();
        positions.reserve(mr.n_syncmers + mr.n_snpmers);
        for(size_t i = 0; i < mr.n_syncmers; i++) {
            positions.push_back(mr.syncmers[i].pos);
        }
        for(size_t i = 0; i < mr.n_snpmers; i++) {
            positions.push_back(mr.snpmers[i].pos);
        }
        if(positions.empty()) {
            continue;
        }

        // Sort + dedup by position. A locus that is both a syncmer and a SNPmer
        // yields one marker; the encoded KmerId is a pure function of position,
        // so deduping by position loses nothing.
        std::sort(positions.begin(), positions.end());
        positions.erase(
            std::unique(positions.begin(), positions.end()), positions.end());

        std::vector<std::pair<uint32_t, KmerId>>& staged =
            findMarkersSimdMinimizersData.stagedMarkers[readId];
        staged.reserve(positions.size());
        for(const uint32_t position: positions) {
            // Drop markers whose encoded k-mer would run past the read end.
            if(uint64_t(position) + uint64_t(k) > baseCount) {
                continue;
            }
            Kmer kmer;
            extractKmer(read, uint64_t(position), uint64_t(k), kmer);
            staged.emplace_back(position, kmer.id(uint64_t(k)));
        }
        totalStaged += staged.size();
    }

    myloasm_read_index_free(&index);

    if(unmatchedReads != 0) {
        performanceLog << timestamp
            << "WARNING: " << unmatchedReads
            << " myloasm reads did not match a dinara read by name and were "
            << "skipped." << endl;
    }
    performanceLog << timestamp
        << "Staged " << totalStaged << " myloasm markers (forward strand)."
        << endl;

    // --- Two-pass CSR store (reuse the simd path's pass 2). -------------------
    // Pass 1: count from the staged sizes (single pass, cheap; no re-sketch).
    markers->beginPass1(2 * readCount);
    markerKmerIds->beginPass1(2 * readCount);
    for(ReadId readId = 0; readId < ReadId(readCount); ++readId) {
        const uint64_t count64 =
            findMarkersSimdMinimizersData.stagedMarkers[readId].size();
        markers->incrementCount(OrientedReadId(readId, 0).getValue(), count64);
        markers->incrementCount(OrientedReadId(readId, 1).getValue(), count64);
        markerKmerIds->incrementCount(OrientedReadId(readId, 0).getValue(), count64);
        markerKmerIds->incrementCount(OrientedReadId(readId, 1).getValue(), count64);
    }

    // Pass 2: store markers + derive strand 1 (identical to the simd path).
    markers->beginPass2();
    markerKmerIds->beginPass2();
    const uint64_t batchSize = 100;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::findMarkersSimdMinimizersPass2, threadCount);
    markers->endPass2(false);
    markerKmerIds->endPass2(false);

    // Release staging memory.
    findMarkersSimdMinimizersData.stagedMarkers.clear();
    findMarkersSimdMinimizersData.stagedMarkers.shrink_to_fit();

    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double(
        (std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin)).count());
    performanceLog << timestamp
        << "Finding markers using myloasm syncmers + SNPmers completed in "
        << tTotal << " s." << endl;
    cout << "Created " << markers->totalSize()
        << " markers using myloasm syncmers + SNPmers." << endl;
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

    // optimization: if we already have them (from the SIMD minimizer path), don't recompute.
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
        ", het " << assemblerInfo->kmerDistributionInfo.coverageHet <<
        ", hom " << assemblerInfo->kmerDistributionInfo.coverageHom << endl;
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
        ", het " << assemblerInfo->kmerDistributionInfo.coverageHet <<
        ", hom " << assemblerInfo->kmerDistributionInfo.coverageHom << endl;
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
void Assembler::applyKmerCountFilter(
    uint64_t minFreq, uint64_t maxFreq, uint64_t threadCount, bool filterPalindromes,
    bool filterRepeatKmers, bool filterLowComplexity)
{
    performanceLog << timestamp << "Filtering markers by frequency [" 
                   << minFreq << ", " << maxFreq << "] and palindromes: " 
                   << (filterPalindromes ? "yes" : "no")
                   << ", repeat k-mers: " << (filterRepeatKmers ? "yes" : "no")
                   << ", low-complexity: " << (filterLowComplexity ? "yes" : "no")
                   << "." << endl;
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
    applyKmerCountFilterData.filterRepeatKmers = filterRepeatKmers;
    applyKmerCountFilterData.filterLowComplexity = filterLowComplexity;

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
    const bool filterRepeatKmers = applyKmerCountFilterData.filterRepeatKmers;
    const bool filterLowComplexity = applyKmerCountFilterData.filterLowComplexity;

    // Short-period tandem-repeat predicate (periods 1-6). Mirrors
    // filterMarkerGraphVerticesByRepeatKmers (--max-anchor-repeat-length).
    const vector<uint64_t> maxAnchorRepeatLength = {6, 4, 4, 4, 4};
    auto isRepeatKmer = [&](const Kmer& kmer0) -> bool {
        for(uint64_t i = 0; i < maxAnchorRepeatLength.size(); i++) {
            const uint64_t period = i + 1;
            const uint64_t maxAllowedCopyNumber = maxAnchorRepeatLength[i];
            uint64_t copies = 0;
            switch(period) {
            case 1: copies = kmer0.countExactRepeatCopies<1>(k); break;
            case 2: copies = kmer0.countExactRepeatCopies<2>(k); break;
            case 3: copies = kmer0.countExactRepeatCopies<3>(k); break;
            case 4: copies = kmer0.countExactRepeatCopies<4>(k); break;
            case 5: copies = kmer0.countExactRepeatCopies<5>(k); break;
            case 6: copies = kmer0.countExactRepeatCopies<6>(k); break;
            default: copies = 0; break;
            }
            if(copies > maxAllowedCopyNumber) return true;
        }
        return false;
    };

    // Low-complexity predicate by distinct sub-k-mer count (lengths 1-3).
    // Mirrors filterMarkerGraphVerticesByDistinctSubkmerCount
    // (--min-anchor-distinct-subkmer-count).
    const vector<uint64_t> minAnchorDistinctSubkmerCount = {4, 12, 24};
    auto isLowComplexity = [&](const Kmer& kmer) -> bool {
        for(uint64_t i = 0; i < minAnchorDistinctSubkmerCount.size(); i++) {
            const uint64_t subKmerLength = i + 1;
            const uint64_t minAllowedCount = minAnchorDistinctSubkmerCount[i];
            if(kmer.count(subKmerLength, k) < minAllowedCount) return true;
        }
        return false;
    };

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

                // Palindromic k-mers can cause ambiguity in directed graph construction
                // because their forward and reverse orientations are identical. They are
                // rejected before the frequency lookup because their frequency cannot
                // change the filtering decision.
                if(applyKmerCountFilterData.filterPalindromes && kmerId == rcKmerId) {
                    continue;
                }

                // Drop short-period tandem-repeat and low-complexity k-mers
                // (same predicates as the marker-graph vertex filters), so these
                // minimizers never seed a marker or a marker-graph vertex.
                // Checked before the frequency lookup: their k-mer content alone
                // decides rejection, independent of coverage.
                if(filterRepeatKmers && isRepeatKmer(kmer)) {
                    continue;
                }
                if(filterLowComplexity && isLowComplexity(kmer)) {
                    continue;
                }

                const KmerId canonical = std::min(kmerId, rcKmerId);
                const uint64_t freq = kmerCounter->getFrequencyFast(canonical);
                
                if(freq >= minF && freq <= maxF) {
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


// =============================================================================
// Remove markers whose canonical k-mer appears more than once in the reference.
// A k-mer that is non-unique in the reference cannot serve as an unambiguous
// anchor for SV detection. This removes such markers from ALL reads (and the
// reference itself) so they don't pollute the inverted index.
//
// The approach:
//   1. Scan reference reads (strand 0), collect all canonical k-mer IDs.
//   2. Sort, find duplicates, build a sorted blacklist of non-unique k-mers.
//   3. Two-pass filter (same pattern as applyKmerCountFilter):
//      - Pass 1: mark which markers to keep (canonical k-mer not in blacklist).
//      - Pass 2: rebuild markers and markerKmerIds without blacklisted entries.
// =============================================================================
void Assembler::removeNonUniqueReferenceMarkers(
    uint64_t referenceReadCount,
    uint64_t threadCount,
    uint64_t maxRefKmerFreq)
{
    const auto tBegin = std::chrono::steady_clock::now();
    performanceLog << timestamp
        << "Removing non-unique reference k-mers (maxFreq="
        << maxRefKmerFreq << ")." << endl;

    checkMarkersAreOpen();
    DINARA_ASSERT(markerKmerIds->isOpen());

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t k = assemblerInfo->k;
    const uint64_t readCount = reads->readCount();

    // =========================================================================
    // Phase 1: Collect all canonical k-mer IDs from reference reads (strand 0).
    // =========================================================================
    vector<KmerId> refKmers;
    for(uint64_t refId = 0; refId < referenceReadCount; ++refId) {
        const auto orientedReadId = OrientedReadId(ReadId(refId), 0);
        const auto kmerIds = (*markerKmerIds)[orientedReadId.getValue()];
        for(uint64_t i = 0; i < kmerIds.size(); ++i) {
            refKmers.push_back(kmerIds[i]);
        }
    }

    cout << "Reference k-mers (strand 0): " << refKmers.size() << endl;

    // =========================================================================
    // Phase 2: Sort and find k-mers exceeding the frequency threshold.
    // =========================================================================
    // A k-mer appearing more than maxRefKmerFreq times in the reference
    // is blacklisted. Default maxRefKmerFreq=1 removes all non-unique
    // k-mers (original behavior). Higher values retain low-copy repeat
    // markers, which is useful for SV detection in VNTR regions where
    // having some markers is better than none.
    sort(refKmers.begin(), refKmers.end());

    vector<KmerId> blacklist;
    size_t i = 0;
    while(i < refKmers.size()) {
        size_t j = i + 1;
        while(j < refKmers.size() && refKmers[j] == refKmers[i]) ++j;
        const uint64_t freq = j - i;
        if(freq > maxRefKmerFreq) {
            blacklist.push_back(refKmers[i]);
        }
        i = j;
    }

    cout << "Reference k-mers exceeding maxFreq=" << maxRefKmerFreq
         << " (blacklisted): " << blacklist.size() << endl;

    // =========================================================================
    // Phase 2b: Adaptive rescue for marker-depleted windows.
    // =========================================================================
    // After blacklisting, some reference windows may have zero remaining
    // markers (common in VNTRs). For those windows, rescue k-mers with
    // freq 2..rescueMaxFreq by removing them from the blacklist. This
    // retains some markers in depleted regions without affecting non-
    // depleted regions.
    if(!blacklist.empty() && maxRefKmerFreq == 1) {
        const uint64_t rescueMaxFreq = 2;
        const uint32_t windowSize = 50;

        // Compute reference marker positions after blacklisting.
        // Use reference read 0 (strand 0).
        const auto refOrientedReadId = OrientedReadId(ReadId(0), 0);
        const auto refMarkerPositions = (*markers)[refOrientedReadId.getValue()];
        const auto refKmerIdsSpan = (*markerKmerIds)[refOrientedReadId.getValue()];
        const uint64_t refLen = reads->getReadRawSequenceLength(ReadId(0));
        const uint32_t nWindows = uint32_t((refLen + windowSize - 1) / windowSize);

        // Count markers per window after blacklisting.
        vector<uint32_t> windowMarkerCount(nWindows, 0);
        for(size_t mi = 0; mi < refMarkerPositions.size(); ++mi) {
            const KmerId kid = refKmerIdsSpan[mi];
            if(std::binary_search(blacklist.begin(), blacklist.end(), kid))
                continue;
            const uint32_t pos = refMarkerPositions[mi].position;
            const uint32_t w = pos / windowSize;
            if(w < nWindows) ++windowMarkerCount[w];
        }

        // Find depleted windows (zero markers after blacklisting).
        uint32_t depletedWindows = 0;
        for(uint32_t w = 0; w < nWindows; ++w) {
            if(windowMarkerCount[w] == 0) ++depletedWindows;
        }

        if(depletedWindows > 0) {
            // Build a set of reference positions in depleted windows.
            // For each blacklisted k-mer, check if it has an occurrence
            // in a depleted window. If so, and its freq <= rescueMaxFreq,
            // rescue it.
            // Build sorted list of rescue-eligible k-mers
            // (freq 2..rescueMaxFreq).
            vector<KmerId> rescueEligible;
            i = 0;
            while(i < refKmers.size()) {
                size_t j2 = i + 1;
                while(j2 < refKmers.size() && refKmers[j2] == refKmers[i]) ++j2;
                const uint64_t freq = j2 - i;
                if(freq > maxRefKmerFreq && freq <= rescueMaxFreq) {
                    rescueEligible.push_back(refKmers[i]);
                }
                i = j2;
            }
            sort(rescueEligible.begin(), rescueEligible.end());

            // Mark windows that are part of a contiguous depleted
            // region (>=5 consecutive depleted windows, i.e. >=250bp).
            // Shorter depleted stretches are likely low-complexity or
            // short repeats where rescue causes more harm than good.
            vector<bool> isVntrDepleted(nWindows, false);
            {
                uint32_t runStart = UINT32_MAX;
                uint32_t runLen = 0;
                for(uint32_t w = 0; w <= nWindows; ++w) {
                    if(w < nWindows && windowMarkerCount[w] == 0) {
                        if(runStart == UINT32_MAX) runStart = w;
                        ++runLen;
                    } else {
                        if(runLen >= 5) {
                            for(uint32_t rw = runStart;
                                rw < runStart + runLen; ++rw)
                                isVntrDepleted[rw] = true;
                        }
                        runStart = UINT32_MAX;
                        runLen = 0;
                    }
                }
            }

            // Check which rescue-eligible k-mers have occurrences
            // in VNTR-depleted windows.
            vector<KmerId> rescueSet;
            for(size_t mi = 0; mi < refMarkerPositions.size(); ++mi) {
                const KmerId kid = refKmerIdsSpan[mi];
                if(!std::binary_search(
                       rescueEligible.begin(),
                       rescueEligible.end(), kid))
                    continue;
                const uint32_t pos = refMarkerPositions[mi].position;
                const uint32_t w = pos / windowSize;
                if(w < nWindows && isVntrDepleted[w]) {
                    rescueSet.push_back(kid);
                }
            }
            sort(rescueSet.begin(), rescueSet.end());
            rescueSet.erase(
                unique(rescueSet.begin(), rescueSet.end()),
                rescueSet.end());

            if(!rescueSet.empty()) {
                // Remove rescued k-mers from the blacklist.
                vector<KmerId> newBlacklist;
                newBlacklist.reserve(blacklist.size());
                for(const auto& kid : blacklist) {
                    if(!std::binary_search(
                           rescueSet.begin(),
                           rescueSet.end(), kid)) {
                        newBlacklist.push_back(kid);
                    }
                }
                cout << "  Rescued " << rescueSet.size()
                     << " k-mers in " << depletedWindows
                     << " depleted windows (rescueMaxFreq="
                     << rescueMaxFreq << ")" << endl;
                blacklist = std::move(newBlacklist);
            }
        }
    }

    if(blacklist.empty()) {
        const auto tEnd = std::chrono::steady_clock::now();
        const double tTotal = 1.e-9 * double(
            std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin).count());
        cout << "No non-unique reference k-mers found. Skipping filter ("
             << tTotal << " s)." << endl;
        performanceLog << timestamp
            << "removeNonUniqueReferenceMarkers completed (no-op) in "
            << tTotal << " s." << endl;
        return;
    }

    // blacklist is already sorted (subset of sorted refKmers).

    // =========================================================================
    // Phase 3: Rebuild markers and markerKmerIds, excluding blacklisted k-mers.
    // =========================================================================
    // Swap current markers to "old" for in-place rebuild.
    const string markersName = markers->getName();
    const string kmerIdsName = markerKmerIds->getName();

    auto oldMarkers = markers;
    if(!markersName.empty()) {
        oldMarkers->rename(markersName + ".old");
    }
    markers = make_shared<MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>>();
    markers->createNew(markersName, largeDataPageSize);

    auto oldMarkerKmerIds = markerKmerIds;
    if(!kmerIdsName.empty()) {
        oldMarkerKmerIds->rename(kmerIdsName + ".old");
    }
    markerKmerIds = make_shared<MemoryMapped::VectorOfVectors<KmerId, uint64_t>>();
    markerKmerIds->createNew(kmerIdsName, largeDataPageSize);

    uint64_t totalBefore = 0;
    uint64_t totalAfter = 0;

    for(ReadId readId = 0; ReadId(readId) < ReadId(readCount); ++readId) {
        const uint64_t readLen = reads->getReadRawSequenceLength(readId);

        // Process strand 0.
        const auto orientedReadId0 = OrientedReadId(readId, 0);
        const auto oldM0 = (*oldMarkers)[orientedReadId0.getValue()];
        const auto oldK0 = (*oldMarkerKmerIds)[orientedReadId0.getValue()];
        const size_t n = oldM0.size();
        totalBefore += n;

        // Determine which markers to keep.
        vector<bool> keep(n);
        size_t nKeep = 0;
        for(size_t i = 0; i < n; ++i) {
            const KmerId kmerId = oldK0[i];
            // Binary search in blacklist.
            keep[i] = !std::binary_search(blacklist.begin(), blacklist.end(), kmerId);
            if(keep[i]) ++nKeep;
        }
        totalAfter += nKeep;

        // Write strand 0 markers.
        markers->appendVector();
        markerKmerIds->appendVector();
        for(size_t i = 0; i < n; ++i) {
            if(keep[i]) {
                markers->append(oldM0[i]);
                markerKmerIds->append(oldK0[i]);
            }
        }

        // Write strand 1 markers (reverse order, reverse complement positions).
        const auto oldK1 = (*oldMarkerKmerIds)[OrientedReadId(readId, 1).getValue()];
        markers->appendVector();
        markerKmerIds->appendVector();
        // Strand 1 markers are in reverse order of strand 0.
        // keep[i] for strand 0 position i corresponds to strand 1 position (n-1-i).
        for(size_t i = n; i > 0; ) {
            --i;
            if(keep[i]) {
                // Recompute strand 1 position from strand 0 position.
                const uint32_t pos0 = oldM0[i].position;
                CompressedMarker cm;
                cm.position = static_cast<uint32_t>(readLen - k - pos0);
                markers->append(cm);
                markerKmerIds->append(oldK1[n - 1 - i]);
            }
        }
    }

    // Remove old data.
    oldMarkers->remove();
    oldMarkerKmerIds->remove();

    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double(
        std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin).count());

    cout << "Marker filtering complete in " << tTotal << " s." << endl;
    cout << "  Markers before: " << totalBefore << " (per strand)" << endl;
    cout << "  Markers after:  " << totalAfter << " (per strand)" << endl;
    cout << "  Removed:        " << (totalBefore - totalAfter) << " markers per strand." << endl;

    performanceLog << timestamp
        << "removeNonUniqueReferenceMarkers completed in " << tTotal << " s. "
        << "Removed " << (totalBefore - totalAfter) << " markers." << endl;
}



// Remove all markers from reads whose marker span covers less than
// minSpanFraction of the read length.
// Span = (lastMarkerPos + k - firstMarkerPos) / readLength.
void Assembler::filterReadsByMarkerSpanCoverage(
    double minSpanFraction,
    uint64_t threadCount)
{
    cout << timestamp << "filterReadsByMarkerSpanCoverage begins "
         << "(minSpanFraction=" << minSpanFraction << ")." << endl;

    checkMarkersAreOpen();
    DINARA_ASSERT(markerKmerIds->isOpen());

    const uint64_t k = assemblerInfo->k;
    const uint64_t readCount = reads->readCount();

    // Determine which reads pass the span coverage threshold.
    vector<bool> keepRead(readCount, true);
    uint64_t discardedCount = 0;
    uint64_t noMarkerCount = 0;

    for(ReadId readId = 0; readId < ReadId(readCount); ++readId) {
        const auto orientedReadId = OrientedReadId(readId, 0);
        const auto readMarkers = (*markers)[orientedReadId.getValue()];
        const uint64_t readLength = reads->getReadRawSequenceLength(readId);

        if(readMarkers.size() == 0) {
            keepRead[readId] = false;
            ++noMarkerCount;
            continue;
        }

        const uint32_t firstPos = readMarkers.front().position;
        const uint32_t lastPos = readMarkers.back().position;
        const uint64_t span = uint64_t(lastPos) + k - uint64_t(firstPos);
        const double fraction = double(span) / double(readLength);

        if(fraction < minSpanFraction) {
            keepRead[readId] = false;
            ++discardedCount;
        }
    }

    cout << timestamp << "filterReadsByMarkerSpanCoverage: "
         << discardedCount << " reads discarded (span < "
         << minSpanFraction << "), "
         << noMarkerCount << " reads had no markers, "
         << (readCount - discardedCount - noMarkerCount) << " reads kept." << endl;

    if(discardedCount == 0 && noMarkerCount == 0) {
        cout << timestamp << "filterReadsByMarkerSpanCoverage: nothing to filter." << endl;
        return;
    }

    // Rebuild markers and markerKmerIds, keeping only passing reads.
    const string markersName = markers->getName();
    const string kmerIdsName = markerKmerIds->getName();

    auto oldMarkers = markers;
    auto oldMarkerKmerIds = markerKmerIds;

    if(!markersName.empty()) {
        oldMarkers->rename(markersName + ".spanfilter.old");
    }
    if(!kmerIdsName.empty()) {
        oldMarkerKmerIds->rename(kmerIdsName + ".spanfilter.old");
    }

    markers = make_shared<MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>>();
    markerKmerIds = make_shared<MemoryMapped::VectorOfVectors<KmerId, uint64_t>>();

    markers->createNew(markersName, largeDataPageSize);
    markerKmerIds->createNew(kmerIdsName, largeDataPageSize);

    // Pass 1: count markers per oriented read.
    const uint64_t orientedReadCount = 2 * readCount;
    markers->beginPass1(orientedReadCount);
    markerKmerIds->beginPass1(orientedReadCount);

    for(ReadId readId = 0; readId < ReadId(readCount); ++readId) {
        for(Strand strand = 0; strand < 2; ++strand) {
            const auto oid = OrientedReadId(readId, strand);
            const uint64_t n = keepRead[readId]
                ? (*oldMarkers)[oid.getValue()].size()
                : 0;
            markers->incrementCount(oid.getValue(), n);
            markerKmerIds->incrementCount(oid.getValue(), n);
        }
    }

    // Pass 2: copy markers for passing reads.
    markers->beginPass2();
    markerKmerIds->beginPass2();

    // store() writes in reverse order (decrements count as index),
    // so iterate backwards to preserve the original marker order.
    for(ReadId readId = 0; readId < ReadId(readCount); ++readId) {
        if(!keepRead[readId]) continue;
        for(Strand strand = 0; strand < 2; ++strand) {
            const auto oid = OrientedReadId(readId, strand);
            const auto oldM = (*oldMarkers)[oid.getValue()];
            const auto oldK = (*oldMarkerKmerIds)[oid.getValue()];
            const uint64_t n = oldM.size();
            for(uint64_t i = n; i > 0; --i) {
                markers->store(oid.getValue(), oldM[i - 1]);
                markerKmerIds->store(oid.getValue(), oldK[i - 1]);
            }
        }
    }

    markers->endPass2(false);
    markerKmerIds->endPass2(false);

    // Clean up old data.
    oldMarkers->remove();
    oldMarkerKmerIds->remove();

    cout << timestamp << "filterReadsByMarkerSpanCoverage done." << endl;
}
