// PngImage.hpp must be included first because of png issues on Ubuntu 16.04.
#include "PngImage.hpp"

// Dinara.
#include "Assembler.hpp"
#include "Alignment.hpp"
#include "AlignmentGraph.hpp"
#include "Align4.hpp"
#include "AssemblerOptions.hpp"
#include "compressAlignment.hpp"
#include "performanceLog.hpp"
#include "ProjectedAlignment.hpp"
#include "Reads.hpp"
#include "HifiasmCigarImport.hpp"
#include "hifiasmCoordinateTransforms.hpp"
#include "overlapClassification.hpp"
#include "span.hpp"
#include "timestamp.hpp"

// Standard libraries.
#include "chrono.hpp"
#include <cmath>
#include <filesystem>
#include "iterator.hpp"
#include "tuple.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>

using namespace dinara;
using namespace std;

void Assembler::deriveChainFromInterval(
    const array<OrientedReadId, 2>& orientedReadIds,
    uint32_t read0Begin, uint32_t read0End,
    uint32_t read1Begin, uint32_t read1End,
    vector< array<uint32_t, 2> >& ordinals) const
{
    ordinals.clear();

    // Markers and their KmerIds share the same oriented indexing (index ==
    // ordinal) and are position-sorted ascending in each oriented frame.
    const auto markers0    = (*markers)[orientedReadIds[0].getValue()];
    const auto markers1    = (*markers)[orientedReadIds[1].getValue()];
    const auto kmerIds0    = (*markerKmerIds)[orientedReadIds[0].getValue()];
    const auto kmerIds1    = (*markerKmerIds)[orientedReadIds[1].getValue()];

    // Collect (KmerId, ordinal) for markers whose position lies in the overlap
    // box on each side. Positions are sorted, so we can stop once past the end.
    // Pair type kept local; sorted by KmerId for the intersection below.
    using KmerOrdinal = pair<KmerId, uint32_t>;
    vector<KmerOrdinal> box0, box1;
    auto collect = [](
        const span<const CompressedMarker>& m,
        const span<const KmerId>& kids,
        uint32_t begin, uint32_t end,
        vector<KmerOrdinal>& out)
    {
        const uint32_t n = uint32_t(m.size());
        for(uint32_t ord = 0; ord < n; ++ord) {
            const uint32_t pos = m[ord].position;
            if(pos < begin) continue;
            if(pos >= end) break;              // sorted: no later marker qualifies
            out.emplace_back(kids[ord], ord);
        }
    };
    collect(markers0, kmerIds0, read0Begin, read0End, box0);
    collect(markers1, kmerIds1, read1Begin, read1End, box1);

    if(box0.empty() || box1.empty()) {
        return;
    }

    // Sort each side by KmerId, then mark KmerIds that occur more than once on a
    // side: within the box such anchors are ambiguous (repeats) and are dropped,
    // exactly like a seed-chaining step ignoring high-multiplicity minimizers.
    auto byKmer = [](const KmerOrdinal& a, const KmerOrdinal& b) {
        return a.first < b.first;
    };
    sort(box0.begin(), box0.end(), byKmer);
    sort(box1.begin(), box1.end(), byKmer);

    // Two-pointer intersection over sorted-by-KmerId lists. Only KmerIds that
    // appear EXACTLY once on each side produce an anchor pair {ord0, ord1}.
    vector< array<uint32_t, 2> > anchors;
    {
        size_t i = 0, j = 0;
        const size_t n0 = box0.size(), n1 = box1.size();
        while(i < n0 && j < n1) {
            if(box0[i].first < box1[j].first) {
                ++i;
            } else if(box1[j].first < box0[i].first) {
                ++j;
            } else {
                // Equal KmerId. Count run length on each side.
                const KmerId key = box0[i].first;
                size_t iEnd = i, jEnd = j;
                while(iEnd < n0 && box0[iEnd].first == key) ++iEnd;
                while(jEnd < n1 && box1[jEnd].first == key) ++jEnd;
                if((iEnd - i) == 1 && (jEnd - j) == 1) {
                    anchors.push_back({box0[i].second, box1[j].second});
                }
                i = iEnd;
                j = jEnd;
            }
        }
    }

    if(anchors.empty()) {
        return;
    }

    // Order anchors by read0 ordinal. Inside hifiasm's validated box the overlap
    // is colinear, so read1 ordinals are then already almost-monotonic; enforce a
    // strictly-increasing guard on read1 to yield a valid chain and drop the rare
    // out-of-order anchor (a mismapped repeat that slipped the multiplicity test).
    sort(anchors.begin(), anchors.end(),
        [](const array<uint32_t, 2>& a, const array<uint32_t, 2>& b) {
            return a[0] < b[0];
        });

    ordinals.reserve(anchors.size());
    uint32_t lastOrd1 = 0;
    bool haveLast = false;
    for(const array<uint32_t, 2>& a : anchors) {
        if(haveLast && a[1] <= lastOrd1) {
            continue;                         // keep read1 strictly increasing
        }
        ordinals.push_back(a);
        lastOrd1 = a[1];
        haveLast = true;
    }
}

void Assembler::computeBaseAlignmentsAndStore(
    const AlignOptions& alignOptions,
    uint64_t threadCount
) {
    const auto tBegin = steady_clock::now();
    const size_t candidateCount = alignmentCandidates.candidates.size();

    cout << timestamp << "Begin computing alignments with evidence for ";
    cout << candidateCount << " alignment candidates." << endl;

    // Check that we have what we need.
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentCandidatesAreOpen();

    // Store parameters so they are accessible to the threads.
    auto& data = computeAlignmentsData;
    data.alignOptions = &alignOptions;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Pick the batch size for computing alignments.
    // Use the same logic as legacy to ensure similar load balancing.
    size_t batchSize = 10;
    if(batchSize > candidateCount/threadCount) {
        batchSize = candidateCount/threadCount;
    }
    if(batchSize == 0) {
        batchSize = 1;
    }

    // Prepare data structures for each thread.
    data.threadAlignmentData.resize(threadCount);
    for(auto& v : data.threadAlignmentData) {
        v.clear();
    }
    data.threadCompressedAlignments.resize(threadCount);
    data.threadCigarStores.resize(threadCount);
    for(auto& store : data.threadCigarStores) {
        store.clear();
    }
    
    // Always resize these to avoid SIGSEGV during aggregation.
    data.threadProjectedAlignmentTime.assign(threadCount, 0.0);
    data.threadCollectionTime.assign(threadCount, 0.0);
    data.threadFilteredByErrorRate.assign(threadCount, 0);
    data.threadFilteredByErrorRateGap.assign(threadCount, 0);
    data.threadFilteredByGapCount.assign(threadCount, 0);

    // If needed for variant clustering, resize this as well.
    if (assemblerInfo->readGraphCreationMethod == 5) {
        data.threadVariantClusteringPositionPairs.resize(threadCount);
    } else {
        data.threadVariantClusteringPositionPairs.clear();
    }
    
    performanceLog << timestamp << "Alignment computation begins (Unified Chaining & Evidence path)." << endl;
    cout << timestamp << "Alignment computation begins (Unified Chaining & Evidence path)." << endl;
    
    setupLoadBalancing(candidateCount, batchSize);
    runThreads(&Assembler::computeBaseAlignmentsAndStoreThreadFunction, threadCount);
    
    performanceLog << timestamp << "Alignment computation completed." << endl;
    cout << timestamp << "Alignment computation completed." << endl;

    // Store the alignments found by each thread into global shared memory.
    performanceLog << timestamp << "Storing the alignment found by each thread." << endl;

    // Pre-size global containers to avoid repeated remaps/reallocations.
    uint64_t totalAlignments = 0;
    uint64_t totalCigarTokens = 0;
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        totalAlignments += data.threadAlignmentData[threadId].size();
        const OverlapCigarStore& localCigarStore = data.threadCigarStores[threadId];
        totalCigarTokens += localCigarStore.tokenCount();
    }

    alignmentData.createNew(largeDataName("AlignmentData"), largeDataPageSize, 0, totalAlignments);
    compressedAlignments.createNew(largeDataName("CompressedAlignments"), largeDataPageSize);

    overlapCigarStore.clear();
    overlapCigarStore.reserve(totalCigarTokens);

    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector<AlignmentData>& threadAlignmentData = data.threadAlignmentData[threadId];
        const size_t idShift = alignmentData.size(); // Current global count serves as offset.
        if(!threadAlignmentData.empty()) {
            // Merge thread-local packed CIGAR arena first to get the offset shift.
            const uint32_t cigarOffsetShift = overlapCigarStore.merge(data.threadCigarStores[threadId]);
            data.threadCigarStores[threadId].clear();

            alignmentData.resize(idShift + threadAlignmentData.size());
            std::copy(
                threadAlignmentData.begin(),
                threadAlignmentData.end(),
                alignmentData.begin() + idShift);
            for(size_t i = 0; i < threadAlignmentData.size(); ++i) {
                alignmentData[idShift + i].info.alignmentId += idShift;
                if(alignmentData[idShift + i].info.cigarOffset != uint32_t(-1)) {
                    alignmentData[idShift + i].info.cigarOffset += cigarOffsetShift;
                }
            }
        }

        const auto& threadCompressedAlignments = *data.threadCompressedAlignments[threadId];
        const auto size = threadCompressedAlignments.size();
        for(size_t i=0; i<size; i++) {
            compressedAlignments.appendVector(
                threadCompressedAlignments[i].begin(),
                threadCompressedAlignments[i].end()
            );
        }

    }

    // Release unused allocated memory.
    alignmentData.unreserve();
    compressedAlignments.unreserve();

    // Create the alignment table (mapping from read to its alignments).
    // This is required for subsequent filtering phases.
    performanceLog << timestamp << "Creating alignment table." << endl;
    computeAlignmentTable();

    // Report CIGAR store memory usage.
    {
        const size_t cigarTotalBytes = overlapCigarStore.memoryUsage();
        auto mb = [](size_t b) { return b / (1024.0 * 1024.0); };
        cout << timestamp << "OverlapCigarStore: "
             << overlapCigarStore.tokenCount() << " tokens, "
             << mb(cigarTotalBytes) << " MB." << endl;
    }

    const auto tEnd = steady_clock::now();
    const double elapsedSeconds = seconds(tEnd - tBegin);
    performanceLog << timestamp << "Done computing alignments. Elapsed time: " << elapsedSeconds << " s." << endl;
    cout << timestamp << "Done computing alignments. Elapsed time: " << elapsedSeconds << " s." << endl;
}




void Assembler::computeBaseAlignmentsAndStoreThreadFunction(size_t threadId) {
    auto& data = computeAlignmentsData;
    const AlignOptions& alignOptions = *data.alignOptions;
    auto& threadAlignmentData = data.threadAlignmentData[threadId];
    const auto& candidates = alignmentCandidates.candidates;
    // Every candidate carries an imported hifiasm CIGAR (candidates and the
    // CIGAR store come from the SAME deduped import list, keyed identically).
    // The marker-ordinal chain is DERIVED per candidate from the hifiasm overlap
    // box (deriveChainFromInterval); there is no separate chaining DP.
    const uint32_t markerK = uint32_t(assemblerInfo->k);
    const bool collectProjectedTiming = (assemblerInfo->readGraphCreationMethod == 5);
    const size_t minAlignedMarkerCount = (alignOptions.minAlignedMarkerCount > 0) ?
        size_t(alignOptions.minAlignedMarkerCount) : 0;
    const double maxErrorRate = alignOptions.maxErrorRate;
    const int64_t dpMatchScore = alignOptions.overlapDpMatchScore;
    const int64_t dpMismatchScore = alignOptions.overlapDpMismatchScore;
    const int64_t dpGapOpen1 = alignOptions.overlapDpGapOpen1;
    const int64_t dpGapExtend1 = alignOptions.overlapDpGapExtend1;
    
    // Initialize compressed alignment storage for this thread.
    data.threadCompressedAlignments[threadId] = make_shared<MemoryMapped::VectorOfVectors<char, uint64_t>>();
    auto& thisThreadCompressedAlignments = *data.threadCompressedAlignments[threadId];
    thisThreadCompressedAlignments.createNew(
        largeDataName("Thread" + to_string(threadId) + "-CompressedAlignments"),
        largeDataPageSize);

    OverlapCigarStore& cigarStore = data.threadCigarStores[threadId];
    string compressedAlignment;
    array<OrientedReadId, 2> orientedReadIds;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t candidateIndex = begin; candidateIndex != end; candidateIndex++) {
            const OrientedReadPair& candidate = candidates[candidateIndex];

            // The marker-ordinal chain is derived from hifiasm's overlap box
            // (deriveChainFromInterval) after the CIGAR walk below. The chain
            // size is only known after that derivation, so the
            // minAlignedMarkerCount support filter is applied there.
            orientedReadIds[0] = OrientedReadId(candidate.readIds[0], 0);
            orientedReadIds[1] = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);
            const array<LongBaseSequenceView, 2> sequenceViews = {
                reads->getRead(orientedReadIds[0].getReadId()),
                reads->getRead(orientedReadIds[1].getReadId())
            };
            const array<span<const CompressedMarker>, 2> markerSpans = {
                (*markers)[orientedReadIds[0].getValue()],
                (*markers)[orientedReadIds[1].getValue()]
            };

            // Compute projected alignment metrics and sparse diffs (mismatches/indels).
            steady_clock::time_point tProjStart;
            if(collectProjectedTiming) {
                tProjStart = steady_clock::now();
            }

            // Reuse hifiasm's base alignment: reframe hifiasm's CIGAR into
            // dinara's read0/read1 canonical frame and walk it clipped to the
            // marker interval. Falls back to A*PA2 (Method QuickRawSparse) only
            // if the imported CIGAR is missing/unusable for this pair.
            bool usedHifiasmCigar = false;
            // The chain is derived from the overlap box (deriveChainFromInterval),
            // so the ctor gets an empty local Alignment that we fill below; under
            // Method::None the ctor does not read it.
            Alignment directAlignment;
            ProjectedAlignment projectedAlignment(
                markerK,
                orientedReadIds,
                sequenceViews,
                directAlignment,
                markerSpans,
                ProjectedAlignment::Method::None,
                dpMatchScore,
                dpMismatchScore,
                dpGapOpen1,
                dpGapExtend1,
                &cigarStore);

            {
                // Candidate readIds are canonical (readIds[0] < readIds[1]); the
                // pair key and strand match the imported-CIGAR store keys.
                const uint64_t pairKey =
                    (uint64_t(candidate.readIds[0]) << 32) | uint64_t(candidate.readIds[1]);
                const HifiasmImportedCigarStore::Record* rec =
                    hifiasmImportedCigarStore.find(pairKey, candidate.isSameStrand);
                if(rec != nullptr && rec->cigarTokenCount > 0) {
                    // Reframe the native hifiasm CIGAR into read0/read1 canonical
                    // frame. read lengths are needed for the reverse-strand /
                    // id-swap cases.
                    const uint32_t qLen =
                        uint32_t(reads->getRead(ReadId(rec->readIdQ)).baseCount);
                    const uint32_t tLen =
                        uint32_t(reads->getRead(ReadId(rec->readIdT)).baseCount);
                    const NormalizedHifiasmCigar norm = normalizeHifiasmCigar(
                        hifiasmImportedCigarStore.tokensOf(*rec),
                        rec->readIdQ, rec->readIdT,
                        rec->qStart, rec->qEnd, rec->tStart, rec->tEnd,
                        qLen, tLen, rec->isSameStrand);

                    // Clip window is the CIGAR's own read0 span (forward coords).
                    // The marker-ordinal chain is derived separately from the
                    // hifiasm overlap interval (deriveChainFromInterval), not from
                    // this walk; the walk only produces statistics and the CIGAR.
                    usedHifiasmCigar = projectedAlignment.constructFromHifiasmCigar(
                        span<const CigarToken>(norm.tokens.data(), norm.tokens.size()),
                        norm.read0Start, norm.read1Start,
                        norm.read0Start, norm.read0End);

                    // Derive the marker-ordinal chain from the overlap box hifiasm
                    // validated. Both intervals are in the oriented frames of
                    // orientedReadIds: read0 forward [read0Start,read0End); read1
                    // in alignment orientation [read1Start,read1End) (norm already
                    // reframed for reverse overlaps). No chaining DP -- the box is
                    // colinear.
                    if(usedHifiasmCigar) {
                        deriveChainFromInterval(
                            orientedReadIds,
                            norm.read0Start, norm.read0End,
                            norm.read1Start, norm.read1End,
                            directAlignment.ordinals);
                    }
                } else if(rec != nullptr) {
                    // Interval-only record (myloasm marker overlap path): the
                    // pair has a validated interval but no base CIGAR. Reframe
                    // the interval into the read0/read1 canonical frame exactly
                    // as above (normalizeHifiasmCigar produces empty tokens here)
                    // and derive the marker chain from it. constructQuickRawSparse
                    // below then builds the per-segment CIGAR from that chain with
                    // A*PA2 -- the same base alignment the hifiasm path gets for
                    // free, but computed from myloasm's interval.
                    const uint32_t qLen =
                        uint32_t(reads->getRead(ReadId(rec->readIdQ)).baseCount);
                    const uint32_t tLen =
                        uint32_t(reads->getRead(ReadId(rec->readIdT)).baseCount);
                    const NormalizedHifiasmCigar norm = normalizeHifiasmCigar(
                        hifiasmImportedCigarStore.tokensOf(*rec),
                        rec->readIdQ, rec->readIdT,
                        rec->qStart, rec->qEnd, rec->tStart, rec->tEnd,
                        qLen, tLen, rec->isSameStrand);
                    deriveChainFromInterval(
                        orientedReadIds,
                        norm.read0Start, norm.read0End,
                        norm.read1Start, norm.read1End,
                        directAlignment.ordinals);
                }
                if(!usedHifiasmCigar) {
                    // No usable CIGAR: build the base alignment from the derived
                    // marker chain (interval-only path) or, if no chain exists,
                    // recompute from scratch. constructQuickRawSparse walks
                    // directAlignment.ordinals; on the interval-only path this is
                    // the chain from deriveChainFromInterval above.
                    projectedAlignment.constructQuickRawSparse();
                }
            }
            if(collectProjectedTiming) {
                data.threadProjectedAlignmentTime[threadId] += seconds(steady_clock::now() - tProjStart);
            }

            // The chain size is only known after deriveChainFromInterval. Drop
            // pairs with no chain or fewer than Align.minAlignedMarkerCount anchors.
            if(directAlignment.ordinals.empty()) {
                continue;
            }
            if(minAlignedMarkerCount > 0 &&
                directAlignment.ordinals.size() < minAlignedMarkerCount) {
                continue;
            }
            
            // Error rate filtering.
            const double projectedErrorRate = projectedAlignment.errorRate();
            if(projectedErrorRate > maxErrorRate) {
                data.threadFilteredByErrorRate[threadId]++;
                continue;
            }

            // Create alignment info summary from the derived-chain alignment.
            AlignmentInfo alignmentInfo(directAlignment,
                uint32_t(markerSpans[0].size()),
                uint32_t(markerSpans[1].size())
            );

            alignmentInfo.errorRate = float(projectedErrorRate);
            alignmentInfo.mismatchCount = uint32_t(projectedAlignment.mismatchCount);
            const double projectedGapErrorRate = projectedAlignment.errorRateGaps();
            alignmentInfo.errorRateGaps = float(projectedGapErrorRate);
            alignmentInfo.gapCount = uint32_t(projectedAlignment.totalIndelBaseCount);
            alignmentInfo.gapEventCount = uint32_t(projectedAlignment.totalGapEventCount);
            alignmentInfo.dpScore = projectedAlignment.totalDpScore;

            AlignmentData thisAlignmentData(candidate, alignmentInfo);
            
            // Store CIGAR boundary positions (marker-based, not extended).
            // read0 is always strand 0 → forward coordinates.
            thisAlignmentData.qs = projectedAlignment.cigarRead0Start;
            thisAlignmentData.qe = projectedAlignment.cigarRead0End;
            // read1 may be strand 1 → convert oriented coords to forward.
            if (candidate.isSameStrand) {
                thisAlignmentData.ts = projectedAlignment.cigarRead1Start;
                thisAlignmentData.te = projectedAlignment.cigarRead1End;
            } else {
                const uint32_t tLen = uint32_t(sequenceViews[1].baseCount);
                thisAlignmentData.ts = tLen - projectedAlignment.cigarRead1End;
                thisAlignmentData.te = tLen - projectedAlignment.cigarRead1Start;
            }

            // Hifiasm-parity EXTENDED coordinates (see AlignmentData::extendedQs
            // doc comment): diagonally extrapolate the tight span above out to
            // read boundaries, matching hifiasm's ma_hit_t convention. Used by
            // any containment/dovetail-type classification (e.g. flagContainedReads),
            // never by consumers that want the real/tight aligned span.
            thisAlignmentData.extendedQs = thisAlignmentData.qs;
            thisAlignmentData.extendedQe = thisAlignmentData.qe;
            thisAlignmentData.extendedTs = thisAlignmentData.ts;
            thisAlignmentData.extendedTe = thisAlignmentData.te;
            extendOverlapToReadBoundaries(
                uint32_t(sequenceViews[0].baseCount),
                uint32_t(sequenceViews[1].baseCount),
                candidate.isSameStrand,
                thisAlignmentData.extendedQs, thisAlignmentData.extendedQe,
                thisAlignmentData.extendedTs, thisAlignmentData.extendedTe);

            thisAlignmentData.hasLargeIndel = projectedAlignment.hasLargeIndel;
            thisAlignmentData.informativeHetSiteCount0 = 0;
            thisAlignmentData.informativeHetSiteCount1 = 0;
            thisAlignmentData.informativeHetSiteScore = 0;
            thisAlignmentData.deleteReasons0 = AlignmentData::DeleteReasonNone;
            thisAlignmentData.deleteReasons1 = AlignmentData::DeleteReasonNone;

            // Per-overlap CIGAR is stored in the OverlapCigarStore arena.
            // alignmentId is kept as the running alignment index (later shifted
            // to a global index during the merge below).
            thisAlignmentData.info.alignmentId = threadAlignmentData.size();
            thisAlignmentData.info.cigarOffset     = projectedAlignment.cigarOffset;
            thisAlignmentData.info.cigarTokenCount = projectedAlignment.cigarTokenCount;

            threadAlignmentData.push_back(thisAlignmentData);
            
            dinara::compress(directAlignment, compressedAlignment);
            thisThreadCompressedAlignments.appendVector(
                compressedAlignment.begin(),
                compressedAlignment.end()
            );
        }
    }
    
    thisThreadCompressedAlignments.unreserve();
}



void Assembler::filterSecondaryAlignmentsPerReadPair(
    uint64_t threadCount,
    bool requireNonRedundantOnBothReads)
{
    // --- Hifiasm-style per-read-pair redundancy filtering ---
    // For each read r0, consider all overlaps to a given partner read (r1).
    // Keep the strongest chains and delete redundant/secondary ones:
    // - score ratio filter (hifiasm -p 0.8)
    // - overlap-on-r0 filter (hifiasm -M 0.5)
    // Optionally also require non-redundancy on the partner read interval (ts/te).

    const auto tBegin = steady_clock::now();
    cout << timestamp << "Filtering secondary alignments per read pair begins." << endl;

    removedSecondaryAlignmentCount = 0;
    removedSecondaryAlignmentBySymmetryOnlyCount = 0;
    filterSecondaryRequireNonRedundantOnBothReads = requireNonRedundantOnBothReads;

    // Run threads.
    const uint64_t readCount = reads->readCount();
    setupLoadBalancing(readCount, 100); // Batch size heuristic.
    runThreads(&Assembler::filterSecondaryAlignmentsPerReadPairThreadFunction, threadCount);

    const auto tEnd = steady_clock::now();
    const double tSeconds = seconds(tEnd - tBegin);
    cout << timestamp << "Filtering secondary alignments per read pair ends in " << tSeconds << " s." << endl;
    cout << timestamp << "Removed " << removedSecondaryAlignmentCount << " redundant alignments (marked isDeleted).";
    if(filterSecondaryRequireNonRedundantOnBothReads) {
        cout << " symmetryOnlyRemoved=" << removedSecondaryAlignmentBySymmetryOnlyCount;
    }
    cout << endl;
}



void Assembler::filterSecondaryAlignmentsPerReadPairThreadFunction(size_t)
{
    uint64_t begin, end;
    uint64_t localRemovedCount = 0;
    uint64_t localSymmetryOnlyRemovedCount = 0;

    // Reuse buffers across reads to avoid repeated allocations.
    struct CandidateInfo {
        uint32_t alignmentId;
        int64_t score;
        uint64_t qs;
        uint64_t qe;
        uint64_t ts;
        uint64_t te;
        bool keptOld = false; // kept by legacy (r0-only) logic
        bool keptNew = false; // kept by current logic (with optional symmetry)
    };
    vector<CandidateInfo> group;
    vector<CandidateInfo> kept;
    group.reserve(16);
    kept.reserve(16);

    // Default batch size from setupLoadBalancing.
    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); r0++) {
            const OrientedReadId orientedR0(r0, 0);
            const auto& table = alignmentTable[orientedR0.getValue()];
            if(table.empty()) {
                continue;
            }

            group.clear();
            kept.clear();

            ReadId currentTarget = invalidReadId;

            auto flushGroup = [&]() {
                if(group.empty()) {
                    return;
                }
                if(group.size() == 1) {
                    // A single alignment for this read pair cannot be redundant.
                    return;
                }

                // Sort group by score descending.
                // Hifiasm tie-break when shared_seed ties: prefer smaller overlapLen.
                sort(group.begin(), group.end(), [](const CandidateInfo& a, const CandidateInfo& b) {
                    if(a.score != b.score) {
                        return a.score > b.score;
                    }
                    const uint64_t lenA = a.qe - a.qs;
                    const uint64_t lenB = b.qe - b.qs;
                    if(lenA != lenB) {
                        return lenA < lenB;
                    }
                    return a.alignmentId < b.alignmentId;
                });

                const int64_t bestScore = group[0].score;
                auto failsScoreRatio = [&](const CandidateInfo& cand) -> bool {
                    // 1) Score ratio filter (hifiasm -p 0.8), done with integer math:
                    // cand.score < bestScore * 0.8  <=>  cand.score * 5 < bestScore * 4
                    return cand.score * 5 < bestScore * 4;
                };
                auto overlapsMoreThanHalf = [&](uint64_t s0, uint64_t e0, uint64_t s1, uint64_t e1, uint64_t len) -> bool {
                    if(len == 0) {
                        return false;
                    }
                    const uint64_t oStart = std::max(s0, s1);
                    const uint64_t oEnd = std::min(e0, e1);
                    if(oEnd <= oStart) {
                        return false;
                    }
                    const uint64_t oLen = oEnd - oStart;
                    return oLen * 2 > len;
                };
                auto overlapsKeptOnR0 = [&](const CandidateInfo& cand, const vector<CandidateInfo>& keptIntervals) -> bool {
                    const uint64_t len = cand.qe - cand.qs;
                    for(const auto& existing : keptIntervals) {
                        if(overlapsMoreThanHalf(cand.qs, cand.qe, existing.qs, existing.qe, len)) {
                            return true;
                        }
                    }
                    return false;
                };
                auto overlapsKeptOnR1 = [&](const CandidateInfo& cand, const vector<CandidateInfo>& keptIntervals) -> bool {
                    const uint64_t len = cand.te - cand.ts;
                    for(const auto& existing : keptIntervals) {
                        if(overlapsMoreThanHalf(cand.ts, cand.te, existing.ts, existing.te, len)) {
                            return true;
                        }
                    }
                    return false;
                };

                // Fast path: legacy behavior (r0-only redundancy).
                if(!filterSecondaryRequireNonRedundantOnBothReads) {
                    for(auto& cand : group) {
                        cand.keptNew = false;
                    }
                    kept.clear();
                    for(auto& cand : group) {
                        if(failsScoreRatio(cand)) {
                            continue;
                        }
                        if(overlapsKeptOnR0(cand, kept)) {
                            continue;
                        }
                        cand.keptNew = true;
                        kept.push_back(cand);
                    }
                    for(const auto& cand : group) {
                        if(!cand.keptNew) {
                            alignmentData[cand.alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
                            localRemovedCount++;
                        }
                    }
                    return;
                }

                // First pass: legacy behavior (r0-only redundancy).
                for(auto& cand : group) {
                    cand.keptOld = false;
                }
                kept.clear();
                for(auto& cand : group) {
                    if(failsScoreRatio(cand)) {
                        continue;
                    }
                    if(overlapsKeptOnR0(cand, kept)) {
                        continue;
                    }
                    cand.keptOld = true;
                    kept.push_back(cand);
                }

                // Second pass: actual behavior (optionally symmetric).
                for(auto& cand : group) {
                    cand.keptNew = false;
                }
                kept.clear();
                for(auto& cand : group) {
                    if(failsScoreRatio(cand)) {
                        continue;
                    }
                    if(overlapsKeptOnR0(cand, kept)) {
                        continue;
                    }
                    if(filterSecondaryRequireNonRedundantOnBothReads && overlapsKeptOnR1(cand, kept)) {
                        continue;
                    }
                    cand.keptNew = true;
                    kept.push_back(cand);
                }

                // Mark deletions and count symmetry-only removals.
                for(const auto& cand : group) {
                    if(!cand.keptNew) {
                        alignmentData[cand.alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
                        localRemovedCount++;
                        if(cand.keptOld) {
                            localSymmetryOnlyRemovedCount++;
                        }
                    }
                }
            };

            // alignmentTable[OrientedReadId(r0,0)] is already sorted by partner OrientedReadId,
            // so partner ReadIds are contiguous (strand 0/1 are adjacent for the same read).
            // Skip the prefix where the partner oriented read id is < orientedR0 (partner read id < r0).
            // This avoids scanning alignments that would be skipped by canonical ordering checks.
            size_t firstIndex = 0;
            {
                size_t lo = 0;
                size_t hi = table.size();
                while(lo < hi) {
                    const size_t mid = (lo + hi) / 2;
                    const uint32_t alignmentId = table[mid];
                    const OrientedReadId other = alignmentData[alignmentId].getOther(orientedR0);
                    if(other < orientedR0) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                firstIndex = lo;
            }

            for(size_t tableIndex = firstIndex; tableIndex < table.size(); ++tableIndex) {
                const uint32_t alignmentId = table[tableIndex];
                const auto& ad = alignmentData[alignmentId];

                // This filter is intended to run after phasing/EC parity:
                // only consider overlaps that are cis on BOTH reads.
                // (Trans overlaps already carry DeleteReasonPhase and should not compete for "best" here.)
                if(!ad.isCisByBothSides()) {
                    continue;
                }

                // Skip alignments already deleted for other reasons. We only want to mark *additional*
                // redundancies among alignments that are still eligible for the graph.
                if(!ad.keptByBothSides()) {
                    continue;
                }

                // Mark self-overlaps as deleted immediately.
                if(ad.readIds[0] == ad.readIds[1]) {
                    alignmentData[alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
                    localRemovedCount++;
                    continue;
                }

                // Canonical ordering: only process from the smaller ReadId side.
                // AlignmentData is expected to be stored with readIds[0] < readIds[1].
                // So for this r0, we only act when r0 == readIds[0].
                if(ad.readIds[0] != r0) {
                    continue;
                }
                const ReadId target = ad.readIds[1];

                if(target != currentTarget) {
                    flushGroup();
                    group.clear();
                    currentTarget = target;
                }

                const auto& info = ad.info;

                // Overlap length on r0.
                const uint64_t overlapLenOnR0 = uint64_t(ad.qe) - uint64_t(ad.qs);

                // Hifiasm keeps the best chain(s) by overlap_region.shared_seed (chaining DP score),
                // not by a base-level alignment DP score. We use a proxy score derived from the
                // marker alignment strength (see AlignmentInfo::hifiasmSharedSeedScoreProxy()).
                (void)overlapLenOnR0;
                const int64_t score = info.hifiasmSharedSeedScoreProxy();

                // Overlap interval on r0 (forward coordinates).
                const uint64_t qs = ad.qs;
                const uint64_t qe = ad.qe;
                const uint64_t ts = ad.ts;
                const uint64_t te = ad.te;

                group.push_back(CandidateInfo{alignmentId, score, qs, qe, ts, te});
            }

            flushGroup();
        }
    }

    removedSecondaryAlignmentCount += localRemovedCount;
    removedSecondaryAlignmentBySymmetryOnlyCount += localSymmetryOnlyRemovedCount;
}



void Assembler::keepOnlyBestAlignmentPerReadPairByDpScore(uint64_t threadCount)
{
    // For each read pair (r0,r1), multiple alignment chains can exist.
    // This keeps only the single best-scoring chain (hifiasm shared_seed proxy score),
    // and deletes the rest as redundant/secondary.
    //
    // This is intentionally aggressive. It can remove legitimate alternative placements
    // (repeats / multi-mapping) and should be used as an experimental clean-up step.
    const auto tBegin = steady_clock::now();
    cout << timestamp << "Keeping only best alignment per read pair by DP score begins." << endl;

    removedSecondaryAlignmentCount = 0;
    removedSecondaryAlignmentBySymmetryOnlyCount = 0;

    const uint64_t readCount = reads->readCount();
    setupLoadBalancing(readCount, 100); // Batch size heuristic.
    runThreads(&Assembler::keepOnlyBestAlignmentPerReadPairByDpScoreThreadFunction, threadCount);

    const auto tEnd = steady_clock::now();
    const double tSeconds = seconds(tEnd - tBegin);
    cout << timestamp << "Keeping only best alignment per read pair by DP score ends in " << tSeconds << " s." << endl;
    cout << timestamp << "Removed " << removedSecondaryAlignmentCount << " alignments (marked DeleteReasonSecondary)." << endl;
}



void Assembler::keepOnlyBestAlignmentPerReadPairByDpScoreThreadFunction(size_t)
{
    uint64_t begin = 0, end = 0;
    uint64_t localRemovedCount = 0;

    // Reuse buffers across reads to avoid repeated allocations.
    struct CandidateInfo {
        uint32_t alignmentId = 0;
        int64_t score = 0;
        ReadId target = invalidReadId;
        uint64_t overlapLen = 0;
    };
    vector<CandidateInfo> group;
    group.reserve(16);

    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); r0++) {
            const OrientedReadId orientedR0(r0, 0);
            const auto& table = alignmentTable[orientedR0.getValue()];
            if(table.empty()) {
                continue;
            }

            group.clear();
            ReadId currentTarget = invalidReadId;

            auto flushGroup = [&]() {
                if(group.size() <= 1) {
                    return;
                }
                // Keep the single best by score.
                size_t bestIndex = 0;
                for(size_t i = 1; i < group.size(); ++i) {
                    if(group[i].score > group[bestIndex].score) {
                        bestIndex = i;
                    } else if(group[i].score == group[bestIndex].score && group[i].overlapLen < group[bestIndex].overlapLen) {
                        // Hifiasm tie-break when shared_seed ties: prefer smaller overlapLen.
                        bestIndex = i;
                    }
                }
                for(size_t i = 0; i < group.size(); ++i) {
                    if(i == bestIndex) continue;
                    alignmentData[group[i].alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
                    ++localRemovedCount;
                }
            };

            // Same canonical prefix skip as filterSecondaryAlignmentsPerReadPair to reduce scanning.
            size_t firstIndex = 0;
            {
                size_t lo = 0;
                size_t hi = table.size();
                while(lo < hi) {
                    const size_t mid = (lo + hi) / 2;
                    const uint32_t alignmentId = table[mid];
                    const OrientedReadId other = alignmentData[alignmentId].getOther(orientedR0);
                    if(other < orientedR0) {
                        lo = mid + 1;
                    } else {
                        hi = mid;
                    }
                }
                firstIndex = lo;
            }

            for(size_t tableIndex = firstIndex; tableIndex < table.size(); ++tableIndex) {
                const uint32_t alignmentId = table[tableIndex];
                const auto& ad = alignmentData[alignmentId];

                // Only consider overlaps that passed phasing (cis on both) and are still eligible.
                if(!ad.isCisByBothSides()) {
                    continue;
                }
                if(!ad.keptByBothSides()) {
                    continue;
                }

                // Mark self-overlaps as deleted immediately.
                if(ad.readIds[0] == ad.readIds[1]) {
                    alignmentData[alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
                    ++localRemovedCount;
                    continue;
                }

                // Canonical ordering: only process from the smaller ReadId side.
                if(ad.readIds[0] != r0) {
                    continue;
                }
                const ReadId target = ad.readIds[1];

                if(target != currentTarget) {
                    flushGroup();
                    group.clear();
                    currentTarget = target;
                }

                const uint64_t overlapLenOnR0 = uint64_t(ad.qe) - uint64_t(ad.qs);
                const int64_t score = ad.info.hifiasmSharedSeedScoreProxy();
                group.push_back(CandidateInfo{alignmentId, score, target, overlapLenOnR0});
            }

            flushGroup();
        }
    }

    removedSecondaryAlignmentCount += localRemovedCount;
}



/**
 * @brief Deduplicate ONT overlap chains to keep one overlap per partner read (hifiasm parity).
 *
 * ## Purpose
 * After lchain+mcopy chaining and base-level alignment, multiple overlap chains may exist
 * for the same (query, partner) read pair. This function implements hifiasm's dedup_chains
 * logic (ecovlp.cpp:2984) to collapse these to a single best overlap per partner.
 *
 * ## Hifiasm Reference: ecovlp.cpp:2984-3030 `dedup_chains(overlap_region_alloc* ol)`
 *
 * ### Hifiasm Algorithm:
 * ```c
 * 1. Sort overlaps by y_id (partner read ID)
 * 2. Group overlaps by y_id
 * 3. For each group with >1 entry:
 *    a. Find best overlap based on:
 *       - Priority 1: Minimize is_match (0=perfect, 1=cis, 2=trans)
 *       - Priority 2: Maximize score = (x_pos_e + 1 - x_pos_s) - 12*non_homopolymer_errors
 *       - Priority 3: Maximize span = x_pos_e + 1 - x_pos_s (tie-breaker)
 *    b. Swap best to position m, increment m
 * 4. Trim ol->length to m (remove non-best in-place)
 * ```
 *
 * ### Dinara Implementation:
 * ```cpp
 * 1. For each query read r0:
 *    a. Process alignments where r0 is readIds[0] (canonical ordering)
 *    b. Group alignments by partner readIds[1]
 *    c. Skip alignments with non-phase deletion reasons
 * 2. For each partner group:
 *    a. Find best alignment based on:
 *       - Priority 1: Minimize isMatchRank (1=cis, 2=trans)
 *       - Priority 2: Maximize score = span - 12*errorCount
 *       - Priority 3: Maximize span (tie-breaker)
 *       - Priority 4: Minimize alignmentId (determinism tie-breaker)
 *    b. Mark non-best as DeleteReasonSecondary (lazy deletion)
 * ```
 *
 * ## Algorithm Equivalence
 *
 * Close to hifiasm's dedup_chains logic, but not byte-for-byte identical.
 *
 * | Aspect | Hifiasm | Dinara | Equivalent? |
 * |--------|---------|--------|-------------|
 * | **Grouping** | By y_id | By readIds[1] (partner) | ✅ Same |
 * | **Priority 1** | Minimize is_match | Minimize stored EC match state | Near-parity |
 * | **Priority 2** | Maximize (span - 12*errors) | Maximize (span - 12*errors) | Approximate |
 * | **Priority 3** | Maximize span | Maximize span | ✅ Same |
 * | **Deletion** | In-place compaction | Mark DeleteReasonSecondary | Behaviorally similar |
 *
 * ### Field Mappings:
 *
 * | Hifiasm Field | Dinara Field | Notes |
 * |---------------|--------------|-------|
 * | `y_id` | `ad.readIds[1]` | Partner read ID |
 * | `x_pos_s, x_pos_e` | `ad.qs, ad.qe` | Query span (when readIds[0] = r0) |
 * | `is_match` | `hifiasmEcMatchState{0,1}` | Stored per read perspective |
 * | error count | `mismatchCount + gapCount` | From the CIGAR walk |
 *
 * ## Why Deduplication is Necessary
 *
 * The lchain+mcopy algorithm (with mcopy_num > 1) can extract multiple high-scoring
 * chains between the same read pair. This is useful during chaining to:
 * 1. Handle repetitive genomic regions (segmental duplications)
 * 2. Capture alternative alignments for ambiguous regions
 * 3. Allow phasing to select the best cis overlap
 *
 * However, after phasing marks cis/trans overlaps, we want only **one** overlap per
 * partner for downstream assembly. This function implements that final deduplication.
 *
 * ## Execution Context
 *
 * This function MUST run after:
 * 1. ✅ `computeBaseAlignmentsAndStore()` - base-level alignment computed
 * 2. ✅ `performHifiasmECParity()` - phasing marks cis/trans (DeleteReasonPhase)
 * 3. ✅ Before final alignment filtering - allows best overlap to survive
 *
 * ## Thread Safety
 *
 * - Parallelized per query read (no conflicts between threads)
 * - Each thread processes disjoint sets of query reads
 * - Atomic counter for removed overlap count
 *
 * @param threadCount Number of threads for parallel processing
 *
 * @complexity O(N*M) where N=read count, M=avg alignments per read
 *             Amortized O(N log M) per read due to binary search for firstIndex
 *
 * @see hifiasm ecovlp.cpp:2984 dedup_chains()
 * @see performHifiasmECParity() for cis/trans phasing
 * @see AlignmentData::DeleteReasonSecondary for lazy deletion marker
 */
void Assembler::deduplicateOntChainsPerPartnerReadHifiasmLike(uint64_t threadCount)
{
    const auto tBegin = steady_clock::now();
    cout << timestamp << "ONT-style dedup (one overlap per partner) begins." << endl;

    removedOntDeduplicatedChainCount = 0;

    const uint64_t readCount = reads->readCount();
    setupLoadBalancing(readCount, 100); // Same heuristic as other per-read filters.
    runThreads(&Assembler::deduplicateOntChainsPerPartnerReadHifiasmLikeThreadFunction, threadCount);

    const auto tEnd = steady_clock::now();
    const double tSeconds = seconds(tEnd - tBegin);
    cout << timestamp << "ONT-style dedup ends in " << tSeconds << " s. Removed "
         << removedOntDeduplicatedChainCount.load() << " redundant overlaps (marked DeleteReasonSecondary)." << endl;
}



/**
 * @brief Worker thread function for ONT chain deduplication per query read.
 *
 * This implements the per-read logic equivalent to hifiasm's dedup_chains.
 * Each thread processes a disjoint subset of query reads.
 *
 * @see deduplicateOntChainsPerPartnerReadHifiasmLike() for algorithm overview
 */
void Assembler::deduplicateOntChainsPerPartnerReadHifiasmLikeThreadFunction(size_t)
{
    uint64_t begin = 0, end = 0;
    uint64_t localRemoved = 0;

    // =========================================================================
    // Per-Alignment Metadata for Group Selection
    // =========================================================================
    struct CandidateInfo {
        uint32_t alignmentId = 0;
        ReadId partner = invalidReadId;
        uint32_t errorCount = 0; // hifiasm-style non-homopolymer error count
        uint32_t span = 0;       // qe - qs (query span)
        uint8_t isMatchRank = 1; // Hifiasm-style: 1=cis, 2=trans
        int64_t score = std::numeric_limits<int64_t>::min(); // span - 12*errorCount
    };
    vector<CandidateInfo> group;
    group.reserve(8);  // Typical case: few overlaps per partner

    // =========================================================================
    // Lambda: Select Best Overlap Per Partner Group (Hifiasm Logic)
    // =========================================================================
    // Implements hifiasm ecovlp.cpp:2994-3014 selection logic:
    //   - Priority 1: Minimize stored hifiasm-style match state
    //   - Priority 2: Maximize score (span - 12*errors)
    //   - Priority 3: Maximize span (tie-breaker)
    //   - Priority 4: Minimize alignmentId (Dinara's determinism tie-breaker)
    auto flushGroup = [&]() {
        if(group.size() <= 1) {
            return;  // No deduplication needed for single overlap
        }

        // -----------------------------------------------------------------
        // Find Best Overlap in Group (Hifiasm Priority Order)
        // -----------------------------------------------------------------
        size_t best = 0;
        for(size_t i = 1; i < group.size(); ++i) {
            const auto& a = group[i];
            const auto& b = group[best];

            // Priority 1: Prefer smaller hifiasm-style match state.
            if(a.isMatchRank != b.isMatchRank) {
                if(a.isMatchRank < b.isMatchRank) best = i;
                continue;
            }

            // Priority 2: Maximize score = span - 12*errors
            // Matches hifiasm: sc = plus - minus; if(sc > mm_sc) sf = 1;
            if(a.score != b.score) {
                if(a.score > b.score) best = i;
                continue;
            }

            // Priority 3: Maximize span (tie-breaker)
            // Matches hifiasm: if((sc == mm_sc) && ((z->x_pos_e+1-z->x_pos_s) > ...)) sf = 1;
            if(a.span != b.span) {
                if(a.span > b.span) best = i;
                continue;
            }

            // Priority 4: Minimize alignmentId (Dinara determinism tie-breaker;
            // hifiasm leaves perfect ties iteration-order dependent).
            if(a.alignmentId < b.alignmentId) {
                best = i;
            }
        }

        // -----------------------------------------------------------------
        // Mark Non-Best as Secondary (Lazy Deletion)
        // -----------------------------------------------------------------
        // Hifiasm swaps best to position m and trims array (in-place deletion).
        // Dinara marks as DeleteReasonSecondary (lazy deletion, removed later).
        for(size_t i = 0; i < group.size(); ++i) {
            if(i == best) continue;  // Keep best overlap
            alignmentData[group[i].alignmentId].addDeleteReasonsBoth(AlignmentData::DeleteReasonSecondary);
            ++localRemoved;
        }
    };

    // =========================================================================
    // Main Loop: Process Each Query Read in Thread's Batch
    // =========================================================================
    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); ++r0) {
            const OrientedReadId orientedR0(r0, 0);
            const auto& table = alignmentTable[orientedR0.getValue()];
            if(table.empty()) {
                continue;  // No alignments for this read
            }

            group.clear();
            ReadId currentPartner = invalidReadId;

            // -----------------------------------------------------------------
            // Binary Search: Skip Prefix Where Partner < r0 (Avoid Duplicates)
            // -----------------------------------------------------------------
            // Dinara stores alignments bidirectionally (both (r0, r1) and (r1, r0)).
            // To avoid processing the same pair twice, we only process pairs where
            // r0 < partner (canonical ordering). Binary search finds the first such entry.
            size_t firstIndex = 0;
            {
                size_t lo = 0;
                size_t hi = table.size();
                while(lo < hi) {
                    const size_t mid = (lo + hi) / 2;
                    const uint32_t alignmentId = table[mid];
                    const OrientedReadId other = alignmentData[alignmentId].getOther(orientedR0);
                    if(other < orientedR0) {
                        lo = mid + 1;  // Partner < r0, skip
                    } else {
                        hi = mid;      // Partner >= r0, could be first
                    }
                }
                firstIndex = lo;
            }

            // -----------------------------------------------------------------
            // Group Alignments by Partner (Hifiasm Groups by y_id)
            // -----------------------------------------------------------------
            for(size_t tableIndex = firstIndex; tableIndex < table.size(); ++tableIndex) {
                const uint32_t alignmentId = table[tableIndex];
                const auto& ad = alignmentData[alignmentId];

                // Canonical ordering: only process groups where r0 is readIds[0]
                // so ad.qs/qe represent span on this query read.
                if(ad.readIds[0] != r0) {
                    continue;
                }
                const ReadId partner = ad.readIds[1];

                // -----------------------------------------------------------------
                // Filter: Skip Overlaps with Non-Phase Deletion Reasons
                // -----------------------------------------------------------------
                // Match hifiasm's ONT behavior: dedup runs after phasing marks
                // cis/trans (is_match), but before overlaps are hard-removed.
                // Therefore we allow phase-deleted (trans) overlaps to participate,
                // but skip overlaps with OTHER deletion reasons (e.g., low quality).
                const auto nonPhase0 = ad.deleteReasons0 & ~AlignmentData::DeleteReasonPhase;
                const auto nonPhase1 = ad.deleteReasons1 & ~AlignmentData::DeleteReasonPhase;
                if(nonPhase0 != AlignmentData::DeleteReasonNone ||
                   nonPhase1 != AlignmentData::DeleteReasonNone) {
                    continue;  // Already filtered, skip
                }

                // -----------------------------------------------------------------
                // New Partner Group: Flush Previous Group
                // -----------------------------------------------------------------
                if(partner != currentPartner) {
                    flushGroup();  // Deduplicate previous partner's group
                    group.clear();
                    currentPartner = partner;
                }

                // -----------------------------------------------------------------
                // Compute Overlap Metrics (Hifiasm Fields)
                // -----------------------------------------------------------------
                const uint32_t span = ad.qe - ad.qs;  // Matches hifiasm: x_pos_e + 1 - x_pos_s

                const uint32_t err =
                    ((ad.info.mismatchCount == invalid<uint32_t>) ? 0U : ad.info.mismatchCount) +
                    ((ad.info.gapCount == invalid<uint32_t>) ? 0U : ad.info.gapCount);

                // Hifiasm-style EC match state: 1=cis, 2=trans.
                // Lower value (cis) is preferred in dedup.
                const uint8_t isMatchRank = ad.getHifiasmEcMatchStateFromReadPerspective(r0);

                // Score formula: span - 12*errors
                // Matches hifiasm: sc = plus - minus (ecovlp.cpp:2997)
                const int64_t score = int64_t(span) - 12 * int64_t(err);

                group.push_back(CandidateInfo{alignmentId, partner, err, span, isMatchRank, score});
            }

            // Flush final partner group for this query read
            flushGroup();
        }
    }

    // Update global removed count (atomic)
    removedOntDeduplicatedChainCount += localRemoved;
}



// For each read pair with multiple chains to the same target on the same
// strand, keep the single best chain (highest hifiasm shared_seed) and delete
// the rest. Mirrors hifiasm's live special_lchain selection; see
// removeMultiChainAlignmentsThreadFunction for details.
void Assembler::removeMultiChainAlignments(uint64_t threadCount)
{
    const auto tBegin = steady_clock::now();
    cout << timestamp << "Removing multi-chain alignments begins." << endl;

    removedMultiChainCount = 0;

    const uint64_t readCount = reads->readCount();
    setupLoadBalancing(readCount, 100);
    runThreads(&Assembler::removeMultiChainAlignmentsThreadFunction, threadCount);

    const auto tEnd = steady_clock::now();
    const double tSeconds = seconds(tEnd - tBegin);
    cout << timestamp << "Removing multi-chain alignments ends in " << tSeconds << " s. Removed "
         << removedMultiChainCount.load() << " alignments from multi-chain pairs." << endl;
}


void Assembler::removeMultiChainAlignmentsThreadFunction(size_t)
{
    uint64_t begin = 0, end = 0;
    uint64_t localRemoved = 0;

    // Candidate chains between r0 and one target, plus the key hifiasm's live
    // selection uses to pick the survivor.
    struct Candidate {
        uint32_t alignmentId;
        int32_t sharedSeed;  // hifiasm overlap_region.shared_seed (chain DP score)
    };
    // Grouped by (target read, strand). hifiasm keeps same-strand (paf[]) and
    // reverse-strand (reverse_paf[]) overlaps in separate stores that never
    // compete, so an inverted-repeat pair legitimately keeps one chain per
    // orientation. Index [0]=same-strand, [1]=reverse-strand.
    std::map<ReadId, array<vector<Candidate>, 2>> byTarget;

    while(getNextBatch(begin, end)) {
        for(ReadId r0 = ReadId(begin); r0 != ReadId(end); ++r0) {
            const OrientedReadId orientedR0(r0, 0);
            const auto& table = alignmentTable[orientedR0.getValue()];
            if(table.empty()) continue;

            byTarget.clear();

            for(size_t tableIndex = 0; tableIndex < table.size(); ++tableIndex) {
                const uint32_t alignmentId = table[tableIndex];
                const auto& ad = alignmentData[alignmentId];

                // Only process each pair once: r0 must be readIds[0].
                if(ad.readIds[0] != r0) continue;

                // Skip self-overlaps.
                if(ad.readIds[0] == ad.readIds[1]) continue;

                // Skip already-deleted alignments.
                if(ad.deleteReasons0 || ad.deleteReasons1) continue;

                // hifiasm's live selection among multiple chains to the same
                // target (special_lchain, anchor.cpp): keep the chain with the
                // highest overlap_region.shared_seed (the minimizer-chain DP
                // score). This is the key hifiasm actually uses; the
                // span - 12*non_homopolymer_errors formula lives only in a
                // disabled (/**...**/) block in ecovlp.cpp and is NOT run.
                const int32_t sharedSeed =
                    (ad.info.sharedSeedScore != invalid<int32_t>) ?
                    ad.info.sharedSeedScore : 0;

                const int strandIdx = ad.isSameStrand ? 0 : 1;
                byTarget[ad.readIds[1]][strandIdx].push_back(
                    Candidate{alignmentId, sharedSeed});
            }

            // For each (target, strand) with >1 chain, keep the highest
            // shared_seed and delete the rest. A legitimate overlap produces one
            // chain per target per strand; multiple chains indicate a repeat
            // region where the chainer found several plausible paths. hifiasm
            // keeps the best-scoring one rather than discarding the pair,
            // preserving coverage. Ties are broken by smaller alignmentId for
            // determinism (hifiasm's internal x_pos_strand tie-break has no
            // dinara equivalent).
            for(auto& [targetId, strandGroups] : byTarget) {
                for(int s = 0; s < 2; ++s) {
                    auto& candidates = strandGroups[s];
                    if(candidates.size() <= 1) continue;

                    size_t best = 0;
                    for(size_t i = 1; i < candidates.size(); ++i) {
                        const auto& c = candidates[i];
                        const auto& b = candidates[best];
                        if(c.sharedSeed > b.sharedSeed ||
                           (c.sharedSeed == b.sharedSeed &&
                            c.alignmentId < b.alignmentId)) {
                            best = i;
                        }
                    }

                    for(size_t i = 0; i < candidates.size(); ++i) {
                        if(i == best) continue;
                        alignmentData[candidates[i].alignmentId].addDeleteReasonsBoth(
                            AlignmentData::DeleteReasonSecondary);
                        ++localRemoved;
                    }
                }
            }
        }
    }

    removedMultiChainCount += localRemoved;
}


