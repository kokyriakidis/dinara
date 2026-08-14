#include "ProjectedAlignment.hpp"
#include "Alignment.hpp"
#include "Assembler.hpp"
#include "Base.hpp"
#include "LongBaseSequence.hpp"
#include "Marker.hpp"
#include "OverlapCigarStore.hpp"
#include "Reads.hpp"
using namespace dinara;

#include "algorithm.hpp"
#include <iostream>




ProjectedAlignment::ProjectedAlignment(
    const Assembler& assembler,
    const array<OrientedReadId, 2>& orientedReadIds,
    const Alignment& alignment,
    Method method,
    int64_t dpMatchScore,
    int64_t dpMismatchScore,
    int64_t dpGapOpen1,
    int64_t dpGapExtend1) :

    ProjectedAlignment(
        uint32_t(assembler.assemblerInfo->k),
        orientedReadIds,
        {
            assembler.getReads().getRead(orientedReadIds[0].getReadId()),
            assembler.getReads().getRead(orientedReadIds[1].getReadId())
        },
        alignment,
        {
            (*assembler.markers)[orientedReadIds[0].getValue()],
            (*assembler.markers)[orientedReadIds[1].getValue()]
        },
        method,
        dpMatchScore,
        dpMismatchScore,
        dpGapOpen1,
        dpGapExtend1)
{
}



ProjectedAlignment::ProjectedAlignment(
    uint32_t k,
    const array<OrientedReadId, 2>& orientedReadIds,
    const array<LongBaseSequenceView, 2>& sequences,
    const Alignment& alignment,
    const array< span<const CompressedMarker>, 2>& markers,
    Method method,
    int64_t dpMatchScore,
    int64_t dpMismatchScore,
    int64_t dpGapOpen1,
    int64_t dpGapExtend1,
    OverlapCigarStore* cigarStore) :
    k(k),
    kHalf(k / 2),
    dpMatchScore(dpMatchScore),
    dpMismatchScore(dpMismatchScore),
    dpGapOpen1(dpGapOpen1),
    dpGapExtend1(dpGapExtend1),
    orientedReadIds(orientedReadIds),
    sequences(sequences),
    alignment(alignment),
    markers(markers),
    cigarStore(cigarStore)
{
    DINARA_ASSERT((k % 2) == 0);

    switch(method) {
        case Method::All:
            constructAll();
            break;
        case Method::QuickRaw:
            constructQuickRaw();
            break;
        case Method::QuickRawSparse:
            constructQuickRawSparse();
            break;
        case Method::None:
            // Caller fills the object explicitly (constructFromHifiasmCigar).
            break;
        default:
            DINARA_ASSERT(0);
    }
}



void ProjectedAlignment::constructAll()
{
    mismatchCountRle = 0;

    // Loop over pairs of consecutive aligned markers (A, B).
    for(uint64_t iB=1; iB<alignment.ordinals.size(); iB++) {
        const uint64_t iA = iB - 1;

        // Get the ordinals of these pair of consecutive aligned markers.
        const array<uint32_t, 2>& ordinalsA = alignment.ordinals[iA];
        const array<uint32_t, 2>& ordinalsB = alignment.ordinals[iB];

        segments.push_back(ProjectedAlignmentSegment(
            kHalf,
            ordinalsA,
            ordinalsB,
            markers));
        ProjectedAlignmentSegment& segment = segments.back();

        // Fill in the base sequences.
        fillSequences(segment);

        // Align them.
        segment.computeAlignment(
            dpMatchScore, dpMismatchScore,
            dpGapOpen1, dpGapExtend1);

        // Same, in RLE.
        segment.fillRleSequences();
        segment.computeRleAlignment();
        mismatchCountRle += segment.mismatchCountRle;
    }



    computeStatistics();
}



// This stores only the following:
// - RLE sequences and RLE alignments for segments for which the RLE sequences
//   of the two oriented reads are different.
// - Total RLE edit distance and total RLE lengths.
// This stores only the raw sequences and alignments for segments for which the raw sequences
// of the two oriented reads are different. Even though it only stores segments with mismatches, 
// the computed total length includes all segments.
void ProjectedAlignment::constructQuickRaw()
{
    // Create the segment outside the loop and reuse it to reduce
    // memory allocation activity.
    ProjectedAlignmentSegment segment;

    // Initialize statistics
    totalLength = {0, 0};
    totalEditDistance = 0;
    mismatchCount = 0;
    totalIndelBaseCount = 0;
    totalGapEventCount = 0;
    totalDpScore = 0;
    hasLargeIndel = false;

    // Hifiasm current overlap scoring uses single-affine parameters.
    // For identical segments we only need the match reward.
    const int64_t dpMatch = dpMatchScore;

    // Loop over pairs of consecutive aligned markers (A, B).
    for(uint64_t iB=1; iB<alignment.ordinals.size(); iB++) {
        const uint64_t iA = iB - 1;

        // Store in the segment the ordinals of these pair of consecutive aligned markers.
        segment.ordinalsA = alignment.ordinals[iA];
        segment.ordinalsB = alignment.ordinals[iB];

        // Store the corresponding positions.
        for(uint64_t i=0; i<2; i++) {
            segment.positionsA[i] = markers[i][segment.ordinalsA[i]].position + kHalf;
            segment.positionsB[i] = markers[i][segment.ordinalsB[i]].position + kHalf;
        }

        // Fill in the base sequences.
        fillSequences(segment);

        // Accumulate total lengths (even for identical sequences)
        for(uint64_t i=0; i<2; i++) {
            totalLength[i] += segment.sequences[i].size();
        }

        // If the raw sequences are the same, don't store the segment.
        if(segment.sequences[0] == segment.sequences[1]) {
            totalDpScore += dpMatch * int64_t(segment.sequences[0].size());
            continue;
        }

        // Align them.
        segment.computeAlignment(
            dpMatchScore, dpMismatchScore,
            dpGapOpen1, dpGapExtend1);

        // Accumulate statistics
        totalEditDistance += segment.editDistance;
        mismatchCount += segment.mismatchCount;
        totalIndelBaseCount += segment.indelBaseCount;
        totalGapEventCount += segment.gapEventCount;
        totalDpScore += segment.dpScore;
        if (segment.hasLargeIndel) hasLargeIndel = true;

        // Store the segment.
        segments.push_back(segment);
    }
}



// Similar to constructQuickRaw, but does not store per-base alignment traces.
// Instead, it stores sparse mismatch and indel events that are sufficient to
// populate AlignedEvidenceStore efficiently.
//
// When cigarStore is non-null, also packs the full per-overlap CIGAR into
// the store using hifiasm-style uint16_t tokens (2-bit op + 14-bit length).
void ProjectedAlignment::constructQuickRawSparse()
{
    sparseMismatches.clear();
    sparseIndels.clear();

    // Initialize statistics.
    totalLength = {0, 0};
    totalEditDistance = 0;
    mismatchCount = 0;
    totalIndelBaseCount = 0;
    totalGapEventCount = 0;
    totalDpScore = 0;
    hasLargeIndel = false;

    // Begin a new CIGAR entry in the store if requested.
    const bool storeCigar = (cigarStore != nullptr);
    if(storeCigar) {
        cigarOffset = cigarStore->beginAlignment();
    }

    // Record the CIGAR's boundary positions (first and last marker pairs).
    // These are in oriented-read coordinates.
    if(alignment.ordinals.size() >= 2) {
        const auto lastIdx = alignment.ordinals.size() - 1;
        cigarRead0Start = markers[0][alignment.ordinals[0][0]].position + kHalf;
        cigarRead1Start = markers[1][alignment.ordinals[0][1]].position + kHalf;
        cigarRead0End   = markers[0][alignment.ordinals[lastIdx][0]].position + kHalf;
        cigarRead1End   = markers[1][alignment.ordinals[lastIdx][1]].position + kHalf;
    }

    // Scoring parameters (hoisted out of the per-segment loop).
    const int64_t match = dpMatchScore;
    const int64_t mismatch = dpMismatchScore;
    auto gapPenalty = [gapOpen1 = dpGapOpen1, gapExtend1 = dpGapExtend1](uint64_t length) -> int64_t {
        return gapOpen1 + gapExtend1 * int64_t(length);
    };

    // Thread-local buffers for ASCII sequences.
    static thread_local vector<uint8_t> asciiSequence0;
    static thread_local vector<uint8_t> asciiSequence1;
    static constexpr uint8_t baseToAscii[4] = {'A', 'C', 'G', 'T'};
    static constexpr auto asciiToBase = []() constexpr {
        array<uint8_t, 256> m{};
        m['A'] = 0; m['C'] = 1; m['G'] = 2; m['T'] = 3;
        return m;
    }();

    // Loop over pairs of consecutive aligned markers (A, B).
    for(uint64_t iB=1; iB<alignment.ordinals.size(); iB++) {
        const uint64_t iA = iB - 1;

        // Compute segment positions directly (no ProjectedAlignmentSegment needed).
        const uint32_t begin0 = markers[0][alignment.ordinals[iA][0]].position + kHalf;
        const uint32_t end0   = markers[0][alignment.ordinals[iB][0]].position + kHalf;
        const uint32_t begin1 = markers[1][alignment.ordinals[iA][1]].position + kHalf;
        const uint32_t end1   = markers[1][alignment.ordinals[iB][1]].position + kHalf;
        const uint32_t len0 = end0 - begin0;
        const uint32_t len1 = end1 - begin1;

        totalLength[0] += len0;
        totalLength[1] += len1;

        // Empty segments: emit indel for CIGAR consistency if needed.
        // read0 is the query, read1 is the target (SAM/PAF convention):
        // query-only bases are an insertion, target-only bases a deletion.
        if(len0 == 0 || len1 == 0) {
            if(storeCigar) {
                if(len0 > 0) cigarStore->pushInsertion(len0);
                else if(len1 > 0) cigarStore->pushDeletion(len1);
            }
            continue;
        }

        // Fill ASCII sequences directly from oriented reads.
        asciiSequence0.resize(len0);
        asciiSequence1.resize(len1);
        for(uint32_t j = 0; j < len0; ++j) {
            asciiSequence0[j] = baseToAscii[getBase(0, begin0 + j).value];
        }
        for(uint32_t j = 0; j < len1; ++j) {
            asciiSequence1[j] = baseToAscii[getBase(1, begin1 + j).value];
        }

        // Identical segments: fast path.
        if(len0 == len1 && memcmp(asciiSequence0.data(), asciiSequence1.data(), len0) == 0) {
            totalDpScore += match * int64_t(len0);
            if(storeCigar) cigarStore->pushMatch(len0);
            continue;
        }

        // Align with A*PA2.
        char* cigar = nullptr;
        size_t cigarLen = 0;
        const int64_t cost = astarpa2_simple(
            asciiSequence0.data(), asciiSequence0.size(),
            asciiSequence1.data(), asciiSequence1.size(),
            (unsigned char**)&cigar, &cigarLen);

        totalEditDistance += cost;

        uint64_t position0 = 0;
        uint64_t position1 = 0;
        size_t currentVal = 0;
        for(size_t i=0; i<cigarLen; i++) {
            const char c = cigar[i];
            if(c >= '0' && c <= '9') {
                currentVal = currentVal * 10 + size_t(c - '0');
                continue;
            }
            if(currentVal == 0) currentVal = 1;

            if(c == 'M' || c == '=' || c == 'X') {
                uint64_t mismatchHere = 0;
                uint32_t runStart = 0;
                bool runIsMismatch = false;

                for(size_t k=0; k<currentVal; ++k) {
                    const uint8_t a0 = asciiSequence0[position0 + k];
                    const uint8_t a1 = asciiSequence1[position1 + k];
                    const bool isMismatch = (a0 != a1);

                    if(isMismatch) {
                        sparseMismatches.push_back(ProjectedAlignmentSparseMismatch{
                            uint32_t(begin0 + uint32_t(position0 + k)),
                            uint32_t(begin1 + uint32_t(position1 + k)),
                            asciiToBase[a0],
                            asciiToBase[a1]
                        });
                        ++mismatchHere;
                    }

                    // Accumulate match/mismatch runs for CIGAR tokens.
                    if(storeCigar) {
                        if(k == 0) {
                            runIsMismatch = isMismatch;
                            runStart = 0;
                        } else if(isMismatch != runIsMismatch) {
                            cigarStore->pushOp(runIsMismatch ? 1 : 0, uint32_t(k) - runStart);
                            runIsMismatch = isMismatch;
                            runStart = uint32_t(k);
                        }
                    }
                }

                if(storeCigar) {
                    cigarStore->pushOp(runIsMismatch ? 1 : 0, uint32_t(currentVal) - runStart);
                }

                mismatchCount += mismatchHere;
                totalDpScore += match * int64_t(currentVal - mismatchHere) + mismatch * int64_t(mismatchHere);
                position0 += currentVal;
                position1 += currentVal;

            } else if(c == 'D') {
                sparseIndels.push_back(ProjectedAlignmentSparseIndel{
                    uint32_t(begin0 + uint32_t(position0)),
                    uint32_t(begin1 + uint32_t(position1)),
                    uint32_t(currentVal), 'D'
                });
                // A*PA2 'D' consumes read0/query, which in SAM/PAF terms is an
                // insertion (query bases absent from the target).
                if(storeCigar) cigarStore->pushInsertion(uint32_t(currentVal));
                position0 += currentVal;
                totalIndelBaseCount += currentVal;
                ++totalGapEventCount;
                totalDpScore -= gapPenalty(uint64_t(currentVal));
                if(currentVal >= 6) hasLargeIndel = true;

            } else if(c == 'I') {
                sparseIndels.push_back(ProjectedAlignmentSparseIndel{
                    uint32_t(begin0 + uint32_t(position0)),
                    uint32_t(begin1 + uint32_t(position1)),
                    uint32_t(currentVal), 'I'
                });
                // A*PA2 'I' consumes read1/target, which in SAM/PAF terms is a
                // deletion (target bases absent from the query).
                if(storeCigar) cigarStore->pushDeletion(uint32_t(currentVal));
                position1 += currentVal;
                totalIndelBaseCount += currentVal;
                ++totalGapEventCount;
                totalDpScore -= gapPenalty(uint64_t(currentVal));
                if(currentVal >= 6) hasLargeIndel = true;

            } else if(c == 'H') {
                // Hard clips: tolerate but ignore.
            } else {
                DINARA_ASSERT(0);
            }

            currentVal = 0;
        }

        astarpa_free_cigar((unsigned char*)cigar);

        DINARA_ASSERT(position0 == asciiSequence0.size());
        DINARA_ASSERT(position1 == asciiSequence1.size());
    }

    // Finalize CIGAR and verify consistency.
    if(storeCigar && cigarOffset != uint32_t(-1)) {
        cigarTokenCount = cigarStore->tokensSince(cigarOffset);
        uint64_t cigarConsumed0 = 0;
        uint64_t cigarConsumed1 = 0;
        cigarStore->forEachOp(cigarOffset, cigarTokenCount, [&](uint8_t op, uint32_t len) {
            // read0 is the query, read1 is the target. Decide consumption
            // through the shared helpers so this check follows the single
            // op convention definition and cannot drift from the producer.
            if(opConsumesQuery(op))  cigarConsumed0 += len;  // read0
            if(opConsumesTarget(op)) cigarConsumed1 += len;  // read1
        });
        DINARA_ASSERT(cigarConsumed0 == totalLength[0]);
        DINARA_ASSERT(cigarConsumed1 == totalLength[1]);
    }
}



bool ProjectedAlignment::constructFromHifiasmCigar(
    span<const CigarToken> normalizedTokens,
    uint32_t read0Anchor,
    uint32_t read1Anchor,
    uint32_t read0Begin,
    uint32_t read0End)
{
    // The hifiasm CIGAR must cover the marker interval on read0; otherwise the
    // caller should fall back to recomputation.
    if(normalizedTokens.empty()) return false;
    if(read0Begin < read0Anchor || read0End < read0Begin) return false;
    {
        // read0 span covered by the tokens = sum of query-consuming ops.
        uint64_t read0Covered = 0;
        for(const CigarToken& t : normalizedTokens) {
            if(opConsumesQuery(t.op())) read0Covered += t.len();
        }
        if(uint64_t(read0End) > uint64_t(read0Anchor) + read0Covered) return false;
    }

    // hifiasm's exported CIGAR labels every column '='(0)/'X'(1)/'I'(2)/'D'(3)
    // per SAM convention, and the aligner emits 'X' for each mismatched column
    // (never a generic 'M' hiding a mismatch inside a match run). So we trust the
    // op labels for every statistic and do NOT re-read bases: no per-column
    // comparison. The sparse per-base diffs are likewise not produced here --
    // the only consumer of ProjectedAlignment::sparse* (the leaf-snarl phasing
    // path) builds its own QuickRawSparse alignment, never this one.
    totalLength = {0, 0};
    totalEditDistance = 0;
    mismatchCount = 0;
    totalIndelBaseCount = 0;
    totalGapEventCount = 0;
    totalDpScore = 0;
    hasLargeIndel = false;

    const bool storeCigar = (cigarStore != nullptr);
    if(storeCigar) {
        cigarOffset = cigarStore->beginAlignment();
    }

    cigarRead0Start = read0Begin;
    cigarRead0End   = read0End;
    // read1 boundaries are resolved from the walk (the read1 position aligned to
    // read0Begin / read0End), mirroring how constructQuickRawSparse derives them
    // from the first/last aligned marker pair.
    cigarRead1Start = read1Anchor;
    cigarRead1End   = read1Anchor;
    bool haveRead1Start = false;

    const int64_t match = dpMatchScore;
    const int64_t mismatch = dpMismatchScore;
    auto gapPenalty = [gapOpen1 = dpGapOpen1, gapExtend1 = dpGapExtend1](uint64_t length) -> int64_t {
        return gapOpen1 + gapExtend1 * int64_t(length);
    };

    // Coalesce runs and walk clipped to [read0Begin, read0End) on read0,
    // replicating OverlapCigarStore::walkRange. op match(0)/mismatch(1) consume
    // both reads; ins(2) consumes read0 (query-only); del(3) consumes read1.
    // Statistics come straight from the op labels (see top-of-function note).
    uint64_t xk = read0Anchor;   // read0 forward position
    uint64_t yk = read1Anchor;   // read1 alignment-orientation position

    // Track read1 boundary aligned to read0End.
    auto noteBoundaries = [&](uint64_t x, uint64_t y) {
        if(!haveRead1Start && x >= read0Begin) {
            cigarRead1Start = uint32_t(y);
            haveRead1Start = true;
        }
    };

    size_t i = 0;
    while(i < normalizedTokens.size()) {
        const uint8_t op = normalizedTokens[i].op();
        uint32_t totalLen = normalizedTokens[i].len();
        size_t peek = i + 1;
        while(peek < normalizedTokens.size() && normalizedTokens[peek].op() == op) {
            totalLen += normalizedTokens[peek].len();
            ++peek;
        }
        i = peek;

        const uint64_t xkEnd = xk + (opConsumesQuery(op)  ? totalLen : 0);
        const uint64_t ykEnd = yk + (opConsumesTarget(op) ? totalLen : 0);

        if(!opConsumesQuery(op)) {
            // Deletion (read1-only). Emitted iff its read0 anchor is in range.
            if(xk >= read0Begin && xk < read0End) {
                noteBoundaries(xk, yk);
                if(storeCigar) cigarStore->pushDeletion(totalLen);
                totalLength[1] += totalLen;
                totalIndelBaseCount += totalLen;
                ++totalGapEventCount;
                totalEditDistance += int64_t(totalLen);
                totalDpScore -= gapPenalty(uint64_t(totalLen));
                if(totalLen >= 6) hasLargeIndel = true;
                cigarRead1End = uint32_t(ykEnd);
            }
        } else if(xkEnd > read0Begin && xk < read0End) {
            // Match/mismatch or insertion overlapping the range: clip to it.
            const uint64_t clipStart = (xk < read0Begin) ? uint64_t(read0Begin) : xk;
            const uint64_t clipEnd   = (xkEnd > read0End) ? uint64_t(read0End) : xkEnd;
            const uint32_t clipLen   = uint32_t(clipEnd - clipStart);
            const uint64_t skipBases = clipStart - xk;
            const uint64_t adjYk = yk + (opConsumesTarget(op) ? skipBases : 0);
            noteBoundaries(clipStart, adjYk);

            if(op == CigarOpIns) {
                // Insertion (read0-only). SAM/PAF 'I'.
                if(storeCigar) cigarStore->pushInsertion(clipLen);
                totalLength[0] += clipLen;
                totalIndelBaseCount += clipLen;
                ++totalGapEventCount;
                totalEditDistance += int64_t(clipLen);
                totalDpScore -= gapPenalty(uint64_t(clipLen));
                if(clipLen >= 6) hasLargeIndel = true;
                cigarRead1End = uint32_t(adjYk);
            } else {
                // Aligned run: op is match(0) or mismatch(1), uniform across the
                // whole run after coalescing. Trust the op label -- a mismatch run
                // is clipLen mismatched columns, a match run is clipLen identical
                // columns -- so no base fetch or per-column comparison is needed.
                const bool runIsMismatch = (op == CigarOpMismatch);
                if(storeCigar && clipLen > 0) {
                    cigarStore->pushOp(runIsMismatch ? 1 : 0, clipLen);
                }
                if(runIsMismatch) {
                    mismatchCount += clipLen;
                    totalEditDistance += int64_t(clipLen);
                    totalDpScore += mismatch * int64_t(clipLen);
                } else {
                    totalDpScore += match * int64_t(clipLen);
                }
                totalLength[0] += clipLen;
                totalLength[1] += clipLen;
                cigarRead1End = uint32_t(adjYk + clipLen);
            }
        }

        if(xkEnd >= read0End && opConsumesQuery(op)) break;
        xk = xkEnd;
        yk = ykEnd;
    }

    if(storeCigar && cigarOffset != uint32_t(-1)) {
        cigarTokenCount = cigarStore->tokensSince(cigarOffset);
        uint64_t consumed0 = 0, consumed1 = 0;
        cigarStore->forEachOp(cigarOffset, cigarTokenCount, [&](uint8_t op, uint32_t len) {
            if(opConsumesQuery(op))  consumed0 += len;
            if(opConsumesTarget(op)) consumed1 += len;
        });
        DINARA_ASSERT(consumed0 == totalLength[0]);
        DINARA_ASSERT(consumed1 == totalLength[1]);
    }

    return true;
}



ProjectedAlignmentSegment::ProjectedAlignmentSegment(
    uint32_t kHalf,
    const array<uint32_t, 2>& ordinalsA,
    const array<uint32_t, 2>& ordinalsB,
    const array< span<const CompressedMarker>, 2>& markers) :
    ordinalsA(ordinalsA),
    ordinalsB(ordinalsB)
{
    for(uint64_t i=0; i<2; i++) {
        positionsA[i] = markers[i][ordinalsA[i]].position + kHalf;
        positionsB[i] = markers[i][ordinalsB[i]].position + kHalf;
    }
}



void ProjectedAlignmentSegment::computeAlignment(
    int64_t dpMatchScore,
    int64_t dpMismatchScore,
    int64_t dpGapOpen1,
    int64_t dpGapExtend1)
{
    const vector<uint8_t>& sequence0 = reinterpret_cast< const vector<uint8_t>& >(sequences[0]);
    const vector<uint8_t>& sequence1 = reinterpret_cast< const vector<uint8_t>& >(sequences[1]);

    if(sequence0 == sequence1) {
        editDistance = 0;
        alignment.resize(sequence0.size());
        fill(alignment.begin(), alignment.end(), make_pair(true, true));
        mismatchCount = 0;
        indelBaseCount = 0;
        gapEventCount = 0;
        // Hifiasm/minimap2-style match reward (HiFi defaults).
        dpScore = dpMatchScore * int64_t(sequence0.size());
        hasLargeIndel = false;
        maxIndelSize = 0;

    } else {
        // Convert sequences to ASCII for A*PA2
        // Use thread_local buffers to avoid allocation
        static thread_local vector<uint8_t> asciiSequence0;
        static thread_local vector<uint8_t> asciiSequence1;
        asciiSequence0.clear();
        asciiSequence1.clear();
        asciiSequence0.reserve(sequence0.size());
        asciiSequence1.reserve(sequence1.size());

        static const char baseToAscii[] = {'A', 'C', 'G', 'T'};

        // Fast conversion using lookup table
        for(uint8_t b : sequence0) {
            asciiSequence0.push_back(b < 4 ? baseToAscii[b] : 'N');
        }
        for(uint8_t b : sequence1) {
            asciiSequence1.push_back(b < 4 ? baseToAscii[b] : 'N');
        }

        char* cigar = nullptr;
        size_t cigarLen = 0;
        int64_t cost = astarpa2_simple(
            asciiSequence0.data(), asciiSequence0.size(),
            asciiSequence1.data(), asciiSequence1.size(),
            (unsigned char**)&cigar, &cigarLen);
        
        // Convert cost to edit distance (A*PA2 cost is usually edit distance if using defaults)
        // But here we just use the returned cost.
        editDistance = cost;

        // Parse CIGAR directly from char*
        alignment.clear();
        
        // Optimization: Pre-calculate total alignment length to avoid reallocations
        size_t totalAlignmentLength = 0;
        size_t currentVal = 0;
        for(size_t i=0; i<cigarLen; i++) {
            char c = cigar[i];
            if(isdigit(c)) {
                currentVal = currentVal * 10 + (c - '0');
            } else {
                if(currentVal == 0) currentVal = 1;
                // All supported operations (M, =, X, I, D, S, N) contribute to alignment vector length
                // H (hard clip) does not.
                if (c != 'H') {
                    totalAlignmentLength += currentVal;
                }
                currentVal = 0;
            }
        }
        
        alignment.resize(totalAlignmentLength);

        // Fill the alignment vector AND compute mismatchCount
        size_t currentIndex = 0;
        currentVal = 0;
        uint64_t position0 = 0;
        uint64_t position1 = 0;
        mismatchCount = 0;
        indelBaseCount = 0;
        gapEventCount = 0;
        dpScore = 0;
        hasLargeIndel = false;
        maxIndelSize = 0;

        const int64_t match = dpMatchScore;
        const int64_t mismatch = dpMismatchScore;
        auto gapPenalty = [gapOpen1 = dpGapOpen1, gapExtend1 = dpGapExtend1](uint64_t length) -> int64_t {
            // Single-affine gap penalty: O + k*E.
            return gapOpen1 + gapExtend1 * int64_t(length);
        };

        for(size_t i=0; i<cigarLen; i++) {
            char c = cigar[i];
            if(isdigit(c)) {
                currentVal = currentVal * 10 + (c - '0');
            } else {
                if(currentVal == 0) currentVal = 1;
                
                pair<bool, bool> op;
                if(c == 'M' || c == '=' || c == 'X') {
                    op = {true, true};
                } else if(c == 'I') {
                    op = {false, true};
                } else if(c == 'D') {
                    op = {true, false};
                } else {
                    op = {true, true}; // Default
                }

                if (c == 'M' || c == '=' || c == 'X' || c == 'I' || c == 'D') {
                    // Fill alignment
                    std::fill(alignment.begin() + currentIndex, alignment.begin() + currentIndex + currentVal, op);
                    currentIndex += currentVal;

                    // Update mismatchCount and positions
                    if (op.first && op.second) { // M, =, X
                        uint64_t mismatchHere = 0;
                        for(size_t k=0; k<currentVal; ++k) {
                            if (sequence0[position0 + k] != sequence1[position1 + k]) {
                                mismatchCount++;
                                mismatchHere++;
                            }
                        }
                        const uint64_t matchHere = uint64_t(currentVal) - mismatchHere;
                        dpScore += match * int64_t(matchHere) + mismatch * int64_t(mismatchHere);
                        position0 += currentVal;
                        position1 += currentVal;
                    } else if (op.first) { // D
                        position0 += currentVal;
                        indelBaseCount += currentVal;
                        gapEventCount++;
                        dpScore -= gapPenalty(uint64_t(currentVal));
                        if(currentVal >= 6) hasLargeIndel = true;
                        if(currentVal > maxIndelSize) maxIndelSize = uint32_t(currentVal);
                    } else if (op.second) { // I
                        position1 += currentVal;
                        indelBaseCount += currentVal;
                        gapEventCount++; // ONE gap event
                        dpScore -= gapPenalty(uint64_t(currentVal));
                        if(currentVal >= 6) hasLargeIndel = true;
                        if(currentVal > maxIndelSize) maxIndelSize = uint32_t(currentVal);
                    }
                }
                currentVal = 0;
            }
        }
        
        astarpa_free_cigar((unsigned char*)cigar);

        DINARA_ASSERT(position0 == sequence0.size());
        DINARA_ASSERT(position1 == sequence1.size());
    }

}



void ProjectedAlignmentSegment::computeRleAlignment()
{
    const vector<uint8_t>& sequence0 = reinterpret_cast< const vector<uint8_t>& >(rleSequences[0]);
    const vector<uint8_t>& sequence1 = reinterpret_cast< const vector<uint8_t>& >(rleSequences[1]);

    if(sequence0 == sequence1) {
        rleEditDistance = 0;
        rleAlignment.resize(sequence0.size());
        fill(rleAlignment.begin(), rleAlignment.end(), make_pair(true, true));
        mismatchCountRle = 0;

    } else {
        // Convert sequences to ASCII for A*PA2
        // Use thread_local buffers to avoid allocation
        static thread_local vector<uint8_t> asciiSequence0;
        static thread_local vector<uint8_t> asciiSequence1;
        asciiSequence0.clear();
        asciiSequence1.clear();
        asciiSequence0.reserve(sequence0.size());
        asciiSequence1.reserve(sequence1.size());

        static const char baseToAscii[] = {'A', 'C', 'G', 'T'};

        // Fast conversion using lookup table
        for(uint8_t b : sequence0) {
            asciiSequence0.push_back(b < 4 ? baseToAscii[b] : 'N');
        }
        for(uint8_t b : sequence1) {
            asciiSequence1.push_back(b < 4 ? baseToAscii[b] : 'N');
        }

        char* cigar = nullptr;
        size_t cigarLen = 0;
        int64_t cost = astarpa2_simple(
            asciiSequence0.data(), asciiSequence0.size(),
            asciiSequence1.data(), asciiSequence1.size(),
            (unsigned char**)&cigar, &cigarLen);
        
        // Convert cost to edit distance (A*PA2 cost is usually edit distance if using defaults)
        // But here we just use the returned cost.
        rleEditDistance = cost;

        // Parse CIGAR directly from char*
        rleAlignment.clear();
        
        // Optimization: Pre-calculate total alignment length to avoid reallocations
        size_t totalAlignmentLength = 0;
        size_t currentVal = 0;
        for(size_t i=0; i<cigarLen; i++) {
            char c = cigar[i];
            if(isdigit(c)) {
                currentVal = currentVal * 10 + (c - '0');
            } else {
                if(currentVal == 0) currentVal = 1;
                // All supported operations (M, =, X, I, D, S, N) contribute to alignment vector length
                // H (hard clip) does not.
                if (c != 'H') {
                    totalAlignmentLength += currentVal;
                }
                currentVal = 0;
            }
        }
        
        rleAlignment.resize(totalAlignmentLength);

        // Fill the alignment vector AND compute mismatchCountRle
        size_t currentIndex = 0;
        currentVal = 0;
        uint64_t position0 = 0;
        uint64_t position1 = 0;
        mismatchCountRle = 0;

        for(size_t i=0; i<cigarLen; i++) {
            char c = cigar[i];
            if(isdigit(c)) {
                currentVal = currentVal * 10 + (c - '0');
            } else {
                if(currentVal == 0) currentVal = 1;
                
                pair<bool, bool> op;
                if(c == 'M' || c == '=' || c == 'X') {
                    op = {true, true};
                } else if(c == 'I') {
                    op = {false, true};
                } else if(c == 'D') {
                    op = {true, false};
                } else {
                    op = {true, true}; // Default
                }

                if (c == 'M' || c == '=' || c == 'X' || c == 'I' || c == 'D') {
                    // Fill alignment
                    std::fill(rleAlignment.begin() + currentIndex, rleAlignment.begin() + currentIndex + currentVal, op);
                    currentIndex += currentVal;

                    // Update mismatchCountRle and positions
                    if (op.first && op.second) { // M, =, X
                        for(size_t k=0; k<currentVal; ++k) {
                            if (sequence0[position0 + k] != sequence1[position1 + k]) {
                                mismatchCountRle++;
                            }
                        }
                        position0 += currentVal;
                        position1 += currentVal;
                    } else if (op.first) { // D
                        position0 += currentVal;
                    } else if (op.second) { // I
                        position1 += currentVal;
                    }
                }
                currentVal = 0;
            }
        }
        
        astarpa_free_cigar((unsigned char*)cigar);

        DINARA_ASSERT(position0 == sequence0.size());
        DINARA_ASSERT(position1 == sequence1.size());
    }

}



void ProjectedAlignment::writeHtml(ostream& html, bool brief) const
{
    html <<
        "<table>"
        "<tr>"
        "<th colspan=6>" << orientedReadIds[0] <<
        "<th colspan=6>" << orientedReadIds[1] <<
        "<th rowspan=2>Edit<br>distance"
        "<th rowspan=2>RLE<br>edit<br>distance"
        "<th rowspan=2>RLE<br>mismatch<br>count"
        "<th rowspan=2 class=left>Alignments"
        "<tr>"
        "<th>OrdinalA"
        "<th>OrdinalB"
        "<th>Skip"
        "<th>PositionA"
        "<th>PositionB"
        "<th>Length"
        "<th>OrdinalA"
        "<th>OrdinalB"
        "<th>Skip"
        "<th>PositionA"
        "<th>PositionB"
        "<th>Length";


    for(const ProjectedAlignmentSegment& segment: segments) {
        if(brief and (segment.editDistance == 0)) {
            continue;
        }
        segment.writeHtml(html);
    }

    html << "</table>";
}



void ProjectedAlignmentSegment::writeHtml(ostream& html) const
{
    html << "<tr>";

    for(uint64_t i=0; i<2; i++) {
        html << "<td class=centered>" << ordinalsA[i];
        html << "<td class=centered>" << ordinalsB[i];
        html << "<td class=centered>" << ordinalsB[i] - ordinalsA[i] - 1;
        html << "<td class=centered>" << positionsA[i];
        html << "<td class=centered>" << positionsB[i];
        html << "<td class=centered>" << positionsB[i] - positionsA[i];
    }

    html << "<td class=centered>" << editDistance;
    html << "<td class=centered>" << rleEditDistance;
    html << "<td class=centered>" << mismatchCountRle;

    html << "<td class=left style='font-family:courier'>";
    writeAlignmentHtml(html);
    html << "<br><br>";
    writeRleAlignmentHtml(html);
}



void ProjectedAlignmentSegment::writeAlignmentHtml(ostream& html) const
{
    const vector<Base>& sequence0 = sequences[0];
    const vector<Base>& sequence1 = sequences[1];

    uint64_t position0 = 0;
    uint64_t position1 = 0;
    std::ostringstream alignment0;
    std::ostringstream alignment1;

    for(const pair<bool, bool>& p: alignment) {
        const bool hasBase0 = p.first;
        const bool hasBase1 = p.second;

        if(hasBase0) {
            alignment0 << sequence0[position0++];
        } else {
            alignment0 << "-";
        }

        if(hasBase1) {
            alignment1 << sequence1[position1++];
        } else {
            alignment1 << "-";
        }

    }

    DINARA_ASSERT(position0 == sequence0.size());
    DINARA_ASSERT(position1 == sequence1.size());

    const string alignment0String = alignment0.str();
    const string alignment1String = alignment1.str();

    for(uint64_t i=0; i<alignment.size(); i++) {
        const bool isDifferent = (alignment0String[i] != alignment1String[i]);
        if(isDifferent) {
            html << "<span style='background-color:pink'>";
        }
        html << alignment0String[i];
        if(isDifferent) {
            html << "</span>";
        }
    }

    html << "<br>";

    for(uint64_t i=0; i<alignment.size(); i++) {
        const bool isDifferent = (alignment0String[i] != alignment1String[i]);
        if(isDifferent) {
            html << "<span style='background-color:pink'>";
        }
        html << alignment1String[i];
        if(isDifferent) {
            html << "</span>";
        }
    }

}



void ProjectedAlignmentSegment::writeRleAlignmentHtml(ostream& html) const
{
    const vector<Base>& sequence0 = rleSequences[0];
    const vector<Base>& sequence1 = rleSequences[1];

    uint64_t position0 = 0;
    uint64_t position1 = 0;
    std::ostringstream alignment0;
    std::ostringstream alignment1;

    for(const pair<bool, bool>& p: rleAlignment) {
        const bool hasBase0 = p.first;
        const bool hasBase1 = p.second;

        if(hasBase0) {
            alignment0 << sequence0[position0++];
        } else {
            alignment0 << "-";
        }

        if(hasBase1) {
            alignment1 << sequence1[position1++];
        } else {
            alignment1 << "-";
        }

    }

    DINARA_ASSERT(position0 == sequence0.size());
    DINARA_ASSERT(position1 == sequence1.size());

    const string alignment0String = alignment0.str();
    const string alignment1String = alignment1.str();

    for(uint64_t i=0; i<rleAlignment.size(); i++) {
        const bool isDifferent = (alignment0String[i] != alignment1String[i]);
        if(isDifferent) {
            html << "<span style='background-color:pink'>";
        }
        html << alignment0String[i];
        if(isDifferent) {
            html << "</span>";
        }
    }

    html << "<br>";

    for(uint64_t i=0; i<rleAlignment.size(); i++) {
        const bool isDifferent = (alignment0String[i] != alignment1String[i]);
        if(isDifferent) {
            html << "<span style='background-color:pink'>";
        }
        html << alignment1String[i];
        if(isDifferent) {
            html << "</span>";
        }
    }

}



void ProjectedAlignment::fillSequences(ProjectedAlignmentSegment& segment) const
{
    for(uint64_t i=0; i<2; i++) {
        vector<Base>& sequence = segment.sequences[i];
        sequence.clear();
        const uint32_t begin = segment.positionsA[i];
        const uint32_t end = segment.positionsB[i];

        // It is possible (due to degenerate marker placement or unexpected marker alignments)
        // that two consecutive aligned marker midpoints map to the same base position, or even
        // that begin/end are not in increasing order. In these cases, this segment contains
        // no bases on this side; allow an empty sequence instead of asserting/crashing.
        if (begin >= end) {
            continue;
        }

        for(uint32_t position=begin; position!=end; position++) {
            sequence.push_back(getBase(i, position));
        }
    }
}




Base ProjectedAlignment::getBase(uint64_t i, uint32_t position) const
{
    const LongBaseSequenceView& sequence = sequences[i];

    if(orientedReadIds[i].getStrand() == 0) {
        return sequence[position];
    } else {
        return sequence[sequence.baseCount - 1 - position].complement();
    }
}



void ProjectedAlignmentSegment::fillRleSequences()
{
    for(uint64_t i=0; i<2; i++) {
        const vector<Base>& sequence = sequences[i];
        vector<Base>& rleSequence = rleSequences[i];
        rleSequence.clear();

        for(const Base b: sequence) {
            if(rleSequence.empty()) {
                rleSequence.push_back(b);
            } else {
                if(rleSequence.back() != b) {
                    rleSequence.push_back(b);
                }
            }
        }
    }
}



void ProjectedAlignment::computeStatistics()
{
    // Compute total lengths.
    totalLength = {0, 0};
    totalLengthRle = {0, 0};
    for(const ProjectedAlignmentSegment& segment: segments) {
        for(uint64_t i=0; i<2; i++) {
            totalLength[i] += segment.sequences[i].size();
            totalLengthRle[i] += segment.rleSequences[i].size();
        }
    }

    // Sanity check on the total lengths.
    for(uint64_t i=0; i<2; i++) {
        DINARA_ASSERT(totalLength[i] == segments.back().positionsB[i] - segments.front().positionsA[i]);
    }

    // Compute total edit distances.
    totalEditDistance = 0;
    totalEditDistanceRle = 0;
    mismatchCount = 0;
    totalIndelBaseCount = 0;
    totalGapEventCount = 0;
    totalDpScore = 0;
    hasLargeIndel = false;
    for(const ProjectedAlignmentSegment& segment: segments) {
        totalEditDistance += segment.editDistance;
        totalEditDistanceRle += segment.rleEditDistance;
        totalIndelBaseCount += segment.indelBaseCount;
        totalGapEventCount += segment.gapEventCount;
        mismatchCount += segment.mismatchCount;
        totalDpScore += segment.dpScore;
        if (segment.hasLargeIndel) hasLargeIndel = true;
    }
}



double ProjectedAlignment::errorRate() const
{
    // One-sided denominator (query/read0 length only), matching hifiasm's
    // non_trim_error_rate which uses tErr/tLen where tLen = query x-length.
    return double(totalEditDistance) / double(totalLength[0]);
}



double ProjectedAlignment::errorRateGaps() const
{
    // One-sided denominator, consistent with errorRate().
    return double(totalIndelBaseCount) / double(totalLength[0]);
}



double ProjectedAlignment::errorRateRle() const
{
    return double(totalEditDistanceRle) / double(totalLengthRle[0] + totalLengthRle[1]);
}



double ProjectedAlignment::Q() const
{
    const double er = errorRate();
    DINARA_ASSERT(er > 0.);
    return -10. * log10(er);
}



double ProjectedAlignment::QRle() const
{
    const double er = errorRateRle();
    DINARA_ASSERT(er > 0.);
    return -10. * log10(er);
}



void ProjectedAlignment::writeStatisticsHtml(ostream& html) const
{
    using std::fixed;
    using std::setprecision;

    html << "<table>";

    // Header line.
    html <<
        "<tr>"
        "<th>"
        "<th>" << orientedReadIds[0] <<
        "<th>" << orientedReadIds[1] <<
        "<th>Total";

    // Length.
    html <<
        "<tr>"
        "<th class=left>Length of aligned portion" <<
        "<td class=centered>" << totalLength[0] <<
        "<td class=centered>" << totalLength[1] <<
        "<td class=centered>" << totalLength[0] + totalLength[1];

    // RLE length.
    html <<
        "<tr>"
        "<th class=left>RLE Length of aligned portion" <<
        "<td class=centered>" << totalLengthRle[0] <<
        "<td class=centered>" << totalLengthRle[1] <<
        "<td class=centered>" << totalLengthRle[0] + totalLengthRle[1];

    // Edit distance.
    html <<
        "<tr>"
        "<th class=left>Edit distance" <<
        "<td colspan=3 class=centered>" << totalEditDistance;

    // RLE edit distance.
    html <<
        "<tr>"
        "<th class=left>RLE edit distance" <<
        "<td colspan=3 class=centered>" << totalEditDistanceRle;

    // Error rate.
    html <<
        "<tr>"
        "<th class=left>Error rate" <<
        "<td colspan=3 class=centered>" << errorRate();

    // RLE error rate.
    html <<
        "<tr>"
        "<th class=left>RLE error rate" <<
        "<td colspan=3 class=centered>" << errorRateRle();

    // Q.
    html <<
        "<tr>"
        "<th class=left>Q (dB)" <<
        "<td colspan=3 class=centered>" << fixed << setprecision(1) << Q();

    // RLE Q.
    html <<
        "<tr>"
        "<th class=left>RLE Q (dB)" <<
        "<td colspan=3 class=centered>" << fixed << setprecision(1) << QRle();

    // RLE mismatch count
    html <<
        "<tr>"
        "<th class=left>RLE mismatch count" <<
        "<td colspan=3 class=centered>" << fixed << setprecision(1) << mismatchCountRle;

    html << "</table>";
}



// Find pairs of mismatching positions in the raw alignments.
void ProjectedAlignment::getMismatchPositions(vector< array<uint32_t, 2> >& mismatchPositions) const
{

    // Start with no mismatches.
    mismatchPositions.clear();

    // Loop over all segments.
    for(const ProjectedAlignmentSegment& segment: segments) {

        // Get the sequences and the alignment for this segment.
        const auto& sequences = segment.sequences;
        const vector< pair<bool, bool> >& alignment = segment.alignment;

        // Loop over the alignment.
        uint32_t positionOffset0 = 0;
        uint32_t positionOffset1 = 0;
        for(const pair<bool, bool>& p: alignment) {
            const bool isBase0 = p.first;
            const bool isBase1 = p.second;

            // If neither is a gap, check if they are the same.
            if(isBase0 and isBase1) {
                const Base base0 = sequences[0][positionOffset0];
                const Base base1 = sequences[1][positionOffset1];

                // If not the same, store these mismatch positions.
                if(base0 != base1) {
                    mismatchPositions.push_back({
                        segment.positionsA[0] + positionOffset0,
                        segment.positionsA[1] + positionOffset1});
                }
            }

            // Increment the position offsets.
            if(isBase0) {
                ++positionOffset0;
            }
            if(isBase1) {
                ++positionOffset1;
            }
        }
        DINARA_ASSERT(positionOffset0 == segment.sequences[0].size());
        DINARA_ASSERT(positionOffset1 == segment.sequences[1].size());
    }
}
