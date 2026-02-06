#include "ProjectedAlignment.hpp"
#include "Alignment.hpp"
#include "Assembler.hpp"
#include "Base.hpp"
#include "LongBaseSequence.hpp"
#include "Marker.hpp"
#include "Reads.hpp"
#include "seqan.hpp"
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
    int64_t dpGapExtend1,
    int64_t dpGapOpen2,
    int64_t dpGapExtend2) :

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
        dpGapExtend1,
        dpGapOpen2,
        dpGapExtend2)
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
    int64_t dpGapOpen2,
    int64_t dpGapExtend2) :
    k(k),
    kHalf(k / 2),
    orientedReadIds(orientedReadIds),
    sequences(sequences),
    alignment(alignment),
    markers(markers),
    dpMatchScore(dpMatchScore),
    dpMismatchScore(dpMismatchScore),
    dpGapOpen1(dpGapOpen1),
    dpGapExtend1(dpGapExtend1),
    dpGapOpen2(dpGapOpen2),
    dpGapExtend2(dpGapExtend2)
{
    DINARA_ASSERT((k % 2) == 0);

    switch(method) {
        case Method::All:
            constructAll();
            break;
        case Method::QuickRle:
            constructQuickRle();
            break;
        case Method::QuickRaw:
            constructQuickRaw();
            break;
        case Method::QuickRawSparse:
            constructQuickRawSparse();
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
            matchScore, mismatchScore, gapScore,
            dpMatchScore, dpMismatchScore,
            dpGapOpen1, dpGapExtend1, dpGapOpen2, dpGapExtend2);

        // Same, in RLE.
        segment.fillRleSequences();
        segment.computeRleAlignment(matchScore, mismatchScore, gapScore);
        mismatchCountRle += segment.mismatchCountRle;
    }



    computeStatistics();
}



// This stores only the following:
// - RLE sequences and RLE alignments for segments for which the RLE sequences
//   of the two oriented reads are different.
// - Total RLE edit distance and total RLE lengths.
void ProjectedAlignment::constructQuickRle()
{
    // Create the segment outside the loop and reuse it to reduce
    // memory allocation activity.
    ProjectedAlignmentSegment segment;

    totalLengthRle = {0, 0};
    totalEditDistanceRle = 0;
    mismatchCountRle = 0;

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

        // Store RLE sequences, without going through the raw sequences for speed.
        // Also increment the total RLE lengths.
        for(uint64_t i=0; i<2; i++) {
            vector<Base>& rleSequence = segment.rleSequences[i];
            rleSequence.clear();
            for(uint32_t position=segment.positionsA[i]; position!=segment.positionsB[i]; position++) {
                const Base b = getBase(i, position);
                if(rleSequence.empty() or b != rleSequence.back()) {
                    rleSequence.push_back(b);
                }
            }
            totalLengthRle[i] += rleSequence.size();
        }

        // If the RLE sequences are the same, there is no contribution to RLE edit distance,
        // and we don't store the segment.
        if(segment.rleSequences[0] == segment.rleSequences[1]) {
            continue;
        }

        // Otherwise, we compute the RLE alignment and store this segment.
        segment.computeRleAlignment(matchScore, mismatchScore, gapScore);
        totalEditDistanceRle += segment.rleEditDistance;
        mismatchCountRle += segment.mismatchCountRle;
        segments.push_back(segment);
    }
}



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
    totalDeletionCount = 0;
    totalGapEventCount = 0;
    totalDpScore = 0;
    hasLargeIndel = false;
    maxIndelSize = 0;

    // Hifiasm/minimap2-style two-piece affine parameters (HiFi defaults).
    // For identical segments we only need the match reward.
    const int64_t dpMatch = dpMatchScore;

    // Loop over pairs of consecutive aligned markers (A, B).
    // Prepend the initial kHalf match (Left Tail)
    // Note: CIGAR generation moved to AlignedEvidenceStore (AssemblerAlign.cpp). Removed here.


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
            matchScore, mismatchScore, gapScore,
            dpMatchScore, dpMismatchScore,
            dpGapOpen1, dpGapExtend1, dpGapOpen2, dpGapExtend2);

        // Accumulate statistics
        totalEditDistance += segment.editDistance;
        mismatchCount += segment.mismatchCount;
        totalDeletionCount += segment.deletionCount;
        totalGapEventCount += segment.gapEventCount;
        totalDpScore += segment.dpScore;
        if (segment.hasLargeIndel) hasLargeIndel = true;
        if (segment.maxIndelSize > maxIndelSize) maxIndelSize = segment.maxIndelSize;

        // Store the segment.
        segments.push_back(segment);
    }

    // Append the final kHalf match (Right Tail)
    // Removed.
}



// Similar to constructQuickRaw, but does not store per-base alignment traces.
// Instead, it stores sparse mismatch and indel events that are sufficient to
// populate AlignedEvidenceStore efficiently.
void ProjectedAlignment::constructQuickRawSparse()
{
    // Reuse a single segment object to reduce allocations.
    ProjectedAlignmentSegment segment;

    sparseMismatches.clear();
    sparseIndels.clear();

    // Initialize statistics.
    totalLength = {0, 0};
    totalEditDistance = 0;
    mismatchCount = 0;
    totalDeletionCount = 0;
    totalGapEventCount = 0;
    totalDpScore = 0;
    hasLargeIndel = false;
    maxIndelSize = 0;

    // Hifiasm/minimap2-style two-piece affine parameters (HiFi defaults).
    // For identical segments we only need the match reward.
    const int64_t dpMatch = dpMatchScore;

    // Thread-local buffers to avoid repeated allocations/conversions.
    // We fill these directly from the oriented read views (via getBase),
    // avoiding the intermediate vector<Base> in ProjectedAlignmentSegment.
    static thread_local vector<uint8_t> asciiSequence0;
    static thread_local vector<uint8_t> asciiSequence1;
    static constexpr uint8_t baseToAscii[4] = {'A', 'C', 'G', 'T'};
    static constexpr array<uint8_t, 256> asciiToBase = []() constexpr {
        array<uint8_t, 256> m{};
        for(size_t i = 0; i < 256; ++i) m[i] = 0;
        m[uint8_t('A')] = 0;
        m[uint8_t('C')] = 1;
        m[uint8_t('G')] = 2;
        m[uint8_t('T')] = 3;
        return m;
    }();

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

        const uint32_t begin0 = segment.positionsA[0];
        const uint32_t end0 = segment.positionsB[0];
        const uint32_t begin1 = segment.positionsA[1];
        const uint32_t end1 = segment.positionsB[1];
        const uint32_t len0 = end0 - begin0;
        const uint32_t len1 = end1 - begin1;

        // Accumulate total lengths (even for identical sequences).
        totalLength[0] += len0;
        totalLength[1] += len1;

        // Empty segments should not happen, but tolerate.
        if(len0 == 0 || len1 == 0) {
            continue;
        }

        // Fill ASCII sequences directly from oriented reads.
        asciiSequence0.resize(len0);
        asciiSequence1.resize(len1);
        for(uint32_t j = 0; j < len0; ++j) {
            const uint8_t b = getBase(0, begin0 + j).value;
            asciiSequence0[j] = baseToAscii[b];
        }
        for(uint32_t j = 0; j < len1; ++j) {
            const uint8_t b = getBase(1, begin1 + j).value;
            asciiSequence1[j] = baseToAscii[b];
        }

        // If the raw sequences are the same, there is no contribution.
        if(asciiSequence0 == asciiSequence1) {
            totalDpScore += dpMatch * int64_t(len0);
            continue;
        }

        // Align them and collect sparse diffs.
        char* cigar = nullptr;
        size_t cigarLen = 0;
        const int64_t cost = astarpa2_simple(
            asciiSequence0.data(), asciiSequence0.size(),
            asciiSequence1.data(), asciiSequence1.size(),
            (unsigned char**)&cigar, &cigarLen);

        const int64_t segEditDistance = cost;
        uint64_t segMismatchCount = 0;
        uint64_t segDeletionCount = 0;
        uint64_t segGapEventCount = 0;
        int64_t segDpScore = 0;
        bool segHasLargeIndel = false;
        uint32_t segMaxIndelSize = 0;

        // Hifiasm/minimap2-style two-piece affine parameters (HiFi defaults):
        // match=+2, mismatch=-4, gapCost(L)=min(O1+E1*L, O2+E2*L) with (O1,E1)=(4,2), (O2,E2)=(24,1).
        const int64_t match = dpMatchScore;
        const int64_t mismatch = dpMismatchScore;
        const int64_t gapOpen1 = dpGapOpen1;
        const int64_t gapExtend1 = dpGapExtend1;
        const int64_t gapOpen2 = dpGapOpen2;
        const int64_t gapExtend2 = dpGapExtend2;
        auto gapPenalty = [&](uint64_t length) -> int64_t {
            const int64_t l = int64_t(length);
            DINARA_ASSERT(l >= 1);
            // Minimap2/ksw2 convention: a gap of length k costs O + k*E.
            const int64_t c1 = gapOpen1 + gapExtend1 * l;
            const int64_t c2 = gapOpen2 + gapExtend2 * l;
            return std::min(c1, c2);
        };

        // Parse CIGAR and collect sparse differences.
        uint64_t position0 = 0;
        uint64_t position1 = 0;
        size_t currentVal = 0;
        for(size_t i=0; i<cigarLen; i++) {
            const char c = cigar[i];
            if(isdigit(c)) {
                currentVal = currentVal * 10 + size_t(c - '0');
                continue;
            }

            if(currentVal == 0) {
                currentVal = 1;
            }

            if(c == 'M' || c == '=' || c == 'X') {
                // Match/mismatch block: record mismatches.
                uint64_t mismatchHere = 0;
                for(size_t k=0; k<currentVal; ++k) {
                    const uint8_t a0 = asciiSequence0[position0 + k];
                    const uint8_t a1 = asciiSequence1[position1 + k];
                    if(a0 != a1) {
                        sparseMismatches.push_back(ProjectedAlignmentSparseMismatch{
                            uint32_t(begin0 + uint32_t(position0 + k)),
                            uint32_t(begin1 + uint32_t(position1 + k)),
                            asciiToBase[a0],
                            asciiToBase[a1]
                        });
                        ++segMismatchCount;
                        ++mismatchHere;
                    }
                }
                const uint64_t matchHere = uint64_t(currentVal) - mismatchHere;
                segDpScore += match * int64_t(matchHere) + mismatch * int64_t(mismatchHere);
                position0 += currentVal;
                position1 += currentVal;

            } else if(c == 'D') {
                // Gap in sequence1 (target).
                sparseIndels.push_back(ProjectedAlignmentSparseIndel{
                    uint32_t(begin0 + uint32_t(position0)),
                    uint32_t(begin1 + uint32_t(position1)),
                    uint32_t(currentVal),
                    'D'
                });
                position0 += currentVal;
                segDeletionCount += currentVal;
                ++segGapEventCount;
                segDpScore -= gapPenalty(uint64_t(currentVal));
                if(currentVal >= 6) {
                    segHasLargeIndel = true;
                }
                if(currentVal > segMaxIndelSize) {
                    segMaxIndelSize = uint32_t(currentVal);
                }

            } else if(c == 'I') {
                // Gap in sequence0 (query).
                sparseIndels.push_back(ProjectedAlignmentSparseIndel{
                    uint32_t(begin0 + uint32_t(position0)),
                    uint32_t(begin1 + uint32_t(position1)),
                    uint32_t(currentVal),
                    'I'
                });
                position1 += currentVal;
                segDeletionCount += currentVal;
                ++segGapEventCount;
                segDpScore -= gapPenalty(uint64_t(currentVal));
                if(currentVal >= 6) {
                    segHasLargeIndel = true;
                }
                if(currentVal > segMaxIndelSize) {
                    segMaxIndelSize = uint32_t(currentVal);
                }

            } else if(c == 'H') {
                // Hard clips do not consume sequence.
                // astarpa2_simple is not expected to generate these, but tolerate.
            } else {
                // Unexpected CIGAR op for A*PA2 simple alignment.
                DINARA_ASSERT(0);
            }

            currentVal = 0;
        }

        astarpa_free_cigar((unsigned char*)cigar);

        DINARA_ASSERT(position0 == asciiSequence0.size());
        DINARA_ASSERT(position1 == asciiSequence1.size());

        // Accumulate statistics.
        totalEditDistance += segEditDistance;
        mismatchCount += segMismatchCount;
        totalDeletionCount += segDeletionCount;
        totalGapEventCount += segGapEventCount;
        totalDpScore += segDpScore;
        if(segHasLargeIndel) {
            hasLargeIndel = true;
        }
        if(segMaxIndelSize > maxIndelSize) {
            maxIndelSize = segMaxIndelSize;
        }
    }
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
    int64_t /* matchScore */,
    int64_t /* mismatchScore */,
    int64_t /* gapScore */,
    int64_t dpMatchScore,
    int64_t dpMismatchScore,
    int64_t dpGapOpen1,
    int64_t dpGapExtend1,
    int64_t dpGapOpen2,
    int64_t dpGapExtend2)
{
    const vector<uint8_t>& sequence0 = reinterpret_cast< const vector<uint8_t>& >(sequences[0]);
    const vector<uint8_t>& sequence1 = reinterpret_cast< const vector<uint8_t>& >(sequences[1]);

    if(sequence0 == sequence1) {
        editDistance = 0;
        alignment.resize(sequence0.size());
        fill(alignment.begin(), alignment.end(), make_pair(true, true));
        mismatchCount = 0;
        deletionCount = 0;
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
        deletionCount = 0;
        gapEventCount = 0;
        dpScore = 0;
        hasLargeIndel = false;
        maxIndelSize = 0;

        const int64_t match = dpMatchScore;
        const int64_t mismatch = dpMismatchScore;
        const int64_t gapOpen1 = dpGapOpen1;
        const int64_t gapExtend1 = dpGapExtend1;
        const int64_t gapOpen2 = dpGapOpen2;
        const int64_t gapExtend2 = dpGapExtend2;
        auto gapPenalty = [&](uint64_t length) -> int64_t {
            const int64_t l = int64_t(length);
            DINARA_ASSERT(l >= 1);
            // Minimap2/ksw2 convention: a gap of length k costs O + k*E.
            const int64_t c1 = gapOpen1 + gapExtend1 * l;
            const int64_t c2 = gapOpen2 + gapExtend2 * l;
            return std::min(c1, c2);
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
                        deletionCount += currentVal;
                        gapEventCount++;
                        dpScore -= gapPenalty(uint64_t(currentVal));
                        if(currentVal >= 6) hasLargeIndel = true;
                        if(currentVal > maxIndelSize) maxIndelSize = uint32_t(currentVal);
                    } else if (op.second) { // I
                        position1 += currentVal;
                        deletionCount += currentVal;
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



void ProjectedAlignmentSegment::computeAlignmentSparse(
    int64_t /* matchScore */,
    int64_t /* mismatchScore */,
    int64_t /* gapScore */,
    int64_t dpMatchScore,
    int64_t dpMismatchScore,
    int64_t dpGapOpen1,
    int64_t dpGapExtend1,
    int64_t dpGapOpen2,
    int64_t dpGapExtend2,
    vector<ProjectedAlignmentSparseMismatch>& sparseMismatches,
    vector<ProjectedAlignmentSparseIndel>& sparseIndels)
{
    const vector<uint8_t>& sequence0 = reinterpret_cast< const vector<uint8_t>& >(sequences[0]);
    const vector<uint8_t>& sequence1 = reinterpret_cast< const vector<uint8_t>& >(sequences[1]);

    // This is only called when the sequences differ.
    DINARA_ASSERT(sequence0 != sequence1);

    // Convert sequences to ASCII for A*PA2 using thread_local buffers to avoid allocations.
    static thread_local vector<uint8_t> asciiSequence0;
    static thread_local vector<uint8_t> asciiSequence1;
    asciiSequence0.clear();
    asciiSequence1.clear();
    asciiSequence0.reserve(sequence0.size());
    asciiSequence1.reserve(sequence1.size());

    static const char baseToAscii[] = {'A', 'C', 'G', 'T'};

    for(uint8_t b : sequence0) {
        asciiSequence0.push_back(b < 4 ? baseToAscii[b] : 'N');
    }
    for(uint8_t b : sequence1) {
        asciiSequence1.push_back(b < 4 ? baseToAscii[b] : 'N');
    }

    char* cigar = nullptr;
    size_t cigarLen = 0;
    const int64_t cost = astarpa2_simple(
        asciiSequence0.data(), asciiSequence0.size(),
        asciiSequence1.data(), asciiSequence1.size(),
        (unsigned char**)&cigar, &cigarLen);

    editDistance = cost;

    // Parse CIGAR and collect sparse differences.
    uint64_t position0 = 0;
    uint64_t position1 = 0;
    mismatchCount = 0;
    deletionCount = 0;
    gapEventCount = 0;
    dpScore = 0;
    hasLargeIndel = false;
    maxIndelSize = 0;

    const int64_t match = dpMatchScore;
    const int64_t mismatch = dpMismatchScore;
    const int64_t gapOpen1 = dpGapOpen1;
    const int64_t gapExtend1 = dpGapExtend1;
    const int64_t gapOpen2 = dpGapOpen2;
    const int64_t gapExtend2 = dpGapExtend2;
    auto gapPenalty = [&](uint64_t length) -> int64_t {
        const int64_t l = int64_t(length);
        DINARA_ASSERT(l >= 1);
        // Minimap2/ksw2 convention: a gap of length k costs O + k*E.
        const int64_t c1 = gapOpen1 + gapExtend1 * l;
        const int64_t c2 = gapOpen2 + gapExtend2 * l;
        return std::min(c1, c2);
    };

    size_t currentVal = 0;
    for(size_t i=0; i<cigarLen; i++) {
        const char c = cigar[i];
        if(isdigit(c)) {
            currentVal = currentVal * 10 + size_t(c - '0');
            continue;
        }

        if(currentVal == 0) {
            currentVal = 1;
        }

        if(c == 'M' || c == '=' || c == 'X') {
            // Match/mismatch block: record mismatches.
            uint64_t mismatchHere = 0;
            for(size_t k=0; k<currentVal; ++k) {
                const uint8_t b0 = sequence0[position0 + k];
                const uint8_t b1 = sequence1[position1 + k];
                if(b0 != b1) {
                    sparseMismatches.push_back(ProjectedAlignmentSparseMismatch{
                        uint32_t(positionsA[0] + uint32_t(position0 + k)),
                        uint32_t(positionsA[1] + uint32_t(position1 + k)),
                        b0,
                        b1
                    });
                    ++mismatchCount;
                    ++mismatchHere;
                }
            }
            const uint64_t matchHere = uint64_t(currentVal) - mismatchHere;
            dpScore += match * int64_t(matchHere) + mismatch * int64_t(mismatchHere);
            position0 += currentVal;
            position1 += currentVal;

        } else if(c == 'D') {
            // Gap in sequence1 (target).
            sparseIndels.push_back(ProjectedAlignmentSparseIndel{
                uint32_t(positionsA[0] + uint32_t(position0)),
                uint32_t(positionsA[1] + uint32_t(position1)),
                uint32_t(currentVal),
                'D'
            });
            position0 += currentVal;
            deletionCount += currentVal;
            ++gapEventCount;
            dpScore -= gapPenalty(uint64_t(currentVal));
            if(currentVal >= 6) {
                hasLargeIndel = true;
            }
            if(currentVal > maxIndelSize) {
                maxIndelSize = uint32_t(currentVal);
            }

        } else if(c == 'I') {
            // Gap in sequence0 (query).
            sparseIndels.push_back(ProjectedAlignmentSparseIndel{
                uint32_t(positionsA[0] + uint32_t(position0)),
                uint32_t(positionsA[1] + uint32_t(position1)),
                uint32_t(currentVal),
                'I'
            });
            position1 += currentVal;
            deletionCount += currentVal;
            ++gapEventCount;
            dpScore -= gapPenalty(uint64_t(currentVal));
            if(currentVal >= 6) {
                hasLargeIndel = true;
            }
            if(currentVal > maxIndelSize) {
                maxIndelSize = uint32_t(currentVal);
            }

        } else if(c == 'H') {
            // Hard clips do not consume sequence.
            // astarpa2_simple is not expected to generate these, but tolerate.
        } else {
            // Unexpected CIGAR op for A*PA2 simple alignment.
            DINARA_ASSERT(0);
        }

        currentVal = 0;
    }

    astarpa_free_cigar((unsigned char*)cigar);

    DINARA_ASSERT(position0 == sequence0.size());
    DINARA_ASSERT(position1 == sequence1.size());
}



void ProjectedAlignmentSegment::computeRleAlignment(
    int64_t /* matchScore */,
    int64_t /* mismatchScore */,
    int64_t /* gapScore */)
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
    totalDeletionCount = 0;
    totalGapEventCount = 0; // Init
    totalDpScore = 0;
    hasLargeIndel = false;
    maxIndelSize = 0;
    for(const ProjectedAlignmentSegment& segment: segments) {
        totalEditDistance += segment.editDistance;
        totalEditDistanceRle += segment.rleEditDistance;
        totalDeletionCount += segment.deletionCount;
        totalGapEventCount += segment.gapEventCount; // Sum
        mismatchCount += segment.mismatchCount;
        totalDpScore += segment.dpScore;
        if (segment.hasLargeIndel) hasLargeIndel = true;
        if (segment.maxIndelSize > maxIndelSize) maxIndelSize = segment.maxIndelSize;
    }
}



double ProjectedAlignment::errorRate() const
{
    return double(totalEditDistance) / double(totalLength[0] + totalLength[1]);
}



double ProjectedAlignment::errorRateGaps() const
{
    return double(totalDeletionCount) / double(totalLength[0] + totalLength[1]);
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
