#pragma once

/******************************************************************

Class ProjectedSAlignment describes the "projection" of
an Alignment in marker space to base space.

******************************************************************/

#include "array.hpp"
#include "cstdint.hpp"
#include "Base.hpp"
#include "invalid.hpp"
#include "iosfwd.hpp"
#include "span.hpp"
#include "utility.hpp"
#include "vector.hpp"
#include <astarpa/astarpa.h>

namespace dinara {
    class ProjectedAlignment;
    class ProjectedAlignmentSegment;

    class Alignment;
    class Assembler;
    class CompressedMarker;
    class LongBaseSequenceView;
    class OrientedReadId;

    struct ProjectedAlignmentSparseMismatch {
        uint32_t position0;
        uint32_t position1;
        uint8_t base0;
        uint8_t base1;
    };

    struct ProjectedAlignmentSparseIndel {
        uint32_t position0;
        uint32_t position1;
        uint32_t length;
        char op; // 'I' or 'D'
    };
}



// Each ProjectedAlignmentSegment describes a pair of consecutive
// aligned markers in the Alignment in marker space.
class dinara::ProjectedAlignmentSegment {
public:

    ProjectedAlignmentSegment(
        uint32_t kHalf,
        const array<uint32_t, 2>& ordinalA,
        const array<uint32_t, 2>& ordinalB,
        const array< span<const CompressedMarker>, 2>& markers);

    ProjectedAlignmentSegment() {}

    // In the arrays below, the index can be 0 or 1 and correspond
    // to the first and second oriented reads in the Alignment.
    // The A and B suffix refer to the left and right markers
    // in the pair of consecutive aligned markers.

    // The ordinals of the pair of consecutive aligned markers
    // described by this ProjectedAlignmentSegment.
    array<uint32_t, 2> ordinalsA;
    array<uint32_t, 2> ordinalsB;

    // The begin/end positions in base space, taken
    // at the midpoints of the two consecutive aligned markers.
    array<uint32_t, 2> positionsA;
    array<uint32_t, 2> positionsB;

    // The corresponding Base sequences.
    array<vector<Base>, 2> sequences;

    // The alignment between the two sequences.
    // See seqan.hpp for its meaning.
    int64_t editDistance = invalid<int64_t>;
    vector< pair<bool, bool> > alignment;
    void computeAlignment(
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore,
        int64_t dpMatchScore,
        int64_t dpMismatchScore,
        int64_t dpGapOpen1,
        int64_t dpGapExtend1,
        int64_t dpGapOpen2,
        int64_t dpGapExtend2);
    void computeAlignmentSparse(
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore,
        int64_t dpMatchScore,
        int64_t dpMismatchScore,
        int64_t dpGapOpen1,
        int64_t dpGapExtend1,
        int64_t dpGapOpen2,
        int64_t dpGapExtend2,
        vector<ProjectedAlignmentSparseMismatch>& sparseMismatches,
        vector<ProjectedAlignmentSparseIndel>& sparseIndels);

    // The Base sequences in RLE represenation.
    array<vector<Base>, 2> rleSequences;
    void fillRleSequences();

    // The alignment between the two RLE sequences.
    // See seqan.hpp for its meaning.
    int64_t rleEditDistance = invalid<int64_t>;
    vector< pair<bool, bool> > rleAlignment;
    void computeRleAlignment(
        int64_t matchScore,
        int64_t mismatchScore,
        int64_t gapScore);

    // The number of mismatches in the RLE alignment.
    uint64_t mismatchCountRle = invalid<uint64_t>;

    // The number of mismatches in the raw alignment.
    uint64_t mismatchCount = invalid<uint64_t>;

    // Total errors excluding those attributed to homopolymer-repeat context
    // using hifiasm's if_is_homopolymer_repeat-style accounting.
    uint64_t nonHomopolymerErrorCount = invalid<uint64_t>;
    
    // The number of deletions in the raw alignment.
    uint64_t deletionCount = invalid<uint64_t>;

    // The number of gap events (indels).
    uint64_t gapEventCount = invalid<uint64_t>;

    // Base-level DP alignment score computed from the CIGAR using hifiasm's current
    // single-affine overlap scoring model.
    int64_t dpScore = invalid<int64_t>;

    // Flag for large indel (>= 6 bases)
    bool hasLargeIndel = false;

    // Largest indel size
    uint32_t maxIndelSize = 0;


    void writeAlignmentHtml(ostream&) const;
    void writeRleAlignmentHtml(ostream&) const;

    void writeHtml(ostream&) const;
};



// Class ProjectedSAlignment describes the "projection" of
// an Alignment in marker space to base space.
// It is a sequence of ProjectedAlignmentSegments as defined above.
class dinara::ProjectedAlignment {
public:
    vector<ProjectedAlignmentSegment> segments;
    vector<ProjectedAlignmentSparseMismatch> sparseMismatches;
    vector<ProjectedAlignmentSparseIndel> sparseIndels;

    enum class Method {
        All,        // Do both RLE and raw alignments, store all segments.
        QuickRle,   // Only do RLE alignments, only store segments where the two RLE sequences differ.
        QuickRaw,   // Only do raw alignments, only store segments where the two raw sequences differ.
        QuickRawSparse, // Only do raw alignments, store sparse diffs (no full alignment trace).
    };

    ProjectedAlignment(
        const Assembler&,
        const array<OrientedReadId, 2>&,
        const Alignment&,
        Method method,
        int64_t dpMatchScore,
        int64_t dpMismatchScore,
        int64_t dpGapOpen1,
        int64_t dpGapExtend1,
        int64_t dpGapOpen2,
        int64_t dpGapExtend2);

    ProjectedAlignment(
        uint32_t k,
        const array<OrientedReadId, 2>&,
        const array<LongBaseSequenceView, 2>&,
        const Alignment&,
        const array< span<const CompressedMarker>, 2>& markers,
        Method method,
        int64_t dpMatchScore,
        int64_t dpMismatchScore,
        int64_t dpGapOpen1,
        int64_t dpGapExtend1,
        int64_t dpGapOpen2,
        int64_t dpGapExtend2);

    void constructAll();
    void constructQuickRle();
    void constructQuickRaw();
    void constructQuickRawSparse();
    
    // Flag to indicate if the projected alignment touches the ends of the markers.
    bool touchesMarkerEnds = false;



    // Marker length and its half.
    uint32_t k;
    uint32_t kHalf;

    // Scoring scheme for edit distance.
    const int64_t matchScore = 0;
    const int64_t mismatchScore = -1;
    const int64_t gapScore = -1;

    // Scoring scheme for overlap DP score.
    // Hifiasm's current overlap path uses single-affine scoring (gapo/gape).
    // We keep the second gap parameter pair in the API for compatibility with
    // existing option plumbing, but overlap DP scoring ignores it.
    const int64_t dpMatchScore;
    const int64_t dpMismatchScore;
    const int64_t dpGapOpen1;
    const int64_t dpGapExtend1;
    const int64_t dpGapOpen2;
    const int64_t dpGapExtend2;

    // The two OrientedReadIds in this alignment.
    const array<OrientedReadId, 2>& orientedReadIds;

    // The base sequences of the reads.
    // These require reverse complementing for OrientedReadIds on strand 1.
    const array<LongBaseSequenceView, 2>& sequences;

    void fillSequences(ProjectedAlignmentSegment&) const;

    // Get the Base at a given position in one of the two oriented reads,
    // doing a reverse complement if necessary.
    Base getBase(uint64_t i, uint32_t position) const;

    // The input Alignment in marker space.
    const Alignment& alignment;

    // The markers for the two oriented reads in this alignment.
    const array< span<const CompressedMarker>, 2>& markers;



    // Statistics for the entire ProjectedAlignment.

    // Number of aligned bases (raw and RLE) in the aligned portions of the two oriented reads.
    array<uint64_t, 2> totalLength;
    array<uint64_t, 2> totalLengthRle;

    // Total edit distance (raw and RLE).
    int64_t totalEditDistance;
    int64_t totalEditDistanceRle;

    // The number of mismatches in the RLE alignment.
    uint64_t mismatchCountRle;

    // The number of mismatches in the raw alignment.
    uint64_t mismatchCount;

    // Total errors excluding homopolymer-repeat-associated errors.
    uint64_t nonHomopolymerErrorCount;
    
    // The number of deletions in the raw alignment.
    uint64_t totalDeletionCount; // Total bases in deletions

    // The number of gap events (indels) in the raw alignment.
    uint64_t totalGapEventCount;

    // Base-level DP score for the entire projected alignment (sum of segment scores).
    int64_t totalDpScore = invalid<int64_t>;

    // Flag for large indel (>= 6 bases)
    bool hasLargeIndel = false;

    // Largest indel size
    uint32_t maxIndelSize = 0;

    void computeStatistics();
    double errorRate() const;
    double errorRateGaps() const;
    double errorRateRle() const;
    double Q() const;
    double QRle() const;

    void writeStatisticsHtml(ostream&) const;
    void writeHtml(ostream&, bool brief) const;

    // Find pairs of mismatching positions in the raw alignments.
    void getMismatchPositions(vector< array<uint32_t, 2> >&) const;
};
