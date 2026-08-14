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

#include "OverlapCigarStore.hpp"

namespace dinara {
    class ProjectedAlignment;
    class ProjectedAlignmentSegment;

    class Alignment;
    class Assembler;
    class CompressedMarker;
    class LongBaseSequenceView;
    class OrientedReadId;
    class OverlapCigarStore;

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
        int64_t dpMatchScore,
        int64_t dpMismatchScore,
        int64_t dpGapOpen1,
        int64_t dpGapExtend1);
    // The Base sequences in RLE represenation.
    array<vector<Base>, 2> rleSequences;
    void fillRleSequences();

    // The alignment between the two RLE sequences.
    // See seqan.hpp for its meaning.
    int64_t rleEditDistance = invalid<int64_t>;
    vector< pair<bool, bool> > rleAlignment;
    void computeRleAlignment();

    // The number of mismatches in the RLE alignment.
    uint64_t mismatchCountRle = invalid<uint64_t>;

    // The number of mismatches in the raw alignment.
    uint64_t mismatchCount = invalid<uint64_t>;

    // Total bases involved in indels (insertions + deletions).
    uint64_t indelBaseCount = invalid<uint64_t>;

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
        QuickRaw,   // Only do raw alignments, only store segments where the two raw sequences differ.
        QuickRawSparse, // Only do raw alignments, store sparse diffs (no full alignment trace).
        None,       // Construct nothing; caller fills the object explicitly
                    // (e.g. constructFromHifiasmCigar). Used to reuse hifiasm's
                    // base alignment instead of recomputing it.
    };

    ProjectedAlignment(
        const Assembler&,
        const array<OrientedReadId, 2>&,
        const Alignment&,
        Method method,
        int64_t dpMatchScore,
        int64_t dpMismatchScore,
        int64_t dpGapOpen1,
        int64_t dpGapExtend1);

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
        OverlapCigarStore* cigarStore = nullptr);

    void constructAll();
    void constructQuickRaw();
    void constructQuickRawSparse();

    // Fill the sparse diffs, statistics and CIGAR from a pre-computed hifiasm
    // CIGAR (already normalized into dinara's read0/read1 canonical frame),
    // clipped to the marker interval [read0Begin, read0End) on read0. This
    // reproduces exactly what constructQuickRawSparse would have written, but
    // walks hifiasm's tokens instead of re-running A*PA2 per segment.
    //
    //   normalizedTokens : CIGAR in read0/read1 canonical frame (op I consumes
    //                      read0, D consumes read1)
    //   read0Anchor      : read0 forward position the tokens start at
    //   read1Anchor      : read1 alignment-orientation position the tokens start at
    //   read0Begin/End   : marker interval on read0 (forward coords) to keep
    //
    // Returns false (leaving the object untouched) when the CIGAR does not span
    // the marker interval, so the caller can fall back to recomputation.
    bool constructFromHifiasmCigar(
        span<const CigarToken> normalizedTokens,
        uint32_t read0Anchor,
        uint32_t read1Anchor,
        uint32_t read0Begin,
        uint32_t read0End);



    // Marker length and its half.
    uint32_t k;
    uint32_t kHalf;

    // Single-affine scoring for overlap DP score (matches hifiasm's gapo/gape).
    const int64_t dpMatchScore;
    const int64_t dpMismatchScore;
    const int64_t dpGapOpen1;
    const int64_t dpGapExtend1;

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

    // Total bases involved in indels (insertions + deletions).
    uint64_t totalIndelBaseCount;

    // The number of gap events (indels) in the raw alignment.
    uint64_t totalGapEventCount;

    // Base-level DP score for the entire projected alignment (sum of segment scores).
    int64_t totalDpScore = invalid<int64_t>;

    // Flag for large indel (>= 6 bases)
    bool hasLargeIndel = false;

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

    // When non-null, constructQuickRawSparse() packs the full per-overlap
    // CIGAR into this store (hifiasm-style uint16_t tokens).
    OverlapCigarStore* cigarStore = nullptr;

    // CIGAR location in the store's arena, populated when cigarStore is non-null.
    uint32_t cigarOffset     = uint32_t(-1);
    uint32_t cigarTokenCount = uint32_t(-1);

    // CIGAR boundary positions (oriented-read coordinates).
    // Set during constructQuickRawSparse.
    uint32_t cigarRead0Start = 0;
    uint32_t cigarRead1Start = 0;
    uint32_t cigarRead0End   = 0;
    uint32_t cigarRead1End   = 0;
};
