#ifndef DINARA_ALIGNMENT_HPP
#define DINARA_ALIGNMENT_HPP

#include "invalid.hpp"
#include "OrientedReadPair.hpp"
#include "ReadId.hpp"

#include "algorithm.hpp"
#include "array.hpp"
#include <cmath>
#include "cstdint.hpp"
#include <stdexcept>
#include "utility.hpp"
#include "vector.hpp"

namespace dinara {

    class Alignment;
    class AlignmentData;
    class AlignmentInfo;
    enum class AlignmentType;
    void reverse(AlignmentType&);

    inline uint32_t computeOverlappingMarkerCount(
        uint32_t markerCount0,
        uint32_t markerCount1,
        int32_t ordinalOffset
    );

    // Classification of the overlap
    enum class CisTransStatus : uint8_t {
        Unknown = 0,
        Cis = 1,
        Trans = 2
    };

    class CompressedMarker;

    namespace MemoryMapped {
        template<class T, class Int> class VectorOfVectors;
    }
}



class dinara::Alignment {
public:

    // The ordinals in each of the two oriented reads of the
    // markers in the alignment.
    vector< array<uint32_t, 2> > ordinals;
    
    // Extended coordinates (computed during chaining for efficiency).
    // These are pre-computed to avoid duplicate marker lookups.
    // Format: half-open intervals [start, end)
    uint32_t qs = 0;  // Query start (extended)
    uint32_t qe = 0;  // Query end (extended)
    uint32_t ts = 0;  // Target start (extended, FORWARD strand for both same/diff)
    uint32_t te = 0;  // Target end (extended, FORWARD strand for both same/diff)
    
    void clear() {
        ordinals.clear();
        qs = qe = ts = te = 0;
    }

    uint32_t maxSkip() const;
    uint32_t maxDrift() const;

    void swap();
    void reverseComplement(uint32_t markerCount0, uint32_t markerCount1);

    void checkStrictlyIncreasing() const;

};



// Enum used to classify an alignment.
enum class dinara::AlignmentType {
    read0IsContained,   // 0 is contained in 1. Draw as 0tee--1.
    read1IsContained,   // 1 is contained in 0. Draw as 1tee--0.
    read0IsBackward,    // No containement, 0 is backward of 1 at both ends. Draw as 0->1.
    read1IsBackward,    // No containement, 1 is backward of 0 at both ends. Draw as 1->0.
    ambiguous           // Draw as 0diamond--diamond1
};
inline void dinara::reverse(AlignmentType& alignmentType)
{
    switch(alignmentType) {
    case AlignmentType::read0IsContained:
        alignmentType = AlignmentType::read1IsContained;
        return;
    case AlignmentType::read1IsContained:
        alignmentType = AlignmentType::read0IsContained;
        return;
    case AlignmentType::read0IsBackward:
        alignmentType = AlignmentType::read1IsBackward;
        return;
    case AlignmentType::read1IsBackward:
        alignmentType = AlignmentType::read0IsBackward;
        return;
    case AlignmentType::ambiguous:
        return;
    default:
        DINARA_ASSERT(0);
    }
}



class dinara::AlignmentInfo {
public:

    // Alignment information for each of the oriented reads in the alignment.
    class Data {
    public:

        Data() : markerCount(0) {}
        Data(
            uint32_t markerCount,
            uint32_t firstOrdinal,
            uint32_t lastOrdinal) :
            markerCount(markerCount),
            firstOrdinal(firstOrdinal),
            lastOrdinal(lastOrdinal) {}


        // Return the number of markers in this oriented read
        // that are on the left of the alignment, before the alignment begins.
        uint32_t leftTrim() const
        {
            return firstOrdinal;
        }

        // Return the number of markers in this oriented read
        // that are on the right of the alignment, after the alignment ends.
        uint32_t rightTrim() const
        {
            return markerCount - 1 - lastOrdinal;
        }

        // Return the number of markers in the range covered by the alignment.
        uint32_t range() const
        {
            return lastOrdinal + 1 - firstOrdinal;
        }

        // Return the number of bases in the range covered by the alignment.
        uint32_t baseRange(
            uint64_t k,
            OrientedReadId,
            const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) const;

        // Sanity check.
        void check() const
        {
            // Use <= to allow for case where the alignment has no markers,
            // in which case everyting is set to zero.
            DINARA_ASSERT(firstOrdinal <= markerCount);
            DINARA_ASSERT(lastOrdinal <= markerCount);
        }

        // Update to reflect reverse complementing of the oriented read.
        void reverseComplement()
        {
            std::swap(firstOrdinal, lastOrdinal);
            firstOrdinal = markerCount - 1 - firstOrdinal;
            lastOrdinal  = markerCount - 1 - lastOrdinal;
        }

        double alignmentCenter() const
        {
            return double(firstOrdinal + lastOrdinal) / 2.;
        }
        uint32_t twiceAlignmentCenter() const
        {
            return firstOrdinal + lastOrdinal;
        }

        double centerPosition() const
        {
            return double(markerCount) / 2.;
        }
        uint32_t twiceCenterPosition() const
        {
            return markerCount;
        }

        // The total number of markers in this oriented read.
        uint32_t markerCount;

        // The ordinal of the first and last marker of this oriented read
        // involved in the alignment.
        uint32_t firstOrdinal;
        uint32_t lastOrdinal;

    };
    array<Data, 2> data;



    // The number of markers in the alignment.
    // This is the same for both oriented reads!
    // It is guaranteed to never be zero for a valid alignment.
    uint32_t markerCount;

    // The minimum, maximum, and average ordinal offset,
    // computed over all aligned markers.
    int32_t minOrdinalOffset;
    int32_t maxOrdinalOffset;
    int32_t averageOrdinalOffset;

    // The maximum ordinal skip between successive aligned markers.
    uint32_t maxSkip;

    // Maximum drift is the maximum absolute value of diagonal shift
    // between successive aligned markers.
    // That, is the maximum absolute value of
    // (ordinals[i][0]-ordinals[i][1]) - (ordinals[i-1][0]-ordinals[i-1][1])
    // over the ordinals of the alignment.
    uint32_t maxDrift;

    // Flag that is set if this alignment is used in the read graph.
    uint8_t isInReadGraph : 1;

    // Uniqueness metric (alignment method 5 only).
    // See Assembler::alignOrientedReads5.
    float uniquenessMetric = std::numeric_limits<float>::signaling_NaN();

    // ProjectedAlignment metrics.
    // Only computed for read graph creation method 4 and 5.
    float errorRateRle = invalid<float>;
    uint32_t mismatchCountRle = invalid<uint32_t>;
    float errorRate = invalid<float>;
    uint32_t mismatchCount = invalid<uint32_t>;
    float errorRateGaps = invalid<float>;
    uint32_t gapCount = invalid<uint32_t>;      // Total gap BASES

    uint32_t gapEventCount = invalid<uint32_t>; // Total gap EVENTS

    // Base-level DP score computed from the base CIGAR using a two-piece affine model
    // (hifiasm/minimap2 HiFi defaults). Populated when ProjectedAlignment is computed.
    int64_t dpScore = invalid<int64_t>;

    // Evidence ID (APES/TASSD index)
    size_t alignmentId = invalid<size_t>;

    void clearFlags()
    {
        isInReadGraph = 0;
    }

    // Approximate overlap score used for hifiasm-style overlap ranking/filtering.
    // This matches the heuristic used elsewhere in Dinara:
    //   score ~= 2*overlapLen - 6*mismatch - 4*gapBases - 4*gapEvents
    // If mismatch/gap metrics are not available, this falls back to 2*overlapLen.
    int64_t hifiasmApproxScore(uint64_t overlapLen) const
    {
        int64_t score = 2 * int64_t(overlapLen);
        const bool hasMetrics =
            (mismatchCount != invalid<uint32_t>) &&
            (gapCount != invalid<uint32_t>);
        if(hasMetrics) {
            score -= 6 * int64_t(mismatchCount);
            score -= 4 * int64_t(gapCount);
            if(gapEventCount != invalid<uint32_t>) {
                score -= 4 * int64_t(gapEventCount);
            }
        }
        return score;
    }

    // Return the stored DP score (computed from the base CIGAR using a two-piece affine model).
    // This is required to be populated anywhere we use hifiasm-style overlap ranking.
    int64_t hifiasmDpScoreOrApprox(uint64_t overlapLen) const
    {
        (void)overlapLen;
        if(dpScore == invalid<int64_t>) {
            throw std::runtime_error(
                "AlignmentInfo::dpScore is not populated. Regenerate alignments with ProjectedAlignment DP scoring enabled.");
        }
        return dpScore;
    }



    // Constructors.
    AlignmentInfo(
        const Alignment& alignment,
        const array<uint32_t, 2>& markerCounts)
    {
        create(alignment, markerCounts);
        clearFlags();
    }
    AlignmentInfo(
        const Alignment& alignment,
        uint32_t markerCount0,
        uint32_t markerCount1)
    {
        create(alignment, array<uint32_t, 2>({markerCount0, markerCount1}));
        clearFlags();
    }
    void create(
        const Alignment&,
        const array<uint32_t, 2>& markerCounts);
    void create(
        const Alignment&,
        uint32_t markerCount0,
        uint32_t markerCount1);
    AlignmentInfo() : markerCount(0)
    {
        clearFlags();
    }



    // Update to reflect a swap the two oriented reads.
    void swap()
    {
        std::swap(data[0], data[1]);
        minOrdinalOffset = -minOrdinalOffset;
        maxOrdinalOffset = -maxOrdinalOffset;
        averageOrdinalOffset = -averageOrdinalOffset;
    }

    // Update to reflect reverse complementing of the two oriented reads.
    void reverseComplement()
    {
        for(size_t i=0; i<2; i++) {
            data[i].reverseComplement();
        }
        const int32_t delta = int32_t(data[0].markerCount) - int32_t(data[1].markerCount);
        minOrdinalOffset =  delta - minOrdinalOffset;
        maxOrdinalOffset =  delta - maxOrdinalOffset;
        averageOrdinalOffset =  delta - averageOrdinalOffset;
    }

    // Some accessors.
    uint32_t leftTrim(size_t i) const {
        return data[i].leftTrim();
    }
    uint32_t rightTrim(size_t i) const {
        return data[i].rightTrim();
    }
    uint32_t range(size_t i) const {    // In markers
        return data[i].range();
    }
    uint32_t baseRange(
        uint64_t k,
        OrientedReadId orientedReadId,
        size_t i,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) const {    // In bases
        return data[i].baseRange(k, orientedReadId, markers);
    }

    // Return the ratio of aligned markers over the alignment range.
    double alignedFraction(size_t i) const
    {
        return double(markerCount) / double(range(i));
    }
    double minAlignedFraction() const
    {
        return min(alignedFraction(0), alignedFraction(1));
    }

    // Compute the left and right trim, expressed in markers.
    // This is the minimum number of markers (over the two oriented reads)
    // that are excluded from the alignment on each side.
    // If the trim is too high, the alignment is suspicious.
    pair<uint32_t, uint32_t> computeTrim() const
    {
        const uint32_t leftTrim  = min(data[0].leftTrim() , data[1].leftTrim() );
        const uint32_t rightTrim = min(data[0].rightTrim(), data[1].rightTrim());
        return make_pair(leftTrim, rightTrim);
    }

    // Find out if this is a containing alignment,
    // that is, if the alignment covers one read
    // entirely, except possibly for up to maxTim
    // markers on each side.
    bool isContaining(uint32_t maxTrim) const {
        for(size_t i=0; i<2; i++) {
            if(leftTrim(i)<=maxTrim && rightTrim(i)<=maxTrim) {
                return true;
            }
        }
        return false;
    }



    // Return true if this alignment has read i unambiguously contained
    // in read 1-i.
    bool isContained(int i, uint32_t maxTrim) const {

        // Figure out if the two reads are fully covered by the alignment
        // (except possibly up to maxTrim markers at eah end).
        array<bool, 2> coversFullRead;
        coversFullRead[0] = (leftTrim(0)<=maxTrim && rightTrim(0)<=maxTrim);
        coversFullRead[1] = (leftTrim(1)<=maxTrim && rightTrim(1)<=maxTrim);

        // Return true only if the first read is fully covered
        // and the second read is not.
        return coversFullRead[i] and not coversFullRead[1-i];
    }



    // Jaccard similarity of the alignment.
    double jaccard() const
    {
        const uint32_t range0 = range(0);
        const uint32_t range1 = range(1);
        const double intersectionRange = double(range0 + range1) / 2.;
        const double unionRange = double(data[0].markerCount + data[1].markerCount) - intersectionRange;
        return intersectionRange / unionRange;
    }


    // Return the offset between the centers of the two oriented reads,
    // as estimated form the alignment.
    // The offset is positive if the center of the second read
    // is to the right of the center of the first read.
    double offsetAtCenter() const
    {
        const double alignmentCenter0 = data[0].alignmentCenter();
        const double alignmentCenter1 = data[1].alignmentCenter();

        const double center0 = data[0].centerPosition();
        const double center1 = data[1].centerPosition();

        return
            (center1 - alignmentCenter1) -
            (center0 - alignmentCenter0);
    }
    int32_t twiceOffsetAtCenter() const
    {
        const uint32_t twiceAlignmentCenter0 = data[0].twiceAlignmentCenter();
        const uint32_t twiceAlignmentCenter1 = data[1].twiceAlignmentCenter();

        const uint32_t twiceCenter0 = data[0].twiceCenterPosition();
        const uint32_t twiceCenter1 = data[1].twiceCenterPosition();

        return
            (int(twiceCenter1) - int(twiceAlignmentCenter1)) -
            (int(twiceCenter0) - int(twiceAlignmentCenter0));
    }


    // Classify this alignment.
    AlignmentType classify(uint32_t maxTrim) const
    {
        // Compute trim.
        const uint32_t leftTrim0  = leftTrim(0);
        const uint32_t leftTrim1  = leftTrim(1);
        const uint32_t rightTrim0 = rightTrim(0);
        const uint32_t rightTrim1 = rightTrim(1);

        // Check for containment.
        const bool isContained0 = (leftTrim0<=maxTrim) && (rightTrim0<=maxTrim);
        const bool isContained1 = (leftTrim1<=maxTrim) && (rightTrim1<=maxTrim);
        if(isContained0 && !isContained1) {
            // 0 is unambiguously contained in 1.
            return AlignmentType::read0IsContained;
        }
        if(isContained1 && !isContained0) {
            // 1 is unambiguously contained in 0.
            return AlignmentType::read1IsContained;
        }
        if(isContained0 && isContained1) {
            // Near complete overlap.
            return AlignmentType::ambiguous;
        }

        // If getting here, no containment found.
        DINARA_ASSERT(!isContained0 && !isContained1);

        // Figure out if one of the two reads is backward at both ends.
        const bool read0IsBackward =
            leftTrim0>maxTrim  && rightTrim0<=maxTrim &&
            leftTrim1<=maxTrim && rightTrim1>=maxTrim;
        const bool read1IsBackward =
            leftTrim1>maxTrim  && rightTrim1<=maxTrim &&
            leftTrim0<=maxTrim && rightTrim0>=maxTrim;
        if(read0IsBackward && !read1IsBackward) {
            return AlignmentType::read0IsBackward;
        }
        if(read1IsBackward && !read0IsBackward) {
            return AlignmentType::read1IsBackward;
        }
        return AlignmentType::ambiguous;
    }

    void write(ostream& s) const
    {
        s << "Alignment with " << markerCount << " aligned markers:\n";
        for(int i=0; i<2; i++) {
            const auto& d = data[i];
            s << "    " << i;
            s << ": first " << d.firstOrdinal;
            s << ", last: " << d.lastOrdinal;
            s << ", total: " << d.markerCount << "\n";
        }
        s << "    Offset at center: " << offsetAtCenter() << endl;
    }
};



class dinara::AlignmentData :
    public dinara::OrientedReadPair {
public:

    using DeleteReasonMask = uint16_t;
    static constexpr DeleteReasonMask DeleteReasonNone        = 0;
    static constexpr DeleteReasonMask DeleteReasonPhase       = 1u << 0; // EC phasing decision
    static constexpr DeleteReasonMask DeleteReasonSecondary   = 1u << 1; // Redundant/secondary per-read-pair filtering
    static constexpr DeleteReasonMask DeleteReasonChemical    = 1u << 2; // ONT chemical chimera masking
    static constexpr DeleteReasonMask DeleteReasonChimeric    = 1u << 3; // Chimeric read filtering
    static constexpr DeleteReasonMask DeleteReasonLocal       = 1u << 4; // ma_hit_sub / local segment filtering
    static constexpr DeleteReasonMask DeleteReasonCoverageCut = 1u << 5; // ma_hit_cut / coverage cuts
    static constexpr DeleteReasonMask DeleteReasonHanging     = 1u << 6; // ma_hit_flt / hanging overlap filter
    static constexpr DeleteReasonMask DeleteReasonContained   = 1u << 7; // ma_hit_contained_advance / contained read removal
    static constexpr DeleteReasonMask DeleteReasonPalindromic = 1u << 8; // palindromic read filtering
    static constexpr DeleteReasonMask DeleteReasonContainedPrune = 1u << 9; // keep best overlap for contained reads

    // The AlignmentInfo computed with the first read on strand 0.
    AlignmentInfo info;

    // Classification of the overlap
    CisTransStatus cisTransStatus = CisTransStatus::Unknown;

    // Flags
    // Number of informative het/SV sites (in query coordinates) covered by this overlap,
    // as computed by performHifiasmECParity, stored separately for each read perspective.
    // informativeHetSiteCount0 corresponds to readIds[0]'s view, informativeHetSiteCount1 to readIds[1]'s view.
    uint32_t informativeHetSiteCount0 = 0;
    uint32_t informativeHetSiteCount1 = 0;
    // Informative-site score for this overlap.
    // This is computed as max(informativeHetSiteCount0, informativeHetSiteCount1) after both sides
    // have been populated by performHifiasmECParity.
    uint32_t informativeHetSiteScore = 0;
    bool coversHetSiteAtLeast(uint32_t minCount) const
    {
        return informativeHetSiteCount0 >= minCount || informativeHetSiteCount1 >= minCount;
    }
    bool coversHetSite() const
    {
        return coversHetSiteAtLeast(1);
    }

    uint32_t getInformativeHetSiteCountFromReadPerspective(ReadId queryReadId) const
    {
        if(readIds[0] == queryReadId) {
            return informativeHetSiteCount0;
        } else if(readIds[1] == queryReadId) {
            return informativeHetSiteCount1;
        } else {
            return 0;
        }
    }

    void updateInformativeHetSiteScore()
    {
        informativeHetSiteScore = std::max(informativeHetSiteCount0, informativeHetSiteCount1);
    }
    
    // Directional deletion reason bitmasks.
    // deleteReasons0: reasons from readIds[0]'s perspective
    // deleteReasons1: reasons from readIds[1]'s perspective
    // An overlap is only kept for graph purposes if BOTH sides have no deletion reasons (conservative AND).
    DeleteReasonMask deleteReasons0 = DeleteReasonNone;
    DeleteReasonMask deleteReasons1 = DeleteReasonNone;

    bool isDeleted0() const { return deleteReasons0 != DeleteReasonNone; }
    bool isDeleted1() const { return deleteReasons1 != DeleteReasonNone; }
    bool keptByBothSides() const { return !isDeleted0() && !isDeleted1(); }
    bool isDeleted() const { return isDeleted0() && isDeleted1(); } // fully deleted (both sides)
    bool isCisByBothSides() const
    {
        return
            ((deleteReasons0 & DeleteReasonPhase) == 0) &&
            ((deleteReasons1 & DeleteReasonPhase) == 0);
    }

    void addDeleteReasonsBoth(DeleteReasonMask reasons)
    {
        deleteReasons0 |= reasons;
        deleteReasons1 |= reasons;
    }
    void addDeleteReasonsFromReadPerspective(ReadId readId, DeleteReasonMask reasons)
    {
        if (readIds[0] == readId) {
            deleteReasons0 |= reasons;
        } else {
            deleteReasons1 |= reasons;
        }
    }
    void clearDeleteReasonsFromReadPerspective(ReadId readId, DeleteReasonMask reasons)
    {
        if (readIds[0] == readId) {
            deleteReasons0 &= ~reasons;
        } else {
            deleteReasons1 &= ~reasons;
        }
    }
    void clearDeleteReasonsBoth(DeleteReasonMask reasons)
    {
        deleteReasons0 &= ~reasons;
        deleteReasons1 &= ~reasons;
    }
    
    bool hasLargeIndel = false;

    // Explicit Coordinates (likely in bases, derived from markers or passed from PAF)
    uint32_t qs = 0; // Query Start
    uint32_t qe = 0; // Query End
    uint32_t ts = 0; // Target Start
    uint32_t te = 0; // Target End

    AlignmentData() {}
    AlignmentData(
        const array<ReadId, 2>& readIds,
        bool isSameStrand,
        const AlignmentInfo& info) :
        OrientedReadPair(readIds, isSameStrand),
        info(info)
    {}
    AlignmentData(
        const OrientedReadPair& orientedReadPair,
        const AlignmentInfo& info) :
        OrientedReadPair(orientedReadPair),
        info(info)
    {}

    // Given an AlignmentData, return its AlignmentInfo,
    // after swapping and/or reverse complementing it
    // to make sure it refers to the given OrientedReadId's,
    // in that order.
    AlignmentInfo orient(OrientedReadId, OrientedReadId) const;

};



// Compute the number of overlapping markers between two
// oriented reads with given number of markers
// and ordinal offset.
inline uint32_t dinara::computeOverlappingMarkerCount(
    uint32_t markerCount0,
    uint32_t markerCount1,
    int32_t ordinalOffset)
{
    const int32_t begin0 = 0;
    const int32_t end0 = int32_t(markerCount0);
    const int32_t begin1 = ordinalOffset;
    const int32_t end1 = begin1 + int32_t(markerCount1);

    const int32_t overlapBegin = max(begin0, begin1);
    const int32_t overlapEnd = min(end0, end1);

    if(overlapEnd < overlapBegin) {
        return 0;
    } else {
        return overlapEnd - overlapBegin;
    }

}

#endif
