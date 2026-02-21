#pragma once

// Shasta2.
#include "Base.hpp"
#include "invalid.hpp"
#include "MarkerKmers.hpp"
#include "ReadId.hpp"
#include "Shasta2Anchors.hpp"

// Standard library.
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "vector.hpp"



namespace dinara {
    class Shasta2LocalAssembly4;

    class Shasta2Anchors;
    class Shasta2AnchorPair;
}




class dinara::Shasta2LocalAssembly4 {
public:

    // This assembles between anchorIdA and anchorIdB of the given Shasta2AnchorPair.
    // This can use all of the OrientedReadIds in the Shasta2AnchorPair
    // and/or in the additionalOrientedReadIds that also appear in the
    // Kmers corresponding to the left and right Anchors.
    // The additionalOrientedReadIds must be sorted.
    // The additionalOrientedReadIds are allowed to contain OrientedReadIds
    // that are also in the Shasta2AnchorPair.
    Shasta2LocalAssembly4(
        const Shasta2Anchors&,
        uint64_t abpoaMaxLength,
        ostream& html,
        bool debug,
        const Shasta2AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds);

    // Assembled sequence and its coverage.
    vector<Base> sequence;
    vector<uint64_t> coverage;

private:

    const Shasta2Anchors& anchors;
    ostream& html;

    // The two anchors of the Shasta2AnchorPair used for this assembly.
    Shasta2AnchorId leftAnchorId;
    Shasta2AnchorId rightAnchorId;

    // MarkerInfos for the Marker Kmers corresponding to the left
    // and right anchors of the Shasta2AnchorPair being assembled.
    vector<MarkerKmers::MarkerInfo> leftMarkerInfos;
    vector<MarkerKmers::MarkerInfo> rightMarkerInfos;
    void fillMarkerInfos();

    // The union of the OrientedReadIds in the Shasta2AnchorPair and
    // the additionalOrientedReadIds.
    vector<OrientedReadId> allOrientedReadIds;
    void gatherAllOrientedReads(
        const Shasta2AnchorPair&,
        const vector<OrientedReadId>& additionalOrientedReadIds);



    // Information on OrientedReadIds that are in
    // the usableOrientedReadIds and also appear on both
    // the left and right anchors.
    class CommonOrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // The ordinals of this OrientedReadId in the left/right anchor.
        uint32_t leftOrdinal = invalid<uint32_t>;
        uint32_t rightOrdinal = invalid<uint32_t>;
        uint32_t ordinalOffset() const
        {
            return rightOrdinal - leftOrdinal;
        }

        // The base positions of this OrientedReadId's marker
        // in the left/right anchor.
        // These are mid positions of the marker in the oriented read sequence.
        uint32_t leftPosition = invalid<uint32_t>;
        uint32_t rightPosition = invalid<uint32_t>;
        uint32_t positionOffset() const
        {
            return rightPosition - leftPosition;
        }

    };
    vector<CommonOrientedReadInfo> commonOrientedReadInfos;
    void gatherCommonOrientedReads();

    // Get the sequence a CommonOrientedReadInfo contributes to the assembly.
    void getSequence(const CommonOrientedReadInfo&, vector<Base>&) const;

    // This does the multiple sequence alignment and stores sequence and coverage.
    void assemble();

    // Html output.
    void writeInput(
        const Shasta2AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds) const;
    void writeAllOrientedReadIds() const;
    void writeCommonOrientedReads() const;
    void writeAssembledSequence() const;
};
