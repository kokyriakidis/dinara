// Dinara.
#include "Shasta2LocalAssembly4.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2AnchorPair.hpp"
#include "Shasta2LocalAssembly5.hpp"
#include "Marker.hpp"
#include "Shasta2PoastaWrapper.hpp"
#include "ReadId.hpp"
#include "Reads.hpp"
using namespace dinara;

// Standard library.
#include "algorithm.hpp"
#include <map>



Shasta2LocalAssembly4::Shasta2LocalAssembly4(
    const Shasta2Anchors& anchors,
    uint64_t abpoaMaxLength,
    ostream& html,
    bool debug,
    const Shasta2AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) :
    anchors(anchors),
    html(html),
    leftAnchorId(anchorPair.anchorIdA),
    rightAnchorId(anchorPair.anchorIdB)
{
    // EXPOSE WHEN CODE STABILIZES.
    // const uint64_t minCoverage = 3;

    DINARA_ASSERT(std::ranges::is_sorted(additionalOrientedReadIds));
    if(html) {
        writeInput(anchorPair, additionalOrientedReadIds);
    }

    gatherAllOrientedReads(anchorPair, additionalOrientedReadIds);
    if(html) {
        writeAllOrientedReadIds();
    }

    fillMarkerInfos();

    gatherCommonOrientedReads();
    if(html) {
        writeCommonOrientedReads();
    }



    // If coverage is too low, use Shasta2LocalAssembly5 instead,
    // which can use OrientedReadIds that appear only on the left
    // or only on the right.
    if(true /* commonOrientedReadInfos.size() < minCoverage */) {
        if(html) {
            html << "<br><br>Switching to Shasta2LocalAssembly5 due to low coverage<hr>";
        }
        Shasta2LocalAssembly5 localAssembly5(
            anchors,
            abpoaMaxLength,
            html,
            debug,
            anchorPair,
            additionalOrientedReadIds);

        // For now ignore the output of Shasta2LocalAssembly5 and continue with Shasta2LocalAssembly4.
        sequence = localAssembly5.sequence;
        coverage = localAssembly5.coverage;
        return;
    }



    assemble();
    if(html) {
        writeAssembledSequence();
    }
}



void Shasta2LocalAssembly4::gatherAllOrientedReads(
    const Shasta2AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds)
{
    std::ranges::set_union(anchorPair.orientedReadIds, additionalOrientedReadIds,
        back_inserter(allOrientedReadIds));

}



void Shasta2LocalAssembly4::fillMarkerInfos()
{
    // Get the Kmers.
    const Kmer leftKmer = anchors.anchorKmer(leftAnchorId);
    const Kmer rightKmer = anchors.anchorKmer(rightAnchorId);

    // Get the MarkerInfos.
    anchors.getMarkerKmers()->get(leftKmer, leftMarkerInfos);
    anchors.getMarkerKmers()->get(rightKmer, rightMarkerInfos);
}



// The CommonOrientedReads are the ones that will be used
// in the assembly. They are the intersection of
// allOrientedReadIds with the OrientedReadIds that
// appear in both the left and right MarkerInfos.
void Shasta2LocalAssembly4::gatherCommonOrientedReads()
{

    const uint32_t kHalf = uint32_t(anchors.k / 2);

    // Find OrientedReadIds that appear in both the left and right markerInfos.
    // An OrientedReadId can appear only once in a marker Kmer corresponding to an Anchor.
    class MarkerInfoPair {
    public:
        OrientedReadId orientedReadId;
        uint32_t leftOrdinal;
        uint32_t rightOrdinal;
    };
    vector<MarkerInfoPair> markerInfoPairs;
    auto itLeft = leftMarkerInfos.begin();
    const auto itLeftEnd = leftMarkerInfos.end();
    auto itRight = rightMarkerInfos.begin();
    const auto itRightEnd = rightMarkerInfos.end();
    while((itLeft != itLeftEnd) and (itRight != itRightEnd)) {
        if(itLeft->orientedReadId < itRight->orientedReadId) {
            ++itLeft;
        } else if(itRight->orientedReadId < itLeft->orientedReadId) {
            ++itRight;
        } else {
            const OrientedReadId orientedReadId = itLeft->orientedReadId;
            DINARA_ASSERT(orientedReadId == itRight->orientedReadId);
            if(itLeft->ordinal < itRight->ordinal) {
                MarkerInfoPair& markerInfoPair = markerInfoPairs.emplace_back();
                markerInfoPair.orientedReadId = orientedReadId;
                markerInfoPair.leftOrdinal = itLeft->ordinal;
                markerInfoPair.rightOrdinal = itRight->ordinal;
            }
            ++itLeft;
            ++itRight;
        }
    }



    // The common OrientedReadIds must appear in both the
    // markerInfoPairs and allOrientedReadIds.
    auto itPairs = markerInfoPairs.begin();
    const auto itPairsEnd = markerInfoPairs.end();
    auto itAll = allOrientedReadIds.begin();
    const auto itAllEnd = allOrientedReadIds.end();
    while((itPairs != itPairsEnd) and (itAll != itAllEnd)) {
        if(itPairs->orientedReadId < *itAll) {
            ++itPairs;
        } else if(*itAll < itPairs->orientedReadId) {
            ++itAll;
        } else {
            const OrientedReadId orientedReadId = itPairs->orientedReadId;
            DINARA_ASSERT(orientedReadId == *itAll);

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

            CommonOrientedReadInfo& commonOrientedReadInfo = commonOrientedReadInfos.emplace_back();
            commonOrientedReadInfo.orientedReadId = orientedReadId;

            commonOrientedReadInfo.leftOrdinal = itPairs->leftOrdinal;
            commonOrientedReadInfo.rightOrdinal = itPairs->rightOrdinal;

            commonOrientedReadInfo.leftPosition = orientedReadMarkers[commonOrientedReadInfo.leftOrdinal].position + kHalf;
            commonOrientedReadInfo.rightPosition = orientedReadMarkers[commonOrientedReadInfo.rightOrdinal].position + kHalf;

            ++itPairs;
            ++itAll;
        }
    }
}



void Shasta2LocalAssembly4::assemble()
{
    DINARA_ASSERT(not commonOrientedReadInfos.empty());

    // Gather distinct sequences and their coverage.
    class DistinctSequence {
    public:
        shared_ptr< vector<Base> > sequencePointer;
        const vector<Base>& sequence() const
        {
            return *sequencePointer;
        }
        uint64_t coverage;

        // Order by decreasing coverage.
        bool operator<(const DistinctSequence& that) {
            return coverage > that.coverage;
        }
    };
    vector<DistinctSequence> distinctSequences;
    vector<Base> orientedReadSequence;
    for(const CommonOrientedReadInfo& info: commonOrientedReadInfos) {
        getSequence(info, orientedReadSequence);

        // See if this sequence is already in our DistinctSequences.
        bool done = false;
        for(DistinctSequence& distinctSequence: distinctSequences) {
            if(distinctSequence.sequence() == orientedReadSequence) {
                ++distinctSequence.coverage;
                done = true;
                break;
            }
        }

        // If we did not have it, create it with coverage 1.
        if(not done) {
            DistinctSequence& distinctSequence = distinctSequences.emplace_back();
            distinctSequence.sequencePointer = make_shared< vector<Base> >(orientedReadSequence);
            distinctSequence.coverage = 1;
        }
    }
    sort(distinctSequences.begin(), distinctSequences.end());



    // Write the DistinctSequences.
    if(html) {
        html <<
            "<h3>Distinct sequences</h3>"
            "<table><tr>"
            "<th>Coverage<th>Length<th>Sequence";
        for(const DistinctSequence& distinctSequence: distinctSequences) {
            const vector<Base>& sequence = distinctSequence.sequence();
            html <<
                "<tr>"
                "<td class=centered>" << distinctSequence.coverage <<
                "<td class=centered>" << sequence.size() <<
                "<td class=left style='font-family:monospace;white-space:nowrap''>";
            std::ranges::copy(sequence, ostream_iterator<Base>(html));
        }
        html << "</table>";
    }



    // If there is a dominant sequence, use it as the consensus.
    const uint64_t maximumCoverage = distinctSequences.front().coverage;
    const uint64_t totalCoverage = commonOrientedReadInfos.size();
    if(maximumCoverage > totalCoverage/2) {
        sequence = distinctSequences.front().sequence();
        coverage.resize(sequence.size(), maximumCoverage);
        if(html) {
            html << "<br>The dominant sequence with coverage " << maximumCoverage <<
                " was used as the consensus.";
            if(maximumCoverage < totalCoverage) {
                html << "<br>Store base coverages are lower bounds.";
            }
        }
        return;
    }



    // If getting here, we have to run the multiple sequence alignment of the distinct sequences.
    vector< pair<vector<Base>, uint64_t> > sequencesWithCoverage;
    for(const DistinctSequence& distinctSequence: distinctSequences) {
         sequencesWithCoverage.emplace_back(distinctSequence.sequence(), distinctSequence.coverage);
    }
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;
    shasta2Poasta(sequencesWithCoverage, consensus, alignment, alignedConsensus);

    // Store assembled sequence
    for(const auto& [base, baseCoverage]: consensus) {
        sequence.push_back(base);
        coverage.push_back(baseCoverage);
    }



    // Write the multiple sequence alignment of the distinct sequences.
    if(html) {
        html <<
            "<h3>Multiple sequence alignment of the distinct sequences</h3>"
            "<table><tr>"
            "<th>Coverage<th>Length<th>Aligned sequence";
        for(uint64_t i=0; i<distinctSequences.size(); i++) {
            const DistinctSequence& distinctSequence = distinctSequences[i];
            const vector<Base>& sequence = distinctSequence.sequence();
            const vector<AlignedBase>& alignmentRow = alignment[i];
            html <<
                "<tr>"
                "<td class=centered>" << distinctSequence.coverage <<
                "<td class=centered>" << sequence.size() <<
                "<td class=left style='font-family:monospace;white-space:nowrap''>";
            for(uint64_t position=0; position<alignmentRow.size(); position++) {
                const AlignedBase base = alignmentRow[position];
                const bool isMismatch = (base != alignedConsensus[position]);
                if(isMismatch) {
                    html << "<span style='background-color:Pink'>";
                }
                html << base;
                if(isMismatch) {
                    html << "</span>";
                }
            }
        }

        // Add a row with aligned consensus.
        html << "<tr><th>Consensus<td class=centered>" << consensus.size() <<
            "<td class=left style='font-family:monospace;white-space:nowrap'>";
        uint64_t consensusPosition = 0;
        for(uint64_t position=0; position<alignedConsensus.size(); position++) {
            const AlignedBase base = alignedConsensus[position];
            if(base.isGap()) {
                html << alignedConsensus[position];
            } else {
                const uint64_t coverage = consensus[consensusPosition].second;
                html << "<span title='Position " << consensusPosition <<
                    ", coverage " << coverage << "'";
                if(coverage != totalCoverage) {
                    html << " style='background-color:pink'";
                }
                html << ">";
                html << alignedConsensus[position];
                html << "</span>";
                ++consensusPosition;
            }
        }
        html << "</table>";
    }
}



// Get the sequence a CommonOrientedReadInfo contributes to the assembly.
void Shasta2LocalAssembly4::getSequence(const CommonOrientedReadInfo& info, vector<Base>& sequence) const
{
    sequence.clear();
    for(uint32_t position=info.leftPosition; position<info.rightPosition; position++) {
        sequence.push_back(anchors.reads.getOrientedReadBase(info.orientedReadId, position));
    }
}



void Shasta2LocalAssembly4::writeInput(
    const Shasta2AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) const
{
    html <<
        "<table>"
        "<tr><th class=left>Left anchor<td class=centered>" << shasta2AnchorIdToString(leftAnchorId) <<
        "<tr><th class=left>Right anchor<td class=centered>" << shasta2AnchorIdToString(rightAnchorId) <<
        "</table>";
    html << "<h3>OrientedReadIds in the Shasta2AnchorPair</h3><table>";
    for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << anchorPair.orientedReadIds.size() << " oriented reads in the Shasta2AnchorPair.";

    html << "<h3>Additional OrientedReadIds</h3>"
        "These are additional OrientedReadIds that can be used in this assembly."
        "<table>";
    for(const OrientedReadId orientedReadId: additionalOrientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << additionalOrientedReadIds.size() << " additional oriented reads.";

}


void Shasta2LocalAssembly4::writeAllOrientedReadIds() const
{
    html << "<h3>All OrientedReadIds</h3>"
        "These are the union of the OrientedReadIds in the Shasta2AnchorPair "
        "and the additional OrientedReadIds"
        "<table>";
    for(const OrientedReadId orientedReadId: allOrientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << allOrientedReadIds.size() << " oriented reads.";
}



void Shasta2LocalAssembly4::writeCommonOrientedReads() const
{
    html <<
        "<h3>Common oriented reads</h3>"
        "These are OrientedReadIds from the above list that are present in both the "
        "left and right anchors, and that visit the right marker k-mer "
        "after (not necessarily immediately after) they visit the left marker k-mer."
        "<table>"
        "<tr>"
        "<th>OrientedReadId"
        "<th>Left<br>ordinal<th>Right<br>ordinal<th>Ordinal<br>offset"
        "<th>Left<br>position<th>Right<br>position<th>Position<br>offset"
        "<th>Sequence";

    vector<Base> sequence;
    for(const CommonOrientedReadInfo& info: commonOrientedReadInfos) {
        getSequence(info, sequence);
        html <<
            "<tr>"
            "<td class=centered>" << info.orientedReadId <<
            "<td class=centered>" << info.leftOrdinal <<
            "<td class=centered>" << info.rightOrdinal <<
            "<td class=centered>" << info.ordinalOffset() <<
            "<td class=centered>" << info.leftPosition <<
            "<td class=centered>" << info.rightPosition <<
            "<td class=centered>" << info.positionOffset() <<
            "<td class=left style='font-family:monospace;white-space: nowrap''>";
        std::ranges::copy(sequence, ostream_iterator<Base>(html));
    }
    html << "</table>" << commonOrientedReadInfos.size() << " common oriented reads.";

}



void Shasta2LocalAssembly4::writeAssembledSequence() const
{
    html <<
        "<h3>Consensus</h3>"
        "<p><span style='font-family:monospace;white-space:nowrap''>"
        ">Shasta2LocalAssembly " << sequence.size() <<
        "<br>";
    for(uint64_t position=0; position<sequence.size(); position++) {
        html << "<span title='Position " << position <<
            ", coverage " << coverage[position] << "'";
        if(coverage[position] != commonOrientedReadInfos.size()) {
            html << " style='background-color:Pink'";
        }
        html << ">";
        html << sequence[position];
        html << "</span>";
    }
    html << "</span><br><br>";

    html <<
        "<table>"
        "<tr><th class=left>Consensus sequence length<td class=left>" << sequence.size() <<
        "<tr><th class=left>Consensus sequence"
        "<td style='font-family:monospace;white-space:nowrap''>";

    for(uint64_t position=0; position<sequence.size(); position++) {
        html << "<span title='Position " << position <<
            ", coverage " << coverage[position] << "'";
        if(coverage[position] != commonOrientedReadInfos.size()) {
            html << " style='background-color:Pink'";
        }
        html << ">";
        html << sequence[position];
        html << "</span>";
    }

    html <<
        "<tr><th class=left >Coverage"
        "<td style='font-family:monospace;white-space:nowrap''>";

    std::map<uint64_t, char> coverageLegend;

    for(uint64_t position=0; position<sequence.size(); position++) {
        const uint64_t cov= coverage[position];
        char coverageCharacter = '*';
        if(cov < 10) {
            coverageCharacter = char(cov) + '0';
        } else if(cov < 36) {
            coverageCharacter = (char(cov) - char(10)) + 'A';
        }
        coverageLegend[cov] = coverageCharacter;
        html << "<span title='Position " << position <<
            ", coverage " << cov << "'";
        if(cov != commonOrientedReadInfos.size()) {
            html << " style='background-color:Pink'";
        }
        html << "'>";
        html << coverageCharacter;
        html << "</span>";
    }

    html << "</table><br><br>";

    // Write the coverage legend.
    html << "<p><table><tr><th>Coverage<th>Symbol";
    for(const auto& p: coverageLegend) {
        html << "<tr><td class=centered>" << p.first << "<td class=centered>" << p.second;
    }
    html << "</table>";


}
