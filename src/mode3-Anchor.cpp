#include "mode3-Anchor.hpp"
#include "Alignment.hpp"
#include "deduplicate.hpp"
#include "html.hpp"
#include "ReadGraph.hpp"
#include "Marker.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"
using namespace dinara;
using namespace mode3;

#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <system_error>
#include <unordered_set>

namespace {
    inline mode3::AnchorId reverseComplementAnchorId(mode3::AnchorId anchorId)
    {
        return anchorId ^ 1ULL;
    }

    inline std::pair<mode3::AnchorId, mode3::AnchorId> canonicalEndpointPair(
        mode3::AnchorId startAnchorId,
        mode3::AnchorId endAnchorId)
    {
        const std::pair<mode3::AnchorId, mode3::AnchorId> p0 = {startAnchorId, endAnchorId};
        const std::pair<mode3::AnchorId, mode3::AnchorId> p1 = {
            reverseComplementAnchorId(endAnchorId),
            reverseComplementAnchorId(startAnchorId)};
        return (p1 < p0) ? p1 : p0;
    }
}

// Explicit instantiation.
#include "MultithreadedObject.tpp"
namespace dinara {
    template class MultithreadedObject<mode3::Anchors>;
}



// This constructor accesses existing Anchors.
Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    DINARA_ASSERT((k %2) == 0);
    kHalf = k / 2;

    anchorMarkerIntervals.accessExistingReadOnly(largeDataName("AnchorMarkerIntervals"));
    anchorSequences.accessExistingReadOnly(largeDataName("AnchorSequences"));
    anchorInfos.accessExistingReadWrite(largeDataName("AnchorInfos"));
    journeys.accessExistingReadOnly(largeDataName("Journeys"));
}



// This constructor creates the Anchors from explicitly provided marker intervals.
Anchors::Anchors(
    const MappedMemoryOwner& mappedMemoryOwner,
    const Reads& reads,
    uint64_t k,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const vector<vector<AnchorMarkerInterval>>& anchorsExplicit,
    uint32_t ordinalOffset,
    uint64_t /* threadCount */) :
    MultithreadedObject<Anchors>(*this),
    MappedMemoryOwner(mappedMemoryOwner),
    reads(reads),
    k(k),
    markers(markers)
{
    performanceLog << timestamp << "Anchor creation from explicit marker intervals begins." << endl;

    DINARA_ASSERT((k % 2) == 0);
    kHalf = k / 2;
    DINARA_ASSERT(reads.representation == 0);
    DINARA_ASSERT(markers.isOpen());

    anchorMarkerIntervals.createNew(
        largeDataName("AnchorMarkerIntervals"),
        largeDataPageSize);
    anchorSequences.createNew(
        largeDataName("AnchorSequences"),
        largeDataPageSize);
    anchorInfos.createNew(
        largeDataName("AnchorInfos"),
        largeDataPageSize);

    for(const auto& anchor : anchorsExplicit) {
        anchorMarkerIntervals.appendVector();
        for(const auto& interval : anchor) {
            anchorMarkerIntervals.append(interval);
        }
        // No explicit sequence for this anchor.
        anchorSequences.appendVector();
    }

    anchorInfos.resize(anchorMarkerIntervals.size());
    for(AnchorId anchorId=0; anchorId<anchorInfos.size(); ++anchorId) {
        anchorInfos[anchorId].ordinalOffset = ordinalOffset;
        anchorInfos[anchorId].componentId = invalid<uint32_t>;
        anchorInfos[anchorId].localAnchorIdInComponent = invalid<uint64_t>;
    }

    performanceLog << timestamp << "Anchor creation from explicit marker intervals ends (" << anchorInfos.size() << " anchors)." << endl;
}



Anchor Anchors::operator[](AnchorId anchorId) const
{
    return anchorMarkerIntervals[anchorId];
}



// This returns the sequence between the midpoint of the first marker of the
// anchor and the midpoint of the second marker.
// When anchorCreationMethod is "FromMarkerKmers", this sequence is empty.
span<const Base> Anchors::anchorSequence(AnchorId anchorId) const
{
    return anchorSequences[anchorId];
}



// This returns the sequence between the beginning of the first marker of the
// anchor and the end of the second marker.
// Get the sequence from the first oriented read of the anchor
// (all oriented reads are guaranteed to have the same anchor sequence).
vector<Base> Anchors::anchorExtendedSequence(AnchorId anchorId) const
{
    // Get the first AnchorMarkerInterval for this Anchor.
    const Anchor anchor = (*this)[anchorId];
    const AnchorMarkerInterval& firstMarkerInterval = anchor.front();

    // Get the OrientedReadId and the ordinals.
    const OrientedReadId orientedReadId = firstMarkerInterval.orientedReadId;
    const uint32_t ordinal0 = firstMarkerInterval.ordinal0;
    const uint32_t ordinal1 = ordinal0 + ordinalOffset(anchorId);

    // Access the markers of this OrientedReadId.
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];

    const CompressedMarker& marker0 = orientedReadMarkers[ordinal0];
    const CompressedMarker& marker1 = orientedReadMarkers[ordinal1];

    const uint32_t begin = marker0.position;
    const uint32_t end = marker1.position + uint32_t(k);

    vector<Base> sequence;
    for(uint32_t position=begin; position!=end; position++) {
        sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
    }

    return sequence;
}



uint64_t Anchors::size() const
{
    return anchorMarkerIntervals.size();
}



void Anchors::check() const
{
    const Anchors& anchors = *this;

    for(AnchorId anchorId=0; anchorId<size(); anchorId++) {
        const Anchor& anchor = anchors[anchorId];
        anchor.check();
    }
}



void Anchor::check() const
{
    const Anchor& anchor = *this;

    for(uint64_t i=1; i<size(); i++) {
        DINARA_ASSERT(anchor[i-1].orientedReadId.getReadId() < anchor[i].orientedReadId.getReadId());
    }
}



// Return the number of common oriented reads between two Anchors.
uint64_t Anchors::countCommon(
    AnchorId anchorId0,
    AnchorId anchorId1,
    bool ignoreNegativeOffsets) const
{
    const Anchors& anchors = *this;
    const Anchor anchor0 = anchors[anchorId0];
    const Anchor anchor1 = anchors[anchorId1];

    return anchor0.countCommon(anchor1, ignoreNegativeOffsets);
}



// Return the number of common oriented reads with another Anchor.
// Oriented reads in each Anchor are sorted and not duplicated.
uint64_t Anchor::countCommon(const Anchor& that, bool ignoreNegativeOffsets) const
{
    const Anchor& anchor0 = *this;
    const Anchor& anchor1 = that;

    auto it0 = anchor0.begin();
    auto it1 = anchor1.begin();

    const auto end0 = anchor0.end();
    const auto end1 = anchor1.end();

    uint64_t count = 0;
    while((it0 != end0) and (it1 != end1)) {
        const OrientedReadId orientedReadId0 = it0->orientedReadId;
        const OrientedReadId orientedReadId1 = it1->orientedReadId;
        if(orientedReadId0 < orientedReadId1) {
            ++it0;
        } else if(orientedReadId1 < orientedReadId0) {
            ++it1;
        } else {
            if(ignoreNegativeOffsets){
                if(it0->ordinal0 < it1->ordinal0) {
                    ++count;
                }
            } else {
                ++count;
            }
            ++it0;
            ++it1;
        }
    }

    return count;
}



void Anchors::analyzeAnchorPair(
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    AnchorPairInfo& info
    ) const
{
    const Anchors& anchors = *this;

    // Prepare for the joint loop over OrientedReadIds of the two Anchors.
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];
    const auto beginA = anchorA.begin();
    const auto beginB = anchorB.begin();
    const auto endA = anchorA.end();
    const auto endB = anchorB.end();

    // Store the total number of OrientedReadIds on the two edges.
    info.totalA = endA - beginA;
    info.totalB = endB - beginB;


    // Joint loop over the MarkerIntervals of the two Anchors,
    // to count the common oreiented reads and compute average offsets.
    info.common = 0;
    int64_t sumTwiceMarkerOffsets = 0;
    int64_t sumTwiceBaseOffsets = 0;
    auto itA = beginA;
    auto itB = beginB;
    while(itA != endA and itB != endB) {

        if(itA->orientedReadId < itB->orientedReadId) {
            ++itA;
            continue;
        }

        if(itB->orientedReadId < itA->orientedReadId) {
            ++itB;
            continue;
        }

        // We found a common OrientedReadId.
        ++info.common;
        const OrientedReadId orientedReadId = itA->orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Compute the offset in markers.
        // This assumes an ordinal offset of 1.
        // If this changes, this code will need some changes.
        const uint32_t ordinalOffsetA = ordinalOffset(anchorIdA);
        const uint32_t ordinalOffsetB = ordinalOffset(anchorIdB);
        const uint32_t ordinalA0 = itA->ordinal0;
        const uint32_t ordinalA1 = ordinalA0 + ordinalOffsetA;
        const uint32_t ordinalB0 = itB->ordinal0;
        const uint32_t ordinalB1 = ordinalB0 + ordinalOffsetB;
        sumTwiceMarkerOffsets += int64_t(ordinalB0) - int64_t(ordinalA0);
        sumTwiceMarkerOffsets += int64_t(ordinalB1) - int64_t(ordinalA1);

        // Compute the offset in bases.
        const int64_t positionA0 = int64_t(orientedReadMarkers[ordinalA0].position);
        const int64_t positionA1 = int64_t(orientedReadMarkers[ordinalA1].position);
        const int64_t positionB0 = int64_t(orientedReadMarkers[ordinalB0].position);
        const int64_t positionB1 = int64_t(orientedReadMarkers[ordinalB1].position);
        sumTwiceBaseOffsets -= positionA0;
        sumTwiceBaseOffsets -= positionA1;
        sumTwiceBaseOffsets += positionB0;
        sumTwiceBaseOffsets += positionB1;

        // Continue the joint loop.
        ++itA;
        ++itB;

    }
    info.onlyA = info.totalA - info.common;
    info.onlyB = info.totalB - info.common;

    // If there are no common reads, this is all we can do.
    if(info.common == 0) {
        info.offsetInMarkers = invalid<int64_t>;
        info.offsetInBases = invalid<int64_t>;
        info.onlyAShort = invalid<uint64_t>;
        info.onlyBShort = invalid<uint64_t>;
        return;
    }

    // Compute the estimated offsets.
    info.offsetInMarkers = int64_t(0.5 * std::round(double(sumTwiceMarkerOffsets) / double(info.common)));
    info.offsetInBases = int64_t(0.5 * std::round(double(sumTwiceBaseOffsets) / double(info.common)));



    // Now do the joint loop again, and count the onlyA and onlyB oriented reads
    // that are too short to appear in the other edge.
    itA = beginA;
    itB = beginB;
    uint64_t onlyACheck = 0;
    uint64_t onlyBCheck = 0;
    info.onlyAShort = 0;
    info.onlyBShort = 0;
    while(true) {
        if(itA == endA and itB == endB) {
            break;
        }

        else if(itB == endB or ((itA!=endA) and (itA->orientedReadId < itB->orientedReadId))) {
            // This oriented read only appears in Anchor A.
            ++onlyACheck;
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge A in this oriented read.
            const uint32_t ordinalA0 = itA->ordinal0;
            const uint32_t ordinalA1 = ordinalA0 + ordinalOffset(anchorIdA);
            const int64_t positionA0 = int64_t(orientedReadMarkers[ordinalA0].position);
            const int64_t positionA1 = int64_t(orientedReadMarkers[ordinalA1].position);

            // Find the hypothetical positions of edge B, assuming the estimated base offset.
            const int64_t positionB0 = positionA0 + info.offsetInBases;
            const int64_t positionB1 = positionA1 + info.offsetInBases;

            // If this ends up outside the read, this counts as onlyAShort.
            if(positionB0 < 0 or positionB1 >= lengthInBases) {
                ++info.onlyAShort;
            }

            ++itA;
            continue;
        }

        else if(itA == endA or ((itB!=endB) and (itB->orientedReadId < itA->orientedReadId))) {
            // This oriented read only appears in Anchor B.
            ++onlyBCheck;
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of edge B in this oriented read.
            const uint32_t ordinalB0 = itB->ordinal0;
            const uint32_t ordinalB1 = ordinalB0 + ordinalOffset(anchorIdB);
            const int64_t positionB0 = int64_t(orientedReadMarkers[ordinalB0].position);
            const int64_t positionB1 = int64_t(orientedReadMarkers[ordinalB1].position);

            // Find the hypothetical positions of edge A, assuming the estimated base offset.
            const int64_t positionA0 = positionB0 - info.offsetInBases;
            const int64_t positionA1 = positionB1 - info.offsetInBases;

            // If this ends up outside the read, this counts as onlyBShort.
            if(positionA0 < 0 or positionA1 >= lengthInBases) {
                ++info.onlyBShort;
            }

            ++itB;
            continue;
        }

        else {
            // This oriented read appears in both edges. In this loop, we
            // don't need to do anything.
            ++itA;
            ++itB;
        }
    }
    DINARA_ASSERT(onlyACheck == info.onlyA);
    DINARA_ASSERT(onlyBCheck == info.onlyB);
}



void Anchors::writeHtml(
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    AnchorPairInfo& info,
    ostream& html) const
{
    const Anchors& anchors = *this;

    // Begin the summary table.
    html <<
        "<table>"
        "<tr><th><th>On<br>anchor A<th>On<br>anchor B";

    // Total.
    html <<
        "<tr><th class=left>Total ";
    writeInformationIcon(html, "The total number of oriented reads on each of the two anchors.");
    html << "<td class=centered>" << info.totalA << "<td class=centered>" << info.totalB;

    // Common.
    html << "<tr><th class=left>Common ";
    writeInformationIcon(html, "The number of common oriented reads between the two anchors.");
    html <<
        "<td class=centered colspan = 2>" << info.common;

    // Only.
    html <<
        "<tr><th class=left>Only ";
    writeInformationIcon(html, "The number of oriented reads that appear in one anchor but not the other.");
    html <<
        "<td class=centered>" << info.onlyA << "<td class=centered>" << info.onlyB;

    // The rest of the summary table can only be written if there are common reads.
    if(info.common > 0) {

        // Only, short.
        html <<
            "<tr><th class=left>Only, short ";
        writeInformationIcon(html, "The number of oriented reads that appear in one anchor only "
            " and are too short to appear on the other anchor, based on the estimated base offset.");
        html <<
            "<td class=centered>" << info.onlyAShort << "<td class=centered>" << info.onlyBShort;

        // Only, missing.
        html <<
            "<tr><th class=left>Only, missing ";
        writeInformationIcon(html, "The number of oriented reads that appear in one anchor only "
            " and are not too short to appear on the other anchor, based on the estimated base offset.");
        html <<
            "<td class=centered>" << info.onlyA - info.onlyAShort << "<td class=centered>" << info.onlyB - info.onlyBShort;
    }

    // End the summary table.
    html << "</table>";

    // Only write out the rest if there are common reads.
    if(info.common == 0) {
        return;
    }

    // Write the table with Jaccard similarities and estimated offsets.
    using std::fixed;
    using std::setprecision;
    html <<
        "<br><table>"
        "<tr><th class=left>Jaccard similarity<td class=centered>" <<
        fixed << setprecision(2) << info.jaccard() <<
        "<tr><th class=left>Corrected Jaccard similarity<td class=centered>" <<
        fixed << setprecision(2) << info.correctedJaccard() <<
        "<tr><th class=left>Estimated offset in markers<td class=centered>" << info.offsetInMarkers <<
        "<tr><th class=left>Estimated offset in bases<td class=centered>" << info.offsetInBases <<
        "</table>";



    // Write the details table.
    html <<
        "<br>In the following table, positions in red are hypothetical, based on the above "
        "estimated base offset."
        "<p><table>";

    // Header row.
    html <<
        "<tr>"
        "<th class=centered rowspan=2>Oriented<br>read id"
        "<th class=centered colspan=2>Length"
        "<th colspan=4>Anchor A"
        "<th colspan=4>Anchor B"
        "<th rowspan=2>Ordinal offset"
        "<th rowspan=2>Base offset"
        "<th rowspan=2>Classification"
        "<tr>"
        "<th>Markers"
        "<th>Bases"
        "<th>Ordinal0"
        "<th>Ordinal1"
        "<th>Position0"
        "<th>Position1"
        "<th>Ordinal0"
        "<th>Ordinal1"
        "<th>Position0"
        "<th>Position1";

    // Prepare for the joint loop over OrientedReadIds of the two anchors.
    const auto markerIntervalsA = anchors[anchorIdA];
    const auto markerIntervalsB = anchors[anchorIdB];
    const auto beginA = markerIntervalsA.begin();
    const auto beginB = markerIntervalsB.begin();
    const auto endA = markerIntervalsA.end();
    const auto endB = markerIntervalsB.end();

    // Joint loop over the AnchorMarkerIntervals of the two Anchors.
    auto itA = beginA;
    auto itB = beginB;
    while(true) {
        if(itA == endA and itB == endB) {
            break;
        }

        else if(itB == endB or ((itA!=endA) and (itA->orientedReadId < itB->orientedReadId))) {
            // This oriented read only appears in Anchor A.
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of Anchor A in this oriented read.
            const uint32_t ordinalA0 = itA->ordinal0;
            const uint32_t ordinalA1 = ordinalA0 + ordinalOffset(anchorIdA);
            const int64_t positionA0 = int64_t(orientedReadMarkers[ordinalA0].position);
            const int64_t positionA1 = int64_t(orientedReadMarkers[ordinalA1].position);

            // Find the hypothetical positions of Anchor B, assuming the estimated base offset.
            const int64_t positionB0 = positionA0 + info.offsetInBases;
            const int64_t positionB1 = positionA1 + info.offsetInBases;
            const bool isShort = positionB0<0 or positionB1 >= lengthInBases;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << ordinalA0 <<
                "<td class=centered>" << ordinalA1 <<
                "<td class=centered>" << positionA0 <<
                "<td class=centered>" << positionA1 <<
                "<td><td>"
                "<td class=centered style='color:Red'>" << positionB0 <<
                "<td class=centered style='color:Red'>" << positionB1 << "<td><td>"
                "<td class=centered>OnlyA, " << (isShort ? "short" : "missing");

            ++itA;
            continue;
        }

        else if(itA == endA or ((itB!=endB) and (itB->orientedReadId < itA->orientedReadId))) {
            // This oriented read only appears in Anchor B.
            const OrientedReadId orientedReadId = itB->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of Anchor B in this oriented read.
            const uint32_t ordinalB0 = itB->ordinal0;
            const uint32_t ordinalB1 = ordinalB0 + ordinalOffset(anchorIdB);
            const int64_t positionB0 = int64_t(orientedReadMarkers[ordinalB0].position);
            const int64_t positionB1 = int64_t(orientedReadMarkers[ordinalB1].position);

            // Find the hypothetical positions of edge A, assuming the estimated base offset.
            const int64_t positionA0 = positionB0 - info.offsetInBases;
            const int64_t positionA1 = positionB1 - info.offsetInBases;
            const bool isShort = positionA0<0 or positionA1 >= lengthInBases;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << lengthInBases <<
                "<td><td>"
                "<td class=centered style='color:Red'>" << positionA0 <<
                "<td class=centered style='color:Red'>" << positionA1 <<
                "<td class=centered>" << ordinalB0 <<
                "<td class=centered>" << ordinalB1 <<
                "<td class=centered>" << positionB0 <<
                "<td class=centered>" << positionB1 << "<td><td>"
                "<td class=centered>OnlyB, " << (isShort ? "short" : "missing");

            ++itB;
            continue;
        }

        else {
            // This oriented read appears in both Anchors.
            const OrientedReadId orientedReadId = itA->orientedReadId;
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const int64_t lengthInBases = int64_t(reads.getReadRawSequenceLength(orientedReadId.getReadId()));

            // Get the positions of Anchor A in this oriented read.
            const uint32_t ordinalA0 = itA->ordinal0;
            const uint32_t ordinalA1 = ordinalA0 + ordinalOffset(anchorIdA);
            const int64_t positionA0 = int64_t(orientedReadMarkers[ordinalA0].position);
            const int64_t positionA1 = int64_t(orientedReadMarkers[ordinalA1].position);

            // Get the positions of Anchor B in this oriented read.
            const uint32_t ordinalB0 = itB->ordinal0;
            const uint32_t ordinalB1 = ordinalB0 + ordinalOffset(anchorIdB);
            const int64_t positionB0 = int64_t(orientedReadMarkers[ordinalB0].position);
            const int64_t positionB1 = int64_t(orientedReadMarkers[ordinalB1].position);

            // Compute estimated offsets.
            const int64_t ordinalOffset = uint64_t(ordinalB1) - uint64_t(ordinalA0);
            const int64_t baseOffset = positionB1 - positionA0;

            html <<
                "<tr><td class=centered>"
                "<a href='exploreRead?readId=" << orientedReadId.getReadId() <<
                "&strand=" << orientedReadId.getStrand() << "'>" << orientedReadId << "</a>"
                "<td class=centered>" << orientedReadMarkers.size() <<
                "<td class=centered>" << lengthInBases <<
                "<td class=centered>" << ordinalA0 <<
                "<td class=centered>" << ordinalA1 <<
                "<td class=centered>" << positionA0 <<
                "<td class=centered>" << positionA1 <<
                "<td class=centered>" << ordinalB0 <<
                "<td class=centered>" << ordinalB1 <<
                "<td class=centered>" << positionB0 <<
                "<td class=centered>" << positionB1 <<
                "<td class=centered>" << ordinalOffset <<
                "<td class=centered>" << baseOffset <<
                "<td class=centered>Common";

            ++itA;
            ++itB;
        }
    }

    // Finish the details table.
    html << "</table>";

}



// Return true if the second Anchor is adjacent to the first one,
// as seen by at least one of the common oriented reads.
// THIS CRITERION MAY BE TOO LOSE.
bool Anchors::areAdjacentAnchors(AnchorId anchorId0, AnchorId anchorId1) const
{
    const auto markerIntervals0 = anchorMarkerIntervals[anchorId0];
    const auto markerIntervals1 = anchorMarkerIntervals[anchorId1];

    // Joint loop over the marker intervals.
    auto it0 = markerIntervals0.begin();
    auto it1 = markerIntervals1.begin();
    const auto end0 = markerIntervals0.end();
    const auto end1 = markerIntervals1.end();
    while(it0 != end0 and it1 != end1) {
        const OrientedReadId orientedReadId0 = it0->orientedReadId;
        const OrientedReadId orientedReadId1 = it1->orientedReadId;

        if(orientedReadId0 < orientedReadId1) {
            ++it0;
        } else if(orientedReadId1 < orientedReadId0) {
            ++it1;
        } else {
            if(it0->ordinal0 + ordinalOffset(anchorId0) == it1->ordinal0) {
                return true;
            }
            ++it0;
            ++it1;
        }
    }

    return false;

}


// Anchors are numbered such that each pair of reverse complemented
// AnchorIds are numbered (n, n+1), where n is even, n = 2*m.
// We represent an AnchorId as a string as follows:
// - AnchorId n is represented as m+
// - AnchorId n+1 is represented as m-
// For example, the reverse complemented pair (150, 151) is represented as (75+, 75-).

string dinara::mode3::anchorIdToString(AnchorId n)
{
    std::ostringstream s;

    const AnchorId m = (n >> 1);
    s << m;

    if(n & 1) {
        s << "-";
    } else {
        s << "+";
    }

    return s.str();
}



AnchorId dinara::mode3::anchorIdFromString(const string& s)
{

    if(s.size() < 2) {
        return invalid<AnchorId>;
    }

    const char cLast = s.back();

    uint64_t lastBit;
    if(cLast == '+') {
        lastBit = 0;
    } else if(cLast == '-') {
        lastBit = 1;
    } else {
        return invalid<AnchorId>;
    }

    const uint64_t m = std::stoul(s.substr(0, s.size()-1));

    return 2* m + lastBit;
}



void Anchors::computeJourneys(uint64_t threadCount, const vector<uint64_t>* keepAnchorWords)
{
    performanceLog << timestamp << "Anchors::computeJourneys begins." << endl;

    const uint64_t orientedReadCount = 2 * reads.readCount();

    // If we are recomputing journeys (potentially with filtering), make sure we are not
    // appending to stale state and that all marker intervals start with invalid positions.
    if (journeys.isOpen()) {
        journeys.remove();
    }

    journeyAnchorFilterData.keepAnchorWords = keepAnchorWords;

    // Reset positionInJourney on all marker intervals. If we filter anchors, some anchors
    // may no longer appear in any journey; their positionInJourney must not retain stale values.
    {
        const uint64_t anchorBatchCount = 1000;
        setupLoadBalancing(size(), anchorBatchCount);
        runThreads(&Anchors::resetPositionInJourneyThreadFunction, threadCount);
    }

    // Pass1: make space for the journeysWithOrdinals.
    journeysWithOrdinals.createNew(largeDataName("tmp-JourneysWithOrdinals"), largeDataPageSize);
    journeysWithOrdinals.beginPass1(orientedReadCount);
    const uint64_t anchorBatchCount = 1000;
    setupLoadBalancing(size(), anchorBatchCount);
    runThreads(&Anchors::computeJourneysThreadFunction1, threadCount);

    // Pass2: store the unsorted journeysWithOrdinals.
    journeysWithOrdinals.beginPass2();
    setupLoadBalancing(size(), anchorBatchCount);
    runThreads(&Anchors::computeJourneysThreadFunction2, threadCount);
    journeysWithOrdinals.endPass2();

    // Pass 3:sort the journeysWithOrdinals and make space for the journeys
    journeys.createNew(largeDataName("Journeys"), largeDataPageSize);
    journeys.beginPass1(orientedReadCount);
    const uint64_t orientedReadBatchCount = 1000;
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Anchors::computeJourneysThreadFunction3, threadCount);

    // Pass 4: copy the sorted journeysWithOrdinals to the journeys.
    journeys.beginPass2();
    setupLoadBalancing(orientedReadCount, orientedReadBatchCount);
    runThreads(&Anchors::computeJourneysThreadFunction4, threadCount);
    journeys.endPass2(false, true);

    journeysWithOrdinals.remove();

    performanceLog << timestamp << "Anchors::computeJourneys ends." << endl;

    writeJourneys();

    journeyAnchorFilterData.keepAnchorWords = nullptr;
}



void Anchors::computeJourneysThreadFunction1(uint64_t /* threadId */)
{
    computeJourneysThreadFunction12(1);
}


void Anchors::writeJourneys() const
{
    const ReadId readCount = reads.readCount();

    ofstream csv("Journeys.csv");
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            csv << orientedReadId << ",";
            const auto journey = journeys[orientedReadId.getValue()];
            for(const AnchorId anchorId: journey) {
                csv << anchorIdToString(anchorId) << ",";
            }
            csv << "\n";
        }
    }
}

namespace {
    struct OverlapEvent {
        uint32_t pos = 0;
        uint8_t type = 0; // 0=end, 1=start
        uint32_t otherOrientedReadValue = 0;
    };

    inline bool getAlignmentIntervalOnOrientedRead(
        const Reads& reads,
        const AlignmentData& ad,
        OrientedReadId orientedReadId,
        uint32_t& start,
        uint32_t& end)
    {
        const ReadId readId = orientedReadId.getReadId();
        if (readId == ad.readIds[0]) {
            start = ad.qs;
            end = ad.qe;
        } else if (readId == ad.readIds[1]) {
            start = ad.ts;
            end = ad.te;
        } else {
            return false;
        }

        const uint32_t readLen = uint32_t(reads.getRead(readId).baseCount);
        if (start > end) return false;
        if (end > readLen) end = readLen;
        if (start > readLen) start = readLen;
        if (end <= start) return false;

        if (orientedReadId.getStrand() == 1) {
            // ad.* are in forward coordinates of the read; invert for reverse-complement strand.
            const uint32_t rcStart = readLen - end;
            const uint32_t rcEnd = readLen - start;
            start = rcStart;
            end = rcEnd;
            if (end <= start) return false;
        }

        return true;
    }
}

void Anchors::writeStableOverlapIntervalsBestAnchors(
    const ReadGraph& readGraph,
    const MemoryMapped::Vector<AlignmentData>& alignmentData,
    uint64_t threadCount,
    const string& outCsvName,
    bool strand0Only,
    uint32_t minIntervalLength)
{
    performanceLog << timestamp << "Anchors::writeStableOverlapIntervalsBestAnchors begins." << endl;

    reads.checkReadsAreOpen();
    DINARA_ASSERT(readGraph.edges.isOpen);
    DINARA_ASSERT(readGraph.connectivity.isOpen());

    stableOverlapIntervalsBestAnchorsData.readGraph = &readGraph;
    stableOverlapIntervalsBestAnchorsData.alignmentData = &alignmentData;
    stableOverlapIntervalsBestAnchorsData.outCsvName = outCsvName;
    stableOverlapIntervalsBestAnchorsData.strand0Only = strand0Only;
    stableOverlapIntervalsBestAnchorsData.minIntervalLength = minIntervalLength;
    stableOverlapIntervalsBestAnchorsData.threadFileNames.assign(threadCount, string());

    const uint64_t workItems = strand0Only ? reads.readCount() : (2 * reads.readCount());
    const uint64_t batchSize = 2000;
    setupLoadBalancing(workItems, batchSize);
    runThreads(&Anchors::writeStableOverlapIntervalsBestAnchorsThreadFunction, threadCount);

    // Merge per-thread files.
    {
        ofstream out(outCsvName);
        out << "OrientedReadId,IntervalBegin,IntervalEnd,IntervalLength,"
               "ActiveOverlapCount,CandidateAnchorCount,ChosenAnchorId,DesiredReads,"
               "CoveredDesiredReads,MissingDesiredReads,ExtraReads,ChosenAnchorCoverage,"
               "ChosenAnchorCount,ChosenAnchorIds,FinalCoveredDesiredReads,FinalMissingDesiredReads,TotalExtraReads\n";

        for (uint64_t t = 0; t < threadCount; ++t) {
            const string& fn = stableOverlapIntervalsBestAnchorsData.threadFileNames[t];
            if (fn.empty()) continue;
            ifstream in(fn);
            if (!in) continue;
            string line;
            bool first = true;
            while (std::getline(in, line)) {
                if (first) { first = false; continue; } // skip header
                out << line << "\n";
            }
            in.close();
            std::error_code ec;
            std::filesystem::remove(fn, ec);
        }
    }

    // Clean up pointers.
    stableOverlapIntervalsBestAnchorsData.readGraph = nullptr;
    stableOverlapIntervalsBestAnchorsData.alignmentData = nullptr;
    stableOverlapIntervalsBestAnchorsData.keepAnchorWords = nullptr;
    stableOverlapIntervalsBestAnchorsData.threadFileNames.clear();

    performanceLog << timestamp << "Anchors::writeStableOverlapIntervalsBestAnchors ends." << endl;
}



vector<uint64_t> Anchors::writeStableOverlapIntervalsBestAnchorsAndCollectKeptAnchors(
    const ReadGraph& readGraph,
    const MemoryMapped::Vector<AlignmentData>& alignmentData,
    uint64_t threadCount,
    const string& outCsvName,
    bool strand0Only,
    uint32_t minIntervalLength,
    bool keepReverseComplementPairs)
{
    // Build a global bitset of anchors chosen as "best" in at least one stable interval.
    // This can be used to filter and recompute journeys afterwards.
    const uint64_t wordCount = (size() + 63ULL) / 64ULL;
    vector<std::atomic<uint64_t>> keepWordsAtomic(wordCount);
    for (uint64_t i = 0; i < wordCount; ++i) {
        keepWordsAtomic[i].store(0ULL, std::memory_order_relaxed);
    }

    stableOverlapIntervalsBestAnchorsData.keepAnchorWords = &keepWordsAtomic;
    writeStableOverlapIntervalsBestAnchors(
        readGraph, alignmentData, threadCount, outCsvName, strand0Only, minIntervalLength);
    stableOverlapIntervalsBestAnchorsData.keepAnchorWords = nullptr;

    vector<uint64_t> keepWords(wordCount);
    for (uint64_t i = 0; i < wordCount; ++i) {
        keepWords[i] = keepWordsAtomic[i].load(std::memory_order_relaxed);
    }

    if (keepReverseComplementPairs) {
        // Reverse-complement AnchorId is `anchorId^1` (toggles the last bit).
        // In the bitset, this just swaps each adjacent bit-pair (0<->1, 2<->3, ...),
        // so we can add the RC anchors in a single post-pass without any atomics.
        static constexpr uint64_t evenMask = 0x5555555555555555ULL; // bits 0,2,4,...
        static constexpr uint64_t oddMask  = 0xAAAAAAAAAAAAAAAAULL; // bits 1,3,5,...
        for (uint64_t i = 0; i < wordCount; ++i) {
            const uint64_t w = keepWords[i];
            keepWords[i] = w | ((w & evenMask) << 1) | ((w & oddMask) >> 1);
        }
    }
    return keepWords;
}



void Anchors::writeStableOverlapIntervalsBestAnchorsThreadFunction(uint64_t threadId)
{
    const auto* readGraphPtr = stableOverlapIntervalsBestAnchorsData.readGraph;
    const auto* alignmentDataPtr = stableOverlapIntervalsBestAnchorsData.alignmentData;
    DINARA_ASSERT(readGraphPtr);
    DINARA_ASSERT(alignmentDataPtr);
    const ReadGraph& readGraph = *readGraphPtr;
    const MemoryMapped::Vector<AlignmentData>& alignmentData = *alignmentDataPtr;

    const string threadFileName =
        stableOverlapIntervalsBestAnchorsData.outCsvName + ".thread" + std::to_string(threadId) + ".csv";
    stableOverlapIntervalsBestAnchorsData.threadFileNames[threadId] = threadFileName;
    ofstream out(threadFileName);
    out << "OrientedReadId,IntervalBegin,IntervalEnd,IntervalLength,"
           "ActiveOverlapCount,CandidateAnchorCount,ChosenAnchorId,DesiredReads,"
           "CoveredDesiredReads,MissingDesiredReads,ExtraReads,ChosenAnchorCoverage,"
           "ChosenAnchorCount,ChosenAnchorIds,FinalCoveredDesiredReads,FinalMissingDesiredReads,TotalExtraReads\n";

    const bool strand0Only = stableOverlapIntervalsBestAnchorsData.strand0Only;
    const uint32_t minIntervalLength = stableOverlapIntervalsBestAnchorsData.minIntervalLength;

    vector<OverlapEvent> events;
    events.reserve(256);
    std::unordered_set<uint32_t> active;
    active.reserve(64);
    std::unordered_set<uint32_t> uncovered;
    uncovered.reserve(64);

    // Loop over all batches assigned to this thread.
    uint64_t begin = 0, end = 0;
    while (getNextBatch(begin, end)) {
        for (uint64_t item = begin; item != end; ++item) {
            OrientedReadId orientedReadId;
            if (strand0Only) {
                orientedReadId = OrientedReadId(ReadId(item), 0);
            } else {
                orientedReadId = OrientedReadId::fromValue(uint32_t(item));
            }

            const auto journey = journeys[orientedReadId.getValue()];
            if (journey.empty()) continue;

            // Build overlap events from incident readGraph edges.
            events.clear();
            const auto inc = readGraph.connectivity[orientedReadId.getValue()];
            events.reserve(2 * inc.size());
            for (const uint32_t edgeId : inc) {
                if (edgeId >= readGraph.edges.size()) continue;
                const ReadGraphEdge& edge = readGraph.edges[edgeId];
                const uint64_t alignmentId = edge.alignmentId;
                if (alignmentId >= alignmentData.size()) continue;
                const AlignmentData& ad = alignmentData[alignmentId];
                if (!ad.keptByBothSides()) continue;

                uint32_t s = 0, e = 0;
                if (!getAlignmentIntervalOnOrientedRead(reads, ad, orientedReadId, s, e)) continue;

                const OrientedReadId other = edge.getOther(orientedReadId);
                const uint32_t otherValue = other.getValue();
                events.push_back({s, 1, otherValue});
                events.push_back({e, 0, otherValue});
            }

            if (events.empty()) continue;

            std::sort(events.begin(), events.end(), [](const OverlapEvent& a, const OverlapEvent& b) {
                if (a.pos != b.pos) return a.pos < b.pos;
                return a.type < b.type; // end(0) before start(1)
            });

            // Precompute anchor mid positions for this journey (monotone increasing).
            const auto orientedMarkers = markers[orientedReadId.getValue()];
            vector<uint32_t> journeyAnchorMid(journey.size());
            for (uint64_t i = 0; i < journey.size(); ++i) {
                const AnchorId anchorId = journey[i];
                const uint32_t ordinal0 = getFirstOrdinal(anchorId, orientedReadId);
                const uint32_t p = orientedMarkers[ordinal0].position;
                journeyAnchorMid[i] = p + uint32_t(kHalf);
            }

                // Sweep stable intervals and pick best anchor.
                active.clear();
                uint64_t journeyIdx = 0;

                // Sweep line over overlap start/end events.
                // After applying all events at position `pos`, the set `active` contains exactly
                // the overlaps that are active over the half-open interval [pos, nextPos).
                uint64_t i = 0;
                while (i < events.size()) {
                    const uint32_t pos = events[i].pos;

                // Apply all end events at pos.
                while (i < events.size() && events[i].pos == pos && events[i].type == 0) {
                    active.erase(events[i].otherOrientedReadValue);
                    ++i;
                }
                // Apply all start events at pos.
                while (i < events.size() && events[i].pos == pos && events[i].type == 1) {
                    active.insert(events[i].otherOrientedReadValue);
                    ++i;
                }

                const uint32_t nextPos = (i < events.size()) ? events[i].pos : pos;
                if (nextPos <= pos) {
                    continue;
                }

                const uint32_t intervalLen = nextPos - pos;
                if (minIntervalLength && intervalLen < minIntervalLength) {
                    continue;
                }
                if (active.empty()) {
                    continue;
                }

                // Advance journeyIdx to first anchor inside the interval.
                while (journeyIdx < journeyAnchorMid.size() && journeyAnchorMid[journeyIdx] < pos) {
                    ++journeyIdx;
                }
                uint64_t jEnd = journeyIdx;
                while (jEnd < journeyAnchorMid.size() && journeyAnchorMid[jEnd] < nextPos) {
                    ++jEnd;
                }
                const uint64_t candidateCount = jEnd - journeyIdx;
                if (candidateCount == 0) {
                    continue;
                }

                const uint32_t desiredReads = uint32_t(active.size()) + 1U; // active overlaps + self
                const uint32_t selfValue = orientedReadId.getValue();

                // Greedy multi-anchor selection within the interval:
                // Many regions will have a single anchor that cleanly contains all desired reads (missing=0).
                // Collapsed/noisy regions often do not: overlap partners can be split across multiple anchors.
                // In that case, keeping just the "best single anchor" can still retain a collapsed anchor.
                //
                // Failsafe: iteratively pick anchors that cover *new* uncovered desired reads until we cover
                // the whole desired set or we can no longer make progress. We still penalize anchors that
                // contain unrelated reads ("extra") to avoid keeping broadly-collapsed anchors.

                // Initialize the desired set and the uncovered desired reads.
                uncovered.clear();
                uncovered.insert(selfValue);
                for (const uint32_t v : active) {
                    uncovered.insert(v);
                }

                vector<AnchorId> chosenAnchors;
                chosenAnchors.reserve(4);
                uint32_t totalExtraReads = 0;

                // We cap the number of chosen anchors per interval to prevent pathological explosion
                // in highly collapsed regions.
                const uint32_t maxChosenAnchors =
                    std::min<uint32_t>(8U, desiredReads); // never need more than desiredReads anchors.

                // Track a "primary" anchor score compatible with the old output columns.
                AnchorId primaryAnchorId = invalid<AnchorId>;
                uint32_t primaryCovered = 0;
                uint32_t primaryMissing = std::numeric_limits<uint32_t>::max();
                uint32_t primaryExtra = std::numeric_limits<uint32_t>::max();
                uint32_t primaryCoverage = 0;

                // Greedy loop.
                for (uint32_t step = 0; step < maxChosenAnchors && !uncovered.empty(); ++step) {
                    AnchorId bestAnchorId = invalid<AnchorId>;
                    uint32_t bestGain = 0;
                    uint32_t bestExtra = std::numeric_limits<uint32_t>::max();
                    uint32_t bestCoveredDesired = 0;
                    uint32_t bestAnchorCoverage = 0;

                    for (uint64_t j = journeyIdx; j < jEnd; ++j) {
                        const AnchorId anchorId = journey[j];
                        bool alreadyChosen = false;
                        for (const AnchorId x : chosenAnchors) {
                            if (x == anchorId) {
                                alreadyChosen = true;
                                break;
                            }
                        }
                        if (alreadyChosen) continue;

                        const Anchor anchor = (*this)[anchorId];

                        // Compute:
                        // - coveredDesired: |anchor ∩ desired|
                        // - gain: |anchor ∩ uncovered|
                        // - extra: |anchor \ desired|
                        uint32_t coveredDesired = 0;
                        uint32_t gain = 0;
                        bool hasSelf = false;
                        for (const auto& mi : anchor) {
                            const uint32_t v = mi.orientedReadId.getValue();
                            const bool inDesired = (v == selfValue) || (active.find(v) != active.end());
                            if (inDesired) {
                                ++coveredDesired;
                                if (uncovered.find(v) != uncovered.end()) {
                                    ++gain;
                                }
                                if (v == selfValue) {
                                    hasSelf = true;
                                }
                            }
                        }
                        if (!hasSelf) continue;
                        if (gain == 0) continue;

                        const uint32_t extra = uint32_t(anchor.size()) - coveredDesired;

                        // Choose anchor with maximum gain, then minimum extra, then maximum coveredDesired.
                        if (gain > bestGain ||
                            (gain == bestGain && extra < bestExtra) ||
                            (gain == bestGain && extra == bestExtra && coveredDesired > bestCoveredDesired)) {
                            bestAnchorId = anchorId;
                            bestGain = gain;
                            bestExtra = extra;
                            bestCoveredDesired = coveredDesired;
                            bestAnchorCoverage = uint32_t(anchor.size());
                        }
                    }

                    if (bestAnchorId == invalid<AnchorId>) {
                        // No candidate can cover any additional desired reads.
                        break;
                    }

                    // Apply selection.
                    chosenAnchors.push_back(bestAnchorId);
                    totalExtraReads += bestExtra;

                    // Update uncovered set by removing desired reads present on the chosen anchor.
                    {
                        const Anchor chosenAnchor = (*this)[bestAnchorId];
                        for (const auto& mi : chosenAnchor) {
                            const uint32_t v = mi.orientedReadId.getValue();
                            uncovered.erase(v);
                        }
                    }

                    // Capture primary-anchor columns (the first chosen anchor).
                    if (step == 0) {
                        primaryAnchorId = bestAnchorId;
                        primaryCovered = bestCoveredDesired;
                        primaryMissing = desiredReads - bestCoveredDesired;
                        primaryExtra = bestExtra;
                        primaryCoverage = bestAnchorCoverage;
                    }

                    // Record chosen anchor in the global keep-bitset (if requested).
                    if (stableOverlapIntervalsBestAnchorsData.keepAnchorWords) {
                        auto& keepWords = *stableOverlapIntervalsBestAnchorsData.keepAnchorWords;
                        const uint64_t word = uint64_t(bestAnchorId) >> 6;
                        const uint64_t bit = uint64_t(bestAnchorId) & 63ULL;
                        if (word < keepWords.size()) {
                            keepWords[word].fetch_or(1ULL << bit, std::memory_order_relaxed);
                        }
                    }
                }

                if (chosenAnchors.empty()) continue;

                const uint32_t finalMissing = uint32_t(uncovered.size());
                const uint32_t finalCovered = desiredReads - finalMissing;

                // Serialize chosen anchors as a semicolon-separated list.
                std::ostringstream chosenList;
                for (uint64_t k = 0; k < chosenAnchors.size(); ++k) {
                    if (k) chosenList << ";";
                    chosenList << anchorIdToString(chosenAnchors[k]);
                }

                out << orientedReadId << "," << pos << "," << nextPos << "," << intervalLen << ","
                    << active.size() << "," << candidateCount << "," << anchorIdToString(primaryAnchorId) << ","
                    << desiredReads << "," << primaryCovered << "," << primaryMissing << "," << primaryExtra << ","
                    << primaryCoverage << ","
                    << chosenAnchors.size() << "," << chosenList.str() << ","
                    << finalCovered << "," << finalMissing << "," << totalExtraReads << "\n";
            }
        }
    }
}



void Anchors::computeJourneysThreadFunction2(uint64_t /* threadId */)
{
    computeJourneysThreadFunction12(2);
}



void Anchors::computeJourneysThreadFunction12(uint64_t pass)
{
    Anchors& anchors = *this;
    const auto* keepWords = anchors.journeyAnchorFilterData.keepAnchorWords;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all AnchorIds in this batch.
        for(AnchorId anchorId=begin; anchorId!=end; anchorId++) {
            if (keepWords) {
                const uint64_t word = uint64_t(anchorId) >> 6;
                const uint64_t bit = uint64_t(anchorId) & 63ULL;
                if (word >= keepWords->size() || (((*keepWords)[word] >> bit) & 1ULL) == 0ULL) {
                    continue;
                }
            }
            Anchor anchor = anchors[anchorId];

            // Loop over the marker intervals of this Anchor.
            for(const auto& anchorMarkerInterval: anchor) {
                const auto orientedReadIdValue = anchorMarkerInterval.orientedReadId.getValue();

                if(pass == 1) {
                    journeysWithOrdinals.incrementCountMultithreaded(orientedReadIdValue);
                } else {
                    journeysWithOrdinals.storeMultithreaded(
                        orientedReadIdValue, {anchorId, anchorMarkerInterval.ordinal0});
                }

            }
        }

    }
}



void Anchors::resetPositionInJourneyThreadFunction(uint64_t /* threadId */)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin = 0, end = 0;
    while (getNextBatch(begin, end)) {
        for (AnchorId anchorId = begin; anchorId != end; ++anchorId) {
            span<AnchorMarkerInterval> markerIntervals = anchorMarkerIntervals[anchorId];
            for (AnchorMarkerInterval& mi : markerIntervals) {
                mi.positionInJourney = invalid<uint32_t>;
            }
        }
    }
}



void Anchors::computeJourneysThreadFunction3(uint64_t /* threadId */)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            auto v = journeysWithOrdinals[orientedReadValue];
            sort(v.begin(), v.end(), OrderPairsBySecondOnly<uint64_t, uint32_t>());
            journeys.incrementCountMultithreaded(orientedReadValue, v.size());
        }
    }
}



void Anchors::computeJourneysThreadFunction4(uint64_t /* threadId */)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all oriented reads assigned to this thread.
        for(uint64_t orientedReadValue=begin; orientedReadValue!=end; orientedReadValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(orientedReadValue));

            // Copy the journeysWithOrdinals to the journeys.
            const auto v = journeysWithOrdinals[orientedReadValue];
            const auto journey = journeys[orientedReadValue];
            DINARA_ASSERT(journey.size() == v.size());
            for(uint64_t i=0; i<v.size(); i++) {
                journey[i] = v[i].first;
            }

            // Store journey information for this oriented read in the marker interval.
            for(uint64_t position=0; position<journey.size(); position++) {
                const AnchorId anchorId = journey[position];
                span<AnchorMarkerInterval> markerIntervals = anchorMarkerIntervals[anchorId];
                bool found = false;
                for(AnchorMarkerInterval& markerInterval: markerIntervals) {
                    if(markerInterval.orientedReadId == orientedReadId) {
                        markerInterval.positionInJourney = uint32_t(position);
                        found = true;
                        break;
                    }
                }
                DINARA_ASSERT(found);
            }
        }
    }
}



// For a given AnchorId, follow the read journeys forward by one step.
// Return a vector of the AnchorIds reached in this way.
// The count vector is the number of oriented reads each of the AnchorIds.
void Anchors::findChildren(
    AnchorId anchorId,
    vector<AnchorId>& children,
    vector<uint64_t>& count,
    uint64_t minCoverage) const
{
    children.clear();
    for(const auto& markerInterval: anchorMarkerIntervals[anchorId]) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto journey = journeys[orientedReadId.getValue()];
        const uint64_t position = markerInterval.positionInJourney;
        if (position == invalid<uint32_t> || position >= journey.size()) {
            continue;
        }
        const uint64_t nextPosition = position + 1;
        if(nextPosition < journey.size()) {
            const AnchorId nextAnchorId = journey[nextPosition];
            children.push_back(nextAnchorId);
        }
    }

    deduplicateAndCountWithThreshold(children, count, minCoverage);
}



// For a given AnchorId, follow the read journeys backward by one step.
// Return a vector of the AnchorIds reached in this way.
// The count vector is the number of oriented reads each of the AnchorIds.
void Anchors::findParents(
    AnchorId anchorId,
    vector<AnchorId>& parents,
    vector<uint64_t>& count,
    uint64_t minCoverage) const
{
    parents.clear();
    for(const auto& markerInterval: anchorMarkerIntervals[anchorId]) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto journey = journeys[orientedReadId.getValue()];
        const uint64_t position = markerInterval.positionInJourney;
        if (position == invalid<uint32_t> || position >= journey.size()) {
            continue;
        }
        if(position > 0) {
            const uint64_t previousPosition = position - 1;
            const AnchorId previousAnchorId = journey[previousPosition];
            parents.push_back(previousAnchorId);
        }
    }

    deduplicateAndCountWithThreshold(parents, count, minCoverage);
}



// Get the first ordinal for the AnchorMarkerInterval corresponding to a
// given AnchorId and OrientedReadId.
// This asserts if the given AnchorId does not contain an AnchorMarkerInterval
// for the requested OrientedReadId.
uint32_t Anchors::getFirstOrdinal(AnchorId anchorId, OrientedReadId orientedReadId) const
{
    for(const auto& markerInterval: anchorMarkerIntervals[anchorId]) {
        if(markerInterval.orientedReadId == orientedReadId) {
            return markerInterval.ordinal0;
        }
    }

    DINARA_ASSERT(0);
}



void Anchors::writeCoverageHistogram() const
{
    vector<uint64_t> histogram;
    for(AnchorId anchorId=0; anchorId<anchorMarkerIntervals.size(); anchorId++) {
        const uint64_t coverage = anchorMarkerIntervals.size(anchorId);
        if(coverage >= histogram.size()) {
            histogram.resize(coverage + 1, 0);
        }
        ++histogram[coverage];
    }

    ofstream csv("AnchorCoverageHistogram.csv");
    csv << "Coverage,Frequency\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        csv << coverage << "," << histogram[coverage] << "\n";
    }
}



// Read following.
void Anchors::followOrientedReads(
    AnchorId anchorId0,
    uint64_t direction,                         // 0 = forward, 1 = backward
    uint64_t minCommonCount,
    double minJaccard,
    double minCorrectedJaccard,
    vector< pair<AnchorId, AnchorPairInfo> >& anchorInfos
    ) const
{
    const Anchor anchor0 = (*this)[anchorId0];

    // Gather all AnchorIds reached by the forward or backward portions of the
    // journeys of the oriented reads on this anchor.
    vector<AnchorId> anchorIds;
    for(const AnchorMarkerInterval& anchorMarkerInterval: anchor0) {
        const OrientedReadId orientedReadId = anchorMarkerInterval.orientedReadId;
        const uint64_t position0 = anchorMarkerInterval.positionInJourney;
        if (position0 == invalid<uint32_t>) {
            continue;
        }
        const auto journey = journeys[orientedReadId.getValue()];
        if (position0 >= journey.size()) {
            continue;
        }

        // Figure out the forward or backward portion of the journey.
        uint64_t begin;
        uint64_t end;
        if(direction == 0) {
            begin = position0 + 1;
            end = journey.size();
        } else {
            begin = 0;
            end = position0;
        }

        // Copy the AnchorIds on this portion of the journey.
        copy(journey.begin() + begin, journey.begin() + end, back_inserter(anchorIds));
    }

    // Only keep the ones we saw at least minCommonCount times.
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(anchorIds, count, minCommonCount);



    // Gather the ones that satisfy our criteria.
    anchorInfos.clear();
    for(const AnchorId anchorId1: anchorIds) {
        AnchorPairInfo info;
        if(direction == 0) {
            analyzeAnchorPair(anchorId0, anchorId1, info);
        } else {
            analyzeAnchorPair(anchorId1, anchorId0, info);
        }
        if(info.common < minCommonCount) {
            continue;
        }
        if(info.common == 0) {
            continue;
        }
        if(info.jaccard() < minJaccard) {
            continue;
        }
        if(info.correctedJaccard() < minCorrectedJaccard) {
            continue;
        }
        anchorInfos.push_back(make_pair(anchorId1, info));
    }


    // Sort them by offset.
    class SortHelper {
    public:
        bool operator()(
            const pair<AnchorId, AnchorPairInfo>& x,
            const pair<AnchorId, AnchorPairInfo>& y
            ) const
        {
            return x.second.offsetInBases < y.second.offsetInBases;
        }
    };
    sort(anchorInfos.begin(), anchorInfos.end(), SortHelper());
}



// For each read, write out the largest gap between adjacent anchors.
// The two oriented reads for a read have the same gaps.
void Anchors::writeAnchorGapsByRead() const
{

    // Open the output csv file and write a header.
    ofstream csv("AnchorGaps.csv");
    csv << "ReadId,Length,Anchor count,Gap\n";

    // Loop over all reads.
    for(ReadId readId=0; readId<reads.readCount(); readId++) {
        const uint64_t readLength = reads.getRead(readId).baseCount;

        // Put in on strand 0. It would have the same gap on strand 1.
        const OrientedReadId orientedReadId(readId, 0);

        // Get the markers and the journey.
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];
        const auto journey = journeys[orientedReadId.getValue()];

        // Loop over adjacent positions in journey.
        uint64_t previousPosition = 0;
        uint64_t maxGap = 0;
        for(uint64_t i=0; i<=journey.size(); i++) {

            uint64_t position = invalid<uint64_t>;
            if(i == journey.size()) {
                position = readLength;
            } else {
                const AnchorId anchorId = journey[i];
                const Anchor anchor = (*this)[anchorId];
                for(const AnchorMarkerInterval& markerInterval: anchor) {
                    if(markerInterval.orientedReadId == orientedReadId) {
                        const uint32_t ordinal = markerInterval.ordinal0;
                        const CompressedMarker& marker = orientedReadMarkers[ordinal];
                        position = marker.position;
                        break;
                    }
                }
            }
            DINARA_ASSERT(position != invalid<uint64_t>);
            DINARA_ASSERT(position >= previousPosition);

            const uint64_t gap = position - previousPosition;
            maxGap = max(maxGap, gap);

            previousPosition = position;
        }
        csv << readId << "," << readLength << "," << journey.size() << "," << maxGap << "\n";
    }
}



void Anchors::writeJourneyEndpoints() const
{
    class PairStats {
    public:
        uint64_t readCount = 0;
        uint64_t anchorCountSum = 0;
        uint64_t terminalSpanSum = 0;
        uint64_t unanchoredStartSum = 0;
        uint64_t unanchoredEndSum = 0;
    };

    // Use strand 0 only so each physical read contributes at most once.
    std::map<std::pair<AnchorId, AnchorId>, PairStats> pairStats;
    uint64_t readCountWithJourneys = 0;
    uint64_t readCountWithTerminalPair = 0;

    ofstream perReadCsv("JourneyEndpoints.csv");
    perReadCsv << "ReadId,AnchorCount,StartAnchor,EndAnchor,CanonicalStartAnchor,CanonicalEndAnchor,"
                  "StartPosition,EndPosition,ReadLength,TerminalSpan,UnanchoredStart,UnanchoredEnd\n";

    for(ReadId readId=0; readId<reads.readCount(); readId++) {
        const OrientedReadId orientedReadId(readId, 0);
        const auto journey = journeys[orientedReadId.getValue()];
        const uint64_t anchorCount = journey.size();
        if(anchorCount == 0) {
            perReadCsv << readId << ",0,,,,,,,,,,\n";
            continue;
        }

        ++readCountWithJourneys;
        const uint64_t readLength = reads.getRead(readId).baseCount;
        const AnchorId startAnchorId = journey.front();
        const AnchorId endAnchorId = journey.back();
        const auto canonicalPair = canonicalEndpointPair(startAnchorId, endAnchorId);

        const auto orientedReadMarkers = markers[orientedReadId.getValue()];
        const uint32_t startOrdinal = getFirstOrdinal(startAnchorId, orientedReadId);
        const uint32_t endOrdinal = getFirstOrdinal(endAnchorId, orientedReadId);
        const uint64_t startPosition = orientedReadMarkers[startOrdinal].position;
        const uint64_t endPosition = orientedReadMarkers[endOrdinal].position;
        DINARA_ASSERT(endPosition >= startPosition);
        DINARA_ASSERT(readLength >= endPosition);

        const uint64_t terminalSpan = endPosition - startPosition;
        const uint64_t unanchoredStart = startPosition;
        const uint64_t unanchoredEnd = readLength - endPosition;

        perReadCsv << readId << "," << anchorCount << ","
            << anchorIdToString(startAnchorId) << ","
            << anchorIdToString(endAnchorId) << ","
            << anchorIdToString(canonicalPair.first) << ","
            << anchorIdToString(canonicalPair.second) << ","
            << startPosition << "," << endPosition << "," << readLength << ","
            << terminalSpan << "," << unanchoredStart << "," << unanchoredEnd << "\n";

        if(anchorCount < 2) {
            continue;
        }

        ++readCountWithTerminalPair;
        PairStats& stats = pairStats[canonicalPair];
        ++stats.readCount;
        stats.anchorCountSum += anchorCount;
        stats.terminalSpanSum += terminalSpan;
        stats.unanchoredStartSum += unanchoredStart;
        stats.unanchoredEndSum += unanchoredEnd;
    }

    vector< pair< pair<AnchorId, AnchorId>, PairStats> > sortedPairStats(
        pairStats.begin(), pairStats.end());
    sort(sortedPairStats.begin(), sortedPairStats.end(),
        [](const auto& x, const auto& y) {
            if(x.second.readCount != y.second.readCount) {
                return x.second.readCount > y.second.readCount;
            }
            return x.first < y.first;
        });

    ofstream pairCsv("JourneyEndpointPairs.csv");
    pairCsv << "CanonicalStartAnchor,CanonicalEndAnchor,ReadCount,"
               "AverageAnchorCount,AverageTerminalSpan,AverageUnanchoredStart,AverageUnanchoredEnd\n";
    for(const auto& p: sortedPairStats) {
        const std::pair<AnchorId, AnchorId>& pair = p.first;
        const PairStats& stats = p.second;
        const double readCountInverse = 1.0 / double(stats.readCount);
        pairCsv
            << anchorIdToString(pair.first) << ","
            << anchorIdToString(pair.second) << ","
            << stats.readCount << ","
            << double(stats.anchorCountSum) * readCountInverse << ","
            << double(stats.terminalSpanSum) * readCountInverse << ","
            << double(stats.unanchoredStartSum) * readCountInverse << ","
            << double(stats.unanchoredEndSum) * readCountInverse << "\n";
    }

    performanceLog << timestamp
        << "Wrote JourneyEndpoints.csv for " << readCountWithJourneys << " reads and "
        << "JourneyEndpointPairs.csv for " << pairStats.size()
        << " canonical endpoint links (from " << readCountWithTerminalPair
        << " reads with at least 2 anchors)." << endl;
}
