// AssemblerMultiSegmentMSA.cpp
//
// Test function for the multi-segment TheseusMSA constructor and align_from.
// Picks one anchor window, builds a multi-segment POA graph from the backbone
// read's inter-anchor segments, then aligns each overlapping read's segments
// using align_from.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"

#include <theseus/heuristics.h>
#include <theseus/penalties.h>
#include <theseus/theseus_msa_aligner.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// Extract the base sequence of an oriented read between two marker ordinals.
// Returns the sequence from the midpoint of ordinalA to the midpoint of ordinalB.
string extractSegmentSequence(
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t k,
    OrientedReadId oid,
    uint32_t ordinalA,
    uint32_t ordinalB)
{
    if(ordinalA >= ordinalB) {
        return {};
    }
    const auto readMarkers = markers[oid.getValue()];
    if(ordinalB >= readMarkers.size()) {
        return {};
    }

    const uint32_t kHalf = uint32_t(k / 2);
    const uint32_t beginPos = readMarkers[ordinalA].position + kHalf;
    const uint32_t endPos   = readMarkers[ordinalB].position + kHalf;
    if(endPos <= beginPos) {
        return {};
    }

    string seq;
    seq.reserve(endPos - beginPos);
    for(uint32_t pos = beginPos; pos < endPos; pos++) {
        seq.push_back(reads.getOrientedReadBase(oid, pos).character());
    }
    return seq;
}

} // anonymous namespace


void Assembler::testMultiSegmentMSA(
    const shared_ptr<Shasta2Anchors>& shasta2Anchors,
    const shared_ptr<Shasta2Journeys>& shasta2Journeys,
    const vector<AnchorWindow>& anchorWindows)
{
    const Reads& readsRef = getReads();
    const auto& markersRef = *markers;
    const uint64_t k = assemblerInfo->k;

    if(anchorWindows.empty()) {
        cout << "testMultiSegmentMSA: no windows." << endl;
        return;
    }

    const AnchorWindow& window = anchorWindows.front();
    const OrientedReadId backboneOid = window.backboneOrientedReadId;
    const auto backboneJourney = (*shasta2Journeys)[backboneOid];

    cout << "testMultiSegmentMSA: window " << window.windowId
         << " backbone " << backboneOid
         << " anchors [" << window.backboneBegin << "," << window.backboneEnd << ")"
         << " reads " << window.readIntervals.size() << endl;

    // Build backbone segments: one segment between each pair of consecutive anchors.
    const uint32_t nBackboneAnchors = window.backboneEnd - window.backboneBegin;
    const uint32_t nSegments = nBackboneAnchors - 1;

    vector<string> segmentStrings;
    segmentStrings.reserve(nSegments);

    // Map from journey position to segment index.
    // Segment i spans from journey position (backboneBegin + i) to (backboneBegin + i + 1).
    for(uint32_t i = 0; i < nSegments; i++) {
        const uint32_t journeyPosLeft  = window.backboneBegin + i;
        const uint32_t journeyPosRight = window.backboneBegin + i + 1;

        const Shasta2AnchorId leftAnchorId  = backboneJourney[journeyPosLeft];
        const Shasta2AnchorId rightAnchorId = backboneJourney[journeyPosRight];

        const uint32_t leftOrdinal  = shasta2Anchors->getOrdinal(leftAnchorId, backboneOid);
        const uint32_t rightOrdinal = shasta2Anchors->getOrdinal(rightAnchorId, backboneOid);

        string seg = extractSegmentSequence(readsRef, markersRef, k, backboneOid, leftOrdinal, rightOrdinal);
        if(seg.empty()) {
            cout << "  segment " << i << " is empty (ordinals " << leftOrdinal
                 << "->" << rightOrdinal << "), skipping window." << endl;
            return;
        }
        segmentStrings.push_back(std::move(seg));
    }

    cout << "  " << nSegments << " backbone segments:";
    for(uint32_t i = 0; i < nSegments; i++) {
        cout << " " << segmentStrings[i].size();
    }
    cout << " bases" << endl;

    // Build string_view vector for the TheseusMSA constructor.
    vector<string_view> segmentViews;
    segmentViews.reserve(nSegments);
    for(const auto& s : segmentStrings) {
        segmentViews.push_back(s);
    }

    // Create the multi-segment MSA.
    theseus::Penalties penalties(0, 2, 3, 1);
    theseus::Heuristics heuristics(false, false);
    vector<theseus::Graph::NodeId> nodeIds;
    theseus::TheseusMSA aligner(penalties, heuristics, segmentViews, nodeIds, 1);

    cout << "  TheseusMSA created with " << nodeIds.size() << " segment nodes." << endl;

    // For each non-backbone read in the window, find which backbone anchors
    // it shares, extract the corresponding sequence, and align_from.
    uint32_t alignedCount = 0;
    uint32_t skippedCount = 0;

    for(const auto& readInterval : window.readIntervals) {
        const OrientedReadId oid = readInterval.orientedReadId;
        if(oid == backboneOid) {
            continue;
        }

        const auto readJourney = (*shasta2Journeys)[oid];
        if(readJourney.empty()) {
            skippedCount++;
            continue;
        }

        // Find the first and last backbone anchor that this read shares.
        // The read's journey positions [readInterval.begin, readInterval.end)
        // correspond to anchors in the read's journey. We need to find which
        // of these anchors are also in the backbone's journey, and map them
        // to segment indices.

        // Build a map from backbone anchor ID to segment boundary index.
        // Backbone anchor at journey position (backboneBegin + i) is boundary i.
        // Segment i spans boundaries i to i+1.
        // We need to find the read's entry boundary and exit boundary.

        uint32_t readEntryBoundary = UINT32_MAX;
        uint32_t readExitBoundary = 0;

        for(uint32_t jp = readInterval.begin; jp < readInterval.end; jp++) {
            if(jp >= readJourney.size()) break;
            const Shasta2AnchorId anchorId = readJourney[jp];

            // Check if this anchor is one of the backbone's anchors.
            for(uint32_t bi = 0; bi <= nSegments; bi++) {
                const uint32_t backboneJp = window.backboneBegin + bi;
                if(backboneJp < backboneJourney.size() && backboneJourney[backboneJp] == anchorId) {
                    readEntryBoundary = min(readEntryBoundary, bi);
                    readExitBoundary  = max(readExitBoundary, bi);
                    break;
                }
            }
        }

        if(readEntryBoundary >= readExitBoundary || readEntryBoundary == UINT32_MAX) {
            skippedCount++;
            continue;
        }

        // Extract the read's sequence between its entry and exit anchors.
        const Shasta2AnchorId entryAnchorId = backboneJourney[window.backboneBegin + readEntryBoundary];
        const Shasta2AnchorId exitAnchorId  = backboneJourney[window.backboneBegin + readExitBoundary];

        const uint32_t entryOrdinal = shasta2Anchors->getOrdinal(entryAnchorId, oid);
        const uint32_t exitOrdinal  = shasta2Anchors->getOrdinal(exitAnchorId, oid);

        if(entryOrdinal == invalid<uint32_t> || exitOrdinal == invalid<uint32_t>) {
            skippedCount++;
            continue;
        }

        string readSeq = extractSegmentSequence(readsRef, markersRef, k, oid, entryOrdinal, exitOrdinal);
        if(readSeq.empty()) {
            skippedCount++;
            continue;
        }

        // Align from the entry boundary's segment node.
        // readEntryBoundary is the boundary index; the segment starting at
        // that boundary is segment index readEntryBoundary, whose node ID
        // is nodeIds[readEntryBoundary].
        if(readEntryBoundary >= nodeIds.size()) {
            skippedCount++;
            continue;
        }

        auto alignment = aligner.align_from(
            readSeq,
            nodeIds[readEntryBoundary],
            1,     // weight
            true,  // is_ends_free (read may not reach the last segment)
            0);    // start_offset

        alignedCount++;
        if(alignedCount <= 5) {
            cout << "  read " << oid
                 << " boundaries [" << readEntryBoundary << "," << readExitBoundary << "]"
                 << " seq " << readSeq.size() << " bases"
                 << " score " << alignment.compute_affine_gap_score(penalties)
                 << endl;
        }
    }

    cout << "  aligned " << alignedCount << " reads, skipped " << skippedCount << endl;

    // Write the MSA and GFA to files.
    {
        const string msaPath = "testMultiSegmentMSA_window" + to_string(window.windowId) + ".fasta";
        ofstream msaFile(msaPath);
        aligner.print_as_msa(msaFile);
        cout << "  MSA written to " << msaPath << endl;
    }
    {
        const string gfaPath = "testMultiSegmentMSA_window" + to_string(window.windowId) + ".gfa";
        ofstream gfaFile(gfaPath);
        aligner.print_as_gfa(gfaFile);
        cout << "  GFA written to " << gfaPath << endl;
    }
}
