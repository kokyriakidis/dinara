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
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <string_view>
#include <unordered_map>
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

    // Report memory before MSA.
    {
        ifstream procStatus("/proc/self/status");
        string line;
        while(getline(procStatus, line)) {
            if(line.find("VmRSS") == 0) {
                cout << "  before MSA: " << line << endl;
            }
        }
    }

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

    // Build read -> sorted backbone boundary info directly from anchor marker
    // intervals. For each backbone boundary anchor, look up all oriented reads
    // that contain it and record the boundary index and marker ordinal.
    // This avoids walking each read's full journey.
    struct BoundaryHit {
        uint32_t boundaryIndex;
        uint32_t ordinal;  // marker ordinal on this read
    };
    unordered_map<uint64_t, vector<BoundaryHit>> readBoundaryHits;

    for(uint32_t bi = 0; bi <= nSegments; bi++) {
        const uint32_t bjp = window.backboneBegin + bi;
        if(bjp >= backboneJourney.size()) break;
        const Shasta2AnchorId anchorId = backboneJourney[bjp];
        const auto anchor = (*shasta2Anchors)[anchorId];
        for(const auto& info : anchor) {
            if(info.orientedReadId == backboneOid) continue;
            readBoundaryHits[info.orientedReadId.getValue()].push_back(
                {bi, info.ordinal});
        }
    }

    // Sort each read's hits by boundary index (they may arrive out of order
    // if the same read appears at multiple backbone anchors).
    for(auto& [readId, hits] : readBoundaryHits) {
        sort(hits.begin(), hits.end(),
            [](const BoundaryHit& a, const BoundaryHit& b) {
                return a.boundaryIndex < b.boundaryIndex;
            });
    }

    // Remove reads with fewer than 2 boundary hits.
    uint32_t skippedReads = 0;
    for(auto it = readBoundaryHits.begin(); it != readBoundaryHits.end(); ) {
        if(it->second.size() < 2) {
            skippedReads++;
            it = readBoundaryHits.erase(it);
        } else {
            ++it;
        }
    }

    cout << "  non-backbone reads with >=2 shared anchors: " << readBoundaryHits.size()
         << ", skipped: " << skippedReads << endl;

    // Prefix sum of backbone segment lengths for computing base spans.
    vector<size_t> segPrefixSum(nSegments + 1, 0);
    for(uint32_t i = 0; i < nSegments; i++) {
        segPrefixSum[i + 1] = segPrefixSum[i] + segmentStrings[i].size();
    }

    // Sort reads by backbone base span (descending) so longer-spanning
    // reads are added to the POA graph first.
    // Use a flat vector of (baseSpan, readIdValue) pairs — sort compares
    // the span directly without any hash lookups.
    vector<pair<size_t, uint64_t>> readsBySpan;
    readsBySpan.reserve(readBoundaryHits.size());
    for(const auto& [readIdValue, hits] : readBoundaryHits) {
        size_t span = segPrefixSum[hits.back().boundaryIndex]
                    - segPrefixSum[hits.front().boundaryIndex];
        readsBySpan.push_back({span, readIdValue});
    }
    sort(readsBySpan.begin(), readsBySpan.end(),
        [](const auto& a, const auto& b) { return a.first > b.first; });

    // For each read, align segments between consecutive shared backbone anchors.
    uint32_t alignedSegments = 0;
    uint32_t alignedReads = 0;
    double totalAlignTime = 0.0;
    double maxAlignTime = 0.0;
    uint32_t maxAlignSeg = 0;
    size_t totalAlignBases = 0;
    int readSeqId = 1;  // 0 is the backbone
    vector<string> msaSeqNames;
    msaSeqNames.push_back(to_string(backboneOid.getValue()));
    vector<uint64_t> msaSeqIds;
    msaSeqIds.push_back(backboneOid.getValue());

    for(const auto& [baseSpan, readIdValue] : readsBySpan) {
        const auto& hits = readBoundaryHits[readIdValue];
        const OrientedReadId oid = OrientedReadId::fromValue(static_cast<ReadId>(readIdValue));

        uint32_t readSegments = 0;
        for(size_t hi = 0; hi + 1 < hits.size(); hi++) {
            const uint32_t prevBoundary = hits[hi].boundaryIndex;
            const uint32_t nextBoundary = hits[hi + 1].boundaryIndex;

            if(nextBoundary <= prevBoundary || prevBoundary >= nodeIds.size()) {
                continue;
            }

            const uint32_t prevOrdinal = hits[hi].ordinal;
            const uint32_t nextOrdinal = hits[hi + 1].ordinal;

            if(nextOrdinal <= prevOrdinal) {
                continue;
            }

            string readSeq = extractSegmentSequence(readsRef, markersRef, k, oid, prevOrdinal, nextOrdinal);
            if(readSeq.empty()) {
                continue;
            }

            // Pass end_node to scope the alignment to the subgraph
            // between the two boundary nodes.
            int endNode = (nextBoundary < nodeIds.size())
                ? static_cast<int>(nodeIds[nextBoundary])
                : -1;  // -1 = sink

            auto t0 = chrono::steady_clock::now();
            auto alignment = aligner.align_from(
                readSeq,
                nodeIds[prevBoundary],
                1,     // weight
                true,  // is_ends_free
                0,     // start_offset
                endNode,
                readSeqId);
            auto t1 = chrono::steady_clock::now();
            double elapsed = chrono::duration<double>(t1 - t0).count();
            totalAlignTime += elapsed;
            totalAlignBases += readSeq.size();
            if(elapsed > maxAlignTime) {
                maxAlignTime = elapsed;
                maxAlignSeg = alignedSegments;
            }

            readSegments++;
            alignedSegments++;

            if(elapsed > 0.1) {
                cout << "  SLOW: read " << oid
                     << " boundaries [" << prevBoundary << "," << nextBoundary << "]"
                     << " seq " << readSeq.size() << " bases"
                     << " took " << elapsed << "s" << endl;
            }

            if(alignedReads < 5 && readSegments == 1) {
                cout << "  read " << oid
                     << " matches=" << hits.size()
                     << " firstSeg boundaries [" << prevBoundary << "," << nextBoundary << "]"
                     << " seq " << readSeq.size() << " bases"
                     << " score " << alignment.compute_affine_gap_score(penalties)
                     << endl;
            }
        }

        if(readSegments > 0) {
            msaSeqNames.push_back(to_string(oid.getValue()));
            msaSeqIds.push_back(oid.getValue());
            alignedReads++;
            readSeqId++;
        }
    }

    cout << "  aligned " << alignedReads << " reads (" << alignedSegments << " segments), skipped " << skippedReads << endl;
    cout << "  total align time: " << totalAlignTime << "s"
         << "  avg: " << (alignedSegments > 0 ? totalAlignTime / alignedSegments * 1000 : 0) << "ms/seg"
         << "  max: " << maxAlignTime << "s (seg#" << maxAlignSeg << ")"
         << "  total bases: " << totalAlignBases << endl;

    // Report memory usage.
    {
        ifstream procStatus("/proc/self/status");
        string line;
        while(getline(procStatus, line)) {
            if(line.find("VmRSS") == 0 || line.find("VmPeak") == 0) {
                cout << "  " << line << endl;
            }
        }
    }

    // Write the MSA and GFA to files.
    {
        const string msaPath = "testMultiSegmentMSA_window" + to_string(window.windowId) + ".fasta";
        ofstream msaFile(msaPath);
        aligner.print_as_msa(msaFile, readSeqId - 1, &msaSeqNames);
        cout << "  MSA written to " << msaPath << endl;
    }
    {
        const string gfaPath = "testMultiSegmentMSA_window" + to_string(window.windowId) + ".gfa";
        ofstream gfaFile(gfaPath);
        aligner.print_as_gfa(gfaFile);
        cout << "  GFA written to " << gfaPath << endl;
    }
}
