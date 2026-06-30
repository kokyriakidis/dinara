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


// Build one multi-segment Theseus MSA for a single anchor window using all
// oriented reads that share at least two of the window's backbone anchors.
// Returns true if an MSA was produced (FASTA + GFA written), false if the
// window was skipped (too few anchors, empty segment, etc.).
// Defined as a member helper below; declared static-like via the Assembler
// method so it can reuse member accessors.
bool Assembler::runOneWindowMultiSegmentMSA(
    const shared_ptr<Shasta2Anchors>& shasta2Anchors,
    const shared_ptr<Shasta2Journeys>& shasta2Journeys,
    const AnchorWindow& window)
{
    const Reads& readsRef = getReads();
    const auto& markersRef = *markers;
    const uint64_t k = assemblerInfo->k;

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
    if(nBackboneAnchors < 2) {
        cout << "  window " << window.windowId << " has < 2 anchors, skipping." << endl;
        return false;
    }
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
            return false;
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
    theseus::Heuristics heuristics;   // defaults: density_drop=off, lag_pruning=off
    vector<theseus::Graph::NodeId> nodeIds;
    theseus::TheseusMSA aligner(penalties, heuristics, segmentViews, nodeIds, 1);

    cout << "  TheseusMSA created with " << nodeIds.size() << " segment nodes." << endl;

    // Build read -> backbone boundary info from anchor membership.
    // For each backbone anchor, look up all reads that contain it and
    // record the boundary index and marker ordinal.
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

    // Sort each read's hits by boundary index.
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

    // Clip boundary hits to pairwise alignment ordinal range.
    // For each read, find its best alignment with the backbone and
    // remove boundary hits outside the alignment's ordinal range.
    // This eliminates spurious hits from transitive closure.
    const auto& clipTable = getAlignmentTable();
    const ReadId backboneReadId = backboneOid.getReadId();
    uint32_t clippedReads = 0;
    for(auto& [readIdValue, hits] : readBoundaryHits) {
        const OrientedReadId oid = OrientedReadId::fromValue(static_cast<ReadId>(readIdValue));
        const ReadId readId = oid.getReadId();

        // Find best alignment with backbone.
        uint32_t bestFirst = 0, bestLast = 0;
        uint32_t bestSpan = 0;
        const auto& aligns = clipTable[oid.getValue()];
        for(uint32_t idx : aligns) {
            const auto& ad = alignmentData[idx];
            ReadId partnerId = (ad.readIds[0] == readId) ? ad.readIds[1] : ad.readIds[0];
            if(partnerId != backboneReadId) continue;

            int targetIdx = (ad.readIds[0] == readId) ? 0 : 1;
            uint32_t firstOrd = ad.info.data[targetIdx].firstOrdinal;
            uint32_t lastOrd  = ad.info.data[targetIdx].lastOrdinal;

            // Flip ordinals if strand mismatch.
            Strand storedStrand = (targetIdx == 0) ? 0
                : (ad.isSameStrand ? 0 : 1);
            if(storedStrand != oid.getStrand()) {
                uint32_t mc = ad.info.data[targetIdx].markerCount;
                uint32_t f = mc - 1 - lastOrd;
                uint32_t l = mc - 1 - firstOrd;
                firstOrd = f;
                lastOrd = l;
            }

            uint32_t span = (lastOrd > firstOrd) ? (lastOrd - firstOrd) : 0;
            if(span > bestSpan) {
                bestSpan = span;
                bestFirst = firstOrd;
                bestLast = lastOrd;
            }
        }

        if(bestSpan > 0) {
            size_t before = hits.size();
            hits.erase(
                std::remove_if(hits.begin(), hits.end(),
                    [&](const BoundaryHit& h) {
                        return h.ordinal < bestFirst || h.ordinal > bestLast;
                    }),
                hits.end());
            if(hits.size() < before) clippedReads++;
        }
    }

    // Re-remove reads that dropped below 2 hits after clipping.
    for(auto it = readBoundaryHits.begin(); it != readBoundaryHits.end(); ) {
        if(it->second.size() < 2) {
            skippedReads++;
            it = readBoundaryHits.erase(it);
        } else {
            ++it;
        }
    }

    // Sort reads by base span descending (longest first).
    vector<pair<uint32_t, uint64_t>> readsBySpan;
    readsBySpan.reserve(readBoundaryHits.size());
    for(const auto& [readIdValue, hits] : readBoundaryHits) {
        const auto readMarkers = markersRef[readIdValue];
        uint32_t firstPos = readMarkers[hits.front().ordinal].position;
        uint32_t lastPos  = readMarkers[hits.back().ordinal].position;
        uint32_t span = (lastPos > firstPos) ? (lastPos - firstPos) : 0;
        readsBySpan.push_back({span, readIdValue});
    }
    sort(readsBySpan.begin(), readsBySpan.end(),
        [](const auto& a, const auto& b) { return a.first > b.first; });

    cout << "  reads with >=2 shared anchors: " << readsBySpan.size()
         << ", skipped: " << skippedReads
         << ", clipped: " << clippedReads << endl;



    // ---------------------------------------------------------------
    // Alignment loop: journey-chained MSA.
    // For each read (longest first), walk its journey interval and
    // find anchors that already have nodes in the graph. Align
    // sub-segments between consecutive known anchors. After each
    // sub-segment, register intermediate journey anchors from the
    // alignment path so subsequent reads can use them as boundaries.
    // ---------------------------------------------------------------
    uint32_t alignedSegments = 0;
    uint32_t alignedReads = 0;
    double totalAlignTime = 0.0;
    double maxAlignTime = 0.0;
    uint32_t maxAlignSeg = 0;
    size_t totalAlignBases = 0;

    vector<double> perReadTime;  // per-read timing
    int readSeqId = 1;  // 0 is the backbone
    vector<string> msaSeqNames;
    msaSeqNames.push_back(to_string(backboneOid.getValue()));
    vector<uint64_t> msaSeqIds;
    msaSeqIds.push_back(backboneOid.getValue());

    for(const auto& [baseSpan, readIdValue] : readsBySpan) {
        const auto& hits = readBoundaryHits[readIdValue];
        const OrientedReadId oid = OrientedReadId::fromValue(static_cast<ReadId>(readIdValue));

        uint32_t readSegments = 0;
        double readTime = 0.0;
        for(size_t hi = 0; hi + 1 < hits.size(); hi++) {
            const uint32_t prevBoundary = hits[hi].boundaryIndex;
            const uint32_t nextBoundary = hits[hi + 1].boundaryIndex;
            if(nextBoundary <= prevBoundary || prevBoundary >= nodeIds.size()) continue;

            const uint32_t prevOrdinal = hits[hi].ordinal;
            const uint32_t nextOrdinal = hits[hi + 1].ordinal;
            if(nextOrdinal <= prevOrdinal) continue;

            string readSeq = extractSegmentSequence(
                readsRef, markersRef, k, oid, prevOrdinal, nextOrdinal);
            if(readSeq.empty()) continue;

            int endNode = (nextBoundary < nodeIds.size())
                ? static_cast<int>(nodeIds[nextBoundary])
                : -1;

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
            readTime += elapsed;
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
        }

        if(readSegments > 0) {
            msaSeqNames.push_back(to_string(oid.getValue()));
            msaSeqIds.push_back(oid.getValue());
            perReadTime.push_back(readTime);
            alignedReads++;
            readSeqId++;
        }
    }

    cout << "  aligned " << alignedReads << " reads (" << alignedSegments << " segments), skipped " << skippedReads << endl;

    // Invariant: every read interval in the window is accounted for as exactly
    // one of: the backbone (the reference, not aligned to itself), an aligned
    // read (>=2 shared anchors, produced >=1 segment), or a skipped read
    // (<2 shared anchors, no inter-anchor segment to align).
    {
        const uint64_t accounted = 1ULL + alignedReads + skippedReads; // +1 backbone
        const uint64_t total = window.readIntervals.size();
        if(accounted != total) {
            cout << "  WARNING: read accounting mismatch: backbone(1) + aligned("
                 << alignedReads << ") + skipped(" << skippedReads << ") = "
                 << accounted << " != readIntervals(" << total << ")" << endl;
        }
        DINARA_ASSERT(accounted == total);
    }

    // Report timing by quartile (reads are already in longest-first order).
    if(perReadTime.size() >= 4) {
        size_t q = perReadTime.size() / 4;
        for(int qi = 0; qi < 4; qi++) {
            size_t start = qi * q;
            size_t end = (qi == 3) ? perReadTime.size() : (qi + 1) * q;
            double sum = 0;
            for(size_t i = start; i < end; i++) sum += perReadTime[i];
            cout << "  Q" << (qi + 1) << " (reads " << start << "-" << (end - 1)
                 << "): " << sum << "s total, " << (sum / (end - start) * 1000) << "ms/read" << endl;
        }
    }
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
        // include_consensus = false: we only want the aligned input rows
        // (backbone + reads), not the majority-voting consensus row.
        aligner.print_as_msa(msaFile, readSeqId - 1, &msaSeqNames,
            /* include_consensus = */ false);
        cout << "  MSA written to " << msaPath << endl;
    }
    {
        const string gfaPath = "testMultiSegmentMSA_window" + to_string(window.windowId) + ".gfa";
        ofstream gfaFile(gfaPath);
        aligner.print_as_gfa(gfaFile);
        cout << "  GFA written to " << gfaPath << endl;
    }
    return true;
}


// Driver: build a per-window all-reads multi-segment MSA for each anchor window.
// For each window it constructs one TheseusMSA seeded by the backbone read and
// aligns every read that shares >=2 of the window's anchors (ends-free), then
// writes the window's MSA (FASTA) and POA graph (GFA).
//
// The number of windows processed is capped by the DINARA_MSA_MAX_WINDOWS
// environment variable (default: 1) to keep the test tractable; set it to 0 to
// process all windows.
void Assembler::testMultiSegmentMSA(
    const shared_ptr<Shasta2Anchors>& shasta2Anchors,
    const shared_ptr<Shasta2Journeys>& shasta2Journeys,
    const vector<AnchorWindow>& anchorWindows)
{
    if(anchorWindows.empty()) {
        cout << "testMultiSegmentMSA: no windows." << endl;
        return;
    }

    uint64_t maxWindows = 1;
    if(const char* env = getenv("DINARA_MSA_MAX_WINDOWS")) {
        maxWindows = strtoull(env, nullptr, 10);  // 0 = all windows
    }

    cout << "testMultiSegmentMSA: " << anchorWindows.size() << " windows available";
    if(maxWindows == 0) {
        cout << ", processing all." << endl;
    } else {
        cout << ", processing up to " << maxWindows << " (set DINARA_MSA_MAX_WINDOWS=0 for all)." << endl;
    }

    uint64_t processed = 0;
    uint64_t produced = 0;
    for(const AnchorWindow& window : anchorWindows) {
        if(maxWindows != 0 && processed >= maxWindows) break;
        processed++;
        if(runOneWindowMultiSegmentMSA(shasta2Anchors, shasta2Journeys, window)) {
            produced++;
        }
    }

    cout << "testMultiSegmentMSA: produced MSA for " << produced
         << " of " << processed << " processed windows." << endl;
}
