// AssemblerMultiSegmentMSA.cpp
//
// Test function for multi-segment MSA using abPOA's subgraph alignment.
// Picks one anchor window, builds a POA graph from the backbone read's
// concatenated segments, then aligns each overlapping read's segments
// using abpoa_align_sequence_to_subgraph.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"

extern "C" {
#include <abpoa/abpoa.h>
}

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
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

// Convert a character base to abPOA's 0-3 encoding (A=0, C=1, G=2, T=3).
inline uint8_t baseToAbpoa(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;
    }
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

    // Concatenate all backbone segments into one sequence for abPOA.
    string backboneSeq;
    size_t totalBackboneBases = 0;
    for(const auto& s : segmentStrings) {
        totalBackboneBases += s.size();
    }
    backboneSeq.reserve(totalBackboneBases);
    for(const auto& s : segmentStrings) {
        backboneSeq += s;
    }

    // Compute exclusive boundary node IDs for abPOA subgraph alignment.
    // After adding the backbone, the graph has:
    //   node 0 = source, nodes 2..N+1 = bases, node 1 = sink.
    // abpoa_align_sequence_to_subgraph(beg, end) aligns to all nodes
    // between beg and end, both excluded.
    //
    // For segment i, the bases occupy nodes [2+cumBases[i], 2+cumBases[i+1]-1].
    // To align to segments [prev..next-1], we need:
    //   excBeg = node before first base of segment prev
    //   excEnd = node after last base of segment next-1
    //
    // We store one exclusive boundary per anchor boundary:
    //   boundaryNodeIds[0] = source (before all bases)
    //   boundaryNodeIds[i] for 0 < i < nSegments:
    //     We need a node that is AFTER the last base of segment i-1
    //     AND BEFORE the first base of segment i.
    //     In the initial linear graph, the last base of segment i-1
    //     has an edge to the first base of segment i. There's no
    //     intermediate node. So we use the last base of segment i-1
    //     as the exclusive boundary — it will be excluded from the
    //     subgraph of segment i, and included in segment i-1's subgraph.
    //   boundaryNodeIds[nSegments] = sink (after all bases)
    vector<int> boundaryNodeIds(nSegments + 1);
    {
        boundaryNodeIds[0] = ABPOA_SRC_NODE_ID;
        size_t cumBases = 0;
        for(uint32_t i = 0; i < nSegments; i++) {
            cumBases += segmentStrings[i].size();
            if(i + 1 < nSegments) {
                // Last base of segment i = node (2 + cumBases - 1).
                boundaryNodeIds[i + 1] = 2 + int(cumBases) - 1;
            } else {
                boundaryNodeIds[i + 1] = ABPOA_SINK_NODE_ID;
            }
        }
    }

    // Initialize abPOA.
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->align_mode = ABPOA_GLOBAL_MODE;
    abpt->gap_mode = ABPOA_AFFINE_GAP;
    abpt->zdrop = -1;       // disable z-drop
    abpt->end_bonus = 0;
    abpt->wb = -1;          // disable adaptive banding
    abpt->ret_cigar = 1;    // need CIGAR for add_graph_alignment
    abpt->out_msa = 0;
    abpt->out_cons = 0;
    abpt->out_gfa = 0;
    abpt->verbose = ABPOA_NONE_VERBOSE;
    abpoa_post_set_para(abpt);

    abpoa_t *ab = abpoa_init();

    // Add the backbone as the first sequence.
    {
        vector<uint8_t> bbEnc(backboneSeq.size());
        for(size_t i = 0; i < backboneSeq.size(); i++) {
            bbEnc[i] = baseToAbpoa(backboneSeq[i]);
        }
        int bbLen = int(bbEnc.size());

        // Reset graph for the backbone length.
        abpoa_reset(ab, abpt, bbLen);

        // Align backbone to empty graph — this creates the initial linear graph.
        abpoa_res_t res;
        res.graph_cigar = nullptr;
        res.n_cigar = 0;
        res.m_cigar = 0;
        abpoa_align_sequence_to_graph(ab, abpt, bbEnc.data(), bbLen, &res);
        abpoa_add_graph_alignment(ab, abpt, bbEnc.data(), nullptr, bbLen,
            nullptr, res, 0, 1, 1);
        if(res.graph_cigar) free(res.graph_cigar);
    }

    cout << "  abPOA graph initialized with " << ab->abg->node_n
         << " nodes (" << totalBackboneBases << " bases)" << endl;
    cout << flush;

    // Build read -> sorted backbone boundary info from anchor marker intervals.
    struct BoundaryHit {
        uint32_t boundaryIndex;
        uint32_t ordinal;
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

    for(auto& [readId, hits] : readBoundaryHits) {
        sort(hits.begin(), hits.end(),
            [](const BoundaryHit& a, const BoundaryHit& b) {
                return a.boundaryIndex < b.boundaryIndex;
            });
    }

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

    // For each read, align segments between consecutive shared backbone anchors
    // using abPOA subgraph alignment.
    uint32_t alignedSegments = 0;
    uint32_t alignedReads = 0;
    double totalAlignTime = 0.0;
    double maxAlignTime = 0.0;
    uint32_t maxAlignSeg = 0;
    size_t totalAlignBases = 0;
    int readSeqId = 1;  // 0 is the backbone

    for(const auto& [readIdValue, hits] : readBoundaryHits) {
        const OrientedReadId oid = OrientedReadId::fromValue(static_cast<ReadId>(readIdValue));

        uint32_t readSegments = 0;
        for(size_t hi = 0; hi + 1 < hits.size(); hi++) {
            const uint32_t prevBoundary = hits[hi].boundaryIndex;
            const uint32_t nextBoundary = hits[hi + 1].boundaryIndex;

            if(nextBoundary <= prevBoundary) {
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

            // Convert to abPOA encoding.
            vector<uint8_t> queryEnc(readSeq.size());
            for(size_t qi = 0; qi < readSeq.size(); qi++) {
                queryEnc[qi] = baseToAbpoa(readSeq[qi]);
            }

            // Subgraph boundaries (both exclusive).
            const int excBeg = boundaryNodeIds[prevBoundary];
            const int excEnd = boundaryNodeIds[nextBoundary];

            auto t0 = chrono::steady_clock::now();

            abpoa_res_t res;
            res.graph_cigar = nullptr;
            res.n_cigar = 0;
            res.m_cigar = 0;

            abpoa_align_sequence_to_subgraph(ab, abpt, excBeg, excEnd,
                queryEnc.data(), int(queryEnc.size()), &res);

            if(res.n_cigar > 0) {
                abpoa_add_subgraph_alignment(ab, abpt, excBeg, excEnd,
                    queryEnc.data(), nullptr, int(queryEnc.size()),
                    nullptr, res, readSeqId, readSeqId + 1, 0);
                if(res.graph_cigar) free(res.graph_cigar);
            }

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
            readSeqId++;

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
                     << " score " << res.best_score
                     << endl;
            }
        }

        if(readSegments > 0) {
            alignedReads++;
        }
    }

    cout << "  aligned " << alignedReads << " reads (" << alignedSegments << " segments), skipped " << skippedReads << endl;
    cout << "  total align time: " << totalAlignTime << "s"
         << "  avg: " << (alignedSegments > 0 ? totalAlignTime / alignedSegments * 1000 : 0) << "ms/seg"
         << "  max: " << maxAlignTime << "s (seg#" << maxAlignSeg << ")"
         << "  total bases: " << totalAlignBases << endl;

    // Write GFA.
    {
        const string gfaPath = "testMultiSegmentMSA_window" + to_string(window.windowId) + ".gfa";
        FILE* gfaFile = fopen(gfaPath.c_str(), "w");
        if(gfaFile) {
            abpoa_generate_gfa(ab, abpt, gfaFile);
            fclose(gfaFile);
            cout << "  GFA written to " << gfaPath << endl;
        }
    }

    abpoa_free(ab);
    abpoa_free_para(abpt);
}
