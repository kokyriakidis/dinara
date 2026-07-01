// AssemblerAbpoaMultiSegmentMSA.cpp
//
// Per-window multi-segment MSA, built on abPOA. The backbone read is laid down
// as a linear partial-order graph (one node per base). Every overlapping read
// is then folded in piece-by-piece between the backbone nodes that correspond
// to the shared anchors it carries, using abPOA's subgraph alignment API
// (abpoa_align_sequence_to_subgraph + abpoa_add_subgraph_alignment).
//
// This replaces the earlier theseus-fork implementation (AssemblerMultiSegment-
// MSA.cpp, excluded from the build). abPOA's banded SIMD DP makes each piece
// alignment O(qlen * band) rather than the fork's scalar WFA O(score^2), which
// removes the quadratic blow-up on reads whose shared anchors are far apart
// (~4.6x faster overall) and avoids the fork's terminal-condition issues.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"

#include <abpoa.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// Extract the base sequence of an oriented read between two marker ordinals as
// abPOA 0-3 base codes (A=0,C=1,G=2,T=3). Mirrors extractSegmentSequence in
// AssemblerMultiSegmentMSA.cpp but emits codes directly (Base::value already is
// the 0-3 encoding abPOA expects).
vector<uint8_t> extractSegmentCodes(
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
    vector<uint8_t> codes;
    codes.reserve(endPos - beginPos);
    for(uint32_t pos = beginPos; pos < endPos; pos++) {
        codes.push_back(reads.getOrientedReadBase(oid, pos).value);
    }
    return codes;
}

} // anonymous namespace


// Build one abPOA multi-segment MSA for a single anchor window using all
// oriented reads that share at least two of the window's backbone anchors.
// Returns true if an MSA was produced, false if the window was skipped.
bool Assembler::runOneWindowAbpoaMultiSegmentMSA(
    const shared_ptr<Shasta2Anchors>& shasta2Anchors,
    const shared_ptr<Shasta2Journeys>& shasta2Journeys,
    const AnchorWindow& window)
{
    const Reads& readsRef = getReads();
    const auto& markersRef = *markers;
    const uint64_t k = assemblerInfo->k;

    const OrientedReadId backboneOid = window.backboneOrientedReadId;
    const auto backboneJourney = (*shasta2Journeys)[backboneOid];

    cout << "testAbpoaMultiSegmentMSA: window " << window.windowId
         << " backbone " << backboneOid
         << " anchors [" << window.backboneBegin << "," << window.backboneEnd << ")"
         << " reads " << window.readIntervals.size() << endl;

    const uint32_t nBackboneAnchors = window.backboneEnd - window.backboneBegin;
    if(nBackboneAnchors < 2) {
        cout << "  window " << window.windowId << " has < 2 anchors, skipping." << endl;
        return false;
    }
    const uint32_t nSegments = nBackboneAnchors - 1;

    // Backbone marker ordinal at each anchor boundary (0..nSegments).
    // anchorOrdinal[bi] is the backbone read's marker ordinal for the anchor at
    // journey position backboneBegin + bi.
    vector<uint32_t> anchorOrdinal(nBackboneAnchors);
    for(uint32_t bi = 0; bi < nBackboneAnchors; bi++) {
        const Shasta2AnchorId anchorId = backboneJourney[window.backboneBegin + bi];
        anchorOrdinal[bi] = shasta2Anchors->getOrdinal(anchorId, backboneOid);
    }

    // Full backbone base sequence as 0-3 codes, plus the base offset of each
    // anchor's MIDPOINT within that sequence.
    //
    // FULL-ANCHOR COVERAGE: anchors are k=50-base markers; using only their
    // midpoint (position + kHalf) as the segment boundary leaves the OUTER kHalf
    // of the two window-terminal anchors with no column at all. To cover the
    // full 50 bases of every anchor, we extend the backbone span by kHalf on
    // each end: it now runs from the START of the first anchor
    // (marker[first].position) to the END of the last anchor
    // (marker[last].position + k), instead of midpoint-to-midpoint. This adds
    // kHalf real columns at each window edge so reads can place the outer halves
    // of their first/last anchors.
    const uint32_t kHalf = uint32_t(k / 2);
    const auto backboneMarkers = markersRef[backboneOid.getValue()];
    const uint32_t backboneBeginPos = backboneMarkers[anchorOrdinal.front()].position;          // first anchor START
    const uint32_t backboneEndPos   = backboneMarkers[anchorOrdinal.back()].position + uint32_t(k); // last anchor END
    if(backboneEndPos <= backboneBeginPos) {
        cout << "  backbone span empty, skipping window." << endl;
        return false;
    }
    vector<uint8_t> backboneCodes;
    backboneCodes.reserve(backboneEndPos - backboneBeginPos);
    for(uint32_t pos = backboneBeginPos; pos < backboneEndPos; pos++) {
        backboneCodes.push_back(readsRef.getOrientedReadBase(backboneOid, pos).value);
    }
    const int backboneLen = static_cast<int>(backboneCodes.size());

    // anchorOffset[bi] = base offset of anchor bi's MIDPOINT into backboneCodes.
    // With the kHalf-extended span, anchorOffset[0] == kHalf (not 0) and
    // anchorOffset[last] == backboneLen - kHalf (not backboneLen).
    vector<int> anchorOffset(nBackboneAnchors);
    for(uint32_t bi = 0; bi < nBackboneAnchors; bi++) {
        const uint32_t pos = backboneMarkers[anchorOrdinal[bi]].position + kHalf;  // midpoint
        anchorOffset[bi] = static_cast<int>(pos - backboneBeginPos);
    }

    cout << "  backbone " << backboneLen << " bases across "
         << nSegments << " segments (full-anchor span, +" << kHalf
         << "bp each end)" << endl;

    // ------------------------------------------------------------------
    // abPOA setup.
    // Scoring chosen to match the theseus path's intent: match=0 implied by
    // POA (abPOA uses positive match), 2/3/1 style mismatch+gap. abPOA needs a
    // positive match score; we use the abPOA defaults (match=2, mismatch=4,
    // affine gap 4/2) which are well tuned for long-read POA. Banding (wb=10)
    // is what keeps each piece alignment O(qlen*band).
    // ------------------------------------------------------------------
    abpoa_t* ab = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    abpt->align_mode = ABPOA_GLOBAL_MODE;
    abpt->gap_mode = ABPOA_AFFINE_GAP;
    abpt->match = 2;
    abpt->mismatch = 4;
    abpt->gap_open1 = 4;
    abpt->gap_ext1 = 2;
    abpt->gap_open2 = 0;       // affine (single-piece) gap
    abpt->gap_ext2 = 0;
    abpt->wb = 10;             // adaptive band; <0 disables banding
    abpt->wf = 0.01;
    abpt->disable_seeding = 1; // we drive the segmentation ourselves
    abpt->progressive_poa = 0;
    abpt->out_msa = 1;         // need RC-MSA output
    abpt->out_cons = 0;
    abpt->ret_cigar = 1;
    abpoa_post_set_para(abpt);  // sets use_read_ids etc. from out_msa

    // Total number of sequences = backbone + all reads sharing >=2 anchors.
    // We need an upper bound up front for read_id bitsets; recount precisely
    // after building boundary hits. Use readIntervals.size() as a safe bound.
    const int totReadBound = static_cast<int>(window.readIntervals.size()) + 1;

    // Seed the backbone as read_id 0. On an empty graph (node_n==2),
    // abpoa_add_subgraph_alignment lays the sequence down as a linear chain and
    // fills qpos_to_node_id with the node id for each base position.
    vector<int> backboneQposToNode(backboneLen, -1);
    {
        abpoa_res_t res;
        res.n_cigar = 0; res.m_cigar = 0; res.graph_cigar = nullptr;
        abpoa_add_subgraph_alignment(
            ab, abpt,
            ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID,
            backboneCodes.data(), nullptr, backboneLen,
            backboneQposToNode.data(), res,
            /* read_id */ 0, /* tot_read_n */ totReadBound,
            /* inc_both_ends */ 1);
    }

    // anchorNode[bi] = abPOA node id of anchor bi's MIDPOINT base.
    // With the kHalf-extended backbone, every anchor midpoint offset lies in
    // [kHalf, backboneLen-kHalf], so ALL anchors map to real interior nodes and
    // are handled uniformly as anchor MATCHes (no SINK special case).
    // backboneNodeAt(off) resolves a backbone base offset to its node id, or to
    // SRC/SINK when the offset falls just before/after the chain (used for the
    // excluded fold boundaries of a read's outer anchor halves).
    vector<int> anchorNode(nBackboneAnchors);
    for(uint32_t bi = 0; bi < nBackboneAnchors; bi++) {
        const int off = anchorOffset[bi];
        anchorNode[bi] = (off >= 0 && off < backboneLen)
            ? backboneQposToNode[off]
            : ABPOA_SINK_NODE_ID;  // defensive; should not occur with extension
    }
    auto backboneNodeAt = [&](int off) -> int {
        if(off < 0) return ABPOA_SRC_NODE_ID;
        if(off >= backboneLen) return ABPOA_SINK_NODE_ID;
        return backboneQposToNode[off];
    };

    // ------------------------------------------------------------------
    // Build read -> backbone boundary hits (identical logic to the theseus
    // path: anchor membership, sort, drop <2, clip to pairwise alignment
    // ordinal range, re-drop <2, sort by base span descending).
    // ------------------------------------------------------------------
    struct BoundaryHit {
        uint32_t boundaryIndex;  // anchor index 0..nSegments
        uint32_t ordinal;        // marker ordinal on this read
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
        if(it->second.size() < 2) { skippedReads++; it = readBoundaryHits.erase(it); }
        else ++it;
    }

    // Clip to pairwise alignment ordinal range.
    const auto& clipTable = getAlignmentTable();
    const ReadId backboneReadId = backboneOid.getReadId();
    for(auto& [readIdValue, hits] : readBoundaryHits) {
        const OrientedReadId oid = OrientedReadId::fromValue(static_cast<ReadId>(readIdValue));
        const ReadId readId = oid.getReadId();
        uint32_t bestFirst = 0, bestLast = 0, bestSpan = 0;
        const auto& aligns = clipTable[oid.getValue()];
        for(uint32_t idx : aligns) {
            const auto& ad = alignmentData[idx];
            ReadId partnerId = (ad.readIds[0] == readId) ? ad.readIds[1] : ad.readIds[0];
            if(partnerId != backboneReadId) continue;
            int targetIdx = (ad.readIds[0] == readId) ? 0 : 1;
            uint32_t firstOrd = ad.info.data[targetIdx].firstOrdinal;
            uint32_t lastOrd  = ad.info.data[targetIdx].lastOrdinal;
            Strand storedStrand = (targetIdx == 0) ? 0 : (ad.isSameStrand ? 0 : 1);
            if(storedStrand != oid.getStrand()) {
                uint32_t mc = ad.info.data[targetIdx].markerCount;
                uint32_t f = mc - 1 - lastOrd;
                uint32_t l = mc - 1 - firstOrd;
                firstOrd = f; lastOrd = l;
            }
            uint32_t span = (lastOrd > firstOrd) ? (lastOrd - firstOrd) : 0;
            if(span > bestSpan) { bestSpan = span; bestFirst = firstOrd; bestLast = lastOrd; }
        }
        if(bestSpan > 0) {
            hits.erase(remove_if(hits.begin(), hits.end(),
                [&](const BoundaryHit& h) {
                    return h.ordinal < bestFirst || h.ordinal > bestLast;
                }), hits.end());
        }
    }
    for(auto it = readBoundaryHits.begin(); it != readBoundaryHits.end(); ) {
        if(it->second.size() < 2) { skippedReads++; it = readBoundaryHits.erase(it); }
        else ++it;
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
         << ", skipped: " << skippedReads << endl;

    // ------------------------------------------------------------------
    // Per-read alignment: fold each read piece (between consecutive shared
    // anchors) into the graph via abPOA subgraph alignment. Each read gets a
    // single read_id; all its pieces accumulate cigars into one whole_res that
    // is added once, so the MSA row is contiguous.
    //
    // abPOA's subgraph API takes EXCLUDED begin/end nodes:
    //   exc_beg | inc_beg ... inc_end | exc_end
    // For a piece spanning anchors a -> b on the backbone, the bases strictly
    // between anchor a and anchor b are aligned. We exclude the node at anchor a
    // (the left pin) and the node at anchor b (the right pin), so beg = node at
    // anchor a, end = node at anchor b.
    // ------------------------------------------------------------------
    uint32_t alignedSegments = 0;
    uint32_t alignedReads = 0;
    double totalAlignTime = 0.0;
    double maxAlignTime = 0.0;
    size_t totalAlignBases = 0;

    vector<double> perReadTime;
    int readSeqId = 1;  // 0 is the backbone

    // abPOA re-runs a full topological sort of the WHOLE graph on every
    // add_subgraph_alignment. Adding each piece separately therefore costs
    // O(graph_nodes) per piece -> O(graph_nodes * pieces), which dominates
    // runtime even though each piece ALIGNMENT is fast. abPOA's own internal
    // aligner avoids this by accumulating all of a read's piece cigars into one
    // whole-read cigar and folding it in with a SINGLE add per read (one sort
    // per read, not per piece). We do the same here.
    //
    // FULL-ANCHOR COVERAGE. Every shared anchor is a k=50-base marker; its
    // MIDPOINT (position + kHalf) is the segment boundary. The read's query span
    // covers the FULL extent of its first and last shared anchors: from the
    // START of the first anchor (rm[firstOrd].position) to the END of the last
    // anchor (rm[lastOrd].position + k). All anchors -- first, interior, last --
    // are emitted as a single-base midpoint MATCH; the outer kHalf of the first
    // and last anchors are aligned as leading/trailing gaps against the
    // kHalf-extended backbone. The stitched cigar is:
    //   <leading gap>          <- outer kHalf of first anchor (before its midpoint)
    //   MATCH(anchorNode[a0])  <- first anchor midpoint
    //   <gap 0>                <- bases between a0 and a1 midpoints
    //   MATCH(anchorNode[a1])
    //   ...
    //   MATCH(anchorNode[am])  <- last anchor midpoint
    //   <trailing gap>         <- outer kHalf of last anchor (after its midpoint)
    //
    // The fold's excluded boundaries (foldBegNode/foldEndNode) are the backbone
    // nodes immediately OUTSIDE the read's covered span (one base before the
    // first covered base, one base after the last). Because those boundaries are
    // now genuinely outside the read, inc_both_ends is 0: membership is recorded
    // only on the read's real covered bases (abPOA records read membership on
    // the SOURCE node of each out-edge; the final synthesized edge records the
    // read's true last base). Each anchor MATCH consumes the read's OWN base, so
    // a sequencing error at an anchor spawns an aligned variant node rather than
    // silently displaying the backbone base. add_subgraph_alignment derives
    // query_id by counting M/I ops, so cigars must cover every query base
    // exactly once, in order.
    auto encodeMatch = [](int nodeId) -> abpoa_cigar_t {
        // CMATCH: node_id << 34 | query_id << 4 | op. query_id is recomputed by
        // the consumer for M ops, so we only need node_id and the op.
        return (static_cast<abpoa_cigar_t>(nodeId) << 34)
             | static_cast<abpoa_cigar_t>(ABPOA_CMATCH);
    };

    const int intK = static_cast<int>(k);

    for(const auto& [baseSpan, readIdValue] : readsBySpan) {
        (void)baseSpan;
        const auto& hits = readBoundaryHits[readIdValue];
        const OrientedReadId oid = OrientedReadId::fromValue(static_cast<ReadId>(readIdValue));
        const auto rm = markersRef[oid.getValue()];

        const uint32_t firstBoundary = hits.front().boundaryIndex;
        const uint32_t lastBoundary  = hits.back().boundaryIndex;
        if(firstBoundary >= nBackboneAnchors || lastBoundary >= nBackboneAnchors) continue;

        // Query span in read coordinates: first-anchor START .. last-anchor END.
        const uint32_t qBegin = rm[hits.front().ordinal].position;          // first anchor start
        const uint32_t qEnd   = rm[hits.back().ordinal].position + intK;    // last anchor end
        if(qEnd <= qBegin) continue;
        const int qlen = static_cast<int>(qEnd - qBegin);

        // Fold boundaries: backbone nodes just OUTSIDE the read's covered span.
        // The read's first covered base aligns to backbone offset
        // (anchorOffset[firstBoundary] - kHalf); its predecessor is the excluded
        // left boundary. Symmetrically on the right. SRC/SINK at the window ends.
        const int firstMidOffset = anchorOffset[firstBoundary];
        const int lastMidOffset  = anchorOffset[lastBoundary];
        const int foldBegNode = backboneNodeAt(firstMidOffset - static_cast<int>(kHalf) - 1);
        const int foldEndNode = backboneNodeAt(lastMidOffset  + static_cast<int>(kHalf) + 1);

        vector<uint8_t> readCodes;
        readCodes.reserve(qlen);
        for(uint32_t pos = qBegin; pos < qEnd; pos++) {
            readCodes.push_back(readsRef.getOrientedReadBase(oid, pos).value);
        }

        vector<abpoa_cigar_t> whole;   // stitched whole-read cigar
        double readTime = 0.0;
        uint32_t readSegments = 0;
        bool readOk = true;

        // Helper to align a gap [begQpos, gapEnd) against subgraph (beg, end)
        // and append its cigars, updating timers.
        auto alignGap = [&](int begNode, int endNode, int begQpos, int gapLen,
                            uint32_t boundaryForLog) {
            if(gapLen <= 0) return;
            abpoa_res_t res;
            res.n_cigar = 0; res.m_cigar = 0; res.graph_cigar = nullptr;
            auto t0 = chrono::steady_clock::now();
            abpoa_align_sequence_to_subgraph(
                ab, abpt, begNode, endNode,
                readCodes.data() + begQpos, gapLen, &res);
            auto t1 = chrono::steady_clock::now();
            double elapsed = chrono::duration<double>(t1 - t0).count();
            for(int ci = 0; ci < res.n_cigar; ci++) whole.push_back(res.graph_cigar[ci]);
            if(res.n_cigar) free(res.graph_cigar);
            totalAlignTime += elapsed;
            readTime += elapsed;
            totalAlignBases += gapLen;
            if(elapsed > maxAlignTime) maxAlignTime = elapsed;
            if(elapsed > 0.1) {
                cout << "  SLOW: read " << oid
                     << " gap near anchor " << boundaryForLog
                     << " seq " << gapLen << " bases took " << elapsed << "s" << endl;
            }
        };

        int begNode = foldBegNode;   // excluded left boundary
        int begQpos = 0;             // query offset of next unaligned base
        for(size_t hi = 0; hi < hits.size(); hi++) {
            const uint32_t boundary = hits[hi].boundaryIndex;
            if(boundary >= nBackboneAnchors) { readOk = false; break; }
            const int pinNode = anchorNode[boundary];

            // Anchor midpoint offset in readCodes.
            const int pinQpos = static_cast<int>(rm[hits[hi].ordinal].position + kHalf) - static_cast<int>(qBegin);
            if(pinQpos < begQpos || pinQpos >= qlen) { readOk = false; break; }

            // Gap of read bases strictly before this anchor midpoint, aligned
            // against (begNode, pinNode). For hi==0 this is the outer kHalf of
            // the first anchor; for hi>0 it is the inter-anchor segment.
            alignGap(begNode, pinNode, begQpos, pinQpos - begQpos, boundary);

            // Anchor midpoint MATCH (read's own base).
            whole.push_back(encodeMatch(pinNode));
            begNode = pinNode;
            begQpos = pinQpos + 1;
            readSegments++;
            alignedSegments++;
        }

        if(readOk) {
            // Trailing gap: outer kHalf of the last anchor, aligned against
            // (lastAnchorNode, foldEndNode).
            alignGap(begNode, foldEndNode, begQpos, qlen - begQpos, lastBoundary);
        }

        if(readOk && readSegments > 0 && !whole.empty()) {
            abpoa_res_t wholeRes;
            wholeRes.n_cigar = static_cast<int>(whole.size());
            wholeRes.m_cigar = static_cast<int>(whole.size());
            wholeRes.graph_cigar = whole.data();
            auto t0 = chrono::steady_clock::now();
            abpoa_add_subgraph_alignment(
                ab, abpt,
                foldBegNode, foldEndNode,
                readCodes.data(), nullptr, qlen,
                nullptr, wholeRes,
                readSeqId, totReadBound, /* inc_both_ends */ 0);
            auto t1 = chrono::steady_clock::now();
            readTime += chrono::duration<double>(t1 - t0).count();
            totalAlignTime += chrono::duration<double>(t1 - t0).count();

            perReadTime.push_back(readTime);
            alignedReads++;
            readSeqId++;
        }
    }

    cout << "  aligned " << alignedReads << " reads ("
         << alignedSegments << " segments), skipped " << skippedReads << endl;

    // Report timing by quartile (reads are already longest-first).
    if(perReadTime.size() >= 4) {
        size_t q = perReadTime.size() / 4;
        for(int qi = 0; qi < 4; qi++) {
            size_t start = qi * q;
            size_t end = (qi == 3) ? perReadTime.size() : (qi + 1) * q;
            double sum = 0;
            for(size_t i = start; i < end; i++) sum += perReadTime[i];
            cout << "  Q" << (qi + 1) << " (reads " << start << "-" << (end - 1)
                 << "): " << sum << "s total, "
                 << (sum / (end - start) * 1000) << "ms/read" << endl;
        }
    }
    cout << "  total align time: " << totalAlignTime << "s"
         << "  avg: " << (alignedSegments > 0 ? totalAlignTime / alignedSegments * 1000 : 0)
         << "ms/seg  max: " << maxAlignTime << "s"
         << "  total bases: " << totalAlignBases << endl;

    // ------------------------------------------------------------------
    // Emit the MSA. We built the graph manually, so abPOA's seq store does not
    // know how many rows exist; set n_seq = backbone + aligned reads. The RC-MSA
    // generator reads per-edge read_id bitsets (populated by add_read_id, on
    // because out_msa=1) to lay each row out. is_rc[]/name[] arrays are
    // pre-sized to CHUNK_READ_N (1024) and zeroed, so this is safe for windows
    // with <=1024 sequences.
    // ------------------------------------------------------------------
    const int totSeq = readSeqId;  // backbone (0) + alignedReads
    ab->abs->n_seq = totSeq;

    {
        const string msaPath = "testAbpoaMultiSegmentMSA_window"
            + to_string(window.windowId) + ".fasta";
        FILE* msaFp = fopen(msaPath.c_str(), "w");
        if(msaFp) {
            // Populate abc (msa_base, msa_len) from the graph, then write it.
            // abpoa_output_rc_msa returns early if msa_len<=0, so the generate
            // step is required for a manually-built graph.
            abpoa_generate_rc_msa(ab, abpt);
            abpoa_output_rc_msa(ab, abpt, msaFp);
            fclose(msaFp);
            cout << "  MSA written to " << msaPath << endl;
        } else {
            cout << "  WARNING: could not open " << msaPath << " for writing" << endl;
        }
    }
    {
        const string gfaPath = "testAbpoaMultiSegmentMSA_window"
            + to_string(window.windowId) + ".gfa";
        FILE* gfaFp = fopen(gfaPath.c_str(), "w");
        if(gfaFp) {
            abpt->out_gfa = 1;
            abpoa_generate_gfa(ab, abpt, gfaFp);
            abpt->out_gfa = 0;
            fclose(gfaFp);
            cout << "  GFA written to " << gfaPath << endl;
        }
    }

    abpoa_free(ab);
    abpoa_free_para(abpt);
    return true;
}


// Driver: build a per-window all-reads abPOA multi-segment MSA for each anchor
// window, mirroring testMultiSegmentMSA so the two engines can be compared.
// Number of windows capped by env DINARA_MSA_MAX_WINDOWS (default 1, 0=all).
void Assembler::testAbpoaMultiSegmentMSA(
    const shared_ptr<Shasta2Anchors>& shasta2Anchors,
    const shared_ptr<Shasta2Journeys>& shasta2Journeys,
    const vector<AnchorWindow>& anchorWindows)
{
    if(anchorWindows.empty()) {
        cout << "testAbpoaMultiSegmentMSA: no windows." << endl;
        return;
    }

    uint64_t maxWindows = 1;
    if(const char* env = getenv("DINARA_MSA_MAX_WINDOWS")) {
        maxWindows = strtoull(env, nullptr, 10);  // 0 = all windows
    }

    cout << "testAbpoaMultiSegmentMSA: " << anchorWindows.size()
         << " windows available";
    if(maxWindows == 0) cout << ", processing all." << endl;
    else cout << ", processing up to " << maxWindows << "." << endl;

    uint64_t processed = 0, produced = 0;
    for(const AnchorWindow& window : anchorWindows) {
        if(maxWindows != 0 && processed >= maxWindows) break;
        processed++;
        if(runOneWindowAbpoaMultiSegmentMSA(shasta2Anchors, shasta2Journeys, window)) {
            produced++;
        }
    }

    cout << "testAbpoaMultiSegmentMSA: produced MSA for " << produced
         << " of " << processed << " processed windows." << endl;
}
