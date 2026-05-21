/// @file PhasingKmeansAlign.cpp
/// @brief Noisy-region abPOA MSA for k-means phasing refinement.
///
/// See PhasingKmeansNoisyMsa.md for the full implementation plan.

#include "PhasingKmeansAlign.hpp"
#include "Assembler.hpp"
#include "Marker.hpp"
#include "Reads.hpp"

#include "abpoa/abpoa.h"

#include <algorithm>
#include <iostream>

using namespace dinara;
using namespace std;

// ============================================================================
// Helpers
// ============================================================================

/// Convert a character base to 0123 encoding (A=0, C=1, G=2, T=3).
static inline uint8_t baseToInt(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;
    }
}

/// Extract the base sequence of an oriented read between two marker ordinals.
/// Returns the sequence from the midpoint of ordinalA to the midpoint of ordinalB,
/// encoded as 0123 (ACGT) for abPOA.
static vector<uint8_t> extractSegmentSeq0123(
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t k,
    OrientedReadId oid,
    uint32_t ordinalA,
    uint32_t ordinalB)
{
    vector<uint8_t> seq;
    if(ordinalA >= ordinalB) return seq;

    const auto readMarkers = markers[oid.getValue()];
    if(ordinalB >= readMarkers.size()) return seq;

    const uint32_t kHalf = uint32_t(k / 2);
    const uint32_t beginPos = readMarkers[ordinalA].position + kHalf;
    const uint32_t endPos   = readMarkers[ordinalB].position + kHalf;
    if(endPos <= beginPos) return seq;

    seq.reserve(endPos - beginPos);
    for(uint32_t pos = beginPos; pos < endPos; pos++) {
        seq.push_back(baseToInt(reads.getOrientedReadBase(oid, pos).character()));
    }
    return seq;
}

// ============================================================================
// collectNoisyReadInfo
// ============================================================================

KmNoisyReadInfo dinara::collectNoisyReadInfo(
    const Assembler& assembler,
    const KmScratchpad& scratch,
    uint32_t regStart, uint32_t regEnd)
{
    KmNoisyReadInfo info;

    const auto& markersRef = *assembler.markers;
    const auto& alignments = assembler.alignmentCandidatesAlignmentsData.alignments;
    const auto& alignmentDataRef = assembler.alignmentData;

    if(scratch.overlaps.empty()) return info;

    // Determine backbone OrientedReadId (always strand 0).
    const auto& firstOv = scratch.overlaps[0];
    const auto& firstAd = alignmentDataRef[firstOv.alignmentId];
    ReadId backboneReadId = firstOv.queryIsRead0
        ? firstAd.readIds[0] : firstAd.readIds[1];
    const OrientedReadId backboneOid(backboneReadId, 0);
    const auto backboneMarkers = markersRef[backboneOid.getValue()];
    const uint32_t bbMarkerCount = uint32_t(backboneMarkers.size());

    const Reads& readsRef = assembler.getReads();
    const uint64_t k = assembler.assemblerInfo->k;
    const uint32_t kHalf = uint32_t(k / 2);

    // ---------------------------------------------------------------
    // Pass 1: For each overlap, find flanking backbone ordinals around
    // the noisy region. Track the widest envelope.
    // ---------------------------------------------------------------

    struct PerReadInfo {
        int overlapIndex;
        int hap;
        uint32_t leftBbOrd;   // backbone ordinal at left anchor
        uint32_t rightBbOrd;  // backbone ordinal at right anchor
        uint32_t leftRdOrd;   // read ordinal at left anchor
        uint32_t rightRdOrd;  // read ordinal at right anchor
        OrientedReadId readOid;
    };
    vector<PerReadInfo> perRead;

    uint32_t globalLeftOrd = bbMarkerCount;
    uint32_t globalRightOrd = 0;

    const uint32_t numOv = uint32_t(scratch.overlaps.size());
    for(uint32_t oi = 0; oi < numOv; oi++) {
        const auto& ov = scratch.overlaps[oi];
        if(ov.alignmentId >= alignments.size()) continue;
        const auto& aln = alignments[ov.alignmentId];
        if(aln.ordinals.empty()) continue;
        const auto& ad = alignmentDataRef[ov.alignmentId];

        // Determine which ordinal slot is backbone, which is the read.
        int bbSlot, rdSlot;
        OrientedReadId readOid;
        if(ov.queryIsRead0) {
            bbSlot = 0; rdSlot = 1;
            readOid = OrientedReadId(ad.readIds[1], ad.isSameStrand ? 0 : 1);
        } else {
            bbSlot = 1; rdSlot = 0;
            readOid = OrientedReadId(ad.readIds[0], 0);
        }

        const bool bbNeedsFlip = (bbSlot == 1) && !ad.isSameStrand;

        // Build a list of (flipped backbone ordinal, j) pairs sorted by
        // backbone position, so the anchor search works regardless of
        // whether bbNeedsFlip reverses the ordinal order.
        struct OrdPos { uint32_t bbOrd; uint32_t bbPos; int j; };
        vector<OrdPos> ordPositions;
        ordPositions.reserve(aln.ordinals.size());
        for(int j = 0; j < int(aln.ordinals.size()); j++) {
            uint32_t bbOrd = aln.ordinals[j][bbSlot];
            if(bbNeedsFlip) bbOrd = bbMarkerCount - 1 - bbOrd;
            if(bbOrd >= bbMarkerCount) continue;
            uint32_t bbPos = backboneMarkers[bbOrd].position;
            ordPositions.push_back({bbOrd, bbPos, j});
        }
        sort(ordPositions.begin(), ordPositions.end(),
             [](const OrdPos& a, const OrdPos& b) { return a.bbPos < b.bbPos; });

        // Find left anchor (last entry with bbPos <= regStart)
        // and right anchor (first entry with bbPos >= regEnd).
        int leftIdx = -1, rightIdx = -1;
        for(int p = 0; p < int(ordPositions.size()); p++) {
            if(ordPositions[p].bbPos <= regStart) leftIdx = p;
            if(ordPositions[p].bbPos >= regEnd && rightIdx == -1) rightIdx = p;
        }

        if(leftIdx < 0 || rightIdx < 0 || leftIdx >= rightIdx) continue;

        // Map back to original ordinal array indices.
        const int leftJ  = ordPositions[leftIdx].j;
        const int rightJ = ordPositions[rightIdx].j;

        // Use the already-flipped backbone ordinals from ordPositions.
        uint32_t leftBbOrd  = ordPositions[leftIdx].bbOrd;
        uint32_t rightBbOrd = ordPositions[rightIdx].bbOrd;
        // leftBbOrd < rightBbOrd guaranteed by sorted order.

        uint32_t leftRdOrd  = aln.ordinals[leftJ][rdSlot];
        uint32_t rightRdOrd = aln.ordinals[rightJ][rdSlot];
        if(leftRdOrd > rightRdOrd) swap(leftRdOrd, rightRdOrd);

        if(leftBbOrd < globalLeftOrd) globalLeftOrd = leftBbOrd;
        if(rightBbOrd > globalRightOrd) globalRightOrd = rightBbOrd;

        perRead.push_back({int(oi), ov.hap, leftBbOrd, rightBbOrd,
                           leftRdOrd, rightRdOrd, readOid});
    }

    if(perRead.empty() || globalLeftOrd >= globalRightOrd) return info;

    // ---------------------------------------------------------------
    // Extract backbone sequence for the widest envelope.
    // ---------------------------------------------------------------

    const uint32_t bbStartPos = backboneMarkers[globalLeftOrd].position + kHalf;
    const uint32_t bbEndPos   = backboneMarkers[globalRightOrd].position + kHalf;
    if(bbEndPos <= bbStartPos) return info;

    const int bbLen = int(bbEndPos - bbStartPos);
    info.backboneSeq.reserve(bbLen);
    for(uint32_t pos = bbStartPos; pos < bbEndPos; pos++) {
        info.backboneSeq.push_back(
            baseToInt(readsRef.getOrientedReadBase(backboneOid, pos).character()));
    }
    info.backboneLen = bbLen;
    info.backboneStartPos = bbStartPos;

    // ---------------------------------------------------------------
    // Pass 2: Extract per-read sequences and compute backbone-relative
    // positions for subgraph alignment.
    // ---------------------------------------------------------------

    for(auto& pr : perRead) {
        // Extract read sequence between its flanking markers.
        vector<uint8_t> readSeq = extractSegmentSeq0123(
            readsRef, markersRef, k, pr.readOid, pr.leftRdOrd, pr.rightRdOrd);
        if(readSeq.empty()) continue;

        // Compute backbone-relative positions for this read's flanking markers.
        const uint32_t readBbLeftPos  = backboneMarkers[pr.leftBbOrd].position + kHalf;
        const uint32_t readBbRightPos = backboneMarkers[pr.rightBbOrd].position + kHalf;
        const int relLeft  = int(readBbLeftPos  - bbStartPos);
        const int relRight = int(readBbRightPos - bbStartPos);
        if(relLeft < 0 || relRight > bbLen || relLeft >= relRight) continue;

        info.overlapIndices.push_back(pr.overlapIndex);
        info.haps.push_back(pr.hap);
        info.seqs.push_back(move(readSeq));
        info.seqLens.push_back(int(info.seqs.back().size()));
        info.bbLeftPos.push_back(relLeft);
        info.bbRightPos.push_back(relRight);
        info.orientedReadIdValues.push_back(pr.readOid.getValue());
    }

    info.nReads = int(info.overlapIndices.size());
    return info;
}

// ============================================================================
// hasBothHaplotypes
// ============================================================================

bool dinara::hasBothHaplotypes(const KmNoisyReadInfo& info, int minHapReads) {
    int hap1 = 0, hap2 = 0;
    for(int i = 0; i < info.nReads; i++) {
        if(info.haps[i] == 1) hap1++;
        else if(info.haps[i] == 2) hap2++;
    }
    return hap1 >= minHapReads && hap2 >= minHapReads;
}

// ============================================================================
// abpoaMsaRun
// ============================================================================

int dinara::abpoaMsaRun(
    const KmNoisyReadInfo& info,
    const vector<int>& readIndices,
    int maxNCons,
    const KmNoisyMsaOptions& opts,
    bool includeBackboneInReads,
    array<KmAbpoaMsaResult, 2>& resultsOut)
{
    resultsOut = {};
    if(readIndices.empty()) return 0;

    const int nReads = int(readIndices.size());
    // Total sequences = backbone (read 0) + nReads.
    const int totalSeqs = 1 + nReads;

    abpoa_t* ab = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    abpt->cons_algrm = ABPOA_MF;
    abpt->sub_aln = 1;
    abpt->inc_path_score = 1;
    abpt->out_cons = 1;
    abpt->out_msa  = 1;
    abpt->max_n_cons = maxNCons;
    abpt->min_freq   = opts.minFreq;
    abpt->match      = opts.match;
    abpt->mismatch   = opts.mismatch;
    abpt->gap_open1  = opts.gapOpen1;
    abpt->gap_ext1   = opts.gapExt1;
    abpt->gap_open2  = opts.gapOpen2;
    abpt->gap_ext2   = opts.gapExt2;
    abpoa_post_set_para(abpt);

    ab->abs->n_seq = totalSeqs;

    // --- Read 0: seed with backbone sequence ---
    // On an empty graph (node_n == 2), abpoa_align_sequence_to_subgraph
    // returns -1 and res.n_cigar stays 0.  abpoa_add_subgraph_alignment
    // detects the empty graph and calls abpoa_add_graph_sequence internally,
    // creating nodes 2..L+1 for a backbone of length L.
    {
        abpoa_res_t res{};
        res.graph_cigar = nullptr;
        res.n_cigar = 0;

        abpoa_align_sequence_to_subgraph(
            ab, abpt, 0, 1,
            const_cast<uint8_t*>(info.backboneSeq.data()),
            info.backboneLen, &res);
        abpoa_add_subgraph_alignment(
            ab, abpt, 0, 1,
            const_cast<uint8_t*>(info.backboneSeq.data()),
            nullptr, info.backboneLen, nullptr,
            res, 0, totalSeqs, 0);
        if(res.n_cigar) free(res.graph_cigar);
    }

    // --- Reads 1..nReads: align to subgraph between flanking markers ---
    // After seeding, backbone position p (0-based) → node ID p + 2.
    // This mapping is stable — subsequent reads add new nodes with higher IDs.
    for(int i = 0; i < nReads; i++) {
        const int ri = readIndices[i];
        const int leftPos  = info.bbLeftPos[ri];
        const int rightPos = info.bbRightPos[ri];

        // Inclusive node IDs for the flanking markers.
        // bbLeftPos is the start of the read's span (inclusive).
        // bbRightPos is the end (exclusive), so the last included base
        // is at bbRightPos - 1.
        const int incBeg = leftPos + 2;
        const int incEnd = rightPos + 2 - 1;
        if(incBeg > incEnd) continue;

        int excBeg = 0, excEnd = 1;
        abpoa_subgraph_nodes(ab, abpt, incBeg, incEnd, &excBeg, &excEnd);

        abpoa_res_t res{};
        res.graph_cigar = nullptr;
        res.n_cigar = 0;

        abpoa_align_sequence_to_subgraph(
            ab, abpt, excBeg, excEnd,
            const_cast<uint8_t*>(info.seqs[ri].data()),
            info.seqLens[ri], &res);
        abpoa_add_subgraph_alignment(
            ab, abpt, excBeg, excEnd,
            const_cast<uint8_t*>(info.seqs[ri].data()),
            nullptr, info.seqLens[ri], nullptr,
            res, i + 1, totalSeqs, 0);
        if(res.n_cigar) free(res.graph_cigar);
    }

    // --- Generate consensus + MSA ---
    abpoa_output(ab, abpt, nullptr);
    abpoa_cons_t* abc = ab->abc;
    int nCons = 0;

    if(abc->n_cons > 0) {
        nCons = abc->n_cons;
        if(nCons > 2) nCons = 2;

        // Backbone MSA row is always at msa_base[0].
        const vector<uint8_t> bbMsaRow(
            abc->msa_base[0], abc->msa_base[0] + abc->msa_len);
        const uint32_t bbStartPos = info.backboneStartPos;

        if(nCons == 2) {
            // Two clusters from de-novo clustering (max_n_cons=2).
            // The backbone (abPOA index 0) is assigned to one cluster.
            // That cluster becomes hap1 (results[0]), the other hap2 (results[1]).

            // Find which cluster contains the backbone.
            int bbCluster = -1;
            for(int ci = 0; ci < 2 && bbCluster < 0; ci++) {
                for(int k = 0; k < abc->clu_n_seq[ci]; k++) {
                    if(abc->clu_read_ids[ci][k] == 0) {
                        bbCluster = ci;
                        break;
                    }
                }
            }
            if(bbCluster < 0) bbCluster = 0; // fallback: assume cluster 0

            // Map abPOA cluster index → output index.
            // bbCluster → results[0] (hap1), other → results[1] (hap2).
            const int outMap[2] = {
                (bbCluster == 0) ? 0 : 1,
                (bbCluster == 0) ? 1 : 0
            };

            for(int ci = 0; ci < 2; ci++) {
                const int oi = outMap[ci]; // output index
                auto& r = resultsOut[oi];
                r.msaLen = abc->msa_len;
                r.backboneMsaRow = bbMsaRow;
                r.backboneStartPos = bbStartPos;
                r.consensusLen = abc->cons_len[ci];
                r.consensusSeq.assign(
                    abc->cons_base[ci],
                    abc->cons_base[ci] + abc->cons_len[ci]);

                // Consensus MSA row at abc->n_seq + ci.
                r.consensusMsaRow.assign(
                    abc->msa_base[abc->n_seq + ci],
                    abc->msa_base[abc->n_seq + ci] + abc->msa_len);

                // Collect sequences in this cluster.
                // Backbone (abPOA index 0) is only included in the hap1
                // cluster (oi == 0) since it's a hap1 read.
                // Map info indices → overlap indices for downstream use.
                for(int k = 0; k < abc->clu_n_seq[ci]; k++) {
                    int abpoaIdx = abc->clu_read_ids[ci][k];
                    if(abpoaIdx < 0 || abpoaIdx > nReads) continue;

                    if(abpoaIdx == 0) {
                        // Backbone — only include in hap1 (oi == 0).
                        if(oi == 0 && includeBackboneInReads) {
                            r.readIndices.push_back(-1);
                            r.readMsaRows.emplace_back(
                                abc->msa_base[0],
                                abc->msa_base[0] + abc->msa_len);
                        }
                    } else {
                        int infoIdx = readIndices[abpoaIdx - 1];
                        r.readIndices.push_back(info.overlapIndices[infoIdx]);
                        r.readMsaRows.emplace_back(
                            abc->msa_base[abpoaIdx],
                            abc->msa_base[abpoaIdx] + abc->msa_len);
                    }
                }
                r.nReads = int(r.readIndices.size());
            }
        } else {
            // Single consensus (per-haplotype path, max_n_cons=1).
            auto& r = resultsOut[0];
            r.msaLen = abc->msa_len;
            r.backboneMsaRow = bbMsaRow;
            r.backboneStartPos = bbStartPos;
            r.consensusLen = abc->cons_len[0];
            r.consensusSeq.assign(
                abc->cons_base[0],
                abc->cons_base[0] + abc->cons_len[0]);

            // Consensus MSA row at abc->n_seq.
            r.consensusMsaRow.assign(
                abc->msa_base[abc->n_seq],
                abc->msa_base[abc->n_seq] + abc->msa_len);

            // Include backbone (if requested) + all haplotype reads.
            const int extraRows = includeBackboneInReads ? 1 : 0;
            r.readIndices.reserve(extraRows + nReads);
            r.readMsaRows.reserve(extraRows + nReads);

            if(includeBackboneInReads) {
                r.readIndices.push_back(-1);
                r.readMsaRows.push_back(bbMsaRow);
            }

            // Haplotype reads — map info indices → overlap indices.
            for(int k = 0; k < nReads; k++) {
                int infoIdx = readIndices[k];
                r.readIndices.push_back(info.overlapIndices[infoIdx]);
                r.readMsaRows.emplace_back(
                    abc->msa_base[k + 1],
                    abc->msa_base[k + 1] + abc->msa_len);
            }
            r.nReads = int(r.readIndices.size());
        }
    }

    abpoa_free_para(abpt);
    abpoa_free(ab);
    return nCons;
}

// ============================================================================
// collectNoisyRegMsa — top-level
// ============================================================================

int dinara::collectNoisyRegMsa(
    const Assembler& assembler,
    const KmScratchpad& scratch,
    uint32_t regStart, uint32_t regEnd,
    const KmNoisyMsaOptions& opts,
    array<KmAbpoaMsaResult, 2>& results)
{
    results = {};

    const uint32_t regLen = regEnd - regStart;
    if(regLen == 0 || int(regLen) > opts.maxNoisyRegLen) {
        return 0;
    }

    KmNoisyReadInfo info = collectNoisyReadInfo(
        assembler, scratch, regStart, regEnd);

    if(info.nReads == 0 || info.nReads > opts.maxNoisyRegCov) {
        return 0;
    }
    if(info.backboneLen == 0) {
        return 0;
    }

    if(hasBothHaplotypes(info, opts.minHapReads)) {
        // Per-haplotype path: separate abPOA per hap, max_n_cons=1.
        // Backbone is included in hap1 reads (h=0) but not hap2 (h=1).
        int nCons = 0;
        for(int h = 0; h < 2; h++) {
            int hap = h + 1;
            vector<int> hapIndices;
            for(int i = 0; i < info.nReads; i++) {
                if(info.haps[i] == hap) hapIndices.push_back(i);
            }
            if(int(hapIndices.size()) < opts.minHapReads) continue;

            const bool includeBb = (h == 0); // backbone is hap1
            array<KmAbpoaMsaResult, 2> hapResults;
            int ret = abpoaMsaRun(info, hapIndices, 1, opts, includeBb, hapResults);
            if(ret > 0) {
                results[h] = move(hapResults[0]);
                nCons++;
            }
        }
        // Like pgphase with_ps_hap: require both haplotypes to succeed.
        if(nCons != 2) return 0;
        return 2;
    } else {
        // Combined fallback: all reads in one POA, max_n_cons=2.
        // Like pgphase no_ps_hap path.
        // Backbone is included in the hap1 cluster (handled inside abpoaMsaRun).
        vector<int> allIndices;
        allIndices.reserve(info.nReads);
        for(int i = 0; i < info.nReads; i++) {
            allIndices.push_back(i);
        }
        return abpoaMsaRun(info, allIndices, 2, opts, true, results);
    }
}
