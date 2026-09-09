/// @file AssemblerPhasing.cpp
/// @brief Hifiasm-parity ONT overlap phasing using OverlapCigarStore.
///
/// Ports hifiasm's rphase_hc pipeline (Correct.cpp:20191) to Dinara.
/// Reads directly from OverlapCigarStore — no AlignedEvidenceStore needed.
///
/// Pipeline per query read:
///   1. Gather overlaps from alignment table
///   2. Sliding-window SNP detection (428 bp windows)
///   3. SNP matrix construction (filter + confirm sites)
///   3b. Adjacent site filter
///   4. DP phasing (longest compatible chain)
///   5. Allele grouping (greedy cis/trans labeling + consistency check)
///   5b. multi_check (weak site promotion)
///   6. Large indel phasing (>=16 bp SV detection + BFS clustering)
///   7. Dedup chains (best overlap per target read)
///
/// @reference hifiasm v0.25.0-r726 (ec9a8b2)

#include "Assembler.hpp"
#include "PhasingTypes.hpp"
#include "Alignment.hpp"
#include "Reads.hpp"
#include "invalid.hpp"
#include "radixSort.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <thread>
#include <vector>

using namespace std;
using namespace dinara;

// ============================================================================
// cigarRead1Start: convert ad.ts (forward) to oriented coordinates for the
// CIGAR cursor. ad.qs/ts now store marker-based (non-extended) positions in
// forward coordinates. The CIGAR's yk is in oriented coordinates, so for RC
// overlaps we convert: oriented_start = targetLen - ad.te.
// ============================================================================

static uint64_t cigarRead1Start(
    const Assembler& assembler,
    const AlignmentData& ad)
{
    if (ad.isSameStrand) {
        return ad.ts;
    }
    const uint32_t targetLen = uint32_t(
        assembler.getReads().getRead(ad.readIds[1]).baseCount);
    return targetLen - ad.te;
}

// ============================================================================
// gatherOverlaps: collect all overlaps for a query read
// ============================================================================

static void gatherOverlaps(
    const Assembler& assembler,
    ReadId queryReadId,
    PhasingScratchpad& scratch)
{
    const auto& alignmentTable = assembler.getAlignmentTable();
    const auto& alignmentData = assembler.alignmentData;

    // alignmentTable is indexed by OrientedReadId::getValue().
    // Use strand 0 (forward) for the query.
    const OrientedReadId orientedQueryId(queryReadId, 0);
    if (orientedQueryId.getValue() >= alignmentTable.size()) return;
    const auto alignmentIds = alignmentTable[orientedQueryId.getValue()];

    for (const uint32_t alignmentId : alignmentIds) {
        const AlignmentData& ad = alignmentData[alignmentId];
        const AlignmentInfo& info = ad.info;

        // Skip deleted overlaps.
        if (ad.isDeleted()) continue;

        // Skip overlaps without CIGAR data.
        if (info.cigarOffset == uint32_t(-1)) continue;
        if (info.cigarTokenCount == 0) continue;

        PhasingOverlap ov;
        ov.alignmentId = alignmentId;
        ov.cigarOffset = info.cigarOffset;
        ov.cigarTokenCount = info.cigarTokenCount;
        ov.errorCount = info.mismatchCount + info.gapCount;
        ov.isMatch = 1; // default cis, matching hifiasm post-alignment state
        ov.strong = 0;

        // Determine which read in the pair is the query.
        // AlignmentData stores readIds[0] < readIds[1] (canonical order).
        // The CIGAR is always read0=query, read1=target in the store.
        if (ad.readIds[0] == queryReadId) {
            ov.queryIsRead0 = 1;
            ov.targetReadId = ad.readIds[1];
            ov.qs = ad.qs;
            ov.qe = ad.qe;
            ov.ts = ad.ts;
            ov.te = ad.te;
            ov.isRev = ad.isSameStrand ? 0 : 1;
        } else {
            ov.queryIsRead0 = 0;
            ov.targetReadId = ad.readIds[0];
            // Swap: query is read1, target is read0.
            ov.qs = ad.ts;
            ov.qe = ad.te;
            ov.ts = ad.qs;
            ov.te = ad.qe;
            ov.isRev = ad.isSameStrand ? 0 : 1;
        }

        scratch.overlaps.push_back(ov);
    }
}

// ============================================================================
// unpackQuerySequence: convert 2-bit packed bases to uint8_t array
// ============================================================================

static void unpackQuerySequence(
    const Assembler& assembler,
    ReadId queryReadId,
    uint32_t queryLen,
    PhasingScratchpad& scratch)
{
    scratch.queryBases.resize(queryLen);
    const auto sequence = assembler.getReads().getRead(queryReadId);
    for (uint32_t i = 0; i < queryLen; i++) {
        scratch.queryBases[i] = sequence[i].value;
    }
}

// ============================================================================
// getBaseAtPosition: look up a base from a read sequence
// ============================================================================

static inline uint8_t getBaseAtPosition(
    const Assembler& assembler,
    ReadId readId,
    uint32_t position,
    bool isReverseComplement)
{
    const auto sequence = assembler.getReads().getRead(readId);
    if (!isReverseComplement) {
        return sequence[position].value;
    } else {
        return sequence[sequence.baseCount - 1 - position].complement().value;
    }
}

// ============================================================================
// ============================================================================
// isPeriodicRepeat: hifiasm hpc_mask_ff equivalent
//
// Checks if position p on the query read is within a periodic repeat
// region of period 1..HPC_RR (4). Searches within a flanking window
// of ±HPC_PL (12) bases around p. For each period r, extends the
// repeat in both directions from p. If the repeat span exceeds
// r * HPC_CC (cutoff), the position is masked.
// ============================================================================

static constexpr uint32_t HPC_RR = 4;   // max repeat period
static constexpr uint32_t HPC_CC = 2;   // cutoff multiplier
static constexpr uint32_t HPC_PL = 12;  // flanking window size

static bool isPeriodicRepeat(
    const uint8_t* seq,
    uint32_t seqLen,
    uint32_t p)
{
    if (seqLen == 0) return false;

    const int64_t sn = int64_t(seqLen);
    const int64_t pp = int64_t(p);

    // Flanking window [s, e) around position p, clamped to [0, seqLen).
    const int64_t s = (pp >= int64_t(HPC_PL)) ? (pp - int64_t(HPC_PL)) : 0;
    const int64_t e = ((pp + int64_t(HPC_PL)) <= sn) ? (pp + int64_t(HPC_PL)) : sn;

    for (uint32_t r = 1; r <= HPC_RR; r++) {
        const int64_t rc = int64_t(r) * int64_t(HPC_CC); // cutoff

        // Scan 0: including p, extend right then left.
        {
            int64_t k = pp + int64_t(r);
            while (k < e && (k - int64_t(r)) >= s && seq[k] == seq[k - r]) k++;
            int64_t ze = k; if (ze > e) ze = e;
            k = pp - 1;
            while (k >= s && (k + int64_t(r)) < e && seq[k] == seq[k + r]) k--;
            int64_t zs = k + 1; if (zs < s) zs = s;
            if ((ze - zs) > int64_t(r) && (ze - zs) >= rc) return true;
        }

        // Scan 1: not including p, right side only.
        {
            int64_t k = pp + int64_t(r) + 1;
            while (k < e && (k - int64_t(r)) >= s && seq[k] == seq[k - r]) k++;
            int64_t zs = pp + 1; if (zs < s) zs = s;
            int64_t ze = k; if (ze > e) ze = e;
            if ((ze - zs) > int64_t(r) && (ze - zs) >= rc) return true;
        }

        // Scan 2: including p, extend left then right.
        {
            int64_t k = pp - int64_t(r);
            while (k >= s && (k + int64_t(r)) < e && seq[k] == seq[k + r]) k--;
            int64_t zs = k + 1; if (zs < s) zs = s;
            k = pp + 1;
            while (k < e && (k - int64_t(r)) >= s && seq[k] == seq[k - r]) k++;
            int64_t ze = k; if (ze > e) ze = e;
            if ((ze - zs) > int64_t(r) && (ze - zs) >= rc) return true;
        }

        // Scan 3: not including p, left side only.
        {
            int64_t k = pp - int64_t(r) - 1;
            while (k >= s && (k + int64_t(r)) < e && seq[k] == seq[k + r]) k--;
            int64_t zs = k + 1; if (zs < s) zs = s;
            int64_t ze = pp; if (ze > e) ze = e;
            if ((ze - zs) > int64_t(r) && (ze - zs) >= rc) return true;
        }
    }

    return false;
}

// ============================================================================
// detectSnpSites: single-pass SNP detection
//
// Fused single-pass approach (hifiasm extract_sub_cigar_hc pattern):
//
// 1. Walk each overlap's CIGAR once. At mismatch positions, increment
//    vote counter AND record (qpos, tpos) per overlap.
// 2. Classify: positions with > OCC_THRES votes become candidates.
// 3. Emit evidence from recorded data — no second CIGAR walk needed.
//    - Alt evidence: from recorded mismatches at candidate positions.
//    - Match evidence: for each candidate position, any overlap whose
//      query range covers it and has no mismatch/indel there.
//
// To handle indels correctly, we also record per-overlap indel ranges
// during the walk. A candidate position inside an indel range for an
// overlap does NOT get match evidence from that overlap.
// ============================================================================

static void detectSnpSites(
    const Assembler& assembler,
    ReadId queryReadId,
    uint32_t queryLen,
    PhasingScratchpad& scratch)
{
    const auto& cigarStore = assembler.getOverlapCigarStore();
    const uint32_t numOv = uint32_t(scratch.overlaps.size());

    // ---- Single-pass CIGAR walk: vote counting + mismatch recording ----
    scratch.flag.assign(queryLen, 0);

    // Per-overlap mismatch records: (qpos, tpos) pairs.
    // Flat vector with per-overlap range indices.
    scratch.mismatchRecords.clear();
    scratch.mismatchRangeBegin.resize(numOv + 1);

    // Per-overlap indel coverage: positions where this overlap has an indel.
    // We track indel ranges as (start, end) pairs per overlap.
    // For efficiency, use a flat vector with per-overlap range indices.
    scratch.indelRecords.clear();
    scratch.indelRangeBegin.resize(numOv + 1);

    for (size_t oi = 0; oi < numOv; oi++) {
        scratch.mismatchRangeBegin[oi] = uint32_t(scratch.mismatchRecords.size());
        scratch.indelRangeBegin[oi] = uint32_t(scratch.indelRecords.size());
        const auto& ov = scratch.overlaps[oi];
        const auto& ad = assembler.alignmentData[ov.alignmentId];
        const bool needsRcConvert =
            (ov.queryIsRead0 == 0) && (ov.isRev != 0);

        cigarStore.forEachOpWithPositions(
            ov.cigarOffset, ov.cigarTokenCount,
            ad.qs, cigarRead1Start(assembler, ad),
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (op == 1) {
                    // Mismatch: increment votes and record positions.
                    if (!needsRcConvert) {
                        const uint32_t qStart_op = ov.queryIsRead0
                            ? uint32_t(xk) : uint32_t(yk);
                        const uint32_t qEnd_op = qStart_op + len;
                        for (uint32_t qpos = qStart_op; qpos < qEnd_op; qpos++) {
                            if (qpos >= queryLen) break;
                            if (scratch.flag[qpos] < 255)
                                scratch.flag[qpos]++;
                            const uint32_t tpos = ov.queryIsRead0
                                ? uint32_t(yk) + (qpos - uint32_t(xk))
                                : uint32_t(xk) + (qpos - uint32_t(yk));
                            scratch.mismatchRecords.push_back({qpos, tpos});
                        }
                    } else {
                        for (uint32_t b = 0; b < len; b++) {
                            const uint32_t qpos_rc = uint32_t(yk) + b;
                            if (qpos_rc >= queryLen) continue;
                            const uint32_t qpos = queryLen - 1 - qpos_rc;
                            if (scratch.flag[qpos] < 255)
                                scratch.flag[qpos]++;
                            const uint32_t tpos = uint32_t(xk) + b;
                            scratch.mismatchRecords.push_back({qpos, tpos});
                        }
                    }
                } else if (op == CigarOpIns || op == CigarOpDel) {
                    // Record query positions consumed by indels. The query read
                    // (xk if queryIsRead0, else yk) is consumed by this indel
                    // iff the op consumes that read's sequence.
                    const bool queryConsumed = ov.queryIsRead0
                        ? opConsumesQuery(op) : opConsumesTarget(op);
                    if (queryConsumed) {
                        // Query positions in CIGAR coordinates.
                        const uint32_t rawStart = ov.queryIsRead0
                            ? uint32_t(xk) : uint32_t(yk);
                        const uint32_t rawEnd = std::min(rawStart + len, queryLen);
                        if (rawStart < rawEnd) {
                            if (!needsRcConvert) {
                                scratch.indelRecords.push_back({rawStart, rawEnd});
                            } else {
                                // RC: reverse the range.
                                const uint32_t fwdStart = queryLen - rawEnd;
                                const uint32_t fwdEnd = queryLen - rawStart;
                                scratch.indelRecords.push_back({fwdStart, fwdEnd});
                            }
                        }
                    }
                }
            });
    }
    scratch.mismatchRangeBegin[numOv] = uint32_t(scratch.mismatchRecords.size());
    scratch.indelRangeBegin[numOv] = uint32_t(scratch.indelRecords.size());

    // Classify: candidate SNPs where vote count > OCC_THRES.
    bool anyCandidates = false;
    for (uint32_t i = 0; i < queryLen; i++) {
        if (scratch.flag[i] > PHASING_OCC_THRES) {
            if (isPeriodicRepeat(scratch.queryBases.data(), queryLen, i)) {
                scratch.flag[i] = 3; // HPC-masked candidate
            } else {
                scratch.flag[i] = 1; // candidate SNP
            }
            anyCandidates = true;
        } else {
            scratch.flag[i] = 0;
        }
    }
    if (!anyCandidates) return;

    // Build sorted candidate position list.
    scratch.candidatePositions.clear();
    for (uint32_t i = 0; i < queryLen; i++) {
        if (scratch.flag[i] != 0) {
            scratch.candidatePositions.push_back(i);
        }
    }
    const uint32_t numCand = uint32_t(scratch.candidatePositions.size());

    // ---- Emit alt evidence from recorded mismatches ----
    for (uint32_t oi = 0; oi < numOv; oi++) {
        const auto& ov = scratch.overlaps[oi];
        const bool targetIsRc =
            (ov.queryIsRead0 != 0) && (ov.isRev != 0);
        const bool needsRcConvert =
            (ov.queryIsRead0 == 0) && (ov.isRev != 0);

        const uint32_t mBegin = scratch.mismatchRangeBegin[oi];
        const uint32_t mEnd = scratch.mismatchRangeBegin[oi + 1];
        for (uint32_t mi = mBegin; mi < mEnd; mi++) {
            const auto& mr = scratch.mismatchRecords[mi];
            if (scratch.flag[mr.qpos] == 0) continue; // not a candidate

            PhasingEvidence ev;
            ev.site = mr.qpos;
            ev.overlapIdx = oi;
            ev.siteIdx = UINT32_MAX;
            ev.isHpc = (scratch.flag[mr.qpos] == 3) ? 1 : 0;
            // For RC path (!queryIsRead0 && isRev), tpos is xk+b with targetIsRc=false.
            // For forward path with queryIsRead0 && isRev, targetIsRc=true.
            ev.base = getBaseAtPosition(
                assembler, ReadId(ov.targetReadId),
                mr.tpos, needsRcConvert ? false : targetIsRc);
            ev.isAlt = 1;
            scratch.evidence.push_back(ev);
        }
    }

    // ---- Emit match evidence from range checks ----
    // Build per-overlap mismatch bitmap for fast lookup.
    scratch.mismatchBitmap.assign(size_t(numOv) * numCand, 0);
    for (uint32_t oi = 0; oi < numOv; oi++) {
        const uint32_t mBegin = scratch.mismatchRangeBegin[oi];
        const uint32_t mEnd = scratch.mismatchRangeBegin[oi + 1];
        for (uint32_t mi = mBegin; mi < mEnd; mi++) {
            const uint32_t qpos = scratch.mismatchRecords[mi].qpos;
            if (scratch.flag[qpos] == 0) continue;
            const auto it = std::lower_bound(
                scratch.candidatePositions.begin(),
                scratch.candidatePositions.end(), qpos);
            if (it != scratch.candidatePositions.end() && *it == qpos) {
                const uint32_t ci = uint32_t(it - scratch.candidatePositions.begin());
                scratch.mismatchBitmap[oi * numCand + ci] = 1;
            }
        }
    }

    // Build per-overlap indel bitmap: mark candidate positions inside indels.
    // indelBitmap[oi * numCand + ci] = 1 if candidate ci is inside an indel for overlap oi.
    scratch.indelBitmap.assign(size_t(numOv) * numCand, 0);
    for (uint32_t oi = 0; oi < numOv; oi++) {
        const uint32_t iBegin = scratch.indelRangeBegin[oi];
        const uint32_t iEnd = scratch.indelRangeBegin[oi + 1];
        for (uint32_t ii = iBegin; ii < iEnd; ii++) {
            const auto& ir = scratch.indelRecords[ii];
            // Find candidates within [ir.qpos, ir.tpos) (tpos used as end).
            auto lo = std::lower_bound(
                scratch.candidatePositions.begin(),
                scratch.candidatePositions.end(), ir.qpos);
            auto hi = std::lower_bound(lo,
                scratch.candidatePositions.end(), ir.tpos);
            for (auto it = lo; it != hi; ++it) {
                const uint32_t ci = uint32_t(it - scratch.candidatePositions.begin());
                scratch.indelBitmap[oi * numCand + ci] = 1;
            }
        }
    }

    for (uint32_t ci = 0; ci < numCand; ci++) {
        const uint32_t qpos = scratch.candidatePositions[ci];
        const uint8_t isHpc = (scratch.flag[qpos] == 3) ? 1 : 0;
        const uint8_t base = scratch.queryBases[qpos];

        for (uint32_t oi = 0; oi < numOv; oi++) {
            const auto& ov = scratch.overlaps[oi];
            // Check if overlap covers this position.
            if (qpos < ov.qs || qpos >= ov.qe) continue;
            // Skip if mismatch or indel at this position.
            if (scratch.mismatchBitmap[oi * numCand + ci]) continue;
            if (scratch.indelBitmap[oi * numCand + ci]) continue;

            PhasingEvidence ev;
            ev.site = qpos;
            ev.overlapIdx = oi;
            ev.siteIdx = UINT32_MAX;
            ev.isHpc = isHpc;
            ev.base = base;
            ev.isAlt = 0;
            scratch.evidence.push_back(ev);
        }
    }
}

// ============================================================================
// buildSnpMatrix: sort evidence, group by site, filter, confirm sites
//
// Port of hifiasm's radix_sort_haplotype_evdience_srt (by site) followed
// by per-site push_info (Correct.cpp:10511) which does radix sort by
// overlapID + fused dedup/count/filter/emit in one pass.
//
// 1. Counting sort all evidence by site position — O(n + maxSite).
// 2. Per site group (hifiasm push_info pattern):
//    a. Counting sort by overlapIdx — O(m + numOverlaps).
//    b. Single fused pass: dedup by overlapIdx, count match/alt/strand.
//    c. Filter: require matchCount >= S_HAP_COV, altCount >= INFOR_COV.
//    d. Emit qualifying evidence with siteIdx assigned, compacting in-place.
// 3. Evidence array contains only qualifying entries on exit.
// ============================================================================

static void buildSnpMatrix(
    PhasingScratchpad& scratch)
{
    if (scratch.evidence.empty()) return;

    const size_t n = scratch.evidence.size();
    const uint32_t numOverlaps = uint32_t(scratch.overlaps.size());

    // --- Counting sort by site position ---
    uint32_t maxSite = 0;
    for (size_t i = 0; i < n; i++) {
        if (scratch.evidence[i].site > maxSite)
            maxSite = scratch.evidence[i].site;
    }

    scratch.countBuf.assign(maxSite + 2, 0);
    for (size_t i = 0; i < n; i++) {
        scratch.countBuf[scratch.evidence[i].site + 1]++;
    }
    for (uint32_t i = 1; i <= maxSite + 1; i++) {
        scratch.countBuf[i] += scratch.countBuf[i - 1];
    }

    scratch.evidenceTmp.resize(n);
    for (size_t i = 0; i < n; i++) {
        scratch.evidenceTmp[scratch.countBuf[scratch.evidence[i].site]++] =
            scratch.evidence[i];
    }
    scratch.evidence.swap(scratch.evidenceTmp);
    // evidence is now sorted by site. evidenceTmp is free for reuse.

    // --- Per-site: fused sort-by-overlapIdx / dedup / count / filter / emit ---
    // Write qualifying evidence compactly into evidenceTmp (output buffer).
    // This mirrors hifiasm's push_info called per site group from rphase_hc.
    scratch.evidenceTmp.clear();
    scratch.evidenceTmp.reserve(n); // output ≤ input; avoid reallocations

    size_t i = 0;
    while (i < n) {
        const uint32_t currentSite = scratch.evidence[i].site;
        const size_t groupBegin = i;
        while (i < n && scratch.evidence[i].site == currentSite) i++;
        const size_t groupEnd = i;
        const size_t groupLen = groupEnd - groupBegin;

        // (a) Sort this site's evidence by overlapIdx.
        //     Hifiasm uses radix sort which falls back to insertion sort
        //     for groups ≤ 64 (RS_MIN_SIZE). In practice all groups are
        //     ≤ 64 (avg ~35), so insertion sort is the fast path.
        if (groupLen > 1) {
            auto* base = scratch.evidence.data() + groupBegin;
            if (groupLen <= 64) {
                // Insertion sort — matches hifiasm's rs_insertsort fallback.
                for (size_t j = 1; j < groupLen; j++) {
                    if (base[j].overlapIdx < base[j - 1].overlapIdx) {
                        PhasingEvidence tmp = base[j];
                        size_t k = j;
                        do {
                            base[k] = base[k - 1];
                            k--;
                        } while (k > 0 && tmp.overlapIdx < base[k - 1].overlapIdx);
                        base[k] = tmp;
                    }
                }
            } else {
                // Counting sort for rare large groups.
                scratch.countBuf.assign(numOverlaps + 1, 0);
                for (size_t j = groupBegin; j < groupEnd; j++) {
                    scratch.countBuf[scratch.evidence[j].overlapIdx + 1]++;
                }
                for (uint32_t k = 1; k <= numOverlaps; k++) {
                    scratch.countBuf[k] += scratch.countBuf[k - 1];
                }
                scratch.sortTmp.resize(groupLen);
                for (size_t j = groupBegin; j < groupEnd; j++) {
                    scratch.sortTmp[
                        scratch.countBuf[scratch.evidence[j].overlapIdx]++] =
                        scratch.evidence[j];
                }
                for (size_t j = 0; j < groupLen; j++) {
                    base[j] = scratch.sortTmp[j];
                }
            }
        }

        // (b) Fused dedup + count pass (hifiasm push_info lines 10517-10533).
        //     Walk sorted-by-overlapIdx evidence. For each overlapIdx run,
        //     keep first entry, accumulate match/alt/strand counts.
        uint32_t matchCount = 0;
        uint32_t fwdStrandRefCount = 0;
        uint32_t altCount[4] = {0, 0, 0, 0};
        uint8_t isHpc = 0;
        const uint8_t queryBase = scratch.queryBases[currentSite];

        // Dedup in-place within evidence[groupBegin..groupEnd).
        size_t dedupEnd = groupBegin;
        size_t m = groupBegin; // start of current overlapIdx run
        for (size_t k = groupBegin + 1; k <= groupEnd; k++) {
            if (k == groupEnd ||
                scratch.evidence[k].overlapIdx !=
                scratch.evidence[m].overlapIdx) {
                // Keep first entry of this overlapIdx run.
                const auto& ev = scratch.evidence[m];
                scratch.evidence[dedupEnd] = ev;
                if (ev.isAlt == 0) {
                    matchCount++;
                    if (scratch.overlaps[ev.overlapIdx].isRev == 0)
                        fwdStrandRefCount++;
                } else {
                    if (ev.base < 4) altCount[ev.base]++;
                }
                if (ev.isHpc) isHpc = 1;
                dedupEnd++;
                m = k;
            }
        }

        // +1 for query read itself (hifiasm: p->occ_0 = 1 + occ_0).
        matchCount += 1;

        // (c) Filter: strand bias, min counts.
        // Reject if all ref-matching overlaps are forward-strand
        // (hifiasm: rev_n == occ_0, Correct.cpp:10567).
        if (fwdStrandRefCount == matchCount - 1) continue;

        // Identify qualifying alt bases and create PhasingSite entries.
        uint32_t baseSiteIdx[4] = {UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX};
        uint32_t lastSiteIdx = UINT32_MAX;

        for (uint8_t b = 0; b < 4; b++) {
            if (b == queryBase) continue;
            if (altCount[b] < 2) continue;
            if (matchCount < PHASING_S_HAP_COV) continue;
            if (altCount[b] < PHASING_INFOR_COV) continue;

            PhasingSite ps;
            ps.site = currentSite;
            ps.queryBase = queryBase;
            ps.altBase = b;
            ps.isHpc = isHpc;
            ps.dpChainId = -1;
            ps.matchCount = matchCount;
            ps.altCount = altCount[b];
            ps.fwdStrandCount = fwdStrandRefCount + 1;
            ps.labelMatchCount = matchCount;
            ps.labelFwdStrandCount = fwdStrandRefCount + 1;
            ps.transConfirmed = 0;
            ps.cisReset = 0;
            ps.promoted = 0;
            ps.evidenceBegin = 0;
            ps.evidenceEnd = 0;

            lastSiteIdx = uint32_t(scratch.sites.size());
            baseSiteIdx[b] = lastSiteIdx;
            scratch.sites.push_back(ps);
        }

        if (lastSiteIdx == UINT32_MAX) continue;

        // (d) Emit qualifying evidence with siteIdx assigned.
        //     Append to evidenceTmp (output buffer). Set evidence ranges.
        //     match → lastSiteIdx, mismatch → baseSiteIdx[base].
        const uint32_t sharedBegin = uint32_t(scratch.evidenceTmp.size());

        for (size_t j = groupBegin; j < dedupEnd; j++) {
            const auto& ev = scratch.evidence[j];
            if (ev.isAlt == 0) {
                PhasingEvidence out = ev;
                out.siteIdx = lastSiteIdx;
                scratch.evidenceTmp.push_back(out);
            } else if (ev.base < 4 && baseSiteIdx[ev.base] != UINT32_MAX) {
                PhasingEvidence out = ev;
                out.siteIdx = baseSiteIdx[ev.base];
                scratch.evidenceTmp.push_back(out);
            }
            // Non-qualifying mismatches → dropped
        }

        const uint32_t sharedEnd = uint32_t(scratch.evidenceTmp.size());

        for (uint8_t b = 0; b < 4; b++) {
            if (baseSiteIdx[b] != UINT32_MAX) {
                scratch.sites[baseSiteIdx[b]].evidenceBegin = sharedBegin;
                scratch.sites[baseSiteIdx[b]].evidenceEnd = sharedEnd;
            }
        }
    }

    // Output is in evidenceTmp. Swap into evidence.
    scratch.evidence.swap(scratch.evidenceTmp);
}

// ============================================================================
// filterAdjacentSites: reject DP-confirmed sites at adjacent positions
//
// Port of hifiasm's adjacent-site filter (Correct.cpp:8855).
//
// In hifiasm this runs inside generate_haplotypes_naive_HiFi, AFTER
// gen_rphase_dp has physically removed DP-rejected sites. It removes
// sites at adjacent positions (p and p+1) from the compacted arrays.
//
// In dinara, DP-rejected sites remain with dpChainId < 0. This function
// runs after runDpPhasing and marks adjacent DP-confirmed sites as
// rejected (dpChainId = -1). Only DP-confirmed sites (dpChainId >= 0)
// are considered when checking adjacency, matching hifiasm's
// post-compaction behavior. labelCisTrans guards on dpChainId in every
// step, so rejected sites are invisible to the labeling logic.
// ============================================================================

static void filterAdjacentSites(
    PhasingScratchpad& scratch)
{
    if (scratch.sites.size() < 2) return;

    const size_t n = scratch.sites.size();

    // Group sites by position. Only groups with at least one
    // DP-confirmed site participate in adjacency checks.
    struct PosGroup { uint32_t pos; size_t begin; size_t end; };
    vector<PosGroup> groups;
    size_t gi = 0;
    while (gi < n) {
        const uint32_t pos = scratch.sites[gi].site;
        size_t gEnd = gi;
        while (gEnd < n && scratch.sites[gEnd].site == pos) gEnd++;
        bool hasConfirmed = false;
        for (size_t i = gi; i < gEnd; i++) {
            if (scratch.sites[i].dpChainId >= 0) {
                hasConfirmed = true;
                break;
            }
        }
        if (hasConfirmed) {
            groups.push_back({pos, gi, gEnd});
        }
        gi = gEnd;
    }

    // Hifiasm (Correct.cpp:8855): if consecutive groups have positions
    // p and p+1, both are removed. We set dpChainId = -1 instead.
    for (size_t g = 0; g + 1 < groups.size(); g++) {
        if (groups[g + 1].pos == groups[g].pos + 1) {
            for (size_t i = groups[g].begin; i < groups[g].end; i++)
                scratch.sites[i].dpChainId = -1;
            for (size_t i = groups[g + 1].begin; i < groups[g + 1].end; i++)
                scratch.sites[i].dpChainId = -1;
        }
    }
}

// ============================================================================
// checkCompatibility: are two SNP sites consistent across shared overlaps?
//
// Port of hifiasm's comput_sc_rphase (Correct.cpp:9600).
//
// Two sites are compatible if overlaps covering both show consistent
// allele assignments: both match the query, or both mismatch.
// Mixed assignments (match at one, mismatch at the other) are skipped.
// Returns true only when evidence exists from both haplotypes (nn0>0, nn1>0).
// ============================================================================

static bool checkCompatibility(
    const PhasingScratchpad& scratch,
    uint32_t siteI,
    uint32_t siteJ)
{
    const auto& si = scratch.sites[siteI];
    const auto& sj = scratch.sites[siteJ];

    if (si.site == sj.site) return false;

    // Walk both evidence ranges (sorted by overlapIdx) in tandem.
    uint32_t pi = si.evidenceBegin;
    uint32_t pj = sj.evidenceBegin;

    // Hifiasm fi/fj classification per evidence entry:
    //   fi=0: match (isAlt==0)
    //   fi=1: mismatch for THIS PhasingSite (isAlt==1 && siteIdx==siteI)
    //   fi=2: mismatch for DIFFERENT PhasingSite at same position
    //
    // Counting rules (comput_sc_rphase):
    //   Both fi=2 with valid siteIdx → treat as fi=fj=0
    //   Either fi=2 (but not both) → incompatible, return false
    //   fi != fj → incompatible, return false
    //   fi == fj → nn[fi]++
    //
    // Compatible iff nn[0] > 0 && nn[1] > 0.
    uint32_t nn0 = 0, nn1 = 0;

    while (pi < si.evidenceEnd && pj < sj.evidenceEnd) {
        const auto& ei = scratch.evidence[pi];
        const auto& ej = scratch.evidence[pj];

        if (ei.overlapIdx < ej.overlapIdx) { pi++; continue; }
        if (ei.overlapIdx > ej.overlapIdx) { pj++; continue; }

        // Classify fi.
        uint8_t fi = 2;
        if (ei.isAlt == 0) {
            fi = 0;
        } else if (ei.siteIdx == siteI) {
            fi = 1;
        }

        // Classify fj.
        uint8_t fj = 2;
        if (ej.isAlt == 0) {
            fj = 0;
        } else if (ej.siteIdx == siteJ) {
            fj = 1;
        }

        // Both fi=2: mismatch for different sites at same position.
        // Treat as both matching if both have valid siteIdx.
        if (fi == 2 && fj == 2 &&
            ei.siteIdx != UINT32_MAX && ej.siteIdx != UINT32_MAX) {
            fi = fj = 0;
        }

        if (fi == 2 || fj == 2) return false;
        if (fi != fj) return false;
        if (fi == 0) nn0++;
        else nn1++;

        pi++;
        pj++;
    }

    return nn0 > 0 && nn1 > 0;
}

// ============================================================================
// isHpcDependent: would this site fail thresholds without HPC evidence?
//
// Port of hifiasm's is_hpc_vec (Correct.cpp:9383).
//
// Temporarily subtracts HPC evidence (isHpc==1) from match/alt counts
// and checks if the site still passes thresholds. Returns true if the
// site depends on HPC evidence to pass (i.e., should be rejected for
// single-site chains).
// ============================================================================

static bool isHpcDependent(
    const PhasingScratchpad& scratch,
    uint32_t siteIdx)
{
    const auto& site = scratch.sites[siteIdx];
    uint32_t mc = site.matchCount;
    uint32_t ac = site.altCount;

    for (uint32_t ei = site.evidenceBegin; ei < site.evidenceEnd; ei++) {
        const auto& ev = scratch.evidence[ei];
        if (!ev.isHpc) continue;
        if (ev.isAlt == 0) {
            if (mc > 0) mc--;
        } else if (ev.siteIdx == siteIdx) {
            if (ac > 0) ac--;
        }
    }

    if (mc < 2 || ac < 2) return true;
    if (mc < PHASING_S_HAP_COV || ac < PHASING_INFOR_COV) return true;
    return false;
}

// ============================================================================
// isStrandBiased: hifiasm is_st_bs macro equivalent
//
// Returns true if a site has strand bias — nearly all ref-matching
// overlaps come from one strand. Hifiasm parameters: st_rate=0.05, st_max=2.
// ============================================================================

static constexpr double PHASING_ST_RATE = 0.05;
static constexpr uint32_t PHASING_ST_MAX = 2;

static bool isStrandBiased(uint32_t matchCount, uint32_t fwdStrandCount)
{
    // fwdStrandCount = forward-strand ref overlaps + 1 (query)
    // matchCount = total ref overlaps + 1 (query)
    // Condition: reverse-strand count <= ST_MAX AND forward fraction >= (1 - ST_RATE)
    if (fwdStrandCount + PHASING_ST_MAX >= matchCount &&
        matchCount * PHASING_ST_RATE + fwdStrandCount >= matchCount) {
        return true;
    }
    return false;
}

/// Check strand bias using immutable counts (for DP and initial labeling).
static bool isStrandBiased(const PhasingSite& site)
{
    return isStrandBiased(site.matchCount, site.fwdStrandCount);
}

/// Check strand bias using mutable label counts (for greedy labeling).
static bool isLabelStrandBiased(const PhasingSite& site)
{
    return isStrandBiased(site.labelMatchCount, site.labelFwdStrandCount);
}

// ============================================================================
// runDpPhasing: longest compatible chain over confirmed SNP sites
//
// Port of hifiasm's gen_rphase_dp + gen_rphase_dp0_single_path
// (Correct.cpp:9648, 9428).
//
// Pre-DP filter: strand-biased sites (is_st_bs) are excluded from the
// DP entirely, matching hifiasm's gen_rphase_dp which physically removes
// them before calling gen_rphase_dp0_single_path.
//
// LIS-style DP: f[i] = max(f[j] + 1) for all j < i where sites i and j
// are compatible. Iterates j from i-1 down to 0 (hifiasm tie-breaking).
//
// Chain confirmation:
//   - Multi-site chains (length > 1): confirmed, but individual sites
//     with matchCount < cc are rejected (hifiasm per-site occ_0 check).
//   - Single-site chains: require !isHpcDependent AND matchCount >= cc.
//   cc = max(6, hetCov * 0.7).
// ============================================================================

static void runDpPhasing(
    PhasingScratchpad& scratch,
    uint32_t hetCov)
{
    const uint32_t n = uint32_t(scratch.sites.size());
    if (n == 0) return;

    const uint32_t cc = max(6U, uint32_t(double(hetCov) * 0.7));

    // Pre-DP filter: mark strand-biased sites for exclusion.
    // Hifiasm (gen_rphase_dp, Correct.cpp:9665) physically removes
    // strand-biased sites before running the DP. We mark them instead
    // and skip them during scoring and chain assignment.
    vector<bool> dpExcluded(n, false);
    for (uint32_t i = 0; i < n; i++) {
        if (isStrandBiased(scratch.sites[i])) {
            dpExcluded[i] = true;
        }
    }

    // Fill incomplete evidence (port of hifiasm fill_incom, Correct.cpp:9597).
    //
    // For each overlap, if it has evidence at site positions A and C but
    // not at intermediate site position B, add a synthetic match entry at
    // B. The assumption: if an overlap covers a site position and has no
    // evidence there, it's a match.
    //
    // Hifiasm does this in gen_rphase_dp (Correct.cpp:9749) before the DP.
    // After adding incomplete entries, it recounts occ_0/occ_1/overlap_num.
    {
        // Collect non-excluded site positions in order.
        vector<uint32_t> activeSitePositions;
        for (uint32_t i = 0; i < n; i++) {
            if (!dpExcluded[i]) {
                activeSitePositions.push_back(scratch.sites[i].site);
            }
        }
        // Deduplicate positions (multi-alt sites share positions).
        dinara::radixSort(activeSitePositions.data(),
            activeSitePositions.data() + activeSitePositions.size(),
            [](uint32_t x) { return x; });
        activeSitePositions.erase(
            unique(activeSitePositions.begin(), activeSitePositions.end()),
            activeSitePositions.end());

        if (activeSitePositions.size() > 1) {
            // Build position → index in activeSitePositions.
            unordered_map<uint32_t, uint32_t> posToIdx;
            for (uint32_t i = 0; i < uint32_t(activeSitePositions.size()); i++) {
                posToIdx[activeSitePositions[i]] = i;
            }

            const uint32_t numOv = uint32_t(scratch.overlaps.size());
            // Per-overlap: last seen position index in activeSitePositions.
            // UINT32_MAX = not yet seen.
            vector<uint32_t> lastPosIdx(numOv, UINT32_MAX);

            // Collect incomplete entries: (overlapIdx, sitePosition).
            vector<pair<uint32_t, uint32_t>> incompleteEntries;

            // Walk evidence sorted by site position (already sorted).
            for (const auto& ev : scratch.evidence) {
                auto it = posToIdx.find(ev.site);
                if (it == posToIdx.end()) continue; // position not active
                uint32_t curPosIdx = it->second;
                uint32_t oi = ev.overlapIdx;

                if (lastPosIdx[oi] == UINT32_MAX) {
                    // First time seeing this overlap — just record.
                    lastPosIdx[oi] = curPosIdx;
                } else if (curPosIdx > lastPosIdx[oi] + 1) {
                    // Gap: add match entries for intermediate positions.
                    for (uint32_t pi = lastPosIdx[oi] + 1; pi < curPosIdx; pi++) {
                        incompleteEntries.push_back({oi, activeSitePositions[pi]});
                    }
                    lastPosIdx[oi] = curPosIdx;
                } else {
                    lastPosIdx[oi] = curPosIdx;
                }
            }

            // Append synthetic match evidence entries.
            for (const auto& [oi, sitePos] : incompleteEntries) {
                PhasingEvidence ev;
                ev.site = sitePos;
                ev.overlapIdx = oi;
                ev.siteIdx = UINT32_MAX; // assigned below during recount
                ev.base = 0;
                ev.isAlt = 0; // match
                ev.isHpc = 0;
                scratch.evidence.push_back(ev);
            }

            if (!incompleteEntries.empty()) {
                // Re-sort evidence by (site, overlapIdx).
                dinara::radixSort(scratch.evidence.data(),
                    scratch.evidence.data() + scratch.evidence.size(),
                    [](const PhasingEvidence& e) -> uint64_t {
                        return (uint64_t(e.site) << 32) | e.overlapIdx;
                    });

                // Assign siteIdx to new entries and rebuild evidenceBegin/End.
                // For synthetic match entries (siteIdx == UINT32_MAX), assign
                // the lastSiteIdx at that position (same as buildSnpMatrix).
                uint32_t si = 0, ei = 0;
                while (si < n && ei < uint32_t(scratch.evidence.size())) {
                    uint32_t pos = scratch.sites[si].site;
                    uint32_t rangeBegin = ei;

                    // Find last site at this position (for lastSiteIdx).
                    uint32_t lastSi = si;
                    while (lastSi + 1 < n &&
                           scratch.sites[lastSi + 1].site == pos) {
                        lastSi++;
                    }

                    // Walk evidence at this position.
                    while (ei < uint32_t(scratch.evidence.size()) &&
                           scratch.evidence[ei].site == pos) {
                        if (scratch.evidence[ei].siteIdx == UINT32_MAX) {
                            // Synthetic match → assign to last site at position.
                            scratch.evidence[ei].siteIdx = lastSi;
                        }
                        ei++;
                    }

                    // Set evidenceBegin/End for all sites at this position.
                    for (uint32_t s = si; s <= lastSi; s++) {
                        scratch.sites[s].evidenceBegin = rangeBegin;
                        scratch.sites[s].evidenceEnd = ei;
                    }
                    si = lastSi + 1;
                }

                // Recount matchCount, fwdStrandCount, altCount from evidence.
                // Hifiasm recounts occ_0, occ_1, overlap_num after fill_incom.
                for (uint32_t i = 0; i < n; i++) {
                    auto& site = scratch.sites[i];
                    uint32_t mc = 0, fwd = 0, ac = 0;
                    for (uint32_t e = site.evidenceBegin; e < site.evidenceEnd; e++) {
                        const auto& ev = scratch.evidence[e];
                        if (ev.isAlt == 0) {
                            mc++;
                            if (scratch.overlaps[ev.overlapIdx].isRev == 0) {
                                fwd++;
                            }
                        } else if (ev.siteIdx == i) {
                            ac++;
                        }
                    }
                    site.matchCount = mc + 1;  // +1 for query
                    site.fwdStrandCount = fwd + 1; // +1 for query
                    site.altCount = ac;
                    // Also update mutable label copies.
                    site.labelMatchCount = site.matchCount;
                    site.labelFwdStrandCount = site.fwdStrandCount;
                }

                // Re-evaluate dpExcluded with updated counts.
                for (uint32_t i = 0; i < n; i++) {
                    dpExcluded[i] = isStrandBiased(scratch.sites[i]);
                }
            }
        }
    }

    scratch.dpScore.assign(n, 1);
    scratch.dpParent.assign(n, -1);

    // DP: iterate j from i-1 down to 0 (hifiasm tie-breaking: last
    // improving j wins, which is the smallest index due to reverse scan).
    for (uint32_t i = 1; i < n; i++) {
        if (dpExcluded[i]) continue;
        int32_t maxF = 1, maxJ = -1;
        for (int32_t j = int32_t(i) - 1; j >= 0; j--) {
            if (dpExcluded[uint32_t(j)]) continue;
            if (!checkCompatibility(scratch, i, uint32_t(j))) continue;
            int32_t sc = scratch.dpScore[uint32_t(j)] + 1;
            if (sc > maxF) {
                maxF = sc;
                maxJ = j;
            }
        }
        scratch.dpScore[i] = maxF;
        scratch.dpParent[i] = maxJ;
    }

    // Sort sites by score descending for greedy chain extraction.
    // Hifiasm: encode (UINT32_MAX - f[i]) << 32 | i, radix sort.
    scratch.dpEndpoints.clear();
    for (uint32_t i = 0; i < n; i++) {
        if (dpExcluded[i]) continue;
        scratch.dpEndpoints.push_back({scratch.dpScore[i], i});
    }
    sort(scratch.dpEndpoints.begin(), scratch.dpEndpoints.end(),
        [](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first > b.first;
            return a.second < b.second;
        });

    // Greedy chain extraction: start from highest-scoring site,
    // follow parent pointers, mark claimed sites.
    scratch.dpChainId.assign(n, -1);
    int32_t chainId = 0;

    // Temporary: collect (chainLength, chainStartInBuffer) pairs.
    vector<pair<uint32_t, uint32_t>> chains;
    vector<uint32_t> chainBuf;

    for (const auto& [score, startIdx] : scratch.dpEndpoints) {
        if (scratch.dpChainId[startIdx] >= 0) continue; // already claimed

        uint32_t bufStart = uint32_t(chainBuf.size());
        int32_t k = int32_t(startIdx);
        while (k >= 0 && scratch.dpChainId[uint32_t(k)] < 0) {
            chainBuf.push_back(uint32_t(k));
            scratch.dpChainId[uint32_t(k)] = INT32_MAX; // temporary claim
            k = scratch.dpParent[uint32_t(k)];
        }
        uint32_t chainLen = uint32_t(chainBuf.size()) - bufStart;
        if (chainLen > 0) {
            chains.push_back({chainLen, bufStart});
        }
    }

    // Reset temporary claims.
    scratch.dpChainId.assign(n, -1);

    // Sort chains by length descending, then by bufStart ascending
    // (hifiasm encodes ((UINT32_MAX - length) << 32) | rn0 and radix
    // sorts, so equal-length chains are ordered by extraction order,
    // i.e. higher-scoring endpoints first).
    sort(chains.begin(), chains.end(),
        [](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first > b.first;
            return a.second < b.second;
        });

    // Assign chain IDs with confirmation logic.
    // Hifiasm processes all chains without a limit.
    for (const auto& [chainLen, bufStart] : chains) {
        int8_t plus = -1; // default: not confirmed

        if (chainLen > 1) {
            plus = 1; // multi-site chain → confirmed
        } else {
            // Single-site chain: require !isHpcDependent AND matchCount >= cc.
            uint32_t idx = chainBuf[bufStart];
            if (!isHpcDependent(scratch, idx) &&
                scratch.sites[idx].matchCount >= cc) {
                plus = 1;
            }
        }

        // Assign chain ID to confirmed sites. For multi-site chains,
        // individual sites with matchCount < cc are still rejected
        // (hifiasm per-site occ_0 check).
        bool anyConfirmed = false;
        for (uint32_t i = 0; i < chainLen; i++) {
            uint32_t idx = chainBuf[bufStart + i];
            if (scratch.sites[idx].matchCount >= cc && plus > 0) {
                scratch.dpChainId[idx] = chainId;
                scratch.sites[idx].dpChainId = chainId;
                anyConfirmed = true;
            }
        }
        if (anyConfirmed) chainId++;
    }
}

// ============================================================================
// labelCisTrans: greedy cis/trans labeling from confirmed SNP sites
//
// Port of hifiasm's generate_haplotypes_naive_HiFi (Correct.cpp:8845).
//
// Six steps matching hifiasm exactly:
//   Step A — Count confirmed mismatches per overlap, build sorted list.
//   Step B — Greedy labeling (sorted by mismatch count desc): mark trans,
//            set transConfirmed for mismatch sites, decrement
//            labelMatchCount/labelFwdStrandCount for match evidence.
//   Step C — Consistency check: flip cis overlaps with mismatches at
//            transConfirmed sites.
//   Step D — Set cisReset for mismatch sites at cis overlaps.
//   Step E — multi_check: set promoted for weak sites in >=2 cis overlaps.
//   Step F — Final loop: set strong flag, flip cis→trans at confirmed
//            sites (transConfirmed && !cisReset, or promoted).
// ============================================================================

static void labelCisTrans(
    PhasingScratchpad& scratch)
{
    const uint32_t numSites = uint32_t(scratch.sites.size());
    const uint32_t numOv = uint32_t(scratch.overlaps.size());
    if (numSites == 0 || numOv == 0) return;

    // Build per-overlap evidence index. Hifiasm sorts evidence globally
    // by overlapID and iterates per-overlap. We build an index instead,
    // keeping our per-site evidence layout intact.
    //
    // overlapEvIdx: indices into scratch.evidence, sorted by overlapIdx.
    // overlapEvBegin[oi], overlapEvEnd[oi]: range for overlap oi.
    vector<uint32_t> overlapEvIdx(scratch.evidence.size());
    iota(overlapEvIdx.begin(), overlapEvIdx.end(), 0U);
    dinara::radixSort(overlapEvIdx.data(),
        overlapEvIdx.data() + overlapEvIdx.size(),
        [&](uint32_t a) { return scratch.evidence[a].overlapIdx; });

    vector<uint32_t> overlapEvBegin(numOv, 0);
    vector<uint32_t> overlapEvEnd(numOv, 0);
    {
        uint32_t pos = 0;
        for (uint32_t oi = 0; oi < numOv; oi++) {
            while (pos < overlapEvIdx.size() &&
                   scratch.evidence[overlapEvIdx[pos]].overlapIdx < oi) pos++;
            overlapEvBegin[oi] = pos;
            while (pos < overlapEvIdx.size() &&
                   scratch.evidence[overlapEvIdx[pos]].overlapIdx == oi) pos++;
            overlapEvEnd[oi] = pos;
        }
    }

    // Helper: count confirmed mismatches for one overlap using current
    // label counts. Mirrors hifiasm's inner loop in Step A/B.
    auto countConfirmedMismatches = [&](uint32_t oi) -> uint32_t {
        uint32_t o = 0;
        for (uint32_t p = overlapEvBegin[oi]; p < overlapEvEnd[oi]; p++) {
            const auto& ev = scratch.evidence[overlapEvIdx[p]];
            if (ev.isAlt == 0) continue; // skip matches
            if (ev.siteIdx >= numSites) continue;
            const auto& s = scratch.sites[ev.siteIdx];
            if (s.dpChainId < 0) continue;
            if (s.labelMatchCount < 2 || s.altCount < 2) continue;
            if (isLabelStrandBiased(s)) continue;
            if (s.labelMatchCount < PHASING_S_HAP_COV ||
                s.altCount < PHASING_INFOR_COV) continue;
            o++;
        }
        return o;
    };

    // ----------------------------------------------------------------
    // Step A: Count confirmed mismatches per overlap, build sorted list.
    //
    // Hifiasm (Correct.cpp:8889): sort evidence by overlapID, group by
    // overlap, count mismatches at confirmed SnpStats, push overlaps
    // with o>0 to snp_srt encoded as ((UINT32_MAX - o) << 32 | l).
    // ----------------------------------------------------------------

    for (auto& ov : scratch.overlaps) {
        ov.confirmedMismatchCount = countConfirmedMismatches(
            uint32_t(&ov - scratch.overlaps.data()));
    }

    // Build sorted list: overlaps with mismatches, sorted by count desc.
    // Encodes (count, overlapIdx) — hifiasm uses (UINT32_MAX - o, l).
    scratch.sortedOverlapIndices.clear();
    for (uint32_t oi = 0; oi < numOv; oi++) {
        if (scratch.overlaps[oi].confirmedMismatchCount > 0) {
            scratch.sortedOverlapIndices.push_back(oi);
        }
    }
    sort(scratch.sortedOverlapIndices.begin(),
         scratch.sortedOverlapIndices.end(),
         [&](uint32_t a, uint32_t b) {
             return scratch.overlaps[a].confirmedMismatchCount >
                    scratch.overlaps[b].confirmedMismatchCount;
         });

    // ----------------------------------------------------------------
    // Step B: Greedy labeling (most mismatches first).
    //
    // Hifiasm (Correct.cpp:8949): for each overlap in sorted order,
    // re-count mismatches (thresholds may have changed from decrements),
    // mark trans, set score=1 for mismatch SnpStats, decrement occ_0
    // on all SnpStats at same position for match evidence, decrement
    // overlap_num for forward-strand matches.
    // ----------------------------------------------------------------

    for (uint32_t sortIdx : scratch.sortedOverlapIndices) {
        // Re-count with current label counts (D1/D2 fix).
        uint32_t o = countConfirmedMismatches(sortIdx);
        if (o == 0) continue;

        auto& ov = scratch.overlaps[sortIdx];
        if (ov.isMatch == 1) ov.isMatch = 2; // trans

        // Walk this overlap's evidence: mark mismatch sites, decrement
        // match counts. Hifiasm physically removes DP-rejected sites
        // before this function, so they never appear. We skip them
        // via dpChainId < 0.
        for (uint32_t p = overlapEvBegin[sortIdx]; p < overlapEvEnd[sortIdx]; p++) {
            const auto& ev = scratch.evidence[overlapEvIdx[p]];
            if (ev.siteIdx >= numSites) continue;
            if (scratch.sites[ev.siteIdx].dpChainId < 0) continue;

            if (ev.isAlt == 1) {
                // Mismatch: set score=1 (Correct.cpp:8979).
                scratch.sites[ev.siteIdx].transConfirmed = 1;
            } else {
                // Match: decrement occ_0 on all sites at same position.
                // Hifiasm walks backwards only from overlapSite
                // (Correct.cpp:8983). Match evidence points to the last
                // site at this position, so backwards covers all.
                const uint32_t pos = scratch.sites[ev.siteIdx].site;
                for (int32_t si = int32_t(ev.siteIdx); si >= 0; si--) {
                    if (scratch.sites[uint32_t(si)].site != pos) break;
                    auto& s = scratch.sites[uint32_t(si)];
                    if (s.dpChainId < 0) continue;
                    if (s.labelMatchCount > 1) s.labelMatchCount--;
                    if (ov.isRev == 0 && s.labelFwdStrandCount > 1) {
                        s.labelFwdStrandCount--;
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------
    // Step C: Consistency check.
    //
    // Hifiasm (Correct.cpp:8998): second pass over sorted overlaps.
    // For each overlap, count mismatches at trans-confirmed sites.
    // If cis and count > 0, flip to trans.
    // (cisReset hasn't been set yet, so transConfirmed alone is correct.)
    // ----------------------------------------------------------------

    for (uint32_t sortIdx : scratch.sortedOverlapIndices) {
        uint32_t o = 0;
        for (uint32_t p = overlapEvBegin[sortIdx]; p < overlapEvEnd[sortIdx]; p++) {
            const auto& ev = scratch.evidence[overlapEvIdx[p]];
            if (ev.isAlt == 0) continue;
            if (ev.siteIdx >= numSites) continue;
            const auto& s = scratch.sites[ev.siteIdx];
            if (s.dpChainId < 0) continue;
            if (s.labelMatchCount < 2 || s.altCount < 2) continue;
            if (isLabelStrandBiased(s)) continue;
            if (s.transConfirmed) o++;
        }
        auto& ov = scratch.overlaps[sortIdx];
        if (ov.isMatch == 1 && o > 0) {
            ov.isMatch = 2;
        }
    }

    // ----------------------------------------------------------------
    // Step D: Mark cis-reset for mismatch sites at cis overlaps.
    //
    // Hifiasm (Correct.cpp:9015): for each cis overlap, set score=-1
    // on its mismatch SnpStats. This ensures only sites with mismatches
    // exclusively at trans overlaps remain trans-confirmed.
    // We set cisReset=1 instead of overwriting transConfirmed, keeping
    // both flags readable.
    // ----------------------------------------------------------------

    for (uint32_t oi = 0; oi < numOv; oi++) {
        if (scratch.overlaps[oi].isMatch != 1) continue;
        for (uint32_t p = overlapEvBegin[oi]; p < overlapEvEnd[oi]; p++) {
            const auto& ev = scratch.evidence[overlapEvIdx[p]];
            if (ev.isAlt == 1 && ev.siteIdx < numSites &&
                scratch.sites[ev.siteIdx].dpChainId >= 0) {
                // Hifiasm: score = -1 (Correct.cpp:9025).
                scratch.sites[ev.siteIdx].cisReset = 1;
            }
        }
    }

    // ----------------------------------------------------------------
    // Step E: multi_check — disabled for ONT (multi_check=0).
    // Hifiasm (Correct.cpp:9033): if(multi_check) { ... }
    // ONT path passes multi_check=0, so this entire block is skipped.

    // ----------------------------------------------------------------
    // Step F: Final loop (Correct.cpp:9085).
    //
    // For each overlap:
    //   trans → strong=1
    //   cis → check evidence at sites with score==1 && occ_0>=2 &&
    //          occ_1>=2 && !is_st_bs. If found: strong=1.
    //          If mismatch: also flip to trans and break.
    //
    // score==1 at this point means transConfirmed was set in Step B
    // and NOT overwritten to -1 in Step D. With multi_check disabled,
    // this is (transConfirmed && !cisReset).
    // ----------------------------------------------------------------

    for (uint32_t oi = 0; oi < numOv; oi++) {
        auto& ov = scratch.overlaps[oi];
        if (ov.isMatch == 2) {
            ov.strong = 1;
        } else if (ov.isMatch == 1) {
            for (uint32_t p = overlapEvBegin[oi]; p < overlapEvEnd[oi]; p++) {
                const auto& ev = scratch.evidence[overlapEvIdx[p]];
                if (ev.siteIdx >= numSites) continue;
                if (ev.isAlt != 0 && ev.isAlt != 1) continue;
                const auto& s = scratch.sites[ev.siteIdx];
                if (s.dpChainId < 0) continue;
                // score == 1: transConfirmed && !cisReset
                if (!(s.transConfirmed && !s.cisReset)) continue;
                if (s.labelMatchCount < 2 || s.altCount < 2) continue;
                if (isLabelStrandBiased(s)) continue;
                ov.strong = 1;
                if (ev.isAlt == 1) {
                    ov.isMatch = 2; // flip cis → trans
                    break;
                }
            }
        }
    }
}

// ============================================================================
// phaseLargeIndels: detect >=16 bp SVs from CIGARs, cluster, phase
//
// Port of hifiasm's rphase_lidel (Correct.cpp:20155).
//
// Only considers cis overlaps (is_match == 1). Four steps:
//   Step 1 — Walk CIGARs and detect contiguous indel regions >= 16 bp.
//   Step 2 — Sort SV events by query position.
//   Step 3 — Cluster events by >= 50% position overlap (BFS connected
//            components), with target-read dedup within each cluster.
//   Step 4 — For each cluster with >= 3 unique targets and >= 3 spanning
//            cis overlaps, label the SV-carrying overlaps as trans.
// ============================================================================

static void phaseLargeIndels(
    const Assembler& assembler,
    ReadId queryReadId,
    uint32_t queryLen,
    PhasingScratchpad& scratch)
{
    const auto& cigarStore = assembler.getOverlapCigarStore();

    // ----------------------------------------------------------------
    // Step 1: Detect SV events from cis overlaps.
    //
    // Walk each cis overlap's CIGAR. Contiguous indel runs (ins/del ops)
    // whose span or base count >= SV_MIN_LEN (16) are recorded as
    // PhasingSvEvents. RC overlaps have their coordinates converted to
    // forward-strand afterward.
    // ----------------------------------------------------------------

    for (uint32_t oi = 0; oi < uint32_t(scratch.overlaps.size()); oi++) {
        const auto& ov = scratch.overlaps[oi];
        if (ov.isMatch != 1) continue; // only cis overlaps

        const auto& ad = assembler.alignmentData[ov.alignmentId];

        // Walk the full CIGAR and detect contiguous error regions.
        uint32_t indelRunStart = UINT32_MAX;
        uint32_t indelRunEnd = 0;
        uint32_t indelRunBases = 0;

        auto flushRun = [&]() {
            if (indelRunStart == UINT32_MAX) return;
            uint32_t runLen = indelRunEnd - indelRunStart;
            if (runLen >= PHASING_SV_MIN_LEN || indelRunBases >= PHASING_SV_MIN_LEN) {
                PhasingSvEvent ev;
                ev.overlapIdx = oi;
                ev.queryPos = indelRunStart;
                ev.queryEnd = indelRunEnd;
                ev.errorBases = indelRunBases;
                scratch.svEvents.push_back(ev);
            }
            indelRunStart = UINT32_MAX;
            indelRunEnd = 0;
            indelRunBases = 0;
        };

        OverlapCigarStore::Cursor cursor;
        cursor.reset(ov.cigarOffset, ov.cigarTokenCount,
            ad.qs, cigarRead1Start(assembler, ad), cigarStore);

        // Track SV events in oriented coordinates during the walk,
        // then convert to forward coordinates afterward.
        const bool qNeedsRc =
            (ov.queryIsRead0 == 0) && (ov.isRev != 0);
        const size_t svEventsBefore = scratch.svEvents.size();

        cigarStore.walkRangeWithCursor(
            cursor, 0, UINT32_MAX,
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                // Map to query position (oriented coordinates).
                uint32_t qpos = ov.queryIsRead0 ? uint32_t(xk) : uint32_t(yk);

                if (op == CigarOpIns || op == CigarOpDel) { // insertion or deletion
                    if (indelRunStart == UINT32_MAX) {
                        indelRunStart = qpos;
                    }
                    // How many query bases this op spans: the query read
                    // (xk if queryIsRead0, else yk) advances iff the op consumes
                    // that read's sequence.
                    const bool queryAdvances = ov.queryIsRead0
                        ? opConsumesQuery(op) : opConsumesTarget(op);
                    indelRunEnd = qpos + (queryAdvances ? len : 0);
                    if (indelRunEnd <= indelRunStart) {
                        indelRunEnd = indelRunStart + 1;
                    }
                    indelRunBases += len;
                } else {
                    // Match or mismatch — flush any pending indel run.
                    flushRun();
                }
            });
        flushRun();

        // Convert RC-coordinate SV events to forward coordinates.
        if (qNeedsRc) {
            for (size_t ei = svEventsBefore; ei < scratch.svEvents.size(); ei++) {
                auto& ev = scratch.svEvents[ei];
                uint32_t fwdStart = queryLen - ev.queryEnd;
                uint32_t fwdEnd = queryLen - ev.queryPos;
                ev.queryPos = fwdStart;
                ev.queryEnd = fwdEnd;
            }
        }
    }

    if (scratch.svEvents.empty()) return;

    // ----------------------------------------------------------------
    // Step 2: Sort SV events by query position.
    // ----------------------------------------------------------------

    dinara::radixSort(scratch.svEvents.data(),
        scratch.svEvents.data() + scratch.svEvents.size(),
        [](const PhasingSvEvent& e) { return e.queryPos; });

    // ----------------------------------------------------------------
    // Step 3: Cluster by position overlap >= 50% with target-read dedup.
    //
    // BFS connected components: events are sorted by position, and we
    // greedily merge events that overlap the current cluster by >= 50%
    // of the smaller span. Within each cluster, deduplicate by target
    // read (keep first per target).
    // ----------------------------------------------------------------

    {
        const size_t nev = scratch.svEvents.size();

        // Assign cluster IDs via greedy merge.
        vector<int32_t> clusterId(nev, -1);
        int32_t nextCluster = 0;

        for (size_t i = 0; i < nev; i++) {
            if (clusterId[i] >= 0) continue;

            // Start a new cluster with event i.
            clusterId[i] = nextCluster;
            uint32_t cStart = scratch.svEvents[i].queryPos;
            uint32_t cEnd = scratch.svEvents[i].queryEnd;

            // BFS: try to add subsequent events.
            bool changed = true;
            while (changed) {
                changed = false;
                for (size_t j = i + 1; j < nev; j++) {
                    if (clusterId[j] >= 0) continue;
                    const auto& ev = scratch.svEvents[j];
                    if (ev.queryPos > cEnd) break; // sorted, no more overlap

                    uint32_t oStart = max(cStart, ev.queryPos);
                    uint32_t oEnd = min(cEnd, ev.queryEnd);
                    int32_t overlap = int32_t(oEnd) - int32_t(oStart);

                    uint32_t cSpan = cEnd - cStart;
                    uint32_t eSpan = ev.queryEnd - ev.queryPos;
                    uint32_t minSpan = min(cSpan, eSpan);

                    if (overlap > 0 && uint32_t(overlap) >= minSpan / 2) {
                        clusterId[j] = nextCluster;
                        cStart = min(cStart, ev.queryPos);
                        cEnd = max(cEnd, ev.queryEnd);
                        changed = true;
                    }
                }
            }

            // Build cluster with target-read dedup.
            PhasingSvCluster cluster;
            cluster.consensusPos = cStart;
            cluster.consensusEnd = cEnd;
            cluster.eventBegin = uint32_t(i);
            cluster.eventCount = 0;

            // Collect unique target reads in this cluster.
            vector<uint32_t> seenTargets;
            for (size_t j = i; j < nev; j++) {
                if (clusterId[j] != nextCluster) continue;
                const uint32_t tid = scratch.overlaps[scratch.svEvents[j].overlapIdx].targetReadId;
                bool dup = false;
                for (uint32_t t : seenTargets) {
                    if (t == tid) { dup = true; break; }
                }
                if (!dup) {
                    seenTargets.push_back(tid);
                    cluster.eventCount++;
                }
            }
            cluster.eventEnd = uint32_t(nev); // will scan by clusterId
            scratch.svClusters.push_back(cluster);
            nextCluster++;
        }

        // ---- Step 4: Label overlaps in qualifying clusters. ----
        // Require >= 3 unique targets in the cluster AND >= 3 cis overlaps
        // spanning the cluster's consensus region.
        for (int32_t ci = 0; ci < nextCluster; ci++) {
            const auto& cluster = scratch.svClusters[ci];
            if (cluster.eventCount < 3) continue; // need >= 3 unique targets (hifiasm c_sz=3)

            // Count how many cis overlaps span this region.
            uint32_t spanCount = 0;
            for (const auto& ov : scratch.overlaps) {
                if (ov.isMatch != 1) continue;
                if (ov.qs <= cluster.consensusPos && ov.qe >= cluster.consensusEnd) {
                    spanCount++;
                }
            }
            if (spanCount < 3) continue;

            // Overlaps with SV events in this cluster → trans.
            for (size_t ei = 0; ei < nev; ei++) {
                if (clusterId[ei] != ci) continue;
                auto& ov = scratch.overlaps[scratch.svEvents[ei].overlapIdx];
                if (ov.isMatch == 1) {
                    ov.isMatch = 2;
                }
            }
        }
    }
}

// ============================================================================
// dedupChains: reduce to one overlap per target read
//
// Port of hifiasm's dedup_chains (ecovlp.cpp:2984).
//
// For each target read, keep the best overlap by:
//   1. Lower is_match wins (cis preferred over trans)
//   2. Higher score wins: score = (qe - qs) - 12 * errors
// ============================================================================

static void dedupChains(
    PhasingScratchpad& scratch)
{
    if (scratch.overlaps.empty()) return;

    // Sort by (targetReadId, isMatch, score desc, span desc).
    // Matches hifiasm dedup_chains: lower isMatch wins (1=cis < 2=trans),
    // then higher score = span - 12*errors, then higher span as tiebreaker.
    sort(scratch.overlaps.begin(), scratch.overlaps.end(),
        [&](const PhasingOverlap& a, const PhasingOverlap& b) {
            if (a.targetReadId != b.targetReadId)
                return a.targetReadId < b.targetReadId;
            if (a.isMatch != b.isMatch)
                return a.isMatch < b.isMatch;
            // Quality-adjusted score: span - 12 * errors.
            int32_t spanA = int32_t(a.qe - a.qs);
            int32_t spanB = int32_t(b.qe - b.qs);
            int32_t scoreA = spanA - 12 * int32_t(a.errorCount);
            int32_t scoreB = spanB - 12 * int32_t(b.errorCount);
            if (scoreA != scoreB)
                return scoreA > scoreB;
            return spanA > spanB;
        });

    // Keep first per target (best by sort order).
    size_t write = 0;
    for (size_t i = 0; i < scratch.overlaps.size(); i++) {
        if (i > 0 &&
            scratch.overlaps[i].targetReadId ==
            scratch.overlaps[i-1].targetReadId) {
            continue;
        }
        scratch.overlaps[write++] = scratch.overlaps[i];
    }
    scratch.overlaps.resize(write);
}

// ============================================================================
// writeResults: write phasing labels back to AlignmentData
// ============================================================================

static void writeResults(
    Assembler& assembler,
    ReadId queryReadId,
    const PhasingScratchpad& scratch)
{
    for (const auto& ov : scratch.overlaps) {
        auto& ad = assembler.alignmentData[ov.alignmentId];

        // Per-read-perspective match state. Each thread writes to a different
        // field (state0 vs state1) based on queryReadId, so no race.
        ad.setHifiasmEcMatchStateFromReadPerspective(queryReadId, ov.isMatch);
    }
}

// ============================================================================
// phaseOverlaps: threaded entry point
//
// Iterates over all reads, running the full phasing pipeline per read.
// Thread-local PhasingScratchpad avoids per-read allocation.
// ============================================================================

void Assembler::phaseOverlaps(uint64_t threadCount)
{
    cout << timestamp << "=== ONT Overlap Phasing Pipeline ===" << endl;

    const uint64_t readCount = getReads().readCount();
    cout << timestamp << "Read count: " << readCount << endl;
    cout << timestamp << "Thread count: " << threadCount << endl;

    if (readCount == 0) {
        cout << timestamp << "No reads to phase." << endl;
        return;
    }

    // Counters for progress reporting.
    atomic<uint64_t> readsProcessed(0);
    atomic<uint64_t> readsWithOverlaps(0);
    atomic<uint64_t> readsWithSites(0);
    atomic<uint64_t> totalCis(0);
    atomic<uint64_t> totalTrans(0);

    // Mutex for debug output (prevents interleaving across threads).
    static std::mutex phasingCoutMutex;

    // Optional debug hook: set DINARA_PHASING_DEBUG_READ to a read ID (uint32)
    // to print detailed per-site and per-overlap accounting.
    ReadId phasingDebugReadId = invalid<ReadId>;
    {
        const char* s = std::getenv("DINARA_PHASING_DEBUG_READ");
        if (s && *s) {
            char* end = nullptr;
            const unsigned long v = std::strtoul(s, &end, 10);
            if (end && end != s && *end == 0
                && v <= std::numeric_limits<uint32_t>::max()
                && v < readCount) {
                phasingDebugReadId = ReadId(uint32_t(v));
                cout << timestamp << "[PHASING-DBG] Enabled for read "
                     << phasingDebugReadId << "-0" << endl;
            } else {
                cout << timestamp << "[PHASING-DBG] Requested read "
                     << s << " is invalid or out of range (readCount="
                     << readCount << ")" << endl;
            }
        }
    }

    // Per-thread timing accumulators (microseconds). Aggregated after join.
    struct alignas(64) ThreadTiming {
        int64_t gather = 0, unpack = 0, detectSnp = 0, buildMatrix = 0;
        int64_t filterAdj = 0, dpPhase = 0, labelCT = 0;
        int64_t largeIndel = 0, dedup = 0, write = 0;
    };
    vector<ThreadTiming> threadTimings(threadCount);

    // Static block scheduling (same pattern as performHifiasmECParity).
    vector<thread> threads;
    uint64_t chunkSize = readCount / threadCount;
    if (chunkSize == 0) chunkSize = 1;

    for (uint64_t t = 0; t < threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t start = t * chunkSize;
            if (start >= readCount) return;
            const uint64_t end = min((t + 1) * chunkSize, readCount);

            PhasingScratchpad scratch;
            ThreadTiming& tt = threadTimings[t];
            using clk = std::chrono::steady_clock;
            auto usec = [](clk::time_point a, clk::time_point b) {
                return std::chrono::duration_cast<std::chrono::microseconds>(b - a).count();
            };

            for (uint64_t rid = start; rid < end; rid++) {
                const ReadId readId(rid);
                scratch.clear();

                const bool isDebugRead =
                    (phasingDebugReadId != invalid<ReadId>
                     && readId == phasingDebugReadId);

                // 1. Gather overlaps.
                auto tp0 = clk::now();
                gatherOverlaps(*this, readId, scratch);
                auto tp1 = clk::now();
                tt.gather += usec(tp0, tp1);
                if (scratch.overlaps.empty()) {
                    if (isDebugRead) {
                        std::lock_guard<std::mutex> lock(phasingCoutMutex);
                        cout << timestamp << "[PHASING-DBG] read=" << readId
                             << " overlaps=0" << endl;
                    }
                    readsProcessed++;
                    continue;
                }
                readsWithOverlaps++;

                // 2. Unpack query sequence.
                tp0 = clk::now();
                const uint32_t queryLen =
                    uint32_t(getReads().getRead(readId).baseCount);
                unpackQuerySequence(*this, readId, queryLen, scratch);
                tp1 = clk::now();
                tt.unpack += usec(tp0, tp1);

                // 3. SNP detection.
                tp0 = clk::now();
                detectSnpSites(*this, readId, queryLen, scratch);
                tp1 = clk::now();
                tt.detectSnp += usec(tp0, tp1);

                // 4. Build SNP matrix (filter + confirm sites).
                tp0 = clk::now();
                buildSnpMatrix(scratch);
                tp1 = clk::now();
                tt.buildMatrix += usec(tp0, tp1);

                // Debug: candidate sites.
                if (isDebugRead) {
                    std::lock_guard<std::mutex> lock(phasingCoutMutex);
                    cout << timestamp << "[PHASING-DBG] read=" << readId
                         << " queryLen=" << queryLen
                         << " overlaps=" << scratch.overlaps.size()
                         << " candidateSites=" << scratch.sites.size()
                         << endl;
                    for (uint32_t si = 0; si < scratch.sites.size(); si++) {
                        const auto& s = scratch.sites[si];
                        cout << timestamp << "[PHASING-DBG]   candidate[" << si
                             << "] pos=" << s.site
                             << " ref=" << "ACGT"[s.queryBase]
                             << " alt=" << "ACGT"[s.altBase]
                             << " refSupport=" << s.matchCount
                             << " altSupport=" << s.altCount
                             << " fwdStrand=" << s.fwdStrandCount
                             << " isHpc=" << int(s.isHpc) << endl;
                    }
                }

                if (!scratch.sites.empty()) {
                    readsWithSites++;

                    // 5. DP phasing.
                    const uint64_t coverageHet =
                        assemblerInfo->kmerDistributionInfo.coverageHet;
                    const uint32_t hetCov = uint32_t(coverageHet / 2);
                    tp0 = clk::now();
                    runDpPhasing(scratch, hetCov);
                    tp1 = clk::now();
                    tt.dpPhase += usec(tp0, tp1);

                    // 5b. Adjacent-site filter (after DP, like hifiasm).
                    // Hifiasm runs this inside generate_haplotypes_naive_HiFi
                    // (Correct.cpp:8855), after gen_rphase_dp.
                    tp0 = clk::now();
                    filterAdjacentSites(scratch);
                    tp1 = clk::now();
                    tt.filterAdj += usec(tp0, tp1);

                    // 6. Allele grouping (cis/trans labeling).
                    tp0 = clk::now();
                    labelCisTrans(scratch);
                    tp1 = clk::now();
                    tt.labelCT += usec(tp0, tp1);

                    // Debug: confirmed sites and overlap labels.
                    if (isDebugRead) {
                        std::lock_guard<std::mutex> lock(phasingCoutMutex);
                        uint32_t confirmed = 0;
                        for (const auto& s : scratch.sites) {
                            if (s.dpChainId >= 0) confirmed++;
                        }
                        cout << timestamp << "[PHASING-DBG] coverageHet="
                             << coverageHet << " hetCov=" << hetCov
                             << " confirmed=" << confirmed
                             << "/" << scratch.sites.size() << endl;
                        for (uint32_t si = 0; si < scratch.sites.size(); si++) {
                            const auto& s = scratch.sites[si];
                            if (s.dpChainId < 0) continue;
                            cout << timestamp << "[PHASING-DBG]   confirmed["
                                 << si << "] pos=" << s.site
                                 << " ref=" << "ACGT"[s.queryBase]
                                 << " alt=" << "ACGT"[s.altBase]
                                 << " refSupport=" << s.matchCount
                                 << " altSupport=" << s.altCount
                                 << " fwdStrand=" << s.fwdStrandCount
                                 << " chain=" << s.dpChainId
                                 << " isHpc=" << int(s.isHpc)
                                 << " labelConfirmed="
                                 << s.isLabelConfirmed() << endl;
                        }
                        uint32_t cisCount = 0, transCount = 0;
                        for (const auto& ov : scratch.overlaps) {
                            if (ov.isMatch == 1) cisCount++;
                            else if (ov.isMatch == 2) transCount++;
                        }
                        cout << timestamp << "[PHASING-DBG] labels cis="
                             << cisCount << " trans=" << transCount
                             << " unlabeled="
                             << (scratch.overlaps.size() - cisCount - transCount)
                             << endl;
                        for (const auto& ov : scratch.overlaps) {
                            const char* label =
                                (ov.isMatch == 1) ? "CIS" :
                                (ov.isMatch == 2) ? "TRANS" : "UNLABELED";
                            cout << timestamp << "[PHASING-DBG]   target="
                                 << ov.targetReadId
                                 << " isRev=" << int(ov.isRev)
                                 << " strong=" << int(ov.strong)
                                 << " confirmedSites=" << ov.confirmedSiteCount
                                 << " mismatches=" << ov.confirmedMismatchCount
                                 << " -> " << label << endl;
                        }
                    }

                    // 7. Large indel phasing.
                    tp0 = clk::now();
                    phaseLargeIndels(*this, readId, queryLen, scratch);
                    tp1 = clk::now();
                    tt.largeIndel += usec(tp0, tp1);
                }

                // 8. Dedup chains (best overlap per target).
                tp0 = clk::now();
                dedupChains(scratch);
                tp1 = clk::now();
                tt.dedup += usec(tp0, tp1);

                // 9. Write results back.
                tp0 = clk::now();
                writeResults(*this, readId, scratch);
                tp1 = clk::now();
                tt.write += usec(tp0, tp1);

                // Count labels for reporting.
                for (const auto& ov : scratch.overlaps) {
                    if (ov.isMatch == 1) totalCis++;
                    else if (ov.isMatch == 2) totalTrans++;
                }

                readsProcessed++;

                // Progress reporting every 10000 reads.
                const uint64_t n = readsProcessed.load();
                if (n % 10000 == 0) {
                    cout << timestamp << "Phased " << n << " / "
                         << readCount << " reads" << endl;
                }
            }
        });
    }

    for (auto& t : threads) t.join();

    // Aggregate and print per-stage timing.
    ThreadTiming total;
    for (const auto& tt : threadTimings) {
        total.gather += tt.gather;
        total.unpack += tt.unpack;
        total.detectSnp += tt.detectSnp;
        total.buildMatrix += tt.buildMatrix;
        total.filterAdj += tt.filterAdj;
        total.dpPhase += tt.dpPhase;
        total.labelCT += tt.labelCT;
        total.largeIndel += tt.largeIndel;
        total.dedup += tt.dedup;
        total.write += tt.write;
    }
    auto ms = [](int64_t us) { return us / 1000; };
    cout << timestamp << "Phasing per-stage timing (ms, sum over " << threadCount << " threads):" << endl;
    cout << timestamp << "  gatherOverlaps:    " << ms(total.gather) << endl;
    cout << timestamp << "  unpackQuery:       " << ms(total.unpack) << endl;
    cout << timestamp << "  detectSnpSites:    " << ms(total.detectSnp) << endl;
    cout << timestamp << "  buildSnpMatrix:    " << ms(total.buildMatrix) << endl;
    cout << timestamp << "  filterAdjacentSites:" << ms(total.filterAdj) << endl;
    cout << timestamp << "  runDpPhasing:      " << ms(total.dpPhase) << endl;
    cout << timestamp << "  labelCisTrans:     " << ms(total.labelCT) << endl;
    cout << timestamp << "  phaseLargeIndels:  " << ms(total.largeIndel) << endl;
    cout << timestamp << "  dedupChains:       " << ms(total.dedup) << endl;
    cout << timestamp << "  writeResults:      " << ms(total.write) << endl;

    cout << timestamp << "Phasing complete." << endl;
    cout << timestamp << "Reads processed: " << readsProcessed.load() << endl;
    cout << timestamp << "Reads with overlaps: " << readsWithOverlaps.load() << endl;
    cout << timestamp << "Reads with SNP sites: " << readsWithSites.load() << endl;
    cout << timestamp << "Total cis labels: " << totalCis.load() << endl;
    cout << timestamp << "Total trans labels: " << totalTrans.load() << endl;
}
