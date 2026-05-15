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
#include "timestamp.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
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
        ov.errorCount = (info.nonHomopolymerErrorCount != uint32_t(-1))
            ? info.nonHomopolymerErrorCount
            : (info.mismatchCount + info.gapCount);
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
// region of period 1..HPC_RR (4). For each period r, extends the
// repeat in both directions from p. If the repeat span exceeds
// r * HPC_CC (cutoff), the position is masked.
// ============================================================================

static constexpr uint32_t HPC_RR = 4;  // max repeat period
static constexpr uint32_t HPC_CC = 2;  // cutoff multiplier

static bool isPeriodicRepeat(
    const uint8_t* seq,
    uint32_t seqLen,
    uint32_t p)
{
    if (seqLen == 0) return false;

    const int64_t sn = int64_t(seqLen);
    const int64_t pp = int64_t(p);

    for (uint32_t r = 1; r <= HPC_RR; r++) {
        const int64_t rc = int64_t(r) * int64_t(HPC_CC); // cutoff

        // Scan 0: including p, extend right then left.
        {
            int64_t k = pp + int64_t(r);
            while (k < sn && (k - int64_t(r)) >= 0 && seq[k] == seq[k - r]) k++;
            int64_t ze = k;
            k = pp - 1;
            while (k >= 0 && (k + int64_t(r)) < sn && seq[k] == seq[k + r]) k--;
            int64_t zs = k + 1;
            if ((ze - zs) > int64_t(r) && (ze - zs) >= rc) return true;
        }

        // Scan 1: not including p, right side only.
        {
            int64_t k = pp + int64_t(r) + 1;
            while (k < sn && (k - int64_t(r)) >= 0 && seq[k] == seq[k - r]) k++;
            int64_t zs = pp + 1;
            int64_t ze = k;
            if ((ze - zs) > int64_t(r) && (ze - zs) >= rc) return true;
        }

        // Scan 2: including p, extend left then right.
        {
            int64_t k = pp - int64_t(r);
            while (k >= 0 && (k + int64_t(r)) < sn && seq[k] == seq[k + r]) k--;
            int64_t zs = k + 1;
            k = pp + 1;
            while (k < sn && (k - int64_t(r)) >= 0 && seq[k] == seq[k - r]) k++;
            int64_t ze = k;
            if ((ze - zs) > int64_t(r) && (ze - zs) >= rc) return true;
        }

        // Scan 3: not including p, left side only.
        {
            int64_t k = pp - int64_t(r) - 1;
            while (k >= 0 && (k + int64_t(r)) < sn && seq[k] == seq[k + r]) k--;
            int64_t zs = k + 1;
            int64_t ze = pp;
            if ((ze - zs) > int64_t(r) && (ze - zs) >= rc) return true;
        }
    }

    return false;
}

// ============================================================================
// detectSnpSites: sliding-window SNP detection
//
// Port of hifiasm's hc_phase_robust_rr (Correct.cpp:10200).
//
// Three stages:
//   Pre-walk — Walk each overlap's CIGAR once, collecting match/mismatch
//              events into scratch.cigarEvents with per-overlap ranges.
//   Pass 1   — For each 428 bp window, count mismatch votes per position.
//              Positions with >= 2 votes become candidate SNPs (unless
//              masked by periodic repeat detection).
//   Pass 2   — Re-scan events at candidate positions, emitting
//              PhasingEvidence entries with the observed base and
//              match/mismatch status.
// ============================================================================

static void detectSnpSites(
    const Assembler& assembler,
    ReadId queryReadId,
    uint32_t queryLen,
    PhasingScratchpad& scratch)
{
    const auto& cigarStore = assembler.getOverlapCigarStore();
    const uint32_t W = PHASING_WINDOW_SIZE;

    // ----------------------------------------------------------------
    // Pre-walk: collect all match/mismatch events per overlap.
    //
    // walkRangeWithCursor filters by xk (read0) coordinates, which only
    // works when queryIsRead0. To handle both orientations uniformly, we
    // do a single full walk per overlap and store events in forward-strand
    // query coordinates. Each overlap's events are contiguous and sorted
    // by qpos, enabling efficient window-based scanning below.
    // ----------------------------------------------------------------

    auto& allEvents = scratch.cigarEvents;
    scratch.cigarEventRanges.resize(scratch.overlaps.size());

    allEvents.reserve(scratch.overlaps.size() * 64); // rough estimate

    for (size_t oi = 0; oi < scratch.overlaps.size(); oi++) {
        auto& ov = scratch.overlaps[oi];
        const auto& ad = assembler.alignmentData[ov.alignmentId];

        scratch.cigarEventRanges[oi].begin = uint32_t(allEvents.size());

        OverlapCigarStore::Cursor cursor;
        cursor.reset(ov.cigarOffset, ov.cigarTokenCount,
            ad.qs, cigarRead1Start(assembler, ad), cigarStore);

        // When queryIsRead0 == 0 and the overlap is RC, yk is in RC
        // coordinates of read1 (the query). queryBases is in forward
        // coordinates, so we must convert: qpos_fwd = queryLen - 1 - qpos_rc.
        const bool queryNeedsRcConvert =
            (ov.queryIsRead0 == 0) && (ov.isRev != 0);

        cigarStore.walkRangeWithCursor(
            cursor, 0, UINT32_MAX,
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (op != 0 && op != 1) return; // skip indels

                for (uint32_t b = 0; b < len; b++) {
                    uint32_t qpos, tpos;
                    if (ov.queryIsRead0) {
                        qpos = uint32_t(xk) + b;
                        tpos = uint32_t(yk) + b;
                    } else {
                        qpos = uint32_t(yk) + b;
                        tpos = uint32_t(xk) + b;
                        if (queryNeedsRcConvert) {
                            if (qpos >= queryLen) continue;
                            qpos = queryLen - 1 - qpos;
                        }
                    }
                    if (qpos >= queryLen) continue;

                    CigarEvent ce;
                    ce.qpos = qpos;
                    ce.tpos = tpos;
                    ce.op = op;
                    allEvents.push_back(ce);
                }
            });

        scratch.cigarEventRanges[oi].end = uint32_t(allEvents.size());

        // For queryIsRead0==0 with RC, the qpos values are reverse-ordered
        // after forward conversion. Sort to enable early-exit in window scans.
        if (queryNeedsRcConvert) {
            sort(allEvents.begin() + scratch.cigarEventRanges[oi].begin,
                 allEvents.begin() + scratch.cigarEventRanges[oi].end,
                 [](const CigarEvent& a, const CigarEvent& b) {
                     return a.qpos < b.qpos;
                 });
        }
    }

    // ----------------------------------------------------------------
    // Windowed SNP detection using pre-collected events.
    // ----------------------------------------------------------------

    for (uint32_t wStart = 0; wStart < queryLen; wStart += W) {
        const uint32_t wEnd = min(wStart + W, queryLen);
        const uint32_t wLen = wEnd - wStart;

        // ---- Pass 1: vote counting ----
        // Count how many overlaps have a mismatch at each position.
        // Positions with > OCC_THRES (1) votes are candidate SNPs.
        scratch.flag.assign(wLen, 0);

        for (size_t oi = 0; oi < scratch.overlaps.size(); oi++) {
            const auto& ov = scratch.overlaps[oi];
            if (ov.qe <= wStart || ov.qs >= wEnd) continue;

            for (uint32_t ei = scratch.cigarEventRanges[oi].begin;
                 ei < scratch.cigarEventRanges[oi].end; ei++) {
                const auto& ce = allEvents[ei];
                if (ce.qpos < wStart) continue;
                if (ce.qpos >= wEnd) break; // events are sorted by qpos
                if (ce.op == 1) { // mismatch
                    uint32_t fi = ce.qpos - wStart;
                    if (scratch.flag[fi] < 255)
                        scratch.flag[fi]++;
                }
            }
        }

        // Classify positions:
        //   flag=0 → not a candidate (too few votes)
        //   flag=1 → candidate SNP
        //   flag=3 → masked by periodic repeat (hifiasm hpc_mask_ff)
        bool anyCandidates = false;
        for (uint32_t i = 0; i < wLen; i++) {
            if (scratch.flag[i] > PHASING_OCC_THRES) {
                const uint32_t absPos = wStart + i;
                if (isPeriodicRepeat(scratch.queryBases.data(), queryLen,
                                     absPos)) {
                    scratch.flag[i] = 3;
                } else {
                    scratch.flag[i] = 1;
                    anyCandidates = true;
                }
            } else {
                scratch.flag[i] = 0;
            }
        }
        if (!anyCandidates) continue;

        // ---- Pass 2: evidence collection ----
        // At candidate positions (flag==1), emit PhasingEvidence entries
        // recording the observed base and whether it matches the query.
        for (size_t oi = 0; oi < scratch.overlaps.size(); oi++) {
            const auto& ov = scratch.overlaps[oi];
            if (ov.qe <= wStart || ov.qs >= wEnd) continue;

            for (uint32_t ei = scratch.cigarEventRanges[oi].begin;
                 ei < scratch.cigarEventRanges[oi].end; ei++) {
                const auto& ce = allEvents[ei];
                if (ce.qpos < wStart) continue;
                if (ce.qpos >= wEnd) break;

                uint32_t fi = ce.qpos - wStart;
                if (scratch.flag[fi] == 0 || scratch.flag[fi] == 3) continue;

                PhasingEvidence ev;
                ev.site = ce.qpos;
                ev.overlapIdx = uint32_t(oi);

                if (ce.op == 0) { // match → query base
                    ev.base = scratch.queryBases[ce.qpos];
                    ev.isAlt = 0;
                } else { // mismatch → look up target base
                    // tpos is in oriented coordinates of the target read.
                    // When queryIsRead0: target is read1, RC if isRev.
                    // When !queryIsRead0: target is read0, always forward.
                    const bool targetIsRc = (ov.queryIsRead0 != 0) && (ov.isRev != 0);
                    ev.base = getBaseAtPosition(
                        assembler,
                        ReadId(ov.targetReadId),
                        ce.tpos, targetIsRc);
                    ev.isAlt = 1;
                }

                ev.isHpc = 0;
                scratch.evidence.push_back(ev);
            }
        }
    }
}

// ============================================================================
// buildSnpMatrix: sort evidence, group by site, filter, confirm sites
//
// Port of hifiasm's SetSnpMatrix + push_info (Correct.cpp:10511).
//
// 1. Sort evidence by (site, overlapIdx).
// 2. Group by site. For each site:
//    - Dedup by overlapIdx (keep first per overlap).
//    - Count matchCount (isAlt==0) and per-base altCount[4].
//    - Require matchCount >= PHASING_S_HAP_COV and altCount >= PHASING_INFOR_COV.
//    - For each qualifying alt base, create a PhasingSite.
// 3. Record evidence ranges per site.
// ============================================================================

static void buildSnpMatrix(
    PhasingScratchpad& scratch)
{
    if (scratch.evidence.empty()) return;

    // Sort by (site, overlapIdx).
    sort(scratch.evidence.begin(), scratch.evidence.end(),
        [](const PhasingEvidence& a, const PhasingEvidence& b) {
            if (a.site != b.site) return a.site < b.site;
            return a.overlapIdx < b.overlapIdx;
        });

    // Dedup: within each (site, overlapIdx) group, keep only the first entry.
    {
        size_t write = 0;
        for (size_t i = 0; i < scratch.evidence.size(); i++) {
            if (i > 0 &&
                scratch.evidence[i].site == scratch.evidence[i-1].site &&
                scratch.evidence[i].overlapIdx == scratch.evidence[i-1].overlapIdx) {
                continue; // skip duplicate
            }
            scratch.evidence[write++] = scratch.evidence[i];
        }
        scratch.evidence.resize(write);
    }

    // Group by site and build PhasingSite entries.
    size_t i = 0;
    while (i < scratch.evidence.size()) {
        const uint32_t site = scratch.evidence[i].site;
        const size_t groupBegin = i;

        // Find end of this site's group.
        while (i < scratch.evidence.size() && scratch.evidence[i].site == site) {
            i++;
        }
        const size_t groupEnd = i;

        // Count match and per-base alt.
        uint32_t matchCount = 0;
        uint32_t fwdStrandRefCount = 0; // forward-strand ref-matching overlaps
        uint32_t altCount[4] = {0, 0, 0, 0};
        uint8_t isHpc = 0;
        uint8_t queryBase = 0;

        for (size_t j = groupBegin; j < groupEnd; j++) {
            const auto& ev = scratch.evidence[j];
            if (j == groupBegin) {
                queryBase = scratch.queryBases[site];
                isHpc = ev.isHpc;
            }
            if (ev.isAlt == 0) {
                matchCount++;
                // Track forward-strand count for strand bias.
                // isRev==0 means same-strand (forward).
                if (scratch.overlaps[ev.overlapIdx].isRev == 0) {
                    fwdStrandRefCount++;
                }
            } else {
                if (ev.base < 4) altCount[ev.base]++;
            }
            if (ev.isHpc) isHpc = 1;
        }

        // Hifiasm adds +1 to matchCount for the query read itself.
        matchCount += 1;

        // Create a PhasingSite for each qualifying alt base (count >= 2).
        // Hifiasm creates separate SnpStats entries per alt allele.
        for (uint8_t b = 0; b < 4; b++) {
            if (b == queryBase) continue;
            if (altCount[b] < 2) continue; // hifiasm: occ_1[i] >= 2

            // Filter: require sufficient evidence on both alleles.
            if (matchCount < PHASING_S_HAP_COV) continue;
            if (altCount[b] < PHASING_INFOR_COV) continue;

            PhasingSite ps;
            ps.site = site;
            ps.queryBase = queryBase;
            ps.altBase = b;
            ps.isHpc = isHpc;
            ps.dpChainId = -1;
            ps.matchCount = matchCount;
            ps.altCount = altCount[b];
            ps.fwdStrandCount = fwdStrandRefCount + 1; // +1 for query read
            ps.evidenceBegin = uint32_t(groupBegin);
            ps.evidenceEnd = uint32_t(groupEnd);

            scratch.sites.push_back(ps);
        }
    }
}

// ============================================================================
// filterAdjacentSites: remove SNP sites at adjacent positions
//
// Port of hifiasm's adjacent-site filter (Correct.cpp:8855).
// If sites exist at positions p and p+1, BOTH are removed.
// This prevents alignment artifacts (e.g., indels manifesting as
// adjacent mismatches) from being treated as het sites.
// ============================================================================

static void filterAdjacentSites(
    PhasingScratchpad& scratch)
{
    if (scratch.sites.size() < 2) return;

    // Sites are already sorted by position (from buildSnpMatrix's
    // sequential processing of sorted evidence).
    // With multi-alt, multiple sites can share the same position.
    // Mark sites for removal if their position is adjacent to another site's.

    const size_t n = scratch.sites.size();
    vector<bool> remove(n, false);

    for (size_t i = 0; i < n; i++) {
        // Check predecessor (different position, adjacent).
        if (i > 0 && scratch.sites[i].site == scratch.sites[i-1].site + 1) {
            remove[i] = true;
            // Also mark all sites at the predecessor position.
            for (size_t j = i - 1; j < n && scratch.sites[j].site == scratch.sites[i].site - 1; ) {
                remove[j] = true;
                if (j == 0) break;
                j--;
            }
        }
        // Check successor (different position, adjacent).
        if (i + 1 < n && scratch.sites[i + 1].site == scratch.sites[i].site + 1) {
            remove[i] = true;
        }
    }

    // Compact.
    size_t write = 0;
    for (size_t i = 0; i < n; i++) {
        if (!remove[i]) {
            scratch.sites[write++] = scratch.sites[i];
        }
    }
    scratch.sites.resize(write);
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

    // Walk both evidence ranges (sorted by overlapIdx) in tandem.
    uint32_t pi = si.evidenceBegin;
    uint32_t pj = sj.evidenceBegin;

    // Count shared overlaps where both sites agree:
    // nn0 = both match query (isAlt==0 at both sites)
    // nn1 = both mismatch query (isAlt==1 at both sites)
    // Hifiasm requires nn0 > 0 && nn1 > 0: evidence from both haplotypes
    // is needed to confirm the sites are on the same haplotype block.
    uint32_t nn0 = 0, nn1 = 0;

    while (pi < si.evidenceEnd && pj < sj.evidenceEnd) {
        const auto& ei = scratch.evidence[pi];
        const auto& ej = scratch.evidence[pj];

        if (ei.overlapIdx < ej.overlapIdx) { pi++; continue; }
        if (ei.overlapIdx > ej.overlapIdx) { pj++; continue; }

        // Same overlap covers both sites.
        // Both match → nn0++. Both mismatch → nn1++.
        // Mixed (one match, one mismatch) → skip (not counted).
        if (ei.isAlt == ej.isAlt) {
            if (ei.isAlt == 0) nn0++;
            else nn1++;
        }

        pi++;
        pj++;
    }

    // Require evidence from both haplotypes.
    return nn0 > 0 && nn1 > 0;
}

// ============================================================================
// runDpPhasing: longest compatible chain over confirmed SNP sites
//
// Port of hifiasm's gen_rphase_dp0_single_path (Correct.cpp:9648).
//
// Standard LIS-style DP: f[i] = max(f[j] + 1) for all j < i where
// sites i and j are compatible. Backtrack to extract chains.
// Sites in chains of length > 1 are confirmed. Length-1 chains require
// stricter thresholds (matchCount >= hetCov * 0.7, minimum 6).
// ============================================================================

static void runDpPhasing(
    PhasingScratchpad& scratch,
    uint32_t hetCov)
{
    const uint32_t n = uint32_t(scratch.sites.size());
    if (n == 0) return;

    scratch.dpScore.assign(n, 1);
    scratch.dpParent.assign(n, -1);

    // DP: O(n^2) but n is typically small (tens to low hundreds of SNP sites).
    for (uint32_t i = 1; i < n; i++) {
        for (uint32_t j = 0; j < i; j++) {
            if (scratch.dpScore[j] + 1 > scratch.dpScore[i]) {
                if (checkCompatibility(scratch, i, j)) {
                    scratch.dpScore[i] = scratch.dpScore[j] + 1;
                    scratch.dpParent[i] = int32_t(j);
                }
            }
        }
    }

    // Extract chains by backtracking from each local maximum.
    scratch.dpChainId.assign(n, -1);
    int8_t chainId = 0;

    // Find chain endpoints: sites where no later site points to them.
    scratch.dpIsEndpoint.assign(n, true);
    for (uint32_t i = 0; i < n; i++) {
        if (scratch.dpParent[i] >= 0) {
            scratch.dpIsEndpoint[uint32_t(scratch.dpParent[i])] = false;
        }
    }

    // Extract chains from longest to shortest.
    scratch.dpEndpoints.clear();
    for (uint32_t i = 0; i < n; i++) {
        if (scratch.dpIsEndpoint[i]) {
            scratch.dpEndpoints.push_back({scratch.dpScore[i], i});
        }
    }
    sort(scratch.dpEndpoints.begin(), scratch.dpEndpoints.end(),
        [](const auto& a, const auto& b) { return a.first > b.first; });

    for (const auto& [score, endIdx] : scratch.dpEndpoints) {
        if (uint32_t(chainId) >= PHASING_MAX_DP_CHAINS) break;
        if (scratch.dpChainId[endIdx] >= 0) continue; // already claimed

        // Backtrack.
        scratch.dpChain.clear();
        int32_t k = int32_t(endIdx);
        while (k >= 0 && scratch.dpChainId[uint32_t(k)] < 0) {
            scratch.dpChain.push_back(uint32_t(k));
            k = scratch.dpParent[uint32_t(k)];
        }

        if (scratch.dpChain.size() < 2) {
            // Single-site chain: apply stricter threshold.
            uint32_t idx = scratch.dpChain[0];
            uint32_t minMatch = max(6U, uint32_t(double(hetCov) * 0.7));
            if (scratch.sites[idx].matchCount < minMatch) continue;
            if (scratch.sites[idx].altCount < PHASING_INFOR_COV) continue;
        }

        // Assign chain ID.
        for (uint32_t idx : scratch.dpChain) {
            scratch.dpChainId[idx] = chainId;
            scratch.sites[idx].dpChainId = chainId;
        }
        chainId++;
    }
}

// ============================================================================
// isStrandBiased: hifiasm is_st_bs macro equivalent
//
// Returns true if a site has strand bias — nearly all ref-matching
// overlaps come from one strand. Hifiasm parameters: st_rate=0.05, st_max=2.
// ============================================================================

static constexpr double PHASING_ST_RATE = 0.05;
static constexpr uint32_t PHASING_ST_MAX = 2;

static bool isStrandBiased(const PhasingSite& site)
{
    // fwdStrandCount = forward-strand ref overlaps + 1 (query)
    // matchCount = total ref overlaps + 1 (query)
    // Condition: reverse-strand count <= ST_MAX AND forward fraction >= (1 - ST_RATE)
    if (site.fwdStrandCount + PHASING_ST_MAX >= site.matchCount &&
        site.matchCount * PHASING_ST_RATE + site.fwdStrandCount >= site.matchCount) {
        return true;
    }
    return false;
}

// ============================================================================
// labelCisTrans: greedy cis/trans labeling from confirmed SNP sites
//
// Port of hifiasm's generate_haplotypes_naive_HiFi (Correct.cpp:8845).
//
// Four phases:
//   Phase 1 — Count confirmed mismatches per overlap (at DP-confirmed sites).
//   Phase 2 — Greedy labeling: overlaps sorted by mismatch count desc;
//             those with mismatches → trans, others → cis.
//   Phase 3 — Consistency check: cis overlaps with mismatches at
//             trans-confirmed sites get flipped to trans.
//   Phase 4 — multi_check: promote weak (sub-threshold) sites that
//             appear across multiple overlaps, flipping additional
//             cis overlaps to trans.
// ============================================================================

static void labelCisTrans(
    PhasingScratchpad& scratch)
{
    // ----------------------------------------------------------------
    // Phase 1: Count confirmed mismatches per overlap.
    //
    // For each DP-confirmed, non-strand-biased site, walk its evidence
    // and tally how many confirmed sites each overlap covers
    // (confirmedSiteCount) and how many of those are mismatches
    // (confirmedMismatchCount).
    // ----------------------------------------------------------------

    for (auto& ov : scratch.overlaps) {
        ov.confirmedSiteCount = 0;
        ov.confirmedMismatchCount = 0;
    }

    for (uint32_t si = 0; si < uint32_t(scratch.sites.size()); si++) {
        const auto& site = scratch.sites[si];
        if (site.dpChainId < 0) continue; // not confirmed by DP
        if (isStrandBiased(site)) continue;

        for (uint32_t ei = site.evidenceBegin; ei < site.evidenceEnd; ei++) {
            const auto& ev = scratch.evidence[ei];
            auto& ov = scratch.overlaps[ev.overlapIdx];
            ov.confirmedSiteCount++;
            if (ev.isAlt) {
                ov.confirmedMismatchCount++;
            }
        }
    }

    // ----------------------------------------------------------------
    // Phase 2: Greedy labeling (most mismatches first).
    //
    // Sort overlaps by confirmedMismatchCount descending. Overlaps with
    // at least one confirmed mismatch → trans (is_match=2). Overlaps
    // covering confirmed sites but with zero mismatches → cis (is_match=1).
    // Overlaps not covering any confirmed site keep their default (cis).
    // ----------------------------------------------------------------

    scratch.sortedOverlapIndices.resize(scratch.overlaps.size());
    iota(scratch.sortedOverlapIndices.begin(),
         scratch.sortedOverlapIndices.end(), 0U);
    sort(scratch.sortedOverlapIndices.begin(),
         scratch.sortedOverlapIndices.end(),
         [&](uint32_t a, uint32_t b) {
             return scratch.overlaps[a].confirmedMismatchCount >
                    scratch.overlaps[b].confirmedMismatchCount;
         });

    for (uint32_t idx : scratch.sortedOverlapIndices) {
        auto& ov = scratch.overlaps[idx];
        if (ov.confirmedSiteCount == 0) {
            ov.isMatch = 1;
            continue;
        }
        if (ov.confirmedMismatchCount > 0) {
            ov.isMatch = 2; // trans
            ov.strong = 1;
        } else {
            ov.isMatch = 1; // cis
        }
    }

    // ----------------------------------------------------------------
    // Phase 3: Consistency check.
    //
    // Build a set of "trans-confirmed" sites: confirmed sites where at
    // least one trans overlap has a mismatch. Then re-check all cis
    // overlaps — if a cis overlap has a mismatch at any trans-confirmed
    // site, flip it to trans.
    // ----------------------------------------------------------------

    // scratch.flag[si] = 1 means site si is trans-confirmed.
    scratch.flag.assign(scratch.sites.size(), 0);

    // Mark trans-confirmed sites.
    for (uint32_t si = 0; si < uint32_t(scratch.sites.size()); si++) {
        const auto& site = scratch.sites[si];
        if (site.dpChainId < 0) continue;
        for (uint32_t ei = site.evidenceBegin; ei < site.evidenceEnd; ei++) {
            const auto& ev = scratch.evidence[ei];
            if (ev.isAlt && scratch.overlaps[ev.overlapIdx].isMatch == 2) {
                scratch.flag[si] = 1;
                break;
            }
        }
    }

    // Flip cis overlaps that mismatch at trans-confirmed sites.
    for (uint32_t si = 0; si < uint32_t(scratch.sites.size()); si++) {
        if (scratch.flag[si] == 0) continue;
        const auto& site = scratch.sites[si];
        for (uint32_t ei = site.evidenceBegin; ei < site.evidenceEnd; ei++) {
            const auto& ev = scratch.evidence[ei];
            auto& ov = scratch.overlaps[ev.overlapIdx];
            if (ov.isMatch == 1 && ev.isAlt) {
                ov.isMatch = 2; // flip cis → trans
            }
        }
    }

    // ----------------------------------------------------------------
    // Phase 4: multi_check — weak site promotion.
    //
    // Catches overlaps that the DP missed because individual sites
    // didn't meet the strict threshold. Steps:
    //   a) Identify "weak sites": evidence positions with occ_0 >= 2
    //      and occ_1 >= 2 that are NOT already DP-confirmed.
    //   b) For each cis overlap, collect its weak mismatch positions.
    //      Require count >= alignLength * 0.04.
    //   c) Apply 32bp proximity filter (remove clustered sites).
    //      Require >= 2 surviving positions per overlap.
    //   d) Positions appearing in >= 2 overlaps get "promoted".
    //   e) Cis overlaps with mismatches at promoted positions → trans.
    // ----------------------------------------------------------------

    {
        // (a) Build set of confirmed site positions for exclusion.
        vector<uint32_t> confirmedPositions;
        for (const auto& site : scratch.sites) {
            if (site.dpChainId >= 0) {
                confirmedPositions.push_back(site.site);
            }
        }
        sort(confirmedPositions.begin(), confirmedPositions.end());

        // Identify weak sites by grouping evidence by position.
        struct WeakSiteInfo {
            uint32_t site;
            uint32_t matchCount;
            uint32_t altCount;
        };
        vector<WeakSiteInfo> weakSites;

        size_t ei = 0;
        while (ei < scratch.evidence.size()) {
            const uint32_t pos = scratch.evidence[ei].site;
            size_t groupEnd = ei;
            uint32_t mc = 0, ac = 0;
            while (groupEnd < scratch.evidence.size() &&
                   scratch.evidence[groupEnd].site == pos) {
                if (scratch.evidence[groupEnd].isAlt) ac++;
                else mc++;
                groupEnd++;
            }
            mc += 1; // +1 for query read

            bool isConfirmed = binary_search(
                confirmedPositions.begin(), confirmedPositions.end(), pos);
            if (!isConfirmed && mc >= 2 && ac >= 2) {
                weakSites.push_back({pos, mc, ac});
            }
            ei = groupEnd;
        }

        if (!weakSites.empty()) {
            // (b) Build per-overlap list of weak mismatch positions.
            vector<uint32_t> weakPositions;
            for (const auto& ws : weakSites) {
                weakPositions.push_back(ws.site);
            }

            const uint32_t numOv = uint32_t(scratch.overlaps.size());
            vector<vector<uint32_t>> overlapWeakMismatches(numOv);
            for (const auto& ev : scratch.evidence) {
                if (!ev.isAlt) continue;
                if (ev.overlapIdx >= numOv) continue;
                if (binary_search(weakPositions.begin(), weakPositions.end(), ev.site)) {
                    overlapWeakMismatches[ev.overlapIdx].push_back(ev.site);
                }
            }

            // (c-d) For each cis overlap, apply density + proximity filters,
            //       then collect promotion candidates.
            vector<uint32_t> promotionCandidates;

            for (uint32_t oi = 0; oi < numOv; oi++) {
                const auto& ov = scratch.overlaps[oi];
                if (ov.isMatch == 2) continue; // already trans

                auto& overlapWeakSites = overlapWeakMismatches[oi];
                if (overlapWeakSites.empty()) continue;

                // Density filter: need >= 4% of alignment length.
                const uint32_t alignLen = ov.qe - ov.qs;
                if (overlapWeakSites.size() < uint32_t(double(alignLen) * 0.04)) continue;

                // 32bp proximity filter: remove sites within 32bp of neighbors.
                sort(overlapWeakSites.begin(), overlapWeakSites.end());
                vector<uint32_t> filtered;
                for (size_t i = 0; i < overlapWeakSites.size(); i++) {
                    bool tooClose = false;
                    if (i > 0 && overlapWeakSites[i] - overlapWeakSites[i-1] < 32)
                        tooClose = true;
                    if (i + 1 < overlapWeakSites.size() &&
                        overlapWeakSites[i+1] - overlapWeakSites[i] < 32)
                        tooClose = true;
                    if (!tooClose) filtered.push_back(overlapWeakSites[i]);
                }

                if (filtered.size() >= 2) {
                    for (uint32_t pos : filtered) {
                        promotionCandidates.push_back(pos);
                    }
                }
            }

            // (d) Positions appearing in >= 2 overlaps get promoted.
            sort(promotionCandidates.begin(), promotionCandidates.end());
            vector<uint32_t> promotedPositions;
            for (size_t i = 0; i < promotionCandidates.size(); ) {
                const uint32_t pos = promotionCandidates[i];
                size_t count = 0;
                while (i < promotionCandidates.size() &&
                       promotionCandidates[i] == pos) {
                    count++;
                    i++;
                }
                if (count >= 2) {
                    promotedPositions.push_back(pos);
                }
            }

            // (e) Flip cis overlaps with mismatches at promoted positions.
            if (!promotedPositions.empty()) {
                for (const auto& ev : scratch.evidence) {
                    if (!ev.isAlt) continue;
                    auto& ov = scratch.overlaps[ev.overlapIdx];
                    if (ov.isMatch != 1) continue;
                    if (binary_search(promotedPositions.begin(),
                                      promotedPositions.end(), ev.site)) {
                        ov.isMatch = 2;
                    }
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

                if (op == 2 || op == 3) { // insertion or deletion
                    if (indelRunStart == UINT32_MAX) {
                        indelRunStart = qpos;
                    }
                    // Compute how many query bases this op spans.
                    // When queryIsRead0: query=xk, op3 (del) advances xk, op2 (ins) doesn't.
                    // When !queryIsRead0: query=yk, op2 (ins) advances yk, op3 (del) doesn't.
                    const bool queryAdvances = ov.queryIsRead0 ? (op == 3) : (op == 2);
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

    sort(scratch.svEvents.begin(), scratch.svEvents.end(),
        [](const PhasingSvEvent& a, const PhasingSvEvent& b) {
            return a.queryPos < b.queryPos;
        });

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

            for (uint64_t rid = start; rid < end; rid++) {
                const ReadId readId(rid);
                scratch.clear();

                // 1. Gather overlaps.
                gatherOverlaps(*this, readId, scratch);
                if (scratch.overlaps.empty()) {
                    readsProcessed++;
                    continue;
                }
                readsWithOverlaps++;

                // 2. Unpack query sequence.
                const uint32_t queryLen =
                    uint32_t(getReads().getRead(readId).baseCount);
                unpackQuerySequence(*this, readId, queryLen, scratch);

                // 3. Sliding-window SNP detection.
                detectSnpSites(*this, readId, queryLen, scratch);

                // 4. Build SNP matrix (filter + confirm sites).
                buildSnpMatrix(scratch);

                // 4b. Remove adjacent sites (positions p and p+1).
                filterAdjacentSites(scratch);

                if (!scratch.sites.empty()) {
                    readsWithSites++;

                    // 5. DP phasing.
                    // hetCov = number of overlaps (proxy for local coverage).
                    const uint32_t hetCov = uint32_t(scratch.overlaps.size());
                    runDpPhasing(scratch, hetCov);

                    // 6. Allele grouping (cis/trans labeling).
                    labelCisTrans(scratch);

                    // 7. Large indel phasing.
                    phaseLargeIndels(*this, readId, queryLen, scratch);
                }

                // 8. Dedup chains (best overlap per target).
                dedupChains(scratch);

                // 9. Write results back.
                writeResults(*this, readId, scratch);

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

    cout << timestamp << "Phasing complete." << endl;
    cout << timestamp << "Reads processed: " << readsProcessed.load() << endl;
    cout << timestamp << "Reads with overlaps: " << readsWithOverlaps.load() << endl;
    cout << timestamp << "Reads with SNP sites: " << readsWithSites.load() << endl;
    cout << timestamp << "Total cis labels: " << totalCis.load() << endl;
    cout << timestamp << "Total trans labels: " << totalTrans.load() << endl;
}
