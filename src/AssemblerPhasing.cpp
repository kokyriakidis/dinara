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
        //   flag=3 → HPC-masked candidate (evidence still collected)
        bool anyCandidates = false;
        for (uint32_t i = 0; i < wLen; i++) {
            if (scratch.flag[i] > PHASING_OCC_THRES) {
                const uint32_t absPos = wStart + i;
                if (isPeriodicRepeat(scratch.queryBases.data(), queryLen,
                                     absPos)) {
                    scratch.flag[i] = 3;
                } else {
                    scratch.flag[i] = 1;
                }
                anyCandidates = true;
            } else {
                scratch.flag[i] = 0;
            }
        }
        if (!anyCandidates) continue;

        // ---- Pass 2: evidence collection ----
        // At candidate positions (flag==1 or flag==3), emit PhasingEvidence
        // entries recording the observed base and whether it matches the
        // query. flag==3 positions are HPC-masked but evidence is still
        // collected (matching hifiasm); the isHpc flag is set on the
        // evidence but is unused in the ONT pipeline.
        for (size_t oi = 0; oi < scratch.overlaps.size(); oi++) {
            const auto& ov = scratch.overlaps[oi];
            if (ov.qe <= wStart || ov.qs >= wEnd) continue;

            for (uint32_t ei = scratch.cigarEventRanges[oi].begin;
                 ei < scratch.cigarEventRanges[oi].end; ei++) {
                const auto& ce = allEvents[ei];
                if (ce.qpos < wStart) continue;
                if (ce.qpos >= wEnd) break;

                uint32_t fi = ce.qpos - wStart;
                if (scratch.flag[fi] == 0) continue;

                PhasingEvidence ev;
                ev.site = ce.qpos;
                ev.overlapIdx = uint32_t(oi);
                ev.siteIdx = UINT32_MAX; // set properly in buildSnpMatrix
                ev.isHpc = (scratch.flag[fi] == 3) ? 1 : 0;

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
    //
    // Hifiasm creates separate SnpStats per qualifying alt base at each
    // position. All SnpStats at the same position share the same evidence
    // range. Each evidence entry has an overlapSite field pointing to the
    // SnpStats it belongs to:
    //   - Match evidence → overlapSite = LAST SnpStats at this position
    //   - Mismatch evidence → overlapSite = SnpStats for its specific base
    //   - Non-qualifying mismatches (base count < 2) → dropped
    //
    // We replicate this with ev.siteIdx on PhasingEvidence.

    const size_t rawEnd = scratch.evidence.size();
    size_t i = 0;
    while (i < rawEnd) {
        const uint32_t currentSite = scratch.evidence[i].site;
        const size_t groupBegin = i;

        // Find end of this site's group.
        while (i < rawEnd && scratch.evidence[i].site == currentSite) {
            i++;
        }
        const size_t groupEnd = i;

        // Count match and per-base alt.
        uint32_t matchCount = 0;
        uint32_t fwdStrandRefCount = 0;
        uint32_t altCount[4] = {0, 0, 0, 0};
        uint8_t isHpc = 0;
        uint8_t queryBase = scratch.queryBases[currentSite];

        for (size_t j = groupBegin; j < groupEnd; j++) {
            const auto& ev = scratch.evidence[j];
            if (ev.isAlt == 0) {
                matchCount++;
                if (scratch.overlaps[ev.overlapIdx].isRev == 0) {
                    fwdStrandRefCount++;
                }
            } else {
                if (ev.base < 4) altCount[ev.base]++;
            }
            if (ev.isHpc) isHpc = 1;
        }

        // +1 for query read (hifiasm: p->occ_0 = 1 + occ_0).
        matchCount += 1;

        // Reject position if ALL matching overlaps are forward-strand
        // (zero reverse-strand). Hifiasm: push_info rejects when
        // rev_n == occ_0 (Correct.cpp:10567). This is a strict subset
        // of the is_st_bs strand-bias check applied later in the DP.
        if (fwdStrandRefCount == matchCount - 1) continue;

        // Identify qualifying alt bases and create PhasingSite entries.
        // Record the siteIdx for each qualifying base.
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
            ps.evidenceBegin = 0; // set below
            ps.evidenceEnd = 0;

            lastSiteIdx = uint32_t(scratch.sites.size());
            baseSiteIdx[b] = lastSiteIdx;
            scratch.sites.push_back(ps);
        }

        if (lastSiteIdx == UINT32_MAX) continue; // no qualifying bases

        // Build shared evidence range: all matches + qualifying mismatches.
        // Set ev.siteIdx exactly like hifiasm sets overlapSite:
        //   match → lastSiteIdx (last SnpStats at this position)
        //   mismatch → baseSiteIdx[base] (SnpStats for its specific base)
        const uint32_t sharedBegin = uint32_t(scratch.evidence.size());

        for (size_t j = groupBegin; j < groupEnd; j++) {
            const auto& ev = scratch.evidence[j];
            if (ev.isAlt == 0) {
                PhasingEvidence out = ev;
                out.siteIdx = lastSiteIdx;
                scratch.evidence.push_back(out);
            } else if (ev.base < 4 && baseSiteIdx[ev.base] != UINT32_MAX) {
                PhasingEvidence out = ev;
                out.siteIdx = baseSiteIdx[ev.base];
                scratch.evidence.push_back(out);
            }
            // Non-qualifying mismatches → dropped
        }

        const uint32_t sharedEnd = uint32_t(scratch.evidence.size());

        // Set shared evidence range on all PhasingSite entries at this position.
        for (uint8_t b = 0; b < 4; b++) {
            if (baseSiteIdx[b] != UINT32_MAX) {
                scratch.sites[baseSiteIdx[b]].evidenceBegin = sharedBegin;
                scratch.sites[baseSiteIdx[b]].evidenceEnd = sharedEnd;
            }
        }
    }

    // Remove the raw evidence (indices 0..rawEnd). Shift tail to front
    // and update all evidenceBegin/evidenceEnd indices.
    if (!scratch.sites.empty()) {
        const size_t tailStart = rawEnd;
        const size_t tailLen = scratch.evidence.size() - tailStart;
        for (size_t j = 0; j < tailLen; j++) {
            scratch.evidence[j] = scratch.evidence[tailStart + j];
        }
        scratch.evidence.resize(tailLen);
        const uint32_t shift = uint32_t(tailStart);
        for (auto& ps : scratch.sites) {
            ps.evidenceBegin -= shift;
            ps.evidenceEnd -= shift;
        }
    }
}

// ============================================================================
// filterAdjacentSites: remove SNP sites at adjacent positions
//
// Port of hifiasm's adjacent-site filter (Correct.cpp:8855).
//
// Removes SNP sites at adjacent positions (p and p+1). Adjacent
// mismatches are typically alignment artifacts — an indel that the
// aligner represents as two consecutive substitutions — rather than
// real heterozygous variants. Keeping them would add noise to the
// phasing DP and labeling.
//
// With multi-alt sites, a single position can have multiple PhasingSite
// entries (one per qualifying alt base). When a position is removed,
// ALL its PhasingSite entries are removed.
//
// Three arrays must stay consistent:
//   sites[]    — compacted, indices shift down
//   evidence[] — entries for removed sites are dropped
//   ev.siteIdx — remapped from old site indices to new ones
//
// Hifiasm does the same compaction (Correct.cpp:8855): it compacts
// snp_stat[], compacts the evidence list[], and rewrites
// overlapSite -= m_off for each surviving evidence entry.
// ============================================================================

static void filterAdjacentSites(
    PhasingScratchpad& scratch)
{
    if (scratch.sites.size() < 2) return;

    const size_t n = scratch.sites.size();
    vector<bool> remove(n, false);

    // --- 1. Identify adjacent positions ---
    //
    // Group sites by position. Sites are sorted by position, and
    // multi-alt entries at the same position are contiguous.
    // If consecutive groups have positions p and p+1, mark both
    // groups for removal.

    struct PosGroup { uint32_t pos; size_t begin; size_t end; };
    vector<PosGroup> groups;
    size_t gi = 0;
    while (gi < n) {
        const uint32_t pos = scratch.sites[gi].site;
        size_t gEnd = gi;
        while (gEnd < n && scratch.sites[gEnd].site == pos) gEnd++;
        groups.push_back({pos, gi, gEnd});
        gi = gEnd;
    }

    for (size_t g = 0; g + 1 < groups.size(); g++) {
        if (groups[g + 1].pos == groups[g].pos + 1) {
            for (size_t i = groups[g].begin; i < groups[g].end; i++)
                remove[i] = true;
            for (size_t i = groups[g + 1].begin; i < groups[g + 1].end; i++)
                remove[i] = true;
        }
    }

    bool anyRemoved = false;
    for (size_t i = 0; i < n; i++) {
        if (remove[i]) { anyRemoved = true; break; }
    }
    if (!anyRemoved) return;

    // --- 2. Build old-to-new site index mapping ---
    //
    // Surviving sites shift down to fill gaps left by removed sites.
    // oldToNew[i] gives the new index for site i, or UINT32_MAX if
    // site i was removed. Used to remap ev.siteIdx in step 4.

    vector<uint32_t> oldToNew(n, UINT32_MAX);
    uint32_t newIdx = 0;
    for (size_t i = 0; i < n; i++) {
        if (!remove[i]) {
            oldToNew[i] = newIdx++;
        }
    }

    // --- 3. Compact sites array ---

    size_t write = 0;
    for (size_t i = 0; i < n; i++) {
        if (!remove[i]) {
            scratch.sites[write++] = scratch.sites[i];
        }
    }
    scratch.sites.resize(write);

    // --- 4. Compact evidence array and remap siteIdx ---
    //
    // Drop evidence entries whose site was removed (oldToNew == UINT32_MAX).
    // For surviving entries, remap ev.siteIdx to the new site index.
    // Evidence order is preserved (sorted by site position, then overlapIdx).

    for (auto& site : scratch.sites) {
        site.evidenceBegin = 0;
        site.evidenceEnd = 0;
    }

    size_t evWrite = 0;
    for (size_t ei = 0; ei < scratch.evidence.size(); ei++) {
        auto& ev = scratch.evidence[ei];

        if (ev.siteIdx < uint32_t(n)) {
            uint32_t newSi = oldToNew[ev.siteIdx];
            if (newSi == UINT32_MAX) continue;
            ev.siteIdx = newSi;
        }

        scratch.evidence[evWrite++] = ev;
    }
    scratch.evidence.resize(evWrite);

    // --- 5. Rebuild evidenceBegin/evidenceEnd on surviving sites ---
    //
    // Both sites and evidence are sorted by position, so we walk them
    // in tandem. All sites at the same position share the same evidence
    // range (set by buildSnpMatrix for multi-alt support).

    if (!scratch.evidence.empty()) {
        uint32_t si = 0;
        uint32_t ei = 0;
        while (si < uint32_t(scratch.sites.size()) &&
               ei < uint32_t(scratch.evidence.size())) {
            uint32_t pos = scratch.sites[si].site;
            uint32_t rangeBegin = ei;
            while (ei < uint32_t(scratch.evidence.size()) &&
                   scratch.evidence[ei].site == pos) {
                ei++;
            }
            while (si < uint32_t(scratch.sites.size()) &&
                   scratch.sites[si].site == pos) {
                scratch.sites[si].evidenceBegin = rangeBegin;
                scratch.sites[si].evidenceEnd = ei;
                si++;
            }
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
    int8_t chainId = 0;

    // Temporary: collect (chainLength, chainStartInBuffer) pairs.
    vector<pair<uint32_t, uint32_t>> chains;
    vector<uint32_t> chainBuf;

    for (const auto& [score, startIdx] : scratch.dpEndpoints) {
        if (scratch.dpChainId[startIdx] >= 0) continue; // already claimed

        uint32_t bufStart = uint32_t(chainBuf.size());
        int32_t k = int32_t(startIdx);
        while (k >= 0 && scratch.dpChainId[uint32_t(k)] < 0) {
            chainBuf.push_back(uint32_t(k));
            scratch.dpChainId[uint32_t(k)] = 127; // temporary claim
            k = scratch.dpParent[uint32_t(k)];
        }
        uint32_t chainLen = uint32_t(chainBuf.size()) - bufStart;
        if (chainLen > 0) {
            chains.push_back({chainLen, bufStart});
        }
    }

    // Reset temporary claims.
    scratch.dpChainId.assign(n, -1);

    // Sort chains by length descending (hifiasm does this).
    sort(chains.begin(), chains.end(),
        [](const auto& a, const auto& b) { return a.first > b.first; });

    // Assign chain IDs with confirmation logic.
    for (const auto& [chainLen, bufStart] : chains) {
        if (uint32_t(chainId) >= PHASING_MAX_DP_CHAINS) break;

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
    sort(overlapEvIdx.begin(), overlapEvIdx.end(),
        [&](uint32_t a, uint32_t b) {
            return scratch.evidence[a].overlapIdx < scratch.evidence[b].overlapIdx;
        });

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
        // match counts.
        for (uint32_t p = overlapEvBegin[sortIdx]; p < overlapEvEnd[sortIdx]; p++) {
            const auto& ev = scratch.evidence[overlapEvIdx[p]];
            if (ev.siteIdx >= numSites) continue;

            if (ev.isAlt == 1) {
                // Mismatch: mark as trans-confirmed.
                // Hifiasm: s->score = 1 for ALL mismatches at trans overlaps.
                scratch.sites[ev.siteIdx].transConfirmed = 1;
            } else {
                // Match: decrement labelMatchCount on ALL sites at this
                // position. Hifiasm walks backwards from overlapSite to
                // find all SnpStats at the same position.
                const uint32_t pos = scratch.sites[ev.siteIdx].site;
                // Walk backwards from ev.siteIdx.
                for (int32_t si = int32_t(ev.siteIdx); si >= 0; si--) {
                    if (scratch.sites[uint32_t(si)].site != pos) break;
                    auto& s = scratch.sites[uint32_t(si)];
                    if (s.labelMatchCount > 1) s.labelMatchCount--;
                    if (ov.isRev == 0 && s.labelFwdStrandCount > 1) {
                        s.labelFwdStrandCount--;
                    }
                }
                // Also walk forwards (ev.siteIdx may not be the last).
                for (uint32_t si = ev.siteIdx + 1; si < numSites; si++) {
                    if (scratch.sites[si].site != pos) break;
                    auto& s = scratch.sites[si];
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
            if (ev.isAlt == 1 && ev.siteIdx < numSites) {
                scratch.sites[ev.siteIdx].cisReset = 1;
            }
        }
    }

    // ----------------------------------------------------------------
    // Step E: multi_check — weak site promotion.
    //
    // Hifiasm (Correct.cpp:9035): for each cis overlap, find mismatches
    // at sites that pass basic thresholds but are NOT confirmed and NOT
    // trans-confirmed (after cisReset). Apply density + 32bp proximity
    // filter. Site indices appearing in >=2 overlaps get promoted.
    // ----------------------------------------------------------------

    {
        // Collect per-overlap weak mismatch siteIdx values.
        vector<uint32_t> promotionCandidates;

        for (uint32_t oi = 0; oi < numOv; oi++) {
            if (scratch.overlaps[oi].isMatch == 2) continue; // skip trans

            // Collect weak mismatch site indices for this overlap.
            vector<uint32_t> weakSiteIndices;
            for (uint32_t p = overlapEvBegin[oi]; p < overlapEvEnd[oi]; p++) {
                const auto& ev = scratch.evidence[overlapEvIdx[p]];
                if (ev.isAlt == 0) continue;
                if (ev.siteIdx >= numSites) continue;
                const auto& s = scratch.sites[ev.siteIdx];
                if (s.labelMatchCount < 2 || s.altCount < 2) continue;
                if (isLabelStrandBiased(s)) continue;
                // Skip confirmed sites (pass full thresholds).
                if (s.labelMatchCount >= PHASING_S_HAP_COV &&
                    s.altCount >= PHASING_INFOR_COV) continue;
                // Skip already confirmed (trans-confirmed and not cis-reset).
                if (s.isLabelConfirmed()) continue;
                weakSiteIndices.push_back(ev.siteIdx);
            }

            uint32_t o = uint32_t(weakSiteIndices.size());
            if (o == 0) continue;

            // Density filter: need >= 4% of alignment length.
            const uint32_t alignLen = scratch.overlaps[oi].qe -
                                      scratch.overlaps[oi].qs;
            if (o < uint32_t(double(alignLen) * 0.04)) continue;

            // Sort by siteIdx (hifiasm sorts overlapSite values).
            sort(weakSiteIndices.begin(), weakSiteIndices.end());

            // 32bp proximity filter: remove sites within 32bp of neighbors.
            // Hifiasm uses site positions looked up via snp_stat[a[i]].site.
            //
            // Note: hifiasm has a stale-pointer bug here (Correct.cpp:9068):
            // the `t` pointer for the last element retains its value from
            // the previous iteration, pointing to the current element
            // itself, so `site + 32 > site` is always true and the last
            // element is always dropped. We replicate this for parity.
            vector<uint32_t> filtered;
            for (size_t i = 0; i < weakSiteIndices.size(); i++) {
                uint32_t pos = scratch.sites[weakSiteIndices[i]].site;
                bool tooClose = false;
                if (i > 0) {
                    uint32_t prevPos = scratch.sites[weakSiteIndices[i-1]].site;
                    if (prevPos + 32 > pos) tooClose = true;
                }
                // Last element: hifiasm always filters it (stale t pointer).
                if (i + 1 < weakSiteIndices.size()) {
                    uint32_t nextPos = scratch.sites[weakSiteIndices[i+1]].site;
                    if (pos + 32 > nextPos) tooClose = true;
                } else {
                    tooClose = true; // replicate hifiasm last-element drop
                }
                if (!tooClose) filtered.push_back(weakSiteIndices[i]);
            }

            if (filtered.size() >= 2) {
                for (uint32_t si : filtered) {
                    promotionCandidates.push_back(si);
                }
            }
        }

        // Promote: siteIdx values appearing in >=2 overlaps.
        if (!promotionCandidates.empty()) {
            sort(promotionCandidates.begin(), promotionCandidates.end());
            for (size_t i = 0; i < promotionCandidates.size(); ) {
                uint32_t si = promotionCandidates[i];
                size_t count = 0;
                while (i < promotionCandidates.size() &&
                       promotionCandidates[i] == si) {
                    count++;
                    i++;
                }
                if (count >= 2) {
                    scratch.sites[si].promoted = 1;
                }
            }
        }
    }

    // ----------------------------------------------------------------
    // Step F: Final loop — set strong flag, apply promotions.
    //
    // Hifiasm (Correct.cpp:9085): for each overlap:
    //   - trans: strong=1
    //   - cis: check evidence at confirmed sites (transConfirmed and not
    //     cisReset, OR promoted). If mismatch found, flip to trans.
    //     If match found, set strong=1.
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
                if (!s.isLabelConfirmed()) continue;
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
                    // Hifiasm: cc = max(6, (het_cov > 0 ? het_cov : hom_cov/n_hap) * 0.7)
                    // het_cov is typically 0; hom_cov = coveragePeak; n_hap = 2.
                    const uint64_t coveragePeak =
                        assemblerInfo->kmerDistributionInfo.coveragePeak;
                    const uint32_t hetCov = uint32_t(coveragePeak / 2);
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
