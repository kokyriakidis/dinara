/// @file AssemblerPhasingKmeans.cpp
/// @brief K-means overlap phasing adapted from pgphase/longcallD.

#include "Assembler.hpp"
#include "PhasingKmeansTypes.hpp"
#include "Alignment.hpp"
#include "Reads.hpp"
#include "invalid.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <numeric>
#include <thread>
#include <vector>

using namespace std;
using namespace dinara;

// ============================================================================
// Helpers
// ============================================================================

static uint64_t kmCigarRead1Start(
    const Assembler& assembler, const AlignmentData& ad)
{
    if (ad.isSameStrand) return ad.ts;
    const uint32_t targetLen = uint32_t(
        assembler.getReads().getRead(ad.readIds[1]).baseCount);
    return targetLen - ad.te;
}

static inline uint8_t kmGetBase(
    const Assembler& assembler, ReadId readId,
    uint32_t position, bool isReverseComplement)
{
    const auto sequence = assembler.getReads().getRead(readId);
    if (!isReverseComplement) return sequence[position].value;
    return sequence[sequence.baseCount - 1 - position].complement().value;
}

static constexpr uint32_t KM_HPC_RR = 4;
static constexpr uint32_t KM_HPC_CC = 2;
static constexpr uint32_t KM_HPC_PL = 12;

static bool kmIsPeriodicRepeat(const uint8_t* seq, uint32_t seqLen, uint32_t p)
{
    if (seqLen == 0) return false;
    const int64_t sn = int64_t(seqLen), pp = int64_t(p);
    const int64_t s = (pp >= int64_t(KM_HPC_PL)) ? (pp - int64_t(KM_HPC_PL)) : 0;
    const int64_t e = ((pp + int64_t(KM_HPC_PL)) <= sn) ? (pp + int64_t(KM_HPC_PL)) : sn;
    for (uint32_t r = 1; r <= KM_HPC_RR; r++) {
        const int64_t rc = int64_t(r) * int64_t(KM_HPC_CC);
        { int64_t k=pp+int64_t(r); while(k<e&&(k-int64_t(r))>=s&&seq[k]==seq[k-r])k++;
          int64_t ze=k; if(ze>e)ze=e; k=pp-1;
          while(k>=s&&(k+int64_t(r))<e&&seq[k]==seq[k+r])k--;
          int64_t zs=k+1; if(zs<s)zs=s;
          if((ze-zs)>int64_t(r)&&(ze-zs)>=rc) return true; }
        { int64_t k=pp+int64_t(r)+1; while(k<e&&(k-int64_t(r))>=s&&seq[k]==seq[k-r])k++;
          int64_t zs=pp+1; if(zs<s)zs=s; int64_t ze=k; if(ze>e)ze=e;
          if((ze-zs)>int64_t(r)&&(ze-zs)>=rc) return true; }
        { int64_t k=pp-int64_t(r); while(k>=s&&(k+int64_t(r))<e&&seq[k]==seq[k+r])k--;
          int64_t zs=k+1; if(zs<s)zs=s; k=pp+1;
          while(k<e&&(k-int64_t(r))>=s&&seq[k]==seq[k-r])k++;
          int64_t ze=k; if(ze>e)ze=e;
          if((ze-zs)>int64_t(r)&&(ze-zs)>=rc) return true; }
        { int64_t k=pp-int64_t(r)-1; while(k>=s&&(k+int64_t(r))<e&&seq[k]==seq[k+r])k--;
          int64_t zs=k+1; if(zs<s)zs=s; int64_t ze=pp; if(ze>e)ze=e;
          if((ze-zs)>int64_t(r)&&(ze-zs)>=rc) return true; }
    }
    return false;
}

// ============================================================================
// Step 1: Gather overlaps
// ============================================================================

static void kmGatherOverlaps(
    const Assembler& assembler, ReadId backboneReadId, KmScratchpad& scratch)
{
    const auto& alignmentTable = assembler.getAlignmentTable();
    const auto& alignmentData = assembler.alignmentData;
    const OrientedReadId orientedId(backboneReadId, 0);
    if (orientedId.getValue() >= alignmentTable.size()) return;
    const auto alignmentIds = alignmentTable[orientedId.getValue()];

    for (const uint32_t alignmentId : alignmentIds) {
        const AlignmentData& ad = alignmentData[alignmentId];
        const AlignmentInfo& info = ad.info;
        if (ad.isDeleted()) continue;
        if (info.cigarOffset == uint32_t(-1) || info.cigarTokenCount == 0) continue;

        KmOverlap ov;
        ov.alignmentId = alignmentId;
        ov.cigarOffset = info.cigarOffset;
        ov.cigarTokenCount = info.cigarTokenCount;
        ov.strong = 0; ov.hap = 0; ov.phaseSet = -1;
        if (ad.readIds[0] == backboneReadId) {
            ov.queryIsRead0 = 1; ov.targetReadId = ad.readIds[1];
            ov.qs = ad.qs; ov.qe = ad.qe; ov.ts = ad.ts; ov.te = ad.te;
            ov.isRev = ad.isSameStrand ? 0 : 1;
        } else {
            ov.queryIsRead0 = 0; ov.targetReadId = ad.readIds[0];
            ov.qs = ad.ts; ov.qe = ad.te; ov.ts = ad.qs; ov.te = ad.qe;
            ov.isRev = ad.isSameStrand ? 0 : 1;
        }
        scratch.overlaps.push_back(ov);
    }
}

// ============================================================================
// Step 2: Parse CIGARs into per-overlap digars (single walk, like pgphase)
// ============================================================================

static void kmParseCigars(
    const Assembler& assembler, ReadId backboneReadId,
    uint32_t backboneLen, KmScratchpad& scratch)
{
    const auto& cigarStore = assembler.getOverlapCigarStore();
    const uint32_t numOv = uint32_t(scratch.overlaps.size());

    scratch.digars.clear();
    scratch.digarBegin.resize(numOv);
    scratch.digarEnd.resize(numOv);

    for (uint32_t oi = 0; oi < numOv; oi++) {
        scratch.digarBegin[oi] = uint32_t(scratch.digars.size());
        const auto& ov = scratch.overlaps[oi];
        const auto& ad = assembler.alignmentData[ov.alignmentId];
        const bool needsRc = (ov.queryIsRead0 == 0) && (ov.isRev != 0);
        const bool qIsR0 = (ov.queryIsRead0 != 0);

        cigarStore.forEachOpWithPositions(
            ov.cigarOffset, ov.cigarTokenCount,
            ad.qs, kmCigarRead1Start(assembler, ad),
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (op == 1) { // Mismatch — one digar per base with actual alt base
                    for (uint32_t b = 0; b < len; b++) {
                        uint32_t bbPos, tPos;
                        if (needsRc) {
                            uint32_t rc = uint32_t(yk) + b;
                            if (rc >= backboneLen) continue;
                            bbPos = backboneLen - 1 - rc;
                            tPos = uint32_t(xk) + b;
                        } else if (qIsR0) {
                            bbPos = uint32_t(xk) + b;
                            tPos = uint32_t(yk) + b;
                        } else {
                            bbPos = uint32_t(yk) + b;
                            tPos = uint32_t(xk) + b;
                        }
                        if (bbPos >= backboneLen) continue;
                        // Read actual alt base from target read.
                        uint8_t altBase = kmGetBase(assembler, ReadId(ov.targetReadId),
                                                    tPos, ov.isRev != 0);
                        scratch.digars.push_back({bbPos, KmVarType::Snp, altBase, 1});
                    }
                } else if (op == 2 || op == 3) { // Ins/Del
                    // Determine if this consumes backbone or target.
                    bool bbConsumed = (op==3 && qIsR0) || (op==2 && !qIsR0);
                    if (bbConsumed) {
                        // Deletion on backbone (or insertion on target that consumes backbone).
                        uint32_t raw = (op==3) ? uint32_t(xk) : uint32_t(yk);
                        uint32_t rawE = (raw + len < backboneLen) ? raw + len : backboneLen;
                        uint32_t s, e;
                        if (needsRc) { s = backboneLen - rawE; e = backboneLen - raw; }
                        else { s = raw; e = rawE; }
                        // One digar for the whole deletion event, anchored at start.
                        if (s < backboneLen)
                            scratch.digars.push_back({s, KmVarType::Deletion, 0, uint16_t(e - s)});
                    } else {
                        // Insertion on backbone (target has extra bases).
                        uint32_t anchor = qIsR0 ? uint32_t(xk) : uint32_t(yk);
                        if (needsRc && anchor > 0) anchor = backboneLen - 1 - anchor;
                        if (anchor < backboneLen)
                            scratch.digars.push_back({anchor, KmVarType::Insertion, 0, uint16_t(len)});
                    }
                }
            });

        // Sort this overlap's digars by backbone position.
        uint32_t dEnd = uint32_t(scratch.digars.size());
        sort(scratch.digars.begin() + scratch.digarBegin[oi],
             scratch.digars.begin() + dEnd);
        scratch.digarEnd[oi] = dEnd;
    }
}

// ============================================================================
// Step 3a: Collect candidate sites from digars (sort + dedup like pgphase)
// ============================================================================

/// Fuzzy insertion comparison (mirrors pgphase exact_comp_var_site_ins).
/// For insertions >= minSvLen at the same position: returns 0 if
/// min(len1,len2) >= max(len1,len2) * 0.8 (within 20% length ratio).
/// For shorter insertions and all other types: exact match required.
static inline int kmCompareKeysFuzzy(const KmVarKey& a, const KmVarKey& b, int minSvLen) {
    uint32_t asp = (a.type == KmVarType::Snp) ? a.pos : (a.pos > 0 ? a.pos - 1 : 0);
    uint32_t bsp = (b.type == KmVarType::Snp) ? b.pos : (b.pos > 0 ? b.pos - 1 : 0);
    if (asp != bsp) return asp < bsp ? -1 : 1;
    if (a.type != b.type) return a.type < b.type ? -1 : 1;
    if (a.refLen != b.refLen) return a.refLen < b.refLen ? -1 : 1;
    if (a.type == KmVarType::Snp) {
        if (a.altLen != b.altLen) return a.altLen < b.altLen ? -1 : 1;
        if (a.altBase != b.altBase) return a.altBase < b.altBase ? -1 : 1;
        return 0;
    }
    if (a.type == KmVarType::Insertion) {
        if (int(a.altLen) < minSvLen) {
            // Short insertion: exact length + base match required.
            if (a.altLen != b.altLen) return a.altLen < b.altLen ? -1 : 1;
            // We don't store insertion sequence, so same length = same.
            return 0;
        }
        // Large insertion: fuzzy length-ratio rule.
        int minL = int(a.altLen) < int(b.altLen) ? int(a.altLen) : int(b.altLen);
        int maxL = int(a.altLen) > int(b.altLen) ? int(a.altLen) : int(b.altLen);
        if (double(minL) >= double(maxL) * 0.8) return 0;
        return int(a.altLen) - int(b.altLen);
    }
    // Deletion: refLen already compared above, that's sufficient.
    return 0;
}

static void kmCollectCandidates(KmScratchpad& scratch, int minSvLen)
{
    // Emit one candidate key per digar event.
    vector<KmVarKey> keys;
    keys.reserve(scratch.digars.size());
    for (const auto& d : scratch.digars) {
        KmVarKey k;
        k.pos = d.pos;
        k.type = d.type;
        k.altBase = (d.type == KmVarType::Snp) ? d.altBase : 0;
        k.refLen = (d.type == KmVarType::Deletion) ? d.len : (d.type == KmVarType::Snp ? 1 : 0);
        k.altLen = (d.type == KmVarType::Insertion) ? d.len : (d.type == KmVarType::Snp ? 1 : 0);
        keys.push_back(k);
    }

    // Sort by exact order.
    sort(keys.begin(), keys.end());

    // Collapse fuzzy large insertions (mirrors pgphase collapse_fuzzy_large_insertions).
    // Then exact dedup for everything else.
    if (!keys.empty()) {
        size_t wi = 1;
        for (size_t ri = 1; ri < keys.size(); ri++) {
            if (kmCompareKeysFuzzy(keys[wi - 1], keys[ri], minSvLen) == 0)
                continue; // merge-equivalent → skip
            if (wi != ri) keys[wi] = keys[ri];
            wi++;
        }
        keys.resize(wi);
    }

    scratch.candidates.clear();
    scratch.candidates.reserve(keys.size());
    for (const auto& k : keys) {
        KmCandidate cand;
        cand.key = k;
        scratch.candidates.push_back(cand);
    }
}

// ============================================================================
// Step 3b: Allele counting via merge walk
// ============================================================================

/// Compare digar to candidate key using fuzzy insertion matching.
/// Used in both allele counting and profile merge walks.
/// Mirrors pgphase lcd_exact_comp_var_site_ins_merge.
static inline int kmCompareDigarToCand(const KmDigarOp& d, const KmVarKey& k, int minSvLen) {
    uint32_t dsp = (d.type == KmVarType::Snp) ? d.pos : (d.pos > 0 ? d.pos - 1 : 0);
    uint32_t ksp = (k.type == KmVarType::Snp) ? k.pos : (k.pos > 0 ? k.pos - 1 : 0);
    if (dsp != ksp) return dsp < ksp ? -1 : 1;
    if (d.type != k.type) return d.type < k.type ? -1 : 1;
    // refLen comparison.
    uint16_t dRefLen = (d.type == KmVarType::Deletion) ? d.len : (d.type == KmVarType::Snp ? 1 : 0);
    if (dRefLen != k.refLen) return dRefLen < k.refLen ? -1 : 1;
    if (d.type == KmVarType::Snp) {
        uint16_t dAltLen = 1;
        if (dAltLen != k.altLen) return dAltLen < k.altLen ? -1 : 1;
        if (d.altBase != k.altBase) return d.altBase < k.altBase ? -1 : 1;
        return 0;
    }
    if (d.type == KmVarType::Insertion) {
        if (int(d.len) < minSvLen) {
            if (d.len != k.altLen) return d.len < k.altLen ? -1 : 1;
            return 0;
        }
        int minL = int(d.len) < int(k.altLen) ? int(d.len) : int(k.altLen);
        int maxL = int(d.len) > int(k.altLen) ? int(d.len) : int(k.altLen);
        if (double(minL) >= double(maxL) * 0.8) return 0;
        return int(d.len) - int(k.altLen);
    }
    return 0;
}

// ============================================================================
// Step 3b: Allele counting only (merge walk, like pgphase collect_allele_counts_from_records)
// ============================================================================

static void kmCountAlleles(KmScratchpad& scratch, int minSvLen)
{
    const uint32_t numCand = uint32_t(scratch.candidates.size());
    const uint32_t numOv = uint32_t(scratch.overlaps.size());
    if (numCand == 0 || numOv == 0) return;

    for (uint32_t ci = 0; ci < numCand; ci++) {
        auto& c = scratch.candidates[ci];
        c.totalCov = 1; c.refCov = 1; c.fwdRef = 1;
    }

    for (uint32_t oi = 0; oi < numOv; oi++) {
        const auto& ov = scratch.overlaps[oi];
        uint32_t di = scratch.digarBegin[oi];
        uint32_t diEnd = scratch.digarEnd[oi];

        uint32_t ciStart = uint32_t(lower_bound(
            scratch.candidates.begin(), scratch.candidates.end(), ov.qs,
            [](const KmCandidate& c, uint32_t p) { return c.key.pos < p; })
            - scratch.candidates.begin());

        for (uint32_t ci = ciStart; ci < numCand && scratch.candidates[ci].key.pos < ov.qe; ci++) {
            auto& c = scratch.candidates[ci];
            c.totalCov++;

            while (di < diEnd && kmCompareDigarToCand(scratch.digars[di], c.key, minSvLen) < 0)
                di++;

            bool isAlt = false;
            if (di < diEnd && kmCompareDigarToCand(scratch.digars[di], c.key, minSvLen) == 0) {
                isAlt = true;
                di++;
            }

            if (isAlt) {
                c.altCov++;
                if (ov.isRev == 0) c.fwdAlt++; else c.revAlt++;
            } else {
                c.refCov++;
                if (ov.isRev == 0) c.fwdRef++; else c.revRef++;
            }
        }
    }

    for (uint32_t ci = 0; ci < numCand; ci++) {
        auto& c = scratch.candidates[ci];
        c.alleleFraction = c.totalCov > 0 ? double(c.altCov) / double(c.totalCov) : 0.0;
    }
}

// ============================================================================
// Step 5: Build per-overlap allele profiles (after classification)
//
// Mirrors pgphase collect_read_var_profile:
// - Skips NON_VAR candidates
// - Detects overlapping variants (deletion spanning a SNP → allele = -1)
// - Separate merge walk per overlap against candidate table
// ============================================================================

/// Check if digar and candidate overlap positionally (pgphase lcd_profile_ovlp_var_site).
/// Returns true if the reference spans of the two events overlap.
/// Positional overlap check (mirrors pgphase lcd_profile_ovlp_var_site).
/// Returns true if the reference spans of digar and candidate overlap.
static inline bool kmDigarOverlapsCand(const KmDigarOp& d, const KmVarKey& k) {
    // Compute reference spans.
    // SNP: refLen=1. Deletion: refLen=len. Insertion: refLen=0.
    uint32_t dRefLen = (d.type == KmVarType::Snp) ? 1 :
                       (d.type == KmVarType::Deletion) ? d.len : 0;
    uint32_t dBeg = d.pos;
    uint32_t dEnd = d.pos + dRefLen;
    uint32_t kBeg = k.pos;
    uint32_t kEnd = k.pos + k.refLen;

    // Both insertions (refLen=0): overlap iff same position.
    if (dRefLen == 0 && k.refLen == 0)
        return dBeg == kBeg;
    // Digar is insertion: overlaps if strictly inside candidate's span.
    if (dRefLen == 0)
        return dBeg > kBeg && dBeg < kEnd;
    // Candidate is insertion: overlaps if strictly inside digar's span.
    if (k.refLen == 0)
        return kBeg > dBeg && kBeg < dEnd;
    // Both have nonzero ref span — standard interval overlap.
    return !(dBeg >= kEnd || kBeg >= dEnd);
}

static void kmBuildOverlapProfiles(KmScratchpad& scratch, int minSvLen)
{
    const uint32_t numCand = uint32_t(scratch.candidates.size());
    const uint32_t numOv = uint32_t(scratch.overlaps.size());
    if (numCand == 0 || numOv == 0) return;

    scratch.overlapProfiles.resize(numOv);

    for (uint32_t oi = 0; oi < numOv; oi++) {
        const auto& ov = scratch.overlaps[oi];
        auto& prof = scratch.overlapProfiles[oi];
        prof.overlapIdx = oi;
        prof.startVarIdx = -1; prof.endVarIdx = -1;
        prof.alleles.clear();

        uint32_t di = scratch.digarBegin[oi];
        uint32_t diEnd = scratch.digarEnd[oi];

        uint32_t ciStart = uint32_t(lower_bound(
            scratch.candidates.begin(), scratch.candidates.end(), ov.qs,
            [](const KmCandidate& c, uint32_t p) { return c.key.pos < p; })
            - scratch.candidates.begin());

        for (uint32_t ci = ciStart; ci < numCand && scratch.candidates[ci].key.pos < ov.qe; ci++) {
            const auto& c = scratch.candidates[ci];

            // Skip NON_VAR candidates (like pgphase skips kLongcalldNonVar).
            if (c.categoryFlag == KM_NON_VAR) continue;

            // Advance digar pointer, skipping digars before this candidate.
            while (di < diEnd && kmCompareDigarToCand(scratch.digars[di], c.key, minSvLen) < 0
                   && !kmDigarOverlapsCand(scratch.digars[di], c.key))
                di++;

            // Determine allele observation.
            int allele;
            if (di < diEnd && kmCompareDigarToCand(scratch.digars[di], c.key, minSvLen) == 0) {
                // Exact match → alt.
                allele = 1;
                di++;
            } else if (di < diEnd && kmDigarOverlapsCand(scratch.digars[di], c.key)) {
                // Digar overlaps candidate but doesn't match exactly.
                // E.g., deletion spanning a SNP position → no observation.
                allele = -1;
            } else {
                // No digar at this position → ref.
                allele = 0;
            }

            // Append to profile.
            if (prof.startVarIdx < 0) prof.startVarIdx = int(ci);
            while (prof.startVarIdx + int(prof.alleles.size()) < int(ci))
                prof.alleles.push_back(-1);
            prof.endVarIdx = int(ci);
            prof.alleles.push_back(allele);
        }

        // Tail: remaining candidates within overlap range that aren't NON_VAR get ref.
        // (Already handled by the for loop — they have no matching digar → allele = 0.)
    }
}

// ============================================================================
// Step 4: Classify candidates and detect noisy regions
// ============================================================================

static void kmClassifyCandidates(
    uint32_t backboneLen, KmScratchpad& scratch, const KmPhasingOptions& opts)
{
    const uint32_t numCand = uint32_t(scratch.candidates.size());
    if (numCand == 0) return;

    for (uint32_t ci = 0; ci < numCand; ci++) {
        auto& c = scratch.candidates[ci];
        if (c.totalCov < int(opts.minDepth) || c.altCov < int(opts.minAltDepth)) {
            c.category = KmVariantCategory::LowCoverage; c.categoryFlag = 0; continue;
        }
        { int fa = c.fwdAlt, ra = c.revAlt;
          if ((fa == 0 || ra == 0) && (fa + ra) >= 4) {
              c.category = KmVariantCategory::StrandBias; c.categoryFlag = 0; continue; } }
        if (c.alleleFraction < opts.minAf) {
            c.category = KmVariantCategory::LowAlleleFraction; c.categoryFlag = 0; continue;
        }
        if (c.alleleFraction > opts.maxAf) {
            c.category = KmVariantCategory::CleanHom; c.categoryFlag = KM_CLEAN_HOM; continue;
        }
        bool isIndel = (c.key.type != KmVarType::Snp);
        if (isIndel) {
            if (kmIsPeriodicRepeat(scratch.backboneBases.data(), backboneLen, c.key.pos)) {
                c.category = KmVariantCategory::RepeatHetIndel; c.categoryFlag = KM_REP_HET_VAR; continue;
            }
            c.category = KmVariantCategory::CleanHetIndel; c.categoryFlag = KM_CLEAN_HET_INDEL; continue;
        }
        if (kmIsPeriodicRepeat(scratch.backboneBases.data(), backboneLen, c.key.pos)) {
            c.category = KmVariantCategory::RepeatHetIndel; c.categoryFlag = KM_REP_HET_VAR; continue;
        }
        c.category = KmVariantCategory::CleanHetSnp; c.categoryFlag = KM_CLEAN_HET_SNP;
    }

    // Noisy region seeds.
    scratch.noisyRegions.clear();
    vector<pair<uint32_t,uint32_t>> seeds;
    for (uint32_t ci = 0; ci < numCand; ci++) {
        if (scratch.candidates[ci].category == KmVariantCategory::RepeatHetIndel)
            seeds.push_back({scratch.candidates[ci].key.pos,
                             scratch.candidates[ci].key.pos + scratch.candidates[ci].key.refLen});
    }
    for (uint32_t ci = 1; ci < numCand; ci++) {
        const auto& prev = scratch.candidates[ci-1];
        const auto& cur = scratch.candidates[ci];
        if (prev.category == KmVariantCategory::LowCoverage || prev.category == KmVariantCategory::NonVariant) continue;
        if (cur.category == KmVariantCategory::LowCoverage || cur.category == KmVariantCategory::NonVariant) continue;
        if (cur.key.pos - prev.key.pos <= 25) {
            if (prev.category != KmVariantCategory::CleanHetSnp || cur.category != KmVariantCategory::CleanHetSnp)
                seeds.push_back({prev.key.pos, cur.key.pos + cur.key.refLen});
        }
    }
    if (seeds.empty()) return;

    sort(seeds.begin(), seeds.end());
    vector<pair<uint32_t,uint32_t>> merged;
    merged.push_back(seeds[0]);
    for (size_t i = 1; i < seeds.size(); i++) {
        auto& last = merged.back();
        if (seeds[i].first <= last.second + opts.noisyRegMergeDis)
            last.second = max(last.second, seeds[i].second);
        else merged.push_back(seeds[i]);
    }
    for (const auto& [s, e] : merged)
        scratch.noisyRegions.push_back({s, e, 0, false});

    size_t ri = 0;
    for (uint32_t ci = 0; ci < numCand && ri < scratch.noisyRegions.size(); ci++) {
        auto& c = scratch.candidates[ci];
        while (ri < scratch.noisyRegions.size() && scratch.noisyRegions[ri].end <= c.key.pos) ri++;
        if (ri < scratch.noisyRegions.size() &&
            c.key.pos >= scratch.noisyRegions[ri].start && c.key.pos < scratch.noisyRegions[ri].end) {
            c.category = KmVariantCategory::NonVariant; c.categoryFlag = KM_NON_VAR;
        }
    }
}

// ============================================================================
// Step 6: K-means phasing — helpers
// ============================================================================

static inline int kmGetAllele(const KmScratchpad& s, uint32_t oi, uint32_t ci) {
    const auto& p = s.overlapProfiles[oi];
    if (p.startVarIdx < 0) return -1;
    int off = int(ci) - p.startVarIdx;
    if (off < 0 || off >= int(p.alleles.size())) return -1;
    return p.alleles[off];
}

static inline int kmAlleSlots(const KmCandidate& c) {
    return c.nUniqAlles > 0 ? c.nUniqAlles : 2;
}

static int kmSelectPivot(const KmScratchpad& s, const vector<uint32_t>& vi) {
    int bSnp = -1, bIndel = -1; int bSnpD = 0, bIndelD = 0;
    for (int i = 0; i < int(vi.size()); i++) {
        const auto& c = s.candidates[vi[i]];
        if (c.category == KmVariantCategory::CleanHetSnp && (bSnp < 0 || c.totalCov > bSnpD))
            { bSnp = i; bSnpD = c.totalCov; }
        else if (c.category == KmVariantCategory::CleanHetIndel && (bIndel < 0 || c.totalCov > bIndelD))
            { bIndel = i; bIndelD = c.totalCov; }
    }
    return bSnp >= 0 ? bSnp : bIndel;
}

static int kmReadToConsScore(KmCandidate& c, int hap, int alle) {
    int vs = (c.category == KmVariantCategory::CleanHetSnp ||
              c.category == KmVariantCategory::CleanHetIndel) ? 2 : 1;
    if (c.hapConsAlle[hap] == -1 && c.hapConsAlle[3-hap] == -1) return 0;
    if (c.hapConsAlle[hap] == -1) c.hapConsAlle[hap] = 1 - c.hapConsAlle[3-hap];
    if (c.hapConsAlle[3-hap] == -1) c.hapConsAlle[3-hap] = 1 - c.hapConsAlle[hap];
    if (c.hapConsAlle[hap] == alle) return vs;
    if (c.hapConsAlle[hap] == -1) return 0;
    return -vs;
}

static void kmUpdateConsAlle(KmCandidate& c, int hap) {
    if (hap == 0) return;
    const auto& prof = c.hapAlleProfile[hap];
    int mx = 0, ma = -1;
    for (int a = 0; a < int(prof.size()); a++)
        if (prof[a] > mx) { mx = prof[a]; ma = a; }
    c.hapConsAlle[hap] = ma;
}

static int kmAssignOverlapHap(KmScratchpad& s, uint32_t oi, uint32_t flags) {
    const auto& prof = s.overlapProfiles[oi];
    if (prof.startVarIdx < 0) return -1;
    int hs[3] = {0,0,0}, nv[3] = {0,0,0};
    for (int ci = prof.startVarIdx; ci <= prof.endVarIdx; ci++) {
        auto& c = s.candidates[ci];
        if ((c.categoryFlag & flags) == 0) continue;
        if (c.isHomopolymerIndel) continue;
        int a = kmGetAllele(s, oi, uint32_t(ci));
        if (a < 0) continue;
        for (int h = 1; h <= 2; h++) {
            int sc = kmReadToConsScore(c, h, a);
            if (sc != 0) nv[h]++;
            if (c.category != KmVariantCategory::CleanHom) hs[h] += sc;
        }
    }
    if (nv[1] == 0 && nv[2] == 0) return -1;
    int mxH = 0, mxS = 0, mnH = 0, mnS = 0;
    for (int h = 1; h <= 2; h++) {
        if (hs[h] > mxS) { mxH = h; mxS = hs[h]; }
        else if (hs[h] < mnS) { mnH = h; mnS = hs[h]; }
    }
    if (mxS == 0 && mnS == 0) return 0;
    return mxS > 0 ? mxH : 3 - mnH;
}

static void kmUpdateProfiles(KmScratchpad& s, uint32_t oi, int hap, uint32_t flags) {
    const auto& prof = s.overlapProfiles[oi];
    if (prof.startVarIdx < 0) return;
    for (int ci = prof.startVarIdx; ci <= prof.endVarIdx; ci++) {
        auto& c = s.candidates[ci];
        if ((c.categoryFlag & flags) == 0) continue;
        int a = kmGetAllele(s, oi, uint32_t(ci));
        if (a < 0 || a >= int(c.hapAlleProfile[1].size())) continue;
        if (hap == 0) {
            c.hapAlleProfile[1][a]++; c.hapAlleProfile[2][a]++;
            kmUpdateConsAlle(c, 1); kmUpdateConsAlle(c, 2);
        } else {
            c.hapAlleProfile[hap][a]++;
            kmUpdateConsAlle(c, hap);
        }
    }
}

// ============================================================================
// Step 6b: Phase-set flip detection
// ============================================================================

static bool kmPhaseSetFlip(KmScratchpad& scratch,
    const vector<uint32_t>& validIdx, const KmPhasingOptions& opts)
{
    const int n = int(validIdx.size());
    vector<int> hetIdx;
    for (int vi = 0; vi < n; vi++) {
        const auto& c = scratch.candidates[validIdx[vi]];
        if (c.hapConsAlle[1] != -1 && c.hapConsAlle[2] != -1 &&
            c.hapConsAlle[1] != c.hapConsAlle[2] && !c.isHomopolymerIndel)
            hetIdx.push_back(vi);
    }

    vector<int> nAgree(n, 0), nConflict(n, 0);
    for (int hi = 1; hi < int(hetIdx.size()); hi++) {
        int vi = hetIdx[hi];
        uint32_t ci = validIdx[vi], prevCi = validIdx[hetIdx[hi-1]];
        const auto& cur = scratch.candidates[ci];
        const auto& prev = scratch.candidates[prevCi];
        for (uint32_t oi = 0; oi < uint32_t(scratch.overlaps.size()); oi++) {
            int hap = scratch.overlaps[oi].hap;
            if (hap == 0) continue;
            int a1 = kmGetAllele(scratch, oi, prevCi);
            int a2 = kmGetAllele(scratch, oi, ci);
            if (a1 < 0 || a2 < 0) continue;
            if (prev.hapConsAlle[hap] == a1 && cur.hapConsAlle[hap] == a2) nAgree[vi]++;
            else if (prev.hapConsAlle[hap] == a1 && cur.hapConsAlle[3-hap] == a2) nConflict[vi]++;
        }
    }

    bool changed = false;
    int flip = 0;
    int32_t ps = int32_t(scratch.candidates[validIdx[0]].key.pos);
    scratch.candidates[validIdx[0]].phaseSet = ps;
    for (int vi = 1; vi < n; vi++) {
        auto& c = scratch.candidates[validIdx[vi]];
        bool isHet = (c.hapConsAlle[1] != -1 && c.hapConsAlle[2] != -1 &&
                      c.hapConsAlle[1] != c.hapConsAlle[2] && !c.isHomopolymerIndel);
        if (isHet) {
            if (nAgree[vi] < int(opts.minSpanningReads) && nConflict[vi] < int(opts.minSpanningReads))
                ps = int32_t(c.key.pos);
            else if (nConflict[vi] > nAgree[vi]) flip ^= 1;
            if (flip == 1) { changed = true; swap(c.hapConsAlle[1], c.hapConsAlle[2]); }
        }
        c.phaseSet = ps;
    }
    return changed;
}

// ============================================================================
// Step 6c: K-means main loop
// ============================================================================

static void kmRunKmeans(KmScratchpad& scratch, const KmPhasingOptions& opts, uint32_t flags)
{
    const uint32_t numOv = uint32_t(scratch.overlaps.size());
    const uint32_t numCand = uint32_t(scratch.candidates.size());
    if (numCand == 0 || numOv == 0) return;

    scratch.validVarIdx.clear();
    for (uint32_t ci = 0; ci < numCand; ci++)
        if ((scratch.candidates[ci].categoryFlag & flags) != 0)
            scratch.validVarIdx.push_back(ci);
    if (scratch.validVarIdx.empty()) return;

    // Init allele profiles.
    for (uint32_t ci : scratch.validVarIdx) {
        auto& c = scratch.candidates[ci];
        int na = kmAlleSlots(c);
        c.hapAlleProfile[0].assign(na, 0);
        c.hapAlleProfile[1].assign(na, 0);
        c.hapAlleProfile[2].assign(na, 0);
        if (c.category == KmVariantCategory::CleanHom) { c.hapConsAlle[1] = 1; c.hapConsAlle[2] = 1; }
        else { c.hapConsAlle[1] = -1; c.hapConsAlle[2] = -1; }
        c.hapConsAlle[0] = (c.refCov >= c.altCov) ? 0 : 1;
    }

    // Phase 1: seed from pivot.
    int pivotVi = kmSelectPivot(scratch, scratch.validVarIdx);
    if (pivotVi < 0) return;
    { auto& pc = scratch.candidates[scratch.validVarIdx[pivotVi]];
      pc.hapConsAlle[1] = 0; pc.hapConsAlle[2] = 1; }

    for (auto& ov : scratch.overlaps) { ov.hap = 0; ov.phaseSet = -1; }
    for (uint32_t oi = 0; oi < numOv; oi++) {
        int hap = kmAssignOverlapHap(scratch, oi, flags);
        if (hap > 0) { scratch.overlaps[oi].hap = hap; kmUpdateProfiles(scratch, oi, hap, flags); }
        else if (hap == 0) kmUpdateProfiles(scratch, oi, 0, flags);
    }
    for (uint32_t ci : scratch.validVarIdx) {
        kmUpdateConsAlle(scratch.candidates[ci], 1);
        kmUpdateConsAlle(scratch.candidates[ci], 2);
    }

    // Phase 2: iterative refinement.
    for (uint32_t iter = 0; iter < opts.maxKmeansIter; iter++) {
        kmPhaseSetFlip(scratch, scratch.validVarIdx, opts);
        for (uint32_t ci : scratch.validVarIdx) {
            auto& c = scratch.candidates[ci]; int na = kmAlleSlots(c);
            c.hapAlleProfile[0].assign(na, 0);
            c.hapAlleProfile[1].assign(na, 0);
            c.hapAlleProfile[2].assign(na, 0);
        }
        bool anyChanged = false;
        for (uint32_t oi = 0; oi < numOv; oi++) {
            int oldH = scratch.overlaps[oi].hap;
            int newH = kmAssignOverlapHap(scratch, oi, flags);
            if (newH > 0) { scratch.overlaps[oi].hap = newH; kmUpdateProfiles(scratch, oi, newH, flags); }
            else { scratch.overlaps[oi].hap = 0; if (newH == 0) kmUpdateProfiles(scratch, oi, 0, flags); }
            if (scratch.overlaps[oi].hap != oldH) anyChanged = true;
        }
        bool consChanged = false;
        for (uint32_t ci : scratch.validVarIdx) {
            auto& c = scratch.candidates[ci];
            int o1 = c.hapConsAlle[1], o2 = c.hapConsAlle[2];
            kmUpdateConsAlle(c, 1); kmUpdateConsAlle(c, 2);
            if (c.hapConsAlle[1] != o1 || c.hapConsAlle[2] != o2) consChanged = true;
        }
        if (!anyChanged && !consChanged) break;
    }

    // Phase 3: assign phase sets to overlaps.
    for (uint32_t oi = 0; oi < numOv; oi++) {
        auto& ov = scratch.overlaps[oi];
        if (ov.hap == 0) continue;
        const auto& prof = scratch.overlapProfiles[oi];
        if (prof.startVarIdx < 0) continue;
        for (int ci = prof.startVarIdx; ci <= prof.endVarIdx; ci++) {
            const auto& c = scratch.candidates[ci];
            if ((c.categoryFlag & flags) == 0) continue;
            if (c.hapConsAlle[1] != -1 && c.hapConsAlle[2] != -1 && c.hapConsAlle[1] != c.hapConsAlle[2])
                { ov.phaseSet = c.phaseSet; break; }
        }
    }

    // Phase 4: fill hapAlt/hapRef.
    for (uint32_t ci : scratch.validVarIdx) {
        auto& c = scratch.candidates[ci];
        if (c.hapConsAlle[1] == 1) { c.hapAlt = 1; c.hapRef = 2; }
        else if (c.hapConsAlle[2] == 1) { c.hapAlt = 2; c.hapRef = 1; }
    }
}

// ============================================================================
// Step 7: Write results
// ============================================================================

static void kmWriteResults(
    Assembler& assembler, ReadId backboneReadId, const KmScratchpad& scratch)
{
    for (const auto& ov : scratch.overlaps) {
        auto& ad = assembler.alignmentData[ov.alignmentId];
        uint8_t matchState = (ov.hap == 2) ? 2 : 1;
        ad.setHifiasmEcMatchStateFromReadPerspective(backboneReadId, matchState);
    }
}

// ============================================================================
// Public entry point
// ============================================================================

void Assembler::phaseOverlapsKmeans(uint64_t threadCount)
{
    cout << timestamp << "=== K-means Overlap Phasing ===" << endl;
    const uint64_t readCount = getReads().readCount();
    cout << timestamp << "Read count: " << readCount << endl;
    if (readCount == 0) { cout << timestamp << "No reads." << endl; return; }

    KmPhasingOptions opts;
    atomic<uint64_t> readsProcessed(0), readsWithOverlaps(0), readsWithSites(0);
    atomic<uint64_t> totalCis(0), totalTrans(0), totalNoisyRegions(0);

    struct alignas(64) TT {
        int64_t gather=0, unpack=0, detect=0, count=0;
        int64_t classify=0, kmeans=0, write=0;
    };
    vector<TT> ttv(threadCount);

    vector<thread> threads;
    uint64_t chunk = readCount / threadCount;
    if (chunk == 0) chunk = 1;

    for (uint64_t t = 0; t < threadCount; t++) {
        threads.emplace_back([&, t]() {
            uint64_t start = t * chunk;
            if (start >= readCount) return;
            uint64_t end = min((t+1)*chunk, readCount);
            KmScratchpad scratch;
            TT& tt = ttv[t];
            using clk = chrono::steady_clock;
            auto us = [](clk::time_point a, clk::time_point b) {
                return chrono::duration_cast<chrono::microseconds>(b-a).count(); };

            for (uint64_t rid = start; rid < end; rid++) {
                ReadId readId(rid);
                scratch.clear();

                auto t0 = clk::now();
                kmGatherOverlaps(*this, readId, scratch);
                auto t1 = clk::now(); tt.gather += us(t0,t1);
                if (scratch.overlaps.empty()) { readsProcessed++; continue; }
                readsWithOverlaps++;

                t0 = clk::now();
                uint32_t bbLen = uint32_t(getReads().getRead(readId).baseCount);
                scratch.backboneBases.resize(bbLen);
                auto seq = getReads().getRead(readId);
                for (uint32_t i = 0; i < bbLen; i++) scratch.backboneBases[i] = seq[i].value;
                t1 = clk::now(); tt.unpack += us(t0,t1);

                t0 = clk::now();
                kmParseCigars(*this, readId, bbLen, scratch);
                t1 = clk::now(); tt.detect += us(t0,t1);

                t0 = clk::now();
                kmCollectCandidates(scratch, opts.minSvLen);
                kmCountAlleles(scratch, opts.minSvLen);
                t1 = clk::now(); tt.count += us(t0,t1);

                t0 = clk::now();
                kmClassifyCandidates(bbLen, scratch, opts);
                t1 = clk::now(); tt.classify += us(t0,t1);
                totalNoisyRegions += scratch.noisyRegions.size();

                uint32_t cleanHet = 0;
                for (const auto& c : scratch.candidates)
                    if (c.category == KmVariantCategory::CleanHetSnp ||
                        c.category == KmVariantCategory::CleanHetIndel) cleanHet++;
                if (cleanHet > 0) readsWithSites++;

                // Build profiles after classification (skips NON_VAR).
                t0 = clk::now();
                kmBuildOverlapProfiles(scratch, opts.minSvLen);
                t1 = clk::now(); tt.count += us(t0,t1);

                t0 = clk::now();
                if (cleanHet > 0) kmRunKmeans(scratch, opts, KM_GERMLINE_CLEAN);
                t1 = clk::now(); tt.kmeans += us(t0,t1);

                t0 = clk::now();
                kmWriteResults(*this, readId, scratch);
                t1 = clk::now(); tt.write += us(t0,t1);

                for (const auto& ov : scratch.overlaps) {
                    if (ov.hap == 1 || ov.hap == 0) totalCis++;
                    else if (ov.hap == 2) totalTrans++;
                }
                readsProcessed++;
                if (readsProcessed.load() % 10000 == 0)
                    cout << timestamp << "K-means phased " << readsProcessed.load()
                         << " / " << readCount << endl;
            }
        });
    }
    for (auto& t : threads) t.join();

    TT total;
    for (const auto& tt : ttv) {
        total.gather += tt.gather; total.unpack += tt.unpack;
        total.detect += tt.detect; total.count += tt.count;
        total.classify += tt.classify;
        total.kmeans += tt.kmeans; total.write += tt.write;
    }
    auto ms = [](int64_t u) { return u/1000; };
    cout << timestamp << "Timing (ms, sum " << threadCount << " threads):" << endl;
    cout << timestamp << "  gather:   " << ms(total.gather) << endl;
    cout << timestamp << "  unpack:   " << ms(total.unpack) << endl;
    cout << timestamp << "  detect:   " << ms(total.detect) << endl;
    cout << timestamp << "  count+prof:" << ms(total.count) << endl;
    cout << timestamp << "  classify: " << ms(total.classify) << endl;
    cout << timestamp << "  kmeans:   " << ms(total.kmeans) << endl;
    cout << timestamp << "  write:    " << ms(total.write) << endl;
    cout << timestamp << "Complete. Reads=" << readsProcessed.load()
         << " withOvlp=" << readsWithOverlaps.load()
         << " withSites=" << readsWithSites.load()
         << " cis=" << totalCis.load()
         << " trans=" << totalTrans.load()
         << " noisyRegs=" << totalNoisyRegions.load() << endl;
}
