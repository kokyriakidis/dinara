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
// Step 2: Walk CIGARs, detect candidate sites
// ============================================================================

static void kmDetectCandidateSites(
    const Assembler& assembler, ReadId backboneReadId,
    uint32_t backboneLen, KmScratchpad& scratch)
{
    const auto& cigarStore = assembler.getOverlapCigarStore();
    const uint32_t numOv = uint32_t(scratch.overlaps.size());

    scratch.mismatchVotes.assign(backboneLen, 0);
    scratch.indelVotes.assign(backboneLen, 0);
    scratch.mismatchRecords.clear();
    scratch.indelRecords.clear();

    for (uint32_t oi = 0; oi < numOv; oi++) {
        const auto& ov = scratch.overlaps[oi];
        const auto& ad = assembler.alignmentData[ov.alignmentId];
        const bool needsRcConvert = (ov.queryIsRead0 == 0) && (ov.isRev != 0);

        cigarStore.forEachOpWithPositions(
            ov.cigarOffset, ov.cigarTokenCount,
            ad.qs, kmCigarRead1Start(assembler, ad),
            [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
                if (op == 1) { // Mismatch
                    if (!needsRcConvert) {
                        const uint32_t qStart = ov.queryIsRead0 ? uint32_t(xk) : uint32_t(yk);
                        for (uint32_t b = 0; b < len; b++) {
                            const uint32_t qpos = qStart + b;
                            if (qpos >= backboneLen) break;
                            if (scratch.mismatchVotes[qpos] < 255) scratch.mismatchVotes[qpos]++;
                            const uint32_t tpos = ov.queryIsRead0 ? uint32_t(yk)+b : uint32_t(xk)+b;
                            scratch.mismatchRecords.push_back({qpos, tpos, oi});
                        }
                    } else {
                        for (uint32_t b = 0; b < len; b++) {
                            const uint32_t qpos_rc = uint32_t(yk) + b;
                            if (qpos_rc >= backboneLen) continue;
                            const uint32_t qpos = backboneLen - 1 - qpos_rc;
                            if (scratch.mismatchVotes[qpos] < 255) scratch.mismatchVotes[qpos]++;
                            scratch.mismatchRecords.push_back({qpos, uint32_t(xk)+b, oi});
                        }
                    }
                } else if (op == 2 || op == 3) { // Ins/Del
                    bool bbConsumed = (op==3 && ov.queryIsRead0) || (op==2 && !ov.queryIsRead0);
                    if (bbConsumed) {
                        const uint32_t rawStart = (op==3) ? uint32_t(xk) : uint32_t(yk);
                        uint32_t rawEnd = min(rawStart + len, backboneLen);
                        if (!needsRcConvert) {
                            for (uint32_t p = rawStart; p < rawEnd; p++)
                                if (scratch.indelVotes[p] < 255) scratch.indelVotes[p]++;
                            scratch.indelRecords.push_back({rawStart, rawEnd, oi, op, uint16_t(len)});
                        } else {
                            uint32_t fS = backboneLen - rawEnd, fE = backboneLen - rawStart;
                            for (uint32_t p = fS; p < fE && p < backboneLen; p++)
                                if (scratch.indelVotes[p] < 255) scratch.indelVotes[p]++;
                            scratch.indelRecords.push_back({fS, fE, oi, op, uint16_t(len)});
                        }
                    } else {
                        uint32_t anchor = ov.queryIsRead0 ? uint32_t(xk) : uint32_t(yk);
                        if (anchor < backboneLen) {
                            if (scratch.indelVotes[anchor] < 255) scratch.indelVotes[anchor]++;
                            scratch.indelRecords.push_back({anchor, anchor+1, oi, op, uint16_t(len)});
                        }
                    }
                }
            });
    }

    // Positions with >= 2 votes become candidates.
    for (uint32_t pos = 0; pos < backboneLen; pos++) {
        if (scratch.mismatchVotes[pos] >= 2 || scratch.indelVotes[pos] >= 2) {
            KmCandidate cand{};
            cand.pos = pos;
            cand.refBase = scratch.backboneBases[pos];
            cand.type = (scratch.mismatchVotes[pos] >= scratch.indelVotes[pos]) ? 0 : 2;
            cand.refLen = 1;
            cand.category = KmVariantCategory::LowCoverage;
            scratch.candidates.push_back(cand);
        }
    }
}

// ============================================================================
// Step 3: Count alleles per candidate site
// ============================================================================

static void kmCountAlleles(
    const Assembler& assembler, ReadId backboneReadId,
    uint32_t backboneLen, KmScratchpad& scratch)
{
    if (scratch.candidates.empty()) return;
    const uint32_t numOv = uint32_t(scratch.overlaps.size());
    const uint32_t numCand = uint32_t(scratch.candidates.size());

    // Position → candidate index.
    vector<uint32_t> posToCand(backboneLen, UINT32_MAX);
    for (uint32_t ci = 0; ci < numCand; ci++)
        posToCand[scratch.candidates[ci].pos] = ci;

    // Per-overlap × candidate: has alt?
    vector<uint8_t> ovHasAlt(size_t(numOv) * numCand, 0);
    for (const auto& mr : scratch.mismatchRecords) {
        if (mr.backbonePos >= backboneLen) continue;
        uint32_t ci = posToCand[mr.backbonePos];
        if (ci != UINT32_MAX) ovHasAlt[size_t(mr.overlapIdx) * numCand + ci] = 1;
    }
    for (const auto& ir : scratch.indelRecords) {
        for (uint32_t p = ir.backboneStart; p < ir.backboneEnd && p < backboneLen; p++) {
            uint32_t ci = posToCand[p];
            if (ci != UINT32_MAX) ovHasAlt[size_t(ir.overlapIdx) * numCand + ci] = 1;
        }
    }

    for (uint32_t ci = 0; ci < numCand; ci++) {
        auto& c = scratch.candidates[ci];
        uint32_t ref = 0, alt = 0, fR = 0, rR = 0, fA = 0, rA = 0, tot = 0;
        for (uint32_t oi = 0; oi < numOv; oi++) {
            const auto& ov = scratch.overlaps[oi];
            if (c.pos < ov.qs || c.pos >= ov.qe) continue;
            tot++;
            if (ovHasAlt[size_t(oi) * numCand + ci]) {
                alt++; if (ov.isRev == 0) fA++; else rA++;
            } else {
                ref++; if (ov.isRev == 0) fR++; else rR++;
            }
        }
        ref++; tot++; fR++; // backbone itself
        c.totalCov = tot; c.refCov = ref; c.altCov = alt;
        c.fwdRef = fR; c.revRef = rR; c.fwdAlt = fA; c.revAlt = rA;
        c.alleleFraction = tot > 0 ? double(alt) / double(tot) : 0.0;
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

    // Pass 1: per-site classification.
    for (uint32_t ci = 0; ci < numCand; ci++) {
        auto& c = scratch.candidates[ci];
        if (c.totalCov < opts.minDepth || c.altCov < opts.minAltDepth) {
            c.category = KmVariantCategory::LowCoverage; c.categoryFlag = 0; continue;
        }
        // Strand bias: all alt on one strand with >= 4 alt reads.
        { int fa = int(c.fwdAlt), ra = int(c.revAlt);
          if ((fa == 0 || ra == 0) && (fa + ra) >= 4) {
              c.category = KmVariantCategory::StrandBias; c.categoryFlag = 0; continue; } }
        if (c.alleleFraction < opts.minAf) {
            c.category = KmVariantCategory::LowAlleleFraction; c.categoryFlag = 0; continue;
        }
        if (c.alleleFraction > opts.maxAf) {
            c.category = KmVariantCategory::CleanHom; c.categoryFlag = KM_CLEAN_HOM; continue;
        }
        if (c.type != 0) { // indel
            if (kmIsPeriodicRepeat(scratch.backboneBases.data(), backboneLen, c.pos)) {
                c.category = KmVariantCategory::RepeatHetIndel; c.categoryFlag = KM_REP_HET_VAR; continue;
            }
            c.category = KmVariantCategory::CleanHetIndel; c.categoryFlag = KM_CLEAN_HET_INDEL; continue;
        }
        // SNP
        if (kmIsPeriodicRepeat(scratch.backboneBases.data(), backboneLen, c.pos)) {
            c.category = KmVariantCategory::RepeatHetIndel; c.categoryFlag = KM_REP_HET_VAR; continue;
        }
        c.category = KmVariantCategory::CleanHetSnp; c.categoryFlag = KM_CLEAN_HET_SNP;
    }

    // Pass 2: noisy region seeds from repeat-het and dense clusters.
    scratch.noisyRegions.clear();
    vector<pair<uint32_t,uint32_t>> seeds;
    for (uint32_t ci = 0; ci < numCand; ci++) {
        if (scratch.candidates[ci].category == KmVariantCategory::RepeatHetIndel)
            seeds.push_back({scratch.candidates[ci].pos, scratch.candidates[ci].pos + scratch.candidates[ci].refLen});
    }
    for (uint32_t ci = 1; ci < numCand; ci++) {
        const auto& prev = scratch.candidates[ci-1];
        const auto& cur = scratch.candidates[ci];
        if (prev.category == KmVariantCategory::LowCoverage || prev.category == KmVariantCategory::NonVariant) continue;
        if (cur.category == KmVariantCategory::LowCoverage || cur.category == KmVariantCategory::NonVariant) continue;
        if (cur.pos - prev.pos <= 25) {
            if (prev.category != KmVariantCategory::CleanHetSnp || cur.category != KmVariantCategory::CleanHetSnp)
                seeds.push_back({prev.pos, cur.pos + cur.refLen});
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

    // Demote candidates inside noisy regions.
    size_t ri = 0;
    for (uint32_t ci = 0; ci < numCand && ri < scratch.noisyRegions.size(); ci++) {
        auto& c = scratch.candidates[ci];
        while (ri < scratch.noisyRegions.size() && scratch.noisyRegions[ri].end <= c.pos) ri++;
        if (ri < scratch.noisyRegions.size() &&
            c.pos >= scratch.noisyRegions[ri].start && c.pos < scratch.noisyRegions[ri].end) {
            c.category = KmVariantCategory::NonVariant; c.categoryFlag = KM_NON_VAR;
        }
    }
}

// ============================================================================
// Step 5: Build per-overlap allele profiles
// ============================================================================

static void kmBuildOverlapProfiles(
    uint32_t backboneLen, KmScratchpad& scratch)
{
    const uint32_t numOv = uint32_t(scratch.overlaps.size());
    const uint32_t numCand = uint32_t(scratch.candidates.size());
    if (numCand == 0 || numOv == 0) return;

    vector<uint32_t> posToCand(backboneLen, UINT32_MAX);
    for (uint32_t ci = 0; ci < numCand; ci++)
        posToCand[scratch.candidates[ci].pos] = ci;

    vector<uint8_t> ovHasAlt(size_t(numOv) * numCand, 0);
    for (const auto& mr : scratch.mismatchRecords) {
        if (mr.backbonePos >= backboneLen) continue;
        uint32_t ci = posToCand[mr.backbonePos];
        if (ci != UINT32_MAX) ovHasAlt[size_t(mr.overlapIdx) * numCand + ci] = 1;
    }
    for (const auto& ir : scratch.indelRecords) {
        for (uint32_t p = ir.backboneStart; p < ir.backboneEnd && p < backboneLen; p++) {
            uint32_t ci = posToCand[p];
            if (ci != UINT32_MAX) ovHasAlt[size_t(ir.overlapIdx) * numCand + ci] = 1;
        }
    }

    scratch.overlapProfiles.resize(numOv);
    for (uint32_t oi = 0; oi < numOv; oi++) {
        auto& prof = scratch.overlapProfiles[oi];
        prof.overlapIdx = oi;
        prof.startVarIdx = -1; prof.endVarIdx = -1;
        prof.alleles.clear();
        const auto& ov = scratch.overlaps[oi];
        for (uint32_t ci = 0; ci < numCand; ci++) {
            const auto& cand = scratch.candidates[ci];
            if (cand.categoryFlag == KM_NON_VAR) continue;
            if (cand.pos < ov.qs || cand.pos >= ov.qe) continue;
            if (prof.startVarIdx < 0) prof.startVarIdx = int(ci);
            while (prof.startVarIdx + int(prof.alleles.size()) < int(ci))
                prof.alleles.push_back(-1);
            prof.endVarIdx = int(ci);
            prof.alleles.push_back(ovHasAlt[size_t(oi) * numCand + ci] ? 1 : 0);
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
    int bSnp = -1, bIndel = -1; uint32_t bSnpD = 0, bIndelD = 0;
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
    int32_t ps = int32_t(scratch.candidates[validIdx[0]].pos);
    scratch.candidates[validIdx[0]].phaseSet = ps;
    for (int vi = 1; vi < n; vi++) {
        auto& c = scratch.candidates[validIdx[vi]];
        bool isHet = (c.hapConsAlle[1] != -1 && c.hapConsAlle[2] != -1 &&
                      c.hapConsAlle[1] != c.hapConsAlle[2] && !c.isHomopolymerIndel);
        if (isHet) {
            if (nAgree[vi] < int(opts.minSpanningReads) && nConflict[vi] < int(opts.minSpanningReads))
                ps = int32_t(c.pos);
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
        int64_t classify=0, profile=0, kmeans=0, write=0;
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
                kmDetectCandidateSites(*this, readId, bbLen, scratch);
                t1 = clk::now(); tt.detect += us(t0,t1);

                t0 = clk::now();
                kmCountAlleles(*this, readId, bbLen, scratch);
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

                t0 = clk::now();
                kmBuildOverlapProfiles(bbLen, scratch);
                t1 = clk::now(); tt.profile += us(t0,t1);

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
        total.classify += tt.classify; total.profile += tt.profile;
        total.kmeans += tt.kmeans; total.write += tt.write;
    }
    auto ms = [](int64_t u) { return u/1000; };
    cout << timestamp << "Timing (ms, sum " << threadCount << " threads):" << endl;
    cout << timestamp << "  gather:   " << ms(total.gather) << endl;
    cout << timestamp << "  unpack:   " << ms(total.unpack) << endl;
    cout << timestamp << "  detect:   " << ms(total.detect) << endl;
    cout << timestamp << "  count:    " << ms(total.count) << endl;
    cout << timestamp << "  classify: " << ms(total.classify) << endl;
    cout << timestamp << "  profile:  " << ms(total.profile) << endl;
    cout << timestamp << "  kmeans:   " << ms(total.kmeans) << endl;
    cout << timestamp << "  write:    " << ms(total.write) << endl;
    cout << timestamp << "Complete. Reads=" << readsProcessed.load()
         << " withOvlp=" << readsWithOverlaps.load()
         << " withSites=" << readsWithSites.load()
         << " cis=" << totalCis.load()
         << " trans=" << totalTrans.load()
         << " noisyRegs=" << totalNoisyRegions.load() << endl;
}
