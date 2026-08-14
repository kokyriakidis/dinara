/// @file AssemblerPhasingKmeans.cpp
/// @brief K-means overlap phasing adapted from pgphase/longcallD.

#include "Assembler.hpp"
#include "PhasingKmeansAlign.hpp"
#include "PhasingKmeansTypes.hpp"
#include "Alignment.hpp"
#include "Reads.hpp"
#include "invalid.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <atomic>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iostream>
#include <mutex>
#include <numeric>
#include <thread>
#include <vector>

using namespace std;
using namespace dinara;

// ============================================================================
// SDUST low-complexity detection (C++ port of heng li's sdust.c from pgphase)
//
// Detects low-complexity regions using symmetric DUST algorithm.
// Input: 2-bit encoded sequence (0=A, 1=C, 2=G, 3=T, >=4 = N/break).
// Output: vector of (start, end) intervals (0-based, half-open).
// Parameters: T=threshold (pgphase default 5), W=window (pgphase default 20).
// ============================================================================

static constexpr int SD_WLEN = 3;
static constexpr int SD_WTOT = (1 << (SD_WLEN << 1));  // 64
static constexpr int SD_WMSK = SD_WTOT - 1;

struct SdustPerfIntv { int start, finish, r, l; };

/// Port of sdust.c shift_window.
static inline void kmSdustShiftWindow(
    int t, deque<int>& w, int T, int W,
    int& L, int& rw, int& rv, int* cw, int* cv)
{
    int s;
    if (int(w.size()) >= W - SD_WLEN + 1) {
        s = w.front(); w.pop_front();
        rw -= --cw[s];
        if (L > int(w.size()))
            --L, rv -= --cv[s];
    }
    w.push_back(t);
    ++L;
    rw += cw[t]++;
    rv += cv[t]++;
    if (cv[t] * 10 > T << 1) {
        do {
            s = w[w.size() - L];
            rv -= --cv[s];
            --L;
        } while (s != t);
    }
}

/// Port of sdust.c save_masked_regions.
static inline void kmSdustSaveMasked(
    vector<pair<uint32_t,uint32_t>>& res, vector<SdustPerfIntv>& P, int start)
{
    if (P.empty() || P.back().start >= start) return;
    auto& p = P.back();
    int saved = 0;
    if (!res.empty()) {
        uint32_t f = res.back().second;
        if (uint32_t(p.start) <= f) {
            saved = 1;
            res.back().second = (f > uint32_t(p.finish)) ? f : uint32_t(p.finish);
        }
    }
    if (!saved)
        res.push_back({uint32_t(p.start), uint32_t(p.finish)});
    int i;
    for (i = int(P.size()) - 1; i >= 0 && P[i].start < start; --i);
    P.resize(i + 1);
}

/// Port of sdust.c find_perfect.
static void kmSdustFindPerfect(
    vector<SdustPerfIntv>& P, const deque<int>& w,
    int T, int start, int L, int rv, const int* cv)
{
    int c[SD_WTOT], r = rv;
    memcpy(c, cv, SD_WTOT * sizeof(int));
    int max_r = 0, max_l = 0;
    for (int i = int(w.size()) - L - 1; i >= 0; --i) {
        int t = w[i];
        r += c[t]++;
        int new_r = r, new_l = int(w.size()) - i - 1;
        if (new_r * 10 > T * new_l) {
            int j;
            for (j = 0; j < int(P.size()) && P[j].start >= i + start; ++j) {
                auto& p = P[j];
                if (max_r == 0 || p.r * max_l > max_r * p.l)
                    max_r = p.r, max_l = p.l;
            }
            if (max_r == 0 || new_r * max_l >= max_r * new_l) {
                max_r = new_r; max_l = new_l;
                SdustPerfIntv pi;
                pi.start = i + start;
                pi.finish = int(w.size()) + (SD_WLEN - 1) + start;
                pi.r = new_r; pi.l = new_l;
                P.insert(P.begin() + j, pi);
            }
        }
    }
}

/// Port of sdust.c sdust_core. Input is already 2-bit encoded (0-3).
static void kmSdust(const uint8_t* seq, int lSeq, int T, int W,
                    vector<pair<uint32_t,uint32_t>>& out)
{
    out.clear();
    if (lSeq < SD_WLEN) return;

    deque<int> w;
    vector<SdustPerfIntv> P;
    int rv = 0, rw = 0, L = 0;
    int cv[SD_WTOT], cw[SD_WTOT];
    memset(cv, 0, sizeof(cv));
    memset(cw, 0, sizeof(cw));
    unsigned t = 0;
    int l = 0;

    for (int i = 0; i <= lSeq; ++i) {
        // Our backbone bases are already 0-3 encoded; treat >=4 as N/break.
        int b = (i < lSeq) ? int(seq[i]) : 4;
        if (b < 4) {
            ++l;
            t = (t << 2 | b) & SD_WMSK;
            if (l >= SD_WLEN) {
                int start = (l - W > 0 ? l - W : 0) + (i + 1 - l);
                kmSdustSaveMasked(out, P, start);
                kmSdustShiftWindow(t, w, T, W, L, rw, rv, cw, cv);
                if (rw * 10 > L * T)
                    kmSdustFindPerfect(P, w, T, start, L, rv, cv);
            }
        } else {
            // pgphase sdust.c: only reset l and t on N/break.
            // The window state (w, L, rv, rw, cv, cw) is NOT reset —
            // old entries age out naturally via shift_window.
            int start = (l - W + 1 > 0 ? l - W + 1 : 0) + (i + 1 - l);
            while (!P.empty()) kmSdustSaveMasked(out, P, start++);
            l = 0; t = 0;
        }
    }
}

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

/// pgphase var_is_homopolymer_pg: checks if flanking reference has a periodic
/// repeat pattern (unit length 1-6, 3 copies) in either direction from the
/// variant boundary.
// kmIsHomopolymer and kmIsRepeatRegion are now inline in PhasingKmeansTypes.hpp.

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
// Step 2: Parse CIGARs into per-overlap digars + noisy region detection
//
// Port of pgphase NoisyRegionBuilder / XidQueue sliding window.
// During CIGAR walk, variant events are pushed into a sliding window.
// When total event weight in the window exceeds max_s (noisyRegMaxXgaps),
// the window span is flagged as a noisy interval for this overlap.
// ============================================================================

/// Sliding window noisy region detector (pgphase XidQueue + xid_push_win).
struct KmNoisyBuilder {
    // Circular queue of events in the window.
    vector<uint32_t> pos;    // backbone position
    vector<uint32_t> lens;   // ref-space length (0 for insertions)
    vector<int>      counts; // event weight (1 for SNP, len for indel)
    int front = 0;
    int rear = -1;
    int totalCount = 0;
    int maxS;   // threshold: noisy when totalCount > maxS
    int win;    // sliding window size in bp

    // Active noisy interval tracking.
    int32_t curStart = -1, curEnd = -1;
    int curQStart = -1, curQEnd = -1;

    // Output.
    vector<KmNoisyRegion>& noisyOut;

    KmNoisyBuilder(int maxSites, int maxSIn, int winIn, vector<KmNoisyRegion>& out)
        : pos(max(1, maxSites)), lens(max(1, maxSites)), counts(max(1, maxSites)),
          maxS(maxSIn), win(winIn), noisyOut(out) {}

    void ensureCapacity() {
        if (size_t(rear + 1) < pos.size()) return;
        size_t nc = pos.size() * 2;
        pos.resize(nc); lens.resize(nc); counts.resize(nc);
    }

    void observe(uint32_t p, uint32_t len, int count) {
        ensureCapacity();
        ++rear;
        pos[rear] = p; lens[rear] = len; counts[rear] = count;
        totalCount += count;

        // Evict events that fell behind the window.
        while (front <= rear &&
               int64_t(pos[front]) + int64_t(lens[front]) - 1 <= int64_t(p) - win) {
            totalCount -= counts[front];
            ++front;
        }

        if (count <= 0) return;
        if (totalCount <= maxS) return;

        int32_t noisyStart = int32_t(pos[front]);
        int32_t noisyEnd = int32_t(pos[rear] + lens[rear]);

        if (curStart == -1) {
            curStart = noisyStart; curEnd = noisyEnd;
            curQStart = front; curQEnd = rear;
            return;
        }
        if (noisyStart <= curEnd) {
            curEnd = noisyEnd; curQEnd = rear;
            return;
        }
        // Flush previous interval.
        int varSize = 0;
        for (int i = curQStart; i <= curQEnd; i++) varSize += counts[i];
        int span = int(curEnd - curStart + 1);
        if (varSize < span) varSize = span;
        noisyOut.push_back({uint32_t(curStart), uint32_t(curEnd), varSize, false});

        curStart = noisyStart; curEnd = noisyEnd;
        curQStart = front; curQEnd = rear;
    }

    void flush() {
        if (curStart == -1) return;
        int varSize = 0;
        for (int i = curQStart; i <= curQEnd; i++) varSize += counts[i];
        int span = int(curEnd - curStart + 1);
        if (varSize < span) varSize = span;
        noisyOut.push_back({uint32_t(curStart), uint32_t(curEnd), varSize, false});
        curStart = -1;
    }
};

static void kmParseCigars(
    const Assembler& assembler, ReadId backboneReadId,
    uint32_t backboneLen, KmScratchpad& scratch, const KmPhasingOptions& opts)
{
    const auto& cigarStore = assembler.getOverlapCigarStore();
    const uint32_t numOv = uint32_t(scratch.overlaps.size());

    scratch.digars.clear();
    scratch.digarBegin.resize(numOv);
    scratch.digarEnd.resize(numOv);
    scratch.overlapNoisyRegions.clear();
    scratch.overlapNoisyBegin.resize(numOv);
    scratch.overlapNoisyEnd.resize(numOv);

    // pgphase: window size depends on technology.
    const int noisyWin = opts.isOnt ? 25 : 100;  // kDefaultNoisyRegSlideWinOnt / Hifi
    const int noisyMaxS = int(opts.noisyRegMaxXgaps);

    for (uint32_t oi = 0; oi < numOv; oi++) {
        scratch.digarBegin[oi] = uint32_t(scratch.digars.size());
        scratch.overlapNoisyBegin[oi] = uint32_t(scratch.overlapNoisyRegions.size());
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
                        // Get the target read's base in the backbone's forward frame.
                        // Read0 is always strand 0; read1's strand depends on isSameStrand.
                        //
                        // When target=read1 (qIsR0=1):
                        //   tPos is read1's oriented position (forward if same-strand,
                        //   RC if opposite). kmGetBase with isRc handles both correctly
                        //   because it converts oriented pos → forward pos + complement.
                        //
                        // When target=read0 (qIsR0=0):
                        //   tPos is read0's forward position (read0 is always strand 0).
                        //   If strands differ (isRev), read0 is on the opposite strand
                        //   from the backbone. We need to complement the base to put it
                        //   in the backbone's frame, but NOT mirror the position (tPos
                        //   is already a valid forward index into read0's sequence).
                        uint8_t altBase;
                        if (qIsR0) {
                            // Target is read1: use kmGetBase which handles oriented positions.
                            altBase = kmGetBase(assembler, ReadId(ov.targetReadId),
                                                tPos, ov.isRev != 0);
                        } else {
                            // Target is read0 (always strand 0): read forward base directly.
                            const auto seq = assembler.getReads().getRead(ReadId(ov.targetReadId));
                            altBase = seq[tPos].value;
                            // If strands differ, complement to match backbone's frame.
                            if (ov.isRev) altBase = uint8_t((~altBase) & 3);
                        }
                        scratch.digars.push_back({bbPos, KmVarType::Snp, altBase, 1, {}});
                    }
                } else if (op == CigarOpIns || op == CigarOpDel) { // Ins/Del
                    // Backbone read (xk if qIsR0, else yk) consumed by this indel
                    // iff the op consumes that read's sequence → deletion on the
                    // backbone; otherwise an insertion at a backbone position.
                    const bool bbConsumed = qIsR0 ? opConsumesQuery(op)
                                                  : opConsumesTarget(op);
                    if (bbConsumed) {
                        uint32_t raw = qIsR0 ? uint32_t(xk) : uint32_t(yk);
                        uint32_t rawE = (raw + len < backboneLen) ? raw + len : backboneLen;
                        uint32_t s, e;
                        if (needsRc) { s = backboneLen - rawE; e = backboneLen - raw; }
                        else { s = raw; e = rawE; }
                        if (s < backboneLen) {
                            scratch.digars.push_back({s, KmVarType::Deletion, 0, uint16_t(e - s), {}});
                        }
                    } else {
                        uint32_t anchor = qIsR0 ? uint32_t(xk) : uint32_t(yk);
                        // For RC, convert boundary (not position): b → backboneLen - b.
                        // pgphase never needs this because BAM CIGARs are reference-forward.
                        if (needsRc) anchor = backboneLen - anchor;
                        if (anchor < backboneLen) {
                            // Extract inserted bases from the non-backbone read
                            // in the backbone's forward frame.
                            std::string insSeq;
                            insSeq.reserve(len);
                            // The non-backbone read's positions for this insertion:
                            //   qIsR0: non-bb is read1, positions yk..yk+len-1
                            //   !qIsR0: non-bb is read0, positions xk..xk+len-1
                            const uint32_t insStart = qIsR0 ? uint32_t(yk) : uint32_t(xk);
                            for (uint32_t b = 0; b < len; b++) {
                                uint8_t base;
                                if (qIsR0) {
                                    base = kmGetBase(assembler, ReadId(ov.targetReadId),
                                                     insStart + b, ov.isRev != 0);
                                } else {
                                    const auto seq = assembler.getReads().getRead(ReadId(ov.targetReadId));
                                    base = seq[insStart + b].value;
                                    if (ov.isRev) base = uint8_t((~base) & 3);
                                }
                                insSeq.push_back("ACGT"[base & 3]);
                            }
                            // For RC overlaps (!qIsR0 && isRev), the bases were
                            // complemented above but are still in read0's forward
                            // order. Reverse to put them in the backbone's frame.
                            if (needsRc) {
                                std::reverse(insSeq.begin(), insSeq.end());
                            }
                            scratch.digars.push_back({anchor, KmVarType::Insertion, 0, uint16_t(len), std::move(insSeq)});
                        }
                    }
                }
            });

        // Sort this overlap's digars by backbone position.
        uint32_t dEnd = uint32_t(scratch.digars.size());
        sort(scratch.digars.begin() + scratch.digarBegin[oi],
             scratch.digars.begin() + dEnd);

        // Run noisy region detection on sorted digars so the sliding window
        // sees events in backbone-position order (required for RC overlaps).
        int maxSites = max(1, int(dEnd - scratch.digarBegin[oi]) * 2 + 8);
        KmNoisyBuilder noisy(maxSites, noisyMaxS, noisyWin, scratch.overlapNoisyRegions);
        for (uint32_t di = scratch.digarBegin[oi]; di < dEnd; di++) {
            const auto& d = scratch.digars[di];
            if (d.type == KmVarType::Snp)
                noisy.observe(d.pos, 1, 1);
            else if (d.type == KmVarType::Deletion)
                noisy.observe(d.pos, d.len, int(d.len));
            else
                noisy.observe(d.pos, 0, int(d.len));
        }
        noisy.flush();
        scratch.digarEnd[oi] = dEnd;
        scratch.overlapNoisyEnd[oi] = uint32_t(scratch.overlapNoisyRegions.size());
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
            // Short insertion: exact length + sequence match required.
            if (a.altLen != b.altLen) return a.altLen < b.altLen ? -1 : 1;
            if (a.altSeq != b.altSeq) return a.altSeq < b.altSeq ? -1 : 1;
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
        k.altSeq = d.altSeq;
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
// Interval helpers (used by noisy region processing and classification)
// ============================================================================

/// Check if position range [qStart, qEnd) overlaps any interval in sorted list.
static bool kmOverlapsAny(const vector<pair<uint32_t,uint32_t>>& intervals,
                           uint32_t qStart, uint32_t qEnd)
{
    auto it = lower_bound(intervals.begin(), intervals.end(), qStart,
        [](const pair<uint32_t,uint32_t>& iv, uint32_t val) { return iv.second <= val; });
    return it != intervals.end() && it->first < qEnd;
}

/// Count how many intervals in sorted list overlap [qStart, qEnd).
static int kmCountOverlaps(const vector<pair<uint32_t,uint32_t>>& intervals,
                            uint32_t qStart, uint32_t qEnd)
{
    auto it = lower_bound(intervals.begin(), intervals.end(), qStart,
        [](const pair<uint32_t,uint32_t>& iv, uint32_t val) { return iv.second <= val; });
    int count = 0;
    while (it != intervals.end() && it->first < qEnd) { count++; ++it; }
    return count;
}

/// Check if [qStart, qEnd) is fully contained in some interval.
static bool kmIsContained(const vector<pair<uint32_t,uint32_t>>& intervals,
                           uint32_t qStart, uint32_t qEnd)
{
    auto it = lower_bound(intervals.begin(), intervals.end(), qStart,
        [](const pair<uint32_t,uint32_t>& iv, uint32_t val) { return iv.second <= val; });
    while (it != intervals.end() && it->first <= qStart) {
        if (it->second >= qEnd) return true;
        ++it;
    }
    return false;
}

/// Extend [start, end) to include overlapping low-complexity intervals.
/// Port of pgphase cr_add_var_to_noisy_cr low_comp extension.
static void kmExtendWithLowComp(const vector<pair<uint32_t,uint32_t>>& lowComp,
                                 uint32_t& start, uint32_t& end)
{
    if (lowComp.empty()) return;
    // Find low-complexity intervals overlapping [start, end).
    auto it = lower_bound(lowComp.begin(), lowComp.end(), start,
        [](const pair<uint32_t,uint32_t>& iv, uint32_t val) { return iv.second <= val; });
    while (it != lowComp.end() && it->first < end) {
        if (it->first < start) start = it->first;
        if (it->second > end) end = it->second;
        ++it;
    }
}

/// Merge overlapping/adjacent intervals (sorted input, merge distance = mergeDis).
static void kmMergeIntervals(vector<pair<uint32_t,uint32_t>>& intervals, uint32_t mergeDis)
{
    if (intervals.size() <= 1) return;
    sort(intervals.begin(), intervals.end());
    vector<pair<uint32_t,uint32_t>> merged;
    merged.push_back(intervals[0]);
    for (size_t i = 1; i < intervals.size(); i++) {
        auto& last = merged.back();
        if (intervals[i].first <= last.second + mergeDis)
            last.second = max(last.second, intervals[i].second);
        else
            merged.push_back(intervals[i]);
    }
    intervals = std::move(merged);
}

// ============================================================================
// Step 2b: Pre-process noisy regions (pgphase pre_process_noisy_regs)
//
// Merges per-overlap noisy intervals into chunk-level noisy regions.
// Filters by read support: requires >= minAltDepth overlaps with noisy
// intervals overlapping each merged region, and noisy ratio >= minAf.
// These chunk-level noisy regions feed into classify pass 2.
// ============================================================================

static void kmPreProcessNoisyRegs(KmScratchpad& scratch, const KmPhasingOptions& opts)
{
    const uint32_t numOv = uint32_t(scratch.overlaps.size());
    scratch.noisyRegions.clear();

    // Collect all per-overlap noisy intervals into a single sorted list.
    vector<pair<uint32_t,uint32_t>> allNoisy;
    for (uint32_t oi = 0; oi < numOv; oi++) {
        for (uint32_t ni = scratch.overlapNoisyBegin[oi]; ni < scratch.overlapNoisyEnd[oi]; ni++) {
            const auto& nr = scratch.overlapNoisyRegions[ni];
            allNoisy.push_back({nr.start, nr.end});
        }
    }
    if (allNoisy.empty()) return;

    // pgphase: cr_extend_noisy_regs_with_low_comp — extend each noisy interval
    // to include overlapping low-complexity (sdust) intervals, then merge.
    for (auto& [s, e] : allNoisy)
        kmExtendWithLowComp(scratch.lowComplexity, s, e);

    // Merge overlapping intervals (pgphase cr_merge).
    kmMergeIntervals(allNoisy, opts.noisyRegMergeDis);

    // For each merged region, count how many overlaps have noisy intervals
    // overlapping it (pgphase: noisy_reg_to_noisy) and how many overlaps
    // span it at all (pgphase: noisy_reg_to_total).
    const int minNoisyReads = int(opts.minAltDepth);
    const float minNoisyRatio = float(opts.minAf);
    // longcallD only checks n_noisy < min_alt_dp and ratio < min_af.

    for (const auto& [mStart, mEnd] : allNoisy) {
        int nTotal = 0;
        int nNoisy = 0;
        for (uint32_t oi = 0; oi < numOv; oi++) {
            const auto& ov = scratch.overlaps[oi];
            // pgphase: any overlap between read span and merged region counts.
            if (ov.qs >= mEnd || ov.qe <= mStart) continue;
            nTotal++;
            // Does this overlap have a noisy interval overlapping the merged region?
            for (uint32_t ni = scratch.overlapNoisyBegin[oi]; ni < scratch.overlapNoisyEnd[oi]; ni++) {
                const auto& nr = scratch.overlapNoisyRegions[ni];
                if (nr.start < mEnd && nr.end > mStart) { nNoisy++; break; }
            }
        }
        // pgphase filter: skip if insufficient support.
        if (nNoisy < minNoisyReads) continue;
        if (nTotal > 0 && float(nNoisy) / float(nTotal) < minNoisyRatio) continue;
        scratch.noisyRegions.push_back({mStart, mEnd, 0, false});
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

static void kmBuildOverlapProfiles(KmScratchpad& scratch, int minSvLen,
    const vector<uint8_t>* skipMask = nullptr)
{
    const uint32_t numCand = uint32_t(scratch.candidates.size());
    const uint32_t numOv = uint32_t(scratch.overlaps.size());
    if (numCand == 0 || numOv == 0) return;

    scratch.overlapProfiles.resize(numOv);

    for (uint32_t oi = 0; oi < numOv; oi++) {
        auto& prof = scratch.overlapProfiles[oi];
        prof.overlapIdx = oi;
        prof.startVarIdx = -1; prof.endVarIdx = -1;
        prof.alleles.clear();

        // Skip overlaps excluded by the mask.
        if (skipMask && !(*skipMask)[oi]) continue;

        const auto& ov = scratch.overlaps[oi];

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
                // No digar at this position.
                // pgphase: during merge walk, always count as ref (allele=0).
                // The noisy region check only applies in the tail loop
                // (after all digars consumed), handled below.
                if (di >= diEnd) {
                    // Tail: all digars consumed. Check per-overlap noisy region.
                    bool inOverlapNoisy = false;
                    for (uint32_t ni = scratch.overlapNoisyBegin[oi]; ni < scratch.overlapNoisyEnd[oi]; ni++) {
                        const auto& nr = scratch.overlapNoisyRegions[ni];
                        if (c.key.pos >= nr.start && c.key.pos < nr.end) {
                            inOverlapNoisy = true; break;
                        }
                    }
                    if (inOverlapNoisy) continue; // pgphase: skip, don't append
                }
                allele = 0;
            }

            // Append to profile.
            if (prof.startVarIdx < 0) prof.startVarIdx = int(ci);
            while (prof.startVarIdx + int(prof.alleles.size()) < int(ci))
                prof.alleles.push_back(-1);
            prof.endVarIdx = int(ci);
            prof.alleles.push_back(allele);
        }
    }
}

// ============================================================================
// Fisher exact test for strand bias (pgphase log_hypergeometric + fisher_exact)
// ============================================================================

/// Log of the hypergeometric PMF term (longcallD `log_hypergeometric`, lgamma only).
/// 2x2 table: row1 (a,b), row2 (c,d).
// kmLogHypergeom and kmFisherExactTwoTail are now inline in PhasingKmeansTypes.hpp.

// ============================================================================
// Step 4: Classify candidates and detect noisy regions
//
// Port of pgphase classify_cand_vars_pgphase (two-pass) +
// post_process_noisy_regs_pgphase + apply_noisy_containment_filter +
// prune_not_candidate_variants.
//
// Backbone bases serve as the reference sequence; sdust low-complexity
// intervals and per-overlap noisy regions (from kmParseCigars) are used
// in place of pgphase's BAM-derived equivalents.
// ============================================================================

/// pgphase classify_variant_initial: local rules only.
/// Returns initial category; does NOT set c.category (caller stores in cats[]).
static KmVariantCategory kmClassifyVariantInitial(
    KmCandidate& c, const uint8_t* bbSeq, uint32_t bbLen,
    const KmPhasingOptions& opts)
{
    c.alleleFraction = c.totalCov > 0 ? double(c.altCov) / double(c.totalCov) : 0.0;

    // pgphase: depth check includes low-quality reads.
    const int depthWithLowQual = c.totalCov + c.lowQualCov;
    if (depthWithLowQual < int(opts.minDepth) || c.altCov < int(opts.minAltDepth))
        return KmVariantCategory::LowCoverage;

    // ONT strand bias: Fisher exact on (fwd_alt, rev_alt, expected, expected).
    if (opts.isOnt) {
        const int fa = c.fwdAlt;
        const int ra = c.revAlt;
        const int expected = (fa + ra) / 2;
        if (expected > 0) {
            const double p = kmFisherExactTwoTail(fa, ra, expected, expected);
            if (p < opts.strandBiasPval)
                return KmVariantCategory::StrandBias;
        }
    }

    if (c.alleleFraction < opts.minAf)
        return KmVariantCategory::LowAlleleFraction;
    if (c.alleleFraction > opts.maxAf)
        return KmVariantCategory::CleanHom;

    // Filter all small indels — nanopore indel noise dominates.
    // Only SNPs and large SVs (>= minSvLen) are used for phasing.
    if (c.key.type == KmVarType::Insertion || c.key.type == KmVarType::Deletion) {
        const int indelLen = (c.key.type == KmVarType::Insertion)
            ? int(c.key.altLen) : int(c.key.refLen);
        if (indelLen >= opts.minSvLen)
            return KmVariantCategory::CleanHetIndel;
        return KmVariantCategory::RepeatHetIndel;
    }
    return KmVariantCategory::CleanHetSnp;
}

/// pgphase category_skipped_for_noisy_flank.
static bool kmCategorySkippedForNoisyFlank(KmVariantCategory c)
{
    return c == KmVariantCategory::LowCoverage ||
           c == KmVariantCategory::StrandBias ||
           c == KmVariantCategory::NonVariant;
}

/// Candidate genomic span on backbone (pgphase variant_genomic_span).
static void kmVariantSpan(const KmVarKey& k, uint32_t& start, uint32_t& end)
{
    if (k.type == KmVarType::Insertion) {
        start = k.pos;
        end = k.pos;  // zero-width on backbone
    } else {
        start = k.pos;
        end = k.pos + k.refLen;
    }
}

static void kmClassifyCandidates(
    uint32_t backboneLen, KmScratchpad& scratch, const KmPhasingOptions& opts)
{
    const uint32_t numCand = uint32_t(scratch.candidates.size());
    if (numCand == 0) return;
    const uint8_t* bbSeq = scratch.backboneBases.data();

    // ── Pass 1: local classification (pgphase classify_variant_initial) ──
    // Fill cats[] and build interval index of all surviving candidates
    // (not LOW_COV, not STRAND_BIAS in ONT mode).
    vector<KmVariantCategory> cats(numCand);
    vector<pair<uint32_t,uint32_t>> varPosIntervals;  // pgphase var_pos_cr

    for (uint32_t ci = 0; ci < numCand; ci++) {
        cats[ci] = kmClassifyVariantInitial(scratch.candidates[ci], bbSeq, backboneLen, opts);

        // Build interval index: skip LOW_COV; skip STRAND_BIAS in ONT mode.
        if (cats[ci] == KmVariantCategory::LowCoverage) continue;
        if (opts.isOnt && cats[ci] == KmVariantCategory::StrandBias) continue;

        uint32_t vs, ve;
        kmVariantSpan(scratch.candidates[ci].key, vs, ve);
        if (ve <= vs) ve = vs + 1;  // insertions: 1bp span for overlap checks
        varPosIntervals.push_back({vs, ve});
    }
    sort(varPosIntervals.begin(), varPosIntervals.end());

    // ── Pass 2: context rules (pgphase classify_cand_vars pass 2) ──
    // chunk_noisy containment → NON_VAR; RepeatHetIndel → noisy seed;
    // dense overlap (n_ov > 1) → noisy seed; LOW_AF → LOW_COV.
    vector<pair<uint32_t,uint32_t>> noisyVarSeeds;  // pgphase noisy_var_cr

    // Build sorted chunk_noisy intervals from pre_process_noisy_regs output.
    vector<pair<uint32_t,uint32_t>> chunkNoisy;
    for (const auto& nr : scratch.noisyRegions)
        chunkNoisy.push_back({nr.start, nr.end});
    sort(chunkNoisy.begin(), chunkNoisy.end());

    for (uint32_t ci = 0; ci < numCand; ci++) {
        KmVariantCategory c = cats[ci];

        // pgphase: skip NON_VAR and STRAND_BIAS before overlap checks.
        if (c == KmVariantCategory::NonVariant || c == KmVariantCategory::StrandBias)
            continue;

        uint32_t vs, ve;
        kmVariantSpan(scratch.candidates[ci].key, vs, ve);
        if (ve <= vs) ve = vs + 1;

        // pgphase: if candidate overlaps chunk_noisy → NON_VAR.
        if (!chunkNoisy.empty() && kmOverlapsAny(chunkNoisy, vs, ve)) {
            cats[ci] = KmVariantCategory::NonVariant;
            continue;
        }

        if (c == KmVariantCategory::LowCoverage)
            continue;

        // RepeatHetIndel → add to noisy seeds (pgphase: cr_add_var_to_noisy_cr, no ratio check).
        // Extend variant span to include overlapping low-complexity intervals.
        if (c == KmVariantCategory::RepeatHetIndel) {
            uint32_t es = vs, ee = ve;
            kmExtendWithLowComp(scratch.lowComplexity, es, ee);
            noisyVarSeeds.push_back({es, ee});
            continue;
        }

        // Dense overlap check: n_ov > 1 in var_pos_cr → add to noisy seeds
        // (pgphase: cr_add_var_to_noisy_cr with check_noisy_reads_ratio=true).
        // Ratio check: of overlaps spanning this position, what fraction have
        // per-overlap noisy intervals here? Only add if ratio >= minAf.
        int nOv = kmCountOverlaps(varPosIntervals, vs, ve);
        if (nOv > 1) {
            int nTotal = 0, nNoisy = 0;
            for (uint32_t oi = 0; oi < uint32_t(scratch.overlaps.size()); oi++) {
                const auto& ov = scratch.overlaps[oi];
                // pgphase: any overlap between read span and variant position.
                if (ov.qs >= ve || ov.qe <= vs) continue;
                nTotal++;
                for (uint32_t ni = scratch.overlapNoisyBegin[oi]; ni < scratch.overlapNoisyEnd[oi]; ni++) {
                    const auto& nr = scratch.overlapNoisyRegions[ni];
                    if (nr.start < ve && nr.end > vs) { nNoisy++; break; }
                }
            }
            if (nTotal > 0 && double(nNoisy) / double(nTotal) >= opts.minAf) {
                uint32_t es = vs, ee = ve;
                kmExtendWithLowComp(scratch.lowComplexity, es, ee);
                noisyVarSeeds.push_back({es, ee});
            }
        }

        // pgphase: LOW_AF → LOW_COV after this loop.
        if (c == KmVariantCategory::LowAlleleFraction)
            cats[ci] = KmVariantCategory::LowCoverage;
    }

    // ── Merge noisy_var_cr seeds into chunk noisy regions ──
    // pgphase: cr_merge2(chunk_noisy, noisy_var_cr, ...) merges pre-existing
    // chunk noisy (from pre_process_noisy_regs) with new seeds from classification.
    if (!noisyVarSeeds.empty()) {
        // Add pre-existing chunk noisy regions to the seed list.
        for (const auto& nr : scratch.noisyRegions)
            noisyVarSeeds.push_back({nr.start, nr.end});
        kmMergeIntervals(noisyVarSeeds, opts.noisyRegMergeDis);
    } else if (!scratch.noisyRegions.empty()) {
        // No new seeds, but keep pre-existing chunk noisy regions.
        for (const auto& nr : scratch.noisyRegions)
            noisyVarSeeds.push_back({nr.start, nr.end});
    }

    // Write merged noisy regions back (pgphase: intervals_from_cr → chunk.noisy_regions).
    scratch.noisyRegions.clear();
    for (const auto& [s, e] : noisyVarSeeds)
        scratch.noisyRegions.push_back({s, e, 0, false});

    // ── Post-process noisy regions: extend boundaries to nearest clean candidates ──
    // (pgphase post_process_noisy_regs_pgphase → collect_noisy_reg_start_end)
    if (!noisyVarSeeds.empty()) {
        const int32_t flankLen = int32_t(opts.noisyRegFlankLen);
        const int nNoisy = int(noisyVarSeeds.size());
        vector<int32_t> startOut(nNoisy), endOut(nNoisy);

        // Find max_left_var_i and min_right_var_i per noisy region.
        vector<int> maxLeftVarI(nNoisy, -1);
        vector<int> minRightVarI(nNoisy, -1);

        int regI = 0, varI = 0;
        while (regI < nNoisy && varI < int(numCand)) {
            if (kmCategorySkippedForNoisyFlank(cats[varI])) { varI++; continue; }
            uint32_t vs, ve;
            kmVariantSpan(scratch.candidates[varI].key, vs, ve);
            uint32_t regStart = noisyVarSeeds[regI].first;
            uint32_t regEnd = noisyVarSeeds[regI].second;
            if (vs > regEnd) {
                if (minRightVarI[regI] == -1) minRightVarI[regI] = varI;
                regI++;
            } else if (ve < regStart) {
                maxLeftVarI[regI] = varI;
                varI++;
            } else {
                varI++;
            }
        }

        for (int ri = 0; ri < nNoisy; ri++) {
            if (maxLeftVarI[ri] == -1) maxLeftVarI[ri] = 0;
            if (minRightVarI[ri] == -1) minRightVarI[ri] = max(0, int(numCand) - 1);

            int32_t regStart = int32_t(noisyVarSeeds[ri].first);
            int32_t regEnd = int32_t(noisyVarSeeds[ri].second);

            // Extend left.
            int32_t curStart = regStart - flankLen;
            for (int v = maxLeftVarI[ri]; v >= 0; v--) {
                if (kmCategorySkippedForNoisyFlank(cats[v])) continue;
                uint32_t vs, ve;
                kmVariantSpan(scratch.candidates[v].key, vs, ve);
                if (int32_t(ve) < curStart - 1) break;
                if (int32_t(vs) - flankLen < curStart)
                    curStart = int32_t(vs) - flankLen;
            }
            startOut[ri] = max(curStart, int32_t(0));

            // Extend right.
            int32_t curEnd = regEnd + flankLen;
            for (int v = minRightVarI[ri]; v < int(numCand); v++) {
                if (kmCategorySkippedForNoisyFlank(cats[v])) continue;
                uint32_t vs, ve;
                kmVariantSpan(scratch.candidates[v].key, vs, ve);
                if (int32_t(vs) > curEnd + 1) break;
                if (int32_t(ve) + flankLen > curEnd)
                    curEnd = int32_t(ve) + flankLen;
            }
            endOut[ri] = min(curEnd, int32_t(backboneLen));
        }

        // Build extended noisy regions and re-merge.
        vector<pair<uint32_t,uint32_t>> extendedNoisy;
        for (int ri = 0; ri < nNoisy; ri++)
            extendedNoisy.push_back({uint32_t(startOut[ri]), uint32_t(endOut[ri])});
        kmMergeIntervals(extendedNoisy, 0);

        // ── Apply noisy containment filter ──
        // (pgphase apply_noisy_containment_filter: candidates fully inside noisy → NON_VAR)
        for (uint32_t ci = 0; ci < numCand; ci++) {
            // pgphase: skip NOT_CAND categories (NON_VAR, LOW_COV, STRAND_BIAS).
            if (cats[ci] == KmVariantCategory::NonVariant ||
                cats[ci] == KmVariantCategory::LowCoverage ||
                cats[ci] == KmVariantCategory::StrandBias)
                continue;
            uint32_t vs, ve;
            kmVariantSpan(scratch.candidates[ci].key, vs, ve);
            if (ve <= vs) ve = vs + 1;
            if (kmIsContained(extendedNoisy, vs, ve))
                cats[ci] = KmVariantCategory::NonVariant;
        }

        // Replace noisy regions with extended+merged result for future MSA.
        scratch.noisyRegions.clear();
        for (const auto& [s, e] : extendedNoisy)
            scratch.noisyRegions.push_back({s, e, 0, false});
    }

    // ── Write final categories and flags ──
    for (uint32_t ci = 0; ci < numCand; ci++) {
        auto& c = scratch.candidates[ci];
        c.category = cats[ci];
        switch (cats[ci]) {
            case KmVariantCategory::CleanHetSnp:    c.categoryFlag = KM_CLEAN_HET_SNP; break;
            case KmVariantCategory::CleanHetIndel:  c.categoryFlag = KM_CLEAN_HET_INDEL; break;
            case KmVariantCategory::CleanHom:       c.categoryFlag = KM_CLEAN_HOM; break;
            case KmVariantCategory::RepeatHetIndel: c.categoryFlag = KM_REP_HET_VAR; break;
            default:                                c.categoryFlag = 0; break;
        }
    }

    // ── Prune non-candidate variants ──
    // (pgphase prune_not_candidate_variants: remove NON_VAR, LOW_COV, STRAND_BIAS)
    {
        vector<KmCandidate> kept;
        kept.reserve(numCand);
        for (auto& c : scratch.candidates) {
            if (c.category == KmVariantCategory::NonVariant ||
                c.category == KmVariantCategory::LowCoverage ||
                c.category == KmVariantCategory::StrandBias)
                continue;
            kept.push_back(std::move(c));
        }
        scratch.candidates = std::move(kept);
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
    // pgphase select_init_var: CleanHetSnp > CleanHetIndel > NoisyCandHet(SNP) > NoisyCandHet(non-hp indel).
    int bSnp = -1, bIndel = -1, bNSnp = -1, bNIndel = -1;
    int bSnpD = 0, bIndelD = 0, bNSnpD = 0, bNIndelD = 0;
    for (int i = 0; i < int(vi.size()); i++) {
        const auto& c = s.candidates[vi[i]];
        if (c.category == KmVariantCategory::CleanHetSnp) {
            if (bSnp < 0 || c.totalCov > bSnpD) { bSnp = i; bSnpD = c.totalCov; }
        } else if (c.category == KmVariantCategory::CleanHetIndel) {
            if (bIndel < 0 || c.totalCov > bIndelD) { bIndel = i; bIndelD = c.totalCov; }
        } else if (c.category == KmVariantCategory::NoisyCandHet) {
            if (c.key.type == KmVarType::Snp) {
                if (bNSnp < 0 || c.totalCov > bNSnpD) { bNSnp = i; bNSnpD = c.totalCov; }
            } else if (!c.isHomopolymerIndel) {
                if (bNIndel < 0 || c.totalCov > bNIndelD) { bNIndel = i; bNIndelD = c.totalCov; }
            }
        }
    }
    if (bSnp >= 0) return bSnp;
    if (bIndel >= 0) return bIndel;
    if (bNSnp >= 0) return bNSnp;
    return bNIndel;
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

static void kmUpdateConsAlle(KmCandidate& c, int hap, bool isOnt = false) {
    if (hap == 0) return;
    const auto& prof = c.hapAlleProfile[hap];
    int mx = 0, ma = -1, total = 0;
    for (int a = 0; a < int(prof.size()); a++) {
        total += prof[a];
        if (prof[a] > mx) { mx = prof[a]; ma = a; }
    }
    // pgphase ONT 67% guard: homopolymer indels need supermajority.
    if (isOnt && c.isHomopolymerIndel && mx < int(total * 0.67))
        ma = -1;
    c.hapConsAlle[hap] = ma;
}

static int kmAssignOverlapHap(KmScratchpad& s, uint32_t oi, uint32_t flags) {
    const auto& prof = s.overlapProfiles[oi];
    if (prof.startVarIdx < 0) return -1;
    int hs[3] = {0,0,0}, nv[3] = {0,0,0};
    for (int ci = prof.startVarIdx; ci <= prof.endVarIdx; ci++) {
        auto& c = s.candidates[ci];
        if ((c.categoryFlag & flags) == 0) continue;
        // pgphase: skip homopolymer indels and NoisyCandHom in hap assignment.
        if (c.isHomopolymerIndel || c.category == KmVariantCategory::NoisyCandHom) continue;
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

/// Phase 1: update profiles AND consensus per overlap (pgphase update_var_hap_profile_cons_alle).
static void kmUpdateProfilesAndCons(KmScratchpad& s, uint32_t oi, int hap, uint32_t flags, bool isOnt) {
    const auto& prof = s.overlapProfiles[oi];
    if (prof.startVarIdx < 0) return;
    for (int ci = prof.startVarIdx; ci <= prof.endVarIdx; ci++) {
        auto& c = s.candidates[ci];
        if ((c.categoryFlag & flags) == 0) continue;
        int a = kmGetAllele(s, oi, uint32_t(ci));
        if (a < 0 || a >= int(c.hapAlleProfile[1].size())) continue;
        if (hap == 0) {
            c.hapAlleProfile[1][a]++; c.hapAlleProfile[2][a]++;
            kmUpdateConsAlle(c, 1, isOnt); kmUpdateConsAlle(c, 2, isOnt);
        } else {
            c.hapAlleProfile[hap][a]++;
            kmUpdateConsAlle(c, hap, isOnt);
        }
    }
}

/// Phase 2: update profiles only, no consensus (pgphase update_var_hap_profile).
static void kmUpdateProfilesOnly(KmScratchpad& s, uint32_t oi, int hap, uint32_t flags) {
    const auto& prof = s.overlapProfiles[oi];
    if (prof.startVarIdx < 0) return;
    for (int ci = prof.startVarIdx; ci <= prof.endVarIdx; ci++) {
        auto& c = s.candidates[ci];
        if ((c.categoryFlag & flags) == 0) continue;
        int a = kmGetAllele(s, oi, uint32_t(ci));
        if (a < 0 || a >= int(c.hapAlleProfile[1].size())) continue;
        if (hap == 0) {
            c.hapAlleProfile[1][a]++; c.hapAlleProfile[2][a]++;
        } else {
            c.hapAlleProfile[hap][a]++;
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
            if (flip == 1) {
                changed = true;
                // pgphase longcallD: swap loop for hap=1..PLOID (=2) does two swaps
                // that cancel out (net identity). We replicate this no-op behavior.
                // The flip only affects `changed` (forces another iteration).
                for (int h = 1; h <= 2; h++) {
                    int tmp = c.hapConsAlle[h];
                    c.hapConsAlle[h] = c.hapConsAlle[3-h];
                    c.hapConsAlle[3-h] = tmp;
                }
            }
        }
        c.phaseSet = ps;
    }
    return changed;
}

// ============================================================================
// Step 6c: K-means main loop
// ============================================================================

void dinara::kmRunKmeans(KmScratchpad& scratch, const KmPhasingOptions& opts, uint32_t flags)
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
    // pgphase: hapAlleProfile[0] is only zeroed on first init (when [1] and [2] are empty).
    for (uint32_t ci : scratch.validVarIdx) {
        auto& c = scratch.candidates[ci];
        int na = kmAlleSlots(c);
        bool firstInit = c.hapAlleProfile[1].empty() && c.hapAlleProfile[2].empty();
        if (firstInit) c.hapAlleProfile[0].assign(na, 0);
        c.hapAlleProfile[1].assign(na, 0);
        c.hapAlleProfile[2].assign(na, 0);
        if (c.category == KmVariantCategory::CleanHom ||
            c.category == KmVariantCategory::NoisyCandHom) { c.hapConsAlle[1] = 1; c.hapConsAlle[2] = 1; }
        else { c.hapConsAlle[1] = -1; c.hapConsAlle[2] = -1; }
        // pgphase get_var_init_max_cov_allele: argmax of allele coverages.
        // ONT guard: homopolymer indels get -1 (no initial consensus).
        int maxCov = 0, maxAlle = -1;
        if (opts.isOnt && c.isHomopolymerIndel) { c.hapConsAlle[0] = -1; continue; }
        int na0 = int(c.alleCovs.empty() ? 2 : c.alleCovs.size());
        for (int a = 0; a < na0; a++) {
            int cov = c.alleCovs.empty() ? (a == 0 ? c.refCov : (a == 1 ? c.altCov : 0))
                                         : c.alleCovs[a];
            if (cov > maxCov) { maxCov = cov; maxAlle = a; }
        }
        c.hapConsAlle[0] = maxAlle;
    }

    // Phase 1: seed from pivot.
    int pivotVi = kmSelectPivot(scratch, scratch.validVarIdx);
    if (pivotVi < 0) return;
    { auto& pc = scratch.candidates[scratch.validVarIdx[pivotVi]];
      pc.hapConsAlle[1] = 0; pc.hapConsAlle[2] = 1; }

    for (auto& ov : scratch.overlaps) { ov.hap = 0; ov.phaseSet = -1; }

    // Phase 1: variant-centric sweep from pivot outward.
    // pgphase: for each variant in sweep order, find overlaps covering it,
    // assign unassigned reads, update profiles+consensus.
    // Sweep order: [pivot, pivot-1, ..., 0, pivot+1, ..., nv-1].
    // HOM variants are skipped in Phase 1.
    const bool isOnt = opts.isOnt;
    const int nv = int(scratch.validVarIdx.size());
    vector<int> sweepOrder(nv);
    sweepOrder[0] = pivotVi;
    for (int vi = pivotVi - 1; vi >= 0; --vi) sweepOrder[pivotVi - vi] = vi;
    for (int vi = pivotVi + 1; vi < nv; ++vi) sweepOrder[vi] = vi;

    for (int idx = 0; idx < nv; idx++) {
        uint32_t ci = scratch.validVarIdx[sweepOrder[idx]];
        const auto& c = scratch.candidates[ci];
        // pgphase: skip HOM in Phase 1 (both CleanHom and NoisyCandHom).
        if (c.category == KmVariantCategory::CleanHom ||
            c.category == KmVariantCategory::NoisyCandHom) continue;

        // Find all overlaps covering this candidate position.
        for (uint32_t oi = 0; oi < numOv; oi++) {
            if (scratch.overlaps[oi].hap != 0) continue; // already assigned
            const auto& prof = scratch.overlapProfiles[oi];
            if (prof.startVarIdx < 0) continue;
            if (int(ci) < prof.startVarIdx || int(ci) > prof.endVarIdx) continue;

            int hap = kmAssignOverlapHap(scratch, oi, flags);
            // pgphase: if hap == -1, seed as hap 1.
            if (hap == -1) hap = 1;
            scratch.overlaps[oi].hap = hap;
            kmUpdateProfilesAndCons(scratch, oi, hap, flags, isOnt);
        }
    }
    // Phase 2: iterative refinement — batch profile update, then batch consensus
    // (pgphase iter_update_var_hap_to_cons_alle: profiles only per read,
    //  then update_var_hap_to_cons_alle for all vars after all reads).
    for (uint32_t iter = 0; iter < opts.maxKmeansIter; iter++) {
        // Phase-set flip detection (pgphase iter_update_var_hap_cons_phase_set).
        bool flipChanged = kmPhaseSetFlip(scratch, scratch.validVarIdx, opts);

        // Save old consensus for convergence check.
        vector<array<int,3>> savedCons(scratch.validVarIdx.size());
        for (size_t vi = 0; vi < scratch.validVarIdx.size(); vi++)
            savedCons[vi] = scratch.candidates[scratch.validVarIdx[vi]].hapConsAlle;

        // Reset profiles (pgphase var_init_hap_to_alle_profile).
        for (uint32_t ci : scratch.validVarIdx) {
            auto& c = scratch.candidates[ci]; int na = kmAlleSlots(c);
            c.hapAlleProfile[0].assign(na, 0);
            c.hapAlleProfile[1].assign(na, 0);
            c.hapAlleProfile[2].assign(na, 0);
        }

        // Assign all overlaps and update profiles only (no consensus).
        // pgphase: if hap == -1, treat as 0 (unassigned) and still update profiles.
        for (uint32_t oi = 0; oi < numOv; oi++) {
            int newH = kmAssignOverlapHap(scratch, oi, flags);
            if (newH == -1) newH = 0;
            scratch.overlaps[oi].hap = newH;
            kmUpdateProfilesOnly(scratch, oi, newH, flags);
        }

        // Batch consensus update after all overlaps (pgphase update_var_hap_to_cons_alle).
        for (uint32_t ci : scratch.validVarIdx) {
            kmUpdateConsAlle(scratch.candidates[ci], 1, isOnt);
            kmUpdateConsAlle(scratch.candidates[ci], 2, isOnt);
        }

        // Convergence check: break when both phase-set flip and consensus unchanged
        // (pgphase: if (c1 == 0 && c2 == 0) break).
        bool consChanged = false;
        for (size_t vi = 0; vi < scratch.validVarIdx.size(); vi++) {
            const auto& cur = scratch.candidates[scratch.validVarIdx[vi]].hapConsAlle;
            if (cur[1] != savedCons[vi][1] || cur[2] != savedCons[vi][2]) { consChanged = true; break; }
        }
        if (!flipChanged && !consChanged) break;
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

    // Phase 4: fill hapAlt/hapRef (pgphase: iterate all candidates, not just valid).
    for (auto& c : scratch.candidates) {
        int c1 = c.hapConsAlle[1];
        int c2 = c.hapConsAlle[2];
        c.hapAlt = 0;
        c.hapRef = 0;
        // pgphase: both -1 → fall back to hapConsAlle[0] (hom_idx).
        if (c1 == -1 && c2 == -1) {
            c1 = c2 = c.hapConsAlle[0];
        }
        if (c1 == -1) c1 = 0; // unknown hap → ref
        if (c2 == -1) c2 = 0;
        const bool h1Alt = (c1 != 0);
        const bool h2Alt = (c2 != 0);
        if (h1Alt && h2Alt)       { c.hapAlt = 3; c.hapRef = 0; }
        else if (h1Alt && !h2Alt) { c.hapAlt = 1; c.hapRef = 2; }
        else if (!h1Alt && h2Alt) { c.hapAlt = 2; c.hapRef = 1; }
        // both ref/unresolved: leave hapAlt/hapRef as 0/0
    }
}

// ============================================================================
// Step 7: Write results
// ============================================================================

void dinara::kmWriteResults(
    Assembler& assembler, ReadId backboneReadId, const KmScratchpad& scratch)
{
    // Check if any overlap was actually phased (hap 1 or 2).
    bool hasPhased = false;
    for (const auto& ov : scratch.overlaps) {
        if (ov.hap == 1 || ov.hap == 2) { hasPhased = true; break; }
    }

    for (const auto& ov : scratch.overlaps) {
        auto& ad = assembler.alignmentData[ov.alignmentId];
        uint8_t matchState;
        if (!hasPhased) {
            // No het sites found — no phasing signal. Leave default (cis=1).
            matchState = 1;
        } else if (ov.hap == 2) {
            matchState = 2; // trans
        } else if (ov.hap == 1) {
            matchState = 1; // cis
        } else {
            matchState = 1; // unassigned — treat as cis
        }
        ad.setHifiasmEcMatchStateFromReadPerspective(backboneReadId, matchState);
    }
}

// ============================================================================
// Step 8: Second-round refinement on cis overlaps
// ============================================================================

/// Iteratively re-cluster cis overlaps to peel off different paralogous copies.
/// Each round recounts alleles from the remaining cis set, reclassifies,
/// and runs k-means. Overlaps that split into hap 2 are peeled off as
/// matchState=3. Repeats until no overlaps split or no het sites remain.
void dinara::kmRefineCis(
    Assembler& assembler, ReadId backboneReadId,
    KmScratchpad& scratch, const KmPhasingOptions& opts,
    uint32_t bbLen, bool dbg)
{
    const uint32_t numOv = uint32_t(scratch.overlaps.size());

    // Track overlap state across refinement rounds.
    // isCis: 1 = in the refinement pool, 0 = excluded (trans or peeled off).
    // Round 1 cis (hap==1) and unassigned (hap==0) are both included —
    // unassigned overlaps lacked phasing signal but could still be from
    // a different paralogous copy.
    // wasTrans: true for overlaps that were trans (hap==2) in round 1.
    vector<uint8_t> isCis(numOv);
    vector<bool> wasTrans(numOv);
    for (uint32_t oi = 0; oi < numOv; oi++) {
        wasTrans[oi] = (scratch.overlaps[oi].hap == 2);
        isCis[oi] = wasTrans[oi] ? 0 : 1;
    }

    const uint32_t maxRefineRounds = 10;

    for (uint32_t round = 0; round < maxRefineRounds; round++) {

        // Recount alleles using only cis overlaps.
        // Reset counts. Backbone counts as ref.
        for (auto& c : scratch.candidates) {
            c.totalCov = 1; c.refCov = 1; c.fwdRef = 1;
            c.altCov = 0; c.fwdAlt = 0; c.revAlt = 0; c.revRef = 0;
        }
        {
            const uint32_t numCand = uint32_t(scratch.candidates.size());
            for (uint32_t oi = 0; oi < numOv; oi++) {
                if (!isCis[oi]) continue;
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
                    while (di < diEnd && kmCompareDigarToCand(scratch.digars[di], c.key, opts.minSvLen) < 0)
                        di++;
                    bool isAlt = (di < diEnd && kmCompareDigarToCand(scratch.digars[di], c.key, opts.minSvLen) == 0);
                    if (isAlt) { c.altCov++; if (ov.isRev == 0) c.fwdAlt++; else c.revAlt++; di++; }
                    else { c.refCov++; if (ov.isRev == 0) c.fwdRef++; else c.revRef++; }
                }
            }
            for (uint32_t ci = 0; ci < numCand; ci++) {
                auto& c = scratch.candidates[ci];
                c.alleleFraction = c.totalCov > 0 ? double(c.altCov) / double(c.totalCov) : 0.0;
            }
        }

        // Reclassify candidates with cis-only counts.
        kmClassifyCandidates(bbLen, scratch, opts);

        uint32_t cleanHet = 0;
        uint32_t cisCount = 0;
        for (const auto& c : scratch.candidates)
            if (c.category == KmVariantCategory::CleanHetSnp ||
                c.category == KmVariantCategory::CleanHetIndel) cleanHet++;
        for (uint32_t oi = 0; oi < numOv; oi++)
            if (isCis[oi]) cisCount++;

        if (dbg)
            cout << "DEBUG refine-cis round " << round
                 << ": cisCount=" << cisCount << " cleanHet=" << cleanHet << endl;

        // No het sites in cis subset — nothing to split.
        if (cleanHet == 0) break;

        // Rebuild profiles for cis-only overlaps.
        kmBuildOverlapProfiles(scratch, opts.minSvLen, &isCis);

        // Reset hap assignments for k-means.
        for (auto& ov : scratch.overlaps) ov.hap = 0;

        // Run k-means on the cis subset.
        kmRunKmeans(scratch, opts, KM_GERMLINE_CLEAN);

        // Count how many cis overlaps split into hap 2.
        uint32_t peeled = 0;
        for (uint32_t oi = 0; oi < numOv; oi++) {
            if (isCis[oi] && scratch.overlaps[oi].hap == 2) {
                isCis[oi] = 0; // peel off
                peeled++;
            }
        }

        if (dbg)
            cout << "DEBUG refine-cis round " << round
                 << ": peeled=" << peeled << endl;

        // Converged — no overlaps split off.
        if (peeled == 0) break;
    }

    // Write final matchStates.
    for (uint32_t oi = 0; oi < numOv; oi++) {
        auto& ad = assembler.alignmentData[scratch.overlaps[oi].alignmentId];
        uint8_t matchState;
        if (wasTrans[oi]) {
            matchState = 2; // trans from round 1
        } else if (isCis[oi]) {
            matchState = 1; // survived all refinement rounds — true cis
        } else {
            matchState = 3; // peeled off during refinement — different copy
        }
        ad.setHifiasmEcMatchStateFromReadPerspective(backboneReadId, matchState);
    }

    if (dbg) {
        int s1 = 0, s2 = 0, s3 = 0;
        for (uint32_t oi = 0; oi < numOv; oi++) {
            const auto& ad = assembler.alignmentData[scratch.overlaps[oi].alignmentId];
            uint8_t ms = ad.getHifiasmEcMatchStateFromReadPerspective(backboneReadId);
            if (ms == 1) s1++;
            else if (ms == 2) s2++;
            else if (ms == 3) s3++;
        }
        cout << "DEBUG refine-cis: FINAL cis=" << s1
             << " trans=" << s2 << " cisDiffCopy=" << s3 << endl;
    }
}

// ============================================================================
// Public entry point
// ============================================================================

void Assembler::phaseOverlapsKmeans(uint64_t threadCount, bool isOnt)
{
    cout << timestamp << "=== K-means Overlap Phasing ===" << endl;
    cout << timestamp << "Using OverlapCigarStore (raw CIGAR parsing)." << endl;
    const uint64_t readCount = getReads().readCount();
    cout << timestamp << "Read count: " << readCount << endl;
    if (readCount == 0) { cout << timestamp << "No reads." << endl; return; }

    // Debug mode: if DINARA_PHASING_DEBUG_READ is set, process only that read
    // with verbose output, single-threaded.
    int64_t debugReadId = -1;
    if (const char* env = getenv("DINARA_PHASING_DEBUG_READ")) {
        debugReadId = atol(env);
        cout << timestamp << "DEBUG MODE: processing only read " << debugReadId << endl;
        threadCount = 1;
    }

    KmPhasingOptions opts;
    opts.isOnt = isOnt;
    if (isOnt) cout << timestamp << "ONT mode: Fisher exact strand bias filter enabled (p < "
                    << opts.strandBiasPval << ")" << endl;
    atomic<uint64_t> readsProcessed(0), readsWithOverlaps(0), readsWithSites(0);
    atomic<uint64_t> totalCis(0), totalTrans(0), totalCisDiffCopy(0), totalNoisyRegions(0);

    struct alignas(64) TT {
        int64_t gather=0, unpack=0, detect=0, count=0;
        int64_t classify=0, kmeans=0, noisyMsa=0, write=0, refine=0;
        // Sub-timers for profiling.
        int64_t unpackSeq=0, unpackSdust=0;
        int64_t countCollect=0, countAlleles=0, countProfiles=0;
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
                if (debugReadId >= 0 && int64_t(rid) != debugReadId) {
                    readsProcessed++; continue;
                }
                ReadId readId(rid);
                scratch.clear();

                const bool dbg = (debugReadId >= 0);
                auto t0 = clk::now();
                kmGatherOverlaps(*this, readId, scratch);
                auto t1 = clk::now(); tt.gather += us(t0,t1);
                if (scratch.overlaps.empty()) {
                    if (dbg) cout << "DEBUG read " << rid << ": no overlaps found." << endl;
                    readsProcessed++; continue;
                }
                readsWithOverlaps++;
                if (dbg) cout << "DEBUG read " << rid << ": " << scratch.overlaps.size() << " overlaps" << endl;

                t0 = clk::now();
                uint32_t bbLen = uint32_t(getReads().getRead(readId).baseCount);
                scratch.backboneBases.resize(bbLen);
                auto seq = getReads().getRead(readId);
                for (uint32_t i = 0; i < bbLen; i++) scratch.backboneBases[i] = seq[i].value;
                auto tMid = clk::now(); tt.unpackSeq += us(t0,tMid);
                // pgphase populate_low_complexity_intervals: sdust on backbone.
                // Disabled until Step 4 (noisy region MSA) is implemented.
                // SDUST is only used to extend noisy region boundaries; without
                // it, noisy regions are slightly smaller but phasing is unaffected.
                // kmSdust(scratch.backboneBases.data(), int(bbLen),
                //         int(opts.sdustThreshold), int(opts.sdustWindow),
                //         scratch.lowComplexity);
                scratch.lowComplexity.clear();
                t1 = clk::now(); tt.unpackSdust += us(tMid,t1);
                tt.unpack += us(t0,t1);
                if (dbg) cout << "DEBUG read " << rid << ": bbLen=" << bbLen
                              << " lowComplexity=" << scratch.lowComplexity.size() << " intervals" << endl;

                t0 = clk::now();
                kmParseCigars(*this, readId, bbLen, scratch, opts);
                t1 = clk::now(); tt.detect += us(t0,t1);
                if (dbg) cout << "DEBUG read " << rid << ": " << scratch.digars.size() << " digars, "
                              << scratch.overlapNoisyRegions.size() << " per-overlap noisy regions" << endl;

                t0 = clk::now();
                kmCollectCandidates(scratch, opts.minSvLen);
                tMid = clk::now(); tt.countCollect += us(t0,tMid);
                kmCountAlleles(scratch, opts.minSvLen);
                t1 = clk::now(); tt.countAlleles += us(tMid,t1);
                tt.count += us(t0,t1);
                if (dbg) cout << "DEBUG read " << rid << ": " << scratch.candidates.size() << " candidates after collect+count" << endl;

                // pgphase step 2.1: pre-process per-overlap noisy regions before classification.
                t0 = clk::now();
                kmPreProcessNoisyRegs(scratch, opts);
                kmClassifyCandidates(bbLen, scratch, opts);
                t1 = clk::now(); tt.classify += us(t0,t1);
                totalNoisyRegions += scratch.noisyRegions.size();

                uint32_t cleanHet = 0, cleanHom = 0, repHet = 0, noisyHet = 0;
                for (const auto& c : scratch.candidates) {
                    if (c.category == KmVariantCategory::CleanHetSnp ||
                        c.category == KmVariantCategory::CleanHetIndel) cleanHet++;
                    else if (c.category == KmVariantCategory::CleanHom) cleanHom++;
                    else if (c.category == KmVariantCategory::RepeatHetIndel) repHet++;
                    else if (c.category == KmVariantCategory::NoisyCandHet) noisyHet++;
                }
                if (cleanHet > 0) readsWithSites++;
                if (dbg) {
                    cout << "DEBUG read " << rid << ": after classify: "
                         << scratch.candidates.size() << " candidates (cleanHet=" << cleanHet
                         << " cleanHom=" << cleanHom << " repHetIndel=" << repHet
                         << " noisyHet=" << noisyHet << "), "
                         << scratch.noisyRegions.size() << " noisy regions" << endl;
                    // Print each candidate
                    for (uint32_t ci = 0; ci < uint32_t(scratch.candidates.size()); ci++) {
                        const auto& c = scratch.candidates[ci];
                        const char* catStr = "?";
                        switch (c.category) {
                            case KmVariantCategory::CleanHetSnp:    catStr = "CleanHetSnp"; break;
                            case KmVariantCategory::CleanHetIndel:  catStr = "CleanHetIndel"; break;
                            case KmVariantCategory::CleanHom:       catStr = "CleanHom"; break;
                            case KmVariantCategory::LowCoverage:    catStr = "LowCov"; break;
                            case KmVariantCategory::LowAlleleFraction: catStr = "LowAF"; break;
                            case KmVariantCategory::StrandBias:     catStr = "StrandBias"; break;
                            case KmVariantCategory::RepeatHetIndel: catStr = "RepHetIndel"; break;
                            case KmVariantCategory::NoisyCandHet:   catStr = "NoisyCandHet"; break;
                            case KmVariantCategory::NoisyCandHom:   catStr = "NoisyCandHom"; break;
                            case KmVariantCategory::NonVariant:     catStr = "NonVar"; break;
                        }
                        const char* typeStr = c.key.type == KmVarType::Snp ? "SNP" :
                                              c.key.type == KmVarType::Insertion ? "INS" : "DEL";
                        cout << "  cand[" << ci << "] pos=" << c.key.pos << " " << typeStr
                             << " refLen=" << c.key.refLen << " altLen=" << c.key.altLen
                             << " cov=" << c.totalCov << " alt=" << c.altCov
                             << " af=" << c.alleleFraction << " cat=" << catStr << endl;
                    }
                }

                // Build profiles after classification (skips NON_VAR).
                t0 = clk::now();
                kmBuildOverlapProfiles(scratch, opts.minSvLen);
                t1 = clk::now(); tt.countProfiles += us(t0,t1);
                tt.count += us(t0,t1);
                if (dbg) {
                    int profiled = 0;
                    for (const auto& p : scratch.overlapProfiles)
                        if (p.startVarIdx >= 0) profiled++;
                    cout << "DEBUG read " << rid << ": " << profiled << "/" << scratch.overlaps.size()
                         << " overlaps have profiles" << endl;
                }

                t0 = clk::now();
                if (cleanHet > 0) kmRunKmeans(scratch, opts, KM_GERMLINE_CLEAN);
                t1 = clk::now(); tt.kmeans += us(t0,t1);

                t0 = clk::now();
                kmWriteResults(*this, readId, scratch);
                t1 = clk::now(); tt.write += us(t0,t1);

                // Step 8: second-round refinement on cis overlaps.
                // Detects different paralogous copies within the cis set.
                t0 = clk::now();
                if (cleanHet > 0)
                    kmRefineCis(*this, readId, scratch, opts, bbLen, dbg);
                t1 = clk::now(); tt.refine += us(t0,t1);

                if (dbg) {
                    int hap1 = 0, hap2 = 0, hap0 = 0;
                    for (const auto& ov : scratch.overlaps) {
                        if (ov.hap == 1) hap1++;
                        else if (ov.hap == 2) hap2++;
                        else hap0++;
                    }
                    cout << "DEBUG read " << rid << ": RESULT hap1(cis)=" << hap1
                         << " hap2(trans)=" << hap2 << " unassigned=" << hap0 << endl;
                }

                // Count final matchStates by reading back from alignmentData.
                for (const auto& ov : scratch.overlaps) {
                    const auto& ad = alignmentData[ov.alignmentId];
                    uint8_t ms = ad.getHifiasmEcMatchStateFromReadPerspective(readId);
                    if (ms == 1) totalCis++;
                    else if (ms == 2) totalTrans++;
                    else if (ms == 3) totalCisDiffCopy++;
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
        total.kmeans += tt.kmeans; total.noisyMsa += tt.noisyMsa;
        total.write += tt.write; total.refine += tt.refine;
        total.unpackSeq += tt.unpackSeq; total.unpackSdust += tt.unpackSdust;
        total.countCollect += tt.countCollect; total.countAlleles += tt.countAlleles;
        total.countProfiles += tt.countProfiles;
    }
    auto ms = [](int64_t u) { return u/1000; };
    cout << timestamp << "Timing (ms, sum " << threadCount << " threads):" << endl;
    cout << timestamp << "  gather:   " << ms(total.gather) << endl;
    cout << timestamp << "  unpack:   " << ms(total.unpack)
         << "  (seq=" << ms(total.unpackSeq) << " sdust=" << ms(total.unpackSdust) << ")" << endl;
    cout << timestamp << "  detect:   " << ms(total.detect) << endl;
    cout << timestamp << "  count+prof:" << ms(total.count)
         << "  (collect=" << ms(total.countCollect) << " alleles=" << ms(total.countAlleles)
         << " profiles=" << ms(total.countProfiles) << ")" << endl;
    cout << timestamp << "  classify: " << ms(total.classify) << endl;
    cout << timestamp << "  kmeans:   " << ms(total.kmeans) << endl;
    cout << timestamp << "  noisyMsa: " << ms(total.noisyMsa) << endl;
    cout << timestamp << "  write:    " << ms(total.write) << endl;
    cout << timestamp << "  refine:   " << ms(total.refine) << endl;
    cout << timestamp << "Complete. Reads=" << readsProcessed.load()
         << " withOvlp=" << readsWithOverlaps.load()
         << " withSites=" << readsWithSites.load()
         << " cis=" << totalCis.load()
         << " trans=" << totalTrans.load()
         << " cisDiffCopy=" << totalCisDiffCopy.load()
         << " noisyRegs=" << totalNoisyRegions.load() << endl;
}
