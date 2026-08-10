/// @file AssemblerWindowKsw2LeafSnarls.cpp
/// @brief Leaf-snarl (multi-column MNP/het-block) detection per anchor window,
///        from INDEPENDENT pairwise ksw2 alignments -- no shared multi-sequence
///        graph, no per-interval fragmentation.
///
/// For each member read sharing anchors with a window's backbone, this builds
/// ONE dinara::KwMemberProfile (WindowHetProfiles.hpp -- the SAME engine-
/// agnostic per-member representation emitHetBubblesFromProfiles and
/// findLeafSnarlsFromProfiles already consume) by running one banded ksw2
/// global alignment per inter-anchor segment of THAT MEMBER ONLY against the
/// backbone, FULLSPAN (every column the member spans gets a real aligned
/// entry, anchor bodies included -- not just inter-anchor gaps). Each
/// member's alignment is entirely independent: nothing accumulates across
/// members, so cost/memory for member i is bounded purely by member i's own
/// overlap span, regardless of how many other members were processed before
/// it or how many intervals the window has. This is deliberately a separate,
/// self-contained copy of the per-segment ksw2 alignment logic in
/// AssemblerWindowKsw2HetSites.cpp (NOT a refactor of it): that file declares
/// its OWN local (anonymous-namespace) KwMemberProfile/KwAlignedCol structs
/// with extra fields (readCBegin/pinBb/insertions/...) for its two-pass
/// consensus feature, which collide by name with the shared
/// dinara::KwMemberProfile from WindowHetProfiles.hpp this function needs --
/// #including that file's helpers here would silently shadow the wrong type.
/// Duplicating the (short, already-proven) per-segment alignment loop avoids
/// that collision entirely and keeps the existing, working ksw2 het-bubble
/// path untouched.
///
/// Verification only: writes nothing to AnchorWindow::hetBubbles or anywhere
/// else. Purely a comparison point against the interval-abPOA-profile and
/// abPOA-shared-graph leaf-snarl detectors on the same windows.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "AssemblerOptions.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "WindowIntervalPoa.hpp"
#include "invalid.hpp"
#include "timestamp.hpp"

extern "C" {
#include "ksw2.h"
}

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace std;
using namespace dinara;

namespace {

// Extra band slack on top of the pin-implied length-diff offset. Mirrors
// AssemblerWindowKsw2HetSites.cpp's kwBandSlack: inter-anchor segments are
// short and nearly diagonal, so a small constant absorbs indel jitter --
// PROVIDED the segment itself is short. A fixed slack is not enough once a
// segment spans thousands of bases (LIS-filtering can legitimately drop many
// intermediate anchors in a marker-sparse stretch, leaving one big gap
// between the surviving pins): several small indels scattered through a long
// segment can each nudge the true alignment off the direct end-to-end
// diagonal, and their CUMULATIVE mid-segment drift can exceed a small fixed
// band even when the segment's own net length difference is small. Once the
// true path falls outside the band, ksw2's banded DP is forced onto a
// suboptimal path for the rest of the segment, producing a cascade of
// artificial "mismatches" that are actually a banding failure, not real
// divergence. Confirmed directly: a 2165bp segment (one member's own segment
// this large, from a real marker-sparse stretch) analyzed with a fixed
// 30bp-slack band showed ~20-40% "mismatch" rates -- an order of magnitude
// above any plausible HiFi/het divergence rate -- which then vanished once
// the band scaled with segment length (see lsnComputeBand). Scaling by a
// small fraction of segment length (not just its net length difference)
// gives the DP room for indel density to accumulate across a long segment,
// while staying a tiny, fixed 30bp for the common short (~50-100bp) case.
constexpr int lsnBandSlack = 30;
constexpr double lsnBandFraction = 0.05;

inline int lsnComputeBand(uint32_t bbSegLen, uint32_t cSegLen) {
    const int lenDiff = abs(int(bbSegLen) - int(cSegLen));
    const int scaled = int(lsnBandFraction * double(max(bbSegLen, cSegLen)));
    return lenDiff + max(lsnBandSlack, scaled);
}

struct LsnPin { uint32_t bbPos; uint32_t cPos; };

// Keep the longest colinear (strictly diagonal-increasing) run of pins: sort
// by backbone position, then take the LIS over read positions. Drops
// repeat-induced off-diagonal pins that would otherwise create bogus
// inter-anchor segments. Identical in intent to lisByCPos in
// AssemblerWindowKsw2HetSites.cpp.
void lsnLisByCPos(vector<LsnPin>& pins) {
    const uint32_t n = uint32_t(pins.size());
    if (n < 2) return;
    vector<uint32_t> tails;
    vector<int32_t> pred(n, -1);
    for (uint32_t i = 0; i < n; i++) {
        uint32_t lo = 0, hi = uint32_t(tails.size());
        while (lo < hi) {
            const uint32_t mid = (lo + hi) / 2;
            if (pins[tails[mid]].cPos < pins[i].cPos) lo = mid + 1;
            else hi = mid;
        }
        if (lo > 0) pred[i] = int32_t(tails[lo - 1]);
        if (lo == uint32_t(tails.size())) tails.push_back(i);
        else tails[lo] = i;
    }
    vector<LsnPin> kept;
    for (int32_t c = int32_t(tails.back()); c >= 0; c = pred[uint32_t(c)])
        kept.push_back(pins[uint32_t(c)]);
    reverse(kept.begin(), kept.end());
    pins.swap(kept);
}

} // anonymous namespace

// Noise tracking now uses the shared IpoaNoiseTracker (WindowIntervalPoa.hpp),
// consolidated there once AssemblerWindowProjectedAlignmentLeafSnarls.cpp
// needed the identical logic -- was a per-file duplicate (LsnNoiseTracker)
// to avoid the KwMemberProfile naming collision described in this file's
// header; that collision only concerned KwMemberProfile/KwAlignedCol, not
// this self-contained tracker, so sharing it is safe.

void Assembler::computeWindowKsw2LeafSnarls(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const AlignOptions& alignOptions,
    uint64_t threadCount) const
{
    if (!markers) {
        cout << timestamp << "computeWindowKsw2LeafSnarls: markers not available. Skipping." << endl;
        return;
    }
    const Reads& rds = getReads();
    const auto& mkrs = *markers;
    const uint64_t k = assemblerInfo->k;

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
        if (threadCount == 0) threadCount = 1;
    }

    cout << timestamp << "computeWindowKsw2LeafSnarls: processing "
         << anchorWindows.size() << " windows with " << threadCount << " threads." << endl;

    // overlapDpMismatchScore is already stored as a negative addend (default
    // -4, symmetric with overlapDpMatchScore's positive +2 -- both are meant
    // to be summed directly into a score, matching how AlignmentInfo::dpScore
    // uses them). ksw2's substitution matrix (mat[]) also expects a SIGNED
    // addend at each entry (positive for match, negative for mismatch, added
    // directly by the DP), not a positive penalty magnitude to subtract.
    // AssemblerWindowKsw2HetSites.cpp negates this value before use
    // (probably copying the sign convention of the gap-open/extend
    // parameters below, which genuinely ARE positive magnitudes ksw2
    // subtracts internally) -- but for the mismatch score that double-negates
    // an already-negative value, turning a -4 penalty into a +4 REWARD for
    // mismatches. Confirmed empirically: using that same formula here
    // produced a 19% aggregate mismatch rate (443/491 members >5%) on real
    // HiFi data, which dropped to a realistic rate once this negation was
    // removed. Do NOT copy the negated form from that file without fixing it
    // there too -- see the note left for the user about this being a
    // pre-existing bug in already-reachable code (DINARA_HET_ENGINE=ksw2).
    const int8_t kswMatch = int8_t(alignOptions.overlapDpMatchScore);
    const int8_t kswMismatch = int8_t(alignOptions.overlapDpMismatchScore);
    const int8_t kswGapO1 = int8_t(alignOptions.overlapDpGapOpen1);
    const int8_t kswGapE1 = int8_t(alignOptions.overlapDpGapExtend1);
    const int8_t kswGapO2 = int8_t(alignOptions.overlapDpGapOpen2);
    const int8_t kswGapE2 = int8_t(alignOptions.overlapDpGapExtend2);
    int8_t mat[25];
    for (int a = 0; a < 5; a++)
        for (int b = 0; b < 5; b++)
            mat[a * 5 + b] = (a == b && a < 4) ? kswMatch : kswMismatch;

    std::atomic<uint64_t> nextWindow{0};
    std::atomic<uint64_t> windowsProcessed{0};
    std::atomic<uint64_t> membersAligned{0};
    std::atomic<uint64_t> snarlsFound{0};

    auto worker = [&]() {
        for (;;) {
            const uint64_t wi = nextWindow.fetch_add(1);
            if (wi >= anchorWindows.size()) break;
            const AnchorWindow& window = anchorWindows[wi];

            const OrientedReadId bbOid = window.backboneOrientedReadId;
            const ReadId bbReadId = bbOid.getReadId();
            const uint32_t bbLen = uint32_t(rds.getRead(bbReadId).baseCount);
            const auto bbJ = journeys[bbOid];

            if (window.backboneEnd <= window.backboneBegin + 1) continue;
            const uint32_t firstAnchorJP = window.backboneBegin;
            const uint32_t lastAnchorJP = window.backboneEnd - 1;
            if (firstAnchorJP >= bbJ.size() || lastAnchorJP >= bbJ.size()) continue;

            const uint32_t firstOrd = anchors.getOrdinal(bbJ[firstAnchorJP], bbOid);
            const uint32_t lastOrd = anchors.getOrdinal(bbJ[lastAnchorJP], bbOid);
            if (firstOrd == invalid<uint32_t> || lastOrd == invalid<uint32_t>) continue;

            const uint32_t windowBbBegin = mkrs[bbOid.getValue()][firstOrd].position;
            const uint32_t windowBbEnd = mkrs[bbOid.getValue()][lastOrd].position + uint32_t(k);
            if (windowBbEnd <= windowBbBegin) continue;

            vector<uint8_t> bbSeqVec(bbLen);
            for (uint32_t i = 0; i < bbLen; i++)
                bbSeqVec[i] = rds.getOrientedReadBase(bbOid, i).value;

            // Live per-read pins, exactly like ksw2DetectHetBubblesInWindow:
            // every anchor a read shares with the backbone inside this window
            // becomes a pin, read from the anchor's own live member list (not
            // the persisted/thinned window.readIntervals sharedPins subset).
            unordered_map<uint64_t, vector<LsnPin>> readPins;
            for (uint32_t jp = window.backboneBegin; jp < window.backboneEnd; jp++) {
                if (jp >= bbJ.size()) break;
                const Shasta2AnchorId aid = bbJ[jp];
                const uint32_t bbOrd = anchors.getOrdinal(aid, bbOid);
                if (bbOrd == invalid<uint32_t> || bbOrd >= mkrs[bbOid.getValue()].size()) continue;
                const uint32_t bbPos = mkrs[bbOid.getValue()][bbOrd].position;
                if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;

                const auto anchor = anchors[aid];
                for (const auto& info : anchor) {
                    const OrientedReadId cOid = info.orientedReadId;
                    if (cOid == bbOid) continue;
                    if (cOid.getValue() >= journeys.size()) continue;
                    const uint32_t cOrd = info.ordinal;
                    if (cOrd >= mkrs[cOid.getValue()].size()) continue;
                    const uint32_t cPos = mkrs[cOid.getValue()][cOrd].position;
                    readPins[cOid.getValue()].push_back(LsnPin{bbPos, cPos});
                }
            }

            vector<KwMemberProfile> profiles;
            profiles.reserve(readPins.size());

            for (auto& [orientedValue, pinsRaw] : readPins) {
                const OrientedReadId cOid = OrientedReadId::fromValue(ReadId(orientedValue));
                if (pinsRaw.size() < 2) continue;

                sort(pinsRaw.begin(), pinsRaw.end(),
                    [](const LsnPin& a, const LsnPin& b) {
                        return a.bbPos < b.bbPos || (a.bbPos == b.bbPos && a.cPos < b.cPos);
                    });
                lsnLisByCPos(pinsRaw);
                if (pinsRaw.size() < 2) continue;

                KwMemberProfile prof;
                prof.oid = cOid;
                prof.bbCovBegin = max(pinsRaw.front().bbPos, windowBbBegin);
                prof.bbCovEnd = min(pinsRaw.back().bbPos + uint32_t(k), windowBbEnd);

                // Noise gate: exclude this member's votes inside locally dense
                // mismatch/indel clusters (see the design note on
                // IpoaNoiseTracker, WindowIntervalPoa.hpp, for why this is
                // required, not optional). Same defaults as
                // ksw2DetectHetBubblesInWindow's HiFi tuning
                // (noisyRegSlideWin=100, noisyRegMaxXgaps=5).
                IpoaNoiseTracker noise(100, 5, prof.noisyRanges);

                const char* seqDebugOid = getenv("DINARA_KSW2_LEAFSNARL_SEQ_DEBUG_OID");
                const bool seqDebugThis = seqDebugOid != nullptr &&
                    cOid == OrientedReadId(string(seqDebugOid));

                // FULLSPAN tiling: segment i covers [pins[i].bbPos, pins[i+1].bbPos)
                // (anchor body included at the segment head), plus a final tail
                // [lastPin.bbPos, lastPin.bbPos+k) so the last anchor body is
                // covered too. Every column the member spans gets a real
                // aligned entry -- no columns silently assumed ref.
                for (size_t pi = 0; pi + 1 < pinsRaw.size(); pi++) {
                    const LsnPin& left = pinsRaw[pi];
                    const LsnPin& right = pinsRaw[pi + 1];
                    const uint32_t bbSegBegin = left.bbPos;
                    const uint32_t bbSegEnd = right.bbPos;
                    const uint32_t cSegBegin = left.cPos;
                    const uint32_t cSegEnd = right.cPos;
                    if (bbSegEnd <= bbSegBegin || cSegEnd <= cSegBegin) continue;
                    if (bbSegEnd <= windowBbBegin || bbSegBegin >= windowBbEnd) continue;

                    const uint32_t bbSegLen = bbSegEnd - bbSegBegin;
                    const uint32_t cSegLen = cSegEnd - cSegBegin;

                    static thread_local vector<uint8_t> query;
                    static thread_local vector<uint8_t> target;
                    query.resize(cSegLen);
                    target.resize(bbSegLen);
                    for (uint32_t i = 0; i < cSegLen; i++)
                        query[i] = rds.getOrientedReadBase(cOid, cSegBegin + i).value;
                    for (uint32_t i = 0; i < bbSegLen; i++)
                        target[i] = rds.getOrientedReadBase(bbOid, bbSegBegin + i).value;

                    const int band = lsnComputeBand(bbSegLen, cSegLen);

                    if (seqDebugThis) {
                        auto dump = [](const vector<uint8_t>& v) {
                            string s; s.reserve(v.size());
                            for (uint8_t b : v) s += "ACGTN"[b < 4 ? b : 4];
                            return s;
                        };
                        cout << "      seqProbe window=" << window.windowId << " oid=" << cOid << " pi=" << pi
                             << " bbSeg=[" << bbSegBegin << "," << bbSegEnd << ") len=" << bbSegLen
                             << " cSeg=[" << cSegBegin << "," << cSegEnd << ") len=" << cSegLen
                             << " band=" << band << "\n"
                             << "        target=" << dump(target) << "\n"
                             << "        query =" << dump(query) << endl;
                    }

                    ksw_extz_t ez;
                    memset(&ez, 0, sizeof(ez));
                    ksw_extd2_sse(
                        nullptr,
                        int(cSegLen), query.data(),
                        int(bbSegLen), target.data(),
                        5, mat,
                        kswGapO1, kswGapE1, kswGapO2, kswGapE2,
                        band, -1, 0, 0, &ez);

                    if (ez.n_cigar > 0 && ez.cigar != nullptr) {
                        uint32_t qpos = cSegBegin;
                        uint32_t tpos = bbSegBegin;
                        for (int ci = 0; ci < ez.n_cigar; ci++) {
                            const uint32_t len = ez.cigar[ci] >> 4;
                            const uint32_t op = ez.cigar[ci] & 0xf;
                            if (op == 0) { // M
                                for (uint32_t j = 0; j < len; j++) {
                                    const uint32_t bbPos = tpos + j;
                                    if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;
                                    const uint32_t rp = qpos + j;
                                    const uint8_t cb = rds.getOrientedReadBase(cOid, rp).value;
                                    prof.alignedCols.push_back(KwAlignedCol{bbPos, rp, cb});
                                    const uint8_t bb = bbSeqVec[bbPos];
                                    if (cb != bb && cb < 4 && bb < 4) {
                                        prof.snps.push_back(KwSnp{bbPos, cb});
                                        noise.observe(bbPos, 1, 1);
                                    }
                                }
                                qpos += len; tpos += len;
                            } else if (op == 1) { // I
                                noise.observe(tpos, 0, int(len));
                                qpos += len;
                            } else if (op == 2) { // D
                                noise.observe(tpos, len, int(len));
                                const uint32_t db = max(tpos, windowBbBegin);
                                const uint32_t de = min(tpos + len, windowBbEnd);
                                if (db < de) prof.deletionRanges.push_back({db, de});
                                tpos += len;
                            } else {
                                tpos += len;
                            }
                        }
                    }
                    if (ez.cigar) free(ez.cigar);
                }
                // Trailing anchor body [lastPin, lastPin+k).
                {
                    const LsnPin& last = pinsRaw.back();
                    const uint32_t bbSegBegin = last.bbPos;
                    const uint32_t bbSegEnd = min(last.bbPos + uint32_t(k), windowBbEnd);
                    const uint32_t cSegBegin = last.cPos;
                    const uint32_t cSegEnd = last.cPos + uint32_t(k);
                    if (bbSegEnd > bbSegBegin && cSegEnd > cSegBegin && bbSegBegin < windowBbEnd) {
                        const uint32_t bbSegLen = bbSegEnd - bbSegBegin;
                        const uint32_t cSegLen = cSegEnd - cSegBegin;
                        static thread_local vector<uint8_t> query;
                        static thread_local vector<uint8_t> target;
                        query.resize(cSegLen);
                        target.resize(bbSegLen);
                        for (uint32_t i = 0; i < cSegLen; i++)
                            query[i] = rds.getOrientedReadBase(cOid, cSegBegin + i).value;
                        for (uint32_t i = 0; i < bbSegLen; i++)
                            target[i] = rds.getOrientedReadBase(bbOid, bbSegBegin + i).value;
                        const int band = lsnComputeBand(bbSegLen, cSegLen);
                        if (seqDebugThis) {
                            auto dump = [](const vector<uint8_t>& v) {
                                string s; s.reserve(v.size());
                                for (uint8_t b : v) s += "ACGTN"[b < 4 ? b : 4];
                                return s;
                            };
                            cout << "      seqProbe window=" << window.windowId << " oid=" << cOid
                                 << " TRAILING"
                                 << " bbSeg=[" << bbSegBegin << "," << bbSegEnd << ") len=" << bbSegLen
                                 << " cSeg=[" << cSegBegin << "," << cSegEnd << ") len=" << cSegLen
                                 << " band=" << band << "\n"
                                 << "        target=" << dump(target) << "\n"
                                 << "        query =" << dump(query) << endl;
                        }
                        ksw_extz_t ez;
                        memset(&ez, 0, sizeof(ez));
                        ksw_extd2_sse(
                            nullptr, int(cSegLen), query.data(), int(bbSegLen), target.data(),
                            5, mat, kswGapO1, kswGapE1, kswGapO2, kswGapE2, band, -1, 0, 0, &ez);
                        if (ez.n_cigar > 0 && ez.cigar != nullptr) {
                            uint32_t qpos = cSegBegin, tpos = bbSegBegin;
                            for (int ci = 0; ci < ez.n_cigar; ci++) {
                                const uint32_t len = ez.cigar[ci] >> 4;
                                const uint32_t op = ez.cigar[ci] & 0xf;
                                if (op == 0) {
                                    for (uint32_t j = 0; j < len; j++) {
                                        const uint32_t bbPos = tpos + j;
                                        if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;
                                        const uint32_t rp = qpos + j;
                                        const uint8_t cb = rds.getOrientedReadBase(cOid, rp).value;
                                        prof.alignedCols.push_back(KwAlignedCol{bbPos, rp, cb});
                                        const uint8_t bb = bbSeqVec[bbPos];
                                        if (cb != bb && cb < 4 && bb < 4) {
                                            prof.snps.push_back(KwSnp{bbPos, cb});
                                            noise.observe(bbPos, 1, 1);
                                        }
                                    }
                                    qpos += len; tpos += len;
                                } else if (op == 1) {
                                    noise.observe(tpos, 0, int(len));
                                    qpos += len;
                                } else if (op == 2) {
                                    noise.observe(tpos, len, int(len));
                                    const uint32_t db = max(tpos, windowBbBegin);
                                    const uint32_t de = min(tpos + len, windowBbEnd);
                                    if (db < de) prof.deletionRanges.push_back({db, de});
                                    tpos += len;
                                } else {
                                    tpos += len;
                                }
                            }
                        }
                        if (ez.cigar) free(ez.cigar);
                    }
                }
                noise.finish();

                sort(prof.alignedCols.begin(), prof.alignedCols.end(),
                    [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos < b.bbPos; });
                prof.alignedCols.erase(
                    unique(prof.alignedCols.begin(), prof.alignedCols.end(),
                        [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos == b.bbPos; }),
                    prof.alignedCols.end());
                if (!prof.noisyRanges.empty()) {
                    sort(prof.noisyRanges.begin(), prof.noisyRanges.end());
                    vector<pair<uint32_t, uint32_t>> merged;
                    merged.reserve(prof.noisyRanges.size());
                    for (const auto& r : prof.noisyRanges) {
                        if (!merged.empty() && r.first <= merged.back().second)
                            merged.back().second = max(merged.back().second, r.second);
                        else merged.push_back(r);
                    }
                    prof.noisyRanges = std::move(merged);
                }
                if (!prof.deletionRanges.empty()) {
                    sort(prof.deletionRanges.begin(), prof.deletionRanges.end());
                    vector<pair<uint32_t, uint32_t>> merged;
                    merged.reserve(prof.deletionRanges.size());
                    for (const auto& r : prof.deletionRanges) {
                        if (!merged.empty() && r.first <= merged.back().second)
                            merged.back().second = max(merged.back().second, r.second);
                        else merged.push_back(r);
                    }
                    prof.deletionRanges = std::move(merged);
                }

                if (getenv("DINARA_KSW2_LEAFSNARL_NOISE_DEBUG") != nullptr) {
                    uint64_t noisyBases = 0;
                    for (const auto& r : prof.noisyRanges) noisyBases += (r.second - r.first);
                    cout << "      noiseProbe window=" << window.windowId << " oid=" << cOid
                         << " alignedCols=" << prof.alignedCols.size()
                         << " snps=" << prof.snps.size()
                         << " cov=[" << prof.bbCovBegin << "," << prof.bbCovEnd << ")"
                         << " noisyRangeCount=" << prof.noisyRanges.size()
                         << " noisyBases=" << noisyBases << endl;
                }
                if (seqDebugThis) {
                    cout << "      snpDump window=" << window.windowId << " oid=" << cOid
                         << " nSnps=" << prof.snps.size() << " first20:";
                    for (size_t i = 0; i < prof.snps.size() && i < 20; i++) {
                        cout << " (" << prof.snps[i].bbPos << ",alt="
                             << "ACGTN"[prof.snps[i].altBase < 4 ? prof.snps[i].altBase : 4] << ")";
                    }
                    cout << endl;
                }

                profiles.push_back(std::move(prof));
                membersAligned.fetch_add(1);
            }

            if (profiles.size() < 2) continue;

            const vector<uint8_t> bbWindowBases(
                bbSeqVec.begin() + windowBbBegin, bbSeqVec.begin() + windowBbEnd);
            const auto snarls = findLeafSnarlsFromProfiles(
                bbOid, profiles, bbWindowBases, windowBbBegin,
                /*minSupport=*/6, /*minVaf=*/0.12,
                /*dropHomopolymer=*/true, /*dropRepeat=*/true);
            snarlsFound.fetch_add(snarls.size());

            for (const LeafSnarl& s : snarls) {
                const int64_t localStart = int64_t(s.start) - int64_t(windowBbBegin);
                const int64_t ctxBegin = max<int64_t>(0, localStart - 10);
                const int64_t ctxEnd = min<int64_t>(int64_t(bbWindowBases.size()), localStart + 12);
                string ctx;
                for (int64_t p = ctxBegin; p < ctxEnd; p++) {
                    const uint8_t b = bbWindowBases[size_t(p)];
                    ctx += "ACGTN"[b < 4 ? b : 4];
                }
                cout << "    windowKsw2LeafSnarl window=" << window.windowId
                     << " bb=" << bbOid
                     << " [" << s.start << "," << s.end << ")"
                     << " ctx=" << ctx
                     << " alleles=" << s.alleles.size();
                for (const LeafSnarlAllele& al : s.alleles) {
                    cout << " {n=" << al.members.size() << " bases=";
                    for (uint8_t b : al.bases) cout << "ACGTN"[b < 4 ? b : 4];
                    cout << "}";
                }
                cout << endl;
            }

            windowsProcessed.fetch_add(1);
        }
    };

    vector<std::thread> threads;
    threads.reserve(threadCount);
    for (uint64_t t = 0; t < threadCount; t++) threads.emplace_back(worker);
    for (auto& th : threads) th.join();

    cout << timestamp << "computeWindowKsw2LeafSnarls: processed " << windowsProcessed.load()
         << " windows, aligned " << membersAligned.load() << " members, found "
         << snarlsFound.load() << " leaf snarls." << endl;
}
