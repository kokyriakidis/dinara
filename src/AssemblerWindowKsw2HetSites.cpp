/// @file AssemblerWindowKsw2HetSites.cpp
/// @brief Banded per-segment het-SNP detection per anchor window using ksw2.
///
/// Self-contained alternative to cigarDetectSnpsInWindow. Instead of reusing
/// the global all-pairs CIGARs, it aligns each member read against the backbone
/// LOCALLY inside the window:
///
///   - The backbone read defines the coordinate system (one column per backbone
///     base in [windowBbBegin, windowBbEnd), oriented backbone frame).
///   - For each member, its persisted shared-anchor pins give exact
///     correspondence points with the backbone. Between consecutive pins we run
///     ONE banded ksw2 global alignment (2-piece affine gap, minimap2 model) of
///     the member's inter-anchor segment against the backbone's inter-anchor
///     segment. The band follows the pin geometry, so every DP is tiny.
///   - The resulting CIGAR is projected onto backbone columns. Mismatched M
///     columns become SNP candidates. (SNPs only for now: indel ops advance the
///     coordinates but emit no variant.)
///
/// Pins work for transitive members too (reads sharing anchors with the
/// backbone but never directly overlapping it), so the window is fully
/// self-contained: no dependence on the global alignment table for variant
/// calling.
///
/// SNP candidates are aggregated, coverage is computed by sweep line, and the
/// same het filters as cigarDetectSnpsInWindow (alt/ref support, allele
/// frequency, Fisher strand bias, homopolymer/repeat context) gate the calls.
/// Output is written into AnchorWindow::hetSnps / cleanHetSnpCount in the same
/// format, so downstream phasing is unchanged.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "AssemblerOptions.hpp"
#include "PhasingKmeansTypes.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "invalid.hpp"
#include "timestamp.hpp"

extern "C" {
#include "ksw2.h"
}

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace dinara;

namespace {

// Minimum reads (incl. backbone) for a window to be considered.
constexpr uint32_t kwMinReadCoverage = 6;
constexpr uint32_t kwMinSnpAltSupport = 3;
constexpr uint32_t kwMinSnpRefSupport = 3;

// Extra band slack added on top of the pin-implied diagonal offset. Inter-anchor
// segments are short and nearly diagonal; a small constant absorbs indel jitter.
constexpr int kwBandSlack = 30;

// A SNP observed on one member at a backbone position.
struct KwSnp {
    uint32_t bbPos;   // oriented backbone base position
    uint8_t altBase;  // 0=A,1=C,2=G,3=T
};

// Per-member outcome: the SNPs it carries, its backbone coverage range, and the
// backbone positions it DELETES (member has no base there). Deleted positions
// must not be counted as ref support for a SNP, so they are excluded from
// effective coverage.
struct KwMemberProfile {
    OrientedReadId oid;
    vector<KwSnp> snps;
    uint32_t bbCovBegin = 0; // oriented backbone coordinates
    uint32_t bbCovEnd = 0;
    // Sorted, non-overlapping [begin, end) backbone ranges deleted by this read.
    vector<pair<uint32_t, uint32_t>> deletionRanges;

    // True if backbone position pos falls inside a deletion in this read.
    bool isDeleted(uint32_t pos) const {
        auto it = upper_bound(deletionRanges.begin(), deletionRanges.end(), pos,
            [](uint32_t p, const pair<uint32_t, uint32_t>& r) { return p < r.first; });
        if (it != deletionRanges.begin()) {
            --it;
            if (pos >= it->first && pos < it->second) return true;
        }
        return false;
    }
};

} // anonymous namespace


// Detect clean het SNPs in an anchor window using banded ksw2 2-piece affine
// alignment of each member's inter-anchor segments against the backbone.
uint32_t Assembler::ksw2DetectSnpsInWindow(
    AnchorWindow& window,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const AlignOptions& alignOptions) const
{
    const Reads& rds = getReads();
    const auto& mkrs = *markers;
    const uint64_t k = assemblerInfo->k;
    const OrientedReadId bbOid = window.backboneOrientedReadId;
    const ReadId bbReadId = bbOid.getReadId();
    const uint32_t bbLen = uint32_t(rds.getRead(bbReadId).baseCount);
    const auto bbJ = journeys[bbOid];

    if (window.backboneEnd <= window.backboneBegin + 1) return 0;

    // Backbone coordinate range for this window (oriented backbone frame).
    // Matches cigarDetectSnpsInWindow: marker.position (no half-k offset) so the
    // backbone sequence index and the het filters share one coordinate system.
    const uint32_t firstAnchorJP = window.backboneBegin;
    const uint32_t lastAnchorJP = window.backboneEnd - 1;
    if (firstAnchorJP >= bbJ.size() || lastAnchorJP >= bbJ.size()) return 0;

    const uint32_t firstOrd = anchors.getOrdinal(bbJ[firstAnchorJP], bbOid);
    const uint32_t lastOrd = anchors.getOrdinal(bbJ[lastAnchorJP], bbOid);
    if (firstOrd == invalid<uint32_t> || lastOrd == invalid<uint32_t>) return 0;

    const uint32_t windowBbBegin = mkrs[bbOid.getValue()][firstOrd].position;
    const uint32_t windowBbEnd = mkrs[bbOid.getValue()][lastOrd].position + uint32_t(k);
    if (windowBbEnd <= windowBbBegin) return 0;

    // Backbone sequence over the FULL read (0-3 encoded). Indexed by oriented
    // backbone base position, so het filters can inspect flanking context.
    vector<uint8_t> bbSeqVec(bbLen);
    for (uint32_t i = 0; i < bbLen; i++)
        bbSeqVec[i] = rds.getOrientedReadBase(bbOid, i).value;
    const uint8_t* bbSeq = bbSeqVec.data();

    // ksw2 scoring: 2-piece affine, reused from AlignOptions (minimap2 model).
    // ksw2 takes positive penalties; AlignOptions stores match positive and
    // mismatch negative, so negate mismatch into a positive penalty.
    const int8_t kswMatch = int8_t(alignOptions.overlapDpMatchScore);
    const int8_t kswMismatch = int8_t(-alignOptions.overlapDpMismatchScore);
    const int8_t kswGapO1 = int8_t(alignOptions.overlapDpGapOpen1);
    const int8_t kswGapE1 = int8_t(alignOptions.overlapDpGapExtend1);
    const int8_t kswGapO2 = int8_t(alignOptions.overlapDpGapOpen2);
    const int8_t kswGapE2 = int8_t(alignOptions.overlapDpGapExtend2);

    // 5x5 score matrix (A,C,G,T,N). N matches nothing.
    int8_t mat[25];
    for (int a = 0; a < 5; a++)
        for (int b = 0; b < 5; b++)
            mat[a * 5 + b] = (a == b && a < 4) ? kswMatch : kswMismatch;

    // Align a member segment against a backbone segment with banded ksw2 global
    // alignment, projecting mismatched columns onto backbone SNP positions.
    // bbSegBegin/cSegBegin are oriented base positions; lengths are in bases.
    auto alignSegment = [&, windowBbBegin, windowBbEnd](
        OrientedReadId cOid,
        uint32_t bbSegBegin, uint32_t bbSegLen,
        uint32_t cSegBegin, uint32_t cSegLen,
        vector<KwSnp>& outSnps,
        vector<pair<uint32_t, uint32_t>>& outDels)
    {
        if (bbSegLen == 0 || cSegLen == 0) return;

        // Build 0-3 encoded query (member) and target (backbone) for this gap.
        static thread_local vector<uint8_t> query;
        static thread_local vector<uint8_t> target;
        query.resize(cSegLen);
        target.resize(bbSegLen);
        for (uint32_t i = 0; i < cSegLen; i++)
            query[i] = rds.getOrientedReadBase(cOid, cSegBegin + i).value;
        for (uint32_t i = 0; i < bbSegLen; i++)
            target[i] = rds.getOrientedReadBase(bbOid, bbSegBegin + i).value;

        // Band from segment length difference plus a small constant slack.
        const int lenDiff = abs(int(bbSegLen) - int(cSegLen));
        const int band = lenDiff + kwBandSlack;

        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ez));
        // Global alignment (no extension/Z-drop), 2-piece affine.
        ksw_extd2_sse(
            nullptr,
            int(cSegLen), query.data(),
            int(bbSegLen), target.data(),
            5, mat,
            kswGapO1, kswGapE1, kswGapO2, kswGapE2,
            band, -1, 0, 0, &ez);

        if (ez.n_cigar > 0 && ez.cigar != nullptr) {
            // Walk the CIGAR. ksw2 op codes: 0=M, 1=I (query/member), 2=D
            // (target/backbone). Query=member, target=backbone.
            uint32_t qpos = cSegBegin;
            uint32_t tpos = bbSegBegin;
            for (int ci = 0; ci < ez.n_cigar; ci++) {
                const uint32_t len = ez.cigar[ci] >> 4;
                const uint32_t op = ez.cigar[ci] & 0xf;
                if (op == 0) { // M: aligned, may be match or mismatch
                    for (uint32_t j = 0; j < len; j++) {
                        const uint32_t bbPos = tpos + j;
                        if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;
                        const uint8_t cb = rds.getOrientedReadBase(cOid, qpos + j).value;
                        const uint8_t bb = bbSeqVec[bbPos];
                        if (cb != bb && cb < 4 && bb < 4) {
                            outSnps.push_back(KwSnp{bbPos, cb});
                        }
                    }
                    qpos += len;
                    tpos += len;
                } else if (op == 1) { // I: insertion in member, no backbone column
                    qpos += len;
                } else if (op == 2) { // D: deletion in member, advances backbone
                    // Record the deleted backbone span (clamped to the window)
                    // so these positions are not miscounted as ref support.
                    const uint32_t db = max(tpos, windowBbBegin);
                    const uint32_t de = min(tpos + len, windowBbEnd);
                    if (db < de) outDels.push_back({db, de});
                    tpos += len;
                } else {
                    // N/other: advance backbone defensively.
                    tpos += len;
                }
            }
        }
        if (ez.cigar) free(ez.cigar);
    };

    // Build per-member profiles by walking persisted shared-anchor pins.
    vector<KwMemberProfile> profiles;
    profiles.reserve(window.readIntervals.size());

    for (size_t ri = 1; ri < window.readIntervals.size(); ri++) {
        const AnchorWindowReadInterval& interval = window.readIntervals[ri];
        const OrientedReadId cOid = interval.orientedReadId;
        if (cOid == bbOid) continue;
        if (cOid.getValue() >= journeys.size()) continue;
        if (interval.sharedPins.size() < 2) continue;

        // Convert persisted pins (journey positions) to base positions on both
        // backbone and member. Pins arrive sorted by backbone journey position.
        struct Pin { uint32_t bbPos; uint32_t cPos; };
        vector<Pin> pins;
        pins.reserve(interval.sharedPins.size());
        for (const AnchorWindowSharedPin& sp : interval.sharedPins) {
            if (sp.backboneJourneyPos >= bbJ.size()) continue;
            const Shasta2AnchorId aid = bbJ[sp.backboneJourneyPos];
            const uint32_t bbOrd = anchors.getOrdinal(aid, bbOid);
            const uint32_t cOrd = anchors.getOrdinal(aid, cOid);
            if (bbOrd == invalid<uint32_t> || cOrd == invalid<uint32_t>) continue;
            if (bbOrd >= mkrs[bbOid.getValue()].size()) continue;
            if (cOrd >= mkrs[cOid.getValue()].size()) continue;
            const uint32_t bbPos = mkrs[bbOid.getValue()][bbOrd].position;
            const uint32_t cPos = mkrs[cOid.getValue()][cOrd].position;
            if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;
            pins.push_back(Pin{bbPos, cPos});
        }
        if (pins.size() < 2) continue;

        // Pins should already be sorted by backbone position; ensure it.
        sort(pins.begin(), pins.end(),
            [](const Pin& a, const Pin& b) { return a.bbPos < b.bbPos; });

        KwMemberProfile prof;
        prof.oid = cOid;
        prof.bbCovBegin = max(pins.front().bbPos, windowBbBegin);
        prof.bbCovEnd = min(pins.back().bbPos + uint32_t(k), windowBbEnd);

        // Align the gap between each consecutive pin pair.
        for (size_t pi = 0; pi + 1 < pins.size(); pi++) {
            const Pin& left = pins[pi];
            const Pin& right = pins[pi + 1];

            // Inter-anchor segment: from the END of the left anchor (+k) to the
            // START of the right anchor, on both backbone and member.
            const uint32_t bbSegBegin = left.bbPos + uint32_t(k);
            const uint32_t bbSegEnd = right.bbPos;
            const uint32_t cSegBegin = left.cPos + uint32_t(k);
            const uint32_t cSegEnd = right.cPos;
            if (bbSegEnd <= bbSegBegin || cSegEnd <= cSegBegin) continue;
            if (bbSegEnd <= windowBbBegin || bbSegBegin >= windowBbEnd) continue;

            alignSegment(cOid,
                bbSegBegin, bbSegEnd - bbSegBegin,
                cSegBegin, cSegEnd - cSegBegin,
                prof.snps, prof.deletionRanges);
        }

        // Sort and merge deletion ranges so isDeleted() (binary search over
        // non-overlapping [begin,end)) is correct.
        if (!prof.deletionRanges.empty()) {
            sort(prof.deletionRanges.begin(), prof.deletionRanges.end());
            vector<pair<uint32_t, uint32_t>> merged;
            merged.reserve(prof.deletionRanges.size());
            for (const auto& r : prof.deletionRanges) {
                if (!merged.empty() && r.first <= merged.back().second) {
                    merged.back().second = max(merged.back().second, r.second);
                } else {
                    merged.push_back(r);
                }
            }
            prof.deletionRanges = move(merged);
        }

        profiles.push_back(move(prof));
    }

    if (profiles.size() + 1 < kwMinReadCoverage) return 0; // +1 for backbone

    // Aggregate SNP candidates: (bbPos, altBase) -> {fwd, rev, total}.
    struct SnpAccum { uint32_t fwd = 0; uint32_t rev = 0; uint32_t total = 0; };
    unordered_map<uint64_t, SnpAccum> snpCounts;
    auto snpKey = [](uint32_t pos, uint8_t alt) -> uint64_t {
        return (uint64_t(pos) << 8) | uint64_t(alt);
    };
    unordered_set<uint32_t> snpPositions;

    for (const auto& prof : profiles) {
        const bool isFwd = (prof.oid.getStrand() == 0);
        for (const KwSnp& s : prof.snps) {
            auto& acc = snpCounts[snpKey(s.bbPos, s.altBase)];
            if (isFwd) acc.fwd++; else acc.rev++;
            acc.total++;
            snpPositions.insert(s.bbPos);
        }
    }

    // Spanning coverage at each SNP position via sweep line over member ranges.
    struct CovEvent { uint32_t pos; int delta; };
    vector<CovEvent> covEvents;
    covEvents.reserve(profiles.size() * 2);
    for (const auto& prof : profiles) {
        if (prof.bbCovBegin < prof.bbCovEnd) {
            covEvents.push_back({prof.bbCovBegin, +1});
            covEvents.push_back({prof.bbCovEnd, -1});
        }
    }
    sort(covEvents.begin(), covEvents.end(),
        [](const CovEvent& a, const CovEvent& b) {
            return a.pos < b.pos || (a.pos == b.pos && a.delta > b.delta);
        });

    vector<uint32_t> sortedSnpPositions(snpPositions.begin(), snpPositions.end());
    sort(sortedSnpPositions.begin(), sortedSnpPositions.end());

    unordered_map<uint32_t, uint32_t> spanningCount;
    spanningCount.reserve(sortedSnpPositions.size());
    {
        int running = 0;
        size_t ei = 0;
        for (uint32_t pos : sortedSnpPositions) {
            while (ei < covEvents.size() && covEvents[ei].pos <= pos) {
                running += covEvents[ei].delta;
                ei++;
            }
            // Backbone itself spans the whole window: +1.
            spanningCount[pos] = uint32_t(max(running, 0)) + 1;
        }
    }

    // Classify SNPs with the same filters as cigarDetectSnpsInWindow.
    KmPhasingOptions opts;
    struct PassingSnp {
        uint32_t pos; uint8_t altBase;
        uint32_t altCov; uint32_t refCov; uint32_t spanning;
    };
    vector<PassingSnp> passingSnps;
    uint32_t cleanHetSnps = 0;

    uint32_t passAltSupport = 0, failRefSupport = 0, failAf = 0;
    uint32_t failStrandBias = 0, failHomopolymer = 0;

    for (uint32_t pos : snpPositions) {
        const uint32_t spanning = spanningCount[pos];
        if (spanning == 0) continue;

        // Effective spanning: exclude members that DELETE this backbone position
        // (they have no base here, so they are neither ref nor alt). The
        // backbone is never deleted, so it stays counted.
        uint32_t delCount = 0;
        for (const auto& prof : profiles) {
            if (pos >= prof.bbCovBegin && pos < prof.bbCovEnd && prof.isDeleted(pos))
                delCount++;
        }
        const uint32_t effSpanning = (spanning > delCount) ? spanning - delCount : 0;
        if (effSpanning == 0) continue;

        for (uint8_t alt = 0; alt < 4; alt++) {
            auto it = snpCounts.find(snpKey(pos, alt));
            if (it == snpCounts.end()) continue;
            const auto& acc = it->second;

            if (acc.total < kwMinSnpAltSupport) continue;
            passAltSupport++;

            const uint32_t refCov = (effSpanning > acc.total) ? effSpanning - acc.total : 0;
            if (refCov < kwMinSnpRefSupport) { failRefSupport++; continue; }

            const double af = double(acc.total) / double(effSpanning);
            if (af < opts.minAf || af > opts.maxAf) { failAf++; continue; }

            const int expected = int(acc.total) / 2;
            if (expected > 0) {
                const double p = kmFisherExactTwoTail(
                    int(acc.fwd), int(acc.rev), expected, expected);
                if (p < opts.strandBiasPval) { failStrandBias++; continue; }
            }

            KmVarKey vkey;
            vkey.pos = pos;
            vkey.type = KmVarType::Snp;
            vkey.altBase = alt;
            vkey.refLen = 1;
            vkey.altLen = 1;
            if (kmIsHomopolymer(bbSeq, bbLen, vkey, 0) ||
                kmIsRepeatRegion(bbSeq, bbLen, vkey, 0)) {
                failHomopolymer++;
                continue;
            }

            passingSnps.push_back({pos, alt, acc.total, refCov, effSpanning});
            cleanHetSnps++;
        }
    }

    cout << "    ksw2DetectSnps bb=" << bbOid
         << " window=[" << windowBbBegin << "," << windowBbEnd << ")"
         << " reads=" << profiles.size()
         << " snpPositions=" << snpPositions.size()
         << " passAltSupport=" << passAltSupport
         << " failRefSupport=" << failRefSupport
         << " failAF=" << failAf
         << " failStrandBias=" << failStrandBias
         << " failHomopolymer=" << failHomopolymer
         << " cleanHetSnps=" << cleanHetSnps << endl;

    // Store passing het SNPs with per-read alt/ref lists (same format as the
    // CIGAR path so downstream phasing is unchanged).
    window.hetSnps.resize(passingSnps.size());
    for (size_t i = 0; i < passingSnps.size(); i++) {
        auto& hs = window.hetSnps[i];
        hs.bbPos = passingSnps[i].pos;
        hs.refBase = bbSeqVec[passingSnps[i].pos];
        hs.altBase = passingSnps[i].altBase;
        hs.altCov = passingSnps[i].altCov;
        hs.refCov = passingSnps[i].refCov;
        hs.spanning = passingSnps[i].spanning;

        // Backbone is ref at every position.
        if (passingSnps[i].pos >= windowBbBegin && passingSnps[i].pos < windowBbEnd)
            hs.refReads.push_back(bbOid);

        for (const auto& prof : profiles) {
            if (passingSnps[i].pos < prof.bbCovBegin ||
                passingSnps[i].pos >= prof.bbCovEnd) continue;
            // A read that deletes this position has no base here: neither ref nor alt.
            if (prof.isDeleted(passingSnps[i].pos)) continue;

            bool hasThisAlt = false;
            bool hasOtherAlt = false;
            for (const KwSnp& s : prof.snps) {
                if (s.bbPos == passingSnps[i].pos) {
                    if (s.altBase == passingSnps[i].altBase) hasThisAlt = true;
                    else hasOtherAlt = true;
                    break;
                }
            }
            if (hasThisAlt) hs.altReads.push_back(prof.oid);
            else if (!hasOtherAlt) hs.refReads.push_back(prof.oid);
        }
    }

    window.cleanHetSnpCount = cleanHetSnps;
    return cleanHetSnps;
}
