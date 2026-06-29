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

// Sliding-window CIGAR-density noise tracker (port of pgphase's XidQueue).
// Each variant event (mismatch=1, insertion=len, deletion=len) is pushed in
// backbone-coordinate order. Within a window of `win` backbone bases, if the
// summed event size exceeds `maxS`, the spanned backbone interval is flagged
// noisy. Contiguous/overlapping noisy spans are merged. The emitted ranges are
// half-open [begin, end) backbone positions.
//
// Differences vs pgphase: no soft-clip handling (ksw2 segments are global, so
// there are no end clips), and coordinates are uint32_t backbone positions.
struct KwNoiseTracker {
    int win;
    int maxS;
    // Pending events in the current window: backbone start pos, ref-length, size.
    vector<uint32_t> pos;
    vector<uint32_t> len;
    vector<int> count;
    size_t front = 0;
    long total = 0;
    // Current open noisy span and the queue indices spanning it.
    long curStart = -1;
    long curEnd = -1;
    vector<pair<uint32_t, uint32_t>>& out; // merged noisy ranges [begin,end)

    KwNoiseTracker(int win_, int maxS_, vector<pair<uint32_t, uint32_t>>& out_)
        : win(win_), maxS(maxS_), out(out_) {}

    // Push one variant event at backbone position p, ref-length l, size c.
    void observe(uint32_t p, uint32_t l, int c) {
        pos.push_back(p);
        len.push_back(l);
        count.push_back(c);
        total += c;
        const size_t rear = pos.size() - 1;

        // Evict events whose ref-end falls before the window's trailing edge.
        while (front <= rear &&
               int64_t(pos[front]) + int64_t(len[front]) - 1 <= int64_t(p) - win) {
            total -= count[front];
            ++front;
        }

        if (c <= 0) return;
        if (total <= maxS) return;

        const long noisyStart = long(pos[front]);
        const long noisyEnd = long(pos[rear]) + long(len[rear]); // half-open

        if (curStart == -1) {
            curStart = noisyStart;
            curEnd = noisyEnd;
            return;
        }
        if (noisyStart <= curEnd) { // overlap/adjacent: extend
            curEnd = max(curEnd, noisyEnd);
            return;
        }
        // Disjoint: flush previous span, open a new one.
        out.push_back({uint32_t(curStart), uint32_t(curEnd)});
        curStart = noisyStart;
        curEnd = noisyEnd;
    }

    // Flush the final open span. Call once after all events.
    void finish() {
        if (curStart == -1) return;
        out.push_back({uint32_t(curStart), uint32_t(curEnd)});
        curStart = -1;
        curEnd = -1;
    }
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
    // Sorted, non-overlapping [begin, end) backbone ranges flagged noisy by the
    // CIGAR-density filter (too many mismatch/indel events locally). A read's
    // SNP votes inside these ranges are not trusted.
    vector<pair<uint32_t, uint32_t>> noisyRanges;

    static bool inRanges(const vector<pair<uint32_t, uint32_t>>& ranges, uint32_t pos) {
        auto it = upper_bound(ranges.begin(), ranges.end(), pos,
            [](uint32_t p, const pair<uint32_t, uint32_t>& r) { return p < r.first; });
        if (it != ranges.begin()) {
            --it;
            if (pos >= it->first && pos < it->second) return true;
        }
        return false;
    }

    // True if backbone position pos falls inside a deletion in this read.
    bool isDeleted(uint32_t pos) const { return inRanges(deletionRanges, pos); }
    // True if backbone position pos falls inside a noisy region in this read.
    bool isNoisy(uint32_t pos) const { return inRanges(noisyRanges, pos); }
};

// One read's allele row across the window's het sites (pgphase ReadVariantProfile).
// alleles[i] is this read's call at hetSnps[i]:
//   0 = ref, 1 = alt, -1 = not covered / non-informative.
struct KwReadVariantProfile {
    OrientedReadId oid;
    int startVarIdx = -1;     // first het column this read informs (-1 = none)
    int endVarIdx = -1;       // last het column this read informs
    vector<int> alleles;      // one entry per het site, column order
};

// Build the row-oriented read x het-site matrix from the column-oriented
// window.hetSnps (altReads / refReads). One row per read that touches at least
// one het site; the backbone read is row 0 (ref at every het site it spans).
// This is the per-window analogue of pgphase collect_read_var_profile and is
// what you inspect to check the het sites for a window.
static vector<KwReadVariantProfile> buildWindowHetSiteMatrix(const AnchorWindow& window) {
    const size_t nSites = window.hetSnps.size();

    // Map oriented-read value -> row index, allocating rows on first sight.
    unordered_map<uint64_t, size_t> rowOf;
    vector<KwReadVariantProfile> rows;

    auto rowFor = [&](OrientedReadId oid) -> KwReadVariantProfile& {
        const uint64_t key = oid.getValue();
        auto it = rowOf.find(key);
        if (it != rowOf.end()) return rows[it->second];
        rowOf.emplace(key, rows.size());
        KwReadVariantProfile p;
        p.oid = oid;
        p.alleles.assign(nSites, -1);
        rows.push_back(move(p));
        return rows.back();
    };

    for (size_t col = 0; col < nSites; col++) {
        const auto& hs = window.hetSnps[col];
        for (const OrientedReadId oid : hs.altReads) rowFor(oid).alleles[col] = 1;
        for (const OrientedReadId oid : hs.refReads) {
            auto& row = rowFor(oid);
            if (row.alleles[col] != 1) row.alleles[col] = 0; // alt wins ties
        }
    }

    // Fill covered span [startVarIdx, endVarIdx] per row.
    for (auto& row : rows) {
        for (size_t i = 0; i < nSites; i++) {
            if (row.alleles[i] < 0) continue;
            if (row.startVarIdx < 0) row.startVarIdx = int(i);
            row.endVarIdx = int(i);
        }
    }

    return rows;
}

// Print the read x het-site matrix for a window. Columns are het sites in
// backbone-position order; '.' = not covered, '0' = ref, '1' = alt.
static void printWindowHetSiteMatrix(const AnchorWindow& window) {
    const size_t nSites = window.hetSnps.size();
    if (nSites == 0) {
        cout << "    hetMatrix bb=" << window.backboneOrientedReadId
             << " window=" << window.windowId << " no het sites" << endl;
        return;
    }

    const vector<KwReadVariantProfile> rows = buildWindowHetSiteMatrix(window);

    cout << "    hetMatrix bb=" << window.backboneOrientedReadId
         << " window=" << window.windowId
         << " sites=" << nSites << " reads=" << rows.size() << endl;

    // Header: het-site backbone positions.
    cout << "      read \\ bbPos:";
    for (const auto& hs : window.hetSnps) cout << ' ' << hs.bbPos;
    cout << endl;

    for (const auto& row : rows) {
        cout << "      " << row.oid << " ";
        for (size_t i = 0; i < nSites; i++) {
            const int a = row.alleles[i];
            cout << (a < 0 ? '.' : char('0' + a));
        }
        cout << "  [" << row.startVarIdx << ".." << row.endVarIdx << "]" << endl;
    }
}

} // anonymous namespace


// Detect clean het SNPs in an anchor window using banded ksw2 2-piece affine
// alignment of each member's inter-anchor segments against the backbone.
uint32_t Assembler::ksw2DetectSnpsInWindow(
    AnchorWindow& window,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const AlignOptions& alignOptions,
    int noisyRegSlideWin,
    int noisyRegMaxXgaps) const
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
        vector<pair<uint32_t, uint32_t>>& outDels,
        KwNoiseTracker& noise)
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
                            // Mismatch: 1 variant base at this backbone position.
                            noise.observe(bbPos, 1, 1);
                        }
                    }
                    qpos += len;
                    tpos += len;
                } else if (op == 1) { // I: insertion in member, no backbone column
                    // Insertion: counts toward local noise density (size = len),
                    // anchored at the current backbone position, ref-length 0.
                    // No SNP/indel variant emitted.
                    noise.observe(tpos, 0, int(len));
                    qpos += len;
                } else if (op == 2) { // D: deletion in member, advances backbone
                    // Deletion: counts toward local noise density (size = len)
                    // over the deleted backbone span.
                    noise.observe(tpos, len, int(len));
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

        // Per-read CIGAR-density noise tracker. Events from all segments are fed
        // in increasing backbone-position order (pins are sorted), so the
        // sliding window over backbone coordinates is well-formed.
        KwNoiseTracker noise(noisyRegSlideWin, noisyRegMaxXgaps, prof.noisyRanges);

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
                prof.snps, prof.deletionRanges, noise);
        }
        noise.finish();

        // noisyRanges are emitted already sorted and merged by the tracker,
        // but adjacent spans flushed separately can still touch; merge them so
        // isNoisy()'s binary search over non-overlapping ranges is correct.
        if (!prof.noisyRanges.empty()) {
            sort(prof.noisyRanges.begin(), prof.noisyRanges.end());
            vector<pair<uint32_t, uint32_t>> merged;
            merged.reserve(prof.noisyRanges.size());
            for (const auto& r : prof.noisyRanges) {
                if (!merged.empty() && r.first <= merged.back().second) {
                    merged.back().second = max(merged.back().second, r.second);
                } else {
                    merged.push_back(r);
                }
            }
            prof.noisyRanges = move(merged);
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
            // A SNP inside this read's locally-noisy region is untrusted
            // (likely an alignment-error cluster, not a real allele): skip it.
            if (prof.isNoisy(s.bbPos)) continue;
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
        // (no base here) or are locally NOISY here (their base is untrusted, so
        // they are counted as neither ref nor alt). The backbone is never
        // deleted/noisy, so it stays counted.
        uint32_t excludedCount = 0;
        for (const auto& prof : profiles) {
            if (pos < prof.bbCovBegin || pos >= prof.bbCovEnd) continue;
            if (prof.isDeleted(pos) || prof.isNoisy(pos))
                excludedCount++;
        }
        const uint32_t effSpanning =
            (spanning > excludedCount) ? spanning - excludedCount : 0;
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

    uint32_t readsWithNoise = 0;
    for (const auto& prof : profiles)
        if (!prof.noisyRanges.empty()) readsWithNoise++;

    cout << "    ksw2DetectSnps bb=" << bbOid
         << " window=[" << windowBbBegin << "," << windowBbEnd << ")"
         << " reads=" << profiles.size()
         << " readsWithNoise=" << readsWithNoise
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
            // A read that deletes this position has no base here, and a read
            // that is locally noisy here has an untrusted base: neither ref nor alt.
            if (prof.isDeleted(passingSnps[i].pos)) continue;
            if (prof.isNoisy(passingSnps[i].pos)) continue;

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

    // Row-oriented read x het-site matrix for inspection (pgphase
    // ReadVariantProfile view). Built from the per-site read lists above.
    printWindowHetSiteMatrix(window);

    return cleanHetSnps;
}
