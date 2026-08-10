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
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace dinara;

namespace {

// Ad-hoc perf counter: total ksw_extd2_sse calls (one per member per
// inter-anchor segment), to compare call-count against the interval-abPOA
// engine's msaRuns (one shared abpoa_msa call per interval, covering ALL
// members touching it). Printed via dinaraPrintKswCallCount() when
// DINARA_HET_TIMING is set (reusing that engine's timing flag for an
// apples-to-apples side-by-side).
inline std::atomic<std::uint64_t> kswAlignCalls{0};

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

// One aligned (CIGAR M) column of a member against the backbone. Captured only
// by the het-bubble path so it can, after SNP calling, reconstruct the k=2
// marker rawPositions (read base position at the bubble's flank columns) and
// verify flank-linearity directly from the member pileup. bbPos is the oriented
// backbone base position; readBase is the member's base there (0-3); readPos is
// the member's ABSOLUTE oriented base position, which is exactly the rawPosition
// abPOA would store for a k=2 anchor member at this column.
struct KwAlignedCol {
    uint32_t bbPos;
    uint32_t readPos;
    uint8_t readBase;
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
    // Read (query) base span [readCBegin, readCEnd) corresponding to the pin
    // range, in the member's own oriented coordinates. Used by the two-pass
    // consensus path to re-extract the read for realignment against the
    // consensus. readCBegin = front pin cPos, readCEnd = back pin cPos + k.
    uint32_t readCBegin = 0;
    uint32_t readCEnd = 0;
    // Pin (shared-anchor) coordinates for the two-pass path: parallel arrays of
    // backbone base position and read base position, sorted by backbone position.
    // Lets pass 2 realign the same short inter-anchor segments against the
    // consensus instead of the raw backbone.
    vector<uint32_t> pinBb;
    vector<uint32_t> pinC;
    // Every aligned (M) column of this member, in increasing backbone-position
    // order. Populated ONLY by the het-bubble path (captureAlignedCols=true) so
    // it can look up the member's read base + absolute read position at any
    // flank backbone column when synthesizing k=2 anchor members. Empty for the
    // flat-SNP path. Sorted by bbPos, one entry per aligned backbone position.
    vector<KwAlignedCol> alignedCols;

    // Read base at backbone position pos (0-3), or 0xff if this member does not
    // align a base there (deleted / outside coverage). Binary search over the
    // sorted alignedCols.
    const KwAlignedCol* colAt(uint32_t pos) const {
        auto it = lower_bound(alignedCols.begin(), alignedCols.end(), pos,
            [](const KwAlignedCol& c, uint32_t p) { return c.bbPos < p; });
        if (it != alignedCols.end() && it->bbPos == pos) return &*it;
        return nullptr;
    }
    // Insertions relative to the backbone: read bases that align to no backbone
    // column, keyed by the backbone position they are inserted BEFORE (the tpos
    // of the ksw2 'I' op). Captured only by the two-pass consensus path so it can
    // build insertion columns. Each entry: (bbPosBefore, inserted read bases).
    // Multiple entries at the same bbPosBefore are possible across segments.
    vector<pair<uint32_t, vector<uint8_t>>> insertions;

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

// Prints the total ksw_extd2_sse call count (one per member per inter-anchor
// segment) accumulated so far, when DINARA_HET_TIMING is set -- a side-by-side
// comparison point against the interval-abPOA engine's own [HetTiming]
// msaRuns count (one shared abpoa_msa call per interval, covering ALL members
// touching it, versus one independent ksw2 call per member per segment here).
// Lives in dinara::main (matching where main.cpp's assemble() calls it,
// since a local extern declaration there resolves within that same nested
// namespace, not global scope) even though it only needs the anonymous-
// namespace counter above it, visible here via normal enclosing-scope lookup.
namespace dinara { namespace main {
void dinaraPrintKswCallCount() {
    if (std::getenv("DINARA_HET_TIMING") == nullptr) return;
    std::cout << "[HetTiming] ksw2 alignSegment (ksw_extd2_sse) calls: "
              << kswAlignCalls.load() << std::endl;
}
} }


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
    // The gap-open/extend parameters below genuinely are positive penalty
    // magnitudes that ksw2 subtracts internally, but the substitution matrix
    // (mat[], used for kswMatch/kswMismatch) is a SIGNED addend matrix: ksw2
    // adds mat[a][b] directly into the DP score, so a mismatch entry must
    // stay negative, matching how AlignmentInfo::dpScore already sums these
    // same AlignOptions fields elsewhere. overlapDpMismatchScore is already
    // stored negative (default -4, symmetric with the positive +2 match
    // score) -- do NOT negate it again (confirmed empirically: an extra
    // negation here turns a -4 penalty into a +4 reward for mismatches,
    // producing a ~19% aggregate mismatch rate on real HiFi data instead of
    // the expected ~0.5%).
    const int8_t kswMatch = int8_t(alignOptions.overlapDpMatchScore);
    const int8_t kswMismatch = int8_t(alignOptions.overlapDpMismatchScore);
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


// Detect clean het bubbles in an anchor window using the SAME banded ksw2
// member pileup as ksw2DetectSnpsInWindow, but emit AnchorWindow::hetBubbles in
// the exact format the abPOA path (testAbpoaMultiSegmentMSA) produces, so it is
// a drop-in alternative that the downstream plan/append/stage passes consume
// unchanged.
//
// This is deliberately a PARALLEL function (not a refactor of the SNP path): the
// flat-SNP path stays byte-for-byte identical for its own separate pipeline.
// The two share only the file-local pileup helpers (KwSnp, KwNoiseTracker,
// KwMemberProfile).
//
// Frame: hetBubble.backboneOffset = bbPos - windowBbBegin, where windowBbBegin
// == abPOA's backboneBeginPos == marker[firstAnchor].position and the span runs
// to marker[lastAnchor].position + k. Both paths therefore index the same
// [firstAnchorStart, lastAnchorEnd) backbone array with the same origin; there
// is NO half-k shift on the SNP offset (kHalf only moves anchor pins, not SNP
// columns).
//
// rawPosition semantics match abPOA exactly: it is the member's ABSOLUTE
// oriented base position at the k=2 anchor's node column. Arms pin at commonPred
// (bbPos-1); the leading hom pins at predPrev (bbPos-2); the trailing hom pins
// at commonSucc (bbPos+1). The backbone read is an identity-frame member of the
// ref arm and both homs (its rawPosition at column c is simply c).
//
// Flank-linearity (faithful reconstruction of abPOA's degree-1 flank guarantee):
// a SNP at column p is accepted only if columns p-2, p-1, p+1, p+2 all lie
// inside the window AND are homozygous over the SNP's spanning members -- no
// spanning member records a mismatch there and none deletes there. This is the
// pileup equivalent of abPOA requiring predPrev->commonPred and
// commonSucc->succNext to be single-in/single-out linear edges. Without it the
// bracketing k=2 homs would not be shared by all entering/leaving reads.
uint32_t Assembler::ksw2DetectHetBubblesInWindow(
    AnchorWindow& window,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const AlignOptions& alignOptions,
    double hetMinVaf,
    uint64_t hetMinSupport,
    bool hetDropHomopolymer,
    bool hetDropRepeat,
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

    window.hetBubbles.clear();
    if (window.backboneEnd <= window.backboneBegin + 1) return 0;

    // Backbone coordinate range for this window. Identical to the abPOA path:
    // START of the first anchor to END of the last anchor (marker.position, no
    // half-k offset), so backboneOffset = bbPos - windowBbBegin matches the
    // abPOA backboneCodes index exactly.
    const uint32_t firstAnchorJP = window.backboneBegin;
    const uint32_t lastAnchorJP = window.backboneEnd - 1;
    if (firstAnchorJP >= bbJ.size() || lastAnchorJP >= bbJ.size()) return 0;

    const uint32_t firstOrd = anchors.getOrdinal(bbJ[firstAnchorJP], bbOid);
    const uint32_t lastOrd = anchors.getOrdinal(bbJ[lastAnchorJP], bbOid);
    if (firstOrd == invalid<uint32_t> || lastOrd == invalid<uint32_t>) return 0;

    const uint32_t windowBbBegin = mkrs[bbOid.getValue()][firstOrd].position;
    const uint32_t windowBbEnd = mkrs[bbOid.getValue()][lastOrd].position + uint32_t(k);
    if (windowBbEnd <= windowBbBegin) return 0;

    vector<uint8_t> bbSeqVec(bbLen);
    for (uint32_t i = 0; i < bbLen; i++)
        bbSeqVec[i] = rds.getOrientedReadBase(bbOid, i).value;
    const uint8_t* bbSeq = bbSeqVec.data();

    // ksw2 scoring: 2-piece affine, reused from AlignOptions (minimap2 model).
    // mat[] is a signed addend matrix (ksw2 adds it directly into the DP
    // score), so the mismatch entry must stay negative -- overlapDpMismatchScore
    // is already stored that way (default -4, symmetric with the positive +2
    // match score); do not negate it again (see the detailed note on the
    // other kswMismatch above in ksw2DetectSnpsInWindow).
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

    // Same banded ksw2 segment alignment as the SNP path, but ALSO records every
    // aligned (M) column into outCols so the bubble pass can recover member read
    // positions and check flank-linearity from the pileup.
    auto alignSegment = [&, windowBbBegin, windowBbEnd](
        OrientedReadId cOid,
        uint32_t bbSegBegin, uint32_t bbSegLen,
        uint32_t cSegBegin, uint32_t cSegLen,
        vector<KwSnp>& outSnps,
        vector<pair<uint32_t, uint32_t>>& outDels,
        vector<KwAlignedCol>& outCols,
        KwNoiseTracker& noise,
        vector<pair<uint32_t, vector<uint8_t>>>* outIns = nullptr)
    {
        if (bbSegLen == 0 || cSegLen == 0) return;
        kswAlignCalls.fetch_add(1, std::memory_order_relaxed);

        static thread_local vector<uint8_t> query;
        static thread_local vector<uint8_t> target;
        query.resize(cSegLen);
        target.resize(bbSegLen);
        for (uint32_t i = 0; i < cSegLen; i++)
            query[i] = rds.getOrientedReadBase(cOid, cSegBegin + i).value;
        for (uint32_t i = 0; i < bbSegLen; i++)
            target[i] = rds.getOrientedReadBase(bbOid, bbSegBegin + i).value;

        const int lenDiff = abs(int(bbSegLen) - int(cSegLen));
        const int band = lenDiff + kwBandSlack;

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
                        const uint8_t bb = bbSeqVec[bbPos];
                        // Record the aligned column (read base + absolute read
                        // position) for later rawPosition / flank lookup.
                        outCols.push_back(KwAlignedCol{bbPos, rp, cb});
                        if (cb != bb && cb < 4 && bb < 4) {
                            outSnps.push_back(KwSnp{bbPos, cb});
                            noise.observe(bbPos, 1, 1);
                        }
                    }
                    qpos += len;
                    tpos += len;
                } else if (op == 1) { // I
                    noise.observe(tpos, 0, int(len));
                    // Capture inserted read bases (two-pass consensus only).
                    if (outIns != nullptr && tpos >= windowBbBegin && tpos <= windowBbEnd) {
                        vector<uint8_t> ins(len);
                        for (uint32_t j = 0; j < len; j++)
                            ins[j] = rds.getOrientedReadBase(cOid, qpos + j).value;
                        outIns->push_back({tpos, std::move(ins)});
                    }
                    qpos += len;
                } else if (op == 2) { // D
                    noise.observe(tpos, len, int(len));
                    const uint32_t db = max(tpos, windowBbBegin);
                    const uint32_t de = min(tpos + len, windowBbEnd);
                    if (db < de) outDels.push_back({db, de});
                    tpos += len;
                } else {
                    tpos += len;
                }
            }
        }
        if (ez.cigar) free(ez.cigar);
    };

    // Build per-read pins LIVE from the window's backbone anchors, exactly like
    // the abPOA path -- NOT from the persisted window.readIntervals/sharedPins
    // subset. Every anchor a read shares with the backbone inside this window
    // becomes a pin, so per-column read coverage matches the abPOA fold. (The
    // persisted sharedPins were thinned by unclaimed-only capture + LIS +
    // duplicate-kmer drop, which starved the pileup and cost recall.)
    //
    // For each backbone anchor in [backboneBegin, backboneEnd) we read its live
    // member list; each non-backbone member contributes one pin (backbone base
    // position, read base position) via the member's own marker ordinal.
    struct Pin { uint32_t bbPos; uint32_t cPos; };
    unordered_map<uint64_t, vector<Pin>> readPins;

    for (uint32_t jp = window.backboneBegin; jp < window.backboneEnd; jp++) {
        if (jp >= bbJ.size()) break;
        const Shasta2AnchorId aid = bbJ[jp];
        const uint32_t bbOrd = anchors.getOrdinal(aid, bbOid);
        if (bbOrd == invalid<uint32_t> || bbOrd >= mkrs[bbOid.getValue()].size())
            continue;
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
            readPins[cOid.getValue()].push_back(Pin{bbPos, cPos});
        }
    }

    // Keep the longest colinear (strictly diagonal-increasing) run of pins per
    // read: sort by backbone position, then take the LIS over read positions.
    // This drops repeat-induced off-diagonal pins that would otherwise create
    // bogus inter-anchor segments (same purpose the persisted sharedPins LIS
    // served, but computed over the FULL live anchor set).
    auto lisByCPos = [](vector<Pin>& pins) {
        const uint32_t n = uint32_t(pins.size());
        if (n < 2) return;
        vector<uint32_t> tails;              // indices into pins
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
        vector<Pin> kept;
        for (int32_t c = int32_t(tails.back()); c >= 0; c = pred[uint32_t(c)])
            kept.push_back(pins[uint32_t(c)]);
        reverse(kept.begin(), kept.end());
        pins.swap(kept);
    };

    // Build per-member profiles from the live pins, capturing aligned columns.
    vector<KwMemberProfile> profiles;
    profiles.reserve(readPins.size());

    for (auto& [orientedValue, pinsRaw] : readPins) {
        // Map key is the full OrientedReadId::getValue() (read<<1 | strand), so
        // reconstruct with fromValue (NOT the ReadId ctor, which is deleted).
        const OrientedReadId cOid = OrientedReadId::fromValue(ReadId(orientedValue));
        if (pinsRaw.size() < 2) continue;

        // Sort by backbone position, then keep the longest colinear run.
        sort(pinsRaw.begin(), pinsRaw.end(),
            [](const Pin& a, const Pin& b) {
                return a.bbPos < b.bbPos || (a.bbPos == b.bbPos && a.cPos < b.cPos);
            });
        lisByCPos(pinsRaw);
        vector<Pin>& pins = pinsRaw;
        if (pins.size() < 2) continue;

        KwMemberProfile prof;
        prof.oid = cOid;
        prof.bbCovBegin = max(pins.front().bbPos, windowBbBegin);
        prof.bbCovEnd = min(pins.back().bbPos + uint32_t(k), windowBbEnd);
        prof.readCBegin = pins.front().cPos;
        prof.readCEnd = pins.back().cPos + uint32_t(k);
        if (getenv("DINARA_KSW2_TWOPASS") != nullptr) {
            prof.pinBb.reserve(pins.size());
            prof.pinC.reserve(pins.size());
            for (const Pin& pn : pins) { prof.pinBb.push_back(pn.bbPos); prof.pinC.push_back(pn.cPos); }
        }

        KwNoiseTracker noise(noisyRegSlideWin, noisyRegMaxXgaps, prof.noisyRanges);

        // Segment tiling. Default (DINARA_KSW2_FULLSPAN unset) keeps the legacy
        // gap-only tiling that aligns [left+k, right) and treats anchor bodies as
        // implicit ref. That silently defaults every read to ref at any column
        // inside a shared-anchor body, so a het SNP that falls in the anchor body
        // for one haplotype is never seen -> alt support collapses and the site
        // is missed.
        //
        // FULLSPAN mode tiles the WHOLE span with no gaps: segment i spans
        // [pins[i].bbPos, pins[i+1].bbPos) on the backbone (anchor body included
        // at the segment head), plus a final [lastPin.bbPos, lastPin.bbPos+k)
        // tail so the last anchor body is covered too. Every read then has a real
        // aligned base at every column it spans, and memberCall reads the actual
        // base instead of assuming ref.
        // Two-pass consensus (DINARA_KSW2_TWOPASS) forces FULLSPAN and captures
        // insertions so a length-changing consensus can be built.
        const bool twoPass = (getenv("DINARA_KSW2_TWOPASS") != nullptr);
        const bool fullSpan = twoPass || (getenv("DINARA_KSW2_FULLSPAN") != nullptr);
        auto* insSink = twoPass ? &prof.insertions : nullptr;

        for (size_t pi = 0; pi + 1 < pins.size(); pi++) {
            const Pin& left = pins[pi];
            const Pin& right = pins[pi + 1];
            const uint32_t bbSegBegin = fullSpan ? left.bbPos : left.bbPos + uint32_t(k);
            const uint32_t bbSegEnd = right.bbPos;
            const uint32_t cSegBegin = fullSpan ? left.cPos : left.cPos + uint32_t(k);
            const uint32_t cSegEnd = right.cPos;
            if (bbSegEnd <= bbSegBegin || cSegEnd <= cSegBegin) continue;
            if (bbSegEnd <= windowBbBegin || bbSegBegin >= windowBbEnd) continue;

            alignSegment(cOid,
                bbSegBegin, bbSegEnd - bbSegBegin,
                cSegBegin, cSegEnd - cSegBegin,
                prof.snps, prof.deletionRanges, prof.alignedCols, noise, insSink);
        }
        // Cover the trailing anchor body [lastPin, lastPin+k) too, so the final
        // k columns are examined rather than assumed ref.
        if (fullSpan && pins.size() >= 1) {
            const Pin& last = pins.back();
            const uint32_t bbSegBegin = last.bbPos;
            const uint32_t bbSegEnd = min(last.bbPos + uint32_t(k), windowBbEnd);
            const uint32_t cSegBegin = last.cPos;
            const uint32_t cSegEnd = last.cPos + uint32_t(k);
            if (bbSegEnd > bbSegBegin && cSegEnd > cSegBegin &&
                bbSegBegin < windowBbEnd)
                alignSegment(cOid,
                    bbSegBegin, bbSegEnd - bbSegBegin,
                    cSegBegin, cSegEnd - cSegBegin,
                    prof.snps, prof.deletionRanges, prof.alignedCols, noise, insSink);
        }
        noise.finish();

        // FULLSPAN can emit the same backbone column from two adjacent segments
        // only at exact segment joints; dedup alignedCols by bbPos (keep first).
        if (fullSpan && prof.alignedCols.size() > 1) {
            sort(prof.alignedCols.begin(), prof.alignedCols.end(),
                [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos < b.bbPos; });
            prof.alignedCols.erase(
                unique(prof.alignedCols.begin(), prof.alignedCols.end(),
                    [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos == b.bbPos; }),
                prof.alignedCols.end());
        }

        // Aligned columns are produced in segment order; sort by backbone
        // position so colAt() binary search is valid. (Segments are already in
        // order but a defensive sort keeps the invariant explicit.)
        sort(prof.alignedCols.begin(), prof.alignedCols.end(),
            [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos < b.bbPos; });

        if (!prof.noisyRanges.empty()) {
            sort(prof.noisyRanges.begin(), prof.noisyRanges.end());
            vector<pair<uint32_t, uint32_t>> merged;
            merged.reserve(prof.noisyRanges.size());
            for (const auto& r : prof.noisyRanges) {
                if (!merged.empty() && r.first <= merged.back().second)
                    merged.back().second = max(merged.back().second, r.second);
                else merged.push_back(r);
            }
            prof.noisyRanges = move(merged);
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
            prof.deletionRanges = move(merged);
        }

        profiles.push_back(move(prof));
    }

    if (profiles.size() + 1 < kwMinReadCoverage) return 0; // +1 for backbone

    // ------------------------------------------------------------------
    // Two-pass consensus realignment (DINARA_KSW2_TWOPASS).
    //
    // Pass 1 above aligned every read to the BACKBONE independently. Because each
    // read solves its own optimal gap placement, indels shift the SNP column
    // read-to-read: the alt allele smears across neighbouring columns and never
    // reaches support. abPOA avoids this by aligning reads to EACH OTHER.
    //
    // This pass approximates that: build a length-changing consensus from the
    // pass-1 pileup (majority base per backbone column + majority insertion
    // columns), then realign each read to the CONSENSUS. Since the consensus
    // already carries the read population's common indel structure, per-read gap
    // ambiguity collapses and columns agree. SNPs are still called vs the
    // backbone base (ref = backbone), so ref/alt semantics are unchanged; only
    // the column frame is stabilised.
    // ------------------------------------------------------------------
    if (getenv("DINARA_KSW2_TWOPASS") != nullptr) {
        const uint32_t winLen = windowBbEnd - windowBbBegin;

        // Per-backbone-column majority base + coverage, and per-insertion-slot
        // (inserted BEFORE backbone position p) the number of reads inserting and
        // the most common inserted sequence.
        vector<array<uint32_t,4>> baseVotes(winLen, {0,0,0,0});
        vector<uint32_t> colCov(winLen, 0);
        // insertion slot p in [0,winLen]: map inserted-seq -> count
        vector<unordered_map<string,uint32_t>> insVotes(winLen + 1);
        vector<uint32_t> insSlotCov(winLen + 1, 0); // reads spanning the slot

        for (const auto& prof : profiles) {
            for (const KwAlignedCol& c : prof.alignedCols) {
                if (c.bbPos < windowBbBegin || c.bbPos >= windowBbEnd) continue;
                if (c.readBase < 4) {
                    baseVotes[c.bbPos - windowBbBegin][c.readBase]++;
                    colCov[c.bbPos - windowBbBegin]++;
                }
            }
            for (const auto& ins : prof.insertions) {
                if (ins.first < windowBbBegin || ins.first > windowBbEnd) continue;
                const uint32_t slot = ins.first - windowBbBegin;
                string s(ins.second.begin(), ins.second.end());
                insVotes[slot][s]++;
            }
        }
        // Reads spanning each insertion slot (slot p sits between col p-1 and p).
        for (const auto& prof : profiles) {
            uint32_t b = max(prof.bbCovBegin, windowBbBegin);
            uint32_t e = min(prof.bbCovEnd, windowBbEnd);
            for (uint32_t p = b; p <= e; p++) insSlotCov[p - windowBbBegin]++;
        }

        // Build consensus sequence + consToBb (backbone pos per consensus col, or
        // -1 for an inserted column) + bbToCons (M-column index per backbone pos).
        vector<uint8_t> consensusSeq;
        vector<int64_t> consToBb;
        vector<uint32_t> bbToCons(winLen, 0);
        consensusSeq.reserve(winLen + winLen / 8);
        consToBb.reserve(consensusSeq.capacity());

        auto emitInsertionSlot = [&](uint32_t slot) {
            if (insVotes[slot].empty()) return;
            // Majority insertion: most common inserted sequence, required on > half
            // of the reads spanning the slot.
            const string* best = nullptr; uint32_t bestN = 0;
            for (const auto& kv : insVotes[slot])
                if (kv.second > bestN) { bestN = kv.second; best = &kv.first; }
            if (best == nullptr) return;
            if (bestN * 2 <= insSlotCov[slot]) return; // not a consensus indel
            for (char ch : *best) {
                consensusSeq.push_back(uint8_t(ch));
                consToBb.push_back(-1);
            }
        };

        // IMPORTANT: substitution columns keep the BACKBONE base, never the
        // majority. At a 50/50 het the majority base is a coin flip between the
        // two haplotypes, so a majority consensus flips haplotype column-to-column
        // and both destroys the SNP signal (alt reads would match consensus) and
        // worsens drift (chimeric reference). We only stabilise the INDEL frame:
        // insert consensus insertion columns, but at every backbone column emit
        // the backbone base so alt reads still mismatch and the SNP is preserved.
        (void)baseVotes; (void)colCov;
        for (uint32_t i = 0; i < winLen; i++) {
            emitInsertionSlot(i); // consensus insertions before this column
            bbToCons[i] = uint32_t(consensusSeq.size());
            consensusSeq.push_back(bbSeqVec[windowBbBegin + i]);
            consToBb.push_back(int64_t(windowBbBegin + i));
        }
        emitInsertionSlot(winLen); // trailing insertions

        // Pass 2: realign each read PER SEGMENT (same short inter-anchor spans as
        // pass 1) against the CONSENSUS subsequence for that segment. Keeping the
        // segments short preserves band validity; the only change vs pass 1 is the
        // target (consensus, which carries the majority insertions) so per-read
        // gap ambiguity collapses onto a shared frame. SNPs are still called vs
        // the backbone base via consToBb.
        for (auto& prof : profiles) {
            if (prof.pinBb.size() < 2) continue;
            const OrientedReadId cOid = prof.oid;

            // Clear pass-1 results; refill from the consensus realignment.
            prof.snps.clear();
            prof.alignedCols.clear();
            prof.deletionRanges.clear();
            prof.noisyRanges.clear();
            KwNoiseTracker noise(noisyRegSlideWin, noisyRegMaxXgaps, prof.noisyRanges);

            // Realign one segment [bbSegBegin,bbSegEnd) of backbone against the
            // corresponding consensus columns [consBegin,consEnd) with the read's
            // [cSegBegin,cSegEnd). Emits SNPs/dels/cols in backbone coordinates.
            auto realignSeg = [&](uint32_t bbSegBegin, uint32_t bbSegEnd,
                                  uint32_t cSegBegin, uint32_t cSegEnd) {
                if (bbSegEnd <= bbSegBegin || cSegEnd <= cSegBegin) return;
                if (bbSegEnd <= windowBbBegin || bbSegBegin >= windowBbEnd) return;
                const uint32_t bb0 = max(bbSegBegin, windowBbBegin);
                const uint32_t bb1 = min(bbSegEnd, windowBbEnd);
                if (bb1 <= bb0) return;
                const uint32_t consBegin = bbToCons[bb0 - windowBbBegin];
                const uint32_t consEnd = bbToCons[(bb1 - 1) - windowBbBegin] + 1;
                if (consEnd <= consBegin) return;
                const uint32_t qLen = cSegEnd - cSegBegin;
                const uint32_t tLen = consEnd - consBegin;
                if (qLen == 0 || tLen == 0) return;

                static thread_local vector<uint8_t> query, target;
                query.resize(qLen);
                target.resize(tLen);
                for (uint32_t i = 0; i < qLen; i++)
                    query[i] = rds.getOrientedReadBase(cOid, cSegBegin + i).value;
                for (uint32_t i = 0; i < tLen; i++)
                    target[i] = consensusSeq[consBegin + i];

                const int band = abs(int(tLen) - int(qLen)) + kwBandSlack;
                ksw_extz_t ez;
                memset(&ez, 0, sizeof(ez));
                ksw_extd2_sse(nullptr, int(qLen), query.data(), int(tLen), target.data(),
                    5, mat, kswGapO1, kswGapE1, kswGapO2, kswGapE2, band, -1, 0, 0, &ez);

                if (ez.n_cigar > 0 && ez.cigar != nullptr) {
                    uint32_t qpos = cSegBegin;
                    uint32_t cpos = consBegin;
                    for (int ci = 0; ci < ez.n_cigar; ci++) {
                        const uint32_t len = ez.cigar[ci] >> 4;
                        const uint32_t op = ez.cigar[ci] & 0xf;
                        if (op == 0) {
                            for (uint32_t j = 0; j < len; j++) {
                                const int64_t bb = consToBb[cpos + j];
                                const uint32_t rp = qpos + j;
                                const uint8_t cb = rds.getOrientedReadBase(cOid, rp).value;
                                if (bb < 0) continue;
                                if (uint32_t(bb) < windowBbBegin || uint32_t(bb) >= windowBbEnd) continue;
                                const uint8_t bbBase = bbSeqVec[bb];
                                prof.alignedCols.push_back(KwAlignedCol{uint32_t(bb), rp, cb});
                                if (cb != bbBase && cb < 4 && bbBase < 4) {
                                    prof.snps.push_back(KwSnp{uint32_t(bb), cb});
                                    noise.observe(uint32_t(bb), 1, 1);
                                }
                            }
                            qpos += len; cpos += len;
                        } else if (op == 1) {
                            const int64_t bb = (cpos < consToBb.size()) ? consToBb[cpos] : -1;
                            if (bb >= 0) noise.observe(uint32_t(bb), 0, int(len));
                            qpos += len;
                        } else if (op == 2) {
                            for (uint32_t j = 0; j < len; j++) {
                                const int64_t bb = consToBb[cpos + j];
                                if (bb >= 0 && uint32_t(bb) >= windowBbBegin && uint32_t(bb) < windowBbEnd) {
                                    prof.deletionRanges.push_back({uint32_t(bb), uint32_t(bb) + 1});
                                    noise.observe(uint32_t(bb), 1, 1);
                                }
                            }
                            cpos += len;
                        } else {
                            cpos += len;
                        }
                    }
                }
                if (ez.cigar) free(ez.cigar);
            };

            // Same FULLSPAN tiling as pass 1 (anchor bodies included), realigned
            // against the consensus.
            for (size_t pi = 0; pi + 1 < prof.pinBb.size(); pi++)
                realignSeg(prof.pinBb[pi], prof.pinBb[pi + 1],
                           prof.pinC[pi], prof.pinC[pi + 1]);
            // Trailing anchor body.
            realignSeg(prof.pinBb.back(), prof.pinBb.back() + uint32_t(k),
                       prof.pinC.back(), prof.pinC.back() + uint32_t(k));
            noise.finish();

            // Dedup aligned columns by bbPos (segment joints can double-emit).
            sort(prof.alignedCols.begin(), prof.alignedCols.end(),
                [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos < b.bbPos; });
            prof.alignedCols.erase(
                unique(prof.alignedCols.begin(), prof.alignedCols.end(),
                    [](const KwAlignedCol& a, const KwAlignedCol& b) { return a.bbPos == b.bbPos; }),
                prof.alignedCols.end());
            // Dedup snps too.
            sort(prof.snps.begin(), prof.snps.end(),
                [](const KwSnp& a, const KwSnp& b) {
                    return a.bbPos < b.bbPos || (a.bbPos == b.bbPos && a.altBase < b.altBase); });
            prof.snps.erase(
                unique(prof.snps.begin(), prof.snps.end(),
                    [](const KwSnp& a, const KwSnp& b) { return a.bbPos == b.bbPos && a.altBase == b.altBase; }),
                prof.snps.end());

            if (!prof.deletionRanges.empty()) {
                sort(prof.deletionRanges.begin(), prof.deletionRanges.end());
                vector<pair<uint32_t,uint32_t>> merged;
                for (const auto& r : prof.deletionRanges) {
                    if (!merged.empty() && r.first <= merged.back().second)
                        merged.back().second = max(merged.back().second, r.second);
                    else merged.push_back(r);
                }
                prof.deletionRanges = move(merged);
            }
            if (!prof.noisyRanges.empty()) {
                sort(prof.noisyRanges.begin(), prof.noisyRanges.end());
                vector<pair<uint32_t,uint32_t>> merged;
                for (const auto& r : prof.noisyRanges) {
                    if (!merged.empty() && r.first <= merged.back().second)
                        merged.back().second = max(merged.back().second, r.second);
                    else merged.push_back(r);
                }
                prof.noisyRanges = move(merged);
            }
        }

        if (getenv("DINARA_KSW2_HET_DEBUG"))
            cout << "      ksw2twopass window=" << window.windowId
                 << " consensusLen=" << consensusSeq.size()
                 << " backboneLen=" << winLen
                 << " insertedCols=" << (consensusSeq.size() - winLen) << endl;
    }

    // Aggregate SNP candidates with strand counts (same as SNP path).
    struct SnpAccum { uint32_t fwd = 0; uint32_t rev = 0; uint32_t total = 0; };
    unordered_map<uint64_t, SnpAccum> snpCounts;
    auto snpKey = [](uint32_t pos, uint8_t alt) -> uint64_t {
        return (uint64_t(pos) << 8) | uint64_t(alt);
    };
    // Diagnostic: DINARA_KSW2_NONOISE=1 disables the CIGAR-density noise filter
    // so we can measure how many alt-haplotype reads it excludes in het clusters.
    const bool disableNoise = (getenv("DINARA_KSW2_NONOISE") != nullptr);

    unordered_set<uint32_t> snpPositions;
    for (const auto& prof : profiles) {
        const bool isFwd = (prof.oid.getStrand() == 0);
        for (const KwSnp& s : prof.snps) {
            if (!disableNoise && prof.isNoisy(s.bbPos)) continue;
            auto& acc = snpCounts[snpKey(s.bbPos, s.altBase)];
            if (isFwd) acc.fwd++; else acc.rev++;
            acc.total++;
            snpPositions.insert(s.bbPos);
        }
    }

    // Spanning coverage via sweep line (same as SNP path).
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
            spanningCount[pos] = uint32_t(max(running, 0)) + 1; // +1 backbone
        }
    }

    // Per-member allele call at a backbone position: -2 = ref (agrees with
    // backbone), a base 0-3 = that alt, or -1 if not covered / deleted / noisy.
    auto memberCall = [disableNoise](const KwMemberProfile& prof, uint32_t pos) -> int {
        if (pos < prof.bbCovBegin || pos >= prof.bbCovEnd) return -1;
        if (prof.isDeleted(pos)) return -1;
        if (!disableNoise && prof.isNoisy(pos)) return -1;
        for (const KwSnp& s : prof.snps) if (s.bbPos == pos) return int(s.altBase);
        return -2; // covered, agrees with backbone (ref)
    };

    // Per-column deletion support (spanning members that delete this backbone
    // position). Needed for the spanning denominator and the flank test, mirroring
    // abPOA's delSupport folded into the deletion skip edge.
    auto delSupportAt = [&](uint32_t pos) -> uint32_t {
        uint32_t d = 0;
        for (const auto& prof : profiles)
            if (pos >= prof.bbCovBegin && pos < prof.bbCovEnd && prof.isDeleted(pos)) d++;
        return d;
    };

    // abPOA-matching per-allele support cutoff. If the caller supplied a nonzero
    // hetMinSupport use it verbatim; otherwise auto-derive from coverageHet with
    // the same rule as testAbpoaMultiSegmentMSA (cut_bd=6, cut_rate=0.7, n_hap=2).
    int minSupport = 6;
    if (hetMinSupport > 0) {
        minSupport = int(hetMinSupport);
    } else {
        constexpr uint64_t cut_bd = 6, cut_rate_num = 7, cut_rate_den = 10, n_hap = 2;
        const uint64_t coverageHet = assemblerInfo.isOpen ?
            assemblerInfo->kmerDistributionInfo.coverageHet : invalid<uint64_t>;
        uint64_t base = 0;
        if (coverageHet != invalid<uint64_t> && coverageHet > 0) base = coverageHet / n_hap;
        uint64_t cc = (base * cut_rate_num) / cut_rate_den;
        if (cc < cut_bd) cc = cut_bd;
        minSupport = int(cc);
    }
    const double minVaf = hetMinVaf;

    // Classify passing SNPs with the SAME gates as detectWindowSnps (abPOA): ref
    // and each alt >= minSupport, alt VAF >= minVaf (no maxAf, no strand bias),
    // homopolymer/STR dropped ONLY when the matching flag is set. Group every
    // passing alt at one column into a single N-armed bubble.
    struct PassingAlt { uint8_t altBase; uint32_t altCov; uint32_t refCov; };
    struct PassingSite { uint32_t pos; uint32_t spanning; uint32_t refCount; vector<PassingAlt> alts; };
    vector<PassingSite> passingSites;

    vector<uint32_t> posList(snpPositions.begin(), snpPositions.end());
    sort(posList.begin(), posList.end());

    for (uint32_t pos : posList) {
        // Per-column base-allele tallies from the member pileup: refCount =
        // members agreeing with backbone (+1 for the backbone read itself),
        // altTotal[b] via snpCounts. spanning = ref + all alts + deletions,
        // exactly like abPOA's totalAlleleSupport + delSupport.
        uint32_t refCount = 1; // backbone is ref
        for (const auto& prof : profiles)
            if (memberCall(prof, pos) == -2) refCount++;

        uint32_t altTotal = 0;
        for (uint8_t b = 0; b < 4; b++) {
            auto it = snpCounts.find(snpKey(pos, b));
            if (it != snpCounts.end()) altTotal += it->second.total;
        }
        const uint32_t del = delSupportAt(pos);
        const uint32_t spanning = refCount + altTotal + del;
        if (spanning == 0) continue;

        if (refCount < uint32_t(minSupport)) continue; // ref allele must be present

        PassingSite site;
        site.pos = pos;
        site.spanning = spanning;
        site.refCount = refCount;

        for (uint8_t alt = 0; alt < 4; alt++) {
            auto it = snpCounts.find(snpKey(pos, alt));
            if (it == snpCounts.end()) continue;
            const uint32_t altCov = it->second.total;
            if (altCov < uint32_t(minSupport)) continue;

            const double af = double(altCov) / double(spanning);
            if (af < minVaf) continue;

            // Homopolymer (unit length 1) and STR (unit length 2..6) contexts are
            // gated independently, matching abPOA. Kept by default.
            KmVarKey vkey;
            vkey.pos = pos; vkey.type = KmVarType::Snp; vkey.altBase = alt;
            vkey.refLen = 1; vkey.altLen = 1;
            const bool inHomopolymer = kmIsRepeatUnitRange(bbSeq, bbLen, vkey, 0, 1, 1);
            const bool inStr         = kmIsRepeatUnitRange(bbSeq, bbLen, vkey, 0, 2, 6);
            if (hetDropHomopolymer && inHomopolymer) continue;
            if (hetDropRepeat && inStr) continue;

            site.alts.push_back({alt, altCov, refCount});
        }
        if (!site.alts.empty()) passingSites.push_back(move(site));
    }

    // Flank-linearity gate (majority-based reconstruction of abPOA's degree-1
    // predPrev/commonPred and commonSucc/succNext test). abPOA tolerates isolated
    // sequencing errors at the flanks -- they branch into low-coverage side nodes
    // that never form a clean bubble, so the MAJORITY (consensus) node keeps its
    // single in/out edge. We mirror that: a flank column breaks linearity only if
    // a competing base allele OR a deletion reaches minSupport there (i.e. the
    // flank is itself a real variant / indel boundary). Below-minSupport noise is
    // absorbed, exactly as the POA absorbs it. Columns p-2,p-1,p+1,p+2 are
    // checked so >=2 linear bases separate accepted SNPs (chainable homs).
    auto flanksLinear = [&](uint32_t pos) -> bool {
        if (pos < windowBbBegin + 2) return false;          // need p-2
        if (pos + 2 >= windowBbEnd) return false;           // need p+2
        const uint32_t flankCols[4] = {pos - 2, pos - 1, pos + 1, pos + 2};
        for (uint32_t c : flankCols) {
            for (uint8_t b = 0; b < 4; b++) {
                auto it = snpCounts.find(snpKey(c, b));
                if (it != snpCounts.end() && it->second.total >= uint32_t(minSupport))
                    return false;                            // real competing allele
            }
            if (delSupportAt(c) >= uint32_t(minSupport)) return false; // real deletion
        }
        return true;
    };

    // Build hetBubbles in backbone order, one per passing site, in the EXACT
    // format testAbpoaMultiSegmentMSA stages (arms + leadHom + hom). The
    // downstream plan/append/stage passes then treat them identically.
    // Per-position diagnostic (DINARA_KSW2_HET_DEBUG=1): for every raw SNP
    // position, dump the base tallies and which gate (if any) it fails, keyed by
    // backboneOffset so it can be cross-referenced against the abPOA log.
    if (getenv("DINARA_KSW2_HET_DEBUG")) {
        for (uint32_t pos : posList) {
            uint32_t refCount = 1;
            for (const auto& prof : profiles)
                if (memberCall(prof, pos) == -2) refCount++;
            uint32_t tally[4] = {0,0,0,0};
            for (uint8_t b = 0; b < 4; b++) {
                auto it = snpCounts.find(snpKey(pos, b));
                if (it != snpCounts.end()) tally[b] = it->second.total;
            }
            const uint32_t del = delSupportAt(pos);
            const uint32_t spanning = refCount + tally[0]+tally[1]+tally[2]+tally[3] + del;
            uint8_t topAlt = 0; uint32_t topCov = 0;
            for (uint8_t b = 0; b < 4; b++) if (tally[b] > topCov) { topCov = tally[b]; topAlt = b; }
            const double af = spanning ? double(topCov)/double(spanning) : 0.0;
            const bool okRef = refCount >= uint32_t(minSupport);
            const bool okAlt = topCov >= uint32_t(minSupport);
            const bool okVaf = af >= minVaf;
            const bool okFlank = flanksLinear(pos);
            cout << "      ksw2dbg window=" << window.windowId
                 << " off=" << (pos - windowBbBegin)
                 << " ref=" << refCount << " topAlt=" << int(topAlt)
                 << " topCov=" << topCov << " del=" << del
                 << " span=" << spanning << " vaf=" << af
                 << " minSup=" << minSupport
                 << (okRef?"":" [failRef]") << (okAlt?"":" [failAlt]")
                 << (okVaf?"":" [failVAF]") << (okFlank?"":" [failFlank]")
                 << endl;
        }
    }

    // One-shot per-profile probe for a single (window,offset): DINARA_KSW2_PROBE
    // = "<windowId>:<backboneOffset>". Explains why members are or aren't counted
    // at that column (bracketed? aligned col? noisy? deleted? base).
    if (const char* probe = getenv("DINARA_KSW2_PROBE")) {
        const string p(probe);
        const auto colon = p.find(':');
        if (colon != string::npos &&
            uint32_t(atoi(p.substr(0, colon).c_str())) == window.windowId) {
            const uint32_t poff = uint32_t(atoi(p.substr(colon + 1).c_str()));
            const uint32_t ppos = windowBbBegin + poff;
            uint32_t nBracket = 0, nAligned = 0, nNoisy = 0, nDel = 0;
            uint32_t baseTally[5] = {0,0,0,0,0};
            for (const auto& prof : profiles) {
                const bool brack = (ppos >= prof.bbCovBegin && ppos < prof.bbCovEnd);
                if (!brack) continue;
                nBracket++;
                if (prof.isNoisy(ppos)) nNoisy++;
                if (prof.isDeleted(ppos)) nDel++;
                const KwAlignedCol* c = prof.colAt(ppos);
                if (c) { nAligned++; baseTally[c->readBase < 4 ? c->readBase : 4]++; }
            }
            // Classify the bracketed-but-not-aligned reads: are they inside a
            // shared-anchor body (a pin covers [pin.bbPos, pin.bbPos+k) around
            // ppos) or in some other gap?
            uint32_t nInAnchorBody = 0, nGapNoAlign = 0;
            for (const auto& prof : profiles) {
                if (ppos < prof.bbCovBegin || ppos >= prof.bbCovEnd) continue;
                if (prof.colAt(ppos)) continue;          // aligned -> already counted
                if (prof.isDeleted(ppos)) continue;      // deletion -> counted
                // Not aligned, not deleted, yet bracketed: does a pin body cover it?
                bool inBody = false;
                for (const KwAlignedCol& c : prof.alignedCols) { (void)c; }
                // We do not store pins post-loop; approximate anchor-body test by
                // checking whether the nearest aligned columns straddle a >=k gap
                // around ppos (i.e. an un-aligned anchor region).
                const KwAlignedCol* lo = nullptr; const KwAlignedCol* hi = nullptr;
                for (const KwAlignedCol& c : prof.alignedCols) {
                    if (c.bbPos < ppos) lo = &c;
                    if (c.bbPos > ppos) { hi = &c; break; }
                }
                if (lo && hi && (hi->bbPos - lo->bbPos) >= uint32_t(k)) inBody = true;
                if (inBody) nInAnchorBody++; else nGapNoAlign++;
            }
            cout << "      ksw2probe window=" << window.windowId << " off=" << poff
                 << " pos=" << ppos << " profiles=" << profiles.size()
                 << " bracketed=" << nBracket << " aligned=" << nAligned
                 << " noisy=" << nNoisy << " deleted=" << nDel
                 << " inAnchorBody=" << nInAnchorBody << " gapNoAlign=" << nGapNoAlign
                 << " bases[A,C,G,T,N]=" << baseTally[0] << "," << baseTally[1]
                 << "," << baseTally[2] << "," << baseTally[3] << "," << baseTally[4]
                 << endl;
        }
    }

    uint32_t emitted = 0;
    for (const PassingSite& site : passingSites) {
        const uint32_t pos = site.pos;
        if (!flanksLinear(pos)) continue;   // no linear flanks -> no clean bubble

        const uint8_t predBase = bbSeqVec[pos - 1];

        AnchorWindow::HetBubble bubble;
        bubble.backboneOffset = pos - windowBbBegin;

        // Collect this site's base-allele member set (ref + kept alts) for the
        // two bracketing homs, and classify each spanning member for the arms.
        vector<OrientedReadId> homMemberReads;   // members spanning via a base allele
        vector<const KwMemberProfile*> homMemberProf;

        // Reference arm (arm 0): backbone + spanning members that agree with
        // backbone (ref) and are not carrying any alt here.
        AnchorWindow::HetAnchor refArm;
        refArm.backboneOffset = bubble.backboneOffset;
        refArm.predBase = predBase;
        refArm.alleleBase = bbSeqVec[pos];
        refArm.isRef = true;
        refArm.members.push_back({bbOid, pos - 1}); // backbone identity frame
        for (const auto& prof : profiles) {
            if (memberCall(prof, pos) != -2) continue;   // -2 == ref/agrees
            const KwAlignedCol* c = prof.colAt(pos - 1);
            if (!c) continue;                            // no recoverable rawPosition
            refArm.members.push_back({prof.oid, c->readPos});
            homMemberReads.push_back(prof.oid);
            homMemberProf.push_back(&prof);
        }
        // The site gate above used raw observation counts (site.refCount),
        // but the arm's real members are re-filtered by colAt (a member
        // without a recoverable rawPosition at pos-1 is dropped), which can
        // shrink the arm well below the count that passed the site gate --
        // in the extreme, down to just the backbone (size 1), which
        // appendHetAnchorPair asserts against (members.size() >= 2). Re-apply
        // minSupport to the PINNED member count so an allele that fell below
        // threshold after pinning does not reach the anchor graph at all
        // (mirrors the equivalent re-check in emitHetBubblesFromProfiles,
        // WindowHetProfiles.hpp, which this function does not otherwise call).
        if (refArm.members.size() >= size_t(minSupport))
            bubble.alleles.push_back(std::move(refArm));

        // Alternate arms.
        for (const PassingAlt& pa : site.alts) {
            AnchorWindow::HetAnchor arm;
            arm.backboneOffset = bubble.backboneOffset;
            arm.predBase = predBase;
            arm.alleleBase = pa.altBase;
            arm.isRef = false;
            for (const auto& prof : profiles) {
                if (memberCall(prof, pos) != int(pa.altBase)) continue;
                const KwAlignedCol* c = prof.colAt(pos - 1);
                if (!c) continue;
                arm.members.push_back({prof.oid, c->readPos});
                homMemberReads.push_back(prof.oid);
                homMemberProf.push_back(&prof);
            }
            // Same post-pinning re-check as the ref arm above: an empty-only
            // guard lets an arm through with as few as 1 member (below what
            // appendHetAnchorPair requires), if pinning dropped most of the
            // site gate's raw altCov.
            if (arm.members.size() < size_t(minSupport)) continue;
            bubble.alleles.push_back(std::move(arm));
        }
        if (bubble.alleles.size() < 2) continue;

        // Leading hom [predPrevBase, predBase] at predPrev (column pos-2). Its
        // members are the backbone plus every base-allele member with a
        // recoverable position there. predBackboneOffset = (pos-2) - windowBbBegin.
        {
            const uint32_t leadOff = (pos - 2) - windowBbBegin;
            AnchorWindow::HetAnchor leadHom;
            leadHom.backboneOffset = leadOff;
            leadHom.predBase = bbSeqVec[pos - 2];
            leadHom.alleleBase = predBase;         // linear next base
            leadHom.isRef = true;
            leadHom.members.push_back({bbOid, pos - 2});
            for (const KwMemberProfile* prof : homMemberProf) {
                const KwAlignedCol* c = prof->colAt(pos - 2);
                if (!c) continue;
                leadHom.members.push_back({prof->oid, c->readPos});
            }
            if (leadHom.members.size() > 1) {      // >1 => at least one non-backbone
                bubble.predBackboneOffset = leadOff;
                bubble.leadHom = std::move(leadHom);
            }
        }

        // Trailing hom [succBase, nextBase] at commonSucc (column pos+1).
        // succBackboneOffset = (pos+1) - windowBbBegin.
        {
            const uint32_t succOff = (pos + 1) - windowBbBegin;
            AnchorWindow::HetAnchor homAnchor;
            homAnchor.backboneOffset = succOff;
            homAnchor.predBase = bbSeqVec[pos + 1];  // succBase
            homAnchor.alleleBase = bbSeqVec[pos + 2]; // nextBase (linear)
            homAnchor.isRef = true;
            homAnchor.members.push_back({bbOid, pos + 1});
            for (const KwMemberProfile* prof : homMemberProf) {
                const KwAlignedCol* c = prof->colAt(pos + 1);
                if (!c) continue;
                homAnchor.members.push_back({prof->oid, c->readPos});
            }
            if (homAnchor.members.size() > 1) {
                bubble.succBackboneOffset = succOff;
                bubble.hom = std::move(homAnchor);
            }
        }

        window.hetBubbles.push_back(std::move(bubble));
        emitted++;
    }

    // NOTE: coincident-hom merging (sharedLeadFromBubble) is intentionally NOT
    // done here, for the same reason documented in emitHetBubblesFromProfiles
    // (WindowHetProfiles.hpp): it must run AFTER planning (main.cpp's Pass 1.5,
    // mergeWindowCoincidentHoms), because a bubble referenced here as "the
    // preceding bubble" can still be dropped by planning (primary-anchor
    // collision, a backward-drift hom member drop, or falling outside any
    // backbone interval) -- and nothing re-validates sharedLeadFromBubble
    // afterward. An earlier version of this function set it right here,
    // before planning ran; confirmed as a real bug (not hypothetical): on
    // real data it produced a bubble whose sharedLeadFromBubble pointed at a
    // since-dropped bubble (plannedInterval == -1), which
    // appendWindowHetAnchors (main.cpp) then asserted against. Leave
    // sharedLeadFromBubble at its default (-1) and let Pass 1.5 set it.

    cout << "    ksw2DetectHetBubbles bb=" << bbOid
         << " window=[" << windowBbBegin << "," << windowBbEnd << ")"
         << " reads=" << profiles.size()
         << " snpPositions=" << snpPositions.size()
         << " passingSites=" << passingSites.size()
         << " hetBubbles=" << emitted << endl;
    // Per-bubble offsets in the SAME frame abPOA prints (backboneOff = offset
    // into backboneCodes = bbPos - windowBbBegin), so the two engines can be
    // diffed column-for-column. windowId ties the line to the abPOA window.
    // Gated: only emitted under DINARA_KSW2_HET_DEBUG to keep normal runs quiet.
    if(getenv("DINARA_KSW2_HET_DEBUG"))
    for(const auto& b : window.hetBubbles) {
        cout << "      ksw2bubble window=" << window.windowId
             << " backboneOff=" << b.backboneOffset
             << " arms=" << b.alleles.size()
             << " ref=" << (b.alleles.empty() ? 0 : int(b.alleles[0].alleleBase));
        for(size_t ai = 1; ai < b.alleles.size(); ai++)
            cout << (ai > 1 ? "," : " alt=") << int(b.alleles[ai].alleleBase);
        cout << endl;
    }

    return emitted;
}
