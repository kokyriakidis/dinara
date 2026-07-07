#pragma once

/// @file WindowHetProfiles.hpp
/// @brief Engine-agnostic per-window het-bubble detection from a member pileup.
///
/// All het-detection engines (ksw2 star alignment and the POA engines) reduce a
/// window to the SAME intermediate representation: a set of per-member profiles
/// (KwMemberProfile) giving, for each read, the backbone columns it aligns, the
/// SNPs it carries, its coverage span, and its deleted/noisy ranges. Given those
/// profiles plus the backbone context, the het gating (min support, VAF,
/// homopolymer/STR, flank-linearity) and the emission of AnchorWindow::hetBubbles
/// (arms + leadHom + hom, in the exact format testAbpoaMultiSegmentMSA stages)
/// are identical. That shared tail lives here so the two engines differ ONLY in
/// how they populate the profiles.

#include "AnchorWindows.hpp"
#include "PhasingKmeansTypes.hpp"
#include "ReadId.hpp"
#include "invalid.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace dinara {

// A SNP observed on one member at a backbone position.
struct KwSnp {
    std::uint32_t bbPos;   // oriented backbone base position
    std::uint8_t altBase;  // 0=A,1=C,2=G,3=T
};

// One aligned (CIGAR/MSA M) column of a member against the backbone. bbPos is
// the oriented backbone base position; readBase is the member's base there
// (0-3); readPos is the member's ABSOLUTE oriented base position, which is
// exactly the rawPosition a k=2 anchor member stores at this column.
struct KwAlignedCol {
    std::uint32_t bbPos;
    std::uint32_t readPos;
    std::uint8_t readBase;
};

// Per-member outcome: SNPs carried, backbone coverage range, deleted positions,
// noisy ranges, and (for het-bubble emission) every aligned column so k=2 anchor
// rawPositions can be recovered at the flank columns.
struct KwMemberProfile {
    OrientedReadId oid;
    std::vector<KwSnp> snps;
    std::uint32_t bbCovBegin = 0; // oriented backbone coordinates
    std::uint32_t bbCovEnd = 0;
    std::vector<KwAlignedCol> alignedCols;  // sorted by bbPos, one per aligned pos

    // Read base + absolute read position at backbone position pos, or nullptr if
    // this member aligns no base there. Binary search over sorted alignedCols.
    const KwAlignedCol* colAt(std::uint32_t pos) const {
        auto it = std::lower_bound(alignedCols.begin(), alignedCols.end(), pos,
            [](const KwAlignedCol& c, std::uint32_t p) { return c.bbPos < p; });
        if (it != alignedCols.end() && it->bbPos == pos) return &*it;
        return nullptr;
    }

    // Sorted, non-overlapping [begin, end) backbone ranges deleted by this read.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> deletionRanges;
    // Sorted, non-overlapping [begin, end) backbone ranges flagged noisy.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> noisyRanges;

    static bool inRanges(const std::vector<std::pair<std::uint32_t, std::uint32_t>>& ranges,
                         std::uint32_t pos) {
        auto it = std::upper_bound(ranges.begin(), ranges.end(), pos,
            [](std::uint32_t p, const std::pair<std::uint32_t, std::uint32_t>& r) {
                return p < r.first; });
        if (it != ranges.begin()) {
            --it;
            if (pos >= it->first && pos < it->second) return true;
        }
        return false;
    }
    bool isDeleted(std::uint32_t pos) const { return inRanges(deletionRanges, pos); }
    bool isNoisy(std::uint32_t pos) const { return inRanges(noisyRanges, pos); }
};

// Emit AnchorWindow::hetBubbles from a completed member pileup. This is the
// shared tail of both engines: SNP aggregation, coverage sweep, per-allele
// support/VAF/homopolymer/STR gating, majority-based flank-linearity, and
// bubble emission (ref arm + alt arms + leadHom + trailing hom), followed by the
// coincident-hom merge. window.hetBubbles is cleared and repopulated.
//
// Inputs:
//   window        : the window; hetBubbles written here.
//   profiles      : per-member profiles built by the engine.
//   bbSeqVec      : full backbone base codes (index = oriented backbone pos).
//   windowBbBegin : backbone base offset of the window start (frame origin).
//   windowBbEnd   : backbone base offset one past the window end.
//   bbOid         : backbone oriented read id (identity frame for arm members).
//   k             : marker k (only used for context; arms are k=2).
//   hetMinVaf/hetMinSupport/hetDropHomopolymer/hetDropRepeat : het gates.
//   coverageHet   : diploid het coverage estimate, or invalid<uint64_t> if
//                   unavailable; used to auto-derive minSupport when
//                   hetMinSupport==0 (same rule as testAbpoaMultiSegmentMSA).
//   enginePrefix  : short tag for debug lines (e.g. "intervalpoa"/"ksw2").
//
// Returns the number of bubbles emitted.
inline std::uint32_t emitHetBubblesFromProfiles(
    AnchorWindow& window,
    const std::vector<KwMemberProfile>& profiles,
    const std::vector<std::uint8_t>& bbSeqVec,
    std::uint32_t windowBbBegin,
    std::uint32_t windowBbEnd,
    OrientedReadId bbOid,
    std::uint64_t k,
    double hetMinVaf,
    std::uint64_t hetMinSupport,
    bool hetDropHomopolymer,
    bool hetDropRepeat,
    std::uint64_t coverageHet,
    const char* enginePrefix)
{
    using std::vector;
    using std::pair;
    using std::uint8_t;
    using std::uint32_t;
    using std::uint64_t;
    using std::unordered_map;
    using std::unordered_set;
    using std::max;
    using std::min;
    using std::sort;

    const uint8_t* bbSeq = bbSeqVec.data();
    const uint32_t bbLen = uint32_t(bbSeqVec.size());

    const bool debug = (std::getenv("DINARA_HET_DEBUG") != nullptr);
    // Diagnostic: disable the per-member noise filter to measure its cost.
    const bool disableNoise = (std::getenv("DINARA_HET_NONOISE") != nullptr);

    // Aggregate SNP candidates with strand counts.
    struct SnpAccum { uint32_t fwd = 0; uint32_t rev = 0; uint32_t total = 0; };
    unordered_map<uint64_t, SnpAccum> snpCounts;
    auto snpKey = [](uint32_t pos, uint8_t alt) -> uint64_t {
        return (uint64_t(pos) << 8) | uint64_t(alt);
    };

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

    // Spanning coverage via sweep line.
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

    // Per-column deletion support (spanning members that delete this position).
    auto delSupportAt = [&](uint32_t pos) -> uint32_t {
        uint32_t d = 0;
        for (const auto& prof : profiles)
            if (pos >= prof.bbCovBegin && pos < prof.bbCovEnd && prof.isDeleted(pos)) d++;
        return d;
    };

    // abPOA-matching per-allele support cutoff. If the caller supplied a nonzero
    // hetMinSupport use it verbatim; otherwise auto-derive from coverageHet with
    // the same rule as testAbpoaMultiSegmentMSA (cut_bd=6, cut_rate=0.7,n_hap=2).
    int minSupport = 6;
    if (hetMinSupport > 0) {
        minSupport = int(hetMinSupport);
    } else {
        constexpr uint64_t cut_bd = 6, cut_rate_num = 7, cut_rate_den = 10, n_hap = 2;
        uint64_t base = 0;
        if (coverageHet != invalid<uint64_t> && coverageHet > 0) base = coverageHet / n_hap;
        uint64_t cc = (base * cut_rate_num) / cut_rate_den;
        if (cc < cut_bd) cc = cut_bd;
        minSupport = int(cc);
    }
    const double minVaf = hetMinVaf;

    // Classify passing SNPs with the SAME gates as detectWindowSnps (abPOA).
    struct PassingAlt { uint8_t altBase; uint32_t altCov; uint32_t refCov; };
    struct PassingSite { uint32_t pos; uint32_t spanning; uint32_t refCount; vector<PassingAlt> alts; };
    vector<PassingSite> passingSites;

    vector<uint32_t> posList(snpPositions.begin(), snpPositions.end());
    sort(posList.begin(), posList.end());

    for (uint32_t pos : posList) {
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

            KmVarKey vkey;
            vkey.pos = pos; vkey.type = KmVarType::Snp; vkey.altBase = alt;
            vkey.refLen = 1; vkey.altLen = 1;
            const bool inHomopolymer = kmIsRepeatUnitRange(bbSeq, bbLen, vkey, 0, 1, 1);
            const bool inStr         = kmIsRepeatUnitRange(bbSeq, bbLen, vkey, 0, 2, 6);
            if (hetDropHomopolymer && inHomopolymer) continue;
            if (hetDropRepeat && inStr) continue;

            site.alts.push_back({alt, altCov, refCount});
        }
        if (!site.alts.empty()) passingSites.push_back(std::move(site));
    }

    // Flank-linearity gate (majority-based reconstruction of abPOA's degree-1
    // predPrev/commonPred and commonSucc/succNext test). A flank column breaks
    // linearity only if a competing base allele OR a deletion reaches minSupport
    // there. Columns p-2,p-1,p+1,p+2 are checked so >=2 linear bases separate
    // accepted SNPs (chainable homs).
    auto flanksLinear = [&](uint32_t pos) -> bool {
        if (pos < windowBbBegin + 2) return false;
        if (pos + 2 >= windowBbEnd) return false;
        const uint32_t flankCols[4] = {pos - 2, pos - 1, pos + 1, pos + 2};
        for (uint32_t c : flankCols) {
            for (uint8_t b = 0; b < 4; b++) {
                auto it = snpCounts.find(snpKey(c, b));
                if (it != snpCounts.end() && it->second.total >= uint32_t(minSupport))
                    return false;
            }
            if (delSupportAt(c) >= uint32_t(minSupport)) return false;
        }
        return true;
    };

    window.hetBubbles.clear();

    // A k=2 anchor member stores a read position p; shasta2 reconstructs the
    // anchor k-mer as read[p..p+2] and requires EVERY member to carry the same
    // two bases. Pinning a member at the predecessor column (colAt(predPos)) is
    // only valid if that read actually spells [kmer0, kmer1] there: the base at
    // predPos equals kmer0, and the very next READ base (predPos's readPos + 1,
    // i.e. no insertion in between) is the aligned base at predPos+1 and equals
    // kmer1. A read with a mismatch/indel at the predecessor otherwise gets
    // pinned into an anchor whose k-mer it does not share, which shasta2 rejects
    // as an inconsistent k-mer. Returns the pred column to pin at, or nullptr.
    auto pinnedKmerCol = [](const KwMemberProfile& prof, uint32_t predPos,
                            uint8_t kmer0, uint8_t kmer1) -> const KwAlignedCol* {
        const KwAlignedCol* c0 = prof.colAt(predPos);
        if (!c0 || c0->readBase != kmer0) return nullptr;
        const KwAlignedCol* c1 = prof.colAt(predPos + 1);
        if (!c1 || c1->readPos != c0->readPos + 1 || c1->readBase != kmer1)
            return nullptr;
        return c0;
    };

    uint32_t emitted = 0;
    for (const PassingSite& site : passingSites) {
        const uint32_t pos = site.pos;
        if (!flanksLinear(pos)) continue;

        const uint8_t predBase = bbSeqVec[pos - 1];

        AnchorWindow::HetBubble bubble;
        bubble.backboneOffset = pos - windowBbBegin;

        vector<const KwMemberProfile*> homMemberProf;

        // Reference arm (arm 0).
        AnchorWindow::HetAnchor refArm;
        refArm.backboneOffset = bubble.backboneOffset;
        refArm.predBase = predBase;
        refArm.alleleBase = bbSeqVec[pos];
        refArm.isRef = true;
        refArm.members.push_back({bbOid, pos - 1});
        for (const auto& prof : profiles) {
            if (memberCall(prof, pos) != -2) continue;
            const KwAlignedCol* c = pinnedKmerCol(prof, pos - 1, predBase, refArm.alleleBase);
            if (!c) continue;
            refArm.members.push_back({prof.oid, c->readPos});
            homMemberProf.push_back(&prof);
        }
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
                const KwAlignedCol* c = pinnedKmerCol(prof, pos - 1, predBase, pa.altBase);
                if (!c) continue;
                arm.members.push_back({prof.oid, c->readPos});
                homMemberProf.push_back(&prof);
            }
            if (arm.members.empty()) continue;
            bubble.alleles.push_back(std::move(arm));
        }
        if (bubble.alleles.size() < 2) continue;

        // Leading hom [predPrevBase, predBase] at predPrev (column pos-2).
        {
            const uint32_t leadOff = (pos - 2) - windowBbBegin;
            AnchorWindow::HetAnchor leadHom;
            leadHom.backboneOffset = leadOff;
            leadHom.predBase = bbSeqVec[pos - 2];
            leadHom.alleleBase = predBase;
            leadHom.isRef = true;
            leadHom.members.push_back({bbOid, pos - 2});
            for (const KwMemberProfile* prof : homMemberProf) {
                const KwAlignedCol* c =
                    pinnedKmerCol(*prof, pos - 2, leadHom.predBase, leadHom.alleleBase);
                if (!c) continue;
                leadHom.members.push_back({prof->oid, c->readPos});
            }
            if (leadHom.members.size() > 1) {
                bubble.predBackboneOffset = leadOff;
                bubble.leadHom = std::move(leadHom);
            }
        }

        // Trailing hom [succBase, nextBase] at commonSucc (column pos+1).
        {
            const uint32_t succOff = (pos + 1) - windowBbBegin;
            AnchorWindow::HetAnchor homAnchor;
            homAnchor.backboneOffset = succOff;
            homAnchor.predBase = bbSeqVec[pos + 1];
            homAnchor.alleleBase = bbSeqVec[pos + 2];
            homAnchor.isRef = true;
            homAnchor.members.push_back({bbOid, pos + 1});
            for (const KwMemberProfile* prof : homMemberProf) {
                const KwAlignedCol* c =
                    pinnedKmerCol(*prof, pos + 1, homAnchor.predBase, homAnchor.alleleBase);
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

    // Coincident-hom merge (mirrors the abPOA Pass 1.5 condition).
    for (size_t i = 1; i < window.hetBubbles.size(); i++) {
        const auto& prev = window.hetBubbles[i - 1];
        auto& cur = window.hetBubbles[i];
        if (cur.predBackboneOffset != 0 && prev.succBackboneOffset != 0 &&
            cur.predBackboneOffset == prev.succBackboneOffset &&
            !cur.leadHom.members.empty() && !prev.hom.members.empty()) {
            cur.sharedLeadFromBubble = int64_t(i - 1);
        }
    }

    if (debug) {
        std::cout << "    " << enginePrefix << "DetectHetBubbles bb=" << bbOid
             << " window=[" << windowBbBegin << "," << windowBbEnd << ")"
             << " reads=" << profiles.size()
             << " snpPositions=" << snpPositions.size()
             << " passingSites=" << passingSites.size()
             << " hetBubbles=" << emitted << std::endl;
        for (const auto& b : window.hetBubbles) {
            std::cout << "      " << enginePrefix << "bubble window=" << window.windowId
                 << " backboneOff=" << b.backboneOffset
                 << " arms=" << b.alleles.size()
                 << " ref=" << (b.alleles.empty() ? 0 : int(b.alleles[0].alleleBase));
            for (size_t ai = 1; ai < b.alleles.size(); ai++)
                std::cout << (ai > 1 ? "," : " alt=") << int(b.alleles[ai].alleleBase);
            std::cout << std::endl;
        }
    }

    return emitted;
}

} // namespace dinara
