/// @file AssemblerWindowCigarMSA.cpp
/// @brief CIGAR-based SNP detection per anchor window.
/// Uses pairwise CIGARs from OverlapCigarStore to build variant profiles
/// against the backbone read, then classifies SNPs for het/hom determination.
/// Transitive reads (sharing anchors but not directly overlapping the backbone)
/// are projected through intermediary reads.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "PhasingKmeansTypes.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "OverlapCigarStore.hpp"
#include "invalid.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace dinara;

namespace {

constexpr uint32_t cwMinReadCoverage = 6;
constexpr uint32_t cwMinSnpAltSupport = 3;
constexpr uint32_t cwMinSnpRefSupport = 3;

// A variant event extracted from a CIGAR against the backbone.
struct CwVariant {
    uint32_t bbPos;       // backbone position
    KmVarType type;       // Snp, Insertion, Deletion
    uint8_t altBase;      // for SNPs: 0=A,1=C,2=G,3=T
    uint16_t len;         // refLen for DEL, altLen for INS, 1 for SNP
    string insSeq;        // inserted bases (INS only)
};

// Per-read variant profile within a window.
struct CwReadProfile {
    OrientedReadId oid;
    vector<CwVariant> variants; // sorted by bbPos
    bool isDirect;              // true if directly overlaps backbone
    uint32_t bbCovBegin = 0;    // backbone coverage range [begin, end)
    uint32_t bbCovEnd = 0;      // in oriented backbone coordinates

    // Backbone positions where this read has a deletion (no bases).
    // Sorted, non-overlapping [begin, end) ranges in oriented backbone frame.
    vector<pair<uint32_t, uint32_t>> deletionRanges;

    // Check if a backbone position falls within a deletion in this read.
    bool isDeleted(uint32_t pos) const {
        // Binary search for the first range whose end > pos.
        auto it = upper_bound(deletionRanges.begin(), deletionRanges.end(), pos,
            [](uint32_t p, const pair<uint32_t, uint32_t>& r) { return p < r.first; });
        if (it != deletionRanges.begin()) {
            --it;
            if (pos >= it->first && pos < it->second) return true;
        }
        return false;
    }
};

// A mapping from intermediary read positions to backbone positions,
// built from the intermediary's CIGAR against the backbone.
// Only match/mismatch positions have a valid mapping.
struct CwPosMap {
    // Sparse map: intermediary position → backbone position.
    // Only populated for match/mismatch ops (op 0 and 1).
    unordered_map<uint32_t, uint32_t> toBb;
};

} // anonymous namespace


// Compute read1 start position for CIGAR walking (matches kmCigarRead1Start).
static inline uint64_t cwCigarRead1Start(
    const Assembler& assembler, const AlignmentData& ad)
{
    if (ad.isSameStrand) return ad.ts;
    const uint32_t targetLen = uint32_t(
        assembler.getReads().getRead(ad.readIds[1]).baseCount);
    return targetLen - ad.te;
}

// Get a base from a read, handling reverse complement orientation.
// Matches kmGetBase from the k-means path.
static inline uint8_t cwGetBaseOriented(
    const Assembler& assembler, ReadId readId,
    uint32_t position, bool isReverseComplement)
{
    const auto sequence = assembler.getReads().getRead(readId);
    if (!isReverseComplement) return sequence[position].value;
    return sequence[sequence.baseCount - 1 - position].complement().value;
}


// Build a position map from an intermediary read's CIGAR against the backbone.
// Maps intermediary-read positions (in backbone's forward frame) to backbone positions.
static CwPosMap cwBuildPosMap(
    const Assembler& assembler,
    ReadId backboneReadId,
    uint32_t backboneLen,
    bool backboneIsStrand1,
    uint32_t alignmentId,
    uint32_t windowBbBegin,
    uint32_t windowBbEnd)
{
    CwPosMap posMap;
    const auto& cigarStore = assembler.getOverlapCigarStore();
    const auto& ad = assembler.alignmentData[alignmentId];
    const auto& info = ad.info;

    if (info.cigarOffset == uint32_t(-1) || info.cigarTokenCount == 0) return posMap;

    const bool bbIsRead0 = (ad.readIds[0] == backboneReadId);
    const ReadId partnerId = bbIsRead0 ? ad.readIds[1] : ad.readIds[0];
    const uint32_t partnerLen = uint32_t(
        assembler.getReads().getRead(partnerId).baseCount);
    const bool isRev = !ad.isSameStrand;
    const bool needsRc = !bbIsRead0 && isRev;
    // When B is read1 and strands differ, yk is in B's RC frame.
    // Normalize to B's forward frame for consistent lookup.
    const bool mirrorPartner = bbIsRead0 && isRev;

    const bool mirrorBb = backboneIsStrand1;

    cigarStore.forEachOpWithPositions(
        info.cigarOffset, info.cigarTokenCount,
        ad.qs, cwCigarRead1Start(assembler, ad),
        [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
            if (op == 0 || op == 1) { // Match or mismatch
                for (uint32_t b = 0; b < len; b++) {
                    uint32_t bbPos, partnerPos;
                    if (needsRc) {
                        uint32_t rc = uint32_t(yk) + b;
                        if (rc >= backboneLen) continue;
                        bbPos = backboneLen - 1 - rc;
                        partnerPos = uint32_t(xk) + b;
                    } else if (bbIsRead0) {
                        bbPos = uint32_t(xk) + b;
                        partnerPos = uint32_t(yk) + b;
                    } else {
                        bbPos = uint32_t(yk) + b;
                        partnerPos = uint32_t(xk) + b;
                    }
                    if (mirrorBb) bbPos = backboneLen - 1 - bbPos;
                    if (mirrorPartner) partnerPos = partnerLen - 1 - partnerPos;
                    if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;
                    posMap.toBb[partnerPos] = bbPos;
                }
            }
            // Insertions/deletions: no 1:1 mapping, skip.
        });

    return posMap;
}


// Parse a transitive read's CIGAR against an intermediary, projecting
// SNPs onto backbone coordinates via the intermediary's position map.
static void cwParseTransitiveCigar(
    const Assembler& assembler,
    ReadId intermediaryReadId,
    uint32_t intermediaryLen,
    uint32_t alignmentId,
    const CwPosMap& posMap,
    vector<CwVariant>& variants,
    uint32_t& covBeginOut,
    uint32_t& covEndOut)
{
    const auto& cigarStore = assembler.getOverlapCigarStore();
    const auto& ad = assembler.alignmentData[alignmentId];
    const auto& info = ad.info;

    if (info.cigarOffset == uint32_t(-1) || info.cigarTokenCount == 0) return;

    const bool intIsRead0 = (ad.readIds[0] == intermediaryReadId);
    const ReadId transitiveReadId = intIsRead0 ? ad.readIds[1] : ad.readIds[0];
    const bool isRev = !ad.isSameStrand;
    const bool needsRc = !intIsRead0 && isRev;

    // Track min/max backbone positions this transitive read maps to.
    uint32_t covMin = UINT32_MAX;
    uint32_t covMax = 0;

    cigarStore.forEachOpWithPositions(
        info.cigarOffset, info.cigarTokenCount,
        ad.qs, cwCigarRead1Start(assembler, ad),
        [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
            // For match and mismatch ops, update coverage from posMap lookups.
            if (op == 0 || op == 1) {
                // Check endpoints of this op for coverage range.
                for (uint32_t b = 0; b < len; b++) {
                    uint32_t intPos;
                    if (needsRc) {
                        uint32_t rc = uint32_t(yk) + b;
                        if (rc >= intermediaryLen) continue;
                        intPos = intermediaryLen - 1 - rc;
                    } else if (intIsRead0) {
                        intPos = uint32_t(xk) + b;
                    } else {
                        intPos = uint32_t(yk) + b;
                    }
                    auto it = posMap.toBb.find(intPos);
                    if (it == posMap.toBb.end()) continue;
                    uint32_t bbPos = it->second;
                    if (bbPos < covMin) covMin = bbPos;
                    if (bbPos + 1 > covMax) covMax = bbPos + 1;

                    // Only extract variants for mismatches.
                    if (op == 1) {
                        uint32_t tPos;
                        if (needsRc) {
                            tPos = uint32_t(xk) + b;
                        } else if (intIsRead0) {
                            tPos = uint32_t(yk) + b;
                        } else {
                            tPos = uint32_t(xk) + b;
                        }

                        uint8_t altBase;
                        if (intIsRead0) {
                            altBase = cwGetBaseOriented(assembler, transitiveReadId, tPos, isRev);
                        } else {
                            const auto seq = assembler.getReads().getRead(transitiveReadId);
                            altBase = seq[tPos].value;
                            if (isRev) altBase = uint8_t((~altBase) & 3);
                        }

                        variants.push_back({bbPos, KmVarType::Snp, altBase, 1, {}});
                    }
                }
            }
            // For transitive projection we only care about SNPs.
            // Indels relative to the intermediary are harder to project
            // and less reliable — skip them.
        });

    // Output coverage range.
    if (covMin < covMax) {
        covBeginOut = covMin;
        covEndOut = covMax;
    } else {
        covBeginOut = 0;
        covEndOut = 0;
    }

    sort(variants.begin(), variants.end(),
        [](const CwVariant& a, const CwVariant& b) { return a.bbPos < b.bbPos; });
}


// ============================================================================
// Parse CIGAR for one overlap against the backbone, extracting variants
// within [windowBbBegin, windowBbEnd) backbone coordinates.
// Follows the same CIGAR walking logic as kmParseCigars.
// ============================================================================
static void cwParseCigarForOverlap(
    const Assembler& assembler,
    ReadId backboneReadId,
    uint32_t backboneLen,
    bool backboneIsStrand1,
    uint32_t alignmentId,
    uint32_t windowBbBegin,
    uint32_t windowBbEnd,
    vector<CwVariant>& variants,
    uint32_t& covBeginOut,
    uint32_t& covEndOut,
    vector<pair<uint32_t, uint32_t>>& deletionRanges)
{
    const auto& cigarStore = assembler.getOverlapCigarStore();
    const auto& ad = assembler.alignmentData[alignmentId];
    const auto& info = ad.info;

    if (info.cigarOffset == uint32_t(-1) || info.cigarTokenCount == 0) return;

    // Determine which read is the backbone (read0 or read1 in the alignment).
    // AlignmentData always stores read0 on strand 0.
    const bool bbIsRead0 = (ad.readIds[0] == backboneReadId);
    const ReadId targetReadId = bbIsRead0 ? ad.readIds[1] : ad.readIds[0];
    const bool isRev = !ad.isSameStrand;
    const bool needsRc = !bbIsRead0 && isRev;

    // When the backbone oriented read is on strand 1, the window coordinates
    // are in the RC frame, but the CIGAR gives positions in read0-forward frame.
    // We convert bbPos to the oriented backbone's frame after extraction.
    const bool mirrorBb = backboneIsStrand1;

    // Track min/max backbone positions touched by this alignment
    // to determine the read's coverage range on the backbone.
    uint32_t covMin = UINT32_MAX;
    uint32_t covMax = 0;

    // Helper: update coverage range from a backbone-consuming op's endpoints.
    auto updateCovRange = [&](uint32_t rawStart, uint32_t rawLen) {
        if (rawLen == 0) return;
        uint32_t rawEnd = min(rawStart + rawLen, backboneLen);
        uint32_t s, e;
        if (needsRc) { s = backboneLen - rawEnd; e = backboneLen - rawStart; }
        else { s = rawStart; e = rawEnd; }
        if (mirrorBb) {
            uint32_t ms = backboneLen - e, me = backboneLen - s;
            s = ms; e = me;
        }
        if (s < covMin) covMin = s;
        if (e > covMax) covMax = e;
    };

    cigarStore.forEachOpWithPositions(
        info.cigarOffset, info.cigarTokenCount,
        ad.qs, cwCigarRead1Start(assembler, ad),
        [&](uint8_t op, uint32_t len, uint64_t xk, uint64_t yk) {
            // Update coverage range for backbone-consuming ops.
            // Op 0/1 (match/mismatch): both reads consumed.
            // Op 3: only read0 (xk) consumed — backbone-consuming iff bbIsRead0.
            // Op 2: only read1 (yk) consumed — backbone-consuming iff !bbIsRead0.
            if (op == 0 || op == 1) {
                uint32_t bbRaw = bbIsRead0 ? uint32_t(xk) : uint32_t(yk);
                updateCovRange(bbRaw, len);
            } else if ((op == 3 && bbIsRead0) || (op == 2 && !bbIsRead0)) {
                uint32_t bbRaw = (op == 3) ? uint32_t(xk) : uint32_t(yk);
                updateCovRange(bbRaw, len);
            }

            if (op == 1) { // Mismatch — one SNP per base
                for (uint32_t b = 0; b < len; b++) {
                    uint32_t bbPos, tPos;
                    if (needsRc) {
                        uint32_t rc = uint32_t(yk) + b;
                        if (rc >= backboneLen) continue;
                        bbPos = backboneLen - 1 - rc;
                        tPos = uint32_t(xk) + b;
                    } else if (bbIsRead0) {
                        bbPos = uint32_t(xk) + b;
                        tPos = uint32_t(yk) + b;
                    } else {
                        bbPos = uint32_t(yk) + b;
                        tPos = uint32_t(xk) + b;
                    }
                    // Convert to oriented backbone frame.
                    if (mirrorBb) bbPos = backboneLen - 1 - bbPos;
                    if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;
                    if (bbPos >= backboneLen) continue;

                    uint8_t altBase;
                    if (bbIsRead0) {
                        altBase = cwGetBaseOriented(assembler, targetReadId, tPos, isRev);
                    } else {
                        const auto seq = assembler.getReads().getRead(targetReadId);
                        altBase = seq[tPos].value;
                        if (isRev) altBase = uint8_t((~altBase) & 3);
                    }
                    // When mirrorBb is set, bbPos was mirrored to the backbone's
                    // RC frame. The alt base must also be complemented to match.
                    if (mirrorBb) altBase = uint8_t((~altBase) & 3);

                    variants.push_back({bbPos, KmVarType::Snp, altBase, 1, {}});
                }
            } else if (op == 2 || op == 3) { // Ins/Del
                bool bbConsumed = (op == 3 && bbIsRead0) || (op == 2 && !bbIsRead0);
                if (bbConsumed) {
                    // Deletion on backbone
                    uint32_t raw = (op == 3) ? uint32_t(xk) : uint32_t(yk);
                    uint32_t rawE = min(raw + len, backboneLen);
                    uint32_t s, e;
                    if (needsRc) { s = backboneLen - rawE; e = backboneLen - raw; }
                    else { s = raw; e = rawE; }
                    if (mirrorBb) {
                        uint32_t ms = backboneLen - e, me = backboneLen - s;
                        s = ms; e = me;
                    }
                    if (s >= windowBbEnd || e <= windowBbBegin) return;
                    s = max(s, windowBbBegin);
                    e = min(e, windowBbEnd);
                    if (s < backboneLen) {
                        variants.push_back({s, KmVarType::Deletion, 0, uint16_t(e - s), {}});
                        deletionRanges.push_back({s, e});
                    }
                } else {
                    // Insertion at backbone position
                    uint32_t anchor = bbIsRead0 ? uint32_t(xk) : uint32_t(yk);
                    if (needsRc) anchor = backboneLen - anchor;
                    if (mirrorBb) anchor = backboneLen - 1 - anchor;
                    if (anchor < windowBbBegin || anchor >= windowBbEnd) return;
                    if (anchor >= backboneLen) return;

                    string insSeq;
                    insSeq.reserve(len);
                    const uint32_t insStart = bbIsRead0 ? uint32_t(yk) : uint32_t(xk);
                    for (uint32_t b = 0; b < len; b++) {
                        uint8_t base;
                        if (bbIsRead0) {
                            base = cwGetBaseOriented(assembler, targetReadId, insStart + b, isRev);
                        } else {
                            const auto seq = assembler.getReads().getRead(targetReadId);
                            base = seq[insStart + b].value;
                            if (isRev) base = uint8_t((~base) & 3);
                        }
                        insSeq.push_back("ACGT"[base & 3]);
                    }
                    if (needsRc) reverse(insSeq.begin(), insSeq.end());
                    variants.push_back({anchor, KmVarType::Insertion, 0, uint16_t(len), move(insSeq)});
                }
            }
        });

    // Clamp coverage range to the window and output.
    if (covMin < covMax) {
        covBeginOut = max(covMin, windowBbBegin);
        covEndOut = min(covMax, windowBbEnd);
    } else {
        covBeginOut = windowBbBegin;
        covEndOut = windowBbBegin;
    }

    sort(variants.begin(), variants.end(),
        [](const CwVariant& a, const CwVariant& b) { return a.bbPos < b.bbPos; });

    // Sort and merge overlapping deletion ranges so isDeleted binary search
    // works correctly (it only checks the range with the largest start <= pos).
    sort(deletionRanges.begin(), deletionRanges.end());
    if (deletionRanges.size() > 1) {
        vector<pair<uint32_t, uint32_t>> merged;
        merged.push_back(deletionRanges[0]);
        for (size_t i = 1; i < deletionRanges.size(); i++) {
            if (deletionRanges[i].first <= merged.back().second) {
                merged.back().second = max(merged.back().second, deletionRanges[i].second);
            } else {
                merged.push_back(deletionRanges[i]);
            }
        }
        deletionRanges = move(merged);
    }
}


// ============================================================================
// Detect clean het SNPs in an anchor window using CIGAR-based variant parsing.
// Returns the number of SNPs passing strand bias and homopolymer/repeat filters.
// ============================================================================
uint32_t Assembler::cigarDetectSnpsInWindow(
    AnchorWindow& window,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys) const
{
    const Reads& rds = getReads();
    const auto& mkrs = *markers;
    const uint64_t k = assemblerInfo->k;
    const OrientedReadId bbOid = window.backboneOrientedReadId;
    const ReadId bbReadId = bbOid.getReadId();
    const uint32_t bbLen = uint32_t(rds.getRead(bbReadId).baseCount);
    const auto bbJ = journeys[bbOid];

    if (window.backboneEnd <= window.backboneBegin + 1) return 0;

    // Compute backbone coordinate range for this window.
    const uint32_t firstAnchorJP = window.backboneBegin;
    const uint32_t lastAnchorJP = window.backboneEnd - 1;
    if (firstAnchorJP >= bbJ.size() || lastAnchorJP >= bbJ.size()) return 0;

    const uint32_t firstOrd = anchors.getOrdinal(bbJ[firstAnchorJP], bbOid);
    const uint32_t lastOrd = anchors.getOrdinal(bbJ[lastAnchorJP], bbOid);
    if (firstOrd == invalid<uint32_t> || lastOrd == invalid<uint32_t>) return 0;

    const uint32_t windowBbBegin = mkrs[bbOid.getValue()][firstOrd].position;
    const uint32_t windowBbEnd = mkrs[bbOid.getValue()][lastOrd].position + uint32_t(k);
    if (windowBbEnd <= windowBbBegin) return 0;

    // Diagnostic: collect unique oriented reads from all backbone anchors.
    unordered_set<uint32_t> anchorReads; // oid.getValue()
    for (uint32_t jp = window.backboneBegin; jp < window.backboneEnd; jp++) {
        if (jp >= bbJ.size()) break;
        const Shasta2AnchorId anchorId = bbJ[jp];
        const auto anchor = anchors[anchorId];
        for (const auto& markerInfo : anchor) {
            if (markerInfo.orientedReadId != bbOid)
                anchorReads.insert(markerInfo.orientedReadId.getValue());
        }
    }

    // Reads from readIntervals (window planning).
    unordered_set<uint32_t> intervalReads;
    for (size_t ri = 1; ri < window.readIntervals.size(); ri++)
        intervalReads.insert(window.readIntervals[ri].orientedReadId.getValue());

    // Build alignment table lookup: for each read in the window, find its
    // alignment with the backbone.
    const auto& alnTable = getAlignmentTable();

    // Collect variant profiles for all reads in the window.
    vector<CwReadProfile> profiles;

    // Map from readId to profile index for transitive projection.
    unordered_map<uint32_t, size_t> readToProfile;

    // Build a lookup from partner ReadId → alignment IDs for the backbone's overlaps.
    // Keeps all chains (including secondary/deduped) so we can pick the one
    // whose backbone coverage overlaps the window region.
    unordered_map<uint32_t, vector<uint32_t>> partnerReadToAlnIds; // ReadId → [alignmentId]

    if (bbOid.getValue() < alnTable.size()) {
        const auto bbAlnIds = alnTable[bbOid.getValue()];
        for (size_t ai = 0; ai < bbAlnIds.size(); ai++) {
            const uint32_t alnId = bbAlnIds[ai];
            const auto& ad = alignmentData[alnId];
            // Accept secondary chains too — they may cover a different region.
            if (ad.isDeleted0() && ad.isDeleted1()) continue;
            if (ad.info.cigarOffset == uint32_t(-1) || ad.info.cigarTokenCount == 0) continue;

            const ReadId partnerId = (ad.readIds[0] == bbReadId)
                ? ad.readIds[1] : ad.readIds[0];
            partnerReadToAlnIds[uint32_t(partnerId)].push_back(alnId);
        }
    }

    // Helper: find the best alignment for a partner that covers the window region.
    // Returns the alignment ID, or uint32_t(-1) if none covers the window.
    auto findBestAlnForWindow = [&](ReadId partnerId) -> uint32_t {
        auto it = partnerReadToAlnIds.find(uint32_t(partnerId));
        if (it == partnerReadToAlnIds.end()) return uint32_t(-1);
        const auto& alnIds = it->second;

        // First try: find one that covers the window.
        for (uint32_t alnId : alnIds) {
            const auto& ad = alignmentData[alnId];
            bool bbIsR0 = (ad.readIds[0] == bbReadId);
            uint32_t bs = bbIsR0 ? ad.qs : ad.ts;
            uint32_t be = bbIsR0 ? ad.qe : ad.te;
            if (be > windowBbBegin && bs < windowBbEnd) return alnId;
        }
        // Fallback: return the first available chain.
        return alnIds.empty() ? uint32_t(-1) : alnIds[0];
    };

    // Process direct overlaps.
    for (size_t ri = 1; ri < window.readIntervals.size(); ri++) {
        const auto& interval = window.readIntervals[ri];
        const OrientedReadId oid = interval.orientedReadId;

        uint32_t bestAlnId = findBestAlnForWindow(oid.getReadId());
        if (bestAlnId == uint32_t(-1)) continue;

        CwReadProfile profile;
        profile.oid = oid;
        profile.isDirect = true;

        cwParseCigarForOverlap(
            *this, bbReadId, bbLen, bbOid.getStrand() == 1, bestAlnId,
            windowBbBegin, windowBbEnd, profile.variants,
            profile.bbCovBegin, profile.bbCovEnd,
            profile.deletionRanges);

        readToProfile[oid.getValue()] = profiles.size();
        profiles.push_back(move(profile));
    }

    // Collect transitive reads: window reads with no direct backbone overlap.
    // For each, find an intermediary read B that is already in profiles and
    // has a pairwise CIGAR with the transitive read C.
    {
        // Reads that need transitive projection.
        vector<size_t> transitiveIntervalIndices;
        for (size_t ri = 1; ri < window.readIntervals.size(); ri++) {
            const OrientedReadId oid = window.readIntervals[ri].orientedReadId;
            if (readToProfile.count(oid.getValue()) == 0)
                transitiveIntervalIndices.push_back(ri);
        }

        // Cache position maps for intermediary reads (built lazily).
        // Key: alignmentId of intermediary-vs-backbone overlap.
        unordered_map<uint32_t, CwPosMap> posMapCache;

        uint32_t transitiveAdded = 0;
        for (size_t ri : transitiveIntervalIndices) {
            const OrientedReadId cOid = window.readIntervals[ri].orientedReadId;
            const ReadId cReadId = cOid.getReadId();

            // Search C's alignment table for an intermediary B that is in profiles.
            if (cOid.getValue() >= alnTable.size()) continue;
            const auto cAlnIds = alnTable[cOid.getValue()];

            bool found = false;
            for (size_t ai = 0; ai < cAlnIds.size() && !found; ai++) {
                const uint32_t cAlnId = cAlnIds[ai];
                const auto& cAd = alignmentData[cAlnId];
                if (!cAd.keptByBothSides()) continue;
                if (cAd.info.cigarOffset == uint32_t(-1) || cAd.info.cigarTokenCount == 0) continue;

                // Find the partner (potential intermediary B).
                const ReadId bReadId = (cAd.readIds[0] == cReadId)
                    ? cAd.readIds[1] : cAd.readIds[0];

                // Check if B is in profiles (has a direct backbone overlap).
                // B could be on either strand in the window.
                OrientedReadId bOid0(bReadId, 0);
                OrientedReadId bOid1(bReadId, 1);
                size_t bProfileIdx = SIZE_MAX;
                auto bit = readToProfile.find(bOid0.getValue());
                if (bit != readToProfile.end()) bProfileIdx = bit->second;
                else {
                    bit = readToProfile.find(bOid1.getValue());
                    if (bit != readToProfile.end()) bProfileIdx = bit->second;
                }
                if (bProfileIdx == SIZE_MAX) continue;

                // Found intermediary B. Find its backbone alignment covering the window.
                uint32_t bBbAlnId = findBestAlnForWindow(bReadId);
                if (bBbAlnId == uint32_t(-1)) continue;

                auto pmIt = posMapCache.find(bBbAlnId);
                if (pmIt == posMapCache.end()) {
                    pmIt = posMapCache.emplace(bBbAlnId,
                        cwBuildPosMap(*this, bbReadId, bbLen,
                                     bbOid.getStrand() == 1, bBbAlnId,
                                     windowBbBegin, windowBbEnd)).first;
                }
                const CwPosMap& posMap = pmIt->second;
                if (posMap.toBb.empty()) continue;

                // Parse C-vs-B CIGAR, projecting SNPs onto backbone coordinates.
                const uint32_t bLen = uint32_t(rds.getRead(bReadId).baseCount);

                CwReadProfile profile;
                profile.oid = cOid;
                profile.isDirect = false;

                cwParseTransitiveCigar(
                    *this, bReadId, bLen, cAlnId, posMap, profile.variants,
                    profile.bbCovBegin, profile.bbCovEnd);

                readToProfile[cOid.getValue()] = profiles.size();
                profiles.push_back(move(profile));
                transitiveAdded++;
                found = true;
            }
        }

        if (transitiveAdded > 0) {
            cout << "      transitive: " << transitiveIntervalIndices.size()
                 << " candidates, " << transitiveAdded << " added" << endl;
        }
    }

    // Anchor-based segment comparison fallback for reads still not profiled.
    // These reads share anchors with the backbone in the window region but
    // have no usable CIGAR (direct or transitive). We compare the raw
    // sequence between consecutive shared anchors.
    {
        // Build backbone anchor set for the window.
        unordered_map<uint64_t, uint32_t> bbAnchorToJP; // anchorId → journey position
        for (uint32_t jp = window.backboneBegin; jp < window.backboneEnd; jp++) {
            if (jp >= bbJ.size()) break;
            bbAnchorToJP[uint64_t(bbJ[jp])] = jp;
        }

        uint32_t anchorFallbackAdded = 0;
        for (size_t ri = 1; ri < window.readIntervals.size(); ri++) {
            const auto& interval = window.readIntervals[ri];
            const OrientedReadId cOid = interval.orientedReadId;
            if (readToProfile.count(cOid.getValue()) != 0) continue;

            if (cOid.getValue() >= journeys.size()) continue;
            const auto cJ = journeys[cOid];

            // Find shared anchors between C and backbone in the window.
            struct SharedAnchor {
                uint32_t bbJP;     // backbone journey position
                uint32_t cOrd;     // C's ordinal for this anchor
                uint32_t bbPos;    // backbone base position
                uint32_t cPos;     // C's base position
            };
            vector<SharedAnchor> shared;

            for (uint32_t pos = interval.begin; pos < interval.end; pos++) {
                if (pos >= cJ.size()) break;
                auto it = bbAnchorToJP.find(uint64_t(cJ[pos]));
                if (it == bbAnchorToJP.end()) continue;

                uint32_t bbJP = it->second;
                uint32_t cOrd = anchors.getOrdinal(cJ[pos], cOid);
                uint32_t bbOrd = anchors.getOrdinal(bbJ[bbJP], bbOid);
                if (cOrd == invalid<uint32_t> || bbOrd == invalid<uint32_t>) continue;

                uint32_t cBasePos = mkrs[cOid.getValue()][cOrd].position;
                uint32_t bbBasePos = mkrs[bbOid.getValue()][bbOrd].position;

                shared.push_back({bbJP, cOrd, bbBasePos, cBasePos});
            }

            // Sort by backbone journey position.
            sort(shared.begin(), shared.end(),
                [](const SharedAnchor& a, const SharedAnchor& b) {
                    return a.bbJP < b.bbJP; });

            // Compare segments between consecutive shared anchors.
            CwReadProfile profile;
            profile.oid = cOid;
            profile.isDirect = false;

            if (shared.size() < 2) {
                // Not enough shared backbone anchors to form a segment —
                // skip this read (no variant or coverage information).
                continue;
            }

            // Coverage range: from first shared anchor to last shared anchor + k,
            // clamped to window.
            profile.bbCovBegin = max(shared.front().bbPos, windowBbBegin);
            profile.bbCovEnd = min(shared.back().bbPos + uint32_t(k), windowBbEnd);

            for (size_t si = 0; si + 1 < shared.size(); si++) {
                const auto& left = shared[si];
                const auto& right = shared[si + 1];

                // Segment on backbone: from left anchor position + k to right anchor position.
                uint32_t bbSegBegin = left.bbPos + uint32_t(k);
                uint32_t bbSegEnd = right.bbPos;
                if (bbSegEnd <= bbSegBegin) continue;

                // Segment on C: from left anchor position + k to right anchor position.
                uint32_t cSegBegin = left.cPos + uint32_t(k);
                uint32_t cSegEnd = right.cPos;
                if (cSegEnd <= cSegBegin) continue;

                uint32_t bbSegLen = bbSegEnd - bbSegBegin;
                uint32_t cSegLen = cSegEnd - cSegBegin;

                // Filter within window coordinate range.
                if (bbSegEnd <= windowBbBegin || bbSegBegin >= windowBbEnd) continue;

                if (bbSegLen == cSegLen) {
                    // Same length: compare base by base (pure SNP detection).
                    for (uint32_t offset = 0; offset < bbSegLen; offset++) {
                        uint32_t bbPos = bbSegBegin + offset;
                        if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;

                        uint8_t bbBase = rds.getOrientedReadBase(bbOid, bbPos).value;
                        uint8_t cBase = rds.getOrientedReadBase(cOid, cSegBegin + offset).value;

                        if (cBase != bbBase) {
                            profile.variants.push_back(
                                {bbPos, KmVarType::Snp, cBase, 1, {}});
                        }
                    }
                } else {
                    // Different length: indel region. Compare prefix from the
                    // left anchor and suffix from the right anchor to detect
                    // SNPs in the non-indel flanks. Each side gets half the
                    // shorter segment to avoid double-counting the middle.
                    uint32_t minLen = min(bbSegLen, cSegLen);
                    uint32_t halfLen = minLen / 2;

                    // Prefix: compare from left anchor forward.
                    for (uint32_t offset = 0; offset < halfLen; offset++) {
                        uint32_t bbPos = bbSegBegin + offset;
                        if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;

                        uint8_t bbBase = rds.getOrientedReadBase(bbOid, bbPos).value;
                        uint8_t cBase = rds.getOrientedReadBase(cOid, cSegBegin + offset).value;

                        if (cBase != bbBase) {
                            profile.variants.push_back(
                                {bbPos, KmVarType::Snp, cBase, 1, {}});
                        }
                    }

                    // Suffix: compare from right anchor backward.
                    for (uint32_t offset = 1; offset <= halfLen; offset++) {
                        uint32_t bbPos = bbSegEnd - offset;
                        if (bbPos < windowBbBegin || bbPos >= windowBbEnd) continue;

                        uint8_t bbBase = rds.getOrientedReadBase(bbOid, bbPos).value;
                        uint8_t cBase = rds.getOrientedReadBase(cOid, cSegEnd - offset).value;

                        if (cBase != bbBase) {
                            profile.variants.push_back(
                                {bbPos, KmVarType::Snp, cBase, 1, {}});
                        }
                    }
                }
            }

            readToProfile[cOid.getValue()] = profiles.size();
            profiles.push_back(move(profile));
            anchorFallbackAdded++;
        }

        if (anchorFallbackAdded > 0) {
            cout << "      anchorFallback: " << anchorFallbackAdded << " added" << endl;
        }
    }

    // Diagnostic: compare anchor reads vs interval reads vs profiled reads.
    {
        unordered_set<uint32_t> profiledReads;
        uint32_t directCount = 0, transitiveCount = 0;
        for (const auto& p : profiles) {
            profiledReads.insert(p.oid.getValue());
            if (p.isDirect) directCount++;
            else transitiveCount++;
        }

        // Anchor reads not in intervals.
        uint32_t anchorNotInInterval = 0;
        for (uint32_t v : anchorReads)
            if (intervalReads.count(v) == 0) anchorNotInInterval++;

        // Interval reads not profiled.
        uint32_t intervalNotProfiled = 0;
        for (uint32_t v : intervalReads)
            if (profiledReads.count(v) == 0) intervalNotProfiled++;

        // Anchor reads not profiled.
        uint32_t anchorNotProfiled = 0;
        for (uint32_t v : anchorReads)
            if (profiledReads.count(v) == 0) anchorNotProfiled++;

        cout << "      coverage: anchorReads=" << anchorReads.size()
             << " intervalReads=" << intervalReads.size()
             << " profiled=" << profiles.size()
             << " (direct=" << directCount
             << " transitive=" << transitiveCount << ")"
             << " anchorNotInInterval=" << anchorNotInInterval
             << " intervalNotProfiled=" << intervalNotProfiled
             << " anchorNotProfiled=" << anchorNotProfiled
             << endl;

    }

    if (profiles.size() + 1 < cwMinReadCoverage) return 0; // +1 for backbone

    // Build backbone sequence for repeat/homopolymer detection.
    vector<uint8_t> bbSeqVec(bbLen);
    for (uint32_t i = 0; i < bbLen; i++)
        bbSeqVec[i] = rds.getOrientedReadBase(bbOid, i).value;
    const uint8_t* bbSeq = bbSeqVec.data();

    KmPhasingOptions opts;

    // Aggregate SNP candidates across all profiles.
    // Key: backbone position → alt base → {fwdCount, revCount, totalCount}
    struct SnpAccum {
        uint32_t fwd = 0;
        uint32_t rev = 0;
        uint32_t total = 0;
    };
    // Map: bbPos → altBase → accumulator
    unordered_map<uint64_t, SnpAccum> snpCounts;

    auto snpKey = [](uint32_t pos, uint8_t alt) -> uint64_t {
        return (uint64_t(pos) << 8) | uint64_t(alt);
    };

    // Track which positions have any SNP.
    unordered_set<uint32_t> snpPositions;

    for (const auto& prof : profiles) {
        for (const auto& v : prof.variants) {
            if (v.type != KmVarType::Snp) continue;
            uint64_t key = snpKey(v.bbPos, v.altBase);
            auto& acc = snpCounts[key];
            if (prof.oid.getStrand() == 0) acc.fwd++;
            else acc.rev++;
            acc.total++;
            snpPositions.insert(v.bbPos);
        }
    }

    // Build sweep-line events from per-read coverage ranges.
    // Each read contributes +1 at bbCovBegin and -1 at bbCovEnd.
    // Sorting and sweeping gives the spanning count at any position.
    struct CovEvent {
        uint32_t pos;
        int delta; // +1 for begin, -1 for end
    };
    vector<CovEvent> covEvents;
    covEvents.reserve(profiles.size() * 2);
    for (const auto& prof : profiles) {
        if (prof.bbCovBegin < prof.bbCovEnd) {
            covEvents.push_back({prof.bbCovBegin, +1});
            covEvents.push_back({prof.bbCovEnd, -1});
        }
    }
    // Sort by position; at the same position, +1 (begin) before -1 (end)
    // so a read starting at pos p is counted as spanning pos p.
    sort(covEvents.begin(), covEvents.end(),
        [](const CovEvent& a, const CovEvent& b) {
            return a.pos < b.pos || (a.pos == b.pos && a.delta > b.delta);
        });

    // Sort SNP positions for sweep.
    vector<uint32_t> sortedSnpPositions(snpPositions.begin(), snpPositions.end());
    sort(sortedSnpPositions.begin(), sortedSnpPositions.end());

    // Sweep: compute spanning read count at each SNP position.
    unordered_map<uint32_t, uint32_t> spanningCount;
    spanningCount.reserve(sortedSnpPositions.size());
    {
        int running = 0;
        size_t ei = 0; // event index
        for (uint32_t pos : sortedSnpPositions) {
            while (ei < covEvents.size() && covEvents[ei].pos <= pos) {
                running += covEvents[ei].delta;
                ei++;
            }
            spanningCount[pos] = uint32_t(max(running, 0));
        }
    }

    // Classify SNPs and collect passing het SNPs for phasing.
    struct PassingSnp {
        uint32_t pos;
        uint8_t altBase;
        uint32_t altCov;
        uint32_t refCov;
        uint32_t spanning;
        uint32_t fwd;
        uint32_t rev;
    };
    vector<PassingSnp> passingSnps;

    uint32_t cleanHetSnps = 0;

    // Filter counters.
    uint32_t passAltSupport = 0;
    uint32_t failRefSupport = 0;
    uint32_t failAf = 0;
    uint32_t failStrandBias = 0;
    uint32_t failHomopolymer = 0;

    for (uint32_t pos : snpPositions) {
        const uint32_t spanning = spanningCount[pos];
        if (spanning == 0) continue;

        // Count reads with a deletion at this position.
        uint32_t delCount = 0;
        for (const auto& prof : profiles) {
            if (pos >= prof.bbCovBegin && pos < prof.bbCovEnd && prof.isDeleted(pos))
                delCount++;
        }
        // Effective spanning: reads that actually have bases at this position.
        const uint32_t effSpanning = (spanning > delCount) ? spanning - delCount : 0;
        if (effSpanning == 0) continue;

        // For each alt allele at this position.
        for (uint8_t alt = 0; alt < 4; alt++) {
            uint64_t key = snpKey(pos, alt);
            auto it = snpCounts.find(key);
            if (it == snpCounts.end()) continue;
            const auto& acc = it->second;

            if (acc.total < cwMinSnpAltSupport) continue;
            passAltSupport++;

            // Ref support = reads with bases at this position minus alt reads.
            const uint32_t refCov = (effSpanning > acc.total) ? effSpanning - acc.total : 0;
            if (refCov < cwMinSnpRefSupport) { failRefSupport++; continue; }

            const double af = double(acc.total) / double(effSpanning);
            if (af < opts.minAf || af > opts.maxAf) { failAf++; continue; }

            // Strand bias.
            const int expected = int(acc.total) / 2;
            if (expected > 0) {
                const double p = kmFisherExactTwoTail(
                    int(acc.fwd), int(acc.rev), expected, expected);
                if (p < opts.strandBiasPval) { failStrandBias++; continue; }
            }

            // Homopolymer / repeat context.
            KmVarKey vkey;
            vkey.pos = pos;
            vkey.type = KmVarType::Snp;
            vkey.altBase = alt;
            vkey.refLen = 1;
            vkey.altLen = 1;
            if (kmIsHomopolymer(bbSeq, bbLen, vkey, 0) ||
                kmIsRepeatRegion(bbSeq, bbLen, vkey, 0)) { failHomopolymer++; continue; }

            passingSnps.push_back({pos, alt, acc.total, refCov, effSpanning,
                                   acc.fwd, acc.rev});
            cleanHetSnps++;
        }
    }

    cout << "    cigarDetectSnps bb=" << bbOid
         << " window=[" << windowBbBegin << "," << windowBbEnd << ")"
         << " reads=" << profiles.size()
         << " snpPositions=" << snpPositions.size()
         << " passAltSupport=" << passAltSupport
         << " failRefSupport=" << failRefSupport
         << " failAF=" << failAf
         << " failStrandBias=" << failStrandBias
         << " failHomopolymer=" << failHomopolymer
         << " cleanHetSnps=" << cleanHetSnps << endl;

    // Store passing het SNPs in the window, with per-read alt/ref lists.
    window.hetSnps.resize(passingSnps.size());
    for (size_t i = 0; i < passingSnps.size(); i++) {
        auto& hs = window.hetSnps[i];
        hs.bbPos = passingSnps[i].pos;
        hs.refBase = bbSeqVec[passingSnps[i].pos];
        hs.altBase = passingSnps[i].altBase;
        hs.altCov = passingSnps[i].altCov;
        hs.refCov = passingSnps[i].refCov;
        hs.spanning = passingSnps[i].spanning;

        // Backbone is always ref at every position.
        if (passingSnps[i].pos >= windowBbBegin && passingSnps[i].pos < windowBbEnd)
            hs.refReads.push_back(bbOid);

        for (const auto& prof : profiles) {
            // Skip reads that don't cover this position.
            if (passingSnps[i].pos < prof.bbCovBegin ||
                passingSnps[i].pos >= prof.bbCovEnd) continue;
            // Skip reads with a deletion at this position.
            if (prof.isDeleted(passingSnps[i].pos)) continue;

            // Check if this read has any SNP variant at this position.
            // A read with a different alt allele is neither ref nor alt
            // for this specific variant — skip it entirely.
            bool hasThisAlt = false;
            bool hasOtherAlt = false;
            for (const auto& v : prof.variants) {
                if (v.type == KmVarType::Snp && v.bbPos == passingSnps[i].pos) {
                    if (v.altBase == passingSnps[i].altBase) {
                        hasThisAlt = true;
                    } else {
                        hasOtherAlt = true;
                    }
                    break;
                }
            }
            if (hasThisAlt) hs.altReads.push_back(prof.oid);
            else if (!hasOtherAlt) hs.refReads.push_back(prof.oid);
            // else: read has a different alt allele — excluded from both lists
        }
    }

    // ── Hifiasm-style DP chaining + transitive closure ──
    // 1. Build per-site evidence: which reads are ref/alt at each passing SNP
    // 2. DP chain to select consistent het SNP sites
    // 3. Transitive closure to classify overlaps as cis/trans
    // ── Hifiasm-style DP chaining + transitive closure ──
    // Port of gen_rphase_dp0_single_path and generate_haplotypes_naive_HiFi.

    const uint32_t numOv = uint32_t(profiles.size());

    // Per-site evidence: for each passing SNP, which reads are ref/alt.
    struct SiteEvidence {
        vector<int8_t> allele; // per-profile: 0=ref, 1=alt, -1=missing
        uint32_t nRef = 0;
        uint32_t nAlt = 0;
    };

    if (cleanHetSnps > 0 && numOv >= 2) {

        // Sort passing SNPs by position.
        sort(passingSnps.begin(), passingSnps.end(),
            [](const PassingSnp& a, const PassingSnp& b) {
                return a.pos < b.pos || (a.pos == b.pos && a.altBase < b.altBase);
            });

        const uint32_t numSites = uint32_t(passingSnps.size());

        // Adjacent-site filter (hifiasm: skip sites at distance 1).
        // A site is dropped if its neighbor is at distance 1 on either side.
        vector<bool> adjDrop(numSites, false);
        for (uint32_t si = 0; si < numSites; si++) {
            if (si > 0 && passingSnps[si].pos == passingSnps[si-1].pos + 1) {
                adjDrop[si] = true;
                adjDrop[si-1] = true;
            }
            if (si + 1 < numSites && passingSnps[si].pos + 1 == passingSnps[si+1].pos) {
                adjDrop[si] = true;
                adjDrop[si+1] = true;
            }
        }

        // Build per-site evidence.
        // profiles[] does NOT include the backbone read. We add it as a
        // virtual entry at index numOv (allele=0 at all covered sites).
        // This matches hifiasm where the query read is part of the evidence.
        const uint32_t bbIdx = numOv; // virtual backbone index
        const uint32_t numEvidence = numOv + 1; // profiles + backbone

        vector<SiteEvidence> sites(numSites);
        for (uint32_t si = 0; si < numSites; si++) {
            sites[si].allele.assign(numEvidence, -1);
            // Backbone is ref at all sites it covers.
            if (passingSnps[si].pos >= windowBbBegin && passingSnps[si].pos < windowBbEnd) {
                sites[si].allele[bbIdx] = 0;
                sites[si].nRef++;
            }
        }

        for (uint32_t oi = 0; oi < numOv; oi++) {
            const auto& prof = profiles[oi];

            // Build set of this read's SNP variants.
            unordered_map<uint64_t, uint8_t> readSnps;
            for (const auto& v : prof.variants) {
                if (v.type != KmVarType::Snp) continue;
                readSnps[snpKey(v.bbPos, v.altBase)] = 1;
            }

            for (uint32_t si = 0; si < numSites; si++) {
                uint32_t pos = passingSnps[si].pos;
                if (pos < prof.bbCovBegin || pos >= prof.bbCovEnd) continue;
                if (prof.isDeleted(pos)) continue;

                uint64_t sk = snpKey(pos, passingSnps[si].altBase);
                if (readSnps.count(sk)) {
                    sites[si].allele[oi] = 1; // alt
                    sites[si].nAlt++;
                } else {
                    // Check for different alt allele at same position.
                    bool hasDiffAlt = false;
                    for (uint8_t a = 0; a < 4; a++) {
                        if (a == passingSnps[si].altBase) continue;
                        if (readSnps.count(snpKey(pos, a))) {
                            hasDiffAlt = true;
                            break;
                        }
                    }
                    if (hasDiffAlt) {
                        sites[si].allele[oi] = -1; // different alt = ambiguous
                    } else {
                        sites[si].allele[oi] = 0; // ref
                        sites[si].nRef++;
                    }
                }
            }
        }

        // ── comput_sc_rphase port ──
        // Two sites are consistent if shared reads agree on allele pattern:
        // need ≥1 read that is ref at both AND ≥1 read that is alt at both.
        // Any read that is ambiguous (allele=-1) at either site is skipped.
        // Any read with fi != fj (ref at one, alt at other) → INT64_MIN.
        auto computScRphase = [&](uint32_t si, uint32_t sj) -> int64_t {
            if (passingSnps[si].pos == passingSnps[sj].pos) return INT64_MIN;

            int nn0 = 0; // reads that are ref at both sites
            int nn1 = 0; // reads that are alt at both sites

            for (uint32_t oi = 0; oi < numEvidence; oi++) {
                int8_t fi = sites[si].allele[oi];
                int8_t fj = sites[sj].allele[oi];
                if (fi < 0 || fj < 0) continue; // missing/ambiguous at either
                if (fi != fj) return INT64_MIN;  // inconsistent
                if (fi == 0) nn0++;
                else nn1++;
            }

            if (nn0 > 0 && nn1 > 0) return 1;
            return INT64_MIN;
        };

        // ── gen_rphase_dp0_single_path port ──
        // DP to find best chain of consistent sites.
        // cc = min ref support threshold.
        // hifiasm: cc = het_cov * cut_rate (cut_rate=0.7), min cut_bd (=6).
        // coveragePeak is homozygous coverage; het_cov = coveragePeak / 2.
        const uint64_t coveragePeak = assemblerInfo->kmerDistributionInfo.coveragePeak;
        const uint64_t hetCov = coveragePeak / 2;
        uint32_t cc = uint32_t(hetCov * 0.7);
        if (cc < 6) cc = 6;

        // Filter out adjacent-dropped sites for DP.
        vector<uint32_t> dpSites;
        for (uint32_t si = 0; si < numSites; si++) {
            if (!adjDrop[si]) dpSites.push_back(si);
        }
        const uint32_t dpN = uint32_t(dpSites.size());

        vector<int32_t> dpF(dpN, 1);   // best score ending here
        vector<int32_t> dpP(dpN, -1);  // predecessor index in dpSites
        vector<bool> dpUsed(dpN, false);

        for (uint32_t i = 0; i < dpN; i++) {
            int32_t maxF = 1;
            int32_t maxJ = -1;
            for (int32_t j = int32_t(i) - 1; j >= 0; j--) {
                int64_t sc = computScRphase(dpSites[i], dpSites[j]);
                if (sc == INT64_MIN) continue;
                int32_t candidate = dpF[j] + int32_t(sc);
                if (candidate > maxF) {
                    maxF = candidate;
                    maxJ = j;
                }
            }
            dpF[i] = maxF;
            dpP[i] = maxJ;
        }

        // Extract disjoint paths by greedy traceback from highest scores.
        // Sort indices by score descending.
        vector<uint32_t> dpOrder(dpN);
        iota(dpOrder.begin(), dpOrder.end(), 0);
        sort(dpOrder.begin(), dpOrder.end(),
            [&](uint32_t a, uint32_t b) { return dpF[a] > dpF[b]; });

        // Greedy path extraction (hifiasm: extract non-overlapping chains).
        struct DpPath {
            vector<uint32_t> siteIndices; // indices into dpSites
        };
        vector<DpPath> paths;

        for (uint32_t idx = 0; idx < dpN; idx++) {
            uint32_t start = dpOrder[idx];
            if (dpUsed[start]) continue;

            DpPath path;
            for (int32_t k = int32_t(start); k >= 0 && !dpUsed[k]; ) {
                path.siteIndices.push_back(k);
                dpUsed[k] = true;
                k = dpP[k];
            }
            if (path.siteIndices.empty()) continue;

            // Reverse so sites are in position order.
            reverse(path.siteIndices.begin(), path.siteIndices.end());
            paths.push_back(std::move(path));
        }

        // Score each path's sites (hifiasm logic):
        // - Chains with length > 1: mark score=1 if occ_0 >= cc, else score=-1
        // - Singletons: mark score=1 only if occ_0 >= cc (and not homopolymer,
        //   already filtered earlier)
        // - Dense clusters (all sites within distance 8): treat as singleton
        //   with stricter filtering.
        vector<int8_t> siteScore(numSites, -1); // -1 = rejected, 1 = validated

        for (const auto& path : paths) {
            const uint32_t pathLen = uint32_t(path.siteIndices.size());

            // Check if all sites in path are within distance 8 of each other
            // (hifiasm: adjacent sites within distance 8 → stricter filtering).
            bool isDense = true;
            if (pathLen > 1) {
                for (uint32_t pi = 1; pi < pathLen; pi++) {
                    uint32_t si0 = dpSites[path.siteIndices[pi-1]];
                    uint32_t si1 = dpSites[path.siteIndices[pi]];
                    if (passingSnps[si1].pos - passingSnps[si0].pos > 8) {
                        isDense = false;
                        break;
                    }
                }
            }

            int8_t plus;
            if (pathLen > 1 && !isDense) {
                plus = 1; // multi-site chain, not dense → validated
            } else {
                plus = -1; // singleton or dense cluster → needs extra check
            }

            for (uint32_t pi = 0; pi < pathLen; pi++) {
                uint32_t si = dpSites[path.siteIndices[pi]];
                if (sites[si].nRef >= cc) {
                    siteScore[si] = plus;
                } else {
                    siteScore[si] = -1;
                }
            }

            // Dense cluster: if any site in a run of consecutive positions
            // is rejected, reject the whole run (hifiasm behavior).
            if (isDense && pathLen > 1) {
                bool anyBad = false;
                for (uint32_t pi = 0; pi < pathLen; pi++) {
                    uint32_t si = dpSites[path.siteIndices[pi]];
                    if (siteScore[si] == -1) { anyBad = true; break; }
                }
                if (anyBad) {
                    for (uint32_t pi = 0; pi < pathLen; pi++) {
                        siteScore[dpSites[path.siteIndices[pi]]] = -1;
                    }
                }
            }
        }

        // Collect validated sites.
        vector<uint32_t> validSites;
        for (uint32_t si = 0; si < numSites; si++) {
            if (siteScore[si] == 1) validSites.push_back(si);
        }

        cout << "    DP chain: " << validSites.size() << " / " << numSites
             << " sites validated (" << paths.size() << " paths, "
             << dpN << " after adj filter)" << endl;

        // ── generate_haplotypes_naive_HiFi port ──
        // Transitive closure to classify reads as cis/trans.

        if (validSites.empty()) {
            // No validated sites — all reads in one cluster.
            window.readClusters.resize(1);
            window.readClusters[0].reserve(numOv);
            for (const auto& prof : profiles) {
                window.readClusters[0].push_back(prof.oid);
            }
            return cleanHetSnps;
        }

        const uint32_t nValid = uint32_t(validSites.size());

        // Step 1: Per-read, count alt alleles at validated sites.
        // Build (overlapIdx, altCount) sorted by altCount descending.
        // This corresponds to hifiasm's snp_srt: sort overlaps by how many
        // alt alleles they carry at validated sites (most informative first).
        struct ReadAltInfo {
            uint32_t oi;
            uint32_t nAlt = 0;   // alt alleles at validated sites
            uint32_t nObs = 0;   // total observations at validated sites
        };
        vector<ReadAltInfo> readAlts(numOv);
        for (uint32_t oi = 0; oi < numOv; oi++) {
            readAlts[oi].oi = oi;
            for (uint32_t vi = 0; vi < nValid; vi++) {
                int8_t a = sites[validSites[vi]].allele[oi];
                if (a < 0) continue;
                readAlts[oi].nObs++;
                if (a == 1) readAlts[oi].nAlt++;
            }
        }

        // Sort by nAlt descending (most informative first).
        // hifiasm: radix_sort_bc64 on ((uint32_t)-1 - o) << 32 | l
        vector<uint32_t> readOrder(numOv);
        iota(readOrder.begin(), readOrder.end(), 0);
        sort(readOrder.begin(), readOrder.end(),
            [&](uint32_t a, uint32_t b) {
                return readAlts[a].nAlt > readAlts[b].nAlt;
            });

        // Per-read: is_match. 0=unset, 1=cis, 2=trans.
        vector<uint8_t> isMatch(numOv, 0);

        // Per-site: score tracking for "not real allele" adjustment.
        // We make a mutable copy of nRef for decrementing.
        vector<uint32_t> siteRefCov(nValid);
        for (uint32_t vi = 0; vi < nValid; vi++) {
            siteRefCov[vi] = sites[validSites[vi]].nRef;
        }

        // Step 2: First pass — seed trans from most-informative overlaps.
        // hifiasm: for each overlap sorted by alt count, if it has alt alleles
        // at validated sites with sufficient support, mark as trans.
        // Then promote those sites (score=1) and decrement ref support
        // ("not real allele" adjustment).
        for (uint32_t ri = 0; ri < numOv; ri++) {
            uint32_t oi = readOrder[ri];
            if (readAlts[oi].nAlt == 0) break; // sorted desc, no more alts
            if (readAlts[oi].nObs == 0) continue;

            // Count alt alleles at validated sites with sufficient support.
            // hifiasm: occ_0 >= s_hap_cov && occ_1 >= infor_cov
            // We use nRef >= 2 && nAlt >= 2 (already filtered in passingSnps).
            uint32_t qualifiedAlts = 0;
            for (uint32_t vi = 0; vi < nValid; vi++) {
                uint32_t si = validSites[vi];
                if (sites[si].allele[oi] != 1) continue;
                if (siteRefCov[vi] >= 2 && sites[si].nAlt >= 2) {
                    qualifiedAlts++;
                }
            }
            if (qualifiedAlts == 0) continue;

            // Mark as trans.
            isMatch[oi] = 2;

            // Promote sites and apply "not real allele" adjustment.
            // hifiasm: for alt sites, set score=1 (promote).
            // For ref sites of this overlap, decrement occ_0 (ref support).
            for (uint32_t vi = 0; vi < nValid; vi++) {
                uint32_t si = validSites[vi];
                int8_t a = sites[si].allele[oi];
                if (a == 1) {
                    siteScore[si] = 1; // promote (already 1, but explicit)
                } else if (a == 0) {
                    // "Not real allele" adjustment: this read is ref at this
                    // site but trans overall, so its ref support is not real.
                    if (siteRefCov[vi] > 1) siteRefCov[vi]--;
                }
            }
        }

        // Step 3: Second pass — propagate.
        // hifiasm: for each overlap sorted by alt count, if still cis (is_match==1)
        // but has alt at a promoted site → mark as trans.
        for (uint32_t ri = 0; ri < numOv; ri++) {
            uint32_t oi = readOrder[ri];
            if (isMatch[oi] == 2) continue; // already trans
            if (readAlts[oi].nObs == 0) continue;

            uint32_t promotedAlts = 0;
            for (uint32_t vi = 0; vi < nValid; vi++) {
                uint32_t si = validSites[vi];
                if (sites[si].allele[oi] == 1 && siteScore[si] == 1) {
                    promotedAlts++;
                }
            }
            if (promotedAlts > 0) {
                isMatch[oi] = 2;
            }
        }

        // Step 4: Reset score for sites where all remaining cis overlaps
        // have ref (hifiasm: if overlap is cis, its alt sites get score=-1).
        // This prevents false promotions from propagating further.
        for (uint32_t vi = 0; vi < nValid; vi++) {
            uint32_t si = validSites[vi];
            bool anyCisAlt = false;
            for (uint32_t oi = 0; oi < numOv; oi++) {
                if (isMatch[oi] != 2 && sites[si].allele[oi] == 1) {
                    anyCisAlt = true;
                    break;
                }
            }
            if (!anyCisAlt) {
                siteScore[si] = -1;
            }
        }

        // Remaining reads with observations → cis.
        for (uint32_t oi = 0; oi < numOv; oi++) {
            if (isMatch[oi] == 0 && readAlts[oi].nObs > 0) {
                isMatch[oi] = 1;
            }
        }

        // Store per-read haplotype (backward compat).
        window.readHaplotypes.resize(numOv);
        uint32_t nCis = 0, nTrans = 0, nUnset = 0;
        for (uint32_t oi = 0; oi < numOv; oi++) {
            window.readHaplotypes[oi].orientedReadId = profiles[oi].oid;
            window.readHaplotypes[oi].hap = isMatch[oi];
            if (isMatch[oi] == 1) nCis++;
            else if (isMatch[oi] == 2) nTrans++;
            else nUnset++;
        }
        cout << "    transitive closure: cis=" << nCis << " trans=" << nTrans
             << " unclassified=" << nUnset << endl;

        // Build clusters: cis, trans, unclassified.
        vector<OrientedReadId> cisCluster, transCluster, uncCluster;
        for (uint32_t oi = 0; oi < numOv; oi++) {
            if (isMatch[oi] == 1)
                cisCluster.push_back(profiles[oi].oid);
            else if (isMatch[oi] == 2)
                transCluster.push_back(profiles[oi].oid);
            else
                uncCluster.push_back(profiles[oi].oid);
        }

        if (!cisCluster.empty())
            window.readClusters.push_back(std::move(cisCluster));
        if (!transCluster.empty())
            window.readClusters.push_back(std::move(transCluster));
        if (!uncCluster.empty())
            window.readClusters.push_back(std::move(uncCluster));

        cout << "    clusters: " << window.readClusters.size();
        for (size_t ci = 0; ci < window.readClusters.size(); ci++) {
            cout << " [" << ci << "]=" << window.readClusters[ci].size();
        }
        cout << endl;
    }

    return cleanHetSnps;
}
