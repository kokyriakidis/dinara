#ifndef DINARA_ANCHOR_WINDOWS_HPP
#define DINARA_ANCHOR_WINDOWS_HPP

#include "ReadId.hpp"
#include "Shasta2Anchors.hpp"
#include "cstdint.hpp"

#include <vector>

namespace dinara {

// A read's interval within an anchor window, expressed as journey positions.
struct AnchorWindowReadInterval {
    OrientedReadId orientedReadId;
    uint32_t begin = 0; // Inclusive position in the oriented read journey.
    uint32_t end = 0;   // Exclusive position in the oriented read journey.
    uint32_t touchedAnchorCount = 0; // Anchors shared with the backbone in this interval.
};

// An alternate path between two consecutive LIS (backbone-shared) anchors,
// formed by the non-shared anchors of a non-direct overlap read.
// The path goes: anchorIdA -> intermediateAnchorIds[0] -> ... -> anchorIdB
// where anchorIdA and anchorIdB are backbone anchors (LIS pillars).
struct AnchorWindowAlternatePath {
    Shasta2AnchorId anchorIdA;  // LIS pillar (backbone anchor) at start.
    Shasta2AnchorId anchorIdB;  // LIS pillar (backbone anchor) at end.
    std::vector<Shasta2AnchorId> intermediateAnchorIds; // Non-backbone anchors between pillars.
};

// An anchor window: a contiguous interval on a backbone read's journey,
// plus the intervals on all other reads that share anchors with the backbone.
struct AnchorWindow {
    uint32_t windowId = 0;
    OrientedReadId backboneOrientedReadId;
    uint32_t backboneBegin = 0; // Inclusive journey position on backbone.
    uint32_t backboneEnd = 0;   // Exclusive journey position on backbone.
    uint32_t claimedAnchorCount = 0;
    std::vector<AnchorWindowReadInterval> readIntervals;

    // Alternate paths from non-direct overlap reads.
    // These form parallel chains between backbone anchors at het sites.
    std::vector<AnchorWindowAlternatePath> alternatePaths;

    // Set by cigarDetectSnpsInWindow: number of clean het SNPs found.
    // 0 means the window is homozygous — only backbone anchors should be kept.
    uint32_t cleanHetSnpCount = 0;

    // Passing het SNP positions (in oriented backbone coordinates).
    struct HetSnp {
        uint32_t bbPos;       // backbone position
        uint8_t altBase;      // alt allele (0=A, 1=C, 2=G, 3=T)
        uint32_t altCov;      // alt read count
        uint32_t refCov;      // ref read count
        uint32_t spanning;    // total reads with bases at this position
        std::vector<OrientedReadId> altReads;  // reads carrying alt
        std::vector<OrientedReadId> refReads;  // reads carrying ref
    };
    std::vector<HetSnp> hetSnps;

    // Per-read haplotype assignment from k-means phasing.
    // hap: 0 = unassigned, 1 = haplotype 1 (cis with backbone), 2 = haplotype 2 (trans).
    struct ReadHaplotype {
        OrientedReadId orientedReadId;
        int hap = 0;
    };
    std::vector<ReadHaplotype> readHaplotypes;
};

} // namespace dinara

#endif // DINARA_ANCHOR_WINDOWS_HPP
