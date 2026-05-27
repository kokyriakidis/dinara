#ifndef DINARA_ANCHOR_WINDOWS_HPP
#define DINARA_ANCHOR_WINDOWS_HPP

#include "ReadId.hpp"
#include "Shasta2Anchors.hpp"
#include "cstdint.hpp"

#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace dinara {

// A read's interval within an anchor window, expressed as journey positions.
struct AnchorWindowReadInterval {
    OrientedReadId orientedReadId;
    uint32_t begin = 0; // Inclusive position in the oriented read journey.
    uint32_t end = 0;   // Exclusive position in the oriented read journey.
    uint32_t touchedAnchorCount = 0; // Anchors shared with the backbone in this interval.

    // Window transitions: which window this read came from and goes to.
    // noWindow means the read starts/ends here (no adjacent window).
    static constexpr uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    uint32_t previousWindow = noWindow; // Window this read was in before this one.
    uint32_t nextWindow = noWindow;     // Window this read goes to after this one.
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

    // Filtered backbone positions: the longest subsequence of journey
    // positions in [backboneBegin, backboneEnd) where every consecutive
    // pair has sufficient common read support. Intra-window edges should
    // be created between consecutive entries in this vector.
    // If empty, use all positions in [backboneBegin, backboneEnd).
    std::vector<uint32_t> filteredBackbonePositions;

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
        uint8_t refBase;      // ref allele (0=A, 1=C, 2=G, 3=T) — backbone base
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

    // Per-read window transition map: groups reads by their (previousWindow, nextWindow)
    // transition pattern. Used for detangling decisions.
    // Key: (previousWindow, nextWindow). noWindow means the read starts/ends here.
    // Value: oriented read IDs following that transition.
    static constexpr uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    std::map<std::pair<uint32_t, uint32_t>, std::vector<OrientedReadId>> transitionReads;

    // Backbone read's transition: which windows the backbone connects from/to.
    uint32_t backbonePreviousWindow = noWindow;
    uint32_t backboneNextWindow = noWindow;

    // Read clusters from phasing.
    // Each cluster is a set of reads that phase together.
    // Cluster 0 contains the backbone read.
    std::vector<std::vector<OrientedReadId>> readClusters;

    // Compute per-cluster transition summaries from readClusters and transitionReads.
    // Returns one map per cluster: (previousWindow, nextWindow) -> read count.
    // Must be called after both readClusters and transitionReads are populated.
    std::vector<std::map<std::pair<uint32_t, uint32_t>, uint64_t>>
    computeClusterTransitions() const {
        // Build reverse lookup: orientedReadId -> (previousWindow, nextWindow).
        std::map<uint32_t, std::pair<uint32_t, uint32_t>> readTransition;
        for(const auto& [key, reads] : transitionReads) {
            for(const auto& oid : reads) {
                readTransition[oid.getValue()] = key;
            }
        }

        std::vector<std::map<std::pair<uint32_t, uint32_t>, uint64_t>> result;
        result.resize(readClusters.size());
        for(uint64_t ci = 0; ci < readClusters.size(); ci++) {
            for(const auto& oid : readClusters[ci]) {
                auto it = readTransition.find(oid.getValue());
                if(it != readTransition.end()) {
                    result[ci][it->second]++;
                }
            }
        }
        return result;
    }

    // True if the cluster at the same index contains unclassified reads
    // (reads that don't span any het sites and belong to both haplotypes).
    std::vector<bool> clusterIsUnclassified;
};

} // namespace dinara

#endif // DINARA_ANCHOR_WINDOWS_HPP
