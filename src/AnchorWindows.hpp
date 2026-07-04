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

// A single anchor a member read shares with the backbone (a "pin").
// Stored as journey positions, which is what window construction already has;
// consumers derive marker ordinals / base positions via getOrdinal as needed.
// Pins for one read are stored in backbone order (ascending backboneJourneyPos).
struct AnchorWindowSharedPin {
    uint32_t readJourneyPos = 0;     // Position in the member's oriented journey.
    uint32_t backboneJourneyPos = 0; // Position in the backbone's journey.
};

// A read's interval within an anchor window, expressed as journey positions.
struct AnchorWindowReadInterval {
    OrientedReadId orientedReadId;
    uint32_t begin = 0; // Inclusive position in the oriented read journey.
    uint32_t end = 0;   // Exclusive position in the oriented read journey.
    uint32_t touchedAnchorCount = 0; // Anchors shared with the backbone in this interval.

    // Backbone anchors this read shares with the backbone (its pins), in
    // backbone order. Populated during window construction so downstream passes
    // (per-segment abPOA MSA, het-site detection) reuse them instead of
    // re-intersecting journeys. Consecutive pins delimit inter-anchor segments.
    std::vector<AnchorWindowSharedPin> sharedPins;

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

    // Base span of the backbone (distance in bases from first to last
    // backbone anchor on the backbone read). Set during window construction.
    uint64_t baseSpan = 0;

    // Endpoint windows: the windows that the backbone read transitions
    // to/from. These are the "real" connections. All other inter-window
    // connections are internal (from non-backbone reads).
    // Set after transitionReads are populated.
    // noWindow means the backbone starts/ends here.
    uint32_t backbonePreviousWindow = AnchorWindowReadInterval::noWindow;
    uint32_t backboneNextWindow = AnchorWindowReadInterval::noWindow;

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

    // Staged het-anchor descriptors from per-window abPOA SNP detection.
    // Each allele of a clean SNP bubble becomes one k=2 anchor: its members are
    // the reads on that allele, each pinned at rawPosition (the read's base
    // position at the bubble's common predecessor). shasta2 re-derives the
    // 2-base k-mer [predBase, alleleBase] from rawPosition, so all members of
    // one allele share the same 2 bases and the two alleles differ in the
    // second base (forming the bubble). Anchor IDs are assigned in a later
    // serial pass that appends them to the primary anchor store.
    struct HetAnchorMember {
        OrientedReadId orientedReadId;
        uint32_t rawPosition;   // read base position at commonPred (predBase)
    };
    struct HetAnchor {
        uint32_t backboneOffset = 0;   // backbone base offset of the SNP column
        uint8_t predBase = 0;          // common predecessor base (0-3)
        uint8_t alleleBase = 0;        // this allele's base (0-3)
        bool isRef = false;            // true = backbone/reference allele
        std::vector<HetAnchorMember> members;
        // Assigned in the post-window append pass (invalid until then).
        Shasta2AnchorId anchorId = std::numeric_limits<Shasta2AnchorId>::max();
    };
    // Grouped per SNP: [0] = ref allele anchor, [1..] = alt allele anchor(s).
    // The predecessor/successor backbone anchors that flank each bubble are the
    // window backbone anchors nearest to backboneOffset; resolved when the
    // anchor graph bubble edges are built.
    struct HetBubble {
        uint32_t backboneOffset = 0;
        std::vector<HetAnchor> alleles;   // [0]=ref, rest=alts
        // Each bubble is bracketed by TWO all-reads k=2 hom anchors so that a
        // backbone anchor only ever connects to a hom (which the backbone read
        // is a member of), never directly to an allele arm. A minority-allele
        // arm's reads are NOT members of the flanking k=50 backbone anchor (they
        // diverge at the 50-mer locus), so a direct backbone->arm edge has an
        // empty read intersection and cannot form a valid AnchorPair. The homs
        // solve this: every read entering/leaving the site passes through them.
        //
        //   leadHom is at predPrev ([predPrevBase, predBase]); all reads
        //   entering the site share it (predPrev->commonPred is linear).
        //   hom (trailing) is at commonSucc ([succBase, nextBase]); all reads
        //   leaving the site share it (commonSucc->succNext is linear).
        //
        // A single interval's chain is thus:
        //   bbA_i -> leadHom_0 -> {alleles_0} -> hom_0
        //         -> leadHom_1 -> {alleles_1} -> hom_1
        //         -> ... -> leadHom_n -> {alleles_n} -> hom_n -> bbA_{i+1}
        // Both backbone-touching edges (bbA_i -> leadHom_0 and hom_n ->
        // bbA_{i+1}) are hom<->backbone, so their intersection is non-empty. All
        // interior edges touch a het/hom anchor and are read-consistent by
        // construction. Every bubble now carries BOTH homs (no last-bubble
        // special case). Empty members => that hom unavailable (bubble dropped).
        uint32_t predBackboneOffset = 0;       // backbone offset of leadHom
        HetAnchor leadHom;
        uint32_t succBackboneOffset = 0;       // backbone offset of trailing hom
        HetAnchor hom;

        // Transient planning field, set by the intra-window edge planner (a
        // pre-pass over the windows) and consumed by the append + staging passes
        // so they agree on exactly which anchors get created and wired (no orphan
        // anchors, no unwired bubbles).
        //   plannedInterval : index i of the backbone interval [bbA_i,
        //                     bbA_{i+1}) that strictly contains this bubble's
        //                     full flank span [predBackboneOffset,
        //                     succBackboneOffset], or -1 if none (bubble
        //                     dropped). Since every bubble is now self-bracketed
        //                     by its own lead/trailing hom, there is no
        //                     last-in-interval special case.
        int32_t plannedInterval = -1;

        // Coincident-hom merge (Pass 1.5). When this bubble's leading hom sits
        // on the SAME backbone column (POA node) as the PRECEDING bubble's
        // trailing hom -- i.e. predBackboneOffset == that bubble's
        // succBackboneOffset, which happens when two accepted SNPs are exactly
        // 3 bp apart -- the two homs are one node with identical members. This
        // field holds the hetBubbles index of that preceding bubble, and its
        // trailing-hom anchor is reused as this bubble's leading hom: the append
        // pass does NOT allocate a separate anchor for this leadHom and the
        // staging pass omits the redundant leadHom step, so the chain is
        // ...arms_prev -> sharedHom -> arms_this... (-1 = not shared).
        int64_t sharedLeadFromBubble = -1;
    };
    std::vector<HetBubble> hetBubbles;

    // Window-local (intra-window) anchor-graph edges, staged on the window so
    // all edges that are fully determined by a single window live in the window
    // structure. The anchor graph constructor replays these directly. Two
    // kinds share this list:
    //   - backbone chain edges (consecutive backbone anchors, + RC mirror):
    //     isHet=false; the constructor builds them via the normal marker-based
    //     Shasta2AnchorPair path (offset is recomputed, the staged offset is
    //     ignored).
    //   - het bubble edges (pred->allele, allele->succ, + RC mirror): isHet=true;
    //     the constructor builds them via the direct read-intersection path
    //     using the staged nominal offset (het anchors are k=2 and have no real
    //     marker ordinal).
    // Inter-window edges are NOT staged here: they depend on multiple windows
    // and are discovered in a separate global pass.
    struct IntraWindowEdge {
        Shasta2AnchorId anchorIdA = std::numeric_limits<Shasta2AnchorId>::max();
        Shasta2AnchorId anchorIdB = std::numeric_limits<Shasta2AnchorId>::max();
        uint32_t offset = 0;   // nominal base offset (used only for het edges)
        bool isHet = false;
    };
    std::vector<IntraWindowEdge> intraWindowEdges;

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

    // Inter-window edges incident to this window.
    // Populated after inter-window edge creation.
    struct InterWindowEdge {
        uint32_t otherWindow;       // The neighbor window (raw ID, may be RC mirror).
        Shasta2AnchorId anchorIdA;  // Last anchor in source window.
        Shasta2AnchorId anchorIdB;  // First anchor in destination window.
        uint64_t readCount;         // Shared read count.
    };
    std::vector<InterWindowEdge> outEdges;  // This window is the source.
    std::vector<InterWindowEdge> inEdges;   // This window is the destination.
};

} // namespace dinara

#endif // DINARA_ANCHOR_WINDOWS_HPP
