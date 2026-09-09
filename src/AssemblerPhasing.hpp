#ifndef DINARA_ASSEMBLER_PHASING_HPP
#define DINARA_ASSEMBLER_PHASING_HPP

#include "ReadId.hpp"
#include "ProjectedAlignment.hpp"
#include "MemoryMappedVector.hpp"
#include <span>
#include <vector>
#include <cstdint>

namespace dinara {

class Assembler;

// Configuration for the phasing algorithm
struct PhasingConfig {
    double heterozygosityRate = 0.001;
    int windowSize = 512;
    int minEvidence = 2;       // Minimum reads supporting an allele to consider it a Het SNP
    int maxDrift = 10;
    
    // Hifiasm Parity Parameters
    int s_hap_cov = 3;         // Minimal support for a valid haplotype (default ~3)
    int infor_cov = 1;         // Minimal support for an informative allele (default 1)
    double st_rate = 0.05;     // Systematic error rate threshold
    int st_max = 2;            // Systematic error absolute count threshold (rarely used)
    uint64_t hom_cov = 0;      // Homozygous peak coverage (Peak Coverage)
    bool is_hpc = false;       // Enable homopolymer checks (usually false for HiFi, but Hifiasm has logic)
    bool is_ont = true;        // ONT mode flag (affects thresholds)
};

// Represents a single piece of evidence (a base at a specific site)
// Similar to Hifiasm's haplotype_evdience
struct HaplotypeEvidence {
    uint32_t site;          // 0-based position on the TARGET read
    uint32_t overlapId;     // Index of the overlap in the local list
    uint8_t base;           // 0=A, 1=C, 2=G, 3=T, 4=Gap
    uint8_t type;           // 0=Match, 1=SNP/Mismatch
    uint8_t isSolid;        // Flag for high-confidence regions
    bool isSameStrand;      // Parity: Tracks strand relative to target (0=Same)
};

// Represents an overlap tailored for phasing, derived from alignment data/CIGAR
struct PhasingOverlap {
    uint32_t alnIdx;            // Index in global alignment table
    uint32_t targetReadId;      // ID of the target/reference read
    uint32_t queryReadId;       // ID of the query read
    uint32_t targetStart;       // Start position on target
    uint32_t targetEnd;         // End position on target
    uint32_t queryStart;        // Start position on query (strand aware)
    uint32_t queryEnd;          // End position on query
    bool queryStrand;           // 0=forward, 1=reverse
    std::vector<uint32_t> cigar; // Encoded CIGAR string
};

class AssemblerPhasing {
public:
    // Main entry point: Filters a list of overlaps for a specific target read
    // by removing those that appear to belong to a different haplotype.
    // Returns indices of overlaps to KEEP.
    static std::vector<uint32_t> filterOverlapsByPhasing(
        const Assembler& assembler,
        ReadId targetReadId,
        const std::vector<PhasingOverlap>& overlaps,
        const PhasingConfig& config = PhasingConfig()
    );

    struct SnpStats {
        uint32_t site;
        uint32_t occ_0; // Ref count
        uint32_t occ_1; // Alt count
        uint32_t overlap_num;
        uint32_t homopolymer_num;     // Parity: Count of evidence inside homopolymers
        uint32_t non_homopolymer_num; // Parity: Count of evidence outside homopolymers
    };
    
    // Returns phasing score per overlap (Hifiasm parity)
    // score > 0 = Alt-supporting (TRANS)
    // score <= 0 = Ref-supporting (CIS)
    struct OverlapPhasingResult {
        int32_t score;           // Positive = TRANS, non-positive = CIS
        bool hasInformativeSite; // True if overlap covers valid SNP site
    };
    static std::vector<OverlapPhasingResult> getOverlapPhasingScores(
        const Assembler& assembler,
        ReadId targetReadId,
        const std::vector<PhasingOverlap>& overlaps,
        const PhasingConfig& config = PhasingConfig()
    );

private:
    // internal helpers
    static void collectHaplotypeEvidence(
        const Assembler& assembler,
        const PhasingOverlap& overlap,
        std::vector<HaplotypeEvidence>& evidenceOut,
        const PhasingConfig& config
    );
    


    // Updated Naive refinement (Stateful)
    static std::vector<uint32_t> refineOverlapsNaive(
        const std::vector<uint32_t>& keptIndices,
        const std::vector<PhasingOverlap>& overlaps,
        const std::vector<HaplotypeEvidence>& evidence,
        std::vector<SnpStats>& snpStats, // Mutable stats
        const PhasingConfig& config
    );

    // SV Phasing (Parity with rphase_lidel)
    static std::vector<uint32_t> filterOverlapsBySV(
        const std::vector<uint32_t>& keptIndices,
        const std::vector<PhasingOverlap>& overlaps,
        const std::vector<HaplotypeEvidence>& evidence,
        std::vector<SnpStats>& snpStats,
        const Assembler& assembler,
        const PhasingConfig& config
    );

    // Returns stats populated during DP preparation.
    static std::vector<uint32_t> generatePhasingDP(
        ReadId targetReadId,
        const std::vector<PhasingOverlap>& overlaps,
        const std::vector<HaplotypeEvidence>& evidence,
        std::vector<SnpStats>& outStats,
        const PhasingConfig& config
    );

};

} // namespace dinara

#endif // DINARA_ASSEMBLER_PHASING_HPP
