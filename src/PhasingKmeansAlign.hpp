#ifndef DINARA_PHASING_KMEANS_ALIGN_HPP
#define DINARA_PHASING_KMEANS_ALIGN_HPP

/// @file PhasingKmeansAlign.hpp
/// @brief Noisy-region abPOA MSA for k-means phasing refinement.
///
/// See PhasingKmeansNoisyMsa.md for the full implementation plan.

#include "PhasingKmeansTypes.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace dinara {

class Assembler;

// ============================================================================
// Per-read data for a noisy region
// ============================================================================

struct KmNoisyReadInfo {
    int nReads = 0;
    std::vector<int>  overlapIndices;  ///< Indices into scratch.overlaps.
    std::vector<int>  haps;            ///< Hap label from initial k-means (1, 2, or 0).

    /// Per-read extracted sequences (0123 = ACGT encoding for abPOA).
    std::vector<std::vector<uint8_t>> seqs;
    std::vector<int> seqLens;

    /// Per-read: backbone base-position of the left flanking marker,
    /// relative to the backbone sequence start (0-based into backboneSeq).
    std::vector<int> bbLeftPos;

    /// Per-read: backbone base-position of the right flanking marker,
    /// relative to the backbone sequence start (0-based into backboneSeq).
    std::vector<int> bbRightPos;

    /// OrientedReadId values (for diagnostics).
    std::vector<uint64_t> orientedReadIdValues;

    /// Backbone sequence spanning the widest envelope (0123 encoding).
    std::vector<uint8_t> backboneSeq;
    int backboneLen = 0;

    /// Absolute backbone base-position of the left edge of the widest envelope.
    uint32_t backboneStartPos = 0;
};

// ============================================================================
// abPOA MSA result for one cluster
// ============================================================================

struct KmAbpoaMsaResult {
    int nReads = 0;                        ///< Number of reads in this cluster.
    std::vector<int> readIndices;          ///< Overlap indices into scratch.overlaps (-1 = backbone).
    std::vector<uint8_t> consensusSeq;     ///< Consensus sequence (0123 encoding).
    int consensusLen = 0;
    int msaLen = 0;                        ///< Number of columns in MSA.

    /// MSA rows for reads in this cluster.
    /// Each row has msaLen columns, 0=A, 1=C, 2=G, 3=T, 4=gap.
    std::vector<std::vector<uint8_t>> readMsaRows;

    /// MSA row for the consensus.
    std::vector<uint8_t> consensusMsaRow;

    /// MSA row for the backbone (read 0 in abPOA).
    /// Available for variant calling (consensus vs backbone).
    std::vector<uint8_t> backboneMsaRow;

    /// Absolute backbone base position of the MSA's left edge.
    /// Corresponds to the widest envelope's left marker midpoint.
    uint32_t backboneStartPos = 0;
};

// ============================================================================
// Options
// ============================================================================

struct KmNoisyMsaOptions {
    int maxNoisyRegLen  = 5000;
    int maxNoisyRegCov  = 200;
    int minHapReads     = 3;

    // abPOA scoring parameters (matching pgphase defaults).
    int match    = 2;
    int mismatch = 4;
    int gapOpen1 = 4;
    int gapExt1  = 2;
    int gapOpen2 = 24;
    int gapExt2  = 1;
    double minFreq = 0.3;
};

// ============================================================================
// Functions
// ============================================================================

/// Collect reads spanning [regStart, regEnd] on the backbone.
/// Extracts per-read sequences between flanking marker ordinals,
/// and the backbone sequence for the widest envelope.
KmNoisyReadInfo collectNoisyReadInfo(
    const Assembler& assembler,
    const KmScratchpad& scratch,
    uint32_t regStart, uint32_t regEnd);

/// Check if we have enough reads from both haplotypes.
bool hasBothHaplotypes(const KmNoisyReadInfo& info, int minHapReads);

/// Run abPOA on a subset of reads. Seeds the POA with the backbone,
/// then aligns each read to the subgraph between its flanking markers.
/// readIndices: which reads from info to include.
/// maxNCons: 1 for per-haplotype, 2 for de-novo clustering.
/// includeBackboneInReads: if true, include the backbone MSA row in
///   readMsaRows with readIndices=-1 (for hap1 where backbone is a member).
/// Returns number of consensus sequences produced (0, 1, or 2).
/// Up to 2 results are written to resultsOut.
int abpoaMsaRun(
    const KmNoisyReadInfo& info,
    const std::vector<int>& readIndices,
    int maxNCons,
    const KmNoisyMsaOptions& opts,
    bool includeBackboneInReads,
    std::array<KmAbpoaMsaResult, 2>& resultsOut);

/// Top-level: collect reads, run MSA (per-haplotype or combined fallback).
/// Returns number of consensus sequences produced (0, 1, or 2).
int collectNoisyRegMsa(
    const Assembler& assembler,
    const KmScratchpad& scratch,
    uint32_t regStart, uint32_t regEnd,
    const KmNoisyMsaOptions& opts,
    std::array<KmAbpoaMsaResult, 2>& results);

/// Outer loop: iterate over all noisy regions, run MSA, extract variants,
/// score reads, and merge into scratch. Port of pgphase collect_noisy_vars_step4.
void kmNoisyMsaStep4(
    const Assembler& assembler,
    KmScratchpad& scratch,
    const KmNoisyMsaOptions& opts);

} // namespace dinara

#endif
