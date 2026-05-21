/// @file PhasingKmeansNoisyMsa.cpp
/// @brief Noisy-region MSA variant calling and re-phasing (pgphase step 4).
///
/// See PhasingKmeansNoisyMsa.md for the full implementation plan.
///
/// For each noisy region:
/// 1. Run abPOA MSA via collectNoisyRegMsa (per-haplotype or combined).
/// 2. Walk the MSA matrix to extract variants (consensus vs backbone).
/// 3. Score all reads at the new variant positions.
/// 4. Merge new variants into the candidate table.
/// 5. Re-run k-means with the expanded candidate set.

#include "PhasingKmeansAlign.hpp"
#include "PhasingKmeansTypes.hpp"
#include "Assembler.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

using namespace dinara;
using namespace std;

// ============================================================================
// Constants
// ============================================================================

/// abPOA gap value in MSA matrix (0=A, 1=C, 2=G, 3=T, 4=gap).
static constexpr uint8_t kGap = 4;

// ============================================================================
// Variant extraction from MSA matrix
// ============================================================================

/// A variant discovered from the MSA consensus vs backbone.
struct NoisyMsaVariant {
    uint32_t backbonePos;  ///< Position on the backbone (absolute).
    KmVarType type;
    uint8_t refBase;       ///< For SNPs: reference base (0-3).
    uint8_t altBase;       ///< For SNPs: alt base (0-3).
    string refSeq;         ///< For indels: ref sequence (char).
    string altSeq;         ///< For indels: alt sequence (char).
    uint32_t msaColStart;  ///< MSA column where this variant starts.
    uint32_t msaColEnd;    ///< MSA column past the end of this variant.
};

/// Build a KmVarKey from a NoisyMsaVariant.
static KmVarKey varToKey(const NoisyMsaVariant& v) {
    KmVarKey key;
    key.pos = v.backbonePos;
    key.type = v.type;
    key.altBase = (v.type == KmVarType::Snp) ? v.altBase : 0;
    if(v.type == KmVarType::Snp) {
        key.refLen = 1; key.altLen = 1;
    } else if(v.type == KmVarType::Insertion) {
        key.refLen = 0; key.altLen = uint16_t(v.altSeq.size());
    } else if(v.type == KmVarType::Deletion) {
        key.refLen = uint16_t(v.refSeq.size()); key.altLen = 0;
    }
    return key;
}

/// Walk the MSA to find positions where the consensus row differs from the
/// backbone row.  Both rows are in 0123/gap encoding from abPOA.
///
/// backboneStartPos: absolute backbone base position of the MSA's left edge.
static vector<NoisyMsaVariant> extractVariantsFromMsa(
    const KmAbpoaMsaResult& result,
    uint32_t backboneStartPos)
{
    vector<NoisyMsaVariant> vars;
    if(result.msaLen == 0) return vars;

    const auto& bbRow = result.backboneMsaRow;
    const auto& consRow = result.consensusMsaRow;
    if(int(bbRow.size()) != result.msaLen) return vars;
    if(int(consRow.size()) != result.msaLen) return vars;

    static const char bases[] = "ACGT";

    // Track backbone position (skip gap columns in backbone).
    uint32_t bbPos = backboneStartPos;

    uint32_t col = 0;
    while(col < uint32_t(result.msaLen)) {
        uint8_t bbVal = bbRow[col];
        uint8_t consVal = consRow[col];

        // Both same → no variant.
        if(bbVal == consVal) {
            if(bbVal != kGap) bbPos++;
            col++;
            continue;
        }

        // SNP: backbone has a base, consensus has a different base.
        if(bbVal < 4 && consVal < 4) {
            NoisyMsaVariant v;
            v.backbonePos = bbPos;
            v.type = KmVarType::Snp;
            v.refBase = bbVal;
            v.altBase = consVal;
            v.msaColStart = col;
            v.msaColEnd = col + 1;
            vars.push_back(v);
            bbPos++;
            col++;
            continue;
        }

        // Insertion: backbone is gap, consensus has a base.
        if(bbVal == kGap && consVal < 4) {
            string insSeq;
            uint32_t insColStart = col;
            while(col < uint32_t(result.msaLen) &&
                  bbRow[col] == kGap && consRow[col] < 4) {
                insSeq += bases[consRow[col]];
                col++;
            }
            if(!insSeq.empty()) {
                NoisyMsaVariant v;
                v.backbonePos = bbPos; // Position before the insertion.
                v.type = KmVarType::Insertion;
                v.altSeq = move(insSeq);
                v.msaColStart = insColStart;
                v.msaColEnd = col;
                vars.push_back(v);
            }
            continue; // Don't increment bbPos — backbone had gaps.
        }

        // Deletion: backbone has a base, consensus is gap.
        if(bbVal < 4 && consVal == kGap) {
            string delSeq;
            uint32_t delStartBbPos = bbPos;
            uint32_t delColStart = col;
            while(col < uint32_t(result.msaLen) &&
                  bbRow[col] < 4 && consRow[col] == kGap) {
                delSeq += bases[bbRow[col]];
                bbPos++;
                col++;
            }
            if(!delSeq.empty()) {
                NoisyMsaVariant v;
                v.backbonePos = delStartBbPos;
                v.type = KmVarType::Deletion;
                v.refSeq = move(delSeq);
                v.msaColStart = delColStart;
                v.msaColEnd = col;
                vars.push_back(v);
            }
            continue;
        }

        // Other cases (both gap, etc.) — skip.
        if(bbVal != kGap) bbPos++;
        col++;
    }

    return vars;
}

// ============================================================================
// Score reads at variant positions
// ============================================================================

/// For each variant, score each read in the MSA: does it match
/// the ref allele (0), the alt allele (1), or neither (-1)?
///
/// Returns scores[varIdx][readIdx] where readIdx indexes into
/// result.readMsaRows / result.readIndices.
static vector<vector<int8_t>> scoreReadsAtVariants(
    const KmAbpoaMsaResult& result,
    const vector<NoisyMsaVariant>& vars)
{
    const int nReads = result.nReads;
    vector<vector<int8_t>> scores(vars.size());
    if(nReads == 0) return scores;

    for(size_t vi = 0; vi < vars.size(); vi++) {
        scores[vi].assign(nReads, -1);
    }

    // For SNPs, we can score directly at the variant's MSA column.
    // For indels, we check the full span of MSA columns.
    for(size_t vi = 0; vi < vars.size(); vi++) {
        const auto& v = vars[vi];

        if(v.type == KmVarType::Snp) {
            uint32_t col = v.msaColStart;
            if(col >= uint32_t(result.msaLen)) continue;

            for(int r = 0; r < nReads; r++) {
                uint8_t readVal = result.readMsaRows[r][col];
                if(readVal == v.refBase) scores[vi][r] = 0;
                else if(readVal == v.altBase) scores[vi][r] = 1;
                // else: gap or other base → stays -1
            }
        }
        else if(v.type == KmVarType::Insertion) {
            // Insertion: backbone is gap, alt is bases.
            // Read has insertion → bases in the gap columns.
            // Read has ref → gaps in the gap columns.
            for(int r = 0; r < nReads; r++) {
                int nBases = 0, nGaps = 0;
                for(uint32_t col = v.msaColStart; col < v.msaColEnd; col++) {
                    uint8_t readVal = result.readMsaRows[r][col];
                    if(readVal < 4) nBases++;
                    else if(readVal == kGap) nGaps++;
                }
                int span = int(v.msaColEnd - v.msaColStart);
                if(nGaps == span) scores[vi][r] = 0;       // all gaps → ref
                else if(nBases == span) scores[vi][r] = 1;  // all bases → alt
                // mixed → stays -1
            }
        }
        else if(v.type == KmVarType::Deletion) {
            // Deletion: backbone has bases, consensus is gap.
            // Read has deletion → gaps in these columns.
            // Read has ref → bases in these columns.
            for(int r = 0; r < nReads; r++) {
                int nBases = 0, nGaps = 0;
                for(uint32_t col = v.msaColStart; col < v.msaColEnd; col++) {
                    uint8_t readVal = result.readMsaRows[r][col];
                    if(readVal < 4) nBases++;
                    else if(readVal == kGap) nGaps++;
                }
                int span = int(v.msaColEnd - v.msaColStart);
                if(nBases == span) scores[vi][r] = 0;      // all bases → ref
                else if(nGaps == span) scores[vi][r] = 1;   // all gaps → alt (deletion)
                // mixed → stays -1
            }
        }
    }

    return scores;
}

// ============================================================================
// Convert NoisyMsaVariant to KmCandidate
// ============================================================================

static KmCandidate variantToCandidate(
    const NoisyMsaVariant& v,
    const vector<int8_t>& readScores,
    KmVariantCategory category)
{
    KmCandidate cand;
    cand.key = varToKey(v);

    // Compute allele counts from all read scores including backbone.
    // The backbone is a real hap1 read and its allele should be counted.
    for(size_t i = 0; i < readScores.size(); i++) {
        int8_t s = readScores[i];
        if(s == 0) { cand.refCov++; cand.totalCov++; }
        else if(s == 1) { cand.altCov++; cand.totalCov++; }
    }

    cand.category = category;
    cand.categoryFlag = kmCategoryToFlag(category);
    return cand;
}

// ============================================================================
// Merge new variants into scratch
// ============================================================================

/// Merge variants from one haplotype's MSA into the scratch candidate table.
/// readIndices maps MSA rows to overlap indices (-1 = backbone, skip).
static void mergeNoisyVariants(
    KmScratchpad& scratch,
    const vector<NoisyMsaVariant>& vars,
    const vector<vector<int8_t>>& readScores,
    const vector<int>& readIndices,
    KmVariantCategory category)
{
    for(size_t vi = 0; vi < vars.size(); vi++) {
        const auto& v = vars[vi];
        const auto& scores = readScores[vi];

        // Skip variants with no alt support.
        int altCount = 0;
        for(size_t i = 0; i < scores.size(); i++) {
            if(scores[i] == 1) altCount++;
        }
        if(altCount == 0) continue;

        // Check if this variant already exists in scratch.candidates.
        KmVarKey key = varToKey(v);
        int candIdx = -1;
        for(int ci = 0; ci < int(scratch.candidates.size()); ci++) {
            if(scratch.candidates[ci].key == key) {
                candIdx = ci;
                break;
            }
        }

        if(candIdx < 0) {
            // New variant — add to candidates.
            KmCandidate cand = variantToCandidate(v, scores, category);
            candIdx = int(scratch.candidates.size());
            scratch.candidates.push_back(cand);
        }

        // Update overlap profiles for reads that were in the MSA.
        // This runs for both new and existing variants, so reads from
        // a second haplotype's MSA still get profiled.
        for(size_t ri = 0; ri < scores.size(); ri++) {
            int oi = readIndices[ri];
            if(oi < 0) continue; // backbone sentinel
            if(oi >= int(scratch.overlapProfiles.size())) continue;

            auto& prof = scratch.overlapProfiles[oi];
            if(prof.startVarIdx < 0) {
                prof.startVarIdx = candIdx;
                prof.endVarIdx = candIdx;
            }
            int off = candIdx - prof.startVarIdx;
            if(off < 0) continue;
            while(int(prof.alleles.size()) <= off) {
                prof.alleles.push_back(-1);
            }
            prof.alleles[off] = scores[ri];
            if(candIdx > prof.endVarIdx) {
                prof.endVarIdx = candIdx;
            }
        }
    }
}

// ============================================================================
// Cross-haplotype scoring
// ============================================================================

/// Build a map from backbone base position to MSA column index.
/// For insertion columns (backbone is gap), no entry is added.
/// Returns a map where bbPosToCol[bbPos - backboneStartPos] = MSA column.
static vector<int> buildBbPosToMsaCol(
    const KmAbpoaMsaResult& result)
{
    // Map backbone-relative position → MSA column.
    // Count backbone bases to determine the size.
    int nBbBases = 0;
    for(int col = 0; col < result.msaLen; col++) {
        if(result.backboneMsaRow[col] != kGap) nBbBases++;
    }

    vector<int> posToCol(nBbBases, -1);
    int bbIdx = 0;
    for(int col = 0; col < result.msaLen; col++) {
        if(result.backboneMsaRow[col] != kGap) {
            if(bbIdx < nBbBases) posToCol[bbIdx] = col;
            bbIdx++;
        }
    }
    return posToCol;
}

/// Score reads from one MSA at variant positions discovered from another MSA.
/// Uses backbone position to find the corresponding MSA column.
static void crossScoreReads(
    const KmAbpoaMsaResult& result,
    const vector<NoisyMsaVariant>& vars,
    vector<vector<int8_t>>& scores)
{
    if(result.nReads == 0 || result.msaLen == 0) return;

    vector<int> posToCol = buildBbPosToMsaCol(result);

    for(size_t vi = 0; vi < vars.size(); vi++) {
        const auto& v = vars[vi];
        int relPos = int(v.backbonePos - result.backboneStartPos);
        if(relPos < 0 || relPos >= int(posToCol.size())) {
            // Out of range — all reads get unknown score.
            for(int r = 0; r < result.nReads; r++)
                scores[vi].push_back(-1);
            continue;
        }

        if(v.type == KmVarType::Snp) {
            int col = posToCol[relPos];
            if(col < 0 || col >= result.msaLen) {
                for(int r = 0; r < result.nReads; r++)
                    scores[vi].push_back(-1);
                continue;
            }

            for(int r = 0; r < result.nReads; r++) {
                uint8_t readVal = result.readMsaRows[r][col];
                int8_t score = -1;
                if(readVal == v.refBase) score = 0;
                else if(readVal == v.altBase) score = 1;
                scores[vi].push_back(score);
            }
        }
        else if(v.type == KmVarType::Insertion) {
            // For insertions, find the insertion columns after the backbone position.
            // In this MSA, the insertion may or may not exist at the same columns.
            // Look for gap columns in the backbone row after posToCol[relPos].
            int anchorCol = posToCol[relPos];
            if(anchorCol < 0) {
                for(int r = 0; r < result.nReads; r++)
                    scores[vi].push_back(-1);
                continue;
            }
            // Count consecutive gap columns in backbone after the anchor.
            int insStart = anchorCol + 1;
            int insEnd = insStart;
            while(insEnd < result.msaLen && result.backboneMsaRow[insEnd] == kGap)
                insEnd++;
            int insSpan = insEnd - insStart;

            for(int r = 0; r < result.nReads; r++) {
                if(insSpan == 0) {
                    // No insertion columns in this MSA → read has ref.
                    scores[vi].push_back(0);
                } else {
                    int nBases = 0, nGaps = 0;
                    for(int col = insStart; col < insEnd; col++) {
                        uint8_t rv = result.readMsaRows[r][col];
                        if(rv < 4) nBases++;
                        else if(rv == kGap) nGaps++;
                    }
                    if(nGaps == insSpan) scores[vi].push_back(0);
                    else if(nBases == insSpan) scores[vi].push_back(1);
                    else scores[vi].push_back(-1);
                }
            }
        }
        else if(v.type == KmVarType::Deletion) {
            // For deletions, check the backbone columns starting at relPos.
            int delLen = int(v.refSeq.size());
            int startCol = posToCol[relPos];
            if(startCol < 0) {
                for(int r = 0; r < result.nReads; r++)
                    scores[vi].push_back(-1);
                continue;
            }
            // Find the next delLen backbone-base columns.
            vector<int> delCols;
            delCols.reserve(delLen);
            for(int p = relPos; p < relPos + delLen && p < int(posToCol.size()); p++) {
                if(posToCol[p] >= 0) delCols.push_back(posToCol[p]);
            }

            for(int r = 0; r < result.nReads; r++) {
                if(int(delCols.size()) != delLen) {
                    scores[vi].push_back(-1);
                    continue;
                }
                int nBases = 0, nGaps = 0;
                for(int col : delCols) {
                    uint8_t rv = result.readMsaRows[r][col];
                    if(rv < 4) nBases++;
                    else if(rv == kGap) nGaps++;
                }
                if(nBases == delLen) scores[vi].push_back(0);
                else if(nGaps == delLen) scores[vi].push_back(1);
                else scores[vi].push_back(-1);
            }
        }
    }
}

// ============================================================================
// Process per-haplotype MSA results
// ============================================================================

/// Process two per-haplotype MSA results: extract variants from both,
/// cross-score reads from each haplotype at all variant positions,
/// then merge into scratch.
static int processTwoHapResults(
    KmScratchpad& scratch,
    const array<KmAbpoaMsaResult, 2>& results)
{
    // Extract variants from both haplotypes.
    vector<NoisyMsaVariant> vars0 = extractVariantsFromMsa(
        results[0], results[0].backboneStartPos);
    vector<NoisyMsaVariant> vars1 = extractVariantsFromMsa(
        results[1], results[1].backboneStartPos);

    // Merge variant lists (union, deduplicated by position+type+allele).
    vector<NoisyMsaVariant> allVars = vars0;
    for(const auto& v : vars1) {
        bool dup = false;
        for(const auto& ev : allVars) {
            if(ev.backbonePos == v.backbonePos && ev.type == v.type &&
               ev.altBase == v.altBase &&
               ev.refSeq == v.refSeq && ev.altSeq == v.altSeq) {
                dup = true;
                break;
            }
        }
        if(!dup) allVars.push_back(v);
    }
    if(allVars.empty()) return 0;

    // Determine category: variants found in both haplotypes are hom,
    // variants found in only one are het.
    vector<KmVariantCategory> categories(allVars.size());
    for(size_t vi = 0; vi < allVars.size(); vi++) {
        KmVarKey key = varToKey(allVars[vi]);
        bool inHap1 = false, inHap2 = false;
        for(const auto& v0 : vars0) {
            if(varToKey(v0) == key) { inHap1 = true; break; }
        }
        for(const auto& v1 : vars1) {
            if(varToKey(v1) == key) { inHap2 = true; break; }
        }
        categories[vi] = (inHap1 && inHap2)
            ? KmVariantCategory::NoisyCandHom
            : KmVariantCategory::NoisyCandHet;
    }

    // Score reads from both MSAs at all variant positions using
    // backbone-position-to-MSA-column mapping (cross-scoring).
    vector<vector<int8_t>> scores0(allVars.size());
    crossScoreReads(results[0], allVars, scores0);

    vector<vector<int8_t>> scores1(allVars.size());
    crossScoreReads(results[1], allVars, scores1);

    // Merge scores and readIndices from both haplotypes.
    for(size_t vi = 0; vi < allVars.size(); vi++) {
        const auto& v = allVars[vi];

        // Combined scores: hap1 reads + hap2 reads.
        vector<int8_t> combinedScores;
        vector<int> combinedReadIndices;
        for(int r = 0; r < int(scores0[vi].size()); r++) {
            combinedScores.push_back(scores0[vi][r]);
            combinedReadIndices.push_back(
                r < int(results[0].readIndices.size())
                    ? results[0].readIndices[r] : -1);
        }
        for(int r = 0; r < int(scores1[vi].size()); r++) {
            combinedScores.push_back(scores1[vi][r]);
            combinedReadIndices.push_back(
                r < int(results[1].readIndices.size())
                    ? results[1].readIndices[r] : -1);
        }

        // Skip variants with no alt support.
        int altCount = 0;
        for(size_t i = 0; i < combinedScores.size(); i++) {
            if(combinedScores[i] == 1) altCount++;
        }
        if(altCount == 0) continue;

        // Check if variant already exists.
        KmVarKey key = varToKey(v);
        int candIdx = -1;
        for(int ci = 0; ci < int(scratch.candidates.size()); ci++) {
            if(scratch.candidates[ci].key == key) {
                candIdx = ci;
                break;
            }
        }

        if(candIdx < 0) {
            KmCandidate cand = variantToCandidate(v, combinedScores, categories[vi]);
            candIdx = int(scratch.candidates.size());
            scratch.candidates.push_back(cand);
        }

        // Update overlap profiles.
        for(size_t ri = 0; ri < combinedScores.size(); ri++) {
            int oi = combinedReadIndices[ri];
            if(oi < 0) continue;
            if(oi >= int(scratch.overlapProfiles.size())) continue;

            auto& prof = scratch.overlapProfiles[oi];
            if(prof.startVarIdx < 0) {
                prof.startVarIdx = candIdx;
                prof.endVarIdx = candIdx;
            }
            int off = candIdx - prof.startVarIdx;
            if(off < 0) continue;
            while(int(prof.alleles.size()) <= off) {
                prof.alleles.push_back(-1);
            }
            prof.alleles[off] = combinedScores[ri];
            if(candIdx > prof.endVarIdx) {
                prof.endVarIdx = candIdx;
            }
        }
    }

    return int(allVars.size());
}

// ============================================================================
// Process one MSA result (single consensus path)
// ============================================================================

/// Extract variants from one KmAbpoaMsaResult, score reads, merge into scratch.
/// Used for the combined fallback (nCons==1) path.
/// Returns number of new variants added.
static int processOneMsaResult(
    KmScratchpad& scratch,
    const KmAbpoaMsaResult& result,
    KmVariantCategory category)
{
    if(result.nReads == 0 || result.msaLen == 0) return 0;

    vector<NoisyMsaVariant> vars = extractVariantsFromMsa(
        result, result.backboneStartPos);
    if(vars.empty()) return 0;

    vector<vector<int8_t>> readScores = scoreReadsAtVariants(
        result, vars);

    mergeNoisyVariants(scratch, vars, readScores,
                       result.readIndices, category);

    return int(vars.size());
}

// ============================================================================
// kmNoisyMsaStep4 — outer loop (pgphase collect_noisy_vars_step4)
// ============================================================================

void dinara::kmNoisyMsaStep4(
    const Assembler& assembler,
    KmScratchpad& scratch,
    const KmNoisyMsaOptions& opts)
{
    if(scratch.noisyRegions.empty()) return;

    int totalNewVars = 0;
    int regionsProcessed = 0;

    for(size_t nri = 0; nri < scratch.noisyRegions.size(); nri++) {
        const auto& reg = scratch.noisyRegions[nri];
        const uint32_t regStart = reg.start;
        const uint32_t regEnd = reg.end;

        if(regEnd <= regStart) continue;
        if(int(regEnd - regStart) > opts.maxNoisyRegLen) continue;

        // Run MSA (per-haplotype or combined fallback).
        array<KmAbpoaMsaResult, 2> results;
        int nCons = collectNoisyRegMsa(
            assembler, scratch, regStart, regEnd, opts, results);

        if(nCons == 0) continue;
        regionsProcessed++;

        if(nCons == 2) {
            // Two haplotype consensuses — extract variants from both,
            // cross-score reads from each haplotype at all variant positions.
            totalNewVars += processTwoHapResults(scratch, results);
        } else if(nCons == 1) {
            // Single consensus (from combined POA or single haplotype).
            // Variants are homozygous candidates — the subsequent k-means
            // re-run will use these profiles to separate haplotypes.
            totalNewVars += processOneMsaResult(
                scratch, results[0],
                KmVariantCategory::NoisyCandHom);
        }
    }

    if(totalNewVars > 0 || regionsProcessed > 0) {
        cout << "Noisy-region MSA: processed " << regionsProcessed
             << "/" << scratch.noisyRegions.size()
             << " regions, recovered " << totalNewVars
             << " new variant sites." << endl;
    }
}
