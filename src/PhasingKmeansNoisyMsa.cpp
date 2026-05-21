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
#include <numeric>
#include <vector>

using namespace dinara;
using namespace std;

// ============================================================================
// Constants
// ============================================================================

static constexpr uint8_t kGap = 4;
static constexpr int kLongIndelLen = 10;
static constexpr float kConsSimilarityThreshold = 0.9f;

// ============================================================================
// Variant extraction from MSA matrix
// ============================================================================

struct NoisyMsaVariant {
    uint32_t backbonePos;
    KmVarType type;
    uint8_t refBase = 0;
    uint8_t altBase = 0;
    string refSeq;
    string altSeq;
    uint32_t msaColStart;
    uint32_t msaColEnd;
};

static KmVarKey varToKey(const NoisyMsaVariant& v)
{
    KmVarKey key;
    key.pos = v.backbonePos;
    key.type = v.type;
    key.altBase = v.altBase;
    if(v.type == KmVarType::Snp) { key.refLen = 1; key.altLen = 1; }
    else if(v.type == KmVarType::Insertion) { key.refLen = 0; key.altLen = uint32_t(v.altSeq.size()); }
    else { key.refLen = uint32_t(v.refSeq.size()); key.altLen = 0; }
    return key;
}

/// Walk the MSA consensus vs backbone rows to extract variants.
/// Port of pgphase make_cand_vars_from_msa / make_cand_vars_from_baln0.
static vector<NoisyMsaVariant> extractVariantsFromMsa(
    const KmAbpoaMsaResult& result, uint32_t backboneStartPos)
{
    vector<NoisyMsaVariant> vars;
    if(result.msaLen == 0) return vars;
    const auto& bbRow = result.backboneMsaRow;
    const auto& consRow = result.consensusMsaRow;
    uint32_t bbPos = backboneStartPos;
    uint32_t col = 0;

    while(col < uint32_t(result.msaLen)) {
        uint8_t bbVal = bbRow[col];
        uint8_t consVal = consRow[col];

        if(bbVal == consVal) {
            if(bbVal != kGap) bbPos++;
            col++;
            continue;
        }

        // SNP: suppress if next column has gap in either row (complex boundary).
        if(bbVal < 4 && consVal < 4) {
            const bool nextBbNonGap = (col + 1 >= uint32_t(result.msaLen)) || bbRow[col + 1] != kGap;
            const bool nextConsNonGap = (col + 1 >= uint32_t(result.msaLen)) || consRow[col + 1] != kGap;
            if(nextBbNonGap && nextConsNonGap) {
                NoisyMsaVariant v;
                v.backbonePos = bbPos; v.type = KmVarType::Snp;
                v.refBase = bbVal; v.altBase = consVal;
                v.msaColStart = col; v.msaColEnd = col + 1;
                vars.push_back(v);
            }
            bbPos++; col++;
            continue;
        }

        // Insertion: backbone gap, consensus base.
        if(bbVal == kGap && consVal < 4) {
            NoisyMsaVariant v;
            v.backbonePos = bbPos; v.type = KmVarType::Insertion;
            v.msaColStart = col;
            while(col < uint32_t(result.msaLen) && bbRow[col] == kGap && consRow[col] < 4) {
                v.altSeq.push_back("ACGT"[consRow[col]]);
                col++;
            }
            v.msaColEnd = col;
            vars.push_back(v);
            continue;
        }

        // Deletion: backbone base, consensus gap.
        if(bbVal < 4 && consVal == kGap) {
            NoisyMsaVariant v;
            v.backbonePos = bbPos; v.type = KmVarType::Deletion;
            v.msaColStart = col;
            while(col < uint32_t(result.msaLen) && bbRow[col] < 4 && consRow[col] == kGap) {
                v.refSeq.push_back("ACGT"[bbRow[col]]);
                bbPos++; col++;
            }
            v.msaColEnd = col;
            vars.push_back(v);
            continue;
        }

        if(bbVal != kGap) bbPos++;
        col++;
    }
    return vars;
}

// ============================================================================
// Partial-cover check
// ============================================================================

/// Check if a read fully covers a variant's MSA columns by verifying
/// the read has non-gap bases at the nearest backbone-base columns flanking
/// the variant span. Mirrors pgphase's full_cover (cover_start && cover_end).
static bool readFullyCoversVariant(
    const vector<uint8_t>& readRow,
    const vector<uint8_t>& bbRow,
    int msaLen,
    uint32_t msaColStart,
    uint32_t msaColEnd)
{
    // Left flank: find nearest backbone-base column before the variant.
    bool leftCover = false;
    for(int col = int(msaColStart) - 1; col >= 0; col--) {
        if(bbRow[col] != kGap) { leftCover = (readRow[col] < 4); break; }
    }
    if(!leftCover && msaColStart == 0) leftCover = true;

    // Right flank: find nearest backbone-base column after the variant.
    bool rightCover = false;
    for(int col = int(msaColEnd); col < msaLen; col++) {
        if(bbRow[col] != kGap) { rightCover = (readRow[col] < 4); break; }
    }
    if(!rightCover && int(msaColEnd) >= msaLen) rightCover = true;

    return leftCover && rightCover;
}

// ============================================================================
// Score reads at variant positions
// ============================================================================

/// For each variant, score each read in the MSA: ref (0), alt (1), or
/// unknown (-1). Reads that don't fully cover the variant get -1.
/// Insertions use consensus-matching with similarity threshold for long indels.
static vector<vector<int8_t>> scoreReadsAtVariants(
    const KmAbpoaMsaResult& result,
    const vector<NoisyMsaVariant>& vars)
{
    const int nReads = result.nReads;
    vector<vector<int8_t>> scores(vars.size());
    if(nReads == 0) return scores;

    for(size_t vi = 0; vi < vars.size(); vi++)
        scores[vi].assign(nReads, -1);

    for(size_t vi = 0; vi < vars.size(); vi++) {
        const auto& v = vars[vi];

        if(v.type == KmVarType::Snp) {
            for(int r = 0; r < nReads; r++) {
                uint8_t readVal = result.readMsaRows[r][v.msaColStart];
                if(readVal == v.refBase) scores[vi][r] = 0;
                else if(readVal == v.altBase) scores[vi][r] = 1;
            }
        }
        else if(v.type == KmVarType::Insertion) {
            int span = int(v.msaColEnd - v.msaColStart);
            for(int r = 0; r < nReads; r++) {
                if(!readFullyCoversVariant(result.readMsaRows[r],
                        result.backboneMsaRow, result.msaLen,
                        v.msaColStart, v.msaColEnd))
                    continue;
                int nMatch = 0, nMismatch = 0, nGaps = 0;
                for(uint32_t col = v.msaColStart; col < v.msaColEnd; col++) {
                    uint8_t rv = result.readMsaRows[r][col];
                    uint8_t cv = result.consensusMsaRow[col];
                    if(rv < 4 && rv == cv) nMatch++;
                    else if(rv < 4) nMismatch++;
                    else if(rv == kGap) nGaps++;
                }
                int nBases = nMatch + nMismatch;
                if(nGaps == span) {
                    scores[vi][r] = 0;  // all gaps -> ref
                } else if(nBases == 0) {
                    continue;  // no bases at all -> unknown
                } else if(span >= kLongIndelLen) {
                    if(nMatch >= int(span * kConsSimilarityThreshold))
                        scores[vi][r] = 1;
                    else if(nBases + nGaps == span)
                        scores[vi][r] = 0;  // fully covered, not similar -> ref
                    // else mixed/partial -> stays -1
                } else {
                    if(nMatch == span && nMismatch == 0)
                        scores[vi][r] = 1;  // exact match -> alt
                    else if(nBases + nGaps == span)
                        scores[vi][r] = 0;  // fully covered, not matching -> ref
                    // else mixed/partial -> stays -1
                }
            }
        }
        else if(v.type == KmVarType::Deletion) {
            int span = int(v.msaColEnd - v.msaColStart);
            for(int r = 0; r < nReads; r++) {
                if(!readFullyCoversVariant(result.readMsaRows[r],
                        result.backboneMsaRow, result.msaLen,
                        v.msaColStart, v.msaColEnd))
                    continue;
                int nBases = 0, nGaps = 0;
                for(uint32_t col = v.msaColStart; col < v.msaColEnd; col++) {
                    uint8_t rv = result.readMsaRows[r][col];
                    if(rv < 4) nBases++;
                    else if(rv == kGap) nGaps++;
                }
                if(nBases == span) scores[vi][r] = 0;
                else if(nGaps == span) scores[vi][r] = 1;
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
    cand.categoryFlag = kmCategoryToFlag(category);
    for(size_t i = 0; i < readScores.size(); i++) {
        int8_t s = readScores[i];
        if(s == 0) { cand.refCov++; cand.totalCov++; }
        else if(s == 1) { cand.altCov++; cand.totalCov++; }
    }
    cand.alleleFraction = cand.totalCov > 0
        ? double(cand.altCov) / double(cand.totalCov) : 0.0;
    return cand;
}

// ============================================================================
// Merge new variants into scratch (pgphase merge_var_profile)
// ============================================================================

struct NewNoisyVar {
    KmCandidate candidate;
    vector<int8_t> scores;
    vector<int> readIndices;
};

/// Sorted-merge new variants into scratch.candidates and rebuild all overlap
/// profiles with remapped indices. Mirrors pgphase merge_var_profile.
/// On collision (same KmVarKey), keeps the old candidate and discards the
/// new variant's read scores (pgphase behavior).
static void mergeNoisyVariants(
    KmScratchpad& scratch,
    vector<NewNoisyVar>& newVars)
{
    if(newVars.empty()) return;

    stable_sort(newVars.begin(), newVars.end(),
        [](const NewNoisyVar& a, const NewNoisyVar& b) {
            return a.candidate.key < b.candidate.key;
        });

    const int nOld = int(scratch.candidates.size());
    const int nNew = int(newVars.size());

    vector<KmCandidate> merged;
    merged.reserve(nOld + nNew);
    vector<int> oldToMerged(nOld, -1);
    vector<int> newToMerged(nNew, -1);

    int oi = 0, ni = 0;
    while(oi < nOld && ni < nNew) {
        const auto& oldKey = scratch.candidates[oi].key;
        const auto& newKey = newVars[ni].candidate.key;
        if(oldKey < newKey) {
            oldToMerged[oi] = int(merged.size());
            merged.push_back(scratch.candidates[oi++]);
        } else if(newKey < oldKey) {
            newToMerged[ni] = int(merged.size());
            merged.push_back(newVars[ni++].candidate);
        } else {
            // Collision: keep old, discard new scores.
            oldToMerged[oi] = int(merged.size());
            merged.push_back(scratch.candidates[oi++]);
            ni++;
        }
    }
    while(oi < nOld) {
        oldToMerged[oi] = int(merged.size());
        merged.push_back(scratch.candidates[oi++]);
    }
    while(ni < nNew) {
        newToMerged[ni] = int(merged.size());
        merged.push_back(newVars[ni++].candidate);
    }

    // Rebuild overlap profiles using old->merged index map.
    const int nOv = int(scratch.overlapProfiles.size());
    vector<KmOverlapProfile> mergedProfiles(nOv);
    for(int o = 0; o < nOv; o++)
        mergedProfiles[o].overlapIdx = scratch.overlapProfiles[o].overlapIdx;

    // Transfer old profile entries.
    for(int o = 0; o < nOv; o++) {
        const auto& oldProf = scratch.overlapProfiles[o];
        if(oldProf.startVarIdx < 0) continue;
        for(int off = 0; off < int(oldProf.alleles.size()); off++) {
            int oldIdx = oldProf.startVarIdx + off;
            if(oldIdx >= nOld) break;
            int mergedIdx = oldToMerged[oldIdx];
            if(mergedIdx < 0) continue;
            auto& mp = mergedProfiles[o];
            if(mp.startVarIdx < 0 || mergedIdx < mp.startVarIdx) {
                if(mp.startVarIdx >= 0) {
                    int shift = mp.startVarIdx - mergedIdx;
                    mp.alleles.insert(mp.alleles.begin(), shift, -1);
                }
                mp.startVarIdx = mergedIdx;
            }
            if(mp.endVarIdx < mergedIdx) mp.endVarIdx = mergedIdx;
            int moff = mergedIdx - mp.startVarIdx;
            while(int(mp.alleles.size()) <= moff) mp.alleles.push_back(-1);
            mp.alleles[moff] = oldProf.alleles[off];
        }
    }

    // Apply new variant read scores.
    for(int nvi = 0; nvi < nNew; nvi++) {
        int mergedIdx = newToMerged[nvi];
        if(mergedIdx < 0) continue;
        const auto& nv = newVars[nvi];
        for(size_t ri = 0; ri < nv.scores.size(); ri++) {
            int ovIdx = nv.readIndices[ri];
            if(ovIdx < 0 || ovIdx >= nOv) continue;
            auto& mp = mergedProfiles[ovIdx];
            if(mp.startVarIdx < 0) {
                mp.startVarIdx = mergedIdx;
                mp.endVarIdx = mergedIdx;
            }
            if(mergedIdx < mp.startVarIdx) {
                int shift = mp.startVarIdx - mergedIdx;
                mp.alleles.insert(mp.alleles.begin(), shift, -1);
                mp.startVarIdx = mergedIdx;
            }
            if(mp.endVarIdx < mergedIdx) mp.endVarIdx = mergedIdx;
            int moff = mergedIdx - mp.startVarIdx;
            while(int(mp.alleles.size()) <= moff) mp.alleles.push_back(-1);
            mp.alleles[moff] = nv.scores[ri];
        }
    }

    for(auto& vi : scratch.validVarIdx) {
        if(int(vi) < nOld) vi = uint32_t(oldToMerged[vi]);
    }

    scratch.candidates = move(merged);
    scratch.overlapProfiles = move(mergedProfiles);
}

// ============================================================================
// Cross-haplotype coverage check (pgphase get_full_cover_from_cons_aln_str)
// ============================================================================

/// Build a map from backbone base position to MSA column index.
static vector<int> buildBbPosToMsaCol(const KmAbpoaMsaResult& result)
{
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

/// Cross-haplotype coverage check. For variants from the OTHER consensus,
/// only checks if the read covers the variant position and assigns ref (0).
/// Mirrors pgphase update_cand_var_profile_from_cons_aln_str21 cross-hap branch.
static int8_t crossCoverageCheck(
    const KmAbpoaMsaResult& result, int readIdx,
    const NoisyMsaVariant& v, const vector<int>& posToCol)
{
    int relPos = int(v.backbonePos - result.backboneStartPos);
    if(relPos < 0 || relPos >= int(posToCol.size())) return -1;
    int col = posToCol[relPos];
    if(col < 0 || col >= result.msaLen) return -1;

    if(v.type == KmVarType::Deletion) {
        int delLen = int(v.refSeq.size());
        vector<int> delCols;
        for(int dc = col; dc < result.msaLen && int(delCols.size()) < delLen; dc++) {
            if(result.backboneMsaRow[dc] != kGap) delCols.push_back(dc);
        }
        if(int(delCols.size()) != delLen) return -1;
        uint32_t delColStart = uint32_t(delCols.front());
        uint32_t delColEnd = uint32_t(delCols.back()) + 1;
        if(!readFullyCoversVariant(result.readMsaRows[readIdx],
                result.backboneMsaRow, result.msaLen, delColStart, delColEnd))
            return -1;
        return 0;
    }

    // SNP or insertion: check if read covers the anchor position.
    uint8_t readVal = result.readMsaRows[readIdx][col];
    return (readVal < 4) ? 0 : -1;
}

// ============================================================================
// Process per-haplotype MSA results
// ============================================================================

/// Two-haplotype path. Port of pgphase update_cand_var_profile_from_cons_aln_str2.
/// Builds unified variant list with varFromCons bitmask (1=hap0, 2=hap1, 3=both).
/// Same-haplotype: full scoring via scoreReadsAtVariants.
/// Cross-haplotype: coverage-only ref via crossCoverageCheck.
static int processTwoHapResults(
    KmScratchpad& scratch,
    const array<KmAbpoaMsaResult, 2>& results)
{
    vector<NoisyMsaVariant> vars0 = extractVariantsFromMsa(
        results[0], results[0].backboneStartPos);
    vector<NoisyMsaVariant> vars1 = extractVariantsFromMsa(
        results[1], results[1].backboneStartPos);
    if(vars0.empty() && vars1.empty()) return 0;

    // Sorted merge with varFromCons bitmask.
    vector<NoisyMsaVariant> allVars;
    vector<int> varFromCons;
    vector<KmVariantCategory> categories;
    allVars.reserve(vars0.size() + vars1.size());
    varFromCons.reserve(vars0.size() + vars1.size());
    categories.reserve(vars0.size() + vars1.size());

    size_t i0 = 0, i1 = 0;
    while(i0 < vars0.size() && i1 < vars1.size()) {
        KmVarKey k0 = varToKey(vars0[i0]);
        KmVarKey k1 = varToKey(vars1[i1]);
        if(k0 < k1) {
            allVars.push_back(vars0[i0++]);
            varFromCons.push_back(1);
            categories.push_back(KmVariantCategory::NoisyCandHet);
        } else if(k1 < k0) {
            allVars.push_back(vars1[i1++]);
            varFromCons.push_back(2);
            categories.push_back(KmVariantCategory::NoisyCandHet);
        } else {
            allVars.push_back(vars0[i0++]);
            varFromCons.push_back(3);
            categories.push_back(KmVariantCategory::NoisyCandHom);
            i1++;
        }
    }
    while(i0 < vars0.size()) {
        allVars.push_back(vars0[i0++]);
        varFromCons.push_back(1);
        categories.push_back(KmVariantCategory::NoisyCandHet);
    }
    while(i1 < vars1.size()) {
        allVars.push_back(vars1[i1++]);
        varFromCons.push_back(2);
        categories.push_back(KmVariantCategory::NoisyCandHet);
    }
    if(allVars.empty()) return 0;

    // Full scoring for same-haplotype variants.
    vector<vector<int8_t>> sameScores0 = scoreReadsAtVariants(results[0], vars0);
    vector<vector<int8_t>> sameScores1 = scoreReadsAtVariants(results[1], vars1);

    // posToCol maps for cross-coverage checks.
    vector<int> posToCol0 = buildBbPosToMsaCol(results[0]);
    vector<int> posToCol1 = buildBbPosToMsaCol(results[1]);

    // Map unified variant index to same-haplotype score index.
    vector<int> sameIdx0(allVars.size(), -1);
    vector<int> sameIdx1(allVars.size(), -1);
    {
        int si0 = 0, si1 = 0;
        for(size_t vi = 0; vi < allVars.size(); vi++) {
            if(varFromCons[vi] & 1) sameIdx0[vi] = si0++;
            if(varFromCons[vi] & 2) sameIdx1[vi] = si1++;
        }
    }

    // Build combined readIndices.
    const int nReads0 = results[0].nReads;
    const int nReads1 = results[1].nReads;
    const int nReadsTotal = nReads0 + nReads1;
    vector<int> combinedReadIndices(nReadsTotal);
    for(int r = 0; r < nReads0; r++)
        combinedReadIndices[r] = r < int(results[0].readIndices.size())
            ? results[0].readIndices[r] : -1;
    for(int r = 0; r < nReads1; r++)
        combinedReadIndices[nReads0 + r] = r < int(results[1].readIndices.size())
            ? results[1].readIndices[r] : -1;

    // Score per variant.
    vector<NewNoisyVar> newVars;
    newVars.reserve(allVars.size());

    for(size_t vi = 0; vi < allVars.size(); vi++) {
        vector<int8_t> varScores(nReadsTotal);

        // Hap0 reads.
        for(int r = 0; r < nReads0; r++) {
            if(varFromCons[vi] & 1) {
                int si = sameIdx0[vi];
                varScores[r] = (si >= 0 && si < int(sameScores0.size()))
                    ? sameScores0[si][r] : -1;
            } else {
                varScores[r] = crossCoverageCheck(results[0], r, allVars[vi], posToCol0);
            }
        }
        // Hap1 reads.
        for(int r = 0; r < nReads1; r++) {
            if(varFromCons[vi] & 2) {
                int si = sameIdx1[vi];
                varScores[nReads0 + r] = (si >= 0 && si < int(sameScores1.size()))
                    ? sameScores1[si][r] : -1;
            } else {
                varScores[nReads0 + r] = crossCoverageCheck(results[1], r, allVars[vi], posToCol1);
            }
        }

        int altCount = 0;
        for(auto s : varScores) { if(s == 1) altCount++; }
        if(altCount == 0) continue;

        NewNoisyVar nv;
        nv.candidate = variantToCandidate(allVars[vi], varScores, categories[vi]);
        nv.scores = move(varScores);
        nv.readIndices = combinedReadIndices;
        newVars.push_back(move(nv));
    }

    int nAdded = int(newVars.size());
    mergeNoisyVariants(scratch, newVars);
    return nAdded;
}

/// Single-consensus path. Port of pgphase update_cand_var_profile_from_cons_aln_str1.
static int processOneMsaResult(
    KmScratchpad& scratch,
    const KmAbpoaMsaResult& result,
    KmVariantCategory category)
{
    if(result.nReads == 0 || result.msaLen == 0) return 0;

    vector<NoisyMsaVariant> vars = extractVariantsFromMsa(
        result, result.backboneStartPos);
    if(vars.empty()) return 0;

    vector<vector<int8_t>> allScores = scoreReadsAtVariants(result, vars);

    vector<NewNoisyVar> newVars;
    newVars.reserve(vars.size());

    for(size_t vi = 0; vi < vars.size(); vi++) {
        vector<int8_t> varScores(result.nReads);
        for(int r = 0; r < result.nReads; r++)
            varScores[r] = allScores[vi][r];

        int altCount = 0;
        for(auto s : varScores) { if(s == 1) altCount++; }
        if(altCount == 0) continue;

        NewNoisyVar nv;
        nv.candidate = variantToCandidate(vars[vi], varScores, category);
        nv.scores = move(varScores);
        nv.readIndices = result.readIndices;
        newVars.push_back(move(nv));
    }

    int nAdded = int(newVars.size());
    mergeNoisyVariants(scratch, newVars);
    return nAdded;
}

// ============================================================================
// sortNoisyRegs (pgphase sort_noisy_regs)
// ============================================================================

static vector<int> sortNoisyRegs(const vector<KmNoisyRegion>& regions)
{
    const int n = int(regions.size());
    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);
    for(int i = 0; i < n; i++) {
        for(int j = i + 1; j < n; j++) {
            const auto& a = regions[idx[i]];
            const auto& b = regions[idx[j]];
            bool doSwap = false;
            if(a.label > b.label) doSwap = true;
            else if(a.label == b.label) {
                if((a.end - a.start) > (b.end - b.start)) doSwap = true;
            }
            if(doSwap) swap(idx[i], idx[j]);
        }
    }
    return idx;
}

// ============================================================================
// collectNoisyVars1 (pgphase collect_noisy_vars1)
// ============================================================================

static int collectNoisyVars1(
    const Assembler& assembler, KmScratchpad& scratch,
    const KmNoisyMsaOptions& opts, int regIdx)
{
    const auto& reg = scratch.noisyRegions[regIdx];
    const uint32_t regStart = reg.start;
    const uint32_t regEnd = reg.end;
    if(regEnd <= regStart) return 0;
    if(int(regEnd - regStart) > opts.maxNoisyRegLen) return 0;

    array<KmAbpoaMsaResult, 2> results;
    int nCons = collectNoisyRegMsa(
        assembler, scratch, regStart, regEnd, opts, results);

    if(nCons == 0) return -1;

    if(nCons == 2)
        return processTwoHapResults(scratch, results);
    if(nCons == 1)
        return processOneMsaResult(scratch, results[0],
                                   KmVariantCategory::NoisyCandHom);
    return 0;
}

// ============================================================================
// kmNoisyMsaStep4 (pgphase collect_noisy_vars_step4)
// ============================================================================

void dinara::kmNoisyMsaStep4(
    const Assembler& assembler,
    KmScratchpad& scratch,
    const KmNoisyMsaOptions& msaOpts,
    const KmPhasingOptions& phasingOpts)
{
    if(scratch.noisyRegions.empty()) return;

    const vector<int> sorted = sortNoisyRegs(scratch.noisyRegions);
    const int nRegs = int(scratch.noisyRegions.size());
    vector<bool> done(nRegs, false);

    int totalNewVars = 0;
    int totalPasses = 0;

    while(true) {
        bool anyDone = false;
        bool anyNewVar = false;

        for(int regIdx : sorted) {
            if(done[regIdx]) continue;
            const int ret = collectNoisyVars1(
                assembler, scratch, msaOpts, regIdx);
            if(ret >= 0) {
                done[regIdx] = true;
                anyDone = true;
                if(ret > 0) {
                    anyNewVar = true;
                    totalNewVars += ret;
                }
            }
        }

        if(anyNewVar)
            kmRunKmeans(scratch, phasingOpts, KM_GERMLINE_ALL);

        totalPasses++;
        if(!anyDone) break;
    }

    if(totalNewVars > 0 || totalPasses > 1) {
        cout << "Noisy-region MSA: " << totalPasses << " pass(es), recovered "
             << totalNewVars << " new variant sites from "
             << nRegs << " regions." << endl;
    }
}
