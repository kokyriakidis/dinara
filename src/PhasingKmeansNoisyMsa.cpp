/// @file PhasingKmeansNoisyMsa.cpp
/// @brief Noisy-region MSA variant calling and re-phasing (pgphase step 4).
///
/// Port of pgphase/longcallD collect_noisy_vars_step4.
///
/// Key design (matching pgphase):
/// - Variant extraction: pairwise alignment of raw backbone vs consensus
///   sequences (edlib NW), NOT direct MSA row comparison. This avoids
///   multi-sequence alignment artifacts that inflate variant counts.
/// - Read scoring: direct MSA column indexing with full_cover gating.
///   Each variant carries msaColStart/msaColEnd mapped back from the
///   pairwise alignment, so reads are scored from the original MSA matrix.

#include "PhasingKmeansAlign.hpp"
#include "PhasingKmeansTypes.hpp"
#include "Assembler.hpp"
#include "edlib.h"

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

static constexpr uint8_t kGap = 5;  // abPOA MSA gap value (abpt->m for nucleotide mode)
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
    else if(v.type == KmVarType::Insertion) {
        key.refLen = 0;
        key.altLen = uint16_t(v.altSeq.size());
        key.altSeq = v.altSeq;
    }
    else { key.refLen = uint16_t(v.refSeq.size()); key.altLen = 0; }
    return key;
}

/// Extract raw (ungapped) sequence from an MSA row.
static vector<uint8_t> extractRawSeq(const vector<uint8_t>& msaRow, int msaLen)
{
    vector<uint8_t> seq;
    seq.reserve(msaLen);
    for(int i = 0; i < msaLen; i++) {
        if(msaRow[i] != kGap) seq.push_back(msaRow[i]);
    }
    return seq;
}

/// Pairwise-align backbone vs consensus sequences using edlib, then extract
/// variants from the alignment. Port of pgphase's approach: WFA ref-vs-cons
/// alignment followed by make_cand_vars_from_baln0.
///
/// Using pairwise alignment instead of MSA row comparison avoids
/// multi-sequence alignment artifacts that inflate variant counts.
static vector<NoisyMsaVariant> extractVariantsFromPairwise(
    const vector<uint8_t>& bbRaw,
    const vector<uint8_t>& consRaw,
    uint32_t backboneStartPos)
{
    vector<NoisyMsaVariant> vars;
    if(bbRaw.empty() || consRaw.empty()) return vars;

    // edlib global (NW) alignment with full path.
    EdlibAlignResult result = edlibAlign(
        reinterpret_cast<const char*>(consRaw.data()), int(consRaw.size()),
        reinterpret_cast<const char*>(bbRaw.data()), int(bbRaw.size()),
        edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, nullptr, 0));

    if(result.status != EDLIB_STATUS_OK || result.alignment == nullptr) {
        edlibFreeAlignResult(result);
        return vars;
    }

    // Build aligned rows from edlib edit operations.
    // edlib: query=cons, target=bb.
    // EDLIB_EDOP_MATCH(0)    = both present, same base
    // EDLIB_EDOP_INSERT(1)   = insertion to target = bb has base, cons has gap
    // EDLIB_EDOP_DELETE(2)   = deletion from target = cons has base, bb has gap
    // EDLIB_EDOP_MISMATCH(3) = both present, different base
    vector<uint8_t> alnRef, alnCons;
    alnRef.reserve(result.alignmentLength);
    alnCons.reserve(result.alignmentLength);
    int ri = 0, ci = 0;
    for(int k = 0; k < result.alignmentLength; k++) {
        unsigned char op = result.alignment[k];
        switch(op) {
            case EDLIB_EDOP_MATCH:
            case EDLIB_EDOP_MISMATCH:
                alnRef.push_back(bbRaw[ri++]);
                alnCons.push_back(consRaw[ci++]);
                break;
            case EDLIB_EDOP_INSERT: // insertion to target: bb base, cons gap
                alnRef.push_back(bbRaw[ri++]);
                alnCons.push_back(kGap);
                break;
            case EDLIB_EDOP_DELETE: // deletion from target: cons base, bb gap
                alnRef.push_back(kGap);
                alnCons.push_back(consRaw[ci++]);
                break;
        }
    }
    edlibFreeAlignResult(result);

    // Extract variants from aligned rows (pgphase make_cand_vars_from_baln0).
    const char* bases = "ACGTN";
    const int alnLen = int(alnRef.size());
    uint32_t bbPos = backboneStartPos;
    int i = 0;

    while(i < alnLen) {
        uint8_t rv = alnRef[i];
        uint8_t cv = alnCons[i];

        // Skip double-gap columns (shouldn't happen with edlib but be safe).
        if(rv == kGap && cv == kGap) { i++; continue; }

        // Match/same.
        if(rv == cv) { if(rv != kGap) bbPos++; i++; continue; }

        // SNP: both non-gap, different.
        if(rv != kGap && cv != kGap) {
            // pgphase: suppress SNP if next column has gap in either row.
            bool nextRefNonGap = (i+1 >= alnLen) || alnRef[i+1] != kGap;
            bool nextConsNonGap = (i+1 >= alnLen) || alnCons[i+1] != kGap;
            if(nextRefNonGap && nextConsNonGap) {
                NoisyMsaVariant v;
                v.backbonePos = bbPos;
                v.type = KmVarType::Snp;
                v.refBase = rv;
                v.altBase = cv;
                vars.push_back(v);
            }
            bbPos++; i++;
            continue;
        }

        // Insertion: ref gap, cons base.
        if(rv == kGap) {
            NoisyMsaVariant v;
            v.backbonePos = bbPos;
            v.type = KmVarType::Insertion;
            while(i < alnLen && alnRef[i] == kGap && alnCons[i] != kGap) {
                v.altSeq.push_back(bases[alnCons[i] < 4 ? alnCons[i] : 4]);
                i++;
            }
            vars.push_back(v);
            continue;
        }

        // Deletion: ref base, cons gap.
        if(cv == kGap) {
            NoisyMsaVariant v;
            v.backbonePos = bbPos;
            v.type = KmVarType::Deletion;
            while(i < alnLen && alnRef[i] != kGap && alnCons[i] == kGap) {
                v.refSeq.push_back(bases[alnRef[i] < 4 ? alnRef[i] : 4]);
                bbPos++; i++;
            }
            vars.push_back(v);
            continue;
        }

        if(rv != kGap) bbPos++;
        i++;
    }
    return vars;
}

/// Top-level variant extraction: extract raw sequences from MSA, pairwise
/// align, and extract variants. Also populates msaColStart/msaColEnd on each
/// variant by mapping backbone positions back to MSA columns (needed for
/// read scoring which still uses the MSA matrix).
static vector<NoisyMsaVariant> extractVariantsFromMsa(
    const KmAbpoaMsaResult& result, uint32_t backboneStartPos)
{
    if(result.msaLen == 0) return {};

    // Extract raw sequences.
    vector<uint8_t> bbRaw = extractRawSeq(result.backboneMsaRow, result.msaLen);
    vector<uint8_t> consRaw = extractRawSeq(result.consensusMsaRow, result.msaLen);

    // Pairwise alignment and variant extraction.
    vector<NoisyMsaVariant> vars = extractVariantsFromPairwise(
        bbRaw, consRaw, backboneStartPos);

    // Build backbone-position-to-MSA-column map for read scoring.
    // We need msaColStart/msaColEnd for each variant so scoreReadsAtVariants
    // can index into the MSA matrix.
    vector<int> bbPosToCol;
    {
        int maxBbPos = 0;
        int bbIdx = 0;
        for(int col = 0; col < result.msaLen; col++) {
            if(result.backboneMsaRow[col] != kGap) bbIdx++;
        }
        maxBbPos = bbIdx;
        bbPosToCol.resize(maxBbPos, -1);
        bbIdx = 0;
        for(int col = 0; col < result.msaLen; col++) {
            if(result.backboneMsaRow[col] != kGap) {
                if(bbIdx < maxBbPos) bbPosToCol[bbIdx] = col;
                bbIdx++;
            }
        }
    }

    // Map variant backbone positions to MSA columns.
    for(auto& v : vars) {
        int relPos = int(v.backbonePos - backboneStartPos);
        if(relPos < 0 || relPos >= int(bbPosToCol.size())) {
            v.msaColStart = v.msaColEnd = 0;
            continue;
        }

        if(v.type == KmVarType::Snp) {
            int col = bbPosToCol[relPos];
            v.msaColStart = uint32_t(col);
            v.msaColEnd = uint32_t(col) + 1;
        }
        else if(v.type == KmVarType::Insertion) {
            // Insertion columns are the gap columns in the backbone row
            // immediately after the anchor position.
            int anchorCol = bbPosToCol[relPos];
            // Walk backwards: insertion is BEFORE the anchor in the MSA.
            // Actually, insertion columns are gap columns in backbone row
            // at or after the anchor. Find the run of backbone-gap columns
            // after the previous backbone base.
            int startCol = anchorCol;
            // Search backwards for insertion gap columns before this anchor.
            // In the MSA, insertions appear as backbone-gap columns before
            // the next backbone base.
            if(relPos > 0) {
                int prevCol = bbPosToCol[relPos - 1];
                startCol = prevCol + 1;
            } else {
                // Insertion before first backbone base.
                startCol = 0;
            }
            // Find the gap run.
            int gapStart = -1, gapEnd = -1;
            for(int col = startCol; col < anchorCol; col++) {
                if(result.backboneMsaRow[col] == kGap) {
                    if(gapStart < 0) gapStart = col;
                    gapEnd = col + 1;
                }
            }
            if(gapStart >= 0) {
                v.msaColStart = uint32_t(gapStart);
                v.msaColEnd = uint32_t(gapEnd);
            } else {
                // Fallback: look after the anchor.
                gapStart = anchorCol + 1;
                gapEnd = gapStart;
                for(int col = gapStart; col < result.msaLen; col++) {
                    if(result.backboneMsaRow[col] == kGap)
                        gapEnd = col + 1;
                    else break;
                }
                v.msaColStart = uint32_t(gapStart);
                v.msaColEnd = uint32_t(gapEnd);
            }
        }
        else if(v.type == KmVarType::Deletion) {
            int delLen = int(v.refSeq.size());
            v.msaColStart = uint32_t(bbPosToCol[relPos]);
            if(relPos + delLen - 1 < int(bbPosToCol.size()))
                v.msaColEnd = uint32_t(bbPosToCol[relPos + delLen - 1]) + 1;
            else
                v.msaColEnd = uint32_t(result.msaLen);
        }
    }

    // Remove variants with invalid MSA column mappings.
    vars.erase(
        remove_if(vars.begin(), vars.end(),
            [](const NoisyMsaVariant& v) {
                return v.msaColStart >= v.msaColEnd;
            }),
        vars.end());

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
                // pgphase full_cover gating: only score reads that span
                // the SNP position (have non-gap flanking bases).
                if(!readFullyCoversVariant(result.readMsaRows[r],
                        result.backboneMsaRow, result.msaLen,
                        v.msaColStart, v.msaColEnd))
                    continue;
                uint8_t readVal = result.readMsaRows[r][v.msaColStart];
                if(readVal == v.altBase) scores[vi][r] = 1;
                else scores[vi][r] = 0;  // ref, other base, or gap at covered pos → ref
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
            for(int r = 0; r < nReads; r++) {
                if(!readFullyCoversVariant(result.readMsaRows[r],
                        result.backboneMsaRow, result.msaLen,
                        v.msaColStart, v.msaColEnd))
                    continue;
                int nNonGap = 0;
                for(uint32_t col = v.msaColStart; col < v.msaColEnd; col++) {
                    if(result.readMsaRows[r][col] != kGap) nNonGap++;
                }
                // pgphase: any non-gap base in deletion region → ref (0),
                // all gaps → alt (1).
                scores[vi][r] = (nNonGap == 0) ? 1 : 0;
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
        vector<int8_t> varScores(nReadsTotal, -1);

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
        vector<int8_t>& varScores = allScores[vi];

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
