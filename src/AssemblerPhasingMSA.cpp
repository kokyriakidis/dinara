/// @file AssemblerPhasingMSA.cpp
/// @brief POA-based overlap phasing using Theseus MSA on anchor windows.
/// Coexists with phaseOverlapsKmeans — both write hifiasmEcMatchState.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "PhasingKmeansTypes.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "TheseusAlignMutex.hpp"
#include "invalid.hpp"
#include "timestamp.hpp"

#include <theseus/graph.h>
#include <theseus/heuristics.h>
#include <theseus/penalties.h>
#include <theseus/theseus_msa_aligner.h>

#include <algorithm>
#include <atomic>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <iostream>
#include <mutex>
#include <numeric>
#include <sstream>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// --- Tuning constants ---
constexpr uint64_t msaMinAnchorCoverage = 6;
constexpr uint64_t msaMaxReadsPerPair = 200;
constexpr uint32_t msaMaxFocalSeqLen = 5000;
constexpr uint32_t msaMinSegmentBases = 500;
constexpr uint64_t msaMinSnpRefSupport = 3;
constexpr uint64_t msaMinSnpAltSupport = 3;
constexpr uint64_t msaMinHomopolymerRun = 3;
constexpr double   msaMinAf = 0.20;
constexpr double   msaMaxAf = 0.80;

// Per-row noisy region detection (mirrors k-means KmNoisyBuilder).
// For each MSA row, a sliding window scans for dense disagreements.
// Per-row noisy intervals are merged into chunk-level regions requiring
// ≥ msaNoisyMinReads rows and ≥ msaNoisyMinRatio of spanning rows.
constexpr int      msaNoisyWindow = 25;       // sliding window size in backbone offsets (ONT)
constexpr int      msaNoisyMaxEvents = 5;     // max event weight before flagging (noisyRegMaxXgaps)
constexpr uint32_t msaNoisyMergeDis = 500;    // merge distance for per-row intervals
constexpr uint32_t msaNoisyMinReads = 3;      // min rows with noisy intervals
constexpr double   msaNoisyMinRatio = 0.20;   // min fraction of spanning rows that are noisy
constexpr uint32_t msaNoisyFlankBp = 10;      // extend each side

// --- Structured variant site ---
struct MsaAltAllele {
    string type;       // "SNP", "INS", "DEL", "MNP"
    string sequence;
    vector<OrientedReadId> reads;
};

struct MsaVariantSite {
    uint32_t backbonePosition;
    string refAllele;
    vector<MsaAltAllele> altAlleles;
    vector<OrientedReadId> refReads;
    uint32_t totalCov = 0;
    uint32_t altCov = 0;
    double af = 0.0;
    bool isInsertionColumn = false; // true for pass 2 (insertion-column) sites
};

// --- Per-read allele observation across all sites in a window ---
struct MsaReadProfile {
    OrientedReadId readId;
    vector<int> alleles;  // one per het site: 0=ref, 1=alt, -1=missing
    int firstSiteIdx = -1;
    int lastSiteIdx = -1;
};

// --- MSA sequence info (mirrors prototype) ---
struct MsaSeqInfo {
    OrientedReadId oid;
    string sequence;
    uint32_t begin = 0;
    uint32_t end = 0;
    bool hasBothAnchors = false;
    char anchorSide = 'B';
};

// --- Counters ---
struct MsaPhasingCounters {
    atomic<uint64_t> windowsProcessed{0};
    atomic<uint64_t> pairsProcessed{0};
    atomic<uint64_t> sitesDetected{0};
    atomic<uint64_t> hetSitesUsed{0};
    atomic<uint64_t> readsPhased{0};
    atomic<uint64_t> cisCount{0};
    atomic<uint64_t> transCount{0};
    atomic<uint64_t> skippedLowCov{0};
    atomic<uint64_t> skippedEmptyFocal{0};
    atomic<uint64_t> skippedLongFocal{0};
};

// --- Helpers ---

static inline bool msaIsBase(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

static string msaExtractSegment(
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t k, OrientedReadId oid,
    uint32_t ordA, uint32_t ordB)
{
    if (ordA >= ordB) return {};
    const auto rm = markers[oid.getValue()];
    if (ordB >= rm.size()) return {};
    const uint32_t kh = uint32_t(k / 2);
    const uint32_t beg = rm[ordA].position + kh;
    const uint32_t end = rm[ordB].position + kh;
    if (end <= beg) return {};
    string s;
    s.reserve(end - beg);
    for (uint32_t p = beg; p < end; p++)
        s.push_back(reads.getOrientedReadBase(oid, p).character());
    return s;
}

static bool msaIsHomopolymerAt(const string& seq, uint32_t pos) {
    if (pos >= seq.size()) return false;
    char base = seq[pos];
    if (!msaIsBase(base)) return false;
    size_t b = pos, e = pos + 1;
    while (b > 0 && seq[b-1] == base) --b;
    while (e < seq.size() && seq[e] == base) ++e;
    return (e - b) >= msaMinHomopolymerRun;
}

static bool msaSnpTouchesHomopolymer(const string& seq, uint32_t pos, char alt) {
    if (msaIsHomopolymerAt(seq, pos)) return true;
    if (pos >= seq.size() || !msaIsBase(alt)) return false;
    string tmp = seq;
    tmp[pos] = alt;
    return msaIsHomopolymerAt(tmp, pos);
}

static string msaAlleleType(const string& ref, const string& alt) {
    if (ref.size() == 1 && alt.size() == 1) return "SNP";
    if (ref.size() == alt.size()) return "MNP";
    if (ref.size() < alt.size()) return "INS";
    return "DEL";
}

static vector<string> msaParseFasta(const string& text) {
    vector<string> seqs;
    string cur;
    istringstream in(text);
    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!cur.empty()) { seqs.push_back(cur); cur.clear(); }
        } else {
            cur += line;
        }
    }
    if (!cur.empty()) seqs.push_back(cur);
    return seqs;
}


// Detect noisy regions from MSA, per-row sliding window (mirrors k-means CIGAR path).
// For each row, scan backbone offsets for dense disagreements using a sliding window.
// Merge per-row noisy intervals into chunk-level regions with read support filtering.
static vector<pair<uint32_t, uint32_t>> msaDetectNoisyRegions(
    const vector<MsaSeqInfo>& seqInfos,
    const vector<string>& alignedSeqs,
    const vector<size_t>& colByOff,
    const vector<size_t>& firstNG,
    const vector<size_t>& lastNG,
    const vector<uint8_t>& hasB,
    uint32_t focalBegin)
{
    vector<pair<uint32_t, uint32_t>> regions;
    if (colByOff.empty() || alignedSeqs.size() < 2) return regions;
    const size_t nOff = colByOff.size();

    // Per-row noisy intervals: flat buffer with per-row ranges.
    vector<pair<uint32_t, uint32_t>> allRowNoisy;
    vector<uint32_t> rowNoisyBegin(alignedSeqs.size());
    vector<uint32_t> rowNoisyEnd(alignedSeqs.size());

    for (size_t r = 1; r < alignedSeqs.size(); r++) {
        rowNoisyBegin[r] = uint32_t(allRowNoisy.size());
        if (!hasB[r]) { rowNoisyEnd[r] = rowNoisyBegin[r]; continue; }

        // Sliding window over backbone offsets for this row.
        // Event queue: positions and weights of disagreements.
        struct Event { uint32_t off; int weight; };
        vector<Event> events;
        int front = 0, totalWeight = 0;
        int32_t curStart = -1, curEnd = -1;

        for (size_t off = 0; off < nOff; off++) {
            size_t c = colByOff[off];
            if (c < firstNG[r] || c > lastNG[r]) continue;

            char ref = alignedSeqs[0][c];
            char b = alignedSeqs[r][c];
            if (!msaIsBase(ref)) continue;

            int w = 0;
            if (b == '-') w = 1;  // deletion
            else if (msaIsBase(b) && b != ref) w = 1;  // substitution

            // Check for insertion: gap columns between this ref col and next.
            if (off + 1 < nOff) {
                size_t nextC = colByOff[off + 1];
                int insLen = 0;
                for (size_t gc = c + 1; gc < nextC; gc++) {
                    if (alignedSeqs[0][gc] == '-' && msaIsBase(alignedSeqs[r][gc]))
                        insLen++;
                }
                if (insLen > 0) w += insLen;
            }

            if (w == 0) continue;

            uint32_t bbPos = focalBegin + uint32_t(off);
            events.push_back({bbPos, w});
            totalWeight += w;

            // Evict events outside the window.
            while (front < int(events.size()) &&
                   int64_t(events[front].off) <= int64_t(bbPos) - msaNoisyWindow) {
                totalWeight -= events[front].weight;
                front++;
            }

            if (totalWeight <= msaNoisyMaxEvents) continue;

            int32_t nStart = int32_t(events[front].off);
            int32_t nEnd = int32_t(bbPos + 1);

            if (curStart == -1) {
                curStart = nStart; curEnd = nEnd;
            } else if (nStart <= curEnd) {
                curEnd = max(curEnd, nEnd);
            } else {
                allRowNoisy.push_back({uint32_t(curStart), uint32_t(curEnd)});
                curStart = nStart; curEnd = nEnd;
            }
        }
        if (curStart != -1)
            allRowNoisy.push_back({uint32_t(curStart), uint32_t(curEnd)});
        rowNoisyEnd[r] = uint32_t(allRowNoisy.size());
    }

    if (allRowNoisy.empty()) return regions;

    // Collect all per-row intervals, merge.
    vector<pair<uint32_t, uint32_t>> merged;
    for (const auto& iv : allRowNoisy)
        merged.push_back(iv);
    sort(merged.begin(), merged.end());
    {
        vector<pair<uint32_t, uint32_t>> tmp;
        tmp.push_back(merged[0]);
        for (size_t i = 1; i < merged.size(); i++) {
            auto& last = tmp.back();
            if (merged[i].first <= last.second + msaNoisyMergeDis)
                last.second = max(last.second, merged[i].second);
            else
                tmp.push_back(merged[i]);
        }
        merged = move(tmp);
    }

    // For each merged region, count how many rows have noisy intervals
    // overlapping it, and how many rows span it.
    for (const auto& [mStart, mEnd] : merged) {
        int nTotal = 0, nNoisy = 0;
        for (size_t r = 1; r < alignedSeqs.size(); r++) {
            if (!hasB[r]) continue;
            // Does this row span the region?
            // Use firstNG/lastNG mapped to backbone offsets.
            // Approximate: check if the row's first/last non-gap columns
            // bracket the region's backbone positions.
            uint32_t rStart = focalBegin;
            uint32_t rEnd = focalBegin + uint32_t(nOff);
            // Find backbone offset of firstNG[r] and lastNG[r].
            for (size_t off = 0; off < nOff; off++) {
                if (colByOff[off] >= firstNG[r]) { rStart = focalBegin + uint32_t(off); break; }
            }
            for (size_t off = nOff; off > 0; off--) {
                if (colByOff[off - 1] <= lastNG[r]) { rEnd = focalBegin + uint32_t(off); break; }
            }
            if (rStart >= mEnd || rEnd <= mStart) continue;
            nTotal++;
            // Does this row have a noisy interval overlapping the merged region?
            for (uint32_t ni = rowNoisyBegin[r]; ni < rowNoisyEnd[r]; ni++) {
                if (allRowNoisy[ni].first < mEnd && allRowNoisy[ni].second > mStart) {
                    nNoisy++;
                    break;
                }
            }
        }
        if (nNoisy < int(msaNoisyMinReads)) continue;
        if (nTotal > 0 && double(nNoisy) / double(nTotal) < msaNoisyMinRatio) continue;
        // Extend flanks.
        uint32_t s = (mStart > msaNoisyFlankBp) ? mStart - msaNoisyFlankBp : 0;
        uint32_t e = mEnd + msaNoisyFlankBp;
        if (!regions.empty() && s <= regions.back().second)
            regions.back().second = max(regions.back().second, e);
        else
            regions.push_back({s, e});
    }

    return regions;
}

// Check if a backbone position falls within any noisy region.
static bool msaIsInNoisyRegion(
    uint32_t pos, const vector<pair<uint32_t, uint32_t>>& noisyRegions)
{
    for (const auto& [s, e] : noisyRegions)
        if (pos >= s && pos < e) return true;
    return false;
}

// Detect variant sites from an MSA of one anchor pair.
// Also detects per-row noisy regions and merges them into chunk-level regions.
static vector<MsaVariantSite> msaDetectVariantSites(
    const vector<MsaSeqInfo>& seqInfos,
    const vector<string>& alignedSeqs,
    uint32_t focalBegin, uint32_t focalEnd,
    vector<pair<uint32_t, uint32_t>>* noisyOut = nullptr)
{
    vector<MsaVariantSite> sites;
    if (focalEnd <= focalBegin) return sites;
    if (alignedSeqs.size() != seqInfos.size() || alignedSeqs.empty()) return sites;
    const size_t colCount = alignedSeqs.front().size();
    for (const auto& s : alignedSeqs)
        if (s.size() != colCount) return sites;
    if (colCount == 0) return sites;

    // Per-row span.
    vector<size_t> firstNG(alignedSeqs.size(), colCount);
    vector<size_t> lastNG(alignedSeqs.size(), 0);
    vector<uint8_t> hasB(alignedSeqs.size(), 0);
    for (size_t r = 0; r < alignedSeqs.size(); r++)
        for (size_t c = 0; c < colCount; c++)
            if (alignedSeqs[r][c] != '-') {
                firstNG[r] = min(firstNG[r], c);
                lastNG[r] = max(lastNG[r], c);
                hasB[r] = 1;
            }

    // Map columns to target offsets.
    vector<size_t> colByOff;
    colByOff.reserve(focalEnd - focalBegin);
    for (size_t c = 0; c < colCount; c++)
        if (alignedSeqs[0][c] != '-')
            colByOff.push_back(c);

    string focalTarget;
    focalTarget.reserve(colByOff.size());
    for (size_t c : colByOff) {
        char b = alignedSeqs[0][c];
        focalTarget.push_back(msaIsBase(b) ? b : 'N');
    }

    // Find dirty offsets.
    vector<uint8_t> isDirty(colByOff.size(), 0);
    for (size_t off = 0; off < colByOff.size(); off++) {
        size_t c = colByOff[off];
        char refBase = alignedSeqs[0][c];
        if (!msaIsBase(refBase)) continue;
        for (size_t r = 1; r < alignedSeqs.size(); r++) {
            if (!hasB[r] || c < firstNG[r] || c > lastNG[r]) continue;
            char b = alignedSeqs[r][c];
            if (b == '-' || (msaIsBase(b) && b != refBase)) { isDirty[off] = 1; break; }
        }
    }

    // Mark insertion gaps on preceding ref column.
    size_t gapStart = SIZE_MAX;
    for (size_t c = 0; c <= colCount; c++) {
        bool isFG = (c < colCount && alignedSeqs[0][c] == '-');
        if (isFG && gapStart == SIZE_MAX) gapStart = c;
        if ((!isFG || c == colCount) && gapStart != SIZE_MAX) {
            size_t ancCol = SIZE_MAX;
            for (size_t j = gapStart; j > 0; --j)
                if (alignedSeqs[0][j-1] != '-') { ancCol = j-1; break; }
            if (ancCol != SIZE_MAX) {
                auto it = lower_bound(colByOff.begin(), colByOff.end(), ancCol);
                if (it != colByOff.end() && *it == ancCol) {
                    size_t off = size_t(it - colByOff.begin());
                    for (size_t r = 1; r < alignedSeqs.size(); r++) {
                        if (!hasB[r] || ancCol < firstNG[r] || ancCol > lastNG[r]) continue;
                        for (size_t j = gapStart; j < c; j++)
                            if (msaIsBase(alignedSeqs[r][j])) { isDirty[off] = 1; goto nextGR; }
                    }
                }
            }
            nextGR: gapStart = SIZE_MAX;
        }
    }

    // Detect per-row noisy regions and merge into chunk-level regions.
    if (noisyOut)
        *noisyOut = msaDetectNoisyRegions(
            seqInfos, alignedSeqs, colByOff, firstNG, lastNG, hasB, focalBegin);

    // Walk dirty runs → variant sites.
    using AKey = tuple<string, string, string>;
    for (size_t bo = 0; bo < isDirty.size(); ) {
        if (!isDirty[bo]) { ++bo; continue; }
        size_t eo = bo;
        while (eo + 1 < isDirty.size() && isDirty[eo + 1]) ++eo;
        size_t bc = colByOff[bo], ec = colByOff[eo];
        for (size_t c = ec + 1; c < colCount && alignedSeqs[0][c] == '-'; c++) ec = c;

        string ref;
        for (size_t c = bc; c <= ec; c++) { char b = alignedSeqs[0][c]; if (msaIsBase(b)) ref.push_back(b); }
        if (ref.empty()) { bo = eo + 1; continue; }

        vector<OrientedReadId> refReads;
        map<AKey, vector<OrientedReadId>> altMap;
        for (size_t r = 0; r < alignedSeqs.size(); r++) {
            if (!hasB[r] || bc < firstNG[r] || ec > lastNG[r]) continue;
            string allele;
            for (size_t c = bc; c <= ec; c++) { char b = alignedSeqs[r][c]; if (msaIsBase(b)) allele.push_back(b); }
            if (allele == ref) refReads.push_back(seqInfos[r].oid);
            else {
                string tp = msaAlleleType(ref, allele);
                AKey key{tp, ref, allele};
                auto& v = altMap[key];
                if (find(v.begin(), v.end(), seqInfos[r].oid) == v.end()) v.push_back(seqInfos[r].oid);
            }
        }

        uint64_t maxAS = 0;
        vector<MsaAltAllele> reportable;
        for (const auto& [key, aReads] : altMap) {
            const auto& [tp, rK, alt] = key; (void)rK;
            if (tp != "SNP" && alt.size() < 2) continue;
            if (tp == "SNP" && msaSnpTouchesHomopolymer(focalTarget, uint32_t(bo), alt[0])) continue;
            reportable.push_back(MsaAltAllele{tp, alt, aReads});
            maxAS = max(maxAS, uint64_t(aReads.size()));
        }
        if (refReads.size() >= msaMinSnpRefSupport && maxAS >= msaMinSnpAltSupport) {
            uint32_t pos = focalBegin + uint32_t(bo);
            uint32_t tot = uint32_t(refReads.size() + maxAS);
            double af = tot > 0 ? double(maxAS) / double(tot) : 0.0;
            sites.push_back(MsaVariantSite{pos, ref, move(reportable), move(refReads), tot, uint32_t(maxAS), af});
        }
        bo = eo + 1;
    }

    // Second pass: detect het sites in insertion columns (where backbone has '-').
    // The backbone doesn't participate here — we look for disagreements among
    // non-backbone reads. The consensus allele among spanning reads serves as "ref".
    {
        size_t insStart = SIZE_MAX;
        for (size_t c = 0; c <= colCount; c++) {
            bool isIns = (c < colCount && alignedSeqs[0][c] == '-');
            if (isIns && insStart == SIZE_MAX) insStart = c;
            if ((!isIns || c == colCount) && insStart != SIZE_MAX) {
                size_t insEnd = c; // exclusive

                // Find the backbone position and base to anchor this insertion site.
                // Use the preceding ref column's backbone offset and base.
                uint32_t anchorPos = focalBegin; // fallback
                char anchorBase = 'N';
                for (size_t j = insStart; j > 0; --j) {
                    if (alignedSeqs[0][j - 1] != '-') {
                        anchorBase = alignedSeqs[0][j - 1];
                        auto it = lower_bound(colByOff.begin(), colByOff.end(), j - 1);
                        if (it != colByOff.end() && *it == j - 1)
                            anchorPos = focalBegin + uint32_t(it - colByOff.begin());
                        break;
                    }
                }

                // Extract alleles from non-backbone reads spanning this insertion.
                map<string, vector<OrientedReadId>> alleleMap;
                for (size_t r = 1; r < alignedSeqs.size(); r++) {
                    if (!hasB[r]) continue;
                    // Check if read spans the insertion region.
                    if (insStart < firstNG[r] || insEnd - 1 > lastNG[r]) continue;
                    string allele;
                    for (size_t j = insStart; j < insEnd; j++) {
                        char b = alignedSeqs[r][j];
                        if (msaIsBase(b)) allele.push_back(b);
                    }
                    // Empty allele = read has all gaps in insertion = no insertion.
                    // Use "-" to represent no-insertion.
                    if (allele.empty()) allele = "-";
                    alleleMap[allele].push_back(seqInfos[r].oid);
                }

                if (alleleMap.size() < 2) { insStart = SIZE_MAX; continue; }

                // Find consensus (most common allele).
                string consensus;
                size_t maxCount = 0;
                for (const auto& [allele, reads_vec] : alleleMap) {
                    if (reads_vec.size() > maxCount) {
                        maxCount = reads_vec.size();
                        consensus = allele;
                    }
                }

                // Build ref/alt from consensus, anchored with the preceding
                // backbone base (VCF-style).  This makes ref/alt comparable
                // strings so prefix/suffix trimming works correctly downstream.
                string anchorStr(1, anchorBase);
                string refStr = (consensus == "-")
                    ? anchorStr
                    : anchorStr + consensus;

                vector<OrientedReadId> refReads2 = alleleMap[consensus];
                vector<MsaAltAllele> reportable2;
                uint64_t maxAS2 = 0;
                for (const auto& [allele, reads_vec] : alleleMap) {
                    if (allele == consensus) continue;
                    if (reads_vec.size() < msaMinSnpAltSupport) continue;
                    // Anchor the alt allele with the backbone base.
                    string anchoredAlt = (allele == "-")
                        ? anchorStr
                        : anchorStr + allele;
                    string tp = msaAlleleType(refStr, anchoredAlt);
                    reportable2.push_back(MsaAltAllele{tp, anchoredAlt, reads_vec});
                    maxAS2 = max(maxAS2, uint64_t(reads_vec.size()));
                }

                if (refReads2.size() >= msaMinSnpRefSupport && maxAS2 >= msaMinSnpAltSupport) {
                    uint32_t tot = uint32_t(refReads2.size() + maxAS2);
                    double af = tot > 0 ? double(maxAS2) / double(tot) : 0.0;
                    MsaVariantSite vs{anchorPos, refStr, move(reportable2),
                        move(refReads2), tot, uint32_t(maxAS2), af, true};
                    sites.push_back(move(vs));
                }

                // Per-column SNP scan within the insertion block.
                // Look at each individual column for single-base disagreements
                // among reads that have bases there.
                // Skip single-column blocks — already handled by block-level scan above.
                if (insEnd - insStart < 2) { insStart = SIZE_MAX; continue; }
                for (size_t j = insStart; j < insEnd; j++) {
                    // Count bases at this column among non-backbone reads.
                    uint32_t baseCounts[4] = {0, 0, 0, 0}; // A C G T
                    vector<OrientedReadId> baseReads[4];
                    for (size_t r = 1; r < alignedSeqs.size(); r++) {
                        if (!hasB[r] || j < firstNG[r] || j > lastNG[r]) continue;
                        char b = alignedSeqs[r][j];
                        int bi = -1;
                        if (b == 'A' || b == 'a') bi = 0;
                        else if (b == 'C' || b == 'c') bi = 1;
                        else if (b == 'G' || b == 'g') bi = 2;
                        else if (b == 'T' || b == 't') bi = 3;
                        if (bi >= 0) {
                            baseCounts[bi]++;
                            baseReads[bi].push_back(seqInfos[r].oid);
                        }
                    }
                    // Find consensus base and check for alt bases.
                    int consBase = -1;
                    uint32_t consCount = 0;
                    for (int bi = 0; bi < 4; bi++) {
                        if (baseCounts[bi] > consCount) {
                            consCount = baseCounts[bi];
                            consBase = bi;
                        }
                    }
                    if (consBase < 0 || consCount < msaMinSnpRefSupport) continue;

                    vector<MsaAltAllele> colAlts;
                    uint64_t colMaxAS = 0;
                    const char baseChars[] = "ACGT";
                    for (int bi = 0; bi < 4; bi++) {
                        if (bi == consBase) continue;
                        if (baseCounts[bi] < msaMinSnpAltSupport) continue;
                        string altStr(1, baseChars[bi]);
                        colAlts.push_back(MsaAltAllele{"SNP", altStr, baseReads[bi]});
                        colMaxAS = max(colMaxAS, uint64_t(baseCounts[bi]));
                    }
                    if (colAlts.empty()) continue;

                    uint32_t colTot = uint32_t(consCount + colMaxAS);
                    double colAf = colTot > 0 ? double(colMaxAS) / double(colTot) : 0.0;
                    string colRef(1, baseChars[consBase]);
                    MsaVariantSite vs{anchorPos, colRef, move(colAlts),
                        baseReads[consBase], colTot, uint32_t(colMaxAS), colAf, true};
                    sites.push_back(move(vs));
                }

                insStart = SIZE_MAX;
            }
        }
    }

    return sites;
}

} // anonymous namespace

// ============================================================================
// Process one anchor window using multi-segment Theseus MSA.
// Builds one POA graph spanning all backbone segments, then aligns each
// read via align_from between its shared anchor nodes. Any read sharing
// >=2 anchors with the backbone participates (not just B-reads).
// ============================================================================
static void msaProcessWindow(
    Assembler& assembler,
    const AnchorWindow& window,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    MsaPhasingCounters& counters)
{
    const Reads& reads = assembler.getReads();
    const auto& mkrs = *assembler.markers;
    const uint64_t k = assembler.assemblerInfo->k;
    const OrientedReadId bbOid = window.backboneOrientedReadId;
    const auto bbJ = journeys[bbOid];
    const uint32_t nBbAnchors = window.backboneEnd - window.backboneBegin;
    if (nBbAnchors < 2) return;
    const uint32_t nSegments = nBbAnchors - 1;

    // Step 1: Build backbone segments.
    const uint32_t kh = uint32_t(k / 2);
    vector<string> segStrings;
    segStrings.reserve(nSegments);
    // Track backbone positions for variant site detection.
    vector<uint32_t> segFocalBegin(nSegments), segFocalEnd(nSegments);

    for (uint32_t i = 0; i < nSegments; i++) {
        uint32_t jp0 = window.backboneBegin + i;
        uint32_t jp1 = window.backboneBegin + i + 1;
        Shasta2AnchorId aL = bbJ[jp0], aR = bbJ[jp1];
        uint32_t oL = anchors.getOrdinal(aL, bbOid);
        uint32_t oR = anchors.getOrdinal(aR, bbOid);
        if (oL == invalid<uint32_t> || oR == invalid<uint32_t>) return;
        string seg = msaExtractSegment(reads, mkrs, k, bbOid, oL, oR);
        if (seg.empty()) return;
        segFocalBegin[i] = mkrs[bbOid.getValue()][oL].position + kh;
        segFocalEnd[i] = mkrs[bbOid.getValue()][oR].position + kh;
        segStrings.push_back(move(seg));
    }

    // Check total backbone span isn't pathological.
    size_t totalBbBases = 0;
    for (const auto& s : segStrings) totalBbBases += s.size();
    if (totalBbBases > msaMaxFocalSeqLen * nSegments) {
        counters.skippedLongFocal++;
        return;
    }

    // Build multi-segment POA graph from backbone.
    vector<string_view> segViews;
    segViews.reserve(nSegments);
    for (const auto& s : segStrings) segViews.push_back(s);

    theseus::Penalties pen(0, 2, 3, 1);
    theseus::Heuristics heur(false, false);
    vector<theseus::Graph::NodeId> nodeIds;
    theseus::TheseusMSA aligner(pen, heur, segViews, nodeIds, 1);

    // Step 2: Find reads sharing anchors with the backbone.
    // Only keep anchor hits that fall within the alignment chain between
    // the read and the backbone. Hits outside the chain are spurious
    // k-mer matches that would corrupt the POA graph.
    struct BoundaryHit { uint32_t boundaryIndex; uint32_t ordinal; };
    unordered_map<uint64_t, vector<BoundaryHit>> readHits;

    // Build read→alignment chain ordinal range lookup.
    // For each direct overlap with the backbone, record the partner read's
    // ordinal range from the alignment chain (matching prototype logic).
    struct ChainRange { uint32_t firstOrd; uint32_t lastOrd; };
    unordered_map<uint64_t, ChainRange> readChainRange;
    {
        const ReadId bbReadId = bbOid.getReadId();
        const auto& at = assembler.getAlignmentTable();
        for (uint32_t idx : at[bbOid.getValue()]) {
            const auto& ad = assembler.alignmentData[idx];
            OrientedReadId partner = ad.getOther(bbOid);
            ReadId partnerReadId = partner.getReadId();

            // Find which index (0 or 1) is the partner read.
            int targetIdx = (ad.readIds[0] == partnerReadId) ? 0 : 1;
            uint32_t firstOrd = ad.info.data[targetIdx].firstOrdinal;
            uint32_t lastOrd  = ad.info.data[targetIdx].lastOrdinal;

            // Flip ordinals if strand mismatch.
            Strand storedStrand = (targetIdx == 0) ? 0
                : (ad.isSameStrand ? Strand(0) : Strand(1));
            if (storedStrand != partner.getStrand()) {
                uint32_t mc = ad.info.data[targetIdx].markerCount;
                uint32_t f = mc - 1 - lastOrd;
                uint32_t l = mc - 1 - firstOrd;
                firstOrd = f;
                lastOrd = l;
            }

            auto it = readChainRange.find(partner.getValue());
            if (it == readChainRange.end() || (lastOrd - firstOrd) > (it->second.lastOrd - it->second.firstOrd)) {
                readChainRange[partner.getValue()] = ChainRange{firstOrd, lastOrd};
            }
        }
    }

    // Collect all anchor hits first (unfiltered).
    for (uint32_t bi = 0; bi <= nSegments; bi++) {
        uint32_t jp = window.backboneBegin + bi;
        if (jp >= bbJ.size()) break;
        Shasta2AnchorId aid = bbJ[jp];
        const auto anc = anchors[aid];
        for (const auto& info : anc) {
            if (info.orientedReadId == bbOid) continue;
            readHits[info.orientedReadId.getValue()].push_back({bi, info.ordinal});
        }
    }

    // Trim hits to chain range: remove hits outside the alignment chain
    // but keep the read if it still has >=2 valid hits.
    for (auto& [rv, hits] : readHits) {
        auto cit = readChainRange.find(rv);
        if (cit != readChainRange.end()) {
            const auto& cr = cit->second;
            hits.erase(
                remove_if(hits.begin(), hits.end(),
                    [&](const BoundaryHit& h) {
                        return h.ordinal < cr.firstOrd || h.ordinal > cr.lastOrd;
                    }),
                hits.end());
        }
    }

    // Sort hits, remove reads with <2 hits.
    for (auto it = readHits.begin(); it != readHits.end(); ) {
        auto& hits = it->second;
        sort(hits.begin(), hits.end(),
            [](const BoundaryHit& a, const BoundaryHit& b) {
                return a.boundaryIndex < b.boundaryIndex; });
        if (hits.size() < 2) it = readHits.erase(it);
        else ++it;
    }

    if (readHits.empty()) { counters.windowsProcessed++; return; }

    // Step 3: Align each read's sub-segments via align_from.
    // Track which reads were successfully aligned for variant detection.
    int readSeqId = 1; // 0 = backbone
    vector<OrientedReadId> alignedReadIds;
    alignedReadIds.push_back(bbOid); // index 0 = backbone

    // Sort reads by span (longest first) for better POA graph quality.
    vector<pair<uint32_t, uint64_t>> readsBySpan;
    readsBySpan.reserve(readHits.size());
    for (const auto& [rv, hits] : readHits) {
        const auto rm = mkrs[rv];
        uint32_t fp = rm[hits.front().ordinal].position;
        uint32_t lp = rm[hits.back().ordinal].position;
        readsBySpan.push_back({(lp > fp) ? (lp - fp) : 0, rv});
    }
    sort(readsBySpan.begin(), readsBySpan.end(),
        [](const auto& a, const auto& b) { return a.first > b.first; });

    for (const auto& [span, rv] : readsBySpan) {
        const auto& hits = readHits[rv];
        OrientedReadId oid = OrientedReadId::fromValue(ReadId(rv));
        bool aligned = false;
        uint32_t segCount = 0;
        uint32_t totalSegBases = 0;
        uint32_t skippedShort = 0;

        // Merge consecutive short inter-anchor spans: walk forward from
        // each anchor hit, skipping intermediate anchors until the
        // accumulated span reaches msaMinSegmentBases.  If no span
        // reaches the threshold, align the longest available span so
        // that short reads still contribute to the MSA.
        size_t hi = 0;
        while (hi + 1 < hits.size()) {
            uint32_t startBI  = hits[hi].boundaryIndex;
            uint32_t startOrd = hits[hi].ordinal;
            if (startBI >= nodeIds.size()) { hi++; continue; }

            // Scan forward for the first anchor that gives a span
            // >= msaMinSegmentBases.  Track the farthest valid anchor
            // in case none reaches the threshold.
            size_t best = hi;
            bool reached = false;
            const auto rm = mkrs[rv];
            uint32_t kh2 = uint32_t(k / 2);
            uint32_t beg = rm[startOrd].position + kh2;
            for (size_t hj = hi + 1; hj < hits.size(); hj++) {
                if (hits[hj].boundaryIndex <= startBI) continue;
                if (hits[hj].ordinal <= startOrd) continue;
                best = hj;
                uint32_t end = rm[hits[hj].ordinal].position + kh2;
                if (end > beg && (end - beg) >= msaMinSegmentBases) {
                    reached = true;
                    break;
                }
            }
            if (best == hi) { hi++; continue; }

            // If we didn't reach the threshold, best is the farthest
            // valid anchor — align whatever we have.
            uint32_t endBI  = hits[best].boundaryIndex;
            uint32_t endOrd = hits[best].ordinal;

            string seg = msaExtractSegment(reads, mkrs, k, oid, startOrd, endOrd);
            if (seg.empty()) { hi = best; continue; }

            if (!reached) skippedShort++;

            int endNode = (endBI < nodeIds.size())
                ? static_cast<int>(nodeIds[endBI]) : -1;
            aligner.align_from(seg, nodeIds[startBI], 1, true, 0, endNode, readSeqId);
            aligned = true;
            segCount++;
            totalSegBases += uint32_t(seg.size());
            hi = best;
        }
        if (aligned) {
            alignedReadIds.push_back(oid);
            readSeqId++;
        }
        cout << "    read " << oid << " hits=" << hits.size()
             << " segs=" << segCount << " bases=" << totalSegBases
             << " skippedShort=" << skippedShort
             << (aligned ? " ALIGNED" : " DROPPED") << endl;
    }

    if (alignedReadIds.size() < msaMinAnchorCoverage) {
        counters.skippedLowCov++;
        counters.windowsProcessed++;
        return;
    }

    // Step 4: Get MSA matrix and detect variant sites.
    auto msaMat = aligner.get_msa_matrix(readSeqId);

    // Convert MSAMatrix to vector<string> for variant detection.
    vector<string> alignedSeqs(msaMat.n_rows);
    for (int r = 0; r < msaMat.n_rows; r++) {
        alignedSeqs[r].resize(msaMat.n_cols);
        for (int c = 0; c < msaMat.n_cols; c++)
            alignedSeqs[r][c] = static_cast<char>(msaMat(r, c));
    }

    // Build MsaSeqInfo for the full MSA (one entry per aligned read).
    vector<MsaSeqInfo> fullSeqInfos;
    fullSeqInfos.reserve(alignedReadIds.size());
    for (const auto& oid : alignedReadIds)
        fullSeqInfos.push_back(MsaSeqInfo{oid, {}, 0, 0, true, 'B'});

    // get_msa_matrix returns num_sequences+1 rows (backbone + reads).
    // Truncate if sizes don't match.
    if (alignedSeqs.size() < fullSeqInfos.size())
        fullSeqInfos.resize(alignedSeqs.size());
    else if (alignedSeqs.size() > fullSeqInfos.size())
        alignedSeqs.resize(fullSeqInfos.size());

    // Write MSA matrix to a text file for inspection.
    {
        ostringstream fnss;
        fnss << "msa_window_" << bbOid << ".txt";
        const string msaFileName = fnss.str();
        ofstream msaFile(msaFileName);
        if (msaFile) {
            const int nRows = int(alignedSeqs.size());
            const int nCols = nRows > 0 ? int(alignedSeqs[0].size()) : 0;
            msaFile << "# MSA matrix: " << nRows << " rows x " << nCols << " cols\n";
            msaFile << "# Row 0 = backbone (" << bbOid << ")\n";
            for (size_t ri = 1; ri < alignedReadIds.size() && ri < alignedSeqs.size(); ri++)
                msaFile << "# Row " << ri << " = " << alignedReadIds[ri] << "\n";
            for (const auto& row : alignedSeqs)
                msaFile << row << "\n";
            cout << "  MSA written to " << msaFileName
                 << " (" << nRows << " rows, " << nCols << " cols)" << endl;
        }
    }

    // Detect variant sites across the full MSA.
    uint32_t fullFocalBegin = segFocalBegin[0];
    uint32_t fullFocalEnd = segFocalEnd[nSegments - 1];

    auto allSites = msaDetectVariantSites(fullSeqInfos, alignedSeqs, fullFocalBegin, fullFocalEnd);
    counters.sitesDetected += allSites.size();
    counters.pairsProcessed += nSegments;

    // Count unique reads across all backbone anchors (before any filtering).
    unordered_set<uint64_t> allAnchorReads;
    for (uint32_t bi = 0; bi <= nSegments; bi++) {
        uint32_t jp = window.backboneBegin + bi;
        if (jp >= bbJ.size()) break;
        const auto anc = anchors[bbJ[jp]];
        for (const auto& info : anc)
            if (info.orientedReadId != bbOid)
                allAnchorReads.insert(info.orientedReadId.getValue());
    }

    // Detailed output for this window.
    cout << "  MSA DETAIL journeyLen=" << bbJ.size()
         << " windowRange=[" << window.backboneBegin << "," << window.backboneEnd << ")"
         << endl;
    cout << "  MSA DETAIL bb=" << bbOid
         << " segs=" << nSegments
         << " alignedReads=" << (alignedReadIds.size() - 1)
         << " msaRows=" << alignedSeqs.size()
         << " seqInfos=" << fullSeqInfos.size()
         << " msaCols=" << (alignedSeqs.empty() ? 0 : alignedSeqs[0].size())
         << " focalBegin=" << fullFocalBegin
         << " focalEnd=" << fullFocalEnd
         << " focalSpan=" << (fullFocalEnd - fullFocalBegin)
         << " variantSites=" << allSites.size()
         << " uniqueReadsInAnchors=" << allAnchorReads.size()
         << endl;
    // Check column consistency.
    if (!alignedSeqs.empty()) {
        bool colConsistent = true;
        for (const auto& s : alignedSeqs)
            if (s.size() != alignedSeqs[0].size()) { colConsistent = false; break; }
        cout << "    colConsistent=" << colConsistent << endl;
        // Count dirty columns (span-aware, matching variant detection).
        vector<size_t> fNG(alignedSeqs.size(), alignedSeqs[0].size());
        vector<size_t> lNG(alignedSeqs.size(), 0);
        for (size_t r = 0; r < alignedSeqs.size(); r++)
            for (size_t c = 0; c < alignedSeqs[0].size(); c++)
                if (alignedSeqs[r][c] != '-') {
                    fNG[r] = min(fNG[r], c);
                    lNG[r] = max(lNG[r], c);
                }
        size_t dirtyCount = 0;
        for (size_t c = 0; c < alignedSeqs[0].size(); c++) {
            char ref = alignedSeqs[0][c];
            if (!msaIsBase(ref)) continue;
            for (size_t r = 1; r < alignedSeqs.size(); r++) {
                if (c < fNG[r] || c > lNG[r]) continue; // outside read span
                char b = alignedSeqs[r][c];
                if (b == '-' || (msaIsBase(b) && b != ref)) { dirtyCount++; break; }
            }
        }
        cout << "    dirtyColumns=" << dirtyCount << endl;
        // For dirty columns, count how many reads span each one.
        vector<uint32_t> dirtyCovs;
        for (size_t c = 0; c < alignedSeqs[0].size(); c++) {
            char ref = alignedSeqs[0][c];
            if (!msaIsBase(ref)) continue;
            bool isDirty = false;
            uint32_t spanCount = 0;
            for (size_t r = 1; r < alignedSeqs.size(); r++) {
                if (c < fNG[r] || c > lNG[r]) continue;
                spanCount++;
                char b = alignedSeqs[r][c];
                if (b == '-' || (msaIsBase(b) && b != ref)) isDirty = true;
            }
            if (isDirty) dirtyCovs.push_back(spanCount);
        }
        if (!dirtyCovs.empty()) {
            sort(dirtyCovs.begin(), dirtyCovs.end());
            cout << "    dirtyCoverage: min=" << dirtyCovs.front()
                 << " median=" << dirtyCovs[dirtyCovs.size()/2]
                 << " max=" << dirtyCovs.back() << endl;
        }
        // Count rows with any non-gap content.
        size_t activeRows = 0;
        for (size_t r = 0; r < alignedSeqs.size(); r++) {
            bool hasContent = false;
            for (char c : alignedSeqs[r])
                if (c != '-') { hasContent = true; break; }
            if (hasContent) activeRows++;
        }
        cout << "    activeRows=" << activeRows << " (of " << alignedSeqs.size() << ")" << endl;
    }
    for (size_t i = 0; i < allSites.size(); i++) {
        const auto& s = allSites[i];
        cout << "    site[" << i << "] pos=" << s.backbonePosition
             << " ref=" << s.refAllele
             << " refCov=" << s.refReads.size()
             << " altCov=" << s.altCov
             << " AF=" << fixed << setprecision(3) << s.af << defaultfloat;
        for (const auto& alt : s.altAlleles)
            cout << " " << alt.type << ":" << alt.sequence << "(" << alt.reads.size() << ")";
        cout << endl;
        // For SNP sites, show which reads carry ref vs alt.
        if (!s.altAlleles.empty()) {
            bool hasSNP = false;
            for (const auto& alt : s.altAlleles)
                if (alt.type == "SNP") { hasSNP = true; break; }
            if (hasSNP) {
                cout << "      ref reads:";
                for (const auto& oid : s.refReads) cout << " " << oid;
                cout << endl;
                for (const auto& alt : s.altAlleles) {
                    if (alt.type != "SNP") continue;
                    cout << "      alt " << alt.sequence << " reads:";
                    for (const auto& oid : alt.reads) cout << " " << oid;
                    cout << endl;
                }
            }
        }
    }

    if (allSites.empty()) { counters.windowsProcessed++; return; }

    // Step 5: Populate KmScratchpad with MSA-derived candidates and profiles,
    // then run the shared k-means pipeline from phaseOverlapsKmeans.

    // Build read→alignmentId map for the backbone.
    unordered_map<uint64_t, uint32_t> readToAlnId;
    {
        const auto& at = assembler.getAlignmentTable();
        for (uint32_t idx : at[bbOid.getValue()]) {
            const auto& ad = assembler.alignmentData[idx];
            OrientedReadId partner = ad.getOther(bbOid);
            readToAlnId[partner.getValue()] = idx;
        }
    }

    // Build KmScratchpad.overlaps — one per aligned read (excluding backbone).
    // Map from OrientedReadId value → overlap index in scratchpad.
    KmScratchpad scratch;
    KmPhasingOptions opts;
    unordered_map<uint64_t, uint32_t> readToOvIdx;
    for (size_t ri = 1; ri < alignedReadIds.size(); ri++) {
        uint64_t rv = alignedReadIds[ri].getValue();
        auto it = readToAlnId.find(rv);
        if (it == readToAlnId.end()) continue;
        KmOverlap ov;
        ov.alignmentId = it->second;
        ov.targetReadId = alignedReadIds[ri].getReadId();
        ov.qs = 0; ov.qe = 0; // Not used for MSA-derived profiles.
        ov.ts = 0; ov.te = 0;
        ov.cigarOffset = 0; ov.cigarTokenCount = 0;
        ov.isRev = 0; ov.queryIsRead0 = 0;
        readToOvIdx[rv] = uint32_t(scratch.overlaps.size());
        scratch.overlaps.push_back(ov);
    }
    if (scratch.overlaps.empty()) { counters.windowsProcessed++; return; }

    // Build numeric backbone sequence for repeat/homopolymer detection.
    // Uses the same 0-3 encoding as the k-means path (A=0, C=1, G=2, T=3).
    // Must use oriented read bases since MSA positions are oriented-read positions.
    const uint32_t bbLen = uint32_t(reads.getRead(bbOid.getReadId()).baseCount);
    vector<uint8_t> bbSeqVec(bbLen);
    for (uint32_t i = 0; i < bbLen; i++)
        bbSeqVec[i] = reads.getOrientedReadBase(bbOid, i).value;
    const uint8_t* bbSeq = bbSeqVec.data();


    // Decompose multiallelic sites into biallelic pairs (pgphase style).
    // Each alt allele with sufficient read support becomes its own candidate.
    // Track (siteIdx, altIdx) per candidate for profile building.
    sort(allSites.begin(), allSites.end(),
        [](const MsaVariantSite& a, const MsaVariantSite& b) {
            return a.backbonePosition < b.backbonePosition; });

    struct CandOrigin { uint32_t siteIdx; uint32_t altIdx; };
    vector<CandOrigin> candOrigins;

    constexpr uint32_t msaHpRunLen = 3;          // min homopolymer run length for SNP filtering
    constexpr uint32_t msaNoisyFlank = 10;       // extend noisy seeds by this many bp

    uint32_t repeatFiltered = 0, strandBiasFiltered = 0;
    uint32_t indelProxFiltered = 0, hpFiltered = 0;
    for (uint32_t si = 0; si < uint32_t(allSites.size()); si++) {
        const auto& site = allSites[si];
        const int refCov = int(site.refReads.size());

        for (uint32_t ai = 0; ai < uint32_t(site.altAlleles.size()); ai++) {
            const auto& alt = site.altAlleles[ai];
            if (alt.reads.size() < msaMinSnpAltSupport) continue;

            KmCandidate cand;

            if (site.isInsertionColumn) {
                // Pass 2 (insertion-column) sites: the variant lives in
                // columns where the backbone has gaps.  key.pos stays at
                // the anchor position (the preceding backbone base) so
                // kmIsHomopolymer/kmIsRepeatRegion check the correct
                // flanking reference context — matching pgphase semantics.
                // The anchor base was prepended to both ref and alt during
                // site construction, so strip it to get the raw content.
                cand.key.pos = site.backbonePosition;
                const string& refA = site.refAllele;
                const string& altA = alt.sequence;
                // Both refA and altA start with the anchor base.
                const uint32_t rawRefLen = uint32_t(refA.size()) - 1; // consensus length
                const uint32_t rawAltLen = uint32_t(altA.size()) - 1; // alt content length

                if (alt.type == "SNP") {
                    cand.key.type = KmVarType::Snp;
                    // Per-column SNPs: alt.sequence is a single base (no anchor prefix).
                    cand.key.altBase = altA.empty() ? 0 :
                        (altA[0] == 'C' ? 1 : altA[0] == 'G' ? 2 :
                         altA[0] == 'T' ? 3 : 0);
                    cand.key.refLen = 1; cand.key.altLen = 1;
                } else if (rawRefLen == 0 && rawAltLen > 0) {
                    // Pure insertion (consensus was "-", reads have bases).
                    cand.key.type = KmVarType::Insertion;
                    cand.key.refLen = 0;
                    cand.key.altLen = uint16_t(rawAltLen);
                    cand.key.altSeq = altA.substr(1); // strip anchor base
                } else if (rawAltLen == 0 && rawRefLen > 0) {
                    // Pure deletion (read has "-", consensus has bases).
                    cand.key.type = KmVarType::Deletion;
                    cand.key.refLen = uint16_t(rawRefLen);
                    cand.key.altLen = 0;
                } else {
                    // Length difference or MNP within insertion columns.
                    if (rawRefLen > rawAltLen) {
                        cand.key.type = KmVarType::Deletion;
                        cand.key.refLen = uint16_t(rawRefLen - rawAltLen);
                        cand.key.altLen = 0;
                    } else if (rawAltLen > rawRefLen) {
                        cand.key.type = KmVarType::Insertion;
                        cand.key.refLen = 0;
                        cand.key.altLen = uint16_t(rawAltLen - rawRefLen);
                        cand.key.altSeq = altA.substr(1 + rawRefLen);
                    } else {
                        // Same length MNP — treat as deletion for filtering.
                        cand.key.type = KmVarType::Deletion;
                        cand.key.refLen = uint16_t(rawRefLen);
                        cand.key.altLen = 0;
                    }
                }
            } else {
                // Pass 1 (backbone dirty-run) sites: trim common prefix/suffix
                // between ref and alt to find the minimal indel motif and its
                // true position on the backbone.
                const string& refA = site.refAllele;
                const string& altA = alt.sequence;
                uint32_t pfx = 0;
                while (pfx < refA.size() && pfx < altA.size() && refA[pfx] == altA[pfx]) pfx++;
                uint32_t sfx = 0;
                while (sfx < refA.size() - pfx && sfx < altA.size() - pfx &&
                       refA[refA.size() - 1 - sfx] == altA[altA.size() - 1 - sfx]) sfx++;
                const uint32_t trimRefLen = uint32_t(refA.size()) - pfx - sfx;
                const uint32_t trimAltLen = uint32_t(altA.size()) - pfx - sfx;

                cand.key.pos = site.backbonePosition + pfx;

                if (alt.type == "SNP") {
                    cand.key.type = KmVarType::Snp;
                    cand.key.altBase = altA.empty() ? 0 :
                        (altA[0] == 'C' ? 1 : altA[0] == 'G' ? 2 :
                         altA[0] == 'T' ? 3 : 0);
                    cand.key.refLen = 1; cand.key.altLen = 1;
                } else if (alt.type == "INS") {
                    cand.key.type = KmVarType::Insertion;
                    cand.key.refLen = 0;
                    cand.key.altLen = uint16_t(trimAltLen);
                    cand.key.altSeq = altA.substr(pfx, trimAltLen);
                } else if (alt.type == "DEL") {
                    cand.key.type = KmVarType::Deletion;
                    cand.key.refLen = uint16_t(trimRefLen);
                    cand.key.altLen = 0;
                } else {
                    // MNP — treat as indel for filtering purposes.
                    cand.key.type = KmVarType::Deletion;
                    cand.key.refLen = uint16_t(trimRefLen);
                    cand.key.altLen = 0;
                }
            }

            const int altCov = int(alt.reads.size());
            const int totalCov = refCov + altCov;
            const double af = totalCov > 0 ? double(altCov) / double(totalCov) : 0.0;

            cand.totalCov = totalCov;
            cand.refCov = refCov;
            cand.altCov = altCov;
            cand.alleleFraction = af;
            cand.alleCovs = {refCov, altCov};
            cand.nUniqAlles = 2;

            // Count fwd/rev strand alt reads for strand bias.
            int fwdAlt = 0, revAlt = 0;
            for (const auto& oid : alt.reads) {
                if (oid.getStrand() == 0) fwdAlt++;
                else revAlt++;
            }
            cand.fwdAlt = fwdAlt;
            cand.revAlt = revAlt;

            // Classify per-pair.
            if (af < opts.minAf) {
                cand.category = KmVariantCategory::LowAlleleFraction;
                cand.categoryFlag = KM_NON_VAR;
            } else if (af > opts.maxAf) {
                cand.category = KmVariantCategory::CleanHom;
                cand.categoryFlag = KM_NON_VAR;
            } else if (cand.key.type == KmVarType::Insertion ||
                       cand.key.type == KmVarType::Deletion) {
                // Match k-means path: all small indels (< minSvLen) are
                // excluded from phasing — nanopore indel noise dominates.
                // Only large SVs (>= minSvLen) get CleanHetIndel.
                const int indelLen = (cand.key.type == KmVarType::Insertion)
                    ? int(cand.key.altLen) : int(cand.key.refLen);
                if (indelLen >= opts.minSvLen) {
                    cand.category = KmVariantCategory::CleanHetIndel;
                    cand.categoryFlag = KM_REP_HET_VAR;
                } else {
                    cand.category = KmVariantCategory::RepeatHetIndel;
                    cand.categoryFlag = KM_REP_HET_VAR;
                    cand.isHomopolymerIndel = true;
                }
                repeatFiltered++;
            } else {
                // SNP — apply additional filters.

                // 1. Strand bias: Fisher exact test on fwd/rev alt counts.
                const int expected = (fwdAlt + revAlt) / 2;
                if (expected > 0) {
                    const double p = kmFisherExactTwoTail(fwdAlt, revAlt, expected, expected);
                    if (p < opts.strandBiasPval) {
                        cand.category = KmVariantCategory::StrandBias;
                        cand.categoryFlag = KM_NON_VAR;
                        strandBiasFiltered++;
                        scratch.candidates.push_back(move(cand));
                        candOrigins.push_back({si, ai});
                        continue;
                    }
                }

                // 2. Homopolymer-adjacent SNP: filter if the SNP position is
                //    flanked by a homopolymer run of length >= msaHpRunLen.
                if (cand.key.pos < bbLen) {
                    uint8_t base = bbSeq[cand.key.pos];
                    // Check forward run.
                    uint32_t fwdRun = 0;
                    for (uint32_t p = cand.key.pos + 1; p < bbLen && bbSeq[p] == base; p++)
                        fwdRun++;
                    // Check backward run.
                    uint32_t bwdRun = 0;
                    for (int64_t p = int64_t(cand.key.pos) - 1; p >= 0 && bbSeq[p] == base; p--)
                        bwdRun++;
                    // Also check if the alt base forms a homopolymer with flanking bases.
                    uint8_t altBase = cand.key.altBase;
                    uint32_t altFwd = 0, altBwd = 0;
                    for (uint32_t p = cand.key.pos + 1; p < bbLen && bbSeq[p] == altBase; p++)
                        altFwd++;
                    for (int64_t p = int64_t(cand.key.pos) - 1; p >= 0 && bbSeq[p] == altBase; p--)
                        altBwd++;

                    if (fwdRun + bwdRun >= msaHpRunLen || altFwd + altBwd >= msaHpRunLen) {
                        cand.category = KmVariantCategory::NonVariant;
                        cand.categoryFlag = KM_NON_VAR;
                        hpFiltered++;
                        scratch.candidates.push_back(move(cand));
                        candOrigins.push_back({si, ai});
                        continue;
                    }
                }

                cand.category = KmVariantCategory::CleanHetSnp;
                cand.categoryFlag = KM_GERMLINE_CLEAN;
            }
            scratch.candidates.push_back(move(cand));
            candOrigins.push_back({si, ai});
        }
    }
    // Noisy seed containment pass (pgphase post_process_noisy_regs + apply_noisy_containment_filter).
    // RepeatHetIndel candidates seed noisy regions. Extend, merge, then filter
    // any CleanHetSnp fully contained within a noisy region.
    {
        vector<pair<uint32_t, uint32_t>> noisySeeds;
        for (const auto& c : scratch.candidates) {
            if (c.category == KmVariantCategory::RepeatHetIndel) {
                uint32_t pos = c.key.pos;
                uint32_t len = max(uint32_t(c.key.refLen), uint32_t(c.key.altLen));
                uint32_t start = (pos > msaNoisyFlank) ? pos - msaNoisyFlank : 0;
                uint32_t end = pos + len + msaNoisyFlank;
                noisySeeds.push_back({start, end});
            }
        }
        if (!noisySeeds.empty()) {
            // Sort and merge.
            sort(noisySeeds.begin(), noisySeeds.end());
            vector<pair<uint32_t, uint32_t>> merged;
            merged.push_back(noisySeeds[0]);
            for (size_t i = 1; i < noisySeeds.size(); i++) {
                if (noisySeeds[i].first <= merged.back().second)
                    merged.back().second = max(merged.back().second, noisySeeds[i].second);
                else
                    merged.push_back(noisySeeds[i]);
            }
            // Filter CleanHetSnp candidates contained in noisy seeds.
            for (auto& c : scratch.candidates) {
                if (c.category != KmVariantCategory::CleanHetSnp) continue;
                for (const auto& [s, e] : merged) {
                    if (c.key.pos >= s && c.key.pos < e) {
                        c.category = KmVariantCategory::NonVariant;
                        c.categoryFlag = KM_NON_VAR;
                        indelProxFiltered++;
                        break;
                    }
                }
            }
        }
    }

    // Print classification summary with per-category site indices.
    uint32_t cleanHet = 0, nLowAf = 0, nCleanHom = 0;
    vector<uint32_t> idxCleanHet, idxRepeat, idxLowAf, idxCleanHom,
                     idxStrandBias, idxIndelProx, idxHpAdj;
    for (uint32_t ci = 0; ci < uint32_t(scratch.candidates.size()); ci++) {
        const auto& c = scratch.candidates[ci];
        if (c.category == KmVariantCategory::CleanHetSnp ||
            c.category == KmVariantCategory::CleanHetIndel) {
            cleanHet++; idxCleanHet.push_back(ci);
        } else if (c.category == KmVariantCategory::RepeatHetIndel) {
            idxRepeat.push_back(ci);
        } else if (c.category == KmVariantCategory::LowAlleleFraction) {
            nLowAf++; idxLowAf.push_back(ci);
        } else if (c.category == KmVariantCategory::CleanHom) {
            nCleanHom++; idxCleanHom.push_back(ci);
        } else if (c.category == KmVariantCategory::StrandBias) {
            idxStrandBias.push_back(ci);
        } else if (c.category == KmVariantCategory::NonVariant) {
            // Could be indelProx or hpAdj — distinguish by checking position.
            // For simplicity, just collect them all.
            idxIndelProx.push_back(ci);
        }
    }
    auto printCandList = [&](const char* label, const vector<uint32_t>& idxs) {
        if (idxs.empty()) return;
        cout << "      " << label << ": ";
        for (size_t j = 0; j < idxs.size(); j++) {
            if (j) cout << ", ";
            const auto& c = scratch.candidates[idxs[j]];
            const auto& o = candOrigins[idxs[j]];
            const auto& alt = allSites[o.siteIdx].altAlleles[o.altIdx];
            cout << "cand[" << idxs[j] << "] pos=" << c.key.pos
                 << " " << alt.type << ":" << alt.sequence
                 << " AF=" << fixed << setprecision(3) << c.alleleFraction;
        }
        cout << endl;
    };
    cout << "    candidates=" << scratch.candidates.size()
         << " (from " << allSites.size() << " sites)"
         << " cleanHet=" << cleanHet
         << " repeatIndel=" << repeatFiltered
         << " strandBias=" << strandBiasFiltered
         << " indelProx=" << indelProxFiltered
         << " hpAdj=" << hpFiltered
         << " lowAF=" << nLowAf
         << " cleanHom=" << nCleanHom << endl;
    printCandList("cleanHet", idxCleanHet);
    printCandList("repeatIndel", idxRepeat);
    printCandList("strandBias", idxStrandBias);
    printCandList("indelProx/hpAdj", idxIndelProx);
    printCandList("lowAF", idxLowAf);
    printCandList("cleanHom", idxCleanHom);
    counters.hetSitesUsed += cleanHet;
    if (cleanHet == 0) { counters.windowsProcessed++; return; }

    // Build KmScratchpad.overlapProfiles — per-overlap allele at each candidate.
    // We bypass kmBuildOverlapProfiles (which needs digars) and populate directly
    // from the MSA variant site read lists.
    const uint32_t numCand = uint32_t(scratch.candidates.size());
    const uint32_t numOv = uint32_t(scratch.overlaps.size());
    scratch.overlapProfiles.resize(numOv);
    for (uint32_t oi = 0; oi < numOv; oi++) {
        scratch.overlapProfiles[oi].overlapIdx = oi;
        scratch.overlapProfiles[oi].startVarIdx = -1;
        scratch.overlapProfiles[oi].endVarIdx = -1;
    }

    // For each candidate (biallelic pair), find which overlaps carry ref/alt.
    // Ref reads at the site get allele=0 for this candidate.
    // Only the alt reads matching this candidate's specific alt allele get allele=1.
    // Reads carrying a different alt at the same site get allele=0 (they're ref
    // with respect to this particular biallelic pair).
    for (uint32_t ci = 0; ci < numCand; ci++) {
        if (scratch.candidates[ci].categoryFlag == KM_NON_VAR) continue;
        const auto& origin = candOrigins[ci];
        const auto& site = allSites[origin.siteIdx];
        const auto& matchingAlt = site.altAlleles[origin.altIdx];

        // Build set of reads carrying the matching alt for fast lookup.
        unordered_set<uint64_t> matchingAltReads;
        for (const auto& oid : matchingAlt.reads)
            matchingAltReads.insert(oid.getValue());

        // Ref reads → allele 0.
        for (const auto& oid : site.refReads) {
            auto it = readToOvIdx.find(oid.getValue());
            if (it == readToOvIdx.end()) continue;
            auto& prof = scratch.overlapProfiles[it->second];
            if (prof.startVarIdx < 0) prof.startVarIdx = int(ci);
            while (prof.startVarIdx + int(prof.alleles.size()) < int(ci))
                prof.alleles.push_back(-1);
            prof.endVarIdx = int(ci);
            prof.alleles.push_back(0);
        }
        // All alt reads at this site.
        for (const auto& alt : site.altAlleles) {
            for (const auto& oid : alt.reads) {
                auto it = readToOvIdx.find(oid.getValue());
                if (it == readToOvIdx.end()) continue;
                auto& prof = scratch.overlapProfiles[it->second];
                if (prof.startVarIdx < 0) prof.startVarIdx = int(ci);
                while (prof.startVarIdx + int(prof.alleles.size()) < int(ci))
                    prof.alleles.push_back(-1);
                prof.endVarIdx = int(ci);
                // Matching alt → allele 1; other alts → allele 0 (ref for this pair).
                prof.alleles.push_back(matchingAltReads.count(oid.getValue()) ? 1 : 0);
            }
        }
    }

    // Step 6: Run shared k-means (round 1).
    kmRunKmeans(scratch, opts, KM_GERMLINE_CLEAN);

    // Step 7: Write initial results.
    const ReadId bbRid = bbOid.getReadId();
    kmWriteResults(assembler, bbRid, scratch);

    // Step 8: Iterative refinement of both cis and trans groups.
    // Each group is independently re-examined: recount alleles within the group,
    // reclassify, rebuild profiles, run k-means, peel off trans reads.
    // This can separate 3+ collapsed copies.

    // Assign each overlap to a group. Groups are numbered starting from 0.
    // After round 1: group 0 = cis (hap 1 or unassigned), group 1 = trans (hap 2).
    vector<int> groupId(numOv);
    for (uint32_t oi = 0; oi < numOv; oi++)
        groupId[oi] = (scratch.overlaps[oi].hap == 2) ? 1 : 0;
    int nextGroupId = 2;

    // Helper: refine a single group. Returns the number of reads peeled off.
    auto refineGroup = [&](int gid, const char* label) -> uint32_t {
        constexpr uint32_t maxRefineRounds = 10;
        uint32_t totalPeeled = 0;

        for (uint32_t round = 0; round < maxRefineRounds; round++) {
            // Count group members.
            uint32_t groupSize = 0;
            for (uint32_t oi = 0; oi < numOv; oi++)
                if (groupId[oi] == gid) groupSize++;
            if (groupSize < 6) break; // need at least 3 ref + 3 alt

            // Recount alleles using only this group's overlaps.
            for (auto& c : scratch.candidates) {
                c.totalCov = 1; c.refCov = 1; c.fwdRef = 1;
                c.altCov = 0; c.fwdAlt = 0; c.revAlt = 0; c.revRef = 0;
            }
            for (uint32_t ci = 0; ci < numCand; ci++) {
                const auto& origin = candOrigins[ci];
                const auto& site = allSites[origin.siteIdx];
                const auto& matchingAlt = site.altAlleles[origin.altIdx];
                unordered_set<uint64_t> matchAltSet;
                for (const auto& oid : matchingAlt.reads)
                    matchAltSet.insert(oid.getValue());

                auto& c = scratch.candidates[ci];
                for (const auto& oid : site.refReads) {
                    auto it = readToOvIdx.find(oid.getValue());
                    if (it == readToOvIdx.end() || groupId[it->second] != gid) continue;
                    c.totalCov++; c.refCov++;
                }
                for (const auto& alt : site.altAlleles) {
                    for (const auto& oid : alt.reads) {
                        auto it = readToOvIdx.find(oid.getValue());
                        if (it == readToOvIdx.end() || groupId[it->second] != gid) continue;
                        c.totalCov++;
                        if (matchAltSet.count(oid.getValue())) c.altCov++;
                        else c.refCov++;
                    }
                }
                c.alleleFraction = c.totalCov > 0 ? double(c.altCov) / double(c.totalCov) : 0.0;
            }

            // Reclassify with group-only counts.
            uint32_t roundCleanHet = 0;
            for (uint32_t ci = 0; ci < numCand; ci++) {
                auto& c = scratch.candidates[ci];
                if (c.alleleFraction < opts.minAf || c.alleleFraction > opts.maxAf) {
                    c.categoryFlag = KM_NON_VAR;
                } else if (c.key.type == KmVarType::Snp) {
                    c.category = KmVariantCategory::CleanHetSnp;
                    c.categoryFlag = KM_GERMLINE_CLEAN;
                    roundCleanHet++;
                } else {
                    c.categoryFlag = KM_REP_HET_VAR;
                }
            }

            cout << "    refine " << label << " round " << round
                 << ": groupSize=" << groupSize
                 << " cleanHet=" << roundCleanHet << endl;

            if (roundCleanHet == 0) break;

            // Rebuild profiles for this group's overlaps.
            scratch.overlapProfiles.clear();
            scratch.overlapProfiles.resize(numOv);
            for (uint32_t oi = 0; oi < numOv; oi++) {
                scratch.overlapProfiles[oi].overlapIdx = oi;
                scratch.overlapProfiles[oi].startVarIdx = -1;
                scratch.overlapProfiles[oi].endVarIdx = -1;
                scratch.overlapProfiles[oi].alleles.clear();
            }
            for (uint32_t ci = 0; ci < numCand; ci++) {
                if (scratch.candidates[ci].categoryFlag == KM_NON_VAR) continue;
                const auto& origin = candOrigins[ci];
                const auto& site = allSites[origin.siteIdx];
                const auto& matchingAlt = site.altAlleles[origin.altIdx];
                unordered_set<uint64_t> matchAltSet;
                for (const auto& oid : matchingAlt.reads)
                    matchAltSet.insert(oid.getValue());

                for (const auto& oid : site.refReads) {
                    auto it = readToOvIdx.find(oid.getValue());
                    if (it == readToOvIdx.end() || groupId[it->second] != gid) continue;
                    auto& prof = scratch.overlapProfiles[it->second];
                    if (prof.startVarIdx < 0) prof.startVarIdx = int(ci);
                    while (prof.startVarIdx + int(prof.alleles.size()) < int(ci))
                        prof.alleles.push_back(-1);
                    prof.endVarIdx = int(ci);
                    prof.alleles.push_back(0);
                }
                for (const auto& alt : site.altAlleles) {
                    for (const auto& oid : alt.reads) {
                        auto it = readToOvIdx.find(oid.getValue());
                        if (it == readToOvIdx.end() || groupId[it->second] != gid) continue;
                        auto& prof = scratch.overlapProfiles[it->second];
                        if (prof.startVarIdx < 0) prof.startVarIdx = int(ci);
                        while (prof.startVarIdx + int(prof.alleles.size()) < int(ci))
                            prof.alleles.push_back(-1);
                        prof.endVarIdx = int(ci);
                        prof.alleles.push_back(matchAltSet.count(oid.getValue()) ? 1 : 0);
                    }
                }
            }

            // Reset hap assignments and run k-means.
            for (auto& ov : scratch.overlaps) ov.hap = 0;
            kmRunKmeans(scratch, opts, KM_GERMLINE_CLEAN);

            // Peel off trans overlaps into a new group.
            uint32_t peeled = 0;
            int newGid = nextGroupId;
            for (uint32_t oi = 0; oi < numOv; oi++) {
                if (groupId[oi] == gid && scratch.overlaps[oi].hap == 2) {
                    groupId[oi] = newGid;
                    peeled++;
                }
            }
            if (peeled > 0) nextGroupId++;
            totalPeeled += peeled;

            cout << "    refine " << label << " round " << round
                 << ": peeled=" << peeled
                 << " -> group " << newGid << endl;

            if (peeled == 0) break;
        }
        return totalPeeled;
    };

    // Refine both initial groups.
    refineGroup(0, "cis");
    refineGroup(1, "trans");

    // Recursively refine any new groups that were created.
    // (Groups created by peeling may themselves contain multiple copies.)
    for (int gid = 2; gid < nextGroupId; gid++) {
        string label = "group" + to_string(gid);
        refineGroup(gid, label.c_str());
    }

    // Write final matchStates.
    // Group 0 = cis (matchState 1), group 1 = trans (matchState 2),
    // groups 2+ = additional copies (matchState 3+, capped at 255).
    for (uint32_t oi = 0; oi < numOv; oi++) {
        auto& ad = assembler.alignmentData[scratch.overlaps[oi].alignmentId];
        int g = groupId[oi];
        uint8_t matchState;
        if (g == 0) matchState = 1;       // cis
        else if (g == 1) matchState = 2;  // trans
        else matchState = 3;              // additional copy
        ad.setHifiasmEcMatchStateFromReadPerspective(bbRid, matchState);
    }

    // Count and log per-window results.
    // Count reads per group.
    map<int, uint32_t> groupCounts;
    for (uint32_t oi = 0; oi < numOv; oi++)
        groupCounts[groupId[oi]]++;

    uint32_t wCis = groupCounts.count(0) ? groupCounts[0] : 0;
    uint32_t wTrans = groupCounts.count(1) ? groupCounts[1] : 0;
    uint32_t wOther = 0;
    for (const auto& [g, cnt] : groupCounts)
        if (g >= 2) wOther += cnt;

    for (const auto& ov : scratch.overlaps) {
        OrientedReadId partner = assembler.alignmentData[ov.alignmentId].getOther(bbOid);
        int g = groupId[&ov - &scratch.overlaps[0]];
        cout << "    read=" << partner << " -> group" << g << endl;
        counters.readsPhased++;
        if (g == 0) counters.cisCount++;
        else if (g == 1) counters.transCount++;
    }
    cout << "  MSA window bb=" << bbOid
         << " segs=" << nSegments
         << " reads=" << scratch.overlaps.size()
         << " sites=" << allSites.size()
         << " het=" << cleanHet
         << " groups=" << nextGroupId
         << " group0(cis)=" << wCis
         << " group1(trans)=" << wTrans;
    if (wOther > 0) cout << " otherGroups=" << wOther;
    cout << endl;
    counters.windowsProcessed++;
}


// ============================================================================
// Public entry point
// ============================================================================
void Assembler::phaseOverlapsMSA(uint64_t threadCount)
{
    cout << timestamp << "=== MSA-Based Overlap Phasing ===" << endl;
    DINARA_ASSERT(shasta2Anchors);
    DINARA_ASSERT(shasta2Journeys);
    DINARA_ASSERT(shasta2Journeys->isOpen());
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();

    if (threadCount == 0)
        threadCount = std::thread::hardware_concurrency();
    threadCount = max<uint64_t>(1, threadCount);

    const uint64_t readCount = reads->readCount();
    cout << timestamp << "Read count: " << readCount << ", threads: " << threadCount << endl;
    if (readCount == 0) return;

    // Build readIdsSortedByLength for window planning.
    vector<ReadId> readIdsSorted(readCount);
    iota(readIdsSorted.begin(), readIdsSorted.end(), ReadId(0));
    sort(readIdsSorted.begin(), readIdsSorted.end(),
        [this](ReadId a, ReadId b) {
            return reads->getRead(a).baseCount > reads->getRead(b).baseCount;
        });

    // Plan anchor windows.
    cout << timestamp << "Planning anchor windows..." << endl;
    vector<AnchorWindow> anchorWindows;
    computeAnchorWindowsClean(
        shasta2Anchors, shasta2Journeys,
        readIdsSorted, anchorWindows, threadCount);
    cout << timestamp << "Planned " << anchorWindows.size() << " anchor windows." << endl;

    if (anchorWindows.empty()) return;

    // Process only the first (largest) window for now.
    MsaPhasingCounters counters;
    const auto begin = chrono::steady_clock::now();
    msaProcessWindow(*this, anchorWindows[0],
        *shasta2Anchors, *shasta2Journeys, counters);

    const auto end = chrono::steady_clock::now();
    double secs = chrono::duration<double>(end - begin).count();

    cout << timestamp << "MSA phasing complete."
         << " windows=" << counters.windowsProcessed.load()
         << " pairs=" << counters.pairsProcessed.load()
         << " sitesDetected=" << counters.sitesDetected.load()
         << " hetSitesUsed=" << counters.hetSitesUsed.load()
         << " readsPhased=" << counters.readsPhased.load()
         << " cis=" << counters.cisCount.load()
         << " trans=" << counters.transCount.load()
         << " skippedLowCov=" << counters.skippedLowCov.load()
         << " skippedLongFocal=" << counters.skippedLongFocal.load()
         << " seconds=" << fixed << setprecision(2) << secs
         << defaultfloat << endl;
}
