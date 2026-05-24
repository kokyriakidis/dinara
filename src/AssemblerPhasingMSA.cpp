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


// Detect variant sites from an MSA of one anchor pair.
static vector<MsaVariantSite> msaDetectVariantSites(
    const vector<MsaSeqInfo>& seqInfos,
    const vector<string>& alignedSeqs,
    uint32_t focalBegin, uint32_t focalEnd)
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
    // Raise the indel length gate for repeat/homopolymer detection.
    // pgphase uses 5, but nanopore homopolymer errors routinely produce
    // indels of 6-10+ bases in long runs. The MSA path lacks the CIGAR-based
    // noisy region detector that compensates in the k-means path.
    constexpr int msaNoisyRegMaxXgaps = 20;

    // Build KmScratchpad.candidates — one per variant site.
    // Sort by position to match k-means expectations.
    sort(allSites.begin(), allSites.end(),
        [](const MsaVariantSite& a, const MsaVariantSite& b) {
            return a.backbonePosition < b.backbonePosition; });

    uint32_t repeatFiltered = 0;
    for (const auto& site : allSites) {
        KmCandidate cand;
        cand.key.pos = site.backbonePosition;
        // Determine variant type from the highest-coverage alt allele.
        if (!site.altAlleles.empty()) {
            size_t bestIdx = 0;
            for (size_t ai = 1; ai < site.altAlleles.size(); ai++)
                if (site.altAlleles[ai].reads.size() > site.altAlleles[bestIdx].reads.size())
                    bestIdx = ai;
            const auto& bestAlt = site.altAlleles[bestIdx];
            if (bestAlt.type == "SNP") {
                cand.key.type = KmVarType::Snp;
                cand.key.altBase = bestAlt.sequence.empty() ? 0 :
                    (bestAlt.sequence[0] == 'C' ? 1 : bestAlt.sequence[0] == 'G' ? 2 :
                     bestAlt.sequence[0] == 'T' ? 3 : 0);
                cand.key.refLen = 1; cand.key.altLen = 1;
            } else if (bestAlt.type == "INS") {
                cand.key.type = KmVarType::Insertion;
                cand.key.refLen = 0;
                cand.key.altLen = uint16_t(bestAlt.sequence.size());
                cand.key.altSeq = bestAlt.sequence;
            } else {
                cand.key.type = KmVarType::Deletion;
                cand.key.refLen = uint16_t(bestAlt.sequence.size());
                cand.key.altLen = 0;
            }
        } else {
            cand.key.type = KmVarType::Snp;
            cand.key.refLen = 1; cand.key.altLen = 1;
        }
        cand.totalCov = int(site.totalCov);
        cand.refCov = int(site.refReads.size());
        cand.altCov = int(site.altCov);
        cand.alleleFraction = site.af;
        cand.alleCovs = {cand.refCov, cand.altCov};
        cand.nUniqAlles = 2;

        // Classify using the same logic as kmClassifyVariantInitial.
        if (site.af < opts.minAf) {
            cand.category = KmVariantCategory::LowAlleleFraction;
            cand.categoryFlag = KM_NON_VAR;
        } else if (site.af > opts.maxAf) {
            cand.category = KmVariantCategory::CleanHom;
            cand.categoryFlag = KM_NON_VAR;
        } else if (cand.key.type == KmVarType::Insertion ||
                   cand.key.type == KmVarType::Deletion) {
            // Filter all indels — nanopore indel noise dominates in MSA.
            // Only SNPs are used for phasing.
            cand.category = KmVariantCategory::RepeatHetIndel;
            cand.categoryFlag = KM_REP_HET_VAR;
            cand.isHomopolymerIndel = true;
            repeatFiltered++;
        } else {
            cand.category = KmVariantCategory::CleanHetSnp;
            cand.categoryFlag = KM_GERMLINE_CLEAN;
        }
        scratch.candidates.push_back(move(cand));
    }
    // Print classification summary with per-category site indices.
    uint32_t cleanHet = 0, nLowAf = 0, nCleanHom = 0;
    vector<uint32_t> idxCleanHet, idxRepeat, idxLowAf, idxCleanHom;
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
        }
    }
    auto printSiteList = [&](const char* label, const vector<uint32_t>& idxs) {
        if (idxs.empty()) return;
        cout << "      " << label << ": ";
        for (size_t j = 0; j < idxs.size(); j++) {
            if (j) cout << ", ";
            cout << "site[" << idxs[j] << "] pos=" << scratch.candidates[idxs[j]].key.pos;
        }
        cout << endl;
    };
    cout << "    classification: cleanHet=" << cleanHet
         << " repeatIndel=" << repeatFiltered
         << " lowAF=" << nLowAf
         << " cleanHom=" << nCleanHom << endl;
    printSiteList("cleanHet", idxCleanHet);
    printSiteList("repeatIndel", idxRepeat);
    printSiteList("lowAF", idxLowAf);
    printSiteList("cleanHom", idxCleanHom);
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

    // For each candidate, find which overlaps carry ref/alt.
    for (uint32_t ci = 0; ci < numCand; ci++) {
        if (scratch.candidates[ci].categoryFlag == KM_NON_VAR) continue;
        const auto& site = allSites[ci];

        // Ref reads.
        for (const auto& oid : site.refReads) {
            auto it = readToOvIdx.find(oid.getValue());
            if (it == readToOvIdx.end()) continue;
            auto& prof = scratch.overlapProfiles[it->second];
            if (prof.startVarIdx < 0) prof.startVarIdx = int(ci);
            while (prof.startVarIdx + int(prof.alleles.size()) < int(ci))
                prof.alleles.push_back(-1);
            prof.endVarIdx = int(ci);
            prof.alleles.push_back(0); // ref
        }
        // Alt reads.
        for (const auto& alt : site.altAlleles) {
            for (const auto& oid : alt.reads) {
                auto it = readToOvIdx.find(oid.getValue());
                if (it == readToOvIdx.end()) continue;
                auto& prof = scratch.overlapProfiles[it->second];
                if (prof.startVarIdx < 0) prof.startVarIdx = int(ci);
                while (prof.startVarIdx + int(prof.alleles.size()) < int(ci))
                    prof.alleles.push_back(-1);
                prof.endVarIdx = int(ci);
                prof.alleles.push_back(1); // alt
            }
        }
    }

    // Dummy digar/noisy region arrays (not used but kmRefineCis may reference).
    scratch.digarBegin.assign(numOv, 0);
    scratch.digarEnd.assign(numOv, 0);
    scratch.overlapNoisyBegin.assign(numOv, 0);
    scratch.overlapNoisyEnd.assign(numOv, 0);

    // Step 6: Run shared k-means.
    kmRunKmeans(scratch, opts, KM_GERMLINE_CLEAN);

    // Step 7: Write results.
    const ReadId bbRid = bbOid.getReadId();
    kmWriteResults(assembler, bbRid, scratch);

    // Step 8: Cis refinement.
    kmRefineCis(assembler, bbRid, scratch, opts, bbLen, false);

    // Count and log per-window results.
    uint32_t wCis = 0, wTrans = 0, wCisDiffCopy = 0, wUnclass = 0;
    for (const auto& ov : scratch.overlaps) {
        const auto& ad = assembler.alignmentData[ov.alignmentId];
        uint8_t ms = ad.getHifiasmEcMatchStateFromReadPerspective(bbRid);
        const char* label = (ms == 1) ? "cis" : (ms == 2) ? "trans" : (ms == 3) ? "cisDiffCopy" : "unclass";
        if (ms == 1) { counters.cisCount++; wCis++; }
        else if (ms == 2) { counters.transCount++; wTrans++; }
        else if (ms == 3) { wCisDiffCopy++; }
        else { wUnclass++; }
        counters.readsPhased++;
        // Find the OrientedReadId for this overlap.
        OrientedReadId partner = ad.getOther(bbOid);
        cout << "    read=" << partner << " -> " << label << endl;
    }
    cout << "  MSA window bb=" << bbOid
         << " segs=" << nSegments
         << " reads=" << scratch.overlaps.size()
         << " sites=" << allSites.size()
         << " het=" << cleanHet
         << " cis=" << wCis
         << " trans=" << wTrans
         << " cisDiffCopy=" << wCisDiffCopy
         << " unclass=" << wUnclass
         << endl;
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
