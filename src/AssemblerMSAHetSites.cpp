// AssemblerMSAHetSites.cpp
//
// Anchor-guided transitive graph MSA for het-site detection.
// For one focal oriented read:
//   1. Walk the focal linear journey and collect every non-focal read that
//      touches one or more focal anchors.
//   2. For each candidate read, keep the longest monotone chain of shared
//      anchors in focal-journey order and in that read's own journey order.
//   3. Extract the candidate subsequence between the chain's outer anchor
//      ordinals and add it to an abPOA graph seeded by the full focal read.
//   4. Write the resulting whole-read graph as GFA.
//   5. Run "vg deconstruct -p focal -a -R" to get snarls/alleles.
//   6. Parse the VCF into VariantEvent entries and store in
//      shasta2VariantEvents[readSlot].

#include "Assembler.hpp"
#include "Base.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"
#include "abpoa/abpoa.h"
// Theseus disabled: headers removed; this file does not use any theseus symbols.

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace dinara;
using namespace std;

// ---------------------------------------------------------------------------
// GFA graph (segment adjacency, used for path reconstruction).
// ---------------------------------------------------------------------------

namespace {

struct GFAGraph {
    map<string, string>         seq;    // segment name → sequence
    map<string, vector<string>> succ;   // successors
    map<string, vector<string>> pred;   // predecessors
    vector<string>              sources; // nodes with no predecessors
};

bool parseGFA(const string& filename, GFAGraph& g)
{
    ifstream in(filename);
    if (!in.is_open()) return false;
    string line;
    while (getline(in, line)) {
        if (line.size() < 2) continue;
        if (line[0] == 'S') {
            istringstream ss(line);
            string tag, name, sq;
            ss >> tag >> name >> sq;
            g.seq[name] = sq;
            g.succ[name];
            g.pred[name];
        } else if (line[0] == 'L') {
            istringstream ss(line);
            string tag, from, fstrand, to, tstrand, cigar;
            ss >> tag >> from >> fstrand >> to >> tstrand >> cigar;
            g.succ[from].push_back(to);
            g.pred[to].push_back(from);
        }
    }
    for (auto& [n, preds] : g.pred)
        if (preds.empty()) g.sources.push_back(n);
    return true;
}

// ---------------------------------------------------------------------------
// Path reconstruction via greedy sequence matching.
//
// Given the compact POA graph and a read sequence, recover the ordered list
// of node names the read traverses. At each branch point, pick the successor
// whose sequence matches the next characters of the read. Returns empty if
// reconstruction fails (read doesn't match any path — should not happen for
// spanning reads that were aligned by theseus).
// ---------------------------------------------------------------------------
vector<string> reconstructPath(const GFAGraph& g,
                                const string&   readSeq)
{
    // Find the unique source node.
    if (g.sources.size() != 1) return {};
    string cur = g.sources[0];

    vector<string> path;
    size_t pos = 0; // current position in readSeq

    set<string> visited;
    while (true) {
        if (!visited.insert(cur).second) return {}; // cycle guard

        const string& nodeSeq = g.seq.at(cur);
        // Check that readSeq matches at this position (allow up to 1 mismatch
        // for sequencing noise on long nodes).
        if (pos + nodeSeq.size() > readSeq.size()) return {};
        path.push_back(cur);
        pos += nodeSeq.size();

        const auto& nexts = g.succ.at(cur);
        if (nexts.empty()) break; // reached sink

        if (nexts.size() == 1) {
            cur = nexts[0];
            continue;
        }

        // Branch: pick the successor whose sequence matches readSeq[pos...].
        string chosen;
        for (const string& nxt : nexts) {
            const string& ns = g.seq.at(nxt);
            if (pos + ns.size() <= readSeq.size() &&
                readSeq.compare(pos, ns.size(), ns) == 0) {
                chosen = nxt;
                break;
            }
        }
        if (chosen.empty()) {
            // No exact match — pick by best prefix match (handles SNP mismatches
            // introduced by sequencing error vs the POA consensus).
            int bestScore = -1;
            for (const string& nxt : nexts) {
                const string& ns = g.seq.at(nxt);
                int matches = 0;
                for (size_t i = 0; i < ns.size() && pos + i < readSeq.size(); i++)
                    if (readSeq[pos + i] == ns[i]) matches++;
                if (matches > bestScore) { bestScore = matches; chosen = nxt; }
            }
        }
        if (chosen.empty()) return {};
        cur = chosen;
    }

    return path;
}

// ---------------------------------------------------------------------------
// Write GFA with embedded P lines (one per read).
// focal_path_name is always placed first so vg deconstruct uses it as ref.
// Returns true on success.
// ---------------------------------------------------------------------------
bool writeGFAWithPaths(
    const string& srcGfaPath,
    const string& dstGfaPath,
    const GFAGraph& g,
    const vector<string>& readSeqs,
    const vector<string>& pathNames)  // pathNames[0] == "focal"
{
    // Copy original S/L lines.
    ifstream src(srcGfaPath);
    if (!src.is_open()) return false;
    ofstream dst(dstGfaPath);
    if (!dst.is_open()) return false;

    string line;
    while (getline(src, line)) {
        if (line.empty()) continue;
        if (line[0] == 'S' || line[0] == 'L')
            dst << line << "\n";
    }

    // Append P lines.
    for (size_t i = 0; i < readSeqs.size(); i++) {
        const vector<string> path = reconstructPath(g, readSeqs[i]);
        if (path.empty()) continue;

        dst << "P\t" << pathNames[i] << "\t";
        for (size_t j = 0; j < path.size(); j++) {
            if (j) dst << ",";
            dst << path[j] << "+";
        }
        dst << "\t*\n";
    }
    return true;
}

// ---------------------------------------------------------------------------
// VCF record from vg deconstruct -n -R.
// ---------------------------------------------------------------------------
struct VcfRecord {
    uint32_t       pos;           // 1-based POS in focal read
    string         id;            // snarl ID (e.g. ">2>5")
    string         ref;
    vector<string> alts;          // may include "*" (star allele from -R)
    int            level;         // LV: 0=top-level
    string         parentSite;    // PS: parent snarl ID (empty if level==0)
    int            parentAllele;  // PA: which allele of parent the ref path is on
                                  //     0  = focal read passes through this nested site
                                  //     >0 = focal read does NOT enter this site
                                  //     -1 = not set (level==0)
    uint32_t       rsStart;       // RS: top-level site start on focal read
    uint32_t       rsEnd;         // RD: top-level site end   on focal read
    vector<int>    genotypes;     // one per sample; -1=missing(.), -2=star(*)
};

// Parse VCF text from vg deconstruct -n -R.
// sampleNamesOut is populated from the #CHROM header line.
vector<VcfRecord> parseVcf(const string& vcfText, vector<string>& sampleNamesOut)
{
    vector<VcfRecord> records;
    istringstream in(vcfText);
    string line;
    bool headerParsed = false;

    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line.size() >= 2 && line[0] == '#' && line[1] == '#') continue;
        if (line[0] == '#') {
            istringstream ss(line);
            vector<string> cols;
            string col;
            while (ss >> col) cols.push_back(col);
            for (size_t i = 9; i < cols.size(); i++)
                sampleNamesOut.push_back(cols[i]);
            headerParsed = true;
            continue;
        }
        if (!headerParsed) continue;

        istringstream ss(line);
        vector<string> cols;
        string col;
        while (ss >> col) cols.push_back(col);
        if (cols.size() < 9) continue;

        VcfRecord rec;
        rec.pos          = (uint32_t)stoul(cols[1]);
        rec.id           = cols[2];
        rec.ref          = cols[3];
        rec.level        = 0;
        rec.parentAllele = -1;
        rec.rsStart      = 0;
        rec.rsEnd        = 0;

        // ALT: comma-separated, may include "*"
        {
            istringstream altss(cols[4]);
            string a;
            while (getline(altss, a, ',')) rec.alts.push_back(a);
        }

        // INFO: parse LV, PS, PA, RS, RD
        {
            istringstream infoss(cols[7]);
            string token;
            while (getline(infoss, token, ';')) {
                if      (token.size() > 3 && token.substr(0,3) == "LV=")
                    rec.level = stoi(token.substr(3));
                else if (token.size() > 3 && token.substr(0,3) == "PS=")
                    rec.parentSite = token.substr(3);
                else if (token.size() > 3 && token.substr(0,3) == "PA=")
                    rec.parentAllele = stoi(token.substr(3));
                else if (token.size() > 3 && token.substr(0,3) == "RS=")
                    rec.rsStart = (uint32_t)stoul(token.substr(3));
                else if (token.size() > 3 && token.substr(0,3) == "RD=")
                    rec.rsEnd   = (uint32_t)stoul(token.substr(3));
            }
        }

        // Genotypes: FORMAT=GT (usually first field; find it explicitly).
        int gtIdx = 0;
        {
            istringstream fmtss(cols[8]);
            string f; int fi = 0;
            while (getline(fmtss, f, ':')) { if (f == "GT") { gtIdx = fi; break; } fi++; }
        }
        for (size_t i = 9; i < cols.size(); i++) {
            // Extract the GT sub-field at gtIdx position within colon-separated FORMAT.
            string gtField = cols[i];
            if (gtIdx > 0) {
                istringstream fs(cols[i]);
                string sub; int fi = 0;
                while (getline(fs, sub, ':')) { if (fi == gtIdx) { gtField = sub; break; } fi++; }
            }
            // Take first allele of potentially phased/unphased GT.
            string allele0 = gtField;
            auto sep = gtField.find_first_of("|/");
            if (sep != string::npos) allele0 = gtField.substr(0, sep);

            int g;
            if (allele0.empty() || allele0 == ".") g = -1;  // missing
            else if (allele0 == "*")               g = -2;  // star allele
            else                                   g = stoi(allele0);
            rec.genotypes.push_back(g);
        }

        records.push_back(std::move(rec));
    }
    return records;
}

// Run "vg deconstruct -p focal -a -R <gfaPath>" and return stdout as string.
// -a: all-snarls — reports all snarls at all levels (LV/PS, renamed from -n
//     in vg ≥ v1.73) and adds PA (which allele of parent contains ref path),
//     RS/RD (top-level parent position on focal read), PL/RL/PR tags.
//     PA=0 means the focal read passes through this nested site (it is on the
//     parent's REF branch). PA>0 means the focal read does not enter this site.
// -R: star-alleles for reads that span the parent but don't traverse a nested
//     site — keeps outer and inner events independent for the phasing DP.
string runVgDeconstruct(const string& gfaPath)
{
    const string cmd = "vg deconstruct -p focal -a -R " + gfaPath + " 2>/dev/null";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "";
    string result;
    char buf[8192];
    while (fgets(buf, sizeof(buf), pipe)) result += buf;
    pclose(pipe);
    return result;
}

// Classify variant type from REF/ALT lengths.
Assembler::MSAVariantType classifyVariant(const string& ref, const string& alt)
{
    using VT = Assembler::MSAVariantType;
    if (ref == alt)                             return VT::COMPLEX;
    if (ref.size() == 1 && alt.size() == 1)    return VT::SNP;
    if (ref.empty() || alt.empty())             return VT::DELETION;
    if (ref.size() < alt.size())                return VT::INSERTION;
    if (ref.size() > alt.size())                return VT::DELETION;
    return VT::MNP;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Assembler::computeMSAHetSites
// ---------------------------------------------------------------------------

void Assembler::computeMSAHetSites(ReadId focalReadId, uint32_t strand)
{
    const OrientedReadId focalRead(focalReadId, strand);
    const Reads& reads = getReads();
    const uint64_t readSlot = focalRead.getValue();
    const uint32_t kmerLen = assemblerInfo.isOpen ? uint32_t(assemblerInfo->k) : 0U;

    if (!markers || !shasta2Anchors || readSlot >= shasta2LinearJourneys.size()) {
        std::cout << "Anchor-guided abPOA graph: read " << focalReadId << "-" << strand
                  << " requires markers, anchors, and a linear journey. Skipping.\n";
        return;
    }

    const auto& linearJourney = shasta2LinearJourneys[readSlot];
    if (linearJourney.size() < 2) {
        std::cout << "Anchor-guided abPOA graph: read " << focalReadId << "-" << strand
                  << " linear journey too short. Skipping.\n";
        return;
    }

    std::cout << "Anchor-guided abPOA graph on read "
              << focalReadId << "-" << strand << ":\n"
              << "  Linear journey anchors: " << linearJourney.size() << "\n";

    if (shasta2VariantEvents.size() <= readSlot) {
        shasta2VariantEvents.resize(readSlot + 1);
    }
    shasta2VariantEvents[readSlot].clear();

    const auto extractReadBases = [&](OrientedReadId oid, uint32_t begin, uint32_t end) -> vector<Base> {
        vector<Base> seq;
        if (end <= begin) {
            return seq;
        }
        seq.reserve(size_t(end - begin));
        for (uint32_t pos = begin; pos < end; ++pos) {
            seq.push_back(reads.getOrientedReadBase(oid, pos));
        }
        return seq;
    };

    auto extractBetweenOrdinals = [&](OrientedReadId oid, uint32_t ordA, uint32_t ordB) -> vector<Base> {
        if (ordA == ordB) {
            return {};
        }
        const uint32_t ordMin = std::min(ordA, ordB);
        const uint32_t ordMax = std::max(ordA, ordB);
        const auto orientedMarkers = (*markers)[oid.getValue()];
        if (ordMax >= orientedMarkers.size()) {
            return {};
        }
        const uint32_t posStart = orientedMarkers[ordMin].position;
        const uint32_t posEnd = orientedMarkers[ordMax].position + kmerLen;
        return extractReadBases(oid, posStart, posEnd);
    };

    vector<uint32_t> focalAnchorBasePos(linearJourney.size(), 0);
    for (size_t i = 0; i < linearJourney.size(); ++i) {
        const Shasta2AnchorId anchorId = linearJourney[i];
        const uint32_t ordinal = shasta2Anchors->getOrdinal(anchorId, focalRead);
        focalAnchorBasePos[i] = (*markers)[focalRead.getValue()][ordinal].position;
    }

    struct SharedAnchorHit {
        Shasta2AnchorId anchorId;
        uint32_t focalPos;
        uint32_t readJourneyPos;
        uint32_t readOrdinal;
    };
    struct AnchoredReadCandidate {
        OrientedReadId orientedReadId;
        bool isDirectSeed = false;
        uint32_t chainLength = 0;
        uint32_t focalBegin = 0;
        uint32_t focalEnd = 0;
        uint32_t readBeginOrdinal = 0;
        uint32_t readEndOrdinal = 0;
        uint32_t readJourneyBegin = 0;
        uint32_t readJourneyEnd = 0;
        uint32_t focalSpan = 0;
        vector<SharedAnchorHit> chainHits;
        vector<Base> sequence;
    };

    constexpr uint32_t minSharedAnchors = 2;
    map<Shasta2AnchorId, uint32_t> focalAnchorPosById;
    for (uint32_t focalIndex = 0; focalIndex < linearJourney.size(); ++focalIndex) {
        focalAnchorPosById[linearJourney[focalIndex]] = focalAnchorBasePos[focalIndex];
    }

    set<OrientedReadId> directOverlapReads;
    if (focalRead.getValue() < alignmentTable.size()) {
        for (const uint32_t alignmentId : alignmentTable[focalRead.getValue()]) {
            if (alignmentId >= alignmentData.size()) {
                continue;
            }
            const AlignmentData& ad = alignmentData[alignmentId];
            if (!ad.keptByBothSides()) {
                continue;
            }

            OrientedReadId o0(ad.readIds[0], 0);
            OrientedReadId o1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
            AlignmentInfo orientedInfo = ad.info;
            if (o0.getReadId() != focalRead.getReadId()) {
                std::swap(o0, o1);
                orientedInfo.swap();
            }
            if (o0.getStrand() != focalRead.getStrand()) {
                o0.flipStrand();
                o1.flipStrand();
                orientedInfo.reverseComplement();
            }
            if (o0 != focalRead) {
                continue;
            }
            if (o1.getValue() >= shasta2LinearJourneys.size()) {
                continue;
            }
            directOverlapReads.insert(o1);
        }
    }

    map<Shasta2AnchorId, uint32_t> mappedAnchorPos = focalAnchorPosById;
    set<OrientedReadId> frontierReads = directOverlapReads;
    set<OrientedReadId> placedReads;
    vector<AnchoredReadCandidate> recruitedReads;
    recruitedReads.reserve(directOverlapReads.size());
    uint32_t directSeedCount = 0;
    uint32_t transitiveSeedCount = 0;

    auto buildCandidate = [&](OrientedReadId oid,
                              const map<Shasta2AnchorId, uint32_t>& anchorPosMap,
                              bool isDirectSeed) -> pair<bool, AnchoredReadCandidate> {
        AnchoredReadCandidate candidate;
        if (oid == focalRead || oid.getValue() >= shasta2LinearJourneys.size()) {
            return {false, std::move(candidate)};
        }
        const auto& readJourney = shasta2LinearJourneys[oid.getValue()];
        if (readJourney.size() < minSharedAnchors) {
            return {false, std::move(candidate)};
        }

        vector<SharedAnchorHit> hits;
        hits.reserve(readJourney.size());
        for (uint32_t readJourneyPos = 0; readJourneyPos < readJourney.size(); ++readJourneyPos) {
            const Shasta2AnchorId anchorId = readJourney[readJourneyPos];
            auto it = anchorPosMap.find(anchorId);
            if (it == anchorPosMap.end()) {
                continue;
            }
            const uint32_t readOrdinal = shasta2Anchors->getOrdinal(anchorId, oid);
            hits.push_back(SharedAnchorHit{
                anchorId,
                it->second,
                readJourneyPos,
                readOrdinal
            });
        }
        if (hits.size() < minSharedAnchors) {
            return {false, std::move(candidate)};
        }

        vector<uint32_t> dp(hits.size(), 1);
        vector<int> prev(hits.size(), -1);
        uint32_t bestLen = 0;
        int bestIdx = -1;
        for (size_t i = 0; i < hits.size(); ++i) {
            for (size_t j = 0; j < i; ++j) {
                if (hits[j].focalPos < hits[i].focalPos &&
                    hits[j].readJourneyPos < hits[i].readJourneyPos &&
                    dp[j] + 1 > dp[i]) {
                    dp[i] = dp[j] + 1;
                    prev[i] = int(j);
                }
            }
            if (dp[i] > bestLen) {
                bestLen = dp[i];
                bestIdx = int(i);
            }
        }
        if (bestLen < minSharedAnchors || bestIdx < 0) {
            return {false, std::move(candidate)};
        }

        vector<int> chainIndexes;
        for (int idx = bestIdx; idx >= 0; idx = prev[size_t(idx)]) {
            chainIndexes.push_back(idx);
        }
        std::reverse(chainIndexes.begin(), chainIndexes.end());

        candidate.orientedReadId = oid;
        candidate.isDirectSeed = isDirectSeed;
        candidate.chainLength = bestLen;
        candidate.chainHits.reserve(chainIndexes.size());
        for (const int idx : chainIndexes) {
            candidate.chainHits.push_back(hits[size_t(idx)]);
        }
        const SharedAnchorHit& first = candidate.chainHits.front();
        const SharedAnchorHit& last = candidate.chainHits.back();
        if (last.focalPos <= first.focalPos || last.readOrdinal == first.readOrdinal) {
            return {false, AnchoredReadCandidate{}};
        }

        candidate.focalBegin = first.focalPos;
        candidate.focalEnd = last.focalPos + kmerLen;
        candidate.readBeginOrdinal = first.readOrdinal;
        candidate.readEndOrdinal = last.readOrdinal;
        candidate.readJourneyBegin = first.readJourneyPos;
        candidate.readJourneyEnd = last.readJourneyPos;
        candidate.focalSpan = (candidate.focalEnd > candidate.focalBegin) ?
            (candidate.focalEnd - candidate.focalBegin) : 0;
        if (candidate.focalSpan == 0) {
            return {false, AnchoredReadCandidate{}};
        }

        candidate.sequence = extractBetweenOrdinals(
            oid, candidate.readBeginOrdinal, candidate.readEndOrdinal);
        if (candidate.sequence.empty()) {
            return {false, AnchoredReadCandidate{}};
        }
        return {true, std::move(candidate)};
    };

    auto addFrontierReadsOnPlacedInterval = [&](const AnchoredReadCandidate& candidate) {
        if (candidate.orientedReadId.getValue() >= shasta2LinearJourneys.size()) {
            return;
        }
        const auto& readJourney = shasta2LinearJourneys[candidate.orientedReadId.getValue()];
        const uint32_t start = std::min(candidate.readJourneyBegin, candidate.readJourneyEnd);
        const uint32_t end = std::min<uint32_t>(
            std::max(candidate.readJourneyBegin, candidate.readJourneyEnd),
            uint32_t(readJourney.size() - 1));
        for (uint32_t readJourneyPos = start; readJourneyPos <= end; ++readJourneyPos) {
            const Shasta2AnchorId anchorId = readJourney[readJourneyPos];
            for (const Shasta2AnchorMarkerInfo& info : (*shasta2Anchors)[anchorId]) {
                if (info.orientedReadId == focalRead || placedReads.contains(info.orientedReadId)) {
                    continue;
                }
                if (info.positionInJourney == invalid<uint32_t>) {
                    continue;
                }
                if (info.orientedReadId.getValue() >= shasta2LinearJourneys.size()) {
                    continue;
                }
                frontierReads.insert(info.orientedReadId);
            }
        }
    };

    auto propagateAnchorMappings = [&](const AnchoredReadCandidate& candidate) {
        if (candidate.orientedReadId.getValue() >= shasta2LinearJourneys.size() ||
            candidate.chainHits.empty()) {
            return;
        }
        const auto& readJourney = shasta2LinearJourneys[candidate.orientedReadId.getValue()];
        const auto& orientedMarkers = (*markers)[candidate.orientedReadId.getValue()];

        for (size_t i = 0; i < candidate.chainHits.size(); ++i) {
            const SharedAnchorHit& hit = candidate.chainHits[i];
            mappedAnchorPos.try_emplace(hit.anchorId, hit.focalPos);
            if (i + 1 >= candidate.chainHits.size()) {
                continue;
            }
            const SharedAnchorHit& next = candidate.chainHits[i + 1];
            if (next.readJourneyPos <= hit.readJourneyPos || next.readOrdinal <= hit.readOrdinal) {
                continue;
            }
            const uint32_t readPos0 = orientedMarkers[hit.readOrdinal].position;
            const uint32_t readPos1 = orientedMarkers[next.readOrdinal].position;
            if (readPos1 <= readPos0 || next.focalPos <= hit.focalPos) {
                continue;
            }
            for (uint32_t readJourneyPos = hit.readJourneyPos + 1;
                 readJourneyPos < next.readJourneyPos; ++readJourneyPos) {
                const Shasta2AnchorId anchorId = readJourney[readJourneyPos];
                if (mappedAnchorPos.contains(anchorId)) {
                    continue;
                }
                const uint32_t readOrdinal = shasta2Anchors->getOrdinal(anchorId, candidate.orientedReadId);
                const uint32_t readPos = orientedMarkers[readOrdinal].position;
                const uint64_t numer = uint64_t(readPos - readPos0) * uint64_t(next.focalPos - hit.focalPos);
                const uint64_t denom = uint64_t(readPos1 - readPos0);
                const uint32_t focalPos = hit.focalPos + uint32_t(numer / denom);
                mappedAnchorPos.emplace(anchorId, focalPos);
            }
        }
    };

    auto candidateBetter = [](const AnchoredReadCandidate& a, const AnchoredReadCandidate& b) {
        if (a.focalSpan != b.focalSpan) {
            return a.focalSpan > b.focalSpan;
        }
        if (a.chainLength != b.chainLength) {
            return a.chainLength > b.chainLength;
        }
        if (a.isDirectSeed != b.isDirectSeed) {
            return a.isDirectSeed && !b.isDirectSeed;
        }
        return a.orientedReadId < b.orientedReadId;
    };

    while (true) {
        bool foundBest = false;
        AnchoredReadCandidate bestCandidate;
        for (const OrientedReadId oid : directOverlapReads) {
            if (placedReads.contains(oid)) {
                continue;
            }
            auto [ok, candidate] = buildCandidate(oid, focalAnchorPosById, true);
            if (!ok) {
                continue;
            }
            if (!foundBest || candidateBetter(candidate, bestCandidate)) {
                bestCandidate = std::move(candidate);
                foundBest = true;
            }
        }
        if (!foundBest) {
            break;
        }
        placedReads.insert(bestCandidate.orientedReadId);
        recruitedReads.push_back(bestCandidate);
        ++directSeedCount;
        propagateAnchorMappings(recruitedReads.back());
        addFrontierReadsOnPlacedInterval(recruitedReads.back());
    }

    while (true) {
        bool foundBest = false;
        AnchoredReadCandidate bestCandidate;
        for (const OrientedReadId oid : frontierReads) {
            if (placedReads.contains(oid)) {
                continue;
            }
            const bool isDirectSeed = directOverlapReads.contains(oid);
            auto [ok, candidate] = buildCandidate(oid, mappedAnchorPos, isDirectSeed);
            if (!ok) {
                continue;
            }
            if (!foundBest || candidateBetter(candidate, bestCandidate)) {
                bestCandidate = std::move(candidate);
                foundBest = true;
            }
        }
        if (!foundBest) {
            break;
        }
        placedReads.insert(bestCandidate.orientedReadId);
        recruitedReads.push_back(bestCandidate);
        if (bestCandidate.isDirectSeed) {
            ++directSeedCount;
        } else {
            ++transitiveSeedCount;
        }
        propagateAnchorMappings(recruitedReads.back());
        addFrontierReadsOnPlacedInterval(recruitedReads.back());
    }

    if (recruitedReads.empty()) {
        std::cout << "  No anchor-supported recruited reads.\n";
        return;
    }

    vector<vector<Base>> seqs;
    vector<OrientedReadId> seqIds;
    vector<string> pathNames;
    vector<uint32_t> overlapSpans;
    seqs.reserve(recruitedReads.size() + 1);
    seqIds.reserve(recruitedReads.size() + 1);
    pathNames.reserve(recruitedReads.size() + 1);
    overlapSpans.reserve(recruitedReads.size() + 1);

    const uint32_t focalReadLength = uint32_t(reads.getRead(focalReadId).baseCount);
    seqs.push_back(extractReadBases(focalRead, 0, focalReadLength));
    seqIds.push_back(focalRead);
    pathNames.push_back("focal");
    overlapSpans.push_back(focalReadLength);

    for (const AnchoredReadCandidate& candidate : recruitedReads) {
        seqs.push_back(candidate.sequence);
        seqIds.push_back(candidate.orientedReadId);
        pathNames.push_back("r" + std::to_string(candidate.orientedReadId.getValue()) +
                            (candidate.isDirectSeed ? "_direct" : "_transitive") +
                            "_chain" + std::to_string(candidate.chainLength) +
                            "_span" + std::to_string(candidate.focalSpan));
        overlapSpans.push_back(candidate.focalSpan);
    }

    if (seqs.size() < 2) {
        std::cout << "  Not enough recruited sequences after anchor extraction.\n";
        return;
    }

    std::cout << "  Recruited reads with monotone shared-anchor chains: "
              << recruitedReads.size() << "\n"
              << "    Direct focal-overlap seeds: " << directSeedCount << "\n"
              << "    Transitive recruited reads: " << transitiveSeedCount << "\n"
              << "  Progressive graph input sequences: " << seqs.size()
              << " (focal + " << (seqs.size() - 1) << " recruited)\n"
              << "  Largest recruited focal span: "
              << (overlapSpans.size() > 1 ? overlapSpans[1] : 0) << " bp\n";

    vector<int> sequenceLengths;
    vector<uint8_t*> sequencePointers;
    vector<char*> sequenceNamePointers;
    sequenceLengths.reserve(seqs.size());
    sequencePointers.reserve(seqs.size());
    sequenceNamePointers.reserve(seqs.size());
    for (size_t i = 0; i < seqs.size(); ++i) {
        sequenceLengths.push_back(int(seqs[i].size()));
        sequencePointers.push_back(reinterpret_cast<uint8_t*>(seqs[i].data()));
        sequenceNamePointers.push_back(pathNames[i].data());
    }

    const string stem = "overlap_abpoa_read" + std::to_string(uint64_t(focalReadId));
    const string gfaPath = stem + ".gfa";

    abpoa_t* ab = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    abpt->align_mode = ABPOA_LOCAL_MODE;
    abpt->out_msa = 0;
    abpt->out_cons = 0;
    abpt->out_gfa = 0;
    abpt->sort_input_seq = 0;
    abpt->progressive_poa = 0;
    abpt->cons_algrm = 1;
    abpt->use_qv = 0;
    abpt->max_n_cons = 1;
    abpoa_post_set_para(abpt);

    const int msaStatus = abpoa_msa(
        ab,
        abpt,
        int(seqs.size()),
        sequenceNamePointers.data(),
        sequenceLengths.data(),
        sequencePointers.data(),
        nullptr,
        nullptr);
    if (msaStatus != 0) {
        abpoa_free(ab);
        abpoa_free_para(abpt);
        std::cout << "  abPOA failed with status " << msaStatus << ".\n";
        return;
    }

    FILE* gfaFile = fopen(gfaPath.c_str(), "w");
    if (gfaFile == nullptr) {
        abpoa_free(ab);
        abpoa_free_para(abpt);
        std::cout << "  Failed to open " << gfaPath << " for writing.\n";
        return;
    }
    abpoa_generate_gfa(ab, abpt, gfaFile);
    fclose(gfaFile);
    abpoa_free(ab);
    abpoa_free_para(abpt);

    const string vcfText = runVgDeconstruct(gfaPath);
    if (vcfText.empty()) {
        std::cout << "  Wrote 1 GFA, but vg deconstruct returned no output.\n";
        return;
    }

    vector<string> sampleNames;
    const vector<VcfRecord> vcfRecords = parseVcf(vcfText, sampleNames);
    if (vcfRecords.empty()) {
        std::cout << "  Wrote 1 GFA, but parsed VCF contained no records.\n";
        return;
    }

    map<string, size_t> pathNameToIdx;
    for (size_t i = 0; i < pathNames.size(); ++i) {
        pathNameToIdx[pathNames[i]] = i;
    }

    uint32_t eventsTotal = 0;

    // Diagnostic counters (printed once at end).
    uint32_t dbg_rawEvents     = 0;  // events before any filter
    uint32_t dbg_filtSupport   = 0;  // dropped by min-support (< 3 reads per allele)
    uint32_t dbg_filtStr       = 0;  // dropped by STR filter
    uint32_t dbg_filtStrand    = 0;  // dropped by strand-bias filter

    // Filtered events for verbose output (populated below).
    struct FilteredEvent {
        VariantEvent ev;
        const char*  reason; // "support" | "STR" | "strand"
    };
    std::vector<FilteredEvent> filteredEvents;

    // --- Convert VCF records → VariantEvents ---
    for (const VcfRecord& rec : vcfRecords) {
        for (size_t ai = 0; ai < rec.alts.size(); ++ai) {
            if (rec.alts[ai] == "*") continue;

            VariantEvent ev;
            ev.segmentIndex = 0;
            ev.focalPosStart = 0;
            ev.focalPosEnd = focalReadLength;
            ev.vcfPos = rec.pos;
            ev.vcfId = rec.id;
            ev.refAllele = rec.ref;
            ev.altAllele = rec.alts[ai];
            ev.level = rec.level;
            ev.parentSiteId = rec.parentSite;
            ev.parentAllele = rec.parentAllele;
            ev.topLevelPosStart = rec.rsStart;
            ev.topLevelPosEnd = rec.rsEnd;
            ev.varType = classifyVariant(rec.ref, rec.alts[ai]);

            const int altGT = int(ai + 1);
            ev.refReads.push_back(seqIds[0]);

            for (size_t si = 0; si < rec.genotypes.size(); ++si) {
                if (si >= sampleNames.size()) break;
                const int gt = rec.genotypes[si];
                if (gt == -1 || gt == -2) continue;

                const string& sname = sampleNames[si];
                auto it = pathNameToIdx.find(sname);
                if (it == pathNameToIdx.end()) continue;
                const OrientedReadId oid = seqIds[it->second];
                if (gt == 0) {
                    ev.refReads.push_back(oid);
                } else if (gt == altGT) {
                    ev.altReads.push_back(oid);
                }
            }

            if (ev.refAllele.empty() || ev.altAllele.empty() ||
                ev.refAllele == ev.altAllele) {
                continue;
            }

            ++dbg_rawEvents;

            if (ev.refReads.size() < 3 || ev.altReads.size() < 3) {
                ++dbg_filtSupport;
                filteredEvents.push_back({ev, "support"});
                continue;
            }

            {
                auto isStrSnarl = [](const string& a, const string& b) -> bool {
                    const size_t na = a.size(), nb = b.size();
                    size_t p = 0;
                    while (p < na && p < nb && a[p] == b[p]) ++p;

                    size_t s = 0;
                    const size_t ra = na - p, rb = nb - p;
                    while (s < ra && s < rb && a[na - 1 - s] == b[nb - 1 - s]) ++s;

                    const size_t ta = ra - s, tb = rb - s;
                    if (ta != 0 && tb != 0) return false;
                    if (ta == 0 && tb == 0) return false;

                    const char* tl = (ta == 0) ? b.data() + p : a.data() + p;
                    const size_t n  = (ta == 0) ? tb : ta;

                    for (size_t r = 1; r <= 4; ++r) {
                        if (n % r != 0) continue;
                        const size_t k = n / r;

                        if (k >= 2) {
                            bool ok = true;
                            for (size_t i = r; i < n && ok; ++i) {
                                if (tl[i] != tl[i - r]) ok = false;
                            }
                            if (ok) return true;
                        } else if (r >= 2) {
                            if (p >= r) {
                                bool match = true;
                                for (size_t j = 0; j < r && match; ++j) {
                                    if (a[p - r + j] != tl[j]) match = false;
                                }
                                if (match) return true;
                            }
                            if (s >= r) {
                                bool match = true;
                                for (size_t j = 0; j < r && match; ++j) {
                                    if (a[na - s + j] != tl[j]) match = false;
                                }
                                if (match) return true;
                            }
                        }
                        {
                            const char unit = tl[0];
                            size_t leftRun = 0;
                            for (size_t j = p; j > 0 && a[j - 1] == unit; --j) {
                                ++leftRun;
                            }
                            if (leftRun >= 2) return true;
                            size_t rightRun = 0;
                            for (size_t j = na - s; j < na && a[j] == unit; ++j) {
                                ++rightRun;
                            }
                            if (rightRun >= 2) return true;
                        }
                    }
                    return false;
                };

                if (isStrSnarl(ev.refAllele, ev.altAllele)) {
                    ++dbg_filtStr;
                    filteredEvents.push_back({ev, "STR"});
                    continue;
                }
            }

            {
                constexpr double stRate = 0.05;
                constexpr uint32_t stMax = 2;

                auto isStrandBiased = [&](const std::vector<OrientedReadId>& supportingReads) -> bool {
                    uint32_t fwd = 0;
                    for (const OrientedReadId& oid : supportingReads) {
                        if (oid.getStrand() == 0) ++fwd;
                    }
                    const uint32_t total = uint32_t(supportingReads.size());
                    const uint32_t rev = total - fwd;
                    if (fwd + stMax >= total &&
                        double(total) * stRate + double(fwd) >= double(total)) {
                        return true;
                    }
                    if (rev + stMax >= total &&
                        double(total) * stRate + double(rev) >= double(total)) {
                        return true;
                    }
                    return false;
                };

                if (isStrandBiased(ev.refReads) || isStrandBiased(ev.altReads)) {
                    ++dbg_filtStrand;
                    filteredEvents.push_back({ev, "strand"});
                    continue;
                }
            }

            shasta2VariantEvents[readSlot].push_back(std::move(ev));
            ++eventsTotal;
        }
    }

    // --- Filter 3: drop adjacent-position events (distance == 1 bp) ---
    // Port of hifiasm's Step 0 adjacent-site removal. Events at consecutive
    // positions (vcfPos_i and vcfPos_j differ by exactly 1) are likely MNP
    // decomposition noise or sequencing-error clusters. Drop any event that
    // has a neighbour (same segmentIndex, sorted by vcfPos) at distance 1.
    // Applied globally across all segments for this read.
    {
        auto& evs = shasta2VariantEvents[readSlot];

        // Sort by (segmentIndex, vcfPos) for adjacency detection.
        std::sort(evs.begin(), evs.end(),
            [](const VariantEvent& a, const VariantEvent& b) {
                if (a.segmentIndex != b.segmentIndex)
                    return a.segmentIndex < b.segmentIndex;
                return a.vcfPos < b.vcfPos;
            });

        // Mark events to drop: any event whose vcfPos is within 1 of a neighbour
        // in the same segment.
        std::vector<bool> drop(evs.size(), false);
        for (size_t i = 0; i < evs.size(); i++) {
            // Check previous event in same segment
            if (i > 0 &&
                evs[i].segmentIndex == evs[i-1].segmentIndex &&
                evs[i].vcfPos <= evs[i-1].vcfPos + 1) {
                drop[i]   = true;
                drop[i-1] = true;
            }
            // Check next event in same segment
            if (i + 1 < evs.size() &&
                evs[i].segmentIndex == evs[i+1].segmentIndex &&
                evs[i+1].vcfPos <= evs[i].vcfPos + 1) {
                drop[i]   = true;
                drop[i+1] = true;
            }
        }

        // Compact: remove dropped events.
        // Guard write != i to avoid self-move: GCC's std::string and std::vector
        // move-assign do not check for self-assignment, so evs[n] = std::move(evs[n])
        // silently empties the strings and vectors.
        size_t write = 0;
        for (size_t i = 0; i < evs.size(); i++) {
            if (!drop[i]) {
                if (write != i) evs[write] = std::move(evs[i]);
                ++write;
            }
        }
        evs.resize(write);
        eventsTotal = (uint32_t)evs.size();
    }

    // --- Summary ---
    std::cout << "  Written 1 GFA, "
              << eventsTotal << " variant event(s) (after all filters)\n"
              << "  Diagnostics:\n"
              << "    Raw events (before filters) : " << dbg_rawEvents    << "\n"
              << "    Dropped by min-support (<3) : " << dbg_filtSupport  << "\n"
              << "    Dropped by STR filter       : " << dbg_filtStr      << "\n"
              << "    Dropped by strand-bias      : " << dbg_filtStrand   << "\n";

    auto printEvent = [](const VariantEvent& ev, const char* tag) {
        const char* typeStr = [&]() -> const char* {
            switch (ev.varType) {
                case MSAVariantType::SNP:      return "SNP";
                case MSAVariantType::INSERTION: return "INS";
                case MSAVariantType::DELETION:  return "DEL";
                case MSAVariantType::MNP:       return "MNP";
                default:                        return "COMPLEX";
            }
        }();
        std::cout << "    " << tag
                  << " seg" << ev.segmentIndex
                  << " LV=" << ev.level
                  << " pos=" << ev.vcfPos
                  << " " << typeStr
                  << " ref=" << ev.refAllele
                  << " alt=" << ev.altAllele
                  << " refReads=" << ev.refReads.size()
                  << " altReads=" << ev.altReads.size();
        if (!ev.parentSiteId.empty())
            std::cout << " parent=" << ev.parentSiteId;
        std::cout << "\n";
    };

    std::cout << "  Passing events:\n";
    for (const auto& ev : shasta2VariantEvents[readSlot])
        printEvent(ev, "PASS");

    std::cout << "  Filtered SNPs:\n";
    for (const auto& fev : filteredEvents)
        if (fev.ev.varType == MSAVariantType::SNP)
            printEvent(fev.ev, fev.reason);
}
