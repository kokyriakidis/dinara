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
#include <theseus/theseus_msa_aligner.h>
#include <theseus/penalties.h>

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



namespace {

string shellQuoteLocal(const string& s)
{
    string quoted = "'";
    for (const char c : s) {
        if (c == '\'') {
            quoted += "'\"'\"'";
        } else {
            quoted += c;
        }
    }
    quoted += "'";
    return quoted;
}

string basesToString(const vector<Base>& seq)
{
    string s;
    s.reserve(seq.size());
    for (const Base& b : seq) {
        s.push_back(b.character());
    }
    return s;
}

bool extractTagValue(const vector<string>& fields, const string& prefix, string& value)
{
    for (const string& field : fields) {
        if (field.rfind(prefix, 0) == 0) {
            value = field.substr(prefix.size());
            return true;
        }
    }
    return false;
}

struct FastGAPafRecord {
    string queryName;
    uint32_t queryStart = 0;
    uint32_t queryEnd = 0;
    char strand = '+';
    string targetName;
    uint32_t targetStart = 0;
    uint32_t targetEnd = 0;
    string cg;
};

bool parseFastGAPafLine(const string& line, FastGAPafRecord& record)
{
    if (line.empty() || line[0] == '#') {
        return false;
    }
    istringstream iss(line);
    vector<string> fields;
    string field;
    while (iss >> field) {
        fields.push_back(field);
    }
    if (fields.size() < 12) {
        return false;
    }

    record.queryName = fields[0];
    record.queryStart = uint32_t(stoul(fields[2]));
    record.queryEnd = uint32_t(stoul(fields[3]));
    record.strand = fields[4].empty() ? '+' : fields[4][0];
    record.targetName = fields[5];
    record.targetStart = uint32_t(stoul(fields[7]));
    record.targetEnd = uint32_t(stoul(fields[8]));
    record.cg.clear();
    extractTagValue(fields, "cg:Z:", record.cg);
    return true;
}

string resolveFastGABinary()
{
    const char* fastgaBinEnv = ::getenv("DINARA_FASTGA_BIN");
    string fastgaBin = (fastgaBinEnv && *fastgaBinEnv) ? string(fastgaBinEnv) : string("FastGA");
    if ((!fastgaBinEnv || !*fastgaBinEnv) && fastgaBin == "FastGA") {
        const char* homeEnv = ::getenv("HOME");
        if (homeEnv && *homeEnv) {
            const string localFastga = string(homeEnv) + "/.local/bin/FastGA";
            if (std::filesystem::exists(localFastga)) {
                fastgaBin = localFastga;
            }
        }
    }
    return fastgaBin;
}

bool envFlagEnabled(const char* name)
{
    const char* value = ::getenv(name);
    if (!value || !*value) {
        return false;
    }
    const string s(value);
    return s == "1" || s == "true" || s == "TRUE" || s == "yes" || s == "YES";
}

} // namespace



void Assembler::benchmarkFastGAOnChainedReadPairs(ReadId focalReadId, uint32_t strand)
{
    using std::chrono::duration_cast;
    using std::chrono::milliseconds;
    using std::chrono::steady_clock;

    const OrientedReadId focalRead(focalReadId, strand);
    const Reads& reads = getReads();
    const uint64_t readSlot = focalRead.getValue();
    const uint32_t kmerLen = assemblerInfo.isOpen ? uint32_t(assemblerInfo->k) : 0U;

    if (!markers || !shasta2Anchors || readSlot >= shasta2LinearJourneys.size()) {
        cout << "FASTGA chained-pair benchmark: read " << focalReadId << "-" << strand
             << " requires markers, anchors, and a linear journey. Skipping." << endl;
        return;
    }

    const auto& focalJourney = shasta2LinearJourneys[readSlot];
    if (focalJourney.size() < 2) {
        cout << "FASTGA chained-pair benchmark: read " << focalReadId << "-" << strand
             << " linear journey too short. Skipping." << endl;
        return;
    }

    string fastgaBin = resolveFastGABinary();
    const char* fastgaThreadsEnv = ::getenv("DINARA_FASTGA_THREADS");
    const uint32_t fastgaThreads = fastgaThreadsEnv ? std::max<uint32_t>(1, uint32_t(stoul(fastgaThreadsEnv))) : 1U;
    const char* fastgaDirEnv = ::getenv("DINARA_FASTGA_DIR");
    const string outputDirectory =
        (fastgaDirEnv && *fastgaDirEnv) ? string(fastgaDirEnv) :
        ("FastGA-read" + to_string(uint64_t(focalReadId)) + "-" + to_string(strand));

    std::filesystem::create_directories(outputDirectory);

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

    auto extractBetweenOrdinals = [&](OrientedReadId oid,
                                      uint32_t ordA,
                                      uint32_t ordB,
                                      uint32_t& posStart,
                                      uint32_t& posEnd) -> vector<Base> {
        posStart = 0;
        posEnd = 0;
        if (ordA == ordB) {
            return {};
        }
        const uint32_t ordMin = min(ordA, ordB);
        const uint32_t ordMax = max(ordA, ordB);
        const auto orientedMarkers = (*markers)[oid.getValue()];
        if (ordMax >= orientedMarkers.size()) {
            return {};
        }
        posStart = orientedMarkers[ordMin].position;
        posEnd = orientedMarkers[ordMax].position + kmerLen;
        return extractReadBases(oid, posStart, posEnd);
    };

    vector<uint32_t> focalAnchorBasePos(focalJourney.size(), 0);
    for (size_t i = 0; i < focalJourney.size(); ++i) {
        const Shasta2AnchorId anchorId = focalJourney[i];
        const uint32_t ordinal = shasta2Anchors->getOrdinal(anchorId, focalRead);
        focalAnchorBasePos[i] = (*markers)[focalRead.getValue()][ordinal].position;
    }

    map<Shasta2AnchorId, uint32_t> focalAnchorPosById;
    for (uint32_t focalIndex = 0; focalIndex < focalJourney.size(); ++focalIndex) {
        focalAnchorPosById[focalJourney[focalIndex]] = focalAnchorBasePos[focalIndex];
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

    struct PairChainHit {
        Shasta2AnchorId anchorId;
        uint32_t focalPos;
        uint32_t sourceJourneyPos;
        uint32_t sourceOrdinal;
        uint32_t targetJourneyPos;
        uint32_t targetOrdinal;
    };
    struct PairTask {
        OrientedReadId sourceReadId;
        OrientedReadId targetReadId;
        bool isDirectSeed = false;
        uint32_t chainLength = 0;
        uint32_t focalBegin = 0;
        uint32_t focalEnd = 0;
        uint32_t focalSpan = 0;
        uint32_t sourceBeginOrdinal = 0;
        uint32_t sourceEndOrdinal = 0;
        uint32_t targetBeginOrdinal = 0;
        uint32_t targetEndOrdinal = 0;
        uint32_t targetJourneyBegin = 0;
        uint32_t targetJourneyEnd = 0;
        uint32_t sourceBeginPos = 0;
        uint32_t sourceEndPos = 0;
        uint32_t targetBeginPos = 0;
        uint32_t targetEndPos = 0;
        vector<PairChainHit> chainHits;
        vector<Base> sourceSequence;
        vector<Base> targetSequence;
    };

    constexpr uint32_t minSharedAnchors = 2;
    map<Shasta2AnchorId, uint32_t> mappedAnchorPos = focalAnchorPosById;
    set<OrientedReadId> placedReads = {focalRead};
    set<OrientedReadId> frontierReads = directOverlapReads;
    vector<PairTask> pairTasks;
    pairTasks.reserve(directOverlapReads.size());

    auto buildPairTask = [&](OrientedReadId sourceReadId,
                             OrientedReadId targetReadId,
                             const map<Shasta2AnchorId, uint32_t>& anchorPosMap,
                             bool isDirectSeed) -> pair<bool, PairTask> {
        PairTask task;
        if (sourceReadId.getValue() >= shasta2LinearJourneys.size() ||
            targetReadId.getValue() >= shasta2LinearJourneys.size() ||
            sourceReadId == targetReadId) {
            return {false, std::move(task)};
        }

        const auto& sourceJourney = shasta2LinearJourneys[sourceReadId.getValue()];
        const auto& targetJourney = shasta2LinearJourneys[targetReadId.getValue()];
        if (sourceJourney.size() < minSharedAnchors || targetJourney.size() < minSharedAnchors) {
            return {false, std::move(task)};
        }

        struct SourceAnchorInfo {
            uint32_t focalPos;
            uint32_t sourceJourneyPos;
            uint32_t sourceOrdinal;
        };
        map<Shasta2AnchorId, SourceAnchorInfo> sourceAnchors;
        for (uint32_t sourceJourneyPos = 0; sourceJourneyPos < sourceJourney.size(); ++sourceJourneyPos) {
            const Shasta2AnchorId anchorId = sourceJourney[sourceJourneyPos];
            auto posIt = anchorPosMap.find(anchorId);
            if (posIt == anchorPosMap.end()) {
                continue;
            }
            sourceAnchors.emplace(
                anchorId,
                SourceAnchorInfo{
                    posIt->second,
                    sourceJourneyPos,
                    shasta2Anchors->getOrdinal(anchorId, sourceReadId)
                });
        }
        if (sourceAnchors.size() < minSharedAnchors) {
            return {false, std::move(task)};
        }

        vector<PairChainHit> hits;
        hits.reserve(targetJourney.size());
        for (uint32_t targetJourneyPos = 0; targetJourneyPos < targetJourney.size(); ++targetJourneyPos) {
            const Shasta2AnchorId anchorId = targetJourney[targetJourneyPos];
            auto it = sourceAnchors.find(anchorId);
            if (it == sourceAnchors.end()) {
                continue;
            }
            hits.push_back(PairChainHit{
                anchorId,
                it->second.focalPos,
                it->second.sourceJourneyPos,
                it->second.sourceOrdinal,
                targetJourneyPos,
                shasta2Anchors->getOrdinal(anchorId, targetReadId)
            });
        }
        if (hits.size() < minSharedAnchors) {
            return {false, std::move(task)};
        }

        vector<uint32_t> dp(hits.size(), 1);
        vector<int> prev(hits.size(), -1);
        uint32_t bestLen = 0;
        int bestIdx = -1;
        for (size_t i = 0; i < hits.size(); ++i) {
            for (size_t j = 0; j < i; ++j) {
                if (hits[j].focalPos < hits[i].focalPos &&
                    hits[j].sourceJourneyPos < hits[i].sourceJourneyPos &&
                    hits[j].targetJourneyPos < hits[i].targetJourneyPos &&
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
            return {false, std::move(task)};
        }

        vector<int> chainIndexes;
        for (int idx = bestIdx; idx >= 0; idx = prev[size_t(idx)]) {
            chainIndexes.push_back(idx);
        }
        reverse(chainIndexes.begin(), chainIndexes.end());

        task.sourceReadId = sourceReadId;
        task.targetReadId = targetReadId;
        task.isDirectSeed = isDirectSeed;
        task.chainLength = bestLen;
        task.chainHits.reserve(chainIndexes.size());
        for (const int idx : chainIndexes) {
            task.chainHits.push_back(hits[size_t(idx)]);
        }

        const PairChainHit& first = task.chainHits.front();
        const PairChainHit& last = task.chainHits.back();
        if (last.focalPos <= first.focalPos ||
            last.sourceOrdinal == first.sourceOrdinal ||
            last.targetOrdinal == first.targetOrdinal) {
            return {false, PairTask{}};
        }

        task.focalBegin = first.focalPos;
        task.focalEnd = last.focalPos + kmerLen;
        task.focalSpan = task.focalEnd - task.focalBegin;
        task.sourceBeginOrdinal = first.sourceOrdinal;
        task.sourceEndOrdinal = last.sourceOrdinal;
        task.targetBeginOrdinal = first.targetOrdinal;
        task.targetEndOrdinal = last.targetOrdinal;
        task.targetJourneyBegin = first.targetJourneyPos;
        task.targetJourneyEnd = last.targetJourneyPos;
        if (task.focalSpan == 0) {
            return {false, PairTask{}};
        }

        task.sourceSequence = extractBetweenOrdinals(
            sourceReadId, task.sourceBeginOrdinal, task.sourceEndOrdinal,
            task.sourceBeginPos, task.sourceEndPos);
        task.targetSequence = extractBetweenOrdinals(
            targetReadId, task.targetBeginOrdinal, task.targetEndOrdinal,
            task.targetBeginPos, task.targetEndPos);
        if (task.sourceSequence.empty() || task.targetSequence.empty()) {
            return {false, PairTask{}};
        }
        return {true, std::move(task)};
    };

    auto propagateAnchorMappings = [&](const PairTask& task) {
        if (task.targetReadId.getValue() >= shasta2LinearJourneys.size()) {
            return;
        }
        const auto& targetJourney = shasta2LinearJourneys[task.targetReadId.getValue()];
        const auto& targetMarkers = (*markers)[task.targetReadId.getValue()];

        for (size_t i = 0; i < task.chainHits.size(); ++i) {
            const PairChainHit& hit = task.chainHits[i];
            mappedAnchorPos.try_emplace(hit.anchorId, hit.focalPos);
            if (i + 1 >= task.chainHits.size()) {
                continue;
            }
            const PairChainHit& next = task.chainHits[i + 1];
            if (next.targetJourneyPos <= hit.targetJourneyPos || next.targetOrdinal <= hit.targetOrdinal) {
                continue;
            }
            const uint32_t targetPos0 = targetMarkers[hit.targetOrdinal].position;
            const uint32_t targetPos1 = targetMarkers[next.targetOrdinal].position;
            if (targetPos1 <= targetPos0 || next.focalPos <= hit.focalPos) {
                continue;
            }
            for (uint32_t targetJourneyPos = hit.targetJourneyPos + 1;
                 targetJourneyPos < next.targetJourneyPos; ++targetJourneyPos) {
                const Shasta2AnchorId anchorId = targetJourney[targetJourneyPos];
                if (mappedAnchorPos.contains(anchorId)) {
                    continue;
                }
                const uint32_t targetOrdinal = shasta2Anchors->getOrdinal(anchorId, task.targetReadId);
                const uint32_t targetPos = targetMarkers[targetOrdinal].position;
                const uint64_t numer =
                    uint64_t(targetPos - targetPos0) * uint64_t(next.focalPos - hit.focalPos);
                const uint64_t denom = uint64_t(targetPos1 - targetPos0);
                mappedAnchorPos.emplace(anchorId, hit.focalPos + uint32_t(numer / denom));
            }
        }
    };

    auto addFrontierReads = [&](const PairTask& task) {
        if (task.targetReadId.getValue() >= shasta2LinearJourneys.size()) {
            return;
        }
        const auto& targetJourney = shasta2LinearJourneys[task.targetReadId.getValue()];
        const uint32_t start = min(task.targetJourneyBegin, task.targetJourneyEnd);
        const uint32_t end = min<uint32_t>(
            max(task.targetJourneyBegin, task.targetJourneyEnd),
            uint32_t(targetJourney.size() - 1));
        for (uint32_t targetJourneyPos = start; targetJourneyPos <= end; ++targetJourneyPos) {
            const Shasta2AnchorId anchorId = targetJourney[targetJourneyPos];
            for (const Shasta2AnchorMarkerInfo& info : (*shasta2Anchors)[anchorId]) {
                if (placedReads.contains(info.orientedReadId) || info.orientedReadId == focalRead) {
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

    auto taskBetter = [](const PairTask& a, const PairTask& b) {
        if (a.focalSpan != b.focalSpan) {
            return a.focalSpan > b.focalSpan;
        }
        if (a.chainLength != b.chainLength) {
            return a.chainLength > b.chainLength;
        }
        if (a.isDirectSeed != b.isDirectSeed) {
            return a.isDirectSeed && !b.isDirectSeed;
        }
        if (a.sourceReadId != b.sourceReadId) {
            return a.sourceReadId < b.sourceReadId;
        }
        return a.targetReadId < b.targetReadId;
    };

    while (true) {
        bool foundBest = false;
        PairTask bestTask;
        for (const OrientedReadId targetReadId : directOverlapReads) {
            if (placedReads.contains(targetReadId)) {
                continue;
            }
            auto [ok, task] = buildPairTask(focalRead, targetReadId, focalAnchorPosById, true);
            if (!ok) {
                continue;
            }
            if (!foundBest || taskBetter(task, bestTask)) {
                bestTask = std::move(task);
                foundBest = true;
            }
        }
        if (!foundBest) {
            break;
        }
        placedReads.insert(bestTask.targetReadId);
        pairTasks.push_back(bestTask);
        propagateAnchorMappings(pairTasks.back());
        addFrontierReads(pairTasks.back());
    }

    while (true) {
        bool foundBest = false;
        PairTask bestTask;
        for (const OrientedReadId targetReadId : frontierReads) {
            if (placedReads.contains(targetReadId)) {
                continue;
            }
            for (const OrientedReadId sourceReadId : placedReads) {
                if (sourceReadId == targetReadId) {
                    continue;
                }
                const bool isDirectSeed = (sourceReadId == focalRead) && directOverlapReads.contains(targetReadId);
                auto [ok, task] = buildPairTask(sourceReadId, targetReadId, mappedAnchorPos, isDirectSeed);
                if (!ok) {
                    continue;
                }
                if (!foundBest || taskBetter(task, bestTask)) {
                    bestTask = std::move(task);
                    foundBest = true;
                }
            }
        }
        if (!foundBest) {
            break;
        }
        placedReads.insert(bestTask.targetReadId);
        pairTasks.push_back(bestTask);
        propagateAnchorMappings(pairTasks.back());
        addFrontierReads(pairTasks.back());
    }

    cout << "FASTGA chained-pair benchmark on read " << focalReadId << "-" << strand << ":" << endl
         << "  Linear journey anchors: " << focalJourney.size() << endl
         << "  Direct focal overlaps: " << directOverlapReads.size() << endl
         << "  Chained read pairs to test: " << pairTasks.size() << endl
         << "  Output directory: " << std::filesystem::absolute(outputDirectory) << endl;

    if (pairTasks.empty()) {
        cout << "  No chained read pairs were identified." << endl;
        return;
    }

    const bool fastgaAvailable =
        std::filesystem::exists(fastgaBin) ||
        (::system(("command -v " + shellQuoteLocal(fastgaBin) + " >/dev/null 2>&1").c_str()) == 0);
    if (!fastgaAvailable) {
        cout << "  FASTGA binary not found: " << fastgaBin
             << " (set DINARA_FASTGA_BIN to override)." << endl;
        return;
    }

    ofstream summary(outputDirectory + "/summary.tsv");
    ofstream events(outputDirectory + "/events.tsv");
    summary << "pairIndex\tsource\ttarget\tdirect\tchainLength\tfocalSpan\tsourceBegin\tsourceEnd\ttargetBegin\ttargetEnd\telapsedMs\texitCode\tpafRecords\tsnpCount\tinsCount\tdelCount\tsvCount\n";
    events << "pairIndex\tsource\ttarget\tqueryName\ttargetName\tqueryPos\ttargetPos\ttype\tsize\n";

    uint64_t totalElapsedMs = 0;
    uint64_t totalPafRecords = 0;
    uint64_t totalSnps = 0;
    uint64_t totalInsertions = 0;
    uint64_t totalDeletions = 0;
    uint64_t totalSvs = 0;
    uint64_t succeededPairs = 0;

    for (size_t pairIndex = 0; pairIndex < pairTasks.size(); ++pairIndex) {
        const PairTask& task = pairTasks[pairIndex];
        const string pairStem = outputDirectory + "/pair" + to_string(pairIndex);
        const string sourceFasta = pairStem + ".source.fa";
        const string targetFasta = pairStem + ".target.fa";
        const string pafPath = pairStem + ".paf";
        const string stderrPath = pairStem + ".stderr";

        {
            ofstream out(sourceFasta);
            out << ">r" << task.sourceReadId.getValue() << "\n"
                << basesToString(task.sourceSequence) << "\n";
        }
        {
            ofstream out(targetFasta);
            out << ">r" << task.targetReadId.getValue() << "\n"
                << basesToString(task.targetSequence) << "\n";
        }

        const string cmd =
            shellQuoteLocal(fastgaBin) +
            " -T" + to_string(fastgaThreads) +
            " -pafxS " +
            shellQuoteLocal(sourceFasta) + " " +
            shellQuoteLocal(targetFasta) +
            " > " + shellQuoteLocal(pafPath) +
            " 2> " + shellQuoteLocal(stderrPath);

        const auto t0 = steady_clock::now();
        const int exitCode = ::system(cmd.c_str());
        const auto elapsedMs = uint64_t(duration_cast<milliseconds>(steady_clock::now() - t0).count());
        totalElapsedMs += elapsedMs;

        uint64_t pafRecords = 0;
        uint64_t snpCount = 0;
        uint64_t insCount = 0;
        uint64_t delCount = 0;
        uint64_t svCount = 0;

        if (exitCode == 0) {
            ++succeededPairs;
            ifstream pafIn(pafPath);
            string line;
            while (getline(pafIn, line)) {
                FastGAPafRecord rec;
                if (!parseFastGAPafLine(line, rec)) {
                    continue;
                }
                ++pafRecords;
                if (rec.cg.empty()) {
                    continue;
                }

                uint32_t queryPos = task.sourceBeginPos + rec.queryStart;
                uint32_t targetPos = task.targetBeginPos + rec.targetStart;
                size_t i = 0;
                while (i < rec.cg.size()) {
                    size_t j = i;
                    while (j < rec.cg.size() && isdigit(rec.cg[j])) {
                        ++j;
                    }
                    if (j == i || j >= rec.cg.size()) {
                        break;
                    }
                    const uint32_t len = uint32_t(stoul(rec.cg.substr(i, j - i)));
                    const char op = rec.cg[j];
                    if (op == '=' || op == 'M') {
                        queryPos += len;
                        targetPos += len;
                    } else if (op == 'X') {
                        for (uint32_t k = 0; k < len; ++k) {
                            events << pairIndex << "\t"
                                   << task.sourceReadId.getValue() << "\t"
                                   << task.targetReadId.getValue() << "\t"
                                   << rec.queryName << "\t"
                                   << rec.targetName << "\t"
                                   << (queryPos + k) << "\t"
                                   << (targetPos + k) << "\tSNP\t1\n";
                        }
                        snpCount += len;
                        queryPos += len;
                        targetPos += len;
                    } else if (op == 'I') {
                        events << pairIndex << "\t"
                               << task.sourceReadId.getValue() << "\t"
                               << task.targetReadId.getValue() << "\t"
                               << rec.queryName << "\t"
                               << rec.targetName << "\t"
                               << queryPos << "\t"
                               << targetPos << "\t"
                               << (len >= 50 ? "SV_INS" : "INS") << "\t"
                               << len << "\n";
                        ++insCount;
                        if (len >= 50) {
                            ++svCount;
                        }
                        queryPos += len;
                    } else if (op == 'D') {
                        events << pairIndex << "\t"
                               << task.sourceReadId.getValue() << "\t"
                               << task.targetReadId.getValue() << "\t"
                               << rec.queryName << "\t"
                               << rec.targetName << "\t"
                               << queryPos << "\t"
                               << targetPos << "\t"
                               << (len >= 50 ? "SV_DEL" : "DEL") << "\t"
                               << len << "\n";
                        ++delCount;
                        if (len >= 50) {
                            ++svCount;
                        }
                        targetPos += len;
                    } else {
                        if (op == 'S' || op == 'H') {
                            queryPos += len;
                        } else if (op == 'N') {
                            targetPos += len;
                        }
                    }
                    i = j + 1;
                }
            }
        }

        totalPafRecords += pafRecords;
        totalSnps += snpCount;
        totalInsertions += insCount;
        totalDeletions += delCount;
        totalSvs += svCount;

        summary << pairIndex << "\t"
                << task.sourceReadId.getValue() << "\t"
                << task.targetReadId.getValue() << "\t"
                << (task.isDirectSeed ? 1 : 0) << "\t"
                << task.chainLength << "\t"
                << task.focalSpan << "\t"
                << task.sourceBeginPos << "\t"
                << task.sourceEndPos << "\t"
                << task.targetBeginPos << "\t"
                << task.targetEndPos << "\t"
                << elapsedMs << "\t"
                << exitCode << "\t"
                << pafRecords << "\t"
                << snpCount << "\t"
                << insCount << "\t"
                << delCount << "\t"
                << svCount << "\n";
    }

    cout << "  FASTGA binary: " << fastgaBin << endl
         << "  FASTGA threads per pair: " << fastgaThreads << endl
         << "  Successful pairs: " << succeededPairs << " / " << pairTasks.size() << endl
         << "  Total wall time across pairs: " << totalElapsedMs << " ms" << endl
         << "  Total PAF records: " << totalPafRecords << endl
         << "  Derived event counts from extended CIGAR:" << endl
         << "    SNPs: " << totalSnps << endl
         << "    Insertions: " << totalInsertions << endl
         << "    Deletions: " << totalDeletions << endl
         << "    SV-sized indels (>=50bp): " << totalSvs << endl
         << "  Wrote: " << std::filesystem::absolute(outputDirectory + "/summary.tsv") << endl
         << "  Wrote: " << std::filesystem::absolute(outputDirectory + "/events.tsv") << endl;
}



void Assembler::alignChainedCandidatesWithFastGA(uint64_t threadCount)
{
    using std::chrono::duration_cast;
    using std::chrono::milliseconds;
    using std::chrono::steady_clock;

    const auto tBegin = steady_clock::now();
    const size_t candidateCount = alignmentCandidates.candidates.size();

    cout << timestamp << "Begin FASTGA alignment of chained candidates for "
         << candidateCount << " candidate pairs." << endl;

    reads->checkReadsAreOpen();
    checkAlignmentCandidatesAreOpen();

    if (alignmentCandidatesAlignmentsData.alignments.size() < candidateCount) {
        throw runtime_error(
            "alignChainedCandidatesWithFastGA requires precomputed chained intervals "
            "for every candidate.");
    }

    const string fastgaBin = resolveFastGABinary();
    const bool fastgaAvailable =
        std::filesystem::exists(fastgaBin) ||
        (::system(("command -v " + shellQuoteLocal(fastgaBin) + " >/dev/null 2>&1").c_str()) == 0);
    if (!fastgaAvailable) {
        throw runtime_error(
            "alignChainedCandidatesWithFastGA could not find the FastGA binary. "
            "Set DINARA_FASTGA_BIN to override.");
    }

    const char* fastgaThreadsEnv = ::getenv("DINARA_FASTGA_THREADS");
    const uint32_t fastgaThreads =
        fastgaThreadsEnv ?
        std::max<uint32_t>(1, uint32_t(stoul(fastgaThreadsEnv))) :
        std::max<uint32_t>(1, uint32_t(threadCount == 0 ? 1 : threadCount));
    const char* outputDirEnv = ::getenv("DINARA_FASTGA_CHAINED_DIR");
    const string outputDirectory =
        (outputDirEnv && *outputDirEnv) ? string(outputDirEnv) : string("FastGA-chained-candidates");
    const bool keepTempFiles = envFlagEnabled("DINARA_FASTGA_KEEP_TEMP");
    const bool writeEvents = envFlagEnabled("DINARA_FASTGA_WRITE_EVENTS");
    const char* maxCandidatesEnv = ::getenv("DINARA_FASTGA_MAX_CANDIDATES");
    const size_t maxCandidates =
        (maxCandidatesEnv && *maxCandidatesEnv) ?
        size_t(stoull(maxCandidatesEnv)) :
        candidateCount;

    std::filesystem::create_directories(outputDirectory);
    const string stderrPath = outputDirectory + "/fastga.stderr";

    ofstream summary(outputDirectory + "/summary.tsv");
    if (!summary) {
        throw runtime_error("Could not open FASTGA summary output in " + outputDirectory);
    }
    summary <<
        "candidateIndex\tread0\tread1\tsameStrand\tqs\tqe\tts\tte\tquerySpan\ttargetSpan\t"
        "elapsedMs\texitCode\tpafRecords\tsnpCount\tinsCount\tdelCount\tsvCount\n";

    unique_ptr<ofstream> eventsFile;
    if (writeEvents) {
        eventsFile = make_unique<ofstream>(outputDirectory + "/events.tsv");
        if (!*eventsFile) {
            throw runtime_error("Could not open FASTGA events output in " + outputDirectory);
        }
        *eventsFile <<
            "candidateIndex\tread0\tread1\tqueryPos\ttargetPos\ttype\tsize\t"
            "queryBase\ttargetBase\n";
    }

    auto extractReadBases = [&](OrientedReadId oid, uint32_t begin, uint32_t end) -> vector<Base> {
        vector<Base> sequence;
        if (end <= begin) {
            return sequence;
        }
        sequence.reserve(size_t(end - begin));
        for (uint32_t pos = begin; pos < end; ++pos) {
            sequence.push_back(reads->getOrientedReadBase(oid, pos));
        }
        return sequence;
    };

    struct PendingSnp {
        uint32_t pos;
        Base base;
    };
    struct PendingIndel {
        uint32_t pos;
        uint32_t len;
        uint8_t type;
    };
    auto sortSnps = [](const PendingSnp& a, const PendingSnp& b) {
        if (a.pos != b.pos) {
            return a.pos < b.pos;
        }
        return a.base < b.base;
    };
    auto sortIndels = [](const PendingIndel& a, const PendingIndel& b) {
        if (a.pos != b.pos) {
            return a.pos < b.pos;
        }
        if (a.type != b.type) {
            return a.type < b.type;
        }
        return a.len < b.len;
    };

    chainedFastgaEvidenceStore.clear();
    chainedFastgaEvidenceStore.index.reserve(std::min(candidateCount, maxCandidates));

    uint64_t totalElapsedMs = 0;
    uint64_t totalPafRecords = 0;
    uint64_t totalSnps = 0;
    uint64_t totalInsertions = 0;
    uint64_t totalDeletions = 0;
    uint64_t totalSvs = 0;
    uint64_t successfulCandidates = 0;

    for (size_t i = 0; i < candidateCount && i < maxCandidates; ++i) {
        const OrientedReadPair& candidate = alignmentCandidates.candidates[i];
        const Alignment& alignment = alignmentCandidatesAlignmentsData.alignments[i];
        const OrientedReadId orientedReadId0(candidate.readIds[0], 0);
        const OrientedReadId orientedReadId1(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);

        const uint32_t readLen0 = uint32_t(reads->getRead(candidate.readIds[0]).baseCount);
        const uint32_t readLen1 = uint32_t(reads->getRead(candidate.readIds[1]).baseCount);
        const uint32_t qs = std::min(alignment.qs, readLen0);
        const uint32_t qe = std::min(alignment.qe, readLen0);
        const uint32_t tsForward = std::min(alignment.ts, readLen1);
        const uint32_t teForward = std::min(alignment.te, readLen1);

        chainedFastgaEvidenceStore.beginAlignment();

        if (qe <= qs || teForward <= tsForward) {
            summary << i << "\t" << candidate.readIds[0] << "\t" << candidate.readIds[1] << "\t"
                    << (candidate.isSameStrand ? 1 : 0) << "\t"
                    << qs << "\t" << qe << "\t" << tsForward << "\t" << teForward << "\t"
                    << (qe > qs ? (qe - qs) : 0) << "\t"
                    << (teForward > tsForward ? (teForward - tsForward) : 0) << "\t"
                    << 0 << "\t-1\t0\t0\t0\t0\t0\n";
            continue;
        }

        const uint32_t targetBegin = candidate.isSameStrand ? tsForward : (readLen1 - teForward);
        const uint32_t targetEnd = candidate.isSameStrand ? teForward : (readLen1 - tsForward);
        vector<Base> querySequence = extractReadBases(orientedReadId0, qs, qe);
        vector<Base> targetSequence = extractReadBases(orientedReadId1, targetBegin, targetEnd);
        if (querySequence.empty() || targetSequence.empty()) {
            summary << i << "\t" << candidate.readIds[0] << "\t" << candidate.readIds[1] << "\t"
                    << (candidate.isSameStrand ? 1 : 0) << "\t"
                    << qs << "\t" << qe << "\t" << tsForward << "\t" << teForward << "\t"
                    << (qe - qs) << "\t" << (teForward - tsForward) << "\t"
                    << 0 << "\t-1\t0\t0\t0\t0\t0\n";
            continue;
        }

        const string pairStem = outputDirectory + "/candidate" + to_string(i);
        const string queryFasta = pairStem + ".query.fa";
        const string targetFasta = pairStem + ".target.fa";
        const string pafPath = pairStem + ".paf";

        {
            ofstream out(queryFasta);
            out << ">read" << candidate.readIds[0] << "\n"
                << basesToString(querySequence) << "\n";
        }
        {
            ofstream out(targetFasta);
            out << ">read" << candidate.readIds[1]
                << (candidate.isSameStrand ? "+\n" : "-\n")
                << basesToString(targetSequence) << "\n";
        }

        const string cmd =
            shellQuoteLocal(fastgaBin) +
            " -T" + to_string(fastgaThreads) +
            " -pafx " +
            shellQuoteLocal(queryFasta) + " " +
            shellQuoteLocal(targetFasta) +
            " > " + shellQuoteLocal(pafPath) +
            " 2>> " + shellQuoteLocal(stderrPath);

        const auto t0 = steady_clock::now();
        const int exitCode = ::system(cmd.c_str());
        const uint64_t elapsedMs =
            uint64_t(duration_cast<milliseconds>(steady_clock::now() - t0).count());
        totalElapsedMs += elapsedMs;

        uint64_t pafRecords = 0;
        uint64_t snpCount = 0;
        uint64_t insCount = 0;
        uint64_t delCount = 0;
        uint64_t svCount = 0;

        if (exitCode == 0) {
            vector<PendingSnp> pendingSnps0;
            vector<PendingSnp> pendingSnps1;
            vector<PendingIndel> pendingIndels0;
            vector<PendingIndel> pendingIndels1;

            ifstream pafIn(pafPath);
            string line;
            while (getline(pafIn, line)) {
                FastGAPafRecord rec;
                if (!parseFastGAPafLine(line, rec) || rec.cg.empty()) {
                    continue;
                }
                ++pafRecords;

                uint32_t queryPos = rec.queryStart;
                uint32_t targetPos = rec.targetStart;
                size_t cgPos = 0;
                while (cgPos < rec.cg.size()) {
                    size_t digitsEnd = cgPos;
                    while (digitsEnd < rec.cg.size() && isdigit(rec.cg[digitsEnd])) {
                        ++digitsEnd;
                    }
                    if (digitsEnd == cgPos || digitsEnd >= rec.cg.size()) {
                        break;
                    }
                    const uint32_t len = uint32_t(stoul(rec.cg.substr(cgPos, digitsEnd - cgPos)));
                    const char op = rec.cg[digitsEnd];

                    if (op == '=' || op == 'M') {
                        queryPos += len;
                        targetPos += len;
                    } else if (op == 'X') {
                        for (uint32_t k = 0; k < len; ++k) {
                            if (queryPos + k >= querySequence.size() || targetPos + k >= targetSequence.size()) {
                                break;
                            }
                            const Base queryBase = querySequence[queryPos + k];
                            const Base targetBaseOriented = targetSequence[targetPos + k];
                            const uint32_t queryAbs = qs + queryPos + k;
                            pendingSnps1.push_back({queryAbs, targetBaseOriented});

                            const uint32_t targetForwardPos =
                                candidate.isSameStrand ?
                                (tsForward + targetPos + k) :
                                (teForward - 1U - (targetPos + k));
                            pendingSnps0.push_back({
                                targetForwardPos,
                                candidate.isSameStrand ? queryBase : queryBase.complement()
                            });

                            if (eventsFile) {
                                *eventsFile << i << "\t" << candidate.readIds[0] << "\t"
                                            << candidate.readIds[1] << "\t"
                                            << queryAbs << "\t" << targetForwardPos
                                            << "\tSNP\t1\t" << queryBase.character() << "\t"
                                            << (candidate.isSameStrand ?
                                                targetBaseOriented.character() :
                                                targetBaseOriented.complement().character())
                                            << "\n";
                            }
                        }
                        snpCount += len;
                        queryPos += len;
                        targetPos += len;
                    } else if (op == 'I') {
                        const uint32_t queryAbs = qs + queryPos;
                        pendingIndels1.push_back({queryAbs, len, 1});

                        const uint32_t targetBoundary =
                            candidate.isSameStrand ?
                            (tsForward + targetPos) :
                            (teForward - targetPos);
                        pendingIndels0.push_back({targetBoundary, len, 0});

                        if (eventsFile) {
                            *eventsFile << i << "\t" << candidate.readIds[0] << "\t"
                                        << candidate.readIds[1] << "\t"
                                        << queryAbs << "\t" << targetBoundary
                                        << "\tDEL\t" << len << "\t.\t.\n";
                        }

                        ++insCount;
                        if (len >= 50) {
                            ++svCount;
                        }
                        queryPos += len;
                    } else if (op == 'D') {
                        const uint32_t queryAbs = qs + queryPos;
                        pendingIndels1.push_back({queryAbs, len, 0});

                        const uint32_t targetForwardPos =
                            candidate.isSameStrand ?
                            (tsForward + targetPos) :
                            (teForward - (targetPos + len));
                        pendingIndels0.push_back({targetForwardPos, len, 1});

                        if (eventsFile) {
                            *eventsFile << i << "\t" << candidate.readIds[0] << "\t"
                                        << candidate.readIds[1] << "\t"
                                        << queryAbs << "\t" << targetForwardPos
                                        << "\tINS\t" << len << "\t.\t.\n";
                        }

                        ++delCount;
                        if (len >= 50) {
                            ++svCount;
                        }
                        targetPos += len;
                    } else if (op == 'S' || op == 'H') {
                        queryPos += len;
                    } else if (op == 'N') {
                        targetPos += len;
                    }

                    cgPos = digitsEnd + 1;
                }
            }

            std::sort(pendingSnps0.begin(), pendingSnps0.end(), sortSnps);
            std::sort(pendingSnps1.begin(), pendingSnps1.end(), sortSnps);
            std::sort(pendingIndels0.begin(), pendingIndels0.end(), sortIndels);
            std::sort(pendingIndels1.begin(), pendingIndels1.end(), sortIndels);

            for (const PendingSnp& snp : pendingSnps0) {
                chainedFastgaEvidenceStore.addSnp0(snp.pos, snp.base.value);
            }
            for (const PendingIndel& indel : pendingIndels0) {
                chainedFastgaEvidenceStore.addIndel0(indel.pos, indel.len, indel.type);
            }
            for (const PendingSnp& snp : pendingSnps1) {
                chainedFastgaEvidenceStore.addSnp1(snp.pos, snp.base.value);
            }
            for (const PendingIndel& indel : pendingIndels1) {
                chainedFastgaEvidenceStore.addIndel1(indel.pos, indel.len, indel.type);
            }

            ++successfulCandidates;
        }

        totalPafRecords += pafRecords;
        totalSnps += snpCount;
        totalInsertions += insCount;
        totalDeletions += delCount;
        totalSvs += svCount;

        summary << i << "\t" << candidate.readIds[0] << "\t" << candidate.readIds[1] << "\t"
                << (candidate.isSameStrand ? 1 : 0) << "\t"
                << qs << "\t" << qe << "\t" << tsForward << "\t" << teForward << "\t"
                << (qe - qs) << "\t" << (teForward - tsForward) << "\t"
                << elapsedMs << "\t" << exitCode << "\t"
                << pafRecords << "\t" << snpCount << "\t"
                << insCount << "\t" << delCount << "\t" << svCount << "\n";

        if (!keepTempFiles) {
            std::error_code ec;
            std::filesystem::remove(queryFasta, ec);
            std::filesystem::remove(targetFasta, ec);
            std::filesystem::remove(pafPath, ec);
        }
    }

    const auto tEnd = steady_clock::now();
    cout << timestamp << "FASTGA chained-candidate benchmark completed." << endl
         << timestamp << "  FASTGA binary: " << fastgaBin << endl
         << timestamp << "  FASTGA threads per pair: " << fastgaThreads << endl
         << timestamp << "  Processed candidates: " << std::min(candidateCount, maxCandidates)
         << " / " << candidateCount << endl
         << timestamp << "  Successful candidates: " << successfulCandidates << endl
         << timestamp << "  Total pair wall time: " << totalElapsedMs << " ms" << endl
         << timestamp << "  Total elapsed wall time: "
         << std::chrono::duration<double>(tEnd - tBegin).count() << " s" << endl
         << timestamp << "  Total PAF records: " << totalPafRecords << endl
         << timestamp << "  Derived event counts:" << endl
         << timestamp << "    SNPs: " << totalSnps << endl
         << timestamp << "    Insertions: " << totalInsertions << endl
         << timestamp << "    Deletions: " << totalDeletions << endl
         << timestamp << "    SV-sized indels (>=50bp): " << totalSvs << endl
         << timestamp << "  Evidence entries stored: " << chainedFastgaEvidenceStore.index.size() << endl
         << timestamp << "  Summary: " << std::filesystem::absolute(outputDirectory + "/summary.tsv") << endl;
    if (eventsFile) {
        cout << timestamp << "  Events: "
             << std::filesystem::absolute(outputDirectory + "/events.tsv") << endl;
    }
}
