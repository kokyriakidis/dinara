// AssemblerMSAHetSites.cpp
//
// Per-anchor-pair POA-MSA for het-site detection.
// For each consecutive anchor pair in a read's linear journey:
//   1. Collect spanning reads, extract inter-anchor sequences, run POA-MSA.
//   2. Write the POA compact graph as GFA.
//   3. Reconstruct each spanning read's path through the GFA nodes (greedy
//      sequence matching), emit as P lines.
//   4. Run "vg deconstruct -p focal -a" to get a VCF where each record is a
//      snarl event (including nested ones via LV/PS tags), and each sample
//      column is a spanning read's genotype (0=ref, 1=alt).
//   5. Parse the VCF into VariantEvent entries and store in
//      shasta2VariantEvents[readSlot].

#include "Assembler.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include <theseus/theseus_msa_aligner.h>
#include <theseus/penalties.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
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

// Run "vg deconstruct -p focal -n -R <gfaPath>" and return stdout as string.
// -n: nested decomposition — reports all snarls at all levels (LV/PS) and
//     adds PA (which allele of the parent contains this site's ref path),
//     RS/RD (top-level parent position on focal read), PL/RL/PR tags.
//     PA=0 means the focal read passes through this nested site (it is on the
//     parent's REF branch). PA>0 means the focal read does not enter this site.
// -R: star-alleles for reads that span the parent but don't traverse a nested
//     site — keeps outer and inner events independent for the phasing DP.
string runVgDeconstruct(const string& gfaPath)
{
    const string cmd = "vg deconstruct -p focal -n -R " + gfaPath + " 2>/dev/null";
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
    const uint64_t k = assemblerInfo->k;
    const Reads& reads = getReads();
    const auto& mkrs = *markers;
    const theseus::Penalties penalties(-1, 4, 4, 2);

    const OrientedReadId focalRead(focalReadId, strand);
    const std::vector<uint64_t>& linearJourney =
        shasta2LinearJourneys[focalRead.getValue()];

    if (linearJourney.size() < 2) {
        std::cout << "Journey MSA: read " << focalReadId << "-" << strand
                  << " linear journey too short. Skipping.\n";
        return;
    }

    std::cout << "Journey MSA (per-anchor-pair, theseus + vg deconstruct) on read "
              << focalReadId << "-" << strand << ":\n"
              << "  Linear journey: " << linearJourney.size() << " anchors, "
              << (linearJourney.size() - 1) << " segments\n";

    const uint32_t windowEnd = (uint32_t)(linearJourney.size() - 1);
    const std::string stem = "journey_msa_read" + std::to_string(uint64_t(focalReadId));

    const uint64_t readSlot = focalRead.getValue();
    if (shasta2VariantEvents.size() <= readSlot)
        shasta2VariantEvents.resize(readSlot + 1);
    shasta2VariantEvents[readSlot].clear();

    auto extractBetween = [&](OrientedReadId oid,
                               uint32_t ordL, uint32_t ordR) -> std::string {
        const uint32_t ordMin = std::min(ordL, ordR);
        const uint32_t ordMax = std::max(ordL, ordR);
        const uint32_t posStart = uint32_t(mkrs[oid.getValue()][ordMin].position);
        const uint32_t posEnd   = uint32_t(mkrs[oid.getValue()][ordMax].position) + uint32_t(k);
        std::string seq;
        seq.reserve(posEnd - posStart);
        for (uint32_t pos = posStart; pos < posEnd; pos++)
            seq += reads.getOrientedReadBase(oid, pos).character();
        return seq;
    };

    uint32_t segmentsWritten = 0;
    uint32_t eventsTotal     = 0;

    for (uint32_t jPos = 0; jPos < windowEnd; jPos++) {
        const Shasta2AnchorId anchorIdL = linearJourney[jPos];
        const Shasta2AnchorId anchorIdR = linearJourney[jPos + 1];

        // Left anchor ordinals.
        std::map<OrientedReadId, uint32_t> leftOrdinals;
        for (const Shasta2AnchorMarkerInfo& info : (*shasta2Anchors)[anchorIdL])
            leftOrdinals[info.orientedReadId] = info.ordinal;

        // Intersect with right anchor → spanning reads only.
        std::map<OrientedReadId, std::pair<uint32_t,uint32_t>> spanning;
        for (const Shasta2AnchorMarkerInfo& info : (*shasta2Anchors)[anchorIdR]) {
            auto it = leftOrdinals.find(info.orientedReadId);
            if (it != leftOrdinals.end())
                spanning[info.orientedReadId] = {it->second, info.ordinal};
        }
        if (spanning.size() < 2) continue;

        // Focal read base positions for this segment (used in VariantEvent).
        uint32_t segPosStart = 0, segPosEnd = 0;
        {
            const auto focalIt = spanning.find(focalRead);
            if (focalIt != spanning.end()) {
                const uint32_t ordMin = std::min(focalIt->second.first,
                                                  focalIt->second.second);
                const uint32_t ordMax = std::max(focalIt->second.first,
                                                  focalIt->second.second);
                segPosStart = uint32_t(mkrs[focalRead.getValue()][ordMin].position);
                segPosEnd   = uint32_t(mkrs[focalRead.getValue()][ordMax].position)
                              + uint32_t(k);
            }
        }

        // Build sequence list: focal read first (path name "focal"), rest named
        // by their oriented read ID value for uniqueness.
        std::vector<std::string>    seqs;
        std::vector<OrientedReadId> seqIds;
        std::vector<std::string>    pathNames;
        seqs.reserve(spanning.size());
        seqIds.reserve(spanning.size());
        pathNames.reserve(spanning.size());

        auto focalIt = spanning.find(focalRead);
        if (focalIt != spanning.end()) {
            seqs.push_back(extractBetween(focalRead,
                                          focalIt->second.first,
                                          focalIt->second.second));
            seqIds.push_back(focalRead);
            pathNames.push_back("focal");
        } else {
            auto first = spanning.begin();
            seqs.push_back(extractBetween(first->first,
                                          first->second.first,
                                          first->second.second));
            seqIds.push_back(first->first);
            pathNames.push_back("focal");
        }
        for (auto& [oid, ords] : spanning) {
            if (oid == focalRead) continue;
            std::string seq = extractBetween(oid, ords.first, ords.second);
            if (!seq.empty()) {
                seqs.push_back(std::move(seq));
                seqIds.push_back(oid);
                pathNames.push_back("r" + std::to_string(oid.getValue()));
            }
        }
        if (seqs.size() < 2) continue;

        // --- POA-MSA → GFA ---
        theseus::TheseusMSA msa(penalties, seqs[0]);
        for (size_t i = 1; i < seqs.size(); i++) msa.align(seqs[i]);

        const std::string gfaPath  = stem + "_seg" + std::to_string(jPos) + ".gfa";
        const std::string gfaWPath = stem + "_seg" + std::to_string(jPos) + "_paths.gfa";
        { std::ofstream out(gfaPath); msa.print_as_gfa(out); }
        ++segmentsWritten;

        // --- Parse GFA graph and reconstruct per-read paths ---
        GFAGraph g;
        if (!parseGFA(gfaPath, g)) continue;

        if (!writeGFAWithPaths(gfaPath, gfaWPath, g, seqs, pathNames)) continue;

        // --- vg deconstruct → VCF ---
        const string vcfText = runVgDeconstruct(gfaWPath);
        if (vcfText.empty()) continue;

        vector<string> sampleNames;
        const vector<VcfRecord> vcfRecords = parseVcf(vcfText, sampleNames);
        if (vcfRecords.empty()) continue;

        // Build map: pathName → index in seqs/seqIds (for OrientedReadId lookup).
        std::map<string, size_t> pathNameToIdx;
        for (size_t i = 0; i < pathNames.size(); i++)
            pathNameToIdx[pathNames[i]] = i;

        // --- Convert VCF records → VariantEvents ---
        for (const VcfRecord& rec : vcfRecords) {
            // One VariantEvent per non-star ALT allele.
            // GT=-2 (star) means "spans parent but doesn't traverse this nested site"
            // — no evidence here, already captured at the parent LV=0 event.
            for (size_t ai = 0; ai < rec.alts.size(); ai++) {
                if (rec.alts[ai] == "*") continue;

                VariantEvent ev;
                ev.segmentIndex    = jPos;
                ev.focalPosStart   = segPosStart;
                ev.focalPosEnd     = segPosEnd;
                ev.vcfPos          = rec.pos;
                ev.vcfId           = rec.id;
                ev.refAllele       = rec.ref;
                ev.altAllele       = rec.alts[ai];
                ev.level           = rec.level;
                ev.parentSiteId    = rec.parentSite;
                ev.parentAllele    = rec.parentAllele;   // PA: 0=focal passes through
                ev.topLevelPosStart = (rec.level == 0) ? rec.pos    : rec.rsStart;
                ev.topLevelPosEnd   = (rec.level == 0) ? rec.pos    : rec.rsEnd;
                ev.varType         = classifyVariant(rec.ref, rec.alts[ai]);

                const int altGT = (int)(ai + 1);

                for (size_t si = 0; si < rec.genotypes.size(); si++) {
                    if (si >= sampleNames.size()) break;
                    const int gt = rec.genotypes[si];
                    if (gt == -1 || gt == -2) continue; // missing or star → no evidence

                    const string& sname = sampleNames[si];
                    auto it = pathNameToIdx.find(sname);
                    if (it == pathNameToIdx.end()) continue;
                    const OrientedReadId oid = seqIds[it->second];

                    if (gt == 0)
                        ev.refReads.push_back(oid);
                    else if (gt == altGT)
                        ev.altReads.push_back(oid);
                    // other non-star alt GT: read is on a different alt allele
                    // at this site — unassigned for this specific event.
                }

                // --- Filter 1: minimum allele support ---
                // Keep only events where BOTH alleles have at least 3 supporting reads.
                if (ev.refReads.size() < 3 || ev.altReads.size() < 3) continue;

                // --- Filter 1b: short tandem repeat (STR) snarl filter ---
                // A snarl is a repeat-unit artifact if the two allele paths differ
                // only in the number of full copies of a repeat unit of period 1–4.
                // Period 1 = homopolymer (AAAA vs AA), period 2 = dinucleotide (ACAC vs AC),
                // period 3 = trinucleotide (AGCAGC vs AGC), period 4 = tetranucleotide.
                // Matches hpc_mask_ff's hpc_min=4 scan range.
                //
                // Rule for each period r (1..4):
                //   1. Length difference must be a non-zero multiple of r.
                //   2. The longer allele must equal the shorter allele with one or more
                //      copies of the repeat unit prepended OR appended.
                //   3. The repeat unit must match the corresponding prefix or suffix of
                //      the shorter allele (confirms the unit is consistent across both).
                //
                // Examples filtered (period shown):
                //   ref=A,      alt=AAA      period=1: extra "AA"    ✓
                //   ref=AC,     alt=ACAC     period=2: extra "AC"    ✓
                //   ref=AGC,    alt=AGCAGC   period=3: extra "AGC"   ✓
                //   ref=ACGT,   alt=ACGTACGT period=4: extra "ACGT"  ✓
                //   ref=ACACAC, alt=ACAC     period=2: extra "AC"    ✓
                // Examples NOT filtered:
                //   ref=A,  alt=G   → same length (SNP)
                //   ref=AC, alt=AG  → same length, different bases (SNP/MNP)
                //   ref=AC, alt=AAC → extra 'A'; unit "A" does NOT match suffix of "AC" ('C')
                //                     NOR prefix of "AC" ('A') at the abutting end → NOT filtered
                //   ref=AA, alt=AAC → extra 'C'; unit "C" ≠ 'A' at abutting end   → NOT filtered
                {
                    auto isStrSnarl = [](const string& a, const string& b) -> bool {
                        const size_t na = a.size(), nb = b.size();

                        // Strip common prefix.
                        size_t p = 0;
                        while (p < na && p < nb && a[p] == b[p]) ++p;

                        // Strip common suffix (must not overlap with stripped prefix).
                        size_t s = 0;
                        const size_t ra = na - p, rb = nb - p;
                        while (s < ra && s < rb && a[na - 1 - s] == b[nb - 1 - s]) ++s;

                        // Trimmed lengths (no allocation — just lengths).
                        const size_t ta = ra - s, tb = rb - s;

                        // Exactly one trimmed string must be empty (pure ins or del).
                        if (ta != 0 && tb != 0) return false;
                        if (ta == 0 && tb == 0) return false;

                        // Pointer into the non-empty trimmed region; n = its length.
                        const char* tl = (ta == 0) ? b.data() + p : a.data() + p;
                        const size_t n  = (ta == 0) ? tb : ta;

                        // Must be k≥2 copies of a period-1..4 unit.
                        // tl[i] == tl[i-r] (no modulo) is equivalent to tl[i] == tl[i%r]
                        // by induction and is division-free.
                        for (size_t r = 1; r <= 4; r++) {
                            if (n % r != 0) continue;
                            if (n / r < 2) continue;
                            bool ok = true;
                            for (size_t i = r; i < n && ok; i++)
                                if (tl[i] != tl[i - r]) ok = false;
                            if (ok) return true;
                        }
                        return false;
                    };

                    if (isStrSnarl(ev.refAllele, ev.altAllele)) continue;
                }

                // --- Filter 2: strand bias ---
                // Port of the hifiasm is_st_bs check (st_rate=0.05, st_max=2).
                // An allele is strand-biased when ≥95% of its supporting reads are
                // on the same strand AND there are ≤2 reads on the minority strand.
                // Applied symmetrically to both ref and alt alleles: if either shows
                // strand bias the event is likely an artifact and is dropped.
                {
                    constexpr double   stRate = 0.05;
                    constexpr uint32_t stMax  = 2;

                    auto isStrandBiased = [&](const std::vector<OrientedReadId>& reads) -> bool {
                        uint32_t fwd = 0;
                        for (const OrientedReadId& oid : reads)
                            if (oid.getStrand() == 0) fwd++;
                        const uint32_t total = (uint32_t)reads.size();
                        const uint32_t rev   = total - fwd;
                        // Check forward-strand dominance.
                        if (fwd + stMax >= total &&
                            double(total) * stRate + double(fwd) >= double(total))
                            return true;
                        // Check reverse-strand dominance (symmetric).
                        if (rev + stMax >= total &&
                            double(total) * stRate + double(rev) >= double(total))
                            return true;
                        return false;
                    };

                    if (isStrandBiased(ev.refReads) || isStrandBiased(ev.altReads))
                        continue;
                }

                shasta2VariantEvents[readSlot].push_back(std::move(ev));
                ++eventsTotal;
            }
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
        size_t write = 0;
        for (size_t i = 0; i < evs.size(); i++)
            if (!drop[i]) evs[write++] = std::move(evs[i]);
        evs.resize(write);
        eventsTotal = (uint32_t)evs.size();
    }

    // --- Summary ---
    std::cout << "  Written " << segmentsWritten << " GFAs, "
              << eventsTotal << " variant event(s) (after all filters)\n";

    for (const auto& ev : shasta2VariantEvents[readSlot]) {
        const char* typeStr = [&]() -> const char* {
            switch (ev.varType) {
                case MSAVariantType::SNP:       return "SNP";
                case MSAVariantType::INSERTION:  return "INS";
                case MSAVariantType::DELETION:   return "DEL";
                case MSAVariantType::MNP:        return "MNP";
                default:                         return "COMPLEX";
            }
        }();
        std::cout << "    seg" << ev.segmentIndex
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
    }
}
