// AssemblerAbpoaMultiSegmentMSA.cpp
//
// Per-window multi-segment MSA, built on abPOA. The backbone read is laid down
// as a linear partial-order graph (one node per base). Every overlapping read
// is then folded in piece-by-piece between the backbone nodes that correspond
// to the shared anchors it carries, using abPOA's subgraph alignment API
// (abpoa_align_sequence_to_subgraph + abpoa_add_subgraph_alignment).
//
// This replaces the earlier theseus-fork implementation (AssemblerMultiSegment-
// MSA.cpp, excluded from the build). abPOA's banded SIMD DP makes each piece
// alignment O(qlen * band) rather than the fork's scalar WFA O(score^2), which
// removes the quadratic blow-up on reads whose shared anchors are far apart
// (~4.6x faster overall) and avoids the fork's terminal-condition issues.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "PhasingKmeansTypes.hpp"  // kmIsHomopolymer (homopolymer/STR context)
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"

#include <abpoa.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// Extract the base sequence of an oriented read between two marker ordinals as
// abPOA 0-3 base codes (A=0,C=1,G=2,T=3). Mirrors extractSegmentSequence in
// AssemblerMultiSegmentMSA.cpp but emits codes directly (Base::value already is
// the 0-3 encoding abPOA expects).
vector<uint8_t> extractSegmentCodes(
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t k,
    OrientedReadId oid,
    uint32_t ordinalA,
    uint32_t ordinalB)
{
    if(ordinalA >= ordinalB) {
        return {};
    }
    const auto readMarkers = markers[oid.getValue()];
    if(ordinalB >= readMarkers.size()) {
        return {};
    }
    const uint32_t kHalf = uint32_t(k / 2);
    const uint32_t beginPos = readMarkers[ordinalA].position + kHalf;
    const uint32_t endPos   = readMarkers[ordinalB].position + kHalf;
    if(endPos <= beginPos) {
        return {};
    }
    vector<uint8_t> codes;
    codes.reserve(endPos - beginPos);
    for(uint32_t pos = beginPos; pos < endPos; pos++) {
        codes.push_back(reads.getOrientedReadBase(oid, pos).value);
    }
    return codes;
}

// One member of a het-anchor allele: an oriented read plus the base position
// (in that read's own strand-0/strand coordinates) at the bubble's common
// predecessor node. This is the rawPosition of the k=2 het anchor marker
// [predBase, alleleBase]: shasta2 re-derives the 2-base k-mer from this
// position, so all members of one allele share the same 2 bases and members of
// different alleles differ in the second base (the bubble).
struct HetAlleleMember {
    OrientedReadId orientedReadId;
    uint32_t rawPosition;             // read base position at commonPred (predBase)
};

// One allele of a detected het site. The reference (backbone) allele and every
// alternate allele that independently clears minSupport/minVAF is represented,
// so a triallelic column becomes a single 3-armed bubble instead of two
// biallelic bubbles that duplicate the reference.
struct WindowSnpAllele {
    char base = 'N';                  // allele base (A/C/G/T)
    uint8_t code = 0;                 // allele base code (0-3)
    int support = 0;                  // reads carrying this allele
    bool isRef = false;               // true for the backbone/reference allele
    vector<OrientedReadId> reads;     // reads on this allele (phasing)
    // Het-anchor members (read + rawPosition at commonPred). Subset of `reads`
    // whose predecessor position was recoverable from the POA.
    vector<HetAlleleMember> members;
};

// A detected clean SNP site (biallelic or multiallelic) in a window POA graph.
struct WindowSnp {
    int msaRank = -1;                 // MSA column (topological rank) of the bubble
    int backboneOffset = -1;          // backbone base offset (into backboneCodes), -1 if backbone gaps here
    int delSupport = 0;               // reads carrying a coexisting deletion (skip edge)
    int spanning = 0;                 // reads spanning the column (VAF denominator)
    bool inHomopolymerOrRepeat = false; // backbone context is a homopolymer/STR run
    // All alleles clearing the gates: [0] = reference, [1..] = alternates. Every
    // entry (including ref) has support >= minSupport; each alt also has
    // vaf >= minVaf. Size >= 2 (ref + >=1 alt) whenever the site is emitted.
    vector<WindowSnpAllele> alleles;

    // Leading hom anchor members: reads at predPrev (the base BEFORE the
    // bubble's common predecessor), rawPosition = predPrev base position,
    // forming the k=2 hom anchor [predPrevBase, predBase]. The flank-linearity
    // test guarantees predPrev->commonPred is linear, so this 2-mer is shared by
    // ALL reads entering the site (every allele passes through predPrev). This
    // is the symmetric counterpart of the trailing hom: it lets the interval's
    // upstream backbone anchor connect to a hom (which the backbone read spans)
    // instead of directly to a minority allele arm (which the backbone k=50
    // anchor does not share reads with). predBackboneOffset is predPrev's
    // backbone offset.
    int predBackboneOffset = -1;
    uint8_t predPrevBase = 0;   // base at predPrev (0-3)
    vector<HetAlleleMember> leadHomMembers;

    // Trailing hom anchor members: reads at commonSucc (rawPosition = succ base
    // position), forming the k=2 hom anchor [succBase, nextBase]. The
    // flank-linearity test guarantees commonSucc->succNext is linear, so this
    // 2-mer is shared by ALL reads spanning the site (all alleles reconverge
    // here). succBackboneOffset is the backbone offset of commonSucc.
    int succBackboneOffset = -1;
    uint8_t succBase = 0;   // base at commonSucc (0-3)
    vector<HetAlleleMember> homMembers;
};

// popcount over a read_ids bitset of `nWords` uint64 words.
inline int bitsetPopcount(const uint64_t* words, int nWords) {
    int c = 0;
    for(int w = 0; w < nWords; w++) c += __builtin_popcountll(words[w]);
    return c;
}

// Collect the read ids set in a bitset as a list of abPOA seq ids.
inline void bitsetToSeqIds(const uint64_t* words, int nWords, vector<int>& out) {
    out.clear();
    for(int w = 0; w < nWords; w++) {
        uint64_t bits = words[w];
        while(bits) {
            const int b = __builtin_ctzll(bits);
            out.push_back(w * 64 + b);
            bits &= bits - 1;
        }
    }
}

// abPOA base code (0-3) -> ASCII. Anything else -> 'N'.
inline char codeToChar(uint8_t c) {
    static const char lut[4] = {'A', 'C', 'G', 'T'};
    return (c < 4) ? lut[c] : 'N';
}

// Detect clean biallelic SNP bubbles in a finished window POA graph.
//
// A clean SNP bubble is a set of aligned nodes (X plus X.aligned_node_id[])
// that occupy the same MSA column but carry different bases, where every allele
// node has exactly one in-edge and one out-edge, and all allele nodes share a
// single common predecessor and a single common successor. This is the classic
// substitution bubble: one source -> {allele nodes} -> one sink.
//
// The backbone (read id 0) allele is the reference. Support for each allele is
// the popcount of the reads passing through the allele node (its out-edge
// read_ids). Per-allele read sets are mapped back to OrientedReadIds through
// seqIdToOrientedRead for downstream phasing.
//
// Filters: distinct non-gap bases, backbone present in the bubble, BOTH the
// reference and the alt allele supported by >= minSupport reads (a het site
// needs both alleles present, not just a strong alt vs the lone backbone read),
// and VAF >= minVaf (VAF = altSupport / spanning).
vector<WindowSnp> detectWindowSnps(
    abpoa_t* ab,
    const vector<int>& backboneQposToNode,
    const vector<uint8_t>& backboneCodes,
    const vector<OrientedReadId>& seqIdToOrientedRead,
    const vector<unordered_map<int, uint32_t>>& seqIdNodeToReadPos,
    int minSupport,
    double minVaf,
    bool dropHomopolymer,
    bool dropRepeat)
{
    vector<WindowSnp> snps;
    int droppedRepeat = 0;  // SNPs suppressed by homopolymer/STR context
    const abpoa_graph_t* abg = ab->abg;
    if(abg == nullptr || abg->node == nullptr) {
        return snps;
    }
    const int nodeN = abg->node_n;
    const int backboneLen = static_cast<int>(backboneCodes.size());

    // Reverse map: backbone node id -> backbone base offset. Node ids are dense
    // and small; use a flat vector sized to nodeN.
    vector<int> nodeToBackboneOffset(nodeN, -1);
    for(int off = 0; off < backboneLen; off++) {
        const int nid = backboneQposToNode[off];
        if(nid >= 0 && nid < nodeN) {
            nodeToBackboneOffset[nid] = off;
        }
    }

    vector<char> visited(nodeN, 0);
    vector<int> tmpSeqIds;

    for(int nid = 0; nid < nodeN; nid++) {
        if(nid == ABPOA_SRC_NODE_ID || nid == ABPOA_SINK_NODE_ID) continue;
        if(visited[nid]) continue;

        const abpoa_node_t& node = abg->node[nid];

        // Gather the allele group: this node plus its aligned nodes.
        vector<int> group;
        group.push_back(nid);
        for(int a = 0; a < node.aligned_node_n; a++) {
            group.push_back(node.aligned_node_id[a]);
        }
        // Mark the whole group visited regardless of outcome.
        for(int g : group) if(g >= 0 && g < nodeN) visited[g] = 1;

        if(group.size() < 2) continue;  // no variation at this column

        // Bubble test: every allele node has exactly one in-edge and one
        // out-edge, and all alleles share a single common predecessor and a
        // single common successor. A coexisting DELETION is allowed: it appears
        // as a skip edge pred->succ (no node), which does not disturb the allele
        // nodes. It is NOT a disqualifier -- the ref/alt substitution is still a
        // real SNP; the deletion reads are simply another allele that must be
        // counted into the spanning depth so VAF is not inflated.
        bool clean = true;
        int commonPred = -1, commonSucc = -1;
        for(int g : group) {
            if(g < 0 || g >= nodeN) { clean = false; break; }
            const abpoa_node_t& an = abg->node[g];
            if(an.in_edge_n != 1 || an.out_edge_n != 1) { clean = false; break; }
            const int pred = an.in_id[0];
            const int succ = an.out_id[0];
            if(commonPred == -1) { commonPred = pred; commonSucc = succ; }
            else if(pred != commonPred || succ != commonSucc) { clean = false; break; }
        }
        if(!clean) continue;

        // Purity test: the ONLY paths through the bubble may be the allele nodes
        // plus at most one deletion skip edge (pred->succ). Any predecessor
        // out-edge or successor in-edge that lands somewhere else means the
        // bubble is genuinely tangled (branch hub, overlapping indel), which we
        // do not report as a SNP. Build the allowed-neighbour set once.
        vector<char> isAllele(nodeN, 0);
        for(int g : group) isAllele[g] = 1;
        const abpoa_node_t& predNode = abg->node[commonPred];
        const abpoa_node_t& succNode = abg->node[commonSucc];

        // Deletion support = reads on the pred->succ skip edge (if present).
        int delSupport = 0;
        bool tangled = false;
        for(int e = 0; e < predNode.out_edge_n; e++) {
            const int tgt = predNode.out_id[e];
            if(tgt == commonSucc) {
                if(predNode.read_ids != nullptr && predNode.read_ids[e] != nullptr)
                    delSupport = bitsetPopcount(predNode.read_ids[e], predNode.read_ids_n);
            } else if(tgt < 0 || tgt >= nodeN || !isAllele[tgt]) {
                tangled = true; break;  // pred branches outside the bubble
            }
        }
        if(tangled) continue;
        for(int e = 0; e < succNode.in_edge_n; e++) {
            const int src = succNode.in_id[e];
            if(src == commonPred) continue;               // deletion skip edge
            if(src < 0 || src >= nodeN || !isAllele[src]) { tangled = true; break; }
        }
        if(tangled) continue;

        // Flank-linearity test. Require the two bases flanking the bubble to be
        // unambiguous (homozygous) so the k=2 het/hom anchors are well-founded:
        //   - the predecessor's own predecessor (predPrev) exists, connects ONLY
        //     to commonPred, and commonPred has a single in-edge. This makes the
        //     het anchor 2-mer [predBase, alleleBase] shared by all reads on the
        //     allele (predBase is linear).
        //   - the successor's own successor (succNext) exists, receives ONLY from
        //     commonSucc, and commonSucc has a single out-edge. This makes the
        //     hom separator 2-mer [succBase, nextBase] shared by all reads
        //     (succBase->nextBase is linear), guaranteeing a clean separator can
        //     be built between this SNP and the next one in the interval.
        // Together these also force >=2 linear bases between any two accepted
        // SNPs, so accepted bubbles are never adjacent and always chainable with
        // a hom in between. SRC(0)/SINK(1) are not valid flanks.
        if(commonPred < 2 || commonSucc < 2) continue;
        if(predNode.in_edge_n != 1) continue;
        const int predPrev = predNode.in_id[0];
        if(predPrev < 2 || predPrev >= nodeN) continue;
        if(abg->node[predPrev].out_edge_n != 1) continue;
        if(abg->node[predPrev].out_id[0] != commonPred) continue;
        if(succNode.out_edge_n != 1) continue;
        const int succNext = succNode.out_id[0];
        if(succNext < 2 || succNext >= nodeN) continue;
        if(abg->node[succNext].in_edge_n != 1) continue;
        if(abg->node[succNext].in_id[0] != commonSucc) continue;

        // Distinct non-gap bases (POA nodes always carry a real base 0-3).
        // Identify the backbone (reference) allele within the group.
        int refNode = -1;
        for(int g : group) {
            if(nodeToBackboneOffset[g] != -1) { refNode = g; break; }
        }
        if(refNode == -1) continue;  // backbone not in bubble -> not a ref SNP

        const uint8_t refCode = abg->node[refNode].base;

        // Per-allele support = popcount of the allele node's out-edge read_ids.
        // With out_edge_n==1 there is exactly one bitset per allele node.
        auto alleleSupport = [&](int g, vector<int>& seqIdsOut) -> int {
            const abpoa_node_t& an = abg->node[g];
            if(an.out_edge_n < 1 || an.read_ids == nullptr || an.read_ids[0] == nullptr) {
                seqIdsOut.clear();
                return 0;
            }
            bitsetToSeqIds(an.read_ids[0], an.read_ids_n, seqIdsOut);
            return static_cast<int>(seqIdsOut.size());
        };

        vector<int> refSeqIds;
        const int refSupport = alleleSupport(refNode, refSeqIds);

        // Collect EVERY distinct non-reference base allele. All alleles at this
        // column (ref + each alt) are emitted together in a SINGLE WindowSnp,
        // so a multiallelic site (ref + >=2 alt bases) becomes one N-armed
        // bubble rather than several biallelic records that duplicate the
        // reference. The total base-allele depth sums all alleles so per-alt
        // VAF (alt.support / spanning) is correct.
        struct AltAllele { int node; uint8_t base; int support; vector<int> seqIds; };
        vector<AltAllele> alts;
        int totalAlleleSupport = 0;
        for(int g : group) {
            vector<int> seqIds;
            const int s = alleleSupport(g, seqIds);
            totalAlleleSupport += s;
            if(g == refNode) continue;
            const uint8_t altCode = abg->node[g].base;
            if(altCode == refCode) continue;  // same base, not a substitution
            alts.push_back({g, altCode, s, std::move(seqIds)});
        }
        if(alts.empty()) continue;  // no distinct alt base

        // Spanning depth folds the deletion reads into the denominator so a
        // ref/alt/del site reports true VAF, not an inflated one. Shared across
        // all alt alleles at this column.
        const int spanning = totalAlleleSupport + delSupport;

        // Homopolymer / short-tandem-repeat context of the backbone at this
        // SNP. kmIsHomopolymer scans repeat units of length 1..6 requiring >=3
        // copies on either flank, so it covers both homopolymers (unit=1) and
        // STRs (unit 2..6). backboneCodes is already numeric (0-3), the exact
        // format the function expects. kmIsRepeatRegion is indel-only and a
        // no-op for SNPs, so it is intentionally not called here.
        // Homopolymer/STR context depends only on the backbone position, so it
        // is the same for every alt allele at this column.
        const int backboneOff = nodeToBackboneOffset[refNode];
        KmVarKey vkey;
        vkey.pos = static_cast<uint32_t>(backboneOff);
        vkey.type = KmVarType::Snp;
        vkey.altBase = alts[0].base;
        vkey.refLen = 1;
        vkey.altLen = 1;
        // Split the repeat-context test into two independent classes so each
        // can be filtered separately: homopolymer = repeat unit length 1,
        // short tandem repeat (STR) = repeat unit length 2..6.
        const uint8_t* bcData = backboneCodes.data();
        const uint32_t bcLen = static_cast<uint32_t>(backboneCodes.size());
        const bool inHomopolymer = kmIsRepeatUnitRange(bcData, bcLen, vkey, 0, 1, 1);
        const bool inStr         = kmIsRepeatUnitRange(bcData, bcLen, vkey, 0, 2, 6);
        const bool inRepeat = inHomopolymer || inStr;
        if(inRepeat) droppedRepeat++;

        // Verification dump (DINARA_SNP_VERIFY=1): print the backbone context
        // around each SNP and an INDEPENDENT recomputation of the repeat rule,
        // flagging any disagreement with kmIsHomopolymer. The independent check
        // mirrors the documented rule (unit len 1..6, >=3 copies on either
        // flank) but is written from scratch so a bug in one won't hide in the
        // other.
        if(getenv("DINARA_SNP_VERIFY")) {
            const int n = static_cast<int>(backboneCodes.size());
            auto at = [&](int p) -> int {
                return (p >= 0 && p < n) ? int(backboneCodes[p]) : -1;
            };
            // Independent repeat test around the SNP at backboneOff.
            // Forward flank starts at backboneOff+1, backward at backboneOff-1
            // (SNP base excluded, matching kmIsHomopolymer's startPos/endPos).
            auto indepRepeat = [&]() -> bool {
                for(int r = 1; r <= 6; r++) {
                    // forward: bases [off+1 .. off+r] repeated >=3x
                    bool fwd = true;
                    for(int c = 1; c < 3 && fwd; c++)
                        for(int j = 0; j < r && fwd; j++) {
                            const int a = at(backboneOff + 1 + j);
                            const int b = at(backboneOff + 1 + c*r + j);
                            if(a < 0 || a != b) fwd = false;
                        }
                    if(fwd) return true;
                    // backward: bases [off-1 .. off-r] repeated >=3x
                    bool bwd = true;
                    for(int c = 1; c < 3 && bwd; c++)
                        for(int j = 0; j < r && bwd; j++) {
                            const int a = at(backboneOff - 1 - j);
                            const int b = at(backboneOff - 1 - c*r - j);
                            if(a < 0 || a != b) bwd = false;
                        }
                    if(bwd) return true;
                }
                return false;
            };
            const bool indep = indepRepeat();
            string ctx;
            for(int p = backboneOff - 18; p <= backboneOff + 18; p++) {
                const int b = at(p);
                char c = (b < 0) ? '.' : codeToChar(uint8_t(b));
                if(p == backboneOff) { ctx += '['; ctx += c; ctx += ']'; }
                else ctx += c;
            }
            // Bubble purity: how many paths actually enter/leave the bubble.
            // For a truly biallelic SNP, pred.out_edge_n and succ.in_edge_n
            // should equal the number of allele nodes. Any excess means an
            // extra path (typically a deletion skip edge) coexists at this
            // column, so the site is really multi-allelic.
            const int nAlleles = static_cast<int>(group.size());
            const int predOut = abg->node[commonPred].out_edge_n;
            const int succIn  = abg->node[commonSucc].in_edge_n;
            const bool extraPath = (predOut > nAlleles) || (succIn > nAlleles);
            cout << "    VERIFY off=" << backboneOff
                 << " " << ctx
                 << " km=" << (inRepeat ? 1 : 0)
                 << " indep=" << (indep ? 1 : 0)
                 << " alleles=" << nAlleles
                 << " predOut=" << predOut
                 << " succIn=" << succIn
                 << (extraPath ? "  <<< EXTRA_PATH(del?)" : "")
                 << (inRepeat != indep ? "  <<< MISMATCH" : "") << endl;
        }

        // Drop by class: homopolymer and STR context are gated independently.
        if(dropHomopolymer && inHomopolymer) continue;
        if(dropRepeat && inStr) continue;

        auto mapReads = [&](const vector<int>& seqIds, vector<OrientedReadId>& out) {
            out.clear();
            for(int sid : seqIds) {
                if(sid >= 0 && sid < static_cast<int>(seqIdToOrientedRead.size())) {
                    out.push_back(seqIdToOrientedRead[sid]);
                }
            }
        };

        // Build het-anchor members for an allele: for each supporting read,
        // recover its base position at the bubble's common predecessor node
        // (the k=2 marker's rawPosition = predBase position). Reads whose
        // predecessor position is not recoverable from the POA are skipped
        // (they still count in support; they just can't seed an anchor member).
        auto mapMembersAt = [&](const vector<int>& seqIds, int atNode,
                                vector<HetAlleleMember>& out) {
            out.clear();
            for(int sid : seqIds) {
                if(sid < 0 || sid >= static_cast<int>(seqIdToOrientedRead.size())) continue;
                if(sid >= static_cast<int>(seqIdNodeToReadPos.size())) continue;
                const auto& nodeMap = seqIdNodeToReadPos[sid];
                const auto it = nodeMap.find(atNode);
                if(it == nodeMap.end()) continue;  // position not recoverable
                out.push_back({seqIdToOrientedRead[sid], it->second});
            }
        };
        auto mapMembers = [&](const vector<int>& seqIds, vector<HetAlleleMember>& out) {
            mapMembersAt(seqIds, commonPred, out);
        };

        // Emit ONE record per column, carrying every allele that clears the
        // gates. For a het site the reference allele and at least one alternate
        // must be independently supported (>= minSupport each; each alt also
        // >= minVaf). Gating only the alt would let through
        // homozygous-alt-vs-backbone columns (refSupport == 1, the lone
        // backbone read) and single-read artifacts, which are not het sites.
        // A multiallelic column (ref + >=2 strong alts) becomes a single
        // N-armed bubble here rather than several biallelic bubbles that would
        // duplicate the reference and be wrongly chained in series.
        if(refSupport < minSupport) continue;

        // Collect the qualifying alt alleles first, so we know the arm count
        // and can skip columns that end up with no usable alt.
        vector<const AltAllele*> keptAlts;
        for(const AltAllele& alt : alts) {
            if(alt.support < minSupport) continue;
            const double vaf = (spanning > 0) ? double(alt.support) / double(spanning) : 0.0;
            if(vaf < minVaf) continue;
            keptAlts.push_back(&alt);
        }
        if(keptAlts.empty()) continue;   // no strong alt: not a het site

        WindowSnp snp;
        snp.msaRank = (abg->node_id_to_msa_rank != nullptr) ? abg->node_id_to_msa_rank[refNode] : -1;
        snp.backboneOffset = backboneOff;
        snp.inHomopolymerOrRepeat = inRepeat;
        snp.delSupport = delSupport;
        snp.spanning = spanning;

        // Reference allele (arm 0).
        {
            WindowSnpAllele refAllele;
            refAllele.base = codeToChar(refCode);
            refAllele.code = refCode;
            refAllele.support = refSupport;
            refAllele.isRef = true;
            mapReads(refSeqIds, refAllele.reads);
            mapMembers(refSeqIds, refAllele.members);
            snp.alleles.push_back(std::move(refAllele));
        }
        // Alternate alleles (arms 1..).
        for(const AltAllele* alt : keptAlts) {
            WindowSnpAllele altAllele;
            altAllele.base = codeToChar(alt->base);
            altAllele.code = alt->base;
            altAllele.support = alt->support;
            altAllele.isRef = false;
            mapReads(alt->seqIds, altAllele.reads);
            mapMembers(alt->seqIds, altAllele.members);
            snp.alleles.push_back(std::move(altAllele));
        }

        // All reads spanning the site (ref + every kept alt) share both flanks:
        // they enter through predPrev and leave through commonSucc (both linear
        // by the flank-linearity test). Build this shared set once for the two
        // bracketing hom anchors.
        vector<int> homSeqIds = refSeqIds;
        for(const AltAllele* alt : keptAlts)
            homSeqIds.insert(homSeqIds.end(), alt->seqIds.begin(), alt->seqIds.end());

        // Leading hom anchor at predPrev: [predPrevBase, predBase]. predPrev is
        // the base BEFORE commonPred; predPrev->commonPred is linear, so this
        // 2-mer is shared by ALL reads entering the site. This lets the upstream
        // backbone anchor connect to a hom instead of directly to a minority
        // allele arm (whose reads the k=50 backbone anchor does not share).
        snp.predBackboneOffset = nodeToBackboneOffset[predPrev];
        snp.predPrevBase = abg->node[predPrev].base;
        mapMembersAt(homSeqIds, predPrev, snp.leadHomMembers);

        // Trailing hom anchor at commonSucc: [succBase, nextBase]. All reads
        // spanning the site reconverge here, so recover each one's position at
        // commonSucc. The flank-linearity test guarantees the next base is
        // linear, so [succBase, nextBase] is shared by all.
        snp.succBackboneOffset = nodeToBackboneOffset[commonSucc];
        snp.succBase = abg->node[commonSucc].base;
        mapMembersAt(homSeqIds, commonSucc, snp.homMembers);

        snps.push_back(std::move(snp));
    }

    // Report in backbone order.
    sort(snps.begin(), snps.end(), [](const WindowSnp& a, const WindowSnp& b) {
        return a.backboneOffset < b.backboneOffset;
    });
    cout << "  homopolymer/STR SNPs "
         << ((dropHomopolymer || dropRepeat) ? "dropped(some): " : "flagged: ")
         << droppedRepeat << endl;
    return snps;
}

} // anonymous namespace


// Build one abPOA multi-segment MSA for a single anchor window using all
// oriented reads that share at least two of the window's backbone anchors.
// Returns true if an MSA was produced, false if the window was skipped.
bool Assembler::runOneWindowAbpoaMultiSegmentMSA(
    const shared_ptr<Shasta2Anchors>& shasta2Anchors,
    const shared_ptr<Shasta2Journeys>& shasta2Journeys,
    AnchorWindow& window,
    std::ostream& out,
    double hetMinVaf,
    uint64_t hetMinSupport,
    bool hetDropHomopolymer,
    bool hetDropRepeat)
{
    const Reads& readsRef = getReads();
    const auto& markersRef = *markers;
    const uint64_t k = assemblerInfo->k;

    const OrientedReadId backboneOid = window.backboneOrientedReadId;
    const auto backboneJourney = (*shasta2Journeys)[backboneOid];

    out << "testAbpoaMultiSegmentMSA: window " << window.windowId
         << " backbone " << backboneOid
         << " anchors [" << window.backboneBegin << "," << window.backboneEnd << ")"
         << " reads " << window.readIntervals.size() << endl;

    const uint32_t nBackboneAnchors = window.backboneEnd - window.backboneBegin;
    if(nBackboneAnchors < 2) {
        out << "  window " << window.windowId << " has < 2 anchors, skipping." << endl;
        return false;
    }
    const uint32_t nSegments = nBackboneAnchors - 1;

    // Backbone marker ordinal at each anchor boundary (0..nSegments).
    // anchorOrdinal[bi] is the backbone read's marker ordinal for the anchor at
    // journey position backboneBegin + bi.
    vector<uint32_t> anchorOrdinal(nBackboneAnchors);
    for(uint32_t bi = 0; bi < nBackboneAnchors; bi++) {
        const Shasta2AnchorId anchorId = backboneJourney[window.backboneBegin + bi];
        anchorOrdinal[bi] = shasta2Anchors->getOrdinal(anchorId, backboneOid);
    }

    // Full backbone base sequence as 0-3 codes, plus the base offset of each
    // anchor's MIDPOINT within that sequence.
    //
    // FULL-ANCHOR COVERAGE: anchors are k=50-base markers; using only their
    // midpoint (position + kHalf) as the segment boundary leaves the OUTER kHalf
    // of the two window-terminal anchors with no column at all. To cover the
    // full 50 bases of every anchor, we extend the backbone span by kHalf on
    // each end: it now runs from the START of the first anchor
    // (marker[first].position) to the END of the last anchor
    // (marker[last].position + k), instead of midpoint-to-midpoint. This adds
    // kHalf real columns at each window edge so reads can place the outer halves
    // of their first/last anchors.
    const uint32_t kHalf = uint32_t(k / 2);
    const auto backboneMarkers = markersRef[backboneOid.getValue()];
    const uint32_t backboneBeginPos = backboneMarkers[anchorOrdinal.front()].position;          // first anchor START
    const uint32_t backboneEndPos   = backboneMarkers[anchorOrdinal.back()].position + uint32_t(k); // last anchor END
    if(backboneEndPos <= backboneBeginPos) {
        out << "  backbone span empty, skipping window." << endl;
        return false;
    }
    vector<uint8_t> backboneCodes;
    backboneCodes.reserve(backboneEndPos - backboneBeginPos);
    for(uint32_t pos = backboneBeginPos; pos < backboneEndPos; pos++) {
        backboneCodes.push_back(readsRef.getOrientedReadBase(backboneOid, pos).value);
    }
    const int backboneLen = static_cast<int>(backboneCodes.size());

    // anchorOffset[bi] = base offset of anchor bi's MIDPOINT into backboneCodes.
    // With the kHalf-extended span, anchorOffset[0] == kHalf (not 0) and
    // anchorOffset[last] == backboneLen - kHalf (not backboneLen).
    vector<int> anchorOffset(nBackboneAnchors);
    for(uint32_t bi = 0; bi < nBackboneAnchors; bi++) {
        const uint32_t pos = backboneMarkers[anchorOrdinal[bi]].position + kHalf;  // midpoint
        anchorOffset[bi] = static_cast<int>(pos - backboneBeginPos);
    }

    out << "  backbone " << backboneLen << " bases across "
         << nSegments << " segments (full-anchor span, +" << kHalf
         << "bp each end)" << endl;

    // ------------------------------------------------------------------
    // abPOA setup.
    // Scoring chosen to match the theseus path's intent: match=0 implied by
    // POA (abPOA uses positive match), 2/3/1 style mismatch+gap. abPOA needs a
    // positive match score; we use the abPOA defaults (match=2, mismatch=4,
    // affine gap 4/2) which are well tuned for long-read POA. Banding (wb=10)
    // is what keeps each piece alignment O(qlen*band).
    // ------------------------------------------------------------------
    abpoa_t* ab = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    abpt->align_mode = ABPOA_GLOBAL_MODE;
    abpt->gap_mode = ABPOA_AFFINE_GAP;
    abpt->match = 2;
    abpt->mismatch = 4;
    abpt->gap_open1 = 4;
    abpt->gap_ext1 = 2;
    abpt->gap_open2 = 0;       // affine (single-piece) gap
    abpt->gap_ext2 = 0;
    abpt->wb = 10;             // adaptive band; <0 disables banding
    abpt->wf = 0.01;
    abpt->disable_seeding = 1; // we drive the segmentation ourselves
    abpt->progressive_poa = 0;
    abpt->out_msa = 1;         // need RC-MSA output
    abpt->out_cons = 0;
    abpt->ret_cigar = 1;
    abpoa_post_set_para(abpt);  // sets use_read_ids etc. from out_msa

    // Total number of sequences = backbone + all reads sharing >=2 anchors.
    // We need an upper bound up front for read_id bitsets; recount precisely
    // after building boundary hits. Use readIntervals.size() as a safe bound.
    const int totReadBound = static_cast<int>(window.readIntervals.size()) + 1;

    // Seed the backbone as read_id 0. On an empty graph (node_n==2),
    // abpoa_add_subgraph_alignment lays the sequence down as a linear chain and
    // fills qpos_to_node_id with the node id for each base position.
    vector<int> backboneQposToNode(backboneLen, -1);
    {
        abpoa_res_t res;
        res.n_cigar = 0; res.m_cigar = 0; res.graph_cigar = nullptr;
        abpoa_add_subgraph_alignment(
            ab, abpt,
            ABPOA_SRC_NODE_ID, ABPOA_SINK_NODE_ID,
            backboneCodes.data(), nullptr, backboneLen,
            backboneQposToNode.data(), res,
            /* read_id */ 0, /* tot_read_n */ totReadBound,
            /* inc_both_ends */ 1);
    }

    // anchorNode[bi] = abPOA node id of anchor bi's MIDPOINT base.
    // With the kHalf-extended backbone, every anchor midpoint offset lies in
    // [kHalf, backboneLen-kHalf], so ALL anchors map to real interior nodes and
    // are handled uniformly as anchor MATCHes (no SINK special case).
    // backboneNodeAt(off) resolves a backbone base offset to its node id, or to
    // SRC/SINK when the offset falls just before/after the chain (used for the
    // excluded fold boundaries of a read's outer anchor halves).
    vector<int> anchorNode(nBackboneAnchors);
    for(uint32_t bi = 0; bi < nBackboneAnchors; bi++) {
        const int off = anchorOffset[bi];
        anchorNode[bi] = (off >= 0 && off < backboneLen)
            ? backboneQposToNode[off]
            : ABPOA_SINK_NODE_ID;  // defensive; should not occur with extension
    }
    auto backboneNodeAt = [&](int off) -> int {
        if(off < 0) return ABPOA_SRC_NODE_ID;
        if(off >= backboneLen) return ABPOA_SINK_NODE_ID;
        return backboneQposToNode[off];
    };

    // ------------------------------------------------------------------
    // Build read -> backbone boundary hits (identical logic to the theseus
    // path: anchor membership, sort, drop <2, clip to pairwise alignment
    // ordinal range, re-drop <2, sort by base span descending).
    // ------------------------------------------------------------------
    struct BoundaryHit {
        uint32_t boundaryIndex;  // anchor index 0..nSegments
        uint32_t ordinal;        // marker ordinal on this read
    };
    unordered_map<uint64_t, vector<BoundaryHit>> readBoundaryHits;

    for(uint32_t bi = 0; bi <= nSegments; bi++) {
        const uint32_t bjp = window.backboneBegin + bi;
        if(bjp >= backboneJourney.size()) break;
        const Shasta2AnchorId anchorId = backboneJourney[bjp];
        const auto anchor = (*shasta2Anchors)[anchorId];
        for(const auto& info : anchor) {
            if(info.orientedReadId == backboneOid) continue;
            readBoundaryHits[info.orientedReadId.getValue()].push_back(
                {bi, info.ordinal});
        }
    }
    for(auto& [readId, hits] : readBoundaryHits) {
        sort(hits.begin(), hits.end(),
            [](const BoundaryHit& a, const BoundaryHit& b) {
                return a.boundaryIndex < b.boundaryIndex;
            });
    }
    uint32_t skippedReads = 0;
    for(auto it = readBoundaryHits.begin(); it != readBoundaryHits.end(); ) {
        if(it->second.size() < 2) { skippedReads++; it = readBoundaryHits.erase(it); }
        else ++it;
    }

    // Clip to pairwise alignment ordinal range.
    const auto& clipTable = getAlignmentTable();
    const ReadId backboneReadId = backboneOid.getReadId();
    for(auto& [readIdValue, hits] : readBoundaryHits) {
        const OrientedReadId oid = OrientedReadId::fromValue(static_cast<ReadId>(readIdValue));
        const ReadId readId = oid.getReadId();
        uint32_t bestFirst = 0, bestLast = 0, bestSpan = 0;
        const auto& aligns = clipTable[oid.getValue()];
        for(uint32_t idx : aligns) {
            const auto& ad = alignmentData[idx];
            ReadId partnerId = (ad.readIds[0] == readId) ? ad.readIds[1] : ad.readIds[0];
            if(partnerId != backboneReadId) continue;
            int targetIdx = (ad.readIds[0] == readId) ? 0 : 1;
            uint32_t firstOrd = ad.info.data[targetIdx].firstOrdinal;
            uint32_t lastOrd  = ad.info.data[targetIdx].lastOrdinal;
            Strand storedStrand = (targetIdx == 0) ? 0 : (ad.isSameStrand ? 0 : 1);
            if(storedStrand != oid.getStrand()) {
                uint32_t mc = ad.info.data[targetIdx].markerCount;
                uint32_t f = mc - 1 - lastOrd;
                uint32_t l = mc - 1 - firstOrd;
                firstOrd = f; lastOrd = l;
            }
            uint32_t span = (lastOrd > firstOrd) ? (lastOrd - firstOrd) : 0;
            if(span > bestSpan) { bestSpan = span; bestFirst = firstOrd; bestLast = lastOrd; }
        }
        if(bestSpan > 0) {
            hits.erase(remove_if(hits.begin(), hits.end(),
                [&](const BoundaryHit& h) {
                    return h.ordinal < bestFirst || h.ordinal > bestLast;
                }), hits.end());
        }
    }
    for(auto it = readBoundaryHits.begin(); it != readBoundaryHits.end(); ) {
        if(it->second.size() < 2) { skippedReads++; it = readBoundaryHits.erase(it); }
        else ++it;
    }

    // Sort reads by base span descending (longest first).
    vector<pair<uint32_t, uint64_t>> readsBySpan;
    readsBySpan.reserve(readBoundaryHits.size());
    for(const auto& [readIdValue, hits] : readBoundaryHits) {
        const auto readMarkers = markersRef[readIdValue];
        uint32_t firstPos = readMarkers[hits.front().ordinal].position;
        uint32_t lastPos  = readMarkers[hits.back().ordinal].position;
        uint32_t span = (lastPos > firstPos) ? (lastPos - firstPos) : 0;
        readsBySpan.push_back({span, readIdValue});
    }
    sort(readsBySpan.begin(), readsBySpan.end(),
        [](const auto& a, const auto& b) { return a.first > b.first; });

    out << "  reads with >=2 shared anchors: " << readsBySpan.size()
         << ", skipped: " << skippedReads << endl;

    // ------------------------------------------------------------------
    // Per-read alignment: fold each read piece (between consecutive shared
    // anchors) into the graph via abPOA subgraph alignment. Each read gets a
    // single read_id; all its pieces accumulate cigars into one whole_res that
    // is added once, so the MSA row is contiguous.
    //
    // abPOA's subgraph API takes EXCLUDED begin/end nodes:
    //   exc_beg | inc_beg ... inc_end | exc_end
    // For a piece spanning anchors a -> b on the backbone, the bases strictly
    // between anchor a and anchor b are aligned. We exclude the node at anchor a
    // (the left pin) and the node at anchor b (the right pin), so beg = node at
    // anchor a, end = node at anchor b.
    // ------------------------------------------------------------------
    uint32_t alignedSegments = 0;
    uint32_t alignedReads = 0;
    double totalAlignTime = 0.0;
    double maxAlignTime = 0.0;
    size_t totalAlignBases = 0;

    vector<double> perReadTime;
    int readSeqId = 1;  // 0 is the backbone

    // Map abPOA read id (readSeqId) -> OrientedReadId, for SNP allele phasing.
    // Index 0 is the backbone; entries 1..alignedReads are filled as reads fold.
    vector<OrientedReadId> seqIdToOrientedRead;
    seqIdToOrientedRead.push_back(backboneOid);  // read id 0 = backbone

    // Per-read node->absolute-base-position map, indexed by abPOA read id.
    // For each read we capture abPOA's qpos_to_node_id (queryPos -> nodeId) and
    // invert it to nodeId -> (read's absolute base position). This lets het
    // anchor generation recover, for any read on an allele, the read base
    // position at the bubble's predecessor node (the k=2 anchor's rawPosition).
    // seqId 0 = backbone: its node<->offset relation is backboneQposToNode, and
    // the absolute read position is backboneBeginPos + backbone offset.
    vector<unordered_map<int, uint32_t>> seqIdNodeToReadPos;
    seqIdNodeToReadPos.emplace_back();  // seqId 0 = backbone
    {
        auto& bbMap = seqIdNodeToReadPos.back();
        bbMap.reserve(backboneLen);
        for(int off = 0; off < backboneLen; off++) {
            const int nid = backboneQposToNode[off];
            if(nid >= 0) bbMap[nid] = backboneBeginPos + uint32_t(off);
        }
    }

    // abPOA re-runs a full topological sort of the WHOLE graph on every
    // add_subgraph_alignment. Adding each piece separately therefore costs
    // O(graph_nodes) per piece -> O(graph_nodes * pieces), which dominates
    // runtime even though each piece ALIGNMENT is fast. abPOA's own internal
    // aligner avoids this by accumulating all of a read's piece cigars into one
    // whole-read cigar and folding it in with a SINGLE add per read (one sort
    // per read, not per piece). We do the same here.
    //
    // FULL-ANCHOR COVERAGE. Every shared anchor is a k=50-base marker; its
    // MIDPOINT (position + kHalf) is the segment boundary. The read's query span
    // covers the FULL extent of its first and last shared anchors: from the
    // START of the first anchor (rm[firstOrd].position) to the END of the last
    // anchor (rm[lastOrd].position + k). All anchors -- first, interior, last --
    // are emitted as a single-base midpoint MATCH; the outer kHalf of the first
    // and last anchors are aligned as leading/trailing gaps against the
    // kHalf-extended backbone. The stitched cigar is:
    //   <leading gap>          <- outer kHalf of first anchor (before its midpoint)
    //   MATCH(anchorNode[a0])  <- first anchor midpoint
    //   <gap 0>                <- bases between a0 and a1 midpoints
    //   MATCH(anchorNode[a1])
    //   ...
    //   MATCH(anchorNode[am])  <- last anchor midpoint
    //   <trailing gap>         <- outer kHalf of last anchor (after its midpoint)
    //
    // The fold's excluded boundaries (foldBegNode/foldEndNode) are the backbone
    // nodes immediately OUTSIDE the read's covered span (one base before the
    // first covered base, one base after the last). Because those boundaries are
    // now genuinely outside the read, inc_both_ends is 0: membership is recorded
    // only on the read's real covered bases (abPOA records read membership on
    // the SOURCE node of each out-edge; the final synthesized edge records the
    // read's true last base). Each anchor MATCH consumes the read's OWN base, so
    // a sequencing error at an anchor spawns an aligned variant node rather than
    // silently displaying the backbone base. add_subgraph_alignment derives
    // query_id by counting M/I ops, so cigars must cover every query base
    // exactly once, in order.
    auto encodeMatch = [](int nodeId) -> abpoa_cigar_t {
        // CMATCH: node_id << 34 | query_id << 4 | op. query_id is recomputed by
        // the consumer for M ops, so we only need node_id and the op.
        return (static_cast<abpoa_cigar_t>(nodeId) << 34)
             | static_cast<abpoa_cigar_t>(ABPOA_CMATCH);
    };

    const int intK = static_cast<int>(k);

    // abPOA allocates the FULL DP matrix for every subgraph alignment
    // (simd_abpoa_realloc: sn * gn * 3 * size, where sn ~ qlen/lanes and gn is
    // the subgraph node span). The adaptive band (wb) limits COMPUTATION but
    // NOT allocation, so a read whose two shared anchors are far apart on the
    // backbone -- a chimera, a repeat mismap, or a large structural variant --
    // forces a matrix of tens of GiB and aborts the whole run via a
    // posix_memalign failure ("[SIMDMalloc] posix_memalign fail!"). Estimate
    // the matrix size before each gap alignment and drop the offending read
    // rather than let abPOA abort. Cap is per-alignment; override with
    // DINARA_ABPOA_MAX_DP_GIB (default 4 GiB, generous for legitimate gaps
    // which are at most a few thousand bases x a few thousand nodes).
    uint64_t maxDpGiB = 4;
    if(const char* e = getenv("DINARA_ABPOA_MAX_DP_GIB")) {
        char* endp = nullptr;
        const unsigned long v = strtoul(e, &endp, 10);
        if(endp != e && v > 0) maxDpGiB = v;
    }
    const uint64_t maxDpBytes = maxDpGiB << 30;
    uint32_t oversizeDroppedReads = 0;

    for(const auto& [baseSpan, readIdValue] : readsBySpan) {
        (void)baseSpan;
        const auto& hits = readBoundaryHits[readIdValue];
        const OrientedReadId oid = OrientedReadId::fromValue(static_cast<ReadId>(readIdValue));
        const auto rm = markersRef[oid.getValue()];

        const uint32_t firstBoundary = hits.front().boundaryIndex;
        const uint32_t lastBoundary  = hits.back().boundaryIndex;
        if(firstBoundary >= nBackboneAnchors || lastBoundary >= nBackboneAnchors) continue;

        // Query span in read coordinates: first-anchor START .. last-anchor END.
        const uint32_t qBegin = rm[hits.front().ordinal].position;          // first anchor start
        const uint32_t qEnd   = rm[hits.back().ordinal].position + intK;    // last anchor end
        if(qEnd <= qBegin) continue;
        const int qlen = static_cast<int>(qEnd - qBegin);

        // Fold boundaries: backbone nodes just OUTSIDE the read's covered span.
        // The read's first covered base aligns to backbone offset
        // (anchorOffset[firstBoundary] - kHalf); its predecessor is the excluded
        // left boundary. Symmetrically on the right. SRC/SINK at the window ends.
        const int firstMidOffset = anchorOffset[firstBoundary];
        const int lastMidOffset  = anchorOffset[lastBoundary];
        const int foldBegNode = backboneNodeAt(firstMidOffset - static_cast<int>(kHalf) - 1);
        const int foldEndNode = backboneNodeAt(lastMidOffset  + static_cast<int>(kHalf) + 1);

        vector<uint8_t> readCodes;
        readCodes.reserve(qlen);
        for(uint32_t pos = qBegin; pos < qEnd; pos++) {
            readCodes.push_back(readsRef.getOrientedReadBase(oid, pos).value);
        }

        vector<abpoa_cigar_t> whole;   // stitched whole-read cigar
        double readTime = 0.0;
        uint32_t readSegments = 0;
        bool readOk = true;

        // Helper to align a gap [begQpos, gapEnd) against subgraph (beg, end)
        // and append its cigars, updating timers.
        auto alignGap = [&](int begNode, int endNode, int begQpos, int gapLen,
                            uint32_t boundaryForLog) -> bool {
            if(gapLen <= 0) return true;
            // Estimate abPOA's DP allocation and skip this read if it would be
            // oversized. gn is the subgraph node span (begNode..endNode in
            // topological index order); the matrix is ~ (gapLen/lanes + 1) * gn
            // * 3 (affine) * 4 bytes (32-bit fallback, the larger case). Using
            // 4 bytes/lane and dropping the /lanes factor is a safe upper bound
            // that never under-counts, so we stay well clear of the real alloc.
            const int begIndex = ab->abg->node_id_to_index[begNode];
            const int endIndex = ab->abg->node_id_to_index[endNode];
            const uint64_t gn = (endIndex >= begIndex)
                ? uint64_t(endIndex - begIndex + 1) : 0;
            const uint64_t estBytes =
                (uint64_t(gapLen) + 1) * gn * 3 * 4;
            if(estBytes > maxDpBytes) {
                out << "  OVERSIZE: read " << oid
                    << " gap near anchor " << boundaryForLog
                    << " would need ~" << (estBytes >> 30)
                    << " GiB (" << gapLen << " bases x " << gn
                    << " nodes); dropping read" << endl;
                return false;
            }
            abpoa_res_t res;
            res.n_cigar = 0; res.m_cigar = 0; res.graph_cigar = nullptr;
            auto t0 = chrono::steady_clock::now();
            abpoa_align_sequence_to_subgraph(
                ab, abpt, begNode, endNode,
                readCodes.data() + begQpos, gapLen, &res);
            auto t1 = chrono::steady_clock::now();
            double elapsed = chrono::duration<double>(t1 - t0).count();
            for(int ci = 0; ci < res.n_cigar; ci++) whole.push_back(res.graph_cigar[ci]);
            if(res.n_cigar) free(res.graph_cigar);
            totalAlignTime += elapsed;
            readTime += elapsed;
            totalAlignBases += gapLen;
            if(elapsed > maxAlignTime) maxAlignTime = elapsed;
            if(elapsed > 0.1) {
                out << "  SLOW: read " << oid
                     << " gap near anchor " << boundaryForLog
                     << " seq " << gapLen << " bases took " << elapsed << "s" << endl;
            }
            return true;
        };

        int begNode = foldBegNode;   // excluded left boundary
        int begQpos = 0;             // query offset of next unaligned base
        for(size_t hi = 0; hi < hits.size(); hi++) {
            const uint32_t boundary = hits[hi].boundaryIndex;
            if(boundary >= nBackboneAnchors) { readOk = false; break; }
            const int pinNode = anchorNode[boundary];

            // Anchor midpoint offset in readCodes.
            const int pinQpos = static_cast<int>(rm[hits[hi].ordinal].position + kHalf) - static_cast<int>(qBegin);
            if(pinQpos < begQpos || pinQpos >= qlen) { readOk = false; break; }

            // Gap of read bases strictly before this anchor midpoint, aligned
            // against (begNode, pinNode). For hi==0 this is the outer kHalf of
            // the first anchor; for hi>0 it is the inter-anchor segment.
            if(!alignGap(begNode, pinNode, begQpos, pinQpos - begQpos, boundary)) {
                readOk = false;
                oversizeDroppedReads++;
                break;
            }

            // Anchor midpoint MATCH (read's own base).
            whole.push_back(encodeMatch(pinNode));
            begNode = pinNode;
            begQpos = pinQpos + 1;
            readSegments++;
            alignedSegments++;
        }

        if(readOk) {
            // Trailing gap: outer kHalf of the last anchor, aligned against
            // (lastAnchorNode, foldEndNode).
            if(!alignGap(begNode, foldEndNode, begQpos, qlen - begQpos, lastBoundary)) {
                readOk = false;
                oversizeDroppedReads++;
            }
        }

        if(readOk && readSegments > 0 && !whole.empty()) {
            abpoa_res_t wholeRes;
            wholeRes.n_cigar = static_cast<int>(whole.size());
            wholeRes.m_cigar = static_cast<int>(whole.size());
            wholeRes.graph_cigar = whole.data();
            // Capture this read's queryPos -> nodeId map so het anchor
            // generation can recover the read's base position at any bubble
            // node. qpos_to_node_id[q] is the node the read's q-th base landed
            // on (-1 for inserted/unaligned bases).
            vector<int> qposToNode(qlen, -1);
            auto t0 = chrono::steady_clock::now();
            abpoa_add_subgraph_alignment(
                ab, abpt,
                foldBegNode, foldEndNode,
                readCodes.data(), nullptr, qlen,
                qposToNode.data(), wholeRes,
                readSeqId, totReadBound, /* inc_both_ends */ 0);
            auto t1 = chrono::steady_clock::now();
            readTime += chrono::duration<double>(t1 - t0).count();
            totalAlignTime += chrono::duration<double>(t1 - t0).count();

            // Invert to nodeId -> absolute read base position (qBegin + q).
            seqIdNodeToReadPos.emplace_back();
            auto& nodeMap = seqIdNodeToReadPos.back();
            nodeMap.reserve(qlen);
            for(int q = 0; q < qlen; q++) {
                const int nid = qposToNode[q];
                if(nid >= 0) nodeMap[nid] = qBegin + uint32_t(q);
            }

            perReadTime.push_back(readTime);
            seqIdToOrientedRead.push_back(oid);  // read id readSeqId = this read
            alignedReads++;
            readSeqId++;
        }
    }

    out << "  aligned " << alignedReads << " reads ("
         << alignedSegments << " segments), skipped " << skippedReads;
    if(oversizeDroppedReads > 0)
        out << ", dropped " << oversizeDroppedReads << " oversize";
    out << endl;

    // Report timing by quartile (reads are already longest-first).
    if(perReadTime.size() >= 4) {
        size_t q = perReadTime.size() / 4;
        for(int qi = 0; qi < 4; qi++) {
            size_t start = qi * q;
            size_t end = (qi == 3) ? perReadTime.size() : (qi + 1) * q;
            double sum = 0;
            for(size_t i = start; i < end; i++) sum += perReadTime[i];
            out << "  Q" << (qi + 1) << " (reads " << start << "-" << (end - 1)
                 << "): " << sum << "s total, "
                 << (sum / (end - start) * 1000) << "ms/read" << endl;
        }
    }
    out << "  total align time: " << totalAlignTime << "s"
         << "  avg: " << (alignedSegments > 0 ? totalAlignTime / alignedSegments * 1000 : 0)
         << "ms/seg  max: " << maxAlignTime << "s"
         << "  total bases: " << totalAlignBases << endl;

    // ------------------------------------------------------------------
    // Emit the MSA. We built the graph manually, so abPOA's seq store does not
    // know how many rows exist; set n_seq = backbone + aligned reads. The RC-MSA
    // generator reads per-edge read_id bitsets (populated by add_read_id, on
    // because out_msa=1) to lay each row out. is_rc[]/name[] arrays are
    // pre-sized to CHUNK_READ_N (1024) and zeroed, so this is safe for windows
    // with <=1024 sequences.
    // ------------------------------------------------------------------
    const int totSeq = readSeqId;  // backbone (0) + alignedReads
    ab->abs->n_seq = totSeq;

    // Populate abc (msa_base, msa_len) from the graph. This also fills
    // node_id_to_msa_rank, which the SNP detector below relies on, so it must
    // run for every window regardless of whether we dump the MSA to disk.
    abpoa_generate_rc_msa(ab, abpt);

    // Per-window FASTA/GFA dumps are a debugging aid and, at genome scale, would
    // write thousands of files into the run directory. Off by default; enable
    // with DINARA_MSA_DUMP=1 (or true/yes/on) to inspect individual windows.
    static const bool dumpMsaFiles = []() {
        const char* env = getenv("DINARA_MSA_DUMP");
        if(env == nullptr) return false;
        const string s(env);
        return s == "1" || s == "true" || s == "yes" || s == "on"
            || s == "TRUE" || s == "YES" || s == "ON";
    }();
    if(dumpMsaFiles) {
        {
            const string msaPath = "testAbpoaMultiSegmentMSA_window"
                + to_string(window.windowId) + ".fasta";
            FILE* msaFp = fopen(msaPath.c_str(), "w");
            if(msaFp) {
                // abpoa_output_rc_msa returns early if msa_len<=0; the generate
                // step above is required for a manually-built graph.
                abpoa_output_rc_msa(ab, abpt, msaFp);
                fclose(msaFp);
                out << "  MSA written to " << msaPath << endl;
            } else {
                out << "  WARNING: could not open " << msaPath << " for writing" << endl;
            }
        }
        {
            const string gfaPath = "testAbpoaMultiSegmentMSA_window"
                + to_string(window.windowId) + ".gfa";
            FILE* gfaFp = fopen(gfaPath.c_str(), "w");
            if(gfaFp) {
                abpt->out_gfa = 1;
                abpoa_generate_gfa(ab, abpt, gfaFp);
                abpt->out_gfa = 0;
                fclose(gfaFp);
                out << "  GFA written to " << gfaPath << endl;
            }
        }
    }

    // ------------------------------------------------------------------
    // Detect clean SNP bubbles in the finished graph. The abpoa_generate_rc_msa
    // call above has populated node_id_to_msa_rank, so columns have ranks.
    // ------------------------------------------------------------------
    {
        // Per-allele support cutoff, derived from dataset coverage using the
        // same rule hifiasm applies at the final het-site stage (see
        // AssemblerHifiasmEC.cpp): cc = max(cut_bd, (coverageHet / n_hap) *
        // cut_rate) with n_hap=2, cut_rate=0.7, cut_bd=6. coverageHet is the
        // k-mer histogram peak (total diploid coverage), so coverageHet/2 is
        // the expected per-haplotype depth and a het allele must carry >=70% of
        // it, floored at 6. When coverageHet is unavailable (histogram not
        // computed) fall back to the floor so we never under-gate.
        // Per-allele support cutoff. If the caller supplied a nonzero
        // hetMinSupport (Assembly.mode3.hetMinSupport) use it verbatim;
        // otherwise auto-derive from coverage as before.
        int minSupport = 6;
        if(hetMinSupport > 0) {
            minSupport = static_cast<int>(hetMinSupport);
        } else {
            constexpr uint64_t cut_bd = 6;
            constexpr uint64_t cut_rate_num = 7;
            constexpr uint64_t cut_rate_den = 10;
            constexpr uint64_t n_hap = 2;
            const uint64_t coverageHet = assemblerInfo.isOpen ?
                assemblerInfo->kmerDistributionInfo.coverageHet : invalid<uint64_t>;
            uint64_t base = 0;
            if(coverageHet != invalid<uint64_t> && coverageHet > 0)
                base = coverageHet / n_hap;
            uint64_t cc = (base * cut_rate_num) / cut_rate_den;
            if(cc < cut_bd) cc = cut_bd;
            minSupport = static_cast<int>(cc);
        }
        const double minVaf = hetMinVaf;
        // Repeat-context SNPs are kept by default (see
        // Assembly.mode3.hetDropHomopolymer / hetDropRepeat): the
        // flank-linearity test already requires a clean homozygous base on each
        // side, so such SNPs are real het sites; dropping them discarded far
        // more true SNPs than it kept. The two classes are gated separately.
        const bool dropHomopolymer = hetDropHomopolymer;
        const bool dropRepeat = hetDropRepeat;
        const vector<WindowSnp> snps = detectWindowSnps(
            ab, backboneQposToNode, backboneCodes, seqIdToOrientedRead,
            seqIdNodeToReadPos, minSupport, minVaf,
            dropHomopolymer, dropRepeat);
        out << "  SNPs detected: " << snps.size()
             << " (minSupport=" << minSupport << ", minVAF=" << minVaf
             << ", dropHomopolymer=" << (dropHomopolymer ? "true" : "false")
             << ", dropRepeat=" << (dropRepeat ? "true" : "false")
             << ")" << endl;
        for(const WindowSnp& snp : snps) {
            const WindowSnpAllele& ref = snp.alleles.front();
            out << "    col rank=" << snp.msaRank
                 << " backboneOff=" << snp.backboneOffset
                 << " " << ref.base << ">";
            // alt bases
            for(size_t ai = 1; ai < snp.alleles.size(); ai++)
                out << (ai > 1 ? "," : "") << snp.alleles[ai].base;
            out << " ref=" << ref.support << " alt=";
            for(size_t ai = 1; ai < snp.alleles.size(); ai++)
                out << (ai > 1 ? "," : "") << snp.alleles[ai].support;
            const int nAlt = int(snp.alleles.size()) - 1;
            const double topVaf = (snp.spanning > 0 && nAlt > 0)
                ? double(snp.alleles[1].support) / double(snp.spanning) : 0.0;
            out << " del=" << snp.delSupport
                 << " span=" << snp.spanning
                 << " VAF=" << topVaf
                 << (nAlt > 1 ? " [multi]" : "")
                 << (snp.delSupport > 0 ? " [+del]" : "")
                 << (snp.inHomopolymerOrRepeat ? " [repeat]" : "")
                 << " arms=" << snp.alleles.size()
                 << " refMembers=" << ref.members.size();
            for(size_t ai = 1; ai < snp.alleles.size(); ai++)
                out << " altMembers[" << ai << "]=" << snp.alleles[ai].members.size();
            out << endl;
        }

        // Stage het-anchor descriptors on the window: one k=2 anchor per allele
        // per SNP record. A biallelic site yields 2 arms (ref + alt); a
        // multiallelic site yields N arms (ref + each strong alt) in a SINGLE
        // bubble. Members carry rawPosition at the bubble's common predecessor;
        // predBase is the backbone base there, the allele base is this allele's
        // base. IDs are assigned in the post-window append pass. Only stage
        // arms that have at least one recoverable member.
        for(const WindowSnp& snp : snps) {
            // predBase: backbone base at (backboneOffset - 1); if the SNP is at
            // offset 0 there is no predecessor base, so skip (cannot form the
            // 2-base marker).
            if(snp.backboneOffset <= 0) continue;
            const uint8_t predBase = backboneCodes[snp.backboneOffset - 1];

            AnchorWindow::HetBubble bubble;
            bubble.backboneOffset = static_cast<uint32_t>(snp.backboneOffset);

            for(const WindowSnpAllele& allele : snp.alleles) {
                if(allele.members.empty()) continue;  // no recoverable member
                AnchorWindow::HetAnchor arm;
                arm.backboneOffset = bubble.backboneOffset;
                arm.predBase = predBase;
                arm.alleleBase = allele.code;
                arm.isRef = allele.isRef;
                for(const HetAlleleMember& m : allele.members)
                    arm.members.push_back({m.orientedReadId, m.rawPosition});
                bubble.alleles.push_back(std::move(arm));
            }
            // Need at least two arms to form a bubble.
            if(bubble.alleles.size() < 2) continue;

            // Leading hom anchor [predPrevBase, predBase] at predPrev. This
            // brackets the bubble upstream so the interval's backbone anchor
            // connects to a hom (shared by every entering read) rather than to a
            // minority allele arm (whose reads the k=50 backbone anchor does not
            // share). Both homs are required for a wired bubble; if the leading
            // hom is unavailable (flank not on a backbone column, or no
            // recoverable members) the plan pass drops the bubble.
            if(snp.predBackboneOffset >= 0 && !snp.leadHomMembers.empty()) {
                bubble.predBackboneOffset = static_cast<uint32_t>(snp.predBackboneOffset);
                AnchorWindow::HetAnchor leadHom;
                leadHom.backboneOffset = bubble.predBackboneOffset;
                leadHom.predBase = snp.predPrevBase;
                // alleleBase is predBase (the linear next base after predPrev).
                leadHom.alleleBase = predBase;
                leadHom.isRef = true;
                for(const HetAlleleMember& m : snp.leadHomMembers)
                    leadHom.members.push_back({m.orientedReadId, m.rawPosition});
                bubble.leadHom = std::move(leadHom);
            }

            // Trailing hom anchor [succBase, nextBase] at commonSucc. predBase
            // of the hom is succBase (the base shared by all reads at
            // commonSucc); alleleBase is the linear next base. Members are all
            // spanning reads. Brackets the bubble downstream.
            if(snp.succBackboneOffset >= 0 && !snp.homMembers.empty()) {
                bubble.succBackboneOffset = static_cast<uint32_t>(snp.succBackboneOffset);
                AnchorWindow::HetAnchor homAnchor;
                homAnchor.backboneOffset = bubble.succBackboneOffset;
                homAnchor.predBase = snp.succBase;
                // nextBase is the backbone base right after commonSucc.
                homAnchor.alleleBase =
                    (snp.succBackboneOffset + 1 < static_cast<int>(backboneCodes.size()))
                    ? backboneCodes[snp.succBackboneOffset + 1] : 0;
                homAnchor.isRef = true;
                for(const HetAlleleMember& m : snp.homMembers)
                    homAnchor.members.push_back({m.orientedReadId, m.rawPosition});
                bubble.hom = std::move(homAnchor);
            }

            window.hetBubbles.push_back(std::move(bubble));
        }
        out << "  staged het bubbles: " << window.hetBubbles.size() << endl;
    }

    abpoa_free(ab);
    abpoa_free_para(abpt);
    return true;
}


// Driver: build a per-window all-reads abPOA multi-segment MSA for each anchor
// window and detect het sites, one window per work item with dynamic load
// balancing (per-window cost varies with read count / backbone length, so a
// dynamic schedule keeps threads busy). Windows are processed independently:
// each worker writes only to its own AnchorWindow (window.hetBubbles), so there
// is no output data race. Het-anchor CREATION is NOT done here -- it grows the
// shared memory-mapped anchor store and must stay in the serial append pass that
// runs after this. By default ALL windows are processed; env
// DINARA_MSA_MAX_WINDOWS=N caps at the first N (0 = all).
void Assembler::testAbpoaMultiSegmentMSA(
    const shared_ptr<Shasta2Anchors>& shasta2Anchors,
    const shared_ptr<Shasta2Journeys>& shasta2Journeys,
    vector<AnchorWindow>& anchorWindows,
    uint64_t threadCount,
    double hetMinVaf,
    uint64_t hetMinSupport,
    bool hetDropHomopolymer,
    bool hetDropRepeat)
{
    if(anchorWindows.empty()) {
        cout << "testAbpoaMultiSegmentMSA: no windows." << endl;
        return;
    }

    // Default: process ALL windows (genome-wide het detection). Set
    // DINARA_MSA_MAX_WINDOWS=N to cap at the first N windows for debugging
    // (0 also means all).
    uint64_t maxWindows = 0;
    if(const char* env = getenv("DINARA_MSA_MAX_WINDOWS")) {
        maxWindows = strtoull(env, nullptr, 10);  // 0 = all windows
    }
    const uint64_t windowEnd = (maxWindows == 0) ?
        anchorWindows.size() : std::min<uint64_t>(maxWindows, anchorWindows.size());

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
        if(threadCount == 0) threadCount = 1;
    }

    cout << "testAbpoaMultiSegmentMSA: " << anchorWindows.size()
         << " windows available, processing " << windowEnd
         << " on " << threadCount << " threads." << endl;

    // Publish the shared inputs/outputs for the thread function.
    abpoaMultiSegmentMSAData.shasta2Anchors = &shasta2Anchors;
    abpoaMultiSegmentMSAData.shasta2Journeys = &shasta2Journeys;
    abpoaMultiSegmentMSAData.anchorWindows = &anchorWindows;
    abpoaMultiSegmentMSAData.windowEnd = windowEnd;
    abpoaMultiSegmentMSAData.processed = 0;
    abpoaMultiSegmentMSAData.produced = 0;
    abpoaMultiSegmentMSAData.hetMinVaf = hetMinVaf;
    abpoaMultiSegmentMSAData.hetMinSupport = hetMinSupport;
    abpoaMultiSegmentMSAData.hetDropHomopolymer = hetDropHomopolymer;
    abpoaMultiSegmentMSAData.hetDropRepeat = hetDropRepeat;

    // One work item per window (batch size 1: per-window cost is high and
    // uneven, so fine-grained dynamic balancing is what we want).
    setupLoadBalancing(windowEnd, 1);
    runThreads(&Assembler::testAbpoaMultiSegmentMSAThreadFunction, threadCount);

    cout << "testAbpoaMultiSegmentMSA: produced MSA for "
         << abpoaMultiSegmentMSAData.produced.load()
         << " of " << abpoaMultiSegmentMSAData.processed.load()
         << " processed windows." << endl;
}


void Assembler::testAbpoaMultiSegmentMSAThreadFunction(size_t)
{
    auto& data = abpoaMultiSegmentMSAData;
    const auto& shasta2Anchors = *data.shasta2Anchors;
    const auto& shasta2Journeys = *data.shasta2Journeys;
    auto& anchorWindows = *data.anchorWindows;

    uint64_t begin = 0, end = 0;
    while(getNextBatch(begin, end)) {
        for(uint64_t i = begin; i != end; i++) {
            AnchorWindow& window = anchorWindows[i];

            // Buffer this window's diagnostics so they flush as one block, not
            // interleaved with other threads' output.
            std::ostringstream buffer;
            const bool ok = runOneWindowAbpoaMultiSegmentMSA(
                shasta2Anchors, shasta2Journeys, window, buffer,
                data.hetMinVaf, data.hetMinSupport,
                data.hetDropHomopolymer, data.hetDropRepeat);

            data.processed.fetch_add(1, std::memory_order_relaxed);
            if(ok) data.produced.fetch_add(1, std::memory_order_relaxed);

            {
                std::lock_guard<std::mutex> lock(mutex);
                cout << buffer.str();
            }
        }
    }
}
