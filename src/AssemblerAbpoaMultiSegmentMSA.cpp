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
#include <chrono>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
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

// A detected clean biallelic SNP in a window POA graph.
struct WindowSnp {
    int msaRank = -1;                 // MSA column (topological rank) of the bubble
    int backboneOffset = -1;          // backbone base offset (into backboneCodes), -1 if backbone gaps here
    char refBase = 'N';               // backbone (reference) allele
    char altBase = 'N';               // alternate allele
    int refSupport = 0;               // reads carrying refBase
    int altSupport = 0;               // reads carrying altBase
    int delSupport = 0;               // reads carrying a coexisting deletion (skip edge)
    int spanning = 0;                 // reads spanning the column (VAF denominator)
    double vaf = 0.0;                  // altSupport / spanning (this alt allele)
    int nAltAlleles = 1;              // # of distinct alt bases at this column (>1 = multiallelic)
    bool inHomopolymerOrRepeat = false; // backbone context is a homopolymer/STR run
    vector<OrientedReadId> refReads;  // reads on the reference allele (phasing)
    vector<OrientedReadId> altReads;  // reads on the alternate allele (phasing)
    // Het-anchor members (read + rawPosition at commonPred) for the two
    // alleles. Parallel to refReads/altReads in the reads they cover, but only
    // include reads whose predecessor position was recoverable from the POA.
    vector<HetAlleleMember> refMembers;
    vector<HetAlleleMember> altMembers;

    // Hom-separator anchor members: reads at commonSucc (rawPosition = succ base
    // position), forming the k=2 hom anchor [succBase, nextBase]. The
    // flank-linearity test guarantees commonSucc->succNext is linear, so this
    // 2-mer is shared by ALL reads spanning the site (both alleles reconverge
    // here). Used to separate consecutive SNPs chained in the same backbone
    // interval. succBackboneOffset is the backbone offset of commonSucc.
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
// Filters: distinct non-gap bases, backbone present in the bubble, alt support
// >= minSupport, and VAF >= minVaf (VAF = altSupport / (refSupport+altSupport)).
vector<WindowSnp> detectWindowSnps(
    abpoa_t* ab,
    const vector<int>& backboneQposToNode,
    const vector<uint8_t>& backboneCodes,
    const vector<OrientedReadId>& seqIdToOrientedRead,
    const vector<unordered_map<int, uint32_t>>& seqIdNodeToReadPos,
    int minSupport,
    double minVaf,
    bool dropRepeatContext)
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

        // Collect EVERY distinct non-reference base allele. A multiallelic site
        // (ref + 2 alt bases) is emitted as separate biallelic records sharing
        // the same column and spanning depth, matching how the rest of the
        // pipeline (KmVarKey, one altBase per record) represents variants. The
        // total base-allele depth sums all alleles so per-alt VAF is correct.
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
        const int nAltAlleles = static_cast<int>(alts.size());

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
        const bool inRepeat = kmIsHomopolymer(
            backboneCodes.data(),
            static_cast<uint32_t>(backboneCodes.size()),
            vkey, /* xid */ 0);
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

        if(dropRepeatContext && inRepeat) continue;

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

        // Emit one biallelic record per distinct alt allele. Each is gated on
        // its OWN support/VAF, so a strong alt is kept even if a second weak alt
        // at the same column is filtered out.
        for(const AltAllele& alt : alts) {
            if(alt.support < minSupport) continue;
            const double vaf = (spanning > 0) ? double(alt.support) / double(spanning) : 0.0;
            if(vaf < minVaf) continue;

            WindowSnp snp;
            snp.msaRank = (abg->node_id_to_msa_rank != nullptr) ? abg->node_id_to_msa_rank[refNode] : -1;
            snp.backboneOffset = backboneOff;
            snp.inHomopolymerOrRepeat = inRepeat;
            snp.refBase = codeToChar(refCode);
            snp.altBase = codeToChar(alt.base);
            snp.refSupport = refSupport;
            snp.altSupport = alt.support;
            snp.delSupport = delSupport;
            snp.spanning = spanning;
            snp.vaf = vaf;
            snp.nAltAlleles = nAltAlleles;
            mapReads(refSeqIds, snp.refReads);
            mapReads(alt.seqIds, snp.altReads);
            mapMembers(refSeqIds, snp.refMembers);
            mapMembers(alt.seqIds, snp.altMembers);

            // Hom-separator anchor at commonSucc: [succBase, nextBase]. All
            // reads spanning the site (ref + this alt) reconverge here, so
            // recover each one's position at commonSucc. rawPosition = succ base
            // position; the flank-linearity test guarantees the next base is
            // linear, so [succBase, nextBase] is shared by all these reads.
            snp.succBackboneOffset = nodeToBackboneOffset[commonSucc];
            snp.succBase = abg->node[commonSucc].base;
            {
                vector<int> homSeqIds = refSeqIds;
                homSeqIds.insert(homSeqIds.end(), alt.seqIds.begin(), alt.seqIds.end());
                mapMembersAt(homSeqIds, commonSucc, snp.homMembers);
            }
            snps.push_back(std::move(snp));
        }
    }

    // Report in backbone order; alt records at the same column stay grouped.
    sort(snps.begin(), snps.end(), [](const WindowSnp& a, const WindowSnp& b) {
        if(a.backboneOffset != b.backboneOffset) return a.backboneOffset < b.backboneOffset;
        return a.altBase < b.altBase;
    });
    cout << "  homopolymer/STR SNPs "
         << (dropRepeatContext ? "dropped: " : "flagged: ")
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
    AnchorWindow& window)
{
    const Reads& readsRef = getReads();
    const auto& markersRef = *markers;
    const uint64_t k = assemblerInfo->k;

    const OrientedReadId backboneOid = window.backboneOrientedReadId;
    const auto backboneJourney = (*shasta2Journeys)[backboneOid];

    cout << "testAbpoaMultiSegmentMSA: window " << window.windowId
         << " backbone " << backboneOid
         << " anchors [" << window.backboneBegin << "," << window.backboneEnd << ")"
         << " reads " << window.readIntervals.size() << endl;

    const uint32_t nBackboneAnchors = window.backboneEnd - window.backboneBegin;
    if(nBackboneAnchors < 2) {
        cout << "  window " << window.windowId << " has < 2 anchors, skipping." << endl;
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
        cout << "  backbone span empty, skipping window." << endl;
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

    cout << "  backbone " << backboneLen << " bases across "
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

    cout << "  reads with >=2 shared anchors: " << readsBySpan.size()
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
                            uint32_t boundaryForLog) {
            if(gapLen <= 0) return;
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
                cout << "  SLOW: read " << oid
                     << " gap near anchor " << boundaryForLog
                     << " seq " << gapLen << " bases took " << elapsed << "s" << endl;
            }
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
            alignGap(begNode, pinNode, begQpos, pinQpos - begQpos, boundary);

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
            alignGap(begNode, foldEndNode, begQpos, qlen - begQpos, lastBoundary);
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

    cout << "  aligned " << alignedReads << " reads ("
         << alignedSegments << " segments), skipped " << skippedReads << endl;

    // Report timing by quartile (reads are already longest-first).
    if(perReadTime.size() >= 4) {
        size_t q = perReadTime.size() / 4;
        for(int qi = 0; qi < 4; qi++) {
            size_t start = qi * q;
            size_t end = (qi == 3) ? perReadTime.size() : (qi + 1) * q;
            double sum = 0;
            for(size_t i = start; i < end; i++) sum += perReadTime[i];
            cout << "  Q" << (qi + 1) << " (reads " << start << "-" << (end - 1)
                 << "): " << sum << "s total, "
                 << (sum / (end - start) * 1000) << "ms/read" << endl;
        }
    }
    cout << "  total align time: " << totalAlignTime << "s"
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

    {
        const string msaPath = "testAbpoaMultiSegmentMSA_window"
            + to_string(window.windowId) + ".fasta";
        FILE* msaFp = fopen(msaPath.c_str(), "w");
        if(msaFp) {
            // Populate abc (msa_base, msa_len) from the graph, then write it.
            // abpoa_output_rc_msa returns early if msa_len<=0, so the generate
            // step is required for a manually-built graph.
            abpoa_generate_rc_msa(ab, abpt);
            abpoa_output_rc_msa(ab, abpt, msaFp);
            fclose(msaFp);
            cout << "  MSA written to " << msaPath << endl;
        } else {
            cout << "  WARNING: could not open " << msaPath << " for writing" << endl;
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
            cout << "  GFA written to " << gfaPath << endl;
        }
    }

    // ------------------------------------------------------------------
    // Detect clean SNP bubbles in the finished graph. abpoa_generate_rc_msa
    // (called above) has populated node_id_to_msa_rank, so columns have ranks.
    // ------------------------------------------------------------------
    {
        const int minSupport = 3;
        const double minVaf = 0.12;
        const bool dropRepeatContext = true;  // drop SNPs in homopolymer/STR runs
        const vector<WindowSnp> snps = detectWindowSnps(
            ab, backboneQposToNode, backboneCodes, seqIdToOrientedRead,
            seqIdNodeToReadPos, minSupport, minVaf, dropRepeatContext);
        cout << "  SNPs detected: " << snps.size()
             << " (minSupport=" << minSupport << ", minVAF=" << minVaf
             << ", dropRepeatContext=" << (dropRepeatContext ? "true" : "false")
             << ")" << endl;
        for(const WindowSnp& snp : snps) {
            cout << "    col rank=" << snp.msaRank
                 << " backboneOff=" << snp.backboneOffset
                 << " " << snp.refBase << ">" << snp.altBase
                 << " ref=" << snp.refSupport
                 << " alt=" << snp.altSupport
                 << " del=" << snp.delSupport
                 << " span=" << snp.spanning
                 << " VAF=" << snp.vaf
                 << (snp.nAltAlleles > 1 ? " [multi]" : "")
                 << (snp.delSupport > 0 ? " [+del]" : "")
                 << (snp.inHomopolymerOrRepeat ? " [repeat]" : "")
                 << " refMembers=" << snp.refMembers.size()
                 << " altMembers=" << snp.altMembers.size() << endl;
        }

        // Stage het-anchor descriptors on the window: one k=2 anchor per allele
        // (ref + this alt) per SNP record. Members carry rawPosition at the
        // bubble's common predecessor; predBase is the backbone base there, the
        // allele base is ref/alt. IDs are assigned in the post-window append
        // pass. Only stage alleles that have at least one recoverable member.
        for(const WindowSnp& snp : snps) {
            if(snp.refMembers.empty() && snp.altMembers.empty()) continue;
            // predBase: backbone base at (backboneOffset - 1); if the SNP is at
            // offset 0 there is no predecessor base, so skip (cannot form the
            // 2-base marker).
            if(snp.backboneOffset <= 0) continue;
            const uint8_t predBase = backboneCodes[snp.backboneOffset - 1];
            const uint8_t refCode = static_cast<uint8_t>(
                snp.refBase == 'A' ? 0 : snp.refBase == 'C' ? 1 :
                snp.refBase == 'G' ? 2 : snp.refBase == 'T' ? 3 : 0);
            const uint8_t altCode = static_cast<uint8_t>(
                snp.altBase == 'A' ? 0 : snp.altBase == 'C' ? 1 :
                snp.altBase == 'G' ? 2 : snp.altBase == 'T' ? 3 : 0);

            AnchorWindow::HetBubble bubble;
            bubble.backboneOffset = static_cast<uint32_t>(snp.backboneOffset);

            AnchorWindow::HetAnchor refAnchor;
            refAnchor.backboneOffset = bubble.backboneOffset;
            refAnchor.predBase = predBase;
            refAnchor.alleleBase = refCode;
            refAnchor.isRef = true;
            for(const HetAlleleMember& m : snp.refMembers)
                refAnchor.members.push_back({m.orientedReadId, m.rawPosition});

            AnchorWindow::HetAnchor altAnchor;
            altAnchor.backboneOffset = bubble.backboneOffset;
            altAnchor.predBase = predBase;
            altAnchor.alleleBase = altCode;
            altAnchor.isRef = false;
            for(const HetAlleleMember& m : snp.altMembers)
                altAnchor.members.push_back({m.orientedReadId, m.rawPosition});

            bubble.alleles.push_back(std::move(refAnchor));
            bubble.alleles.push_back(std::move(altAnchor));

            // Hom separator anchor [succBase, nextBase] at commonSucc. predBase
            // of the hom is succBase (the base shared by all reads at
            // commonSucc); alleleBase is the linear next base. Members are all
            // spanning reads. Used to chain/separate consecutive SNPs.
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
        cout << "  staged het bubbles: " << window.hetBubbles.size() << endl;
    }

    abpoa_free(ab);
    abpoa_free_para(abpt);
    return true;
}


// Driver: build a per-window all-reads abPOA multi-segment MSA for each anchor
// window, mirroring testMultiSegmentMSA so the two engines can be compared.
// Number of windows capped by env DINARA_MSA_MAX_WINDOWS (default 1, 0=all).
void Assembler::testAbpoaMultiSegmentMSA(
    const shared_ptr<Shasta2Anchors>& shasta2Anchors,
    const shared_ptr<Shasta2Journeys>& shasta2Journeys,
    vector<AnchorWindow>& anchorWindows)
{
    if(anchorWindows.empty()) {
        cout << "testAbpoaMultiSegmentMSA: no windows." << endl;
        return;
    }

    uint64_t maxWindows = 1;
    if(const char* env = getenv("DINARA_MSA_MAX_WINDOWS")) {
        maxWindows = strtoull(env, nullptr, 10);  // 0 = all windows
    }

    cout << "testAbpoaMultiSegmentMSA: " << anchorWindows.size()
         << " windows available";
    if(maxWindows == 0) cout << ", processing all." << endl;
    else cout << ", processing up to " << maxWindows << "." << endl;

    uint64_t processed = 0, produced = 0;
    for(AnchorWindow& window : anchorWindows) {
        if(maxWindows != 0 && processed >= maxWindows) break;
        processed++;
        if(runOneWindowAbpoaMultiSegmentMSA(shasta2Anchors, shasta2Journeys, window)) {
            produced++;
        }
    }

    cout << "testAbpoaMultiSegmentMSA: produced MSA for " << produced
         << " of " << processed << " processed windows." << endl;
}
