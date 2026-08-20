// Transcribe abPOA local topology into the Shasta2AnchorGraph.
//
// For each anchor-graph edge A -> B whose two-sided read coverage is at least
// minCommonForHet, we run an abPOA MSA over the reads' inter-anchor sequences.
// The MSA column matrix (msa_base[row][col], row order == anchorPair
// orientedReadIds because sort_input_seq=0) is a materialization of the POA
// graph in rank order. We read that topology and transcribe it into the anchor
// graph's own vocabulary:
//
//   * A column where every covered read carries the SAME non-gap base is a
//     strict convergence point (a "hom" position). The edge's own endpoints A
//     and B are the outer convergence points by construction (all reads share
//     both). Interior convergence columns separate consecutive bubbles.
//
//   * A maximal run of columns between two convergence points where reads
//     disagree is a bubble. Reads are grouped by their exact subsequence over
//     the run: each distinct group is an allele arm. An arm that clears the
//     per-allele support floor becomes a het anchor (appendHetAnchorPair); the
//     original mixed edge is replaced by flank -> arm -> flank edges. A read
//     whose subsequence over the run is empty (it skips the run: a deletion
//     allele) takes the direct flank -> flank edge, so the deletion IS the
//     short path -- no arm anchor for it.
//
// This unifies SNPs (a one-column bubble) and big indels (a multi-column
// bubble, one side long / one side empty) with no size special-case: indel
// length falls out of the run length, "both alleles pass" is just counting
// reads per group.
//
// Two phases:
//   1. Planning (parallel, read-only): per edge, build an EdgePlan describing
//      the bubbles and their arms. Touches only the graph for reads; allocates
//      no anchors, mutates nothing.
//   2. Apply (serial): appendHetAnchorPair for each arm (grows the anchor
//      store, so must be serial), add the matching graph vertices in lockstep
//      (vertex_descriptor == AnchorId), wire the fwd + RC flank->arm->flank
//      edges, and disable the original edge (+ its RC mirror).
//
// Iteration 1 scope: to keep the point-anchor (k=2) representation provably
// round-trippable, we only transcribe edges with exactly ONE interior bubble
// that is flanked on both sides by the edge endpoints (a clean single-site
// edge: one SNP, or one >=15bp indel, or one small indel). Edges with multiple
// interior bubbles, or a bubble touching an end, are counted and left
// unchanged for a later iteration. This matches the "sparse now, revisit if too
// coarse" decision.

#include "Shasta2AnchorGraph.hpp"

#include "Base.hpp"
#include "Reads.hpp"
#include "Shasta2AnchorPair.hpp"
#include "Shasta2Anchors.hpp"
#include "WindowIntervalPoa.hpp"   // IpoaAbHandle

#include <boost/graph/iteration_macros.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <thread>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// Longest-row length above which abpoa's DP can blow up memory on an
// anchor-less span; such edges are skipped. Mirrors DINARA_IPOA_MAX_LEN on the
// window het path.
int hetOnGraphMaxLen()
{
    static const int cap = [] {
        const char* e = std::getenv("DINARA_IPOA_MAX_LEN");
        if(e == nullptr) return 16384;
        return std::atoi(e);
    }();
    return cap;
}

// Per-allele support floor: an arm must carry at least this many reads to
// become an anchor. Coverage-relative: max(2, covered/5). appendHetAnchorPair
// itself requires >= 2 members.
uint64_t alleleFloor(uint64_t coveredRows)
{
    const uint64_t frac = coveredRows / 5;
    return frac < 2 ? 2 : frac;
}

// One planned arm of one bubble: the reads carrying this allele and, for each,
// the raw read position that pins its het anchor (the first base of the
// allele's subsequence, in the read's own coordinate frame). Empty members ==
// the deletion allele: no anchor, reads take the direct flank->flank edge.
struct ArmPlan {
    vector<pair<OrientedReadId, uint32_t>> members;   // (read, rawPosition)
    bool isDeletion = false;                          // true => no arm anchor
};

// One planned bubble on an edge: the two flanking anchor ids (here always the
// edge endpoints in iteration 1) and the allele arms.
struct BubblePlan {
    Shasta2AnchorId flankA = invalid<Shasta2AnchorId>;
    Shasta2AnchorId flankB = invalid<Shasta2AnchorId>;
    vector<ArmPlan> arms;
};

// The transcription plan for one edge. Empty bubbles => nothing to do.
struct EdgePlan {
    Shasta2AnchorGraph::edge_descriptor edge;
    vector<BubblePlan> bubbles;
};

// Result of analyzing one edge's MSA matrix.
enum class AnalyzeResult {
    NoBubble,        // fully homozygous span (all columns convergent)
    Planned,         // exactly one clean interior bubble -> plan produced
    DeferMultiSite,  // >1 interior bubble (iteration 2)
    DeferEndBubble,  // a bubble touches column 0 or the last column
    DeferComplex     // other unsupported shape
};

// Analyze the row-major MSA matrix (rows == orientedReadIds order) for one edge
// and, on the clean single-bubble case, fill `out`. `msa[r]` is row r's aligned
// bases (value 0..3, or 4 for gap), all of length ncols. `posA[r]` is read r's
// absolute base position of anchorA (so aligned column maps back to a raw read
// position by counting that row's non-gap bases up to the column). `total` is
// the row count.
AnalyzeResult analyzeEdgeMsa(
    const vector<vector<uint8_t>>& msa,
    const vector<uint32_t>& posA,
    const vector<OrientedReadId>& rowRead,   // rowRead[r] == read of MSA row r
    uint64_t total,
    BubblePlan& out)
{
    if(msa.empty()) return AnalyzeResult::NoBubble;
    const uint64_t ncols = msa.front().size();
    if(ncols == 0) return AnalyzeResult::NoBubble;

    // A column is convergent iff every row carries the SAME non-gap base
    // (strict: no gaps, single base value across all rows).
    auto isConvergent = [&](uint64_t c) -> bool {
        uint8_t b = msa[0][c];
        if(b >= 4) return false;
        for(uint64_t r=1; r<total; r++) {
            if(msa[r][c] != b) return false;
        }
        return true;
    };

    // Find maximal runs of non-convergent columns (candidate bubbles).
    vector<pair<uint64_t, uint64_t>> runs;   // [begin, end) column ranges
    uint64_t c = 0;
    while(c < ncols) {
        if(isConvergent(c)) { c++; continue; }
        uint64_t begin = c;
        while(c < ncols && !isConvergent(c)) c++;
        runs.push_back({begin, c});
    }

    if(runs.empty()) return AnalyzeResult::NoBubble;
    if(runs.size() > 1) return AnalyzeResult::DeferMultiSite;

    const auto [rbegin, rend] = runs.front();
    // Iteration 1 requires the bubble to be strictly interior: at least one
    // convergent column on each side (so both flanks are the edge endpoints and
    // every read is defined before/after the run).
    if(rbegin == 0 || rend == ncols) return AnalyzeResult::DeferEndBubble;

    // Precompute, per row, the raw read position at column rbegin: posA[r] plus
    // the number of that row's non-gap bases strictly before rbegin.
    vector<uint32_t> rawAtRunStart(total);
    for(uint64_t r=0; r<total; r++) {
        uint32_t nonGap = 0;
        for(uint64_t col=0; col<rbegin; col++) {
            if(msa[r][col] < 4) nonGap++;
        }
        rawAtRunStart[r] = posA[r] + nonGap;
    }

    // Group rows by their exact subsequence (non-gap bases) over [rbegin, rend).
    // Key: the concatenated non-gap base values. Empty key == deletion allele.
    std::map<vector<uint8_t>, vector<uint64_t>> groups;   // allele -> row indices
    for(uint64_t r=0; r<total; r++) {
        vector<uint8_t> allele;
        for(uint64_t col=rbegin; col<rend; col++) {
            const uint8_t b = msa[r][col];
            if(b < 4) allele.push_back(b);
        }
        groups[allele].push_back(r);
    }

    // Build arms for groups that clear the per-allele floor.
    const uint64_t floor = alleleFloor(total);
    out = BubblePlan();
    uint64_t passing = 0;
    for(const auto& [allele, rows] : groups) {
        if(rows.size() < floor) continue;
        passing++;
        ArmPlan arm;
        if(allele.empty()) {
            arm.isDeletion = true;   // no anchor; reads take flank->flank
        } else {
            arm.members.reserve(rows.size());
            for(const uint64_t r: rows) {
                // Pin the arm anchor at the first base of the allele on this
                // read: rawAtRunStart[r] is the read's raw position of its first
                // non-gap base at/after the run start, which for a member of a
                // non-empty allele is exactly that allele's first base.
                arm.members.push_back({rowRead[r], rawAtRunStart[r]});
            }
        }
        out.arms.push_back(std::move(arm));
    }

    // Need at least two passing alleles for a real bubble; otherwise this is
    // effectively homozygous at the supported level (noise minority filtered).
    if(passing < 2) return AnalyzeResult::NoBubble;

    return AnalyzeResult::Planned;
}

// Build the abPOA MSA row matrix for one edge. Returns false (edge skipped) if
// fewer than 2 non-empty rows, or the longest row exceeds maxLenCap. On success
// fills `msa` (row-major, each row length == ncols, values 0..3 or 4=gap),
// `posA` (per surviving row, read's absolute anchorA base position), and
// `rowRead` (per surviving row, the OrientedReadId). Rows with an empty
// inter-anchor span are dropped, so surviving rows may be fewer than coverage;
// their read/posA are kept in sync.
bool buildEdgeMsa(
    IpoaAbHandle& ah,
    const Shasta2AnchorPair& anchorPair,
    const vector<pair<Shasta2AnchorPair::Positions, Shasta2AnchorPair::Positions>>& positions,
    const vector<vector<Base>>& sequences,
    int maxLenCap,
    vector<vector<uint8_t>>& msa,
    vector<uint32_t>& posA,
    vector<OrientedReadId>& rowRead)
{
    msa.clear();
    posA.clear();
    rowRead.clear();

    // Assemble abPOA code rows for non-empty sequences, tracking read + posA.
    vector<vector<uint8_t>> codeStore;
    vector<int> seqLens;
    vector<uint8_t*> seqPtrs;
    int maxLen = 0;
    for(uint64_t i=0; i<sequences.size(); i++) {
        const vector<Base>& seq = sequences[i];
        if(seq.empty()) continue;
        codeStore.emplace_back();
        vector<uint8_t>& codes = codeStore.back();
        codes.resize(seq.size());
        for(uint64_t j=0; j<seq.size(); j++) {
            codes[j] = seq[j].value;
        }
        if(int(codes.size()) > maxLen) maxLen = int(codes.size());
        posA.push_back(positions[i].first.basePosition);
        rowRead.push_back(anchorPair.orientedReadIds[i]);
    }

    if(codeStore.size() < 2) return false;
    if(maxLenCap > 0 && maxLen > maxLenCap) return false;

    const int nSeq = int(codeStore.size());
    seqLens.resize(nSeq);
    seqPtrs.resize(nSeq);
    for(int r=0; r<nSeq; r++) {
        seqLens[r] = int(codeStore[r].size());
        seqPtrs[r] = codeStore[r].data();
    }

    abpoa_reset(ah.ab, ah.abpt, maxLen > 0 ? maxLen : 1);
    abpoa_msa(ah.ab, ah.abpt, nSeq, nullptr, seqLens.data(),
              seqPtrs.data(), nullptr, nullptr);

    const abpoa_cons_t* abc = ah.ab->abc;
    if(abc == nullptr || abc->msa_len <= 0 || abc->n_seq != nSeq) return false;

    const uint64_t ncols = uint64_t(abc->msa_len);
    msa.resize(nSeq);
    for(int r=0; r<nSeq; r++) {
        msa[r].resize(ncols);
        for(uint64_t col=0; col<ncols; col++) {
            const uint8_t code = abc->msa_base[r][col];
            msa[r][col] = (code < 4) ? code : uint8_t(4);
        }
    }
    return true;
}

} // namespace

// Serial apply phase. Consumes the per-thread EdgePlans, creates het anchors,
// grows the graph vertex set in lockstep (vertex_descriptor == AnchorId), wires
// the fwd + RC flank->arm->flank edges, and disables the original edges.
// Defined at file scope (not in the anonymous namespace) so it can use the
// graph's private helpers via the public API only.
namespace dinara {
void applyHetBubblePlans(
    Shasta2AnchorGraph& graph,
    Shasta2Anchors& anchors,
    const vector<vector<EdgePlan>>& threadPlans,
    HetOnGraphResult& result);
}

void dinara::applyHetBubblePlans(
    Shasta2AnchorGraph& graph,
    Shasta2Anchors& anchors,
    const vector<vector<EdgePlan>>& threadPlans,
    HetOnGraphResult& result)
{
    // Helper: keep graph vertices in lockstep with the anchor store. After an
    // appendHetAnchorPair the store has 2 new anchors (canonical even, RC odd);
    // add exactly enough vertices so vertex_descriptor == AnchorId still holds.
    auto syncVertices = [&]() {
        const uint64_t anchorCount = anchors.size();
        while(num_vertices(graph) < anchorCount) {
            boost::add_vertex(graph);
        }
    };

    // Helper: add a directed edge idA -> idB carrying the intersection of the
    // two anchors' member reads, after dropping non-forward offsets. Mirrors the
    // window path's addHetEdge. Returns true if an edge was added.
    auto addHetEdge = [&](Shasta2AnchorId idA, Shasta2AnchorId idB) -> bool {
        const Shasta2Anchor anchorA = anchors[idA];
        const Shasta2Anchor anchorB = anchors[idB];
        Shasta2AnchorPair edgePair;
        edgePair.anchorIdA = idA;
        edgePair.anchorIdB = idB;
        auto itA = anchorA.begin();
        auto itB = anchorB.begin();
        while(itA != anchorA.end() && itB != anchorB.end()) {
            if(itA->orientedReadId < itB->orientedReadId) { ++itA; continue; }
            if(itB->orientedReadId < itA->orientedReadId) { ++itB; continue; }
            edgePair.orientedReadIds.push_back(itA->orientedReadId);
            ++itA; ++itB;
        }
        edgePair.removeNegativeOffsets(anchors);
        if(edgePair.orientedReadIds.empty()) return false;

        // Offset: het anchors are synthetic k=2 markers with no real marker
        // ordinals, so getAverageOffset (which indexes each read's marker array
        // by ordinal) would read out of bounds. Compute the average base-offset
        // directly from anchor member positions instead, which are valid for
        // both primary and het anchors. All kept reads are strictly forward
        // (positionB > positionA) after removeNegativeOffsets.
        vector<pair<uint32_t, uint32_t>> anchorPositions;
        edgePair.getAnchorPositions(anchors, anchorPositions);
        uint64_t offsetSum = 0;
        for(const auto& p: anchorPositions) {
            offsetSum += uint64_t(p.second - p.first);
        }
        const uint32_t nominalOffset = anchorPositions.empty() ? 0 :
            uint32_t(offsetSum / anchorPositions.size());

        Shasta2AnchorGraph::edge_descriptor e;
        bool added = false;
        boost::tie(e, added) = boost::add_edge(
            idA, idB,
            Shasta2AnchorGraphEdge(edgePair, nominalOffset, graph.nextEdgeId++),
            graph);
        graph[e].useForAssembly = true;
        return true;
    };

    // Build the RC member list for an arm: flip strand and mirror the raw
    // position, matching appendHetAnchorPair's own RC formula so the two agree.
    // We do NOT call appendHetAnchorPair for the RC arm separately -- that
    // function already appends BOTH strands and returns the canonical id, with
    // the RC at canonicalId+1. So here we only need the canonical member list.

    for(const vector<EdgePlan>& plans: threadPlans) {
        for(const EdgePlan& plan: plans) {
            bool anyArmCreated = false;

            for(const BubblePlan& bubble: plan.bubbles) {
                result.bubblesTranscribed++;

                for(const ArmPlan& arm: bubble.arms) {
                    if(arm.isDeletion) {
                        // Deletion allele: reads take the direct flank->flank
                        // edge. That edge already exists (it is the original
                        // edge) but will be disabled below; re-add a dedicated
                        // flank->flank edge carrying just... the intersection of
                        // the two flanks, which is all reads. Simplest correct
                        // behavior: add flankA->flankB via addHetEdge.
                        if(addHetEdge(bubble.flankA, bubble.flankB)) {
                            result.deletionEdgesAdded++;
                            anyArmCreated = true;
                        }
                        continue;
                    }
                    if(arm.members.size() < 2) continue;   // appendHetAnchorPair needs >=2

                    const Shasta2AnchorId armId =
                        anchors.appendHetAnchorPair(arm.members);
                    syncVertices();
                    result.hetAnchorsCreated++;

                    // Wire fwd: flankA -> arm -> flankB.
                    if(addHetEdge(bubble.flankA, armId))    result.armEdgesAdded++;
                    if(addHetEdge(armId, bubble.flankB))    result.armEdgesAdded++;

                    // Wire RC: flankB^1 -> arm^1 -> flankA^1. appendHetAnchorPair
                    // placed the RC arm at armId+1; the flanks' RC ids are
                    // flank^1 (canonical/RC pairing 2i/2i+1).
                    const Shasta2AnchorId armRc = armId + 1;
                    const Shasta2AnchorId flankArc = bubble.flankA ^ 1ULL;
                    const Shasta2AnchorId flankBrc = bubble.flankB ^ 1ULL;
                    if(uint64_t(flankBrc) < anchors.size() &&
                       uint64_t(flankArc) < anchors.size()) {
                        if(addHetEdge(flankBrc, armRc)) result.armEdgesAdded++;
                        if(addHetEdge(armRc, flankArc)) result.armEdgesAdded++;
                    }

                    anyArmCreated = true;
                }
            }

            if(anyArmCreated) {
                graph.disableEdge(plan.edge);
                result.originalEdgesDisabled++;
            }
        }
    }
}



HetOnGraphResult dinara::transcribeHetBubbles(
    Shasta2AnchorGraph& graph,
    Shasta2Anchors& anchors,
    uint64_t minCommonForHet,
    uint64_t threadCount)
{
    HetOnGraphResult result;
    const int maxLenCap = hetOnGraphMaxLen();
    const auto t0 = std::chrono::steady_clock::now();

    // Materialize edge descriptors so threads can index them.
    vector<Shasta2AnchorGraph::edge_descriptor> edgeList;
    edgeList.reserve(num_edges(graph));
    BGL_FORALL_EDGES(e, graph, Shasta2AnchorGraphBaseClass) {
        edgeList.push_back(e);
    }
    const uint64_t edgeCount = edgeList.size();
    result.edgesTotal = edgeCount;

    if(threadCount == 0) threadCount = std::thread::hardware_concurrency();
    if(threadCount == 0) threadCount = 1;
    if(threadCount > edgeCount && edgeCount > 0) threadCount = edgeCount;

    // ------------------------------------------------------------------
    // Phase 1: planning (parallel, read-only). Each thread analyzes edges off a
    // shared atomic counter and appends completed EdgePlans to a private list.
    // ------------------------------------------------------------------
    std::atomic<uint64_t> nextEdge{0};
    vector<vector<EdgePlan>> threadPlans(threadCount);
    vector<HetOnGraphResult> threadStats(threadCount);

    auto planner = [&](uint64_t threadId) {
        HetOnGraphResult& st = threadStats[threadId];
        vector<EdgePlan>& plans = threadPlans[threadId];

        IpoaAbHandle ah;
        vector<pair<Shasta2AnchorPair::Positions, Shasta2AnchorPair::Positions>> positions;
        vector<vector<Base>> sequences;
        vector<vector<uint8_t>> msa;
        vector<uint32_t> posA;
        vector<OrientedReadId> rowRead;

        for(;;) {
            const uint64_t ei = nextEdge.fetch_add(1);
            if(ei >= edgeCount) break;

            const Shasta2AnchorGraph::edge_descriptor e = edgeList[ei];
            const Shasta2AnchorGraphEdge& edge = graph[e];
            const Shasta2AnchorPair& anchorPair = edge.anchorPair;

            const uint64_t coverage = anchorPair.size();
            if(coverage < minCommonForHet) { st.edgesSkippedCoverage++; continue; }
            st.edgesConsidered++;

            anchorPair.get(anchors, positions, sequences);

            // Lossless identical early-out: byte-identical non-empty rows can
            // yield no bubble.
            {
                bool anyNonEmpty = false;
                const vector<Base>* first = nullptr;
                bool allIdentical = true;
                for(const vector<Base>& seq: sequences) {
                    if(seq.empty()) continue;
                    if(!anyNonEmpty) { first = &seq; anyNonEmpty = true; continue; }
                    if(seq != *first) { allIdentical = false; break; }
                }
                // Only skip when there are also no empty rows (an empty vs
                // non-empty split is itself a deletion bubble).
                bool anyEmpty = false;
                for(const vector<Base>& seq: sequences) {
                    if(seq.empty()) { anyEmpty = true; break; }
                }
                if(anyNonEmpty && allIdentical && !anyEmpty) {
                    st.edgesSkippedIdentical++;
                    continue;
                }
            }

            if(!buildEdgeMsa(ah, anchorPair, positions, sequences,
                             maxLenCap, msa, posA, rowRead)) {
                st.edgesSkippedLen++;
                continue;
            }
            st.edgesMsad++;

            BubblePlan bubble;
            const AnalyzeResult ar =
                analyzeEdgeMsa(msa, posA, rowRead, rowRead.size(), bubble);
            switch(ar) {
            case AnalyzeResult::NoBubble:
                break;
            case AnalyzeResult::DeferMultiSite:
                st.edgesDeferredMultiSite++;
                break;
            case AnalyzeResult::DeferEndBubble:
                st.edgesDeferredEndBubble++;
                break;
            case AnalyzeResult::DeferComplex:
                st.edgesDeferredComplex++;
                break;
            case AnalyzeResult::Planned: {
                bubble.flankA = anchorPair.anchorIdA;
                bubble.flankB = anchorPair.anchorIdB;
                EdgePlan plan;
                plan.edge = e;
                plan.bubbles.push_back(std::move(bubble));
                plans.push_back(std::move(plan));
                st.edgesPlanned++;
                break;
            }
            }
        }
    };

    {
        vector<std::thread> threads;
        threads.reserve(threadCount);
        for(uint64_t t=0; t<threadCount; t++) threads.emplace_back(planner, t);
        for(auto& th: threads) th.join();
    }

    // Merge planning stats.
    for(const HetOnGraphResult& st: threadStats) {
        result.edgesConsidered        += st.edgesConsidered;
        result.edgesMsad              += st.edgesMsad;
        result.edgesSkippedCoverage   += st.edgesSkippedCoverage;
        result.edgesSkippedLen        += st.edgesSkippedLen;
        result.edgesSkippedIdentical  += st.edgesSkippedIdentical;
        result.edgesPlanned           += st.edgesPlanned;
        result.edgesDeferredMultiSite += st.edgesDeferredMultiSite;
        result.edgesDeferredEndBubble += st.edgesDeferredEndBubble;
        result.edgesDeferredComplex   += st.edgesDeferredComplex;
    }

    // ------------------------------------------------------------------
    // Phase 2: apply (serial). appendHetAnchorPair grows the anchor store and
    // add_vertex must keep vertex_descriptor == AnchorId, so this cannot be
    // parallelized.
    // ------------------------------------------------------------------
    applyHetBubblePlans(graph, anchors, threadPlans, result);

    const auto t1 = std::chrono::steady_clock::now();
    result.elapsedSeconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    return result;
}
