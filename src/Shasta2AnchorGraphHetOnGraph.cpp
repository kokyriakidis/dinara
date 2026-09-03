// Detect abPOA local topology on the Shasta2AnchorGraph and turn it into new
// het/hom anchors.
//
// For each anchor-graph edge A -> B whose two-sided read coverage is at least
// minCommonForHet, we run an abPOA MSA over the reads' inter-anchor sequences.
// The MSA column matrix (msa_base[row][col], row order == anchorPair
// orientedReadIds because sort_input_seq=0) is a materialization of the POA
// graph in rank order. We read that topology in the anchor graph's own
// vocabulary:
//
//   * A column where every covered read carries the SAME non-gap base is a
//     strict convergence point. The edge's own endpoints A and B are the
//     outer convergence points by construction (all reads share both).
//
//   * A maximal run of columns between two convergence points where reads
//     disagree is a candidate site. Reads are grouped by their exact
//     subsequence over the run: each distinct group is an allele. A group
//     that clears the per-allele support floor is a passing allele; a run is
//     a REAL site only when >=2 alleles pass. A read whose subsequence over
//     the run is empty (it skips the run: a deletion allele) needs no new
//     anchor -- see below.
//
// This unifies SNPs (a one-column site) and indels (a multi-column site, one
// side long / one side empty) with no size special-case.
//
// What this file does NOT do: wire anything into the graph. Detection only
// appends a new anchor per passing non-deletion allele
// (Shasta2Anchors::appendHetAnchorPair) for each real site; it never touches
// the anchor graph or disables anything. The deletion allele and any read not
// in a passing allele (including every read on edges with no real site at
// all) get no new anchor -- they simply have nothing inserted at that point
// in their journey, so after the caller rebuilds journeys and the anchor
// graph from scratch (Shasta2Journeys::rebuildAfterNewAnchors, then a fresh
// Shasta2AnchorGraph), such a read naturally takes the direct flank->flank
// edge. This replaces an earlier design that manually wired flank->arm->flank
// (+ RC) edges in place: once new anchors just become ordinary entries in
// each read's journey, the (already-verified) journey->anchor-graph builder
// produces the correct topology for any number of sites per edge with no
// special-casing, and reads outside every passing allele get the "direct"
// edge for free instead of needing an explicit isDeletion/non-member case.
//
// A run that fails the >=2-passing-allele floor is NOT a veto on the rest of
// the edge -- it is simply not a site (absorbed into the surrounding
// homozygous background), and every other run on the edge is evaluated
// independently. Runs closer together than minConvergentGap() are merged
// before this analysis, to absorb abPOA-fragmented single variants (see its
// comment for the empirical motivation).
//
// Position-tie safety: a new anchor's stored position is (raw read position
// of its first base) + hetAnchorKHalf(); the edge's own flank anchors store
// (marker start + k/2). Shasta2AnchorPair::assertNoNegativeOffsets crashes if
// two anchors adjacent in some read's rebuilt journey have non-increasing
// positions for that read, so a new anchor's position must never tie or
// cross its neighbors'. Requiring the surviving chain of real sites to be
// strictly interior (first site doesn't touch column 0, last doesn't touch
// the last column -- see analyzeEdgeMsa) guarantees at least one convergent
// column of margin on each side, which is enough: Shasta2AnchorPair::get's
// inter-anchor sequence runs from anchorA's own position to anchorB's own
// position exclusive, so a site's raw position is always > anchorA's own
// position (by at least hetAnchorKHalf) and, given the required trailing
// margin, always < anchorB's own position for every member read.

#include "Shasta2AnchorGraph.hpp"

#include "Base.hpp"
#include "Reads.hpp"
#include "Shasta2AnchorPair.hpp"
#include "Shasta2Anchors.hpp"
#include "WindowIntervalPoa.hpp"   // IpoaAbHandle

#include <boost/graph/iteration_macros.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <mutex>
#include <thread>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// Investigation aid (DINARA_HET_MULTISITE_DEBUG=1): tallies the shape of
// multi-run edges (how many separate runs they contain, how far apart, and
// why a run failed to become a real site) so the merge/absorb policy can be
// tuned against what the data actually looks like. Thread-safe (mutex-
// protected); negligible overhead when off.
bool multiSiteDebugEnabled()
{
    static const bool on = (std::getenv("DINARA_HET_MULTISITE_DEBUG") != nullptr);
    return on;
}
std::mutex multiSiteDebugMutex;
std::map<uint64_t, uint64_t> multiSiteRunCountHistogram;   // runs.size() -> count
std::map<uint64_t, uint64_t> multiSiteGapHistogram;        // gap (in columns) between consecutive runs, bucketed
void recordMultiSiteShape(const vector<pair<uint64_t, uint64_t>>& runs)
{
    std::lock_guard<std::mutex> lock(multiSiteDebugMutex);
    multiSiteRunCountHistogram[runs.size()]++;
    for(uint64_t i = 0; i + 1 < runs.size(); i++) {
        const uint64_t gap = runs[i + 1].first - runs[i].second;
        // Bucket gaps geometrically so the histogram stays small.
        uint64_t bucket = 0;
        uint64_t g = gap;
        while(g > 0) { bucket++; g >>= 1; }
        multiSiteGapHistogram[bucket]++;
    }
}

// For a run (on a multi-run edge) that failed the >=2-passing-alleles check:
// how many distinct exact-subsequence groups it split into and how dominant
// the largest group was, to tell apart "genuinely too many alleles" from "one
// real consensus plus scattered noise elsewhere in the window".
std::map<uint64_t, uint64_t> multiSiteFailGroupCountHistogram;       // distinct alleles -> occurrences
std::map<uint64_t, uint64_t> multiSiteFailLargestGroupFracHistogram; // floor(10*largest/total) -> occurrences
void recordMultiSiteFailure(uint64_t numGroups, uint64_t largestGroupSize, uint64_t total)
{
    std::lock_guard<std::mutex> lock(multiSiteDebugMutex);
    multiSiteFailGroupCountHistogram[numGroups]++;
    const uint64_t frac10 = (total > 0) ? (largestGroupSize * 10) / total : 0;
    multiSiteFailLargestGroupFracHistogram[frac10]++;
}

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

// Minimum number of convergent columns required between two candidate sites
// to keep them as separate sites; runs separated by fewer columns are merged
// into one combined run instead. Empirically motivated: on a real test
// region, over half of all inter-site gaps on multi-run edges were 1-3
// columns, and edges with dozens of "sites" a few columns apart are far more
// consistent with one messy, indel-shifted divergent region (which abPOA's
// strict all-rows-identical convergence test fragments) than with that many
// genuinely independent heterozygous sites. Not a biological law -- a
// tunable heuristic; override with DINARA_HET_MIN_CONVERGENT_GAP.
uint64_t minConvergentGap()
{
    static const uint64_t gap = [] () -> uint64_t {
        const char* e = std::getenv("DINARA_HET_MIN_CONVERGENT_GAP");
        if(e == nullptr) return 10;
        return uint64_t(std::atoll(e));
    }();
    return gap;
}

// One new anchor to create: an allele's member list ((read, raw position of
// the allele's first base) pairs). The deletion allele and any allele that
// doesn't clear the support floor need no anchor and never appear here.
using ArmMembers = vector<pair<OrientedReadId, uint32_t>>;

// Everything found on one edge: every real site's non-deletion alleles,
// flattened. Site boundaries do not matter beyond this function -- creating
// each anchor and letting the caller rebuild journeys from scratch places it
// correctly regardless of how many other sites share the edge.
struct EdgeHetSites {
    vector<ArmMembers> arms;
    uint64_t siteCount = 0;   // number of real (>=2 passing alleles) sites found
};

// Result of analyzing one edge's MSA matrix.
enum class AnalyzeResult {
    NoBubble,        // no real site found (fully homozygous, or every
                     // non-convergent run was noise below the allele floor)
    Planned,         // >=1 real site found; outSites has the new anchors
    DeferEndBubble,  // the surviving chain of real sites would touch a span end
    DeferComplex     // reserved, currently unused
};

// Analyze the row-major MSA matrix (rows == anchorPair orientedReadIds order,
// via buildEdgeMsa) for one edge. `msa[r]` is row r's aligned bases (value
// 0..3, or 4 for gap), all of length ncols. `posA[r]` is read r's absolute
// base position of anchorA (so an aligned column maps back to a raw read
// position by counting that row's non-gap bases up to the column). `total`
// is the row count. See the file header for the overall algorithm and the
// position-tie safety argument behind the strictly-interior check below.
AnalyzeResult analyzeEdgeMsa(
    const vector<vector<uint8_t>>& msa,
    const vector<uint32_t>& posA,
    const vector<OrientedReadId>& rowRead,
    uint64_t total,
    EdgeHetSites& outSites)
{
    outSites.arms.clear();
    outSites.siteCount = 0;
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

    // Find maximal runs of non-convergent columns (candidate sites).
    vector<pair<uint64_t, uint64_t>> rawRuns;   // [begin, end) column ranges
    {
        uint64_t c = 0;
        while(c < ncols) {
            if(isConvergent(c)) { c++; continue; }
            uint64_t begin = c;
            while(c < ncols && !isConvergent(c)) c++;
            rawRuns.push_back({begin, c});
        }
    }
    if(rawRuns.empty()) return AnalyzeResult::NoBubble;

    // Merge runs separated by fewer than minConvergentGap() convergent
    // columns into one combined run -- see minConvergentGap()'s comment.
    vector<pair<uint64_t, uint64_t>> runs;
    runs.push_back(rawRuns.front());
    const uint64_t gapFloor = minConvergentGap();
    for(uint64_t i = 1; i < rawRuns.size(); i++) {
        const uint64_t gap = rawRuns[i].first - runs.back().second;
        if(gap < gapFloor) {
            runs.back().second = rawRuns[i].second;
        } else {
            runs.push_back(rawRuns[i]);
        }
    }

    if(multiSiteDebugEnabled() && runs.size() > 1) recordMultiSiteShape(runs);

    // Prefix non-gap counts per row, computed once, so converting any column
    // into a raw read position is O(1) regardless of how many runs there are
    // (avoids O(runs * ncols) work per row on many-site edges).
    vector<vector<uint32_t>> prefixNonGap(total, vector<uint32_t>(ncols + 1, 0));
    for(uint64_t r = 0; r < total; r++) {
        for(uint64_t col = 0; col < ncols; col++) {
            prefixNonGap[r][col + 1] = prefixNonGap[r][col] + (msa[r][col] < 4 ? 1 : 0);
        }
    }
    auto rawPositionAtColumn = [&](uint64_t r, uint64_t col) -> uint32_t {
        return posA[r] + prefixNonGap[r][col];
    };

    // Analyze each (merged) run independently. A run that fails the >=2-
    // passing-allele floor is simply not a site -- absorbed into the
    // surrounding background, not a veto -- so every other run on this edge
    // is unaffected (see file header).
    const uint64_t floor = alleleFloor(total);
    vector<pair<uint64_t, uint64_t>> passingRuns;   // [begin, end) of real sites only
    for(uint64_t ri = 0; ri < runs.size(); ri++) {
        const auto [rbegin, rend] = runs[ri];

        // Group rows by their exact subsequence (non-gap bases) over
        // [rbegin, rend). Empty key == deletion allele.
        std::map<vector<uint8_t>, vector<uint64_t>> groups;
        for(uint64_t r=0; r<total; r++) {
            vector<uint8_t> allele;
            for(uint64_t col=rbegin; col<rend; col++) {
                const uint8_t b = msa[r][col];
                if(b < 4) allele.push_back(b);
            }
            groups[allele].push_back(r);
        }

        uint64_t passing = 0;
        vector<ArmMembers> siteArms;
        for(const auto& [allele, rows] : groups) {
            if(rows.size() < floor) continue;
            passing++;
            if(allele.empty()) continue;   // deletion allele: no anchor needed
            ArmMembers arm;
            arm.reserve(rows.size());
            for(const uint64_t r: rows) {
                // Pin the arm anchor at the first base of the allele on this
                // read: the raw position of its first non-gap base at/after
                // the run start, which for a member of a non-empty allele is
                // exactly that allele's first base.
                arm.push_back({rowRead[r], rawPositionAtColumn(r, rbegin)});
            }
            siteArms.push_back(std::move(arm));
        }

        if(passing < 2) {
            if(multiSiteDebugEnabled() && runs.size() > 1) {
                uint64_t largest = 0;
                for(const auto& [allele, rows] : groups) {
                    if(rows.size() > largest) largest = rows.size();
                }
                recordMultiSiteFailure(groups.size(), largest, total);
            }
            continue;
        }

        passingRuns.push_back({rbegin, rend});
        outSites.siteCount++;
        for(ArmMembers& arm: siteArms) outSites.arms.push_back(std::move(arm));
    }

    if(passingRuns.empty()) return AnalyzeResult::NoBubble;

    // The chain of real sites must be strictly interior: the first can't
    // start at column 0, the last can't end at the last column. See the file
    // header for why this is required for position-tie safety, not just a
    // scope choice.
    if(passingRuns.front().first == 0 || passingRuns.back().second == ncols) {
        outSites.arms.clear();
        outSites.siteCount = 0;
        return AnalyzeResult::DeferEndBubble;
    }

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



HetOnGraphResult dinara::transcribeHetBubbles(
    const Shasta2AnchorGraph& graph,
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
    // Phase 1: detection (parallel, read-only). Each thread analyzes edges
    // off a shared atomic counter and appends this edge's real-site alleles
    // to a private list. The graph is never touched -- see the file header
    // for why no wiring happens here.
    // ------------------------------------------------------------------
    std::atomic<uint64_t> nextEdge{0};
    vector<vector<ArmMembers>> threadArms(threadCount);
    vector<HetOnGraphResult> threadStats(threadCount);

    auto planner = [&](uint64_t threadId) {
        HetOnGraphResult& st = threadStats[threadId];
        vector<ArmMembers>& arms = threadArms[threadId];

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

            EdgeHetSites sites;
            const AnalyzeResult ar = analyzeEdgeMsa(msa, posA, rowRead, rowRead.size(), sites);
            switch(ar) {
            case AnalyzeResult::NoBubble:
                break;
            case AnalyzeResult::DeferEndBubble:
                st.edgesDeferredEndBubble++;
                break;
            case AnalyzeResult::DeferComplex:
                st.edgesDeferredComplex++;
                break;
            case AnalyzeResult::Planned:
                st.edgesPlanned++;
                if(sites.siteCount > 1) st.edgesPlannedMultiSite++;
                st.sitesTranscribed += sites.siteCount;
                for(ArmMembers& arm: sites.arms) arms.push_back(std::move(arm));
                break;
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
        result.edgesPlannedMultiSite  += st.edgesPlannedMultiSite;
        result.sitesTranscribed       += st.sitesTranscribed;
        result.edgesDeferredEndBubble += st.edgesDeferredEndBubble;
        result.edgesDeferredComplex   += st.edgesDeferredComplex;
    }

    if(multiSiteDebugEnabled()) {
        std::lock_guard<std::mutex> lock(multiSiteDebugMutex);
        cout << "Multi-run edge shape: run-count histogram (runs.size() -> edges):" << endl;
        for(const auto& [runCount, edges]: multiSiteRunCountHistogram) {
            cout << "  " << runCount << " runs: " << edges << " edges" << endl;
        }
        cout << "Multi-run edge shape: inter-run gap histogram (log2 bucket of column gap -> occurrences):" << endl;
        for(const auto& [bucket, occurrences]: multiSiteGapHistogram) {
            const uint64_t lo = (bucket == 0) ? 0 : (1ULL << (bucket - 1));
            const uint64_t hi = (1ULL << bucket) - 1;
            cout << "  [" << lo << "-" << hi << "] columns: " << occurrences << " gaps" << endl;
        }
        cout << "Absorbed runs (merged run, <2 passing alleles): distinct-allele-group histogram:" << endl;
        for(const auto& [numGroups, occurrences]: multiSiteFailGroupCountHistogram) {
            cout << "  " << numGroups << " distinct groups: " << occurrences << " (run,edge) absorptions" << endl;
        }
        cout << "Absorbed runs: largest-group-fraction histogram (floor(10*largest/total) -> occurrences):" << endl;
        for(const auto& [frac10, occurrences]: multiSiteFailLargestGroupFracHistogram) {
            cout << "  largest group covers [" << frac10*10 << "-" << (frac10+1)*10 << "]% of rows: " << occurrences << " absorptions" << endl;
        }
    }

    // ------------------------------------------------------------------
    // Phase 2: apply (serial). appendHetAnchorPair grows the anchor store, so
    // this cannot be parallelized. No graph mutation: the caller rebuilds
    // journeys and the anchor graph from scratch once every new anchor has
    // been appended (see file header).
    // ------------------------------------------------------------------
    for(vector<ArmMembers>& arms: threadArms) {
        for(ArmMembers& arm: arms) {
            if(arm.size() < 2) continue;   // appendHetAnchorPair requires >= 2
            anchors.appendHetAnchorPair(arm);
            result.hetAnchorsCreated++;
        }
    }

    const auto t1 = std::chrono::steady_clock::now();
    result.elapsedSeconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    return result;
}
