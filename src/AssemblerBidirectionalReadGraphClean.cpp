/*******************************************************************************

AssemblerBidirectionalReadGraphClean.cpp
========================================

String-graph-style cleaning operations that work directly on the
BidirectionalReadGraph (BRG).

The BRG stores undirected edges (one per alignment, one vertex per physical
read).  String-graph cleaning algorithms (transitive reduction, tip cutting,
weak-arc removal) need directed arcs with extension lengths.

Strategy:
  1. Build a temporary in-memory directed-arc view from BRG edges, using
     computeMaHit2ArcLikeHifiasm to derive direction + extension lengths.
  2. Run the cleaning algorithms on the directed view (marking arcs deleted).
  3. Sync deletions back to BRG edge isDeleted flags.

This avoids creating a separate persistent StringGraph structure while
providing full hifiasm-parity cleaning semantics.

*******************************************************************************/

// Dinara.
#include "Assembler.hpp"
#include "timestamp.hpp"
#include "chrono.hpp"
#include "Reads.hpp"

// Standard library.
#include "algorithm.hpp"
#include <queue>
#include <numeric>

using namespace dinara;
using namespace std;



// ============================================================================
// Anonymous namespace: temporary directed-arc structures and helpers
// ============================================================================
namespace {

// A directed arc derived from a BRG edge.
// Stored in RC-twin pairs at consecutive indices: arcId ^ 1 is the twin.
struct BrgDirectedArc {
    uint32_t from;         // OrientedReadId::getValue()
    uint32_t to;           // OrientedReadId::getValue()
    uint32_t len;          // Non-overlap extension length (hifiasm asg_arc_len)
    uint32_t overlapLen;   // Overlap length (hifiasm p->ol)
    uint32_t brgEdgeId;    // Back-pointer to BidirectionalReadGraph::edges[]
    uint8_t  del;          // Soft-delete flag
};

// Result of computeMaHit2ArcLikeHifiasm.
struct ArcFromAlignmentResult {
    uint32_t fromVertex = uint32_t(invalidReadId);
    uint32_t toVertex = uint32_t(invalidReadId);
    uint32_t len = 0;
    uint32_t overlapLen = 0;
};

// Hifiasm-parity arc computation: alignment coordinates → directed arc.
// Duplicated from AssemblerStringGraph.cpp to keep this file self-contained.
inline bool computeMaHit2ArcForBrg(
    ReadId qn,
    ReadId tn,
    bool rev,
    uint32_t qs,
    uint32_t qe,
    uint32_t ql,
    uint32_t ts,
    uint32_t te,
    uint32_t tl,
    uint32_t maxHang,
    float maxHangRate,
    uint32_t minOverlapLen,
    ArcFromAlignmentResult& out)
{
    if (qe <= qs || te <= ts) return false;
    if (qe > ql || te > tl) return false;

    const int32_t tl5 = rev ? int32_t(tl) - int32_t(te) : int32_t(ts);
    const int32_t tl3 = rev ? int32_t(ts) : int32_t(tl) - int32_t(te);

    const int32_t qsI = int32_t(qs);
    const int32_t qeI = int32_t(qe);
    const int32_t qlI = int32_t(ql);
    const int32_t tSpan = int32_t(te) - int32_t(ts);
    const int32_t qSpan = qeI - qsI;
    if (qSpan <= 0 || tSpan <= 0) return false;

    const int32_t qRightHang = qlI - qeI;
    if (qRightHang < 0 || tl5 < 0 || tl3 < 0) return false;

    const int32_t ext5 = std::min(qsI, tl5);
    const int32_t ext3 = std::min(qRightHang, tl3);
    const int32_t alignedQ = qSpan + ext5 + ext3;
    const int32_t alignedT = tSpan + ext5 + ext3;

    if (ext5 > int32_t(maxHang) || ext3 > int32_t(maxHang)) return false;
    if (float(qSpan) < float(alignedQ) * maxHangRate) return false;
    if (float(tSpan) < float(alignedT) * maxHangRate) return false;

    // Containment checks.
    if (qsI <= tl5 && qRightHang <= tl3) return false;
    if (qsI >= tl5 && qRightHang >= tl3) return false;

    if (alignedQ < int32_t(minOverlapLen) || alignedT < int32_t(minOverlapLen)) return false;

    uint32_t uEnd = 0;
    uint32_t vEnd = 0;
    int32_t l = 0;
    if (qsI > tl5) {
        uEnd = 0;
        vEnd = rev ? 1U : 0U;
        l = qsI - tl5;
    } else {
        uEnd = 1;
        vEnd = rev ? 0U : 1U;
        l = (qlI - qeI) - tl3;
    }
    if (l < 0 || uint32_t(l) > ql) return false;

    out.fromVertex = OrientedReadId(qn, Strand(uEnd)).getValue();
    out.toVertex   = OrientedReadId(tn, Strand(vEnd)).getValue();
    out.len        = uint32_t(l);
    out.overlapLen = ql - uint32_t(l);
    return true;
}


// ============================================================================
// BrgDirectedView — temporary directed graph derived from BRG
// ============================================================================
class BrgDirectedView {
public:
    vector<BrgDirectedArc>          arcs;
    vector<vector<uint32_t>>        outgoing;  // outgoing[orientedReadId] → arc indices
    vector<vector<uint32_t>>        incoming;  // incoming[orientedReadId] → arc indices
    vector<uint8_t>                 readDeleted; // readDeleted[readId]
    uint32_t                        vertexCount = 0;  // 2 * readCount
    uint64_t                        nReads = 0;

    // Build directed arcs from BRG edges + alignment data.
    void buildFromBrg(
        const BidirectionalReadGraph& brg,
        const MemoryMapped::Vector<AlignmentData>& alignmentData,
        const Reads& reads);

    // Rebuild adjacency lists to exclude deleted arcs, sort outgoing by len.
    void rebuildAdjacency();

    // Sync arc deletions back to BRG edges.
    void syncDeletionsToBrg(BidirectionalReadGraph& brg) const;

    // --- Cleaning algorithms ---
    uint64_t transitiveReduce(uint32_t fuzz);
    uint64_t cutTips(uint32_t maxShortTipReads);
    uint64_t cutWeakArcs(
        uint32_t maxExtReads,
        double lenRatio,
        uint32_t minDiff,
        const MemoryMapped::Vector<AlignmentData>& alignmentData,
        const MemoryMapped::VectorOfVectors<uint32_t, uint32_t>& alignmentTable);
    uint64_t removeSingleNodeBubbles(uint32_t maxShortTipReads);
    uint64_t breakShortCycles(uint32_t maxCycleReads);
    uint64_t dropShortOverlaps(double dropRatio, uint32_t minOverlapLen);

private:
    // Helpers matching hifiasm semantics.
    uint32_t countOutNonDel(uint32_t v) const;
    bool firstOutNonDel(uint32_t v, uint32_t& arcIdOut) const;
    void symmetrizeArcDeletion();
    void deleteRead(ReadId readId);

    // Hifiasm parity: asg_is_utg_end
    int isUtgEnd(uint32_t v, uint32_t* nextVertex) const;
    // Hifiasm parity: asg_check_unambi1
    uint32_t checkUnambiguous1(uint32_t v) const;
    // Hifiasm parity: asg_topocut_aux
    int topocutAux(uint32_t v, int maxExtReads) const;
};


// -----------------------------------------------------------------------
// BrgDirectedView::buildFromBrg
// -----------------------------------------------------------------------
void BrgDirectedView::buildFromBrg(
    const BidirectionalReadGraph& brg,
    const MemoryMapped::Vector<AlignmentData>& alignmentData,
    const Reads& reads)
{
    nReads = brg.readCount();
    vertexCount = uint32_t(2 * nReads);

    readDeleted.assign(nReads, 0);
    arcs.clear();

    // Reserve: each BRG edge may produce 2 directed arcs (or 0 if rejected).
    arcs.reserve(2 * brg.edges.size());

    static constexpr uint32_t maSgMaxHang = 1000;
    static constexpr float    maSgMaxHangRate = 0.8f;
    static constexpr uint32_t maSgMinOverlapLen = 50;

    for (uint64_t edgeId = 0; edgeId < brg.edges.size(); ++edgeId) {
        const auto& edge = brg.edges[edgeId];
        if (edge.isDeleted) continue;

        const uint64_t alId = edge.alignmentId;
        if (alId >= alignmentData.size()) continue;

        const AlignmentData& ad = alignmentData[alId];
        const ReadId qn = ad.readIds[0];
        const ReadId tn = ad.readIds[1];
        const bool rev = !ad.isSameStrand;

        const uint32_t ql = uint32_t(reads.getReadRawSequenceLength(qn));
        const uint32_t tl = uint32_t(reads.getReadRawSequenceLength(tn));
        if (ql == 0 || tl == 0) continue;

        ArcFromAlignmentResult arc;
        const bool ok = computeMaHit2ArcForBrg(
            qn, tn, rev,
            ad.qs, ad.qe, ql,
            ad.ts, ad.te, tl,
            maSgMaxHang, maSgMaxHangRate, maSgMinOverlapLen,
            arc);
        if (!ok) continue;

        // Store arc pair at consecutive indices (arcId^1 is RC twin).
        BrgDirectedArc a;
        a.from       = arc.fromVertex;
        a.to         = arc.toVertex;
        a.len        = arc.len;
        a.overlapLen = arc.overlapLen;
        a.brgEdgeId  = uint32_t(edgeId);
        a.del        = 0;
        arcs.push_back(a);

        BrgDirectedArc b;
        b.from       = a.to ^ 1U;
        b.to         = a.from ^ 1U;
        b.len        = a.len;
        b.overlapLen = a.overlapLen;
        b.brgEdgeId  = uint32_t(edgeId);
        b.del        = 0;
        arcs.push_back(b);
    }

    // Build initial adjacency.
    rebuildAdjacency();
}


// -----------------------------------------------------------------------
// BrgDirectedView::rebuildAdjacency
// -----------------------------------------------------------------------
void BrgDirectedView::rebuildAdjacency()
{
    outgoing.assign(vertexCount, {});
    incoming.assign(vertexCount, {});

    // Mark arcs incident to deleted reads.
    if (!readDeleted.empty()) {
        for (uint64_t arcId = 0; arcId < arcs.size(); ++arcId) {
            auto& a = arcs[arcId];
            if (a.del) continue;
            const uint32_t fromRead = a.from >> 1U;
            const uint32_t toRead   = a.to >> 1U;
            if ((fromRead < readDeleted.size() && readDeleted[fromRead]) ||
                (toRead   < readDeleted.size() && readDeleted[toRead])) {
                a.del = 1;
            }
        }
    }

    for (uint64_t arcId = 0; arcId < arcs.size(); ++arcId) {
        const auto& a = arcs[arcId];
        if (a.del) continue;
        outgoing[a.from].push_back(uint32_t(arcId));
        incoming[a.to].push_back(uint32_t(arcId));
    }

    // Sort outgoing adjacency by len (required by transitive reduction).
    for (uint32_t v = 0; v < vertexCount; ++v) {
        std::stable_sort(outgoing[v].begin(), outgoing[v].end(),
            [&](uint32_t aId, uint32_t bId) {
                return arcs[aId].len < arcs[bId].len;
            });
    }
}


// -----------------------------------------------------------------------
// BrgDirectedView::syncDeletionsToBrg
// -----------------------------------------------------------------------
void BrgDirectedView::syncDeletionsToBrg(BidirectionalReadGraph& brg) const
{
    // An arc is deleted ⟹ mark its source BRG edge as deleted.
    // Also handle per-read deletion: mark all incident BRG edges.
    for (const auto& a : arcs) {
        if (a.del && a.brgEdgeId < brg.edges.size()) {
            brg.edges[a.brgEdgeId].isDeleted = 1;
        }
    }

    // Per-read deletion: mark all edges incident to deleted reads.
    for (ReadId readId = 0; readId < ReadId(nReads); ++readId) {
        if (readDeleted[readId]) {
            for (const uint32_t edgeId : brg.connectivity[readId]) {
                brg.edges[edgeId].isDeleted = 1;
            }
        }
    }
}


// -----------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------
uint32_t BrgDirectedView::countOutNonDel(uint32_t v) const
{
    uint32_t n = 0;
    for (const uint32_t arcId : outgoing[v]) {
        if (!arcs[arcId].del) ++n;
    }
    return n;
}

bool BrgDirectedView::firstOutNonDel(uint32_t v, uint32_t& arcIdOut) const
{
    for (const uint32_t arcId : outgoing[v]) {
        if (!arcs[arcId].del) {
            arcIdOut = arcId;
            return true;
        }
    }
    return false;
}

void BrgDirectedView::symmetrizeArcDeletion()
{
    DINARA_ASSERT((arcs.size() & 1ULL) == 0);
    for (uint64_t arcId = 0; arcId < arcs.size(); arcId += 2) {
        const uint8_t del = uint8_t(arcs[arcId].del | arcs[arcId ^ 1ULL].del);
        arcs[arcId].del = del;
        arcs[arcId ^ 1ULL].del = del;
    }
}

void BrgDirectedView::deleteRead(ReadId readId)
{
    if (readId < ReadId(readDeleted.size())) {
        readDeleted[readId] = 1;
    }
    const uint32_t v0 = uint32_t(readId) << 1U;
    const uint32_t v1 = v0 ^ 1U;
    for (const uint32_t arcId : outgoing[v0]) {
        arcs[arcId].del = 1;
    }
    for (const uint32_t arcId : outgoing[v1]) {
        arcs[arcId].del = 1;
    }
}

// Hifiasm parity: asg_is_utg_end.
// Examines v^1's outgoing arcs to determine unitig end type.
int BrgDirectedView::isUtgEnd(uint32_t v, uint32_t* nextVertex) const
{
    constexpr int ASG_ET_MERGEABLE  = 0;
    constexpr int ASG_ET_TIP        = 1;
    constexpr int ASG_ET_MULTI_OUT  = 2;
    constexpr int ASG_ET_MULTI_NEI  = 3;

    const uint32_t v1 = v ^ 1U;
    const uint32_t nv = countOutNonDel(v1);
    if (nv == 0) return ASG_ET_TIP;
    if (nv > 1)  return ASG_ET_MULTI_OUT;

    uint32_t arcId = 0;
    firstOutNonDel(v1, arcId);
    const uint32_t to = arcs[arcId].to;
    if (nextVertex) *nextVertex = to;

    const uint32_t w = to ^ 1U;
    const uint32_t nw = countOutNonDel(w);
    if (nw != 1) return ASG_ET_MULTI_NEI;
    return ASG_ET_MERGEABLE;
}

uint32_t BrgDirectedView::checkUnambiguous1(uint32_t v) const
{
    uint32_t next = std::numeric_limits<uint32_t>::max();
    uint32_t count = 0;
    for (const uint32_t arcId : outgoing[v]) {
        if (arcs[arcId].del) continue;
        next = arcs[arcId].to;
        if (++count > 1) break;
    }
    return (count == 1) ? next : std::numeric_limits<uint32_t>::max();
}

int BrgDirectedView::topocutAux(uint32_t v, int maxExtReads) const
{
    int nExt = 1;
    for (; nExt < maxExtReads && v != std::numeric_limits<uint32_t>::max(); ++nExt) {
        if (checkUnambiguous1(v ^ 1U) == std::numeric_limits<uint32_t>::max()) {
            --nExt;
            break;
        }
        v = checkUnambiguous1(v);
    }
    return nExt;
}


// ============================================================================
// Transitive Reduction
// ============================================================================
// Hifiasm parity: asg_arc_del_trans.
// For each vertex v, if v→w→x is reachable within len(v→x) + fuzz,
// and v→x is a direct arc, then v→x is transitive and can be deleted.
uint64_t BrgDirectedView::transitiveReduce(uint32_t fuzz)
{
    vector<uint8_t> mark(vertexCount, 0);
    uint64_t reduced = 0;

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (readDeleted[v >> 1U]) {
            for (const uint32_t arcId : outgoing[v]) {
                if (!arcs[arcId].del) { arcs[arcId].del = 1; ++reduced; }
            }
            continue;
        }

        const auto& out = outgoing[v];
        if (out.empty()) continue;

        // Mark all direct neighbors.
        for (const uint32_t arcId : out) {
            if (!arcs[arcId].del) mark[arcs[arcId].to] = 1;
        }

        // Find longest arc len (outgoing is sorted by len, last non-deleted).
        uint32_t longestLen = 0;
        for (auto it = out.rbegin(); it != out.rend(); ++it) {
            if (!arcs[*it].del) { longestLen = arcs[*it].len; break; }
        }
        const uint32_t L = longestLen + fuzz;

        // For each neighbor w of v, check if w reaches any other neighbor x of v
        // within the length bound.
        for (const uint32_t arcIdVw : out) {
            if (arcs[arcIdVw].del) continue;
            const uint32_t w = arcs[arcIdVw].to;
            if (mark[w] != 1) continue;

            for (const uint32_t arcIdWx : outgoing[w]) {
                if (arcs[arcIdWx].del) continue;
                if (arcs[arcIdWx].len + arcs[arcIdVw].len > L) break; // sorted by len
                const uint32_t x = arcs[arcIdWx].to;
                if (mark[x]) mark[x] = 2; // x is transitively reachable
            }
        }

        // Delete transitively reachable arcs.
        for (const uint32_t arcId : out) {
            if (arcs[arcId].del) continue;
            const uint32_t w = arcs[arcId].to;
            if (mark[w] == 2) {
                arcs[arcId].del = 1;
                ++reduced;
            }
            mark[w] = 0;
        }
    }

    if (reduced) {
        // Post-process: symmetrize + cleanup + remove multi/asymmetric arcs.
        symmetrizeArcDeletion();

        // Remove multi-arcs (keep one v→w, delete duplicates).
        {
            vector<uint32_t> count(vertexCount, 0);
            for (uint32_t v = 0; v < vertexCount; ++v) {
                for (auto it = outgoing[v].rbegin(); it != outgoing[v].rend(); ++it) {
                    if (arcs[*it].del) continue;
                    ++count[arcs[*it].to];
                }
                for (auto it = outgoing[v].rbegin(); it != outgoing[v].rend(); ++it) {
                    if (arcs[*it].del) continue;
                    if (--count[arcs[*it].to] != 0) {
                        arcs[*it].del = 1;
                    }
                }
            }
        }

        // Remove asymmetric arcs: delete u→v if (v^1)→(u^1) is absent.
        rebuildAdjacency();
        for (uint32_t u = 0; u < vertexCount; ++u) {
            for (const uint32_t arcId : outgoing[u]) {
                if (arcs[arcId].del) continue;
                const uint32_t vv = arcs[arcId].to;
                const uint32_t rcFrom = vv ^ 1U;
                const uint32_t rcTo   = u ^ 1U;
                bool found = false;
                for (const uint32_t rcArcId : outgoing[rcFrom]) {
                    if (arcs[rcArcId].del) continue;
                    if (arcs[rcArcId].to == rcTo) { found = true; break; }
                }
                if (!found) arcs[arcId].del = 1;
            }
        }

        rebuildAdjacency();
    }

    return reduced;
}


// ============================================================================
// Tip Cutting
// ============================================================================
// Hifiasm parity: asg_arc_cut_tips.
// Tips are short dead-end unitigs. A vertex v is a tip start if outdegree(v^1)==0
// (no incoming arcs). Walk the mergeable chain; if shorter than maxShortTipReads,
// delete all reads in the tip.
uint64_t BrgDirectedView::cutTips(uint32_t maxShortTipReads)
{
    if (maxShortTipReads == 0) return 0;

    constexpr int ASG_ET_MERGEABLE = 0;

    auto outdegree0 = [&](uint32_t v) -> bool {
        return countOutNonDel(v) == 0;
    };

    // Phase 1: gather candidate tips (sorted by unitig length, shortest first).
    vector<uint64_t> candidates;
    candidates.reserve(vertexCount / 8 + 16);

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (readDeleted[v >> 1U]) continue;
        if (!outdegree0(v ^ 1U)) continue; // tip start: no arcs into v

        uint32_t w = v;
        uint32_t kv = 1;
        uint32_t i = 0;
        for (i = 0; i < maxShortTipReads; ++i) {
            uint32_t nextVertex = 0;
            if (isUtgEnd(w ^ 1U, &nextVertex) != ASG_ET_MERGEABLE) break;
            w = nextVertex;
            ++kv;
        }
        if (i < maxShortTipReads) {
            candidates.push_back((uint64_t(kv) << 32) | uint64_t(v));
        }
    }

    std::sort(candidates.begin(), candidates.end());

    // Phase 2: cut tips shortest first, re-checking each candidate.
    uint64_t cutCount = 0;
    for (const uint64_t key : candidates) {
        const uint32_t v = uint32_t(key);
        if (readDeleted[v >> 1U]) continue;
        if (!outdegree0(v ^ 1U)) continue;

        vector<uint32_t> path;
        path.push_back(v);
        uint32_t w = v;
        uint32_t i = 0;
        for (i = 0; i < maxShortTipReads; ++i) {
            uint32_t nextVertex = 0;
            if (isUtgEnd(w ^ 1U, &nextVertex) != ASG_ET_MERGEABLE) break;
            w = nextVertex;
            path.push_back(w);
        }
        if (i < maxShortTipReads) {
            for (const uint32_t vv : path) {
                deleteRead(ReadId(vv >> 1U));
            }
            ++cutCount;
        }
    }

    if (cutCount) {
        symmetrizeArcDeletion();
        rebuildAdjacency();
    }
    return cutCount;
}


// ============================================================================
// Weak Arc Cutting
// ============================================================================
// Hifiasm parity: asg_arc_cut_weak (ONT path).
// For vertices with ≥2 outgoing arcs, arcs with overlapLen significantly
// less than the best are "weak". If the strong neighbor overlaps the weak
// neighbor (dedup check) and topology permits, delete the weak arc.
uint64_t BrgDirectedView::cutWeakArcs(
    uint32_t maxExtReads,
    double lenRatio,
    uint32_t minDiff,
    const MemoryMapped::Vector<AlignmentData>& alignmentData,
    const MemoryMapped::VectorOfVectors<uint32_t, uint32_t>& alignmentTable)
{
    if (maxExtReads == 0) return 0;

    auto hasAlignmentBetweenOrientedReads = [&](OrientedReadId a, OrientedReadId b) -> bool {
        const span<const uint32_t> section = alignmentTable[a.getValue()];
        auto it = std::lower_bound(
            section.begin(), section.end(), b,
            [&](uint32_t alId, const OrientedReadId& target) {
                const OrientedReadId other = alignmentData[alId].getOther(a);
                return other < target;
            });
        if (it == section.end()) return false;
        return alignmentData[*it].getOther(a) == b;
    };

    auto hasAnyAlignmentBetweenReads = [&](ReadId a, ReadId b) -> bool {
        for (Strand sa = 0; sa < 2; ++sa) {
            for (Strand sb = 0; sb < 2; ++sb) {
                if (hasAlignmentBetweenOrientedReads(OrientedReadId(a, sa), OrientedReadId(b, sb)))
                    return true;
            }
        }
        return false;
    };

    // Phase 1: gather weak-arc candidates.
    vector<uint64_t> candidates;
    candidates.reserve(1024);

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (readDeleted[v >> 1U]) continue;

        uint32_t maxOl = 0;
        uint32_t outCount = 0;
        uint32_t maxArcId = std::numeric_limits<uint32_t>::max();
        for (const uint32_t arcId : outgoing[v]) {
            const auto& a = arcs[arcId];
            if (a.del) continue;
            ++outCount;
            if (a.overlapLen > maxOl) {
                maxOl = a.overlapLen;
                maxArcId = arcId;
            }
        }
        if (outCount < 2) continue;
        if (maxArcId == std::numeric_limits<uint32_t>::max()) continue;

        const ReadId strongRead = ReadId(arcs[maxArcId].to >> 1U);
        for (const uint32_t arcId : outgoing[v]) {
            const auto& a = arcs[arcId];
            if (a.del || arcId == maxArcId) continue;
            if (a.overlapLen >= maxOl) continue;
            if (uint32_t(double(a.overlapLen)) + minDiff > maxOl) continue;
            if (double(a.overlapLen) > double(maxOl) * lenRatio) continue;

            const ReadId weakRead = ReadId(a.to >> 1U);
            if (!hasAnyAlignmentBetweenReads(strongRead, weakRead)) continue;

            candidates.push_back((uint64_t(a.overlapLen) << 32) | uint64_t(arcId));
        }
    }

    std::sort(candidates.begin(), candidates.end());

    // Phase 2: evaluate candidates and delete.
    uint64_t cut = 0;
    for (const uint64_t key : candidates) {
        const uint32_t arcId = uint32_t(key);
        if (arcId >= arcs.size() || arcs[arcId].del) continue;

        const uint32_t v  = arcs[arcId].from;
        const uint32_t to = arcs[arcId].to;
        const uint32_t w  = to ^ 1U;
        if (readDeleted[v >> 1U] || readDeleted[w >> 1U]) continue;

        const uint32_t kv = countOutNonDel(v);
        const uint32_t kw = countOutNonDel(w);
        if (kv <= 1 && kw <= 1) continue;

        const uint32_t twinId = arcId ^ 1U;
        if (twinId >= arcs.size()) continue;
        const auto& ve = arcs[arcId];
        const auto& we = arcs[twinId];
        const uint32_t mmOl = std::min(ve.overlapLen, we.overlapLen);

        // Check v-side best alternative.
        uint32_t vOlMax = 0;
        if (kv >= 2) {
            for (const uint32_t altId : outgoing[v]) {
                const auto& a = arcs[altId];
                if (a.del || a.to == ve.to) continue;
                if (a.overlapLen <= ve.overlapLen) continue;
                if (!hasAnyAlignmentBetweenReads(ReadId(a.to >> 1U), ReadId(ve.to >> 1U))) continue;
                vOlMax = std::max(vOlMax, a.overlapLen);
            }
            if (vOlMax == 0) continue;
            if (double(mmOl) > double(vOlMax) * lenRatio) continue;
            if (mmOl + minDiff > vOlMax) continue;
        }

        // Check w-side best alternative.
        uint32_t wOlMax = 0;
        if (kw >= 2) {
            for (const uint32_t altId : outgoing[w]) {
                const auto& a = arcs[altId];
                if (a.del || a.to == we.to) continue;
                if (a.overlapLen <= we.overlapLen) continue;
                if (!hasAnyAlignmentBetweenReads(ReadId(a.to >> 1U), ReadId(we.to >> 1U))) continue;
                wOlMax = std::max(wOlMax, a.overlapLen);
            }
            if (wOlMax == 0) continue;
            if (double(mmOl) > double(wOlMax) * lenRatio) continue;
            if (mmOl + minDiff > wOlMax) continue;
        }

        // Topology criterion.
        bool toDel = false;
        if (kv > 1 && kw > 1) {
            toDel = true;
        } else if (kw == 1) {
            if (topocutAux(w ^ 1U, int(maxExtReads)) < int(maxExtReads)) toDel = true;
        } else if (kv == 1) {
            if (topocutAux(v ^ 1U, int(maxExtReads)) < int(maxExtReads)) toDel = true;
        }

        if (toDel) {
            arcs[arcId].del = 1;
            arcs[twinId].del = 1;
            ++cut;
        }
    }

    if (cut) {
        symmetrizeArcDeletion();
        rebuildAdjacency();
    }
    return cut;
}


// ============================================================================
// Single-Node Bubble Removal
// ============================================================================
// Simplified hifiasm asg_arc_del_single_node_directly.
// Removes the weaker branch of a 2-branch bubble that reconverges in one step:
//   v → a → t
//   v → b → t
uint64_t BrgDirectedView::removeSingleNodeBubbles(uint32_t maxShortTipReads)
{
    uint64_t reduced = 0;

    auto incomingNonDelCount = [&](uint32_t v) -> uint32_t {
        uint32_t n = 0;
        for (const uint32_t arcId : incoming[v]) {
            if (!arcs[arcId].del) ++n;
        }
        return n;
    };

    auto hasIncomingFrom = [&](uint32_t v, uint32_t from) -> bool {
        for (const uint32_t arcId : incoming[v]) {
            if (arcs[arcId].del) continue;
            if (arcs[arcId].from == from) return true;
        }
        return false;
    };

    auto uniqueOutTo = [&](uint32_t v, uint32_t& to) -> bool {
        uint32_t arcId = 0;
        if (!firstOutNonDel(v, arcId)) return false;
        if (countOutNonDel(v) != 1) return false;
        to = arcs[arcId].to;
        return true;
    };

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (readDeleted[v >> 1U]) continue;

        // Need exactly 2 outgoing non-deleted arcs.
        array<uint32_t, 2> outArcIds{};
        uint32_t outCount = 0;
        for (const uint32_t arcId : outgoing[v]) {
            if (arcs[arcId].del) continue;
            if (outCount < 2) outArcIds[outCount] = arcId;
            ++outCount;
        }
        if (outCount != 2) continue;

        const uint32_t arcVa = outArcIds[0];
        const uint32_t arcVb = outArcIds[1];
        const uint32_t a = arcs[arcVa].to;
        const uint32_t b = arcs[arcVb].to;
        if (a == b) continue;

        if (incomingNonDelCount(a) != 1 || incomingNonDelCount(b) != 1) continue;
        if (!hasIncomingFrom(a, v) || !hasIncomingFrom(b, v)) continue;

        uint32_t tA = 0, tB = 0;
        if (!uniqueOutTo(a, tA) || !uniqueOutTo(b, tB) || tA != tB) continue;

        if (maxShortTipReads < 2) continue;

        // Score: prefer branch with more overlap.
        uint32_t arcAtoT = 0, arcBtoT = 0;
        firstOutNonDel(a, arcAtoT);
        firstOutNonDel(b, arcBtoT);
        const uint32_t scoreA = arcs[arcVa].overlapLen + arcs[arcAtoT].overlapLen;
        const uint32_t scoreB = arcs[arcVb].overlapLen + arcs[arcBtoT].overlapLen;

        const bool removeB = (scoreA >= scoreB);
        const uint32_t delArcVx   = removeB ? arcVb   : arcVa;
        const uint32_t delArcXtoT = removeB ? arcBtoT : arcAtoT;

        arcs[delArcVx].del = 1;
        arcs[delArcVx ^ 1U].del = 1;
        arcs[delArcXtoT].del = 1;
        arcs[delArcXtoT ^ 1U].del = 1;
        reduced += 4;
    }

    if (reduced) {
        symmetrizeArcDeletion();
        rebuildAdjacency();
    }
    return reduced;
}


// ============================================================================
// Short Cycle Breaking
// ============================================================================
// Break simple short cycles (all vertices 1-in/1-out) by deleting one edge.
uint64_t BrgDirectedView::breakShortCycles(uint32_t maxCycleReads)
{
    if (maxCycleReads == 0) return 0;

    auto outDeg = [&](uint32_t v) -> uint32_t { return countOutNonDel(v); };
    auto inDeg  = [&](uint32_t v) -> uint32_t {
        uint32_t n = 0;
        for (const uint32_t arcId : incoming[v]) {
            if (!arcs[arcId].del) ++n;
        }
        return n;
    };
    auto uniqueOutArc = [&](uint32_t v, uint32_t& arcIdOut) -> bool {
        if (outDeg(v) != 1) return false;
        return firstOutNonDel(v, arcIdOut);
    };

    vector<uint8_t> visited(vertexCount, 0);
    uint64_t broken = 0;

    for (uint32_t start = 0; start < vertexCount; ++start) {
        if (visited[start]) continue;
        if (readDeleted[start >> 1U]) { visited[start] = 1; continue; }
        if (outDeg(start) != 1 || inDeg(start) != 1) { visited[start] = 1; continue; }

        vector<uint32_t> path;
        uint32_t v = start;
        uint32_t steps = 0;
        bool cycle = false;

        while (true) {
            if (steps++ > maxCycleReads) break;
            if (visited[v] && v != start) break;
            path.push_back(v);

            uint32_t arcId = 0;
            if (!uniqueOutArc(v, arcId)) break;
            const uint32_t next = arcs[arcId].to;

            if (next == start) {
                cycle = true;
                path.push_back(next);
                break;
            }
            if (outDeg(next) != 1 || inDeg(next) != 1) break;
            v = next;
        }

        for (const uint32_t u : path) visited[u] = 1;

        if (!cycle || path.size() <= 2 || path.size() - 1 > maxCycleReads) continue;

        const uint32_t prev = path[path.size() - 2];
        uint32_t arcPrev = 0;
        if (!uniqueOutArc(prev, arcPrev)) continue;
        if (arcs[arcPrev].to != start) continue;

        arcs[arcPrev].del = 1;
        arcs[arcPrev ^ 1U].del = 1;
        broken += 2;
    }

    if (broken) {
        symmetrizeArcDeletion();
        rebuildAdjacency();
    }
    return broken;
}


// ============================================================================
// Drop Short Overlaps
// ============================================================================
// For each vertex with ≥2 outgoing arcs, drop arcs whose overlapLen is
// below dropRatio × maxOverlapLen.
uint64_t BrgDirectedView::dropShortOverlaps(double dropRatio, uint32_t minOverlapLen)
{
    uint64_t removed = 0;
    vector<uint32_t> arcIds;
    arcIds.reserve(64);

    for (uint32_t v = 0; v < vertexCount; ++v) {
        if (readDeleted[v >> 1U]) continue;
        if (outgoing[v].size() < 2) continue;

        arcIds.clear();
        uint32_t maxOl = 0;
        for (const uint32_t arcId : outgoing[v]) {
            if (arcs[arcId].del) continue;
            arcIds.push_back(arcId);
            maxOl = std::max(maxOl, arcs[arcId].overlapLen);
        }
        if (arcIds.size() < 2) continue;

        const uint32_t thres = std::max<uint32_t>(
            uint32_t(double(maxOl) * dropRatio + 0.499), minOverlapLen);

        // Sort by overlapLen descending.
        std::sort(arcIds.begin(), arcIds.end(), [&](uint32_t aId, uint32_t bId) {
            const auto& aa = arcs[aId];
            const auto& bb = arcs[bId];
            if (aa.overlapLen != bb.overlapLen) return aa.overlapLen > bb.overlapLen;
            if (aa.len != bb.len) return aa.len < bb.len;
            return aa.to < bb.to;
        });

        for (size_t i = 1; i < arcIds.size(); ++i) {
            const uint32_t arcId = arcIds[i];
            if (arcs[arcId].overlapLen < thres) {
                arcs[arcId].del = 1;
                arcs[arcId ^ 1U].del = 1;
                removed += 2;
            }
        }
    }

    if (removed) {
        symmetrizeArcDeletion();
        rebuildAdjacency();
    }
    return removed;
}


} // anonymous namespace



// ============================================================================
// Assembler public methods
// ============================================================================

// ---------------------------------------------------------------------------
// reduceBidirectionalReadGraphTransitive
// ---------------------------------------------------------------------------
uint64_t Assembler::reduceBidirectionalReadGraphTransitive(uint32_t gapFuzz)
{
    reads->checkReadsAreOpen();
    checkBidirectionalReadGraphIsOpen();
    checkAlignmentDataAreOpen();

    const auto t0 = steady_clock::now();
    cout << timestamp << "BRG transitive reduction begins (gapFuzz=" << gapFuzz << ")." << endl;

    BrgDirectedView view;
    view.buildFromBrg(bidirectionalReadGraph, alignmentData, *reads);

    cout << timestamp << "  Built directed view: " << view.arcs.size() << " directed arcs from "
         << bidirectionalReadGraph.edges.size() << " BRG edges." << endl;

    const uint64_t reduced = view.transitiveReduce(gapFuzz);
    cout << timestamp << "  Transitively reduced " << reduced << " arcs." << endl;

    view.syncDeletionsToBrg(bidirectionalReadGraph);

    cout << timestamp << "BRG transitive reduction complete in "
         << seconds(steady_clock::now() - t0) << " s." << endl;
    return reduced;
}


// ---------------------------------------------------------------------------
// cutBidirectionalReadGraphTips
// ---------------------------------------------------------------------------
uint64_t Assembler::cutBidirectionalReadGraphTips(uint32_t maxShortTipReads)
{
    reads->checkReadsAreOpen();
    checkBidirectionalReadGraphIsOpen();
    checkAlignmentDataAreOpen();

    const auto t0 = steady_clock::now();
    cout << timestamp << "BRG tip cutting begins (maxShortTipReads=" << maxShortTipReads << ")." << endl;

    BrgDirectedView view;
    view.buildFromBrg(bidirectionalReadGraph, alignmentData, *reads);
    const uint64_t cut = view.cutTips(maxShortTipReads);

    cout << timestamp << "  Cut " << cut << " tips." << endl;

    view.syncDeletionsToBrg(bidirectionalReadGraph);

    cout << timestamp << "BRG tip cutting complete in "
         << seconds(steady_clock::now() - t0) << " s." << endl;
    return cut;
}


// ---------------------------------------------------------------------------
// cutBidirectionalReadGraphWeakArcs
// ---------------------------------------------------------------------------
uint64_t Assembler::cutBidirectionalReadGraphWeakArcs(
    uint32_t maxExtReads,
    double lenRatio,
    uint32_t minDiff)
{
    reads->checkReadsAreOpen();
    checkBidirectionalReadGraphIsOpen();
    checkAlignmentDataAreOpen();

    const auto t0 = steady_clock::now();
    cout << timestamp << "BRG weak-arc cutting begins (maxExtReads=" << maxExtReads
         << ", lenRatio=" << lenRatio << ", minDiff=" << minDiff << ")." << endl;

    BrgDirectedView view;
    view.buildFromBrg(bidirectionalReadGraph, alignmentData, *reads);
    const uint64_t cut = view.cutWeakArcs(maxExtReads, lenRatio, minDiff,
                                           alignmentData, alignmentTable);

    cout << timestamp << "  Cut " << cut << " weak arcs." << endl;

    view.syncDeletionsToBrg(bidirectionalReadGraph);

    cout << timestamp << "BRG weak-arc cutting complete in "
         << seconds(steady_clock::now() - t0) << " s." << endl;
    return cut;
}


// ---------------------------------------------------------------------------
// cleanBidirectionalReadGraphInitial
// ---------------------------------------------------------------------------
// Initial cleaning pass: transitive reduction + tip cutting.
// Matches hifiasm Steps 9–11 pipeline semantics.
void Assembler::cleanBidirectionalReadGraphInitial(
    uint32_t gapFuzz,
    uint32_t maxShortTipReads)
{
    reads->checkReadsAreOpen();
    checkBidirectionalReadGraphIsOpen();
    checkAlignmentDataAreOpen();

    const auto t0 = steady_clock::now();
    cout << timestamp << "BRG initial cleaning begins (gapFuzz=" << gapFuzz
         << ", maxShortTipReads=" << maxShortTipReads << ")." << endl;

    // Build the directed view once and run all cleaning stages sequentially.
    BrgDirectedView view;
    view.buildFromBrg(bidirectionalReadGraph, alignmentData, *reads);

    cout << timestamp << "  Built directed view: " << view.arcs.size() << " directed arcs from "
         << bidirectionalReadGraph.edges.size() << " BRG edges." << endl;

    // Step 10: Transitive reduction.
    const uint64_t reduced = view.transitiveReduce(gapFuzz);
    cout << timestamp << "  Transitively reduced " << reduced << " arcs." << endl;

    // Step 11: Cut short tips.
    if (maxShortTipReads > 0) {
        const uint64_t cut = view.cutTips(maxShortTipReads);
        cout << timestamp << "  Cut " << cut << " tips." << endl;
    }

    // Sync all deletions back to BRG.
    view.syncDeletionsToBrg(bidirectionalReadGraph);

    // Count surviving edges.
    uint64_t survivingEdges = 0;
    for (uint64_t i = 0; i < bidirectionalReadGraph.edges.size(); ++i) {
        if (!bidirectionalReadGraph.edges[i].isDeleted) ++survivingEdges;
    }

    cout << timestamp << "BRG initial cleaning complete: "
         << survivingEdges << " / " << bidirectionalReadGraph.edges.size()
         << " edges survive (" << seconds(steady_clock::now() - t0) << " s)." << endl;
}


// ---------------------------------------------------------------------------
// cleanBidirectionalReadGraphIterative
// ---------------------------------------------------------------------------
// Iterative cleaning: cycles + bubbles + tips + overlap-ratio drops over
// multiple rounds with increasing stringency.
// Matches hifiasm ul_clean_gfa semantics.
void Assembler::cleanBidirectionalReadGraphIterative(
    uint32_t cleanRounds,
    double minDropRate,
    double maxDropRate,
    uint32_t maxShortTipReads)
{
    reads->checkReadsAreOpen();
    checkBidirectionalReadGraphIsOpen();
    checkAlignmentDataAreOpen();

    if (cleanRounds == 0) return;
    if (minDropRate <= 0. || maxDropRate <= 0. || minDropRate > maxDropRate || maxDropRate >= 1.) {
        throw runtime_error("cleanBidirectionalReadGraphIterative: invalid drop-rate range.");
    }

    const auto t0 = steady_clock::now();
    cout << timestamp << "BRG iterative cleaning begins: " << cleanRounds
         << " rounds, drop ratio " << minDropRate << " → " << maxDropRate << "." << endl;

    const double step = (cleanRounds == 1) ? 0. : (maxDropRate - minDropRate) / double(cleanRounds - 1);
    double dropRatio = minDropRate;

    for (uint32_t round = 0; round < cleanRounds; ++round, dropRatio += step) {
        if (dropRatio > maxDropRate) dropRatio = maxDropRate;

        cout << timestamp << "  Round " << (round + 1) << "/" << cleanRounds
             << " (dropRatio=" << dropRatio << ")." << endl;

        // Build a fresh directed view from current BRG state.
        BrgDirectedView view;
        view.buildFromBrg(bidirectionalReadGraph, alignmentData, *reads);

        // Pre-clean: cycles + bubbles + tips.
        uint64_t totalCycles = 0, totalBubbles = 0, totalTips = 0;
        for (uint32_t iter = 0; iter < 10; ++iter) {
            const uint64_t cycles  = view.breakShortCycles(100);
            const uint64_t bubbles = view.removeSingleNodeBubbles(maxShortTipReads);
            const uint64_t tips    = view.cutTips(maxShortTipReads);
            totalCycles  += cycles;
            totalBubbles += bubbles;
            totalTips    += tips;
            if (cycles == 0 && bubbles == 0 && tips == 0) break;
        }
        cout << timestamp << "    Pre-clean: broke " << totalCycles << " cycles, removed "
             << totalBubbles << " bubble arcs, cut " << totalTips << " tips." << endl;

        // Drop short overlaps.
        const uint64_t dropped = view.dropShortOverlaps(dropRatio, 0);
        cout << timestamp << "    Dropped " << dropped << " short-overlap arcs." << endl;

        // Post-clean: cycles + bubbles + tips again.
        uint64_t totalCycles2 = 0, totalBubbles2 = 0, totalTips2 = 0;
        for (uint32_t iter = 0; iter < 10; ++iter) {
            const uint64_t cycles  = view.breakShortCycles(100);
            const uint64_t bubbles = view.removeSingleNodeBubbles(maxShortTipReads);
            const uint64_t tips    = view.cutTips(maxShortTipReads);
            totalCycles2  += cycles;
            totalBubbles2 += bubbles;
            totalTips2    += tips;
            if (cycles == 0 && bubbles == 0 && tips == 0) break;
        }
        cout << timestamp << "    Post-clean: broke " << totalCycles2 << " cycles, removed "
             << totalBubbles2 << " bubble arcs, cut " << totalTips2 << " tips." << endl;

        // Sync this round's deletions back to BRG.
        view.syncDeletionsToBrg(bidirectionalReadGraph);
    }

    // Count surviving edges.
    uint64_t survivingEdges = 0;
    for (uint64_t i = 0; i < bidirectionalReadGraph.edges.size(); ++i) {
        if (!bidirectionalReadGraph.edges[i].isDeleted) ++survivingEdges;
    }

    cout << timestamp << "BRG iterative cleaning complete: "
         << survivingEdges << " / " << bidirectionalReadGraph.edges.size()
         << " edges survive (" << seconds(steady_clock::now() - t0) << " s)." << endl;
}
