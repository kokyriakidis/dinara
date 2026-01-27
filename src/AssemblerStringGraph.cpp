// Dinara.
#include "Assembler.hpp"
#include "timestamp.hpp"
#include "Reads.hpp"

using namespace dinara;
using namespace std;

namespace {
    struct ArcFromAlignmentResult {
        uint32_t fromVertex = uint32_t(invalidReadId);
        uint32_t toVertex = uint32_t(invalidReadId);
        uint32_t len = 0;
        uint32_t overlapLen = 0;
    };

    inline bool computeMaHit2ArcLikeHifiasm(
        ReadId qn,
        ReadId tn,
        bool rev,
        uint32_t qs,
        uint32_t qe,
        uint32_t ql,
        uint32_t ts,
        uint32_t te,
        uint32_t tl,
        ArcFromAlignmentResult& out)
    {
        if (qe <= qs || te <= ts) {
            return false;
        }
        if (qe > ql || te > tl) {
            return false;
        }

        const int32_t tl5 = rev ? int32_t(tl) - int32_t(te) : int32_t(ts);
        const int32_t tl3 = rev ? int32_t(ts) : int32_t(tl) - int32_t(te);

        const int32_t qsI = int32_t(qs);
        const int32_t qeI = int32_t(qe);
        const int32_t qlI = int32_t(ql);

        // Mirrors hifiasm's end selection:
        //   if (qs > tl5): query-to-target overlap => uEnd=0, vEnd=rev, l=qs-tl5
        //   else:          target-to-query overlap => uEnd=1, vEnd=!rev, l=(ql-qe)-tl3
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
        if (l < 0 || uint32_t(l) > ql) {
            return false;
        }

        const uint32_t from = OrientedReadId(qn, Strand(uEnd)).getValue();
        const uint32_t to = OrientedReadId(tn, Strand(vEnd)).getValue();
        out.fromVertex = from;
        out.toVertex = to;
        out.len = uint32_t(l);
        out.overlapLen = ql - uint32_t(l);
        return true;
    }
}


void Assembler::createStringGraphUsingSelectedAlignments(const vector<bool>& keepAlignment)
{
    checkAlignmentDataAreOpen();
    reads->checkReadsAreOpen();

    if (keepAlignment.size() != alignmentData.size()) {
        throw runtime_error("createStringGraphUsingSelectedAlignments: keepAlignment size mismatch.");
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t vertexCount = 2 * readCount;

    // Recreate from scratch.
    stringGraph.remove();
    stringGraph.arcs.createNew(largeDataName("StringGraphArcs"), largeDataPageSize);
    stringGraph.outgoing.createNew(largeDataName("StringGraphOutgoing"), largeDataPageSize);
    stringGraph.incoming.createNew(largeDataName("StringGraphIncoming"), largeDataPageSize);
    stringGraph.readDeleted.createNew(largeDataName("StringGraphReadDeleted"), largeDataPageSize);
    stringGraph.readDeleted.reserveAndResize(readCount);
    std::fill(stringGraph.readDeleted.begin(), stringGraph.readDeleted.end(), uint8_t(0));

    // Reserve a lower bound to reduce reallocations (each kept overlap yields two arcs).
    const size_t keptCount = size_t(count(keepAlignment.begin(), keepAlignment.end(), true));
    stringGraph.arcs.reserve(2 * keptCount);

    for (uint64_t alignmentId = 0; alignmentId < alignmentData.size(); ++alignmentId) {
        if (!keepAlignment[alignmentId]) {
            continue;
        }
        const AlignmentData& ad = alignmentData[alignmentId];

        const ReadId qn = ad.readIds[0];
        const ReadId tn = ad.readIds[1];
        const bool rev = !ad.isSameStrand;

        uint32_t ql = 0;
        uint32_t tl = 0;
        if (validReadIntervals.empty()) {
            ql = uint32_t(reads->getReadRawSequenceLength(qn));
            tl = uint32_t(reads->getReadRawSequenceLength(tn));
        } else {
            if (qn >= validReadIntervals.size() || tn >= validReadIntervals.size()) continue;
            const auto& rq = validReadIntervals[qn];
            const auto& rt = validReadIntervals[tn];
            if (rq.isDeleted || rt.isDeleted) {
                if (qn < readCount) stringGraph.readDeleted[qn] = 1;
                if (tn < readCount) stringGraph.readDeleted[tn] = 1;
                continue;
            }
            ql = rq.end - rq.start;
            tl = rt.end - rt.start;
        }
        if (ql == 0 || tl == 0) continue;

        ArcFromAlignmentResult arc;
        if (!computeMaHit2ArcLikeHifiasm(qn, tn, rev, ad.qs, ad.qe, ql, ad.ts, ad.te, tl, arc)) {
            continue;
        }

        StringGraphArc a;
        a.from = arc.fromVertex;
        a.to = arc.toVertex;
        a.len = arc.len;
        a.overlapLen = arc.overlapLen;
        a.alignmentId = alignmentId;
        a.del = 0;
        stringGraph.arcs.push_back(a);

        // Reverse-complement twin: (to^1) -> (from^1).
        StringGraphArc b = a;
        b.from = a.to ^ 1U;
        b.to = a.from ^ 1U;
        stringGraph.arcs.push_back(b);
    }

    stringGraph.unreserve();

    // Build outgoing adjacency (like hifiasm asg.idx/asg_arc_a).
    stringGraph.outgoing.beginPass1(vertexCount);
    stringGraph.incoming.beginPass1(vertexCount);
    for (uint64_t arcId = 0; arcId < stringGraph.arcs.size(); ++arcId) {
        const StringGraphArc& a = stringGraph.arcs[arcId];
        if (a.del) continue;
        DINARA_ASSERT(a.from < vertexCount);
        stringGraph.outgoing.incrementCount(a.from);
        stringGraph.incoming.incrementCount(a.to);
    }
    stringGraph.outgoing.beginPass2();
    stringGraph.incoming.beginPass2();
    for (uint64_t arcId = 0; arcId < stringGraph.arcs.size(); ++arcId) {
        const StringGraphArc& a = stringGraph.arcs[arcId];
        if (a.del) continue;
        stringGraph.outgoing.store(a.from, uint32_t(arcId));
        stringGraph.incoming.store(a.to, uint32_t(arcId));
    }
    stringGraph.outgoing.endPass2();
    stringGraph.incoming.endPass2();
    stringGraph.unreserve();
}



void Assembler::accessStringGraph()
{
    stringGraph.arcs.accessExistingReadOnly(largeDataName("StringGraphArcs"));
    stringGraph.outgoing.accessExistingReadOnly(largeDataName("StringGraphOutgoing"));
    try {
        stringGraph.incoming.accessExistingReadOnly(largeDataName("StringGraphIncoming"));
    } catch(...) {
        // Backward compatibility: older datasets may not have incoming adjacency persisted.
    }
    try {
        stringGraph.readDeleted.accessExistingReadOnly(largeDataName("StringGraphReadDeleted"));
    } catch(...) {
        // Backward compatibility: older datasets may not have read deletion state persisted.
    }
}

void Assembler::accessStringGraphReadWrite()
{
    stringGraph.arcs.accessExistingReadWrite(largeDataName("StringGraphArcs"));
    stringGraph.outgoing.accessExistingReadWrite(largeDataName("StringGraphOutgoing"));
    try {
        stringGraph.incoming.accessExistingReadWrite(largeDataName("StringGraphIncoming"));
    } catch(...) {
        // Backward compatibility: older datasets may not have incoming adjacency persisted.
    }
    try {
        stringGraph.readDeleted.accessExistingReadWrite(largeDataName("StringGraphReadDeleted"));
    } catch(...) {
        // Backward compatibility: older datasets may not have read deletion state persisted.
    }
}

void Assembler::checkStringGraphIsOpen() const
{
    if (!stringGraph.arcs.isOpen) {
        throw runtime_error("String graph arcs are not accessible.");
    }
    if (!stringGraph.outgoing.isOpen()) {
        throw runtime_error("String graph outgoing adjacency is not accessible.");
    }
}
