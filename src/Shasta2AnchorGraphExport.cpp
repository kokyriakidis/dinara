// Export dinara's Shasta2AnchorGraph in shasta2-compatible binary format.
//
// shasta2's --external-anchor-graph-name expects a file in the
// MemoryMapped::Vector<char> format, containing a
// boost::archive::binary_oarchive of shasta2::AnchorGraph.
//
// We build a shasta2::AnchorGraphBaseClass (the boost::adjacency_list
// with shasta2::AnchorGraphEdge properties), wrap it in a struct that
// mirrors AnchorGraph's serialize() pattern (using base_object<>),
// and write the result in the MemoryMapped file format.

#include "Shasta2AnchorGraph.hpp"

// shasta2 types.
#include "shasta2/AnchorGraph.hpp"
#include "shasta2/AnchorPair.hpp"

// Boost.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

// Standard library.
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace dinara;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ofstream;
using std::ostringstream;


namespace {

// MemoryMapped::Vector<char> file header (must match shasta2's layout exactly).
static constexpr size_t kMagicNumber = 0xa3756fd4b5d8bcc1ULL;
static constexpr size_t kHeaderSize  = 4096;
static constexpr size_t kPageSize    = 4096;

struct MmapHeader {
    size_t headerSize;
    size_t objectSize;
    size_t objectCount;
    size_t pageSize;
    size_t pageCount;
    size_t fileSize;
    size_t capacity;
    size_t magicNumber;
    char   padding[kHeaderSize - 8 * sizeof(size_t)];
};
static_assert(sizeof(MmapHeader) == kHeaderSize);


// Write serialized data in MemoryMapped::Vector<char> file format.
void writeMmapFile(const string& fileName, const string& data)
{
    const size_t objectCount = data.size();
    const size_t totalBytes  = kHeaderSize + objectCount;
    const size_t pageCount   = (totalBytes + kPageSize - 1) / kPageSize;
    const size_t fileSize    = pageCount * kPageSize;

    MmapHeader header;
    std::memset(&header, 0, sizeof(header));
    header.headerSize  = kHeaderSize;
    header.objectSize  = sizeof(char);
    header.objectCount = objectCount;
    header.pageSize    = kPageSize;
    header.pageCount   = pageCount;
    header.fileSize    = fileSize;
    header.capacity    = (fileSize - kHeaderSize) / sizeof(char);
    header.magicNumber = kMagicNumber;

    ofstream out(fileName, std::ios::binary);
    if(!out) {
        throw std::runtime_error("Cannot open " + fileName + " for writing.");
    }
    out.write(reinterpret_cast<const char*>(&header), kHeaderSize);
    out.write(data.data(), data.size());

    // Pad to page-aligned fileSize.
    const size_t paddingBytes = fileSize - kHeaderSize - data.size();
    if(paddingBytes > 0) {
        vector<char> zeros(paddingBytes, 0);
        out.write(zeros.data(), zeros.size());
    }
}


// Wrapper that mirrors shasta2::AnchorGraph's serialization structure.
// Upstream AnchorGraph::serialize() does:
//     ar & base_object<AnchorGraphBaseClass>(*this);
//     ar & orientedReadIds;
// Each AnchorGraphEdge stores anchorIdA/anchorIdB plus [begin,end) indices into
// this graph-level orientedReadIds vector (rather than an inline AnchorPair), so
// the wrapper must carry and serialize that vector for the archive to match what
// shasta2's --external-anchor-graph-name loader expects.
struct AnchorGraphWrapper : public shasta2::AnchorGraphBaseClass {
    vector<shasta2::OrientedReadId> orientedReadIds;
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int) {
        ar & boost::serialization::base_object<shasta2::AnchorGraphBaseClass>(*this);
        ar & orientedReadIds;
    }
};

} // anonymous namespace


void Shasta2AnchorGraph::saveForShasta2(
    const string& fileName, const Shasta2Anchors& anchors) const
{
    const auto& dinaraGraph = *this;
    const uint64_t nVertices = boost::num_vertices(dinaraGraph);

    // Build a wrapper graph with shasta2 types.
    AnchorGraphWrapper shastaGraph;
    // Add vertices (AnchorGraphBaseClass uses vecS, so add_vertex is sequential).
    for(uint64_t i = 0; i < nVertices; i++) {
        boost::add_vertex(shastaGraph);
    }

    // Export-time monotonicity verification, bound to the exact edge set being
    // serialized. The Shasta2AnchorGraph constructor already verifies edges, but
    // that runs BEFORE trimBackbones and on the in-memory graph; this pass runs
    // on the precise set written to disk, so nothing between construction and
    // export can slip a backward/equal-position edge into the file shasta2 loads.
    // For each exported edge A -> B and each read on the pair, the read's stored
    // position at B must strictly exceed its position at A (shasta2 assembles by
    // read-following in increasing position and asserts positionB > positionA).
    uint64_t verifiedEdges = 0;
    uint64_t verifiedReadSteps = 0;

    uint64_t edgeCount = 0;
    uint64_t skippedEdgeCount = 0;
    BGL_FORALL_EDGES(e, dinaraGraph, Shasta2AnchorGraphBaseClass) {
        const auto& dEdge = dinaraGraph[e];

        // Only export active edges.
        if(!dEdge.useForAssembly) {
            ++skippedEdgeCount;
            continue;
        }

        const auto src = boost::source(e, dinaraGraph);
        const auto tgt = boost::target(e, dinaraGraph);

        // Verify this edge is forward-monotonic on every shared read before it
        // is written. Both anchor member lists are sorted by OrientedReadId, so
        // shared reads are found by a linear merge. A read on the edge's anchor
        // pair that is missing from an anchor is skipped here (shasta2's
        // gatherOrientedReads likewise skips it); only reads present on BOTH
        // anchors constrain assembly order.
        {
            const Shasta2Anchor anchorA = anchors[dEdge.anchorPair.anchorIdA];
            const Shasta2Anchor anchorB = anchors[dEdge.anchorPair.anchorIdB];
            auto itA = anchorA.begin();
            auto itB = anchorB.begin();
            const auto endA = anchorA.end();
            const auto endB = anchorB.end();
            while(itA != endA && itB != endB) {
                if(itA->orientedReadId < itB->orientedReadId) { ++itA; continue; }
                if(itB->orientedReadId < itA->orientedReadId) { ++itB; continue; }
                // Shared read: exported ordering must be strictly forward. Every
                // anchor exports position - k/2 with the SAME uniform export
                // shift (1 for k=2, 0 for k=0), so the shift cancels regardless
                // of k; compare stored positions directly.
                if(!(itB->position > itA->position)) {
                    throw runtime_error(
                        "Shasta2 anchor-graph export failed: backward edge "
                        "anchor " +
                        shasta2AnchorIdToString(dEdge.anchorPair.anchorIdA) +
                        " -> anchor " +
                        shasta2AnchorIdToString(dEdge.anchorPair.anchorIdB) +
                        " on oriented read " +
                        std::to_string(itA->orientedReadId.getValue()) +
                        " has non-increasing position (" +
                        std::to_string(itA->position) + " -> " +
                        std::to_string(itB->position) + ").");
                }
                ++verifiedReadSteps;
                ++itA;
                ++itB;
            }
            ++verifiedEdges;
        }

        // Append this edge's oriented reads to the shared vector and record the
        // [begin,end) range (upstream's AnchorGraph::addEdge layout). The edge no
        // longer carries an offset; shasta2 recomputes offsets during assembly.
        const uint64_t begin = shastaGraph.orientedReadIds.size();
        for(const auto& rid : dEdge.anchorPair.orientedReadIds) {
            shastaGraph.orientedReadIds.push_back(
                shasta2::OrientedReadId::fromValue(rid.getValue()));
        }
        const uint64_t end = shastaGraph.orientedReadIds.size();

        shasta2::AnchorGraphEdge shastaEdge(
            dEdge.anchorPair.anchorIdA,
            dEdge.anchorPair.anchorIdB,
            begin, end,
            dEdge.id,
            /* useForAssembly */ true);

        boost::add_edge(src, tgt, shastaEdge, shastaGraph);
        ++edgeCount;
    }

    cout << "Export monotonicity check: verified " << verifiedEdges
         << " edges forward across " << verifiedReadSteps
         << " shared-read steps." << endl;

    cout << "Exporting AnchorGraph for shasta2: "
         << nVertices << " vertices, "
         << edgeCount << " active edges ("
         << skippedEdgeCount << " disabled edges skipped)." << endl;

    // Serialize with the same archive structure as shasta2::AnchorGraph.
    ostringstream oss;
    {
        boost::archive::binary_oarchive archive(oss);
        archive << shastaGraph;
    }

    writeMmapFile(fileName, oss.str());

    cout << "Wrote shasta2-compatible AnchorGraph to " << fileName
         << " (" << oss.str().size() << " bytes serialized)." << endl;

    // Round-trip verification: deserialize and compare.
    {
        const string& serialized = oss.str();
        std::istringstream iss(serialized);
        AnchorGraphWrapper roundTrip;
        {
            boost::archive::binary_iarchive archive(iss);
            archive >> roundTrip;
        }

        const uint64_t rtVertices = boost::num_vertices(roundTrip);
        const uint64_t rtEdges = boost::num_edges(roundTrip);
        bool mismatch = false;

        if(rtVertices != nVertices) {
            cout << "ROUND-TRIP MISMATCH: vertices " << rtVertices << " vs " << nVertices << endl;
            mismatch = true;
        }
        if(rtEdges != edgeCount) {
            cout << "ROUND-TRIP MISMATCH: edges " << rtEdges << " vs " << edgeCount << endl;
            mismatch = true;
        }
        if(roundTrip.orientedReadIds.size() != shastaGraph.orientedReadIds.size()) {
            cout << "ROUND-TRIP MISMATCH: orientedReadIds "
                 << roundTrip.orientedReadIds.size() << " vs "
                 << shastaGraph.orientedReadIds.size() << endl;
            mismatch = true;
        }

        // Compare each edge: anchor ids, index range, and the reads in that range.
        uint64_t edgeIdx = 0;
        uint64_t badEdges = 0;
        auto [origIt, origEnd] = boost::edges(shastaGraph);
        auto [rtIt, rtEnd] = boost::edges(roundTrip);
        for(; origIt != origEnd && rtIt != rtEnd; ++origIt, ++rtIt, ++edgeIdx) {
            const auto& origEdge = shastaGraph[*origIt];
            const auto& rtEdge = roundTrip[*rtIt];

            bool bad = false;
            if(origEdge.anchorIdA != rtEdge.anchorIdA ||
               origEdge.anchorIdB != rtEdge.anchorIdB) {
                bad = true;
            }
            const uint64_t origN = origEdge.orientedReadIdsEnd - origEdge.orientedReadIdsBegin;
            const uint64_t rtN = rtEdge.orientedReadIdsEnd - rtEdge.orientedReadIdsBegin;
            if(origN != rtN) {
                bad = true;
            } else {
                for(uint64_t i = 0; i < origN; i++) {
                    const auto o = shastaGraph.orientedReadIds[origEdge.orientedReadIdsBegin + i];
                    const auto r = roundTrip.orientedReadIds[rtEdge.orientedReadIdsBegin + i];
                    if(o.getValue() != r.getValue()) { bad = true; break; }
                }
            }
            if(bad) {
                ++badEdges;
                if(badEdges <= 5) {
                    cout << "ROUND-TRIP EDGE MISMATCH #" << edgeIdx
                         << ": orig " << origEdge.anchorIdA << "->" << origEdge.anchorIdB
                         << " (" << origN << " reads)"
                         << " vs rt " << rtEdge.anchorIdA << "->" << rtEdge.anchorIdB
                         << " (" << rtN << " reads)" << endl;
                }
            }
        }
        if(badEdges > 0) {
            cout << "ROUND-TRIP: " << badEdges << " edges differ!" << endl;
            mismatch = true;
        }
        if(!mismatch) {
            cout << "Round-trip verification passed: " << rtVertices << " vertices, "
                 << rtEdges << " edges, "
                 << roundTrip.orientedReadIds.size() << " oriented reads match." << endl;
        }
    }
}
