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
// AnchorGraph::serialize() does: ar & base_object<AnchorGraphBaseClass>(*this)
// This adds class tracking metadata that a standalone AnchorGraphBaseClass
// serialization does not include. Without this wrapper, the archive is
// 5 bytes shorter and shasta2's deserializer reads misaligned data.
struct AnchorGraphWrapper : public shasta2::AnchorGraphBaseClass {
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int) {
        ar & boost::serialization::base_object<shasta2::AnchorGraphBaseClass>(*this);
    }
};

} // anonymous namespace


void Shasta2AnchorGraph::saveForShasta2(const string& fileName) const
{
    const auto& dinaraGraph = *this;
    const uint64_t nVertices = boost::num_vertices(dinaraGraph);

    // Build a wrapper graph with shasta2 types.
    AnchorGraphWrapper shastaGraph;
    // Add vertices (AnchorGraphBaseClass uses vecS, so add_vertex is sequential).
    for(uint64_t i = 0; i < nVertices; i++) {
        boost::add_vertex(shastaGraph);
    }

    uint64_t edgeCount = 0;
    BGL_FORALL_EDGES(e, dinaraGraph, Shasta2AnchorGraphBaseClass) {
        const auto src = boost::source(e, dinaraGraph);
        const auto tgt = boost::target(e, dinaraGraph);
        const auto& dEdge = dinaraGraph[e];

        // Convert dinara AnchorPair -> shasta2 AnchorPair.
        shasta2::AnchorPair shastaPair;
        shastaPair.anchorIdA = dEdge.anchorPair.anchorIdA;
        shastaPair.anchorIdB = dEdge.anchorPair.anchorIdB;
        shastaPair.orientedReadIds.resize(dEdge.anchorPair.orientedReadIds.size());
        for(size_t i = 0; i < dEdge.anchorPair.orientedReadIds.size(); i++) {
            shastaPair.orientedReadIds[i] =
                shasta2::OrientedReadId::fromValue(
                    dEdge.anchorPair.orientedReadIds[i].getValue());
        }

        shasta2::AnchorGraphEdge shastaEdge(shastaPair, dEdge.offset, dEdge.id);
        shastaEdge.useForAssembly = dEdge.useForAssembly;

        boost::add_edge(src, tgt, shastaEdge, shastaGraph);
        ++edgeCount;
    }

    cout << "Exporting AnchorGraph for shasta2: "
         << nVertices << " vertices, "
         << edgeCount << " edges." << endl;

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

        // Compare each edge.
        uint64_t edgeIdx = 0;
        uint64_t badEdges = 0;
        auto [origIt, origEnd] = boost::edges(shastaGraph);
        auto [rtIt, rtEnd] = boost::edges(roundTrip);
        for(; origIt != origEnd && rtIt != rtEnd; ++origIt, ++rtIt, ++edgeIdx) {
            const auto& origPair = shastaGraph[*origIt].anchorPair;
            const auto& rtPair = roundTrip[*rtIt].anchorPair;

            bool bad = false;
            if(origPair.anchorIdA != rtPair.anchorIdA ||
               origPair.anchorIdB != rtPair.anchorIdB) {
                bad = true;
            }
            if(origPair.orientedReadIds.size() != rtPair.orientedReadIds.size()) {
                bad = true;
            } else {
                for(size_t i = 0; i < origPair.orientedReadIds.size(); i++) {
                    if(origPair.orientedReadIds[i].getValue() != rtPair.orientedReadIds[i].getValue()) {
                        bad = true;
                        break;
                    }
                }
            }
            if(bad) {
                ++badEdges;
                if(badEdges <= 5) {
                    cout << "ROUND-TRIP EDGE MISMATCH #" << edgeIdx
                         << ": orig " << origPair.anchorIdA << "->" << origPair.anchorIdB
                         << " (" << origPair.orientedReadIds.size() << " reads)"
                         << " vs rt " << rtPair.anchorIdA << "->" << rtPair.anchorIdB
                         << " (" << rtPair.orientedReadIds.size() << " reads)" << endl;
                }
            }
        }
        if(badEdges > 0) {
            cout << "ROUND-TRIP: " << badEdges << " edges differ!" << endl;
            mismatch = true;
        }
        if(!mismatch) {
            cout << "Round-trip verification passed: " << rtVertices << " vertices, "
                 << rtEdges << " edges match." << endl;
        }
    }
}
