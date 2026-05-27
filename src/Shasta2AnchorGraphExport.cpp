// Export dinara's Shasta2AnchorGraph in shasta2-compatible binary format.
//
// shasta2's --external-anchor-graph-name expects a file in the
// MemoryMapped::Vector<char> format, containing a
// boost::archive::binary_oarchive of shasta2::AnchorGraph.
//
// We build a shasta2::AnchorGraphBaseClass (the boost::adjacency_list
// with shasta2::AnchorGraphEdge properties), serialize it with the
// same archive structure that shasta2::AnchorGraph uses, and write
// the result in the MemoryMapped file format.

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

} // anonymous namespace


void Shasta2AnchorGraph::saveForShasta2(const string& fileName) const
{
    const auto& dinaraGraph = *this;
    const uint64_t nVertices = boost::num_vertices(dinaraGraph);

    // Build a shasta2::AnchorGraphBaseClass and populate it.
    shasta2::AnchorGraphBaseClass shastaGraph(nVertices);

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

    // Serialize using the same pattern as shasta2::AnchorGraph::save(ostream&).
    // shasta2::AnchorGraph::save does: archive << *this
    // shasta2::AnchorGraph::serialize does: ar & base_object<AnchorGraphBaseClass>(*this)
    //
    // For non-polymorphic, non-exported types, boost binary archives
    // don't encode class names. The archive is purely structural:
    // the adjacency_list serialization with AnchorGraphEdge properties.
    // Serializing the base class directly produces identical bytes.
    ostringstream oss;
    {
        boost::archive::binary_oarchive archive(oss);
        archive << shastaGraph;
    }

    writeMmapFile(fileName, oss.str());

    cout << "Wrote shasta2-compatible AnchorGraph to " << fileName
         << " (" << oss.str().size() << " bytes serialized)." << endl;
}
