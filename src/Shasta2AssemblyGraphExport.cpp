// Export dinara's Shasta2AssemblyGraph as a shasta2-compatible AnchorGraph.
//
// Each assembly graph edge step becomes one anchor graph edge.
// This allows shasta2 to consume the cleaned assembly graph
// (after tip removal and superbubble popping) as an external anchor graph.

#include "Shasta2AssemblyGraph.hpp"
#include "Shasta2Anchors.hpp"

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

    const size_t paddingBytes = fileSize - kHeaderSize - data.size();
    if(paddingBytes > 0) {
        vector<char> zeros(paddingBytes, 0);
        out.write(zeros.data(), zeros.size());
    }
}


// Wrapper that mirrors shasta2::AnchorGraph's serialization structure.
struct AnchorGraphWrapper : public shasta2::AnchorGraphBaseClass {
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int) {
        ar & boost::serialization::base_object<shasta2::AnchorGraphBaseClass>(*this);
    }
};

} // anonymous namespace


void Shasta2AssemblyGraph::saveForShasta2(const string& fileName) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;

    // Get anchor count for vertex creation.
    DINARA_ASSERT(anchorsPointer != 0);
    const uint64_t nVertices = anchorsPointer->size();

    // Build a wrapper graph with shasta2 types.
    AnchorGraphWrapper shastaGraph;
    for(uint64_t i = 0; i < nVertices; i++) {
        boost::add_vertex(shastaGraph);
    }

    // Each assembly graph edge step becomes one anchor graph edge.
    uint64_t edgeCount = 0;
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];

        for(const Shasta2AssemblyGraphEdgeStep& step : edge) {
            const Shasta2AnchorPair& dPair = step.anchorPair;

            // Convert dinara AnchorPair -> shasta2 AnchorPair.
            shasta2::AnchorPair shastaPair;
            shastaPair.anchorIdA = dPair.anchorIdA;
            shastaPair.anchorIdB = dPair.anchorIdB;
            shastaPair.orientedReadIds.resize(dPair.orientedReadIds.size());
            for(size_t i = 0; i < dPair.orientedReadIds.size(); i++) {
                shastaPair.orientedReadIds[i] =
                    shasta2::OrientedReadId::fromValue(
                        dPair.orientedReadIds[i].getValue());
            }

            shasta2::AnchorGraphEdge shastaEdge(shastaPair, step.offset, edgeCount);
            shastaEdge.useForAssembly = true;

            boost::add_edge(
                uint64_t(dPair.anchorIdA),
                uint64_t(dPair.anchorIdB),
                shastaEdge,
                shastaGraph);
            ++edgeCount;
        }
    }

    cout << "Exporting AssemblyGraph as AnchorGraph for shasta2: "
         << nVertices << " vertices, "
         << edgeCount << " edges." << endl;

    // Serialize.
    ostringstream oss;
    {
        boost::archive::binary_oarchive archive(oss);
        archive << shastaGraph;
    }

    writeMmapFile(fileName, oss.str());

    cout << "Wrote shasta2-compatible AnchorGraph to " << fileName
         << " (" << oss.str().size() << " bytes serialized)." << endl;
}
