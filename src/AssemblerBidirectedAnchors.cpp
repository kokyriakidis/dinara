// Assembler method for BRG-native anchor creation.

#include "Assembler.hpp"
#include "mode3-BidirectedAnchor.hpp"
#include "timestamp.hpp"

#include <filesystem>
#include <stdexcept>

using namespace dinara;


shared_ptr<mode3::BidirectedAnchors> Assembler::createBidirectedAnchors(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t minEdgeCoverage,
    uint64_t threadCount)
{
    cout << timestamp << "createBidirectedAnchors begins." << endl;

    auto bidirectedAnchors = make_shared<mode3::BidirectedAnchors>(
        *this,
        getReads(),
        assemblerInfo->k,
        *markers,
        markerGraph,
        minAnchorCoverage,
        maxAnchorCoverage,
        threadCount);

    bidirectedAnchors->computeJourneys(threadCount);
    bidirectedAnchors->computeEdges(threadCount, minEdgeCoverage);

    cout << timestamp << "createBidirectedAnchors ends with "
         << bidirectedAnchors->size() << " anchors." << endl;

    this->bidirectedAnchors = bidirectedAnchors;
    return bidirectedAnchors;
}


void Assembler::saveAnchorGraph(const std::string& fileName, bool includePaths) const
{
    if (!bidirectedAnchors) {
        throw std::runtime_error("Assembler::saveAnchorGraph called but bidirectedAnchors is not set.");
    }
    bidirectedAnchors->writeGfa(fileName, includePaths);
    cout << timestamp << "Saved anchor graph to "
         << std::filesystem::absolute(fileName).string()
         << " (" << bidirectedAnchors->size() << " segments)." << endl;
}



void Assembler::accessBidirectedAnchors()
{
    bidirectedAnchors = make_shared<mode3::BidirectedAnchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers);
}
