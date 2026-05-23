
// Dinara.
#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2AnchorPair.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"
using namespace dinara;
using namespace std;

// Standard libraries.
#include "chrono.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <unordered_set>
#include <vector>



void Assembler::testAnchorWindowsCleanLongestRead(
    uint64_t threadCount)
{
    cout << timestamp << "testAnchorWindowsCleanLongestRead begins." << endl;
    const auto t0 = steady_clock::now();

    DINARA_ASSERT(shasta2Anchors);
    DINARA_ASSERT(shasta2Journeys);
    DINARA_ASSERT(shasta2Journeys->isOpen());
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();

    const uint64_t readCount = reads->readCount();
    DINARA_ASSERT(readCount > 0);

    // Find the longest read.
    ReadId longestReadId = 0;
    uint64_t longestLength = 0;
    for(ReadId readId = 0; readId < readCount; readId++) {
        const uint64_t len = reads->getRead(readId).baseCount;
        if(len > longestLength) {
            longestLength = len;
            longestReadId = readId;
        }
    }

    cout << timestamp << "Longest read: " << longestReadId
         << " length=" << longestLength << " bases" << endl;

    // Check that the longest read has a journey.
    const OrientedReadId longestOid(longestReadId, 0);
    if(longestOid.getValue() >= shasta2Journeys->size()) {
        cout << timestamp << "Longest read has no journey. Aborting test." << endl;
        return;
    }
    const auto longestJourney = (*shasta2Journeys)[longestOid];
    cout << timestamp << "Longest read journey has " << longestJourney.size()
         << " anchors." << endl;
    if(longestJourney.empty()) {
        cout << timestamp << "Longest read journey is empty. Aborting test." << endl;
        return;
    }

    // Run computeAnchorWindowsClean with only the longest read as backbone candidate.
    vector<ReadId> singleRead = {longestReadId};
    vector<AnchorWindow> anchorWindows;
    computeAnchorWindowsClean(
        shasta2Anchors,
        shasta2Journeys,
        singleRead,
        anchorWindows,
        threadCount);

    cout << timestamp << "computeAnchorWindowsClean produced "
         << anchorWindows.size() << " window(s) for the longest read." << endl;

    // Collect kept anchor IDs from Window 0 only (the longest read's window).
    // The other windows are from touched reads that were re-pushed onto the
    // heap and are not relevant for this test.
    unordered_set<Shasta2AnchorId> keptAnchorSet;
    uint64_t totalReadIntervals = 0;
    if(!anchorWindows.empty()) {
        const AnchorWindow& window = anchorWindows[0];
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];
        for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
            keptAnchorSet.insert(backboneJourney[pos]);
        }
        totalReadIntervals = window.readIntervals.size();
    }

    cout << timestamp << "Kept anchors: " << keptAnchorSet.size()
         << " (out of " << shasta2Anchors->size() << " total)"
         << " totalReadIntervals=" << totalReadIntervals << endl;

    // Print per-window diagnostics.
    for(uint64_t i = 0; i < anchorWindows.size(); i++) {
        const AnchorWindow& w = anchorWindows[i];
        cout << timestamp << "  Window " << i
             << ": backbone=" << w.backboneOrientedReadId
             << " [" << w.backboneBegin << "," << w.backboneEnd << ")"
             << " backboneAnchors=" << (w.backboneEnd - w.backboneBegin)
             << " claimedAnchors=" << w.claimedAnchorCount
             << " reads=" << w.readIntervals.size()
             << endl;
    }

    // Build a restricted anchor graph using only the kept anchors.
    // Instead of constructing the full Shasta2AnchorGraph (expensive),
    // we iterate only over kept anchors and find children among them.
    cout << timestamp << "Building restricted anchor graph from kept anchors..." << endl;

    // For each kept anchor, find children (minEdgeCoverage=1 to include
    // transitions supported by even a single read) and keep only those
    // that are also kept.
    struct RestrictedEdge {
        Shasta2AnchorId anchorIdA;
        Shasta2AnchorId anchorIdB;
        uint64_t coverage;
    };
    vector<RestrictedEdge> restrictedEdges;

    vector<Shasta2AnchorId> children;
    vector<uint64_t> counts;
    for(const Shasta2AnchorId anchorIdA : keptAnchorSet) {
        shasta2Anchors->findChildren(*shasta2Journeys, anchorIdA, children, counts, 1);
        DINARA_ASSERT(children.size() == counts.size());
        for(uint64_t i = 0; i < children.size(); i++) {
            const Shasta2AnchorId anchorIdB = children[i];
            if(keptAnchorSet.find(anchorIdB) == keptAnchorSet.end()) {
                continue;
            }
            restrictedEdges.push_back(RestrictedEdge{
                anchorIdA, anchorIdB, counts[i]});
        }
    }

    cout << timestamp << "Restricted anchor graph: "
         << keptAnchorSet.size() << " vertices, "
         << restrictedEdges.size() << " edges." << endl;

    // Write GFA.
    const string gfaFileName = "LongestReadAnchorGraph.gfa";
    ofstream gfa(gfaFileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + gfaFileName + " for writing.");
    }
    gfa << "H\tVN:Z:1.0\n";

    // Write vertices for kept anchors.
    for(const Shasta2AnchorId anchorId : keptAnchorSet) {
        gfa << "S\t" << anchorId << "\t*\tLN:i:1\n";
    }

    // Write edges.
    for(const RestrictedEdge& edge : restrictedEdges) {
        gfa << "L\t" << edge.anchorIdA << "\t+\t"
            << edge.anchorIdB << "\t+\t0M"
            << "\tRC:i:" << edge.coverage
            << "\n";
    }

    cout << timestamp << "Wrote " << gfaFileName << endl;

    const auto t1 = steady_clock::now();
    const double elapsedSeconds = seconds(t1 - t0);
    cout << timestamp << "testAnchorWindowsCleanLongestRead ends."
         << " seconds=" << fixed << setprecision(2) << elapsedSeconds
         << defaultfloat << endl;
}
