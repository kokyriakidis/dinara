
// Dinara.
#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"
using namespace dinara;
using namespace std;

// Standard libraries.
#include "chrono.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <numeric>
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

    // Sort all reads by length (longest first).
    vector<ReadId> readIdsSortedByLength(readCount);
    std::iota(readIdsSortedByLength.begin(), readIdsSortedByLength.end(), ReadId(0));
    std::sort(readIdsSortedByLength.begin(), readIdsSortedByLength.end(),
        [&](ReadId a, ReadId b) {
            return reads->getRead(a).baseCount
                 > reads->getRead(b).baseCount;
        });

    cout << timestamp << "Longest read: " << readIdsSortedByLength[0]
         << " length=" << reads->getRead(readIdsSortedByLength[0]).baseCount
         << " bases" << endl;

    // Run computeAnchorWindowsClean with all reads.
    vector<AnchorWindow> anchorWindows;
    computeAnchorWindowsClean(
        shasta2Anchors,
        shasta2Journeys,
        readIdsSortedByLength,
        anchorWindows,
        threadCount);

    cout << timestamp << "computeAnchorWindowsClean produced "
         << anchorWindows.size() << " windows." << endl;

    // Write GFA: each window is a chain of its backbone anchors.
    const string gfaFileName = "AnchorWindowsClean.gfa";
    ofstream gfa(gfaFileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + gfaFileName + " for writing.");
    }
    gfa << "H\tVN:Z:1.0\n";

    uint64_t totalVertices = 0;
    uint64_t totalEdges = 0;

    for(const AnchorWindow& window : anchorWindows) {
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];

        // Write vertices: backbone anchors in this window.
        for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
            gfa << "S\t" << backboneJourney[pos] << "\t*\tLN:i:1\n";
            ++totalVertices;
        }

        // Write edges: consecutive pairs in the backbone journey.
        for(uint32_t pos = window.backboneBegin; pos + 1 < window.backboneEnd; pos++) {
            gfa << "L\t" << backboneJourney[pos] << "\t+\t"
                << backboneJourney[pos + 1] << "\t+\t0M\n";
            ++totalEdges;
        }
    }

    cout << timestamp << "Wrote " << gfaFileName
         << ": " << anchorWindows.size() << " chains, "
         << totalVertices << " vertices, "
         << totalEdges << " edges." << endl;

    const auto t1 = steady_clock::now();
    const double elapsedSeconds = seconds(t1 - t0);
    cout << timestamp << "testAnchorWindowsCleanLongestRead ends."
         << " seconds=" << fixed << setprecision(2) << elapsedSeconds
         << defaultfloat << endl;
}
