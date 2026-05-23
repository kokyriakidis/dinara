
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
#include <map>
#include <numeric>
#include <set>
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
    const uint64_t anchorCount = shasta2Anchors->size();

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

    // Build anchorId -> windowId map from backbone anchors.
    // Only backbone anchors are "kept" (they define the window chain).
    const uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    vector<uint32_t> anchorToWindow(anchorCount, noWindow);
    for(uint32_t windowId = 0; windowId < uint32_t(anchorWindows.size()); windowId++) {
        const AnchorWindow& window = anchorWindows[windowId];
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];
        for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
            const Shasta2AnchorId anchorId = backboneJourney[pos];
            anchorToWindow[uint64_t(anchorId)] = windowId;
        }
    }

    // Find inter-window connecting edges by walking read journeys.
    // For each read, walk its journey and track which window each anchor
    // belongs to. When the window changes, record a connecting edge from
    // the last backbone anchor in the previous window to the first backbone
    // anchor in the next window.
    // Use a set to deduplicate edges and count coverage.
    struct ConnectingEdge {
        Shasta2AnchorId anchorIdA;  // Last backbone anchor in window A.
        Shasta2AnchorId anchorIdB;  // First backbone anchor in window B.
        bool operator<(const ConnectingEdge& other) const {
            if(anchorIdA != other.anchorIdA) return anchorIdA < other.anchorIdA;
            return anchorIdB < other.anchorIdB;
        }
    };
    map<ConnectingEdge, uint64_t> connectingEdgeCoverage;

    const uint64_t orientedReadCount = 2 * readCount;
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
        if(oid.getValue() >= shasta2Journeys->size()) continue;
        const auto journey = (*shasta2Journeys)[oid];
        if(journey.empty()) continue;

        // Walk the journey, tracking the last backbone anchor seen in the
        // current window.
        uint32_t currentWindow = noWindow;
        Shasta2AnchorId lastAnchorInCurrentWindow = 0;

        for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
            const Shasta2AnchorId anchorId = journey[pos];
            if(uint64_t(anchorId) >= anchorCount) continue;
            const uint32_t windowId = anchorToWindow[uint64_t(anchorId)];
            if(windowId == noWindow) continue;

            if(windowId != currentWindow) {
                // Window transition.
                if(currentWindow != noWindow) {
                    // Add connecting edge from last anchor in previous window
                    // to this anchor (first in new window).
                    ConnectingEdge edge{lastAnchorInCurrentWindow, anchorId};
                    ++connectingEdgeCoverage[edge];
                }
                currentWindow = windowId;
            }
            lastAnchorInCurrentWindow = anchorId;
        }
    }

    cout << timestamp << "Found " << connectingEdgeCoverage.size()
         << " inter-window connecting edges." << endl;

    // Write GFA.
    const string gfaFileName = "AnchorWindowsClean.gfa";
    ofstream gfa(gfaFileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + gfaFileName + " for writing.");
    }
    gfa << "H\tVN:Z:1.0\n";

    uint64_t totalVertices = 0;
    uint64_t totalIntraEdges = 0;

    // Write window chains.
    for(const AnchorWindow& window : anchorWindows) {
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];

        // Vertices.
        for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
            gfa << "S\t" << backboneJourney[pos] << "\t*\tLN:i:1\n";
            ++totalVertices;
        }

        // Intra-window edges: consecutive backbone pairs.
        for(uint32_t pos = window.backboneBegin; pos + 1 < window.backboneEnd; pos++) {
            gfa << "L\t" << backboneJourney[pos] << "\t+\t"
                << backboneJourney[pos + 1] << "\t+\t0M\n";
            ++totalIntraEdges;
        }
    }

    // Write inter-window connecting edges.
    uint64_t totalInterEdges = 0;
    for(const auto& [edge, coverage] : connectingEdgeCoverage) {
        gfa << "L\t" << edge.anchorIdA << "\t+\t"
            << edge.anchorIdB << "\t+\t0M"
            << "\tRC:i:" << coverage << "\n";
        ++totalInterEdges;
    }

    cout << timestamp << "Wrote " << gfaFileName
         << ": " << anchorWindows.size() << " chains, "
         << totalVertices << " vertices, "
         << totalIntraEdges << " intra-window edges, "
         << totalInterEdges << " inter-window edges." << endl;

    const auto t1 = steady_clock::now();
    const double elapsedSeconds = seconds(t1 - t0);
    cout << timestamp << "testAnchorWindowsCleanLongestRead ends."
         << " seconds=" << fixed << setprecision(2) << elapsedSeconds
         << defaultfloat << endl;
}
