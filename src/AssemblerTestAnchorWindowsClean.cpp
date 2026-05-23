
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

    // Build anchorId -> windowId and anchorId -> backbone position maps.
    // Only backbone anchors are "kept" (they define the window chain).
    const uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    vector<uint32_t> anchorToWindow(anchorCount, noWindow);
    vector<uint32_t> anchorToBackbonePos(anchorCount, 0);
    for(uint32_t windowId = 0; windowId < uint32_t(anchorWindows.size()); windowId++) {
        const AnchorWindow& window = anchorWindows[windowId];
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];
        for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
            const Shasta2AnchorId anchorId = backboneJourney[pos];
            anchorToWindow[uint64_t(anchorId)] = windowId;
            anchorToBackbonePos[uint64_t(anchorId)] = pos;
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
        // current window. Only move forward in backbone order within a window.
        uint32_t currentWindow = noWindow;
        Shasta2AnchorId lastAnchorInCurrentWindow = 0;
        uint32_t lastBackbonePosInCurrentWindow = 0;

        for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
            const Shasta2AnchorId anchorId = journey[pos];
            if(uint64_t(anchorId) >= anchorCount) continue;
            const uint32_t windowId = anchorToWindow[uint64_t(anchorId)];
            if(windowId == noWindow) continue;

            const uint32_t backbonePos = anchorToBackbonePos[uint64_t(anchorId)];

            if(windowId == currentWindow) {
                // Same window — only update if moving forward in backbone order.
                if(backbonePos > lastBackbonePosInCurrentWindow) {
                    lastAnchorInCurrentWindow = anchorId;
                    lastBackbonePosInCurrentWindow = backbonePos;
                }
            } else {
                // Window transition.
                if(currentWindow != noWindow) {
                    ConnectingEdge edge{lastAnchorInCurrentWindow, anchorId};
                    ++connectingEdgeCoverage[edge];
                }
                currentWindow = windowId;
                lastAnchorInCurrentWindow = anchorId;
                lastBackbonePosInCurrentWindow = backbonePos;
            }
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

    uint64_t totalAltEdges = 0;

    // Write window chains and alternate paths.
    for(const AnchorWindow& window : anchorWindows) {
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];

        // Backbone vertices.
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

        // Alternate path vertices and edges.
        for(const AnchorWindowAlternatePath& altPath : window.alternatePaths) {
            for(const Shasta2AnchorId mid : altPath.intermediateAnchorIds) {
                gfa << "S\t" << mid << "\t*\tLN:i:1\n";
                ++totalVertices;
            }
            // Chain: anchorIdA -> intermediates -> anchorIdB.
            Shasta2AnchorId prev = altPath.anchorIdA;
            for(const Shasta2AnchorId mid : altPath.intermediateAnchorIds) {
                gfa << "L\t" << prev << "\t+\t" << mid << "\t+\t0M\n";
                ++totalAltEdges;
                prev = mid;
            }
            gfa << "L\t" << prev << "\t+\t" << altPath.anchorIdB << "\t+\t0M\n";
            ++totalAltEdges;
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
         << totalAltEdges << " alternate-path edges, "
         << totalInterEdges << " inter-window edges." << endl;

    // Write CSV for Bandage coloring: each anchor gets a color based on its window.
    const string csvFileName = "AnchorWindowsClean.csv";
    ofstream csv(csvFileName);
    if(!csv) {
        throw runtime_error("Cannot open " + csvFileName + " for writing.");
    }
    csv << "Name,Color\n";

    // Generate distinct colors for each window using HSL with fixed S and L.
    const uint32_t windowCount = uint32_t(anchorWindows.size());
    for(uint32_t windowId = 0; windowId < windowCount; windowId++) {
        const AnchorWindow& window = anchorWindows[windowId];
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];

        // HSL to RGB conversion with S=0.7, L=0.5.
        const double hue = (360.0 * windowId) / windowCount;
        const double s = 0.7;
        const double l = 0.5;
        const double c = (1.0 - std::abs(2.0 * l - 1.0)) * s;
        const double x = c * (1.0 - std::abs(std::fmod(hue / 60.0, 2.0) - 1.0));
        const double m = l - c / 2.0;
        double r1, g1, b1;
        if(hue < 60)       { r1 = c; g1 = x; b1 = 0; }
        else if(hue < 120) { r1 = x; g1 = c; b1 = 0; }
        else if(hue < 180) { r1 = 0; g1 = c; b1 = x; }
        else if(hue < 240) { r1 = 0; g1 = x; b1 = c; }
        else if(hue < 300) { r1 = x; g1 = 0; b1 = c; }
        else               { r1 = c; g1 = 0; b1 = x; }
        const int r = int((r1 + m) * 255);
        const int g = int((g1 + m) * 255);
        const int b = int((b1 + m) * 255);

        // Write each backbone anchor with this color.
        for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
            csv << backboneJourney[pos] << ","
                << "#" << hex << setfill('0')
                << setw(2) << r << setw(2) << g << setw(2) << b
                << dec << "\n";
        }
    }

    cout << timestamp << "Wrote " << csvFileName << endl;

    const auto t1 = steady_clock::now();
    const double elapsedSeconds = seconds(t1 - t0);
    cout << timestamp << "testAnchorWindowsCleanLongestRead ends."
         << " seconds=" << fixed << setprecision(2) << elapsedSeconds
         << defaultfloat << endl;
}
