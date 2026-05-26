
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
    uint64_t threadCount,
    uint64_t minInterWindowCoverage)
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

    // Run CIGAR-based SNP detection on each window.
    {
        cout << timestamp << "Running CIGAR-based SNP detection on "
             << anchorWindows.size() << " windows..." << endl;
        const auto t0 = steady_clock::now();
        uint64_t hetWindows = 0;
        uint64_t totalSnps = 0;
        for(AnchorWindow& window : anchorWindows) {
            window.cleanHetSnpCount = cigarDetectSnpsInWindow(
                window, *shasta2Anchors, *shasta2Journeys);
            if(window.cleanHetSnpCount > 0) {
                hetWindows++;
                totalSnps += window.cleanHetSnpCount;
            }
        }
        const auto t1 = steady_clock::now();
        const double secs = seconds(t1 - t0);
        cout << timestamp << "CIGAR-based SNP detection complete."
             << " hetWindows=" << hetWindows
             << " homWindows=" << (anchorWindows.size() - hetWindows)
             << " totalCleanHetSnps=" << totalSnps
             << " seconds=" << fixed << setprecision(2) << secs
             << defaultfloat << endl;
    }

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
    // For each read, when it crosses from window A to window B, record
    // the (lastAnchorInA, firstAnchorInB) pair. Collect all such pairs
    // per window pair and emit edges for all of them.
    struct AnchorPairKey {
        Shasta2AnchorId anchorIdA;
        Shasta2AnchorId anchorIdB;
        bool operator<(const AnchorPairKey& o) const {
            if(anchorIdA != o.anchorIdA) return anchorIdA < o.anchorIdA;
            return anchorIdB < o.anchorIdB;
        }
    };
    // For each window pair, collect all candidate (anchorA, anchorB) pairs
    // with the number of reads that suggested each pair.
    std::map<std::pair<uint32_t, uint32_t>,
             std::map<AnchorPairKey, uint32_t>> windowPairCandidates;

    const uint64_t orientedReadCount = 2 * readCount;
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
        if(oid.getValue() >= shasta2Journeys->size()) continue;
        const auto journey = (*shasta2Journeys)[oid];
        if(journey.empty()) continue;

        uint32_t currentWindow = noWindow;
        Shasta2AnchorId lastAnchorInCurrentWindow = 0;

        for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
            const Shasta2AnchorId anchorId = journey[pos];
            if(uint64_t(anchorId) >= anchorCount) continue;
            const uint32_t windowId = anchorToWindow[uint64_t(anchorId)];
            if(windowId == noWindow) continue;

            if(windowId == currentWindow) {
                lastAnchorInCurrentWindow = anchorId;
            } else {
                if(currentWindow != noWindow) {
                    auto key = make_pair(currentWindow, windowId);
                    AnchorPairKey apk{lastAnchorInCurrentWindow, anchorId};
                    windowPairCandidates[key][apk]++;
                }
                currentWindow = windowId;
                lastAnchorInCurrentWindow = anchorId;
            }
        }
    }

    // Count total candidate edges across all window pairs.
    uint64_t totalCandidateEdges = 0;
    for(const auto& [windowPair, candidates] : windowPairCandidates) {
        totalCandidateEdges += candidates.size();
    }
    cout << timestamp << "Found " << windowPairCandidates.size()
         << " window pairs with " << totalCandidateEdges
         << " candidate inter-window edges." << endl;

    // Count common oriented reads between two anchors using two-pointer merge.
    // Both anchors' marker info spans are sorted by orientedReadId.
    auto commonReadCount = [&](Shasta2AnchorId anchorIdA, Shasta2AnchorId anchorIdB) -> uint64_t {
        const auto a = (*shasta2Anchors)[anchorIdA];
        const auto b = (*shasta2Anchors)[anchorIdB];
        uint64_t count = 0;
        uint64_t i = 0, j = 0;
        while(i < a.size() && j < b.size()) {
            const auto oidA = a[i].orientedReadId;
            const auto oidB = b[j].orientedReadId;
            if(oidA == oidB) {
                ++count;
                ++i;
                ++j;
            } else if(oidA < oidB) {
                ++i;
            } else {
                ++j;
            }
        }
        return count;
    };

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
    unordered_set<uint64_t> emittedVertices;

    // Write window chains and alternate paths.
    for(const AnchorWindow& window : anchorWindows) {
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];

        // Backbone vertices.
        for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
            if(emittedVertices.insert(uint64_t(backboneJourney[pos])).second) {
                gfa << "S\t" << backboneJourney[pos] << "\t*\tLN:i:1\n";
                ++totalVertices;
            }
        }

        // Intra-window edges: consecutive backbone pairs.
        for(uint32_t pos = window.backboneBegin; pos + 1 < window.backboneEnd; pos++) {
            const Shasta2AnchorId idA = backboneJourney[pos];
            const Shasta2AnchorId idB = backboneJourney[pos + 1];
            gfa << "L\t" << idA << "\t+\t"
                << idB << "\t+\t0M"
                << "\tRC:i:" << commonReadCount(idA, idB) << "\n";
            ++totalIntraEdges;
        }

        // Alternate path vertices and edges — only for het windows.
        if(window.cleanHetSnpCount > 0) {
            for(const AnchorWindowAlternatePath& altPath : window.alternatePaths) {
                for(const Shasta2AnchorId mid : altPath.intermediateAnchorIds) {
                    if(emittedVertices.insert(uint64_t(mid)).second) {
                        gfa << "S\t" << mid << "\t*\tLN:i:1\n";
                        ++totalVertices;
                    }
                }
                // Chain: anchorIdA -> intermediates -> anchorIdB.
                Shasta2AnchorId prev = altPath.anchorIdA;
                for(const Shasta2AnchorId mid : altPath.intermediateAnchorIds) {
                    gfa << "L\t" << prev << "\t+\t" << mid << "\t+\t0M"
                        << "\tRC:i:" << commonReadCount(prev, mid) << "\n";
                    ++totalAltEdges;
                    prev = mid;
                }
                gfa << "L\t" << prev << "\t+\t" << altPath.anchorIdB << "\t+\t0M"
                    << "\tRC:i:" << commonReadCount(prev, altPath.anchorIdB) << "\n";
                ++totalAltEdges;
            }
        }
    }

    // Write inter-window connecting edges: pick the candidate with the
    // highest commonReadCount, skip if below minInterWindowCoverage.
    uint64_t totalInterEdges = 0;
    for(const auto& [windowPair, candidates] : windowPairCandidates) {
        // Find the candidate with the highest RC.
        const AnchorPairKey* bestApk = nullptr;
        uint64_t bestRc = 0;
        for(const auto& [apk, count] : candidates) {
            const uint64_t rc = commonReadCount(apk.anchorIdA, apk.anchorIdB);
            if(rc > bestRc) {
                bestRc = rc;
                bestApk = &apk;
            }
        }
        if(bestApk && bestRc >= minInterWindowCoverage) {
            if(emittedVertices.insert(uint64_t(bestApk->anchorIdA)).second) {
                gfa << "S\t" << bestApk->anchorIdA << "\t*\tLN:i:1\n";
                ++totalVertices;
            }
            if(emittedVertices.insert(uint64_t(bestApk->anchorIdB)).second) {
                gfa << "S\t" << bestApk->anchorIdB << "\t*\tLN:i:1\n";
                ++totalVertices;
            }
            gfa << "L\t" << bestApk->anchorIdA << "\t+\t"
                << bestApk->anchorIdB << "\t+\t0M"
                << "\tRC:i:" << bestRc << "\n";
            ++totalInterEdges;
        }
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

        // Write each backbone anchor and its RC mirror with this color.
        for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
            csv << backboneJourney[pos] << ","
                << "#" << hex << setfill('0')
                << setw(2) << r << setw(2) << g << setw(2) << b
                << dec << "\n";
            const uint64_t rcId = uint64_t(backboneJourney[pos]) ^ 1ULL;
            if(rcId < anchorCount) {
                csv << rcId << ","
                    << "#" << hex << setfill('0')
                    << setw(2) << r << setw(2) << g << setw(2) << b
                    << dec << "\n";
            }
        }
    }

    cout << timestamp << "Wrote " << csvFileName << endl;

    const auto t1 = steady_clock::now();
    const double elapsedSeconds = seconds(t1 - t0);
    cout << timestamp << "testAnchorWindowsCleanLongestRead ends."
         << " seconds=" << fixed << setprecision(2) << elapsedSeconds
         << defaultfloat << endl;
}


// Write AnchorWindowsClean.gfa and .csv from pre-computed windows.
void Assembler::writeAnchorWindowsCleanGfa(
    const vector<AnchorWindow>& anchorWindows,
    uint64_t minInterWindowCoverage)
{
    DINARA_ASSERT(shasta2Anchors);
    DINARA_ASSERT(shasta2Journeys);

    const uint64_t readCount = reads->readCount();
    const uint64_t anchorCount = shasta2Anchors->size();

    // Build anchorId -> windowId and anchorId -> backbone position maps.
    // For each original window W (windowId), also create a mirror RC window
    // (windowId + windowCount) whose backbone anchors are the RC of the originals.
    const uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    const uint32_t windowCount = uint32_t(anchorWindows.size());
    vector<uint32_t> anchorToWindow(anchorCount, noWindow);
    vector<uint32_t> anchorToBackbonePos(anchorCount, 0);
    for(uint32_t windowId = 0; windowId < windowCount; windowId++) {
        const AnchorWindow& window = anchorWindows[windowId];
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];

        // Use filtered backbone positions if available.
        const auto& positions = window.filteredBackbonePositions;
        if(!positions.empty()) {
            for(const uint32_t pos : positions) {
                const Shasta2AnchorId anchorId = backboneJourney[pos];
                anchorToWindow[uint64_t(anchorId)] = windowId;
                anchorToBackbonePos[uint64_t(anchorId)] = pos;
                const uint64_t rcAid = uint64_t(anchorId) ^ 1ULL;
                if(rcAid < anchorCount) {
                    anchorToWindow[rcAid] = windowId + windowCount;
                    anchorToBackbonePos[rcAid] = pos;
                }
            }
        } else {
            for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                const Shasta2AnchorId anchorId = backboneJourney[pos];
                anchorToWindow[uint64_t(anchorId)] = windowId;
                anchorToBackbonePos[uint64_t(anchorId)] = pos;
                const uint64_t rcAid = uint64_t(anchorId) ^ 1ULL;
                if(rcAid < anchorCount) {
                    anchorToWindow[rcAid] = windowId + windowCount;
                    anchorToBackbonePos[rcAid] = pos;
                }
            }
        }
    }

    // Find inter-window connecting edges by walking read journeys.
    // Collect all candidate (anchorA, anchorB) pairs per window pair.
    struct AnchorPairKey {
        Shasta2AnchorId anchorIdA;
        Shasta2AnchorId anchorIdB;
        bool operator<(const AnchorPairKey& o) const {
            if(anchorIdA != o.anchorIdA) return anchorIdA < o.anchorIdA;
            return anchorIdB < o.anchorIdB;
        }
    };
    std::map<std::pair<uint32_t, uint32_t>,
             std::map<AnchorPairKey, uint32_t>> windowPairCandidates;

    const uint64_t orientedReadCount = 2 * readCount;
    for(uint64_t oidValue = 0; oidValue < orientedReadCount; oidValue++) {
        const OrientedReadId oid = OrientedReadId::fromValue(ReadId(oidValue));
        if(oid.getValue() >= shasta2Journeys->size()) continue;
        const auto journey = (*shasta2Journeys)[oid];
        if(journey.empty()) continue;

        uint32_t currentWindow = noWindow;
        Shasta2AnchorId lastAnchorInCurrentWindow = 0;

        for(uint32_t pos = 0; pos < uint32_t(journey.size()); pos++) {
            const Shasta2AnchorId anchorId = journey[pos];
            if(uint64_t(anchorId) >= anchorCount) continue;
            const uint32_t windowId = anchorToWindow[uint64_t(anchorId)];
            if(windowId == noWindow) continue;

            if(windowId == currentWindow) {
                lastAnchorInCurrentWindow = anchorId;
            } else {
                if(currentWindow != noWindow) {
                    auto key = make_pair(currentWindow, windowId);
                    AnchorPairKey apk{lastAnchorInCurrentWindow, anchorId};
                    windowPairCandidates[key][apk]++;
                }
                currentWindow = windowId;
                lastAnchorInCurrentWindow = anchorId;
            }
        }
    }

    uint64_t totalCandidateEdges = 0;
    for(const auto& [windowPair, candidates] : windowPairCandidates) {
        totalCandidateEdges += candidates.size();
    }
    cout << timestamp << "Found " << windowPairCandidates.size()
         << " window pairs with " << totalCandidateEdges
         << " candidate inter-window edges." << endl;

    // Count common oriented reads between two anchors.
    auto commonReadCount = [&](Shasta2AnchorId anchorIdA, Shasta2AnchorId anchorIdB) -> uint64_t {
        const auto a = (*shasta2Anchors)[anchorIdA];
        const auto b = (*shasta2Anchors)[anchorIdB];
        uint64_t count = 0;
        uint64_t i = 0, j = 0;
        while(i < a.size() && j < b.size()) {
            const auto oidA = a[i].orientedReadId;
            const auto oidB = b[j].orientedReadId;
            if(oidA == oidB) { ++count; ++i; ++j; }
            else if(oidA < oidB) { ++i; }
            else { ++j; }
        }
        return count;
    };

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
    unordered_set<uint64_t> emittedVertices;

    // Helper to emit a vertex if not already emitted.
    auto emitVertex = [&](Shasta2AnchorId anchorId) {
        if(emittedVertices.insert(uint64_t(anchorId)).second) {
            gfa << "S\t" << anchorId << "\t*\tLN:i:1\n";
            ++totalVertices;
        }
    };

    // Helper to emit a link.
    auto emitLink = [&](Shasta2AnchorId idA, Shasta2AnchorId idB, uint64_t& counter) {
        gfa << "L\t" << idA << "\t+\t"
            << idB << "\t+\t0M"
            << "\tRC:i:" << commonReadCount(idA, idB) << "\n";
        ++counter;
    };

    for(const AnchorWindow& window : anchorWindows) {
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];

        // Use filtered backbone positions if available.
        const auto& positions = window.filteredBackbonePositions;
        const bool hasFiltered = !positions.empty();

        // Backbone vertices and their RC mirrors.
        if(hasFiltered) {
            for(const uint32_t pos : positions) {
                emitVertex(backboneJourney[pos]);
                const Shasta2AnchorId rcId = Shasta2AnchorId(uint64_t(backboneJourney[pos]) ^ 1ULL);
                if(uint64_t(rcId) < anchorCount) emitVertex(rcId);
            }
        } else {
            for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                emitVertex(backboneJourney[pos]);
                const Shasta2AnchorId rcId = Shasta2AnchorId(uint64_t(backboneJourney[pos]) ^ 1ULL);
                if(uint64_t(rcId) < anchorCount) emitVertex(rcId);
            }
        }

        // Intra-window edges: consecutive filtered backbone pairs + RC mirrors.
        if(hasFiltered) {
            for(uint64_t i = 0; i + 1 < positions.size(); i++) {
                const Shasta2AnchorId idA = backboneJourney[positions[i]];
                const Shasta2AnchorId idB = backboneJourney[positions[i + 1]];
                emitLink(idA, idB, totalIntraEdges);
                const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(idA) ^ 1ULL);
                const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(idB) ^ 1ULL);
                if(uint64_t(rcA) < anchorCount && uint64_t(rcB) < anchorCount) {
                    emitLink(rcB, rcA, totalIntraEdges);
                }
            }
        } else {
            for(uint32_t pos = window.backboneBegin; pos + 1 < window.backboneEnd; pos++) {
                const Shasta2AnchorId idA = backboneJourney[pos];
                const Shasta2AnchorId idB = backboneJourney[pos + 1];
                emitLink(idA, idB, totalIntraEdges);
                const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(idA) ^ 1ULL);
                const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(idB) ^ 1ULL);
                if(uint64_t(rcA) < anchorCount && uint64_t(rcB) < anchorCount) {
                    emitLink(rcB, rcA, totalIntraEdges);
                }
            }
        }

        // Alternate path vertices and edges — only for het windows.
        // Also create RC mirror edges.
        if(window.cleanHetSnpCount > 0) {
            for(const AnchorWindowAlternatePath& altPath : window.alternatePaths) {
                // Forward chain: collect vertices and edges.
                vector<Shasta2AnchorId> forwardChain;
                forwardChain.push_back(altPath.anchorIdA);
                for(const Shasta2AnchorId mid : altPath.intermediateAnchorIds) {
                    emitVertex(mid);
                    const Shasta2AnchorId rcMid = Shasta2AnchorId(uint64_t(mid) ^ 1ULL);
                    if(uint64_t(rcMid) < anchorCount) emitVertex(rcMid);
                    forwardChain.push_back(mid);
                }
                forwardChain.push_back(altPath.anchorIdB);

                // Forward edges.
                for(uint64_t i = 0; i + 1 < forwardChain.size(); i++) {
                    emitLink(forwardChain[i], forwardChain[i+1], totalAltEdges);
                }
                // RC mirror edges (reversed).
                for(uint64_t i = forwardChain.size() - 1; i > 0; i--) {
                    const Shasta2AnchorId rcA = Shasta2AnchorId(uint64_t(forwardChain[i]) ^ 1ULL);
                    const Shasta2AnchorId rcB = Shasta2AnchorId(uint64_t(forwardChain[i-1]) ^ 1ULL);
                    if(uint64_t(rcA) < anchorCount && uint64_t(rcB) < anchorCount) {
                        emitLink(rcA, rcB, totalAltEdges);
                    }
                }
            }
        }
    }

    // Inter-window connecting edges: pick the candidate with the
    // highest commonReadCount, skip if below minInterWindowCoverage.
    uint64_t totalInterEdges = 0;
    for(const auto& [windowPair, candidates] : windowPairCandidates) {
        const AnchorPairKey* bestApk = nullptr;
        uint64_t bestRc = 0;
        for(const auto& [apk, count] : candidates) {
            const uint64_t rc = commonReadCount(apk.anchorIdA, apk.anchorIdB);
            if(rc > bestRc) {
                bestRc = rc;
                bestApk = &apk;
            }
        }
        if(bestApk && bestRc >= minInterWindowCoverage) {
            if(emittedVertices.insert(uint64_t(bestApk->anchorIdA)).second) {
                gfa << "S\t" << bestApk->anchorIdA << "\t*\tLN:i:1\n";
                ++totalVertices;
            }
            if(emittedVertices.insert(uint64_t(bestApk->anchorIdB)).second) {
                gfa << "S\t" << bestApk->anchorIdB << "\t*\tLN:i:1\n";
                ++totalVertices;
            }
            gfa << "L\t" << bestApk->anchorIdA << "\t+\t"
                << bestApk->anchorIdB << "\t+\t0M"
                << "\tRC:i:" << bestRc << "\n";
            ++totalInterEdges;
        }
    }

    cout << timestamp << "Wrote " << gfaFileName
         << ": " << anchorWindows.size() << " chains, "
         << totalVertices << " vertices, "
         << totalIntraEdges << " intra-window edges, "
         << totalAltEdges << " alternate-path edges, "
         << totalInterEdges << " inter-window edges." << endl;

    // Write CSV for Bandage coloring.
    const string csvFileName = "AnchorWindowsClean.csv";
    ofstream csv(csvFileName);
    if(!csv) {
        throw runtime_error("Cannot open " + csvFileName + " for writing.");
    }
    csv << "Name,Color\n";

    for(uint32_t windowId = 0; windowId < windowCount; windowId++) {
        const AnchorWindow& window = anchorWindows[windowId];
        const OrientedReadId backboneOid = window.backboneOrientedReadId;
        const auto backboneJourney = (*shasta2Journeys)[backboneOid];

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

        // Use filtered backbone positions if available.
        const auto& positions = window.filteredBackbonePositions;
        auto emitCsvEntry = [&](uint32_t pos) {
            csv << backboneJourney[pos] << ","
                << "#" << hex << setfill('0')
                << setw(2) << r << setw(2) << g << setw(2) << b
                << dec << "\n";
            const uint64_t rcId = uint64_t(backboneJourney[pos]) ^ 1ULL;
            if(rcId < anchorCount) {
                csv << rcId << ","
                    << "#" << hex << setfill('0')
                    << setw(2) << r << setw(2) << g << setw(2) << b
                    << dec << "\n";
            }
        };

        if(!positions.empty()) {
            for(const uint32_t pos : positions) {
                emitCsvEntry(pos);
            }
        } else {
            for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                emitCsvEntry(pos);
            }
        }
    }

    cout << timestamp << "Wrote " << csvFileName << endl;
}
