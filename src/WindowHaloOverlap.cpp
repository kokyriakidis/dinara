// Read-only diagnostic: measure pairwise window halo overlaps.

#include "WindowHaloOverlap.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace dinara;
using std::cout;
using std::endl;
using std::vector;

void dinara::reportWindowHaloOverlaps(
    const vector<vector<uint32_t>>& windowHalos,
    uint64_t anchorCount)
{
    cout << timestamp << "reportWindowHaloOverlaps begins." << endl;

    const uint64_t windowCount = windowHalos.size();
    if(windowCount == 0) {
        cout << "=== Window halo overlap diagnostic ===\n"
             << "No window halos (halo build disabled or empty).\n"
             << "=== end window halo overlap diagnostic ===" << endl;
        return;
    }

    // Inverted index: anchorId -> windows whose halo contains it.
    // Built compactly so anchors in no halo cost nothing beyond the slot.
    vector<vector<uint32_t>> anchorToWindows(anchorCount);
    uint64_t totalHaloEntries = 0;
    for(uint32_t w = 0; w < windowCount; w++) {
        for(const uint32_t aid : windowHalos[w]) {
            if(aid < anchorCount) {
                anchorToWindows[aid].push_back(w);
                ++totalHaloEntries;
            }
        }
    }

    // Tally pairwise overlaps (shared halo anchors) from co-occurrence.
    // Key: ordered pair (w0 < w1). Value: number of shared halo anchors.
    std::map<std::pair<uint32_t, uint32_t>, uint64_t> pairOverlap;
    uint64_t sharedAnchorCount = 0;     // anchors in >= 2 halos
    uint64_t maxOwnersForAnchor = 0;    // promiscuity of the busiest anchor
    for(uint64_t aid = 0; aid < anchorCount; aid++) {
        const vector<uint32_t>& owners = anchorToWindows[aid];
        if(owners.size() < 2) continue;
        ++sharedAnchorCount;
        maxOwnersForAnchor = std::max<uint64_t>(maxOwnersForAnchor, owners.size());
        for(size_t i = 0; i < owners.size(); i++) {
            for(size_t j = i + 1; j < owners.size(); j++) {
                const uint32_t a = owners[i];
                const uint32_t b = owners[j];
                const auto key = (a < b) ?
                    std::make_pair(a, b) : std::make_pair(b, a);
                ++pairOverlap[key];
            }
        }
    }

    // Per-window: number of distinct other windows it overlaps, and its max
    // single-partner overlap.
    vector<uint32_t> windowPartnerCount(windowCount, 0);
    vector<uint64_t> windowMaxOverlap(windowCount, 0);
    for(const auto& [pair, overlap] : pairOverlap) {
        ++windowPartnerCount[pair.first];
        ++windowPartnerCount[pair.second];
        windowMaxOverlap[pair.first] =
            std::max(windowMaxOverlap[pair.first], overlap);
        windowMaxOverlap[pair.second] =
            std::max(windowMaxOverlap[pair.second], overlap);
    }

    // Overlap-size distribution (bucketed).
    auto bucket = [](uint64_t n) -> const char* {
        if(n == 0) return "0";
        if(n == 1) return "1";
        if(n <= 4) return "2-4";
        if(n <= 16) return "5-16";
        if(n <= 64) return "17-64";
        return "65+";
    };
    std::map<std::string, uint64_t> overlapSizeHist;
    uint64_t maxPairOverlap = 0;
    for(const auto& [pair, overlap] : pairOverlap) {
        (void)pair;
        ++overlapSizeHist[bucket(overlap)];
        maxPairOverlap = std::max(maxPairOverlap, overlap);
    }

    // Partner-count distribution (how many neighbors each window overlaps).
    std::map<std::string, uint64_t> partnerCountHist;
    uint64_t isolatedWindows = 0;
    for(uint32_t w = 0; w < windowCount; w++) {
        const uint32_t pc = windowPartnerCount[w];
        if(pc == 0) ++isolatedWindows;
        ++partnerCountHist[bucket(pc)];
    }

    // Mean halo size.
    const double meanHaloSize = double(totalHaloEntries) / double(windowCount);

    cout << "=== Window halo overlap diagnostic (read-only) ===\n";
    cout << "Windows:                       " << windowCount << "\n";
    cout << "Mean halo size (anchors):      "
         << std::fixed << std::setprecision(1) << meanHaloSize
         << std::defaultfloat << "\n";
    cout << "Anchors shared by >=2 halos:   " << sharedAnchorCount << "\n";
    cout << "Busiest anchor owners:         " << maxOwnersForAnchor << "\n";
    cout << "Overlapping window pairs:      " << pairOverlap.size() << "\n";
    cout << "Largest pair overlap:          " << maxPairOverlap << "\n";
    cout << "-- Overlap size per window pair (shared halo anchors) --\n";
    for(const auto& [b, count] : overlapSizeHist) {
        cout << "  " << b << " anchor(s): " << count << " pair(s)\n";
    }
    cout << "-- Per-window partner count (how many windows it overlaps) --\n";
    cout << "  isolated (0 partners): " << isolatedWindows << " window(s)"
         << "   <- no halo overlap with any neighbor\n";
    for(const auto& [b, count] : partnerCountHist) {
        cout << "  " << b << " partner(s): " << count << " window(s)\n";
    }
    cout << "   A healthy halo model wants most windows overlapping a SMALL\n";
    cout << "   number of neighbors with a SUBSTANTIAL shared-anchor count\n";
    cout << "   (specific overlap). Many partners or tiny overlaps => the\n";
    cout << "   overlap signal is promiscuous/weak.\n";
    cout << "=== end window halo overlap diagnostic ===" << endl;

    cout << timestamp << "reportWindowHaloOverlaps ends." << endl;
}
