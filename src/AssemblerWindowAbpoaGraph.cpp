// AssemblerWindowAbpoaGraph.cpp
//
// Per-window progressive abPOA graph construction.
//
// For each anchor window:
//   1. Seed an abPOA graph with the backbone read's base sequence (the
//      sequence spanning the window's backbone anchors). The backbone is
//      always sequence 0, so it becomes the spine of the POA graph.
//   2. For every other member read in the window, find the backbone anchors
//      it shares with the backbone. The first and last shared anchor define
//      the read's "anchor interval" with the backbone. The read is then
//      aligned only to the backbone subgraph spanning that interval (via
//      abpoa_align_sequence_to_subgraph), and added progressively.
//   3. Write the resulting graph as GFA (one file per window) for diagnostics.
//
// The abPOA graph / MSA produced here is the substrate for later het-site
// detection (a separate pass can read the GFA or re-run with MSA output).
//
// This mirrors the verified progressive-subgraph pattern in
// PhasingKmeansAlign.cpp::abpoaMsaRun and the shared-anchor derivation in
// AssemblerAnchorWindowsClean.cpp::createWindow.

#include "Assembler.hpp"
#include "AnchorWindows.hpp"
#include "Base.hpp"
#include "Marker.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"

#include "abpoa/abpoa.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// The interval a single ORIENTED read (read + strand, never just a ReadId)
// shares with the window backbone. The same physical read in its two
// orientations is two distinct OrientedReadIds and yields two distinct
// intervals here.
//
// The interval is expressed in two coordinate spaces, both derived from the
// same set of shared anchors:
//   - Backbone space (base coords relative to the backbone sequence start):
//     [backboneBegin, backboneEnd). Selects the abPOA graph slice.
//   - Read space (this oriented read's marker ordinals):
//     [readOrdinalBegin, readOrdinalEnd). Selects the bases to extract.
struct OrientedReadInterval {
    OrientedReadId orientedReadId;
    int backboneBegin = 0;       // Inclusive backbone base position of first shared anchor.
    int backboneEnd = 0;         // Exclusive backbone base position past last shared anchor.
    uint32_t readOrdinalBegin = 0; // First shared-anchor marker ordinal on this oriented read.
    uint32_t readOrdinalEnd = 0;   // Last shared-anchor marker ordinal on this oriented read.
    vector<uint8_t> seq;         // Read bases over [readOrdinalBegin, readOrdinalEnd) (0123).
};

inline uint8_t baseToInt(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;
    }
}

} // anonymous namespace

void Assembler::computeWindowAbpoaGraphs(
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const string& outputPrefix,
    uint64_t threadCount) const
{
    if(!markers) {
        cout << timestamp << "computeWindowAbpoaGraphs: markers not available. Skipping." << endl;
        return;
    }
    const Reads& reads = getReads();
    const uint32_t kHalf = assemblerInfo.isOpen ? uint32_t(assemblerInfo->k / 2) : 0U;

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
        if(threadCount == 0) threadCount = 1;
    }

    cout << timestamp << "computeWindowAbpoaGraphs: processing "
         << anchorWindows.size() << " windows with " << threadCount
         << " threads (prefix=\"" << outputPrefix << "\")." << endl;

    std::atomic<uint64_t> nextWindow{0};
    std::atomic<uint64_t> windowsWritten{0};
    std::atomic<uint64_t> windowsSkipped{0};
    std::atomic<uint64_t> readsAligned{0};

    // Extract the base sequence of an oriented read between two marker
    // ordinals (midpoint to midpoint), encoded as 0123. Mirrors
    // PhasingKmeansAlign.cpp::extractSegmentSeq0123.
    const auto extractSeq0123 = [&](OrientedReadId oid, uint32_t ordA, uint32_t ordB) -> vector<uint8_t> {
        vector<uint8_t> seq;
        if(ordA >= ordB) return seq;
        const auto readMarkers = (*markers)[oid.getValue()];
        if(ordB >= readMarkers.size()) return seq;
        const uint32_t beginPos = readMarkers[ordA].position + kHalf;
        const uint32_t endPos   = readMarkers[ordB].position + kHalf;
        if(endPos <= beginPos) return seq;
        seq.reserve(endPos - beginPos);
        for(uint32_t pos = beginPos; pos < endPos; pos++) {
            seq.push_back(baseToInt(reads.getOrientedReadBase(oid, pos).character()));
        }
        return seq;
    };

    auto worker = [&]() {
        // Each thread owns its own abPOA instances — no shared mutable state.
        for(;;) {
            const uint64_t wi = nextWindow.fetch_add(1);
            if(wi >= anchorWindows.size()) break;
            const AnchorWindow& window = anchorWindows[wi];

            const OrientedReadId backboneOid = window.backboneOrientedReadId;

            // Ordered backbone journey positions forming the spine.
            const auto backboneJourney = journeys[backboneOid];
            vector<uint32_t> bbPositions = window.filteredBackbonePositions;
            if(bbPositions.empty()) {
                for(uint32_t pos = window.backboneBegin; pos < window.backboneEnd; pos++) {
                    bbPositions.push_back(pos);
                }
            }
            if(bbPositions.size() < 2) {
                windowsSkipped.fetch_add(1);
                continue;
            }

            // Backbone base coordinates: midpoint of first and last backbone
            // anchor. The backbone sequence spans [firstMid, lastMid).
            const uint32_t firstOrdinal = anchors.getOrdinal(
                backboneJourney[bbPositions.front()], backboneOid);
            const uint32_t lastOrdinal = anchors.getOrdinal(
                backboneJourney[bbPositions.back()], backboneOid);
            const auto bbMarkers = (*markers)[backboneOid.getValue()];
            if(lastOrdinal >= bbMarkers.size() || firstOrdinal >= lastOrdinal) {
                windowsSkipped.fetch_add(1);
                continue;
            }
            const uint32_t bbStartBase = bbMarkers[firstOrdinal].position + kHalf;
            const uint32_t bbEndBase   = bbMarkers[lastOrdinal].position + kHalf;
            if(bbEndBase <= bbStartBase) {
                windowsSkipped.fetch_add(1);
                continue;
            }

            // Backbone sequence (0123).
            vector<uint8_t> backboneSeq;
            backboneSeq.reserve(bbEndBase - bbStartBase);
            for(uint32_t pos = bbStartBase; pos < bbEndBase; pos++) {
                backboneSeq.push_back(baseToInt(reads.getOrientedReadBase(backboneOid, pos).character()));
            }
            const int backboneLen = int(backboneSeq.size());
            if(backboneLen < 2) {
                windowsSkipped.fetch_add(1);
                continue;
            }

            // Map backbone anchorId -> base position relative to backbone start.
            // Only positions that lie within [bbStartBase, bbEndBase] map to a
            // valid backbone node (0-based into backboneSeq).
            std::unordered_map<uint64_t, int> backboneAnchorToBasePos;
            backboneAnchorToBasePos.reserve(bbPositions.size() * 2);
            for(const uint32_t pos : bbPositions) {
                const Shasta2AnchorId aid = backboneJourney[pos];
                const uint32_t ord = anchors.getOrdinal(aid, backboneOid);
                if(ord >= bbMarkers.size()) continue;
                const uint32_t mid = bbMarkers[ord].position + kHalf;
                if(mid < bbStartBase || mid > bbEndBase) continue;
                backboneAnchorToBasePos[uint64_t(aid)] = int(mid - bbStartBase);
            }

            // Compute the oriented read interval for each non-backbone member:
            // find its shared backbone anchors -> backbone-space interval +
            // read-space ordinal interval + extracted read sequence. Keyed by
            // OrientedReadId, so a read's two strands produce two intervals.
            vector<OrientedReadInterval> members;
            members.reserve(window.readIntervals.size());

            for(const AnchorWindowReadInterval& ri : window.readIntervals) {
                const OrientedReadId oid = ri.orientedReadId;
                if(oid == backboneOid) continue;
                if(oid.getValue() >= journeys.size()) continue;
                const auto journey = journeys[oid];
                if(journey.empty()) continue;

                // Shared backbone anchors, recording read ordinal + backbone
                // base position, in this oriented read's journey order.
                int backboneBegin = -1, backboneEnd = -1;
                uint32_t readOrdinalBegin = UINT32_MAX, readOrdinalEnd = 0;
                bool any = false;
                for(uint32_t readPos = 0; readPos < uint32_t(journey.size()); readPos++) {
                    const uint64_t aid = uint64_t(journey[readPos]);
                    auto it = backboneAnchorToBasePos.find(aid);
                    if(it == backboneAnchorToBasePos.end()) continue;
                    const int bbBase = it->second;
                    const uint32_t readOrd = anchors.getOrdinal(Shasta2AnchorId(aid), oid);
                    if(!any) {
                        backboneBegin = bbBase;
                        backboneEnd = bbBase;
                        readOrdinalBegin = readOrd;
                        readOrdinalEnd = readOrd;
                        any = true;
                    } else {
                        backboneBegin = std::min(backboneBegin, bbBase);
                        backboneEnd = std::max(backboneEnd, bbBase);
                        readOrdinalBegin = std::min(readOrdinalBegin, readOrd);
                        readOrdinalEnd = std::max(readOrdinalEnd, readOrd);
                    }
                }
                if(!any || readOrdinalBegin >= readOrdinalEnd) continue;

                OrientedReadInterval interval;
                interval.orientedReadId = oid;
                interval.backboneBegin = backboneBegin;
                // Exclusive right bound; clamp to backbone length.
                interval.backboneEnd = std::min(backboneEnd + 1, backboneLen);
                interval.readOrdinalBegin = readOrdinalBegin;
                interval.readOrdinalEnd = readOrdinalEnd;
                if(interval.backboneBegin >= interval.backboneEnd) continue;
                interval.seq = extractSeq0123(oid, readOrdinalBegin, readOrdinalEnd);
                if(interval.seq.size() < 2) continue;
                members.push_back(std::move(interval));
            }

            // Add reads to the MSA from the biggest shared interval to the
            // smallest. The order reads are added to a POA graph matters:
            // adding the longest-overlapping reads first builds the most
            // reliable spine before shorter reads are layered on. Interval
            // size is the backbone span [backboneBegin, backboneEnd); ties broken
            // by extracted read length, then read id for determinism.
            std::sort(members.begin(), members.end(),
                [](const OrientedReadInterval& a, const OrientedReadInterval& b) {
                    const int spanA = a.backboneEnd - a.backboneBegin;
                    const int spanB = b.backboneEnd - b.backboneBegin;
                    if(spanA != spanB) return spanA > spanB;
                    if(a.seq.size() != b.seq.size()) return a.seq.size() > b.seq.size();
                    return a.orientedReadId < b.orientedReadId;
                });

            // Total sequences = backbone + members.
            const int nMembers = int(members.size());
            const int totalSeqs = 1 + nMembers;

            // --- Run progressive abPOA ---
            abpoa_t* ab = abpoa_init();
            abpoa_para_t* abpt = abpoa_init_para();
            abpt->align_mode = ABPOA_GLOBAL_MODE;
            abpt->out_msa  = 1;
            abpt->out_cons = 1;
            abpt->out_gfa  = 0;
            abpt->sort_input_seq = 0;   // keep backbone as seq 0 (the spine).
            abpt->progressive_poa = 0;
            abpt->max_n_cons = 1;
            abpt->use_qv = 0;
            abpoa_post_set_para(abpt);

            ab->abs->n_seq = totalSeqs;

            // Seed backbone as sequence 0. On an empty graph,
            // abpoa_add_subgraph_alignment creates nodes 2..L+1; thereafter
            // backbone base position p maps to node id p + 2.
            {
                abpoa_res_t res{};
                res.graph_cigar = nullptr;
                res.n_cigar = 0;
                abpoa_align_sequence_to_subgraph(
                    ab, abpt, 0, 1, backboneSeq.data(), backboneLen, &res);
                abpoa_add_subgraph_alignment(
                    ab, abpt, 0, 1, backboneSeq.data(), nullptr,
                    backboneLen, nullptr, res, 0, totalSeqs, 0);
                if(res.n_cigar) free(res.graph_cigar);
            }

            // Align each member to the backbone subgraph spanning its anchor
            // interval [backboneBegin, backboneEnd).
            for(int i = 0; i < nMembers; i++) {
                const OrientedReadInterval& interval = members[size_t(i)];

                // Inclusive backbone node IDs for the flanking positions.
                const int incBeg = interval.backboneBegin + 2;
                const int incEnd = interval.backboneEnd + 2 - 1;
                if(incBeg > incEnd) continue;

                int excBeg = 0, excEnd = 1;
                abpoa_subgraph_nodes(ab, abpt, incBeg, incEnd, &excBeg, &excEnd);

                abpoa_res_t res{};
                res.graph_cigar = nullptr;
                res.n_cigar = 0;
                abpoa_align_sequence_to_subgraph(
                    ab, abpt, excBeg, excEnd,
                    const_cast<uint8_t*>(interval.seq.data()),
                    int(interval.seq.size()), &res);
                abpoa_add_subgraph_alignment(
                    ab, abpt, excBeg, excEnd,
                    const_cast<uint8_t*>(interval.seq.data()), nullptr,
                    int(interval.seq.size()), nullptr, res, i + 1, totalSeqs, 0);
                if(res.n_cigar) free(res.graph_cigar);
                readsAligned.fetch_add(1);
            }

            // Generate consensus + MSA so the graph is complete, then GFA.
            abpoa_output(ab, abpt, nullptr);

            const string gfaPath = outputPrefix + "window" +
                std::to_string(uint64_t(window.windowId)) + ".gfa";
            FILE* gfaFile = fopen(gfaPath.c_str(), "w");
            if(gfaFile != nullptr) {
                abpoa_generate_gfa(ab, abpt, gfaFile);
                fclose(gfaFile);
                windowsWritten.fetch_add(1);
            } else {
                windowsSkipped.fetch_add(1);
            }

            abpoa_free(ab);
            abpoa_free_para(abpt);
        }
    };

    vector<std::thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t = 0; t < threadCount; t++) {
        threads.emplace_back(worker);
    }
    for(auto& th : threads) {
        th.join();
    }

    cout << timestamp << "computeWindowAbpoaGraphs: wrote " << windowsWritten.load()
         << " GFAs, skipped " << windowsSkipped.load() << " windows, aligned "
         << readsAligned.load() << " member reads." << endl;
}
