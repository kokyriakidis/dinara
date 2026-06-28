// AssemblerWindowAbpoaGraph.cpp
//
// Per-window, per-segment (inter-anchor) abPOA MSA construction.
//
// For each anchor window we build ONE abPOA graph:
//   1. Seed the graph with the backbone read's base sequence spanning the
//      window's backbone anchors. The backbone is sequence 0 (the spine);
//      backbone base position p maps to graph node id p + 2. Seeding an empty
//      graph does no DP, so the backbone may be arbitrarily long.
//   2. For every other member read, find the backbone anchors it shares with
//      the backbone (its "pins"). Between each pair of CONSECUTIVE shared
//      anchors, align only that inter-anchor segment of the read to the
//      backbone subgraph spanning the same two anchors. All segments of one
//      member use the same abPOA read_id, so the member occupies a single MSA
//      row that threads through its shared anchors.
//   3. After all members are added, extract the row-column MSA matrix
//      (abc->msa_base) and write it as a CSV per window.
//
// Why per-segment: abPOA's global alignment DP is O(graph_span x query_len).
// Aligning a whole long-read window in one shot makes a span x span matrix
// (~150k x 150k => ~128 GiB) and aborts. Splitting at shared anchors bounds
// each DP to the (small) inter-anchor gap, so there is NO size cap and no
// member is dropped for size.
//
// Intervals are computed per OrientedReadId (read + strand), never per ReadId;
// a read's two strands are distinct OrientedReadIds with distinct pins.
//
// Anchor/sequence handling mirrors AssemblerAnchorWindowsClean.cpp::createWindow
// and PhasingKmeansAlign.cpp::abpoaMsaRun.

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
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// A shared backbone anchor for one oriented member read.
struct AnchorPin {
    int backbonePos;       // 0-based base position into the backbone sequence.
    uint32_t readOrdinal;  // Marker ordinal of this anchor on the member read.
};

// A member oriented read and its shared-anchor pins (in backbone order).
struct OrientedReadMember {
    OrientedReadId orientedReadId;
    vector<AnchorPin> pins;   // Sorted ascending by backbonePos.
    int backboneSpan = 0;     // pins.back - pins.front; used for add order.
};

inline uint8_t baseToInt(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4; // N
    }
}

// abPOA encodes bases 0..3 = ACGT, 4 = N, and the gap value equals abpt->m.
inline char intToBaseOrGap(uint8_t v, int gapValue) {
    if(int(v) == gapValue) return '-';
    switch(v) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
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
    std::atomic<uint64_t> membersAligned{0};
    std::atomic<uint64_t> segmentsAligned{0};
    std::atomic<uint64_t> segmentsSkipped{0};

    // Extract the base sequence of an oriented read between two marker ordinals
    // (midpoint to midpoint), encoded as 0123/4. Mirrors the extraction used in
    // PhasingKmeansAlign.cpp / createWindow.
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
            // anchor. The backbone sequence spans [bbStartBase, bbEndBase).
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

            // Backbone sequence (0123/4).
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

            // Consume the shared-anchor pins persisted on each member during
            // window construction (readJourneyPos, backboneJourneyPos), and
            // convert them to (backbone base position, read marker ordinal).
            // No journey re-intersection here — that work was done once in
            // createWindow.
            vector<OrientedReadMember> members;
            members.reserve(window.readIntervals.size());

            for(const AnchorWindowReadInterval& ri : window.readIntervals) {
                const OrientedReadId oid = ri.orientedReadId;
                if(oid == backboneOid) continue;
                if(oid.getValue() >= journeys.size()) continue;
                if(ri.sharedPins.size() < 2) continue;
                const auto readMarkers = (*markers)[oid.getValue()];

                vector<AnchorPin> pins;
                pins.reserve(ri.sharedPins.size());
                for(const AnchorWindowSharedPin& sp : ri.sharedPins) {
                    // Backbone base position: midpoint of the anchor on the
                    // backbone, relative to backbone sequence start. Must lie
                    // within the seeded backbone span.
                    if(sp.backboneJourneyPos >= uint32_t(backboneJourney.size())) continue;
                    const Shasta2AnchorId aid = backboneJourney[sp.backboneJourneyPos];
                    const uint32_t bbOrd = anchors.getOrdinal(aid, backboneOid);
                    if(bbOrd >= bbMarkers.size()) continue;
                    const uint32_t bbMid = bbMarkers[bbOrd].position + kHalf;
                    if(bbMid < bbStartBase || bbMid >= bbEndBase) continue;
                    const int backbonePos = int(bbMid - bbStartBase);

                    // Read marker ordinal of the same anchor on this oriented read.
                    const uint32_t readOrd = anchors.getOrdinal(aid, oid);
                    if(readOrd >= readMarkers.size()) continue;

                    pins.push_back(AnchorPin{backbonePos, readOrd});
                }
                if(pins.size() < 2) continue;

                // Pins arrive in backbone order from construction; re-sort
                // defensively in case any were dropped above.
                std::sort(pins.begin(), pins.end(),
                    [](const AnchorPin& a, const AnchorPin& b) {
                        return a.backbonePos < b.backbonePos;
                    });

                OrientedReadMember m;
                m.orientedReadId = oid;
                m.pins = std::move(pins);
                m.backboneSpan = m.pins.back().backbonePos - m.pins.front().backbonePos;
                if(m.backboneSpan <= 0) continue;
                members.push_back(std::move(m));
            }

            // Add members from the biggest backbone span to the smallest, so the
            // longest-overlapping reads shape the graph spine before shorter
            // ones. Ties broken by oriented read id for determinism.
            std::sort(members.begin(), members.end(),
                [](const OrientedReadMember& a, const OrientedReadMember& b) {
                    if(a.backboneSpan != b.backboneSpan) return a.backboneSpan > b.backboneSpan;
                    return a.orientedReadId < b.orientedReadId;
                });

            const int nMembers = int(members.size());
            const int totalSeqs = 1 + nMembers; // backbone + members.

            // --- Build the per-window abPOA graph ---
            abpoa_t* ab = abpoa_init();
            abpoa_para_t* abpt = abpoa_init_para();
            abpt->align_mode = ABPOA_GLOBAL_MODE;
            // sub_aln: align reads to a subgraph (node-range restricted). Without
            // it abpoa_align_sequence_to_subgraph ignores the node bounds.
            abpt->sub_aln = 1;
            abpt->inc_path_score = 1;
            abpt->out_msa  = 1;   // populate abc->msa_base (row-column MSA).
            abpt->out_cons = 1;   // also compute a consensus row.
            abpt->out_gfa  = 0;
            abpt->sort_input_seq = 0;   // keep backbone as seq 0 (the spine).
            abpt->progressive_poa = 0;
            abpt->max_n_cons = 1;
            abpt->use_qv = 0;
            abpoa_post_set_para(abpt);

            ab->abs->n_seq = totalSeqs;

            // Seed backbone as read_id 0. On an empty graph,
            // abpoa_add_subgraph_alignment creates nodes 2..L+1 directly (no DP),
            // so backbone base position p maps to node id p + 2.
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

            // Add each member segment-by-segment between consecutive shared
            // anchors. All segments of one member use the same read_id, so the
            // member fills a single MSA row.
            for(int i = 0; i < nMembers; i++) {
                const OrientedReadMember& m = members[size_t(i)];
                const int readId = i + 1;
                bool anySegment = false;

                for(size_t p = 0; p + 1 < m.pins.size(); p++) {
                    const AnchorPin& a = m.pins[p];
                    const AnchorPin& b = m.pins[p + 1];

                    // Require strictly increasing read ordinal (forward,
                    // monotone) and backbone position. Tangled/duplicate pins
                    // are skipped.
                    if(b.readOrdinal <= a.readOrdinal) {
                        segmentsSkipped.fetch_add(1);
                        continue;
                    }
                    if(b.backbonePos <= a.backbonePos) {
                        segmentsSkipped.fetch_add(1);
                        continue;
                    }

                    // Backbone subgraph for this segment: nodes covering
                    // backbone bases [a.backbonePos, b.backbonePos).
                    const int incBeg = a.backbonePos + 2;
                    const int incEnd = b.backbonePos + 2 - 1;
                    if(incBeg > incEnd) {
                        segmentsSkipped.fetch_add(1);
                        continue;
                    }

                    // Read bases for this segment: [midpoint(a), midpoint(b)).
                    vector<uint8_t> segSeq = extractSeq0123(m.orientedReadId, a.readOrdinal, b.readOrdinal);
                    if(segSeq.empty()) {
                        segmentsSkipped.fetch_add(1);
                        continue;
                    }

                    int excBeg = 0, excEnd = 1;
                    abpoa_subgraph_nodes(ab, abpt, incBeg, incEnd, &excBeg, &excEnd);

                    abpoa_res_t res{};
                    res.graph_cigar = nullptr;
                    res.n_cigar = 0;
                    abpoa_align_sequence_to_subgraph(
                        ab, abpt, excBeg, excEnd,
                        segSeq.data(), int(segSeq.size()), &res);
                    // inc_both_ends = 1: register this read_id at the segment's
                    // BEGIN anchor node. abPOA always excludes the end node, so
                    // for adjacent segments sharing an anchor this fills the
                    // internal anchor columns and makes the member's MSA row
                    // thread continuously through its shared anchors (with 0,
                    // every internal anchor column would be a gap for the read).
                    abpoa_add_subgraph_alignment(
                        ab, abpt, excBeg, excEnd,
                        segSeq.data(), nullptr,
                        int(segSeq.size()), nullptr, res, readId, totalSeqs, 1);
                    if(res.n_cigar) free(res.graph_cigar);
                    segmentsAligned.fetch_add(1);
                    anySegment = true;
                }
                if(anySegment) membersAligned.fetch_add(1);
            }

            // --- Extract the row-column MSA matrix and write CSV ---
            abpoa_generate_rc_msa(ab, abpt);
            const abpoa_cons_t* abc = ab->abc;

            bool wrote = false;
            if(abc != nullptr && abc->msa_base != nullptr && abc->msa_len > 0 && abc->n_seq > 0) {
                const int gapValue = abpt->m; // abPOA stores gaps as value m (=5).
                const string csvPath = outputPrefix + "window" +
                    std::to_string(uint64_t(window.windowId)) + ".msa.csv";
                FILE* f = fopen(csvPath.c_str(), "w");
                if(f != nullptr) {
                    // Header: row label, then one column per MSA column.
                    fprintf(f, "row,orientedReadId,role");
                    for(int c = 0; c < abc->msa_len; c++) fprintf(f, ",c%d", c);
                    fprintf(f, "\n");

                    // Row 0 = backbone; rows 1..nMembers = members in add order.
                    for(int rowReadId = 0; rowReadId < abc->n_seq; rowReadId++) {
                        const char* role = (rowReadId == 0) ? "backbone" : "member";
                        uint64_t oidVal;
                        if(rowReadId == 0) {
                            oidVal = uint64_t(backboneOid.getValue());
                        } else {
                            oidVal = uint64_t(members[size_t(rowReadId - 1)].orientedReadId.getValue());
                        }
                        fprintf(f, "%d,%llu,%s", rowReadId, (unsigned long long)oidVal, role);
                        const uint8_t* rowBases = abc->msa_base[rowReadId];
                        for(int c = 0; c < abc->msa_len; c++) {
                            fprintf(f, ",%c", intToBaseOrGap(rowBases[c], gapValue));
                        }
                        fprintf(f, "\n");
                    }

                    // Append the consensus row(s) if present.
                    for(int consI = 0; consI < abc->n_cons; consI++) {
                        const int row = abc->n_seq + consI;
                        fprintf(f, "%d,-,consensus", row);
                        const uint8_t* rowBases = abc->msa_base[row];
                        for(int c = 0; c < abc->msa_len; c++) {
                            fprintf(f, ",%c", intToBaseOrGap(rowBases[c], gapValue));
                        }
                        fprintf(f, "\n");
                    }
                    fclose(f);
                    wrote = true;
                }
            }
            if(wrote) windowsWritten.fetch_add(1);
            else windowsSkipped.fetch_add(1);

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
         << " MSA CSVs, skipped " << windowsSkipped.load() << " windows, aligned "
         << membersAligned.load() << " members over " << segmentsAligned.load()
         << " segments (" << segmentsSkipped.load() << " segments skipped)." << endl;
}
