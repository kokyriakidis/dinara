// AssemblerDirectOverlapMSA.cpp
//
// Diagnostic: build a single multi-segment Theseus MSA for one focal read
// using all its direct overlaps from alignmentTable.
// The focal read's journey defines the backbone segments; overlapping reads
// are aligned via align_from using shared anchor positions as boundaries.
// Reads with only one shared anchor are aligned ends-free from that anchor.
// Reports timing and coverage to evaluate feasibility of per-read MSA.

#include "Assembler.hpp"
#include "Reads.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "timestamp.hpp"

#include <theseus/heuristics.h>
#include <theseus/penalties.h>
#include <theseus/theseus_msa_aligner.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

string extractSegment(
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t k,
    OrientedReadId oid,
    uint32_t ordinalA,
    uint32_t ordinalB)
{
    if(ordinalA >= ordinalB) {
        return {};
    }
    const auto readMarkers = markers[oid.getValue()];
    if(ordinalB >= readMarkers.size()) {
        return {};
    }

    const uint32_t kHalf = uint32_t(k / 2);
    const uint32_t beginPos = readMarkers[ordinalA].position + kHalf;
    const uint32_t endPos   = readMarkers[ordinalB].position + kHalf;
    if(endPos <= beginPos) {
        return {};
    }

    string seq;
    seq.reserve(endPos - beginPos);
    for(uint32_t pos = beginPos; pos < endPos; pos++) {
        seq.push_back(reads.getOrientedReadBase(oid, pos).character());
    }
    return seq;
}

// Extract sequence from a marker ordinal extending in one direction,
// bounded by a maximum length.
string extractSegmentToEnd(
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    uint64_t k,
    OrientedReadId oid,
    uint32_t ordinal,
    bool toRight,
    uint32_t maxLength)
{
    const auto readMarkers = markers[oid.getValue()];
    if(ordinal >= readMarkers.size()) {
        return {};
    }

    const uint32_t kHalf = uint32_t(k / 2);
    const uint32_t anchorPos = readMarkers[ordinal].position + kHalf;
    const uint32_t readLength = uint32_t(reads.getRead(oid.getReadId()).baseCount);

    uint32_t beginPos, endPos;
    if(toRight) {
        beginPos = anchorPos;
        endPos = min(anchorPos + maxLength, readLength);
    } else {
        endPos = anchorPos;
        beginPos = (anchorPos > maxLength) ? (anchorPos - maxLength) : 0;
    }

    if(endPos <= beginPos) {
        return {};
    }

    string seq;
    seq.reserve(endPos - beginPos);
    for(uint32_t pos = beginPos; pos < endPos; pos++) {
        seq.push_back(reads.getOrientedReadBase(oid, pos).character());
    }
    return seq;
}

} // anonymous namespace


void Assembler::testDirectOverlapMSA(
    const shared_ptr<Shasta2Anchors>& shasta2Anchors,
    const shared_ptr<Shasta2Journeys>& shasta2Journeys,
    ReadId focalReadId)
{
    // Reset Theseus timing counters so we only measure this function.
    theseus::TheseusMSA::reset_timing_counters();

    const auto totalBegin = chrono::steady_clock::now();

    const Reads& readsRef = getReads();
    const auto& markersRef = *markers;
    const uint64_t k = assemblerInfo->k;
    const OrientedReadId focalOid(focalReadId, 0);

    cout << timestamp << "[DirectOverlapMSA] Focal read " << focalOid
         << " length " << readsRef.getRead(focalReadId).baseCount << " bp" << endl;

    // Get the focal read's journey.
    if(focalOid.getValue() >= shasta2Journeys->size()) {
        cout << timestamp << "[DirectOverlapMSA] Focal read has no journey." << endl;
        return;
    }
    const auto focalJourney = (*shasta2Journeys)[focalOid];
    if(focalJourney.size() < 2) {
        cout << timestamp << "[DirectOverlapMSA] Focal journey too short: "
             << focalJourney.size() << " anchors." << endl;
        return;
    }

    cout << timestamp << "[DirectOverlapMSA] Focal journey: "
         << focalJourney.size() << " anchors." << endl;

    // Build backbone segments between consecutive journey anchors.
    const uint32_t nAnchors = uint32_t(focalJourney.size());
    const uint32_t nSegments = nAnchors - 1;

    vector<string> segmentStrings;
    segmentStrings.reserve(nSegments);

    for(uint32_t i = 0; i < nSegments; i++) {
        const Shasta2AnchorId leftAnchorId  = focalJourney[i];
        const Shasta2AnchorId rightAnchorId = focalJourney[i + 1];

        const uint32_t leftOrdinal  = shasta2Anchors->getOrdinal(leftAnchorId, focalOid);
        const uint32_t rightOrdinal = shasta2Anchors->getOrdinal(rightAnchorId, focalOid);

        string seg = extractSegment(readsRef, markersRef, k, focalOid, leftOrdinal, rightOrdinal);
        if(seg.empty()) {
            cout << timestamp << "[DirectOverlapMSA] Empty backbone segment " << i
                 << " (ordinals " << leftOrdinal << "->" << rightOrdinal << "), aborting." << endl;
            return;
        }
        segmentStrings.push_back(std::move(seg));
    }

    uint64_t totalBackboneBases = 0;
    for(const auto& s : segmentStrings) {
        totalBackboneBases += s.size();
    }
    cout << timestamp << "[DirectOverlapMSA] " << nSegments << " backbone segments, "
         << totalBackboneBases << " total bases." << endl;

    // Build string_view vector for the multi-segment constructor.
    vector<string_view> segmentViews;
    segmentViews.reserve(nSegments);
    for(const auto& s : segmentStrings) {
        segmentViews.push_back(s);
    }

    // Create the multi-segment MSA.
    const auto msaConstructBegin = chrono::steady_clock::now();
    theseus::Penalties penalties(0, 2, 3, 1);
    theseus::Heuristics heuristics(false, false);
    vector<theseus::Graph::NodeId> nodeIds;
    theseus::TheseusMSA aligner(penalties, heuristics, segmentViews, nodeIds, 1);
    const auto msaConstructEnd = chrono::steady_clock::now();
    const double constructSeconds = chrono::duration<double>(msaConstructEnd - msaConstructBegin).count();

    cout << timestamp << "[DirectOverlapMSA] TheseusMSA created with "
         << nodeIds.size() << " segment nodes in "
         << constructSeconds << "s." << endl;

    // Build a map from Shasta2AnchorId -> boundary index for the focal journey.
    unordered_map<uint64_t, uint32_t> anchorToBoundary;
    anchorToBoundary.reserve(nAnchors);
    for(uint32_t i = 0; i < nAnchors; i++) {
        anchorToBoundary[uint64_t(focalJourney[i])] = i;
    }

    // Gather all direct overlaps from alignmentTable.
    const auto discoveryBegin = chrono::steady_clock::now();
    const auto& aTable = getAlignmentTable();
    const auto focalAlignments = aTable[focalOid.getValue()];

    // For each overlapping read, find which backbone anchors it shares.
    struct BoundaryHit {
        uint32_t boundaryIndex;
        uint32_t ordinal;  // marker ordinal on the overlapping read
    };

    struct OverlapReadInfo {
        OrientedReadId oid;
        vector<BoundaryHit> hits;
        uint32_t baseSpan = 0;
    };

    vector<OverlapReadInfo> overlapReads;
    unordered_set<uint64_t> seenReads;

    uint32_t totalAlignments = 0;
    uint32_t skippedDeleted = 0;
    uint32_t skippedDuplicate = 0;
    uint32_t skippedNoJourney = 0;
    uint32_t readsWithZeroHits = 0;
    uint32_t readsWithOneHit = 0;
    uint32_t readsWithTwoOrMoreHits = 0;

    for(uint32_t alignmentId : focalAlignments) {
        ++totalAlignments;
        const AlignmentData& ad = alignmentData[alignmentId];

        // Skip deleted alignments (non-phase reasons).
        constexpr AlignmentData::DeleteReasonMask nonPhase =
            ~AlignmentData::DeleteReasonPhase;
        if((ad.deleteReasons0 & nonPhase) || (ad.deleteReasons1 & nonPhase)) {
            ++skippedDeleted;
            continue;
        }

        // Get the partner oriented read.
        const OrientedReadId partnerOid = ad.getOther(focalOid);
        const uint64_t partnerValue = partnerOid.getValue();

        // Skip duplicates.
        if(!seenReads.insert(partnerValue).second) {
            ++skippedDuplicate;
            continue;
        }

        // Skip if partner has no journey.
        if(partnerValue >= shasta2Journeys->size()) {
            ++skippedNoJourney;
            continue;
        }
        const auto partnerJourney = (*shasta2Journeys)[partnerOid];
        if(partnerJourney.empty()) {
            ++skippedNoJourney;
            continue;
        }

        // Walk the partner's journey and find anchors shared with the focal backbone.
        vector<BoundaryHit> hits;
        for(uint32_t jp = 0; jp < partnerJourney.size(); jp++) {
            const Shasta2AnchorId anchorId = partnerJourney[jp];
            auto it = anchorToBoundary.find(uint64_t(anchorId));
            if(it != anchorToBoundary.end()) {
                const uint32_t ordinal = shasta2Anchors->getOrdinal(anchorId, partnerOid);
                if(ordinal != invalid<uint32_t>) {
                    hits.push_back({it->second, ordinal});
                }
            }
        }

        if(hits.empty()) {
            ++readsWithZeroHits;
            continue;
        }

        // Sort by boundary index.
        sort(hits.begin(), hits.end(),
            [](const BoundaryHit& a, const BoundaryHit& b) {
                return a.boundaryIndex < b.boundaryIndex;
            });

        // Compute base span for sorting.
        const auto partnerMarkers = markersRef[partnerOid.getValue()];
        uint32_t firstPos = partnerMarkers[hits.front().ordinal].position;
        uint32_t lastPos  = partnerMarkers[hits.back().ordinal].position;
        uint32_t span = (lastPos > firstPos) ? (lastPos - firstPos) : 0;

        if(hits.size() == 1) {
            ++readsWithOneHit;
        } else {
            ++readsWithTwoOrMoreHits;
        }

        overlapReads.push_back({partnerOid, std::move(hits), span});
    }

    // Sort by base span descending (longest overlaps first).
    sort(overlapReads.begin(), overlapReads.end(),
        [](const OverlapReadInfo& a, const OverlapReadInfo& b) {
            return a.baseSpan > b.baseSpan;
        });

    const auto discoveryEnd = chrono::steady_clock::now();
    const double discoverySeconds = chrono::duration<double>(discoveryEnd - discoveryBegin).count();

    cout << timestamp << "[DirectOverlapMSA] Overlap discovery:"
         << " time=" << fixed << setprecision(6) << discoverySeconds << "s"
         << defaultfloat
         << " totalAlignments=" << totalAlignments
         << " skippedDeleted=" << skippedDeleted
         << " skippedDuplicate=" << skippedDuplicate
         << " skippedNoJourney=" << skippedNoJourney
         << " readsWithZeroHits=" << readsWithZeroHits
         << " readsWithOneHit=" << readsWithOneHit
         << " readsWithTwoOrMoreHits=" << readsWithTwoOrMoreHits
         << " totalOverlapReads=" << overlapReads.size()
         << endl;

    // Align each overlapping read's segments into the POA graph.
    // Detailed timing: separate sequence extraction from align_from calls.
    const auto alignBegin = chrono::steady_clock::now();
    uint32_t alignedReads = 0;
    uint32_t alignedMultiAnchorReads = 0;
    uint32_t alignedSingleAnchorReads = 0;
    uint32_t alignedSegments = 0;
    uint64_t totalAlignBases = 0;
    double maxSegmentTime = 0.0;
    uint32_t skippedEmptySegments = 0;
    uint32_t skippedBadOrdinals = 0;

    // Phase timing accumulators.
    double totalExtractSeconds = 0.0;
    double totalAlignFromSeconds = 0.0;

    // Segment length buckets for align_from time distribution.
    // Buckets: [0,100), [100,500), [500,1000), [1000,2000), [2000,5000), [5000,+inf)
    struct TimeBucket {
        uint32_t count = 0;
        double totalSeconds = 0.0;
        uint64_t totalBases = 0;
    };
    constexpr uint32_t bucketBounds[] = {100, 500, 1000, 2000, 5000};
    constexpr uint32_t nBuckets = 6;
    TimeBucket buckets[nBuckets];

    auto getBucket = [&](uint64_t seqLen) -> uint32_t {
        for(uint32_t b = 0; b < nBuckets - 1; b++) {
            if(seqLen < bucketBounds[b]) return b;
        }
        return nBuckets - 1;
    };

    // Precompute cumulative backbone segment lengths for distance lookups.
    // cumulativeLen[i] = total backbone bases from boundary 0 to boundary i.
    vector<uint64_t> cumulativeLen(nAnchors, 0);
    for(uint32_t i = 0; i < nSegments; i++) {
        cumulativeLen[i + 1] = cumulativeLen[i] + segmentStrings[i].size();
    }

    // Minimum segment length for emitting an align_from call.
    // Consecutive shared anchors closer than this on the backbone are merged.
    constexpr uint32_t minMergedSegmentLength = 500;

    // Max extension length for single-anchor reads.
    const uint32_t avgSegmentLength = uint32_t(totalBackboneBases / max<uint32_t>(nSegments, 1));
    const uint32_t maxExtensionLength = avgSegmentLength * 5;

    uint32_t mergedAwayAnchors = 0;
    uint32_t emittedMergedSegments = 0;

    int readSeqId = 1;  // 0 is the backbone
    vector<string> msaSeqNames;
    msaSeqNames.push_back(focalOid.getString());
    vector<uint64_t> msaSeqIds;
    msaSeqIds.push_back(focalOid.getValue());

    for(const auto& overlapRead : overlapReads) {
        const OrientedReadId oid = overlapRead.oid;
        const auto& hits = overlapRead.hits;

        uint32_t readSegments = 0;

        if(hits.size() >= 2) {
            // Multi-anchor read: merge consecutive short segments.
            // Walk hits and accumulate backbone distance. Emit a segment
            // when the accumulated distance exceeds minMergedSegmentLength
            // or we reach the last hit.
            size_t segStart = 0;  // index into hits for current merged segment start

            for(size_t hi = 1; hi < hits.size(); hi++) {
                const uint32_t startBoundary = hits[segStart].boundaryIndex;
                const uint32_t curBoundary = hits[hi].boundaryIndex;

                if(curBoundary <= startBoundary) {
                    continue;
                }

                uint64_t backboneDist = cumulativeLen[curBoundary] - cumulativeLen[startBoundary];
                bool isLast = (hi + 1 == hits.size());

                if(backboneDist >= minMergedSegmentLength || isLast) {
                    // Emit merged segment from segStart to hi.
                    const uint32_t prevBoundary = hits[segStart].boundaryIndex;
                    const uint32_t nextBoundary = hits[hi].boundaryIndex;
                    const uint32_t prevOrdinal = hits[segStart].ordinal;
                    const uint32_t nextOrdinal = hits[hi].ordinal;

                    if(nextBoundary > prevBoundary && prevBoundary < nodeIds.size()
                       && nextOrdinal > prevOrdinal) {

                        auto tExtract0 = chrono::steady_clock::now();
                        string readSeq = extractSegment(
                            readsRef, markersRef, k, oid, prevOrdinal, nextOrdinal);
                        auto tExtract1 = chrono::steady_clock::now();
                        totalExtractSeconds += chrono::duration<double>(tExtract1 - tExtract0).count();

                        if(readSeq.empty()) {
                            skippedEmptySegments++;
                        } else {
                            int endNode = (nextBoundary < nodeIds.size())
                                ? static_cast<int>(nodeIds[nextBoundary])
                                : -1;

                            auto t0 = chrono::steady_clock::now();
                            aligner.align_from(
                                readSeq,
                                nodeIds[prevBoundary],
                                1,     // weight
                                true,  // is_ends_free
                                0,     // start_offset
                                endNode,
                                readSeqId);
                            auto t1 = chrono::steady_clock::now();
                            double elapsed = chrono::duration<double>(t1 - t0).count();
                            totalAlignFromSeconds += elapsed;
                            maxSegmentTime = max(maxSegmentTime, elapsed);
                            totalAlignBases += readSeq.size();

                            uint32_t b = getBucket(readSeq.size());
                            buckets[b].count++;
                            buckets[b].totalSeconds += elapsed;
                            buckets[b].totalBases += readSeq.size();

                            // Count merged-away intermediate anchors.
                            uint32_t intermediateAnchors = uint32_t(hi - segStart - 1);
                            mergedAwayAnchors += intermediateAnchors;
                            if(intermediateAnchors > 0) {
                                emittedMergedSegments++;
                            }

                            readSegments++;
                            alignedSegments++;
                        }
                    } else {
                        skippedBadOrdinals++;
                    }

                    segStart = hi;
                }
            }

            if(readSegments > 0) {
                alignedMultiAnchorReads++;
            }
        } else {
            // Single-anchor read: align ends-free from the shared anchor.
            const uint32_t boundary = hits[0].boundaryIndex;
            const uint32_t ordinal = hits[0].ordinal;

            if(boundary < nodeIds.size()) {
                const auto partnerMarkers = markersRef[oid.getValue()];
                const uint32_t anchorPos = partnerMarkers[ordinal].position;
                const uint32_t readLength = uint32_t(readsRef.getRead(oid.getReadId()).baseCount);

                // Extend toward whichever side has more sequence.
                uint32_t rightExtent = readLength - anchorPos;
                uint32_t leftExtent = anchorPos;
                bool toRight = (rightExtent >= leftExtent);

                auto tExtract0 = chrono::steady_clock::now();
                string readSeq = extractSegmentToEnd(
                    readsRef, markersRef, k, oid, ordinal, toRight, maxExtensionLength);
                auto tExtract1 = chrono::steady_clock::now();
                totalExtractSeconds += chrono::duration<double>(tExtract1 - tExtract0).count();

                if(!readSeq.empty()) {
                    auto t0 = chrono::steady_clock::now();
                    aligner.align_from(
                        readSeq,
                        nodeIds[boundary],
                        1,     // weight
                        true,  // is_ends_free
                        0,     // start_offset
                        -1,    // end_node: sink (ends-free)
                        readSeqId);
                    auto t1 = chrono::steady_clock::now();
                    double elapsed = chrono::duration<double>(t1 - t0).count();
                    totalAlignFromSeconds += elapsed;
                    maxSegmentTime = max(maxSegmentTime, elapsed);
                    totalAlignBases += readSeq.size();

                    uint32_t b = getBucket(readSeq.size());
                    buckets[b].count++;
                    buckets[b].totalSeconds += elapsed;
                    buckets[b].totalBases += readSeq.size();

                    readSegments++;
                    alignedSegments++;
                    alignedSingleAnchorReads++;
                }
            }
        }

        if(readSegments > 0) {
            msaSeqNames.push_back(oid.getString());
            msaSeqIds.push_back(oid.getValue());
            alignedReads++;
            readSeqId++;
        }
    }

    const auto alignEnd = chrono::steady_clock::now();
    const double alignSeconds = chrono::duration<double>(alignEnd - alignBegin).count();

    const auto totalEnd = chrono::steady_clock::now();
    const double totalSeconds = chrono::duration<double>(totalEnd - totalBegin).count();

    cout << timestamp << "[DirectOverlapMSA] Alignment complete."
         << " alignedReads=" << alignedReads
         << " multiAnchorReads=" << alignedMultiAnchorReads
         << " singleAnchorReads=" << alignedSingleAnchorReads
         << " alignedSegments=" << alignedSegments
         << " totalAlignBases=" << totalAlignBases
         << " skippedEmptySegments=" << skippedEmptySegments
         << " skippedBadOrdinals=" << skippedBadOrdinals
         << " mergedAwayAnchors=" << mergedAwayAnchors
         << " emittedMergedSegments=" << emittedMergedSegments
         << defaultfloat << endl;

    cout << timestamp << "[DirectOverlapMSA] Timing breakdown:"
         << " constructPOA=" << fixed << setprecision(6) << constructSeconds << "s"
         << " extractSeq=" << totalExtractSeconds << "s"
         << " alignFrom=" << totalAlignFromSeconds << "s"
         << " loopOverhead=" << (alignSeconds - totalExtractSeconds - totalAlignFromSeconds) << "s"
         << " totalAlignLoop=" << alignSeconds << "s"
         << " totalSeconds=" << totalSeconds << "s"
         << defaultfloat << endl;

    cout << timestamp << "[DirectOverlapMSA] Per-segment avg:"
         << " extract=" << fixed << setprecision(4)
         << (alignedSegments > 0 ? totalExtractSeconds / alignedSegments * 1e6 : 0) << "us"
         << " alignFrom=" << (alignedSegments > 0 ? totalAlignFromSeconds / alignedSegments * 1e6 : 0) << "us"
         << " maxAlignFrom=" << maxSegmentTime * 1e3 << "ms"
         << defaultfloat << endl;

    // Print segment length buckets.
    {
        const char* bucketLabels[] = {"<100bp", "100-500bp", "500-1000bp", "1000-2000bp", "2000-5000bp", ">=5000bp"};
        cout << timestamp << "[DirectOverlapMSA] align_from time by segment length:" << endl;
        for(uint32_t b = 0; b < nBuckets; b++) {
            if(buckets[b].count == 0) continue;
            cout << timestamp << "  " << bucketLabels[b]
                 << ": n=" << buckets[b].count
                 << " totalTime=" << fixed << setprecision(6) << buckets[b].totalSeconds << "s"
                 << " avgTime=" << setprecision(4) << (buckets[b].totalSeconds / buckets[b].count * 1e6) << "us"
                 << " totalBases=" << buckets[b].totalBases
                 << " avgBases=" << (buckets[b].totalBases / buckets[b].count)
                 << defaultfloat << endl;
        }
    }

    // Print Theseus internal timing breakdown.
    {
        auto tc = theseus::TheseusMSA::get_timing_counters();
        cout << timestamp << "[DirectOverlapMSA] Theseus internals:"
             << " subgraph=" << fixed << setprecision(6) << (tc.subgraph_ns / 1e9) << "s"
             << " newAlignment=" << (tc.new_alignment_ns / 1e9) << "s"
             << " dpLoop=" << (tc.dp_loop_ns / 1e9) << "s"
             << " poaUpdate=" << (tc.poa_update_ns / 1e9) << "s"
             << " calls=" << tc.align_calls
             << " totalScore=" << tc.total_score
             << " avgScore=" << (tc.align_calls > 0 ? tc.total_score / tc.align_calls : 0)
             << defaultfloat << endl;
        theseus::TheseusMSA::reset_timing_counters();
    }

    // Write MSA output.
    {
        const string msaPath = "DirectOverlapMSA_read" + to_string(focalReadId) + ".fasta";
        ofstream msaFile(msaPath);
        aligner.print_as_msa(msaFile, readSeqId, &msaSeqNames);
        cout << timestamp << "[DirectOverlapMSA] MSA written to " << msaPath
             << " (" << readSeqId << " sequences)." << endl;
    }

    // Memory report.
    {
        ifstream procStatus("/proc/self/status");
        string line;
        while(getline(procStatus, line)) {
            if(line.find("VmRSS") == 0 || line.find("VmPeak") == 0) {
                cout << timestamp << "[DirectOverlapMSA] " << line << endl;
            }
        }
    }
}
