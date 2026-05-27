// Split-read classification for SV detection.
//
// After DP chaining with multi-chain extraction, multiple chains may exist
// for the same read-vs-reference pair. This step classifies each chain as:
//   - primary:       best-scoring chain for this read
//   - supplementary: covers a different part of the read (split-read SV signal)
//   - secondary:     overlaps a primary on the read (alternative mapping, not SV)
//
// The logic follows minimap2's mm_set_parent: chains are processed in
// descending score order. A chain becomes a new primary if its query span
// is mostly non-overlapping with all existing primaries. Otherwise it's
// secondary to the primary it overlaps most.

#include "Assembler.hpp"
#include "Alignment.hpp"
#include "AlignmentCandidates.hpp"
#include "Reads.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <vector>

using namespace dinara;
using namespace std;


void Assembler::classifySplitAlignments(
    uint64_t referenceReadCount,
    double maskLevel,
    const string& outputFileName)
{
    performanceLog << timestamp
        << "Classifying split alignments." << endl;
    const auto tBegin = std::chrono::steady_clock::now();

    const auto& candidates = alignmentCandidates.candidates;
    const auto& alignments = alignmentCandidatesAlignmentsData.alignments;
    const uint64_t n = candidates.size();
    DINARA_ASSERT(alignments.size() == n);

    // Group alignments by read pair. Since mcopy can produce multiple
    // entries for the same OrientedReadPair, we need to find groups.
    // Build an index: for each alignment, store (readId, alignmentIndex).
    // We group by the non-reference read ID.

    struct ChainInfo {
        uint64_t alignmentIndex;
        ReadId readId;          // The short read (non-reference).
        ReadId refId;           // The reference read.
        bool isSameStrand;
        uint32_t qs, qe;       // Query (read) span.
        uint32_t ts, te;        // Target (reference) span.
        uint32_t score;         // Chain score (number of ordinal pairs).
        // Classification result.
        enum Class { Primary, Supplementary, Secondary } classification;
        int32_t parentIndex;    // Index within the read's chain group (-1 = self).
    };

    vector<ChainInfo> chains(n);
    for(uint64_t i = 0; i < n; ++i) {
        const auto& cand = candidates[i];
        const auto& al = alignments[i];

        // In OrientedReadPair, readIds[0] < readIds[1].
        // The reference is read 0..referenceReadCount-1.
        ReadId readId, refId;
        if(cand.readIds[0] < ReadId(referenceReadCount)) {
            refId = cand.readIds[0];
            readId = cand.readIds[1];
        } else {
            readId = cand.readIds[0];
            refId = cand.readIds[1];
        }

        chains[i].alignmentIndex = i;
        chains[i].readId = readId;
        chains[i].refId = refId;
        chains[i].isSameStrand = cand.isSameStrand;
        chains[i].qs = al.qs;
        chains[i].qe = al.qe;
        chains[i].ts = al.ts;
        chains[i].te = al.te;
        chains[i].score = uint32_t(al.ordinals.size());
        chains[i].classification = ChainInfo::Primary;
        chains[i].parentIndex = -1;
    }

    // Sort by (readId, descending score) so we process best chains first.
    vector<uint64_t> order(n);
    for(uint64_t i = 0; i < n; ++i) order[i] = i;
    sort(order.begin(), order.end(), [&](uint64_t a, uint64_t b) {
        if(chains[a].readId != chains[b].readId)
            return chains[a].readId < chains[b].readId;
        return chains[a].score > chains[b].score;
    });

    // Process each read's chains.
    uint64_t nPrimary = 0, nSupplementary = 0, nSecondary = 0;
    uint64_t groupStart = 0;

    while(groupStart < n) {
        const ReadId currentRead = chains[order[groupStart]].readId;
        uint64_t groupEnd = groupStart + 1;
        while(groupEnd < n && chains[order[groupEnd]].readId == currentRead) {
            ++groupEnd;
        }

        // Indices within this group that are primaries (in processing order).
        vector<uint64_t> primaries;

        for(uint64_t gi = groupStart; gi < groupEnd; ++gi) {
            auto& chain = chains[order[gi]];

            if(primaries.empty()) {
                // First chain for this read is always primary.
                chain.classification = ChainInfo::Primary;
                chain.parentIndex = -1;
                primaries.push_back(order[gi]);
                ++nPrimary;
                continue;
            }

            // Check query overlap with each existing primary.
            const uint32_t si = chain.qs;
            const uint32_t ei = chain.qe;
            const uint32_t chainLen = ei - si;
            if(chainLen == 0) {
                chain.classification = ChainInfo::Secondary;
                chain.parentIndex = 0;
                ++nSecondary;
                continue;
            }

            // Compute total coverage of this chain's query span by existing primaries.
            // Collect covered intervals, merge, compute uncovered length.
            struct Interval { uint32_t s, e; };
            vector<Interval> covered;
            for(uint64_t pi : primaries) {
                const auto& prim = chains[pi];
                const uint32_t sj = prim.qs;
                const uint32_t ej = prim.qe;
                // Clip to this chain's span.
                if(ej <= si || sj >= ei) continue;
                uint32_t cs = (sj < si) ? si : sj;
                uint32_t ce = (ej > ei) ? ei : ej;
                covered.push_back({cs, ce});
            }

            uint32_t coveredLen = 0;
            if(!covered.empty()) {
                // Sort and merge intervals.
                sort(covered.begin(), covered.end(),
                    [](const Interval& a, const Interval& b) { return a.s < b.s; });
                uint32_t ms = covered[0].s, me = covered[0].e;
                for(size_t ci = 1; ci < covered.size(); ++ci) {
                    if(covered[ci].s <= me) {
                        me = max(me, covered[ci].e);
                    } else {
                        coveredLen += me - ms;
                        ms = covered[ci].s;
                        me = covered[ci].e;
                    }
                }
                coveredLen += me - ms;
            }

            const double overlapFrac = double(coveredLen) / double(chainLen);

            if(overlapFrac <= maskLevel) {
                // Mostly non-overlapping with existing primaries → supplementary.
                chain.classification = ChainInfo::Supplementary;
                chain.parentIndex = -1;
                primaries.push_back(order[gi]);
                ++nSupplementary;
            } else {
                // Overlapping → secondary.
                // Find the primary with the most overlap.
                int32_t bestParent = 0;
                uint32_t bestOverlap = 0;
                for(size_t pi = 0; pi < primaries.size(); ++pi) {
                    const auto& prim = chains[primaries[pi]];
                    uint32_t sj = prim.qs, ej = prim.qe;
                    if(ej <= si || sj >= ei) continue;
                    uint32_t ol = min(ei, ej) - max(si, sj);
                    if(ol > bestOverlap) {
                        bestOverlap = ol;
                        bestParent = int32_t(pi);
                    }
                }
                chain.classification = ChainInfo::Secondary;
                chain.parentIndex = bestParent;
                ++nSecondary;
            }
        }

        groupStart = groupEnd;
    }

    // Write output TSV.
    ofstream tsv(outputFileName);
    if(!tsv) {
        throw runtime_error("Cannot open output file: " + outputFileName);
    }
    tsv << "alignment_id\tread_id\tread_name\tref_read_id\tstrand\t"
           "query_start\tquery_end\tref_start\tref_end\t"
           "n_anchors\tclassification\n";

    for(uint64_t i = 0; i < n; ++i) {
        const auto& chain = chains[i];
        const auto readNameSpan = reads->getReadName(chain.readId);
        const string readName(readNameSpan.data(), readNameSpan.size());

        const char* classStr = "primary";
        if(chain.classification == ChainInfo::Supplementary) classStr = "supplementary";
        else if(chain.classification == ChainInfo::Secondary) classStr = "secondary";

        tsv << chain.alignmentIndex << '\t'
            << chain.readId << '\t'
            << readName << '\t'
            << chain.refId << '\t'
            << (chain.isSameStrand ? '+' : '-') << '\t'
            << chain.qs << '\t'
            << chain.qe << '\t'
            << chain.ts << '\t'
            << chain.te << '\t'
            << chain.score << '\t'
            << classStr << '\n';
    }

    tsv.close();

    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            tEnd - tBegin).count());

    cout << "Split-read classification completed in " << tTotal << " s." << endl;
    cout << "  Total chains: " << n << endl;
    cout << "  Primary: " << nPrimary << endl;
    cout << "  Supplementary (split-read SV): " << nSupplementary << endl;
    cout << "  Secondary: " << nSecondary << endl;
    cout << "  Output: " << outputFileName << endl;

    // Summary of reads with supplementary alignments (split reads).
    // These are the SV-carrying reads.
    uint64_t nSplitReads = 0;
    ReadId prevRead = ReadId(-1);
    for(uint64_t i = 0; i < n; ++i) {
        if(chains[order[i]].classification == ChainInfo::Supplementary &&
           chains[order[i]].readId != prevRead) {
            ++nSplitReads;
            prevRead = chains[order[i]].readId;
        }
    }
    cout << "  Reads with split alignments: " << nSplitReads << endl;

    performanceLog << timestamp
        << "Split-read classification completed in " << tTotal << " s." << endl;
}


// =============================================================================
// Build a Theseus MSA for each reference sequence using chain anchors.
//
// For each reference read:
//   1. Collect all unique anchor positions (marker ordinals on the reference)
//      from all chains involving that reference.
//   2. Sort them to define backbone segment boundaries.
//   3. Extract reference subsequences between consecutive boundaries.
//   4. Create a TheseusMSA with the multi-segment constructor.
//   5. For each chain, extract the short read's subsequence between
//      consecutive shared anchors and call align_from.
//   6. Output the MSA and consensus.
//
// The ordinals in each Alignment are {refOrdinal, readOrdinal} because
// readIdA < readIdB is enforced during chaining, and the reference (read 0)
// is always readIdA.
// =============================================================================
#ifdef DINARA_HAVE_THESEUS
#include <theseus/heuristics.h>
#include <theseus/penalties.h>
#include <theseus/theseus_msa_aligner.h>
#endif

void Assembler::buildSvMSA(
    uint64_t referenceReadCount,
    const string& outputPrefix)
{
#ifndef DINARA_HAVE_THESEUS
    throw runtime_error("buildSvMSA requires the theseus library (not available in this build).");
#else
    const auto tBegin = std::chrono::steady_clock::now();
    performanceLog << timestamp << "Building SV MSA." << endl;

    const auto& candidates = alignmentCandidates.candidates;
    const auto& alignments = alignmentCandidatesAlignmentsData.alignments;
    const uint64_t n = candidates.size();
    DINARA_ASSERT(alignments.size() == n);

    const Reads& readsRef = getReads();
    const auto& markersRef = *markers;
    const uint64_t k = assemblerInfo->k;
    const uint32_t kHalf = uint32_t(k / 2);

    // =========================================================================
    // Step 1: Group chains by reference read ID.
    // =========================================================================
    // For each chain, determine which is the reference and which is the read.
    // readIds[0] < readIds[1], and reference reads are 0..referenceReadCount-1.
    // ordinals[i] = {ordinalA, ordinalB} where A is readIds[0] (reference).

    struct ChainEntry {
        uint64_t chainIndex;
        ReadId refId;
        ReadId readId;
        bool isSameStrand;
    };

    // Map from refId -> list of chain entries.
    vector<vector<ChainEntry>> chainsByRef(referenceReadCount);

    for(uint64_t i = 0; i < n; ++i) {
        const auto& cand = candidates[i];
        const auto& al = alignments[i];
        if(al.ordinals.size() < 2) continue;  // Need at least 2 anchors for a segment.

        ReadId refId, readId;
        if(cand.readIds[0] < ReadId(referenceReadCount)) {
            refId = cand.readIds[0];
            readId = cand.readIds[1];
        } else {
            // Both are short reads — skip (shouldn't happen in svanchors).
            continue;
        }

        chainsByRef[refId].push_back({i, refId, readId, cand.isSameStrand});
    }

    // =========================================================================
    // Step 2: For each reference read, build MSA.
    // =========================================================================
    uint64_t totalMSAs = 0;
    uint64_t totalAlignedReads = 0;
    uint64_t totalAlignedSegments = 0;

    for(ReadId refId = 0; refId < ReadId(referenceReadCount); ++refId) {
        const auto& chainsForRef = chainsByRef[refId];
        if(chainsForRef.empty()) continue;

        const OrientedReadId refOid(refId, 0);
        const auto refMarkers = markersRef[refOid.getValue()];
        const uint32_t refLength = uint32_t(readsRef.getRead(refId).baseCount);

        // Collect all unique reference marker ordinals across all chains.
        vector<uint32_t> allRefOrdinals;
        for(const auto& ce : chainsForRef) {
            const auto& al = alignments[ce.chainIndex];
            for(const auto& ord : al.ordinals) {
                allRefOrdinals.push_back(ord[0]);  // Reference ordinal.
            }
        }

        // Sort and deduplicate.
        sort(allRefOrdinals.begin(), allRefOrdinals.end());
        allRefOrdinals.erase(
            unique(allRefOrdinals.begin(), allRefOrdinals.end()),
            allRefOrdinals.end());

        if(allRefOrdinals.size() < 2) continue;

        // Build a map from reference ordinal -> boundary index.
        unordered_map<uint32_t, uint32_t> ordinalToBoundary;
        ordinalToBoundary.reserve(allRefOrdinals.size());
        for(uint32_t bi = 0; bi < uint32_t(allRefOrdinals.size()); ++bi) {
            ordinalToBoundary[allRefOrdinals[bi]] = bi;
        }

        const uint32_t nBoundaries = uint32_t(allRefOrdinals.size());
        const uint32_t nSegments = nBoundaries - 1;

        // -----------------------------------------------------------------
        // Step 3: Extract backbone segments from the reference.
        // -----------------------------------------------------------------
        // Each segment spans from the midpoint of marker at ordinal[i]
        // to the midpoint of marker at ordinal[i+1].
        vector<string> segmentStrings;
        segmentStrings.reserve(nSegments);
        bool badSegment = false;

        for(uint32_t si = 0; si < nSegments; ++si) {
            const uint32_t ordLeft = allRefOrdinals[si];
            const uint32_t ordRight = allRefOrdinals[si + 1];

            if(ordLeft >= refMarkers.size() || ordRight >= refMarkers.size()) {
                badSegment = true;
                break;
            }

            const uint32_t posLeft = refMarkers[ordLeft].position + kHalf;
            const uint32_t posRight = refMarkers[ordRight].position + kHalf;

            if(posRight <= posLeft) {
                badSegment = true;
                break;
            }

            string seg;
            seg.reserve(posRight - posLeft);
            for(uint32_t pos = posLeft; pos < posRight; ++pos) {
                seg.push_back(readsRef.getOrientedReadBase(refOid, pos).character());
            }
            segmentStrings.push_back(std::move(seg));
        }

        if(badSegment || segmentStrings.empty()) continue;

        // -----------------------------------------------------------------
        // Step 4: Create TheseusMSA with multi-segment constructor.
        // -----------------------------------------------------------------
        vector<string_view> segmentViews;
        segmentViews.reserve(nSegments);
        for(const auto& s : segmentStrings) {
            segmentViews.push_back(s);
        }

        theseus::Penalties penalties(0, 2, 3, 1);  // match=0, mismatch=2, gap_open=3, gap_extend=1
        theseus::Heuristics heuristics(false, false);
        vector<theseus::Graph::NodeId> nodeIds;
        theseus::TheseusMSA aligner(penalties, heuristics, segmentViews, nodeIds, 1);

        // -----------------------------------------------------------------
        // Step 5: Align each chain's read segments into the MSA.
        // -----------------------------------------------------------------
        int seqId = 1;  // 0 is the backbone (reference).
        vector<string> seqNames;
        seqNames.push_back(string("ref_") + to_string(uint32_t(refId)));

        for(const auto& ce : chainsForRef) {
            const auto& al = alignments[ce.chainIndex];
            const ReadId readId = ce.readId;
            const bool isSameStrand = ce.isSameStrand;

            // Determine the oriented read ID for the short read.
            // If same strand, the read is on strand 0. If different strand, strand 1.
            const Strand readStrand = isSameStrand ? 0 : 1;
            const OrientedReadId readOid(readId, readStrand);
            const auto readMarkers = markersRef[readOid.getValue()];

            // Walk the chain's ordinals. Each consecutive pair of ordinals
            // that map to backbone boundaries defines a segment to align.
            // ordinals are {refOrdinal, readOrdinal}.
            const auto& ordinals = al.ordinals;

            // Map each ordinal to its boundary index.
            vector<pair<uint32_t, uint32_t>> boundaryAndReadOrdinal;
            for(const auto& ord : ordinals) {
                auto it = ordinalToBoundary.find(ord[0]);
                if(it != ordinalToBoundary.end()) {
                    boundaryAndReadOrdinal.push_back({it->second, ord[1]});
                }
            }

            if(boundaryAndReadOrdinal.size() < 2) continue;

            // Sort by boundary index (should already be sorted if chain is collinear).
            sort(boundaryAndReadOrdinal.begin(), boundaryAndReadOrdinal.end());

            const auto readNameSpan = readsRef.getReadName(readId);
            const string readName(readNameSpan.data(), readNameSpan.size());
            seqNames.push_back(readName);

            uint32_t alignedSegs = 0;
            for(size_t j = 0; j + 1 < boundaryAndReadOrdinal.size(); ++j) {
                const uint32_t bLeft = boundaryAndReadOrdinal[j].first;
                const uint32_t bRight = boundaryAndReadOrdinal[j + 1].first;
                const uint32_t readOrdLeft = boundaryAndReadOrdinal[j].second;
                const uint32_t readOrdRight = boundaryAndReadOrdinal[j + 1].second;

                if(bRight <= bLeft) continue;
                if(readOrdRight <= readOrdLeft) continue;
                if(readOrdLeft >= readMarkers.size() || readOrdRight >= readMarkers.size()) continue;

                // Extract read subsequence between the two anchor positions.
                const uint32_t readPosLeft = readMarkers[readOrdLeft].position + kHalf;
                const uint32_t readPosRight = readMarkers[readOrdRight].position + kHalf;
                if(readPosRight <= readPosLeft) continue;

                string readSeg;
                readSeg.reserve(readPosRight - readPosLeft);
                for(uint32_t pos = readPosLeft; pos < readPosRight; ++pos) {
                    readSeg.push_back(readsRef.getOrientedReadBase(readOid, pos).character());
                }

                if(readSeg.empty()) continue;

                // Determine start and end nodes in the graph.
                if(bLeft >= nodeIds.size()) continue;
                int endNode = (bRight < nodeIds.size())
                    ? static_cast<int>(nodeIds[bRight])
                    : -1;

                aligner.align_from(
                    readSeg,
                    nodeIds[bLeft],
                    1,      // weight
                    true,   // is_ends_free
                    0,      // start_offset
                    endNode,
                    seqId);

                ++alignedSegs;
                ++totalAlignedSegments;
            }

            if(alignedSegs > 0) {
                ++seqId;
                ++totalAlignedReads;
            } else {
                seqNames.pop_back();
            }
        }

        // -----------------------------------------------------------------
        // Step 6: Output MSA and consensus.
        // -----------------------------------------------------------------
        if(seqId <= 1) continue;  // Only backbone, no reads aligned.

        const string msaFileName = outputPrefix + "_ref" + to_string(uint32_t(refId)) + ".msa";
        {
            ofstream msaOut(msaFileName);
            if(msaOut) {
                aligner.print_as_msa(msaOut, seqId, &seqNames);
            }
        }

        const string consensusFileName = outputPrefix + "_ref" + to_string(uint32_t(refId)) + ".consensus";
        {
            ofstream consOut(consensusFileName);
            if(consOut) {
                vector<int> consensusWeights;
                string consensusSeq;
                string consensusGapped;
                aligner.majority_voting_consensus(consensusWeights, consensusSeq, consensusGapped);
                consOut << ">ref" << uint32_t(refId) << "_consensus\n"
                        << consensusSeq << "\n";
            }
        }

        ++totalMSAs;
        cout << "  Reference " << refId << ": " << nBoundaries << " anchor boundaries, "
             << nSegments << " segments, " << chainsForRef.size() << " chains, "
             << (seqId - 1) << " reads aligned. Output: " << msaFileName << endl;
    }

    const auto tEnd = std::chrono::steady_clock::now();
    const double tTotal = 1.e-9 * double(
        std::chrono::duration_cast<std::chrono::nanoseconds>(tEnd - tBegin).count());

    cout << "SV MSA completed in " << tTotal << " s." << endl;
    cout << "  MSAs built: " << totalMSAs << endl;
    cout << "  Total reads aligned: " << totalAlignedReads << endl;
    cout << "  Total segments aligned: " << totalAlignedSegments << endl;

    performanceLog << timestamp
        << "buildSvMSA completed in " << tTotal << " s. "
        << totalMSAs << " MSAs, " << totalAlignedReads << " reads." << endl;
#endif // DINARA_HAVE_THESEUS
}
