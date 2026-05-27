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

    // Store classifications for use by buildSvMSA.
    chainClassifications.resize(n);
    for(uint64_t i = 0; i < n; ++i) {
        switch(chains[i].classification) {
            case ChainInfo::Primary:
                chainClassifications[i] = ChainClassification::Primary;
                break;
            case ChainInfo::Supplementary:
                chainClassifications[i] = ChainClassification::Supplementary;
                break;
            case ChainInfo::Secondary:
                chainClassifications[i] = ChainClassification::Secondary;
                break;
        }
    }

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

        // Skip secondary chains — they overlap a primary on the query
        // and would double-align the same read region.
        if(i < chainClassifications.size()
           && chainClassifications[i] == ChainClassification::Secondary) {
            continue;
        }

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
        // Step 5: Group chains by read, classify SV type, sort, then align.
        // -----------------------------------------------------------------
        // Group chains by readId so we can classify each read's SV type
        // from its set of chains (primary + supplementary).
        //
        // SV type classification per read:
        //   - Deletion: read has supplementary chains and the gap between
        //     chains on the reference is larger than on the query.
        //   - Insertion: gap on the query is larger than on the reference.
        //   - Reference-like: only a primary chain (no supplementary).
        //
        // Alignment order: deletions first (create shortcut paths in the
        // graph), then insertions (create longer alternative paths), then
        // reference-like reads (reinforce the backbone).

        enum class SvType { Deletion = 0, Inversion = 1, Insertion = 2, ReferenceLike = 3 };

        struct ReadGroup {
            ReadId readId;
            SvType svType;
            int64_t svSize;  // |refSpan - readSpan|, larger = bigger SV.
            uint32_t breakpointRefPos;  // Reference position (bp) of the SV breakpoint.
            int32_t clusterId;  // SV cluster assignment (-1 = unclustered).
            vector<size_t> chainIndicesInRef;  // Indices into chainsForRef.
        };

        // Group chains by readId.
        unordered_map<uint32_t, vector<size_t>> readToChainIndices;
        for(size_t ci = 0; ci < chainsForRef.size(); ++ci) {
            readToChainIndices[uint32_t(chainsForRef[ci].readId)].push_back(ci);
        }

        vector<ReadGroup> readGroups;
        readGroups.reserve(readToChainIndices.size());

        for(auto& [rid, indices] : readToChainIndices) {
            ReadGroup rg;
            rg.readId = ReadId(rid);
            rg.chainIndicesInRef = std::move(indices);
            rg.svSize = 0;
            rg.breakpointRefPos = 0;
            rg.clusterId = -1;

            // Check for inversion: chains on both strands.
            const auto& firstCeClassify = chainsForRef[rg.chainIndicesInRef[0]];
            const Strand primaryStrand = firstCeClassify.isSameStrand ? 0 : 1;
            bool hasPrimaryStrand = false;
            bool hasOppositeStrand = false;
            for(size_t ci : rg.chainIndicesInRef) {
                const Strand s = chainsForRef[ci].isSameStrand ? 0 : 1;
                if(s == primaryStrand) hasPrimaryStrand = true;
                else hasOppositeStrand = true;
            }

            // Compute per-chain reference spans for breakpoint detection.
            // Each chain's ref extent: [minRefPos, maxRefPos].
            struct ChainRefSpan {
                uint32_t minRefPos;
                uint32_t maxRefPos;
                size_t chainIdxInRef;
            };
            vector<ChainRefSpan> chainSpans;
            for(size_t ci : rg.chainIndicesInRef) {
                const auto& al = alignments[chainsForRef[ci].chainIndex];
                uint32_t cMinRef = UINT32_MAX, cMaxRef = 0;
                for(const auto& ord : al.ordinals) {
                    if(ord[0] < refMarkers.size()) {
                        const uint32_t pos = refMarkers[ord[0]].position;
                        cMinRef = std::min(cMinRef, pos);
                        cMaxRef = std::max(cMaxRef, pos);
                    }
                }
                if(cMinRef <= cMaxRef) {
                    chainSpans.push_back({cMinRef, cMaxRef, ci});
                }
            }
            // Sort chains by reference start position.
            sort(chainSpans.begin(), chainSpans.end(),
                [](const ChainRefSpan& a, const ChainRefSpan& b) {
                    return a.minRefPos < b.minRefPos;
                });

            if(hasPrimaryStrand && hasOppositeStrand) {
                // Inversion: chains on both strands.
                // svSize = total reference span of the inverted chains.
                // breakpoint = start of the first inverted chain on the reference.
                uint32_t invMinRef = UINT32_MAX, invMaxRef = 0;
                for(size_t ci : rg.chainIndicesInRef) {
                    const Strand s = chainsForRef[ci].isSameStrand ? 0 : 1;
                    if(s == primaryStrand) continue;  // Only measure inverted chains.
                    const auto& al = alignments[chainsForRef[ci].chainIndex];
                    for(const auto& ord : al.ordinals) {
                        if(ord[0] < refMarkers.size()) {
                            const uint32_t pos = refMarkers[ord[0]].position;
                            invMinRef = std::min(invMinRef, pos);
                            invMaxRef = std::max(invMaxRef, pos);
                        }
                    }
                }
                if(invMinRef != UINT32_MAX && invMaxRef > invMinRef) {
                    rg.svSize = int64_t(invMaxRef) - int64_t(invMinRef);
                    rg.breakpointRefPos = invMinRef;
                }
                rg.svType = SvType::Inversion;
            } else {
                // No inversion — classify by examining gaps both
                // BETWEEN chains and WITHIN chains.
                //
                // Between chains: for each pair of consecutive chains
                // (sorted by ref pos), compare refGap vs readGap.
                //
                // Within chains: for each pair of consecutive anchors
                // in a chain, compare the ref gap vs read gap. This
                // catches insertions where the chain doesn't split
                // (the DP bridges the insertion with a gap penalty).

                const OrientedReadId classifyOid(rg.readId, primaryStrand);
                const auto classifyReadMarkers = markersRef[classifyOid.getValue()];

                // Collect all anchor pairs (refPos, readPos) across all
                // primary-strand chains, sorted by refPos.
                struct AnchorPoint {
                    uint32_t refPos;
                    uint32_t readPos;
                };
                vector<AnchorPoint> allAnchors;

                for(size_t ci : rg.chainIndicesInRef) {
                    const Strand s = chainsForRef[ci].isSameStrand ? 0 : 1;
                    if(s != primaryStrand) continue;
                    const auto& al = alignments[chainsForRef[ci].chainIndex];
                    for(const auto& ord : al.ordinals) {
                        if(ord[0] < refMarkers.size() && ord[1] < classifyReadMarkers.size()) {
                            allAnchors.push_back({
                                refMarkers[ord[0]].position,
                                classifyReadMarkers[ord[1]].position
                            });
                        }
                    }
                }

                // Sort by reference position.
                sort(allAnchors.begin(), allAnchors.end(),
                    [](const AnchorPoint& a, const AnchorPoint& b) {
                        return a.refPos < b.refPos;
                    });

                if(allAnchors.size() >= 2) {
                    // Compute the diagonal offset (readPos - refPos) at
                    // each anchor. An SV causes a shift in the diagonal:
                    //   Deletion: diagonal decreases (read loses sequence)
                    //   Insertion: diagonal increases (read gains sequence)
                    //
                    // Find the maximum diagonal change between any two
                    // anchors. This catches both between-chain gaps and
                    // within-chain insertions/deletions.

                    // Compute diagonal at each anchor.
                    vector<int64_t> diagonals(allAnchors.size());
                    for(size_t ai = 0; ai < allAnchors.size(); ++ai) {
                        diagonals[ai] = int64_t(allAnchors[ai].readPos)
                                      - int64_t(allAnchors[ai].refPos);
                    }

                    // Find the maximum diagonal change.
                    // This is the max(diag) - min(diag) with the
                    // constraint that min comes before max (insertion)
                    // or max comes before min (deletion).
                    int64_t minDiag = diagonals[0];
                    int64_t maxDiag = diagonals[0];
                    size_t minDiagIdx = 0;
                    size_t maxDiagIdx = 0;

                    // Track running min and the max rise from it (insertion).
                    int64_t runningMin = diagonals[0];
                    size_t runningMinIdx = 0;
                    int64_t bestRise = 0;  // max(diag[j] - diag[i]) for j > i
                    uint32_t bestRiseBreakpoint = 0;

                    // Track running max and the max drop from it (deletion).
                    int64_t runningMax = diagonals[0];
                    size_t runningMaxIdx = 0;
                    int64_t bestDrop = 0;  // max(diag[i] - diag[j]) for j > i
                    uint32_t bestDropBreakpoint = 0;

                    for(size_t ai = 1; ai < allAnchors.size(); ++ai) {
                        // Check for insertion (diagonal rises).
                        const int64_t rise = diagonals[ai] - runningMin;
                        if(rise > bestRise) {
                            bestRise = rise;
                            // Breakpoint is where the diagonal starts rising.
                            bestRiseBreakpoint = allAnchors[runningMinIdx].refPos;
                        }
                        if(diagonals[ai] < runningMin) {
                            runningMin = diagonals[ai];
                            runningMinIdx = ai;
                        }

                        // Check for deletion (diagonal drops).
                        const int64_t drop = runningMax - diagonals[ai];
                        if(drop > bestDrop) {
                            bestDrop = drop;
                            bestDropBreakpoint = allAnchors[runningMaxIdx].refPos;
                        }
                        if(diagonals[ai] > runningMax) {
                            runningMax = diagonals[ai];
                            runningMaxIdx = ai;
                        }
                    }

                    // Minimum SV size to report (filter noise from
                    // marker resolution).
                    const int64_t minSvSize = 20;

                    if(bestRise > bestDrop && bestRise > minSvSize) {
                        rg.svType = SvType::Insertion;
                        rg.svSize = bestRise;
                        rg.breakpointRefPos = bestRiseBreakpoint;
                    } else if(bestDrop > minSvSize) {
                        rg.svType = SvType::Deletion;
                        rg.svSize = bestDrop;
                        rg.breakpointRefPos = bestDropBreakpoint;
                    } else {
                        rg.svType = SvType::ReferenceLike;
                    }
                } else {
                    rg.svType = SvType::ReferenceLike;
                }
            }

            readGroups.push_back(std::move(rg));
        }

        // -----------------------------------------------------------------
        // Step 5b: Cluster SV-carrying reads by breakpoint position and
        // SV size, following the hifiasm detectSVSites pattern.
        //
        // Within each SV type (Deletion, Insertion, Inversion), reads
        // whose breakpoints fall within SV_WINDOW bp and whose SV sizes
        // are within SV_SIZE_RATIO of each other are assigned to the
        // same cluster. ReferenceLike reads get clusterId = -1.
        // -----------------------------------------------------------------
        constexpr uint32_t SV_WINDOW = 50;
        constexpr double SV_SIZE_RATIO = 0.20;
        constexpr uint32_t SV_MIN_SUPPORT = 1;  // Minimum reads per cluster.

        // Separate indices by SV type for independent clustering.
        vector<size_t> svIndices;  // Indices into readGroups for SV reads.
        for(size_t i = 0; i < readGroups.size(); ++i) {
            if(readGroups[i].svType != SvType::ReferenceLike) {
                svIndices.push_back(i);
            }
        }

        // Sort SV reads by (svType, breakpointRefPos) for clustering.
        sort(svIndices.begin(), svIndices.end(),
            [&readGroups](size_t a, size_t b) {
                if(readGroups[a].svType != readGroups[b].svType)
                    return static_cast<int>(readGroups[a].svType)
                         < static_cast<int>(readGroups[b].svType);
                return readGroups[a].breakpointRefPos < readGroups[b].breakpointRefPos;
            });

        int32_t nextClusterId = 0;

        for(size_t ii = 0; ii < svIndices.size(); ) {
            // Start a new cluster seed from svIndices[ii].
            const size_t seedIdx = svIndices[ii];
            const SvType seedType = readGroups[seedIdx].svType;
            const uint32_t seedPos = readGroups[seedIdx].breakpointRefPos;
            const int64_t seedSize = readGroups[seedIdx].svSize;

            // Gather all reads within the window that match type and size.
            vector<size_t> clusterMembers;
            size_t jj = ii;
            while(jj < svIndices.size()) {
                const size_t idx = svIndices[jj];
                if(readGroups[idx].svType != seedType) break;
                if(readGroups[idx].breakpointRefPos >= seedPos + SV_WINDOW) break;

                // Size tolerance check.
                const int64_t sz = readGroups[idx].svSize;
                const int64_t refSz = (seedSize > 0) ? seedSize : 1;
                const int64_t diff = std::abs(sz - seedSize);
                if(diff <= int64_t(std::abs(refSz) * SV_SIZE_RATIO) || seedSize == 0) {
                    clusterMembers.push_back(idx);
                }
                ++jj;
            }

            if(clusterMembers.size() >= SV_MIN_SUPPORT) {
                for(size_t idx : clusterMembers) {
                    readGroups[idx].clusterId = nextClusterId;
                }
                ++nextClusterId;
            }

            // Advance past the window.
            ii = (jj > ii + 1) ? jj : ii + 1;
        }

        const int32_t totalClusters = nextClusterId;

        // Sort: deletions first (largest first), then inversions, then
        // insertions (largest first), then reference-like. Within each
        // SV type, larger SVs come first so they establish the graph
        // paths before smaller ones.
        sort(readGroups.begin(), readGroups.end(),
            [](const ReadGroup& a, const ReadGroup& b) {
                if(a.svType != b.svType)
                    return static_cast<int>(a.svType) < static_cast<int>(b.svType);
                // Within same type, larger SV first.
                return a.svSize > b.svSize;
            });

        // Now align reads in the sorted order.
        int seqId = 1;  // 0 is the backbone (reference).
        vector<string> seqNames;
        seqNames.push_back(string("ref_") + to_string(uint32_t(refId)));

        for(const auto& rg : readGroups) {
            const ReadId readId = rg.readId;
            const auto readNameSpan = readsRef.getReadName(readId);
            string readName(readNameSpan.data(), readNameSpan.size());

            // Tag with cluster ID if assigned.
            if(rg.clusterId >= 0) {
                readName += "_C" + to_string(rg.clusterId);
            }

            if(rg.chainIndicesInRef.empty()) continue;

            // Split chains by strand. Same-strand and opposite-strand
            // chains have ordinals on different marker arrays and must
            // be processed separately.
            // Strand groups: [0] = strand matching primary, [1] = opposite.
            const auto& firstCe = chainsForRef[rg.chainIndicesInRef[0]];
            const Strand primaryStrand = firstCe.isSameStrand ? 0 : 1;

            // Collect {boundaryIndex, readOrdinal, strand} per strand group.
            struct StrandGroup {
                Strand strand;
                vector<pair<uint32_t, uint32_t>> boundaryAndReadOrdinal;
            };
            StrandGroup strandGroups[2];
            strandGroups[0].strand = primaryStrand;
            strandGroups[1].strand = 1 - primaryStrand;

            for(size_t ci : rg.chainIndicesInRef) {
                const auto& ce = chainsForRef[ci];
                const auto& al = alignments[ce.chainIndex];
                const Strand thisStrand = ce.isSameStrand ? 0 : 1;
                int groupIdx = (thisStrand == primaryStrand) ? 0 : 1;

                for(const auto& ord : al.ordinals) {
                    auto it = ordinalToBoundary.find(ord[0]);
                    if(it != ordinalToBoundary.end()) {
                        strandGroups[groupIdx].boundaryAndReadOrdinal.push_back(
                            {it->second, ord[1]});
                    }
                }
            }

            bool anyAligned = false;

            // Process each strand group (primary strand first, then inverted).
            for(int gi = 0; gi < 2; ++gi) {
                auto& sg = strandGroups[gi];
                if(sg.boundaryAndReadOrdinal.size() < 2) continue;

                sort(sg.boundaryAndReadOrdinal.begin(), sg.boundaryAndReadOrdinal.end());
                sg.boundaryAndReadOrdinal.erase(
                    unique(sg.boundaryAndReadOrdinal.begin(), sg.boundaryAndReadOrdinal.end()),
                    sg.boundaryAndReadOrdinal.end());

                const OrientedReadId sgOid(readId, sg.strand);
                const auto sgMarkers = markersRef[sgOid.getValue()];

                const uint32_t sgBMin = sg.boundaryAndReadOrdinal.front().first;
                const uint32_t sgBMax = sg.boundaryAndReadOrdinal.back().first;
                if(sgBMax <= sgBMin) continue;
                if(sgBMin >= nodeIds.size()) continue;

                const string strandTag = (gi == 0) ? "" : "_INV";

                // Phase 1: Anchored segments for this strand group.
                for(size_t j = 0; j + 1 < sg.boundaryAndReadOrdinal.size(); ++j) {
                    const uint32_t bLeft = sg.boundaryAndReadOrdinal[j].first;
                    const uint32_t bRight = sg.boundaryAndReadOrdinal[j + 1].first;
                    const uint32_t readOrdLeft = sg.boundaryAndReadOrdinal[j].second;
                    const uint32_t readOrdRight = sg.boundaryAndReadOrdinal[j + 1].second;

                    if(bRight != bLeft + 1) continue;
                    if(readOrdRight <= readOrdLeft) continue;
                    if(readOrdLeft >= sgMarkers.size() || readOrdRight >= sgMarkers.size()) continue;

                    const uint32_t readPosLeft = sgMarkers[readOrdLeft].position + kHalf;
                    const uint32_t readPosRight = sgMarkers[readOrdRight].position + kHalf;
                    if(readPosRight <= readPosLeft) continue;

                    string readSeg;
                    readSeg.reserve(readPosRight - readPosLeft);
                    for(uint32_t pos = readPosLeft; pos < readPosRight; ++pos) {
                        readSeg.push_back(readsRef.getOrientedReadBase(sgOid, pos).character());
                    }

                    if(readSeg.empty()) continue;
                    if(bLeft >= nodeIds.size()) continue;
                    int endNode = (bRight < nodeIds.size())
                        ? static_cast<int>(nodeIds[bRight])
                        : -1;

                    seqNames.push_back(readName + "_seg" + to_string(j) + "_P1" + strandTag);

                    aligner.align_from(
                        readSeg,
                        nodeIds[bLeft],
                        1,      // weight
                        true,   // is_ends_free
                        0,      // start_offset
                        endNode,
                        seqId);

                    ++seqId;
                    anyAligned = true;
                    ++totalAlignedSegments;
                }

                // Phase 2: Full-span alignment for this strand group.
                const uint32_t readOrdMin = sg.boundaryAndReadOrdinal.front().second;
                const uint32_t readOrdMax = sg.boundaryAndReadOrdinal.back().second;
                if(readOrdMax > readOrdMin
                   && readOrdMin < sgMarkers.size()
                   && readOrdMax < sgMarkers.size()) {

                    const uint32_t readPosLeft = sgMarkers[readOrdMin].position + kHalf;
                    const uint32_t readPosRight = sgMarkers[readOrdMax].position + kHalf;

                    if(readPosRight > readPosLeft) {
                        string readSeq;
                        readSeq.reserve(readPosRight - readPosLeft);
                        for(uint32_t pos = readPosLeft; pos < readPosRight; ++pos) {
                            readSeq.push_back(readsRef.getOrientedReadBase(sgOid, pos).character());
                        }

                        if(!readSeq.empty()) {
                            int endNode = (sgBMax < nodeIds.size())
                                ? static_cast<int>(nodeIds[sgBMax])
                                : -1;

                            seqNames.push_back(readName + "_span_P2" + strandTag);

                            aligner.align_from(
                                readSeq,
                                nodeIds[sgBMin],
                                1,      // weight
                                true,   // is_ends_free
                                0,      // start_offset
                                endNode,
                                seqId);

                            ++seqId;
                            anyAligned = true;
                            ++totalAlignedSegments;
                        }
                    }
                }
            }

            if(anyAligned) {
                ++totalAlignedReads;
            }
        }

        // -----------------------------------------------------------------
        // Phase 3: Indirect alignment of insertion-internal reads via
        // a read graph BFS.
        //
        // Build a lightweight read graph from all-vs-all chains:
        //   Nodes = reads (short reads only, not reference)
        //   Edges = chains between reads (with shared ordinal pairs)
        //
        // For each placed read, build a map: readOrdinal -> boundaryIndex
        // (from its reference chains). Then BFS outward from placed reads.
        // When an unplaced read U is reached via a placed/already-placed
        // read P, map P's ordinals through the P→boundary map to find
        // U's backbone boundaries. This handles transitive cases: U1→U2→P.
        // -----------------------------------------------------------------
        {
            // Read graph: adjacency list.
            // Edge = {neighborReadId, chainIndex, iAmReadA}.
            struct ReadGraphEdge {
                uint32_t neighborReadId;
                uint64_t chainIndex;
                bool iAmReadA;  // true if this node is readIds[0] in the chain.
            };
            unordered_map<uint32_t, vector<ReadGraphEdge>> readGraph;

            // Build graph from all read-vs-read chains.
            for(uint64_t i = 0; i < n; ++i) {
                const auto& cand = candidates[i];
                const auto& al = alignments[i];
                if(al.ordinals.size() < 2) continue;

                // Skip secondaries.
                if(i < chainClassifications.size()
                   && chainClassifications[i] == ChainClassification::Secondary) {
                    continue;
                }

                const ReadId idA = cand.readIds[0];
                const ReadId idB = cand.readIds[1];

                // Skip read-vs-reference pairs.
                if(idA < ReadId(referenceReadCount)) continue;

                readGraph[uint32_t(idA)].push_back({uint32_t(idB), i, true});
                readGraph[uint32_t(idB)].push_back({uint32_t(idA), i, false});
            }

            // For each placed read, build ordinal -> boundary map.
            // Key: (readId, readOrdinal) -> boundaryIndex.
            // Use a flat map: readId -> (ordinal -> boundary).
            unordered_map<uint32_t,
                unordered_map<uint32_t, uint32_t>> readOrdinalToBoundary;

            for(const auto& ce : chainsForRef) {
                const auto& al = alignments[ce.chainIndex];
                auto& ordMap = readOrdinalToBoundary[uint32_t(ce.readId)];
                for(const auto& ord : al.ordinals) {
                    auto it = ordinalToBoundary.find(ord[0]);
                    if(it != ordinalToBoundary.end()) {
                        ordMap[ord[1]] = it->second;
                    }
                }
            }

            // Track placed reads.
            unordered_set<uint32_t> placedReadIds;
            for(const auto& rg : readGroups) {
                placedReadIds.insert(uint32_t(rg.readId));
            }

            // BFS from placed reads outward.
            queue<uint32_t> bfsQueue;
            for(uint32_t rid : placedReadIds) {
                bfsQueue.push(rid);
            }

            uint64_t indirectAligned = 0;

            while(!bfsQueue.empty()) {
                uint32_t currentId = bfsQueue.front();
                bfsQueue.pop();

                auto graphIt = readGraph.find(currentId);
                if(graphIt == readGraph.end()) continue;

                // currentId must have an ordinal->boundary map.
                auto mapIt = readOrdinalToBoundary.find(currentId);
                if(mapIt == readOrdinalToBoundary.end()) continue;
                const auto& currentOrdMap = mapIt->second;

                for(const auto& edge : graphIt->second) {
                    if(placedReadIds.count(edge.neighborReadId)) continue;

                    const auto& al = alignments[edge.chainIndex];
                    // ordinals: {ordA, ordB}.
                    // If edge.iAmReadA, currentId is A, neighbor is B.
                    const uint32_t myOrdIdx = edge.iAmReadA ? 0 : 1;
                    const uint32_t neighborOrdIdx = edge.iAmReadA ? 1 : 0;

                    // Map through: my ordinal -> boundary -> neighbor ordinal.
                    vector<pair<uint32_t, uint32_t>> neighborBoundaries;
                    for(const auto& ord : al.ordinals) {
                        auto bIt = currentOrdMap.find(ord[myOrdIdx]);
                        if(bIt != currentOrdMap.end()) {
                            neighborBoundaries.push_back({bIt->second, ord[neighborOrdIdx]});
                        }
                    }

                    if(neighborBoundaries.size() < 2) continue;

                    sort(neighborBoundaries.begin(), neighborBoundaries.end());
                    neighborBoundaries.erase(
                        unique(neighborBoundaries.begin(), neighborBoundaries.end()),
                        neighborBoundaries.end());

                    const uint32_t bMinN = neighborBoundaries.front().first;
                    const uint32_t bMaxN = neighborBoundaries.back().first;
                    if(bMaxN <= bMinN) continue;
                    if(bMinN >= nodeIds.size()) continue;

                    const uint32_t nOrdMin = neighborBoundaries.front().second;
                    const uint32_t nOrdMax = neighborBoundaries.back().second;
                    if(nOrdMax <= nOrdMin) continue;

                    // Determine strand from the chain's isSameStrand.
                    const auto& cand = candidates[edge.chainIndex];
                    // The neighbor's strand relative to the reference depends
                    // on the chain of edges. For simplicity, try strand 0 first.
                    Strand nStrand = 0;
                    auto nMarkers = markersRef[OrientedReadId(ReadId(edge.neighborReadId), 0).getValue()];
                    if(nOrdMin >= nMarkers.size() || nOrdMax >= nMarkers.size()) {
                        nStrand = 1;
                        nMarkers = markersRef[OrientedReadId(ReadId(edge.neighborReadId), 1).getValue()];
                        if(nOrdMin >= nMarkers.size() || nOrdMax >= nMarkers.size()) continue;
                    }

                    const OrientedReadId nOid(ReadId(edge.neighborReadId), nStrand);
                    const uint32_t nPosLeft = nMarkers[nOrdMin].position + kHalf;
                    const uint32_t nPosRight = nMarkers[nOrdMax].position + kHalf;
                    if(nPosRight <= nPosLeft) continue;

                    string nSeq;
                    nSeq.reserve(nPosRight - nPosLeft);
                    for(uint32_t pos = nPosLeft; pos < nPosRight; ++pos) {
                        nSeq.push_back(readsRef.getOrientedReadBase(nOid, pos).character());
                    }
                    if(nSeq.empty()) continue;

                    int endNodeN = (bMaxN < nodeIds.size())
                        ? static_cast<int>(nodeIds[bMaxN])
                        : -1;

                    const auto nNameSpan = readsRef.getReadName(ReadId(edge.neighborReadId));
                    seqNames.push_back(string(nNameSpan.data(), nNameSpan.size()));

                    aligner.align_from(
                        nSeq,
                        nodeIds[bMinN],
                        1,      // weight
                        true,   // is_ends_free
                        0,      // start_offset
                        endNodeN,
                        seqId);

                    ++seqId;
                    ++indirectAligned;
                    ++totalAlignedReads;
                    ++totalAlignedSegments;

                    // Mark as placed and propagate its ordinal->boundary map
                    // so further BFS steps can use it.
                    placedReadIds.insert(edge.neighborReadId);
                    auto& newOrdMap = readOrdinalToBoundary[edge.neighborReadId];
                    for(const auto& nb : neighborBoundaries) {
                        newOrdMap[nb.second] = nb.first;
                    }
                    bfsQueue.push(edge.neighborReadId);
                }
            }

            if(indirectAligned > 0) {
                cout << "    Indirectly aligned " << indirectAligned
                     << " insertion-internal reads via read graph BFS." << endl;
            }

            // ---------------------------------------------------------
            // Phase 4: Detect insertions from the read graph.
            //
            // Reads fully inside an insertion have no reference chains
            // but chain to other reads. We detect insertions by finding
            // connected components of unanchored reads that bridge
            // between left-anchored and right-anchored reads.
            //
            // For each edge (chain) between two reads, we compute the
            // "overhang": how far the neighbor read extends beyond the
            // overlap region, in base pairs. Walking a path from a
            // left-anchored read through unanchored reads to a
            // right-anchored read, the sum of overhangs gives the
            // total inserted sequence length.
            // ---------------------------------------------------------
            {
                // For each read, record its reference anchor range (bp).
                struct RefAnchorRange {
                    uint32_t minRefPos;
                    uint32_t maxRefPos;
                };
                unordered_map<uint32_t, RefAnchorRange> readRefRange;

                for(const auto& ce : chainsForRef) {
                    const auto& al = alignments[ce.chainIndex];
                    auto it = readRefRange.find(uint32_t(ce.readId));
                    uint32_t rMin = (it != readRefRange.end()) ? it->second.minRefPos : UINT32_MAX;
                    uint32_t rMax = (it != readRefRange.end()) ? it->second.maxRefPos : 0;
                    for(const auto& ord : al.ordinals) {
                        if(ord[0] < refMarkers.size()) {
                            const uint32_t pos = refMarkers[ord[0]].position;
                            rMin = std::min(rMin, pos);
                            rMax = std::max(rMax, pos);
                        }
                    }
                    readRefRange[uint32_t(ce.readId)] = {rMin, rMax};
                }

                // Identify unanchored reads (in the read graph but no ref chains).
                unordered_set<uint32_t> unanchoredReads;
                for(const auto& [rid, edges] : readGraph) {
                    if(readRefRange.find(rid) == readRefRange.end()) {
                        unanchoredReads.insert(rid);
                    }
                }

                // Helper: compute the overhang (bp) of the neighbor read
                // beyond the overlap region. Returns both left and right
                // overhangs so the caller can pick the correct direction.
                struct OverhangResult {
                    int64_t leftOverhang;   // bp before overlap start
                    int64_t rightOverhang;  // bp after overlap end
                    uint32_t overlapMinPos; // neighbor's overlap start (bp)
                    uint32_t overlapMaxPos; // neighbor's overlap end (bp)
                    uint32_t readMinPos;    // neighbor's first marker pos
                    uint32_t readMaxPos;    // neighbor's last marker pos
                };

                auto computeOverhangs = [&](const ReadGraphEdge& edge) -> OverhangResult {
                    OverhangResult result = {0, 0, 0, 0, 0, 0};
                    const auto& al = alignments[edge.chainIndex];
                    if(al.ordinals.size() < 2) return result;

                    const uint32_t neighborOrdIdx = edge.iAmReadA ? 1 : 0;

                    // Find the overlap region in the neighbor's ordinal space.
                    uint32_t overlapMinOrd = UINT32_MAX, overlapMaxOrd = 0;
                    for(const auto& ord : al.ordinals) {
                        overlapMinOrd = std::min(overlapMinOrd, ord[neighborOrdIdx]);
                        overlapMaxOrd = std::max(overlapMaxOrd, ord[neighborOrdIdx]);
                    }

                    // Get neighbor's markers (strand 0).
                    const auto nMarkers = markersRef[
                        OrientedReadId(ReadId(edge.neighborReadId), 0).getValue()];
                    const uint32_t nTotal = uint32_t(nMarkers.size());
                    if(nTotal == 0 || overlapMaxOrd >= nTotal || overlapMinOrd >= nTotal)
                        return result;

                    result.readMinPos = nMarkers[0].position;
                    result.readMaxPos = nMarkers[nTotal - 1].position;
                    result.overlapMinPos = nMarkers[overlapMinOrd].position;
                    result.overlapMaxPos = nMarkers[overlapMaxOrd].position;
                    result.rightOverhang = int64_t(result.readMaxPos) - int64_t(result.overlapMaxPos);
                    result.leftOverhang = int64_t(result.overlapMinPos) - int64_t(result.readMinPos);

                    return result;
                };

                if(!unanchoredReads.empty()) {
                    // Find connected components of unanchored reads.
                    unordered_map<uint32_t, int32_t> componentAssignment;
                    int32_t nextComponent = 0;

                    for(uint32_t seed : unanchoredReads) {
                        if(componentAssignment.count(seed)) continue;

                        // BFS to find the component.
                        queue<uint32_t> q;
                        q.push(seed);
                        componentAssignment[seed] = nextComponent;
                        vector<uint32_t> component;

                        while(!q.empty()) {
                            uint32_t cur = q.front();
                            q.pop();
                            component.push_back(cur);

                            auto gIt = readGraph.find(cur);
                            if(gIt == readGraph.end()) continue;
                            for(const auto& edge : gIt->second) {
                                if(unanchoredReads.count(edge.neighborReadId)
                                   && !componentAssignment.count(edge.neighborReadId)) {
                                    componentAssignment[edge.neighborReadId] = nextComponent;
                                    q.push(edge.neighborReadId);
                                }
                            }
                        }

                        // Collect all anchored reads that neighbor this component,
                        // with their ref positions and the connecting edge.
                        struct AnchoredNeighborInfo {
                            uint32_t readId;
                            uint32_t refMidPos;
                            uint32_t refMinPos;
                            uint32_t refMaxPos;
                            uint32_t connectingUnanchoredId;
                            uint64_t chainIndex;
                        };
                        vector<AnchoredNeighborInfo> anchoredNeighbors;

                        for(uint32_t uid : component) {
                            auto gIt = readGraph.find(uid);
                            if(gIt == readGraph.end()) continue;
                            for(const auto& edge : gIt->second) {
                                auto rrIt = readRefRange.find(edge.neighborReadId);
                                if(rrIt == readRefRange.end()) continue;
                                anchoredNeighbors.push_back({
                                    edge.neighborReadId,
                                    (rrIt->second.minRefPos + rrIt->second.maxRefPos) / 2,
                                    rrIt->second.minRefPos,
                                    rrIt->second.maxRefPos,
                                    uid,
                                    edge.chainIndex
                                });
                            }
                        }

                        if(anchoredNeighbors.size() < 2) {
                            ++nextComponent;
                            continue;
                        }

                        // Sort by ref position.
                        sort(anchoredNeighbors.begin(), anchoredNeighbors.end(),
                            [](const AnchoredNeighborInfo& a, const AnchoredNeighborInfo& b) {
                                return a.refMidPos < b.refMidPos;
                            });

                        // Deduplicate by readId.
                        {
                            vector<AnchoredNeighborInfo> deduped;
                            unordered_set<uint32_t> seen;
                            for(const auto& an : anchoredNeighbors) {
                                if(seen.insert(an.readId).second) {
                                    deduped.push_back(an);
                                }
                            }
                            anchoredNeighbors = std::move(deduped);
                        }

                        if(anchoredNeighbors.size() < 2) {
                            ++nextComponent;
                            continue;
                        }

                        // Left flank = anchored reads with smallest ref pos.
                        // Right flank = anchored reads with largest ref pos.
                        const uint32_t leftRefMax = anchoredNeighbors.front().refMaxPos;
                        const uint32_t rightRefMin = anchoredNeighbors.back().refMinPos;

                        // Ref gap between flanks.
                        const int64_t refGap = int64_t(rightRefMin) - int64_t(leftRefMax);

                        // Split into left and right sets.
                        const uint32_t splitPos = (leftRefMax + rightRefMin) / 2;
                        unordered_set<uint32_t> leftAnchored, rightAnchored;
                        for(const auto& an : anchoredNeighbors) {
                            if(an.refMidPos <= splitPos) leftAnchored.insert(an.readId);
                            else rightAnchored.insert(an.readId);
                        }

                        if(leftAnchored.empty() || rightAnchored.empty()) {
                            ++nextComponent;
                            continue;
                        }

                        // Estimate insertion size using coverage-based method.
                        //
                        // The unanchored reads tile the insertion. With
                        // coverage C, the insertion size ≈ total bases
                        // in unanchored reads / C. We estimate C from
                        // the overall dataset: totalBases / refLength.
                        //
                        // Also compute a tiling estimate: for each chain
                        // between reads in the component, the overlap
                        // region tells us how much sequence is shared.
                        // The unique contribution of each read is
                        // readLength - maxOverlap. Sum these for the
                        // tiling path length.

                        // Method 1: Coverage-based estimate.
                        uint64_t componentTotalBases = 0;
                        for(uint32_t uid : component) {
                            const auto uMarkers = markersRef[
                                OrientedReadId(ReadId(uid), 0).getValue()];
                            if(uMarkers.size() >= 2) {
                                componentTotalBases +=
                                    uMarkers[uMarkers.size() - 1].position
                                    - uMarkers[0].position;
                            }
                        }

                        // Estimate coverage from the dataset.
                        // Total bases / reference length.
                        const double estCoverage = std::max(1.0,
                            double(readsRef.getTotalBaseCount())
                            / double(refMarkers.empty() ? 1u :
                                uint32_t(refMarkers[refMarkers.size() - 1].position)));

                        const int64_t coverageEstimate =
                            int64_t(double(componentTotalBases) / estCoverage);

                        // Method 2: Tiling estimate.
                        // For each read in the component, find its maximum
                        // overlap with any other read in the component
                        // (from chain coordinates). The unique extension
                        // is readSpan - maxOverlap.
                        int64_t tilingEstimate = 0;
                        for(uint32_t uid : component) {
                            const auto uMarkers = markersRef[
                                OrientedReadId(ReadId(uid), 0).getValue()];
                            if(uMarkers.size() < 2) continue;
                            const int64_t readSpan =
                                int64_t(uMarkers[uMarkers.size() - 1].position)
                                - int64_t(uMarkers[0].position);

                            // Find max overlap with any neighbor in the component.
                            int64_t maxOverlap = 0;
                            auto gIt = readGraph.find(uid);
                            if(gIt != readGraph.end()) {
                                for(const auto& edge : gIt->second) {
                                    // Only consider edges within the component
                                    // or to anchored neighbors.
                                    bool inComponent = (componentAssignment.count(edge.neighborReadId)
                                        && componentAssignment[edge.neighborReadId] == nextComponent);
                                    if(!inComponent) continue;

                                    auto oh = computeOverhangs(edge);
                                    const int64_t overlapSpan =
                                        int64_t(oh.overlapMaxPos) - int64_t(oh.overlapMinPos);
                                    maxOverlap = std::max(maxOverlap, overlapSpan);
                                }
                            }

                            const int64_t uniqueExtension = readSpan - maxOverlap;
                            if(uniqueExtension > 0) {
                                tilingEstimate += uniqueExtension;
                            }
                        }

                        // Use the average of both estimates, but prefer
                        // the coverage estimate for larger components
                        // (more statistically robust).
                        const int64_t bestInsertionSize = (component.size() >= 5)
                            ? coverageEstimate
                            : std::max(coverageEstimate, tilingEstimate);
                        const uint32_t bestBreakpoint = leftRefMax;
                        const int bestPathReads = int(component.size());

                        // Filter: for a real insertion, the left and right
                        // anchored reads should be close on the reference
                        // (the insertion adds sequence at a single point).
                        // Allow refGap up to ~500bp to account for marker
                        // resolution and short-read mapping uncertainty.
                        // Filters for a real insertion:
                        // 1. Left and right anchored reads should be close
                        //    on the reference (insertion at a single point).
                        // 2. Component should have enough reads for a
                        //    reliable estimate (at least 3).
                        // 3. Insertion size should be positive and meaningful.
                        const bool validBreakpoint = (std::abs(refGap) < 500);
                        const bool enoughReads = (component.size() >= 3);

                        if(bestInsertionSize > 20 && validBreakpoint && enoughReads) {
                            cout << "    Insertion detected via read graph: "
                                 << "size=" << bestInsertionSize << "bp, "
                                 << "breakpoint=" << bestBreakpoint << ", "
                                 << "component=" << component.size() << " reads, "
                                 << "path=" << bestPathReads << " intermediate reads, "
                                 << "refGap=" << refGap << "bp."
                                 << endl;
                        }

                        ++nextComponent;
                    }

                    if(nextComponent > 0) {
                        cout << "    Read-graph insertion scan: "
                             << unanchoredReads.size() << " unanchored reads in "
                             << nextComponent << " components." << endl;
                    }
                }
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

        // -----------------------------------------------------------------
        // Step 7: Output per-cluster SV summary.
        // -----------------------------------------------------------------
        if(totalClusters > 0) {
            const string clusterFileName = outputPrefix + "_ref"
                + to_string(uint32_t(refId)) + ".sv_clusters.tsv";
            ofstream clusterOut(clusterFileName);
            if(clusterOut) {
                clusterOut << "cluster_id\tsv_type\tnum_reads\tbreakpoint_pos\t"
                           << "mean_sv_size\tmin_sv_size\tmax_sv_size\n";

                for(int32_t cid = 0; cid < totalClusters; ++cid) {
                    uint32_t count = 0;
                    int64_t sumSize = 0;
                    int64_t minSize = INT64_MAX;
                    int64_t maxSize = INT64_MIN;
                    uint64_t sumPos = 0;
                    SvType cType = SvType::ReferenceLike;

                    for(const auto& rg : readGroups) {
                        if(rg.clusterId != cid) continue;
                        ++count;
                        sumSize += rg.svSize;
                        sumPos += rg.breakpointRefPos;
                        minSize = std::min(minSize, rg.svSize);
                        maxSize = std::max(maxSize, rg.svSize);
                        cType = rg.svType;
                    }

                    if(count == 0) continue;

                    const char* typeStr = "UNKNOWN";
                    switch(cType) {
                        case SvType::Deletion:  typeStr = "DEL"; break;
                        case SvType::Inversion: typeStr = "INV"; break;
                        case SvType::Insertion: typeStr = "INS"; break;
                        case SvType::ReferenceLike: typeStr = "REF"; break;
                    }

                    clusterOut << cid << "\t"
                               << typeStr << "\t"
                               << count << "\t"
                               << (sumPos / count) << "\t"
                               << (sumSize / int64_t(count)) << "\t"
                               << minSize << "\t"
                               << maxSize << "\n";
                }
            }
            cout << "    SV clusters: " << totalClusters
                 << " (output: " << clusterFileName << ")" << endl;
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
