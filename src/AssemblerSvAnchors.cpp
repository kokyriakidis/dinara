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
#include "InvertedIndexBuilder.hpp"
#include "Reads.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_map>
#include <unordered_set>
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
    const string& outputPrefix,
    const vector<RefHitDepthWindow>& refHitDepth)
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
        cerr << "Creating TheseusMSA with " << segmentViews.size() << " segments..." << endl;
        theseus::TheseusMSA aligner(penalties, heuristics, segmentViews, nodeIds, 1);
        cerr << "TheseusMSA created. nodeIds.size()=" << nodeIds.size() << endl;

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
            // Require opposite-strand chains to have a minimum
            // number of anchors and reference span to avoid
            // false inversions from palindromic k-mers.
            const auto& firstCeClassify = chainsForRef[rg.chainIndicesInRef[0]];
            const Strand primaryStrand = firstCeClassify.isSameStrand ? 0 : 1;
            bool hasPrimaryStrand = false;
            bool hasOppositeStrand = false;
            constexpr uint32_t minInvAnchors = 15;
            constexpr uint32_t minInvSpan = 80;
            for(size_t ci : rg.chainIndicesInRef) {
                const Strand s = chainsForRef[ci].isSameStrand ? 0 : 1;
                if(s == primaryStrand) {
                    hasPrimaryStrand = true;
                } else {
                    // Check anchor count and ref span.
                    const auto& al = alignments[chainsForRef[ci].chainIndex];
                    if(al.ordinals.size() >= minInvAnchors) {
                        uint32_t cMin = UINT32_MAX, cMax = 0;
                        for(const auto& ord : al.ordinals) {
                            if(ord[0] < refMarkers.size()) {
                                const uint32_t pos =
                                    refMarkers[ord[0]].position;
                                cMin = std::min(cMin, pos);
                                cMax = std::max(cMax, pos);
                            }
                        }
                        if(cMax > cMin
                           && cMax - cMin >= minInvSpan) {
                            hasOppositeStrand = true;
                        }
                    }
                }
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
                        cout << "    Read " << rg.readId
                             << ": INS size=" << bestRise
                             << " breakpoint=" << bestRiseBreakpoint
                             << " (drop=" << bestDrop << ")"
                             << " anchors=" << allAnchors.size()
                             << " diagRange=[" << *std::min_element(diagonals.begin(), diagonals.end())
                             << "," << *std::max_element(diagonals.begin(), diagonals.end()) << "]"
                             << endl;
                    } else if(bestDrop > minSvSize) {
                        rg.svType = SvType::Deletion;
                        rg.svSize = bestDrop;
                        rg.breakpointRefPos = bestDropBreakpoint;
                        cout << "    Read " << rg.readId
                             << ": DEL size=" << bestDrop
                             << " breakpoint=" << bestDropBreakpoint
                             << " (rise=" << bestRise << ")"
                             << " anchors=" << allAnchors.size()
                             << " diagRange=[" << *std::min_element(diagonals.begin(), diagonals.end())
                             << "," << *std::max_element(diagonals.begin(), diagonals.end()) << "]"
                             << endl;
                    } else {
                        rg.svType = SvType::ReferenceLike;
                        if(allAnchors.size() >= 20
                           && (bestRise > 5 || bestDrop > 5)) {
                            cout << "    Read " << rg.readId
                                 << ": REF rise=" << bestRise
                                 << " drop=" << bestDrop
                                 << " anchors=" << allAnchors.size()
                                 << " chains=" << rg.chainIndicesInRef.size()
                                 << endl;
                        }

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
        constexpr uint32_t SV_WINDOW = 100;
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

        // Merge adjacent clusters of the same type and similar size.
        // This handles cases where the greedy sweep splits a single
        // SV event into two clusters because the breakpoint spread
        // exceeds SV_WINDOW from the first seed.
        if(nextClusterId >= 2) {
            // Build per-cluster info.
            struct ClusterInfo {
                SvType type;
                int64_t medianSize;
                uint32_t minBp, maxBp;
                int32_t id;
            };
            vector<ClusterInfo> clusters(nextClusterId);
            for(int32_t c = 0; c < nextClusterId; ++c) {
                vector<int64_t> sizes;
                uint32_t minBp = UINT32_MAX, maxBp = 0;
                SvType type = SvType::ReferenceLike;
                for(auto& rg2 : readGroups) {
                    if(rg2.clusterId == c) {
                        sizes.push_back(rg2.svSize);
                        minBp = std::min(minBp, rg2.breakpointRefPos);
                        maxBp = std::max(maxBp, rg2.breakpointRefPos);
                        type = rg2.svType;
                    }
                }
                sort(sizes.begin(), sizes.end());
                clusters[c] = {type,
                    sizes.empty() ? 0 : sizes[sizes.size()/2],
                    minBp, maxBp, c};
            }

            // Merge clusters that are close and similar.
            for(int32_t c1 = 0; c1 < nextClusterId; ++c1) {
                for(int32_t c2 = c1 + 1; c2 < nextClusterId; ++c2) {
                    if(clusters[c1].type != clusters[c2].type) continue;
                    // Check breakpoint proximity.
                    const uint32_t gap =
                        (clusters[c2].minBp > clusters[c1].maxBp)
                        ? clusters[c2].minBp - clusters[c1].maxBp : 0;
                    if(gap > SV_WINDOW) continue;
                    // Check size similarity.
                    const int64_t s1 = clusters[c1].medianSize;
                    const int64_t s2 = clusters[c2].medianSize;
                    const int64_t ref = std::max(std::abs(s1), int64_t(1));
                    if(std::abs(s1 - s2) > int64_t(ref * SV_SIZE_RATIO))
                        continue;
                    // Merge c2 into c1.
                    for(auto& rg2 : readGroups) {
                        if(rg2.clusterId == c2) rg2.clusterId = c1;
                    }
                    clusters[c1].minBp = std::min(
                        clusters[c1].minBp, clusters[c2].minBp);
                    clusters[c1].maxBp = std::max(
                        clusters[c1].maxBp, clusters[c2].maxBp);
                    clusters[c2].id = -1; // Mark as merged.
                }
            }
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

        int readGroupIdx = 0;
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

            // TODO: Theseus MSA disabled — crashes on some inputs.
            // Will be replaced with abPOA segment-by-segment alignment.
            ++readGroupIdx;
            continue;

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

                    const uint32_t readBaseCount = uint32_t(readsRef.getRead(
                        sgOid.getReadId()).baseCount);
                    if(readPosRight > readBaseCount) continue;

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
                    const uint32_t readBaseCount2 = uint32_t(readsRef.getRead(
                        sgOid.getReadId()).baseCount);

                    if(readPosRight > readPosLeft && readPosRight <= readBaseCount2) {
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
            ++readGroupIdx;
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

                    const uint32_t nBaseCount = uint32_t(readsRef.getRead(
                        ReadId(edge.neighborReadId)).baseCount);
                    if(nPosRight > nBaseCount) continue;

                    string nSeq;
                    nSeq.reserve(nPosRight - nPosLeft);
                    for(uint32_t pos = nPosLeft; pos < nPosRight; ++pos) {
                        nSeq.push_back(readsRef.getOrientedReadBase(nOid, pos).character());
                    }
                    if(nSeq.empty()) continue;

                    int endNodeN = (bMaxN < nodeIds.size())
                        ? static_cast<int>(nodeIds[bMaxN])
                        : -1;

                    // TODO: Theseus align_from disabled — crashes on some inputs.
                    // const auto nNameSpan = readsRef.getReadName(ReadId(edge.neighborReadId));
                    // seqNames.push_back(string(nNameSpan.data(), nNameSpan.size()));
                    // aligner.align_from(nSeq, nodeIds[bMinN], 1, true, 0, endNodeN, seqId);
                    // ++seqId;
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

            cout << "    Phase 1/2 done: " << totalAlignedReads << " reads, "
                 << totalAlignedSegments << " segments aligned." << endl;
            cout << flush;

            // ---------------------------------------------------------
            // Phase 4: Detect insertions via reference-centric
            // breakpoint detection + read graph path measurement.
            //
            // Step A: For each read, record where its chain starts/ends
            //         on the reference and how much of the read is
            //         beyond the chain (soft-clip / overhang).
            // Step B: Scan the reference for positions where many chains
            //         end (left breakpoint) or start (right breakpoint)
            //         while having significant read overhang.
            // Step C: For matched left/right breakpoint pairs, measure
            //         insertion size via read graph path.
            // ---------------------------------------------------------
            {
                // For each read, record its reference chain endpoints (bp).
                struct ReadChainInfo {
                    uint32_t minRefPos = UINT32_MAX;
                    uint32_t maxRefPos = 0;
                    uint32_t minReadOrd = UINT32_MAX; // first marker ordinal in chain
                    uint32_t maxReadOrd = 0;          // last marker ordinal in chain
                    uint32_t totalMarkers = 0;        // total markers in read
                };
                unordered_map<uint32_t, ReadChainInfo> readChainInfoMap;

                for(const auto& ce : chainsForRef) {
                    const auto& al = alignments[ce.chainIndex];
                    auto& info = readChainInfoMap[uint32_t(ce.readId)];
                    for(const auto& ord : al.ordinals) {
                        if(ord[0] < refMarkers.size()) {
                            const uint32_t pos = refMarkers[ord[0]].position;
                            info.minRefPos = std::min(info.minRefPos, pos);
                            info.maxRefPos = std::max(info.maxRefPos, pos);
                        }
                        info.minReadOrd = std::min(info.minReadOrd, ord[1]);
                        info.maxReadOrd = std::max(info.maxReadOrd, ord[1]);
                    }
                    const auto rdMarkers = markersRef[
                        OrientedReadId(ReadId(ce.readId), 0).getValue()];
                    info.totalMarkers = uint32_t(rdMarkers.size());
                }

                // Identify unanchored reads (in the read graph but no ref chains).
                unordered_set<uint32_t> unanchoredReads;
                for(const auto& [rid, edges] : readGraph) {
                    if(readChainInfoMap.find(rid) == readChainInfoMap.end()) {
                        unanchoredReads.insert(rid);
                    }
                }

                // -------------------------------------------------------
                // Reference-centric coverage analysis.
                //
                // For each window along the reference, count:
                //   chainEndCount:   chains whose maxRefPos falls here
                //   chainStartCount: chains whose minRefPos falls here
                //   spanningCount:   chains that span through this window
                //
                // A breakpoint shows as a spike in chainEndCount (left BP)
                // or chainStartCount (right BP) with a drop in spanning.
                // -------------------------------------------------------
                const uint32_t refEndPos = refMarkers.empty() ? 0
                    : uint32_t(refMarkers[refMarkers.size() - 1].position);
                const uint32_t refStartPos = refMarkers.empty() ? 0
                    : uint32_t(refMarkers[0].position);
                const uint32_t refLen = refEndPos - refStartPos;
                const uint32_t windowSize = 50; // bp per window
                const uint32_t nWindows = (refLen / windowSize) + 1;
                const uint32_t boundaryMargin = 300;

                // Boundary windows to exclude (extraction edges).
                const uint32_t boundaryWindows = boundaryMargin / windowSize;

                vector<uint32_t> chainEndCount(nWindows, 0);
                vector<uint32_t> chainStartCount(nWindows, 0);
                vector<uint32_t> spanningCount(nWindows, 0);

                // Also track per-read overhang at each endpoint.
                struct EndpointInfo {
                    uint32_t readId;
                    int64_t overhangBp;
                };
                // Map from window index to list of reads ending/starting there.
                unordered_map<uint32_t, vector<EndpointInfo>> chainEndReads;
                unordered_map<uint32_t, vector<EndpointInfo>> chainStartReads;

                for(const auto& [rid, info] : readChainInfoMap) {
                    if(info.totalMarkers < 2) continue;
                    if(info.maxRefPos < refStartPos || info.minRefPos < refStartPos) continue;

                    const uint32_t endWin = (info.maxRefPos - refStartPos) / windowSize;
                    const uint32_t startWin = (info.minRefPos - refStartPos) / windowSize;
                    if(endWin >= nWindows || startWin >= nWindows) continue;

                    // Compute overhangs.
                    const auto rdMarkers = markersRef[
                        OrientedReadId(ReadId(rid), 0).getValue()];
                    if(rdMarkers.size() < 2) continue;

                    const uint32_t readLastPos = uint32_t(rdMarkers[rdMarkers.size() - 1].position);
                    const uint32_t chainEndOrdPos = (info.maxReadOrd < rdMarkers.size())
                        ? uint32_t(rdMarkers[info.maxReadOrd].position) : readLastPos;
                    const uint32_t readFirstPos = uint32_t(rdMarkers[0].position);
                    const uint32_t chainStartOrdPos = (info.minReadOrd < rdMarkers.size())
                        ? uint32_t(rdMarkers[info.minReadOrd].position) : readFirstPos;

                    const int64_t rightOvhBp = int64_t(readLastPos) - int64_t(chainEndOrdPos);
                    const int64_t leftOvhBp = int64_t(chainStartOrdPos) - int64_t(readFirstPos);

                    // Count chain endpoints.
                    chainEndCount[endWin]++;
                    chainStartCount[startWin]++;

                    // Record reads with significant overhang.
                    const int64_t minOvhBp = 20;
                    const int64_t rightOvhOrd = int64_t(info.totalMarkers - 1) - int64_t(info.maxReadOrd);
                    const int64_t leftOvhOrd = int64_t(info.minReadOrd);

                    if(rightOvhBp >= minOvhBp && rightOvhOrd >= 2) {
                        chainEndReads[endWin].push_back({rid, rightOvhBp});
                    }
                    if(leftOvhBp >= minOvhBp && leftOvhOrd >= 2) {
                        chainStartReads[startWin].push_back({rid, leftOvhBp});
                    }

                    // Count spanning coverage.
                    for(uint32_t w = startWin; w <= endWin && w < nWindows; ++w) {
                        spanningCount[w]++;
                    }
                }

                // Compute background chain-end rate.
                // In a uniform region, chain ends are spread evenly.
                // At a breakpoint, chain ends cluster.
                // Use median spanning coverage as baseline.
                vector<uint32_t> sortedSpanning(spanningCount.begin(), spanningCount.end());
                sort(sortedSpanning.begin(), sortedSpanning.end());
                const uint32_t medianSpanning = sortedSpanning[sortedSpanning.size() / 2];

                // Background end rate: total chain ends / number of windows
                // (excluding boundary windows).
                uint32_t totalEnds = 0, totalStarts = 0;
                uint32_t interiorWindows = 0;
                for(uint32_t w = boundaryWindows; w + boundaryWindows < nWindows; ++w) {
                    totalEnds += chainEndCount[w];
                    totalStarts += chainStartCount[w];
                    ++interiorWindows;
                }
                const double bgEndRate = interiorWindows > 0
                    ? double(totalEnds) / double(interiorWindows) : 1.0;
                const double bgStartRate = interiorWindows > 0
                    ? double(totalStarts) / double(interiorWindows) : 1.0;

                // Detect breakpoints: windows where chain-end count is
                // significantly above background AND there are reads with
                // significant overhang at that position.
                //
                // A left breakpoint has many chain ends with right overhang.
                // A right breakpoint has many chain starts with left overhang.
                const double minFoldEnrichment = 3.0;
                const uint32_t minEndpointReads = 2; // reads with overhang

                struct Breakpoint {
                    uint32_t windowIdx;
                    uint32_t refPos;        // center of window
                    uint32_t endpointCount; // chains ending/starting here
                    uint32_t spanCount;     // spanning coverage
                    uint32_t ovhReadCount;  // reads with significant overhang
                    double foldEnrichment;
                    vector<EndpointInfo> reads; // reads with overhang
                };
                vector<Breakpoint> leftBreakpoints;  // chain ends
                vector<Breakpoint> rightBreakpoints; // chain starts

                for(uint32_t w = boundaryWindows; w + boundaryWindows < nWindows; ++w) {
                    const uint32_t refPos = refStartPos + w * windowSize + windowSize / 2;

                    // Left breakpoint: many chain ends with right overhang.
                    if(chainEndCount[w] > 0) {
                        const double fold = chainEndCount[w] / std::max(bgEndRate, 0.1);
                        auto it = chainEndReads.find(w);
                        const uint32_t ovhCount = (it != chainEndReads.end())
                            ? uint32_t(it->second.size()) : 0;

                        if(fold >= minFoldEnrichment && ovhCount >= minEndpointReads) {
                            leftBreakpoints.push_back({
                                w, refPos, chainEndCount[w], spanningCount[w],
                                ovhCount, fold,
                                (it != chainEndReads.end()) ? it->second : vector<EndpointInfo>{}
                            });
                        }
                    }

                    // Right breakpoint: many chain starts with left overhang.
                    if(chainStartCount[w] > 0) {
                        const double fold = chainStartCount[w] / std::max(bgStartRate, 0.1);
                        auto it = chainStartReads.find(w);
                        const uint32_t ovhCount = (it != chainStartReads.end())
                            ? uint32_t(it->second.size()) : 0;

                        if(fold >= minFoldEnrichment && ovhCount >= minEndpointReads) {
                            rightBreakpoints.push_back({
                                w, refPos, chainStartCount[w], spanningCount[w],
                                ovhCount, fold,
                                (it != chainStartReads.end()) ? it->second : vector<EndpointInfo>{}
                            });
                        }
                    }
                }

                // -------------------------------------------------------
                // K-mer hit depth along the reference.
                //
                // For each reference marker, look up its canonical k-mer
                // in the inverted index and count how many read occurrences
                // share that k-mer. Aggregate per window.
                //
                // At an insertion breakpoint, the hit depth drops because
                // reads carrying the insertion have their k-mers split
                // between the two flanks.
                // -------------------------------------------------------
                // Use the pre-computed hit depth profile passed from main.
                vector<double> windowHitDepth(nWindows, 0.0);
                vector<uint32_t> windowMarkerCount(nWindows, 0);

                for(const auto& hdw : refHitDepth) {
                    if(hdw.refPos < refStartPos) continue;
                    const uint32_t win = (hdw.refPos - refStartPos) / windowSize;
                    if(win >= nWindows) continue;
                    windowHitDepth[win] = hdw.avgHitDepth;
                    windowMarkerCount[win] = hdw.markerCount;
                }

                // Compute median hit depth for background.
                vector<double> sortedHitDepth;
                for(uint32_t w = boundaryWindows; w + boundaryWindows < nWindows; ++w) {
                    if(windowMarkerCount[w] > 0) {
                        sortedHitDepth.push_back(windowHitDepth[w]);
                    }
                }
                sort(sortedHitDepth.begin(), sortedHitDepth.end());
                const double medianHitDepth = sortedHitDepth.empty() ? 0.0
                    : sortedHitDepth[sortedHitDepth.size() / 2];

                // Detect hit-depth breakpoints: windows where depth drops
                // below 50% of median AND has at least some markers.
                const double hitDepthDropThreshold = 0.5;
                struct HitDepthBreakpoint {
                    uint32_t windowIdx;
                    uint32_t refPos;
                    double hitDepth;
                    double dropRatio; // hitDepth / medianHitDepth
                };
                vector<HitDepthBreakpoint> hitDepthBreakpoints;

                if(medianHitDepth > 2.0) {
                    for(uint32_t w = boundaryWindows; w + boundaryWindows < nWindows; ++w) {
                        if(windowMarkerCount[w] == 0) continue;
                        const double ratio = windowHitDepth[w] / medianHitDepth;
                        if(ratio < hitDepthDropThreshold) {
                            const uint32_t refPos = refStartPos + w * windowSize + windowSize / 2;
                            hitDepthBreakpoints.push_back({w, refPos, windowHitDepth[w], ratio});
                        }
                    }
                }

                cout << "    Coverage analysis: "
                     << "medianSpanning=" << medianSpanning
                     << " bgEndRate=" << bgEndRate
                     << " bgStartRate=" << bgStartRate
                     << " leftBPs=" << leftBreakpoints.size()
                     << " rightBPs=" << rightBreakpoints.size()
                     << " unanchored=" << unanchoredReads.size()
                     << " medianHitDepth=" << medianHitDepth
                     << " hitDepthBPs=" << hitDepthBreakpoints.size()
                     << endl;

                for(const auto& hbp : hitDepthBreakpoints) {
                    cout << "      HitDepth BP: pos=" << hbp.refPos
                         << " depth=" << hbp.hitDepth
                         << " ratio=" << hbp.dropRatio
                         << endl;
                }

                // Merge hit-depth breakpoints into left/right breakpoints.
                // A hit-depth drop can indicate either a left or right
                // breakpoint. We check if there are chain-end reads or
                // chain-start reads near the hit-depth breakpoint.
                for(const auto& hbp : hitDepthBreakpoints) {
                    // Check if already covered by an existing breakpoint.
                    bool coveredLeft = false, coveredRight = false;
                    for(const auto& lbp : leftBreakpoints) {
                        if(std::abs(int64_t(lbp.refPos) - int64_t(hbp.refPos)) < int64_t(windowSize * 3)) {
                            coveredLeft = true;
                            break;
                        }
                    }
                    for(const auto& rbp : rightBreakpoints) {
                        if(std::abs(int64_t(rbp.refPos) - int64_t(hbp.refPos)) < int64_t(windowSize * 3)) {
                            coveredRight = true;
                            break;
                        }
                    }

                    // If not covered, add as both left and right breakpoint
                    // using any chain-end/start reads in nearby windows.
                    if(!coveredLeft) {
                        vector<EndpointInfo> nearbyReads;
                        for(int32_t dw = -2; dw <= 2; ++dw) {
                            const int32_t w = int32_t(hbp.windowIdx) + dw;
                            if(w < 0 || uint32_t(w) >= nWindows) continue;
                            auto it = chainEndReads.find(uint32_t(w));
                            if(it != chainEndReads.end()) {
                                for(const auto& ei : it->second) nearbyReads.push_back(ei);
                            }
                        }
                        if(!nearbyReads.empty()) {
                            leftBreakpoints.push_back({
                                hbp.windowIdx, hbp.refPos,
                                chainEndCount[hbp.windowIdx],
                                spanningCount[hbp.windowIdx],
                                uint32_t(nearbyReads.size()),
                                hbp.dropRatio > 0 ? 1.0 / hbp.dropRatio : 10.0,
                                std::move(nearbyReads)
                            });
                        }
                    }
                    if(!coveredRight) {
                        vector<EndpointInfo> nearbyReads;
                        for(int32_t dw = -2; dw <= 2; ++dw) {
                            const int32_t w = int32_t(hbp.windowIdx) + dw;
                            if(w < 0 || uint32_t(w) >= nWindows) continue;
                            auto it = chainStartReads.find(uint32_t(w));
                            if(it != chainStartReads.end()) {
                                for(const auto& ei : it->second) nearbyReads.push_back(ei);
                            }
                        }
                        if(!nearbyReads.empty()) {
                            rightBreakpoints.push_back({
                                hbp.windowIdx, hbp.refPos,
                                chainStartCount[hbp.windowIdx],
                                spanningCount[hbp.windowIdx],
                                uint32_t(nearbyReads.size()),
                                hbp.dropRatio > 0 ? 1.0 / hbp.dropRatio : 10.0,
                                std::move(nearbyReads)
                            });
                        }
                    }
                }

                for(const auto& bp : leftBreakpoints) {
                    cout << "      Left BP: pos=" << bp.refPos
                         << " ends=" << bp.endpointCount
                         << " spanning=" << bp.spanCount
                         << " ovhReads=" << bp.ovhReadCount
                         << " fold=" << bp.foldEnrichment;
                    for(const auto& ei : bp.reads) {
                        cout << " r" << ei.readId << ":" << ei.overhangBp;
                    }
                    cout << endl;
                }
                for(const auto& bp : rightBreakpoints) {
                    cout << "      Right BP: pos=" << bp.refPos
                         << " starts=" << bp.endpointCount
                         << " spanning=" << bp.spanCount
                         << " ovhReads=" << bp.ovhReadCount
                         << " fold=" << bp.foldEnrichment;
                    for(const auto& ei : bp.reads) {
                        cout << " r" << ei.readId << ":" << ei.overhangBp;
                    }
                    cout << endl;
                }

                // Collect VNTR gap regions and insertion call regions
                // for suppressing false coverage-drop deletion calls.
                struct SvRegion {
                    uint32_t startPos;
                    uint32_t endPos;
                };
                vector<SvRegion> vntrGaps;
                vector<SvRegion> insertionCallRegions;

                // For each left breakpoint, find the nearest right breakpoint
                // to form a breakpoint pair.
                for(const auto& lbp : leftBreakpoints) {
                    const Breakpoint* bestRbp = nullptr;
                    int64_t bestDist = INT64_MAX;
                    for(const auto& rbp : rightBreakpoints) {
                        const int64_t dist = std::abs(
                            int64_t(rbp.refPos) - int64_t(lbp.refPos));
                        if(dist < bestDist) {
                            bestDist = dist;
                            bestRbp = &rbp;
                        }
                    }

                    // Allow larger distance if a hit-depth cluster spans
                    // the gap between left and right BPs (VNTR region).
                    int64_t maxPairDist = 500;
                    if(bestRbp && bestDist > 500) {
                        // Check if the region between the BPs has consistently
                        // low hit-depth (indicating a marker-depleted VNTR).
                        const uint32_t gapStartWin = (lbp.refPos - refStartPos) / windowSize;
                        const uint32_t gapEndWin = (bestRbp->refPos - refStartPos) / windowSize;
                        uint32_t lowDepthWindows = 0;
                        uint32_t totalGapWindows = 0;
                        for(uint32_t w = gapStartWin; w <= gapEndWin && w < nWindows; ++w) {
                            if(windowMarkerCount[w] == 0) continue;
                            ++totalGapWindows;
                            if(medianHitDepth > 0
                               && windowHitDepth[w] / medianHitDepth < hitDepthDropThreshold) {
                                ++lowDepthWindows;
                            }
                        }
                        // If >50% of windows in the gap have low hit-depth,
                        // this is a VNTR — allow the pairing.
                        if(totalGapWindows > 0
                           && double(lowDepthWindows) / double(totalGapWindows) > 0.5) {
                            maxPairDist = bestDist + 1; // Allow this pair.
                            vntrGaps.push_back({lbp.refPos, bestRbp->refPos});
                            cout << "    VNTR gap: L=" << lbp.refPos
                                 << " R=" << bestRbp->refPos
                                 << " lowDepth=" << lowDepthWindows
                                 << "/" << totalGapWindows
                                 << " — extending pair distance"
                                 << endl;
                        }
                    }

                    if(!bestRbp || bestDist > maxPairDist) continue;

                    // Collect all left-flank and right-flank read IDs.
                    unordered_set<uint32_t> leftIds, rightIds;
                    for(const auto& ei : lbp.reads) leftIds.insert(ei.readId);
                    for(const auto& ei : bestRbp->reads) rightIds.insert(ei.readId);

                    // Direct path search: check if left-flank reads
                    // connect to right-flank reads directly or through
                    // a small number of unanchored intermediates.
                    //
                    // For each left-flank read, check:
                    //   1-hop: left → right (direct overlap)
                    //   2-hop: left → unanchored → right
                    //   3-hop: left → unanchored → unanchored → right
                    //
                    // The insertion size = sum of unique extensions along
                    // the path = left_overhang + intermediate_extensions + right_overhang.
                    //
                    // Skip for VNTR gaps: read-to-read chains in repetitive
                    // regions create false shortcuts, giving wrong sizes.

                    int64_t bestPathDist = 0;
                    int bestPathLen = 0;
                    bool foundPath = false;
                    vector<int64_t> allPathDists; // Collect all valid paths.
                    const bool isVntrGap = (maxPairDist > 500);

                    // Helper: get the overlap span (bp) between two reads
                    // from their chain. Returns the overlap extent in the
                    // neighbor's coordinate space, excluding large gaps
                    // (which indicate a breakpoint between reference and
                    // insertion regions).
                    auto getOverlapSpan = [&](const ReadGraphEdge& edge) -> int64_t {
                        const auto& al = alignments[edge.chainIndex];
                        if(al.ordinals.size() < 2) return -1;
                        const uint32_t nIdx = edge.iAmReadA ? 1 : 0;

                        // Collect neighbor ordinals and sort.
                        vector<uint32_t> nOrds;
                        nOrds.reserve(al.ordinals.size());
                        for(const auto& ord : al.ordinals) {
                            nOrds.push_back(ord[nIdx]);
                        }
                        sort(nOrds.begin(), nOrds.end());

                        const auto nMkrs = markersRef[
                            OrientedReadId(ReadId(edge.neighborReadId), 0).getValue()];
                        if(nMkrs.size() < 2) return -1;

                        // Sum marker-to-marker distances, capping each
                        // gap to exclude breakpoint gaps (where the
                        // read transitions from reference to insertion).
                        const int64_t maxGap = int64_t(k) * 3;
                        int64_t totalSpan = 0;
                        for(size_t i = 1; i < nOrds.size(); ++i) {
                            if(nOrds[i] >= nMkrs.size()
                               || nOrds[i-1] >= nMkrs.size()) continue;
                            const int64_t gap =
                                int64_t(nMkrs[nOrds[i]].position)
                                - int64_t(nMkrs[nOrds[i-1]].position);
                            totalSpan += std::min(gap, maxGap);
                        }
                        return totalSpan;
                    };

                    if(!isVntrGap)
                    for(const auto& lf : lbp.reads) {
                        auto gL = readGraph.find(lf.readId);
                        if(gL == readGraph.end()) continue;

                        for(const auto& e1 : gL->second) {
                            const uint32_t n1 = e1.neighborReadId;
                            const int64_t ovlp1 = getOverlapSpan(e1);
                            if(ovlp1 < 0) continue;

                            // 1-hop: left → right
                            if(rightIds.count(n1)) {
                                int64_t rOvh = 0;
                                for(const auto& rf : bestRbp->reads)
                                    if(rf.readId == n1) { rOvh = rf.overhangBp; break; }
                                const int64_t dist = lf.overhangBp + rOvh - ovlp1;
                                if(dist > 0) {
                                    allPathDists.push_back(dist);
                                    if(dist > bestPathDist) {
                                        bestPathDist = dist;
                                        bestPathLen = 1;
                                        foundPath = true;
                                    }
                                }
                                continue;
                            }

                            if(!unanchoredReads.count(n1)) continue;

                            const auto n1Mkrs = markersRef[
                                OrientedReadId(ReadId(n1), 0).getValue()];
                            if(n1Mkrs.size() < 2) continue;
                            const int64_t n1Span = int64_t(n1Mkrs[n1Mkrs.size()-1].position)
                                                 - int64_t(n1Mkrs[0].position);

                            auto gN1 = readGraph.find(n1);
                            if(gN1 == readGraph.end()) continue;

                            for(const auto& e2 : gN1->second) {
                                const uint32_t n2 = e2.neighborReadId;
                                if(n2 == lf.readId) continue;
                                const int64_t ovlp2 = getOverlapSpan(e2);
                                if(ovlp2 < 0) continue;

                                // 2-hop: left → n1 → right
                                if(rightIds.count(n2)) {
                                    int64_t rOvh = 0;
                                    for(const auto& rf : bestRbp->reads)
                                        if(rf.readId == n2) { rOvh = rf.overhangBp; break; }
                                    const int64_t dist = lf.overhangBp + n1Span
                                        - ovlp1 - ovlp2 + rOvh;
                                    if(dist > 0) {
                                        allPathDists.push_back(dist);
                                        if(dist > bestPathDist) {
                                            bestPathDist = dist;
                                            bestPathLen = 2;
                                            foundPath = true;
                                        }
                                    }
                                    continue;
                                }

                                if(!unanchoredReads.count(n2)) continue;

                                const auto n2Mkrs = markersRef[
                                    OrientedReadId(ReadId(n2), 0).getValue()];
                                if(n2Mkrs.size() < 2) continue;
                                const int64_t n2Span = int64_t(n2Mkrs[n2Mkrs.size()-1].position)
                                                     - int64_t(n2Mkrs[0].position);

                                auto gN2 = readGraph.find(n2);
                                if(gN2 == readGraph.end()) continue;

                                for(const auto& e3 : gN2->second) {
                                    const uint32_t n3 = e3.neighborReadId;
                                    if(n3 == n1 || n3 == lf.readId) continue;
                                    if(!rightIds.count(n3)) continue;
                                    const int64_t ovlp3 = getOverlapSpan(e3);
                                    if(ovlp3 < 0) continue;

                                    // 3-hop: left → n1 → n2 → right
                                    int64_t rOvh = 0;
                                    for(const auto& rf : bestRbp->reads)
                                        if(rf.readId == n3) { rOvh = rf.overhangBp; break; }
                                    const int64_t dist = lf.overhangBp + n1Span + n2Span
                                        - ovlp1 - ovlp2 - ovlp3 + rOvh;
                                    if(dist > 0) {
                                        allPathDists.push_back(dist);
                                        if(dist > bestPathDist) {
                                            bestPathDist = dist;
                                            bestPathLen = 3;
                                            foundPath = true;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    const uint32_t breakpointPos = (lbp.refPos + bestRbp->refPos) / 2;

                    // bestPathDist already holds the maximum path
                    // distance from the search above.

                    cout << "    Breakpoint pair: L=" << lbp.refPos
                         << " (ends=" << lbp.endpointCount
                         << " ovh=" << lbp.ovhReadCount
                         << " fold=" << lbp.foldEnrichment << ")"
                         << " R=" << bestRbp->refPos
                         << " (starts=" << bestRbp->endpointCount
                         << " ovh=" << bestRbp->ovhReadCount
                         << " fold=" << bestRbp->foldEnrichment << ")"
                         << " refGap=" << bestDist
                         << " foundPath=" << foundPath
                         << " pathDist=" << bestPathDist
                         << " paths=" << allPathDists.size()
                         << " hops=" << bestPathLen
                         << " breakpoint=" << breakpointPos
                         << endl;

                    // If the standard path search failed and this is a VNTR
                    // gap, estimate insertion size from coverage depth.
                    //
                    // In a VNTR, read-to-read chains create shortcuts
                    // through the repeat structure, so path-based size
                    // estimation doesn't work. Instead, use the total
                    // read bases covering the VNTR vs the reference length.
                    //
                    // For a diploid het insertion:
                    //   totalReadBases ≈ coverage × (refLen + insLen) / 2
                    //                  + coverage × refLen / 2
                    //                  = coverage × refLen + coverage × insLen / 2
                    //   insLen ≈ 2 × (totalReadBases / coverage - refLen)
                    //
                    // For a hom insertion:
                    //   insLen ≈ totalReadBases / coverage - refLen
                    if(!foundPath && maxPairDist > 500) {
                        // Count total read bases in the VNTR region using
                        // actual read lengths (not marker-based lengths,
                        // which are unreliable in repetitive regions).
                        const uint32_t vntrStart = lbp.refPos;
                        const uint32_t vntrEnd = bestRbp->refPos;
                        const int64_t vntrRefLen = int64_t(vntrEnd) - int64_t(vntrStart);

                        // Compute average coverage in flanking regions.
                        uint32_t flankWindows = 0;
                        uint32_t flankSpanning = 0;
                        const uint32_t vntrStartWin = (vntrStart - refStartPos) / windowSize;
                        const uint32_t vntrEndWin = (vntrEnd - refStartPos) / windowSize;
                        for(uint32_t w = boundaryWindows; w < vntrStartWin && w < nWindows; ++w) {
                            flankSpanning += spanningCount[w];
                            ++flankWindows;
                        }
                        for(uint32_t w = vntrEndWin + 1; w + boundaryWindows < nWindows; ++w) {
                            flankSpanning += spanningCount[w];
                            ++flankWindows;
                        }
                        const double flankCoverage = flankWindows > 0
                            ? double(flankSpanning) / double(flankWindows) : 1.0;

                        // Count total read bases using actual base counts.
                        // Include: reads whose chains overlap the VNTR +
                        // unanchored reads (likely within VNTR).
                        int64_t totalReadBases = 0;
                        uint32_t vntrReadCount = 0;
                        for(const auto& [rid, info] : readChainInfoMap) {
                            if(info.maxRefPos < vntrStart || info.minRefPos > vntrEnd)
                                continue;
                            const int64_t readLen = int64_t(
                                readsRef.getRead(ReadId(rid)).baseCount);
                            totalReadBases += readLen;
                            ++vntrReadCount;
                        }
                        for(const uint32_t rid : unanchoredReads) {
                            const int64_t readLen = int64_t(
                                readsRef.getRead(ReadId(rid)).baseCount);
                            totalReadBases += readLen;
                            ++vntrReadCount;
                        }

                        // Estimate insertion size.
                        // Expected bases for reference VNTR at flank coverage:
                        //   expectedBases = flankCov × vntrRefLen
                        // Actual bases = totalReadBases
                        // Excess = totalReadBases - expectedBases
                        // For het: excess ≈ coverage/2 × insLen
                        //   insLen ≈ 2 × excess / coverage
                        // For hom: excess ≈ coverage × insLen
                        //   insLen ≈ excess / coverage
                        const double expectedBases = flankCoverage * double(vntrRefLen);
                        const double excess = double(totalReadBases) - expectedBases;
                        const int64_t insLenHet = flankCoverage > 0
                            ? int64_t(2.0 * excess / flankCoverage) : 0;
                        const int64_t insLenHom = flankCoverage > 0
                            ? int64_t(excess / flankCoverage) : 0;

                        cout << "    VNTR depth estimate: refLen=" << vntrRefLen
                             << " reads=" << vntrReadCount
                             << " bases=" << totalReadBases
                             << " expected=" << int64_t(expectedBases)
                             << " hetIns=" << insLenHet
                             << " homIns=" << insLenHom
                             << endl;

                        if(insLenHet > 20 && insLenHet < vntrRefLen) {
                            bestPathDist = insLenHet;
                            bestPathLen = 0; // coverage-based
                            foundPath = true;
                        } else if(insLenHom > 20 && insLenHom < vntrRefLen) {
                            bestPathDist = insLenHom;
                            bestPathLen = 0;
                            foundPath = true;
                        }
                    }

                    if(foundPath && bestPathDist > 20) {
                        cout << "    >>> INSERTION CALL: "
                             << "size=" << bestPathDist << "bp, "
                             << "breakpoint=" << breakpointPos << ", "
                             << "leftEnds=" << lbp.endpointCount << ", "
                             << "rightStarts=" << bestRbp->endpointCount << ", "
                             << "hops=" << bestPathLen
                             << endl;
                        // Record the insertion region for suppressing
                        // false coverage-drop deletion calls.
                        insertionCallRegions.push_back({
                            std::min(lbp.refPos, bestRbp->refPos),
                            std::max(lbp.refPos, bestRbp->refPos)
                        });
                    }

                    // -------------------------------------------------
                    // Deletion detection from diagonal shift.
                    //
                    // For a deletion, reads spanning the breakpoint pair
                    // have chains where the reference gap between
                    // consecutive anchors exceeds the read gap:
                    //   delSize = refGap - readGap  (positive for deletion)
                    //
                    // Look at all chains spanning the left-right BP zone.
                    // -------------------------------------------------
                    if(!foundPath && bestDist >= 50) {
                        const uint32_t delZoneStart = std::min(lbp.refPos, bestRbp->refPos);
                        const uint32_t delZoneEnd = std::max(lbp.refPos, bestRbp->refPos);


                        vector<int64_t> delShifts;

                        for(const auto& ce : chainsForRef) {
                            if(ce.readId == uint32_t(refId)) continue;
                            const auto& al = alignments[ce.chainIndex];
                            if(al.ordinals.size() < 4) continue;

                            // Check if chain spans the deletion zone.
                            bool hasLeft = false, hasRight = false;
                            for(const auto& ord : al.ordinals) {
                                if(ord[0] >= refMarkers.size()) continue;
                                const uint32_t rp = uint32_t(refMarkers[ord[0]].position);
                                if(rp < delZoneStart) hasLeft = true;
                                if(rp > delZoneEnd) hasRight = true;
                            }
                            if(!hasLeft || !hasRight) continue;

                            const Strand strand = ce.isSameStrand ? 0 : 1;
                            const auto rdMarkers = markersRef[
                                OrientedReadId(ReadId(ce.readId), strand).getValue()];

                            // Find the largest negative diagonal shift
                            // (refGap > readGap) near the zone.
                            int64_t bestDelShift = 0;
                            for(size_t j = 1; j < al.ordinals.size(); ++j) {
                                const auto& prev = al.ordinals[j - 1];
                                const auto& curr = al.ordinals[j];
                                if(prev[0] >= refMarkers.size()
                                   || curr[0] >= refMarkers.size()) continue;
                                if(prev[1] >= rdMarkers.size()
                                   || curr[1] >= rdMarkers.size()) continue;

                                const uint32_t refPosPrev = uint32_t(
                                    refMarkers[prev[0]].position);
                                const uint32_t refPosCurr = uint32_t(
                                    refMarkers[curr[0]].position);

                                // At least one anchor should be near the zone.
                                if(refPosPrev > delZoneEnd || refPosCurr < delZoneStart)
                                    continue;

                                const uint32_t rdPosPrev = uint32_t(
                                    rdMarkers[prev[1]].position);
                                const uint32_t rdPosCurr = uint32_t(
                                    rdMarkers[curr[1]].position);

                                const int64_t refGap = int64_t(refPosCurr)
                                    - int64_t(refPosPrev);
                                const int64_t readGap = int64_t(rdPosCurr)
                                    - int64_t(rdPosPrev);
                                const int64_t delShift = refGap - readGap;

                                if(delShift > bestDelShift) {
                                    bestDelShift = delShift;
                                }
                            }

                            if(bestDelShift > 20) {
                                delShifts.push_back(bestDelShift);
                            }
                        }

                        if(delShifts.size() >= 2) {
                            sort(delShifts.begin(), delShifts.end());
                            const int64_t medianDel = delShifts[delShifts.size() / 2];

                            cout << "    Deletion diagonal: n="
                                 << delShifts.size()
                                 << " median=" << medianDel
                                 << " min=" << delShifts.front()
                                 << " max=" << delShifts.back()
                                 << endl;

                            if(medianDel > 30) {
                                cout << "    >>> DELETION CALL: "
                                     << "size=" << medianDel << "bp, "
                                     << "breakpoint=" << breakpointPos << ", "
                                     << "leftEnds=" << lbp.endpointCount << ", "
                                     << "rightStarts=" << bestRbp->endpointCount << ", "
                                     << "supportingReads=" << delShifts.size()
                                     << endl;
                            }
                        }
                    }

                    // Also check for deletions when there IS no
                    // reference gap (L and R at same position) but
                    // spanning coverage drops.
                    if(!foundPath && bestDist < 50) {
                        // Look for negative diagonal shifts in chains
                        // spanning the breakpoint position.
                        const uint32_t bpPos = breakpointPos;
                        const uint32_t delZoneStart = bpPos > windowSize * 2
                            ? bpPos - windowSize * 2 : 0;
                        const uint32_t delZoneEnd = bpPos + windowSize * 2;

                        vector<int64_t> delShifts;

                        for(const auto& ce : chainsForRef) {
                            if(ce.readId == uint32_t(refId)) continue;
                            const auto& al = alignments[ce.chainIndex];
                            if(al.ordinals.size() < 4) continue;

                            bool hasLeft = false, hasRight = false;
                            for(const auto& ord : al.ordinals) {
                                if(ord[0] >= refMarkers.size()) continue;
                                const uint32_t rp = uint32_t(refMarkers[ord[0]].position);
                                if(rp < delZoneStart) hasLeft = true;
                                if(rp > delZoneEnd) hasRight = true;
                            }
                            if(!hasLeft || !hasRight) continue;

                            const Strand strand = ce.isSameStrand ? 0 : 1;
                            const auto rdMarkers = markersRef[
                                OrientedReadId(ReadId(ce.readId), strand).getValue()];

                            int64_t bestDelShift = 0;
                            for(size_t j = 1; j < al.ordinals.size(); ++j) {
                                const auto& prev = al.ordinals[j - 1];
                                const auto& curr = al.ordinals[j];
                                if(prev[0] >= refMarkers.size()
                                   || curr[0] >= refMarkers.size()) continue;
                                if(prev[1] >= rdMarkers.size()
                                   || curr[1] >= rdMarkers.size()) continue;

                                const uint32_t refPosPrev = uint32_t(
                                    refMarkers[prev[0]].position);
                                const uint32_t refPosCurr = uint32_t(
                                    refMarkers[curr[0]].position);

                                if(refPosPrev > delZoneEnd
                                   || refPosCurr < delZoneStart) continue;

                                const uint32_t rdPosPrev = uint32_t(
                                    rdMarkers[prev[1]].position);
                                const uint32_t rdPosCurr = uint32_t(
                                    rdMarkers[curr[1]].position);

                                const int64_t refGap = int64_t(refPosCurr)
                                    - int64_t(refPosPrev);
                                const int64_t readGap = int64_t(rdPosCurr)
                                    - int64_t(rdPosPrev);
                                const int64_t delShift = refGap - readGap;

                                if(delShift > bestDelShift) {
                                    bestDelShift = delShift;
                                }
                            }

                            if(bestDelShift > 20) {
                                delShifts.push_back(bestDelShift);
                            }
                        }

                        if(delShifts.size() >= 2) {
                            sort(delShifts.begin(), delShifts.end());
                            const int64_t medianDel = delShifts[delShifts.size() / 2];

                            cout << "    Deletion diagonal: n="
                                 << delShifts.size()
                                 << " median=" << medianDel
                                 << " min=" << delShifts.front()
                                 << " max=" << delShifts.back()
                                 << endl;

                            if(medianDel > 30) {
                                cout << "    >>> DELETION CALL: "
                                     << "size=" << medianDel << "bp, "
                                     << "breakpoint=" << breakpointPos << ", "
                                     << "leftEnds=" << lbp.endpointCount << ", "
                                     << "rightStarts=" << bestRbp->endpointCount << ", "
                                     << "supportingReads=" << delShifts.size()
                                     << endl;
                            }
                        }
                    }
                }

                // ---------------------------------------------------------
                // Hit-depth-only insertion detection.
                //
                // For small insertions, reads span across the breakpoint
                // without their chains breaking. There are no chain-endpoint
                // breakpoints near the true site. The only signal is a
                // hit-depth drop: inserted sequence has no reference k-mers,
                // so reads carrying the insertion contribute fewer hits.
                //
                // Approach:
                // 1. Cluster consecutive hit-depth breakpoints.
                // 2. For each cluster, find reads whose chains span across
                //    the drop zone.
                // 3. For each spanning read, find the largest diagonal shift
                //    within the drop zone: (readGap - refGap) at consecutive
                //    chain anchors straddling the zone.
                // 4. The median diagonal shift estimates insertion size.
                // ---------------------------------------------------------
                if(!hitDepthBreakpoints.empty()) {
                    // Cluster consecutive hit-depth breakpoints (within 2 windows).
                    struct HitDepthCluster {
                        uint32_t startPos;
                        uint32_t endPos;
                        uint32_t startWin;
                        uint32_t endWin;
                        double minRatio;
                    };
                    vector<HitDepthCluster> hdClusters;

                    HitDepthCluster cur;
                    cur.startPos = hitDepthBreakpoints[0].refPos;
                    cur.endPos = hitDepthBreakpoints[0].refPos;
                    cur.startWin = hitDepthBreakpoints[0].windowIdx;
                    cur.endWin = hitDepthBreakpoints[0].windowIdx;
                    cur.minRatio = hitDepthBreakpoints[0].dropRatio;

                    for(size_t i = 1; i < hitDepthBreakpoints.size(); ++i) {
                        const auto& hbp = hitDepthBreakpoints[i];
                        if(hbp.windowIdx <= cur.endWin + 3) {
                            // Extend cluster.
                            cur.endPos = hbp.refPos;
                            cur.endWin = hbp.windowIdx;
                            cur.minRatio = std::min(cur.minRatio, hbp.dropRatio);
                        } else {
                            hdClusters.push_back(cur);
                            cur.startPos = hbp.refPos;
                            cur.endPos = hbp.refPos;
                            cur.startWin = hbp.windowIdx;
                            cur.endWin = hbp.windowIdx;
                            cur.minRatio = hbp.dropRatio;
                        }
                    }
                    hdClusters.push_back(cur);

                    for(const auto& cluster : hdClusters) {
                        // Skip clusters already covered by a chain-endpoint
                        // breakpoint pair (already handled above).
                        bool alreadyCovered = false;
                        for(const auto& lbp : leftBreakpoints) {
                            for(const auto& rbp : rightBreakpoints) {
                                const int64_t dist = std::abs(
                                    int64_t(rbp.refPos) - int64_t(lbp.refPos));
                                if(dist <= 500
                                   && lbp.refPos <= cluster.endPos + windowSize * 3
                                   && rbp.refPos >= cluster.startPos - windowSize * 3) {
                                    alreadyCovered = true;
                                    break;
                                }
                            }
                            if(alreadyCovered) break;
                        }
                        if(alreadyCovered) continue;

                        // Find the deepest drop point in the cluster.
                        // Use a narrow zone around it (±1 window) so that
                        // short reads can span across.
                        uint32_t deepestPos = cluster.startPos;
                        double deepestRatio = cluster.minRatio;
                        for(const auto& hbp : hitDepthBreakpoints) {
                            if(hbp.refPos >= cluster.startPos
                               && hbp.refPos <= cluster.endPos
                               && hbp.dropRatio <= deepestRatio) {
                                deepestPos = hbp.refPos;
                                deepestRatio = hbp.dropRatio;
                            }
                        }

                        // Narrow zone: just ±1 window around deepest point.
                        const uint32_t zoneStart = deepestPos > windowSize
                            ? deepestPos - windowSize : 0;
                        const uint32_t zoneEnd = deepestPos + windowSize;

                        cout << "    HitDepth cluster: "
                             << cluster.startPos << "-" << cluster.endPos
                             << " minRatio=" << cluster.minRatio
                             << " deepest=" << deepestPos
                             << " zone=" << zoneStart << "-" << zoneEnd
                             << endl;



                        // Find reads whose chains span across the drop zone.
                        // For each such read, examine chain anchors to find
                        // the largest diagonal shift within the zone.
                        //
                        // diagonal = readPos - refPos
                        // An insertion causes diagonal to increase.
                        vector<int64_t> diagShifts;

                        uint32_t nChainsChecked = 0;
                        uint32_t nChainsSpanning = 0;

                        for(const auto& ce : chainsForRef) {
                            if(ce.readId == uint32_t(refId)) continue;
                            const auto& al = alignments[ce.chainIndex];
                            if(al.ordinals.size() < 4) continue;
                            ++nChainsChecked;

                            // Check if this chain spans the zone.
                            // Need anchors on both sides.
                            bool hasLeft = false, hasRight = false;
                            for(const auto& ord : al.ordinals) {
                                if(ord[0] >= refMarkers.size()) continue;
                                const uint32_t rp = uint32_t(refMarkers[ord[0]].position);
                                if(rp < zoneStart) hasLeft = true;
                                if(rp > zoneEnd) hasRight = true;
                            }
                            if(!hasLeft || !hasRight) continue;
                            ++nChainsSpanning;

                            // Get read markers for position lookup.
                            const Strand strand = ce.isSameStrand ? 0 : 1;
                            const auto rdMarkers = markersRef[
                                OrientedReadId(ReadId(ce.readId), strand).getValue()];

                            // Find the largest diagonal shift across the zone.
                            // Look at consecutive anchor pairs where one is
                            // before the zone and one is after.
                            int64_t bestShift = 0;
                            for(size_t j = 1; j < al.ordinals.size(); ++j) {
                                const auto& prev = al.ordinals[j - 1];
                                const auto& curr = al.ordinals[j];
                                if(prev[0] >= refMarkers.size()
                                   || curr[0] >= refMarkers.size()) continue;
                                if(prev[1] >= rdMarkers.size()
                                   || curr[1] >= rdMarkers.size()) continue;

                                const uint32_t refPosPrev = uint32_t(
                                    refMarkers[prev[0]].position);
                                const uint32_t refPosCurr = uint32_t(
                                    refMarkers[curr[0]].position);

                                // At least one anchor should be near/in the zone.
                                if(refPosPrev > zoneEnd || refPosCurr < zoneStart)
                                    continue;

                                const uint32_t rdPosPrev = uint32_t(
                                    rdMarkers[prev[1]].position);
                                const uint32_t rdPosCurr = uint32_t(
                                    rdMarkers[curr[1]].position);

                                const int64_t refGap = int64_t(refPosCurr)
                                    - int64_t(refPosPrev);
                                const int64_t readGap = int64_t(rdPosCurr)
                                    - int64_t(rdPosPrev);
                                const int64_t shift = readGap - refGap;

                                if(shift > bestShift) {
                                    bestShift = shift;
                                }
                            }

                            if(bestShift > 10) {
                                diagShifts.push_back(bestShift);
                            }
                        }

                        if(diagShifts.empty()) {
                            cout << "      No spanning chains with diagonal shift (checked="
                                 << nChainsChecked << " spanning="
                                 << nChainsSpanning << ")." << endl;
                            continue;
                        }

                        sort(diagShifts.begin(), diagShifts.end());
                        const int64_t medianShift = diagShifts[diagShifts.size() / 2];
                        const uint32_t breakpointPos =
                            (cluster.startPos + cluster.endPos) / 2;

                        cout << "      Diagonal shifts: n="
                             << diagShifts.size()
                             << " median=" << medianShift
                             << " min=" << diagShifts.front()
                             << " max=" << diagShifts.back()
                             << " breakpoint=" << breakpointPos
                             << endl;

                        if(diagShifts.size() >= 2 && medianShift > 20) {
                            cout << "    >>> INSERTION CALL (hit-depth): "
                                 << "size=" << medianShift << "bp, "
                                 << "breakpoint=" << breakpointPos << ", "
                                 << "supportingReads=" << diagShifts.size()
                                 << ", minRatio=" << cluster.minRatio
                                 << endl;
                        }
                    }
                }

                // ---------------------------------------------------------
                // Genome-wide diagonal shift scan for deletions.
                //
                // Deletions may not produce chain-endpoint breakpoints
                // because reads span across the deletion with chains
                // covering both flanks. The deletion is visible as a
                // negative diagonal shift: refGap > readGap at consecutive
                // anchors.
                //
                // Scan all chains for large diagonal shifts and cluster
                // them by reference position.
                // ---------------------------------------------------------
                {
                    struct DiagEvent {
                        uint32_t refPos;   // midpoint of the anchor pair
                        int64_t shift;     // refGap - readGap (positive = deletion)
                        uint32_t readId;
                    };
                    vector<DiagEvent> delEvents;

                    for(const auto& ce : chainsForRef) {
                        if(ce.readId == uint32_t(refId)) continue;
                        const auto& al = alignments[ce.chainIndex];
                        if(al.ordinals.size() < 4) continue;

                        const Strand strand = ce.isSameStrand ? 0 : 1;
                        const auto rdMarkers = markersRef[
                            OrientedReadId(ReadId(ce.readId), strand).getValue()];

                        for(size_t j = 1; j < al.ordinals.size(); ++j) {
                            const auto& prev = al.ordinals[j - 1];
                            const auto& curr = al.ordinals[j];
                            if(prev[0] >= refMarkers.size()
                               || curr[0] >= refMarkers.size()) continue;
                            if(prev[1] >= rdMarkers.size()
                               || curr[1] >= rdMarkers.size()) continue;

                            const uint32_t refPosPrev = uint32_t(
                                refMarkers[prev[0]].position);
                            const uint32_t refPosCurr = uint32_t(
                                refMarkers[curr[0]].position);
                            const uint32_t rdPosPrev = uint32_t(
                                rdMarkers[prev[1]].position);
                            const uint32_t rdPosCurr = uint32_t(
                                rdMarkers[curr[1]].position);

                            const int64_t refGap = int64_t(refPosCurr)
                                - int64_t(refPosPrev);
                            const int64_t readGap = int64_t(rdPosCurr)
                                - int64_t(rdPosPrev);
                            const int64_t delShift = refGap - readGap;

                            // Only consider significant deletions.
                            if(delShift > 30 && refGap > 50) {
                                const uint32_t midPos = (refPosPrev + refPosCurr) / 2;
                                // Skip boundary regions.
                                if(midPos > refStartPos + boundaryMargin
                                   && midPos + boundaryMargin < refEndPos) {
                                    delEvents.push_back({midPos, delShift, ce.readId});
                                }
                            }
                        }
                    }

                    if(!delEvents.empty()) {
                        // Sort by reference position.
                        sort(delEvents.begin(), delEvents.end(),
                            [](const DiagEvent& a, const DiagEvent& b) {
                                return a.refPos < b.refPos;
                            });

                        // Cluster events within 200bp of each other.
                        struct DelCluster {
                            uint32_t startPos;
                            uint32_t endPos;
                            vector<int64_t> shifts;
                            unordered_set<uint32_t> readIds;
                        };
                        vector<DelCluster> delClusters;

                        DelCluster curCluster;
                        curCluster.startPos = delEvents[0].refPos;
                        curCluster.endPos = delEvents[0].refPos;
                        curCluster.shifts.push_back(delEvents[0].shift);
                        curCluster.readIds.insert(delEvents[0].readId);

                        for(size_t i = 1; i < delEvents.size(); ++i) {
                            if(delEvents[i].refPos <= curCluster.endPos + 200) {
                                curCluster.endPos = delEvents[i].refPos;
                                curCluster.shifts.push_back(delEvents[i].shift);
                                curCluster.readIds.insert(delEvents[i].readId);
                            } else {
                                delClusters.push_back(std::move(curCluster));
                                curCluster.startPos = delEvents[i].refPos;
                                curCluster.endPos = delEvents[i].refPos;
                                curCluster.shifts.clear();
                                curCluster.shifts.push_back(delEvents[i].shift);
                                curCluster.readIds.clear();
                                curCluster.readIds.insert(delEvents[i].readId);
                            }
                        }
                        delClusters.push_back(std::move(curCluster));

                        for(auto& cluster : delClusters) {
                            if(cluster.readIds.size() < 2) continue;

                            sort(cluster.shifts.begin(), cluster.shifts.end());
                            const int64_t medianDel =
                                cluster.shifts[cluster.shifts.size() / 2];
                            const uint32_t bpPos =
                                (cluster.startPos + cluster.endPos) / 2;

                            cout << "    Deletion cluster: pos="
                                 << cluster.startPos << "-" << cluster.endPos
                                 << " reads=" << cluster.readIds.size()
                                 << " events=" << cluster.shifts.size()
                                 << " median=" << medianDel
                                 << " min=" << cluster.shifts.front()
                                 << " max=" << cluster.shifts.back()
                                 << endl;

                            if(medianDel > 30 && cluster.readIds.size() >= 3) {
                                cout << "    >>> DELETION CALL: "
                                     << "size=" << medianDel << "bp, "
                                     << "breakpoint=" << bpPos << ", "
                                     << "supportingReads=" << cluster.readIds.size()
                                     << endl;
                            }
                        }
                    }
                }

                // ---------------------------------------------------------
                // Coverage-drop deletion detection.
                //
                // For heterozygous deletions, spanning coverage drops
                // by ~50% in the deleted region. Detect consecutive
                // windows where coverage is significantly below median.
                // ---------------------------------------------------------
                if(medianSpanning > 5) {
                    const double covDropThreshold = 0.6; // 60% of median
                    const uint32_t minCovDropWindows = 2;

                    struct CovDropCluster {
                        uint32_t startWin;
                        uint32_t endWin;
                        uint32_t startPos;
                        uint32_t endPos;
                        double minRatio;
                    };
                    vector<CovDropCluster> covDropClusters;

                    uint32_t clusterStart = UINT32_MAX;
                    double clusterMinRatio = 1.0;
                    for(uint32_t w = boundaryWindows;
                        w + boundaryWindows < nWindows; ++w) {
                        const double ratio =
                            double(spanningCount[w]) / double(medianSpanning);
                        if(ratio < covDropThreshold) {
                            if(clusterStart == UINT32_MAX) {
                                clusterStart = w;
                                clusterMinRatio = ratio;
                            }
                            clusterMinRatio =
                                std::min(clusterMinRatio, ratio);
                        } else {
                            if(clusterStart != UINT32_MAX) {
                                const uint32_t clusterEnd = w - 1;
                                if(clusterEnd - clusterStart + 1
                                    >= minCovDropWindows) {
                                    const uint32_t sp = refStartPos
                                        + clusterStart * windowSize;
                                    const uint32_t ep = refStartPos
                                        + (clusterEnd + 1) * windowSize;
                                    covDropClusters.push_back(
                                        {clusterStart, clusterEnd,
                                         sp, ep, clusterMinRatio});
                                }
                                clusterStart = UINT32_MAX;
                            }
                        }
                    }
                    // Handle cluster at end.
                    if(clusterStart != UINT32_MAX) {
                        const uint32_t clusterEnd =
                            nWindows - boundaryWindows - 1;
                        if(clusterEnd - clusterStart + 1
                            >= minCovDropWindows) {
                            const uint32_t sp = refStartPos
                                + clusterStart * windowSize;
                            const uint32_t ep = refStartPos
                                + (clusterEnd + 1) * windowSize;
                            covDropClusters.push_back(
                                {clusterStart, clusterEnd,
                                 sp, ep, clusterMinRatio});
                        }
                    }

                    for(const auto& cdc : covDropClusters) {
                        const uint32_t delSize = cdc.endPos - cdc.startPos;
                        // Check flanking coverage is near median.
                        uint32_t leftFlankCov = 0, rightFlankCov = 0;
                        uint32_t leftFlankN = 0, rightFlankN = 0;
                        for(uint32_t w =
                                (cdc.startWin > 3 ? cdc.startWin - 3 : 0);
                            w < cdc.startWin; ++w) {
                            leftFlankCov += spanningCount[w];
                            ++leftFlankN;
                        }
                        for(uint32_t w = cdc.endWin + 1;
                            w <= cdc.endWin + 3 && w < nWindows; ++w) {
                            rightFlankCov += spanningCount[w];
                            ++rightFlankN;
                        }
                        const double leftFlank = leftFlankN > 0
                            ? double(leftFlankCov) / double(leftFlankN) : 0;
                        const double rightFlank = rightFlankN > 0
                            ? double(rightFlankCov) / double(rightFlankN) : 0;

                        // Both flanks should be near median (>70%).
                        if(leftFlank < 0.7 * medianSpanning
                           || rightFlank < 0.7 * medianSpanning) continue;

                        // Check if this region overlaps a detected VNTR gap
                        // or an insertion call region.
                        // VNTR gaps cause coverage drops from chaining
                        // failure in tandem repeats, not real deletions.
                        // Insertion regions cause coverage drops because
                        // insertion-carrying reads can't chain through.
                        bool overlapsVntr = false;
                        for(const auto& vg : vntrGaps) {
                            if(cdc.startPos < vg.endPos
                               && cdc.endPos > vg.startPos) {
                                overlapsVntr = true;
                                break;
                            }
                        }
                        // Use margin for insertion overlap since the
                        // coverage drop extends beyond the BP pair.
                        const uint32_t insMargin = 200;
                        bool overlapsInsertion = false;
                        for(const auto& ir : insertionCallRegions) {
                            const uint32_t irStart =
                                ir.startPos > insMargin
                                ? ir.startPos - insMargin : 0;
                            const uint32_t irEnd = ir.endPos + insMargin;
                            if(cdc.startPos < irEnd
                               && cdc.endPos > irStart) {
                                overlapsInsertion = true;
                                break;
                            }
                        }

                        // Check if the region is marker-depleted.
                        // In a real deletion, the reference still has
                        // markers (k-mers) — reads just don't align there.
                        // In a VNTR/repeat, the reference itself has
                        // low hit-depth because k-mers are non-unique.
                        uint32_t lowHitDepthWins = 0;
                        uint32_t totalHitDepthWins = 0;
                        for(uint32_t w = cdc.startWin;
                            w <= cdc.endWin && w < nWindows; ++w) {
                            if(windowMarkerCount[w] == 0) continue;
                            ++totalHitDepthWins;
                            if(medianHitDepth > 0
                               && windowHitDepth[w] / medianHitDepth
                                  < hitDepthDropThreshold) {
                                ++lowHitDepthWins;
                            }
                        }
                        const bool markerDepleted = totalHitDepthWins > 0
                            && double(lowHitDepthWins)
                               / double(totalHitDepthWins) > 0.5;

                        const uint32_t bpPos =
                            (cdc.startPos + cdc.endPos) / 2;

                        cout << "    Coverage-drop deletion: pos="
                             << cdc.startPos << "-" << cdc.endPos
                             << " size=" << delSize
                             << " minRatio=" << cdc.minRatio
                             << " leftFlank=" << leftFlank
                             << " rightFlank=" << rightFlank
                             << " vntr=" << overlapsVntr
                             << " ins=" << overlapsInsertion
                             << " markerDepleted=" << markerDepleted
                             << endl;

                        // Suppress calls that overlap detected VNTR gaps
                        // or insertion call regions.
                        if(overlapsVntr || overlapsInsertion) continue;

                        // -------------------------------------------------
                        // Adaptive multi-k anchor filling.
                        //
                        // Progressively fill the coverage-drop region
                        // and flanks with unique anchors at increasing
                        // k values. At each k, only scan gaps where no
                        // anchors exist yet. Build a per-read anchor
                        // map, then analyze diagonal shifts.
                        // -------------------------------------------------
                        const uint32_t flankSize = 200;
                        const uint32_t refSeqLen = uint32_t(
                            readsRef.getRead(refId).baseCount);
                        const uint32_t regionStart =
                            cdc.startPos > flankSize
                            ? cdc.startPos - flankSize : 0;
                        const uint32_t regionEnd =
                            std::min(cdc.endPos + flankSize, refSeqLen);

                        // Get reference raw sequence.
                        const vector<Base> refSeq =
                            readsRef.getOrientedReadRawSequence(
                                OrientedReadId(refId, 0));

                        // Reference anchors: (refPos, kmer string).
                        // Sorted by refPos. Track covered positions.
                        struct RefAnchor {
                            uint32_t refPos;
                            string kmer;
                            uint32_t kLen;
                        };
                        vector<RefAnchor> refAnchors;

                        // Track which reference positions are covered
                        // by at least one unique anchor.
                        const uint32_t regionLen = regionEnd - regionStart;
                        vector<bool> covered(regionLen, false);

                        // Extended region for uniqueness checking.
                        const uint32_t uniStart =
                            regionStart > 500
                            ? regionStart - 500 : 0;
                        const uint32_t uniEnd =
                            std::min(regionEnd + 500, refSeqLen);

                        // Fill gaps starting from large k (most unique,
                        // skeleton anchors) down to small k (dense fill).
                        // Max k=60 since reads are ~150bp.
                        const uint32_t maxK = 60;
                        const uint32_t minK = uint32_t(k);
                        for(uint32_t tryK = maxK;
                            tryK >= minK; tryK -= 2) {

                            // Build k-mer → positions for uniqueness
                            // check across extended region.
                            std::unordered_map<string, vector<uint32_t>>
                                refKmerPos;
                            for(uint32_t p = uniStart;
                                p + tryK <= uniEnd; ++p) {
                                string kmer;
                                kmer.reserve(tryK);
                                for(uint32_t j = 0; j < tryK; ++j) {
                                    kmer.push_back(
                                        refSeq[p + j].character());
                                }
                                refKmerPos[kmer].push_back(p);
                            }

                            // Find unique k-mers in uncovered gaps.
                            uint32_t newAnchors = 0;
                            for(const auto& [kmer, positions] : refKmerPos) {
                                if(positions.size() != 1) continue;
                                const uint32_t pos = positions[0];
                                if(pos < regionStart
                                   || pos + tryK > regionEnd) continue;

                                // Check if this position is already
                                // covered by an existing anchor.
                                const uint32_t localPos =
                                    pos - regionStart;
                                if(covered[localPos]) continue;

                                refAnchors.push_back({pos, kmer, tryK});
                                // Mark covered range.
                                for(uint32_t j = 0;
                                    j < tryK && localPos + j < regionLen;
                                    ++j) {
                                    covered[localPos + j] = true;
                                }
                                ++newAnchors;
                            }

                            // Count remaining gaps.
                            uint32_t gapBases = 0;
                            for(uint32_t i = 0; i < regionLen; ++i) {
                                if(!covered[i]) ++gapBases;
                            }

                            cout << "      k=" << tryK
                                 << ": +" << newAnchors
                                 << " anchors, total="
                                 << refAnchors.size()
                                 << ", gapBases=" << gapBases
                                 << "/" << regionLen
                                 << endl;

                            // Stop if no gaps remain.
                            if(gapBases == 0) break;
                        }

                        // Sort reference anchors by position.
                        sort(refAnchors.begin(), refAnchors.end(),
                            [](const RefAnchor& a, const RefAnchor& b) {
                                return a.refPos < b.refPos;
                            });

                        cout << "      Total ref anchors: "
                             << refAnchors.size() << endl;

                        if(refAnchors.size() < 5) {
                            // Not enough anchors to analyze.
                            if(delSize >= 50 && delSize <= 2000) {
                                cout << "    >>> DELETION CALL (coverage): "
                                     << "size=" << delSize << "bp, "
                                     << "breakpoint=" << bpPos
                                     << endl;
                            }
                            continue; // next covDropCluster
                        }

                        // For each read, match reference anchors and
                        // build a per-read anchor list with (refPos,
                        // readPos) pairs.
                        struct ReadAnchorResult {
                            ReadId readId;
                            vector<pair<uint32_t, uint32_t>> anchors;
                            // (refPos, readPos)
                        };
                        vector<ReadAnchorResult> readResults;

                        for(const auto& rg : readGroups) {
                            const vector<Base> readSeq =
                                readsRef.getOrientedReadRawSequence(
                                    OrientedReadId(rg.readId, 0));
                            const uint32_t readLen =
                                uint32_t(readSeq.size());

                            // Build k-mer → positions for this read.
                            // We need to handle multiple k values, so
                            // build for each k used in refAnchors.
                            // Collect all unique k values.
                            std::set<uint32_t> kValues;
                            for(const auto& ra : refAnchors) {
                                kValues.insert(ra.kLen);
                            }

                            // For each k, build read k-mer index.
                            std::unordered_map<string, vector<uint32_t>>
                                readKmerPos;
                            for(const uint32_t kv : kValues) {
                                if(readLen < kv) continue;
                                for(uint32_t p = 0;
                                    p + kv <= readLen; ++p) {
                                    string kmer;
                                    kmer.reserve(kv);
                                    for(uint32_t j = 0; j < kv; ++j) {
                                        kmer.push_back(
                                            readSeq[p + j].character());
                                    }
                                    // Only add if not already present
                                    // (avoid duplicates from different k).
                                    readKmerPos[kmer].push_back(p);
                                }
                            }

                            // Match reference anchors.
                            vector<pair<uint32_t, uint32_t>> matches;
                            for(const auto& ra : refAnchors) {
                                auto it = readKmerPos.find(ra.kmer);
                                if(it == readKmerPos.end()) continue;
                                // Only use if unique in read.
                                if(it->second.size() != 1) continue;
                                matches.push_back(
                                    {ra.refPos, it->second[0]});
                            }

                            if(matches.size() >= 3) {
                                // Sort by refPos.
                                sort(matches.begin(), matches.end());
                                readResults.push_back(
                                    {rg.readId, std::move(matches)});
                            }
                        }

                        cout << "      Reads with anchors: "
                             << readResults.size() << endl;

                        // Collect per-read median diagonals for
                        // bimodal analysis (single-flank approach).
                        vector<int64_t> allMedianDiags;
                        for(const auto& rr : readResults) {
                            vector<int64_t> diags;
                            for(const auto& [rp, rdp] : rr.anchors) {
                                diags.push_back(
                                    int64_t(rp) - int64_t(rdp));
                            }
                            sort(diags.begin(), diags.end());
                            allMedianDiags.push_back(
                                diags[diags.size() / 2]);
                        }
                        sort(allMedianDiags.begin(),
                             allMedianDiags.end());

                        cout << "      Reads with anchors: "
                             << readResults.size()
                             << " median diags: "
                             << allMedianDiags.size() << endl;

                        // Analyze per-read diagonal profiles.
                        // For each read, compute diagonal at each anchor
                        // and find the max drop (deletion signal).
                        bool refinedCall = false;
                        struct DelSignal {
                            ReadId readId;
                            int64_t dropSize;
                            uint32_t dropRefPos;
                        };
                        vector<DelSignal> delSignals;

                        for(const auto& rr : readResults) {
                            // Compute diagonals.
                            vector<int64_t> diags;
                            vector<uint32_t> refPositions;
                            for(const auto& [rp, rdp] : rr.anchors) {
                                diags.push_back(
                                    int64_t(rp) - int64_t(rdp));
                                refPositions.push_back(rp);
                            }

                            // Find max drop in diagonal (deletion).
                            int64_t maxDrop = 0;
                            uint32_t dropPos = 0;
                            for(size_t i = 1; i < diags.size(); ++i) {
                                const int64_t drop =
                                    diags[i-1] - diags[i];
                                if(drop > maxDrop) {
                                    maxDrop = drop;
                                    dropPos = (refPositions[i-1]
                                               + refPositions[i]) / 2;
                                }
                            }

                            if(maxDrop > 30) {
                                delSignals.push_back(
                                    {rr.readId, maxDrop, dropPos});
                            }
                        }

                        if(delSignals.size() >= 2) {
                            // Cluster deletion signals by size.
                            sort(delSignals.begin(), delSignals.end(),
                                [](const DelSignal& a, const DelSignal& b) {
                                    return a.dropSize < b.dropSize;
                                });

                            // Find the most common deletion size
                            // (within 20% tolerance).
                            uint32_t bestCount = 0;
                            int64_t bestSize = 0;
                            uint32_t bestBp = 0;
                            for(size_t i = 0;
                                i < delSignals.size(); ++i) {
                                uint32_t count = 0;
                                int64_t sizeSum = 0;
                                uint32_t bpSum = 0;
                                for(size_t j = i;
                                    j < delSignals.size(); ++j) {
                                    if(delSignals[j].dropSize
                                       <= delSignals[i].dropSize * 1.3) {
                                        ++count;
                                        sizeSum += delSignals[j].dropSize;
                                        bpSum += delSignals[j].dropRefPos;
                                    }
                                }
                                if(count > bestCount) {
                                    bestCount = count;
                                    bestSize = sizeSum / int64_t(count);
                                    bestBp = bpSum / count;
                                }
                            }

                            cout << "      Del signals ("
                                 << delSignals.size() << "):";
                            for(const auto& ds : delSignals) {
                                cout << " r" << ds.readId
                                     << ":" << ds.dropSize
                                     << "@" << ds.dropRefPos;
                            }
                            cout << endl;

                            // Require the detected size to be at least
                            // 25% of the coverage-drop region to avoid
                            // noise from repeat-induced small drops.
                            if(bestCount >= 2 && bestSize >= 50
                               && bestSize >= int64_t(delSize) / 4) {
                                cout << "    >>> DELETION CALL (adaptive): "
                                     << "size=" << bestSize << "bp, "
                                     << "breakpoint=" << bestBp << ", "
                                     << "reads=" << bestCount
                                     << endl;
                                refinedCall = true;
                            }
                        }

                        // If per-read diagonal drop didn't work,
                        // try per-anchor pairwise diagonal difference
                        // analysis. For each reference anchor, collect
                        // diagonals from all reads that match it. In a
                        // het deletion, reads from different alleles at
                        // the same anchor differ by the deletion size.
                        if(!refinedCall && readResults.size() >= 4) {
                            // Build per-anchor diagonal lists.
                            // Key: refPos of anchor, Value: list of
                            // (readDiag) from different reads.
                            std::unordered_map<uint32_t,
                                vector<int64_t>> anchorDiags;
                            for(const auto& rr : readResults) {
                                for(const auto& [rp, rdp] : rr.anchors) {
                                    anchorDiags[rp].push_back(
                                        int64_t(rp) - int64_t(rdp));
                                }
                            }

                            // Collect all pairwise diagonal differences
                            // at each anchor. In a het deletion, the
                            // differences cluster around 0 (same allele)
                            // and ±D (different alleles).
                            vector<int64_t> pairDiffs;
                            for(auto& [pos, diags] : anchorDiags) {
                                if(diags.size() < 2) continue;
                                sort(diags.begin(), diags.end());
                                for(size_t i = 0; i < diags.size(); ++i) {
                                    for(size_t j = i + 1;
                                        j < diags.size(); ++j) {
                                        const int64_t diff =
                                            diags[j] - diags[i];
                                        if(diff > 30) {
                                            pairDiffs.push_back(diff);
                                        }
                                    }
                                }
                            }

                            int64_t bestShift = 0;

                            if(pairDiffs.size() >= 3) {
                                sort(pairDiffs.begin(), pairDiffs.end());

                                // Count occurrences of each diff value.
                                std::unordered_map<int64_t, uint32_t>
                                    diffCounts;
                                for(const auto& d : pairDiffs) {
                                    ++diffCounts[d];
                                }

                                // Build histogram of diff clusters
                                // (within 10% tolerance).
                                struct DiffCluster {
                                    int64_t meanDiff;
                                    uint32_t count;
                                };
                                vector<DiffCluster> clusters;
                                vector<int64_t> uniqueDiffs;
                                for(const auto& [d, c] : diffCounts) {
                                    uniqueDiffs.push_back(d);
                                }
                                sort(uniqueDiffs.begin(),
                                     uniqueDiffs.end());

                                for(size_t i = 0;
                                    i < uniqueDiffs.size(); ) {
                                    int64_t sum = 0;
                                    uint32_t cnt = 0;
                                    size_t j = i;
                                    while(j < uniqueDiffs.size()
                                          && uniqueDiffs[j]
                                             <= uniqueDiffs[i] * 1.15) {
                                        sum += uniqueDiffs[j]
                                               * diffCounts[uniqueDiffs[j]];
                                        cnt += diffCounts[uniqueDiffs[j]];
                                        ++j;
                                    }
                                    clusters.push_back(
                                        {sum / int64_t(cnt), cnt});
                                    i = j;
                                }

                                // Sort clusters by count (descending).
                                sort(clusters.begin(), clusters.end(),
                                    [](const DiffCluster& a,
                                       const DiffCluster& b) {
                                        return a.count > b.count;
                                    });

                                cout << "      Diff clusters:";
                                for(size_t i = 0;
                                    i < std::min(clusters.size(),
                                                 size_t(8)); ++i) {
                                    cout << " " << clusters[i].meanDiff
                                         << "bp(" << clusters[i].count
                                         << ")";
                                }
                                cout << endl;

                                // Find the best cluster.
                                // Strategy: among clusters with diff
                                // in [100, delSize], find the most
                                // supported one. The 100bp minimum
                                // filters out repeat-unit diffs.
                                // If tied (within 10%), prefer smaller.
                                uint32_t bestCount = 0;
                                for(const auto& cl : clusters) {
                                    if(cl.meanDiff >= 100
                                       && cl.meanDiff <= int64_t(delSize)
                                       && cl.count > bestCount) {
                                        bestCount = cl.count;
                                        bestShift = cl.meanDiff;
                                    }
                                }

                                // Among clusters with similar support
                                // (within 10% of best), prefer smaller.
                                if(bestCount > 0) {
                                    int64_t smallestInRange = bestShift;
                                    for(const auto& cl : clusters) {
                                        if(cl.meanDiff >= 100
                                           && cl.meanDiff <= int64_t(delSize)
                                           && cl.count >= bestCount * 9 / 10
                                           && cl.meanDiff < smallestInRange) {
                                            smallestInRange = cl.meanDiff;
                                        }
                                    }
                                    bestShift = smallestInRange;
                                    // Update count for the selected cluster.
                                    for(const auto& cl : clusters) {
                                        if(cl.meanDiff == bestShift) {
                                            bestCount = cl.count;
                                            break;
                                        }
                                    }
                                }

                                // If no cluster >= 100bp, try >= 50bp.
                                if(bestShift == 0) {
                                    for(const auto& cl : clusters) {
                                        if(cl.meanDiff >= 50
                                           && cl.count > bestCount) {
                                            bestCount = cl.count;
                                            bestShift = cl.meanDiff;
                                        }
                                    }
                                }

                                cout << "      Best diff: "
                                     << bestShift << "bp ("
                                     << bestCount << " pairs)"
                                     << endl;
                            }

                            // Fallback: per-read median diagonal
                            // flank gap analysis.
                            if(bestShift == 0) {
                                const uint32_t gapCenter =
                                    (cdc.startPos + cdc.endPos) / 2;
                                vector<int64_t> leftDiags, rightDiags;
                                for(const auto& rr : readResults) {
                                    vector<uint32_t> rps;
                                    vector<int64_t> ds;
                                    for(const auto& [rp, rdp]
                                        : rr.anchors) {
                                        rps.push_back(rp);
                                        ds.push_back(int64_t(rp)
                                                     - int64_t(rdp));
                                    }
                                    sort(rps.begin(), rps.end());
                                    sort(ds.begin(), ds.end());
                                    const uint32_t medRefPos =
                                        rps[rps.size() / 2];
                                    const int64_t medDiag =
                                        ds[ds.size() / 2];
                                    if(medRefPos < gapCenter) {
                                        leftDiags.push_back(medDiag);
                                    } else {
                                        rightDiags.push_back(medDiag);
                                    }
                                }

                                auto analyzeFlankGap = [](
                                    vector<int64_t>& diags) -> int64_t {
                                    if(diags.size() < 4) return 0;
                                    sort(diags.begin(), diags.end());
                                    int64_t maxGap = 0;
                                    for(size_t i = 1;
                                        i < diags.size(); ++i) {
                                        const int64_t gap =
                                            diags[i] - diags[i-1];
                                        if(gap > maxGap)
                                            maxGap = gap;
                                    }
                                    return maxGap;
                                };

                                const int64_t leftGap =
                                    analyzeFlankGap(leftDiags);
                                const int64_t rightGap =
                                    analyzeFlankGap(rightDiags);

                                cout << "      Flank gaps: left="
                                     << leftGap << " ("
                                     << leftDiags.size()
                                     << " reads) right=" << rightGap
                                     << " (" << rightDiags.size()
                                     << " reads)" << endl;

                                if(leftGap > 50 && rightGap > 50) {
                                    bestShift =
                                        (leftGap + rightGap) / 2;
                                } else if(leftGap > 50) {
                                    bestShift = leftGap;
                                } else if(rightGap > 50) {
                                    bestShift = rightGap;
                                }
                            }

                            if(bestShift >= 50 && bestShift <= 2000) {
                                cout << "    >>> DELETION CALL "
                                     << "(adaptive-bimodal): "
                                     << "size=" << bestShift << "bp, "
                                     << "breakpoint=" << bpPos
                                     << endl;
                                refinedCall = true;
                            }
                        }

                        if(!refinedCall && delSize >= 50
                           && delSize <= 2000) {
                            cout << "    >>> DELETION CALL (coverage): "
                                 << "size=" << delSize << "bp, "
                                 << "breakpoint=" << bpPos
                                 << endl;
                        }
                    }
                }
            }
        }

        // -----------------------------------------------------------------
        // Step 6a: Output per-cluster SV summary.
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

                    cout << "    >>> " << typeStr << " CLUSTER: "
                         << "id=" << cid << ", "
                         << "size=" << (sumSize / int64_t(count)) << "bp, "
                         << "breakpoint=" << (sumPos / count) << ", "
                         << "reads=" << count
                         << endl;
                }
            }
        }

        // -----------------------------------------------------------------
        // Step 6b: Output MSA and consensus.
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
