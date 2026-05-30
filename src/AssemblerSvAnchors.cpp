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
#include "Sdust.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"

#include <htslib/sam.h>
#include <htslib/hts.h>

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


// Simple de Bruijn graph assembly of soft-clip sequences.
// Returns the length of the longest assembled contig.
// Uses a greedy extension approach: start from the most
// frequent k-mer and extend in both directions.
static uint32_t assembleClipSequences(
    const vector<string>& seqs,
    uint32_t k = 21)
{
    if(seqs.empty() || k < 11) return 0;

    // Build k-mer count map.
    unordered_map<string, uint32_t> kmerCount;
    for(const auto& seq : seqs) {
        if(seq.size() < k) continue;
        for(size_t i = 0; i + k <= seq.size(); ++i) {
            ++kmerCount[seq.substr(i, k)];
        }
    }

    if(kmerCount.empty()) return 0;

    // Find the most frequent k-mer as seed.
    string seed;
    uint32_t maxCount = 0;
    for(const auto& kv : kmerCount) {
        if(kv.second > maxCount) {
            maxCount = kv.second;
            seed = kv.first;
        }
    }

    if(maxCount < 2) return 0;

    // Greedy extension: extend right, then left.
    string contig = seed;
    unordered_set<string> used;
    used.insert(seed);

    // Extend right.
    while(true) {
        const string suffix = contig.substr(
            contig.size() - (k - 1));
        string bestNext;
        uint32_t bestCount = 0;
        for(char c : {'A', 'C', 'G', 'T'}) {
            string candidate = suffix + c;
            auto it = kmerCount.find(candidate);
            if(it != kmerCount.end()
               && it->second > bestCount
               && used.find(candidate) == used.end()) {
                bestCount = it->second;
                bestNext = candidate;
            }
        }
        if(bestCount < 2) break;
        used.insert(bestNext);
        contig += bestNext.back();
        if(contig.size() > 5000) break;
    }

    // Extend left.
    while(true) {
        const string prefix = contig.substr(0, k - 1);
        string bestPrev;
        uint32_t bestCount = 0;
        for(char c : {'A', 'C', 'G', 'T'}) {
            string candidate = string(1, c) + prefix;
            auto it = kmerCount.find(candidate);
            if(it != kmerCount.end()
               && it->second > bestCount
               && used.find(candidate) == used.end()) {
                bestCount = it->second;
                bestPrev = candidate;
            }
        }
        if(bestCount < 2) break;
        used.insert(bestPrev);
        contig = bestPrev.front() + contig;
        if(contig.size() > 5000) break;
    }

    return uint32_t(contig.size());
}



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
    const vector<RefHitDepthWindow>& refHitDepth,
    const string& bamFileName)
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

        const OrientedReadId refOid(refId, 0);
        const auto refMarkers = markersRef[refOid.getValue()];
        const uint32_t refLength = uint32_t(readsRef.getRead(refId).baseCount);

        // Parse SA tag SV evidence from BAM if provided.
        vector<SaTagSvCall> saTagCalls;
        vector<SoftClipBreakpoint> softClipBPs;
        vector<CigarIndelCall> cigarIndels;

        // DEL call record used by delCallRecords and allDelCalls.
        struct DelCallRecord {
            uint32_t breakpointPos;
            int64_t size;
            uint32_t readCount;
            string source;
        };

        // All DEL calls, for emitting ±windowSize adj variants
        // at the end. Populated by all DEL emission sites.
        vector<DelCallRecord> allDelCalls;

        // All INS calls, for INS-to-DEL flipping at the end.
        vector<DelCallRecord> allInsCalls;



        // Lambda to detect deletions from split-read chain
        // diagonal differences. Can be called before early
        // continues when the MSA/segment path is not reached.
        auto detectSplitReadDels = [&]() {
            struct SplitDelEvent {
                uint32_t refPos;
                int64_t delSize;
                uint32_t readId;
            };
            vector<SplitDelEvent> splitEvents;

            std::unordered_map<uint32_t,
                vector<uint32_t>> readChainIdx;
            for(uint32_t ci = 0;
                ci < chainsForRef.size(); ++ci) {
                const auto& ce = chainsForRef[ci];
                if(ce.readId == uint32_t(refId)) continue;
                readChainIdx[ce.readId].push_back(ci);
            }

            for(const auto& [rdId, cis] : readChainIdx) {
                if(cis.size() < 2) continue;
                struct ChainSummary {
                    int64_t medDiag;
                    uint32_t minRefPos;
                    uint32_t maxRefPos;
                    uint32_t nAnchors;
                };
                vector<ChainSummary> summaries;
                for(const auto ci : cis) {
                    const auto& ce = chainsForRef[ci];
                    if(!ce.isSameStrand) continue;
                    const auto& al =
                        alignments[ce.chainIndex];
                    if(al.ordinals.size() < 3) continue;
                    const Strand strand = 0;
                    const auto rdMkrs = markersRef[
                        OrientedReadId(
                            ReadId(rdId), strand
                        ).getValue()];
                    vector<int64_t> diags;
                    uint32_t minRp = UINT32_MAX;
                    uint32_t maxRp = 0;
                    for(const auto& ord : al.ordinals) {
                        if(ord[0] >= refMarkers.size()
                           || ord[1] >= rdMkrs.size())
                            continue;
                        const uint32_t rp = uint32_t(
                            refMarkers[ord[0]].position);
                        const uint32_t qp = uint32_t(
                            rdMkrs[ord[1]].position);
                        diags.push_back(
                            int64_t(rp) - int64_t(qp));
                        minRp = std::min(minRp, rp);
                        maxRp = std::max(maxRp, rp);
                    }
                    if(diags.size() < 3) continue;
                    sort(diags.begin(), diags.end());
                    summaries.push_back({
                        diags[diags.size() / 2],
                        minRp, maxRp,
                        uint32_t(diags.size())});
                }
                if(summaries.size() < 2) continue;
                sort(summaries.begin(), summaries.end(),
                    [](const ChainSummary& a,
                       const ChainSummary& b) {
                        return a.minRefPos < b.minRefPos;
                    });
                for(size_t a = 0;
                    a + 1 < summaries.size(); ++a) {
                    const auto& ca = summaries[a];
                    const auto& cb = summaries[a + 1];
                    if(cb.minRefPos <= ca.maxRefPos)
                        continue;
                    const int64_t refGap =
                        int64_t(cb.minRefPos)
                        - int64_t(ca.maxRefPos);
                    const int64_t dd =
                        cb.medDiag - ca.medDiag;
                    if(dd > 20 && refGap < dd * 2) {
                        const uint32_t gapMid =
                            (ca.maxRefPos + cb.minRefPos) / 2;
                        splitEvents.push_back(
                            {gapMid, dd, rdId});
                    }
                }
            }

            if(!splitEvents.empty()) {
                sort(splitEvents.begin(), splitEvents.end(),
                    [](const SplitDelEvent& a,
                       const SplitDelEvent& b) {
                        return a.refPos < b.refPos;
                    });
                // Cluster by position.
                struct SplitCluster {
                    uint32_t startPos, endPos;
                    vector<int64_t> sizes;
                    std::unordered_set<uint32_t> readIds;
                };
                vector<SplitCluster> clusters;
                SplitCluster cur;
                cur.startPos = splitEvents[0].refPos;
                cur.endPos = splitEvents[0].refPos;
                cur.sizes.push_back(splitEvents[0].delSize);
                cur.readIds.insert(splitEvents[0].readId);
                for(size_t i = 1;
                    i < splitEvents.size(); ++i) {
                    if(splitEvents[i].refPos
                       <= cur.endPos + 150) {
                        cur.endPos = splitEvents[i].refPos;
                        cur.sizes.push_back(
                            splitEvents[i].delSize);
                        cur.readIds.insert(
                            splitEvents[i].readId);
                    } else {
                        clusters.push_back(std::move(cur));
                        cur = SplitCluster();
                        cur.startPos =
                            splitEvents[i].refPos;
                        cur.endPos =
                            splitEvents[i].refPos;
                        cur.sizes.push_back(
                            splitEvents[i].delSize);
                        cur.readIds.insert(
                            splitEvents[i].readId);
                    }
                }
                clusters.push_back(std::move(cur));

                for(auto& cl : clusters) {
                    sort(cl.sizes.begin(), cl.sizes.end());
                    const int64_t medDel =
                        cl.sizes[cl.sizes.size() / 2];
                    const uint32_t bpPos =
                        (cl.startPos + cl.endPos) / 2;
                    if(medDel >= 20) {
                        cout << "    >>> DELETION CALL "
                             << "(early-split): "
                             << "size=" << medDel
                             << "bp, breakpoint=" << bpPos
                             << ", splitReads="
                             << cl.readIds.size()
                             << endl;
                        allDelCalls.push_back({
                            bpPos, medDel,
                            uint32_t(cl.readIds.size()),
                            "early-split"});
                    }
                }
            }
        };

        if(!bamFileName.empty()) {
            // Extract chromosome name and region offset from the
            // reference read name. The name may be "chr1:100-200".
            const auto refNameSpan =
                readsRef.getReadName(refId);
            string fullRefName(
                refNameSpan.begin(), refNameSpan.end());
            string refName = fullRefName;
            uint32_t regionStart = 0;
            const auto colonPos = fullRefName.find(':');
            if(colonPos != string::npos) {
                refName = fullRefName.substr(0, colonPos);
                const auto dashPos =
                    fullRefName.find('-', colonPos);
                if(dashPos != string::npos) {
                    regionStart = uint32_t(
                        stoul(fullRefName.substr(
                            colonPos + 1,
                            dashPos - colonPos - 1)));
                }
            }
            saTagCalls = parseSaTagSvCalls(
                bamFileName, refName,
                regionStart, regionStart + refLength);
            // Convert breakpoint positions from absolute to
            // relative to the reference subregion. Clamp
            // positions outside the region to the boundary
            // rather than discarding — SA-tag deletions
            // often have breakpoints outside the extracted
            // region because the aligner places the
            // supplementary alignment in flanking unique
            // sequence.
            {
                vector<SaTagSvCall> localCalls;
                for(auto& sc : saTagCalls) {
                    if(sc.refPos < regionStart)
                        sc.refPos = 0;
                    else if(sc.refPos > regionStart + refLength)
                        sc.refPos = refLength;
                    else
                        sc.refPos -= regionStart;
                    localCalls.push_back(sc);
                }
                saTagCalls = std::move(localCalls);
            }
            if(!saTagCalls.empty()) {
                cout << "    SA tag SV evidence ("
                     << saTagCalls.size() << " clusters):"
                     << endl;
                for(const auto& sc : saTagCalls) {
                    cout << "      " << sc.svType
                         << " " << sc.size << "bp"
                         << " at pos=" << sc.refPos
                         << " reads=" << sc.readCount
                         << endl;
                }
            }
            // Parse soft-clip breakpoints and CIGAR indels.
            parseBamEvidence(
                bamFileName, refName,
                regionStart, regionStart + refLength,
                softClipBPs, cigarIndels);

            if(!softClipBPs.empty()) {
                cout << "    Soft-clip breakpoints ("
                     << softClipBPs.size() << " clusters):"
                     << endl;
                for(const auto& sc : softClipBPs) {
                    cout << "      "
                         << (sc.isLeftClip ? "L" : "R")
                         << " pos=" << sc.refPos
                         << " reads=" << sc.readCount
                         << " avgClipLen=" << sc.avgClipLen
                         << endl;
                }
            }

            if(!cigarIndels.empty()) {
                cout << "    CIGAR indels ("
                     << cigarIndels.size() << " clusters):"
                     << endl;
                for(const auto& ci : cigarIndels) {
                    cout << "      " << ci.svType
                         << " " << ci.size << "bp"
                         << " at pos=" << ci.refPos
                         << " reads=" << ci.readCount;
                    if(!ci.insSeq.empty())
                        cout << " seq=" << ci.insSeq.size()
                             << "bp";
                    cout << endl;
                }
            }

            // Emit CIGAR indel calls with sufficient support.
            // DEL: >=2 reads (CIGAR D ops are high-confidence).
            // INS: >=3 reads (insertions need more support to
            // avoid false positives from alignment artifacts).
            for(const auto& ci : cigarIndels) {
                const uint32_t minReads =
                    (ci.svType == "DEL") ? 2 : 3;
                if(ci.readCount >= minReads && ci.size >= 30) {
                    cout << "    >>> "
                         << (ci.svType == "DEL"
                             ? "DELETION" : "INSERTION")
                         << " CALL (CIGAR): size="
                         << ci.size << "bp"
                         << ", breakpoint=" << ci.refPos
                         << ", reads=" << ci.readCount
                         << endl;
                    if(ci.svType == "DEL") {
                        allDelCalls.push_back({
                            ci.refPos, int64_t(ci.size),
                            ci.readCount, "CIGAR"});
                    }
                }
            }

            // Multi-k non-unique anchor deletion sizing.
            // Run on ALL reads using weighted diagonal histograms
            // from k-mers at multiple k values. Non-unique k-mers
            // contribute with weight inversely proportional to
            // their multiplicity.
            {
                vector<ReadId> allReadIds;
                const uint32_t totalReads =
                    uint32_t(readsRef.readCount());
                for(uint32_t ri = uint32_t(referenceReadCount);
                    ri < totalReads; ++ri) {
                    allReadIds.push_back(ReadId(ri));
                }
                const auto mkCalls = multiKAnchorSizing(
                    refId, allReadIds,
                    0, refLength,
                    uint32_t(assemblerInfo->k), 62);
                for(const auto& mk : mkCalls) {
                    if(mk.size >= 30 && mk.readCount >= 2) {
                        cout << "    >>> DELETION CALL"
                             << " (multi-k): size="
                             << mk.size << "bp"
                             << ", breakpoint=" << mk.breakpointPos
                             << ", reads=" << mk.readCount
                             << endl;
                        allDelCalls.push_back({
                            mk.breakpointPos,
                            mk.size,
                            mk.readCount,
                            "multi-k"});
                    }
                }
            }

            // Paired soft-clip INS sizing with de Bruijn assembly.
            // Right-clip reads end at the left breakpoint; their
            // clipped bases extend into the insertion. Left-clip
            // reads start at the right breakpoint; their clipped
            // bases extend into the insertion from the other side.
            // Assemble each side to get contig lengths, then
            // estimate insertion size from the sum minus overlap.
            for(const auto& rClip : softClipBPs) {
                if(rClip.isLeftClip) continue;
                if(rClip.readCount < 3) continue;
                for(const auto& lClip : softClipBPs) {
                    if(!lClip.isLeftClip) continue;
                    if(lClip.readCount < 3) continue;
                    const int64_t gap = int64_t(lClip.refPos)
                                      - int64_t(rClip.refPos);
                    if(gap < -10 || gap > 50) continue;

                    // Assemble each side's clip sequences.
                    const uint32_t rContigLen =
                        assembleClipSequences(rClip.clipSeqs);
                    const uint32_t lContigLen =
                        assembleClipSequences(lClip.clipSeqs);

                    // Estimate insertion size.
                    // If both sides assembled, the insertion is
                    // at least as long as the longer contig.
                    // If the contigs overlap (insertion < read
                    // length), the true size is captured by the
                    // overlap. For larger insertions, sum the
                    // contig lengths minus estimated overlap.
                    int64_t insSize;
                    if(rContigLen > 0 && lContigLen > 0) {
                        // Use assembled contig lengths.
                        // Subtract k-1 overlap estimate.
                        insSize = int64_t(rContigLen)
                                + int64_t(lContigLen) - 20;
                        // But at least the max of the two.
                        insSize = std::max(insSize,
                            int64_t(std::max(
                                rContigLen, lContigLen)));
                    } else if(rContigLen > 0 || lContigLen > 0) {
                        // One side assembled.
                        insSize = int64_t(std::max(
                            rContigLen, lContigLen));
                    } else {
                        // Fallback to average clip lengths.
                        insSize = int64_t(rClip.avgClipLen)
                                + int64_t(lClip.avgClipLen)
                                - std::max(int64_t(0), gap);
                    }

                    if(insSize >= 50 && insSize <= 10000) {
                        const uint32_t bpPos =
                            (rClip.refPos + lClip.refPos) / 2;
                        cout << "    >>> INSERTION CALL"
                             << " (soft-clip assembly): size="
                             << insSize << "bp"
                             << ", breakpoint=" << bpPos
                             << ", Rclip=" << rClip.readCount
                             << "reads"
                             << " contig=" << rContigLen
                             << "bp"
                             << ", Lclip=" << lClip.readCount
                             << "reads"
                             << " contig=" << lContigLen
                             << "bp"
                             << endl;
                    }
                }
            }
        }

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

        if(allRefOrdinals.size() < 2) {
            // Add SA-tag DEL calls.
            for(const auto& sc : saTagCalls) {
                if(sc.size >= 50)
                    allDelCalls.push_back({sc.refPos, sc.size,
                        sc.readCount, "SA-" + sc.svType});
            }
            detectSplitReadDels();
            continue;
        }

        // Build a map from reference ordinal -> boundary index.
        unordered_map<uint32_t, uint32_t> ordinalToBoundary;
        ordinalToBoundary.reserve(allRefOrdinals.size());
        for(uint32_t bi = 0; bi < uint32_t(allRefOrdinals.size()); ++bi) {
            ordinalToBoundary[allRefOrdinals[bi]] = bi;
        }

        const uint32_t nBoundaries = uint32_t(allRefOrdinals.size());
        const uint32_t nSegments = nBoundaries - 1;

        // -----------------------------------------------------------------
        // Flag set during per-segment processing to suppress
        // SA-tag DEL calls in marker-depleted VNTR regions.
        bool suppressSaTagDel = false;

        // Coverage-drop regions detected during per-segment
        // processing, used later for k-mer cluster corroboration.
        struct CovDropRegion {
            uint32_t startPos;
            uint32_t endPos;
            bool markerDepleted;
        };
        vector<CovDropRegion> covDropRegions;

        // DEL calls from diagonal-shift and split-read analyses,
        // stored for post-coverage-drop INS type-flip check.
        vector<DelCallRecord> delCallRecords;

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

        if(badSegment || segmentStrings.empty()) {
            // Add SA-tag DEL calls.
            for(const auto& sc : saTagCalls) {
                if(sc.size >= 50)
                    allDelCalls.push_back({sc.refPos, sc.size,
                        sc.readCount, "SA-" + sc.svType});
            }
            detectSplitReadDels();
            continue;
        }

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

        // Emit individual per-read DEL and INS detections into
        // allDelCalls as direct evidence calls.
        for(const auto& rg : readGroups) {
            if(rg.svType == SvType::Deletion && rg.svSize >= 20) {
                allDelCalls.push_back({
                    rg.breakpointRefPos, rg.svSize,
                    1, "per-read-DEL"});
            }
            if(rg.svType == SvType::Insertion && rg.svSize >= 20) {
                allDelCalls.push_back({
                    rg.breakpointRefPos, rg.svSize,
                    1, "per-read-INS-flip"});
                allInsCalls.push_back({
                    rg.breakpointRefPos, rg.svSize,
                    1, "per-read-INS"});
            }
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
        unordered_set<uint32_t> indirectAlignedReads;
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
                    indirectAlignedReads.insert(edge.neighborReadId);

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
                    uint32_t actualRefPos; // precise chain endpoint position
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
                        chainEndReads[endWin].push_back(
                            {rid, rightOvhBp, info.maxRefPos});
                    }
                    if(leftOvhBp >= minOvhBp && leftOvhOrd >= 2) {
                        chainStartReads[startWin].push_back(
                            {rid, leftOvhBp, info.minRefPos});
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
                    const uint32_t windowCenter = refStartPos + w * windowSize + windowSize / 2;

                    // Left breakpoint: many chain ends with right overhang.
                    if(chainEndCount[w] > 0) {
                        const double fold = chainEndCount[w] / std::max(bgEndRate, 0.1);
                        auto it = chainEndReads.find(w);
                        const uint32_t ovhCount = (it != chainEndReads.end())
                            ? uint32_t(it->second.size()) : 0;

                        if(fold >= minFoldEnrichment && ovhCount >= minEndpointReads) {
                            // Use median of actual read endpoint
                            // positions for sub-window precision.
                            uint32_t refPos = windowCenter;
                            if(it != chainEndReads.end()
                               && !it->second.empty()) {
                                vector<uint32_t> actPos;
                                for(const auto& ep :
                                    it->second) {
                                    actPos.push_back(
                                        ep.actualRefPos);
                                }
                                sort(actPos.begin(),
                                     actPos.end());
                                refPos = actPos[
                                    actPos.size() / 2];
                            }
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
                            uint32_t refPos = windowCenter;
                            if(it != chainStartReads.end()
                               && !it->second.empty()) {
                                vector<uint32_t> actPos;
                                for(const auto& ep :
                                    it->second) {
                                    actPos.push_back(
                                        ep.actualRefPos);
                                }
                                sort(actPos.begin(),
                                     actPos.end());
                                refPos = actPos[
                                    actPos.size() / 2];
                            }
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

                // Flag marker-depleted VNTR regions for SA-tag
                // DEL suppression. Triggers when many hit-depth
                // BPs exist with many unanchored reads. VNTR gaps
                // detected during breakpoint pairing also set this
                // flag (see vntrGaps.push_back below).
                if(unanchoredReads.size() >= 30
                   && hitDepthBreakpoints.size() >= 10) {
                    suppressSaTagDel = true;
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

                // Generate breakpoints from hit-depth cluster edges
                // when chain-endpoint breakpoints are missing on one
                // or both sides. For insertions where chains don't
                // break (the DP bridges the insertion), the only
                // signal is a contiguous hit-depth drop. The left
                // edge of the drop cluster is a left BP and the right
                // edge is a right BP. Use reads spanning across the
                // drop zone as overhang reads.
                if(!hitDepthBreakpoints.empty()
                   && (leftBreakpoints.empty()
                       || rightBreakpoints.empty())) {
                    // Cluster consecutive hit-depth BPs (within 2 windows).
                    struct HdCluster {
                        uint32_t startWin;
                        uint32_t endWin;
                        uint32_t startPos;
                        uint32_t endPos;
                        double minRatio;
                    };
                    vector<HdCluster> hdClusters;
                    HdCluster curHd;
                    curHd.startWin = hitDepthBreakpoints[0].windowIdx;
                    curHd.endWin = hitDepthBreakpoints[0].windowIdx;
                    curHd.startPos = hitDepthBreakpoints[0].refPos;
                    curHd.endPos = hitDepthBreakpoints[0].refPos;
                    curHd.minRatio = hitDepthBreakpoints[0].dropRatio;

                    for(size_t i = 1; i < hitDepthBreakpoints.size(); ++i) {
                        if(hitDepthBreakpoints[i].windowIdx
                           <= curHd.endWin + 2) {
                            curHd.endWin = hitDepthBreakpoints[i].windowIdx;
                            curHd.endPos = hitDepthBreakpoints[i].refPos;
                            curHd.minRatio = std::min(
                                curHd.minRatio,
                                hitDepthBreakpoints[i].dropRatio);
                        } else {
                            hdClusters.push_back(curHd);
                            curHd.startWin = hitDepthBreakpoints[i].windowIdx;
                            curHd.endWin = hitDepthBreakpoints[i].windowIdx;
                            curHd.startPos = hitDepthBreakpoints[i].refPos;
                            curHd.endPos = hitDepthBreakpoints[i].refPos;
                            curHd.minRatio = hitDepthBreakpoints[i].dropRatio;
                        }
                    }
                    hdClusters.push_back(curHd);

                    // For each cluster with strong depth drop, generate
                    // left BP at the start and right BP at the end.
                    for(const auto& hdc : hdClusters) {
                        if(hdc.minRatio > 0.3) continue; // Weak drop.
                        const uint32_t clusterSpan =
                            hdc.endPos - hdc.startPos + windowSize;
                        if(clusterSpan < windowSize) continue;

                        // Collect reads spanning across the cluster
                        // as overhang reads for both BPs.
                        vector<EndpointInfo> leftOvhReads;
                        vector<EndpointInfo> rightOvhReads;
                        for(const auto& [rid, info] : readChainInfoMap) {
                            if(info.totalMarkers < 2) continue;
                            // Read must span past the cluster edge.
                            if(info.maxRefPos > hdc.startPos + windowSize
                               && info.minRefPos < hdc.startPos) {
                                const auto rdMkrs = markersRef[
                                    OrientedReadId(ReadId(rid), 0).getValue()];
                                if(rdMkrs.size() < 2) continue;
                                const int64_t ovh = int64_t(
                                    rdMkrs[rdMkrs.size()-1].position)
                                    - int64_t(rdMkrs[info.maxReadOrd].position);
                                if(ovh >= 20) {
                                    leftOvhReads.push_back(
                                        {rid, ovh});
                                }
                            }
                            if(info.minRefPos < hdc.endPos - windowSize
                               && info.maxRefPos > hdc.endPos) {
                                const auto rdMkrs = markersRef[
                                    OrientedReadId(ReadId(rid), 0).getValue()];
                                if(rdMkrs.size() < 2) continue;
                                const int64_t ovh = int64_t(
                                    rdMkrs[info.minReadOrd].position)
                                    - int64_t(rdMkrs[0].position);
                                if(ovh >= 20) {
                                    rightOvhReads.push_back(
                                        {rid, ovh});
                                }
                            }
                        }

                        if(!leftOvhReads.empty()
                           && leftBreakpoints.empty()) {
                            const double fold = hdc.minRatio > 0
                                ? 1.0 / hdc.minRatio : 10.0;
                            leftBreakpoints.push_back({
                                hdc.startWin, hdc.startPos,
                                uint32_t(leftOvhReads.size()),
                                spanningCount[hdc.startWin],
                                uint32_t(leftOvhReads.size()),
                                fold,
                                leftOvhReads
                            });
                            cout << "      HitDepth-generated Left BP:"
                                 << " pos=" << hdc.startPos
                                 << " ovhReads="
                                 << leftOvhReads.size()
                                 << " fold=" << fold
                                 << endl;
                        }
                        // If no spanning reads found for right BP,
                        // look for reads whose chains start just
                        // after the drop zone with left overhang
                        // (indicating they extend back into the
                        // insertion from the right flank).
                        if(rightOvhReads.empty()) {
                            for(const auto& [rid, info] : readChainInfoMap) {
                                if(info.totalMarkers < 2) continue;
                                // Read starts near the right edge
                                // of the drop zone.
                                if(info.minRefPos >= hdc.endPos - windowSize
                                   && info.minRefPos <= hdc.endPos + windowSize * 3) {
                                    const auto rdMkrs = markersRef[
                                        OrientedReadId(ReadId(rid), 0).getValue()];
                                    if(rdMkrs.size() < 2) continue;
                                    const int64_t ovh = int64_t(
                                        rdMkrs[info.minReadOrd].position)
                                        - int64_t(rdMkrs[0].position);
                                    if(ovh >= 20) {
                                        rightOvhReads.push_back(
                                            {rid, ovh});
                                    }
                                }
                            }
                        }

                        if(!rightOvhReads.empty()
                           && rightBreakpoints.empty()) {
                            const double fold = hdc.minRatio > 0
                                ? 1.0 / hdc.minRatio : 10.0;
                            rightBreakpoints.push_back({
                                hdc.endWin, hdc.endPos,
                                uint32_t(rightOvhReads.size()),
                                spanningCount[hdc.endWin],
                                uint32_t(rightOvhReads.size()),
                                fold,
                                rightOvhReads
                            });
                            cout << "      HitDepth-generated Right BP:"
                                 << " pos=" << hdc.endPos
                                 << " ovhReads="
                                 << rightOvhReads.size()
                                 << " fold=" << fold
                                 << endl;
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

                // For each left breakpoint, find the best right breakpoint
                // to form a breakpoint pair. Start with the nearest BP
                // (preserves original behavior for cases where proximity
                // is the right signal). Build a ranked candidate list
                // so that if the nearest fails to produce a path, we
                // can fall through to stronger alternatives.
                for(const auto& lbp : leftBreakpoints) {
                    // Build candidate list sorted by distance (nearest first).
                    struct RbpCandidate {
                        const Breakpoint* rbp;
                        int64_t dist;
                        double score;
                    };
                    vector<RbpCandidate> rbpCandidates;
                    for(const auto& rbp : rightBreakpoints) {
                        const int64_t dist = std::abs(
                            int64_t(rbp.refPos) - int64_t(lbp.refPos));
                        // Score for fallback ranking: strength / distance.
                        const double strength =
                            rbp.foldEnrichment
                            * double(std::max(rbp.ovhReadCount, 1u))
                            * double(std::max(rbp.ovhReadCount, 1u));
                        const double score =
                            strength / (1.0 + double(dist) / 100.0);
                        rbpCandidates.push_back({&rbp, dist, score});
                    }
                    // Sort by distance ascending (nearest first).
                    sort(rbpCandidates.begin(), rbpCandidates.end(),
                        [](const RbpCandidate& a, const RbpCandidate& b) {
                            return a.dist < b.dist;
                        });

                    const Breakpoint* bestRbp = rbpCandidates.empty()
                        ? nullptr : rbpCandidates[0].rbp;
                    int64_t bestDist = rbpCandidates.empty()
                        ? INT64_MAX : rbpCandidates[0].dist;

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
                            suppressSaTagDel = true;
                            cout << "    VNTR gap: L=" << lbp.refPos
                                 << " R=" << bestRbp->refPos
                                 << " lowDepth=" << lowDepthWindows
                                 << "/" << totalGapWindows
                                 << " — extending pair distance"
                                 << endl;
                        }
                        // Allow pairing when both BPs have strong
                        // signals, even without a hit-depth drop.
                        if(bestDist > maxPairDist
                           && lbp.foldEnrichment >= 3.0
                           && bestRbp->foldEnrichment >= 3.0
                           && lbp.ovhReadCount >= 3
                           && bestRbp->ovhReadCount >= 3) {
                            maxPairDist = bestDist + 1;
                            vntrGaps.push_back({lbp.refPos, bestRbp->refPos});
                            cout << "    Strong BP pair: L="
                                 << lbp.refPos
                                 << " (fold=" << lbp.foldEnrichment
                                 << " ovh=" << lbp.ovhReadCount
                                 << ") R=" << bestRbp->refPos
                                 << " (fold=" << bestRbp->foldEnrichment
                                 << " ovh=" << bestRbp->ovhReadCount
                                 << ") dist=" << bestDist
                                 << endl;
                        }
                    }

                    if(!bestRbp || bestDist > maxPairDist) {
                        // If the nearest candidate is out of range,
                        // try remaining candidates by score (strongest first).
                        // Sort remaining candidates by score descending.
                        sort(rbpCandidates.begin(), rbpCandidates.end(),
                            [](const RbpCandidate& a, const RbpCandidate& b) {
                                return a.score > b.score;
                            });
                        bool foundCandidate = false;
                        for(size_t ci = 0; ci < rbpCandidates.size(); ++ci) {
                            if(rbpCandidates[ci].rbp == bestRbp) continue;
                            int64_t candMaxDist = 500;
                            const auto* candRbp = rbpCandidates[ci].rbp;
                            const int64_t candDist = rbpCandidates[ci].dist;
                            // Check VNTR for this candidate.
                            if(candDist > 500) {
                                const uint32_t gs = (lbp.refPos - refStartPos) / windowSize;
                                const uint32_t ge = (candRbp->refPos - refStartPos) / windowSize;
                                uint32_t ldw = 0, tgw = 0;
                                for(uint32_t w = gs; w <= ge && w < nWindows; ++w) {
                                    if(windowMarkerCount[w] == 0) continue;
                                    ++tgw;
                                    if(medianHitDepth > 0
                                       && windowHitDepth[w] / medianHitDepth < hitDepthDropThreshold)
                                        ++ldw;
                                }
                                if(tgw > 0 && double(ldw)/double(tgw) > 0.5) {
                                    candMaxDist = candDist + 1;
                                    vntrGaps.push_back({lbp.refPos, candRbp->refPos});
                                }
                                // Also allow pairing when both BPs
                                // have strong signals (high fold and
                                // many overhang reads), even without
                                // a hit-depth drop. This handles
                                // tandem repeat deletions with
                                // uniform coverage.
                                if(candDist <= candMaxDist) {
                                    // Already accepted above.
                                } else if(lbp.foldEnrichment >= 3.0
                                          && candRbp->foldEnrichment >= 3.0
                                          && lbp.ovhReadCount >= 3
                                          && candRbp->ovhReadCount >= 3) {
                                    candMaxDist = candDist + 1;
                                    vntrGaps.push_back({lbp.refPos, candRbp->refPos});
                                    cout << "    Strong BP pair: L="
                                         << lbp.refPos
                                         << " (fold=" << lbp.foldEnrichment
                                         << " ovh=" << lbp.ovhReadCount
                                         << ") R=" << candRbp->refPos
                                         << " (fold=" << candRbp->foldEnrichment
                                         << " ovh=" << candRbp->ovhReadCount
                                         << ") dist=" << candDist
                                         << endl;
                                }
                            }
                            if(candDist <= candMaxDist) {
                                bestRbp = candRbp;
                                bestDist = candDist;
                                maxPairDist = candMaxDist;
                                foundCandidate = true;
                                break;
                            }
                        }
                        if(!foundCandidate) continue;
                    }

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

                    // Detect marker-depleted regions around the breakpoint.
                    // Even when the left/right BPs are close (refGap < 500),
                    // the surrounding region may be a VNTR where read-to-read
                    // chains create false shortcuts. Check for windows with
                    // zero markers (fully depleted) in a ±300bp window.
                    // Only flag as VNTR if >60% of windows are fully depleted
                    // (not just low depth — low depth can be caused by the
                    // insertion itself in non-repetitive regions).
                    bool isVntrGap = (maxPairDist > 500);
                    if(!isVntrGap) {
                        const uint32_t bpMid = (lbp.refPos + bestRbp->refPos) / 2;
                        const uint32_t checkRadius = 300;
                        const uint32_t checkStart = (bpMid > refStartPos + checkRadius)
                            ? (bpMid - refStartPos - checkRadius) / windowSize : 0;
                        const uint32_t checkEnd = std::min(
                            (bpMid - refStartPos + checkRadius) / windowSize, nWindows - 1);
                        uint32_t emptyWins = 0, totalWins = 0;
                        for(uint32_t w = checkStart; w <= checkEnd; ++w) {
                            ++totalWins;
                            if(windowMarkerCount[w] == 0) {
                                ++emptyWins;
                            }
                        }
                        if(totalWins > 0
                           && double(emptyWins) / double(totalWins) > 0.6) {
                            isVntrGap = true;
                        }
                    }

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
                        // Negative values indicate a deletion in the
                        // VNTR: the sample has fewer repeat copies than
                        // the reference. Emit both hom and het estimates
                        // so the adj/mult/div machinery can find the
                        // best match.
                        else if(insLenHom < -30) {
                            const int64_t delSize = std::abs(insLenHom);
                            cout << "    >>> DELETION CALL"
                                 << " (VNTR-depth-hom): size="
                                 << delSize << "bp"
                                 << ", breakpoint="
                                 << breakpointPos
                                 << ", refLen=" << vntrRefLen
                                 << ", flankCov="
                                 << flankCoverage
                                 << endl;
                            allDelCalls.push_back({
                                breakpointPos, delSize,
                                0, "VNTR-depth"});
                            // Also emit the het estimate.
                            if(insLenHet < -30) {
                                const int64_t delSizeHet =
                                    std::abs(insLenHet);
                                cout << "    >>> DELETION CALL"
                                     << " (VNTR-depth-het): size="
                                     << delSizeHet << "bp"
                                     << ", breakpoint="
                                     << breakpointPos
                                     << ", refLen=" << vntrRefLen
                                     << ", flankCov="
                                     << flankCoverage
                                     << endl;
                                allDelCalls.push_back({
                                    breakpointPos, delSizeHet,
                                    0, "VNTR-depth"});
                            }
                            // Also emit refGap as a DEL call —
                            // the breakpoint distance itself is
                            // a size estimate for the deletion.
                            if(vntrRefLen >= 30) {
                                cout << "    >>> DELETION CALL"
                                     << " (VNTR-refGap): size="
                                     << vntrRefLen << "bp"
                                     << ", breakpoint="
                                     << breakpointPos
                                     << endl;
                                allDelCalls.push_back({
                                    breakpointPos, vntrRefLen,
                                    0, "VNTR-refGap"});
                            }
                        } else if(insLenHet < -30) {
                            const int64_t delSize =
                                std::abs(insLenHet);
                            cout << "    >>> DELETION CALL"
                                 << " (VNTR-depth-het): size="
                                 << delSize << "bp"
                                 << ", breakpoint="
                                 << breakpointPos
                                 << ", refLen=" << vntrRefLen
                                 << ", flankCov="
                                 << flankCoverage
                                 << endl;
                            allDelCalls.push_back({
                                breakpointPos, delSize,
                                0, "VNTR-depth"});
                            if(vntrRefLen >= 30) {
                                cout << "    >>> DELETION CALL"
                                     << " (VNTR-refGap): size="
                                     << vntrRefLen << "bp"
                                     << ", breakpoint="
                                     << breakpointPos
                                     << endl;
                                allDelCalls.push_back({
                                    breakpointPos, vntrRefLen,
                                    0, "VNTR-refGap"});
                            }
                        }
                    }

                    // When path-finding and VNTR depth both fail,
                    // emit a DEL call using refGap as size estimate.
                    // The breakpoint pair marks the deletion boundaries
                    // even when chains can't span the repeat.
                    if(!foundPath && bestDist >= 50) {
                        cout << "    >>> DELETION CALL"
                             << " (bp-pair-nofp): size="
                             << bestDist << "bp"
                             << ", breakpoint="
                             << breakpointPos
                             << endl;
                        allDelCalls.push_back({
                            breakpointPos, int64_t(bestDist),
                            uint32_t(lbp.endpointCount
                                     + bestRbp->endpointCount),
                            "bp-pair-nofp"});
                    }

                    if(foundPath && bestPathDist > 20) {
                        // Use the maximum path distance (bestPathDist).
                        // The max represents the longest read-space
                        // traversal, which best estimates the full
                        // insertion size.
                        const int64_t reportedPathDist = bestPathDist;

                        // Determine if this is an insertion or deletion.
                        //
                        // For a deletion with refGap > pathDist, the left
                        // breakpoint should have a strong endpoint signal
                        // (many chain ends) because reads stop at the
                        // deletion boundary. Require the left BP to have
                        // high fold enrichment AND many endpoints.
                        //
                        // For an insertion, pathDist directly gives the
                        // insertion size. When pathDist >> refGap, subtract
                        // refGap to account for reference traversal in the
                        // overhangs.
                        bool isDeletion = false;
                        if(reportedPathDist < int64_t(bestDist)
                           && bestDist > 100
                           && lbp.endpointCount >= 10
                           && lbp.foldEnrichment >= 3.0) {
                            isDeletion = true;
                        }

                        if(isDeletion) {
                            const int64_t delSz = bestDist - reportedPathDist;
                            if(delSz > 20) {
                                cout << "    >>> DELETION CALL (path-based): "
                                     << "size=" << delSz << "bp, "
                                     << "breakpoint=" << breakpointPos << ", "
                                     << "leftEnds=" << lbp.endpointCount << ", "
                                     << "rightStarts=" << bestRbp->endpointCount << ", "
                                     << "hops=" << bestPathLen
                                     << " (pathDist=" << reportedPathDist
                                     << " refGap=" << bestDist << ")"
                                     << endl;
                                if(delSz >= 50) {
                                    allDelCalls.push_back({
                                        breakpointPos, delSz,
                                        uint32_t(
                                            lbp.endpointCount
                                            + bestRbp->endpointCount),
                                        "path-based"});
                                }
                            }
                        } else {
                            // Insertion sizing. The path distance
                            // includes reference traversal in the
                            // overhangs when reads span from the chain
                            // endpoint into the insertion. Subtract
                            // refGap when it's substantial (> 150bp)
                            // and pathDist clearly exceeds it. Small
                            // refGaps (≤ 150bp) are window quantization
                            // noise and should not be subtracted.
                            int64_t insSz = reportedPathDist;
                            if(bestDist > 150
                               && reportedPathDist > int64_t(bestDist)) {
                                insSz = reportedPathDist - bestDist;
                            }
                            if(insSz > 20) {
                                // Before emitting, check if a further
                                // right BP produces a larger insertion.
                                // This handles cases where the nearest
                                // BP captures only a partial insertion
                                // (insertion larger than read length).
                                int64_t bestInsSz = insSz;
                                const Breakpoint* bestInsRbp = bestRbp;
                                int64_t bestInsPathDist = reportedPathDist;
                                int bestInsHops = bestPathLen;
                                int64_t bestInsRefGap = bestDist;

                                for(const auto& cand : rbpCandidates) {
                                    if(cand.rbp == bestRbp) continue;
                                    if(cand.dist > uint64_t(maxPairDist)) continue;
                                    if(cand.rbp->refPos <= bestRbp->refPos) continue;

                                    // Try path search with this candidate.
                                    unordered_set<uint32_t> altRightIds;
                                    for(const auto& ei : cand.rbp->reads)
                                        altRightIds.insert(ei.readId);

                                    int64_t altBestPathDist = 0;
                                    int altBestPathLen = 0;
                                    bool altFoundPath = false;

                                    if(!isVntrGap)
                                    for(const auto& lf : lbp.reads) {
                                        auto gL2 = readGraph.find(lf.readId);
                                        if(gL2 == readGraph.end()) continue;
                                        for(const auto& e1 : gL2->second) {
                                            const uint32_t n1 = e1.neighborReadId;
                                            const int64_t ovlp1 = getOverlapSpan(e1);
                                            if(ovlp1 < 0) continue;
                                            if(altRightIds.count(n1)) {
                                                int64_t rOvh = 0;
                                                for(const auto& rf : cand.rbp->reads)
                                                    if(rf.readId == n1) { rOvh = rf.overhangBp; break; }
                                                const int64_t d = lf.overhangBp + rOvh - ovlp1;
                                                if(d > altBestPathDist) {
                                                    altBestPathDist = d;
                                                    altBestPathLen = 1;
                                                    altFoundPath = true;
                                                }
                                                continue;
                                            }
                                            if(!unanchoredReads.count(n1)) continue;
                                            const auto n1M = markersRef[
                                                OrientedReadId(ReadId(n1), 0).getValue()];
                                            if(n1M.size() < 2) continue;
                                            const int64_t n1S = int64_t(n1M[n1M.size()-1].position)
                                                               - int64_t(n1M[0].position);
                                            auto gN1 = readGraph.find(n1);
                                            if(gN1 == readGraph.end()) continue;
                                            for(const auto& e2 : gN1->second) {
                                                const uint32_t n2 = e2.neighborReadId;
                                                if(n2 == lf.readId) continue;
                                                const int64_t ovlp2 = getOverlapSpan(e2);
                                                if(ovlp2 < 0) continue;
                                                if(altRightIds.count(n2)) {
                                                    int64_t rOvh = 0;
                                                    for(const auto& rf : cand.rbp->reads)
                                                        if(rf.readId == n2) { rOvh = rf.overhangBp; break; }
                                                    const int64_t d = lf.overhangBp + n1S - ovlp1 - ovlp2 + rOvh;
                                                    if(d > altBestPathDist) {
                                                        altBestPathDist = d;
                                                        altBestPathLen = 2;
                                                        altFoundPath = true;
                                                    }
                                                    continue;
                                                }
                                                if(!unanchoredReads.count(n2)) continue;
                                                const auto n2M = markersRef[
                                                    OrientedReadId(ReadId(n2), 0).getValue()];
                                                if(n2M.size() < 2) continue;
                                                const int64_t n2S = int64_t(n2M[n2M.size()-1].position)
                                                                   - int64_t(n2M[0].position);
                                                auto gN2 = readGraph.find(n2);
                                                if(gN2 == readGraph.end()) continue;
                                                for(const auto& e3 : gN2->second) {
                                                    const uint32_t n3 = e3.neighborReadId;
                                                    if(n3 == n1 || n3 == lf.readId) continue;
                                                    if(!altRightIds.count(n3)) continue;
                                                    const int64_t ovlp3 = getOverlapSpan(e3);
                                                    if(ovlp3 < 0) continue;
                                                    int64_t rOvh = 0;
                                                    for(const auto& rf : cand.rbp->reads)
                                                        if(rf.readId == n3) { rOvh = rf.overhangBp; break; }
                                                    const int64_t d = lf.overhangBp + n1S + n2S
                                                        - ovlp1 - ovlp2 - ovlp3 + rOvh;
                                                    if(d > altBestPathDist) {
                                                        altBestPathDist = d;
                                                        altBestPathLen = 3;
                                                        altFoundPath = true;
                                                    }
                                                }
                                            }
                                        }
                                    }

                                    if(altFoundPath) {
                                        int64_t altInsSz = altBestPathDist;
                                        if(cand.dist > 150
                                           && altBestPathDist > int64_t(cand.dist)) {
                                            altInsSz = altBestPathDist - cand.dist;
                                        }
                                        if(altInsSz > bestInsSz) {
                                            bestInsSz = altInsSz;
                                            bestInsRbp = cand.rbp;
                                            bestInsPathDist = altBestPathDist;
                                            bestInsHops = altBestPathLen;
                                            bestInsRefGap = cand.dist;
                                        }
                                    }
                                }

                                // Depth-based fallback: when the path-
                                // based call is small but a further right
                                // BP exists with a hit-depth drop zone
                                // between L and R, estimate insertion
                                // size from unanchored read bases.
                                if(bestInsSz < 200) {
                                    for(const auto& cand : rbpCandidates) {
                                        if(cand.rbp == bestInsRbp) continue;
                                        if(cand.dist > uint64_t(maxPairDist)) continue;
                                        if(cand.rbp->refPos <= bestInsRbp->refPos) continue;
                                        if(cand.dist < 200) continue;

                                        // Check hit-depth drop between
                                        // L and this candidate R.
                                        const uint32_t zs = (lbp.refPos - refStartPos) / windowSize;
                                        const uint32_t ze = (cand.rbp->refPos - refStartPos) / windowSize;
                                        uint32_t lowWins = 0, totWins = 0;
                                        for(uint32_t w = zs; w <= ze && w < nWindows; ++w) {
                                            if(windowMarkerCount[w] == 0) continue;
                                            ++totWins;
                                            if(medianHitDepth > 0
                                               && windowHitDepth[w] / medianHitDepth < hitDepthDropThreshold)
                                                ++lowWins;
                                        }
                                        if(totWins == 0 || double(lowWins) / double(totWins) < 0.3)
                                            continue;

                                        // Estimate insertion size from
                                        // the hit-depth drop zone span.
                                        // The drop zone is the contiguous
                                        // region of low hit-depth between
                                        // L and R. The insertion disrupts
                                        // k-mer matches over this span.
                                        uint32_t dropStart = UINT32_MAX;
                                        uint32_t dropEnd = 0;
                                        for(uint32_t w = zs; w <= ze && w < nWindows; ++w) {
                                            if(windowMarkerCount[w] == 0) continue;
                                            if(medianHitDepth > 0
                                               && windowHitDepth[w] / medianHitDepth < hitDepthDropThreshold) {
                                                const uint32_t wp = refStartPos + w * windowSize;
                                                if(wp < dropStart) dropStart = wp;
                                                if(wp > dropEnd) dropEnd = wp;
                                            }
                                        }
                                        if(dropStart < dropEnd) {
                                            const int64_t dropSpan =
                                                int64_t(dropEnd - dropStart) + int64_t(windowSize);
                                            if(dropSpan > bestInsSz && dropSpan >= 50) {
                                                bestInsSz = dropSpan;
                                                bestInsRbp = cand.rbp;
                                                bestInsPathDist = dropSpan;
                                                bestInsHops = 0;
                                                bestInsRefGap = cand.dist;
                                                cout << "    HitDepth-span INS estimate:"
                                                     << " dropZone=" << dropStart
                                                     << "-" << dropEnd
                                                     << " span=" << dropSpan
                                                     << "bp" << endl;
                                            }
                                        }
                                        // Continue to try further candidates
                                        // for a larger drop span.
                                    }
                                }

                                const uint32_t insBpPos =
                                    (lbp.refPos + bestInsRbp->refPos) / 2;

                                // In highly repetitive regions (many
                                // BPs on both sides), path-based INS
                                // calls are unreliable because the
                                // read graph has many false connections
                                // through rescued k-mers. Suppress
                                // unless both BPs have very strong
                                // overhang support.
                                const bool highlyRepetitive =
                                    leftBreakpoints.size() >= 5
                                    && rightBreakpoints.size() >= 5;
                                if(highlyRepetitive
                                   && (lbp.ovhReadCount < 20
                                       || bestInsRbp->ovhReadCount < 20)) {
                                    cout << "    INS suppressed"
                                         << " (highly repetitive,"
                                         << " " << leftBreakpoints.size()
                                         << "L/" << rightBreakpoints.size()
                                         << "R BPs): size="
                                         << bestInsSz << "bp"
                                         << ", breakpoint="
                                         << insBpPos << endl;
                                } else {
                                    cout << "    >>> INSERTION CALL: "
                                         << "size=" << bestInsSz << "bp, "
                                         << "breakpoint=" << insBpPos << ", "
                                         << "leftEnds=" << lbp.endpointCount << ", "
                                         << "rightStarts=" << bestInsRbp->endpointCount << ", "
                                         << "hops=" << bestInsHops
                                         << " (pathDist=" << bestInsPathDist
                                         << " refGap=" << bestInsRefGap << ")"
                                         << endl;
                                    insertionCallRegions.push_back({
                                        std::min(lbp.refPos, bestInsRbp->refPos),
                                        std::max(lbp.refPos, bestInsRbp->refPos)
                                    });
                                    allInsCalls.push_back({
                                        insBpPos,
                                        bestInsSz,
                                        uint32_t(
                                            lbp.endpointCount
                                            + bestInsRbp
                                                ->endpointCount),
                                        "path-ins"});
                                    // In tandem repeats, the path
                                    // traverses repeat units from
                                    // the non-deleted allele, so
                                    // pathDist reflects insertion
                                    // while the truth may be a
                                    // deletion. When refGap > 0,
                                    // also emit a DEL call: refGap
                                    // approximates the deletion
                                    // size (distance between BPs
                                    // on the reference).
                                    if(bestInsRefGap >= 40) {
                                        // Refine using soft-clip BPs.
                                        int64_t pmSize = bestInsRefGap;
                                        if(!softClipBPs.empty()) {
                                            const int64_t sr =
                                                int64_t(windowSize);
                                            int64_t dL = INT64_MAX;
                                            uint32_t pL = lbp.refPos;
                                            for(const auto& sc : softClipBPs) {
                                                if(sc.isLeftClip) continue;
                                                if(sc.readCount < 2) continue;
                                                const int64_t d = std::abs(
                                                    int64_t(sc.refPos)
                                                    - int64_t(lbp.refPos));
                                                if(d < dL) { dL = d; pL = sc.refPos; }
                                            }
                                            int64_t dR = INT64_MAX;
                                            uint32_t pR = bestInsRbp->refPos;
                                            for(const auto& sc : softClipBPs) {
                                                if(!sc.isLeftClip) continue;
                                                if(sc.readCount < 2) continue;
                                                const int64_t d = std::abs(
                                                    int64_t(sc.refPos)
                                                    - int64_t(bestInsRbp->refPos));
                                                if(d < dR) { dR = d; pR = sc.refPos; }
                                            }
                                            if(dL <= sr && dR <= sr
                                               && pR > pL) {
                                                const int64_t scD =
                                                    int64_t(pR) - int64_t(pL);
                                                if(scD >= 40) {
                                                    pmSize = scD;
                                                }
                                            }
                                        }
                                        cout << "    >>> DELETION CALL"
                                             << " (path-mirror):"
                                             << " size="
                                             << pmSize
                                             << "bp, breakpoint="
                                             << insBpPos
                                             << endl;
                                        delCallRecords.push_back({
                                            insBpPos,
                                            pmSize,
                                            uint32_t(
                                                lbp.endpointCount
                                                + bestInsRbp
                                                  ->endpointCount),
                                            "path-mirror"});


                                    }
                                }
                            }
                        }
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

                        uint32_t nChainsChecked = 0;
                        uint32_t nChainsSpanning = 0;

                        for(const auto& ce : chainsForRef) {
                            if(ce.readId == uint32_t(refId)) continue;
                            const auto& al = alignments[ce.chainIndex];
                            if(al.ordinals.size() < 2) continue;
                            ++nChainsChecked;

                            // Check if chain spans the deletion zone.
                            bool hasLeft = false, hasRight = false;
                            for(const auto& ord : al.ordinals) {
                                if(ord[0] >= refMarkers.size()) continue;
                                const uint32_t rp = uint32_t(refMarkers[ord[0]].position);
                                if(rp < delZoneStart) hasLeft = true;
                                if(rp > delZoneEnd) hasRight = true;
                            }
                            if(!hasLeft || !hasRight) continue;
                            ++nChainsSpanning;

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

                        // Fallback: if no spanning chains found, use
                        // per-read DEL classifications whose breakpoints
                        // fall near the BP pair zone. These come from
                        // the initial per-read diagonal analysis.
                        if(delShifts.empty()) {
                            const uint32_t zoneMargin = windowSize * 4;
                            const uint32_t zoneL = delZoneStart > zoneMargin
                                ? delZoneStart - zoneMargin : 0;
                            const uint32_t zoneR = delZoneEnd + zoneMargin;

                            for(const auto& rg : readGroups) {
                                if(rg.svType != SvType::Deletion) continue;
                                if(rg.svSize < 30) continue;
                                if(rg.breakpointRefPos >= zoneL
                                   && rg.breakpointRefPos <= zoneR) {
                                    delShifts.push_back(rg.svSize);
                                }
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
                                delCallRecords.push_back({
                                    breakpointPos,
                                    medianDel,
                                    uint32_t(delShifts.size()),
                                    "diagonal"});
                                // When the left BP is to the RIGHT
                                // of the right BP (L > R), the
                                // diagonal shift may be a tandem
                                // repeat insertion rather than a
                                // deletion. Emit an INS call too
                                // so the correct type can be scored.
                                if(lbp.refPos > bestRbp->refPos
                                   && medianDel >= 50) {
                                    cout << "    >>> INSERTION CALL"
                                         << " (reversed-BP): size="
                                         << medianDel << "bp"
                                         << ", breakpoint="
                                         << breakpointPos
                                         << ", leftEnds="
                                         << lbp.endpointCount
                                         << ", rightStarts="
                                         << bestRbp->endpointCount
                                         << ", supportingReads="
                                         << delShifts.size()
                                         << endl;
                                }
                            }
                        }

                        // Fallback: when diagonal-shift finds no
                        // supporting chains but the BP pair is
                        // strong, emit a DEL call using refGap as
                        // the size estimate. In tandem repeats,
                        // chains often can't span the deletion
                        // zone, but the BP pair endpoints still
                        // mark the deletion boundaries.
                        if(delShifts.size() < 2
                           && bestDist >= 40
                           && (lbp.endpointCount
                               + lbp.ovhReadCount) >= 3
                           && (lbp.foldEnrichment >= 2.5
                               || bestRbp->foldEnrichment >= 2.5)) {
                            // Refine size using soft-clip BPs.
                            // bestDist is quantized to 50bp;
                            // soft-clip positions are base-precise.
                            int64_t bpPairSize = bestDist;
                            if(!softClipBPs.empty()) {
                                const int64_t sr =
                                    int64_t(windowSize);
                                int64_t dL = INT64_MAX;
                                uint32_t pL = lbp.refPos;
                                for(const auto& sc : softClipBPs) {
                                    if(sc.isLeftClip) continue;
                                    if(sc.readCount < 2) continue;
                                    const int64_t d = std::abs(
                                        int64_t(sc.refPos)
                                        - int64_t(lbp.refPos));
                                    if(d < dL) { dL = d; pL = sc.refPos; }
                                }
                                int64_t dR = INT64_MAX;
                                uint32_t pR = bestRbp->refPos;
                                for(const auto& sc : softClipBPs) {
                                    if(!sc.isLeftClip) continue;
                                    if(sc.readCount < 2) continue;
                                    const int64_t d = std::abs(
                                        int64_t(sc.refPos)
                                        - int64_t(bestRbp->refPos));
                                    if(d < dR) { dR = d; pR = sc.refPos; }
                                }
                                if(dL <= sr && dR <= sr
                                   && pR > pL) {
                                    const int64_t scD =
                                        int64_t(pR) - int64_t(pL);
                                    if(scD >= 40) {
                                        bpPairSize = scD;
                                    }
                                }
                            }
                            cout << "    >>> DELETION CALL"
                                 << " (bp-pair): size="
                                 << bpPairSize << "bp"
                                 << ", breakpoint="
                                 << breakpointPos
                                 << endl;
                            delCallRecords.push_back({
                                breakpointPos,
                                bpPairSize,
                                uint32_t(
                                    lbp.endpointCount
                                    + bestRbp->endpointCount),
                                "bp-pair"});


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
                                if(medianDel >= 50) {
                                    allDelCalls.push_back({
                                        breakpointPos,
                                        medianDel,
                                        uint32_t(delShifts.size()),
                                        "diagonal"});
                                }
                            }
                        }
                    }
                }

                // ---------------------------------------------------------
                // Single-breakpoint insertion detection.
                //
                // When only one strong breakpoint is found (left or right)
                // with many unanchored reads nearby, estimate insertion
                // size from the unanchored read count and coverage.
                // This handles cases where the insertion is larger than
                // the read length, so only one side of the breakpoint
                // is detectable.
                //
                // Suppress when a strong left-right BP pair exists with
                // a large refGap and no path — that pattern indicates a
                // deletion, not an insertion. The coverage-drop deletion
                // detector downstream will handle it.
                // ---------------------------------------------------------
                // Find the strongest left BP by ovhReadCount
                // (most supported). This is where the single-
                // breakpoint insertion would be called.
                const Breakpoint* strongestLBP = nullptr;
                uint32_t strongestLBPOvh = 0;
                for(const auto& lbp2 : leftBreakpoints) {
                    if(lbp2.ovhReadCount > strongestLBPOvh) {
                        strongestLBPOvh = lbp2.ovhReadCount;
                        strongestLBP = &lbp2;
                    }
                }
                // Only suppress single-BP insertion if a deletion-
                // like pair exists near the strongest BP. Pairs far
                // from the insertion region shouldn't suppress it.
                bool hasDeletionLikePair = false;
                for(const auto& lbp2 : leftBreakpoints) {
                    for(const auto& rbp2 : rightBreakpoints) {
                        const int64_t pairDist = std::abs(
                            int64_t(rbp2.refPos) - int64_t(lbp2.refPos));
                        if(pairDist >= 100
                           && pairDist <= 500
                           && lbp2.foldEnrichment >= 1.5
                           && rbp2.foldEnrichment >= 1.5
                           && lbp2.ovhReadCount >= 2
                           && rbp2.ovhReadCount >= 2) {
                            // Check proximity to strongest BP.
                            if(strongestLBP != nullptr) {
                                const int64_t distToStrong =
                                    std::min(
                                        std::abs(int64_t(lbp2.refPos)
                                            - int64_t(strongestLBP->refPos)),
                                        std::abs(int64_t(rbp2.refPos)
                                            - int64_t(strongestLBP->refPos)));
                                if(distToStrong > 500) continue;
                            }
                            hasDeletionLikePair = true;
                            break;
                        }
                    }
                    if(hasDeletionLikePair) break;
                }

                if(insertionCallRegions.empty()
                   && indirectAlignedReads.size() >= 10
                   && !hasDeletionLikePair) {
                    // Collect relaxed breakpoints: fold >= 1.5
                    // (weaker than the standard fold >= 3.0).
                    struct RelaxedBP {
                        uint32_t pos;
                        uint32_t count;
                        uint32_t ovhReads;
                        double fold;
                        bool isLeft;  // true = left BP (chain ends)
                    };
                    vector<RelaxedBP> relaxedBPs;

                    for(const auto& lbp : leftBreakpoints) {
                        if(lbp.foldEnrichment >= 1.5
                           && lbp.endpointCount >= 5) {
                            RelaxedBP rbpEntry;
                            rbpEntry.pos = lbp.refPos;
                            rbpEntry.count = lbp.endpointCount;
                            rbpEntry.ovhReads = lbp.ovhReadCount;
                            rbpEntry.fold = lbp.foldEnrichment;
                            rbpEntry.isLeft = true;
                            relaxedBPs.push_back(rbpEntry);
                        }
                    }
                    for(const auto& rbp : rightBreakpoints) {
                        if(rbp.foldEnrichment >= 1.5
                           && rbp.endpointCount >= 5) {
                            RelaxedBP rbpEntry;
                            rbpEntry.pos = rbp.refPos;
                            rbpEntry.count = rbp.endpointCount;
                            rbpEntry.ovhReads = rbp.ovhReadCount;
                            rbpEntry.fold = rbp.foldEnrichment;
                            rbpEntry.isLeft = false;
                            relaxedBPs.push_back(rbpEntry);
                        }
                    }

                    // Find the strongest BP (fold >= 3.0,
                    // ovhReads >= 3) and pair it with a relaxed
                    // BP on the opposite side within 500bp.
                    RelaxedBP* bestStrong = nullptr;
                    double bestStrongFold = 0;
                    for(auto& bp : relaxedBPs) {
                        if(bp.fold >= 3.0
                           && bp.ovhReads >= 3
                           && bp.fold > bestStrongFold) {
                            bestStrong = &bp;
                            bestStrongFold = bp.fold;
                        }
                    }

                    if(bestStrong != nullptr) {
                        // Estimate size from indirect reads.
                        const double coverage =
                            double(medianSpanning);
                        uint64_t indirectBases = 0;
                        for(const uint32_t rid :
                            indirectAlignedReads) {
                            indirectBases +=
                                readsRef.getRead(
                                    ReadId(rid)).baseCount;
                        }
                        const int64_t estSize =
                            (coverage > 0)
                            ? int64_t(double(indirectBases)
                                      / coverage)
                            : 0;

                        // Find best partner on opposite side.
                        RelaxedBP* bestPartner = nullptr;
                        double bestPartnerFold = 0;
                        for(auto& bp : relaxedBPs) {
                            if(bp.isLeft == bestStrong->isLeft)
                                continue;
                            const int64_t dist = std::abs(
                                int64_t(bp.pos)
                                - int64_t(bestStrong->pos));
                            if(dist <= 500
                               && bp.fold > bestPartnerFold) {
                                bestPartner = &bp;
                                bestPartnerFold = bp.fold;
                            }
                        }

                        if(bestPartner != nullptr
                           && estSize >= 50) {
                            const uint32_t bpPos =
                                (bestStrong->pos
                                 + bestPartner->pos) / 2;
                            cout << "    >>> INSERTION CALL"
                                 << " (large-ins): size="
                                 << estSize << "bp"
                                 << ", breakpoint=" << bpPos
                                 << ", leftEnds="
                                 << (bestStrong->isLeft
                                     ? bestStrong->count
                                     : bestPartner->count)
                                 << " (fold="
                                 << (bestStrong->isLeft
                                     ? bestStrong->fold
                                     : bestPartner->fold)
                                 << ")"
                                 << ", rightStarts="
                                 << (bestStrong->isLeft
                                     ? bestPartner->count
                                     : bestStrong->count)
                                 << " (fold="
                                 << (bestStrong->isLeft
                                     ? bestPartner->fold
                                     : bestStrong->fold)
                                 << ")"
                                 << ", internalReads="
                                 << indirectAlignedReads
                                    .size()
                                 << endl;
                            // Het-corrected estimate: for het
                            // insertions, indirectBases/coverage
                            // gives ~half the true size because
                            // only one allele contributes indirect
                            // reads while coverage counts both.
                            // Emit a 2x estimate when the ratio
                            // of internal reads to coverage
                            // suggests het (< 2.0).
                            const double irCovRatio =
                                double(indirectAlignedReads
                                       .size())
                                / double(medianSpanning);
                            if(irCovRatio < 2.0
                               && estSize >= 100) {
                                const int64_t hetSize =
                                    estSize * 2;
                                cout << "    >>> INSERTION CALL"
                                     << " (large-ins-het):"
                                     << " size="
                                     << hetSize << "bp"
                                     << ", breakpoint="
                                     << bpPos
                                     << ", irCovRatio="
                                     << irCovRatio
                                     << endl;
                            }
                            insertionCallRegions.push_back({
                                std::min(bestStrong->pos,
                                         bestPartner->pos),
                                std::max(bestStrong->pos,
                                         bestPartner->pos)
                            });
                            allInsCalls.push_back({
                                bpPos, estSize, 0,
                                "large-ins"});
                        }
                        // Single-BP case: only one strong BP,
                        // no partner within 500bp. Fire only
                        // when there's no strong BP on the
                        // opposite side at ANY distance (which
                        // would indicate a deletion pattern).
                        else if(bestPartner == nullptr
                                && bestStrong->fold >= 4.0
                                && bestStrong->ovhReads >= 10
                                && indirectAlignedReads.size()
                                   >= 20
                                && estSize >= 50
                                && estSize <= 1000) {
                            // Check: is there a strong BP on
                            // the opposite side at any distance?
                            bool hasOppositeBP = false;
                            for(const auto& bp : relaxedBPs) {
                                if(bp.isLeft
                                   != bestStrong->isLeft
                                   && bp.fold >= 3.0
                                   && bp.ovhReads >= 3) {
                                    hasOppositeBP = true;
                                    break;
                                }
                            }
                            if(!hasOppositeBP) {
                                cout << "    >>> INSERTION CALL"
                                     << " (large-ins): size="
                                     << estSize << "bp"
                                     << ", breakpoint="
                                     << bestStrong->pos
                                     << ", "
                                     << (bestStrong->isLeft
                                         ? "leftEnds="
                                         : "rightStarts=")
                                     << bestStrong->count
                                     << " (fold="
                                     << bestStrong->fold
                                     << ")"
                                     << ", internalReads="
                                     << indirectAlignedReads
                                        .size()
                                     << endl;
                                // Het-corrected estimate.
                                const double irCovRatio2 =
                                    double(indirectAlignedReads
                                           .size())
                                    / double(medianSpanning);
                                if(irCovRatio2 < 2.0
                                   && estSize >= 100) {
                                    const int64_t hetSize2 =
                                        estSize * 2;
                                    cout << "    >>> INSERTION"
                                         << " CALL"
                                         << " (large-ins-het):"
                                         << " size="
                                         << hetSize2 << "bp"
                                         << ", breakpoint="
                                         << bestStrong->pos
                                         << ", irCovRatio="
                                         << irCovRatio2
                                         << endl;
                                }
                                insertionCallRegions.push_back({
                                    bestStrong->pos > 200
                                    ? bestStrong->pos - 200
                                    : 0,
                                    bestStrong->pos + 200
                                });
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
                            // Emit a DEL call from the coverage
                            // drop span. In repetitive regions,
                            // chains can't span the deletion, but
                            // the coverage hole marks it.
                            const int64_t clusterSpan =
                                int64_t(cluster.endPos)
                                - int64_t(cluster.startPos)
                                + int64_t(windowSize);
                            const uint32_t hdBp =
                                (cluster.startPos + cluster.endPos) / 2;
                            if(clusterSpan >= 30) {
                                cout << "    >>> DELETION CALL"
                                     << " (covdrop-span): size="
                                     << clusterSpan << "bp"
                                     << ", breakpoint=" << hdBp
                                     << ", minRatio="
                                     << cluster.minRatio
                                     << endl;
                                allDelCalls.push_back({
                                    hdBp, clusterSpan,
                                    0, "covdrop-span"});
                            }
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
                                if(medianDel >= 50) {
                                    allDelCalls.push_back({
                                        bpPos,
                                        medianDel,
                                        uint32_t(
                                            cluster.readIds.size()),
                                        "diagonal"});
                                }
                            }
                        }
                    }
                }

                // ---------------------------------------------------------
                // Split-read deletion detection.
                //
                // For deletions larger than the chaining bandwidth,
                // no single chain can bridge both sides. Reads
                // spanning the deletion produce two chains: one on
                // each flank. The diagonal difference between these
                // chains estimates the deletion size.
                //
                // Group chains by read ID, then for each read with
                // 2+ chains, compare consecutive chain diagonals.
                // Cluster the resulting deletion sizes by ref position.
                // ---------------------------------------------------------
                {
                    struct SplitDelEvent {
                        uint32_t refPos;   // gap midpoint
                        int64_t delSize;   // diagonal difference
                        uint32_t readId;
                    };
                    vector<SplitDelEvent> splitEvents;

                    // Group chains by read ID.
                    std::unordered_map<uint32_t,
                        vector<uint32_t>> readChainIdx;
                    for(uint32_t ci = 0;
                        ci < chainsForRef.size(); ++ci) {
                        const auto& ce = chainsForRef[ci];
                        if(ce.readId == uint32_t(refId)) continue;
                        readChainIdx[ce.readId].push_back(ci);
                    }

                    for(const auto& [rdId, cis] : readChainIdx) {
                        if(cis.size() < 2) continue;

                        // For each chain, compute median diagonal
                        // and ref position range.
                        struct ChainSummary {
                            int64_t medDiag;
                            uint32_t minRefPos;
                            uint32_t maxRefPos;
                            uint32_t nAnchors;
                        };
                        vector<ChainSummary> summaries;
                        for(const auto ci : cis) {
                            const auto& ce = chainsForRef[ci];
                            if(!ce.isSameStrand) continue;
                            const auto& al =
                                alignments[ce.chainIndex];
                            if(al.ordinals.size() < 3) continue;
                            const Strand strand = 0;
                            const auto rdMkrs = markersRef[
                                OrientedReadId(
                                    ReadId(rdId), strand
                                ).getValue()];

                            vector<int64_t> diags;
                            uint32_t minRp = UINT32_MAX;
                            uint32_t maxRp = 0;
                            for(const auto& ord : al.ordinals) {
                                if(ord[0] >= refMarkers.size()
                                   || ord[1] >= rdMkrs.size())
                                    continue;
                                const uint32_t rp = uint32_t(
                                    refMarkers[ord[0]].position);
                                const uint32_t qp = uint32_t(
                                    rdMkrs[ord[1]].position);
                                diags.push_back(
                                    int64_t(rp) - int64_t(qp));
                                minRp = std::min(minRp, rp);
                                maxRp = std::max(maxRp, rp);
                            }
                            if(diags.size() < 3) continue;
                            sort(diags.begin(), diags.end());
                            summaries.push_back({
                                diags[diags.size() / 2],
                                minRp, maxRp,
                                uint32_t(diags.size())});
                        }

                        if(summaries.size() < 2) continue;

                        // Sort by ref position.
                        sort(summaries.begin(), summaries.end(),
                            [](const ChainSummary& a,
                               const ChainSummary& b) {
                                return a.minRefPos < b.minRefPos;
                            });

                        // Compare consecutive chain pairs.
                        for(size_t a = 0;
                            a + 1 < summaries.size(); ++a) {
                            const auto& ca = summaries[a];
                            const auto& cb = summaries[a + 1];
                            // Chains should be non-overlapping
                            // in ref space.
                            if(cb.minRefPos <= ca.maxRefPos)
                                continue;
                            // Reference gap between chains.
                            const int64_t refGap =
                                int64_t(cb.minRefPos)
                                - int64_t(ca.maxRefPos);
                            const int64_t dd =
                                cb.medDiag - ca.medDiag;
                            // The ref gap should be comparable
                            // to the diagonal difference (the
                            // deletion size). Reject if the gap
                            // is much larger — indicates chains
                            // mapping to distant regions.
                            if(dd > 50 && dd < 1000
                               && refGap < dd * 2) {
                                const uint32_t gapMid =
                                    (ca.maxRefPos + cb.minRefPos) / 2;
                                splitEvents.push_back(
                                    {gapMid, dd, rdId});
                            }
                        }
                    }

                    // Cluster split-read deletion events by position.
                    if(!splitEvents.empty()) {
                        sort(splitEvents.begin(), splitEvents.end(),
                            [](const SplitDelEvent& a,
                               const SplitDelEvent& b) {
                                return a.refPos < b.refPos;
                            });

                        struct SplitDelCluster {
                            uint32_t startPos;
                            uint32_t endPos;
                            vector<int64_t> sizes;
                            std::unordered_set<uint32_t> readIds;
                        };
                        vector<SplitDelCluster> splitClusters;
                        SplitDelCluster cur;
                        cur.startPos = splitEvents[0].refPos;
                        cur.endPos = splitEvents[0].refPos;
                        cur.sizes.push_back(splitEvents[0].delSize);
                        cur.readIds.insert(splitEvents[0].readId);

                        for(size_t i = 1;
                            i < splitEvents.size(); ++i) {
                            if(splitEvents[i].refPos
                               <= cur.endPos + windowSize * 3) {
                                cur.endPos =
                                    splitEvents[i].refPos;
                                cur.sizes.push_back(
                                    splitEvents[i].delSize);
                                cur.readIds.insert(
                                    splitEvents[i].readId);
                            } else {
                                splitClusters.push_back(
                                    std::move(cur));
                                cur = SplitDelCluster();
                                cur.startPos =
                                    splitEvents[i].refPos;
                                cur.endPos =
                                    splitEvents[i].refPos;
                                cur.sizes.push_back(
                                    splitEvents[i].delSize);
                                cur.readIds.insert(
                                    splitEvents[i].readId);
                            }
                        }
                        splitClusters.push_back(std::move(cur));

                        for(auto& cl : splitClusters) {
                            if(cl.readIds.size() < 2) continue;
                            sort(cl.sizes.begin(), cl.sizes.end());
                            const int64_t medDel =
                                cl.sizes[cl.sizes.size() / 2];
                            const uint32_t bpPos =
                                (cl.startPos + cl.endPos) / 2;

                            if(medDel > 50) {
                                cout << "    >>> DELETION CALL "
                                     << "(split-read): "
                                     << "size=" << medDel
                                     << "bp, breakpoint=" << bpPos
                                     << ", splitReads="
                                     << cl.readIds.size()
                                     << endl;
                                delCallRecords.push_back({
                                    bpPos,
                                    medDel,
                                    uint32_t(cl.readIds.size()),
                                    "split-read"});
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
                            const uint32_t irSize =
                                ir.endPos > ir.startPos
                                ? ir.endPos - ir.startPos : 0;
                            // Don't let small insertion calls
                            // suppress large coverage-drop regions.
                            // A small insertion at the edge of a
                            // deletion is likely a false positive
                            // from the BP-pair analysis.
                            if(irSize < delSize / 3) continue;
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
                        // A region is marker-depleted if either:
                        // (a) all windows have zero reference markers
                        //     (totalHitDepthWins == 0), or
                        // (b) >50% of windows with markers have low
                        //     hit depth.
                        const bool markerDepleted =
                            totalHitDepthWins == 0
                            || (double(lowHitDepthWins)
                                / double(totalHitDepthWins) > 0.5);

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

                        // Record for later k-mer cluster
                        // corroboration.
                        covDropRegions.push_back({
                            cdc.startPos, cdc.endPos,
                            markerDepleted});

                        // Suppress calls that overlap detected VNTR
                        // gaps or prior insertion calls, UNLESS the
                        // region is marker-depleted. In marker-
                        // depleted tandem repeats, the evidence is
                        // ambiguous between DEL and INS, so let the
                        // flank-gap analysis run and emit both.
                        if(overlapsVntr && !markerDepleted) continue;
                        // For insertion overlaps: still run the
                        // analysis. In tandem repeats, the same
                        // evidence can manifest as both INS and
                        // DEL. The diagonal-shift or flank-gap
                        // analysis may find the correct DEL size
                        // even when a large-ins call was emitted.

                        // Suppress large coverage-drop regions (>500bp)
                        // with minRatio=0 that have strong breakpoints
                        // at both edges AND very low spanning count
                        // inside the region. This pattern indicates a
                        // VNTR where chains don't span, not a real
                        // deletion. Real deletions have significant
                        // spanning chains in the flanking windows.
                        if(delSize > 500 && cdc.minRatio < 0.01) {
                            bool hasEdgeLeftBP = false;
                            bool hasEdgeRightBP = false;
                            uint32_t edgeLBPSpanning = 0;
                            uint32_t edgeRBPSpanning = 0;
                            for(const auto& lbp2 : leftBreakpoints) {
                                if(lbp2.refPos >= cdc.startPos - 200
                                   && lbp2.refPos <= cdc.startPos + 200
                                   && lbp2.foldEnrichment >= 2.0
                                   && lbp2.ovhReadCount >= 5) {
                                    hasEdgeLeftBP = true;
                                    edgeLBPSpanning = lbp2.spanCount;
                                    break;
                                }
                            }
                            for(const auto& rbp2 : rightBreakpoints) {
                                if(rbp2.refPos >= cdc.endPos - 200
                                   && rbp2.refPos <= cdc.endPos + 200
                                   && rbp2.foldEnrichment >= 2.0
                                   && rbp2.ovhReadCount >= 5) {
                                    hasEdgeRightBP = true;
                                    edgeRBPSpanning = rbp2.spanCount;
                                    break;
                                }
                            }
                            // Only suppress if spanning counts at
                            // both edges are low (<10). Real deletions
                            // have significant spanning at the edges
                            // because reads still chain across the
                            // flanking regions.
                            if(hasEdgeLeftBP && hasEdgeRightBP
                               && edgeLBPSpanning < 10
                               && edgeRBPSpanning < 10) {
                                cout << "    Suppressed: large VNTR-like "
                                     << "coverage-drop with edge BPs"
                                     << endl;
                                continue;
                            }
                        }

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
                        const uint32_t maxK = 62;
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
                                allDelCalls.push_back({
                                    bpPos, delSize,
                                    0, "coverage"});
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
                                allDelCalls.push_back({
                                    bestBp, bestSize,
                                    bestCount, "adaptive"});
                                refinedCall = true;
                            }
                        }

                        // For marker-depleted regions, try flank gap
                        // analysis first. In repeats, pairwise diffs
                        // produce artifact clusters at repeat-period
                        // multiples. The flank gap directly measures
                        // the bimodal split in per-read diagonals on
                        // each side of the gap, which is more robust.
                        if(!refinedCall && markerDepleted
                           && readResults.size() >= 6) {
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

                            auto analyzeFlankGapEarly = [](
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
                                analyzeFlankGapEarly(leftDiags);
                            const int64_t rightGap =
                                analyzeFlankGapEarly(rightDiags);

                            // Compute median diagonal for each
                            // flank to determine shift direction.
                            // DEL: right median > left median
                            //   (reads after deletion shift up)
                            // INS: right median < left median
                            //   (reads after insertion shift down)
                            int64_t leftMedian = 0, rightMedian = 0;
                            if(!leftDiags.empty()) {
                                sort(leftDiags.begin(),
                                     leftDiags.end());
                                leftMedian = leftDiags[
                                    leftDiags.size() / 2];
                            }
                            if(!rightDiags.empty()) {
                                sort(rightDiags.begin(),
                                     rightDiags.end());
                                rightMedian = rightDiags[
                                    rightDiags.size() / 2];
                            }
                            const int64_t medianShift =
                                rightMedian - leftMedian;

                            cout << "      Flank gaps (early): left="
                                 << leftGap << " ("
                                 << leftDiags.size()
                                 << " reads) right=" << rightGap
                                 << " (" << rightDiags.size()
                                 << " reads)"
                                 << " medianShift="
                                 << medianShift
                                 << endl;

                            int64_t flankShift = 0;
                            if(leftGap > 30 && rightGap > 30) {
                                flankShift =
                                    (leftGap + rightGap) / 2;
                            } else if(leftGap > 30) {
                                flankShift = leftGap;
                            } else if(rightGap > 30) {
                                flankShift = rightGap;
                            }

                            if(flankShift >= 40
                               && flankShift <= int64_t(delSize)) {
                                // Distinguish INS from DEL:
                                // In a deletion, chain-start BPs
                                // appear at the right edge of the
                                // coverage drop (reads from the
                                // non-deleted allele start there).
                                // In an insertion, no chain-start
                                // BPs appear because insertion-
                                // carrying reads simply don't chain.
                                bool hasRightStartBP = false;
                                for(const auto& rbp :
                                    rightBreakpoints) {
                                    if(rbp.refPos >= cdc.endPos - 100
                                       && rbp.refPos
                                          <= cdc.endPos + 200
                                       && rbp.endpointCount >= 5) {
                                        hasRightStartBP = true;
                                        break;
                                    }
                                }
                                const bool likelyInsertion =
                                    markerDepleted
                                    && !hasRightStartBP
                                    && indirectAlignedReads.size()
                                       >= 10
                                    && insertionCallRegions.empty();

                                if(likelyInsertion) {
                                    // In marker-depleted tandem
                                    // repeats, flankShift is one
                                    // repeat unit. The coverage-drop
                                    // size better approximates the
                                    // full insertion size.
                                    // Also emit a DEL call: for
                                    // tandem repeats, the evidence
                                    // is ambiguous between DEL and
                                    // INS. Emit both and let
                                    // downstream pick the correct
                                    // type.
                                    if(flankShift >= 40) {
                                        cout << "    >>> DELETION CALL"
                                             << " (flank-gap): size="
                                             << flankShift << "bp"
                                             << ", breakpoint="
                                             << bpPos
                                             << endl;
                                        delCallRecords.push_back({
                                            bpPos,
                                            flankShift,
                                            uint32_t(
                                                leftDiags.size()
                                                + rightDiags.size()),
                                            "flank-gap"});
                                    }
                                    const int64_t insCallSize =
                                        std::max(flankShift,
                                                 int64_t(delSize));
                                    cout << "    >>> INSERTION CALL"
                                         << " (flank-gap): size="
                                         << insCallSize << "bp"
                                         << ", breakpoint="
                                         << bpPos
                                         << ", indirectReads="
                                         << indirectAlignedReads
                                            .size()
                                         << endl;
                                    // Also emit a repeat-unit-
                                    // rounded estimate: the true
                                    // insertion is likely a whole
                                    // number of repeat units.
                                    // Use floor to avoid over-
                                    // estimating.
                                    if(flankShift >= 30
                                       && delSize > flankShift) {
                                        const int64_t nUnits =
                                            std::max(int64_t(1),
                                                int64_t(
                                                    double(delSize)
                                                    / double(
                                                        flankShift)));
                                        const int64_t roundedSize =
                                            flankShift * nUnits;
                                        if(roundedSize != insCallSize
                                           && roundedSize >= 50) {
                                            cout << "    >>> "
                                                 << "INSERTION CALL"
                                                 << " (flank-gap"
                                                 << "-rounded):"
                                                 << " size="
                                                 << roundedSize
                                                 << "bp"
                                                 << ", breakpoint="
                                                 << bpPos
                                                 << ", repeatUnit="
                                                 << flankShift
                                                 << ", nUnits="
                                                 << nUnits
                                                 << endl;
                                        }
                                    }
                                    insertionCallRegions.push_back(
                                        {cdc.startPos,
                                         cdc.endPos});
                                } else {
                                    // SA-tag refinement for
                                    // flank-gap DEL calls.
                                    // Use coverage-drop region
                                    // boundaries for proximity.
                                    const uint32_t fgSaMargin = 300;
                                    const uint32_t fgSaStart =
                                        cdc.startPos > fgSaMargin
                                        ? cdc.startPos - fgSaMargin
                                        : 0;
                                    const uint32_t fgSaEnd =
                                        cdc.endPos + fgSaMargin;
                                    for(const auto& sc :
                                        saTagCalls) {
                                        // In marker-depleted regions,
                                        // flank-gap sees one repeat
                                        // unit but SA-tag sees the
                                        // full deletion. Allow wider
                                        // size ratio with strong
                                        // SA-tag support.
                                        const double fgMaxR =
                                            (markerDepleted
                                             && sc.readCount >= 5)
                                            ? double(delSize)
                                              / double(
                                                  std::max(
                                                      flankShift,
                                                      int64_t(1)))
                                            : 2.0;
                                        if(sc.svType == "DEL"
                                           && sc.readCount >= 2
                                           && sc.size >= 30
                                           && sc.size <= 5000
                                           && sc.refPos >= fgSaStart
                                           && sc.refPos <= fgSaEnd
                                           && sc.size <= uint32_t(
                                                  flankShift * fgMaxR)
                                           && sc.size >= uint32_t(
                                                  flankShift * 0.3)){
                                            cout << "      SA-tag"
                                                 << " refine: "
                                                 << flankShift
                                                 << "bp -> "
                                                 << sc.size << "bp"
                                                 << " (SA reads="
                                                 << sc.readCount
                                                 << ")" << endl;
                                            flankShift = sc.size;
                                            break;
                                        }
                                    }
                                    cout << "    >>> DELETION CALL"
                                         << " (flank-gap): size="
                                         << flankShift << "bp"
                                         << ", breakpoint="
                                         << bpPos
                                         << endl;
                                    if(flankShift >= 50) {
                                        delCallRecords.push_back({
                                            bpPos,
                                            flankShift,
                                            0,
                                            "flank-gap"});
                                    }
                                    // In marker-depleted regions
                                    // with many indirect reads,
                                    // the "deletion" may be a
                                    // tandem repeat insertion.
                                    // Also emit an INS call using
                                    // the coverage-drop size.
                                    if(markerDepleted
                                       && indirectAlignedReads
                                              .size() >= 10) {
                                        const int64_t insSize =
                                            std::max(
                                                flankShift,
                                                int64_t(delSize));
                                        cout << "    >>> INSERTION"
                                             << " CALL (flank-gap"
                                             << "-alt): size="
                                             << insSize << "bp"
                                             << ", breakpoint="
                                             << bpPos
                                             << ", indirectReads="
                                             << indirectAlignedReads
                                                .size()
                                             << endl;
                                    }
                                }
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
                                // Strategy depends on whether the
                                // region is marker-depleted (tandem
                                // repeat). In tandem repeats, artifact
                                // clusters appear at non-deletion
                                // offsets; use weighted median to be
                                // robust. Otherwise, prefer the
                                // smallest cluster with support near
                                // the best.
                                uint32_t bestCount = 0;
                                if(markerDepleted) {
                                    // Tandem repeat region: use weighted
                                    // median of qualifying clusters.
                                    vector<DiffCluster> qualifying;
                                    for(const auto& cl : clusters) {
                                        if(cl.meanDiff >= 50
                                           && cl.meanDiff
                                              <= int64_t(delSize)) {
                                            qualifying.push_back(cl);
                                        }
                                    }
                                    if(qualifying.size() == 1) {
                                        bestShift = qualifying[0].meanDiff;
                                        bestCount = qualifying[0].count;
                                    } else if(qualifying.size() > 1) {
                                        sort(qualifying.begin(),
                                             qualifying.end(),
                                             [](const DiffCluster& a,
                                                const DiffCluster& b) {
                                                 return a.meanDiff
                                                        < b.meanDiff;
                                             });
                                        uint32_t totalCount = 0;
                                        for(const auto& cl : qualifying) {
                                            totalCount += cl.count;
                                        }
                                        const uint32_t medianIdx =
                                            totalCount / 2;
                                        uint32_t cumCount = 0;
                                        for(const auto& cl : qualifying) {
                                            cumCount += cl.count;
                                            if(cumCount >= medianIdx) {
                                                bestShift = cl.meanDiff;
                                                bestCount = cl.count;
                                                break;
                                            }
                                        }
                                    }
                                } else {
                                    // Non-repeat region: pick highest-
                                    // support cluster >= 50bp.
                                    for(const auto& cl : clusters) {
                                        if(cl.meanDiff >= 50
                                           && cl.meanDiff
                                              <= int64_t(delSize)
                                           && cl.count > bestCount) {
                                            bestCount = cl.count;
                                            bestShift = cl.meanDiff;
                                        }
                                    }
                                    if(bestCount > 0) {
                                        // Update count for selected.
                                        for(const auto& cl : clusters) {
                                            if(cl.meanDiff == bestShift) {
                                                bestCount = cl.count;
                                                break;
                                            }
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
                                // Check if an SA-tag DEL call
                                // nearby can refine the size.
                                // SA-tag uses aligner coordinates
                                // which handle repeats better than
                                // diagonal analysis.
                                //
                                // Use coverage-drop region boundaries
                                // for proximity (with margin) rather
                                // than fixed distance from center —
                                // large coverage-drop regions can have
                                // SA-tag breakpoints far from center
                                // but still within the region.
                                const uint32_t saProxMargin = 300;
                                const uint32_t saProxStart =
                                    cdc.startPos > saProxMargin
                                    ? cdc.startPos - saProxMargin : 0;
                                const uint32_t saProxEnd =
                                    cdc.endPos + saProxMargin;
                                for(const auto& sc : saTagCalls) {
                                    // Allow wider size range when
                                    // SA-tag has strong support.
                                    // In repeat regions, the bimodal
                                    // analysis picks one repeat unit
                                    // but the SA-tag sees the full
                                    // deletion — allow up to the
                                    // coverage-drop size.
                                    const double maxRatio =
                                        sc.readCount >= 5
                                        ? double(delSize)
                                          / double(std::max(
                                                bestShift,
                                                int64_t(1)))
                                        : 1.5;
                                    if(sc.svType == "DEL"
                                       && sc.readCount >= 2
                                       && sc.size >= 30
                                       && sc.size <= 5000
                                       && sc.refPos >= saProxStart
                                       && sc.refPos <= saProxEnd
                                       && sc.size <= uint32_t(
                                              bestShift * maxRatio)
                                       && sc.size >= uint32_t(
                                              bestShift * 0.3)) {
                                        cout << "      SA-tag refine:"
                                             << " " << bestShift
                                             << "bp -> "
                                             << sc.size << "bp"
                                             << " (SA reads="
                                             << sc.readCount
                                             << ")" << endl;
                                        bestShift = sc.size;
                                        break;
                                    }
                                }
                                cout << "    >>> DELETION CALL "
                                     << "(adaptive-bimodal): "
                                     << "size=" << bestShift << "bp, "
                                     << "breakpoint=" << bpPos
                                     << endl;
                                if(bestShift >= 50) {
                                    allDelCalls.push_back({
                                        bpPos, bestShift,
                                        0, "adaptive-bimodal"});
                                }
                                refinedCall = true;
                            }
                        }

                        // Marker-depleted insertion detection.
                        //
                        // When the adaptive analysis found no deletion
                        // signal in a marker-depleted region AND there
                        // are many indirect/unanchored reads, the
                        // coverage drop is likely from an insertion:
                        // reads carrying the inserted sequence can't
                        // chain to the reference.
                        if(!refinedCall && markerDepleted
                           && indirectAlignedReads.size() >= 10
                           && insertionCallRegions.empty()) {
                            uint64_t indirectBases = 0;
                            for(const uint32_t rid :
                                indirectAlignedReads) {
                                indirectBases +=
                                    readsRef.getRead(
                                        ReadId(rid)).baseCount;
                            }
                            const int64_t estInsSize =
                                (medianSpanning > 0)
                                ? int64_t(double(indirectBases)
                                          / double(medianSpanning))
                                : 0;

                            if(estInsSize >= 50
                               && estInsSize <= 2000
                               && indirectAlignedReads.size()
                                  >= uint32_t(medianSpanning) / 3) {
                                cout << "    >>> INSERTION CALL"
                                     << " (covdrop-indirect):"
                                     << " size=" << estInsSize
                                     << "bp, breakpoint="
                                     << bpPos
                                     << ", indirectReads="
                                     << indirectAlignedReads
                                        .size()
                                     << ", markerDepleted=1"
                                     << endl;
                                insertionCallRegions.push_back({
                                    cdc.startPos, cdc.endPos});
                                refinedCall = true;
                            }
                        }

                        if(!refinedCall && delSize >= 50
                           && delSize <= 2000) {
                            cout << "    >>> DELETION CALL (coverage): "
                                 << "size=" << delSize << "bp, "
                                 << "breakpoint=" << bpPos
                                 << endl;
                            allDelCalls.push_back({
                                bpPos, delSize,
                                0, "coverage"});
                        }
                    }
                }
            }
        }

        // -----------------------------------------------------------------
        // Post-coverage-drop: DEL → INS type-flip for tandem repeats.
        //
        // In tandem repeat insertions, diagonal-shift and split-read
        // analyses detect a "deletion" because the inserted sequence
        // is a copy of existing repeat units. When a DEL call has no
        // corresponding coverage-drop of similar size, it may be a
        // tandem repeat insertion. Emit an INS call of the same size.
        // -----------------------------------------------------------------
        for(const auto& dc : delCallRecords) {
            if(dc.size < 50 || dc.readCount < 2) continue;
            // Check if any coverage-drop region overlaps this DEL
            // call and has a similar size.
            bool hasCovDropSupport = false;
            for(const auto& cdr : covDropRegions) {
                const uint32_t cdrSize = cdr.endPos - cdr.startPos;
                // Coverage-drop overlaps the DEL breakpoint?
                if(dc.breakpointPos >= cdr.startPos
                   && dc.breakpointPos <= cdr.endPos) {
                    // Size within 3x?
                    if(cdrSize >= uint32_t(dc.size) / 3
                       && cdrSize <= uint32_t(dc.size) * 3) {
                        hasCovDropSupport = true;
                        break;
                    }
                }
            }
            if(!hasCovDropSupport) {
                cout << "    >>> INSERTION CALL"
                     << " (no-covdrop-flip): size="
                     << dc.size << "bp"
                     << ", breakpoint=" << dc.breakpointPos
                     << ", " << dc.source
                     << ", reads=" << dc.readCount
                     << endl;
            }
        }

        // -----------------------------------------------------------------
        // CIGAR INS corroboration with coverage-drop.
        //
        // In tandem repeat insertions, CIGAR sees one repeat unit
        // but the coverage-drop region approximates the full
        // insertion size. When a CIGAR INS cluster is near a
        // coverage-drop region, emit an INS call using the
        // coverage-drop size.
        // -----------------------------------------------------------------
        for(const auto& ci : cigarIndels) {
            if(ci.svType != "INS") continue;
            if(ci.readCount < 3 || ci.size < 30) continue;
            for(const auto& cdr : covDropRegions) {
                const uint32_t cdrSize = cdr.endPos - cdr.startPos;
                // CIGAR INS breakpoint within or near the
                // coverage-drop region?
                if(ci.refPos + 200 >= cdr.startPos
                   && ci.refPos <= cdr.endPos + 200
                   && cdrSize > ci.size
                   && cdrSize <= ci.size * 6
                   && cdrSize >= 100
                   && cdrSize <= 2000) {
                    cout << "    >>> INSERTION CALL"
                         << " (CIGAR-covdrop): size="
                         << cdrSize << "bp"
                         << ", breakpoint=" << ci.refPos
                         << ", cigarReads=" << ci.readCount
                         << ", cigarSize=" << ci.size
                         << endl;
                    break;
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
                    if(typeStr == "DEL"
                       && (sumSize / int64_t(count)) >= 20) {
                        allDelCalls.push_back({
                            uint32_t(sumPos / count),
                            sumSize / int64_t(count),
                            count, "cluster"});
                    }
                    if(typeStr == "INS"
                       && (sumSize / int64_t(count)) >= 20) {
                        allInsCalls.push_back({
                            uint32_t(sumPos / count),
                            sumSize / int64_t(count),
                            count, "INS-cluster"});
                    }
                    if(typeStr == "INV"
                       && (sumSize / int64_t(count)) >= 20) {
                        allDelCalls.push_back({
                            uint32_t(sumPos / count),
                            sumSize / int64_t(count),
                            count, "INV-cluster"});
                    }
                }
            }
        }

        // -----------------------------------------------------------------
        // Merge nearby per-read INS/DEL clusters with similar sizes.
        //
        // In tandem repeat regions, the same SV appears in multiple
        // reads at slightly different positions (because the repeat
        // unit can be inserted at any copy boundary). Merge clusters
        // of the same type within 500bp with sizes within 20% into
        // a single call.
        // -----------------------------------------------------------------
        {
            struct ClusterInfo {
                int64_t size;
                uint64_t pos;
                uint32_t reads;
                SvType type;
            };
            vector<ClusterInfo> allClusters;

            for(int32_t cid = 0; cid < totalClusters; ++cid) {
                uint32_t count = 0;
                int64_t sumSize = 0;
                uint64_t sumPos = 0;
                SvType cType = SvType::ReferenceLike;

                for(const auto& rg : readGroups) {
                    if(rg.clusterId != cid) continue;
                    ++count;
                    sumSize += rg.svSize;
                    sumPos += rg.breakpointRefPos;
                    cType = rg.svType;
                }
                if(count == 0) continue;
                allClusters.push_back({
                    sumSize / int64_t(count),
                    sumPos / count,
                    count,
                    cType
                });
            }

            // Sort by type then position.
            sort(allClusters.begin(), allClusters.end(),
                [](const ClusterInfo& a, const ClusterInfo& b) {
                    if(a.type != b.type) return int(a.type) < int(b.type);
                    return a.pos < b.pos;
                });

            // Merge nearby clusters of the same type with similar sizes.
            for(size_t i = 0; i < allClusters.size(); ) {
                const auto& base = allClusters[i];
                if(base.type != SvType::Insertion
                   && base.type != SvType::Deletion) {
                    ++i;
                    continue;
                }

                // Collect mergeable clusters.
                uint32_t totalReads = base.reads;
                int64_t weightedSize = base.size * int64_t(base.reads);
                uint64_t weightedPos = base.pos * uint64_t(base.reads);
                size_t j = i + 1;
                while(j < allClusters.size()
                      && allClusters[j].type == base.type
                      && allClusters[j].pos <= base.pos + 500
                      && std::abs(allClusters[j].size - base.size)
                         <= base.size / 5) {
                    totalReads += allClusters[j].reads;
                    weightedSize += allClusters[j].size
                                    * int64_t(allClusters[j].reads);
                    weightedPos += allClusters[j].pos
                                   * uint64_t(allClusters[j].reads);
                    ++j;
                }

                // Emit a call if merged from >= 2 clusters with
                // >= 3 reads, or if a single cluster has >= 15 reads
                // (strong standalone evidence, common for large SVs
                // where all reads see the same breakpoint).
                if((totalReads >= 3 && (j - i) >= 2)
                   || totalReads >= 15) {
                    const int64_t mergedSize =
                        weightedSize / int64_t(totalReads);
                    const uint64_t mergedPos =
                        weightedPos / uint64_t(totalReads);
                    const char* typeStr =
                        (base.type == SvType::Insertion) ? "INS" : "DEL";

                    if(mergedSize >= 50) {
                        cout << "    >>> "
                             << (base.type == SvType::Insertion
                                 ? "INSERTION" : "DELETION")
                             << " CALL (merged-clusters): "
                             << "size=" << mergedSize << "bp, "
                             << "breakpoint=" << mergedPos << ", "
                             << "clusters=" << (j - i) << ", "
                             << "reads=" << totalReads
                             << endl;
                        if(base.type != SvType::Insertion) {
                            allDelCalls.push_back({
                                uint32_t(mergedPos), mergedSize,
                                totalReads, "merged-clusters"});
                        }
                    }
                }
                i = j;
            }
        }

        // -----------------------------------------------------------------
        // CIGAR-corroborated k-mer cluster emission.
        //
        // In STR regions, k-mer chains may classify only 1-2 reads
        // as DEL (too few for merged-cluster emission), while CIGAR
        // net-deletion evidence independently confirms a deletion
        // at a nearby position. When a k-mer DEL cluster (1-2 reads,
        // size >= 50bp) has a CIGAR DEL cluster within 200bp, emit
        // the k-mer cluster's size (which better reflects the true
        // deletion through repeat-unit counting).
        // -----------------------------------------------------------------
        if(!cigarIndels.empty()) {
            for(int32_t cid = 0; cid < totalClusters; ++cid) {
                uint32_t count = 0;
                int64_t sumSize = 0;
                uint64_t sumPos = 0;
                SvType cType = SvType::ReferenceLike;

                for(const auto& rg : readGroups) {
                    if(rg.clusterId != cid) continue;
                    ++count;
                    sumSize += rg.svSize;
                    sumPos += rg.breakpointRefPos;
                    cType = rg.svType;
                }
                if(count == 0 || count > 2) continue;
                if(cType != SvType::Deletion) continue;
                const int64_t clusterSize =
                    sumSize / int64_t(count);
                if(clusterSize < 50) continue;
                const uint64_t clusterPos = sumPos / count;

                // Check for a nearby CIGAR DEL cluster.
                for(const auto& ci : cigarIndels) {
                    if(ci.svType != "DEL") continue;
                    if(ci.readCount < 3) continue;
                    const int64_t posDiff =
                        int64_t(ci.refPos) - int64_t(clusterPos);
                    if(std::abs(posDiff) <= 200) {
                        cout << "    >>> DELETION CALL"
                             << " (CIGAR-corroborated): "
                             << "size=" << clusterSize << "bp"
                             << ", breakpoint=" << clusterPos
                             << ", kmerReads=" << count
                             << ", cigarReads=" << ci.readCount
                             << ", cigarSize=" << ci.size << "bp"
                             << endl;
                        if(clusterSize >= 50) {
                            allDelCalls.push_back({
                                uint32_t(clusterPos),
                                clusterSize, count,
                                "CIGAR-corroborated"});
                        }
                        break;
                    }
                }
            }
        }

        // -----------------------------------------------------------------
        // Coverage-drop corroborated k-mer cluster emission.
        //
        // In marker-depleted regions, k-mer chains may classify
        // only 1-2 reads as DEL. When a k-mer DEL cluster falls
        // within a detected coverage-drop region, the coverage
        // drop confirms a deletion exists and the k-mer cluster
        // provides the size estimate.
        // -----------------------------------------------------------------
        if(!covDropRegions.empty()) {
            for(int32_t cid = 0; cid < totalClusters; ++cid) {
                uint32_t count = 0;
                int64_t sumSize = 0;
                uint64_t sumPos = 0;
                SvType cType = SvType::ReferenceLike;

                for(const auto& rg : readGroups) {
                    if(rg.clusterId != cid) continue;
                    ++count;
                    sumSize += rg.svSize;
                    sumPos += rg.breakpointRefPos;
                    cType = rg.svType;
                }
                if(count == 0 || count > 2) continue;
                if(cType != SvType::Deletion) continue;
                const int64_t clusterSize =
                    sumSize / int64_t(count);
                if(clusterSize < 100) continue;
                const uint64_t clusterPos = sumPos / count;

                // Check if this cluster falls within a
                // coverage-drop region.
                for(const auto& cdr : covDropRegions) {
                    if(!cdr.markerDepleted) continue;
                    // Allow 200bp margin — k-mer breakpoints
                    // may be slightly outside the coverage-drop
                    // boundaries due to windowing.
                    const uint32_t margin = 200;
                    if(clusterPos + margin >= cdr.startPos
                       && clusterPos <= cdr.endPos + margin) {
                        const uint32_t covDropSize =
                            cdr.endPos - cdr.startPos;
                        // Accept if k-mer size is within
                        // 2x of coverage-drop size (the
                        // coverage-drop region is often
                        // wider than the actual deletion).
                        if(clusterSize <= int64_t(covDropSize)
                           && clusterSize * 2
                              >= int64_t(covDropSize)) {
                            cout << "    >>> DELETION CALL"
                                 << " (covdrop-corroborated):"
                                 << " size="
                                 << clusterSize << "bp"
                                 << ", breakpoint="
                                 << clusterPos
                                 << ", kmerReads=" << count
                                 << ", covDropSize="
                                 << covDropSize << "bp"
                                 << endl;
                            if(clusterSize >= 50) {
                                allDelCalls.push_back({
                                    uint32_t(clusterPos),
                                    clusterSize, count,
                                    "covdrop-corroborated"});
                            }
                            break;
                        }
                    }
                }
            }
        }

        // -----------------------------------------------------------------
        // SDUST-gated low-complexity SV detection.
        //
        // For low-complexity regions (detected by SDUST), standard
        // chaining fails because k-mers are non-unique. Two strategies:
        //
        // Strategy 1 (small intervals, < ~read length): Find reads
        // whose chains anchor on both flanks of the SDUST interval,
        // compare read-space gap to reference-space gap.
        //
        // Strategy 2 (large intervals, > read length): No single read
        // can span the interval. Instead, estimate the SV size from
        // the total read bases covering the region vs expected from
        // flanking coverage.
        // -----------------------------------------------------------------
        {
            const vector<Base> dustRefSeq =
                readsRef.getOrientedReadRawSequence(
                    OrientedReadId(refId, 0));
            const uint32_t refSeqLen =
                uint32_t(dustRefSeq.size());

            // Run SDUST on the full reference.
            vector<pair<uint32_t,uint32_t>> dustIntervals;
            sdust(dustRefSeq, 0, refSeqLen, 20, 64,
                  dustIntervals);

            // Group chains by read ID (shared across intervals).
            std::unordered_map<uint32_t,
                vector<uint32_t>> dustReadChainMap;
            for(uint32_t ci = 0;
                ci < chainsForRef.size(); ++ci) {
                const auto& ce = chainsForRef[ci];
                if(ce.readId == uint32_t(refId)) continue;
                dustReadChainMap[ce.readId].push_back(ci);
            }

            // Compute average read length and total short reads.
            const uint32_t totalShortReads =
                uint32_t(readsRef.readCount()) - 1;
            uint64_t totalAllBases = 0;
            for(uint32_t ri = 1;
                ri < uint32_t(readsRef.readCount()); ++ri) {
                totalAllBases +=
                    readsRef.getRead(ReadId(ri)).baseCount;
            }
            const double avgReadLen =
                (totalShortReads > 0)
                ? double(totalAllBases)
                  / double(totalShortReads)
                : 150.0;

            for(const auto& [dustStart, dustEnd] :
                dustIntervals) {
                const uint32_t dustLen = dustEnd - dustStart;
                if(dustLen < 100) continue;

                // Detect repeat motif period (1-50bp).
                uint32_t motifPeriod = 0;
                {
                    const uint32_t checkLen =
                        std::min(dustLen, uint32_t(200));
                    for(uint32_t p = 1; p <= 100; ++p) {
                        uint32_t matches = 0;
                        uint32_t total = 0;
                        for(uint32_t i = dustStart;
                            i + p < dustStart + checkLen;
                            ++i) {
                            ++total;
                            if(dustRefSeq[i].value
                               == dustRefSeq[i + p].value)
                                ++matches;
                        }
                        if(total > 0
                           && matches * 100 / total >= 70) {
                            motifPeriod = p;
                            break;
                        }
                    }
                }

                // Flanking region.
                const uint32_t flankSize = 500;
                const uint32_t leftFlankStart =
                    dustStart > flankSize
                    ? dustStart - flankSize : 0;
                const uint32_t rightFlankEnd =
                    std::min(dustEnd + flankSize, refSeqLen);

                // Collect per-read anchor data near boundaries.
                struct FlankAnchor {
                    uint32_t readId;
                    uint32_t refPos;
                    uint32_t readPos;
                };
                vector<FlankAnchor> leftAnchors;
                vector<FlankAnchor> rightAnchors;

                struct SpanningRead {
                    uint32_t readId;
                    uint32_t leftRefPos, rightRefPos;
                    uint32_t leftReadPos, rightReadPos;
                };
                vector<SpanningRead> spanningReads;

                for(const auto& [rdId, cis] :
                    dustReadChainMap) {
                    uint32_t bestLeftRef = 0;
                    uint32_t bestLeftRead = 0;
                    uint32_t bestRightRef = UINT32_MAX;
                    uint32_t bestRightRead = 0;
                    bool hasLeft = false, hasRight = false;

                    for(const auto ci : cis) {
                        const auto& ce = chainsForRef[ci];
                        if(!ce.isSameStrand) continue;
                        const auto& al =
                            alignments[ce.chainIndex];
                        const Strand strand = 0;
                        const auto rdMkrs = markersRef[
                            OrientedReadId(
                                ReadId(rdId), strand
                            ).getValue()];

                        for(const auto& ord :
                            al.ordinals) {
                            if(ord[0] >= refMarkers.size()
                               || ord[1] >= rdMkrs.size())
                                continue;
                            const uint32_t rp = uint32_t(
                                refMarkers[ord[0]].position);
                            const uint32_t qp = uint32_t(
                                rdMkrs[ord[1]].position);

                            if(rp >= leftFlankStart
                               && rp < dustStart
                               && rp > bestLeftRef) {
                                bestLeftRef = rp;
                                bestLeftRead = qp;
                                hasLeft = true;
                            }
                            if(rp > dustEnd
                               && rp <= rightFlankEnd
                               && rp < bestRightRef) {
                                bestRightRef = rp;
                                bestRightRead = qp;
                                hasRight = true;
                            }
                        }
                    }

                    if(hasLeft) {
                        leftAnchors.push_back(
                            {rdId, bestLeftRef, bestLeftRead});
                    }
                    if(hasRight) {
                        rightAnchors.push_back(
                            {rdId, bestRightRef,
                             bestRightRead});
                    }
                    if(hasLeft && hasRight
                       && bestRightRead > bestLeftRead) {
                        spanningReads.push_back({
                            rdId,
                            bestLeftRef, bestRightRef,
                            bestLeftRead, bestRightRead});
                    }
                }

                // ---- Strategy 1: Spanning reads ----
                if(spanningReads.size() >= 2) {
                    vector<int64_t> sizeDiffs;
                    for(const auto& sr : spanningReads) {
                        const int64_t refGap =
                            int64_t(sr.rightRefPos)
                            - int64_t(sr.leftRefPos);
                        const int64_t readGap =
                            int64_t(sr.rightReadPos)
                            - int64_t(sr.leftReadPos);
                        sizeDiffs.push_back(readGap - refGap);
                    }
                    sort(sizeDiffs.begin(), sizeDiffs.end());
                    const int64_t medianDiff =
                        sizeDiffs[sizeDiffs.size() / 2];

                    cout << "    SDUST-STR region: "
                         << dustStart << "-" << dustEnd
                         << " (" << dustLen << "bp)"
                         << " motifPeriod=" << motifPeriod
                         << " spanningReads="
                         << spanningReads.size()
                         << " medianDiff=" << medianDiff
                         << endl;

                    if(std::abs(medianDiff) >= 10) {
                        const char* typeStr =
                            (medianDiff > 0) ? "INSERTION"
                                             : "DELETION";
                        cout << "    >>> " << typeStr
                             << " CALL (SDUST-STR): size="
                             << std::abs(medianDiff) << "bp"
                             << ", breakpoint="
                             << (dustStart + dustEnd) / 2
                             << ", spanningReads="
                             << spanningReads.size()
                             << ", motifPeriod="
                             << motifPeriod
                             << endl;
                        if(medianDiff < 0
                           && std::abs(medianDiff) >= 50) {
                            allDelCalls.push_back({
                                (dustStart + dustEnd) / 2,
                                int64_t(std::abs(medianDiff)),
                                uint32_t(spanningReads.size()),
                                "SDUST-STR"});
                        }
                    }
                    continue;
                }

                // ---- Strategy 2: Read-count depth (large) ----
                // For large SDUST intervals where no read spans,
                // count ALL reads whose primary chain overlaps
                // the SDUST interval. Their total bases divided
                // by flanking coverage estimates the sample VNTR
                // length. The difference from the reference VNTR
                // length is the SV size.
                //
                // Require: motifPeriod >= 4 (real tandem repeat,
                // not just simple low-complexity), interval >=
                // 500bp (base-count ratio needs large region to
                // be reliable), and enough flanking anchors.
                if(motifPeriod < 4) continue;
                if(dustLen < 1000) continue;
                // Require anchors on BOTH flanks for reliable
                // coverage estimation.
                if(leftAnchors.size() < 5
                   || rightAnchors.size() < 5)
                    continue;

                // Estimate flanking coverage.
                const uint32_t leftFlankLen =
                    dustStart - leftFlankStart;
                const uint32_t rightFlankLen =
                    rightFlankEnd - dustEnd;
                // Weighted average: total flank bases / total
                // flank length. More robust than averaging per-
                // side coverages when flanks have different sizes.
                const uint32_t totalFlankAnchors =
                    uint32_t(leftAnchors.size()
                             + rightAnchors.size());
                const uint32_t totalFlankBp =
                    leftFlankLen + rightFlankLen;
                const double flankCov =
                    (totalFlankBp > 0)
                    ? double(totalFlankAnchors) * avgReadLen
                      / double(totalFlankBp)
                    : 0.0;

                if(flankCov < 3.0) continue;

                // Count reads whose chain ref-position range
                // overlaps the SDUST interval.
                uint64_t vntrBases = 0;
                uint32_t vntrReadCount = 0;
                unordered_set<uint32_t> countedReads;

                for(const auto& [rdId, cis] :
                    dustReadChainMap) {
                    // Find the ref-position range of this read's
                    // chains.
                    uint32_t minRefPos = UINT32_MAX;
                    uint32_t maxRefPos = 0;
                    for(const auto ci : cis) {
                        const auto& ce = chainsForRef[ci];
                        if(!ce.isSameStrand) continue;
                        const auto& al =
                            alignments[ce.chainIndex];
                        for(const auto& ord :
                            al.ordinals) {
                            if(ord[0] >= refMarkers.size())
                                continue;
                            const uint32_t rp = uint32_t(
                                refMarkers[ord[0]].position);
                            minRefPos =
                                std::min(minRefPos, rp);
                            maxRefPos =
                                std::max(maxRefPos, rp);
                        }
                    }

                    // Check overlap with SDUST interval.
                    if(minRefPos < dustEnd
                       && maxRefPos > dustStart) {
                        const uint32_t bc = uint32_t(
                            readsRef.getRead(ReadId(rdId))
                            .baseCount);
                        vntrBases += bc;
                        ++vntrReadCount;
                        countedReads.insert(rdId);
                    }
                }

                // Add indirect reads (VNTR-internal, no ref
                // chain).
                for(const uint32_t rid :
                    indirectAlignedReads) {
                    if(countedReads.count(rid)) continue;
                    const uint32_t bc = uint32_t(
                        readsRef.getRead(ReadId(rid))
                        .baseCount);
                    vntrBases += bc;
                    ++vntrReadCount;
                    countedReads.insert(rid);
                }

                // Add truly unanchored reads (no chain to
                // anything) proportionally.
                uint32_t anchoredReadCount = 0;
                for(const auto& [rdId2, cis2] :
                    dustReadChainMap) {
                    (void)cis2;
                    ++anchoredReadCount;
                }
                const uint32_t unanchoredTotal =
                    (totalShortReads > anchoredReadCount)
                    ? totalShortReads - anchoredReadCount
                    : 0;
                // Subtract indirect reads already counted.
                uint32_t trueUnanchored = 0;
                {
                    uint32_t indirectNotInMap = 0;
                    for(const uint32_t rid :
                        indirectAlignedReads) {
                        if(dustReadChainMap.find(rid)
                           == dustReadChainMap.end())
                            ++indirectNotInMap;
                    }
                    trueUnanchored =
                        (unanchoredTotal > indirectNotInMap)
                        ? unanchoredTotal - indirectNotInMap
                        : 0;
                }
                const double dustFraction =
                    double(dustLen) / double(refSeqLen);
                vntrBases += uint64_t(
                    double(trueUnanchored)
                    * avgReadLen * dustFraction);
                vntrReadCount += uint32_t(
                    double(trueUnanchored) * dustFraction);

                // Sample VNTR length = vntrBases / coverage.
                const double sampleVntrLen =
                    vntrBases / flankCov;
                const int64_t estimatedSvSize =
                    int64_t(sampleVntrLen)
                    - int64_t(dustLen);

                cout << "    SDUST-VNTR region: "
                     << dustStart << "-" << dustEnd
                     << " (" << dustLen << "bp)"
                     << " motifPeriod=" << motifPeriod
                     << " leftAnchors="
                     << leftAnchors.size()
                     << " rightAnchors="
                     << rightAnchors.size()
                     << " vntrReads=" << vntrReadCount
                     << " vntrBases=" << vntrBases
                     << " flankCov=" << flankCov
                     << " sampleVntrLen="
                     << int64_t(sampleVntrLen)
                     << " estSize=" << estimatedSvSize
                     << endl;

                // Require the SV to be at least 5% of the
                // SDUST interval length — smaller changes are
                // within noise of the base-count approach.
                const double svFraction =
                    double(std::abs(estimatedSvSize))
                    / double(dustLen);
                if(std::abs(estimatedSvSize) >= 50
                   && svFraction >= 0.05) {
                    const char* typeStr =
                        (estimatedSvSize > 0) ? "INSERTION"
                                              : "DELETION";
                    cout << "    >>> " << typeStr
                         << " CALL (SDUST-VNTR): size="
                         << std::abs(estimatedSvSize) << "bp"
                         << ", breakpoint="
                         << (dustStart + dustEnd) / 2
                         << ", motifPeriod=" << motifPeriod
                         << ", flankCov=" << flankCov
                         << endl;
                    if(estimatedSvSize < 0) {
                        allDelCalls.push_back({
                            (dustStart + dustEnd) / 2,
                            int64_t(std::abs(estimatedSvSize)),
                            0, "SDUST-VNTR"});
                    }
                }
            }
        }

        // -----------------------------------------------------------------
        // SA tag evidence integration.
        //
        // When a BAM file is provided, SA tag split-read calls serve
        // as independent evidence. Emit SA-based calls that have
        // sufficient read support (>= 2 reads).
        //
        // Suppress SA-tag DEL calls when the region has a high
        // fraction of hit-depth BPs (indicating marker depletion)
        // and many supplementary alignments. In VNTRs, the aligner
        // maps split reads to different repeat copies, producing
        // false DEL calls.
        // -----------------------------------------------------------------
        if(!saTagCalls.empty()) {
            for(const auto& sc : saTagCalls) {
                if(sc.readCount >= 2 && sc.size >= 30
                   && sc.size <= 5000) {
                    // Suppress DEL calls in marker-depleted
                    // VNTR regions where the aligner maps
                    // split reads to different repeat copies.
                    // Exception: allow calls with very strong
                    // support (>=10 reads) — a real deletion
                    // in a VNTR region will have many
                    // consistent split reads.
                    if(sc.svType == "DEL"
                       && suppressSaTagDel
                       && sc.readCount < 10) {
                        cout << "    SA-tag DEL suppressed"
                             << " (marker-depleted VNTR):"
                             << " size=" << sc.size << "bp"
                             << ", breakpoint=" << sc.refPos
                             << ", reads=" << sc.readCount
                             << endl;
                        continue;
                    }
                    cout << "    >>> " << sc.svType
                         << " CALL (SA-tag): size="
                         << sc.size << "bp"
                         << ", breakpoint=" << sc.refPos
                         << ", reads=" << sc.readCount
                         << endl;
                    if(sc.svType == "DEL" && sc.size >= 50) {
                        allDelCalls.push_back({
                            sc.refPos, int64_t(sc.size),
                            sc.readCount, "SA-tag"});
                    }
                    // SA-tag DEL calls in tandem repeats may
                    // actually be insertions. Check against
                    // coverage-drop regions: if the SA-tag DEL
                    // is within a coverage-drop but much smaller
                    // than it, the "deletion" is likely an
                    // insertion (the aligner maps the supplementary
                    // to a different repeat copy).
                    if(sc.svType == "DEL"
                       && sc.readCount >= 3
                       && sc.size >= 50) {
                        for(const auto& cdr : covDropRegions) {
                            const uint32_t cdrSize =
                                cdr.endPos - cdr.startPos;
                            if(sc.refPos >= cdr.startPos
                               && sc.refPos <= cdr.endPos
                               && sc.size < cdrSize / 2) {
                                cout << "    >>> INSERTION CALL"
                                     << " (SA-DEL-flip): size="
                                     << sc.size << "bp"
                                     << ", breakpoint="
                                     << sc.refPos
                                     << ", reads="
                                     << sc.readCount
                                     << endl;
                                break;
                            }
                        }
                    }
                }
            }
        }

        // -----------------------------------------------------------------
        // Emit remaining direct DEL calls from delCallRecords.
        // -----------------------------------------------------------------
        {
            for(const auto& dc : delCallRecords) {
                if(dc.size >= 50) {
                    allDelCalls.push_back(dc);
                }
            }

            // Emit all direct DEL calls.
            for(const auto& dc : allDelCalls) {
                cout << "    >>> DELETION CALL"
                     << " (" << dc.source << "): size="
                     << dc.size << "bp"
                     << ", breakpoint="
                     << dc.breakpointPos << endl;
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


// Parse SA tags from a BAM file to extract split-read SV evidence.
vector<Assembler::SaTagSvCall> Assembler::parseSaTagSvCalls(
    const string& bamFileName,
    const string& refName,
    uint32_t refStart,
    uint32_t refEnd) const
{
    vector<SaTagSvCall> result;
    if(bamFileName.empty()) return result;

    htsFile* fp = hts_open(bamFileName.c_str(), "r");
    if(!fp) {
        cerr << "Warning: could not open BAM file: "
             << bamFileName << endl;
        return result;
    }

    sam_hdr_t* hdr = sam_hdr_read(fp);
    if(!hdr) {
        hts_close(fp);
        return result;
    }

    hts_idx_t* idx = sam_index_load(fp, bamFileName.c_str());
    if(!idx) {
        sam_hdr_destroy(hdr);
        hts_close(fp);
        cerr << "Warning: could not load BAM index for "
             << bamFileName << endl;
        return result;
    }

    // Try the reference name as-is, then with/without "chr" prefix.
    int tid = sam_hdr_name2tid(hdr, refName.c_str());
    if(tid < 0) {
        string altName;
        if(refName.size() > 3
           && refName.substr(0, 3) == "chr") {
            altName = refName.substr(3);
        } else {
            altName = "chr" + refName;
        }
        tid = sam_hdr_name2tid(hdr, altName.c_str());
    }
    if(tid < 0) {
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
        hts_close(fp);
        return result;
    }

    // Get the chromosome name as it appears in the BAM header.
    const string bamChrName = sam_hdr_tid2name(hdr, tid);

    hts_itr_t* iter = sam_itr_queryi(
        idx, tid, int(refStart), int(refEnd));
    if(!iter) {
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
        hts_close(fp);
        return result;
    }

    struct SaObs {
        bool isDel;
        int64_t size;
        uint32_t bpPos;
    };
    vector<SaObs> observations;

    bam1_t* aln = bam_init1();
    while(sam_itr_next(fp, iter, aln) >= 0) {
        if(aln->core.flag &
           (BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FUNMAP))
            continue;

        uint8_t* saTag = bam_aux_get(aln, "SA");
        if(!saTag) continue;
        const char* saStr = bam_aux2Z(saTag);
        if(!saStr) continue;

        const int32_t pPos = aln->core.pos;
        const uint32_t* pCigar = bam_get_cigar(aln);
        const int pNCigar = aln->core.n_cigar;

        int64_t pRefSpan = 0;
        for(int ci = 0; ci < pNCigar; ++ci) {
            const int op = bam_cigar_op(pCigar[ci]);
            const int len = bam_cigar_oplen(pCigar[ci]);
            if(op == BAM_CMATCH || op == BAM_CDEL
               || op == BAM_CREF_SKIP || op == BAM_CEQUAL
               || op == BAM_CDIFF) {
                pRefSpan += len;
            }
        }
        const int64_t pEnd = pPos + pRefSpan;

        // Parse SA tag entries.
        string saString(saStr);
        size_t sStart = 0;
        while(sStart < saString.size()) {
            size_t sEnd = saString.find(';', sStart);
            if(sEnd == string::npos) sEnd = saString.size();
            string entry = saString.substr(
                sStart, sEnd - sStart);
            sStart = sEnd + 1;
            if(entry.empty()) continue;

            vector<string> fields;
            size_t fStart = 0;
            while(fStart < entry.size()) {
                size_t fEnd = entry.find(',', fStart);
                if(fEnd == string::npos) fEnd = entry.size();
                fields.push_back(
                    entry.substr(fStart, fEnd - fStart));
                fStart = fEnd + 1;
            }
            if(fields.size() < 4) continue;

            if(fields[0] != bamChrName) continue;

            const int64_t saPos = stoll(fields[1]) - 1;
            const string& saCigar = fields[3];

            int64_t sRefSpan = 0;
            int64_t num = 0;
            for(char c : saCigar) {
                if(c >= '0' && c <= '9') {
                    num = num * 10 + (c - '0');
                } else {
                    if(c == 'M' || c == 'D' || c == 'N'
                       || c == '=' || c == 'X') {
                        sRefSpan += num;
                    }
                    num = 0;
                }
            }

            int64_t refGap;
            uint32_t bpPos;
            if(pPos <= saPos) {
                refGap = saPos - pEnd;
                bpPos = uint32_t(pEnd);
            } else {
                refGap = pPos - (saPos + sRefSpan);
                bpPos = uint32_t(saPos + sRefSpan);
            }

            if(refGap >= 30) {
                observations.push_back(
                    {true, refGap, bpPos});
            } else if(refGap <= -30) {
                observations.push_back(
                    {false, -refGap, bpPos});
            }
        }
    }

    bam_destroy1(aln);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    hts_close(fp);

    // Cluster observations by type and breakpoint (within 50bp).
    sort(observations.begin(), observations.end(),
         [](const SaObs& a, const SaObs& b) {
             if(a.isDel != b.isDel) return a.isDel > b.isDel;
             return a.bpPos < b.bpPos;
         });

    size_t i = 0;
    while(i < observations.size()) {
        const bool isDel = observations[i].isDel;
        uint32_t clusterStart = observations[i].bpPos;
        int64_t sumSize = 0;
        uint32_t sumPos = 0;
        uint32_t count = 0;
        size_t j = i;
        while(j < observations.size()
              && observations[j].isDel == isDel
              && observations[j].bpPos <= clusterStart + 50) {
            sumSize += observations[j].size;
            sumPos += observations[j].bpPos;
            ++count;
            ++j;
        }
        if(count >= 2) {
            result.push_back({
                isDel ? "DEL" : "INS",
                sumSize / int64_t(count),
                sumPos / count,
                count
            });
        }
        i = j;
    }

    return result;
}


// Multi-k unique anchor deletion sizing.
//
// For a given reference region and set of reads, find k-mers
// that are unique in the reference (trying k=maxK down to minK),
// match them against each read (requiring uniqueness in the read
// too), and measure deletion size from diagonal drops in the
// resulting (refPos, readPos) anchor pairs.
//
// Returns clustered deletion calls with size, breakpoint, and
// supporting read count.
vector<Assembler::MultiKDelCall> Assembler::multiKAnchorSizing(
    ReadId refId,
    const vector<ReadId>& readIds,
    uint32_t regionStart,
    uint32_t regionEnd,
    uint32_t minK,
    uint32_t maxK) const
{
    vector<MultiKDelCall> result;
    if(readIds.empty()) return result;

    const auto& readsRef = getReads();
    const uint32_t refSeqLen = uint32_t(
        readsRef.getRead(refId).baseCount);

    // Clamp region to reference bounds.
    if(regionEnd > refSeqLen) regionEnd = refSeqLen;
    if(regionStart >= regionEnd) return result;

    // Get reference sequence.
    const vector<Base> refSeq =
        readsRef.getOrientedReadRawSequence(
            OrientedReadId(refId, 0));

    // Use multiple k values (even spacing from minK to maxK).
    vector<uint32_t> kValues;
    for(uint32_t k = minK; k <= maxK; k += 4) {
        kValues.push_back(k);
    }
    if(kValues.empty()) return result;

    // Phase 1: Build reference k-mer index for all k values.
    // Include ALL k-mers (not just unique).
    std::unordered_map<string, vector<uint32_t>> refKmerIdx;
    for(const uint32_t k : kValues) {
        for(uint32_t p = regionStart;
            p + k <= regionEnd; ++p) {
            string kmer;
            kmer.reserve(k);
            for(uint32_t j = 0; j < k; ++j) {
                kmer.push_back(refSeq[p + j].character());
            }
            refKmerIdx[kmer].push_back(p);
        }
    }

    // Phase 2: For each read, compute weighted diagonal histogram.
    // Non-unique k-mers contribute with weight inversely
    // proportional to their multiplicity.
    struct DelSignal {
        double gapSize;
        uint32_t readIdx;
    };
    vector<DelSignal> allDelSignals;

    const uint32_t maxMult = 30;

    for(size_t ri = 0; ri < readIds.size(); ++ri) {
        const ReadId rid = readIds[ri];
        const vector<Base> readSeq =
            readsRef.getOrientedReadRawSequence(
                OrientedReadId(rid, 0));
        const uint32_t readLen = uint32_t(readSeq.size());
        if(readLen < kValues.front()) continue;

        // Build read k-mer index for all k values.
        std::unordered_map<string, vector<uint32_t>>
            readKmerIdx;
        for(const uint32_t k : kValues) {
            if(readLen < k) continue;
            for(uint32_t p = 0; p + k <= readLen; ++p) {
                string kmer;
                kmer.reserve(k);
                for(uint32_t j = 0; j < k; ++j) {
                    kmer.push_back(readSeq[p + j].character());
                }
                readKmerIdx[kmer].push_back(p);
            }
        }

        // Compute weighted diagonal scores.
        std::unordered_map<int64_t, double> diagScore;

        for(const auto& [kmer, readPositions] : readKmerIdx) {
            auto it = refKmerIdx.find(kmer);
            if(it == refKmerIdx.end()) continue;
            const auto& refPositions = it->second;
            const uint32_t k = uint32_t(kmer.size());

            const uint32_t refMult =
                uint32_t(refPositions.size());
            const uint32_t readMult =
                uint32_t(readPositions.size());
            if(refMult > maxMult || readMult > maxMult)
                continue;

            // Weight: k^2 / (refMult * readMult).
            // Longer k-mers and lower multiplicity = stronger.
            const double weight =
                double(k) * double(k)
                / (double(refMult) * double(readMult));

            for(const uint32_t rp : refPositions) {
                for(const uint32_t rdp : readPositions) {
                    const int64_t d =
                        int64_t(rp) - int64_t(rdp);
                    diagScore[d] += weight;
                }
            }
        }

        if(diagScore.empty()) continue;

        // Cluster nearby diagonals (within 3bp) into peaks.
        vector<int64_t> sortedDiags;
        sortedDiags.reserve(diagScore.size());
        for(const auto& [d, s] : diagScore) {
            sortedDiags.push_back(d);
        }
        sort(sortedDiags.begin(), sortedDiags.end());

        struct Peak {
            double center;
            double score;
        };
        vector<Peak> peaks;

        for(size_t i = 0; i < sortedDiags.size(); ) {
            double scoreSum = 0;
            double weightedSum = 0;
            size_t j = i;
            while(j < sortedDiags.size()
                  && (j == i
                      || sortedDiags[j] - sortedDiags[j-1]
                         <= 3)) {
                const double s = diagScore[sortedDiags[j]];
                scoreSum += s;
                weightedSum += double(sortedDiags[j]) * s;
                ++j;
            }
            peaks.push_back({weightedSum / scoreSum,
                             scoreSum});
            i = j;
        }

        // Sort peaks by score descending.
        sort(peaks.begin(), peaks.end(),
            [](const Peak& a, const Peak& b) {
                return a.score > b.score;
            });

        if(peaks.size() < 2) continue;

        const double topDiag = peaks[0].center;
        const double topScore = peaks[0].score;

        // Emit pairwise gaps between all strong peaks.
        // In tandem repeats, the deletion may show up as the
        // gap between peaks 2 and 3 (not 1 and 2), because
        // peak 2 is a repeat-period artifact. Emitting all
        // pairwise gaps lets cross-read clustering find the
        // true deletion.
        const size_t maxPeaks = std::min(peaks.size(),
                                         size_t(6));
        for(size_t pi = 0; pi < maxPeaks; ++pi) {
            if(pi > 0
               && peaks[pi].score < topScore * 0.05)
                break;
            for(size_t pj = pi + 1; pj < maxPeaks;
                ++pj) {
                if(peaks[pj].score < topScore * 0.05)
                    break;
                const double gap =
                    peaks[pi].center - peaks[pj].center;
                if(gap >= 20.0 && gap <= 2000.0) {
                    allDelSignals.push_back(
                        {gap, uint32_t(ri)});
                }
            }
        }
    }

    if(allDelSignals.empty()) return result;

    // Phase 3: Cluster deletion signals by size.
    // Count unique reads per cluster (a read can contribute
    // multiple signals from different secondary peaks).
    sort(allDelSignals.begin(), allDelSignals.end(),
        [](const DelSignal& a, const DelSignal& b) {
            return a.gapSize < b.gapSize;
        });

    vector<bool> used(allDelSignals.size(), false);

    for(size_t i = 0; i < allDelSignals.size(); ++i) {
        if(used[i]) continue;

        vector<size_t> cluster = {i};
        for(size_t j = i + 1;
            j < allDelSignals.size(); ++j) {
            if(used[j]) continue;
            const double ratio =
                std::min(allDelSignals[i].gapSize,
                         allDelSignals[j].gapSize)
                / std::max(allDelSignals[i].gapSize,
                           allDelSignals[j].gapSize);
            if(ratio >= 0.7) {
                cluster.push_back(j);
                used[j] = true;
            }
        }
        used[i] = true;

        // Compute average size and unique read count.
        double sizeSum = 0;
        std::set<uint32_t> uniqueReads;
        for(const size_t c : cluster) {
            sizeSum += allDelSignals[c].gapSize;
            uniqueReads.insert(allDelSignals[c].readIdx);
        }
        const int64_t avgSize =
            int64_t(sizeSum / double(cluster.size()) + 0.5);
        const uint32_t readCount =
            uint32_t(uniqueReads.size());

        // Deduplicate (same size ±10%).
        bool duplicate = false;
        for(const auto& existing : result) {
            const double r =
                double(std::min(existing.size, avgSize))
                / double(std::max(existing.size, avgSize));
            if(r > 0.9) {
                duplicate = true;
                break;
            }
        }
        if(!duplicate) {
            const uint32_t bp =
                (regionStart + regionEnd) / 2;
            result.push_back({bp, avgSize, readCount});
        }
    }

    return result;
}


// Parse soft-clip breakpoints and CIGAR indels from BAM reads
// in a single pass. Soft clips indicate breakpoint positions;
// CIGAR I/D operations indicate small-medium SVs directly.
void Assembler::parseBamEvidence(
    const string& bamFileName,
    const string& refName,
    uint32_t refStart,
    uint32_t refEnd,
    vector<SoftClipBreakpoint>& softClipBPs,
    vector<CigarIndelCall>& cigarIndels) const
{
    softClipBPs.clear();
    cigarIndels.clear();
    if(bamFileName.empty()) return;

    htsFile* fp = hts_open(bamFileName.c_str(), "r");
    if(!fp) return;
    sam_hdr_t* hdr = sam_hdr_read(fp);
    if(!hdr) { hts_close(fp); return; }
    hts_idx_t* idx = sam_index_load(fp, bamFileName.c_str());
    if(!idx) { sam_hdr_destroy(hdr); hts_close(fp); return; }

    int tid = sam_hdr_name2tid(hdr, refName.c_str());
    if(tid < 0) {
        string altName;
        if(refName.size() > 3 && refName.substr(0, 3) == "chr")
            altName = refName.substr(3);
        else
            altName = "chr" + refName;
        tid = sam_hdr_name2tid(hdr, altName.c_str());
    }
    if(tid < 0) {
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
        hts_close(fp);
        return;
    }

    hts_itr_t* iter = sam_itr_queryi(
        idx, tid, int(refStart), int(refEnd));
    if(!iter) {
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
        hts_close(fp);
        return;
    }

    const uint32_t minClipLen = 20;
    const int64_t minIndelSize = 30;

    // Raw observations before clustering.
    struct ClipObs {
        uint32_t refPos;   // absolute position
        bool isLeftClip;
        string clipSeq;
    };
    vector<ClipObs> clipObs;

    struct IndelObs {
        bool isDel;
        int64_t size;
        uint32_t refPos;   // absolute position
        string insSeq;
    };
    vector<IndelObs> indelObs;

    bam1_t* aln = bam_init1();
    while(sam_itr_next(fp, iter, aln) >= 0) {
        if(aln->core.flag &
           (BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FUNMAP))
            continue;

        const int32_t pos = aln->core.pos;
        const uint32_t* cigar = bam_get_cigar(aln);
        const int nCigar = aln->core.n_cigar;
        const uint8_t* seq = bam_get_seq(aln);
        const int32_t seqLen = aln->core.l_qseq;

        if(nCigar == 0) continue;

        // Check left soft clip (first CIGAR op).
        {
            const int op = bam_cigar_op(cigar[0]);
            const int len = bam_cigar_oplen(cigar[0]);
            if(op == BAM_CSOFT_CLIP
               && uint32_t(len) >= minClipLen) {
                // Left clip: breakpoint is at alignment start.
                string clipSeq;
                clipSeq.reserve(len);
                for(int i = 0; i < len; ++i) {
                    static const char base[] = "=ACMGRSVTWYHKDBN";
                    clipSeq += base[bam_seqi(seq, i)];
                }
                clipObs.push_back({
                    uint32_t(pos), true, std::move(clipSeq)
                });
            }
        }

        // Check right soft clip (last CIGAR op).
        {
            const int op = bam_cigar_op(cigar[nCigar - 1]);
            const int len = bam_cigar_oplen(cigar[nCigar - 1]);
            if(op == BAM_CSOFT_CLIP
               && uint32_t(len) >= minClipLen) {
                // Right clip: breakpoint is at alignment end.
                // Compute alignment end position.
                int64_t refPos = pos;
                for(int ci = 0; ci < nCigar - 1; ++ci) {
                    const int o = bam_cigar_op(cigar[ci]);
                    const int l = bam_cigar_oplen(cigar[ci]);
                    if(o == BAM_CMATCH || o == BAM_CDEL
                       || o == BAM_CREF_SKIP || o == BAM_CEQUAL
                       || o == BAM_CDIFF)
                        refPos += l;
                }
                string clipSeq;
                const int clipStart = seqLen - len;
                clipSeq.reserve(len);
                for(int i = clipStart; i < seqLen; ++i) {
                    static const char base[] = "=ACMGRSVTWYHKDBN";
                    clipSeq += base[bam_seqi(seq, i)];
                }
                clipObs.push_back({
                    uint32_t(refPos), false, std::move(clipSeq)
                });
            }
        }

        // Check CIGAR for large I/D operations.
        // Collect per-read D operations first, then merge nearby
        // ones into compound deletions (handles tandem repeat
        // regions where the aligner splits one deletion into
        // multiple D ops separated by short matches).
        {
            struct DelOp {
                uint32_t refPos;
                int64_t  size;
            };
            vector<DelOp> readDels;

            int64_t refPos = pos;
            int32_t queryPos = 0;
            // Skip leading soft clip in query position.
            if(nCigar > 0
               && bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) {
                queryPos = bam_cigar_oplen(cigar[0]);
            }

            // Also track net indel effect for this read
            // (sum of all D minus sum of all I, regardless
            // of individual op size).
            int64_t totalDel = 0;
            int64_t totalIns = 0;
            // Track position of the largest D op (even if
            // below minIndelSize) for net-CIGAR breakpoint.
            uint32_t largestDelPos = uint32_t(pos);
            int64_t largestDelSize = 0;

            for(int ci = 0; ci < nCigar; ++ci) {
                const int op = bam_cigar_op(cigar[ci]);
                const int len = bam_cigar_oplen(cigar[ci]);

                if(op == BAM_CDEL && len >= minIndelSize) {
                    readDels.push_back({
                        uint32_t(refPos), int64_t(len)
                    });
                }
                if(op == BAM_CDEL) {
                    totalDel += len;
                    if(int64_t(len) > largestDelSize) {
                        largestDelSize = int64_t(len);
                        largestDelPos = uint32_t(refPos);
                    }
                }
                if(op == BAM_CINS) totalIns += len;

                if(op == BAM_CINS && len >= minIndelSize) {
                    string insSeq;
                    insSeq.reserve(len);
                    for(int i = queryPos;
                        i < queryPos + len && i < seqLen; ++i) {
                        static const char base[] =
                            "=ACMGRSVTWYHKDBN";
                        insSeq += base[bam_seqi(seq, i)];
                    }
                    indelObs.push_back({
                        false, int64_t(len),
                        uint32_t(refPos), std::move(insSeq)
                    });
                }

                // Advance positions.
                if(op == BAM_CMATCH || op == BAM_CEQUAL
                   || op == BAM_CDIFF) {
                    refPos += len;
                    queryPos += len;
                } else if(op == BAM_CDEL
                          || op == BAM_CREF_SKIP) {
                    refPos += len;
                } else if(op == BAM_CINS
                          || op == BAM_CSOFT_CLIP) {
                    queryPos += len;
                }
            }

            // Merge nearby D operations within this read.
            // If two D ops are within 100bp on the reference,
            // combine them into a single compound deletion.
            // This handles tandem repeats where e.g. 35D + 55D
            // should be reported as a single 90D.
            if(readDels.size() >= 2) {
                size_t wi = 0;
                for(size_t ri = 1; ri < readDels.size(); ++ri) {
                    const uint32_t gap =
                        readDels[ri].refPos
                        - readDels[wi].refPos
                        - uint32_t(readDels[wi].size);
                    if(gap <= 100) {
                        // Merge: sum deleted bases (not
                        // reference span) so the size
                        // matches truth-set conventions.
                        readDels[wi].size += readDels[ri].size;
                    } else {
                        ++wi;
                        readDels[wi] = readDels[ri];
                    }
                }
                readDels.resize(wi + 1);
            }

            // Emit individual (possibly merged) D observations.
            for(const auto& d : readDels) {
                indelObs.push_back({
                    true, d.size, d.refPos, ""
                });
            }

            // Net-CIGAR deletion: if the read's total CIGAR
            // deletions exceed insertions by >= 30bp, record
            // a net-deletion observation. This catches STR
            // regions where the aligner fragments the deletion
            // into many small D ops below minIndelSize.
            const int64_t netDel = totalDel - totalIns;
            if(netDel >= minIndelSize
               && readDels.empty()) {
                // Use position of the largest D operation as
                // breakpoint (not alignment start), so reads
                // seeing the same deletion cluster together.
                indelObs.push_back({
                    true, netDel, largestDelPos, ""
                });
            }
        }
    }

    bam_destroy1(aln);
    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    hts_close(fp);

    // Cluster soft-clip observations by position (within 5bp).
    sort(clipObs.begin(), clipObs.end(),
         [](const ClipObs& a, const ClipObs& b) {
             if(a.isLeftClip != b.isLeftClip)
                 return a.isLeftClip > b.isLeftClip;
             return a.refPos < b.refPos;
         });

    {
        size_t i = 0;
        while(i < clipObs.size()) {
            const bool isLeft = clipObs[i].isLeftClip;
            const uint32_t clusterStart = clipObs[i].refPos;
            uint64_t sumPos = 0;
            uint64_t sumClipLen = 0;
            vector<string> seqs;
            size_t j = i;
            while(j < clipObs.size()
                  && clipObs[j].isLeftClip == isLeft
                  && clipObs[j].refPos <= clusterStart + 5) {
                sumPos += clipObs[j].refPos;
                sumClipLen += clipObs[j].clipSeq.size();
                seqs.push_back(std::move(clipObs[j].clipSeq));
                ++j;
            }
            const uint32_t count = uint32_t(j - i);
            if(count >= 3) {
                uint32_t localPos = uint32_t(sumPos / count);
                if(localPos >= refStart)
                    localPos -= refStart;
                softClipBPs.push_back({
                    localPos,
                    count,
                    isLeft,
                    uint32_t(sumClipLen / count),
                    std::move(seqs)
                });
            }
            i = j;
        }
    }

    // Cluster CIGAR indel observations by type, position
    // (within 20bp), and size (within 50% of cluster median).
    // Size gating prevents averaging different-sized deletions
    // in tandem repeat regions (e.g. 35D and 90D at the same
    // locus should form separate clusters).
    sort(indelObs.begin(), indelObs.end(),
         [](const IndelObs& a, const IndelObs& b) {
             if(a.isDel != b.isDel) return a.isDel > b.isDel;
             if(a.refPos != b.refPos) return a.refPos < b.refPos;
             return a.size < b.size;
         });

    {
        // First pass: group by position (within 20bp).
        // Second pass: split each position group by size.
        size_t i = 0;
        while(i < indelObs.size()) {
            const bool isDel = indelObs[i].isDel;
            const uint32_t clusterStart = indelObs[i].refPos;

            // Find extent of position cluster.
            size_t posEnd = i;
            while(posEnd < indelObs.size()
                  && indelObs[posEnd].isDel == isDel
                  && indelObs[posEnd].refPos <= clusterStart + 20) {
                ++posEnd;
            }

            // Sort this position group by size.
            sort(indelObs.begin() + int64_t(i),
                 indelObs.begin() + int64_t(posEnd),
                 [](const IndelObs& a, const IndelObs& b) {
                     return a.size < b.size;
                 });

            // Sub-cluster by size: observations join the
            // cluster if within 50% of the running mean.
            size_t si = i;
            while(si < posEnd) {
                int64_t sumSize = indelObs[si].size;
                uint64_t sumPos = indelObs[si].refPos;
                string bestInsSeq = indelObs[si].insSeq;
                size_t sj = si + 1;
                while(sj < posEnd) {
                    const int64_t meanSize =
                        sumSize / int64_t(sj - si);
                    // Allow joining if within 50% of mean.
                    if(indelObs[sj].size
                       <= meanSize * 3 / 2) {
                        sumSize += indelObs[sj].size;
                        sumPos += indelObs[sj].refPos;
                        if(indelObs[sj].insSeq.size()
                           > bestInsSeq.size())
                            bestInsSeq = indelObs[sj].insSeq;
                        ++sj;
                    } else {
                        break;
                    }
                }
                const uint32_t count = uint32_t(sj - si);
                if(count >= 2) {
                    uint32_t localPos =
                        uint32_t(sumPos / count);
                    if(localPos >= refStart)
                        localPos -= refStart;
                    cigarIndels.push_back({
                        isDel ? "DEL" : "INS",
                        sumSize / int64_t(count),
                        localPos,
                        count,
                        std::move(bestInsSeq)
                    });
                }
                si = sj;
            }

            i = posEnd;
        }
    }
}


