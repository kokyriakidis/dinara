// Dinara.
#include "Assembler.hpp"
#include "Reads.hpp"
#include "chrono.hpp"
#include "hifiasmECInternals.hpp"
#include "timestamp.hpp"

// Standard library.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

struct GlobalReadSite {
    uint32_t siteId = 0;
    uint32_t position = 0;
    uint8_t allele = 0;
};

struct GlobalReadSiteIndex {
    // Global het sites per read, sorted by position in the read's forward coordinate space.
    // Used to iterate over sites within an alignment's [qs, qe) window.
    // The target-side concordance check no longer needs a second per-site index:
    // it is done by querying the pre-computed alignedEvidenceStore SNP streams directly.
    vector< vector<GlobalReadSite> > sitesByReadPos;
    uint64_t keptSiteCount = 0;
};

inline bool positionOrder(const GlobalReadSite& a, const GlobalReadSite& b)
{
    if (a.position != b.position) {
        return a.position < b.position;
    }
    return a.siteId < b.siteId;
}

inline bool siteOrder(const GlobalReadSite& a, const GlobalReadSite& b)
{
    if (a.siteId != b.siteId) {
        return a.siteId < b.siteId;
    }
    return a.position < b.position;
}

inline GlobalReadSiteIndex buildGlobalReadSiteIndexFromClusters(
    const Assembler& assembler,
    const Assembler::GlobalMismatchSiteClusters& clusters,
    uint32_t minAlleleSupport,
    uint32_t minAllelesWithSupport)
{
    GlobalReadSiteIndex out;

    const Reads& reads = assembler.getReads();
    const uint64_t readCount = reads.readCount();
    out.sitesByReadPos.assign(readCount, vector<GlobalReadSite>{});

    const uint32_t siteCount = clusters.clusterMemberOffsets.empty() ?
        0U :
        uint32_t(clusters.clusterMemberOffsets.size() - 1);
    if (siteCount == 0) {
        return out;
    }

    for (uint32_t siteId = 0; siteId < siteCount; siteId++) {
        const uint64_t begin = clusters.clusterMemberOffsets[siteId];
        const uint64_t end = clusters.clusterMemberOffsets[siteId + 1];

        unordered_map<uint32_t, array<uint32_t, 4> > countByRead;
        unordered_map<uint32_t, array<uint32_t, 4> > minPosByRead;
        countByRead.reserve(size_t(end - begin + 8));
        minPosByRead.reserve(size_t(end - begin + 8));

        for (uint64_t i = begin; i < end; i++) {
            const auto& node = clusters.nodes[clusters.clusterMembers[i]];
            const uint32_t rid = uint32_t(node.first);
            const uint32_t pos = node.second;
            const uint8_t base =
                reads.getOrientedReadBase(OrientedReadId(node.first, 0), pos).value;
            if (base >= 4) {
                continue;
            }

            auto countIt = countByRead.find(rid);
            if (countIt == countByRead.end()) {
                countIt = countByRead.insert({rid, array<uint32_t, 4>{0, 0, 0, 0}}).first;
                minPosByRead.insert({
                    rid,
                    array<uint32_t, 4>{
                        numeric_limits<uint32_t>::max(),
                        numeric_limits<uint32_t>::max(),
                        numeric_limits<uint32_t>::max(),
                        numeric_limits<uint32_t>::max()}
                });
            }
            countIt->second[base]++;
            auto& minPos = minPosByRead[rid];
            if (pos < minPos[base]) {
                minPos[base] = pos;
            }
        }

        if (countByRead.empty()) {
            continue;
        }

        struct Winner {
            uint32_t rid = 0;
            uint32_t pos = 0;
            uint8_t allele = 0;
        };
        vector<Winner> winners;
        winners.reserve(countByRead.size());
        array<uint32_t, 4> winnerAlleleCounts{0, 0, 0, 0};

        for (const auto& [rid, counts] : countByRead) {
            uint32_t bestAllele = 0;
            uint32_t bestCount = 0;
            bool tied = false;
            for (uint32_t allele = 0; allele < 4; allele++) {
                if (counts[allele] > bestCount) {
                    bestCount = counts[allele];
                    bestAllele = allele;
                    tied = false;
                } else if (counts[allele] == bestCount && counts[allele] > 0) {
                    tied = true;
                }
            }
            if (bestCount == 0 || tied) {
                continue;
            }

            const auto minPosIt = minPosByRead.find(rid);
            if (minPosIt == minPosByRead.end()) {
                continue;
            }
            const uint32_t bestPos = minPosIt->second[bestAllele];
            if (bestPos == numeric_limits<uint32_t>::max()) {
                continue;
            }

            winners.push_back(Winner{
                rid,
                bestPos,
                uint8_t(bestAllele)});
            winnerAlleleCounts[bestAllele]++;
        }

        if (winners.empty()) {
            continue;
        }

        uint32_t supportedAlleles = 0;
        for (uint32_t allele = 0; allele < 4; allele++) {
            if (winnerAlleleCounts[allele] >= minAlleleSupport) {
                supportedAlleles++;
            }
        }
        if (supportedAlleles < minAllelesWithSupport) {
            continue;
        }

        out.keptSiteCount++;
        for (const auto& w : winners) {
            if (w.rid >= out.sitesByReadPos.size()) {
                continue;
            }
            out.sitesByReadPos[w.rid].push_back(GlobalReadSite{
                siteId,
                w.pos,
                w.allele});
        }
    }

    for (uint64_t rid = 0; rid < readCount; rid++) {
        auto& byPos = out.sitesByReadPos[rid];
        sort(byPos.begin(), byPos.end(), positionOrder);
    }

    return out;
}

struct CandidateGlobal {
    uint32_t alignmentId = invalid<uint32_t>;
    uint32_t qs = 0;
    uint32_t qe = 0;
    uint32_t ts = 0;
    uint32_t te = 0;
    uint32_t targetId = 0;
    bool targetIsRc = false;
    uint8_t isMatch = 1;
};

} // namespace


void Assembler::performGlobalSiteECParity(uint64_t threadCount)
{
    cout << timestamp << "=== Global-Site EC Pipeline (Experimental) ===" << endl;

    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    if (!alignmentTable.isOpen()) {
        throw runtime_error("performGlobalSiteECParity requires alignmentTable to be open.");
    }

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    if (threadCount == 0) {
        threadCount = 1;
    }

    AlignOptions clusteringAlignOptions{};
    clusteringAlignOptions.maxErrorRate = 1.0;

    const auto tBeginAll = steady_clock::now();
    const auto tBuildSitesBegin = steady_clock::now();
    const GlobalMismatchSiteClusters clusters = clusterMismatchingPositionsIntoGlobalHetSites(
        clusteringAlignOptions,
        threadCount,
        false,
        false);
    // Global-site EC is primarily a phasing/consistency filter. For low-coverage or small
    // integration scenarios we still want to keep biallelic sites even if each allele is
    // supported by only a single read (the later decision step still requires multiple
    // shared sites, see minSharedSitesForDecision).
    const GlobalReadSiteIndex globalIndex = buildGlobalReadSiteIndexFromClusters(
        *this,
        clusters,
        1, // minAlleleSupport
        2  // minAllelesWithSupport
    );
    const double tBuildSites = seconds(steady_clock::now() - tBuildSitesBegin);

    const uint64_t readCount = reads->readCount();
    vector<thread> threads;
    const uint64_t chunkSize = max<uint64_t>(1, readCount / threadCount);

    struct alignas(64) Timing {
        double gatherCandidates = 0.;
        double phaseUsingGlobalSites = 0.;
        double finalizeFlags = 0.;
        uint64_t readsVisited = 0;
        uint64_t readsWithAlignments = 0;
        uint64_t readsWithCandidates = 0;
    };
    vector<Timing> timings(threadCount);

    static constexpr uint32_t minInformativeSitesPerRead = 2;
    static constexpr uint32_t minSharedSitesForDecision = 2;

    for (uint64_t t = 0; t < threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t start = t * chunkSize;
            const uint64_t end = (t == threadCount - 1) ? readCount : min<uint64_t>((t + 1) * chunkSize, readCount);
            Timing& timing = timings[t];

            // Per-thread scratch pad reused across reads.
            HifiasmECScratchPad scratch;

            // Temp per-site accumulators (reused across reads).
            struct SiteAcc {
                uint32_t occ_0 = 1;       // starts at 1 (query itself counts as ref)
                uint32_t fwd_ref_cov = 1;  // starts at 1 (query is fwd)
                std::array<uint32_t, 4> altCount{};
            };
            struct TempEv {
                uint32_t overlapID;
                uint32_t ri;       // index into querySitesByPos
                uint8_t  type;     // 0=ref, 1=alt
                uint8_t  misBase;  // alt base index (0-3); 4 = unknown
            };
            vector<SiteAcc> siteAcc;
            vector<TempEv>  tempEvidence;
            vector<int32_t> siteToRow;
            vector<uint8_t> siteDominantAlt;

            for (uint64_t readId = start; readId < end; readId++) {
                timing.readsVisited++;

                const OrientedReadId orientedReadId(ReadId(readId), 0);
                if (orientedReadId.getValue() >= alignmentTable.size()) {
                    continue;
                }
                const span<const uint32_t> alignments = alignmentTable[orientedReadId.getValue()];
                if (alignments.empty()) {
                    continue;
                }
                timing.readsWithAlignments++;

                // ------------------------------------------------------------------
                // Phase 1: gather candidates (unchanged from old code).
                // ------------------------------------------------------------------
                const auto tGatherBegin = steady_clock::now();
                scratch.clear();
                scratch.candidates.reserve(alignments.size());

                for (const uint32_t alignmentId : alignments) {
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if (!ad.keptByBothSides()) {
                        continue;
                    }

                    OrientedReadId o0(ad.readIds[0], 0);
                    OrientedReadId o1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
                    AlignmentInfo orientedInfo = ad.info;
                    if (o0.getReadId() != orientedReadId.getReadId()) {
                        std::swap(o0, o1);
                        orientedInfo.swap();
                    }
                    DINARA_ASSERT(o0.getReadId() == orientedReadId.getReadId());

                    if (o0.getStrand() != orientedReadId.getStrand()) {
                        o0.flipStrand();
                        o1.flipStrand();
                        orientedInfo.reverseComplement();
                    }
                    DINARA_ASSERT(o0 == orientedReadId);

                    uint32_t qsCore = 0;
                    uint32_t qeCore = 0;
                    uint32_t tsCoreOriented = 0;
                    uint32_t teCoreOriented = 0;
                    bool fromMarkers = false;
                    const ReadId queryReadId = ReadId(readId);
                    const uint32_t kmerLen = assemblerInfo.isOpen ? uint32_t(assemblerInfo->k) : 0U;

                    if (markers && assemblerInfo.isOpen) {
                        const auto m0 = (*markers)[o0.getValue()];
                        const auto m1 = (*markers)[o1.getValue()];
                        if (!m0.empty() && !m1.empty()) {
                            qsCore = m0[orientedInfo.data[0].firstOrdinal].position;
                            qeCore = m0[orientedInfo.data[0].lastOrdinal].position + kmerLen;
                            tsCoreOriented = m1[orientedInfo.data[1].firstOrdinal].position;
                            teCoreOriented = m1[orientedInfo.data[1].lastOrdinal].position + kmerLen;
                            fromMarkers = true;
                        }
                    }
                    if (!fromMarkers) {
                        if (ad.readIds[0] == queryReadId) {
                            qsCore = ad.qs;
                            qeCore = ad.qe;
                            tsCoreOriented = ad.ts;
                            teCoreOriented = ad.te;
                        } else {
                            qsCore = ad.ts;
                            qeCore = ad.te;
                            tsCoreOriented = ad.qs;
                            teCoreOriented = ad.qe;
                        }
                    }

                    const uint32_t tLen = uint32_t(reads->getRead(o1.getReadId()).baseCount);
                    uint32_t tsFwd = tsCoreOriented;
                    uint32_t teFwd = teCoreOriented;
                    if (fromMarkers && o1.getStrand() != 0) {
                        tsFwd = tLen - teCoreOriented;
                        teFwd = tLen - tsCoreOriented;
                    }

                    CandidateEC ce;
                    ce.alignmentId = alignmentId;
                    ce.qs = qsCore;
                    ce.qe = qeCore;
                    ce.ts = tsFwd;
                    ce.te = teFwd;
                    ce.targetId = uint32_t(o1.getReadId());
                    ce.isRev = (o1.getStrand() != 0);
                    ce.is_match = 1;
                    ce.strong = 0;
                    scratch.candidates.push_back(ce);
                }
                timing.gatherCandidates += seconds(steady_clock::now() - tGatherBegin);

                if (scratch.candidates.empty()) {
                    continue;
                }
                timing.readsWithCandidates++;

                // ------------------------------------------------------------------
                // Phase 2: build snpStats + hapEvidence from global het sites.
                //
                // For each global het site at position p (in query's forward coords):
                //   - Query the pre-computed SNP stream for each candidate that covers p.
                //   - A hit means the candidate carries an alt base → type=1, misBase=base.
                //   - No hit means concordant → type=0.
                //   - After iterating all candidates, pick the dominant alt base and emit
                //     one SnpStats row per site (occ_1 >= 2 required, same as detectHetSites).
                //   - fwd_ref_cov = concordant FW-strand overlaps + 1 (for the query itself).
                // ------------------------------------------------------------------
                const auto tPhaseBegin = steady_clock::now();

                const vector<GlobalReadSite>& querySitesByPos = globalIndex.sitesByReadPos[readId];
                const bool isInformativeRead = querySitesByPos.size() >= minInformativeSitesPerRead;

                if (querySitesByPos.empty() || !isInformativeRead) {
                    // No het sites for this read; nothing to phase — write back cis.
                    const ReadId queryReadId = ReadId(readId);
                    for (const CandidateEC& ce : scratch.candidates) {
                        alignmentData[ce.alignmentId].clearDeleteReasonsFromReadPerspective(
                            queryReadId, AlignmentData::DeleteReasonPhase);
                    }
                    continue;
                }

                const size_t nSites = querySitesByPos.size();
                const size_t nCands = scratch.candidates.size();

                siteAcc.assign(nSites, SiteAcc{});
                tempEvidence.clear();
                tempEvidence.reserve(nSites * 4);

                for (size_t c = 0; c < nCands; c++) {
                    const CandidateEC& ce = scratch.candidates[c];
                    if (ce.qe <= ce.qs) continue;

                    const AlignmentData& ad = alignmentData[ce.alignmentId];
                    const size_t evidenceId = ad.info.alignmentId;
                    if (evidenceId == invalid<size_t>) continue;
                    const uint32_t evidenceId32 = uint32_t(evidenceId);
                    const bool rIsQuery = (ad.readIds[0] == ReadId(readId));

                    const auto qBeginIt = lower_bound(
                        querySitesByPos.begin(), querySitesByPos.end(),
                        uint32_t(ce.qs),
                        [](const GlobalReadSite& a, uint32_t pos) { return a.position < pos; });
                    const auto qEndIt = lower_bound(
                        querySitesByPos.begin(), querySitesByPos.end(),
                        uint32_t(ce.qe),
                        [](const GlobalReadSite& a, uint32_t pos) { return a.position < pos; });

                    for (auto it = qBeginIt; it != qEndIt; ++it) {
                        const size_t ri = size_t(it - querySitesByPos.begin());
                        const uint32_t p = it->position;

                        uint8_t altBase = 4;
                        bool hasMismatch = false;
                        auto snpCb = [&](uint32_t, uint8_t base) {
                            hasMismatch = true;
                            altBase = base;
                        };
                        if (rIsQuery) {
                            alignedEvidenceStore.forEachSnp1InRange(evidenceId32, p, p + 1, snpCb);
                        } else {
                            alignedEvidenceStore.forEachSnp0InRange(evidenceId32, p, p + 1, snpCb);
                        }

                        if (hasMismatch && altBase < 4) {
                            siteAcc[ri].altCount[altBase]++;
                            tempEvidence.push_back({uint32_t(c), uint32_t(ri), 1, altBase});
                        } else {
                            siteAcc[ri].occ_0++;
                            if (!ce.isRev) siteAcc[ri].fwd_ref_cov++;
                            tempEvidence.push_back({uint32_t(c), uint32_t(ri), 0, 4});
                        }
                    }
                }

                // Build SnpStats rows (one per site with dominant alt occ_1 >= 2).
                siteToRow.assign(nSites, -1);
                siteDominantAlt.assign(nSites, 4);
                scratch.snpStats.clear();

                for (size_t ri = 0; ri < nSites; ri++) {
                    const SiteAcc& acc = siteAcc[ri];
                    uint8_t bestBase = 4;
                    uint32_t bestCount = 0;
                    for (uint8_t b = 0; b < 4; b++) {
                        if (acc.altCount[b] > bestCount) {
                            bestCount = acc.altCount[b];
                            bestBase = b;
                        }
                    }
                    if (bestCount < 2) continue;  // require occ_1 >= 2

                    siteToRow[ri] = int32_t(scratch.snpStats.size());
                    siteDominantAlt[ri] = bestBase;

                    SnpStats stat;
                    stat.site        = querySitesByPos[ri].position;
                    stat.occ_0       = acc.occ_0;
                    stat.occ_1       = bestCount;
                    stat.fwd_ref_cov = acc.fwd_ref_cov;
                    stat.refBase     = 'N';  // not used by gen_rphase_dp
                    stat.altBase     = "ACGT"[bestBase];
                    stat.is_homopolymer = 0;
                    stat.score       = -1;
                    stat.dpScore     = 0;
                    scratch.snpStats.push_back(stat);
                }

                if (scratch.snpStats.empty()) {
                    // No sites with enough alt support; write back cis.
                    const ReadId queryReadId = ReadId(readId);
                    for (const CandidateEC& ce : scratch.candidates) {
                        alignmentData[ce.alignmentId].clearDeleteReasonsFromReadPerspective(
                            queryReadId, AlignmentData::DeleteReasonPhase);
                    }
                    continue;
                }

                // Build hapEvidence remapped to row indices.
                scratch.hapEvidence.clear();
                scratch.hapEvidence.reserve(tempEvidence.size());
                for (const TempEv& te : tempEvidence) {
                    const int32_t row = siteToRow[te.ri];
                    if (row < 0) continue;
                    // For alt entries keep only the dominant base.
                    if (te.type == 1 && te.misBase != siteDominantAlt[te.ri]) continue;

                    HaplotypeEvidence ev;
                    ev.overlapID   = te.overlapID;
                    ev.site        = querySitesByPos[te.ri].position;
                    ev.overlapSite = uint32_t(row);
                    ev.type        = te.type;
                    ev.misBase     = te.misBase < 4 ? te.misBase : 0;
                    ev.hp          = false;
                    scratch.hapEvidence.push_back(ev);
                }
                sort(scratch.hapEvidence.begin(), scratch.hapEvidence.end());

                // No insertion holes — DP will skip hole processing.
                scratch.insertionOffsets.clear();
                scratch.insertionIntervals.clear();
                scratch.insertionBaseCount.clear();

                timing.phaseUsingGlobalSites += seconds(steady_clock::now() - tPhaseBegin);

                // ------------------------------------------------------------------
                // Phase 3: run the hifiasm DP + transitive-closure pipeline.
                // ------------------------------------------------------------------
                gen_rphase_dp(*this, scratch);
                generate_haplotypes_naive_HiFi(*this, scratch);

                // ------------------------------------------------------------------
                // Phase 4: write back DeleteReasonPhase flags.
                // ------------------------------------------------------------------
                const auto tFinalizeBegin = steady_clock::now();
                const ReadId queryReadId = ReadId(readId);
                for (size_t c = 0; c < scratch.candidates.size(); c++) {
                    const CandidateEC& ce = scratch.candidates[c];
                    AlignmentData& ad = alignmentData[ce.alignmentId];

                    const bool keep = (ce.is_match == 1);
                    if (keep) {
                        ad.clearDeleteReasonsFromReadPerspective(queryReadId, AlignmentData::DeleteReasonPhase);
                    } else {
                        ad.addDeleteReasonsFromReadPerspective(queryReadId, AlignmentData::DeleteReasonPhase);
                    }

                    // Track how many DP-retained sites were shared with this candidate.
                    const uint32_t informativeCount = uint32_t(scratch.snpStats.size());
                    if (ad.readIds[0] == queryReadId) {
                        ad.informativeHetSiteCount0 = informativeCount;
                    } else if (ad.readIds[1] == queryReadId) {
                        ad.informativeHetSiteCount1 = informativeCount;
                    }
                }
                timing.finalizeFlags += seconds(steady_clock::now() - tFinalizeBegin);
            }
        });
    }

    for (auto& th : threads) {
        th.join();
    }

    // =========================================================================
    // Transitive phasing via BFS 2-coloring of the phasing graph.
    //
    // The per-thread direct pass labeled each alignment cis or trans from each
    // read's perspective via DeleteReasonPhase flags. Here we build a graph
    // over reads with those edge labels and propagate haplotype assignments
    // (color 0 / color 1) by BFS so that reads never directly aligned can
    // still be phased through intermediate reads.
    //
    // Edge label from a read's perspective (using that side's flags):
    //   isTrans{0,1}    — read voted this alignment trans (DeleteReasonPhase set)
    //   hasCisEvidence  — voted cis with enough shared sites (not trans,
    //                     but informativeHetSiteCount >= minSharedSitesForDecision)
    //   unknown         — insufficient shared sites → skip edge in BFS
    //
    // After 2-coloring, apply the color result to every alignment where both
    // endpoints have a consistent (non-contradicted) color:
    //   same color  → cis  → clear DeleteReasonPhase from both sides
    //   diff color  → trans → set   DeleteReasonPhase on both sides
    // =========================================================================
    const auto tBfsBegin = steady_clock::now();
    {
        // Whether an alignment was deleted for reasons OTHER than phasing
        // (containment, chemical, etc.). Such alignments are excluded from the
        // BFS entirely; they were not candidates in the direct pass and must
        // not participate in color propagation.
        auto deletedForNonPhaseReasons = [](const AlignmentData& ad) -> bool {
            constexpr AlignmentData::DeleteReasonMask nonPhase =
                ~AlignmentData::DeleteReasonPhase;
            return ((ad.deleteReasons0 & nonPhase) != 0) ||
                   ((ad.deleteReasons1 & nonPhase) != 0);
        };

        // Edge label helpers from each side's perspective.
        auto isTrans0 = [](const AlignmentData& ad) -> bool {
            return (ad.deleteReasons0 & AlignmentData::DeleteReasonPhase) != 0;
        };
        auto hasCisEvidence0 = [&](const AlignmentData& ad) -> bool {
            return !isTrans0(ad) &&
                   (ad.informativeHetSiteCount0 >= minSharedSitesForDecision);
        };
        auto isTrans1 = [](const AlignmentData& ad) -> bool {
            return (ad.deleteReasons1 & AlignmentData::DeleteReasonPhase) != 0;
        };
        auto hasCisEvidence1 = [&](const AlignmentData& ad) -> bool {
            return !isTrans1(ad) &&
                   (ad.informativeHetSiteCount1 >= minSharedSitesForDecision);
        };

        static constexpr uint8_t NO_COLOR = 255;
        vector<uint8_t> readColor(readCount, NO_COLOR);
        vector<bool>    inconsistent(readCount, false);

        // BFS using a vector as a queue (index-based, safe across push_back).
        // Each entry stores (readId, depth).  We propagate at most 2 hops
        // (depth 0 = seed, depth 1 = direct neighbours, depth 2 = neighbours
        // of neighbours).  Reads at depth 2 get a colour but do NOT enqueue
        // their own neighbours, so nothing beyond hop 2 is coloured.
        vector<pair<uint32_t, uint8_t>> bfsQueue;
        bfsQueue.reserve(4096);

        for (uint32_t seed = 0; seed < uint32_t(readCount); ++seed) {
            if (readColor[seed] != NO_COLOR || inconsistent[seed]) continue;

            // Only start a component from a read that has at least one labeled
            // edge; reads with no phasing signal stay uncolored.
            const OrientedReadId seedOriented(ReadId(seed), 0);
            if (seedOriented.getValue() >= alignmentTable.size()) continue;
            const span<const uint32_t> seedAlns =
                alignmentTable[seedOriented.getValue()];

            bool hasLabeledEdge = false;
            for (const uint32_t alnId : seedAlns) {
                const AlignmentData& ad = alignmentData[alnId];
                if (deletedForNonPhaseReasons(ad)) continue;
                const bool seedIsR0 = (ad.readIds[0] == ReadId(seed));
                if (seedIsR0 ? (isTrans0(ad) || hasCisEvidence0(ad))
                             : (isTrans1(ad) || hasCisEvidence1(ad))) {
                    hasLabeledEdge = true;
                    break;
                }
            }
            if (!hasLabeledEdge) continue;

            // Seed this component with color 0 and BFS outward (max 2 hops).
            readColor[seed] = 0;
            bfsQueue.clear();
            bfsQueue.push_back({seed, 0});

            for (size_t qi = 0; qi < bfsQueue.size(); ++qi) {
                const uint32_t cur   = bfsQueue[qi].first;
                const uint8_t  depth = bfsQueue[qi].second;
                if (inconsistent[cur]) continue;
                const uint8_t curColor = readColor[cur];

                // At depth 2 we already have the neighbour-of-neighbour
                // coloured; don't go any deeper.
                if (depth >= 2) continue;

                const OrientedReadId curOriented(ReadId(cur), 0);
                if (curOriented.getValue() >= alignmentTable.size()) continue;
                const span<const uint32_t> curAlns =
                    alignmentTable[curOriented.getValue()];

                for (const uint32_t alnId : curAlns) {
                    const AlignmentData& ad = alignmentData[alnId];
                    if (deletedForNonPhaseReasons(ad)) continue;

                    const bool curIsR0 = (ad.readIds[0] == ReadId(cur));
                    const uint32_t neighbor = curIsR0 ?
                        uint32_t(ad.readIds[1]) : uint32_t(ad.readIds[0]);

                    // Determine edge label from cur's perspective.
                    bool edgeIsTrans, edgeHasLabel;
                    if (curIsR0) {
                        edgeIsTrans  = isTrans0(ad);
                        edgeHasLabel = edgeIsTrans || hasCisEvidence0(ad);
                    } else {
                        edgeIsTrans  = isTrans1(ad);
                        edgeHasLabel = edgeIsTrans || hasCisEvidence1(ad);
                    }
                    if (!edgeHasLabel) continue;

                    // trans edge flips the color; cis edge preserves it.
                    const uint8_t neighborExpected =
                        edgeIsTrans ? uint8_t(1u - curColor) : curColor;

                    if (readColor[neighbor] == NO_COLOR) {
                        readColor[neighbor] = neighborExpected;
                        // Only enqueue if we haven't yet reached depth 2.
                        bfsQueue.push_back({neighbor, uint8_t(depth + 1)});
                    } else if (readColor[neighbor] != neighborExpected) {
                        // Contradiction: two paths assign different colors.
                        // Mark the read inconsistent; do not change its color
                        // and do not propagate further from it.
                        inconsistent[neighbor] = true;
                    }
                }
            }
        }

        // Apply transitive coloring: override direct-pass flags for every
        // alignment where both endpoints have a consistent, non-NO_COLOR color.
        uint64_t transitivelyPhased = 0;
        for (AlignmentData& ad : alignmentData) {
            if (deletedForNonPhaseReasons(ad)) continue;
            const uint32_t r0 = uint32_t(ad.readIds[0]);
            const uint32_t r1 = uint32_t(ad.readIds[1]);
            if (readColor[r0] == NO_COLOR || readColor[r1] == NO_COLOR) continue;
            if (inconsistent[r0] || inconsistent[r1]) continue;

            ++transitivelyPhased;
            if (readColor[r0] == readColor[r1]) {
                // Cis: same haplotype — clear trans flags from both sides.
                ad.clearDeleteReasonsBoth(AlignmentData::DeleteReasonPhase);
            } else {
                // Trans: opposite haplotypes — set trans flags on both sides.
                ad.addDeleteReasonsBoth(AlignmentData::DeleteReasonPhase);
            }
        }

        uint64_t colored = 0, inconsistentCount = 0;
        for (uint32_t r = 0; r < uint32_t(readCount); ++r) {
            if (readColor[r] != NO_COLOR) ++colored;
            if (inconsistent[r])          ++inconsistentCount;
        }
        const double tBfs = seconds(steady_clock::now() - tBfsBegin);
        cout << timestamp << "  transitive 2-color: " << tBfs << " s"
             << " (colored=" << colored
             << ", inconsistent=" << inconsistentCount
             << ", transitivelyPhased=" << transitivelyPhased << ")" << endl;
    }

    for (AlignmentData& ad : alignmentData) {
        ad.updateInformativeHetSiteScore();
    }

    Timing total;
    for (const auto& t : timings) {
        total.gatherCandidates += t.gatherCandidates;
        total.phaseUsingGlobalSites += t.phaseUsingGlobalSites;
        total.finalizeFlags += t.finalizeFlags;
        total.readsVisited += t.readsVisited;
        total.readsWithAlignments += t.readsWithAlignments;
        total.readsWithCandidates += t.readsWithCandidates;
    }

    const double tAll = seconds(steady_clock::now() - tBeginAll);
    cout << timestamp << "Global-site EC timings:" << endl;
    cout << timestamp << "  build global sites: " << tBuildSites << " s"
         << " (clusters=" << clusters.clusterRepresentatives.size()
         << ", keptSites=" << globalIndex.keptSiteCount << ")" << endl;
    cout << timestamp << "  gather candidates:  " << total.gatherCandidates << " s" << endl;
    cout << timestamp << "  phase from globals: " << total.phaseUsingGlobalSites << " s" << endl;
    cout << timestamp << "  finalize flags:     " << total.finalizeFlags << " s" << endl;
    cout << timestamp << "Global-site EC wall time: " << tAll << " s"
         << " (reads=" << total.readsVisited
         << ", withAlignments=" << total.readsWithAlignments
         << ", withCandidates=" << total.readsWithCandidates << ")" << endl;
}
