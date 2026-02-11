// Dinara.
#include "Assembler.hpp"
#include "Reads.hpp"
#include "chrono.hpp"
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
    vector< vector<GlobalReadSite> > sitesByReadPos;
    vector< vector<GlobalReadSite> > sitesByReadSite;
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
    out.sitesByReadSite.assign(readCount, vector<GlobalReadSite>{});

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

        auto& bySite = out.sitesByReadSite[rid];
        bySite = byPos;
        sort(bySite.begin(), bySite.end(), siteOrder);
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

    const bool debugCheckBases = (::getenv("DINARA_EC_GLOBALSITE_ASSERT_BASES") != nullptr);

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

                vector<CandidateGlobal> candidates;
                candidates.reserve(alignments.size());

                const auto tGatherBegin = steady_clock::now();
                for (const uint32_t alignmentId : alignments) {
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if (ad.isDeleted()) {
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

                    CandidateGlobal c;
                    c.alignmentId = alignmentId;
                    c.qs = qsCore;
                    c.qe = qeCore;
                    c.ts = tsFwd;
                    c.te = teFwd;
                    c.targetId = uint32_t(o1.getReadId());
                    c.targetIsRc = (o1.getStrand() != 0);
                    c.isMatch = 1;
                    candidates.push_back(c);
                }
                timing.gatherCandidates += seconds(steady_clock::now() - tGatherBegin);

                if (candidates.empty()) {
                    continue;
                }
                timing.readsWithCandidates++;

                const vector<GlobalReadSite>& querySitesByPos = globalIndex.sitesByReadPos[readId];
                const bool isInformativeRead = querySitesByPos.size() >= minInformativeSitesPerRead;
                vector<uint32_t> informativeCountByCandidate(candidates.size(), 0);

                const auto tPhaseBegin = steady_clock::now();
                for (size_t c = 0; c < candidates.size(); c++) {
                    CandidateGlobal& candidate = candidates[c];
                    if (candidate.qe <= candidate.qs || candidate.te <= candidate.ts) {
                        continue;
                    }
                    if (querySitesByPos.empty()) {
                        continue;
                    }
                    if (candidate.targetId >= globalIndex.sitesByReadSite.size()) {
                        continue;
                    }
                    const auto& targetSitesBySite = globalIndex.sitesByReadSite[candidate.targetId];
                    if (targetSitesBySite.empty()) {
                        continue;
                    }

                    const auto qBeginIt = lower_bound(
                        querySitesByPos.begin(),
                        querySitesByPos.end(),
                        candidate.qs,
                        [](const GlobalReadSite& a, uint32_t pos) { return a.position < pos; });
                    const auto qEndIt = lower_bound(
                        querySitesByPos.begin(),
                        querySitesByPos.end(),
                        candidate.qe,
                        [](const GlobalReadSite& a, uint32_t pos) { return a.position < pos; });

                    uint32_t sharedSites = 0;
                    uint32_t concordantAlleles = 0;
                    uint32_t discordantAlleles = 0;

                    static const uint8_t complementBase[4] = {3, 2, 1, 0};
                    for (auto it = qBeginIt; it != qEndIt; ++it) {
                        const uint32_t siteId = it->siteId;
                        const auto tIt = lower_bound(
                            targetSitesBySite.begin(),
                            targetSitesBySite.end(),
                            siteId,
                            [](const GlobalReadSite& a, uint32_t s) { return a.siteId < s; });
                        if (tIt == targetSitesBySite.end() || tIt->siteId != siteId) {
                            continue;
                        }
                        if (tIt->position < candidate.ts || tIt->position >= candidate.te) {
                            continue;
                        }

                        if (debugCheckBases) {
                            const uint8_t qBase = reads->getOrientedReadBase(
                                OrientedReadId(ReadId(readId), 0), it->position).value;
                            if (qBase < 4) {
                                DINARA_ASSERT(qBase == it->allele);
                            }
                            const uint8_t tBase = reads->getOrientedReadBase(
                                OrientedReadId(ReadId(candidate.targetId), 0), tIt->position).value;
                            if (tBase < 4) {
                                DINARA_ASSERT(tBase == tIt->allele);
                            }
                        }

                        sharedSites++;
                        uint8_t targetAllele = tIt->allele;
                        if (candidate.targetIsRc && targetAllele < 4) {
                            targetAllele = complementBase[targetAllele];
                        }
                        if (targetAllele == it->allele) {
                            concordantAlleles++;
                        } else {
                            discordantAlleles++;
                        }
                    }

                    informativeCountByCandidate[c] = sharedSites;
                    if (
                        isInformativeRead &&
                        sharedSites >= minSharedSitesForDecision &&
                        discordantAlleles > concordantAlleles) {
                        candidate.isMatch = 2;
                    }
                }
                timing.phaseUsingGlobalSites += seconds(steady_clock::now() - tPhaseBegin);

                const auto tFinalizeBegin = steady_clock::now();
                const ReadId queryReadId = ReadId(readId);
                for (size_t c = 0; c < candidates.size(); c++) {
                    const CandidateGlobal& candidate = candidates[c];
                    AlignmentData& ad = alignmentData[candidate.alignmentId];

                    const bool keep = !isInformativeRead || (candidate.isMatch == 1);
                    if (keep) {
                        ad.clearDeleteReasonsFromReadPerspective(queryReadId, AlignmentData::DeleteReasonPhase);
                    } else {
                        ad.addDeleteReasonsFromReadPerspective(queryReadId, AlignmentData::DeleteReasonPhase);
                    }

                    const uint32_t informativeCount = informativeCountByCandidate[c];
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
