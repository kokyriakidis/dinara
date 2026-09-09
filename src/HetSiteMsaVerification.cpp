#include "HetSiteMsaVerification.hpp"

#include "Reads.hpp"
#include "timestamp.hpp"

#include "abpoa/abpoa.h"

#include <algorithm>
#include <atomic>
#include <iostream>
#include <mutex>
#include <thread>
#include <unordered_map>

using namespace dinara;
using std::cout;
using std::endl;



// ============================================================================
// The pure half.
// ============================================================================

MsaSiteVerdict dinara::verifyColumnAgreement(
    const uint8_t* const* msa,
    int nSeq,
    int msaLen,
    int gapValue,
    const vector<uint32_t>& snpOffsetInRow,
    double minAgreementFraction)
{
    MsaSiteVerdict verdict;
    if(nSeq < 2 || msaLen <= 0 || int(snpOffsetInRow.size()) != nSeq) {
        verdict.reason = MsaSiteVerdict::Reason::tooFewMembers;
        return verdict;
    }

    // Walk each row, counting non-gap entries, to find the column holding that
    // row's SNP base. A row whose sequence is shorter than its recorded offset
    // (possible if the interval was clipped) simply does not vote.
    vector<int> columnOfRow(nSeq, -1);
    for(int row = 0; row < nSeq; row++) {
        uint32_t seen = 0;
        for(int c = 0; c < msaLen; c++) {
            if(int(msa[row][c]) == gapValue) {
                continue;
            }
            if(seen == snpOffsetInRow[size_t(row)]) {
                columnOfRow[size_t(row)] = c;
                break;
            }
            ++seen;
        }
    }

    // Majority column among the rows that placed.
    std::unordered_map<int, uint32_t> votes;
    for(const int c : columnOfRow) {
        if(c >= 0) {
            ++votes[c];
        }
    }
    if(votes.empty()) {
        verdict.reason = MsaSiteVerdict::Reason::msaFailed;
        return verdict;
    }
    int bestColumn = -1;
    uint32_t bestVotes = 0;
    for(const auto& [c, n] : votes) {
        // Ties broken on the smaller column so the result does not depend on
        // hash order -- the same determinism requirement as everywhere else.
        if(n > bestVotes || (n == bestVotes && c < bestColumn)) {
            bestVotes = n;
            bestColumn = c;
        }
    }

    uint32_t placed = 0;
    for(const int c : columnOfRow) {
        if(c >= 0) {
            ++placed;
        }
    }

    verdict.column = bestColumn;
    verdict.agreeing = bestVotes;
    verdict.total = placed;

    if(double(bestVotes) < minAgreementFraction * double(placed)) {
        verdict.reason = MsaSiteVerdict::Reason::columnsDisagree;
        return verdict;
    }

    // The agreed column must still look het: at least two distinct bases with
    // at least 2 supporting rows. A column that collapses to one allele under
    // the MSA is exactly the artifact this pass exists to catch.
    uint32_t baseCount[4] = {0, 0, 0, 0};
    for(int row = 0; row < nSeq; row++) {
        const uint8_t v = msa[row][bestColumn];
        if(int(v) == gapValue || v >= 4) {
            continue;
        }
        ++baseCount[v];
    }
    uint32_t alleles = 0;
    for(const uint32_t n : baseCount) {
        if(n >= 2) {
            ++alleles;
        }
    }
    verdict.allelesAtColumn = alleles;
    if(alleles < 2) {
        verdict.reason = MsaSiteVerdict::Reason::notBiallelic;
        return verdict;
    }

    verdict.verified = true;
    verdict.reason = MsaSiteVerdict::Reason::ok;
    return verdict;
}



// ============================================================================
// The driver.
// ============================================================================

namespace {

// (position, anchorId) per oriented read, ascending by position. Built once and
// shared by every thread.
using AnchorPositions = vector<vector<std::pair<uint32_t, Shasta2AnchorId>>>;

AnchorPositions buildAnchorPositions(
    const Shasta2Anchors& anchors, uint64_t orientedReadCount)
{
    AnchorPositions byOrientedRead(orientedReadCount);
    const uint64_t anchorCount = anchors.size();
    for(Shasta2AnchorId id = 0; id < anchorCount; id++) {
        for(const Shasta2AnchorMarkerInfo& mi : anchors[id]) {
            const uint64_t v = mi.orientedReadId.getValue();
            if(v < orientedReadCount) {
                byOrientedRead[v].push_back({mi.position, id});
            }
        }
    }
    for(auto& v : byOrientedRead) {
        std::sort(v.begin(), v.end());
    }
    return byOrientedRead;
}

// Anchor ids within [position - flank, position) on this oriented read.
void anchorsBefore(
    const vector<std::pair<uint32_t, Shasta2AnchorId>>& sorted,
    uint32_t position, uint32_t flank,
    vector<Shasta2AnchorId>& out)
{
    out.clear();
    const uint32_t low = (position > flank) ? (position - flank) : 0u;
    auto it = std::lower_bound(sorted.begin(), sorted.end(),
        std::make_pair(low, Shasta2AnchorId(0)));
    for(; it != sorted.end() && it->first < position; ++it) {
        out.push_back(it->second);
    }
}

// Anchor ids within (position, position + flank] on this oriented read.
void anchorsAfter(
    const vector<std::pair<uint32_t, Shasta2AnchorId>>& sorted,
    uint32_t position, uint32_t flank,
    vector<Shasta2AnchorId>& out)
{
    out.clear();
    auto it = std::upper_bound(sorted.begin(), sorted.end(),
        std::make_pair(position, std::numeric_limits<Shasta2AnchorId>::max()));
    for(; it != sorted.end() && it->first <= position + flank; ++it) {
        out.push_back(it->second);
    }
}

// Position of anchorId on this oriented read, or invalid if absent.
bool positionOf(
    const vector<std::pair<uint32_t, Shasta2AnchorId>>& sorted,
    Shasta2AnchorId anchorId, uint32_t& positionOut)
{
    for(const auto& [position, id] : sorted) {
        if(id == anchorId) {
            positionOut = position;
            return true;
        }
    }
    return false;
}

// Run abPOA over `sequences` (each already 0..3 encoded, oriented) and hand the
// row-major MSA matrix to `consume`. Returns false if abPOA produced nothing.
template<class Consume>
bool runMsa(const vector<vector<uint8_t>>& sequences, Consume&& consume)
{
    const int nSeq = int(sequences.size());
    if(nSeq < 2) {
        return false;
    }

    abpoa_t* ab = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    abpt->out_msa = 1;
    abpt->out_cons = 0;
    abpt->disable_seeding = 1;
    abpoa_post_set_para(abpt);

    vector<int> lengths(nSeq);
    vector<uint8_t*> seqs(nSeq);
    // abPOA takes non-const pointers but does not modify the input; the copies
    // live in `sequences`, which outlives this call.
    for(int i = 0; i < nSeq; i++) {
        lengths[size_t(i)] = int(sequences[size_t(i)].size());
        seqs[size_t(i)] = const_cast<uint8_t*>(sequences[size_t(i)].data());
    }

    abpoa_msa(ab, abpt, nSeq, nullptr, lengths.data(), seqs.data(),
        nullptr, nullptr);

    bool ok = false;
    abpoa_cons_t* abc = ab->abc;
    if(abc != nullptr && abc->msa_len > 0 && abc->msa_base != nullptr) {
        consume(abc->msa_base, nSeq, int(abc->msa_len), int(abpt->m));
        ok = true;
    }

    abpoa_free(ab);
    abpoa_free_para(abpt);
    return ok;
}

} // namespace



uint64_t dinara::msaVerifyHetSites(
    vector<Assembler::CigarSnpSite>& sites,
    const Reads& reads,
    const Shasta2Anchors& anchors,
    uint32_t flankBases,
    double minAgreementFraction,
    uint64_t threadCount,
    vector<uint64_t>* reasonCountsOut)
{
    if(sites.empty()) {
        return 0;
    }
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    if(threadCount == 0) {
        threadCount = 1;
    }

    const uint64_t orientedReadCount = 2 * reads.readCount();
    const AnchorPositions anchorPositions =
        buildAnchorPositions(anchors, orientedReadCount);

    vector<uint8_t> keep(sites.size(), 1);
    const uint64_t reasonCount = 6;
    vector<vector<uint64_t>> reasonsByThread(threadCount, vector<uint64_t>(reasonCount, 0));

    const uint64_t chunk = (sites.size() + threadCount - 1) / threadCount;
    vector<std::thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t = 0; t < threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = std::min<uint64_t>(sites.size(), (t + 1) * chunk);
            vector<uint64_t>& reasons = reasonsByThread[t];

            // Hoisted scratch: this loop runs tens of thousands of times.
            vector<Shasta2AnchorId> before, after;
            std::unordered_map<Shasta2AnchorId, uint32_t> leftVotes, rightVotes;
            vector<std::pair<OrientedReadId, uint32_t>> members;
            vector<vector<uint8_t>> sequences;
            vector<uint32_t> snpOffset;

            for(uint64_t siteIndex = begin; siteIndex < end; siteIndex++) {
                const Assembler::CigarSnpSite& site = sites[siteIndex];

                members.clear();
                for(const auto& allele : site.alleles) {
                    for(const auto& m : allele) {
                        members.push_back(m);
                    }
                }
                if(members.size() < 2) {
                    keep[siteIndex] = 0;
                    ++reasons[uint64_t(MsaSiteVerdict::Reason::tooFewMembers)];
                    continue;
                }

                // Anchors bracketing the SNP, chosen by QUORUM rather than
                // by unanimity. Requiring every member to share both bounds
                // sounds safer but is not: one read missing one anchor empties
                // the intersection, and measured on E821 that rejected 4207 of
                // 4542 sites for "no shared anchors" while only ONE site
                // actually failed the homology test. That measures the
                // requirement, not the data. So take the anchor each side that
                // the most members carry; members lacking it simply do not
                // participate, exactly like a read clipped by the interval.
                leftVotes.clear();
                rightVotes.clear();
                bool viable = true;
                for(const auto& [orientedReadId, position] : members) {
                    const uint64_t v = orientedReadId.getValue();
                    if(v >= orientedReadCount) { viable = false; break; }
                    const auto& sorted = anchorPositions[v];
                    anchorsBefore(sorted, position, flankBases, before);
                    anchorsAfter(sorted, position, flankBases, after);
                    for(const Shasta2AnchorId id : before) ++leftVotes[id];
                    for(const Shasta2AnchorId id : after) ++rightVotes[id];
                }
                if(!viable || leftVotes.empty() || rightVotes.empty()) {
                    keep[siteIndex] = 0;
                    ++reasons[uint64_t(MsaSiteVerdict::Reason::noSharedAnchors)];
                    continue;
                }

                // Most-carried anchor on each side, ties broken on the smaller
                // id so the choice does not depend on hash order.
                const auto pickBest = [](
                    const std::unordered_map<Shasta2AnchorId, uint32_t>& votes)
                {
                    Shasta2AnchorId best = 0;
                    uint32_t bestVotes = 0;
                    bool have = false;
                    for(const auto& [id, n] : votes) {
                        if(!have || n > bestVotes || (n == bestVotes && id < best)) {
                            best = id; bestVotes = n; have = true;
                        }
                    }
                    return best;
                };
                const Shasta2AnchorId leftAnchor = pickBest(leftVotes);
                const Shasta2AnchorId rightAnchor = pickBest(rightVotes);

                // Extract each member's substring between those anchors, and
                // where its SNP base sits inside it.
                sequences.clear();
                snpOffset.clear();
                for(const auto& [orientedReadId, position] : members) {
                    const auto& sorted = anchorPositions[orientedReadId.getValue()];
                    uint32_t leftPos = 0, rightPos = 0;
                    if(!positionOf(sorted, leftAnchor, leftPos)) continue;
                    if(!positionOf(sorted, rightAnchor, rightPos)) continue;
                    if(!(leftPos <= position && position < rightPos)) continue;

                    vector<uint8_t> seq;
                    seq.reserve(rightPos - leftPos);
                    for(uint32_t p = leftPos; p < rightPos; p++) {
                        seq.push_back(uint8_t(
                            reads.getOrientedReadBase(orientedReadId, p).value));
                    }
                    snpOffset.push_back(position - leftPos);
                    sequences.push_back(std::move(seq));
                }
                if(sequences.size() < 2) {
                    keep[siteIndex] = 0;
                    ++reasons[uint64_t(MsaSiteVerdict::Reason::tooFewMembers)];
                    continue;
                }

                MsaSiteVerdict verdict;
                const bool ran = runMsa(sequences,
                    [&](uint8_t** msa, int nSeq, int msaLen, int gapValue) {
                        verdict = verifyColumnAgreement(
                            msa, nSeq, msaLen, gapValue, snpOffset,
                            minAgreementFraction);
                    });
                if(!ran) {
                    keep[siteIndex] = 0;
                    ++reasons[uint64_t(MsaSiteVerdict::Reason::msaFailed)];
                    continue;
                }
                ++reasons[uint64_t(verdict.reason)];
                if(!verdict.verified) {
                    keep[siteIndex] = 0;
                }
            }
        });
    }
    for(std::thread& thread : threads) {
        thread.join();
    }

    vector<uint64_t> reasons(reasonCount, 0);
    for(const auto& r : reasonsByThread) {
        for(uint64_t i = 0; i < reasonCount; i++) {
            reasons[i] += r[i];
        }
    }
    if(reasonCountsOut != nullptr) {
        *reasonCountsOut = reasons;
    }

    vector<Assembler::CigarSnpSite> kept;
    kept.reserve(sites.size());
    for(uint64_t i = 0; i < sites.size(); i++) {
        if(keep[i]) {
            kept.push_back(std::move(sites[i]));
        }
    }
    const uint64_t rejected = sites.size() - kept.size();
    sites.swap(kept);
    return rejected;
}
