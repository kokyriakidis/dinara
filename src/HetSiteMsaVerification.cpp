#include "HetSiteMsaVerification.hpp"

#include "Reads.hpp"
#include "span.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"

#include "abpoa/abpoa.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <array>
#include <limits>
#include <thread>
#include <unordered_map>

using namespace dinara;

// Below this many placed rows the biallelic and partition tests have too little
// to work with: two alleles need >= 2 rows each, so anything under 4 can only
// ever fail or pass for want of evidence rather than on the evidence.
static constexpr uint32_t defaultMinPlacedRows = 4;



const char* dinara::msaVerdictReasonName(MsaSiteVerdict::Reason reason)
{
    switch(reason) {
    case MsaSiteVerdict::Reason::ok:                return "confirmed";
    case MsaSiteVerdict::Reason::noSharedAnchors:   return "no shared bracketing anchors";
    case MsaSiteVerdict::Reason::tooFewMembers:     return "too few members placed";
    case MsaSiteVerdict::Reason::msaFailed:         return "MSA failed";
    case MsaSiteVerdict::Reason::columnsDisagree:   return "members landed in different columns";
    case MsaSiteVerdict::Reason::notBiallelic:      return "agreed column not biallelic";
    case MsaSiteVerdict::Reason::partitionMismatch: return "column splits reads differently than the arms";
    }
    return "unknown";
}



// ============================================================================
// The pure half.
// ============================================================================

MsaSiteVerdict dinara::verifyColumnAgreement(
    const uint8_t* const* msa,
    int nSeq,
    int msaLen,
    int gapValue,
    const vector<uint32_t>& snpOffsetInRow,
    const vector<uint8_t>& armOfRow,
    double minAgreementFraction,
    uint32_t minPlacedRows)
{
    MsaSiteVerdict verdict;
    if(nSeq < 2 || msaLen <= 0 ||
        int(snpOffsetInRow.size()) != nSeq || int(armOfRow.size()) != nSeq) {
        verdict.reason = MsaSiteVerdict::Reason::tooFewMembers;
        return verdict;
    }

    // Walk each row, counting NON-GAP entries, to find the column holding that
    // row's SNP base. A row whose sequence is shorter than its recorded offset
    // (the interval clipped it) simply does not place, and does not vote.
    vector<int> columnOfRow(nSeq, -1);
    uint32_t placed = 0;
    for(int row = 0; row < nSeq; row++) {
        uint32_t seen = 0;
        for(int c = 0; c < msaLen; c++) {
            if(int(msa[row][c]) == gapValue) {
                continue;
            }
            if(seen == snpOffsetInRow[size_t(row)]) {
                columnOfRow[size_t(row)] = c;
                ++placed;
                break;
            }
            ++seen;
        }
    }
    verdict.total = placed;

    if(placed < minPlacedRows) {
        verdict.reason = MsaSiteVerdict::Reason::tooFewMembers;
        return verdict;
    }

    // Majority column among the rows that placed.
    std::unordered_map<int, uint32_t> votes;
    for(const int c : columnOfRow) {
        if(c >= 0) {
            ++votes[c];
        }
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
    verdict.column = bestColumn;
    verdict.agreeing = bestVotes;

    if(double(bestVotes) < minAgreementFraction * double(placed)) {
        verdict.reason = MsaSiteVerdict::Reason::columnsDisagree;
        return verdict;
    }

    // The agreed column must still look het. Only rows that PLACED vote here:
    // a row whose SNP base sits elsewhere is not evidence about this column,
    // and counting it was letting clipped rows tip the allele census.
    std::array<uint32_t, 4> baseCount{0, 0, 0, 0};
    for(int row = 0; row < nSeq; row++) {
        if(columnOfRow[size_t(row)] < 0) {
            continue;
        }
        const uint8_t v = msa[row][bestColumn];
        if(int(v) == gapValue || v >= 4) {
            continue;
        }
        ++baseCount[v];
    }
    // The two best-supported bases, each needing >= 2 rows.
    int allele0 = -1, allele1 = -1;
    for(int b = 0; b < 4; b++) {
        if(baseCount[size_t(b)] < 2) {
            continue;
        }
        if(allele0 < 0 || baseCount[size_t(b)] > baseCount[size_t(allele0)]) {
            allele1 = allele0;
            allele0 = b;
        } else if(allele1 < 0 || baseCount[size_t(b)] > baseCount[size_t(allele1)]) {
            allele1 = b;
        }
    }
    verdict.allelesAtColumn = uint32_t((allele0 >= 0 ? 1 : 0) + (allele1 >= 0 ? 1 : 0));
    if(allele1 < 0) {
        verdict.reason = MsaSiteVerdict::Reason::notBiallelic;
        return verdict;
    }

    // Does the column split the reads the same way the ARMS do? Both arms can
    // be biallelic and still disagree about which reads carry which allele, and
    // phasing consumes the partition, not the allele count. Score both possible
    // arm->allele bijections and take the better one.
    std::array<std::array<uint32_t, 2>, 2> table{};   // [arm][alleleIndex]
    for(int row = 0; row < nSeq; row++) {
        if(columnOfRow[size_t(row)] < 0) {
            continue;
        }
        const uint8_t arm = armOfRow[size_t(row)];
        if(arm > 1) {
            continue;                       // only the two main arms vote
        }
        const uint8_t v = msa[row][bestColumn];
        if(int(v) == allele0) {
            ++table[arm][0];
        } else if(int(v) == allele1) {
            ++table[arm][1];
        }
    }
    const uint32_t straight = table[0][0] + table[1][1];
    const uint32_t swapped  = table[0][1] + table[1][0];
    verdict.partitionAgreeing = std::max(straight, swapped);
    verdict.partitionTotal =
        table[0][0] + table[0][1] + table[1][0] + table[1][1];
    if(verdict.partitionTotal < minPlacedRows) {
        verdict.reason = MsaSiteVerdict::Reason::tooFewMembers;
        return verdict;
    }
    if(double(verdict.partitionAgreeing) <
        minAgreementFraction * double(verdict.partitionTotal)) {
        verdict.reason = MsaSiteVerdict::Reason::partitionMismatch;
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

// The anchor index: for each oriented read, its anchors as (position,
// anchorId) ascending by position, stored as a flat CSR rather than a vector
// per read.
//
// This is the expensive half of the pass and it does NOT scale with the number
// of sites -- it walks every anchor member in the assembly (16M on a 6 Mb
// fixture), so on a large genome it grows with the genome while the
// verification itself grows only with the site count. Measured single-threaded
// it was 0.52 s of a 1.12 s pass, so it runs on every thread: a counting pass
// sizes each read's slice exactly, a fill pass scatters into it, and the
// per-read sorts are independent.
class AnchorIndex {
public:
    vector<uint64_t> begin_;                                  // size n+1
    vector<std::pair<uint32_t, Shasta2AnchorId>> data_;

    span<const std::pair<uint32_t, Shasta2AnchorId>> operator[](uint64_t v) const
    {
        return span<const std::pair<uint32_t, Shasta2AnchorId>>(
            data_.data() + begin_[v], data_.data() + begin_[v + 1]);
    }
};

AnchorIndex buildAnchorIndex(
    const Shasta2Anchors& anchors, uint64_t orientedReadCount, uint64_t threadCount)
{
    AnchorIndex index;
    index.begin_.assign(orientedReadCount + 1, 0);
    const uint64_t anchorCount = anchors.size();

    const auto runOverAnchors = [&](auto&& body) {
        const uint64_t chunk = (anchorCount + threadCount - 1) / threadCount;
        vector<std::thread> threads;
        threads.reserve(threadCount);
        for(uint64_t t = 0; t < threadCount; t++) {
            threads.emplace_back([&, t]() {
                body(t * chunk, std::min(anchorCount, (t + 1) * chunk));
            });
        }
        for(std::thread& thread : threads) thread.join();
    };

    // Counting pass. One relaxed atomic per oriented read: O(orientedReadCount)
    // rather than O(threadCount * orientedReadCount), which matters once the
    // read count is large.
    {
        vector<std::atomic<uint64_t>> counts(orientedReadCount);
        for(std::atomic<uint64_t>& c : counts) c.store(0, std::memory_order_relaxed);
        runOverAnchors([&](uint64_t begin, uint64_t end) {
            for(Shasta2AnchorId id = begin; id < end; id++) {
                for(const Shasta2AnchorMarkerInfo& mi : anchors[id]) {
                    const uint64_t v = mi.orientedReadId.getValue();
                    if(v < orientedReadCount) {
                        counts[v].fetch_add(1, std::memory_order_relaxed);
                    }
                }
            }
        });
        uint64_t total = 0;
        for(uint64_t v = 0; v < orientedReadCount; v++) {
            index.begin_[v] = total;
            total += counts[v].load(std::memory_order_relaxed);
        }
        index.begin_[orientedReadCount] = total;
        index.data_.resize(total);
    }

    // Fill pass. Entries land inside a read's slice in nondeterministic order,
    // which does not matter: the slice is sorted immediately afterwards.
    {
        vector<std::atomic<uint64_t>> cursor(orientedReadCount);
        for(uint64_t v = 0; v < orientedReadCount; v++) {
            cursor[v].store(index.begin_[v], std::memory_order_relaxed);
        }
        runOverAnchors([&](uint64_t begin, uint64_t end) {
            for(Shasta2AnchorId id = begin; id < end; id++) {
                for(const Shasta2AnchorMarkerInfo& mi : anchors[id]) {
                    const uint64_t v = mi.orientedReadId.getValue();
                    if(v < orientedReadCount) {
                        index.data_[cursor[v].fetch_add(1, std::memory_order_relaxed)] =
                            {mi.position, id};
                    }
                }
            }
        });
    }

    // Sort each read's slice; reads are independent.
    {
        const uint64_t chunk = (orientedReadCount + threadCount - 1) / threadCount;
        vector<std::thread> threads;
        threads.reserve(threadCount);
        for(uint64_t t = 0; t < threadCount; t++) {
            threads.emplace_back([&, t]() {
                const uint64_t begin = t * chunk;
                const uint64_t end = std::min(orientedReadCount, (t + 1) * chunk);
                for(uint64_t v = begin; v < end; v++) {
                    std::sort(index.data_.data() + index.begin_[v],
                              index.data_.data() + index.begin_[v + 1]);
                }
            });
        }
        for(std::thread& thread : threads) thread.join();
    }

    return index;
}

// One candidate bounding anchor, accumulated across the members that carry it.
class Candidate {
public:
    uint32_t votes = 0;
    uint64_t distanceSum = 0;   // summed |anchor position - SNP position|
};

// Rank candidates: most members first, then closest to the SNP (a tighter
// interval is a cheaper alignment and gives ambiguity less room), then smallest
// id so the choice never depends on hash order.
Shasta2AnchorId pickBoundingAnchor(
    const std::unordered_map<Shasta2AnchorId, Candidate>& candidates)
{
    Shasta2AnchorId best = 0;
    const Candidate* bestCandidate = nullptr;
    for(const auto& [id, candidate] : candidates) {
        if(bestCandidate == nullptr) {
            best = id; bestCandidate = &candidate; continue;
        }
        if(candidate.votes != bestCandidate->votes) {
            if(candidate.votes > bestCandidate->votes) { best = id; bestCandidate = &candidate; }
            continue;
        }
        // Compare mean distance without dividing: a*d2 < b*d1 <=> d1/v1 < d2/v2.
        const uint64_t lhs = candidate.distanceSum * uint64_t(bestCandidate->votes);
        const uint64_t rhs = bestCandidate->distanceSum * uint64_t(candidate.votes);
        if(lhs != rhs) {
            if(lhs < rhs) { best = id; bestCandidate = &candidate; }
            continue;
        }
        if(id < best) { best = id; bestCandidate = &candidate; }
    }
    return best;
}

// One abPOA instance, reused across every site a thread handles.
//
// abpoa_init/abpoa_free allocate substantial structures, and doing that per
// site meant four allocations per alignment for thousands of very short
// alignments. abpoa_reset is exactly the supported way to reuse a graph, so
// each thread keeps one of these for its whole batch.
class MsaRunner {
public:
    MsaRunner()
    {
        ab = abpoa_init();
        abpt = abpoa_init_para();
        abpt->out_msa = 1;
        abpt->out_cons = 0;
        abpt->disable_seeding = 1;
        abpoa_post_set_para(abpt);
    }
    ~MsaRunner()
    {
        abpoa_free(ab);
        abpoa_free_para(abpt);
    }
    MsaRunner(const MsaRunner&) = delete;
    MsaRunner& operator=(const MsaRunner&) = delete;

    // Align `sequences` and hand the row-major MSA matrix to `consume`.
    // False if abPOA produced nothing usable.
    template<class Consume>
    bool run(const vector<vector<uint8_t>>& sequences, Consume&& consume)
    {
        const int nSeq = int(sequences.size());
        if(nSeq < 2) {
            return false;
        }

        lengths.resize(size_t(nSeq));
        seqs.resize(size_t(nSeq));
        int maxLength = 0;
        // abPOA takes non-const pointers but does not modify the input; the
        // sequences outlive this call.
        for(int i = 0; i < nSeq; i++) {
            const int length = int(sequences[size_t(i)].size());
            lengths[size_t(i)] = length;
            seqs[size_t(i)] = const_cast<uint8_t*>(sequences[size_t(i)].data());
            maxLength = std::max(maxLength, length);
        }

        abpoa_reset(ab, abpt, maxLength);
        abpoa_msa(ab, abpt, nSeq, nullptr, lengths.data(), seqs.data(),
            nullptr, nullptr);

        const abpoa_cons_t* const abc = ab->abc;
        if(abc == nullptr || abc->msa_len <= 0 || abc->msa_base == nullptr) {
            return false;
        }
        consume(abc->msa_base, nSeq, int(abc->msa_len), int(abpt->m));
        return true;
    }

private:
    abpoa_t* ab = nullptr;
    abpoa_para_t* abpt = nullptr;
    vector<int> lengths;
    vector<uint8_t*> seqs;
};

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

    const auto tBegin = std::chrono::steady_clock::now();
    const uint64_t orientedReadCount = 2 * reads.readCount();
    const AnchorIndex anchorIndex =
        buildAnchorIndex(anchors, orientedReadCount, threadCount);
    const auto tIndexed = std::chrono::steady_clock::now();

    vector<uint8_t> keep(sites.size(), 1);
    const uint64_t reasonCount = MsaSiteVerdict::reasonCount;
    vector<vector<uint64_t>> reasonsByThread(
        threadCount, vector<uint64_t>(reasonCount, 0));

    const uint64_t chunk = (sites.size() + threadCount - 1) / threadCount;
    vector<std::thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t = 0; t < threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = std::min<uint64_t>(sites.size(), (t + 1) * chunk);
            vector<uint64_t>& reasons = reasonsByThread[t];

            // One abPOA instance for this thread's whole batch, plus hoisted
            // scratch: this loop runs tens of thousands of times.
            MsaRunner msaRunner;
            std::unordered_map<Shasta2AnchorId, Candidate> leftCandidates, rightCandidates;
            // Per member, the positions of the candidate anchors on THAT read,
            // so the extraction pass never rescans the read's anchor list.
            vector<std::unordered_map<Shasta2AnchorId, uint32_t>> positionsByMember;
            vector<std::pair<OrientedReadId, uint32_t>> members;
            vector<uint8_t> armOfMember;
            vector<vector<uint8_t>> sequences;
            vector<uint32_t> snpOffset;
            vector<uint8_t> armOfRow;

            for(uint64_t siteIndex = begin; siteIndex < end; siteIndex++) {
                const Assembler::CigarSnpSite& site = sites[siteIndex];

                // Flatten the arms, remembering which arm each member came
                // from: that assignment is what the partition check verifies.
                members.clear();
                armOfMember.clear();
                for(uint64_t arm = 0; arm < site.alleles.size(); arm++) {
                    for(const auto& m : site.alleles[arm]) {
                        members.push_back(m);
                        armOfMember.push_back(uint8_t(std::min<uint64_t>(arm, 255)));
                    }
                }
                if(members.size() < defaultMinPlacedRows) {
                    keep[siteIndex] = 0;
                    ++reasons[uint64_t(MsaSiteVerdict::Reason::tooFewMembers)];
                    continue;
                }

                // Candidate bounding anchors, gathered by QUORUM rather than by
                // unanimity. Requiring every member to share both bounds sounds
                // safer but is not: one read missing one anchor empties the
                // intersection, and measured on E821 that rejected 4207 of 4542
                // sites for "no shared anchors" while exactly ONE failed the
                // homology test -- it measured the requirement, not the data.
                // Members lacking the chosen bounds simply do not participate,
                // exactly like a read clipped by the interval.
                leftCandidates.clear();
                rightCandidates.clear();
                positionsByMember.assign(members.size(), {});
                bool viable = true;
                for(uint64_t i = 0; i < members.size(); i++) {
                    const auto& [orientedReadId, position] = members[i];
                    const uint64_t v = orientedReadId.getValue();
                    if(v >= orientedReadCount) { viable = false; break; }
                    const auto sorted = anchorIndex[v];

                    const uint32_t low = (position > flankBases) ? (position - flankBases) : 0u;
                    auto it = std::lower_bound(sorted.begin(), sorted.end(),
                        std::make_pair(low, Shasta2AnchorId(0)));
                    for(; it != sorted.end() && it->first <= position + flankBases; ++it) {
                        const uint32_t anchorPosition = it->first;
                        const Shasta2AnchorId id = it->second;
                        // Record the position on THIS read once, so extraction
                        // below is a hash lookup instead of a linear rescan.
                        positionsByMember[i].emplace(id, anchorPosition);
                        if(anchorPosition < position) {
                            Candidate& c = leftCandidates[id];
                            ++c.votes;
                            c.distanceSum += position - anchorPosition;
                        } else if(anchorPosition > position) {
                            Candidate& c = rightCandidates[id];
                            ++c.votes;
                            c.distanceSum += anchorPosition - position;
                        }
                    }
                }
                if(!viable || leftCandidates.empty() || rightCandidates.empty()) {
                    keep[siteIndex] = 0;
                    ++reasons[uint64_t(MsaSiteVerdict::Reason::noSharedAnchors)];
                    continue;
                }

                const Shasta2AnchorId leftAnchor = pickBoundingAnchor(leftCandidates);
                const Shasta2AnchorId rightAnchor = pickBoundingAnchor(rightCandidates);

                // Extract each member's substring between those anchors, where
                // its SNP base sits inside it, and which arm it belongs to.
                sequences.clear();
                snpOffset.clear();
                armOfRow.clear();
                for(uint64_t i = 0; i < members.size(); i++) {
                    const auto& [orientedReadId, position] = members[i];
                    const auto& byId = positionsByMember[i];
                    const auto leftIt = byId.find(leftAnchor);
                    const auto rightIt = byId.find(rightAnchor);
                    if(leftIt == byId.end() || rightIt == byId.end()) continue;
                    const uint32_t leftPos = leftIt->second;
                    const uint32_t rightPos = rightIt->second;
                    if(!(leftPos <= position && position < rightPos)) continue;

                    vector<uint8_t> seq;
                    seq.reserve(rightPos - leftPos);
                    for(uint32_t p = leftPos; p < rightPos; p++) {
                        seq.push_back(uint8_t(
                            reads.getOrientedReadBase(orientedReadId, p).value));
                    }
                    snpOffset.push_back(position - leftPos);
                    armOfRow.push_back(armOfMember[i]);
                    sequences.push_back(std::move(seq));
                }
                if(sequences.size() < defaultMinPlacedRows) {
                    keep[siteIndex] = 0;
                    ++reasons[uint64_t(MsaSiteVerdict::Reason::tooFewMembers)];
                    continue;
                }

                MsaSiteVerdict verdict;
                const bool ran = msaRunner.run(sequences,
                    [&](uint8_t** msa, int nSeq, int msaLen, int gapValue) {
                        verdict = verifyColumnAgreement(
                            msa, nSeq, msaLen, gapValue, snpOffset, armOfRow,
                            minAgreementFraction, defaultMinPlacedRows);
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

    const auto tVerified = std::chrono::steady_clock::now();
    const auto seconds = [](auto a, auto b) {
        return 1.e-9 * double(std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count());
    };
    performanceLog << timestamp << "MSA het verification: "
        << seconds(tBegin, tIndexed) << " s building the anchor index, "
        << seconds(tIndexed, tVerified) << " s verifying "
        << sites.size() << " sites." << endl;

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
