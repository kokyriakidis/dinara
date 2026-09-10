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
    case MsaSiteVerdict::Reason::noisyNeighbourhood: return "neighbouring columns also disagree (smear, not a point variant)";
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
    uint32_t minPlacedRows,
    uint32_t flankColumns,
    uint32_t noisyFlankColumns)
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

    // INVARIANT CHECK, not a filter. Once every member contributes a window
    // and all of them place in one column, both arms are present by
    // construction and the CIGAR pass already proved their bases differ -- so
    // this cannot fail unless the FRAME is wrong. It fires zero times as
    // shipped, and that zero is the signal: under the shared-anchor windowing
    // this replaced, it fired 4157 times because the frame was starving one arm
    // of the site. Keep it as the canary for the next framing regression.
    //
    // Only rows that PLACED vote here:
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

    // INVARIANT CHECK, not a filter, for the same reason as the biallelic test
    // above -- the arms were assigned from the same bases this reads back, so
    // it fires zero times as shipped. Worth keeping anyway: it is the ONLY
    // guard anywhere on the read PARTITION, which is what phasing consumes and
    // is not measured by any truth comparison in this tree. If the arms and the
    // alignment ever disagree about which read carries which allele, nothing
    // else notices.
    //
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
    // No separate size check here: partitionTotal counts rows that placed AND
    // carry one of the two alleles, and the biallelic gate above already
    // required two alleles with >= 2 rows each, so it cannot be below 4 once
    // that gate has passed.
    if(double(verdict.partitionAgreeing) <
        minAgreementFraction * double(verdict.partitionTotal)) {
        verdict.reason = MsaSiteVerdict::Reason::partitionMismatch;
        return verdict;
    }

    // Does the site LOOK like a point variant, or like a smear?
    //
    // A misaligned difference spreads disagreement across several neighbouring
    // columns; a clean SNP is one disagreeing column among agreeing ones. This
    // is hifiasm's "a snp very close to another is not a real snp"
    // (generate_haplotypes_DP) asked in alignment space, where an intervening
    // indel cannot distort the distance between two sites.
    //
    // It buys precision and costs a LOT of recall, because the premise is only
    // half true for a diploid: a het SNP's neighbours are frequently OTHER het
    // SNPs, since heterozygous variants cluster in haplotype-divergent regions.
    // So it also finds the regions richest in real variants. Measured against
    // the HG002 chr12 track (see the header for the table), rejecting at 4, 2
    // or 1 noisy flanking columns costs 8, 17 and 33 true variants per false
    // positive removed. Off by default for that reason; it is exposed because
    // the right point on that curve is a policy choice, not a fact.
    if(noisyFlankColumns > 0) {
        uint32_t noisy = 0;
        const int lo = std::max(0, bestColumn - int(flankColumns));
        const int hi = std::min(msaLen - 1, bestColumn + int(flankColumns));
        for(int c = lo; c <= hi; c++) {
            if(c == bestColumn) {
                continue;
            }
            std::array<uint32_t, 4> n{0, 0, 0, 0};
            for(int row = 0; row < nSeq; row++) {
                if(columnOfRow[size_t(row)] < 0) continue;
                const uint8_t v = msa[row][c];
                if(int(v) == gapValue || v >= 4) continue;
                ++n[v];
            }
            uint32_t alleles = 0;
            for(const uint32_t k : n) {
                if(k >= 2) ++alleles;
            }
            if(alleles >= 2) ++noisy;
        }
        verdict.noisyFlankColumns = noisy;
        if(noisy >= noisyFlankColumns) {
            verdict.reason = MsaSiteVerdict::Reason::noisyNeighbourhood;
            return verdict;
        }
    }

    verdict.verified = true;
    verdict.reason = MsaSiteVerdict::Reason::ok;
    return verdict;
}



// ============================================================================
// The driver.
// ============================================================================

namespace {

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
        // Disable adaptive banding. A band saves time by assuming the
        // alignment stays near the diagonal, which is the assumption that
        // breaks in the noisy, indel-rich stretches these sites live in.
        // longcallD sets the same for its realignment.
        abpt->wb = -1;
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
    uint32_t flankBases,
    double minAgreementFraction,
    uint32_t flankColumns,
    uint32_t noisyFlankColumns,
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

    vector<uint8_t> keep(sites.size(), 1);
    const uint64_t reasonCount = MsaSiteVerdict::reasonCount;
    vector<vector<uint64_t>> reasonsByThread(
        threadCount, vector<uint64_t>(reasonCount, 0));

    // DYNAMIC scheduling, one site at a time.
    //
    // Cost per site scales with its member count and window length, coverage
    // varies along the genome, and collapseDuplicateLoci leaves sites ordered
    // by read -- so a static partition hands contiguous, cost-correlated blocks
    // to each thread. Measured, that matters less than it sounds: 8.79 s vs
    // 9.54 s on four threads, and nothing at twenty (3.82 vs 3.80). Kept
    // because it is strictly more robust to a skewed site distribution and one
    // relaxed atomic per site is nothing against ~15 ms of POA.
    //
    // Do not read the thread scaling here as poor. This box has TEN physical
    // cores with two threads each, so 4 -> 20 threads is 2.5x more compute, not
    // 5x; the measured 2.3-2.5x is the hardware limit, not an inefficiency.
    std::atomic<uint64_t> nextSite{0};
    vector<std::thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t = 0; t < threadCount; t++) {
        threads.emplace_back([&, t]() {
            vector<uint64_t>& reasons = reasonsByThread[t];

            // One abPOA instance for this thread's whole batch, plus hoisted
            // scratch: this loop runs tens of thousands of times.
            MsaRunner msaRunner;
            vector<std::pair<OrientedReadId, uint32_t>> members;
            vector<uint8_t> armOfMember;
            vector<vector<uint8_t>> sequences;
            vector<uint32_t> snpOffset;
            vector<uint8_t> armOfRow;

            for(;;) {
                const uint64_t siteIndex =
                    nextSite.fetch_add(1, std::memory_order_relaxed);
                if(siteIndex >= sites.size()) {
                    break;
                }
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

                // PER-READ windows: each member contributes flankBases either
                // side of its OWN copy of the SNP. No shared bound, and no
                // anchor involved in the windowing at all.
                //
                // Requiring a shared bounding anchor was self-defeating. An
                // anchor is a k=50 marker stored at its MIDPOINT, so it spans
                // about +/-25 bases, and two reads share it only if they share
                // all 50 -- which reads carrying different alleles cannot do
                // when the window covers the SNP. Every anchor within ~k/2 of
                // the site is therefore arm-specific by construction. Bounding
                // by the nearest shared anchor picks one only a single
                // haplotype holds and starves the MSA of the other arm: that
                // variant rejected 4157 of 4542 sites as "not biallelic", which
                // was the frame's failure being reported as the data's.
                // Bounding by a well-attended anchor avoided that only because
                // attendance is a proxy for "outside the k/2 window", which is
                // why its results drifted with the search radius (1406 sites
                // confirmed at radius 50 against 3998 at 400).
                //
                // POA does not need shared endpoints -- it aligns ragged ends
                // perfectly well -- so dropping the requirement dissolves the
                // conflict: every member participates and no site is lost for
                // want of a bound two haplotypes happen to share.
                sequences.clear();
                snpOffset.clear();
                armOfRow.clear();
                bool viable = true;
                for(uint64_t i = 0; i < members.size(); i++) {
                    const auto& [orientedReadId, position] = members[i];
                    const uint64_t v = orientedReadId.getValue();
                    if(v >= orientedReadCount) { viable = false; break; }

                    const uint64_t readLength =
                        reads.getReadRawSequenceLength(orientedReadId.getReadId());
                    const uint32_t leftPos =
                        (position > flankBases) ? (position - flankBases) : 0u;
                    const uint32_t rightPos = uint32_t(std::min<uint64_t>(
                        readLength, uint64_t(position) + flankBases + 1));
                    if(rightPos <= position) continue;

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
                if(not viable) {
                    keep[siteIndex] = 0;
                    ++reasons[uint64_t(MsaSiteVerdict::Reason::noSharedAnchors)];
                    continue;
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
                            minAgreementFraction, defaultMinPlacedRows,
                            flankColumns, noisyFlankColumns);
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
        << seconds(tBegin, tVerified) << " s for "
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
