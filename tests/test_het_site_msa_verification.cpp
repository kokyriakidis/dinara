// Tests for the pure half of the experimental MSA het-site verification
// (src/HetSiteMsaVerification.cpp).
//
// verifyColumnAgreement answers three questions the pairwise CIGAR pass cannot:
// do the members' SNP bases land in the same COLUMN (homology), is that column
// BIALLELIC (a real variant), and does its read split match the ARM split the
// CIGAR pass assigned (usable for phasing). The coordinate walk -- offsets
// count NON-GAP entries, columns include gaps -- and the arm bijection are
// where this gets silently wrong, so both are pinned here on hand-built
// matrices. No abPOA, no Assembler.

#include "../external/catch2/catch.hpp"

#include "../src/HetSiteMsaVerification.hpp"

#include <vector>

using namespace dinara;

namespace {

constexpr uint8_t A = 0, C = 1, G = 2, T = 3;
constexpr uint8_t gap = 5;       // abPOA's gap value for a 5-letter alphabet.
constexpr uint32_t minPlaced = 4;

// Build a row-pointer view over hand-written rows, as abPOA hands them over.
class Matrix {
public:
    std::vector<std::vector<uint8_t>> rows;
    std::vector<const uint8_t*> ptrs;

    explicit Matrix(std::vector<std::vector<uint8_t>> r) : rows(std::move(r))
    {
        for(const auto& row: rows) ptrs.push_back(row.data());
    }
    const uint8_t* const* data() const { return ptrs.data(); }
    int nSeq() const { return int(rows.size()); }
    int len() const { return rows.empty() ? 0 : int(rows[0].size()); }
};

MsaSiteVerdict run(const Matrix& m,
    const std::vector<uint32_t>& offsets,
    const std::vector<uint8_t>& arms,
    double minAgreement = 0.9,
    uint32_t rejectAtNoisyFlank = 0)      // 0 = neighbourhood test off, as shipped
{
    return verifyColumnAgreement(
        m.data(), m.nSeq(), m.len(), gap, offsets, arms, minAgreement, minPlaced,
        /*flankColumns*/ 10, rejectAtNoisyFlank);
}

} // namespace



TEST_CASE("MSA verification confirms a site whose members share one column",
          "[msahet]")
{
    // Four ungapped rows, SNP at offset 2 in every one. Two carry G, two T, and
    // the arms agree with the bases -- a clean biallelic column.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, T, T, A},
        {A, C, T, T, A},
    });

    const auto v = run(m, {2, 2, 2, 2}, {0, 0, 1, 1});

    CHECK(v.verified);
    CHECK(v.reason == MsaSiteVerdict::Reason::ok);
    CHECK(v.column == 2);
    CHECK(v.agreeing == 4);
    CHECK(v.total == 4);
    CHECK(v.allelesAtColumn == 2);
    CHECK(v.partitionAgreeing == 4);
    CHECK(v.partitionTotal == 4);
}


TEST_CASE("MSA verification rejects a site whose members land in different columns",
          "[msahet]")
{
    // This is the artifact the pass exists to catch. Rows 2 and 3 carry an
    // insertion column, so their third base sits at column 3 while rows 0 and 1
    // have theirs at column 2. The pairwise CIGARs called this one site; the
    // MSA says the bases are not the same locus.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, gap, C, T, A},
        {A, gap, C, T, A},
    });

    const auto v = run(m, {2, 2, 2, 2}, {0, 0, 1, 1});

    CHECK_FALSE(v.verified);
    CHECK(v.reason == MsaSiteVerdict::Reason::columnsDisagree);
    CHECK(v.agreeing == 2);      // a 2-2 split reaches no majority at 0.9
    CHECK(v.total == 4);
}


TEST_CASE("MSA verification tolerates a minority landing elsewhere", "[msahet]")
{
    // One row of five disagreeing is 80% agreement: kept at a 0.75 threshold,
    // rejected at 0.9. The threshold is the whole knob, so both directions are
    // pinned.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, T, T, A},
        {A, C, T, T, A},
        {A, gap, C, T, A},
    });
    const std::vector<uint32_t> offsets{2, 2, 2, 2, 2};
    const std::vector<uint8_t> arms{0, 0, 1, 1, 1};

    CHECK(run(m, offsets, arms, 0.75).verified);
    CHECK_FALSE(run(m, offsets, arms, 0.9).verified);
}


TEST_CASE("MSA verification rejects a column that is not biallelic", "[msahet]")
{
    // Every member agrees on the column AND on the base. Under the MSA this is
    // not a het site at all -- the pairwise disagreement was misplacement.
    // Measured on E821 this is the largest rejection class, and the one that
    // was corrupting phasing.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, G, T, A},
    });

    const auto v = run(m, {2, 2, 2, 2}, {0, 0, 1, 1});

    CHECK_FALSE(v.verified);
    CHECK(v.reason == MsaSiteVerdict::Reason::notBiallelic);
    CHECK(v.column == 2);
    CHECK(v.allelesAtColumn == 1);
}


TEST_CASE("MSA verification needs two reads supporting each allele", "[msahet]")
{
    // Three G and one T: the T is a singleton, so the column carries only one
    // allele with >= 2 support. A single divergent read is an error, not a
    // haplotype.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, T, T, A},
    });

    const auto v = run(m, {2, 2, 2, 2}, {0, 0, 0, 1});

    CHECK_FALSE(v.verified);
    CHECK(v.reason == MsaSiteVerdict::Reason::notBiallelic);
    CHECK(v.allelesAtColumn == 1);
}


TEST_CASE("MSA verification rejects a column that splits reads unlike the arms",
          "[msahet]")
{
    // Biallelic, one column, every member placed -- and still useless. The MSA
    // says {row0,row1} carry G and {row2,row3} carry T, while the CIGAR pass
    // put row0 and row2 in one arm. Both partitions cannot be right, and
    // phasing consumes the arms at face value, so this site is worse than none.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, T, T, A},
        {A, C, T, T, A},
    });

    const auto v = run(m, {2, 2, 2, 2}, {0, 1, 0, 1});

    CHECK_FALSE(v.verified);
    CHECK(v.reason == MsaSiteVerdict::Reason::partitionMismatch);
    CHECK(v.allelesAtColumn == 2);        // it IS biallelic
    CHECK(v.partitionAgreeing == 2);      // 2 of 4 either way round
    CHECK(v.partitionTotal == 4);
}


TEST_CASE("MSA verification accepts a partition whose arm labels are swapped",
          "[msahet]")
{
    // Arm 0 vs arm 1 is a label, not an identity: the detector's arm order
    // carries no meaning. Only the SPLIT has to match, so the check scores both
    // bijections and takes the better. Rejecting this would throw away half of
    // all correct sites.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, T, T, A},
        {A, C, T, T, A},
    });

    const auto v = run(m, {2, 2, 2, 2}, {1, 1, 0, 0});

    CHECK(v.verified);
    CHECK(v.partitionAgreeing == 4);
    CHECK(v.partitionTotal == 4);
}


TEST_CASE("MSA verification counts offsets in non-gap entries, not columns",
          "[msahet]")
{
    // The offset is an index into the READ's bases; leading gaps must not
    // consume it. Getting this backwards silently shifts every site.
    Matrix m({
        {gap, gap, A, C, G, T},
        {gap, gap, A, C, G, T},
        {gap, gap, A, C, T, T},
        {gap, gap, A, C, T, T},
    });

    const auto v = run(m, {2, 2, 2, 2}, {0, 0, 1, 1});   // third base = G/T

    CHECK(v.verified);
    CHECK(v.column == 4);        // not 2: two leading gap columns
    CHECK(v.allelesAtColumn == 2);
}


TEST_CASE("MSA verification ignores a row whose sequence is too short", "[msahet]")
{
    // A row clipped by the interval never reaches its offset. It must not vote
    // and must not crash -- it simply does not participate. The other four rows
    // are left as a clean 2-2 split so this isolates the non-placement: an
    // earlier version of this test used three rows, where dropping the short
    // one also broke biallelism, and so tested two things at once.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, T, T, A},
        {A, C, T, T, A},
        {A, gap, gap, gap, gap},   // one base only; offset 2 unreachable
    });

    const auto v = run(m, {2, 2, 2, 2, 2}, {0, 0, 1, 1, 1});

    CHECK(v.total == 4);         // the short row did not place
    CHECK(v.agreeing == 4);
    CHECK(v.verified);
}


TEST_CASE("MSA verification refuses too few placed rows", "[msahet]")
{
    // Two alleles need >= 2 rows each, so under four placed rows the biallelic
    // and partition tests can only pass or fail for want of evidence. Three
    // perfectly consistent rows are still not enough to conclude anything.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, T, T, A},
    });

    const auto v = run(m, {2, 2, 2}, {0, 0, 1});

    CHECK_FALSE(v.verified);
    CHECK(v.reason == MsaSiteVerdict::Reason::tooFewMembers);
    CHECK(v.total == 3);
}


TEST_CASE("MSA verification refuses fewer than two rows", "[msahet]")
{
    Matrix m({{A, C, G, T, A}});

    const auto v = run(m, {2}, {0});

    CHECK_FALSE(v.verified);
    CHECK(v.reason == MsaSiteVerdict::Reason::tooFewMembers);
}


TEST_CASE("MSA verification can reject a site whose neighbourhood is also noisy",
          "[msahet]")
{
    // Column 2 is the site and columns 1 and 3 disagree too -- the signature of
    // a smeared misalignment rather than a point variant. Off by default
    // (rejectAtNoisyFlank 0) because a het SNP's neighbours are often other het
    // SNPs; this pins the mechanism, not the policy.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, G, T, C, A},
        {A, G, T, C, A},
    });
    const std::vector<uint32_t> offsets{2, 2, 2, 2};
    const std::vector<uint8_t> arms{0, 0, 1, 1};

    CHECK(run(m, offsets, arms, 0.9, /*rejectAtNoisyFlank*/ 0).verified);

    const auto strict = run(m, offsets, arms, 0.9, /*rejectAtNoisyFlank*/ 2);
    CHECK_FALSE(strict.verified);
    CHECK(strict.reason == MsaSiteVerdict::Reason::noisyNeighbourhood);
    CHECK(strict.noisyFlankColumns == 2);      // columns 1 and 3
}
