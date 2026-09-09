// Tests for the pure half of the experimental MSA het-site verification
// (src/HetSiteMsaVerification.cpp).
//
// verifyColumnAgreement answers one question: mapped into a single MSA frame,
// do all the member reads' SNP bases land in the SAME column? That is the
// homology check the pairwise CIGAR pass cannot make, so the coordinate walk
// -- offsets count NON-GAP entries, columns include gaps -- is the part worth
// pinning down. These build the matrix by hand; no abPOA, no Assembler.

#include "../external/catch2/catch.hpp"

#include "../src/HetSiteMsaVerification.hpp"

#include <vector>

using namespace dinara;

namespace {

constexpr uint8_t A = 0, C = 1, G = 2, T = 3;
constexpr uint8_t gap = 5;   // abPOA's gap value for a 5-letter alphabet.

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

} // namespace



TEST_CASE("MSA verification confirms a site whose members share one column",
          "[msahet]")
{
    // Four ungapped rows, SNP at offset 2 in every one. Two carry G, two T --
    // a clean biallelic column.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, T, T, A},
        {A, C, T, T, A},
    });
    const std::vector<uint32_t> offsets{2, 2, 2, 2};

    const auto v = verifyColumnAgreement(m.data(), m.nSeq(), m.len(), gap, offsets, 0.9);

    CHECK(v.verified);
    CHECK(v.reason == MsaSiteVerdict::Reason::ok);
    CHECK(v.column == 2);
    CHECK(v.agreeing == 4);
    CHECK(v.total == 4);
    CHECK(v.allelesAtColumn == 2);
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
    const std::vector<uint32_t> offsets{2, 2, 2, 2};

    const auto v = verifyColumnAgreement(m.data(), m.nSeq(), m.len(), gap, offsets, 0.9);

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

    CHECK(verifyColumnAgreement(m.data(), m.nSeq(), m.len(), gap, offsets, 0.75).verified);
    CHECK_FALSE(verifyColumnAgreement(m.data(), m.nSeq(), m.len(), gap, offsets, 0.9).verified);
}


TEST_CASE("MSA verification rejects a column that is not biallelic", "[msahet]")
{
    // Every member agrees on the column AND on the base. Under the MSA this is
    // not a het site at all -- the pairwise disagreement was misplacement.
    Matrix m({
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, G, T, A},
        {A, C, G, T, A},
    });
    const std::vector<uint32_t> offsets{2, 2, 2, 2};

    const auto v = verifyColumnAgreement(m.data(), m.nSeq(), m.len(), gap, offsets, 0.9);

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
    const std::vector<uint32_t> offsets{2, 2, 2, 2};

    const auto v = verifyColumnAgreement(m.data(), m.nSeq(), m.len(), gap, offsets, 0.9);

    CHECK_FALSE(v.verified);
    CHECK(v.reason == MsaSiteVerdict::Reason::notBiallelic);
    CHECK(v.allelesAtColumn == 1);
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
    const std::vector<uint32_t> offsets{2, 2, 2, 2};   // third base = G/T

    const auto v = verifyColumnAgreement(m.data(), m.nSeq(), m.len(), gap, offsets, 0.9);

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
    const std::vector<uint32_t> offsets{2, 2, 2, 2, 2};

    const auto v = verifyColumnAgreement(m.data(), m.nSeq(), m.len(), gap, offsets, 0.9);

    CHECK(v.total == 4);         // the short row did not place
    CHECK(v.agreeing == 4);
    CHECK(v.verified);
}


TEST_CASE("MSA verification refuses fewer than two rows", "[msahet]")
{
    Matrix m({{A, C, G, T, A}});
    const std::vector<uint32_t> offsets{2};

    const auto v = verifyColumnAgreement(m.data(), m.nSeq(), m.len(), gap, offsets, 0.9);

    CHECK_FALSE(v.verified);
    CHECK(v.reason == MsaSiteVerdict::Reason::tooFewMembers);
}
