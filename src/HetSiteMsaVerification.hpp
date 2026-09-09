#ifndef DINARA_HET_SITE_MSA_VERIFICATION_HPP
#define DINARA_HET_SITE_MSA_VERIFICATION_HPP

// EXPERIMENTAL: verify CIGAR-detected het sites with an MSA over the interval
// between two shared anchors.
//
// WHAT THIS CHECKS THAT NOTHING ELSE DOES
//
// detectCigarSnpSites works from PAIRWISE hifiasm CIGARs. Its self-check
// confirms that a mismatch column's two bases really differ -- but not that
// they are HOMOLOGOUS. A pairwise alignment that places a mismatch one base off
// (routine next to an indel or inside a short repeat) yields a column whose
// bases genuinely differ while describing different loci, and that self-check
// passes it. Nor is agreement across pairs independent evidence: alignments
// sharing an indel context tend to share the same misplacement, so one artifact
// can be counted many times.
//
// An MSA puts every covering read in ONE coordinate frame, so the question
// "are these bases the same locus?" becomes answerable: map each member's SNP
// base into the MSA and see whether they all land in the same column.
//
// WHY ANCHORS ARE THE RIGHT INTERVAL
//
// Every member read shares the bounding anchors' k-mers, so the substrings
// between them are homologous by construction and comparable in length -- a
// well-posed MSA. This is what an arbitrary sequence window cannot guarantee,
// and why bounding by shared anchors is the load-bearing part of the design
// rather than an implementation detail.
//
// LIMITS
//
// Inside a short tandem repeat the MSA is itself ambiguous: abPOA picks *an*
// alignment that need not be more correct than the pairwise one. This raises
// confidence where a consistent frame exists; it does not rescue the cases
// where none does. It is a filter, so it can only ever REMOVE sites -- the risk
// is bounded and shows up as a recall cost, which is measurable against a truth
// set.

#include "Assembler.hpp"
#include "Shasta2Anchors.hpp"
#include "cstdint.hpp"
#include "vector.hpp"

namespace dinara {

class Reads;

// Per-site outcome, for reporting and for tests.
class MsaSiteVerdict {
public:
    bool verified = false;
    int column = -1;              // agreed MSA column, or -1
    uint32_t agreeing = 0;        // rows whose SNP base landed on that column
    uint32_t total = 0;           // rows placed in the MSA at all
    uint32_t allelesAtColumn = 0; // distinct bases there with >= 2 support

    // Why a site failed, for the funnel report. Stable strings, not prose.
    enum class Reason {
        ok,
        noSharedAnchors,   // no anchor pair brackets the site in every member
        tooFewMembers,     // fewer than 2 members survived interval extraction
        msaFailed,         // abPOA produced no usable matrix
        columnsDisagree,   // members' SNP bases landed in different columns
        notBiallelic       // agreed column is not cleanly biallelic
    };
    Reason reason = Reason::noSharedAnchors;
};


// The pure half: given an MSA matrix and, per row, the offset of that row's SNP
// base within its own sequence, decide whether every row's SNP base lands in
// the same column.
//
// `msa` is row-major with `nSeq` rows of `msaLen` entries, encoded as abPOA
// encodes them (0=A 1=C 2=G 3=T, >=4 non-informative, gapValue = gap). Rows are
// members only -- there is no backbone row, because the interval is bounded by
// shared anchors rather than built around a spine.
//
// `snpOffsetInRow[row]` counts NON-GAP entries: the SNP base is the
// snpOffsetInRow[row]-th base of that row's sequence.
//
// A site is verified when at least `minAgreementFraction` of placed rows land
// on the majority column AND that column carries at least two distinct bases
// with >= 2 supporting rows (a het site that collapses to one allele under MSA
// was a pairwise artifact).
//
// Kept free of abPOA, Reads and Assembler so it can be tested on hand-built
// matrices; the coordinate logic is where this gets subtly wrong.
MsaSiteVerdict verifyColumnAgreement(
    const uint8_t* const* msa,
    int nSeq,
    int msaLen,
    int gapValue,
    const vector<uint32_t>& snpOffsetInRow,
    double minAgreementFraction);


// The driver: verify every site and erase the ones that fail, returning the
// number erased. Sites are independent, so this runs on `threadCount` threads.
//
// `flankBases` is how far to look on each side of the SNP for a bounding anchor
// shared by all members; a larger value gives the MSA more context and a better
// chance of a shared pair, at the cost of a longer alignment.
uint64_t msaVerifyHetSites(
    vector<Assembler::CigarSnpSite>& sites,
    const Reads& reads,
    const Shasta2Anchors& anchors,
    uint32_t flankBases,
    double minAgreementFraction,
    uint64_t threadCount,
    vector<uint64_t>* reasonCountsOut = nullptr);

}

#endif
