#ifndef DINARA_HET_SITE_MSA_VERIFICATION_HPP
#define DINARA_HET_SITE_MSA_VERIFICATION_HPP

// EXPERIMENTAL: verify CIGAR-detected het sites with an MSA over the interval
// between two anchors the member reads share.
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
// An MSA puts every covering read in ONE coordinate frame, so three questions
// become answerable, in increasing order of what they buy:
//
//   1. Do the members' SNP bases land in the same COLUMN?   (homology)
//   2. Is that column actually BIALLELIC?                   (a real variant)
//   3. Does the column's read split match the ARM split
//      the CIGAR pass assigned?                             (usable for phasing)
//
// (3) is the one that matters downstream. Two arms can both be biallelic and
// still disagree about WHICH reads carry which allele, and a site like that is
// worse than no site at all: phasing takes the partition at face value. A site
// that passes (1) and (2) but fails (3) is a confidently wrong haplotype call.
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
    // Rows whose (arm, base) matches the best arm->allele bijection, over the
    // rows that placed and carry one of the two alleles.
    uint32_t partitionAgreeing = 0;
    uint32_t partitionTotal = 0;

    // Why a site failed, for the funnel report. Stable and ordered so the
    // report reads as a funnel; keep in sync with reasonName().
    enum class Reason {
        ok,
        noSharedAnchors,    // no bracketing anchor carried by enough members
        tooFewMembers,      // too few members survived interval extraction
        msaFailed,          // abPOA produced no usable matrix
        columnsDisagree,    // members' SNP bases landed in different columns
        notBiallelic,       // agreed column is not cleanly biallelic
        partitionMismatch   // column splits the reads differently than the arms
    };
    static constexpr uint64_t reasonCount = 7;
    Reason reason = Reason::noSharedAnchors;
};

// Human-readable name for a reason, for the funnel report.
const char* msaVerdictReasonName(MsaSiteVerdict::Reason);


// The pure half: given an MSA matrix, the offset of each row's SNP base within
// its own sequence, and the arm the CIGAR pass assigned that row, decide
// whether the site survives.
//
// `msa` is row-major with `nSeq` rows of `msaLen` entries, encoded as abPOA
// encodes them (0=A 1=C 2=G 3=T, >=4 non-informative, gapValue = gap). Rows are
// members only -- there is no backbone row, because the interval is bounded by
// shared anchors rather than built around a spine.
//
// `snpOffsetInRow[row]` counts NON-GAP entries: the SNP base is the
// snpOffsetInRow[row]-th base of that row's sequence. `armOfRow[row]` is 0 or 1.
//
// A site survives when, in order:
//   - at least `minPlacedRows` rows place at all (below that the biallelic and
//     partition tests have too little to work with to mean anything);
//   - at least `minAgreementFraction` of placed rows land on one column;
//   - that column carries two distinct bases with >= 2 supporting rows;
//   - at least `minAgreementFraction` of the rows carrying one of those two
//     bases sit on the arm the best bijection assigns.
//
// Kept free of abPOA, Reads and Assembler so it can be tested on hand-built
// matrices; the coordinate logic is where this gets subtly wrong.
MsaSiteVerdict verifyColumnAgreement(
    const uint8_t* const* msa,
    int nSeq,
    int msaLen,
    int gapValue,
    const vector<uint32_t>& snpOffsetInRow,
    const vector<uint8_t>& armOfRow,
    double minAgreementFraction,
    uint32_t minPlacedRows);


// The driver: verify every site and erase the ones that fail, returning the
// number erased. Sites are independent, so this runs on `threadCount` threads.
//
// `flankBases` is the half-width of the window each member contributes, centred
// on its own copy of the site. No anchors are involved: see the note above.
uint64_t msaVerifyHetSites(
    vector<Assembler::CigarSnpSite>& sites,
    const Reads& reads,
    uint32_t flankBases,
    double minAgreementFraction,
    uint64_t threadCount,
    vector<uint64_t>* reasonCountsOut = nullptr);

}

#endif
