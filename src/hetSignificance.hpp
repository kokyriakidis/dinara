#ifndef DINARA_HET_SIGNIFICANCE_HPP
#define DINARA_HET_SIGNIFICANCE_HPP

// Shared significance test for "is this minority allele real, or is it just
// sequencing error misreading the dominant allele?".
//
// Used by both het-site detectors -- the per-edge abPOA one
// (Shasta2AnchorGraphHetOnGraph.cpp) and the CIGAR-driven one
// (AssemblerCigarSnpSites.cpp). It lives here rather than being duplicated
// because the threshold is a decision, not an implementation detail: two
// copies of it would drift and the two detectors would silently disagree
// about what counts as a site.

#include "cstdint.hpp"

#include <boost/math/distributions/binomial.hpp>

namespace dinara {

    // Significance level, matching myloasm's SNPmer caller exactly. A
    // conventional threshold, deliberately not exposed as a tunable.
    constexpr double hetSignificance = 0.05;

    // One-sided binomial tail: with `trials` reads carrying the dominant
    // allele and an assumed per-read error rate, how likely is it that
    // `successes` reads would misread it the SAME way by chance? A value at or
    // below hetSignificance means the count is unlikely to be pure noise -- a
    // real second haplotype, rather than a guess from a fixed count or
    // fraction floor.
    //
    // P(X >= k) = 1 - P(X <= k-1); computed via boost::math's complement idiom
    // rather than 1-cdf(k-1) for numerical stability at small p-values.
    // Returns 1.0 (never significant) for trials == 0 or successes == 0.
    inline double binomialTailPValue(
        uint64_t trials, uint64_t successes, double errorRate)
    {
        if(trials == 0 || successes == 0) return 1.0;
        const boost::math::binomial_distribution<double> dist(
            double(trials), errorRate);
        return boost::math::cdf(
            boost::math::complement(dist, double(successes) - 1.0));
    }

}

#endif
