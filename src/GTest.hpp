#pragma once

// G-test (log-likelihood ratio test) for tangle matrices.
// Adapted from shasta2's implementation.
// https://en.wikipedia.org/wiki/G-test
//
// Given an observed tangle matrix (entrances x exits), enumerates all
// possible connectivity hypotheses and computes the G statistic for
// each. The hypothesis with the lowest G is the best fit.

#include <cstdint>
#include <iosfwd>
#include <vector>

namespace dinara {
    class GTest;
}

class dinara::GTest {
public:
    // Construct and run the G-test on the given tangle matrix.
    // epsilon: mixing weight for the random (uniform) component.
    //   Expected = (1-epsilon)*ideal + epsilon*random.
    // onlyConsiderInjective: skip hypotheses where an entrance maps
    //   to multiple exits or vice versa.
    // onlyConsiderPermutation: only consider 1-to-1 mappings.
    GTest(
        const std::vector<std::vector<uint64_t>>& tangleMatrix,
        double epsilon,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation);

    GTest(
        const std::vector<std::vector<double>>& tangleMatrix,
        double epsilon,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation);

    bool success = false;

    class Hypothesis {
    public:
        std::vector<std::vector<bool>> connectivityMatrix;
        double G;

        Hypothesis(
            const std::vector<std::vector<bool>>& connectivityMatrix,
            double G) :
            connectivityMatrix(connectivityMatrix),
            G(G)
            {}

        // Sort by G (lowest = best fit).
        bool operator<(const Hypothesis& that) const {
            return G < that.G;
        }
    };
    std::vector<Hypothesis> hypotheses;

    // Check if each entrance maps to exactly one exit.
    static bool isForwardInjective(const std::vector<std::vector<bool>>& connectivityMatrix);

    // Check if each exit maps to exactly one entrance.
    static bool isBackwardInjective(const std::vector<std::vector<bool>>& connectivityMatrix);

private:
    void run(
        const std::vector<std::vector<double>>& tangleMatrix,
        double epsilon,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation);
};
