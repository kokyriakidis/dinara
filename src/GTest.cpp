// G-test (log-likelihood ratio test) for tangle matrices.
// Adapted from shasta2's implementation.

#include "GTest.hpp"
#include "DINARA_ASSERT.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace dinara;
using std::vector;



GTest::GTest(
    const vector<vector<double>>& tangleMatrix,
    double epsilon,
    bool onlyConsiderInjective,
    bool onlyConsiderPermutation)
{
    run(tangleMatrix, epsilon, onlyConsiderInjective, onlyConsiderPermutation);
}



GTest::GTest(
    const vector<vector<uint64_t>>& tangleMatrixInteger,
    double epsilon,
    bool onlyConsiderInjective,
    bool onlyConsiderPermutation)
{
    vector<vector<double>> tangleMatrixDouble;
    for(const auto& row : tangleMatrixInteger) {
        vector<double> rowDouble;
        for(const uint64_t v : row) {
            rowDouble.push_back(double(v));
        }
        tangleMatrixDouble.push_back(rowDouble);
    }
    run(tangleMatrixDouble, epsilon, onlyConsiderInjective, onlyConsiderPermutation);
}



void GTest::run(
    const vector<vector<double>>& tangleMatrix,
    double epsilon,
    bool onlyConsiderInjective,
    bool onlyConsiderPermutation)
{
    const uint64_t entranceCount = tangleMatrix.size();
    DINARA_ASSERT(entranceCount > 0);

    const uint64_t exitCount = tangleMatrix.front().size();
    DINARA_ASSERT(exitCount > 0);
    for(uint64_t i = 0; i < entranceCount; i++) {
        DINARA_ASSERT(tangleMatrix[i].size() == exitCount);
    }

    // Limit to 16 entries max (4x4 matrix).
    const uint64_t totalEntries = entranceCount * exitCount;
    if(totalEntries > 16) {
        success = false;
        return;
    }

    // Compute coverage sums per entrance (row) and per exit (column).
    double totalCoverage = 0;
    vector<double> entranceCoverage(entranceCount, 0.);
    vector<double> exitCoverage(exitCount, 0.);
    for(uint64_t i = 0; i < entranceCount; i++) {
        for(uint64_t j = 0; j < exitCount; j++) {
            const double v = tangleMatrix[i][j];
            totalCoverage += v;
            entranceCoverage[i] += v;
            exitCoverage[j] += v;
        }
    }

    // Random (uniform) tangle matrix: expected under no connectivity structure.
    vector<vector<double>> randomMatrix(entranceCount, vector<double>(exitCount, 0.));
    for(uint64_t i = 0; i < entranceCount; i++) {
        for(uint64_t j = 0; j < exitCount; j++) {
            randomMatrix[i][j] =
                entranceCoverage[i] * exitCoverage[j] / totalCoverage;
        }
    }

    // Temporary matrices.
    vector<vector<double>> idealA(entranceCount, vector<double>(exitCount));
    vector<vector<double>> idealB(entranceCount, vector<double>(exitCount));
    vector<vector<double>> idealMatrix(entranceCount, vector<double>(exitCount));
    vector<vector<double>> expectedMatrix(entranceCount, vector<double>(exitCount));
    vector<vector<bool>> connectivityMatrix(entranceCount, vector<bool>(exitCount));

    // Enumerate all 2^(entranceCount*exitCount) connectivity hypotheses.
    const uint64_t N = 1ULL << totalEntries;
    hypotheses.clear();

    for(uint64_t connInt = 0; connInt < N; connInt++) {

        // Build connectivity matrix from bits.
        uint64_t mask = 1;
        for(uint64_t i = 0; i < entranceCount; i++) {
            for(uint64_t j = 0; j < exitCount; j++) {
                connectivityMatrix[i][j] = ((connInt & mask) != 0);
                mask <<= 1;
            }
        }

        // Filter by injectivity/permutation constraints.
        const bool fwdInj = isForwardInjective(connectivityMatrix);
        const bool bwdInj = isBackwardInjective(connectivityMatrix);
        if(onlyConsiderInjective && !(fwdInj || bwdInj)) continue;
        if(onlyConsiderPermutation && !(fwdInj && bwdInj)) continue;

        // Compute ideal tangle matrix A: distribute each entrance's
        // coverage equally among its connected exits.
        for(uint64_t i = 0; i < entranceCount; i++) {
            uint64_t nonZero = 0;
            for(uint64_t j = 0; j < exitCount; j++) {
                if(connectivityMatrix[i][j]) ++nonZero;
            }
            if(nonZero == 0) {
                for(uint64_t j = 0; j < exitCount; j++) idealA[i][j] = 0.;
            } else {
                const double val = entranceCoverage[i] / double(nonZero);
                for(uint64_t j = 0; j < exitCount; j++) {
                    idealA[i][j] = connectivityMatrix[i][j] ? val : 0.;
                }
            }
        }

        // Compute ideal tangle matrix B: distribute each exit's
        // coverage equally among its connected entrances.
        for(uint64_t j = 0; j < exitCount; j++) {
            uint64_t nonZero = 0;
            for(uint64_t i = 0; i < entranceCount; i++) {
                if(connectivityMatrix[i][j]) ++nonZero;
            }
            if(nonZero == 0) {
                for(uint64_t i = 0; i < entranceCount; i++) idealB[i][j] = 0.;
            } else {
                const double val = exitCoverage[j] / double(nonZero);
                for(uint64_t i = 0; i < entranceCount; i++) {
                    idealB[i][j] = connectivityMatrix[i][j] ? val : 0.;
                }
            }
        }

        // Average the two ideal matrices.
        for(uint64_t i = 0; i < entranceCount; i++) {
            for(uint64_t j = 0; j < exitCount; j++) {
                idealMatrix[i][j] = 0.5 * (idealA[i][j] + idealB[i][j]);
            }
        }

        // Expected = (1-epsilon)*ideal + epsilon*random.
        for(uint64_t i = 0; i < entranceCount; i++) {
            for(uint64_t j = 0; j < exitCount; j++) {
                expectedMatrix[i][j] =
                    (1. - epsilon) * idealMatrix[i][j] +
                    epsilon * randomMatrix[i][j];
            }
        }

        // G statistic: 2 * sum(observed * ln(observed/expected)).
        // Converted to decibels (multiply by 10, use log10).
        double G = 0.;
        for(uint64_t i = 0; i < entranceCount; i++) {
            for(uint64_t j = 0; j < exitCount; j++) {
                const double observed = tangleMatrix[i][j];
                if(observed > 0.) {
                    const double expected = expectedMatrix[i][j];
                    G += observed * log10(observed / expected);
                }
            }
        }
        G *= 20.; // Factor of 2 and convert to decibels.

        hypotheses.emplace_back(connectivityMatrix, G);
    }

    std::sort(hypotheses.begin(), hypotheses.end());
    success = !hypotheses.empty();
}



bool GTest::isForwardInjective(const vector<vector<bool>>& m)
{
    for(uint64_t i = 0; i < m.size(); i++) {
        uint64_t count = 0;
        for(uint64_t j = 0; j < m[i].size(); j++) {
            if(m[i][j]) ++count;
        }
        if(count != 1) return false;
    }
    return true;
}



bool GTest::isBackwardInjective(const vector<vector<bool>>& m)
{
    if(m.empty()) return false;
    const uint64_t exitCount = m.front().size();
    for(uint64_t j = 0; j < exitCount; j++) {
        uint64_t count = 0;
        for(uint64_t i = 0; i < m.size(); i++) {
            if(m[i][j]) ++count;
        }
        if(count != 1) return false;
    }
    return true;
}
