#include "Shasta2GTest.hpp"

#include "DINARA_ASSERT.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>

using namespace dinara;
using namespace std;

Shasta2GTest::Shasta2GTest(
    const vector< vector<double> >& tangleMatrix,
    const double epsilon,
    const bool onlyConsiderInjective,
    const bool onlyConsiderPermutation)
{
    run(tangleMatrix, epsilon, onlyConsiderInjective, onlyConsiderPermutation);
}

Shasta2GTest::Shasta2GTest(
    const vector< vector<uint64_t> >& tangleMatrixInteger,
    const double epsilon,
    const bool onlyConsiderInjective,
    const bool onlyConsiderPermutation)
{
    vector< vector<double> > tangleMatrixDouble;
    for(const vector<uint64_t>& rowInteger: tangleMatrixInteger) {
        vector<double> rowDouble;
        for(const uint64_t valueInteger: rowInteger) {
            rowDouble.push_back(double(valueInteger));
        }
        tangleMatrixDouble.push_back(rowDouble);
    }

    run(tangleMatrixDouble, epsilon, onlyConsiderInjective, onlyConsiderPermutation);
}

void Shasta2GTest::run(
    const vector< vector<double> >& tangleMatrix,
    const double epsilon,
    const bool onlyConsiderInjective,
    const bool onlyConsiderPermutation)
{
    const uint64_t entranceCount = tangleMatrix.size();
    DINARA_ASSERT(entranceCount > 0);

    const uint64_t exitCount = tangleMatrix.front().size();
    DINARA_ASSERT(exitCount > 0);
    for(uint64_t i=0; i<entranceCount; i++) {
        DINARA_ASSERT(tangleMatrix[i].size() == exitCount);
    }

    const uint64_t totalTangleMatrixEntryCount = entranceCount * exitCount;
    if(totalTangleMatrixEntryCount > 16) {
        success = false;
        return;
    }

    double totalCommonCoverage = 0.;
    vector<double> entranceCommonCoverage(entranceCount, 0.);
    vector<double> exitCommonCoverage(exitCount, 0.);
    for(uint64_t i=0; i<entranceCount; i++) {
        for(uint64_t j=0; j<exitCount; j++) {
            const double coverage = tangleMatrix[i][j];
            totalCommonCoverage += coverage;
            entranceCommonCoverage[i] += coverage;
            exitCommonCoverage[j] += coverage;
        }
    }

    vector< vector<double> > randomTangleMatrix(entranceCount, vector<double>(exitCount, 0.));
    for(uint64_t i=0; i<entranceCount; i++) {
        for(uint64_t j=0; j<exitCount; j++) {
            randomTangleMatrix[i][j] =
                entranceCommonCoverage[i] *
                exitCommonCoverage[j] /
                totalCommonCoverage;
        }
    }

    vector< vector<double> > idealTangleMatrixA(entranceCount, vector<double>(exitCount));
    vector< vector<double> > idealTangleMatrixB(entranceCount, vector<double>(exitCount));
    vector< vector<double> > idealTangleMatrix(entranceCount, vector<double>(exitCount));
    vector< vector<double> > expectedTangleMatrix(entranceCount, vector<double>(exitCount));
    vector< vector<bool> > connectivityMatrix(entranceCount, vector<bool>(exitCount));

    const uint64_t N = 1ULL << totalTangleMatrixEntryCount;
    hypotheses.clear();
    for(uint64_t connectivityInteger=0; connectivityInteger<N; connectivityInteger++) {

        uint64_t mask = 1;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                connectivityMatrix[iEntrance][iExit] = ((connectivityInteger & mask) != 0);
                mask = mask << 1;
            }
        }

        const bool hypothesisIsForwardInjective = isForwardInjective(connectivityMatrix);
        const bool hypothesisIsBackwardInjective = isBackwardInjective(connectivityMatrix);
        const bool isInjective = hypothesisIsForwardInjective || hypothesisIsBackwardInjective;
        const bool isPermutation = hypothesisIsForwardInjective && hypothesisIsBackwardInjective;
        if(onlyConsiderInjective && (!isInjective)) {
            continue;
        }
        if(onlyConsiderPermutation && (!isPermutation)) {
            continue;
        }

        for(uint64_t i=0; i<entranceCount; i++) {
            uint64_t nonZeroCount = 0;
            for(uint64_t j=0; j<exitCount; j++) {
                if(connectivityMatrix[i][j]) {
                    ++nonZeroCount;
                }
            }
            if(nonZeroCount == 0) {
                for(uint64_t j=0; j<exitCount; j++) {
                    idealTangleMatrixA[i][j] = 0.;
                }
            } else {
                const double value = entranceCommonCoverage[i] / double(nonZeroCount);
                for(uint64_t j=0; j<exitCount; j++) {
                    if(connectivityMatrix[i][j]) {
                        idealTangleMatrixA[i][j] = value;
                    } else {
                        idealTangleMatrixA[i][j] = 0.;
                    }
                }
            }
        }

        for(uint64_t j=0; j<exitCount; j++) {
            uint64_t nonZeroCount = 0;
            for(uint64_t i=0; i<entranceCount; i++) {
                if(connectivityMatrix[i][j]) {
                    ++nonZeroCount;
                }
            }
            if(nonZeroCount == 0) {
                for(uint64_t i=0; i<entranceCount; i++) {
                    idealTangleMatrixB[i][j] = 0.;
                }
            } else {
                const double value = exitCommonCoverage[j] / double(nonZeroCount);
                for(uint64_t i=0; i<entranceCount; i++) {
                    if(connectivityMatrix[i][j]) {
                        idealTangleMatrixB[i][j] = value;
                    } else {
                        idealTangleMatrixB[i][j] = 0.;
                    }
                }
            }
        }

        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            for(uint64_t iExit=0; iExit<exitCount; iExit++) {
                idealTangleMatrix[iEntrance][iExit] =
                    0.5 * (idealTangleMatrixA[iEntrance][iExit] + idealTangleMatrixB[iEntrance][iExit]);
            }
        }

        for(uint64_t i=0; i<entranceCount; i++) {
            for(uint64_t j=0; j<exitCount; j++) {
                expectedTangleMatrix[i][j] =
                    (1. - epsilon) * idealTangleMatrix[i][j] +
                    epsilon * randomTangleMatrix[i][j];
            }
        }

        double G = 0.;
        for(uint64_t i=0; i<entranceCount; i++) {
            for(uint64_t j=0; j<exitCount; j++) {
                const double actualCoverage = tangleMatrix[i][j];
                if(actualCoverage > 0.) {
                    const double expectedCoverage = expectedTangleMatrix[i][j];
                    G += actualCoverage * log10(actualCoverage / expectedCoverage);
                }
            }
        }
        G *= 20.;

        hypotheses.emplace_back(Hypothesis(connectivityMatrix, G));
    }

    sort(hypotheses.begin(), hypotheses.end());
    success = !hypotheses.empty();
}

void Shasta2GTest::writeHtml(ostream& html) const
{
    html <<
        "<h3>G test</h3>"
        "<table><tr>"
        "<th>Connectivity<br>matrix"
        "<th>G<br>(dB)";

    for(const auto& hypothesis: hypotheses) {
        const vector< vector<bool> >& connectivityMatrix = hypothesis.connectivityMatrix;

        html << "<tr><td style='display: flex; align-items: center; justify-content: center;'>";

        html << "<table>";
        for(uint64_t i=0; i<connectivityMatrix.size(); i++) {
            html << "<tr>";
            for(uint64_t j=0; j<connectivityMatrix[i].size(); j++) {
                html << "<td class=centered>" << hypothesis.connectivityMatrix[i][j];
            }
        }
        html << "</table>";

        html << "<td class=centered>" << fixed << setprecision(1) << hypothesis.G;
    }

    html << "</table>";
}

bool Shasta2GTest::isForwardInjective(const vector< vector<bool> >& connectivityMatrix)
{
    const uint64_t entranceCount = connectivityMatrix.size();
    const uint64_t exitCount = connectivityMatrix.front().size();

    for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
        uint64_t count = 0;
        for(uint64_t iExit=0; iExit<exitCount; iExit++) {
            if(connectivityMatrix[iEntrance][iExit]) {
                ++count;
            }
        }
        if(count != 1) {
            return false;
        }
    }
    return true;
}

bool Shasta2GTest::isBackwardInjective(const vector< vector<bool> >& connectivityMatrix)
{
    const uint64_t entranceCount = connectivityMatrix.size();
    const uint64_t exitCount = connectivityMatrix.front().size();

    for(uint64_t iExit=0; iExit<exitCount; iExit++) {
        uint64_t count = 0;
        for(uint64_t iEntrance=0; iEntrance<entranceCount; iEntrance++) {
            if(connectivityMatrix[iEntrance][iExit]) {
                ++count;
            }
        }
        if(count != 1) {
            return false;
        }
    }
    return true;
}
