#pragma once

#include "invalid.hpp"

#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "vector.hpp"

namespace dinara {
    class Shasta2GTest;
}

class dinara::Shasta2GTest {
public:
    Shasta2GTest(
        const vector< vector<uint64_t> >& tangleMatrix,
        double epsilon,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation);

    Shasta2GTest(
        const vector< vector<double> >& tangleMatrix,
        double epsilon,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation);

    bool success = false;

    void writeHtml(ostream&) const;

    class Hypothesis {
    public:
        vector< vector<bool> > connectivityMatrix;
        double G = invalid<double>;

        Hypothesis(
            const vector< vector<bool> >& connectivityMatrix,
            double G) :
            connectivityMatrix(connectivityMatrix),
            G(G)
        {}

        Hypothesis() {}

        bool operator<(const Hypothesis& that) const
        {
            return G < that.G;
        }
    };

    vector<Hypothesis> hypotheses;

    static bool isForwardInjective(const vector< vector<bool> >& connectivityMatrix);
    static bool isBackwardInjective(const vector< vector<bool> >& connectivityMatrix);

private:
    void run(
        const vector< vector<double> >& tangleMatrix,
        double epsilon,
        bool onlyConsiderInjective,
        bool onlyConsiderPermutation);
};
