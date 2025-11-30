#include "SimpleConsensusCaller.hpp"
#include "Coverage.hpp"
using namespace dinara;


Consensus SimpleConsensusCaller::operator()(
    const Coverage& coverage) const
{
    const AlignedBase base = coverage.mostFrequentBase();
    const size_t repeatCount = coverage.mostFrequentRepeatCount(base);
    return Consensus(base, repeatCount);
}

