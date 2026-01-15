#ifndef DINARA_MINIMIZER_CHECKER_HPP
#define DINARA_MINIMIZER_CHECKER_HPP

#include "KmerChecker.hpp"
#include <vector>
#include <algorithm>

namespace dinara {
    class MinimizerChecker;
}

class dinara::MinimizerChecker : public dinara::KmerChecker {
public:
    // Sorted vector of allowed KmerIds for efficient binary search.
    std::vector<KmerId> allowedKmers;

    // Construct from a vector of KmerIds. 
    // The input vector does not need to be sorted; we sort it here.
    MinimizerChecker(std::vector<KmerId>& kmers) {
        allowedKmers.swap(kmers);
        std::sort(allowedKmers.begin(), allowedKmers.end());
    }

    // Check if a KmerId is a marker.
    bool isMarker(KmerId kmerId) const override {
        return std::binary_search(allowedKmers.begin(), allowedKmers.end(), kmerId);
    }
};

#endif
