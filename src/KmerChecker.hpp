#ifndef DINARA_KMER_CHECKER_HPP
#define DINARA_KMER_CHECKER_HPP

// Dinara.
#include "dinaraTypes.hpp"

namespace dinara {
    class KmerChecker;
    class HashedKmerChecker;
}



// The KmerChecker is an abstract class that knows how to find
// out if a k-mer is a marker.
// All implementations must guarantee that if a KmerId if a marker
// its reverse complement is also a marker.
class dinara::KmerChecker {
public:
    virtual bool isMarker(KmerId) const = 0;
};

#endif
