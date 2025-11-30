#ifndef DINARA_KMER_CHECKER_FACTORY_HPP
#define DINARA_KMER_CHECKER_FACTORY_HPP

// Dinara.
#include "KmerChecker.hpp"
#include "memory.hpp"

namespace dinara {
    class KmerCheckerFactory;

    class KmerChecker;
    class KmersOptions;
    class Reads;
    class MappedMemoryOwner;

}



// The KmerCheckerFactory knows how to create the appropriate
// type of KmerChecker for the options used.
class dinara::KmerCheckerFactory {
public:

    static shared_ptr<KmerChecker> createNew(
        const KmersOptions&,
        uint64_t threadCount,
        const Reads&,
        const MappedMemoryOwner&);

    static shared_ptr<KmerChecker> createFromBinaryData(
        uint64_t k,
        uint64_t generationMethod,
        const Reads&,
        const MappedMemoryOwner&);
};

#endif

