#include "KmerChecker.hpp"
#include "KmerCounter.hpp"
#include <memory>

namespace dinara {
    class FrequencyChecker;
}

class dinara::FrequencyChecker : public dinara::KmerChecker {
public:
    std::shared_ptr<dinara::KmerCounter> kmerCounter;
    uint64_t minFreq;
    uint64_t maxFreq;

    FrequencyChecker(
        std::shared_ptr<dinara::KmerCounter> kmerCounter,
        uint64_t minFreq,
        uint64_t maxFreq) :
        kmerCounter(kmerCounter),
        minFreq(minFreq),
        maxFreq(maxFreq) {}

    // Check if a KmerId is a marker based on its global frequency.
    bool isMarker(KmerId kmerId) const override {
        // Check frequency directly from KmerCounter (O(1) lookup).
        uint64_t frequency = kmerCounter->getFrequency(kmerId);
        return (frequency >= minFreq && frequency <= maxFreq);
    }
};

#endif
