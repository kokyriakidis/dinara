#pragma once

#include <boost/pending/disjoint_sets.hpp>

#include "vector.hpp"

namespace dinara {
    class Shasta2DisjointSets;
}

class dinara::Shasta2DisjointSets {
public:
    explicit Shasta2DisjointSets(uint64_t n);

    void initializeDisconnected();
    void initializeDisconnected(uint64_t i);

    void unionSet(uint64_t, uint64_t);
    void link(uint64_t, uint64_t);
    uint64_t findSet(uint64_t);

    void gatherComponents(
        uint64_t minComponentSize,
        vector< vector<uint64_t> >&);

private:
    vector<uint64_t> rank;
    vector<uint64_t> parent;
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets;
};

