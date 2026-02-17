#include "Shasta2DisjointSets.hpp"
#include "utility.hpp"
#include "orderPairs.hpp"

using namespace dinara;
using namespace std;

Shasta2DisjointSets::Shasta2DisjointSets(uint64_t n) :
    rank(n),
    parent(n),
    disjointSets(&rank[0], &parent[0])
{
    initializeDisconnected();
}

void Shasta2DisjointSets::initializeDisconnected()
{
    const uint64_t n = rank.size();
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }
}

void Shasta2DisjointSets::initializeDisconnected(uint64_t i)
{
    disjointSets.make_set(i);
}

void Shasta2DisjointSets::unionSet(uint64_t i, uint64_t j)
{
    disjointSets.union_set(i, j);
}

void Shasta2DisjointSets::link(uint64_t i, uint64_t j)
{
    disjointSets.link(i, j);
}

uint64_t Shasta2DisjointSets::findSet(uint64_t i)
{
    return disjointSets.find_set(i);
}

void Shasta2DisjointSets::gatherComponents(
    uint64_t minComponentSize,
    vector< vector<uint64_t> >& components)
{
    const uint64_t n = rank.size();
    vector< vector<uint64_t> > allComponents(n);
    for(uint64_t i=0; i<n; i++) {
        const uint64_t componentId = findSet(i);
        allComponents[componentId].push_back(i);
    }

    vector< pair<uint64_t, uint64_t> > componentTable;
    for(uint64_t i=0; i<n; i++) {
        const vector<uint64_t>& component = allComponents[i];
        if(component.size() >= minComponentSize) {
            componentTable.emplace_back(i, component.size());
        }
    }

    std::ranges::sort(componentTable, OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    components.clear();
    for(const auto& p: componentTable) {
        const uint64_t i = p.first;
        const vector<uint64_t>& component = allComponents[i];
        components.emplace_back(component);
    }
}
