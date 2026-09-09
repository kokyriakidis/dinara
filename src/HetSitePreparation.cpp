#include "HetSitePreparation.hpp"

#include <algorithm>

using namespace dinara;



uint64_t dinara::collapseDuplicateLoci(vector<Assembler::CigarSnpSite>& sites)
{
    // Sort keys are computed ONCE per site, not inside the comparator: both
    // scan every member of a site, so evaluating them per comparison would
    // rescan the members O(log n) times each for nothing.
    vector<uint64_t> order(sites.size());
    vector<pair<uint64_t, pair<uint64_t, uint32_t>>> key(sites.size());
    for(uint64_t i = 0; i < sites.size(); i++) {
        order[i] = i;
        uint64_t support = 0;
        pair<uint64_t, uint32_t> first{~0ULL, 0};
        for(const auto& allele: sites[i].alleles) {
            support += allele.size();
            for(const auto& m: allele) {
                const pair<uint64_t, uint32_t> k{m.first.getValue(), m.second};
                if(k < first) first = k;
            }
        }
        key[i] = {support, first};
    }
    std::sort(order.begin(), order.end(), [&](uint64_t x, uint64_t y) {
        if(key[x].first != key[y].first) {
            return key[x].first > key[y].first;      // better supported first
        }
        return key[x].second < key[y].second;        // then deterministic
    });

    std::unordered_map<uint64_t, uint64_t> claim;    // occurrence -> site
    vector<bool> dropped(sites.size(), false);
    uint64_t duplicateSites = 0;
    for(const uint64_t i: order) {
        bool clash = false;
        for(const auto& allele: sites[i].alleles) {
            for(const auto& m: allele) {
                const uint64_t occurrence =
                    (uint64_t(m.first.getValue()) << 32) | uint64_t(m.second);
                if(claim.count(occurrence)) { clash = true; break; }
            }
            if(clash) break;
        }
        if(clash) { dropped[i] = true; ++duplicateSites; continue; }
        for(const auto& allele: sites[i].alleles) {
            for(const auto& m: allele) {
                const uint64_t occurrence =
                    (uint64_t(m.first.getValue()) << 32) | uint64_t(m.second);
                claim.emplace(occurrence, i);
            }
        }
    }

    // Rebuild in the DETERMINISTIC order, not in the order detection happened
    // to append. detectCigarSnpSites fills its output from parallel workers
    // under a mutex, so that order depends on thread scheduling -- and it
    // becomes the het anchor ids, which are the last tie-break both here and in
    // the journey rebuild. Two runs on identical input once resolved ties
    // differently and produced different graphs because of exactly this.
    vector<uint64_t> keepOrder;
    keepOrder.reserve(sites.size() - duplicateSites);
    for(uint64_t i = 0; i < sites.size(); i++) {
        if(!dropped[i]) keepOrder.push_back(i);
    }
    std::sort(keepOrder.begin(), keepOrder.end(), [&](uint64_t x, uint64_t y) {
        return key[x].second < key[y].second;
    });
    vector<Assembler::CigarSnpSite> kept;
    kept.reserve(keepOrder.size());
    for(const uint64_t i: keepOrder) {
        kept.push_back(std::move(sites[i]));
    }
    sites.swap(kept);

    return duplicateSites;
}



uint64_t dinara::dropAlreadyAnchoredArmMembers(
    vector<Assembler::CigarSnpSite>& sites,
    const std::unordered_map<uint64_t, Shasta2AnchorId>& occupied,
    uint32_t hetKHalf)
{
    uint64_t removed = 0;
    for(Assembler::CigarSnpSite& site: sites) {
        for(auto& members: site.alleles) {
            const size_t before = members.size();
            members.erase(std::remove_if(members.begin(), members.end(),
                [&](const pair<OrientedReadId, uint32_t>& m) {
                    return occupied.count(
                        (uint64_t(m.first.getValue()) << 32) |
                        uint64_t(m.second + hetKHalf)) != 0;
                }), members.end());
            removed += before - members.size();
        }
    }
    return removed;
}



std::unordered_map<uint64_t, Shasta2AnchorId> dinara::buildOccupiedPositions(
    const Shasta2Anchors& anchors)
{
    std::unordered_map<uint64_t, Shasta2AnchorId> occupied;
    const uint64_t anchorCount = anchors.size();
    for(Shasta2AnchorId id = 0; id < anchorCount; id++) {
        for(const Shasta2AnchorMarkerInfo& mi : anchors[id]) {
            occupied.emplace(
                (uint64_t(mi.orientedReadId.getValue()) << 32) |
                uint64_t(mi.position), id);
        }
    }
    return occupied;
}
