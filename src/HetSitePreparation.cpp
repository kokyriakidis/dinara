#include "HetSitePreparation.hpp"

#include <algorithm>
#include <atomic>
#include <thread>

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
    const AnchorOccupancy& occupied)
{
    uint64_t removed = 0;
    for(Assembler::CigarSnpSite& site: sites) {
        for(auto& members: site.alleles) {
            const size_t before = members.size();
            members.erase(std::remove_if(members.begin(), members.end(),
                [&](const pair<OrientedReadId, uint32_t>& m) {
                    return occupied.contains(m.first.getValue(), m.second);
                }), members.end());
            removed += before - members.size();
        }
    }
    return removed;
}



AnchorOccupancy dinara::AnchorOccupancy::forTesting(
    uint64_t orientedReadCount,
    const vector<std::pair<uint64_t, uint32_t>>& entries)
{
    AnchorOccupancy occupancy;
    occupancy.begin_.assign(orientedReadCount + 1, 0);
    for(const auto& [v, position] : entries) {
        static_cast<void>(position);
        if(v < orientedReadCount) ++occupancy.begin_[v + 1];
    }
    for(uint64_t v = 0; v < orientedReadCount; v++) {
        occupancy.begin_[v + 1] += occupancy.begin_[v];
    }
    occupancy.positions_.resize(occupancy.begin_[orientedReadCount]);
    vector<uint64_t> cursor(occupancy.begin_.begin(), occupancy.begin_.end() - 1);
    for(const auto& [v, position] : entries) {
        if(v < orientedReadCount) occupancy.positions_[cursor[v]++] = position;
    }
    for(uint64_t v = 0; v < orientedReadCount; v++) {
        std::sort(occupancy.positions_.begin() + std::ptrdiff_t(occupancy.begin_[v]),
                  occupancy.positions_.begin() + std::ptrdiff_t(occupancy.begin_[v + 1]));
    }
    return occupancy;
}


AnchorOccupancy dinara::buildOccupiedPositions(
    const Shasta2Anchors& anchors, uint64_t threadCount)
{
    if(threadCount == 0) threadCount = std::thread::hardware_concurrency();
    if(threadCount == 0) threadCount = 1;

    const uint64_t anchorCount = anchors.size();
    uint64_t orientedReadCount = 0;
    for(Shasta2AnchorId id = 0; id < anchorCount; id++) {
        for(const Shasta2AnchorMarkerInfo& mi : anchors[id]) {
            orientedReadCount = std::max<uint64_t>(
                orientedReadCount, uint64_t(mi.orientedReadId.getValue()) + 1);
        }
    }

    AnchorOccupancy occupancy;
    occupancy.begin_.assign(orientedReadCount + 1, 0);

    const auto runOverAnchors = [&](auto&& body) {
        const uint64_t chunk = (anchorCount + threadCount - 1) / threadCount;
        vector<std::thread> threads;
        threads.reserve(threadCount);
        for(uint64_t t = 0; t < threadCount; t++) {
            threads.emplace_back([&, t]() {
                body(t * chunk, std::min(anchorCount, (t + 1) * chunk));
            });
        }
        for(std::thread& thread : threads) thread.join();
    };

    // Count, prefix-sum, fill, sort -- the same shape as the marker CSR. Fill
    // order inside a read's slice is nondeterministic and irrelevant: the slice
    // is sorted immediately afterwards.
    {
        vector<std::atomic<uint64_t>> counts(orientedReadCount);
        for(std::atomic<uint64_t>& c : counts) c.store(0, std::memory_order_relaxed);
        runOverAnchors([&](uint64_t begin, uint64_t end) {
            for(Shasta2AnchorId id = begin; id < end; id++) {
                for(const Shasta2AnchorMarkerInfo& mi : anchors[id]) {
                    counts[mi.orientedReadId.getValue()].fetch_add(
                        1, std::memory_order_relaxed);
                }
            }
        });
        uint64_t total = 0;
        for(uint64_t v = 0; v < orientedReadCount; v++) {
            occupancy.begin_[v] = total;
            total += counts[v].load(std::memory_order_relaxed);
        }
        occupancy.begin_[orientedReadCount] = total;
        occupancy.positions_.resize(total);
    }
    {
        vector<std::atomic<uint64_t>> cursor(orientedReadCount);
        for(uint64_t v = 0; v < orientedReadCount; v++) {
            cursor[v].store(occupancy.begin_[v], std::memory_order_relaxed);
        }
        runOverAnchors([&](uint64_t begin, uint64_t end) {
            for(Shasta2AnchorId id = begin; id < end; id++) {
                for(const Shasta2AnchorMarkerInfo& mi : anchors[id]) {
                    occupancy.positions_[cursor[mi.orientedReadId.getValue()].fetch_add(
                        1, std::memory_order_relaxed)] = mi.position;
                }
            }
        });
    }
    {
        const uint64_t chunk = (orientedReadCount + threadCount - 1) / threadCount;
        vector<std::thread> threads;
        threads.reserve(threadCount);
        for(uint64_t t = 0; t < threadCount; t++) {
            threads.emplace_back([&, t]() {
                const uint64_t b = t * chunk;
                const uint64_t e = std::min(orientedReadCount, (t + 1) * chunk);
                for(uint64_t v = b; v < e; v++) {
                    std::sort(
                        occupancy.positions_.begin() + std::ptrdiff_t(occupancy.begin_[v]),
                        occupancy.positions_.begin() + std::ptrdiff_t(occupancy.begin_[v + 1]));
                }
            });
        }
        for(std::thread& thread : threads) thread.join();
    }

    return occupancy;
}
