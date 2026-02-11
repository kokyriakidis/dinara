// Dinara.
#include "Assembler.hpp"
#include "Kmer.hpp"
#include "KmerCounter.hpp"
#include "Reads.hpp"
#include "extractKmer.hpp"
#include "findMarkerId.hpp"
#include "timestamp.hpp"
#include "mode3-Anchor.hpp"

// Standard library.
#include "algorithm.hpp"
#include <atomic>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>
#include <optional>
#include <thread>
#include <unordered_map>
#include <unordered_set>

using namespace dinara;
using namespace std;

namespace {
    inline MarkerGraph::VertexId asVertexId(MarkerGraph::CompressedVertexId x)
    {
        return MarkerGraph::VertexId(uint64_t(x));
    }

    using Interval = dinara::mode3::AnchorMarkerInterval;

    inline Interval reverseComplementInterval(
        Interval interval,
        const dinara::MemoryMapped::VectorOfVectors<dinara::CompressedMarker, uint64_t>& markers)
    {
        interval.orientedReadId.flipStrand();
        const uint64_t markerCount = markers.size(interval.orientedReadId.getValue());
        interval.ordinal0 = uint32_t(markerCount) - 1 - interval.ordinal0;
        return interval;
    }

    inline vector<Interval> reverseComplementAnchor(
        const vector<Interval>& anchor,
        const dinara::MemoryMapped::VectorOfVectors<dinara::CompressedMarker, uint64_t>& markers)
    {
        vector<Interval> rc = anchor;
        for(auto& interval : rc) {
            interval = reverseComplementInterval(interval, markers);
        }
        std::sort(rc.begin(), rc.end(), [](const Interval& a, const Interval& b) {
            return a.orientedReadId < b.orientedReadId;
        });
        return rc;
    }

    Kmer getMarkerKmer(
        const dinara::Reads& reads,
        const dinara::MemoryMapped::VectorOfVectors<dinara::CompressedMarker, uint64_t>& markers,
        uint64_t k,
        dinara::OrientedReadId orientedReadId,
        uint32_t ordinal)
    {
        const dinara::ReadId readId = orientedReadId.getReadId();
        const dinara::Strand strand = orientedReadId.getStrand();
        const dinara::LongBaseSequenceView readSequence = reads.getRead(readId);

        if(strand == 0) {
            const auto orientedMarkers = markers[orientedReadId.getValue()];
            const uint32_t position = orientedMarkers[ordinal].position;
            Kmer kmer;
            extractKmer(readSequence, position, k, kmer);
            return kmer;
        }

        // Strand 1: use the corresponding marker on strand 0 and reverse-complement.
        const dinara::OrientedReadId orientedReadId0(readId, 0);
        const auto orientedMarkers0 = markers[orientedReadId0.getValue()];
        const uint32_t markerCount0 = uint32_t(orientedMarkers0.size());
        const uint32_t ordinal0 = markerCount0 - 1 - ordinal;
        const uint32_t position0 = orientedMarkers0[ordinal0].position;
        Kmer kmer0;
        extractKmer(readSequence, position0, k, kmer0);
        return kmer0.reverseComplement(k);
    }

    bool mapMarkerOrdinalByOffsetAndKmer(
        const dinara::Reads& reads,
        const dinara::MemoryMapped::VectorOfVectors<dinara::CompressedMarker, uint64_t>& markers,
        uint64_t k,
        const dinara::AlignmentInfo& info,
        dinara::OrientedReadId orientedReadId0,
        uint32_t ordinal0,
        dinara::OrientedReadId orientedReadId1,
        uint32_t& ordinal1Out,
        const Kmer& seedKmer,
        uint32_t maxSearchRadius,
        uint32_t maxOffsetRange)
    {
        (void)orientedReadId0;

        if(ordinal0 < info.data[0].firstOrdinal || ordinal0 > info.data[0].lastOrdinal) {
            return false;
        }

        // ord0 - ord1 = offset  =>  ord1 = ord0 - offset
        int64_t lo = int64_t(ordinal0) - int64_t(info.maxOrdinalOffset);
        int64_t hi = int64_t(ordinal0) - int64_t(info.minOrdinalOffset);
        if(lo > hi) {
            std::swap(lo, hi);
        }

        // Clamp to aligned range on read1.
        lo = std::max<int64_t>(lo, int64_t(info.data[1].firstOrdinal));
        hi = std::min<int64_t>(hi, int64_t(info.data[1].lastOrdinal));
        if(lo > hi) {
            return false;
        }

        if(uint64_t(hi - lo) > maxOffsetRange) {
            return false;
        }

        const int64_t center = int64_t(ordinal0) - int64_t(info.averageOrdinalOffset);

        auto tryOrdinal1 = [&](int64_t candidate) -> bool {
            if(candidate < lo || candidate > hi) {
                return false;
            }
            const uint32_t o1 = uint32_t(candidate);
            const Kmer kmer1 = getMarkerKmer(reads, markers, k, orientedReadId1, o1);
            if(kmer1 == seedKmer) {
                ordinal1Out = o1;
                return true;
            }
            return false;
        };

        // Try center, then expand out.
        if(tryOrdinal1(center)) {
            return true;
        }
        for(uint32_t d=1; d<=maxSearchRadius; ++d) {
            if(tryOrdinal1(center - int64_t(d))) {
                return true;
            }
            if(tryOrdinal1(center + int64_t(d))) {
                return true;
                    }
                }
        return false;
    }

    struct OverlapEvent {
        uint32_t ordinal = 0;
        int8_t delta = 0; // +1 start, -1 end (end = lastOrdinal+1)
        uint32_t edgeId = 0;
    };

    struct OverlapInterval {
        uint32_t start = 0;
        uint32_t end = 0; // one past last
    };

    inline bool findOverlapIntervalIndex(
        const vector<OverlapInterval>& intervals,
        uint32_t ordinal,
        uint32_t& indexOut)
    {
        if(intervals.empty()) {
            return false;
        }
        auto it = std::upper_bound(
            intervals.begin(),
            intervals.end(),
            ordinal,
            [](uint32_t value, const OverlapInterval& x) {
                return value < x.start;
            });
        if(it == intervals.begin()) {
            return false;
        }
        --it;
        if(ordinal >= it->end) {
            return false;
        }
        indexOut = uint32_t(it - intervals.begin());
        return true;
    }

    struct ArticulationResult {
        vector<uint8_t> isArticulation;
    };

    ArticulationResult findArticulationPoints(const vector<vector<uint32_t>>& adj)
    {
        const uint32_t n = uint32_t(adj.size());
        ArticulationResult result;
        result.isArticulation.assign(n, 0);
        if(n == 0) {
            return result;
        }

        vector<int32_t> disc(n, -1);
        vector<int32_t> low(n, -1);
        vector<int32_t> parent(n, -1);
        int32_t time = 0;

        std::function<void(uint32_t)> dfs = [&](uint32_t u) {
            disc[u] = low[u] = time++;
            uint32_t childCount = 0;
            for(const uint32_t v : adj[u]) {
                if(disc[v] == -1) {
                    parent[v] = int32_t(u);
                    ++childCount;
                    dfs(v);
                    low[u] = std::min(low[u], low[v]);

                    if(parent[u] == -1 && childCount > 1) {
                        result.isArticulation[u] = 1;
                    }
                    if(parent[u] != -1 && low[v] >= disc[u]) {
                        result.isArticulation[u] = 1;
                    }
                } else if(int32_t(v) != parent[u]) {
                    low[u] = std::min(low[u], disc[v]);
                }
            }
        };

        for(uint32_t i=0; i<n; ++i) {
            if(disc[i] == -1) {
                dfs(i);
            }
        }
        return result;
    }

    vector<vector<uint32_t>> connectedComponents(
        const vector<vector<uint32_t>>& adj,
        const vector<uint8_t>& removed)
    {
        const uint32_t n = uint32_t(adj.size());
        vector<uint8_t> visited(n, 0);
        vector<vector<uint32_t>> comps;
        comps.reserve(n);

        vector<uint32_t> stack;
        stack.reserve(n);

        for(uint32_t i=0; i<n; ++i) {
            if(removed[i] || visited[i]) {
                continue;
            }
            visited[i] = 1;
            stack.clear();
            stack.push_back(i);
            comps.emplace_back();
            auto& comp = comps.back();

            while(!stack.empty()) {
                const uint32_t u = stack.back();
                stack.pop_back();
                comp.push_back(u);
                for(const uint32_t v : adj[u]) {
                    if(removed[v] || visited[v]) {
                        continue;
                    }
                    visited[v] = 1;
                    stack.push_back(v);
                }
            }
        }
        return comps;
    }

    uint32_t countCommonSortedNeighbors(const vector<uint32_t>& a, const vector<uint32_t>& b)
    {
        uint32_t common = 0;
        size_t i = 0;
        size_t j = 0;
        while(i < a.size() && j < b.size()) {
            const uint32_t x = a[i];
            const uint32_t y = b[j];
            if(x == y) {
                ++common;
                ++i;
                ++j;
            } else if(x < y) {
                ++i;
            } else {
                ++j;
            }
        }
        return common;
    }

    // Peel a quasi-clique core from a candidate group.
    // Removes nodes that do not have enough neighbors *within the group*.
    // This is robust against a small number of bridging/chimeric reads.
    //
    // tau is represented as a rational num/denom (for speed and determinism):
    // keep node u only if deg_in_group(u) >= ceil( (num/denom) * (k-1) ),
    // where k is the current group size.
    vector<uint32_t> peelQuasiCliqueCore(
        const vector<vector<uint32_t>>& adj,
        const vector<uint32_t>& group,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        uint32_t tauNum,
        uint32_t tauDen)
    {
        const uint32_t n = uint32_t(adj.size());
        const uint32_t k0 = uint32_t(group.size());
        if(k0 == 0) {
            return {};
        }
        if(k0 == 1) {
            if(minAnchorCoverage <= 1 && maxAnchorCoverage >= 1) {
                return group;
            } else {
                return {};
            }
        }

        vector<uint8_t> inGroup(n, 0);
        for(const uint32_t u : group) {
            DINARA_ASSERT(u < n);
            inGroup[u] = 1;
        }

        uint32_t k = k0;
        vector<uint32_t> deg(n, 0);
        for(const uint32_t u : group) {
            uint32_t d = 0;
            for(const uint32_t v : adj[u]) {
                if(inGroup[v]) {
                    ++d;
                }
            }
            deg[u] = d;
        }

        auto minDegForSize = [&](uint32_t kk) -> uint32_t {
            if(kk <= 1) {
                return 0;
            }
            // ceil(tau * (kk-1)) with tau = tauNum/tauDen
            const uint32_t rhs = kk - 1;
            return (tauNum * rhs + tauDen - 1) / tauDen;
        };

        vector<uint32_t> queue;
        queue.reserve(k0);

        uint32_t minDeg = minDegForSize(k);
        for(const uint32_t u : group) {
            if(inGroup[u] && deg[u] < minDeg) {
                queue.push_back(u);
            }
        }

        // Iteratively remove nodes with insufficient internal degree.
        while(!queue.empty()) {
            const uint32_t u = queue.back();
            queue.pop_back();
            if(!inGroup[u]) {
                continue;
            }
            if(deg[u] >= minDeg) {
                continue;
            }

            // Remove u.
            inGroup[u] = 0;
            --k;

            if(k < minAnchorCoverage) {
                return {};
            }

            // The threshold increases/decreases with k.
            minDeg = minDegForSize(k);

            // Update neighbors.
            for(const uint32_t v : adj[u]) {
                if(!inGroup[v]) {
                    continue;
                }
                if(deg[v] > 0) {
                    --deg[v];
                }
                if(deg[v] < minDeg) {
                    queue.push_back(v);
                }
            }
        }

        vector<uint32_t> core;
        core.reserve(k);
        for(const uint32_t u : group) {
            if(inGroup[u]) {
                core.push_back(u);
            }
        }

        if(core.size() < minAnchorCoverage || core.size() > maxAnchorCoverage) {
            return {};
        }

        return core;
    }

    vector<vector<uint32_t>> peelGroups(
        const vector<vector<uint32_t>>& adj,
        const vector<vector<uint32_t>>& groups,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        // Default quasi-clique density.
        // Use 4/5 (0.8) to allow for some missing edges after readGraph filtering.
        constexpr uint32_t tauNum = 4;
        constexpr uint32_t tauDen = 5;

        vector<vector<uint32_t>> out;
        out.reserve(groups.size());
        for(const auto& g : groups) {
            auto core = peelQuasiCliqueCore(adj, g, minAnchorCoverage, maxAnchorCoverage, tauNum, tauDen);
            if(!core.empty()) {
                out.push_back(std::move(core));
            }
        }
        return out;
    }

    // Greedy dense-core extraction.
    // This is used as a fallback when topology-based splitting (articulation/triangles/Jaccard)
    // fails to separate two dense regions connected by multiple bridge reads.
    vector<vector<uint32_t>> greedyDenseCoreSplit(
        const vector<vector<uint32_t>>& adj,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        const uint32_t n = uint32_t(adj.size());
        if(n < 2*minAnchorCoverage) {
            return {};
        }

        // Precompute neighbor lists are assumed sorted/unique by callers.
        vector<uint8_t> alive(n, 1);
        uint32_t aliveCount = n;

        // Local clustering coefficient score per node, computed as:
        // cc(u) = (2 * edges among alive neighbors of u) / (deg(u)*(deg(u)-1)).
        // We avoid floating point by comparing the fraction with numerator=sumCommon,
        // denominator=deg*(deg-1), where sumCommon is Σ_{v in N(u)} |N(u)∩N(v)|.
        auto clusteringScore = [&](uint32_t u, vector<uint32_t>& neighborsAlive) -> pair<uint64_t, uint64_t> {
            neighborsAlive.clear();
            for(const uint32_t v : adj[u]) {
                if(alive[v]) {
                    neighborsAlive.push_back(v);
                }
            }
            const uint64_t deg = neighborsAlive.size();
            if(deg < 2) {
                return {0, 1};
            }
            uint64_t sumCommon = 0;
            for(const uint32_t v : neighborsAlive) {
                sumCommon += countCommonSortedNeighbors(neighborsAlive, adj[v]);
            }
            const uint64_t denom = deg * (deg - 1);
            return {sumCommon, denom};
        };

        vector<uint32_t> remaining;
        remaining.reserve(n);
        for(uint32_t i=0; i<n; ++i) remaining.push_back(i);

        vector<vector<uint32_t>> cores;
        cores.reserve(4);

        while(aliveCount >= minAnchorCoverage) {
            // Pick the best seed among alive nodes: maximize clustering coefficient, then degree.
            uint32_t best = invalid<uint32_t>;
            uint64_t bestNum = 0;
            uint64_t bestDen = 1;
            uint32_t bestDeg = 0;
            vector<uint32_t> neighborsAlive;
            neighborsAlive.reserve(128);
            for(uint32_t u=0; u<n; ++u) {
                if(!alive[u]) continue;
                // Degree among alive nodes.
                uint32_t du = 0;
                for(const uint32_t v : adj[u]) {
                    if(alive[v]) ++du;
                }
                if(du < minAnchorCoverage - 1) continue;
                const auto [num, den] = clusteringScore(u, neighborsAlive);
                const bool better =
                    (best == invalid<uint32_t>) ||
                    (num * bestDen > bestNum * den) ||
                    (num * bestDen == bestNum * den && du > bestDeg);
                if(better) {
                    best = u;
                    bestNum = num;
                    bestDen = den;
                    bestDeg = du;
                }
            }
            if(best == invalid<uint32_t>) {
                break;
            }

            // Candidate = seed + its alive neighbors.
            vector<uint32_t> cand;
            cand.reserve(adj[best].size() + 1);
            cand.push_back(best);
            for(const uint32_t v : adj[best]) {
                if(alive[v]) cand.push_back(v);
            }
            std::sort(cand.begin(), cand.end());
            cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

            // Peel a quasi-clique core from candidate.
            auto core = peelQuasiCliqueCore(adj, cand, minAnchorCoverage, maxAnchorCoverage, 4, 5);
            if(core.empty()) {
                // Seed not useful. Remove it and continue.
                alive[best] = 0;
                --aliveCount;
                continue;
            }

            // Remove core nodes from alive set.
            for(const uint32_t u : core) {
                if(alive[u]) {
                    alive[u] = 0;
                    --aliveCount;
                }
            }
            cores.push_back(std::move(core));

            // Stop if we already have at least 2 cores and too few nodes remain to form another.
            if(cores.size() >= 2 && aliveCount < minAnchorCoverage) {
                break;
            }
        }

        if(cores.size() >= 2) {
            return cores;
        }
        return {};
    }

    struct MclMatrix {
        uint32_t n = 0;
        vector<double> a; // row-major, size n*n
        explicit MclMatrix(uint32_t n = 0) : n(n), a(size_t(n) * size_t(n), 0.0) {}
        inline double& operator()(uint32_t r, uint32_t c) {
            return a[size_t(r) * size_t(n) + size_t(c)];
        }
        inline double operator()(uint32_t r, uint32_t c) const {
            return a[size_t(r) * size_t(n) + size_t(c)];
        }
    };

    vector<vector<uint32_t>> symmetrizeAdj(const vector<vector<uint32_t>>& adj)
    {
        const uint32_t n = uint32_t(adj.size());
        vector<vector<uint32_t>> und(n);
        for(uint32_t u=0; u<n; ++u) {
            und[u].reserve(adj[u].size());
        }
        for(uint32_t u=0; u<n; ++u) {
            for(const uint32_t v : adj[u]) {
                if(v >= n || v == u) {
                    continue;
                }
                und[u].push_back(v);
                und[v].push_back(u);
            }
        }
        for(uint32_t u=0; u<n; ++u) {
            auto& nbr = und[u];
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
        }
        return und;
    }

    double graphDensityUndirected(const vector<vector<uint32_t>>& undAdj)
    {
        const uint32_t n = uint32_t(undAdj.size());
        if(n < 2) {
            return 0.0;
        }
        uint64_t sumDeg = 0;
        for(uint32_t u=0; u<n; ++u) {
            sumDeg += undAdj[u].size();
        }
        const double e = double(sumDeg) / 2.0;
        return (2.0 * e) / (double(n) * double(n - 1));
    }

    double averageLocalClusteringCoefficientUndirected(const vector<vector<uint32_t>>& undAdj)
    {
        const uint32_t n = uint32_t(undAdj.size());
        if(n < 3) {
            return 0.0;
        }

        // Build a dense boolean adjacency for triangle counting.
        vector<uint8_t> mat(size_t(n) * size_t(n), 0);
        auto at = [&](uint32_t r, uint32_t c) -> uint8_t& {
            return mat[size_t(r) * size_t(n) + size_t(c)];
        };
        for(uint32_t u=0; u<n; ++u) {
            for(const uint32_t v : undAdj[u]) {
                at(u, v) = 1;
            }
        }

        double sum = 0.0;
        uint32_t count = 0;
        for(uint32_t u=0; u<n; ++u) {
            const uint32_t deg = uint32_t(undAdj[u].size());
            if(deg < 2) {
                continue;
            }
            uint64_t triEdges = 0;
            const auto& nbr = undAdj[u];
            for(size_t i=0; i<nbr.size(); ++i) {
                const uint32_t v = nbr[i];
                for(size_t j=i+1; j<nbr.size(); ++j) {
                    const uint32_t w = nbr[j];
                    if(at(v, w)) {
                        ++triEdges;
                    }
                }
            }
            const double possible = double(deg) * double(deg - 1) / 2.0;
            sum += double(triEdges) / possible;
            ++count;
        }
        if(count == 0) {
            return 0.0;
        }
        return sum / double(count);
    }

    void normalizeColumns(MclMatrix& m)
    {
        const uint32_t n = m.n;
        for(uint32_t c=0; c<n; ++c) {
            double s = 0.0;
            for(uint32_t r=0; r<n; ++r) {
                s += m(r, c);
            }
            if(s == 0.0) {
                // Ensure column-stochastic.
                for(uint32_t r=0; r<n; ++r) {
                    m(r, c) = 0.0;
                }
                m(c, c) = 1.0;
                s = 1.0;
            }
            const double inv = 1.0 / s;
            for(uint32_t r=0; r<n; ++r) {
                m(r, c) *= inv;
            }
        }
    }

    MclMatrix multiply(const MclMatrix& a, const MclMatrix& b)
    {
        DINARA_ASSERT(a.n == b.n);
        const uint32_t n = a.n;
        MclMatrix out(n);
        for(uint32_t i=0; i<n; ++i) {
            for(uint32_t k=0; k<n; ++k) {
                const double aik = a(i, k);
                if(aik == 0.0) {
                    continue;
                }
                for(uint32_t j=0; j<n; ++j) {
                    out(i, j) += aik * b(k, j);
                }
            }
        }
        return out;
    }

    void inflateAndPrune(MclMatrix& m, double inflation, double pruneThreshold)
    {
        const uint32_t n = m.n;
        for(uint32_t c=0; c<n; ++c) {
            for(uint32_t r=0; r<n; ++r) {
                double x = m(r, c);
                if(x <= 0.0) {
                    continue;
                }
                x = std::pow(x, inflation);
                if(x < pruneThreshold) {
                    x = 0.0;
                }
                m(r, c) = x;
            }
        }
    }

    double maxAbsDiff(const MclMatrix& a, const MclMatrix& b)
    {
        DINARA_ASSERT(a.n == b.n);
        double d = 0.0;
        for(size_t i=0; i<a.a.size(); ++i) {
            d = std::max(d, std::abs(a.a[i] - b.a[i]));
        }
        return d;
    }

    vector<vector<uint32_t>> mclClusterUndirected(
        const vector<vector<uint32_t>>& undAdj,
        double inflation,
        uint32_t maxIterations)
    {
        const uint32_t n = uint32_t(undAdj.size());
        if(n == 0) {
            return {};
        }
        if(n == 1) {
            return { {0} };
        }

        // Initialize matrix with self-loops and adjacency.
        MclMatrix m(n);
        for(uint32_t i=0; i<n; ++i) {
            m(i, i) = 1.0;
        }
        for(uint32_t u=0; u<n; ++u) {
            for(const uint32_t v : undAdj[u]) {
                if(v < n) {
                    m(v, u) = 1.0;
                }
            }
        }
        normalizeColumns(m);

        constexpr double pruneThreshold = 1e-4;
        constexpr double convergenceTol = 1e-3;

        for(uint32_t it=0; it<maxIterations; ++it) {
            MclMatrix next = multiply(m, m);          // expansion (power 2)
            inflateAndPrune(next, inflation, pruneThreshold);
            normalizeColumns(next);
            const double delta = maxAbsDiff(m, next);
            m = std::move(next);
            if(delta < convergenceTol) {
                break;
            }
        }

        // Extract clusters by following the strongest attractor per node until convergence.
        vector<uint32_t> argmax(n, 0);
        for(uint32_t c=0; c<n; ++c) {
            uint32_t best = 0;
            double bestVal = m(0, c);
            for(uint32_t r=1; r<n; ++r) {
                const double v = m(r, c);
                if(v > bestVal) {
                    bestVal = v;
                    best = r;
                }
            }
            argmax[c] = best;
        }

        auto attractor = [&](uint32_t x) -> uint32_t {
            for(uint32_t step=0; step<n; ++step) {
                const uint32_t y = argmax[x];
                if(y == x) {
                    return x;
                }
                x = y;
            }
            return x;
        };

        unordered_map<uint32_t, vector<uint32_t>> groups;
        groups.reserve(n);
        for(uint32_t i=0; i<n; ++i) {
            groups[attractor(i)].push_back(i);
        }

        vector<vector<uint32_t>> out;
        out.reserve(groups.size());
        for(auto& kv : groups) {
            auto& g = kv.second;
            std::sort(g.begin(), g.end());
            out.push_back(std::move(g));
        }

        // Deterministic order.
        std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) {
            if(a.size() != b.size()) return a.size() > b.size();
            return a < b;
        });

        return out;
    }

    struct SmallBitset {
        static constexpr uint32_t MaxWords = 4; // 256 bits
        uint32_t wordCount = 0;
        std::array<uint64_t, MaxWords> w{};

        SmallBitset() = default;
        explicit SmallBitset(uint32_t bits) { reset(bits); }

        void reset(uint32_t bits)
        {
            wordCount = (bits + 63U) / 64U;
            for(uint32_t i=0; i<MaxWords; ++i) w[i] = 0;
        }

        inline bool empty() const
        {
            for(uint32_t i=0; i<wordCount; ++i) {
                if(w[i]) return false;
            }
            return true;
        }

        inline void set(uint32_t i)
        {
            w[i >> 6] |= (1ULL << (i & 63U));
        }

        inline void clear(uint32_t i)
        {
            w[i >> 6] &= ~(1ULL << (i & 63U));
        }

        inline bool test(uint32_t i) const
        {
            return (w[i >> 6] >> (i & 63U)) & 1ULL;
        }

        inline uint32_t popcount() const
        {
            uint32_t c = 0;
            for(uint32_t i=0; i<wordCount; ++i) {
                c += uint32_t(__builtin_popcountll(w[i]));
            }
            return c;
        }

        inline SmallBitset operator&(const SmallBitset& other) const
        {
            SmallBitset out;
            out.wordCount = wordCount;
            for(uint32_t i=0; i<wordCount; ++i) {
                out.w[i] = w[i] & other.w[i];
            }
            return out;
        }

        inline SmallBitset operator|(const SmallBitset& other) const
        {
            SmallBitset out;
            out.wordCount = wordCount;
            for(uint32_t i=0; i<wordCount; ++i) {
                out.w[i] = w[i] | other.w[i];
            }
            return out;
        }

        inline SmallBitset operator-(const SmallBitset& other) const
        {
            SmallBitset out;
            out.wordCount = wordCount;
            for(uint32_t i=0; i<wordCount; ++i) {
                out.w[i] = w[i] & ~other.w[i];
            }
            return out;
        }

        inline void operator|=(const SmallBitset& other)
        {
            for(uint32_t i=0; i<wordCount; ++i) {
                w[i] |= other.w[i];
            }
        }
    };

    template<class F>
    inline void forEachSetBit(const SmallBitset& s, F&& f)
    {
        for(uint32_t wi=0; wi<s.wordCount; ++wi) {
            uint64_t x = s.w[wi];
            while(x) {
                const uint32_t b = uint32_t(__builtin_ctzll(x));
                const uint32_t idx = wi * 64U + b;
                f(idx);
                x &= x - 1;
            }
        }
    }

    vector<SmallBitset> enumerateMaximalCliques(
        const vector<SmallBitset>& nbr,
        uint32_t nodeCount,
        uint32_t maxCliques,
        bool& exceeded)
    {
        exceeded = false;
        vector<SmallBitset> cliques;
        cliques.reserve(32);

        SmallBitset R(nodeCount), P(nodeCount), X(nodeCount);
        for(uint32_t i=0; i<nodeCount; ++i) {
            P.set(i);
        }

        std::function<void(const SmallBitset&, SmallBitset&, SmallBitset&)> rec =
            [&](const SmallBitset& Rcur, SmallBitset& Pcur, SmallBitset& Xcur) {
                if(exceeded) return;
                if(Pcur.empty() && Xcur.empty()) {
                    cliques.push_back(Rcur);
                    if(cliques.size() >= maxCliques) {
                        exceeded = true;
                    }
                    return;
                }

                // Choose pivot u from P ∪ X to maximize |P ∩ N(u)|.
                SmallBitset unionPX = Pcur | Xcur;
                uint32_t pivot = invalid<uint32_t>;
                uint32_t best = 0;
                forEachSetBit(unionPX, [&](uint32_t u) {
                    const uint32_t c = (Pcur & nbr[u]).popcount();
                    if(c > best) {
                        best = c;
                        pivot = u;
                    }
                });
                if(pivot == invalid<uint32_t>) {
                    pivot = 0;
                }

                SmallBitset candidates = Pcur - nbr[pivot];
                vector<uint32_t> candList;
                candList.reserve(candidates.popcount());
                forEachSetBit(candidates, [&](uint32_t v) { candList.push_back(v); });

                for(const uint32_t v : candList) {
                    if(exceeded) return;
                    if(!Pcur.test(v)) {
                        continue;
                    }
                    SmallBitset Rnext = Rcur;
                    Rnext.set(v);
                    SmallBitset Pnext = Pcur & nbr[v];
                    SmallBitset Xnext = Xcur & nbr[v];
                    rec(Rnext, Pnext, Xnext);
                    Pcur.clear(v);
                    Xcur.set(v);
                }
            };

        rec(R, P, X);
        return cliques;
    }

    vector<vector<uint32_t>> splitVertexByCliqueCover(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        const vector<uint8_t>& isCore,
        uint32_t attachMinSupport,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage,
        bool& exploded)
    {
        exploded = false;
        const uint32_t n = uint32_t(adjAll.size());
        if(n == 0) {
            return {};
        }

        const auto& baseAdj = hasAnyCisEdge ? adjCis : adjAll;
        vector<vector<uint32_t>> und = symmetrizeAdj(baseAdj);

        vector<uint32_t> coreNodes;
        coreNodes.reserve(n);
        for(uint32_t u=0; u<n; ++u) {
            if(u < isCore.size() && isCore[u]) {
                coreNodes.push_back(u);
            }
        }
        if(coreNodes.size() < 2) {
            coreNodes.clear();
            for(uint32_t u=0; u<n; ++u) {
                coreNodes.push_back(u);
            }
        }

        const uint32_t cn = uint32_t(coreNodes.size());
        if(cn < 2 || cn > 256) {
            return {};
        }

        vector<int32_t> toCore(n, -1);
        for(uint32_t i=0; i<cn; ++i) {
            toCore[coreNodes[i]] = int32_t(i);
        }

        vector<SmallBitset> nbr(cn);
        for(uint32_t i=0; i<cn; ++i) {
            nbr[i].reset(cn);
        }

        uint64_t sumDeg = 0;
        for(uint32_t i=0; i<cn; ++i) {
            const uint32_t u = coreNodes[i];
            for(const uint32_t v : und[u]) {
                const int32_t j = (v < n) ? toCore[v] : -1;
                if(j >= 0 && uint32_t(j) != i) {
                    nbr[i].set(uint32_t(j));
                }
            }
            sumDeg += nbr[i].popcount();
        }

        // If the core is almost a clique, just keep it (will be peelled later if needed).
        {
            const double e = double(sumDeg) / 2.0;
            const double density = (2.0 * e) / (double(cn) * double(cn - 1));
            if(density >= 0.97) {
                vector<uint32_t> all;
                all.reserve(n);
                for(uint32_t u=0; u<n; ++u) all.push_back(u);
                if(all.size() >= minAnchorCoverage && all.size() <= maxAnchorCoverage) {
                    return {all};
                }
            }
        }

        // Fast path: if the core graph has multiple connected components and each component
        // is a clique, we can use those components directly as disjoint core cliques.
        vector<SmallBitset> bigCliques;
        {
            vector<uint8_t> visited(cn, 0);
            vector<uint32_t> stack;
            stack.reserve(cn);
            for(uint32_t start=0; start<cn; ++start) {
                if(visited[start]) continue;
                visited[start] = 1;
                stack.clear();
                stack.push_back(start);
                vector<uint32_t> comp;
                while(!stack.empty()) {
                    const uint32_t u = stack.back();
                    stack.pop_back();
                    comp.push_back(u);
                    forEachSetBit(nbr[u], [&](uint32_t v) {
                        if(v >= cn) return;
                        if(!visited[v]) {
                            visited[v] = 1;
                            stack.push_back(v);
                        }
                    });
                }
                if(comp.size() < 2) {
                    continue;
                }
                SmallBitset compSet(cn);
                for(const uint32_t u : comp) compSet.set(u);
                bool isClique = true;
                const uint32_t k = uint32_t(comp.size());
                for(const uint32_t u : comp) {
                    const uint32_t d = (nbr[u] & compSet).popcount();
                    if(d != k - 1) {
                        isClique = false;
                        break;
                    }
                }
                if(isClique) {
                    bigCliques.push_back(compSet);
                }
            }
        }

        if(bigCliques.empty()) {
            constexpr uint32_t maxCliques = 256;
            bool exceeded = false;
            auto cliques = enumerateMaximalCliques(nbr, cn, maxCliques, exceeded);
            if(exceeded) {
                exploded = true;
                return {};
            }
            bigCliques.reserve(cliques.size());
            for(const auto& c : cliques) {
                if(c.popcount() >= 2) {
                    bigCliques.push_back(c);
                }
            }
        }
        if(bigCliques.empty()) {
            return {};
        }

        std::sort(bigCliques.begin(), bigCliques.end(),
            [&](const SmallBitset& a, const SmallBitset& b) {
                const uint32_t sa = a.popcount();
                const uint32_t sb = b.popcount();
                if(sa != sb) return sa > sb;
                uint32_t fa = invalid<uint32_t>;
                uint32_t fb = invalid<uint32_t>;
                forEachSetBit(a, [&](uint32_t i) { if(fa == invalid<uint32_t>) fa = i; });
                forEachSetBit(b, [&](uint32_t i) { if(fb == invalid<uint32_t>) fb = i; });
                return fa < fb;
            });

        SmallBitset used(cn);
        used.reset(cn);
        vector<SmallBitset> chosen;
        chosen.reserve(8);
        for(const auto& c : bigCliques) {
            if((c & used).empty()) {
                chosen.push_back(c);
                used |= c;
            }
        }
        if(chosen.empty()) {
            return {};
        }
        // If we could only pick one disjoint clique and it does not cover all core nodes,
        // don't proceed: we'd risk dropping an entire alternate region/haplotype.
        // Let the more robust fallback splitter handle these cases.
        if(chosen.size() == 1 && used.popcount() < cn) {
            return {};
        }

        const uint32_t gCount = uint32_t(chosen.size());
        vector<vector<uint32_t>> groups(gCount);
        vector<int32_t> groupOfNode(n, -1);
        vector<uint32_t> coreSize(gCount, 0);
        for(uint32_t gi=0; gi<gCount; ++gi) {
            forEachSetBit(chosen[gi], [&](uint32_t ci) {
                const uint32_t u = coreNodes[ci];
                groupOfNode[u] = int32_t(gi);
                groups[gi].push_back(u);
                ++coreSize[gi];
            });
        }

        // Attach remaining core nodes based on core-core support.
        constexpr double coreAttachFrac = 0.90;
        vector<uint32_t> counts(gCount);
        for(uint32_t ci=0; ci<cn; ++ci) {
            const uint32_t u = coreNodes[ci];
            if(groupOfNode[u] != -1) continue;
            std::fill(counts.begin(), counts.end(), 0);
            for(const uint32_t v : und[u]) {
                if(v >= n) continue;
                if(!(v < isCore.size() && isCore[v])) continue;
                const int32_t gi = groupOfNode[v];
                if(gi >= 0) ++counts[uint32_t(gi)];
            }
            uint32_t bestGi = std::numeric_limits<uint32_t>::max();
            uint32_t bestCount = 0;
            uint32_t secondCount = 0;
            for(uint32_t gi=0; gi<gCount; ++gi) {
                const uint32_t c = counts[gi];
                if(c > bestCount) {
                    secondCount = bestCount;
                    bestCount = c;
                    bestGi = gi;
                } else if(c == bestCount && c != 0) {
                    secondCount = bestCount;
                } else if(c > secondCount) {
                    secondCount = c;
                }
            }
            if(bestGi == std::numeric_limits<uint32_t>::max()) continue;
            if(secondCount == bestCount) continue;
            uint32_t required = attachMinSupport;
            if(coreSize[bestGi] >= 4) {
                required = std::max(required, uint32_t(std::ceil(coreAttachFrac * double(coreSize[bestGi]))));
            }
            if(bestCount < required) continue;
            if(groups[bestGi].size() >= maxAnchorCoverage) continue;
            groupOfNode[u] = int32_t(bestGi);
            groups[bestGi].push_back(u);
            ++coreSize[bestGi];
        }

        // Attach non-core nodes (typically contained reads) using support to core nodes only.
        constexpr double nonCoreAttachFrac = 0.80;
        for(uint32_t u=0; u<n; ++u) {
            if(u < isCore.size() && isCore[u]) continue;
            if(groupOfNode[u] != -1) continue;

            std::fill(counts.begin(), counts.end(), 0);
            uint32_t degToCores = 0;
            for(const uint32_t v : und[u]) {
                if(v >= n) continue;
                if(!(v < isCore.size() && isCore[v])) continue;
                const int32_t gi = groupOfNode[v];
                if(gi >= 0) {
                    ++counts[uint32_t(gi)];
                    ++degToCores;
                }
            }

            uint32_t bestGi = std::numeric_limits<uint32_t>::max();
            uint32_t bestCount = 0;
            uint32_t secondCount = 0;
            for(uint32_t gi=0; gi<gCount; ++gi) {
                const uint32_t c = counts[gi];
                if(c > bestCount) {
                    secondCount = bestCount;
                    bestCount = c;
                    bestGi = gi;
                } else if(c == bestCount && c != 0) {
                    secondCount = bestCount;
                } else if(c > secondCount) {
                    secondCount = c;
                }
            }
            if(bestGi == std::numeric_limits<uint32_t>::max()) continue;
            if(secondCount == bestCount) continue;
            // Ambiguity/bridge check: if support is close for multiple clusters, do not attach.
            // Integer math for bestCount/secondCount <= 1.25  <=>  bestCount*4 <= secondCount*5.
            if(secondCount != 0 && bestCount * 4 <= secondCount * 5) {
                continue;
            }

            uint32_t required = attachMinSupport;
            if(coreSize[bestGi] >= 6) {
                required = std::max(required, uint32_t(std::ceil(nonCoreAttachFrac * double(coreSize[bestGi]))));
            }
            if(degToCores >= 6) {
                required = std::max(required, uint32_t((degToCores + 3) / 4));
                required = std::min(required, uint32_t(4));
            }

            if(bestCount < required) continue;
            if(groups[bestGi].size() >= maxAnchorCoverage) continue;
            groupOfNode[u] = int32_t(bestGi);
            groups[bestGi].push_back(u);
        }

        vector<vector<uint32_t>> kept;
        kept.reserve(groups.size());
        for(auto& g : groups) {
            std::sort(g.begin(), g.end());
            g.erase(std::unique(g.begin(), g.end()), g.end());
            if(g.size() >= minAnchorCoverage && g.size() <= maxAnchorCoverage) {
                kept.push_back(std::move(g));
            }
        }
        if(kept.empty()) {
            return {};
        }

        return peelGroups(und, kept, minAnchorCoverage, maxAnchorCoverage);
    }

    vector<vector<uint32_t>> splitVertexByOverlapSupport(
        const vector<vector<uint32_t>>& adj,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        const uint32_t n = uint32_t(adj.size());
        if(n == 0) {
            return {};
        }
        if(n == 1) {
            if(minAnchorCoverage <= 1 && maxAnchorCoverage >= 1) {
                return { {0} };
            } else {
                return {};
            }
        }

        // Drop isolated nodes (no overlap support) if doing so keeps a valid anchor.
        vector<uint8_t> removedIsolated(n, 0);
        uint32_t keptAfterIsolated = 0;
        for(uint32_t i=0; i<n; ++i) {
            if(adj[i].empty()) {
                removedIsolated[i] = 1;
            } else {
                ++keptAfterIsolated;
            }
        }
        if(keptAfterIsolated < n &&
           keptAfterIsolated >= minAnchorCoverage && keptAfterIsolated <= maxAnchorCoverage) {
            return { [&]{
                vector<uint32_t> all;
                all.reserve(keptAfterIsolated);
                for(uint32_t i=0; i<n; ++i) {
                    if(!removedIsolated[i]) all.push_back(i);
                }
                return all;
            }() };
        }

        // If already disconnected, keep components.
        {
            vector<uint8_t> none(n, 0);
            const auto comps = connectedComponents(adj, none);
            if(comps.size() > 1) {
                vector<vector<uint32_t>> kept;
                for(const auto& comp : comps) {
                    if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                        kept.push_back(comp);
                    }
                }
                if(!kept.empty()) {
                    auto peeled = peelGroups(adj, kept, minAnchorCoverage, maxAnchorCoverage);
                    if(!peeled.empty()) {
                        return peeled;
                    }
                }
            }
        }

        // Remove articulation points to avoid "single-bridge" collapses.
        const auto art = findArticulationPoints(adj);
        vector<uint8_t> removed = art.isArticulation;
        const auto compsNoArt = connectedComponents(adj, removed);
        vector<vector<uint32_t>> keptNoArt;
        for(const auto& comp : compsNoArt) {
            if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                keptNoArt.push_back(comp);
            }
        }
        if(!keptNoArt.empty()) {
            auto peeled = peelGroups(adj, keptNoArt, minAnchorCoverage, maxAnchorCoverage);
            if(peeled.size() >= 2) {
                return peeled;
            }
        }

        // If the vertex is held together by "weak" cross-edges with divergent neighbor sets,
        // split using only edges that have at least one common neighbor (triangle support).
        // This catches multi-bridge cases that are not single articulation points.
        {
            vector<vector<uint32_t>> strongAdj(n);
            for(uint32_t u=0; u<n; ++u) {
                strongAdj[u].reserve(adj[u].size());
            }
            for(uint32_t u=0; u<n; ++u) {
                for(const uint32_t v : adj[u]) {
                    if(v <= u) {
                        continue;
                    }
                    const uint32_t common = countCommonSortedNeighbors(adj[u], adj[v]);
                    if(common >= 1) {
                        strongAdj[u].push_back(v);
                        strongAdj[v].push_back(u);
                    }
                }
            }
            for(uint32_t u=0; u<n; ++u) {
                auto& nbr = strongAdj[u];
                std::sort(nbr.begin(), nbr.end());
                nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
            }

            vector<uint8_t> none(n, 0);
            const auto compsStrong = connectedComponents(strongAdj, none);
            vector<vector<uint32_t>> keptStrong;
            for(const auto& comp : compsStrong) {
                if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                    keptStrong.push_back(comp);
                }
            }
            if(!keptStrong.empty()) {
                auto peeled = peelGroups(strongAdj, keptStrong, minAnchorCoverage, maxAnchorCoverage);
                if(peeled.size() >= 2) {
                    return peeled;
                }
            }
        }

        // If we still have a single component, try a stricter edge filter based on neighbor-set similarity.
        // This targets cases where 2+ bridge/chimeric reads connect two dense groups:
        // triangle support alone can keep those bridge edges, but the Jaccard similarity of neighbor sets
        // is typically low for bridge edges.
        {
            // Heuristic guard: Jaccard is not meaningful for tiny graphs.
            if(n >= 6) {
                vector<vector<uint32_t>> jaccAdj(n);
                for(uint32_t u=0; u<n; ++u) {
                    jaccAdj[u].reserve(adj[u].size());
                }

                constexpr double minJaccard = 0.6;
                constexpr uint32_t minCommon = 2;

                for(uint32_t u=0; u<n; ++u) {
                    for(const uint32_t v : adj[u]) {
                        if(v <= u) {
                            continue;
                        }
                        const uint32_t common = countCommonSortedNeighbors(adj[u], adj[v]);
                        const uint32_t du = uint32_t(adj[u].size());
                        const uint32_t dv = uint32_t(adj[v].size());
                        const uint32_t uni = du + dv - common;
                        if(uni == 0) {
                            continue;
                        }
                        const double j = double(common) / double(uni);
                        if(common >= minCommon || j >= minJaccard) {
                            jaccAdj[u].push_back(v);
                            jaccAdj[v].push_back(u);
                        }
                    }
                }
                for(uint32_t u=0; u<n; ++u) {
                    auto& nbr = jaccAdj[u];
                    std::sort(nbr.begin(), nbr.end());
                    nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
                }

                vector<uint8_t> none(n, 0);
                const auto compsJ = connectedComponents(jaccAdj, none);
                vector<vector<uint32_t>> keptJ;
                for(const auto& comp : compsJ) {
                    if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                        keptJ.push_back(comp);
                    }
                }
                if(!keptJ.empty()) {
                    auto peeled = peelGroups(jaccAdj, keptJ, minAnchorCoverage, maxAnchorCoverage);
                    if(peeled.size() >= 2) {
                        return peeled;
                    }
                }
            }
        }

        // Greedy attempt to extract 2+ dense cores even when the topology-based split did not.
        // This is robust to multiple bridge reads.
        {
            const auto cores = greedyDenseCoreSplit(adj, minAnchorCoverage, maxAnchorCoverage);
            if(cores.size() >= 2) {
                return cores;
            }
        }

        // As a last attempt, peel a dense core from the full vertex.
        // This can drop bridge reads even when we cannot confidently split into 2+ anchors.
        {
            vector<uint32_t> all;
            all.reserve(n);
            for(uint32_t i=0; i<n; ++i) {
                all.push_back(i);
            }
            auto core = peelQuasiCliqueCore(adj, all, minAnchorCoverage, maxAnchorCoverage, 4, 5);
            if(!core.empty() && core.size() < all.size()) {
                return {std::move(core)};
            }
        }

        // Fall back to keeping all nodes if it yields a valid anchor.
        if(n >= minAnchorCoverage && n <= maxAnchorCoverage) {
            vector<uint32_t> all;
            all.reserve(n);
            for(uint32_t i=0; i<n; ++i) {
                all.push_back(i);
            }
            return {all};
        }

        // Otherwise, give up.
        return {};
    }

    // Like splitVertexByOverlapSupport, but uses phasing information when available:
    // - If there are at least two non-trivial connected components in the CIS-only overlap graph,
    //   we treat them as the core clusters.
    // - Reads without CIS edges (typically no het sites / unphased) are attached to exactly one
    //   core cluster based on the number of CIS edges they have to that cluster.
    // - This prevents unphased reads from "gluing" two haplotype clusters together via triangles.
    vector<vector<uint32_t>> splitVertexByOverlapSupportWithPhasing(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        if(!hasAnyCisEdge) {
            return splitVertexByOverlapSupport(adjAll, minAnchorCoverage, maxAnchorCoverage);
        }

        const uint32_t n = uint32_t(adjAll.size());
        vector<uint8_t> none(n, 0);
        const auto cisComps = connectedComponents(adjCis, none);

        vector<vector<uint32_t>> core;
        core.reserve(cisComps.size());
        for(const auto& comp : cisComps) {
            if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                core.push_back(comp);
            }
        }

        // If we don't have enough CIS structure to split, fall back to topology-only splitting.
        if(core.size() < 2) {
            return splitVertexByOverlapSupport(adjAll, minAnchorCoverage, maxAnchorCoverage);
        }

        vector<int32_t> assigned(n, -1);
        for(uint32_t ci=0; ci<uint32_t(core.size()); ++ci) {
            for(const uint32_t u : core[ci]) {
                assigned[u] = int32_t(ci);
            }
        }

        // Attach remaining nodes to exactly one core component using CIS support counts.
        vector<uint32_t> counts(core.size());
        for(uint32_t u=0; u<n; ++u) {
            if(assigned[u] != -1) {
                continue;
            }
            std::fill(counts.begin(), counts.end(), 0);
            for(const uint32_t v : adjCis[u]) {
                const int32_t ci = assigned[v];
                if(ci >= 0) {
                    ++counts[uint32_t(ci)];
                }
            }
            uint32_t bestCi = std::numeric_limits<uint32_t>::max();
            uint32_t bestCount = 0;
            uint32_t secondCount = 0;
            for(uint32_t ci=0; ci<counts.size(); ++ci) {
                const uint32_t c = counts[ci];
                if(c > bestCount) {
                    secondCount = bestCount;
                    bestCount = c;
                    bestCi = ci;
                } else if(c == bestCount && c != 0) {
                    // Tie.
                    secondCount = bestCount;
                } else if(c > secondCount) {
                    secondCount = c;
                }
            }
            if(bestCount == 0) {
                continue; // no CIS support -> leave unassigned
            }
            if(secondCount == bestCount) {
                continue; // ambiguous -> leave unassigned
            }
            if(core[bestCi].size() >= maxAnchorCoverage) {
                continue;
            }
            assigned[u] = int32_t(bestCi);
            core[bestCi].push_back(u);
        }

        // Return the split only if we still have at least two valid groups.
        vector<vector<uint32_t>> kept;
        kept.reserve(core.size());
        for(auto& comp : core) {
            if(comp.size() >= minAnchorCoverage && comp.size() <= maxAnchorCoverage) {
                kept.push_back(comp);
            }
        }
        if(kept.size() >= 2) {
            // Peel within each returned group to drop weakly supported members.
            return peelGroups(adjAll, kept, minAnchorCoverage, maxAnchorCoverage);
        }

        return splitVertexByOverlapSupport(adjAll, minAnchorCoverage, maxAnchorCoverage);
    }

    // Split using a designated set of "core" nodes (typically non-contained reads).
    // We first split the induced subgraph on the core nodes, then attach non-core nodes
    // to exactly one core cluster based on overlap-support edges. Ambiguous non-core nodes
    // (ties or insufficient support) are dropped.
    vector<vector<uint32_t>> splitVertexByOverlapSupportWithCoreMask(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        const vector<uint8_t>& isCore,
        uint32_t coreMinSize,
        uint32_t attachMinSupport,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        const uint32_t n = uint32_t(adjAll.size());
        if(n == 0 || isCore.size() != n) {
            return splitVertexByOverlapSupportWithPhasing(
                adjAll, adjCis, hasAnyCisEdge, minAnchorCoverage, maxAnchorCoverage);
        }

        vector<uint32_t> coreNodes;
        coreNodes.reserve(n);
        for(uint32_t i=0; i<n; ++i) {
            if(isCore[i]) {
                coreNodes.push_back(i);
            }
        }

        // Self-tuning: core splitting can be useful even when the number of non-contained reads
        // is small. We only require at least 2 core nodes (and at least one non-core node),
        // and let the final minAnchorCoverage decide whether the split is kept.
        if(coreNodes.size() < std::max<uint32_t>(2, coreMinSize) || coreNodes.size() >= n) {
            return splitVertexByOverlapSupportWithPhasing(
                adjAll, adjCis, hasAnyCisEdge, minAnchorCoverage, maxAnchorCoverage);
        }

        vector<int32_t> oldToCore(n, -1);
        for(uint32_t i=0; i<uint32_t(coreNodes.size()); ++i) {
            oldToCore[coreNodes[i]] = int32_t(i);
        }

        const uint32_t cn = uint32_t(coreNodes.size());
        vector<vector<uint32_t>> coreAdjAll(cn);
        vector<vector<uint32_t>> coreAdjCis(cn);
        bool coreHasAnyCisEdge = false;

        for(uint32_t ii=0; ii<cn; ++ii) {
            const uint32_t u = coreNodes[ii];
            auto& outAll = coreAdjAll[ii];
            for(const uint32_t v : adjAll[u]) {
                const int32_t jj = oldToCore[v];
                if(jj >= 0) {
                    outAll.push_back(uint32_t(jj));
                }
            }

            if(hasAnyCisEdge) {
                auto& outCis = coreAdjCis[ii];
                for(const uint32_t v : adjCis[u]) {
                    const int32_t jj = oldToCore[v];
                    if(jj >= 0) {
                        outCis.push_back(uint32_t(jj));
                        coreHasAnyCisEdge = true;
                    }
                }
            }
        }

        for(uint32_t i=0; i<cn; ++i) {
            auto& nbr = coreAdjAll[i];
            std::sort(nbr.begin(), nbr.end());
            nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
            auto& nbrC = coreAdjCis[i];
            std::sort(nbrC.begin(), nbrC.end());
            nbrC.erase(std::unique(nbrC.begin(), nbrC.end()), nbrC.end());
        }

        // Split only the core graph. Use a very permissive minimum coverage (1) so that
        // 2+ small core clusters (even singletons) can be proposed; contained reads will
        // then be attached, and the final minAnchorCoverage filter will decide.
        const uint64_t minCoreCoverage = 1;
        auto coreGroups = splitVertexByOverlapSupportWithPhasing(
            coreAdjAll,
            coreAdjCis,
            coreHasAnyCisEdge,
            minCoreCoverage,
            maxAnchorCoverage);
        if(coreGroups.size() < 2) {
            return splitVertexByOverlapSupportWithPhasing(
                adjAll, adjCis, hasAnyCisEdge, minAnchorCoverage, maxAnchorCoverage);
        }

        // Map core groups back to original indices and build an assignment array.
        vector<int32_t> groupOfNode(n, -1);
        vector<vector<uint32_t>> groups;
        groups.reserve(coreGroups.size());
        for(uint32_t gi=0; gi<uint32_t(coreGroups.size()); ++gi) {
            vector<uint32_t> g;
            g.reserve(coreGroups[gi].size());
            for(const uint32_t coreIdx : coreGroups[gi]) {
                const uint32_t u = coreNodes[coreIdx];
                g.push_back(u);
                groupOfNode[u] = int32_t(gi);
            }
            groups.push_back(std::move(g));
        }

        // Attach non-core nodes to exactly one core cluster using overlap support.
        // Use adjAll to avoid relying on phasing for contained reads.
        vector<uint32_t> counts(groups.size());
        for(uint32_t u=0; u<n; ++u) {
            if(isCore[u]) {
                continue;
            }
            std::fill(counts.begin(), counts.end(), 0);
            for(const uint32_t v : adjAll[u]) {
                const int32_t gi = (v < n) ? groupOfNode[v] : -1;
                if(gi >= 0) {
                    ++counts[uint32_t(gi)];
                }
            }

            uint32_t degToCore = 0;
            for(uint32_t gi=0; gi<counts.size(); ++gi) {
                degToCore += counts[gi];
            }

            uint32_t bestGi = std::numeric_limits<uint32_t>::max();
            uint32_t bestCount = 0;
            uint32_t secondCount = 0;
            for(uint32_t gi=0; gi<counts.size(); ++gi) {
                const uint32_t c = counts[gi];
                if(c > bestCount) {
                    secondCount = bestCount;
                    bestCount = c;
                    bestGi = gi;
                } else if(c == bestCount && c != 0) {
                    secondCount = bestCount; // tie
                } else if(c > secondCount) {
                    secondCount = c;
                }
            }

            // Self-tuning attachment threshold: for high-degree bridge candidates, require more support
            // to avoid attaching broadly connected contained reads.
            uint32_t dynamicMinSupport = attachMinSupport;
            if(degToCore >= 4) {
                // ceil(degToCore/4), capped.
                dynamicMinSupport = std::max(dynamicMinSupport, uint32_t((degToCore + 3) / 4));
                dynamicMinSupport = std::min(dynamicMinSupport, uint32_t(3));
            }
            if(bestCount < dynamicMinSupport) {
                continue;
            }
            if(secondCount == bestCount) {
                continue; // ambiguous bridge: connects similarly to multiple clusters
            }
            if(bestGi == std::numeric_limits<uint32_t>::max()) {
                continue;
            }
            if(groups[bestGi].size() >= maxAnchorCoverage) {
                continue;
            }
            groupOfNode[u] = int32_t(bestGi);
            groups[bestGi].push_back(u);
        }

        // Keep only valid anchors and peel within each group.
        vector<vector<uint32_t>> kept;
        kept.reserve(groups.size());
        for(auto& g : groups) {
            std::sort(g.begin(), g.end());
            g.erase(std::unique(g.begin(), g.end()), g.end());
            if(g.size() >= minAnchorCoverage && g.size() <= maxAnchorCoverage) {
                kept.push_back(std::move(g));
            }
        }
        if(kept.size() >= 2) {
            return peelGroups(adjAll, kept, minAnchorCoverage, maxAnchorCoverage);
        }

        return splitVertexByOverlapSupportWithPhasing(
            adjAll, adjCis, hasAnyCisEdge, minAnchorCoverage, maxAnchorCoverage);
    }
}

#if DINARA_TESTING
namespace dinara::testing {
    vector<vector<uint32_t>> splitVertexByOverlapSupportForTesting(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        return splitVertexByOverlapSupportWithPhasing(
            adjAll, adjCis, hasAnyCisEdge, minAnchorCoverage, maxAnchorCoverage);
    }

    vector<vector<uint32_t>> mclClusterForTesting(
        const vector<vector<uint32_t>>& adj,
        double inflation,
        uint32_t maxIterations)
    {
        const auto und = symmetrizeAdj(adj);
        return mclClusterUndirected(und, inflation, maxIterations);
    }

    vector<vector<uint32_t>> splitVertexByOverlapSupportWithCoreMaskForTesting(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        const vector<uint8_t>& isCore,
        uint32_t coreMinSize,
        uint32_t attachMinSupport,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        return splitVertexByOverlapSupportWithCoreMask(
            adjAll,
            adjCis,
            hasAnyCisEdge,
            isCore,
            coreMinSize,
            attachMinSupport,
            minAnchorCoverage,
            maxAnchorCoverage);
    }

    vector<vector<uint32_t>> autoSplitVertexForTesting(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        const vector<uint8_t>& isCore,
        bool useNonContainedCores,
        uint32_t coreMinSize,
        uint32_t attachMinSupport,
        bool useMclSecondary,
        uint32_t mclMinVertexSize,
        double mclInflation,
        uint32_t mclMaxIterations,
        double suspiciousMaxDensity,
        double suspiciousMaxAverageClustering,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        const uint32_t n = uint32_t(adjAll.size());
        vector<vector<uint32_t>> groups;

        if(useNonContainedCores) {
            groups = splitVertexByOverlapSupportWithCoreMask(
                adjAll,
                adjCis,
                hasAnyCisEdge,
                isCore,
                coreMinSize,
                attachMinSupport,
                minAnchorCoverage,
                maxAnchorCoverage);
        } else {
            groups = splitVertexByOverlapSupportWithPhasing(
                adjAll,
                adjCis,
                hasAnyCisEdge,
                minAnchorCoverage,
                maxAnchorCoverage);
        }

        if(groups.size() < 2 && useMclSecondary && n >= mclMinVertexSize) {
            const auto& baseAdj = hasAnyCisEdge ? adjCis : adjAll;
            const auto undAdj = symmetrizeAdj(baseAdj);
            const double density = graphDensityUndirected(undAdj);
            const double avgClustering = averageLocalClusteringCoefficientUndirected(undAdj);
            if(density <= suspiciousMaxDensity && avgClustering <= suspiciousMaxAverageClustering) {
                auto clusters = mclClusterUndirected(undAdj, mclInflation, mclMaxIterations);
                vector<vector<uint32_t>> kept;
                kept.reserve(clusters.size());
                for(auto& c : clusters) {
                    if(c.size() >= minAnchorCoverage && c.size() <= maxAnchorCoverage) {
                        kept.push_back(std::move(c));
                    }
                }
                if(kept.size() >= 2) {
                    auto peeled = peelGroups(adjAll, kept, minAnchorCoverage, maxAnchorCoverage);
                    if(peeled.size() >= 2) {
                        groups = std::move(peeled);
                    }
                }
            }
        }

        return groups;
    }

    std::pair<vector<vector<uint32_t>>, bool> autoSplitVertexWithMclTriedFlagForTesting(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        const vector<uint8_t>& isCore,
        bool useNonContainedCores,
        uint32_t coreMinSize,
        uint32_t attachMinSupport,
        bool useMclSecondary,
        uint32_t mclMinVertexSize,
        double mclInflation,
        uint32_t mclMaxIterations,
        double suspiciousMaxDensity,
        double suspiciousMaxAverageClustering,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        const uint32_t n = uint32_t(adjAll.size());
        vector<vector<uint32_t>> groups;
        bool mclTried = false;

        if(useNonContainedCores) {
            groups = splitVertexByOverlapSupportWithCoreMask(
                adjAll,
                adjCis,
                hasAnyCisEdge,
                isCore,
                coreMinSize,
                attachMinSupport,
                minAnchorCoverage,
                maxAnchorCoverage);
        } else {
            groups = splitVertexByOverlapSupportWithPhasing(
                adjAll,
                adjCis,
                hasAnyCisEdge,
                minAnchorCoverage,
                maxAnchorCoverage);
        }

        if(groups.size() < 2 && useMclSecondary && n >= mclMinVertexSize) {
            const auto& baseAdj = hasAnyCisEdge ? adjCis : adjAll;
            const auto undAdj = symmetrizeAdj(baseAdj);
            const double density = graphDensityUndirected(undAdj);
            const double avgClustering = averageLocalClusteringCoefficientUndirected(undAdj);
            if(density <= suspiciousMaxDensity && avgClustering <= suspiciousMaxAverageClustering) {
                mclTried = true;
                auto clusters = mclClusterUndirected(undAdj, mclInflation, mclMaxIterations);
                vector<vector<uint32_t>> kept;
                kept.reserve(clusters.size());
                for(auto& c : clusters) {
                    if(c.size() >= minAnchorCoverage && c.size() <= maxAnchorCoverage) {
                        kept.push_back(std::move(c));
                    }
                }
                if(kept.size() >= 2) {
                    auto peeled = peelGroups(adjAll, kept, minAnchorCoverage, maxAnchorCoverage);
                    if(peeled.size() >= 2) {
                        groups = std::move(peeled);
                    }
                }
            }
        }

        return {groups, mclTried};
    }

    vector<vector<uint32_t>> splitVertexByCliqueCoverForTesting(
        const vector<vector<uint32_t>>& adjAll,
        const vector<vector<uint32_t>>& adjCis,
        bool hasAnyCisEdge,
        const vector<uint8_t>& isCore,
        uint32_t attachMinSupport,
        uint64_t minAnchorCoverage,
        uint64_t maxAnchorCoverage)
    {
        bool exploded = false;
        auto groups = splitVertexByCliqueCover(
            adjAll, adjCis, hasAnyCisEdge, isCore, attachMinSupport, minAnchorCoverage, maxAnchorCoverage, exploded);
        // For tests, treat "exploded" as no split.
        if(exploded) {
            return {};
        }
        return groups;
    }
}
#endif



shared_ptr<mode3::Anchors> Assembler::createAnchorsFromMarkerGraphVerticesAtOverlapEvents(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount)
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkReadGraphIsOpen();
    if(!alignmentData.isOpen) {
        throw runtime_error("Alignment data are not accessible.");
    }
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    const auto& mgVertices = markerGraph.vertices();

    vector<vector<MarkerGraphVertexId>> threadSelected(threadCount);

    uint64_t chunk = orientedReadCount / threadCount;
    if(chunk == 0) chunk = 1;

    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = (t == threadCount - 1) ? orientedReadCount : min(orientedReadCount, (t+1) * chunk);
            auto& out = threadSelected[t];
            out.reserve((end - begin) * 4);

            vector<uint32_t> eventOrdinals;
            eventOrdinals.reserve(256);

            for(uint64_t v=begin; v<end; ++v) {
                const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
                const uint32_t markerCount = uint32_t(markers->size(orientedReadId.getValue()));
                if(markerCount == 0) {
                    continue;
                }

                eventOrdinals.clear();

                // Collect start/end events from read-graph overlaps incident to this oriented read.
                for(const uint32_t edgeId : readGraph.connectivity[orientedReadId.getValue()]) {
                    const ReadGraphEdge& edge = readGraph.edges[edgeId];
                    if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                        continue;
                    }
                    const uint64_t alignmentId = edge.alignmentId;
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if(ad.isDeleted()) {
                        continue;
                    }
                    if(!ad.info.isInReadGraph) {
                        continue;
                    }
                    const OrientedReadId other = edge.getOther(orientedReadId);
                    const AlignmentInfo info = ad.orient(orientedReadId, other);

                    const uint32_t first = info.data[0].firstOrdinal;
                    const uint32_t last = info.data[0].lastOrdinal;
                    if(first < markerCount) {
                        eventOrdinals.push_back(first);
                    }
                    const uint32_t afterLast = last + 1;
                    if(afterLast < markerCount) {
                        eventOrdinals.push_back(afterLast);
                    }
                }

                if(eventOrdinals.empty()) {
                    continue;
                }

                std::sort(eventOrdinals.begin(), eventOrdinals.end());
                eventOrdinals.erase(std::unique(eventOrdinals.begin(), eventOrdinals.end()), eventOrdinals.end());

                MarkerGraphVertexId lastSelected = MarkerGraph::invalidVertexId;
                for(const uint32_t ordinal : eventOrdinals) {
                    if(ordinal >= markerCount) {
                        continue;
                    }
                    const MarkerId markerId = getMarkerId(orientedReadId, ordinal);
                    const auto compressedVertex = markerGraph.vertexTable[markerId];
                    if(compressedVertex == MarkerGraph::invalidCompressedVertexId) {
                        continue;
                    }
                    const MarkerGraphVertexId vertexId = asVertexId(compressedVertex);
                    if(vertexId >= mgVertices.size()) {
                        continue;
                    }
                    const uint64_t cov = mgVertices.size(vertexId);
                    if(cov < minAnchorCoverage || cov > maxAnchorCoverage) {
                        continue;
                    }

                    const MarkerGraphVertexId rc = markerGraph.reverseComplementVertex[vertexId];
                    const MarkerGraphVertexId canonical = min(vertexId, rc);
                    if(canonical == lastSelected) {
                        continue;
                    }
                    out.push_back(canonical);
                    lastSelected = canonical;
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    // Merge and deduplicate selected canonical vertices.
    vector<MarkerGraphVertexId> selected;
    {
        size_t total = 0;
        for(const auto& v : threadSelected) total += v.size();
        selected.reserve(total);
        for(auto& v : threadSelected) {
            selected.insert(selected.end(), v.begin(), v.end());
        }
    }
    std::sort(selected.begin(), selected.end());
    selected.erase(std::unique(selected.begin(), selected.end()), selected.end());

    cout << timestamp << "Selected " << selected.size()
         << " marker graph vertices from overlap events for anchors." << endl;

    return make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        markerGraph,
        selected,
        minAnchorCoverage,
        maxAnchorCoverage,
        threadCount);
}



shared_ptr<mode3::Anchors> Assembler::createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount,
    bool enableColinearityPeeling,
    double minDominantFractionToPeel)
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkReadGraphIsOpen();
    if(!alignmentData.isOpen) {
        throw runtime_error("Alignment data are not accessible.");
    }
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    const auto& mgVertices = markerGraph.vertices();

    struct ChainEntry {
        uint32_t ordinal = 0;
        MarkerGraphVertexId vertexId = MarkerGraph::invalidVertexId;
    };

    struct Event {
        uint32_t ordinal = 0;
        int8_t delta = 0;  // +1 start, -1 end (end = lastOrdinal+1)
        ReadId otherOrientedReadValue = invalidReadId;
    };

    // Optional: per oriented read, store the chosen canonical vertex chain as (ordinal, vertexId).
    // This is used to derive local colinearity/context signatures for vertex peeling/splitting.
    vector<vector<ChainEntry>> chosenVertexChains;
    if(enableColinearityPeeling) {
        chosenVertexChains.resize(orientedReadCount);
    }

    vector<vector<MarkerGraphVertexId>> threadSelected(threadCount);

    uint64_t chunk = orientedReadCount / threadCount;
    if(chunk == 0) chunk = 1;

    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = (t == threadCount - 1) ? orientedReadCount : min(orientedReadCount, (t+1) * chunk);
            auto& out = threadSelected[t];
            out.reserve((end - begin) * 2);

            vector<Event> events;
            events.reserve(512);

            std::unordered_map<MarkerGraphVertexId, bool> duplicateReadIdCache;
            duplicateReadIdCache.reserve(4096);

            // Track active overlap neighbors while sweeping event ordinals.
            std::unordered_map<ReadId, uint32_t> activeNeighborCount;
            activeNeighborCount.reserve(256);
            vector<pair<MarkerId, MarkerId>> neighborRanges;
            neighborRanges.reserve(256);

            // Cache per-segment "explained neighbor" counts for candidate vertices.
            std::unordered_map<MarkerGraphVertexId, uint32_t> explainedCache;
            explainedCache.reserve(128);

            for(uint64_t v=begin; v<end; ++v) {
                const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
                const uint32_t markerCount = uint32_t(markers->size(orientedReadId.getValue()));
                if(markerCount == 0) {
                    continue;
                }

                vector<ChainEntry>* chainPointer = nullptr;
                if(enableColinearityPeeling) {
                    chainPointer = &chosenVertexChains[v];
                    chainPointer->clear();
                }

                events.clear();

                // Collect start/end events from read-graph overlaps incident to this oriented read.
                for(const uint32_t edgeId : readGraph.connectivity[orientedReadId.getValue()]) {
                    const ReadGraphEdge& edge = readGraph.edges[edgeId];
                    if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                        continue;
                    }
                    const uint64_t alignmentId = edge.alignmentId;
                    const AlignmentData& ad = alignmentData[alignmentId];
                    // Prefer AND semantics: only use overlaps that both reads keep.
                    if(!ad.keptByBothSides()) {
                        continue;
                    }
                    // If cis/trans is populated, skip trans edges when defining segments.
                    if(ad.cisTransStatus == CisTransStatus::Trans) {
                        continue;
                    }
                    if(!ad.info.isInReadGraph) {
                        continue;
                    }

                    const OrientedReadId other = edge.getOther(orientedReadId);
                    const AlignmentInfo info = ad.orient(orientedReadId, other);

                    const uint32_t first = info.data[0].firstOrdinal;
                    const uint32_t last = info.data[0].lastOrdinal;
                    if(first < markerCount) {
                        events.push_back({first, +1, other.getValue()});
                    }
                    const uint32_t afterLast = last + 1;
                    if(afterLast < markerCount) {
                        events.push_back({afterLast, -1, other.getValue()});
                    }
                }

                if(events.empty()) {
                    continue;
                }

                std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
                    if(a.ordinal != b.ordinal) {
                        return a.ordinal < b.ordinal;
                    }
                    // Order does not matter for correctness (we apply all events at the same ordinal),
                    // but deterministic ordering helps reproducibility.
                    if(a.delta != b.delta) {
                        return a.delta > b.delta;
                    }
                    return a.otherOrientedReadValue < b.otherOrientedReadValue;
                });

                activeNeighborCount.clear();
                int32_t activeOverlapCount = 0;
                MarkerGraphVertexId lastSelected = MarkerGraph::invalidVertexId;

                // Sweep over distinct event ordinals.
                for(size_t i=0; i<events.size(); ) {
                    const uint32_t segmentStart = events[i].ordinal;

                    // Apply all events at this ordinal before processing the segment starting here.
                    while(i < events.size() && events[i].ordinal == segmentStart) {
                        const Event& e = events[i];
                        activeOverlapCount += e.delta;
                        const ReadId otherValue = e.otherOrientedReadValue;
                        if(otherValue != invalidReadId) {
                            if(e.delta > 0) {
                                ++activeNeighborCount[otherValue];
                            } else {
                                auto it = activeNeighborCount.find(otherValue);
                                if(it != activeNeighborCount.end()) {
                                    if(it->second > 1) {
                                        --it->second;
                                    } else {
                                        activeNeighborCount.erase(it);
                                    }
                                }
                            }
                        }
                        ++i;
                    }

                    const uint32_t segmentEnd = (i < events.size()) ? events[i].ordinal : markerCount;
                    if(activeOverlapCount <= 0) {
                        continue;
                    }
                    if(segmentEnd <= segmentStart) {
                        continue;
                    }

                    // Build marker-id ranges for currently active overlap neighbors.
                    // Each oriented read corresponds to a contiguous range of global MarkerIds.
                    neighborRanges.clear();
                    neighborRanges.reserve(activeNeighborCount.size());
                    for(const auto& p : activeNeighborCount) {
                        const ReadId otherValue = p.first;
                        if(otherValue == orientedReadId.getValue()) {
                            continue;
                        }
                        const OrientedReadId other = OrientedReadId::fromValue(otherValue);
                        const uint32_t otherMarkerCount = uint32_t(markers->size(otherValue));
                        if(otherMarkerCount == 0) {
                            continue;
                        }
                        const MarkerId beginMarkerId = getMarkerId(other, 0);
                        neighborRanges.push_back({beginMarkerId, beginMarkerId + otherMarkerCount});
                    }
                    std::sort(neighborRanges.begin(), neighborRanges.end(),
                        [](const auto& a, const auto& b) { return a.first < b.first; });

                    explainedCache.clear();

                    // Scan all ordinals in this segment and select the "best" marker graph vertex.
                    // Canonicalize by RC to avoid duplicates.
                    // In addition to coverage, score candidates by how well they explain the
                    // currently-active overlap neighborhood of this read (read-local context).
                    const uint32_t mid = segmentStart + (segmentEnd - segmentStart) / 2;

                    uint32_t bestAnyExplained = 0;
                    uint64_t bestAnyCov = 0;
                    MarkerGraphVertexId bestAny = MarkerGraph::invalidVertexId;
                    uint32_t bestAnyDistanceToMid = std::numeric_limits<uint32_t>::max();
                    uint32_t bestAnyOrdinal = std::numeric_limits<uint32_t>::max();

                    uint32_t bestCleanExplained = 0;
                    uint64_t bestCleanCov = 0;
                    MarkerGraphVertexId bestClean = MarkerGraph::invalidVertexId;
                    uint32_t bestCleanDistanceToMid = std::numeric_limits<uint32_t>::max();
                    uint32_t bestCleanOrdinal = std::numeric_limits<uint32_t>::max();

                    MarkerGraphVertexId previous = MarkerGraph::invalidVertexId;

                    auto isBetterCandidate =
                        [](uint32_t explained,
                           uint64_t cov,
                           uint32_t distanceToMid,
                           MarkerGraphVertexId canonical,
                           uint32_t bestExplained,
                           uint64_t bestCov,
                           uint32_t bestDistanceToMid,
                           MarkerGraphVertexId bestCanonical) -> bool
                        {
                            if(canonical == MarkerGraph::invalidVertexId) {
                                return false;
                            }
                            if(bestCanonical == MarkerGraph::invalidVertexId) {
                                return true;
                            }
                            if(explained != bestExplained) {
                                return explained > bestExplained;
                            }
                            if(cov != bestCov) {
                                return cov > bestCov;
                            }
                            if(distanceToMid != bestDistanceToMid) {
                                return distanceToMid < bestDistanceToMid;
                            }
                            return canonical < bestCanonical;
                        };

                    auto countExplainedNeighbors = [&](MarkerGraphVertexId vertexId) -> uint32_t {
                        auto it = explainedCache.find(vertexId);
                        if(it != explainedCache.end()) {
                            return it->second;
                        }
                        uint32_t explained = 0;
                        if(!neighborRanges.empty()) {
                            const span<const MarkerId> vertexMarkerIds = mgVertices[vertexId];
                            for(const auto& r : neighborRanges) {
                                const MarkerId beginId = r.first;
                                const MarkerId endId = r.second;
                                auto posIt = std::lower_bound(vertexMarkerIds.begin(), vertexMarkerIds.end(), beginId);
                                if(posIt != vertexMarkerIds.end() && *posIt < endId) {
                                    ++explained;
                                }
                            }
                        }
                        explainedCache.insert({vertexId, explained});
                        return explained;
                    };

                    for(uint32_t ordinal = segmentStart; ordinal < segmentEnd; ++ordinal) {
                        const MarkerId markerId = getMarkerId(orientedReadId, ordinal);
                        const auto compressedVertex = markerGraph.vertexTable[markerId];
                        if(compressedVertex == MarkerGraph::invalidCompressedVertexId) {
                            continue;
                        }
                        const MarkerGraphVertexId vertexId = asVertexId(compressedVertex);
                        if(vertexId >= mgVertices.size()) {
                            continue;
                        }
                        const uint64_t cov = mgVertices.size(vertexId);
                        if(cov < minAnchorCoverage || cov > maxAnchorCoverage) {
                            continue;
                        }

                        const MarkerGraphVertexId rc = markerGraph.reverseComplementVertex[vertexId];
                        const MarkerGraphVertexId canonical = min(vertexId, rc);
                        if(canonical == previous) {
                            continue;
                        }
                        previous = canonical;

                        const uint32_t distanceToMid = (ordinal > mid) ? (ordinal - mid) : (mid - ordinal);

                        const uint32_t explained = countExplainedNeighbors(vertexId);

                        // Update best-any (no cleanliness constraint).
                        if(isBetterCandidate(
                               explained, cov, distanceToMid, canonical,
                               bestAnyExplained, bestAnyCov, bestAnyDistanceToMid, bestAny)) {
                            bestAnyExplained = explained;
                            bestAnyCov = cov;
                            bestAny = canonical;
                            bestAnyDistanceToMid = distanceToMid;
                            bestAnyOrdinal = ordinal;
                        }

                        // Update best-clean only if this vertex has no duplicate ReadIds.
                        if(isBetterCandidate(
                               explained, cov, distanceToMid, canonical,
                               bestCleanExplained, bestCleanCov, bestCleanDistanceToMid, bestClean)) {
                            auto dupIt = duplicateReadIdCache.find(canonical);
                            bool hasDuplicateReadIds = false;
                            if(dupIt != duplicateReadIdCache.end()) {
                                hasDuplicateReadIds = dupIt->second;
                            } else {
                                hasDuplicateReadIds = markerGraph.vertexHasDuplicateReadIds(canonical, *markers);
                                duplicateReadIdCache.insert({canonical, hasDuplicateReadIds});
                            }
                            if(!hasDuplicateReadIds) {
                                bestCleanExplained = explained;
                                bestCleanCov = cov;
                                bestClean = canonical;
                                bestCleanDistanceToMid = distanceToMid;
                                bestCleanOrdinal = ordinal;
                                // Only safe to early-exit if we found a clean vertex at maximum allowed coverage.
                                if(bestCleanCov == maxAnchorCoverage) {
                                    break;
                                }
                            }
                        }
                    }

                    // Prefer the candidate that best matches the active overlap neighborhood.
                    // Use cleanliness (no duplicate ReadIds) as a tie-breaker, not as an absolute rule.
                    MarkerGraphVertexId chosen = bestAny;
                    if(bestClean != MarkerGraph::invalidVertexId) {
                        if(bestAny == MarkerGraph::invalidVertexId ||
                            bestCleanExplained > bestAnyExplained ||
                            (bestCleanExplained == bestAnyExplained)) {
                            chosen = bestClean;
                        }
                    }
                    if(chosen == MarkerGraph::invalidVertexId) {
                        continue;
                    }
                    if(chosen == lastSelected) {
                        continue;
                    }
                    out.push_back(chosen);
                    if(chainPointer) {
                        uint32_t chosenOrdinal = bestAnyOrdinal;
                        if(chosen == bestClean) {
                            chosenOrdinal = bestCleanOrdinal;
                        }
                        if(chosenOrdinal != std::numeric_limits<uint32_t>::max()) {
                            chainPointer->push_back(ChainEntry{chosenOrdinal, chosen});
                        }
                    }
                    lastSelected = chosen;
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    // Merge and deduplicate selected canonical vertices.
    vector<MarkerGraphVertexId> selected;
    {
        size_t total = 0;
        for(const auto& v : threadSelected) total += v.size();
        selected.reserve(total);
        for(auto& v : threadSelected) {
            selected.insert(selected.end(), v.begin(), v.end());
        }
    }
    std::sort(selected.begin(), selected.end());
    selected.erase(std::unique(selected.begin(), selected.end()), selected.end());

    cout << timestamp << "Selected " << selected.size()
         << " marker graph vertices (best per overlap interval) for anchors." << endl;

    if(!enableColinearityPeeling) {
        return make_shared<mode3::Anchors>(
            MappedMemoryOwner(*this),
            getReads(),
            assemblerInfo->k,
            *markers,
            markerGraph,
            selected,
            minAnchorCoverage,
            maxAnchorCoverage,
            threadCount);
    }

    if(minDominantFractionToPeel < 0.) {
        minDominantFractionToPeel = 0.;
    }
    if(minDominantFractionToPeel > 1.) {
        minDominantFractionToPeel = 1.;
    }

    struct CtxKey {
        MarkerGraphVertexId a = MarkerGraph::invalidVertexId;
        MarkerGraphVertexId b = MarkerGraph::invalidVertexId; // a <= b
    };
    struct CtxKeyHash {
        size_t operator()(const CtxKey& k) const
        {
            // Hash-combine (boost-like).
            const uint64_t x = k.a;
            const uint64_t y = k.b;
            uint64_t h = x + 0x9e3779b97f4a7c15ULL;
            h ^= y + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
            return size_t(h);
        }
    };
    struct CtxKeyEq {
        bool operator()(const CtxKey& x, const CtxKey& y) const
        {
            return x.a == y.a && x.b == y.b;
        }
    };

    constexpr double maxUnclassifiedFractionToAct = 0.10;
    constexpr double minTwoClusterExplainedFractionToSplit = 0.95;
    constexpr double minSecondFractionToSplit = 0.25;

    vector<vector<vector<Interval>>> threadAnchors(threadCount);
    vector<pair<uint64_t, uint64_t>> threadSplitPeelCounts(threadCount, {0, 0});

    const uint64_t selectedCount = uint64_t(selected.size());
    uint64_t selectedChunk = (threadCount == 0) ? selectedCount : (selectedCount / threadCount);
    if(selectedChunk == 0) {
        selectedChunk = 1;
    }

    threads.clear();
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; ++t) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * selectedChunk;
            const uint64_t end =
                (t == threadCount - 1) ? selectedCount : min(selectedCount, (t + 1) * selectedChunk);
            auto& outAnchors = threadAnchors[t];
            outAnchors.reserve((end - begin) * 2ULL);

            uint64_t splitCount = 0;
            uint64_t peelCount = 0;

            vector<Interval> vertexIntervals;
            vector<Interval> filteredA;
            vector<Interval> filteredB;
            vector<CtxKey> ctxKeys;
            std::unordered_map<CtxKey, uint32_t, CtxKeyHash, CtxKeyEq> ctxCounts;
            ctxCounts.reserve(64);

            for(uint64_t i=begin; i<end; ++i) {
                const MarkerGraphVertexId vertexId = selected[i];
                const auto vertexMarkerIds = mgVertices[vertexId];
                if(vertexMarkerIds.empty()) {
                    continue;
                }

                vertexIntervals.clear();
                vertexIntervals.reserve(vertexMarkerIds.size());
                for(const MarkerId markerId : vertexMarkerIds) {
                    OrientedReadId orientedReadId;
                    uint32_t ordinal0;
                    tie(orientedReadId, ordinal0) = dinara::findMarkerId(markerId, *markers);
                    vertexIntervals.push_back(Interval{orientedReadId, ordinal0});
                }

                std::sort(vertexIntervals.begin(), vertexIntervals.end(),
                    [](const Interval& a, const Interval& b) {
                        if(a.orientedReadId != b.orientedReadId) {
                            return a.orientedReadId < b.orientedReadId;
                        }
                        return a.ordinal0 < b.ordinal0;
                    });
                vertexIntervals.erase(
                    std::unique(vertexIntervals.begin(), vertexIntervals.end(),
                        [](const Interval& a, const Interval& b) {
                            return a.orientedReadId == b.orientedReadId;
                        }),
                    vertexIntervals.end());

                if(vertexIntervals.size() < minAnchorCoverage || vertexIntervals.size() > maxAnchorCoverage) {
                    continue;
                }

                bool didSplitOrPeel = false;

                // Compute per-read local context signature for this vertex using the per-read chain.
                ctxCounts.clear();
                ctxKeys.clear();
                ctxKeys.reserve(vertexIntervals.size());
                size_t classified = 0;
                size_t unclassified = 0;
                for(const Interval& interval : vertexIntervals) {
                    const ReadId orientedReadValue = interval.orientedReadId.getValue();
                    const auto& chain = chosenVertexChains[orientedReadValue];

                    CtxKey key;
                    if(!chain.empty()) {
                        const auto it = std::lower_bound(
                            chain.begin(), chain.end(), interval.ordinal0,
                            [](const ChainEntry& e, uint32_t ordinal) { return e.ordinal < ordinal; });
                        const size_t pos = size_t(it - chain.begin());

                        // If this vertex appears in the chain near this ordinal, use its immediate neighbors.
                        const bool matchAt = (pos < chain.size() && chain[pos].vertexId == vertexId);
                        const bool matchBefore = (pos > 0 && chain[pos - 1].vertexId == vertexId);

                        MarkerGraphVertexId prev = MarkerGraph::invalidVertexId;
                        MarkerGraphVertexId next = MarkerGraph::invalidVertexId;
                        if(matchAt && matchBefore) {
                            // Ambiguous (two occurrences straddling ordinal). Leave unclassified.
                        } else if(matchAt || matchBefore) {
                            const size_t occ = matchAt ? pos : (pos - 1);
                            if(occ > 0) {
                                prev = chain[occ - 1].vertexId;
                            }
                            if(occ + 1 < chain.size()) {
                                next = chain[occ + 1].vertexId;
                            }
                        } else {
                            // Vertex not present in the chain at this location: use the surrounding chain entries.
                            if(pos > 0) {
                                prev = chain[pos - 1].vertexId;
                            }
                            if(pos < chain.size()) {
                                next = chain[pos].vertexId;
                            }
                        }

                        if(prev != MarkerGraph::invalidVertexId || next != MarkerGraph::invalidVertexId) {
                            key.a = prev;
                            key.b = next;
                            if(key.a > key.b) {
                                std::swap(key.a, key.b);
                            }
                        }
                    }

                    ctxKeys.push_back(key);
                    if(key.a == MarkerGraph::invalidVertexId && key.b == MarkerGraph::invalidVertexId) {
                        ++unclassified;
                    } else {
                        ++classified;
                        ++ctxCounts[key];
                    }
                }

                if(!ctxCounts.empty() && classified >= size_t(minAnchorCoverage)) {
                    const double unclassifiedFraction = double(unclassified) / double(vertexIntervals.size());
                    if(unclassifiedFraction <= maxUnclassifiedFractionToAct) {
                        // Find the top two context signatures by read count.
                        CtxKey bestKey;
                        CtxKey secondKey;
                        uint32_t bestCount = 0;
                        uint32_t secondCount = 0;
                        for(const auto& p : ctxCounts) {
                            const uint32_t count = p.second;
                            if(count > bestCount) {
                                secondCount = bestCount;
                                secondKey = bestKey;
                                bestCount = count;
                                bestKey = p.first;
                            } else if(count > secondCount) {
                                secondCount = count;
                                secondKey = p.first;
                            }
                        }

                        if(bestCount > 0) {
                            const double bestFraction = double(bestCount) / double(classified);
                            const double secondFraction = double(secondCount) / double(classified);
                            const double twoClusterFraction = double(bestCount + secondCount) / double(classified);

                            // Split when extremely clear: two dominant signatures explain almost everything.
                            if(secondCount >= minAnchorCoverage &&
                                secondFraction >= minSecondFractionToSplit &&
                                twoClusterFraction >= minTwoClusterExplainedFractionToSplit) {

                                filteredA.clear();
                                filteredB.clear();
                                filteredA.reserve(bestCount);
                                filteredB.reserve(secondCount);
                                for(size_t j=0; j<vertexIntervals.size(); ++j) {
                                    const CtxKey& k = ctxKeys[j];
                                    if(k.a == bestKey.a && k.b == bestKey.b) {
                                        filteredA.push_back(vertexIntervals[j]);
                                    } else if(k.a == secondKey.a && k.b == secondKey.b) {
                                        filteredB.push_back(vertexIntervals[j]);
                                    }
                                }

                                if(filteredA.size() >= minAnchorCoverage && filteredB.size() >= minAnchorCoverage) {
                                    outAnchors.push_back(filteredA);
                                    outAnchors.push_back(reverseComplementAnchor(filteredA, *markers));
                                    outAnchors.push_back(filteredB);
                                    outAnchors.push_back(reverseComplementAnchor(filteredB, *markers));
                                    ++splitCount;
                                    didSplitOrPeel = true;
                                }
                            }

                            // If we didn't split, we can still peel to the dominant signature when clear.
                            if(!didSplitOrPeel &&
                                bestFraction >= minDominantFractionToPeel &&
                                bestCount < vertexIntervals.size()) {
                                filteredA.clear();
                                filteredA.reserve(bestCount);
                                for(size_t j=0; j<vertexIntervals.size(); ++j) {
                                    const CtxKey& k = ctxKeys[j];
                                    if(k.a == bestKey.a && k.b == bestKey.b) {
                                        filteredA.push_back(vertexIntervals[j]);
                                    }
                                }
                                if(filteredA.size() >= minAnchorCoverage) {
                                    outAnchors.push_back(filteredA);
                                    outAnchors.push_back(reverseComplementAnchor(filteredA, *markers));
                                    ++peelCount;
                                    didSplitOrPeel = true;
                                }
                            }
                        }
                    }
                }

                if(!didSplitOrPeel) {
                    outAnchors.push_back(vertexIntervals);
                    outAnchors.push_back(reverseComplementAnchor(vertexIntervals, *markers));
                }
            }

            threadSplitPeelCounts[t] = {splitCount, peelCount};
        });
    }

    for(auto& th : threads) {
        th.join();
    }

    uint64_t splitTotal = 0;
    uint64_t peelTotal = 0;
    size_t totalExplicit = 0;
    for(uint64_t t=0; t<threadCount; ++t) {
        splitTotal += threadSplitPeelCounts[t].first;
        peelTotal += threadSplitPeelCounts[t].second;
        totalExplicit += threadAnchors[t].size();
    }

    vector<vector<Interval>> anchorsExplicit;
    anchorsExplicit.reserve(totalExplicit);
    for(uint64_t t=0; t<threadCount; ++t) {
        for(auto& anchor : threadAnchors[t]) {
            anchorsExplicit.push_back(std::move(anchor));
        }
    }

    if(splitTotal > 0 || peelTotal > 0) {
        cout << timestamp << "Colinearity vertex peeling: split=" << splitTotal
             << ", peeled=" << peelTotal << endl;
    }

    return make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        anchorsExplicit,
        /*ordinalOffset*/ 0,
        threadCount);
}



shared_ptr<mode3::Anchors> Assembler::createAnchorsFromMarkerGraphVerticesBestPerOverlapIntervalDecomposed(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount)
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkReadGraphIsOpen();
    if(!alignmentData.isOpen) {
        throw runtime_error("Alignment data are not accessible.");
    }
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    const uint64_t orientedReadCount = 2 * readCount;
    const auto& mgVertices = markerGraph.vertices();

    struct Event {
        uint32_t ordinal = 0;
        int32_t delta = 0; // +1 start, -1 end (end = lastOrdinal+1)
    };

    // Phase 1: select canonical marker graph vertices (same logic as best-per-interval).
    vector<vector<MarkerGraphVertexId>> threadSelected(threadCount);
    uint64_t chunk = orientedReadCount / threadCount;
    if(chunk == 0) chunk = 1;

    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = (t == threadCount - 1) ? orientedReadCount : min(orientedReadCount, (t+1) * chunk);
            auto& out = threadSelected[t];
            out.reserve((end - begin) * 2);

            vector<Event> events;
            events.reserve(512);
            vector<Event> merged;
            merged.reserve(512);

            std::unordered_map<MarkerGraphVertexId, bool> duplicateReadIdCache;
            duplicateReadIdCache.reserve(4096);

            for(uint64_t v=begin; v<end; ++v) {
                const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(v));
                const uint32_t markerCount = uint32_t(markers->size(orientedReadId.getValue()));
                if(markerCount == 0) {
                    continue;
                }

                events.clear();
                for(const uint32_t edgeId : readGraph.connectivity[orientedReadId.getValue()]) {
                    const ReadGraphEdge& edge = readGraph.edges[edgeId];
                    if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                        continue;
                    }
                    const uint64_t alignmentId = edge.alignmentId;
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if(!ad.keptByBothSides()) {
                        continue;
                    }
                    if(!ad.info.isInReadGraph) {
                        continue;
                    }

                    const OrientedReadId other = edge.getOther(orientedReadId);
                    const AlignmentInfo info = ad.orient(orientedReadId, other);

                    const uint32_t first = info.data[0].firstOrdinal;
                    const uint32_t last = info.data[0].lastOrdinal;
                    if(first < markerCount) {
                        events.push_back({first, +1});
                    }
                    const uint32_t afterLast = last + 1;
                    if(afterLast < markerCount) {
                        events.push_back({afterLast, -1});
                    }
                }

                if(events.empty()) {
                    continue;
                }

                std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
                    return a.ordinal < b.ordinal;
                });

                merged.clear();
                for(const auto& e : events) {
                    if(merged.empty() || merged.back().ordinal != e.ordinal) {
                        merged.push_back(e);
                    } else {
                        merged.back().delta += e.delta;
                    }
                }

                int32_t active = 0;
                MarkerGraphVertexId lastSelected = MarkerGraph::invalidVertexId;

                for(size_t i=0; i<merged.size(); ++i) {
                    active += merged[i].delta;
                    const uint32_t segmentStart = merged[i].ordinal;
                    const uint32_t segmentEnd = (i + 1 < merged.size()) ? merged[i + 1].ordinal : markerCount;
                    if(active <= 0) {
                        continue;
                    }
                    if(segmentEnd <= segmentStart) {
                        continue;
                    }

                    const uint32_t mid = segmentStart + (segmentEnd - segmentStart) / 2;
                    uint64_t bestAnyCov = 0;
                    MarkerGraphVertexId bestAny = MarkerGraph::invalidVertexId;
                    uint32_t bestAnyDistanceToMid = std::numeric_limits<uint32_t>::max();

                    uint64_t bestCleanCov = 0;
                    MarkerGraphVertexId bestClean = MarkerGraph::invalidVertexId;
                    uint32_t bestCleanDistanceToMid = std::numeric_limits<uint32_t>::max();

                    MarkerGraphVertexId previous = MarkerGraph::invalidVertexId;
                    for(uint32_t ordinal = segmentStart; ordinal < segmentEnd; ++ordinal) {
                        const MarkerId markerId = getMarkerId(orientedReadId, ordinal);
                        const auto compressedVertex = markerGraph.vertexTable[markerId];
                        if(compressedVertex == MarkerGraph::invalidCompressedVertexId) {
                            continue;
                        }
                        const MarkerGraphVertexId vertexId = asVertexId(compressedVertex);
                        if(vertexId >= mgVertices.size()) {
                            continue;
                        }
                        const uint64_t cov = mgVertices.size(vertexId);
                        if(cov < minAnchorCoverage || cov > maxAnchorCoverage) {
                            continue;
                        }

                        const MarkerGraphVertexId rc = markerGraph.reverseComplementVertex[vertexId];
                        const MarkerGraphVertexId canonical = min(vertexId, rc);
                        if(canonical == previous) {
                            continue;
                        }
                        previous = canonical;

                        const uint32_t distanceToMid = (ordinal > mid) ? (ordinal - mid) : (mid - ordinal);
                        if(cov > bestAnyCov ||
                           (cov == bestAnyCov && distanceToMid < bestAnyDistanceToMid) ||
                           (cov == bestAnyCov && distanceToMid == bestAnyDistanceToMid && canonical < bestAny)) {
                            bestAnyCov = cov;
                            bestAny = canonical;
                            bestAnyDistanceToMid = distanceToMid;
                        }

                        if(cov > bestCleanCov ||
                           (cov == bestCleanCov && distanceToMid < bestCleanDistanceToMid) ||
                           (cov == bestCleanCov && distanceToMid == bestCleanDistanceToMid && canonical < bestClean)) {
                            auto dupIt = duplicateReadIdCache.find(canonical);
                            bool hasDuplicateReadIds = false;
                            if(dupIt != duplicateReadIdCache.end()) {
                                hasDuplicateReadIds = dupIt->second;
                            } else {
                                hasDuplicateReadIds = markerGraph.vertexHasDuplicateReadIds(canonical, *markers);
                                duplicateReadIdCache.insert({canonical, hasDuplicateReadIds});
                            }
                            if(!hasDuplicateReadIds) {
                                bestCleanCov = cov;
                                bestClean = canonical;
                                bestCleanDistanceToMid = distanceToMid;
                                if(bestCleanCov == maxAnchorCoverage) {
                                    break;
                                }
                            }
                        }
                    }

                    const MarkerGraphVertexId chosen =
                        (bestClean != MarkerGraph::invalidVertexId) ? bestClean : bestAny;
                    if(chosen == MarkerGraph::invalidVertexId) {
                        continue;
                    }
                    if(chosen == lastSelected) {
                        continue;
                    }
                    out.push_back(chosen);
                    lastSelected = chosen;
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    vector<MarkerGraphVertexId> selected;
    {
        size_t total = 0;
        for(const auto& v : threadSelected) total += v.size();
        selected.reserve(total);
        for(auto& v : threadSelected) {
            selected.insert(selected.end(), v.begin(), v.end());
        }
    }
    std::sort(selected.begin(), selected.end());
    selected.erase(std::unique(selected.begin(), selected.end()), selected.end());

    cout << timestamp << "Selected " << selected.size()
         << " marker graph vertices (best per overlap interval) for anchor decomposition." << endl;

    // Phase 2: decompose each selected vertex using overlap support among its oriented reads.
    vector<vector<vector<Interval>>> threadAnchors(threadCount);
    uint64_t vChunk = selected.size() / threadCount;
    if(vChunk == 0) vChunk = 1;

    threads.clear();
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * vChunk;
            const uint64_t end = (t == threadCount - 1) ? selected.size() : min(selected.size(), (t + 1) * vChunk);
            auto& outAnchors = threadAnchors[t];
            outAnchors.reserve((end - begin) * 2);

            std::unordered_map<OrientedReadId::Int, uint32_t> index;
            index.reserve(256);

            vector<vector<uint32_t>> adjAll;
            adjAll.reserve(256);
            vector<vector<uint32_t>> adjCis;
            adjCis.reserve(256);

            vector<Interval> vertexIntervals;
            vertexIntervals.reserve(256);

            for(uint64_t i=begin; i<end; ++i) {
                const MarkerGraphVertexId vertexId = selected[i];
                const MarkerGraphVertexId rcVertexId = markerGraph.reverseComplementVertex[vertexId];
                if(vertexId > rcVertexId) {
                    continue; // not canonical
                }

                const auto vertexMarkerIds = mgVertices[vertexId];
                if(vertexMarkerIds.size() < minAnchorCoverage) {
                    continue;
                }

                vertexIntervals.clear();
                // Keep at most one marker per oriented read.
                OrientedReadId::Int lastOrientedReadValue = std::numeric_limits<OrientedReadId::Int>::max();
                for(const MarkerId markerId : vertexMarkerIds) {
                    OrientedReadId orientedReadId;
                    uint32_t ordinal0;
                    tie(orientedReadId, ordinal0) = dinara::findMarkerId(markerId, *markers);
                    const OrientedReadId::Int orientedReadValue = orientedReadId.getValue();
                    if(orientedReadValue == lastOrientedReadValue) {
                        continue;
                    }
                    lastOrientedReadValue = orientedReadValue;
                    vertexIntervals.emplace_back(orientedReadId, ordinal0);
                }

                if(vertexIntervals.size() < minAnchorCoverage) {
                    continue;
                }

                // Build an overlap-support graph among the oriented reads in this vertex.
                index.clear();
                for(uint32_t u=0; u<uint32_t(vertexIntervals.size()); ++u) {
                    index.emplace(vertexIntervals[u].orientedReadId.getValue(), u);
                }

                const uint32_t n = uint32_t(vertexIntervals.size());
                adjAll.assign(n, {});
                adjCis.assign(n, {});
                bool hasAnyCisEdge = false;
                for(uint32_t u=0; u<n; ++u) {
                    const OrientedReadId orientedReadId = vertexIntervals[u].orientedReadId;
                    for(const uint32_t edgeId : readGraph.connectivity[orientedReadId.getValue()]) {
                        const ReadGraphEdge& edge = readGraph.edges[edgeId];
                        if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                            continue;
                        }
                        const uint64_t alignmentId = edge.alignmentId;
                        const AlignmentData& ad = alignmentData[alignmentId];
                        if(!ad.keptByBothSides()) {
                            continue;
                        }
                        if(!ad.info.isInReadGraph) {
                            continue;
                        }
                        if(ad.cisTransStatus == CisTransStatus::Trans) {
                            continue;
                        }
                        const OrientedReadId other = edge.getOther(orientedReadId);
                        const auto it = index.find(other.getValue());
                        if(it == index.end()) {
                            continue;
                        }
                        const uint32_t vIdx = it->second;
                        if(vIdx == u) {
                            continue;
                        }
                        const uint32_t ordinalU = vertexIntervals[u].ordinal0;
                        const uint32_t ordinalV = vertexIntervals[vIdx].ordinal0;
                        const AlignmentInfo info = ad.orient(orientedReadId, other);
                        if(ordinalU < info.data[0].firstOrdinal || ordinalU > info.data[0].lastOrdinal) {
                            continue;
                        }
                        if(ordinalV < info.data[1].firstOrdinal || ordinalV > info.data[1].lastOrdinal) {
                            continue;
                        }
                        adjAll[u].push_back(vIdx);
                        if(ad.cisTransStatus == CisTransStatus::Cis) {
                            adjCis[u].push_back(vIdx);
                            hasAnyCisEdge = true;
                        }
                    }
                }
                for(uint32_t u=0; u<n; ++u) {
                    auto& nbr = adjAll[u];
                    std::sort(nbr.begin(), nbr.end());
                    nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
                }
                for(uint32_t u=0; u<n; ++u) {
                    auto& nbr = adjCis[u];
                    std::sort(nbr.begin(), nbr.end());
                    nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
                }

                const auto groups = splitVertexByOverlapSupportWithPhasing(
                    adjAll,
                    adjCis,
                    hasAnyCisEdge,
                    minAnchorCoverage,
                    maxAnchorCoverage);
                for(const auto& group : groups) {
                    vector<Interval> anchor;
                    anchor.reserve(group.size());
                    for(const uint32_t idxInVertex : group) {
                        anchor.push_back(vertexIntervals[idxInVertex]);
                    }
                    std::sort(anchor.begin(), anchor.end(), [](const Interval& a, const Interval& b) {
                        return a.orientedReadId < b.orientedReadId;
                    });
                    if(anchor.size() < minAnchorCoverage || anchor.size() > maxAnchorCoverage) {
                        continue;
                    }
                    outAnchors.push_back(anchor);
                    outAnchors.push_back(reverseComplementAnchor(anchor, *markers));
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    vector<vector<Interval>> anchorsExplicit;
    {
        size_t total = 0;
        for(const auto& v : threadAnchors) total += v.size();
        anchorsExplicit.reserve(total);
        for(auto& v : threadAnchors) {
            anchorsExplicit.insert(anchorsExplicit.end(), v.begin(), v.end());
        }
    }

    cout << timestamp << "Constructed " << anchorsExplicit.size()
         << " explicit anchors (including reverse complements) after decomposition." << endl;

    return make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        anchorsExplicit,
        /*ordinalOffset*/ 0,
        threadCount);
}



shared_ptr<mode3::Anchors> Assembler::createAnchorsFromMarkerGraphVerticesSplitUsingReadGraph(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    const Mode3AssemblyOptions& mode3Options,
    uint64_t threadCount)
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkReadGraphIsOpen();
    if(!alignmentData.isOpen) {
        throw runtime_error("Alignment data are not accessible.");
    }
    checkMarkerGraphVerticesAreAvailable();
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const auto& mgVertices = markerGraph.vertices();
    const uint64_t vertexCount = mgVertices.size();

    vector<vector<vector<Interval>>> threadAnchors(threadCount);

    std::atomic<uint64_t> canonicalVertices{0};
    std::atomic<uint64_t> candidateVertices{0};
    std::atomic<uint64_t> splitVertices{0};
    std::atomic<uint64_t> cliqueTriedVertices{0};
    std::atomic<uint64_t> cliqueUsedVertices{0};
    std::atomic<uint64_t> cliqueSplitVertices{0};
    std::atomic<uint64_t> cliqueExplodedVertices{0};
    std::atomic<uint64_t> mclTriedVertices{0};
    std::atomic<uint64_t> mclSplitVertices{0};
    std::atomic<uint64_t> emittedAnchors{0};

    uint64_t chunk = vertexCount / threadCount;
    if(chunk == 0) chunk = 1;

    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = (t == threadCount - 1) ? vertexCount : min(vertexCount, (t+1) * chunk);
            auto& outAnchors = threadAnchors[t];
            outAnchors.reserve((end - begin) / 8);

            vector<Interval> vertexIntervals;
            vertexIntervals.reserve(64);

            vector<vector<uint32_t>> adjAll;
            vector<vector<uint32_t>> adjCis;
            unordered_map<uint32_t, uint32_t> index;
            index.reserve(128);

            for(MarkerGraphVertexId vertexId=begin; vertexId<end; ++vertexId) {
                const MarkerGraphVertexId rcVertexId = markerGraph.reverseComplementVertex[vertexId];
                if(vertexId > rcVertexId) {
                    continue;
                }
                ++canonicalVertices;

                const auto vertexMarkerIds = mgVertices[vertexId];
                if(vertexMarkerIds.size() < minAnchorCoverage) {
                    continue;
                }
                if(vertexMarkerIds.size() > maxAnchorCoverage) {
                    // Fast reject: if marker count already exceeds max anchor coverage,
                    // the deduplicated read coverage can still be <= max, but this saves work
                    // for very large vertices that will almost always be excluded.
                    // Keep a small slack to avoid rejecting legitimate cases.
                    if(vertexMarkerIds.size() > maxAnchorCoverage * 3) {
                        continue;
                    }
                }
                if(markerGraph.vertexHasDuplicateReadIds(vertexId, *markers)) {
                    continue;
                }

                vertexIntervals.clear();
                for(const MarkerId markerId : vertexMarkerIds) {
                    OrientedReadId orientedReadId;
                    uint32_t ordinal0;
                    tie(orientedReadId, ordinal0) = dinara::findMarkerId(markerId, *markers);
                    vertexIntervals.push_back(Interval{orientedReadId, ordinal0});
                }

                std::sort(vertexIntervals.begin(), vertexIntervals.end(),
                    [](const Interval& a, const Interval& b) {
                        if(a.orientedReadId != b.orientedReadId) {
                            return a.orientedReadId < b.orientedReadId;
                        }
                        return a.ordinal0 < b.ordinal0;
                    });
                vertexIntervals.erase(
                    std::unique(vertexIntervals.begin(), vertexIntervals.end(),
                        [](const Interval& a, const Interval& b) {
                            return a.orientedReadId == b.orientedReadId;
                        }),
                    vertexIntervals.end());

                const uint32_t n = uint32_t(vertexIntervals.size());
                if(n < minAnchorCoverage || n > maxAnchorCoverage) {
                    continue;
                }
                ++candidateVertices;

                index.clear();
                index.reserve(n * 2);
                for(uint32_t i=0; i<n; ++i) {
                    index[vertexIntervals[i].orientedReadId.getValue()] = i;
                }

                adjAll.assign(n, {});
                adjCis.assign(n, {});
                bool hasAnyCisEdge = false;

                for(uint32_t u=0; u<n; ++u) {
                    const OrientedReadId orientedReadId = vertexIntervals[u].orientedReadId;
                    for(const uint32_t edgeId : readGraph.connectivity[orientedReadId.getValue()]) {
                        const ReadGraphEdge& edge = readGraph.edges[edgeId];
                        if(edge.crossesStrands || edge.hasInconsistentAlignment) {
                            continue;
                        }
                        const uint64_t alignmentId = edge.alignmentId;
                        const AlignmentData& ad = alignmentData[alignmentId];
                        if(!ad.info.isInReadGraph) {
                            continue;
                        }
                        const OrientedReadId other = edge.getOther(orientedReadId);
                        const auto it = index.find(other.getValue());
                        if(it == index.end()) {
                            continue;
                        }
                        const uint32_t vIdx = it->second;
                        if(vIdx == u) {
                            continue;
                        }

                        // Only use overlaps that cover both marker ordinals.
                        const uint32_t ordinalU = vertexIntervals[u].ordinal0;
                        const uint32_t ordinalV = vertexIntervals[vIdx].ordinal0;
                        const AlignmentInfo info = ad.orient(orientedReadId, other);
                        if(ordinalU < info.data[0].firstOrdinal || ordinalU > info.data[0].lastOrdinal) {
                            continue;
                        }
                        if(ordinalV < info.data[1].firstOrdinal || ordinalV > info.data[1].lastOrdinal) {
                            continue;
                        }

                        adjAll[u].push_back(vIdx);
                        if(ad.cisTransStatus == CisTransStatus::Cis) {
                            adjCis[u].push_back(vIdx);
                            hasAnyCisEdge = true;
                        }
                    }
                }

                for(uint32_t u=0; u<n; ++u) {
                    auto& nbr = adjAll[u];
                    std::sort(nbr.begin(), nbr.end());
                    nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
                }
                for(uint32_t u=0; u<n; ++u) {
                    auto& nbr = adjCis[u];
                    std::sort(nbr.begin(), nbr.end());
                    nbr.erase(std::unique(nbr.begin(), nbr.end()), nbr.end());
                }

                auto groups = [&]() -> vector<vector<uint32_t>> {
                    const auto& opt = mode3Options.vertexSplitOptions;

                    // Core mask: non-contained reads are treated as core evidence by default.
                    vector<uint8_t> isCore(n, 1);
                    if(opt.useNonContainedCores) {
                        for(uint32_t u=0; u<n; ++u) {
                            const ReadId rid = vertexIntervals[u].orientedReadId.getReadId();
                            if(reads->getFlags(rid).isContained) {
                                isCore[u] = 0;
                            }
                        }
                    }

                    // First try a clique-cover splitter on the overlap-support graph.
                    ++cliqueTriedVertices;
                    bool exploded = false;
                    auto cliqueGroups = splitVertexByCliqueCover(
                        adjAll,
                        adjCis,
                        hasAnyCisEdge,
                        isCore,
                        opt.attachMinSupport,
                        minAnchorCoverage,
                        maxAnchorCoverage,
                        exploded);
                    if(exploded) {
                        ++cliqueExplodedVertices;
                    }
                    if(!cliqueGroups.empty()) {
                        ++cliqueUsedVertices;
                        if(cliqueGroups.size() >= 2) {
                            ++cliqueSplitVertices;
                        }
                        return cliqueGroups;
                    }

                    // Fall back to the existing robust splitter.
                    if(opt.useNonContainedCores) {
                        return splitVertexByOverlapSupportWithCoreMask(
                            adjAll,
                            adjCis,
                            hasAnyCisEdge,
                            isCore,
                            opt.coreMinSize,
                            opt.attachMinSupport,
                            minAnchorCoverage,
                            maxAnchorCoverage);
                    } else {
                        return splitVertexByOverlapSupportWithPhasing(
                            adjAll,
                            adjCis,
                            hasAnyCisEdge,
                            minAnchorCoverage,
                            maxAnchorCoverage);
                    }
                }();

                // Optional secondary splitter: Markov Clustering (MCL) for suspicious vertices that
                // remain a single cluster after the default bridge-removal + peeling logic.
                if(groups.size() < 2 && mode3Options.vertexSplitOptions.useMclSecondary) {
                    const auto& opt = mode3Options.vertexSplitOptions;
                    if(n >= opt.mclMinVertexSize) {
                        const auto& baseAdj = hasAnyCisEdge ? adjCis : adjAll;
                        const auto undAdj = symmetrizeAdj(baseAdj);
                        const double density = graphDensityUndirected(undAdj);
                        const double avgClustering = averageLocalClusteringCoefficientUndirected(undAdj);
                        if(density <= opt.suspiciousMaxDensity &&
                            avgClustering <= opt.suspiciousMaxAverageClustering) {

                            ++mclTriedVertices;
                            auto clusters = mclClusterUndirected(
                                undAdj,
                                opt.mclInflation,
                                opt.mclMaxIterations);

                            vector<vector<uint32_t>> kept;
                            kept.reserve(clusters.size());
                            for(auto& c : clusters) {
                                if(c.size() >= minAnchorCoverage && c.size() <= maxAnchorCoverage) {
                                    kept.push_back(std::move(c));
                                }
                            }

                            if(kept.size() >= 2) {
                                auto peeled = peelGroups(adjAll, kept, minAnchorCoverage, maxAnchorCoverage);
                                if(peeled.size() >= 2) {
                                    groups = std::move(peeled);
                                    ++mclSplitVertices;
                                }
                            }
                        }
                    }
                }

                if(groups.size() >= 2) {
                    ++splitVertices;
                }

                for(const auto& group : groups) {
                    vector<Interval> anchor;
                    anchor.reserve(group.size());
                    for(const uint32_t idxInVertex : group) {
                        anchor.push_back(vertexIntervals[idxInVertex]);
                    }
                    std::sort(anchor.begin(), anchor.end(), [](const Interval& a, const Interval& b) {
                        return a.orientedReadId < b.orientedReadId;
                    });
                    if(anchor.size() < minAnchorCoverage || anchor.size() > maxAnchorCoverage) {
                        continue;
                    }
                    outAnchors.push_back(anchor);
                    outAnchors.push_back(reverseComplementAnchor(anchor, *markers));
                    emittedAnchors += 2;
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    vector<vector<Interval>> anchorsExplicit;
    {
        size_t total = 0;
        for(const auto& v : threadAnchors) total += v.size();
        anchorsExplicit.reserve(total);
        for(auto& v : threadAnchors) {
            anchorsExplicit.insert(anchorsExplicit.end(), v.begin(), v.end());
        }
    }

    cout << timestamp << "Vertex-split anchors: canonicalVertices=" << canonicalVertices.load() <<
        ", candidates=" << candidateVertices.load() <<
        ", splitVertices=" << splitVertices.load() <<
        ", cliqueTriedVertices=" << cliqueTriedVertices.load() <<
        ", cliqueUsedVertices=" << cliqueUsedVertices.load() <<
        ", cliqueSplitVertices=" << cliqueSplitVertices.load() <<
        ", cliqueExplodedVertices=" << cliqueExplodedVertices.load() <<
        ", mclTriedVertices=" << mclTriedVertices.load() <<
        ", mclSplitVertices=" << mclSplitVertices.load() <<
        ", anchorsEmitted=" << anchorsExplicit.size() << " (including reverse complements)." << endl;

    return make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        anchorsExplicit,
        /*ordinalOffset*/ 0,
        threadCount);
}



shared_ptr<mode3::Anchors> Assembler::createAnchorsFromOverlapsBestPerOverlapInterval(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount)
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    computeMarkerKmerIds(threadCount);
    checkReadGraphIsOpen();
    if(!alignmentData.isOpen) {
        throw runtime_error("Alignment data are not accessible.");
    }

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    // Only generate anchors from strand-0 oriented reads (read+).
    // We always emit the reverse-complement anchor explicitly, so processing strand-1 reads would
    // create duplicate anchors and slow down anchor generation.
    const uint64_t orientedReadCount = readCount;
    const uint64_t k = assemblerInfo->k;
    // Same defaults as Shasta2 (--max-anchor-repeat-length).
    const vector<uint64_t> maxAnchorRepeatLength = {8, 3, 3, 3, 3};

    // Experimental: If we have informative het coverage (set by AssemblerHifiasmEC),
    // only generate anchors from reads that cover at least one informative het site,
    // and only using their cis overlaps (trans overlaps are excluded by readGraph filtering,
    // but we also explicitly skip those if cisTransStatus is populated).
    bool restrictToInformativeHetReads = false;
    vector<uint8_t> readHasInformativeHet(readCount, 0);
    {
        // We treat the dataset as "informative-het mode" if ANY alignment is flagged as
        // covering an informative het/SV site. This is independent of read graph construction
        // and avoids missing directional phasing cases that were filtered out by conservative AND.
        for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); ++alignmentId) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if(!ad.coversHetSite()) {
                continue;
            }
            restrictToInformativeHetReads = true;
            if(ad.readIds[0] < readCount) {
                readHasInformativeHet[ad.readIds[0]] = 1;
            }
            if(ad.readIds[1] < readCount) {
                readHasInformativeHet[ad.readIds[1]] = 1;
            }
        }
    }
    if(restrictToInformativeHetReads) {
        cout << timestamp << "Using informative-het-only mode for overlap-only anchors: "
             << "only reads covering an informative het site will contribute anchors, using cis overlaps only." << endl;
    } else {
        cout << timestamp << "No informative het coverage detected in the filtered read graph; "
             << "generating overlap-only anchors from all reads." << endl;
    }

    if((not markerKmerIds) or (not markerKmerIds->isOpen())) {
        throw runtime_error("MarkerKmerIds are required for FromOverlapsBestPerOverlapInterval (duplicate ReadId filter).");
    }

    // Precompute the set of canonical marker k-mers that have duplicate ReadIds in MarkerKmerIds
    // (meaning the same canonical k-mer appears more than once in at least one read).
    // These are problematic as anchor seeds because the seed k-mer is not unique within some read.
    vector<vector<KmerId>> duplicateCanonicalKmerIdsByThread(threadCount);
    {
        vector<thread> dupThreads;
        dupThreads.reserve(threadCount);

        uint64_t chunk2 = readCount / threadCount;
        if(chunk2 == 0) chunk2 = 1;

        for(uint64_t t=0; t<threadCount; t++) {
            dupThreads.emplace_back([&, t]() {
                const uint64_t begin = t * chunk2;
                const uint64_t end = (t == threadCount - 1) ? readCount : min(readCount, (t+1) * chunk2);

                auto& duplicates = duplicateCanonicalKmerIdsByThread[t];
                duplicates.clear();
                duplicates.reserve(1024);

                unordered_set<KmerId, KmerIdHasher> seen;
                seen.reserve(4096);

                for(uint64_t r=begin; r<end; ++r) {
                    const ReadId readId = ReadId(r);
                    const OrientedReadId or0(readId, 0);
                    const OrientedReadId or1(readId, 1);
                    const auto kmerIds0 = (*markerKmerIds)[or0.getValue()];
                    const auto kmerIds1 = (*markerKmerIds)[or1.getValue()];
                    const uint32_t markerCount = uint32_t(kmerIds0.size());
                    if(markerCount == 0) {
                        continue;
                    }
                    DINARA_ASSERT(kmerIds1.size() == markerCount);

                    seen.clear();
                    if(seen.bucket_count() < size_t(markerCount) * 2ULL) {
                        seen.reserve(size_t(markerCount) * 2ULL);
                    }

                    for(uint32_t ordinal0=0; ordinal0<markerCount; ++ordinal0) {
                        const uint32_t ordinal1 = markerCount - 1U - ordinal0;
                        const KmerId id0 = kmerIds0[ordinal0];
                        const KmerId id1 = kmerIds1[ordinal1];
                        const KmerId canonical = (id0 <= id1) ? id0 : id1;
                        const auto inserted = seen.insert(canonical).second;
                        if(not inserted) {
                            duplicates.push_back(canonical);
                        }
                    }
                }
            });
        }
        for(auto& th : dupThreads) {
            th.join();
        }
    }

    vector<KmerId> duplicateCanonicalKmerIds;
    {
        size_t total = 0;
        for(const auto& v : duplicateCanonicalKmerIdsByThread) total += v.size();
        duplicateCanonicalKmerIds.reserve(total);
        for(auto& v : duplicateCanonicalKmerIdsByThread) {
            duplicateCanonicalKmerIds.insert(duplicateCanonicalKmerIds.end(), v.begin(), v.end());
        }
        sort(duplicateCanonicalKmerIds.begin(), duplicateCanonicalKmerIds.end());
        duplicateCanonicalKmerIds.erase(
            unique(duplicateCanonicalKmerIds.begin(), duplicateCanonicalKmerIds.end()),
            duplicateCanonicalKmerIds.end());
    }

    unordered_set<KmerId, KmerIdHasher> canonicalKmerIdsWithDuplicateReadIds;
    canonicalKmerIdsWithDuplicateReadIds.reserve(duplicateCanonicalKmerIds.size() * 2ULL + 1ULL);
    for(const KmerId& id : duplicateCanonicalKmerIds) {
        canonicalKmerIdsWithDuplicateReadIds.insert(id);
    }
    cout << timestamp << "Identified " << canonicalKmerIdsWithDuplicateReadIds.size()
         << " canonical marker k-mers with duplicate ReadIds (will not be used as anchor seeds)." << endl;

    struct CandidateAnchor {
        ReadId seedReadId;
        uint32_t overlapIntervalIndex = 0;
        uint32_t intervalStart = 0;
        uint32_t intervalEnd = 0; // one past last
        uint32_t seedOrdinal = 0;
        uint32_t support = 0;
        vector<Interval> anchor;
    };

    // For each read (strand 0), store the overlap-event intervals (where the active overlap set is constant and non-empty).
    vector<vector<OverlapInterval>> overlapIntervalsPerRead(readCount);

    vector<vector<CandidateAnchor>> threadCandidates(threadCount);

    uint64_t chunk = orientedReadCount / threadCount;
    if(chunk == 0) chunk = 1;

    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = (t == threadCount - 1) ? orientedReadCount : min(orientedReadCount, (t+1) * chunk);
            auto& outCandidates = threadCandidates[t];
            outCandidates.reserve(end - begin);

            vector<OverlapEvent> events;
            events.reserve(512);
            vector<uint32_t> activeEdgeIds;
            activeEdgeIds.reserve(256);

            vector<uint32_t> segmentEdgeIds;
            segmentEdgeIds.reserve(256);

            for(uint64_t v=begin; v<end; ++v) {
                const OrientedReadId orientedReadId0(ReadId(v), 0);
                if(reads->getFlags(orientedReadId0.getReadId()).isContained) {
                    continue;
                }
                const uint32_t markerCount0 = uint32_t(markers->size(orientedReadId0.getValue()));
                if(markerCount0 == 0) {
                    continue;
                }
                const OrientedReadId orientedReadId0rc(ReadId(v), 1);
                const auto orientedReadKmerIds0 = (*markerKmerIds)[orientedReadId0.getValue()];
                const auto orientedReadKmerIds1 = (*markerKmerIds)[orientedReadId0rc.getValue()];
                DINARA_ASSERT(orientedReadKmerIds0.size() == markerCount0);
                DINARA_ASSERT(orientedReadKmerIds1.size() == markerCount0);

                auto& outIntervals = overlapIntervalsPerRead[v];
                outIntervals.clear();

                events.clear();
                if(restrictToInformativeHetReads) {
                    if(v >= readHasInformativeHet.size() || !readHasInformativeHet[v]) {
                        continue;
                    }
                }

                // Collect start/end events from filtered readGraph overlaps.
                for(const uint32_t edgeId : readGraph.connectivity[orientedReadId0.getValue()]) {
                    const ReadGraphEdge& edge = readGraph.edges[edgeId];
                    if(edge.crossesStrands) {
                        continue;
                    }
                    if(edge.hasInconsistentAlignment) {
                        continue;
                    }
                    const uint64_t alignmentId = edge.alignmentId;
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if(!ad.keptByBothSides()) {
                        continue;
                    }
                    if(!ad.info.isInReadGraph) {
                        continue;
                    }
                    if(restrictToInformativeHetReads) {
                        // Use directional phasing decisions:
                        // if this read marked the overlap as trans, it will have DeleteReasonPhase set
                        // from this read's perspective (even if the other read kept it).
                        const ReadId seedReadId = orientedReadId0.getReadId();
                        const AlignmentData::DeleteReasonMask reasons =
                            (ad.readIds[0] == seedReadId) ? ad.deleteReasons0 : ad.deleteReasons1;
                        if(reasons & AlignmentData::DeleteReasonPhase) {
                            continue;
                        }
                    }

                    const OrientedReadId orientedReadId1 = edge.getOther(orientedReadId0);
                    const AlignmentInfo info = ad.orient(orientedReadId0, orientedReadId1);
                    const uint32_t first = info.data[0].firstOrdinal;
                    const uint32_t last = info.data[0].lastOrdinal;
                    if(first < markerCount0) {
                        events.push_back({first, +1, edgeId});
                    }
                    const uint32_t afterLast = last + 1;
                    if(afterLast < markerCount0) {
                        events.push_back({afterLast, -1, edgeId});
                    }
                }

                if(events.empty()) {
                    continue;
                }

                std::sort(events.begin(), events.end(), [](const OverlapEvent& a, const OverlapEvent& b) {
                    if(a.ordinal != b.ordinal) {
                        return a.ordinal < b.ordinal;
                    }
                    // Deterministic: process end events before start events at the same ordinal.
                    return a.delta < b.delta;
                });

                auto deactivate = [&](uint32_t edgeIdToRemove) {
                    const auto it = std::find(activeEdgeIds.begin(), activeEdgeIds.end(), edgeIdToRemove);
                    if(it != activeEdgeIds.end()) {
                        *it = activeEdgeIds.back();
                        activeEdgeIds.pop_back();
                    }
                };

                auto shouldSkipKmerDueToRepeats = [&](const Kmer& kmer0) -> bool {
                    for(uint64_t i=0; i<maxAnchorRepeatLength.size(); i++) {
                        const uint64_t period = i + 1;
                        const uint64_t maxAllowedCopyNumber = maxAnchorRepeatLength[i];
                        uint64_t copies = 0;
                        switch(period) {
                        case 1: copies = kmer0.countExactRepeatCopies<1>(k); break;
                        case 2: copies = kmer0.countExactRepeatCopies<2>(k); break;
                        case 3: copies = kmer0.countExactRepeatCopies<3>(k); break;
                        case 4: copies = kmer0.countExactRepeatCopies<4>(k); break;
                        case 5: copies = kmer0.countExactRepeatCopies<5>(k); break;
                        case 6: copies = kmer0.countExactRepeatCopies<6>(k); break;
                        default:
                            // Shasta2 only supports up to period 6.
                            copies = 0;
                            break;
                        }
                        if(copies > maxAllowedCopyNumber) {
                            return true;
                        }
                    }
                    return false;
                };

                // Build an anchor for this seed marker by mapping the seed k-mer across all active overlaps.
                // Returns the full (unique-read) support count. The stored anchor is capped to maxAnchorCoverage
                // to keep memory bounded, but we still scan all overlaps to measure support accurately.
                auto buildAnchorAtSeed = [&](uint32_t seedOrdinal, const Kmer& seedKmer,
                    const vector<uint32_t>& edgeIds, vector<Interval>& anchorOut) -> uint32_t
                {
                    anchorOut.clear();
                    anchorOut.reserve(64);
                    anchorOut.emplace_back(orientedReadId0, seedOrdinal);

                    std::unordered_set<ReadId> usedReadIds;
                    usedReadIds.reserve(128);
                    usedReadIds.insert(orientedReadId0.getReadId());
                    uint32_t fullSupport = 1;
                    const uint32_t storeLimit = uint32_t(maxAnchorCoverage);

                    for(const uint32_t activeEdgeId : edgeIds) {
                        const ReadGraphEdge& edge = readGraph.edges[activeEdgeId];
                        if(edge.crossesStrands) {
                            continue;
                        }
                        if(edge.hasInconsistentAlignment) {
                            continue;
                        }
                        const uint64_t alignmentId = edge.alignmentId;
                        const AlignmentData& ad = alignmentData[alignmentId];
                        if(!ad.keptByBothSides()) {
                            continue;
                        }
                        if(!ad.info.isInReadGraph) {
                            continue;
                        }

                        const OrientedReadId orientedReadId1 = edge.getOther(orientedReadId0);
                        const ReadId readId1 = orientedReadId1.getReadId();
                        if(usedReadIds.contains(readId1)) {
                            continue;
                        }

                        const AlignmentInfo info = ad.orient(orientedReadId0, orientedReadId1);
                        uint32_t ordinal1 = 0;
                        if(!mapMarkerOrdinalByOffsetAndKmer(
                            *reads,
                            *markers,
                            k,
                            info,
                            orientedReadId0,
                            seedOrdinal,
                            orientedReadId1,
                            ordinal1,
                            seedKmer,
                            /*maxSearchRadius*/ 8,
                            /*maxOffsetRange*/ 32)) {
                            continue;
                        }

                        ++fullSupport;
                        if(anchorOut.size() < storeLimit) {
                            anchorOut.emplace_back(orientedReadId1, ordinal1);
                        }
                        usedReadIds.insert(readId1);
                    }
                    return fullSupport;
                };

                vector<Interval> tmpAnchor;
                vector<Interval> bestAnchor;

                activeEdgeIds.clear();
                size_t ei = 0;
                while(ei < events.size()) {
                    const uint32_t ordinal = events[ei].ordinal;

                    // Apply all events at this ordinal.
                    while(ei < events.size() && events[ei].ordinal == ordinal) {
                        if(events[ei].delta > 0) {
                            activeEdgeIds.push_back(events[ei].edgeId);
                        } else {
                            deactivate(events[ei].edgeId);
                        }
                        ++ei;
                    }

                    const uint32_t nextOrdinal = (ei < events.size()) ? events[ei].ordinal : markerCount0;
                    if(activeEdgeIds.empty()) {
                        continue;
                    }
                    if(nextOrdinal <= ordinal) {
                        continue;
                    }

                    const uint32_t segmentStart = ordinal;
                    const uint32_t segmentEnd = nextOrdinal;

                    // This overlap-event segment has a constant active overlap set.
                    // We further split long segments into shorter anchor intervals to avoid large gaps.
                    constexpr uint32_t maxIntervalMarkers = 200;
                    const uint32_t segmentLen = segmentEnd - segmentStart;
                    const bool splitSegment = (segmentLen <= maxIntervalMarkers);

                    // Early exit: even with perfect k-mer mapping, we cannot exceed 1 + (# unique reads in active edges).
                    segmentEdgeIds = activeEdgeIds;
                    std::sort(segmentEdgeIds.begin(), segmentEdgeIds.end(), [&](uint32_t a, uint32_t b) {
                        const ReadId ra = readGraph.edges[a].getOther(orientedReadId0).getReadId();
                        const ReadId rb = readGraph.edges[b].getOther(orientedReadId0).getReadId();
                        if(ra != rb) {
                            return ra < rb;
                        }
                        return a < b;
                    });
                    {
                        ReadId prev = invalid<ReadId>;
                        uint32_t uniqueOtherReads = 0;
                        for(const uint32_t edgeId : segmentEdgeIds) {
                            const ReadId r = readGraph.edges[edgeId].getOther(orientedReadId0).getReadId();
                            if(r != prev) {
                                ++uniqueOtherReads;
                                prev = r;
                            }
                        }
                        const uint32_t maxPossible = std::min<uint32_t>(uint32_t(maxAnchorCoverage), 1 + uniqueOtherReads);
                        if(maxPossible < minAnchorCoverage) {
                            continue;
                        }
                    }

                    // Select one anchor per interval. If the overlap-event segment is long, do not split it:
                    // pick the best seed across the entire segment.
                    for(uint32_t intervalStart = segmentStart; intervalStart < segmentEnd; intervalStart += (splitSegment ? maxIntervalMarkers : segmentLen)) {
                        const uint32_t intervalEnd = splitSegment ? std::min(segmentEnd, intervalStart + maxIntervalMarkers) : segmentEnd;
                        const uint32_t intervalLen = intervalEnd - intervalStart;
                        if(intervalLen == 0) {
                            continue;
                        }

                        outIntervals.push_back({intervalStart, intervalEnd});
                        const uint32_t overlapIntervalIndex = uint32_t(outIntervals.size() - 1);

                        // Pick the seed marker with maximum support by scanning all marker ordinals in the interval.
                        vector<uint32_t> seedCandidates;
                        seedCandidates.reserve(intervalLen);
                        for(uint32_t o = intervalStart; o < intervalEnd; ++o) {
                            seedCandidates.push_back(o);
                        }

                        uint32_t bestSeed = invalid<uint32_t>;
                        uint32_t bestScore = 0;        // capped to maxAnchorCoverage
                        uint32_t bestFullSupport = 0;  // full unique-read support (can exceed maxAnchorCoverage)
                        bestAnchor.clear();

                        auto seedBetter = [&](uint32_t a, uint32_t b) -> bool {
                            // Prefer the seed closer to the interval center (more robust), then lower ordinal for determinism.
                            const int64_t center = int64_t(intervalStart) + int64_t(intervalLen) / 2;
                            const int64_t da = (int64_t(a) >= center) ? (int64_t(a) - center) : (center - int64_t(a));
                            const int64_t db = (int64_t(b) >= center) ? (int64_t(b) - center) : (center - int64_t(b));
                            if(da != db) {
                                return da < db;
                            }
                            return a < b;
                        };

                        for(const uint32_t seedOrdinal : seedCandidates) {
                            if(seedOrdinal < intervalStart || seedOrdinal >= intervalEnd) {
                                continue;
                            }

                            // Skip this seed if its canonical marker k-mer has duplicate ReadIds globally
                            // (it appears more than once in at least one read).
                            const uint32_t seedOrdinalRc = markerCount0 - 1U - seedOrdinal;
                            const KmerId id0 = orientedReadKmerIds0[seedOrdinal];
                            const KmerId id1 = orientedReadKmerIds1[seedOrdinalRc];
                            const KmerId canonicalId = (id0 <= id1) ? id0 : id1;
                            if(canonicalKmerIdsWithDuplicateReadIds.contains(canonicalId)) {
                                continue;
                            }

                            const Kmer seedKmer = getMarkerKmer(*reads, *markers, k, orientedReadId0, seedOrdinal);
                            if(shouldSkipKmerDueToRepeats(seedKmer)) {
                                continue;
                            }
                            const uint32_t fullSupport = buildAnchorAtSeed(seedOrdinal, seedKmer, segmentEdgeIds, tmpAnchor);
                            if(fullSupport < minAnchorCoverage) {
                                continue;
                            }
                            const uint32_t score = std::min<uint32_t>(fullSupport, uint32_t(maxAnchorCoverage));

                            // Prefer higher score, then lower full support (avoids extremely high-multiplicity seeds
                            // when a max-coverage seed exists), then tie-break by proximity to interval center.
                            if(score > bestScore ||
                               (score == bestScore && (bestSeed == invalid<uint32_t> || fullSupport < bestFullSupport)) ||
                               (score == bestScore && fullSupport == bestFullSupport &&
                                   bestSeed != invalid<uint32_t> && seedBetter(seedOrdinal, bestSeed))) {
                                bestScore = score;
                                bestSeed = seedOrdinal;
                                bestFullSupport = fullSupport;
                                bestAnchor = tmpAnchor;
                            }
                        }

                        if(bestSeed == invalid<uint32_t> || bestScore < minAnchorCoverage) {
                            continue;
                        }

                        std::sort(bestAnchor.begin(), bestAnchor.end(), [](const Interval& a, const Interval& b) {
                            return a.orientedReadId < b.orientedReadId;
                        });

                        // Cap stored membership to maxAnchorCoverage to keep candidates bounded.
                        if(bestAnchor.size() > maxAnchorCoverage) {
                            bestAnchor.resize(maxAnchorCoverage);
                        }

                        CandidateAnchor candidate;
                        candidate.seedReadId = ReadId(v);
                        candidate.overlapIntervalIndex = overlapIntervalIndex;
                        candidate.intervalStart = intervalStart;
                        candidate.intervalEnd = intervalEnd;
                        candidate.seedOrdinal = bestSeed;
                        candidate.support = bestScore;
                        candidate.anchor = std::move(bestAnchor);
                        outCandidates.push_back(std::move(candidate));
                    }
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    vector<CandidateAnchor> candidates;
    {
        size_t total = 0;
        for(const auto& v : threadCandidates) total += v.size();
        candidates.reserve(total);
        for(auto& v : threadCandidates) {
            candidates.insert(candidates.end(), v.begin(), v.end());
        }
    }

    cout << timestamp << "Constructed " << candidates.size()
         << " candidate anchors from overlaps (one per overlap-event interval, with long-interval splitting)." << endl;

    // Diagnostics: candidate support distribution (support is capped to maxAnchorCoverage).
    {
        uint64_t count = 0;
        uint64_t sum = 0;
        uint32_t minSupport = std::numeric_limits<uint32_t>::max();
        uint32_t maxSupport = 0;
        vector<uint64_t> hist(uint64_t(maxAnchorCoverage) + 1, 0);
        for(const CandidateAnchor& c : candidates) {
            ++count;
            sum += c.support;
            minSupport = std::min(minSupport, c.support);
            maxSupport = std::max(maxSupport, c.support);
            if(c.support <= maxAnchorCoverage) {
                ++hist[c.support];
            }
        }
        if(count) {
            cout << timestamp << "[DIAG] Candidate anchor support (capped): "
                 << "min=" << minSupport << " max=" << maxSupport
                 << " mean=" << double(sum) / double(count)
                 << " (n=" << count << ")." << endl;
        }
    }

    // Select anchors with cross-read overlap-interval claiming:
    // if a read already has an anchor in a given overlap-event interval, skip selecting another anchor for that interval.
    // This reduces duplicates and spreads anchors across overlap-covered regions.
    vector<vector<uint8_t>> intervalClaimed(readCount);
    for(uint64_t v=0; v<readCount; ++v) {
        intervalClaimed[v].assign(overlapIntervalsPerRead[v].size(), 0);
    }

    std::sort(candidates.begin(), candidates.end(), [](const CandidateAnchor& a, const CandidateAnchor& b) {
        if(a.support != b.support) {
            return a.support > b.support; // stronger anchors first
        }
        if(a.seedReadId != b.seedReadId) {
            return a.seedReadId < b.seedReadId;
        }
        if(a.intervalStart != b.intervalStart) {
            return a.intervalStart < b.intervalStart;
        }
        if(a.intervalEnd != b.intervalEnd) {
            return a.intervalEnd < b.intervalEnd;
        }
        return a.seedOrdinal < b.seedOrdinal;
    });

    vector<vector<Interval>> selected;
    selected.reserve(candidates.size());

    auto claimIntervalForMarker = [&](const Interval& interval) -> std::optional<pair<ReadId, uint32_t>> {
        const ReadId readId = interval.orientedReadId.getReadId();
        if(uint64_t(readId) >= readCount) {
            return std::nullopt;
        }
        uint32_t idx = 0;
        if(!findOverlapIntervalIndex(overlapIntervalsPerRead[readId], interval.ordinal0, idx)) {
            return std::nullopt;
        }
        return pair<ReadId, uint32_t>(readId, idx);
    };

    auto selectCandidate = [&](const CandidateAnchor& candidate, bool allowClaimedInRepair) -> bool {
        (void)allowClaimedInRepair;
        const uint64_t v = uint64_t(candidate.seedReadId);
        if(v >= readCount) {
            return false;
        }
        if(candidate.overlapIntervalIndex >= overlapIntervalsPerRead[v].size()) {
            return false;
        }
        if(intervalClaimed[v][candidate.overlapIntervalIndex]) {
            return false;
        }

        vector<Interval> anchor;
        anchor.reserve(candidate.anchor.size());

        vector<pair<ReadId, uint32_t>> toClaim;
        toClaim.reserve(candidate.anchor.size());

        // Always include the seed marker, so we anchor this interval on the seed read.
        anchor.emplace_back(OrientedReadId(candidate.seedReadId, 0), candidate.seedOrdinal);
        toClaim.emplace_back(candidate.seedReadId, candidate.overlapIntervalIndex);

        // Add additional read intervals from this candidate, preferring unclaimed overlap-intervals for claiming.
        // Already-claimed intervals are allowed to contribute to the anchor (they improve coverage),
        // but they are not re-claimed. This prevents coverage collapse as claiming progresses.
        vector<pair<Interval, pair<ReadId, uint32_t>>> claimedDeferred;
        claimedDeferred.reserve(candidate.anchor.size());

        for(const Interval& interval : candidate.anchor) {
            const ReadId readId = interval.orientedReadId.getReadId();
            if(readId == candidate.seedReadId) {
                continue;
            }

            const auto keyOpt = claimIntervalForMarker(interval);
            if(!keyOpt) {
                continue;
            }
            const auto [rid, idx] = *keyOpt;
            if(intervalClaimed[rid][idx]) {
                claimedDeferred.push_back({interval, {rid, idx}});
                continue;
            }
            anchor.push_back(interval);
            toClaim.push_back({rid, idx});
            if(anchor.size() >= maxAnchorCoverage) {
                break;
            }
        }

        // Fill remaining capacity with already-claimed intervals (no new claims) to improve coverage.
        if(anchor.size() < maxAnchorCoverage) {
            for(const auto& x : claimedDeferred) {
                anchor.push_back(x.first);
                if(anchor.size() >= maxAnchorCoverage) {
                    break;
                }
            }
        }

        if(anchor.size() < minAnchorCoverage || anchor.size() > maxAnchorCoverage) {
            return false;
        }

        std::sort(anchor.begin(), anchor.end(), [](const Interval& a, const Interval& b) {
            return a.orientedReadId < b.orientedReadId;
        });

        selected.push_back(anchor);

        // Claim all unclaimed overlap-event intervals used by this anchor.
        for(const auto& [rid, idx] : toClaim) {
            intervalClaimed[rid][idx] = 1;
        }
        return true;
    };

    // Pass 1: strict (only uses unclaimed intervals on all reads).
    uint64_t selectedStrict = 0;
    for(const CandidateAnchor& candidate : candidates) {
        if(selectCandidate(candidate, /*allowClaimedInRepair*/ false)) {
            ++selectedStrict;
        }
    }

    // Pass 2 (repair) disabled for now.
    // It can increase coverage by allowing already-claimed intervals on other reads,
    // but it changes the strict "one anchor claims one interval" behavior.
    uint64_t selectedRepair = 0;

    cout << timestamp << "Selected " << selected.size()
         << " anchors (strict=" << selectedStrict << ", repair=" << selectedRepair << ")." << endl;

    // Diagnostics: selected anchor size distribution.
    {
        uint64_t count = 0;
        uint64_t sum = 0;
        uint64_t minSize = std::numeric_limits<uint64_t>::max();
        uint64_t maxSize = 0;
        vector<uint64_t> hist(uint64_t(maxAnchorCoverage) + 1, 0);
        for(const auto& a : selected) {
            const uint64_t n = a.size();
            ++count;
            sum += n;
            minSize = std::min(minSize, n);
            maxSize = std::max(maxSize, n);
            if(n <= maxAnchorCoverage) {
                ++hist[n];
            }
        }
        if(count) {
            cout << timestamp << "[DIAG] Selected anchor size: "
                 << "min=" << minSize << " max=" << maxSize
                 << " mean=" << double(sum) / double(count)
                 << " (n=" << count << ")." << endl;
        }
    }

    uint64_t totalIntervals = 0;
    uint64_t claimedIntervals = 0;
    for(uint64_t v=0; v<readCount; ++v) {
        totalIntervals += intervalClaimed[v].size();
        for(const uint8_t x : intervalClaimed[v]) {
            claimedIntervals += (x != 0);
        }
    }
    if(totalIntervals) {
        cout << timestamp << "Claimed " << claimedIntervals << " / " << totalIntervals
             << " overlap intervals (" << double(claimedIntervals) / double(totalIntervals) << ")." << endl;
    }

    // Emit explicit reverse complements.
    vector<vector<Interval>> anchorsExplicit;
    anchorsExplicit.reserve(selected.size() * 2);
    for(const auto& anchor : selected) {
        anchorsExplicit.push_back(anchor);
        anchorsExplicit.push_back(reverseComplementAnchor(anchor, *markers));
    }

    // Make output deterministic.
    auto anchorLessLex = [](const vector<Interval>& a, const vector<Interval>& b) -> bool {
        const size_t n = std::min(a.size(), b.size());
        for(size_t i=0; i<n; ++i) {
            if(a[i].orientedReadId != b[i].orientedReadId) {
                return a[i].orientedReadId < b[i].orientedReadId;
            }
            if(a[i].ordinal0 != b[i].ordinal0) {
                return a[i].ordinal0 < b[i].ordinal0;
            }
        }
        return a.size() < b.size();
    };
    std::sort(anchorsExplicit.begin(), anchorsExplicit.end(), [&](const vector<Interval>& a, const vector<Interval>& b) {
        if(a.size() != b.size()) {
            return a.size() > b.size();
        }
        return anchorLessLex(a, b);
    });
    anchorsExplicit.erase(std::unique(anchorsExplicit.begin(), anchorsExplicit.end(),
        [&](const vector<Interval>& a, const vector<Interval>& b) {
            return a.size() == b.size() && !anchorLessLex(a, b) && !anchorLessLex(b, a);
        }), anchorsExplicit.end());

    cout << timestamp << "Selected " << anchorsExplicit.size()
         << " anchors (including reverse complements) after deduplication." << endl;

    auto anchors = make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        anchorsExplicit,
        /*ordinalOffset*/ 0,
        threadCount);

    // These anchors are used for mode3-style exploration (and sometimes for downstream processing).
    // The local anchor graph and other navigation features require journeys.
    anchors->computeJourneys(threadCount);
    return anchors;
}



// ============================================================================
// createAnchorsFromOverlapsBestPerOverlapIntervalBidirectional
// ============================================================================
// BidirectionalReadGraph-aware variant of createAnchorsFromOverlapsBestPerOverlapInterval.
//
// Uses bidirectionalReadGraph (one vertex per physical read, one edge per alignment)
// instead of the strand-doubled readGraph.  Orientation-aware traversal via
// edge.traverse(readId, strand) replaces the legacy crossesStrands gate:
//
//   Old: readGraph.connectivity[orientedReadId.getValue()]
//        edge.crossesStrands → skip
//        edge.getOther(orientedReadId) → OrientedReadId
//
//   New: bidirectionalReadGraph.connectivity[readId]
//        edge.isDeleted → skip  (no crossesStrands concept)
//        edge.traverse(readId, 0) → (ReadId, Strand) → OrientedReadId
//
// This preserves cross-strand overlaps (inversions, segdups) that the legacy
// graph would discard, letting anchors span inversion breakpoints.
shared_ptr<mode3::Anchors> Assembler::createAnchorsFromOverlapsBestPerOverlapIntervalBidirectional(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount)
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    computeMarkerKmerIds(threadCount);
    checkBidirectionalReadGraphIsOpen();
    if(!alignmentData.isOpen) {
        throw runtime_error("Alignment data are not accessible.");
    }

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    const uint64_t readCount = reads->readCount();
    // Only generate anchors from strand-0 oriented reads (read+).
    // We always emit the reverse-complement anchor explicitly.
    const uint64_t orientedReadCount = readCount;
    const uint64_t k = assemblerInfo->k;
    const vector<uint64_t> maxAnchorRepeatLength = {8, 3, 3, 3, 3};

    // Informative het restriction (same logic as the strand-doubled version).
    bool restrictToInformativeHetReads = false;
    vector<uint8_t> readHasInformativeHet(readCount, 0);
    {
        for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); ++alignmentId) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if(!ad.coversHetSite()) {
                continue;
            }
            restrictToInformativeHetReads = true;
            if(ad.readIds[0] < readCount) {
                readHasInformativeHet[ad.readIds[0]] = 1;
            }
            if(ad.readIds[1] < readCount) {
                readHasInformativeHet[ad.readIds[1]] = 1;
            }
        }
    }
    if(restrictToInformativeHetReads) {
        cout << timestamp << "BRG: Using informative-het-only mode for overlap-only anchors." << endl;
    } else {
        cout << timestamp << "BRG: No informative het coverage detected; generating overlap-only anchors from all reads." << endl;
    }

    if((not markerKmerIds) or (not markerKmerIds->isOpen())) {
        throw runtime_error("MarkerKmerIds are required for FromOverlapsBestPerOverlapIntervalBidirectional.");
    }

    // Precompute canonical marker k-mers with duplicate ReadIds.
    vector<vector<KmerId>> duplicateCanonicalKmerIdsByThread(threadCount);
    {
        vector<thread> dupThreads;
        dupThreads.reserve(threadCount);

        uint64_t chunk2 = readCount / threadCount;
        if(chunk2 == 0) chunk2 = 1;

        for(uint64_t t=0; t<threadCount; t++) {
            dupThreads.emplace_back([&, t]() {
                const uint64_t begin = t * chunk2;
                const uint64_t end = (t == threadCount - 1) ? readCount : min(readCount, (t+1) * chunk2);

                auto& duplicates = duplicateCanonicalKmerIdsByThread[t];
                duplicates.clear();
                duplicates.reserve(1024);

                unordered_set<KmerId, KmerIdHasher> seen;
                seen.reserve(4096);

                for(uint64_t r=begin; r<end; ++r) {
                    const ReadId readId = ReadId(r);
                    const OrientedReadId or0(readId, 0);
                    const OrientedReadId or1(readId, 1);
                    const auto kmerIds0 = (*markerKmerIds)[or0.getValue()];
                    const auto kmerIds1 = (*markerKmerIds)[or1.getValue()];
                    const uint32_t markerCount = uint32_t(kmerIds0.size());
                    if(markerCount == 0) {
                        continue;
                    }
                    DINARA_ASSERT(kmerIds1.size() == markerCount);

                    seen.clear();
                    if(seen.bucket_count() < size_t(markerCount) * 2ULL) {
                        seen.reserve(size_t(markerCount) * 2ULL);
                    }

                    for(uint32_t ordinal0=0; ordinal0<markerCount; ++ordinal0) {
                        const uint32_t ordinal1 = markerCount - 1U - ordinal0;
                        const KmerId id0 = kmerIds0[ordinal0];
                        const KmerId id1 = kmerIds1[ordinal1];
                        const KmerId canonical = (id0 <= id1) ? id0 : id1;
                        const auto inserted = seen.insert(canonical).second;
                        if(not inserted) {
                            duplicates.push_back(canonical);
                        }
                    }
                }
            });
        }
        for(auto& th : dupThreads) {
            th.join();
        }
    }

    vector<KmerId> duplicateCanonicalKmerIds;
    {
        size_t total = 0;
        for(const auto& v : duplicateCanonicalKmerIdsByThread) total += v.size();
        duplicateCanonicalKmerIds.reserve(total);
        for(auto& v : duplicateCanonicalKmerIdsByThread) {
            duplicateCanonicalKmerIds.insert(duplicateCanonicalKmerIds.end(), v.begin(), v.end());
        }
        sort(duplicateCanonicalKmerIds.begin(), duplicateCanonicalKmerIds.end());
        duplicateCanonicalKmerIds.erase(
            unique(duplicateCanonicalKmerIds.begin(), duplicateCanonicalKmerIds.end()),
            duplicateCanonicalKmerIds.end());
    }

    unordered_set<KmerId, KmerIdHasher> canonicalKmerIdsWithDuplicateReadIds;
    canonicalKmerIdsWithDuplicateReadIds.reserve(duplicateCanonicalKmerIds.size() * 2ULL + 1ULL);
    for(const KmerId& id : duplicateCanonicalKmerIds) {
        canonicalKmerIdsWithDuplicateReadIds.insert(id);
    }
    cout << timestamp << "BRG: Identified " << canonicalKmerIdsWithDuplicateReadIds.size()
         << " canonical marker k-mers with duplicate ReadIds." << endl;

    struct CandidateAnchor {
        ReadId seedReadId;
        uint32_t overlapIntervalIndex = 0;
        uint32_t intervalStart = 0;
        uint32_t intervalEnd = 0;
        uint32_t seedOrdinal = 0;
        uint32_t support = 0;
        vector<Interval> anchor;
    };

    vector<vector<OverlapInterval>> overlapIntervalsPerRead(readCount);
    vector<vector<CandidateAnchor>> threadCandidates(threadCount);

    uint64_t chunk = orientedReadCount / threadCount;
    if(chunk == 0) chunk = 1;

    // ========================================================================
    // Main parallel loop: for each physical read (strand 0), sweep-line over
    // BRG edges to find overlap-event intervals and select anchor seeds.
    // ========================================================================
    vector<thread> threads;
    threads.reserve(threadCount);
    for(uint64_t t=0; t<threadCount; t++) {
        threads.emplace_back([&, t]() {
            const uint64_t begin = t * chunk;
            const uint64_t end = (t == threadCount - 1) ? orientedReadCount : min(orientedReadCount, (t+1) * chunk);
            auto& outCandidates = threadCandidates[t];
            outCandidates.reserve(end - begin);

            vector<OverlapEvent> events;
            events.reserve(512);
            vector<uint32_t> activeEdgeIds;
            activeEdgeIds.reserve(256);

            vector<uint32_t> segmentEdgeIds;
            segmentEdgeIds.reserve(256);

            for(uint64_t v=begin; v<end; ++v) {
                const ReadId readId0 = ReadId(v);
                const OrientedReadId orientedReadId0(readId0, 0);
                if(reads->getFlags(readId0).isContained) {
                    continue;
                }
                const uint32_t markerCount0 = uint32_t(markers->size(orientedReadId0.getValue()));
                if(markerCount0 == 0) {
                    continue;
                }
                const OrientedReadId orientedReadId0rc(readId0, 1);
                const auto orientedReadKmerIds0 = (*markerKmerIds)[orientedReadId0.getValue()];
                const auto orientedReadKmerIds1 = (*markerKmerIds)[orientedReadId0rc.getValue()];
                DINARA_ASSERT(orientedReadKmerIds0.size() == markerCount0);
                DINARA_ASSERT(orientedReadKmerIds1.size() == markerCount0);

                auto& outIntervals = overlapIntervalsPerRead[v];
                outIntervals.clear();

                events.clear();
                if(restrictToInformativeHetReads) {
                    if(v >= readHasInformativeHet.size() || !readHasInformativeHet[v]) {
                        continue;
                    }
                }

                // ----------------------------------------------------------
                // Collect start/end events from BidirectionalReadGraph edges.
                // ----------------------------------------------------------
                // BRG: connectivity is indexed by ReadId (not OrientedReadId::getValue()).
                // Each edge stores isSameStrand instead of crossesStrands.
                // We use edge.traverse(readId0, 0) to derive the neighbour's OrientedReadId.
                for(const uint32_t edgeId : bidirectionalReadGraph.connectivity[readId0]) {
                    const BidirectionalReadGraphEdge& edge = bidirectionalReadGraph.edges[edgeId];
                    if(edge.isDeleted) {
                        continue;
                    }
                    if(edge.hasInconsistentAlignment) {
                        continue;
                    }
                    const uint64_t alignmentId = edge.alignmentId;
                    const AlignmentData& ad = alignmentData[alignmentId];
                    if(!ad.keptByBothSides()) {
                        continue;
                    }
                    if(!ad.info.isInReadGraph) {
                        continue;
                    }
                    if(restrictToInformativeHetReads) {
                        const AlignmentData::DeleteReasonMask reasons =
                            (ad.readIds[0] == readId0) ? ad.deleteReasons0 : ad.deleteReasons1;
                        if(reasons & AlignmentData::DeleteReasonPhase) {
                            continue;
                        }
                    }

                    // Orientation-aware traversal: derive the neighbour's OrientedReadId.
                    const auto [readId1, strand1] = edge.traverse(readId0, Strand(0));
                    const OrientedReadId orientedReadId1(readId1, strand1);

                    const AlignmentInfo info = ad.orient(orientedReadId0, orientedReadId1);
                    const uint32_t first = info.data[0].firstOrdinal;
                    const uint32_t last = info.data[0].lastOrdinal;
                    if(first < markerCount0) {
                        events.push_back({first, +1, edgeId});
                    }
                    const uint32_t afterLast = last + 1;
                    if(afterLast < markerCount0) {
                        events.push_back({afterLast, -1, edgeId});
                    }
                }

                if(events.empty()) {
                    continue;
                }

                std::sort(events.begin(), events.end(), [](const OverlapEvent& a, const OverlapEvent& b) {
                    if(a.ordinal != b.ordinal) {
                        return a.ordinal < b.ordinal;
                    }
                    return a.delta < b.delta;
                });

                auto deactivate = [&](uint32_t edgeIdToRemove) {
                    const auto it = std::find(activeEdgeIds.begin(), activeEdgeIds.end(), edgeIdToRemove);
                    if(it != activeEdgeIds.end()) {
                        *it = activeEdgeIds.back();
                        activeEdgeIds.pop_back();
                    }
                };

                auto shouldSkipKmerDueToRepeats = [&](const Kmer& kmer0) -> bool {
                    for(uint64_t i=0; i<maxAnchorRepeatLength.size(); i++) {
                        const uint64_t period = i + 1;
                        const uint64_t maxAllowedCopyNumber = maxAnchorRepeatLength[i];
                        uint64_t copies = 0;
                        switch(period) {
                        case 1: copies = kmer0.countExactRepeatCopies<1>(k); break;
                        case 2: copies = kmer0.countExactRepeatCopies<2>(k); break;
                        case 3: copies = kmer0.countExactRepeatCopies<3>(k); break;
                        case 4: copies = kmer0.countExactRepeatCopies<4>(k); break;
                        case 5: copies = kmer0.countExactRepeatCopies<5>(k); break;
                        case 6: copies = kmer0.countExactRepeatCopies<6>(k); break;
                        default:
                            copies = 0;
                            break;
                        }
                        if(copies > maxAllowedCopyNumber) {
                            return true;
                        }
                    }
                    return false;
                };

                // Build an anchor for this seed marker using BRG edges.
                auto buildAnchorAtSeed = [&](uint32_t seedOrdinal, const Kmer& seedKmer,
                    const vector<uint32_t>& edgeIds, vector<Interval>& anchorOut) -> uint32_t
                {
                    anchorOut.clear();
                    anchorOut.reserve(64);
                    anchorOut.emplace_back(orientedReadId0, seedOrdinal);

                    std::unordered_set<ReadId> usedReadIds;
                    usedReadIds.reserve(128);
                    usedReadIds.insert(readId0);
                    uint32_t fullSupport = 1;
                    const uint32_t storeLimit = uint32_t(maxAnchorCoverage);

                    for(const uint32_t activeEdgeId : edgeIds) {
                        const BidirectionalReadGraphEdge& edge = bidirectionalReadGraph.edges[activeEdgeId];
                        if(edge.isDeleted) {
                            continue;
                        }
                        if(edge.hasInconsistentAlignment) {
                            continue;
                        }
                        const uint64_t alignmentId = edge.alignmentId;
                        const AlignmentData& ad = alignmentData[alignmentId];
                        if(!ad.keptByBothSides()) {
                            continue;
                        }
                        if(!ad.info.isInReadGraph) {
                            continue;
                        }

                        // BRG orientation-aware traversal.
                        const auto [readId1, strand1] = edge.traverse(readId0, Strand(0));
                        const OrientedReadId orientedReadId1(readId1, strand1);
                        if(usedReadIds.contains(readId1)) {
                            continue;
                        }

                        const AlignmentInfo info = ad.orient(orientedReadId0, orientedReadId1);
                        uint32_t ordinal1 = 0;
                        if(!mapMarkerOrdinalByOffsetAndKmer(
                            *reads,
                            *markers,
                            k,
                            info,
                            orientedReadId0,
                            seedOrdinal,
                            orientedReadId1,
                            ordinal1,
                            seedKmer,
                            /*maxSearchRadius*/ 8,
                            /*maxOffsetRange*/ 32)) {
                            continue;
                        }

                        ++fullSupport;
                        if(anchorOut.size() < storeLimit) {
                            anchorOut.emplace_back(orientedReadId1, ordinal1);
                        }
                        usedReadIds.insert(readId1);
                    }
                    return fullSupport;
                };

                vector<Interval> tmpAnchor;
                vector<Interval> bestAnchor;

                activeEdgeIds.clear();
                size_t ei = 0;
                while(ei < events.size()) {
                    const uint32_t ordinal = events[ei].ordinal;

                    while(ei < events.size() && events[ei].ordinal == ordinal) {
                        if(events[ei].delta > 0) {
                            activeEdgeIds.push_back(events[ei].edgeId);
                        } else {
                            deactivate(events[ei].edgeId);
                        }
                        ++ei;
                    }

                    const uint32_t nextOrdinal = (ei < events.size()) ? events[ei].ordinal : markerCount0;
                    if(activeEdgeIds.empty()) {
                        continue;
                    }
                    if(nextOrdinal <= ordinal) {
                        continue;
                    }

                    const uint32_t segmentStart = ordinal;
                    const uint32_t segmentEnd = nextOrdinal;

                    constexpr uint32_t maxIntervalMarkers = 200;
                    const uint32_t segmentLen = segmentEnd - segmentStart;
                    const bool splitSegment = (segmentLen <= maxIntervalMarkers);

                    // Early exit: count unique other reads from BRG edges.
                    segmentEdgeIds = activeEdgeIds;
                    std::sort(segmentEdgeIds.begin(), segmentEdgeIds.end(), [&](uint32_t a, uint32_t b) {
                        // BRG: getOther returns ReadId directly (no OrientedReadId unwrapping).
                        const ReadId ra = bidirectionalReadGraph.edges[a].getOther(readId0);
                        const ReadId rb = bidirectionalReadGraph.edges[b].getOther(readId0);
                        if(ra != rb) {
                            return ra < rb;
                        }
                        return a < b;
                    });
                    {
                        ReadId prev = invalid<ReadId>;
                        uint32_t uniqueOtherReads = 0;
                        for(const uint32_t edgeId : segmentEdgeIds) {
                            const ReadId r = bidirectionalReadGraph.edges[edgeId].getOther(readId0);
                            if(r != prev) {
                                ++uniqueOtherReads;
                                prev = r;
                            }
                        }
                        const uint32_t maxPossible = std::min<uint32_t>(uint32_t(maxAnchorCoverage), 1 + uniqueOtherReads);
                        if(maxPossible < minAnchorCoverage) {
                            continue;
                        }
                    }

                    for(uint32_t intervalStart = segmentStart; intervalStart < segmentEnd; intervalStart += (splitSegment ? maxIntervalMarkers : segmentLen)) {
                        const uint32_t intervalEnd = splitSegment ? std::min(segmentEnd, intervalStart + maxIntervalMarkers) : segmentEnd;
                        const uint32_t intervalLen = intervalEnd - intervalStart;
                        if(intervalLen == 0) {
                            continue;
                        }

                        outIntervals.push_back({intervalStart, intervalEnd});
                        const uint32_t overlapIntervalIndex = uint32_t(outIntervals.size() - 1);

                        vector<uint32_t> seedCandidates;
                        seedCandidates.reserve(intervalLen);
                        for(uint32_t o = intervalStart; o < intervalEnd; ++o) {
                            seedCandidates.push_back(o);
                        }

                        uint32_t bestSeed = invalid<uint32_t>;
                        uint32_t bestScore = 0;
                        uint32_t bestFullSupport = 0;
                        bestAnchor.clear();

                        auto seedBetter = [&](uint32_t a, uint32_t b) -> bool {
                            const int64_t center = int64_t(intervalStart) + int64_t(intervalLen) / 2;
                            const int64_t da = (int64_t(a) >= center) ? (int64_t(a) - center) : (center - int64_t(a));
                            const int64_t db = (int64_t(b) >= center) ? (int64_t(b) - center) : (center - int64_t(b));
                            if(da != db) {
                                return da < db;
                            }
                            return a < b;
                        };

                        for(const uint32_t seedOrdinal : seedCandidates) {
                            if(seedOrdinal < intervalStart || seedOrdinal >= intervalEnd) {
                                continue;
                            }

                            const uint32_t seedOrdinalRc = markerCount0 - 1U - seedOrdinal;
                            const KmerId id0 = orientedReadKmerIds0[seedOrdinal];
                            const KmerId id1 = orientedReadKmerIds1[seedOrdinalRc];
                            const KmerId canonicalId = (id0 <= id1) ? id0 : id1;
                            if(canonicalKmerIdsWithDuplicateReadIds.contains(canonicalId)) {
                                continue;
                            }

                            const Kmer seedKmer = getMarkerKmer(*reads, *markers, k, orientedReadId0, seedOrdinal);
                            if(shouldSkipKmerDueToRepeats(seedKmer)) {
                                continue;
                            }
                            const uint32_t fullSupport = buildAnchorAtSeed(seedOrdinal, seedKmer, segmentEdgeIds, tmpAnchor);
                            if(fullSupport < minAnchorCoverage) {
                                continue;
                            }
                            const uint32_t score = std::min<uint32_t>(fullSupport, uint32_t(maxAnchorCoverage));

                            if(score > bestScore ||
                               (score == bestScore && (bestSeed == invalid<uint32_t> || fullSupport < bestFullSupport)) ||
                               (score == bestScore && fullSupport == bestFullSupport &&
                                   bestSeed != invalid<uint32_t> && seedBetter(seedOrdinal, bestSeed))) {
                                bestScore = score;
                                bestSeed = seedOrdinal;
                                bestFullSupport = fullSupport;
                                bestAnchor = tmpAnchor;
                            }
                        }

                        if(bestSeed == invalid<uint32_t> || bestScore < minAnchorCoverage) {
                            continue;
                        }

                        std::sort(bestAnchor.begin(), bestAnchor.end(), [](const Interval& a, const Interval& b) {
                            return a.orientedReadId < b.orientedReadId;
                        });

                        if(bestAnchor.size() > maxAnchorCoverage) {
                            bestAnchor.resize(maxAnchorCoverage);
                        }

                        CandidateAnchor candidate;
                        candidate.seedReadId = readId0;
                        candidate.overlapIntervalIndex = overlapIntervalIndex;
                        candidate.intervalStart = intervalStart;
                        candidate.intervalEnd = intervalEnd;
                        candidate.seedOrdinal = bestSeed;
                        candidate.support = bestScore;
                        candidate.anchor = std::move(bestAnchor);
                        outCandidates.push_back(std::move(candidate));
                    }
                }
            }
        });
    }
    for(auto& th : threads) {
        th.join();
    }

    vector<CandidateAnchor> candidates;
    {
        size_t total = 0;
        for(const auto& v : threadCandidates) total += v.size();
        candidates.reserve(total);
        for(auto& v : threadCandidates) {
            candidates.insert(candidates.end(), v.begin(), v.end());
        }
    }

    cout << timestamp << "BRG: Constructed " << candidates.size()
         << " candidate anchors from overlaps." << endl;

    // Diagnostics.
    {
        uint64_t count = 0;
        uint64_t sum = 0;
        uint32_t minSupport = std::numeric_limits<uint32_t>::max();
        uint32_t maxSupport = 0;
        for(const CandidateAnchor& c : candidates) {
            ++count;
            sum += c.support;
            minSupport = std::min(minSupport, c.support);
            maxSupport = std::max(maxSupport, c.support);
        }
        if(count) {
            cout << timestamp << "[DIAG] BRG candidate anchor support: "
                 << "min=" << minSupport << " max=" << maxSupport
                 << " mean=" << double(sum) / double(count)
                 << " (n=" << count << ")." << endl;
        }
    }

    // Select anchors with cross-read overlap-interval claiming.
    vector<vector<uint8_t>> intervalClaimed(readCount);
    for(uint64_t v=0; v<readCount; ++v) {
        intervalClaimed[v].assign(overlapIntervalsPerRead[v].size(), 0);
    }

    std::sort(candidates.begin(), candidates.end(), [](const CandidateAnchor& a, const CandidateAnchor& b) {
        if(a.support != b.support) {
            return a.support > b.support;
        }
        if(a.seedReadId != b.seedReadId) {
            return a.seedReadId < b.seedReadId;
        }
        if(a.intervalStart != b.intervalStart) {
            return a.intervalStart < b.intervalStart;
        }
        if(a.intervalEnd != b.intervalEnd) {
            return a.intervalEnd < b.intervalEnd;
        }
        return a.seedOrdinal < b.seedOrdinal;
    });

    vector<vector<Interval>> selected;
    selected.reserve(candidates.size());

    auto claimIntervalForMarker = [&](const Interval& interval) -> std::optional<pair<ReadId, uint32_t>> {
        const ReadId readId = interval.orientedReadId.getReadId();
        if(uint64_t(readId) >= readCount) {
            return std::nullopt;
        }
        uint32_t idx = 0;
        if(!findOverlapIntervalIndex(overlapIntervalsPerRead[readId], interval.ordinal0, idx)) {
            return std::nullopt;
        }
        return pair<ReadId, uint32_t>(readId, idx);
    };

    auto selectCandidate = [&](const CandidateAnchor& candidate) -> bool {
        const uint64_t v = uint64_t(candidate.seedReadId);
        if(v >= readCount) {
            return false;
        }
        if(candidate.overlapIntervalIndex >= overlapIntervalsPerRead[v].size()) {
            return false;
        }
        if(intervalClaimed[v][candidate.overlapIntervalIndex]) {
            return false;
        }

        vector<Interval> anchor;
        anchor.reserve(candidate.anchor.size());

        vector<pair<ReadId, uint32_t>> toClaim;
        toClaim.reserve(candidate.anchor.size());

        anchor.emplace_back(OrientedReadId(candidate.seedReadId, 0), candidate.seedOrdinal);
        toClaim.emplace_back(candidate.seedReadId, candidate.overlapIntervalIndex);

        vector<pair<Interval, pair<ReadId, uint32_t>>> claimedDeferred;
        claimedDeferred.reserve(candidate.anchor.size());

        for(const Interval& interval : candidate.anchor) {
            const ReadId readId = interval.orientedReadId.getReadId();
            if(readId == candidate.seedReadId) {
                continue;
            }

            const auto keyOpt = claimIntervalForMarker(interval);
            if(!keyOpt) {
                continue;
            }
            const auto [rid, idx] = *keyOpt;
            if(intervalClaimed[rid][idx]) {
                claimedDeferred.push_back({interval, {rid, idx}});
                continue;
            }
            anchor.push_back(interval);
            toClaim.push_back({rid, idx});
            if(anchor.size() >= maxAnchorCoverage) {
                break;
            }
        }

        if(anchor.size() < maxAnchorCoverage) {
            for(const auto& x : claimedDeferred) {
                anchor.push_back(x.first);
                if(anchor.size() >= maxAnchorCoverage) {
                    break;
                }
            }
        }

        if(anchor.size() < minAnchorCoverage || anchor.size() > maxAnchorCoverage) {
            return false;
        }

        std::sort(anchor.begin(), anchor.end(), [](const Interval& a, const Interval& b) {
            return a.orientedReadId < b.orientedReadId;
        });

        selected.push_back(anchor);

        for(const auto& [rid, idx] : toClaim) {
            intervalClaimed[rid][idx] = 1;
        }
        return true;
    };

    uint64_t selectedCount = 0;
    for(const CandidateAnchor& candidate : candidates) {
        if(selectCandidate(candidate)) {
            ++selectedCount;
        }
    }

    cout << timestamp << "BRG: Selected " << selected.size() << " anchors." << endl;

    // Diagnostics: selected anchor size distribution.
    {
        uint64_t count = 0;
        uint64_t sum = 0;
        uint64_t minSize = std::numeric_limits<uint64_t>::max();
        uint64_t maxSize = 0;
        for(const auto& a : selected) {
            const uint64_t n = a.size();
            ++count;
            sum += n;
            minSize = std::min(minSize, n);
            maxSize = std::max(maxSize, n);
        }
        if(count) {
            cout << timestamp << "[DIAG] BRG selected anchor size: "
                 << "min=" << minSize << " max=" << maxSize
                 << " mean=" << double(sum) / double(count)
                 << " (n=" << count << ")." << endl;
        }
    }

    uint64_t totalIntervals = 0;
    uint64_t claimedIntervals = 0;
    for(uint64_t v=0; v<readCount; ++v) {
        totalIntervals += intervalClaimed[v].size();
        for(const uint8_t x : intervalClaimed[v]) {
            claimedIntervals += (x != 0);
        }
    }
    if(totalIntervals) {
        cout << timestamp << "BRG: Claimed " << claimedIntervals << " / " << totalIntervals
             << " overlap intervals (" << double(claimedIntervals) / double(totalIntervals) << ")." << endl;
    }

    // Emit explicit reverse complements.
    vector<vector<Interval>> anchorsExplicit;
    anchorsExplicit.reserve(selected.size() * 2);
    for(const auto& anchor : selected) {
        anchorsExplicit.push_back(anchor);
        anchorsExplicit.push_back(reverseComplementAnchor(anchor, *markers));
    }

    // Deterministic sort and dedup.
    auto anchorLessLex = [](const vector<Interval>& a, const vector<Interval>& b) -> bool {
        const size_t n = std::min(a.size(), b.size());
        for(size_t i=0; i<n; ++i) {
            if(a[i].orientedReadId != b[i].orientedReadId) {
                return a[i].orientedReadId < b[i].orientedReadId;
            }
            if(a[i].ordinal0 != b[i].ordinal0) {
                return a[i].ordinal0 < b[i].ordinal0;
            }
        }
        return a.size() < b.size();
    };
    std::sort(anchorsExplicit.begin(), anchorsExplicit.end(), [&](const vector<Interval>& a, const vector<Interval>& b) {
        if(a.size() != b.size()) {
            return a.size() > b.size();
        }
        return anchorLessLex(a, b);
    });
    anchorsExplicit.erase(std::unique(anchorsExplicit.begin(), anchorsExplicit.end(),
        [&](const vector<Interval>& a, const vector<Interval>& b) {
            return a.size() == b.size() && !anchorLessLex(a, b) && !anchorLessLex(b, a);
        }), anchorsExplicit.end());

    cout << timestamp << "BRG: Selected " << anchorsExplicit.size()
         << " anchors (including reverse complements) after deduplication." << endl;

    auto anchors = make_shared<mode3::Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        anchorsExplicit,
        /*ordinalOffset*/ 0,
        threadCount);

    anchors->computeJourneys(threadCount);
    return anchors;
}
