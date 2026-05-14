// Dinara.
#include "Assembler.hpp"
#include "Alignment.hpp"
#include "ProjectedAlignment.hpp"
#include "compressAlignment.hpp"
#include "DINARA_ASSERT.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"

// Standard library.
#include <algorithm>
#include <cstdint>
#include <deque>
#include <limits>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace dinara;
using namespace std;

namespace {

// Disjoint-set union over dynamically created mismatch nodes.
// Used to cluster transitive mismatch links across read overlaps.
class DynamicDisjointSets {
public:
    void reserve(uint64_t n) {
        parent.reserve(size_t(n));
        rank.reserve(size_t(n));
    }

    uint64_t add() {
        const uint64_t id = parent.size();
        parent.push_back(id);
        rank.push_back(0);
        return id;
    }

    uint64_t find(uint64_t x) {
        uint64_t root = x;
        while (parent[root] != root) {
            root = parent[root];
        }
        while (parent[x] != x) {
            const uint64_t next = parent[x];
            parent[x] = root;
            x = next;
        }
        return root;
    }

    void unite(uint64_t a, uint64_t b) {
        a = find(a);
        b = find(b);
        if (a == b) {
            return;
        }
        if (rank[a] < rank[b]) {
            std::swap(a, b);
        }
        parent[b] = a;
        if (rank[a] == rank[b]) {
            rank[a]++;
        }
    }

    uint64_t size() const {
        return parent.size();
    }

private:
    vector<uint64_t> parent;
    vector<uint8_t> rank;
};

// Key for a mismatch node identified by (readId, position).
inline uint64_t packReadPosKey(ReadId readId, uint32_t position)
{
    return (uint64_t(readId) << 32) | uint64_t(position);
}

// Key for propagated membership identified by (siteId, readId).
inline uint64_t packSiteReadKey(uint32_t siteId, ReadId readId)
{
    return (uint64_t(siteId) << 32) | uint64_t(readId);
}

// Decode delta-coded SNP evidence tokens into absolute positions.
// Hop tokens advance the running position but do not represent mismatches.
inline bool nextMismatchPosition(
    span<const SnpEvidence> tokens,
    size_t& tokenIndex,
    uint32_t& runningPos,
    uint32_t& mismatchPos)
{
    while (tokenIndex < tokens.size()) {
        const SnpEvidence ev = tokens[tokenIndex++];
        runningPos += ev.delta();
        if (ev.isHop()) {
            continue;
        }
        mismatchPos = runningPos;
        return true;
    }
    return false;
}

template<class F>
inline bool forEachMismatchPair(
    span<const SnpEvidence> qTokens,
    span<const SnpEvidence> tTokens,
    bool reverse,
    uint32_t reserveHint,
    vector<uint32_t>& reverseScratch,
    F&& f)
{
    if (!reverse) {
        size_t qTokenIndex = 0;
        size_t tTokenIndex = 0;
        uint32_t qRunningPos = 0;
        uint32_t tRunningPos = 0;
        uint32_t qMismatchPos = 0;
        uint32_t tMismatchPos = 0;
        while (true) {
            const bool qOk = nextMismatchPosition(qTokens, qTokenIndex, qRunningPos, qMismatchPos);
            const bool tOk = nextMismatchPosition(tTokens, tTokenIndex, tRunningPos, tMismatchPos);
            if (qOk != tOk) {
                return false;
            }
            if (!qOk) {
                return true;
            }
            f(qMismatchPos, tMismatchPos);
        }
    }

    reverseScratch.clear();
    if (reserveHint > reverseScratch.capacity()) {
        reverseScratch.reserve(reserveHint);
    }
    size_t tTokenIndex = 0;
    uint32_t tRunningPos = 0;
    uint32_t tMismatchPos = 0;
    while (nextMismatchPosition(tTokens, tTokenIndex, tRunningPos, tMismatchPos)) {
        reverseScratch.push_back(tMismatchPos);
    }

    size_t qTokenIndex = 0;
    uint32_t qRunningPos = 0;
    uint32_t qMismatchPos = 0;
    uint64_t mismatchIndex = 0;
    while (nextMismatchPosition(qTokens, qTokenIndex, qRunningPos, qMismatchPos)) {
        if (mismatchIndex >= reverseScratch.size()) {
            return false;
        }
        f(qMismatchPos, reverseScratch[reverseScratch.size() - 1 - mismatchIndex]);
        mismatchIndex++;
    }

    return mismatchIndex == reverseScratch.size();
}

class FlatNodeIdTable {
public:
    void reserve(uint64_t expectedSize)
    {
        // Keep load factor <= ~0.7.
        uint64_t capacity = 8;
        while (capacity * 7 < expectedSize * 10) {
            capacity *= 2;
        }
        rehash(capacity);
    }

    template<class CreateNodeId>
    uint64_t getOrCreate(uint64_t key, CreateNodeId&& createNodeId)
    {
        if (occupied.empty()) {
            rehash(8);
        }
        if ((used + 1) * 10 >= occupied.size() * 7) {
            rehash(occupied.size() * 2);
        }

        const uint64_t mask = occupied.size() - 1;
        uint64_t i = hash(key) & mask;
        while (true) {
            if (!occupied[i]) {
                occupied[i] = 1;
                keys[i] = key;
                values[i] = createNodeId();
                used++;
                return values[i];
            }
            if (keys[i] == key) {
                return values[i];
            }
            i = (i + 1) & mask;
        }
    }

private:
    vector<uint64_t> keys;
    vector<uint64_t> values;
    vector<uint8_t> occupied;
    uint64_t used = 0;

    static uint64_t hash(uint64_t x)
    {
        // splitmix64
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    void rehash(uint64_t newCapacity)
    {
        DINARA_ASSERT((newCapacity & (newCapacity - 1)) == 0);
        vector<uint64_t> oldKeys = std::move(keys);
        vector<uint64_t> oldValues = std::move(values);
        vector<uint8_t> oldOccupied = std::move(occupied);

        keys.assign(newCapacity, 0);
        values.assign(newCapacity, 0);
        occupied.assign(newCapacity, 0);
        used = 0;

        if (oldOccupied.empty()) {
            return;
        }

        const uint64_t mask = newCapacity - 1;
        for (size_t j = 0; j < oldOccupied.size(); j++) {
            if (!oldOccupied[j]) {
                continue;
            }
            const uint64_t key = oldKeys[j];
            uint64_t i = hash(key) & mask;
            while (occupied[i]) {
                i = (i + 1) & mask;
            }
            occupied[i] = 1;
            keys[i] = key;
            values[i] = oldValues[j];
            used++;
        }
    }
};

struct MapResult {
    bool ok = false;
    bool hole = false;
    uint32_t mappedPos = 0;
};

struct AlignmentAnchors {
    uint32_t first0 = 0;
    uint32_t first1 = 0;
    uint32_t last0 = 0;
    uint32_t last1 = 0;
    Strand strand1 = 0;
    bool ok = false;
};

// Shared per-thread cache for alignment anchor materialization.
// Reused across global-het traversal/verification/propagation paths.
struct ThreadLocalAlignmentCache {
    const Assembler* assembler = nullptr;
    vector<AlignmentAnchors> anchors;
    vector<uint8_t> anchorsComputed;
    vector<uint32_t> alignmentSeenEpoch;
    uint32_t epoch = 1;
};

inline AlignmentAnchors computeAlignmentAnchors(
    const Assembler& assembler,
    const AlignmentData& ad,
    uint32_t kHalf)
{
    AlignmentAnchors a;
    if (ad.info.markerCount == 0) {
        return a;
    }

    const OrientedReadId orientedReadId0(ad.readIds[0], 0);
    const OrientedReadId orientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
    a.strand1 = orientedReadId1.getStrand();

    const auto markers0 = (*assembler.markers)[orientedReadId0.getValue()];
    const auto markers1 = (*assembler.markers)[orientedReadId1.getValue()];
    if (markers0.empty() || markers1.empty()) {
        return a;
    }

    const uint32_t ord0First = ad.info.data[0].firstOrdinal;
    const uint32_t ord0Last = ad.info.data[0].lastOrdinal;
    const uint32_t ord1First = ad.info.data[1].firstOrdinal;
    const uint32_t ord1Last = ad.info.data[1].lastOrdinal;
    if (ord0First >= markers0.size() || ord0Last >= markers0.size() ||
        ord1First >= markers1.size() || ord1Last >= markers1.size()) {
        return a;
    }

    a.first0 = markers0[ord0First].position + kHalf;
    a.last0 = markers0[ord0Last].position + kHalf;

    const uint32_t pos1FirstOriented = markers1[ord1First].position + kHalf;
    const uint32_t pos1LastOriented = markers1[ord1Last].position + kHalf;
    if (a.strand1 == 0) {
        a.first1 = pos1FirstOriented;
        a.last1 = pos1LastOriented;
    } else {
        const uint32_t len1 = uint32_t(assembler.getReads().getRead(ad.readIds[1]).baseCount);
        if (pos1FirstOriented >= len1 || pos1LastOriented >= len1) {
            return a;
        }
        a.first1 = (len1 - 1U) - pos1FirstOriented; // Rightmost in forward coords.
        a.last1 = (len1 - 1U) - pos1LastOriented;   // Leftmost in forward coords.
    }

    a.ok = true;
    return a;
}

inline ThreadLocalAlignmentCache& getThreadLocalAlignmentCache(
    const Assembler& assembler,
    size_t alignmentCount)
{
    static thread_local ThreadLocalAlignmentCache cache;
    if (cache.assembler != &assembler || cache.anchors.size() != alignmentCount) {
        cache.assembler = &assembler;
        cache.anchors.assign(alignmentCount, AlignmentAnchors{});
        cache.anchorsComputed.assign(alignmentCount, 0);
        cache.alignmentSeenEpoch.assign(alignmentCount, 0);
        cache.epoch = 1;
    }
    return cache;
}

inline uint32_t beginAlignmentScanEpoch(ThreadLocalAlignmentCache& cache)
{
    cache.epoch++;
    if (cache.epoch == 0) {
        std::fill(cache.alignmentSeenEpoch.begin(), cache.alignmentSeenEpoch.end(), 0);
        cache.epoch = 1;
    }
    return cache.epoch;
}

inline const AlignmentAnchors& getCachedAlignmentAnchors(
    ThreadLocalAlignmentCache& cache,
    const Assembler& assembler,
    uint64_t alignmentId,
    const AlignmentData& ad,
    uint32_t kHalf)
{
    DINARA_ASSERT(alignmentId < cache.anchors.size());
    if (!cache.anchorsComputed[alignmentId]) {
        cache.anchors[alignmentId] = computeAlignmentAnchors(assembler, ad, kHalf);
        cache.anchorsComputed[alignmentId] = 1;
    }
    return cache.anchors[alignmentId];
}

// Map a position across an overlap using only sparse indels, in forward read coordinates.
// The anchor positions define the diagonal reference, and "dir" is +1 for same-strand overlaps
// and -1 for reverse-complement overlaps (forward coords run in opposite directions).
inline MapResult mapPositionUsingSparseIndels(
    uint32_t fromPos,
    uint32_t anchorFrom,
    uint32_t anchorTo,
    int dir,
    span<const IndelEvidence> indels)
{
    // Mapping is only valid within the region covered by the evidence origin.
    if (fromPos < anchorFrom) {
        return {};
    }

    int64_t mapped =
        int64_t(anchorTo) + int64_t(dir) * (int64_t(fromPos) - int64_t(anchorFrom));

    for (const IndelEvidence ev : indels) {
        const uint32_t p = ev.pos();
        if (p > fromPos) {
            break;
        }
        const uint32_t len = ev.len();
        if (ev.isInsertion()) {
            mapped += int64_t(dir) * int64_t(len);
        } else {
            const uint64_t end = uint64_t(p) + uint64_t(len);
            if (uint64_t(fromPos) < end) {
                return MapResult{false, true, 0};
            }
            mapped -= int64_t(dir) * int64_t(len);
        }
    }

    if (mapped < 0 || mapped > int64_t(std::numeric_limits<uint32_t>::max())) {
        return {};
    }
    return MapResult{true, false, uint32_t(mapped)};
}

class FlatAssignmentTable {
public:
    using Key = uint64_t;
    static constexpr Key emptyKey = std::numeric_limits<Key>::max();

    class ObserveResult {
    public:
        bool insertedKey = false;
        bool conflict = false;
    };

    uint64_t size() const {
        return used;
    }

    void reserve(uint64_t expectedSize) {
        // Keep load factor <= ~0.7.
        uint64_t capacity = 8;
        while (capacity * 7 < expectedSize * 10) {
            capacity *= 2;
        }
        rehash(capacity);
    }

    // Observe an assignment for (siteId, readId) -> position.
    // Keeps up to two competing position candidates with support counts.
    ObserveResult observe(Key key, uint32_t position)
    {
        if (keys.empty()) {
            rehash(8);
        }
        if ((used + 1) * 10 >= keys.size() * 7) {
            rehash(keys.size() * 2);
        }

        const uint64_t mask = keys.size() - 1;
        uint64_t i = hash(key) & mask;
        while (true) {
            if (keys[i] == emptyKey) {
                keys[i] = key;
                position0[i] = position;
                position1[i] = 0;
                support0[i] = 1;
                support1[i] = 0;
                hasPosition1[i] = 0;
                used++;
                return ObserveResult{true, false};
            }
            if (keys[i] == key) {
                if (position0[i] == position) {
                    support0[i]++;
                    return ObserveResult{};
                }
                if (hasPosition1[i] && position1[i] == position) {
                    support1[i]++;
                    return ObserveResult{};
                }

                ObserveResult out;
                out.conflict = true;
                if (!hasPosition1[i]) {
                    hasPosition1[i] = 1;
                    position1[i] = position;
                    support1[i] = 1;
                    return out;
                }

                // If we already have two candidates and observe a third one, retain two strong
                // candidates with a small-space majority heuristic.
                if (support0[i] <= support1[i]) {
                    if (support0[i] == 1) {
                        position0[i] = position;
                        support0[i] = 1;
                        // Do not enqueue propagation for third+ low-support churn.
                        // This replacement is only for final disambiguation bookkeeping.
                    } else {
                        support0[i]--;
                    }
                } else {
                    if (support1[i] == 1) {
                        position1[i] = position;
                        support1[i] = 1;
                        // Do not enqueue propagation for third+ low-support churn.
                        // This replacement is only for final disambiguation bookkeeping.
                    } else {
                        support1[i]--;
                    }
                }
                return out;
            }
            i = (i + 1) & mask;
        }
    }

    template<class F>
    void forEachOccupiedChosen(F&& f) const
    {
        for (size_t i = 0; i < keys.size(); i++) {
            if (keys[i] != emptyKey) {
                uint32_t chosenPos = position0[i];
                uint32_t otherPos = 0;
                uint32_t otherSupport = 0;
                if (hasPosition1[i]) {
                    if (support1[i] > support0[i] ||
                        (support1[i] == support0[i] && position1[i] < position0[i])) {
                        chosenPos = position1[i];
                        otherPos = position0[i];
                        otherSupport = support0[i];
                    } else {
                        otherPos = position1[i];
                        otherSupport = support1[i];
                    }
                }
                f(keys[i], chosenPos, otherPos, otherSupport);
            }
        }
    }

private:
    vector<Key> keys;
    vector<uint32_t> position0;
    vector<uint32_t> position1;
    vector<uint32_t> support0;
    vector<uint32_t> support1;
    vector<uint8_t> hasPosition1;
    uint64_t used = 0;

    static uint64_t hash(uint64_t x)
    {
        // splitmix64
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    void rehash(uint64_t newCapacity)
    {
        DINARA_ASSERT((newCapacity & (newCapacity - 1)) == 0);
        vector<Key> oldKeys = std::move(keys);
        vector<uint32_t> oldPosition0 = std::move(position0);
        vector<uint32_t> oldPosition1 = std::move(position1);
        vector<uint32_t> oldSupport0 = std::move(support0);
        vector<uint32_t> oldSupport1 = std::move(support1);
        vector<uint8_t> oldHasPosition1 = std::move(hasPosition1);

        keys.assign(newCapacity, emptyKey);
        position0.assign(newCapacity, 0);
        position1.assign(newCapacity, 0);
        support0.assign(newCapacity, 0);
        support1.assign(newCapacity, 0);
        hasPosition1.assign(newCapacity, 0);
        used = 0;

        if (!oldKeys.empty()) {
            const uint64_t mask = keys.size() - 1;
            for (size_t j = 0; j < oldKeys.size(); j++) {
                const Key key = oldKeys[j];
                if (key == emptyKey) {
                    continue;
                }
                uint64_t i = hash(key) & mask;
                while (keys[i] != emptyKey) {
                    i = (i + 1) & mask;
                }
                keys[i] = key;
                position0[i] = oldPosition0[j];
                position1[i] = oldPosition1[j];
                support0[i] = oldSupport0[j];
                support1[i] = oldSupport1[j];
                hasPosition1[i] = oldHasPosition1[j];
                used++;
            }
        }
    }
};

inline void materializeMismatchClustersFromComponents(
    DynamicDisjointSets& dsets,
    Assembler::GlobalMismatchSiteClusters& result,
    const Reads& readsRef)
{
    if (dsets.size() == 0) {
        result.clusterMemberOffsets = {0};
        return;
    }

    // Group nodes by DSU representative.
    vector<pair<uint64_t, uint64_t>> rootAndNode;
    rootAndNode.reserve(result.nodes.size());
    for (uint64_t nodeId = 0; nodeId < result.nodes.size(); nodeId++) {
        rootAndNode.push_back({dsets.find(nodeId), nodeId});
    }
    sort(rootAndNode.begin(), rootAndNode.end());

    result.clusterMemberOffsets.clear();
    result.clusterRepresentatives.clear();
    result.clusterMembers.clear();
    result.alleleCounts.clear();
    result.alleleCountsAreTransitive.clear();
    result.transitiveSiteMemberCounts.clear();
    result.clusterMemberOffsets.push_back(0);

    size_t i = 0;
    while (i < rootAndNode.size()) {
        const uint64_t root = rootAndNode[i].first;
        size_t j = i + 1;
        while (j < rootAndNode.size() && rootAndNode[j].first == root) {
            j++;
        }

        // Keep only non-singleton components (singletons are not mismatch-linked sites).
        if (j - i >= 2) {
            result.clusterRepresentatives.push_back(root);
            array<uint32_t, 4> counts{0, 0, 0, 0};
            for (size_t k = i; k < j; k++) {
                const uint64_t nodeId = rootAndNode[k].second;
                result.clusterMembers.push_back(nodeId);

                const auto& node = result.nodes[nodeId];
                const uint8_t base =
                    readsRef.getOrientedReadBase(OrientedReadId(node.first, 0), node.second).value;
                DINARA_ASSERT(base < 4);
                counts[base]++;
            }
            result.alleleCounts.push_back(counts);
            result.alleleCountsAreTransitive.push_back(0);
            result.transitiveSiteMemberCounts.push_back(uint32_t(j - i)); // mismatch-only count
            result.clusterMemberOffsets.push_back(uint64_t(result.clusterMembers.size()));
        }
        i = j;
    }
}

} // namespace


// Build global mismatch-site clusters across all accepted overlaps.
// Mismatch nodes are (readId, position), and connected components become global sites.
Assembler::GlobalMismatchSiteClusters Assembler::clusterMismatchingPositionsIntoGlobalHetSites(
    const AlignOptions& alignOptions,
    uint64_t threadCount,
    bool includeDeletedAlignments,
    bool readGraphOnly
) const
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    (void)threadCount; // Global-all clustering path remains single-threaded.

    // Typical datasets have several mismatch nodes per overlap. Reserve aggressively
    // to avoid repeated rehash/reallocation in the hot clustering loop.
    const uint64_t nodeReserveHint = alignmentData.size() * 6ULL + 1024ULL;

    DynamicDisjointSets dsets;
    dsets.reserve(nodeReserveHint);
    FlatNodeIdTable nodeIdByKey;
    nodeIdByKey.reserve(nodeReserveHint);

    GlobalMismatchSiteClusters result;
    result.nodes.reserve(size_t(nodeReserveHint));

    uint64_t totalMismatchEdges = 0;
    uint64_t skippedMismatchedCounts = 0;

    auto getOrAddNode = [&](ReadId readId, uint32_t pos) -> uint64_t {
        const uint64_t key = packReadPosKey(readId, pos);
        return nodeIdByKey.getOrCreate(key, [&]() {
            const uint64_t nodeId = dsets.add();
            DINARA_ASSERT(nodeId == result.nodes.size());
            result.nodes.push_back({readId, pos});
            return nodeId;
        });
    };

    // Fast path: use the already-stored mismatch evidence streams.
    // These store the same mismatches in both read coordinates (but without explicit pairing),
    // and the order is consistent with the underlying sparse mismatch list, so we can pair
    // by mismatch index (reversing the target order if the overlap is reverse-complemented).
    if (alignedEvidenceStore.index.size() >= alignmentData.size()) {
        vector<uint32_t> reverseScratch;
        reverseScratch.reserve(1024);

        for (uint64_t alignmentId = 0; alignmentId < alignmentData.size(); alignmentId++) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if (!includeDeletedAlignments && ad.isDeleted()) {
                continue;
            }
            if (readGraphOnly && !ad.info.isInReadGraph) {
                continue;
            }
            if (ad.info.errorRate > float(alignOptions.maxErrorRate)) {
                continue;
            }
            if (ad.info.mismatchCount == 0) {
                continue;
            }

            const uint32_t evidenceId = uint32_t(ad.info.alignmentId);
            const span<const SnpEvidence> qTokens = alignedEvidenceStore.getSnps1(evidenceId);
            const span<const SnpEvidence> tTokens = alignedEvidenceStore.getSnps0(evidenceId);
            const ReadId r0 = ad.readIds[0];
            const ReadId r1 = ad.readIds[1];
            uint64_t mismatchPairs = 0;
            const bool ok = forEachMismatchPair(
                qTokens,
                tTokens,
                !ad.isSameStrand,
                ad.info.mismatchCount,
                reverseScratch,
                [&](uint32_t pos0, uint32_t pos1) {
                    const uint64_t n0 = getOrAddNode(r0, pos0);
                    const uint64_t n1 = getOrAddNode(r1, pos1);
                    dsets.unite(n0, n1);
                    mismatchPairs++;
                });
            if (!ok) {
                skippedMismatchedCounts++;
                continue;
            }
            totalMismatchEdges += mismatchPairs;
        }

    } else {
        // Slow fallback: recompute sparse mismatches from stored marker-space alignments.
        if (!compressedAlignments.isOpen()) {
            throw runtime_error("clusterMismatchingPositionsIntoGlobalHetSites requires compressedAlignments to be open if evidence is unavailable.");
        }

        Alignment alignment;
        for (uint64_t alignmentId = 0; alignmentId < alignmentData.size(); alignmentId++) {
            const AlignmentData& ad = alignmentData[alignmentId];
            if (!includeDeletedAlignments && ad.isDeleted()) {
                continue;
            }
            if (readGraphOnly && !ad.info.isInReadGraph) {
                continue;
            }

            if (ad.info.errorRate > float(alignOptions.maxErrorRate)) {
                continue;
            }

            const span<const char> compressedAlignment = compressedAlignments[alignmentId];
            dinara::decompress(compressedAlignment, alignment);

            const array<OrientedReadId, 2> orientedReadIds = {
                OrientedReadId(ad.readIds[0], 0),
                OrientedReadId(ad.readIds[1], ad.isSameStrand ? 0 : 1)
            };

            const ProjectedAlignment projectedAlignment(
                *this,
                orientedReadIds,
                alignment,
                ProjectedAlignment::Method::QuickRawSparse,
                alignOptions.overlapDpMatchScore,
                alignOptions.overlapDpMismatchScore,
                alignOptions.overlapDpGapOpen1,
                alignOptions.overlapDpGapExtend1
            );

            const ReadId r0 = ad.readIds[0];
            const ReadId r1 = ad.readIds[1];
            const uint32_t len1 = uint32_t(reads->getRead(r1).baseCount);
            const bool r1Rev = (orientedReadIds[1].getStrand() == 1);

            for (const auto& m : projectedAlignment.sparseMismatches) {
                const uint32_t pos0 = m.position0;
                const uint32_t pos1Oriented = m.position1;
                if (r1Rev) {
                    DINARA_ASSERT(pos1Oriented < len1);
                }
                const uint32_t pos1 = r1Rev ? ((len1 - 1U) - pos1Oriented) : pos1Oriented;

                const uint64_t n0 = getOrAddNode(r0, pos0);
                const uint64_t n1 = getOrAddNode(r1, pos1);
                dsets.unite(n0, n1);
                totalMismatchEdges++;
            }
        }
    }

    materializeMismatchClustersFromComponents(dsets, result, getReads());

    cout << timestamp << "Global mismatch-site clustering: nodes=" << result.nodes.size()
         << " edges=" << totalMismatchEdges
         << " clusters=" << result.clusterRepresentatives.size();
    if (skippedMismatchedCounts) {
        cout << " skippedInconsistent=" << skippedMismatchedCounts;
    }
    cout << endl;

    return result;
}



// Restrict mismatch-site clustering to the read-graph component reachable from seedReadId.
Assembler::GlobalMismatchSiteClusters Assembler::clusterMismatchingPositionsIntoGlobalHetSitesReachableFromRead(
    ReadId seedReadId,
    const AlignOptions& alignOptions,
    uint64_t threadCount,
    uint64_t maxReadsToProcess,
    uint64_t maxAlignmentsToProcess,
    bool includeDeletedAlignments,
    bool readGraphOnly
) const
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    if (!alignmentTable.isOpen()) {
        throw runtime_error("clusterMismatchingPositionsIntoGlobalHetSitesReachableFromRead requires alignmentTable to be open.");
    }
    if (alignedEvidenceStore.index.size() < alignmentData.size()) {
        throw runtime_error("clusterMismatchingPositionsIntoGlobalHetSitesReachableFromRead requires alignedEvidenceStore to be populated.");
    }
    if (seedReadId >= reads->readCount()) {
        throw runtime_error("clusterMismatchingPositionsIntoGlobalHetSitesReachableFromRead: invalid seedReadId.");
    }

    if (threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    threadCount = std::max<uint64_t>(1, threadCount);

    // For reachable clustering, initial footprint should still be sizeable to avoid
    // growth churn in highly connected components.
    const uint64_t nodeReserveHint =
        (maxAlignmentsToProcess ? (maxAlignmentsToProcess * 4ULL + 1024ULL) : 65536ULL);

    vector<uint8_t> readVisited(reads->readCount(), 0);
    vector<ReadId> queue;
    queue.reserve(1024);
    size_t queueHead = 0;

    readVisited[seedReadId] = 1;
    queue.push_back(seedReadId);
    uint64_t readsProcessed = 0;

    vector<uint8_t> alignmentProcessed(alignmentData.size(), 0);
    uint64_t alignmentsProcessed = 0;
    bool hitReadLimit = false;
    bool hitAlignmentLimit = false;

    vector<uint64_t> selectedAlignmentIds;
    selectedAlignmentIds.reserve(16384);

    while (queueHead < queue.size()) {
        const ReadId currentReadId = queue[queueHead++];
        readsProcessed++;
        if (maxReadsToProcess && readsProcessed > maxReadsToProcess) {
            hitReadLimit = true;
            break;
        }

        const OrientedReadId orientedReadId(currentReadId, 0);
        const span<const uint32_t> section = alignmentTable[orientedReadId.getValue()];

        for (const uint32_t alignmentId32 : section) {
            const uint64_t alignmentId = alignmentId32;
            if (alignmentId >= alignmentProcessed.size()) {
                continue;
            }
            if (alignmentProcessed[alignmentId]) {
                continue;
            }
            alignmentProcessed[alignmentId] = 1;
            alignmentsProcessed++;
            if (maxAlignmentsToProcess && alignmentsProcessed > maxAlignmentsToProcess) {
                hitAlignmentLimit = true;
                break;
            }

            const AlignmentData& ad = alignmentData[alignmentId];
            if (!includeDeletedAlignments && ad.isDeleted()) {
                continue;
            }
            if (readGraphOnly && !ad.info.isInReadGraph) {
                continue;
            }
            if (ad.info.errorRate > float(alignOptions.maxErrorRate)) {
                continue;
            }
            if (ad.info.mismatchCount == 0) {
                continue;
            }
            const uint32_t evidenceId = uint32_t(ad.info.alignmentId);
            if (evidenceId >= alignedEvidenceStore.index.size()) {
                continue;
            }

            // Expand BFS frontier by reads participating in this overlap.
            for (int i = 0; i < 2; i++) {
                const ReadId r = ad.readIds[i];
                if (r >= readVisited.size()) {
                    continue;
                }
                if (!readVisited[r]) {
                    readVisited[r] = 1;
                    queue.push_back(r);
                }
            }
            selectedAlignmentIds.push_back(alignmentId);
        }

        if (hitAlignmentLimit) {
            break;
        }
    }

    DynamicDisjointSets dsets;
    dsets.reserve(nodeReserveHint);
    FlatNodeIdTable nodeIdByKey;
    nodeIdByKey.reserve(nodeReserveHint);

    GlobalMismatchSiteClusters result;
    result.nodes.reserve(size_t(nodeReserveHint));

    uint64_t totalMismatchEdges = 0;
    uint64_t skippedMismatchedCounts = 0;

    auto getOrAddNode = [&](ReadId readId, uint32_t pos) -> uint64_t {
        const uint64_t key = packReadPosKey(readId, pos);
        return nodeIdByKey.getOrCreate(key, [&]() {
            const uint64_t nodeId = dsets.add();
            DINARA_ASSERT(nodeId == result.nodes.size());
            result.nodes.push_back({readId, pos});
            return nodeId;
        });
    };

    auto processAlignmentIntoDsu = [&](uint64_t alignmentId, vector<uint32_t>& reverseScratch) {
        const AlignmentData& ad = alignmentData[alignmentId];
        const uint32_t evidenceId = uint32_t(ad.info.alignmentId);
        const span<const SnpEvidence> qTokens = alignedEvidenceStore.getSnps1(evidenceId);
        const span<const SnpEvidence> tTokens = alignedEvidenceStore.getSnps0(evidenceId);
        const ReadId r0 = ad.readIds[0];
        const ReadId r1 = ad.readIds[1];
        uint64_t mismatchPairs = 0;
        const bool ok = forEachMismatchPair(
            qTokens,
            tTokens,
            !ad.isSameStrand,
            ad.info.mismatchCount,
            reverseScratch,
            [&](uint32_t pos0, uint32_t pos1) {
                const uint64_t n0 = getOrAddNode(r0, pos0);
                const uint64_t n1 = getOrAddNode(r1, pos1);
                dsets.unite(n0, n1);
                mismatchPairs++;
            });
        if (!ok) {
            skippedMismatchedCounts++;
            return;
        }
        totalMismatchEdges += mismatchPairs;
    };

    // BFS discovery is naturally serial; mismatch decoding/pair extraction is parallelizable.
    // We parallelize extraction into packed-key edges and keep DSU/node assignment serial.
    if (threadCount <= 1 || selectedAlignmentIds.size() < 1024) {
        vector<uint32_t> reverseScratch;
        reverseScratch.reserve(1024);
        for (const uint64_t alignmentId : selectedAlignmentIds) {
            processAlignmentIntoDsu(alignmentId, reverseScratch);
        }
    } else {
        const uint64_t workerCount = std::min<uint64_t>(threadCount, selectedAlignmentIds.size());
        vector<vector<pair<uint64_t, uint64_t>>> edgesByWorker(workerCount);
        vector<uint64_t> skippedByWorker(workerCount, 0);
        vector<thread> workers;
        workers.reserve(size_t(workerCount));
        const uint64_t chunkSize = (selectedAlignmentIds.size() + workerCount - 1) / workerCount;

        for (uint64_t workerId = 0; workerId < workerCount; workerId++) {
            workers.emplace_back([&, workerId]() {
                const uint64_t begin = workerId * chunkSize;
                const uint64_t end = std::min<uint64_t>(selectedAlignmentIds.size(), begin + chunkSize);
                auto& localEdges = edgesByWorker[workerId];
                vector<uint32_t> reverseScratch;
                reverseScratch.reserve(1024);

                for (uint64_t i = begin; i < end; i++) {
                    const uint64_t alignmentId = selectedAlignmentIds[i];
                    const AlignmentData& ad = alignmentData[alignmentId];
                    const uint32_t evidenceId = uint32_t(ad.info.alignmentId);
                    const span<const SnpEvidence> qTokens = alignedEvidenceStore.getSnps1(evidenceId);
                    const span<const SnpEvidence> tTokens = alignedEvidenceStore.getSnps0(evidenceId);
                    const ReadId r0 = ad.readIds[0];
                    const ReadId r1 = ad.readIds[1];
                    const bool ok = forEachMismatchPair(
                        qTokens,
                        tTokens,
                        !ad.isSameStrand,
                        ad.info.mismatchCount,
                        reverseScratch,
                        [&](uint32_t pos0, uint32_t pos1) {
                            localEdges.push_back({
                                packReadPosKey(r0, pos0),
                                packReadPosKey(r1, pos1)});
                        });
                    if (!ok) {
                        skippedByWorker[workerId]++;
                    }
                }
            });
        }
        for (auto& worker : workers) {
            worker.join();
        }

        for (uint64_t workerId = 0; workerId < workerCount; workerId++) {
            skippedMismatchedCounts += skippedByWorker[workerId];
            const auto& localEdges = edgesByWorker[workerId];
            totalMismatchEdges += localEdges.size();
            for (const auto& edge : localEdges) {
                const uint64_t key0 = edge.first;
                const uint64_t key1 = edge.second;
                const uint64_t n0 = getOrAddNode(ReadId(key0 >> 32), uint32_t(key0));
                const uint64_t n1 = getOrAddNode(ReadId(key1 >> 32), uint32_t(key1));
                dsets.unite(n0, n1);
            }
        }
    }

    if (dsets.size() == 0) {
        result.clusterMemberOffsets = {0};
        cout << timestamp << "Global mismatch-site clustering (reachable): nodes=0 edges=0 clusters=0" << endl;
        return result;
    }

    materializeMismatchClustersFromComponents(dsets, result, getReads());

    cout << timestamp << "Global mismatch-site clustering (reachable): seed=" << seedReadId
         << " readsProcessed=" << readsProcessed
         << " alignmentsProcessed=" << alignmentsProcessed
         << " nodes=" << result.nodes.size()
         << " edges=" << totalMismatchEdges
         << " clusters=" << result.clusterRepresentatives.size();
    if (skippedMismatchedCounts) {
        cout << " skippedInconsistent=" << skippedMismatchedCounts;
    }
    if (hitReadLimit) {
        cout << " hitReadLimit=" << maxReadsToProcess;
    }
    if (hitAlignmentLimit) {
        cout << " hitAlignmentLimit=" << maxAlignmentsToProcess;
    }
    cout << endl;

    return result;
}



// Traverse read-graph overlaps from a seed (read,position), map through sparse indels,
// and aggregate transitive site members plus allele counts.
Assembler::TransitiveHetSiteCoverage Assembler::gatherTransitiveHetSiteCoverage(
    ReadId seedReadId,
    uint32_t seedPosition,
    const AlignOptions& alignOptions,
    uint64_t maxNodesToVisit,
    uint64_t maxAlignmentsToScan,
    bool includeDeletedAlignments
) const
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    if (alignedEvidenceStore.index.size() < alignmentData.size()) {
        throw runtime_error("gatherTransitiveHetSiteCoverage requires alignedEvidenceStore to be populated.");
    }

    TransitiveHetSiteCoverage out;
    if (seedReadId >= reads->readCount()) {
        return out;
    }
    const uint64_t seedLen = reads->getRead(seedReadId).baseCount;
    if (uint64_t(seedPosition) >= seedLen) {
        return out;
    }

    unordered_map<uint32_t, uint32_t> positionByRead;
    positionByRead.reserve(1024);
    vector< pair<ReadId, uint32_t> > queue;
    queue.reserve(1024);
    size_t queueHead = 0;

    positionByRead.insert({seedReadId, seedPosition});
    queue.push_back({seedReadId, seedPosition});

    ThreadLocalAlignmentCache& cache = getThreadLocalAlignmentCache(*this, alignmentData.size());
    const uint32_t scanEpoch = beginAlignmentScanEpoch(cache);
    const Reads& readsRef = getReads();
    const uint32_t kHalf = uint32_t(assemblerInfo->k / 2);

    while (queueHead < queue.size()) {
        if (maxNodesToVisit && uint64_t(queue.size()) >= maxNodesToVisit) {
            out.hitNodeLimit = true;
            break;
        }

        const auto [currentReadId, currentPos] = queue[queueHead++];
        const OrientedReadId orientedReadId(currentReadId, 0);
        const span<const uint32_t> section = alignmentTable[orientedReadId.getValue()];

        for (const uint32_t alignmentId32 : section) {
            if (maxAlignmentsToScan && out.alignmentsScanned >= maxAlignmentsToScan) {
                out.hitAlignmentLimit = true;
                break;
            }

            const uint64_t alignmentId = alignmentId32;
            if (alignmentId >= alignmentData.size() || alignmentId >= cache.alignmentSeenEpoch.size()) {
                continue;
            }
            if (cache.alignmentSeenEpoch[alignmentId] == scanEpoch) {
                continue;
            }
            cache.alignmentSeenEpoch[alignmentId] = scanEpoch;
            out.alignmentsScanned++;

            const AlignmentData& ad = alignmentData[alignmentId];
            if (!includeDeletedAlignments && ad.isDeleted()) {
                continue;
            }
            if (ad.info.errorRate > float(alignOptions.maxErrorRate)) {
                continue;
            }

            const uint32_t evidenceId = uint32_t(ad.info.alignmentId);
            if (evidenceId >= alignedEvidenceStore.index.size()) {
                continue;
            }

            const AlignmentAnchors& anchors =
                getCachedAlignmentAnchors(cache, *this, alignmentId, ad, kHalf);
            if (!anchors.ok) {
                continue;
            }

            const int dir = ad.isSameStrand ? 1 : -1;

            ReadId otherReadId = invalidReadId;
            MapResult mapped;
            if (currentReadId == ad.readIds[0]) {
                // Map read0 -> read1 (use indels stored on read0 coordinate system).
                if (currentPos < ad.qs || currentPos >= ad.qe) {
                    continue;
                }
                otherReadId = ad.readIds[1];
                mapped = mapPositionUsingSparseIndels(
                    currentPos,
                    anchors.first0,
                    anchors.first1,
                    dir,
                    alignedEvidenceStore.getIndels1(evidenceId));
            } else if (currentReadId == ad.readIds[1]) {
                // Map read1 -> read0 (use indels stored on read1 coordinate system).
                if (currentPos < ad.ts || currentPos >= ad.te) {
                    continue;
                }
                otherReadId = ad.readIds[0];
                const uint32_t anchorFrom = (anchors.strand1 == 0) ? anchors.first1 : anchors.last1;
                const uint32_t anchorTo = (anchors.strand1 == 0) ? anchors.first0 : anchors.last0;
                mapped = mapPositionUsingSparseIndels(
                    currentPos,
                    anchorFrom,
                    anchorTo,
                    dir,
                    alignedEvidenceStore.getIndels0(evidenceId));
            } else {
                continue;
            }

            if (mapped.hole) {
                out.mappingHoles++;
                continue;
            }
            if (!mapped.ok) {
                continue;
            }

            const uint64_t otherLen = reads->getRead(otherReadId).baseCount;
            if (uint64_t(mapped.mappedPos) >= otherLen) {
                continue;
            }

            const auto it = positionByRead.find(otherReadId);
            if (it == positionByRead.end()) {
                if (maxNodesToVisit && uint64_t(queue.size()) + 1 >= maxNodesToVisit) {
                    out.hitNodeLimit = true;
                    continue;
                }
                positionByRead.insert({otherReadId, mapped.mappedPos});
                queue.push_back({otherReadId, mapped.mappedPos});
            } else {
                if (it->second != mapped.mappedPos) {
                    out.mappingConflicts++;
                }
            }
        }

        if (out.hitAlignmentLimit) {
            break;
        }
    }

    out.members = queue;
    for (const auto& node : out.members) {
        const ReadId readId = node.first;
        const uint32_t pos = node.second;
        const uint8_t base = readsRef.getOrientedReadBase(OrientedReadId(readId, 0), pos).value;
        if (base < 4) {
            out.alleleCounts[base]++;
        }
    }

    return out;
}


// Debug verifier for one global site: remap expected members independently and report
// agreement, holes, and mapping failures.
Assembler::GlobalHetSitePositionVerificationStats Assembler::debugVerifyGlobalHetSitePositionsUsingReadGraph(
    const GlobalMismatchSiteClusters& clusters,
    const GlobalHetSiteAlleleMembers& members,
    uint32_t siteId,
    const AlignOptions& alignOptions,
    uint64_t maxNodesToVisit,
    uint64_t maxAlignmentsToScan,
    bool includeDeletedAlignments
) const
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    if (!alignmentTable.isOpen()) {
        throw runtime_error("debugVerifyGlobalHetSitePositionsUsingReadGraph requires alignmentTable to be open.");
    }
    if (alignedEvidenceStore.index.size() < alignmentData.size()) {
        throw runtime_error("debugVerifyGlobalHetSitePositionsUsingReadGraph requires alignedEvidenceStore to be populated.");
    }

    GlobalHetSitePositionVerificationStats out;
    out.siteId = siteId;
    if (siteId >= members.offsets.size()) {
        return out;
    }
    if (siteId + 1 >= clusters.clusterMemberOffsets.size()) {
        return out;
    }

    // Expected per-read positions for this siteId (from the propagated member table).
    unordered_map<ReadId, uint32_t> expected;
    {
        const auto& off = members.offsets[siteId];
        const uint64_t total = off[4] - off[0];
        expected.reserve(size_t(total * 10 / 7 + 8));
        for (int allele = 0; allele < 4; allele++) {
            for (uint64_t i = off[allele]; i < off[allele + 1]; i++) {
                expected.insert({members.members[i].readId, members.members[i].position});
            }
        }
    }
    out.expectedMembers = expected.size();

    // Multi-source BFS from mismatch members of this site.
    unordered_map<ReadId, uint32_t> reached;
    reached.reserve(expected.size() ? expected.size() : 128);
    deque< pair<ReadId, uint32_t> > queue;
    {
        const uint64_t begin = clusters.clusterMemberOffsets[siteId];
        const uint64_t end = clusters.clusterMemberOffsets[siteId + 1];
        for (uint64_t i = begin; i < end; i++) {
            const auto& node = clusters.nodes[clusters.clusterMembers[i]];
            if (reached.insert({node.first, node.second}).second) {
                queue.push_back({node.first, node.second});
            }
        }
    }

    const uint32_t kHalf = uint32_t(assemblerInfo->k / 2);
    ThreadLocalAlignmentCache& cache = getThreadLocalAlignmentCache(*this, alignmentData.size());

    while (!queue.empty()) {
        if (maxNodesToVisit && reached.size() >= maxNodesToVisit) {
            out.hitNodeLimit = true;
            break;
        }

        const auto [readId, pos] = queue.front();
        queue.pop_front();

        const OrientedReadId orientedReadId(readId, 0);
        const span<const uint32_t> section = alignmentTable[orientedReadId.getValue()];

        for (const uint32_t alignmentId32 : section) {
            if (maxAlignmentsToScan && out.checkedMappings >= maxAlignmentsToScan) {
                out.hitAlignmentLimit = true;
                break;
            }

            const uint64_t alignmentId = alignmentId32;
            if (alignmentId >= alignmentData.size()) {
                continue;
            }
            const AlignmentData& ad = alignmentData[alignmentId];
            if (!includeDeletedAlignments && ad.isDeleted()) {
                continue;
            }
            if (!ad.info.isInReadGraph) {
                continue;
            }
            if (ad.info.errorRate > float(alignOptions.maxErrorRate)) {
                continue;
            }

            const uint32_t evidenceId = uint32_t(ad.info.alignmentId);
            if (evidenceId >= alignedEvidenceStore.index.size()) {
                continue;
            }

            const AlignmentAnchors& a =
                getCachedAlignmentAnchors(cache, *this, alignmentId, ad, kHalf);
            if (!a.ok) {
                continue;
            }

            const int dir = ad.isSameStrand ? 1 : -1;

            ReadId otherReadId = invalidReadId;
            MapResult mapped;
            if (readId == ad.readIds[0]) {
                const uint32_t fromBegin = std::min(a.first0, a.last0);
                const uint32_t fromEnd = std::max(a.first0, a.last0);
                if (pos < fromBegin || pos >= fromEnd) {
                    continue;
                }
                otherReadId = ad.readIds[1];
                mapped = mapPositionUsingSparseIndels(
                    pos,
                    a.first0,
                    a.first1,
                    dir,
                    alignedEvidenceStore.getIndels1(evidenceId));
            } else if (readId == ad.readIds[1]) {
                const uint32_t fromBegin = std::min(a.first1, a.last1);
                const uint32_t fromEnd = std::max(a.first1, a.last1);
                if (pos < fromBegin || pos >= fromEnd) {
                    continue;
                }
                otherReadId = ad.readIds[0];
                const uint32_t anchorFrom = (a.strand1 == 0) ? a.first1 : a.last1;
                const uint32_t anchorTo = (a.strand1 == 0) ? a.first0 : a.last0;
                mapped = mapPositionUsingSparseIndels(
                    pos,
                    anchorFrom,
                    anchorTo,
                    dir,
                    alignedEvidenceStore.getIndels0(evidenceId));
            } else {
                continue;
            }

            out.checkedMappings++;

            if (mapped.hole) {
                out.mappingHoles++;
                continue;
            }
            if (!mapped.ok) {
                out.mappingFailures++;
                continue;
            }

            const uint64_t otherLen = reads->getRead(otherReadId).baseCount;
            if (uint64_t(mapped.mappedPos) >= otherLen) {
                continue;
            }

            if (const auto it = expected.find(otherReadId); it != expected.end()) {
                if (it->second != mapped.mappedPos) {
                    out.mismatchedPositions++;
                }
            }

            if (reached.insert({otherReadId, mapped.mappedPos}).second) {
                queue.push_back({otherReadId, mapped.mappedPos});
            }
        }

        if (out.hitAlignmentLimit) {
            break;
        }
    }

    // Count how many expected members were reached.
    for (const auto& p : expected) {
        if (reached.find(p.first) != reached.end()) {
            out.reachedMembers++;
        }
    }

    return out;
}



vector<int8_t> Assembler::computeReadGraphStrandsFromSeed(
    ReadId seedReadId,
    uint64_t& conflicts,
    bool includeDeletedAlignments
) const
{
    reads->checkReadsAreOpen();
    checkAlignmentDataAreOpen();
    if (!alignmentTable.isOpen()) {
        throw runtime_error("computeReadGraphStrandsFromSeed requires alignmentTable to be open.");
    }

    conflicts = 0;
    vector<int8_t> strandByRead(reads->readCount(), int8_t(-1));
    if (seedReadId >= reads->readCount()) {
        return strandByRead;
    }

    deque<ReadId> q;
    strandByRead[uint64_t(seedReadId)] = 0;
    q.push_back(seedReadId);

    while (!q.empty()) {
        const ReadId r = q.front();
        q.pop_front();

        const OrientedReadId orientedReadId(r, 0);
        const span<const uint32_t> section = alignmentTable[orientedReadId.getValue()];
        for (const uint32_t alignmentId32 : section) {
            const uint64_t alignmentId = alignmentId32;
            if (alignmentId >= alignmentData.size()) {
                continue;
            }
            const AlignmentData& ad = alignmentData[alignmentId];
            if (!includeDeletedAlignments && ad.isDeleted()) {
                continue;
            }
            if (!ad.info.isInReadGraph) {
                continue;
            }

            ReadId other = invalidReadId;
            if (ad.readIds[0] == r) {
                other = ad.readIds[1];
            } else if (ad.readIds[1] == r) {
                other = ad.readIds[0];
            } else {
                continue;
            }

            const int8_t delta = ad.isSameStrand ? int8_t(0) : int8_t(1);
            const int8_t implied = int8_t(strandByRead[uint64_t(r)] ^ delta);
            int8_t& dst = strandByRead[uint64_t(other)];
            if (dst == -1) {
                dst = implied;
                q.push_back(other);
            } else if (dst != implied) {
                conflicts++;
            }
        }
    }

    return strandByRead;
}



// Re-orient propagated members into a focal strand frame.
Assembler::GlobalHetSiteOrientedAlleleMembers Assembler::orientGlobalHetSiteAlleleMembers(
    const GlobalHetSiteAlleleMembers& forwardMembers,
    const vector<int8_t>& strandByRead
) const
{
    reads->checkReadsAreOpen();

    GlobalHetSiteOrientedAlleleMembers out;
    out.offsets = forwardMembers.offsets;
    out.members.clear();

    if (forwardMembers.offsets.empty()) {
        return out;
    }

    static constexpr uint8_t complementBase[4] = {3, 2, 1, 0};
    const uint32_t siteCount = uint32_t(forwardMembers.offsets.size());
    const uint64_t readCount = reads->readCount();
    vector<uint32_t> readLengthByRead(readCount, 0);
    for (uint64_t r = 0; r < readCount; r++) {
        readLengthByRead[r] = uint32_t(getReads().getRead(ReadId(r)).baseCount);
    }

    // First pass: count members per site per oriented allele.
    vector< array<uint32_t, 4> > counts(siteCount, array<uint32_t, 4>{0, 0, 0, 0});
    for (uint32_t siteId = 0; siteId < siteCount; siteId++) {
        const auto& off = forwardMembers.offsets[siteId];
        for (int allele = 0; allele < 4; allele++) {
            for (uint64_t i = off[allele]; i < off[allele + 1]; i++) {
                const auto& m = forwardMembers.members[i];
                const ReadId rid = m.readId;
                const uint32_t posFwd = m.position;
                if (uint64_t(rid) >= readLengthByRead.size()) {
                    continue;
                }
                const uint32_t len = readLengthByRead[uint64_t(rid)];
                if (len == 0 || posFwd >= len) {
                    continue;
                }
                const int8_t s = (uint64_t(rid) < strandByRead.size()) ? strandByRead[uint64_t(rid)] : int8_t(-1);
                const uint8_t b = (s == 1) ? complementBase[uint8_t(allele)] : uint8_t(allele);
                counts[siteId][b]++;
            }
        }
    }

    // Prefix sums to offsets.
    out.offsets.assign(siteCount, array<uint64_t, 5>{0, 0, 0, 0, 0});
    uint64_t total = 0;
    for (uint32_t siteId = 0; siteId < siteCount; siteId++) {
        out.offsets[siteId][0] = total;
        total += counts[siteId][0];
        out.offsets[siteId][1] = total;
        total += counts[siteId][1];
        out.offsets[siteId][2] = total;
        total += counts[siteId][2];
        out.offsets[siteId][3] = total;
        total += counts[siteId][3];
        out.offsets[siteId][4] = total;
    }

    out.members.assign(total, GlobalHetSiteOrientedAlleleMembers::Member{});
    vector< array<uint32_t, 4> > writeIndex(siteCount, array<uint32_t, 4>{0, 0, 0, 0});

    // Second pass: write members into oriented-allele buckets.
    for (uint32_t siteId = 0; siteId < siteCount; siteId++) {
        const auto& off = forwardMembers.offsets[siteId];
        for (int allele = 0; allele < 4; allele++) {
            for (uint64_t i = off[allele]; i < off[allele + 1]; i++) {
                const auto& m = forwardMembers.members[i];
                const ReadId rid = m.readId;
                const uint32_t posFwd = m.position;
                if (uint64_t(rid) >= readLengthByRead.size()) {
                    continue;
                }
                const uint32_t len = readLengthByRead[uint64_t(rid)];
                const int8_t s = (uint64_t(rid) < strandByRead.size()) ? strandByRead[uint64_t(rid)] : int8_t(-1);
                const Strand strand = (s == 1) ? 1 : 0;
                const uint8_t b = (strand == 1) ? complementBase[uint8_t(allele)] : uint8_t(allele);

                if (len == 0 || posFwd >= len) {
                    continue;
                }
                const uint32_t posOriented = (strand == 1) ? ((len - 1U) - posFwd) : posFwd;

                const uint64_t begin = out.offsets[siteId][b];
                const uint32_t w = writeIndex[siteId][b]++;
                const uint64_t idx = begin + w;
                if (idx >= out.members.size()) {
                    continue;
                }
                out.members[idx] = GlobalHetSiteOrientedAlleleMembers::Member{
                    OrientedReadId(rid, strand),
                    posOriented
                };
            }
        }
    }

    return out;
}



// Recompute cluster allele counts using transitive read-graph coverage.
void Assembler::updateGlobalMismatchSiteClusterAlleleCountsWithTransitiveCoverage(
    GlobalMismatchSiteClusters& clusters,
    const AlignOptions& alignOptions,
    uint64_t maxNodesToVisitPerSite,
    uint64_t maxAlignmentsToScanPerSite,
    ReadId restrictToClustersInvolvingRead,
    uint64_t maxClustersToUpdate,
    bool includeDeletedAlignments
) const
{
    if (clusters.alleleCountsAreTransitive.size() != clusters.alleleCounts.size()) {
        clusters.alleleCountsAreTransitive.assign(clusters.alleleCounts.size(), 0);
    }
    if (clusters.transitiveSiteMemberCounts.size() != clusters.alleleCounts.size()) {
        clusters.transitiveSiteMemberCounts.assign(clusters.alleleCounts.size(), 0);
    }

    uint64_t updated = 0;
    for (size_t clusterId = 0;
        clusterId + 1 < clusters.clusterMemberOffsets.size();
        clusterId++) {

        if (maxClustersToUpdate && updated >= maxClustersToUpdate) {
            break;
        }
        if (clusterId < clusters.alleleCountsAreTransitive.size() &&
            clusters.alleleCountsAreTransitive[clusterId]) {
            continue;
        }

        const uint64_t begin = clusters.clusterMemberOffsets[clusterId];
        const uint64_t end = clusters.clusterMemberOffsets[clusterId + 1];
        if (begin >= end) {
            continue;
        }

        ReadId seedReadId = clusters.nodes[clusters.clusterMembers[begin]].first;
        uint32_t seedPos = clusters.nodes[clusters.clusterMembers[begin]].second;

        if (restrictToClustersInvolvingRead != invalidReadId) {
            bool found = false;
            for (uint64_t i = begin; i < end; i++) {
                const auto& node = clusters.nodes[clusters.clusterMembers[i]];
                if (node.first == restrictToClustersInvolvingRead) {
                    seedReadId = node.first;
                    seedPos = node.second;
                    found = true;
                    break;
                }
            }
            if (!found) {
                continue;
            }
        } else {
            // Choose a seed node that minimizes overlap scanning work (heuristic):
            // pick the member whose read has the smallest alignmentTable[read-strand0] degree.
            if (alignmentTable.isOpen()) {
                size_t bestDegree = std::numeric_limits<size_t>::max();
                for (uint64_t i = begin; i < end; i++) {
                    const auto& node = clusters.nodes[clusters.clusterMembers[i]];
                    const OrientedReadId r(node.first, 0);
                    const size_t degree = alignmentTable[r.getValue()].size();
                    if (degree < bestDegree) {
                        bestDegree = degree;
                        seedReadId = node.first;
                        seedPos = node.second;
                    }
                }
            }
        }

        const auto coverage = gatherTransitiveHetSiteCoverage(
            seedReadId,
            seedPos,
            alignOptions,
            maxNodesToVisitPerSite,
            maxAlignmentsToScanPerSite,
            includeDeletedAlignments
        );
        clusters.alleleCounts[clusterId] = coverage.alleleCounts;
        clusters.alleleCountsAreTransitive[clusterId] = 1;
        clusters.transitiveSiteMemberCounts[clusterId] = uint32_t(coverage.members.size());
        updated++;
    }
}



// Propagate global-site memberships over read-graph overlaps.
// Per (site,read) policy:
// - keep up to two competing mapped positions with support counts;
// - drop if a hole is observed;
// - drop if competing tracked positions imply different bases.
Assembler::GlobalHetSiteAlleleMembers Assembler::computeGlobalHetSiteAlleleMembersUsingReadGraph(
    const GlobalMismatchSiteClusters& clusters,
    const AlignOptions& alignOptions,
    uint64_t maxPendingTasks,
    bool includeDeletedAlignments,
    ReadId seedReadId
) const
{
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    if (!alignmentTable.isOpen()) {
        throw runtime_error("computeGlobalHetSiteAlleleMembersUsingReadGraph requires alignmentTable to be open.");
    }
    if (alignedEvidenceStore.index.size() < alignmentData.size()) {
        throw runtime_error("computeGlobalHetSiteAlleleMembersUsingReadGraph requires alignedEvidenceStore to be populated.");
    }

    GlobalHetSiteAlleleMembers out;
    const uint32_t siteCount =
        clusters.clusterMemberOffsets.empty() ?
        0U :
        uint32_t(clusters.clusterMemberOffsets.size() - 1);
    if (siteCount == 0) {
        out.offsets.clear();
        out.members.clear();
        return out;
    }

    const uint32_t kHalf = uint32_t(assemblerInfo->k / 2);
    ThreadLocalAlignmentCache& cache = getThreadLocalAlignmentCache(*this, alignmentData.size());

    struct SitePos {
        uint32_t position = 0;
        uint32_t siteId = 0;
    };

    const uint64_t readCount = reads->readCount();
    vector<uint32_t> readLengthByRead(readCount, 0);
    for (uint64_t r = 0; r < readCount; r++) {
        readLengthByRead[r] = uint32_t(reads->getRead(ReadId(r)).baseCount);
    }
    vector< vector<SitePos> > pendingByRead(readCount);
    vector<uint8_t> pendingIsSorted(readCount, 1);
    vector<uint8_t> readInQueue(readCount, 0);
    deque<ReadId> readQueue;
    uint64_t enqueuedTaskCount = 0;
    // Strict no-gap policy only needs hole presence (not counts).
    unordered_set<uint64_t> holeKeys;
    holeKeys.reserve(uint64_t(clusters.clusterMembers.size() / 4 + 1024));
    const Reads& readsRef = getReads();

    FlatAssignmentTable table;
    table.reserve(uint64_t(clusters.clusterMembers.size()) * 4 + 1024); // heuristic

    // Phase 1 helper: register a discovered (site,read)->position assignment.
    auto enqueueAssignment = [&](ReadId readId, uint32_t siteId, uint32_t position) {
        if (maxPendingTasks && enqueuedTaskCount >= maxPendingTasks) {
            return;
        }
        const uint64_t key = packSiteReadKey(siteId, readId);
        const auto observed = table.observe(key, position);
        if (observed.conflict) {
            out.mappingConflicts++;
        }
        // Propagate only when this (site,read) key is first discovered.
        // Alternate candidate positions still contribute support/conflict bookkeeping,
        // but do not fan out additional propagation branches.
        if (!observed.insertedKey) {
            return;
        }
        out.propagatedAssignments++;

        vector<SitePos>& pending = pendingByRead[uint64_t(readId)];
        if (pendingIsSorted[uint64_t(readId)] && !pending.empty()) {
            const SitePos& last = pending.back();
            if (last.position > position ||
                (last.position == position && last.siteId > siteId)) {
                pendingIsSorted[uint64_t(readId)] = 0;
            }
        }
        pending.push_back(SitePos{position, siteId});
        enqueuedTaskCount++;

        if (!readInQueue[uint64_t(readId)]) {
            readInQueue[uint64_t(readId)] = 1;
            readQueue.push_back(readId);
        }
    };

    // Phase 1: seed propagation from mismatch cluster members (or focal seed only).
    auto seedAssignments = [&]() {
        for (uint32_t siteId = 0; siteId < siteCount; siteId++) {
            const uint64_t begin = clusters.clusterMemberOffsets[siteId];
            const uint64_t end = clusters.clusterMemberOffsets[siteId + 1];
            if (seedReadId != invalidReadId) {
                uint32_t bestPos = std::numeric_limits<uint32_t>::max();
                for (uint64_t i = begin; i < end; i++) {
                    const auto& node = clusters.nodes[clusters.clusterMembers[i]];
                    if (node.first == seedReadId && node.second < bestPos) {
                        bestPos = node.second;
                    }
                }
                if (bestPos != std::numeric_limits<uint32_t>::max()) {
                    enqueueAssignment(seedReadId, siteId, bestPos);
                }
            } else {
                for (uint64_t i = begin; i < end; i++) {
                    const auto& node = clusters.nodes[clusters.clusterMembers[i]];
                    enqueueAssignment(node.first, siteId, node.second);
                }
            }
            if (maxPendingTasks && enqueuedTaskCount >= maxPendingTasks) {
                break;
            }
        }
    };

    // Phase 2: multi-source propagation on read-graph overlaps.
    auto propagateAssignmentsOnReadGraph = [&]() {
        while (!readQueue.empty()) {
            const ReadId readId = readQueue.front();
            readQueue.pop_front();
            readInQueue[uint64_t(readId)] = 0;

            vector<SitePos>& pending = pendingByRead[uint64_t(readId)];
            if (pending.empty()) {
                continue;
            }
            if (!pendingIsSorted[uint64_t(readId)]) {
                sort(pending.begin(), pending.end(),
                    [](const SitePos& a, const SitePos& b) {
                        if (a.position != b.position) {
                            return a.position < b.position;
                        }
                        return a.siteId < b.siteId;
                    });
                pendingIsSorted[uint64_t(readId)] = 1;
            }

            const OrientedReadId orientedReadId(readId, 0);
            const span<const uint32_t> section = alignmentTable[orientedReadId.getValue()];

            for (const uint32_t alignmentId32 : section) {
                const uint64_t alignmentId = alignmentId32;
                if (alignmentId >= alignmentData.size()) {
                    continue;
                }
                const AlignmentData& ad = alignmentData[alignmentId];
                if (!ad.info.isInReadGraph) {
                    continue;
                }
                if (!includeDeletedAlignments && ad.isDeleted()) {
                    continue;
                }
                if (ad.info.errorRate > float(alignOptions.maxErrorRate)) {
                    continue;
                }

                const uint32_t evidenceId = uint32_t(ad.info.alignmentId);
                if (evidenceId >= alignedEvidenceStore.index.size()) {
                    continue;
                }

                const AlignmentAnchors& a =
                    getCachedAlignmentAnchors(cache, *this, alignmentId, ad, kHalf);
                if (!a.ok) {
                    continue;
                }

                ReadId otherReadId = invalidReadId;
                uint32_t intervalBegin = 0;
                uint32_t intervalEnd = 0;
                uint32_t anchorFrom = 0;
                uint32_t anchorTo = 0;
                int dir = ad.isSameStrand ? 1 : -1;
                span<const IndelEvidence> indels;

                if (readId == ad.readIds[0]) {
                    otherReadId = ad.readIds[1];
                    intervalBegin = ad.qs;
                    intervalEnd = ad.qe;
                    anchorFrom = a.first0;
                    anchorTo = a.first1;
                    indels = alignedEvidenceStore.getIndels1(evidenceId); // on read0 coords
                } else if (readId == ad.readIds[1]) {
                    otherReadId = ad.readIds[0];
                    intervalBegin = ad.ts;
                    intervalEnd = ad.te;
                    anchorFrom = (a.strand1 == 0) ? a.first1 : a.last1;
                    anchorTo = (a.strand1 == 0) ? a.first0 : a.last0;
                    indels = alignedEvidenceStore.getIndels0(evidenceId); // on read1 coords
                } else {
                    continue;
                }
                if (intervalBegin >= intervalEnd) {
                    continue;
                }

                // Find pending entries in this overlap interval.
                const auto lb = lower_bound(
                    pending.begin(),
                    pending.end(),
                    intervalBegin,
                    [](const SitePos& x, uint32_t v) { return x.position < v; });
                const auto ub = lower_bound(
                    pending.begin(),
                    pending.end(),
                    intervalEnd,
                    [](const SitePos& x, uint32_t v) { return x.position < v; });
                if (lb == ub) {
                    continue;
                }

                // Map selected positions in one pass over sparse indels.
                int64_t netShift = 0;
                size_t indelIndex = 0;
                bool pendingDeletion = false;
                uint64_t pendingDeletionEnd = 0;
                int64_t pendingDeletionShift = 0; // negative len

                const uint64_t otherLen = readLengthByRead[uint64_t(otherReadId)];

                for (auto it = lb; it != ub; ++it) {
                    const uint32_t fromPos = it->position;
                    const uint32_t siteId = it->siteId;

                    if (fromPos < anchorFrom) {
                        continue;
                    }

                    while (pendingDeletion && uint64_t(fromPos) >= pendingDeletionEnd) {
                        netShift += pendingDeletionShift;
                        pendingDeletion = false;
                    }

                    while (indelIndex < indels.size() && indels[indelIndex].pos() <= fromPos) {
                        const IndelEvidence ev = indels[indelIndex++];
                        const uint32_t p = ev.pos();
                        const uint32_t len = ev.len();
                        if (ev.isInsertion()) {
                            netShift += int64_t(len);
                        } else {
                            const uint64_t end = uint64_t(p) + uint64_t(len);
                            if (uint64_t(fromPos) < end) {
                                pendingDeletion = true;
                                pendingDeletionEnd = end;
                                pendingDeletionShift = -int64_t(len);
                                break;
                            } else {
                                netShift -= int64_t(len);
                            }
                        }
                    }

                    if (pendingDeletion && uint64_t(fromPos) < pendingDeletionEnd) {
                        out.mappingHoles++;
                        const uint64_t key = packSiteReadKey(siteId, otherReadId);
                        holeKeys.insert(key);
                        continue;
                    }

                    const int64_t mapped =
                        int64_t(anchorTo) +
                        int64_t(dir) * ((int64_t(fromPos) - int64_t(anchorFrom)) + netShift);
                    if (mapped < 0 || mapped >= int64_t(otherLen)) {
                        continue;
                    }
                    enqueueAssignment(otherReadId, siteId, uint32_t(mapped));
                }

                if (maxPendingTasks && enqueuedTaskCount >= maxPendingTasks) {
                    break;
                }
            }

            pending.clear();
            pendingIsSorted[uint64_t(readId)] = 1;

            if (maxPendingTasks && enqueuedTaskCount >= maxPendingTasks) {
                break;
            }
        }
    };

    // Phase 3: materialize final per-site per-allele members after consistency filtering.
    auto materializeMembers = [&]() {
        out.offsets.assign(siteCount, array<uint64_t, 5>{0, 0, 0, 0, 0});
        vector< array<uint32_t, 4> > counts(siteCount, array<uint32_t, 4>{0, 0, 0, 0});

        struct StagedMember {
            uint32_t sid = 0;
            ReadId rid = invalidReadId;
            uint32_t position = 0;
            uint8_t allele = 0;
        };
        vector<StagedMember> staged;
        staged.reserve(size_t(table.size()));

        table.forEachOccupiedChosen([&](
            uint64_t key,
            uint32_t position,
            uint32_t otherPosition,
            uint32_t otherSupport) {
            const uint32_t sid = uint32_t(key >> 32);
            const ReadId rid = ReadId(uint32_t(key & 0xffffffffU));
            if (sid >= siteCount) {
                return;
            }
            if (uint64_t(rid) >= readLengthByRead.size()) {
                return;
            }
            const uint32_t readLength = readLengthByRead[uint64_t(rid)];
            if (position >= readLength) {
                return;
            }

            // Strict no-gap policy: if this (site,read) was ever observed as a hole
            // during propagation, do not emit a base member for it.
            if (holeKeys.find(key) != holeKeys.end()) {
                out.mappingConflicts++;
                return;
            }

            const uint8_t base = readsRef.getOrientedReadBase(OrientedReadId(rid, 0), position).value;
            if (base >= 4) {
                out.mappingConflicts++;
                return;
            }
            // Require exactly one consistent base across tracked competing positions.
            if (otherSupport > 0 && otherPosition < readLength) {
                const uint8_t otherBase =
                    readsRef.getOrientedReadBase(OrientedReadId(rid, 0), otherPosition).value;
                if (otherBase >= 4 || otherBase != base) {
                    out.mappingConflicts++;
                    return;
                }
            }

            counts[sid][base]++;
            staged.push_back(StagedMember{sid, rid, position, base});
        });

        uint64_t total = 0;
        for (uint32_t sid = 0; sid < siteCount; sid++) {
            out.offsets[sid][0] = total;
            total += counts[sid][0];
            out.offsets[sid][1] = total;
            total += counts[sid][1];
            out.offsets[sid][2] = total;
            total += counts[sid][2];
            out.offsets[sid][3] = total;
            total += counts[sid][3];
            out.offsets[sid][4] = total;
        }

        out.members.assign(total, GlobalHetSiteAlleleMembers::Member{});
        vector< array<uint32_t, 4> > writeIndex(siteCount, array<uint32_t, 4>{0, 0, 0, 0});

        for (const auto& m : staged) {
            const uint64_t begin = out.offsets[m.sid][m.allele];
            const uint32_t w = writeIndex[m.sid][m.allele]++;
            const uint64_t idx = begin + w;
            if (idx >= out.members.size()) {
                continue;
            }
            out.members[idx] = GlobalHetSiteAlleleMembers::Member{m.rid, m.position};
        }
    };

    seedAssignments();
    propagateAssignmentsOnReadGraph();
    materializeMembers();

    return out;
}


Assembler::GlobalHetSiteReadIndex Assembler::buildFilteredGlobalHetSiteReadIndex(
    const GlobalHetSiteAlleleMembers& members,
    uint32_t minAlleleSupport,
    uint32_t minAllelesWithMinSupport
) const
{
    reads->checkReadsAreOpen();

    GlobalHetSiteReadIndex out;
    const uint32_t siteCount = uint32_t(members.offsets.size());
    const uint64_t readCount = reads->readCount();

    out.siteAlleleCounts.assign(siteCount, array<uint32_t, 4>{0, 0, 0, 0});
    out.sitePassesFilter.assign(siteCount, 0);
    out.sitesByRead.assign(readCount, vector<GlobalHetSiteReadIndex::ReadSite>{});
    if (siteCount == 0 || readCount == 0) {
        return out;
    }

    vector<uint32_t> readLengthByRead(readCount, 0);
    for (uint64_t r = 0; r < readCount; r++) {
        readLengthByRead[r] = uint32_t(reads->getRead(ReadId(r)).baseCount);
    }

    for (uint32_t siteId = 0; siteId < siteCount; siteId++) {
        const auto& off = members.offsets[siteId];
        const array<uint32_t, 4> counts{
            uint32_t(off[1] - off[0]),
            uint32_t(off[2] - off[1]),
            uint32_t(off[3] - off[2]),
            uint32_t(off[4] - off[3])
        };
        out.siteAlleleCounts[siteId] = counts;

        uint32_t supportedAlleles = 0;
        for (int allele = 0; allele < 4; allele++) {
            if (counts[allele] >= minAlleleSupport) {
                supportedAlleles++;
            }
        }
        if (supportedAlleles < minAllelesWithMinSupport) {
            continue;
        }

        out.sitePassesFilter[siteId] = 1;
        out.keptSiteCount++;

        for (int allele = 0; allele < 4; allele++) {
            const uint64_t b0 = off[allele];
            const uint64_t b1 = off[allele + 1];
            for (uint64_t i = b0; i < b1; i++) {
                const auto& member = members.members[i];
                const uint64_t rid = uint64_t(member.readId);
                if (rid >= readCount || member.position >= readLengthByRead[rid]) {
                    out.droppedInvalidReadSiteCount++;
                    continue;
                }
                out.sitesByRead[rid].push_back(GlobalHetSiteReadIndex::ReadSite{
                    siteId,
                    member.position,
                    uint8_t(allele)
                });
            }
        }
    }

    for (uint64_t rid = 0; rid < readCount; rid++) {
        auto& readSites = out.sitesByRead[rid];
        if (readSites.empty()) {
            continue;
        }

        sort(readSites.begin(), readSites.end(),
            [](const GlobalHetSiteReadIndex::ReadSite& a, const GlobalHetSiteReadIndex::ReadSite& b) {
                if (a.siteId != b.siteId) {
                    return a.siteId < b.siteId;
                }
                if (a.readPosition != b.readPosition) {
                    return a.readPosition < b.readPosition;
                }
                return a.allele < b.allele;
            });

        size_t write = 0;
        size_t begin = 0;
        while (begin < readSites.size()) {
            size_t end = begin + 1;
            bool consistent = true;
            const uint32_t refPosition = readSites[begin].readPosition;
            const uint8_t refAllele = readSites[begin].allele;
            const uint32_t siteId = readSites[begin].siteId;

            while (end < readSites.size() && readSites[end].siteId == siteId) {
                if (readSites[end].readPosition != refPosition ||
                    readSites[end].allele != refAllele) {
                    consistent = false;
                }
                end++;
            }

            if (consistent) {
                readSites[write++] = readSites[begin];
            } else {
                out.droppedAmbiguousReadSiteCount++;
            }
            begin = end;
        }
        readSites.resize(write);

        sort(readSites.begin(), readSites.end(),
            [](const GlobalHetSiteReadIndex::ReadSite& a, const GlobalHetSiteReadIndex::ReadSite& b) {
                if (a.readPosition != b.readPosition) {
                    return a.readPosition < b.readPosition;
                }
                return a.siteId < b.siteId;
            });
    }

    return out;
}
