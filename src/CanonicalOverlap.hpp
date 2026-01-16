#ifndef DINARA_CANONICAL_OVERLAP_HPP
#define DINARA_CANONICAL_OVERLAP_HPP

// CanonicalOverlap: Per-Read overlap storage matching Hifiasm efficiency.
// 40 bytes per overlap, indexed by ReadId (not OrientedReadId).
// Canonical ordering: qn < tn always.

#include <cstdint>
#include <vector>
#include <algorithm>

namespace dinara {

// Compact overlap representation (40 bytes)
// Combines essential fields from AlignmentData + AlignmentInfo
struct CanonicalOverlap {
    // Read IDs (8 bytes) - Canonical: qn < tn
    uint32_t qn;            // Query ReadId
    uint32_t tn;            // Target ReadId
    
    // Coordinates (16 bytes)
    uint32_t qs;            // Query Start
    uint32_t qe;            // Query End
    uint32_t ts;            // Target Start
    uint32_t te;            // Target End
    
    // Marker ordinal info (8 bytes) - Condensed AlignmentInfo
    uint16_t markerCount;   // Number of aligned markers
    uint16_t firstOrdinal0; // First marker ordinal for read0
    uint16_t lastOrdinal0;  // Last marker ordinal for read0
    uint16_t firstOrdinal1; // First marker ordinal for read1 (last = first + markerCount - 1 roughly)
    
    // Bitfields (4 bytes)
    uint32_t rev        : 1;  // 0=same strand, 1=reverse complement
    uint32_t del        : 1;  // Deleted by filtering
    uint32_t is_match   : 2;  // 0=unknown, 1=CIS, 2=TRANS
    uint32_t strong     : 1;  // Has informative SNP evidence
    uint32_t lg_indel   : 1;  // Has large indel
    uint32_t cc         : 26; // Chain count / alignment length
    
    // Reserved/padding (4 bytes)
    uint32_t cigarIdx;      // Index into separate phasingCigars VectorOfVectors
    
    // === Helper Methods ===
    
    // Phasing status
    bool isUnknown() const { return is_match == 0; }
    bool isCis() const { return is_match == 1; }
    bool isTrans() const { return is_match == 2; }
    
    // Deletion status (by filtering)
    bool isDeleted() const { return del == 1; }
    void setDeleted(bool v) { del = v ? 1 : 0; }
    
    // Phasing classification
    void markCis() { is_match = 1; }
    void markTrans() { is_match = 2; }
    void markStrong() { strong = 1; }
    
    // Coordinate helpers
    uint32_t queryLength() const { return qe - qs; }
    uint32_t targetLength() const { return te - ts; }
    
    // Left/right unaligned tails (for containment detection)
    // Note: Requires read lengths to be passed in
    uint32_t leftTailQ(uint32_t qLen) const { return qs; }
    uint32_t rightTailQ(uint32_t qLen) const { return qLen - qe; }
    uint32_t leftTailT(uint32_t tLen) const { return ts; }
    uint32_t rightTailT(uint32_t tLen) const { return tLen - te; }
};

static_assert(sizeof(CanonicalOverlap) == 40, "CanonicalOverlap must be 40 bytes");


// Per-Read overlap index (Hifiasm sources/reverse_sources equivalent)
class OverlapIndex {
public:
    // sources[qn] = overlaps where this read is query (canonical: qn < tn)
    std::vector<std::vector<CanonicalOverlap>> sources;
    
    // reverseIndices[tn] = {(qn, idx)} pairs for overlaps where tn is target
    // This allows finding overlaps involving tn without duplicating data
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> reverseIndices;
    
    // Total overlap count
    uint64_t totalOverlaps = 0;
    
    // Initialize for N reads
    void resize(uint32_t readCount) {
        sources.resize(readCount);
        reverseIndices.resize(readCount);
        totalOverlaps = 0;
    }
    
    // Clear all data
    void clear() {
        for (auto& v : sources) v.clear();
        for (auto& v : reverseIndices) v.clear();
        totalOverlaps = 0;
    }
    
    // Add an overlap (ensures canonical ordering: qn < tn)
    void addOverlap(CanonicalOverlap ov) {
        // Ensure canonical ordering
        if (ov.qn > ov.tn) {
            // Swap to make qn < tn
            std::swap(ov.qn, ov.tn);
            std::swap(ov.qs, ov.ts);
            std::swap(ov.qe, ov.te);
            std::swap(ov.firstOrdinal0, ov.firstOrdinal1);
            // Note: lastOrdinal1 is not stored, but markerCount is same
        }
        
        uint32_t idx = (uint32_t)sources[ov.qn].size();
        sources[ov.qn].push_back(ov);
        reverseIndices[ov.tn].push_back({ov.qn, idx});
        totalOverlaps++;
    }
    
    // Get overlap by (qn, index)
    CanonicalOverlap& getOverlap(uint32_t qn, uint32_t idx) {
        return sources[qn][idx];
    }
    const CanonicalOverlap& getOverlap(uint32_t qn, uint32_t idx) const {
        return sources[qn][idx];
    }
    
    // Get all overlaps where readId is QUERY (qn)
    std::vector<CanonicalOverlap>& getOverlapsAsQuery(uint32_t readId) {
        return sources[readId];
    }
    const std::vector<CanonicalOverlap>& getOverlapsAsQuery(uint32_t readId) const {
        return sources[readId];
    }
    
    // Iterate overlaps where readId is TARGET (tn)
    // Callback: void(CanonicalOverlap&)
    template<typename Func>
    void forEachOverlapAsTarget(uint32_t readId, Func&& func) {
        for (auto& [qn, idx] : reverseIndices[readId]) {
            func(sources[qn][idx]);
        }
    }
    
    // Get count of overlaps for a read (as query + as target)
    uint64_t getOverlapCount(uint32_t readId) const {
        return sources[readId].size() + reverseIndices[readId].size();
    }
    
    // Sort overlaps by target for each source (useful for some algorithms)
    void sortByTarget() {
        for (auto& vec : sources) {
            std::sort(vec.begin(), vec.end(), 
                [](const CanonicalOverlap& a, const CanonicalOverlap& b) {
                    return a.tn < b.tn;
                });
        }
        // Rebuild reverse indices after sorting
        for (auto& v : reverseIndices) v.clear();
        for (uint32_t qn = 0; qn < sources.size(); qn++) {
            for (uint32_t idx = 0; idx < sources[qn].size(); idx++) {
                uint32_t tn = sources[qn][idx].tn;
                reverseIndices[tn].push_back({qn, idx});
            }
        }
    }
    
    // Count non-deleted overlaps
    uint64_t countActive() const {
        uint64_t count = 0;
        for (const auto& vec : sources) {
            for (const auto& ov : vec) {
                if (!ov.isDeleted()) count++;
            }
        }
        return count;
    }
    
    // Count TRANS overlaps
    uint64_t countTrans() const {
        uint64_t count = 0;
        for (const auto& vec : sources) {
            for (const auto& ov : vec) {
                if (ov.isTrans()) count++;
            }
        }
        return count;
    }
};

} // namespace dinara

#endif // DINARA_CANONICAL_OVERLAP_HPP
