#pragma once
// Internal types and function declarations shared between AssemblerHifiasmEC.cpp
// and AssemblerGlobalSiteEC.cpp.  NOT part of the public API.

#include <array>
#include <cstdint>
#include <vector>

namespace dinara {

class Assembler;

// ---------------------------------------------------------------------------
// Candidate overlap descriptor for a single query read.
// ---------------------------------------------------------------------------
struct CandidateEC {
    uint32_t alignmentId;
    uint64_t qs, qe, ts, te;
    uint32_t targetId;
    bool isRev;
    uint8_t is_match = 1;
    uint8_t strong = 0;
};

// ---------------------------------------------------------------------------
// Per-allele statistics for a single SNP row.
// ---------------------------------------------------------------------------
struct SnpStats {
    uint32_t site;
    uint32_t occ_0;
    uint32_t occ_1;
    uint32_t fwd_ref_cov;
    char refBase;
    char altBase;
    uint8_t is_homopolymer;
    int score;
    int dpScore;

    SnpStats() {
        site = uint32_t(-1);
        occ_0 = occ_1 = 0;
        fwd_ref_cov = 0;
        refBase = altBase = 0;
        is_homopolymer = 0;
        score = 0;
        dpScore = 0;
    }
};

// ---------------------------------------------------------------------------
// Per-overlap evidence entry for a single site.
// ---------------------------------------------------------------------------
struct HaplotypeEvidence {
    uint32_t overlapID;
    uint32_t site;
    uint32_t overlapSite;
    uint8_t type;
    uint8_t misBase;
    bool hp;

    bool operator<(const HaplotypeEvidence& other) const {
        if (site != other.site) return site < other.site;
        return overlapID < other.overlapID;
    }
};

// ---------------------------------------------------------------------------
// Sweep-line event for coverage computation.
// ---------------------------------------------------------------------------
struct SweepEvent {
    size_t siteIdx;
    uint32_t candIdx;
    bool isEnd;
    bool operator<(const SweepEvent& other) const {
        if (siteIdx != other.siteIdx) return siteIdx < other.siteIdx;
        return isEnd < other.isEnd;
    }
};

// ---------------------------------------------------------------------------
struct RawSV {
    uint32_t overlapID;
    uint32_t site;
    int64_t size;
};

// ---------------------------------------------------------------------------
struct RphaseDpTiming {
    double buildAltAnyBits = 0.;
    double buildRefBits = 0.;
    double buildHpBits = 0.;
    double transitions = 0.;
    double extractPaths = 0.;
};

// ---------------------------------------------------------------------------
// Thread-local scratchpad for the parity EC pipeline.
// ---------------------------------------------------------------------------
struct HifiasmECScratchPad {
    struct GapInterval {
        uint32_t begin;
        uint32_t end;
    };

    std::vector<CandidateEC> candidates;
    std::vector<SnpStats> snpStats;
    std::vector<HaplotypeEvidence> hapEvidence;
    std::vector<SnpStats> svStats;
    std::vector<HaplotypeEvidence> svEvidence;

    std::vector<uint32_t> uniqueSites;
    std::vector<int32_t> diffTotal;
    std::vector<int32_t> diffFwd;
    std::vector<uint32_t> siteTotalCov;
    std::vector<uint32_t> siteFwdCov;

    std::vector<uint32_t> insertionOffsets;
    std::vector<GapInterval> insertionIntervals;
    std::vector<uint32_t> insertionBaseCount;

    std::vector<int> validIndices;
    std::vector<int64_t> f;
    std::vector<int> p;
    std::vector<int> indexMap;
    std::vector<uint64_t> flatBits;
    std::vector<uint64_t> flatAnyBits;
    std::vector<uint64_t> flatHpBits;
    std::vector<SweepEvent> events;
    std::vector<SweepEvent> gapEvents;
    std::vector<uint64_t> active;
    std::vector<uint64_t> gapped;

    std::vector<RawSV> rawSVs;
    std::vector<size_t> svIndices;
    std::vector<uint8_t> unpackedRead;
    std::vector<uint8_t> covered;
    std::vector<uint32_t> path;
    std::vector<uint64_t> supportBits;
    std::vector<uint64_t> conflictBits;

    void clear() {
        candidates.clear();
        snpStats.clear();
        hapEvidence.clear();
        svStats.clear();
        svEvidence.clear();
        uniqueSites.clear();
        diffTotal.clear();
        diffFwd.clear();
        siteTotalCov.clear();
        siteFwdCov.clear();
        insertionOffsets.clear();
        insertionIntervals.clear();
        insertionBaseCount.clear();
        validIndices.clear();
        f.clear();
        p.clear();
        indexMap.clear();
        flatBits.clear();
        flatAnyBits.clear();
        flatHpBits.clear();
        events.clear();
        gapEvents.clear();
        active.clear();
        gapped.clear();
        rawSVs.clear();
        svIndices.clear();
        unpackedRead.clear();
        covered.clear();
        path.clear();
        supportBits.clear();
        conflictBits.clear();
    }
};

} // namespace dinara

// ---------------------------------------------------------------------------
// Functions implemented in AssemblerHifiasmEC.cpp.
// These live at global scope in that TU (via "using namespace dinara;").
// ---------------------------------------------------------------------------
void gen_rphase_dp(
    dinara::Assembler& assembler,
    dinara::HifiasmECScratchPad& scratch,
    dinara::RphaseDpTiming* timing = nullptr);

void generate_haplotypes_naive_HiFi(
    dinara::Assembler& assembler,
    dinara::HifiasmECScratchPad& scratch);
