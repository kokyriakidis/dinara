/**
 * @file test_integration.cpp
 * @brief Full integration tests for marker generation and filtering using real Assembler.
 * 
 * These tests verify the complete dinara marker pipeline:
 * 1. Read loading from FASTQ
 * 2. Marker generation using SIMD closed syncmers
 * 3. K-mer counting from marker KmerIds
 * 4. Marker filtering by frequency
 */

#include "../external/catch2/catch.hpp"

// Dinara headers
#include "../src/Assembler.hpp"
#include "../src/Reads.hpp"
#include "../src/KmerCounter.hpp"
#include "../src/ProjectedAlignment.hpp"
#include "../src/DINARA_ASSERT.hpp"
#include "../src/markerAccessFunctions.hpp"

// Standard library
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <string>
#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>
#include <tuple>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

using namespace dinara;

namespace {
class ScopedSilenceIo {
public:
    ScopedSilenceIo()
        : nullStream("/dev/null")
        , oldCoutBuf(std::cout.rdbuf())
        , oldCerrBuf(std::cerr.rdbuf())
    {
        std::cout.rdbuf(nullStream.rdbuf());
        std::cerr.rdbuf(nullStream.rdbuf());
    }

    ~ScopedSilenceIo()
    {
        std::cout.rdbuf(oldCoutBuf);
        std::cerr.rdbuf(oldCerrBuf);
    }

    ScopedSilenceIo(const ScopedSilenceIo&) = delete;
    ScopedSilenceIo& operator=(const ScopedSilenceIo&) = delete;

private:
    std::ofstream nullStream;
    std::streambuf* oldCoutBuf;
    std::streambuf* oldCerrBuf;
};

template<class F>
decltype(auto) withSilencedIo(F&& f)
{
    ScopedSilenceIo silence;
    return std::forward<F>(f)();
}

class ScopedCurrentPath {
public:
    explicit ScopedCurrentPath(const fs::path& p) : old(fs::current_path())
    {
        fs::current_path(p);
    }
    ~ScopedCurrentPath() { fs::current_path(old); }

    ScopedCurrentPath(const ScopedCurrentPath&) = delete;
    ScopedCurrentPath& operator=(const ScopedCurrentPath&) = delete;

private:
    fs::path old;
};

template<class F>
decltype(auto) withSilencedIoInDir(const fs::path& dir, F&& f)
{
    ScopedCurrentPath scoped(dir);
    return withSilencedIo(std::forward<F>(f));
}

fs::path makeUniqueTempDir(const std::string& prefix)
{
    static std::atomic<uint64_t> counter{0};
    const auto base = fs::temp_directory_path();

    for (int attempt = 0; attempt < 100; ++attempt) {
        const auto suffix =
            std::to_string(::getpid()) + "_" +
            std::to_string(counter.fetch_add(1)) + "_" +
            std::to_string(uint64_t(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::path dir = base / (prefix + suffix);
        std::error_code ec;
        if (fs::create_directory(dir, ec)) {
            return dir;
        }
    }
    throw std::runtime_error("Failed to create unique temp directory under " + base.string());
}

std::string toString(LongBaseSequenceView view)
{
    std::string s;
    s.resize(view.baseCount);
    for (size_t i = 0; i < view.baseCount; ++i) {
        static const char map[] = {'A', 'C', 'G', 'T'};
        s[i] = map[view[i].value & 3U];
    }
    return s;
}

template<class AlignmentDataContainer>
uint32_t findAlignmentEvidenceId(
    const AlignmentDataContainer& alignmentData,
    ReadId a,
    ReadId b,
    std::optional<bool> expectedSameStrand = std::nullopt)
{
    for (const auto& ad : alignmentData) {
        if (ad.isDeleted()) continue;
        const bool matchesPair =
            (ad.readIds[0] == a && ad.readIds[1] == b) ||
            (ad.readIds[0] == b && ad.readIds[1] == a);
        if (!matchesPair) continue;
        if (expectedSameStrand && ad.isSameStrand != *expectedSameStrand) continue;
        return uint32_t(ad.info.alignmentId);
    }
    return invalid<uint32_t>;
}

template<class AlignmentDataContainer>
const AlignmentData* findAlignmentDataPtr(
    const AlignmentDataContainer& alignmentData,
    ReadId a,
    ReadId b,
    std::optional<bool> expectedSameStrand = std::nullopt)
{
    for (const auto& ad : alignmentData) {
        if (ad.isDeleted()) continue;
        const bool matchesPair =
            (ad.readIds[0] == a && ad.readIds[1] == b) ||
            (ad.readIds[0] == b && ad.readIds[1] == a);
        if (!matchesPair) continue;
        if (expectedSameStrand && ad.isSameStrand != *expectedSameStrand) continue;
        return &ad;
    }
    return nullptr;
}
} // namespace

// =============================================================================
// HELPER: Generate random DNA sequence
// =============================================================================
std::string randomSequence(size_t length, uint32_t seed = 42) {
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dis(0, 3);
    std::string seq(length, 'N');
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (size_t i = 0; i < length; i++) {
        seq[i] = bases[dis(gen)];
    }
    return seq;
}

// =============================================================================
// FIXTURE: Creates temp directory with test FASTQ and initializes Assembler
// =============================================================================
class AssemblerIntegrationFixture {
public:
    fs::path testDir;
    fs::path fastqPath;
    std::unique_ptr<Assembler> assembler;
    
    AssemblerIntegrationFixture() {
        testDir = makeUniqueTempDir("dinara_integration_");
        fastqPath = testDir / "test_reads.fastq";
    }
    
    ~AssemblerIntegrationFixture() {
        assembler.reset();
        try {
            fs::remove_all(testDir);
        } catch (...) {}
    }
    
    void createFastq(const std::vector<std::string>& sequences) {
        std::ofstream out(fastqPath);
        if (!out) {
            throw std::runtime_error("Failed to open " + fastqPath.string());
        }
        for (size_t i = 0; i < sequences.size(); i++) {
            out << "@read_" << i << "\n";
            out << sequences[i] << "\n";
            out << "+\n";
            out << std::string(sequences[i].size(), '~') << "\n";
        }
    }
    
    void initAssembler() {
        withSilencedIoInDir(testDir, [&] {
            assembler = std::make_unique<Assembler>(
                testDir.string() + "/",
                true,   // createNew
                0,      // readRepresentation = raw
                4096    // pageSize
            );
        });
    }
    
    void loadReads(uint64_t minReadLength = 0) {
        withSilencedIoInDir(testDir, [&] {
            assembler->addReads(fastqPath.string(), minReadLength, true, 1);
            assembler->computeReadIdsSortedByName();
        });
    }
    
    void generateMarkers(int k = 16, int s = 4) {
        withSilencedIoInDir(testDir, [&] { assembler->findMarkersSimdClosedSyncmers(1, k, s); });
    }
    
    void countKmers() {
        withSilencedIoInDir(testDir, [&] { assembler->countKmersFromMarkerKmerIds(1); });
    }
    
    void applyFilter(uint64_t minFreq, uint64_t maxFreq) {
        withSilencedIoInDir(testDir, [&] { assembler->applyKmerCountFilter(minFreq, maxFreq, 1); });
    }
    
    void findCandidates() {
        withSilencedIoInDir(testDir, [&] { assembler->findAlignmentCandidatesInvertedIndex(0.1, 100, 1); });
    }

    // Granular pipeline for testing
    void buildIndex() {
        withSilencedIoInDir(testDir, [&] { assembler->buildInvertedIndex(1); });
    }

    void chainCandidates(double maxDriftRate = 0.1, uint64_t maxChainLimit = 100) {
        withSilencedIoInDir(testDir, [&] { assembler->chainAlignmentCandidates(maxDriftRate, maxChainLimit, 1); });
    }
    
    void computeAlignments() {
        AlignOptions options;
        // Initialize all required fields for projected alignment on precomputed chains
        options.alignMethod = 6;
        options.maxSkip = 100;
        options.maxDrift = 100;
        options.maxTrim = 10000;
        options.minAlignedMarkerCount = 4;
        options.minAlignedFraction = 0.0;
        options.maxMarkerFrequency = 1000;
        options.matchScore = 3;
        options.mismatchScore = -1;
        options.gapScore = -1;
        options.downsamplingFactor = 0.1;
        options.bandExtend = 10;
        options.maxBand = 1000;
        options.sameChannelReadAlignmentSuppressDeltaThreshold = 0;
        options.suppressContainments = false;
        options.align4DeltaX = 200;
        options.align4DeltaY = 10;
        options.align4MinEntryCountPerCell = 10;
        options.align4MaxDistanceFromBoundary = 100;
        options.align5DriftRateTolerance = 0.02;
        options.align5MinBandExtend = 10;
        options.maxErrorRate = 0.3; // Higher for synthetic tests
        withSilencedIoInDir(testDir, [&] { assembler->computeAlignmentsWithEvidence(options, 1); });
    }

    struct PafOverlapSpec {
        std::string qName;
        std::string tName;
        bool sameStrand;
        uint32_t qStart;
        uint32_t qEnd;
        uint32_t tStart;
        uint32_t tEnd;
    };

    // Create a PAF file with overlap information
    std::string createPafFile(const std::vector<PafOverlapSpec>& overlaps) {
        std::string pafPath = (testDir / "overlaps.paf").string();
        std::ofstream out(pafPath);
        if (!out) {
            throw std::runtime_error("Failed to open " + pafPath);
        }

        const auto& reads = assembler->getReads();

        // PAF format: qName qLen qStart qEnd strand tName tLen tStart tEnd matches alignLen mapQ
        for (const auto& ov : overlaps) {
            const ReadId qId = reads.getReadId(ov.qName);
            const ReadId tId = reads.getReadId(ov.tName);
            const uint64_t qLen = reads.getReadRawSequenceLength(qId);
            const uint64_t tLen = reads.getReadRawSequenceLength(tId);

            if (ov.qEnd > qLen || ov.tEnd > tLen || ov.qStart >= ov.qEnd || ov.tStart >= ov.tEnd) {
                throw std::runtime_error("Invalid PAF overlap coordinates for " + ov.qName + " vs " + ov.tName);
            }

            const uint32_t qSpan = ov.qEnd - ov.qStart;
            const uint32_t tSpan = ov.tEnd - ov.tStart;
            const uint32_t alignLen = std::min(qSpan, tSpan);
            const uint32_t matches = alignLen; // synthetic exact overlap

            out << ov.qName << "\t" << qLen << "\t" << ov.qStart << "\t" << ov.qEnd << "\t"
                << (ov.sameStrand ? "+" : "-") << "\t"
                << ov.tName << "\t" << tLen << "\t" << ov.tStart << "\t" << ov.tEnd << "\t"
                << matches << "\t" << alignLen << "\t60\n";
        }
        return pafPath;
    }

    void importPafCandidates(const std::string& pafPath) {
        withSilencedIoInDir(testDir, [&] { assembler->importAlignmentCandidatesFromPaf(pafPath); });
    }

    void chainPafCandidates(double maxDriftRate = 0.1, uint64_t maxChainLimit = 100) {
        withSilencedIoInDir(testDir, [&] { assembler->chainPafCandidates(maxDriftRate, maxChainLimit, 1); });
    }
};

// =============================================================================
// HELPER: Reverse complement sequence
// =============================================================================
std::string reverseComplement(const std::string& seq) {
    std::string rc = seq;
    std::reverse(rc.begin(), rc.end());
    for (auto& b : rc) {
        if (b == 'A') b = 'T';
        else if (b == 'C') b = 'G';
        else if (b == 'G') b = 'C';
        else if (b == 'T') b = 'A';
    }
    return rc;
}

namespace dinara::hifiasmEcTestHooks {
    int64_t scoreLinkStrict(
        const std::vector<uint64_t>& flatBits,
        const std::vector<uint64_t>& flatAnyBits,
        size_t siteI,
        size_t siteJ,
        size_t nWords);
    void resetDpSameSiteComparisons();
    uint64_t getDpSameSiteComparisons();
}

// =============================================================================
// HELPER: Toggle base (for SNPs)
// =============================================================================
char otherBase(char b) {
    if (b == 'A') return 'C';
    return 'A';
}

// =============================================================================
// TEST: Basic Read Loading
// =============================================================================
TEST_CASE("Integration: Assembler read loading", "[integration][reads]") {
    AssemblerIntegrationFixture fixture;
    
    std::string seq1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";  // 64 bases
    std::string seq2 = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";  // 64 bases
    
    fixture.createFastq({seq1, seq2});
    fixture.initAssembler();
    fixture.loadReads();
    
    const auto& reads = fixture.assembler->getReads();
    REQUIRE(reads.readCount() == 2);
    CHECK(toString(reads.getRead(0)) == seq1);
    CHECK(toString(reads.getRead(1)) == seq2);
    CHECK(std::string(reads.getReadName(0).begin(), reads.getReadName(0).end()) == "read_0");
    CHECK(std::string(reads.getReadName(1).begin(), reads.getReadName(1).end()) == "read_1");
}

// =============================================================================
// TEST: Marker Generation
// =============================================================================
TEST_CASE("Integration: Marker generation from syncmers", "[integration][markers]") {
    AssemblerIntegrationFixture fixture;
    
    // 100 base sequence - should produce several markers
    std::string seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
                      "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";
    
    fixture.createFastq({seq});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(15, 5);
    
    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    // Read 0 strand 0
    const auto read0Strand0 = OrientedReadId(0, 0).getValue();
    const auto read0Strand1 = OrientedReadId(0, 1).getValue();
    const size_t count0 = markers->size(read0Strand0);
    const size_t count1 = markers->size(read0Strand1);
    REQUIRE(count0 > 0);
    CHECK(count0 == count1);

    // Markers are sorted and within bounds.
    const auto readMarkers = (*markers)[read0Strand0];
    for (size_t i = 0; i < readMarkers.size(); ++i) {
        CHECK(uint32_t(readMarkers[i].position) + 15 <= seq.size());
        if (i > 0) {
            CHECK(uint32_t(readMarkers[i - 1].position) <= uint32_t(readMarkers[i].position));
        }
    }
}

// =============================================================================
// TEST: K-mer Counting
// =============================================================================
TEST_CASE("Integration: K-mer counting", "[integration][counting]") {
    AssemblerIntegrationFixture fixture;
    
    // Create reads with varied k-mers to satisfy getHistogramHistogram logic
    // We need more repeated k-mers than unique ones to create a valley in the histogram.
    std::string seqSH = randomSequence(200, 123);
    std::string seqLG = randomSequence(1000, 456); 
    
    fixture.createFastq({seqLG, seqLG, seqSH});  // 2x Long (Freq 2), 1x Short (Freq 1)
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(15, 5);
    fixture.countKmers();
    
    auto& kmerCounter = *(fixture.assembler->kmerCounter);
    REQUIRE(kmerCounter.kmerIdFrequencies.size() > 0);

    bool foundFreq1 = false;
    bool foundFreq2Plus = false;
    uint64_t maxFreq = 0;
    for (auto& kv : kmerCounter.kmerIdFrequencies) {
        const uint64_t freq = kv.second;
        maxFreq = std::max(maxFreq, freq);
        if (freq == 1) foundFreq1 = true;
        if (freq >= 2) foundFreq2Plus = true;
    }
    CHECK(foundFreq1);      // seqSH contributes rare markers
    CHECK(foundFreq2Plus);  // seqLG duplicated
    CHECK(maxFreq >= 2);
}

// =============================================================================
// TEST: Marker Filtering
// =============================================================================
TEST_CASE("Integration: Marker filtering by frequency", "[integration][filtering]") {
    AssemblerIntegrationFixture fixture;
    
    // Create reads with varied k-mers to satisfy getHistogramInfo valley search
    // seqLG (repeated twice) creates Freq(2). seqSH (repeated once) creates Freq(1).
    // Ensure Freq(2) > Freq(1) by making seqLG much longer.
    std::string seqSH = randomSequence(200, 789);
    std::string seqLG = randomSequence(1000, 101112);
    
    fixture.createFastq({seqLG, seqLG, seqSH});
    fixture.initAssembler();
    fixture.loadReads();
    const uint64_t k = 15;
    fixture.generateMarkers(int(k), 5);
    fixture.countKmers();
    
    auto* markersBeforePtr = fixture.assembler->markers.get();
    REQUIRE(markersBeforePtr != nullptr);

    const uint64_t markersBefore = markersBeforePtr->totalSize();
    const auto read0s0 = markersBeforePtr->size(OrientedReadId(0, 0).getValue());
    const auto read2s0 = markersBeforePtr->size(OrientedReadId(2, 0).getValue());
    
    // Filter: keep only k-mers with frequency >= 2
    fixture.applyFilter(2, 1000);
    
    auto* markersAfterPtr = fixture.assembler->markers.get();
    REQUIRE(markersAfterPtr != nullptr);
    const uint64_t markersAfter = markersAfterPtr->totalSize();
    const auto read0s0After = markersAfterPtr->size(OrientedReadId(0, 0).getValue());
    const auto read2s0After = markersAfterPtr->size(OrientedReadId(2, 0).getValue());
    
    CHECK(markersAfter <= markersBefore);
    REQUIRE(markersAfter > 0);
    REQUIRE(read0s0After > 0);
    CHECK(read0s0After <= read0s0);

    // The unique read (read_2) should lose essentially all its markers under minFreq=2.
    CHECK(read2s0After == 0);

    // Every remaining marker must have k-mer frequency within the requested bounds.
    const auto& reads = fixture.assembler->getReads();
    auto& kmerCounter = *(fixture.assembler->kmerCounter);
    kmerCounter.buildFrequencyLUT();

    for (uint32_t orientedReadIdValue = 0; orientedReadIdValue < markersAfterPtr->size(); ++orientedReadIdValue) {
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
        const auto orientedMarkers = (*markersAfterPtr)[orientedReadIdValue];
        for (uint32_t ordinal = 0; ordinal < orientedMarkers.size(); ++ordinal) {
            const Kmer kmer = dinara::getOrientedReadMarkerKmer(
                orientedReadId, ordinal, k, reads, *markersAfterPtr);
            const KmerId id = KmerId(kmer.id(k));
            const KmerId rc = KmerId(kmer.reverseComplement(k).id(k));
            const KmerId canonical = std::min(id, rc);
            const uint64_t freq = kmerCounter.getFrequencyFast(canonical);
            CHECK(freq >= 2);
            CHECK(freq <= 1000);
        }
    }
}

// =============================================================================
// TEST: Forward/Reverse Complement Symmetry
// =============================================================================
TEST_CASE("Integration: F/RC marker symmetry", "[integration][symmetry]") {
    AssemblerIntegrationFixture fixture;
    
    std::string seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    
    fixture.createFastq({seq});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(15, 5);
    
    auto* markers = fixture.assembler->markers.get();
    
    size_t strand0 = markers->size(OrientedReadId(0, 0).getValue());
    size_t strand1 = markers->size(OrientedReadId(0, 1).getValue());
    
    CHECK(strand0 == strand1);
}
// =============================================================================
// TEST: Projected Alignment and Evidence Storage
// =============================================================================
TEST_CASE("Integration: Projected alignment and evidence storage", "[integration][evidence]") {
    AssemblerIntegrationFixture fixture;
    
    // Create base sequence [prefix][variant_region][suffix]
    std::string prefix = randomSequence(1000, 111);
    std::string mid = randomSequence(1000, 222);
    std::string suffix = randomSequence(1000, 333);
    
    // Read 0: Pure sequence [prefix][mid][suffix]
    std::string seq0 = prefix + mid + suffix;
    
    // Read 1: SNP at mid[500] (Target-projection view)
    std::string seq1 = seq0;
    char targetBase = seq1[1500];
    seq1[1500] = otherBase(targetBase);
    
    // Read 2: 50bp deletion in mid
    std::string seq2 = prefix + mid.substr(0, 475) + mid.substr(525) + suffix;
    
    // Read 3: RC of seq1 (Test F-R orientation)
    std::string seq3 = reverseComplement(seq1);
    
    fixture.createFastq({seq0, seq1, seq2, seq3});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000); 

    fixture.findCandidates();
    fixture.computeAlignments();

    const auto& store = fixture.assembler->alignedEvidenceStore;
    const auto& reads = fixture.assembler->getReads();
    const auto& alignmentData = fixture.assembler->alignmentData;
    CAPTURE(store.index.size());
    CAPTURE(alignmentData.size());
    REQUIRE(store.index.size() > 0);
    REQUIRE(alignmentData.size() > 0);

    // --- SNP (read_0 vs read_1, same strand) ---
    const AlignmentData* ad01 = findAlignmentDataPtr(alignmentData, ReadId(0), ReadId(1), true);
    REQUIRE(ad01 != nullptr);
    const uint32_t ev01 = uint32_t(ad01->info.alignmentId);
    const OrientedReadId q01(ad01->readIds[0], 0);
    const OrientedReadId t01(ad01->readIds[1], ad01->isSameStrand ? 0 : 1);
    const uint8_t expectedQ01 = reads.getOrientedReadBase(q01, 1500).value;
    const uint8_t expectedT01 = reads.getOrientedReadBase(t01, 1500).value;

    std::vector<std::pair<uint32_t, uint8_t>> s0_01;
    store.forEachSnp0InRange(ev01, 1490, 1510, [&](uint32_t pos, uint8_t base) { s0_01.push_back({pos, base}); });
    std::vector<std::pair<uint32_t, uint8_t>> s1_01;
    store.forEachSnp1InRange(ev01, 1490, 1510, [&](uint32_t pos, uint8_t base) { s1_01.push_back({pos, base}); });
    REQUIRE(s0_01.size() == 1);
    REQUIRE(s1_01.size() == 1);
    CHECK(s0_01[0].first == 1500);
    CHECK(s1_01[0].first == 1500);
    CHECK(s0_01[0].second == expectedQ01);
    CHECK(s1_01[0].second == expectedT01);

    // --- Deletion (read_0 vs read_2) ---
    const uint32_t ev02 = findAlignmentEvidenceId(alignmentData, ReadId(0), ReadId(2));
    REQUIRE(ev02 != invalid<uint32_t>);
    const auto indels0 = store.getIndels0(ev02);
    const auto indels1 = store.getIndels1(ev02);
    REQUIRE_FALSE((indels0.empty() && indels1.empty()));

    // The injected deletion was 50bp starting at read_0 coordinate 1000+475=1475.
    uint32_t sumInWindow0 = 0;
    uint32_t sumInWindow1 = 0;
    for (const auto& ev : indels0) if (ev.pos() >= 1450 && ev.pos() <= 1550) sumInWindow0 += ev.len();
    for (const auto& ev : indels1) if (ev.pos() >= 1450 && ev.pos() <= 1550) sumInWindow1 += ev.len();
    CHECK(std::max(sumInWindow0, sumInWindow1) >= 45);

    // --- F-R orientation (read_0 vs read_3, different strands) ---
    const AlignmentData* ad03 = nullptr;
    for (const auto& ad : alignmentData) {
        if (ad.isDeleted()) continue;
        if (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(3) && ad.isSameStrand == false) {
            ad03 = &ad;
            break;
        }
    }
    if (!ad03) {
        ad03 = findAlignmentDataPtr(alignmentData, ReadId(0), ReadId(3), false);
    }
    REQUIRE(ad03 != nullptr);
    CHECK(ad03->isSameStrand == false);
    const uint32_t ev03 = uint32_t(ad03->info.alignmentId);
    const OrientedReadId q03(ad03->readIds[0], 0);
    const OrientedReadId t03(ad03->readIds[1], ad03->isSameStrand ? 0 : 1);
    const uint8_t expectedQ03 = reads.getOrientedReadBase(q03, 1500).value;
    const uint8_t expectedT03 = reads.getOrientedReadBase(t03, 1500).value;
    static const uint8_t complementBase[4] = {3, 2, 1, 0};
    const uint32_t tRawLen03 = uint32_t(reads.getReadRawSequenceLength(ad03->readIds[1]));
    const uint32_t expectedS0Pos03 = tRawLen03 - 1U - 1500U;

    std::vector<std::pair<uint32_t, uint8_t>> s0_03;
    store.forEachSnp0InRange(ev03, 1490, 1510, [&](uint32_t pos, uint8_t base) { s0_03.push_back({pos, base}); });
    std::vector<std::pair<uint32_t, uint8_t>> s1_03;
    store.forEachSnp1InRange(ev03, 1490, 1510, [&](uint32_t pos, uint8_t base) { s1_03.push_back({pos, base}); });
    REQUIRE(s0_03.size() == 1);
    REQUIRE(s1_03.size() == 1);
    CHECK(s0_03[0].first == expectedS0Pos03);
    CHECK(s1_03[0].first == 1500);
    CHECK(s0_03[0].second == complementBase[expectedQ03]);
    CHECK(s1_03[0].second == expectedT03);
}

TEST_CASE("Integration: Overlap-event anchors select vertices uniquely", "[integration][anchors][events]") {
    AssemblerIntegrationFixture fixture;

    // Two reads with enough markers. We will manually construct:
    // - a single read-graph overlap (and its RC twin)
    // - a minimal marker graph with two canonical vertices (and their RC vertices)
    // such that the overlap start/end events map to those vertices.
    fixture.createFastq({randomSequence(4000, 501), randomSequence(4000, 502)});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);

    withSilencedIoInDir(fixture.testDir, [&] {
        // -----------------------
        // Alignment + ReadGraph
        // -----------------------
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        AlignmentInfo info;
        info.alignmentId = 0;
        info.isInReadGraph = 1;
        info.markerCount = lastOrdinal + 1 - firstOrdinal;
        info.data[0] = AlignmentInfo::Data(m0, firstOrdinal, lastOrdinal);
        info.data[1] = AlignmentInfo::Data(m1, firstOrdinal, lastOrdinal);
        info.minOrdinalOffset = 0;
        info.maxOrdinalOffset = 0;
        info.averageOrdinalOffset = 0;
        info.maxSkip = 0;
        info.maxDrift = 0;

        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);

        fixture.assembler->readGraph.edges.resize(2);
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0p, r1p};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[0] = e;
        }
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0m, r1m};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[1] = e;
        }

        // connectivity is indexed by orientedReadIdValue = readId*2 + strand.
        // Index 0: r0+, 1: r0-, 2: r1+, 3: r1-.
        auto& conn = fixture.assembler->readGraph.connectivity;
        conn.appendVector(); // 0
        conn.append(0);
        conn.appendVector(); // 1
        conn.append(1);
        conn.appendVector(); // 2
        conn.append(0);
        conn.appendVector(); // 3
        conn.append(1);

        // -----------------------
        // Minimal MarkerGraph
        // -----------------------
        auto& mg = fixture.assembler->markerGraph;
        mg.constructVertices();
        mg.vertices().createNew("", 4096);
        mg.vertexTable.createNew("", 4096);
        mg.vertexTable.resize(markers->totalSize());
        std::fill(mg.vertexTable.begin(), mg.vertexTable.end(), MarkerGraph::invalidCompressedVertexId);

        // Two canonical vertices and their reverse complements:
        //  vertex0 <-> vertex2, vertex1 <-> vertex3.
        const MarkerGraphVertexId v0 = 0;
        const MarkerGraphVertexId v1 = 1;
        const MarkerGraphVertexId v2 = 2;
        const MarkerGraphVertexId v3 = 3;

        mg.reverseComplementVertex.createNew("", 4096);
        mg.reverseComplementVertex.resize(4);
        mg.reverseComplementVertex[v0] = v2;
        mg.reverseComplementVertex[v2] = v0;
        mg.reverseComplementVertex[v1] = v3;
        mg.reverseComplementVertex[v3] = v1;

        const MarkerId r0p_first = fixture.assembler->getMarkerId(r0p, firstOrdinal);
        const MarkerId r1p_first = fixture.assembler->getMarkerId(r1p, firstOrdinal);
        const MarkerId r0p_after = fixture.assembler->getMarkerId(r0p, lastOrdinal + 1);
        const MarkerId r1p_after = fixture.assembler->getMarkerId(r1p, lastOrdinal + 1);

        const uint32_t rcFirstOrdinal0 = uint32_t(m0) - 1 - firstOrdinal;
        const uint32_t rcFirstOrdinal1 = uint32_t(m1) - 1 - firstOrdinal;
        const uint32_t rcAfterOrdinal0 = uint32_t(m0) - 1 - (lastOrdinal + 1);
        const uint32_t rcAfterOrdinal1 = uint32_t(m1) - 1 - (lastOrdinal + 1);

        const MarkerId r0m_first = fixture.assembler->getMarkerId(r0m, rcFirstOrdinal0);
        const MarkerId r1m_first = fixture.assembler->getMarkerId(r1m, rcFirstOrdinal1);
        const MarkerId r0m_after = fixture.assembler->getMarkerId(r0m, rcAfterOrdinal0);
        const MarkerId r1m_after = fixture.assembler->getMarkerId(r1m, rcAfterOrdinal1);

        // Vertex 0: overlap "start" event markers.
        {
            vector<MarkerId> ids = {r0p_first, r1p_first};
            std::sort(ids.begin(), ids.end());
            mg.vertices().appendVector();
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(v0));
            }
        }
        // Vertex 1: overlap "end" event markers (last+1).
        {
            vector<MarkerId> ids = {r0p_after, r1p_after};
            std::sort(ids.begin(), ids.end());
            mg.vertices().appendVector();
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(v1));
            }
        }
        // Vertex 2: RC of vertex 0.
        {
            vector<MarkerId> ids = {r0m_first, r1m_first};
            std::sort(ids.begin(), ids.end());
            mg.vertices().appendVector();
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(v2));
            }
        }
        // Vertex 3: RC of vertex 1.
        {
            vector<MarkerId> ids = {r0m_after, r1m_after};
            std::sort(ids.begin(), ids.end());
            mg.vertices().appendVector();
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(v3));
            }
        }
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromMarkerGraphVerticesAtOverlapEvents(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() == 4); // 2 canonical vertices * (forward+RC)

    // Verify that each (forward, rc) pair is consecutive and consistent.
    for(mode3::AnchorId a = 0; a < anchors->size(); a += 2) {
        const auto forward = (*anchors)[a];
        const auto reverse = (*anchors)[a + 1];
        REQUIRE(forward.size() == reverse.size());
        // The reverse-complement anchor is stored sorted by OrientedReadId (like all anchors).
        std::vector<dinara::mode3::AnchorMarkerInterval> expectedReverse;
        expectedReverse.reserve(forward.size());
        for(const auto& f : forward) {
            auto expected = f;
            expected.orientedReadId.flipStrand();
            const uint64_t markerCount = fixture.assembler->markers->size(expected.orientedReadId.getValue());
            expected.ordinal0 = uint32_t(markerCount) - 1 - f.ordinal0;
            expectedReverse.push_back(expected);
        }
        std::sort(expectedReverse.begin(), expectedReverse.end(),
            [](const dinara::mode3::AnchorMarkerInterval& x, const dinara::mode3::AnchorMarkerInterval& y) {
                return x.orientedReadId < y.orientedReadId;
            });

        for(size_t i = 0; i < reverse.size(); ++i) {
            CHECK(reverse[i].orientedReadId == expectedReverse[i].orientedReadId);
            CHECK(reverse[i].ordinal0 == expectedReverse[i].ordinal0);
        }
    }
}

TEST_CASE("Integration: Best-per-interval anchors pick the max-coverage vertex inside the interval", "[integration][anchors][events][best]") {
    AssemblerIntegrationFixture fixture;

    // Four reads. We manually construct readGraph overlaps so that read 0 has one active overlap interval,
    // and construct a marker graph where:
    // - The highest-coverage interior vertex contains both strands for read 0 (duplicate ReadId) and must be skipped.
    // - The next-best interior vertex is clean and should be selected.
    fixture.createFastq({
        randomSequence(4000, 601),
        randomSequence(4000, 602),
        randomSequence(4000, 603),
        randomSequence(4000, 604),
    });
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const OrientedReadId r2p(ReadId(2), 0);
    const OrientedReadId r2m(ReadId(2), 1);
    const OrientedReadId r3p(ReadId(3), 0);
    const OrientedReadId r3m(ReadId(3), 1);

    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    const uint32_t m2 = uint32_t(markers->size(r2p.getValue()));
    const uint32_t m3 = uint32_t(markers->size(r3p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);
    REQUIRE(m2 > 16);
    REQUIRE(m3 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);
    REQUIRE(lastOrdinal + 1 < m2);
    REQUIRE(lastOrdinal + 1 < m3);

    const uint32_t interiorBadOrdinal = 2;
    const uint32_t interiorGoodOrdinal = 3;
    REQUIRE(interiorBadOrdinal >= firstOrdinal);
    REQUIRE(interiorBadOrdinal < lastOrdinal + 1);
    REQUIRE(interiorGoodOrdinal >= firstOrdinal);
    REQUIRE(interiorGoodOrdinal < lastOrdinal + 1);

    withSilencedIoInDir(fixture.testDir, [&] {
        // -----------------------
        // AlignmentData: (0,1), (0,2), (0,3) all same-strand, same ordinal interval.
        // -----------------------
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        auto fillAlignment = [&](size_t idx, ReadId otherReadId) {
            AlignmentInfo info;
            info.alignmentId = idx;
            info.isInReadGraph = 1;
            info.markerCount = lastOrdinal + 1 - firstOrdinal;

            info.data[0] = AlignmentInfo::Data(m0, firstOrdinal, lastOrdinal);
            uint32_t otherMarkerCount = 0;
            if(otherReadId == ReadId(1)) otherMarkerCount = m1;
            if(otherReadId == ReadId(2)) otherMarkerCount = m2;
            if(otherReadId == ReadId(3)) otherMarkerCount = m3;
            info.data[1] = AlignmentInfo::Data(otherMarkerCount, firstOrdinal, lastOrdinal);

            info.minOrdinalOffset = 0;
            info.maxOrdinalOffset = 0;
            info.averageOrdinalOffset = 0;
            info.maxSkip = 0;
            info.maxDrift = 0;

            AlignmentData ad(OrientedReadPair(ReadId(0), otherReadId, true), info);
            fixture.assembler->alignmentData[idx] = ad;
        };
        fillAlignment(0, ReadId(1));
        fillAlignment(1, ReadId(2));
        fillAlignment(2, ReadId(3));

        // -----------------------
        // ReadGraph edges: each alignment yields an edge for + strand and an edge for - strand.
        // We only need connectivity for read 0 (both strands) but keep it consistent for endpoints too.
        // -----------------------
        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);
        fixture.assembler->readGraph.edges.resize(6);

        auto setEdge = [&](size_t edgeIndex, OrientedReadId a, OrientedReadId b, uint64_t alignmentId) {
            ReadGraphEdge e;
            e.orientedReadIds = {a, b};
            e.alignmentId = alignmentId;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[edgeIndex] = e;
        };
        setEdge(0, r0p, r1p, 0);
        setEdge(1, r0m, r1m, 0);
        setEdge(2, r0p, r2p, 1);
        setEdge(3, r0m, r2m, 1);
        setEdge(4, r0p, r3p, 2);
        setEdge(5, r0m, r3m, 2);

        auto& conn = fixture.assembler->readGraph.connectivity;
        // orientedReadCount = 8 (4 reads * 2 strands).
        // Build vectors in order 0..7.
        conn.appendVector(); // 0 r0+
        conn.append(0); conn.append(2); conn.append(4);
        conn.appendVector(); // 1 r0-
        conn.append(1); conn.append(3); conn.append(5);

        conn.appendVector(); // 2 r1+
        conn.append(0);
        conn.appendVector(); // 3 r1-
        conn.append(1);

        conn.appendVector(); // 4 r2+
        conn.append(2);
        conn.appendVector(); // 5 r2-
        conn.append(3);

        conn.appendVector(); // 6 r3+
        conn.append(4);
        conn.appendVector(); // 7 r3-
        conn.append(5);

        // -----------------------
        // Minimal MarkerGraph:
        // - boundary at firstOrdinal maps to low-coverage vertex v0 (2 markers: r0,r1)
        // - boundary at afterLast maps to low-coverage vertex v2 (2 markers: r0,r1)
        // - interior ordinal maps to:
        //     v4: highest coverage but *invalid* (contains r0+ and r0- -> duplicate ReadId)
        //     v6: next-best coverage and valid (unique ReadIds)
        // RC pairs are (0<->1), (2<->3), (4<->5), (6<->7).
        // -----------------------
        auto& mg = fixture.assembler->markerGraph;
        mg.constructVertices();
        mg.vertices().createNew("", 4096);
        mg.vertexTable.createNew("", 4096);
        mg.vertexTable.resize(markers->totalSize());
        std::fill(mg.vertexTable.begin(), mg.vertexTable.end(), MarkerGraph::invalidCompressedVertexId);

        mg.reverseComplementVertex.createNew("", 4096);
        mg.reverseComplementVertex.resize(8);
        mg.reverseComplementVertex[0] = 1;
        mg.reverseComplementVertex[1] = 0;
        mg.reverseComplementVertex[2] = 3;
        mg.reverseComplementVertex[3] = 2;
        mg.reverseComplementVertex[4] = 5;
        mg.reverseComplementVertex[5] = 4;
        mg.reverseComplementVertex[6] = 7;
        mg.reverseComplementVertex[7] = 6;

        const MarkerId r0_first = fixture.assembler->getMarkerId(r0p, firstOrdinal);
        const MarkerId r1_first = fixture.assembler->getMarkerId(r1p, firstOrdinal);
        const MarkerId r0_after = fixture.assembler->getMarkerId(r0p, lastOrdinal + 1);
        const MarkerId r1_after = fixture.assembler->getMarkerId(r1p, lastOrdinal + 1);

        const MarkerId r0_bad = fixture.assembler->getMarkerId(r0p, interiorBadOrdinal);
        const MarkerId r1_bad = fixture.assembler->getMarkerId(r1p, interiorBadOrdinal);
        const MarkerId r2_bad = fixture.assembler->getMarkerId(r2p, interiorBadOrdinal);
        const MarkerId r3_bad = fixture.assembler->getMarkerId(r3p, interiorBadOrdinal);

        const uint32_t rcBadOrdinal0 = uint32_t(m0) - 1 - interiorBadOrdinal;
        const MarkerId r0m_bad = fixture.assembler->getMarkerId(r0m, rcBadOrdinal0);

        const MarkerId r0_good = fixture.assembler->getMarkerId(r0p, interiorGoodOrdinal);
        const MarkerId r1_good = fixture.assembler->getMarkerId(r1p, interiorGoodOrdinal);
        const MarkerId r2_good = fixture.assembler->getMarkerId(r2p, interiorGoodOrdinal);
        const MarkerId r3_good = fixture.assembler->getMarkerId(r3p, interiorGoodOrdinal);

        // Build 8 vertices in order.
        mg.vertices().clear();

        // Vertex 0: start boundary (low coverage, 2 markers).
        mg.vertices().appendVector();
        {
            std::vector<MarkerId> ids = {r0_first, r1_first};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(0));
            }
        }

        // Vertex 1: rc partner (unused content for this test).
        mg.vertices().appendVector();

        // Vertex 2: end boundary (low coverage, 2 markers).
        mg.vertices().appendVector();
        {
            std::vector<MarkerId> ids = {r0_after, r1_after};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(2));
            }
        }

        // Vertex 3: rc partner (unused content).
        mg.vertices().appendVector();

        // Vertex 4: interior (highest coverage but invalid due to duplicate ReadId: r0+ and r0-).
        mg.vertices().appendVector();
        {
            std::vector<MarkerId> ids = {r0_bad, r0m_bad, r1_bad, r2_bad, r3_bad};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(4));
            }
        }

        // Vertex 5: rc partner (unused content).
        mg.vertices().appendVector();

        // Vertex 6: interior (next-best, valid, unique ReadIds).
        mg.vertices().appendVector();
        {
            std::vector<MarkerId> ids = {r0_good, r1_good, r2_good, r3_good};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(6));
            }
        }

        // Vertex 7: rc partner for vertex 6 (valid).
        mg.vertices().appendVector();
        {
            const uint32_t rcGoodOrdinal0 = uint32_t(m0) - 1 - interiorGoodOrdinal;
            const uint32_t rcGoodOrdinal1 = uint32_t(m1) - 1 - interiorGoodOrdinal;
            const uint32_t rcGoodOrdinal2 = uint32_t(m2) - 1 - interiorGoodOrdinal;
            const uint32_t rcGoodOrdinal3 = uint32_t(m3) - 1 - interiorGoodOrdinal;

            const MarkerId r0m_good = fixture.assembler->getMarkerId(r0m, rcGoodOrdinal0);
            const MarkerId r1m_good = fixture.assembler->getMarkerId(r1m, rcGoodOrdinal1);
            const MarkerId r2m_good = fixture.assembler->getMarkerId(r2m, rcGoodOrdinal2);
            const MarkerId r3m_good = fixture.assembler->getMarkerId(r3m, rcGoodOrdinal3);

            std::vector<MarkerId> ids = {r0m_good, r1m_good, r2m_good, r3m_good};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(7));
            }
        }
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() == 2); // one selected vertex => (forward + rc)

    const auto a0 = (*anchors)[0];
    REQUIRE(a0.coverage() == 4);
    // All plus-strand read ids appear at the interior ordinal.
    {
        std::vector<OrientedReadId> ids;
        for(const auto& mi : a0) {
            ids.push_back(mi.orientedReadId);
        }
        std::sort(ids.begin(), ids.end(), [](const OrientedReadId& a, const OrientedReadId& b) {
            return a.getValue() < b.getValue();
        });
        REQUIRE(ids.size() == 4);
        CHECK(ids[0] == OrientedReadId(ReadId(0), 0));
        CHECK(ids[1] == OrientedReadId(ReadId(1), 0));
        CHECK(ids[2] == OrientedReadId(ReadId(2), 0));
        CHECK(ids[3] == OrientedReadId(ReadId(3), 0));
    }
}

TEST_CASE("Integration: Best-per-interval anchors fall back to duplicate-ReadId vertex when no clean vertex exists", "[integration][anchors][events][best][fallback]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(4000, 701), randomSequence(4000, 702)});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);

    const uint32_t interiorOrdinal = 2;
    const uint32_t rcInterior0 = uint32_t(m0) - 1 - interiorOrdinal;

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        AlignmentInfo info;
        info.alignmentId = 0;
        info.isInReadGraph = 1;
        info.markerCount = lastOrdinal + 1 - firstOrdinal;
        info.data[0] = AlignmentInfo::Data(m0, firstOrdinal, lastOrdinal);
        info.data[1] = AlignmentInfo::Data(m1, firstOrdinal, lastOrdinal);
        info.minOrdinalOffset = 0;
        info.maxOrdinalOffset = 0;
        info.averageOrdinalOffset = 0;
        info.maxSkip = 0;
        info.maxDrift = 0;

        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);
        fixture.assembler->readGraph.edges.resize(2);
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0p, r1p};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[0] = e;
        }
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0m, r1m};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[1] = e;
        }
        auto& conn = fixture.assembler->readGraph.connectivity;
        conn.appendVector(); // r0+
        conn.append(0);
        conn.appendVector(); // r0-
        conn.append(1);
        conn.appendVector(); // r1+
        conn.append(0);
        conn.appendVector(); // r1-
        conn.append(1);

        // Marker graph with only one interior vertex in range, but it has duplicate ReadId (r0+ and r0-).
        auto& mg = fixture.assembler->markerGraph;
        mg.constructVertices();
        mg.vertices().createNew("", 4096);
        mg.vertexTable.createNew("", 4096);
        mg.vertexTable.resize(markers->totalSize());
        std::fill(mg.vertexTable.begin(), mg.vertexTable.end(), MarkerGraph::invalidCompressedVertexId);

        mg.reverseComplementVertex.createNew("", 4096);
        mg.reverseComplementVertex.resize(2);
        mg.reverseComplementVertex[0] = 1;
        mg.reverseComplementVertex[1] = 0;

        mg.vertices().clear();
        mg.vertices().appendVector();
        {
            const MarkerId r0p_mid = fixture.assembler->getMarkerId(r0p, interiorOrdinal);
            const MarkerId r0m_mid = fixture.assembler->getMarkerId(r0m, rcInterior0);
            const MarkerId r1p_mid = fixture.assembler->getMarkerId(r1p, interiorOrdinal);
            std::vector<MarkerId> ids = {r0p_mid, r0m_mid, r1p_mid};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(0));
            }
        }
        mg.vertices().appendVector(); // rc partner unused
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromMarkerGraphVerticesBestPerOverlapInterval(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() == 2); // selected vertex => (forward + rc)

    const auto a0 = (*anchors)[0];
    // Must contain both strands for read 0 (duplicate ReadId), because no clean alternative exists.
    bool hasR0p = false;
    bool hasR0m = false;
    for(const auto& mi : a0) {
        if(mi.orientedReadId == r0p) hasR0p = true;
        if(mi.orientedReadId == r0m) hasR0m = true;
    }
    CHECK(hasR0p);
    CHECK(hasR0m);
}

TEST_CASE("Integration: Best-per-interval anchor decomposition splits a vertex bridged by a single read", "[integration][anchors][events][best][decompose]") {
    AssemblerIntegrationFixture fixture;

    // Five reads. We construct a single marker graph vertex containing markers from all five,
    // but only one read (read 4) connects the two groups via overlaps. The decomposed anchor
    // method must split this into two anchors (plus their reverse complements).
    fixture.createFastq({
        randomSequence(4000, 801),
        randomSequence(4000, 802),
        randomSequence(4000, 803),
        randomSequence(4000, 804),
        randomSequence(4000, 805),
    });
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const OrientedReadId r2p(ReadId(2), 0);
    const OrientedReadId r2m(ReadId(2), 1);
    const OrientedReadId r3p(ReadId(3), 0);
    const OrientedReadId r3m(ReadId(3), 1);
    const OrientedReadId r4p(ReadId(4), 0);
    const OrientedReadId r4m(ReadId(4), 1);

    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    const uint32_t m2 = uint32_t(markers->size(r2p.getValue()));
    const uint32_t m3 = uint32_t(markers->size(r3p.getValue()));
    const uint32_t m4 = uint32_t(markers->size(r4p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);
    REQUIRE(m2 > 16);
    REQUIRE(m3 > 16);
    REQUIRE(m4 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    const uint32_t interiorOrdinal = 2;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);
    REQUIRE(lastOrdinal + 1 < m2);
    REQUIRE(lastOrdinal + 1 < m3);
    REQUIRE(lastOrdinal + 1 < m4);

    withSilencedIoInDir(fixture.testDir, [&] {
        // Alignments:
        // - Group A: (0,1)
        // - Group B: (2,3)
        // - Bridge read 4 overlaps to all in both groups.
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(6);

        auto fillAlignment = [&](size_t idx, ReadId a, ReadId b, uint32_t ma, uint32_t mb) {
            AlignmentInfo info;
            info.alignmentId = idx;
            info.isInReadGraph = 1;
            info.markerCount = lastOrdinal + 1 - firstOrdinal;
            info.data[0] = AlignmentInfo::Data(ma, firstOrdinal, lastOrdinal);
            info.data[1] = AlignmentInfo::Data(mb, firstOrdinal, lastOrdinal);
            info.minOrdinalOffset = 0;
            info.maxOrdinalOffset = 0;
            info.averageOrdinalOffset = 0;
            info.maxSkip = 0;
            info.maxDrift = 0;

            AlignmentData ad(OrientedReadPair(a, b, true), info);
            fixture.assembler->alignmentData[idx] = ad;
        };

        fillAlignment(0, ReadId(0), ReadId(1), m0, m1);
        fillAlignment(1, ReadId(2), ReadId(3), m2, m3);
        fillAlignment(2, ReadId(0), ReadId(4), m0, m4);
        fillAlignment(3, ReadId(1), ReadId(4), m1, m4);
        fillAlignment(4, ReadId(2), ReadId(4), m2, m4);
        fillAlignment(5, ReadId(3), ReadId(4), m3, m4);

        // ReadGraph edges: each alignment yields an edge for + and an edge for -.
        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);
        fixture.assembler->readGraph.edges.resize(12);

        auto setEdge = [&](size_t edgeIndex, OrientedReadId a, OrientedReadId b, uint64_t alignmentId) {
            ReadGraphEdge e;
            e.orientedReadIds = {a, b};
            e.alignmentId = alignmentId;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[edgeIndex] = e;
        };
        // Alignment 0: (0,1)
        setEdge(0, r0p, r1p, 0);
        setEdge(1, r0m, r1m, 0);
        // Alignment 1: (2,3)
        setEdge(2, r2p, r3p, 1);
        setEdge(3, r2m, r3m, 1);
        // Alignment 2: (0,4)
        setEdge(4, r0p, r4p, 2);
        setEdge(5, r0m, r4m, 2);
        // Alignment 3: (1,4)
        setEdge(6, r1p, r4p, 3);
        setEdge(7, r1m, r4m, 3);
        // Alignment 4: (2,4)
        setEdge(8, r2p, r4p, 4);
        setEdge(9, r2m, r4m, 4);
        // Alignment 5: (3,4)
        setEdge(10, r3p, r4p, 5);
        setEdge(11, r3m, r4m, 5);

        auto& conn = fixture.assembler->readGraph.connectivity;
        // orientedReadCount = 10 (5 reads * 2 strands). Build vectors in order 0..9.
        conn.appendVector(); // 0 r0+
        conn.append(0); conn.append(4);
        conn.appendVector(); // 1 r0-
        conn.append(1); conn.append(5);

        conn.appendVector(); // 2 r1+
        conn.append(0); conn.append(6);
        conn.appendVector(); // 3 r1-
        conn.append(1); conn.append(7);

        conn.appendVector(); // 4 r2+
        conn.append(2); conn.append(8);
        conn.appendVector(); // 5 r2-
        conn.append(3); conn.append(9);

        conn.appendVector(); // 6 r3+
        conn.append(2); conn.append(10);
        conn.appendVector(); // 7 r3-
        conn.append(3); conn.append(11);

        conn.appendVector(); // 8 r4+
        conn.append(4); conn.append(6); conn.append(8); conn.append(10);
        conn.appendVector(); // 9 r4-
        conn.append(5); conn.append(7); conn.append(9); conn.append(11);

        // MarkerGraph: one vertex that contains all five reads on the + strand at interiorOrdinal.
        auto& mg = fixture.assembler->markerGraph;
        mg.constructVertices();
        mg.vertices().createNew("", 4096);
        mg.vertexTable.createNew("", 4096);
        mg.vertexTable.resize(markers->totalSize());
        std::fill(mg.vertexTable.begin(), mg.vertexTable.end(), MarkerGraph::invalidCompressedVertexId);

        mg.reverseComplementVertex.createNew("", 4096);
        mg.reverseComplementVertex.resize(2);
        mg.reverseComplementVertex[0] = 1;
        mg.reverseComplementVertex[1] = 0;

        mg.vertices().clear();
        mg.vertices().appendVector();
        {
            const MarkerId r0_mid = fixture.assembler->getMarkerId(r0p, interiorOrdinal);
            const MarkerId r1_mid = fixture.assembler->getMarkerId(r1p, interiorOrdinal);
            const MarkerId r2_mid = fixture.assembler->getMarkerId(r2p, interiorOrdinal);
            const MarkerId r3_mid = fixture.assembler->getMarkerId(r3p, interiorOrdinal);
            const MarkerId r4_mid = fixture.assembler->getMarkerId(r4p, interiorOrdinal);
            std::vector<MarkerId> ids = {r0_mid, r1_mid, r2_mid, r3_mid, r4_mid};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(0));
            }
        }
        mg.vertices().appendVector(); // rc partner unused
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromMarkerGraphVerticesBestPerOverlapIntervalDecomposed(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() == 4); // two decomposed anchors => (forward + rc) each

    bool foundGroupA = false;
    bool foundGroupB = false;
    for(dinara::mode3::AnchorId anchorId=0; anchorId<anchors->size(); ++anchorId) {
        const auto a = (*anchors)[anchorId];
        bool hasR0 = false, hasR1 = false, hasR2 = false, hasR3 = false, hasR4 = false;
        for(const auto& mi : a) {
            if(mi.orientedReadId == r0p) hasR0 = true;
            if(mi.orientedReadId == r1p) hasR1 = true;
            if(mi.orientedReadId == r2p) hasR2 = true;
            if(mi.orientedReadId == r3p) hasR3 = true;
            if(mi.orientedReadId == r4p) hasR4 = true;
        }
        if(hasR0) {
            CHECK(hasR1);
            CHECK_FALSE(hasR2);
            CHECK_FALSE(hasR3);
            CHECK_FALSE(hasR4);
            foundGroupA = true;
        }
        if(hasR2) {
            CHECK(hasR3);
            CHECK_FALSE(hasR0);
            CHECK_FALSE(hasR1);
            CHECK_FALSE(hasR4);
            foundGroupB = true;
        }
    }
    CHECK(foundGroupA);
    CHECK(foundGroupB);
}

TEST_CASE("Integration: Best-per-interval anchor decomposition splits a vertex connected by multiple weak cross-edges", "[integration][anchors][events][best][decompose][multibridge]") {
    AssemblerIntegrationFixture fixture;

    // Six reads. Two tight clusters (0,1,2) and (3,4,5) with two weak cross-edges (0-3 and 1-4).
    // This graph has no articulation point, but cross-edges have divergent neighbor sets and should be dropped.
    fixture.createFastq({
        randomSequence(4000, 901),
        randomSequence(4000, 902),
        randomSequence(4000, 903),
        randomSequence(4000, 904),
        randomSequence(4000, 905),
        randomSequence(4000, 906),
    });
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const OrientedReadId r2p(ReadId(2), 0);
    const OrientedReadId r2m(ReadId(2), 1);
    const OrientedReadId r3p(ReadId(3), 0);
    const OrientedReadId r3m(ReadId(3), 1);
    const OrientedReadId r4p(ReadId(4), 0);
    const OrientedReadId r4m(ReadId(4), 1);
    const OrientedReadId r5p(ReadId(5), 0);
    const OrientedReadId r5m(ReadId(5), 1);

    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    const uint32_t m2 = uint32_t(markers->size(r2p.getValue()));
    const uint32_t m3 = uint32_t(markers->size(r3p.getValue()));
    const uint32_t m4 = uint32_t(markers->size(r4p.getValue()));
    const uint32_t m5 = uint32_t(markers->size(r5p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);
    REQUIRE(m2 > 16);
    REQUIRE(m3 > 16);
    REQUIRE(m4 > 16);
    REQUIRE(m5 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    const uint32_t interiorOrdinal = 2;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);
    REQUIRE(lastOrdinal + 1 < m2);
    REQUIRE(lastOrdinal + 1 < m3);
    REQUIRE(lastOrdinal + 1 < m4);
    REQUIRE(lastOrdinal + 1 < m5);

    withSilencedIoInDir(fixture.testDir, [&] {
        // Alignments for cluster A.
        // Alignments for cluster B.
        // Two cross edges: (0,3) and (1,4).
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(8);

        auto fillAlignment = [&](size_t idx, ReadId a, ReadId b, uint32_t ma, uint32_t mb) {
            AlignmentInfo info;
            info.alignmentId = idx;
            info.isInReadGraph = 1;
            info.markerCount = lastOrdinal + 1 - firstOrdinal;
            info.data[0] = AlignmentInfo::Data(ma, firstOrdinal, lastOrdinal);
            info.data[1] = AlignmentInfo::Data(mb, firstOrdinal, lastOrdinal);
            info.minOrdinalOffset = 0;
            info.maxOrdinalOffset = 0;
            info.averageOrdinalOffset = 0;
            info.maxSkip = 0;
            info.maxDrift = 0;

            AlignmentData ad(OrientedReadPair(a, b, true), info);
            fixture.assembler->alignmentData[idx] = ad;
        };

        // A clique: 0-1, 0-2, 1-2
        fillAlignment(0, ReadId(0), ReadId(1), m0, m1);
        fillAlignment(1, ReadId(0), ReadId(2), m0, m2);
        fillAlignment(2, ReadId(1), ReadId(2), m1, m2);
        // B clique: 3-4, 3-5, 4-5
        fillAlignment(3, ReadId(3), ReadId(4), m3, m4);
        fillAlignment(4, ReadId(3), ReadId(5), m3, m5);
        fillAlignment(5, ReadId(4), ReadId(5), m4, m5);
        // Cross edges
        fillAlignment(6, ReadId(0), ReadId(3), m0, m3);
        fillAlignment(7, ReadId(1), ReadId(4), m1, m4);

        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);
        fixture.assembler->readGraph.edges.resize(16);

        auto setEdge = [&](size_t edgeIndex, OrientedReadId a, OrientedReadId b, uint64_t alignmentId) {
            ReadGraphEdge e;
            e.orientedReadIds = {a, b};
            e.alignmentId = alignmentId;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[edgeIndex] = e;
        };

        // For each alignment, add + and -.
        setEdge(0, r0p, r1p, 0); setEdge(1, r0m, r1m, 0);
        setEdge(2, r0p, r2p, 1); setEdge(3, r0m, r2m, 1);
        setEdge(4, r1p, r2p, 2); setEdge(5, r1m, r2m, 2);
        setEdge(6, r3p, r4p, 3); setEdge(7, r3m, r4m, 3);
        setEdge(8, r3p, r5p, 4); setEdge(9, r3m, r5m, 4);
        setEdge(10, r4p, r5p, 5); setEdge(11, r4m, r5m, 5);
        setEdge(12, r0p, r3p, 6); setEdge(13, r0m, r3m, 6);
        setEdge(14, r1p, r4p, 7); setEdge(15, r1m, r4m, 7);

        auto& conn = fixture.assembler->readGraph.connectivity;
        // orientedReadCount = 12 (6 reads * 2 strands) in order 0..11.
        conn.appendVector(); // r0+
        conn.append(0); conn.append(2); conn.append(12);
        conn.appendVector(); // r0-
        conn.append(1); conn.append(3); conn.append(13);

        conn.appendVector(); // r1+
        conn.append(0); conn.append(4); conn.append(14);
        conn.appendVector(); // r1-
        conn.append(1); conn.append(5); conn.append(15);

        conn.appendVector(); // r2+
        conn.append(2); conn.append(4);
        conn.appendVector(); // r2-
        conn.append(3); conn.append(5);

        conn.appendVector(); // r3+
        conn.append(6); conn.append(8); conn.append(12);
        conn.appendVector(); // r3-
        conn.append(7); conn.append(9); conn.append(13);

        conn.appendVector(); // r4+
        conn.append(6); conn.append(10); conn.append(14);
        conn.appendVector(); // r4-
        conn.append(7); conn.append(11); conn.append(15);

        conn.appendVector(); // r5+
        conn.append(8); conn.append(10);
        conn.appendVector(); // r5-
        conn.append(9); conn.append(11);

        // MarkerGraph: one vertex contains all six reads on + strand at interiorOrdinal.
        auto& mg = fixture.assembler->markerGraph;
        mg.constructVertices();
        mg.vertices().createNew("", 4096);
        mg.vertexTable.createNew("", 4096);
        mg.vertexTable.resize(markers->totalSize());
        std::fill(mg.vertexTable.begin(), mg.vertexTable.end(), MarkerGraph::invalidCompressedVertexId);

        mg.reverseComplementVertex.createNew("", 4096);
        mg.reverseComplementVertex.resize(2);
        mg.reverseComplementVertex[0] = 1;
        mg.reverseComplementVertex[1] = 0;

        mg.vertices().clear();
        mg.vertices().appendVector();
        {
            const MarkerId r0_mid = fixture.assembler->getMarkerId(r0p, interiorOrdinal);
            const MarkerId r1_mid = fixture.assembler->getMarkerId(r1p, interiorOrdinal);
            const MarkerId r2_mid = fixture.assembler->getMarkerId(r2p, interiorOrdinal);
            const MarkerId r3_mid = fixture.assembler->getMarkerId(r3p, interiorOrdinal);
            const MarkerId r4_mid = fixture.assembler->getMarkerId(r4p, interiorOrdinal);
            const MarkerId r5_mid = fixture.assembler->getMarkerId(r5p, interiorOrdinal);
            std::vector<MarkerId> ids = {r0_mid, r1_mid, r2_mid, r3_mid, r4_mid, r5_mid};
            std::sort(ids.begin(), ids.end());
            for(const MarkerId id : ids) {
                mg.vertices().append(id);
                mg.vertexTable[id] = MarkerGraph::CompressedVertexId(uint64_t(0));
            }
        }
        mg.vertices().appendVector(); // rc partner unused
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromMarkerGraphVerticesBestPerOverlapIntervalDecomposed(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() == 4); // two anchors (A and B) + RC each

    bool foundA = false;
    bool foundB = false;
    for(dinara::mode3::AnchorId anchorId=0; anchorId<anchors->size(); ++anchorId) {
        const auto a = (*anchors)[anchorId];
        bool has0=false, has1=false, has2=false, has3=false, has4=false, has5=false;
        for(const auto& mi : a) {
            if(mi.orientedReadId == r0p) has0 = true;
            if(mi.orientedReadId == r1p) has1 = true;
            if(mi.orientedReadId == r2p) has2 = true;
            if(mi.orientedReadId == r3p) has3 = true;
            if(mi.orientedReadId == r4p) has4 = true;
            if(mi.orientedReadId == r5p) has5 = true;
        }
        if(has0 || has1 || has2) {
            CHECK(has0);
            CHECK(has1);
            CHECK(has2);
            CHECK_FALSE(has3);
            CHECK_FALSE(has4);
            CHECK_FALSE(has5);
            foundA = true;
        }
        if(has3 || has4 || has5) {
            CHECK_FALSE(has0);
            CHECK_FALSE(has1);
            CHECK_FALSE(has2);
            CHECK(has3);
            CHECK(has4);
            CHECK(has5);
            foundB = true;
        }
    }
    CHECK(foundA);
    CHECK(foundB);
}

TEST_CASE("Integration: Overlap-only best-per-interval anchors work without marker graph vertices", "[integration][anchors][overlaps][best]") {
    AssemblerIntegrationFixture fixture;

    // Two identical reads => markers (k-mers) match at the same ordinals.
    const std::string seq = randomSequence(4000, 1001);
    fixture.createFastq({seq, seq});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);

    auto* markers = fixture.assembler->markers.get();
    REQUIRE(markers != nullptr);

    const OrientedReadId r0p(ReadId(0), 0);
    const OrientedReadId r0m(ReadId(0), 1);
    const OrientedReadId r1p(ReadId(1), 0);
    const OrientedReadId r1m(ReadId(1), 1);
    const uint32_t m0 = uint32_t(markers->size(r0p.getValue()));
    const uint32_t m1 = uint32_t(markers->size(r1p.getValue()));
    REQUIRE(m0 > 16);
    REQUIRE(m1 > 16);

    const uint32_t firstOrdinal = 1;
    const uint32_t lastOrdinal = 3;
    REQUIRE(lastOrdinal + 1 < m0);
    REQUIRE(lastOrdinal + 1 < m1);

    withSilencedIoInDir(fixture.testDir, [&] {
        // Alignment + ReadGraph (same as the overlap-event test, but we do not create a marker graph).
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        AlignmentInfo info;
        info.alignmentId = 0;
        info.isInReadGraph = 1;
        info.markerCount = lastOrdinal + 1 - firstOrdinal;
        info.data[0] = AlignmentInfo::Data(m0, firstOrdinal, lastOrdinal);
        info.data[1] = AlignmentInfo::Data(m1, firstOrdinal, lastOrdinal);
        info.minOrdinalOffset = 0;
        info.maxOrdinalOffset = 0;
        info.averageOrdinalOffset = 0;
        info.maxSkip = 0;
        info.maxDrift = 0;

        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->readGraph.edges.createNew("", 4096);
        fixture.assembler->readGraph.connectivity.createNew("", 4096);
        fixture.assembler->readGraph.edges.resize(2);
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0p, r1p};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[0] = e;
        }
        {
            ReadGraphEdge e;
            e.orientedReadIds = {r0m, r1m};
            e.alignmentId = 0;
            e.crossesStrands = 0;
            e.hasInconsistentAlignment = 0;
            fixture.assembler->readGraph.edges[1] = e;
        }

        auto& conn = fixture.assembler->readGraph.connectivity;
        conn.appendVector(); // r0+
        conn.append(0);
        conn.appendVector(); // r0-
        conn.append(1);
        conn.appendVector(); // r1+
        conn.append(0);
        conn.appendVector(); // r1-
        conn.append(1);
    });

    auto anchors = withSilencedIoInDir(fixture.testDir, [&] {
        return fixture.assembler->createAnchorsFromOverlapsBestPerOverlapInterval(
            /*minAnchorCoverage*/ 2,
            /*maxAnchorCoverage*/ 100,
            /*threadCount*/ 2);
    });
    REQUIRE(anchors);
    REQUIRE(anchors->size() >= 2); // at least one anchor + its RC anchor

    // Verify we can find at least one anchor containing both reads on strand 0 and its RC on strand 1.
    bool foundForward = false;
    bool foundReverse = false;
    for(uint64_t i=0; i<anchors->size(); ++i) {
        const auto a = (*anchors)[i];
        if(a.coverage() != 2) {
            continue;
        }
        if(a[0].orientedReadId == r0p && a[1].orientedReadId == r1p) {
            foundForward = true;
        }
        if(a[0].orientedReadId == r0m && a[1].orientedReadId == r1m) {
            foundReverse = true;
        }
        if(foundForward && foundReverse) {
            break;
        }
    }
    CHECK(foundForward);
    CHECK(foundReverse);
}

TEST_CASE("Integration: ma_hit_flt keeps contained overlaps", "[integration][hifiasm][filter]") {
    AssemblerIntegrationFixture fixture;

    // 3 reads: 0 is contained in 1, and 2 has an internal overlap to 1.
    fixture.createFastq({randomSequence(1000, 1), randomSequence(1500, 2), randomSequence(1000, 3)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(2);

        // Alignment 0: QCONT (query read 0 is fully contained in target read 1).
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 1000;
            ad.ts = 250; ad.te = 1250;
            fixture.assembler->alignmentData[0] = ad;
        }

        // Alignment 1: internal match (should be removed by ma_hit_flt).
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(2), ReadId(1), true), info);
            ad.qs = 200; ad.qe = 800;   // internal on query
            ad.ts = 100; ad.te = 700;   // internal on target
            fixture.assembler->alignmentData[1] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();

        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->applyCoverageCuts(50, 1);
        fixture.assembler->filterHangingOverlaps(1000, 0.8, 50, 1);
    });

    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted());
    CHECK(fixture.assembler->alignmentData[1].isDeleted());
    CHECK((fixture.assembler->alignmentData[1].deleteReasons0 & AlignmentData::DeleteReasonHanging) != 0);
    CHECK((fixture.assembler->alignmentData[1].deleteReasons1 & AlignmentData::DeleteReasonHanging) != 0);
}

TEST_CASE("Integration: ma_hit_cut skips deleted endpoints", "[integration][hifiasm][filter]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 11), randomSequence(1000, 12)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        AlignmentInfo info;
        info.alignmentId = 0;
        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
        ad.qs = 0; ad.qe = 1000;
        ad.ts = 0; ad.te = 1000;
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->computeAlignmentTableForTesting();

        // With min_dp=2, a single overlap yields depth 1, so both reads get deleted in ma_hit_sub.
        fixture.assembler->filterLocalSegments(2, 1);

        // Hifiasm parity: ma_hit_cut skips overlaps if either endpoint read is deleted; it does not force-delete the edge.
        fixture.assembler->applyCoverageCuts(50, 1);
    });

    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted());
    CHECK(fixture.assembler->alignmentData[0].deleteReasons0 == AlignmentData::DeleteReasonNone);
    CHECK(fixture.assembler->alignmentData[0].deleteReasons1 == AlignmentData::DeleteReasonNone);
}

TEST_CASE("Integration: ma_hit_sub selects the max-depth interval", "[integration][hifiasm][filter][sub]") {
    AssemblerIntegrationFixture fixture;

    // Read 0 has depth>=2 only on [400,600); reads 1 and 2 have only depth 1 everywhere.
    fixture.createFastq({randomSequence(1000, 201), randomSequence(1000, 202), randomSequence(1000, 203)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(2);

        // Overlap 0: [0,600) on read 0.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 600;
            ad.ts = 0; ad.te = 600;
            fixture.assembler->alignmentData[0] = ad;
        }

        // Overlap 1: [400,1000) on read 0.
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(2), true), info);
            ad.qs = 400; ad.qe = 1000;
            ad.ts = 0; ad.te = 600;
            fixture.assembler->alignmentData[1] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(2, 1);
    });

    const auto v0 = fixture.assembler->getValidReadIntervalForTesting(ReadId(0));
    CHECK_FALSE(v0.isDeleted);
    CHECK(v0.start == 400);
    CHECK(v0.end == 600);

    const auto v1 = fixture.assembler->getValidReadIntervalForTesting(ReadId(1));
    const auto v2 = fixture.assembler->getValidReadIntervalForTesting(ReadId(2));
    CHECK(v1.isDeleted);
    CHECK(v2.isDeleted);
}

TEST_CASE("Integration: ma_hit_cut trims and normalizes overlaps to valid segments", "[integration][hifiasm][filter][cut]") {
    AssemblerIntegrationFixture fixture;

    // Build a case where read 0's valid region is [100,900] and read 1's valid region is [200,800].
    // The overlap (0,1) is full-length in raw coordinates and must be trimmed + normalized.
    fixture.createFastq({
        randomSequence(1000, 301),
        randomSequence(1000, 302),
        randomSequence(1000, 303),
        randomSequence(1000, 304),
    });
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        // Full-length overlap between reads 0 and 1.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 1000;
            ad.ts = 0; ad.te = 1000;
            fixture.assembler->alignmentData[0] = ad;
        }

        // Second overlap to give read 0 depth>=2 only on [100,900].
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(2), true), info);
            ad.qs = 100; ad.qe = 900;
            ad.ts = 0; ad.te = 800;
            fixture.assembler->alignmentData[1] = ad;
        }

        // Second overlap to give read 1 depth>=2 only on [200,800].
        {
            AlignmentInfo info;
            info.alignmentId = 2;
            AlignmentData ad(OrientedReadPair(ReadId(1), ReadId(3), true), info);
            ad.qs = 200; ad.qe = 800;
            ad.ts = 0; ad.te = 600;
            fixture.assembler->alignmentData[2] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(2, 1);

        const auto v0 = fixture.assembler->getValidReadIntervalForTesting(ReadId(0));
        const auto v1 = fixture.assembler->getValidReadIntervalForTesting(ReadId(1));
        REQUIRE_FALSE(v0.isDeleted);
        REQUIRE_FALSE(v1.isDeleted);
        REQUIRE(v0.start == 100);
        REQUIRE(v0.end == 900);
        REQUIRE(v1.start == 200);
        REQUIRE(v1.end == 800);

        fixture.assembler->applyCoverageCuts(50, 1);
    });

    const AlignmentData& a01 = fixture.assembler->alignmentData[0];
    REQUIRE_FALSE(a01.isDeleted());
    CHECK(a01.deleteReasons0 == AlignmentData::DeleteReasonNone);
    CHECK(a01.deleteReasons1 == AlignmentData::DeleteReasonNone);

    // After trimming and normalization:
    // - read0 segment is [100,900], read1 segment is [200,800]
    // - intersection in raw coordinates is [200,800]
    // - normalized: on read0 => [100,700], on read1 => [0,600]
    CHECK(a01.qs == 100);
    CHECK(a01.qe == 700);
    CHECK(a01.ts == 0);
    CHECK(a01.te == 600);
}

TEST_CASE("Integration: ma_hit_contained_advance does not RC-map swapped intervals", "[integration][hifiasm][filter]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 41), randomSequence(1000, 42)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        // One reverse-strand overlap. When viewed from read 1's perspective, the query/target
        // intervals must be swapped without RC-mapping (hifiasm set_reverse_overlap parity).
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), false), info);
            ad.qs = 0; ad.qe = 900;
            ad.ts = 0; ad.te = 900;
            fixture.assembler->alignmentData[0] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();

        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->applyCoverageCuts(50, 1);
        fixture.assembler->filterHangingOverlaps(1000, 0.8, 50, 1);
        fixture.assembler->removeContainedReads(1000, 0.8, 50, 1);
    });

    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted());
    CHECK((fixture.assembler->alignmentData[0].deleteReasons0 & AlignmentData::DeleteReasonContained) == 0);
    CHECK((fixture.assembler->alignmentData[0].deleteReasons1 & AlignmentData::DeleteReasonContained) == 0);
}

TEST_CASE("Integration: try_rescue_overlaps clears directional phase deletions on consensus span", "[integration][hifiasm][filter][rescue]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({
        randomSequence(1000, 501),
        randomSequence(1000, 502),
        randomSequence(1000, 503),
        randomSequence(1000, 504),
        randomSequence(1000, 505),
    });
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(4);

        // Four conflict overlaps: read 0 deletes (phase) but other endpoints keep.
        // All intervals share the consensus intersection [300,650).
        const std::array<std::pair<uint32_t, uint32_t>, 4> qIntervals = {{
            {0, 700},
            {100, 900},
            {200, 800},
            {300, 650},
        }};
        for (uint32_t i = 0; i < 4; ++i) {
            AlignmentInfo info;
            info.alignmentId = i;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1 + i), true), info);
            ad.qs = qIntervals[i].first;
            ad.qe = qIntervals[i].second;
            ad.ts = 0;
            ad.te = qIntervals[i].second - qIntervals[i].first;
            ad.addDeleteReasonsFromReadPerspective(ReadId(0), AlignmentData::DeleteReasonPhase);
            fixture.assembler->alignmentData[i] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->rescuePhasedOverlaps(4, 1);
    });

    for (uint32_t i = 0; i < 4; ++i) {
        const auto& ad = fixture.assembler->alignmentData[i];
        CHECK_FALSE(ad.isDeleted0());
        CHECK_FALSE(ad.isDeleted1());
        CHECK(ad.deleteReasons0 == AlignmentData::DeleteReasonNone);
        CHECK(ad.deleteReasons1 == AlignmentData::DeleteReasonNone);
    }
}

TEST_CASE("Integration: ma_hit_contained_advance compresses containment chains", "[integration][hifiasm][filter][contained][chain]") {
    AssemblerIntegrationFixture fixture;

    // Chain: 0 contained in 1, 1 contained in 2. Keep (2,3) as a dovetail so container stays alive.
    fixture.createFastq({
        randomSequence(500, 601),
        randomSequence(800, 602),
        randomSequence(1200, 603),
        randomSequence(1200, 604),
    });
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        // 0 in 1.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 500;
            ad.ts = 100; ad.te = 600;
            fixture.assembler->alignmentData[0] = ad;
        }
        // 1 in 2.
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(1), ReadId(2), true), info);
            ad.qs = 0; ad.qe = 800;
            ad.ts = 200; ad.te = 1000;
            fixture.assembler->alignmentData[1] = ad;
        }
        // 2 dovetail 3 (kept, and NOT a containment).
        {
            AlignmentInfo info;
            info.alignmentId = 2;
            AlignmentData ad(OrientedReadPair(ReadId(2), ReadId(3), true), info);
            ad.qs = 0; ad.qe = 1100;
            ad.ts = 100; ad.te = 1200;
            fixture.assembler->alignmentData[2] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->applyCoverageCuts(50, 1);
        fixture.assembler->filterHangingOverlaps(1000, 0.8, 50, 1);
        fixture.assembler->removeContainedReads(1000, 0.8, 50, 1);
    });

    CHECK(fixture.assembler->getContainmentRootForTesting(ReadId(0)) == ReadId(2));
    CHECK(fixture.assembler->getContainmentRootForTesting(ReadId(1)) == ReadId(2));
    CHECK(fixture.assembler->getContainmentRootForTesting(ReadId(2)) == ReadId(invalidReadId));

    // Contained overlaps are deleted with the contained reason.
    CHECK((fixture.assembler->alignmentData[0].deleteReasons0 & AlignmentData::DeleteReasonContained) != 0);
    CHECK((fixture.assembler->alignmentData[0].deleteReasons1 & AlignmentData::DeleteReasonContained) != 0);
    CHECK((fixture.assembler->alignmentData[1].deleteReasons0 & AlignmentData::DeleteReasonContained) != 0);
    CHECK((fixture.assembler->alignmentData[1].deleteReasons1 & AlignmentData::DeleteReasonContained) != 0);

    // The dovetail edge (2,3) remains.
    CHECK_FALSE(fixture.assembler->alignmentData[2].isDeleted());
}

TEST_CASE("Integration: StringGraph arcs follow hifiasm ma_hit2arc end selection", "[integration][stringgraph][arcs]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 701), randomSequence(1000, 702)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        // Construct a dovetail overlap in the query-to-target sense (qs > tl5, rev=0):
        // query interval [200,1000), target interval [0,800).
        AlignmentInfo info;
        info.alignmentId = 0;
        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
        ad.qs = 200; ad.qe = 1000;
        ad.ts = 0;   ad.te = 800;
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->applyCoverageCuts(50, 1);
        fixture.assembler->filterHangingOverlaps(1000, 0.8, 50, 1);

        std::vector<bool> keep(1, true);
        fixture.assembler->createStringGraphUsingSelectedAlignments(keep);
    });

    REQUIRE(fixture.assembler->stringGraph.arcs.size() == 2);

    const auto& a = fixture.assembler->stringGraph.arcs[0];
    const auto& b = fixture.assembler->stringGraph.arcs[1];

    // Expected: uEnd=0, vEnd=0 => 0+ -> 1+
    CHECK(a.from == OrientedReadId(ReadId(0), 0).getValue());
    CHECK(a.to == OrientedReadId(ReadId(1), 0).getValue());
    CHECK(a.len == 200);
    CHECK(a.overlapLen == 800);

    // Twin is (to^1)->(from^1).
    CHECK(b.from == (a.to ^ 1U));
    CHECK(b.to == (a.from ^ 1U));
    CHECK(b.len == a.len);
    CHECK(b.overlapLen == a.overlapLen);

    // Adjacency lists include both arcs.
    {
        std::vector<uint32_t> outFrom;
        for (const uint32_t arcId : fixture.assembler->stringGraph.outgoing[a.from]) {
            outFrom.push_back(arcId);
        }
        std::ranges::sort(outFrom);
        CHECK(outFrom == std::vector<uint32_t>({0}));

        std::vector<uint32_t> inTo;
        for (const uint32_t arcId : fixture.assembler->stringGraph.incoming[a.to]) {
            inTo.push_back(arcId);
        }
        std::ranges::sort(inTo);
        CHECK(inTo == std::vector<uint32_t>({0}));

        std::vector<uint32_t> outRcFrom;
        for (const uint32_t arcId : fixture.assembler->stringGraph.outgoing[b.from]) {
            outRcFrom.push_back(arcId);
        }
        std::ranges::sort(outRcFrom);
        CHECK(outRcFrom == std::vector<uint32_t>({1}));

        std::vector<uint32_t> inRcTo;
        for (const uint32_t arcId : fixture.assembler->stringGraph.incoming[b.to]) {
            inRcTo.push_back(arcId);
        }
        std::ranges::sort(inRcTo);
        CHECK(inRcTo == std::vector<uint32_t>({1}));
    }
}

TEST_CASE("Integration: StringGraph arcs use target forward coordinates with rev flag (no RC-mapping)", "[integration][stringgraph][arcs][rev]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 711), randomSequence(1000, 712)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        // Reverse-strand overlap where target interval is still in target's forward coordinates.
        // Pick values that trigger the qs>tl5 branch in hifiasm's ma_hit2arc for rev=1:
        // Choose values that survive ma_hit_flt internal-fraction checks and still satisfy qs>tl5:
        // tl5 = tl - te = 10, qs = 100 => qs>tl5 and ext5=min(qs,tl5)=10.
        AlignmentInfo info;
        info.alignmentId = 0;
        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), false), info);
        ad.qs = 100; ad.qe = 1000;
        ad.ts = 0;   ad.te = 990;
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->applyCoverageCuts(50, 1);
        fixture.assembler->filterHangingOverlaps(1000, 0.8, 50, 1);
    });

    std::vector<bool> keep(1, true);
    fixture.assembler->createStringGraphUsingSelectedAlignments(keep);

    CAPTURE(fixture.assembler->alignmentData.size());
    CAPTURE(fixture.assembler->alignmentData[0].qs);
    CAPTURE(fixture.assembler->alignmentData[0].qe);
    CAPTURE(fixture.assembler->alignmentData[0].ts);
    CAPTURE(fixture.assembler->alignmentData[0].te);
    CAPTURE(fixture.assembler->alignmentData[0].isSameStrand);
    const auto v0 = fixture.assembler->getValidReadIntervalForTesting(ReadId(0));
    const auto v1 = fixture.assembler->getValidReadIntervalForTesting(ReadId(1));
    CAPTURE(v0.start);
    CAPTURE(v0.end);
    CAPTURE(v0.isDeleted);
    CAPTURE(v1.start);
    CAPTURE(v1.end);
    CAPTURE(v1.isDeleted);

    REQUIRE(fixture.assembler->stringGraph.arcs.size() == 2);
    const auto& a = fixture.assembler->stringGraph.arcs[0];

    // Expected from hifiasm ma_hit2arc for rev=1 and qs>tl5: uEnd=0, vEnd=1 => 0+ -> 1-
    CHECK(a.from == OrientedReadId(ReadId(0), 0).getValue());
    CHECK(a.to == OrientedReadId(ReadId(1), 1).getValue());
    CHECK(a.len == 90);
    CHECK(a.overlapLen == 910);

    // Incoming/outgoing adjacency matches the arc endpoints.
    {
        std::vector<uint32_t> outFrom;
        for (const uint32_t arcId : fixture.assembler->stringGraph.outgoing[a.from]) {
            outFrom.push_back(arcId);
        }
        std::ranges::sort(outFrom);
        CHECK(outFrom == std::vector<uint32_t>({0}));

        std::vector<uint32_t> inTo;
        for (const uint32_t arcId : fixture.assembler->stringGraph.incoming[a.to]) {
            inTo.push_back(arcId);
        }
        std::ranges::sort(inTo);
        CHECK(inTo == std::vector<uint32_t>({0}));
    }
}

TEST_CASE("Integration: StringGraph transitive reduction removes transitive arcs (hifiasm asg_arc_del_trans)", "[integration][stringgraph][clean][transitive]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 721), randomSequence(1000, 722), randomSequence(1000, 723)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        // Create arcs:
        // 0+ -> 1+ len=100
        // 0+ -> 2+ len=200 (transitive via 1+)
        // 1+ -> 2+ len=50
        // Use qs branch (qs > tl5 with tl5=ts=0) so len = qs and overlapLen = ql - qs.
        auto makeAlignment = [&](uint64_t alignmentId, ReadId qn, ReadId tn, uint32_t len) {
            AlignmentInfo info;
            info.alignmentId = alignmentId;
            AlignmentData ad(OrientedReadPair(qn, tn, true), info);
            ad.qs = len;
            ad.qe = 1000;
            ad.ts = 0;
            ad.te = 1000 - len;
            fixture.assembler->alignmentData[alignmentId] = ad;
        };

        makeAlignment(0, ReadId(0), ReadId(1), 100);
        makeAlignment(1, ReadId(0), ReadId(2), 200);
        makeAlignment(2, ReadId(1), ReadId(2), 50);

        fixture.assembler->computeAlignmentTableForTesting();

        std::vector<bool> keep(3, true);
        fixture.assembler->createStringGraphUsingSelectedAlignments(keep);

        REQUIRE(fixture.assembler->stringGraph.arcs.size() == 6);

        // Run only transitive reduction (skip tip cutting by passing maxShortTipReads=0).
        fixture.assembler->cleanStringGraphInitialHifiasm(/*gapFuzz*/0, /*maxShortTipReads*/0);
    });

    // alignmentId=1 creates arcId=2 (and its twin arcId=3). It should be deleted as transitive.
    CHECK(fixture.assembler->stringGraph.arcs[2].del == 1);
    CHECK(fixture.assembler->stringGraph.arcs[3].del == 1);
    CHECK(fixture.assembler->stringGraph.arcs[0].del == 0);
    CHECK(fixture.assembler->stringGraph.arcs[4].del == 0);

    // After cleanup, outgoing adjacency from 0+ should only contain arcId=0 (0+->1+).
    const uint32_t v0 = OrientedReadId(ReadId(0), 0).getValue();
    std::vector<uint32_t> out0;
    for (const uint32_t arcId : fixture.assembler->stringGraph.outgoing[v0]) {
        out0.push_back(arcId);
    }
    CHECK(out0 == std::vector<uint32_t>({0}));
}

TEST_CASE("Integration: StringGraph cuts short tips by read count (hifiasm asg_cut_tip)", "[integration][stringgraph][clean][tip]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 731), randomSequence(1000, 732), randomSequence(1000, 733)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(2);

        // Chain: 0+ -> 1+ -> 2+.
        auto makeAlignment = [&](uint64_t alignmentId, ReadId qn, ReadId tn, uint32_t len) {
            AlignmentInfo info;
            info.alignmentId = alignmentId;
            AlignmentData ad(OrientedReadPair(qn, tn, true), info);
            ad.qs = len;
            ad.qe = 1000;
            ad.ts = 0;
            ad.te = 1000 - len;
            fixture.assembler->alignmentData[alignmentId] = ad;
        };
        makeAlignment(0, ReadId(0), ReadId(1), 100);
        makeAlignment(1, ReadId(1), ReadId(2), 100);

        fixture.assembler->computeAlignmentTableForTesting();

        std::vector<bool> keep(2, true);
        fixture.assembler->createStringGraphUsingSelectedAlignments(keep);

        REQUIRE(fixture.assembler->stringGraph.readDeleted.isOpen);
        REQUIRE(fixture.assembler->stringGraph.readDeleted.size() == 3);

        // With maxShortTipReads=2, the 3-read chain is NOT removed.
        fixture.assembler->cleanStringGraphInitialHifiasm(/*gapFuzz*/0, /*maxShortTipReads*/2);
        CHECK(fixture.assembler->stringGraph.readDeleted[ReadId(0)] == 0);
        CHECK(fixture.assembler->stringGraph.readDeleted[ReadId(1)] == 0);
        CHECK(fixture.assembler->stringGraph.readDeleted[ReadId(2)] == 0);

        // With maxShortTipReads=3, the 3-read tip unitig is removed.
        fixture.assembler->cleanStringGraphInitialHifiasm(/*gapFuzz*/0, /*maxShortTipReads*/3);
    });

    CHECK(fixture.assembler->stringGraph.readDeleted[ReadId(0)] == 1);
    CHECK(fixture.assembler->stringGraph.readDeleted[ReadId(1)] == 1);
    CHECK(fixture.assembler->stringGraph.readDeleted[ReadId(2)] == 1);

    // All arcs should be deleted.
    for (uint64_t arcId = 0; arcId < fixture.assembler->stringGraph.arcs.size(); ++arcId) {
        CHECK(fixture.assembler->stringGraph.arcs[arcId].del == 1);
    }
}

TEST_CASE("Integration: StringGraph pre-clean removes simple 1-step bubbles (topology only)", "[integration][stringgraph][clean][preclean]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({
        randomSequence(1000, 741), randomSequence(1000, 742), randomSequence(1000, 743),
        randomSequence(1000, 744), randomSequence(1000, 745), randomSequence(1000, 746),
        randomSequence(1000, 747)
    });
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(7);

        auto makeAlignment = [&](uint64_t alignmentId, ReadId qn, ReadId tn, uint32_t len) {
            AlignmentInfo info;
            info.alignmentId = alignmentId;
            AlignmentData ad(OrientedReadPair(qn, tn, true), info);
            ad.qs = len;
            ad.qe = 1000;
            ad.ts = 0;
            ad.te = 1000 - len;
            fixture.assembler->alignmentData[alignmentId] = ad;
        };

        // Long chain into v=3 so it isn't a tip (prevents tip cutting from deleting the bubble source):
        // 0 -> 1 -> 2 -> 3
        makeAlignment(0, ReadId(0), ReadId(1), 100);
        makeAlignment(1, ReadId(1), ReadId(2), 100);
        makeAlignment(2, ReadId(2), ReadId(3), 100);

        // Bubble at v=3: branches to 4 and 5 reconverge at 6.
        // Strong branch: len=100 (ol=900)
        // Weak branch:   len=400 (ol=600)
        makeAlignment(3, ReadId(3), ReadId(4), 100);  // 3->4
        makeAlignment(4, ReadId(3), ReadId(5), 400);  // 3->5 (weak)
        makeAlignment(5, ReadId(4), ReadId(6), 100);  // 4->6
        makeAlignment(6, ReadId(5), ReadId(6), 400);  // 5->6 (weak)

        fixture.assembler->computeAlignmentTableForTesting();
        std::vector<bool> keep(7, true);
        fixture.assembler->createStringGraphUsingSelectedAlignments(keep);

        // Run pre-clean (bubble removal + tip cutting).
        fixture.assembler->cleanStringGraphPreCleanHifiasm(/*maxShortTipReads*/3);
    });

    // alignmentId=4 (3->5) => arcId=8,9 should be deleted by bubble removal.
    CHECK(fixture.assembler->stringGraph.arcs[8].del == 1);
    CHECK(fixture.assembler->stringGraph.arcs[9].del == 1);

    // alignmentId=6 (5->6) => arcId=12,13 should be deleted by bubble removal.
    CHECK(fixture.assembler->stringGraph.arcs[12].del == 1);
    CHECK(fixture.assembler->stringGraph.arcs[13].del == 1);

    // Strong branch remains.
    CHECK(fixture.assembler->stringGraph.arcs[6].del == 0);   // 3->4
    CHECK(fixture.assembler->stringGraph.arcs[10].del == 0);  // 4->6
}

TEST_CASE("Integration: StringGraph drop-short-overlaps prunes weak overlaps by ratio", "[integration][stringgraph][clean][drop]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 751), randomSequence(1000, 752), randomSequence(1000, 753)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        // v=0 has two outgoing arcs:
        // 0->1 strong (len=100 => ol=900)
        // 0->2 weak   (len=700 => ol=300)
        // Also add 1->2 so the weak edge isn't needed.
        auto makeAlignment = [&](uint64_t alignmentId, ReadId qn, ReadId tn, uint32_t len) {
            AlignmentInfo info;
            info.alignmentId = alignmentId;
            AlignmentData ad(OrientedReadPair(qn, tn, true), info);
            ad.qs = len;
            ad.qe = 1000;
            ad.ts = 0;
            ad.te = 1000 - len;
            fixture.assembler->alignmentData[alignmentId] = ad;
        };
        makeAlignment(0, ReadId(0), ReadId(1), 100);
        makeAlignment(1, ReadId(0), ReadId(2), 700);
        makeAlignment(2, ReadId(1), ReadId(2), 100);

        fixture.assembler->computeAlignmentTableForTesting();
        std::vector<bool> keep(3, true);
        fixture.assembler->createStringGraphUsingSelectedAlignments(keep);

        // Drop with ratio 0.5: threshold ~ 450 => ol=300 should be deleted.
        fixture.assembler->cleanStringGraphDropShortOverlaps(0.5, /*minOverlapLen*/0);
    });

    // alignmentId=1 creates arcId=2 and twin arcId=3.
    CHECK(fixture.assembler->stringGraph.arcs[2].del == 1);
    CHECK(fixture.assembler->stringGraph.arcs[3].del == 1);
    CHECK(fixture.assembler->stringGraph.arcs[0].del == 0);
    CHECK(fixture.assembler->stringGraph.arcs[4].del == 0);

    // After cleanup, outgoing adjacency from 0+ should only contain arcId=0.
    const uint32_t v0 = OrientedReadId(ReadId(0), 0).getValue();
    std::vector<uint32_t> out0;
    for (const uint32_t arcId : fixture.assembler->stringGraph.outgoing[v0]) {
        out0.push_back(arcId);
    }
    CHECK(out0 == std::vector<uint32_t>({0}));
}

TEST_CASE("Integration: StringGraph breaks simple short cycles", "[integration][stringgraph][clean][cycle]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 771), randomSequence(1000, 772), randomSequence(1000, 773)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        auto makeAlignment = [&](uint64_t alignmentId, ReadId qn, ReadId tn, uint32_t len) {
            AlignmentInfo info;
            info.alignmentId = alignmentId;
            AlignmentData ad(OrientedReadPair(qn, tn, true), info);
            ad.qs = len;
            ad.qe = 1000;
            ad.ts = 0;
            ad.te = 1000 - len;
            fixture.assembler->alignmentData[alignmentId] = ad;
        };

        // 0 -> 1 -> 2 -> 0 cycle.
        makeAlignment(0, ReadId(0), ReadId(1), 100);
        makeAlignment(1, ReadId(1), ReadId(2), 100);
        makeAlignment(2, ReadId(2), ReadId(0), 100);

        fixture.assembler->computeAlignmentTableForTesting();
        std::vector<bool> keep(3, true);
        fixture.assembler->createStringGraphUsingSelectedAlignments(keep);

        // Break short cycles (<=100 reads).
        const uint64_t broken = fixture.assembler->cleanStringGraphBreakShortCycles(100);
        CHECK(broken > 0);
    });

    // The cycle should be broken by deleting exactly one overlap's arc pair on this simple case.
    // At least one of the three alignment-derived arc pairs must be deleted.
    const bool pair0Deleted = fixture.assembler->stringGraph.arcs[0].del && fixture.assembler->stringGraph.arcs[1].del;
    const bool pair1Deleted = fixture.assembler->stringGraph.arcs[2].del && fixture.assembler->stringGraph.arcs[3].del;
    const bool pair2Deleted = fixture.assembler->stringGraph.arcs[4].del && fixture.assembler->stringGraph.arcs[5].del;
    CHECK((pair0Deleted || pair1Deleted || pair2Deleted));
}

TEST_CASE("Integration: StringGraph explorer can list active vertices from arcs", "[integration][stringgraph][http]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 761), randomSequence(1000, 762)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(1);

        AlignmentInfo info;
        info.alignmentId = 0;
        AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
        ad.qs = 100; ad.qe = 1000;
        ad.ts = 0;   ad.te = 900;
        fixture.assembler->alignmentData[0] = ad;

        fixture.assembler->computeAlignmentTableForTesting();
        std::vector<bool> keep(1, true);
        fixture.assembler->createStringGraphUsingSelectedAlignments(keep);
    });

    // Request without readId, but with listActive=on should not error and should list at least one vertex.
    std::ostringstream html;
    fixture.assembler->exploreStringGraph(std::vector<std::string>{
        "exploreStringGraph", "listActive", "on", "activeLimit", "10"}, html);
    const std::string s = html.str();
    CHECK(s.find("Active vertices") != std::string::npos);
    CHECK(s.find("exploreStringGraph?readId=") != std::string::npos);
    const bool hasAtLeastOneVertex =
        (s.find("0-0") != std::string::npos) || (s.find("1-0") != std::string::npos);
    CHECK(hasAtLeastOneVertex);
}

TEST_CASE("Integration: UnitigGraph construction matches hifiasm ma_ug_gen boundaries", "[integration][unitiggraph][build]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({
        randomSequence(1000, 801),
        randomSequence(1000, 802),
        randomSequence(1000, 803),
        randomSequence(1000, 804),
        randomSequence(1000, 805),
    });
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(4);

        // Build a "fork into a merge" pattern at read 2:
        //   0 -> 1 -> 2
        //   3 -> 4 -> 2   (realized via (2,4) overlap whose RC twin is 4+ -> 2+)
        //
        // This makes read 2 have multiple incoming arcs, so unitigs stop before read 2.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 200; ad.qe = 1000;
            ad.ts = 0;   ad.te = 800;
            fixture.assembler->alignmentData[0] = ad;
        }
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(1), ReadId(2), true), info);
            ad.qs = 200; ad.qe = 1000;
            ad.ts = 0;   ad.te = 800;
            fixture.assembler->alignmentData[1] = ad;
        }
        {
            AlignmentInfo info;
            info.alignmentId = 2;
            AlignmentData ad(OrientedReadPair(ReadId(3), ReadId(4), true), info);
            ad.qs = 200; ad.qe = 1000;
            ad.ts = 0;   ad.te = 800;
            fixture.assembler->alignmentData[2] = ad;
        }
        {
            AlignmentInfo info;
            info.alignmentId = 3;
            // Overlap stored as (2,4). Choose coordinates that produce arc 2- -> 4-
            // so that the RC twin becomes 4+ -> 2+ (second incoming arc into read 2).
            AlignmentData ad(OrientedReadPair(ReadId(2), ReadId(4), true), info);
            ad.qs = 0;   ad.qe = 800;
            ad.ts = 200; ad.te = 1000;
            fixture.assembler->alignmentData[3] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->applyCoverageCuts(50, 1);
        fixture.assembler->filterHangingOverlaps(1000, 0.8, 50, 1);

        std::vector<bool> keep(4, true);
        fixture.assembler->createStringGraphUsingSelectedAlignments(keep);
        fixture.assembler->createUnitigGraphFromStringGraph();
    });

    // Expect 3 unitigs: [0,1], [3,4], and [2] (in some orientation).
    REQUIRE(fixture.assembler->unitigGraph.unitigs.size() == 3);

    const auto findUnitigContaining = [&](uint32_t stringVertex) -> uint32_t {
        for (uint32_t u = 0; u < fixture.assembler->unitigGraph.unitigs.size(); ++u) {
            const auto path = fixture.assembler->unitigGraph.unitigPaths[u];
            for (const uint64_t x : path) {
                if (uint32_t(x >> 32) == stringVertex) {
                    return u;
                }
            }
        }
        return uint32_t(invalidReadId);
    };

    const uint32_t v0p = OrientedReadId(ReadId(0), 0).getValue();
    const uint32_t v1p = OrientedReadId(ReadId(1), 0).getValue();
    const uint32_t v2m = OrientedReadId(ReadId(2), 1).getValue(); // unitig for read 2 starts at 2- (ma_ug_gen rule)
    const uint32_t v3p = OrientedReadId(ReadId(3), 0).getValue();
    const uint32_t v4p = OrientedReadId(ReadId(4), 0).getValue();

    const uint32_t u01 = findUnitigContaining(v0p);
    const uint32_t u34 = findUnitigContaining(v3p);
    const uint32_t u2 = findUnitigContaining(v2m);
    REQUIRE(u01 != uint32_t(invalidReadId));
    REQUIRE(u34 != uint32_t(invalidReadId));
    REQUIRE(u2 != uint32_t(invalidReadId));
    CHECK(u01 == findUnitigContaining(v1p));
    CHECK(u34 == findUnitigContaining(v4p));

    // There should be unitig-graph arcs corresponding to the boundary overlaps into read 2.
    // Specifically, one arc should be derived from the string arc 1+ -> 2+.
    uint64_t stringArc12 = std::numeric_limits<uint64_t>::max();
    for (uint64_t arcId = 0; arcId < fixture.assembler->stringGraph.arcs.size(); ++arcId) {
        const auto& arc = fixture.assembler->stringGraph.arcs[arcId];
        if (arc.del) continue;
        if (arc.from == v1p && arc.to == OrientedReadId(ReadId(2), 0).getValue()) {
            stringArc12 = arcId;
            break;
        }
    }
    REQUIRE(stringArc12 != std::numeric_limits<uint64_t>::max());

    bool foundUnitigArc12 = false;
    for (uint64_t arcId = 0; arcId < fixture.assembler->unitigGraph.arcs.size(); ++arcId) {
        const auto& ua = fixture.assembler->unitigGraph.arcs[arcId];
        if (ua.del) continue;
        if (ua.stringArcId != stringArc12) continue;
        const uint32_t fromUnitig = ua.from >> 1U;
        const uint32_t toUnitig = ua.to >> 1U;
        CHECK(fromUnitig == u01);
        CHECK(toUnitig == u2);
        foundUnitigArc12 = true;
        break;
    }
    CHECK(foundUnitigArc12);
}

TEST_CASE("Integration: UnitigGraph transitive reduction removes transitive arcs (hifiasm asg_arc_del_trans)", "[integration][unitiggraph][clean][transitive]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 811)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        // Construct a tiny unitig graph directly:
        // 0+ -> 1+ -> 2+ and 0+ -> 2+ (transitive).
        auto& ug = fixture.assembler->unitigGraph;
        ug.remove();
        ug.unitigs.createNew("", 4096);
        ug.unitigPaths.createNew("", 4096);
        ug.arcs.createNew("", 4096);
        ug.outgoing.createNew("", 4096);
        ug.incoming.createNew("", 4096);
        ug.unitigDeleted.createNew("", 4096);

        ug.unitigs.resize(3);
        ug.unitigDeleted.resize(3);
        std::fill(ug.unitigDeleted.begin(), ug.unitigDeleted.end(), uint8_t(0));

        // Minimal path placeholders.
        ug.unitigPaths.appendVector(std::vector<uint64_t>({}));
        ug.unitigPaths.appendVector(std::vector<uint64_t>({}));
        ug.unitigPaths.appendVector(std::vector<uint64_t>({}));

        auto addArcPair = [&](uint32_t from, uint32_t to, uint32_t len, uint32_t overlapLen) {
            UnitigGraphArc a;
            a.from = from;
            a.to = to;
            a.len = len;
            a.overlapLen = overlapLen;
            a.del = 0;
            ug.arcs.push_back(a);
            UnitigGraphArc b = a;
            b.from = to ^ 1U;
            b.to = from ^ 1U;
            ug.arcs.push_back(b);
        };

        const uint32_t v0p = (0U << 1) | 0U;
        const uint32_t v1p = (1U << 1) | 0U;
        const uint32_t v2p = (2U << 1) | 0U;
        addArcPair(v0p, v1p, 100, 900);
        addArcPair(v1p, v2p, 100, 900);
        addArcPair(v0p, v2p, 210, 790);

        const uint32_t vertexCount = 6;
        ug.outgoing.beginPass1(vertexCount);
        ug.incoming.beginPass1(vertexCount);
        for (uint64_t arcId = 0; arcId < ug.arcs.size(); ++arcId) {
            ug.outgoing.incrementCount(ug.arcs[arcId].from);
            ug.incoming.incrementCount(ug.arcs[arcId].to);
        }
        ug.outgoing.beginPass2();
        ug.incoming.beginPass2();
        for (uint64_t arcId = 0; arcId < ug.arcs.size(); ++arcId) {
            ug.outgoing.store(ug.arcs[arcId].from, uint32_t(arcId));
            ug.incoming.store(ug.arcs[arcId].to, uint32_t(arcId));
        }
        ug.outgoing.endPass2();
        ug.incoming.endPass2();
    });

    fixture.assembler->cleanUnitigGraphInitialHifiasm(/*gapFuzz*/0, /*maxShortTipUnitigs*/0);

    // The transitive arc 0+ -> 2+ is the third arc pair (arcId 4 and 5).
    CHECK(fixture.assembler->unitigGraph.arcs.size() == 6);
    CHECK(fixture.assembler->unitigGraph.arcs[4].from == ((0U << 1) | 0U));
    CHECK(fixture.assembler->unitigGraph.arcs[4].to == ((2U << 1) | 0U));
    CHECK(fixture.assembler->unitigGraph.arcs[4].del == 1);
    CHECK(fixture.assembler->unitigGraph.arcs[5].del == 1);

    // Outgoing from 0+ should now only include the arc to 1+.
    std::vector<uint32_t> out;
    for (const uint32_t arcId : fixture.assembler->unitigGraph.outgoing[(0U << 1) | 0U]) {
        out.push_back(arcId);
    }
    std::ranges::sort(out);
    CHECK(out == std::vector<uint32_t>({0}));
}

TEST_CASE("Integration: detect_chimeric_reads deletes edges for simple chimera", "[integration][hifiasm][filter]") {
    AssemblerIntegrationFixture fixture;

    fixture.createFastq({randomSequence(1000, 21), randomSequence(1000, 22), randomSequence(1000, 23)});
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(3);

        // For read 0 (length 1000), create a left anchor (qs==0) and right anchor (qe==rLen),
        // with a gap between them -> simple chimera.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 400;
            ad.ts = 0; ad.te = 400;
            fixture.assembler->alignmentData[0] = ad;
        }
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(2), true), info);
            ad.qs = 600; ad.qe = 1000;
            ad.ts = 0; ad.te = 400;
            fixture.assembler->alignmentData[1] = ad;
        }
        // Control overlap not involving read 0.
        {
            AlignmentInfo info;
            info.alignmentId = 2;
            AlignmentData ad(OrientedReadPair(ReadId(1), ReadId(2), true), info);
            ad.qs = 0; ad.qe = 1000;
            ad.ts = 0; ad.te = 1000;
            fixture.assembler->alignmentData[2] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();
        fixture.assembler->filterLocalSegments(0, 1);
        fixture.assembler->detectChimericReads(1);
    });

    CHECK(fixture.assembler->alignmentData[0].isDeleted());
    CHECK(fixture.assembler->alignmentData[1].isDeleted());
    CHECK_FALSE(fixture.assembler->alignmentData[2].isDeleted());
    CHECK((fixture.assembler->alignmentData[0].deleteReasons0 & AlignmentData::DeleteReasonChimeric) != 0);
    CHECK((fixture.assembler->alignmentData[0].deleteReasons1 & AlignmentData::DeleteReasonChimeric) != 0);
    CHECK((fixture.assembler->alignmentData[1].deleteReasons0 & AlignmentData::DeleteReasonChimeric) != 0);
    CHECK((fixture.assembler->alignmentData[1].deleteReasons1 & AlignmentData::DeleteReasonChimeric) != 0);
}

TEST_CASE("Integration: ONT chemical arc mask deletes overlaps for low-depth reads", "[integration][hifiasm][filter][ont][chemical]") {
    AssemblerIntegrationFixture fixture;

    // 4 reads. Read 0 has two full-length overlaps (min depth >=2) => not flagged.
    // Read 3 has two non-overlapping intervals => min depth == 0 => flagged.
    fixture.createFastq({
        randomSequence(2000, 31),
        randomSequence(2000, 32),
        randomSequence(2000, 33),
        randomSequence(2000, 34),
    });
    fixture.initAssembler();
    fixture.loadReads();

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->alignmentData.createNew("", 4096);
        fixture.assembler->alignmentData.resize(4);

        // Read 0 overlaps (0,1) and (0,2) spanning the full read.
        {
            AlignmentInfo info;
            info.alignmentId = 0;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(1), true), info);
            ad.qs = 0; ad.qe = 2000;
            ad.ts = 0; ad.te = 2000;
            fixture.assembler->alignmentData[0] = ad;
        }
        {
            AlignmentInfo info;
            info.alignmentId = 1;
            AlignmentData ad(OrientedReadPair(ReadId(0), ReadId(2), true), info);
            ad.qs = 0; ad.qe = 2000;
            ad.ts = 0; ad.te = 2000;
            fixture.assembler->alignmentData[1] = ad;
        }

        // Read 3 has two disjoint overlaps (gap in the middle).
        {
            AlignmentInfo info;
            info.alignmentId = 2;
            AlignmentData ad(OrientedReadPair(ReadId(1), ReadId(3), true), info);
            // Interval on read 3: [0,800)
            ad.qs = 0; ad.qe = 2000;
            ad.ts = 0; ad.te = 800;
            fixture.assembler->alignmentData[2] = ad;
        }
        {
            AlignmentInfo info;
            info.alignmentId = 3;
            AlignmentData ad(OrientedReadPair(ReadId(2), ReadId(3), true), info);
            // Interval on read 3: [1200,2000)
            ad.qs = 0; ad.qe = 2000;
            ad.ts = 1200; ad.te = 2000;
            fixture.assembler->alignmentData[3] = ad;
        }

        fixture.assembler->computeAlignmentTableForTesting();

        // Use a small flank and cov=1 to make the test robust on short reads.
        fixture.assembler->applyOntChemicalArcMask(1, 0, 0.02, 1);
    });

    // Overlaps incident to read 3 are chemically deleted.
    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted());
    CHECK_FALSE(fixture.assembler->alignmentData[1].isDeleted());
    CHECK(fixture.assembler->alignmentData[2].isDeleted());
    CHECK(fixture.assembler->alignmentData[3].isDeleted());
    CHECK((fixture.assembler->alignmentData[2].deleteReasons0 & AlignmentData::DeleteReasonChemical) != 0);
    CHECK((fixture.assembler->alignmentData[2].deleteReasons1 & AlignmentData::DeleteReasonChemical) != 0);
    CHECK((fixture.assembler->alignmentData[3].deleteReasons0 & AlignmentData::DeleteReasonChemical) != 0);
    CHECK((fixture.assembler->alignmentData[3].deleteReasons1 & AlignmentData::DeleteReasonChemical) != 0);
}

TEST_CASE("Integration: reverse overlap coordinates do not underflow", "[integration][hifiasm][coords]") {
    AssemblerIntegrationFixture fixture;

    const std::string seq = randomSequence(5000, 4242);
    const std::string seqRc = reverseComplement(seq);

    fixture.createFastq({seq, seqRc});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    const AlignmentData* ad = findAlignmentDataPtr(fixture.assembler->alignmentData, ReadId(0), ReadId(1), false);
    REQUIRE(ad != nullptr);
    REQUIRE_FALSE(ad->isSameStrand);

    const uint32_t len0 = uint32_t(fixture.assembler->getReads().getReadRawSequenceLength(ReadId(0)));
    const uint32_t len1 = uint32_t(fixture.assembler->getReads().getReadRawSequenceLength(ReadId(1)));

    CHECK(ad->qs <= ad->qe);
    CHECK(ad->qe <= len0);
    CHECK(ad->ts <= ad->te);
    CHECK(ad->te <= len1);
}

TEST_CASE("Large Matching Regions and SNP Chaining", "[evidence][delta]") {
    AssemblerIntegrationFixture fixture;
    
    // 1. Create reads of 40kb
    std::string base = randomSequence(40000, 123);
    std::string mutant = base;
    // Single SNP at 35,000
    mutant[35000] = (base[35000] == 'A') ? 'C' : 'A';
    
    fixture.createFastq({base, mutant});
    fixture.initAssembler();
    fixture.loadReads();

    // 2. Marker generation
    fixture.generateMarkers(32, 10);
    fixture.countKmers();
    fixture.applyFilter(0, 100);

    // 3. Build inverted index (Phases 1-4)
    fixture.buildIndex();

    // 4. Run DP chaining (Phase 5)
    fixture.chainCandidates(0.1, 100);

    // 5. Compute Projected Alignments on the chained candidates
    fixture.computeAlignments();

    const auto& store = fixture.assembler->alignedEvidenceStore;
    const auto& alignmentData = fixture.assembler->alignmentData;
    REQUIRE(store.index.size() > 0);
    REQUIRE(alignmentData.size() > 0);

    const uint32_t evidenceId = findAlignmentEvidenceId(alignmentData, ReadId(0), ReadId(1), true);
    REQUIRE(evidenceId != invalid<uint32_t>);

    // A single SNP at 35k should decode correctly (with hop tokens as needed).
    std::vector<uint32_t> positions;
    store.forEachSnp0InRange(evidenceId, 34990, 35010, [&](uint32_t pos, uint8_t /*base*/) {
        positions.push_back(pos);
    });
    REQUIRE(positions.size() == 1);
    CHECK(positions[0] == 35000);

    const auto tokens = store.getSnps0(evidenceId);
    const uint32_t expectedHops = 35000 / SnpEvidence::MAX_DELTA;
    REQUIRE(tokens.size() == expectedHops + 1);
    uint32_t hopCount = 0;
    for (const auto& t : tokens) if (t.isHop()) ++hopCount;
    CHECK(hopCount == expectedHops);
}

TEST_CASE("AlignedEvidenceStore delta encoding and range decoding", "[evidence][delta][checkpoints]") {
    using dinara::AlignedEvidenceStore;
    using dinara::SnpEvidence;

    // This is a focused correctness test for the delta-encoded SNP token stream:
    // - Handles deltas larger than 14 bits using hop tokens.
    // - Handles delta exactly MAX_DELTA using a hop + delta(0) SNP token.
    // - Checkpointed range decoding returns correct absolute positions.

    AlignedEvidenceStore store;
    store.reserve(4, 0, 0);

    // Alignment 0: exact MAX_DELTA.
    store.beginAlignment();
    store.addSnp0(SnpEvidence::MAX_DELTA, 2);
    store.addSnp1(SnpEvidence::MAX_DELTA, 3);
    {
        auto tokens = store.getSnps0(0);
        REQUIRE(tokens.size() == 2);
        CHECK(tokens[0].delta() == SnpEvidence::MAX_DELTA);
        CHECK(tokens[0].isHop());
        CHECK(tokens[1].delta() == 0);
        CHECK_FALSE(tokens[1].isHop());

        std::vector<uint32_t> positions;
        store.forEachSnp0InRange(0, 0, SnpEvidence::MAX_DELTA + 1, [&](uint32_t pos, uint8_t /*base*/) {
            positions.push_back(pos);
        });
        REQUIRE(positions.size() == 1);
        CHECK(positions[0] == SnpEvidence::MAX_DELTA);

        std::vector<uint32_t> positions1;
        store.forEachSnp1InRange(0, 0, SnpEvidence::MAX_DELTA + 1, [&](uint32_t pos, uint8_t /*base*/) {
            positions1.push_back(pos);
        });
        REQUIRE(positions1.size() == 1);
        CHECK(positions1[0] == SnpEvidence::MAX_DELTA);
    }

    // Alignment 1: huge delta requiring multiple hops.
    store.beginAlignment();
    const uint32_t hugePos = 50000;
    store.addSnp0(hugePos, 1);
    store.addSnp1(hugePos, 2);
    {
        auto tokens = store.getSnps0(1);
        // Expect floor(hugePos / MAX_DELTA) hops + 1 SNP token.
        const uint32_t expectedHops = hugePos / SnpEvidence::MAX_DELTA;
        REQUIRE(tokens.size() == expectedHops + 1);
        for (uint32_t i = 0; i < expectedHops; ++i) {
            CHECK(tokens[i].isHop());
            CHECK(tokens[i].delta() == SnpEvidence::MAX_DELTA);
        }
        CHECK_FALSE(tokens[expectedHops].isHop());

        std::vector<uint32_t> positions;
        store.forEachSnp0InRange(1, 0, hugePos + 1, [&](uint32_t pos, uint8_t /*base*/) {
            positions.push_back(pos);
        });
        REQUIRE(positions.size() == 1);
        CHECK(positions[0] == hugePos);

        // Range query that should include it.
        positions.clear();
        store.forEachSnp0InRange(1, hugePos, hugePos + 1, [&](uint32_t pos, uint8_t /*base*/) {
            positions.push_back(pos);
        });
        REQUIRE(positions.size() == 1);
        CHECK(positions[0] == hugePos);

        // Range query that should exclude it.
        positions.clear();
        store.forEachSnp0InRange(1, hugePos + 1, hugePos + 100, [&](uint32_t /*pos*/, uint8_t /*base*/) {
            positions.push_back(0);
        });
        CHECK(positions.empty());

        std::vector<uint32_t> positions1;
        store.forEachSnp1InRange(1, hugePos, hugePos + 1, [&](uint32_t pos, uint8_t /*base*/) {
            positions1.push_back(pos);
        });
        REQUIRE(positions1.size() == 1);
        CHECK(positions1[0] == hugePos);
    }

    // Alignment 2: dense tokens to exercise checkpoints and seeking.
    store.beginAlignment();
    for (uint32_t p = 0; p < 200; ++p) {
        store.addSnp0(p, uint8_t(p & 3));
    }
    {
        std::vector<uint32_t> positions;
        store.forEachSnp0InRange(2, 150, 155, [&](uint32_t pos, uint8_t /*base*/) {
            positions.push_back(pos);
        });
        REQUIRE(positions.size() == 5);
        CHECK(positions[0] == 150);
        CHECK(positions[1] == 151);
        CHECK(positions[2] == 152);
        CHECK(positions[3] == 153);
        CHECK(positions[4] == 154);
    }
}

TEST_CASE("AlignedEvidenceStore indel evidence handles small and large positions", "[evidence][indel]") {
    using dinara::AlignedEvidenceStore;

    AlignedEvidenceStore store;
    store.reserve(1, 0, 0);
    store.beginAlignment();

    // Store a small and a large indel position (absolute positions, monotonic).
    store.addIndel0(100, 25, 1);     // deletion
    store.addIndel0(50000, 30, 0);   // insertion

    auto indels = store.getIndels0(0);
    REQUIRE(indels.size() == 2);
    CHECK(indels[0].pos() == 100);
    CHECK(indels[0].len() == 25);
    CHECK(indels[0].isDeletion());

    CHECK(indels[1].pos() == 50000);
    CHECK(indels[1].len() == 30);
    CHECK(indels[1].isInsertion());

    // Mirror the range-scan logic used in detectSVSites (lower_bound + early stop).
    auto it = std::lower_bound(indels.begin(), indels.end(), uint32_t(49900),
        [](const dinara::IndelEvidence& e, uint32_t value) { return e.pos() < value; });
    REQUIRE(it != indels.end());
    CHECK(it->pos() == 50000);
}

// =============================================================================
// TEST: PAF Import and Chaining
// =============================================================================
TEST_CASE("Integration: PAF import and chaining", "[integration][paf][chaining]") {
    AssemblerIntegrationFixture fixture;
    
    // Create overlapping reads with shared sequence regions.
    // Read 0 and Read 1 share a large overlap region.
    std::string sharedRegion = randomSequence(800, 444);
    std::string prefix0 = randomSequence(200, 555);
    std::string suffix1 = randomSequence(200, 666);
    
    std::string seq0 = prefix0 + sharedRegion;  // read_0: [prefix0][shared]
    std::string seq1 = sharedRegion + suffix1;  // read_1: [shared][suffix1]
    
    // Read 2 is unrelated (should not be in PAF)
    std::string seq2 = randomSequence(1000, 777);
    
    fixture.createFastq({seq0, seq1, seq2});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);

    // Create PAF file with the known exact overlap between read_0 and read_1:
    // read_0: [prefix0][shared] => overlap is [200,1000)
    // read_1: [shared][suffix1] => overlap is [0,800)
    const std::string pafPath = fixture.createPafFile({
        AssemblerIntegrationFixture::PafOverlapSpec{
            "read_0", "read_1", true,
            uint32_t(prefix0.size()), uint32_t(seq0.size()),
            0U, uint32_t(sharedRegion.size())
        }
    });

    fixture.buildIndex();
    fixture.importPafCandidates(pafPath);
    const size_t importedCount = fixture.assembler->alignmentCandidates.candidates.size();
    CAPTURE(importedCount);
    REQUIRE(importedCount >= 1);

    bool foundPairImported = false;
    for (size_t i = 0; i < importedCount; ++i) {
        const auto& c = fixture.assembler->alignmentCandidates.candidates[i];
        if ((c.readIds[0] == 0 && c.readIds[1] == 1) || (c.readIds[0] == 1 && c.readIds[1] == 0)) {
            foundPairImported = true;
            CHECK(c.isSameStrand == true);
        }
    }
    REQUIRE(foundPairImported);

    fixture.chainPafCandidates(0.1, 100);
    const size_t chainedCount = fixture.assembler->alignmentCandidates.candidates.size();
    const size_t precomputedAlignments = fixture.assembler->alignmentCandidatesAlignmentsData.alignments.size();
    CAPTURE(chainedCount);
    CAPTURE(precomputedAlignments);
    REQUIRE(chainedCount > 0);
    CHECK(precomputedAlignments == chainedCount);

    // Compare that inverted-index discovery also finds the same pair.
    AssemblerIntegrationFixture fixture2;
    fixture2.createFastq({seq0, seq1, seq2});
    fixture2.initAssembler();
    fixture2.loadReads();
    fixture2.generateMarkers(16, 5);
    fixture2.countKmers();
    fixture2.applyFilter(1, 1000);
    fixture2.findCandidates();
    bool foundPairInverted = false;
    for (const auto& c : fixture2.assembler->alignmentCandidates.candidates) {
        if ((c.readIds[0] == 0 && c.readIds[1] == 1) || (c.readIds[0] == 1 && c.readIds[1] == 0)) {
            foundPairInverted = true;
            break;
        }
    }
    CHECK(foundPairInverted);

    // Chained PAF candidates can be used for alignment computation.
    fixture.computeAlignments();
    const size_t alignmentCount = fixture.assembler->alignmentData.size();
    CAPTURE(alignmentCount);
    REQUIRE(alignmentCount >= 1);

    const AlignmentData* ad01 = findAlignmentDataPtr(fixture.assembler->alignmentData, ReadId(0), ReadId(1), true);
    REQUIRE(ad01 != nullptr);
    const uint32_t evidenceId = uint32_t(ad01->info.alignmentId);

    // Alignment should cover most of the true overlap.
    CHECK(uint32_t(ad01->qe - ad01->qs) >= 700);
    CHECK(uint32_t(ad01->te - ad01->ts) >= 700);

    // This overlap has no injected variants, so evidence streams should be empty (or at least empty in the overlapped range).
    std::vector<uint32_t> snpPos;
    fixture.assembler->alignedEvidenceStore.forEachSnp0InRange(
        evidenceId, uint32_t(prefix0.size()), uint32_t(seq0.size()), [&](uint32_t pos, uint8_t) { snpPos.push_back(pos); });
    CHECK(snpPos.empty());

    const auto indels0 = fixture.assembler->alignedEvidenceStore.getIndels0(evidenceId);
    const auto indels1 = fixture.assembler->alignedEvidenceStore.getIndels1(evidenceId);
    const uint32_t ovBegin = uint32_t(prefix0.size());
    const uint32_t ovEnd = uint32_t(seq0.size());
    for (const auto& ev : indels0) {
        const bool outside = (ev.pos() < ovBegin) || (ev.pos() >= ovEnd);
        CHECK(outside);
    }
    for (const auto& ev : indels1) {
        const bool outside = (ev.pos() < ovBegin) || (ev.pos() >= ovEnd);
        CHECK(outside);
    }
}

TEST_CASE("Integration: AlignedEvidenceStore stores SNP and indel evidence", "[integration][evidence][variants]") {
    AssemblerIntegrationFixture fixture;

    // Create reads with deterministic variants.
    std::string prefix = randomSequence(500, 123);
    std::string mid = randomSequence(200, 789);
    std::string suffix = randomSequence(500, 456);

    std::string seq0 = prefix + mid + suffix;
    std::string seq1 = seq0;

    // Inject 5 SNPs in seq1.
    const std::array<uint32_t, 5> snpSites{510U, 540U, 570U, 600U, 630U};
    for (uint32_t site : snpSites) {
        seq1[site] = otherBase(seq1[site]);
    }

    // Inject a 15bp deletion in seq1 near the middle (affects indel evidence).
    seq1.erase(650, 15);

    fixture.createFastq({seq0, seq1});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(10, 3);
    fixture.countKmers();
    fixture.applyFilter(1, 100);
    fixture.findCandidates();
    fixture.computeAlignments();

    const auto& store = fixture.assembler->alignedEvidenceStore;
    const auto& alignmentData = fixture.assembler->alignmentData;
    REQUIRE(alignmentData.size() >= 1);

    const uint32_t evidenceId = findAlignmentEvidenceId(alignmentData, ReadId(0), ReadId(1), true);
    REQUIRE(evidenceId != invalid<uint32_t>);

    // Verify each injected SNP is present and projects a base from read_0 in target coordinates.
    const auto& reads = fixture.assembler->getReads();
    for (uint32_t site : snpSites) {
        std::vector<uint8_t> bases;
        store.forEachSnp0InRange(evidenceId, site, site + 1, [&](uint32_t pos, uint8_t base) {
            CHECK(pos == site);
            bases.push_back(base);
        });
        REQUIRE(bases.size() == 1);
        CHECK(bases[0] == reads.getOrientedReadBase(OrientedReadId(ReadId(0), 0), site).value);
    }

    // Verify the 15bp event exists in indel evidence near the expected locus.
    const auto indels0 = store.getIndels0(evidenceId);
    const auto indels1 = store.getIndels1(evidenceId);
    REQUIRE_FALSE((indels0.empty() && indels1.empty()));
    uint32_t totalIndelLen = 0;
    for (const auto& indel : indels0) if (indel.pos() >= 635 && indel.pos() <= 665) totalIndelLen += indel.len();
    for (const auto& indel : indels1) if (indel.pos() >= 635 && indel.pos() <= 665) totalIndelLen += indel.len();
    CHECK(totalIndelLen >= 15);
}

TEST_CASE("Integration: AlignedEvidenceStore handles large SNP deltas and mixed variants", "[integration][evidence][delta][indel]") {
    AssemblerIntegrationFixture fixture;

    // Create 20kb reads to test large SNP deltas (> 16383) with an indel.
    std::string seq0 = randomSequence(20000, 111);
    std::string seq1 = seq0;

    // Inject SNP at 100.
    seq1[100] = otherBase(seq1[100]);

    // Inject SNP at 18000 (far from the first SNP).
    seq1[18000] = otherBase(seq1[18000]);

    // Inject a 10bp deletion at 10000 (on target).
    seq1.erase(10000, 10);

    fixture.createFastq({seq0, seq1});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(10, 3);
    fixture.countKmers();
    fixture.applyFilter(1, 100);
    fixture.findCandidates();
    fixture.computeAlignments();

    const auto& store = fixture.assembler->alignedEvidenceStore;
    const auto& alignmentData = fixture.assembler->alignmentData;
    REQUIRE(alignmentData.size() >= 1);

    const uint32_t evidenceId = findAlignmentEvidenceId(alignmentData, ReadId(0), ReadId(1), true);
    REQUIRE(evidenceId != invalid<uint32_t>);

    // The far SNP should force hop tokens in the delta stream.
    const auto tokens = store.getSnps0(evidenceId);
    bool foundHop = false;
    for (const auto& t : tokens) if (t.isHop()) { foundHop = true; break; }
    CHECK(foundHop);

    // Verify we can find the near SNP at 100, and a far SNP around 18000 (may shift by indel in target coords).
    std::vector<uint32_t> positions;
    store.forEachSnp0InRange(evidenceId, 0, 20050, [&](uint32_t pos, uint8_t /*base*/) { positions.push_back(pos); });
    CHECK(std::find(positions.begin(), positions.end(), 100U) != positions.end());
    const bool hasFar =
        (std::find(positions.begin(), positions.end(), 17990U) != positions.end()) ||
        (std::find(positions.begin(), positions.end(), 18000U) != positions.end());
    CHECK(hasFar);

    // Verify the 10bp deletion exists near 10000.
    const auto indels0 = store.getIndels0(evidenceId);
    const auto indels1 = store.getIndels1(evidenceId);
    uint32_t totalMixedIndelLen = 0;
    for (const auto& indel : indels0) if (indel.pos() >= 9990 && indel.pos() <= 10010) totalMixedIndelLen += indel.len();
    for (const auto& indel : indels1) if (indel.pos() >= 9990 && indel.pos() <= 10010) totalMixedIndelLen += indel.len();
    CHECK(totalMixedIndelLen >= 10);
}

TEST_CASE("Integration: filterSecondaryAlignmentsPerReadPair removes redundant duplicates", "[integration][filtering][besthit]") {
    AssemblerIntegrationFixture fixture;

    // Two near-identical reads to ensure a single strong overlap candidate is found.
    std::string seq0 = randomSequence(3000, 4242);
    std::string seq1 = seq0;
    // Add a tiny number of SNPs to avoid accidental perfect-repeat corner cases, but keep overlap strong.
    seq1[100] = otherBase(seq1[100]);
    seq1[2500] = otherBase(seq1[2500]);

    fixture.createFastq({seq0, seq1});
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);

    // Create chained candidates first.
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 100);

    auto& candidates = fixture.assembler->alignmentCandidates.candidates;
    auto& chainedAlignments = fixture.assembler->alignmentCandidatesAlignmentsData.alignments;
    REQUIRE(candidates.size() == chainedAlignments.size());
    REQUIRE(candidates.size() >= 1);

    // Find the 0-1 candidate.
    uint64_t baseIndex = invalid<uint64_t>;
    for (uint64_t i = 0; i < candidates.size(); ++i) {
        const auto& c = candidates[i];
        if (c.readIds[0] == ReadId(0) && c.readIds[1] == ReadId(1)) {
            baseIndex = i;
            break;
        }
    }
    REQUIRE(baseIndex != invalid<uint64_t>);

    // Duplicate the same (read0,read1) candidate/alignment twice.
    const auto baseCand = candidates[baseIndex];
    const auto baseAln = chainedAlignments[baseIndex];
    candidates.push_back(baseCand);
    chainedAlignments.push_back(baseAln);
    candidates.push_back(baseCand);
    chainedAlignments.push_back(baseAln);

    REQUIRE(candidates.size() == chainedAlignments.size());

    // Compute base-space alignments and evidence; this will create multiple AlignmentData
    // entries for the same read pair.
    fixture.computeAlignments();

    auto countPair = [&](bool onlyActive) -> uint64_t {
        uint64_t n = 0;
        for (const auto& ad : fixture.assembler->alignmentData) {
            const bool isPair =
                (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(1)) ||
                (ad.readIds[0] == ReadId(1) && ad.readIds[1] == ReadId(0));
            if (!isPair) continue;
            if (onlyActive && ad.isDeleted()) continue;
            ++n;
        }
        return n;
    };

    const uint64_t beforeTotal = countPair(false);
    const uint64_t beforeActive = countPair(true);
    REQUIRE(beforeTotal >= 3);
    REQUIRE(beforeActive >= 3);

    fixture.assembler->filterSecondaryAlignmentsPerReadPair(1);

    const uint64_t afterTotal = countPair(false);
    const uint64_t afterActive = countPair(true);
    CHECK(afterTotal == beforeTotal);
    CHECK(afterActive == 1);

    // Ensure redundant alignments are annotated as secondary on both sides.
    uint64_t secondaryDeleted = 0;
    uint64_t secondaryKept = 0;
    for (const auto& ad : fixture.assembler->alignmentData) {
        const bool isPair =
            (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(1)) ||
            (ad.readIds[0] == ReadId(1) && ad.readIds[1] == ReadId(0));
        if (!isPair) continue;
        const bool hasSecondary0 = (ad.deleteReasons0 & AlignmentData::DeleteReasonSecondary) != 0;
        const bool hasSecondary1 = (ad.deleteReasons1 & AlignmentData::DeleteReasonSecondary) != 0;
        if (ad.isDeleted()) {
            CHECK(hasSecondary0);
            CHECK(hasSecondary1);
            ++secondaryDeleted;
        } else if (ad.keptByBothSides()) {
            CHECK_FALSE(hasSecondary0);
            CHECK_FALSE(hasSecondary1);
            ++secondaryKept;
        }
    }
    CHECK(secondaryKept == 1);
    CHECK(secondaryDeleted >= 2);
}

TEST_CASE("Integration: gen_rphase_dp skips same-site transitions", "[integration][dp][parity]") {
    using dinara::hifiasmEcTestHooks::getDpSameSiteComparisons;
    using dinara::hifiasmEcTestHooks::resetDpSameSiteComparisons;

    resetDpSameSiteComparisons();

    AssemblerIntegrationFixture fixture;

    // Build a small multiallelic site on read_0 at position P:
    // - 3 overlaps support alt=C
    // - 3 overlaps support alt=G
    // - 2 overlaps support ref (query base)
    // Each alt has at least one reverse-oriented overlap to avoid strand-bias filtering.
    std::string seq0 = randomSequence(2000, 91231);
    const size_t P = 777;
    seq0[P] = 'A';

    auto makeVariant = [&](char altBase) {
        std::string s = seq0;
        s[P] = altBase;
        return s;
    };

    std::vector<std::string> seqs;
    seqs.push_back(seq0);

    // Ref-supporting overlaps.
    seqs.push_back(seq0);
    seqs.push_back(reverseComplement(seq0));

    // Alt=C (2 forward, 1 reverse)
    seqs.push_back(makeVariant('C'));
    seqs.push_back(makeVariant('C'));
    seqs.push_back(reverseComplement(makeVariant('C')));

    // Alt=G (2 forward, 1 reverse)
    seqs.push_back(makeVariant('G'));
    seqs.push_back(makeVariant('G'));
    seqs.push_back(reverseComplement(makeVariant('G')));

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);

    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    // Run parity EC to trigger gen_rphase_dp.
    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    // With the same-site transition skip, the DP should never attempt to score
    // links between different alleles of the same genomic position.
    CHECK(getDpSameSiteComparisons() == 0);
}

TEST_CASE("Integration: comput_sc_rphase_strict requires both haplotype supports", "[integration][dp][parity][score]") {
    using dinara::hifiasmEcTestHooks::scoreLinkStrict;

    // Two sites, one 64-bit word of candidates.
    const size_t nWords = 1;
    const size_t nSites = 2;
    std::vector<uint64_t> flatBits(nSites * 2 * nWords, 0);
    std::vector<uint64_t> flatAnyBits(nSites * nWords, 0);

    auto setRef = [&](size_t site, uint32_t cand) {
        flatBits[site * 2 * nWords + 2 * 0] |= (1ULL << (cand & 63));
    };
    auto setAlt = [&](size_t site, uint32_t cand) {
        flatBits[site * 2 * nWords + 2 * 0 + 1] |= (1ULL << (cand & 63));
    };

    SECTION("Missing ref-ref support => no link") {
        // alt-alt reads (support hap1 across both)
        setAlt(0, 2); setAlt(1, 2);
        setAlt(0, 3); setAlt(1, 3);
        // ref-only at site 0 (do not cover site 1)
        setRef(0, 0);
        setRef(0, 1);
        // ref-only at site 1 (do not cover site 0)
        setRef(1, 4);
        setRef(1, 5);

        CHECK(scoreLinkStrict(flatBits, flatAnyBits, 0, 1, nWords) == INT64_MIN);
    }

    SECTION("Missing alt-alt support => no link") {
        flatBits.assign(nSites * 2 * nWords, 0);
        // ref-ref reads (support hap0 across both)
        setRef(0, 0); setRef(1, 0);
        setRef(0, 1); setRef(1, 1);
        // alt-only at site 0 (do not cover site 1)
        setAlt(0, 2);
        setAlt(0, 3);
        // alt-only at site 1 (do not cover site 0)
        setAlt(1, 4);
        setAlt(1, 5);

        CHECK(scoreLinkStrict(flatBits, flatAnyBits, 0, 1, nWords) == INT64_MIN);
    }

    SECTION("Both ref-ref and alt-alt support => link score is finite") {
        flatBits.assign(nSites * 2 * nWords, 0);
        // ref-ref
        setRef(0, 0); setRef(1, 0);
        // alt-alt
        setAlt(0, 1); setAlt(1, 1);
        // extra non-covering evidence at one site should not break link
        setRef(0, 2);
        setAlt(1, 3);

        const int64_t sc = scoreLinkStrict(flatBits, flatAnyBits, 0, 1, nWords);
        CHECK(sc != INT64_MIN);
    }

    SECTION("Third allele at only one site => link is rejected") {
        flatBits.assign(nSites * 2 * nWords, 0);
        flatAnyBits.assign(nSites * nWords, 0);

        // Ensure both haplotype supports would exist without the third-allele issue.
        setRef(0, 0); setRef(1, 0);
        setAlt(0, 1); setAlt(1, 1);

        // Candidate 2: "Other" at site 0, Ref at site 1 => should invalidate the entire link.
        flatAnyBits[0] |= (1ULL << 2);
        setRef(1, 2);

        CHECK(scoreLinkStrict(flatBits, flatAnyBits, 0, 1, nWords) == INT64_MIN);
    }

    SECTION("Third allele at both sites counts as ref-ref support (rareRef)") {
        flatBits.assign(nSites * 2 * nWords, 0);
        flatAnyBits.assign(nSites * nWords, 0);

        // Candidate 0: other at both sites => treated as ref-ref support.
        flatAnyBits[0] |= (1ULL << 0);
        flatAnyBits[1] |= (1ULL << 0);

        // Candidate 1: alt-alt support.
        setAlt(0, 1); setAlt(1, 1);

        const int64_t sc = scoreLinkStrict(flatBits, flatAnyBits, 0, 1, nWords);
        CHECK(sc != INT64_MIN);
    }
}

TEST_CASE("Integration: high coveragePeak suppresses DP site retention (singleton site)", "[integration][hifiasm][ec][parity][cc]") {
    AssemblerIntegrationFixture fixture;

    // Build exactly one informative SNP site on read_0 (singleton chain):
    // totalCov=5 (>=5), occ_1=3, occ_0=3 (refCov=2 + query).
    std::string seq0 = randomSequence(2000, 616161);
    const size_t P = 777;
    seq0[P] = 'A';
    auto altC = [&] {
        std::string s = seq0;
        s[P] = 'C';
        return s;
    };

    std::vector<std::string> seqs;
    seqs.push_back(seq0);                 // read_0 (query)
    seqs.push_back(seq0);                 // ref overlap 1
    seqs.push_back(reverseComplement(seq0)); // ref overlap 2 (reverse)
    seqs.push_back(altC());               // alt overlap 1
    seqs.push_back(altC());               // alt overlap 2
    seqs.push_back(reverseComplement(altC())); // alt overlap 3 (reverse)

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    // Force a very high coveragePeak (cc becomes very large in DP).
    // With DP enabled, hifiasm compacts the SNP list to score==1 sites (occ_0>=cc),
    // so this should eliminate SNP sites before trans-closure runs (no overlaps removed).
    fixture.assembler->assemblerInfo->kmerDistributionInfo.coveragePeak = 200;

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    auto alignmentHasSnpAtQueryPos = [&](const AlignmentData& ad, uint32_t pos) -> bool {
        if (ad.info.alignmentId == invalid<size_t>) return false;
        const uint32_t evidenceId = uint32_t(ad.info.alignmentId);
        const bool queryIsRead0 = (ad.readIds[0] == ReadId(0));
        const uint32_t qs = queryIsRead0 ? ad.qs : ad.ts;
        const uint32_t qe = queryIsRead0 ? ad.qe : ad.te;
        if (!(qs <= pos && qe > pos)) return false;
        bool found = false;
        if (queryIsRead0) {
            fixture.assembler->alignedEvidenceStore.forEachSnp1InRange(
                evidenceId, pos, pos + 1,
                [&](uint32_t p, uint8_t /*b*/) { if (p == pos) found = true; }
            );
        } else {
            fixture.assembler->alignedEvidenceStore.forEachSnp0InRange(
                evidenceId, pos, pos + 1,
                [&](uint32_t p, uint8_t /*b*/) { if (p == pos) found = true; }
            );
        }
        return found;
    };

    bool sawAltAtP = false;
    for (const auto& ad : fixture.assembler->alignmentData) {
        if (ad.readIds[0] != ReadId(0) && ad.readIds[1] != ReadId(0)) continue;
        const bool hasP = alignmentHasSnpAtQueryPos(ad, uint32_t(P));
        const bool deletedFromRead0 = (ad.readIds[0] == ReadId(0)) ? ad.isDeleted0() : ad.isDeleted1();
        if (hasP) sawAltAtP = true;
        CHECK_FALSE(deletedFromRead0);
        CHECK_FALSE(ad.coversHetSite);
    }
    REQUIRE(sawAltAtP);
}

TEST_CASE("Integration: high coveragePeak suppresses DP site retention (chained sites)", "[integration][hifiasm][ec][parity][cc]") {
    AssemblerIntegrationFixture fixture;

    // Two SNP sites on read_0 that form a valid DP chain (alt-alt and ref-ref support exist).
    std::string seq0 = randomSequence(2000, 717171);
    const size_t P1 = 500;
    const size_t P2 = 1200;
    seq0[P1] = 'A';
    seq0[P2] = 'A';

    auto altC2 = [&] {
        std::string s = seq0;
        s[P1] = 'C';
        s[P2] = 'C';
        return s;
    };

    std::vector<std::string> seqs;
    seqs.push_back(seq0); // read_0 query

    // Ref-ref overlaps (2 forward, 1 reverse)
    seqs.push_back(seq0);
    seqs.push_back(seq0);
    seqs.push_back(reverseComplement(seq0));

    // Alt-alt overlaps (2 forward, 1 reverse)
    seqs.push_back(altC2());
    seqs.push_back(altC2());
    seqs.push_back(reverseComplement(altC2()));

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    // Make cc very large in DP. With DP enabled, no SNP sites survive into trans-closure.
    fixture.assembler->assemblerInfo->kmerDistributionInfo.coveragePeak = 200;

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    auto alignmentHasSnpAtQueryPos = [&](const AlignmentData& ad, uint32_t pos) -> bool {
        if (ad.info.alignmentId == invalid<size_t>) return false;
        const uint32_t evidenceId = uint32_t(ad.info.alignmentId);
        const bool queryIsRead0 = (ad.readIds[0] == ReadId(0));
        const uint32_t qs = queryIsRead0 ? ad.qs : ad.ts;
        const uint32_t qe = queryIsRead0 ? ad.qe : ad.te;
        if (!(qs <= pos && qe > pos)) return false;
        bool found = false;
        if (queryIsRead0) {
            fixture.assembler->alignedEvidenceStore.forEachSnp1InRange(
                evidenceId, pos, pos + 1,
                [&](uint32_t p, uint8_t /*b*/) { if (p == pos) found = true; }
            );
        } else {
            fixture.assembler->alignedEvidenceStore.forEachSnp0InRange(
                evidenceId, pos, pos + 1,
                [&](uint32_t p, uint8_t /*b*/) { if (p == pos) found = true; }
            );
        }
        return found;
    };

    bool sawAlt = false;
    for (const auto& ad : fixture.assembler->alignmentData) {
        if (ad.readIds[0] != ReadId(0) && ad.readIds[1] != ReadId(0)) continue;
        const bool hasP1 = alignmentHasSnpAtQueryPos(ad, uint32_t(P1));
        const bool hasP2 = alignmentHasSnpAtQueryPos(ad, uint32_t(P2));
        const bool deletedFromRead0 = (ad.readIds[0] == ReadId(0)) ? ad.isDeleted0() : ad.isDeleted1();
        if (hasP1 || hasP2) sawAlt = true;
        CHECK_FALSE(deletedFromRead0);
        CHECK_FALSE(ad.coversHetSite);
    }
    REQUIRE(sawAlt);
}

TEST_CASE("Integration: homopolymer context can suppress singleton DP sites (target HP mask)", "[integration][hifiasm][ec][parity][hpc]") {
    AssemblerIntegrationFixture fixture;

    // Create one SNP site on read_0 at position P where:
    // - query has 'C' (breaks a homopolymer)
    // - alt overlaps have a 6-A homopolymer spanning P and show 'A' at P (HP-suspect)
    // In hifiasm's DP (no-QV) path, singleton sites can be rejected if overlaps are
    // homopolymer-suspect (is_hpc_vec behavior). This suppresses the SNP site before
    // trans-closure runs, so no overlaps are filtered here.
    std::string seq0 = randomSequence(2000, 818181);
    const size_t P = 777;

    // Make the region around P an A-run in the alt reads: positions [P-3, P+2] are 'A'.
    for (int d = -3; d <= 2; ++d) {
        seq0[size_t(int(P) + d)] = 'A';
    }
    // Break the run on the query at the SNP site.
    seq0[P] = 'C';

    auto altA = [&] {
        std::string s = seq0;
        s[P] = 'A'; // restores a 6-A run around P
        return s;
    };

    std::vector<std::string> seqs;
    seqs.push_back(seq0);                    // read_0 query (ref base = C)
    // Ref overlaps: 5 total so occ_0 = refCov + query >= 6 (needed for hifiasm cc cut_bd=6)
    seqs.push_back(seq0);                    // ref overlap 1
    seqs.push_back(reverseComplement(seq0)); // ref overlap 2 (reverse)
    seqs.push_back(seq0);                    // ref overlap 3
    seqs.push_back(seq0);                    // ref overlap 4
    seqs.push_back(reverseComplement(seq0)); // ref overlap 5 (reverse)
    // Alt overlaps (2 forward + 1 reverse) => occ_1=3, totalCov=8
    seqs.push_back(altA());
    seqs.push_back(altA());
    seqs.push_back(reverseComplement(altA()));

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    // cc is hifiasm-style (cut_rate=0.7, cut_bd=6). For small coverage peaks it floors to 6;
    // the extra ref overlaps above ensure the singleton passes occ_0>=cc so the HP check is exercised.
    fixture.assembler->assemblerInfo->kmerDistributionInfo.coveragePeak = 4;

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    auto alignmentHasSnpAtQueryPos = [&](const AlignmentData& ad, uint32_t pos) -> bool {
        if (ad.info.alignmentId == invalid<size_t>) return false;
        const uint32_t evidenceId = uint32_t(ad.info.alignmentId);
        const bool queryIsRead0 = (ad.readIds[0] == ReadId(0));
        const uint32_t qs = queryIsRead0 ? ad.qs : ad.ts;
        const uint32_t qe = queryIsRead0 ? ad.qe : ad.te;
        if (!(qs <= pos && qe > pos)) return false;
        bool found = false;
        if (queryIsRead0) {
            fixture.assembler->alignedEvidenceStore.forEachSnp1InRange(
                evidenceId, pos, pos + 1,
                [&](uint32_t p, uint8_t /*b*/) { if (p == pos) found = true; }
            );
        } else {
            fixture.assembler->alignedEvidenceStore.forEachSnp0InRange(
                evidenceId, pos, pos + 1,
                [&](uint32_t p, uint8_t /*b*/) { if (p == pos) found = true; }
            );
        }
        return found;
    };

    bool sawAltAtP = false;
    for (const auto& ad : fixture.assembler->alignmentData) {
        if (ad.readIds[0] != ReadId(0) && ad.readIds[1] != ReadId(0)) continue;
        const bool hasP = alignmentHasSnpAtQueryPos(ad, uint32_t(P));
        const bool deletedFromRead0 = (ad.readIds[0] == ReadId(0)) ? ad.isDeleted0() : ad.isDeleted1();
        if (hasP) sawAltAtP = true;
        CHECK_FALSE(deletedFromRead0);
        CHECK_FALSE(ad.coversHetSite);
    }
    REQUIRE(sawAltAtP);
}

TEST_CASE("Integration: performHifiasmECParity ignores singleton mismatch evidence for informative-site coverage", "[integration][hifiasm][ec][parity]") {
    AssemblerIntegrationFixture fixture;

    // Query read (0): baseline.
    const size_t readLen = 4000;
    std::string seq0 = randomSequence(readLen, 12345);

    // Construct a robust, informative biallelic SNP with plenty of support,
    // and mixed orientation, to survive strict filters.
    const size_t P = 1500;
    const size_t Q = 500;
    // Ensure the region around P is NOT flagged as HP-suspect by the ONT-style
    // hpc_mask_ff predicate (HPC_PL=12, HPC_RR=4, HPC_CC=2), otherwise the
    // singleton is_hpc_vec step could legitimately reject the site.
    auto ontHpcMask = [&](const std::string& s, size_t pos) -> bool {
        const int64_t len = (int64_t)s.size();
        const int64_t p = (int64_t)pos;
        const int64_t hpc_len = 12;
        const int64_t hpc_min = 4;
        const int64_t hpc_cut = 2;
        const int64_t e = (p + hpc_len <= len) ? (p + hpc_len) : len;
        const int64_t beg = (p >= hpc_len) ? (p - hpc_len) : 0;
        for (int64_t r = 1; r <= hpc_min; ++r) {
            const int64_t rc = r * hpc_cut;
            int64_t k;
            // Including p, forward check
            for (k = p + r; (k < e) && ((k - r) >= beg) && (s[size_t(k)] == s[size_t(k - r)]); ++k) {}
            int64_t ze = k;
            for (k = p - 1; (k >= beg) && ((k + r) < e) && (s[size_t(k)] == s[size_t(k + r)]); --k) {}
            int64_t zs = k + 1;
            if (((ze - zs) > r) && ((ze - zs) >= rc)) return true;

            // Excluding p, forward check
            for (k = p + r + 1; (k < e) && ((k - r) >= beg) && (s[size_t(k)] == s[size_t(k - r)]); ++k) {}
            ze = k; zs = p + 1;
            if (((ze - zs) > r) && ((ze - zs) >= rc)) return true;

            // Including p, reverse check
            for (k = p - r; (k >= beg) && ((k + r) < e) && (s[size_t(k)] == s[size_t(k + r)]); --k) {}
            zs = k + 1;
            for (k = p + 1; (k < e) && ((k - r) >= beg) && (s[size_t(k)] == s[size_t(k - r)]); ++k) {}
            ze = k;
            if (((ze - zs) > r) && ((ze - zs) >= rc)) return true;

            // Excluding p, reverse check
            for (k = p - r - 1; (k >= beg) && ((k + r) < e) && (s[size_t(k)] == s[size_t(k + r)]); --k) {}
            zs = k + 1; ze = p;
            if (((ze - zs) > r) && ((ze - zs) >= rc)) return true;
        }
        return false;
    };

    const size_t windowHalf = 30; // Must be >= HPC_PL
    REQUIRE(P >= windowHalf);
    REQUIRE(P + windowHalf < seq0.size());

    std::string safeWindow;
    const size_t center = windowHalf;
    for (uint32_t attempt = 0; attempt < 5000; ++attempt) {
        safeWindow = randomSequence(2 * windowHalf + 1, 9000 + attempt);
        if (!ontHpcMask(safeWindow, center)) {
            const std::string safeWindowRc = reverseComplement(safeWindow);
            if (!ontHpcMask(safeWindowRc, safeWindowRc.size() - 1 - center)) {
                break;
            }
        }
    }
    REQUIRE(safeWindow.size() == 2 * windowHalf + 1);

    seq0.replace(P - windowHalf, safeWindow.size(), safeWindow);
    REQUIRE_FALSE(ontHpcMask(seq0, P));
    REQUIRE_FALSE(ontHpcMask(reverseComplement(seq0), (readLen - 1 - P)));

    std::string seqAlt = seq0;
    seqAlt[P] = otherBase(seqAlt[P]);

    // Build reads:
    // - 3 forward alt
    // - 3 reverse alt (stored as RC, aligned as reverse strand)
    // - 3 forward ref (one of them is "noisy ref" with a singleton mismatch at Q)
    // - 3 reverse ref (stored as RC)
    std::vector<std::string> seqs;
    seqs.push_back(seq0); // read_0

    std::vector<uint32_t> altReadIndices;
    std::vector<uint32_t> refReadIndices;

    auto addForward = [&](const std::string& s) {
        const uint32_t idx = uint32_t(seqs.size());
        seqs.push_back(s);
        return idx;
    };
    auto addReverse = [&](const std::string& s) {
        const uint32_t idx = uint32_t(seqs.size());
        seqs.push_back(reverseComplement(s));
        return idx;
    };

    for (int i = 0; i < 3; ++i) altReadIndices.push_back(addForward(seqAlt));
    for (int i = 0; i < 3; ++i) altReadIndices.push_back(addReverse(seqAlt));

    // Forward refs.
    refReadIndices.push_back(addForward(seq0));
    refReadIndices.push_back(addForward(seq0));
    std::string seqNoisyRef = seq0;
    seqNoisyRef[Q] = otherBase(seqNoisyRef[Q]);
    const uint32_t noisyRefIndex = addForward(seqNoisyRef);
    refReadIndices.push_back(noisyRefIndex);

    // Reverse refs.
    for (int i = 0; i < 3; ++i) refReadIndices.push_back(addReverse(seq0));

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);

    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    auto countPairAlignments = [&](ReadId a, ReadId b) -> uint64_t {
        uint64_t n = 0;
        for (const auto& ad : fixture.assembler->alignmentData) {
            const bool isPair =
                (ad.readIds[0] == a && ad.readIds[1] == b) ||
                (ad.readIds[0] == b && ad.readIds[1] == a);
            if (isPair) ++n;
        }
        return n;
    };
    REQUIRE(countPairAlignments(ReadId(0), ReadId(altReadIndices.front())) >= 1);
    REQUIRE(countPairAlignments(ReadId(0), ReadId(noisyRefIndex)) >= 1);

    // Sanity-check that the evidence store actually contains the engineered SNP at P
    // for an alt overlap, and the singleton mismatch at Q for the noisy ref overlap.
    auto baseToInt = [&](char c) -> uint8_t {
        if (c == 'A') return 0;
        if (c == 'C') return 1;
        if (c == 'G') return 2;
        if (c == 'T') return 3;
        return 4;
    };

    // Find whether there exists an alignment for (read_0, other) that covers `pos`
    // and has a SNP at `pos` in the query coordinate, and return the strand flag
    // used by parity EC (`isRev = !isSameStrand`).
    auto findSnpSupportAt = [&](uint32_t otherReadIdx, uint32_t pos, bool& isRevOut) -> bool {
        for (const auto& ad : fixture.assembler->alignmentData) {
            const bool isPair =
                (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(otherReadIdx)) ||
                (ad.readIds[0] == ReadId(otherReadIdx) && ad.readIds[1] == ReadId(0));
            if (!isPair) continue;
            if (ad.info.alignmentId == invalid<size_t>) continue;

            const bool queryIsRead0 = (ad.readIds[0] == ReadId(0));
            const uint32_t qs = queryIsRead0 ? ad.qs : ad.ts;
            const uint32_t qe = queryIsRead0 ? ad.qe : ad.te;
            if (!(qs <= pos && qe > pos)) continue;

            bool found = false;
            uint8_t base = 4;
            if (queryIsRead0) {
                fixture.assembler->alignedEvidenceStore.forEachSnp1InRange(
                    uint32_t(ad.info.alignmentId), pos, pos + 1,
                    [&](uint32_t p, uint8_t b) { if (p == pos) { found = true; base = b; } }
                );
            } else {
                fixture.assembler->alignedEvidenceStore.forEachSnp0InRange(
                    uint32_t(ad.info.alignmentId), pos, pos + 1,
                    [&](uint32_t p, uint8_t b) { if (p == pos) { found = true; base = b; } }
                );
            }
            if (found) {
                REQUIRE(base < 4);
                REQUIRE(base != baseToInt(seqs[0][pos]));
                isRevOut = !ad.isSameStrand;
                return true;
            }
        }
        return false;
    };

    // Identify which overlaps actually contribute the engineered mismatch at P.
    std::vector<uint32_t> altSupportAtP;
    bool hasAltFwd = false;
    bool hasAltRev = false;
    for (uint32_t idx : altReadIndices) {
        bool isRev = false;
        if (findSnpSupportAt(idx, uint32_t(P), isRev)) {
            altSupportAtP.push_back(idx);
            hasAltRev |= isRev;
            hasAltFwd |= !isRev;
        }
    }
    REQUIRE(altSupportAtP.size() >= 3);
    REQUIRE(hasAltFwd);
    REQUIRE(hasAltRev);

    // The noisy ref overlap must contain the singleton mismatch at Q and must NOT
    // contain the engineered SNP at P.
    {
        bool isRev = false;
        REQUIRE(findSnpSupportAt(noisyRefIndex, uint32_t(Q), isRev));
        CHECK_FALSE(findSnpSupportAt(noisyRefIndex, uint32_t(P), isRev));
    }

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    // In this synthetic setup, the only potential informative SNP is at position P.
    // Verify that `coversHetSite` matches that position (no off-by-one / wrong projection).
    auto alignmentHasSnpAtQueryPos = [&](const AlignmentData& ad, uint32_t pos) -> bool {
        if (ad.info.alignmentId == invalid<size_t>) return false;
        const uint32_t evidenceId = uint32_t(ad.info.alignmentId);

        const bool queryIsRead0 = (ad.readIds[0] == ReadId(0));
        const uint32_t qs = queryIsRead0 ? ad.qs : ad.ts;
        const uint32_t qe = queryIsRead0 ? ad.qe : ad.te;
        if (!(qs <= pos && qe > pos)) return false;

        bool found = false;
        if (queryIsRead0) {
            fixture.assembler->alignedEvidenceStore.forEachSnp1InRange(
                evidenceId, pos, pos + 1,
                [&](uint32_t p, uint8_t /*b*/) { if (p == pos) found = true; }
            );
        } else {
            fixture.assembler->alignedEvidenceStore.forEachSnp0InRange(
                evidenceId, pos, pos + 1,
                [&](uint32_t p, uint8_t /*b*/) { if (p == pos) found = true; }
            );
        }
        return found;
    };

    bool foundAnyAltMarked = false;
    for (const auto& ad : fixture.assembler->alignmentData) {
        if (ad.readIds[0] != ReadId(0) && ad.readIds[1] != ReadId(0)) continue;
        const bool hasSnpAtP = alignmentHasSnpAtQueryPos(ad, uint32_t(P));
        if (ad.coversHetSite) {
            CHECK(hasSnpAtP);
        }
        if (hasSnpAtP && ad.coversHetSite) {
            foundAnyAltMarked = true;
        }
    }
    REQUIRE(foundAnyAltMarked);

    // The noisy ref overlap has only a singleton mismatch at Q (pos != P), so it must not be
    // mapped to an informative SNP row at P and must not get coversHetSite.
    bool noisyHasCoversHetSite = false;
    for (const auto& ad : fixture.assembler->alignmentData) {
        const bool isPair =
            (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(noisyRefIndex)) ||
            (ad.readIds[0] == ReadId(noisyRefIndex) && ad.readIds[1] == ReadId(0));
        if (!isPair) continue;
        noisyHasCoversHetSite |= ad.coversHetSite;
        CHECK_FALSE(alignmentHasSnpAtQueryPos(ad, uint32_t(P)));
        CHECK(alignmentHasSnpAtQueryPos(ad, uint32_t(Q)));
    }
    CHECK_FALSE(noisyHasCoversHetSite);
}

TEST_CASE("Integration: performHifiasmECParity removes trans overlaps (ALT-supporting) like hifiasm", "[integration][hifiasm][ec][parity][trans]") {
    AssemblerIntegrationFixture fixture;

    const size_t readLen = 2500;
    std::string seq0 = randomSequence(readLen, 424242);
    const size_t P = 1200;

    // Ensure the region around P is NOT flagged by the ONT-style hpc_mask_ff predicate
    // for both the reference and alternative base, and for both strands.
    auto ontHpcMask = [&](const std::string& s, size_t pos) -> bool {
        const int64_t len = (int64_t)s.size();
        const int64_t p = (int64_t)pos;
        const int64_t hpc_len = 12;
        const int64_t hpc_min = 4;
        const int64_t hpc_cut = 2;
        const int64_t e = (p + hpc_len <= len) ? (p + hpc_len) : len;
        const int64_t beg = (p >= hpc_len) ? (p - hpc_len) : 0;
        for (int64_t r = 1; r <= hpc_min; ++r) {
            const int64_t rc = r * hpc_cut;
            int64_t k;
            for (k = p + r; (k < e) && ((k - r) >= beg) && (s[size_t(k)] == s[size_t(k - r)]); ++k) {}
            int64_t ze = k;
            for (k = p - 1; (k >= beg) && ((k + r) < e) && (s[size_t(k)] == s[size_t(k + r)]); --k) {}
            int64_t zs = k + 1;
            if (((ze - zs) > r) && ((ze - zs) >= rc)) return true;

            for (k = p + r + 1; (k < e) && ((k - r) >= beg) && (s[size_t(k)] == s[size_t(k - r)]); ++k) {}
            ze = k; zs = p + 1;
            if (((ze - zs) > r) && ((ze - zs) >= rc)) return true;

            for (k = p - r; (k >= beg) && ((k + r) < e) && (s[size_t(k)] == s[size_t(k + r)]); --k) {}
            zs = k + 1;
            for (k = p + 1; (k < e) && ((k - r) >= beg) && (s[size_t(k)] == s[size_t(k - r)]); ++k) {}
            ze = k;
            if (((ze - zs) > r) && ((ze - zs) >= rc)) return true;

            for (k = p - r - 1; (k >= beg) && ((k + r) < e) && (s[size_t(k)] == s[size_t(k + r)]); --k) {}
            zs = k + 1; ze = p;
            if (((ze - zs) > r) && ((ze - zs) >= rc)) return true;
        }
        return false;
    };

    const size_t windowHalf = 30;
    REQUIRE(P >= windowHalf);
    REQUIRE(P + windowHalf < seq0.size());

    std::string safeWindow;
    for (uint32_t attempt = 0; attempt < 10000; ++attempt) {
        safeWindow = randomSequence(2 * windowHalf + 1, 99000 + attempt);
        std::string safeWindowAlt = safeWindow;
        safeWindowAlt[windowHalf] = otherBase(safeWindowAlt[windowHalf]);

        if (ontHpcMask(safeWindow, windowHalf)) continue;
        if (ontHpcMask(safeWindowAlt, windowHalf)) continue;

        const std::string rc = reverseComplement(safeWindow);
        const std::string rcAlt = reverseComplement(safeWindowAlt);
        const size_t rcPos = rc.size() - 1 - windowHalf;
        if (ontHpcMask(rc, rcPos)) continue;
        if (ontHpcMask(rcAlt, rcPos)) continue;
        break;
    }
    REQUIRE(safeWindow.size() == 2 * windowHalf + 1);

    seq0.replace(P - windowHalf, safeWindow.size(), safeWindow);
    REQUIRE_FALSE(ontHpcMask(seq0, P));

    std::string seqAlt = seq0;
    seqAlt[P] = otherBase(seqAlt[P]);
    REQUIRE_FALSE(ontHpcMask(seqAlt, P));

    // Build reads:
    // With hifiasm-parity cc (cut_bd=6), we need enough ref support so the singleton
    // SNP can satisfy occ_0 >= 6 (refCov+query >= 6). Use 5 ref overlaps and 3 alt overlaps.
    std::vector<std::string> seqs;
    seqs.push_back(seq0); // read_0 query
    const uint32_t ref1 = uint32_t(seqs.size()); seqs.push_back(seq0);
    const uint32_t ref2 = uint32_t(seqs.size()); seqs.push_back(seq0);
    const uint32_t ref3 = uint32_t(seqs.size()); seqs.push_back(reverseComplement(seq0));
    const uint32_t ref4 = uint32_t(seqs.size()); seqs.push_back(seq0);
    const uint32_t ref5 = uint32_t(seqs.size()); seqs.push_back(reverseComplement(seq0));

    const uint32_t alt1 = uint32_t(seqs.size()); seqs.push_back(seqAlt);
    const uint32_t alt2 = uint32_t(seqs.size()); seqs.push_back(seqAlt);
    const uint32_t alt3 = uint32_t(seqs.size()); seqs.push_back(reverseComplement(seqAlt));

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    // cc is computed like hifiasm (cut_rate=0.7, cut_bd=6) from coveragePeak; for low peaks
    // it still floors to 6, so coverage must come from overlaps (already ensured above).
    fixture.assembler->assemblerInfo->kmerDistributionInfo.coveragePeak = 4;

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    auto expectDeletedFromRead0 = [&](uint32_t other, bool shouldDelete) {
        const auto& ads = fixture.assembler->alignmentData;
        const AlignmentData* found = nullptr;
        for (const auto& ad : ads) {
            const bool matchesPair =
                (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(other)) ||
                (ad.readIds[0] == ReadId(other) && ad.readIds[1] == ReadId(0));
            if (!matchesPair) continue;
            found = &ad;
            break;
        }
        REQUIRE(found != nullptr);

        const bool deletedFromRead0 = (found->readIds[0] == ReadId(0)) ? found->isDeleted0() : found->isDeleted1();
        CHECK(deletedFromRead0 == shouldDelete);
        CHECK(found->coversHetSite == shouldDelete); // only ALT-supporting overlaps should have mismatch evidence at P
    };

    expectDeletedFromRead0(ref1, false);
    expectDeletedFromRead0(ref2, false);
    expectDeletedFromRead0(ref3, false);
    expectDeletedFromRead0(ref4, false);
    expectDeletedFromRead0(ref5, false);

    expectDeletedFromRead0(alt1, true);
    expectDeletedFromRead0(alt2, true);
    expectDeletedFromRead0(alt3, true);
}

TEST_CASE("Integration: performHifiasmECParity drops adjacent SNP sites like hifiasm", "[integration][hifiasm][ec][parity][adjacent]") {
    AssemblerIntegrationFixture fixture;

    const size_t readLen = 2500;
    std::string seq0 = randomSequence(readLen, 515151);
    const size_t P = 1200;
    const size_t P2 = P + 1;
    REQUIRE(P2 < seq0.size());

    // Engineer a clean adjacent-variant pair.
    seq0[P] = 'A';
    seq0[P2] = 'A';
    std::string seqAlt = seq0;
    seqAlt[P] = 'C';
    seqAlt[P2] = 'C';

    // Build reads:
    // With hifiasm-parity cc (cut_bd=6), ensure there is enough ref support so sites
    // would otherwise validate (then adjacency filtering is the reason no overlaps are removed).
    std::vector<std::string> seqs;
    seqs.push_back(seq0); // read_0 query
    seqs.push_back(seq0);
    seqs.push_back(seq0);
    seqs.push_back(reverseComplement(seq0));
    seqs.push_back(seq0);
    seqs.push_back(reverseComplement(seq0));
    seqs.push_back(seqAlt);
    seqs.push_back(seqAlt);
    seqs.push_back(reverseComplement(seqAlt));

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    fixture.assembler->assemblerInfo->kmerDistributionInfo.coveragePeak = 4;

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    // Adjacent sites should be dropped at the HiFi overlap-marking stage, so
    // the read becomes non-informative and no overlaps are removed.
    for (const auto& ad : fixture.assembler->alignmentData) {
        if (ad.isDeleted()) continue;
        if (ad.readIds[0] != ReadId(0) && ad.readIds[1] != ReadId(0)) continue;
        const bool deletedFromRead0 = (ad.readIds[0] == ReadId(0)) ? ad.isDeleted0() : ad.isDeleted1();
        CHECK_FALSE(deletedFromRead0);
        CHECK_FALSE(ad.coversHetSite);
    }
}

TEST_CASE("Integration: performHifiasmECParity multi_check rescues weak sites and removes trans overlaps", "[integration][hifiasm][ec][parity][multicheck]") {
    AssemblerIntegrationFixture fixture;

    const size_t readLen = 1200;
    std::string seq0 = randomSequence(readLen, 616161);

    // Create many weak SNPs supported by exactly 2 ALT overlaps out of 5 total overlaps
    // (occ_1=2, totalCov=5 => below informative threshold 3), plus 2 isolated sites
    // that survive the multi_check +/-32bp neighbor filter.
    std::vector<size_t> variantSites;
    variantSites.push_back(50);   // isolated
    variantSites.push_back(1100); // isolated
    for (size_t p = 300; p < 610; p += 5) { // dense cluster (will be discarded by +/-32 filter)
        variantSites.push_back(p);
    }

    for (const size_t p : variantSites) {
        seq0[p] = 'A';
    }
    std::string seqAlt = seq0;
    for (const size_t p : variantSites) {
        seqAlt[p] = 'C';
    }

    // Reads: query + 5 ref overlaps + 2 alt overlaps = totalCov=7 at all variant sites.
    // Ensure occ_0 = refCov + query >= 6 so DP (cc cut_bd=6) keeps these sites and
    // multi_check can run on the DP-retained SNP list.
    std::vector<std::string> seqs;
    seqs.push_back(seq0); // read_0 query
    
    const uint32_t ref1 = uint32_t(seqs.size()); seqs.push_back(seq0);
    const uint32_t ref2 = uint32_t(seqs.size()); seqs.push_back(seq0);
    const uint32_t ref3 = uint32_t(seqs.size()); seqs.push_back(reverseComplement(seq0)); // reverse
    const uint32_t ref4 = uint32_t(seqs.size()); seqs.push_back(seq0);
    const uint32_t ref5 = uint32_t(seqs.size()); seqs.push_back(reverseComplement(seq0)); // reverse

    const uint32_t alt1 = uint32_t(seqs.size()); seqs.push_back(seqAlt);
    const uint32_t alt2 = uint32_t(seqs.size()); seqs.push_back(seqAlt);

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    // cc is hifiasm-style (cut_bd=6), but multi_check is independent of DP validation and
    // should still rescue repeated weak sites.
    fixture.assembler->assemblerInfo->kmerDistributionInfo.coveragePeak = 4;

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    auto expectDeletedFromRead0 = [&](uint32_t other, bool shouldDelete) {
        const auto& ads = fixture.assembler->alignmentData;
        const AlignmentData* found = nullptr;
        for (const auto& ad : ads) {
            const bool matchesPair =
                (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(other)) ||
                (ad.readIds[0] == ReadId(other) && ad.readIds[1] == ReadId(0));
            if (!matchesPair) continue;
            found = &ad;
            break;
        }
        REQUIRE(found != nullptr);

        const bool deletedFromRead0 = (found->readIds[0] == ReadId(0)) ? found->isDeleted0() : found->isDeleted1();
        CHECK(deletedFromRead0 == shouldDelete);
    };

    // Multi_check should validate weak sites and then remove ALT overlaps (TRANS).
    expectDeletedFromRead0(ref1, false);
    expectDeletedFromRead0(ref2, false);
    expectDeletedFromRead0(ref3, false);
    expectDeletedFromRead0(ref4, false);
    expectDeletedFromRead0(ref5, false);
    expectDeletedFromRead0(alt1, true);
    expectDeletedFromRead0(alt2, true);
}

TEST_CASE("Integration: multi_check does not use ref support from trans overlaps", "[integration][hifiasm][ec][parity][multicheck][ciscounts]") {
    AssemblerIntegrationFixture fixture;

    const size_t readLen = 1600;
    std::string seq0 = randomSequence(readLen, 717272);

    // Two strong sites (X,Y) split overlaps into ALT (trans) and REF (cis).
    const size_t X = 400;
    const size_t Y = 450;
    seq0[X] = 'A';
    seq0[Y] = 'A';

    // Many weak sites that are ALT only on 2 overlaps, REF on trans overlaps (3 overlaps).
    // Globally these sites have occ_1=2 and occ_0 boosted by trans ref support; after hifiasm-style
    // trans-closure "not real allele" decrements, occ_0 drops to 1 (query only) and these sites must
    // NOT be rescued by multi_check.
    std::vector<size_t> weakSites;
    weakSites.push_back(900);   // isolated (keep outside the strong-site helper overlaps)
    weakSites.push_back(1200);  // isolated
    for (size_t p = 520; p < 720; p += 5) weakSites.push_back(p); // dense cluster to satisfy up=0.04

    std::string seqCis = seq0;
    for (size_t p : weakSites) {
        seq0[p] = 'A';
        seqCis[p] = 'C'; // ALT only on cis overlaps
    }
    // Ensure (X,Y) differ only on trans overlaps.
    seqCis[X] = 'A';
    seqCis[Y] = 'A';

    std::string seqTrans = seq0;
    seqTrans[X] = 'C'; // ALT at validated site
    seqTrans[Y] = 'C'; // ALT at validated site

    std::vector<std::string> seqs;
    seqs.push_back(seq0); // read_0 query

    // Add extra REF overlaps at the strong sites (X,Y) using contained prefix reads.
    // They overlap the strong sites but do not cover the weak-site cluster (>=520),
    // so they raise occ_0 at X/Y without contributing REF support at weak sites.
    const size_t strongPrefixLen = 500; // covers X=400,Y=450 but excludes weakSites cluster starting at 520
    REQUIRE(strongPrefixLen <= seq0.size());
    const std::string seqPrefix = seq0.substr(0, strongPrefixLen);
    seqs.push_back(seqPrefix);
    seqs.push_back(seqPrefix);
    seqs.push_back(reverseComplement(seqPrefix));

    // Two overlaps that carry ALT at weak sites, REF at (X,Y).
    const uint32_t cis1 = uint32_t(seqs.size()); seqs.push_back(seqCis);
    const uint32_t cis2 = uint32_t(seqs.size()); seqs.push_back(seqCis);

    // TRANS overlaps (5) - ALT at X, REF at weak sites.
    // Having 5 trans overlaps ensures weak sites start with occ_0 = 6 (refCov+query),
    // so they survive DP but are removed from multi_check after occ_0 decrements to 1.
    const uint32_t tr1 = uint32_t(seqs.size()); seqs.push_back(seqTrans);
    const uint32_t tr2 = uint32_t(seqs.size()); seqs.push_back(seqTrans);
    const uint32_t tr3 = uint32_t(seqs.size()); seqs.push_back(reverseComplement(seqTrans)); // reverse
    const uint32_t tr4 = uint32_t(seqs.size()); seqs.push_back(seqTrans);
    const uint32_t tr5 = uint32_t(seqs.size()); seqs.push_back(reverseComplement(seqTrans)); // reverse

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    fixture.assembler->assemblerInfo->kmerDistributionInfo.coveragePeak = 4;

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    // TRANS overlaps must be removed (due to validated site X).
    auto deletedFromRead0 = [&](uint32_t other) -> bool {
        for (const auto& ad : fixture.assembler->alignmentData) {
            const bool matchesPair =
                (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(other)) ||
                (ad.readIds[0] == ReadId(other) && ad.readIds[1] == ReadId(0));
            if (!matchesPair) continue;
            return (ad.readIds[0] == ReadId(0)) ? ad.isDeleted0() : ad.isDeleted1();
        }
        REQUIRE(false);
        return false;
    };

    CHECK_FALSE(deletedFromRead0(cis1));
    CHECK_FALSE(deletedFromRead0(cis2));
    CHECK(deletedFromRead0(tr1));
    CHECK(deletedFromRead0(tr2));
    CHECK(deletedFromRead0(tr3));
    CHECK(deletedFromRead0(tr4));
    CHECK(deletedFromRead0(tr5));
}

TEST_CASE("Integration: performHifiasmECParity excludes gap overlaps from SNP-site coverage", "[integration][hifiasm][ec][parity][gap][coverage]") {
    AssemblerIntegrationFixture fixture;

    const size_t readLen = 2200;
    const size_t P = 1100;
    const size_t delLen = 10;
    REQUIRE(P + delLen < readLen);

    // Build a sequence where the deleted segment (starting at P) is unique, to make
    // gap placement deterministic. Enforce the SNP base first so uniqueness is checked
    // against the actual sequence used for alignment.
    std::string seq0;
    for (uint32_t attempt = 0; attempt < 500; ++attempt) {
        seq0 = randomSequence(readLen, 888000 + attempt);
        seq0[P] = 'A';
        const std::string delSeq = seq0.substr(P, delLen);
        const size_t first = seq0.find(delSeq);
        const size_t second = seq0.find(delSeq, first + 1);
        if (first == P && second == std::string::npos) break;
    }
    REQUIRE(seq0.size() == readLen);
    std::string seqAlt = seq0;
    seqAlt[P] = 'C';

    // One overlap has a deletion relative to the query (so the query has an insertion / gap on coverage).
    std::string seqGap = seq0;
    seqGap.erase(P, delLen);

    // Reads:
    // - 3 ALT overlaps (mismatch at P) => occ_1=3
    // - 1 REF overlap
    // - 1 GAP overlap: does NOT provide a base at P and must not count toward ref coverage.
    std::vector<std::string> seqs;
    seqs.push_back(seq0);                    // read_0 query
    seqs.push_back(seqAlt);                  // read_1 alt
    seqs.push_back(seqAlt);                  // read_2 alt
    seqs.push_back(seqAlt);                  // read_3 alt
    seqs.push_back(seq0);                    // read_4 ref
    seqs.push_back(seqGap);                  // read_5 gap (shorter)

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    // Verify we have an insertion (gap) event on the query side for the (read_0, read_5) overlap that spans P.
    const AlignmentData* gapAd = nullptr;
    for (const auto& ad : fixture.assembler->alignmentData) {
        const bool matchesPair =
            (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(5)) ||
            (ad.readIds[0] == ReadId(5) && ad.readIds[1] == ReadId(0));
        if (!matchesPair) continue;
        gapAd = &ad;
        break;
    }
    REQUIRE(gapAd != nullptr);
    REQUIRE(gapAd->info.alignmentId != invalid<size_t>);

    const auto& store = fixture.assembler->alignedEvidenceStore;
    const uint32_t evidenceId = uint32_t(gapAd->info.alignmentId);
    span<const IndelEvidence> gapIndels =
        (gapAd->readIds[0] == ReadId(0)) ? store.getIndels1(evidenceId) : store.getIndels0(evidenceId);
    bool spansP = false;
    for (const auto& e : gapIndels) {
        // In Dinara's query-coordinate indel stream, type==1 corresponds to query bases
        // that have no aligned partner (gap in the other read), i.e. a coverage hole.
        if (!e.isDeletion()) continue;
        const uint32_t b = e.pos();
        const uint32_t epos = b + e.len();
        if (b <= P && P < epos) {
            spansP = true;
            break;
        }
    }
    REQUIRE(spansP);

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    // With correct gap-aware coverage, P is not informative and the read is non-informative,
    // so no overlaps are deleted and none are marked as covering het sites.
    for (const auto& ad : fixture.assembler->alignmentData) {
        if (ad.readIds[0] != ReadId(0) && ad.readIds[1] != ReadId(0)) continue;
        const bool deletedFromRead0 = (ad.readIds[0] == ReadId(0)) ? ad.isDeleted0() : ad.isDeleted1();
        CHECK_FALSE(deletedFromRead0);
        CHECK_FALSE(ad.coversHetSite);
    }
}

TEST_CASE("Integration: performHifiasmECParity filters SV/large-indel overlaps like hifiasm", "[integration][hifiasm][ec][parity][sv]") {
    AssemblerIntegrationFixture fixture;

    const size_t readLen = 3000;
    const size_t P = 1500;
    const size_t delLen = 30; // >= SV_MIN_LEN (20)
    REQUIRE(P + delLen < readLen);

    // Make the deleted segment unique so the indel anchors at the intended position.
    std::string seq0;
    for (uint32_t attempt = 0; attempt < 500; ++attempt) {
        seq0 = randomSequence(readLen, 991000 + attempt);
        const std::string delSeq = seq0.substr(P, delLen);
        const size_t first = seq0.find(delSeq);
        const size_t second = seq0.find(delSeq, first + 1);
        if (first == P && second == std::string::npos) break;
    }
    REQUIRE(seq0.size() == readLen);

    // ALT allele overlaps: delete a 30bp segment (query has an insertion relative to these overlaps).
    std::string seqAlt = seq0;
    seqAlt.erase(P, delLen);

    std::vector<std::string> seqs;
    seqs.push_back(seq0); // read_0 query

    // 3 ALT overlaps (one reverse) and 3 REF overlaps (one reverse) to satisfy strand-balance checks.
    const uint32_t alt1 = uint32_t(seqs.size()); seqs.push_back(seqAlt);
    const uint32_t alt2 = uint32_t(seqs.size()); seqs.push_back(seqAlt);
    const uint32_t alt3 = uint32_t(seqs.size()); seqs.push_back(reverseComplement(seqAlt)); // reverse

    const uint32_t ref1 = uint32_t(seqs.size()); seqs.push_back(seq0);
    const uint32_t ref2 = uint32_t(seqs.size()); seqs.push_back(seq0);
    const uint32_t ref3 = uint32_t(seqs.size()); seqs.push_back(reverseComplement(seq0)); // reverse

    fixture.createFastq(seqs);
    fixture.initAssembler();
    fixture.loadReads();
    fixture.generateMarkers(16, 5);
    fixture.countKmers();
    fixture.applyFilter(1, 1000);
    fixture.buildIndex();
    fixture.chainCandidates(0.1, 200);
    fixture.computeAlignments();

    auto requireLargeIndelEvidenceNearP = [&](uint32_t other) {
        const AlignmentData* adPtr = nullptr;
        for (const auto& ad : fixture.assembler->alignmentData) {
            const bool matchesPair =
                (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(other)) ||
                (ad.readIds[0] == ReadId(other) && ad.readIds[1] == ReadId(0));
            if (!matchesPair) continue;
            adPtr = &ad;
            break;
        }
        REQUIRE(adPtr != nullptr);
        REQUIRE(adPtr->info.alignmentId != invalid<size_t>);

        const auto& store = fixture.assembler->alignedEvidenceStore;
        const uint32_t evidenceId = uint32_t(adPtr->info.alignmentId);
        span<const IndelEvidence> indels =
            (adPtr->readIds[0] == ReadId(0)) ? store.getIndels1(evidenceId) : store.getIndels0(evidenceId);

        uint32_t sumLenInWindow = 0;
        for (const auto& e : indels) {
            const uint32_t pos = e.pos();
            if (pos + 100 < P) continue;
            if (pos > P + 100) continue;
            sumLenInWindow += e.len();
        }
        CAPTURE(other);
        CAPTURE(sumLenInWindow);
        REQUIRE(sumLenInWindow >= delLen);
    };

    requireLargeIndelEvidenceNearP(alt1);
    requireLargeIndelEvidenceNearP(alt2);
    requireLargeIndelEvidenceNearP(alt3);

    withSilencedIoInDir(fixture.testDir, [&] {
        fixture.assembler->performHifiasmECParity(1);
    });

    auto deletedFromRead0 = [&](uint32_t other) -> bool {
        for (const auto& ad : fixture.assembler->alignmentData) {
            const bool matchesPair =
                (ad.readIds[0] == ReadId(0) && ad.readIds[1] == ReadId(other)) ||
                (ad.readIds[0] == ReadId(other) && ad.readIds[1] == ReadId(0));
            if (!matchesPair) continue;
            return (ad.readIds[0] == ReadId(0)) ? ad.isDeleted0() : ad.isDeleted1();
        }
        REQUIRE(false);
        return false;
    };

    // SV allele overlaps should be filtered (trans), ref overlaps kept.
    CHECK(deletedFromRead0(alt1));
    CHECK(deletedFromRead0(alt2));
    CHECK(deletedFromRead0(alt3));

    CHECK_FALSE(deletedFromRead0(ref1));
    CHECK_FALSE(deletedFromRead0(ref2));
    CHECK_FALSE(deletedFromRead0(ref3));
}
