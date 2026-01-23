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
    const std::string safeWindow = "ACGTTGCACTGATCGTACGTA"; // 21bp
    REQUIRE(P >= 10);
    REQUIRE(P + 10 < seq0.size());
    seq0.replace(P - 10, safeWindow.size(), safeWindow);

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
