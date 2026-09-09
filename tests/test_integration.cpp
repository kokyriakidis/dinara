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
#include "../src/MarkerGraph.hpp"
#include "../src/mode3-DirectedAnchorGraph.hpp"

// Standard library
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <string>
#include <iostream>
#include <limits>
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

AlignOptions makeDefaultAlignOptionsForIntegration()
{
    AlignOptions options;
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
    options.maxErrorRate = 0.3;
    // Base-level overlap DP scoring defaults (hifiasm/minimap2 parity).
    options.overlapDpMatchScore = 2;
    options.overlapDpMismatchScore = -4;
    options.overlapDpGapOpen1 = 4;
    options.overlapDpGapExtend1 = 2;
    options.overlapDpGapOpen2 = 24;
    options.overlapDpGapExtend2 = 1;
    return options;
}

} // namespace

TEST_CASE("MarkerGraph vertex coverage histogram counts canonical vertices", "[markerGraph][histogram]")
{
    MarkerGraph markerGraph;
    markerGraph.constructVertices();
    markerGraph.vertices().createNew("", 4096);

    // Create 4 vertices with coverages: [3,3,1,2].
    markerGraph.vertices().appendVector(3);
    markerGraph.vertices().appendVector(3);
    markerGraph.vertices().appendVector(1);
    markerGraph.vertices().appendVector(2);

    // Reverse-complement pairs: (0,1) and (2,3).
    markerGraph.reverseComplementVertex.createNew("", 4096);
    markerGraph.reverseComplementVertex.resize(4);
    markerGraph.reverseComplementVertex[0] = 1;
    markerGraph.reverseComplementVertex[1] = 0;
    markerGraph.reverseComplementVertex[2] = 3;
    markerGraph.reverseComplementVertex[3] = 2;

    const auto all = markerGraph.computeVertexCoverageHistogram(false);
    REQUIRE(all.size() >= 4);
    CHECK(all[1] == 1);
    CHECK(all[2] == 1);
    CHECK(all[3] == 2);

    const auto canonical = markerGraph.computeVertexCoverageHistogram(true);
    REQUIRE(canonical.size() >= 4);
    CHECK(canonical[1] == 1); // vertex 2
    CHECK(canonical[2] == 0);
    CHECK(canonical[3] == 1); // vertex 0
}

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
        // The simd-minimizers backend requires the window length l = k + w - 1
        // to be odd (k and w same parity). Nudge w up by one when needed so the
        // various (k, s) combinations used by these tests stay valid.
        const int w = ((k + s) % 2 == 0) ? s : s + 1;
        withSilencedIoInDir(testDir, [&] {
            assembler->findMarkersSimdMinimizers(1, k, w, /*useHifiasm*/ false);
        });
    }
    
    void countKmers() {
        withSilencedIoInDir(testDir, [&] { assembler->countKmersFromMarkerKmerIds(1); });
    }
    
    void applyFilter(uint64_t minFreq, uint64_t maxFreq) {
        withSilencedIoInDir(testDir, [&] { assembler->applyKmerCountFilter(minFreq, maxFreq, 1); });
    }
};

TEST_CASE("ReadGraph.filterSecondaryRequireNonRedundantOnBothReads removes symmetric-only overlaps",
    "[integration][readgraph][secondary][symmetry]")
{
    AssemblerIntegrationFixture fixture;

    const std::string r0 = randomSequence(3000, 9301);
    const std::string r1 = randomSequence(3000, 9302);
    fixture.createFastq({r0, r1});
    fixture.initAssembler();
    fixture.loadReads();

    fixture.assembler->alignmentData.createNew("", 4096);

    AlignmentInfo infoA;
    infoA.dpScore = 1000;
    infoA.mismatchCount = 0;
    infoA.gapCount = 0;
    infoA.gapEventCount = 0;
    AlignmentData adA(array<ReadId, 2>{ReadId(0), ReadId(1)}, true, infoA);
    adA.qs = 0;
    adA.qe = 1000;
    adA.ts = 500;
    adA.te = 1500;
    fixture.assembler->alignmentData.push_back(adA);

    AlignmentInfo infoB = infoA;
    infoB.dpScore = 900; // within 80% of best
    AlignmentData adB(array<ReadId, 2>{ReadId(0), ReadId(1)}, true, infoB);
    adB.qs = 1600;
    adB.qe = 2600;
    // Highly redundant on the partner read interval, but non-overlapping on r0.
    adB.ts = 520;
    adB.te = 1520;
    fixture.assembler->alignmentData.push_back(adB);

    fixture.assembler->computeAlignmentTableForTesting();

    // Legacy behavior (r0-only redundancy): both are kept.
    fixture.assembler->filterSecondaryAlignmentsPerReadPair(1, false);
    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted0());
    CHECK_FALSE(fixture.assembler->alignmentData[1].isDeleted0());

    // Symmetric behavior: second overlap is removed due to partner-read redundancy only.
    fixture.assembler->filterSecondaryAlignmentsPerReadPair(1, true);
    CHECK_FALSE(fixture.assembler->alignmentData[0].isDeleted0());
    CHECK(fixture.assembler->alignmentData[1].isDeleted0());
    CHECK(fixture.assembler->getRemovedSecondaryAlignmentBySymmetryOnlyCountForTesting() == 1);
}

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




















namespace {

uint64_t addDagSimpleNode(
    dinara::mode3::DirectedAnchorGraph& dag,
    uint64_t markerId,
    uint64_t lengthBp,
    double coverage)
{
    dinara::mode3::DagNodeInfo info;
    info.anchorChain.push_back(markerId);
    info.lengthBp = lengthBp;
    info.coverage = coverage;
    info.removed = false;
    return dag.addNode(info);
}

bool dagHasOutEdge(
    const dinara::mode3::DirectedAnchorGraph& dag,
    dinara::mode3::DagNodeId from,
    dinara::mode3::DagNodeId to)
{
    const auto& out = dag.getOutEdges(from);
    return std::find(out.begin(), out.end(), to) != out.end();
}

uint64_t findDagSegmentByAnchorChain(
    const dinara::mode3::DirectedAnchorGraph& dag,
    const std::vector<dinara::mode3::DagNodeId>& chain)
{
    for(uint64_t segId : dag.getActiveNodeIds()) {
        if(dag.getNode(segId).anchorChain == chain) {
            return segId;
        }
    }
    return std::numeric_limits<uint64_t>::max();
}

} // namespace

TEST_CASE("DirectedAnchorGraph unitigifyAll collapses MBG-style linear chains and rewrites paths",
    "[integration][dag][unitigify][mbg]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segR = addDagSimpleNode(dag, 90, 7, 5.0);
    const uint64_t segS = addDagSimpleNode(dag, 91, 7, 5.0);
    const uint64_t segP = addDagSimpleNode(dag, 100, 8, 10.0);
    const uint64_t segA = addDagSimpleNode(dag, 101, 10, 20.0);
    const uint64_t segB = addDagSimpleNode(dag, 102, 12, 30.0);
    const uint64_t segC = addDagSimpleNode(dag, 103, 14, 40.0);
    const uint64_t segQ = addDagSimpleNode(dag, 104, 9, 50.0);
    const uint64_t segT = addDagSimpleNode(dag, 110, 7, 5.0);
    const uint64_t segU = addDagSimpleNode(dag, 111, 7, 5.0);

    dag.addEdge(fwdNodeId(segR), fwdNodeId(segP), 1);
    dag.addEdge(fwdNodeId(segS), fwdNodeId(segP), 1);
    dag.addEdge(fwdNodeId(segP), fwdNodeId(segA), 2);
    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 3);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segC), 4);
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segQ), 1);
    dag.addEdge(fwdNodeId(segQ), fwdNodeId(segT), 1);
    dag.addEdge(fwdNodeId(segQ), fwdNodeId(segU), 1);

    const uint64_t p0 = dag.addPath({
        fwdNodeId(segP),
        fwdNodeId(segA),
        fwdNodeId(segB),
        fwdNodeId(segC),
        fwdNodeId(segQ)
    }, 1);
    const uint64_t p1 = dag.addPath({
        fwdNodeId(segB),
        fwdNodeId(segC),
        fwdNodeId(segQ)
    }, 1);
    const uint64_t p2 = dag.addPath({
        fwdNodeId(segP),
        fwdNodeId(segA),
        fwdNodeId(segB)
    }, 1);

    REQUIRE(dag.nodeCount() == 9);
    REQUIRE(dag.pathCount() == 3);

    dag.unitigifyAll();

    REQUIRE(dag.nodeExists(segP) == false);
    REQUIRE(dag.nodeExists(segA) == false);
    REQUIRE(dag.nodeExists(segB) == false);
    REQUIRE(dag.nodeExists(segC) == false);
    REQUIRE(dag.nodeExists(segQ) == false);
    REQUIRE(dag.nodeExists(segR));
    REQUIRE(dag.nodeExists(segS));
    REQUIRE(dag.nodeExists(segT));
    REQUIRE(dag.nodeExists(segU));
    REQUIRE(dag.nodeCount() == 5);
    REQUIRE(dag.pathCount() == 3);

    uint64_t mergedSeg = std::numeric_limits<uint64_t>::max();
    for(uint64_t segId : dag.getActiveNodeIds()) {
        const auto& node = dag.getNode(segId);
        if(node.anchorChain == std::vector<DagNodeId>{100, 101, 102, 103, 104}) {
            mergedSeg = segId;
            break;
        }
    }
    REQUIRE(mergedSeg != std::numeric_limits<uint64_t>::max());

    const auto& mergedInfo = dag.getNode(mergedSeg);
    REQUIRE(mergedInfo.lengthBp == (8 + (10 - 2) + (12 - 3) + (14 - 4) + (9 - 1)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(segR), fwdNodeId(mergedSeg)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(segS), fwdNodeId(mergedSeg)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedSeg), fwdNodeId(segT)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedSeg), fwdNodeId(segU)));

    REQUIRE(dag.getPath(p0) ==
        std::vector<DagNodeId>{fwdNodeId(mergedSeg)});
    REQUIRE(dag.getPath(p1) ==
        std::vector<DagNodeId>{fwdNodeId(mergedSeg)});
    REQUIRE(dag.getPath(p2) ==
        std::vector<DagNodeId>{fwdNodeId(mergedSeg)});
}

TEST_CASE("DirectedAnchorGraph unitigifyAll handles circular unitigs in MBG style",
    "[integration][dag][unitigify][mbg][circular]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segA = addDagSimpleNode(dag, 201, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 202, 11, 11.0);
    const uint64_t segC = addDagSimpleNode(dag, 203, 12, 12.0);

    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 2);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segC), 2);
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segA), 2);

    const uint64_t p0 = dag.addPath({
        fwdNodeId(segA),
        fwdNodeId(segB),
        fwdNodeId(segC),
        fwdNodeId(segA)
    }, 1);

    REQUIRE(dag.nodeCount() == 3);
    dag.unitigifyAll();
    REQUIRE(dag.nodeCount() == 1);
    REQUIRE(dag.pathCount() == 1);

    const auto active = dag.getActiveNodeIds();
    REQUIRE(active.size() == 1);
    const uint64_t mergedSeg = active.front();

    const auto& rewritten = dag.getPath(p0);
    REQUIRE(rewritten.empty() == false);
    for(DagNodeId n : rewritten) {
        REQUIRE(segmentOf(n) == mergedSeg);
    }
}

TEST_CASE("DirectedAnchorGraph unitigifyAll leaves non-unitigifiable branch nodes unchanged",
    "[integration][dag][unitigify][mbg][no-merge]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segA = addDagSimpleNode(dag, 401, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 402, 11, 11.0);
    const uint64_t segC = addDagSimpleNode(dag, 403, 12, 12.0);
    const uint64_t segD = addDagSimpleNode(dag, 404, 13, 13.0);

    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 2);
    dag.addEdge(fwdNodeId(segA), fwdNodeId(segC), 2);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segD), 2);
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segD), 2);

    const uint64_t p0 = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segB), fwdNodeId(segD)
    }, 1);
    const uint64_t p1 = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segC), fwdNodeId(segD)
    }, 1);

    dag.unitigifyAll();

    REQUIRE(dag.nodeCount() == 4);
    REQUIRE(dag.nodeExists(segA));
    REQUIRE(dag.nodeExists(segB));
    REQUIRE(dag.nodeExists(segC));
    REQUIRE(dag.nodeExists(segD));

    REQUIRE(dag.getPath(p0) == std::vector<DagNodeId>{
        fwdNodeId(segA), fwdNodeId(segB), fwdNodeId(segD)
    });
    REQUIRE(dag.getPath(p1) == std::vector<DagNodeId>{
        fwdNodeId(segA), fwdNodeId(segC), fwdNodeId(segD)
    });

    for(uint64_t segId : dag.getActiveNodeIds()) {
        REQUIRE(dag.getNode(segId).anchorChain.size() == 1);
    }
}

TEST_CASE("DirectedAnchorGraph unitigifyAll merges disjoint chains and is idempotent",
    "[integration][dag][unitigify][mbg][idempotent]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segA = addDagSimpleNode(dag, 501, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 502, 12, 12.0);
    const uint64_t segC = addDagSimpleNode(dag, 503, 14, 14.0);
    const uint64_t segD = addDagSimpleNode(dag, 504, 16, 16.0);
    const uint64_t segE = addDagSimpleNode(dag, 505, 9, 9.0);

    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 3);
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segD), 4);

    const uint64_t p0 = dag.addPath({fwdNodeId(segA), fwdNodeId(segB)}, 1);
    const uint64_t p1 = dag.addPath({fwdNodeId(segC), fwdNodeId(segD)}, 1);

    dag.unitigifyAll();

    REQUIRE(dag.nodeCount() == 3);
    REQUIRE(dag.nodeExists(segE));
    REQUIRE_FALSE(dag.nodeExists(segA));
    REQUIRE_FALSE(dag.nodeExists(segB));
    REQUIRE_FALSE(dag.nodeExists(segC));
    REQUIRE_FALSE(dag.nodeExists(segD));

    const uint64_t mergedAB =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{501, 502});
    const uint64_t mergedCD =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{503, 504});
    REQUIRE(mergedAB != std::numeric_limits<uint64_t>::max());
    REQUIRE(mergedCD != std::numeric_limits<uint64_t>::max());

    REQUIRE(dag.getPath(p0) == std::vector<DagNodeId>{fwdNodeId(mergedAB)});
    REQUIRE(dag.getPath(p1) == std::vector<DagNodeId>{fwdNodeId(mergedCD)});

    const uint64_t nodesAfterFirst = dag.nodeCount();
    const uint64_t edgesAfterFirst = dag.edgeCount();
    const auto path0AfterFirst = dag.getPath(p0);
    const auto path1AfterFirst = dag.getPath(p1);

    dag.unitigifyAll();

    REQUIRE(dag.nodeCount() == nodesAfterFirst);
    REQUIRE(dag.edgeCount() == edgesAfterFirst);
    REQUIRE(dag.getPath(p0) == path0AfterFirst);
    REQUIRE(dag.getPath(p1) == path1AfterFirst);
}

TEST_CASE("DirectedAnchorGraph unitigifyAll preserves MBG endpoint rewiring corner cases",
    "[integration][dag][unitigify][mbg][rewire]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segA = addDagSimpleNode(dag, 601, 20, 8.0);
    const uint64_t segB = addDagSimpleNode(dag, 602, 21, 8.0);
    const uint64_t segC = addDagSimpleNode(dag, 603, 22, 8.0);
    const uint64_t segX = addDagSimpleNode(dag, 604, 10, 5.0);
    const uint64_t segY = addDagSimpleNode(dag, 605, 10, 5.0);

    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 3);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segC), 4);
    dag.addEdge(revNodeId(segA), revNodeId(segC), 5);  // rc(first)->rc(last)
    dag.addEdge(revNodeId(segA), fwdNodeId(segA), 2);  // rc(first)->first
    dag.addEdge(revNodeId(segA), fwdNodeId(segX), 1);  // rc(first)->outside
    dag.addEdge(fwdNodeId(segC), revNodeId(segC), 2);  // last->rc(last)
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segY), 1);  // last->outside

    const uint64_t pFwd = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segB), fwdNodeId(segC)
    }, 1);
    const uint64_t pRev = dag.addPath({
        revNodeId(segC), revNodeId(segB), revNodeId(segA)
    }, 1);

    dag.unitigifyAll();

    const uint64_t mergedABC =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{601, 602, 603});
    REQUIRE(mergedABC != std::numeric_limits<uint64_t>::max());

    REQUIRE(dag.getPath(pFwd) == std::vector<DagNodeId>{fwdNodeId(mergedABC)});
    REQUIRE(dag.getPath(pRev) == std::vector<DagNodeId>{revNodeId(mergedABC)});

    REQUIRE(dagHasOutEdge(dag, revNodeId(mergedABC), revNodeId(mergedABC)));
    REQUIRE(dagHasOutEdge(dag, revNodeId(mergedABC), fwdNodeId(mergedABC)));
    REQUIRE(dagHasOutEdge(dag, revNodeId(mergedABC), fwdNodeId(segX)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedABC), revNodeId(mergedABC)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedABC), fwdNodeId(segY)));
}

TEST_CASE("DirectedAnchorGraph unitigifyAll rewrites internal and reverse-entry path orientations",
    "[integration][dag][unitigify][mbg][path-rewrite]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    const uint64_t segA = addDagSimpleNode(dag, 701, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 702, 11, 11.0);
    const uint64_t segC = addDagSimpleNode(dag, 703, 12, 12.0);
    const uint64_t segZ = addDagSimpleNode(dag, 704, 8, 6.0);
    const uint64_t segW = addDagSimpleNode(dag, 705, 8, 6.0);

    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 2);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segC), 2);
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segW), 1);
    dag.addEdge(fwdNodeId(segZ), revNodeId(segC), 1);

    const uint64_t p0 = dag.addPath({fwdNodeId(segB), fwdNodeId(segC)}, 1);
    const uint64_t p1 = dag.addPath({revNodeId(segB), revNodeId(segA)}, 1);
    const uint64_t p2 = dag.addPath({
        fwdNodeId(segZ), revNodeId(segC), revNodeId(segB)
    }, 1);
    const uint64_t p3 = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segB), fwdNodeId(segC)
    }, 1);

    dag.unitigifyAll();

    const uint64_t mergedABC =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{701, 702, 703});
    REQUIRE(mergedABC != std::numeric_limits<uint64_t>::max());

    REQUIRE(dag.getPath(p0) == std::vector<DagNodeId>{fwdNodeId(mergedABC)});
    REQUIRE(dag.getPath(p1) == std::vector<DagNodeId>{revNodeId(mergedABC)});
    REQUIRE(dag.getPath(p2) == std::vector<DagNodeId>{
        fwdNodeId(segZ), revNodeId(mergedABC)
    });
    REQUIRE(dag.getPath(p3) == std::vector<DagNodeId>{fwdNodeId(mergedABC)});

    const auto& crossing = dag.getPathsCrossingNode(mergedABC);
    auto hasPath = [&](uint64_t pathIdx) {
        for(const auto& occ : crossing) {
            if(occ.pathIdx == pathIdx) return true;
        }
        return false;
    };
    REQUIRE(hasPath(p0));
    REQUIRE(hasPath(p1));
    REQUIRE(hasPath(p2));
    REQUIRE(hasPath(p3));
}

// Regression test: batch unitigifyAll must preserve edges between
// adjacent chains separated by a hub. Two chains [A,B] and [C,D]
// connected through hub H (H has degree > 1 in both directions)
// should produce merged nodes M_AB and M_CD with edges M_AB→H and H→M_CD.
TEST_CASE("DirectedAnchorGraph unitigifyAll preserves cross-chain edges",
    "[integration][dag][unitigify][mbg][cross-chain]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    // Hub P has out-degree 2 → can't be part of any chain.
    const uint64_t segP = addDagSimpleNode(dag, 800, 10, 5.0);
    const uint64_t segX = addDagSimpleNode(dag, 806, 8, 4.0);
    // Chain 1: A→B
    const uint64_t segA = addDagSimpleNode(dag, 801, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 802, 12, 12.0);
    // Hub H between chains.
    // H needs in-degree >= 2 (so B can't extend into H)
    // and out-degree >= 2 (so C can't extend backward into H).
    const uint64_t segH = addDagSimpleNode(dag, 808, 10, 6.0);
    const uint64_t segZ = addDagSimpleNode(dag, 809, 8, 4.0);  // H→Z
    const uint64_t segW = addDagSimpleNode(dag, 810, 8, 4.0);  // W→H
    // Chain 2: C→D
    const uint64_t segC = addDagSimpleNode(dag, 803, 14, 14.0);
    const uint64_t segD = addDagSimpleNode(dag, 804, 11, 11.0);
    // Hub Q has in-degree 2.
    const uint64_t segQ = addDagSimpleNode(dag, 805, 10, 5.0);
    const uint64_t segY = addDagSimpleNode(dag, 807, 8, 4.0);

    dag.addEdge(fwdNodeId(segP), fwdNodeId(segA), 2);
    dag.addEdge(fwdNodeId(segP), fwdNodeId(segX), 1);  // P out-degree 2
    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 3);   // chain 1 internal
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segH), 2);   // chain 1 → hub
    dag.addEdge(fwdNodeId(segW), fwdNodeId(segH), 1);   // H in-degree 2
    dag.addEdge(fwdNodeId(segH), fwdNodeId(segC), 2);   // hub → chain 2
    dag.addEdge(fwdNodeId(segH), fwdNodeId(segZ), 1);   // H out-degree 2
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segD), 3);   // chain 2 internal
    dag.addEdge(fwdNodeId(segD), fwdNodeId(segQ), 2);
    dag.addEdge(fwdNodeId(segY), fwdNodeId(segQ), 1);   // Q in-degree 2

    // Path crossing both chains through hub.
    const uint64_t p0 = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segB),
        fwdNodeId(segH),
        fwdNodeId(segC), fwdNodeId(segD)
    }, 1);

    dag.unitigifyAll();

    // Find merged nodes by anchor chain.
    const uint64_t mergedAB =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{801, 802});
    const uint64_t mergedCD =
        findDagSegmentByAnchorChain(dag, std::vector<DagNodeId>{803, 804});
    REQUIRE(mergedAB != std::numeric_limits<uint64_t>::max());
    REQUIRE(mergedCD != std::numeric_limits<uint64_t>::max());

    // Edges through hub must exist.
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedAB), fwdNodeId(segH)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(segH), fwdNodeId(mergedCD)));

    // External edges preserved.
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(segP), fwdNodeId(mergedAB)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedCD), fwdNodeId(segQ)));

    // Path must be rewritten to use both merged nodes.
    REQUIRE(dag.getPath(p0) == std::vector<DagNodeId>{
        fwdNodeId(mergedAB), fwdNodeId(segH), fwdNodeId(mergedCD)
    });
}

// Regression test: batch unitigifyAll must perform secondary merges.
// When chain [A,B] and chain [C,D] are adjacent (B→C) with both
// endpoints having degree 1 after merging, the merged nodes M1→M2
// should themselves be merged into a single node M3.
TEST_CASE("DirectedAnchorGraph unitigifyAll performs secondary merges",
    "[integration][dag][unitigify][mbg][secondary-merge]")
{
    using namespace dinara::mode3;

    DirectedAnchorGraph dag;

    // Hub P with out-degree 2 prevents A from extending backward.
    const uint64_t segP = addDagSimpleNode(dag, 900, 10, 5.0);
    const uint64_t segX = addDagSimpleNode(dag, 906, 8, 4.0);
    // Chain 1: A→B (P has out-degree 2, so A can't extend backward past P)
    const uint64_t segA = addDagSimpleNode(dag, 901, 10, 10.0);
    const uint64_t segB = addDagSimpleNode(dag, 902, 12, 12.0);
    // Chain 2: C→D (Q has in-degree 2, so D can't extend forward past Q)
    const uint64_t segC = addDagSimpleNode(dag, 903, 14, 14.0);
    const uint64_t segD = addDagSimpleNode(dag, 904, 11, 11.0);
    // Hub Q with in-degree 2 prevents D from extending forward.
    const uint64_t segQ = addDagSimpleNode(dag, 905, 10, 5.0);
    const uint64_t segY = addDagSimpleNode(dag, 907, 8, 4.0);

    dag.addEdge(fwdNodeId(segP), fwdNodeId(segA), 2);
    dag.addEdge(fwdNodeId(segP), fwdNodeId(segX), 1);  // P out-degree 2
    dag.addEdge(fwdNodeId(segA), fwdNodeId(segB), 3);
    dag.addEdge(fwdNodeId(segB), fwdNodeId(segC), 2);   // cross-chain
    dag.addEdge(fwdNodeId(segC), fwdNodeId(segD), 3);
    dag.addEdge(fwdNodeId(segD), fwdNodeId(segQ), 2);
    dag.addEdge(fwdNodeId(segY), fwdNodeId(segQ), 1);   // Q in-degree 2

    // Path spanning all four chain members.
    const uint64_t p0 = dag.addPath({
        fwdNodeId(segA), fwdNodeId(segB),
        fwdNodeId(segC), fwdNodeId(segD)
    }, 1);

    dag.unitigifyAll();

    // After batch merge: M_AB and M_CD both have degree 1 connecting
    // to each other. The secondary pass should merge them into M_ABCD.
    const uint64_t mergedABCD =
        findDagSegmentByAnchorChain(dag,
            std::vector<DagNodeId>{901, 902, 903, 904});
    REQUIRE(mergedABCD != std::numeric_limits<uint64_t>::max());

    // The path should reference the fully-merged node.
    REQUIRE(dag.getPath(p0) == std::vector<DagNodeId>{
        fwdNodeId(mergedABCD)
    });

    // External edges must reach the final merged node.
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(segP), fwdNodeId(mergedABCD)));
    REQUIRE(dagHasOutEdge(dag, fwdNodeId(mergedABCD), fwdNodeId(segQ)));
}


// =============================================================================
// Test: .fastq.gz loading produces the same reads as .fastq
// =============================================================================
TEST_CASE("ReadLoader handles gzip-compressed fastq files", "[readloader][gz]") {
    // Create a small fastq in memory, write both plain and gzipped versions.
    auto dir = makeUniqueTempDir("dinara_gz_test_");
    auto plainPath = dir / "reads.fastq";
    auto gzPath = dir / "reads.fastq.gz";

    // Write plain fastq.
    {
        std::ofstream out(plainPath);
        for(int i = 0; i < 100; i++) {
            std::string seq(200, 'A');
            // Make each read slightly different.
            for(int j = 0; j < 200; j++) {
                seq[j] = "ACGT"[(i * 7 + j * 3) % 4];
            }
            out << "@read_" << i << "\n" << seq << "\n+\n" << std::string(200, '~') << "\n";
        }
    }

    // Create gzipped version using system gzip.
    {
        std::string cmd = "gzip -k " + plainPath.string();
        REQUIRE(std::system(cmd.c_str()) == 0);
        REQUIRE(fs::exists(gzPath));
    }

    // Load from plain fastq.
    uint64_t countPlain = 0;
    {
        auto subdir = dir / "plain";
        fs::create_directories(subdir);
        std::string prefix = subdir.string() + "/";
        withSilencedIoInDir(subdir, [&] {
            Assembler assembler(prefix, true, 0, 4096);
            assembler.addReads(plainPath.string(), 0, false, 1);
            countPlain = assembler.getReads().readCount();
        });
    }

    // Load from gzipped fastq.
    uint64_t countGz = 0;
    {
        auto subdir = dir / "gz";
        fs::create_directories(subdir);
        std::string prefix = subdir.string() + "/";
        withSilencedIoInDir(subdir, [&] {
            Assembler assembler(prefix, true, 0, 4096);
            assembler.addReads(gzPath.string(), 0, false, 1);
            countGz = assembler.getReads().readCount();
        });
    }

    REQUIRE(countPlain == 100);
    REQUIRE(countGz == 100);
    REQUIRE(countPlain == countGz);

    fs::remove_all(dir);
}


TEST_CASE("ReadLoader: GIAB fastq.gz matches plain fastq", "[readloader][gz][giab]") {
    // Use the real GIAB test file if available.
    fs::path plainPath = fs::path(__FILE__).parent_path() / "GIAB_HG002_PAW70337_RAW_chr1_15-15.4.fastq";
    if(!fs::exists(plainPath)) {
        WARN("Skipping: GIAB test file not found");
        return;
    }

    auto dir = makeUniqueTempDir("dinara_gz_giab_");

    // Create gzipped version.
    fs::path gzPath = dir / "reads.fastq.gz";
    {
        std::string cmd = "gzip -c " + plainPath.string() + " > " + gzPath.string();
        REQUIRE(std::system(cmd.c_str()) == 0);
        REQUIRE(fs::exists(gzPath));
    }

    uint64_t countPlain = 0;
    {
        auto subdir = dir / "plain";
        fs::create_directories(subdir);
        std::string prefix = subdir.string() + "/";
        withSilencedIoInDir(subdir, [&] {
            Assembler assembler(prefix, true, 0, 4096);
            assembler.addReads(plainPath.string(), 0, false, 1);
            countPlain = assembler.getReads().readCount();
        });
    }

    uint64_t countGz = 0;
    {
        auto subdir = dir / "gz";
        fs::create_directories(subdir);
        std::string prefix = subdir.string() + "/";
        withSilencedIoInDir(subdir, [&] {
            Assembler assembler(prefix, true, 0, 4096);
            assembler.addReads(gzPath.string(), 0, false, 1);
            countGz = assembler.getReads().readCount();
        });
    }

    REQUIRE(countPlain > 0);
    REQUIRE(countPlain == countGz);
    INFO("GIAB: " << countPlain << " reads from both plain and gz");

    fs::remove_all(dir);
}
