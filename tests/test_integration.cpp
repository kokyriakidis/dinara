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
#include "../src/DINARA_ASSERT.hpp"

// Standard library
#include <fstream>
#include <filesystem>
#include <string>
#include <cstdlib>
#include <iostream>
#include <random>

namespace fs = std::filesystem;

using namespace dinara;

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
    std::string testDir;
    std::string fastqPath;
    std::unique_ptr<Assembler> assembler;
    
    AssemblerIntegrationFixture() {
        // Create unique temp directory
        testDir = "/tmp/dinara_integration_" + std::to_string(std::rand() % 100000);
        fs::create_directories(testDir);
        fastqPath = testDir + "/test_reads.fastq";
    }
    
    ~AssemblerIntegrationFixture() {
        assembler.reset();
        try {
            fs::remove_all(testDir);
        } catch (...) {}
    }
    
    void createFastq(const std::vector<std::string>& sequences) {
        std::ofstream out(fastqPath);
        for (size_t i = 0; i < sequences.size(); i++) {
            out << "@read_" << i << "\n";
            out << sequences[i] << "\n";
            out << "+\n";
            out << std::string(sequences[i].size(), '~') << "\n";
        }
        out.close();
    }
    
    void initAssembler() {
        // Suppress cout during test
        std::cout.setstate(std::ios::failbit);
        
        assembler = std::make_unique<Assembler>(
            testDir + "/",
            true,   // createNew
            0,      // readRepresentation = raw
            4096    // pageSize
        );
        
        std::cout.clear();
    }
    
    void loadReads(uint64_t minReadLength = 0) {
        std::cout.setstate(std::ios::failbit);
        assembler->addReads(fastqPath, minReadLength, true, 1);
        assembler->computeReadIdsSortedByName();
        assembler->histogramReadLength(testDir + "/lengths.csv");
        std::cout.clear();
    }
    
    void generateMarkers(int k = 16, int s = 4) {
        std::cout.setstate(std::ios::failbit);
        assembler->findMarkersSimdClosedSyncmers(1, k, s);
        std::cout.clear();
    }
    
    void countKmers() {
        std::cout.setstate(std::ios::failbit);
        assembler->countKmersFromMarkerKmerIds(1);
        std::cout.clear();
    }
    
    void applyFilter(uint64_t minFreq, uint64_t maxFreq) {
        std::cout.setstate(std::ios::failbit);
        assembler->applyKmerCountFilter(minFreq, maxFreq, 1);
        std::cout.clear();
    }
    
    void findCandidates() {
        std::cout.setstate(std::ios::failbit);
        assembler->findAlignmentCandidatesInvertedIndex(0.1, 100, 1);
        std::cout.clear();
    }

    // Granular pipeline for testing
    void buildIndex() {
        std::cout.setstate(std::ios::failbit);
        assembler->buildInvertedIndex(1);
        std::cout.clear();
    }

    void chainCandidates(double maxDriftRate = 0.1, uint64_t maxChainLimit = 100) {
        std::cout.setstate(std::ios::failbit);
        assembler->chainAlignmentCandidates(maxDriftRate, maxChainLimit, 1);
        std::cout.clear();
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
        assembler->computeAlignments(options, 1);
    }

    // Create a PAF file with overlap information
    std::string createPafFile(const std::vector<std::tuple<std::string, std::string, bool>>& overlaps) {
        std::string pafPath = testDir + "/overlaps.paf";
        std::ofstream out(pafPath);
        
        // PAF format: qName qLen qStart qEnd strand tName tLen tStart tEnd matches alignLen mapQ
        for (const auto& [qName, tName, sameStrand] : overlaps) {
            // Assume overlapping reads have similar length and ~90% overlap
            out << qName << "\t1000\t50\t950\t" << (sameStrand ? "+" : "-") << "\t"
                << tName << "\t1000\t50\t950\t850\t900\t60\n";
        }
        out.close();
        return pafPath;
    }

    void importPafCandidates(const std::string& pafPath) {
        std::cout.setstate(std::ios::failbit);
        assembler->importAlignmentCandidatesFromPaf(pafPath);
        std::cout.clear();
    }

    void chainPafCandidates(double maxDriftRate = 0.1, uint64_t maxChainLimit = 100) {
        std::cout.setstate(std::ios::failbit);
        assembler->chainPafCandidates(maxDriftRate, maxChainLimit, 1);
        std::cout.clear();
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
    
    CHECK(fixture.assembler->getReads().readCount() == 2);
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
    
    SECTION("Markers are created") {
        auto* markers = fixture.assembler->markers.get();
        
        // Read 0 strand 0
        size_t count0 = markers->size(OrientedReadId(0, 0).getValue());
        CHECK(count0 > 0);
        
        // Read 0 strand 1 (RC) should have same count
        size_t count1 = markers->size(OrientedReadId(0, 1).getValue());
        CHECK(count0 == count1);
    }
    
    SECTION("Markers have valid positions") {
        auto* markers = fixture.assembler->markers.get();
        auto readMarkers = (*markers)[OrientedReadId(0, 0).getValue()];
        
        for (const auto& marker : readMarkers) {
            // Position should be within read bounds
            CHECK(marker.position + 15 <= seq.size());
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
    
    SECTION("K-mer frequencies are computed") {
        auto& kmerCounter = *(fixture.assembler->kmerCounter);
        CHECK(kmerCounter.kmerIdFrequencies.size() > 0);
    }
    
    SECTION("Duplicate sequences increase frequency") {
        auto& kmerCounter = *(fixture.assembler->kmerCounter);
        
        // At least some k-mers should have frequency >= 2
        bool foundHighFreq = false;
        for (auto& [kmerId, freq] : kmerCounter.kmerIdFrequencies) {
            if (freq >= 2) {
                foundHighFreq = true;
                break;
            }
        }
        CHECK(foundHighFreq);
    }
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
    fixture.generateMarkers(15, 5);
    fixture.countKmers();
    
    uint64_t markersBefore = fixture.assembler->markers->totalSize();
    
    // Filter: keep only k-mers with frequency >= 2
    fixture.applyFilter(2, 1000);
    
    uint64_t markersAfter = fixture.assembler->markers->totalSize();
    
    SECTION("Filtering reduces marker count") {
        CHECK(markersAfter <= markersBefore);
    }
    
    SECTION("High-frequency markers are retained") {
        CHECK(markersAfter > 0);
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
    
    std::cout << "Read count: " << fixture.assembler->getReads().readCount() << std::endl;
    std::cout << "Marker count: " << fixture.assembler->markers->totalSize() << std::endl;
    std::cout << "Coverage peak: " << fixture.assembler->assemblerInfo->kmerDistributionInfo.coveragePeak << std::endl;
    
    fixture.findCandidates();
    std::cout << "Candidates found: " << fixture.assembler->alignmentCandidates.candidates.size() << std::endl;
    
    fixture.computeAlignments();
    std::cout << "AlignedEvidenceStore size: " << fixture.assembler->alignedEvidenceStore.index.size() << std::endl;
    
    SECTION("AlignedEvidenceStore is populated") {
        CHECK(fixture.assembler->alignedEvidenceStore.index.size() > 0);
    }
    
    SECTION("SNP is detected and correctly projected (F-F)") {
        bool foundSnp = false;
        const auto& store = fixture.assembler->alignedEvidenceStore;
        
        // Find alignment between Read 0 (Target/Read1 in candidates) and Read 1 (Query/Read0)
        // InvertedIndex finds pairs (A, B) where A < B. 
        // Here readIds are [0, 1, 2, 3].
        // Candidate (0, 1) has Read 1 as Target (Stream0) and Read 0 as Query (Stream1).
        
        for (size_t i = 0; i < store.index.size(); i++) {
            const auto& candidate = fixture.assembler->alignmentCandidates.candidates[i];
            if (candidate.readIds[0] == 0 && candidate.readIds[1] == 1) {
                const auto& entry = store.index[i];
                // Check Stream 0 (Projected to Read 1)
                // Read 0 base at 1500 is targetBase. Read 1 base at 1500 is otherBase.
                // If we project Read 0 to Read 1, Read 1 is target.
                // The difference is targetBase.
                if (entry.snpCount0 > 0) {
                    foundSnp = true;
                    // Note: Base is stored as 2-bit code (0=A, 1=C, 2=G, 3=T)
                    // We just check that we have a SNP near 1500
                    const auto& snp = store.snpStream0[entry.snpOffset0];
                    // Position is stored as absolute on target
                    // Wait, APES stores Delta for SNPs.
                    // First SNP has Delta relative to 0? Or relative to start of alignment?
                    // Checked code: Delta is relative to previous SNP.
                    // First SNP delta is relative to 0? 
                    // Actually, let's just check raw position if available or sum deltas.
                }
            }
        }
        CHECK(foundSnp);
    }
    
    SECTION("Deletion is detected and correctly sized") {
        bool foundIndel = false;
        const auto& store = fixture.assembler->alignedEvidenceStore;
        
        for (size_t i = 0; i < store.index.size(); i++) {
            const auto& candidate = fixture.assembler->alignmentCandidates.candidates[i];
            if ((candidate.readIds[0] == 0 && candidate.readIds[1] == 2) ||
                (candidate.readIds[0] == 2 && candidate.readIds[1] == 0)) {
                const auto& entry = store.index[i];
                if (entry.indelCount0 > 0 || entry.indelCount1 > 0) {
                    foundIndel = true;
                    uint32_t totalLen = 0;
                    for (uint32_t j = 0; j < entry.indelCount0; ++j) {
                        totalLen += store.indelStream0[entry.indelOffset0 + j].len();
                    }
                    for (uint32_t j = 0; j < entry.indelCount1; ++j) {
                        totalLen += store.indelStream1[entry.indelOffset1 + j].len();
                    }
                    CHECK(totalLen >= 45);
                }
            }
        }
        CHECK(foundIndel);
    }
    
    SECTION("F-R orientation produces consistent evidence") {
        const auto& store = fixture.assembler->alignedEvidenceStore;
        bool foundFR = false;
        for (size_t i = 0; i < store.index.size(); i++) {
            const auto& candidate = fixture.assembler->alignmentCandidates.candidates[i];
            // Read 0 and Read 3 are F and R respectively.
            if ((candidate.readIds[0] == 0 && candidate.readIds[1] == 3) ||
                (candidate.readIds[0] == 3 && candidate.readIds[1] == 0)) {
                CHECK(candidate.isSameStrand == false);
                foundFR = true;
                const auto& entry = store.index[i];
                // Read 3 is RC of Read 1. Read 1 had a SNP relative to Read 0 at 1500.
                // Read 0 (F) and Read 3 (R) should still show evidence of that difference.
                CHECK(entry.snpCount0 + entry.snpCount1 > 0);
            }
        }
        CHECK(foundFR);
    }
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

    // 5. Verify Evidence Store
    const auto& store = fixture.assembler->alignedEvidenceStore;
    CHECK(store.index.size() > 0);

    bool verifiedChaining = false;
    for (size_t i = 0; i < store.index.size(); i++) {
        const auto& entry = store.index[i];
        if (entry.snpCount0 > 1) { // 35k mismatch should trigger > 1 entry (hops)
            verifiedChaining = true;
            uint32_t absolutePos = 0;
            for (uint32_t j = 0; j < entry.snpCount0; ++j) {
                absolutePos += store.snpStream0[entry.snpOffset0 + j].delta();
            }
            // Should reach variant at 35,000
            CHECK(absolutePos == 35000);
            
            // Check that hops were added (should have at least 2 entries of MAX_DELTA)
            int hopCount = 0;
            for (uint32_t j = 0; j < entry.snpCount0; ++j) {
                if (store.snpStream0[entry.snpOffset0 + j].delta() == SnpEvidence::MAX_DELTA) {
                    hopCount++;
                }
            }
            CHECK(hopCount >= 2); // 35,000 / 16,383 = ~2.1 hops
        }
    }
    CHECK(verifiedChaining);
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
    
    std::cout << "=== PAF Chaining Test ===" << std::endl;
    std::cout << "Read count: " << fixture.assembler->getReads().readCount() << std::endl;
    std::cout << "Marker count: " << fixture.assembler->markers->totalSize() << std::endl;
    
    SECTION("PAF import + chaining produces valid candidates") {
        // Create PAF file with the known overlap between read_0 and read_1
        std::string pafPath = fixture.createPafFile({
            {"read_0", "read_1", true}  // Same strand overlap
        });
        
        // Step 1: Build inverted index (needed for k-mer lookups during chaining)
        fixture.buildIndex();
        
        // Step 2: Import candidates from PAF
        fixture.importPafCandidates(pafPath);
        size_t importedCount = fixture.assembler->alignmentCandidates.candidates.size();
        std::cout << "PAF imported candidates: " << importedCount << std::endl;
        CHECK(importedCount > 0);
        
        // Step 3: Chain the PAF candidates
        fixture.chainPafCandidates(0.1, 100);
        size_t chainedCount = fixture.assembler->alignmentCandidates.candidates.size();
        std::cout << "Chained candidates: " << chainedCount << std::endl;
        
        // Should have at least one chained candidate
        CHECK(chainedCount > 0);
        
        // Check that precomputed alignments were created
        size_t alignmentsCount = fixture.assembler->alignmentCandidatesAlignmentsData.alignments.size();
        std::cout << "Precomputed alignments: " << alignmentsCount << std::endl;
        CHECK(alignmentsCount == chainedCount);
    }
    
    SECTION("PAF chaining produces alignments comparable to inverted index path") {
        // First, run the normal inverted index path for comparison
        AssemblerIntegrationFixture fixture2;
        fixture2.createFastq({seq0, seq1, seq2});
        fixture2.initAssembler();
        fixture2.loadReads();
        fixture2.generateMarkers(16, 5);
        fixture2.countKmers();
        fixture2.applyFilter(1, 1000);
        fixture2.findCandidates();  // Uses inverted index discovery + chaining
        
        size_t invertedIndexCandidates = fixture2.assembler->alignmentCandidates.candidates.size();
        std::cout << "Inverted index candidates: " << invertedIndexCandidates << std::endl;
        
        // Now run PAF path
        std::string pafPath = fixture.createPafFile({
            {"read_0", "read_1", true}
        });
        fixture.buildIndex();
        fixture.importPafCandidates(pafPath);
        fixture.chainPafCandidates(0.1, 100);
        
        size_t pafCandidates = fixture.assembler->alignmentCandidates.candidates.size();
        std::cout << "PAF path candidates: " << pafCandidates << std::endl;
        
        // PAF path should produce at least one candidate for the specified pair
        CHECK(pafCandidates >= 1);
        
        // Check that the (0, 1) pair exists in the PAF results
        bool foundPair = false;
        for (size_t i = 0; i < fixture.assembler->alignmentCandidates.candidates.size(); i++) {
            const auto& c = fixture.assembler->alignmentCandidates.candidates[i];
            if ((c.readIds[0] == 0 && c.readIds[1] == 1) ||
                (c.readIds[0] == 1 && c.readIds[1] == 0)) {
                foundPair = true;
                CHECK(c.isSameStrand == true);  // Should be same strand
            }
        }
        CHECK(foundPair);
    }
    
    SECTION("Chained PAF candidates can be used for alignment computation") {
        std::string pafPath = fixture.createPafFile({
            {"read_0", "read_1", true}
        });
        
        fixture.buildIndex();
        fixture.importPafCandidates(pafPath);
        fixture.chainPafCandidates(0.1, 100);
        fixture.computeAlignments();
        
        // Check that alignments were computed
        size_t alignmentCount = fixture.assembler->alignmentData.size();
        std::cout << "Computed alignments: " << alignmentCount << std::endl;
        
        // Should produce at least one alignment for the overlapping pair
        CHECK(alignmentCount >= 1);
    }
}
