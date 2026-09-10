// Tests for CIGAR-driven het-site detection (AssemblerCigarSnpSites.cpp).
//
// Every case here corresponds to a bug that actually reached a run: a site
// counted once per covering read instead of once, an arm holding one read on
// both strands (which aborts the external-anchor export), a mismatch column
// whose bases do not actually differ. They are unit tests rather than fixture
// runs because none of that needs real data -- the detector reads only the
// CIGAR store and the read sequences, and both can be built by hand.
//
// Note what is NOT covered, and why: duplicate-locus collapse and the
// "member already sits in an anchor" filter live inside main.cpp's assemble(),
// not behind any callable interface, so they cannot be reached from a test at
// all. That is a structural gap, not an oversight.

#include "../external/catch2/catch.hpp"

#include "../src/Assembler.hpp"
#include "../src/HifiasmImportedCigarStore.hpp"
#include "../src/OverlapCigarStore.hpp"
#include "../src/Reads.hpp"

#include <filesystem>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace dinara;
namespace fs = std::filesystem;

namespace {

fs::path makeSnpTestDir()
{
    static uint64_t counter = 0;
    const fs::path dir = fs::temp_directory_path() /
        ("dinara_snp_test_" + std::to_string(::getpid()) + "_" + std::to_string(counter++));
    fs::create_directories(dir);
    return dir;
}

// A set of reads that are identical except at one base, plus an all-against-all
// set of full-length same-strand overlaps whose CIGARs say exactly where they
// differ. This is the smallest input that exercises the real detection path.
class SnpScenario {
public:
    fs::path dir;
    std::unique_ptr<Assembler> assembler;
    std::vector<std::string> sequences;

    explicit SnpScenario(const std::vector<std::string>& seqs) : sequences(seqs)
    {
        dir = makeSnpTestDir();
        const fs::path fastq = dir / "reads.fastq";
        {
            std::ofstream out(fastq);
            for(size_t i = 0; i < seqs.size(); i++) {
                out << "@read_" << i << "\n" << seqs[i] << "\n+\n"
                    << std::string(seqs[i].size(), '~') << "\n";
            }
        }
        const auto cwd = fs::current_path();
        fs::current_path(dir);
        assembler = std::make_unique<Assembler>(dir.string() + "/", true, 0, 4096);
        assembler->addReads(fastq.string(), 0, true, 1);
        assembler->computeReadIdsSortedByName();
        fs::current_path(cwd);
    }

    ~SnpScenario()
    {
        assembler.reset();
        try { fs::remove_all(dir); } catch(...) {}
    }

    // Add a full-length, same-strand overlap between two reads whose sequences
    // have equal length. The CIGAR is derived from the sequences themselves, so
    // the tokens cannot disagree with the bases the detector will read back --
    // which is the invariant the detector's own self-checks assert.
    void addOverlap(uint64_t pairKey, ReadId q, ReadId t)
    {
        const std::string& a = sequences[q];
        const std::string& b = sequences[t];
        REQUIRE(a.size() == b.size());

        std::vector<uint16_t> tokens;
        uint8_t runOp = 0xff;
        uint16_t runLen = 0;
        for(size_t i = 0; i < a.size(); i++) {
            const uint8_t op = (a[i] == b[i]) ? uint8_t(CigarOpMatch) : uint8_t(CigarOpMismatch);
            if(op == runOp) {
                ++runLen;
            } else {
                if(runOp != 0xff) tokens.push_back(CigarToken(runOp, runLen).data);
                runOp = op;
                runLen = 1;
            }
        }
        if(runOp != 0xff) tokens.push_back(CigarToken(runOp, runLen).data);

        assembler->hifiasmImportedCigarStore.add(
            pairKey, /*isSameStrand*/ true,
            span<const uint16_t>(tokens.data(), tokens.size()),
            uint32_t(q), uint32_t(t),
            0, uint32_t(a.size()), 0, uint32_t(b.size()));
    }

    void addAllPairs()
    {
        uint64_t key = 0;
        for(ReadId q = 0; q < ReadId(sequences.size()); q++) {
            for(ReadId t = q + 1; t < ReadId(sequences.size()); t++) {
                addOverlap(key++, q, t);
            }
        }
    }

    // Run detection with every statistical filter made permissive, so a failure
    // points at the detection mechanics rather than at threshold tuning. The
    // filters have their own tuning evidence; these tests are about whether the
    // machinery finds the site once and assigns the right reads to each allele.
    std::vector<Assembler::CigarSnpSite> detect(uint64_t alleleFloor = 2)
    {
        std::vector<Assembler::CigarSnpSite> sites;
        const auto cwd = fs::current_path();
        fs::current_path(dir);
        // The detector reports its funnel to cout, which is useful in a run and
        // noise in a test suite. Capture it rather than silencing it globally,
        // so a failure can still be diagnosed by printing `captured`.
        std::ostringstream captured;
        std::streambuf* const savedCout = std::cout.rdbuf(captured.rdbuf());
        assembler->detectCigarSnpSites(
            /*minDisagreeCount*/ 2,
            /*minDisagreeFraction*/ 0.0,
            /*maxDisagreeFraction*/ 1.0,
            /*hetErrorRate*/ 0.025,
            /*strandBiasPValue*/ 0.0,     // 0 disables: p < 0 is never true
            /*siteMinPurity*/ 0.0,
            /*siteMinAltDominance*/ 0.0,
            /*minAlleleFraction*/ 0.0,
            /*filterHomopolymer*/ false,
            /*filterStr*/ false,
            /*minHomopolymerRun*/ 5,
            /*alleleCoverageRate*/ 0.0,
            /*alleleCoverageFloor*/ alleleFloor,
            /*ploidy*/ 2,
            &sites,
            /*threadCount*/ 1);
        std::cout.rdbuf(savedCout);
        fs::current_path(cwd);
        return sites;
    }
};

// A 120-base backbone with no homopolymer or short-repeat structure around the
// variable site, so context filters have nothing to catch even when enabled.
std::string backbone()
{
    static const std::string unit = "ACGTAGCTTGACCATGGACTGCATCGTAAGCTGACTTGCA";
    return unit + unit + unit;   // 120 bases
}

std::vector<std::string> twoAlleleReads(size_t nA, size_t nC, size_t snpPos)
{
    std::vector<std::string> seqs;
    for(size_t i = 0; i < nA + nC; i++) {
        std::string s = backbone();
        s[snpPos] = (i < nA) ? 'A' : 'C';
        seqs.push_back(s);
    }
    return seqs;
}

} // namespace



TEST_CASE("CIGAR SNP detection finds a planted het site once", "[snp][cigar]")
{
    // Six reads, three carrying each allele, all overlapping each other. The
    // site is visible from every read, so a detector without deduplication
    // would report it six times.
    const size_t snpPos = 61;
    SnpScenario s(twoAlleleReads(3, 3, snpPos));
    s.addAllPairs();

    const auto sites = s.detect();

    REQUIRE(sites.size() == 1);
    REQUIRE(sites[0].alleles.size() == 2);

    // Each arm holds exactly the reads carrying that allele, and every member
    // is recorded at the SNP position in its own read.
    std::vector<size_t> armSizes{sites[0].alleles[0].size(), sites[0].alleles[1].size()};
    std::sort(armSizes.begin(), armSizes.end());
    CHECK(armSizes[0] == 3);
    CHECK(armSizes[1] == 3);

    for(const auto& arm: sites[0].alleles) {
        for(const auto& [orientedReadId, position]: arm) {
            CHECK(position == snpPos);
        }
    }
}


TEST_CASE("CIGAR SNP detection assigns each read to the arm it actually carries",
          "[snp][cigar]")
{
    const size_t snpPos = 61;
    SnpScenario s(twoAlleleReads(2, 4, snpPos));
    s.addAllPairs();

    const auto sites = s.detect();
    REQUIRE(sites.size() == 1);
    REQUIRE(sites[0].alleles.size() == 2);

    // Reads 0-1 carry 'A', reads 2-5 carry 'C'. Whichever arm is which, the
    // partition must follow the sequence, not the arrival order.
    for(const auto& arm: sites[0].alleles) {
        std::set<char> basesInArm;
        for(const auto& [orientedReadId, position]: arm) {
            basesInArm.insert(s.sequences[orientedReadId.getReadId()][position]);
        }
        CHECK(basesInArm.size() == 1);   // an arm is one allele, never mixed
    }

    std::vector<size_t> armSizes{sites[0].alleles[0].size(), sites[0].alleles[1].size()};
    std::sort(armSizes.begin(), armSizes.end());
    CHECK(armSizes[0] == 2);
    CHECK(armSizes[1] == 4);
}


TEST_CASE("CIGAR SNP detection reports no site when all reads agree", "[snp][cigar]")
{
    // Same construction, one allele. Every CIGAR is a single match run, so no
    // position ever accumulates a disagreeing partner.
    SnpScenario s(twoAlleleReads(6, 0, 61));
    s.addAllPairs();
    CHECK(s.detect().empty());
}


TEST_CASE("CIGAR SNP detection needs two disagreeing partners, not one",
          "[snp][cigar]")
{
    // Five reads carry 'A' and one carries 'C'. The lone read disagrees with
    // everyone, but no read sees TWO partners disagreeing with it at that
    // position except the singleton itself -- and the singleton's own arm has
    // one member, below the arm floor. This is the boundary the candidate rule
    // (>= 2 disagreeing) is meant to sit on.
    SnpScenario s(twoAlleleReads(5, 1, 61));
    s.addAllPairs();

    const auto sites = s.detect(/*alleleFloor*/ 2);
    for(const auto& site: sites) {
        for(const auto& arm: site.alleles) {
            CHECK(arm.size() >= 2);   // never emit a single-read allele
        }
    }
}


TEST_CASE("CIGAR SNP detection finds two independent sites separately",
          "[snp][cigar]")
{
    // Two variable positions far apart, with a different read partition at
    // each, so neither can be explained by the other.
    const size_t posA = 41, posB = 91;
    std::vector<std::string> seqs;
    for(size_t i = 0; i < 6; i++) {
        std::string s = backbone();
        s[posA] = (i < 3) ? 'A' : 'C';
        s[posB] = (i % 2 == 0) ? 'G' : 'T';
        seqs.push_back(s);
    }
    SnpScenario s(seqs);
    s.addAllPairs();

    const auto sites = s.detect();
    REQUIRE(sites.size() == 2);

    std::set<uint32_t> positions;
    for(const auto& site: sites) {
        REQUIRE(site.alleles.size() == 2);
        REQUIRE_FALSE(site.alleles[0].empty());
        positions.insert(site.alleles[0].front().second);
    }
    CHECK(positions == std::set<uint32_t>{uint32_t(posA), uint32_t(posB)});
}


TEST_CASE("CIGAR SNP arms never hold one ReadId on both strands", "[snp][cigar]")
{
    // Shasta2Anchors rejects an anchor containing the same ReadId twice, and an
    // arm is exactly what becomes an anchor. Arms are collected by walking every
    // overlap of the owning read with no per-read bookkeeping, so a partner
    // reachable through two overlaps can be appended twice -- which is how this
    // aborted a real export.
    const size_t snpPos = 61;
    SnpScenario s(twoAlleleReads(3, 3, snpPos));
    s.addAllPairs();
    // The same pair again under a different key: a second overlap between reads
    // that already overlap.
    s.addOverlap(/*pairKey*/ 1000, 0, 3);
    s.addOverlap(/*pairKey*/ 1001, 1, 4);

    const auto sites = s.detect();
    for(const auto& site: sites) {
        for(const auto& arm: site.alleles) {
            std::set<ReadId> seen;
            for(const auto& [orientedReadId, position]: arm) {
                const ReadId r = orientedReadId.getReadId();
                CHECK(seen.insert(r).second);   // each ReadId at most once
            }
        }
    }
}
