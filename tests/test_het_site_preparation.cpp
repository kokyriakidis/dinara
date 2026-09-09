// Tests for the two steps between detection and anchor creation
// (src/HetSitePreparation.cpp).
//
// These were inline in main.cpp's assemble() until they were extracted, which
// is why they had no coverage despite deciding which sites survive and which
// reads join an arm. Both are pure functions, so the tests build their input
// directly and never touch an Assembler.

#include "../external/catch2/catch.hpp"

#include "../src/HetSitePreparation.hpp"

#include <set>
#include <vector>

using namespace dinara;

namespace {

using Member = std::pair<OrientedReadId, uint32_t>;

OrientedReadId fwd(uint32_t readId) { return OrientedReadId(readId, 0); }
OrientedReadId rev(uint32_t readId) { return OrientedReadId(readId, 1); }

Assembler::CigarSnpSite makeSite(
    const std::vector<Member>& armA,
    const std::vector<Member>& armB)
{
    Assembler::CigarSnpSite site;
    site.alleles.push_back(vector<Member>(armA.begin(), armA.end()));
    site.alleles.push_back(vector<Member>(armB.begin(), armB.end()));
    return site;
}

uint64_t totalMembers(const vector<Assembler::CigarSnpSite>& sites)
{
    uint64_t n = 0;
    for(const auto& s: sites) for(const auto& a: s.alleles) n += a.size();
    return n;
}

// The set of member occurrences a site covers, which is what identifies a locus.
std::set<std::pair<uint64_t, uint32_t>> occurrencesOf(const Assembler::CigarSnpSite& s)
{
    std::set<std::pair<uint64_t, uint32_t>> out;
    for(const auto& a: s.alleles) {
        for(const auto& m: a) out.insert({m.first.getValue(), m.second});
    }
    return out;
}

} // namespace



TEST_CASE("collapseDuplicateLoci keeps distinct loci", "[hetprep]")
{
    // Two sites sharing no member occurrence are different loci, whatever
    // else they have in common -- here the same reads at different positions.
    vector<Assembler::CigarSnpSite> sites;
    sites.push_back(makeSite({{fwd(0), 100}, {fwd(1), 100}},
                             {{fwd(2), 100}, {fwd(3), 100}}));
    sites.push_back(makeSite({{fwd(0), 500}, {fwd(1), 500}},
                             {{fwd(2), 500}, {fwd(3), 500}}));

    CHECK(collapseDuplicateLoci(sites) == 0);
    CHECK(sites.size() == 2);
}


TEST_CASE("collapseDuplicateLoci drops a second detection of one locus", "[hetprep]")
{
    // Same locus seen by two owners with slightly different covering sets: the
    // second has one member fewer. They share occurrences, so they are the same
    // locus, and the better-supported one must be the survivor.
    vector<Assembler::CigarSnpSite> sites;
    sites.push_back(makeSite({{fwd(0), 100}, {fwd(1), 100}},
                             {{fwd(2), 100}, {fwd(3), 100}}));      // 4 members
    sites.push_back(makeSite({{fwd(0), 100}, {fwd(1), 100}},
                             {{fwd(2), 100}}));                      // 3 members

    CHECK(collapseDuplicateLoci(sites) == 1);
    REQUIRE(sites.size() == 1);
    CHECK(totalMembers(sites) == 4);      // the better-supported one survived
}


TEST_CASE("collapseDuplicateLoci keeps the better-supported site regardless of order",
          "[hetprep]")
{
    // Detection output order is thread-scheduling dependent, so the result must
    // not depend on which duplicate was appended first.
    const auto big   = makeSite({{fwd(0), 7}, {fwd(1), 7}, {fwd(4), 7}},
                                {{fwd(2), 7}, {fwd(3), 7}});          // 5
    const auto small = makeSite({{fwd(0), 7}, {fwd(1), 7}},
                                {{fwd(2), 7}});                       // 3

    vector<Assembler::CigarSnpSite> bigFirst{big, small};
    vector<Assembler::CigarSnpSite> smallFirst{small, big};

    CHECK(collapseDuplicateLoci(bigFirst) == 1);
    CHECK(collapseDuplicateLoci(smallFirst) == 1);
    REQUIRE(bigFirst.size() == 1);
    REQUIRE(smallFirst.size() == 1);
    CHECK(totalMembers(bigFirst) == 5);
    CHECK(totalMembers(smallFirst) == 5);
    CHECK(occurrencesOf(bigFirst[0]) == occurrencesOf(smallFirst[0]));
}


TEST_CASE("collapseDuplicateLoci output order is deterministic", "[hetprep]")
{
    // Survivors are returned sorted by their smallest member occurrence, not in
    // arrival order, because the index becomes the het anchor id and anchor ids
    // are the last tie-break in the journey rebuild. Two runs that disagreed
    // here produced different graphs from identical input.
    const auto atFive = makeSite({{fwd(9), 5}, {fwd(8), 5}}, {{fwd(7), 5}});
    const auto atOne  = makeSite({{fwd(3), 1}, {fwd(2), 1}}, {{fwd(1), 1}});

    vector<Assembler::CigarSnpSite> a{atFive, atOne};
    vector<Assembler::CigarSnpSite> b{atOne, atFive};
    collapseDuplicateLoci(a);
    collapseDuplicateLoci(b);

    REQUIRE(a.size() == 2);
    REQUIRE(b.size() == 2);
    for(size_t i = 0; i < a.size(); i++) {
        CHECK(occurrencesOf(a[i]) == occurrencesOf(b[i]));
    }
}


TEST_CASE("collapseDuplicateLoci treats the two strands of a read as distinct",
          "[hetprep]")
{
    // An occurrence is (OrientedReadId, position). The same ReadId on opposite
    // strands is a different occurrence and must not merge two real loci.
    vector<Assembler::CigarSnpSite> sites;
    sites.push_back(makeSite({{fwd(0), 42}, {fwd(1), 42}}, {{fwd(2), 42}}));
    sites.push_back(makeSite({{rev(0), 42}, {rev(1), 42}}, {{rev(2), 42}}));

    CHECK(collapseDuplicateLoci(sites) == 0);
    CHECK(sites.size() == 2);
}


TEST_CASE("dropAlreadyAnchoredArmMembers removes only occupied positions",
          "[hetprep]")
{
    // Site positions and stored anchor positions are the same frame: a het
    // anchor is a zero-length marker at the SNP base itself. An earlier version
    // offset one side by hetKHalf; the mismatch matched nothing and removed 192
    // of the wrong members while leaving every real collision in place, so the
    // frames being identical is the property under test.
    std::unordered_map<uint64_t, Shasta2AnchorId> occupied;
    const auto occupy = [&](OrientedReadId r, uint32_t position) {
        occupied.emplace((uint64_t(r.getValue()) << 32) | uint64_t(position), 0);
    };
    occupy(fwd(1), 100);      // read 1 at raw 100 is already anchored
    occupy(fwd(3), 100);

    vector<Assembler::CigarSnpSite> sites;
    sites.push_back(makeSite({{fwd(0), 100}, {fwd(1), 100}},
                             {{fwd(2), 100}, {fwd(3), 100}}));

    CHECK(dropAlreadyAnchoredArmMembers(sites, occupied) == 2);
    REQUIRE(sites.size() == 1);
    REQUIRE(sites[0].alleles.size() == 2);
    CHECK(sites[0].alleles[0].size() == 1);
    CHECK(sites[0].alleles[1].size() == 1);
    CHECK(sites[0].alleles[0].front().first == fwd(0));
    CHECK(sites[0].alleles[1].front().first == fwd(2));
}


TEST_CASE("dropAlreadyAnchoredArmMembers ignores a different position on the same read",
          "[hetprep]")
{
    // Occupancy is per (read, position), not per read: a read anchored
    // elsewhere is still a usable member here.
    std::unordered_map<uint64_t, Shasta2AnchorId> occupied;
    occupied.emplace((uint64_t(fwd(1).getValue()) << 32) | uint64_t(999), 0);

    vector<Assembler::CigarSnpSite> sites;
    sites.push_back(makeSite({{fwd(0), 100}, {fwd(1), 100}}, {{fwd(2), 100}}));

    CHECK(dropAlreadyAnchoredArmMembers(sites, occupied) == 0);
    CHECK(totalMembers(sites) == 3);
}


TEST_CASE("dropAlreadyAnchoredArmMembers can empty an arm", "[hetprep]")
{
    // Every member of one arm already anchored leaves that arm empty. The
    // function does not decide what that means -- the arm floor downstream
    // does -- but it must not silently keep the members.
    std::unordered_map<uint64_t, Shasta2AnchorId> occupied;
    for(uint32_t r: {2u, 3u}) {
        occupied.emplace((uint64_t(fwd(r).getValue()) << 32) | uint64_t(100), 0);
    }

    vector<Assembler::CigarSnpSite> sites;
    sites.push_back(makeSite({{fwd(0), 100}, {fwd(1), 100}},
                             {{fwd(2), 100}, {fwd(3), 100}}));

    CHECK(dropAlreadyAnchoredArmMembers(sites, occupied) == 2);
    REQUIRE(sites.size() == 1);
    CHECK(sites[0].alleles[0].size() == 2);
    CHECK(sites[0].alleles[1].empty());
}
