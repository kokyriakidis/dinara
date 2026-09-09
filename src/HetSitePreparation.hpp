#ifndef DINARA_HET_SITE_PREPARATION_HPP
#define DINARA_HET_SITE_PREPARATION_HPP

// The two steps between "detectCigarSnpSites produced sites" and "those sites
// become anchors".
//
// They lived inline in main.cpp's assemble(), which put real decisions -- which
// sites survive, which reads join an arm -- behind no callable interface and so
// out of reach of any test. Both are pure functions of their inputs, so they
// belong here.
//
//   collapseDuplicateLoci   the same locus detected twice, because ownership
//                           elected different owners from divergent covering
//                           sets. Detected by shared member OCCURRENCE.
//
//   dropAlreadyAnchoredArmMembers
//                           a member whose (read, position) an existing anchor
//                           already holds needs no het anchor: that anchor
//                           necessarily isolates its allele.
//
// Both are declared on plain vectors rather than on Assembler so a test can
// call them with hand-built input.

#include "Assembler.hpp"
#include "Shasta2Anchors.hpp"
#include "cstdint.hpp"
#include "vector.hpp"

#include <unordered_map>

namespace dinara {

// Drop every site that is a second detection of a locus another site already
// claims, keeping the better-supported one.
//
// Two sites are the same locus exactly when they share a member OCCURRENCE --
// the same (OrientedReadId, position) -- because one read cannot carry two
// different alleles at one base. No coordinate mapping is needed.
//
// The keeper is the site with the most members over both arms, ties broken on
// the smallest (OrientedReadId, position) so the result does not depend on the
// order detection happened to produce (which is thread-scheduling dependent).
// Survivors are returned in that same deterministic order, because their index
// becomes the het anchor id.
//
// Returns the number of sites dropped.
uint64_t collapseDuplicateLoci(vector<Assembler::CigarSnpSite>& sites);


// Remove arm members whose (read, position) is already held by an existing
// anchor, and report how many were removed.
//
// A marker anchor groups reads sharing a k=50 k-mer centred on its position,
// and that k-mer spans the SNP base -- so reads carrying different alleles were
// never in the same marker anchor. An anchor already at a member's position
// therefore isolates that member's allele, necessarily and not just usually,
// and building a het anchor over it would clone it and then compete with it for
// a position shasta2 lets only one anchor hold.
//
// `occupied` maps (orientedReadId.getValue() << 32 | position) to nothing --
// only membership is used. Build it with buildOccupiedPositions().
//
// Site positions and stored anchor positions are in the SAME frame: a het
// anchor is a zero-length marker whose stored position is the SNP base itself
// (see Shasta2Anchors::appendHetAnchorPair), so no offset is applied here. This
// used to take a hetKHalf argument for the old 2-base het marker; getting that
// offset wrong matched nothing and removed the wrong members, so the frames are
// now identical by construction rather than by a caller-supplied constant.
uint64_t dropAlreadyAnchoredArmMembers(
    vector<Assembler::CigarSnpSite>& sites,
    const std::unordered_map<uint64_t, Shasta2AnchorId>& occupied);


// (orientedRead, storedPosition) -> the anchor holding it, over every anchor
// currently in the store.
std::unordered_map<uint64_t, Shasta2AnchorId> buildOccupiedPositions(
    const Shasta2Anchors& anchors);

}

#endif
