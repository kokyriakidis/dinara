#pragma once

// HetAnchorK.hpp
//
// Experimental toggle for the het/hom SNP anchor marker length k.
//
// Default (k=2): each het/hom anchor is a 2-base marker. A member read stores
// the position of the SNP's PREDECESSOR base (pos-1) and shasta2 re-derives the
// 2-mer [predBase, alleleBase]; the stored midpoint is predReadPos + k/2 =
// predReadPos + 1 = the SNP read position. Export subtracts k/2 = 1 uniformly.
//
// Experimental (k=0, DINARA_HET_K=0): each het/hom anchor is a zero-length
// position marker. A member read stores the EXACT SNP base position and there
// is no 2-mer to keep consistent, so the k-mer-adjacency pin guard is relaxed.
// The stored midpoint is snpReadPos + k/2 = snpReadPos + 0 = the SNP read
// position -- IDENTICAL to the k=2 stored midpoint. Export subtracts k/2 = 0.
//
// Because the stored midpoint is the SNP read position in BOTH modes, every
// consumer that works in the stored-position frame (anchor ordering, primary
// collision, backward-edge drop, the monotonicity verifier) is unchanged; only
// (a) the emission pin column, (b) hetKHalf, and (c) the export shift switch
// with k. The shasta2 external-anchor export must be invoked with --k equal to
// hetAnchorK() so the loader re-adds the matching k/2.
//
// Shasta2 developer (Paolo Carnevali) confirmed k=0 is valid: the local-assembly
// interval math [x + k/2, ...) holds with k/2 = 0, so a k=0 anchor position is
// simply the first base of / first base after a local assembly, and using the
// het allele's exact position is fine.
//
// This is EXPERIMENTAL and unproven. The default stays k=2 so the working path
// is untouched; set DINARA_HET_K=0 to opt in.

#include <cstdint>
#include <cstdlib>

namespace dinara {

// Marker length k for het/hom anchors. Read once from DINARA_HET_K; only 0 and
// 2 are supported (any other value falls back to the k=2 default).
inline uint32_t hetAnchorK() {
    static const uint32_t k = [] {
        const char* e = std::getenv("DINARA_HET_K");
        if(e != nullptr && e[0] == '0' && e[1] == '\0') return uint32_t(0);
        return uint32_t(2);
    }();
    return k;
}

// Half marker length (k/2) used as the stored-midpoint offset and the uniform
// export shift. 1 for k=2, 0 for k=0.
inline uint32_t hetAnchorKHalf() {
    return hetAnchorK() / 2u;
}

}
