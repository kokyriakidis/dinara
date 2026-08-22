#pragma once

// Snpmer het-anchor producer for the external Shasta2 anchor export.
//
// myloasm detects, per read, the SNPmer markers: k-mers whose middle base marks
// a heterozygous site. Two reads carrying the same allele of the same site share
// the marker's CANONICAL k-mer key. This producer groups all reads by canonical
// key (one key = one allele of one site), filters paralogs/duplicates, and
// appends each allele as a Shasta2 k=2 het anchor via
// Shasta2Anchors::appendHetAnchorPair -- so snpmer anchors are exported in
// Shasta2ExternalAnchors alongside the primary and window-detected het anchors.
//
// This is the "no window" direct-append path: snpmer sites carry no backbone
// frame or bracketing hom anchors, so they are NOT wired as anchor-graph bubble
// edges; they are export-only additional anchors.
//
// The marker->member conversion, with snpCol = pos + mid the forward SNP
// column (mid = (k-1)/2):
//   rawPosition = snpCol - 1                    predBase in the read's frame
//   member      = (OrientedReadId(readId, 0), rawPosition)
// appendHetAnchorPair stores rawPosition + hetK/2 = snpCol (the allele column)
// and generates the reverse-complement (strand 1) placement itself.
//
// Every member is placed in the read's OWN sequencing frame (strand 0). We do
// NOT pre-orient members by a per-SNP canonical strand: a read has one physical
// orientation, but two SNPs on it can each have their canonical allele readable
// on OPPOSITE strands, so per-SNP strands would put both a strand-0 and a
// strand-1 member on one read. appendHetAnchorPair mirrors each member, so the
// RC of the strand-1 member would land on strand 0 on top of a neighboring
// strand-0 canonical member -- two anchors at the same oriented-read position,
// which crashes shasta2's journey builder. Allele separation needs no
// orientation: it comes from the strand-invariant canonical key (allele0 and
// allele1 have distinct keys).
//
// Reads are matched myloasm<->dinara BY NAME (myloasm re-parses the input files
// and may index a different subset/order than dinara's store).

#include <cstdint>
#include <string>
#include <vector>

namespace dinara {
    class Assembler;

    // Run the snpmer producer over the given input files and append the
    // resulting het anchors to assembler.shasta2Anchors. Returns the number of
    // allele anchors appended (2 per biallelic site kept). Must be called AFTER
    // all primary anchors and any window het anchors are appended (so
    // hetAnchorFirstId still marks the primary/het boundary) and BEFORE
    // writeExternalAnchors.
    //
    // minMembersPerAllele: alleles with fewer members are skipped
    //   (appendHetAnchorPair requires >= 2; a lone read is a spurious branch).
    // requireBiallelic: if true, only sites where BOTH alleles survive filtering
    //   are emitted (a real het bubble); if false, a surviving single allele is
    //   also emitted.
    uint64_t appendSnpmerHetAnchors(
        Assembler& assembler,
        const std::vector<std::string>& inputFileNames,
        uint64_t threadCount,
        uint64_t minMembersPerAllele,
        bool requireBiallelic);
}
