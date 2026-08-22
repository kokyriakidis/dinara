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
// The marker->member conversion (validated on real data at 100% allele and
// oriented-2-base-marker consistency; see exp_coord_validate in the
// hifiasmCandidates submodule):
//   canonMid    = (key >> (2*mid)) & 3            middle (allele) base, mid=(k-1)/2
//   o           = (rawFwd[pos+mid] == canonMid) ? 0 : 1     strand reading canonical
//   rawPosition = o ? (L - (pos+mid) - 2) : (pos + mid - 1) predBase, oriented frame
//   member      = (OrientedReadId(readId, o), rawPosition)
// appendHetAnchorPair stores rawPosition + hetK/2, i.e. the SNP (allele) column.
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
