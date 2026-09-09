#ifndef DINARA_ABPOA_HANDLE_HPP
#define DINARA_ABPOA_HANDLE_HPP

// RAII wrapper around an abPOA aligner plus its parameter block.
//
// This lived in WindowIntervalPoa.hpp until the window subsystem was removed.
// It has nothing to do with windows -- it is just the abPOA setup every caller
// needs -- so it was extracted here rather than deleted with its old home.

#include <abpoa.h>

namespace dinara {

// abPOA handle wrapper. Thread-local; reset per interval by abpoa_reset. Long-
// read-tuned affine scoring, adaptive banding, minimizer seeding + partitioning
// (so wide gaps become many tiny window POAs, not one huge DP), row order
// preserved so row 0 is always the backbone.
struct IpoaAbHandle {
    abpoa_t* ab = nullptr;
    abpoa_para_t* abpt = nullptr;
    IpoaAbHandle() {
        ab = abpoa_init();
        abpt = abpoa_init_para();
        abpt->align_mode = ABPOA_GLOBAL_MODE;
        abpt->gap_mode = ABPOA_AFFINE_GAP;
        abpt->match = 2;
        abpt->mismatch = 4;
        abpt->gap_open1 = 4;
        abpt->gap_ext1 = 2;
        abpt->gap_open2 = 0;
        abpt->gap_ext2 = 0;
        abpt->wb = 10;
        abpt->wf = 0.01;
        // Minimizer-based seeding + partitioning (shasta2 LocalAssembly values).
        // With disable_seeding=1 AND progressive_poa=0, abpoa takes the
        // abpoa_poa() branch: ONE monolithic O(qlen * graphNodes) DP over the
        // whole interval. On wide anchor-less gaps (up to ~4.7 kb x ~50 reads)
        // that single DP was ~500 ms and dominated het time. Leaving seeding on
        // (as shasta2 does) makes abpoa take the abpoa_anchor_poa() branch: it
        // finds shared minimizers, partitions each sequence into small windows
        // at those anchors, and runs a tiny DP per window -- exactly the
        // sub-tiling we'd otherwise hand-roll, done inside abpoa. Row order is
        // still input order (abpoa_anchor_poa adds each read at its ORIGINAL
        // index read_id = read_id_map[_i], so msa_base[0] stays the backbone).
        abpt->disable_seeding = 0;
        abpt->w = 6;
        abpt->k = 9;
        abpt->min_w = 10;
        abpt->progressive_poa = 0;
        abpt->sort_input_seq = 0;   // keep row 0 = backbone
        abpt->out_msa = 1;
        abpt->out_cons = 0;
        abpt->ret_cigar = 1;
        abpoa_post_set_para(abpt);
    }
    ~IpoaAbHandle() { if (ab) abpoa_free(ab); if (abpt) abpoa_free_para(abpt); }
};

}

#endif
