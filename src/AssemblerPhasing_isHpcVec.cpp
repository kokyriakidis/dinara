#include "AssemblerPhasing.hpp"

using namespace dinara;

// Helper for homopolymer filtering
// Corresponds to Hifiasm's is_hpc_vec
bool isHpcVec(AssemblerPhasing::SnpStats& s, const PhasingConfig& config) {
    // If we had homopolymer-tagged evidence, we would subtract it here.
    // s.occ_0 -= s.homopolymer_evidence_ref;
    // s.occ_1 -= s.homopolymer_evidence_alt;
    
    // Check if remaining (non-HP) evidence is sufficient
    bool insufficient = false;
    
    // Logic from Hifiasm Line 9396:
    // if((ai->occ_0 < 2 || ai->occ_1 < 2) || (!(ai->occ_0 >= asm_opt.s_hap_cov && ai->occ_1 >= asm_opt.infor_cov))) f = 1;

    if (s.occ_0 < 2 || s.occ_1 < 2) insufficient = true;
    else if (!(s.occ_0 >= (uint32_t)config.s_hap_cov && s.occ_1 >= (uint32_t)config.infor_cov)) insufficient = true;
    
    // Restore counts if we modified them (we didn't here, but keeping pattern)
    // s.occ_0 += ...
    
    return insufficient; // Returns true if it fails the check (is "HPC-like" or weak)
}
