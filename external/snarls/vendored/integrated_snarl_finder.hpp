// Vendored from vg (https://github.com/vgteam/vg)
// Stripped: removed find_snarls/find_snarls_parallel (protobuf/SnarlManager),
//           replaced snarls.hpp with snarl_finder_base.hpp,
//           removed bdsg overlay dependency (takes RankedHandleGraph directly).
// License: MIT

#ifndef VG_INTEGRATED_SNARL_FINDER_HPP_INCLUDED
#define VG_INTEGRATED_SNARL_FINDER_HPP_INCLUDED

#include "snarl_finder_base.hpp"

#include <functional>
#include <vector>
#include <unordered_map>
#include <utility>

namespace vg {

using namespace std;

class IntegratedSnarlFinder : public HandleGraphSnarlFinder {
private:
    class MergedAdjacencyGraph;

    const std::unordered_map<nid_t, size_t> extra_node_weight;

    void traverse_computed_decomposition(MergedAdjacencyGraph& cactus,
        const MergedAdjacencyGraph& forest,
        vector<pair<size_t, vector<handle_t>>>& longest_paths,
        unordered_map<handle_t, handle_t>& towards_deepest_leaf,
        vector<pair<size_t, handle_t>>& longest_cycles,
        unordered_map<handle_t, handle_t>& next_along_cycle,
        const function<void(handle_t)>& begin_chain, const function<void(handle_t)>& end_chain,
        const function<void(handle_t)>& begin_snarl, const function<void(handle_t)>& end_snarl) const;

public:
    /// Construct a snarl finder on a RankedHandleGraph.
    /// The graph must outlive this object.
    IntegratedSnarlFinder(const RankedHandleGraph& graph,
        const std::unordered_map<nid_t, size_t>& extra_node_weight = {});

    virtual ~IntegratedSnarlFinder() = default;

    void traverse_decomposition(
        const function<void(handle_t)>& begin_chain,
        const function<void(handle_t)>& end_chain,
        const function<void(handle_t)>& begin_snarl,
        const function<void(handle_t)>& end_snarl) const override;
};

}

#endif
