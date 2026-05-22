// Minimal base class for snarl finders, replacing vg's snarls.hpp.
// No protobuf, no SnarlManager — just the traverse_decomposition interface.

#ifndef DINARA_SNARL_FINDER_BASE_HPP_INCLUDED
#define DINARA_SNARL_FINDER_BASE_HPP_INCLUDED

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <functional>

namespace vg {

using handlegraph::HandleGraph;
using handlegraph::RankedHandleGraph;
using handlegraph::handle_t;
using handlegraph::nid_t;

/// Minimal base class that declares the traverse_decomposition interface.
/// Replaces HandleGraphSnarlFinder from vg's snarls.hpp without pulling
/// in protobuf or SnarlManager.
class HandleGraphSnarlFinder {
protected:
    const HandleGraph* graph;

public:
    HandleGraphSnarlFinder(const HandleGraph* graph) : graph(graph) {}
    virtual ~HandleGraphSnarlFinder() = default;

    /// Visit all snarls and chains via nested begin/end callbacks.
    /// The first call is always begin_chain (the root chain).
    /// Level-0 top snarls are those whose begin_snarl fires at depth 1.
    virtual void traverse_decomposition(
        const std::function<void(handle_t)>& begin_chain,
        const std::function<void(handle_t)>& end_chain,
        const std::function<void(handle_t)>& begin_snarl,
        const std::function<void(handle_t)>& end_snarl) const = 0;
};

} // namespace vg

#endif
