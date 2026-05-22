// Adapter that wraps a Dinara graph as a handlegraph::RankedHandleGraph
// for use with IntegratedSnarlFinder.
//
// This is a skeleton — fill in the methods when integrating with a
// specific Dinara graph type (ReadGraph, AssemblyGraph2, etc.).

#ifndef DINARA_HANDLE_GRAPH_HPP_INCLUDED
#define DINARA_HANDLE_GRAPH_HPP_INCLUDED

#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <cassert>
#include <string>
#include <vector>
#include <functional>

namespace dinara {

using handlegraph::handle_t;
using handlegraph::nid_t;
using handlegraph::edge_t;

/// Adapter that presents a Dinara graph as a RankedHandleGraph.
///
/// Node IDs are expected to be dense 1-based integers (1..N).
/// Each node has two orientations (forward/reverse).
/// Edges connect oriented node sides.
///
/// Subclass or modify this to wrap your specific graph type.
class DinaraHandleGraph : public handlegraph::RankedHandleGraph {
public:
    // ── Data to populate before calling the snarl finder ──

    struct Node {
        uint64_t length = 0;  // sequence length in bases
    };

    struct Edge {
        nid_t from;
        bool from_reverse;
        nid_t to;
        bool to_reverse;
    };

    std::vector<Node> nodes;  // 0-indexed, node i has ID i+1
    std::vector<Edge> edges;

    // ── HandleGraph interface ──

    bool has_node(nid_t node_id) const override {
        return node_id >= 1 && node_id <= static_cast<nid_t>(nodes.size());
    }

    handle_t get_handle(const nid_t& node_id, bool is_reverse = false) const override {
        return handlegraph::number_bool_packing::pack(static_cast<uint64_t>(node_id), is_reverse);
    }

    nid_t get_id(const handle_t& handle) const override {
        return static_cast<nid_t>(handlegraph::number_bool_packing::unpack_number(handle));
    }

    bool get_is_reverse(const handle_t& handle) const override {
        return handlegraph::number_bool_packing::unpack_bit(handle);
    }

    handle_t flip(const handle_t& handle) const override {
        return handlegraph::number_bool_packing::toggle_bit(handle);
    }

    size_t get_length(const handle_t& handle) const override {
        nid_t id = get_id(handle);
        assert(id >= 1 && id <= static_cast<nid_t>(nodes.size()));
        return nodes[id - 1].length;
    }

    std::string get_sequence(const handle_t& /*handle*/) const override {
        // Snarl finder doesn't need actual sequences.
        return "";
    }

    size_t get_node_count() const override {
        return nodes.size();
    }

    nid_t min_node_id() const override {
        return nodes.empty() ? 0 : 1;
    }

    nid_t max_node_id() const override {
        return static_cast<nid_t>(nodes.size());
    }

    bool follow_edges_impl(const handle_t& handle, bool go_left,
        const std::function<bool(const handle_t&)>& iteratee) const override {
        // go_left=false: follow edges leaving the right side of handle
        // go_left=true:  follow edges entering the left side of handle
        //                (equivalently, leaving the right side of flip(handle))
        nid_t id = get_id(handle);
        bool rev = get_is_reverse(handle);

        for (const auto& e : edges) {
            handle_t target;
            bool match = false;

            if (!go_left) {
                // Looking for edges from (id, rev) going right
                if (e.from == id && e.from_reverse == rev) {
                    target = get_handle(e.to, e.to_reverse);
                    match = true;
                } else if (e.to == id && e.to_reverse != rev) {
                    // Reverse edge: to-side reversed matches from-side
                    target = get_handle(e.from, !e.from_reverse);
                    match = true;
                }
            } else {
                // Looking for edges into (id, rev) from the left
                if (e.to == id && e.to_reverse == rev) {
                    target = get_handle(e.from, e.from_reverse);
                    match = true;
                } else if (e.from == id && e.from_reverse != rev) {
                    target = get_handle(e.to, !e.to_reverse);
                    match = true;
                }
            }

            if (match) {
                if (!iteratee(target)) return false;
            }
        }
        return true;
    }

    bool for_each_handle_impl(
        const std::function<bool(const handle_t&)>& iteratee,
        bool /*parallel*/ = false) const override {
        for (nid_t i = 1; i <= static_cast<nid_t>(nodes.size()); i++) {
            if (!iteratee(get_handle(i, false))) return false;
        }
        return true;
    }

    // ── RankedHandleGraph interface ──
    // Dense 1-based ranks matching node IDs.

    size_t id_to_rank(const nid_t& node_id) const override {
        return static_cast<size_t>(node_id);
    }

    nid_t rank_to_id(const size_t& rank) const override {
        return static_cast<nid_t>(rank);
    }

    // handle_to_rank and rank_to_handle use the default implementations
    // from RankedHandleGraph which delegate to id_to_rank/rank_to_id.
};

} // namespace dinara

#endif
