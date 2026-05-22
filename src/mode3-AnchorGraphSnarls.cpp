// Snarl decomposition of the AnchorGraph using vg's IntegratedSnarlFinder.

#include "mode3-AnchorGraph.hpp"
#include "dinara_handle_graph.hpp"
#include "integrated_snarl_finder.hpp"
#include "timestamp.hpp"

#include <boost/graph/iteration_macros.hpp>

#include <iostream>

using namespace dinara;
using namespace mode3;
using std::cout;
using std::endl;

void AnchorGraph::findSnarls(
    vector<pair<uint64_t, uint64_t>>& topSnarls) const
{
    const AnchorGraph& g = *this;
    topSnarls.clear();

    // Build a DinaraHandleGraph from the AnchorGraph.
    // Node IDs are localAnchorId + 1 (1-based).
    const uint64_t n = vertexDescriptors.size();
    dinara::DinaraHandleGraph hg;
    hg.nodes.resize(n);

    // Set node lengths. Use 1 as a placeholder — the snarl finder
    // uses lengths only for tie-breaking, not for correctness.
    for (uint64_t i = 0; i < n; i++) {
        hg.nodes[i].length = 1;
    }

    // Add edges. Each directed edge A -> B in the AnchorGraph becomes
    // an edge from node (A+1) forward to node (B+1) forward in the handle graph.
    BGL_FORALL_EDGES_T(e, g, AnchorGraph) {
        if (g[e].isNonTransitiveReductionEdge) continue;
        const auto& sv = g[boost::source(e, g)];
        const auto& tv = g[boost::target(e, g)];
        hg.edges.push_back({
            static_cast<handlegraph::nid_t>(sv.localAnchorId + 1), false,
            static_cast<handlegraph::nid_t>(tv.localAnchorId + 1), false
        });
    }

    cout << timestamp << "findSnarls: " << n << " nodes, "
         << hg.edges.size() << " edges." << endl;

    // Run the snarl finder.
    vg::IntegratedSnarlFinder finder(hg);

    // Traverse the decomposition. Track depth to identify level-0 top snarls.
    // Structure: root chain [ snarl [ chain [ snarl ... ] ] ]
    // Depth 1 = inside root chain. Depth 2 = inside a top snarl.
    int depth = 0;
    handlegraph::handle_t snarlStart;

    finder.traverse_decomposition(
        [&](handlegraph::handle_t) { depth++; },   // begin_chain
        [&](handlegraph::handle_t) { depth--; },   // end_chain
        [&](handlegraph::handle_t h) {              // begin_snarl
            depth++;
            if (depth == 2) {
                snarlStart = h;
            }
        },
        [&](handlegraph::handle_t h) {              // end_snarl
            if (depth == 2) {
                // h is the boundary node at the end of this top snarl.
                uint64_t id0 = hg.get_id(snarlStart) - 1;  // back to 0-based local anchor id
                uint64_t id1 = hg.get_id(h) - 1;
                topSnarls.push_back({id0, id1});
            }
            depth--;
        }
    );

    cout << timestamp << "findSnarls: found " << topSnarls.size()
         << " top-level snarls." << endl;
}
