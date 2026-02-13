// Assembler method for Verkko-style directed anchor graph resolution.
// NOTE: The primary pipeline now builds the DAG from mode3::Anchors
// directly in main.cpp. This standalone entry point is kept for
// compatibility but requires pre-built mode3::Anchors (stored externally).

#include "Assembler.hpp"
#include "mode3-DirectedAnchorGraph.hpp"
#include "timestamp.hpp"

using namespace dinara;


void Assembler::runDirectedAnchorGraphResolution()
{
    cout << timestamp << "runDirectedAnchorGraphResolution begins." << endl;

    // This standalone method is no longer the primary entry point.
    // The expanded pipeline in main.cpp builds mode3::Anchors,
    // computes journeys, and feeds them directly to
    // DirectedAnchorGraph::buildFromAnchors().
    throw runtime_error(
        "runDirectedAnchorGraphResolution: this standalone entry point "
        "is deprecated. Use the expanded pipeline in main.cpp which builds "
        "the DAG directly from mode3::Anchors.");
}
