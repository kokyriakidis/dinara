// Shasta2.
#include "Shasta2AnchorGraph.hpp"

#include <boost/graph/iteration_macros.hpp>

#include <fstream>
#include <set>

using namespace dinara;
using namespace std;

void Shasta2AnchorGraph::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }
    gfa << "H\tVN:Z:1.0\n";

    // Write vertices (Anchors)
    BGL_FORALL_VERTICES(v, *this, Shasta2AnchorGraph) {
        // Use '*' for sequence, assign length = 1 as placeholder
        gfa << "S\t" << v << "\t*\tLN:i:1\n";
    }

    // Write edges
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        const auto& edge = (*this)[e];
        if(!edge.useForAssembly) {
            continue;
        }

        // Determine orientation from the first supporting read.
        string sourceOrientation = "+";
        string targetOrientation = "+";

        if(!edge.anchorPair.orientedReadIds.empty()) {
            if(edge.anchorPair.orientedReadIds[0].getStrand() == 1) {
                sourceOrientation = "-";
                targetOrientation = "-";
            }
        }

        gfa << "L\t" << source(e, *this) << "\t" << sourceOrientation << "\t"
            << target(e, *this) << "\t" << targetOrientation << "\t0M\n";
    }
}



void Shasta2AnchorGraph::writeBubbleFinderGraph(const string& fileName, bool useForAssemblyOnly) const
{
    ofstream graphFile(fileName);
    if(!graphFile) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }

    std::set<vertex_descriptor> verticesInOutput;
    uint64_t edgeCount = 0;
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        const auto& edge = (*this)[e];
        if(useForAssemblyOnly && !edge.useForAssembly) {
            continue;
        }
        verticesInOutput.insert(source(e, *this));
        verticesInOutput.insert(target(e, *this));
        ++edgeCount;
    }

    graphFile << verticesInOutput.size() << " " << edgeCount << "\n";
    BGL_FORALL_EDGES(e, *this, Shasta2AnchorGraph) {
        const auto& edge = (*this)[e];
        if(useForAssemblyOnly && !edge.useForAssembly) {
            continue;
        }
        graphFile << source(e, *this) << " " << target(e, *this) << "\n";
    }
}
