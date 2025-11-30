#ifndef DINARA_ASSEMBLE_MARKER_GRAPH_PATH_HPP
#define DINARA_ASSEMBLE_MARKER_GRAPH_PATH_HPP

#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "span.hpp"

namespace dinara {

    class AssembledSegment;
    class Reads;

    void assembleMarkerGraphPath(
        uint64_t readRepresentation,
        uint64_t k,
        const Reads& reads,
        const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
        const MarkerGraph&,
        const span<const MarkerGraph::EdgeId>& markerGraphPath,
        bool storeCoverageData,
        AssembledSegment& assembledSegment);

}



#endif
