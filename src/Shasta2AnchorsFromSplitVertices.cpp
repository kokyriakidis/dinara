// Shasta2AnchorsFromSplitVertices.cpp
//
// This functionality has been removed.
// The file is kept as a stub to avoid build system changes.

#include "Shasta2AnchorsFromSplitVertices.hpp"

#include <stdexcept>

using namespace dinara;
using namespace std;

std::shared_ptr<Shasta2Anchors> dinara::createShasta2AnchorsFromSplitVertices(
    const MappedMemoryOwner&,
    const Reads&,
    uint64_t,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>&,
    const MarkerGraph&,
    const ReadGraph&,
    uint64_t,
    uint64_t,
    uint64_t,
    const MemoryMapped::Vector<AlignmentData>*,
    const MemoryMapped::VectorOfVectors<uint32_t, uint32_t>*)
{
    throw runtime_error("createShasta2AnchorsFromSplitVertices has been removed.");
}
