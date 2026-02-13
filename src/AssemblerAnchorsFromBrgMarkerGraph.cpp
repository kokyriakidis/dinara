// ============================================================================
// Anchor creation from BRG-aware (self-RC) marker graph vertices.
//
// When marker graph vertices are created from the BRG with forward/RC kmer
// merging (createMarkerGraphVerticesFromBrg), each vertex is self-RC: it
// contains markers from both strands of the same physical reads.
//
// The standard vertex-based anchor code (mode3-AnchorFromMarkerGraph.cpp)
// rejects such vertices because they have "duplicate ReadIds" (same ReadId,
// different strands).
//
// This file provides createAnchorsFromBrgMarkerGraphVertices(), which handles
// self-RC vertices correctly by extracting only strand-0 markers for the
// forward anchor and deriving the RC anchor by flipping.  Non-self-RC vertices
// (e.g. palindromic kmers where both map to the same vertex, or vertices that
// happen to have separate RC partners) are handled with the standard logic.
// ============================================================================

#include "Assembler.hpp"
#include "findMarkerId.hpp"
#include "Marker.hpp"
#include "MarkerGraph.hpp"
#include "mode3-Anchor.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
#include "timestamp.hpp"

using namespace dinara;
using namespace mode3;
using namespace std;


shared_ptr<Anchors> Assembler::createAnchorsFromBrgMarkerGraphVertices(
    uint64_t minAnchorCoverage,
    uint64_t maxAnchorCoverage,
    uint64_t threadCount)
{
    performanceLog << timestamp
        << "Begin createAnchorsFromBrgMarkerGraphVertices." << endl;

    // Check prerequisites.
    reads->checkReadsAreOpen();
    checkMarkersAreOpen();
    DINARA_ASSERT(markerGraph.vertices().isOpen());
    DINARA_ASSERT(markerGraph.vertexTable.isOpen);
    DINARA_ASSERT(markerGraph.reverseComplementVertex.isOpen);

    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    using VertexId = MarkerGraph::VertexId;
    using Interval = AnchorMarkerInterval;

    const auto& mgVertices = markerGraph.vertices();
    const VertexId vertexCount = markerGraph.vertexCount();

    // We will collect forward + RC anchor pairs into this vector.
    vector<vector<Interval>> anchorsExplicit;
    anchorsExplicit.reserve(vertexCount);  // upper bound

    uint64_t selfRcCount = 0;
    uint64_t normalCount = 0;
    uint64_t skippedCoverage = 0;
    uint64_t skippedDuplicate = 0;

    for(VertexId vertexId = 0; vertexId < vertexCount; ++vertexId) {

        const VertexId rcVertexId = markerGraph.reverseComplementVertex[vertexId];

        // Only process canonical vertices (vertexId <= rcVertexId).
        if(vertexId > rcVertexId) {
            continue;
        }

        const auto vertexMarkerIds = mgVertices[vertexId];
        const bool isSelfRc = (vertexId == rcVertexId);

        if(isSelfRc) {
            // ---- Self-RC vertex: extract strand-0 markers only ----

            vector<Interval> fwdAnchor;
            fwdAnchor.reserve(vertexMarkerIds.size() / 2 + 1);

            for(const MarkerId markerId : vertexMarkerIds) {
                OrientedReadId orientedReadId;
                uint32_t ordinal;
                tie(orientedReadId, ordinal) = findMarkerId(markerId);
                if(orientedReadId.getStrand() == 0) {
                    fwdAnchor.push_back(Interval(orientedReadId, ordinal));
                }
            }

            // Sort by OrientedReadId.
            sort(fwdAnchor.begin(), fwdAnchor.end(),
                [](const Interval& a, const Interval& b) {
                    return a.orientedReadId < b.orientedReadId;
                });

            // Deduplicate by OrientedReadId (keep first occurrence).
            fwdAnchor.erase(
                unique(fwdAnchor.begin(), fwdAnchor.end(),
                    [](const Interval& a, const Interval& b) {
                        return a.orientedReadId == b.orientedReadId;
                    }),
                fwdAnchor.end());

            // Check coverage.
            if(fwdAnchor.size() < minAnchorCoverage ||
               fwdAnchor.size() > maxAnchorCoverage) {
                ++skippedCoverage;
                continue;
            }

            // Build RC anchor by flipping.
            vector<Interval> rcAnchor;
            rcAnchor.reserve(fwdAnchor.size());
            for(const auto& interval : fwdAnchor) {
                OrientedReadId rcOrientedReadId = interval.orientedReadId;
                rcOrientedReadId.flipStrand();
                const uint64_t markerCount = markers->size(rcOrientedReadId.getValue());
                const uint32_t rcOrdinal = uint32_t(markerCount) - 1 - interval.ordinal0;
                rcAnchor.push_back(Interval(rcOrientedReadId, rcOrdinal));
            }
            sort(rcAnchor.begin(), rcAnchor.end(),
                [](const Interval& a, const Interval& b) {
                    return a.orientedReadId < b.orientedReadId;
                });

            anchorsExplicit.push_back(std::move(fwdAnchor));
            anchorsExplicit.push_back(std::move(rcAnchor));
            ++selfRcCount;

        } else {
            // ---- Standard non-self-RC vertex pair ----

            // Check coverage on the whole vertex.
            if(vertexMarkerIds.size() < minAnchorCoverage ||
               vertexMarkerIds.size() > maxAnchorCoverage) {
                ++skippedCoverage;
                continue;
            }

            // Check for duplicate ReadIds.
            if(markerGraph.vertexHasDuplicateReadIds(vertexId, *markers)) {
                ++skippedDuplicate;
                continue;
            }

            // Gather all markers.
            vector<Interval> fwdAnchor;
            fwdAnchor.reserve(vertexMarkerIds.size());
            for(const MarkerId markerId : vertexMarkerIds) {
                OrientedReadId orientedReadId;
                uint32_t ordinal;
                tie(orientedReadId, ordinal) = findMarkerId(markerId);
                fwdAnchor.push_back(Interval(orientedReadId, ordinal));
            }

            // Sort by OrientedReadId.
            sort(fwdAnchor.begin(), fwdAnchor.end(),
                [](const Interval& a, const Interval& b) {
                    return a.orientedReadId < b.orientedReadId;
                });

            // Deduplicate by OrientedReadId.
            fwdAnchor.erase(
                unique(fwdAnchor.begin(), fwdAnchor.end(),
                    [](const Interval& a, const Interval& b) {
                        return a.orientedReadId == b.orientedReadId;
                    }),
                fwdAnchor.end());

            // Re-check coverage.
            if(fwdAnchor.size() < minAnchorCoverage ||
               fwdAnchor.size() > maxAnchorCoverage) {
                ++skippedCoverage;
                continue;
            }

            // Build RC anchor by flipping.
            vector<Interval> rcAnchor;
            rcAnchor.reserve(fwdAnchor.size());
            for(const auto& interval : fwdAnchor) {
                OrientedReadId rcOrientedReadId = interval.orientedReadId;
                rcOrientedReadId.flipStrand();
                const uint64_t markerCount = markers->size(rcOrientedReadId.getValue());
                const uint32_t rcOrdinal = uint32_t(markerCount) - 1 - interval.ordinal0;
                rcAnchor.push_back(Interval(rcOrientedReadId, rcOrdinal));
            }
            sort(rcAnchor.begin(), rcAnchor.end(),
                [](const Interval& a, const Interval& b) {
                    return a.orientedReadId < b.orientedReadId;
                });

            anchorsExplicit.push_back(std::move(fwdAnchor));
            anchorsExplicit.push_back(std::move(rcAnchor));
            ++normalCount;
        }
    }

    cout << timestamp << "BRG marker graph anchor creation: "
         << selfRcCount << " self-RC vertices, "
         << normalCount << " standard vertex pairs, "
         << skippedCoverage << " skipped (coverage), "
         << skippedDuplicate << " skipped (duplicate ReadIds)." << endl;
    cout << timestamp << "Total anchors (including RC): "
         << anchorsExplicit.size() << endl;

    // Create Anchors using the explicit-intervals constructor.
    auto anchors = make_shared<Anchors>(
        MappedMemoryOwner(*this),
        getReads(),
        assemblerInfo->k,
        *markers,
        anchorsExplicit,
        /*ordinalOffset*/ 0,
        threadCount);

    anchors->computeJourneys(threadCount);

    performanceLog << timestamp
        << "End createAnchorsFromBrgMarkerGraphVertices." << endl;

    return anchors;
}
