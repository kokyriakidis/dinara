#pragma once

#include "Shasta2AssemblyGraph.hpp"

#include "span.hpp"

namespace dinara {
    class Shasta2AssemblyGraphPostprocessor;
}

// Shasta2AssemblyGraph functionality needed only during postprocessing.
// It is used in the http server and in the Python API.
class dinara::Shasta2AssemblyGraphPostprocessor : public Shasta2AssemblyGraph {
public:
    Shasta2AssemblyGraphPostprocessor(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        const Shasta2AssemblyGraphOptions&,
        const string& assemblyStage);

    // Map from vertex id to vertex_descriptor.
    std::map<uint64_t, vertex_descriptor> vertexMap;

    // Map from edge id to edge_descriptor.
    std::map<uint64_t, edge_descriptor> edgeMap;

    // Python-callable function to get an id given an edge-descriptor.
    uint64_t getId(edge_descriptor e) const
    {
        return (*this)[e].id;
    }

    // Annotations of where each AnchorId is used in the current state
    // of Shasta2AssemblyGraph.
    class Annotation {
    public:
        Shasta2AnchorId anchorId;
        vertex_descriptor v = null_vertex();
        edge_descriptor e;
        uint64_t step = invalid<uint64_t>;
        bool isAnchorIdA;

        bool operator<(const Annotation& that) const
        {
            return anchorId < that.anchorId;
        }

        Annotation(Shasta2AnchorId anchorId) :
            anchorId(anchorId)
        {}
        Annotation(Shasta2AnchorId anchorId, vertex_descriptor v) :
            anchorId(anchorId),
            v(v)
        {}
        Annotation(
            Shasta2AnchorId anchorId,
            edge_descriptor e,
            uint64_t step,
            bool isAnchorIdA) :
            anchorId(anchorId),
            e(e),
            step(step),
            isAnchorIdA(isAnchorIdA)
        {}
    };

    // The vector of annotations is kept sorted so we can do searches.
    vector<Annotation> annotations;
    void computeAnnotations();
    span<const Annotation> getAnnotations(Shasta2AnchorId) const;

    // Return true if an AnchorId has one or more vertex annotations.
    bool hasVertexAnnotation(Shasta2AnchorId) const;

    // Find the edges that an AnchorId has annotations for.
    void findAnnotationEdges(Shasta2AnchorId, vector<edge_descriptor>&) const;
};

