// Shasta2
#include "Shasta2AssemblyGraphPostprocessor.hpp"
#include "deduplicate.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>



Shasta2AssemblyGraphPostprocessor::Shasta2AssemblyGraphPostprocessor(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const Shasta2AssemblyGraphOptions& options,
    const string& assemblyStage) :
    Shasta2AssemblyGraph(anchors, journeys, options, assemblyStage)
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;

    // Fill in the vertex map.
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        vertexMap.insert(make_pair(assemblyGraph[v].id, v));
    }

    // Fill in the edge map.
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        edgeMap.insert(make_pair(assemblyGraph[e].id, e));
    }

    computeAnnotations();
}


void Shasta2AssemblyGraphPostprocessor::computeAnnotations()
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;
    annotations.clear();

    // Vertices.
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        annotations.emplace_back(assemblyGraph[v].anchorId, v);
    }

    // Edges.
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
        for(uint64_t step=0; step<edge.size(); step++) {
            const Shasta2AnchorPair& anchorPair = edge[step].anchorPair;
            annotations.emplace_back(anchorPair.anchorIdA, e, step, true);
            annotations.emplace_back(anchorPair.anchorIdB, e, step, false);
        }
    }

    sort(annotations.begin(), annotations.end());
}



span<const Shasta2AssemblyGraphPostprocessor::Annotation>
    Shasta2AssemblyGraphPostprocessor::getAnnotations(Shasta2AnchorId anchorId) const
{
    const auto p = std::equal_range(annotations.begin(), annotations.end(), Annotation(anchorId));
    return span<const Annotation>(p.first, p.second);
}



bool Shasta2AssemblyGraphPostprocessor::hasVertexAnnotation(Shasta2AnchorId anchorId) const
{
    const auto annotations = getAnnotations(anchorId);
    for(const Annotation& annotation: annotations) {
        if(annotation.v != null_vertex()) {
            return true;
        }
    }
    return false;
}



void Shasta2AssemblyGraphPostprocessor::findAnnotationEdges(
    Shasta2AnchorId anchorId,
    vector<edge_descriptor>& edges) const
{
    edges.clear();

    const auto annotations = getAnnotations(anchorId);
    for(const Annotation& annotation: annotations) {
        if(annotation.v == null_vertex()) {
            edges.push_back(annotation.e);
        }
    }
    deduplicate(edges);
}

