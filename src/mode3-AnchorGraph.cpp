// Dinara.
#include "Alignment.hpp"
#include "mode3-AnchorGraph.hpp"
#include "deduplicate.hpp"
#include "longestPath.hpp"
#include "Marker.hpp"
#include "MurmurHash2.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
#include "weightedShuffle.hpp"
using namespace dinara;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"
#include <limits>
#include <map>
#include <queue>
#include <unordered_set>



// Create the AnchorGraph and its vertices and edges given a vector of AnchorIds.
AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    span<const AnchorId> anchorIds,
    uint64_t minEdgeCoverage,
    const MemoryMapped::Vector<AlignmentData>* alignmentData,
    const MemoryMapped::VectorOfVectors<uint32_t, uint32_t>* alignmentTable) :
    anchorIds(anchorIds),
    anchorsPointer(&anchors),
    alignmentDataPointer(alignmentData),
    alignmentTablePointer(alignmentTable)
{
    buildPhaseDeletedPartnerTable();

    // Check that the AnchorIds are sorted and distinct.
    for(uint64_t i=1; i<anchorIds.size(); i++) {
        DINARA_ASSERT(anchorIds[i-1] < anchorIds[i]);
    }

    // Create the vertices.
    for(uint64_t localAnchorId=0; localAnchorId<anchorIds.size(); localAnchorId++) {
        AnchorGraphVertex vertex;
        vertex.localAnchorId = localAnchorId;
        const vertex_descriptor v = add_vertex(vertex, *this);
        vertexDescriptors.push_back(v);
    }

    // Create the edges.
    vector<AnchorId> children;
    vector<uint64_t> counts;
    for(uint64_t localAnchorId0=0; localAnchorId0<anchorIds.size(); localAnchorId0++) {
        const AnchorId anchorId0 = anchorIds[localAnchorId0];
        anchors.findChildren(anchorId0, children, counts, minEdgeCoverage);
        const uint64_t n = children.size();
        DINARA_ASSERT(n == counts.size());
        for(uint64_t i=0; i<n; i++) {
            const AnchorId anchorId1 = children[i];
            const uint64_t coverage = counts[i];
            const AnchorId localAnchorId1 = anchors.getLocalAnchorIdInComponent(anchorId1);

            AnchorPairInfo info;
            anchors.analyzeAnchorPair(anchorId0, anchorId1, info);

            addEdgeFromLocalAnchorIds(localAnchorId0, localAnchorId1, info, coverage);
        }
    }
}



// Create the AnchorGraph for ALL anchors (not per-component).
// In this case, localAnchorId == anchorId for every anchor.
AnchorGraph::AnchorGraph(
    const Anchors& anchors,
    uint64_t minEdgeCoverage,
    const MemoryMapped::Vector<AlignmentData>* alignmentData,
    const MemoryMapped::VectorOfVectors<uint32_t, uint32_t>* alignmentTable) :
    anchorsPointer(&anchors),
    alignmentDataPointer(alignmentData),
    alignmentTablePointer(alignmentTable)
{
    buildPhaseDeletedPartnerTable();

    const uint64_t anchorCount = anchors.size();

    // Build the full set of anchor IDs.
    ownedAnchorIds.resize(anchorCount);
    for(uint64_t i = 0; i < anchorCount; i++) {
        ownedAnchorIds[i] = i;
    }
    anchorIds = span<const AnchorId>(ownedAnchorIds);

    // Create the vertices.
    for(uint64_t localAnchorId = 0; localAnchorId < anchorCount; localAnchorId++) {
        AnchorGraphVertex vertex;
        vertex.localAnchorId = localAnchorId;
        const vertex_descriptor v = add_vertex(vertex, *this);
        vertexDescriptors.push_back(v);
    }

    // Create the edges.
    vector<AnchorId> children;
    vector<uint64_t> counts;
    for(uint64_t anchorId0 = 0; anchorId0 < anchorCount; anchorId0++) {
        anchors.findChildren(anchorId0, children, counts, minEdgeCoverage);
        const uint64_t n = children.size();
        DINARA_ASSERT(n == counts.size());
        for(uint64_t i = 0; i < n; i++) {
            const AnchorId anchorId1 = children[i];
            const uint64_t coverage = counts[i];

            AnchorPairInfo info;
            anchors.analyzeAnchorPair(anchorId0, anchorId1, info);

            addEdgeFromLocalAnchorIds(anchorId0, anchorId1, info, coverage);
        }
    }
}



void AnchorGraph::addEdgeFromLocalAnchorIds(
    uint64_t localAnchorId0,
    uint64_t localAnchorId1,
    const AnchorPairInfo& info,
    uint64_t coverage)
{
    boost::add_edge(
        vertexDescriptors[localAnchorId0],
        vertexDescriptors[localAnchorId1],
        AnchorGraphEdge(info, coverage), *this);
}



// Write a AnchorGraph in graphviz format.
void AnchorGraph::writeGraphviz(
    const string& name,
    const AnchorGraphDisplayOptions& options,
    const Anchors& anchors) const
{
    ofstream out(name + ".dot");

    const AnchorGraph& graph = *this;
    out << "digraph " << name << " {\n";

    BGL_FORALL_VERTICES(v, graph, AnchorGraph) {
        out << getAnchorId(v);

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "[";
        }

        if(options.labels) {
            out << "label=\"";
            out << getAnchorId(v) << "\\n" << anchors[getAnchorId(v)].coverage();
            out << "\" ";
        }

        if(options.tooltips) {
            out << "tooltip=\"";
            out << getAnchorId(v);
            out << "\" ";
        }

        if(options.labels or options.tooltips or options.colorVertices) {
            out << "]";
        }
        out << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, AnchorGraph) {
        const AnchorGraphEdge& edge = graph[e];
        if(not options.showNonTransitiveReductionEdges and edge.isNonTransitiveReductionEdge) {
            continue;
        }
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        out <<
            getAnchorId(v0) << "->" <<
            getAnchorId(v1);

        if(edge.isNonTransitiveReductionEdge or options.labels or options.tooltips or options.colorEdges) {
            out << " [";
        }

        if(edge.isNonTransitiveReductionEdge) {
            out << "style=dashed ";
        }

        if(options.tooltips) {
            out <<
                "tooltip=\"" <<
                getAnchorId(v0) << "->" <<
                getAnchorId(v1) << " ";
            if(edge.coverage != invalid<uint64_t>) {
                out << edge.coverage << "/";
            }
            out <<
                edge.info.common << " " <<
                std::fixed << std::setprecision(2) << edge.info.correctedJaccard() << " " <<
                edge.info.offsetInBases << "\" ";
        }

        if(options.labels) {
            out <<
                "label=\"";
            if(edge.coverage != invalid<uint64_t>) {
                out << edge.coverage << "/";
            }
            out <<
                edge.info.common << "\\n" <<
                std::fixed << std::setprecision(2) << edge.info.correctedJaccard() << "\\n" <<
                edge.info.offsetInBases << "\" ";

        }

        // Color.
        if(options.colorEdges) {
            const double correctedJaccard = edge.info.correctedJaccard();
            if(correctedJaccard <= options.redJ) {
                out << " color=red ";
            } else if(correctedJaccard >= options.greenJ) {
                out << " color=green ";
            } else {
                const double hue = (correctedJaccard - options.redJ) / (3. * (options.greenJ - options.redJ));
                out << " color=\"" << hue << ",1,1\" ";
            }
        }

        if(edge.isNonTransitiveReductionEdge or options.labels or options.tooltips or options.colorEdges) {
            out << "]";
        }
        out << ";\n";
    }

    out << "}\n";
}



void AnchorGraph::writeEdgeCoverageHistogram(const string& fileName) const
{
    const AnchorGraph& primaryGraph = *this;

    // Create a histogram indexed by histogram[coverage][commonCount].
    vector< vector<uint64_t> > histogram;

    // Loop over all edges.
    BGL_FORALL_EDGES(e, primaryGraph, AnchorGraph) {
        const AnchorGraphEdge& edge = primaryGraph[e];
        const uint64_t coverage = edge.coverage;
        const uint64_t commonCount = edge.info.common;
        DINARA_ASSERT(coverage <= commonCount);

        // Increment the histogram, making space as necessary.
        if(coverage >= histogram.size()) {
            histogram.resize(coverage + 1);
        }
        vector<uint64_t>& h = histogram[coverage];
        if(commonCount >= h.size()) {
            h.resize(commonCount + 1, 0);
        }
        ++h[commonCount];
    }

    // Write out the histogram.
    ofstream csv(fileName);
    csv << "Coverage,Common count,Loss,Frequency\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        const vector<uint64_t>& h = histogram[coverage];
        for(uint64_t commonCount=0; commonCount<h.size(); commonCount++) {
            const uint64_t frequency = h[commonCount];

            if(frequency > 0) {
                const uint64_t loss = commonCount - coverage;
                csv << coverage << ",";
                csv << commonCount << ",";
                csv << loss << ",";
                csv << frequency << "\n";
            }
        }
    }
}



void AnchorGraph::writeEdgeDetails(
    const string& fileName,
    const Anchors& anchors) const
{
    const AnchorGraph& anchorGraph = *this;

    ofstream csv(fileName);
    csv << "AnchorId0,AnchorId1,Coverage,Common,Common with positive offset,Offset,\n";


    // Loop over all edges.
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        const AnchorGraphEdge& edge = anchorGraph[e];
        const uint64_t coverage = edge.coverage;
        const uint64_t commonCount = edge.info.common;
        DINARA_ASSERT(coverage <= commonCount);

        const vertex_descriptor v0 = source(e, anchorGraph);
        const vertex_descriptor v1 = target(e, anchorGraph);
        const AnchorId anchorId0 = getAnchorId(v0);
        const AnchorId anchorId1 = getAnchorId(v1);

        const uint64_t commonCountPositiveOffset = anchors.countCommon(anchorId0, anchorId1, true);
        DINARA_ASSERT(commonCountPositiveOffset <= commonCount);

        csv << anchorIdToString(anchorId0) << ",";
        csv << anchorIdToString(anchorId1) << ",";
        csv << coverage << ",";
        csv << commonCount << ",";
        csv << commonCountPositiveOffset << ",";
        csv << edge.info.offsetInBases << ",";
        csv << "\n";

    }

}



void AnchorGraph::removeNegativeOffsetEdges()
{
    AnchorGraph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, AnchorGraph) {
        const AnchorGraphEdge& edge = graph[e];
        DINARA_ASSERT(edge.info.common > 0);
        if(edge.info.offsetInBases < 0) {
            edgesToBeRemoved.push_back(e);
        }
    }

    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



// Remove cross-edges.
// This removes an edge v0->v1 if the following are all true:
// - It is not marked as removed by transitive reduction.
// - Its coverage is at most lowCoverageThreshold.
// - Its estimated offset is at least minOffset.
// - v0 has at least one out-edge with coverage at least highCoverageThreshold
//   (ignoring edges marked as removed by transitive reduction).
// - v1 has at least one in-edge with coverage at least highCoverageThreshold.
//   (ignoring edges marked as removed by transitive reduction).
void AnchorGraph::removeCrossEdges(
    uint64_t lowCoverageThreshold,
    uint64_t highCoverageThreshold,
    uint64_t minOffset,
    bool debug)
{
    AnchorGraph& graph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, AnchorGraph) {
        const AnchorGraphEdge& edge = graph[e];

        // If it is marked as removed by transitive reduction, skip it.
        if(edge.isNonTransitiveReductionEdge) {
            continue;
        }

        // Check coverage.
        if(edge.coverage > lowCoverageThreshold) {
            continue;
        }

        // Check estimated offset.
        if(edge.info.offsetInBases < int64_t(minOffset)) {
            continue;
        }

        // Check out-edges of v0.
        const vertex_descriptor v0 = source(e, graph);
        bool v0HasStrongOutEdge = false;
        BGL_FORALL_OUTEDGES(v0, e0, graph, AnchorGraph) {
            // If it is marked as removed by transitive reduction, ignore it.
            if(graph[e0].isNonTransitiveReductionEdge) {
                continue;
            }
            if(graph[e0].coverage >= highCoverageThreshold) {
                v0HasStrongOutEdge = true;
                break;
            }
        }
        if(not v0HasStrongOutEdge) {
            continue;
        }

        // Check in-edges of v1.
        const vertex_descriptor v1 = target(e, graph);
        bool v1HasStrongOutEdge = false;
        BGL_FORALL_INEDGES(v1, e1, graph, AnchorGraph) {
            // If it is marked as removed by transitive reduction, ignore it.
            if(graph[e1].isNonTransitiveReductionEdge) {
                continue;
            }
            if(graph[e1].coverage >= highCoverageThreshold) {
                v1HasStrongOutEdge = true;
                break;
            }
        }
        if(not v1HasStrongOutEdge) {
            continue;
        }

        // If all above checks passed, this edge will be removed.
        edgesToBeRemoved.push_back(e);
        if(debug) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            cout << "Removing cross edge " <<
                getAnchorId(v0) << "->" <<
                getAnchorId(v1) << endl;
        }
    }

    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



// Remove edges for which loss = (commonCount - coverage) / commonCount > maxLoss
void AnchorGraph::removeWeakEdges(double maxLoss, bool debug)
{
    AnchorGraph& graph = *this;

    // Find the edges we are going to remove.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, AnchorGraph) {
        const AnchorGraphEdge& edge = graph[e];
        const double loss = double(edge.info.common - edge.coverage) / double(edge.info.common);
        if(loss > maxLoss) {
            edgesToBeRemoved.push_back(e);

            if(debug) {
                const vertex_descriptor v0 = source(e, graph);
                const vertex_descriptor v1 = target(e, graph);
                cout << "Removing weak edge " <<
                    getAnchorId(v0) << "->" <<
                    getAnchorId(v1) << ", loss " << loss << endl;
            }
        }
    }



    // Remove the edges we found.
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

}



void AnchorGraph::buildPhaseDeletedPartnerTable()
{
    phaseDeletedPartnersByRead.clear();
    hasPhaseDeletedPartnerCache = false;

    if(alignmentDataPointer == nullptr || alignmentTablePointer == nullptr) {
        return;
    }

    const auto& alignmentData = *alignmentDataPointer;
    const auto& alignmentTable = *alignmentTablePointer;
    if(!alignmentData.isOpen || !alignmentTable.isOpen()) {
        return;
    }

    phaseDeletedPartnersByRead.resize(alignmentTable.size());
    for(uint64_t orientedReadValue=0; orientedReadValue<alignmentTable.size(); orientedReadValue++) {
        auto& partners = phaseDeletedPartnersByRead[orientedReadValue];
        const OrientedReadId orientedReadId = OrientedReadId::fromValue(uint32_t(orientedReadValue));

        for(const uint32_t alignmentId: alignmentTable[orientedReadValue]) {
            if(alignmentId >= alignmentData.size()) {
                continue;
            }
            const AlignmentData& ad = alignmentData[alignmentId];

            AlignmentData::DeleteReasonMask reasons = AlignmentData::DeleteReasonNone;
            if(ad.readIds[0] == orientedReadId.getReadId()) {
                reasons = ad.deleteReasons0;
            } else if(ad.readIds[1] == orientedReadId.getReadId()) {
                reasons = ad.deleteReasons1;
            } else {
                continue;
            }

            if((reasons & AlignmentData::DeleteReasonPhase) == 0) {
                continue;
            }
            partners.push_back(ad.getOther(orientedReadId));
        }
        deduplicate(partners);
    }

    hasPhaseDeletedPartnerCache = true;
}



void AnchorGraph::getSupportingOrientedReads(
    edge_descriptor e,
    vector<OrientedReadId>& supportingReads) const
{
    supportingReads.clear();
    if(anchorsPointer == nullptr) {
        return;
    }

    const AnchorGraph& anchorGraph = *this;
    const auto v0 = source(e, anchorGraph);
    const auto v1 = target(e, anchorGraph);
    const AnchorId anchorId0 = getAnchorId(v0);
    const AnchorId anchorId1 = getAnchorId(v1);

    const auto markerIntervals = (*anchorsPointer)[anchorId0];
    for(const auto& markerInterval: markerIntervals) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto journey = anchorsPointer->journeys[orientedReadId.getValue()];
        const uint64_t position = markerInterval.positionInJourney;
        if(position == invalid<uint32_t> || position >= journey.size()) {
            continue;
        }
        const uint64_t nextPosition = position + 1;
        if(nextPosition < journey.size() && journey[nextPosition] == anchorId1) {
            supportingReads.push_back(orientedReadId);
        }
    }

    deduplicate(supportingReads);
}



bool AnchorGraph::hasPhaseDeletedConflict(
    const vector<OrientedReadId>& firstSet,
    const vector<OrientedReadId>& secondSet) const
{
    if(!hasPhaseDeletedPartnerCache || secondSet.empty()) {
        return false;
    }

    std::unordered_set<uint32_t> secondSetValues;
    secondSetValues.reserve(secondSet.size());
    for(const OrientedReadId orientedReadId: secondSet) {
        secondSetValues.insert(orientedReadId.getValue());
    }

    for(const OrientedReadId orientedReadId: firstSet) {
        const uint64_t orientedReadValue = orientedReadId.getValue();
        if(orientedReadValue >= phaseDeletedPartnersByRead.size()) {
            continue;
        }
        const auto& partners = phaseDeletedPartnersByRead[orientedReadValue];
        for(const OrientedReadId partnerOrientedReadId: partners) {
            if(secondSetValues.contains(partnerOrientedReadId.getValue())) {
                return true;
            }
        }
    }

    return false;
}



bool AnchorGraph::keepCandidateEdgeByPhasing(
    edge_descriptor candidateEdge,
    const vector<edge_descriptor>& alternatePathEdges) const
{
    if(!hasPhaseDeletedPartnerCache) {
        return false;
    }

    vector<OrientedReadId> candidateReads;
    getSupportingOrientedReads(candidateEdge, candidateReads);
    if(candidateReads.empty()) {
        return false;
    }

    vector<OrientedReadId> alternateReads;
    vector<OrientedReadId> pathEdgeReads;
    for(const edge_descriptor pathEdge: alternatePathEdges) {
        getSupportingOrientedReads(pathEdge, pathEdgeReads);
        copy(pathEdgeReads.begin(), pathEdgeReads.end(), back_inserter(alternateReads));
    }
    deduplicate(alternateReads);
    if(alternateReads.empty()) {
        return false;
    }

    vector<OrientedReadId> candidateOnly;
    vector<OrientedReadId> alternateOnly;
    std::set_difference(
        candidateReads.begin(), candidateReads.end(),
        alternateReads.begin(), alternateReads.end(),
        back_inserter(candidateOnly));
    std::set_difference(
        alternateReads.begin(), alternateReads.end(),
        candidateReads.begin(), candidateReads.end(),
        back_inserter(alternateOnly));

    if(candidateOnly.empty() || alternateOnly.empty()) {
        return false;
    }

    return
        hasPhaseDeletedConflict(candidateOnly, alternateOnly) ||
        hasPhaseDeletedConflict(alternateOnly, candidateOnly);
}



bool AnchorGraph::containsEdge(
    const vector<edge_descriptor>& edges,
    edge_descriptor edge) const
{
    for(const edge_descriptor e: edges) {
        if(e == edge) {
            return true;
        }
    }
    return false;
}



bool AnchorGraph::findAlternatePath(
    edge_descriptor candidateEdge,
    uint64_t maxDistance,
    vector<edge_descriptor>& alternatePathEdges,
    bool onlyHigherCoverage) const
{
    alternatePathEdges.clear();

    const AnchorGraph& anchorGraph = *this;
    const uint64_t candidateCoverage = anchorGraph[candidateEdge].coverage;
    const vertex_descriptor sourceVertex = source(candidateEdge, anchorGraph);
    const vertex_descriptor targetVertex = target(candidateEdge, anchorGraph);

    std::queue<vertex_descriptor> q;
    q.push(sourceVertex);

    std::map<vertex_descriptor, uint64_t> distance;
    distance.insert(make_pair(sourceVertex, 0));
    std::map<vertex_descriptor, edge_descriptor> parentEdge;

    while(not q.empty()) {
        const vertex_descriptor vA = q.front();
        q.pop();

        const auto itDistanceA = distance.find(vA);
        DINARA_ASSERT(itDistanceA != distance.end());
        const uint64_t distanceA = itDistanceA->second;
        const uint64_t distanceB = distanceA + 1;

        BGL_FORALL_OUTEDGES(vA, eAB, anchorGraph, AnchorGraph) {
            if(eAB == candidateEdge) {
                continue;
            }
            const AnchorGraphEdge& edgeAB = anchorGraph[eAB];
            if(edgeAB.isNonTransitiveReductionEdge) {
                continue;
            }
            if(onlyHigherCoverage && edgeAB.coverage <= candidateCoverage) {
                continue;
            }

            const vertex_descriptor vB = target(eAB, anchorGraph);
            if(vB == targetVertex) {
                alternatePathEdges.push_back(eAB);
                vertex_descriptor vTrace = vA;
                while(vTrace != sourceVertex) {
                    auto itParent = parentEdge.find(vTrace);
                    if(itParent == parentEdge.end()) {
                        alternatePathEdges.clear();
                        return false;
                    }
                    const edge_descriptor parent = itParent->second;
                    alternatePathEdges.push_back(parent);
                    vTrace = source(parent, anchorGraph);
                }
                reverse(alternatePathEdges.begin(), alternatePathEdges.end());
                return not alternatePathEdges.empty();
            }

            if(distance.contains(vB)) {
                continue;
            }
            if(distanceB < maxDistance) {
                q.push(vB);
                distance.insert(make_pair(vB, distanceB));
                parentEdge.insert(make_pair(vB, eAB));
            }
        }
    }

    return false;
}



void AnchorGraph::getIntervalVertices(
    vertex_descriptor sourceVertex,
    vertex_descriptor targetVertex,
    uint64_t maxDistance,
    vector<vertex_descriptor>& intervalVertices) const
{
    intervalVertices.clear();
    const AnchorGraph& anchorGraph = *this;

    std::queue<vertex_descriptor> qForward;
    qForward.push(sourceVertex);
    std::map<vertex_descriptor, uint64_t> distanceForward;
    distanceForward.insert(make_pair(sourceVertex, 0));
    while(not qForward.empty()) {
        const vertex_descriptor vA = qForward.front();
        qForward.pop();
        const uint64_t distanceA = distanceForward.find(vA)->second;
        const uint64_t distanceB = distanceA + 1;
        BGL_FORALL_OUTEDGES(vA, eAB, anchorGraph, AnchorGraph) {
            if(anchorGraph[eAB].isNonTransitiveReductionEdge) {
                continue;
            }
            const vertex_descriptor vB = target(eAB, anchorGraph);
            if(distanceForward.contains(vB)) {
                continue;
            }
            if(distanceB <= maxDistance) {
                distanceForward.insert(make_pair(vB, distanceB));
                qForward.push(vB);
            }
        }
    }

    std::queue<vertex_descriptor> qBackward;
    qBackward.push(targetVertex);
    std::map<vertex_descriptor, uint64_t> distanceBackward;
    distanceBackward.insert(make_pair(targetVertex, 0));
    while(not qBackward.empty()) {
        const vertex_descriptor vA = qBackward.front();
        qBackward.pop();
        const uint64_t distanceA = distanceBackward.find(vA)->second;
        const uint64_t distanceB = distanceA + 1;
        BGL_FORALL_INEDGES(vA, eBA, anchorGraph, AnchorGraph) {
            if(anchorGraph[eBA].isNonTransitiveReductionEdge) {
                continue;
            }
            const vertex_descriptor vB = source(eBA, anchorGraph);
            if(distanceBackward.contains(vB)) {
                continue;
            }
            if(distanceB <= maxDistance) {
                distanceBackward.insert(make_pair(vB, distanceB));
                qBackward.push(vB);
            }
        }
    }

    for(const auto& p: distanceForward) {
        const vertex_descriptor v = p.first;
        const uint64_t dForward = p.second;
        const auto itBackward = distanceBackward.find(v);
        if(itBackward == distanceBackward.end()) {
            continue;
        }
        const uint64_t dBackward = itBackward->second;
        if(dForward + dBackward <= maxDistance) {
            intervalVertices.push_back(v);
        }
    }
}



bool AnchorGraph::isAcyclicInterval(
    const vector<vertex_descriptor>& intervalVertices,
    vector<vertex_descriptor>& topologicalOrder) const
{
    topologicalOrder.clear();
    if(intervalVertices.empty()) {
        return true;
    }

    const AnchorGraph& anchorGraph = *this;
    std::map<vertex_descriptor, uint64_t> vertexIndex;
    for(uint64_t i=0; i<intervalVertices.size(); i++) {
        vertexIndex.insert(make_pair(intervalVertices[i], i));
    }

    vector<uint64_t> inDegree(intervalVertices.size(), 0);
    for(const vertex_descriptor v: intervalVertices) {
        BGL_FORALL_OUTEDGES(v, eAB, anchorGraph, AnchorGraph) {
            if(anchorGraph[eAB].isNonTransitiveReductionEdge) {
                continue;
            }
            const vertex_descriptor vB = target(eAB, anchorGraph);
            if(not vertexIndex.contains(vB)) {
                continue;
            }
            ++inDegree[vertexIndex.find(vB)->second];
        }
    }

    std::queue<vertex_descriptor> q;
    for(const auto& p: vertexIndex) {
        if(inDegree[p.second] == 0) {
            q.push(p.first);
        }
    }

    while(not q.empty()) {
        const vertex_descriptor vA = q.front();
        q.pop();
        topologicalOrder.push_back(vA);

        BGL_FORALL_OUTEDGES(vA, eAB, anchorGraph, AnchorGraph) {
            if(anchorGraph[eAB].isNonTransitiveReductionEdge) {
                continue;
            }
            const vertex_descriptor vB = target(eAB, anchorGraph);
            const auto itB = vertexIndex.find(vB);
            if(itB == vertexIndex.end()) {
                continue;
            }
            const uint64_t iB = itB->second;
            DINARA_ASSERT(inDegree[iB] > 0);
            --inDegree[iB];
            if(inDegree[iB] == 0) {
                q.push(vB);
            }
        }
    }

    return topologicalOrder.size() == intervalVertices.size();
}



vector<AnchorGraph::edge_descriptor> AnchorGraph::computeBestHaplotypePathInInterval(
    vertex_descriptor sourceVertex,
    vertex_descriptor targetVertex,
    const vector<vertex_descriptor>& topologicalOrder,
    const vector<edge_descriptor>& intervalEdges,
    const std::map<vertex_descriptor, uint64_t>& vertexIndex,
    const vector<array<int64_t, 2> >& edgeHaplotypeScores,
    uint64_t haplotype,
    const vector<bool>& disallowedInternalVertex) const
{
    DINARA_ASSERT(haplotype < 2);
    const AnchorGraph& anchorGraph = *this;
    const int64_t minusInfinity = std::numeric_limits<int64_t>::min() / 4;

    vector<int64_t> bestScore(vertexIndex.size(), minusInfinity);
    vector<edge_descriptor> predecessor(vertexIndex.size());
    vector<bool> hasPredecessor(vertexIndex.size(), false);

    const uint64_t sourceIndex = vertexIndex.find(sourceVertex)->second;
    const uint64_t targetIndex = vertexIndex.find(targetVertex)->second;
    bestScore[sourceIndex] = 0;

    for(const vertex_descriptor vA: topologicalOrder) {
        const uint64_t iA = vertexIndex.find(vA)->second;
        if(bestScore[iA] == minusInfinity) {
            continue;
        }
        if(vA != sourceVertex && vA != targetVertex && disallowedInternalVertex[iA]) {
            continue;
        }

        BGL_FORALL_OUTEDGES(vA, eAB, anchorGraph, AnchorGraph) {
            if(anchorGraph[eAB].isNonTransitiveReductionEdge) {
                continue;
            }
            if(not containsEdge(intervalEdges, eAB)) {
                continue;
            }

            const vertex_descriptor vB = target(eAB, anchorGraph);
            const auto itB = vertexIndex.find(vB);
            if(itB == vertexIndex.end()) {
                continue;
            }
            const uint64_t iB = itB->second;
            if(vB != sourceVertex && vB != targetVertex && disallowedInternalVertex[iB]) {
                continue;
            }

            uint64_t edgeIndex = invalid<uint64_t>;
            for(uint64_t i=0; i<intervalEdges.size(); i++) {
                if(intervalEdges[i] == eAB) {
                    edgeIndex = i;
                    break;
                }
            }
            if(edgeIndex == invalid<uint64_t>) {
                continue;
            }

            const int64_t weight =
                edgeHaplotypeScores[edgeIndex][haplotype] -
                edgeHaplotypeScores[edgeIndex][1 - haplotype];
            const int64_t candidateScore = bestScore[iA] + weight;
            if(candidateScore > bestScore[iB]) {
                bestScore[iB] = candidateScore;
                predecessor[iB] = eAB;
                hasPredecessor[iB] = true;
            }
        }
    }

    if(sourceVertex != targetVertex && !hasPredecessor[targetIndex]) {
        return {};
    }

    vector<edge_descriptor> path;
    vertex_descriptor v = targetVertex;
    while(v != sourceVertex) {
        const uint64_t i = vertexIndex.find(v)->second;
        if(!hasPredecessor[i]) {
            return {};
        }
        const edge_descriptor e = predecessor[i];
        path.push_back(e);
        v = source(e, anchorGraph);
    }
    reverse(path.begin(), path.end());
    return path;
}



void AnchorGraph::resolvePhasedSuperbubbles(uint64_t maxDistance, bool debug)
{
    AnchorGraph& anchorGraph = *this;
    if(!hasPhaseDeletedPartnerCache) {
        if(debug) {
            cout << "resolvePhasedSuperbubbles skipped: no phased partner cache." << endl;
        }
        return;
    }

    vector<edge_descriptor> edgesToInspect;
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        if(!anchorGraph[e].isNonTransitiveReductionEdge) {
            edgesToInspect.push_back(e);
        }
    }

    uint64_t candidateCount = 0;
    uint64_t resolvedCount = 0;

    for(const edge_descriptor candidateEdge: edgesToInspect) {
        if(anchorGraph[candidateEdge].isNonTransitiveReductionEdge) {
            continue;
        }

        vector<edge_descriptor> alternatePathEdges;
        if(!findAlternatePath(candidateEdge, maxDistance, alternatePathEdges, false)) {
            continue;
        }
        if(!keepCandidateEdgeByPhasing(candidateEdge, alternatePathEdges)) {
            continue;
        }
        ++candidateCount;

        const vertex_descriptor sourceVertex = source(candidateEdge, anchorGraph);
        const vertex_descriptor targetVertex = target(candidateEdge, anchorGraph);

        vector<vertex_descriptor> intervalVertices;
        getIntervalVertices(sourceVertex, targetVertex, maxDistance, intervalVertices);
        if(intervalVertices.size() < 2) {
            continue;
        }

        std::map<vertex_descriptor, uint64_t> vertexIndex;
        for(uint64_t i=0; i<intervalVertices.size(); i++) {
            vertexIndex.insert(make_pair(intervalVertices[i], i));
        }
        if(!vertexIndex.contains(sourceVertex) || !vertexIndex.contains(targetVertex)) {
            continue;
        }

        vector<edge_descriptor> intervalEdges;
        for(const vertex_descriptor vA: intervalVertices) {
            BGL_FORALL_OUTEDGES(vA, eAB, anchorGraph, AnchorGraph) {
                if(anchorGraph[eAB].isNonTransitiveReductionEdge) {
                    continue;
                }
                const vertex_descriptor vB = target(eAB, anchorGraph);
                if(vertexIndex.contains(vB)) {
                    intervalEdges.push_back(eAB);
                }
            }
        }
        if(intervalEdges.empty()) {
            continue;
        }

        vector<edge_descriptor> pathHaplotype0;
        vector<edge_descriptor> pathHaplotype1;

        vector<vertex_descriptor> topologicalOrder;
        const bool isAcyclic = isAcyclicInterval(intervalVertices, topologicalOrder);
        if(!isAcyclic) {
            pathHaplotype0 = {candidateEdge};
            pathHaplotype1 = alternatePathEdges;
        } else {
            vector<vector<OrientedReadId> > edgeReads(intervalEdges.size());
            vector<OrientedReadId> allReads;
            for(uint64_t i=0; i<intervalEdges.size(); i++) {
                getSupportingOrientedReads(intervalEdges[i], edgeReads[i]);
                copy(edgeReads[i].begin(), edgeReads[i].end(), back_inserter(allReads));
            }
            deduplicate(allReads);
            if(allReads.empty()) {
                continue;
            }

            std::map<uint32_t, uint64_t> readIndex;
            for(uint64_t i=0; i<allReads.size(); i++) {
                readIndex.insert(make_pair(allReads[i].getValue(), i));
            }

            vector<vector<uint64_t> > conflicts(allReads.size());
            bool hasConflictEvidence = false;
            for(uint64_t i=0; i<allReads.size(); i++) {
                const uint64_t orientedReadValue = allReads[i].getValue();
                if(orientedReadValue >= phaseDeletedPartnersByRead.size()) {
                    continue;
                }
                for(const OrientedReadId partner: phaseDeletedPartnersByRead[orientedReadValue]) {
                    const auto it = readIndex.find(partner.getValue());
                    if(it == readIndex.end()) {
                        continue;
                    }
                    const uint64_t j = it->second;
                    if(i == j) {
                        continue;
                    }
                    conflicts[i].push_back(j);
                    hasConflictEvidence = true;
                }
                deduplicate(conflicts[i]);
            }
            if(!hasConflictEvidence) {
                continue;
            }

            vector<int8_t> readHaplotype(allReads.size(), int8_t(-1));
            std::queue<uint64_t> q;
            for(uint64_t i=0; i<allReads.size(); i++) {
                if(readHaplotype[i] != -1) {
                    continue;
                }
                readHaplotype[i] = 0;
                q.push(i);
                while(not q.empty()) {
                    const uint64_t v = q.front();
                    q.pop();
                    for(const uint64_t w: conflicts[v]) {
                        const int8_t required = int8_t(1 - readHaplotype[v]);
                        if(readHaplotype[w] == -1) {
                            readHaplotype[w] = required;
                            q.push(w);
                        }
                    }
                }
            }

            vector<array<int64_t, 2> > edgeHaplotypeScores(intervalEdges.size(), {0, 0});
            for(uint64_t i=0; i<intervalEdges.size(); i++) {
                for(const OrientedReadId orientedReadId: edgeReads[i]) {
                    const auto it = readIndex.find(orientedReadId.getValue());
                    if(it == readIndex.end()) {
                        continue;
                    }
                    const int8_t hap = readHaplotype[it->second];
                    if(hap == 0) {
                        ++edgeHaplotypeScores[i][0];
                    } else if(hap == 1) {
                        ++edgeHaplotypeScores[i][1];
                    }
                }
            }

            vector<bool> disallowedInternalVertex(intervalVertices.size(), false);
            pathHaplotype0 = computeBestHaplotypePathInInterval(
                sourceVertex,
                targetVertex,
                topologicalOrder,
                intervalEdges,
                vertexIndex,
                edgeHaplotypeScores,
                0,
                disallowedInternalVertex);
            if(pathHaplotype0.empty()) {
                pathHaplotype0 = {candidateEdge};
            }

            for(const edge_descriptor e: pathHaplotype0) {
                const vertex_descriptor v0 = source(e, anchorGraph);
                const vertex_descriptor v1 = target(e, anchorGraph);
                if(v0 != sourceVertex && v0 != targetVertex) {
                    disallowedInternalVertex[vertexIndex.find(v0)->second] = true;
                }
                if(v1 != sourceVertex && v1 != targetVertex) {
                    disallowedInternalVertex[vertexIndex.find(v1)->second] = true;
                }
            }

            pathHaplotype1 = computeBestHaplotypePathInInterval(
                sourceVertex,
                targetVertex,
                topologicalOrder,
                intervalEdges,
                vertexIndex,
                edgeHaplotypeScores,
                1,
                disallowedInternalVertex);
            if(pathHaplotype1.empty()) {
                pathHaplotype1 = alternatePathEdges;
            }
            if(pathHaplotype1.empty()) {
                continue;
            }
        }

        vector<edge_descriptor> keptEdges = pathHaplotype0;
        copy(pathHaplotype1.begin(), pathHaplotype1.end(), back_inserter(keptEdges));

        for(const edge_descriptor e: intervalEdges) {
            anchorGraph[e].isNonTransitiveReductionEdge = !containsEdge(keptEdges, e);
        }
        ++resolvedCount;

        if(debug) {
            cout << "Resolved phased superbubble "
                << getAnchorId(sourceVertex) << "->" << getAnchorId(targetVertex)
                << " intervalVertices=" << intervalVertices.size()
                << " intervalEdges=" << intervalEdges.size()
                << " keep0=" << pathHaplotype0.size()
                << " keep1=" << pathHaplotype1.size() << endl;
        }
    }

    cout << "resolvePhasedSuperbubbles: " << resolvedCount <<
        " resolved out of " << candidateCount << " phased candidates." << endl;
}



// Transitive reduction.
// Processes edges in order of increasing coverage.
// For each edge at a given coverage level, does BFS from source to target
// using only edges with higher coverage that are not already flagged,
// within maxDistance hops. If an alternate path exists, the edge is flagged
// as isNonTransitiveReductionEdge.
void AnchorGraph::transitiveReduction(
    uint64_t maxEdgeCoverage,
    uint64_t maxDistance)
{
    AnchorGraph& anchorGraph = *this;
    cout << "AnchorGraph transitive reduction begins." << endl;

    // Initially make sure all edges are not flagged.
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        anchorGraph[e].isNonTransitiveReductionEdge = false;
    }

    // Loop over edge coverage.
    // At each iteration we only consider edges with this coverage.
    vector<edge_descriptor> edgesToProcess;
    vector<edge_descriptor> edgesToRemove;
    for(uint64_t edgeCoverage = 1; edgeCoverage <= maxEdgeCoverage; edgeCoverage++) {

        // Gather edges with this coverage.
        edgesToProcess.clear();
        BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
            if(anchorGraph[e].coverage == edgeCoverage) {
                edgesToProcess.push_back(e);
            }
        }

        // If there are none, there is nothing to do.
        if(edgesToProcess.empty()) {
            continue;
        }

        // Loop over all edges with this coverage.
        edgesToRemove.clear();
        for(const edge_descriptor e: edgesToProcess) {
            if(transitiveReductionCanRemove(e, maxDistance)) {
                edgesToRemove.push_back(e);
            }
        }

        // Flag edges removed at this iteration.
        for(const edge_descriptor e: edgesToRemove) {
            anchorGraph[e].isNonTransitiveReductionEdge = true;
        }
        cout << "Edge coverage " << edgeCoverage <<
            ": processed " << edgesToProcess.size() <<
            " edges and flagged " << edgesToRemove.size() << endl;
    }
    cout << "AnchorGraph transitive reduction ends." << endl;

    uint64_t keptCount = 0;
    BGL_FORALL_EDGES(e, anchorGraph, AnchorGraph) {
        if(not anchorGraph[e].isNonTransitiveReductionEdge) {
            ++keptCount;
        }
    }
    cout << keptCount << " edges kept out of " <<
        num_edges(anchorGraph) << " total." << endl;
}



// BFS helper for transitive reduction.
// Returns true if an alternate path exists from source(e) to target(e)
// using only edges with higher coverage that are not already flagged,
// within maxDistance hops.
bool AnchorGraph::transitiveReductionCanRemove(
    edge_descriptor e,
    uint64_t maxDistance) const
{
    const AnchorGraph& anchorGraph = *this;
    const uint64_t edgeCoverage = anchorGraph[e].coverage;

    const vertex_descriptor v0 = source(e, anchorGraph);
    const vertex_descriptor v1 = target(e, anchorGraph);

    // Do a forward BFS starting at v0, using edges
    // not flagged as isNonTransitiveReductionEdge
    // with coverage greater than edgeCoverage
    // and with maximum distance (number of edges)
    // equal to maxDistance.
    // If we encounter v1, the edge is removable unless phasing says to keep it.
    std::queue<vertex_descriptor> q;
    q.push(v0);

    // A map to store vertices already encountered and their distance from v0.
    std::map<vertex_descriptor, uint64_t> m;
    m.insert(make_pair(v0, 0));
    std::map<vertex_descriptor, edge_descriptor> parentEdge;

    bool foundAlternatePath = false;

    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor vA = q.front();
        q.pop();
        const auto itA = m.find(vA);
        DINARA_ASSERT(itA != m.end());
        const uint64_t distanceA = itA->second;
        const uint64_t distanceB = distanceA + 1;

        // Loop over its out-edges that are not flagged
        // and have sufficient coverage.
        BGL_FORALL_OUTEDGES(vA, eAB, anchorGraph, AnchorGraph) {
            const AnchorGraphEdge& edgeAB = anchorGraph[eAB];
            if(edgeAB.isNonTransitiveReductionEdge) {
                continue;
            }

            // Only use edges with higher coverage for the BFS.
            if(edgeAB.coverage <= edgeCoverage) {
                continue;
            }

            // If we reached v1, check if phased read evidence says this is a real bubble.
            const vertex_descriptor vB = target(eAB, anchorGraph);
            if(vB == v1) {
                foundAlternatePath = true;

                vector<edge_descriptor> alternatePathEdges;
                alternatePathEdges.push_back(eAB);
                vertex_descriptor vTrace = vA;
                while(vTrace != v0) {
                    auto itParent = parentEdge.find(vTrace);
                    if(itParent == parentEdge.end()) {
                        alternatePathEdges.clear();
                        break;
                    }
                    const edge_descriptor parent = itParent->second;
                    alternatePathEdges.push_back(parent);
                    vTrace = source(parent, anchorGraph);
                }
                reverse(alternatePathEdges.begin(), alternatePathEdges.end());

                if(!alternatePathEdges.empty() &&
                    keepCandidateEdgeByPhasing(e, alternatePathEdges)) {
                    return false;
                }
                continue;
            }

            // If we already encountered vB, don't do anything.
            if(m.contains(vB)) {
                continue;
            }

            if(distanceB < maxDistance) {
                q.push(vB);
                m.insert(make_pair(vB, distanceB));
                parentEdge.insert(make_pair(vB, eAB));
            }
        }
    }

    return foundAlternatePath;
}
