// Shasta.
#include "Shasta2RestrictedAnchorGraph.hpp"
#include "shasta2/algorithm.hpp"
#include "shasta2/stdexcept.hpp"
#include "shasta2/approximateTopologicalSort.hpp"
#include "deduplicate.hpp"
#include "dominatorTree.hpp"
#include "shasta2/findReachableVertices.hpp"
#include "Shasta2Journeys.hpp"
#include "shasta2/longestPath.hpp"
#include "Marker.hpp"
#include "Shasta2TangleMatrix1.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <filesystem>

// Standard library.
#include "fstream.hpp"



Shasta2RestrictedAnchorGraph::Shasta2RestrictedAnchorGraph(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const Shasta2TangleMatrix1& tangleMatrix1,
    uint64_t iEntrance,
    uint64_t iExit,
    ostream& html)
{
    constructFromTangleMatrix1(anchors, journeys, tangleMatrix1, iEntrance, iExit, html);
}



void Shasta2RestrictedAnchorGraph::constructFromTangleMatrix1(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const Shasta2TangleMatrix1& tangleMatrix1,
    uint64_t iEntrance,
    uint64_t iExit,
    ostream& html)
{
    using Graph = Shasta2RestrictedAnchorGraph;
    Graph& graph = *this;

    fillJourneyPortions(journeys, tangleMatrix1, iEntrance, iExit, html);
    gatherAllAnchorIds(journeys);
    fillJourneyPortionsAnchorIndexes(journeys);
    gatherTransitions(html);

    // Create the CycleAvoider.
    cycleAvoider = new shasta2::CycleAvoider<Shasta2RestrictedAnchorGraph>(*this);


    // A vector<bool>, indexed by anchorIndex, to contain anchorIndexes
    // for vertices we already added.
    vector<bool> verticesAdded(anchorMap.size(), false);

    // Add the entrance and exit vertices.
    const Shasta2AssemblyGraph::edge_descriptor entrance = tangleMatrix1.entrances[iEntrance];
    const Shasta2AssemblyGraph::edge_descriptor exit = tangleMatrix1.exits[iExit];
    const Shasta2AnchorId anchorId0 = tangleMatrix1.assemblyGraph[entrance].back().anchorPair.anchorIdB;
    const Shasta2AnchorId anchorId1 = tangleMatrix1.assemblyGraph[exit].front().anchorPair.anchorIdA;
    DINARA_ASSERT(anchorId0 != anchorId1);
    const vertex_descriptor v0 = addVertex(anchorId0);
    const vertex_descriptor v1 = addVertex(anchorId1);
    const uint64_t anchorIndex0 = getAnchorIndex(anchorId0);
    const uint64_t anchorIndex1 = getAnchorIndex(anchorId1);
    verticesAdded[anchorIndex0] = true;
    verticesAdded[anchorIndex1] = true;

    /*
    cout << tangleMatrix1.assemblyGraph[entrance].id << " " << tangleMatrix1.assemblyGraph[exit].id <<
        " " << allAnchorIds.size() << endl;
    */
    if(html) {
        html << "<br>Total number of distinct AnchorIds on the JourneyPortions is " << anchorMap.size();
    }

    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(anchorMap.size());
    vector<uint64_t> parent(anchorMap.size());
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<anchorMap.size(); i++) {
        disjointSets.make_set(i);
    }



    // Add edges in order of decreasing coverage.
    vector<uint64_t> anchorIndexesToAdd;
    for(uint64_t coverage=transitions.size()-1; coverage>0; coverage--) {
        const vector<Transition> transitionsToAdd = transitions[coverage];

        // Gather the anchorIndexes for the vertices that we need to add.
        anchorIndexesToAdd.clear();
        for(const Transition& transition: transitionsToAdd) {
            const uint64_t anchorIndex0 = transition.anchorIndex0;
            if(not verticesAdded[anchorIndex0]) {
                anchorIndexesToAdd.push_back(anchorIndex0);
                verticesAdded[anchorIndex0] = true;
            }
            const uint64_t anchorIndex1 = transition.anchorIndex1;
            if(not verticesAdded[anchorIndex1]) {
                anchorIndexesToAdd.push_back(anchorIndex1);
                verticesAdded[anchorIndex1] = true;
            }
        }

        // Add these vertices.
        for(const uint64_t anchorIndex: anchorIndexesToAdd) {
            const Shasta2AnchorId anchorId = anchorIds[anchorIndex];
            addVertex(anchorId);
            // cout << "Added vertex with Shasta2AnchorId " << anchorId << endl;
        }

        // Now we can add an edge for each of these transitions.
        for(const Transition& transition: transitionsToAdd) {
            const uint64_t anchorIndex0 = transition.anchorIndex0;
            const uint64_t anchorIndex1 = transition.anchorIndex1;
            const Shasta2AnchorId anchorId0 = anchorIds[anchorIndex0];
            const Shasta2AnchorId anchorId1 = anchorIds[anchorIndex1];
            // cout << "Possible edge between AnchorIds " << anchorId0 << " " << anchorId1 << endl;
            const vertex_descriptor v0 = getExistingVertex(anchorId0);
            const vertex_descriptor v1 = getExistingVertex(anchorId1);
            DINARA_ASSERT(v0 != null_vertex());
            DINARA_ASSERT(v1 != null_vertex());

            edge_descriptor e;
            bool edgeWasAdded = false;
            tie(e, edgeWasAdded) = cycleAvoider->addEdge(v0, v1);
            if(edgeWasAdded) {
                Shasta2RestrictedAnchorGraphEdge& edge = graph[e];
                edge.anchorPair.anchorIdA = anchorId0;
                edge.anchorPair.anchorIdB = anchorId1;
                disjointSets.union_set(anchorIndex0, anchorIndex1);
            }
        }

        // If the entrance and exit are different connected components,
        // continue the loop over coverage.
        if(disjointSets.find_set(anchorIndex0) != disjointSets.find_set(anchorIndex1)) {
            continue;
        }

        // If the exit is reachable from the entrance, stop the loop over coverage.
        if(shasta2::isReachable(graph, v0, v1, 0)) {
            break;
        }

        // If still not reachable and coverage is 1, we can't reduce coverage further.
        if(coverage == 1) {
            delete cycleAvoider;
            cycleAvoider = 0;
            throw Unreachable();
        }
    }



    // Destroy the CycleAvoider.
    delete cycleAvoider;
    cycleAvoider = 0;



    // Now fill in the OrientedReadIds of each vertex and edge.
    for(uint64_t j=0; j<journeyPortions.size(); j++) {
        const Shasta2JourneyPortion& journeyPortion = journeyPortions[j];
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;
        const vector<uint64_t>& journeyPortionAnchorIndexes = journeyPortionsAnchorIndexes[j];

        // Vertices.
        for(const uint64_t anchorIndex: journeyPortionAnchorIndexes) {
            const uint64_t anchorId = anchorIds[anchorIndex];
            const vertex_descriptor v = getVertex(anchorId);
            if(v != null_vertex()) {
                graph[v].orientedReadIds.push_back(orientedReadId);
            }
        }

        // Edges.
        for(uint64_t i1=1; i1<journeyPortionAnchorIndexes.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const uint64_t anchorIndex0 = journeyPortionAnchorIndexes[i0];
            const uint64_t anchorIndex1 = journeyPortionAnchorIndexes[i1];
            const uint64_t anchorId0 = anchorIds[anchorIndex0];
            const vertex_descriptor v0 = getVertex(anchorId0);
            if(v0 == null_vertex()) {
                continue;
            }
            const uint64_t anchorId1 = anchorIds[anchorIndex1];
            const vertex_descriptor v1 = getVertex(anchorId1);
            if(v1 == null_vertex()) {
                continue;
            }
            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(v0, v1, graph);
            if(edgeExists) {
                graph[e].anchorPair.orientedReadIds.push_back(orientedReadId);
            }
        }
    }

    // Fill in the offsets.
    BGL_FORALL_EDGES(e, graph, Graph) {
        Shasta2RestrictedAnchorGraphEdge& edge = graph[e];
        edge.offset = edge.anchorPair.getAverageOffset(anchors);
    }

    // Remove unreachable portions.
    keepBetween(anchorId0, anchorId1);
}



// Gather all the distinct AnchorIds that appear in the JourneyPortions
// and store them in the anchorMap.
void Shasta2RestrictedAnchorGraph::gatherAllAnchorIds(const Shasta2Journeys& journeys)
{
    // Initialize the anchorMap.
    uint64_t totalJourneyLength = 0;
    for(const Shasta2JourneyPortion&  journeyPortion:journeyPortions) {
        totalJourneyLength += journeyPortion.end - journeyPortion.begin;
    }
    anchorMap.initialize(totalJourneyLength);


    uint64_t nextAnchorIndex = 0;
    for(const Shasta2JourneyPortion&  journeyPortion:journeyPortions) {
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;
        const Shasta2Journey journey = journeys[orientedReadId];
        const uint32_t begin = journeyPortion.begin;
        const uint32_t end = journeyPortion.end;
        for(uint32_t i=begin; i<end; i++) {
            const Shasta2AnchorId anchorId = journey[i];
            auto p = anchorMap.insertNewOrGetExisting(make_pair(anchorId, AnchorInformation()));
            AnchorInformation& anchorInformation = p->second;
            if(anchorInformation.anchorIndex == invalid<uint64_t>) {
                anchorInformation.anchorIndex = nextAnchorIndex++;
                anchorIds.push_back(anchorId);
                /*
                cout << "Added Shasta2AnchorId " << anchorId <<
                    " with anchorIndex " << anchorInformation.anchorIndex << endl;
                */
            }
        }
    }
}



void Shasta2RestrictedAnchorGraph::fillJourneyPortionsAnchorIndexes(const Shasta2Journeys& journeys)
{
    const uint64_t orientedReadCount = journeyPortions.size();
    journeyPortionsAnchorIndexes.clear();
    journeyPortionsAnchorIndexes.resize(orientedReadCount);

    // Loop over the oriented reads.
    for(uint64_t j=0; j<orientedReadCount; j++) {
        const Shasta2JourneyPortion& journeyPortion = journeyPortions[j];
        vector<uint64_t>& journeyPortionAnchorIndexes =journeyPortionsAnchorIndexes[j];

        // Loop over the journey portion of this oriented read.
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;
        const Shasta2Journey journey = journeys[orientedReadId];
        const uint32_t begin = journeyPortion.begin;
        const uint32_t end = journeyPortion.end;
        journeyPortionAnchorIndexes.reserve(end - begin);
        for(uint32_t i=begin; i<end; i++) {
            const Shasta2AnchorId anchorId = journey[i];
            journeyPortionAnchorIndexes.push_back(getAnchorIndex(anchorId));
        }
    }
}



// Gather all transitions(anchorIndex0, anchorIndex1) for consecutive
// anchors in the journey portions. The number of times each
// transition appears in the journeys is its coverage.
// Store the transitions in a vector indexed by coverage.
void Shasta2RestrictedAnchorGraph::gatherTransitions(ostream& html)
{
    // First, gather all transitions and store the anchorIndexes1
    // in a vector indexed by anchorIndex0.
    vector< vector<uint64_t> > anchorIndexes1(anchorMap.size());
    for(const vector<uint64_t>& journeyPortionAnchorIndexes: journeyPortionsAnchorIndexes) {
        for(uint64_t i1=1; i1<journeyPortionAnchorIndexes.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const uint64_t anchorIndex0 = journeyPortionAnchorIndexes[i0];
            const uint64_t anchorIndex1 = journeyPortionAnchorIndexes[i1];
            DINARA_ASSERT(anchorIndex0 < anchorMap.size());
            DINARA_ASSERT(anchorIndex1 < anchorMap.size());
            anchorIndexes1[anchorIndex0].push_back(anchorIndex1);
        }
    }



    // Now for each anchorIndex0 deduplicate and count the anchorIndexes1.
    // Store its transitions according to its coverage.
    transitions.clear();
    vector<uint64_t> coverage;
    for(uint64_t anchorIndex0=0; anchorIndex0<anchorMap.size(); anchorIndex0++) {
        vector<uint64_t>& v = anchorIndexes1[anchorIndex0];
        deduplicateAndCount(v, coverage);

        for(uint64_t i=0; i<v.size(); i++) {
            const uint64_t c = coverage[i];
            if(c >= transitions.size()) {
                transitions.resize(c + 1);
            }
            transitions[c].push_back(Transition({anchorIndex0, v[i]}));
        }
    }

    if(html) {
        html << "<h3>Number of journey transitions by coverage</h3>"
            "<table><tr><th>Coverage<th>Frequency";
        for(uint64_t coverage=0; coverage<transitions.size(); coverage++) {
            const uint64_t frequency = transitions[coverage].size();
            if(frequency) {
                html << "<tr><td class=centered>" << coverage <<
                    "<td class=centered>" << frequency;
            }
        }
        html << "</table>";
    }

    if(transitions.empty()) {
        throw NoTransitions();
    }
}





// Fill the journey portions using a Shasta2TangleMatrix1.
void Shasta2RestrictedAnchorGraph::fillJourneyPortions(
    const Shasta2Journeys& journeys,
    const Shasta2TangleMatrix1& tangleMatrix1,
    uint64_t iEntrance,
    uint64_t iExit,
    ostream& html)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double drift = 0.2;

    using OrientedReadInfo = Shasta2TangleMatrix1::OrientedReadInfo;

    DINARA_ASSERT(tangleMatrix1.assemblyGraph.getAnchorsPointer() != nullptr);
    const Shasta2Anchors& anchors = *tangleMatrix1.assemblyGraph.getAnchorsPointer();
    const auto& markers = anchors.markers;

    const vector<OrientedReadInfo>& entranceOrientedReadInfos = tangleMatrix1.entranceOrientedReadInfos[iEntrance];
    const vector<OrientedReadInfo>& exitOrientedReadInfos = tangleMatrix1.exitOrientedReadInfos[iExit];



    // Loop over the common oriented reads to estimate an offset.
    auto itEntrance = entranceOrientedReadInfos.begin();
    auto itExit = exitOrientedReadInfos.begin();
    const auto itEntranceEnd = entranceOrientedReadInfos.end();
    const auto itExitEnd = exitOrientedReadInfos.end();
    uint64_t offsetSum = 0;
    uint64_t offsetCount = 0;
    while((itEntrance!=itEntranceEnd) and (itExit!=itExitEnd)) {
        if(itEntrance->orientedReadId < itExit->orientedReadId) {
            ++itEntrance;
            continue;
        }
        if(itExit->orientedReadId < itEntrance->orientedReadId) {
            ++itExit;
            continue;
        }
        const OrientedReadId orientedReadId = itEntrance->orientedReadId;
        DINARA_ASSERT(orientedReadId == itExit->orientedReadId);

        if(not tangleMatrix1.goesBackward(orientedReadId)) {
            const Shasta2Journey journey = journeys[orientedReadId];
            const uint32_t entrancePositionInJourney = itEntrance->positionInJourney;
            const uint32_t exitPositionInJourney = itExit->positionInJourney;
            const Shasta2AnchorId entranceAnchorId = journey[entrancePositionInJourney];
            const Shasta2AnchorId exitAnchorId = journey[exitPositionInJourney];
            const uint32_t entranceOrdinal = anchors.getOrdinal(entranceAnchorId, orientedReadId);
            const uint32_t exitOrdinal = anchors.getOrdinal(exitAnchorId, orientedReadId);
            const auto orientedReadMarkers = markers[orientedReadId.getValue()];
            const uint64_t offset =
                orientedReadMarkers[exitOrdinal].position -
                orientedReadMarkers[entranceOrdinal].position;
            offsetSum += offset;
            ++offsetCount;
        }

        ++itEntrance;
        ++itExit;
    }
    const double averageOffset = double(offsetSum) / double(offsetCount);
    const uint32_t maxOffset = uint32_t(std::round(averageOffset * (1. + drift)));



    // Joint loop over the oriented reads of the entrance and exit.
    itEntrance = entranceOrientedReadInfos.begin();
    itExit = exitOrientedReadInfos.begin();
    while(true) {

        // If both iterators are at their end, we are done.
        if((itEntrance==itEntranceEnd) and (itExit==itExitEnd)) {
            break;
        }



        // Case where the OrientedReadId appears only in the entrance.
        // In this case we use a journey portion that begins at the
        // journey position stored in the OrientedReadInfo and
        // ends at a distance maxOffset after that.
        if(
            (itEntrance != itEntranceEnd) and
            (itExit==itExitEnd or (itEntrance->orientedReadId < itExit->orientedReadId))) {
            const OrientedReadId orientedReadId = itEntrance->orientedReadId;
            if(not tangleMatrix1.goesBackward(orientedReadId)) {
                const Shasta2Journey journey = journeys[orientedReadId];
                const uint32_t begin = itEntrance->positionInJourney;
                uint32_t end = uint32_t(journey.size());

                // End the Shasta2JourneyPortion at a distance maxOffset after begin.
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];
                uint32_t beginPosition = invalid<uint32_t>;
                for(uint32_t positionInJourney=begin; positionInJourney<end; positionInJourney++) {
                    const Shasta2AnchorId anchorId = journey[positionInJourney];
                    const uint32_t ordinal = anchors.getOrdinal(anchorId, orientedReadId);
                    const uint32_t position = orientedReadMarkers[ordinal].position;
                    if(positionInJourney == begin) {
                        beginPosition = position;
                    } else {
                        const uint32_t offset = position - beginPosition;
                        if(offset >maxOffset) {
                            end = positionInJourney;
                            break;
                        }
                    }
                }

                journeyPortions.emplace_back(orientedReadId, begin, end);
                // cout << "Entrance only " << orientedReadId << endl;
            }
            ++itEntrance;
        }



        // Case where the OrientedReadId appears only in the exit.
        // In this case we use a journey portion that begins at the
        // beginning of the journey and ends at the
        // journey position stored in the OrientedReadInfo.
        else if(
            (itExit != itExitEnd) and
            (itEntrance==itEntranceEnd or (itExit->orientedReadId < itEntrance->orientedReadId))) {
            const OrientedReadId orientedReadId = itExit->orientedReadId;
            if(not tangleMatrix1.goesBackward(orientedReadId)) {
                const Shasta2Journey journey = journeys[orientedReadId];
                uint32_t begin = 0;
                const uint32_t end = itExit->positionInJourney + 1;

                // Begin the Shasta2JourneyPortion at a distance maxOffset before begin.
                const auto orientedReadMarkers = markers[orientedReadId.getValue()];
                uint32_t endPosition = invalid<uint32_t>;
                for(uint32_t positionInJourney=end-1; /* Check later */ ; positionInJourney--) {
                    const Shasta2AnchorId anchorId = journey[positionInJourney];
                    const uint32_t ordinal = anchors.getOrdinal(anchorId, orientedReadId);
                    const uint32_t position = orientedReadMarkers[ordinal].position;
                    if(positionInJourney == end - 1) {
                        endPosition = position;
                    } else {
                        const uint32_t offset = endPosition - position;
                        if(offset > maxOffset) {
                            begin = positionInJourney + 1;
                            break;
                        }
                    }

                    if(positionInJourney == 0) {
                        break;
                    }
                }

                journeyPortions.emplace_back(orientedReadId, begin, end);
                // cout << "Exit only " << orientedReadId << endl;
            }
            ++itExit;
        }



        // Case where the OrientedReadId appears in both the entrance and the exit.
        // In this case we use a journey portion that begins at the
        // journey position stored in the entrance OrientedReadInfo and
        // ends at the journey position stored in the exit OrientedReadInfo.
        else if(
            (itEntrance != itEntranceEnd) and
            (itExit!=itExitEnd) and
            (itEntrance->orientedReadId == itExit->orientedReadId)) {
            const OrientedReadId orientedReadId = itEntrance->orientedReadId;
            if(not tangleMatrix1.goesBackward(orientedReadId)) {
                const uint32_t begin = itEntrance->positionInJourney;
                const uint32_t end = itExit->positionInJourney + 1;
                journeyPortions.emplace_back(orientedReadId, begin, end);
                // cout << "Both " << orientedReadId << endl;
            }
            ++itEntrance;
            ++itExit;
        }

        // All cases have been covered.
        else {
            DINARA_ASSERT(0);
        }
    }

    if(html) {
        html << "<h4>Shasta2Journey portions used in this Shasta2RestrictedAnchorGraph</h4>"
            "<table><tr>"
            "<th>OrientedReadId<th>Begin<th>End<th>Length";
        for(const Shasta2JourneyPortion journeyPortion: journeyPortions) {
            html << "<tr>"
                "<td class=centered>" << journeyPortion.orientedReadId <<
                "<td class=centered>" << journeyPortion.begin <<
                "<td class=centered>" << journeyPortion.end <<
                "<td class=centered>" << journeyPortion.end - journeyPortion.begin;
        }
        html << "</table>";
    }
}



// Create a new vertex and add it to the vertexTable, without resorting
// the vertexTable.
Shasta2RestrictedAnchorGraph::vertex_descriptor Shasta2RestrictedAnchorGraph::addVertex(Shasta2AnchorId anchorId)
{

    // Add the vertex and update the anchorMap.
    vertex_descriptor v = null_vertex();
    if(cycleAvoider) {
        v = cycleAvoider->addVertex();
        (*this)[v].anchorId = anchorId;
    }
    else {
        v = add_vertex(Shasta2RestrictedAnchorGraphVertex(anchorId), *this);
    }

    // Store it in the anchorMap.
    auto p = anchorMap.getExisting(anchorId);
    DINARA_ASSERT(p);
    AnchorInformation& anchorInformation = p->second;
    DINARA_ASSERT(anchorInformation.v == null_vertex());
    anchorInformation.v = v;


    return v;
}



void Shasta2RestrictedAnchorGraph::writeGraphviz(
    const string& fileName,
    const vector<Shasta2AnchorId>& highlightVertices) const
{
    ofstream dot(fileName);
    writeGraphviz(dot, highlightVertices);
}



void Shasta2RestrictedAnchorGraph::writeGraphviz(
    ostream& dot,
    const vector<Shasta2AnchorId>& highlightVertices) const
{
    using Graph = Shasta2RestrictedAnchorGraph;
    const Graph& graph = *this;

    dot << "digraph Shasta2RestrictedAnchorGraph {\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Shasta2AnchorId anchorId = graph[v].anchorId;
        dot <<
            "\"" << shasta2AnchorIdToString(anchorId) << "\""
            " [label=\"" << shasta2AnchorIdToString(anchorId) << "\\n" << graph[v].orientedReadIds.size() << "\"";
        if(find(highlightVertices.begin(), highlightVertices.end(), anchorId) != highlightVertices.end()) {
            dot << " style=filled fillcolor=pink";
        }
        dot << "];\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Shasta2AnchorId anchorId0 = graph[v0].anchorId;
        const Shasta2AnchorId anchorId1 = graph[v1].anchorId;
        const Shasta2RestrictedAnchorGraphEdge& edge = graph[e];
        const uint64_t coverage = edge.anchorPair.size();
        const uint64_t offset = edge.offset;

        dot << "\"" << shasta2AnchorIdToString(anchorId0) << "\"->\"" <<
            shasta2AnchorIdToString(anchorId1) << "\""
            "["
            "penwidth=" << std::setprecision(2) << 0.5 * double(coverage) <<
            " label=\"" << coverage << "\\n" << offset << "\"";

        if(edge.isOptimalPathEdge) {
            dot << " color=green";
        }

        dot << "];\n";
    }

    dot << "}\n";
}



void Shasta2RestrictedAnchorGraph::writeHtml(
    ostream& html,
    const vector<Shasta2AnchorId>& highlightVertices) const
{
    // Write it in graphviz format in a temporary location.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName =
        (std::filesystem::temp_directory_path() / (uuid + ".dot")).string();
    writeGraphviz(dotFileName, highlightVertices);

    // Keep html output local to avoid extra binary dependencies.
    html << "<p>Wrote Graphviz file " << dotFileName << ".";
}



// Only keep vertices that are forward reachable from the
// vertex at anchorId0 and backward reachable from the vertex at anchorId1.
// To permit future optimizations, we don't really remove
// the vertices - we only disconnect them for the rest of the graph.
void Shasta2RestrictedAnchorGraph::keepBetween(Shasta2AnchorId anchorId0, Shasta2AnchorId anchorId1)
{
    using Graph = Shasta2RestrictedAnchorGraph;
    Graph& graph = *this;

    const vertex_descriptor v0 = getExistingVertex(anchorId0);
    const vertex_descriptor v1 = getExistingVertex(anchorId1);

    std::set<vertex_descriptor> reachableVertices;
    vector<vertex_descriptor> verticesToBeRemoved;

    // Remove vertices that are not forward reachable from v0.
    // To permit future optimizations, we don't really remove
    // the vertices - we only disconnect them for the rest of the graph.
    shasta2::findReachableVertices(graph, v0, 0, reachableVertices);
    DINARA_ASSERT(reachableVertices.contains(v1));
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not reachableVertices.contains(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        // vertexMap.erase(graph[v].anchorId);
        boost::clear_vertex(v, graph);
        // boost::remove_vertex(v, graph);
    }
    reachableVertices.clear();
    verticesToBeRemoved.clear();

    // Remove vertices that are not backward reachable from v1.
    // To permit future optimizations, we don't really remove
    // the vertices - we only disconnect them for the rest of the graph.
    shasta2::findReachableVertices(graph, v1, 1, reachableVertices);
    DINARA_ASSERT(reachableVertices.contains(v0));
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not reachableVertices.contains(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        // vertexMap.erase(graph[v].anchorId);
        boost::clear_vertex(v, graph);
        // boost::remove_vertex(v, graph);
    }
}




void Shasta2RestrictedAnchorGraph::approximateTopologicalSort()
{
    using Graph = Shasta2RestrictedAnchorGraph;
    Graph& graph = *this;

    // Gather the edges with their coverage.
    vector< pair<edge_descriptor, uint64_t> > edgeTable;
    BGL_FORALL_EDGES(e, graph, Graph)
    {
        const uint64_t coverage = graph[e].anchorPair.size();
        edgeTable.push_back(make_pair(e, coverage));
    }

    // Sort by decreasing coverage.
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());
    vector<edge_descriptor> edgesSortedByDecreasingCoverage;
    for(const auto& p: edgeTable) {
        edgesSortedByDecreasingCoverage.push_back(p.first);
    }

    // Do the approximate topological sort.
    shasta2::approximateTopologicalSort(graph, edgesSortedByDecreasingCoverage);
}



// Find the optimal assembly path.
void Shasta2RestrictedAnchorGraph::findOptimalPath(
    Shasta2AnchorId anchorId0,
    Shasta2AnchorId anchorId1,
    vector<edge_descriptor>& optimalPath)
{
    using Graph = Shasta2RestrictedAnchorGraph;
    Graph& graph = *this;

    // Find the vertices corresponding to anchorId0 and anchorId1.
    const vertex_descriptor v0 = getExistingVertex(anchorId0);
    const vertex_descriptor v1 = getExistingVertex(anchorId1);

    // Compute the longest path.
    shasta2::longestPath(graph, optimalPath);
    if(
        optimalPath.empty() or
        (source(optimalPath.front(), graph) != v0) or
        (target(optimalPath.back(), graph) != v1)) {
        ofstream dot("RestrictedAnchorGraphFailure.dot");
        writeGraphviz(dot, vector<Shasta2AnchorId>{anchorId0, anchorId1});
    }
    DINARA_ASSERT(not optimalPath.empty());
    DINARA_ASSERT(source(optimalPath.front(), graph) == v0);
    DINARA_ASSERT(target(optimalPath.back(), graph) == v1);


    // Set the isOptimalPathEdge flags.
    for(const edge_descriptor e: optimalPath) {
        graph[e].isOptimalPathEdge = true;
    }


}



// Write a table showing which OrientedReadIds are in each vertex.
// Vertices are written out in rank order.
void Shasta2RestrictedAnchorGraph::writeOrientedReadsInVertices(ostream& html) const
{
    using Graph = Shasta2RestrictedAnchorGraph;
    const Graph& graph = *this;

    // Gather the vertices in rank order.
    vector< pair<vertex_descriptor, uint64_t> > verticesSortedByRank;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        verticesSortedByRank.push_back(make_pair(v, graph[v].rank));
    }
    std::ranges::sort(verticesSortedByRank, OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());

    // Write the table header.
    html << "<table><tr><th>Oriented<br>read<br<id>";
    for(const auto& p: verticesSortedByRank) {
        const vertex_descriptor v = p.first;
        const Shasta2AnchorId anchorId = graph[v].anchorId;
        html << "<th>" << shasta2AnchorIdToString(anchorId);
    }



    // Write one row for each OrientedReadId.
    for(const Shasta2JourneyPortion& journeyPortion: journeyPortions) {
        const OrientedReadId orientedReadId = journeyPortion.orientedReadId;

        html << "<tr><th>" << orientedReadId;

        for(const auto& p: verticesSortedByRank) {
            const vertex_descriptor v = p.first;
            const bool vertexContainsOrientedReadId =
                std::binary_search(
                    graph[v].orientedReadIds.begin(),
                    graph[v].orientedReadIds.end(),
                    orientedReadId);

            html << "<td";
            if(vertexContainsOrientedReadId) {
                html << " style='background-color:LightGreen'";
            }
            html << ">" << int(vertexContainsOrientedReadId);
        }
    }


    html << "</table>";
}



void Shasta2RestrictedAnchorGraph::writeAnchorMap() const
{
    cout << "AnchorMap:" << endl;
    DINARA_ASSERT(0);
}



void Shasta2RestrictedAnchorGraph::writeAnchorIds() const
{
    cout << "AnchorIds:" << endl;
    for(uint64_t anchorIndex=0; anchorIndex<anchorIds.size(); anchorIndex++) {
    }
    DINARA_ASSERT(0);
}
