// Dinara.
#include "Shasta2LocalAssembly5.hpp"
#include "Shasta2AbpoaWrapper.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2AnchorPair.hpp"
#include "approximateTopologicalSort.hpp"
#include "CondensedGraph.hpp"
#include "deduplicate.hpp"
#include "dominatorTree.hpp"
#include "findReachableVertices.hpp"
#include "graphvizToHtml.hpp"
#include "Kmer.hpp"
#include "longestPath.hpp"
#include "MarkerKmers.hpp"
#include "Marker.hpp"
#include "orderPairs.hpp"
#include "Shasta2PoastaWrapper.hpp"
#include "Reads.hpp"
#include "tmpDirectory.hpp"
#include "transitiveReduction.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "algorithm.hpp"
#include <cmath>
#include "fstream.hpp"
#include <queue>



Shasta2LocalAssembly5::Shasta2LocalAssembly5(
    const Shasta2Anchors& anchors,
    uint64_t abpoaMaxLength,
    ostream& html,
    [[maybe_unused]] bool debug,
    const Shasta2AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) :
    anchors(anchors),
    abpoaMaxLength(abpoaMaxLength),
    html(html),
    leftAnchorId(anchorPair.anchorIdA),
    rightAnchorId(anchorPair.anchorIdB)
{
    // Sanity check.
    DINARA_ASSERT(std::ranges::is_sorted(additionalOrientedReadIds));

    // Summarize the input for this local assembly.
    if(html) {
        writeInput(anchorPair, additionalOrientedReadIds);
    }

    // Gather the oriented reads that could possibly participate in this local assembly.
    gatherAllOrientedReads(anchorPair, additionalOrientedReadIds);
    if(false) {
        writeAllOrientedReadIds();
    }

    // Gather the MarkerInfos on the left and right marker
    // for this local assembly.
    fillMarkerInfos();
    if(false) {
        writeMarkerInfos();
    }

    // Gather the oriented reads that will actually be used in this local assembly.
    gatherOrientedReads();

    // Estimate the base offset for this local assembly using the oriented reads
    // that appear on both the left and right markers.
    estimateOffset();

    // Fill in the LocalRegions of all oriented reads.
    // These are the oriented read regions that will be used in this local assembly.
    fillLocalRegions();

    // Gather the Kmers that will be used in this local assembly.
    gatherKmers();
    if(html) {
        writeKmers();
        writeOrientedReads();
    }

    // Create the initial graph.
    createGraph();
    removeInaccessibleVertices();
    computeDominatorTree();
    if(html) {
        html << "<h3>Initial Shasta2LocalAssembly5 graph</h3>";
        writeGraph();
    }

    // Graph cleanup.
    removeLowCoverageEdges();
    DINARA_ASSERT(isReachable(*this, vLeft, vRight, 0));
    removeInaccessibleVertices();

    // Compute the assembly path.
    computeAssemblyPath();

    if(html) {
         html << "<h3>Final Shasta2LocalAssembly5 graph</h3>";
         writeGraph();
    }

    // Assemble sequence.
    assemble();
    if(html) {
        writeAssembledSequence();
    }

}



void Shasta2LocalAssembly5::gatherAllOrientedReads(
    const Shasta2AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds)
{
    std::ranges::set_union(anchorPair.orientedReadIds, additionalOrientedReadIds,
        back_inserter(allOrientedReadIds));

}



void Shasta2LocalAssembly5::fillMarkerInfos()
{
    // Get the Kmers.
    leftKmer = anchors.anchorKmer(leftAnchorId);
    rightKmer = anchors.anchorKmer(rightAnchorId);

    // Get the MarkerInfos.
    anchors.markerKmers.get(leftKmer, leftMarkerInfos);
    anchors.markerKmers.get(rightKmer, rightMarkerInfos);
}



// This gathers the oriented reads that will be used used in this local assembly.
// These are all OrientedReadIds that are in allOrientedReadIds and
// also on the left or right Marker.
// This can be done more efficiently but the code would be less readable.
void Shasta2LocalAssembly5::gatherOrientedReads()
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    // Comparator for MarkerInfos that takes into account just the OrientedReadId.
    // This is safe because the left and right Markers are Shasta2Anchors,
    // and so an OrientedReadId can appear more than once.
    class MarkerInfoComparator {
    public:
        bool operator() (const MarkerKmers::MarkerInfo& x, const MarkerKmers::MarkerInfo& y) const
        {
            return x.orientedReadId < y.orientedReadId;
        }
    };
    MarkerInfoComparator markerInfoComparator;



    // Loop over allOrientedReadIds.
    for(const OrientedReadId orientedReadId: allOrientedReadIds) {

        // Look for left and right MarkerInfos containing this oriented read.
        const auto itLeft = std::lower_bound(leftMarkerInfos.begin(), leftMarkerInfos.end(),
            MarkerKmers::MarkerInfo(orientedReadId, 0), markerInfoComparator);
        const auto itRight = std::lower_bound(rightMarkerInfos.begin(), rightMarkerInfos.end(),
            MarkerKmers::MarkerInfo(orientedReadId, 0), markerInfoComparator);

        const bool isOnLeftMarker =
            (itLeft != leftMarkerInfos.end()) and
            (itLeft->orientedReadId == orientedReadId);
        const bool isOnRightMarker =
            (itRight != rightMarkerInfos.end()) and
            (itRight->orientedReadId == orientedReadId);

        // We can use it if it is on at least one of the left and right markers.
        // If on both, we also have to check that the right ordinal is greater than the left ordinal
        // (order violation).
        if(isOnLeftMarker or isOnRightMarker) {

            // Check for order violation.
            const bool isOnBothMarkers = isOnLeftMarker and isOnRightMarker;
            const bool orderViolation = isOnBothMarkers and (itRight->ordinal <= itLeft->ordinal);
            if(not orderViolation) {

                // Ok, we can use this oriented read.
                const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
                OrientedReadInfo& orientedReadInfo = orientedReadInfos.emplace_back();
                orientedReadInfo.orientedReadId = orientedReadId;
                orientedReadInfo.isOnLeftMarker = isOnLeftMarker;
                orientedReadInfo.isOnRightMarker = isOnRightMarker;
                if(isOnLeftMarker) {
                    orientedReadInfo.leftOrdinal = itLeft->ordinal;
                    orientedReadInfo.leftPosition = orientedReadMarkers[itLeft->ordinal].position + kHalf;
                }
                if(isOnRightMarker) {
                    orientedReadInfo.rightOrdinal = itRight->ordinal;
                    orientedReadInfo.rightPosition = orientedReadMarkers[itRight->ordinal].position + kHalf;
                }
            }
        }
    }
}



void Shasta2LocalAssembly5::estimateOffset()
{
    uint64_t sum = 0;
    uint64_t n = 0;
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        if(orientedReadInfo.isOnBothMarkers()) {
            sum += orientedReadInfo.positionOffset();
            ++n;
        }
    }

    DINARA_ASSERT(n > 0);
    offset = uint32_t(std::round(double(sum) / double(n)));

    if(html) {
        html << "<br>Estimated offset using " << n <<
            " oriented reads on the left and right markers is " << offset << " bases.";
    }
}



void Shasta2LocalAssembly5::fillLocalRegions()
{
    const uint32_t length = uint32_t(std::round((1. + drift) * double(offset)));

    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        orientedReadInfo.fillLocalRegion(anchors, length);
    }
}



void Shasta2LocalAssembly5::OrientedReadInfo::fillLocalRegion(
    const Shasta2Anchors& anchors,
    uint32_t length)
{
    const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    if(isOnBothMarkers()) {
        localRegion.firstOrdinal = leftOrdinal;
        localRegion.lastOrdinal = rightOrdinal;
    }

    else if(isOnLeftMarker) {

        localRegion.firstOrdinal = leftOrdinal;
        localRegion.lastOrdinal = leftOrdinal;

        for(uint32_t lastOrdinal=localRegion.firstOrdinal; lastOrdinal<orientedReadMarkers.size(); ++lastOrdinal) {
            if(orientedReadMarkers[lastOrdinal].position + kHalf - leftPosition < length) {
                localRegion.lastOrdinal = lastOrdinal;
            } else {
                break;
            }
        }
    }

    else if(isOnRightMarker) {
        localRegion.lastOrdinal = rightOrdinal;

        if(localRegion.lastOrdinal == 0) {
            localRegion.firstOrdinal = 0;
        } else {

            // Move left by at least length bases.
            for(localRegion.firstOrdinal = localRegion.lastOrdinal - 1;
                /* Check later */ ; localRegion.firstOrdinal--) {
                if(rightPosition - orientedReadMarkers[localRegion.firstOrdinal].position + kHalf >= length) {
                    break;
                }
                if(localRegion.firstOrdinal == 0) {
                    break;
                }
            }
        }

    }

    else {
        DINARA_ASSERT(0);
    }
    DINARA_ASSERT(localRegion.firstOrdinal < orientedReadMarkers.size());
    DINARA_ASSERT(localRegion.lastOrdinal < orientedReadMarkers.size());


    localRegion.firstPosition = orientedReadMarkers[localRegion.firstOrdinal].position + kHalf;
    localRegion.lastPosition = orientedReadMarkers[localRegion.lastOrdinal].position + kHalf;


    // Fill in the LocalRegion Kmers.
    for(uint32_t ordinal = localRegion.firstOrdinal;
        ordinal <= localRegion.lastOrdinal; ordinal++) {
        localRegion.kmers.push_back(anchors.markers.getKmer(orientedReadId, ordinal));
    }

    // Fill in nonUniqueKmers.
    localRegion.computeNonUniqueKmers();
}



void Shasta2LocalAssembly5::OrientedReadInfo::LocalRegion::computeNonUniqueKmers()
{
    nonUniqueKmers = kmers;
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(nonUniqueKmers, count, 2UL);
}



void Shasta2LocalAssembly5::writeInput(
    const Shasta2AnchorPair& anchorPair,
    const vector<OrientedReadId>& additionalOrientedReadIds) const
{
    html <<
        "<table>"
        "<tr><th class=left>Left anchor<td class=centered>" << shasta2AnchorIdToString(leftAnchorId) <<
        "<tr><th class=left>Right anchor<td class=centered>" << shasta2AnchorIdToString(rightAnchorId) <<
        "</table>";
    html << "<h3>OrientedReadIds in the Shasta2AnchorPair</h3><table>";
    for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << anchorPair.orientedReadIds.size() << " oriented reads in the Shasta2AnchorPair.";

    html << "<h3>Additional OrientedReadIds</h3>"
        "These are additional OrientedReadIds that can be used in this assembly."
        "<table>";
    for(const OrientedReadId orientedReadId: additionalOrientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << additionalOrientedReadIds.size() << " additional oriented reads.";

}



void Shasta2LocalAssembly5::writeAllOrientedReadIds() const
{
    html << "<h3>All OrientedReadIds</h3>"
        "These are the union of the OrientedReadIds in the Shasta2AnchorPair "
        "and the additional OrientedReadIds"
        "<table>";
    for(const OrientedReadId orientedReadId: allOrientedReadIds) {
        html << "<tr><td class=centered>" << orientedReadId;
    }
    html << "</table>" << allOrientedReadIds.size() << " oriented reads.";
}



void Shasta2LocalAssembly5::writeMarkerInfos() const
{
    writeMarkerInfos("Left", leftMarkerInfos);
    writeMarkerInfos("Right", rightMarkerInfos);
}



void Shasta2LocalAssembly5::writeMarkerInfos(
    const string& side,
    const vector<MarkerKmers::MarkerInfo>& markerInfos) const
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    html << "<h3>" << side << " marker</h3>"
        "<table>"
        "<tr><th>OrientedReadId<th>Ordinal<th>Position";
    for(const MarkerKmers::MarkerInfo& markerInfo: markerInfos) {
        const auto orientedReadMarkers = anchors.markers[markerInfo.orientedReadId.getValue()];
        html <<
            "<tr>"
            "<td class=centered>" << markerInfo.orientedReadId <<
            "<td class=centered>" << markerInfo.ordinal <<
            "<td class=centered>" << orientedReadMarkers[markerInfo.ordinal].position + kHalf;
    }
    html << "</table>" << markerInfos.size() << " oriented reads.";

}



void Shasta2LocalAssembly5::writeOrientedReads() const
{
    html <<
        "<h3>Oriented reads used in this local assembly</h3>"
        "<table><tr>"
        "<th>OrientedReadId"
        "<th>Marker<br>Count"
        "<th>Left<br>marker<br>ordinal"
        "<th>Right<br>marker<br>ordinal"
        "<th>Left<br>to right<br>ordinal<br>offset"
        "<th>Local<br>region<br>first<br>ordinal"
        "<th>Local<br>region<br>last<br>ordinal"
        "<th>Local<br>region<br>ordinal<br>offset"
        "<th>Local<br>region<br>number<br>of markers"
        "<th>Number<br>of markers<br>for assembly"
        "<th>Read<br>length"
        "<th>Left<br>marker<br>position"
        "<th>Right<br>marker<br>position"
        "<th>Left<br>to right<br>position<br>offset"
        "<th>Local<br>region<br>first<br>position"
        "<th>Local<br>region<br>last<br>position"
        "<th>Local<br>region<br>position<br>offset<br>(length)";

    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        const uint64_t markerCount = anchors.markers[orientedReadInfo.orientedReadId.getValue()].size();
        const uint64_t baseCount = anchors.reads.getRead(orientedReadInfo.orientedReadId.getReadId()).baseCount;


        html <<
            "<tr><td class=centered>" << orientedReadInfo.orientedReadId <<
            "<td class=centered>" << markerCount;



        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftMarker) {
            html << orientedReadInfo.leftOrdinal;
        }

        html << "<td class=centered>";
        if(orientedReadInfo.isOnRightMarker) {
            html << orientedReadInfo.rightOrdinal;
        }

        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftMarker and orientedReadInfo.isOnRightMarker) {
            html << orientedReadInfo.rightOrdinal - orientedReadInfo.leftOrdinal;
        }

        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.firstOrdinal;
        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.lastOrdinal;
        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.ordinalOffset();
        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.kmers.size();
        html << "<td class=centered>";
        html << orientedReadInfo.ordinalsForAssembly.size();

        html << "<td class=centered>" << baseCount;

        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftMarker) {
            html << orientedReadInfo.leftPosition;
        }

        html << "<td class=centered>";
        if(orientedReadInfo.isOnRightMarker) {
            html << orientedReadInfo.rightPosition;
        }

        html << "<td class=centered>";
        if(orientedReadInfo.isOnLeftMarker and orientedReadInfo.isOnRightMarker) {
            html << orientedReadInfo.rightPosition - orientedReadInfo.leftPosition;
        }

        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.firstPosition;
        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.lastPosition;
        html << "<td class=centered>";
        html << orientedReadInfo.localRegion.positionOffset();
    }

    html << "</table>" << orientedReadInfos.size() << " oriented reads.";

}



void Shasta2LocalAssembly5::gatherKmers()
{
    // Gather all Kmers that are not unique
    // in the LocalRegion of one or more reads.
    // These wil be excluded from the assembly.
    vector<Kmer> allNonUniqueKmers;
    for(const OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        std::ranges::copy(orientedReadInfo.localRegion.nonUniqueKmers, back_inserter(allNonUniqueKmers));
    }
    deduplicate(allNonUniqueKmers);

    if(html) {
        html << "<br>Found " << allNonUniqueKmers.size() <<
            " marker k-mers that appear more than once in the local rethis local assembly graph generation.";
    }


    // Gather the Kmers that are unique in the LocalRegions of all oriented reads
    // and count how many times each of them appears.
    kmers.clear();
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        for(uint32_t ordinal = orientedReadInfo.localRegion.firstOrdinal;
            ordinal <= orientedReadInfo.localRegion.lastOrdinal; ordinal++) {

            const Kmer& kmer = orientedReadInfo.localRegion.getKmer(ordinal);
            if(not binary_search(allNonUniqueKmers.begin(), allNonUniqueKmers.end(), kmer)) {
                kmers.push_back(kmer);
            }
        }
    }
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(kmers, count, minVertexCoverage);


    // Make sure the left and right Kmers are included.
    const bool addLeftKmer = not binary_search(kmers.begin(), kmers.end(), leftKmer);
    const bool addRightKmer = not binary_search(kmers.begin(), kmers.end(), rightKmer);
    if(addLeftKmer) {
        kmers.push_back(leftKmer);
    }
    if(addRightKmer) {
        kmers.push_back(rightKmer);
    }
    if(addLeftKmer or addRightKmer) {
        sort(kmers.begin(), kmers.end());
    }
    DINARA_ASSERT(binary_search(kmers.begin(), kmers.end(), leftKmer));
    DINARA_ASSERT(binary_search(kmers.begin(), kmers.end(), rightKmer));

#if 0
    cout << "AAA " << kmers.size() << endl;
    for(uint64_t kmerId=0; kmerId<kmers.size(); kmerId++) {
        cout << kmerId << " ";
        kmers[kmerId].write(cout, anchors.k);
        cout << endl;
    }
#endif


    // Fill in the ordinalsForAssembly in the OrientedReadInfos.
    for(OrientedReadInfo& orientedReadInfo: orientedReadInfos) {
        for(uint32_t ordinal = orientedReadInfo.localRegion.firstOrdinal;
            ordinal <= orientedReadInfo.localRegion.lastOrdinal; ordinal++) {

            const Kmer& kmer = orientedReadInfo.localRegion.getKmer(ordinal);
            if(binary_search(kmers.begin(), kmers.end(), kmer)) {
                orientedReadInfo.ordinalsForAssembly.push_back(ordinal);
            }
        }
    }
}



// Create the graph using the current ordinalsForAssembly.
void Shasta2LocalAssembly5::createGraph()
{
    using Graph = Shasta2LocalAssembly5;
    Graph& graph = *this;

    DINARA_ASSERT(binary_search(kmers.begin(), kmers.end(), leftKmer));
    DINARA_ASSERT(binary_search(kmers.begin(), kmers.end(), rightKmer));

    // Clear first, just in case.
    graph.clear();

    // Generate vertices.
    // Each Kmer generates a vertex.
    vector<vertex_descriptor> vertexTable;
    vLeft = null_vertex();
    vRight = null_vertex();
    for(uint64_t kmerId=0; kmerId<kmers.size(); kmerId++) {
        const Kmer& kmer = kmers[kmerId];
        const vertex_descriptor v = add_vertex(Shasta2LocalAssembly5Vertex(kmerId), graph);
        vertexTable.push_back(v);
        if(kmer == leftKmer) {
            vLeft = v;
        }
        if(kmer == rightKmer) {
            vRight = v;
        }
    }
    DINARA_ASSERT(vLeft != null_vertex());
    DINARA_ASSERT(vRight != null_vertex());



    // Fill in the infos of the vertices.
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadInfos.size(); orientedReadIndex++) {
        const OrientedReadInfo& orientedReadInfo = orientedReadInfos[orientedReadIndex];
        for(uint32_t ordinal: orientedReadInfo.ordinalsForAssembly) {
            const Kmer& kmer = orientedReadInfo.localRegion.getKmer(ordinal);
            const auto it = std::lower_bound(kmers.begin(), kmers.end(), kmer);

           DINARA_ASSERT(it != kmers.end());
            DINARA_ASSERT(*it == kmer);
            const vertex_descriptor v = vertexTable[it - kmers.begin()];
            Shasta2LocalAssembly5Vertex& vertex = graph[v];
            auto& info = vertex.infos.emplace_back();
            info.orientedReadIndex = orientedReadIndex;
            info.ordinal = ordinal;
        }
    }


    // Generate the edges.
    for(uint64_t orientedReadIndex=0; orientedReadIndex<orientedReadInfos.size(); orientedReadIndex++) {
        const OrientedReadInfo& orientedReadInfo = orientedReadInfos[orientedReadIndex];
        for(uint64_t i1=1; i1<orientedReadInfo.ordinalsForAssembly.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const uint32_t ordinal0 = orientedReadInfo.ordinalsForAssembly[i0];
            const uint32_t ordinal1 = orientedReadInfo.ordinalsForAssembly[i1];
            const Kmer& kmer0 = orientedReadInfo.localRegion.getKmer(ordinal0);
            const Kmer& kmer1 = orientedReadInfo.localRegion.getKmer(ordinal1);
            const auto it0 = std::lower_bound(kmers.begin(), kmers.end(), kmer0);
            DINARA_ASSERT(it0 != kmers.end());
            DINARA_ASSERT(*it0 == kmer0);
            const auto it1 = std::lower_bound(kmers.begin(), kmers.end(), kmer1);
            DINARA_ASSERT(it1 != kmers.end());
            DINARA_ASSERT(*it1 == kmer1);
            const vertex_descriptor v0 = vertexTable[it0 - kmers.begin()];
            const vertex_descriptor v1 = vertexTable[it1 - kmers.begin()];
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            if(not edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
                DINARA_ASSERT(edgeExists);
            }
            Shasta2LocalAssembly5Edge& edge = graph[e];
            auto& info = edge.infos.emplace_back();
            info.orientedReadIndex = orientedReadIndex;
            info.ordinal0 = ordinal0;
            info.ordinal1 = ordinal1;
        }
    }



    // Fill in the extendedInfos of the edges.
    BGL_FORALL_EDGES(e, graph, Graph) {
        Shasta2LocalAssembly5Edge& edge = graph[e];

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        const Shasta2LocalAssembly5Vertex& vertex0 = graph[v0];
        const Shasta2LocalAssembly5Vertex& vertex1 = graph[v1];

        const vector<Shasta2LocalAssembly5Vertex::Info>& infos0 = vertex0.infos;
        const vector<Shasta2LocalAssembly5Vertex::Info>& infos1 = vertex1.infos;

        // Joint loop over the infos of the two vertices.
        auto it0 = infos0.begin();
        auto it1 = infos1.begin();
        const auto end0 = infos0.end();
        const auto end1 = infos1.end();
        while((it0 != end0) and (it1 != end1)) {
            if(it0->orientedReadIndex < it1->orientedReadIndex) {
                ++it0;
            }
            else if(it1->orientedReadIndex < it0->orientedReadIndex) {
                ++it1;
            } else {
                if(it0->ordinal < it1->ordinal) {
                    Shasta2LocalAssembly5Edge::Info& info = edge.extendedInfos.emplace_back();
                    info.orientedReadIndex = it0->orientedReadIndex;
                    info.ordinal0 = it0->ordinal;
                    info.ordinal1 = it1->ordinal;
                }
                ++it0;
                ++it1;
            }
        }
    }



    if(html) {
        html << "<br>The Shasta2LocalAssembly5 graph has " << num_vertices(graph) <<
            " vertices and " << num_edges(graph) << " edges.";
    }

}



void Shasta2LocalAssembly5::writeGraphviz(const string& fileName)
{
    ofstream dot(fileName);
    writeGraphviz(dot);

}



void Shasta2LocalAssembly5::writeGraphviz(ostream& dot)
{
    using Graph = Shasta2LocalAssembly5;
    Graph& graph = *this;

    // To make the output more readable, we output vertices in
    // approximate topological order.
    vector< pair<edge_descriptor, uint64_t> > edgesWithCoverage;
    BGL_FORALL_EDGES(e, graph, Graph) {
        edgesWithCoverage.push_back(make_pair(e, graph[e].infos.size()));
    }
    sort(edgesWithCoverage.begin(), edgesWithCoverage.end(),
        OrderPairsBySecondOnlyGreater<edge_descriptor, uint64_t>());
    vector<edge_descriptor> edgesOrderedByDecreasingCoverage;
    for(const auto& [e, ignore]: edgesWithCoverage) {
        edgesOrderedByDecreasingCoverage.push_back(e);
    }
    approximateTopologicalSort(graph, edgesOrderedByDecreasingCoverage);
    vector< pair<vertex_descriptor, uint64_t> > verticesOrderedByRank;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        verticesOrderedByRank.push_back(make_pair(v, graph[v].rank));
    }
    sort(verticesOrderedByRank.begin(), verticesOrderedByRank.end(),
        OrderPairsBySecondOnly<vertex_descriptor, uint64_t>());


    dot << "digraph Shasta2LocalAssembly5 {\n";

    for(const auto& [v, ignore]: verticesOrderedByRank) {
        const Shasta2LocalAssembly5Vertex& vertex = graph[v];
        dot << vertex.kmerId;

        // Begin vertex attributes.
        dot << "[";

        // Label.
        dot << "label=\"V" << vertex.kmerId << "\\n" << vertex.infos.size() << "\"";

        // Tooltip.
        dot << " tooltip=\"";
        kmers[vertex.kmerId].write(dot, anchors.k);
        dot << "\"";

        // Color.
        if(v== vLeft) {
            dot << " style=filled fillcolor=LightGreen";
        }
        else if(v== vRight) {
            dot << " style=filled fillcolor=LightPink";
        }
        else if(vertex.isOnDominatorTreePath) {
            dot << " style=filled fillcolor=LightBlue";
        }

        // End vertex attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, Graph) {
        const Shasta2LocalAssembly5Edge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const uint64_t coverage = edge.infos.size();
        const uint64_t extendedCoverage = edge.extendedInfos.size();
        const uint32_t offset = edgeOffset(e);

        dot << graph[v0].kmerId << "->" << graph[v1].kmerId;

        // Begin edge attributes.
        dot << "[";

        // Label.
        dot << "label=\"" << coverage << "/" << extendedCoverage <<
            "\\n" << offset << "\"";

        // Tooltip.
        dot << " tooltip=\"" << coverage << "/" << extendedCoverage <<
            "\\n" << offset << "\"";

        // Thickness.
        dot << " penwidth=" << std::fixed << std::setprecision(2) << 0.3 * double(extendedCoverage);

        // Color.
        if(edge.isOnAssemblyPath) {
            dot << " color=Green";
        }

        // End edge attributes.
        dot << "]";

        // End the line for this edge.
        dot << ";\n";
    }

    dot << "}\n";
}



void Shasta2LocalAssembly5::writeGraph()
{
    // Write it in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    writeGraphviz(dotFileName);


    // Display it in html in svg format.
    const double timeout = 120.;
    const string options = "-Nshape=rectangle -Gbgcolor=gray95";
    html << "<p>";
    try {
        graphvizToHtml(dotFileName, "dot", timeout, options, html);
    } catch(const std::exception& exception) {
        html << "<br>Error during graph layout: " << exception.what();
    }

}



// Remove vertices that are not forward accessible from vLeft
// and backward accessible from vRight.
void Shasta2LocalAssembly5::removeInaccessibleVertices()
{
    using Graph = Shasta2LocalAssembly5;
    Graph& graph = *this;

    std::set<vertex_descriptor> leftAccessibleVertices;
    findReachableVertices(graph, vLeft, 0, leftAccessibleVertices);

    std::set<vertex_descriptor> rightAccessibleVertices;
    findReachableVertices(graph, vRight, 1, rightAccessibleVertices);

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const bool keep =  leftAccessibleVertices.contains(v) and rightAccessibleVertices.contains(v);
        if(not keep) {
            verticesToBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }

}



// Compute a dominator tree starting at vLeft
// and the dominator tree path from vLeft to vRight.
// To deal with possible cycles, we have to compute
// the dominator tree on the condensed graph.
void Shasta2LocalAssembly5::computeDominatorTree()
{
    using Graph = Shasta2LocalAssembly5;
    Graph& graph = *this;

    const bool debug = false;

    // Create the CondensedGraph.
    std::map<Graph::vertex_descriptor, CondensedGraph::vertex_descriptor> vertexMap;
    CondensedGraph condensedGraph =
        createCondensedGraph<Graph, CondensedGraph>(graph, vertexMap);

    // Compute the dominator tree.
    // This sets the dominator field of the vertices
    // reachable by leftAnchorVertex.
    lengauer_tarjan_dominator_tree_general(condensedGraph, vertexMap[vLeft]);

    // Compute the dominator tree path on the CondensedGraph.
    vector<CondensedGraph::vertex_descriptor> condensedDominatorPath;
    CondensedGraph::vertex_descriptor cv = vertexMap[vRight];
    while(cv != CondensedGraph::null_vertex()) {
        condensedDominatorPath.push_back(cv);
        cv = condensedGraph[cv].dominator;
    }
    std::ranges::reverse(condensedDominatorPath);

    if(debug) {
        cout << "Dominator tree path on the condensed graph:";
        for(const CondensedGraph::vertex_descriptor cv: condensedDominatorPath) {
            const vector<vertex_descriptor>& vertices = condensedGraph[cv].vertices;
            if(vertices.size() == 1) {
                const vertex_descriptor v = vertices.front();
                cout << " V" << graph[v].kmerId;
            } else {
                cout << " (" << vertices.size() << " vertices)";
            }
        }
        cout << endl;
    }


    // Sanity checks.
    {
        // The first vertex of the condensedDominatorPath must contain vLeft.
        const CondensedGraph::vertex_descriptor cv = condensedDominatorPath.front();
        const vector<vertex_descriptor>& vertices = condensedGraph[cv].vertices;
        DINARA_ASSERT(std::ranges::find(vertices, vLeft) != vertices.end());
    }
    {
        // The last vertex of the condensedDominatorPath must contain vRight.
        const CondensedGraph::vertex_descriptor cv = condensedDominatorPath.back();
        const vector<vertex_descriptor>& vertices = condensedGraph[cv].vertices;
        DINARA_ASSERT(std::ranges::find(vertices, vRight) != vertices.end());
    }



    // Create the dominator tree path in the Shasta2LocalAssembly5 graph.
    dominatorTreePath.clear();
    dominatorTreePath.push_back(vLeft);
    graph[vLeft].isOnDominatorTreePath = true;
    for(uint64_t i=1; i<condensedDominatorPath.size()-1; i++) {
        const CondensedGraph::vertex_descriptor cv = condensedDominatorPath[i];
        const vector<vertex_descriptor>& vertices = condensedGraph[cv].vertices;
        if(vertices.size() == 1) {
            const vertex_descriptor v = vertices.front();
            dominatorTreePath.push_back(v);
            graph[v].isOnDominatorTreePath = true;
        }
    }
    dominatorTreePath.push_back(vRight);
    graph[vRight].isOnDominatorTreePath = true;

}



// The assembly path is computed by joining together
// partial assembly paths between adjacent vertices
// of the dominator tree path.
void Shasta2LocalAssembly5::removeLowCoverageEdges()
{
    using Graph = Shasta2LocalAssembly5;
    Graph& graph = *this;

    // We use color for BFSs.
    BGL_FORALL_VERTICES(v, graph, Graph) {
        graph[v].color = 0;
    }

    // Remove low coverage edges, using a different threshold
    // for each leg of the dominator tree path.
    for(uint64_t i1=1; i1<dominatorTreePath.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const vertex_descriptor v0 = dominatorTreePath[i0];
        const vertex_descriptor v1 = dominatorTreePath[i1];
        removeLowCoverageEdges(v0, v1);
    }
}



// Given two vertices which are adjacent in the dominator tree path, vA and vB,
// remove low coverage edges in-between without destroying reachability
// of vB from vA.
void Shasta2LocalAssembly5::removeLowCoverageEdges(
    vertex_descriptor vA,
    vertex_descriptor vB)
{
    using Graph = Shasta2LocalAssembly5;
    Graph& graph = *this;
    const bool debug = false;

    if(debug) {
        cout << "Working on assembly path portion between V" <<
            graph[vA].kmerId << " and V" << graph[vB].kmerId << endl;
    }

    // Find all the vertices "between" vA and vB.
    // Use a BFS that starts at vA and stops at vB.
    // This sets the color of these vertices to 1.
    vector<vertex_descriptor> inBetweenVertices;
    std::queue<vertex_descriptor> q;
    q.push(vA);
    graph[vA].color = 1;
    graph[vB].color = 1;
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();
        BGL_FORALL_OUTEDGES(v0, e, graph, Graph) {
            const vertex_descriptor v1 = target(e, graph);
            if(v1 != vB) {
                if(graph[v1].color == 0) {
                    inBetweenVertices.push_back(v1);
                    graph[v1].color = 1;
                    q.push(v1);
                }
            }
        }

    }

    if(debug) {
        cout << "This portion of the graph contains " << inBetweenVertices.size() << " vertices:" << endl;
        for(const vertex_descriptor v: inBetweenVertices) {
            cout << "V" << graph[v].kmerId << " ";
        }
        cout << endl;
    }



    // A vertex predicate that select only vertices with color>0.
    class VertexPredicate {
    public:
        bool operator()(const vertex_descriptor& v) const
        {
            return (*graph)[v].color > 0;
        }
        VertexPredicate(const Graph& graph) : graph(&graph) {}
        VertexPredicate() : graph(0) {}
        const Graph* graph;
    };
    const VertexPredicate vertexPredicate(graph);



    // An edge predicate that select only edges with coverage >= minCoverage.
    uint64_t minCoverage = 0;
    class EdgePredicate {
    public:
        bool operator()(const edge_descriptor& e) const
        {
            return (*graph)[e].extendedInfos.size() >= *minCoveragePointer;
        }
        EdgePredicate(const Graph& graph, uint64_t* minCoveragePointer) :
            graph(&graph),
            minCoveragePointer(minCoveragePointer) {}
        EdgePredicate() : graph(0) {}
        const Graph* graph;
        uint64_t* minCoveragePointer;
    };
    EdgePredicate edgePredicate(graph, &minCoverage);

    // A filtered graph that includes all vertices with color>0 and
    // all the edges with coverage >= minCoverage.
    using FilteredGraph = boost::filtered_graph<Graph, EdgePredicate, VertexPredicate>;
    FilteredGraph filteredGraph(graph, edgePredicate, vertexPredicate);



    // Find the highest edge coverage threshold that ensures reachability.
    for(minCoverage=1; ; ++minCoverage) {
        if(not isReachable(filteredGraph, vA, vB, 0)) {
            --minCoverage;
            break;
        }
    }
    if(debug) {
        cout << "Edge coverage threshold for reachability is " << minCoverage << endl;
    }

    // Remove edges with lower coverage.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_OUTEDGES(vA, e, graph, Graph) {
        if(graph[e].extendedInfos.size() < minCoverage) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const vertex_descriptor v: inBetweenVertices) {
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            if(graph[e].extendedInfos.size() < minCoverage) {
                edgesToBeRemoved.push_back(e);
            }
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        if(debug) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            cout << "Removing edge " << graph[v0].kmerId << "->" << graph[v1].kmerId << endl;
        }
        boost::remove_edge(e, graph);
    }


    // Reset the colors.
    graph[vA].color = 0;
    graph[vB].color = 0;
    for(const vertex_descriptor v: inBetweenVertices) {
        graph[v].color = 0;
    }

}



void Shasta2LocalAssembly5::removeIsolatedVertices()
{
    using Graph = Shasta2LocalAssembly5;
    Graph& graph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if((in_degree(v, graph) == 0) and (out_degree(v, graph)==0)) {
            verticesToBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::remove_vertex(v, graph);
    }
}



void Shasta2LocalAssembly5::writeKmers() const
{
    html <<
        "<h3>Kmers</h3>"
        "<table><tr><th>Vertex<th>Kmer";

    for(uint64_t vertexId=0; vertexId<kmers.size(); vertexId++) {
        html <<
            "<tr><td class=centered>" << vertexId <<
            "<td style='font-family:monospace'>";
        kmers[vertexId].write(html, anchors.k);
    }

    html << "</table>";
}



// Compute position offset for an edge by averaging the
// position offsets of the participating orientedReads.
uint32_t Shasta2LocalAssembly5::edgeOffset(edge_descriptor e) const
{
    using Graph = Shasta2LocalAssembly5;
    const Graph& graph = *this;

    uint64_t offsetSum = 0;
    for(const Shasta2LocalAssembly5Edge::Info& info: graph[e].infos) {
        const OrientedReadInfo& orientedReadInfo = orientedReadInfos[info.orientedReadIndex];
        const OrientedReadId orientedReadId = orientedReadInfo.orientedReadId;
        const uint32_t ordinal0 = info.ordinal0;
        const uint32_t ordinal1 = info.ordinal1;
        const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];

        offsetSum += (orientedReadMarkers[ordinal1].position - orientedReadMarkers[ordinal0].position);
    }

    return uint32_t(std::round(double(offsetSum) / double(graph[e].infos.size())));
}



void Shasta2LocalAssembly5::computeAssemblyPath()
{
    using Graph = Shasta2LocalAssembly5;
    Graph& graph = *this;

    std::map<Graph::vertex_descriptor, CondensedGraph::vertex_descriptor> vertexMap;
    const CondensedGraph condensedGraph =
        createCondensedGraph<Graph, CondensedGraph>(graph, vertexMap);

    // Compute the longest path in the CondensedGraph.
    // This will always work because the CondensedGraph is guaranteed to be acyclic.
    vector<CondensedGraph::edge_descriptor> longestPathCondensed;
    longestPath(condensedGraph, longestPathCondensed);

    // Sanity checks.
    {
        // The longest path must begin at vLeft or at a strongly connected component
        // that contains vLeft.
        const CondensedGraph::edge_descriptor ce = longestPathCondensed.front();
        const vertex_descriptor cv = source(ce, condensedGraph);
        const vector<vertex_descriptor> vertices = condensedGraph[cv].vertices;
        DINARA_ASSERT(std::ranges::find(vertices, vLeft) != vertices.end());
    }
    {
        // The longest path must end at vRight or at a strongly connected component
        // that contains vRight.
        const CondensedGraph::edge_descriptor ce = longestPathCondensed.back();
        const vertex_descriptor cv = target(ce, condensedGraph);
        const vector<vertex_descriptor> vertices = condensedGraph[cv].vertices;
        DINARA_ASSERT(std::ranges::find(vertices, vRight) != vertices.end());
    }


    // Construct the vertices of the longest path in the Shasta2LocalAssembly5 graph.
    // This begins at vLeft, ends at vRight, and its internal vertices
    // are the internal vertices of the longestPathCondensed that correspond
    // to a single vertex of the Shasta2LocalAssembly5 graph
    // (not a non-trivial strongly connected component).
    vector<vertex_descriptor> assemblyPathVertices;
    assemblyPathVertices.push_back(vLeft);
    for(uint64_t i=1; i<longestPathCondensed.size(); i++) {
        const CondensedGraph::edge_descriptor e = longestPathCondensed[i];
        const CondensedGraph::vertex_descriptor v = source(e, condensedGraph);
        const vector<vertex_descriptor>& vertices = condensedGraph[v].vertices;
        if(vertices.size() == 1) {
            const vertex_descriptor v = vertices.front();
            assemblyPathVertices.push_back(v);
        }
    }
    assemblyPathVertices.push_back(vRight);

    if(false) {
        html << "<br>Assembly path vertices:";
        for(const vertex_descriptor v: assemblyPathVertices) {
            html << " " << graph[v].kmerId;
        }
    }


    // Now we create the assembly path.
    // In most cases there is an edge between successive vertices of the assembly path.
    assemblyPath.clear();
    for(uint64_t i1=1; i1<assemblyPathVertices.size(); i1++) {
        const uint64_t i0 = i1 - 1;

        const vertex_descriptor v0 = assemblyPathVertices[i0];
        const vertex_descriptor v1 = assemblyPathVertices[i1];

        auto[e, edgeExists] = boost::edge(v0, v1, graph);
        if(edgeExists) {
            assemblyPath.push_back(e);
            graph[e].isOnAssemblyPath = true;
        } else {

            // Use a shortest path between v0 and v1.
            /*
            cout << "Looking for a shortest path between V" <<
                graph[v0].kmerId << " and V" << graph[v1].kmerId <<
                "." << endl;
            */
            std::queue<vertex_descriptor> q;
            q.push(v0);
            std::map<vertex_descriptor, vertex_descriptor> predecessorMap;
            bool found = false;
            while(not q.empty()) {
                const vertex_descriptor vA = q.front();
                // cout << "Dequeued " << graph[vA].kmerId << endl;
                q.pop();
                BGL_FORALL_OUTEDGES(vA, e, graph, Graph) {
                    const vertex_descriptor vB = target(e, graph);
                    // cout << "Found " << graph[vB].kmerId << endl;
                    if(not predecessorMap.contains(vB)) {
                        // cout << "Queued " << graph[vB].kmerId << endl;
                        predecessorMap.insert({vB, vA});
                        q.push(vB);
                        if(vB == v1) {
                            found = true;
                        }
                    }
                }
                if(found) {
                    break;
                }

            }
            DINARA_ASSERT(found);

            // Walk back the predecessorMap to find the shortest path.
            vector<edge_descriptor> path;
            vertex_descriptor v = v1;
            while(v != v0) {
                const vertex_descriptor vPrevious = predecessorMap[v];
                auto[e, edgeExists] = edge(vPrevious, v, graph);
                DINARA_ASSERT(edgeExists);
                path.push_back(e);
                v = vPrevious;
            }
            std::ranges::reverse(path);

            /*
            cout << "Found this shortest path between V" <<
                graph[v0].kmerId << " and V" << graph[v1].kmerId <<
                ":" << endl;
            for(const edge_descriptor e: path) {
                const vertex_descriptor vA = source(e, graph);
                const vertex_descriptor vB = target(e, graph);
                cout << "V" << graph[vA].kmerId << "->V" << graph[vB].kmerId << endl;

            }
            */

            // Add this shortest path to the assembly path.
            for(const edge_descriptor e: path) {
                assemblyPath.push_back(e);
                graph[e].isOnAssemblyPath = true;
            }

        }
    }

}



// Assemble all sequence.
void Shasta2LocalAssembly5::assemble()
{
    for(const edge_descriptor e: assemblyPath) {
        assemble(e);
    }
}



// Assemble the sequence contributed by a single edge
// of the assembly path.
void Shasta2LocalAssembly5::assemble(edge_descriptor e)
{
    const uint32_t kHalf = uint32_t(anchors.k / 2);

    using Graph = Shasta2LocalAssembly5;
    Graph& graph = *this;
    const Shasta2LocalAssembly5Edge& edge = graph[e];

    const uint64_t sequenceBegin = sequence.size();

    if(html) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        html << "<h3>Assembly details for edge V" << graph[v0].kmerId <<
            "&#x27a1;V" << graph[v1].kmerId <<
            ", coverage " << edge.extendedInfos.size() << "</h3>\n"
            "<table><tr>"
            "<th>Oriented<br>read id"
            "<th>Left<br>ordinal"
            "<th>Right<br>ordinal"
            "<th>Ordinal<br>offset"
            "<th>Left<br>position"
            "<th>Right<br>position"
            "<th>Sequence<br>length"
            "<th class=left>Sequence"
            "\n";

        for(const Shasta2LocalAssembly5Edge::Info& edgeInfo: edge.extendedInfos) {
            const OrientedReadInfo& orientedReadInfo = orientedReadInfos[edgeInfo.orientedReadIndex];
            const OrientedReadId orientedReadId = orientedReadInfo.orientedReadId;

            const uint32_t ordinal0 = edgeInfo.ordinal0;
            const uint32_t ordinal1 = edgeInfo.ordinal1;
            DINARA_ASSERT(ordinal1 > ordinal0);
            const uint32_t ordinalOffset = ordinal1 - ordinal0;

            const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
            const uint32_t position0 = orientedReadMarkers[ordinal0].position + kHalf;
            const uint32_t position1 = orientedReadMarkers[ordinal1].position + kHalf;
            DINARA_ASSERT(position1 > position0);
            const uint32_t positionOffset = position1 - position0;

            html <<
                "<tr>"
                "<td class=centered>" << orientedReadId <<
                "<td class=centered>" << ordinal0 <<
                "<td class=centered>" << ordinal1 <<
                "<td class=centered>" << ordinalOffset <<
                "<td class=centered>" << position0 <<
                "<td class=centered>" << position1 <<
                "<td class=centered>" << positionOffset <<
                "<td class=left style='font-family:monospace;white-space:nowrap'>";
            for(uint32_t position=position0; position!=position1; position++) {
                html << anchors.reads.getOrientedReadBase(orientedReadId, position);
            }

            html << "\n";
        }
        html << "</table>\n";
    }



    // Gather distinct sequences and their coverage.
    class DistinctSequence {
    public:
        shared_ptr< vector<Base> > sequencePointer;
        const vector<Base>& sequence() const
        {
            return *sequencePointer;
        }
        uint64_t coverage;

        // Order by decreasing coverage.
        bool operator<(const DistinctSequence& that) {
            return coverage > that.coverage;
        }
    };
    vector<DistinctSequence> distinctSequences;
    vector<Base> orientedReadSequence;

    for(const Shasta2LocalAssembly5Edge::Info& edgeInfo: edge.extendedInfos) {
        const OrientedReadInfo& orientedReadInfo = orientedReadInfos[edgeInfo.orientedReadIndex];
        const OrientedReadId orientedReadId = orientedReadInfo.orientedReadId;

        const uint32_t ordinal0 = edgeInfo.ordinal0;
        const uint32_t ordinal1 = edgeInfo.ordinal1;
        DINARA_ASSERT(ordinal1 > ordinal0);

        const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
        const uint32_t position0 = orientedReadMarkers[ordinal0].position + kHalf;
        const uint32_t position1 = orientedReadMarkers[ordinal1].position + kHalf;
        DINARA_ASSERT(position1 > position0);

        orientedReadSequence.clear();
        for(uint32_t position=position0; position!=position1; position++) {
            const Base base = anchors.reads.getOrientedReadBase(orientedReadId, position);
            orientedReadSequence.push_back(base);
        }

        // See if this sequence is already in our DistinctSequences.
        bool done = false;
        for(DistinctSequence& distinctSequence: distinctSequences) {
            if(distinctSequence.sequence() == orientedReadSequence) {
                ++distinctSequence.coverage;
                done = true;
                break;
            }
        }

        // If we did not have it, create it with coverage 1.
        if(not done) {
            DistinctSequence& distinctSequence = distinctSequences.emplace_back();
            distinctSequence.sequencePointer = make_shared< vector<Base> >(orientedReadSequence);
            distinctSequence.coverage = 1;
        }
    }
    sort(distinctSequences.begin(), distinctSequences.end());



    // Write the DistinctSequences.
    if(html) {
        html <<
            "<h4>Distinct sequences</h4>"
            "<table><tr>"
            "<th>Coverage<th>Length<th>Sequence";
        for(const DistinctSequence& distinctSequence: distinctSequences) {
            const vector<Base>& sequence = distinctSequence.sequence();
            html <<
                "<tr>"
                "<td class=centered>" << distinctSequence.coverage <<
                "<td class=centered>" << sequence.size() <<
                "<td class=left style='font-family:monospace;white-space:nowrap''>";
            std::ranges::copy(sequence, ostream_iterator<Base>(html));
        }
        html << "</table>";
    }



    // If there is a dominant sequence, use it as the consensus.
    const uint64_t maximumCoverage = distinctSequences.front().coverage;
    const uint64_t totalCoverage = edge.extendedInfos.size();
    if(maximumCoverage > totalCoverage/2) {
        const vector<Base>& dominantSequence = distinctSequences.front().sequence();
        std::ranges::copy(dominantSequence, back_inserter(sequence));
        for(uint64_t i=0; i<dominantSequence.size(); i++) {
            coverage.push_back(maximumCoverage);
        }
        if(html) {
            html << "<br>The dominant sequence with coverage " << maximumCoverage <<
                " was used as the consensus.";
            if(maximumCoverage < totalCoverage) {
                html << "<br>Stored base coverages are lower bounds.";
            }
        }
        const uint64_t sequenceEnd = sequence.size();
        html << "<br>This edge assembled positions [" << sequenceBegin <<
            "," << sequenceEnd <<
            ") of the sequence of this local assembly.";
        return;
    }



    // If getting here, we have to run the multiple sequence alignment of the distinct sequences.
    vector< pair<vector<Base>, uint64_t> > sequencesWithCoverage;
    uint64_t maxLength = 0;
    for(const DistinctSequence& distinctSequence: distinctSequences) {
         sequencesWithCoverage.emplace_back(distinctSequence.sequence(), distinctSequence.coverage);
         maxLength = max(maxLength, distinctSequence.sequence().size());
    }
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;
    if(maxLength < abpoaMaxLength) {
        shasta2Abpoa(sequencesWithCoverage, consensus, alignment, alignedConsensus);
    } else {
        shasta2Poasta(sequencesWithCoverage, consensus, alignment, alignedConsensus);
    }

    // Store assembled sequence
    for(const auto& [base, baseCoverage]: consensus) {
        sequence.push_back(base);
        coverage.push_back(baseCoverage);
    }



    // Write the multiple sequence alignment of the distinct sequences.
    if(html) {
        html <<
            "<h4>Multiple sequence alignment of the distinct sequences</h4>"
            "<table><tr>"
            "<th>Coverage<th>Length<th>Aligned sequence";
        for(uint64_t i=0; i<distinctSequences.size(); i++) {
            const DistinctSequence& distinctSequence = distinctSequences[i];
            const vector<Base>& sequence = distinctSequence.sequence();
            const vector<AlignedBase>& alignmentRow = alignment[i];
            html <<
                "<tr>"
                "<td class=centered>" << distinctSequence.coverage <<
                "<td class=centered>" << sequence.size() <<
                "<td class=left style='font-family:monospace;white-space:nowrap'>";
            for(uint64_t i=0; i<alignmentRow.size(); i++) {
                const AlignedBase base = alignmentRow[i];
                const bool isMismatch = (base != alignedConsensus[i]);
                if(isMismatch) {
                    html << "<span style='background-color:Pink'>";
                }
                html << base;
                if(isMismatch) {
                    html << "</span>";
                }
            }
        }

        // Add a row with aligned consensus.
        html << "<tr><th>Consensus<td class=centered>" << consensus.size() <<
            "<td class=left style='font-family:monospace;white-space:nowrap'>";
        uint64_t consensusPosition = 0;
        for(uint64_t i=0; i<alignedConsensus.size(); i++) {
            const AlignedBase base = alignedConsensus[i];
            if(base.isGap()) {
                html << alignedConsensus[i];
            } else {
                const uint64_t coverage = consensus[consensusPosition].second;
                html << "<span title='Position " << consensusPosition <<
                    ", coverage " << coverage << "'";
                if(coverage != edge.extendedInfos.size()) {
                    html << " style='background-color:pink'";
                }
                html << ">";
                html << alignedConsensus[i];
                html << "</span>";
                ++consensusPosition;
            }
        }
        html << "</table>";



        // Write consensus and its coverage.
        html <<
            "<h4>Consensus (" << consensus.size() << " bases)</h4>"
            "<table>";

        // Consensus.
        html <<
            "<tr><th class=left>Consensus<td class=left style='font-family:monospace;white-space:nowrap'>";
        for(uint64_t position=0; position<consensus.size(); position++) {
            html << "<span title='Position " << position <<
                ", coverage " << consensus[position].second << "'";
            if(consensus[position].second != edge.extendedInfos.size()) {
                html << " style='background-color:Pink'";
            }
            html << ">";
            html << consensus[position].first;
            html << "</span>";
        }

        // Consensus coverage.
        std::map<uint64_t, char> coverageLegend;
        html <<
            "<tr><th class=left>Coverage<td class=left style='font-family:monospace;white-space:nowrap'>";
        for(uint64_t position=0; position<consensus.size(); position++) {
            const uint64_t coverage = consensus[position].second;
            char c = '*';
            if(coverage < 10) {
                c = char(coverage) + '0';
            } else if(coverage < 36) {
                c = (char(coverage) - char(10)) + 'A';
            }
            coverageLegend[coverage] = c;
            html << "<span title='Position " << position <<
                ", coverage " << consensus[position].second << "'";
            if(consensus[position].second != edge.extendedInfos.size()) {
                html << " style='background-color:Pink'";
            }
            html << "'>";
            html << c;
            html << "</span>";
        }
        html << "</table>\n";

        // Coverage legend.
        html << "<br><table><tr><th>Coverage<th>Symbol";
        for(const auto&[coverage, character]: coverageLegend) {
            html << "<tr><td class=centered>" << coverage <<
                "<td class=centered>" << character;
        }
        html << "</table>";


    }

    const uint64_t sequenceEnd = sequence.size();
    html << "<br>This edge assembled positions [" << sequenceBegin <<
        "," << sequenceEnd <<
        ") of the sequence of this local assembly.";
}






void Shasta2LocalAssembly5::writeAssembledSequence() const
{
    html <<
        "<h4>Assembled sequence</h4>";

    html <<
        "<table>"
        "<tr><th class=left>Consensus sequence length<td class=left>" << sequence.size() <<
        "<tr><th class=left>Consensus sequence"
        "<td style='font-family:monospace;white-space:nowrap''>";

    for(uint64_t position=0; position<sequence.size(); position++) {
        const Base b = sequence[position];
        html << "<span title='" << position << "'>" << b << "</span>";
    }

    html <<
        "<tr><th class=left >Coverage"
        "<td style='font-family:monospace;white-space:nowrap''>";

    std::map<uint64_t, char> coverageLegend;

    for(uint64_t position=0; position<sequence.size(); position++) {
        const uint64_t coverageThisPosition = coverage[position];
        const char c = (coverageThisPosition < 10) ? char(coverageThisPosition + '0') : char(coverageThisPosition - 10 + 'A');
        coverageLegend.insert({coverageThisPosition, c});

        html << c;

    }

    html << "</table><br><br>";

    // Write the coverage legend.
    html << "<p><table><tr><th>Coverage<th>Symbol";
    for(const auto& p: coverageLegend) {
        html << "<tr><td class=centered>" << p.first << "<td class=centered>" << p.second;
    }
    html << "</table>";


}
