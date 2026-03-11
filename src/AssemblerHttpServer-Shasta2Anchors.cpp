#include "Assembler.hpp"
#include "Shasta2AssemblyGraphPostprocessor.hpp"
#include "Shasta2AnchorGraph.hpp"
#include "Shasta2AnchorPair.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "Shasta2LocalAnchorGraph.hpp"
#include "Shasta2LocalReadAnchorGraph.hpp"
#include "graphvizToHtml.hpp"
#include "deduplicate.hpp"
#include "LocalReadGraph.hpp"
#include "platformDependent.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/tokenizer.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

using namespace dinara;
using namespace std;

namespace {

void accessShasta2HttpData(Assembler& assembler)
{
    const MappedMemoryOwner shasta2Owner = assembler.shasta2MappedMemoryOwner();

    if(!assembler.shasta2Anchors) {
        assembler.shasta2Anchors = make_shared<Shasta2Anchors>(
            "",
            shasta2Owner,
            assembler.getReads(),
            assembler.assemblerInfo->k,
            *assembler.markers);
    }

    if(!assembler.shasta2Journeys) {
        assembler.shasta2Journeys = make_shared<Shasta2Journeys>(shasta2Owner);
    }

    if(!assembler.shasta2AnchorGraph) {
        try {
            assembler.shasta2AnchorGraph = make_shared<Shasta2AnchorGraph>(
                shasta2Owner,
                "Shasta2AnchorGraphAfterTransitiveReduction");
        } catch(const exception&) {
            try {
                assembler.shasta2AnchorGraph = make_shared<Shasta2AnchorGraph>(
                    shasta2Owner,
                    "Shasta2AnchorGraph");
            } catch(const exception&) {
                assembler.shasta2AnchorGraph.reset();
            }
        }
    }
}

} // namespace

void Assembler::exploreShasta2Anchor(const vector<string>& request, ostream& html)
{
    try {
        accessShasta2HttpData(*this);
    } catch(const exception&) {
    }
    if(!shasta2Anchors || !shasta2Journeys) {
        html << "<p>Shasta2 anchors/journeys not available.";
        return;
    }
    const Shasta2Anchors& anchors = *shasta2Anchors;
    const Shasta2Journeys& journeys = *shasta2Journeys;

    const uint64_t k = assemblerInfo->k;

    // Get the request parameters.
    string anchorIdString;
    const bool anchorIdStringIsPresent = HttpServer::getParameterValue(request, "anchorIdString", anchorIdString);
    boost::trim(anchorIdString);
    if(not anchorIdStringIsPresent) {
        uint64_t anchorId = invalid<uint64_t>;
        if(HttpServer::getParameterValue(request, "anchorId", anchorId)) {
            anchorIdString = to_string(anchorId);
        }
    }

    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    // Begin the form.
    html <<
        "<h2>Shasta2Anchor information</h2>"
        "<form>"
        "<table>"

    // Shasta2AnchorId
        "<tr><th class=left>Shasta2Anchor id<td class=centered>"
        "<input type=text name=anchorIdString required style='text-align:center'";
    if(anchorIdStringIsPresent) {
        html << " value='" << anchorIdString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors.size() - 1 << " followed by + or -.'>";

    // Assembly stage for annotations.
    html <<
        "<tr title='Leave blank for no annotations'>"
        "<th class=left>Assembly stage for annotations"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center'";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=10>";

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Show anchor details'> "
        "</form>";



    // If the anchor id missing or invalid, stop here.
    if(not anchorIdStringIsPresent) {
        return;
    }
    const Shasta2AnchorId anchorId = shasta2AnchorIdFromString(anchorIdString);

    if((anchorId == invalid<Shasta2AnchorId>) or (anchorId >= anchors.size())) {
        html << "<p>Invalid anchor id. Must be a number between 0 and " <<
            anchors.size() - 1 << " followed by + or -.";
        return;
    }


    html << "<h1>Shasta2Anchor " << anchorIdString << "</h1>";

    const auto markerInfos = anchors[anchorId];
    const uint64_t coverage = markerInfos.size();
    const vector<Base> kmerSequence = anchors.anchorKmerSequence(anchorId);


    vector<Shasta2AnchorId> parents;
    vector<uint64_t> parentsCoverage;
    anchors.findParents(journeys, anchorId, parents, parentsCoverage);

    vector<Shasta2AnchorId> children;
    vector<uint64_t> childrenCoverage;
    anchors.findChildren(journeys, anchorId, children, childrenCoverage);

    // Write a summary table.
    html <<
        "<table>"
        "<tr><th class=left>Coverage<td class=centered>" << coverage <<
        "<tr><th class=left>K-mer sequence<td class=centered style='font-family:monospace'>";
    copy(kmerSequence.begin(), kmerSequence.end(), ostream_iterator<Base>(html));

    html <<
        "<tr><th class=left>Parent anchors"
        "<td class=centered>";
    for(uint64_t i=0; i<parents.size(); i++) {
        const string parentAnchorIdString = shasta2AnchorIdToString(parents[i]);
        if(i != 0) {
            html << "<br>";
        }
        html <<
            "<a href='exploreShasta2Anchor?anchorIdString=" << HttpServer::urlEncode(parentAnchorIdString) << "'>" <<
            parentAnchorIdString << "</a> coverage " << parentsCoverage[i] <<
            " <a href='exploreShasta2AnchorPair2?anchorIdAString=" << HttpServer::urlEncode(parentAnchorIdString) <<
            "&anchorIdBString=" << HttpServer::urlEncode(anchorIdString) <<
            "&adjacentInJourney=on'>pair</a>";

    }

    html <<
        "<tr><th class=left>Children anchors"
        "<td class=centered>";
    for(uint64_t i=0; i<children.size(); i++) {
        const string childAnchorIdString = shasta2AnchorIdToString(children[i]);
        if(i != 0) {
            html << "<br>";
        }
        html <<
            "<a href='exploreShasta2Anchor?anchorIdString=" << HttpServer::urlEncode(childAnchorIdString) << "'>" <<
            childAnchorIdString << "</a> coverage " << childrenCoverage[i] <<
            " <a href='exploreShasta2AnchorPair2?anchorIdAString=" << HttpServer::urlEncode(anchorIdString) <<
            "&anchorIdBString=" << HttpServer::urlEncode(childAnchorIdString) <<
            "&adjacentInJourney=on'>pair</a>";

    }
    html << "</table>";



    // Assembly graph annotations, if requested.
    if(not assemblyStage.empty()) {

        const Shasta2AssemblyGraphOptions options;
        const Shasta2AssemblyGraphPostprocessor& assemblyGraph =
            getShasta2AssemblyGraph(assemblyStage, options);
        const auto annotations = assemblyGraph.getAnnotations(anchorId);

        html << "<h2>Assembly graph annotations at assembly stage " << assemblyStage << "</h2>";

        if(annotations.empty()) {
            html << "This Shasta2Anchor is not referenced in assembly stage " << assemblyStage;
        } else {
            html << "<ul>";

            for(const auto& annotation: annotations) {
                html << "<li>";

                if(annotation.v == Shasta2AssemblyGraph::null_vertex()) {
                    // This Shasta2AnchorId is used in a step.
                    const Shasta2AssemblyGraphEdge& edge = assemblyGraph[annotation.e];
                    const string segmentUrl = "exploreShasta2SegmentSteps?assemblyStage=" + assemblyStage +
                        "&segmentName=" + to_string(edge.id);
                    const string stepUrl = "exploreShasta2SegmentStep?assemblyStage=" + assemblyStage +
                        "&segmentName=" + to_string(edge.id) + "&stepId=" + to_string(annotation.step);
                    html <<
                        "Segment <a href='" << segmentUrl << "'>" << edge.id << "</a>"
                        ", step <a href='" << stepUrl << "'>" << annotation.step << "</a>"
                        " of " << edge.size() <<
                        ", " <<
                        (annotation.isAnchorIdA ? "first" : "second") <<
                        " anchor.";

                } else {

                    // This Shasta2AnchorId is used in a vertex.
                    const Shasta2AssemblyGraph::vertex_descriptor v = annotation.v;
                    html << "Assembly graph vertex with";

                    // Incoming segments.
                    if(in_degree(v, assemblyGraph) == 0) {
                        html << " no incoming segments";
                    } else {
                        html << " incoming segments";
                        BGL_FORALL_INEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
                            const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
                            const string segmentUrl = "exploreShasta2SegmentSteps?assemblyStage=" + assemblyStage +
                                "&segmentName=" + to_string(edge.id);
                            html << " <a href='" << segmentUrl << "'>" << edge.id << "</a>";
                        }
                        html << ",";
                    }

                    // Outgoing segments.
                    if(out_degree(v, assemblyGraph) == 0) {
                        html << " no outgoing segments";
                    } else {
                        html << " outgoing segments";
                        BGL_FORALL_OUTEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
                            const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
                            const string segmentUrl = "exploreShasta2SegmentSteps?assemblyStage=" + assemblyStage +
                                "&segmentName=" + to_string(edge.id);
                            html << " <a href='" << segmentUrl << "'>" << edge.id << "</a>";
                        }
                        html << ".";
                    }
                }
            }

            html << "</ul>";
        }
    }



    std::map<pair<Shasta2AnchorId, Shasta2AnchorId>, uint64_t> tangleMatrix;

    // Write the marker intervals of this Shasta2Anchor.
    html <<
        "<h2>Marker intervals</h2>"
        "<table>"
        "<tr>"
        "<th>Index"
        "<th>Oriented<br>read<br>id"
        "<th>Position<br>in<br>journey"
        "<th>Ordinal"
        "<th>Position"
        "<th>Journey"
        "<th>Previous<br>anchor<br>in journey"
        "<th>Next<br>anchor<br>in journey";

    // Loop over the marker intervals.
    for(uint64_t i=0; i<coverage; i++) {
        const Shasta2AnchorMarkerInfo& markerInfo = markerInfos[i];
        const OrientedReadId orientedReadId = markerInfo.orientedReadId;
        const auto journey = journeys[orientedReadId];

        const uint32_t ordinal = markerInfo.ordinal;

        const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
        const uint32_t position = orientedReadMarkers[ordinal].position;

        Shasta2AnchorId previousAnchorInJourney = invalid<Shasta2AnchorId>;
        if(markerInfo.positionInJourney > 0) {
            previousAnchorInJourney = journey[markerInfo.positionInJourney - 1];
        }
        Shasta2AnchorId nextAnchorInJourney = invalid<Shasta2AnchorId>;
        if(markerInfo.positionInJourney < journey.size() - 1) {
            nextAnchorInJourney = journey[markerInfo.positionInJourney + 1];
        }

        html <<
            "<tr>"
            "<td class=centered>" << i;

        // The OrientedReadId is written with an hyperlink that will
        // display the portion of the oriented read around this Shasta2Anchor.
        const string url =
            "exploreReadSequence?"
            "readId=" + to_string(orientedReadId.getReadId()) +
            "&strand=" + to_string(orientedReadId.getStrand()) +
            "&beginPosition=" + to_string((position > 2 * k) ? (position - 2 * k) : 0) +
            "&endPosition=" + to_string(position + 3 * k - 1);
        const string journeyUrl =
            "exploreShasta2Journey?"
            "readId=" + to_string(orientedReadId.getReadId()) +
            "&strand=" + to_string(orientedReadId.getStrand()) +
            "&beginPosition=" + to_string((position > 10 * k) ? (position - 10 * k) : 0) +
            "&endPosition=" + to_string(position + 11 * k);
        html <<
            "<td class=centered>" <<
            "<a href='" << url << "'>" <<
            orientedReadId << "</a>";

       html <<
            "<td class=centered>" << markerInfo.positionInJourney <<
            "<td class=centered>" << ordinal <<
            "<td class=centered>" << position <<
            "<td class=centered><a href='" << journeyUrl << "'>journey</a>";

       // Previous anchor in journey.
       html << "<td class=centered>";
       if(previousAnchorInJourney != invalid<Shasta2AnchorId>) {
           const string previousAnchorIdString = shasta2AnchorIdToString(previousAnchorInJourney);
           html << "<a href='exploreShasta2Anchor?anchorIdString=" <<
               HttpServer::urlEncode(previousAnchorIdString) << "'>" <<
               previousAnchorIdString << "</a>";
       }

       // Next anchor in journey.
       html << "<td class=centered>";
       if(nextAnchorInJourney != invalid<Shasta2AnchorId>) {
           const string nextAnchorIdString = shasta2AnchorIdToString(nextAnchorInJourney);
           html << "<a href='exploreShasta2Anchor?anchorIdString=" <<
               HttpServer::urlEncode(nextAnchorIdString) << "'>" <<
               nextAnchorIdString << "</a>";
       }

       auto it = tangleMatrix.find(make_pair(previousAnchorInJourney, nextAnchorInJourney));
       if(it == tangleMatrix.end()) {
           tangleMatrix.insert(make_pair(make_pair(previousAnchorInJourney, nextAnchorInJourney), 1));
       } else {
           ++(it->second);
       }
    }
    html << "</table>";



    html << "<h2>Tangle matrix</h2>";
    html << "<table><tr><th>In<th>Out<th>Coverage";
    for(const auto& p: tangleMatrix) {
        const Shasta2AnchorId previousAnchorInJourney = p.first.first;
        const Shasta2AnchorId nextAnchorInJourney = p.first.second;
        const uint64_t coverage = p.second;

        html << "<tr><td class=centered>";
        if(previousAnchorInJourney != invalid<Shasta2AnchorId>) {
            html << shasta2AnchorIdToString(previousAnchorInJourney);
        }

        html << "<td class=centered>";
        if(nextAnchorInJourney != invalid<Shasta2AnchorId>) {
            html << shasta2AnchorIdToString(nextAnchorInJourney);
        }

        html << "<td class=centered>" << coverage;
    }
}



void Assembler::exploreShasta2AnchorPair2(const vector<string>& request, ostream& html)
{
    try {
        accessShasta2HttpData(*this);
    } catch(const exception&) {
    }
    if(!shasta2Anchors || !shasta2Journeys) {
        html << "<p>Shasta2 anchors/journeys not available.";
        return;
    }
    const Shasta2Anchors& anchors = *shasta2Anchors;
    const Shasta2Journeys& journeys = *shasta2Journeys;

    // Get the parameters for the request
    string anchorIdAString;
    const bool anchorIdAStringIsPresent = HttpServer::getParameterValue(request, "anchorIdAString", anchorIdAString);
    boost::trim(anchorIdAString);
    string anchorIdBString;
    const bool anchorIdBStringIsPresent = HttpServer::getParameterValue(request, "anchorIdBString", anchorIdBString);
    boost::trim(anchorIdBString);

    string adjacentInJourneyString;
    const bool adjacentInJourney = HttpServer::getParameterValue(request,
        "adjacentInJourney", adjacentInJourneyString);



    // Write the form.
    html << "<form><table>";

    html <<
        "<tr><th class=left>Shasta2Anchor A"
        "<td class=centered><input type=text name=anchorIdAString required style='text-align:center'";
    if(anchorIdAStringIsPresent) {
        html << " value='" << anchorIdAString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors.size() - 1 << " followed by + or -.'><br>";

    html <<
        "<tr><th class=left>Shasta2Anchor B"
        "<td class=centered><input type=text name=anchorIdBString required style='text-align:center'";
    if(anchorIdBStringIsPresent) {
        html << " value='" << anchorIdBString + "'";
    }
    html <<
        " size=8 title='Enter an anchor id between 0 and " <<
        anchors.size() - 1 << " followed by + or -.'><br>"

        "<tr><th>Adjacent in journey"
        "<td class=centered><input type=checkbox name=adjacentInJourney" <<
        (adjacentInJourney ? " checked" : "") <<
        ">"

        "</table>"
        "<input type=submit value='Get anchor pair information'>"
        "</form>";



    // Check the AnchorIds
    if(not (anchorIdAStringIsPresent and anchorIdBStringIsPresent)) {
        return;
    }
    const Shasta2AnchorId anchorIdA = shasta2AnchorIdFromString(anchorIdAString);
    const Shasta2AnchorId anchorIdB = shasta2AnchorIdFromString(anchorIdBString);

    if((anchorIdA == invalid<Shasta2AnchorId>) or (anchorIdA >= anchors.size())) {
        html << "<p>Invalid anchor id " << anchorIdAString << ". Must be a number between 0 and " <<
            anchors.size() - 1 << " followed by + or -.";
        return;
    }

    if((anchorIdB == invalid<Shasta2AnchorId>) or (anchorIdB >= anchors.size())) {
        html << "<p>Invalid anchor id " << anchorIdBString << " .Must be a number between 0 and " <<
            anchors.size() - 1 << " followed by + or -.";
        return;
    }

    if(anchorIdA == anchorIdB) {
        html << "Specify two distinct anchors.";
        return;
    }



    // Create a Shasta2AnchorPair from these two anchors.
    const Shasta2AnchorPair anchorPair(anchors, anchorIdA, anchorIdB, adjacentInJourney);

    // Output to html.
    html << "<h2>Shasta2Anchor pair</h2>";
    anchorPair.writeAllHtml(html, anchors, journeys);
}



void Assembler::exploreShasta2Journey(const vector<string>& request, ostream& html)
{
    try {
        accessShasta2HttpData(*this);
    } catch(const exception&) {
    }
    if(!shasta2Anchors || !shasta2Journeys) {
        html << "<p>Shasta2 anchors/journeys not available.";
        return;
    }
    const Shasta2Anchors& anchors = *shasta2Anchors;
    const Shasta2Journeys& journeys = *shasta2Journeys;

    html <<
        "<h2>Oriented read journey</h2>"
        "<p>The journey of an oriented read is the sequence of anchors it visits.";


    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = HttpServer::getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = HttpServer::getParameterValue(request, "strand", strand);
    uint32_t beginPosition = 0;
    const bool beginPositionIsPresent = HttpServer::getParameterValue(request, "beginPosition", beginPosition);
    uint32_t endPosition = 0;
    const bool endPositionIsPresent = HttpServer::getParameterValue(request, "endPosition", endPosition);

    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<th class=left>Numeric read id"
        "<td><input type=text name=readId" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " title='Enter a read id between 0 and " << anchors.reads.readCount()-1 << "'>"

        "<tr>"
        "<th class=left>Strand"
        "<td>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);

    html <<
        "<tr>"
        "<th class=left>Begin position"
        "<td><input type=text name=beginPosition"
        " title='Leave blank to begin display at beginning of read.'";
    if(beginPositionIsPresent) {
        html << " value=" << beginPosition;
    }
    html << ">";

    html <<
        "<tr>"
        "<th class=left>End position"
        "<td><input type=text name=endPosition"
        " title='Leave blank to end display at end of read.'";
    if(endPositionIsPresent) {
        html << " value=" << endPosition;
    }
    html << ">";


    html <<
        "</table>"
        "<input type=submit value='Display'>"
        "</form>";

    if(not readIdIsPresent) {
        html << "Specify a numeric read id.";
        return;
    }

    // If the strand is missing, stop here.
    if(not strandIsPresent) {
        return;
    }

    // Sanity checks.
    if(readId >= anchors.reads.readCount()) {
        html << "<p>Invalid read id.";
        return;
    }
    if(strand!=0 && strand!=1) {
        html << "<p>Invalid strand.";
        return;
    }

    // Adjust the position range, if necessary.
    if(!beginPositionIsPresent) {
        beginPosition = 0;
    }
    if(!endPositionIsPresent) {
        endPosition = uint32_t(anchors.reads.getReadRawSequenceLength(readId));
    }
    if(endPosition <= beginPosition) {
        html << "<p>Invalid choice of begin and end position.";
        return;
    }

    // Access the information we need.
    const OrientedReadId orientedReadId(readId, strand);
    DINARA_ASSERT(journeys.isOpen());
    DINARA_ASSERT(journeys.size() == 2 * anchors.reads.readCount());
    const span<const Shasta2AnchorId> journey = journeys[orientedReadId];

    // Page title.
    html << "<h2>Shasta2Journey of oriented read " << orientedReadId << "</h2>";

    // Begin the main table.
    html <<
        "<table><tr>"
        "<th>Position<br>in journey"
        "<th>Shasta2Anchor"
        "<th>Shasta2Anchor<br>coverage"
        "<th>Marker<br>ordinal"
        "<th>Marker<br>position"
        "<th>Previous<br>anchor"
        "<th>Next<br>anchor";

    // Loop over the anchors in the journey of this oriented read.
    for(uint64_t positionInJourney=0; positionInJourney<journey.size(); positionInJourney++) {
        const Shasta2AnchorId anchorId = journey[positionInJourney];
        const uint64_t anchorCoverage = anchors[anchorId].coverage();
        const string anchorIdString = shasta2AnchorIdToString(anchorId);

        const uint64_t ordinal = anchors.getOrdinal(anchorId, orientedReadId);

        const auto orientedReadMarkers = anchors.markers[orientedReadId.getValue()];
        const uint32_t position = orientedReadMarkers[ordinal].position;
        const Shasta2AnchorId previousAnchorId =
            (positionInJourney > 0) ? journey[positionInJourney - 1] : invalid<Shasta2AnchorId>;
        const Shasta2AnchorId nextAnchorId =
            ((positionInJourney + 1) < journey.size()) ? journey[positionInJourney + 1] : invalid<Shasta2AnchorId>;

        if(position < beginPosition) {
            continue;
        }
        if(position >= endPosition) {
            continue;
        }

        html <<
            "<tr>"
            "<td class=centered>" << positionInJourney <<
            "<td class=centered>" <<
            "<a href='exploreShasta2Anchor?anchorIdString=" << HttpServer::urlEncode(anchorIdString) << "'>" <<
            anchorIdString << "</a>"
            "<td class=centered>" << anchorCoverage <<
            "<td class=centered>" << ordinal <<
            "<td class=centered>" << position;

        html << "<td class=centered>";
        if(previousAnchorId != invalid<Shasta2AnchorId>) {
            const string previousAnchorIdString = shasta2AnchorIdToString(previousAnchorId);
            html << "<a href='exploreShasta2Anchor?anchorIdString=" <<
                HttpServer::urlEncode(previousAnchorIdString) << "'>" <<
                previousAnchorIdString << "</a>";
        }

        html << "<td class=centered>";
        if(nextAnchorId != invalid<Shasta2AnchorId>) {
            const string nextAnchorIdString = shasta2AnchorIdToString(nextAnchorId);
            html << "<a href='exploreShasta2Anchor?anchorIdString=" <<
                HttpServer::urlEncode(nextAnchorIdString) << "'>" <<
                nextAnchorIdString << "</a>";
        }
    }

    html << "</table>";
}



void Assembler::exploreShasta2LocalAnchorGraph(
    const vector<string>& request,
    ostream& html)
{
    try {
        accessShasta2HttpData(*this);
    } catch(const exception&) {
    }
    if(!shasta2Anchors || !shasta2Journeys || !shasta2AnchorGraph) {
        html << "<p>Shasta2 anchors/journeys/anchor graph not available.";
        return;
    }
    const Shasta2Anchors& anchors = *shasta2Anchors;
    const Shasta2Journeys& journeys = *shasta2Journeys;
    static_cast<void>(journeys); // Not currently needed here, but kept for upstream parity.

    // Get the options that control graph creation.
    string anchorIdsString;
    HttpServer::getParameterValue(request, "anchorIdsString", anchorIdsString);
    boost::trim(anchorIdsString);

    uint64_t distance = 10;
    HttpServer::getParameterValue(request, "distance", distance);

    uint64_t minCoverage = 1;
    HttpServer::getParameterValue(request, "minCoverage", minCoverage);

    string includeEdgesNotMarkedForAssemblyString;
    bool includeEdgesNotMarkedForAssembly = HttpServer::getParameterValue(request,
        "includeEdgesNotMarkedForAssembly", includeEdgesNotMarkedForAssemblyString);


    // Get the options that control graph display.
    const Shasta2LocalAnchorGraphDisplayOptions displayOptions(request);



    // Start the form.
    html << "<form><table>";

    // Form items for options that control graph creation.
    html <<
        "<tr title='Enter comma or space separated anchor ids, each a number between 0 and " <<
        anchors.size() - 1 << ".'>"
        "<th class=left>Starting anchor ids"
        "<td class=centered><input type=text name=anchorIdsString style='text-align:center' required";
    if(not anchorIdsString.empty()) {
        html << " value='" << anchorIdsString + "'";
    }
    html <<
        " size=8 title='Enter comma separated anchor ids, each between 0 and " <<
        anchors.size() - 1 << ".'>";

    html << "<tr>"
        "<th class=left>Distance"
        "<td class=centered>"
        "<input type=text name=distance style='text-align:center' required size=8 value=" <<
        distance << ">";

    html << "<tr>"
        "<th class=left>Minimum coverage"
        "<td class=centered>"
        "<input type=text name=minCoverage style='text-align:center' required size=8 value=" <<
        minCoverage << ">";

    html << "<tr>"
        "<th class=left>Include edges not marked for use in assembly"
        "<td class=centered>"
        "<input type=checkbox name=includeEdgesNotMarkedForAssembly" <<
        (includeEdgesNotMarkedForAssembly ? " checked" : "") <<
        ">";

    // Form items for options that control graph display.
    displayOptions.writeForm(html);

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Create local anchor graph'>"
        "</form>";



    // If the anchor id are missing, stop here.
    if(anchorIdsString.empty()) {
        return;
    }


    // Extract the AnchorIds.
    vector<Shasta2AnchorId> anchorIds;
    boost::tokenizer< boost::char_separator<char> > tokenizer(anchorIdsString, boost::char_separator<char>(", "));

    for(const string& anchorIdString: tokenizer) {
        const Shasta2AnchorId anchorId = shasta2AnchorIdFromString(anchorIdString);

        if((anchorId == invalid<Shasta2AnchorId>) or (anchorId >= anchors.size())) {
            html << "<p>Invalid anchor id " << anchorIdString << ". Must be a number between 0 and " <<
                anchors.size() - 1 << ".";
            return;
        }

        anchorIds.push_back(anchorId);
    }
    deduplicate(anchorIds);

    // Access the Shasta2AnchorGraph we are going to use.
    const Shasta2AnchorGraph& anchorGraph = *shasta2AnchorGraph;


    // If needed, get the AssemblyGraph for this assembly stage.
    const Shasta2AssemblyGraphPostprocessor* assemblyGraphPointer = 0;
    if(displayOptions.vertexColoring == "byAssemblyAnnotations") {
        const Shasta2AssemblyGraphOptions options;
        const Shasta2AssemblyGraphPostprocessor& assemblyGraph =
            getShasta2AssemblyGraph(displayOptions.assemblyStage, options);
        assemblyGraphPointer = &assemblyGraph;
    }



    // Create the Shasta2LocalAnchorGraph starting from these AnchorIds and moving
    // away up to the specified distance.
    Shasta2LocalAnchorGraph graph(
        anchors,
        anchorGraph,
        anchorIds,
        distance,
		minCoverage,
        not includeEdgesNotMarkedForAssembly);

    html << "<h1>Local anchor graph</h1>";
    html << "The local anchor graph has " << num_vertices(graph) <<
         " vertices and " << num_edges(graph) << " edges.";

    // Write it to html.
    graph.writeHtml(html, displayOptions, assemblyGraphPointer);

}



void Assembler::exploreShasta2LocalReadAnchorGraph(
    const vector<string>& request,
    ostream& html)
{
    try {
        accessShasta2HttpData(*this);
    } catch(const exception&) {
    }
    if(!shasta2Anchors || !shasta2Journeys) {
        html << "<p>Shasta2 anchors/journeys not available.";
        return;
    }
    const Shasta2Anchors& anchors = *shasta2Anchors;
    const Shasta2Journeys& journeys = *shasta2Journeys;
    Shasta2LocalReadAnchorGraph graph(anchors, journeys, request, html);
}



void Assembler::exploreShasta2LocalReadGraph(
    const vector<string>& request,
    ostream& html)
{
    try {
        checkReadGraphIsOpen();
    } catch(const exception& e) {
        html << "<p>" << e.what();
        return;
    }

    // Reuse the standard read graph viewer.
    exploreReadGraph(request, html);
}
