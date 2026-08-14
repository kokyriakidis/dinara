// Dinara.
#include "Assembler.hpp"
#include "compressAlignment.hpp"
#include "ConsensusCaller.hpp"
#include "Coverage.hpp"
#include "hsv.hpp"
#include "InducedAlignment.hpp"
#include "MarkerConnectivityGraph.hpp"
#include "MurmurHash2.hpp"
#include "platformDependent.hpp"
#include "Reads.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"
#include "filesystem.hpp"
#include "fstream.hpp"
#include "iterator.hpp"
#include <queue>




void Assembler::exploreMarkerGraphVertex(const vector<string>& request, ostream& html)
{
    // Get the vertex id.
    MarkerGraph::VertexId vertexId = 0;
    const bool vertexIdIsPresent = getParameterValue(request, "vertexId", vertexId);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show details for marker graph vertex'> "
        "<input type=text name=vertexId required" <<
        (vertexIdIsPresent ? (" value=" + to_string(vertexId)) : "") <<
        " size=8 title='Enter a vertex id between 0 and " << markerGraph.vertexCount()-1 << "'>";
    html << "</form>";

    // If the vertex id missing or invalid, stop here.
    if(!vertexIdIsPresent || !vertexIdIsPresent) {
        return;
    }
    if(vertexId >= markerGraph.vertexCount()) {
        html << "<p>Invalid vertex id. Must be less than " << markerGraph.vertexCount() << ".";
        return;
    }

    // Access the markers of this vertex.
    span<MarkerId> markerIds = markerGraph.getVertexMarkerIds(vertexId);
    const size_t markerCount = markerIds.size();
    DINARA_ASSERT(markerCount > 0);

    // Get the marker sequence.
    const KmerId kmerId = getMarkerGraphVertexKmerId(vertexId);
    const size_t k = assemblerInfo->k;
    const Kmer kmer(kmerId, k);



    // Extract the information we need.
    vector<OrientedReadId> orientedReadIds(markerCount);
    vector<uint32_t> ordinals(markerCount);
    vector< vector<uint8_t> > repeatCounts(markerCount, vector<uint8_t>(k));
    for(size_t j=0; j<markerCount; j++) {
        const MarkerId markerId = markerIds[j];
        const CompressedMarker& marker = markers->begin()[markerId];
        tie(orientedReadIds[j], ordinals[j]) = findMarkerId(markerId);

        // Get the repeat count for this marker at each of the k positions.
        for(size_t i=0; i<k; i++) {
            Base base;
            if(assemblerInfo->readRepresentation == 1) {
                tie(base, repeatCounts[j][i]) =
                    reads->getOrientedReadBaseAndRepeatCount(orientedReadIds[j], uint32_t(marker.position+i));
            } else {
                base = reads->getOrientedReadBase(orientedReadIds[j], uint32_t(marker.position+i));
                repeatCounts[j][i] = 1;
            }
            DINARA_ASSERT(base == kmer[i]);
        }
    }



    // Find all the repeat counts represented.
    std::set<size_t> repeatCountsSet;
    for(const auto& v: repeatCounts) {
        for(const auto r: v) {
            repeatCountsSet.insert(r);
        }
    }



    // Compute consensus repeat counts at each of the k positions.
    const bool consensusIsAvailable = markerGraph.vertexRepeatCounts.isOpen;
    vector<size_t> consensusRepeatCounts(k);
    if(consensusIsAvailable) {
        const auto storedConsensusRepeatCounts =
            markerGraph.vertexRepeatCounts.begin() + k * vertexId;
        for(size_t i=0; i<k; i++) {

            Coverage coverage;
            for(size_t j=0; j<markerCount; j++) {
                coverage.addRead(
                    AlignedBase(kmer[i]),
                    orientedReadIds[j].getStrand(),
                    repeatCounts[j][i]);
            }

            const Consensus consensus = (*consensusCaller)(coverage);
            DINARA_ASSERT(Base(consensus.base) == kmer[i]);
            consensusRepeatCounts[i] = consensus.repeatCount;

            // Check that this repeat count agrees with what was
            // computed during the assembly.
            if(consensusRepeatCounts[i] != storedConsensusRepeatCounts[i]) {
                html << "<p><b>Stored consensus repeat counts do not agree with "
                    "the values computed on the fly.</b>" << endl;
            }
        }
    }




    // Compute concordant and discordant coverage at each position.
    vector<size_t> concordantCoverage(k, 0);
    vector<size_t> discordantCoverage(k, 0);
    if(consensusIsAvailable) {
        for(size_t i=0; i<k; i++) {
            for(size_t j=0; j<markerCount; j++) {
                if(repeatCounts[j][i] == consensusRepeatCounts[i]) {
                    ++concordantCoverage[i];
                } else {
                    ++discordantCoverage[i];
                }
            }
        }
    }


    // Page title.
    const string titleUrl =
        "exploreMarkerGraph0?vertexId=" + to_string(vertexId) +
        "&maxDistance=3"
        "&useWeakEdges=on"
        "&usePrunedEdges=on"
        "&useSuperBubbleEdges=on"
        "&useLowCoverageCrossEdges=on"
        "&sizePixels=800"
        "&timeout=30";
    html << "<h1>Marker graph vertex <a href='" << titleUrl << "'> "<< vertexId << "</a></h1>";



    // Table with summary information for this vertex.
    html <<
        "<table>"
        "<tr><th class=left>Coverage<td class=centered>" << markerCount <<
        "<tr><th class=left>Marker sequence (run-length)" <<
        "<td class=centered style='font-family:monospace'>";
    kmer.write(html, assemblerInfo->k);


    if(consensusIsAvailable) {
        // Write a row with consensus repeat counts.
        html <<
            "<tr><th class=left>Consensus repeat counts"
            "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<k; i++) {
            const size_t repeatCount = consensusRepeatCounts[i];
            if(repeatCount < 10) {
                html << repeatCount;
            } else {
                html << "*";
            }
        }

        // Write a row with the consensus raw sequence.
        html <<
            "<tr><th class=left>Consensus raw sequence"
            "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<k; i++) {
            const Base base = kmer[i];
            const size_t repeatCount = consensusRepeatCounts[i];
            for(size_t r=0; r<repeatCount; r++) {
                html << base;
            }
        }
    }



    // Write a row with outgoing edges.
    html <<
        "<tr><th class=left>Edges starting at this vertex<td class=centered>";
    for(const MarkerGraph::EdgeId edgeId: markerGraph.edgesBySource[vertexId]) {
        html << "<a href='exploreMarkerGraphEdge?edgeId=" << edgeId << "'";
        if(markerGraph.edges[edgeId].wasRemoved()) {
            html << " style='color:#a8b9ea'";
        }
        html << ">" << edgeId << "</a> ";
    }
    html <<
        "<tr><th class=left>Edges ending at this vertex<td class=centered>";
    for(const MarkerGraph::EdgeId edgeId: markerGraph.edgesByTarget[vertexId]) {
        html << "<a href='exploreMarkerGraphEdge?edgeId=" << edgeId << "'";
        if(markerGraph.edges[edgeId].wasRemoved()) {
            html << " style='color:#a8b9ea'";
        }
        html << ">" << edgeId << "</a> ";
    }
    const MarkerGraph::VertexId vertexIdRc = markerGraph.reverseComplementVertex[vertexId];
    html << "<tr><th class=left>Reverse complement vertex<td class=centered>"
        "<a href='exploreMarkerGraphVertex?vertexId=" << vertexIdRc << "'>" <<
        vertexIdRc << "</a> ";


    html << "</table>";



    // Table with one row for each marker.
    html <<
        "<h3>Markers of this vertex</h3>"
        "<table>"
        "<tr>"
        "<th>Oriented<br>read"
        "<th title='Ordinal of the marker in the oriented read'>Ordinal"
        "<th>Repeat<br>counts";
    for(size_t j=0; j<markerIds.size(); j++) {
        const OrientedReadId orientedReadId = orientedReadIds[j];
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();
        const uint32_t ordinal = ordinals[j];

        // Oriented read id.
        html <<
            "<tr>"
            "<td class=centered>"
            "<a href='exploreRead"
            "?readId=" << readId <<
            "&strand=" << strand << "'>" <<
            orientedReadId << "</a>";

        // Marker ordinal.
        html <<
            "<td class=centered>" <<
            "<a href='exploreRead"
            "?readId=" << readId <<
            "&strand=" << strand <<
            "&amp;showMarkers=on"
            "&highlightMarker=" << ordinal <<
            "'>" <<
            ordinal << "</a>";

        // Repeat counts.
        html << "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<k; i++) {
            const uint8_t repeatCount = repeatCounts[j][i];
            if(repeatCount < 10) {
                html << int(repeatCount);
            } else {
                html << "*";
            }
        }
    }


    if(consensusIsAvailable) {
        // Write a row with consensus repeat counts.
        html <<
            "<tr><th colspan=2 class=left>Consensus repeat counts"
            "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<k; i++) {
            const size_t repeatCount = consensusRepeatCounts[i];
            if(repeatCount < 10) {
                html << repeatCount;
            } else {
                html << "*";
            }
        }
    }



    // Write a row with the marker sequence (run-length).
    html <<
        "<tr><th colspan=2 class=left>Run-length sequence"
        "<td class=centered style='font-family:monospace'>";
    kmer.write(html, assemblerInfo->k);



    // Write rows with coverage information for each represented repeat value.
    if(consensusIsAvailable) {
        for(size_t repeatCount: repeatCountsSet) {
            html <<
                "<tr><th colspan=2 class=left>Coverage for repeat count " << repeatCount <<
                "<td class=centered style='font-family:monospace'>";
            for(size_t i=0; i<k; i++) {

                // Compute coverage for this repeat count, at this position.
                size_t coverage = 0;
                for(size_t j=0; j<markerCount; j++) {
                    if(repeatCounts[j][i] == repeatCount) {
                        coverage++;
                    }
                }

                // Write it out.
                if(coverage == 0) {
                    html << ".";
                } else if(coverage < 10) {
                    html << coverage;
                } else {
                    html << "*";
                }
            }
        }



        // Write a row with concordant coverage.
        html <<
            "<tr><th colspan=2 class=left>Concordant repeat count coverage"
            "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<assemblerInfo->k; i++) {
            const size_t coverage = concordantCoverage[i];
            if(coverage == 0) {
                html << ".";
            } else if(coverage < 10) {
                html << coverage;
            } else {
                html << "*";
            }
        }



        // Write a row with discordant coverage.
        html <<
            "<tr><th colspan=2 class=left>Discordant repeat count coverage"
            "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<assemblerInfo->k; i++) {
            const size_t coverage = discordantCoverage[i];
            if(coverage == 0) {
                html << ".";
            } else if(coverage < 10) {
                html << coverage;
            } else {
                html << "*";
            }
        }



        // Write a row with the consensus raw sequence.
        html <<
            "<tr><th colspan=2 class=left>Consensus raw sequence"
            "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<k; i++) {
            const Base base = kmer[i];
            const size_t repeatCount = consensusRepeatCounts[i];
            for(size_t r=0; r<repeatCount; r++) {
                html << base;
            }
        }
    }

    html << "</table>";
}






void Assembler::exploreMarkerGraphEdge(const vector<string>& request, ostream& html)
{
    // Get the edge id.
    MarkerGraph::EdgeId edgeId = 0;
    const bool edgeIdIsPresent = getParameterValue(request, "edgeId", edgeId);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show details for marker graph edge'> "
        "<input type=text name=edgeId required" <<
        (edgeIdIsPresent ? (" value=" + to_string(edgeId)) : "") <<
        " size=8 title='Enter an edge id between 0 and " << markerGraph.edges.size()-1 << "'>";
    html << "</form>";

    // If the edge id missing or invalid, stop here.
    if(!edgeIdIsPresent) {
        return;
    }
    if(edgeId >= markerGraph.edges.size()) {
        html << "<p>Invalid edge id. Must be less than " << markerGraph.edges.size() << ".";
        return;
    }

    // Access the edge.
    const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
    const array<MarkerGraph::VertexId, 2> vertexIds = {edge.source, edge.target};
    const size_t markerCount = markerGraph.edgeCoverage(edgeId);

    // The marker intervals of this edge.
    const span<MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
    DINARA_ASSERT(markerIntervals.size() == markerCount);

    html << "<h1>Marker graph edge " << edgeId << "</h1>";

    // Table to summarize this edge.
    html <<
        "<table style='display:block;white-space:nowrap;'>"
        "<tr><th class=left>Source vertex<td class=centered>"
        "<a href='exploreMarkerGraphVertex?vertexId=" << vertexIds[0] << "'>" << vertexIds[0] << "</a>"
        "<tr><th class=left>Target vertex<td class=centered>"
        "<a href='exploreMarkerGraphVertex?vertexId=" << vertexIds[1] << "'>" << vertexIds[1] << "</a>"
        "<tr><th class=left>Coverage<td class=centered>" << markerCount;

    const MarkerGraph::EdgeId edgeIdRc = markerGraph.reverseComplementEdge[edgeId];
    html << "<tr><th class=left>Reverse complement edge<td class=centered>"
        "<a href='exploreMarkerGraphEdge?edgeId=" << edgeIdRc << "'>" <<
        edgeIdRc << "</a> ";

    html << "</table>";

    // Table of the marker intervals of this edge.
    html <<
        "<h3>Marker intervals</h3>"
        "<table>"
        "<tr><th>Oriented<br>read<th>Left<br>ordinal<th>Right<br>ordinal";
    for(size_t j=0; j<markerCount; j++) {
        const MarkerInterval& markerInterval = markerIntervals[j];
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();
        html <<
            "<tr><td class=centered>"
            "<a href='exploreRead?readId=" << readId <<
            "&strand=" << strand << "'>" << orientedReadId << "</a>"
            "<td class=centered>" << markerInterval.ordinals[0] <<
            "<td class=centered>" << markerInterval.ordinals[1];
    }
    html << "</table>";
}



void Assembler::exploreMarkerGraphInducedAlignment(
    const vector<string>& request,
    ostream& html)
{
    html <<
        "<h1>Display the induced alignment matrix of two oriented reads</h1>"
        "<p>The marker graph induces an effective alignment between each pair "
        "of oriented reads which can be obtained by following each of the oriented reads "
        "in the marker graph. Aligned markers are those that are on the same vertex. "
        "The induced alignment matrix of two oriented reads <i>x</i> and <i>y</i> "
        "with <i>n<sub>x</sub></i> and <i>n<sub>y</sub></i> markers is an "
        "<i>n<sub>x</sub></i>&times;<i>n<sub>y</sub></i> matrix. "
        "Element <i>ij</i> of the matrix is 1 if marker <i>i</i> of <i>x</i> "
        "and marker <i>j</i> of <i>y</i> "
        "are on the same marker graph vertex and 0 otherwise.";

    // Get the read ids and strands from the request.
    ReadId readId0 = 0;
    const bool readId0IsPresent = getParameterValue(request, "readId0", readId0);
    Strand strand0 = 0;
    const bool strand0IsPresent = getParameterValue(request, "strand0", strand0);
    ReadId readId1 = 0;
    const bool readId1IsPresent = getParameterValue(request, "readId1", readId1);
    Strand strand1 = 0;
    const bool strand1IsPresent = getParameterValue(request, "strand1", strand1);
    string ordinalType = "ordinals";
    getParameterValue(request, "ordinalType", ordinalType);

    // Write the form.
    html <<
        "<p>Display the induced alignment matrix of these two reads:"
        "<form>"
        "<input type=text name=readId0 required size=8 " <<
        (readId0IsPresent ? "value="+to_string(readId0) : "") <<
        " title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand0", strand0IsPresent && strand0==0, strand0IsPresent && strand0==1);
    html <<
        "<br><input type=text name=readId1 required size=8 " <<
        (readId1IsPresent ? "value="+to_string(readId1) : "") <<
        " title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand1", strand1IsPresent && strand1==0, strand1IsPresent && strand1==1);

    html <<
        "<p>Plot alignment matrix using "
        "<input type=radio required name=ordinalType value='ordinals'" <<
        (ordinalType == "ordinals" ? " checked=on" : "") <<
        ">ordinals "
        "<input type=radio required name=ordinalType value='compressedOrdinals'" <<
        (ordinalType == "compressedOrdinals" ? " checked=on" : "") <<
        ">compressed ordinals"
        "<p><input type=submit value='Display induced alignment'></form>";


     // If the readId's or strand's are missing, stop here.
     if(!readId0IsPresent || !strand0IsPresent || !readId1IsPresent || !strand1IsPresent) {
         return;
     }

     // Compute the induced alignment.
     const OrientedReadId orientedReadId0(readId0, strand0);
     const OrientedReadId orientedReadId1(readId1, strand1);
     InducedAlignment inducedAlignment;
     computeInducedAlignment(orientedReadId0, orientedReadId1, inducedAlignment);
     fillCompressedOrdinals(orientedReadId0, orientedReadId1, inducedAlignment);
     html << "The induced alignment has " << inducedAlignment.data.size() <<
         " marker pairs.";

     // Write the alignment matrix to a png file.
     inducedAlignment.writePngImage(
         uint32_t(markers->size(orientedReadId0.getValue())),
         uint32_t(markers->size(orientedReadId1.getValue())),
         ordinalType == "compressedOrdinals",
         "Alignment.png");

     // Create a base64 version of the png file.
     const string command = "base64 Alignment.png > Alignment.png.base64";
     ::system(command.c_str());

     // Write the image to html.
     html <<
         "<p><img id=\"alignmentMatrix\" onmousemove=\"updateTitle(event)\" "
         "src=\"data:image/png;base64,";
         ifstream png("Alignment.png.base64");
         html << png.rdbuf();
         html << "\"/>"
             "<script>"
             "function updateTitle(e)"
             "{"
             "    var element = document.getElementById(\"alignmentMatrix\");"
             "    var rectangle = element.getBoundingClientRect();"
             "    var x = e.clientX - Math.round(rectangle.left);"
             "    var y = e.clientY - Math.round(rectangle.top);"
             "    element.title = " <<
             "\"" << orientedReadId0 << " marker \" + x + \", \" + "
             "\"" << orientedReadId1 << " marker \" + y;"
             "}"
             "</script>";



     // Write the induced alignment in a table.
     html <<
         "<p><table><tr><th>Vertex"
         "<th>Ordinal<br> in " << orientedReadId0 <<
         "<th>Ordinal<br> in " << orientedReadId1 <<
         "<th>Compressed<br>ordinal<br> in " << orientedReadId0 <<
         "<th>Compressed<br>ordinal<br> in " << orientedReadId1;
    for(const InducedAlignmentData& d: inducedAlignment.data) {
         html <<
             "<tr><td class=centered>" << d.vertexId <<
             "<td class=centered>" << d.ordinal0 <<
             "<td class=centered>" << d.ordinal1 <<
             "<td class=centered>" << d.compressedOrdinal0 <<
             "<td class=centered>" << d.compressedOrdinal1;
     }
     html << "</table>";



}



void Assembler::exploreMarkerCoverage(
    const vector<string>& request,
    ostream&html)
{
    html <<
        "<h1>Marker coverage of an oriented read</h1>"
        "<p>For each marker of a read, marker coverage is defined "
        "as the coverage (number of oriented reads) for the "
        "marker graph vertex associated with that marker, "
        "or zero if there is no marker graph vertex "
        "associated with the marker.";

    // Get the parameters from the request.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);
    uint32_t firstOrdinal = 0;
    getParameterValue(request, "firstOrdinal", firstOrdinal);
    uint32_t lastOrdinal = 0;
    getParameterValue(request, "lastOrdinal", lastOrdinal);
    int width = 600;
    getParameterValue(request, "width", width);
    int height = 400;
    getParameterValue(request, "height", height);

    // Write the form.
    html <<
        "<form><table>"
        "<tr><td>Read id<td class=centered>"
        "<input type=text name=readId required style='text-align:center'" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " size=8 title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
        "<tr><td>Strand<td class=centered>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);
    html <<
        "<tr><td>First ordinal<td class=centered>"
        "<input type=text name=firstOrdinal style='text-align:center' size=8 value='" << firstOrdinal << "'>"
        "<tr><td>Last ordinal<br>(0 for unlimited)<td class=centered>"
        "<input type=text name=lastOrdinal style='text-align:center' size=8 value='" << lastOrdinal << "'>"
        "<tr><td>Plot width<td class=centered>"
        "<input type=text name=width style='text-align:center' size=8 value='" << width << "'>"
        "<tr><td>Plot height<td class=centered>"
        "<input type=text name=height style='text-align:center' size=8 value='" << height << "'>"
        "</table><input type=submit value='Plot'></form>";

    // If the readId or strand are missing, stop here.
    if(!readIdIsPresent || !strandIsPresent) {
        return;
    }

    const OrientedReadId orientedReadId(readId, strand);
    html << "<h2>Marker coverage of oriented read " << orientedReadId << "</h2>";

    std::ostringstream gnuplotCommands;
    gnuplotCommands <<
        "set border linewidth 1\n"
        "set xtics out nomirror\n"
        "set mxtics 10\n"
        "set ytics out nomirror\n"
        "set grid xtics mxtics ytics linestyle 1 linewidth 1 linecolor rgb '#e0e0e0'\n"
        "plot '-' with points pointtype 7 pointsize 0.5 linecolor rgb '#0000ff' notitle\n";

    const uint32_t markerCount = uint32_t(markers->size(orientedReadId.getValue()));
    if(lastOrdinal == 0) {
        lastOrdinal = markerCount - 1;
    }
    DINARA_ASSERT(lastOrdinal >= firstOrdinal);
    for(uint32_t ordinal=firstOrdinal; ordinal<=lastOrdinal; ordinal++) {
        const MarkerGraph::VertexId vertexId =
            getGlobalMarkerGraphVertex(orientedReadId, ordinal);
        if(vertexId == MarkerGraph::invalidCompressedVertexId) {
            gnuplotCommands << ordinal << " " << "0\n";
        } else {
            const uint64_t coverage = markerGraph.vertexCoverage(vertexId);
            gnuplotCommands << ordinal << " " << coverage << "\n";
        }
    }

    gnuplotCommands << "e\n";
    writeGnuPlotPngToHtml(html, width, height, gnuplotCommands.str());
}



// This shows a table that follows a reads and it alignments in the marker graph.
void Assembler::followReadInMarkerGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the request parameters.
    ReadId readId0 = 0;
    const bool readId0IsPresent = getParameterValue(request, "readId", readId0);
    Strand strand0 = 0;
    const bool strand0IsPresent = getParameterValue(request, "strand", strand0);
    uint32_t firstOrdinal = 0;
    const bool firstOrdinalIsPresent = getParameterValue(request, "firstOrdinal", firstOrdinal);
    uint32_t lastOrdinal = 0;
    const bool lastOrdinalIsPresent = getParameterValue(request, "lastOrdinal", lastOrdinal);
    string whichAlignments = "ReadGraphAlignments";
    getParameterValue(request, "whichAlignments", whichAlignments);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Follow this oriented read and its alignments in the marker graph'> "
        "<br>Read <input type=text name=readId required" <<
        (readId0IsPresent ? (" value=" + to_string(readId0)) : "") <<
        " size=8 title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand", strand0IsPresent && strand0==0, strand0IsPresent && strand0==1);
    html <<
        " First marker ordinal: <input type=text name=firstOrdinal required" <<
        (firstOrdinalIsPresent ? (" value=" + to_string(firstOrdinal)) : "") << ">"
        " Last marker ordinal: <input type=text name=lastOrdinal required" <<
        (lastOrdinalIsPresent ? (" value=" + to_string(lastOrdinal)) : "") << ">";
    html << "<br><input type=radio name=whichAlignments value=AllAlignments" <<
        (whichAlignments=="AllAlignments" ? " checked=checked" : "") << "> All alignments";
    html << "<br><input type=radio name=whichAlignments value=ReadGraphAlignments" <<
        (whichAlignments=="ReadGraphAlignments" ? " checked=checked" : "") <<
        "> Only alignments used in the read graph.";
    html << "</form>";

    // If a required parameter is missing, stop here.
    if(!readId0IsPresent || !strand0IsPresent || !firstOrdinalIsPresent || !lastOrdinalIsPresent) {
        return;
    }
    const OrientedReadId orientedReadId0(readId0, strand0);
    const uint32_t markerCount0 = uint32_t(markers->size(orientedReadId0.getValue()));
    const uint32_t ordinal0Begin = firstOrdinal;
    const uint32_t ordinal0End = min(markerCount0, lastOrdinal + 1);




    // Loop over alignment involving this oriented read, as stored in the
    // alignment table.
    Alignment alignment;
    vector<OrientedReadId> orientedReadIds1;
    vector<bool> isInReadGraph;
    vector< vector<uint32_t> > alignedOrdinals1Matrix; // alignedOrdinals1Matrix[i][ordinal0] = ordinal1;
    const auto alignmentTable0 = alignmentTable[orientedReadId0.getValue()];
    for(const auto alignmentId: alignmentTable0) {
        const AlignmentData& ad = alignmentData[alignmentId];

        // If this alignment is not in the read graph and only read graph alignments
        // were requested, skip it.
        if((whichAlignments=="ReadGraphAlignments") and (not ad.info.isInReadGraph)) {
            continue;
        }

        // The alignment is stored with its first read on strand 0.
        OrientedReadId alignmentOrientedReadId0(ad.readIds[0], 0);
        OrientedReadId alignmentOrientedReadId1(ad.readIds[1],
            ad.isSameStrand ? 0 : 1);

        // Access the alignment and decompress it.
        const span<char> compressedAlignment = compressedAlignments[alignmentId];
        const span<const char> constCompressedAlignment(compressedAlignment.begin(), compressedAlignment.end());
        decompress(constCompressedAlignment, alignment);
        DINARA_ASSERT(alignment.ordinals.size() == ad.info.markerCount);

        // Swap the reads, if necessary.
        bool swapReads = false;
        if(alignmentOrientedReadId0.getReadId() != orientedReadId0.getReadId()) {
            swap(alignmentOrientedReadId0, alignmentOrientedReadId1);
            swapReads = true;
        }

        // Reverse complement, if necessary.
        bool reverseComplement = false;
        if(alignmentOrientedReadId0 != orientedReadId0) {
            alignmentOrientedReadId0.flipStrand();
            alignmentOrientedReadId1.flipStrand();
            reverseComplement = true;
        }
        DINARA_ASSERT(alignmentOrientedReadId0 == orientedReadId0);
        const OrientedReadId orientedReadId1 = alignmentOrientedReadId1;
        const uint32_t markerCount1 = uint32_t(markers->size(orientedReadId1.getValue()));
        orientedReadIds1.push_back(orientedReadId1);
        isInReadGraph.push_back(ad.info.isInReadGraph);

        // Store aligned ordinals of orientedReadId1.
        alignedOrdinals1Matrix.resize(orientedReadIds1.size());
        vector<uint32_t>& alignedOrdinals1 = alignedOrdinals1Matrix.back();
        alignedOrdinals1.resize(markerCount0, std::numeric_limits<uint32_t>::max());
        for(const auto& ordinals: alignment.ordinals) {
            uint32_t ordinal0 = ordinals[0];
            uint32_t ordinal1 = ordinals[1];
            if(swapReads) {
                swap(ordinal0, ordinal1);
            }
            if(reverseComplement) {
                ordinal0 = markerCount0 - 1 - ordinal0;
                ordinal1 = markerCount1 - 1 - ordinal1;
            }
            DINARA_ASSERT(alignedOrdinals1[ordinal0] == std::numeric_limits<uint32_t>::max());
            alignedOrdinals1[ordinal0] = ordinal1;
        }


    }



    // Write the page header.
    html << "<h1>Follow oriented read " << orientedReadId0 <<
        " and its alignments in the marker graph</h1>"
        "<p>This follows oriented read " << orientedReadId0 <<
        " and its alignments in the marker graph.";
    if(whichAlignments=="ReadGraphAlignments") {
        html << " You selected to only display alignments that "
            "correspond to a read graph edge.";
    } else {
        html << " Alignments that correspond to a read graph edge are displayed "
            "with a light blue background.";
    }


    // Write the table header.
    html << "<table><tr><th colspan=2 style='background-color:Beige'>" << orientedReadId0;
    for(uint64_t i=0; i<orientedReadIds1.size(); i++) {
        html << "<th";
        if(isInReadGraph[i]) {
            html << " style='background-color:LightCyan'";
        }
        html << " colspan=2>" << orientedReadIds1[i];
    }
    html << "<tr>"
        "<th style='background-color:Beige'>Marker<br>ordinal"
        "<th style='background-color:Beige'>Marker<br>graph<br>vertex";
    for(uint64_t i=0; i<orientedReadIds1.size(); i++) {
        html << "<th";
        if(isInReadGraph[i]) {
            html << " style='background-color:LightCyan'";
        }
        html << ">";
        html << "Marker<br>ordinal<th";
        if(isInReadGraph[i]) {
            html << " style='background-color:LightCyan'";
        }
        html << ">Marker<br>graph<br>vertex";
    }



    // Write a table row for each marker in orientedReadId0.
    for(uint32_t ordinal0=ordinal0Begin; ordinal0<ordinal0End; ordinal0++) {
        const MarkerId markerId0 = getMarkerId(orientedReadId0, ordinal0);
        const MarkerGraph::CompressedVertexId vertexId0 = markerGraph.vertexTable[markerId0];

        // Write ordinal and marker graph vertex for orientedReadId0.
        html << "<tr><td class=centered title='" << orientedReadId0 << " ordinal'"
            " style='background-color:Beige'>" << ordinal0 <<
            "<td class=centered title='" << orientedReadId0 << " marker graph vertex'";
        if(vertexId0 == MarkerGraph::invalidCompressedVertexId) {
            html << " style='background-color:Beige'>";
        }
        else {
            html << " style='background-color:Aquamarine'>" <<
                "<a href='exploreMarkerGraphVertex?vertexId=" << vertexId0 << "'>" << vertexId0 << "</a>";
        }


        // Loop over aligned reads.
        for(uint64_t i=0; i<orientedReadIds1.size(); i++) {
            const OrientedReadId orientedReadId1 = orientedReadIds1[i];
            const vector<uint32_t>& alignedOrdinals1 = alignedOrdinals1Matrix[i];

            // Marker ordinal.
            html << "<td class=centered title='" << orientedReadId1 << " ordinal'";
            if(isInReadGraph[i]) {
                html << " style='background-color:LightCyan'";
            }
            html << ">";
            const uint32_t ordinal1 = alignedOrdinals1[ordinal0];
            if(ordinal1 != std::numeric_limits<uint32_t>::max()) {
                html << ordinal1;
            }



            // Marker graph vertex.
            if(ordinal1 == std::numeric_limits<uint32_t>::max()) {

                // There is no aligned ordinal.
                html << "<td class=centered title='" << orientedReadId1 << " marker graph vertex'";
                if(isInReadGraph[i]) {
                    html << " style='background-color:LightCyan'";
                }
                html << ">";

            } else {

                // There is an aligned ordinal. Look for a marker graph vertex.
                const MarkerId markerId1 = getMarkerId(orientedReadId1, ordinal1);
                const MarkerGraph::CompressedVertexId vertexId1 = markerGraph.vertexTable[markerId1];
                html << "<td class=centered title='" << orientedReadId1 << " marker graph vertex'";
                if(vertexId1 == MarkerGraph::invalidCompressedVertexId) {

                    // There is no marker graph vertex.
                    if(isInReadGraph[i]) {
                        html << " style='background-color:LightCyan'";
                    }
                    html << ">";

                } else {

                    // There is a marker graph vertex.
                    if(vertexId0 != MarkerGraph::invalidCompressedVertexId) {
                        // Color based on whether or not it is the same as vertexId0.
                        if(vertexId1 == vertexId0) {
                            html << "style='background-color:Aquamarine'";
                        } else {
                            html << "style='background-color:LightPink'";
                        }
                    } else {
                        if(isInReadGraph[i]) {
                            html << " style='background-color:LightCyan'";
                        }
                    }
                    html << ">" <<
                        "<a href='exploreMarkerGraphVertex?vertexId=" << vertexId1 << "'>" << vertexId1 << "</a>";
                }
            }
        }
    }

    // Finish the table.
    html << "</table>";

}



void Assembler::exploreMarkerConnectivity(
    const vector<string>& request,
    ostream& html)
{
    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);
    uint32_t ordinal = 0;
    const bool ordinalIsPresent = getParameterValue(request, "ordinal", ordinal);
    string whichAlignments = "ReadGraphAlignments";
    getParameterValue(request, "whichAlignments", whichAlignments);
    string labelsString;
    const bool labels = getParameterValue(request, "labels", labelsString);
    double timeout = 30;
    getParameterValue(request, "timeout", timeout);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Explore connectivity of this marker'> "
        "<br>Read <input type=text name=readId required" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " size=8 title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);
    html <<
        "<br>Marker ordinal <input type=text name=ordinal required";
    if(ordinalIsPresent) {
        html << " value=" << ordinal;
    }
    html <<
        ">"
        "<br><input type=radio name=whichAlignments value=AllAlignments" <<
        (whichAlignments=="AllAlignments" ? " checked=checked" : "") << "> Use all alignments";
    html << "<br><input type=radio name=whichAlignments value=ReadGraphAlignments" <<
        (whichAlignments=="ReadGraphAlignments" ? " checked=checked" : "") <<
        "> Only use alignments in the read graph.";
    html << "<br><input type=checkbox name=labels" <<
         (labels ? " checked" : "") <<
         "> Labels"
        "<br>Timeout (seconds) for graph layout"
        " <input type=text required name=timeout size=8 style='text-align:center'" <<
        " value='" << timeout <<
        "'>"
        "</form>";
    const bool useReadGraphAlignmentsOnly = (whichAlignments == "ReadGraphAlignments");

    // If the required parameters are missing, stop here.
    if(not(readIdIsPresent and strandIsPresent and ordinalIsPresent)) {
        return;
    }
    const OrientedReadId orientedReadId(readId, strand);

    // Check the ordinal.
    const uint64_t markerCount = markers->size(orientedReadId.getValue());
    if(ordinal >= markerCount) {
        html << "<p>" << orientedReadId << " has " << markerCount << " markers.";
        return;
    }



    // Create an undirected graph in which each vertex represents a marker
    // and aligned markers are joined by an edge.
    MarkerConnectivityGraph graph;
    createMarkerConnectivityGraph(orientedReadId, ordinal, useReadGraphAlignmentsOnly, graph);



    // Count how many times each oriented read appears.
    std::map<OrientedReadId, uint64_t> frequencyMap;
    BGL_FORALL_VERTICES(v, graph, MarkerConnectivityGraph) {
        const OrientedReadId orientedReadId = graph[v].first;
        ++frequencyMap[orientedReadId];
    }

    html << "<br>The marker connectivity graph has " <<
        num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges.";


    // Write the graph out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    ofstream dotFile(dotFileName);
    dotFile << "graph MarkerConnectivity {\n";
    BGL_FORALL_VERTICES(v, graph, MarkerConnectivityGraph) {
        const MarkerDescriptor markerDescriptor = graph[v];
        const OrientedReadId orientedReadId1 = markerDescriptor.first;
        const uint32_t ordinal1 = markerDescriptor.second;
        dotFile << "\"" << orientedReadId1 << "-" << ordinal1 << "\"";
        if(labels) {
            dotFile <<
                " [label=\"" << orientedReadId1 << "\\n" << ordinal1 <<
                "\"";
            if(frequencyMap[orientedReadId1] != 1) {
                dotFile << " style=filled fillcolor=pink";
            } else {
                dotFile << " style=filled fillcolor=cornsilk";
            }
            dotFile << "]";
        }
        dotFile << ";\n";
    }
    BGL_FORALL_EDGES(e, graph, MarkerConnectivityGraph) {
        const auto v0 = source(e, graph);
        const auto v1 = target(e, graph);
        const MarkerDescriptor markerDescriptor0 = graph[v0];
        const MarkerDescriptor markerDescriptor1 = graph[v1];
        dotFile
            << "\"" << markerDescriptor0.first << "-" << markerDescriptor0.second << "\"--"
            << "\"" << markerDescriptor1.first << "-" << markerDescriptor1.second << "\";\n";
    }
    dotFile << "}\n";
    dotFile.close();



    // Use graphviz to render it to svg.
    const string command = timeoutCommand() + " " + to_string(int(timeout)) + " sfdp -O -T svg " + dotFileName +
        ( labels ? " -Goverlap=false -Gsplines=true -Gsmoothing=triangle" :
            " -Nshape=point -Gsize=10 -Gratio=expand -Epenwidth=0.4");
    const int commandStatus = ::system(command.c_str());
    if(WIFEXITED(commandStatus)) {
        const int exitStatus = WEXITSTATUS(commandStatus);
        if(exitStatus == 124) {
            html << "<p>Timeout for graph layout exceeded.";
            std::filesystem::remove(dotFileName);
            return;
        }
        else if(exitStatus!=0) {
            // filesystem::remove(dotFileName);
            throw runtime_error("Error " + to_string(exitStatus) + " running graph layout command: " + command);
        }
    } else if(WIFSIGNALED(commandStatus)) {
        const int signalNumber = WTERMSIG(commandStatus);
        throw runtime_error("Signal " + to_string(signalNumber) + " while running graph layout command: " + command);
    } else {
        throw runtime_error("Abnormal status " + to_string(commandStatus) + " while running graph layout command: " + command);
    }

    // Remove the .dot file.
    std::filesystem::remove(dotFileName);

    // Buttons to resize the svg locally.
    const int sizePixels = 800;
    addScaleSvgButtons(html, sizePixels);
    html << "<br>Found " << num_vertices(graph) << " markers.";

    // Display the svg file.
    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << "<div id=svgDiv style='display:none'>"; // Make it invisible until after we scale it.
    html << svgFile.rdbuf();
    svgFile.close();

    // Scale to desired size, then make it visible.
    html <<
        "</div>"
        "<script>"
        "var svgElement = document.getElementsByTagName('svg')[0];"
        "svgElement.setAttribute('width', " << sizePixels << ");"
        "document.getElementById('svgDiv').setAttribute('style', 'display:block');"
        "</script>";

    // Remove the .svg file.
    std::filesystem::remove(svgFileName);
}

