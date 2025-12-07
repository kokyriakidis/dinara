// Dinara.
#include "Assembler.hpp"
#include "ClusterGraph.hpp"
#include "timestamp.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/algorithm/string.hpp>

// Standard library.
#include "fstream.hpp"
#include <filesystem>
#include <limits>
#include <sstream>



// Create the ClusterGraph from variant clustering data and save it.
void Assembler::createClusterGraph(uint64_t minEdgeCoverage)
{
    cout << timestamp << "Creating ClusterGraph with minEdgeCoverage=" << minEdgeCoverage << endl;
    
    clusterGraph = std::make_shared<ClusterGraph>(*this, minEdgeCoverage);
    
    cout << timestamp << "ClusterGraph created with " 
         << clusterGraph->vertexCount() << " vertices and "
         << clusterGraph->edgeCount() << " edges "
         << "(" << clusterGraph->rcPairEdgeCount() << " RC pair edges)." << endl;
    
    // Save to binary file.
    const string fileName = largeDataFileNamePrefix.empty() ? 
        "ClusterGraph.bin" : (largeDataFileNamePrefix + "ClusterGraph.bin");
    clusterGraph->save(fileName);
    cout << timestamp << "ClusterGraph saved to " << fileName << endl;
}



// Load the ClusterGraph from binary file.
void Assembler::loadClusterGraph()
{
    const string fileName = largeDataFileNamePrefix.empty() ? 
        "ClusterGraph.bin" : (largeDataFileNamePrefix + "ClusterGraph.bin");
    
    if(!std::filesystem::exists(fileName)) {
        throw runtime_error("ClusterGraph file not found: " + fileName);
    }
    
    clusterGraph = std::make_shared<ClusterGraph>();
    clusterGraph->load(fileName);
    
    cout << timestamp << "ClusterGraph loaded with " 
         << clusterGraph->vertexCount() << " vertices and "
         << clusterGraph->edgeCount() << " edges." << endl;
}



// Access the ClusterGraph (load if not exists).
ClusterGraph& Assembler::getClusterGraph()
{
    if(!clusterGraph) {
        loadClusterGraph();
    }
    return *clusterGraph;
}



// HTTP server function to explore the cluster graph.
void Assembler::exploreClusterGraph(const vector<string>& request, ostream& html)
{
    // Get graph creation parameters.
    uint64_t minEdgeCoverage = 0;  // Default to 0 to show all edges
    HttpServer::getParameterValue(request, "minEdgeCoverage", minEdgeCoverage);
    
    string clusterIdsString;
    HttpServer::getParameterValue(request, "clusterIds", clusterIdsString);
    boost::trim(clusterIdsString);
    
    uint64_t startClusterId = std::numeric_limits<uint64_t>::max();
    HttpServer::getParameterValue(request, "startClusterId", startClusterId);
    
    uint64_t distance = 2;
    HttpServer::getParameterValue(request, "distance", distance);
    
    uint64_t sizeThreshold = 1000;
    HttpServer::getParameterValue(request, "sizeThreshold", sizeThreshold);
    
    // Get display options.
    const ClusterGraphDisplayOptions displayOptions(request);

    // Write the form.
    html << "<h1>Cluster Graph Explorer</h1>";
    html << "<form><table>";
    
    // Graph creation options.
    html << "<tr><th class=left>Minimum edge coverage"
         << "<td><input type=text name=minEdgeCoverage value=" << minEdgeCoverage 
         << " size=8 style='text-align:center'>";
    
    html << "<tr><th class=left>Cluster IDs (comma separated, leave empty to use start cluster + distance)"
         << "<td><input type=text name=clusterIds value=\"" << clusterIdsString 
         << "\" size=40>";
    
    html << "<tr><th class=left>Start cluster ID (used if cluster IDs is empty)"
         << "<td><input type=text name=startClusterId value=\"";
    if(startClusterId != std::numeric_limits<uint64_t>::max()) {
        html << startClusterId;
    }
    html << "\" size=12 style='text-align:center'>";
    
    html << "<tr><th class=left>Distance (number of hops from start cluster)"
         << "<td><input type=text name=distance value=" << distance 
         << " size=8 style='text-align:center'>";
    
    html << "<tr><th class=left>Max vertices (0 for unlimited)"
         << "<td><input type=text name=sizeThreshold value=" << sizeThreshold
         << " size=8 style='text-align:center'>";
    
    // Display options.
    displayOptions.writeForm(html);
    
    html << "</table><input type=submit value='Create Cluster Graph'></form>";

    // Try to load the pre-built ClusterGraph.
    const string clusterGraphFile = largeDataFileNamePrefix.empty() ? 
        "ClusterGraph.bin" : (largeDataFileNamePrefix + "ClusterGraph.bin");
    
    if(!std::filesystem::exists(clusterGraphFile)) {
        html << "<p>ClusterGraph file not found: " << clusterGraphFile << ".<br>"
             << "You need to create it first during assembly or by running createClusterGraph.";
        return;
    }
    
    // Load the full graph if not already loaded.
    if(!clusterGraph) {
        try {
            loadClusterGraph();
        } catch(const std::exception& e) {
            html << "<p>Error loading ClusterGraph: " << e.what();
            return;
        }
    }
    
    html << "<p>Full ClusterGraph has " << clusterGraph->vertexCount() 
         << " vertices and " << clusterGraph->edgeCount() << " edges.";

    // Extract the subgraph based on user parameters.
    std::unique_ptr<ClusterGraph> graph;
    
    if(!clusterIdsString.empty()) {
        // Parse cluster IDs from comma-separated string.
        std::vector<uint64_t> clusterIds;
        std::istringstream iss(clusterIdsString);
        std::string token;
        while(std::getline(iss, token, ',')) {
            boost::trim(token);
            if(!token.empty()) {
                try {
                    clusterIds.push_back(std::stoull(token));
                } catch(...) {
                    html << "<p>Invalid cluster ID: " << token;
                    return;
                }
            }
        }
        
        if(clusterIds.empty()) {
            html << "<p>No valid cluster IDs provided.";
            return;
        }
        
        graph = std::make_unique<ClusterGraph>(
            clusterGraph->extractSubgraph(clusterIds, minEdgeCoverage));
        
    } else if(startClusterId != std::numeric_limits<uint64_t>::max()) {
        // Use start cluster + distance.
        graph = std::make_unique<ClusterGraph>(
            clusterGraph->extractNeighborhood(startClusterId, distance, minEdgeCoverage));
        
        if(graph->vertexCount() == 0) {
            html << "<p>Cluster " << startClusterId << " not found in the graph.";
            return;
        }
        
        html << "<p>Showing neighborhood of cluster " << startClusterId 
             << " within distance " << distance << ".";
    } else {
        // Show the full graph (filtered by minEdgeCoverage).
        if(minEdgeCoverage > 0) {
            // Need to filter edges - create a copy with all clusters.
            graph = std::make_unique<ClusterGraph>(
                clusterGraph->extractSubgraph(clusterGraph->clusterIds, minEdgeCoverage));
        } else {
            // No filtering needed, but we need a copy to display.
            graph = std::make_unique<ClusterGraph>(*clusterGraph);
        }
    }

    // Check size threshold.
    if(sizeThreshold > 0 && graph->vertexCount() > sizeThreshold) {
        html << "<p>Graph has " << graph->vertexCount() << " vertices, "
             << "which exceeds the size threshold of " << sizeThreshold << ". "
             << "Increase the threshold or specify specific cluster IDs.";
        return;
    }

    html << "<h2>Cluster Graph</h2>";
    html << "<p>The cluster graph has " << graph->vertexCount() 
         << " vertices and " << graph->edgeCount() << " edges "
         << "(" << graph->rcPairEdgeCount() << " RC pair edges shown as dashed magenta lines).";

    // Write the interactive graph.
    graph->writeHtml(html, displayOptions);
}

