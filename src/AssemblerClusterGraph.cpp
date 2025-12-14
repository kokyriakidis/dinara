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
#include <iomanip>
#include <limits>
#include <sstream>



// Create the ClusterGraph from variant clustering data and save it.
void Assembler::createClusterGraph(uint64_t minEdgeCoverage)
{
    cout << timestamp << "Creating ClusterGraph with minEdgeCoverage=" << minEdgeCoverage << endl;
    
    clusterGraph = std::make_shared<ClusterGraph>(*this, minEdgeCoverage);

    // Update edge coverages to count all shared reads (not just consecutive).
    cout << timestamp << "Updating edge coverages (shared read count)..." << endl;
    clusterGraph->updateEdgeCoverages(*this);

    // Prune edges that violate reverse complement symmetry.
    cout << timestamp << "Pruning asymmetric edges..." << endl;
    clusterGraph->pruneAsymmetricEdges();

    // Simplify graph with transitive reduction.
    cout << timestamp << "Applying transitive reduction..." << endl;
    // clusterGraph->transitiveReduction();

    // Resolve forks using local read following.
    cout << timestamp << "Resolving forks with read following..." << endl;
    // clusterGraph->resolveForks();
    
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



// HTTP server function to explore a single variant cluster (shows reads and alleles).
void Assembler::exploreVariantCluster(const vector<string>& request, ostream& html)
{
    // Get the cluster ID.
    uint64_t clusterIdx = std::numeric_limits<uint64_t>::max();
    HttpServer::getParameterValue(request, "clusterIdx", clusterIdx);
    
    html << "<h1>Variant Cluster Explorer</h1>";
    
    // Write form for entering cluster ID.
    html << "<form><table>";
    html << "<tr><th class=left>Cluster ID"
         << "<td><input type=text name=clusterIdx value=\"";
    if(clusterIdx != std::numeric_limits<uint64_t>::max()) {
        html << clusterIdx;
    }
    html << "\" size=12 style='text-align:center'>";
    html << "</table><input type=submit value='Explore Cluster'></form>";
    
    // If no cluster specified, return.
    if(clusterIdx == std::numeric_limits<uint64_t>::max()) {
        return;
    }
    
    // Check if variant clustering data is available.
    if(!variantClusteringMembersByRepIdx.isOpen()) {
        html << "<p>Variant clustering data not available. "
             << "Run variant clustering first or access the data.";
        return;
    }
    
    // Check if cluster exists.
    if(clusterIdx >= variantClusteringMembersByRepIdx.size()) {
        html << "<p>Cluster " << clusterIdx << " not found (max cluster ID: " 
             << variantClusteringMembersByRepIdx.size() - 1 << ").";
        return;
    }
    
    const auto& members = variantClusteringMembersByRepIdx[clusterIdx];
    if(members.empty()) {
        html << "<p>Cluster " << clusterIdx << " has no members (may not be a valid cluster).";
        return;
    }
    
    const bool hasAlleleData = variantClusteringPositionPairAlleles.isOpen;
    
    // Show cluster summary.
    html << "<h2>Cluster " << clusterIdx << "</h2>";
    html << "<p>Coverage: " << members.size() << " position pairs.";
    
    // Check if valid cluster.
    if(variantClusteringValidClusters.isOpen && clusterIdx < variantClusteringValidClusters.size()) {
        if(variantClusteringValidClusters[clusterIdx]) {
            html << "<br>Status: <span style='color:green'>Valid cluster</span>";
        } else {
            html << "<br>Status: <span style='color:red'>Invalid cluster</span>";
        }
    }
    
    // Count alleles.
    std::array<uint32_t, 5> alleleCounts = {0};
    std::map<OrientedReadId, std::vector<std::pair<uint32_t, uint8_t>>> readData;  // position, allele
    
    for(uint64_t memberIdx : members) {
        const auto& pp = variantClusteringPositionPairs[memberIdx];
        uint8_t allele = 255;  // Unknown.
        if(hasAlleleData && memberIdx < variantClusteringPositionPairAlleles.size()) {
            allele = variantClusteringPositionPairAlleles[memberIdx];
            if(allele < 5) {
                alleleCounts[allele]++;
            }
        }
        readData[pp.first].push_back({pp.second, allele});
    }
    
    // Show allele counts.
    if(hasAlleleData) {
        html << "<h3>Allele Counts</h3>";
        html << "<table><tr><th>Allele<th>Count<th>Fraction</tr>";
        const char* baseNames[] = {"A", "C", "G", "T", "Gap"};
        uint64_t totalWithAllele = 0;
        for(int a = 0; a < 5; a++) {
            totalWithAllele += alleleCounts[a];
        }
        for(int a = 0; a < 5; a++) {
            if(alleleCounts[a] > 0) {
                double fraction = totalWithAllele > 0 ? 
                    double(alleleCounts[a]) / double(totalWithAllele) : 0.0;
                html << "<tr><td>" << baseNames[a] 
                     << "<td class=centered>" << alleleCounts[a]
                     << "<td class=centered>" << std::fixed << std::setprecision(2) << fraction;
            }
        }
        html << "</table>";
    }
    
    // Show reads.
    html << "<h3>Reads (" << readData.size() << " unique reads)</h3>";
    html << "<table><tr><th>Read<th>Strand<th>Positions";
    if(hasAlleleData) {
        html << "<th>Alleles";
    }
    html << "</tr>";
    
    for(const auto& rd : readData) {
        const OrientedReadId orientedReadId = rd.first;
        const auto& positions = rd.second;
        
        html << "<tr>";
        html << "<td><a href='exploreRead?readId=" << orientedReadId.getReadId() 
             << "&strand=" << orientedReadId.getStrand() << "'>" 
             << orientedReadId << "</a>";
        html << "<td class=centered>" << orientedReadId.getStrand();
        
        // Show positions.
        html << "<td class=centered>";
        bool first = true;
        for(const auto& pos : positions) {
            if(!first) html << ", ";
            first = false;
            html << pos.first;
        }
        
        // Show alleles.
        if(hasAlleleData) {
            html << "<td class=centered>";
            first = true;
            for(const auto& pos : positions) {
                if(!first) html << ", ";
                first = false;
                if(pos.second < 5) {
                    const char* baseNames[] = {"A", "C", "G", "T", "-"};
                    html << baseNames[pos.second];
                } else {
                    html << "?";
                }
            }
        }
        html << "</tr>";
    }
    html << "</table>";
    
    // Link to explore related clusters.
    html << "<h3>Related Clusters</h3>";
    html << "<p><a href='exploreClusterGraph?startClusterId=" << clusterIdx 
         << "&distance=1'>Show neighborhood (distance 1)</a>";
    html << "<br><a href='exploreClusterGraph?startClusterId=" << clusterIdx 
         << "&distance=2'>Show neighborhood (distance 2)</a>";
}



// HTTP server function to explore multiple variant clusters together.
void Assembler::exploreVariantClusters(const vector<string>& request, ostream& html)
{
    // Get the cluster IDs.
    string clusterIdsString;
    HttpServer::getParameterValue(request, "clusterIds", clusterIdsString);
    boost::trim(clusterIdsString);
    
    html << "<h1>Multi-Cluster Variant Explorer</h1>";
    
    // Write form.
    html << "<form><table>";
    html << "<tr><th class=left>Cluster IDs (comma separated)"
         << "<td><input type=text name=clusterIds value=\"" << clusterIdsString 
         << "\" size=60>";
    html << "</table><input type=submit value='Analyze Clusters'></form>";
    
    if(clusterIdsString.empty()) {
        html << "<p>Enter cluster IDs to analyze.";
        html << "<p><b>Tip:</b> You can select clusters from the "
             << "<a href='exploreClusterGraph'>Cluster Graph</a> and paste them here.";
        return;
    }
    
    // Parse cluster IDs.
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
    
    // Check if variant clustering data is available.
    if(!variantClusteringMembersByRepIdx.isOpen()) {
        html << "<p>Variant clustering data not available.";
        return;
    }
    
    const bool hasAlleleData = variantClusteringPositionPairAlleles.isOpen;
    
    // Collect all reads across all clusters.
    // Map: OrientedReadId -> Map<ClusterId, vector<pair<position, allele>>>
    std::map<OrientedReadId, std::map<uint64_t, std::vector<std::pair<uint32_t, uint8_t>>>> allReadData;
    std::map<uint64_t, std::array<uint32_t, 5>> clusterAlleleCounts;
    
    for(uint64_t clusterIdx : clusterIds) {
        if(clusterIdx >= variantClusteringMembersByRepIdx.size()) {
            html << "<p>Cluster " << clusterIdx << " not found.";
            continue;
        }
        
        const auto& members = variantClusteringMembersByRepIdx[clusterIdx];
        clusterAlleleCounts[clusterIdx] = {0};
        
        for(uint64_t memberIdx : members) {
            const auto& pp = variantClusteringPositionPairs[memberIdx];
            uint8_t allele = 255;
            if(hasAlleleData && memberIdx < variantClusteringPositionPairAlleles.size()) {
                allele = variantClusteringPositionPairAlleles[memberIdx];
                if(allele < 5) {
                    clusterAlleleCounts[clusterIdx][allele]++;
                }
            }
            allReadData[pp.first][clusterIdx].push_back({pp.second, allele});
        }
    }
    
    html << "<h2>Summary</h2>";
    html << "<p>Analyzing " << clusterIds.size() << " clusters with " 
         << allReadData.size() << " unique reads.";
    
    // Show per-cluster summaries.
    html << "<h3>Cluster Summaries</h3>";
    html << "<table><tr><th>Cluster<th>Coverage";
    if(hasAlleleData) {
        html << "<th>Alleles";
    }
    html << "<th>Actions</tr>";
    
    for(uint64_t clusterIdx : clusterIds) {
        if(clusterIdx >= variantClusteringMembersByRepIdx.size()) continue;
        
        const auto& members = variantClusteringMembersByRepIdx[clusterIdx];
        html << "<tr><td class=centered>" << clusterIdx;
        html << "<td class=centered>" << members.size();
        
        if(hasAlleleData) {
            html << "<td>";
            const auto& counts = clusterAlleleCounts[clusterIdx];
            const char* baseNames[] = {"A", "C", "G", "T", "-"};
            bool first = true;
            for(int a = 0; a < 5; a++) {
                if(counts[a] > 0) {
                    if(!first) html << " ";
                    first = false;
                    html << baseNames[a] << ":" << counts[a];
                }
            }
        }
        
        html << "<td><a href='exploreVariantCluster?clusterIdx=" << clusterIdx << "'>Details</a>";
        html << "</tr>";
    }
    html << "</table>";
    
    // Show read overlap matrix.
    html << "<h3>Read Overlap Matrix</h3>";
    html << "<p>Number of reads shared between clusters:";
    
    // Compute overlaps.
    std::map<std::pair<uint64_t, uint64_t>, uint64_t> overlaps;
    for(const auto& readEntry : allReadData) {
        const auto& clusterMap = readEntry.second;
        std::vector<uint64_t> readClusters;
        for(const auto& cm : clusterMap) {
            readClusters.push_back(cm.first);
        }
        for(size_t i = 0; i < readClusters.size(); i++) {
            for(size_t j = i; j < readClusters.size(); j++) {
                uint64_t c1 = readClusters[i];
                uint64_t c2 = readClusters[j];
                if(c1 > c2) std::swap(c1, c2);
                overlaps[{c1, c2}]++;
            }
        }
    }
    
    html << "<table><tr><th>";
    for(uint64_t c : clusterIds) {
        html << "<th>" << c;
    }
    html << "</tr>";
    
    for(uint64_t c1 : clusterIds) {
        html << "<tr><th>" << c1;
        for(uint64_t c2 : clusterIds) {
            uint64_t k1 = std::min(c1, c2);
            uint64_t k2 = std::max(c1, c2);
            auto it = overlaps.find({k1, k2});
            if(it != overlaps.end()) {
                html << "<td class=centered>" << it->second;
            } else {
                html << "<td class=centered>0";
            }
        }
        html << "</tr>";
    }
    html << "</table>";
    
    // Show all reads with their cluster memberships.
    html << "<h3>Read Details</h3>";
    html << "<p>Reads and their alleles in each cluster:";
    html << "<table><tr><th>Read<th>Strand";
    for(uint64_t c : clusterIds) {
        html << "<th>C" << c;
    }
    html << "</tr>";
    
    for(const auto& readEntry : allReadData) {
        const OrientedReadId orientedReadId = readEntry.first;
        const auto& clusterMap = readEntry.second;
        
        html << "<tr>";
        html << "<td><a href='exploreRead?readId=" << orientedReadId.getReadId() 
             << "&strand=" << orientedReadId.getStrand() << "'>" 
             << orientedReadId << "</a>";
        html << "<td class=centered>" << orientedReadId.getStrand();
        
        for(uint64_t clusterIdx : clusterIds) {
            html << "<td class=centered>";
            auto it = clusterMap.find(clusterIdx);
            if(it != clusterMap.end()) {
                const auto& positions = it->second;
                if(hasAlleleData) {
                    bool first = true;
                    for(const auto& pos : positions) {
                        if(!first) html << ",";
                        first = false;
                        if(pos.second < 5) {
                            const char* baseNames[] = {"A", "C", "G", "T", "-"};
                            html << baseNames[pos.second];
                        } else {
                            html << "?";
                        }
                    }
                } else {
                    html << positions.size();
                }
            } else {
                html << "-";
            }
        }
        html << "</tr>";
    }
    html << "</table>";
    
    // Link to cluster graph.
    html << "<h3>Visualize</h3>";
    html << "<p><a href='exploreClusterGraph?clusterIds=" << clusterIdsString 
         << "'>View as Cluster Graph</a>";
}