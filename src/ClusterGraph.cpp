// Dinara.
#include "ClusterGraph.hpp"
#include "Assembler.hpp"
#include "computeLayout.hpp"
#include "deduplicate.hpp"
#include "html.hpp"
#include "HttpServer.hpp"
#include "orderPairs.hpp"
#include "platformDependent.hpp"
#include "Reads.hpp"
#include "runCommandWithTimeout.hpp"
#include "timestamp.hpp"
using namespace dinara;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <map>
#include <set>



// Create the ClusterGraph from the Assembler's variant clustering data.
// Only includes valid clusters.
// The graph is directed: edges go from lower position cluster to higher position cluster.
ClusterGraph::ClusterGraph(
    const Assembler& assembler,
    uint64_t minEdgeCoverage)
{
    // Collect all valid cluster IDs.
    const uint64_t clusterCount = assembler.variantClusteringValidClusters.size();
    for(uint64_t clusterId = 0; clusterId < clusterCount; clusterId++) {
        if(assembler.variantClusteringValidClusters[clusterId]) {
            clusterIds.push_back(clusterId);
        }
    }

    // Sort cluster IDs.
    std::sort(clusterIds.begin(), clusterIds.end());

    // Check if allele data is available.
    const bool hasAlleleData = assembler.variantClusteringPositionPairAlleles.isOpen;

    // Create vertices.
    for(uint64_t localClusterId = 0; localClusterId < clusterIds.size(); localClusterId++) {
        const uint64_t clusterId = clusterIds[localClusterId];
        ClusterGraphVertex vertex(clusterId, localClusterId);
        
        // Get coverage and compute strand balance.
        uint64_t strand0Count = 0;
        uint64_t strand1Count = 0;
        std::array<uint32_t, 5> alleleCounts = {0};
        
        for(uint64_t memberIdx : assembler.variantClusteringMembersByRepIdx[clusterId]) {
            const auto& pp = assembler.variantClusteringPositionPairs[memberIdx];
            if(pp.first.getStrand() == 0) {
                strand0Count++;
            } else {
                strand1Count++;
            }
            
            // Count alleles (only if allele data is available).
            if(hasAlleleData && memberIdx < assembler.variantClusteringPositionPairAlleles.size()) {
                const uint8_t allele = assembler.variantClusteringPositionPairAlleles[memberIdx];
                if(allele < 5) {
                    alleleCounts[allele]++;
                }
            }
        }
        
        vertex.coverage = strand0Count + strand1Count;
        if(vertex.coverage > 0) {
            vertex.strand0Fraction = double(strand0Count) / double(vertex.coverage);
        }
        
        // Count valid alleles.
        if(hasAlleleData) {
            for(int a = 0; a < 5; a++) {
                if(alleleCounts[a] >= 5) { // minAlleleCoverage
                    vertex.alleleCount++;
                }
            }
        }
        
        const vertex_descriptor v = add_vertex(vertex, *this);
        vertexDescriptors.push_back(v);
        clusterIdToVertex[clusterId] = v;
    }

    // Build edges by following reads.
    buildEdgesFromReads(assembler, minEdgeCoverage);
    
    // Find and connect reverse-complement cluster pairs.
    findAndConnectRcPairs(assembler);
}



// Create the ClusterGraph for a specific set of cluster IDs.
// The graph is directed: edges go from lower position cluster to higher position cluster.
ClusterGraph::ClusterGraph(
    const Assembler& assembler,
    const std::vector<uint64_t>& inputClusterIds,
    uint64_t minEdgeCoverage)
{
    // Copy and sort cluster IDs.
    clusterIds = inputClusterIds;
    std::sort(clusterIds.begin(), clusterIds.end());
    
    // Remove duplicates.
    clusterIds.erase(std::unique(clusterIds.begin(), clusterIds.end()), clusterIds.end());

    // Check if allele data is available.
    const bool hasAlleleData = assembler.variantClusteringPositionPairAlleles.isOpen;

    // Create vertices.
    for(uint64_t localClusterId = 0; localClusterId < clusterIds.size(); localClusterId++) {
        const uint64_t clusterId = clusterIds[localClusterId];
        ClusterGraphVertex vertex(clusterId, localClusterId);
        
        // Get coverage and compute strand balance.
        if(clusterId < assembler.variantClusteringMembersByRepIdx.size()) {
            uint64_t strand0Count = 0;
            uint64_t strand1Count = 0;
            std::array<uint32_t, 5> alleleCounts = {0};
            
            for(uint64_t memberIdx : assembler.variantClusteringMembersByRepIdx[clusterId]) {
                const auto& pp = assembler.variantClusteringPositionPairs[memberIdx];
                if(pp.first.getStrand() == 0) {
                    strand0Count++;
                } else {
                    strand1Count++;
                }
                
                // Count alleles (only if allele data is available).
                if(hasAlleleData && memberIdx < assembler.variantClusteringPositionPairAlleles.size()) {
                    const uint8_t allele = assembler.variantClusteringPositionPairAlleles[memberIdx];
                    if(allele < 5) {
                        alleleCounts[allele]++;
                    }
                }
            }
            
            vertex.coverage = strand0Count + strand1Count;
            if(vertex.coverage > 0) {
                vertex.strand0Fraction = double(strand0Count) / double(vertex.coverage);
            }
            
            // Count valid alleles.
            if(hasAlleleData) {
                for(int a = 0; a < 5; a++) {
                    if(alleleCounts[a] >= 5) {
                        vertex.alleleCount++;
                    }
                }
            }
        }
        
        const vertex_descriptor v = add_vertex(vertex, *this);
        vertexDescriptors.push_back(v);
        clusterIdToVertex[clusterId] = v;
    }

    // Build edges by following reads.
    buildEdgesFromReads(assembler, minEdgeCoverage);
    
    // Find and connect reverse-complement cluster pairs.
    findAndConnectRcPairs(assembler);
}



// Save the ClusterGraph to a binary file.
void ClusterGraph::save(const std::string& fileName) const
{
    std::ofstream out(fileName, std::ios::binary);
    if(!out) {
        throw std::runtime_error("Cannot open file for writing: " + fileName);
    }
    
    const ClusterGraph& graph = *this;
    
    // Write number of vertices.
    uint64_t numVertices = num_vertices(graph);
    out.write(reinterpret_cast<const char*>(&numVertices), sizeof(numVertices));
    
    // Write vertex data.
    BGL_FORALL_VERTICES(v, graph, ClusterGraph) {
        const ClusterGraphVertex& vertex = graph[v];
        out.write(reinterpret_cast<const char*>(&vertex.clusterId), sizeof(vertex.clusterId));
        out.write(reinterpret_cast<const char*>(&vertex.localClusterId), sizeof(vertex.localClusterId));
        out.write(reinterpret_cast<const char*>(&vertex.coverage), sizeof(vertex.coverage));
        out.write(reinterpret_cast<const char*>(&vertex.alleleCount), sizeof(vertex.alleleCount));
        out.write(reinterpret_cast<const char*>(&vertex.strand0Fraction), sizeof(vertex.strand0Fraction));
        out.write(reinterpret_cast<const char*>(&vertex.rcPartnerClusterId), sizeof(vertex.rcPartnerClusterId));
    }
    
    // Write number of edges.
    uint64_t numEdges = num_edges(graph);
    out.write(reinterpret_cast<const char*>(&numEdges), sizeof(numEdges));
    
    // Write edge data (source clusterId, target clusterId, edge properties).
    BGL_FORALL_EDGES(e, graph, ClusterGraph) {
        const ClusterGraphEdge& edge = graph[e];
        uint64_t sourceClusterId = graph[source(e, graph)].clusterId;
        uint64_t targetClusterId = graph[target(e, graph)].clusterId;
        
        out.write(reinterpret_cast<const char*>(&sourceClusterId), sizeof(sourceClusterId));
        out.write(reinterpret_cast<const char*>(&targetClusterId), sizeof(targetClusterId));
        out.write(reinterpret_cast<const char*>(&edge.coverage), sizeof(edge.coverage));
        out.write(reinterpret_cast<const char*>(&edge.averageOffset), sizeof(edge.averageOffset));
        out.write(reinterpret_cast<const char*>(&edge.isRcPairEdge), sizeof(edge.isRcPairEdge));
    }
}



// Load the ClusterGraph from a binary file.
void ClusterGraph::load(const std::string& fileName)
{
    std::ifstream in(fileName, std::ios::binary);
    if(!in) {
        throw std::runtime_error("Cannot open file for reading: " + fileName);
    }
    
    // Clear existing data.
    clear();
    clusterIds.clear();
    vertexDescriptors.clear();
    clusterIdToVertex.clear();
    
    // Read number of vertices.
    uint64_t numVertices;
    in.read(reinterpret_cast<char*>(&numVertices), sizeof(numVertices));
    
    // Read vertex data and create vertices.
    for(uint64_t i = 0; i < numVertices; i++) {
        ClusterGraphVertex vertex;
        in.read(reinterpret_cast<char*>(&vertex.clusterId), sizeof(vertex.clusterId));
        in.read(reinterpret_cast<char*>(&vertex.localClusterId), sizeof(vertex.localClusterId));
        in.read(reinterpret_cast<char*>(&vertex.coverage), sizeof(vertex.coverage));
        in.read(reinterpret_cast<char*>(&vertex.alleleCount), sizeof(vertex.alleleCount));
        in.read(reinterpret_cast<char*>(&vertex.strand0Fraction), sizeof(vertex.strand0Fraction));
        in.read(reinterpret_cast<char*>(&vertex.rcPartnerClusterId), sizeof(vertex.rcPartnerClusterId));
        
        const vertex_descriptor v = add_vertex(vertex, *this);
        clusterIds.push_back(vertex.clusterId);
        vertexDescriptors.push_back(v);
        clusterIdToVertex[vertex.clusterId] = v;
    }
    
    // Read number of edges.
    uint64_t numEdges;
    in.read(reinterpret_cast<char*>(&numEdges), sizeof(numEdges));
    
    // Read edge data and create edges.
    for(uint64_t i = 0; i < numEdges; i++) {
        uint64_t sourceClusterId, targetClusterId;
        ClusterGraphEdge edge;
        
        in.read(reinterpret_cast<char*>(&sourceClusterId), sizeof(sourceClusterId));
        in.read(reinterpret_cast<char*>(&targetClusterId), sizeof(targetClusterId));
        in.read(reinterpret_cast<char*>(&edge.coverage), sizeof(edge.coverage));
        in.read(reinterpret_cast<char*>(&edge.averageOffset), sizeof(edge.averageOffset));
        in.read(reinterpret_cast<char*>(&edge.isRcPairEdge), sizeof(edge.isRcPairEdge));
        
        auto itSource = clusterIdToVertex.find(sourceClusterId);
        auto itTarget = clusterIdToVertex.find(targetClusterId);
        if(itSource != clusterIdToVertex.end() && itTarget != clusterIdToVertex.end()) {
            boost::add_edge(itSource->second, itTarget->second, edge, *this);
        }
    }
}



// Extract a subgraph containing only clusters within 'distance' hops from startClusterId.
ClusterGraph ClusterGraph::extractNeighborhood(
    uint64_t startClusterId, 
    uint64_t distance,
    uint64_t minEdgeCoverage) const
{
    const ClusterGraph& graph = *this;
    
    // Check if start cluster exists.
    auto itStart = clusterIdToVertex.find(startClusterId);
    if(itStart == clusterIdToVertex.end()) {
        return ClusterGraph();  // Return empty graph.
    }
    
    // BFS to find all clusters within 'distance'.
    std::set<uint64_t> reachableClusters;
    std::vector<uint64_t> currentLevel;
    std::vector<uint64_t> nextLevel;
    
    currentLevel.push_back(startClusterId);
    reachableClusters.insert(startClusterId);
    
    for(uint64_t d = 0; d < distance && !currentLevel.empty(); d++) {
        nextLevel.clear();
        for(uint64_t clusterId : currentLevel) {
            auto it = clusterIdToVertex.find(clusterId);
            if(it == clusterIdToVertex.end()) continue;
            
            vertex_descriptor v = it->second;
            
            // Check out-edges.
            BGL_FORALL_OUTEDGES(v, e, graph, ClusterGraph) {
                // Skip edges below coverage threshold (but always include RC pair edges).
                if(!graph[e].isRcPairEdge && graph[e].coverage < minEdgeCoverage) continue;
                
                uint64_t neighborId = graph[target(e, graph)].clusterId;
                if(reachableClusters.find(neighborId) == reachableClusters.end()) {
                    reachableClusters.insert(neighborId);
                    nextLevel.push_back(neighborId);
                }
            }
            
            // Check in-edges (for bidirectional traversal).
            BGL_FORALL_INEDGES(v, e, graph, ClusterGraph) {
                if(!graph[e].isRcPairEdge && graph[e].coverage < minEdgeCoverage) continue;
                
                uint64_t neighborId = graph[source(e, graph)].clusterId;
                if(reachableClusters.find(neighborId) == reachableClusters.end()) {
                    reachableClusters.insert(neighborId);
                    nextLevel.push_back(neighborId);
                }
            }
        }
        currentLevel = std::move(nextLevel);
    }
    
    // Convert to vector.
    std::vector<uint64_t> clusterIdsToInclude(reachableClusters.begin(), reachableClusters.end());
    
    return extractSubgraph(clusterIdsToInclude, minEdgeCoverage);
}



// Extract a subgraph containing only the specified cluster IDs.
ClusterGraph ClusterGraph::extractSubgraph(
    const std::vector<uint64_t>& clusterIdsToInclude,
    uint64_t minEdgeCoverage) const
{
    const ClusterGraph& graph = *this;
    ClusterGraph subgraph;
    
    // Build set for fast lookup.
    std::set<uint64_t> includeSet(clusterIdsToInclude.begin(), clusterIdsToInclude.end());
    
    // Create vertices.
    uint64_t localClusterId = 0;
    for(uint64_t clusterId : clusterIdsToInclude) {
        auto it = clusterIdToVertex.find(clusterId);
        if(it == clusterIdToVertex.end()) continue;
        
        // Copy vertex data.
        ClusterGraphVertex vertex = graph[it->second];
        vertex.localClusterId = localClusterId++;
        
        const vertex_descriptor v = add_vertex(vertex, subgraph);
        subgraph.clusterIds.push_back(clusterId);
        subgraph.vertexDescriptors.push_back(v);
        subgraph.clusterIdToVertex[clusterId] = v;
    }
    
    // Create edges (only between included vertices, filtered by coverage).
    BGL_FORALL_EDGES(e, graph, ClusterGraph) {
        const ClusterGraphEdge& edge = graph[e];
        
        // Apply coverage filter (but always include RC pair edges).
        if(!edge.isRcPairEdge && edge.coverage < minEdgeCoverage) continue;
        
        uint64_t sourceClusterId = graph[source(e, graph)].clusterId;
        uint64_t targetClusterId = graph[target(e, graph)].clusterId;
        
        // Check if both endpoints are in the subgraph.
        if(includeSet.find(sourceClusterId) == includeSet.end()) continue;
        if(includeSet.find(targetClusterId) == includeSet.end()) continue;
        
        auto itSource = subgraph.clusterIdToVertex.find(sourceClusterId);
        auto itTarget = subgraph.clusterIdToVertex.find(targetClusterId);
        if(itSource != subgraph.clusterIdToVertex.end() && 
           itTarget != subgraph.clusterIdToVertex.end()) {
            boost::add_edge(itSource->second, itTarget->second, edge, subgraph);
        }
    }
    
    return subgraph;
}



// Creates directed edges: source cluster has lower position, target cluster has higher position.
void ClusterGraph::buildEdgesFromReads(
    const Assembler& assembler,
    uint64_t minEdgeCoverage)
{
    // Build a map from (sourceClusterId, targetClusterId) -> (coverage, sumOffset).
    // Directed: source has lower position on read, target has higher position.
    std::map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, int64_t>> edgeMap;

    const auto& positionPairs = assembler.variantClusteringPositionPairs;
    
    // Build a lookup table from position pair index to cluster ID.
    // We do this by iterating through variantClusteringMembersByRepIdx.
    const uint64_t positionPairCount = positionPairs.size();
    std::vector<uint64_t> positionPairToClusterId(positionPairCount, std::numeric_limits<uint64_t>::max());
    
    for(uint64_t clusterId = 0; clusterId < assembler.variantClusteringMembersByRepIdx.size(); clusterId++) {
        const auto& members = assembler.variantClusteringMembersByRepIdx[clusterId];
        for(uint64_t memberIdx : members) {
            if(memberIdx < positionPairCount) {
                positionPairToClusterId[memberIdx] = clusterId;
            }
        }
    }

    // Iterate over all reads.
    const uint64_t readCount = assembler.getReads().readCount();
    
    for(ReadId readId = 0; readId < readCount; readId++) {
        // Process both strands.
        for(Strand strand = 0; strand < 2; strand++) {
            OrientedReadId orientedReadId(readId, strand);
            
            // Find position pairs for this read.
            auto it = std::lower_bound(
                positionPairs.begin(), 
                positionPairs.end(), 
                std::make_pair(orientedReadId, uint32_t(0)));

            // Collect clusters for this read, sorted by position.
            struct ClusterOnRead {
                uint64_t clusterId;
                uint32_t position;
            };
            std::vector<ClusterOnRead> clustersOnRead;
            
            while(it != positionPairs.end() && it->first == orientedReadId) {
                uint64_t index = it - positionPairs.begin();
                uint64_t clusterId = positionPairToClusterId[index];
                
                // Only include clusters that are in this graph.
                if(clusterId != std::numeric_limits<uint64_t>::max() &&
                   clusterIdToVertex.find(clusterId) != clusterIdToVertex.end()) {
                    clustersOnRead.push_back({clusterId, it->second});
                }
                ++it;
            }

            // Sort by position.
            std::sort(clustersOnRead.begin(), clustersOnRead.end(),
                [](const ClusterOnRead& a, const ClusterOnRead& b) {
                    return a.position < b.position;
                });

            // Remove duplicate cluster appearances (keep first occurrence).
            std::vector<ClusterOnRead> uniqueClusters;
            std::set<uint64_t> seenClusters;
            for(const auto& c : clustersOnRead) {
                if(seenClusters.find(c.clusterId) == seenClusters.end()) {
                    uniqueClusters.push_back(c);
                    seenClusters.insert(c.clusterId);
                }
            }

            // Create directed edges between consecutive clusters.
            // Direction: lower position -> higher position (already sorted by position).
            for(size_t i = 1; i < uniqueClusters.size(); i++) {
                // Source cluster is at lower position, target cluster is at higher position.
                const uint64_t sourceClusterId = uniqueClusters[i-1].clusterId;
                const uint64_t targetClusterId = uniqueClusters[i].clusterId;
                const int64_t offset = int64_t(uniqueClusters[i].position) - int64_t(uniqueClusters[i-1].position);
                
                // Store directed edge (source, target) based on position order.
                auto key = std::make_pair(sourceClusterId, targetClusterId);
                
                auto& edgeData = edgeMap[key];
                edgeData.first++;  // coverage
                edgeData.second += offset;  // sum of offsets (always positive since sorted by position)
            }
        }
    }

    // Create edges from the map.
    for(const auto& [key, data] : edgeMap) {
        const uint64_t coverage = data.first;
        if(coverage < minEdgeCoverage) {
            continue;
        }
        
        const int64_t averageOffset = data.second / int64_t(coverage);
        
        // key.first = source cluster, key.second = target cluster
        auto itSource = clusterIdToVertex.find(key.first);
        auto itTarget = clusterIdToVertex.find(key.second);
        
        if(itSource != clusterIdToVertex.end() && itTarget != clusterIdToVertex.end()) {
            addEdge(itSource->second, itTarget->second, coverage, averageOffset);
        }
    }
}



void ClusterGraph::addEdge(
    vertex_descriptor v0,
    vertex_descriptor v1,
    uint64_t coverage,
    int64_t averageOffset)
{
    ClusterGraphEdge edge(coverage, averageOffset);
    boost::add_edge(v0, v1, edge, *this);
}



// Find and connect reverse-complement cluster pairs.
// For each cluster, find a member's RC position pair and look up its cluster.
void ClusterGraph::findAndConnectRcPairs(const Assembler& assembler)
{
    const auto& positionPairs = assembler.variantClusteringPositionPairs;
    const auto& membersByRepIdx = assembler.variantClusteringMembersByRepIdx;
    auto& disjointSets = *assembler.variantClusteringDisjointSets;  // Non-const: find() modifies internal state
    
    std::set<std::pair<uint64_t, uint64_t>> addedRcEdges;
    
    for(const auto& [clusterId, v] : clusterIdToVertex) {
        if(clusterId >= membersByRepIdx.size()) continue;
        if(membersByRepIdx[clusterId].empty()) continue;
        
        // Take the first member of this cluster.
        uint64_t memberIdx = membersByRepIdx[clusterId][0];
        const auto& pp = positionPairs[memberIdx];
        
        ReadId readId = pp.first.getReadId();
        Strand strand = pp.first.getStrand();
        uint32_t position = pp.second;
        
        // Find the RC: same read, opposite strand, complementary position.
        Strand rcStrand = 1 - strand;
        OrientedReadId rcOrientedReadId(readId, rcStrand);
        uint64_t readLength = assembler.getReads().getReadRawSequenceLength(readId);
        uint32_t rcPosition = uint32_t(readLength - 1 - position);
        
        // Binary search for the RC position pair.
        // positionPairs is sorted by (OrientedReadId, position).
        auto searchKey = std::make_pair(rcOrientedReadId, rcPosition);
        auto it = std::lower_bound(positionPairs.begin(), positionPairs.end(), searchKey);
        
        if(it == positionPairs.end() || *it != searchKey) continue;
        
        uint64_t rcMemberIdx = it - positionPairs.begin();
        
        // Find the cluster (representative) for the RC member.
        uint64_t rcClusterId = disjointSets.find(rcMemberIdx);
        
        if(rcClusterId == clusterId) continue;  // Same cluster (shouldn't happen)
        
        // Check if the RC cluster is in this graph.
        auto itRc = clusterIdToVertex.find(rcClusterId);
        if(itRc == clusterIdToVertex.end()) continue;
        
        // Create edge key (smaller ID first to avoid duplicates).
        auto key = std::make_pair(std::min(clusterId, rcClusterId), std::max(clusterId, rcClusterId));
        
        if(addedRcEdges.find(key) == addedRcEdges.end()) {
            auto itV0 = clusterIdToVertex.find(key.first);
            auto itV1 = clusterIdToVertex.find(key.second);
            
            if(itV0 != clusterIdToVertex.end() && itV1 != clusterIdToVertex.end()) {
                // Add RC pair edge.
                ClusterGraphEdge rcEdge = ClusterGraphEdge::createRcPairEdge();
                boost::add_edge(itV0->second, itV1->second, rcEdge, *this);
                
                // Update vertex RC partner info.
                (*this)[itV0->second].rcPartnerClusterId = key.second;
                (*this)[itV1->second].rcPartnerClusterId = key.first;
                
                addedRcEdges.insert(key);
            }
        }
    }
}



// Count RC pair edges.
uint64_t ClusterGraph::rcPairEdgeCount() const
{
    uint64_t count = 0;
    const ClusterGraph& graph = *this;
    BGL_FORALL_EDGES(e, graph, ClusterGraph) {
        if(graph[e].isRcPairEdge) {
            count++;
        }
    }
    return count;
}



// Write edge details to a CSV file.
void ClusterGraph::writeEdgeDetails(const std::string& fileName) const
{
    std::ofstream out(fileName);
    out << "ClusterId0,ClusterId1,Coverage,AverageOffset\n";
    
    const ClusterGraph& graph = *this;
    BGL_FORALL_EDGES(e, graph, ClusterGraph) {
        const ClusterGraphEdge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        
        out << graph[v0].clusterId << ","
            << graph[v1].clusterId << ","
            << edge.coverage << ","
            << edge.averageOffset << "\n";
    }
}



// Write vertex details to a CSV file.
void ClusterGraph::writeVertexDetails(const std::string& fileName) const
{
    std::ofstream out(fileName);
    out << "ClusterId,LocalClusterId,Coverage,AlleleCount,InDegree,OutDegree\n";
    
    const ClusterGraph& graph = *this;
    BGL_FORALL_VERTICES(v, graph, ClusterGraph) {
        const ClusterGraphVertex& vertex = graph[v];
        
        out << vertex.clusterId << ","
            << vertex.localClusterId << ","
            << vertex.coverage << ","
            << vertex.alleleCount << ","
            << in_degree(v, graph) << ","
            << out_degree(v, graph) << "\n";
    }
}



// Remove edges with coverage below threshold.
// RC pair edges are never removed (they have isRcPairEdge = true).
void ClusterGraph::removeWeakEdges(uint64_t minCoverage)
{
    ClusterGraph& graph = *this;
    
    // Collect edges to remove.
    std::vector<edge_descriptor> edgesToRemove;
    BGL_FORALL_EDGES(e, graph, ClusterGraph) {
        // Never remove RC pair edges.
        if(graph[e].isRcPairEdge) {
            continue;
        }
        if(graph[e].coverage < minCoverage) {
            edgesToRemove.push_back(e);
        }
    }
    
    // Remove them.
    for(const auto& e : edgesToRemove) {
        boost::remove_edge(e, graph);
    }
}



// ClusterGraphDisplayOptions constructor from HTTP request.
ClusterGraphDisplayOptions::ClusterGraphDisplayOptions(const std::vector<std::string>& request)
{
    HttpServer::getParameterValue(request, "sizePixels", sizePixels);
    HttpServer::getParameterValue(request, "edgeLengthScale", edgeLengthScale);
    HttpServer::getParameterValue(request, "timeout", timeout);
    HttpServer::getParameterValue(request, "vertexSize", vertexSize);
    
    std::string vertexSizeByCoverageString;
    vertexSizeByCoverage = HttpServer::getParameterValue(request, 
        "vertexSizeByCoverage", vertexSizeByCoverageString);
    
    std::string vertexLabelsString;
    vertexLabels = HttpServer::getParameterValue(request, 
        "vertexLabels", vertexLabelsString);
    
    HttpServer::getParameterValue(request, "vertexColoring", vertexColoring);
    HttpServer::getParameterValue(request, "edgeThickness", edgeThickness);
    
    std::string edgeLabelsString;
    edgeLabels = HttpServer::getParameterValue(request, 
        "edgeLabels", edgeLabelsString);
    
    HttpServer::getParameterValue(request, "edgeColoring", edgeColoring);
    HttpServer::getParameterValue(request, "arrowSize", arrowSize);
    HttpServer::getParameterValue(request, "lowCoverage", lowCoverage);
    HttpServer::getParameterValue(request, "highCoverage", highCoverage);
}



void ClusterGraphDisplayOptions::writeForm(std::ostream& html) const
{
    html <<
        "<tr>"
        "<th title='Graphics size in pixels.'>"
        "Graphics size in pixels"
        "<td class=centered><input type=text required name=sizePixels size=8 style='text-align:center'" <<
        " value='" << sizePixels << "'>";

    html <<
        "<tr>"
        "<th>Layout"
        "<td class=left>"
        "<input type=text name=edgeLengthScale style='text-align:center' required size=6 value=" <<
        edgeLengthScale << "> Edge length scale (larger = more spread out)"
        "<br><input type=text name=timeout style='text-align:center' required size=6 value=" <<
        timeout << "> Timeout (seconds)";

    html <<
        "<tr>"
        "<th>Vertices"
        "<td class=left>"
        "<input type=text name=vertexSize style='text-align:center' required size=6 value=" <<
        vertexSize << "> Vertex size (arbitrary units)"
        "<br><input type=checkbox name=vertexSizeByCoverage" <<
        (vertexSizeByCoverage ? " checked" : "") <<
        "> Size proportional to coverage"
        "<hr>"
        "<input type=checkbox name=vertexLabels" <<
        (vertexLabels ? " checked" : "") << "> Labels (dot layout only)"
        "<hr>"
        "<b>Vertex coloring</b>"
        "<br><input type=radio required name=vertexColoring value='black'" <<
        (vertexColoring == "black" ? " checked=on" : "") << ">Black"
        "<br><input type=radio required name=vertexColoring value='byAlleleCount'" <<
        (vertexColoring == "byAlleleCount" ? " checked=on" : "") << ">By allele count";

    html <<
        "<tr>"
        "<th>Edges"
        "<td class=left>"
        "<b>Edge coloring</b>"
        "<br><input type=radio required name=edgeColoring value='black'" <<
        (edgeColoring == "black" ? " checked=on" : "") << ">Black"
        "<br><input type=radio required name=edgeColoring value='byCoverage'" <<
        (edgeColoring == "byCoverage" ? " checked=on" : "") << ">By coverage"
        "<hr>"
        "<b>Edge graphics</b>"
        "<br><input type=text name=edgeThickness style='text-align:center' required size=6 value=" <<
        edgeThickness << "> Thickness (arbitrary units)"
        "<br><input type=text name=arrowSize style='text-align:center' required size=6 value=" <<
        arrowSize << "> Arrow size (arbitrary units)"
        "<hr>"
        "<input type=checkbox name=edgeLabels" <<
        (edgeLabels ? " checked" : "") << "> Labels (dot layout only)"
        "<hr>"
        "<b>Coverage coloring thresholds</b>"
        "<br><input type=text name=lowCoverage style='text-align:center' required size=6 value=" <<
        lowCoverage << "> Low coverage (red)"
        "<br><input type=text name=highCoverage style='text-align:center' required size=6 value=" <<
        highCoverage << "> High coverage (green)";
}



// Write the graph in Graphviz format to a stream.
void ClusterGraph::writeGraphviz(
    std::ostream& out,
    const ClusterGraphDisplayOptions& options) const
{
    const ClusterGraph& graph = *this;
    out << "digraph ClusterGraph {\n";
    out << "rankdir=LR;\n";
    out << "node [shape=point width=" << (0.05 * options.vertexSize) << "];\n";
    out << "edge [arrowsize=" << options.arrowSize << "];\n";

    BGL_FORALL_VERTICES(v, graph, ClusterGraph) {
        const ClusterGraphVertex& vertex = graph[v];
        out << vertex.clusterId;
        out << " [";
        
        // Tooltip.
        out << "tooltip=\"Cluster " << vertex.clusterId;
        out << ", coverage " << vertex.coverage;
        out << ", " << vertex.alleleCount << " alleles\"";
        
        // Size.
        if(options.vertexSizeByCoverage) {
            double size = 0.05 * options.vertexSize * std::sqrt(double(vertex.coverage));
            out << " width=" << size;
        }
        
        // Labels (for dot layout).
        if(options.vertexLabels) {
            out << " shape=rectangle label=\"C" << vertex.clusterId << "\\ncov=" << vertex.coverage << "\"";
        }
        
        // Vertex coloring.
        if(options.vertexColoring == "byAlleleCount") {
            if(vertex.alleleCount >= 2) {
                out << " color=\"green\" fillcolor=\"lightgreen\" style=filled";
            } else if(vertex.alleleCount == 1) {
                out << " color=\"orange\" fillcolor=\"lightyellow\" style=filled";
            } else {
                out << " color=\"gray\" fillcolor=\"lightgray\" style=filled";
            }
        }
        
        // URL to explore this cluster.
        out << " URL=\"exploreVariantCluster?clusterIdx=" << vertex.clusterId << "\"";
        
        out << "];\n";
    }

    BGL_FORALL_EDGES(e, graph, ClusterGraph) {
        const ClusterGraphEdge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        
        out << graph[v0].clusterId << " -> " << graph[v1].clusterId;
        out << " [";
        
        // Tooltip.
        out << "tooltip=\"Coverage " << edge.coverage << ", offset " << edge.averageOffset << "\"";
        
        // Label.
        if(options.edgeLabels) {
            out << " label=\"" << edge.coverage << "\"";
        }
        
        // Thickness.
        double penwidth = options.edgeThickness * std::sqrt(double(edge.coverage));
        out << " penwidth=" << penwidth;
        
        // Edge coloring.
        if(options.edgeColoring == "byCoverage") {
            double fraction = 0.0;
            if(options.highCoverage > options.lowCoverage) {
                fraction = double(edge.coverage - options.lowCoverage) / 
                           double(options.highCoverage - options.lowCoverage);
                fraction = std::max(0.0, std::min(1.0, fraction));
            }
            // Red to green gradient.
            uint8_t r = static_cast<uint8_t>(255 * (1.0 - fraction));
            uint8_t g = static_cast<uint8_t>(255 * fraction);
            char color[10];
            snprintf(color, sizeof(color), "#%02x%02x00", r, g);
            out << " color=\"" << color << "\"";
        }
        
        out << "];\n";
    }

    out << "}\n";
}



// Write the graph in Graphviz format to a file.
void ClusterGraph::writeGraphviz(
    const std::string& fileName,
    const ClusterGraphDisplayOptions& options) const
{
    std::ofstream out(fileName);
    writeGraphviz(out, options);
}



// Compute the graph layout using Graphviz.
void ClusterGraph::computeLayout(const ClusterGraphDisplayOptions& options)
{
    layout.clear();
    
    // Build additional Graphviz options.
    std::string additionalOptions;
    
    // K parameter: ideal edge length (controls spacing in sfdp, fdp, neato)
    if(options.edgeLengthScale != 1.0) {
        double kValue = options.edgeLengthScale;
        additionalOptions = "-Gk=" + std::to_string(kValue);
    }
    
    computeLayoutGraphviz(
        *this,
        "sfdp",
        options.timeout,
        layout,
        additionalOptions);
}



void ClusterGraph::computeLayoutBoundingBox()
{
    boundingBox.xMin = std::numeric_limits<double>::max();
    boundingBox.xMax = std::numeric_limits<double>::lowest();
    boundingBox.yMin = boundingBox.xMin;
    boundingBox.yMax = boundingBox.xMax;
    
    for(const auto& p: layout) {
        const std::array<double, 2>& xy = p.second;
        const double x = xy[0];
        const double y = xy[1];
        boundingBox.xMin = std::min(boundingBox.xMin, x);
        boundingBox.xMax = std::max(boundingBox.xMax, x);
        boundingBox.yMin = std::min(boundingBox.yMin, y);
        boundingBox.yMax = std::max(boundingBox.yMax, y);
    }
}



void ClusterGraph::Box::makeSquare()
{
    if(xSize() > ySize()) {
        const double delta = (xSize() - ySize()) / 2.;
        yMin -= delta;
        yMax += delta;
    } else {
        const double delta = (ySize() - xSize()) / 2.;
        xMin -= delta;
        xMax += delta;
    }
}



void ClusterGraph::Box::extend(double factor)
{
    const double ext = factor * std::max(xSize(), ySize());
    xMin -= ext;
    xMax += ext;
    yMin -= ext;
    yMax += ext;
}



void ClusterGraph::writeVertices(
    std::ostream& html, 
    const ClusterGraphDisplayOptions& options) const
{
    const ClusterGraph& graph = *this;
    
    const double scalingFactor = 0.002;  // sfdp scaling factor
    
    html << "\n<g id='vertices' style='stroke:none'>";
    
    BGL_FORALL_VERTICES(v, graph, ClusterGraph) {
        const ClusterGraphVertex& vertex = graph[v];
        
        // Get position.
        const auto it = layout.find(v);
        if(it == layout.end()) continue;
        const double x = it->second[0];
        const double y = it->second[1];
        
        // Compute radius.
        double radius = scalingFactor * options.vertexSize;
        if(options.vertexSizeByCoverage) {
            radius *= std::sqrt(double(vertex.coverage));
        }
        
        // Choose color.
        std::string color = "black";
        if(options.vertexColoring == "byAlleleCount") {
            if(vertex.alleleCount >= 2) {
                color = "green";
            } else if(vertex.alleleCount == 1) {
                color = "orange";
            } else {
                color = "gray";
            }
        }
        
        html << "\n<circle id='C" << vertex.clusterId << "'"
             << " cx='" << x << "' cy='" << y << "' r='" << radius << "'"
             << " fill='" << color << "'"
             << " onclick=\"window.location='exploreVariantCluster?clusterIdx=" << vertex.clusterId << "'\""
             << " style='cursor:pointer'>"
             << "<title>Cluster " << vertex.clusterId 
             << ", coverage " << vertex.coverage
             << ", " << vertex.alleleCount << " alleles</title>"
             << "</circle>";
    }
    
    html << "\n</g>";
}



void ClusterGraph::writeEdges(
    std::ostream& html,
    const ClusterGraphDisplayOptions& options) const
{
    const ClusterGraph& graph = *this;
    
    const double scalingFactor = 0.002;  // sfdp scaling factor
    
    // Draw edges.
    html << "\n<g id='edges'>";
    
    BGL_FORALL_EDGES(e, graph, ClusterGraph) {
        const ClusterGraphEdge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        
        // Get positions.
        const auto it0 = layout.find(v0);
        const auto it1 = layout.find(v1);
        if(it0 == layout.end() || it1 == layout.end()) continue;
        
        const double x0 = it0->second[0];
        const double y0 = it0->second[1];
        const double x1 = it1->second[0];
        const double y1 = it1->second[1];
        
        // Handle RC pair edges differently.
        if(edge.isRcPairEdge) {
            // RC pair edges: dashed magenta lines.
            double thickness = scalingFactor * options.edgeThickness * 2.0;  // Fixed thickness
            html << "\n<line x1='" << x0 << "' y1='" << y0 
                 << "' x2='" << x1 << "' y2='" << y1 
                 << "' stroke='magenta' stroke-width='" << thickness 
                 << "' stroke-dasharray='" << (thickness * 3) << "," << (thickness * 2) << "'>"
                 << "<title>RC Pair Edge</title>"
                 << "</line>";
            continue;
        }
        
        // Compute thickness for normal edges.
        double thickness = scalingFactor * options.edgeThickness * std::sqrt(double(edge.coverage));
        
        // Choose color.
        std::string color = "black";
        if(options.edgeColoring == "byCoverage") {
            double fraction = 0.0;
            if(options.highCoverage > options.lowCoverage) {
                fraction = double(edge.coverage - options.lowCoverage) / 
                           double(options.highCoverage - options.lowCoverage);
                fraction = std::max(0.0, std::min(1.0, fraction));
            }
            uint8_t r = static_cast<uint8_t>(255 * (1.0 - fraction));
            uint8_t g = static_cast<uint8_t>(255 * fraction);
            char colorBuf[10];
            snprintf(colorBuf, sizeof(colorBuf), "#%02x%02x00", r, g);
            color = colorBuf;
        }
        
        html << "\n<line x1='" << x0 << "' y1='" << y0 
             << "' x2='" << x1 << "' y2='" << y1 
             << "' stroke='" << color << "' stroke-width='" << thickness << "'>"
             << "<title>Coverage " << edge.coverage 
             << ", offset " << edge.averageOffset << "</title>"
             << "</line>";
    }
    
    html << "\n</g>";
    
    // Draw arrows (skip RC pair edges).
    html << "\n<g id='arrows'>";
    
    BGL_FORALL_EDGES(e, graph, ClusterGraph) {
        const ClusterGraphEdge& edge = graph[e];
        
        // No arrows for RC pair edges.
        if(edge.isRcPairEdge) continue;
        
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        
        const auto it0 = layout.find(v0);
        const auto it1 = layout.find(v1);
        if(it0 == layout.end() || it1 == layout.end()) continue;
        
        const double x0 = it0->second[0];
        const double y0 = it0->second[1];
        const double x1 = it1->second[0];
        const double y1 = it1->second[1];
        
        std::string color = "black";
        if(options.edgeColoring == "byCoverage") {
            double fraction = 0.0;
            if(options.highCoverage > options.lowCoverage) {
                fraction = double(edge.coverage - options.lowCoverage) / 
                           double(options.highCoverage - options.lowCoverage);
                fraction = std::max(0.0, std::min(1.0, fraction));
            }
            uint8_t r = static_cast<uint8_t>(255 * (1.0 - fraction));
            uint8_t g = static_cast<uint8_t>(255 * fraction);
            char colorBuf[10];
            snprintf(colorBuf, sizeof(colorBuf), "#%02x%02x00", r, g);
            color = colorBuf;
        }
        
        // Arrow at 70% of the way to target.
        const double relativeArrowLength = 0.3;
        const double x2 = (1. - relativeArrowLength) * x1 + relativeArrowLength * x0;
        const double y2 = (1. - relativeArrowLength) * y1 + relativeArrowLength * y0;
        
        html << "\n<line x1='" << x1 << "' y1='" << y1 
             << "' x2='" << x2 << "' y2='" << y2 
             << "' stroke='" << color << "' stroke-width='" << (0.2 * scalingFactor * options.arrowSize * double(edge.coverage)) << "' />";
    }
    
    html << "\n</g>";
}



void ClusterGraph::writeSvgControls(std::ostream& html) const
{
    html << "<p><table>";

    // Add drag and zoom.
    addSvgDragAndZoom(html);

    // Buttons to change vertex size.
    html << R"stringDelimiter(
    <tr><th class=left>Vertex size<td>
    <button type='button' onClick='changeVertexSize(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeVertexSize(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeVertexSize(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeVertexSize(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeVertexSize(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeVertexSize(10.)' style='width:3em'>+++</button>
        <script>
        function changeVertexSize(factor)
        {
            var vertexGroup = document.getElementById('vertices');
            var vertices = vertexGroup.getElementsByTagName('circle');
            for(i=0; i<vertices.length; i++) {
                v = vertices[i];
                v.setAttribute('r', factor * v.getAttribute('r'));
            }
        }
        </script>
        )stringDelimiter";

    // Buttons to change edge thickness.
    html << R"stringDelimiter(
    <tr><th class=left>Edge thickness<td>
    <button type='button' onClick='changeThickness(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeThickness(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeThickness(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeThickness(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeThickness(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeThickness(10.)' style='width:3em'>+++</button>
        <script>
        function changeThickness(factor)
        {
            var edgeGroup = document.getElementById('edges');
            var edges = edgeGroup.getElementsByTagName('line');
            for(i=0; i<edges.length; i++) {
                e = edges[i];
                e.setAttribute('stroke-width', factor * e.getAttribute('stroke-width'));
            }

            var arrowsGroup = document.getElementById('arrows');
            var arrows = arrowsGroup.getElementsByTagName('line');
            for(i=0; i<arrows.length; i++) {
                a = arrows[i];
                a.setAttribute('stroke-width', factor * a.getAttribute('stroke-width'));
            }
        }
        </script>
        )stringDelimiter";

    // Zoom buttons.
    html << R"stringDelimiter(
        <tr title='Or use the mouse wheel.'><th class=left>Zoom<td>
        <button type='button' onClick='zoomSvg(0.1)' style='width:3em'>---</button>
        <button type='button' onClick='zoomSvg(0.5)' style='width:3em'>--</button>
        <button type='button' onClick='zoomSvg(0.8)' style='width:3em'>-</button>
        <button type='button' onClick='zoomSvg(1.25)' style='width:3em'>+</button>
        <button type='button' onClick='zoomSvg(2.)' style='width:3em'>++</button>
        <button type='button' onClick='zoomSvg(10.)' style='width:3em'>+++</button>
    )stringDelimiter";

    // Buttons to highlight a cluster and zoom to a cluster.
    html << R"stringDelimiter(
        <tr><td colspan=2>
        <button onClick='highlightCluster()'>Highlight</button>
        <button onClick='zoomToCluster()'>Zoom to</button>cluster
        <input id=selectedClusterId type=text size=10 style='text-align:center'>
    <script>
    function zoomToCluster()
    {
        var clusterId = document.getElementById("selectedClusterId").value;
        zoomToGivenCluster(clusterId);
    }
    function zoomToGivenCluster(clusterId)
    {
        var element = document.getElementById("C" + clusterId);
        if(!element) {
            alert("Cluster " + clusterId + " not found");
            return;
        }
        var box = element.getBBox();
        var xCenter = box.x + 0.5 * box.width;
        var yCenter = box.y + 0.5 * box.height;

        var enlargeFactor = 5.;
        var size = enlargeFactor * Math.max(box.width, box.height);
        var factor = size / width;
        width = size;
        height = size;
        x = xCenter - 0.5 * size;
        y = yCenter - 0.5 * size;
        var svg = document.querySelector('svg');
        svg.setAttribute('viewBox', `${x} ${y} ${size} ${size}`);
        ratio = size / svg.getBoundingClientRect().width;
        svg.setAttribute('font-size', svg.getAttribute('font-size') * factor);
    }
    function highlightCluster()
    {
        var clusterId = document.getElementById("selectedClusterId").value;
        var element = document.getElementById("C" + clusterId);
        if(!element) {
            alert("Cluster " + clusterId + " not found");
            return;
        }
        element.style.fill = "Magenta";
    }
    </script>
    )stringDelimiter";

    html << "</table>";

    // Scroll down to the svg.
    const std::string svgId = "ClusterGraph";
    html <<
        "<script>"
        "document.getElementById('" << svgId << "').scrollIntoView({block:'center'});"
        "</script>";

    html <<
        "<p>Use Ctrl+Click to pan."
        "<p>Use Ctrl-Wheel or the above buttons to zoom.";
}



// Write the graph as interactive HTML/SVG.
void ClusterGraph::writeHtml(
    std::ostream& html,
    const ClusterGraphDisplayOptions& options)
{
    // Use scientific notation because svg does not accept floating points
    // ending with a decimal point.
    html << std::scientific;

    // Compute layout.
    computeLayout(options);
    computeLayoutBoundingBox();

    Box viewportBox = boundingBox;
    viewportBox.extend(0.05);
    viewportBox.makeSquare();

    // Begin the svg.
    const std::string svgId = "ClusterGraph";
    html <<
        "\n<br><div style='display:inline-block;vertical-align:top;'>"
        "<svg id='" << svgId <<
        "' width='" << options.sizePixels <<
        "' height='" << options.sizePixels <<
        "' viewbox='" << viewportBox.xMin << " " << viewportBox.yMin << " " <<
        viewportBox.xSize() << " " <<
        viewportBox.ySize() << "'"
        " style='background-color:#f0f0f0'"
        ">\n";

    // Write the edges first so they don't obscure the vertices.
    writeEdges(html, options);

    // Write the vertices.
    writeVertices(html, options);

    // Finish the svg.
    html << "</svg></div>";

    // Side panel with controls.
    html << "<div style='display:inline-block;margin-left:20px'>";
    writeSvgControls(html);
    html << "</div>";
}

