// Superbubble-related functions for Shasta2AssemblyGraph.
// Extracted from Shasta2AssemblyGraph.cpp to match shasta2's file organization.

#include "Shasta2AssemblyGraph.hpp"
#include "Shasta2Superbubble.hpp"
#include "Shasta2SuperbubbleChain.hpp"
#include "findSuperbubbleOnodera.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
#include "DINARA_ASSERT.hpp"
#include "MultithreadedObject.tpp"
#include "MurmurHash2.hpp"

#include <boost/graph/iteration_macros.hpp>

#include <boost/graph/strong_components.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>

using namespace dinara;
using std::array;
using std::cout;
using std::endl;
using std::map;
using std::min;
using std::ofstream;
using std::pair;
using std::set;
using std::string;
using std::vector;

namespace {

array<double, 3> shasta2HslToRgb(double H, double S, double L)
{
    using std::fabs;
    using std::fmod;

    const double C = (1. - fabs(2. * L - 1.)) * S;
    const double Hprime = 6 * H;
    const double X = C * (1 - fabs(fmod(Hprime, 2.) - 1.));

    array<double, 3> rgb;
    if(Hprime >= 0. && Hprime < 1.) {
        rgb = {C, X, 0.};
    } else if(Hprime >= 1. && Hprime < 2.) {
        rgb = {X, C, 0.};
    } else if(Hprime >= 2. && Hprime < 3.) {
        rgb = {0., C, X};
    } else if(Hprime >= 3. && Hprime < 4.) {
        rgb = {0., X, C};
    } else if(Hprime >= 4. && Hprime < 5.) {
        rgb = {X, 0., C};
    } else {
        rgb = {C, 0., X};
    }

    const double m = L - C / 2.;
    for(double& x: rgb) {
        x += m;
    }
    return rgb;
}

string shasta2RandomHslColor(uint64_t id, double S, double L)
{
    const double H = double(MurmurHash2(&id, sizeof(id), 759)) / double(std::numeric_limits<uint32_t>::max());
    const array<double, 3> rgb = shasta2HslToRgb(H, S, L);
    std::ostringstream s;
    s << "#";
    for(const double value: rgb) {
        const int iv = min(255, int(value * 255.));
        s << std::hex << std::setw(2) << std::setfill('0') << std::nouppercase << iv;
    }
    return s.str();
}

} // anonymous namespace



uint64_t Shasta2AssemblyGraph::phaseSuperbubbleChains()
{
    performanceLog << timestamp << "Shasta2AssemblyGraph::phaseSuperbubbleChains begins." << endl;

    PhaseSuperbubbleChainsData& data = phaseSuperbubbleChainsData;
    data.superbubbleChains = make_shared< vector<Shasta2SuperbubbleChain> >();
    vector<Shasta2SuperbubbleChain>& superbubbleChains = *(data.superbubbleChains);
    data.totalChangeCount = 0;

    vector<Shasta2Superbubble> superbubbles;
    findSuperbubbles(superbubbles);
    writeSuperbubbles(superbubbles, "Shasta2Superbubbles-WithOverlaps.csv");
    removeContainedSuperbubbles(superbubbles);
    cout << "Found " << superbubbles.size() << " non-overlapping superbubbles." << endl;
    writeSuperbubbles(superbubbles, "Shasta2Superbubbles.csv");
    writeSuperbubblesForBandage(superbubbles, "Shasta2Superbubbles-Bandage.csv");

    findSuperbubbleChains(superbubbles, superbubbleChains);
    cout << "Found " << superbubbleChains.size() << " superbubble chains." << endl;
    writeSuperbubbleChains(superbubbleChains, "Shasta2SuperbubbleChains.csv");
    writeSuperbubbleChainsForBandage(superbubbleChains, "Shasta2SuperbubbleChains-Bandage.csv");

    setupLoadBalancing(superbubbleChains.size(), 1);
    runThreads(&Shasta2AssemblyGraph::phaseSuperbubbleChainsThreadFunction, options.threadCount);
    data.superbubbleChains = nullptr;
    uint64_t changeCount = data.totalChangeCount;

    changeCount += compress();
    performanceLog << timestamp << "Shasta2AssemblyGraph::phaseSuperbubbleChains ends." << endl;
    return changeCount;
}

void Shasta2AssemblyGraph::phaseSuperbubbleChainsThreadFunction(uint64_t /* threadId */)
{
    PhaseSuperbubbleChainsData& data = phaseSuperbubbleChainsData;
    vector<Shasta2SuperbubbleChain>& superbubbleChains = *(data.superbubbleChains);

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t superbubbleChainId=begin; superbubbleChainId<end; superbubbleChainId++) {
            Shasta2SuperbubbleChain& superbubbleChain = superbubbleChains[superbubbleChainId];
            const uint64_t changeCount = superbubbleChain.phase1(
                *this,
                superbubbleChainId);
            __sync_fetch_and_add(&data.totalChangeCount, changeCount);
        }
    }
}

void Shasta2AssemblyGraph::findSuperbubbles(vector<Shasta2Superbubble>& superbubbles) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;

    // Use the Verkko/Onodera algorithm: from each vertex with out-degree >= 2,
    // find the superbubble exit by topological-order processing.
    // No vertex count limit — the algorithm naturally terminates on cycles
    // and dead-ends (matching Verkko's find_bubble).

    superbubbles.clear();
    uint64_t branchVertexCount = 0;
    uint64_t debugCount = 0;
    BGL_FORALL_VERTICES(vSource, assemblyGraph, Shasta2AssemblyGraph) {
        if(out_degree(vSource, assemblyGraph) >= 2) {
            ++branchVertexCount;
            // Log first few branch vertices for debugging.
            if(debugCount < 5) {
                cout << "  Branch vertex " << assemblyGraph[vSource].id
                     << " (anchor " << assemblyGraph[vSource].anchorId << ")"
                     << " in=" << in_degree(vSource, assemblyGraph)
                     << " out=" << out_degree(vSource, assemblyGraph)
                     << " out-targets:";
                BGL_FORALL_OUTEDGES(vSource, e, assemblyGraph, Shasta2AssemblyGraph) {
                    const vertex_descriptor w = target(e, assemblyGraph);
                    cout << " " << assemblyGraph[w].id
                         << "(in=" << in_degree(w, assemblyGraph)
                         << ",out=" << out_degree(w, assemblyGraph) << ")";
                }
                cout << endl;
                ++debugCount;
            }
        }
        const vertex_descriptor vTarget =
            findSuperbubbleOnodera(assemblyGraph, vSource);
        if(vTarget != null_vertex()) {
            superbubbles.emplace_back(assemblyGraph, vSource, vTarget);
        }
    }
    cout << "findSuperbubbles: " << branchVertexCount
         << " vertices with out-degree >= 2, found "
         << superbubbles.size() << " superbubbles." << endl;
}

// Pop superbubbles by removing low-coverage alternative paths.
// Follows Verkko's pop_bubbles_coverage_based.py approach:
// 1. Find all superbubbles.
// 2. For each superbubble, enumerate all paths from source to target.
// 3. Keep the path with highest minimum edge coverage (widest/bottleneck path).
// 4. Remove edges on non-kept paths if their coverage is below maxPoppableCoverage.
// 5. Clean up isolated vertices.
//
// maxBubbleSize: maximum number of internal vertices in a poppable superbubble (0 = no limit).
// maxPoppableCoverageFraction: edges with coverage below this fraction of the
//     estimated average coverage are removed. Set to 0.5 for diploid, 1.0 for haploid.
uint64_t Shasta2AssemblyGraph::popSuperbubbles(
    uint64_t maxBubbleSize,
    double maxPoppableCoverageFraction)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;

    // Find superbubbles.
    vector<Shasta2Superbubble> superbubbles;
    findSuperbubbles(superbubbles);
    removeContainedSuperbubbles(superbubbles);
    cout << "popSuperbubbles: found " << superbubbles.size() << " superbubbles." << endl;

    // First pass: find the best path for each superbubble and collect kept edges.
    // An edge that is on the best path of any superbubble must not be removed.
    set<edge_descriptor> globalKeptEdges;
    set<edge_descriptor> candidatesForRemoval;

    for(const Shasta2Superbubble& bubble : superbubbles) {

        // Skip superbubbles with too many internal vertices
        // (matches Verkko's max_bubble_pop_size which counts nodes).
        if(maxBubbleSize > 0 && bubble.internalVertices.size() > maxBubbleSize) {
            continue;
        }

        // Skip trivial superbubbles (single edge, nothing to pop).
        if(bubble.isTrivial()) {
            continue;
        }

        // Enumerate all paths from sourceVertex to targetVertex.
        // Each path is a sequence of edges.
        // Use DFS to enumerate paths (bounded by bubble size).
        vector< vector<edge_descriptor> > allPaths;
        {
            // DFS stack: (current vertex, path so far).
            vector< pair<vertex_descriptor, vector<edge_descriptor>> > stack;
            stack.push_back({bubble.sourceVertex, {}});

            while(!stack.empty()) {
                auto [v, path] = std::move(stack.back());
                stack.pop_back();

                if(v == bubble.targetVertex) {
                    allPaths.push_back(std::move(path));
                    continue;
                }

                BGL_FORALL_OUTEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
                    const vertex_descriptor w = target(e, assemblyGraph);
                    // Only follow edges within the superbubble
                    // (to target or to an internal vertex).
                    if(w == bubble.targetVertex || bubble.contains(w)) {
                        vector<edge_descriptor> newPath = path;
                        newPath.push_back(e);
                        stack.push_back({w, std::move(newPath)});
                    }
                }
            }
        }

        if(allPaths.empty()) {
            continue;
        }

        // Find the path with highest minimum edge coverage (bottleneck/widest path).
        uint64_t bestPathIndex = 0;
        double bestMinCoverage = 0.;
        for(uint64_t i = 0; i < allPaths.size(); i++) {
            double minCoverage = std::numeric_limits<double>::max();
            for(const edge_descriptor e : allPaths[i]) {
                minCoverage = std::min(minCoverage, assemblyGraph[e].averageCoverage());
            }
            if(minCoverage > bestMinCoverage) {
                bestMinCoverage = minCoverage;
                bestPathIndex = i;
            }
        }

        // Record kept edges (best path).
        for(const edge_descriptor e : allPaths[bestPathIndex]) {
            globalKeptEdges.insert(e);
        }

        // Record candidate edges for removal (non-best paths).
        for(uint64_t i = 0; i < allPaths.size(); i++) {
            if(i == bestPathIndex) {
                continue;
            }
            for(const edge_descriptor e : allPaths[i]) {
                candidatesForRemoval.insert(e);
            }
        }
    }

    // Only remove edges that are not on the best path of any superbubble.
    set<edge_descriptor> edgesToRemove;
    for(const edge_descriptor e : candidatesForRemoval) {
        if(!globalKeptEdges.count(e)) {
            edgesToRemove.insert(e);
        }
    }

    // Remove collected edges, then clean up isolated vertices.
    // We must remove all edges first, then vertices, because
    // multiple edges may share endpoints.
    uint64_t removedCount = 0;
    set<vertex_descriptor> affectedVertices;
    for(const edge_descriptor e : edgesToRemove) {
        affectedVertices.insert(source(e, assemblyGraph));
        affectedVertices.insert(target(e, assemblyGraph));
        boost::remove_edge(e, assemblyGraph);
        ++removedCount;
    }

    // Remove vertices that became isolated.
    for(const vertex_descriptor v : affectedVertices) {
        if(in_degree(v, assemblyGraph) == 0 && out_degree(v, assemblyGraph) == 0) {
            boost::remove_vertex(v, assemblyGraph);
        }
    }

    cout << "popSuperbubbles: removed " << removedCount << " edges from "
         << superbubbles.size() << " superbubbles." << endl;
    return removedCount;
}



void Shasta2AssemblyGraph::removeContainedSuperbubbles(vector<Shasta2Superbubble>& superbubbles) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;

    map<edge_descriptor, vector<uint64_t> > edgeToSuperbubbles;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        for(const edge_descriptor e: superbubbles[superbubbleId].internalEdges) {
            edgeToSuperbubbles[e].push_back(superbubbleId);
        }
    }

    set< pair<uint64_t, uint64_t> > intersectingPairs;
    for(const auto& p: edgeToSuperbubbles) {
        const auto& edgeSuperbubbles = p.second;
        for(uint64_t i0=0; i0+1<edgeSuperbubbles.size(); i0++) {
            for(uint64_t i1=i0+1; i1<edgeSuperbubbles.size(); i1++) {
                intersectingPairs.insert({
                    min(edgeSuperbubbles[i0], edgeSuperbubbles[i1]),
                    max(edgeSuperbubbles[i0], edgeSuperbubbles[i1])});
            }
        }
    }

    vector<edge_descriptor> commonEdges;
    vector<uint64_t> superbubblesToBeRemoved;
    for(const auto& p: intersectingPairs) {
        const auto& internalEdges0 = superbubbles[p.first].internalEdges;
        const auto& internalEdges1 = superbubbles[p.second].internalEdges;

        commonEdges.clear();
        std::set_intersection(
            internalEdges0.begin(), internalEdges0.end(),
            internalEdges1.begin(), internalEdges1.end(),
            back_inserter(commonEdges),
            orderById);

        if(commonEdges.size() == internalEdges0.size()) {
            superbubblesToBeRemoved.push_back(p.first);
        } else if(commonEdges.size() == internalEdges1.size()) {
            superbubblesToBeRemoved.push_back(p.second);
        } else {
            const uint64_t superbubbleId0 = p.first;
            const uint64_t superbubbleId1 = p.second;
            cout << "Superbubbles " << superbubbleId0 << " and " << superbubbleId1 << " intersect." << endl;
            cout << "Superbubbles " << superbubbleId0 << " has " << internalEdges0.size() << " internal edges:";
            for(const edge_descriptor e: internalEdges0) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
            cout << "Superbubbles " << superbubbleId1 << " has " << internalEdges1.size() << " internal edges:";
            for(const edge_descriptor e: internalEdges1) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << endl;
            cout << "Found " << commonEdges.size() << " common edges." << endl;
            DINARA_ASSERT(0);
        }
    }
    deduplicate(superbubblesToBeRemoved);

    vector<Shasta2Superbubble> filtered;
    filtered.reserve(superbubbles.size());
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(!binary_search(
            superbubblesToBeRemoved.begin(),
            superbubblesToBeRemoved.end(),
            superbubbleId)) {
            filtered.push_back(superbubbles[superbubbleId]);
        }
    }
    superbubbles.swap(filtered);
}

void Shasta2AssemblyGraph::findSuperbubbleChains(
    const vector<Shasta2Superbubble>& superbubbles,
    vector<Shasta2SuperbubbleChain>& superbubbleChains) const
{
    std::map<vertex_descriptor, vector<uint64_t> > mapBySource;
    std::map<vertex_descriptor, vector<uint64_t> > mapByTarget;
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        const Shasta2Superbubble& superbubble = superbubbles[superbubbleId];
        mapBySource[superbubble.sourceVertex].push_back(superbubbleId);
        mapByTarget[superbubble.targetVertex].push_back(superbubbleId);
    }

    for(const auto& p: mapBySource) {
        DINARA_ASSERT(p.second.size() == 1);
    }
    for(const auto& p: mapByTarget) {
        DINARA_ASSERT(p.second.size() == 1);
    }

    vector<bool> wasUsed(superbubbles.size(), false);
    vector<uint64_t> forward;
    vector<uint64_t> backward;

    superbubbleChains.clear();
    for(uint64_t superbubbleId=0; superbubbleId<superbubbles.size(); superbubbleId++) {
        if(wasUsed[superbubbleId]) {
            continue;
        }

        forward.clear();
        vertex_descriptor v = superbubbles[superbubbleId].targetVertex;
        while(true) {
            const auto it = mapBySource.find(v);
            if(it == mapBySource.end()) {
                break;
            }
            const vector<uint64_t>& nextVector = it->second;
            DINARA_ASSERT(nextVector.size() == 1);
            const uint64_t nextSuperbubbleId = nextVector.front();
            forward.push_back(nextSuperbubbleId);
            v = superbubbles[nextSuperbubbleId].targetVertex;
        }

        backward.clear();
        v = superbubbles[superbubbleId].sourceVertex;
        while(true) {
            const auto it = mapByTarget.find(v);
            if(it == mapByTarget.end()) {
                break;
            }
            const vector<uint64_t>& previousVector = it->second;
            DINARA_ASSERT(previousVector.size() == 1);
            const uint64_t previousSuperbubbleId = previousVector.front();
            backward.push_back(previousSuperbubbleId);
            v = superbubbles[previousSuperbubbleId].sourceVertex;
        }

        superbubbleChains.emplace_back();
        Shasta2SuperbubbleChain& superbubbleChain = superbubbleChains.back();
        reverse(backward.begin(), backward.end());
        for(const uint64_t id: backward) {
            wasUsed[id] = true;
            superbubbleChain.push_back(superbubbles[id]);
        }
        superbubbleChain.push_back(superbubbles[superbubbleId]);
        wasUsed[superbubbleId] = true;
        for(const uint64_t id: forward) {
            wasUsed[id] = true;
            superbubbleChain.push_back(superbubbles[id]);
        }
    }
}



// Find the non-trivial strongly connected components.
// Each component is stored with vertices sorted to permit binary searches.
void Shasta2AssemblyGraph::findStrongComponents(
    vector< vector<vertex_descriptor> >& strongComponents) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;

    // Map the vertices to integers.
    uint64_t vertexIndex = 0;
    std::map<vertex_descriptor, uint64_t> vertexMap;
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        vertexMap.insert({v, vertexIndex++});
    }

    // Compute strong components.
    std::map<vertex_descriptor, uint64_t> componentMap;
    boost::strong_components(
        assemblyGraph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexMap)));

    // Gather the vertices in each strong component.
    std::map<uint64_t, vector<vertex_descriptor> > componentVertices;
    for(const auto& p: componentMap) {
        componentVertices[p.second].push_back(p.first);
    }

    strongComponents.clear();
    for(const auto& p: componentVertices) {
        const vector<vertex_descriptor>& component = p.second;
        if(component.size() > 1) {
            strongComponents.push_back(component);
            sort(strongComponents.back().begin(), strongComponents.back().end());
        }
    }
}



// This creates a csv file that can be loaded in bandage to see
// the strongly connected components.
void Shasta2AssemblyGraph::colorStrongComponents() const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;

    vector< vector<vertex_descriptor> > strongComponents;
    findStrongComponents(strongComponents);

    ofstream csv("StrongComponents.csv");
    csv << "Id,Color,Component\n";

    // Loop over the non-trivial strongly connected components.
    for(uint64_t i=0; i<strongComponents.size(); i++) {
        const vector<vertex_descriptor>& strongComponent = strongComponents[i];

        for(const vertex_descriptor v0: strongComponent) {
            BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, Shasta2AssemblyGraph) {
                const vertex_descriptor v1 = target(e, assemblyGraph);
                if(std::binary_search(strongComponent.begin(), strongComponent.end(), v1)) {
                    csv << assemblyGraph[e].id << ",Green," << i << "\n";
                }
            }
        }
    }
}

void Shasta2AssemblyGraph::writeSuperbubbles(
    const vector<Shasta2Superbubble>& superbubbles,
    const string& fileName) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;
    ofstream csv(fileName);
    csv << "Superbubble id,Type,Source ploidy,Target ploidy,Internal edges\n";

    for(uint64_t id=0; id<superbubbles.size(); id++) {
        const Shasta2Superbubble& superbubble = superbubbles[id];
        csv << id << ",";

        if(superbubble.isBubble()) {
            const uint64_t ploidy = superbubble.ploidy();
            if(ploidy == 1) {
                csv << "Edge,";
            } else {
                csv << "Bubble,";
            }
        } else {
            csv << "Superbubble,";
        }

        csv << superbubble.sourcePloidy() << ",";
        csv << superbubble.targetPloidy() << ",";
        for(const edge_descriptor e: superbubble.internalEdges) {
            csv << assemblyGraph[e].id << ",";
        }
        csv << "\n";
    }
}

void Shasta2AssemblyGraph::writeSuperbubblesForBandage(
    const vector<Shasta2Superbubble>& superbubbles,
    const string& fileName) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;
    ofstream csv(fileName);
    csv << "Segment,Color\n";
    for(uint64_t id=0; id<superbubbles.size(); id++) {
        const Shasta2Superbubble& superbubble = superbubbles[id];
        const string color = shasta2RandomHslColor(id, 0.75, 0.5);
        for(const edge_descriptor e: superbubble.internalEdges) {
            csv << assemblyGraph[e].id << "," << color << "\n";
        }
    }
}

void Shasta2AssemblyGraph::writeSuperbubbleChains(
    const vector<Shasta2SuperbubbleChain>& superbubbleChains,
    const string& fileName) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;
    ofstream csv(fileName);
    csv << "ChainId,Position,Type,Source ploidy,Target ploidy,"
        "Internal vertices count,Internal edges count,Internal edges\n";

    for(uint64_t chainId=0; chainId<superbubbleChains.size(); chainId++) {
        const Shasta2SuperbubbleChain& superbubbleChain = superbubbleChains[chainId];
        for(uint64_t position=0; position<superbubbleChain.size(); position++) {
            const Shasta2Superbubble& superbubble = superbubbleChain[position];
            csv << chainId << ",";
            csv << position << ",";

            if(superbubble.isBubble()) {
                const uint64_t ploidy = superbubble.ploidy();
                if(ploidy == 1) {
                    csv << "Edge,";
                } else {
                    csv << "Bubble,";
                }
            } else {
                csv << "Superbubble,";
            }

            csv << superbubble.sourcePloidy() << ",";
            csv << superbubble.targetPloidy() << ",";
            csv << superbubble.internalVertices.size() << ",";
            csv << superbubble.internalEdges.size() << ",";
            for(const edge_descriptor e: superbubble.internalEdges) {
                csv << assemblyGraph[e].id << ",";
            }
            csv << "\n";
        }
    }
}

void Shasta2AssemblyGraph::writeSuperbubbleChainsForBandage(
    const vector<Shasta2SuperbubbleChain>& superbubbleChains,
    const string& fileName) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;
    ofstream csv(fileName);
    csv << "Segment,Color,Chain,Position\n";

    for(uint64_t chainId=0; chainId<superbubbleChains.size(); chainId++) {
        const Shasta2SuperbubbleChain& superbubbleChain = superbubbleChains[chainId];
        const string color = shasta2RandomHslColor(chainId, 0.75, 0.5);

        for(uint64_t position=0; position<superbubbleChain.size(); position++) {
            const Shasta2Superbubble& superbubble = superbubbleChain[position];
            for(const edge_descriptor e: superbubble.internalEdges) {
                csv <<
                    assemblyGraph[e].id << "," <<
                    color << "," <<
                    chainId << "," <<
                    position << "\n";
            }
        }
    }
}
