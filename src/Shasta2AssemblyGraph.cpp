#include <boost/pending/disjoint_sets.hpp>

#include "Shasta2AssemblyGraph.hpp"
#include "Shasta2ReadFollowing.hpp"
#include "Shasta2AreSimilarSequences.hpp"
#include "Shasta2RestrictedAnchorGraph.hpp"
#include "Shasta2DisjointSets.hpp"
#include "Shasta2LocalAssembly4.hpp"
#include "Shasta2Superbubble.hpp"
#include "Shasta2SuperbubbleChain.hpp"
#include "Shasta2TangleMatrix1.hpp"
#include "performanceLog.hpp"
#include "timestamp.hpp"
#include "DINARA_ASSERT.hpp"
#include "MultithreadedObject.tpp"
namespace dinara {
    template class MultithreadedObject<Shasta2AssemblyGraph>;
}
#include "deduplicate.hpp"
#include "findConvergingVertex.hpp"
#include "findLinearChains.hpp"
#include "MurmurHash2.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/strong_components.hpp>

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <limits>
#include <numeric>
#include <queue>
#include <ranges>
#include <set>
#include <sstream>
#include <unordered_set>

using namespace dinara;
using namespace std;

namespace {

uint64_t computeN50(vector<uint64_t> lengths)
{
    if(lengths.empty()) {
        return 0;
    }

    sort(lengths.begin(), lengths.end(), greater<uint64_t>());
    const uint64_t total = accumulate(lengths.begin(), lengths.end(), uint64_t(0));
    const uint64_t threshold = (total + 1) / 2;

    uint64_t sum = 0;
    for(const uint64_t length: lengths) {
        sum += length;
        if(sum >= threshold) {
            return length;
        }
    }

    return 0;
}

class Shasta2AnchorGraphEdgePredicate {
public:
    bool operator()(const Shasta2AnchorGraph::edge_descriptor& e) const
    {
        return (*anchorGraph)[e].useForAssembly;
    }

    explicit Shasta2AnchorGraphEdgePredicate(const Shasta2AnchorGraph& anchorGraph) :
        anchorGraph(&anchorGraph)
    {}

    Shasta2AnchorGraphEdgePredicate() = default;

private:
    const Shasta2AnchorGraph* anchorGraph = nullptr;
};

uint64_t countCommonOrientedReadIds(
    const vector<OrientedReadId>& x,
    const vector<OrientedReadId>& y)
{
    uint64_t n = 0;
    auto xIt = x.begin();
    auto yIt = y.begin();
    while((xIt != x.end()) && (yIt != y.end())) {
        if(*xIt < *yIt) {
            ++xIt;
        } else if(*yIt < *xIt) {
            ++yIt;
        } else {
            ++n;
            ++xIt;
            ++yIt;
        }
    }
    return n;
}

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

} // namespace

uint64_t Shasta2AssemblyGraphEdge::offset() const
{
    uint64_t x = 0;
    for(const auto& step: *this) {
        x += step.offset;
    }
    return x;
}

void Shasta2AssemblyGraphEdge::getSequence(vector<dinara::Base>& sequence) const
{
    sequence.clear();
    for(const auto& step: *this) {
        copy(step.sequence.begin(), step.sequence.end(), back_inserter(sequence));
    }
}

uint64_t Shasta2AssemblyGraphEdge::sequenceLength() const
{
    DINARA_ASSERT(wasAssembled);

    uint64_t length = 0;
    for(const auto& step: *this) {
        length += step.sequence.size();
    }
    return length;
}

uint64_t Shasta2AssemblyGraphEdge::length() const
{
    return wasAssembled ? sequenceLength() : offset();
}

double Shasta2AssemblyGraphEdge::averageCoverage() const
{
    uint64_t sum = 0;
    for(const Shasta2AssemblyGraphEdgeStep& step: *this) {
        sum += step.anchorPair.orientedReadIds.size();
    }

    return double(sum) / double(size());
}

double Shasta2AssemblyGraphEdge::lengthWeightedAverageCoverage() const
{
    uint64_t sum0 = 0;
    uint64_t sum1 = 0;
    for(const Shasta2AssemblyGraphEdgeStep& step: *this) {
        const uint64_t length = wasAssembled ? step.sequence.size() : step.offset;
        const uint64_t coverage = step.anchorPair.orientedReadIds.size();
        sum0 += length;
        sum1 += length * coverage;
    }

    return double(sum1) / double(sum0);
}

double Shasta2AssemblyGraphEdgeAverageCoverage(const Shasta2AssemblyGraphEdge& edge)
{
    uint64_t sum = 0;
    for(const Shasta2AssemblyGraphEdgeStep& step: edge) {
        sum += step.anchorPair.orientedReadIds.size();
    }

    return double(sum) / double(edge.size());
}

double Shasta2AssemblyGraphEdgeLengthWeightedAverageCoverage(const Shasta2AssemblyGraphEdge& edge)
{
    uint64_t sum0 = 0;
    uint64_t sum1 = 0;
    for(const Shasta2AssemblyGraphEdgeStep& step: edge) {
        const uint64_t length = edge.wasAssembled ? step.sequence.size() : step.offset;
        const uint64_t coverage = step.anchorPair.orientedReadIds.size();
        sum0 += length;
        sum1 += length * coverage;
    }

    return double(sum1) / double(sum0);
}

Shasta2AssemblyGraph::Shasta2AssemblyGraph(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const Shasta2AnchorGraph& anchorGraph,
    const Shasta2AssemblyGraphOptions& options) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<Shasta2AssemblyGraph>(*this),
    orderById(*this),
    anchorsPointer(&anchors),
    journeysPointer(&journeys),
    options(options)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;

    // Create a filtered AnchorGraph containing only edges marked useForAssembly.
    using FilteredAnchorGraph = boost::filtered_graph<
        Shasta2AnchorGraph,
        Shasta2AnchorGraphEdgePredicate>;
    const FilteredAnchorGraph filteredAnchorGraph(
        anchorGraph,
        Shasta2AnchorGraphEdgePredicate(anchorGraph));

    // Find linear chains of edges in the filtered AnchorGraph.
    vector< list<Shasta2AnchorGraph::edge_descriptor> > chains;
    findLinearChains(filteredAnchorGraph, 1, chains);

    // Generate vertices at chain endpoints.
    map<Shasta2AnchorId, vertex_descriptor> vertexMap;
    for(const auto& chain: chains) {
        if(chain.empty()) {
            continue;
        }

        const Shasta2AnchorId anchorId0 = anchorGraph[chain.front()].anchorPair.anchorIdA;
        const Shasta2AnchorId anchorId1 = anchorGraph[chain.back()].anchorPair.anchorIdB;

        if(!vertexMap.contains(anchorId0)) {
            const vertex_descriptor v0 = add_vertex(
                Shasta2AssemblyGraphVertex(anchorId0, nextVertexId++),
                assemblyGraph);
            vertexMap.insert(make_pair(anchorId0, v0));
        }

        if(!vertexMap.contains(anchorId1)) {
            const vertex_descriptor v1 = add_vertex(
                Shasta2AssemblyGraphVertex(anchorId1, nextVertexId++),
                assemblyGraph);
            vertexMap.insert(make_pair(anchorId1, v1));
        }
    }
    DINARA_ASSERT(vertexMap.size() == num_vertices(assemblyGraph));

    // Generate one AssemblyGraph edge for each chain.
    for(const auto& chain: chains) {
        if(chain.empty()) {
            continue;
        }

        const Shasta2AnchorId anchorId0 = anchorGraph[chain.front()].anchorPair.anchorIdA;
        const Shasta2AnchorId anchorId1 = anchorGraph[chain.back()].anchorPair.anchorIdB;
        const vertex_descriptor v0 = vertexMap.at(anchorId0);
        const vertex_descriptor v1 = vertexMap.at(anchorId1);

        edge_descriptor e;
        tie(e, ignore) = add_edge(
            v0,
            v1,
            Shasta2AssemblyGraphEdge(nextEdgeId++),
            assemblyGraph);
        Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];

        for(const Shasta2AnchorGraph::edge_descriptor eA: chain) {
            const Shasta2AnchorGraphEdge& edgeA = anchorGraph[eA];
            edge.emplace_back(edgeA.anchorPair, edgeA.offset);
        }
    }

    check();

    cout << "The Shasta2AssemblyGraph has " << num_vertices(*this)
        << " vertices and " << num_edges(*this) << " edges." << endl;
}



// Construct from anchor graph with window information.
// Like the base constructor, but also:
// - Copies the anchorToWindow mapping from the anchor graph
// - Stores a pointer to the anchor windows
// - Builds anchorChain for each edge (full chain of AnchorIds like shasta2's Chain)
// - Populates windowSequence for each edge and windowId for each vertex
Shasta2AssemblyGraph::Shasta2AssemblyGraph(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const Shasta2AnchorGraph& anchorGraph,
    const vector<AnchorWindow>& anchorWindows,
    const Shasta2AssemblyGraphOptions& options) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<Shasta2AssemblyGraph>(*this),
    orderById(*this),
    anchorsPointer(&anchors),
    journeysPointer(&journeys),
    options(options),
    anchorWindowsPointer(&anchorWindows)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;

    // Copy the anchor-to-window mapping from the anchor graph.
    anchorToWindow = anchorGraph.anchorToWindow;
    windowCount = anchorGraph.windowCount;

    // Create a filtered AnchorGraph containing only edges marked useForAssembly.
    using FilteredAnchorGraph = boost::filtered_graph<
        Shasta2AnchorGraph,
        Shasta2AnchorGraphEdgePredicate>;
    const FilteredAnchorGraph filteredAnchorGraph(
        anchorGraph,
        Shasta2AnchorGraphEdgePredicate(anchorGraph));

    // Find linear chains of edges in the filtered AnchorGraph.
    vector< list<Shasta2AnchorGraph::edge_descriptor> > chains;
    findLinearChains(filteredAnchorGraph, 1, chains);

    // Generate vertices at chain endpoints.
    map<Shasta2AnchorId, vertex_descriptor> vertexMap;
    for(const auto& chain: chains) {
        if(chain.empty()) {
            continue;
        }

        const Shasta2AnchorId anchorId0 = anchorGraph[chain.front()].anchorPair.anchorIdA;
        const Shasta2AnchorId anchorId1 = anchorGraph[chain.back()].anchorPair.anchorIdB;

        if(!vertexMap.contains(anchorId0)) {
            const vertex_descriptor v0 = add_vertex(
                Shasta2AssemblyGraphVertex(anchorId0, nextVertexId++),
                assemblyGraph);
            vertexMap.insert(make_pair(anchorId0, v0));
        }

        if(!vertexMap.contains(anchorId1)) {
            const vertex_descriptor v1 = add_vertex(
                Shasta2AssemblyGraphVertex(anchorId1, nextVertexId++),
                assemblyGraph);
            vertexMap.insert(make_pair(anchorId1, v1));
        }
    }
    DINARA_ASSERT(vertexMap.size() == num_vertices(assemblyGraph));

    // Generate one AssemblyGraph edge for each chain.
    // Also build the anchorChain (full sequence of AnchorIds).
    for(const auto& chain: chains) {
        if(chain.empty()) {
            continue;
        }

        const Shasta2AnchorId anchorId0 = anchorGraph[chain.front()].anchorPair.anchorIdA;
        const Shasta2AnchorId anchorId1 = anchorGraph[chain.back()].anchorPair.anchorIdB;
        const vertex_descriptor v0 = vertexMap.at(anchorId0);
        const vertex_descriptor v1 = vertexMap.at(anchorId1);

        edge_descriptor e;
        tie(e, ignore) = add_edge(
            v0,
            v1,
            Shasta2AssemblyGraphEdge(nextEdgeId++),
            assemblyGraph);
        Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];

        // Build the anchor chain and steps simultaneously.
        for(const Shasta2AnchorGraph::edge_descriptor eA: chain) {
            const Shasta2AnchorGraphEdge& edgeA = anchorGraph[eA];
            edge.emplace_back(edgeA.anchorPair, edgeA.offset);

            // Add the source anchor of this step to the chain.
            // (The target of the last step is added after the loop.)
            edge.anchorChain.push_back(edgeA.anchorPair.anchorIdA);
        }
        // Add the final anchor (target of the last step).
        edge.anchorChain.push_back(anchorGraph[chain.back()].anchorPair.anchorIdB);
    }

    check();

    // Populate window annotations for all vertices and edges.
    populateWindowAnnotations();

    cout << "The Shasta2AssemblyGraph has " << num_vertices(*this)
        << " vertices and " << num_edges(*this) << " edges"
        << " (with window info: " << windowCount << " windows)." << endl;
}



// Populate anchorChain and windowSequence for all edges,
// and windowId for all vertices, using the current anchorToWindow mapping.
void Shasta2AssemblyGraph::populateWindowAnnotations()
{
    if(windowCount == 0) {
        return;
    }

    Shasta2AssemblyGraph& assemblyGraph = *this;

    // Populate vertex windowId.
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        assemblyGraph[v].windowId = getWindowId(assemblyGraph[v].anchorId);
    }

    // Populate edge anchorChain (if not already set) and windowSequence.
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];

        // Build anchorChain from steps if not already populated.
        if(edge.anchorChain.empty() && !edge.empty()) {
            for(const Shasta2AssemblyGraphEdgeStep& step : edge) {
                edge.anchorChain.push_back(step.anchorPair.anchorIdA);
            }
            edge.anchorChain.push_back(edge.back().anchorPair.anchorIdB);
        }

        // Build windowSequence from anchorChain.
        edge.windowSequence.clear();
        for(const Shasta2AnchorId anchorId : edge.anchorChain) {
            const uint32_t w = getWindowId(anchorId);
            if(w == noWindow) {
                continue;
            }
            if(edge.windowSequence.empty() || edge.windowSequence.back() != w) {
                edge.windowSequence.push_back(w);
            }
        }
    }
}



Shasta2AssemblyGraph::Shasta2AssemblyGraph(
    const Shasta2Anchors& anchors,
    const Shasta2Journeys& journeys,
    const Shasta2AssemblyGraphOptions& options,
    const string& stage) :
    MappedMemoryOwner(anchors),
    MultithreadedObject<Shasta2AssemblyGraph>(*this),
    orderById(*this),
    anchorsPointer(&anchors),
    journeysPointer(&journeys),
    options(options)
{
    loadStage(stage);
}

Shasta2AssemblyGraph::Shasta2AssemblyGraph(
    const MappedMemoryOwner& mappedMemoryOwner,
    const string& name) :
    MappedMemoryOwner(mappedMemoryOwner),
    MultithreadedObject<Shasta2AssemblyGraph>(*this),
    orderById(*this)
{
    load(name);
}

void Shasta2AssemblyGraph::simplifyAndAssemble()
{
    writePerformanceStatistics("Shasta2AssemblyGraph::simplifyAndAssemble begins");

    const uint64_t minComponentN50 = 100000;

    writeIntermediateStageIfRequested("A");

    for(uint64_t iteration=0; iteration<options.simplifyMaxIterationCount; iteration++) {
        uint64_t changeCount = 0;

        changeCount += bubbleCleanup();
        changeCount += compress();
        writeIntermediateStageIfRequested("B" + to_string(iteration));

        changeCount += phaseSuperbubbleChains();
        writeIntermediateStageIfRequested("C" + to_string(iteration));

        if(changeCount == 0) {
            break;
        }
    }

    findAndConnectAssemblyPaths();
    writeIntermediateStageIfRequested("D");

    removeIsolatedVertices();
    removeLowN50Components(minComponentN50);
    writeIntermediateStageIfRequested("E");

    assembleAll();
    write("Final");
    writeFasta("Final");

    writePerformanceStatistics("Shasta2AssemblyGraph::simplifyAndAssemble ends");
}

uint64_t Shasta2AssemblyGraph::bubbleCleanup()
{
    uint64_t modifiedCount = 0;

    vector< pair<vertex_descriptor, vertex_descriptor> > excludeList;
    while(true) {
        const uint64_t modifiedCountThisIteration = bubbleCleanupIteration(excludeList);
        if(modifiedCountThisIteration == 0) {
            break;
        }
        modifiedCount += modifiedCountThisIteration;
    }

    return modifiedCount;
}

void Shasta2AssemblyGraph::findBubbles(vector<Bubble>& bubbles) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;
    bubbles.clear();

    map<vertex_descriptor, vector<edge_descriptor> > m;
    BGL_FORALL_VERTICES(v0, assemblyGraph, Shasta2AssemblyGraph) {
        m.clear();
        BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, Shasta2AssemblyGraph) {
            const vertex_descriptor v1 = target(e, assemblyGraph);
            m[v1].push_back(e);
        }

        for(const auto& p: m) {
            const vertex_descriptor v1 = p.first;
            const vector<edge_descriptor>& edges = p.second;
            if(edges.size() > 1) {
                Bubble bubble;
                bubble.v0 = v0;
                bubble.v1 = v1;
                bubble.edges = edges;
                std::ranges::sort(bubble.edges, orderById);
                bubbles.push_back(bubble);
            }
        }
    }
}

uint64_t Shasta2AssemblyGraph::bubbleCleanupIteration(
    vector< pair<vertex_descriptor, vertex_descriptor> >& excludeList)
{
    performanceLog << timestamp << "Bubble cleanup iteration begins." << endl;
    Shasta2AssemblyGraph& assemblyGraph = *this;

    // Find all bubbles.
    vector<Bubble> allBubbles;
    findBubbles(allBubbles);

    // Find candidate bubbles.
    vector<Bubble> candidateBubbles;
    for(const Bubble& bubble: allBubbles) {
        if(std::ranges::contains(excludeList, make_pair(bubble.v0, bubble.v1))) {
            continue;
        }

        bool hasLongBranch = false;
        for(const edge_descriptor e: bubble.edges) {
            if(assemblyGraph[e].offset() > options.bubbleCleanupMaxBubbleLength) {
                hasLongBranch = true;
                break;
            }
        }

        if(!hasLongBranch) {
            candidateBubbles.push_back(bubble);
        }
    }

    // Assemble sequence for all the edges of these bubbles.
    edgesToBeAssembled.clear();
    for(const Bubble& bubble: candidateBubbles) {
        for(const edge_descriptor e: bubble.edges) {
            if(not assemblyGraph[e].wasAssembled) {
                edgesToBeAssembled.push_back(e);
            }
        }
    }
    assemble();

    uint64_t modifiedCount = 0;
    for(const Bubble& bubble: candidateBubbles) {
        if(bubbleCleanup(bubble)) {
            ++modifiedCount;
        }
    }

    for(const Bubble& bubble: candidateBubbles) {
        excludeList.push_back(make_pair(bubble.v0, bubble.v1));
    }
    std::ranges::sort(excludeList);

    performanceLog << timestamp << "Bubble cleanup iteration ends." << endl;
    return modifiedCount;
}

bool Shasta2AssemblyGraph::analyzeBubble(
    const Bubble& bubble,
    const vector<uint64_t> minRepeatCount,
    vector< pair<uint64_t, uint64_t> >& similarPairs) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;
    using dinara::Base;

    DINARA_ASSERT(bubble.edges.size() > 1);

    // Gather the assembled sequences of all bubble branches.
    vector< vector<Base> > sequences;
    for(const edge_descriptor e: bubble.edges) {
        sequences.emplace_back();
        DINARA_ASSERT(assemblyGraph[e].wasAssembled);
        assemblyGraph[e].getSequence(sequences.back());
    }

    // Loop over pairs of bubble branches.
    ostream html(0);
    similarPairs.clear();
    for(uint64_t i0=0; i0+1<bubble.edges.size(); i0++) {
        const vector<Base>& sequence0 = sequences[i0];
        for(uint64_t i1=i0+1; i1<bubble.edges.size(); i1++) {
            const vector<Base>& sequence1 = sequences[i1];
            if(shasta2AreSimilarSequences(sequence0, sequence1, minRepeatCount, html)) {
                similarPairs.emplace_back(i0, i1);
            }
        }
    }

    return false;
}

bool Shasta2AssemblyGraph::bubbleCleanup(const Bubble& bubble)
{
    const vector<uint64_t> minRepeatCount = {0, 2, 2, 2, 2, 2, 2};

    DINARA_ASSERT(anchorsPointer);
    const Shasta2Anchors& anchors = *anchorsPointer;

    Shasta2AssemblyGraph& assemblyGraph = *this;
    const uint64_t ploidy = bubble.edges.size();
    DINARA_ASSERT(ploidy > 1);

    // Find similar branch pairs by comparing assembled branch sequences.
    vector< pair<uint64_t, uint64_t> > similarPairs;
    analyzeBubble(bubble, minRepeatCount, similarPairs);
    if(similarPairs.empty()) {
        return false;
    }

    // Find oriented reads that appear in each branch.
    vector< vector<OrientedReadId> > allOrientedReadIds(ploidy);
    for(uint64_t i=0; i<ploidy; i++) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[bubble.edges[i]];
        for(const auto& step: edge) {
            copy(
                step.anchorPair.orientedReadIds.begin(),
                step.anchorPair.orientedReadIds.end(),
                back_inserter(allOrientedReadIds[i]));
        }
        deduplicate(allOrientedReadIds[i]);
    }

    // Build the bubble anchor pair and usable oriented reads per branch.
    const Shasta2AnchorId anchorId0 = assemblyGraph[bubble.v0].anchorId;
    const Shasta2AnchorId anchorId1 = assemblyGraph[bubble.v1].anchorId;
    const Shasta2AnchorPair bubbleAnchorPair(anchors, anchorId0, anchorId1, false);

    vector< vector<OrientedReadId> > usableOrientedReadIds(ploidy);
    for(uint64_t i=0; i<ploidy; i++) {
        set_intersection(
            bubbleAnchorPair.orientedReadIds.begin(),
            bubbleAnchorPair.orientedReadIds.end(),
            allOrientedReadIds[i].begin(),
            allOrientedReadIds[i].end(),
            back_inserter(usableOrientedReadIds[i]));
    }

    // Compute groups of similar branches.
    vector< vector<uint64_t> > branchGroups;
    if(ploidy == 2) {
        DINARA_ASSERT(similarPairs.size() == 1);
        branchGroups.push_back({0, 1});
    } else if(similarPairs.size() == (ploidy * (ploidy - 1)) / 2) {
        branchGroups.resize(1);
        for(uint64_t i=0; i<ploidy; i++) {
            branchGroups.front().push_back(i);
        }
    } else {
        vector<uint64_t> rank(ploidy);
        vector<uint64_t> parent(ploidy);
        boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
        for(uint64_t i=0; i<ploidy; i++) {
            disjointSets.make_set(i);
        }

        for(const auto& p: similarPairs) {
            disjointSets.union_set(p.first, p.second);
        }

        vector< vector<uint64_t> > groups(ploidy);
        for(uint64_t i=0; i<ploidy; i++) {
            groups[disjointSets.find_set(i)].push_back(i);
        }

        for(const auto& group: groups) {
            if(not group.empty()) {
                branchGroups.push_back(group);
            }
        }
    }

    for(const auto& branchGroup: branchGroups) {
        if(branchGroup.size() == 1) {
            continue;
        }

        Shasta2AnchorPair newAnchorPair;
        newAnchorPair.anchorIdA = anchorId0;
        newAnchorPair.anchorIdB = anchorId1;
        for(const uint64_t i: branchGroup) {
            copy(
                usableOrientedReadIds[i].begin(),
                usableOrientedReadIds[i].end(),
                back_inserter(newAnchorPair.orientedReadIds));
        }
        deduplicate(newAnchorPair.orientedReadIds);

        if(newAnchorPair.size() < options.bubbleCleanupMinCommonCount) {
            continue;
        }

        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(
            bubble.v0,
            bubble.v1,
            Shasta2AssemblyGraphEdge(nextEdgeId++),
            assemblyGraph);
        Shasta2AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
        edgeNew.emplace_back(newAnchorPair, newAnchorPair.getAverageOffset(anchors));
        edgeNew.anchorChain = {anchorId0, anchorId1};
        // windowSequence will be populated by populateWindowAnnotations if needed.
        if(windowCount > 0) {
            const uint32_t w0 = getWindowId(anchorId0);
            const uint32_t w1 = getWindowId(anchorId1);
            if(w0 != noWindow) edgeNew.windowSequence.push_back(w0);
            if(w1 != noWindow && (edgeNew.windowSequence.empty() || edgeNew.windowSequence.back() != w1)) {
                edgeNew.windowSequence.push_back(w1);
            }
        }

        for(const uint64_t i: branchGroup) {
            boost::remove_edge(bubble.edges[i], assemblyGraph);
        }
    }

    DINARA_ASSERT(out_degree(bubble.v0, assemblyGraph) > 0);
    DINARA_ASSERT(in_degree(bubble.v1, assemblyGraph) > 0);

    return true;
}



// Compute oriented read journeys in the Shasta2AssemblyGraph.
void Shasta2AssemblyGraph::computeJourneys()
{
    DINARA_ASSERT(anchorsPointer);
    DINARA_ASSERT(journeysPointer);
    const Shasta2Anchors& anchors = *anchorsPointer;
    const Shasta2Journeys& journeys = *journeysPointer;

    const uint64_t orientedReadCount = journeys.size();
    Shasta2AssemblyGraph& assemblyGraph = *this;

    // Compute uncompressed journeys for all oriented reads.
    class AssemblyGraphJourneyEntry {
    public:
        edge_descriptor e;
        uint64_t stepId;
        uint64_t positionInJourneyA;
        uint64_t positionInJourneyB;

        bool operator<(const AssemblyGraphJourneyEntry& that) const
        {
            return positionInJourneyA + positionInJourneyB <
                that.positionInJourneyA + that.positionInJourneyB;
        }
    };
    vector< vector<AssemblyGraphJourneyEntry> > assemblyGraphJourneys(orientedReadCount);

    // Loop over Shasta2AssemblyGraph edges.
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];

        // Loop over steps of this edge.
        for(uint64_t stepId=0; stepId<edge.size(); stepId++) {
            const Shasta2AssemblyGraphEdgeStep& step = edge[stepId];

            // Locate the anchors for this step.
            const Shasta2AnchorPair& anchorPair = step.anchorPair;
            const Shasta2AnchorId anchorIdA = anchorPair.anchorIdA;
            const Shasta2AnchorId anchorIdB = anchorPair.anchorIdB;
            const Shasta2Anchor anchorA = anchors[anchorIdA];
            const Shasta2Anchor anchorB = anchors[anchorIdB];

            // Loop over OrientedReadIds of this step.
            auto itA = anchorA.begin();
            auto itB = anchorB.begin();
            for(const OrientedReadId orientedReadId: anchorPair.orientedReadIds) {

                // Locate this OrientedReadId in the two anchors.
                for(; (itA != anchorA.end()) and (itA->orientedReadId != orientedReadId); ++itA) {}
                DINARA_ASSERT(itA != anchorA.end());
                DINARA_ASSERT(itA->orientedReadId == orientedReadId);
                const Shasta2AnchorMarkerInfo& infoA = *itA;
                for(; (itB != anchorB.end()) and (itB->orientedReadId != orientedReadId); ++itB) {}
                DINARA_ASSERT(itB != anchorB.end());
                const Shasta2AnchorMarkerInfo& infoB = *itB;
                DINARA_ASSERT(itB->orientedReadId == orientedReadId);

                const uint32_t positionInJourneyA = infoA.positionInJourney;
                const uint32_t positionInJourneyB = infoB.positionInJourney;

                AssemblyGraphJourneyEntry entry;
                entry.e = e;
                entry.stepId = stepId;
                entry.positionInJourneyA = positionInJourneyA;
                entry.positionInJourneyB = positionInJourneyB;

                assemblyGraphJourneys[orientedReadId.getValue()].push_back(entry);
            }
        }
    }

    // Sort the journeys.
    for(vector<AssemblyGraphJourneyEntry>& v: assemblyGraphJourneys) {
        sort(v.begin(), v.end());
    }

    // Create the compressed journeys.
    // Here, we only consider transitions between Shasta2AssemblyGraph edges.
    // Compressed journeys consisting of only one edge are considered empty.
    compressedJourneys.clear();
    compressedJourneys.resize(orientedReadCount);
    for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
        const vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney = assemblyGraphJourneys[orientedReadIdValue];
        vector<edge_descriptor>& compressedJourney = compressedJourneys[orientedReadIdValue];

        for(uint64_t i1=0; i1<assemblyGraphJourney.size(); i1++) {
            const AssemblyGraphJourneyEntry& entry1 = assemblyGraphJourney[i1];
            const edge_descriptor e1 = entry1.e;

            if(i1 == 0) {
                compressedJourney.push_back(e1);
            } else {
                const uint64_t i0 = i1 - 1;
                const AssemblyGraphJourneyEntry& entry0 = assemblyGraphJourney[i0];
                const edge_descriptor e0 = entry0.e;
                if(e1 != e0) {
                    compressedJourney.push_back(e1);
                }
            }
        }
        if(compressedJourney.size() == 1) {
            compressedJourney.clear();
        }
    }

    // Write out the compressed journeys.
    {
        ofstream csv("AssemblyGraphCompressedJourneys.csv");
        for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
            const vector<edge_descriptor>& compressedJourney = compressedJourneys[orientedReadIdValue];
            csv << orientedReadId << ",";
            for(const edge_descriptor e: compressedJourney) {
                csv << assemblyGraph[e].id << ",";
            }
            csv << "\n";
        }
    }

    {
        ofstream csv("AssemblyGraphJourneys.csv");
        csv << "OrientedReadId,Segment,Step,PositionInJourneyA,PositionInJourneyB\n";
        for(ReadId orientedReadIdValue=0; orientedReadIdValue<orientedReadCount; orientedReadIdValue++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(orientedReadIdValue);
            const vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney = assemblyGraphJourneys[orientedReadIdValue];
            for(const AssemblyGraphJourneyEntry& entry: assemblyGraphJourney) {
                csv << orientedReadId << ",";
                csv << assemblyGraph[entry.e].id << ",";
                csv << entry.stepId << ",";
                csv << entry.positionInJourneyA << ",";
                csv << entry.positionInJourneyB << "\n";
            }
        }
    }
}


uint64_t Shasta2AssemblyGraph::compress()
{
    Shasta2AssemblyGraph& assemblyGraph = *this;
    uint64_t compressCount = 0;

    vector< list<edge_descriptor> > chains;
    findLinearChains(assemblyGraph, 2, chains);

    for(const auto& chain: chains) {
        DINARA_ASSERT(chain.size() > 1);

        const edge_descriptor e0 = chain.front();
        const edge_descriptor e1 = chain.back();
        const vertex_descriptor v0 = source(e0, assemblyGraph);
        const vertex_descriptor v1 = target(e1, assemblyGraph);

        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(
            v0,
            v1,
            Shasta2AssemblyGraphEdge(nextEdgeId++),
            assemblyGraph);
        Shasta2AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];

        // Concatenate the steps and anchor chains of all the edges in the chain.
        bool firstEdge = true;
        for(const edge_descriptor e: chain) {
            const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
            copy(edge.begin(), edge.end(), back_inserter(edgeNew));

            // Merge anchor chains: skip the first anchor of subsequent edges
            // (it's the same as the last anchor of the previous edge).
            if(!edge.anchorChain.empty()) {
                if(firstEdge) {
                    edgeNew.anchorChain = edge.anchorChain;
                } else {
                    // Skip first element (duplicate of previous chain's last).
                    edgeNew.anchorChain.insert(
                        edgeNew.anchorChain.end(),
                        edge.anchorChain.begin() + 1,
                        edge.anchorChain.end());
                }
            }

            // Merge window sequences: append, removing consecutive duplicates.
            for(const uint32_t w : edge.windowSequence) {
                if(edgeNew.windowSequence.empty() || edgeNew.windowSequence.back() != w) {
                    edgeNew.windowSequence.push_back(w);
                }
            }

            firstEdge = false;
        }

        if(compressDebugLevel >= 1) {
            cout << "Compress " << assemblyGraph[chain.front()].id << "..." <<
                assemblyGraph[chain.back()].id <<
                " into " << edgeNew.id << endl;
        }
        if(compressDebugLevel >= 2) {
            cout << "Compress";
            for(const edge_descriptor e: chain) {
                cout << " " << assemblyGraph[e].id;
            }
            cout << " into " << edgeNew.id << endl;
        }
        if(compressDebugLevel >= 3) {
            uint64_t stepCount = 0;
            for(const edge_descriptor e: chain) {
                const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
                const uint64_t stepBegin = stepCount;
                const uint64_t stepEnd = stepBegin + edge.size();
                cout << edge.id << " " << edge.size() << " steps become " <<
                    edgeNew.id << " steps " << stepBegin << "-" << stepEnd << endl;
                stepCount = stepEnd;
            }
        }

        // Now we can remove the edges of the chain and its internal vertices.
        bool isFirst = true;
        for(const edge_descriptor e: chain) {
            if(isFirst) {
                isFirst = false;
            } else {
                const vertex_descriptor v = source(e, assemblyGraph);
                boost::clear_vertex(v, assemblyGraph);
                boost::remove_vertex(v, assemblyGraph);
            }
        }

        ++compressCount;
    }

    return compressCount;
}



// Remove low-coverage dead-end edges (tips).
// A tip is an edge where one endpoint is a dead end (no other connections
// on that side). Adapted from MBG's tryRemoveTip.
//
// For a tip edge v0->v1 with dead end at v0:
// - The tip edge must have coverage <= maxRemovableCoverage and length <= maxRemovableLength.
// - v1 must have at least one outgoing edge (the side away from the tip)
//   with coverage >= minSafeCoverage, going to a vertex other than v0.
//
// Symmetric logic for dead end at v1 (check v0's incoming edges).
// Remove short tip chains (like hifiasm's asg_arc_cut_tips).
// A tip is removed if its total window count <= maxTipWindows
// AND its total length <= maxTipLength bp. Process shortest first.
// When a tip is removed, its RC mirror chain is also removed
// (like hifiasm's asg_seq_del which deletes both orientations).
uint64_t Shasta2AssemblyGraph::removeShortTips(uint32_t maxTipWindows, uint64_t maxTipLength)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;

    // Build anchorId -> vertex_descriptor map for RC mirror lookup.
    map<Shasta2AnchorId, vertex_descriptor> anchorToVertex;
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        anchorToVertex[assemblyGraph[v].anchorId] = v;
    }

    // Count distinct windows across edges.
    auto countDistinctWindows = [](const vector<Shasta2AssemblyGraphEdge*>& edges) -> uint32_t {
        set<uint32_t> windows;
        for(const auto* edge : edges) {
            for(const uint32_t w : edge->windowSequence) {
                windows.insert(w);
            }
        }
        return uint32_t(windows.size());
    };

    // Walk a tip chain from a dead-end vertex. Returns the edges and vertices.
    // walkForward=true: v has in-degree 0, walk via out-edges.
    // walkForward=false: v has out-degree 0, walk via in-edges.
    auto walkTipChain = [&](vertex_descriptor v, bool walkForward,
        vector<edge_descriptor>& chainEdges,
        vector<vertex_descriptor>& chainVertices)
    {
        chainEdges.clear();
        chainVertices.clear();

        if(walkForward) {
            vertex_descriptor w = v;
            chainVertices.push_back(w);
            while(true) {
                if(out_degree(w, assemblyGraph) != 1) break;
                const edge_descriptor e = *out_edges(w, assemblyGraph).first;
                chainEdges.push_back(e);
                const vertex_descriptor next = target(e, assemblyGraph);
                chainVertices.push_back(next);
                if(in_degree(next, assemblyGraph) != 1) break;
                w = next;
            }
        } else {
            vertex_descriptor w = v;
            chainVertices.push_back(w);
            while(true) {
                if(in_degree(w, assemblyGraph) != 1) break;
                const edge_descriptor e = *in_edges(w, assemblyGraph).first;
                chainEdges.push_back(e);
                const vertex_descriptor prev = source(e, assemblyGraph);
                chainVertices.push_back(prev);
                if(out_degree(prev, assemblyGraph) != 1) break;
                w = prev;
            }
        }
    };

    // Remove a chain (edges + isolated vertices) and its RC mirror.
    // Like hifiasm's asg_seq_del which deletes the node (both orientations)
    // and all arcs to/from it.
    //
    // The last vertex in chainVertices is the branch point — it is NOT
    // an internal node and must not be fully cleared. Only internal
    // vertices (dead end + middle nodes) get all their edges removed.
    auto removeChainAndRc = [&](
        const vector<edge_descriptor>& chainEdges,
        const vector<vertex_descriptor>& chainVertices)
    {
        // Internal vertices = all chain vertices except the last one
        // (the branch point). These correspond to hifiasm's deleted nodes.
        set<vertex_descriptor> internalVertices;
        if(chainVertices.size() > 1) {
            for(uint64_t i = 0; i + 1 < chainVertices.size(); i++) {
                internalVertices.insert(chainVertices[i]);
            }
        } else if(chainVertices.size() == 1) {
            // Single vertex tip — the vertex itself is internal.
            internalVertices.insert(chainVertices[0]);
        }

        // Find RC mirror vertices for internal vertices.
        set<vertex_descriptor> rcInternalVertices;
        for(const vertex_descriptor cv : internalVertices) {
            const Shasta2AnchorId rcAnchorId =
                Shasta2AnchorId(uint64_t(assemblyGraph[cv].anchorId) ^ 1ULL);
            auto it = anchorToVertex.find(rcAnchorId);
            if(it != anchorToVertex.end()) {
                rcInternalVertices.insert(it->second);
            }
        }

        // Remove forward chain edges.
        for(const edge_descriptor e : chainEdges) {
            boost::remove_edge(e, assemblyGraph);
        }

        // Clear all edges from RC internal vertices (like hifiasm's asg_seq_del).
        for(const vertex_descriptor rv : rcInternalVertices) {
            boost::clear_vertex(rv, assemblyGraph);
        }

        // Remove isolated forward internal vertices.
        for(const vertex_descriptor cv : internalVertices) {
            if(in_degree(cv, assemblyGraph) == 0 && out_degree(cv, assemblyGraph) == 0) {
                anchorToVertex.erase(assemblyGraph[cv].anchorId);
                boost::remove_vertex(cv, assemblyGraph);
            }
        }

        // Remove isolated RC internal vertices.
        for(const vertex_descriptor rv : rcInternalVertices) {
            if(in_degree(rv, assemblyGraph) == 0 && out_degree(rv, assemblyGraph) == 0) {
                anchorToVertex.erase(assemblyGraph[rv].anchorId);
                boost::remove_vertex(rv, assemblyGraph);
            }
        }
    };

    // Collect tip candidates.
    struct TipCandidate {
        uint32_t totalWindows;
        uint64_t totalLength;
        vertex_descriptor startVertex;
        bool walkForward;
    };
    vector<TipCandidate> candidates;

    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        // Dead end on incoming side (in-degree 0, out-degree >= 1).
        if(in_degree(v, assemblyGraph) == 0 && out_degree(v, assemblyGraph) >= 1) {
            vector<edge_descriptor> chainEdges;
            vector<vertex_descriptor> chainVertices;
            walkTipChain(v, true, chainEdges, chainVertices);

            if(!chainEdges.empty()) {
                uint64_t totalLength = 0;
                vector<Shasta2AssemblyGraphEdge*> edgePtrs;
                for(const edge_descriptor e : chainEdges) {
                    totalLength += assemblyGraph[e].length();
                    edgePtrs.push_back(&assemblyGraph[e]);
                }
                const uint32_t totalWindows = countDistinctWindows(edgePtrs);
                if(totalWindows <= maxTipWindows && totalLength <= maxTipLength) {
                    candidates.push_back({totalWindows, totalLength, v, true});
                }
            }
        }

        // Dead end on outgoing side (out-degree 0, in-degree >= 1).
        if(out_degree(v, assemblyGraph) == 0 && in_degree(v, assemblyGraph) >= 1) {
            vector<edge_descriptor> chainEdges;
            vector<vertex_descriptor> chainVertices;
            walkTipChain(v, false, chainEdges, chainVertices);

            if(!chainEdges.empty()) {
                uint64_t totalLength = 0;
                vector<Shasta2AssemblyGraphEdge*> edgePtrs;
                for(const edge_descriptor e : chainEdges) {
                    totalLength += assemblyGraph[e].length();
                    edgePtrs.push_back(&assemblyGraph[e]);
                }
                const uint32_t totalWindows = countDistinctWindows(edgePtrs);
                if(totalWindows <= maxTipWindows && totalLength <= maxTipLength) {
                    candidates.push_back({totalWindows, totalLength, v, false});
                }
            }
        }
    }

    // Sort shortest first (like hifiasm's radix_sort_srt64).
    sort(candidates.begin(), candidates.end(),
        [](const TipCandidate& a, const TipCandidate& b) {
            if(a.totalWindows != b.totalWindows) return a.totalWindows < b.totalWindows;
            return a.totalLength < b.totalLength;
        });

    // Process each candidate. Re-walk to verify it's still a tip
    // (earlier removals may have changed the graph).
    uint64_t removedCount = 0;

    for(const TipCandidate& candidate : candidates) {
        const vertex_descriptor v = candidate.startVertex;

        // Re-check dead-end condition.
        if(candidate.walkForward) {
            if(in_degree(v, assemblyGraph) != 0 || out_degree(v, assemblyGraph) == 0) continue;
        } else {
            if(out_degree(v, assemblyGraph) != 0 || in_degree(v, assemblyGraph) == 0) continue;
        }

        // Re-walk and re-check thresholds.
        vector<edge_descriptor> chainEdges;
        vector<vertex_descriptor> chainVertices;
        walkTipChain(v, candidate.walkForward, chainEdges, chainVertices);

        if(chainEdges.empty()) continue;

        uint64_t totalLength = 0;
        vector<Shasta2AssemblyGraphEdge*> edgePtrs;
        for(const edge_descriptor e : chainEdges) {
            totalLength += assemblyGraph[e].length();
            edgePtrs.push_back(&assemblyGraph[e]);
        }
        const uint32_t totalWindows = countDistinctWindows(edgePtrs);
        if(totalWindows > maxTipWindows || totalLength > maxTipLength) continue;

        // Remove the tip chain and its RC mirror.
        removeChainAndRc(chainEdges, chainVertices);
        ++removedCount;
    }

    // Log all edges with their vertex degrees (for debugging).
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        cout << "  Edge " << edge.id
             << ": src=" << assemblyGraph[v0].anchorId
             << "(in=" << in_degree(v0, assemblyGraph)
             << ",out=" << out_degree(v0, assemblyGraph) << ")"
             << " -> tgt=" << assemblyGraph[v1].anchorId
             << "(in=" << in_degree(v1, assemblyGraph)
             << ",out=" << out_degree(v1, assemblyGraph) << ")"
             << " len=" << edge.length()
             << " ws=[";
        for(uint64_t i = 0; i < edge.windowSequence.size(); i++) {
            if(i > 0) cout << ",";
            cout << edge.windowSequence[i];
        }
        cout << "]" << endl;
    }

    // Log remaining tips that were NOT removed (for debugging).
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        bool isDeadEnd = false;
        bool forward = false;
        if(in_degree(v, assemblyGraph) == 0 && out_degree(v, assemblyGraph) >= 1) {
            isDeadEnd = true; forward = true;
        }
        if(out_degree(v, assemblyGraph) == 0 && in_degree(v, assemblyGraph) >= 1) {
            isDeadEnd = true; forward = false;
        }
        if(!isDeadEnd) continue;

        vector<edge_descriptor> chainEdges;
        vector<vertex_descriptor> chainVertices;
        walkTipChain(v, forward, chainEdges, chainVertices);
        if(chainEdges.empty()) continue;

        uint64_t totalLength = 0;
        vector<Shasta2AssemblyGraphEdge*> edgePtrs;
        for(const edge_descriptor e : chainEdges) {
            totalLength += assemblyGraph[e].length();
            edgePtrs.push_back(&assemblyGraph[e]);
        }
        const uint32_t totalWindows = countDistinctWindows(edgePtrs);

        cout << "  Remaining tip: vertex " << assemblyGraph[v].anchorId
             << " direction=" << (forward ? "fwd" : "bwd")
             << " edges=" << chainEdges.size()
             << " windows=" << totalWindows
             << " length=" << totalLength
             << " edgeIds=";
        for(const edge_descriptor e : chainEdges) {
            cout << assemblyGraph[e].id << " ";
        }
        cout << endl;
    }

    cout << "removeShortTips removed " << removedCount
         << " tips with <= " << maxTipWindows << " windows"
         << " and <= " << maxTipLength << " bp." << endl;
    return removedCount;
}



uint64_t Shasta2AssemblyGraph::removeLowCoverageTips(
    double maxRemovableCoverage,
    double minSafeCoverage,
    uint64_t maxRemovableLength)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;
    uint64_t removedCount = 0;

    // Collect edges to remove (can't modify graph while iterating).
    vector<edge_descriptor> edgesToRemove;

    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];

        // Check coverage and length thresholds.
        if(edge.averageCoverage() > maxRemovableCoverage) {
            continue;
        }
        if(edge.length() > maxRemovableLength) {
            continue;
        }

        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        // Check if this is a tip: one end must be a dead end.
        // A dead end at v0 means in_degree(v0)==0 and out_degree(v0)==1.
        // A dead end at v1 means out_degree(v1)==0 and in_degree(v1)==1.
        const bool deadEndAtSource =
            (in_degree(v0, assemblyGraph) == 0) && (out_degree(v0, assemblyGraph) == 1);
        const bool deadEndAtTarget =
            (out_degree(v1, assemblyGraph) == 0) && (in_degree(v1, assemblyGraph) == 1);

        if(!deadEndAtSource && !deadEndAtTarget) {
            continue;
        }

        // MBG's tryRemoveTip checks three things for each neighbor:
        // 1. The neighbor node must have coverage >= minSafeCoverage (hard abort).
        // 2. The transition to the neighbor must have coverage <= maxRemovableCoverage (hard abort).
        // 3. The neighbor must have a safe edge on its other side (not back to tip).
        //
        // In our directed graph:
        // - "neighbor nodes" = edges on the other side of the non-dead-end vertex.
        // - We don't have separate transition coverage, so check 2 is covered by
        //   the tip edge's own coverage check above.
        // - Check 1: ALL edges on the other side must have coverage >= minSafeCoverage.
        // - Check 3: at least one of those edges must lead to a vertex with a safe
        //   connection further into the graph.

        bool canRemove = false;

        if(deadEndAtSource) {
            // v1 is the non-dead-end vertex. The tip arrives on v1's incoming side.
            // "Neighbor nodes" = v1's outgoing edges (the other side).

            // Check 1: ALL outgoing edges from v1 must have coverage >= minSafeCoverage.
            bool allNeighborsSafe = true;
            if(out_degree(v1, assemblyGraph) == 0) {
                allNeighborsSafe = false; // No neighbors on the other side.
            }
            BGL_FORALL_OUTEDGES(v1, eOther, assemblyGraph, Shasta2AssemblyGraph) {
                if(assemblyGraph[eOther].averageCoverage() < minSafeCoverage) {
                    allNeighborsSafe = false;
                    break;
                }
            }

            // Check 3: At least one outgoing edge from v1 must lead to a vertex
            // that has a safe connection further (not back to v0).
            if(allNeighborsSafe) {
                BGL_FORALL_OUTEDGES(v1, eOther, assemblyGraph, Shasta2AssemblyGraph) {
                    const vertex_descriptor vOtherTarget = target(eOther, assemblyGraph);
                    if(vOtherTarget == v0) continue;
                    if(assemblyGraph[eOther].averageCoverage() >= minSafeCoverage) {
                        canRemove = true;
                        break;
                    }
                }
            }
        }

        if(deadEndAtTarget) {
            // v0 is the non-dead-end vertex. The tip leaves from v0's outgoing side.
            // "Neighbor nodes" = v0's incoming edges (the other side).

            // Check 1: ALL incoming edges to v0 must have coverage >= minSafeCoverage.
            bool allNeighborsSafe = true;
            if(in_degree(v0, assemblyGraph) == 0) {
                allNeighborsSafe = false;
            }
            BGL_FORALL_INEDGES(v0, eOther, assemblyGraph, Shasta2AssemblyGraph) {
                if(assemblyGraph[eOther].averageCoverage() < minSafeCoverage) {
                    allNeighborsSafe = false;
                    break;
                }
            }

            // Check 3: At least one incoming edge to v0 must come from a vertex
            // that has a safe connection further (not from v1).
            if(allNeighborsSafe) {
                BGL_FORALL_INEDGES(v0, eOther, assemblyGraph, Shasta2AssemblyGraph) {
                    const vertex_descriptor vOtherSource = source(eOther, assemblyGraph);
                    if(vOtherSource == v1) continue;
                    if(assemblyGraph[eOther].averageCoverage() >= minSafeCoverage) {
                        canRemove = true;
                        break;
                    }
                }
            }
        }

        if(canRemove) {
            edgesToRemove.push_back(e);
        }
    }

    // Remove the collected edges and their dead-end vertices.
    for(const edge_descriptor e: edgesToRemove) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        boost::remove_edge(e, assemblyGraph);
        ++removedCount;

        // Remove the dead-end vertex if it's now isolated.
        if(in_degree(v0, assemblyGraph) == 0 && out_degree(v0, assemblyGraph) == 0) {
            boost::remove_vertex(v0, assemblyGraph);
        }
        if(in_degree(v1, assemblyGraph) == 0 && out_degree(v1, assemblyGraph) == 0) {
            boost::remove_vertex(v1, assemblyGraph);
        }
    }

    cout << "removeLowCoverageTips removed " << removedCount << " tips." << endl;
    return removedCount;
}



// Remove low-coverage crosslink edges at branching vertices.
// Adapted from MBG's tryRemoveCrosslinks.
//
// At a vertex v with out-degree >= 2:
// - There must be at least one outgoing edge with coverage >= minSafeCoverage.
// - A low-coverage outgoing edge e (coverage <= maxRemovableCoverage) is removable if:
//   (a) The target vertex of e has out-degree >= 2 (it has connections on the
//       other side, away from v — MBG's edges[reverse(edge)].size() >= 2).
//   (b) The target vertex has at least one outgoing edge (not back to v)
//       with coverage >= minSafeCoverage.
//
// Symmetric logic for incoming edges at vertices with in-degree >= 2.
uint64_t Shasta2AssemblyGraph::removeLowCoverageCrosslinks(
    double maxRemovableCoverage,
    double minSafeCoverage)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;
    uint64_t removedCount = 0;

    vector<edge_descriptor> edgesToRemove;

    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {

        // Check outgoing edges from v.
        if(out_degree(v, assemblyGraph) >= 2) {
            // There must be at least one safe outgoing edge.
            bool hasSafeOutgoing = false;
            BGL_FORALL_OUTEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
                if(assemblyGraph[e].averageCoverage() >= minSafeCoverage) {
                    hasSafeOutgoing = true;
                    break;
                }
            }

            if(hasSafeOutgoing) {
                BGL_FORALL_OUTEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
                    if(assemblyGraph[e].averageCoverage() > maxRemovableCoverage) {
                        continue;
                    }
                    // The target vertex must have connections on its OTHER side
                    // (outgoing, away from v).
                    const vertex_descriptor vTarget = target(e, assemblyGraph);
                    if(out_degree(vTarget, assemblyGraph) < 2) {
                        continue;
                    }
                    // The target must have a safe outgoing edge (not back to v).
                    bool targetHasSafeOnOtherSide = false;
                    BGL_FORALL_OUTEDGES(vTarget, eOther, assemblyGraph, Shasta2AssemblyGraph) {
                        const vertex_descriptor vOtherTarget = target(eOther, assemblyGraph);
                        if(vOtherTarget == v) continue;
                        if(assemblyGraph[eOther].averageCoverage() >= minSafeCoverage) {
                            targetHasSafeOnOtherSide = true;
                            break;
                        }
                    }
                    if(targetHasSafeOnOtherSide) {
                        edgesToRemove.push_back(e);
                    }
                }
            }
        }

        // Check incoming edges to v.
        if(in_degree(v, assemblyGraph) >= 2) {
            bool hasSafeIncoming = false;
            BGL_FORALL_INEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
                if(assemblyGraph[e].averageCoverage() >= minSafeCoverage) {
                    hasSafeIncoming = true;
                    break;
                }
            }

            if(hasSafeIncoming) {
                BGL_FORALL_INEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
                    if(assemblyGraph[e].averageCoverage() > maxRemovableCoverage) {
                        continue;
                    }
                    // The source vertex must have connections on its OTHER side
                    // (incoming, away from v).
                    const vertex_descriptor vSource = source(e, assemblyGraph);
                    if(in_degree(vSource, assemblyGraph) < 2) {
                        continue;
                    }
                    // The source must have a safe incoming edge (not from v).
                    bool sourceHasSafeOnOtherSide = false;
                    BGL_FORALL_INEDGES(vSource, eOther, assemblyGraph, Shasta2AssemblyGraph) {
                        const vertex_descriptor vOtherSource = source(eOther, assemblyGraph);
                        if(vOtherSource == v) continue;
                        if(assemblyGraph[eOther].averageCoverage() >= minSafeCoverage) {
                            sourceHasSafeOnOtherSide = true;
                            break;
                        }
                    }
                    if(sourceHasSafeOnOtherSide) {
                        edgesToRemove.push_back(e);
                    }
                }
            }
        }
    }

    // Deduplicate (an edge could be flagged from both its source and target).
    std::ranges::sort(edgesToRemove, orderById);
    edgesToRemove.erase(
        std::unique(edgesToRemove.begin(), edgesToRemove.end()),
        edgesToRemove.end());

    for(const edge_descriptor e: edgesToRemove) {
        boost::remove_edge(e, assemblyGraph);
        ++removedCount;
    }

    cout << "removeLowCoverageCrosslinks removed " << removedCount << " crosslinks." << endl;
    return removedCount;
}



// Remove edges and vertices with estimated copy number 0.
// Adapted from MBG's cleanComponent.
//
// Copy number is estimated as round(edgeCoverage / averageCoverage).
// If estimatedAverageCoverage is 0, it is auto-estimated from long edges (>= 50kb).
//
// Unlike a naive global filter, this validates each connected component:
// for every vertex, the sum of incoming edge copy counts must equal the sum
// of outgoing edge copy counts (flow conservation). Only components where
// the copy number model fits cleanly are cleaned.
uint64_t Shasta2AssemblyGraph::cleanByCopyNumber(
    double estimatedAverageCoverage)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;

    // Auto-estimate average coverage from long edges if not provided.
    if(estimatedAverageCoverage == 0.) {
        const uint64_t minLengthForEstimate = 50000;
        double coverageSum = 0.;
        double lengthSum = 0.;
        BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
            const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
            const uint64_t edgeLength = edge.length();
            if(edgeLength >= minLengthForEstimate) {
                coverageSum += edge.averageCoverage() * double(edgeLength);
                lengthSum += double(edgeLength);
            }
        }
        if(lengthSum > 0.) {
            estimatedAverageCoverage = coverageSum / lengthSum;
        } else {
            // Fall back to all edges.
            BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
                const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
                const uint64_t edgeLength = edge.length();
                coverageSum += edge.averageCoverage() * double(edgeLength);
                lengthSum += double(edgeLength);
            }
            if(lengthSum > 0.) {
                estimatedAverageCoverage = coverageSum / lengthSum;
            } else {
                cout << "cleanByCopyNumber: no edges to estimate coverage from." << endl;
                return 0;
            }
        }
    }

    cout << "cleanByCopyNumber: estimated average coverage = " << estimatedAverageCoverage << endl;

    // Helper: compute copy number for a given coverage.
    auto getCopyNumber = [&](double coverage) -> uint64_t {
        return uint64_t((coverage + estimatedAverageCoverage / 2.) / estimatedAverageCoverage);
    };

    // Find connected components (ignoring edge direction).
    // Build vertex index map.
    map<vertex_descriptor, uint64_t> vertexIndexMap;
    vector<vertex_descriptor> vertexTable;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        vertexTable.push_back(v);
        vertexIndexMap[v] = vertexIndex++;
    }

    const uint64_t n = vertexTable.size();
    if(n == 0) {
        cout << "cleanByCopyNumber: no vertices." << endl;
        return 0;
    }

    // Union-Find for connected components.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i = 0; i < n; i++) {
        disjointSets.make_set(i);
    }
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        disjointSets.union_set(vertexIndexMap[v0], vertexIndexMap[v1]);
    }

    // Group vertices by component.
    map<uint64_t, vector<uint64_t>> components;
    for(uint64_t i = 0; i < n; i++) {
        components[disjointSets.find_set(i)].push_back(i);
    }

    // Process each component.
    uint64_t totalEdgesRemoved = 0;
    uint64_t totalVerticesRemoved = 0;

    for(const auto& [componentId, componentVertexIndices]: components) {

        // Collect all edges in this component.
        set<edge_descriptor> componentEdges;
        for(const uint64_t vIdx: componentVertexIndices) {
            const vertex_descriptor v = vertexTable[vIdx];
            BGL_FORALL_OUTEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
                componentEdges.insert(e);
            }
        }

        // Validate the component: check that copy numbers are consistent.
        // For each edge, coverage must fit an integer copy count within ±0.4.
        // For each vertex, sum of incoming copy counts must equal sum of outgoing.
        bool componentValid = true;

        // Check edge copy numbers are well-defined.
        map<edge_descriptor, uint64_t> edgeCopyNumbers;
        for(const edge_descriptor e: componentEdges) {
            const double coverage = assemblyGraph[e].averageCoverage();
            const uint64_t copyNum = getCopyNumber(coverage);
            const double normalizedCov = coverage / estimatedAverageCoverage;
            // Check that coverage is within ±0.4 of the integer copy number.
            if(copyNum <= 1 && std::abs(normalizedCov - double(copyNum)) > 0.4) {
                componentValid = false;
                break;
            }
            edgeCopyNumbers[e] = copyNum;
        }

        if(!componentValid) {
            continue;
        }

        // Check flow conservation at each vertex in the component.
        for(const uint64_t vIdx: componentVertexIndices) {
            const vertex_descriptor v = vertexTable[vIdx];

            uint64_t inCopySum = 0;
            BGL_FORALL_INEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
                if(edgeCopyNumbers.contains(e)) {
                    inCopySum += edgeCopyNumbers[e];
                }
            }

            uint64_t outCopySum = 0;
            BGL_FORALL_OUTEDGES(v, e, assemblyGraph, Shasta2AssemblyGraph) {
                if(edgeCopyNumbers.contains(e)) {
                    outCopySum += edgeCopyNumbers[e];
                }
            }

            // At internal vertices (not sources/sinks), flow must be conserved.
            // At sources (in_degree==0) or sinks (out_degree==0), skip the check.
            if(in_degree(v, assemblyGraph) > 0 && out_degree(v, assemblyGraph) > 0) {
                if(inCopySum != outCopySum) {
                    componentValid = false;
                    break;
                }
            }
        }

        if(!componentValid) {
            continue;
        }

        // Component is valid. Remove edges and vertices with copy number 0.
        vector<edge_descriptor> edgesToRemove;
        for(const auto& [e, copyNum]: edgeCopyNumbers) {
            if(copyNum == 0) {
                edgesToRemove.push_back(e);
            }
        }

        for(const edge_descriptor e: edgesToRemove) {
            boost::remove_edge(e, assemblyGraph);
            ++totalEdgesRemoved;
        }
    }

    // Remove any vertices that became isolated.
    vector<vertex_descriptor> verticesToRemove;
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        if(in_degree(v, assemblyGraph) == 0 && out_degree(v, assemblyGraph) == 0) {
            verticesToRemove.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToRemove) {
        boost::remove_vertex(v, assemblyGraph);
        ++totalVerticesRemoved;
    }

    cout << "cleanByCopyNumber removed " << totalEdgesRemoved
         << " edges and " << totalVerticesRemoved << " vertices." << endl;
    return totalEdgesRemoved;
}



// Remove weak stalks (dead-end linear chains leading to branch points).
// Ported from Shasta2AnchorGraph::cutWeakStalksLeadingToBranch.
//
// A stalk starts at a tip vertex (in-degree 0 or out-degree 0),
// follows a linear chain, and if it reaches a branch point
// with the read union across all edge steps still <= maxReadCount, it is cut.
// If the chain reaches a dead end (both in-degree 0 and out-degree 0),
// it is not cut.
uint64_t Shasta2AssemblyGraph::cutWeakStalks(uint64_t maxReadCount)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;

    // Compute the read union across a set of edges.
    // Returns false (and stops early) if the union exceeds the threshold.
    auto readUnionWithinThreshold = [&](
        const vector<edge_descriptor>& edges,
        uint64_t threshold) -> bool {
        std::unordered_set<uint64_t> readValues;
        readValues.reserve(threshold + 1);
        for(const edge_descriptor e : edges) {
            const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
            for(const Shasta2AssemblyGraphEdgeStep& step : edge) {
                for(const OrientedReadId& orientedReadId : step.anchorPair.orientedReadIds) {
                    readValues.insert(orientedReadId.getValue());
                    if(readValues.size() > threshold) {
                        return false;
                    }
                }
            }
        }
        return true;
    };

    // Try to collect a weak stalk starting from vStart.
    // forward=true: vStart has in-degree 0, walk forward.
    // forward=false: vStart has out-degree 0, walk backward.
    auto tryCollectWeakStalk = [&](
        vertex_descriptor vStart,
        bool forward,
        vector<edge_descriptor>& candidateEdgesToCut)
    {
        vector<edge_descriptor> chainEdges;
        vertex_descriptor current = vStart;

        if(forward) {
            // Source tip: in-degree must be 0, out-degree must be 1.
            if(in_degree(current, assemblyGraph) != 0) return;
            if(out_degree(current, assemblyGraph) != 1) return;
        } else {
            // Sink tip: out-degree must be 0, in-degree must be 1.
            if(out_degree(current, assemblyGraph) != 0) return;
            if(in_degree(current, assemblyGraph) != 1) return;
        }

        bool shouldCut = false;

        while(true) {
            if(forward) {
                const uint64_t outDeg = out_degree(current, assemblyGraph);

                // Dead end — don't cut.
                if(outDeg == 0) {
                    break;
                }

                // Branch point — check if we should cut.
                if(outDeg > 1) {
                    shouldCut = !chainEdges.empty() &&
                        readUnionWithinThreshold(chainEdges, maxReadCount);
                    break;
                }

                // Single out-edge — extend the chain.
                edge_descriptor e;
                BGL_FORALL_OUTEDGES(current, eOut, assemblyGraph, Shasta2AssemblyGraph) {
                    e = eOut;
                }
                const vertex_descriptor next = target(e, assemblyGraph);

                // If next is a merge point (in-degree > 1), cut here.
                if(in_degree(next, assemblyGraph) > 1) {
                    chainEdges.push_back(e);
                    shouldCut = readUnionWithinThreshold(chainEdges, maxReadCount);
                    break;
                }

                chainEdges.push_back(e);

                // Check threshold early.
                if(!readUnionWithinThreshold(chainEdges, maxReadCount)) {
                    break;
                }

                current = next;
            } else {
                const uint64_t inDeg = in_degree(current, assemblyGraph);

                // Dead end — don't cut.
                if(inDeg == 0) {
                    break;
                }

                // Branch point — check if we should cut.
                if(inDeg > 1) {
                    shouldCut = !chainEdges.empty() &&
                        readUnionWithinThreshold(chainEdges, maxReadCount);
                    break;
                }

                // Single in-edge — extend the chain.
                edge_descriptor e;
                BGL_FORALL_INEDGES(current, eIn, assemblyGraph, Shasta2AssemblyGraph) {
                    e = eIn;
                }
                const vertex_descriptor prev = source(e, assemblyGraph);

                // If prev is a fork (out-degree > 1), cut here.
                if(out_degree(prev, assemblyGraph) > 1) {
                    chainEdges.push_back(e);
                    shouldCut = readUnionWithinThreshold(chainEdges, maxReadCount);
                    break;
                }

                chainEdges.push_back(e);

                // Check threshold early.
                if(!readUnionWithinThreshold(chainEdges, maxReadCount)) {
                    break;
                }

                current = prev;
            }
        }

        if(shouldCut) {
            candidateEdgesToCut.insert(
                candidateEdgesToCut.end(),
                chainEdges.begin(),
                chainEdges.end());
        }
    };

    // Collect all weak stalks.
    vector<edge_descriptor> candidateEdgesToCut;
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        tryCollectWeakStalk(v, true, candidateEdgesToCut);
        tryCollectWeakStalk(v, false, candidateEdgesToCut);
    }

    // Remove collected edges and clean up isolated vertices.
    // Use a set to avoid double-removal.
    set<edge_descriptor> edgesToRemove(
        candidateEdgesToCut.begin(), candidateEdgesToCut.end());

    uint64_t cutCount = 0;
    set<vertex_descriptor> affectedVertices;
    for(const edge_descriptor e : edgesToRemove) {
        affectedVertices.insert(source(e, assemblyGraph));
        affectedVertices.insert(target(e, assemblyGraph));
        boost::remove_edge(e, assemblyGraph);
        ++cutCount;
    }

    for(const vertex_descriptor v : affectedVertices) {
        if(in_degree(v, assemblyGraph) == 0 && out_degree(v, assemblyGraph) == 0) {
            boost::remove_vertex(v, assemblyGraph);
        }
    }

    cout << "cutWeakStalks removed " << cutCount
         << " edges (maxReadCount=" << maxReadCount << ")." << endl;
    return cutCount;
}



void Shasta2AssemblyGraph::findAndConnectAssemblyPaths()
{
    writePerformanceStatistics("Shasta2AssemblyGraph::findAndConnectAssemblyPaths begins");

    vector< vector<edge_descriptor> > assemblyPaths;
    findAssemblyPaths(assemblyPaths);
    connectAssemblyPaths(assemblyPaths);

    writePerformanceStatistics("Shasta2AssemblyGraph::findAndConnectAssemblyPaths ends");
}

bool Shasta2AssemblyGraph::canConnect(
    const edge_descriptor e0,
    const edge_descriptor e1) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;
    const vertex_descriptor v0 = target(e0, assemblyGraph);
    const vertex_descriptor v1 = source(e1, assemblyGraph);
    if(v0 == v1) {
        return true;
    }
    const Shasta2AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const Shasta2AnchorId anchorId1 = assemblyGraph[v1].anchorId;
    if(anchorId0 == anchorId1) {
        return true;
    }

    DINARA_ASSERT(anchorsPointer);
    DINARA_ASSERT(journeysPointer);

    ostream html(0);
    const Shasta2TangleMatrix1 tangleMatrix(assemblyGraph, {e0}, {e1}, html);

    try {
        ostream html(0);
        Shasta2RestrictedAnchorGraph restrictedAnchorGraph(
            *anchorsPointer,
            *journeysPointer,
            tangleMatrix,
            0,
            0,
            html);
        vector<Shasta2RestrictedAnchorGraph::edge_descriptor> path;
        restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, path);
        uint64_t minCoverage = std::numeric_limits<uint64_t>::max();
        for(const Shasta2RestrictedAnchorGraph::edge_descriptor e: path) {
            const auto& edge = restrictedAnchorGraph[e];
            minCoverage = min(minCoverage, edge.anchorPair.size());
        }
        if(minCoverage == 0) {
            return false;
        }
    } catch(const std::exception& e) {
        (void)e;
        return false;
    }

    return true;
}

void Shasta2AssemblyGraph::findAssemblyPaths(
    vector< vector<edge_descriptor> >& assemblyPaths) const
{
    Shasta2ReadFollowing::Graph readFollowingGraph(*this);
    readFollowingGraph.findPaths(assemblyPaths);
}

void Shasta2AssemblyGraph::connectAssemblyPaths(
    const vector< vector<edge_descriptor> >& assemblyPaths)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;
    DINARA_ASSERT(anchorsPointer);
    DINARA_ASSERT(journeysPointer);
    const Shasta2Anchors& anchors = *anchorsPointer;
    const Shasta2Journeys& journeys = *journeysPointer;
    const bool debug = false;

    set<edge_descriptor> pathInitialSegments;
    set<edge_descriptor> pathFinalSegments;
    for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {
        DINARA_ASSERT(assemblyPath.size() > 1);
        const edge_descriptor e0 = assemblyPath.front();
        const edge_descriptor e1 = assemblyPath.back();
        DINARA_ASSERT(e0 != e1);
        DINARA_ASSERT(not pathInitialSegments.contains(e0));
        DINARA_ASSERT(not pathFinalSegments.contains(e0));
        DINARA_ASSERT(not pathInitialSegments.contains(e1));
        DINARA_ASSERT(not pathFinalSegments.contains(e1));
        pathInitialSegments.insert(e0);
        pathFinalSegments.insert(e1);
    }

    map<edge_descriptor, uint64_t> pathInternalSegments;
    for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {
        for(uint64_t i=1; i<assemblyPath.size()-1; i++) {
            const edge_descriptor e = assemblyPath[i];
            DINARA_ASSERT(not pathInitialSegments.contains(e));
            DINARA_ASSERT(not pathFinalSegments.contains(e));
            const auto it = pathInternalSegments.find(e);
            if(it == pathInternalSegments.end()) {
                pathInternalSegments.insert({e, 1});
            } else {
                ++it->second;
            }
        }
    }

    for(const vector<edge_descriptor>& assemblyPath: assemblyPaths) {
        vector<edge_descriptor> newAssemblyPath;
        for(uint64_t i=0; i<assemblyPath.size(); i++) {
            const edge_descriptor e = assemblyPath[i];
            const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            const Shasta2AnchorId anchorId0 = assemblyGraph[v0].anchorId;
            const Shasta2AnchorId anchorId1 = assemblyGraph[v1].anchorId;

            vertex_descriptor v0New = null_vertex();
            vertex_descriptor v1New = null_vertex();
            uint64_t idNew = invalid<uint64_t>;
            if(i == 0) {
                v0New = v0;
                v1New = add_vertex(Shasta2AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
                idNew = edge.id;
            } else if(i == assemblyPath.size() - 1) {
                v0New = add_vertex(Shasta2AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
                v1New = v1;
                idNew = edge.id;
            } else {
                v0New = add_vertex(Shasta2AssemblyGraphVertex(anchorId0, nextVertexId++), assemblyGraph);
                v1New = add_vertex(Shasta2AssemblyGraphVertex(anchorId1, nextVertexId++), assemblyGraph);
                if(pathInternalSegments[e] == 1) {
                    idNew = edge.id;
                } else {
                    idNew = nextEdgeId++;
                }
            }

            edge_descriptor eNew;
            tie(eNew, ignore) = add_edge(v0New, v1New, Shasta2AssemblyGraphEdge(idNew), assemblyGraph);
            Shasta2AssemblyGraphEdge& edgeNew = assemblyGraph[eNew];
            edgeNew = edge;
            edgeNew.id = idNew;
            newAssemblyPath.push_back(eNew);
        }

        // Add bridge edges between consecutive segments in the path.
        for(uint64_t i1=1; i1<newAssemblyPath.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const edge_descriptor e0 = newAssemblyPath[i0];
            const edge_descriptor e1 = newAssemblyPath[i1];

            const vertex_descriptor v0 = target(e0, assemblyGraph);
            const vertex_descriptor v1 = source(e1, assemblyGraph);
            const Shasta2AnchorId anchorId0 = assemblyGraph[v0].anchorId;
            const Shasta2AnchorId anchorId1 = assemblyGraph[v1].anchorId;

            edge_descriptor eBridge;
            tie(eBridge, ignore) = add_edge(
                v0,
                v1,
                Shasta2AssemblyGraphEdge(nextEdgeId++),
                assemblyGraph);
            Shasta2AssemblyGraphEdge& bridgeEdge = assemblyGraph[eBridge];

            if(anchorId0 != anchorId1) {
                ostream html(0);
                const Shasta2TangleMatrix1 tangleMatrix(
                    assemblyGraph,
                    vector<edge_descriptor>(1, e0),
                    vector<edge_descriptor>(1, e1),
                    html);

                Shasta2RestrictedAnchorGraph restrictedAnchorGraph(
                    anchors, journeys, tangleMatrix, 0, 0, html);
                vector<Shasta2RestrictedAnchorGraph::edge_descriptor> longestPath;
                restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

                for(const Shasta2RestrictedAnchorGraph::edge_descriptor re: longestPath) {
                    const auto& rEdge = restrictedAnchorGraph[re];
                    if(rEdge.anchorPair.size() == 0) {
                        bridgeEdge.clear();
                        DINARA_ASSERT(0);
                    }
                    bridgeEdge.emplace_back(rEdge.anchorPair, rEdge.offset);
                    bridgeEdge.anchorChain.push_back(rEdge.anchorPair.anchorIdA);
                }
                if(!longestPath.empty()) {
                    bridgeEdge.anchorChain.push_back(
                        restrictedAnchorGraph[longestPath.back()].anchorPair.anchorIdB);
                }
            }
        }
    }

    // Remove old edges that were part of one or more assembly paths.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        if(pathInitialSegments.contains(e) ||
            pathFinalSegments.contains(e) ||
            pathInternalSegments.contains(e)) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, assemblyGraph);
    }

    // Compress the linear chains we created.
    uint64_t oldCompressDebugLevel = compressDebugLevel;
    if(debug) {
        compressDebugLevel = 1;
    }
    compress();
    if(debug) {
        compressDebugLevel = oldCompressDebugLevel;
    }
}



// The detangling process can generate empty edges (edges without steps).
// This removes them by collapsing the vertices they join.
void Shasta2AssemblyGraph::removeEmptyEdges()
{
    Shasta2AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    // We need to find groups of vertices that need to be collapsed together.
    // Usually it will be just two vertices to be collapsed together,
    // but it could be larger groups if there are adjacent empty edges.
    // So we need to find the connected components generated by the empty edges.

    // Map vertices to integer.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    vector<vertex_descriptor> vertexTable;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
        vertexTable.push_back(v);
    }

    // Initialize the disjoint sets data structure.
    const uint64_t n = vertexIndexMap.size();
    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Loop over the empty edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        if(assemblyGraph[e].empty()) {
            const vertex_descriptor v0 = source(e, assemblyGraph);
            const vertex_descriptor v1 = target(e, assemblyGraph);
            const uint64_t i0 = vertexIndexMap[v0];
            const uint64_t i1 = vertexIndexMap[v1];
            disjointSets.union_set(i0, i1);
            edgesToBeRemoved.push_back(e);
        }
    }

    // Gather the vertices in each connected component.
    vector< vector<uint64_t> > groups(n);
    for(uint64_t i=0; i<n; i++) {
        groups[disjointSets.find_set(i)].push_back(i);
    }

    // Collapse each group into a single vertex.
    for(const vector<uint64_t>& group: groups) {
        if(group.size() < 2) {
            continue;
        }
        if(debug) {
            cout << "Found a group of " << group.size() << " vertices to be collapsed:";
            for(const uint64_t vertexIndex: group) {
                const vertex_descriptor v = vertexTable[vertexIndex];
                cout << " " << assemblyGraph[v].id << "/" << shasta2AnchorIdToString(assemblyGraph[v].anchorId);
            }
            cout << endl;
        }

        // Create the collapsed vertex.
        const Shasta2AnchorId anchorId = assemblyGraph[vertexTable[group.front()]].anchorId;
        const vertex_descriptor vNew =
            add_vertex(Shasta2AssemblyGraphVertex(anchorId, assemblyGraph.nextVertexId++), assemblyGraph);

        // For each vertex in this group, reroute all incoming/outgoing non-empty edges
        // to/from the new vertex.
        for(const uint64_t vertexIndex: group) {
            const vertex_descriptor v = vertexTable[vertexIndex];

            // Loop over non-empty incoming edges.
            BGL_FORALL_INEDGES(v, eOld, assemblyGraph, Shasta2AssemblyGraph) {
                Shasta2AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
                if(oldEdge.empty()) {
                    continue;
                }
                const vertex_descriptor u = source(eOld, assemblyGraph);

                // Create the new edge, with target vNew.
                edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(u, vNew, assemblyGraph);
                Shasta2AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
                newEdge.id = oldEdge.id;
                newEdge.swapSteps(oldEdge);

                edgesToBeRemoved.push_back(eOld);
            }

            // Loop over non-empty outgoing edges.
            BGL_FORALL_OUTEDGES(v, eOld, assemblyGraph, Shasta2AssemblyGraph) {
                Shasta2AssemblyGraphEdge& oldEdge = assemblyGraph[eOld];
                if(oldEdge.empty()) {
                    continue;
                }
                const vertex_descriptor u = target(eOld, assemblyGraph);

                // Create the new edge, with source vNew.
                edge_descriptor eNew;
                tie(eNew, ignore) = add_edge(vNew, u, assemblyGraph);
                Shasta2AssemblyGraphEdge& newEdge = assemblyGraph[eNew];
                newEdge.id = oldEdge.id;
                newEdge.swapSteps(oldEdge);

                edgesToBeRemoved.push_back(eOld);
            }
        }
    }

    deduplicate(edgesToBeRemoved);
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, assemblyGraph);
    }

    // Remove any vertices that were left isolated.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        if((in_degree(v, assemblyGraph) == 0) and (out_degree(v, assemblyGraph) == 0)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::remove_vertex(v, assemblyGraph);
    }
}

void Shasta2AssemblyGraph::removeIsolatedVertices()
{
    Shasta2AssemblyGraph& assemblyGraph = *this;
    vector<vertex_descriptor> toBeRemoved;

    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        if((in_degree(v, assemblyGraph) == 0) && (out_degree(v, assemblyGraph) == 0)) {
            toBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: toBeRemoved) {
        boost::remove_vertex(v, assemblyGraph);
    }
}

void Shasta2AssemblyGraph::removeLowN50Components(uint64_t minN50)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;
    const bool debug = false;

    uint64_t vertexIndex = 0;
    vector<vertex_descriptor> vertexTable;
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        vertexTable.push_back(v);
        vertexIndexMap.insert({v, vertexIndex++});
    }

    Shasta2DisjointSets disjointSets(vertexIndexMap.size());
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        disjointSets.unionSet(vertexIndexMap[v0], vertexIndexMap[v1]);
    }

    vector< vector<uint64_t> > components;
    disjointSets.gatherComponents(1, components);
    if(debug) {
        cout << "Found " << components.size() << " connected components of the Shasta2AssemblyGraph." << endl;
    }

    for(const vector<uint64_t>& component: components) {
        vector<uint64_t> lengths;
        for(const uint64_t vIndex0: component) {
            const vertex_descriptor v0 = vertexTable[vIndex0];
            BGL_FORALL_OUTEDGES(v0, e, assemblyGraph, Shasta2AssemblyGraph) {
                const vertex_descriptor v1 = target(e, assemblyGraph);
                const uint64_t vIndex1 = vertexIndexMap[v1];
                DINARA_ASSERT(disjointSets.findSet(vIndex0) == disjointSets.findSet(vIndex1));
                lengths.push_back(assemblyGraph[e].length());
            }
        }

        const uint64_t totalLength = std::accumulate(lengths.begin(), lengths.end(), uint64_t(0));
        std::ranges::sort(lengths, std::greater<uint64_t>());
        uint64_t cumulativeLength = 0;
        uint64_t n50 = 0;
        for(const uint64_t length: lengths) {
            cumulativeLength += length;
            if(2 * cumulativeLength >= totalLength) {
                n50 = length;
                break;
            }
        }

        const bool keep = (n50 >= minN50);
        if(debug) {
            cout << (keep ? "Keeping" : "Discarding") <<
                " a connected component with " << lengths.size() <<
                " segments, total length " << totalLength <<
                ", N50 " << n50 << endl;
        }

        if(not keep) {
            for(const uint64_t vIndex: component) {
                const vertex_descriptor v = vertexTable[vIndex];
                boost::clear_vertex(v, assemblyGraph);
                boost::remove_vertex(v, assemblyGraph);
            }
        }
    }
}

void Shasta2AssemblyGraph::assembleAll()
{
    writePerformanceStatistics("Shasta2AssemblyGraph::assembleAll begins");

    Shasta2AssemblyGraph& assemblyGraph = *this;
    edgesToBeAssembled.clear();
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        edgesToBeAssembled.push_back(e);
    }
    assemble();
    edgesToBeAssembled.clear();

    writePerformanceStatistics("Shasta2AssemblyGraph::assembleAll ends");
}

void Shasta2AssemblyGraph::assemble(edge_descriptor e)
{
    edgesToBeAssembled.clear();
    edgesToBeAssembled.push_back(e);
    assemble();
}

void Shasta2AssemblyGraph::assemble()
{
    performanceLog << timestamp << "Sequence assembly begins for " << edgesToBeAssembled.size() <<
        " assembly graph edges." << endl;
    Shasta2AssemblyGraph& assemblyGraph = *this;

    stepsToBeAssembled.clear();
    for(const edge_descriptor e: edgesToBeAssembled) {
        Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
        for(uint64_t i=0; i<edge.size(); i++) {
            stepsToBeAssembled.push_back(make_pair(e, i));
        }
    }

    const uint64_t batchCount = 1;
    setupLoadBalancing(stepsToBeAssembled.size(), batchCount);
    runThreads(&Shasta2AssemblyGraph::assembleThreadFunction, options.threadCount);

    for(const edge_descriptor e: edgesToBeAssembled) {
        assemblyGraph[e].wasAssembled = true;
    }

    edgesToBeAssembled.clear();
    stepsToBeAssembled.clear();

    performanceLog << timestamp << "Sequence assembly ends." << endl;
}

void Shasta2AssemblyGraph::assembleThreadFunction(uint64_t /* threadId */)
{
    Shasta2AssemblyGraph& assemblyGraph = *this;

    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        for(uint64_t j=begin; j!=end; j++) {
            const auto& p = stepsToBeAssembled[j];
            const edge_descriptor e = p.first;
            const uint64_t i = p.second;
            Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
            DINARA_ASSERT(i < edge.size());
            assembleStep(e, i);
        }
    }
}

void Shasta2AssemblyGraph::clearSequence()
{
    Shasta2AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.wasAssembled = false;
        for(Shasta2AssemblyGraphEdgeStep& step: edge) {
            step.sequence.clear();
            step.sequence.shrink_to_fit();
        }
    }
}

void Shasta2AssemblyGraph::assembleStep(edge_descriptor e, uint64_t i)
{
    DINARA_ASSERT(anchorsPointer);

    Shasta2AssemblyGraph& assemblyGraph = *this;
    Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
    DINARA_ASSERT(i < edge.size());

    Shasta2AssemblyGraphEdgeStep& step = edge[i];
    step.sequence.clear();
    if(step.anchorPair.anchorIdA == step.anchorPair.anchorIdB) {
        return;
    }

    vector<OrientedReadId> additionalOrientedReadIds;
    if(i > 0) {
        copy(
            edge[i - 1].anchorPair.orientedReadIds.begin(),
            edge[i - 1].anchorPair.orientedReadIds.end(),
            back_inserter(additionalOrientedReadIds));
    }
    if(i + 1 < edge.size()) {
        copy(
            edge[i + 1].anchorPair.orientedReadIds.begin(),
            edge[i + 1].anchorPair.orientedReadIds.end(),
            back_inserter(additionalOrientedReadIds));
    }
    deduplicate(additionalOrientedReadIds);

    ostream html(0);
    Shasta2LocalAssembly4 localAssembly(
        *anchorsPointer,
        options.abpoaMaxLength,
        html,
        false,
        edge[i].anchorPair,
        additionalOrientedReadIds);
    step.sequence = localAssembly.sequence;
}

vector<dinara::Base> Shasta2AssemblyGraph::consensusSequence(
    const vector< vector<dinara::Base> >& sequences)
{
    if(sequences.empty()) {
        return {};
    }

    // Prefer the most frequent exact sequence. Break ties by longer sequence.
    map<vector<dinara::Base>, uint64_t> sequenceCounts;
    for(const auto& sequence: sequences) {
        ++sequenceCounts[sequence];
    }

    uint64_t bestCount = 0;
    const vector<dinara::Base>* best = nullptr;
    for(const auto& p: sequenceCounts) {
        const vector<dinara::Base>& sequence = p.first;
        const uint64_t count = p.second;
        if(best == nullptr || count > bestCount ||
            (count == bestCount && sequence.size() > best->size())) {
            best = &sequence;
            bestCount = count;
        }
    }

    if(best == nullptr) {
        return {};
    }
    return *best;
}

void Shasta2AssemblyGraph::writePerformanceStatistics(const string& message) const
{
    std::lock_guard<std::mutex> lock(mutex);
    performanceLog << timestamp << message << ". " <<
        num_vertices(*this) << " vertices, " <<
        num_edges(*this) << " edges." << endl;
}

void Shasta2AssemblyGraph::write(const string& stage)
{
    cout << "Stage " << stage << ": " <<
        num_vertices(*this) << " vertices, " <<
        num_edges(*this) << " edges. Next edge id is " << nextEdgeId << "." << endl;

    saveStage(stage);
    writeGfa("Shasta2Assembly-" + stage + ".gfa");
    writeGraphviz("Shasta2Assembly-" + stage + ".dot");
    writeCsv("Shasta2Assembly-" + stage + ".csv");
}

void Shasta2AssemblyGraph::writeIntermediateStageIfRequested(const string& stage)
{
    if(!options.writeIntermediateAssemblyStages) {
        return;
    }

    write(stage);
}

void Shasta2AssemblyGraph::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    if(!gfa) {
        throw runtime_error("Cannot open " + fileName + " for writing.");
    }

    writeGfa(gfa);
}

void Shasta2AssemblyGraph::writeGfa(ostream& gfa) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;

    gfa << "H\tVN:Z:1.0\n";

    vector<dinara::Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
        const double coverage = Shasta2AssemblyGraphEdgeLengthWeightedAverageCoverage(edge);

        gfa << "S\t";
        gfa << edge.id << "\t";

        if(edge.wasAssembled) {
            edge.getSequence(sequence);
            copy(sequence.begin(), sequence.end(), ostream_iterator<dinara::Base>(gfa));
            const uint64_t length = sequence.size();
            gfa << "\tLN:i:" << length;
            gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(length)));
        } else {
            if(edge.empty()) {
                gfa << "*\tLN:i:0";
            } else {
                const uint64_t offset = edge.offset();
                gfa << "*\tLN:i:" << offset;
                gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(offset)));
            }
        }

        // Add window sequence tag if available.
        if(!edge.windowSequence.empty()) {
            gfa << "\tws:Z:";
            for(uint64_t i = 0; i < edge.windowSequence.size(); i++) {
                if(i > 0) gfa << ",";
                gfa << edge.windowSequence[i];
            }
        }

        // Add anchor chain length tag.
        if(!edge.anchorChain.empty()) {
            gfa << "\tac:i:" << edge.anchorChain.size();
        }

        gfa << "\n";
    }

    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        BGL_FORALL_INEDGES(v, e0, assemblyGraph, Shasta2AssemblyGraph) {
            const uint64_t id0 = assemblyGraph[e0].id;
            BGL_FORALL_OUTEDGES(v, e1, assemblyGraph, Shasta2AssemblyGraph) {
                const uint64_t id1 = assemblyGraph[e1].id;
                gfa << "L\t" << id0 << "\t+\t" << id1 << "\t+\t*\n";
            }
        }
    }
}

void Shasta2AssemblyGraph::writeFasta(const string& stage) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;
    ofstream fasta("Shasta2Assembly-" + stage + ".fasta");

    vector<dinara::Base> sequence;
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
        edge.getSequence(sequence);
        fasta << ">" << edge.id << "\n";
        copy(sequence.begin(), sequence.end(), ostream_iterator<dinara::Base>(fasta));
        fasta << "\n";
    }
}

void Shasta2AssemblyGraph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}

void Shasta2AssemblyGraph::writeGraphviz(ostream& dot) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;

    dot << "digraph Shasta2AssemblyGraph {\n";

    BGL_FORALL_VERTICES(v, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphVertex& vertex = assemblyGraph[v];
        dot <<
            vertex.id <<
            " [label=\"" << shasta2AnchorIdToString(vertex.anchorId);
        if(vertex.windowId != Shasta2AssemblyGraphVertex::noWindow) {
            dot << "\\nW" << vertex.windowId;
        }
        dot << "\\n" << vertex.id << "\"];\n";
    }

    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);
        const Shasta2AssemblyGraphVertex& vertex0 = assemblyGraph[v0];
        const Shasta2AssemblyGraphVertex& vertex1 = assemblyGraph[v1];
        dot <<
            vertex0.id << "->" <<
            vertex1.id <<
            " [label=\"" << edge.id << "\\n" <<
            (edge.wasAssembled ? edge.sequenceLength() : edge.offset()) <<
            "\\n" << edge.size() << " steps";
        if(!edge.windowSequence.empty()) {
            dot << "\\nW:";
            for(uint64_t i = 0; i < edge.windowSequence.size(); i++) {
                if(i > 0) dot << ",";
                dot << edge.windowSequence[i];
            }
        }
        dot << "\"];\n";
    }

    dot << "}\n";
}

void Shasta2AssemblyGraph::writeCsv(const string& fileName) const
{
    ofstream csv(fileName);
    writeCsv(csv);
}

void Shasta2AssemblyGraph::writeCsv(ostream& csv) const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;
    csv << "Segment,Number of steps,Anchor chain length,Average coverage,Estimated length,Actual length,Window sequence\n";
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t coverage = uint64_t(std::round(Shasta2AssemblyGraphEdgeAverageCoverage(edge)));
        csv <<
            edge.id << "," <<
            edge.size() << "," <<
            edge.anchorChain.size() << "," <<
            coverage << "," <<
            edge.offset() << ",";
        if(edge.wasAssembled) {
            csv << edge.sequenceLength();
        }
        csv << ",";
        for(uint64_t i = 0; i < edge.windowSequence.size(); i++) {
            if(i > 0) csv << " ";
            csv << edge.windowSequence[i];
        }
        csv << "\n";
    }
}

void Shasta2AssemblyGraph::check() const
{
    const Shasta2AssemblyGraph& assemblyGraph = *this;

    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
        DINARA_ASSERT(!edge.empty());

        const vertex_descriptor v0 = source(e, assemblyGraph);
        const vertex_descriptor v1 = target(e, assemblyGraph);

        const Shasta2AnchorId anchorId0 = assemblyGraph[v0].anchorId;
        const Shasta2AnchorId anchorId1 = assemblyGraph[v1].anchorId;

        DINARA_ASSERT(edge.front().anchorPair.anchorIdA == anchorId0);
        DINARA_ASSERT(edge.back().anchorPair.anchorIdB == anchorId1);

        for(uint64_t i1=1; i1<edge.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            DINARA_ASSERT(edge[i0].anchorPair.anchorIdB == edge[i1].anchorPair.anchorIdA);
        }

        // Validate anchor chain consistency with steps, if populated.
        if(!edge.anchorChain.empty()) {
            DINARA_ASSERT(edge.anchorChain.size() == edge.size() + 1);
            DINARA_ASSERT(edge.anchorChain.front() == anchorId0);
            DINARA_ASSERT(edge.anchorChain.back() == anchorId1);
            for(uint64_t i = 0; i < edge.size(); i++) {
                DINARA_ASSERT(edge.anchorChain[i] == edge[i].anchorPair.anchorIdA);
                DINARA_ASSERT(edge.anchorChain[i + 1] == edge[i].anchorPair.anchorIdB);
            }
        }
    }
}

void Shasta2AssemblyGraph::save(ostream& s) const
{
    boost::archive::binary_oarchive archive(s);
    archive << *this;
}

void Shasta2AssemblyGraph::load(istream& s)
{
    boost::archive::binary_iarchive archive(s);
    archive >> *this;
}

void Shasta2AssemblyGraph::save(const string& name) const
{
    if(largeDataFileNamePrefix.empty()) {
        return;
    }

    std::ostringstream s;
    save(s);
    const string dataString = s.str();

    MemoryMapped::Vector<char> data;
    data.createNew(largeDataName(name), largeDataPageSize);
    data.resize(dataString.size());
    copy(dataString.begin(), dataString.end(), data.begin());
}

void Shasta2AssemblyGraph::load(const string& name)
{
    MemoryMapped::Vector<char> data;
    try {
        data.accessExistingReadOnly(largeDataName(name));
    } catch(std::exception&) {
        throw runtime_error(name + " is not available.");
    }

    const string dataString(data.begin(), data.size());
    std::istringstream s(dataString);
    try {
        load(s);
    } catch(std::exception& e) {
        throw runtime_error(string("Error reading ") + name + ": " + e.what());
    }
}

void Shasta2AssemblyGraph::saveStage(const string& stage) const
{
    save("Shasta2AssemblyGraph-" + stage);
}

void Shasta2AssemblyGraph::loadStage(const string& stage)
{
    load("Shasta2AssemblyGraph-" + stage);
}
