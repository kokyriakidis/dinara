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

        for(const uint64_t i: branchGroup) {
            boost::remove_edge(bubble.edges[i], assemblyGraph);
        }
    }

    DINARA_ASSERT(out_degree(bubble.v0, assemblyGraph) > 0);
    DINARA_ASSERT(in_degree(bubble.v1, assemblyGraph) > 0);

    return true;
}

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

    const auto hasSelfEdge = [&assemblyGraph](const vertex_descriptor v) {
        bool edgeExists = false;
        tie(ignore, edgeExists) = boost::edge(v, v, assemblyGraph);
        return edgeExists;
    };

    vector< pair<vertex_descriptor, vertex_descriptor> > forwardPairs;
    BGL_FORALL_VERTICES(vSource, assemblyGraph, Shasta2AssemblyGraph) {
        if(hasSelfEdge(vSource)) {
            continue;
        }
        const vertex_descriptor vTarget =
            findConvergingVertexGeneral(assemblyGraph, vSource, options.findSuperbubblesMaxDistance);
        if(vTarget != null_vertex() && !hasSelfEdge(vTarget)) {
            forwardPairs.emplace_back(vSource, vTarget);
        }
    }
    sort(forwardPairs.begin(), forwardPairs.end());

    const boost::reverse_graph<Shasta2AssemblyGraph> reverseAssemblyGraph(assemblyGraph);
    vector< pair<vertex_descriptor, vertex_descriptor> > backwardPairs;
    BGL_FORALL_VERTICES(vTarget, assemblyGraph, Shasta2AssemblyGraph) {
        if(hasSelfEdge(vTarget)) {
            continue;
        }
        const vertex_descriptor vSource =
            findConvergingVertexGeneral(reverseAssemblyGraph, vTarget, options.findSuperbubblesMaxDistance);
        if(vSource != null_vertex() && !hasSelfEdge(vSource)) {
            backwardPairs.emplace_back(vSource, vTarget);
        }
    }
    sort(backwardPairs.begin(), backwardPairs.end());

    vector< pair<vertex_descriptor, vertex_descriptor> > bidirectionalPairs;
    set_intersection(
        forwardPairs.begin(),
        forwardPairs.end(),
        backwardPairs.begin(),
        backwardPairs.end(),
        back_inserter(bidirectionalPairs));

    superbubbles.clear();
    for(const auto& p: bidirectionalPairs) {
        superbubbles.emplace_back(assemblyGraph, p.first, p.second);
    }
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

        // Concatenate the steps of all the edges in the chain.
        for(const edge_descriptor e: chain) {
            const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
            copy(edge.begin(), edge.end(), back_inserter(edgeNew));
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
            gfa << "\n";

        } else {
            if(edge.empty()) {
                gfa << "*\tLN:i:0\n";
            } else {
                const uint64_t offset = edge.offset();
                gfa << "*\tLN:i:" << offset;
                gfa << "\tRC:i:" << uint64_t(std::round(coverage * double(offset)));
                gfa << "\n";
            }
        }
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
            " [label=\"" << shasta2AnchorIdToString(vertex.anchorId) <<
            "\\n" << vertex.id << "\"];\n";
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
            "\\n" << edge.size() <<
            "\"];\n";
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
    csv << "Segment,Number of steps,Average coverage,Estimated length,Actual length\n";
    BGL_FORALL_EDGES(e, assemblyGraph, Shasta2AssemblyGraph) {
        const Shasta2AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t coverage = uint64_t(std::round(Shasta2AssemblyGraphEdgeAverageCoverage(edge)));
        csv <<
            edge.id << "," <<
            edge.size() << "," <<
            coverage << "," <<
            edge.offset() << ",";
        if(edge.wasAssembled) {
            csv << edge.sequenceLength();
        }
        csv << ",\n";
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
