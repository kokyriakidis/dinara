#pragma once

#include "Shasta2AnchorGraph.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "Base.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <limits>
#include <mutex>

#include "memory.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace dinara {

        class Shasta2AssemblyGraph;
        class Shasta2AssemblyGraphVertex;
        class Shasta2AssemblyGraphEdge;
        class Shasta2AssemblyGraphEdgeStep;
        class Shasta2AssemblyGraphOptions;
        class Shasta2Superbubble;
        class Shasta2SuperbubbleChain;

        using Shasta2AssemblyGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            Shasta2AssemblyGraphVertex,
            Shasta2AssemblyGraphEdge>;
}

class dinara::Shasta2AssemblyGraphVertex {
public:
    Shasta2AnchorId anchorId = invalid<Shasta2AnchorId>;
    uint64_t id = invalid<uint64_t>;

    // Window that this anchor belongs to (noWindow if unmapped).
    static constexpr uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    uint32_t windowId = noWindow;

    Shasta2AssemblyGraphVertex(Shasta2AnchorId anchorId, uint64_t id) :
        anchorId(anchorId),
        id(id)
    {}
    Shasta2AssemblyGraphVertex() {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorId;
        ar & id;
        ar & windowId;
    }
};

class dinara::Shasta2AssemblyGraphEdgeStep {
public:
    Shasta2AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    vector<Base> sequence;

    Shasta2AssemblyGraphEdgeStep(
        const Shasta2AnchorPair& anchorPair,
        uint64_t offset) :
        anchorPair(anchorPair),
        offset(offset)
    {}
    Shasta2AssemblyGraphEdgeStep() {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorPair;
        ar & offset;
        ar & sequence;
    }
};

class dinara::Shasta2AssemblyGraphEdge : public vector<Shasta2AssemblyGraphEdgeStep> {
public:
    uint64_t id = invalid<uint64_t>;
    bool wasAssembled = false;

    // Full chain of AnchorIds along this edge, like shasta2's Chain.
    // Includes the source and target vertex anchors.
    // anchorChain[0] == source vertex anchorId,
    // anchorChain.back() == target vertex anchorId.
    // For a single-step edge (one anchor pair A->B), anchorChain = {A, B}.
    // For a multi-step edge, anchorChain = {A0, A1, ..., An} where
    // step[i] connects anchorChain[i] to anchorChain[i+1].
    vector<Shasta2AnchorId> anchorChain;

    // Sequence of normalized window IDs traversed by this edge's anchor chain.
    // Consecutive duplicates are removed. noWindow entries are omitted.
    // This gives the window-level path of the edge.
    static constexpr uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    vector<uint32_t> windowSequence;

    Shasta2AssemblyGraphEdge(uint64_t id = invalid<uint64_t>) :
        id(id)
    {}

    uint64_t offset() const;
    uint64_t sequenceLength() const;
    void getSequence(vector<Base>&) const;
    uint64_t length() const;
    double averageCoverage() const;
    double lengthWeightedAverageCoverage() const;
    void swapSteps(Shasta2AssemblyGraphEdge& that)
    {
        vector<Shasta2AssemblyGraphEdgeStep>::swap(that);
    }

    Shasta2AnchorId firstAnchorId() const
    {
        return front().anchorPair.anchorIdA;
    }
    Shasta2AnchorId lastAnchorId() const
    {
        return back().anchorPair.anchorIdB;
    }

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<Shasta2AssemblyGraphEdgeStep> >(*this);
        ar & id;
        ar & wasAssembled;
        ar & anchorChain;
        ar & windowSequence;
    }
};

class dinara::Shasta2AssemblyGraphOptions {
public:
    uint64_t simplifyMaxIterationCount = 3;
    uint64_t threadCount = 1;
    bool writeIntermediateAssemblyStages = true;
    uint64_t bubbleCleanupMaxBubbleLength = 10000;
    uint64_t bubbleCleanupMinCommonCount = 6;
    uint64_t phasingDistance = 12;
    uint64_t phasingMinCommonCount = 6;
    double detangleEpsilon = 0.05;
    double detangleMaxLogP = 30.;
    double detangleMinLogPDelta = 10.;
    uint64_t findSuperbubblesMaxDistance = 10;
    uint64_t representativeRegionStepCount = 10;
    uint64_t readFollowingMinCommonCount = 6;
    double readFollowingMinCorrectedJaccard = 0.7;
    uint32_t readFollowingPruneLength = 100000;
    uint64_t readFollowingSegmentLengthThreshold = 500000;
    uint64_t readFollowingMinPathLength = 2;
    uint64_t readFollowingMaxPathCount = 0; // 0 means no limit.
    uint64_t abpoaMaxLength = 1000;
};

class dinara::Shasta2AssemblyGraph :
    public Shasta2AssemblyGraphBaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<Shasta2AssemblyGraph> {
public:
    Shasta2AssemblyGraph(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        const Shasta2AnchorGraph&,
        const Shasta2AssemblyGraphOptions& = Shasta2AssemblyGraphOptions());

    // Construct from anchor graph with window information.
    // This stores the anchor-to-window mapping and populates
    // anchorChain, windowSequence, and vertex windowId fields.
    Shasta2AssemblyGraph(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        const Shasta2AnchorGraph&,
        const vector<AnchorWindow>& anchorWindows,
        const Shasta2AssemblyGraphOptions& = Shasta2AssemblyGraphOptions());

    Shasta2AssemblyGraph(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        const Shasta2AssemblyGraphOptions&,
        const string& stage);

    Shasta2AssemblyGraph(
        const MappedMemoryOwner&,
        const string& name);

    uint64_t nextVertexId = 0;
    uint64_t nextEdgeId = 0;

    void simplifyAndAssemble();
    uint64_t compress();
    uint64_t compressDebugLevel = 0; // 1=minimal, 2=compact, 3=detailed.
    uint64_t bubbleCleanup();
    uint64_t phaseSuperbubbleChains();

    // Graph cleaning steps, to be called after the initial compress()
    // and before simplifyAndAssemble().
    // Adapted from MBG (Multiplex de Bruijn Graph assembler).
    uint64_t removeLowCoverageTips(
        double maxRemovableCoverage = 3.,
        double minSafeCoverage = 10.,
        uint64_t maxRemovableLength = 10000);
    uint64_t removeLowCoverageCrosslinks(
        double maxRemovableCoverage = 2.,
        double minSafeCoverage = 10.);
    uint64_t cleanByCopyNumber(
        double estimatedAverageCoverage = 0.); // 0 means auto-estimate

    // Remove short tip chains (like hifiasm's asg_arc_cut_tips).
    // A tip is removed if its total window count is <= maxTipWindows
    // AND its total length is <= maxTipLength bp.
    // Shorter tips are processed first so their removal can expose longer ones.
    uint64_t removeShortTips(uint32_t maxTipWindows = 3, uint64_t maxTipLength = 90000);

    // Pop superbubbles by removing low-coverage alternative paths.
    // Follows Verkko's approach (pop_bubbles_coverage_based.py).
    uint64_t popSuperbubbles(
        uint64_t maxBubbleSize = 10,
        double maxPoppableCoverageFraction = 0.5);

    // Remove parallel edges between the same pair of vertices.
    // When multiple edges connect the same source to the same target,
    // keep the one with highest average coverage and remove the rest.
    uint64_t removeParallelEdges();

    // Remove weak stalks (dead-end linear chains leading to branch points).
    // Ported from Shasta2AnchorGraph::cutWeakStalksLeadingToBranch.
    // A stalk starts at a tip vertex (in-degree 0 or out-degree 0),
    // follows a linear chain, and if it reaches a branch point
    // with the read union across all edges still <= maxReadCount, it is cut.
    uint64_t cutWeakStalks(uint64_t maxReadCount);
    void phaseSuperbubbleChainsThreadFunction(uint64_t threadId);

    // Convert the assembly graph to a bidirected anchor graph.
    BidirectedAnchorGraph toBidirected(
        const std::vector<AnchorWindow>& anchorWindows,
        const Shasta2Journeys& journeys) const;

    // Compute compressed journeys in the Shasta2AssemblyGraph.
    void computeJourneys();

    vector< vector<edge_descriptor> > compressedJourneys;



    void findAndConnectAssemblyPaths();
    bool canConnect(edge_descriptor, edge_descriptor) const;
    void removeEmptyEdges();
    void removeIsolatedVertices();
    void removeLowN50Components(uint64_t minN50);
    void assembleAll();
    void writeIntermediateStageIfRequested(const string& stage);

    // Export the assembly graph as a shasta2-compatible AnchorGraph.
    // Each assembly graph edge step becomes one anchor graph edge.
    void saveForShasta2(const string& fileName) const;

    bool hasSelfEdge(vertex_descriptor v) const
    {
        bool edgeExists = false;
        tie(ignore, edgeExists) = boost::edge(v, v, *this);
        return edgeExists;
    }

    class OrderById {
    public:
        explicit OrderById(const Shasta2AssemblyGraph& assemblyGraph) :
            assemblyGraph(assemblyGraph)
        {}
        const Shasta2AssemblyGraph& assemblyGraph;
        bool operator()(vertex_descriptor x, vertex_descriptor y) const
        {
            return assemblyGraph[x].id < assemblyGraph[y].id;
        }
        bool operator()(edge_descriptor x, edge_descriptor y) const
        {
            return assemblyGraph[x].id < assemblyGraph[y].id;
        }
        using EdgePair = pair<edge_descriptor, edge_descriptor>;
        bool operator()(const EdgePair& x, const EdgePair& y) const
        {
            if(assemblyGraph[x.first].id < assemblyGraph[y.first].id) {
                return true;
            }
            if(assemblyGraph[x.first].id > assemblyGraph[y.first].id) {
                return false;
            }
            return assemblyGraph[x.second].id < assemblyGraph[y.second].id;
        }
    };
    const OrderById orderById;

    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
    void writeFasta(const string& stage) const;
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
    void writeCsv(const string& fileName) const;
    void writeCsv(ostream&) const;

    // Diagnostic (read-only): report candidate bridges between maximal 1-1
    // linear segments, using read window-paths. A segment is an edge of this
    // graph; its endpoint windows are the head (windowSequence.front()) and
    // tail (windowSequence.back()). For each read, walk its window path
    // (readWindows, raw window IDs) and record every transition from one
    // segment's tail window to another segment's head window, accumulating
    // distinct physical-read support. Reports how many segment pairs reads
    // bridge and the per-tail target multiplicity (1 target vs many). Adds no
    // edges. windowReads/readWindows come from the source Shasta2AnchorGraph.
    void reportSegmentBridges(
        const std::map<uint32_t, std::set<uint32_t>>& windowReads,
        const std::map<uint32_t, vector<uint32_t>>& readWindows) const;

    void write(const string& stage);
    void check() const;
    void clearSequence();

    void findStrongComponents(vector< vector<vertex_descriptor> >&) const;
    void colorStrongComponents() const;

    // Window-level information.
    // Populated by the constructor that takes anchorWindows.

    // Per-anchor window assignment (anchorId -> normalized windowId).
    // noWindow means unmapped.
    static constexpr uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    vector<uint32_t> anchorToWindow;
    uint32_t windowCount = 0;

    // Pointer to the anchor windows (not owned, must outlive this object).
    const vector<AnchorWindow>* anchorWindowsPointer = nullptr;

    // Populate anchorChain and windowSequence for all edges,
    // and windowId for all vertices, using the current anchorToWindow mapping.
    // Called automatically by the window-aware constructor.
    // Can also be called after graph modifications (compress, etc.)
    // to refresh the window annotations.
    void populateWindowAnnotations();

    // Get the normalized window ID for an anchor.
    uint32_t getWindowId(Shasta2AnchorId anchorId) const
    {
        const uint64_t aid = uint64_t(anchorId);
        if(aid >= anchorToWindow.size()) return noWindow;
        const uint32_t w = anchorToWindow[aid];
        if(w == noWindow) return noWindow;
        return (w >= windowCount) ? (w - windowCount) : w;
    }

    // Check if window information is available.
    bool hasWindowInfo() const { return windowCount > 0; }

    class PhaseSuperbubbleChainsData {
    public:
        shared_ptr< vector<Shasta2SuperbubbleChain> > superbubbleChains;
        uint64_t totalChangeCount = 0;
    };
    PhaseSuperbubbleChainsData phaseSuperbubbleChainsData;

    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<Shasta2AssemblyGraphBaseClass>(*this);
        ar & nextVertexId;
        ar & nextEdgeId;
        ar & anchorToWindow;
        ar & windowCount;
    }
    void save(ostream&) const;
    void load(istream&);

    void save(const string& name = "Shasta2AssemblyGraph") const;
    void load(const string& name);
    void saveStage(const string& stage) const;
    void loadStage(const string& stage);
    const Shasta2AssemblyGraphOptions& getOptions() const
    {
        return options;
    }
    const Shasta2Anchors* getAnchorsPointer() const
    {
        return anchorsPointer;
    }
    const Shasta2Journeys* getJourneysPointer() const
    {
        return journeysPointer;
    }
    std::mutex& getMutex()
    {
        return mutex;
    }
    mutable std::mutex mutex;

private:
    class Bubble {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        vector<edge_descriptor> edges;
    };
    uint64_t bubbleCleanupIteration(vector< pair<vertex_descriptor, vertex_descriptor> >& excludeList);
    void findBubbles(vector<Bubble>&) const;
    bool bubbleCleanup(const Bubble&);
    bool analyzeBubble(
        const Bubble&,
        const vector<uint64_t> minRepeatCount,
        vector< pair<uint64_t, uint64_t> >& similarPairs) const;
    void findSuperbubbles(vector<Shasta2Superbubble>&) const;
    void removeContainedSuperbubbles(vector<Shasta2Superbubble>&) const;
    void findSuperbubbleChains(
        const vector<Shasta2Superbubble>&,
        vector<Shasta2SuperbubbleChain>&) const;
    void findAssemblyPaths(vector< vector<edge_descriptor> >&) const;
    void connectAssemblyPaths(const vector< vector<edge_descriptor> >&);
    void writeSuperbubbles(const vector<Shasta2Superbubble>&, const string& fileName) const;
    void writeSuperbubblesForBandage(const vector<Shasta2Superbubble>&, const string& fileName) const;
    void writeSuperbubbleChains(const vector<Shasta2SuperbubbleChain>&, const string& fileName) const;
    void writeSuperbubbleChainsForBandage(const vector<Shasta2SuperbubbleChain>&, const string& fileName) const;
    void writePerformanceStatistics(const string& message) const;
    void assemble(edge_descriptor);
    void assemble();
    void assembleThreadFunction(uint64_t threadId);
    void assembleStep(edge_descriptor, uint64_t);
    static vector<dinara::Base> consensusSequence(
        const vector< vector<dinara::Base> >& sequences);

    [[maybe_unused]] const Shasta2Anchors* anchorsPointer = 0;
    [[maybe_unused]] const Shasta2Journeys* journeysPointer = 0;
    Shasta2AssemblyGraphOptions options;
    vector<edge_descriptor> edgesToBeAssembled;
    vector< pair<edge_descriptor, uint64_t> > stepsToBeAssembled;
};
