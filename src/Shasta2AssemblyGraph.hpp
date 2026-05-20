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

    Shasta2AssemblyGraphVertex(Shasta2AnchorId anchorId, uint64_t id) :
        anchorId(anchorId),
        id(id)
    {}
    Shasta2AssemblyGraphVertex() {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorId;
        ar & id;
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

    // Pop superbubbles by removing low-coverage alternative paths.
    // Follows Verkko's approach (pop_bubbles_coverage_based.py).
    uint64_t popSuperbubbles(
        uint64_t maxBubbleSize = 10,
        double maxPoppableCoverageFraction = 0.5);
    void phaseSuperbubbleChainsThreadFunction(uint64_t threadId);

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
    void write(const string& stage);
    void check() const;
    void clearSequence();

    void findStrongComponents(vector< vector<vertex_descriptor> >&) const;
    void colorStrongComponents() const;

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
