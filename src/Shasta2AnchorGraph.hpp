#pragma once

// Shasta2AnchorGraph.hpp


#include "Reads.hpp"
#include "Shasta2AnchorPair.hpp"
#include "Shasta2Anchors.hpp"
#include "Shasta2Journeys.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/base_object.hpp>

// Standard library.
#include <map>
#include <set>
#include "utility.hpp"
#include "vector.hpp"

namespace dinara {

        class Shasta2AnchorGraph;
        class Shasta2AnchorGraphEdge;
        class Shasta2Anchors;

        // Result of the POA-topology detection pass (transcribeHetBubbles):
        // per-edge abPOA is used to find real sites (SNPs and indels) on each
        // edge, each transcribed into new het anchors (one per passing
        // non-deletion allele). See Shasta2AnchorGraphHetOnGraph.cpp's file
        // header for why this only appends anchors and never touches the graph.
        struct HetOnGraphResult {
            uint64_t edgesTotal = 0;        // all edges in the graph
            uint64_t edgesConsidered = 0;   // coverage >= minCommonForHet, MSA attempted
            uint64_t edgesMsad = 0;         // MSA actually produced (>=2 non-empty rows)
            uint64_t edgesSkippedMirror = 0; // this edge's RC mirror is processed instead
            uint64_t edgesSkippedCoverage = 0;
            uint64_t edgesSkippedLen = 0;   // skipped by maxLen guard
            uint64_t edgesSkippedIdentical = 0; // all read sequences identical (no MSA)
            uint64_t edgesPlanned = 0;      // edges with >=1 real site found
            uint64_t edgesPlannedMultiSite = 0; // of edgesPlanned, those with >1 real site
            uint64_t edgesDeferredEndBubble = 0; // the real-site chain would touch a span end
            uint64_t edgesDeferredComplex = 0;   // other unsupported shapes (currently unused)
            uint64_t sitesTranscribed = 0;  // total real sites found across all edges
            uint64_t hetAnchorsCreated = 0; // new allele-arm anchors appended
            double elapsedSeconds = 0.0;
        };

        // Detect abPOA local topology (SNPs and indels) on anchor-graph edges
        // and append a new het anchor for each passing non-deletion allele at
        // each real site (Shasta2Anchors::appendHetAnchorPair). For every edge
        // whose two-sided coverage is at least minCommonForHet, an abPOA MSA is
        // run over the reads' inter-anchor sequences; where reads diverge into
        // >=2 alleles clearing a one-sided binomial significance test against
        // hetErrorRate (myloasm-style; see the .cpp file's binomialTailPValue),
        // each non-deletion allele becomes a new anchor. This mutates the
        // anchor store but NEVER the graph passed in -- the caller must
        // rebuild journeys (Shasta2Journeys::rebuildAfterNewAnchors) and then a
        // fresh Shasta2AnchorGraph from them before the new anchors take
        // effect; see the .cpp file header for why. Detection is parallelized
        // over threadCount (0 = hardware concurrency); anchor creation is
        // serial. Returns counts for reporting.
        HetOnGraphResult transcribeHetBubbles(
            const Shasta2AnchorGraph&,
            Shasta2Anchors&,
            uint64_t minCommonForHet,
            double hetErrorRate,
            uint64_t threadCount = 0);

        using Shasta2AnchorGraphBaseClass = boost::adjacency_list<
            boost::listS,
            boost::vecS,
            boost::bidirectionalS,
            boost::no_property,
            Shasta2AnchorGraphEdge>;
}

// An edge's oriented-read-ids live in Shasta2AnchorGraph::orientedReadIds
// (one arena shared by every edge in the graph, exactly like upstream
// shasta2's AnchorGraph), not embedded per edge -- each edge stores only the
// begin/end indices of its own slice. This is a storage-layout choice, not a
// construction-rule difference: a graph with hundreds of thousands of edges
// previously gave each one its own independently heap-allocated
// vector<OrientedReadId>, which is real fragmentation and per-vector
// overhead at that scale. Get a usable Shasta2AnchorPair back out via
// Shasta2AnchorGraph::getAnchorPair(edge_descriptor), which reconstructs one
// (with its own owned copy of the read-id slice, same as before) on demand.
class dinara::Shasta2AnchorGraphEdge {
public:
    Shasta2AnchorId anchorIdA = invalid<Shasta2AnchorId>;
    Shasta2AnchorId anchorIdB = invalid<Shasta2AnchorId>;

    // Begin/end indexes in Shasta2AnchorGraph::orientedReadIds for the
    // OrientedReadIds that belong to this edge.
    uint64_t orientedReadIdsBegin = invalid<uint64_t>;
    uint64_t orientedReadIdsEnd = invalid<uint64_t>;

    uint64_t offset = invalid<uint64_t>;
    uint64_t id = invalid<uint64_t>;
    bool useForAssembly = false;

    // Inter-window edge attributes.
    // supportingSpanPrev/Next: base-level span of the selected read
    // in the source (prev) and destination (next) window.
    // 0 for intra-window edges.
    uint64_t supportingSpanPrev = 0;
    uint64_t supportingSpanNext = 0;
    // sharedReadCount: number of distinct reads that touch both windows
    // (not just the reads using this specific anchor pair).
    uint64_t sharedReadCount = 0;

    Shasta2AnchorGraphEdge(
        Shasta2AnchorId anchorIdA,
        Shasta2AnchorId anchorIdB,
        uint64_t orientedReadIdsBegin,
        uint64_t orientedReadIdsEnd,
        uint64_t offset,
        uint64_t id) :
        anchorIdA(anchorIdA),
        anchorIdB(anchorIdB),
        orientedReadIdsBegin(orientedReadIdsBegin),
        orientedReadIdsEnd(orientedReadIdsEnd),
        offset(offset),
        id(id)
    {}

    Shasta2AnchorGraphEdge() {}

    uint64_t coverage() const {return orientedReadIdsEnd - orientedReadIdsBegin;}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorIdA;
        ar & anchorIdB;
        ar & orientedReadIdsBegin;
        ar & orientedReadIdsEnd;
        ar & offset;
        ar & id;
        ar & useForAssembly;
        ar & supportingSpanPrev;
        ar & supportingSpanNext;
        ar & sharedReadCount;
    }
};



class dinara::Shasta2AnchorGraph :
    public Shasta2AnchorGraphBaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<Shasta2AnchorGraph> {
public:
    using AnchorPairKey = std::pair<Shasta2AnchorId, Shasta2AnchorId>;

    // Construct the AnchorGraph from the Journeys using the same edge
    // creation rule as mode3::AnchorGraph: for each anchor, call
    // Shasta2Anchors::findChildren and create one edge per child that
    // satisfies minEdgeCoverage. The threadCount parameter is accepted
    // for API compatibility but is not used.
    Shasta2AnchorGraph(
        const Shasta2Anchors&,
        const Shasta2Journeys&,
        uint64_t minEdgeCoverage,
        uint64_t threadCount);

    // Default constructor (empty graph).
    Shasta2AnchorGraph() : MultithreadedObject<Shasta2AnchorGraph>(*this) {}

    // Constructor from binary data.
    Shasta2AnchorGraph(const MappedMemoryOwner&, const string& name);

    uint64_t nextEdgeId = 0;

    // To reduce memory fragmentation, the OrientedReadIds of all edges are
    // stored together in this vector (one arena for the whole graph), rather
    // than each edge owning its own vector -- see Shasta2AnchorGraphEdge's
    // comment. Shasta2AnchorGraphEdge::orientedReadIdsBegin/End index into it.
    vector<OrientedReadId> orientedReadIds;

    // Reconstruct a usable Shasta2AnchorPair for an edge (its own owned copy
    // of the read-id slice, exactly as if it had never left one). Mirrors
    // shasta2's own AnchorGraph::getAnchorPair.
    Shasta2AnchorPair getAnchorPair(edge_descriptor e) const
    {
        const Shasta2AnchorGraphEdge& edge = (*this)[e];
        return Shasta2AnchorPair(
            edge.anchorIdA,
            edge.anchorIdB,
            span<const OrientedReadId>(
                orientedReadIds.begin() + edge.orientedReadIdsBegin,
                orientedReadIds.begin() + edge.orientedReadIdsEnd));
    }

    // Append edgeOrientedReadIds to the shared arena and add the edge.
    // Centralizes nextEdgeId assignment and useForAssembly, matching what
    // every construction call site used to do by hand.
    edge_descriptor addEdge(
        Shasta2AnchorId anchorIdA,
        Shasta2AnchorId anchorIdB,
        const vector<OrientedReadId>& edgeOrientedReadIds,
        uint64_t offset,
        bool useForAssembly = true)
    {
        const uint64_t begin = orientedReadIds.size();
        orientedReadIds.insert(orientedReadIds.end(),
            edgeOrientedReadIds.begin(), edgeOrientedReadIds.end());
        const uint64_t end = orientedReadIds.size();
        Shasta2AnchorGraph& anchorGraph = *this;
        edge_descriptor e;
        bool added = false;
        boost::tie(e, added) = boost::add_edge(
            anchorIdA, anchorIdB,
            Shasta2AnchorGraphEdge(anchorIdA, anchorIdB, begin, end, offset, nextEdgeId++),
            anchorGraph);
        anchorGraph[e].useForAssembly = useForAssembly;
        return e;
    }

    void transitiveReduction(
        uint64_t transitiveReductionMaxEdgeCoverage,
        uint64_t maxDistance);
    uint64_t cutWeakStalksLeadingToBranch(
        const Shasta2Anchors& anchors,
        uint64_t maxTipReadCount);
    // Disable an edge and its RC mirror (dst^1 -> src^1).
    void disableEdge(edge_descriptor e);

    // Remove het/hom allele-arm tips: het/hom anchors (id >= hetAnchorFirstId)
    // that end up with active edges on only one side. The staging guarantees
    // every arm is hom-flanked on both sides, but addHetEdge can drop one of an
    // arm's two flank edges (empty read intersection or non-forward offset),
    // leaving the arm hanging. This disables the arm's surviving edge(s) and
    // cascades (a stranded hom can in turn become one-sided), so no het/hom
    // anchor is exported connected on only one side. Backbone/primary anchors
    // are never touched (their legitimate one-sided ends are window/telomere
    // boundaries handled by trimBackbones). Returns the number of edges
    // disabled.
    uint64_t removeHetArmTips(const Shasta2Anchors& anchors);

    // Remove edges between a window and its RC counterpart.
    uint64_t removeRcWindowConnections();

    // Window-level transitive reduction: if A→B→C exists and A→C also exists,
    // remove A→C (the direct edge is redundant).
    uint64_t windowTransitiveReduction();

private:
    bool transitiveReductionCanRemove(edge_descriptor, uint64_t transitiveReductionMaxDistance) const;
public:

    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object<Shasta2AnchorGraphBaseClass>(*this);
        ar & orientedReadIds;
    }
    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    void save(const string& name) const;
    void load(const string& name);

    // Write the graph to GFA format.
    void writeGfa(const string& fileName) const;

    // Write Bandage color CSV.
    void writeCsv(const string& fileName) const;

    // Save binary data (dinara's own format).
    void saveAnchorGraph(const string& name = "Shasta2AnchorGraph") const { save(name); }

    // Export in shasta2-compatible MemoryMapped::Vector<char> format.
    // The output file can be passed to shasta2 via --external-anchor-graph-name.
    // Verifies every serialized edge is forward-monotonic on its shared reads
    // (throws on a violation) so the exported graph can never trip shasta2's
    // LocalAssembly positionB > positionA assertion.
    //
    // dropMap: the SAME per-canonical-anchor member drop set passed to
    // writeExternalAnchors (journey position-tie resolution). A read dropped
    // from an anchor is no longer a member of that anchor in the exported set,
    // so it must also be removed from every edge incident to that anchor;
    // otherwise shasta2's AnchorPair::getAverageOffset -- which requires each
    // edge oriented read to be a common member of BOTH endpoint anchors --
    // asserts (it == orientedReadIds.end()). Edges left with zero reads after
    // filtering are skipped (a zero-read AnchorPair divides by zero).
    void saveForShasta2(
        const string& fileName,
        const Shasta2Anchors& anchors,
        const Shasta2Anchors::ExternalAnchorDropMap* dropMap = nullptr) const;

    // Per-anchor window assignment (populated by the anchor-window constructor).
    // Maps anchorId -> windowId. noWindow means unmapped.
    // Includes RC mirror windows (windowId >= windowCount).
    static constexpr uint32_t noWindow = std::numeric_limits<uint32_t>::max();
    vector<uint32_t> anchorToWindow;
    uint32_t windowCount = 0;

    // Per-window: distinct oriented read IDs that touch the window.
    std::map<uint32_t, std::set<uint32_t>> windowReads;

    // Per-read: ordered window sequence (consecutive duplicates removed).
    // Gives the order in which the read visits windows.
    std::map<uint32_t, vector<uint32_t>> readWindows;

};
