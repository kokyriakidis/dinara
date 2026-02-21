#pragma once

// Dinara.
#include "Base.hpp"
#include "DINARA_ASSERT.hpp"
#include "invalid.hpp"
#include "Kmer.hpp"
#include "MarkerKmers.hpp"
#include "ReadId.hpp"
#include "Shasta2Anchors.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "span.hpp"
#include "vector.hpp"



namespace dinara {
    class Shasta2LocalAssembly5;
    class Shasta2LocalAssembly5Vertex;
    class Shasta2LocalAssembly5Edge;

    using Shasta2LocalAssembly5BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Shasta2LocalAssembly5Vertex,
        Shasta2LocalAssembly5Edge>;

    class Shasta2Anchors;
    class Shasta2AnchorPair;
    class Base;
    class Marker;
}



class dinara::Shasta2LocalAssembly5Vertex {
public:

    // The id of the Kmer that corresponds to this vertex.
    // This is the index of this Kmer in the Shasta2LocalAssembly5::kmers vector.
    uint64_t kmerId;

    class Info {
    public:

        // Index in orientedReadInfos.
        uint64_t orientedReadIndex;

        // Ordinal.
        uint32_t ordinal;
    };
    vector<Info> infos;

    bool isOnDominatorTreePath = false;

    // These are used for approximate topological ordering.
    // The color is also used by computeAssemblyPath.
    uint64_t color = invalid<uint64_t>;
    uint64_t rank = invalid<uint64_t>;

    Shasta2LocalAssembly5Vertex(uint64_t kmerId = invalid<uint64_t>) :
        kmerId(kmerId) {}
};



// For edge v0->v1 we store:
// - vector<Info> infos: using all the reads that visit v1
//   immediately after visiting v0.
// - vector<Info> extendedInfos: using all the reads that visit v1
//   after visiting v0, but not immediately after.
class dinara::Shasta2LocalAssembly5Edge {
public:
    class Info {
    public:

        // Index in orientedReadInfos.
        uint64_t orientedReadIndex;

        // Ordinal.
        uint32_t ordinal0;

        // Ordinal.
        uint32_t ordinal1;
    };
    vector<Info> infos;
    vector<Info> extendedInfos;

    // This is used for approximate topological ordering.
    // which is only used by writeGraphviz.
    bool isDagEdge = false;

    bool isOnAssemblyPath = false;

    // Default constructor.
    Shasta2LocalAssembly5Edge() {}

};



class dinara::Shasta2LocalAssembly5 : public Shasta2LocalAssembly5BaseClass {
public:
    // Hide boost::adjacency_list::Base;
    using Base = dinara::Base;

    // This assembles between anchorIdA and anchorIdB of the given Shasta2AnchorPair.
    // This can use all of the OrientedReadIds in the Shasta2AnchorPair
    // and/or in the additionalOrientedReadIds that also appear in the
    // Kmers corresponding to the left OR right Anchors.
    // In contrast, Shasta2LocalAssembly4 uses oriented reads that appear in the
    // Kmers corresponding to the left AND right Anchors.
    // At sufficiently high coverage, Shasta2LocalAssembly4 is a better choice
    // because it is simpler and faster. Shasta2LocalAssembly5 is better
    // at low coverage, when it is important to use as many reads as possible,
    // The additionalOrientedReadIds must be sorted.
    // The additionalOrientedReadIds are allowed to contain OrientedReadIds
    // that are also in the Shasta2AnchorPair.
    Shasta2LocalAssembly5(
        const Shasta2Anchors&,
        uint64_t abpoaMaxLength,
        ostream& html,
        bool debug,
        const Shasta2AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds);

    // Assembled sequence and its coverage.
    vector<Base> sequence;
    vector<uint64_t> coverage;

private:

    // EXPOSE WHEN CODE STABILIZES.
    const double drift = 0.1;
    const uint64_t minVertexCoverage = 3;

    const Shasta2Anchors& anchors;
    uint64_t abpoaMaxLength;
    ostream& html;

    // The two anchors of the Shasta2AnchorPair used for this assembly.
    Shasta2AnchorId leftAnchorId;
    Shasta2AnchorId rightAnchorId;

    // The corresponding Kmers.
    Kmer leftKmer;
    Kmer rightKmer;

    // MarkerInfos for the Marker Kmers corresponding to the left
    // and right anchors of the Shasta2AnchorPair being assembled.
    vector<MarkerKmers::MarkerInfo> leftMarkerInfos;
    vector<MarkerKmers::MarkerInfo> rightMarkerInfos;
    void fillMarkerInfos();

    // The union of the OrientedReadIds in the Shasta2AnchorPair and
    // the additionalOrientedReadIds.
    vector<OrientedReadId> allOrientedReadIds;
    void gatherAllOrientedReads(
        const Shasta2AnchorPair&,
        const vector<OrientedReadId>& additionalOrientedReadIds);



    // The oriented reads used in this local assembly.
    // These are all OrientedReadIds that are in allOrientedReadIds and
    // also on the left or right Marker.
    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        // Whether this OrientedReadId appear in the left and right Markers.
        bool isOnLeftMarker = false;
        bool isOnRightMarker = false;
        bool isOnBothMarkers() const
        {
            return isOnLeftMarker and isOnRightMarker;
        }

        // The ordinals of this OrientedReadId in the left/right Markers, if any.
        uint32_t leftOrdinal = invalid<uint32_t>;
        uint32_t rightOrdinal = invalid<uint32_t>;
        uint32_t ordinalOffset() const
        {
            DINARA_ASSERT(isOnLeftMarker);
            DINARA_ASSERT(isOnRightMarker);
            DINARA_ASSERT(rightOrdinal > leftOrdinal);
            return rightOrdinal - leftOrdinal;
        }

        // The base positions of this OrientedReadId's marker
        // in the left/right Markers, if any.
        // These are positions in the oriented read sequence
        // of the mid point of the marker.
        uint32_t leftPosition = invalid<uint32_t>;
        uint32_t rightPosition = invalid<uint32_t>;
        uint32_t positionOffset() const
        {
            DINARA_ASSERT(isOnLeftMarker);
            DINARA_ASSERT(isOnRightMarker);
            DINARA_ASSERT(rightPosition > leftPosition);
            return rightPosition - leftPosition;
        }



        // The region of this oriented read that can be used in this local assembly.
        // The positions are positions of the marker midpoints.
        class LocalRegion {
        public:
            uint32_t firstOrdinal = invalid<uint32_t>;
            uint32_t lastOrdinal = invalid<uint32_t>;
            uint32_t firstPosition = invalid<uint32_t>;
            uint32_t lastPosition = invalid<uint32_t>;

            uint32_t ordinalOffset() const
            {
                return lastOrdinal - firstOrdinal;
            }

            uint32_t positionOffset() const
            {
                return lastPosition - firstPosition;
            }

            // The Kmers for all ordinals in [firstOrdinalForAssembly, lastOrdinalForAssembly],
            // in that order.
            vector<Kmer> kmers;

            const Kmer& getKmer(uint32_t ordinal) const
            {
                DINARA_ASSERT(ordinal >= firstOrdinal);
                DINARA_ASSERT(ordinal <= lastOrdinal);
                return kmers[ordinal - firstOrdinal];
            }

            // The Kmers that appear more than once in this LocalRegion, sorted by Kmer.
            vector<Kmer> nonUniqueKmers;
            void computeNonUniqueKmers();
        };
        LocalRegion localRegion;
        void fillLocalRegion(const Shasta2Anchors&, uint32_t length);


        // The ordinals of the Marker Kmers to be used for assembly.
        // Initially these are all the ordinals in [firstOrdinalForAssembly, lastOrdinalForAssembly],
        // except for the ones that correspond to Kmers that are non-unique in one or more reads.
        // Later, we remove ordinals that causes cycles in the graph.
        // These are sorted.
        vector<uint32_t> ordinalsForAssembly;

    };
    vector<OrientedReadInfo> orientedReadInfos;
    void gatherOrientedReads();

    uint32_t offset;
    void estimateOffset();

    // This fills the LocalRegion of all OrientedReadInfos.
    void fillLocalRegions();

    // The Kmers that will be used in this assembly.
    vector<Kmer> kmers;

    // Gather the Kmers that will be used in this assembly
    // and fill in the OrientedReadInfos;
    void gatherKmers();


    // Create the graph using the current ordinalsForAssembly.
    vertex_descriptor vLeft;
    vertex_descriptor vRight;
    void createGraph();

    // Compute position offset for an edge by averaging the
    // position offsets of the participating orientedReads.
    uint32_t edgeOffset(edge_descriptor) const;

    // Remove vertices that are not forward accessible from vLeft
    // and backward accessible from vRight.
    void removeInaccessibleVertices();

    // Compute a dominator tree starting at vLeft and the
    // dominator tree path between vLeft and vRight.
    void computeDominatorTree();
    vector<vertex_descriptor> dominatorTreePath;

    // Given two vertices which are adjacent in the dominator tree path, vA and vB,
    // remove low coverage edges in-between without destroying reachability
    // of vB from vA.
    void removeLowCoverageEdges(vertex_descriptor vA, vertex_descriptor vB);
    void removeLowCoverageEdges();

    // Compute the assembly path.
    void computeAssemblyPath();
    vector<edge_descriptor> assemblyPath;

    void removeIsolatedVertices();

    // Assemble sequence.
    void assemble();
    void assemble(edge_descriptor);


    // Html output.
    void writeInput(
        const Shasta2AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds) const;
    void writeAllOrientedReadIds() const;
    void writeMarkerInfos() const;
    void writeMarkerInfos(const string& side, const vector<MarkerKmers::MarkerInfo>&) const;
    void writeOrientedReads() const;
    void writeGraphviz(const string& fileName);
    void writeGraphviz(ostream&);
    void writeGraph();
    void writeKmers() const;
    void writeAssembledSequence() const;


    // A condensed version of the Shasta2LocalAssembly5 graph, used
    // to compute the assembly path. Each strongly connected component
    // of the Shasta2LocalAssembly5 graph is collapsed into a single vertex.
    // The CondensedGraph is guaranteed to be acyclic. See:
    // https://en.wikipedia.org/wiki/Strongly_connected_component#Definitions
    // https://cp-algorithms.com/graph/strongly-connected-components.html
    class CondensedGraphVertex;
    using CondensedGraphBaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        CondensedGraphVertex>;
    class CondensedGraphVertex {
    public:
        vector<vertex_descriptor> vertices;
        CondensedGraphBaseClass::vertex_descriptor dominator = null_vertex();
    };
    class CondensedGraph : public CondensedGraphBaseClass {
    };

};
