#pragma once

#include "Shasta2AnchorPair.hpp"
#include "Shasta2Anchors.hpp"
#include "Base.hpp"
#include "Kmer.hpp"
#include "Marker.hpp"
#include "ReadId.hpp"
#include "invalid.hpp"

#include <boost/graph/adjacency_list.hpp>

#include "iosfwd.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace dinara {
    class Shasta2LocalAssembly3;
    class Shasta2LocalAssembly3Vertex;
    class Shasta2LocalAssembly3Edge;

    using Markers = MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>;

    using Shasta2LocalAssembly3BaseClass = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Shasta2LocalAssembly3Vertex,
        Shasta2LocalAssembly3Edge>;
}

class dinara::Shasta2LocalAssembly3Vertex {
public:
    uint64_t kmerIndex;

    class Data {
    public:
        uint64_t orientedReadIndex;
        uint32_t ordinal;
    };
    vector<Data> data;
    uint64_t coverage() const
    {
        return data.size();
    }

    Shasta2LocalAssembly3Vertex(uint64_t kmerIndex) : kmerIndex(kmerIndex) {}

    Shasta2LocalAssembly3BaseClass::vertex_descriptor dominator =
        Shasta2LocalAssembly3BaseClass::null_vertex();
    bool isOnDominatorTreePath = false;
};

class dinara::Shasta2LocalAssembly3Edge {
public:
    class Data {
    public:
        uint64_t orientedReadIndex;
        uint32_t ordinal0;
        uint32_t ordinal1;
    };
    vector<Data> data;
    uint64_t coverage() const
    {
        return data.size();
    }
};

class dinara::Shasta2LocalAssembly3 : public Shasta2LocalAssembly3BaseClass {
public:
    using Base = dinara::Base;

    Shasta2LocalAssembly3(
        const Shasta2Anchors&,
        uint64_t abpoaMaxLength,
        ostream& html,
        bool debug,
        const Shasta2AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds);

public:
    Shasta2AnchorId leftAnchorId;
    Shasta2AnchorId rightAnchorId;

    class OrientedReadInfo {
    public:
        OrientedReadId orientedReadId;

        bool isOnAnchorPair = false;
        bool isOnLeftAnchor = false;
        bool isOnRightAnchor = false;
        bool isOnBothAnchors() const
        {
            return isOnLeftAnchor and isOnRightAnchor;
        }

        uint32_t leftOrdinal = invalid<uint32_t>;
        uint32_t rightOrdinal = invalid<uint32_t>;
        uint32_t ordinalOffset() const
        {
            DINARA_ASSERT(isOnLeftAnchor);
            DINARA_ASSERT(isOnRightAnchor);
            return rightOrdinal - leftOrdinal;
        }

        uint32_t leftPosition = invalid<uint32_t>;
        uint32_t rightPosition = invalid<uint32_t>;
        uint32_t positionOffset() const
        {
            DINARA_ASSERT(isOnLeftAnchor);
            DINARA_ASSERT(isOnRightAnchor);
            return rightPosition - leftPosition;
        }

        uint32_t firstOrdinalForAssembly = invalid<uint32_t>;
        uint32_t lastOrdinalForAssembly = invalid<uint32_t>;
        void fillFirstLastOrdinalForAssembly(const Markers&, uint32_t length);
        uint32_t firstPositionForAssembly(const Markers&) const;
        uint32_t lastPositionForAssembly(const Markers&) const;

        class OrientedReadKmerInfo {
        public:
            Kmer kmer;
            uint64_t kmerIndex;
        };
        vector<OrientedReadKmerInfo> orientedReadKmerInfos;
        void fillOrientedReadKmers(const Shasta2Anchors&);

        bool operator<(const OrientedReadInfo& that) const
        {
            return orientedReadId < that.orientedReadId;
        }
    };

    vector<OrientedReadInfo> orientedReadInfos;

    void gatherOrientedReadsOnAnchorPair(
        const Shasta2Anchors&,
        const Shasta2AnchorPair&);

    uint32_t offset;
    void estimateOffset();

    void gatherAdditionalOrientedReads(
        const Shasta2Anchors&,
        const Shasta2AnchorPair&,
        const vector<OrientedReadId>& additionalOrientedReadIds,
        double drift);

    void fillFirstLastOrdinalForAssembly(const Markers&, double drift);
    void fillOrientedReadKmers(const Shasta2Anchors&);

    vector<Kmer> kmers;
    void gatherKmers(const Shasta2Anchors&);
    uint64_t getKmerIndex(const Kmer&) const;

    uint64_t leftAnchorKmerIndex;
    uint64_t rightAnchorKmerIndex;

    vector<vertex_descriptor> vertexMap;
    void createVertices();
    void createEdges();
    vertex_descriptor leftAnchorVertex;
    vertex_descriptor rightAnchorVertex;

    uint32_t edgePositionOffset(edge_descriptor, const Markers&) const;

    void computeDominatorTree();
    vector<vertex_descriptor> dominatorTreePath;

    void assemble(
        const Shasta2Anchors&,
        uint64_t abpoaMaxLength,
        ostream& html,
        bool debug);
    void assemble(
        const Shasta2Anchors&,
        vertex_descriptor,
        vertex_descriptor,
        uint64_t abpoaMaxLength,
        ostream& html,
        bool debug);

    class AssemblyInfo {
    public:
        OrientedReadId orientedReadId;
        uint32_t ordinal0;
        uint32_t ordinal1;
        uint32_t position0;
        uint32_t position1;
        vector<Base> sequence;
    };

    vector<Base> sequence;
    vector<uint64_t> coverage;

    void writeInput(
        ostream& html,
        bool debug,
        const Shasta2AnchorPair& anchorPair,
        const vector<OrientedReadId>& additionalOrientedReadIds) const;
    void writeOrientedReads(const Shasta2Anchors&, ostream& html) const;
    void writeKmers(ostream& html, uint64_t k) const;
    void writeOrientedReadKmers(ostream& html) const;
    void writeAssemblyInfos(ostream& html, const vector<AssemblyInfo>&) const;
    void writeGraphviz(const string& fileName, const Markers&) const;
    void writeGraphviz(ostream&, const Markers&) const;
    void writeHtml(ostream&, const Markers&) const;
    void writeConsensus(
        ostream& html,
        const vector< pair<Base, uint64_t> >& consensus,
        uint64_t maxCoverage) const;
    void writeAlignment(
        ostream& html,
        const vector< pair<vector<Base>, uint64_t> >& sequenceWithCoverage,
        const vector< pair<Base, uint64_t> >& consensus,
        const vector< vector<AlignedBase> >& alignment,
        const vector<AlignedBase>& alignedConsensus,
        const vector<AssemblyInfo>&) const;
};
