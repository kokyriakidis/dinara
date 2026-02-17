#include "Shasta2SuperbubbleChain.hpp"

#include "Shasta2GTest.hpp"
#include "Shasta2PhasingGraph.hpp"
#include "Shasta2RestrictedAnchorGraph.hpp"
#include "Shasta2Tangle1.hpp"
#include "Shasta2TangleMatrix1.hpp"
#include "DINARA_ASSERT.hpp"
#include "deduplicate.hpp"

#include <boost/graph/iteration_macros.hpp>

#include <algorithm>
#include <limits>
#include <mutex>
#include <set>

using namespace dinara;
using namespace std;

uint64_t Shasta2SuperbubbleChain::phase1(
    Shasta2AssemblyGraph& assemblyGraph,
    const uint64_t superbubbleChainId) const
{
    (void)superbubbleChainId;
    if((assemblyGraph.getAnchorsPointer() == nullptr) || (assemblyGraph.getJourneysPointer() == nullptr)) {
        return 0;
    }

    const uint64_t phasingDistance = assemblyGraph.getOptions().phasingDistance;

    Shasta2PhasingGraph phasingGraph;
    for(uint64_t position=0; position<size(); position++) {
        const Shasta2Superbubble& bubble = at(position);
        if(!bubble.isBubble() || bubble.isTrivial()) {
            continue;
        }
        phasingGraph.addVertex(position);
    }

    for(uint64_t position0=0; position0<size(); position0++) {
        const Shasta2Superbubble& bubble0 = at(position0);
        if(!bubble0.isBubble() || bubble0.isTrivial()) {
            continue;
        }

        uint64_t n0 = 0;
        for(uint64_t position1=position0+1; position1<size(); position1++) {
            const Shasta2Superbubble& bubble1 = at(position1);
            if(!bubble1.isBubble() || bubble1.isTrivial()) {
                continue;
            }
            if(n0 > phasingDistance) {
                continue;
            }
            ++n0;

            ostream html(0);
            Shasta2TangleMatrix1 tangleMatrix(
                assemblyGraph,
                bubble0.internalEdges,
                bubble1.internalEdges,
                html);

            const Shasta2GTest gTest(
                tangleMatrix.tangleMatrix,
                assemblyGraph.getOptions().detangleEpsilon,
                true,
                true);
            if(!gTest.success) {
                continue;
            }

            const auto& bestHypothesis = gTest.hypotheses.front();
            const double bestG = bestHypothesis.G;
            if(bestG > assemblyGraph.getOptions().detangleMaxLogP) {
                continue;
            }
            if(gTest.hypotheses.size() > 1) {
                const double secondBestG = gTest.hypotheses[1].G;
                if(secondBestG - bestG < assemblyGraph.getOptions().detangleMinLogPDelta) {
                    continue;
                }
            }

            if(!(
                Shasta2GTest::isForwardInjective(bestHypothesis.connectivityMatrix) &&
                Shasta2GTest::isBackwardInjective(bestHypothesis.connectivityMatrix)
                )) {
                continue;
            }

            bool coverageCheckFailed = false;
            for(uint64_t iEntrance=0; iEntrance<tangleMatrix.entrances.size(); iEntrance++) {
                const Shasta2AssemblyGraph::vertex_descriptor v0 = target(tangleMatrix.entrances[iEntrance], assemblyGraph);
                const Shasta2AnchorId anchorId0 = assemblyGraph[v0].anchorId;
                for(uint64_t iExit=0; iExit<tangleMatrix.exits.size(); iExit++) {
                    if(!bestHypothesis.connectivityMatrix[iEntrance][iExit]) {
                        continue;
                    }

                    const Shasta2AssemblyGraph::vertex_descriptor v1 = source(tangleMatrix.exits[iExit], assemblyGraph);
                    const Shasta2AnchorId anchorId1 = assemblyGraph[v1].anchorId;
                    if(anchorId1 == anchorId0) {
                        continue;
                    }

                    try {
                        ostream html2(0);
                        Shasta2RestrictedAnchorGraph restrictedAnchorGraph(
                            *assemblyGraph.getAnchorsPointer(),
                            *assemblyGraph.getJourneysPointer(),
                            tangleMatrix,
                            iEntrance,
                            iExit,
                            html2);

                        vector<Shasta2RestrictedAnchorGraph::edge_descriptor> longestPath;
                        restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

                        uint64_t minCoverage = numeric_limits<uint64_t>::max();
                        for(const auto e: longestPath) {
                            const auto& edge = restrictedAnchorGraph[e];
                            minCoverage = min(minCoverage, edge.anchorPair.size());
                        }
                        if(minCoverage == 0) {
                            coverageCheckFailed = true;
                        }
                    } catch(const std::exception&) {
                        coverageCheckFailed = true;
                    }

                    if(coverageCheckFailed) {
                        break;
                    }
                }
                if(coverageCheckFailed) {
                    break;
                }
            }
            if(coverageCheckFailed) {
                continue;
            }

            phasingGraph.addEdge(position0, position1, bestHypothesis);
        }
    }

    phasingGraph.removeLowDegreeVertices(1);
    if(num_vertices(phasingGraph) < 2) {
        return 0;
    }

    phasingGraph.computeConnectedComponents();
    phasingGraph.findLongestPaths();

    uint64_t changeCount = 0;

    lock_guard<mutex> lock(assemblyGraph.getMutex());

    vector<Shasta2AssemblyGraph::vertex_descriptor> removedVertices;
    for(uint64_t componentId=0; componentId<phasingGraph.longestPaths.size(); componentId++) {
        const vector<Shasta2PhasingGraph::edge_descriptor>& longestPath = phasingGraph.longestPaths[componentId];
        for(const Shasta2PhasingGraph::edge_descriptor e: longestPath) {
            const Shasta2PhasingGraph::vertex_descriptor v0 = source(e, phasingGraph);
            const Shasta2PhasingGraph::vertex_descriptor v1 = target(e, phasingGraph);

            const uint64_t position0 = phasingGraph[v0].position;
            const uint64_t position1 = phasingGraph[v1].position;

            const Shasta2Superbubble& bubble0 = at(position0);
            const Shasta2Superbubble& bubble1 = at(position1);

            vector<Shasta2AssemblyGraph::vertex_descriptor> tangleVertices;
            tangleVertices.push_back(bubble0.targetVertex);
            tangleVertices.push_back(bubble1.sourceVertex);
            for(uint64_t position=position0+1; position<position1; position++) {
                const Shasta2Superbubble& superbubble = at(position);
                tangleVertices.push_back(superbubble.sourceVertex);
                tangleVertices.push_back(superbubble.targetVertex);
                copy(
                    superbubble.internalVertices.begin(),
                    superbubble.internalVertices.end(),
                    back_inserter(tangleVertices));
            }

            sort(
                tangleVertices.begin(),
                tangleVertices.end(),
                [&assemblyGraph](
                    const Shasta2AssemblyGraph::vertex_descriptor a,
                    const Shasta2AssemblyGraph::vertex_descriptor b) {
                    return assemblyGraph[a].id < assemblyGraph[b].id;
                });
            tangleVertices.erase(
                unique(
                    tangleVertices.begin(),
                    tangleVertices.end(),
                    [&assemblyGraph](
                        const Shasta2AssemblyGraph::vertex_descriptor a,
                        const Shasta2AssemblyGraph::vertex_descriptor b) {
                        return assemblyGraph[a].id == assemblyGraph[b].id;
                    }),
                tangleVertices.end());

            bool vertexWasRemoved = false;
            for(const auto tv: tangleVertices) {
                const auto it = find_if(
                    removedVertices.begin(),
                    removedVertices.end(),
                    [&assemblyGraph, tv](const Shasta2AssemblyGraph::vertex_descriptor x) {
                        return assemblyGraph[x].id == assemblyGraph[tv].id;
                    });
                if(it != removedVertices.end()) {
                    vertexWasRemoved = true;
                    break;
                }
            }
            if(vertexWasRemoved) {
                continue;
            }

            Shasta2Tangle1 tangle(assemblyGraph, tangleVertices);
            DINARA_ASSERT(tangle.tangleMatrix().entrances.size());

            const Shasta2GTest::Hypothesis& bestHypothesis = phasingGraph[e].bestHypothesis;
            const vector< vector<bool> >& connectivityMatrix = bestHypothesis.connectivityMatrix;
            for(uint64_t iEntrance=0; iEntrance<tangle.tangleMatrix().entrances.size(); iEntrance++) {
                for(uint64_t iExit=0; iExit<tangle.tangleMatrix().exits.size(); iExit++) {
                    if(connectivityMatrix[iEntrance][iExit]) {
                        DINARA_ASSERT(tangle.addConnectPair(iEntrance, iExit));
                    }
                }
            }

            tangle.detangle();
            ++changeCount;

            for(const auto rv: tangle.removedVertices) {
                removedVertices.push_back(rv);
            }
            sort(
                removedVertices.begin(),
                removedVertices.end(),
                [&assemblyGraph](
                    const Shasta2AssemblyGraph::vertex_descriptor a,
                    const Shasta2AssemblyGraph::vertex_descriptor b) {
                    return assemblyGraph[a].id < assemblyGraph[b].id;
                });
            removedVertices.erase(
                unique(
                    removedVertices.begin(),
                    removedVertices.end(),
                    [&assemblyGraph](
                        const Shasta2AssemblyGraph::vertex_descriptor a,
                        const Shasta2AssemblyGraph::vertex_descriptor b) {
                        return assemblyGraph[a].id == assemblyGraph[b].id;
                    }),
                removedVertices.end());
        }
    }

    return changeCount;
}
