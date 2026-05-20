// Standalone superbubble finder for GFA files.
// Uses the Onodera et al. 2013 algorithm (same as Verkko/rukki).
//
// Usage: findSuperbubbleOnoderaGfa input.gfa [maxCount]
//
// Reads a GFA1 file, builds a directed graph using Verkko's convention
// (each segment becomes two oriented nodes: >name and <name),
// finds all superbubbles, and writes them to stdout.

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;



struct GfaVertex {
    string name; // e.g. ">seg1" or "<seg1"
};

using GfaGraph = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::bidirectionalS, GfaVertex>;
using vertex_descriptor = GfaGraph::vertex_descriptor;



// Onodera et al. 2013 superbubble detection.
static vertex_descriptor findSuperbubble(
    const GfaGraph& graph,
    vertex_descriptor vStart,
    uint64_t maxCount)
{
    // Start vertex must have out-degree >= 2 (excluding self-loops).
    {
        uint64_t nonSelfOutDegree = 0;
        BGL_FORALL_OUTEDGES(vStart, e, graph, GfaGraph) {
            if(target(e, graph) != vStart) {
                ++nonSelfOutDegree;
            }
        }
        if(nonSelfOutDegree < 2) {
            return GfaGraph::null_vertex();
        }
    }

    unordered_set<vertex_descriptor> visited;
    unordered_map<vertex_descriptor, uint64_t> remainingIncoming;
    vector<vertex_descriptor> ready;
    uint64_t notReadyCount = 0;

    ready.push_back(vStart);
    remainingIncoming[vStart] = 0;

    while(!ready.empty()) {
        if(maxCount > 0 && visited.size() + notReadyCount > maxCount) {
            return GfaGraph::null_vertex();
        }

        const vertex_descriptor v = ready.back();
        ready.pop_back();
        visited.insert(v);

        if(out_degree(v, graph) == 0) {
            return GfaGraph::null_vertex();
        }

        BGL_FORALL_OUTEDGES(v, e, graph, GfaGraph) {
            const vertex_descriptor w = target(e, graph);

            if(w == v) return GfaGraph::null_vertex();
            if(w == vStart) return GfaGraph::null_vertex();

            if(remainingIncoming.find(w) == remainingIncoming.end()) {
                notReadyCount++;
                remainingIncoming[w] = in_degree(w, graph);
            }

            auto& rem = remainingIncoming[w];
            rem--;

            if(rem == 0) {
                ready.push_back(w);
                notReadyCount--;
            }
        }

        if(ready.size() == 1 && notReadyCount == 0) {
            const vertex_descriptor t = ready.back();

            // Reject if exit has edge back to start.
            BGL_FORALL_OUTEDGES(t, e, graph, GfaGraph) {
                if(target(e, graph) == vStart) {
                    return GfaGraph::null_vertex();
                }
            }

            return t;
        }
    }

    return GfaGraph::null_vertex();
}



int main(int argc, char* argv[])
{
    if(argc < 2) {
        cerr << "Usage: " << argv[0] << " input.gfa [maxCount]" << endl;
        return 1;
    }

    const string gfaPath = argv[1];
    const uint64_t maxCount = (argc >= 3) ? std::stoull(argv[2]) : 0;

    ifstream gfa(gfaPath);
    if(!gfa) {
        cerr << "Cannot open " << gfaPath << endl;
        return 1;
    }

    GfaGraph graph;
    unordered_map<string, vertex_descriptor> nameToVertex;

    auto getOrCreate = [&](const string& name) -> vertex_descriptor {
        auto it = nameToVertex.find(name);
        if(it != nameToVertex.end()) return it->second;
        const vertex_descriptor v = add_vertex(graph);
        graph[v].name = name;
        nameToVertex[name] = v;
        return v;
    };

    // Verkko convention: >name = forward, <name = reverse.
    auto orientedName = [](const string& seg, const string& orient) -> string {
        return (orient == "+" ? ">" : "<") + seg;
    };

    auto revName = [](const string& name) -> string {
        return (name[0] == '>' ? "<" : ">") + name.substr(1);
    };

    string line;
    while(std::getline(gfa, line)) {
        if(line.empty()) continue;

        std::istringstream iss(line);
        string recordType;
        iss >> recordType;

        if(recordType == "S") {
            string name;
            iss >> name;
            getOrCreate(">" + name);
            getOrCreate("<" + name);
        } else if(recordType == "L") {
            string seg1, orient1, seg2, orient2;
            iss >> seg1 >> orient1 >> seg2 >> orient2;

            // Forward edge.
            string fromName = orientedName(seg1, orient1);
            string toName = orientedName(seg2, orient2);
            add_edge(getOrCreate(fromName), getOrCreate(toName), graph);

            // Reverse-complement edge.
            string rcFrom = revName(toName);
            string rcTo = revName(fromName);
            add_edge(getOrCreate(rcFrom), getOrCreate(rcTo), graph);
        }
    }

    cerr << "Graph: " << num_vertices(graph) << " oriented nodes, "
         << num_edges(graph) << " edges" << endl;

    // Find superbubbles.
    uint64_t count = 0;
    BGL_FORALL_VERTICES(vSource, graph, GfaGraph) {
        const vertex_descriptor vTarget = findSuperbubble(graph, vSource, maxCount);
        if(vTarget != GfaGraph::null_vertex()) {
            cout << graph[vSource].name << "\t" << graph[vTarget].name << endl;
            count++;
        }
    }

    cerr << "Found " << count << " superbubbles" << endl;
    return 0;
}
