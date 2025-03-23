#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <fstream>
#include <unordered_map>
#include <ctime>      // For clock()
#include <iomanip>    // For setprecision
using namespace std;
using Graph = vector<set<int>>;

int largest_clique_size = 0;
int total_maximal_cliques = 0;
map<int, int> clique_size_distribution;

void CLIQUE(const Graph& graph, set<int>&& P, set<int>&& R, set<int>&& X) {
    if (P.empty() && X.empty()) {
        int clique_size = R.size();
        largest_clique_size = max(largest_clique_size, clique_size);
        total_maximal_cliques++;
        clique_size_distribution[clique_size]++;
        return;
    }
    
    set<int> P_union_X;
    merge(P.begin(), P.end(),
          X.begin(), X.end(),
          inserter(P_union_X, P_union_X.end()));
    
    if (!P_union_X.empty()) {
        int u = *P_union_X.begin();
        for (int v : P_union_X) {
            if (graph[v].size() > graph[u].size()) {
                u = v;
            }
        }
        
        set<int> P_minus_N_u;
        set_difference(P.begin(), P.end(),
                      graph[u].begin(), graph[u].end(),
                      inserter(P_minus_N_u, P_minus_N_u.end()));
        
        for (int v : P_minus_N_u) {
            set<int> new_P, new_X;
            const set<int>& N_v = graph[v];
            
            set_intersection(P.begin(), P.end(),
                            N_v.begin(), N_v.end(),
                            inserter(new_P, new_P.end()));
            
            set_intersection(X.begin(), X.end(),
                            N_v.begin(), N_v.end(),
                            inserter(new_X, new_X.end()));
            
            set<int> new_R(R);
            new_R.insert(v);
            
            CLIQUE(graph, 
                  move(new_P), 
                  move(new_R), 
                  move(new_X));
            
            P.erase(v);
            X.insert(v);
        }
    }
}

void findAllMaximalCliques(const Graph& graph) {
    set<int> P;
    for (int i = 0; i < graph.size(); i++) {
        P.insert(i);
    }
    CLIQUE(graph, move(P), {}, {});
}

int main() {
    ifstream inputFile("as-skitter (1).txt");
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }
    
    // Node mapping structures
    unordered_map<int, int> nodeToIndex;
    vector<set<int>> graph;
    int currentIndex = 0;
    
    // Read and process edges
    int reported_n, m;
    inputFile >> reported_n >> m;  // Read and ignore reported node count
    
    // Start timing
    clock_t start = clock();
    
    for (int i = 0; i < m; i++) {
        int u, v;
        inputFile >> u >> v;
        
        // Map nodes to indices
        auto mapNode = [&](int node) {
            if (!nodeToIndex.count(node)) {
                nodeToIndex[node] = currentIndex++;
                graph.emplace_back();
            }
            return nodeToIndex[node];
        };
        
        int uIndex = mapNode(u);
        int vIndex = mapNode(v);
        
        // Add undirected edges
        graph[uIndex].insert(vIndex);
        graph[vIndex].insert(uIndex);
    }
    inputFile.close();
    
    // Resize graph in case some nodes had no edges
    graph.resize(nodeToIndex.size());
    
    findAllMaximalCliques(graph);
    
    // End timing
    clock_t end = clock();
    double duration = double(end - start) / CLOCKS_PER_SEC;
    
    cout << "1. Largest clique size: " << largest_clique_size << endl;
    cout << "2. Total maximal cliques: " << total_maximal_cliques << endl;
    cout << "3. Size distribution:" << endl;
    for (const auto& pair : clique_size_distribution) {
        cout << " Size " << pair.first << ": " << pair.second << endl;
    }
    
    // Print runtime
    cout << fixed << setprecision(4);
    cout << "4. Algorithm runtime: " << duration << " seconds" << endl;
    
    return 0;
}