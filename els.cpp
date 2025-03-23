#include <bits/stdc++.h>
#include <fstream>
#include <ctime> // Include ctime for clock function
using namespace std;

using Graph = vector<unordered_set<int>>;
int largest_clique_size = 0;
int total_maximal_cliques = 0;
map<int, int> clique_size_distribution;

// Function to compute the degeneracy ordering
vector<int> degeneracyOrdering(const Graph& graph) {
    int n = graph.size();
    vector<int> degree(n);
    vector<unordered_set<int>> bin(n);
    vector<int> vert(n);
    vector<int> pos(n);
    
    for (int i = 0; i < n; ++i) {
        degree[i] = graph[i].size();
        bin[degree[i]].insert(i);
        vert[i] = i;
    }
    
    int degeneracy = 0;
    int j = 0;
    for (int i = 0; i < n; ++i) {
        while (bin[j].empty()) ++j;
        degeneracy = max(degeneracy, j);
        int v = *bin[j].begin();
        bin[j].erase(v);
        pos[v] = i;
        
        for (int u : graph[v]) {
            if (pos[u] > i) {
                bin[degree[u]].erase(u);
                --degree[u];
                bin[degree[u]].insert(u);
            }
        }
    }
    
    vector<int> ordering(n);
    for (int i = 0; i < n; ++i) {
        ordering[pos[i]] = i;
    }
    return ordering;
}

void BronKerboschPivot(const Graph& graph, unordered_set<int>& P, unordered_set<int> R, unordered_set<int> X) {
    if (P.empty() && X.empty()) {
        // Report maximal clique
        int clique_size = R.size();
        largest_clique_size = max(largest_clique_size, clique_size);
        total_maximal_cliques++;
        clique_size_distribution[clique_size]++;
        return;
    }
    
    unordered_set<int> P_union_X;
    P_union_X.insert(P.begin(), P.end());
    P_union_X.insert(X.begin(), X.end());
    
    int u = *P_union_X.begin();
    for (int v : P_union_X) {
        if (graph[v].size() > graph[u].size()) {
            u = v;
        }
    }
    
    unordered_set<int> P_minus_N_u;
    for (int v : P) {
        if (graph[u].find(v) == graph[u].end()) {
            P_minus_N_u.insert(v);
        }
    }
    
    for (int v : P_minus_N_u) {
        unordered_set<int> new_P, new_X;
        for (int w : P) if (graph[v].find(w) != graph[v].end()) new_P.insert(w);
        for (int w : X) if (graph[v].find(w) != graph[v].end()) new_X.insert(w);
        
        R.insert(v);
        BronKerboschPivot(graph, new_P, R, new_X);
        R.erase(v);
        
        P.erase(v);
        X.insert(v);
    }
}

void BronKerboschDegeneracy(const Graph& graph) {
    vector<int> ordering = degeneracyOrdering(graph);
    int n = graph.size();
    
    for (int i = 0; i < n; ++i) {
        int vi = ordering[i];
        unordered_set<int> P, X;
        
        for (int vj : graph[vi]) {
            if (find(ordering.begin() + i + 1, ordering.end(), vj) != ordering.end()) {
                P.insert(vj);
            } else if (find(ordering.begin(), ordering.begin() + i, vj) != ordering.begin() + i) {
                X.insert(vj);
            }
        }
        
        unordered_set<int> R = {vi};
        BronKerboschPivot(graph, P, R, X);
    }
}

int main() {
    clock_t start_time = clock(); // Start timer
    
    ifstream inputFile("input.txt");
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }

    // Node mapping structures
    unordered_map<int, int> nodeToIndex;
    Graph graph;
    int currentIndex = 0;

    // Read and process edges
    int reported_n, m;
    inputFile >> reported_n >> m; // Read and ignore reported node count

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

    BronKerboschDegeneracy(graph);

    clock_t end_time = clock(); // End timer
    double duration = (double)(end_time - start_time) / CLOCKS_PER_SEC; // Calculate duration in seconds
    
    cout << "1. Largest size of the clique: " << largest_clique_size << endl;
    cout << "2. Total number of maximal cliques: " << total_maximal_cliques << endl;
    cout << "3. Distribution of different size cliques:" << endl;

    for (const auto& pair : clique_size_distribution) {
        cout << "   Size " << pair.first << ": " << pair.second << endl;
    }

    cout << "Execution time: " << duration << " seconds" << endl;

    return 0;
}
