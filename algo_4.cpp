#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <functional>
#include <limits>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <chrono>
#include <unordered_set>
#include <list>
#include <numeric>
#include <stack>
using namespace std;

struct Graph {
    int n;                         
    vector<vector<int>> adj;        // adjacency
    vector<int> originalIds;        // mapping to original vertex ids
    
    Graph(int N = 0) : n(N), adj(N) {
        originalIds.resize(N);
        for (int i = 0; i < N; i++) originalIds[i] = i;
    }
    
    void addEdge(int u, int v){
        if(u==v) return;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    /* induced subgraph over keep[ ] */
    Graph induced(const vector<char>& keep) const {
        vector<int> id(n, -1); int cnt = 0;
        for(int i=0; i<n; ++i) if(keep[i]) id[i] = cnt++;
        
        Graph H(cnt);
        // Transfer original IDs
        for(int i=0; i<n; ++i) {
            if(keep[i]) {
                H.originalIds[id[i]] = originalIds[i];
            }
        }
        
        // Add edges
        for(int u=0; u<n; ++u) {
            if(keep[u]) {
                for(int v: adj[u]) {
                    if(keep[v] && u < v) {
                        H.addEdge(id[u], id[v]);
                    }
                }
            }
        }
        
        return H;
    }
    
    long long edgeCount() const {
        long long s = 0;
        for(auto &v:adj) s += v.size();
        return s/2;
    }
    
    // Print nodes sorted by original ids
    void printNodes() const {
        vector<int> nodes(n);
        for(int i=0; i<n; i++) {
            nodes[i] = originalIds[i];
        }
        sort(nodes.begin(), nodes.end());
        for(int i=0; i<n; i++) {
            cout << nodes[i];
            if(i < n-1) cout << " ";
        }
        cout << endl;
    }
};

/* ---------- h-clique counting utility functions ---------- */
// Utility function to count triangles in a graph
long long countTriangles(const Graph& G) {
    long long count = 0;
    for(int u = 0; u < G.n; ++u) {
        unordered_set<int> neighbors(G.adj[u].begin(), G.adj[u].end());
        for(int v : G.adj[u]) {
            if(v > u) { // Count each triangle once
                for(int w : G.adj[v]) {
                    if(w > v && neighbors.count(w)) { // Ensure u < v < w to count each triangle once
                        count++;
                    }
                }
            }
        }
    }
    return count;
}

// Function to count h-cliques (generalized)
long long countHCliques(const Graph& G, int h) {
    if(h == 1) {
        return G.n;  // For h=1, cliques are just vertices
    } else if(h == 2) {
        return G.edgeCount(); // For h=2, cliques are edges
    } else if(h == 3) {
        return countTriangles(G); // For h=3, cliques are triangles
    } else {
        // For h > 3, implement a more general clique counting algorithm
        // This is a simplified placeholder for demonstration
        // A more efficient implementation would be needed for larger h values
        return 0; // Return 0 for unsupported h values
    }
}

/* ---------- Dinic max-flow (double capacities) ---------- */
struct Dinic {
    struct Edge{ int to, rev; double cap; };
    int N,s,t;  vector<vector<Edge>> G; vector<int> level, it;
    Dinic(int n):N(n),G(n),level(n),it(n){}
    void addEdge(int u,int v,double c){
        Edge a{v,(int)G[v].size(),c};
        Edge b{u,(int)G[u].size(),0};
        G[u].push_back(a); G[v].push_back(b);
    }
    bool bfs(){
        fill(level.begin(), level.end(), -1);
        queue<int> q; level[s]=0; q.push(s);
        while(!q.empty()){
            int v=q.front(); q.pop();
            for(auto &e:G[v]) if(e.cap>1e-9 && level[e.to]<0){
                level[e.to]=level[v]+1; q.push(e.to);
            }
        }
        return level[t]>=0;
    }
    double dfs(int v,double f){
        if(v==t || f<1e-9) return f;
        for(int &i=it[v]; i<(int)G[v].size(); ++i){
            Edge &e=G[v][i];
            if(e.cap>1e-9 && level[v]<level[e.to]){
                double ret = dfs(e.to, min(f,e.cap));
                if(ret>1e-9){
                    e.cap -= ret;
                    G[e.to][e.rev].cap += ret;
                    return ret;
                }
            }
        }
        return 0;
    }
    double maxflow(int S,int T){
        s=S; t=T; double flow=0;
        while(bfs()){
            fill(it.begin(), it.end(),0);
            while(double f=dfs(s,1e100)) flow+=f;
        }
        return flow;
    }
    /* vertices reachable from s in residual graph */
    vector<char> minCut(){
        vector<char> vis(N,0);
        queue<int> q; q.push(s); vis[s]=1;
        while(!q.empty()){
            int v=q.front(); q.pop();
            for(auto &e:G[v]) if(e.cap>1e-9 && !vis[e.to]){
                vis[e.to]=1; q.push(e.to);
            }
        }
        return vis;
    }
};

/* ---------- core decomposition (for h=2, edge density) ----------*/
// For edge density, we can use degree-based core decomposition
struct CoreInfo {
    vector<int> core;          // core number of each vertex
    int kmax;                  // max core number
    double bestDensity;        // ρ''  (highest residual density)
    int bestK;                 // k'' = ceil(ρ'')
};

CoreInfo coreDecompose(const Graph& G, int h) {
    int n = G.n;
    vector<int> deg(n);
    for(int i = 0; i < n; ++i) deg[i] = G.adj[i].size();
    int maxDeg = 0;
    for(int i = 0; i < n; ++i) maxDeg = max(maxDeg, deg[i]);
    
    vector<vector<int>> bucket(maxDeg + 1);
    for(int i = 0; i < n; ++i) bucket[deg[i]].push_back(i);

    vector<int> core(n, 0);
    vector<char> alive(n, 1);
    int processed = 0;
    
    vector<int> vOrder; // Order of vertex removal
    
    for(int curd = 0; curd <= maxDeg; ++curd) {
        while(!bucket[curd].empty()) {
            int v = bucket[curd].back(); bucket[curd].pop_back();
            if(!alive[v] || deg[v] != curd) continue;
            
            alive[v] = 0; 
            core[v] = curd; 
            vOrder.push_back(v);
            processed++;
            
            for(int nb: G.adj[v]) if(alive[nb]) {
                int d = deg[nb];
                deg[nb]--;
                // Move from bucket[d] to bucket[d-1]
                bucket[d-1].push_back(nb);
            }
        }
    }
    
    // Calculate densities in reverse order of vertex removal
    int kmax = 0;
    for(int i = 0; i < n; i++) {
        kmax = max(kmax, core[i]);
    }
    
    // Calculate best residual density
    double bestRho = 0;
    int bestK = 0;
    
    // Create a graph that initially has all vertices
    vector<char> keepAll(n, 1);
    Graph currentG = G;
    long long edges = currentG.edgeCount();
    
    // Remove vertices in order and recalculate density
    for(int i = 0; i < vOrder.size(); i++) {
        int remaining = n - i;
        if(remaining > 0) {
            double rho = (double)edges / remaining;
            if(rho > bestRho) {
                bestRho = rho;
                bestK = ceil(rho);
            }
        }
        
        // Remove the next vertex
        if(i < vOrder.size()) {
            int v = vOrder[i];
            keepAll[v] = 0;
            
            // Update edge count
            for(int nb: G.adj[v]) {
                if(keepAll[nb]) {
                    edges--;
                }
            }
        }
    }
    
    return {core, kmax, bestRho, bestK};
}

/* ---------- connected components utility ---------- */
vector<vector<int>> connectedComponents(const Graph& G) {
    int n = G.n; 
    vector<char> vis(n, 0); 
    vector<vector<int>> comps;
    
    for(int i = 0; i < n; ++i) {
        if(!vis[i]) {
            vector<int> comp; 
            stack<int> st; 
            st.push(i); 
            vis[i] = 1;
            
            while(!st.empty()) {
                int v = st.top(); 
                st.pop(); 
                comp.push_back(v);
                
                for(int nb : G.adj[v]) {
                    if(!vis[nb]) {
                        vis[nb] = 1; 
                        st.push(nb);
                    }
                }
            }
            
            comps.push_back(std::move(comp));
        }
    }
    
    return comps;
}

/* ---------- build flow network for edge density (h=2) ---------- */
Dinic buildFlowNetwork(const Graph& C, double alpha) {
    int n = C.n;
    long long m = C.edgeCount();
    int N = 2 + n;               // s=0, t=1, vertices offset +2
    Dinic D(N);
    
    // Source to each vertex with capacity m
    for(int v = 0; v < n; ++v) {
        D.addEdge(0, v+2, m);  // s→v
    }
    
    // Each vertex to sink with capacity m + 2*alpha - deg(v)
    for(int v = 0; v < n; ++v) {
        double cap = m + 2*alpha - C.adj[v].size();
        D.addEdge(v+2, 1, cap);  // v→t
    }
    
    // Bidirectional edges between vertices with capacity 1
    for(int u = 0; u < n; ++u) {
        for(int v : C.adj[u]) {
            if(u < v) {
                D.addEdge(u+2, v+2, 1);
                D.addEdge(v+2, u+2, 1);
            }
        }
    }
    
    return D;
}

/* ---------- binary search on one connected component ---------- */
pair<double, vector<char>> densestOnComponent(const Graph& C, double l, double u, int h) {
    if(C.n <= 1) return {0, vector<char>(C.n, 1)};
    
    double eps = 1.0 / (C.n * (C.n-1));
    vector<char> bestKeep(C.n, 0);
    double bestAlpha = 0;
    
    while(u - l >= eps) {
        double alpha = (l + u) / 2.0;
        
        // Build the flow network for this alpha
        Dinic D = buildFlowNetwork(C, alpha);
        D.maxflow(0, 1);
        auto reach = D.minCut();  // reach[v]==1 → in S
        
        // Create the subgraph defined by S\{s}
        vector<char> keep(C.n, 0);
        int kept = 0;
        for(int i = 0; i < C.n; ++i) {
            // +2 because the vertices in the flow network start at index 2
            if(reach[i+2]) {
                keep[i] = 1;
                kept++;
            }
        }
        
        if(kept == 0) {     // S = {s}
            u = alpha;   // infeasible
        } else {
            l = alpha;   // feasible, remember best
            bestKeep = keep;
            bestAlpha = alpha;
        }
    }
    
    return {bestAlpha, bestKeep};
}

/* ---------- CoreExact (Algorithm 4 for h=2 edge density) ---------- */
Graph coreExact(const Graph& G, int h) {
    // Skip if graph is empty
    if(G.n == 0) return G;
    
    /* step-1: core decomposition & residual densest lower-bound */
    CoreInfo info = coreDecompose(G, h);
    int kmax = info.kmax;
    int kpp = info.bestK;           // k''         (line 2)
    double l = info.bestDensity;    // ρ'' (lower) (line 3)
    double u = kmax;                // upper       (line 3)
    
    if(kpp == 0) kpp = 1;  // Ensure kpp is at least 1
    
    /* build (k'',Ψ)-core (here k''-core) */
    vector<char> keep(G.n, 0);
    for(int v = 0; v < G.n; ++v) {
        if(info.core[v] >= kpp) {
            keep[v] = 1;
        }
    }
    
    Graph Gk = G.induced(keep);  // line 4
    
    // If k''-core is empty, return original graph
    if(Gk.n == 0) return G;
    
    auto compsIdx = connectedComponents(Gk);
    
    Graph best(0);  // Empty graph initially
    double bestDensity = 0;
    
    /* iterate components (line 5) */
    for(auto &idxList : compsIdx) {
        vector<char> keep2(Gk.n, 0);
        for(int v : idxList) keep2[v] = 1;
        Graph C = Gk.induced(keep2);
        
        // Line 6: if l > k'', refine C as C ∩ ([l], Ψ)-core
        if(l > kpp) {
            vector<char> refinedKeep(C.n, 0);
            int ceilL = ceil(l);
            for(int v = 0; v < C.n; ++v) {
                int originalV = C.originalIds[v];
                int originalIndex = -1;
                for(int i = 0; i < G.n; ++i) {
                    if(G.originalIds[i] == originalV) {
                        originalIndex = i;
                        break;
                    }
                }
                if(originalIndex != -1 && info.core[originalIndex] >= ceilL) {
                    refinedKeep[v] = 1;
                }
            }
            C = C.induced(refinedKeep);
            if(C.n == 0) continue;  // Skip empty component
        }
        
        
        auto result = densestOnComponent(C, l, u, h);
        double alpha = result.first;
        auto& subgraphKeep = result.second;
        
        
        // Create subgraph U = S\{s} (line 19)
        Graph U = C.induced(subgraphKeep);
        
        if(U.n > 0) {
            // For h=2, density is |E|/|V|
            double rhoU = (double)U.edgeCount() / U.n;
            
            // Line 20: if ρ(G[U], Ψ) > ρ(D, Ψ) then D ← G[U]
            if(rhoU > bestDensity) {
                bestDensity = rhoU;
                best = U;
            }
        }
    }
    
    // If we didn't find a better subgraph, return original graph
    if(best.n == 0) return G;
    
    return best;
}

/* ---------- Greedy peeling algorithm for comparison ---------- */
// This is often used as a fast approximation or to test correctness
Graph greedyDensest(const Graph& G) {
    int n = G.n;
    vector<int> deg(n);
    for(int i = 0; i < n; ++i) deg[i] = G.adj[i].size();
    
    vector<char> removed(n, 0);
    vector<pair<double, vector<char>>> candidates; // (density, keep)
    
    long long edges = G.edgeCount();
    int remaining = n;
    
    // Start with all vertices
    vector<char> keep(n, 1);
    double density = (double)edges / n;
    candidates.push_back({density, keep});
    
    // Repeatedly remove the vertex with minimum degree
    for(int iter = 0; iter < n - 1; iter++) {
        // Find minimum degree vertex
        int minDeg = n;
        int minVertex = -1;
        for(int v = 0; v < n; v++) {
            if(!removed[v] && deg[v] < minDeg) {
                minDeg = deg[v];
                minVertex = v;
            }
        }
        
        if(minVertex == -1) break;
        
        // Remove the vertex
        removed[minVertex] = 1;
        remaining--;
        
        // Update edge count and degrees
        for(int nb : G.adj[minVertex]) {
            if(!removed[nb]) {
                edges--;
                deg[nb]--;
            }
        }
        
        // Calculate new density
        if(remaining > 0) {
            double newDensity = (double)edges / remaining;
            
            // Update keep vector and add to candidates
            keep[minVertex] = 0;
            candidates.push_back({newDensity, keep});
        }
    }
    
    // Find the subgraph with maximum density
    double maxDensity = 0;
    vector<char> bestKeep;
    
    for (auto& candidate : candidates) {
        double d = candidate.first;
        auto& k = candidate.second;
        if (d > maxDensity) {
            maxDensity = d;
            bestKeep = k;
        }
    }
    
    
    return G.induced(bestKeep);
}

/* ---------- I/O helpers ---------- */
pair<int, Graph> readEdgeList(istream& in) {
    int h, numVertices, numEdges;
    in >> h;  // Read h parameter
    in >> numVertices >> numEdges;
    
    Graph G(numVertices);
    
    for(int i = 0; i < numEdges; i++) {
        int u, v;
        in >> u >> v;
        G.addEdge(u, v);
    }
    
    return {h, G};
}

/* ---------- helper to generate all h-cliques ---------- */
void enumerateCliques(const Graph& G, int h, vector<vector<int>>& cliques, vector<int> current = {}, int start = 0) {
    if (current.size() == h) {
        cliques.push_back(current);
        return;
    }
    for (int i = start; i < G.n; ++i) {
        bool valid = true;
        for (int v : current) {
            if (find(G.adj[i].begin(), G.adj[i].end(), v) == G.adj[i].end()) {
                valid = false;
                break;
            }
        }
        if (valid) {
            current.push_back(i);
            enumerateCliques(G, h, cliques, current, i + 1);
            current.pop_back();
        }
    }
}

/* ---------- main ---------- */
int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int h;
    Graph G;
    
    if(argc > 1) {
        ifstream fin(argv[1]);
        if(!fin) { cerr << "cannot open " << argv[1] << "\n"; return 1; }
        tie(h, G) = readEdgeList(fin);
    } else {
        tie(h, G) = readEdgeList(cin);
    }

    auto start = chrono::high_resolution_clock::now();

    // STEP 1: Enumerate all h-cliques of G first
    vector<vector<int>> cliquesH;
    if (h >= 2 && h <= 5) {  // You can restrict to h <= 5 to avoid explosion
        enumerateCliques(G, h, cliquesH);
    } else {
        cerr << "Unsupported h=" << h << ", counting skipped.\n";
    }

    // STEP 2: Find densest subgraph
    Graph D1 = coreExact(G, h);
    Graph D2 = greedyDensest(G);
    
    double density1 = D1.n > 0 ? (double)D1.edgeCount() / D1.n : 0;
    double density2 = D2.n > 0 ? (double)D2.edgeCount() / D2.n : 0;
    
    Graph D = (density1 >= density2) ? D1 : D2;

    // STEP 3: Calculate h-clique count inside D
    unordered_set<int> selectedVertices;
    for (int i = 0; i < D.n; ++i) {
        selectedVertices.insert(D.originalIds[i]);
    }

    int cliqueCountInD = 0;
    for (const auto& clique : cliquesH) {
        bool allInside = true;
        for (int v : clique) {
            if (selectedVertices.find(v) == selectedVertices.end()) {
                allInside = false;
                break;
            }
        }
        if (allInside)
            cliqueCountInD++;
    }

    double rhoOpt = selectedVertices.empty() ? 0.0 : (double)cliqueCountInD / selectedVertices.size();

    auto end = chrono::high_resolution_clock::now();
    double elapsed_sec = chrono::duration<double>(end - start).count();

    // STEP 4: Output
    cerr << "h=" << h << " V=" << G.n << "  E=" << G.edgeCount() << "\n";
    cout << "Densest Subgraph Nodes:" << endl;
    D.printNodes();
    cout << "Nodes in densest subgraph: " << selectedVertices.size() << endl;
    cout << "h-Cliques in densest subgraph: " << cliqueCountInD << endl;
    cout << "h-Clique Density: " << rhoOpt << endl;
    cout << "Total execution time: " << elapsed_sec << " seconds" << endl;
    
    return 0;
}