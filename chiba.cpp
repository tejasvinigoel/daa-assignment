#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <map>
#include <climits>
#include <unordered_set>
#include <numeric>

using namespace std;

using Graph = vector<set<int>>;

int largest_clique_size = 0;
int total_maximal_cliques = 0;
map<int, int> clique_size_distribution;
set<set<int>> unique_cliques; // Track unique cliques

vector<int> T, S;

void print_clique(const set<int> &clique, const vector<int> &original_order)
{
    // Removed printing of intermediate steps
    // cout << "Clique: {";
    // for (auto it = clique.begin(); it != clique.end(); ++it)
    // {
    //     if (it != clique.begin())
    //         cout << ", ";
    //     cout << original_order[*it];
    // }
    // cout << "}\n";
}

void UPDATE(int i, const set<int> &C, const Graph &graph, const vector<int> &original_order)
{
    if (i == graph.size())
    {
        if (C.size() >= 2 && unique_cliques.find(C) == unique_cliques.end())
        {
            // Maximality check against all existing cliques
            bool is_maximal = true;
            for (const auto &existing : unique_cliques)
            {
                if (includes(existing.begin(), existing.end(), C.begin(), C.end()))
                {
                    is_maximal = false;
                    break;
                }
            }

            if (is_maximal)
            {
                // Remove all subsets of this new clique
                vector<set<int>> to_remove;
                for (const auto &existing : unique_cliques)
                {
                    if (includes(C.begin(), C.end(), existing.begin(), existing.end()))
                    {
                        to_remove.push_back(existing);
                    }
                }
                for (const auto &s : to_remove)
                {
                    unique_cliques.erase(s);
                }

                unique_cliques.insert(C);
                int size = C.size();
                largest_clique_size = max(largest_clique_size, size);
                clique_size_distribution[size]++;
                total_maximal_cliques = unique_cliques.size();
                // Removed intermediate printing of clique:
                // print_clique(C, original_order);
            }
        }
        return;
    }

    const int current = i;
    const auto &N_i = graph[current];

    // Case 1: C is subset of N(i)
    if (all_of(C.begin(), C.end(), [&](int v)
               { return N_i.count(v); }))
    {
        UPDATE(i + 1, C, graph, original_order);
    }

    // Prepare C ∩ N(i) and C - N(i)
    set<int> C_intersect_Ni, C_diff_Ni;
    for (int v : C)
    {
        if (N_i.count(v))
            C_intersect_Ni.insert(v);
        else
            C_diff_Ni.insert(v);
    }

    // Track modifications for T and S
    vector<pair<int, int>> T_mods, S_mods;

    // Update T[y] = |N(y) ∩ C ∩ N(i)|
    for (int x : C_intersect_Ni)
    {
        for (int y : graph[x])
        {
            if (y != current && !C.count(y))
            {
                T_mods.emplace_back(y, T[y]);
                T[y]++;
            }
        }
    }

    // Update S[y] = |N(y) ∩ (C - N(i))|
    for (int x : C_diff_Ni)
    {
        for (int y : graph[x])
        {
            if (!C.count(y))
            {
                S_mods.emplace_back(y, S[y]);
                S[y]++;
            }
        }
    }

    // Strict maximality test
    bool is_maximal = true;
    for (int y : N_i)
    {
        if (y < current && !C.count(y) && T[y] == (int)C_intersect_Ni.size())
        {
            is_maximal = false;
            break;
        }
    }

    // Lexicographic test with exact bounds
    if (is_maximal)
    {
        vector<int> C_diff_sorted(C_diff_Ni.begin(), C_diff_Ni.end());
        sort(C_diff_sorted.begin(), C_diff_sorted.end());
        const int p = C_diff_sorted.size();

        for (int k = 0; k < p && is_maximal; ++k)
        {
            int jk = C_diff_sorted[k];
            for (int y : graph[jk])
            {
                if (y < current && !C.count(y) && T[y] == (int)C_intersect_Ni.size())
                {
                    if (y >= jk)
                    {
                        S[y]--;
                    }
                    else
                    {
                        int first_j = -1;
                        for (int j : C_diff_sorted)
                        {
                            if (y < j)
                            {
                                first_j = j;
                                break;
                            }
                        }
                        if (first_j == jk)
                        {
                            int j_prev = k > 0 ? C_diff_sorted[k - 1] : -1;
                            if (S[y] + k >= p && y >= j_prev)
                            {
                                is_maximal = false;
                                break;
                            }
                        }
                    }
                }
            }
        }

        if (is_maximal && C_intersect_Ni.empty())
        {
            int jp = C_diff_sorted.empty() ? -1 : C_diff_sorted.back();
            for (int y : C)
            {
                if (y < current && T[y] == 0 && S[y] == 0)
                {
                    if (jp < y || jp < current - 1)
                    {
                        is_maximal = false;
                        break;
                    }
                }
            }
        }
    }

    // Restore T and S
    for (const auto &element : T_mods)
        T[element.first] = element.second;
    for (const auto &element : S_mods)
        S[element.first] = element.second;

    // Case 2: Expand with current vertex if maximal
    if (is_maximal)
    {
        set<int> new_C = C_intersect_Ni;
        new_C.insert(current);
        UPDATE(i + 1, new_C, graph, original_order);
    }

    // Process without current vertex
    UPDATE(i + 1, C, graph, original_order);
}

void CLIQUE(const Graph &graph, const vector<int> &original_order)
{
    const int n = graph.size();
    T.assign(n, 0);
    S.assign(n, 0);
    UPDATE(0, {}, graph, original_order);
}

int main()
{
    ifstream inputFile("input.txt");
    // ... (rest of main remains the same until after CLIQUE call)

    if (!inputFile.is_open())
    {
        cerr << "Error: Could not open input file." << endl;
        return 1;
    }

    int n, m;
    inputFile >> n >> m;

    vector<pair<int, int>> edges;
    int min_id = INT_MAX;
    int max_id = 0;
    for (int i = 0; i < m; i++)
    {
        int u, v;
        inputFile >> u >> v;
        edges.push_back(make_pair(u, v));
        min_id = min(min_id, min(u, v));
        max_id = max(max_id, max(u, v));
    }
    inputFile.close();

    const bool convert = (min_id == 1);
    const int graph_size = max_id + (convert ? 0 : 1);

    // Build original graph
    Graph graph(graph_size);
    for (const auto &edge : edges)
    {
        int u = edge.first;
        int v = edge.second;
        if (convert)
        {
            u--;
            v--;
        }
        graph[u].insert(v);
        graph[v].insert(u);
    }

    // Sort vertices by degree
    vector<int> vertices(graph_size);
    iota(vertices.begin(), vertices.end(), 0);
    sort(vertices.begin(), vertices.end(), [&](int a, int b)
         { return graph[a].size() < graph[b].size(); });

    // Create new index mapping
    vector<int> new_index(graph_size);
    for (int i = 0; i < graph_size; ++i)
    {
        new_index[vertices[i]] = i;
    }

    // Build sorted adjacency list
    Graph sorted_graph(graph_size);
    for (int u = 0; u < graph_size; ++u)
    {
        for (int v : graph[u])
        {
            sorted_graph[new_index[u]].insert(new_index[v]);
        }
    }

    // Removed processing print: cout << "Processing cliques...\n";
    CLIQUE(sorted_graph, vertices);

    // After CLIQUE call, update distribution from unique_cliques
    clique_size_distribution.clear();
    for (const auto &clique : unique_cliques)
    {
        clique_size_distribution[clique.size()]++;
    }
    largest_clique_size = clique_size_distribution.rbegin()->first;
    total_maximal_cliques = unique_cliques.size();

    cout << "\nFinal Results:\n";
    cout << "1. Largest size of the clique: " << largest_clique_size << "\n"
         << "2. Total number of maximal cliques: " << total_maximal_cliques << "\n"
         << "3. Distribution of different size cliques:\n";
    for (const auto &pair : clique_size_distribution)
    {
        cout << "   Size " << pair.first << ": " << pair.second << "\n";
    }

    return 0;
}