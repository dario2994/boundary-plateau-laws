#ifndef boundary_taylor_graph
#define boundary_taylor_graph
#include "utils.cpp"

// Vertices are 0-indexed.
// The order of edges matters; since this will be used to represent an
// embedding of a planar graph.
struct Graph {
    int n;
    vector<vector<int>> edges;

    Graph(): n(0), edges({}) {}

    void add_vertex() {
        edges.push_back({});
        n++;
    }

    // Adds edge u->v, so that it's the first edge for u and the last for v.
    void add_edge(int u, int v) {
        edges[u].insert(edges[u].begin(), v);
        edges[v].push_back(u);
    }

    int degree(int v) const { return edges[v].size(); }

    vector<int> dfs_order;
    vector<bool> visited;
    void _dfs(int v, int father) {
        visited[v] = true;
        int i0 = 0;
        while (edges[v][i0] != father) i0++;
        for (int k = 0; k < (int)edges[v].size(); k++) {
            int u = get(edges[v], i0 + k);
            dfs_order.push_back(v);
            dfs_order.push_back(u);
            if (visited[u]) continue;
            _dfs(u, v);
        }
    }

    // Return a vector containing:
    //  x_0, y_0, x_1, y_1, x_2, y_2, ...
    // where (x_i, y_i) are the edges visited by the dfs starting from the edge
    // (u, v).
    // The edges are considered even if that's a dead-end for the dfs (i.e.,
    // if the final node was already visited).
    vector<int> ordered_dfs(int u, int v) {
        visited.resize(n, false);
        fill(visited.begin(), visited.end(), false);
        dfs_order.clear();
        _dfs(u, v);
        return dfs_order;
    }

    void pretty_print() {
        for (int v = 0; v < n; v++) {
            cout << v << ": ";
            for (int u: edges[v]) cout << u << " ";
            cout << endl;
        }
    }
};

#endif
