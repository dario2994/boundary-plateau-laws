#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;


struct CombinatorialMap {
    struct Dart {
        int from, to;
        int id;
    };

    int V, E, F; // Vertices, Edges, Faces
    vector<Dart> edges;
    
    vector<int> _rev;
    Dart Rev(Dart e) const { return edges[_rev[e.id]]; }
    
    pair<Dart,Dart> _AddEdge(int u, int v) {
        assert(0 <= u and u < V and 0 <= v and v < V);
    
        Dart d1;
        d1.from = u, d1.to = v, d1.id = edges.size();
        edges.push_back(d1);
        Dart d2;
        d2.from = v, d2.to = u, d2.id = edges.size();
        edges.push_back(d2);
        
        _rev[d1.id] = d2.id, _rev[d2.id] = d1.id;
        
        return {d1, d2};
    }
    
    vector<int> _next, _prev;
    Dart Next(Dart e) const { return edges[_next[e.id]]; }
    Dart Prev(Dart e) const { return edges[_prev[e.id]]; }
    Dart Phi(Dart e) const { return edges[_next[_rev[e.id]]; } 
    
    void _AssignNext(Dart a, Dart b) {
        assert(a.from == b.from);
        _next[a.id] = b.id;
        _prev[b.id] = a.id;
    }
    
    // Initialize it as two faces with boundary given by _V edges.
    // _V = 1 is possible and means that the boundary of both faces
    // is a loop.
    CombinatorialMap(int _V) {
        assert(1 <= _V);
        V = _V;
        E = V;
        F = 2;
        vector<pair<Dart,Dart>> cycle;
        for (int i = 0; i < V; i++) {
            cycle.push_back(_AddEdge(i, (i+1) % V));
        }
        for (int i = 0; i < V; i++) {
            auto e1 = cycle[i];
            auto e2 = cycle[(i+1)% V];
            _AssignNext(e1.second, e2.first);
        }
    }
    
    // Attach a path consisting of len edges (if len > 0, then len-1
    // vertices are added to the graph). This creates a new face.
    // The first edge must be Next(Rev(bef)) and the last edge must be Rev(Prev(aft)).
    // 
    // For this function to keep the assumption that the combinatorial
    // map is induced by a spherical embedding, bef and aft shall belong to
    // the boundary of an existing face. This condition is never checked. 
    // Maybe we shall add a check?
    //
    //
    // ACHTUNG: It might be that bef.to == aft.from.
    //          This special case requires additional care.
    void AddPath(int len, const Dart& bef, const Dart& aft) {
        E += len;
        F += 1;
        
        const Dart& bef_next = Next(Rev(bef));
        const Dart& aft_prev = Prev(aft);
        
        vector<int> path_v(len+1);
        path_v[0] = bef.to;
        path_v[len] = aft.from;
        for (int i = 1; i < len; i++) {
            path_v[i] = V;
            V++;
        }
        
        vector<pair<Dart, Dart>> path_e;
        for (int i = 0; i < len; i++) 
            path_e.push_back(_AddEdge(path_v[i], path_v[i+1]));
        
        _AssignNext(Rev(bef), path_e[0].first);
        _AssignNext(path_e[0].first, bef_next);
        _AssignNext(aft_prev, path_e[len-1].second);
        _AssignNext(path_e[len-1].second, Rev(aft));
        for (int i = 0; i < len-1; i++)
            _AssignNext(path_e[i].second, path_e[i+1].first);
        
        // Taking care of bef.to == aft.from    
        if (bef.to == aft.from) {
            assert(bef_next.id == aft.id);
            _AssignNext(path_e[0].first, path_e[len-1].second);
            
        }
    }
};

void _dfs(const CombinatorialMap& A, int e, vector<bool>& visited, 
          vector<int>& ord) {
    if (visited[e]) return;
    visited[e] = true;
    ord.push_back(e);
    _dfs(A, A._rev[e], visited, ord);
    _dfs(A, A._next[e], visited, ord);
}

vector<int> OrderedVisit(const CombinatorialMap& A, int e) {
    assert(0 <= e and e < 2*A.E);
    vector<bool> visited(2*A.E, false);
    vector<int> ord;
    _dfs(A, e, visited, ord);
    return ord;
}

// Check whether the map iso between the darts of A and B induces
// an isomorphism of the combinatorial map.
// It is sufficient to check that it commutes with Rev and Next.
bool CheckIsomorphism(const CombinatorialMap& A,
                      const CombinatorialMap& B,
                      const vector<int>& iso) {
    int E = A.E;
    assert(A.E == E and B.E == E and ssize(iso) == 2*E);
    for (int i = 0; i < 2*E; i++) {
        if (B._rev[iso[i]] != iso[A._rev[i]]) return false;
        if (B._next[iso[i]] != iso[A._next[i]]) return false;
    }
    return true;
}


// Returns, if it exists, a bijection between the darts that is an
// isomorphism of the combinatorial maps.
// The ids of the darts are used for a concise representation.
vector<int> FindIsomorphism(const CombinatorialMap& A,
                            const CombinatorialMap& B) {
    if (A.V != B.V or A.E != B.E or A.F != B.F)
        return {{}, {}};
    
    int E = A.E;
    vector<int> ordA = OrderedVisit(A, 0);
    for (int e = 0; e < 2*E; e++) {
        vector<int> ordB = OrderedVisit(B, e);
        vector<int> iso(2*E);
        for (int i = 0; i < 2*E; i++) iso[ordA[i]] = ordB[i];
        
        if (CheckIsomorphism(A, B, iso))
            return iso;
    }
    return {};
}

vector<vector<ULL>> MatrixPow(const vector<vector<ULL>>& A,
                              const vector<vector<ULL>>& B) {
    int n = A.size();
    assert(ssize(A[0]) == n and ssize(B) == n and ssize(B[0]) == n);
    vector<vector<ULL>> C(n, vector<ULL>(n, 0));
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++) C[i][j] += A[i][k] * B[k][j];
    return C;
}

ULL IsomorphismSignature(const CombinatorialMap& A) {
    static const ULL p_stay = 424'142'591;
    static const ULL p_rev = 984'689'391;
    static const ULL p_next = 123'456'789;
    
    vector<vector<ULL>> mat(2*A.E, vector<ULL>(2*A.E, 0));
    for (int i = 0; i < 2*A.E; i++) {
        mat[i][i] += p_stay;
        mat[i][A._rev[i]] += p_rev;
        mat[i][A._next[i]] += p_next;
    }
    for (int _ = 0; _ < 15; _++) mat = MatrixPow(mat, mat);
    ULL trace = 0;
    for (int i = 0; i < 2*A.E; i++) trace += mat[i][i];
    return trace;
}

int main() {
    cout << "Ciao" << endl;
    
    
    //~ PlanarEmbedding X;
    //~ AreIsomorphic(X, 0, 0, X, 0, 0); 
}
