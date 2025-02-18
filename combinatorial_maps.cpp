#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

struct Dart {
    int from, to;
    int id;
    
    bool operator !=(const Dart& other) {
        return id != other.id;
    }
};

bool operator ==(const Dart& d1, const Dart& d2) {
    return d1.id == d2.id;
}

bool operator !=(const Dart& d1, const Dart& d2) {
    return d1.id != d2.id;
}

ostream& operator<<(ostream& out, const Dart& d) {
    out << "(" << d.id << ": " << d.from << "->" << d.to << ")";
    return out;
}

struct CombinatorialMap {
    int V, E, F; // Vertices, Edges, Faces
    vector<int> degree; // vertex degree
    // Two "copies" of each edge, one for each direction. 
    // Even loops are present twice.
    vector<Dart> edges;
    
    vector<int> _rev;
    vector<int> _next, _prev;
    
    
    pair<Dart,Dart> _AddEdge(int u, int v) {
        assert(0 <= u and u < V and 0 <= v and v < V);
        degree[u]++, degree[v]++;
        E++;
    
        Dart d1;
        d1.from = u, d1.to = v, d1.id = edges.size();
        edges.push_back(d1);
        Dart d2;
        d2.from = v, d2.to = u, d2.id = edges.size();
        edges.push_back(d2);
        
        _rev.push_back(d2.id);
        _rev.push_back(d1.id);
        
        _next.push_back(-1), _next.push_back(-1);
        _prev.push_back(-1), _prev.push_back(-1);
        
        return {d1, d2};
    }
    
    void _AssignNext(Dart a, Dart b) {
        assert(a.from == b.from);
        _next[a.id] = b.id;
        _prev[b.id] = a.id;
    }
    
    // Two vertices connected by an edge.
    //
    // ACHTUNG: We implicitly assume the graph does not consist of a 
    // single vertex and some loops.
    CombinatorialMap() {
        V = 2;
        degree.resize(V, 0);
        F = 1;
        E = 0;
        auto edge = _AddEdge(0, 1);
        _AssignNext(edge.first, edge.first);
        _AssignNext(edge.second, edge.second);
    }
    
    Dart Rev(Dart e) const { return edges[_rev[e.id]]; }
    // Next is counterclockwise
    Dart Next(Dart e) const { return edges[_next[e.id]]; }
    Dart Prev(Dart e) const { return edges[_prev[e.id]]; }
    // Phi traverses the boundary of a face clockwise.
    Dart Phi(Dart e) const { return edges[_next[_rev[e.id]]]; } 
    
    // Create a new vertex and link it (with a new edge) to an existing one.
    // 
    // The new edge is connected so that Next(bef) == new edge.
    pair<Dart,Dart> AddVertex(Dart bef) {
        V++;
        degree.push_back(0);
        
        const Dart& bef_next = Next(bef);
        pair<Dart,Dart> new_edge = _AddEdge(bef.from, V-1);
        
        _AssignNext(new_edge.second, new_edge.second);
        _AssignNext(bef, new_edge.first);
        _AssignNext(new_edge.first, bef_next);
        
        return new_edge;
    }
    
    // Link with an edge two existing vertices.
    // phi(bef) = new_edge 
    // phi(new_edge) = aft
    // 
    // This works (i.e., maintains the property that the combinatorial
    // maps represent a spherical embedding) only if bef and aft
    // belong to the boundary of the same face.
    // It may be bef = aft.
    pair<Dart,Dart> AddEdge(Dart bef, Dart aft) {
        // We check that aft = phi^k(bef) for some k >= 0.
        auto it = bef;
        while (it != aft) {
            it = Phi(it);
            assert(it != bef);
        }
        
        // We Create the new edge
        pair<Dart, Dart> new_edge = _AddEdge(bef.to, aft.from);
        
        // We update the next function
        const Dart& bef_next = Next(Rev(bef));
        const Dart& aft_prev = Prev(aft);
        
        _AssignNext(Rev(bef), new_edge.first);
        if (bef.to != aft.from) {
            _AssignNext(new_edge.first, bef_next);
            _AssignNext(aft_prev, new_edge.second);
        } else {
            _AssignNext(new_edge.first, new_edge.second);
        }
        _AssignNext(new_edge.second, aft);
        
        // We have added a new face in the process.
        F++;
        
        return new_edge;
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
        return {};
    
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

ULL Trace(const vector<vector<ULL>>& mat) {
    ULL trace = 0;
    for (int i = 0; i < ssize(mat); i++) trace += mat[i][i];
    return trace;
}

// Returns a square matrix 2E x 2E. Let M_A be the matrix.
// If f is a bijection from the darts of A to the darts of B that
// is a isomorphism of combinatorial maps, then M_B is obtained
// permuting the rows/columns of M_A according to f.
// In particular, tr(M_A) = tr(M_B).
vector<vector<ULL>> MatrixSignature(const CombinatorialMap& A) {
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
    return mat;
}

vector<CombinatorialMap> Evolve(const CombinatorialMap& A) {
    vector<CombinatorialMap> res;
    
    int leaf = -1;
    for (int v = 0; v < A.V; v++) if (A.degree[v] == 1) leaf = v;
    
    for (Dart d: A.edges) {
        if (leaf != -1 and d.from != leaf) continue;
        if (A.degree[d.from] >= 3) continue;
        CombinatorialMap A2 = A;
        A2.AddVertex(d);
        res.push_back(A2);
    }
    
    for (Dart d: A.edges) {
        if (leaf != -1 and d.to != leaf) continue;
        if (A.degree[d.to] >= 3) continue;

        Dart d2 = d;
        while (1) {
            bool to_consider = !(A.degree[d2.from] >= 3) 
                           and !(d.to == d2.from and A.degree[d.to] >= 2);
            
            if (to_consider) {
                CombinatorialMap A2 = A;
                A2.AddEdge(d, d2);
                res.push_back(A2);
            }
            
            d2 = A.Phi(d2);
            if (d2 == d) break;
        }
    }
    
    return res;
}


void RemoveDuplicates(vector<CombinatorialMap>& V) {
    map<ULL, vector<CombinatorialMap>> signature2list;
    for (const CombinatorialMap& A: V) {
        signature2list[Trace(MatrixSignature(A))].push_back(A);
    }
    V.clear();
    
    for (const auto& pp: signature2list) {
        const auto& vec = pp.second;
        vector<CombinatorialMap> unique_vec;
        for (const CombinatorialMap& A: vec) {
            bool to_add = true;
            for (const CombinatorialMap& B: unique_vec) {
                if (!FindIsomorphism(A, B).empty()) {
                    to_add = false;
                    break;
                }
            }
            if (to_add) {
                unique_vec.push_back(A);
                V.push_back(A);
            }
        }
    }
}

vector<CombinatorialMap> Evolve(const vector<CombinatorialMap>& V) {
    vector<CombinatorialMap> new_V;
    for (const CombinatorialMap& A: V) {
        const vector<CombinatorialMap>& vec = Evolve(A);
        for (const CombinatorialMap& B: vec) new_V.push_back(B);
    }
    RemoveDuplicates(new_V);
    return new_V;
}

int main() {
    vector<CombinatorialMap> vec(1);
    
    for (int i = 2; i <= 6; i++) {
        vec = Evolve(vec);
        cout << i << ": " << vec.size() << endl;
    }
}
