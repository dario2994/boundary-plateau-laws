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
    // Phi traverses the boundary of a (bounded) face clockwise.
    int _phi(int id) const { return _next[_rev[id]]; }
    Dart Phi(Dart e) const { return edges[_phi(e.id)]; } 
    
    
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

struct PlanarMap:CombinatorialMap {
    // An arbitrary dart on the boundary of the unbounded face.
    Dart boundary;
    
    PlanarMap(): CombinatorialMap() {
        boundary = edges[0];
    }
    
    // Returns a list of all darts on the boundary of the unbounded face.
    vector<Dart> BoundaryEdges() const {
        vector<Dart> vec = {boundary};
        while (Phi(vec.back()) != boundary)
            vec.push_back(Phi(vec.back()));
        
        return vec;
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

// Returns the list of the darts encountered when visiting A starting
// from e. The visit is a dfs on the graph with vertices given by the
// darts and edges:
//  e -> rev[e]
//  e -> next[e]
// First, the visit traverses rev and then next.
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

// Returns, if it exists, an isomorphism between A and B that maps dart_A
// into dart_B.
// The ids of the darts are used for a concise representation.
vector<int> FindIsomorphism(const CombinatorialMap& A,
                            const CombinatorialMap& B,
                            Dart dart_A,
                            Dart dart_B) {
    int E = A.E;
    
    vector<int> ordA = OrderedVisit(A, dart_A.id);
    vector<int> ordB = OrderedVisit(B, dart_B.id);
    vector<int> iso(2*E);
    for (int i = 0; i < 2*E; i++) iso[ordA[i]] = ordB[i];
    
    if (CheckIsomorphism(A, B, iso))
        return iso;
    return {};
}


// Returns, if it exists, a bijection between the darts that is an
// isomorphism of the combinatorial maps.
// The ids of the darts are used for a concise representation.
vector<int> FindIsomorphism(const CombinatorialMap& A,
                            const CombinatorialMap& B) {
    if (A.V != B.V or A.E != B.E or A.F != B.F)
        return {};

    for (Dart d: B.edges) {
        const vector<int>& iso = FindIsomorphism(A, B, A.edges[0], d);
        if (!iso.empty())
            return iso;
    }
    return {};
}

vector<int> FindIsomorphism(const PlanarMap& A,
                            const PlanarMap& B) {
    if (A.V != B.V or A.E != B.E or A.F != B.F)
        return {};
    vector<Dart> boundary = B.BoundaryEdges();
    for (Dart d: boundary) {
        const vector<int>& iso = FindIsomorphism(A, B, A.edges[0], d);
        if (!iso.empty())
            return iso;
    }
    return {};
}

vector<vector<ULL>> MatrixMultiplication(const vector<vector<ULL>>& A,
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
    for (int _ = 0; _ < 15; _++) mat = MatrixMultiplication(mat, mat);
    return mat;
}

ULL IsomorphismSignature(const CombinatorialMap& A) {
    return Trace(MatrixSignature(A));
}

ULL IsomorphismSignature(const PlanarMap& A) {
    const auto& mat = MatrixSignature(A);
    ULL res = Trace(mat);
    vector<Dart> boundary = A.BoundaryEdges();
    for (Dart d: boundary) res += 131'657'115 * mat[d.id][d.id];
    return res;
}


// Checks whether A can eventually be transformed (by attaching edges
// on the "outside") into a interesting geodesic net.
bool Validate(const PlanarMap& A) {
    if (A.E > 32 or A.F > 12)
        return false;
    
    set<int> boundary_vertices;
    const vector<Dart>& boundary = A.BoundaryEdges();
    for (const Dart& d: boundary)
        boundary_vertices.insert(d.from);
    
    for (int v = 0; v < A.V; v++) {
        if (!boundary_vertices.count(v) and A.degree[v] != 3) 
            return false;

        if (A.degree[v] > 3)
            return false;
    }
    
    multiset<int> faces_sides;
    vector<bool> visited(2*A.E);
    for (int i = 0; i < 2*A.E; i++) {
        if (visited[i]) continue;
        Dart curr = A.edges[i];
        int cnt = 0;
        while (!visited[curr.id]) {
            visited[curr.id] = true;
            cnt++;
            curr = A.Phi(curr);
        }
        faces_sides.insert(cnt);
    }
    
    auto it = faces_sides.find(ssize(boundary));
    assert(it != faces_sides.end());
    faces_sides.erase(it);
    
    if (faces_sides.empty()) return true;
    if (*faces_sides.rbegin() > 7) return false;
    if (ssize(faces_sides) <= 2) return true;
    if (*prev(prev(faces_sides.end())) > 5) return false;
    
    return true;
}


vector<PlanarMap> Evolve(const PlanarMap& A) {
    vector<PlanarMap> res;
    
    // If A has a leaf, the next edge we add will be adjacent to a leaf.
    bool leaf_exists = false;
    for (int v = 0; v < A.V; v++) leaf_exists |= A.degree[v] == 1;
    
    // Edges belonging to the boundary of the unbounded face of A.
    vector<Dart> boundary = {A.boundary};
    while (A.Phi(boundary.back()) != A.boundary)
        boundary.push_back(A.Phi(boundary.back()));
    
    // Adding a new vertex "outside" of A.
    for (Dart d: boundary) {
        d = A.Rev(d);
        if (leaf_exists and A.degree[d.from] != 1) continue;
        
        // This check can be removed as it is enforced by the validator.
        //~ if (A.degree[d.from] >= 3) continue;
        
        PlanarMap A2 = A;
        A2.AddVertex(d);
        if (Validate(A2))
            res.push_back(A2);
    }
    
    
    for (Dart d1: boundary) for (Dart d2: boundary) {
        if (leaf_exists and A.degree[d1.to] != 1 and A.degree[d2.from] != 1)
            continue;
        
        // These checks can be removed as they are enforced by
        // the validator.
        //~ if (A.degree[d1.to] >= 3 or A.degree[d2.from] >= 3) continue;
        //~ if (d1.to == d2.from and A.degree[d1.to] >= 2) continue;
        
        PlanarMap A2 = A;
        pair<Dart, Dart> new_edge = A2.AddEdge(d1, d2);
        A2.boundary = new_edge.first;
        if (Validate(A2))
            res.push_back(A2);
        if (A.Phi(d1) == d2) {
            A2.boundary = new_edge.second;
            if (Validate(A2))
                res.push_back(A2);
        }
    }
    
    return res;
}


// This function supports both CombinatorialMap and PlanarMap.
template <typename T>
void RemoveDuplicates(vector<T>& V) {
    map<ULL, vector<T>> signature2list;
    for (const T& A: V) {
        signature2list[IsomorphismSignature(A)].push_back(A);
    }
    V.clear();
    
    for (const auto& pp: signature2list) {
        const auto& vec = pp.second;
        vector<T> unique_vec;
        for (const T& A: vec) {
            bool to_add = true;
            for (const T& B: unique_vec) {
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

vector<PlanarMap> Evolve(const vector<PlanarMap>& V) {
    vector<PlanarMap> new_V;
    for (const PlanarMap& A: V) {
        const vector<PlanarMap>& vec = Evolve(A);
        for (const PlanarMap& B: vec) new_V.push_back(B);
    }
    RemoveDuplicates(new_V);
    return new_V;
}


int main() {
    vector<PlanarMap> vec(1);
    
    for (int i = 2; i <= 10; i++) {
        vec = Evolve(vec);
        cout << i << ": " << vec.size() << endl;
    }
}
