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
    int _V, _E, _F; // Vertices, Edges, Faces
    vector<int> _degrees; // vertex degree
    // Two "copies" of each edge, one for each direction. 
    // Even loops are present twice.
    vector<Dart> _darts;
    
    int V() const { return _V; }
    int E() const { return _E; }
    int F() const { return _F; }
    int degree(int i) const { return _degrees[i]; }
    const Dart& dart(int i) const { return _darts[i]; }
    const vector<Dart>& darts() const { return _darts; }
    
    vector<int> _rev;
    vector<int> _next, _prev;
    
    
    // Adds an edge between u and v without properly taking care of
    // _next and _prev.
    // This takes only of the graph-theoretic component of a CombinatorialMap.
    // The user is responsible for properly updating _next and _prev
    // via _AssignNext.
    // Returns the new dart u->v and its reverse v->u.
    pair<Dart,Dart> _AddEdge(int u, int v) {
        assert(0 <= u and u < _V and 0 <= v and v < _V);
        _degrees[u]++, _degrees[v]++;
        _E++;
    
        Dart d1;
        d1.from = u, d1.to = v, d1.id = _darts.size();
        _darts.push_back(d1);
        Dart d2;
        d2.from = v, d2.to = u, d2.id = _darts.size();
        _darts.push_back(d2);
        
        _rev.push_back(d2.id);
        _rev.push_back(d1.id);
        
        _next.push_back(-1), _next.push_back(-1);
        _prev.push_back(-1), _prev.push_back(-1);
        
        return {d1, d2};
    }
    
    // Assigns Next(a) = b. (updates also Prev appropriately).
    void _AssignNext(const Dart& a, const Dart& b) {
        assert(a.from == b.from);
        _next[a.id] = b.id;
        _prev[b.id] = a.id;
    }
    
    const Dart& Rev(const Dart& e) const { return _darts[_rev[e.id]]; }
    // Next is counterclockwise
    const Dart& Next(const Dart& e) const { return _darts[_next[e.id]]; }
    const Dart& Prev(const Dart& e) const { return _darts[_prev[e.id]]; }
    // Phi traverses the boundary of a (bounded) face clockwise.
    int _phi(int id) const { return _next[_rev[id]]; }
    const Dart& Phi(const Dart& e) const { return _darts[_phi(e.id)]; } 
    
    
    // Create a new vertex and link it (with a new edge) to an existing one.
    // 
    // The new edge is connected so that Phi(bef) = new edge.
    pair<Dart,Dart> AddVertex(Dart bef) {
        _V++;
        _degrees.push_back(0);
        
        Dart bef_rev = Rev(bef);
        Dart bef_phi = Phi(bef);
        
        //~ cout << bef_rev << endl;
        
        pair<Dart,Dart> new_edge = _AddEdge(bef.to, V()-1);
        //~ cout << bef_rev << endl;
        
        //~ cout << new_edge.first << " " << new_edge.second << endl;
        
        _AssignNext(new_edge.second, new_edge.second);
        //~ cout << bef_rev << " " << new_edge.first << endl;
        _AssignNext(bef_rev, new_edge.first);
        _AssignNext(new_edge.first, bef_phi);
        
        return new_edge;
    }
    
    // Link with an edge two existing vertices.
    // Phi(bef) = new_edge 
    // Phi(new_edge) = aft
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
        _F++;
        
        return new_edge;
    }
    
    vector<vector<Dart>> GetFaces() const {
        vector<vector<Dart>> faces;
        vector<bool> visited(2*E());
        for (int i = 0; i < 2*E(); i++) {
            if (visited[i]) continue;
            Dart curr = dart(i);
            vector<Dart> face;
            while (!visited[curr.id]) {
                visited[curr.id] = true;
                face.push_back(curr);
                curr = Phi(curr);
            }
            faces.push_back(face);
        }
        return faces;
    }
    
    // Runs a dfs on the Darts Graph.
    void _dfs(int d, vector<bool>& visited, vector<int>& ord) const {
        if (visited[d]) return;
        visited[d] = true;
        ord.push_back(d);
        _dfs(_rev[d], visited, ord);
        _dfs(_next[d], visited, ord);
    }


    // Returns the list of the darts encountered when visiting the darts
    // graph, starting from d. 
    // The visit is a depth first search on the graph with vertices
    // given by the darts and edges:
    //  d -> rev[d]
    //  d -> next[d]
    // First, the visit traverses rev and then next.
    //
    // The ids of the darts are used for a concise representation.
    vector<int> OrderedVisit(const Dart& d) const {
        vector<bool> visited(2*E(), false);
        vector<int> ord;
        _dfs(d.id, visited, ord);
        return ord;
    }
};

// Starting from either an edge or a loop, one can generate any
// valid combinatorial map using AddVertex or AddEdge.

// Two vertices connected by an edge.
CombinatorialMap OneEdgeConstructor() {
    CombinatorialMap A;
    A._V = 2;
    A._degrees.resize(A.V(), 0);
    A._E = 0;
    A._F = 1;
    auto edge = A._AddEdge(0, 1);
    A._AssignNext(edge.first, edge.first);
    A._AssignNext(edge.second, edge.second);
    return A;
}

// One vertex with a loop.
CombinatorialMap OneLoopConstructor() {
    CombinatorialMap A;
    A._V = 1;
    A._degrees.resize(A.V(), 0);
    A._E = 0;
    A._F = 2;
    auto edge = A._AddEdge(0, 0);
    A._AssignNext(edge.first, edge.second);
    A._AssignNext(edge.second, edge.first);
    return A;
}

struct PlanarMap:CombinatorialMap {
    // An arbitrary dart on the boundary of the unbounded face.
    Dart boundary;
    
    PlanarMap(const CombinatorialMap& A): CombinatorialMap(A) {
        boundary = dart(0);
    }
    
    // Returns a list of all darts on the boundary of the unbounded face.
    vector<Dart> BoundaryEdges() const {
        vector<Dart> vec = {boundary};
        while (Phi(vec.back()) != boundary)
            vec.push_back(Phi(vec.back()));
        
        return vec;
    }
};

ostream& operator<<(ostream& out, const CombinatorialMap& A) {
    out << "{V = " << A.V();
    out << ", E = " << A.E();
    out << ", F = " << A.F();
    for (Dart d: A.darts()) 
        out << ", " << d;
    out << ", next = [";
    for (Dart d: A.darts())
        out << A.Next(d).id << ", ";
    out << "]}";
    return out;
}

ostream& operator<<(ostream& out, const PlanarMap& A) {
    out << "(" << static_cast<const CombinatorialMap&>(A) << ", boundary = [";
    for (Dart d: A.BoundaryEdges())
        out << d << ", ";
    out << "])";
    return out;
}

// Check whether the map iso between the darts of A and B induces
// an isomorphism of the combinatorial map.
// It is sufficient to check that it commutes with Rev and Next.
bool CheckIsomorphism(const CombinatorialMap& A,
                      const CombinatorialMap& B,
                      const vector<int>& iso) {
    int E = A.E();
    assert(A.E() == E and B.E() == E and ssize(iso) == 2*E);
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
    int E = A.E();
    
    vector<int> ordA = A.OrderedVisit(dart_A);
    vector<int> ordB = B.OrderedVisit(dart_B);
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
    if (A.V() != B.V() or A.E() != B.E() or A.F() != B.F())
        return {};

    for (Dart d: B.darts()) {
        const vector<int>& iso = FindIsomorphism(A, B, A.dart(0), d);
        if (!iso.empty())
            return iso;
    }
    return {};
}

vector<int> FindIsomorphism(const PlanarMap& A,
                            const PlanarMap& B) {
    if (A.V() != B.V() or A.E() != B.E() or A.F() != B.F())
        return {};
    vector<Dart> boundary = B.BoundaryEdges();
    for (Dart d: boundary) {
        const vector<int>& iso = FindIsomorphism(A, B, A.boundary, d);
        if (!iso.empty())
            return iso;
    }
    return {};
}

template<typename T>
bool AreIsomorphic(const T& A, const T& B) {
    return !FindIsomorphism(A, B).empty();
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
vector<vector<ULL>> MatrixSignature(const CombinatorialMap& A) {
    static const ULL p_stay = 424'142'591;
    static const ULL p_rev = 984'689'391;
    static const ULL p_next = 123'456'789;
    
    vector<vector<ULL>> mat(2*A.E(), vector<ULL>(2*A.E(), 0));
    for (Dart d: A.darts()) {
        int i = d.id;
        mat[i][i] += p_stay;
        mat[i][A.Rev(d).id] += p_rev;
        mat[i][A.Next(d).id] += p_next;
    }
    for (int _ = 0; _ < 8; _++) mat = MatrixMultiplication(mat, mat);
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
// on the "outside") into a valid geodesic net.
bool ValidateSubgraph(const PlanarMap& A) {
    if (A.E() > 30 or A.F() > 12)
        return false;
    
    int leaves_num = 0;
    for (int v = 0; v < A.V(); v++) if (A.degree(v) == 1) leaves_num++;
    assert(leaves_num <= 1);
    
    set<int> boundary_vertices;
    const vector<Dart>& boundary = A.BoundaryEdges();
    for (const Dart& d: boundary)
        boundary_vertices.insert(d.from);
    
    for (int v = 0; v < A.V(); v++) {
        if (!boundary_vertices.count(v) and A.degree(v) != 3) 
            return false;

        if (A.degree(v) > 3)
            return false;
    }
    
    // Computing the number of sides of all internal faces.
    const auto& faces = A.GetFaces();
    multiset<int> faces_sides;
    for (const auto& face: faces) faces_sides.insert(face.size());

    auto boundary_it = faces_sides.find(ssize(boundary));
    assert(boundary_it != faces_sides.end());
    faces_sides.erase(boundary_it); // Discarding the external face.
    
    // We check that: 
    // - there are at most two faces with strictly more than 5 sides.
    // - all faces have at most 7 sides.
    auto big = faces_sides.end();
    if (faces_sides.empty()) return true;
    if (*prev(big) > 7) return false;
    if (ssize(faces_sides) <= 2) return true;
    if (*prev(prev(prev(big))) > 5) return false;
    
    // It shall hold \sum (6-sides) = 12 over all faces.
    // We compute the internal contribution, which we already know,
    // as well as a lower bound on the external contribution.
    // Then, we check that the sum of the two is <= 12.
    int internal_contribution = 0;
    for (int s: faces_sides) internal_contribution += 6-s;
    
    // First, we determine a lower bound on the number of darts that
    // need to be added to complete the net.
    int missing_degrees = 0;
    for (int v: boundary_vertices) missing_degrees += 3-A.degree(v);
    
    // Then, we determine a lower bound on the number of darts which, 
    // in the final configuration, are not the boundary of a, currently,
    // internal face.
    int external_sides = ssize(boundary) + missing_degrees;
    int external_contribution = -external_sides;
    
    // Given the lower bound external_sides, we determine a lower bound
    // on the number of faces that will be added to reach the final
    // configuration.
    int external_faces = 0;
    int big_faces = (*prev(big) > 5)
                   +(*prev(prev(big)) > 5);
    while (external_sides > 0) {
        external_faces++;
        if (big_faces < 2) {
            big_faces++;
            external_sides -= 7;
        } else external_sides -= 5;
    }
    external_contribution += 6 * external_faces;
    
    // If it would be possible to have one more large face we must take
    // that into account to have a correct lower bound.
    if (big_faces < 2) external_contribution -= 1;
    
    if (internal_contribution + external_contribution > 12) 
        return false;
    
    return true;
}

// Checks whether A is a valid geodesic net.
bool ValidateNet(const CombinatorialMap& A) {    
    if (6*A.F() - 2*A.E() != 12) return false;

    for (int v = 0; v < A.V(); v++) {
        if (A.degree(v) != 3) 
            return false;
    }
    
    // TODO: Make a member in CombinatorialMap that returns
    //       the list of faces.
    const auto& faces = A.GetFaces();
    multiset<int> faces_sides;
    for (const auto& face: faces) faces_sides.insert(face.size());
        
    // We check that: 
    // - there are at most two faces with strictly more than 5 sides.
    // - all faces have at most 7 sides.
    auto big = faces_sides.end();
    assert(!faces_sides.empty());
    if (*prev(big) > 7) return false;
    if (ssize(faces_sides) > 2 and *prev(prev(prev(big))) > 5) return false;
        
    return true;
}


vector<PlanarMap> Evolve(const PlanarMap& A) {
    const vector<Dart>& boundary = A.BoundaryEdges();
    vector<PlanarMap> res;
    
    int leaf = -1;
    for (int v = 0; v < A.V(); v++) if (A.degree(v) == 1) leaf = v;
    
    auto add_to_res = [&](const PlanarMap& B) {
        if (leaf == -1 or B.degree(leaf) > 1)
            res.push_back(B);
    };
    
    for (Dart bef: boundary) {
        PlanarMap A2 = A;
        A2.AddVertex(bef);
        add_to_res(A2);
    }
        
    for (Dart bef: boundary) for (Dart aft: boundary) {
        PlanarMap A2 = A;
        pair<Dart, Dart> new_edge = A2.AddEdge(bef, aft);
        A2.boundary = new_edge.first;
        add_to_res(A2);
        if (bef.to == aft.from) {
            A2.boundary = new_edge.second;
            add_to_res(A2);
        }
    }
    
    return res;
}

template <typename T>
struct MaintainUnique {
    map<ULL, vector<T>> signature2list;
    
    void Add(const T& x) {
        bool is_new = true;
        vector<T>& vec = signature2list[IsomorphismSignature(x)];
        for (const T& y: vec) {
            if (AreIsomorphic(x, y)) {
                is_new = false;
                break;
            }
        }
        if (is_new)
            vec.push_back(x);
        assert(ssize(vec) < 30);
        //~ if (ssize(vec) > 1) {
            //~ auto A = vec[0];
            //~ auto B = vec[1];
            //~ int E = A.E;
            
            //~ cout << A << endl;
            //~ cout << B << endl;
            
            //~ vector<int> ordA = OrderedVisit(A, 8);
            //~ vector<int> ordB = OrderedVisit(B, 10);
            
            //~ cout << "ordA = ";
            //~ for (int i = 0; i < 2*E; i++)
                //~ cout << ordA[i] << ", ";
            //~ cout << endl;
            //~ cout << "ordB = ";
            //~ for (int i = 0; i < 2*E; i++)
                //~ cout << ordB[i] << ", ";
            //~ cout << endl;
            
            //~ vector<int> iso(2*E);
            //~ for (int i = 0; i < 2*E; i++) iso[ordA[i]] = ordB[i];
    
            //~ cout << "Is iso? " << CheckIsomorphism(A, B, iso) << endl;
            
            //~ exit(0);
        //~ }
    }
    
    vector<T> GetUnique() {
        vector<T> res;
        for (const auto& pp: signature2list) for (const T& x: pp.second) {
            res.push_back(x);
        }
        return res;
    }
};

vector<PlanarMap> Evolve(const vector<PlanarMap>& V) {
    MaintainUnique<PlanarMap> mu;
    for (const PlanarMap& A: V) {
        const vector<PlanarMap>& vec = Evolve(A);
        for (const PlanarMap& B: vec) {
            if (ValidateSubgraph(B))
                mu.Add(B);
        }
    }
    return mu.GetUnique();
}


int main() {
    vector<PlanarMap> vec = {OneLoopConstructor()};
    
    MaintainUnique<CombinatorialMap> nets;
    
    for (int E = 2; E <= 30; E++) {
        vec = Evolve(vec);
        
        // Adding a cycle of length E.
        PlanarMap cycle = OneEdgeConstructor();
        for (int k = 2; k < E; k++) {
            cycle.AddVertex(cycle.dart(2*(k-2)));
        }
        cycle.AddEdge(cycle.dart(2*(E-2)), cycle.dart(0));
        if (ValidateSubgraph(cycle)) 
            vec.push_back(cycle);
        
        for (const auto& A: vec) if (ValidateNet(A)) nets.Add(A);
        
        cout << E << ": " << vec.size() << " " << nets.GetUnique().size() << endl;
        
        if (vec.size() > 1e6) {
            break;
        }
    }
}


// 16: 39434 92

// 11:42, 21 282132 216
