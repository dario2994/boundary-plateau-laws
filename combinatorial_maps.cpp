#ifndef COMBINATORIAL_MAPS
#define COMBINATORIAL_MAPS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "utils.cpp"

struct Dart {
    int from, to;
    int id;
    
    bool operator ==(const Dart& other) const {
        return id == other.id;
    }
    
    bool operator !=(const Dart& other) const {
        return id != other.id;
    }
};

ostream& operator<<(ostream& out, const Dart& d) {
    out << "(" << d.id << ": " << d.from << "->" << d.to << ")";
    return out;
}


// Class representing a connected CombinatorialMap with at least one
// edge.
// Both loops and repeated edges are supported.
struct CombinatorialMap {
  private:
    // The default constructor is private as any user should construct
    // a CombinatorialMap via constructOneEdge() or constructOneLoop()
    CombinatorialMap() {}
  public:
    friend CombinatorialMap constructOneEdge();
    friend CombinatorialMap constructOneLoop();
    friend CombinatorialMap parseCombinatorialMap(istream& in);
    
  private:
    int _V, _E, _F; // Vertices, Edges, Faces
    vector<int> _degrees; // vertex degree
    // Two "copies" of each edge, one for each direction. 
    // Even loops are present twice.
    vector<Dart> _darts;
    
    vector<int> _rev;
    vector<int> _next, _prev;
    int _phi(int id) const { return _next[_rev[id]]; }
  
  public:  
    int V() const { return _V; }
    int E() const { return _E; }
    int F() const { return _F; }
    int degree(int i) const { return _degrees[i]; }
    const Dart& dart(int i) const { return _darts[i]; }
    const vector<Dart>& darts() const { return _darts; }
    
    const Dart& rev(const Dart& e) const { return _darts[_rev[e.id]]; }
    // next is counterclockwise
    const Dart& next(const Dart& e) const { return _darts[_next[e.id]]; }
    const Dart& prev(const Dart& e) const { return _darts[_prev[e.id]]; }
    // phi traverses the boundary of a (bounded) face clockwise.
    const Dart& phi(const Dart& e) const { return _darts[_phi(e.id)]; } 
    const Dart& phi(const Dart& e, int k) const { 
        if (k == 1)
            return phi(e);
        return phi(phi(e), k-1);
    }
    
  private:
    void checkDart(const Dart& d) const {
        assert(0 <= d.id and d.id < 2*_E);
        assert(_darts[d.id].from == d.from);
        assert(_darts[d.id].to == d.to);
    }
    
  
    // Adds an edge between u and v without properly taking care of
    // _next and _prev.
    // This changes only the graph-theoretic component of a CombinatorialMap.
    // The user is responsible for properly updating _next and _prev
    // via _assignNext.
    // Returns the new dart u->v and its reverse v->u.
    pair<Dart,Dart> createGraphEdge(int u, int v) {
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
    
    // Assigns next(a) = b. (updates also prev appropriately).
    void assignNext(const Dart& a, const Dart& b) {
        assert(a.from == b.from);
        _next[a.id] = b.id;
        _prev[b.id] = a.id;
    }
    
  public:
    // Create a new vertex and link it (with a new edge) to an existing one.
    // 
    // The new edge is connected so that phi(bef) = new edge.
    pair<Dart,Dart> addVertex(Dart bef) {
        checkDart(bef);
        
        _V++;
        _degrees.push_back(0);
        
        Dart bef_rev = rev(bef);
        Dart bef_phi = phi(bef);
        
        pair<Dart,Dart> new_edge = createGraphEdge(bef.to, V()-1);
        assignNext(new_edge.second, new_edge.second);
        assignNext(bef_rev, new_edge.first);
        assignNext(new_edge.first, bef_phi);
        
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
    pair<Dart,Dart> addEdge(Dart bef, Dart aft) {
        checkDart(bef);
        checkDart(aft);
        
        // We check that aft = phi^k(bef) for some k >= 0, which is
        // equivalent to the fact that bef and aft belong to the
        // boundary of the same face.
        auto it = bef;
        while (it != aft) {
            it = phi(it);
            assert(it != bef);
        }
        
        // We Create the new edge
        pair<Dart, Dart> new_edge = createGraphEdge(bef.to, aft.from);
        
        // We update the next function
        const Dart& bef_next = next(rev(bef));
        const Dart& aft_prev = prev(aft);
        
        assignNext(rev(bef), new_edge.first);
        if (bef.to != aft.from) {
            assignNext(new_edge.first, bef_next);
            assignNext(aft_prev, new_edge.second);
        } else {
            assignNext(new_edge.first, new_edge.second);
        }
        assignNext(new_edge.second, aft);
        
        // We have added a new face in the process.
        _F++;
        
        return new_edge;
    }
    
    // Replaces the edge associated to d with two edges linked by
    // a vertex. It does not change the structure of the combinatorial
    // map.
    // This function is very involved for two reasons:
    // 1. There are a few pesky cases (for example when one of the two
    //    endpoints of d is a leaf, or if d is a loop).
    // 2. The required operation invalidates references to darts.
    pair<Dart, Dart> edgeSubdivision(Dart d) {
        checkDart(d);
        // We use id to refer to darts as references are going to be 
        // invalidated once we call createGraphEdge
        int d_id = d.id;
        int e_id = rev(d).id;
        
        Dart prev_e = prev(dart(e_id));
        Dart next_e = next(dart(e_id));
        
        // Creating a new vertex.
        _V++;
        _degrees.push_back(0);
        
        // x is the central vertex
        int x = _V-1;
        int y = dart(d_id).to;
        
        // We update e and d to be adjacent to x.
        _darts[d_id].to = x;
        _darts[e_id].from = x;
        _degrees[x] = 1;
        _degrees[y]--;
        
        // We create an edge from x to y.
        pair<Dart, Dart> pp = createGraphEdge(x, y);
        
        // We fix prev/next in x and y.
        assignNext(dart(e_id), pp.first);
        assignNext(pp.first, dart(e_id));
        
        if (_degrees[y] > 1) {
            assignNext(prev_e, pp.second);
            assignNext(pp.second, next_e);
        } else {
            assignNext(pp.second, pp.second);
        }
        
        
        return {dart(d_id), pp.first};
    }
    
    vector<Dart> face(Dart d) const {
        checkDart(d);
        
        Dart curr = d;
        vector<Dart> face;
        do {
            face.push_back(curr);
            curr = phi(curr);
        } while (curr != d);
        return face;
    }
    
    vector<vector<Dart>> faces() const {
        vector<vector<Dart>> faces;
        vector<bool> visited(2*E());
        for (int i = 0; i < 2*E(); i++) {
            if (visited[i]) continue;
            Dart curr = dart(i);
            vector<Dart> face;
            while (!visited[curr.id]) {
                visited[curr.id] = true;
                face.push_back(curr);
                curr = phi(curr);
            }
            faces.push_back(face);
        }
        return faces;
    }
    
    // A function that asserts if the class does not satisfy one of
    // the assumptions that should be satisfied by a well-structured
    // combinatorial map.
    // The functions available to the user to create and modify a
    // combinatorial map (constructOneEdge, constructOneLoop, 
    // AddVertex, addEdge, edgeSubdivision) maintain the validity
    // of all the assumptions.
    void checkConsistency() const {
        // First, we check the consistency of the graph structure.
        assert(_V > 0);
        assert(_E > 0);
        assert(ssize(_darts) == 2*_E);
        assert(ssize(_degrees) == _V);
        
        // Checking that degrees are correct and edges are reasonable.
        vector<int> new_degrees(_V, 0);
        for (const Dart& d: _darts) {
            assert(0 <= d.from and d.from < _V);
            assert(0 <= d.to and d.to < _V);
            new_degrees[d.from]++;
        }
        assert(_degrees == new_degrees);
        
        // Checking connectedness. (we use a DSU)
        vector<int> father(_V);
        for (int v = 0; v < _V; v++) father[v] = v;
        function<int(int)> ancestor = [&](int v) {
            if (v != father[v]) {
                father[v] = ancestor(father[v]);
            }
            return father[v];
        };
        for (const Dart& d: _darts) {
            int x = d.from, y = d.to;
            x = ancestor(x), y = ancestor(y);
            father[x] = y;
        }
        int root = ancestor(0);
        for (int v = 0; v < _V; v++) 
            assert(ancestor(v) == root);
        
        // We move to checking the consistency of the combinatorial
        // map structure.
        assert(ssize(_rev) == 2*_E);
        assert(ssize(_next) == 2*_E);
        assert(ssize(_prev) == 2*_E);
        
        // Checking that _rev, _next, _prev belong to [0, 2*_E).
        for (int i = 0; i < 2*_E; i++) {
            assert(0 <= _rev[i] and _rev[i] < 2*_E);
            assert(0 <= _next[i] and _next[i] < 2*_E);
            assert(0 <= _prev[i] and _prev[i] < 2*_E);
        }
        
        // Checking that rev is consistent.
        for (const Dart& d: _darts) {
            const Dart& e = rev(d);
            assert(d.from == e.to);
            assert(d.to == e.from);
            assert(rev(e) == d);
        }
        
        // Checking that prev(next) = id.
        for (const Dart& d: _darts) {
            assert(prev(next(d)) == d);
        }
        
        // Checking that, for every vertex, next is a cycle over
        // all darts stemming from that vertex.
        for (const Dart& d: _darts) {
            int q = _degrees[d.from];
            Dart curr = next(d);
            for (int i = 0; i < q-1; i++) {
                assert(curr != d);
                curr = next(curr);
            }
            assert(curr == d);
        }
        
        // Checking that the number of faces is correct.
        assert(ssize(faces()) == _F);
    }
    
  private:
    // Runs a dfs on the Darts Graph.
    void _dfs(int d, vector<bool>& visited, vector<int>& ord) const {
        if (visited[d]) return;
        visited[d] = true;
        ord.push_back(d);
        _dfs(_rev[d], visited, ord);
        _dfs(_next[d], visited, ord);
    }

  public:
    // Returns the list of the darts encountered when visiting the darts
    // graph, starting from d. 
    // The visit is a depth first search on the graph with vertices
    // given by the darts and edges:
    //  d -> rev[d]
    //  d -> next[d]
    // First, the visit traverses rev and then next.
    //
    // The ids of the darts are used for a concise representation.
    vector<int> orderedVisit(const Dart& d) const {
        checkDart(d);
        
        vector<bool> visited(2*E(), false);
        vector<int> ord;
        _dfs(d.id, visited, ord);
        return ord;
    }
    
  protected:
    
    // Returns the adjacency matrix of the darts graph, giving different
    // weights to rev and next edges.
    // It is a 2E x 2E square matrix. 
    // Let M_A be the matrix associated to the combinatorial map A.
    // If f is a bijection from the darts of A to the darts of B that
    // is a isomorphism of combinatorial maps, then M_B is obtained
    // permuting the rows/columns of M_A according to f.
    vector<vector<ULL>> matrixSignature() const {
        static const ULL p_rev = 984'689'395;
        static const ULL p_next = 123'456'789;
        
        vector<vector<ULL>> mat(2*E(), vector<ULL>(2*E(), 0));
        for (Dart d: darts()) {
            int i = d.id;
            mat[i][rev(d).id] += p_rev;
            mat[i][next(d).id] += p_next;
        }
        return mat;
    }
    
  public:

    // IsomorphismsSignature computes a large power of MatrixSignature, 
    // and then applies a nonlinear hash function.
    //
    // One may think that even much simpler implementations of 
    // isomorphismSignature could have a similar ability to tell 
    // whether two combinatorial maps are isomorphic. 
    // For example, why not simply returning the trace?
    // Well, because it does not work.
    // The current implementations are the simplest we could find that
    // don't have any collisions during the execution of the algorithm.
    ULL isomorphismSignature() const {
        const auto& mat = largeMatrixPower(matrixSignature());
        
        ULL res = 0;
        for (const Dart& d: darts()) {
            res += mat[d.id][rev(d).id] ^ mat[d.id][next(d).id];
        }
        
        return res;
    }
    
    void serialize(ostream& out) const {
        out << V() << " " << E() << " " << F() << endl;
        for (int i = 0; i < 2*E(); i++) {
            Dart d = dart(i);
            out << d.from << " " << d.to << endl;
        }
        for (int i = 0; i < 2*E(); i++) {
            Dart d = dart(i);
            out << rev(d).id << " ";
        }
        out << endl;
        
        for (int i = 0; i < 2*E(); i++) {
            Dart d = dart(i);
            out << next(d).id << " ";
        }
        out << endl; 
    }
};

// Starting from either an edge or a loop, one can generate any
// valid combinatorial map using AddVertex or addEdge.

// Two vertices connected by an edge.
CombinatorialMap constructOneEdge() {
    CombinatorialMap A;
    A._V = 2;
    A._degrees.resize(A.V(), 0);
    A._E = 0;
    A._F = 1;
    auto edge = A.createGraphEdge(0, 1);
    A.assignNext(edge.first, edge.first);
    A.assignNext(edge.second, edge.second);
    return A;
}

// One vertex with a loop.
CombinatorialMap constructOneLoop() {
    CombinatorialMap A;
    A._V = 1;
    A._degrees.resize(A.V(), 0);
    A._E = 0;
    A._F = 2;
    auto edge = A.createGraphEdge(0, 0);
    A.assignNext(edge.first, edge.second);
    A.assignNext(edge.second, edge.first);
    return A;
}

CombinatorialMap parseCombinatorialMap(istream& in) {
    CombinatorialMap A;
    in >> A._V >> A._E >> A._F;
    int n = 2*A._E;
    
    A._degrees.resize(A._V, 0);
    
    A._darts.resize(n);
    for (int i = 0; i < n; i++) {
        Dart d;
        d.id = i;
        in >> d.from >> d.to;
        A._darts[i] = d;
        A._degrees[d.from]++;
    }
    A._rev.resize(n);
    A._next.resize(n);
    A._prev.resize(n);
    
    for (int i = 0; i < n; i++) in >> A._rev[i];
    for (int i = 0; i < n; i++) {
        in >> A._next[i];
        A._prev[A._next[i]] = i;
    }
    return A;
}


ostream& operator<<(ostream& out, const CombinatorialMap& A) {
    out << "V = " << A.V();
    out << ", E = " << A.E();
    out << ", F = " << A.F();
    out << "\n";
    out << "Darts: ";
    for (Dart d: A.darts()) 
        out << ", " << d;
    out << "\n";
    out << "Next:\n";
    for (Dart d: A.darts())
        out << d << " ----> " << A.next(d) << "\n";
    return out;
}



// Functions to determine whether two combinatorial maps are isomorphic.


// Check whether the map iso between the darts of A and B induces
// an isomorphism of the combinatorial map.
// It is sufficient to check that it commutes with rev and next.
bool checkIsomorphism(const CombinatorialMap& A,
                      const CombinatorialMap& B,
                      const vector<int>& iso) {
    int E = A.E();
    assert(A.E() == E and B.E() == E and ssize(iso) == 2*E);
    for (int i = 0; i < 2*E; i++) {
        const Dart& a = A.dart(i);
        const Dart& b = B.dart(iso[i]);
        if (B.rev(b).id != iso[A.rev(a).id]) return false;
        if (B.next(b).id != iso[A.next(a).id]) return false;
    }
    return true;
}

// Returns, if it exists, an isomorphism between A and B that maps dart_A
// into dart_B.
// The ids of the darts are used for a concise representation.
vector<int> findIsomorphism(const CombinatorialMap& A,
                            const CombinatorialMap& B,
                            Dart dart_A,
                            Dart dart_B) {
    int E = A.E();
    
    vector<int> ordA = A.orderedVisit(dart_A);
    vector<int> ordB = B.orderedVisit(dart_B);
    vector<int> iso(2*E);
    for (int i = 0; i < 2*E; i++) iso[ordA[i]] = ordB[i];
    
    if (checkIsomorphism(A, B, iso))
        return iso;
    return {};
}


// Returns, if it exists, a bijection between the darts that is an
// isomorphism of the combinatorial maps.
// The ids of the darts are used for a concise representation.
vector<int> findIsomorphism(const CombinatorialMap& A,
                            const CombinatorialMap& B) {
    if (A.V() != B.V() or A.E() != B.E() or A.F() != B.F())
        return {};

    for (Dart d: B.darts()) {
        const vector<int>& iso = findIsomorphism(A, B, A.dart(0), d);
        if (!iso.empty())
            return iso;
    }
    return {};
}


#endif
