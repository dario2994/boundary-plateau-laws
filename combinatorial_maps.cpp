#ifndef COMBINATORIAL_MAPS
#define COMBINATORIAL_MAPS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

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
    // a CombinatorialMap via OneEdgeConstructor() or OneLoopConstructor()
    CombinatorialMap() {}
  public:
    friend CombinatorialMap OneEdgeConstructor();
    friend CombinatorialMap OneLoopConstructor();
    friend CombinatorialMap ParserConstructor(istream& in);
    
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
    
    const Dart& Rev(const Dart& e) const { return _darts[_rev[e.id]]; }
    // Next is counterclockwise
    const Dart& Next(const Dart& e) const { return _darts[_next[e.id]]; }
    const Dart& Prev(const Dart& e) const { return _darts[_prev[e.id]]; }
    // Phi traverses the boundary of a (bounded) face clockwise.
    const Dart& Phi(const Dart& e) const { return _darts[_phi(e.id)]; } 
    
  private:
    void CheckDart(const Dart& d) const {
        assert(0 <= d.id and d.id < 2*_E);
        assert(_darts[d.id].from == d.from);
        assert(_darts[d.id].to == d.to);
    }
    
  
    // Adds an edge between u and v without properly taking care of
    // _next and _prev.
    // This changes only the graph-theoretic component of a CombinatorialMap.
    // The user is responsible for properly updating _next and _prev
    // via _AssignNext.
    // Returns the new dart u->v and its reverse v->u.
    pair<Dart,Dart> CreateGraphEdge(int u, int v) {
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
    void AssignNext(const Dart& a, const Dart& b) {
        assert(a.from == b.from);
        _next[a.id] = b.id;
        _prev[b.id] = a.id;
    }
    
  public:
    // Create a new vertex and link it (with a new edge) to an existing one.
    // 
    // The new edge is connected so that Phi(bef) = new edge.
    pair<Dart,Dart> AddVertex(Dart bef) {
        CheckDart(bef);
        
        _V++;
        _degrees.push_back(0);
        
        Dart bef_rev = Rev(bef);
        Dart bef_phi = Phi(bef);
        
        pair<Dart,Dart> new_edge = CreateGraphEdge(bef.to, V()-1);
        AssignNext(new_edge.second, new_edge.second);
        AssignNext(bef_rev, new_edge.first);
        AssignNext(new_edge.first, bef_phi);
        
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
        CheckDart(bef);
        CheckDart(aft);
        
        // We check that aft = phi^k(bef) for some k >= 0, which is
        // equivalent to the fact that bef and aft belong to the
        // boundary of the same face.
        auto it = bef;
        while (it != aft) {
            it = Phi(it);
            assert(it != bef);
        }
        
        // We Create the new edge
        pair<Dart, Dart> new_edge = CreateGraphEdge(bef.to, aft.from);
        
        // We update the next function
        const Dart& bef_next = Next(Rev(bef));
        const Dart& aft_prev = Prev(aft);
        
        AssignNext(Rev(bef), new_edge.first);
        if (bef.to != aft.from) {
            AssignNext(new_edge.first, bef_next);
            AssignNext(aft_prev, new_edge.second);
        } else {
            AssignNext(new_edge.first, new_edge.second);
        }
        AssignNext(new_edge.second, aft);
        
        // We have added a new face in the process.
        _F++;
        
        return new_edge;
    }
    
    // Replaces the edge associated to d with two edges linked by
    // a vertex. It does not change the structure of the combinatorial
    // map.
    // It handles also the case of d being a loop.
    pair<Dart, Dart> EdgeSubdivision(Dart _d) {
        CheckDart(_d);
        
        Dart& d = _darts[_d.id];
        Dart& e = _darts[_rev[d.id]];
        Dart prev_e = Prev(e);
        Dart next_e = Next(e);
        
        // Creating a new vertex.
        _V++;
        _degrees.push_back(0);
        
        // x is the central vertex
        int x = _V-1;
        int y = d.to;
        
        // We update e and d to be adjacent to x.
        d.to = x;
        e.from = x;
        _degrees[x] = 1;
        _degrees[y]--;
        
        // We create an edge from x to y.
        pair<Dart, Dart> pp = CreateGraphEdge(x, y);
        
        exit(0);
        
        // We fix prev/next in x and y.
        AssignNext(e, pp.first);
        AssignNext(pp.first, e);
        AssignNext(prev_e, pp.second);
        AssignNext(pp.second, next_e);
        
        
        return {d, pp.first};
    }
    
    vector<Dart> GetFace(Dart d) const {
        CheckDart(d);
        
        Dart curr = d;
        vector<Dart> face;
        do {
            face.push_back(curr);
            curr = Phi(curr);
        } while (curr != d);
        return face;
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
    
    // A function that asserts if the class does not satisfy one of
    // the assumptions that should be satisfied by a well-structured
    // combinatorial map.
    // The functions available to the user to create and modify a
    // combinatorial map (OneEdgeConstructor, OneLoopConstructor, 
    // AddVertex, AddEdge, EdgeSubdivision) maintain the validity
    // of all the assumptions.
    void CheckConsistency() const {
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
        
        // Checking that Rev is consistent.
        for (const Dart& d: _darts) {
            const Dart& e = Rev(d);
            assert(d.from == e.to);
            assert(d.to == e.from);
            assert(Rev(e) == d);
        }
        
        // Checking that prev(next) = id.
        for (const Dart& d: _darts) {
            assert(Prev(Next(d)) == d);
        }
        
        // Checking that, for every vertex, Next is a cycle over
        // all darts stemming from that vertex.
        for (const Dart& d: _darts) {
            int q = _degrees[d.from];
            Dart curr = Next(d);
            for (int i = 0; i < q-1; i++) {
                assert(curr != d);
                curr = Next(curr);
            }
            assert(curr == d);
        }
        
        // Checking that the number of faces is correct.
        assert(ssize(GetFaces()) == _F);
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
    vector<int> OrderedVisit(const Dart& d) const {
        CheckDart(d);
        
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
    auto edge = A.CreateGraphEdge(0, 1);
    A.AssignNext(edge.first, edge.first);
    A.AssignNext(edge.second, edge.second);
    return A;
}

// One vertex with a loop.
CombinatorialMap OneLoopConstructor() {
    CombinatorialMap A;
    A._V = 1;
    A._degrees.resize(A.V(), 0);
    A._E = 0;
    A._F = 2;
    auto edge = A.CreateGraphEdge(0, 0);
    A.AssignNext(edge.first, edge.second);
    A.AssignNext(edge.second, edge.first);
    return A;
}

CombinatorialMap ParserConstructor(istream& in) {
    CombinatorialMap A;
    in >> A._V >> A._E >> A._F;
    int n = 2*A._E;
    
    A._darts.resize(n);
    for (int i = 0; i < n; i++) {
        Dart d;
        d.id = i;
        in >> d.from >> d.to;
        A._darts.push_back(d);
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

class PlanarMap: public CombinatorialMap {
  private:
    // An arbitrary dart on the boundary of the unbounded face.
    Dart boundary_dart;
  public:
    PlanarMap(const CombinatorialMap& A): CombinatorialMap(A) {
        boundary_dart = dart(0);
    }
    
    // Returns a list of all darts on the boundary of the unbounded face.
    vector<Dart> BoundaryDarts() const {
        vector<Dart> vec = {boundary_dart};
        while (Phi(vec.back()) != boundary_dart)
            vec.push_back(Phi(vec.back()));
        
        return vec;
    }
    
    // We need to override the definition of EdgeSubdivision as it
    // may invalidate boundary_dart.
    pair<Dart, Dart> EdgeSubdivision(Dart d) {
        vector<Dart> boundary_darts = BoundaryDarts();
        bool on_boundary = false;
        bool on_boundary_rev = false;
        for (const Dart& e: boundary_darts) {
            on_boundary |= d == e;
            on_boundary_rev |= Rev(d) == e;
        }
        
        pair<Dart,Dart> pp = CombinatorialMap::EdgeSubdivision(d);
        if (on_boundary)
            boundary_dart = pp.first;
        if (on_boundary_rev)
            boundary_dart = pp.second;
        
        return pp;
    }
    
    void SetBoundaryDart(Dart d) {
        boundary_dart = d;
    }
    const Dart& GetBoundaryDart() const { return boundary_dart; }
};

ostream& operator<<(ostream& out, const PlanarMap& A) {
    out << "(" << static_cast<const CombinatorialMap&>(A) << ", boundary = [";
    for (Dart d: A.BoundaryDarts())
        out << d << ", ";
    out << "])";
    return out;
}


#endif
