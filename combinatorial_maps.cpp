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
    // Adds an edge between u and v without properly taking care of
    // _next and _prev.
    // This changes only the graph-theoretic component of a CombinatorialMap.
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
    
  public:
    // Create a new vertex and link it (with a new edge) to an existing one.
    // 
    // The new edge is connected so that Phi(bef) = new edge.
    pair<Dart,Dart> AddVertex(Dart bef) {
        _V++;
        _degrees.push_back(0);
        
        Dart bef_rev = Rev(bef);
        Dart bef_phi = Phi(bef);
        
        pair<Dart,Dart> new_edge = _AddEdge(bef.to, V()-1);
        _AssignNext(new_edge.second, new_edge.second);
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
    
    vector<Dart> GetFace(Dart d) const {
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

ostream& operator<<(ostream& out, const PlanarMap& A) {
    out << "(" << static_cast<const CombinatorialMap&>(A) << ", boundary = [";
    for (Dart d: A.BoundaryEdges())
        out << d << ", ";
    out << "])";
    return out;
}


#endif
