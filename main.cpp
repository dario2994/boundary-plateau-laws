#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "isomorphisms.cpp"


// Checks whether A is a valid geodesic net.
bool ValidateNet(const CombinatorialMap& A) {    
    if (6*A.F() - 2*A.E() != 12) return false;

    for (int v = 0; v < A.V(); v++) {
        if (A.degree(v) != 3) 
            return false;
    }
    
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

// Checks whether A can eventually be transformed (by attaching edges
// on the "outside") into a valid geodesic net.
// It (often) gives false positives.
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
    
    vector<CombinatorialMap> res = nets.GetUnique();
    
    ofstream human("list_combinatorial_maps_human.txt");
    ofstream computer("list_combinatorial_maps_computer.txt");
    
    computer << ssize(res) << endl;
    
    for (int _ = 0; _ < ssize(res); _++) {
        const CombinatorialMap& A = res[_];
        human << "Candidate " << _+1 << endl;
        human << A << endl << endl;
        
        computer << A.V() << " " << A.E() << " " << A.F() << endl;
        for (int i = 0; i < 2*A.E(); i++) {
            Dart d = A.dart(i);
            computer << d.from << " " << d.to << endl;
        }
        for (int i = 0; i < 2*A.E(); i++) {
            Dart d = A.dart(i);
            computer << A.Rev(d).id << " ";
        }
        computer << endl;
        
        for (int i = 0; i < 2*A.E(); i++) {
            Dart d = A.dart(i);
            computer << A.Next(d).id << " ";
        }
        computer << endl;    
    }
}

// 13: 7517 45, 2.789s

// 16: 39434 92

// 20: 200793 156, 6m6.114s

