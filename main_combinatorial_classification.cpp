#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "utils.cpp"
#include "combinatorial_maps.cpp"
#include "planar_maps.cpp"


// Checks whether A is a valid geodesic net.
bool validateNet(const CombinatorialMap& A) {    
    // A must be embeddable on the sphere.
    if (A.V() - A.E() + A.F() != 2) return false;

    // A must be 3-regular.
    for (int v = 0; v < A.V(); v++) {
        if (A.degree(v) != 3) 
            return false;
    }
    
    const auto& faces = A.getFaces();
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


int computeExternalContribution(int big_faces, int darts) {
    int res = -darts;
    while (darts > 0) {
        res += 6;
        darts -= (big_faces > 0) ? 7 : 5;
        big_faces--;
    }
    return res;
}

int memo_external_contribution[3][80];
void precomputeExternalContribution() {
    for (int big_faces = 0; big_faces <= 2; big_faces++) {
        for (int darts = 70; darts >= 0; darts--) {
            int& z = memo_external_contribution[big_faces][darts];
            z = 1000;
            if (darts <= 62)
                z = memo_external_contribution[big_faces][darts+6];
            z = min(z, computeExternalContribution(big_faces, darts));
        }
    }
}

// Checks whether A can eventually be transformed (by attaching edges
// on the "outside") into a valid geodesic net.
// It (often) gives false positives.
bool validateSubgraph(const PlanarMap& A) {
    if (A.E() > 30 or A.F() > 12)
        return false;
    
    int leaves_num = 0;
    for (int v = 0; v < A.V(); v++) if (A.degree(v) == 1) leaves_num++;
    assert(leaves_num <= 1);
    
    set<int> boundary_vertices;
    const vector<Dart>& boundary = A.boundaryDarts();
    for (const Dart& d: boundary)
        boundary_vertices.insert(d.from);
    
    for (int v = 0; v < A.V(); v++) {
        if (!boundary_vertices.count(v) and A.degree(v) != 3) 
            return false;

        if (A.degree(v) > 3)
            return false;
    }
    
    // Computing the number of sides of all internal faces.
    const auto& faces = A.getFaces();
    multiset<int> faces_sides;
    for (const auto& face: faces) faces_sides.insert(face.size());

    auto boundary_it = faces_sides.find(ssize(boundary));
    assert(boundary_it != faces_sides.end());
    faces_sides.erase(boundary_it); // Discarding the external face.
    
    // We check that: 
    // - there are at most two faces with strictly more than 5 sides.
    // - all faces have at most 7 sides.
    auto big = faces_sides.end();
    assert(!faces_sides.empty());
    if (*prev(big) > 7) return false;
    if (ssize(faces_sides) >= 3 and *prev(prev(prev(big))) > 5) return false;
    
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
    if (missing_degrees % 2) missing_degrees += 3;
    if (A.E() + missing_degrees/2 > 30) return false;
    
    // Then, we determine a lower bound on the number of darts which, 
    // in the final configuration, are not (currently) the boundary of 
    // an internal face.
    int external_sides = ssize(boundary) + missing_degrees;
    
    // We compute the number of faces with at least 6 sides already
    // present.
    int big_faces = (*prev(big) > 5)
                   +(ssize(faces_sides) >= 2 and *prev(prev(big)) > 5);
    
    // Given the lower bound external_sides, we determine a lower bound
    // on \sum (6-sides) over the faces that will be added to reach 
    // the final configuration.
    int external_contribution = 
        memo_external_contribution[2-big_faces][external_sides];

    if (internal_contribution + external_contribution > 12) 
        return false;
    
    return true;
}


vector<PlanarMap> evolve(const PlanarMap& A) {
    const vector<Dart>& boundary = A.boundaryDarts();
    vector<PlanarMap> res;
    
    int leaf = -1;
    for (int v = 0; v < A.V(); v++) if (A.degree(v) == 1) leaf = v;
    
    auto add_to_res = [&](const PlanarMap& B, Dart new_edge) {
        if (B.degree(new_edge.from) > 3 or B.degree(new_edge.to) > 3)
            return;
        
        if (leaf != -1 and B.degree(leaf) == 1)
            return;
        
        res.push_back(B);
        return;
    };
    
    for (Dart bef: boundary) {
        PlanarMap A2 = A;
        const auto& new_edge = A2.addVertex(bef);
        add_to_res(A2, new_edge.first);
    }
        
    for (Dart bef: boundary) for (Dart aft: boundary) {
        PlanarMap A2 = A;
        pair<Dart, Dart> new_edge = A2.addEdge(bef, aft);
        A2.setBoundaryDart(new_edge.first);
        add_to_res(A2, new_edge.first);
        if (bef.to == aft.from) {
            A2.setBoundaryDart(new_edge.second);
            add_to_res(A2, new_edge.second);
        }
    }
    
    return res;
}

vector<PlanarMap> evolve(const vector<PlanarMap>& V) {
    CollectionUpToIsomorphism<PlanarMap> mu;
    for (const PlanarMap& A: V) {
        const vector<PlanarMap>& vec = evolve(A);
        for (const PlanarMap& B: vec) {
            if (validateSubgraph(B))
                mu.add(B);
        }
    }
    return mu.getList();
}


int main() {
    precomputeExternalContribution();
    
    vector<PlanarMap> vec = {constructOneLoop()};
    
    CollectionUpToIsomorphism<CombinatorialMap> nets;
    
    cout << "Edges = k\tPlanar maps with k edges\tTwice-derived nets with <= k edges" << endl;
    
    for (int E = 2; E <= 30; E++) {
        vec = evolve(vec);
        
        // Adding a cycle of length E.
        PlanarMap cycle = constructOneEdge();
        for (int k = 2; k < E; k++) {
            cycle.addVertex(cycle.dart(2*(k-2)));
        }
        cycle.addEdge(cycle.dart(2*(E-2)), cycle.dart(0));
        if (validateSubgraph(cycle)) 
            vec.push_back(cycle);
        
        for (const auto& A: vec) {
            if (validateNet(A)) {
                A.checkConsistency();
                nets.add(A);
            }
        }
        
        cout << E << "\t\t" << vec.size() << "\t\t\t\t" << nets.getList().size() << endl;
    }
    
    // We save the twice-derived nets into two files.
    // One of the files is optimized for being read by a human, the
    // other one is optimized for being parsed by a computer.
    vector<CombinatorialMap> res = nets.getList();
    ofstream human("data/derived_twice_human.txt");
    ofstream computer("data/derived_twice_computer.txt");
    
    computer << ssize(res) << endl;
    
    for (int _ = 0; _ < ssize(res); _++) {
        const CombinatorialMap& A = res[_];
        human << "Candidate " << _ << endl;
        human << A << endl << endl;
        
        A.serialize(computer);
    }
}

//~ Sul laptop di Camillo, una precedente versione (ma forse equivalente):
//~
//~ 2: 2 0
//~ 3: 6 2
//~ 4: 8 2
//~ 5: 16 2
//~ 6: 42 6
//~ 7: 89 6
//~ 8: 203 6
//~ 9: 500 18
//~ 10: 1182 18
//~ 11: 2665 18
//~ 12: 6184 45
//~ 13: 13157 45
//~ 14: 26871 45
//~ 15: 54319 92
//~ 16: 98182 92
//~ 17: 187961 92
//~ 18: 312123 156
//~ 19: 548530 156
//~ 20: 854319 156
//~ 21: 1353619 216
//~ 22: 1935872 216
//~ 23: 2593212 216
//~ 24: 2609982 276
//~ 25: 1392917 276
//~ 26: 460832 276
//~ 27: 107835 311
//~ 28: 18496 311
//~ 29: 2135 311
//~ 30: 170 333

//~ real    187m27.562s
//~ user    186m56.867s
//~ sys    0m3.095s
