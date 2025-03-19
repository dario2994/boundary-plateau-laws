#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "derived_once_nets.cpp"


vector<CombinatorialMap> parseDerivedTwiceNets() {
    ifstream in("data/derived_twice_computer.txt");
    int cnt;
    in >> cnt;
    vector<CombinatorialMap> res;
    for (int i = 0; i < cnt; i++) {
        CombinatorialMap A = parseCombinatorialMap(in);
        A.checkConsistency();
        res.push_back(A);
    }
    return res;
}

vector<DerivedOnceNet> addPin(const DerivedOnceNet& A) {
    vector<DerivedOnceNet> res;
    
    for (const vector<Dart>& face: A.faces()) {
        // We add a pin inside face.
        // A pin may generate one or two concave angles.
        for (const Dart& d: face) {
            DerivedOnceNet A2 = A;
            pair<Dart, Dart> pp = A2.edgeSubdivision(d);
            res.push_back(A2); // One edge coming out of the pin.
            for (const Dart& e: A2.face(pp.first)) {
                DerivedOnceNet A3 = A2;
                A3.edgeSubdivision(e);
                res.push_back(A3); // Two edges coming out of the pin.
            }
        }
    }
    
    return res;
}

bool validateDerivedOnceNet(const DerivedOnceNet& A) {
    // Check that convex faces have distinct edges and the number
    // of sides belongs to [2, 5].
    for (const vector<Dart>& face: A.faces()) {
        if (!A.isConvexFace(face[0]))
            continue;
        if (ssize(face) < 2 or ssize(face) > 5)
            return false;
        set<int> rev_sides;
        for (const Dart& d: face) {
            if (rev_sides.count(d.id))
                return false;
            rev_sides.insert(A.rev(d).id);
        }
    }
    
    // Check that any face has at most 2 obtuse angles.
    int total_area = 0;
    for (const vector<Dart>& face: A.faces()) {
        int cnt = 0;
        for (const Dart& d: face) {
            if (A.isObtuse(d)) 
                cnt++;
        }
        if (cnt > 2)
            return false;
            
        int area = 6 - ssize(face) + 2*cnt;
        if (area <= 0 or area > 12)
            return false;
        
        total_area += area;
    }
    assert(total_area == 12);
    
    // Two adjacent convex faces have constraints on their sizes.
    for (const Dart& d: A.darts()) {
        const Dart& e = A.rev(d);
        if (!A.isConvexFace(d) or !A.isConvexFace(e))
            continue;
        int d_sides = A.face(d).size();
        int e_sides = A.face(e).size();
        if (d_sides > e_sides) 
            swap(d_sides, e_sides);
        if (d_sides == 2 and e_sides != 2)
            return false;
        if (d_sides == 3 and e_sides == 5)
            return false;
    }
    
    // Two obtuse angles cannot be adjacent, unless they belong
    // to a triangle.
    for (const Dart& d: A.darts()) {
        if (A.isObtuse(d) and A.isObtuse(A.phi(d))
            and A.face(d).size() != 3)
            return false;
    }
    
    return true;
}

// Returns a copy of the original net, with some fake vertices added
// so that there are no more repeated edges.
// The list of integers returned represents the list of vertices that
// were added.
pair<DerivedOnceNet, vector<int>> removeRepeatedEdges(DerivedOnceNet A) {
    vector<int> new_vertices;
    while (1) {
        bool changed = false;
        set<pair<int,int>> edges;
        for (const Dart& d: A.darts()) {
            if (edges.count({d.from, d.to})) {
                A.edgeSubdivision(d);
                changed = true;
                new_vertices.push_back(A.V()-1);
                break;
            }
            edges.insert({d.from, d.to});
        }
        if (!changed) break;
    }
    return {A, new_vertices};
}

int main() {
    vector<CombinatorialMap> derived_twice = 
        parseDerivedTwiceNets();
    cout << "There are " << derived_twice.size()
         << " derived-twice net candidates." << endl << endl;
    
    CollectionUpToIsomorphism<DerivedOnceNet> derived_once;
    
    auto validate_and_add = [&](const DerivedOnceNet& A) {
        if (validateDerivedOnceNet(A))
            derived_once.add(A);
    };
    
    int it = 1;
    for (const CombinatorialMap& A: derived_twice) {
        vector<DerivedOnceNet> pin1 = addPin(A);
        for (const DerivedOnceNet& A2: pin1) {
            vector<DerivedOnceNet> pin2 = addPin(A2);
            for (const DerivedOnceNet& A3: pin2)
                validate_and_add(A3);
        }
        if (it % 10 == 0) {
            cout << "Processed successfully " << it 
                 << " derived-twice net candidates, which generated "
                 << derived_once.getList().size() 
                 << " derived-once net candidates." << endl;
        }
        it++;
    }
    cout << endl;
    
    const vector<DerivedOnceNet> derived_once_vec = derived_once.getList();
    
    cout << "We found " << derived_once_vec.size() 
         << " derived-once nets candidates." << endl;
         
    ofstream computer("data/derived_once_computer.txt");
    ofstream human("data/derived_once_human.txt");
    ofstream for_visualization("data/derived_once_for_visualization.txt");
    computer << derived_once_vec.size() << endl;
    for_visualization << derived_once_vec.size() << endl;
    for (int _ = 0; _ < ssize(derived_once_vec); _++) {
        const DerivedOnceNet& A = derived_once_vec[_];
        
        human << "Candidate " << _ << endl;
        human << A << endl;
        
        A.serialize(computer);
        
        // For visualization, we use the python library matplotlib.
        // This library does not support combinatorial maps with
        // repeated edges, so we create a version of the derived-once
        // net that does not contain any repeat edge.
        auto pp = removeRepeatedEdges(A);
        pp.first.serialize(for_visualization);
        for_visualization << pp.second.size() << "\n";
        for (int v: pp.second)
            for_visualization << v << " ";
        for_visualization << "\n";
    }
}
