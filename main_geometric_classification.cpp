#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "derived_once_nets.cpp"


vector<CombinatorialMap> ParseDerivedTwiceNetCandidates() {
    ifstream in("data/derived_twice_computer.txt");
    int cnt;
    in >> cnt;
    vector<CombinatorialMap> res;
    for (int i = 0; i < cnt; i++) {
        CombinatorialMap A = ParserConstructor(in);
        A.CheckConsistency();
        res.push_back(A);
    }
    return res;
}

vector<DerivedOnceNet> AddPin(const DerivedOnceNet& A) {
    vector<DerivedOnceNet> res;
    
    for (const vector<Dart>& face: A.GetFaces()) {
        // We add a pin inside face.
        // A pin may generate one or two concave angles.
        for (const Dart& d: face) {
            DerivedOnceNet A2 = A;
            pair<Dart, Dart> pp = A2.EdgeSubdivision(d);
            res.push_back(A2); // One edge coming out of the pin.
            for (const Dart& e: A2.GetFace(pp.first)) {
                DerivedOnceNet A3 = A2;
                A3.EdgeSubdivision(e);
                res.push_back(A3); // Two edges coming out of the pin.
            }
        }
    }
    
    return res;
}

bool ValidateDerivedOnceNet(const DerivedOnceNet& A) {
    // Check that convex faces have distinct edges and the number
    // of sides belongs to [2, 5].
    for (const vector<Dart>& face: A.GetFaces()) {
        if (!A.is_convex_face(face[0]))
            continue;
        if (ssize(face) < 2 or ssize(face) > 5)
            return false;
        set<int> rev_sides;
        for (const Dart& d: face) {
            if (rev_sides.count(d.id))
                return false;
            rev_sides.insert(A.Rev(d).id);
        }
    }
    
    // Check that any face has at most 2 obtuse angles.
    int total_area = 0;
    for (const vector<Dart>& face: A.GetFaces()) {
        int cnt = 0;
        for (const Dart& d: face) {
            if (A.is_obtuse(d)) 
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
        const Dart& e = A.Rev(d);
        if (!A.is_convex_face(d) or !A.is_convex_face(e))
            continue;
        int d_sides = A.GetFace(d).size();
        int e_sides = A.GetFace(e).size();
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
        if (A.is_obtuse(d) and A.is_obtuse(A.Phi(d))
            and A.GetFace(d).size() != 3)
            return false;
    }
    
    return true;
}

bool ContainsBiangle(const DerivedOnceNet& A) {
    int cnt = 0;
    for (const vector<Dart>& face: A.GetFaces()) {
        if (!A.is_convex_face(face[0]))
            continue;
        cnt += ssize(face) == 2;
    }
    return cnt > 0;
}

bool ContainsTriangle(const DerivedOnceNet& A) {
    int cnt = 0;
    for (const vector<Dart>& face: A.GetFaces()) {
        if (!A.is_convex_face(face[0]))
            continue;
        cnt += ssize(face) == 3;
    }
    return cnt > 0;
}

int main() {
    vector<CombinatorialMap> derived_twice = 
        ParseDerivedTwiceNetCandidates();
    cout << "There are " << derived_twice.size()
         << " derived-twice net candidates." << endl << endl;
    
    CollectionUpToIsomorphism<DerivedOnceNet> derived_once;
    
    auto validate_and_add = [&](const DerivedOnceNet& A) {
        if (ValidateDerivedOnceNet(A))
            derived_once.Add(A);
    };
    
    int it = 1;
    for (const CombinatorialMap& A: derived_twice) {
        //~ validate_and_add(A);
        vector<DerivedOnceNet> pin1 = AddPin(A);
        for (const DerivedOnceNet& A2: pin1) {
            vector<DerivedOnceNet> pin2 = AddPin(A2);
            for (const DerivedOnceNet& A3: pin2)
                validate_and_add(A3);
        }
        if (it % 10 == 0) {
            cout << "Processed successfully " << it 
                 << " derived-twice net candidates, which generated "
                 << derived_once.GetList().size() 
                 << " derived-once net candidates." << endl;
        }
        it++;
    }
    cout << endl;
    
    const vector<DerivedOnceNet> derived_once_vec = derived_once.GetList();
    
    int biangle = 0, triangle = 0;
    for (const DerivedOnceNet& A: derived_once_vec) {
        biangle += ContainsBiangle(A);
        triangle += ContainsTriangle(A);
    }
    cout << biangle << " " << triangle << endl;
    //~ return 0;
    
    cout << "We found " << derived_once_vec.size() 
         << " derived-once nets candidates." << endl;
         
    ofstream computer("derived_once_computer.txt");
    ofstream human("derived_once_human.txt");
    computer << derived_once_vec.size() << endl;
    for (int _ = 0; _ < ssize(derived_once_vec); _++) {
        const DerivedOnceNet& A = derived_once_vec[_];
        A.Serialize(computer);
        
        human << "Candidate " << _+1 << endl;
        human << A << endl;
        
        if (ContainsBiangle(A)) {
            cout << _+1 << " ";
        }
    }
}
