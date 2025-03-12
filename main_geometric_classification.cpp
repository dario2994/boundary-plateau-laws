#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "isomorphisms.cpp"

class DerivedOnceNets : public CombinatorialMap {
  private:
    // If i belongs to _obtuse_angles, then the angle between
    // dart(i) and Phi(dart(i)) is 240 degrees.
    set<int> _obtuse_angles;
    
  public:
    DerivedOnceNets(const CombinatorialMap& A): CombinatorialMap(A) {}
    
    pair<Dart, Dart> AddEdge(Dart bef, Dart aft) = delete;
    pair<Dart, Dart> AddVertex(Dart bef, Dart aft) = delete;
    
    pair<Dart, Dart> EdgeSubdivision(Dart d) {
        bool d_obtuse = _obtuse_angles.count(d.id);
        bool d_rev_obtuse = _obtuse_angles.count(Rev(d).id);
        
        _obtuse_angles.erase(d.id);
        _obtuse_angles.erase(Rev(d).id);
        
        pair<Dart,Dart> pp = CombinatorialMap::EdgeSubdivision(d);
        
        // The newly created angle is obtuse.
        _obtuse_angles.insert(pp.first.id);
        
        if (d_obtuse)
            _obtuse_angles.insert(pp.second.id);
        if (d_rev_obtuse)
            _obtuse_angles.insert(Rev(pp.first).id);
        
        return pp;
    }
    
    bool is_obtuse(const Dart& d) const { return _obtuse_angles.count(d.id); }
    vector<Dart> obtuse_angles() const {
        vector<Dart> res;
        for (int i: _obtuse_angles)
            res.push_back(dart(i));
        return res;
    }
    
    bool is_convex_face(const Dart& d) const {
        vector<Dart> face = GetFace(d);
        for (const Dart& e: face) {
            if (is_obtuse(e))
                return false;
        }
        return true;
    }
};

// Returns, if it exists, a bijection between the darts that is an
// isomorphism of the DerivedOnceNets.
// The ids of the darts are used for a concise representation.
vector<int> FindIsomorphism(const DerivedOnceNets& A,
                            const DerivedOnceNets& B) {
    if (A.V() != B.V() or A.E() != B.E() or A.F() != B.F())
        return {};
        
    if (A.obtuse_angles().size() != B.obtuse_angles().size())
        return {};

    for (Dart d: B.darts()) {
        const vector<int>& iso = FindIsomorphism(A, B, A.dart(0), d);
        if (iso.empty())
            continue;
        
        bool works = true;
        for (Dart e: A.obtuse_angles())
            works &= B.is_obtuse(B.dart(iso[e.id]));
        
        if (works) 
            return iso;
    }
    return {};
}

ULL IsomorphismSignature(const DerivedOnceNets& A) {
    const auto& mat = LargePower(MatrixSignature(A));
    
    ULL res = 0;
    for (const Dart& d: A.darts()) {
        res += mat[d.id][A.Next(d).id]^mat[d.id][A.Rev(d).id];
    }
    
    vector<Dart> obtuse_angles = A.obtuse_angles();
    for (Dart d: obtuse_angles) {
        ULL tmp = 0;
        for (Dart e: A.darts())
            tmp ^= mat[d.id][e.id];
        res += tmp;
    }
    return res;
}


vector<CombinatorialMap> ParseDerivedTwiceNetCandidates() {
    ifstream in("list_combinatorial_maps_computer.txt");
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

vector<DerivedOnceNets> AddPin(const DerivedOnceNets& A) {
    vector<DerivedOnceNets> res;
    
    for (const vector<Dart>& face: A.GetFaces()) {
        // We add a pin inside face.
        // A pin may generate one or two concave angles.
        for (const Dart& d: face) {
            DerivedOnceNets A2 = A;
            pair<Dart, Dart> pp = A2.EdgeSubdivision(d);
            res.push_back(A2); // One edge coming out of the pin.
            for (const Dart& e: A2.GetFace(pp.first)) {
                DerivedOnceNets A3 = A2;
                A3.EdgeSubdivision(e);
                res.push_back(A3); // Two edges coming out of the pin.
            }
        }
    }
    
    return res;
}

bool ValidateDerivedOnceNet(const DerivedOnceNets& A) {
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

int main() {
    //~ CombinatorialMap A = OneEdgeConstructor();
    //~ A.AddEdge(A.dart(0), A.dart(1));
    //~ A.AddVertex(A.dart(2));
    //~ A.CheckConsistency();
    //~ A.EdgeSubdivision(A.dart(0));
    //~ A.EdgeSubdivision(A.dart(1));
    //~ A.EdgeSubdivision(A.dart(2));
    //~ A.EdgeSubdivision(A.dart(3));
    //~ A.CheckConsistency();
    
    vector<CombinatorialMap> vec = ParseDerivedTwiceNetCandidates();
    
    MaintainUnique<DerivedOnceNets> nets;
    
    for (const CombinatorialMap& A: vec) {
        vector<DerivedOnceNets> pin1 = AddPin(A);
        for (const DerivedOnceNets& A2: pin1) {
            vector<DerivedOnceNets> pin2 = AddPin(A2);
            for (const DerivedOnceNets& A3: pin2)
                nets.Add(A3);
        }
    }
}
