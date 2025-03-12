#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "isomorphisms.cpp"

class DerivedOnceNets : public CombinatorialMap {
  private:
    // If i belongs to obtuse_angles, then the angle between
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
    
    bool is_obtuse(Dart d) const { return _obtuse_angles.count(d.id); }
    vector<Dart> obtuse_angles() const {
        vector<Dart> res;
        for (int i: _obtuse_angles)
            res.push_back(dart(i));
        return res;
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

int main() {
    CombinatorialMap A = OneEdgeConstructor();
    A.AddEdge(A.dart(0), A.dart(1));
    A.AddVertex(A.dart(2));
    A.CheckConsistency();
    A.EdgeSubdivision(A.dart(0));
    A.EdgeSubdivision(A.dart(1));
    A.EdgeSubdivision(A.dart(2));
    A.EdgeSubdivision(A.dart(3));
    A.CheckConsistency();
    
    ParseDerivedTwiceNetCandidates();
}
