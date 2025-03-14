#ifndef DERIVED_ONCE_NETS
#define DERIVED_ONCE_NETS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "combinatorial_maps.cpp"
#include "utils.cpp"

class DerivedOnceNet : public CombinatorialMap {
  private:
    // If i belongs to _obtuse_angles, then the angle between
    // dart(i) and Phi(dart(i)) is 240 degrees.
    set<int> _obtuse_angles;
    
  public:
    DerivedOnceNet(const CombinatorialMap& A): CombinatorialMap(A) {}
    friend DerivedOnceNet ParseDerivedOnceNet(istream& in);
    
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
        
    ULL IsomorphismSignature() const {
        const auto& mat = FastLargeMatrixPower(MatrixSignature());
        
        ULL res = 0;
        for (const Dart& d: darts()) {
            res += mat[d.id][Next(d).id]^mat[d.id][Rev(d).id];
        }
        
        for (const Dart& d: obtuse_angles()) {
            ULL tmp = 0;
            for (Dart e: darts())
                tmp ^= mat[d.id][e.id];
            res += tmp;
        }
        return res;
    }
    
    void Serialize(ostream& out) const {
        CombinatorialMap::Serialize(out);
        
        out << _obtuse_angles.size() << "\n";
        for (int id: _obtuse_angles)
            out << id << " ";
        out << "\n";
    }
};

ostream& operator<<(ostream& out, const DerivedOnceNet& A) {
    out << static_cast<const CombinatorialMap&>(A) << "\n";
    out << "Obtuse angles: " << "\n";
    
    vector<Dart> obtuse_angles;
    for (const Dart& d: A.darts()) {
        if (A.is_obtuse(d))
            obtuse_angles.push_back(d);
    }
    
    for (const Dart& d: obtuse_angles) {
        out << d << " " << A.Phi(d) << "\n";
    }
    
    return out;
}

// Returns, if it exists, a bijection between the darts that is an
// isomorphism of the DerivedOnceNets.
// The ids of the darts are used for a concise representation.
vector<int> FindIsomorphism(const DerivedOnceNet& A,
                            const DerivedOnceNet& B) {
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

DerivedOnceNet ParseDerivedOnceNet(istream& in) {
    DerivedOnceNet A(ParseCombinatorialMap(in));
    int obtuse_angles_num;
    in >> obtuse_angles_num;
    
    for (int i = 0; i < obtuse_angles_num; i++) {
        int id;
        in >> id;
        A._obtuse_angles.insert(id);
    }
    return A;
}


#endif
