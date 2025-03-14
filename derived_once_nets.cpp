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
    friend DerivedOnceNet parseDerivedOnceNet(istream& in);
    
    pair<Dart, Dart> addEdge(Dart bef, Dart aft) = delete;
    pair<Dart, Dart> addVertex(Dart bef, Dart aft) = delete;
    
    pair<Dart, Dart> edgeSubdivision(Dart d) {
        bool d_obtuse = _obtuse_angles.count(d.id);
        bool d_rev_obtuse = _obtuse_angles.count(rev(d).id);
        
        _obtuse_angles.erase(d.id);
        _obtuse_angles.erase(rev(d).id);
        
        pair<Dart,Dart> pp = CombinatorialMap::edgeSubdivision(d);
        
        // The newly created angle is obtuse.
        _obtuse_angles.insert(pp.first.id);
        
        if (d_obtuse)
            _obtuse_angles.insert(pp.second.id);
        if (d_rev_obtuse)
            _obtuse_angles.insert(rev(pp.first).id);
        
        return pp;
    }
    
    bool isObtuse(const Dart& d) const { return _obtuse_angles.count(d.id); }
    
    vector<Dart> obtuseAngles() const {
        vector<Dart> res;
        for (int i: _obtuse_angles)
            res.push_back(dart(i));
        return res;
    }
    
    bool isConvexFace(const Dart& d) const {
        vector<Dart> face = getFace(d);
        for (const Dart& e: face) {
            if (isObtuse(e))
                return false;
        }
        return true;
    }
        
    ULL isomorphismSignature() const {
        const auto& mat = FastLargeMatrixPower(matrixSignature());
        
        ULL res = 0;
        for (const Dart& d: darts()) {
            res += mat[d.id][next(d).id]^mat[d.id][rev(d).id];
        }
        
        for (const Dart& d: obtuseAngles()) {
            ULL tmp = 0;
            for (Dart e: darts())
                tmp ^= mat[d.id][e.id];
            res += tmp;
        }
        return res;
    }
    
    void serialize(ostream& out) const {
        CombinatorialMap::serialize(out);
        
        out << _obtuse_angles.size() << "\n";
        for (int id: _obtuse_angles)
            out << id << " ";
        out << "\n";
    }
};

ostream& operator<<(ostream& out, const DerivedOnceNet& A) {
    out << static_cast<const CombinatorialMap&>(A) << "\n";
    out << "Obtuse angles: " << "\n";
    
    for (const Dart& d: A.obtuseAngles()) {
        out << d << " " << A.phi(d) << "\n";
    }
    
    return out;
}

// Returns, if it exists, a bijection between the darts that is an
// isomorphism of the DerivedOnceNets.
// The ids of the darts are used for a concise representation.
vector<int> findIsomorphism(const DerivedOnceNet& A,
                            const DerivedOnceNet& B) {
    if (A.V() != B.V() or A.E() != B.E() or A.F() != B.F())
        return {};
        
    if (A.obtuseAngles().size() != B.obtuseAngles().size())
        return {};

    for (Dart d: B.darts()) {
        const vector<int>& iso = findIsomorphism(A, B, A.dart(0), d);
        if (iso.empty())
            continue;
        
        bool works = true;
        for (Dart e: A.obtuseAngles())
            works &= B.isObtuse(B.dart(iso[e.id]));
        
        if (works) 
            return iso;
    }
    return {};
}

DerivedOnceNet parseDerivedOnceNet(istream& in) {
    DerivedOnceNet A(parseCombinatorialMap(in));
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
