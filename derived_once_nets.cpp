#ifndef DERIVED_ONCE_NETS
#define DERIVED_ONCE_NETS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "combinatorial_maps.cpp"
#include "utils.cpp"

// Class representing a derived-once net.
// It is a combinatorial map such that all vertices have degree 
// either 2 or 3. 
// Moreover, for all vertices of degree 2, we know which side has
// angle 120 and which side has angle 240.

class DerivedOnceNet : public CombinatorialMap {
  private:
    // If i belongs to _concave_angles, then the angle between
    // dart(i) and phi(dart(i)) is 240 degrees.
    vector<Dart> _concave_angles;
    
  public:
    DerivedOnceNet(const CombinatorialMap& A): CombinatorialMap(A) {}
    friend DerivedOnceNet parseDerivedOnceNet(istream& in);
    
    pair<Dart, Dart> addEdge(Dart bef, Dart aft) = delete;
    pair<Dart, Dart> addVertex(Dart bef, Dart aft) = delete;
    
    bool isConcaveAngle(const Dart& d) const {
        return find(_concave_angles.begin(), _concave_angles.end(), d) 
               != _concave_angles.end();
    }
    
    pair<Dart, Dart> edgeSubdivision(Dart d) {
        bool d_obtuse = isConcaveAngle(d);
        bool d_rev_obtuse = isConcaveAngle(rev(d));
        
        _concave_angles.erase(
            remove(_concave_angles.begin(), _concave_angles.end(), d), _concave_angles.end());
        _concave_angles.erase(
            remove(_concave_angles.begin(), _concave_angles.end(), rev(d)), _concave_angles.end());
        
        pair<Dart,Dart> pp = CombinatorialMap::edgeSubdivision(d);
        
        // The newly created angle is obtuse.
        _concave_angles.push_back(pp.first);
        
        if (d_obtuse)
            _concave_angles.push_back(pp.second);
        if (d_rev_obtuse)
            _concave_angles.push_back(rev(pp.first));
        
        return pp;
    }
    
    const vector<Dart>& concaveAngles() const {
        return _concave_angles;
    }
    
    bool isConvexFace(const Dart& d) const {
        vector<Dart> face = CombinatorialMap::face(d);
        for (const Dart& e: face) {
            if (isConcaveAngle(e))
                return false;
        }
        return true;
    }
    
    int faceArea(const Dart& d) const {
        const vector<Dart>& f = face(d);
        int res = 6 - ssize(f);
        for (const Dart& d: f) {
            if (isConcaveAngle(d))
                res += 2;
        }
        return res;
    }
        
    ULL isomorphismSignature() const {
        const auto& mat = largeMatrixPower(matrixSignature());
        
        ULL res = 0;
        for (const Dart& d: darts()) {
            res += mat[d.id][next(d).id]^mat[d.id][rev(d).id];
        }
        
        for (const Dart& d: concaveAngles()) {
            ULL tmp = 0;
            for (Dart e: darts())
                tmp ^= mat[d.id][e.id];
            res += tmp;
        }
        return res;
    }
    
    void serialize(ostream& out) const {
        CombinatorialMap::serialize(out);
        
        out << concaveAngles().size() << "\n";
        for (const Dart& d: concaveAngles())
            out << d.id << " ";
        out << "\n";
    }
};

ostream& operator<<(ostream& out, const DerivedOnceNet& A) {
    out << static_cast<const CombinatorialMap&>(A) << "\n";
    out << "Obtuse angles: " << "\n";
    
    for (const Dart& d: A.concaveAngles()) {
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
        
    if (A.concaveAngles().size() != B.concaveAngles().size())
        return {};

    for (Dart d: B.darts()) {
        const vector<int>& iso = findIsomorphism(A, B, A.dart(0), d);
        if (iso.empty())
            continue;
        
        bool works = true;
        for (Dart e: A.concaveAngles())
            works &= B.isConcaveAngle(B.dart(iso[e.id]));
        
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
        A._concave_angles.push_back(A.dart(id));
    }
    return A;
}


#endif
