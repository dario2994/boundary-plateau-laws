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
    friend DerivedOnceNet symmetrize(const DerivedOnceNet&);
    
    pair<Dart, Dart> addEdge(Dart bef, Dart aft) = delete;
    pair<Dart, Dart> addVertex(Dart bef, Dart aft) = delete;
    
    bool isConcaveAngle(const Dart& d) const {
        return find(_concave_angles.begin(), _concave_angles.end(), d) 
               != _concave_angles.end();
    }
    
    int concaveAngleId(const Dart& d) const {
        auto it = find(_concave_angles.begin(), _concave_angles.end(), d);
        assert(it != _concave_angles.end());
        return it - _concave_angles.begin();
    }
    
    double phiAngle(const Dart& d) const {
        return isConcaveAngle(d) ? M_PI/3 : -M_PI/3;
    }
    
    double nextAngle(const Dart& d) const {
        return isConcaveAngle(rev(d)) ? 4*M_PI/3 : 2*M_PI/3;
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
        for (const Dart& e: f) {
            if (isConcaveAngle(e))
                res += 2;
        }
        return res;
    }
    
    // Returns {# sides, # concave angles} of the face of d.
    pair<int,int> faceType(const Dart& d) const {
        const vector<Dart>& f = face(d);
        int concave_num = 0;
        for (const Dart& e: f) concave_num += isConcaveAngle(e);
        
        return {ssize(f), concave_num};
    }
    
    // Returns the angles (in radiants) of the face, starting from the
    // one just after d.
    vector<double> faceAngles(const Dart& d) const {
        vector<double> angles;
        for (const Dart& e: face(d)) {
            double angle = isConcaveAngle(e) ? 4*M_PI/3 : 2*M_PI/3;
            angles.push_back(angle);
        }
        return angles;
    }
    
    // Returns the angles (in degrees) of the face, starting from the
    // one just after d.
    vector<int> faceAnglesDeg(const Dart& d) const {
        vector<int> angles;
        for (const Dart& e: face(d)) {
            int angle = isConcaveAngle(e) ? 240 : 120;
            angles.push_back(angle);
        }
        return angles;
    }
    
    // Whether d and e belong to the same face.
    bool sameFace(const Dart& d, const Dart& e) const {
        const vector<Dart>& f = face(d);
        for (const Dart& d1: f)
            if (d1 == e)
                return true;
        return false;
    }
    
    
    bool isConvexFace(const vector<Dart>& f) const {
        return isConvexFace(f[0]);
    }
    
    int faceArea(const vector<Dart>& f) const {
        return faceArea(f[0]);
    }
    
    pair<int,int> faceType(const vector<Dart>& f) const {
        return faceType(f[0]);
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

DerivedOnceNet symmetrize(const DerivedOnceNet& A) {
    vector<Dart> sym_concave_angles;
    for (const Dart& d: A._concave_angles)
        sym_concave_angles.push_back(A.rev(A.phi(d)));
    
    DerivedOnceNet B(symmetrize(static_cast<CombinatorialMap>(A)));
    B._concave_angles = sym_concave_angles;
    return B;
}


#endif
