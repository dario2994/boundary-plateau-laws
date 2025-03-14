#ifndef PLANAR_MAPS
#define PLANAR_MAPS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "combinatorial_maps.cpp"
#include "utils.cpp"


class PlanarMap: public CombinatorialMap {
  private:
    // An arbitrary dart on the boundary of the unbounded face.
    Dart _boundary_dart;
  public:
    PlanarMap(const CombinatorialMap& A): CombinatorialMap(A) {
        _boundary_dart = dart(0);
    }
    
    // Returns a list of all darts on the boundary of the unbounded face.
    vector<Dart> boundaryDarts() const {
        vector<Dart> vec = {_boundary_dart};
        while (phi(vec.back()) != _boundary_dart)
            vec.push_back(phi(vec.back()));
        
        return vec;
    }
    
    // We need to override the definition of EdgeSubdivision as it
    // may invalidate boundary_dart.
    // TODO: Consider erasing it, as we don't use it.
    //       Though, in principle, this is a function that may be useful.
    pair<Dart, Dart> edgeSubdivision(Dart d) {
        vector<Dart> boundary_darts = boundaryDarts();
        bool on_boundary = false;
        bool on_boundary_rev = false;
        for (const Dart& e: boundary_darts) {
            on_boundary |= d == e;
            on_boundary_rev |= rev(d) == e;
        }
        
        pair<Dart,Dart> pp = CombinatorialMap::edgeSubdivision(d);
        if (on_boundary)
            _boundary_dart = pp.first;
        if (on_boundary_rev)
            _boundary_dart = rev(pp.first);
        
        return pp;
    }
    
    void setBoundaryDart(Dart d) { _boundary_dart = d; }
    const Dart& getBoundaryDart() const { return _boundary_dart; }
    
    ULL isomorphismSignature() const {
        const auto& mat = FastLargeMatrixPower(matrixSignature());
        
        ULL res = 0;
        for (const Dart& d: darts()) {
            res += mat[d.id][next(d).id]^mat[d.id][rev(d).id];
        }
        
        vector<Dart> boundary = boundaryDarts();
        for (Dart d: boundary) {
            ULL tmp = 0;
            for (Dart e: darts())
                tmp ^= mat[d.id][e.id];
            res += tmp;
        }
        return res;
    }
};

ostream& operator<<(ostream& out, const PlanarMap& A) {
    out << "(" << static_cast<const CombinatorialMap&>(A) << ", boundary = [";
    for (Dart d: A.boundaryDarts())
        out << d << ", ";
    out << "])";
    return out;
}


vector<int> findIsomorphism(const PlanarMap& A,
                            const PlanarMap& B) {
    if (A.V() != B.V() or A.E() != B.E() or A.F() != B.F())
        return {};
    vector<Dart> boundary = B.boundaryDarts();
    for (Dart d: boundary) {
        const vector<int>& iso = findIsomorphism(A, B, A.getBoundaryDart(), d);
        if (!iso.empty())
            return iso;
    }
    return {};
}


#endif
