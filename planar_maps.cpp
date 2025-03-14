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
    Dart boundary_dart;
  public:
    PlanarMap(const CombinatorialMap& A): CombinatorialMap(A) {
        boundary_dart = dart(0);
    }
    
    // Returns a list of all darts on the boundary of the unbounded face.
    vector<Dart> BoundaryDarts() const {
        vector<Dart> vec = {boundary_dart};
        while (Phi(vec.back()) != boundary_dart)
            vec.push_back(Phi(vec.back()));
        
        return vec;
    }
    
    // We need to override the definition of EdgeSubdivision as it
    // may invalidate boundary_dart.
    // TODO: Consider erasing it, as we don't use it.
    //       Though, in principle, this is a function that may be useful.
    pair<Dart, Dart> EdgeSubdivision(Dart d) {
        vector<Dart> boundary_darts = BoundaryDarts();
        bool on_boundary = false;
        bool on_boundary_rev = false;
        for (const Dart& e: boundary_darts) {
            on_boundary |= d == e;
            on_boundary_rev |= Rev(d) == e;
        }
        
        pair<Dart,Dart> pp = CombinatorialMap::EdgeSubdivision(d);
        if (on_boundary)
            boundary_dart = pp.first;
        if (on_boundary_rev)
            boundary_dart = Rev(pp.first);
        
        return pp;
    }
    
    void SetBoundaryDart(Dart d) {
        boundary_dart = d;
    }
    const Dart& GetBoundaryDart() const { return boundary_dart; }
    
    ULL IsomorphismSignature() const {
        const auto& mat = FastLargeMatrixPower(MatrixSignature());
        
        ULL res = 0;
        for (const Dart& d: darts()) {
            res += mat[d.id][Next(d).id]^mat[d.id][Rev(d).id];
        }
        
        vector<Dart> boundary = BoundaryDarts();
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
    for (Dart d: A.BoundaryDarts())
        out << d << ", ";
    out << "])";
    return out;
}


vector<int> FindIsomorphism(const PlanarMap& A,
                            const PlanarMap& B) {
    if (A.V() != B.V() or A.E() != B.E() or A.F() != B.F())
        return {};
    vector<Dart> boundary = B.BoundaryDarts();
    for (Dart d: boundary) {
        const vector<int>& iso = FindIsomorphism(A, B, A.GetBoundaryDart(), d);
        if (!iso.empty())
            return iso;
    }
    return {};
}


#endif
