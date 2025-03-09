#ifndef ISOMORPHISMS
#define ISOMORPHISMS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"

// Check whether the map iso between the darts of A and B induces
// an isomorphism of the combinatorial map.
// It is sufficient to check that it commutes with Rev and Next.
bool CheckIsomorphism(const CombinatorialMap& A,
                      const CombinatorialMap& B,
                      const vector<int>& iso) {
    int E = A.E();
    assert(A.E() == E and B.E() == E and ssize(iso) == 2*E);
    for (int i = 0; i < 2*E; i++) {
        const Dart& a = A.dart(i);
        const Dart& b = B.dart(iso[i]);
        if (B.Rev(b).id != iso[A.Rev(a).id]) return false;
        if (B.Next(b).id != iso[A.Next(a).id]) return false;
    }
    return true;
}

// Returns, if it exists, an isomorphism between A and B that maps dart_A
// into dart_B.
// The ids of the darts are used for a concise representation.
vector<int> FindIsomorphism(const CombinatorialMap& A,
                            const CombinatorialMap& B,
                            Dart dart_A,
                            Dart dart_B) {
    int E = A.E();
    
    vector<int> ordA = A.OrderedVisit(dart_A);
    vector<int> ordB = B.OrderedVisit(dart_B);
    vector<int> iso(2*E);
    for (int i = 0; i < 2*E; i++) iso[ordA[i]] = ordB[i];
    
    if (CheckIsomorphism(A, B, iso))
        return iso;
    return {};
}


// Returns, if it exists, a bijection between the darts that is an
// isomorphism of the combinatorial maps.
// The ids of the darts are used for a concise representation.
vector<int> FindIsomorphism(const CombinatorialMap& A,
                            const CombinatorialMap& B) {
    if (A.V() != B.V() or A.E() != B.E() or A.F() != B.F())
        return {};

    for (Dart d: B.darts()) {
        const vector<int>& iso = FindIsomorphism(A, B, A.dart(0), d);
        if (!iso.empty())
            return iso;
    }
    return {};
}

vector<int> FindIsomorphism(const PlanarMap& A,
                            const PlanarMap& B) {
    if (A.V() != B.V() or A.E() != B.E() or A.F() != B.F())
        return {};
    vector<Dart> boundary = B.BoundaryEdges();
    for (Dart d: boundary) {
        const vector<int>& iso = FindIsomorphism(A, B, A.boundary, d);
        if (!iso.empty())
            return iso;
    }
    return {};
}

template<typename T>
bool AreIsomorphic(const T& A, const T& B) {
    return !FindIsomorphism(A, B).empty();
}

vector<vector<ULL>> MatrixMultiplication(const vector<vector<ULL>>& A,
                                         const vector<vector<ULL>>& B) {
    int n = A.size();
    assert(ssize(A[0]) == n and ssize(B) == n and ssize(B[0]) == n);
    vector<vector<ULL>> C(n, vector<ULL>(n, 0));
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++) C[i][j] += A[i][k] * B[k][j];
    return C;
}

ULL Trace(const vector<vector<ULL>>& mat) {
    ULL trace = 0;
    for (int i = 0; i < ssize(mat); i++) trace += mat[i][i];
    return trace;
}

// Returns a square matrix 2E x 2E. Let M_A be the matrix.
// If f is a bijection from the darts of A to the darts of B that
// is a isomorphism of combinatorial maps, then M_B is obtained
// permuting the rows/columns of M_A according to f.
vector<vector<ULL>> MatrixSignature(const CombinatorialMap& A) {
    static const ULL p_stay = 424'142'591;
    static const ULL p_rev = 984'689'391;
    static const ULL p_next = 123'456'789;
    
    vector<vector<ULL>> mat(2*A.E(), vector<ULL>(2*A.E(), 0));
    for (Dart d: A.darts()) {
        int i = d.id;
        mat[i][i] += p_stay;
        mat[i][A.Rev(d).id] += p_rev;
        mat[i][A.Next(d).id] += p_next;
    }
    for (int _ = 0; _ < 8; _++) mat = MatrixMultiplication(mat, mat);
    return mat;
}

ULL IsomorphismSignature(const CombinatorialMap& A) {
    return Trace(MatrixSignature(A));
}

ULL IsomorphismSignature(const PlanarMap& A) {
    const auto& mat = MatrixSignature(A);
    ULL res = Trace(mat);
    vector<Dart> boundary = A.BoundaryEdges();
    for (Dart d: boundary) res += 131'657'115 * mat[d.id][d.id];
    return res;
}


// Class to maintain a collection of objects without duplicates.
// Two objects are considered duplicates if they are isomorphic.
// The underlying class T must have: IsomorphismSignature and
// AreIsomorphic.
// Both CombinatorialMap and PlanarMap are valid types for T.
template <typename T>
struct MaintainUnique {
    map<ULL, vector<T>> signature2list;
    
    void Add(const T& x) {
        bool is_new = true;
        vector<T>& vec = signature2list[IsomorphismSignature(x)];
        for (const T& y: vec) {
            if (AreIsomorphic(x, y)) {
                is_new = false;
                break;
            }
        }
        if (is_new)
            vec.push_back(x);
        assert(ssize(vec) < 30);
        //~ if (ssize(vec) > 1) {
            //~ auto A = vec[0];
            //~ auto B = vec[1];
            //~ int E = A.E;
            
            //~ cout << A << endl;
            //~ cout << B << endl;
            
            //~ vector<int> ordA = OrderedVisit(A, 8);
            //~ vector<int> ordB = OrderedVisit(B, 10);
            
            //~ cout << "ordA = ";
            //~ for (int i = 0; i < 2*E; i++)
                //~ cout << ordA[i] << ", ";
            //~ cout << endl;
            //~ cout << "ordB = ";
            //~ for (int i = 0; i < 2*E; i++)
                //~ cout << ordB[i] << ", ";
            //~ cout << endl;
            
            //~ vector<int> iso(2*E);
            //~ for (int i = 0; i < 2*E; i++) iso[ordA[i]] = ordB[i];
    
            //~ cout << "Is iso? " << CheckIsomorphism(A, B, iso) << endl;
            
            //~ exit(0);
        //~ }
    }
    
    vector<T> GetUnique() {
        vector<T> res;
        for (const auto& pp: signature2list) for (const T& x: pp.second) {
            res.push_back(x);
        }
        return res;
    }
};


#endif
