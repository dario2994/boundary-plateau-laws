#ifndef UTILS
#define UTILS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;


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


vector<vector<ULL>> FastLargeMatrixPower(vector<vector<ULL>> mat) {
    for (int _ = 0; _ < 8; _++) mat = MatrixMultiplication(mat, mat);
    return mat;
}

// Class to maintain a collection of objects without duplicates, when
// to isomorphic objects are considered equal.
// The underlying class T must have: IsomorphismSignature and
// AreIsomorphic.
template <typename T>
struct CollectionUpToIsomorphism {
    map<ULL, vector<T>> signature2list;
    
    void Add(const T& x) {
        bool is_new = true;
        vector<T>& vec = signature2list[x.IsomorphismSignature()];
        for (const T& y: vec) {
            if (AreIsomorphic(x, y)) {
                is_new = false;
                break;
            }
        }
        if (is_new)
            vec.push_back(x);
        // The following assert is a safety measure to check that
        // the signature used is efficient.
        assert(ssize(vec) < 10);
    }
    
    vector<T> GetList() {
        vector<T> res;
        for (const auto& pp: signature2list) for (const T& x: pp.second) {
            res.push_back(x);
        }
        return res;
    }
};


#endif
