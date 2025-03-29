#ifndef METRIZATION
#define METRIZATION

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "derived_once_nets.cpp"
#include "spherical_geometry.cpp"


// Returns a vector label, so that label[d.id] == label[e.id]
// implies that the lengths of darts d and e must be the same.
vector<int> identifyEqualSides(const DerivedOnceNet& A) {
    vector<int> labels;
    
    // Init labels.
    for (Dart d: A.darts()) {
        assert(d.id == labels.size());
        labels.push_back(d.id);
    }
    
    // Cheap implementation of DSU.
    function<int(int)> ancestor = [&](int x) {
        if (labels[x] != x) labels[x] = ancestor(labels[x]);
        return labels[x];
    };
    
    function<void(Dart,Dart)> set_equal = [&](Dart d, Dart e) {
        int x = ancestor(d.id), y = ancestor(e.id);
        labels[x] = y;
    };
    
    function<bool(Dart,Dart)> are_equal = [&](Dart d, Dart e) {
        int x = ancestor(d.id), y = ancestor(e.id);
        return x == y;
    };
    
    // Identification rules
    auto convex_biangle = [&](Dart d) {
        if (A.faceType(d) != pair<int,int>{2, 0})
            return;
        set_equal(d, A.phi(d));
    };
    
    auto convex_triangle = [&](Dart d) {
        if (A.faceType(d) != pair<int,int>{3, 0})
            return;
        set_equal(d, A.phi(d));
        set_equal(d, A.phi(A.phi(d)));
    };
    
    auto convex_quadrilateral = [&](Dart d) {
        if (A.faceType(d) != pair<int,int>{4, 0})
            return;
        set_equal(d, A.phi(A.phi(d)));
        set_equal(A.phi(d), A.phi(A.phi(A.phi(d))));
    };
    
    auto pentagon = [&](Dart d) {
        if (ssize(A.face(d)) != 5)
            return;
        
        // Checking symmetry wrt d of concave angles.
        bool sym1 = A.isConcaveAngle(d) == A.isConcaveAngle(A.phi(d, 4));
        bool sym2 = A.isConcaveAngle(A.phi(d)) == A.isConcaveAngle(A.phi(d, 3));
        if (!sym1 or !sym2)
            return;
        
        // Do something.
        if (are_equal(A.phi(d), A.phi(d, 4)))
            set_equal(A.phi(d, 2), A.phi(d, 3));
        if (are_equal(A.phi(d, 2), A.phi(d, 3)))
            set_equal(A.phi(d), A.phi(d, 4));
    };
    
    // Finalize labels.
    for (Dart d: A.darts())
        ancestor(d.id);
    
    return labels;
}


#endif
