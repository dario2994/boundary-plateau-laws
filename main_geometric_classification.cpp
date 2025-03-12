#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "isomorphisms.cpp"

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
}
