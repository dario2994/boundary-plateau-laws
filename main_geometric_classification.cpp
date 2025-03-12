#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "isomorphisms.cpp"

vector<CombinatorialMap> ParseDerivedTwiceNetCandidates() {
    ifstream in("list_combinatorial_maps_computer.txt");
    int cnt;
    in >> cnt;
    vector<CombinatorialMap> res;
    for (int i = 0; i < cnt; i++) {
        CombinatorialMap A = ParserConstructor(in);
        A.CheckConsistency();
        res.push_back(A);
    }
    return res;
}

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
    
    ParseDerivedTwiceNetCandidates();
}
