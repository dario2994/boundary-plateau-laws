#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "derived_once_nets.cpp"


vector<DerivedOnceNet> ParseDerivedOnceNets() {
    ifstream in("data/derived_once_computer.txt");
    int cnt;
    in >> cnt;
    vector<DerivedOnceNet> res;
    for (int i = 0; i < cnt; i++) {
        DerivedOnceNet A = ParseDerivedOnceNet(in);
        A.CheckConsistency();
        res.push_back(A);
    }
    return res;
}

bool ContainsBiangle(const DerivedOnceNet& A) {
    int cnt = 0;
    for (const vector<Dart>& face: A.GetFaces()) {
        if (!A.is_convex_face(face[0]))
            continue;
        cnt += ssize(face) == 2;
    }
    return cnt > 0;
}

bool ContainsTriangle(const DerivedOnceNet& A) {
    int cnt = 0;
    for (const vector<Dart>& face: A.GetFaces()) {
        if (!A.is_convex_face(face[0]))
            continue;
        cnt += ssize(face) == 3;
    }
    return cnt > 0;
}

int main() {
    vector<DerivedOnceNet> derived_once = 
        ParseDerivedOnceNets();
    cout << "There are " << derived_once.size()
         << " derived-once net candidates." << endl << endl;
        
    int biangle = 0, triangle = 0;
    for (const DerivedOnceNet& A: derived_once) {
        biangle += ContainsBiangle(A);
        triangle += ContainsTriangle(A);
    }
    cout << biangle << " " << triangle << endl;
}
