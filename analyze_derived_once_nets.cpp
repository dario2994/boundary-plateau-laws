#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "derived_once_nets.cpp"


vector<DerivedOnceNet> parseDerivedOnceNets() {
    ifstream in("data/derived_once_computer.txt");
    int cnt;
    in >> cnt;
    vector<DerivedOnceNet> res;
    for (int i = 0; i < cnt; i++) {
        DerivedOnceNet A = parseDerivedOnceNet(in);
        A.checkConsistency();
        res.push_back(A);
    }
    return res;
}

bool onePin(const DerivedOnceNet& A) {
    int cnt = 0;
    for (const vector<Dart>& face: A.faces()) {
        cnt += !A.isConvexFace(face[0]);
    }
    return cnt == 1;
}

bool containsBiangle(const DerivedOnceNet& A) {
    int cnt = 0;
    for (const vector<Dart>& face: A.faces()) {
        if (!A.isConvexFace(face[0]))
            continue;
        cnt += ssize(face) == 2;
    }
    return cnt > 0;
}

bool containsTriangle(const DerivedOnceNet& A) {
    int cnt = 0;
    for (const vector<Dart>& face: A.faces()) {
        if (!A.isConvexFace(face[0]))
            continue;
        cnt += ssize(face) == 3;
    }
    return cnt > 0;
}

int faceArea(const DerivedOnceNet& A, Dart d) {
    const vector<Dart>& face = A.face(d);
    int res = 6 - ssize(face);
    for (const Dart& d: face) {
        if (A.isObtuse(d))
            res += 2;
    }
    return res;
}

bool sameFace(const DerivedOnceNet& A, Dart d, Dart e) {
    const vector<Dart>& face = A.face(d);
    for (const Dart& e1: face) {
        if (e1 == e)
            return true;
    }
    return false;
}


bool isPesceMagro(int id, const DerivedOnceNet& A) {
    auto faces = A.faces();
    vector<pair<int,bool>> types;
    for (auto face: faces) {
        types.push_back({ssize(face), A.isConvexFace(face[0])});
    }
    sort(types.begin(), types.end());
    
    vector<pair<int,bool>> hope = {{2, 1}, {4, 1}, {5, 0}, {5, 0}};
    return types == hope;
}



int main() {
    vector<DerivedOnceNet> derived_once = 
        parseDerivedOnceNets();
    cout << "There are " << derived_once.size()
         << " derived-once net candidates." << endl << endl;
        
    
    int with_biangle = 0;
    for (int i = 0; i < ssize(derived_once); i++) {
        const DerivedOnceNet& A = derived_once[i];
        with_biangle += containsBiangle(A);
        if (isPesceMagro(i, A)) {
            cout << i << " ";
        }
    }
    cout << endl;
    cout << with_biangle << endl;
}
