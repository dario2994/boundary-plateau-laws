#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "derived_once_nets.cpp"


int countObtuseAngles(const DerivedOnceNet& A, Dart d) {
    int res = 0;
    for (Dart e: A.face(d)) {
        res += A.isConcaveAngle(e);
    }
    return res;
}

vector<DerivedOnceNet> parseDerivedOnceNets(string in_file) {
    ifstream in(in_file);
    int cnt;
    in >> cnt;
    vector<DerivedOnceNet> res;
    for (int i = 0; i < cnt; i++) {
        const DerivedOnceNet& A = parseDerivedOnceNet(in);
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


bool sameFace(const DerivedOnceNet& A, Dart d, Dart e) {
    const vector<Dart>& face = A.face(d);
    for (const Dart& e1: face) {
        if (e1 == e)
            return true;
    }
    return false;
}

bool isPesceMagro(const DerivedOnceNet& A) {
    auto faces = A.faces();
    if (ssize(faces)!= 3) return false;
        vector<pair<int,bool>> types;
    for (auto face: faces) {
        if ((A.isConvexFace(face[0])==true) and (ssize(face)!=4)) return false;
        if ((A.isConvexFace(face[0])==false) and (ssize(face)!=3)) return false;
    }
    return true;
}

bool containsConcaveQuadrilateral(const DerivedOnceNet& A) {
    for (auto face: A.faces()) {
        if (ssize(face) == 4 and !A.isConvexFace(face[0]))
            return true;
    }
    return false;
}

bool containsConcaveTriangle(const DerivedOnceNet& A) {
    for (auto face: A.faces()) {
        if (ssize(face) == 3 and !A.isConvexFace(face[0]))
            return true;
    }
    return false;
}

bool containsConcavePentagon(const DerivedOnceNet& A) {
    for (auto face: A.faces()) {
        if (ssize(face) == 5 and !A.isConvexFace(face[0]))
            return true;
    }
    return false;
}

bool containsAdjacentConcaveAngles(const DerivedOnceNet& A) {
    for (const Dart& d: A.darts()) {
        if (A.isConcaveAngle(d) and A.isConcaveAngle(A.phi(d)) 
            and A.face(d).size() != 3)
            return true;
    }
    return false;
}

bool containsHexagonsTwoObtuse(const DerivedOnceNet& A) {
    for (auto face: A.faces()) {
        if (ssize(face) == 6 and countObtuseAngles(A, face[0]) == 2)
            return true;
    }
    return false;
}

int main() {
    vector<DerivedOnceNet> derived_once = 
        parseDerivedOnceNets("data/derived_once_computer.txt");
    int count=0;
    for (auto A: derived_once) {
        if (isPesceMagro(A)==true) cout << count;
        count++; 
    }
}
