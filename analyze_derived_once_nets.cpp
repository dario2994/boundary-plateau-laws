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
    for (const vector<Dart>& face: A.getFaces()) {
        cnt += !A.isConvexFace(face[0]);
    }
    return cnt == 1;
}

bool containsBiangle(const DerivedOnceNet& A) {
    int cnt = 0;
    for (const vector<Dart>& face: A.getFaces()) {
        if (!A.isConvexFace(face[0]))
            continue;
        cnt += ssize(face) == 2;
    }
    return cnt > 0;
}

bool containsTriangle(const DerivedOnceNet& A) {
    int cnt = 0;
    for (const vector<Dart>& face: A.getFaces()) {
        if (!A.isConvexFace(face[0]))
            continue;
        cnt += ssize(face) == 3;
    }
    return cnt > 0;
}

int faceArea(const DerivedOnceNet& A, Dart d) {
    const vector<Dart>& face = A.getFace(d);
    int res = 6 - ssize(face);
    for (const Dart& d: face) {
        if (A.isObtuse(d))
            res += 2;
    }
    return res;
}

bool sameFace(const DerivedOnceNet& A, Dart d, Dart e) {
    const vector<Dart>& face = A.getFace(d);
    for (const Dart& e1: face) {
        if (e1 == e)
            return true;
    }
    return false;
}


bool criterion_biangle(int id, const DerivedOnceNet& A) {
    vector<vector<Dart>> faces = A.getFaces();
    
    for (const vector<Dart>& face: faces) {
        if (ssize(face) != 2)
            continue;
        
        if (!A.isConvexFace(face[0]))
            return false;
            
        if (sameFace(A, A.rev(face[0]), A.rev(face[1])))
            return false;
            
        if (A.isConvexFace(A.rev(face[0])) 
            or A.isConvexFace(A.rev(face[1])))
            return false;
            
        if (faceArea(A, A.rev(face[0])) + faceArea(A, A.rev(face[1])) < 3)
            return false;
            
        for (const Dart& d: A.darts()) {
            if (A.isConvexFace(d) and A.isConvexFace(A.rev(d))
                and ssize(A.getFace(d)) == 3
                and ssize(A.getFace(A.rev(d))) == 4)
                return false;
            if (d != face[0] and d != face[1] and ssize(A.getFace(d)) == 2)
                return false;
        }
    }
    return true;
}

bool criterion_simple_boundary(int id, const DerivedOnceNet& A) {
    vector<vector<Dart>> faces = A.getFaces();
    
    for (const Dart& d: A.darts()) {
        if (sameFace(A, d, A.rev(d)))
            return false;
    }
    return true;
}

bool criterion_concave_angles_far(int id, const DerivedOnceNet& A) {
    for (const vector<Dart>& face: A.getFaces()) {
        for (const Dart& d: face) {
            if (!A.isObtuse(d)) continue;
            
            int dist = 1;
            Dart curr = A.phi(d);
            while (curr != d) {
                if (A.isObtuse(curr)) {
                    assert(dist != 1);
                    if (5 < ssize(face) - dist)
                        return false;
                }
                curr = A.phi(curr);
                dist++;
            }
        }
    }
    return true;
}

bool isPesceMagro(int id, const DerivedOnceNet& A) {
    auto faces = A.getFaces();
    vector<pair<int,bool>> types;
    for (auto face: faces) {
        types.push_back({ssize(face), A.isConvexFace(face[0])});
    }
    sort(types.begin(), types.end());
    // 
    vector<pair<int,bool>> hope = {{2, 1}, {4, 1}, {5, 0}, {5, 0}};
    return types == hope;
}


int main() {
    vector<DerivedOnceNet> derived_once = 
        parseDerivedOnceNets();
    cout << "There are " << derived_once.size()
         << " derived-once net candidates." << endl << endl;
        
    int biangle = 0;
    int excluded = 0;
    int non_s1 = 0;
    int bomba = 0;
    for (int i = 0; i < ssize(derived_once); i++) {
        const DerivedOnceNet& A = derived_once[i];
        biangle += containsBiangle(A);
        if (!criterion_biangle(i, A))
            excluded++;
        else if (!criterion_simple_boundary(i, A))
            non_s1++;
        else if (!criterion_concave_angles_far(i, A))
            bomba++;
        else if (isPesceMagro(i, A)) {
            cout << i << " ";
        }
        //~ if (onePin(A)) {
            //~ cout << "One pin: " << i << endl;
        //~ }
    }
    cout << endl;
    cout << excluded << " " << non_s1 << " " << bomba << endl;
}
