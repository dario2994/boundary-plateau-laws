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


bool criterionBiangle(int id, const DerivedOnceNet& A) {
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

bool criterionSimpleBoundary(int id, const DerivedOnceNet& A) {
    vector<vector<Dart>> faces = A.getFaces();
    
    for (const Dart& d: A.darts()) {
        if (sameFace(A, d, A.rev(d)))
            return false;
    }
    return true;
}

bool criterionConcaveAnglesFar(int id, const DerivedOnceNet& A) {
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



int main() {
    vector<DerivedOnceNet> derived_once = 
        parseDerivedOnceNets();
    cout << "There are " << derived_once.size()
         << " derived-once net candidates." << endl << endl;
         
    // The order of the criteria matters, as a criterion may implicitly
    // assume that all the previous ones pass.
    vector<pair<string, function<bool(int,const DerivedOnceNet&)>>> criteria = {  
        {"biangle\t\t", criterionBiangle},
        {"simple boundary\t", criterionSimpleBoundary},
        {"concave angles are far", criterionConcaveAnglesFar}
    };
    
    map<string, vector<int>> excluded_by; 
    vector<int> survivors;
    
    for (int id = 0; id < ssize(derived_once); id++) {
        const DerivedOnceNet& A = derived_once[id];
        bool survives = true;
        for (auto pp: criteria) {
            string name = pp.first;
            auto criterion = pp.second;
            if (!criterion(id, A)) {
                survives = false;
                excluded_by[name].push_back(id);
                break;
            }
        }
        if (survives) 
            survivors.push_back(id);
    }
    
    cout << "Criterion\t\t\t# Excluded" << endl;
    for (auto pp: criteria) {
        string name = pp.first;
        cout << name << "\t\t" << excluded_by[name].size() << endl;
    }
    cout << endl;
    
    cout << "# survivors = " << survivors.size() << endl;
    
    ofstream out("data/derived_once_survivors.txt");
    out << survivors.size() << endl;
    for (int x: survivors)
        out << x << " ";
    out << endl;
}
