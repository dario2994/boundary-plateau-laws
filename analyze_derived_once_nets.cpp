#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "derived_once_nets.cpp"
#include "metrization.cpp"


int countObtuseAngles(const DerivedOnceNet& A, Dart d) {
    int res = 0;
    for (Dart e: A.face(d)) {
        res += A.isConcaveAngle(e);
    }
    return res;
}

int sumAreaConvexFaces(const DerivedOnceNet& A) {
    int res = 0;
    for (auto face: A.faces()) {
        if (A.isConvexFace(face))
            res += A.faceArea(face);
    }
    return res;
}

int pinContribution(const DerivedOnceNet& A, Dart d) {
    if (A.faceType(d).second != 2)
        return 0;
    //~ assert(A.degree(d.from) == 3 and A.degree(d.to) == 3);
    int dist = -1;
    Dart curr = d;
    while (1) {
        if (A.isConcaveAngle(curr)) {
            if (dist == -1) dist = 1;
            else break;
        } else if (dist >= 0) dist++;
        
        curr = A.phi(curr);
    }
    int contribution = max(4-dist, 0);
    return contribution;
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

bool containsTriangleQuadrilateral(const DerivedOnceNet& A) {
    for (Dart d:A.darts()) {
        if (A.faceType(d)==pair<int,int> {3,0} and
            A.faceType(A.rev(d))==pair<int,int> {4,0}) 
            return true;
    }
    return false;       
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

bool containsPentagonTwoObtuse(const DerivedOnceNet& A) {
    for (auto face: A.faces()) {
        if (A.faceType(face[0])==pair<int,int> {5,2})
            return true;
    }
    return false;    
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

bool containsPentagonOneObtuse(const DerivedOnceNet& A) {
    for (auto face: A.faces()) {
        if (A.faceType(face[0])==pair<int,int> {5,1})
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

bool containsAdjacentQuadrilaterals(const DerivedOnceNet& A) {
    pair <int,int> quad={4,0};
    for (Dart d: A.darts()) {
        if ((A.faceType(d)==quad) 
        and (A.faceType(A.rev(d))==quad)) return true;
    }
    return false;
}

bool containsDoubleBlade(const DerivedOnceNet& A) {
    pair <int,int> quad={4,0};
    pair <int,int> pent={5,0};
    for (Dart d: A.darts()) {
        if ((A.faceType(d)!=quad) 
            or (A.faceType(A.rev(d))!=quad)) continue;
        if ((A.faceType(A.rev(A.phi(d)))==pent) and
            (A.faceType(A.rev(A.phi(A.rev(d))))==pent)) 
            return true;
    }
    return false;
}

bool findPentagonalPrism(const DerivedOnceNet& A) {
    for (auto face:A.faces()) {
        Dart d=face[0];
        if ((A.faceType(d)!= pair<int,int> {7,2}) and
            (A.faceType(d)!= pair<int,int> {4,0})) return false;
    }    
    return true;
}

bool findDodecahedron(const DerivedOnceNet& A) {
    for (auto face:A.faces()) {
        Dart d=face[0];
        if ((A.faceType(d)!= pair<int,int> {8,2}) and
            (A.faceType(d)!= pair<int,int> {5,0})) return false;
    }    
    return true;
}

bool findHeppes7(const DerivedOnceNet& A) {
    for (auto face:A.faces()) {
        Dart d=face[0];
        if ((A.faceType(d)!= pair<int,int> {7,2}) and
            (A.faceType(d)!= pair<int,int> {5,0})) return false;
    }    
    return true;
}

bool containsTripleQuadrilateral(const DerivedOnceNet &A) {
    pair<int,int> quad = {4,0}; 
    for (Dart d:A.darts()) {
        if ((A.faceType(d)!=quad) or (A.faceType(A.rev(d))!= quad))
            continue;
        if (A.faceType(A.rev(A.phi(A.phi(d))))==quad) return true;
    }
    return false;
}

bool findHeppes5(const DerivedOnceNet& A) {
    for (auto face:A.faces()) {
        Dart d=face[0];
        if ((A.faceType(d)!= pair<int,int> {6,2}) and
            (A.faceType(d)!= pair<int,int> {5,0})) return false;
    }    
    return true;
}

bool findFalange(const DerivedOnceNet& A) {
    if (containsTripleQuadrilateral(A)==false) return false;
    for (Dart d: A.darts()) {
        if ((A.faceType(d)!=pair<int,int> {4,0}) and
            (A.faceType(d)!=pair<int,int> {5,1}) and
            (A.faceType(d)!=pair<int,int> {6,1}) and
            (A.faceType(d)!=pair<int,int> {7,1})) return false;
    }
    return true;
}

bool containsPentagon(const DerivedOnceNet& A) {
    for (Dart d:A.darts()) {
        if (A.faceType(d)==pair <int,int> {5,0}) return true;
    }
    return false;
}

bool analyzeconjectureHeptagon(const DerivedOnceNet& A) {
    for (Dart d:A.darts()) {
        if ((A.faceType(d)!=pair <int,int> {7,2}) or
            (A.faceType(A.rev(d))!=pair<int,int> {5,0}))
            continue;
        Dart e=A.phi(A.phi(d));
        Dart f=A.phi(A.phi(e));
        if ((A.isConcaveAngle(e)==true) and (A.isConcaveAngle(f)==true)) return true;
    }
    return false;
}

int main() {
    vector<DerivedOnceNet> derived_once = 
        parseDerivedOnceNets("data/derived_once_computer.txt");
        
    ifstream in("data/derived_once_survivors.txt");
    int cnt_survivors;
    in >> cnt_survivors;
    vector<int> survivors(cnt_survivors);
    for (int i = 0; i < cnt_survivors; i++)
        in >> survivors[i];
    in.close();
    
    cout << "There are " << derived_once.size()
         << " derived-once net candidates." << endl;
    cout << "There are " << survivors.size()
         << " survivors." << endl;
    
    //~ MetricDerivedOnceNet mdon(derived_once[1476]);
    //~ cout << mdon.isConsistent() << endl;
    
    //~ return 0;        
    
    int cnt_interesting = 0;
    for (int id: survivors) {
        const DerivedOnceNet& A = derived_once[id];
        if (containsTriangle(A)==false) continue;
        cout << id << endl;
        cnt_interesting++;
    }
    cout << endl;
    cout << "Ci sono "<< cnt_interesting << " reti con triangolo" << endl;
}

