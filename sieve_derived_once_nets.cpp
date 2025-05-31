#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "derived_once_nets.cpp"
//#include "metrization.cpp"


vector<DerivedOnceNet> parseDerivedOnceNets(string in_file) {
    ifstream in(in_file);
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


bool containsBiangle(const DerivedOnceNet& A) {
    for (const vector<Dart>& f: A.faces()) {
        if (A.faceType(f) == pair<int,int>{2, 0})
            return true;
    }
    return false;
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


// Whether the net may contain an edge >= pi.
bool mayHaveEdgesGePi(const DerivedOnceNet& A) {
    if (containsBiangle(A))
        return true;

    for (Dart d: A.darts()) {
        Dart e = A.rev(d);
        if (A.isConvexFace(d) or A.isConvexFace(e))
            continue;
        
        if (A.faceArea(d) + A.faceArea(e) <= 
                6 + pinContribution(A, d) + pinContribution(A, e))
            continue;

        return true;
    }
    return false;
}

// Whether the net may contain an edge > pi.
bool mayHaveEdgesGtPi(const DerivedOnceNet& A) {
    if (containsBiangle(A))
        return false;
    return mayHaveEdgesGePi(A);
}

int sumAreaConvexFaces(const DerivedOnceNet& A) {
    int res = 0;
    for (auto face: A.faces()) {
        if (A.isConvexFace(face))
            res += A.faceArea(face);
    }
    return res;
}

bool criterionBiangle(const DerivedOnceNet& A) {
    vector<vector<Dart>> faces = A.faces();
    
    for (const vector<Dart>& face: faces) {
        if (ssize(face) != 2)
            continue;
        
        if (!A.isConvexFace(face))
            return false;
            
        if (A.sameFace(A.rev(face[0]), A.rev(face[1])))
            return false;
            
        if (A.isConvexFace(A.rev(face[0])) 
            or A.isConvexFace(A.rev(face[1])))
            return false;

        int strange_area = sumAreaConvexFaces(A) - 4;
        strange_area += pinContribution(A, A.rev(face[0]));
        strange_area += pinContribution(A, A.rev(face[1]));
        
        if (strange_area > 5)
            return false;
            
        for (const Dart& d: A.darts()) {
            if (A.faceType(d) == pair<int,int>{3, 0} 
                    and A.faceType(A.rev(d)) == pair<int,int>{4, 0})
                return false;
            if (d != face[0] and d != face[1] and ssize(A.face(d)) == 2)
                return false;
        }
    }
    return true;
}

bool criterionSimpleBoundary(const DerivedOnceNet& A) {
    vector<vector<Dart>> faces = A.faces();
    
    for (const Dart& d: A.darts()) {
        if (A.sameFace(d, A.rev(d)))
            return false;
    }
    return true;
}

bool criterionAdjacentConcaveAngles(const DerivedOnceNet& A) {
    for (const Dart& d: A.darts()) {
        if (A.isConcaveAngle(d) and A.isConcaveAngle(A.phi(d)) 
                and ssize(A.face(d)) != 3)
            return false;
    }
    return true;
}

bool criterionConcaveAnglesFar(const DerivedOnceNet& A) {
    for (const vector<Dart>& face: A.faces()) {
        for (const Dart& d: face) {
            if (!A.isConcaveAngle(d)) continue;
            
            int dist = 1;
            Dart curr = A.phi(d);
            while (curr != d) {
                if (A.isConcaveAngle(curr)) {
                    //~ assert(dist != 1);
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

bool criterionConvexFacesIntersection(const DerivedOnceNet& A) {
    for (Dart d: A.darts()) {
        if (!A.isConvexFace(d) or !A.isConvexFace(A.rev(d)))
            continue;
        vector<Dart> face2 = A.face(A.rev(d));
        int cnt = 0;
        for (Dart e: face2)
            cnt += A.sameFace(d, e);
        if (cnt > 1)
            return false;
    }
    return true;
}

bool criterionQuadrilateralsAreConvex(const DerivedOnceNet& A) {
    for (auto face: A.faces()) {
        if (A.faceType(face).first == 4 and !A.isConvexFace(face))
            return false;
    }
    
    return true;
}

bool isPesceMagro(const DerivedOnceNet& A) {
    auto faces = A.faces();
    if (ssize(faces)!= 3) return false;
    for (auto face: faces) {
        if (A.isConvexFace(face[0]) and ssize(face) != 4) return false;
        if (!A.isConvexFace(face[0]) and ssize(face) != 3) return false;
    }
    return true;
}

bool criterionTrianglesAreConvex(const DerivedOnceNet& A) {
    if (isPesceMagro(A)) return true;
    for (Dart d: A.darts()) {
        if (A.faceType(d) == pair<int,int>{3,1}) return false;
    }
    return true;        
}


bool criterionConsecutiveVerticesDegree2(const DerivedOnceNet& A) {
    int consecutive_deg2 = 0;
    for (Dart d: A.darts()) {
        int curr = 0;
        while (A.degree(A.phi(d, curr).from) == 2)
            curr++;
        consecutive_deg2 = max(consecutive_deg2, curr);
    }
    if (consecutive_deg2 == 3)
        return sumAreaConvexFaces(A) <= 7;
    if (consecutive_deg2 >= 4)
        return sumAreaConvexFaces(A) <= 5 and !containsBiangle(A);
    
    return true; 
}



int main() {
    vector<DerivedOnceNet> derived_once = 
        parseDerivedOnceNets("data/derived_once_computer.txt");
    cout << "There are " << derived_once.size()
         << " derived-once net candidates." << endl << endl;
         
    // The order of the criteria matters, as a criterion may implicitly
    // assume that all the previous ones pass.
    vector<pair<string, function<bool(const DerivedOnceNet&)>>> criteria = {  
        {"biangle\t\t\t", criterionBiangle},
        {"simple boundary\t\t", criterionSimpleBoundary},
        {"adjacent concave angles\t", criterionAdjacentConcaveAngles},
        {"concave angles are far\t", criterionConcaveAnglesFar},
        {"convex faces intersection", criterionConvexFacesIntersection},
        {"quadrilaterals are convex", criterionQuadrilateralsAreConvex},
        {"triangles are convex\t", criterionTrianglesAreConvex},
        {"consecutive degree 2 vertices", criterionConsecutiveVerticesDegree2},
    };
    
    map<string, vector<int>> excluded_by; 
    vector<int> survivors;

    for (int id = 0; id < ssize(derived_once); id++) {
        const DerivedOnceNet& A = derived_once[id];
        bool survives = true;
        for (auto pp: criteria) {
            string name = pp.first;
            auto criterion = pp.second;
            if (!criterion(A)) {
                survives = false;
                excluded_by[name].push_back(id);
                break;
            }
        }

        if (survives) 
            survivors.push_back(id);
    }
    
    cout << "Criterion\t\t\t\t# Excluded" << endl;
    for (auto pp: criteria) {
        string name = pp.first;
        cout << name << "\t\t" << excluded_by[name].size() << endl;
    }
    
    cout << "# survivors = " << survivors.size() << endl;
    
    ofstream out("data/derived_once_survivors.txt");
    out << survivors.size() << endl;
    for (int id: survivors)
        out << id << " ";
    out << endl;
    out.close();
    
    
    cout << "Nets which may have edges strictly longer than pi: " << endl;
    for (int id: survivors)
        if (mayHaveEdgesGtPi(derived_once[id]))
            cout << id << endl;
}
