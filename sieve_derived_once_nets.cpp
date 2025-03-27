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

int countObtuseAngles(const DerivedOnceNet& A, Dart d) {
    int res = 0;
    for (Dart e: A.face(d)) {
        res += A.isConcaveAngle(e);
    }
    return res;
}

int countPins(const DerivedOnceNet& A) {
    int cnt = 0;
    for (const vector<Dart>& face: A.faces()) {
        cnt += !A.isConvexFace(face[0]);
    }
    return cnt;
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

int pinContribution(const DerivedOnceNet& A, Dart d) {
    if (countObtuseAngles(A, d) == 1)
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

// Whether the associated net may contain an edge strictly longer than pi.
bool containsLongEdges(const DerivedOnceNet& A) {
    if (containsBiangle(A))
        return false;

    for (Dart d: A.darts()) {
        Dart e = A.rev(d);
        if (A.isConvexFace(d) or A.isConvexFace(e))
            continue;
        //~ auto face_d = A.face(d);
        //~ auto face_e = A.face(e);
        if (A.faceArea(d) + A.faceArea(e) <= 6 + pinContribution(A, d) + pinContribution(A, e))
            continue;
        return true;
        
    }
    return false;
}

int sumAreaConvexFaces(const DerivedOnceNet& A) {
    int res = 0;
    for (auto face: A.faces()) {
        if (A.isConvexFace(face[0]))
            res += A.faceArea(face[0]);
    }
    return res;
}


bool criterionBiangle(const DerivedOnceNet& A) {
    vector<vector<Dart>> faces = A.faces();
    
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
        
        int strange_area = sumAreaConvexFaces(A) - 4;
        strange_area += pinContribution(A, A.rev(face[0]));
        strange_area += pinContribution(A, A.rev(face[1]));
        
        if (strange_area > 5)
            return false;
            
        for (const Dart& d: A.darts()) {
            if (A.isConvexFace(d) and A.isConvexFace(A.rev(d))
                and ssize(A.face(d)) == 3
                and ssize(A.face(A.rev(d))) == 4)
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
        if (sameFace(A, d, A.rev(d)))
            return false;
    }
    return true;
}

bool criterionAdjacentConcaveAngles(const DerivedOnceNet& A) {
    if (containsLongEdges(A))
        return true; // We cannot apply this criterion.
    for (const Dart& d: A.darts()) {
        if (A.isConcaveAngle(d) and A.isConcaveAngle(A.phi(d)) and A.face(d).size() != 3)
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

bool criterionTwoAdjacentSquares(const DerivedOnceNet& A) {
    for (Dart d: A.darts()) {
        Dart e = A.rev(d);
        if (ssize(A.face(d)) != 4 or ssize(A.face(e)) != 4)
            continue;
        Dart d1 = A.rev(A.phi(A.phi(d)));
        Dart e1 = A.rev(A.phi(A.phi(A.phi(e))));
        
        assert(d1 != e1);
        if (A.isConcaveAngle(d1) and A.isConcaveAngle(e1)) 
            return false;
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
            cnt += sameFace(A, d, e);
        if (cnt > 1)
            return false;
    }
    return true;
}

bool criterionFourPentagons(const DerivedOnceNet& A) {
    for (const vector<Dart>& face: A.faces()) {
        if (ssize(face) != 4 or !A.isConvexFace(face[0]))
            continue;
        bool all_pentagons = true;
        for (Dart d: face) {
            auto pentagon = A.face(A.rev(d));
            all_pentagons &= ssize(pentagon) == 5 and A.isConvexFace(pentagon[0]);
        }
        if (all_pentagons)
            return false;
            
        //~ for (Dart d: face) {
            //~ Dart e = A.rev(A.phi(A.phi(A.rev(d))));
            //~ if (A.isConvexFace(e))
                //~ return false;
        //~ }
    }
    return true;
}

bool criterionThreeRectangles(const DerivedOnceNet& A) {
    for (const Dart& d: A.darts()) {
        if (A.degree(d.from) != 3)
            continue;
        Dart a = d; 
        Dart b = A.next(a); 
        Dart c = A.next(b);
        bool works = ssize(A.face(a)) == 4;
        works &= ssize(A.face(b)) == 4;
        works &= ssize(A.face(c)) == 4;
        
        works &= A.isConvexFace(a);
        works &= A.isConvexFace(b);
        
        works &= A.isConvexFace(c) or
            (A.isConcaveAngle(c) and A.isConcaveAngle(A.phi(A.phi(c))));
        
        if (works)
            return false;
    }
    return true;
}

bool criterionRectangleLongSide(const DerivedOnceNet& A) {
    for (auto face: A.faces()) {
        if (ssize(face) != 4 or countObtuseAngles(A, face[0]) != 1)
            continue;
        bool works = false;
        for (Dart d: face) {
            if (A.degree(d.from) < 3 or A.degree(d.to) < 3)
                continue;
            if (A.isConvexFace(A.rev(d)))
                continue;
            int weird_area = sumAreaConvexFaces(A);
            weird_area += pinContribution(A, d);
            weird_area += pinContribution(A, A.rev(d));
            if (weird_area <= 5)
                works = true;
        }
        if (!works)
            return false;
    }
    
    return true;
}

bool criterionTriangleRectangle(const DerivedOnceNet& A) {
    for (Dart d: A.darts()) {
        Dart e = A.rev(d);
        if (!A.isConvexFace(d) or !A.isConvexFace(e))
            continue;
        if (ssize(A.face(d)) != 3 or ssize(A.face(e)) != 4)
            continue;
        if (A.degree(A.phi(d).to) == 2)
            return false;
    }
    return true;
}

bool criterionQuadrilateralOnePin(const DerivedOnceNet& A) {
    for (auto face: A.faces()) {
        if (ssize(face) == 4 and countObtuseAngles(A, face[0]) == 1)
            return false;
    }
    return true;
}

bool criterionPentagonThreeQuadrilaterals(const DerivedOnceNet& A) {
    for (Dart d: A.darts()) {
        if (ssize(A.face(d)) != 5 or !A.isConvexFace(d))
            continue;
        Dart e = A.phi(d);
        Dart f = A.phi(e);
        bool works = true;
        for (Dart z: {A.rev(d), A.rev(e), A.rev(f)}) {
            works &= ssize(A.face(z)) == 4;
            works &= A.isConvexFace(z);
        }
        if (works)
            return false;
    }
    return true;
}

bool criterionTriangleRectangle2(const DerivedOnceNet& A) {
    for (Dart d: A.darts()) {
        Dart e = A.rev(d);
        if (ssize(A.face(d)) != 3 or !A.isConvexFace(d))
            continue;
        if (ssize(A.face(e)) != 4 or !A.isConvexFace(e))
            continue;
        
        for (Dart f: {A.phi(e), A.phi(A.phi(e)), A.phi(A.phi(A.phi(e)))}) {
            Dart g = A.rev(f);
            if (!A.isConvexFace(g))
                continue;
            if (ssize(A.face(g)) != 4)
                return false;
            if (f != A.phi(A.phi(e)))
                return false;
        }
    }
    return true;
}

bool criterionTwoQuadrilateralsFourPentagons(const DerivedOnceNet& A) {
    for (Dart d: A.darts()) {
        Dart e = A.rev(d);
        if (ssize(A.face(d)) != 4 or !A.isConvexFace(d))
            continue;
        if (ssize(A.face(e)) != 4 or !A.isConvexFace(e))
            continue;
            
        bool works = true;
        for (Dart _: {A.phi(d), A.phi(e)}) {
            Dart z = _;
            for (int i = 0; i < 3; i++) {
                Dart f = A.rev(z);
                works &= ssize(A.face(f)) == 5;
                works &= A.isConvexFace(f);
                z = A.phi(z);
            }
        }
        
        if (works)
            return false;
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
    vector<pair<string, function<bool(const DerivedOnceNet&)>>> criteria = {  
        {"biangle\t\t\t", criterionBiangle},
        {"simple boundary\t\t", criterionSimpleBoundary},
        {"adjacent concave angles\t", criterionAdjacentConcaveAngles},
        {"concave angles are far\t", criterionConcaveAnglesFar},
        {"two adjacent squares\t", criterionTwoAdjacentSquares},
        {"convex faces intersection", criterionConvexFacesIntersection},
        {"four pentagons\t\t", criterionFourPentagons},
        {"three rectangles\t", criterionThreeRectangles},
        {"rectangle long side\t", criterionRectangleLongSide},
        {"triangle rectangle\t", criterionTriangleRectangle},
        {"quadrilateral one pin\t", criterionQuadrilateralOnePin},
        {"pentagon 3 quadrilaterals", criterionPentagonThreeQuadrilaterals},
        {"triangle rectangle 2\t", criterionTriangleRectangle2},
        {"2 quadrilaterals 4 pentagons", criterionTwoQuadrilateralsFourPentagons}
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
    cout << endl;
    
    cout << "# survivors = " << survivors.size() << endl;
    
    ofstream out("data/derived_once_survivors.txt");
    out << survivors.size() << endl;
    for (int id: survivors)
        out << id << " ";
    out << endl;
    out.close();
}
