#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "derived_once_nets.cpp"
#include "metrization.cpp"


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

bool sameFace(const DerivedOnceNet& A, Dart d, Dart e) {
    const vector<Dart>& face = A.face(d);
    for (const Dart& e1: face) {
        if (e1 == e)
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

// Whether the associated net may contain an edge strictly longer than pi.
bool containsLongEdges(const DerivedOnceNet& A) {
    if (containsBiangle(A))
        return false;

    for (Dart d: A.darts()) {
        Dart e = A.rev(d);
        if (A.isConvexFace(d) or A.isConvexFace(e))
            continue;
        
        if (A.faceArea(d) + A.faceArea(e) <= 
                6 + pinContribution(A, d) + pinContribution(A, e))
            continue;

        //aggiunto da Camillo, vedi Lemma 4.7(iv)

        if ((A.faceArea(d)<=pinContribution(A,d)+4) 
                or (A.faceArea(e)<=pinContribution(A,e)+4))
            continue;

        return true;
        
    }
    return false;
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
            
        if (sameFace(A, A.rev(face[0]), A.rev(face[1])))
            return false;
            
        if (A.isConvexFace(A.rev(face[0])) 
            or A.isConvexFace(A.rev(face[1])))
            return false;
        
        // TODO: Do it as in containsLongEdges.
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
        if (sameFace(A, d, A.rev(d)))
            return false;
    }
    return true;
}

bool criterionAdjacentConcaveAngles(const DerivedOnceNet& A) {
    if (containsLongEdges(A))
        return true; // We cannot apply this criterion.
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
            cnt += sameFace(A, d, e);
        if (cnt > 1)
            return false;
    }
    return true;
}

bool criterionFourPentagons(const DerivedOnceNet& A) {
    for (const vector<Dart>& face: A.faces()) {
        if (A.faceType(face) != pair<int,int>{4, 0})
            continue;
        bool all_pentagons = true;
        for (Dart d: face) {
            auto pentagon = A.face(A.rev(d));
            all_pentagons &= A.faceType(A.rev(d)) == pair<int,int>{5, 0};
        }
        if (all_pentagons)
            return false;
    }
    return true;
}

bool criterionThreeQuadrilaterals(const DerivedOnceNet& A) {
    for (const Dart& d: A.darts()) {
        if (A.degree(d.from) != 3)
            continue;
        Dart a = d; 
        Dart b = A.next(a); 
        Dart c = A.next(b);
        if (A.faceType(a) == pair<int,int>{4, 0}
                and A.faceType(b) == pair<int,int>{4, 0}
                and A.faceType(c) == pair<int,int>{4, 0})
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


bool criterionTriangleQuadrilateral(const DerivedOnceNet& A) {
    for (Dart d: A.darts()) {
        Dart e = A.rev(d);
        if (A.faceType(d) != pair<int,int>{3, 0} 
            or A.faceType(e) != pair<int,int>{4, 0})
            continue;
        if (A.degree(A.phi(d).to) == 2)
            return false;
    }
    return true;
}

bool criterionPentagonThreeQuadrilaterals(const DerivedOnceNet& A) {
    for (Dart d: A.darts()) {
        if (A.faceType(d) != pair<int,int>{5, 0})
            continue;
        Dart e = A.phi(d);
        Dart f = A.phi(e);
        bool works = true;
        for (Dart z: {A.rev(d), A.rev(e), A.rev(f)}) {
            works &= A.faceType(z) == pair<int,int>{4, 0};
        }
        if (works)
            return false;
    }
    return true;
}

bool criterionTriangleQuadrilateral2(const DerivedOnceNet& A) {
    // TODO: Should be polished.
    for (Dart d: A.darts()) {
        Dart e = A.rev(d);
        if (A.faceType(d) != pair<int,int>{3, 0}
            or A.faceType(e) != pair<int,int>{4, 0})
            continue;
        Dart opposite = A.phi(A.phi(e));
        if (A.isConvexFace(A.rev(A.phi(e)))) return false;
        if (A.isConvexFace(A.rev(A.phi(opposite)))) return false;
        if (!A.isConvexFace(A.rev(opposite))) continue;
        if (ssize(A.face(A.rev(opposite))) != 4) return false;
        opposite = A.phi(A.phi(A.rev(opposite)));
        if (!A.isConvexFace(A.rev(opposite))) continue;
        if (ssize(A.face(A.rev(opposite))) != 4) return false;
        opposite = A.phi(A.phi(A.rev(opposite)));
        if (A.isConvexFace(A.rev(opposite))) return false;
    }
    return true;
}

bool criterionTwoQuadrilateralsThreePentagons(const DerivedOnceNet& A) {
    for (Dart d: A.darts()) {
        Dart e = A.rev(d);
        if (A.faceType(d) != pair<int,int>{4, 0}
            or A.faceType(e) != pair<int,int>{4, 0})
            continue;
            
        int cnt = 0;
        Dart z = A.phi(d);
        for (int i = 0; i < 3; i++) {
            Dart f = A.rev(z);
            if (A.faceType(f) == pair<int,int>{5, 0}) cnt++;
            z = A.phi(z);
        }
        
        z = A.phi(A.phi(e));
        if (A.faceType(A.rev(z))== pair<int,int>{5, 0}) cnt++;

        if (cnt >= 3)
            return false;
    }
    return true;
}

bool criterionConsecutiveVerticesDegree2(const DerivedOnceNet& A) {
    int consecutive_deg2 = 0;
    for (Dart d: A.darts()) {
        // TODO: This should be improved.
        //       Should we simplify the logic? Is it ok not to check
        //       that these are distinct vertices.
        // Camillo: I think it is OK, if they are not three distinct
        //        vertices then there would be a concave biangle,
        //        which is excluded. 
        if (A.degree(d.from) == 2 and
            A.degree(d.to) == 2 and
            A.degree(A.phi(d).to) == 2) {
            consecutive_deg2 = max(consecutive_deg2, 
                                   3 + (A.degree(A.phi(d, 2).to) == 2));
        }
    }
    if (consecutive_deg2 == 3)
        return sumAreaConvexFaces(A) <= 7;
    if (consecutive_deg2 == 4)
        return sumAreaConvexFaces(A) <= 5 and !containsBiangle(A);
    
    return true; 
}

bool criterionPent1(const DerivedOnceNet& A) {
    pair <int,int> pent1= {5,1};
    pair <int,int> quad = {5,0};
    for (Dart d:A.darts()) {
        if ((A.degree(d.from)!=2) or A.faceType(d)!=pent1) continue;
        Dart opp=A.phi(A.phi(d));
        if (A.faceType(A.rev(opp))!=quad) continue;
        if ((A.faceType(A.rev(A.phi(d)))==quad) or 
            (A.faceType(A.rev(A.phi(opp)))==quad)) return false;
    }
    return true;
}

//bool criterioSperimentale(const DerivedOnceNet& A) {
//    pair <int,int> quad = {4,0};
//    pair <int,int> pent = {5,0};
//    for (Dart d:A.darts()) {
//        if ((A.faceType(d)!=quad) or (A.faceType(A.rev(d))!=quad))
//            continue;
//        Dart e=A.rev(A.phi(d));
//        if (A.faceType(e)!=pent) continue;
//        if ((A.degree((A.phi(A.phi(d))).to)==2) and 
//            (A.degree((A.phi(A.phi(e))).to)==2)) return false;
//        if ((A.degree((A.phi(A.phi(A.rev(d)))).to)==2) and
//            (A.degree((A.phi(A.phi(A.phi(e)))).to)==2)) return false;        
//    }
//    return true;
//}

bool criterionMetrization(const DerivedOnceNet& A) {
    return MetricDerivedOnceNet(A).isConsistent();
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
        {"four pentagons\t\t", criterionFourPentagons},
        {"quadrilaterals are convex", criterionQuadrilateralsAreConvex},
        {"triangles are convex\t", criterionTrianglesAreConvex},
        {"three quadrilaterals\t", criterionThreeQuadrilaterals},
        {"triangle quadrilateral\t", criterionTriangleQuadrilateral},
        {"pentagon 3 quadrilaterals", criterionPentagonThreeQuadrilaterals},
        {"triangle quadrilateral 2", criterionTriangleQuadrilateral2},
        {"2 quadrilaterals 3 pentagons", criterionTwoQuadrilateralsThreePentagons},
        {"consecutive degree 2 vertices", criterionConsecutiveVerticesDegree2},
        {"metrization\t\t", criterionMetrization},
        //{"criterio sperimentale\t",criterioSperimentale},
        {"pentagon {5,1} + quadrilaterals", criterionPent1},
    };
    
    map<string, vector<int>> excluded_by; 
    vector<int> survivors;
    
    auto parse_subdirectories = [&](string directory_path) {
        vector<int> subdirs;

        for (const auto& entry: filesystem::directory_iterator(directory_path)) {
            if (entry.is_directory()) {
                string filename = entry.path().filename().string();
                int num = stoi(filename);
                subdirs.push_back(num);
            } 
        }
        
        return subdirs;
    }; 
    vector<int> blacklisted = parse_subdirectories("adhoc-blacklist/");
    vector<int> whitelisted = parse_subdirectories("adhoc-whitelist/");
    
    //~ vector<DerivedOnceNet> derived_once_old = 
        //~ parseDerivedOnceNets("data/derived_once_computer_old.txt");
        
    //~ auto findNet = [&](const DerivedOnceNet& A) {
        //~ auto B = symmetrize(A);
        //~ bool found = false;
        //~ for (int new_id = 0; new_id < ssize(derived_once); new_id++) {
            //~ if (areIsomorphic(A, derived_once[new_id])
                //~ or areIsomorphic(B, derived_once[new_id])) {
                //~ return new_id;
            //~ }
        //~ }
        //~ assert(0);
        //~ return -1;
    //~ };
        
    // Renaming whitelist and blacklist if the numbering changes.
    //~ cout << "BLACK: " << endl;
    //~ for (int old_id: blacklisted) {
        //~ int new_id = findNet(derived_once_old[old_id]);
        //~ cout << old_id << ": " << new_id << endl;
        //~ filesystem::rename("adhoc-blacklist/" + to_string(old_id),
                           //~ "adhoc-blacklist/" + to_string(new_id));
    //~ }
    
    //~ cout << "WHITE: " << endl;
    //~ for (int old_id: whitelisted) {
        //~ int new_id = findNet(derived_once_old[old_id]);
        //~ cout << old_id << ": " << new_id << endl;
        //~ filesystem::rename("adhoc-whitelist/" + to_string(old_id),
                           //~ "adhoc-whitelist/" + to_string(new_id));
    //~ }

    vector<int> excluded_white;
    vector<int> excluded_black;

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
        if (survives) {
            for (int id2: blacklisted) if (id2 == id) {
                survives = false;
                excluded_black.push_back(id);
            }
        }

        if (survives) {
            for (int id2: whitelisted) if (id2 == id) {
                survives = false;
                excluded_white.push_back(id);
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
    
    cout << "Blacklisted ad hoc\t\t\t" << excluded_black.size() << endl;
    cout << "Whitelisted\t\t\t\t" << excluded_white.size() << endl;
    cout << endl;    
    
    cout << "# survivors = " << survivors.size() << endl;
    
    ofstream out("data/derived_once_survivors.txt");
    out << survivors.size() << endl;
    for (int id: survivors)
        out << id << " ";
    out << endl;
    out.close();
    
    //~ for (int x: excluded_by["three quadrilaterals\t"])
        //~ cout << x << ", ";
}

// 788
