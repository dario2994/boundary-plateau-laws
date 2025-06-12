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

int main() {
    vector<DerivedOnceNet> derived_once =
        parseDerivedOnceNets("data/derived_once_computer.txt");
    cout << "There are " << derived_once.size()
         << " derived-once net candidates." << endl << endl;

    // The order of the criteria matters, as a criterion may implicitly
    // assume that all the previous ones pass.
    vector<pair<string, function<bool(const DerivedOnceNet&)>>> criteria = {
        // Only one net excluded by this criterion would not be excluded by
        // the following criteria. We keep it because it is useful to
        // know that all nets survived have this property.
        {"simple boundary\t\t", criterionSimpleBoundary},
        // The following two criteria are crucial. They exclude many nets and
        // some of the nets they exclude are either not succesfully excluded
        // by embed_or_exclude or it would require a lot of time to exclude
        // them with embed_or_exclude.
        {"biangle\t\t\t", criterionBiangle},
        {"concave angles are far\t", criterionConcaveAnglesFar}
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
}
