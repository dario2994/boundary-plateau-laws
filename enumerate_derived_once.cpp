#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "combinatorial_maps.cpp"
#include "derived_once_nets.cpp"


// Check that the number of edges of convex faces belongs to [2, 5].
bool checkConvexFaces(const DerivedOnceNet& A) {
    for (const vector<Dart>& face: A.faces()) {
        if (!A.isConvexFace(face))
            continue;
        if (ssize(face) < 2 or ssize(face) > 5)
            return false;
    }
    return true;
}

// Check that any face has at most 2 concave angles and its area is positive.
bool checkNonConvexFaces(const DerivedOnceNet& A) {
    for (const vector<Dart>& face: A.faces()) {
        if (A.isConvexFace(face))
            continue;

        int cnt_concave = A.faceType(face).second;
        if (cnt_concave > 2)
            return false;

        int area = A.faceArea(face);
        if (area <= 0 or area > 12)
            return false;
    }

    return true;
}

// Two adjacent convex faces have constraints on the number of sides.
bool checkAdjacentConvexFaces(const DerivedOnceNet& A) {
    for (const Dart& d: A.darts()) {
        const Dart& e = A.rev(d);
        if (!A.isConvexFace(d) or !A.isConvexFace(e))
            continue;
        int d_sides = A.face(d).size();
        int e_sides = A.face(e).size();
        if (d_sides > e_sides)
            swap(d_sides, e_sides);
        if (d_sides == 2 and e_sides != 2)
            return false;
        if (d_sides == 3 and e_sides == 5)
            return false;
    }
    return true;
}

// TODO: This code is repeated identical in metrized_nets.cpp.
//
// Consider the face containing d.
// Returns the boundary distance, avoiding d, between the two concave
// angles on the boundary of the face.
// Assumptions: The face containing d has two concave angles.
int concaveAnglesBoundaryDistanceAvoiding(const DerivedOnceNet& A, Dart d) {
    if (A.faceType(d).second != 2)
        return 0;

    int dist = -1;
    Dart curr = d;
    while (1) {
        if (A.isConcaveAngle(curr)) {
            if (dist == -1) dist = 1;
            else break;
        } else if (dist >= 0) dist++;

        curr = A.phi(curr);
    }
    return dist;
}

int pinContribution(const DerivedOnceNet& A, Dart d) {
    return max(4-concaveAnglesBoundaryDistanceAvoiding(A, d), 0);
}

int sumAreaConvexFaces(const DerivedOnceNet& A) {
    int res = 0;
    for (auto face: A.faces()) {
        if (A.isConvexFace(face))
            res += A.faceArea(face);
    }
    return res;
}

bool checkBiangle(const DerivedOnceNet& A) {
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

bool checkSimpleBoundary(const DerivedOnceNet& A) {
    vector<vector<Dart>> faces = A.faces();

    for (const Dart& d: A.darts()) {
        if (A.sameFace(d, A.rev(d)))
            return false;
    }
    return true;
}

bool checkConcaveAnglesFar(const DerivedOnceNet& A) {
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

bool validateDerivedOnceNet(const DerivedOnceNet& A) {
    // Exactly two non-convex faces.
    int cnt_non_convex_faces = 0;
    for (const auto& face: A.faces())
        cnt_non_convex_faces += !A.isConvexFace(face);
    assert(cnt_non_convex_faces == 2);

    if (!checkConvexFaces(A))
        return false;

    if (!checkNonConvexFaces(A))
        return false;

    if (!checkAdjacentConvexFaces(A))
        return false;

    if (!checkBiangle(A))
        return false;

    if (!checkSimpleBoundary(A))
        return false;

    if (!checkConcaveAnglesFar(A))
        return false;

    return true;
}

// Returns a copy of the original net, with some fake vertices added
// so that there are no more repeated edges.
// The list of integers returned represents the list of vertices that
// were added.
pair<DerivedOnceNet, vector<int>> removeRepeatedEdges(DerivedOnceNet A) {
    vector<int> new_vertices;
    while (1) {
        bool changed = false;
        set<pair<int,int>> edges;
        for (const Dart& d: A.darts()) {
            if (edges.count({d.from, d.to})) {
                A.edgeSubdivision(d);
                changed = true;
                new_vertices.push_back(A.V()-1);
                break;
            }
            edges.insert({d.from, d.to});
        }
        if (!changed) break;
    }
    return {A, new_vertices};
}

vector<DerivedOnceNet> addPin(const DerivedOnceNet& A) {
    vector<DerivedOnceNet> res;

    for (const vector<Dart>& face: A.faces()) {
        if (!A.isConvexFace(face))
            continue;
        // We add a pin inside face.
        // A pin may generate one or two concave angles.
        for (const Dart& d: face) {
            DerivedOnceNet A2 = A;
            pair<Dart, Dart> pp = A2.edgeSubdivision(d);
            res.push_back(A2); // One edge coming out of the pin.
            for (const Dart& e: A2.face(pp.first)) {
                DerivedOnceNet A3 = A2;
                A3.edgeSubdivision(e);
                res.push_back(A3); // Two edges coming out of the pin.
            }
        }
    }

    return res;
}

int main() {
    vector<CombinatorialMap> derived_twice_candidates =
        parseCombinatorialMaps("data/derived_twice_computer.txt");
    cout << "There are " << derived_twice_candidates.size()
         << " derived-twice spherical graph candidates." << endl << endl;

    CollectionUpToIsomorphism<DerivedOnceNet> derived_once_candidates;

    for (int i = 0; i < ssize(derived_twice_candidates); i++) {
        const CombinatorialMap& A = derived_twice_candidates[i];
        vector<DerivedOnceNet> pin1 = addPin(A);
        for (const DerivedOnceNet& A2: pin1) {
            vector<DerivedOnceNet> pin2 = addPin(A2);
            for (const DerivedOnceNet& A3: pin2)
                if (validateDerivedOnceNet(A3))
                    derived_once_candidates.add(A3);
        }
        if (i > 0 and i % 10 == 0) {
            cout << "Processed " << i
                 << " derived-twice net candidates, which generated "
                 << derived_once_candidates.getList().size()
                 << " derived-once net candidates." << endl;
        }
    }
    cout << endl;

    // Removing symmetric duplicates.
    vector<DerivedOnceNet> res = derived_once_candidates.getList();
    vector<DerivedOnceNet> res_up_to_symmetry;
    for (int i = 0; i < ssize(res); i++) {
        const DerivedOnceNet& sym = symmetrize(res[i]);
        bool is_new = true;
        for (int j = 0; j < i and is_new; j++) {
            is_new &= !areIsomorphic(sym, res[j]);
        }
        if (is_new)
            res_up_to_symmetry.push_back(sym);
    }
    res = res_up_to_symmetry;


    cout << "We found " << res.size()
         << " derived-once nets candidates." << endl;

    ofstream computer("data/derived_once_computer.txt");
    ofstream human("data/derived_once_human.txt");
    ofstream for_visualization("data/derived_once_for_visualization.txt");
    computer << res.size() << endl;
    for_visualization << res.size() << endl;
    for (int _ = 0; _ < ssize(res); _++) {
        const DerivedOnceNet& A = res[_];

        human << "Candidate " << _ << endl;
        human << A << endl;

        A.serialize(computer);

        // For visualization, we use the python library matplotlib.
        // This library does not support combinatorial maps with
        // repeated edges, so we create a version of the derived-once
        // net that does not contain any repeat edge.
        auto pp = removeRepeatedEdges(A);
        pp.first.serialize(for_visualization);
        for_visualization << pp.second.size() << "\n";
        for (int v: pp.second)
            for_visualization << v << " ";
        for_visualization << "\n";
    }
}