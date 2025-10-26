#ifndef SIMPLICIAL_COMPLEXES
#define SIMPLICIAL_COMPLEXES

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "utils.cpp"
#include "interval_geometry.cpp"


// Class representing an abstract simplicial complex with at least one vertex.
// In essence, this is a family of subsets of {0, 1, ..., V-1} which contains
// all singletons and is closed with respect to subset.
// On top of this, every vertex is mapped to a point in R^3.
struct ImmersedSimplicialComplex {
  private:
    int _V; // Number of vertices. The vertices are indexed 0, 1, ..., V-1.
    vector<array<double, 3>> points;
    vector<set<int>> faces; // All faces, not only maximal ones.

    void generateSubsetClosure() {
    set<set<int>> new_faces;
    for (set<int> face_set: faces) {
        int k = ssize(face_set);
        vector<int> face_vec(face_set.begin(), face_set.end());
        for (int bb = 0; bb < (1<<k); bb++) {
            set<int> subface;
            for (int i = 0; i < k; i++) {
                if (bb&(1<<i))
                    subface.insert(face_vec[i]);
            }
            new_faces.insert(subface);
        }
    }
    faces = vector<set<int>>(new_faces.begin(), new_faces.end());
}

  public:
    ImmersedSimplicialComplex(int V, vector<set<int>> faces): _V(V), points(V), faces(faces) {
        for (const auto& face: faces)
            for (int v: face)
                assert(0 <= v and v < V);
        generateSubsetClosure();
    }
    ImmersedSimplicialComplex(vector<array<double,3>> _points,
                              vector<set<int>> faces): ImmersedSimplicialComplex(ssize(_points), faces) {
        points = _points;
    }

    // friend ImmersedSimplicialComplex parseImmersedSimplicialComplex(istream&);
    friend ostream& operator<<(ostream&, const ImmersedSimplicialComplex&);

    int V() const { return _V; }
    double area() {
        double total_area = 0;
        for (auto face: faces) {
            if (ssize(face) != 3)
                continue;
            vector<array<double, 3>> vertices;
            for (int v: face)
                vertices.push_back(points[v]);
            array<double,3> triangle_vec = (vertices[2]-vertices[0])^(vertices[1]-vertices[0]);
            double triangle_area = 0.5 * norm(triangle_vec);
            total_area += triangle_area;
        }
        return total_area;
    }

    bool isCollapsible() {
        int F = ssize(faces);
        map<set<int>, int> face2id;
        for (const auto& face: faces)
            face2id[face] = ssize(face2id);

        vector<vector<int>> boundary(F);
        vector<vector<int>> coboundary(F);
        for (const auto& face: faces) {
            int id = face2id[face];
            for (int v: face) {
                set<int> subface = face;
                subface.erase(v);
                assert(face2id.count(subface));
                int sub_id = face2id[subface];
                boundary[id].push_back(sub_id);
                coboundary[sub_id].push_back(id);
            }
        }
        queue<int> free_faces;
        vector<bool> removed(F, false);
        for (int id = 0; id < F; id++) {
            if (ssize(coboundary[id]) == 1)
                free_faces.push(id);
        }

        auto remove_face = [&](int id) {
            assert(!removed[id]);
            assert(coboundary[id].empty());
            // cout << "Removing: ";
            // for (int v: faces[id])
            //     cout << v << " ";
            // cout << endl;
            removed[id] = true;
            for (int sub_id: boundary[id]) {
                vector<int>& vec = coboundary[sub_id];
                vec.erase(remove(vec.begin(), vec.end(), id), vec.end());
                if (ssize(vec) == 1)
                    free_faces.push(sub_id);
            }
        };

        while (!free_faces.empty()) {
            int id = free_faces.front();
            free_faces.pop();
            if (removed[id] or coboundary[id].empty())
                continue;
            assert(ssize(coboundary[id]) == 1);
            int super_id = coboundary[id][0];
            remove_face(super_id);
            remove_face(id);
        }

        for (int i = 0; i < F; i++)
            if (!removed[i])
                return false;
        return true;
    }

    // bool containsSegment(array<double, 3> A, array<double, 3> B) {
    //     return false;
        // TODO: This shall be implemented!
        //       Before implementing it, I should come up with a precise
        //       definition.
    // }
};


ImmersedSimplicialComplex parseImmersedSimplicialComplex(istream& in) {
    int V, F;
    in >> V >> F;
    vector<array<double, 3>> points(V);
    for (int i = 0; i < V; i++)
        in >> points[i][0] >> points[i][1] >> points[i][2];
    vector<set<int>> faces;
    for (int i = 0; i < F; i++) {
        int k;
        in >> k;
        set<int> face;
        for (int j = 0; j < k; j++) {
            int v;
            in >> v;
            face.insert(v);
        }
        faces.push_back(face);
    }

    return ImmersedSimplicialComplex(points, faces);
}

ostream& operator<<(ostream& out, const ImmersedSimplicialComplex& X) {
    out << "V = " << X.V() << ", faces = " << (ssize(X.faces)-X.V()-1) << endl;
    for (auto face: X.faces) {
        if (ssize(face) <= 1)
            continue;
        for (int v: face) cout << v << " ";
        out << endl;
    }
    return out;
}

void testIsCollapsible() {
    ImmersedSimplicialComplex boundary_triangle(3, {{0, 1}, {1, 2}, {2, 0}});
    ImmersedSimplicialComplex triangle(3, {{0, 1}, {1, 2}, {2, 0}, {0, 1, 2}});
    ImmersedSimplicialComplex boundary_tetrahedron(4, {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}});
    ImmersedSimplicialComplex tetrahedron(4, {{0, 1, 2, 3}});
    ImmersedSimplicialComplex annulus(6, {
        {0, 3, 1},
        {3, 1, 4},
        {1, 4, 2},
        {4, 2, 5},
        {2, 5, 0},
        {5, 0, 3}
    });
    ImmersedSimplicialComplex filled_annulus(6, {
        {0, 3, 1},
        {3, 1, 4},
        {1, 4, 2},
        {4, 2, 5},
        {2, 5, 0},
        {5, 0, 3},
        {0, 1, 2}
    });

    cout << "boundary_triangle: " << boundary_triangle.isCollapsible() << endl;
    cout << "triangle: " << triangle.isCollapsible() << endl;
    cout << "boundary_tetrahedron: " << boundary_tetrahedron.isCollapsible() << endl;
    cout << "tetrahedron: " << tetrahedron.isCollapsible() << endl;
    cout << "annulus: " << annulus.isCollapsible() << endl;
    cout << "filled_annulus: " << filled_annulus.isCollapsible() << endl;
}

int main() {
    // testIsCollapsible();
    ifstream foo_in("todel.txt");
    ImmersedSimplicialComplex X = parseImmersedSimplicialComplex(foo_in);
    cout << X.V() << endl;
    cout << X.isCollapsible() << endl;
}




#endif
