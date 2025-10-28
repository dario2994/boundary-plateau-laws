#ifndef SIMPLICIAL_COMPLEXES
#define SIMPLICIAL_COMPLEXES

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "utils.cpp"
#include "interval_geometry.cpp"

// Distance of the point P to the line through A and B.
// Assumption: |A-B| is not much smaller than |A|, |B|.
double distanceFromLine(const array<double, 3>& A,
                        const array<double, 3>& B,
                        const array<double, 3>& P) {
    assert(norm(A-B) > 0.01 * max(norm(A), norm(B)));
    double parallelogram_area = norm((B-A) ^ (P-A));
    double line_segment_length = norm(B-A);

    return parallelogram_area / line_segment_length;
}

// Distance of the segment CD to the line through A and B.
// Assumption: |A-B| is not much smaller than |A|, |B|.
double distanceFromLine(const array<double, 3>& A,
                        const array<double, 3>& B,
                        const array<double, 3>& C,
                        const array<double, 3>& D) {
    assert(norm(A-B) > 0.01 * max(norm(A), norm(B)));

    array<double, 3> x = (B - A) ^ (D - C);
    array<double, 3> y = (B - A) ^ (A - C);

    double s_unclamped = (norm(x) < 1e-5) ? 0.0 : (x * y) / (x * x);
    double s = max(0.0, min(1.0, s_unclamped));

    array<double, 3> closest_point = C + (D - C) * s;

    return distanceFromLine(A, B, closest_point);
}

void testSegmentDistanceFromLine() {
    struct TestCase {
        string name;
        array<double, 3> A, B, C, D;
        double expected;
    };

    vector<TestCase> tests = {
        {"Intersection",
         {0.0, 0.0, 0.5}, {0.0, 0.0, -0.5},
         {0.5, 0.0, 0.0}, {-0.5, 0.0, 0.0}, 0.0},
        {"Closest interior",
         {0.0, 0.0, 0.5}, {0.0, 0.0, -0.5},
         {0.5, 0.5, 0.0}, {-0.5, 0.5, 0.0}, 0.5},
        {"Closest is C (s=0)",
         {0.0, 0.0, 0.5}, {0.0, 0.0, -0.5},
         {0.2, 0.5, 0.0}, {0.8, 0.5, 0.0}, sqrt(0.29)},
        {"Closest is D (s=1)",
         {0.0, 0.0, 0.5}, {0.0, 0.0, -0.5},
         {-0.8, 0.5, 0.0}, {-0.2, 0.5, 0.0}, sqrt(0.29)},
        {"Parallel lines",
         {0.5, 0.0, 0.0}, {-0.5, 0.0, 0.0},
         {0.5, 0.8, 0.0}, {-0.5, 0.8, 0.0}, 0.8}
    };


    for (const auto& test : tests) {
        double actual = distanceFromLine(test.A, test.B, test.C, test.D);
        assert(abs(actual - test.expected) < 1e-5);
    }
}

// Distance of the triangle CDE to the line through A and B.
// Assumption: |A-B| is not much smaller than |A|, |B|.
double distanceFromLine(const array<double, 3>& A,
                        const array<double, 3>& B,
                        const array<double, 3>& C,
                        const array<double, 3>& D,
                        const array<double, 3>& E) {
    assert(norm(A-B) > 0.01 * max(norm(A), norm(B)));
    double boundary_distance = min({
        distanceFromLine(A, B, C, D),
        distanceFromLine(A, B, D, E),
        distanceFromLine(A, B, E, C),
    });

    // The answer is either boundary_distance or it is 0.
    // Determining whether is zero is annoying.
    // We use the Moller-Trumbore algorithm.
    // The implementation is freely inspired (i.e., copied) from wikipedia.
    array<double, 3> CD = D - C;
    array<double, 3> CE = E - C;
    array<double, 3> CA = A - C;

    array<double, 3> U = (B-A) ^ CE;
    array<double, 3> V = CA ^ CD;

    double det = CD * U;

    // There is a plane that "almost" contains all points.
    if (abs(det) < 1e-5)
        return boundary_distance;

    double inv_det = 1.0 / det;
    double u = inv_det * (CA * U);
    double v = inv_det * ((B-A) * V);

    // The intersection between the line and the plane containing the triangle
    // falls outside of the triangle.
    if (u <= 0 or v <= 0 or 1 <= u + v)
        return boundary_distance;

    return 0.0;
}


#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <cassert>
#include <cmath> // For sqrt
#include <algorithm> // For min

// --- Assume your corrected function and its dependencies are here ---
// (e.g., operator-, operator^, norm, dot-product,
// and the 4-argument distanceFromLine for segments)
//
// double distanceFromLine(const array<double, 3>& A, ...);
// double distanceFromLine(const array<double, 3>& A, ...);
//
// --- End of assumed code ---


void testTriangleDistanceFromLine() {
    struct TestCase {
        std::string name;
        std::array<double, 3> A, B, C, D, E;
        double expected;
    };

    // Let the line AB be the z-axis for simplicity in most tests
    std::array<double, 3> line_A = {0.0, 0.0, 1.0};
    std::array<double, 3> line_B = {0.0, 0.0, -1.0};

    std::vector<TestCase> tests = {
        {"Intersection",
         line_A, line_B, // z-axis
         {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {-1.0, -1.0, 0.0}, // Triangle in z=0 plane, contains origin
         0.0},

        {"Closest to Edge Interior",
         line_A, line_B, // z-axis
         {1.0, 1.0, 0.0}, {1.0, -1.0, 0.0}, {1.0, 0.0, 1.0}, // Triangle in x=1 plane
         1.0}, // Closest point is (1,0,0) on the interior of edge C-D

        {"Closest to Vertex C",
         line_A, line_B, // z-axis
         {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {1.0, 2.0, 0.0}, // Triangle in z=0, all in x>0, y>0
         sqrt(2.0)}, // Closest point is (1,1,0), which is vertex C

        {"Parallel Line",
         line_A, line_B, // z-axis
         {2.0, 0.0, 0.0}, {2.0, 1.0, 0.0}, {2.0, 0.0, 1.0}, // Triangle in x=2 plane
         2.0}, // Closest edge is C-E, distance is 2.0

        {"Degenerate Triangle (Miss)",
         line_A, line_B, // z-axis
         {1.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {1.0, 0.0, 2.0}, // Collinear points in x=1
         1.0}, // Should fall back to boundary, min dist is 1.0

        {"Degenerate Triangle (Hit)",
         line_A, line_B, // z-axis
         {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, // Collinear points on x-axis
         0.0} // Should fall back to boundary, segment D-E contains origin
    };


    for (const auto& test : tests) {
        double actual = distanceFromLine(test.A, test.B, test.C, test.D, test.E);
        assert(abs(actual - test.expected) < 1e-5);
    }
}


// Class representing an abstract simplicial complex with at least one vertex.
// In essence, this is a family of subsets of {0, 1, ..., V-1} which contains
// all singletons and is closed with respect to subset.
// On top of this, every vertex is mapped to a point in R^3.
struct ImmersedSimplicialComplex {
  private:
    // The vertices are indexed 0, 1, ..., V-1.
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
    int V() const { return ssize(points); }

    ImmersedSimplicialComplex(int _V, vector<set<int>> faces): points(_V), faces(faces) {
        for (const auto& face: faces)
            for (int v: face)
                assert(0 <= v and v < V());
        generateSubsetClosure();
    }
    ImmersedSimplicialComplex(vector<array<double,3>> _points,
                              vector<set<int>> faces): ImmersedSimplicialComplex(ssize(_points), faces) {
        points = _points;
    }

    friend ostream& operator<<(ostream&, const ImmersedSimplicialComplex&);

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

    bool containsSegment(const array<double, 3>& A,
                         const array<double, 3>& B) {
        const double eps = 1e-5;

        // The vertices are mapped on the line AB.
        vector<bool> on_the_line(V(), false);
        vector<int> points_on_line;
        for (int i = 0; i < V(); i++) {
            double dist = distanceFromLine(A, B, points[i]);
            if (dist > 100*eps)
                continue;
            // We expect a quantitative gap between points on the line and points
            // out of the line.
            if (dist > eps)
                return false;
            points_on_line.push_back(i);
            on_the_line[i] = true;
        }

        // We save the edges of the simplicial complex in a convenient
        // datastructure to speed up the computation.
        set<pair<int,int>> edges;
        for (auto face: faces) {
            if (ssize(face) != 2)
                continue;
            int x = *face.begin();
            int y = *face.rbegin();
            edges.insert({x, y});
            edges.insert({y, x});
        }

        // Sorting the vertices by the order in which they are mapped on AB,
        // we find a path (using edges of the simplicial complex).
        auto position_on_line = [&](int i) {
            return (B-A) * points[i];
        };

        sort(points_on_line.begin(), points_on_line.end(), [&](int i, int j) {
            return position_on_line(i) < position_on_line(j);
        });

        if (ssize(points_on_line) < 2)
            return false;
        for (int i = 0; i < ssize(points_on_line)-1; i++)
            if (!edges.count({points_on_line[i], points_on_line[i+1]}))
                return false;

        // The first vertex is mapped onto A.
        if (norm(points[points_on_line[0]]-A) > eps)
            return false;

        // The last vertex is mapped onto B.
        if (norm(points[points_on_line.back()]-B) > eps)
            return false;
        cout << "#Points on segments = " << ssize(points_on_line) << endl;

        // Any face that has only vertices in the sequence is either a
        // singleton or an edge of the path.
        for (int i = 0; i < ssize(points_on_line); i++) {
            for (int j = i+2; j < ssize(points_on_line); j++) {
                int u = points_on_line[i];
                int v = points_on_line[j];
                if (edges.count({u, v}))
                    return false;
            }
        }

        // All faces (edges and triangles) that don't have any vertex in the
        // sequence do not intersect the line.
        for (auto face: faces) {
            if (ssize(face) <= 1)
                continue;
            bool has_vertex_on_line = false;
            for (int v: face) {
                if (on_the_line[v])
                    has_vertex_on_line = true;
            }
            if (has_vertex_on_line)
                continue;

            assert(ssize(face) <= 3);
            if (ssize(face) == 2) {
                array<double, 3> C = points[*face.begin()];
                array<double, 3> D = points[*face.rbegin()];
                if (distanceFromLine(A, B, C, D) < 100*eps)
                    return false;
            }
            if (ssize(face) == 3) {
                vector<int> ids(face.begin(), face.end());
                if (distanceFromLine(A, B,
                                     points[ids[0]], points[ids[1]], points[ids[2]]) < 100*eps)
                    return false;
            }

        }

        return true;
    }
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
    testSegmentDistanceFromLine();
    testTriangleDistanceFromLine();
    // return 0;

    // testIsCollapsible();
    ifstream foo_in("todel.txt");
    ImmersedSimplicialComplex X = parseImmersedSimplicialComplex(foo_in);
    cout << X.V() << endl;
    cout << X.isCollapsible() << endl;
    cout << X.area() << endl;

    array<double, 3> pin1 = {0.647757, -0.761847, 0};
    array<double, 3> pin2 = {-pin1[0], -pin1[1], -pin1[2]};
    cout << X.containsSegment(pin1, pin2) << endl;
}




#endif
