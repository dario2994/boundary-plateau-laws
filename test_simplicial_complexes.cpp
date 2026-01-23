#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "simplicial_complexes.cpp"

void testSegmentDistanceFromLine() {
    struct TestCase {
        string name;
        pt A, B, C, D;
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



void testTriangleDistanceFromLine() {
    struct TestCase {
        string name;
        pt A, B, C, D, E;
        double expected;
    };

    // Let the line AB be the z-axis for simplicity in most tests
    pt line_A = {0.0, 0.0, 1.0};
    pt line_B = {0.0, 0.0, -1.0};

    vector<TestCase> tests = {
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
    testIsCollapsible();
}
