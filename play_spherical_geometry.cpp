#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "spherical_geometry.cpp"
//~ #include "spherical_faces.cpp"



// Saves a collection of geodesic segments to file (in a format suitable
// for visualization).
void savePicture(vector<GeodesicSegment> vec, string out_file) {
    ofstream out(out_file);
    out << ssize(vec) << endl;
    for (const GeodesicSegment& side: vec) {
        side.serialize(out);
        out << endl;
    }
    out.close();
}

void unitTest() {
    for (double A = 0.1; A < 2*M_PI; A += 0.1) 
    for (double B = 0.1; B < 2*M_PI; B += 0.1)
    for (double C = 0.1; C < 2*M_PI; C += 0.1) {
        auto sides = computeTriangleSides(A, B, C);
        if (ssize(sides) == 3)
            continue;
        for (auto tri: sides) {
            TangentVector tv;
            tv = tv.exp(tri[2]).rot(B-M_PI);
            tv = tv.exp(tri[0]).rot(C-M_PI);
            tv = tv.exp(tri[1]).rot(A-M_PI);
            assert(norm(tv.P-TangentVector().P) < eps);
            assert(norm(tv.V-TangentVector().V) < eps);
        }
    }
}

vector<GeodesicSegment> findFace(vector<double> angles, vector<double> sides, 
                        vector<int> missing_sides) {
    vector<vector<double>> candidates = 
        completeGeneralizedPolygon(angles, sides, missing_sides);
    if (ssize(candidates) == 3)
        return {}; // Cannot be determined.
        
    vector<vector<GeodesicSegment>> valid_faces;
    for (const auto& one_candidate: candidates) {
        auto poly = constructGeneralizedPolygon(angles, one_candidate);
        if (poly.empty())
            continue;
        assert(ssize(poly) == ssize(angles));
        
        if (isValidPolygon(poly) and checkPinProperties(poly))
            valid_faces.push_back(poly);
    }
    assert(ssize(valid_faces) <= 1);
    if (valid_faces.empty())
        return {};
        
    return valid_faces[0];
}

// Some values of sides may be negative. 
// Those are then randomly generated.
void searchFaceRealization(vector<double> angles, vector<double> sides) {
    static random_device rd;  // Non-deterministic random number generator
    static mt19937 gen(rd()); // Mersenne Twister PRNG
    static uniform_real_distribution<double> rnd_side(0.02, 2*M_PI-0.02);

    int n = angles.size();
    assert(3 <= n and n <= 9);
    for (int _ = 0; _ < 100'000; _++) {
        vector<double> sides2 = sides;
        vector<int> missing_sides;
        for (int i = 0; i < n; i++) {
            if (sides2[i] < 0) {
                if (ssize(missing_sides) < 3)
                    missing_sides.push_back(i);
                sides2[i] = rnd_side(gen);
            }
        }
        assert(ssize(missing_sides) == 3);
        
        auto poly = findFace(angles, sides2, missing_sides);
        if (!poly.empty()) {
            TangentVector a = poly[0].rot(-M_PI/3);
            TangentVector b = poly[2].rot(-M_PI/3);
            poly.push_back({a, 3});
            poly.push_back({b, 3});
            
            for (int i = 0; i < n; i++)
                cout << poly[i].len << " ";
            cout << endl;
            
            //~ cout << poly[1].len - poly[2].len << endl;
            
            savePicture(poly, "todel.txt");
            return;
        }
    }
}

int main() {
    double a120 = M_PI*0.66666666666666666666666;
    searchFaceRealization({a120, a120, a120, a120, a120}, {0.01, 1.85, -1, -1, -1});
}
