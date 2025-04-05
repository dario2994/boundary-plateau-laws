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
    //~ cout.precision(10);
    //~ for (int i = 0; i < ssize(sides); i++)
        //~ cout << sides[i] << " ";
    //~ cout << endl;
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
    //~ assert(ssize(valid_faces) <= 1);
    if (valid_faces.empty())
        return {};
        
    return valid_faces[0];
}

// Some values of sides may be negative. 
// Those are then randomly generated.
void searchFaceRealization(vector<double> angles, vector<double> sides) {
    static random_device rd;  // Non-deterministic random number generator
    static mt19937 gen(rd()); // Mersenne Twister PRNG
    static uniform_real_distribution<double> rnd_side(0.1, 2*M_PI-0.1);

    int n = angles.size();
    assert(3 <= n and n <= 9);
    for (int _ = 0; _ < 1'000'000; _++) {
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
            for (int i = 0; i < n; i++)
                cout << poly[i].len << " ";
            cout << endl;
            
            savePicture(poly, "todel.txt");
            return;
        }
    }
}

int main() {
    //~ TangentVector a(pt{0.9921950484, 0.1246955727, 0}, 
                    //~ pt{-0.06234778636, 0.4960975242, -0.8660254038});
    //~ TangentVector b(pt{-0.4809434217, 0.3948067189, -0.7828288956}, 
                    //~ pt{-0.8700617272, -0.3250090607, 0.3706233955});
    
    //~ pt equator_a = a.equatorNormal();
    //~ pt equator_b = b.equatorNormal();
    //~ cout << equator_a << " " << equator_b << endl;
    
    //~ double len = intersectWithEquator(a, equator_b);
    //~ cout << len << endl;
    
    double a120 = M_PI*0.66666666666666666666666;
    double a240 = 2*a120;
    
    //~ cout << "A" << endl;
    //~ searchFaceRealization({a120, a120, 2*a120, a120,a120, 2*a120}, 
                          //~ {-1, -1, -1, -1, -1, -1});
    
    //~ cout << "B" << endl;
    //~ searchFaceRealization({a120/2, a120, 2*a120}, 
                          //~ {-1, -1, -1});
    
    //~ cout << "C" << endl;
    //~ searchFaceRealization({a120, a120, a120, 2*a120}, 
                          //~ {-1, -1, -1, -1});                      
    
    double square = 2*atan(sqrt(0.5));
    cout << square << endl;
    searchFaceRealization({a120, a120, a120, a120, a240}, {-1, -1, M_PI+0.9, -1, -1});
    
    //~ auto angles = {a120, a120, a120, 2*a120, a120, 2*a120};
    //~ auto x = completeGeneralizedPolygon(angles, 
                          //~ {-1, -1, -1, 1.064049786, 0.1251442919, 3.080912009}, {0, 1, 2});
    //~ for (auto y: x) {
        //~ constructGeneralizedPolygon(angles, y);
    //~ }
}
