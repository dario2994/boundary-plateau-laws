#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "spherical_geometry.cpp"
#include "spherical_faces.cpp"



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
        for (auto tri: sides) {
            TangentVector tv;
            tv = tv.exp(tri[2]).rotate(B-M_PI);
            tv = tv.exp(tri[0]).rotate(C-M_PI);
            tv = tv.exp(tri[1]).rotate(A-M_PI);
            assert(norm(tv.P-TangentVector().P) < eps);
            assert(norm(tv.V-TangentVector().V) < eps);
        }
    }
}

// Some values of sides may be negative. 
// Those are then randomly generated.
void searchFaceRealization(vector<int> angles, vector<double> sides) {
    static random_device rd;  // Non-deterministic random number generator
    static mt19937 gen(rd()); // Mersenne Twister PRNG
    static uniform_real_distribution<double> rnd_side(0.02, 2*M_PI-0.02);

    int n = angles.size();
    assert(3 <= n and n <= 9);
    for (int _ = 0; _ < 10'000'000; _++) {
        vector<double> sides2 = sides;
        for (int i = 0; i < n-3; i++) {
            if (sides2[i] < 0)
                sides2[i] = rnd_side(gen);
        }
        auto polys = constructFace(angles, sides2);
        if (!polys.empty()) {
            auto poly = polys[0];
            for (int i = 0; i < n; i++)
                cout << poly[i].len << " ";
            cout << endl;
            
            savePicture(poly, "todel.txt");
            return;
        }
    }
}

void searchFaceBounds(vector<int> angles) {
    static random_device rd;  // Non-deterministic random number generator
    static mt19937 gen(rd()); // Mersenne Twister PRNG
    static uniform_real_distribution<double> rnd_side(0.02, 2*M_PI-0.02);
    
    int n = angles.size();
    assert(3 <= n and n <= 9);
    
    vector<double> _min(n, 1000);
    vector<double> _max(n, 0);
    
    int cnt = 0;
    for (int _ = 0; _ < 1e9; _++) {
        if (_ > 2e6 and cnt == 0)
            break;
        if (_ >= 1e7 and cnt >= 1e4)
            break;
        vector<double> sides(n-3, -1);
        for (int i = 0; i < n-3; i++)
            sides[i] = rnd_side(gen);
        auto polys = constructFace(angles, sides);
        //~ if (polys.size() >= 2) {
            //~ vector<GeodesicSegment> vec = polys[0];
            //~ for (auto gs: polys[1])
                //~ vec.push_back(gs);
            //~ savePicture(vec, "todel.txt");
            //~ exit(0);
        //~ }
        for (auto poly: polys) {
            cnt += 1;
            for (int i = 0; i < n; i++) {
                double l = poly[i].len;
                _min[i] = min(_min[i], l);
                _max[i] = max(_max[i], l);
            }
        }
    }
    
    cout << "Found " << cnt << " candidates with angles: ";
    for (int z: angles) cout << z << " ";
    cout << endl;
    if (cnt != 0) {
        for (int i = 0; i < n; i++)
            cout << "[" << _min[i] << ", " << _max[i] << "]" << endl;
    }
    cout << endl;
}

int main() {
    for (vector<int> angles: vector<vector<int>>{
            {120,120,120}, 
            {120,120,120,120}, 
            {120,120,120,120,120},
            {120,120,240},
            {120,240,240},
            {120,120,120,240},
            {120,120,240,240},
            {120,240,120,240}
    })
        searchFaceBounds(angles);
    
    for (vector<int> angles: vector<vector<int>>{
            {120,120,120,120,120},
            {120,120,120,120,240},
            {120,120,120,240,240},
            {120,120,240,120,240},
            {120,120,120,120,120,240},
            {120,120,120,120,240,240},
            {120,120,120,240,120,240},
            {120,120,240,120,120,240}
    })
        searchFaceBounds(angles);
        
    for (vector<int> angles: vector<vector<int>>{
            {120,120,120,120,120,120,240},
            {120,120,120,120,120,240,240},
            {120,120,120,120,240,120,240},
            {120,120,120,240,120,120,240},
            {120,120,120,120,120,240,120,240},
            {120,120,120,120,240,120,120,240},
            {120,120,120,240,120,120,120,240},
            {120,120,120,120,120,240,120,120,240},
            {120,120,120,120,240,120,120,120,240}
    })
        searchFaceBounds(angles);
}
