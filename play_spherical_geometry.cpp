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
            tv = tv.exp(tri[2]).rot(B-M_PI);
            tv = tv.exp(tri[0]).rot(C-M_PI);
            tv = tv.exp(tri[1]).rot(A-M_PI);
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
        vector<int> missing_sides;
        for (int i = 0; i < n; i++) {
            if (sides2[i] < 0) {
                if (ssize(missing_sides) < 3)
                    missing_sides.push_back(i);
                sides2[i] = rnd_side(gen);
            }
        }
        assert(ssize(missing_sides) == 3);
        
        auto polys = constructFace(angles, sides2, missing_sides);
        if (!polys.empty()) {
            auto poly = polys[0];
            for (int i = 0; i < n; i++)
                cout << poly[i].len << " ";
            cout << endl;
            
            cout << poly[1].len - poly[2].len << endl;
            
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
        vector<double> sides(n, -1);
        for (int i = 0; i < n; i++)
            sides[i] = rnd_side(gen);
        auto polys = constructFace(angles, sides, {0, 1, 2});
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
    searchFaceRealization({120, 120, 120}, {-1, -1, -1});
    searchFaceRealization({120, 120, 120, 120}, {2*atan(sqrt(0.5)), -1, -1, -1});
    searchFaceRealization({120, 120, 120, 120, 240}, {0.9153972, 0.9153972, -1, -1, -1});
    searchFaceRealization({120, 120, 120, 120, 120}, {0.72972766, 0.72972766, -1, -1, -1});
    
    cout.precision(10);
    for (double l = 0.72972766; l < 0.7297277; l += 1e-8) {
        searchFaceRealization({120, 120, 120, 120, 120}, {0.72972766, 0.72972766, -1, -1, -1});
    }
    
    //~ searchFaceRealization({120, 120, 240, 120, 120, 240}, {-1, 1, -1, -1, -1, -1});
    //~ double a120 = 120*M_PI/180;
    
    //~ SimpleMap sm = SimpleMap().rot(a120-M_PI).exp(1).rot(a120-M_PI);
    //~ cout << sm.alpha - M_PI << endl;
    //~ return 0;
    //~ auto res = completeGeneralizedPolygon({a120, a120, a120, a120}, {1, 1, 1, 1}, {1, 2, 3});
    //~ constructPolygon({120, 120, 120, 120}, {1, 1, 1, 1}, {1, 2, 3});
    //~ for (int i = 0; i < 4; i++) cout << res[0][i] << endl;
    
    //~ for (vector<int> angles: {120, 120, 120, 120, 240})
        //~ searchFaceBounds({120, 120, 120, 120, 240});
    
    //~ function<vector<GeodesicSegment>(double)> f = [&](double len) {
        //~ double a = M_PI/3;
        //~ auto polys = 
            //~ constructPolygon({2*a, a, 2*a, 4*a}, {len, 0, 0, 0}, {1, 2, 3});
        //~ if (polys.empty())
            //~ return vector<GeodesicSegment>();
        //~ for (auto poly: polys) {
            //~ if (validatePolygon(poly))
                //~ return poly;
        //~ }
        //~ return vector<GeodesicSegment>();
    //~ };
    
    //~ for (double len = 0.1; len < M_PI; len += 0.1)
        //~ cout << len << ": " << f(len).size() << endl;
}
