#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "spherical_geometry.cpp"

void savePolygon(vector<GeodesicSegment> poly, string out_file) {
    ofstream out(out_file);
    out << ssize(poly) << endl;
    for (const GeodesicSegment& side: poly) {
        side.serialize(out);
        out << endl;
    }
    out.close();
}

bool validatePolygon(const vector<GeodesicSegment>& vec) {
    for (int i = 0; i < ssize(vec); i++) {
        if (vec[i].len < 0.01 or vec[i].len > 2*M_PI-0.01)
            return false;
        int j = (i + 1) % ssize(vec);
        if (vec[i].len >= M_PI and vec[j].len >= M_PI)
            return false;
        for (int k = i+2; k < ssize(vec); k++) {
            if (i == 0 and k == ssize(vec)-1)
                continue;
            if (areIntersecting(vec[i], vec[k]))
                return false;
        }
    }
    return true;
}

vector<vector<GeodesicSegment>> completePolygon(vector<GeodesicSegment>& vec, 
                     double angle1, double angle2, double angle3) {
    auto completion_candidates = 
        completePolygon(vec.back(), vec[0].rotate(-(angle3-M_PI)), 
                        angle1, angle2);
        
    vector<vector<GeodesicSegment>> res;
    for (array<GeodesicSegment, 3> three_sides: completion_candidates) {
        auto vec2 = vec;
        vec2.pop_back();
        for (int i = 0; i < 3; i++) 
            vec2.push_back(three_sides[i]);
        if (validatePolygon(vec2))
            res.push_back(vec2);
    }
    
    return res;
}

// Given n angles and and n-3 sides, it returns all (nonintersecting)
// polygons with this data.
vector<vector<GeodesicSegment>> constructPolygon(vector<double> angles, 
                                                 vector<double> sides) {
    int n = angles.size();
    assert(ssize(sides) == n-3);
    vector<GeodesicSegment> vec = {GeodesicSegment()};
    for (int i = 0; i < n-3; i++) {
        vec.back().len = sides[i];
        TangentVector tv = vec.back().exp(sides[i])
                                     .rotate(angles[i]-M_PI);
        vec.push_back(GeodesicSegment(tv, 0));
    }
    
    return completePolygon(vec, angles[n-3], angles[n-2], angles[n-1]);
}

// angles[i] \in \{120, 240\}
vector<GeodesicSegment> constructFace(vector<int> angles,
                                      vector<double> sides) {
    static const double eps = 1e-6;
    
    int n = ssize(angles);
    assert(n >= 3);
    vector<double> angles_rad(n);
    for (int i = 0; i < n; i++) {
        assert(angles[i] == 120 or angles[i] == 240);
        angles_rad[i] = angles[i] * (M_PI/180);
    }
    auto polys = constructPolygon(angles_rad, sides);
    if (polys.empty())
        return {};
    assert(polys.size() == 1);
    auto poly = polys[0];
    
    vector<int> big_vertices;
    for (int i = 0; i < n; i++) {
        if (angles[i] == 240)
            big_vertices.push_back((i+1) % n);
    }
    
    assert(ssize(big_vertices) <= 2);
    
    if (big_vertices.size() == 1) {
        int u = big_vertices[0];
        TangentVector tu = poly[u].rotate(-(2*M_PI)/3);
        poly.push_back({tu, 0.1});
    } else if (big_vertices.size() == 2) {
        int u = big_vertices[0];
        int v = big_vertices[1];
        TangentVector tu = poly[u].rotate(-(2*M_PI)/3);
        TangentVector tv = poly[v].rotate(-(2*M_PI)/3);
        double len_u = intersectWithEquator(tu, tv.equatorNormal());
        double len_v = intersectWithEquator(tv, tu.equatorNormal());
        
        GeodesicSegment gsu(tu, len_u);
        GeodesicSegment gsv(tv, len_v);
        
        // The angle in the pin belongs to [120,240].
        tu = tu.exp(len_u);
        tv = tv.exp(len_v);
        double theta = tu.inverseRotate(-tv.V);
        //~ if (min(theta, 2*M_PI-theta) > M_PI/3)
            //~ return {};
        
        poly.push_back(gsu);
        poly.push_back(gsv);
        
        if (norm(gsu.to() - gsv.to()) > eps)
            return {};
        
        // The links to the pin stay inside the face.
        for (int i = 0; i < n; i++) {
            if ((i+1) % n != u and i != u and areIntersecting(poly[i], gsu))
                return {};
            if ((i+1) % n != v and i != v and areIntersecting(poly[i], gsv))
                return {};
        }
    }
    
    return poly;
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

void searchFaceRealization(vector<int> angles) {
    static random_device rd;  // Non-deterministic random number generator
    static mt19937 gen(rd()); // Mersenne Twister PRNG
    static uniform_real_distribution<double> rnd_side(0.02, 2*M_PI-0.02);

    
    int n = angles.size();
    assert(3 <= n and n <= 9);
    for (int _ = 0; _ < 10'000'000; _++) {
        vector<double> sides(n-3, -1);
        for (int i = 0; i < n-3; i++)
            sides[i] = rnd_side(gen);
        auto poly = constructFace(angles, sides);
        if (!poly.empty()) {
            for (int i = 0; i < n; i++)
                cout << poly[i].len << " ";
            cout << endl;
            savePolygon(poly, "todel.txt");
            return;
        }
    }
}

int main() {
    vector<int> angles = {120, 120, 120, 120, 240, 240};
    searchFaceRealization(angles);
}
