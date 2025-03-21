#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "spherical_geometry.cpp"

void savePolygon(vector<GeodesicSegment> poly, string out_file) {
    ofstream out(out_file);
    out << ssize(poly) << endl;
    for (const GeodesicSegment& side: poly) {
        pt P = side.P;
        pt v = side.v;
        out << P.x << " " << P.y << " " << P.z << " | ";
        out << v.x << " " << v.y << " " << v.z << " | ";
        out << side.len << endl;
    }
    out.close();
}

bool completePolygon(vector<GeodesicSegment>& vec, 
                     double angle1, double angle2, double angle3) {
    auto missing_sides = completePolygon(vec.back(), vec[0], angle1, angle2, angle3);
    if (missing_sides[0].len < -0.5)
        return false;
    
    vec.pop_back();
    for (int i = 0; i < 3; i++) 
        vec.push_back(missing_sides[i]);
    
    for (int i = 0; i < ssize(vec); i++) {
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

vector<GeodesicSegment> constructPolygon(vector<double> angles, 
                                         vector<double> sides) {
    int n = angles.size();
    assert(ssize(sides) == n-3);
    vector<GeodesicSegment> vec = {GeodesicSegment()};
    for (int i = 0; i < n-3; i++) {
        vec.back().len = sides[i];
        TangentVector tv = vec.back();
        tv.exp(sides[i]);
        tv.rotate(angles[i]-M_PI);
        vec.push_back(GeodesicSegment(tv, 0));
    }
    if (completePolygon(vec, angles[n-3], angles[n-2], angles[n-1]))
        return vec;
    return {};
}


int main() {
    //~ vector<bool> big_angle = {0, 0, 0};
    //~ vector<double> len = {1, 2, 1};
    
    TangentVector tv1;
    double a120 = M_PI*0.666666;
    double a240 = M_PI*1.333333;
    vector<double> angles = {a120, a120, a120, a120};
    vector<double> sides = {1.5};
    vector<GeodesicSegment> poly = constructPolygon(angles, sides);
    assert(!poly.empty());
    savePolygon(poly, "todel.txt");
    
    //~ vector<bool> big_angle = {0, 0, 0, 1, 0, 0};
    //~ vector<double> len = {1.14, 1.7, 0.7, 0.2, 0.8, 0.93};
    //~ assert(ssize(big_angle) == ssize(len));
    //~ vector<TangentVector> poly = {TangentVector()};

    //~ for (int i = 0; i < ssize(len)-1; i++) {
        //~ TangentVector t = poly[i];
        //~ t.exp(len[i]);
        //~ if (big_angle[i])
            //~ t.rotate(M_PI * 0.33333333333);
        //~ else 
            //~ t.rotate(M_PI * -0.3333333333);
        //~ poly.push_back(t);
    //~ }
    
    //~ ofstream out("todel.txt");
    //~ out << ssize(poly) << endl;
    //~ for (int i = 0; i < ssize(poly); i++) {
        //~ pt P = poly[i].P;
        //~ pt v = poly[i].v;
        //~ out << P.x << " " << P.y << " " << P.z << " | ";
        //~ out << v.x << " " << v.y << " " << v.z << " | ";
        //~ out << len[i] << endl;
    //~ }
    //~ out.close();
}
