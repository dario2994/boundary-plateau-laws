#ifndef SPHERICAL_FACES
#define SPHERICAL_FACES

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "spherical_geometry.cpp"


// Check that all edges have a minimum length (to avoid precision issues)
// and that they are non-intersecting.
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

// Given n angles and and n-3 sides, it returns all (nonintersecting)
// polygons with this data.
// The tangent vector identifying the direction of the first edge from
// the first point can be given.
vector<vector<GeodesicSegment>> constructPolygon(
        vector<double> angles, vector<double> sides,
        vector<int> missing_sides,
        TangentVector initial_dir=TangentVector()) {
    int n = angles.size();
    assert(ssize(sides) == n);
    assert(n >= 3);
    assert(ssize(missing_sides) == 3);
    
    auto candidates = completeGeneralizedPolygon(angles, sides, missing_sides);
    
    vector<vector<GeodesicSegment>> res;
    for (auto ss: candidates) {
        vector<GeodesicSegment> vec = {{initial_dir, 0}};
        for (int i = 0; i < n; i++) {
            vec.back().len = ss[i];
            TangentVector tv = vec.back().exp(ss[i])
                                         .rot(angles[i]-M_PI);
            vec.push_back(GeodesicSegment(tv, 0));
        }
        assert(norm(vec[0].P-vec[n].P) < eps and norm(vec[0].V-vec[n].V) < eps);
        vec.pop_back();
        if (validatePolygon(vec))
            res.push_back(vec);
    }
    return res;
}

// Constructs all the faces with the given data that
// satisfies all the properties of a valid face of a net.
//
// In some rare cases there are two candidates. In the vast majority of
// cases, there is either one or two candidates.
// 
// angles[i] \in \{120, 240\}
// ssize(angles) == ssize(sides) + 3
vector<vector<GeodesicSegment>> constructFace(
        vector<int> angles, vector<double> sides, 
        vector<int> missing_sides,
        TangentVector initial_dir=TangentVector()) {
    int n = ssize(angles);
    assert(n >= 3);
    
    vector<double> angles_rad(n);
    for (int i = 0; i < n; i++) {
        assert(angles[i] == 120 or angles[i] == 240);
        angles_rad[i] = angles[i] * (M_PI/180);
    }
    auto candidates = 
        constructPolygon(angles_rad, sides, missing_sides, initial_dir);
    assert(candidates.size() <= 2);
    
    vector<vector<GeodesicSegment>> polys;
    
    vector<int> big_vertices;
    for (int i = 0; i < n; i++) {
        if (angles[i] == 240)
            big_vertices.push_back((i+1) % n);
    }
    
    assert(ssize(big_vertices) <= 2);

    for (auto poly: candidates) {
        bool works = true;
        if (big_vertices.size() == 1) {
            int u = big_vertices[0];
            TangentVector tu = poly[u].rot(-(2*M_PI)/3);
            poly.push_back({tu, 0.1});
        } else if (big_vertices.size() == 2) {
            int u = big_vertices[0];
            int v = big_vertices[1];
            TangentVector tu = poly[u].rot(-(2*M_PI)/3);
            TangentVector tv = poly[v].rot(-(2*M_PI)/3);
            double len_u = intersectWithEquator(tu, tv.equatorNormal());
            double len_v = intersectWithEquator(tv, tu.equatorNormal());
            
            GeodesicSegment gsu(tu, len_u);
            GeodesicSegment gsv(tv, len_v);
            
            // The angle in the pin belongs to [120,240].
            tu = tu.exp(len_u);
            tv = tv.exp(len_v);
            double theta = tu.inverseRot(-tv.V);
            if (min(theta, 2*M_PI-theta) > M_PI/3)
                works = false;
            
            poly.push_back(gsu);
            poly.push_back(gsv);
            
            if (norm(gsu.to() - gsv.to()) > eps)
                works = false;
            
            // The links to the pin stay inside the face.
            for (int i = 0; i < n; i++) {
                if ((i+1) % n != u and i != u and areIntersecting(poly[i], gsu))
                    works = false;
                if ((i+1) % n != v and i != v and areIntersecting(poly[i], gsv))
                    works = false;
            }
        }
        
        if (works)
            polys.push_back(poly);
    }
    
    assert(polys.size() <= 1);
    
    return polys;
}

#endif
