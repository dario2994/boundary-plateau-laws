#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "spherical_geometry.cpp"


vector<double> completeFace(vector<double> angles, vector<double> sides) {
    int n = ssize(angles);
    vector<int> missing;
    for (int i = 0; i < n; i++) {
        if (sides[i] < 0)
            missing.push_back(i);
    }
    
    vector<vector<double>> candidates = 
        completeGeneralizedPolygon(angles, sides, missing);
    assert(ssize(candidates) != 3);
        
    vector<vector<double>> valid_faces;
    for (const auto& one_candidate: candidates) {
        auto poly = constructGeneralizedPolygon(angles, one_candidate);
        if (poly.empty())
            continue;
        assert(ssize(poly) == ssize(angles));
        
        if (isValidPolygon(poly) and checkPinProperties(poly))
            valid_faces.push_back(one_candidate);
    }
    if (valid_faces.empty())
        return {};
    assert(ssize(valid_faces) == 1);
    return valid_faces[0];
}

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
        
        if (isValidPolygon(poly))
            valid_faces.push_back(poly);
        //~ if (isValidPolygon(poly) and checkPinProperties(poly))
            //~ valid_faces.push_back(poly);
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
    double a120 = M_PI*0.66666666666666666666666;
    double a240 = 2*a120;
    
    //~ searchFaceRealization({a120, a120, a120, a240, a120, a240}, {-1, -1, -1, -1, -1, -1});
    
    
    auto f = [&](double x, bool save)->pair<double,double> {
        vector<GeodesicSegment> picture;
        auto add_to_picture = [&](vector<GeodesicSegment> gs) {
            for (auto x: gs)
                picture.push_back(x);
        };
        
        auto new_face = [&](vector<double> angles, 
                            vector<double> sides, 
                            TangentVector tv_init)->vector<GeodesicSegment> {
            sides = 
                completeFace(angles, sides);
            if (sides.empty())
                return {};
            vector<GeodesicSegment> res = 
                constructGeneralizedPolygon(angles, sides, tv_init);
            add_to_picture(res);
            
            return res;
        };
        
        auto face1 = new_face({a120, a120, a120, a120}, 
                              {x, -1, -1, -1}, 
                              TangentVector());
        if (face1.empty())
            return {-1, -1};
        
        double y = face1[1].len;
        //~ cout << y << endl;
        
        auto face2 = new_face({a120, a120, a120, a120},
                              {x, -1, -1, -1},
                              face1[0].rev());
        if (face2.empty())
            return {-1, -1};
            
        auto face3 = new_face({a120, a120, a120, a120, a120},
                              {y, y, -1, -1, -1},
                              face1[1].rev());
        if (face3.empty())
            return {-1, -1};
        
        auto face4 = new_face({a120, a120, a120, a120, a120},
                              {y, y, -1, -1, -1},
                              face2[1].rev());
        if (face4.empty())
            return {-1, -1};                      
            
        auto face5 = new_face({a120, a120, a120, a120, a120},
                              {x, face3[4].len, -1, -1, -1},
                              face1[2].rev());
        if (face5.empty())
            return {-1, -1};
            
        if (save)
            savePicture(picture, "todel.txt");
        
        double len1 = face4[2].len;
        double len2 = face5[4].len;
        
        return {len1, len2};
    };
    
    // Configuration that we already understand.
    //~ for (double x = 1.45; x < 1.47; x += 0.0001) {
        //~ auto pp = f(x, false);
        //~ if (pp.first < 0)
            //~ continue;
        //~ double y = pp.first, z = pp.second;
        //~ if (abs(y-z) > 0.01)
            //~ continue;
        //~ cout << x << " " << y << " " << z << " " << y-z << endl;
    //~ }
    
    //~ f(1.4626, true);
    
    // Another configuration.
    
    for (double x = 0.1; x < 2*M_PI; x += 0.001) {
        auto pp = f(x, false);
        if (pp.first < 0)
            continue;
        double y = pp.first, z = pp.second;
        //~ if (abs(y-z) > 0.01)
            //~ continue;
        cout << x << " " << y << " " << z << " " << y-z << endl;
    }
    
}
