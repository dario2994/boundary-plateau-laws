#ifndef METRIZATION
#define METRIZATION

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "derived_once_nets.cpp"
#include "spherical_geometry.cpp"


int DEBUGGING;

// Returns a vector label, so that label[d.id] == label[e.id]
// implies that the lengths of darts d and e must be the same.
vector<int> identifyEqualSides(const DerivedOnceNet& A) {
    DEBUGGING = -A.E();
    vector<int> labels;
    
    // Init labels.
    for (Dart d: A.darts()) {
        assert(d.id == ssize(labels));
        labels.push_back(d.id);
    }
    
    // Cheap implementation of DSU.
    function<int(int)> ancestor = [&](int x) {
        if (labels[x] != x) labels[x] = ancestor(labels[x]);
        return labels[x];
    };
    
    function<void(Dart,Dart)> set_equal = [&](Dart d, Dart e) {
        int x = ancestor(d.id), y = ancestor(e.id);
        if (x != y) {
            labels[x] = y;
            DEBUGGING++;
        }
    };
    
    function<bool(Dart,Dart)> are_equal = [&](Dart d, Dart e) {
        int x = ancestor(d.id), y = ancestor(e.id);
        return x == y;
    };
    
    // Identification rules
    auto reverse_dart = [&](Dart d) {
        set_equal(d, A.rev(d));
    };
    
    auto convex_biangle = [&](Dart d) {
        if (A.faceType(d) != pair<int,int>{2, 0})
            return;
        set_equal(d, A.phi(d));
    };
    
    auto convex_triangle = [&](Dart d) {
        if (A.faceType(d) != pair<int,int>{3, 0})
            return;
        set_equal(d, A.phi(d));
        set_equal(d, A.phi(A.phi(d)));
    };
    
    auto convex_quadrilateral = [&](Dart d) {
        if (A.faceType(d) != pair<int,int>{4, 0})
            return;
        set_equal(d, A.phi(A.phi(d)));
        set_equal(A.phi(d), A.phi(A.phi(A.phi(d))));
    };
    
    auto pentagon = [&](Dart d) {
        if (ssize(A.face(d)) != 5)
            return;
        
        // Checking symmetry wrt d of concave angles.
        bool sym1 = A.isConcaveAngle(d) == A.isConcaveAngle(A.phi(d, 4));
        bool sym2 = A.isConcaveAngle(A.phi(d)) == A.isConcaveAngle(A.phi(d, 3));
        if (!sym1 or !sym2)
            return;
        
        // Do something.
        if (are_equal(A.phi(d), A.phi(d, 4)))
            set_equal(A.phi(d, 2), A.phi(d, 3));
        if (are_equal(A.phi(d, 2), A.phi(d, 3)))
            set_equal(A.phi(d), A.phi(d, 4));
    };
    
    while (1) {
        vector<int> bef = labels;
        for (Dart d: A.darts()) {
            reverse_dart(d);
            convex_biangle(d);
            convex_triangle(d);
            convex_quadrilateral(d);
            pentagon(d);
        }
        // Nothing changed during this iteration.
        if (bef == labels)
            break;
    }
    
    // Finalize labels.
    for (Dart d: A.darts())
        ancestor(d.id);
    
    return labels;
}

struct DartMetricData {
    bool known;
    double len;
    DartMetricData(): known(false) {}
    DartMetricData(double _len): known(true), len(_len) {}
    bool set(double _len, double eps = 1e-4) { 
        if (known) {
            return abs(len-_len) < eps;
        }
        known = true;
        len = _len; 
        return true;
    }
    
    double get() const {
        return known ? len : -1;
    }
};

vector<DartMetricData> determineSideLengths(const DerivedOnceNet& A) {
    static const vector<double> regular_convex_polygons = {-1, -1, 
        M_PI, // biangle
        acos(-1.0/3), // triangle
        2*atan(sqrt(0.5)), // quadrilateral
        acos(sqrt(5)/3)}; // pentagon
    
    vector<int> label = identifyEqualSides(A);
    vector<DartMetricData> lengths(2*A.E());
    
    auto get_length = [&](Dart d) {
        return lengths[label[d.id]].get();
    };
    
    auto set_length = [&](Dart d, double len) {
        return lengths[label[d.id]].set(len);
    };
    
    // TODO:
    // We should be able to mix the "equal sides" info with the "some
    // sides are determined info". Currently we don't do it.
    auto handle_regular_faces = [&](Dart d) {
        if (!A.isConvexFace(d))
            return true;
        int n = ssize(A.face(d));
        assert(2 <= n and n <= 5);
        double len = regular_convex_polygons[n];
        bool is_regular = true;
        for (int i = 0; i < n; i++) 
            is_regular &= label[d.id] == label[A.phi(d, i).id];
        if (!is_regular)
            return true;
        for (int i = 0; i < n; i++)
            if (!set_length(A.phi(d, i), len))
                return false;
        return true;
    };
    
    for (Dart d: A.darts())
        if (!handle_regular_faces(d))
            return {};
    
    auto complete_face = [&](Dart d) {
        vector<double> sides;
        vector<int> missing_sides;
        vector<double> angles;
        
        int n = A.faceType(d).first;
        for (int i = 0; i < n; i++) {
            Dart e = A.phi(d, i);
            sides.push_back(get_length(e));
            if (get_length(e) < 0)
                missing_sides.push_back(i);
            if (A.isConcaveAngle(e))
                angles.push_back(240 * (M_PI/180));
            else
                angles.push_back(120 * (M_PI / 180));
        } 
        
        assert(ssize(sides) == n);
        assert(ssize(angles) == n);
        
        if (n <= 3)
            return true;
        
        if (ssize(missing_sides) > 3) // Cannot be determined.
            return true;
        
        while (ssize(missing_sides) < 3) {
            int it = 0;
            while (it < ssize(missing_sides) and missing_sides[it] == it)
                it++;
            missing_sides.push_back(it);
            sort(missing_sides.begin(), missing_sides.end());
        }
        
        vector<vector<double>> candidates = 
            completeGeneralizedPolygon(angles, sides, missing_sides);
            
        if (ssize(candidates) == 3)
            return true; // Cannot be determined.
        
        vector<vector<double>> valid_faces;
        for (const auto& one_candidate: candidates) {
            auto poly = constructGeneralizedPolygon(angles, one_candidate);
            if (poly.empty())
                continue;
            assert(ssize(poly) == n);
            if (isValidPolygon(poly) and checkPinProperties(poly))
                valid_faces.push_back(one_candidate);
        }
        assert(ssize(valid_faces) <= 1);
        if (valid_faces.empty())
            return false;
        
        sides = valid_faces[0];
        for (int i = 0; i < n; i++)
            if (!set_length(A.phi(d, i), sides[i]))
                return false;
        return true;
    };
    
    // TODO: Repeating this 10 times is arbitrary, but should be ok.
    for (int i = 0; i < 10; i++) 
        for (Dart d: A.darts())
            if (!complete_face(d))
                return {};
    
    vector<DartMetricData> res(2*A.E());
    for (Dart d: A.darts())
        res[d.id] = lengths[label[d.id]];
        
    for (Dart d: A.darts()) {
        if (A.faceType(d) != pair<int,int>{5, 0})
            continue;
        if (res[d.id].known and res[d.id].len > regular_convex_polygons[3])
            return {};
    }
        
    return res;
}



vector<DerivedOnceNet> parseDerivedOnceNets(string in_file) {
    ifstream in(in_file);
    int cnt;
    in >> cnt;
    vector<DerivedOnceNet> res;
    for (int i = 0; i < cnt; i++) {
        const DerivedOnceNet& A = parseDerivedOnceNet(in);
        A.checkConsistency();
        res.push_back(A);
    }
    return res;
}

int main() {
    vector<DerivedOnceNet> derived_once = 
        parseDerivedOnceNets("data/derived_once_computer.txt");
    
    ifstream in("data/derived_once_survivors.txt");
    int cnt_survivors;
    in >> cnt_survivors;
    vector<int> survivors(cnt_survivors);
    for (int i = 0; i < cnt_survivors; i++)
        in >> survivors[i];
    in.close();
    
    int great = 0;
    map<int,vector<int>> measure_effectiveness;
    for (int id: survivors) {
        const DerivedOnceNet& A = derived_once[id];
        vector<DartMetricData> labels = determineSideLengths(A);
        if (labels.empty()) {
            great++;
            continue;
        }
        int cnt = 0;
        for (auto z: labels)
            cnt += z.known;
        measure_effectiveness[cnt].push_back(id);
        
        if (!labels.empty() and cnt == 0)
            cout << id << " ";
    }
    cout << endl;
    
    cout << great << endl;
    for (auto pp: measure_effectiveness)
        cout << pp.first << " " << pp.second.size() << endl;
    
    for (int x: measure_effectiveness[6])
        cout << x << " ";
    cout << endl;
}


#endif
