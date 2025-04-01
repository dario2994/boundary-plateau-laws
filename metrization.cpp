#ifndef METRIZATION
#define METRIZATION

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "derived_once_nets.cpp"
#include "spherical_geometry.cpp"

// Information on the length of a dart.
// If consistent = false, we have discovered an inconsistency. In this
// case the internal state shall be considered broken.
// Otherwise, the internal state contains the following data:
//  - If known = true, then the length is known to be len.
//  - The length belongs to [lower_bound, upper_bound].
struct MetricDart {
    bool consistent;
    bool known;
    double len;
    double lower_bound, upper_bound; // Strict inequalities.
    MetricDart(): 
        consistent(true), known(false), lower_bound(0), upper_bound(2*M_PI) {}
    MetricDart(double _len): MetricDart() {
        known = true;
        len = _len;
    }
    
    void setLen(double _len, double eps = 1e-5) { 
        if (known) consistent &=  abs(len-_len) < eps;
        else if (_len <= lower_bound or _len >= upper_bound) consistent = false;
        else {
            known = true;
            len = _len; 
        }
    }
    
    double getLen() const {
        return known ? len : -1;
    }
    
    void setLowerBound(double lb) {
        assert(0 <= lb and lb <= 2*M_PI);
        if (known) consistent &= len > lb;
        consistent &= lb < upper_bound;
        lower_bound = max(lb, lower_bound);
    }
    
    void setUpperBound(double ub) {
        assert(0 <= ub and ub <= 2*M_PI);
        if (known) consistent &= len < ub;
        consistent &= lower_bound < ub;
        upper_bound = min(ub, upper_bound);
    }
};

// Keep track of return value.
MetricDart merge(MetricDart d, MetricDart e) {
    if (e.known)
        d.setLen(e.len);
    d.setLowerBound(e.lower_bound);
    d.setUpperBound(e.upper_bound);
    
    return d;
}

class MetricDerivedOnceNet : public DerivedOnceNet {
  private:
    // We identify darts if we can show that they have the same length.
    // _label is the low-level vector used to implement the associated
    // DSU structure.
    mutable vector<int> _label;
    vector<MetricDart> _metric_darts;
    bool consistent;
    
    // Implementation of DSU.
    int ancestor(int x) const {
        if (_label[x] != x) 
            _label[x] = ancestor(_label[x]);
        return _label[x];
    }
    
    void identifyDarts(const Dart& d, const Dart& e) {
        int x = ancestor(d.id), y = ancestor(e.id);
        if (x != y) {
            _label[x] = y;
            _metric_darts[y] = merge(_metric_darts[x], _metric_darts[y]);
        }
    }
    
    bool equalLength(const Dart& d, const Dart& e) const { 
        return ancestor(d.id) == ancestor(e.id); 
    }
    
    MetricDart& metricDart(const Dart& d) { 
        return _metric_darts[ancestor(d.id)];
    }
    
    const MetricDart& metricDart(const Dart& d) const { 
        return _metric_darts[ancestor(d.id)];
    }
    
    void setLength(const Dart& d, double len) { 
        metricDart(d).setLen(len); 
    }
    
    double getLength(const Dart& d) const { 
        return metricDart(d).getLen();
    }
    
    void findRegularFaces(const Dart& d) {
        static const vector<double> regular_sides = {-1, -1, 
            M_PI, // biangle
            acos(-1.0/3), // triangle
            2*atan(sqrt(0.5)), // quadrilateral
            acos(sqrt(5)/3)}; // pentagon
        if (!isConvexFace(d))
            return;
            
        bool is_regular = true;
        for (const Dart& e: face(d))
            is_regular &= equalLength(d, e);
        if (!is_regular)
            return;
        
        
        int n = ssize(face(d));
        for (const Dart& e: face(d))
            setLength(d, regular_sides[n]);
    }
    
    void findSymmetries(const Dart& d) {
        if (!isConvexFace(d))
            return;
        
        int n = ssize(face(d));
        if (n <= 3) { // biangles and triangles are regular.
            for (int i = 1; i < n; i++)
                identifyDarts(d, phi(d, i));
        } else if (n == 4) { // quadrilaterals are symmetric.
            identifyDarts(d, phi(d, 2));
            identifyDarts(phi(d, 1), phi(d, 3));
            return;
        } else if (n == 5) { // TODO: Maybe handle symmetric pentagons?
            if (equalLength(d, phi(d))) 
                identifyDarts(phi(d, 2), phi(d, 4));
            if (equalLength(d, phi(d, 2))) 
                identifyDarts(phi(d, 3), phi(d, 4));
        }
    }
    
    void determineFace(const Dart& d) {
        int n = ssize(face(d));
        if (n <= 2)
            return;
        
        vector<double> angles = faceAngles(d);        
        vector<double> sides;
        vector<int> missing_sides;
        for (int i = 0; i < n; i++) {
            const Dart& e = phi(d, i);
            sides.push_back(getLength(e));
            if (sides.back() < 0)
                missing_sides.push_back(i);
        }
        
        if (ssize(missing_sides) > 3) // Cannot be determined.
            return;
        
        // Faking some missing sides until we have exactly three.
        for (int i = 0; i < n and ssize(missing_sides) < 3; i++) {
            if (sides[i] < 0)
                missing_sides.push_back(i);
        }
        
        vector<vector<double>> candidates = 
            completeGeneralizedPolygon(angles, sides, missing_sides);
            
        if (ssize(candidates) == 3)
            return; // Cannot be determined.
        
        vector<vector<double>> valid_faces;
        for (const auto& candidate: candidates) {
            auto poly = constructGeneralizedPolygon(angles, candidate);
            if (poly.empty())
                continue;

            if (isValidPolygon(poly) and checkPinProperties(poly))
                valid_faces.push_back(candidate);
        }
        // It turns out that there is always at most one candidate, even
        // though there could be two a priori.
        assert(ssize(valid_faces) <= 1);
        if (valid_faces.empty()) {
            consistent = false;
            return;
        }
        
        sides = valid_faces[0];
        for (int i = 0; i < n; i++)
            setLength(phi(d, i), sides[i]);
    }
    
    void checkMetricDartsConsistency() {
        for (const Dart& d: darts())
            consistent &= metricDart(d).consistent;
    }
    
  public:
    MetricDerivedOnceNet(const DerivedOnceNet& A): DerivedOnceNet(A) {
        consistent = true;
        // Init labels.
        for (const Dart& d: darts()) {
            assert(d.id == ssize(_label));
            _label.push_back(d.id);
            _metric_darts.push_back(MetricDart());
        }
        
        for (const Dart& d: darts())
            identifyDarts(d, rev(d));
            
        for (int _ = 0; _ < 5; _++)
            for (const Dart& d: darts())
                findSymmetries(d);
        
        if (!consistent)
            return;
        
        for (const Dart& d: darts())
            findRegularFaces(d);
        
        if (!consistent)
            return;
        
        for (int _ = 0; _ < 10 and consistent; _++)
            for (const Dart& d: darts())
                determineFace(d);
    }
};


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
        vector<double> angles = A.faceAngles(d);
        
        int n = A.faceType(d).first;
        for (int i = 0; i < n; i++) {
            Dart e = A.phi(d, i);
            sides.push_back(get_length(e));
            if (get_length(e) < 0)
                missing_sides.push_back(i);
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
