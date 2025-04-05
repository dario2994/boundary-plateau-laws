#ifndef METRIZATION
#define METRIZATION

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

#include "derived_once_nets.cpp"
#include "spherical_geometry.cpp"

static const vector<double> REGULAR_SIDES = {-1, -1, 
    M_PI, // biangle
    acos(-1.0/3), // triangle
    2*atan(sqrt(0.5)), // quadrilateral
    acos(sqrt(5)/3)}; // pentagon

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
    
    void setLen(double _len) { 
        if (_len < 1e-6)
            consistent = false;
        
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
    if (!e.consistent)
        d.consistent = false;
    if (e.known)
        d.setLen(e.len);
    d.setLowerBound(e.lower_bound);
    d.setUpperBound(e.upper_bound);
    
    return d;
}


// TODO: This shall expose some decent APIs.
class MetricDerivedOnceNet : public DerivedOnceNet {
  private:
    // We identify darts if we can show that they have the same length.
    // _label is the low-level vector used to implement the associated
    // DSU structure.
    mutable vector<int> _label;
    vector<MetricDart> _metric_darts;
    vector<GeodesicSegment> _embedded_darts;
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
    
    void findRegularFaces(const Dart& d) {
        if (!isConvexFace(d))
            return;
            
        bool is_regular = true;
        for (const Dart& e: face(d))
            is_regular &= equalLength(d, e);
        if (!is_regular)
            return;
        
        
        int n = ssize(face(d));
        for (const Dart& e: face(d))
            metricDart(d).setLen(REGULAR_SIDES[n]);
    }
    
    void findSymmetries(const Dart& d, map<int,vector<Dart>>& partition) {
        if (!isConvexFace(d))
            return;
            
        int h = (ancestor(d.id)+1) * 1000;
        
        int n = ssize(face(d));
        if (n <= 3) { // biangles and triangles are regular.
            for (int i = 0; i < n; i++)
                partition[n].push_back(d);
        } else if (n == 4) { // quadrilaterals are symmetric.
            partition[4 + h].push_back(phi(d));
            partition[4 + h].push_back(phi(d, 3));
            if (equalLength(d, phi(d))) partition[4].push_back(d);
        } else if (n == 5) { // TODO: Maybe handle non-convex pentagons?
            if (equalLength(d, phi(d))) {
                partition[5 + h].push_back(phi(d, 2));
                partition[5 + h].push_back(phi(d, 4));
            }
            if (equalLength(d, phi(d, 2))) {
                partition[6 + h].push_back(phi(d, 3));
                partition[6 + h].push_back(phi(d, 4));
            }
            if (equalLength(d, phi(d)) and equalLength(d, phi(d, 2)))
                partition[5].push_back(d);
        }
    }
    
    void generalSideBounds(const Dart& d) {
        if (faceType(d) == pair<int,int>{5, 0}) {
            for (const Dart& e: face(d))
                metricDart(d).setUpperBound(REGULAR_SIDES[3]);
        }
        if (faceType(d) == pair<int,int>{5, 1} and isConcaveAngle(d)) {
            for (int i: {0, 1, 2, 4}) // TODO: We prove a bit less
                metricDart(phi(d, i)).setUpperBound(REGULAR_SIDES[3]);
            metricDart(phi(d, 3)).setLowerBound(REGULAR_SIDES[3]);
        }
        if (faceType(d) == pair<int,int>{5, 0} and equalLength(d, phi(d))) {
            metricDart(d).setUpperBound(REGULAR_SIDES[4]);
            metricDart(phi(d)).setUpperBound(REGULAR_SIDES[4]);
        }
    }
    
    // This part of the code does not exclude a single net.
    void quadrilateralSideBounds(const Dart& d) {
        if (faceType(d) != pair<int,int>{4, 0})
            return;
        double l = metricDart(d).lower_bound;
        double u = metricDart(d).upper_bound;
        if (0 < l and l < M_PI)
            metricDart(phi(d)).setUpperBound(2*atan(0.5 / tan(l/2)) + eps);
        if (0 < u and u < M_PI)
            metricDart(phi(d)).setLowerBound(2*atan(0.5 / tan(u/2)) - eps);
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
            sides.push_back(metricDart(e).getLen());
            if (sides.back() < 0)
                missing_sides.push_back(i);
        }
        
        if (ssize(missing_sides) > 3) // Cannot be determined.
            return;
        
        // Faking some missing sides until we have exactly three.
        for (int i = 0; i < n and ssize(missing_sides) < 3; i++) {
            if (sides[i] >= 0)
                missing_sides.push_back(i);
        }
        assert(ssize(missing_sides) == 3);
        
        vector<vector<double>> candidates = 
            completeGeneralizedPolygon(angles, sides, missing_sides);
            
        if (ssize(candidates) == 3) {
            auto is_missing = [&](int id) {
                id %= n;
                bool found = false;
                for (int x: missing_sides)
                    found |= x == id;
                return found;
            };
            // The lengths cannot be determined, but if the missing
            // sides are adjacent, the central one is long either
            // 0 or pi.
            // The 0 case leads to consistent = false.
            int first_missing = -1;
            for (int i = 0; i < n; i++) {
                if (is_missing(i) and is_missing(i+1) and is_missing(i+2))
                    first_missing = i;
            }
            if (first_missing == -1)
                return;
            
            double sum_angles = angles[first_missing] + 
                                angles[(first_missing + 1) % n];
            if (abs(sum_angles - 2*M_PI) < 0.1) {
                consistent = false;
            } else {
                metricDart(phi(d, first_missing+1)).setLen(M_PI);
            }
            return;
        }
        
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
            metricDart(phi(d, i)).setLen(sides[i]);
    }
    
    void checkMetricDartsConsistency() {
        for (const Dart& d: darts())
            consistent &= metricDart(d).consistent;
    }
    
    // Temporary storage for the dfs of embedDart
    vector<Dart> _visited_darts;
    void embedDart(const Dart& d, TangentVector tv) {
        MetricDart& md = metricDart(d);
        GeodesicSegment& ed = _embedded_darts[d.id];
        
        if (!md.known)
            return;
        if (ed.len >= 0) { // Check == tv
            if (norm(ed.P - tv.P) > eps or
                norm(ed.V - tv.V) > eps)
                md.consistent = false;
            
            return;
        }
        _visited_darts.push_back(d);
        ed = {tv, md.len};
        TangentVector rev_d = tv.exp(md.len);
        rev_d.V = -rev_d.V;
        embedDart(rev(d), rev_d);
        
        if (degree(d.from) == 3) {
            embedDart(next(d), tv.rot(-2*M_PI/3));
            embedDart(next(next(d)), tv.rot(2*M_PI/3));
        } else if (degree(d.from == 2)) {
            if (isConcaveAngle(rev(d))) embedDart(next(d), tv.rot(2*M_PI/3));
            else embedDart(next(d), tv.rot(-2*M_PI/3));
        }
    } 
    
  public:
    bool isConsistent() { 
        checkMetricDartsConsistency();
        return consistent; 
    }
  
    MetricDerivedOnceNet(const DerivedOnceNet& A): DerivedOnceNet(A) {
        consistent = true;
        // Init labels.
        for (const Dart& d: darts()) {
            assert(d.id == ssize(_label));
            _label.push_back(d.id);
            _metric_darts.push_back(MetricDart());
        }
        
        // TODO: This code needs to be cleaned.
        
        for (const Dart& d: darts())
            identifyDarts(d, rev(d));
            
        for (int _ = 0; _ < 10; _++) {
            map<int,vector<Dart>> partition;
            for (const Dart& d: darts())
                findSymmetries(d, partition);
            for (auto pp: partition)
                for (Dart d: pp.second)
                    identifyDarts(pp.second[0], d);
        }
        
        if (!isConsistent())
            return;
        
        for (const Dart& d: darts())
            findRegularFaces(d);
        
        if (!isConsistent())
            return;
        
        for (int _ = 0; _ < 12 and isConsistent(); _++)
            for (const Dart& d: darts())
                determineFace(d);
        
        if (!isConsistent())
            return;
            
        for (const Dart& d: darts())
            generalSideBounds(d);
        
        if (!isConsistent())
            return;
            
        for (int _ = 0; _ < 3 and isConsistent(); _++)
            for (const Dart& d: darts())
                quadrilateralSideBounds(d);
        
        if (!isConsistent())
            return;
                
        _embedded_darts.resize(2*E());
                
        for (Dart uff: darts()) {
            if (_embedded_darts[uff.id].len > 0 or !metricDart(uff).known)
                continue;
            _visited_darts.clear();
            embedDart(uff, TangentVector());
            
            if (!isConsistent())
                return;
            
            // Check non-intersection of the embedded darts.
            for (Dart x: _visited_darts) for (Dart y: _visited_darts) {
                if (x.from > x.to)
                    continue;
                if (y.from > y.to)
                    continue;
                if (x.from == y.from and x.to == y.to)
                    continue;
                GeodesicSegment& ex = _embedded_darts[x.id];
                GeodesicSegment& ey = _embedded_darts[y.id];
                if (x.from == y.from or x.from == y.to or
                    x.to == y.from or x.to == y.to) {
                    if (ex.len > M_PI + eps and ey.len > M_PI + eps)
                        consistent = false;
                } else if (areIntersecting(ex, ey))
                    consistent = false;
            } // Shall we also print that?
            
            if (!isConsistent())
                return;
            
            vector<pt> equators_with_pin;
            for (Dart d: _visited_darts) { 
                if (isConcaveAngle(d)) {
                    GeodesicSegment ed = _embedded_darts[rev(d).id];
                    equators_with_pin.push_back(ed.rot(-2*M_PI/3).equatorNormal());
                }
            }
            
            vector<pt> pins;
            for (int i = 1; i < ssize(equators_with_pin); i++) {
                pt eq0 = equators_with_pin[0];
                pt eq1 = equators_with_pin[i];
                if (norm(eq0-eq1) > eps and norm(eq0+eq1) > eps) {
                    pt pin = normalize(eq0^eq1);
                    pins = {pin, -pin};
                    break;
                }
            }
            if (!pins.empty()) {
                for (pt eq: equators_with_pin) {
                    if (!areOrthogonal(pins[0], eq)) {
                        consistent = false;
                    }
                }
            }
            
            auto is_face_determined = [&](Dart d) {
                bool determined = true;
                for (Dart e: face(d))
                    determined &= metricDart(e).known;
                return determined;
            };
            
            auto is_good_region = [&](Dart d) {
                return is_face_determined(d) and isConvexFace(d);
            };
            
            auto good_boundary = [&](Dart d) {
                return is_good_region(d) and !is_good_region(rev(d));
            };
            
            // TODO: For some reason, we are never able to determine the pin positions! This is rather strange.
            assert(pins.empty());
            
            //~ cout << ssize(_visited_darts) << " " << ssize(darts()) << endl;
            //~ assert(ssize(_visited_darts) == ssize(darts()));
            
            for (Dart d: _visited_darts) { 
                if (isConcaveAngle(d) and isConvexFace(rev(d)) 
                    and is_face_determined(rev(d))) {
                    //~ cout << d << endl;
                    bool intersection = false;
                    TangentVector ed = _embedded_darts[rev(d).id].rot(M_PI/3);
                    //~ cout << _embedded_darts[rev(d).id] << endl;
                    //~ cout << _embedded_darts[next(rev(d)).id] << endl;
                    //~ cout << ed << endl;
                     // TODO: Bad! This is ugly, but doing it properly
                     //       requires some thought.
                     //       Needs to be fixed asap.
                    GeodesicSegment gs(ed, M_PI-1e-3);
                    for (Dart e: _visited_darts) {
                        if (e.from == d.to or e.to == d.to)
                            continue;
                        GeodesicSegment gs2 = _embedded_darts[e.id];
                        gs2.len += 1e-3;
                        
                        //~ if (good_boundary(e))
                            //~ cout << e << endl;
                        
                        if (good_boundary(e) and areIntersecting(gs2, gs))
                            intersection = true;
                    }
                    if (!intersection)
                        consistent = false;
                }
            }
        }
    }
};


#endif
