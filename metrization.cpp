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
    acos(sqrt(5)/3)}; // pentagon, 0.72972765622696...
    
// The length a = 58.25684287.. given by Taylor at page 502.
// If x = tan(a/2), then 8x^6+39x^4-4=0. 
static const double QUAD_PENT = 1.0167737169750890006208651870679862;

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
    
    void setLen(double _len) { 
        if (_len < 1e-6)
            consistent = false;
        
        if (known) consistent &=  abs(len-_len) < eps;
        else {
            known = true;
            len = _len; 
            setLowerBound(_len-eps);
            setUpperBound(_len+eps);
        }
    }
};


vector<int> findEqualLengths(const DerivedOnceNet& A) {
    int n = 2*A.E();
    
    // BEGIN: Implementation of DSU.
    vector<int> label(n);
    for (int i = 0; i < n; i++) label[i] = i;
    
    function<int(int)> ancestor = [&](int x) {
        return (label[x] == x) ? x : (label[x] = ancestor(label[x]));
    };
    
    auto are_equal = [&](const Dart& d, const Dart& e) {
        return ancestor(d.id) == ancestor(e.id);
    };
    
    auto set_equal = [&](const Dart& d, const Dart& e) {
        int x = ancestor(d.id), y = ancestor(e.id);
        if (x != y) 
            label[x] = y;
    };
    // END: Implementation of DSU.
    
    // The ideas are all contained in the following function.
    // We employ multiple simple strategy to prove that two different
    // darts shall have the same length.
    auto find_simmetries = [&](Dart d, map<int,vector<Dart>>& partition) {
        if (!A.isConvexFace(d))
            return;
            
        int h = (ancestor(d.id) + 1) * 1000;
        
        int n = A.faceType(d).first;
        if (n <= 3) { // biangles and triangles are regular.
            for (int i = 0; i < n; i++)
                partition[n].push_back(d);
        } else if (n == 4) { // quadrilaterals are symmetric.
            partition[4 + h].push_back(A.phi(d));
            partition[4 + h].push_back(A.phi(d, 3));
            // square
            if (are_equal(d, A.phi(d))) partition[4].push_back(d);
        } else if (n == 5) { // TODO: Maybe handle non-convex pentagons?
            if (are_equal(d, A.phi(d))) {
                partition[5 + h].push_back(A.phi(d, 2));
                partition[5 + h].push_back(A.phi(d, 4));
            }
            if (are_equal(d, A.phi(d, 2))) {
                partition[6 + h].push_back(A.phi(d, 3));
                partition[6 + h].push_back(A.phi(d, 4));
            }
            if (are_equal(d, A.phi(d)) and are_equal(d, A.phi(d, 2)))
                partition[5].push_back(d);
        }
    };
    
    for (const Dart& d: A.darts())
        set_equal(d, A.rev(d));
    
    while (1) {
        vector<int> bef_label = label;
        map<int,vector<Dart>> partition;
        for (const Dart& d: A.darts())
            find_simmetries(d, partition);
        
        for (auto _: partition) {
            vector<Dart> equivalence_class = _.second;
            const Dart& repr = equivalence_class[0];
            for (const Dart& d: equivalence_class)
                set_equal(repr, d);
        }
        
        if (label == bef_label)
            break;
    }
    
    for (int i = 0; i < n; i++)
        ancestor(i);
    
    return label;
}


// TODO: This shall expose some decent APIs.
class MetricDerivedOnceNet : public DerivedOnceNet {
  private:
    // If consistent = false, then we have proven that the net cannot
    // be embedded in the sphere.
    bool consistent;
    
    // If _label[d.id] == _label[e.id], then d and e have necessarily
    // the same length.
    vector<int> _label;
    vector<MetricDart> _metric_darts;
    
    bool equalLength(const Dart& d, const Dart& e) const {
        return _label[d.id] == _label[e.id]; 
    }
    
    MetricDart& metricDart(const Dart& d) { 
        return _metric_darts[_label[d.id]];
    }
    
    const MetricDart& metricDart(const Dart& d) const { 
        return _metric_darts[_label[d.id]];
    }
    
    // If _len_darts[d.id] > 0, then we are able to determine that
    // d has length _len_darts[d.id].
    // It contains the same information as metricDart(d).len, but it is
    // conceptually easier to access.
    vector<double> _len_darts;
    
    double dartLen(const Dart& d) const { return _len_darts[d.id]; }
    bool isLenKnown(const Dart& d) const { return _len_darts[d.id] > 0; }
    
    void determineRegularFace(const Dart& d) {
        if (!isConvexFace(d))
            return;
            
        bool is_regular = true;
        for (const Dart& e: face(d))
            is_regular &= equalLength(d, e);
        if (!is_regular)
            return;
        
        int n = ssize(face(d));
        metricDart(d).setLen(REGULAR_SIDES[n]);
    }
    
    void determineFace(const Dart& d) {
        int n = ssize(face(d));
        if (n <= 2)
            return;
        
        vector<double> angles = faceAngles(d);        
        vector<double> sides;
        vector<int> missing;
        vector<int> known;
        for (int i = 0; i < n; i++) {
            const MetricDart& e = metricDart(phi(d, i));
            sides.push_back(e.len);
            if (e.known)
                known.push_back(i);
            else
                missing.push_back(i);
        }
        
        if (ssize(missing) > 3) // Cannot be determined.
            return;
    
        // Faking some missing sides until we have exactly three.
        while (ssize(missing) < 3) {
            missing.push_back(known.back());
            known.pop_back();
        }
        assert(ssize(missing) == 3);
        
        vector<vector<double>> candidates = 
            completeGeneralizedPolygon(angles, sides, missing);
            
        if (ssize(candidates) == 3) {
            auto is_missing = [&](int id) {
                return find(missing.begin(), missing.end(), id % n)
                       != missing.end();
            };
            // The lengths cannot be determined, but if the three missing
            // sides are adjacent, the central one is long either
            // 0 or pi.
            // The 0 case leads to consistent = false.
            int pos = -1;
            for (int i = 0; i < n; i++) {
                if (is_missing(i) and is_missing(i+1) and is_missing(i+2))
                    pos = i;
            }
            
            if (pos == -1)
                return;
            
            double sum_angles = angles[pos] + angles[(pos + 1) % n];
            if (abs(sum_angles - 2*M_PI) < 0.1)
                consistent = false;
            else
                metricDart(phi(d, pos+1)).setLen(M_PI);
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
    
    void enforceSideInequalities(const Dart& d) {
        if (faceType(d) == pair<int,int>{5, 0}) {
            metricDart(d).setUpperBound(REGULAR_SIDES[3]-eps);
        }
        if (faceType(d) == pair<int,int>{5, 1} and isConcaveAngle(d)) {
            for (int i: {0, 1, 2, 4}) // TODO: We prove a bit less
                metricDart(phi(d, i)).setUpperBound(REGULAR_SIDES[3] - eps);
            metricDart(phi(d, 3)).setLowerBound(REGULAR_SIDES[3] + eps);
        }
        if (faceType(d) == pair<int,int>{5, 0} and equalLength(d, phi(d))) {
            metricDart(d).setUpperBound(REGULAR_SIDES[4] - eps);
            metricDart(phi(d)).setUpperBound(REGULAR_SIDES[4] - eps);
        }
        // This part of the code does not exclude a single net.
        if (faceType(d) == pair<int,int>{4, 0}) {
            double l = metricDart(d).lower_bound;
            double u = metricDart(d).upper_bound;
            if (0 < l and l < M_PI)
                metricDart(phi(d)).setUpperBound(2*atan(0.5 / tan(l/2)) + eps);
            if (0 < u and u < M_PI)
                metricDart(phi(d)).setLowerBound(2*atan(0.5 / tan(u/2)) - eps);
        }
        // It seems like Conjecture 21.2 is not useful.
        //~ if (faceAnglesDeg(d) == vector<int>{120, 120, 240, 120, 240}) {
            //~ metricDart(d).setUpperBound(REGULAR_SIDES[3]-eps);
            //~ metricDart(phi(d, 2)).setUpperBound(REGULAR_SIDES[3]-eps);
        //~ }
    }
    
    void embedDart(const Dart& d, TangentVector tv, 
                   vector<GeodesicSegment>& embedded_darts) {
        GeodesicSegment& ed = embedded_darts[d.id];
        
        if (!isLenKnown(d))
            return;
        if (ed.len >= 0) { // Check ed == tv
            consistent &= norm(ed.P - tv.P) < eps 
                          and norm(ed.V - tv.V) < eps;
            return;
        }
        
        double len = dartLen(d);
        ed = {tv, len};
        TangentVector rev_d = tv.exp(len);
        rev_d.V = -rev_d.V;
        embedDart(rev(d), rev_d, embedded_darts);
        
        if (degree(d.from) == 3) {
            embedDart(next(d), tv.rot(-2*M_PI/3), embedded_darts);
            embedDart(next(next(d)), tv.rot(2*M_PI/3), embedded_darts);
        } else if (degree(d.from == 2)) {
            if (isConcaveAngle(rev(d))) 
                embedDart(next(d), tv.rot(2*M_PI/3), embedded_darts);
            else 
                embedDart(next(d), tv.rot(-2*M_PI/3), embedded_darts);
        }
    } 
    
  public:
    bool isConsistent() { 
        for (const Dart& d: darts())
            consistent &= metricDart(d).consistent;
        return consistent; 
    }
  
    MetricDerivedOnceNet(const DerivedOnceNet& A): DerivedOnceNet(A) {
        consistent = true;
        
        // Step 1: Identify equal sides.
        _label = findEqualLengths(*this);
        _metric_darts.resize(2*E());
        
        // Step 2: Determine the length of some sides.
        for (const Dart& d: darts())
            determineRegularFace(d);
        
        for (int _ = 0; _ < 12 and isConsistent(); _++)
            for (const Dart& d: darts())
                determineFace(d);
                
        if (!isConsistent())
            return;
            
        // Step 3: Enforce all inequalities we have proven on sides
        //         depending on the type of the face they belong to.
        for (int _ = 0; _ < 10 and isConsistent(); _++)
            for (const Dart& d: darts())
                enforceSideInequalities(d);
            
        if (!isConsistent())
            return;
            
        // Step 4: Embed a part of the net (i.e., the darts whose length
        //         we were able to determine in step 2) in the sphere.
        _len_darts.resize(2*E(), -1);
        for (const Dart& d: darts()) {
            if (metricDart(d).known) {
                _len_darts[d.id] = metricDart(d).len;
            }
        }
        
        vector<bool> visited_darts(2*E(), false);
        vector<vector<GeodesicSegment>> embedded_components;
        vector<vector<Dart>> active_components;
        
        for (const Dart& d: darts()) {
            if (visited_darts[d.id] or dartLen(d) < 0)
                continue;
            vector<GeodesicSegment> embedded_darts(2*E());
            vector<Dart> active_darts;
            embedDart(d, TangentVector(), embedded_darts);
            for (int i = 0; i < 2*E(); i++) {
                if (embedded_darts[i].len > 0) {
                    visited_darts[i] = true;
                    active_darts.push_back(dart(i));
                }
            }
            
            if (!isConsistent())
                return;
            embedded_components.push_back(embedded_darts);
            active_components.push_back(active_darts);
        }
        
        // Step 5: Check that all embedded components are non-self-intersecting.
        for (int _ = 0; _ < ssize(embedded_components); _++) {
            const vector<GeodesicSegment>& embedded_darts = embedded_components[_];
            const vector<Dart>& active_darts = active_components[_];
            
            for (const Dart& x: active_darts) for (const Dart& y: active_darts) {
                if (x.from > x.to)
                    continue;
                if (y.from > y.to)
                    continue;
                if (x.from == y.from and x.to == y.to)
                    continue;
                const GeodesicSegment& ex = embedded_darts[x.id];
                const GeodesicSegment& ey = embedded_darts[y.id];
                if (x.from == y.from or x.from == y.to or
                    x.to == y.from or x.to == y.to) {
                    if (ex.len > M_PI + eps and ey.len > M_PI + eps)
                        consistent = false;
                } else if (areIntersecting(ex, ey))
                    consistent = false;
            }
            
            if (!isConsistent())
                return;
        }
        
        // Step 6: Look for inconsistencies involving the two pins.
        for (int _ = 0; _ < ssize(embedded_components); _++) {
            const vector<GeodesicSegment>& embedded_darts = embedded_components[_];
            const vector<Dart>& active_darts = active_components[_];
               
            vector<pt> equators_with_pin;
            for (const Dart& d: active_darts) { 
                if (isConcaveAngle(d)) {
                    GeodesicSegment ed = embedded_darts[rev(d).id];
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
            
            // TODO: For some reason, we are never able to determine the pin positions! This is rather strange.
            assert(pins.empty());
            
            if (!pins.empty()) {
                for (pt eq: equators_with_pin) 
                    consistent &= areOrthogonal(pins[0], eq);
            }
            
            auto is_face_determined = [&](const Dart& d) {
                bool determined = true;
                for (Dart e: face(d))
                    determined &= dartLen(e) > 0;
                return determined;
            };
            
            auto is_good_region = [&](const Dart& d) {
                return is_face_determined(d) and isConvexFace(d);
            };
            
            auto good_boundary = [&](const Dart& d) {
                return is_good_region(d) and !is_good_region(rev(d));
            };
            
            for (const Dart& d: active_darts) { 
                if (isConcaveAngle(d) and isConvexFace(rev(d)) 
                    and is_face_determined(rev(d))) {
                    bool intersection = false;
                    TangentVector ed = embedded_darts[rev(d).id].rot(M_PI/3);
                     // TODO: Bad! This is ugly, but doing it properly
                     //       requires some thought.
                     //       Needs to be fixed asap.
                    GeodesicSegment gs(ed, M_PI-1e-3);
                    for (const Dart& e: active_darts) {
                        if (e.from == d.to or e.to == d.to)
                            continue;
                        GeodesicSegment gs2 = embedded_darts[e.id];
                        gs2.len += 1e-3;
                        
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
