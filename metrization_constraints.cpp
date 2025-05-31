#ifndef METRIZATION_CONSTRAINTS
#define METRIZATION_CONSTRAINTS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "interval_arithmetic.cpp"
#include "quaternions.cpp"
#include "derived_once_nets.cpp"
#include "metrized_nets.cpp"
#include "utils.cpp"

// Utility function to enforce the constraint on the sides of a 
// quadrilateral.
// 
// I is strictly inside [0, \pi].
// One side of a quadrilateral belongs to I.
// The other side belongs to the interval returned by this function.
Interval quadrilateralConstraint(Interval I) {
    if (I.empty())
        return Interval(0, -1);
    const static double eps = 1e-10;
    static auto f = [](double x) { return 2*atan(0.5 / tan(x/2)); };
    return Interval(f(I.r) - eps, f(I.l) + eps);
}

// Utility function to enforce the constraint on the sides of a 
// pentagon.
// 
// I, J are strictly inside [0, \ell_T].
// Let ABCDE be a spherical pentagon with angles 120.
// If AB \in I, BC\in J, then DE belongs to the interval returned by 
// this function.
Interval pentagon1Constraint(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
    
    static auto f = [](double x, double y) {
        return 0.5*(1.0/3 + 
            cos(x) + cos(y) + cos(x)*cos(y) - 2*sin(x)*sin(y));
    };
    
    double high = f(I.l, J.l);
    double low = f(I.r, J.r);
    
    return acos(Interval(low, high)); 
}

// Utility function to enforce the constraint on the sides of a 
// pentagon.
// 
// I, J are strictly inside [0, \ell_T].
// Let ABCDE be a spherical pentagon with angles 120.
// AB \in I, DE\in J.
// Then BC belongs to the interval returned by this function.
Interval pentagon2Constraint(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
        
    static auto get_theta = [](double x) {
        return atan(2*sin(x) / (1+cos(x)));
    };
    
    static auto get_lambda = [](double x) {
        return 2*sin(x) / sin(get_theta(x));
    };
    
    static auto f = [&](double x, double z) {
        return (2*cos(z)-1.0/3-cos(x))/get_lambda(x);
    };
    
    double low = f(I.l, J.r);
    double high = f(I.r, J.l);
    
    // Applying acos on [low, high] is wrong here, hence we do the 
    // computation explicitly. TODO: Check.
    
    if (low > 1 or high < -1)
        return Interval(0, -1); // empty
    
    double l = (high < 1) ?  acos(high) : 0;
    l -= get_theta(I.r);
    double r = (low > -1) ? acos(low) : M_PI;
    r -= get_theta(I.l);
    Interval K(l, r);
    
    const static double ellT = acos(-1.0/3);
    K.intersect(Interval(0, ellT));
    
    return K; 
}


bool refineMetrizationConvexContinent(MetrizedNet& net, 
                                      const vector<Dart>& initial_darts) {
    StopWatch _("refineMetrizationConvexContinent");
        
    const DerivedOnceNet& A = net.A;
    
    vector<bool> to_propagate(2*A.E(), false);
    queue<Dart> propagation_queue; 
    for (const Dart& d: initial_darts) {
        int id = d.id;
        propagation_queue.push(d);
        propagation_queue.push(A.rev(d));
        to_propagate[id] = true;
        to_propagate[id+1] = true;
    }
        
    bool empty_interval = false;
    auto update_length = [&](const Dart& d, const Interval& I) {
        if (net.length(d).empty()) {
            empty_interval = true;
            return;
        }
        if (I.empty()) {
            net.length(d) = I;
            empty_interval = true;
            return;
        }
        double bef = net.length(d).len();
        net.length(d).intersect(I);
        double aft = net.length(d).len();
        
        bool updated = bef > 1e-12 and aft < 0.9*bef;
        if (updated) {
            if (!to_propagate[d.id]) {
                to_propagate[d.id] = true;
                propagation_queue.push(d);
            }
            if (!to_propagate[A.rev(d).id]) {
                to_propagate[A.rev(d).id] = true;
                propagation_queue.push(A.rev(d));
            }
        }
    };
    
    while (!empty_interval and !propagation_queue.empty()) {
        const Dart d = propagation_queue.front();
        propagation_queue.pop();
        
        if (!to_propagate[d.id] or !A.isConvexFace(d))
            continue;
        to_propagate[d.id] = false;
        
        const vector<Dart>& f = A.face(d);
        
        if (ssize(f) == 4) { // Quadrilateral
            update_length(f[2], net.length(f[0]));
            const Interval& I = quadrilateralConstraint(net.length(f[0]));
            update_length(f[1], I);
            update_length(f[3], I);
        }
        if (ssize(f) == 5) { // Pentagon
            update_length(f[3], pentagon1Constraint(net.length(f[0]), 
                                                    net.length(f[1])));
            update_length(f[2], pentagon1Constraint(net.length(f[0]), 
                                                    net.length(f[4])));
            update_length(f[4], pentagon2Constraint(net.length(f[0]), 
                                                    net.length(f[2])));
            update_length(f[1], pentagon2Constraint(net.length(f[0]), 
                                                    net.length(f[3])));
            update_length(f[3], pentagon2Constraint(net.length(f[2]), 
                                                    net.length(f[0])));
            update_length(f[2], pentagon2Constraint(net.length(f[3]), 
                                                    net.length(f[0])));
        }
    }
    
    return net.totalError() >= 0;
}


template<typename T>
T operator*(const array<T, 3>& X, const array<T, 3>& Y) {
    return X[0]*Y[0] + X[1]*Y[1] + X[2]*Y[2];
}

template<typename T>
array<T, 3> operator-(const array<T, 3>& X, const array<T, 3>& Y) {
    return {X[0]-Y[0], X[1]-Y[1], X[2]-Y[2]};
}

template<typename T>
array<T, 3> operator^(const array<T, 3>& X, const array<T, 3>& Y) {
    return {X[1]*Y[2]-X[2]*Y[1], X[2]*Y[0]-X[0]*Y[2], X[0]*Y[1]-X[1]*Y[0]};
}

template <typename T>
T det(const array<T, 3>& A, const array<T, 3>& B, const array<T, 3>& C) {
    return A * (B ^ C);
}


// Checks that the equators corresponding to any three edges going to
// a pin are concurrent.
bool canPinsBeAntipodal(const MetrizedNet& net) {
    StopWatch _("canPinsBeAntipodal");
    
    int q = ssize(net.A.concaveAngles());
    for (int i = 0; i < q; i++) 
        for (int j = i+1; j < q; j++)
            for (int k = j+1; k < q; k++) {
                Dart root = net.A.phi(net.A.concaveAngles()[i]);
                
                vector<array<Interval, 3>> equators;
                for (int z: {i, j, k}) {
                    Dart d = net.A.concaveAngles()[z];
                    auto isom = net.getIsometryToPin(root, d);
                    equators.push_back(applyIsometry(isom, {0, 0, 1}));
                }                
                if (!det(equators[0], equators[1], equators[2]).containsZero())
                    return false;
            }
    
    return true;
}

// Checks that the angles between edges meeting at a pin belongs to
// [120, 240].
bool checkPinAngles(const MetrizedNet& net) {
    StopWatch _("checkPinAngles");
    
    const DerivedOnceNet& A = net.A;
    
    for (const auto& f: A.faces()) {
        if (A.faceType(f).second != 2) // Two pins
            continue;
        
        vector<Dart> concave_angles;
        for (Dart d:  f) {
            if (A.isConcaveAngle(d))
                concave_angles.push_back(d);
        }
        assert(ssize(concave_angles) == 2);
        
        Dart root = A.phi(concave_angles[0]);
        auto iso0 = net.getIsometryToPin(root, concave_angles[0]);
        auto iso1 = net.getIsometryToPin(root, concave_angles[1]);
        
        auto eq0 = applyIsometry(iso0, {0, 0, 1});
        auto eq1 = applyIsometry(iso1, {0, 0, 1});
        
        Interval scal_prod = eq0[0]*eq1[0] + 
                             eq0[1]*eq1[1] + 
                             eq0[2]*eq1[2];
        
        // scal_prod <= -1/2 is equivalent to the condition on
        // the angle at the pin.
        if (scal_prod.l > -0.5 + 1e-8)
            return false;
    }
    return true;
}

// If (cos x, sin(x)) \in IxJ \cup (-I)x(-J) with 0 <= x <= pi/2, 
// then x belongs to the interval returned by this function.
Interval intersectWithCircle(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
        
    // We shall return the interval of 0 <= l < pi such that:
    // (cos l/2, sin l/2) belongs to IxJ or -(IxJ).
    auto intersectWithCircleQuadrant = [](Interval X, Interval Y) {
        X.intersect(Interval(0, 1));
        Y.intersect(Interval(0, 1));
        if (X.empty() or Y.empty())
            return Interval(0, -1);
        return intersection(Interval(acos(X.r), acos(X.l)),
                            Interval(asin(Y.l), asin(Y.r)));
    };
    
    Interval inter1 = intersectWithCircleQuadrant(I, J);
    Interval inter2 = intersectWithCircleQuadrant(-I, -J);
    
    if (!inter1.empty() and !inter2.empty())
        return Interval(0, 10);
    return inter1.empty() ? inter2 : inter1;
}

// Assume that you know the lengths of k-1 sides of a non-convex face
// with k sides. This function determines (robustly) whether they
// can form the non-convex face (i.e., the angles are correct) and
// enforces the right length of the remaining side.
bool refineMetrizationNonConvexFaces(MetrizedNet& net) {
    StopWatch _("refineMetrizationNonConvexFaces");
    
    const DerivedOnceNet& A = net.A;
    
    for (Dart d: A.darts()) {
        if (A.isConvexFace(d))
            continue;
        int q = ssize(A.face(d));
        
        // Run the check only for the (almost) most undetermined side.
        bool is_maximum_dart = true;
        for (int i = 1; i < q; i++) {
            if (net.length(A.phi(d, i)).len() > 1.1 * net.length(d).len())
                is_maximum_dart = false;
        }
        if (!is_maximum_dart)
            continue;
        
        vector<PathStep> path;
        for (int i = 0; i < q; i++) {
            Dart e = A.phi(d, i);
            if (i > 0)
                path.emplace_back(e);
            path.emplace_back(PathStep('R', A.phiAngle(e)));
        }
        Quaternion<Interval> T = net.isometryFromPath(path);
        // T is, up to sign, (cos(l/2), 0, 0, -sin(l/2)), where
        // l is the length of d.
        if (!T[1].containsZero() or !T[2].containsZero()) {
            net.length(d) = Interval(0, -1);
            return false;
        }
        Interval I = 2*intersectWithCircle(T[0], -T[3]);
        net.length(d).intersect(I);
        if (net.length(d).empty())
            return false;
    }
    return true;
}

// T is an approximate isometry.
// Consider the directed geodesic a starting from (e_1, e_2) and the 
// directed geodesic b starting from (Te_1, Te_2).
//
// If this function returns l, it is guaranteed that we will intersect
// the equator defined by b going along a for at most length l.
// If halfEquator is true, then it is also guaranteed that the intersection
// happens in the "first half" of the equator (starting from b). If the
// first intersection would be on the other half of the equator, this
// function would return a large number.
double intersectWithEquator(const Quaternion<Interval>& T,
                            bool halfEquator) {
    //~ StopWatch _("farthestIntersection");
    
    auto Te3 = applyIsometry(T, {0, 0, 1}); // For efficiency.
    
    if (halfEquator) {
        // sign1 determines on which side of a is Te_1.
        int sign1 = sign(applyIsometry(T, {1, 0, 0})[2]);
        // sign2 determines on which side of b is e_1.
        int sign2 = sign(Te3[0]);
        
        // The first intersection happens in less than pi on both geodesics
        // if and only if sign1 and sign2 have different signs.
        if (sign1 * sign2 != -1)
            return 10;
    }
    
    // Te_3^e_3 is, up to scaling, the intersection between a and b.
    // Notice that Te_3^e_3 = (X, Y, 0).
    Interval X = Te3[1];
    Interval Y = -Te3[0];
    
    if (Y.containsZero())
        return 10;
        
    // We identify the *first* intersection.
    if (Y.l < 0) X = -X, Y = -Y;
    
    double theta = max(atan2(Y.l, X.l), atan2(Y.r, X.l));
    return theta;
}

bool checkNonIntersection(const Quaternion<Interval>& T, 
                          Interval& I, Interval& J) {
    //~ StopWatch _("checkNonIntersection2");
    
    double lenI = intersectWithEquator(T, true);
    if (lenI >= I.r)
        return true;
    double lenJ = intersectWithEquator(conj(T), true);
    
    if (I.l >= lenI)
        J.intersect(Interval(0, lenJ));
    if (J.l >= lenJ)
        I.intersect(Interval(0, lenI));
        
    return !I.empty() and !J.empty();
}

// Checks that non-convex faces are not self-intersecting.
bool checkNonIntersection(MetrizedNet& net) {
    StopWatch _("checkNonIntersection");
    
    const DerivedOnceNet& A = net.A;
    for (Dart d: A.darts()) {
        if (A.isConvexFace(d))
            continue;

        // We don't use net.isometryFromPath to improve the performance.
        double sum_theta = 0;
        sum_theta += M_PI;
        const vector<Dart>& f = A.face(d);
        Quaternion<Interval> iso(1);
        
        for (int i = 0; i < ssize(f)-2; i++) {
            sum_theta += A.phiAngle(f[i]);
            if (i != 0) {
                Dart e = f[i+1];
                bool ans = checkNonIntersection(
                        iso * basicRotation<Interval>(sum_theta),
                        net.length(d),
                        net.length(e));
                if (!ans)
                    return false;
            }
            iso *= basicIsometry<Interval>(sum_theta, net.length(f[i+1]));
            
            double total_error = 
                iso[0].len() + iso[1].len() + iso[2].len() + iso[3].len();
            if (total_error > 1)
                break;
        }
    }
    return true;
}

// Returns true if P belongs strictly to the given convex face, where
// we assume that the dart root is (e_1, e_2).
// If up_to_symmetry is true, then this returns true if either P or -P
// belongs strictly to the convex face. 
bool strictlyInsideConvexFace(const MetrizedNet& net, 
                        array<Interval, 3> P, 
                        const Dart& root,
                        const vector<Dart>& face,
                        bool up_to_symmetry = false) {
    assert(net.A.isConvexFace(face));
    
    vector<int> signs;
    for (const Dart& d: face) {
        const auto& Z = applyIsometry(net.getIsometry(root, d), {0, 0, 1});
        int s = sign(P * Z);
        signs.push_back(s);
        
        if (s == 0 or s != signs[0])
            break; // For performance.
    }
    bool same_sign = true;
    for (int s: signs) 
        same_sign &= s != 0 and s == signs[0];
    
    if (!same_sign)
        return false;
        
    return up_to_symmetry or signs[0] == -1;
}

// Checks that no pins fall strictly inside a convex face.
bool checkPinOutsideConvexFaces(const MetrizedNet& net) {
    StopWatch _("checkPinOutsideConvexFaces");
    
    int cnt = ssize(net.A.concaveAngles());
    for (int i = 0; i < cnt; i++) for (int j = i+1; j < cnt; j++) {
        const Dart& root = net.A.concaveAngles()[i];
        const Dart& other = net.A.concaveAngles()[j];
        const auto& X = applyIsometry(net.getIsometryToPin(root, root), {0, 0, 1});
        const auto& Y = applyIsometry(net.getIsometryToPin(root, other), {0, 0, 1});
        const auto& P = X^Y;
        for (const auto& f: net.A.faces()) {
            if (!net.A.isConvexFace(f))
                continue;
            
            // Check that none of the two intersections of X and Y lie
            // in the convex face f.
            if (strictlyInsideConvexFace(net, P, root, f, true))
                return false;
        }
    }
    return true;
}

// If it returns true, the given points belong to a strict hemisphere.
// Complexity: O(n^4), where n = ssize(vec).
bool belongToStrictHemisphere(vector<array<Interval, 3>> vec) {
    auto subroutine = [&]() {
        int n = ssize(vec);
        for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
        for (int k = j+1; k < n; k++) {
            auto normal = (vec[j]-vec[i])^(vec[k]-vec[i]);
            vector<int> signs;
            for (auto q: vec) {
                int s = sign(q * normal);
                signs.push_back(s);
            }
            bool same_sign = true;
            for (int s: signs) 
                same_sign &= s != 0 and s == signs[0];
            if (same_sign)
                return true;
        }
        return false;
    };
    if (subroutine())
        return true;
    
    // If all points in vec are coplanar, subroutine returns false.
    // To solve this issue, we add a "random" point to vec and we
    // call subroutine again.
    vec.push_back({0.84395, 0.19283, 0.96591});
    return subroutine();
}

// For each concave angles, it returns an upperbound on the length of
// the edge going from that angles into the pin.
//
// Consider an "edge to pin" s. 
// Observe that s cannot intersect any edge of the net. 
// Furthermore, consider the equator generated by another edge t, in the
// same face, that is adjacent to the other "edge to pin" (if it is present).
// Then, we know on which side of t the pin must be. If the vertex
// identified by s and the pin are on same side of t, then s cannot 
// intersect the equator described by t.
// In this way, we have determined the maximum length of s.
vector<double> maxDistanceToPin(const MetrizedNet& net) {
    const DerivedOnceNet& A = net.A;
    
    vector<double> len_to_pin;
    for (const Dart& to_pin: A.concaveAngles()) {
        double len = 10;
        for (const Dart& d: A.darts()) {
            const auto& T = net.getIsometryToPin(d, to_pin);
            
            if (intersectWithEquator(T, true) < net.length(d).l)
                len = min(len, intersectWithEquator(conj(T), true));
            
            if (!A.isConcaveAngle(d) and !A.isConcaveAngle(A.rev(d)))
                continue;
            
            if (!A.sameFace(d, to_pin))
                continue;
            
            if (to_pin.to == d.from or to_pin.to == d.to)
                continue;
            
            //~ assert(applyIsometry(T, {1, 0, 0})[2].r < 0);
            if (applyIsometry(T, {1, 0, 0})[2].r < 0)
                len = min(len, intersectWithEquator(conj(T), false));
        }
        len_to_pin.push_back(len);
    }
    
    return len_to_pin;
}

// This is a check useful only when either: 
//   - All the equators identifying the pins coincide, 
//   - At least one pin falls on an edge.
//
// Consider two "edges to pins" going to different pins.
// Let l1, l2 be the upper bounds on their lenghts. 
// And let q1, q2 be the edges extended to have length l1 and l2.
// Since p\in q1 and -p\in q2 (where p is a pin), q1 and q2 do not 
// belong to a strict hemisphere. 
//
// We check that they do not belong to a strict hemisphere, instead of 
// checking that q1 and -q2 intersect because:
//  - It involves less special cases.
//  - It is sufficient for our purposes.
bool checkPinsSpecial(const MetrizedNet& net) {
    const DerivedOnceNet& A = net.A;
    
    vector<double> len_to_pin = maxDistanceToPin(net);

    for (int i = 0; i < ssize(A.concaveAngles()); i++) {
        for (int j = i+1; j < ssize(A.concaveAngles()); j++) {
            Dart d1 = A.concaveAngles()[i], d2 = A.concaveAngles()[j];
            if (A.sameFace(d1, d2))
                continue;
            double l1 = len_to_pin[i], l2 = len_to_pin[j];
            
            if (l1 >= M_PI-1e-5 or l2 >= M_PI-1e-5)
                return true;
            
            const Dart& root = A.phi(d1);
            vector<PathStep> path1 = net.getPathToPin(root, d1);
            vector<PathStep> path2 = net.getPathToPin(root, d2);
            vector<PathStep> path1_end = path1;
            path1_end.push_back(PathStep('E', l1));
            vector<PathStep> path2_end = path2;
            path2_end.push_back(PathStep('E', l2));
            
            vector<array<Interval, 3>> vec;
            for (auto p: {path1, path2, path1_end, path2_end}) {
                auto pt = applyIsometry(net.isometryFromPath(p), {1, 0, 0});
                vec.push_back(pt);
            }
            
            if (belongToStrictHemisphere(vec))
                return false;
        }
    }
    return true;
}

// If there is exactly one non-convex face, then one of the pins
// must coincide with a degree-3 vertex of the net.
bool checkPinIsVertex(const MetrizedNet& net) {
    const DerivedOnceNet& A = net.A;
    
    if (A.concaveAngles().empty())
        return true;
    
    bool one_concave_face = true;
    for (const Dart& to_pin: A.concaveAngles())
        one_concave_face &= A.sameFace(to_pin, A.concaveAngles()[0]);
    
    if (!one_concave_face)
        return true;
        
    vector<double> len_to_pin = maxDistanceToPin(net);
    
    for (const Dart& root: A.darts()) {
        if (A.degree(root.from) < 3)
            continue;
        bool can_be_pin = true;
        for (int i = 0; i < ssize(A.concaveAngles()); i++) {
            const Dart& to_pin = A.concaveAngles()[i];
            
            auto T = net.getIsometryToPin(root, to_pin);
            auto equator = applyIsometry(T, {0, 0, 1});
            can_be_pin &= (equator * array<Interval, 3>{1, 0, 0}).containsZero();
            
            if (!can_be_pin)
                break;
            
            double l = len_to_pin[i];
            
            if (l >= M_PI-1e-5)
                continue;
            
            vector<PathStep> path_begin = net.getPathToPin(root, to_pin);
            vector<PathStep> path_end = path_begin;
            path_end.push_back(PathStep('E', l));
            
            vector<array<Interval, 3>> vec;
            for (auto p: {{}, path_begin, path_end}) {
                auto pt = applyIsometry(net.isometryFromPath(p), {1, 0, 0});
                vec.push_back(pt);
            }
            
            can_be_pin &= !belongToStrictHemisphere(vec);
            
        }
        if (can_be_pin) 
            return true;
    }
    return false;
}


// Checks that there are no two distinct darts d, e such that:
//  1. The faces corresponding to d and e are convex and with the 
//     same number of sides.
//  2. The lengths of the sides of the face of d are very close to the
//     lengths of the sides of the face of e.
//  3. The direction of the d is very close to the direction of e.
//
// This check is sufficient to rule out a few embedding with "overlapping"
// convex faces. (as in a "double cover")
bool checkFacesNotIdentical(const MetrizedNet& net) {
    auto interval_is_point = [&](const Interval& I, double z) {
        const double eps = 0.001;
        return z - eps <= I.l and I.r <= z + eps;
    };
    
    const DerivedOnceNet& A = net.A;
    
    for (const Dart& d: A.darts()) for (const Dart& e: A.darts()) {
        if (d == e) 
            continue;
            
        bool special = d.from == 11 and d.to == 4;
        special &= e.from == 5 and e.to == 1;
        
        if (!A.isConvexFace(d) or !A.isConvexFace(e))
            continue;
        if (ssize(A.face(d)) != ssize(A.face(e)))
            continue;
        
        bool identical = true;
        int q = ssize(A.face(d));
        for (int i = 0; i < q; i++) {
            const Dart& d2 = A.phi(d, i);
            const Dart& e2 = A.phi(e, i);
            const Interval& I = net.length(d2);
            const Interval& J = net.length(e2);
            if (abs(J.r-I.l) > 0.001)
                identical = false;
            if (abs(I.r-J.l) > 0.001)
                identical = false;
        }
        
        if (!identical)
            continue;
        
        Quaternion<Interval> T = net.getIsometry(d, e);
        if (!interval_is_point(T[0], 1) and !interval_is_point(T[0], -1))
            continue;
        
        for (int i = 1; i < 3; i++)
            identical &= interval_is_point(T[i], 0);
        if (identical)
            return false;
    }
    return true;
}

// No vertex can belong to the interior of a convex face.
bool checkVerticesOutsideConvexFaces(const MetrizedNet& net) {
    for (const Dart& root: net.A.darts()) {
        for (const auto& face: net.A.faces()) {
            if (!net.A.isConvexFace(face))
                continue;
            if (strictlyInsideConvexFace(net, array<Interval,3>{1, 0, 0}, root, face))
                return false;
        }
    }
    return true;
}




// Implementation of the lower bounds on edge lengths that follow from
// the stability of the net.
bool checkLengthLowerbound(MetrizedNet& net) {
    static const auto fun = [](double x, double y) {
        return min(0.01, 0.1*min(x, y)) - 1e-10;
    };
    
    const DerivedOnceNet& A = net.A;
    
    auto min_adjacent = [&](const Dart& d)->double {
        vector<double> vals = {net.length(A.next(d)).l};
        if (A.degree(d.from) == 3)
            vals.push_back(net.length(A.next(A.next(d))).l);
        return *min(vals.begin(), vals.end());
    };
    
    for (const Dart& d: A.darts()) {
        Interval& Z = net.length(d);
        double low = fun(min_adjacent(d), 
                         min_adjacent(A.rev(d)));
        Z.l = max(Z.l, low);
        if (Z.empty()) 
            return false;
    }
    
    for (const Dart& d: A.darts()) {
        if (A.degree(d.from) != 2)
            continue;
        
        
        Interval& X = net.length(d);
        Interval& Y = net.length(A.next(d));
        double low = fun(min_adjacent(A.rev(d)), 
                         min_adjacent(A.rev(A.next(d))));
        if (X.r < low)
            Y.l = max(Y.l, low);
        if (Y.r < low)
            X.l = max(X.l, low);
        
        if (X.empty() or Y.empty())
            return false;
    }
    
    return true;
}

#endif
