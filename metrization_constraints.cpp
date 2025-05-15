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


// Enforces the constraint that the three numbers in sides correspond
// to the three sides of a spherical triangle with angles 120.
void enforceTriangleConstraint(vector<Interval*> sides) {
    assert(ssize(sides) == 3);
    
    static const double L = acos(-1.0/3);
    static const Interval I(L-1e-9, L+1e-9); 
    for (int i = 0; i < 3; i++) {
        sides[i]->intersect(I);
    }
}

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

// Enforces the constraint that the four numbers in sides correspond
// to the four sides of a spherical triangle with angles 120.
//
// Assumption: all sides are subsets of [0, pi].
void enforceQuadrilateralConstraint(vector<Interval*> sides) {
    assert(ssize(sides) == 4);
    
    sides[0]->intersect(*sides[2]);
    sides[2]->intersect(*sides[0]);
    
    sides[1]->intersect(*sides[3]);
    sides[3]->intersect(*sides[1]);
    
    Interval even = quadrilateralConstraint(*sides[1]);
    Interval odd = quadrilateralConstraint(*sides[0]);
    
    sides[0]->intersect(even);
    sides[2]->intersect(even);
    sides[1]->intersect(odd);
    sides[3]->intersect(odd);
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
    
    // Using acos on [low, high] is not necessarily correct here, hence
    // we do the computation explicitly.
    
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


// Enforces the constraint that the five numbers in sides correspond
// to the five sides of a spherical triangle with angles 120.
//
// Assumption: all sides are subsets of [0, \ell_T].
void enforcePentagonConstraint(vector<Interval*> sides) {
    assert(ssize(sides) == 5);
    
    auto s = [&](int id) { return sides[id % 5]; };
        
    auto compute_len = [&]() {
        double res = 0;
        for (int i = 0; i < 5; i++) {
            if (s(i)->empty())
                return -1.0;
            res += sides[i]->len();
        }
        return res;
    };
    
    double curr_len = compute_len();
    while (curr_len >= 1e-8) {
        for (int i = 0; i < 5; i++) {
            s(i+3)->intersect(pentagon1Constraint(*s(i), *s(i+1)));
            s(i)->intersect(pentagon2Constraint(*s(i+1), *s(i+3)));
            s(i)->intersect(pentagon2Constraint(*s(i+4), *s(i+2)));
        }
        double new_len = compute_len();
        if (new_len > 0.99*curr_len)
            break;
        curr_len = new_len;
    }
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

template <typename T>
T det(const array<T, 3>& A, const array<T, 3>& B, const array<T, 3>& C) {
    return A[0]*(B[1]*C[2] - B[2]*C[1])
         + A[1]*(B[2]*C[0] - B[0]*C[2])
         + A[2]*(B[0]*C[1] - B[1]*C[0]);
}

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
                    auto isom = net.getPathToPin(root, d);
                    equators.push_back(applyIsometry(isom, {0, 0, 1}));
                }                
                if (!det(equators[0], equators[1], equators[2]).containsZero())
                    return false;
            }
    
    return true;
}


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
        auto iso0 = net.getPathToPin(root, concave_angles[0]);
        auto iso1 = net.getPathToPin(root, concave_angles[1]);
        
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


bool refineMetrizationConvexContinent(MetrizedNet& net) {
    StopWatch _("refineMetrizationConvexContinent");
    
    const DerivedOnceNet& A = net.A;
    
    while (totalLen(net.lengths()) >= 1e-8) {
        double len_before = totalLen(net.lengths());
        for (const vector<Dart>& f: A.faces()) {
            if (!A.isConvexFace(f))
                continue;
            
            vector<Interval*> sides;
            for (const Dart& d: f)
                sides.push_back(&net.length(d));
            
            assert(3 <= ssize(f) and ssize(f) <= 5);    
            
            if (ssize(f) == 3)
                enforceTriangleConstraint(sides);
            if (ssize(f) == 4)
                enforceQuadrilateralConstraint(sides);    
            if (ssize(f) == 5) 
                enforcePentagonConstraint(sides);
            
            for (const Dart& d: f) 
                net.length(A.rev(d)).intersect(net.length(d));
        }
        if (totalLen(net.lengths()) > 0.99 * len_before)
            break;
    }
    
    return totalLen(net.lengths()) >= 0;
}

bool refineMetrizationNonConvexFaces(MetrizedNet& net) {
    StopWatch _("refineMetrizationNonConvexFaces");
    
    const DerivedOnceNet& A = net.A;
    
    for (Dart d: A.darts()) {
        if (A.isConvexFace(d))
            continue;
        int q = ssize(A.face(d));
        
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
            path.emplace_back(A.phiAngle(e));
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
        net.length(A.rev(d)).intersect(I);
        if (net.length(d).empty())
            return false;
    }
    return true;
}

// T is an approximate isometry.
// Consider the directed geodesic a starting from (e_1, e_2) and the 
// directed geodesic b starting from (Te_1, Te_2).
// 
// This function checks whether the first intersection of a and b
// happens after less then pi distance from both starting points.
// If this is the case, it returns the maximum distance of the 
// intersection point along a. Otherwise it returns a large number.
double farthestIntersection(const Quaternion<Interval>& T) {
    //~ StopWatch _("farthestIntersection");
    
    auto Te3 = applyIsometry(T, {0, 0, 1}); // For efficiency.
    
    // sign1 determines on which side of a is Te_1.
    int sign1 = sign(applyIsometry(T, {1, 0, 0})[2]);
    // sign2 determines on which side of b is e_1.
    int sign2 = sign(Te3[0]);
    
    // The first intersection happens in less than pi on both geodesics
    // if and only if sign1 and sign2 have different signs.
    if (sign1 * sign2 != -1)
        return 10;
    
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
    
    double lenI = farthestIntersection(T);
    double lenJ = farthestIntersection(conj(T));
    
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
        vector<PathStep> path;
        path.emplace_back(M_PI);
        int q = ssize(A.face(d));
        for (int i = 0; i < q-2; i++) {
            path.emplace_back(A.phiAngle(A.phi(d, i)));
            if (i != 0) {
                Dart e = A.phi(d, i+1);
                bool ans = checkNonIntersection(
                        net.isometryFromPath(path),
                        net.length(d),
                        net.length(e));
                net.length(A.rev(d)).intersect(net.length(d));
                net.length(A.rev(e)).intersect(net.length(e));
                if (!ans)
                    return false;
            }
            path.emplace_back(A.phi(d, i+1));
        }
    }
    return true;
}


bool checkPinOutsideConvexFaces(const MetrizedNet& net) {
    StopWatch _("checkPinOutsideConvexFaces");
    
    int cnt = ssize(net.A.concaveAngles());
    for (int i = 0; i < cnt; i++) for (int j = i+1; j < cnt; j++) {
        const Dart& root = net.A.concaveAngles()[i];
        const Dart& other = net.A.concaveAngles()[j];
        const auto& X = applyIsometry(net.getPathToPin(root, root), {0, 0, 1});
        const auto& Y = applyIsometry(net.getPathToPin(root, other), {0, 0, 1});
        for (const auto& f: net.A.faces()) {
            if (!net.A.isConvexFace(f))
                continue;
            // Check that none of the two intersections of X and Y lie
            // in the convex face f.
            
            vector<int> signs;
            for (const Dart& d: f) {
                const auto& Z = applyIsometry(net.getPath(root, d), {0, 0, 1});
                int s = sign(det(X, Y, Z));
                signs.push_back(s);
                
                if (s == 0 or s != signs[0])
                    break; // For performance.
            }
            bool is_inside = true;
            for (int s: signs) {
                is_inside &= s != 0 and s == signs[0];
            }
            if (is_inside)
                return false;
        }
    }
    return true;
}

#endif
