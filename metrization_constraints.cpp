#ifndef METRIZATION_CONSTRAINTS
#define METRIZATION_CONSTRAINTS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "interval_arithmetic.cpp"
#include "quaternions.cpp"
#include "derived_once_nets.cpp"


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

struct PathStep {
    char type;
    double rotation;
    Dart d;
    
    PathStep(double _rotation): type('R'), rotation(_rotation) {}
    PathStep(Dart _d): type('E'), d(_d) {}
    
    double cost(const vector<Interval>& lengths) {
        return 1e-10 + ((type == 'R') ? 0.0 : lengths[d.id].len());
    }
};

// Given a sequence of operations (i.e., a path), it returns the 
// quaternion representing the isometry that maps (e_1, e_2) into the 
// result of the operations starting from (e_1, e_2).
//
// The operations can be:
//  - {'R', \theta}: Apply \rot_\theta.
//  - {'E', l}: Apply \exp_l.
Quaternion<Interval> isometryFromPath(const vector<PathStep>& path, 
                                      const vector<Interval>& lengths) {
    double sum_theta = 0;
    Quaternion<Interval> res(1);
    for (auto op: path) {
        if (op.type == 'R')
            sum_theta += op.rotation;
        if (op.type == 'E')
            res *= basicIsometry<Interval>(sum_theta, lengths[op.d.id]);
    }
    res *= basicRotation<Interval>(sum_theta);
    
    return res;
}

vector<vector<PathStep>> pathToConcaveAngles(
        const DerivedOnceNet& A, const vector<Interval>& lengths) {
    // Step 0: Find the vertices with concave angles.
    vector<int> concave_angles;
    for (Dart d: A.darts())
        if (A.isConcaveAngle(d))
            concave_angles.push_back(d.to);
            
    if (concave_angles.empty())
        return {};
        
    //~ if (ssize(concave_angles) <= 2)
        //~ return {}; // We save some time.
    
    
    // Step 1: Find the best "center".
    //         TODO: Now we just select arbitrarily one of the 
    //               concave angles.
    int v = concave_angles[0];
        
    // Step 2: Run Dijkstra from the chosen dart to the concave angles.
    Dart root;
    for (Dart d: A.darts())
        if (d.from == v)
            root = d;
    
    using pdd = pair<double, Dart>;
    vector<double> distance(2*A.E(), 1e9);
    vector<char> father(2*A.E(), '-'); // n = next, r = rev
    priority_queue<pdd, vector<pdd>, greater<pdd>> q;
    
    distance[root.id] = 0;
    q.push(pdd{0.0, root});
    while (!q.empty()) {
        auto pp = q.top();
        q.pop();
        Dart d = pp.second;
        double dist = pp.first;
        if (dist > distance[d.id])
            continue;
        
        vector<pdd> adj = {{dist + 1e-10, A.next(d)}, 
                           {dist + lengths[d.id].len() + 1e-10, A.rev(d)}};
        for (auto [dist2, d2]: adj) {
            if (dist2 >= distance[d2.id]) 
                continue;
            distance[d2.id] = dist2;
            father[d2.id] = (d.from == d2.from) ? 'n' : 'r';
            q.push({dist2, d2});
        }
    }
    assert(father[root.id] == '-');
    
    // Step 3: Reconstruct the paths from the chosen dart to the edges
    //         going into the pin.
    vector<vector<PathStep>> res;
    for (int goal: concave_angles) {
        double cost = 1e9;
        Dart last_dart;
        for (Dart d: A.darts()) {
            if (d.from == goal and cost > distance[d.id]) {
                cost = distance[d.id];
                last_dart = d;
            }
        }
        Dart curr = last_dart;
        vector<PathStep> path;
        
        // To compress consecutive rotations in a single one.
        auto append_to_path = [&](PathStep op) {
            if (!path.empty() and path.back().type == 'R' and op.type == 'R')
                path.back().rotation += op.rotation;
            else
                path.push_back(op);
        };
        
        // Rotation to go in the direction of the pin.
        double angle_to_pin = A.isConcaveAngle(A.rev(curr)) ? 2*M_PI/3 : M_PI/3;
        append_to_path(angle_to_pin);
        // We go from last_dart to root.
        while (father[curr.id] != '-') {
            if (father[curr.id] == 'n') { // Rotation by 120 or 240.
                double angle = A.isConcaveAngle(A.rev(A.prev(curr))) ? -2*M_PI/3 : 2*M_PI/3;
                append_to_path(angle);
                curr = A.prev(curr);
            } else if (father[curr.id] == 'r') { // Exponential + Rotation by 180.
                append_to_path(M_PI);
                curr = A.rev(curr);
                append_to_path(curr);
            } else assert(0);
        }
        reverse(path.begin(), path.end());
        res.push_back(path);
    }
    
    return res;
}

template <typename T>
T det(const array<T, 3>& A, const array<T, 3>& B, const array<T, 3>& C) {
    return A[0]*(B[1]*C[2] - B[2]*C[1])
         + A[1]*(B[2]*C[0] - B[0]*C[2])
         + A[2]*(B[0]*C[1] - B[1]*C[0]);
}

// Determines, approximately, whether there exists a plane that
// intersects all the given cubes.
// If it returns false, than it cannot. If it returns yes, it could
// be a false positive. 
//
// Implementing a precise check (without false positives) would most
// likely require some linear programming. It would be demanding to
// implement and possibly very slow.
bool areContainedInPlane(vector<array<Interval, 3>> cubes) {
    int q = ssize(cubes);
    for (int i = 0; i < q; i++) 
        for (int j = 0; j < q; j++)
            for (int k = j+1; k < q; k++) {
                if (i == j or i == k)
                    continue;
                if (!det(cubes[i], cubes[j], cubes[k]).containsZero())
                    return false;
            }
    return true;
}

template<typename T>
array<T, 3> applyIsometry(const Quaternion<T>& S, 
                          const array<T, 3>& vec) {
    Quaternion<T> vec2(0, vec[0], vec[1], vec[2]);
    Quaternion<T> res = S*vec2*conj(S);
    return {res[1], res[2], res[3]};
}

bool canPinsBeAntipodal(const DerivedOnceNet& A, 
                        const vector<Interval>& lengths) {
    vector<vector<PathStep>> paths = pathToConcaveAngles(A, lengths);
    vector<Quaternion<Interval>> isometries;
    for (auto path: paths)
        isometries.push_back(isometryFromPath(path, lengths));
    
    vector<array<Interval, 3>> equators;
    for (auto q: isometries)
        equators.push_back(applyIsometry(q, {0, 0, 1}));
    
    return areContainedInPlane(equators);
}


bool checkPinAngles(const DerivedOnceNet& A,
                    const vector<Interval>& lengths) {
    for (const auto& f: A.faces()) {
        if (A.faceType(f).second != 2) // Two pins
            continue;
        
        vector<int> concave_angles;
        for (int i = 0; i < ssize(f); i++) {
            if (A.isConcaveAngle(f[i]))
                concave_angles.push_back(i);
        }
        assert(ssize(concave_angles) == 2);
        vector<PathStep> path;
        path.emplace_back(2*M_PI/3);
        for (int i = concave_angles[0]+1; i <= concave_angles[1]; i++) {
            path.emplace_back(f[i]);
            path.emplace_back(-M_PI/3);
        }
        // T transforms the (directed) edge to one pin to the directed 
        // edge to the other pin.
        Quaternion<Interval> T = isometryFromPath(path, lengths);
        array<Interval, 3> pin_edge = applyIsometry(T, {0, 0, 1});
        
        // scal_prod <= -1/2 is equivalent to the condition on
        // the angle at the pin.
        Interval scal_prod = pin_edge[2];
        //~ cout << scal_prod << endl;
        if (scal_prod.l > -0.5 + 1e-8)
            return false;
    }
    return true;
}


bool refineMetrizationConvexContinent(const DerivedOnceNet& A, 
                                      vector<Interval>& lengths) {
    while (totalLen(lengths) >= 1e-8) {
        double len_before = totalLen(lengths);
        for (const vector<Dart>& f: A.faces()) {
            if (!A.isConvexFace(f))
                continue;
            
            vector<Interval*> sides;
            for (const Dart& d: f)
                sides.push_back(&lengths[d.id]);
            
            assert(3 <= ssize(f) and ssize(f) <= 5);    
            
            if (ssize(f) == 3)
                enforceTriangleConstraint(sides);
            if (ssize(f) == 4)
                enforceQuadrilateralConstraint(sides);    
            if (ssize(f) == 5) 
                enforcePentagonConstraint(sides);
            
            for (const Dart& d: f) 
                lengths[A.rev(d).id].intersect(lengths[d.id]);
        }
        if (totalLen(lengths) > 0.99 * len_before)
            break;
    }
    
    return totalLen(lengths) >= 0;
}

bool refineMetrizationNonConvexFaces(const DerivedOnceNet& A, 
                                     vector<Interval>& lengths) {
    for (int _ = 0; _ < 2; _++) { // We run it twice without a reason.
        for (Dart d: A.darts()) {
            if (A.isConvexFace(d))
                continue;
            int q = ssize(A.face(d));
            vector<PathStep> path;
            for (int i = 0; i < q; i++) {
                Dart e = A.phi(d, i);
                if (i > 0)
                    path.emplace_back(e);
                path.emplace_back(A.phiAngle(e));
            }
            Quaternion<Interval> T = isometryFromPath(path, lengths);
            // T is, up to sign, (cos(l/2), 0, 0, -sin(l/2)), where
            // l is the length of d.
            if (!T[1].containsZero() or !T[2].containsZero()) {
                lengths[d.id] = Interval(0, -1);
                return false;
            }
            Interval I = 2*intersectWithCircle(T[0], -T[3]);
            lengths[d.id].intersect(I);
            lengths[A.rev(d).id].intersect(I);
            if (lengths[d.id].empty())
                return false;
        }
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
double farthestIntersection(Quaternion<Interval> T) {
    // sign1 determines on which side of a is Te_1.
    Interval sign1 = applyIsometry(T, {1, 0, 0})[2];
    // sign2 determines on which side of b is e_1.
    Interval sign2 = applyIsometry(T, {0, 0, 1})[0];
    
    // The first intersection happens in less than pi on both geodesics
    // if and only if sign1 and sign2 have different signs.
    if (sign1.containsZero() or sign2.containsZero())
        return 10;
    if (sign1.l > 0 and sign2.l > 0)
        return 10;
    if (sign1.l < 0 and sign2.l < 0)
        return 10;
    
    // Te_3^e_3 is, up to scaling, the intersection between a and b.
    // Notice that Te_3^e_3 = (X, Y, 0).
    Interval X = applyIsometry(T, {0, 0, 1})[1];
    Interval Y = -applyIsometry(T, {0, 0, 1})[0];
    
    if (Y.containsZero())
        return 10;
        
    // We identify the *first* intersection.
    if (Y.l < 0) X = -X, Y = -Y;
    
    double theta = max(atan2(Y.l, X.l), atan2(Y.r, X.l));
    return theta;
}

bool checkNonIntersection(Quaternion<Interval> T, Interval& I, Interval& J) {
    double lenI = farthestIntersection(T);
    double lenJ = farthestIntersection(conj(T));
    
    //~ if (lenI <= 5 and lenJ <= 5)
        //~ cout << lenI << " " << lenJ << endl;
    
    if (I.l >= lenI)
        J.intersect(Interval(0, lenJ));
    if (J.l >= lenJ)
        I.intersect(Interval(0, lenI));
        
    return !I.empty() and !J.empty();
}

// Checks that non-convex faces have are not self-intersecting.
bool checkNonIntersection(const DerivedOnceNet& A, vector<Interval>& lengths) {
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
                        isometryFromPath(path, lengths),
                        lengths[d.id],
                        lengths[e.id]);
                lengths[A.rev(d).id].intersect(lengths[d.id]);
                lengths[A.rev(e).id].intersect(lengths[e.id]);
                if (!ans)
                    return false;
            }
            path.emplace_back(A.phi(d, i+1));
        }
    }
    return true;
}

#endif
