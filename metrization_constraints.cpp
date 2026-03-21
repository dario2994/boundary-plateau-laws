#ifndef METRIZATION_CONSTRAINTS
#define METRIZATION_CONSTRAINTS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "marked_combinatorial_maps.cpp"
#include "utils.cpp"
#include "interval_arithmetic.cpp"
#include "interval_geometry.cpp"
#include "metrized_nets.cpp"



// Utility function to enforce the constraint on the sides of a
// quadrilateral.
//
// I is strictly inside [0, \pi].
// One side of a quadrilateral belongs to I.
// The other side belongs to the interval returned by this function.
Interval quadrilateralConstraint(Interval I) {
    if (I.empty())
        return EMPTY_INTERVAL;
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
        return EMPTY_INTERVAL;

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
        return EMPTY_INTERVAL;

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
        return EMPTY_INTERVAL;

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

    const MarkedCombinatorialMap& A = net.A;

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


// Checks that the equators corresponding to any three edges going to
// a pin are concurrent.
bool checkPinsAntipodality(const MetrizedNet& net) {
    StopWatch _("checkPinsAntipodality");

    int q = ssize(net.A.concaveAngles());
    for (int i = 0; i < q; i++)
        for (int j = i+1; j < q; j++)
            for (int k = j+1; k < q; k++) {
                Dart root = net.A.concaveAngles()[i];

                vector<ApproximatePoint> equators;
                for (int z: {i, j, k}) {
                    Dart d = net.A.concaveAngles()[z];
                    // Using getIsometryFromPinToPin instead of getIsometryToPin
                    // does not change the result, but dramatically changes
                    // the precision of the computation.
                    // Our choice guarantees that equators[0] = (0, 0, 1).
                    // This reduces the computation of a 3x3 determinant to a
                    // 2x2 determinant.
                    // In practice, this is necessary to "solve" net 28 in a
                    // reasonable amount of time.
                    auto isom = net.getIsometryFromPinToPin(root, d);
                    equators.push_back(isom({0, 0, 1}));
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

    const MarkedCombinatorialMap& A = net.A;

    for (const auto& f: A.faces()) {
        if (A.faceType(f).second != 2) // The face does not have two concave angles.
            continue;

        vector<Dart> concave_angles;
        for (Dart d:  f) {
            if (A.isConcaveAngle(d))
                concave_angles.push_back(d);
        }
        assert(ssize(concave_angles) == 2);

        ApproximateIsometry iso = net.getIsometryFromPinToPin(concave_angles[0], concave_angles[1]);
        Interval scal_prod = iso({0, 0, 1})[2];

        // scal_prod <= -1/2 is equivalent to the condition on the angle at the pin.
        if (scal_prod.l > -0.5 + 1e-8)
            return false;
    }
    return true;
}

// Assume that you know the lengths of k-1 sides of a non-convex face
// with k sides. This function determines (robustly) whether they
// can form the non-convex face (i.e., the angles are correct) and
// enforces the right length of the remaining side.
bool refineMetrizationNonConvexFaces(MetrizedNet& net) {
    StopWatch _("refineMetrizationNonConvexFaces");

    const MarkedCombinatorialMap& A = net.A;

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
        ApproximateIsometry T = net.isometryFromPath(path);
        // Let l be the length of d.
        // T(e_1) = (cos l, -sin l, 0)
        // T(e_2) = (sin l,  cos l, 0)
        // T(e_3) = (0,          0, 0)
        ApproximatePoint T0 = T({1, 0, 0});
        ApproximatePoint T1 = T({0, 1, 0});
        ApproximatePoint T2 = T({0, 0, 1});
        if (!T0[2].containsZero() or !T1[2].containsZero()
              or !T2[0].containsZero() or !T2[1].containsZero()
              or !(T2[2]-1).containsZero()) {
            net.length(d) = EMPTY_INTERVAL;
            return false;
        }
        Interval constraint1 = intersectWithCircle(T0[0], -T0[1], net.length(d));
        Interval constraint2 = intersectWithCircle(T1[1],  T1[0], net.length(d));
        Interval total_constraint = intersection(constraint1, constraint2);
        net.length(d).intersect(total_constraint);
        if (net.length(d).empty())
            return false;
    }
    return true;
}

// Checks that non-convex faces are not self-intersecting.
bool checkNonIntersection(MetrizedNet& net) {
    StopWatch _("checkNonIntersection");

    const MarkedCombinatorialMap& A = net.A;
    for (Dart d: A.darts()) {
        if (A.isConvexFace(d))
            continue;

        // We don't use net.isometryFromPath to improve the performance.
        double sum_theta = 0;
        sum_theta += M_PI;
        const vector<Dart>& f = A.face(d);
        ApproximateIsometry iso(1, 0);

        for (int i = 0; i < ssize(f)-2; i++) {
            sum_theta += A.phiAngle(f[i]);
            if (i != 0) {
                Dart e = f[i+1];
                int ans = areIntersecting(
                        iso * basicRotation(sum_theta),
                        net.length(d),
                        net.length(e));
                if (ans == 1)
                    return false;
            }
            iso *= basicTranslation(sum_theta, net.length(f[i+1]));

            // We stop checking when the uncertainty gets so large
            // that we would not be able to get any reasonable answer
            // from now on.
            if (iso.getUncertainty() > 0.5)
                break;
        }
    }
    return true;
}

// The face must be convex.
// Set the dart root to (e_1, e_2).
// According to the data in net, this function returns:
// - 1 if it can prove that the point P is (strictly) inside the face.
// - 0 if it can prove that the point P is not (strictly) inside the face.
// - INDETERMINATE if it cannot prove anything.
//
// If up_to_symmetry is true, then this returns 1 if it can show that either
// P or -P is inside the face, 0 if it can show that both are not inside,
// and INDETERMINATE if it cannot prove anything.
int isInsideConvexFace(const MetrizedNet& net,
                       ApproximatePoint P,
                       const Dart& root,
                       const vector<Dart>& face,
                       bool up_to_symmetry = false) {
    assert(net.A.isConvexFace(face));

    vector<int> signs;
    for (const Dart& d: face) {
        const auto& Z = net.getIsometry(root, d)({0, 0, 1});
        int s = sign(P * Z);
        signs.push_back(s);

        if (s == INDETERMINATE)
            return INDETERMINATE;

        if (s != signs[0])
            return 0;
    }

    return up_to_symmetry or signs[0] == -1;
}

// Checks that no pins fall strictly inside a convex face.
bool checkPinOutsideConvexFaces(const MetrizedNet& net) {
    StopWatch _("checkPinOutsideConvexFaces");

    int cnt = ssize(net.A.concaveAngles());
    for (int i = 0; i < cnt; i++) for (int j = i+1; j < cnt; j++) {
        const Dart& root = net.A.concaveAngles()[i];
        const Dart& other = net.A.concaveAngles()[j];

        const auto& X = net.getIsometryToPin(root, root)({0, 0, 1});
        const auto& Y = net.getIsometryToPin(root, other)({0, 0, 1});
        const auto& P = X^Y;
        for (const auto& f: net.A.faces()) {
            if (!net.A.isConvexFace(f))
                continue;

            // Check that none of the two intersections of X and Y lie
            // in the convex face f.
            if (isInsideConvexFace(net, P, root, f, true) == 1)
                return false;
        }
    }
    return true;
}

// For each concave angle, it returns an upperbound on the length of
// the edge going from that angle into the pin.
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
    const MarkedCombinatorialMap& A = net.A;

    vector<double> len_to_pin;
    for (const Dart& to_pin: A.concaveAngles()) {
        double len = 10;
        for (const Dart& d: A.darts()) {
            const auto& T = net.getIsometryToPin(d, to_pin);

            int x = halfEquatorsIntersection(T);
            if (x == 1 or x == -1) {
                if (firstIntersectionWithEquator(T).r < net.length(d).l)
                    len = min(len, firstIntersectionWithEquator(inverse(T)).r);
            }

            if (!A.isConcaveAngle(d) and !A.isConcaveAngle(A.rev(d)))
                continue;

            if (!A.sameFace(d, to_pin))
                continue;

            if (to_pin.to == d.from or to_pin.to == d.to)
                continue;

            if (T({1, 0, 0})[2].r < 0)
                len = min(len, firstIntersectionWithEquator(inverse(T)).r);
        }
        len_to_pin.push_back(len);
    }

    return len_to_pin;
}

// This is a check useful (almost) only when all the equators identifying the
// pins coincide (as in this case checkPinsAntipodality does not work).
//
// Consider two "edges to pins" going to different pins.
// Let l1, l2 be the upper bounds on their lenghts.
// And let q1, q2 be the edges extended to have length l1 and l2.
// Since p\in q1 and -p\in q2 (where p is a pin), it must hold
// dist(q1, q2) + l1 + l2 >= pi.
//
// Why aren't we directly checking whether q1 and -q2 intersect? Because this
// kind of check is very hard in degenerate cases and this function is interesting
// almost only in degenerate cases (e.g., when q1, q2 belong to the same
// equator).
bool checkPinsAntipodality2(const MetrizedNet& net) {
    StopWatch _("checkPinsAntipodality2");

    auto eucl_to_sphere = [&](double d) { return 2*asin(d/2); };
    const MarkedCombinatorialMap& A = net.A;

    vector<double> len_to_pin = maxDistanceToPin(net);

    for (int i = 0; i < ssize(A.concaveAngles()); i++) {
        for (int j = i+1; j < ssize(A.concaveAngles()); j++) {
            Dart d1 = A.concaveAngles()[i], d2 = A.concaveAngles()[j];
            if (A.sameFace(d1, d2))
                continue;
            double l1 = len_to_pin[i], l2 = len_to_pin[j];

            if (l1 >= M_PI-1e-5 or l2 >= M_PI-1e-5)
                return true;

            ApproximatePoint start1 = {1, 0, 0}, end1 = {cos(l1), sin(l1), 0};
            vector<PathStep> path = net.getPathFromPinToPin(d1, d2);
            ApproximatePoint start2 = net.isometryFromPath(path)({1, 0, 0});
            path.push_back(PathStep('E', l2));
            ApproximatePoint end2 = net.isometryFromPath(path)({1, 0, 0});

            for (auto p1: {start1, end1}) for (auto p2: {start2, end2}) {
                Interval eucl_dist = norm(p1-p2);
                if (eucl_dist.r > 2-0.1)
                    continue;
                double max_sphere_dist = eucl_to_sphere(eucl_dist.r);
                if (l1 + l2 + max_sphere_dist < M_PI - 1e-5)
                    return false;
            }
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

    const MarkedCombinatorialMap& A = net.A;

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
        double low = 0.5 * fun(min_adjacent(A.rev(d)),
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

// For 100 times, we randomly pick a geodesic segment g of length 1 on
// the sphere.
// For each face F, we check that the the (signed) total intersection of
// g with \partial F belongs to {-1, 0, 1}.
// Moreover,
//  - If g(0) and g(1) belong to the same face, then all intersection
//    values shall be 0.
//  - If g(0) and g(1) belong to different faces, then exactly one face
//    shall yield value -1 and exactly one face shall yield value -1.
bool checkFacesAreBoundaries(const MetrizedNet& net) {
    StopWatch _("checkFacesAreBoundaries");

    for (int _ = 0; _ < 100; _++) {
        ApproximateIsometry T = randomIsometry();
        int pos = 0, neg = 0;
        bool ambiguous = false;
        const Dart& root = net.A.dart(0);

        for (const auto& f: net.A.faces()) {
            int cnt = 0;
            for (const Dart& d: f) {
                int x = intersectionSign(T*net.getIsometry(root, d),
                                         Interval(1), net.length(d));
                if (x == INDETERMINATE) {
                    ambiguous = true;
                    break;
                }
                cnt += x;
            }
            if (ambiguous)
                break;
            if (abs(cnt) > 1)
                return false;
            pos += cnt == 1;
            neg += cnt == 1;
        }
        if (ambiguous)
            continue;

        if (pair<int,int>{pos, neg} != pair<int,int>{0, 0} and
            pair<int,int>{pos, neg} != pair<int,int>{1, 1})
            return false;
    }
    return true;
}

#endif
