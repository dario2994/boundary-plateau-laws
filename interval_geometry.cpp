#ifndef INTERVAL_GEOMETRY
#define INTERVAL_GEOMETRY

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "interval_arithmetic.cpp"
#include "quaternions.cpp"

double deg2rad(int deg) {
    return deg * pi.mid() / 180;
}


// Element of R^3.
typedef array<double, 3> Point;

// Parallelepiped in R^3 with axis-aligned sides, i.e., the product of three intervals.
typedef array<Interval, 3> ApproximatePoint;

ostream& operator<<(ostream& out, Point p) {
    out << p[0] << " " << p[1] << " " << p[2];
    return out;
}

// Class representing a ball in SO(3) endowed with the operator norm.
// The ball is described by its center, an isometry, and its radius (a positive
// real value).
//
// The center is represented by a quaternion q.
// A (nonzero) quaternion q represents the isometry v -> qvq^{-1}, through the
// identification R^3 \cong <i, j, k> \subseteq Quaternions.
//
// We enforce that q has norm equal to 1 (up to precision error), as this
// allows us to use q^{-1} = \bar q.
class ApproximateIsometry {
  private:
    Quaternion<double> q;
    double uncertainty;

    // Make interval with uncertainty
    Interval miwu(double mid) const {
        return Interval(mid - uncertainty, mid + uncertainty);
    }
  public:
    ApproximateIsometry(const Quaternion<double>& _q, double _uncertainty) : q(_q), uncertainty(_uncertainty) {
        double n2 = norm2(q);
        assert(1-1e-10 <= n2 and n2 <= 1 + 1e-10);
        assert(0 < uncertainty);
    }

    // Let A be an approximate isometry.
    // It is guaranteed that \{a*v: a\in A\} \subseteq A(v).
    ApproximatePoint operator()(const Point& v) const {
        Quaternion<double> v4(0, v[0], v[1], v[2]);
        Quaternion<double> res = q * v4 * conj(q);
        return {miwu(res[1]), miwu(res[2]), miwu(res[3])};
    }

    ApproximateIsometry& operator *=(const ApproximateIsometry& S) {
        q *= S.q;
        uncertainty += S.uncertainty;
        return *this;
    }
    // Let A, B be two approximate isometries.
    // It is guaranteed that \{a*b: a\in A, b\in B\} \subseteq A*B.
    friend ApproximateIsometry operator *(const ApproximateIsometry&, const ApproximateIsometry&);
    // It is guaranteed that \{a^{-1}: a\in A\} \subseteq inverse(A).
    friend ApproximateIsometry inverse(const ApproximateIsometry&);

    double getUncertainty() { return uncertainty; }
};

// (S, T) -> S*T is 1-Lip in both arguments (wrt the operator norm).
ApproximateIsometry operator *(const ApproximateIsometry& S, const ApproximateIsometry& T) {
    return {S.q * T.q, S.uncertainty + T.uncertainty + 1e-10};
}

// The inverse is an isometry of SO(3) (wrt the operator norm).
ApproximateIsometry inverse(const ApproximateIsometry& S) {
    return {conj(S.q), S.uncertainty + 1e-10};
}

// Let v = (e_1, e_2).
// This function returns (a superset of) R(v, theta).
ApproximateIsometry basicRotation(int theta) {
    double theta2 = 0.5 * deg2rad(theta);
    return ApproximateIsometry({cos(theta2), sin(theta2), 0, 0}, 1e-10);
}

// Let v = (e_1, e_2).
// This function returns (a superset of) E(\rot_\theta v, len)
ApproximateIsometry basicTranslation(int theta_deg, Interval len) {
    double theta_rad = deg2rad(theta_deg);
    Interval len2 = 0.5 * len;
    return ApproximateIsometry({cos(len2.mid()), 0,
            -sin(theta_rad) * sin(len2.mid()), cos(theta_rad) * sin(len2.mid())}, len.len() / 2 + 1e-10);
}


// Return an ApproximateIsometry with very small uncertainty.
ApproximateIsometry randomIsometry() {
    static default_random_engine engine(123);
    uniform_real_distribution<double> dist(-1.0, 1.0);

    double a = dist(engine), b = dist(engine), c = dist(engine), d = dist(engine);
    double z = sqrt(a*a + b*b + c*c + d*d);
    if (z < 0.1)
        return randomIsometry();
    a /= z, b /= z, c /= z, d /= z;
    return ApproximateIsometry({a, b, c, d}, 1e-10);
}


// (superset of the) Intersection of the unit circle with a rectangle.
// It returns a superset of
// \{l \in domain: (cos l, sin l) \in XxY \}
//
// Assumption: domain \subseteq [0, 2pi].
Interval intersectWithCircle(Interval X, Interval Y, const Interval& domain) {
    X.intersect({-1, 1});
    Y.intersect({-1, 1});
    if (X.empty() or Y.empty())
        return EMPTY_INTERVAL;
    Interval I = acos(X);
    Interval J = asin(Y);
    vector<Interval> lI = {I, 2*pi - I};
    vector<Interval> lJ = {pi - J};
    if (J.containsZero()) {
        lJ.emplace_back(0, J.r);
        lJ.emplace_back((2*pi + J).l, (2*pi).r);
    } else if (J.l > 0) {
        lJ.push_back(J);
    } else if (J.r < 0) {
        lJ.push_back(2*pi + J);
    }
    vector<Interval> lK;
    for (const Interval& i: lI) for (const Interval& j: lJ) {
        Interval k = intersection(intersection(i, j), domain);
        if (!k.empty())
            lK.push_back(k);
    }
    return hull(lK);
}

template<typename T>
array<T, 3> operator*(const array<T, 3>& X, T lambda) {
    return {X[0] * lambda, X[1] * lambda, X[2] * lambda};
}

template<typename T>
array<T, 3> operator/(const array<T, 3>& X, T lambda) {
    return X * (T(1)/lambda);
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
array<T, 3> operator+(const array<T, 3>& X, const array<T, 3>& Y) {
    return {X[0]+Y[0], X[1]+Y[1], X[2]+Y[2]};
}

template<typename T>
array<T, 3> operator-(const array<T, 3>& X) {
    return {-X[0], -X[1], -X[2]};
}

template<typename T>
array<T, 3> operator^(const array<T, 3>& X, const array<T, 3>& Y) {
    return {X[1]*Y[2]-X[2]*Y[1], X[2]*Y[0]-X[0]*Y[2], X[0]*Y[1]-X[1]*Y[0]};
}

template <typename T>
T det(const array<T, 3>& A, const array<T, 3>& B, const array<T, 3>& C) {
    return A * (B ^ C);
}

template <typename T>
T norm(const array<T, 3>& X) {
    return sqrt(square(X[0]) + square(X[1]) + square(X[2]));
}

double norm(const Point& X) {
    return sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
}

Point make_deterministic(ApproximatePoint p) {
    return {p[0].mid(), p[1].mid(), p[2].mid()};
}


// T is an approximate isometry.
// Let a be the oriented half-equator starting from (e_1, e_2).
// Let b be the oriented half-equator starting from (T'e_1, T'e_2), where
// T' is an isometry in T.
//
// The return values of this function are:
//  0: a and b are guaranteed not to intersect.
//  1 or -1: a and b are guaranteed to intersect. The sign depends on
//           the sign of the intersection.
//  INDETERMINATE: This function does not guarantee anything.
int halfEquatorsIntersection(const ApproximateIsometry& T) {
    // sign1 determines on which side of a is Te_1.
    int sign1 = sign(T({1, 0, 0})[2]);
    // sign2 determines on which side of b is e_1.
    int sign2 = sign(T({0, 0, 1})[0]);

    if (sign1 == 0 or sign2 == 0)
        return INDETERMINATE;

    if (sign1 * sign2 == 1)
        return 0; // The two half-equators do not intersect.

    // Either 1 or -1 depending on the orientation of the intersection.
    return sign1;
}

// T is an approximate isometry.
// Let a be the oriented geodesic starting from (e_1, e_2).
// Let b be the equator going through (T'e_1, T'e_2), where T' is an isometry in T.
//
// This function returns a nonempty interval I \subseteq [0, pi].
// It is guaranteed that the first intersection of a with b happens at
// exp((e_1, e_2), t) with t\in I.
Interval firstIntersectionWithEquator(const ApproximateIsometry& T) {
    auto Te3 = T({0, 0, 1});

    // Te_3^e_3 is, up to scaling, the intersection between a and b.
    // Notice that Te_3^e_3 = (X, Y, 0).
    Interval X = Te3[1];
    Interval Y = -Te3[0];

    if (Y.containsZero())
        return Interval(0, pi.r);

    // We identify the *first* intersection.
    if (Y.l < 0) X = -X, Y = -Y;

    double min_theta = min(atan2(Y.l, X.r), atan2(Y.r, X.r));
    double max_theta = max(atan2(Y.l, X.l), atan2(Y.r, X.l));
    return Interval(min_theta, max_theta);
}

// T is an approximate isometry.
// I, J are two intervals contained in [0, pi.r].
// Let a be a geodesic segment with length in I starting from (e_1, e_2).
// Let b be a geodesic segment with length in J starting from (T'e_1, T'e_2),
// where T' is an isometry in T.
//
// The return values of this function are:
//  0: a and b are guaranteed not to intersect.
//  1 or -1: a and b are guaranteed to intersect. The sign depends on
//           the sign of the intersection.
//  INDETERMINATE: This function does not guarantee anything.
int internalIntersectionSign(const ApproximateIsometry& T,
                             const Interval& I, const Interval& J) {
    assert(0 <= I.l and I.r <= pi.r);
    assert(0 <= J.l and J.r <= pi.r);

    int x = halfEquatorsIntersection(T);
    if (x == 0)
        return 0;
    if (x == INDETERMINATE)
        return INDETERMINATE;

    Interval I2 = firstIntersectionWithEquator(T);
    Interval J2 = firstIntersectionWithEquator(inverse(T));

    if (I2.l > I.r or J2.l > J.r)
        return 0;

    if (I2.r < I.l and J2.r < J.l)
        return x;

    return INDETERMINATE;
}

// T is an approximate isometry.
// I, J are two intervals contained in [0, 2pi].
// Let a be a geodesic segment with length in I starting from (e_1, e_2).
// Let b be a geodesic segment with length in J starting from (T'e_1, T'e_2),
// where T' is an isometry in T.
//
// The return values of this function are:
//  0: a and b are guaranteed not to intersect or to have exactly two intersections (in the latter
//     case, the 0 represents the signed sum of the two intersections points).
//  1 or -1: a and b are guaranteed to intersect. The sign depends on
//           the sign of the intersection.
//  INDETERMINATE: This function does not guarantee anything.
int intersectionSign(const ApproximateIsometry& T,
                     const Interval& I, const Interval& J) {
    int ans = 0;
    vector<pair<Interval, Interval>> halves = {{I, J}};
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
        Interval I2 = intersection(Interval(0, pi.r), i == 0 ? I : I-pi);
        Interval J2 = intersection(Interval(0, pi.r), j == 0 ? J : J-pi);
        ApproximateIsometry S = T;
        ApproximateIsometry k({0, 0, 0, 1}, 1e-10);
        if (i == 1)
            S = k * S;
        if (j == 1)
            S = S * k;
        int curr = internalIntersectionSign(S, I2, J2);
        if (ans != INDETERMINATE and curr != INDETERMINATE)
            ans += curr;
        else
            ans = INDETERMINATE;
    }
    if (ans == INDETERMINATE)
        return INDETERMINATE;
    return ans;
}

// T is an approximate isometry.
// I, J are two intervals contained in [0, 2pi].
// Let a be a geodesic segment with length in I starting from (e_1, e_2).
// Let b be a geodesic segment with length in J starting from (T'e_1, T'e_2),
// where T' is an isometry in T.
//
// The return values of this function are:
//  0: a and b are guaranteed not to intersect.
//  1: a and b are guaranteed to intersect.
//  INDETERMINATE: This function does not guarantee anything.
int areIntersecting(const ApproximateIsometry& T,
                    Interval& I, Interval& J) {
    if (I.r >= pi.l or J.r >= pi.l)
        return INDETERMINATE; // We do not handle these cases.
    int x = halfEquatorsIntersection(T);
    if (x == INDETERMINATE)
        return INDETERMINATE;
    if (x == 0)
        return 0;

    double lenI = firstIntersectionWithEquator(T).r;
    if (lenI >= I.r) // We put this check here to speed up the execution.
        return 0;
    double lenJ = firstIntersectionWithEquator(inverse(T)).r;
    if (lenJ >= J.r)
        return 0;

    if (lenI <= I.l and lenJ <= J.l)
        return 1;

    return INDETERMINATE;
}



#endif
