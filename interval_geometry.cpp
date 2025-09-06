#ifndef INTERVAL_GEOMETRY
#define INTERVAL_GEOMETRY

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "interval_arithmetic.cpp"
#include "quaternions.cpp"


// A quaternion with q with unit norm induces an isometry on R^3 via
// v -> qv\bar q, through the identification R^3 \cong <i, j, k> \subseteq Quaternions.
class ApproximateIsometry {
  private:
    Quaternion<double> q;
    double uncertainty;

    // Make interval with uncertainty
    Interval miwu(double mid) const {
        return Interval(mid - uncertainty, mid + uncertainty);
    }
  public:
    ApproximateIsometry(const Quaternion<double>& _q, double _uncertainty) : q(_q), uncertainty(_uncertainty) {}

    array<Interval, 3> operator()(array<double, 3> vec) const {
        Quaternion<double> vec4(0, vec[0], vec[1], vec[2]);
        Quaternion<double> res = q * vec4 * conj(q);
        return {miwu(res[1]), miwu(res[2]), miwu(res[3])};
    }

    // Quaternion<Interval> getQuaternion() const { return q; }

    ApproximateIsometry& operator *=(const ApproximateIsometry& S) {
        q *= S.q;
        uncertainty += S.uncertainty;
        return *this;
    }
    friend ApproximateIsometry operator *(const ApproximateIsometry&, const ApproximateIsometry&);
    friend ApproximateIsometry inverse(const ApproximateIsometry&);

    double getUncertainty() { return uncertainty; }
};

// Here we are implicitly using that (S, T) -> S*T is 1-Lip in both its
// arguments when SO(3) is endowed with the operator norm.
ApproximateIsometry operator *(const ApproximateIsometry& S, const ApproximateIsometry& T) {
    return {S.q * T.q, S.uncertainty + T.uncertainty};
}

// Here we are implicitly using that the inverse is an isometry of SO(3)
// endowed with the operator norm.
ApproximateIsometry inverse(const ApproximateIsometry& S) {
    return {conj(S.q), S.uncertainty};
}

// Set v = (e_1, e_2).
// This function returns R(v, theta).
ApproximateIsometry basicRotation(double theta) {
    double theta2 = 0.5 * theta;
    return ApproximateIsometry({cos(theta2), sin(theta2), 0, 0}, 0);
}

// Set v = (e_1, e_2).
// This function returns E(\rot_\theta v, len)
ApproximateIsometry basicTranslation(double theta, Interval len) {
    Interval len2 = 0.5 * len;
    return ApproximateIsometry({cos(len2.mid()), 0,
            -sin(theta) * sin(len2.mid()), cos(theta) * sin(len2.mid())}, len.len() / 2);
}


// TODO: This function shall guarantee that it returns an isometry (i.e.,
//       a quaternion with norm 1 is contained in the 4-dim box).
//       Should this simply use basicRotation and basicTranslation?
ApproximateIsometry randomIsometry() {
    static default_random_engine engine(123);
    uniform_real_distribution<double> dist(-1.0, 1.0);

    double a = dist(engine), b = dist(engine), c = dist(engine), d = dist(engine);
    double z = sqrt(a*a + b*b + c*c + d*d);
    if (z < 0.1)
        return randomIsometry();
    a /= z, b /= z, c /= z, d /= z;
    return ApproximateIsometry({a, b, c, d}, 0);
}


// (superset of the) Intersection of the unit circle with a rectangle.
// It returns a superset of
// \{l \in domain: (cos l, sin l) \in XxY \}
// Assumption: domain \subseteq [0, 2pi].
//
// TODO: This shall be tested. It can be robustly tested as follows:
//  1. Generate a two random intervals X, Y in [-10, 10].
//  2. Generate a random 0 <= theta <= 2*pi. If cos(theta) \in X
//     and sin(\theta) \in Y, then call this function and check that theta
//     belongs to the returned interval.
Interval intersectWithCircle(Interval X, Interval Y, const Interval& domain) {
    X.intersect({-1, 1});
    Y.intersect({-1, 1});
    if (X.empty() or Y.empty())
        return EMPTY_INTERVAL;
    Interval I = {acos(X.r), acos(X.l)}; // Subset of [0, pi].
    Interval J = {asin(Y.l), asin(Y.r)}; // Subset of [-pi/2, pi/2].
    vector<Interval> lI = {I, {2*M_PI-I.r, 2*M_PI-I.l}};
    vector<Interval> lJ = {{M_PI-J.r, M_PI-J.l}};
    if (J.containsZero()) {
        lJ.emplace_back(0, J.r);
        lJ.emplace_back(2*M_PI + J.l, 2*M_PI);
    } else if (J.l > 0) {
        lJ.push_back(J);
    } else if (J.r < 0) {
        lJ.emplace_back(2*M_PI + J.l, 2*M_PI + J.r);
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
    return sqrt(pow(X[0], 2) + pow(X[1], 2) + pow(X[2], 2));
}

array<double, 3> make_deterministic(array<Interval, 3> p) {
    return {p[0].mid(), p[1].mid(), p[2].mid()};
}


// T is an approximate isometry.
// Let a be the oriented half-equator starting from (e_1, e_2).
// Let b be the oriented half-equator starting from (Te_1, Te_2).
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
// Let b be the equator going through (Te_1, Te_2).
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
        return Interval(0, M_PI);

    // We identify the *first* intersection.
    if (Y.l < 0) X = -X, Y = -Y;

    double min_theta = min(atan2(Y.l, X.r), atan2(Y.r, X.r));
    double max_theta = max(atan2(Y.l, X.l), atan2(Y.r, X.l));
    return Interval(min_theta, max_theta);
}

// T is an approximate isometry.
// Let a be the geodesic segment with length I starting from (e_1, e_2).
// Let b be the geodesic segment with length J starting from (Te_1, Te_2).
// We assume that I, J \subseteq [0, pi].
//
// a and b may have zero or one intersections. Moreover, since we don't know b exactly,
// we may not be able to tell whether they intersect.
//
// The return values of this function are:
//  0: a and b are guaranteed not to intersect.
//  1 or -1: a and b are guaranteed to intersect. The sign depends on
//           the sign of the intersection.
//  INDETERMINATE: This function does not guarantee anything.
int internalIntersectionSign(const ApproximateIsometry& T,
                             const Interval& I, const Interval& J) {
    assert(0 <= I.l and I.r <= M_PI+1e-10);
    assert(0 <= J.l and J.r <= M_PI+1e-10);

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
// Let a be the geodesic segment with length I starting from (e_1, e_2).
// Let b be the geodesic segment with length J starting from (Te_1, Te_2).
// We assume that I, J \subseteq [0, 2pi].
//
// a and b may have zero or one or two intersections. Moreover, since we don't know them exactly,
// we may not be able to tell how many intersections they have.
//
// The return values of this function are:
//  0: a and b are guaranteed not to intersect or to have exactly two intersections (in the latter
//     case, the 0 represents the signed sum of the two intersections points).
//  1 or -1: a and b are guaranteed to intersect. The sign depends on
//           the sign of the intersection.
//  INDETERMINATE: This function does not guarantee anything.
int intersectionSign(const ApproximateIsometry& T,
                     const Interval& I, const Interval& J) {
    const Interval& pi = Interval(0, M_PI);
    int ans = 0;
    vector<pair<Interval, Interval>> halves = {{I, J}};
    for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) {
        Interval I2 = intersection(pi, i == 0 ? I : I-M_PI);
        Interval J2 = intersection(pi, j == 0 ? J : J-M_PI);
        ApproximateIsometry S = T;
        ApproximateIsometry k({0, 0, 0, 1}, 0);
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


// Consider the two geodesic segments:
// - One starting from (e_1, e_2) with length I.
// - One starting from (Te_1, Te_2) with length J.
// This function returns false if they necessarily intersect.
// (it can change I and J to satisfy the non-intersection requirement)
//
// If the intersection happens after at least pi along one of the two
// geodesic segments, it is not detected.
// The guarantee that "false -> they intersect" is not broken, but it
// returns true even if the segments are clearly intersecting.
//
// TODO: Do we really need to allow this function to change I and J?
//       It seems like it makes the logic more complicated, without
//       necessarily improving the performance.
//       If we change it, it becomes easier (though still somewhat messy) to support the
//       intersections past pi.
bool checkNonIntersection(const ApproximateIsometry& T,
                          Interval& I, Interval& J) {
    int x = halfEquatorsIntersection(T);
    if (x != 1 and x != -1)
        return true;
    double lenI = firstIntersectionWithEquator(T).r;
    if (lenI >= I.r) // for performance
        return true;
    double lenJ = firstIntersectionWithEquator(inverse(T)).r;

    if (I.l >= lenI)
        J.intersect(Interval(0, lenJ));
    if (J.l >= lenJ)
        I.intersect(Interval(0, lenI));

    return !I.empty() and !J.empty();
}

// If it returns true, the given points belong to a strict hemisphere.
// Complexity: O(n^4), where n = ssize(vec).
//
// TODO: This needs to be tested. I don't know how to test it.
// TODO: Also the naming is unsafe, as it seems to imply that if it returns
//       false then they cannot belong to the same Hemisphere.
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
                same_sign &= s != INDETERMINATE and s == signs[0];
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



#endif
