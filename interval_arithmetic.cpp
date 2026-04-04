#ifndef INTERVAL_ARITHMETIC
#define INTERVAL_ARITHMETIC

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

static const int INDETERMINATE = -9832521;

// Struct maintaining a closed interval [l, r].
// In many of our applications, the interval is a subset of [0, pi].
struct Interval {
    double l, r;

    Interval(): l(-1e9), r(1e9) {}
    Interval(double _l, double _r): l(_l), r(_r) {}
    Interval(double x): l(x), r(x) {}

    bool empty() const { return r < l; }

    void intersect(Interval I) {
        l = max(l, I.l);
        r = min(r, I.r);
    }

    double len() const { return empty() ? 0.0 : (r-l); }
    double mid() const {
        assert(!empty());
        return (l + r) / 2;
    }
    bool contains(double x) const { return l <= x and x <= r; }
    bool containsZero() const { return contains(0); }
};

const Interval EMPTY_INTERVAL(0, -1);
const Interval pi(nextafter(numbers::pi, -INFINITY),
                  nextafter(numbers::pi, INFINITY));

Interval paddedInterval(double l, double r) {
    return Interval(l - 1e-10, r + 1e-10);
}

Interval paddedInterval(double x) {
    return paddedInterval(x, x);
}

// I < J: i < j for all i in I and all j in J.
bool operator <(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return true;
    return I.r < J.l;
}

ostream& operator <<(ostream& out, const Interval& I) {
    out << "[" << I.l << ", " << I.r << "]";
    return out;
}

// Returns whether I \subseteq J.
bool isSubset(const Interval& I, const Interval& J) {
    if (I.empty())
        return true;
    if (J.empty())
        return false;
    return J.l <= I.l and I.r <= J.r;
}

// Returns the smallest interval containing all the given values.
Interval hull(const vector<double>& vec) {
    if (vec.empty())
        return EMPTY_INTERVAL;
    return Interval(*min_element(vec.begin(), vec.end()),
                    *max_element(vec.begin(), vec.end()));
}

// Returns the smallest interval containing all the given intervals.
Interval hull(const vector<Interval>& vec) {
    vector<double> vals;
    for (const Interval& I: vec) {
        if (!I.empty()) {
            vals.push_back(I.l);
            vals.push_back(I.r);
        }
    }
    return hull(vals);
}

Interval intersection(Interval I, const Interval& J) {
    I.intersect(J);
    return I;
}

// I+J
Interval operator +(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return EMPTY_INTERVAL;
    return Interval(I.l + J.l, I.r + J.r);
}

Interval& operator +=(Interval& I, const Interval& J) {
    I = I + J;
    return I;
}


// I-J
Interval operator -(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return EMPTY_INTERVAL;
    return Interval(I.l - J.r, I.r - J.l);
}

// I*J
Interval operator *(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return EMPTY_INTERVAL;
    vector<double> vec = {I.l*J.l, I.l*J.r, I.r*J.l, I.r*J.r};
    return hull(vec);
}

// -I
Interval operator -(Interval I) {
    if (I.empty())
        return EMPTY_INTERVAL;
    return Interval(-I.r, -I.l);
}

// I / J
Interval operator /(const Interval& I, const Interval& J) {
    assert(!J.containsZero());
    if (I.empty() or J.empty())
        return EMPTY_INTERVAL;
    return I * Interval(1.0 / J.r, 1.0 / J.l);
}

// {x²: x\in I}
Interval square(Interval I) {
    if (I.empty())
        return EMPTY_INTERVAL; // empty
    double a = I.l * I.l;
    double b = I.r * I.r;
    if (I.containsZero())
        return Interval(0, max(a, b));
    return Interval(min(a, b), max(a, b));
}

// sqrt(I)
Interval sqrt(Interval I) {
    if (I.empty())
        return EMPTY_INTERVAL;
    assert(I.l >= 0);
    return Interval(sqrt(I.l), sqrt(I.r));
}


// If cos(theta) \in I, then theta belongs to the returned interval.
// Returned interval is a subset of [0, pi.r]
Interval acos(Interval I) {
    I.intersect(Interval(-1, 1));
    if (I.empty())
        return EMPTY_INTERVAL;
    double l = acos(I.r);
    double r = acos(I.l);

    const static double eps = 1e-10;
    Interval J(l-eps, r+eps);
    J.intersect(Interval(0, pi.r));
    return J;
}

// If sin(theta) \in I, then theta belongs to the returned interval.
// Returned interval is a subset of [-pi.r/2, pi.r/2]
Interval asin(Interval I) {
    I.intersect(Interval(-1, 1));
    if (I.empty())
        return EMPTY_INTERVAL;

    // asin is strictly increasing, so l maps to l and r maps to r
    double l = asin(I.l);
    double r = asin(I.r);

    const static double eps = 1e-10;
    Interval J(l - eps, r + eps);
    J.intersect(Interval(-pi.r / 2.0, pi.r / 2.0));
    return J;
}


// Returns whether x+2pik belongs to I for some integer k.
// 0 <= x < 2pi.
bool containsMod2pi(Interval I, double x) {
    if (I.empty())
        return false;
    double lower = (I.l - x) / (2*M_PI);
    double upper = (I.r - x) / (2*M_PI);

    int k_min = static_cast<int>(ceil(lower));
    int k_max = static_cast<int>(floor(upper));

    return k_min <= k_max;
}

// sin(I)
Interval sin(const Interval& I) {
    if (I.empty())
        return EMPTY_INTERVAL;
    double sl = sin(I.l), sr = sin(I.r);

    double l = containsMod2pi(I, 3*M_PI/2) ? -1 : min(sl, sr);
    double r = containsMod2pi(I, M_PI/2) ? 1 : max(sl, sr);

    return Interval(l, r);
}

// cos(I)
Interval cos(const Interval& I) {
    if (I.empty())
        return EMPTY_INTERVAL;
    double sl = cos(I.l), sr = cos(I.r);

    double l = containsMod2pi(I, M_PI) ? -1 : min(sl, sr);
    double r = containsMod2pi(I, 0) ? 1 : max(sl, sr);

    return Interval(l, r);
}

// Returns 1, -1, INDETERMINATE if respectively I is an interval of positive values,
// an interval of negative values, or it contains zero.
int sign(const Interval& I) {
    if (I.l <= 0 and I.r >= 0)
        return INDETERMINATE;
    return (I.l > 0) - (I.r < 0);
}

#endif
