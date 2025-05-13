#ifndef INTERVAL_ARITHMETIC
#define INTERVAL_ARITHMETIC

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

// Struct maintaining a closed interval [l, r].
// In many of our applications, the interval is a subset of [0, pi].
struct Interval {
    double l, r;
    
    Interval(): l(-1e9), r(1e9) {}
    Interval(double _l, double _r): l(_l), r(_r) {}
    Interval(double x): l(x-1e-10), r(x+1e-10) {}
    
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
    bool containsZero() const { return l <= 0 and 0 <= r; } // TODO: Remove.
};

ostream& operator <<(ostream& out, const Interval& I) {
    out << "[" << I.l << ", " << I.r << "]";
    return out;
}

// Returns the smallest interval containing all the given intervals.
Interval hull(vector<double> vec) { 
    return Interval(*min_element(vec.begin(), vec.end()), 
                    *max_element(vec.begin(), vec.end()));
}

Interval intersection(Interval I, const Interval& J) {
    I.intersect(J);
    return I;
}

// If 0 does not belong to I, returns 1/I.
Interval inverse(const Interval I) {
    assert(!I.containsZero());
    return Interval(1.0 / I.r, 1.0 / I.l);
}

// I+J
Interval operator +(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
    return Interval(I.l + J.l, I.r + J.r);
}

Interval& operator +=(Interval& I, const Interval& J) {
    I = I + J;
    return I;
}


// I-J
Interval operator -(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
    return Interval(I.l - J.r, I.r - J.l);
}

// I*J
Interval operator *(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
    vector<double> vec = {I.l*J.l, I.l*J.r, I.r*J.l, I.r*J.r};
    return hull(vec);
}

// -I
Interval operator-(Interval I) {
    if (I.empty())
        return Interval(0, -1);
    return Interval(-I.r, -I.l);
}


// If cos(theta) \in I, then theta belongs to the returned interval.
// Returned interval is a subset of [0, pi]
Interval acos(Interval I) {
    I.intersect(Interval(-1, 1));
    if (I.empty())
        return Interval(0, -1); // empty
    double l = acos(I.r);
    double r = acos(I.l);
    
    const static double eps = 1e-7;
    Interval J(l-eps, r+eps);
    J.intersect(Interval(0, M_PI));
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
        return Interval(0, -1);
    double sl = sin(I.l), sr = sin(I.r);
    
    double l = containsMod2pi(I, 3*M_PI/2) ? -1 : min(sl, sr);
    double r = containsMod2pi(I, M_PI/2) ? 1 : max(sl, sr);
    
    return Interval(l, r);
}

// cos(I)
Interval cos(const Interval& I) {
    if (I.empty())
        return Interval(0, -1);
    double sl = cos(I.l), sr = cos(I.r);
    
    double l = containsMod2pi(I, M_PI) ? -1 : min(sl, sr);
    double r = containsMod2pi(I, 0) ? 1 : max(sl, sr);
    
    return Interval(l, r);
}

double totalLen(const vector<Interval>& vec) {
    double res = 0;
    for (const Interval& I: vec) {
        if (I.empty())
            return -1;
        res += I.len();
    }
    return res;
}

int sign(const Interval& I) {
    return (I.l > 0) - (I.r < 0);
}

#endif
