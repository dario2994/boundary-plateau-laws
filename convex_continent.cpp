#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

struct Interval {
    double l, r;
    
    Interval(): l(-1e9), r(1e9) {}
    Interval(double _l, double _r): l(_l), r(_r) {}
    
    bool empty() { return r < l; }
    
    void intersect(Interval I) {
        l = max(l, I.l);
        r = min(r, I.r);
    }
    
    double len() { return empty() ? 0.0 : (r-l); }
};

// I is strictly inside [0, \pi].
// One side of a quadrilateral belongs to I.
// The other side belongs to the interval returned by this function.
Interval quadrilateralConstraint(Interval I) {
    if (I.empty())
        return Interval(0, -1);
    const static double eps = 1e-8;
    static auto f = [](double x) { return 2*atan(0.5 / tan(x/2)); };
    return Interval(f(I.r) - eps, f(I.l) + eps);
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


// I, J are strictly inside [0, \ell_T].
// If AB \in I, BC\in J, then DE belongs to the interval returned by 
// this function.
Interval pentagon1Constraint(Interval I, Interval J) {
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

// I, J are strictly inside [0, \ell_T].
// AB \in I, DE\in J.
// Then BC belongs to the interval returned by this function.
Interval pentagon2Constraint(Interval I, Interval J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
    
    static auto g = [](double x) {
        double theta = atan(2*sin(x) / (1+cos(x)));
        double lambda = 2*sin(x) / sin(theta);
        return pair<double,double>{theta, lambda};
    };
    
    static auto f = [&](double x, double z) {
        auto pp = g(x);
        double lambda = pp.second;
        return (2*cos(z)-1.0/3-cos(x))/lambda;
    };
    
    double low = f(I.r, J.l);
    double high = f(I.l, J.r);
    
    Interval K = acos(Interval(low, high));
    K.l -= g(I.l).first;
    K.r -= g(I.r).first;
    
    const static double ellT = acos(-1.0/3);
    K.intersect(Interval(0, ellT));
    
    return K; 
}


int main() {
    
}

