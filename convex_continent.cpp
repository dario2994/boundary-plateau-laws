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

struct Quaternion{
    double q[4];
    Quaternion(): q{0,0,0,0} {}
    Quaternion(double a, double b, double c, double d): q{a,b,c,d} {}
    Quaternion(double r): q{r, 0, 0, 0} {}
    
    double& operator [](int i) { return q[i]; }
    const double& operator [](int i) const { return q[i]; }
};

Quaternion operator +(const Quaternion& x, const Quaternion& y) {
    Quaternion z;
    for (int i = 0; i < 4; i++)
        z[i] = x[i] + y[i];
    return z; 
}

Quaternion operator -(const Quaternion& x, const Quaternion& y) {
    Quaternion z;
    for (int i = 0; i < 4; i++)
        z[i] = x[i] - y[i];
    return z; 
}

Quaternion operator *(const Quaternion& x, const Quaternion& y) {
    Quaternion z;
    z[0] = x[0]*y[0] - x[1]*y[1] - x[2]*y[2] - x[3]*y[3];
    z[1] = x[0]*y[1] + x[1]*y[0] + x[2]*y[3] - x[3]*y[2];
    z[2] = x[2]*y[0] + x[0]*y[2] + x[3]*y[1] - x[1]*y[3];
    z[3] = x[3]*y[0] + x[0]*y[3] + x[1]*y[2] - x[2]*y[1];
    return z;
}

Quaternion& operator *=(Quaternion& x, const Quaternion& y) {
    x = x * y;
    return x;
}

Quaternion conj(const Quaternion& z) {
    return Quaternion(z[0], -z[1], -z[2], -z[3]);
}

double abs(const Quaternion& z) {
    double res = 0;
    for (int i = 0; i < 4; i++)
        res += pow(z[i], 2);
    return res;
}


ostream& operator <<(ostream& out, const Quaternion& q) {
    out << q[0] << " + i" << q[1] << " + j" << q[2] << " + k" << q[3];
    return out;
}

void checkGeneralizedPolygon(
    const vector<double>& angles, const vector<double>& sides) {
    int k = ssize(angles);
    vector<double> theta(k+1);
    theta[0] = 0;
    for (int i = 0; i < k; i++) theta[i+1] = theta[i] + angles[i] - M_PI;
    
    Quaternion res = 1;
    for (int i = 0; i < k; i++) {
        Quaternion q(0, 0, -sin(theta[i]), cos(theta[i]));
        res *= 1 + tan(sides[i]/2) * q;
    }
    cout << res << endl;
    
    cout << res[1]/res[0] << " " << tan(-theta[k]/2) << endl;
}

#include "spherical_geometry.cpp"


int main() {
    double a120 = 2*M_PI / 3;
    double a240 = 2*a120;
    vector<double> angles = {a120, a120, a120, a120, a120, a240};
    for (int _ = 0; _ < 100; _++) {
        auto maybe = completeGeneralizedPolygon(angles, 
            {(rand() % 100)/50.0, (rand() % 100)/50.0, (rand() % 100)/50.0, 0, 0, 0}, {3, 4, 5});
        if (maybe.empty())
            continue;
        
        vector<double> sides = maybe[0];
        //~ cout << sides[0] << " " << sides[1] << " " << sides[2] << " " << sides[3] << endl;
        checkGeneralizedPolygon(angles, sides);
    }
}

