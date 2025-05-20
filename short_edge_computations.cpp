#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;
#include "spherical_geometry.cpp"


struct Configuration {
    double L;
    double delta;
    ConvexSphericalTriangle APO, AQO;

    Configuration(double _L, double _delta): L(_L), delta(_delta) {
        assert(0 < L and L < 0.8);
        
        APO = ConvexSphericalTriangle(
            {{'a', delta/2}, {'c', L}, {'B', 2*M_PI/3}});
        AQO = ConvexSphericalTriangle(
            {{'b', APO.get('b')}, {'B', 2*M_PI/3}, {'C', M_PI/2 - APO.get('C')}});
    }
    
    double cost() const {
        double q1 = 2*L;
        double q2 = AQO.get('a') + 2*AQO.get('c');
        double q3 = 2 * (APO.area() + AQO.area());
        
        return (q1-q2)/2 - q3;
    }
};


// Useful for heuristics.
double argmin(function<double(double)> f, double l, double r) {
    for (int _ = 0; _ < 100; _++) {
        double m1 = (2*l+r)/3;
        double m2 = (2*r+l)/3;
        if (f(m1) < f(m2))
            r = m2;
        else
            l = m1;
    }
    return l;
}

double argmin_naive(function<double(double)> f, double l, double r) {
    double ans = 0;
    double f_ans = 1e18;
    for (double m = l; m < r; m += (r-l)/3000) {
        double f_m = f(m);
        if (f_m < f_ans) {
            f_ans = f_m;
            ans = m;
        }
    }
    return ans;
}

int main() {
    double L = 0.2;
    //~ double delta = 0.02;
    //~ double L = argmin_naive([&](double z) { return -Configuration(delta, z).cost(); }, delta*3, 0.5);
    //~ double L = 0.02;
    //~ double delta = 1.0/250;
    cout << L << endl;
    Configuration C(L, L/10);
    //~ cout << C.AQO.get('c') << " " << C.AQO.get('a') << endl;
    //~ cout << C.APO.area() << " " << C.AQO.area() << endl;
    cout << C.cost() << endl;
}

