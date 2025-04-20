#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;
#include "spherical_geometry.cpp"


struct Configuration {
    double delta;
    double L;
    ConvexSphericalTriangle APO, AQO;

    Configuration(double _delta, double _L): delta(_delta), L(_L) {
        assert(0 < L and L < 0.8);
        //~ double cos_AO = cos(L)*cos(delta/2) + cos(2*M_PI/3)*sin(L)*sin(delta/2);
        //~ AO = acos(cos_AO);
        
        //~ double sin_AOP = sin(L)*sin(2*M_PI/3)/sin(AO);
        //~ assert(sin_AOP > sin(M_PI/6));
        //~ AOP = asin(sin_AOP);
        
        
        //~ QOA = M_PI/2 - AOP;
        
        //~ double sin_AQ = 2/sqrt(3) * sin(AO) * sin(QOA);
        //~ AQ = asin(sin_AQ);
        
        //~ double tan_QO2 = tan((AO-AQ)/2) * sin((2*M_PI/3 + QOA)/2)
                                            //~ /
                                          //~ sin((2*M_PI/3 - QOA)/2);
        //~ QO = 2*atan(tan_QO2);
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

int main() {
    double L = 0.5;
    double delta = 0.7;
    //~ double L = 0.02;
    //~ double delta = 1.0/250;
    Configuration C(delta, L);
    cout << C.cost() << endl;
}

