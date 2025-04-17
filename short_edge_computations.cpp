#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;


double I(double a, double b, double c) {
    assert(0 < a and a < b);
    assert(0 < c);
    double z = (b-a) / c;
    a /= z, b /= z;
    return pow(z, 2) * (pow(1+b*b, 1.5) - pow(1+a*a, 1.5)) / 3;
}

double J(double a, double b) {
    assert(0 < a and a < b);
    return (pow(b, 2) - pow(a, 2)) / 2;
}

double perturbation_cost(double r0, double r1, double r2, double r3,
                         double l0, double l1, double l2, double t) {
    assert(0 < r0 and r0 < r1 and r1 < r2 and r2 < r3);
    assert(l1 > max(l0, l2));
    assert(l2 < l0);
    
    return l2 * J(r1, r2) + l1 * (I(r0, r1, t) + I(r2, r3, t)) - l0 * J(r0, r3);
}


struct Configuration {
    double delta;
    double L;
    double AO, AOP, QOA, AQ, QO;
    Configuration(double _delta, double _L): delta(_delta), L(_L) {
        assert(0 < L and L < 0.8);
        double cos_AO = cos(L)*cos(delta/2) + cos(2*M_PI/3)*sin(L)*sin(delta/2);
        AO = acos(cos_AO);
        
        double sin_AOP = sin(L)*sin(2*M_PI/3)/sin(AO);
        assert(sin_AOP > sin(M_PI/6));
        AOP = asin(sin_AOP);
        
        
        QOA = M_PI/2 - AOP;
        
        double sin_AQ = 2/sqrt(3) * sin(AO) * sin(QOA);
        AQ = asin(sin_AQ);
        
        double tan_QO2 = tan((AO-AQ)/2) * sin((2*M_PI/3 + QOA)/2)
                                            /
                                          sin((2*M_PI/3 - QOA)/2);
        QO = 2*atan(tan_QO2);
    }
    
    double cost(double r0, double r1, double r2, double r3) const {
        double l0 = 2*L;
        double l1 = 2 * AO;
        double l2 = QO + 2*AQ;
        double t = delta/2 + QO;
        
        double x = l0-l2;
        double y = 2*sqrt(2) * l1 * t;
        //~ cout << x << " > " << y << endl;
        cout << x - y << endl;
        
        //~ cout << l0 << " " << l1 << " " << l2 << " " << t << endl;
        
        return perturbation_cost(r0, r1, r2, r3, l0, l1, l2, t);
        
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
    //~ double L = 0.12;
    //~ double delta = 0.02;
    double L = 0.02;
    double delta = 1.0/250;
    Configuration C(delta, L);
    cout << C.cost(0.7, 1, 10, 12.5) << endl;
    
    
    // Useful for heuristics.
    auto optimize_stuff = [C](double r2) {
        double r0 = 0.5;
        double r1 = 1;
        double r3 = 3;
    
        for (int _ = 0; _ < 5; _++) {
            r3 = argmin([r0, r1, r2, C](double x) { 
                    return C.cost(r0, r1, r2, x);
            },
                        r2*(1.01), r2*5);
            r0 = argmin([r1, r2, r3, C](double x) { 
                return C.cost(x, r1, r2, r3); 
            },
                        0.1*r1, r1);
        }
        cout << r0 << " " << r1 << " " << r2 << " " << r3 << endl;
        return C.cost(r0, r1, r2, r3);
    };
}

