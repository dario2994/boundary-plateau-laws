#ifndef SPHERICAL_GEOMETRY
#define SPHERICAL_GEOMETRY

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;


struct pt {
    double x, y, z;
    pt(): x(0), y(0), z(0) {}
    pt(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}
};

ostream& operator <<(ostream& out, const pt& P) {
    out << "(" << P.x << ", " << P.y << ", " << P.z << ")";
    return out;
}

pt operator +(const pt& P, const pt& Q) {
    return {P.x+Q.x, P.y+Q.y, P.z+Q.z};
}

pt operator -(const pt& P, const pt& Q) {
    return {P.x-Q.x, P.y-Q.y, P.z-Q.z};
}

pt& operator +=(pt& P, const pt& Q) {
    P.x += Q.x, P.y += Q.y, P.z += Q.z;
    return P;
}

pt& operator -=(pt& P, const pt& Q) {
    P.x -= Q.x, P.y -= Q.y, P.z -= Q.z;
    return P;
}

pt operator *(const pt& P, double lambda) {
    return {P.x*lambda, P.y*lambda, P.y*lambda};
}

pt operator *(double lambda, const pt& P) {
    return {P.x*lambda, P.y*lambda, P.y*lambda};
}

pt operator /(const pt& P, double lambda) {
    return P * (1/lambda);
}

double operator *(const pt& P, const pt& Q) {
    return P.x*Q.x + P.y*Q.y + P.z*Q.z;
}

pt operator ^(const pt& P, const pt& Q) {
    return {P.y*Q.z - P.z*Q.y,
            P.z*Q.x - P.x*Q.z,
            P.x*Q.y - P.y*Q.x};
}

double norm2(const pt& P) { return P*P; }
double norm(const pt& P) { return sqrt(norm2(P)); }
double dist2(const pt& P, const pt& Q) { return norm2(P-Q); }
double dist(const pt& P, const pt& Q) { return norm(P-Q); }
pt normalize(const pt& P) { return P / norm(P); }

bool onSphere(const pt& P) {
    static const double eps = 1e-8;
    double z = norm2(P);
    return 1 - eps <= z and z <= 1 + eps;
}

bool areOrthogonal(const pt& P, const pt& Q) {
    static const double eps = 1e-8;
    return abs(P*Q) <= eps;
}

struct TangentVector {
    // Position, velocity.
    pt P, v;
    
    void check() {
        assert(onSphere(P));
        assert(onSphere(v));
        assert(areOrthogonal(P, v));
    }
    
    TangentVector(pt _P, pt _v): P(_P), v(_v) {
        check();
    }
    
    void rotate(double theta) {
        pt u = P^v;
        v = cos(theta)*v + sin(theta)*u;
        check();
    }
    void exp(double theta) {
        P = cos(theta)*P + sin(theta)*v;
        check();
    }
};



#endif
