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

pt operator-(const pt& P) {
    return {-P.x, -P.y, -P.z};
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
    return {P.x*lambda, P.y*lambda, P.z*lambda};
}

pt operator *(double lambda, const pt& P) {
    return {P.x*lambda, P.y*lambda, P.z*lambda};
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
    
    TangentVector(): P(1, 0, 0), v(0, 1, 0) { check(); }
    
    TangentVector(pt _P, pt _v): P(_P), v(_v) {
        check();
    }
    
    pt equatorNormal() { return P^v; }
    
    void rotate(double theta) {
        v = cos(theta)*v + sin(theta)*equatorNormal();
        check();
    }
    void exp(double theta) {
        pt P2 = cos(theta)*P + sin(theta)*v;
        v = cos(theta)*v - sin(theta)*P;
        P = P2;
        
        check();
    }
};

bool onEquator(pt p, TangentVector tv) {
    static const double eps = 1e-6;
    return abs(tv.equatorNormal() * p) < eps;
}

double invertExp(TangentVector tv, pt q) {
    assert(onEquator(q, tv));
    double theta = atan2(q*tv.v, q*tv.P);
    if (theta < 0) theta += 2*M_PI;
    return theta;
}

// Returns the smallest positive length such that tv.exp(len).P is 
// orthogonal to q.
// Assumption: tv.P ^ tv.v is independent from q.
double intersectWithEquator(TangentVector tv, pt q) {
    double theta = atan2(-(tv.P*q), tv.v*q);
    while (theta < 0) theta += M_PI;
    return theta;
}

struct GeodesicSegment:TangentVector {
    double len;
    GeodesicSegment(TangentVector tv, double _len): 
            TangentVector(tv), len(_len) {
        assert(0 < len and len < 2*M_PI);
    }
    GeodesicSegment(pt P, pt v, double _len): 
            GeodesicSegment(TangentVector(P, v), len) {}
    
    pt from() { return P; }
    pt to() { 
        TangentVector tmp(*this);
        tmp.exp(len);
        return tmp.P;
    }
};

// Returns whether Q belongs (approximately) to the geodesic segment a.
// Assumption: Q and the geodesic segment belong (approximately) to
// the same plane.
bool onGeodesicSegment(pt Q, GeodesicSegment a) {
    assert(onEquator(Q, a));
    static const double eps = 1e-6;
    
    double theta = atan2(Q * a.v, Q * a.from());
    for (int i = -3; i <= 3; i++) {
        double t = theta + i * (2*M_PI);
        if (-eps <= t and t <= a.len + eps)
            return true;
    }
    return false;
}

// If it returns false, the geodesic segments are guaranteed to be
// disjoint.
// If it returns true, the geodesic segments are guaranteed to have
// distance very small. The function may return true even if the two
// geodesic segments are disjoint (but very close to intersecting).
bool areIntersecting(GeodesicSegment a, GeodesicSegment b) {
    static const double eps = 1e-6;
    
    pt pi_a = a.equatorNormal();
    pt pi_b = b.equatorNormal();
    
    if (norm(pi_a - pi_b) < eps or norm(pi_a + pi_b) < eps) { 
        // Almost on the same equator.
        return onGeodesicSegment(a.from(), b) or
               onGeodesicSegment(a.to(), b) or
               onGeodesicSegment(b.from(), a) or
               onGeodesicSegment(b.to(), a);
    } 
    // On different equators, here the problem is more stable.
    // The drawback is that taking care of arcs longer than pi is 
    // painful.
    double maybe_len_a = intersectWithEquator(a, pi_b);
    if (M_PI - eps <= maybe_len_a)
        maybe_len_a -= M_PI;

    for (double len = maybe_len_a; len <= a.len+eps; len += M_PI) {
        GeodesicSegment a2 = a;
        a2.len = len;
        if (onGeodesicSegment(a2.to(), b))
            return true;
    }
    return false;
}

// Returns a, given A, B, C.
// All angles are not 0 or pi.
double sphericalCosineRule(double A, double B, double C) {
    double cosa = (cos(A) + cos(B)*cos(C)) / (sin(B) * sin(C));
    if (cosa < -1 or cosa > 1)
        return -1;
    return acos(cosa);
}

// We start in vt1.
// We run exp for len1
// We turn by angle1.
// We are in vt2.
// We run exp for len2
// We turn by angle2.
// We are in vt3.
// We run exp for len3
// We turn by angle3.
// We are in vt4
array<double, 3> completePolygon(TangentVector vt1, TangentVector vt4,
        double angle1, double angle2, double angle3) {
    vt4.rotate(-angle3);
    // The assumption is that they are not on (approximately) the
    // same equator.
    if (onEquator(vt4.P, vt1) and onEquator(vt4.v, vt1))
        return {-1, -1, -1};
        
    double X1 = intersectWithEquator(vt1, vt4.equatorNormal());
    TangentVector vt1_bis = vt1;
    vt1_bis.exp(X1);
    
    vt4.v = -vt4.v;
    double X4 = invertExp(vt4, vt1_bis.P);
    TangentVector vt4_bis = vt4;
    vt4_bis.exp(X4);
    vt4.v = -vt4.v;
    vt4_bis.v = -vt4_bis.v;
    
    // Now, vt1_bis and vt4_bis are centered in the same point.
    // We shall start from vt1_bis, go for len1, turn by angle1,
    // go for len2, turn by angle2, go for len3, end at vt4_bis.
    double angle0 = ?
}


#endif
