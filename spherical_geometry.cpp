#ifndef SPHERICAL_GEOMETRY
#define SPHERICAL_GEOMETRY

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

static const double eps = 1e-7;

// Clip angle to [0, 2pi).
double normalizeAngle(double theta) {
    while (theta < 0) theta += 2*M_PI;
    while (theta >= 2*M_PI) theta -= 2*M_PI;
    return theta;
}

// Class representing a point in R^3.
struct pt {
    double x, y, z;
    pt(): x(0), y(0), z(0) {}
    pt(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}
    void serialize(ostream& out) const {
        out << x << " " << y << " " << z;
    }
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
    double z = norm2(P);
    return 1 - eps <= z and z <= 1 + eps;
}

bool areOrthogonal(const pt& P, const pt& Q) {
    return abs(P*Q) <= eps;
}

// Class representing a unit tangent vector to the sphere.
struct TangentVector {
    // Position, velocity.
    // |P| = |V| = 1.
    // P and V are orthogonal.
    pt P, V;
    
    void check() const {
        assert(onSphere(P));
        assert(onSphere(V));
        assert(areOrthogonal(P, V));
    }
    
    TangentVector(): P(1, 0, 0), V(0, 1, 0) { check(); }
    
    TangentVector(pt _P, pt _V): P(_P), V(_V) { check(); }
    
    // Returns a unit vector orthogonal to the equator identified by
    // this tangent vector.
    pt equatorNormal() const { return P^V; }
    
    // Whether Q belongs to the equator identified by this tangent vector.
    // Assumption: |Q|=1.
    bool onEquator(pt Q) const {
        return abs(equatorNormal() * Q) < eps;
    }
    
    // Rotate the direction of V by theta.
    TangentVector rot(double theta) const {
        return {P, cos(theta)*V + sin(theta)*equatorNormal()};
    }
    
    // Apply the exponential map for time theta.
    // This changes both P and V.
    TangentVector exp(double theta) const {
        return {cos(theta)*P + sin(theta)*V, 
                cos(theta)*V - sin(theta)*P};
    }
    
    // Find theta such that exp(theta).P = Q.
    // Assumption: |Q|=1, onEquator(Q).
    double inverseExp(pt Q) const {
        assert(onEquator(Q));
        double theta = atan2(Q * V, Q * P);
        return normalizeAngle(theta);
    }
    
    // Find theta such that rot(theta).V = U
    // Assumption: |U|=1, areOrthogonal(U, P)
    double inverseRot(pt U) const {
        double theta = atan2(U * equatorNormal(), U * V);
        return normalizeAngle(theta);
    }
    
    void serialize(ostream& out) const {
        P.serialize(out);
        out << " ";
        V.serialize(out);
    }
};

// Returns the smallest positive length such that tv.exp(len).P is 
// orthogonal to equator_normal.
// Assumption: 
//     |equator_normal|=1
//     tv.equatorNormal() is independent from equator_normal.
double intersectWithEquator(TangentVector tv, pt equator_normal) {
    // This is a "bad" way to take care of the instability.
    if (norm2(equator_normal^tv.equatorNormal()) < eps)
        return 0;
    double theta = tv.inverseExp(normalize(equator_normal^tv.equatorNormal()));
    if (theta > M_PI) theta -= M_PI;
    return theta;
}

double distance(pt P, pt Q) {
    if (norm2(P-Q) < eps) return norm(P-Q);
    if (norm2(P+Q) < eps) return M_PI - norm(P+Q);
    pt equator = normalize(P^Q);
    TangentVector tv(P, P^equator);
    double len = tv.inverseExp(Q);
    return min(len, 2*M_PI - len);
}

// Class representing a spherical geodesic arc with length in (0, 2pi). 
struct GeodesicSegment:TangentVector {
    double len;
    
    GeodesicSegment(): TangentVector(), len(-1) {}
    
    GeodesicSegment(TangentVector tv, double _len): 
            TangentVector(tv), len(_len) {
        assert(0 <= len and len < 2*M_PI);
    }
    GeodesicSegment(pt P, pt v, double _len): 
            GeodesicSegment(TangentVector(P, v), _len) {}
    
    pt from() const { return P; }
    pt to() const { return exp(len).P; }
    
    // Whether Q belongs to the geodesic segment.
    // Assumption: |Q|=1, onEquator(Q).
    bool contains(pt Q) const {
        double theta = inverseExp(Q);
        return theta >= (2*M_PI - eps) or theta <= len + eps;
    }
    
    void serialize(ostream& out) const {
        TangentVector::serialize(out);
        out << " ";
        out << len;
    }
};

// If it returns false, the geodesic segments are guaranteed to be
// disjoint.
// If it returns true, the geodesic segments are guaranteed to have
// distance very small. The function may return true even if the two
// geodesic segments are disjoint (but very close to intersecting).
bool areIntersecting(GeodesicSegment a, GeodesicSegment b) {
    pt equator_a = a.equatorNormal();
    pt equator_b = b.equatorNormal();
    
    if (norm(equator_a - equator_b) < eps 
        or norm(equator_a + equator_b) < eps) { 
        // Almost on the same equator.
        return b.contains(a.from()) or
               b.contains(a.to()) or
               a.contains(b.from()) or
               a.contains(b.to());
    } 
    // On different equators, here the problem is more stable.
    // The drawback is that taking care of arcs longer than pi is 
    // painful.
    double maybe_len_a = intersectWithEquator(a, equator_b);
    if (M_PI - eps <= maybe_len_a)
        maybe_len_a -= M_PI;

    for (double len = maybe_len_a; len <= a.len+eps; len += M_PI) {
        if (b.contains(a.exp(len).P))
            return true;
    }
    return false;
}

// Whether there exists a spherical convex triangle with angles A, B, C.
// Assumption: 0 < A, B, C < 2pi.
bool validateConvexTriangle(double A, double B, double C) {
    bool condition1 = M_PI + 2*min({A, B, C}) > A + B + C;
    bool condition2 = M_PI < A + B + C and A + B + C < 3*M_PI;
    auto is_convex = [&](double theta) { return 0 < theta and theta < M_PI; };
    bool condition3 = is_convex(A) and is_convex(B) and is_convex(C);
    
    return condition1 and condition2 and condition3;
}

// Returns a, given A, B, C.
// Assumption: validateConvexTriangle(A, B, C).
double sphericalCosineRule(double A, double B, double C) {
    assert(validateConvexTriangle(A, B, C));
    double cosa = (cos(A) + cos(B)*cos(C)) / (sin(B) * sin(C));
    assert(-1 <= cosa and cosa <=1);
    return acos(cosa);
}

// Returns the possible side lengths of a generalized spherical triangle
// with angles A, B, C.
// Assumption: A, B, C \in (0, 2\pi)\setminus\{\pi\}.
vector<array<double,3>> computeTriangleSides(double A, double B, double C) {
    int cnt = (A < M_PI) + (B < M_PI) + (C < M_PI);
    
    if (cnt % 2 == 0)
        return computeTriangleSides(2*M_PI - A, 2*M_PI - B, 2*M_PI - C);
    
    if (cnt == 3) {
        if (!validateConvexTriangle(A, B, C))
            return {};
        double a = sphericalCosineRule(A, B, C);
        double b = sphericalCosineRule(B, C, A);
        double c = sphericalCosineRule(C, A, B);
        
        return {{a, b, c}, {2*M_PI-a, 2*M_PI-b, 2*M_PI-c}};
    }
    
    // Exactly one in A, B, C is lower than pi.
    // Which one? The one corresponding to id according to A->0, B->1, C->2.
    int id = (B < M_PI) + 2*(C < M_PI);
    auto res = 
        computeTriangleSides(fmod(A, M_PI), fmod(B, M_PI), fmod(C, M_PI));
    for (auto& tri: res) tri[id] = 2*M_PI - tri[id];
    
    return res;
}

// Any transformation of the unit tangent of S^2 that commutes with
// SO(R^3).
struct SimpleMap {
    double t1, alpha, t2;
    SimpleMap(): t1(0), alpha(0), t2(0) {}
    SimpleMap(double _t1, double _alpha, double _t2):
        t1(_t1), alpha(_alpha), t2(_t2) {}
        
    TangentVector operator()(const TangentVector& tv) const {
        return tv.exp(t1).rot(alpha).exp(t2);
    }
    
    SimpleMap exp(double t) const {
        return {t1, alpha, t2+t};
    }
    
    SimpleMap rot(double beta) const {
        TangentVector v0;
        TangentVector v3 = (*this)(v0).rot(beta);
        pt goal_equator = v3.equatorNormal();
        
        // This one is the only unstable operation.
        double new_t1 = intersectWithEquator(v0, v3.equatorNormal());
        TangentVector v1 = v0.exp(new_t1);
        assert(v1.onEquator(goal_equator));
        pt half_dir = goal_equator^v1.P;
        double new_alpha = v1.inverseRot(half_dir);
        TangentVector v2 = v1.rot(new_alpha);
        double new_t2 = v2.inverseExp(v3.P);
        
        return {new_t1, new_alpha, new_t2};
    }
};

// Finds all the possible side lengths that coincide with the given one
// apart from the three indexes missing_sides, and which create a 
// generalized spherical polygon when coupled with angles.
// 
// sides[0], angles0, sides[1], angles[1], ..., sides[n-1], angles[n-1].
vector<vector<double>> completePolygon(
        const vector<double>& angles, vector<double> sides,
        vector<int> missing_sides) {
    int n = ssize(angles);
    assert(n >= 3);
    assert(ssize(sides) == n);
    assert(ssize(missing_sides) == 3);
    
    sort(missing_sides.begin(), missing_sides.end());
    assert(0 <= missing_sides[0]
           and missing_sides[0] < missing_sides[1]
           and missing_sides[1] < missing_sides[2]
           and missing_sides[2] < n);
        
    vector<SimpleMap> strange_angles(3);
    for (int i = 0; i < 3; i++) {
        SimpleMap sm = SimpleMap().rot(angles[missing_sides[i]]);
        for (int j = missing_sides[i]+1; 
                j != missing_sides[i+1]; 
                j = (j+1) % n) {
            sm = sm.exp(sides[j]).rot(angles[j]);
        }
        strange_angles[i] = sm;
    }
    vector<array<double, 3>> tri_lengths = computeTriangleSides(
        strange_angles[0].alpha,
        strange_angles[1].alpha,
        strange_angles[2].alpha);
    
    vector<vector<double>> res;
    for (auto tri: tri_lengths) {
        for (int i = 0; i < 3; i++) {
            sides[missing_sides[i]] = 
                tri[i] - strange_angles[i].t1 - strange_angles[i].t2;
        }
        res.push_back(sides);
    }
    return res;
}

// We complete a polygon by finding the lengths of the last three
// sides.
//
// We start in vt1.
// We run exp for len1
// We turn by angle1-pi.
// We are in vt2.
// We run exp for len2
// We turn by angle2-pi.
// We are in vt3.
// We run exp for len3
// We turn by angle3-pi.
// We are in vt4
vector<array<GeodesicSegment, 3>> completePolygon(
        const TangentVector& vt1, TangentVector vt4,
        double angle1, double angle2) {    
    vt4.V = -vt4.V;
    
    // The assumption is that they are not on (approximately) the
    // same equator.
    if (vt1.onEquator(vt4.P) and vt1.onEquator(vt4.V))
        return {};
        
    double X1 = intersectWithEquator(vt1, vt4.equatorNormal());
    const TangentVector& vt1_pole = vt1.exp(X1);
    
    double X3 = vt4.inverseExp(vt1_pole.P);
    TangentVector vt4_pole = vt4.exp(X3);
    
    assert(norm(vt1_pole.P - vt4_pole.P) < eps);
    
    // Now, vt1_pole and vt4_pole are centered in the same point.
    // They form an angle = angle0.
    double angle0 = vt4_pole.inverseRot(vt1_pole.V);
    vector<array<double, 3>> side_lengths = 
        computeTriangleSides(angle0, angle1, angle2);
    
    vector<array<GeodesicSegment, 3>> res;
    for (auto tri: side_lengths) {
        double len1 = normalizeAngle(tri[2] + X1);
        double len2 = tri[0];
        double len3 = normalizeAngle(tri[1] + X3);
        GeodesicSegment g1(vt1, len1);
        GeodesicSegment g2(g1.exp(g1.len).rot(angle1 - M_PI), len2);
        GeodesicSegment g3(g2.exp(g2.len).rot(angle2 - M_PI), len3);
        res.push_back(array<GeodesicSegment, 3>{g1, g2, g3});
    }
    return res;
}


#endif
