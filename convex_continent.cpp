#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ULL;

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
    
    double len() { return empty() ? 0.0 : (r-l); }
    bool containsZero() { return l <= 0 and 0 <= r; }
};

ostream& operator <<(ostream& out, const Interval& I) {
    out << "[" << I.l << ", " << I.r << "]";
    return out;
}

Interval hull(vector<double> vec) { 
    return Interval(*min_element(vec.begin(), vec.end()), 
                    *max_element(vec.begin(), vec.end()));
}

// Interval arithmetic.
Interval operator +(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
    return Interval(I.l + J.l, I.r + J.r);
}

Interval operator -(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
    return Interval(I.l - J.r, I.r - J.l);
}

Interval operator *(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
    vector<double> vec = {I.l*J.l, I.l*J.r, I.r*J.l, I.r*J.r};
    return hull(vec);
}

Interval operator-(Interval I) {
    if (I.empty())
        return Interval(0, -1);
    return Interval(-I.r, -I.l);
}

void enforceTriangleConstraint(vector<Interval*> sides) {
    assert(ssize(sides) == 3);
    
    static const double L = acos(-1.0/3);
    static const Interval I(L-1e-9, L+1e-9); 
    for (int i = 0; i < 3; i++) {
        sides[i]->intersect(I);
    }
}

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

void enforceQuadrilateralConstraint(vector<Interval*> sides) {
    assert(ssize(sides) == 4);
    
    sides[0]->intersect(*sides[2]);
    sides[2]->intersect(*sides[0]);
    
    sides[1]->intersect(*sides[3]);
    sides[3]->intersect(*sides[1]);
    
    Interval even = quadrilateralConstraint(*sides[1]);
    Interval odd = quadrilateralConstraint(*sides[0]);
    
    sides[0]->intersect(even);
    sides[2]->intersect(even);
    sides[1]->intersect(odd);
    sides[3]->intersect(odd);
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
Interval pentagon1Constraint(const Interval& I, const Interval& J) {
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
Interval pentagon2Constraint(const Interval& I, const Interval& J) {
    if (I.empty() or J.empty())
        return Interval(0, -1);
        
    static auto get_theta = [](double x) {
        return atan(2*sin(x) / (1+cos(x)));
    };
    
    static auto get_lambda = [](double x) {
        return 2*sin(x) / sin(get_theta(x));
    };
    
    static auto f = [&](double x, double z) {
        return (2*cos(z)-1.0/3-cos(x))/get_lambda(x);
    };
    
    double low = f(I.l, J.r);
    double high = f(I.r, J.l);
    
    //~ cout << low << " " << high << endl;
    
    // Using acos on [low, high] is not necessarily correct here, hence
    // we do the computation explicitly.
    
    if (low > 1 or high < -1)
        return Interval(0, -1); // empty
    
    double l = (high < 1) ?  acos(high) : 0;
    l -= get_theta(I.r);
    double r = (low > -1) ? acos(low) : M_PI;
    r -= get_theta(I.l);
    //~ cout << l << " " << r << endl;
    Interval K(l, r);
    
    const static double ellT = acos(-1.0/3);
    K.intersect(Interval(0, ellT));
    
    //~ cout << I << " " << J << " " << K << endl;
    
    return K; 
}

void enforcePentagonConstraint(vector<Interval*> sides) {
    auto s = [&](int id) {
        return sides[id % 5];
    };
    
    auto print_state = [&]() {
        cout << "enforcePentagonConstraint" << endl;
        for (int i = 0; i < 5; i++)
            cout << *s(i) << " ";
        cout << endl;
    };
    
    //~ print_state();
    
    auto compute_len = [&]() {
        double res = 0;
        for (int i = 0; i < 5; i++) {
            if (s(i)->empty())
                return -1.0;
            res += sides[i]->len();
        }
        return res;
    };
    double curr_len = compute_len();
    while (curr_len >= 1e-8) {
        for (int i = 0; i < 5; i++) {
            s(i+3)->intersect(pentagon1Constraint(*s(i), *s(i+1)));
            s(i)->intersect(pentagon2Constraint(*s(i+1), *s(i+3)));
            s(i)->intersect(pentagon2Constraint(*s(i+4), *s(i+2)));
        }
        //~ print_state();
        double new_len = compute_len();
        if (new_len > 0.99*curr_len)
            break;
        curr_len = new_len;
    }
    return;
}

template<typename T>
struct Quaternion{
    T q[4];
    Quaternion(): q{T(0),T(0),T(0),T(0)} {}
    Quaternion(T a, T b, T c, T d): q{a,b,c,d} {}
    Quaternion(T r): q{r, T(0), T(0), T(0)} {}
    
    T& operator [](int i) { return q[i]; }
    const T& operator [](int i) const { return q[i]; }
};

template<typename T>
Quaternion<T> operator +(const Quaternion<T>& x, const Quaternion<T>& y) {
    Quaternion<T> z;
    for (int i = 0; i < 4; i++)
        z[i] = x[i] + y[i];
    return z; 
}

template<typename T>
Quaternion<T> operator -(const Quaternion<T>& x, const Quaternion<T>& y) {
    Quaternion<T> z;
    for (int i = 0; i < 4; i++)
        z[i] = x[i] - y[i];
    return z; 
}

template<typename T>
Quaternion<T> operator *(const Quaternion<T>& x, const Quaternion<T>& y) {
    Quaternion<T> z;
    z[0] = x[0]*y[0] - x[1]*y[1] - x[2]*y[2] - x[3]*y[3];
    z[1] = x[0]*y[1] + x[1]*y[0] + x[2]*y[3] - x[3]*y[2];
    z[2] = x[2]*y[0] + x[0]*y[2] + x[3]*y[1] - x[1]*y[3];
    z[3] = x[3]*y[0] + x[0]*y[3] + x[1]*y[2] - x[2]*y[1];
    return z;
}

template<typename T>
Quaternion<T>& operator *=(Quaternion<T>& x, const Quaternion<T>& y) {
    x = x * y;
    return x;
}

template<typename T>
Quaternion<T> conj(const Quaternion<T>& z) {
    return Quaternion<T>(z[0], -z[1], -z[2], -z[3]);
}

template<typename T>
T abs(const Quaternion<T>& z) {
    T res = 0;
    for (int i = 0; i < 4; i++)
        res += z[i]*z[i];
    return res;
}


template<typename T>
ostream& operator <<(ostream& out, const Quaternion<T>& q) {
    out << q[0] << " + i" << q[1] << " + j" << q[2] << " + k" << q[3];
    return out;
}

void checkGeneralizedPolygon(
    const vector<double>& angles, const vector<double>& sides) {
    int k = ssize(angles);
    vector<double> theta(k+1);
    theta[0] = 0;
    for (int i = 0; i < k; i++) theta[i+1] = theta[i] + angles[i] - M_PI;
    
    Quaternion<Interval> res = Interval(1);
    for (int i = 0; i < k; i++) {
        Quaternion<Interval> q(0, 0, -sin(theta[i]), cos(theta[i]));
        double l = tan(sides[i]/2);
        Interval Il(l-0.001, l+0.001);
        res *= Quaternion<Interval>(1) + Quaternion<Interval>(Il) * q;
    }
    cout << res << endl;
    
    cout << res[1]*cos(-theta[k]/2)-res[0]*sin(-theta[k]/2) << endl;
    //~ cout << res[1]/res[0] << " " << tan(-theta[k]/2) << endl;
}

Interval tan(const Interval& I) {
    if (I.empty())
        return Interval(0, -1);
    assert(-M_PI/2+1e-5 < I.l and I.r < M_PI/2-1e-5);
    return Interval(tan(I.l), tan(I.r));
}

Interval atan(const Interval& I) {
    if (I.empty())
        return Interval(0, -1);
    assert(0 <= I.l);
    return Interval(atan(I.l), atan(I.r));
}


// Returns the range of positive t such that line(1, tx, ty)
// intersects I x J x K. 
Interval subproblemInterval(double x, double y, 
        const Interval& I, const Interval& J, const Interval& K) {
    if (!(y*J - x*K).containsZero())
        return Interval(0, -1);
    return Interval(0, 1e6); // TODO: Implement this function!
}

// We "find" the length of the last side.
Interval findGeneralizedPolygonMissingSide(const vector<double>& angles,
                                           vector<Interval> sides) {
    int k = ssize(angles);
    assert(ssize(sides) == k-1);
    //~ bool bad_situa = false;
    //~ for (Interval I: sides)
        //~ if (I.len() >= 0.3)
            //~ bad_situa = true;
    
    //~ if (bad_situa)
        //~ return Interval(-10, 10);
    
    vector<double> theta(k+1);
    theta[0] = 0;
    for (int i = 0; i < k; i++) theta[i+1] = theta[i] + angles[i] - M_PI;
    
    Quaternion<Interval> res(cos(theta[k]/2), -sin(theta[k]/2), 0, 0);
    for (int i = 0; i < k-1; i++) {
        Quaternion<Interval> q(0, 0, -sin(theta[i]), cos(theta[i]));
        //~ cout << sides[i] << endl;
        Interval Il = tan(0.5 * sides[i]);
        res *= Quaternion<Interval>(1) - Quaternion<Interval>(Il) * q;
    }
    // 1 + tan(len/2) * (0, 0, -sin(theta[k-1]), cos(theta[k-1]))
    
    //~ if (bad_situa and !res[1].containsZero()) {
        //~ cout << "WTF" << endl;
        //~ cout << res[1] << endl;
    //~ }
    
    if (!res[1].containsZero())
        return Interval(0, -1); // empty
    
    double x = -sin(theta[k-1]), y = cos(theta[k-1]);
    
    Interval I = subproblemInterval(x, y, res[0], res[2], res[3]);
    return atan(I) * 2;
}

#include "derived_once_nets.cpp"

vector<DerivedOnceNet> parseDerivedOnceNets(string in_file) {
    ifstream in(in_file);
    int cnt;
    in >> cnt;
    vector<DerivedOnceNet> res;
    for (int i = 0; i < cnt; i++) {
        DerivedOnceNet A = parseDerivedOnceNet(in);
        A.checkConsistency();
        res.push_back(A);
    }
    return res;
}


// Determines using only convex faces.
bool isDetermined(const DerivedOnceNet& A, vector<bool> determined) {
    auto is_face_determined = [&](Dart d) {
        if (!A.isConvexFace(d))
            return false;
        int missing = 0;
        for (Dart e: A.face(d)) {
            missing += !determined[e.id];
        }
        return 0 < missing and missing <= 3;
    };
    auto determine_face = [&](Dart d) {
        for (Dart e: A.face(d)) {
            determined[e.id] = 1;
            determined[A.rev(e).id] = 1;
        }
    };
    while (1) {
        bool changed = false;
        for (Dart d: A.darts()) {
            if (is_face_determined(d)) {
                changed = true;
                determine_face(d);
            }
        }
        if (!changed)
            break;
    }
    int missing = 0;
    for (Dart d: A.darts()) missing += A.isConvexFace(d) and !determined[d.id];
    return missing == 0;
}

vector<int> dartsToFix(const DerivedOnceNet& A) {
    int n = A.E();
    if (isDetermined(A, vector<bool>(2*n, false)))
        return {};
        
    for (int ans = 1; ans <= 6; ans++) {
        long long bb = 0;
        while (bb < (1ll<<n)) {
            if (__builtin_popcount(bb) == ans) {
                vector<bool> vec(2*n, false);
                for (int i = 0; i < n; i++) {
                    if (bb&(1ll<<i)) {
                        vec[2*i] = 1;
                        vec[2*i+1] = 1;
                    }
                }
                if (isDetermined(A, vec)) {
                    vector<int> res;
                    for (int i = 0; i < n; i++)
                        if (bb&(1ll<<i))
                            res.push_back(2*i);
                    return res;
                }
                bb += 1ll<<__builtin_ctzll(bb);
            } else bb++;
        }
    }
    assert(0);
    return {-1};
}

#include "spherical_geometry.cpp"

bool metrizeConvexContinent(const DerivedOnceNet& A, vector<Interval>& lengths) {
    auto print_status = [&]() {
        cout << "State: " << endl;
        for (int i = 0; i < ssize(lengths); i += 2)
            cout << A.dart(i) << " " << lengths[i] << endl;
        cout << endl;
    };
    
    auto compute_len = [&]() {
        double res = 0;
        for (int i = 0; i < 2*A.E(); i++) {
            if (lengths[i].empty())
                return -1.0; // impossible
            res += lengths[i].len();
        }
        return res;
    };
    double curr_len = compute_len();
    while (curr_len >= 1e-8) {
        for (auto f: A.faces()) {
            if (!A.isConvexFace(f))
                continue;
            
            vector<Interval*> sides;
            for (const Dart& d: f)
                sides.push_back(&lengths[d.id]);
            if (ssize(f) == 3)
                enforceTriangleConstraint(sides);
            if (ssize(f) == 4)
                enforceQuadrilateralConstraint(sides);    
            if (ssize(f) == 5) {
                enforcePentagonConstraint(sides);
            }
            
            for (const Dart& d: f) 
                lengths[A.rev(d).id].intersect(lengths[d.id]);
        }
        //~ print_status();
        double new_len = compute_len();
        if (new_len > 0.99*curr_len)
            break;
        curr_len = new_len;
    }
    
    if (curr_len < 0)
        return false;
    
    
    // Now, we check the non-convex faces.
    for (Dart d: A.darts()) {
        if (A.isConvexFace(d))
            continue;
        int q = ssize(A.face(d));
        vector<double> angles;
        vector<Interval> sides;
        for (int i = 0; i < q; i++) {
            Dart e = A.phi(d, i);
            angles.push_back(A.isConcaveAngle(e) ? 4*M_PI/3 : 2*M_PI/3);
            if (i < q-1)
                sides.push_back(lengths[e.id]);
        }
        Dart e = A.phi(d, q-1);
        Interval& I = lengths[e.id];
        I.intersect(findGeneralizedPolygonMissingSide(angles, sides));
        lengths[A.rev(e).id].intersect(I);
        if (I.empty())
            return false;
    }
    
    return true;
}

int boh_cnt = 0;

void metrizeNet(const DerivedOnceNet& A, 
                vector<int> fixed_darts, vector<Interval> fixed_lengths,
                vector<Interval> lengths) {
    auto print_status = [&]() {
        cout << "Fixed: " << endl;
        int q = ssize(fixed_darts);
        for (int i = 0; i < q; i++)
            cout << A.dart(fixed_darts[i]) << ": " 
                 << fixed_lengths[i] << "; ";
        cout << endl;
        cout << "State: " << endl;
        for (int i = 0; i < ssize(lengths); i += 2)
            cout << A.dart(i) << " " << lengths[i] << endl;
        cout << endl;
    };
    
    assert(ssize(fixed_darts) == ssize(fixed_lengths));
    int q = ssize(fixed_darts);
    
    for (int i = 0; i < q; i++)
        lengths[fixed_darts[i]].intersect(fixed_lengths[i]);
    
    for (Dart d: A.darts())
        lengths[A.rev(d).id].intersect(lengths[d.id]);
    
    bool works = metrizeConvexContinent(A, lengths);
    
    if (!works) 
        return;
    
    if (fixed_lengths[0].len() < 0.01) {
        boh_cnt++;
        cout << boh_cnt << endl;
        print_status();
        return; // Too small.
    }
    
    for (int bb = 0; bb < (1<<q); bb++) {
        vector<Interval> new_lengths(q);
        for (int i = 0; i < q; i++) {
            const Interval& I = fixed_lengths[i];
            double m = (I.l + I.r)/2;
            if (bb&(1<<i))
                new_lengths[i] = Interval(I.l, m);
            else
                new_lengths[i] = Interval(m, I.r);
        }
        metrizeNet(A, fixed_darts, new_lengths, lengths);
    }
}


int main() {
    //~ Interval I(1.81, 1.817);
    //~ Interval J(0.85, 0.97);
    //~ cout << pentagon2Constraint(I, J) << endl;
    
    //~ return 0;
    
    
    auto net = parseDerivedOnceNets("data/derived_once_computer.txt")[2];
    vector<int> dart_ids = dartsToFix(net);
    for (int id: dart_ids) {
        cout << net.dart(id) << endl;
    }
    
    vector<Interval> lengths(2*net.E(), Interval(0, M_PI));
    
    for (Dart d: net.darts()) {
        // TODO: Shall be proven
        Interval& I = lengths[d.id];
        I.intersect(Interval(0.0001, M_PI-0.0001));
        if (net.isConvexFace(d)) {
            I.intersect(Interval(0.01, M_PI-0.001));
            int q = ssize(net.face(d));
            assert(q != 2);
            if (q == 5)
                I.intersect(Interval(0, acos(-1.0/3)-0.0001));
        }
    }
    for (Dart d: net.darts())
        lengths[d.id].intersect(lengths[net.rev(d).id]);
    
    metrizeNet(net, 
            dart_ids, vector<Interval>(ssize(dart_ids), Interval(0, M_PI)),
            lengths);
}

