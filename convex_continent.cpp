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
    
    double len() const { return empty() ? 0.0 : (r-l); }
    bool contains(double x) const { return l <= x and x <= r; }
    bool containsZero() const { return l <= 0 and 0 <= r; } // TODO: Remove.
};

ostream& operator <<(ostream& out, const Interval& I) {
    out << "[" << I.l << ", " << I.r << "]";
    return out;
}

Interval hull(vector<double> vec) { 
    return Interval(*min_element(vec.begin(), vec.end()), 
                    *max_element(vec.begin(), vec.end()));
}

Interval intersection(Interval I, const Interval& J) {
    I.intersect(J);
    return I;
}

Interval inverse(const Interval I) {
    assert(!I.containsZero());
    return Interval(1.0 / I.r, 1.0 / I.l);
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
    
    // Using acos on [low, high] is not necessarily correct here, hence
    // we do the computation explicitly.
    
    if (low > 1 or high < -1)
        return Interval(0, -1); // empty
    
    double l = (high < 1) ?  acos(high) : 0;
    l -= get_theta(I.r);
    double r = (low > -1) ? acos(low) : M_PI;
    r -= get_theta(I.l);
    Interval K(l, r);
    
    const static double ellT = acos(-1.0/3);
    K.intersect(Interval(0, ellT));
    
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

Interval sin(const Interval& I) {
    if (I.empty())
        return Interval(0, -1);
    double sl = sin(I.l), sr = sin(I.r);
    
    double l = containsMod2pi(I, 3*M_PI/2) ? -1 : min(sl, sr);
    double r = containsMod2pi(I, M_PI/2) ? 1 : max(sl, sr);
    
    return Interval(l, r);
}

Interval cos(const Interval& I) {
    if (I.empty())
        return Interval(0, -1);
    double sl = cos(I.l), sr = cos(I.r);
    
    double l = containsMod2pi(I, M_PI) ? -1 : min(sl, sr);
    double r = containsMod2pi(I, 0) ? 1 : max(sl, sr);
    
    return Interval(l, r);
}

// We employ the identification R^3 \cong <i, j, k> \subseteq Quaternions.
// Then, any quaternion q with |q|=1 induces an isometry via v -> qv\bar q.
 
// Set v = (e_1, e_2).
// This function returns R(v, theta).
template<typename T> 
Quaternion<T> basicRotation(T theta) {
    T theta2 = 0.5 * theta;
    return Quaternion<T>(cos(theta2), sin(theta2), 0, 0);
}

// This function returns E(\rot_\theta v, len)
template<typename T> 
Quaternion<T> basicIsometry(T theta, T len) {
    T len2 = 0.5 * len;
    return Quaternion<T>(cos(len2), 0, 
            -sin(theta) * sin(len2), cos(theta) * sin(len2)); 
}

Interval findGeneralizedPolygonMissingSide2(const vector<double>& angles,
                                           vector<Interval> sides) {
    int k = ssize(angles);
    assert(ssize(sides) == k-1);
        
    vector<double> theta(k+1);
    theta[0] = 0;
    for (int i = 0; i < k; i++) theta[i+1] = theta[i] + angles[i] - M_PI;
    
    Quaternion<Interval> res(1);
    for (int i = k-2; i >= 0; i--) {
        Quaternion<Interval> isom = basicIsometry<Interval>(theta[i], sides[i]);
        res *= conj(isom);
        //~ Quaternion<Interval> q(0, 0, -sin(theta[i]), cos(theta[i]));
        //~ Interval angle = 0.5 * sides[i];
        //~ res *= Quaternion<Interval>(cos(angle)) - Quaternion<Interval>(sin(angle)) * q;
    }
    res *= basicRotation<Interval>(-theta[k]);
    // Up to sign, this is basicIsometry(theta[k-1], len) =
    // cos(len/2) + sin(len/2) * (0, 0, -sin(theta[k-1]), cos(theta[k-1]))
    
    if (!res[1].containsZero())
        return Interval(0, -1); // empty
    
    double x = -sin(theta[k-1]), y = cos(theta[k-1]);
    
    Interval I = res[0];
    Interval J(-2, 2);
    if (abs(x) > 0.01)
        J.intersect(res[2] * (1/x));
    else if (!res[2].containsZero())
        return Interval(0, -1);
    
    if (abs(y) > 0.01)
        J.intersect(res[3] * (1/y));
    else if (!res[3].containsZero())
        return Interval(0, -1);
             
    if (I.empty() or J.empty())
        return Interval(0, -1);
    
    if (I.len() > 0.9 or J.len() > 0.9)
        return Interval(0, 10);
        
    // (cos l/2, sin l/2) belongs to IxJ or -(IxJ).
    // With 0 <= l < pi.
    
    auto intersectWithCircleQuadrant = [](Interval X, Interval Y) {
        X.intersect(Interval(0, 1));
        Y.intersect(Interval(0, 1));
        if (X.empty() or Y.empty())
            return Interval(0, -1);
        return intersection(Interval(acos(X.r), acos(X.l)),
                            Interval(asin(Y.l), asin(Y.r)));
    };
    
    Interval inter1 = intersectWithCircleQuadrant(I, J);
    Interval inter2 = intersectWithCircleQuadrant((-1)*I, (-1)*J);
    
    assert(inter1.empty() or inter2.empty());
    Interval inter = inter1.empty() ? inter2 : inter1;
    
    return 2 * inter;
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

// Returns the smallest possible set of darts such that once their 
// length is chosen, the lengths of all the sides of convex faces
// can be uniquely determined through simple rules.
vector<Dart> findFreeParametersConvex(const DerivedOnceNet& A) {
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
                    vector<Dart> res;
                    for (int i = 0; i < n; i++)
                        if (bb&(1ll<<i))
                            res.push_back(A.dart(2*i));
                    return res;
                }
                bb += 1ll<<__builtin_ctzll(bb);
            } else bb++;
        }
    }
    assert(0);
    return {};
}

vector<Dart> findFreeParametersNonConvex(const DerivedOnceNet& A) {
    vector<Dart> res;
    for (auto f: A.faces()) {
        if (A.isConvexFace(f))
            continue;
        for (Dart d: f) {
            if (!A.isConvexFace(A.rev(d)))
                res.push_back(d);
        }
        break;
    }
    
    if (!res.empty())
        res.pop_back();

    return res;
}

vector<Dart> findFreeParameters(const DerivedOnceNet& A) {
    vector<Dart> convex = findFreeParametersConvex(A);
    vector<Dart> non_convex = findFreeParametersNonConvex(A);
    
    vector<Dart> joined = convex;
    for (Dart d: non_convex)
        joined.push_back(d);
    return joined;
}

#include "spherical_geometry.cpp"

// Given a sequence of operations, it returns the quaternion representing
// the isometry that maps (e_1, e_2) into the result of the operations 
// starting from (e_1, e_2).
//
// The operations can be:
//  - {'R', \theta}: Apply \rot_\theta.
//  - {'E', l}: Apply \exp_l.
template<typename T>
Quaternion<T> isometryFromOperations(vector<pair<char, T>> ops) {
    T sum_theta = 0;
    Quaternion<T> res = 1;
    for (auto op: ops) {
        if (op.first == 'R')
            sum_theta += op.second;
        if (op.first == 'E')
            res *= basicIsometry(sum_theta, op.second);
    }
    res *= basicRotation(sum_theta);
    
    return res;
}


template <typename T>
vector<vector<pair<char, T>>> pathToConcaveAngles(
        const DerivedOnceNet& A, const vector<Interval>& lengths) {
    // Step 0: Find the vertices with concave angles.
    vector<int> concave_angles;
    for (Dart d: A.darts())
        if (A.isConcaveAngle(d))
            concave_angles.push_back(d.to);
            
    if (concave_angles.empty())
        return {};
    
    
    // Step 1: Find the best "center".
    //         TODO: Now we just select arbitrarily one of the 
    //               concave angles.
    int v = concave_angles[0];
        
    // Step 2: Run Dijkstra from the chosen dart to the concave angles.
    Dart root;
    for (Dart d: A.darts())
        if (d.from == v)
            root = d;
    using pdd = pair<double, Dart>;
    vector<double> distance(2*A.E(), 1e9);
    priority_queue<pdd, vector<pdd>, greater<pdd>> q;
    
    distance[root.id] = 0;
    q.push({0, root});
    while (!q.empty()) {
        auto pp = q.top();
        q.pop();
        Dart d = pp.second;
        double dist = pp.first;
        if (dist > distance[d.id])
            continue;
        
        vector<pdd> adj = {{dist + 1e-10, A.next(d)}, 
                           {dist + lengths[d.id].len() + 1e-10, A.phi(d)}};
        for (auto [dist2, d2]: adj) {
            if (dist2 >= distance[d2.id]) 
                continue;
            distance[d2.id] = dist2;
            q.push({dist2, d2});
        }
    }
    
    // Step 3: Reconstruct the paths from the chosen dart to the edges
    //         going into the pin.
}

bool refineMetrization(const DerivedOnceNet& A, vector<Interval>& lengths) {
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
            
            assert(3 <= ssize(f) and ssize(f) <= 5);    
            
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
        double new_len = compute_len();
        if (new_len > 0.99*curr_len)
            break;
        curr_len = new_len;
    }
    
    if (curr_len < 0)
        return false;
    
    // Now, we check the non-convex faces.
    for (int _ = 0; _ < 2; _++) { // We run it twice.
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
            I.intersect(findGeneralizedPolygonMissingSide2(angles, sides));
            lengths[A.rev(e).id].intersect(I);
            if (I.empty())
                return false;
        }
    }
    
    return true;
}

int boh_cnt = 0;

bool searchApproximateMetrization(const DerivedOnceNet& A, 
                const vector<Dart>& fixed_darts, vector<Interval> fixed_lengths,
                vector<Interval> lengths) {
    auto print_status = [&]() {
        cout << "Fixed: " << endl;
        int q = ssize(fixed_darts);
        for (int i = 0; i < q; i++)
            cout << fixed_darts[i] << ": " 
                 << fixed_lengths[i] << "; ";
        cout << endl;
        //~ cout << "State: " << endl;
        //~ for (int i = 0; i < ssize(lengths); i += 2)
            //~ cout << A.dart(i) << " " << lengths[i] << endl;
        //~ cout << endl;
    };
    
    assert(ssize(fixed_darts) == ssize(fixed_lengths));
    int q = ssize(fixed_darts);
    
    for (int i = 0; i < q; i++)
        lengths[fixed_darts[i].id].intersect(fixed_lengths[i]);
    
    for (Dart d: A.darts())
        lengths[A.rev(d).id].intersect(lengths[d.id]);
    
    if (!refineMetrization(A, lengths))
        return false;
    
    if (fixed_lengths[0].len() < 0.001) 
        return true; // We found a valid approximate metrization.
    
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
        if (searchApproximateMetrization(A, fixed_darts, new_lengths, lengths))
            return true;
    }
    
    return false;
}

bool isApproximateMetrizationPossible(const DerivedOnceNet& A) {
    vector<Dart> fixed_darts = findFreeParameters(A);
    
    //~ cout << "The " << ssize(fixed_darts) << " fixed darts are: ";
    //~ for (Dart d: fixed_darts) {
        //~ cout << d << " ";
    //~ }
    //~ cout << endl;
    cout << ssize(fixed_darts) << " ";
    cout.flush();
    
    int cnt = 0;
    for (Dart d: A.darts())
        if (A.isConcaveAngle(d))
            cnt++;
    if (ssize(fixed_darts) >= 5 or (cnt == 4 and ssize(fixed_darts) >= 4)) { // slow and non-effective
        cout << "SKIPPED";
        return true;
    }
    
    vector<Interval> lengths(2*A.E(), Interval(0, M_PI));
    
    for (Dart d: A.darts()) {
        Interval& I = lengths[d.id];
        if (A.isConvexFace(d)) {
            I.intersect(Interval(0.02, M_PI-0.01));
            int q = ssize(A.face(d));
            assert(q != 2);
            if (q == 5)
                I.intersect(Interval(0, acos(-1.0/3)-0.0001));
        }
    }
    for (Dart d: A.darts())
        lengths[d.id].intersect(lengths[A.rev(d).id]);
        
    return  searchApproximateMetrization(A, fixed_darts, 
            vector<Interval>(ssize(fixed_darts), Interval(0, M_PI)),
            lengths);
}


int main() {
    // 1s with the tangent implementation
    // 2, 3, 5, 8, 12, 19, 31, 83
    
    ifstream in("data/derived_once_survivors.txt");
    int cnt_survivors;
    in >> cnt_survivors;
    vector<int> survivors(cnt_survivors);
    for (int i = 0; i < cnt_survivors; i++)
        in >> survivors[i];
    in.close();
    
    for (int id: survivors) {
        //~ if (id < 600)
            //~ continue;
        cout << id << " ";
        auto net = parseDerivedOnceNets("data/derived_once_computer.txt")[id];
        bool possible = isApproximateMetrizationPossible(net);
        if (possible) cout << "POSSIBLE";
        cout << endl;
    }
}


// Correct output.
//~ 2 3 
//~ 5 2 
//~ 7 4 SKIPPEDPOSSIBLE
//~ 8 3 
//~ 13 3 
//~ 19 3 POSSIBLE
//~ 31 2 
//~ 33 4 SKIPPEDPOSSIBLE
//~ 49 3 
//~ 54 3 
//~ 57 1 
//~ 62 3 POSSIBLE

