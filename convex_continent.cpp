#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "interval_arithmetic.cpp"
#include "quaternions.cpp"


// Enforces the constraint that the three numbers in sides correspond
// to the three sides of a spherical triangle with angles 120.
void enforceTriangleConstraint(vector<Interval*> sides) {
    assert(ssize(sides) == 3);
    
    static const double L = acos(-1.0/3);
    static const Interval I(L-1e-9, L+1e-9); 
    for (int i = 0; i < 3; i++) {
        sides[i]->intersect(I);
    }
}

// Utility function to enforce the constraint on the sides of a 
// quadrilateral.
// 
// I is strictly inside [0, \pi].
// One side of a quadrilateral belongs to I.
// The other side belongs to the interval returned by this function.
Interval quadrilateralConstraint(Interval I) {
    if (I.empty())
        return Interval(0, -1);
    const static double eps = 1e-10;
    static auto f = [](double x) { return 2*atan(0.5 / tan(x/2)); };
    return Interval(f(I.r) - eps, f(I.l) + eps);
}

// Enforces the constraint that the four numbers in sides correspond
// to the four sides of a spherical triangle with angles 120.
//
// Assumption: all sides are subsets of [0, pi].
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


// Utility function to enforce the constraint on the sides of a 
// pentagon.
// 
// I, J are strictly inside [0, \ell_T].
// Let ABCDE be a spherical pentagon with angles 120.
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

// Utility function to enforce the constraint on the sides of a 
// pentagon.
// 
// I, J are strictly inside [0, \ell_T].
// Let ABCDE be a spherical pentagon with angles 120.
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


// Enforces the constraint that the five numbers in sides correspond
// to the five sides of a spherical triangle with angles 120.
//
// Assumption: all sides are subsets of [0, \ell_T].
void enforcePentagonConstraint(vector<Interval*> sides) {
    assert(ssize(sides) == 5);
    
    auto s = [&](int id) { return sides[id % 5]; };
        
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
        double new_len = compute_len();
        if (new_len > 0.99*curr_len)
            break;
        curr_len = new_len;
    }
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
        
    // We shall return the interval of 0 <= l < pi such that:
    // (cos l/2, sin l/2) belongs to IxJ or -(IxJ).
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

struct PathStep {
    char type;
    double rotation;
    Dart d;
    
    PathStep(double _rotation): type('R'), rotation(_rotation) {}
    PathStep(Dart _d): type('E'), d(_d) {}
    
    double cost(const vector<Interval>& lengths) {
        return 1e-10 + ((type == 'R') ? 0.0 : lengths[d.id].len());
    }
};

// Given a sequence of operations (i.e., a path), it returns the 
// quaternion representing the isometry that maps (e_1, e_2) into the 
// result of the operations starting from (e_1, e_2).
//
// The operations can be:
//  - {'R', \theta}: Apply \rot_\theta.
//  - {'E', l}: Apply \exp_l.
Quaternion<Interval> isometryFromPath(const vector<PathStep>& path, 
                                      const vector<Interval>& lengths) {
    double sum_theta = 0;
    Quaternion<Interval> res(1);
    for (auto op: path) {
        if (op.type == 'R')
            sum_theta += op.rotation;
        if (op.type == 'E')
            res *= basicIsometry<Interval>(sum_theta, lengths[op.d.id]);
    }
    res *= basicRotation<Interval>(sum_theta);
    
    return res;
}

vector<vector<PathStep>> pathToConcaveAngles(
        const DerivedOnceNet& A, const vector<Interval>& lengths) {
    // Step 0: Find the vertices with concave angles.
    vector<int> concave_angles;
    for (Dart d: A.darts())
        if (A.isConcaveAngle(d))
            concave_angles.push_back(d.to);
            
    if (concave_angles.empty())
        return {};
        
    //~ if (ssize(concave_angles) <= 2)
        //~ return {}; // We save some time.
    
    
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
    vector<char> father(2*A.E(), '-'); // n = next, r = rev
    priority_queue<pdd, vector<pdd>, greater<pdd>> q;
    
    distance[root.id] = 0;
    q.push(pdd{0.0, root});
    while (!q.empty()) {
        auto pp = q.top();
        q.pop();
        Dart d = pp.second;
        double dist = pp.first;
        if (dist > distance[d.id])
            continue;
        
        vector<pdd> adj = {{dist + 1e-10, A.next(d)}, 
                           {dist + lengths[d.id].len() + 1e-10, A.rev(d)}};
        for (auto [dist2, d2]: adj) {
            if (dist2 >= distance[d2.id]) 
                continue;
            distance[d2.id] = dist2;
            father[d2.id] = (d.from == d2.from) ? 'n' : 'r';
            q.push({dist2, d2});
        }
    }
    assert(father[root.id] == '-');
    
    // Step 3: Reconstruct the paths from the chosen dart to the edges
    //         going into the pin.
    vector<vector<PathStep>> res;
    for (int goal: concave_angles) {
        double cost = 1e9;
        Dart last_dart;
        for (Dart d: A.darts()) {
            if (d.from == goal and cost > distance[d.id]) {
                cost = distance[d.id];
                last_dart = d;
            }
        }
        Dart curr = last_dart;
        vector<PathStep> path;
        
        // To compress consecutive rotations in a single one.
        auto append_to_path = [&](PathStep op) {
            if (!path.empty() and path.back().type == 'R' and op.type == 'R')
                path.back().rotation += op.rotation;
            else
                path.push_back(op);
        };
        
        // Rotation to go in the direction of the pin.
        double angle_to_pin = A.isConcaveAngle(A.rev(curr)) ? 2*M_PI/3 : M_PI/3;
        append_to_path(angle_to_pin);
        // We go from last_dart to root.
        while (father[curr.id] != '-') {
            if (father[curr.id] == 'n') { // Rotation by 120 or 240.
                double angle = A.isConcaveAngle(A.rev(A.prev(curr))) ? -2*M_PI/3 : 2*M_PI/3;
                append_to_path(angle);
                curr = A.prev(curr);
            } else if (father[curr.id] == 'r') { // Exponential + Rotation by 180.
                append_to_path(M_PI);
                curr = A.rev(curr);
                append_to_path(curr);
            } else assert(0);
        }
        reverse(path.begin(), path.end());
        res.push_back(path);
    }
    
    return res;
}

template <typename T>
T det(const array<T, 3>& A, const array<T, 3>& B, const array<T, 3>& C) {
    return A[0]*(B[1]*C[2] - B[2]*C[1])
         + A[1]*(B[2]*C[0] - B[0]*C[2])
         + A[2]*(B[0]*C[1] - B[1]*C[0]);
}

// Determines, approximately, whether there exists a plane that
// intersects all the given cubes.
// If it returns false, than it cannot. If it returns yes, it could
// be a false positive. 
//
// Implementing a precise check (without false positives) would most
// likely require some linear programming. It would be demanding to
// implement and possibly very slow.
bool areContainedInPlane(vector<array<Interval, 3>> cubes) {
    int q = ssize(cubes);
    for (int i = 0; i < q; i++) 
        for (int j = 0; j < q; j++)
            for (int k = j+1; k < q; k++) {
                if (i == j or i == k)
                    continue;
                if (!det(cubes[i], cubes[j], cubes[k]).containsZero())
                    return false;
            }
    return true;
}

bool canPinsBeAntipodal(const DerivedOnceNet& A, 
                        const vector<Interval>& lengths) {
    vector<vector<PathStep>> paths = pathToConcaveAngles(A, lengths);
    vector<Quaternion<Interval>> isometries;
    for (auto path: paths)
        isometries.push_back(isometryFromPath(path, lengths));
    
    vector<array<Interval, 3>> equators;
    for (auto q: isometries) {
        Quaternion<Interval> v(0, 0, 0, 1);
        Quaternion<Interval> equator = q * v * conj(q);
        equators.push_back({equator[1], equator[2], equator[3]});
    }
    
    return areContainedInPlane(equators);
}

// Returns the time elapsed in nanoseconds from 1 January 1970, at 00:00:00.
long long get_time() {
    return chrono::duration_cast<chrono::nanoseconds>(
        chrono::steady_clock::now().time_since_epoch())
        .count();
}
long long begin_time;


bool refineMetrizationConvexContinent(const DerivedOnceNet& A, 
                                      vector<Interval>& lengths) {
    while (totalLen(lengths) >= 1e-8) {
        double len_before = totalLen(lengths);
        for (const vector<Dart>& f: A.faces()) {
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
            if (ssize(f) == 5) 
                enforcePentagonConstraint(sides);
            
            for (const Dart& d: f) 
                lengths[A.rev(d).id].intersect(lengths[d.id]);
        }
        if (totalLen(lengths) > 0.99 * len_before)
            break;
    }
    
    return totalLen(lengths) >= 0;
}

bool refineMetrizationNonConvexFaces(const DerivedOnceNet& A, 
                                     vector<Interval>& lengths) {
    for (int _ = 0; _ < 2; _++) { // We run it twice without a reason.
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


bool searchApproximateMetrization(const DerivedOnceNet& A, 
                                  const vector<Dart>& fixed_darts, 
                                  vector<Interval> lengths) {
    // If the computation takes too long, we skip.
    if (get_time() - begin_time > 5*60ll * 1e9) {
        cout << "skipped, ";
        return true;
    }
    
    // If we reached the desired precision of the metrization, we 
    // return.
    if (totalLen(lengths) < 1e-3) 
        return true;
    
    auto print_status = [&]() {
        cout << "State: " << endl;
        for (int i = 0; i < ssize(lengths); i += 2)
            cout << A.dart(i) << " " << lengths[i] << endl;
        cout << endl;
    };
    
    
    for (Dart d: A.darts())
        lengths[A.rev(d).id].intersect(lengths[d.id]);
    
    if (!refineMetrizationConvexContinent(A, lengths))
        return false;
    
    if (!refineMetrizationNonConvexFaces(A, lengths))
        return false;
    
    if (!canPinsBeAntipodal(A, lengths)) 
        return false;
    
    int q = ssize(fixed_darts);
    for (int bb = 0; bb < (1<<q); bb++) {
        vector<Interval> new_lengths = lengths;
        for (int i = 0; i < q; i++) {
            Interval& I = new_lengths[fixed_darts[i].id];
            double m = (I.l + I.r)/2;
            if (bb&(1<<i))
                I = Interval(I.l, m);
            else
                I = Interval(m, I.r);
        }
        if (searchApproximateMetrization(A, fixed_darts, new_lengths))
            return true;
    }
    
    return false;
}

bool isApproximateMetrizationPossible(const DerivedOnceNet& A) {
    vector<Dart> fixed_darts = findFreeParameters(A);
    
    cout << "# of parameters = " << ssize(fixed_darts) << ", ";
    cout.flush();
    
    int cnt = 0;
    for (Dart d: A.darts())
        if (A.isConcaveAngle(d))
            cnt++;
    
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
        
    return  searchApproximateMetrization(A, fixed_darts, lengths);
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
        cout << "net id = " << id << ", ";
        auto net = parseDerivedOnceNets("data/derived_once_computer.txt")[id];
        
        begin_time = get_time();
        bool possible = isApproximateMetrizationPossible(net);
        double duration_s = (get_time() - begin_time) / 1e9;
        
        if (possible) cout << "POSSIBLE";
        cout << " (time = " << duration_s << ")";
        cout << endl;
    }
}

// On 1351: It detects possible (correctly) in 1m43s.
//          We shall speed it up.

