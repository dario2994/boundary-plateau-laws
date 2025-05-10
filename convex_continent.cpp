#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "interval_arithmetic.cpp"
#include "quaternions.cpp"
#include "derived_once_nets.cpp"
#include "metrization_constraints.cpp"


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
        for (Dart e: A.face(d))
            missing += !determined[e.id];
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
    for (Dart d: A.darts()) 
        missing += A.isConvexFace(d) and !determined[d.id];
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

// Returns the time elapsed in nanoseconds from 1 January 1970, at 00:00:00.
long long get_time() {
    return chrono::duration_cast<chrono::nanoseconds>(
        chrono::steady_clock::now().time_since_epoch())
        .count();
}
long long begin_time;

bool searchApproximateMetrization(const DerivedOnceNet& A, 
                                  const vector<Dart>& fixed_darts, 
                                  vector<Interval> lengths) {
    // If the computation takes too long, we skip.
    if (get_time() - begin_time > 7*60*60ll * 1e9) {
        cout << "skipped, ";
        return true;
    }
    // If we reached the desired precision of the metrization, we 
    // return.
    if (totalLen(lengths) < 1e-5) 
        return true;
    
    auto print_status = [&]() {
        cout << "State: " << endl;
        for (int i = 0; i < ssize(lengths); i += 2)
            cout << A.dart(i) << " " << lengths[i] << endl;
        cout << endl;
    };
    
    
    //~ if (lengths[fixed_darts[0].id].len() > 0.1) {
        //~ cout << int((get_time() - begin_time) / 1e9) << endl;
        //~ print_status();
    //~ }
    
    
    for (Dart d: A.darts())
        lengths[A.rev(d).id].intersect(lengths[d.id]);
    
    if (!refineMetrizationConvexContinent(A, lengths))
        return false;
    
    if (!refineMetrizationNonConvexFaces(A, lengths))
        return false;
    
    if (!canPinsBeAntipodal(A, lengths)) 
        return false;
    
    if (!checkPinAngles(A, lengths))
        return false;
        
    if (!checkNonIntersection(A, lengths))
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
    ifstream in("data/derived_once_survivors.txt");
    int cnt_survivors;
    in >> cnt_survivors;
    vector<int> survivors(cnt_survivors);
    for (int i = 0; i < cnt_survivors; i++)
        in >> survivors[i];
    in.close();
    
    for (int id: {113}) {//{19,62,113,121,247,278,292,312,330,331}) {
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

