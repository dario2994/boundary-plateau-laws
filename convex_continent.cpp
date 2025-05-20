#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "interval_arithmetic.cpp"
#include "quaternions.cpp"
#include "derived_once_nets.cpp"
#include "metrization_constraints.cpp"
#include "metrized_nets.cpp"


#include "spherical_geometry.cpp"

void savePicture(vector<pair<pt, string>> vertices, 
                 vector<GeodesicSegment> edges, string out_file) {
    ofstream out(out_file);
    out << ssize(vertices) << " " << ssize(edges) << endl;
    for (auto pp: vertices) {
        pp.first.serialize(out);
        out << " " << pp.second << endl;
    }
    for (const GeodesicSegment& side: edges) {
        side.serialize(out);
        out << endl;
    }
    out.close();
}

// TODO: Translate this in the quaternion language.
pair<vector<pair<pt, string>>, vector<GeodesicSegment>> 
        embedNet(const MetrizedNet& net) {
    const DerivedOnceNet& A = net.A;
    
    int n = 2*A.E();
    vector<GeodesicSegment> embedding(n);
    vector<bool> visited(n, false);
    function<void(Dart)> dfs = [&](Dart d) {
        GeodesicSegment& curr = embedding[d.id];
        double len = net.length(d).mid();
        curr.len = len;
        visited[d.id] = true;
        // Assuming that rev != next. True if there are no self-loops.
        if (!visited[A.next(d).id]) {
            embedding[A.next(d).id] = GeodesicSegment(curr.rot(A.nextAngle(d)), 0);
            dfs(A.next(d));
        }
        if (!visited[A.rev(d).id]) {
            embedding[A.rev(d).id] = GeodesicSegment(curr.exp(len).rot(M_PI), 0);
            dfs(A.rev(d));
        }
    };
    embedding[0] = GeodesicSegment(TangentVector(), 0);
    dfs(A.dart(0));
    
    vector<pair<pt, string>> vertices(A.V());
    for (Dart d: A.darts()) {
        vertices[d.from] = {embedding[d.id].P, to_string(d.from)};
    }
    
    return {vertices, embedding};
}



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
// It returns the length of the longest implication chain that is
// necessary to determine a length.
// The length of an implication chain depends on the "pentagonal rules" 
// applied.
int lengthsImplications(const DerivedOnceNet& A, 
                        const vector<Dart>& fixed_darts) {
    static const int INF = 1'000'000;
    vector<int> vec(2*A.E(), INF);
    for (const Dart& d: fixed_darts) {
        vec[d.id] = 1;
        vec[A.rev(d).id] = 1;
    }
    
    bool updated = true;
    auto determine_edge = [&](const Dart& d, int val) {
        for (int id: {d.id, A.rev(d).id}) {
            if (vec[id] > val) {
                vec[id] = val;
                updated = true;
            }
        }
    };
    
    while (updated) {
        updated = false;
        for (const Dart& d: A.darts()) {
            if (!A.isConvexFace(d))
                continue;
            const vector<Dart>& f = A.face(d);
            if (ssize(f) <= 3) for (const Dart& e: f) determine_edge(e, 0);
            
            if (ssize(f) == 4) {
                int m = INF;
                for (const Dart& e: f) m = min(m, vec[e.id]);
                for (const Dart& e: f) determine_edge(e, m);
            }
            
            if (ssize(f) == 5) {
                vector<int> x;
                for (int i = 0; i < 5; i++)
                    x.push_back(vec[f[i].id]);
                determine_edge(f[3], max(x[0], x[1]) + 1);
                determine_edge(f[4], max(x[0], x[2]) + 1);
                determine_edge(f[3], max(x[0], x[2]) + 1);
                
                //~ determine_edge(f[3], x[0] + x[1] + 1);
                //~ determine_edge(f[4], x[0] + x[2] + 1);
                //~ determine_edge(f[3], x[0] + x[2] + 1);
            }
        }
    }
    
    int outcome = 0;
    for (const Dart& d: A.darts()) {
        if (!A.isConvexFace(d) and !A.isConvexFace(A.rev(d)))
            continue;
        outcome = max(outcome, vec[d.id]);
    }
    
    if (outcome == INF) 
        return -1;
    
    return outcome;
}

// Returns the smallest possible set of darts such that once their 
// length is chosen, the lengths of all the sides of convex faces
// can be uniquely determined through simple rules.
// Among the valid possible sets, it returns one that determines
// all the lengths in a "numerically robust" way.
vector<Dart> findFreeParametersConvex(const DerivedOnceNet& A) {
    int n = A.E();
    if (lengthsImplications(A, {}) >= 0)
        return {};
        
    for (int ans = 1; ans <= 6; ans++) {
        int least_outcome = 1e9;
        vector<Dart> res;
        long long bb = 0;
        while (bb < (1ll<<n)) {
            if (__builtin_popcountll(bb) == ans) {
                vector<Dart> fixed_darts;
                for (int i = 0; i < n; i++)
                    if (bb&(1ll<<i)) 
                        fixed_darts.push_back(A.dart(2*i));
                
                int outcome = lengthsImplications(A, fixed_darts);
                if (outcome >= 0 and outcome < least_outcome) {
                    least_outcome = outcome;
                    res = fixed_darts;
                }
                bb += 1ll<<__builtin_ctzll(bb);
            } else bb++;
        }
        if (!res.empty()) {
            //~ cout << least_outcome << endl;
            return res;
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
    
    // This is a way to guarantee that if res consists of three adjacent
    // edges, the central one is not chosen.
    // This is heuristically convenient.
    //~ sort(res.begin(), res.end(), [&](Dart d) {
        //~ return -(A.degree(d.from) + A.degree(d.to));
    //~ });
    
    if (!res.empty())
        res.pop_back();

    return res;
}

// TODO: Keeping these globally is horrible.
long long begin_time;
int current_net_id;

// Return values:
// 0 no
// 1 yes
// 2 maybe
int searchApproximateMetrization(MetrizedNet net, 
                                 const vector<Dart>& convex_fixed_darts,
                                 const vector<Dart>& non_convex_fixed_darts) {
    const DerivedOnceNet& A = net.A;
    
    auto print_status = [&]() {
        cout << "State: " << net.totalError() << endl;
        for (int i = 0; i < 2*A.E(); i += 2) {
            Dart d = A.dart(i);
            cout << d << " " << net.length(d) << "\t" << 
                    net.length(d).len() << endl;
        }
        cout << endl;
    };
    
    //~ if (net.totalError() < 10) {
        //~ if (get_time() % 30000 == 0) {
            //~ print_status();
            //~ cout << endl;
        //~ }
        //~ boring_counter++;
        //~ return 0;
    //~ }
    
    // If the computation takes too long, we skip.
    if (get_time() - begin_time > 10ll * 1e9) 
        return 2;
    // If we reached the desired precision of the metrization, we 
    // return.
    if (net.totalError() < 1e-5) {
        
        // The following two checks are run only here because:
        //  - They handle rare pathological situations.
        //  - They are almost always implied by other checks.
        //  - They are slow.
        
        if (!checkPinsSpecial(net))
            return 0;
            
        //~ if (!checkFacesNotIdentical(net))
            //~ return 0;

        auto pp = embedNet(net);
        savePicture(pp.first, pp.second, "data/metrization/" + to_string(current_net_id) + ".txt");
        //~ for (const Dart& d: A.darts())
            //~ if (net.length(d).mid() < 0.1)
                //~ cout << d << ": " << net.length(d).mid() << "; ";
        return 1;
    }
    
    // TODO: We shall memorize the outcome of 
    //  1. checkPinAngles (actually this is already very fast)
    //  2. checkNonIntersection
    //  3. checkPinOutsideConvexFaces
    
    if (!refineMetrizationConvexContinent(net, convex_fixed_darts))
        return 0;
    
    if (!refineMetrizationNonConvexFaces(net))
        return 0;
    
    if (!canPinsBeAntipodal(net)) 
        return 0;
    
    if (!checkPinAngles(net))
        return 0;
        
    if (!checkNonIntersection(net))
        return 0;
        
    if (!checkPinOutsideConvexFaces(net))
        return 0;
        
    //~ if (!checkLengthLowerbound(net))
        //~ return 0;
        
    double max_error_boundary_continent = 0;
    for (const Dart& d: A.darts()) {
        if (A.isConvexFace(d) == A.isConvexFace(A.rev(d)))
            continue;
        double l = net.length(d).len();
        max_error_boundary_continent = max(max_error_boundary_continent, l);
    }
    
    vector<Dart> fixed_darts = convex_fixed_darts;
    double max_non_convex = 0;
    for (const Dart& d: non_convex_fixed_darts) {
        double l = net.length(d).len();
        max_non_convex = max(max_non_convex, l);
    }
    if (max_non_convex > max_error_boundary_continent) {
        for (const Dart& d: non_convex_fixed_darts)
            fixed_darts.push_back(d);
    }
        
    int q = ssize(fixed_darts);
    for (int bb = 0; bb < (1<<q); bb++) {
        auto new_net = net;
        for (int i = 0; i < q; i++) {
            Interval& I = new_net.length(fixed_darts[i]);
            double m = (I.l + I.r)/2;
            if (bb&(1<<i))
                I = Interval(I.l, m);
            else
                I = Interval(m, I.r);
        }

        // This is the only external overhead necessary for 
        // caching minimal paths.
        // It could be part of a custom copy constructor, but it
        // would be even more annoying.
        new_net.incrementVersion();
        int outcome = searchApproximateMetrization(new_net, 
            convex_fixed_darts, non_convex_fixed_darts);
        if (outcome > 0)
            return outcome;
    }
    
    return 0;
}

// Return values:
// 0 no
// 1 yes
// 2 maybe
int isApproximateMetrizationPossible(const DerivedOnceNet& A) {
    vector<Dart> convex = findFreeParametersConvex(A);
    vector<Dart> non_convex = findFreeParametersNonConvex(A);
    
    cout << "# of parameters = " << ssize(convex) << " " 
                                 << ssize(non_convex) << ", ";
    
    //~ for (const Dart& d: convex)
        //~ cout << d << " ";
    cout.flush();
    
    return searchApproximateMetrization(MetrizedNet(A), convex, non_convex);
}

int main() {
    ifstream in("data/derived_once_survivors.txt");
    int cnt_survivors;
    in >> cnt_survivors;
    vector<int> survivors(cnt_survivors);
    for (int i = 0; i < cnt_survivors; i++)
        in >> survivors[i];
    in.close();
    
    // Doubled: 2202, 2287
    
    vector<int> new_survivors;
    for (int id: {278, 603, 746, 948, 1017, 1027, 1046, 1120, 1598, 1659, 1680, 1839, 2053, 2474, 2513}) {
        cout << "net id = " << setw(4) << id << ":\t";
        current_net_id = id;
        auto A = parseDerivedOnceNets("data/derived_once_computer.txt")[id];
        
        begin_time = get_time();
        int outcome = isApproximateMetrizationPossible(A);
        double duration_s = (get_time() - begin_time) / 1e9;
        
        cout.precision(2);
        cout << "(time = " << fixed << setw(6) << duration_s << ") ";
        if (outcome == 1) cout << "POSSIBLE";
        if (outcome == 2) cout << "SKIPPED";
        if (outcome != 0)
            new_survivors.push_back(id);
        cout << endl;
    }
    
    StopWatch::printState();
    
    return 0;
    
    //~ cout << cnt1 << " " << cnt2 << " " << cnt3 << " " << cnt4 << endl;
    
    cout << "New survivors num = " << new_survivors.size() << endl;
    for (int id: new_survivors)
        cout << id << " ";
    cout << endl;
    
    return 0;
    
    ofstream out("data/derived_once_survivors.txt");
    out << new_survivors.size() << endl;
    for (int id: new_survivors)
        out << id << " ";
    out << endl;
    out.close();
}

// Useful regexes
//~ net id = (.+)\:\t# of parameters = ., \(time = .+\) \n
//~ net id = +([0-9]+)\:\t# of parameters = ., \(time = .+\)
//~ (POSSIBLE|SKIPPED)\n
//~ ([0-9]+ SKIPPED|POSSIBLE)\n

