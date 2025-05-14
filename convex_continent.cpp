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

pair<vector<pair<pt, string>>, vector<GeodesicSegment>> 
        embedNet(const MetrizedNet& net) {
    const vector<Interval>& lengths = net.lengths;
    const DerivedOnceNet& A = net.A;
    
    int n = ssize(lengths);
    vector<GeodesicSegment> embedding(n);
    vector<bool> visited(n, false);
    function<void(Dart)> dfs = [&](Dart d) {
        GeodesicSegment& curr = embedding[d.id];
        double len = lengths[d.id].mid();
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

long long begin_time;
int current_net_id;


// Return values:
// 0 no
// 1 yes
// 2 maybe
int searchApproximateMetrization(MetrizedNet net, 
                                 const vector<Dart>& fixed_darts) {
    const DerivedOnceNet& A = net.A;
    
    // If the computation takes too long, we skip.
    if (get_time() - begin_time > 10*60ll * 1e9) 
        return 2;
    // If we reached the desired precision of the metrization, we 
    // return.
    if (totalLen(net.lengths) < 1e-5) {
        checkPinAngles(net);
        auto pp = embedNet(net);
        savePicture(pp.first, pp.second, "data/metrization/" + to_string(current_net_id) + ".txt");
        return 1;
    }
    
    auto print_status = [&]() {
        cout << "State: " << endl;
        for (int i = 0; i < ssize(net.lengths); i += 2)
            cout << A.dart(i) << " " << net.lengths[i] << endl;
        cout << endl;
    };
    
    
    for (Dart d: A.darts())
        net.lengths[A.rev(d).id].intersect(net.lengths[d.id]);
    
    if (!refineMetrizationConvexContinent(net))
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
    
    int q = ssize(fixed_darts);
    for (int bb = 0; bb < (1<<q); bb++) {
        vector<Interval> new_lengths = net.lengths;
        for (int i = 0; i < q; i++) {
            Interval& I = new_lengths[fixed_darts[i].id];
            double m = (I.l + I.r)/2;
            if (bb&(1<<i))
                I = Interval(I.l, m);
            else
                I = Interval(m, I.r);
        }
        auto new_net = net;
        new_net.lengths = new_lengths;
        // This is the only external overhead necessary for the
        // caching mechanis of minimal paths to work.
        // It could be part of a custom copy constructor, but it
        // would be even more annoying.
        new_net.incrementVersion();
        int outcome = searchApproximateMetrization(new_net, fixed_darts);
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
    vector<Dart> fixed_darts = findFreeParameters(A);
    
    cout << "# of parameters = " << ssize(fixed_darts) << ", ";
    cout.flush();
    
    return searchApproximateMetrization(MetrizedNet(A), fixed_darts);
}

int main() {
    ifstream in("data/derived_once_survivors.txt");
    int cnt_survivors;
    in >> cnt_survivors;
    vector<int> survivors(cnt_survivors);
    for (int i = 0; i < cnt_survivors; i++)
        in >> survivors[i];
    in.close();
    
    
    vector<int> new_survivors;
    for (int id: survivors) {
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

