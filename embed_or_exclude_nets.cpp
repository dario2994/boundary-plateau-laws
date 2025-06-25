#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "interval_arithmetic.cpp"
#include "quaternions.cpp"
#include "derived_once_nets.cpp"
#include "metrization_constraints.cpp"
#include "metrized_nets.cpp"

typedef array<double, 3> pt;
ostream& operator<<(ostream& out, pt p) {
    out << p[0] << " " << p[1] << " " << p[2];
    return out;
}

void savePicture(vector<pair<pt, string>> vertices,
                 vector<tuple<pt, pt, double>> edges, string out_file) {
    ofstream out(out_file);
    out << ssize(vertices) << " " << ssize(edges) << endl;
    for (auto pp: vertices) {
        out << pp.first << " " << pp.second << endl;
    }
    for (const auto& side: edges) {
        auto [p, v, len] = side;
        out << p << " " << v << " " << len << endl;
    }
    out.close();
}

pair<vector<pair<pt, string>>, vector<tuple<pt, pt, double>>>
        embedNet(const MetrizedNet& net) {
    const DerivedOnceNet& A = net.A;

    vector<pair<pt, string>> vertex_embedding(A.V());
    vector<tuple<pt, pt, double>> dart_embedding(2 * A.E());
    Dart root = A.dart(0);
    for (Dart d: A.darts()) {
        ApproximateIsometry q = net.getIsometry(root, d);
        vertex_embedding[d.from] = {
            make_deterministic(q({1, 0, 0})),
            to_string(d.from)
        };
        dart_embedding[d.id] = {
            make_deterministic(q({1, 0, 0})),
            make_deterministic(q({0, 1, 0})),
            net.length(d).mid()
        };
    }

    return {vertex_embedding, dart_embedding};
}


// Determines the lengths of sides using only the constraints given by
// convex faces. It assumes that the lengths of fixed_darts are known
// a priori.
//
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
        if (!res.empty())
            return res;
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


// Takes an approximately metrized once-derived net and returns:
//  0: No refinement of the given approximate metrization can be realized by
//     a stable cone.
//  1: The current approximate metrization is precise enough (error 1e-5) and
//     the metrized net has antipodal pins, correct angles, and no side
//     extremely short.
//  INDETERMINATE: The program stopped because the execution was taking
//                 too long.
pair<int, optional<MetrizedNet>> searchApproximateMetrization(
        MetrizedNet net,
        const vector<Dart>& convex_fixed_darts,
        const vector<Dart>& non_convex_fixed_darts,
        long long begin_time,
        bool run_pathological_checks,
        const optional<MetrizedNet>& forbidden_metrization) {
    const DerivedOnceNet& A = net.A;

    // If the computation takes too long, we skip.
    // The slowest net, which is 72, requires 38 hours.
    static const long long one_second = 1e9;
    static const long long one_minute = 60*one_second;
    [[maybe_unused]] static const long long one_hour = 60*one_minute;
    if (get_time() - begin_time > 10*one_minute)
        return {INDETERMINATE, nullopt};

    if (net.totalError() < -0.5)
        return {0, nullopt};

    if (forbidden_metrization.has_value()) {
        bool is_refinement = true;
        for (const Dart& d: A.darts()) {
            is_refinement &= isSubset(net.length(d),
                                      forbidden_metrization->length(d));
        }
        if (is_refinement)
            return {0, nullopt};
    }

    if (!refineMetrizationConvexContinent(net, convex_fixed_darts))
        return {0, nullopt};

    if (!refineMetrizationNonConvexFaces(net))
        return {0, nullopt};

    if (!canPinsBeAntipodal(net))
        return {0, nullopt};

    if (!checkPinAngles(net))
        return {0, nullopt};

    if (!checkNonIntersection(net))
        return {0, nullopt};

    if (!checkPinOutsideConvexFaces(net))
        return {0, nullopt};

    if (!checkLengthLowerbound(net))
        return {0, nullopt};

    bool run_now_pathological_checks = net.totalError() < 1e-2
                                       and run_pathological_checks;
    if (net.totalError() < 1e-5)
        run_now_pathological_checks = true;

    if (run_now_pathological_checks) {
        // The following checks are run only here because:
        //  - They handle rare pathological situations.
        //  - They are almost always implied by other checks.
        //  - They are slow.
        if (!checkVerticesOutsideConvexFaces(net))
            return {0, nullopt};

        if (!checkFacesNotIdentical(net))
            return {0, nullopt};

        if (!checkPinIsVertex(net))
            return {0, nullopt};

        if (!checkPinsSpecial(net))
            return {0, nullopt};

        if (!checkFacesAreBoundaries(net))
            return {0, nullopt};
    }

    // If we reached the desired precision of the metrization, we
    // end the search.
    if (net.totalError() < 1e-5)
        return {1, net};


    // Computes the maximum error on the lengths at the boundary of
    // the continent made of convex faces.
    double max_error_boundary_continent = 0;
    for (const Dart& d: A.darts()) {
        if (A.isConvexFace(d) == A.isConvexFace(A.rev(d)))
            continue;
        double l = net.length(d).len();
        max_error_boundary_continent = max(max_error_boundary_continent, l);
    }

    // Computes the maximum error among the non_convex_fixed_darts.
    double max_non_convex = 0;
    for (const Dart& d: non_convex_fixed_darts) {
        double l = net.length(d).len();
        max_non_convex = max(max_non_convex, l);
    }

    // Decide whether we want to refine (i.e., half) the precision
    // only on convex_fixed_darts or also on non_convex_fixed_darts.
    vector<Dart> fixed_darts = convex_fixed_darts;
    if (max_non_convex > max_error_boundary_continent) {
        for (const Dart& d: non_convex_fixed_darts)
            fixed_darts.push_back(d);
    }

    // For all possible halvings of the confidence interval for the lengths
    // of fixed_darts, we call recursively the search.
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
        auto outcome = searchApproximateMetrization(new_net,
            convex_fixed_darts, non_convex_fixed_darts, begin_time, run_pathological_checks, forbidden_metrization);
        if (outcome.first != 0)
            return outcome;
    }

    return {0, nullopt};
}

// Takes a once-derived net and returns:
//  0: The net cannot be realized by a stable cone.
//  1: The net admits an approximate embedding (up to precision 1e-5) with
//     antipodal pins, correct angles, and no side extremely short.
//  INDETERMINATE: The program was not able to show either of the previous
//                 two possibilities.
pair<int, optional<MetrizedNet>> isApproximateMetrizationPossible(
        const DerivedOnceNet& A, bool run_pathological_checks, const optional<MetrizedNet>& forbidden_metrization) {
    vector<Dart> convex = findFreeParametersConvex(A);
    vector<Dart> non_convex = findFreeParametersNonConvex(A);

    clog << "# of parameters = " << ssize(convex) << " "
                                 << ssize(non_convex) << ", ";
    clog.flush();

    auto net = MetrizedNet(A);
    if (!refineMetrizationConvexContinent(net, A.darts()))
        return {0, nullopt};

    return searchApproximateMetrization(
            net, convex, non_convex, get_time(), run_pathological_checks, forbidden_metrization);
}

void saveMetrizedNet(const MetrizedNet& net, string picture, string data) {
    auto pp = embedNet(net);
    savePicture(pp.first, pp.second, picture);
    ofstream metrized_net_out(data);
    net.serialize(metrized_net_out);
}

int main() {
    ifstream in("data/derived_once_survivors.txt");
    int cnt_survivors;
    in >> cnt_survivors;
    vector<int> survivors(cnt_survivors);
    for (int i = 0; i < cnt_survivors; i++)
        in >> survivors[i];
    in.close();

    set<int> possible_nets = {
        72,     278,    434,    521,    691,
        936,    949,    1011,   1034,   1223,
        1351,   1405,   1514,   1829,   2017,
        2119,   2192,   2312,   2454,   2676,
        2700,   2769,   2799,   2891,   3048,
        3222,   3233,   3295,   3338};
    // Slow impossible nets: 729, 1002
    // Slow possible nets: 72, 521, 949, 1829
    vector<int> new_survivors;
    for (int id: survivors) {
        if (!possible_nets.count(id))
            continue;
        if (set<int>{72, 521, 949, 1829}.count(id))
            continue;
        clog << "net id = " << setw(4) << id << ":\t";
        clog.flush();
        auto A = parseDerivedOnceNets("data/derived_once_computer.txt")[id];

        // We skip pathological checks for net 1002 because it completes
        // in roughly 2 hours without them, but takes 10 hours with them.
        bool run_pathological_checks = !possible_nets.count(id) and id != 1002;

        optional<MetrizedNet> forbidden_metrization;
        string metrized_net_filename = "data/metrized_nets/" + to_string(id) + ".txt";
        if (possible_nets.count(id) and
                filesystem::exists(metrized_net_filename)) {
            ifstream metrized_net_in(metrized_net_filename);
            MetrizedNet net2 = parseMetrizedNet(metrized_net_in);
            forbidden_metrization.emplace(net2);

            for (const Dart& d: A.darts()) {
                const Interval I = net2.length(d);
                assert(!I.empty());
                double m = I.mid();
                double len = I.len();
                assert(len < 1e-5);
                len = 2e-4;
                forbidden_metrization->length(d) = Interval(m-len/2, m + len/2);
            }
        }

        long long begin_time = get_time();
        auto [outcome, net] =
            isApproximateMetrizationPossible(A, run_pathological_checks, forbidden_metrization);
        double duration_s = (get_time() - begin_time) / 1e9;

        clog.precision(2);
        clog << "(time = " << fixed << setw(6) << duration_s << ") ";

        if (outcome == INDETERMINATE) {
            clog << "SKIPPED" << endl;
            new_survivors.push_back(id);
            continue;
        }

        if (outcome == 1)
            assert(net.has_value());

        if (forbidden_metrization.has_value()) {
            if (outcome == 0)
                clog << "RIGID" << endl;
            else {
                clog << "NOT RIGID" << endl;
                saveMetrizedNet(*net, "picture_wrong.txt", "data_wrong.txt");
            }
            new_survivors.push_back(id);
            continue;
        }

        assert(outcome == (int)possible_nets.count(id));

        if (outcome == 1) {
            saveMetrizedNet(*net, "data/metrized_nets/" + to_string(id) + "_visualization.txt", "data/metrized_nets/" + to_string(id) + ".txt");
            clog << "POSSIBLE";
            new_survivors.push_back(id);
        }
        clog << endl;
    }

    StopWatch::printState();

    ofstream out("data/derived_once_metric_survivors.txt");
    out << new_survivors.size() << endl;
    for (int id: new_survivors)
        out << id << " ";
    out << endl;
    out.close();
}