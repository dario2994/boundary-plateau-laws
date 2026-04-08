#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "marked_combinatorial_maps.cpp"
#include "interval_arithmetic.cpp"
#include "metrized_nets.cpp"
#include "metrization_constraints.cpp"

// Move the following two functions somewhere else. Maybe in metrized_nets.
void savePicture(vector<pair<Point, string>> vertices,
                 vector<tuple<Point, Point, double>> edges, string out_file) {
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

pair<vector<pair<Point, string>>, vector<tuple<Point, Point, double>>>
        embedNet(const MetrizedNet& net) {
    const MarkedCombinatorialMap& A = net.A;

    vector<pair<Point, string>> vertex_embedding(A.V());
    vector<tuple<Point, Point, double>> dart_embedding(2 * A.E());
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
int lengthsImplications(const MarkedCombinatorialMap& A,
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
// can be uniquely determined through simple inference rules.
// Among the valid possible sets, it returns one that determines
// all the lengths in a "numerically robust" way.
vector<Dart> findFreeParametersConvex(const MarkedCombinatorialMap& A) {
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

vector<Dart> findFreeParametersNonConvex(const MarkedCombinatorialMap& A) {
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

// Parameters for searchApproximateMetrization.
struct SAMParams {
    vector<Dart> convex_fixed_darts;
    vector<Dart> non_convex_fixed_darts;
    long long begin_time;
    bool run_pathological_checks;
    optional<MetrizedNet> forbidden_metrization;
    bool special_rigidity_736;
};

// Takes an approximately metrized once-derived net and returns:
//  Nothing: No refinement of the given approximate metrization is the link
//           of a stable cone.
//  A net: A refinement of the given approximate metrization that is precise enough
//         (error 1e-5), has antipodal pins, correct angles, and no extremely
//         short side.
optional<MetrizedNet> searchApproximateMetrization(
        MetrizedNet net,
        const SAMParams& params) {
    const MarkedCombinatorialMap& A = net.A;

    if (params.special_rigidity_736) {
        net.length(A.dart(16)).intersect(net.length(A.dart(12)));
        net.length(A.dart(12)).intersect(net.length(A.dart(16)));
    }


    // TODO: Make this more reasonable.
    if (net.totalError() < -0.5)
        return nullopt;

    if (params.forbidden_metrization.has_value()) {
        bool is_refinement = true;
        for (const Dart& d: A.darts()) {
            is_refinement &= isSubset(
                net.length(d),
                params.forbidden_metrization->length(d));
        }
        if (is_refinement)
            return nullopt;
    }

    if (!refineMetrizationConvexContinent(net, params.convex_fixed_darts))
        return nullopt;

    if (!refineMetrizationNonConvexFaces(net))
        return nullopt;

    if (!checkPinsAntipodality(net))
        return nullopt;

    if (!checkPinAngles(net))
        return nullopt;

    if (!checkNonIntersection(net))
        return nullopt;

    if (!checkPinOutsideConvexFaces(net))
        return nullopt;

    if (!checkLengthLowerbound(net))
        return nullopt;

    bool run_now_pathological_checks = net.totalError() < 1e-2
                                       and params.run_pathological_checks;
    if (net.totalError() < 1e-5)
        run_now_pathological_checks = true;

    if (run_now_pathological_checks) {
        // The following checks are run only here because:
        //  - They handle rare pathological situations.
        //  - They are _almost_ always implied by other checks.
        //  - They are slow.
        if (!checkPinsAntipodality2(net))
            return nullopt;

        if (!checkFacesAreBoundaries(net))
            return nullopt;
    }

    // If we reached the desired precision of the metrization, we
    // end the search.
    if (net.totalError() < 1e-5)
        return net;


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
    for (const Dart& d: params.non_convex_fixed_darts) {
        double l = net.length(d).len();
        max_non_convex = max(max_non_convex, l);
    }

    // Decide whether we want to refine (i.e., half) the precision
    // only on convex_fixed_darts or also on non_convex_fixed_darts.
    vector<Dart> fixed_darts = params.convex_fixed_darts;
    if (max_non_convex > max_error_boundary_continent) {
        for (const Dart& d: params.non_convex_fixed_darts)
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

        auto maybe_res = searchApproximateMetrization(new_net, params);
        if (maybe_res)
            return maybe_res;
    }

    return nullopt;
}

void saveMetrizedNet(const MetrizedNet& net, string picture, string data) {
    auto pp = embedNet(net);
    savePicture(pp.first, pp.second, picture);
    ofstream metrized_net_out(data);
    net.serialize(metrized_net_out);
}

void relaxMetrization(MetrizedNet& net, double len) {
    for (const Dart& d: net.A.darts()) {
        const Interval I = net.length(d);
        assert(!I.empty());
        double m = I.mid();
        net.length(d) = Interval(m-len/2, m + len/2);
    }
}

// Return values:
// 0: Impossible.
// 1: Possible.
// 2: Possible and rigid.
int processNet(int net_id) {
    static const set<int> possible_nets = {
        28, 99, 158, 182, 250,
        319, 325, 346, 361, 421,
        468, 484, 520, 614, 681,
        712, 736, 779, 827, 909,
        920, 938, 1029, 1096, 1100,
        1122};

    static const set<int> flexible_nets = {827,   1029};

    auto A = parseMarkedCombinatorialMaps("data/once_derived_computer.txt")[net_id];

    SAMParams params;

    params.run_pathological_checks = !possible_nets.count(net_id);
    params.convex_fixed_darts = findFreeParametersConvex(A);
    params.non_convex_fixed_darts = findFreeParametersNonConvex(A);
    params.special_rigidity_736 = false;
    params.forbidden_metrization = nullopt;

    MetrizedNet net(A);

    // Necessary to take care of nets that are metrically fixed obviously.
    if (!refineMetrizationConvexContinent(net, A.darts()))
        return 0;

    if (possible_nets.count(net_id)) {
        // Existence
        params.begin_time = get_time();
        auto res = searchApproximateMetrization(net, params);
        assert(res.has_value());
        MetrizedNet net2 = *res;
        saveMetrizedNet(net2,
            "data/metrized_nets/" + to_string(net_id) + "_visualization.txt",
            "data/metrized_nets/" + to_string(net_id) + ".txt");
        if (flexible_nets.count(net_id))
            return 1;

        // Rigidity
        relaxMetrization(net2, 2e-4);
        // For net 736, we prove rigidity by first proving that all lengths
        // but two are rigid.
        // Then, with this information, we prove algebraically that the
        // two missing lengths must be equal.
        // Finally, we rerun the rigidity-algorithm with this additional
        // information to obtain the desired rigidity.
        if (net_id == 736) {
            net2.length(A.dart(16)) = Interval(0, 10);
            net2.length(A.dart(12)) = Interval(0, 10);
        }
        params.forbidden_metrization.emplace(net2);
        params.begin_time = get_time();
        auto res2 = searchApproximateMetrization(net, params);
        assert(!res2.has_value());

        if (net_id == 736) {
            params.special_rigidity_736 = true;
            params.begin_time = get_time();
            auto res3 = searchApproximateMetrization(net, params);
            assert(!res3.has_value());
            params.forbidden_metrization.reset();
            auto res4 = searchApproximateMetrization(net, params);
            assert(res4.has_value());
            saveMetrizedNet(*res4,
                "data/metrized_nets/" + to_string(net_id) + "_visualization.txt",
                "data/metrized_nets/" + to_string(net_id) + ".txt");
            return 2;
        }
        return 2;
    } else {
        params.begin_time = get_time();
        auto res = searchApproximateMetrization(net, params);
        assert(!res.has_value());
        return 0;
    }
}

// TODO: Move in a different file.
void testingSearchApproximateMetrization(int test_num = -1) {
    struct TestCase {
        int net_id;
        bool run_pathological_checks;
        optional<vector<Interval>> initial_metrization;
        bool expected_outcome;
        double expected_duration; // seconds;
        optional<vector<double>> expected_metrization;
    };
    vector<TestCase> cases;

    cases.emplace_back(99, false, nullopt, 1, 0.06, vector<double>{1.627925,1.627925,0.882358,0.882358,1.627925,1.627925,0.882358,0.882358,0.882358,0.882358,0.882358,0.882358,0.145808,0.145808,0.632887,0.632887,0.145810,0.145810,1.627925,1.627925,0.882358,0.882358,1.627925,1.627925,0.882358,0.882358,0.145810,0.145810,1.093351,1.093351,1.093353,1.093353,0.632888,0.632888,0.145808,0.145808,0.323185,0.323185,1.093353,1.093353,1.093351,1.093351,0.632884,0.632884,0.632884,0.632884});

    cases.emplace_back(100, true, nullopt, 0, 0.3, nullopt);

    cases.emplace_back(11, true, nullopt, 0, 0.2, nullopt);

    cases.emplace_back(39, true, nullopt, 0, 0.33, nullopt);

    cases.emplace_back(712, false, nullopt, 1, 0.13, vector<double>{1.462616,1.462616,1.462616,1.462616,1.016773,1.016773,1.016773,1.016773,1.462616,1.462616,0.236659,0.236659,1.016775,1.016775,1.462614,1.462614,0.236655,0.236655,1.462614,1.462614,1.016775,1.016775,1.462614,1.462614,1.016773,1.016773,1.016775,1.016775,1.016775,1.016775,1.016773,1.016773});

    cases.emplace_back(356, true, nullopt, 0, 10, nullopt);

    cases.emplace_back(45, true, nullopt, 0, 4.5, nullopt);

    cases.emplace_back(66, true, nullopt, 0, 20, nullopt);

    cases.emplace_back(1100, false, nullopt, 1, 0.05, vector<double>{1.230960,1.230960,1.230960,1.230960,1.230959,1.230959,1.230961,1.230961,1.230960,1.230960,1.230960,1.230960,1.230959,1.230959,1.230960,1.230960,1.230959,1.230959,1.230958,1.230958});

    cases.emplace_back(1047, true, nullopt, 0, 10, nullopt);

    for (int test_id = 0; test_id < ssize(cases); test_id++) {
        if (test_num != -1 and test_id != test_num)
            continue;
        const TestCase& tc = cases[test_id];
        cout << "Test " << test_id << ", net " << tc.net_id << ": ";
        cout.flush();
        auto A = parseMarkedCombinatorialMaps("data/once_derived_computer.txt")[tc.net_id];

        SAMParams params;

        params.run_pathological_checks = tc.run_pathological_checks;
        params.convex_fixed_darts = findFreeParametersConvex(A);
        params.non_convex_fixed_darts = findFreeParametersNonConvex(A);
        params.special_rigidity_736 = false;
        params.forbidden_metrization = nullopt;

        MetrizedNet net(A);

        if (tc.initial_metrization) {
            assert(ssize(*tc.initial_metrization) == ssize(params.convex_fixed_darts));
            for (int i = 0; i < ssize(*tc.initial_metrization); i++) {
                Dart d = params.convex_fixed_darts[i];
                net.length(d).intersect(tc.initial_metrization.value()[i]);
            }
        }

        StopWatch stopwatch;
        auto maybe_metrized_net = searchApproximateMetrization(net, params);
        double elapsed_time = stopwatch.get();

        assert(maybe_metrized_net.has_value() == tc.expected_outcome);
        cout << "outcome = " << maybe_metrized_net.has_value() << ", ";
        cout << "elapsed = "<< elapsed_time;
        double time_diff = (tc.expected_duration - elapsed_time) / tc.expected_duration;
        if (time_diff > 0.1)
            cout << " (improved from " << tc.expected_duration << ")";
        if (time_diff < -0.1)
            cout << " (degraded from " << tc.expected_duration << ")";

        if (maybe_metrized_net.has_value()) {
            const double acceptable_total_error = 1e-5;
            assert(maybe_metrized_net->totalError() < acceptable_total_error);
            if (tc.expected_metrization) {
                for (Dart d: A.darts()) {
                    double expected_len = (*tc.expected_metrization)[d.id];
                    Interval outcome_len = maybe_metrized_net->length(d);
                    const double acceptable_error = 1e-4;
                    assert(outcome_len.l - acceptable_error < expected_len);
                    assert(expected_len < outcome_len.r + acceptable_error);
                }
            } else {
                vector<double> lens;
                for (auto d: maybe_metrized_net->A.darts())
                    lens.push_back(maybe_metrized_net->length(d).mid());
                cout << " Info to add to the test: ";
                cout << "vector<double>{";
                cout.precision(6);
                cout << fixed << lens[0];
                for (int i = 1; i < ssize(lens); i++)
                    cout << fixed << "," << lens[i];
                cout << "}";
            }
        }
        cout << endl;
    }

    cout << endl;
    StopWatch::printState();
}

int main() {
    testingSearchApproximateMetrization();
    return 0;
    ifstream in("data/once_derived_computer.txt");
    int cnt_once_derived;
    in >> cnt_once_derived;
    vector<int> once_derived(cnt_once_derived);
    for (int i = 0; i < cnt_once_derived; i++)
        in >> once_derived[i];
    in.close();

    vector<int> embeddable_once_derived;

    for (int id = 0; id < ssize(once_derived); id++) {
        clog << "net id = " << setw(4) << id << ":\t";
        clog.flush();

        long long begin_time = get_time();
        int outcome = processNet(id);
        double duration_s = (get_time() - begin_time) / 1e9;

        clog.precision(2);
        clog << "(time = " << fixed << setw(10) << duration_s << ") ";

        if (outcome == 0)
            clog << "IMPOSSIBLE";
        else if (outcome == 1)
            clog << "POSSIBLE";
        else if (outcome == 2)
            clog << "RIGID";
        else assert(0);
        clog << endl;

        if (outcome != 0)
            embeddable_once_derived.push_back(id);
    }
    clog << endl; // TODO: Use clog everywhere?

    ofstream out("data/embeddable_nets.txt");
    out << embeddable_once_derived.size() << endl;
    for (int id: embeddable_once_derived)
        out << id << " ";
    out << endl;
    out.close();
}
