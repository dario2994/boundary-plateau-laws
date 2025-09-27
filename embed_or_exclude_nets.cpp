#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "derived_once_nets.cpp"
#include "interval_arithmetic.cpp"
#include "metrized_nets.cpp"
#include "metrization_constraints.cpp"

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

vector<pair<long long, MetrizedNet>> virtual_stack;
long long begin_timestamp = 0;
long long last_timestamp = 0;
const long long timestamp_diff = 10*60; // Ten minutes.
const long long interesting_duration = 60; // One minute.

// Parameters for searchApproximateMetrization.
struct SAMParams {
    vector<Dart> convex_fixed_darts;
    vector<Dart> non_convex_fixed_darts;
    long long begin_time;
    bool run_pathological_checks;
    optional<MetrizedNet> forbidden_metrization;
    bool special_rigidity_2192;
};


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
        const SAMParams& params) {
    if (get_time() - last_timestamp > timestamp_diff * 1e9) {
        last_timestamp = get_time(); // Even if nothing is printed.
        int it = -1;
        for (int i = 0; i < ssize(virtual_stack); i++) {
            if (get_time() - virtual_stack[i].first > interesting_duration * 1e9)
                it = i;
        }
        for (int i = 0; i <= it; i++) {
            long long duration = get_time() - virtual_stack[i].first;
            duration /= 1e9;
            duration /= 60; // In minutes.
            auto zet = virtual_stack[i].second;
            cout << duration << "min: ";
            for (auto d: params.convex_fixed_darts)
                cout << zet.length(d) << " ";
            cout << endl;
        }
        if (it != -1) {
            auto zet = virtual_stack[it].second;
            for (auto d: net.A.darts())
                if (d.from < d.to)
                    cout << d << ": " << zet.length(d) << endl;
            cout << endl;
        }
    }

    const DerivedOnceNet& A = net.A;

    if (params.special_rigidity_2192) {
        net.length(A.dart(16)).intersect(net.length(A.dart(12)));
        net.length(A.dart(12)).intersect(net.length(A.dart(16)));
    }

    // If the computation takes too long, we skip.
    // The slowest net, which is 72, requires TODO hours.
    static const long long one_second = 1e9;
    static const long long one_minute = 60*one_second;
    [[maybe_unused]] static const long long one_hour = 60*one_minute;
    if (get_time() - params.begin_time > 1000*one_hour)
        return {INDETERMINATE, nullopt};

    if (net.totalError() < -0.5)
        return {0, nullopt};

    if (params.forbidden_metrization.has_value()) {
        bool is_refinement = true;
        for (const Dart& d: A.darts()) {
            is_refinement &= isSubset(
                net.length(d),
                params.forbidden_metrization->length(d));
        }
        if (is_refinement)
            return {0, nullopt};
    }

    if (!refineMetrizationConvexContinent(net, params.convex_fixed_darts))
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
                                       and params.run_pathological_checks;
    if (net.totalError() < 1e-5)
        run_now_pathological_checks = true;

    if (run_now_pathological_checks) {
        // The following checks are run only here because:
        //  - They handle rare pathological situations.
        //  - They are _almost_ always implied by other checks.
        //  - They are slow.

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

        // This is the only external overhead necessary for
        // caching minimal paths.
        // It could be part of a custom copy constructor, but it
        // would be even more annoying.
        new_net.incrementVersion();
        virtual_stack.emplace_back(get_time(), new_net);
        auto outcome = searchApproximateMetrization(new_net, params);
        virtual_stack.pop_back();
        if (outcome.first != 0)
            return outcome;
    }

    return {0, nullopt};
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
// INDETERMINATE: The computation was taking too long and it was interrupted
//                before it could prove anything.
int processNet(int net_id) {
    static const set<int> possible_nets = {
        72,     278,    434,    521,    691,
        936,    949,    1011,   1034,   1223,
        1351,   1405,   1514,   1829,   2017,
        2119,   2192,   2312,   2454,   2676,
        2700,   2769,   2799,   2891,   3048,
        3222,   3233,   3295,   3338};

    static const set<int> flexible_nets = {2454,   3048};

    auto A = parseDerivedOnceNets("data/derived_once_computer.txt")[net_id];

    SAMParams params;

    // We skip pathological checks for net 1002 because it completes
    // in roughly 2 hours without them, but takes 10 hours with them.
    params.run_pathological_checks = !possible_nets.count(net_id)
                                     and net_id != 1002;
    params.convex_fixed_darts = findFreeParametersConvex(A);
    params.non_convex_fixed_darts = findFreeParametersNonConvex(A);
    params.special_rigidity_2192 = false;
    params.forbidden_metrization = nullopt;

    MetrizedNet net(A);

    // Necessary to take care of nets that are metrically fixed obviously.
    if (!refineMetrizationConvexContinent(net, A.darts()))
        return 0;

    if (possible_nets.count(net_id)) {
        // The following lines use an existing metrization to skip the
        // existence proof.
        // string metrized_net_filename = "data/metrized_nets/" + to_string(net_id) + ".txt";
        // if (possible_nets.count(net_id) and
        //         filesystem::exists(metrized_net_filename)) {
        //     ifstream metrized_net_in(metrized_net_filename);
        //     MetrizedNet net2 = parseMetrizedNet(metrized_net_in);
        // }

        // Existence
        params.begin_time = get_time();
        auto res = searchApproximateMetrization(net, params);
        assert(res.first == 1 or res.first == INDETERMINATE);
        if (res.first == INDETERMINATE)
            return INDETERMINATE;
        assert(res.second.has_value());
        MetrizedNet net2 = *res.second;
        saveMetrizedNet(net2,
            "data/metrized_nets/" + to_string(net_id) + "_visualization.txt",
            "data/metrized_nets/" + to_string(net_id) + ".txt");
        if (flexible_nets.count(net_id))
            return 1;

        // Rigidity
        relaxMetrization(net2, 2e-4);
        // For net 2192, we prove rigidity by first proving that all lengths
        // but two are rigid.
        // Then, with this information, we prove algebraically that the
        // two missing lengths must be equal.
        // Finally, we rerun the rigidity-algorithm with this additional
        // information to finally obtain the desired rigidity.
        if (net_id == 2192) {
            net2.length(A.dart(16)) = Interval(0, 10);
            net2.length(A.dart(12)) = Interval(0, 10);
        }
        params.forbidden_metrization.emplace(net2);
        params.begin_time = get_time();
        int res2 = searchApproximateMetrization(net, params).first;
        assert(res2 == 0 or res2 == INDETERMINATE);
        if (res2 == INDETERMINATE)
            return INDETERMINATE; // We might also decide to return 1.
        if (net_id == 2192) {
            params.special_rigidity_2192 = true;
            params.begin_time = get_time();
            int res3 = searchApproximateMetrization(net, params).first;
            assert(res3 == INDETERMINATE or res3 == 0);
            if (res3 == INDETERMINATE)
                return INDETERMINATE;
            return 2;
        }
        return 2;
    } else {
        params.begin_time = get_time();
        int res = searchApproximateMetrization(net, params).first;
        assert(res == 0 or res == INDETERMINATE);
        return res;
    }
}


void testingSearchApproximateMetrization(int test_num = -1) {
    struct TestCase {
        int net_id;
        bool run_pathological_checks;
        optional<vector<Interval>> initial_metrization;
        int expected_outcome;
        double expected_duration; // seconds;
        optional<vector<double>> expected_metrization;
    };
    vector<TestCase> cases;

    cases.emplace_back(278, false, nullopt, 1, 0.06, vector<double>{1.627925,1.627925,0.882358,0.882358,1.627925,1.627925,0.882358,0.882358,0.882358,0.882358,0.882358,0.882358,0.145808,0.145808,0.632887,0.632887,0.145810,0.145810,1.627925,1.627925,0.882358,0.882358,1.627925,1.627925,0.882358,0.882358,0.145810,0.145810,1.093351,1.093351,1.093353,1.093353,0.632888,0.632888,0.145808,0.145808,0.323185,0.323185,1.093353,1.093353,1.093351,1.093351,0.632884,0.632884,0.632884,0.632884});

    cases.emplace_back(279, true, nullopt, 0, 0.3, nullopt);

    cases.emplace_back(19, true, nullopt, 0, 0.44, nullopt);

    cases.emplace_back(105, true, nullopt, 0, 0.33, nullopt);

    cases.emplace_back(2119, false, nullopt, 1, 0.13, vector<double>{1.462616,1.462616,1.462616,1.462616,1.016773,1.016773,1.016773,1.016773,1.462616,1.462616,0.236659,0.236659,1.016775,1.016775,1.462614,1.462614,0.236655,0.236655,1.462614,1.462614,1.016775,1.016775,1.462614,1.462614,1.016773,1.016773,1.016775,1.016775,1.016775,1.016775,1.016773,1.016773});

    cases.emplace_back(1025, true, nullopt, 0, 11, nullopt);

    cases.emplace_back(115, true, nullopt, 0, 4.5, nullopt);

    cases.emplace_back(3233, false, nullopt, 1, 0.05, vector<double>{1.230960,1.230960,1.230960,1.230960,1.230959,1.230959,1.230961,1.230961,1.230960,1.230960,1.230960,1.230960,1.230959,1.230959,1.230960,1.230960,1.230959,1.230959,1.230958,1.230958});

    // cases.emplace_back(3091, true, nullopt, 0, 200, nullopt); // SLOW

    for (int test_id = 0; test_id < ssize(cases); test_id++) {
        if (test_num != -1 and test_id != test_num)
            continue;
        const TestCase& tc = cases[test_id];
        cout << "Test " << test_id << ", net " << tc.net_id << ": ";
        cout.flush();
        auto A = parseDerivedOnceNets("data/derived_once_computer.txt")[tc.net_id];

        SAMParams params;

        params.run_pathological_checks = tc.run_pathological_checks;
        params.convex_fixed_darts = findFreeParametersConvex(A);
        params.non_convex_fixed_darts = findFreeParametersNonConvex(A);
        params.special_rigidity_2192 = false;
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
        auto [outcome, maybe_metrized_net] = searchApproximateMetrization(net, params);
        double elapsed_time = stopwatch.get();

        assert(outcome == tc.expected_outcome);
        cout << "outcome = " << outcome << ", ";
        cout << "elapsed = "<< elapsed_time;
        double time_diff = (tc.expected_duration - elapsed_time) / tc.expected_duration;
        if (time_diff > 0.1)
            cout << " (improved from " << tc.expected_duration << ")";
        if (time_diff < -0.1)
            cout << " (degraded from " << tc.expected_duration << ")";

        if (outcome == 1) {
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

    ifstream in("data/derived_once_survivors.txt");
    int cnt_survivors;
    in >> cnt_survivors;
    vector<int> survivors(cnt_survivors);
    for (int i = 0; i < cnt_survivors; i++)
        in >> survivors[i];
    in.close();

    // Slow impossible nets: 729, 1002
    // Slow possible nets: 72, 521, 949, 1829, 2192
    vector<int> new_survivors;
    for (int id: {72}) {
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
        else if (outcome == INDETERMINATE)
            clog << "INDETERMINATE";
        else assert(0);
        clog << endl;

        if (outcome != 0)
            new_survivors.push_back(id);
    }
    clog << endl;

    StopWatch::printState();

    ofstream out("data/derived_once_metric_survivors.txt");
    out << new_survivors.size() << endl;
    for (int id: new_survivors)
        out << id << " ";
    out << endl;
    out.close();
}