#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "marked_combinatorial_maps.cpp"
#include "interval_arithmetic.cpp"
#include "metrized_nets.cpp"
#include "metrization_constraints.cpp"
#include "search_approximate_metrization.cpp"


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
}
