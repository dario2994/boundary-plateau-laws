#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "marked_combinatorial_maps.cpp"
#include "interval_arithmetic.cpp"
#include "metrized_nets.cpp"
#include "metrization_constraints.cpp"
#include "search_approximate_metrization.cpp"

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

    static const set<int> flexible_nets = {827, 1029};

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
        relaxMetrization(net2, 2e-4); // TODO: Is this constant good enough?
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

int main() {
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
