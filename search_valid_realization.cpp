#ifndef SEARCH_VALID_REALIZATION
#define SEARCH_VALID_REALIZATION

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "marked_combinatorial_maps.cpp"
#include "interval_arithmetic.cpp"
#include "metrized_nets.cpp"
#include "metrization_constraints.cpp"


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

// Returns a subset of the darts such that them and their reverse belong to the
// boundary of non-convex faces.
// If we know the lengths of all edges on the boundary of a convex face
// and also the edges in the returned subset, then the length of all edges
// can be determined with simple inference rules.
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
    optional<ApproximatelyMetrizedNet> forbidden_metrization;
    bool special_rigidity_736;
};

// Takes an approximately metrized net and returns:
// Nothing: We can prove that no refinement of the given approximately metrized
//          net is the link of a minimizing cone.
// A metrized net: A realization of the given approximately metrized net that
//                 is (approximately) a stationary geodesic net with two
//                 antipodal pins without short edges.
optional<ApproximatelyMetrizedNet> searchValidRealization(
        ApproximatelyMetrizedNet net,
        const SAMParams& params) {
    const MarkedCombinatorialMap& A = net.A;

    if (params.special_rigidity_736) {
        net.length(A.dart(16)).intersect(net.length(A.dart(12)));
        net.length(A.dart(12)).intersect(net.length(A.dart(16)));
    }


    if (net.isImpossible())
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

        auto maybe_res = searchValidRealization(new_net, params);
        if (maybe_res)
            return maybe_res;
    }

    return nullopt;
}

#endif