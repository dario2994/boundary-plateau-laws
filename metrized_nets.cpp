#ifndef METRIZED_NETS
#define METRIZED_NETS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "marked_combinatorial_maps.cpp"
#include "interval_arithmetic.cpp"
#include "interval_geometry.cpp"
#include "utils.cpp"

bool containsBiangle(const MarkedCombinatorialMap& A) {
    for (const vector<Dart>& f: A.faces()) {
        if (A.faceType(f) == pair<int,int>{2, 0})
            return true;
    }
    return false;
}

int pinContribution(const MarkedCombinatorialMap& A, Dart d) {
    if (A.faceType(d).second != 2)
        return 0;

    int dist = -1;
    Dart curr = d;
    while (1) {
        if (A.isConcaveAngle(curr)) {
            if (dist == -1) dist = 1;
            else break;
        } else if (dist >= 0) dist++;

        curr = A.phi(curr);
    }
    int contribution = max(4-dist, 0);
    return contribution;
}


// If this function returns true, then we can prove that the length of d is <= pi.
bool isShorterThanPi(const MarkedCombinatorialMap& A, const Dart& d) {
    if (containsBiangle(A))
        return true;

    Dart e = A.rev(d);
    if (A.isConvexFace(d) or A.isConvexFace(e))
        return true;

    if (A.faceArea(d) + A.faceArea(e) <=
        6 + pinContribution(A, d) + pinContribution(A, e))
        return true;

    return false;
}

// A step in a path on a sphere made of rotations and translations.
//
// - PathStep('R', t): Corresponds to a rotation by t (\rot_t in the paper).
// - PathStep('E', t): Corresponds to the exponential for time t
//   (\exp_t in the paper).
//
// Additionally, there's also PathStep(Dart) which corresponds
// to PathStep('E', Dart.len()).
// This is a trick necessary to speed up some computation.
struct PathStep {
    char type;
    int deg;
    double len;
    Dart d;

    static PathStep rot(int _deg) {
        PathStep step; step.type = 'R'; step.deg = _deg; return step;
    }
    // TODO: Sarebbe ideale cancellare questo del tutto!
    // Viene usato soltanto una volta.
    static PathStep exp(double _len) {
        PathStep step; step.type = 'E'; step.len = _len; return step;
    }
    static PathStep exp(Dart _d) {
        PathStep step; step.type = 'D'; step.d = _d; return step;
    }

    double cost(const vector<Interval>& lengths) {
        return 1e-10 + ((type == 'D') ? lengths[d.id].len() : 0.0 );
    }
};

ostream& operator<<(ostream& out, const PathStep& x) {
    if (x.type == 'D') {
        out << "D(" << x.d.from << "->" << x.d.to << ")";
    } else if (x.type == 'R') {
        out << "R(" << x.deg << ")";
    } else if (x.type == 'E') {
        out << "E(" << x.len << ")";
    } else assert(0);

    return out;
}

// To compress consecutive rotations in a single one.
void appendToPath(vector<PathStep>& path, const PathStep& op) {
    if (path.empty() or path.back().type != op.type or op.type == 'D')
        path.push_back(op);
    else if (op.type == 'R')
        path.back().deg += op.deg;
    else if (op.type == 'E')
        path.back().len += op.len;
}

class MetrizedNet {
  public:
    const MarkedCombinatorialMap A;

  private:
    vector<Interval> _lengths;

    // Cache for the optimal paths (i.e., the ones with smallest total error)
    // between darts.
    mutable vector<vector<vector<PathStep>>> optimal_path;


    void initOptimalPaths() {
        optimal_path.clear();
        optimal_path.resize(2*A.E());
    }

  public:

    MetrizedNet(const MetrizedNet& other) : A(other.A), _lengths(other._lengths) {
        initOptimalPaths();
    }

    // Deleting the copy-assignment operator (since our copy-logic is not the default one).
    MetrizedNet& operator=(const MetrizedNet& other) = delete;

    Interval& length(const Dart& d) {
        return _lengths[d.id>>1];
    }
    const Interval& length(const Dart& d) const {
        return _lengths[d.id>>1];
    }

    MetrizedNet(const MarkedCombinatorialMap& _A): A(_A) {
        initOptimalPaths();

        _lengths.resize(A.E(), Interval(0, (2*pi).r));

        for (Dart d: A.darts()) {
            Interval& I = length(d);

            if (isShorterThanPi(A, d))
                I.intersect(Interval(0, pi.r));

            if (A.isConvexFace(d)) {
                int q = ssize(A.face(d));
                assert(2 <= q and q <= 5);

                static const double ell_T = acos(-1.0/3);

                if (q == 2)
                    I.intersect(pi);
                if (q == 3)
                    I.intersect(paddedInterval(ell_T));
                if (q == 4)
                    I.intersect(Interval(0.01, pi.r - 1e-3));
                if (q == 5)
                    I.intersect(Interval(0.01, ell_T - 1e-3));
            }
        }
    }

    double totalError() const {
        double res = 0;
        for (const Interval& I: _lengths) {
            if (I.empty())
                return -1;
            res += I.len();
        }
        return res;
    }

    bool isImpossible() const {
        for (const Interval& I: _lengths)
            if (I.empty())
                return true;
        return false;
    }

    ApproximateIsometry isometryFromPath(const vector<PathStep>& path) const {
        StopWatch _("isometryFromPath");

        int sum_theta = 0;
        ApproximateIsometry res(1, 1e-10);
        for (auto op: path) {
            if (op.type == 'R')
                sum_theta += op.deg;
            else if (op.type == 'E')
                res *= basicTranslation(sum_theta, op.len);
            else if (op.type == 'D')
                res *= basicTranslation(sum_theta, length(op.d));
            else assert(0);
        }
        res *= basicRotation(sum_theta);

        return res;
    }

  private:

    void computeOptimalPaths(const Dart& root) const {
        if (!optimal_path[root.id].empty())
            return;
        StopWatch _("computeOptimalPaths"); // TODO: Remove all the StopWatch instances.
        // Run Dijkstra from the root.
        // The cost of traversing a dart is the length of the corresponding
        // interval.
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
                               {dist + length(d).len() + 1e-10, A.rev(d)}};
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

        optimal_path[root.id].clear();
        for (const Dart& d: A.darts()) {
            Dart curr = d;
            vector<PathStep> path;

            // We go from last_dart to root.
            while (father[curr.id] != '-') {
                if (father[curr.id] == 'n') { // Rotation by 120 or 240.
                    curr = A.prev(curr);
                    appendToPath(path, PathStep::rot(A.nextAngle(curr)));
                } else if (father[curr.id] == 'r') { // Exponential + Rotation by 180.
                    appendToPath(path, PathStep::rot(180));
                    curr = A.rev(curr);
                    appendToPath(path, PathStep::exp(curr));
                } else assert(0);
            }
            reverse(path.begin(), path.end());
            optimal_path[root.id].push_back(path);
        }
    }

  public:
    vector<PathStep> getPath(const Dart& from, const Dart& to) const {
        computeOptimalPaths(from);
        return optimal_path[from.id][to.id];
    }

    vector<PathStep> getPathToPin(const Dart& from, const Dart& to) const {
        assert(A.isConcaveAngle(to));
        vector<PathStep> path = getPath(from, A.rev(to));
        appendToPath(path, PathStep::rot(120));
        return path;
    }

    vector<PathStep> getPathFromPinToPin(const Dart& from, const Dart& to) const {
        assert(A.isConcaveAngle(from));
        assert(A.isConcaveAngle(to));
        vector<PathStep> path = getPath(A.rev(from), A.rev(to));
        reverse(path.begin(), path.end());
        appendToPath(path, PathStep::rot(-120));
        reverse(path.begin(), path.end());
        appendToPath(path, PathStep::rot(120));
        return path;
    }

    ApproximateIsometry getIsometry(const Dart& from, const Dart& to) const {
        return isometryFromPath(getPath(from, to));
    }

    ApproximateIsometry getIsometryToPin(const Dart& from, const Dart& to) const {
        return isometryFromPath(getPathToPin(from, to));
    }

    ApproximateIsometry getIsometryFromPinToPin(const Dart& from, const Dart& to) const {
        return isometryFromPath(getPathFromPinToPin(from, to));
    }

    void serialize(ostream& out) const {
        A.serialize(out);

        for (const Dart& d: A.darts()) {
            Interval I = length(d);
            out << std::format("{:.17g} {:.17g}\n", I.l, I.r);
        }
    }

};

MetrizedNet parseMetrizedNet(istream& in) {
    MarkedCombinatorialMap A = parseMarkedCombinatorialMap(in);
    MetrizedNet net(A);
    for (Dart d: A.darts()) {
        double l, r;
        in >> l >> r;
        net.length(d) = Interval(l, r);
    }
    return net;
}

#endif
