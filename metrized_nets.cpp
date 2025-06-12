#ifndef METRIZED_NETS
#define METRIZED_NETS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "derived_once_nets.cpp"
#include "interval_arithmetic.cpp"
#include "quaternions.cpp"
#include "interval_geometry.cpp"
#include "utils.cpp"

bool containsBiangle(const DerivedOnceNet& A) {
    for (const vector<Dart>& f: A.faces()) {
        if (A.faceType(f) == pair<int,int>{2, 0})
            return true;
    }
    return false;
}

int pinContribution(const DerivedOnceNet& A, Dart d) {
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
bool isShorterThanPi(const DerivedOnceNet& A, const Dart& d) {
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


struct PathStep {
    char type;
    double t;
    Dart d;

    PathStep(char _type, double _t): type(_type), t(_t) {
        assert(type == 'R' or type == 'E');
    }
    PathStep(Dart _d): type('D'), d(_d) {}

    double cost(const vector<Interval>& lengths) {
        return 1e-10 + ((type == 'D') ? lengths[d.id].len() : 0.0 );
    }
};

ostream& operator<<(ostream& out, const PathStep& x) {
    if (x.type == 'D') {
        out << "D(" << x.d.from << "->" << x.d.to << ")";
    } else if (x.type == 'R') {
        int angle = round(x.t * 360 / (2*M_PI));
        out << "R(" << angle << ")";
    } else if (x.type == 'E') {
        out << "E(" << x.t << ")";
    } else assert(0);

    return out;
}

class MetrizedNet {
  public:
    const DerivedOnceNet& A;

  private:
    vector<Interval> _lengths;
    // TODO: If these are too slow, we can use raw pointers.
    // Defined as shared_ptr to avoid the cost of copying.
    shared_ptr<vector<vector<vector<PathStep>>>> optimal_path;
    shared_ptr<vector<vector<vector<PathStep>>>> to_pins;
    // Indexed by initial dart.
    shared_ptr<long long> current_version;
    shared_ptr<vector<long long>> last_update;

  public:
    Interval& length(const Dart& d) {
        return _lengths[d.id>>1];
    }
    const Interval& length(const Dart& d) const {
        return _lengths[d.id>>1];
    }

    MetrizedNet(const DerivedOnceNet& _A): A(_A),
            optimal_path(make_shared<vector<vector<vector<PathStep>>>>()),
            to_pins(make_shared<vector<vector<vector<PathStep>>>>()),
            current_version(make_shared<long long>()),
            last_update(make_shared<vector<long long>>()) {
        optimal_path->resize(2*A.E());
        to_pins->resize(2*A.E());
        *current_version = 1;
        last_update->resize(2*A.E(), 0);

        _lengths.resize(A.E(), Interval(0, 2*M_PI));

        for (Dart d: A.darts()) {
            Interval& I = length(d);

            if (isShorterThanPi(A, d))
                I.intersect(Interval(0, M_PI));

            if (A.isConvexFace(d)) {
                int q = ssize(A.face(d));
                assert(2 <= q and q <= 5);

                static const double ell_T = acos(-1.0/3);

                if (q == 2)
                    I.intersect(Interval(M_PI));
                if (q == 3)
                    I.intersect(Interval(ell_T));
                if (q == 4)
                    I.intersect(Interval(0.01, M_PI-0.001));
                if (q == 5)
                    I.intersect(Interval(0, ell_T-1e-5));
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

    void incrementVersion() const {
        (*current_version)++;
    }

    ApproximateIsometry isometryFromPath(const vector<PathStep>& path) const {
        StopWatch _("isometryFromPath");

        double sum_theta = 0;
        ApproximateIsometry res(1);
        for (auto op: path) {
            if (op.type == 'R')
                sum_theta += op.t;
            else if (op.type == 'E')
                res *= basicTranslation(sum_theta, op.t);
            else if (op.type == 'D')
                res *= basicTranslation(sum_theta, length(op.d));
            else assert(0);
        }
        res *= basicRotation(sum_theta);

        return res;
    }

  private:
    bool rebuildPolicy(const  Dart& root) const {
        long long curr = *current_version;
        long long& last = (*last_update)[root.id];
        bool outcome = curr - last >= 100 or 2 * last <= curr or curr <= 2;
        if (outcome)
            last = curr;
        return outcome;
    }

    void rebuildEmbedding(const Dart& root) const {
        if (!rebuildPolicy(root))
            return;
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

        // To compress consecutive rotations in a single one.
        auto append_to_path = [&](vector<PathStep>& path, PathStep op) {
            if (!path.empty() and path.back().type == op.type
                and op.type != 'D')
                path.back().t += op.t;
            else
                path.push_back(op);
        };

        (*optimal_path)[root.id].clear();
        for (const Dart& d: A.darts()) {
            Dart curr = d;
            vector<PathStep> path;

            // We go from last_dart to root.
            while (father[curr.id] != '-') {
                if (father[curr.id] == 'n') { // Rotation by 120 or 240.
                    curr = A.prev(curr);
                    append_to_path(path, PathStep('R', A.nextAngle(curr)));
                } else if (father[curr.id] == 'r') { // Exponential + Rotation by 180.
                    append_to_path(path, PathStep('R', M_PI));
                    curr = A.rev(curr);
                    append_to_path(path, curr);
                } else assert(0);
            }
            reverse(path.begin(), path.end());
            (*optimal_path)[root.id].push_back(path);
        }

        // If d belongs to dart_to_pin, then rotating d by 120 yields
        // a dart going to a pin.
        (*to_pins)[root.id].clear();
        for (Dart d: A.concaveAngles()) {
            assert(A.isConcaveAngle(d));
            vector<PathStep> path = (*optimal_path)[root.id][A.rev(d).id];
            append_to_path(path, PathStep('R', 2*M_PI/3));
            (*to_pins)[root.id].push_back(path);
        }
    }

  public:
    vector<PathStep> getPath(const Dart& from, const Dart& to) const {
        rebuildEmbedding(from);
        return (*optimal_path)[from.id][to.id];
    }

    vector<PathStep> getPathToPin(const Dart& from, const Dart& to) const {
        rebuildEmbedding(from);
        return (*to_pins)[from.id][A.concaveAngleId(to)];
    }

    ApproximateIsometry getIsometry(const Dart& from, const Dart& to) const {
        rebuildEmbedding(from);
        return isometryFromPath((*optimal_path)[from.id][to.id]);
    }

    ApproximateIsometry getIsometryToPin(const Dart& from, const Dart& to) const {
        rebuildEmbedding(from);
        return isometryFromPath((*to_pins)[from.id][A.concaveAngleId(to)]);
    }

};

#endif
