#ifndef METRIZED_NETS
#define METRIZED_NETS

#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "derived_once_nets.cpp"
#include "interval_arithmetic.cpp"
#include "quaternions.cpp"
#include "utils.cpp"


struct PathStep {
    char type;
    double rotation;
    Dart d;
    
    PathStep(double _rotation): type('R'), rotation(_rotation) {}
    PathStep(Dart _d): type('E'), d(_d) {}
    
    double cost(const vector<Interval>& lengths) {
        return 1e-10 + ((type == 'R') ? 0.0 : lengths[d.id].len());
    }
};

ostream& operator<<(ostream& out, const PathStep& x) {
    if (x.type == 'R') {
        int angle = round(x.rotation * 360 / (2*M_PI));
        out << "R(" << angle << ")";
    } else {
        out << "E(" << x.d.from << "->" << x.d.to << ")"; 
    }
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
    MetrizedNet(const DerivedOnceNet& _A): A(_A),
            optimal_path(make_shared<vector<vector<vector<PathStep>>>>()),
            to_pins(make_shared<vector<vector<vector<PathStep>>>>()),
            current_version(make_shared<long long>()),
            last_update(make_shared<vector<long long>>()) {
        optimal_path->resize(2*A.E());
        to_pins->resize(2*A.E());
        *current_version = 1;
        last_update->resize(2*A.E(), 0);
        
        _lengths.resize(2*A.E(), Interval(0, M_PI));
        
        for (Dart d: A.darts()) {
            Interval& I = _lengths[d.id];
            if (A.isConvexFace(d)) {
                I.intersect(Interval(0.02, M_PI-0.01));
                int q = ssize(A.face(d));
                assert(q != 2);
                if (q == 5)
                    I.intersect(Interval(0, acos(-1.0/3)-0.0001));
            }
        }
        for (Dart d: A.darts())
            _lengths[d.id].intersect(_lengths[A.rev(d).id]);
    }
    
    // This should be used.
    Interval& length(const Dart& d) {
        return _lengths[d.id];
    }
    const Interval& length(const Dart& d) const {
        return _lengths[d.id];
    }
    
    const vector<Interval>& lengths() const { return _lengths; }
    
    void incrementVersion() const {
        (*current_version)++;
    }
    
    Quaternion<Interval> isometryFromPath(const vector<PathStep>& path) const {
        StopWatch _("isometryFromPath");
        
        double sum_theta = 0;
        Quaternion<Interval> res(1);
        for (auto op: path) {
            if (op.type == 'R')
                sum_theta += op.rotation;
            if (op.type == 'E')
                res *= basicIsometry<Interval>(sum_theta, _lengths[op.d.id]);
        }
        res *= basicRotation<Interval>(sum_theta);
        
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
                               {dist + _lengths[d.id].len() + 1e-10, A.rev(d)}};
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
            if (!path.empty() and path.back().type == 'R' and op.type == 'R')
                path.back().rotation += op.rotation;
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
                    append_to_path(path, A.nextAngle(curr));
                } else if (father[curr.id] == 'r') { // Exponential + Rotation by 180.
                    append_to_path(path, M_PI);
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
            append_to_path(path, 2*M_PI/3);
            (*to_pins)[root.id].push_back(path);
        }
    }
    
  public:
    Quaternion<Interval> getPath(const Dart& from, const Dart& to) const {
        rebuildEmbedding(from);
        return isometryFromPath((*optimal_path)[from.id][to.id]);
    }
    
    Quaternion<Interval> getPathToPin(const Dart& from, const Dart& to) const {
        rebuildEmbedding(from);
        return isometryFromPath((*to_pins)[from.id][A.concaveAngleId(to)]);
    }
    
};

#endif
