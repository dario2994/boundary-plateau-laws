#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "interval_arithmetic.cpp"
#include "quaternions.cpp"
#include "derived_once_nets.cpp"
#include "metrized_nets.cpp"
#include "interval_geometry.cpp"

typedef array<double, 3> pt;

ostream& operator<<(ostream& out, pt p) {
    out << p[0] << " " << p[1] << " " << p[2];
    return out;
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

MetrizedNet parseMetrizedNet(const DerivedOnceNet& A, string in_file) {
    MetrizedNet net(A);
    ifstream in(in_file);
    assert(in.is_open());

    for (Dart d: A.darts()) {
        double l, r;
        in >> l >> r;
        net.length(d) = Interval(l, r);
    }
    return net;
}


vector<pt> getVertexEmbedding(const MetrizedNet& net) {
    const DerivedOnceNet& A = net.A;
    const Dart root = A.dart(0);

    vector<pt> res(A.V());
    for (Dart d: A.darts()) {
        ApproximateIsometry q = net.getIsometry(root, d);
        res[d.from] = make_deterministic(q({1, 0, 0}));
    }

    return res;
}

// Adds a pin.
// If connected_to = {d}, then phi(d) will be the dart from d.to to the pin.
// If connected_to = {d, e}, then
//   phi(d) = z_1 and phi(z_1) = z_2, phi(z_2) = e,
// where z_1, z_2 are the darts (d.to->pin) and (pin->e.from).
int addPin(CombinatorialMap& A, vector<Dart> connected_to) {
    assert(1 <= ssize(connected_to) and ssize(connected_to) <= 2);
    Dart d1 = A.addVertex(connected_to[0]).first;
    if (connected_to.size() == 2)
        A.addEdge(d1, connected_to[1]);

    return d1.to;
}

vector<Dart> getFaceDarts(const CombinatorialMap& A, vector<int> vertices) {
    for (Dart d: A.darts()) {
        if (d.from == vertices[0] and d.to == vertices[1]) {
            bool works = true;
            for (int i = 0; i < ssize(vertices); i++) {
                works &= A.phi(d, i).from == vertices[i];
            }
            works &= A.phi(d, ssize(vertices)).from == vertices[0];
            if (works)
                return A.face(d);
        }
    }
    assert(0);
    return {};
}

vector<Dart> getConcaveAngles(const DerivedOnceNet& A,
                              const vector<Dart>& face) {
    vector<Dart> res;
    for (Dart d: face)
        if (A.isConcaveAngle(d))
            res.push_back(d);
    return res;
}

// It returns a "heuristic" of the cost of putting a pin (with given
// coordinates) in a given face.
// Returns the average of the (Euclidean) distance from the pin to the
// concave angles.
// It is a naive heuristic used to determine the (concave) face a pin belongs to.
double pinAverageDistance(const DerivedOnceNet& A, vector<pt> embedding,
                   pt pin, const vector<Dart>& face) {
    vector<Dart> concave_angles = getConcaveAngles(A, face);
    assert(!concave_angles.empty());

    double res = 0;
    for (Dart concave_angle: concave_angles) {
        pt P = embedding[concave_angle.to];
        res += pow(norm(P-pin), 2);
    }
    return sqrt(res) / ssize(concave_angles);
}



void saveSurfaceEvolverFormat(const vector<pt>& coords,
                              int pin_id,
                              const vector<bool>& is_vertex_fixed,
                              const vector<vector<int>>& faces,
                              string out_file) {
    pt line = coords[pin_id];

    map<pair<int,int>, int> edge2id;
    auto get_edge = [&](int a, int b) {
        assert(a != b);
        assert(0 <= a and a < ssize(coords));
        assert(0 <= b and b < ssize(coords));

        if (!edge2id.count({a, b})) {
            int id = 1 + edge2id.size() / 2;
            edge2id[{a, b}] = id;
            edge2id[{b, a}] = -id;
        }
        return edge2id[{a, b}];
    };

    auto is_vertex_on_line = [&](int v) {
        return norm(coords[v] ^ line) < 0.001;
    };

    auto find_pivot_index = [&](vector<int> face) {
        for (int i = 0; i < ssize(face); i++)
            if (is_vertex_on_line(face[i]))
                return i;
        return 0;
    };

    vector<array<int, 3>> triangles;

    auto triangulate = [&](vector<int> face) {
        assert(ssize(face) >= 3);
        int n = ssize(face);
        int p = find_pivot_index(face);
        for (int i = 0; i < ssize(face); i++) {
            int j = (i+1) % n;
            if (i == p or j == p)
                continue;
            triangles.push_back({
                get_edge(face[p], face[i]),
                get_edge(face[i], face[j]),
                get_edge(face[j], face[p])
            });
        }
    };

    for (vector<int> face: faces)
        triangulate(face);

    ofstream out(out_file);

    out << "boundary 1 parameters 1" << endl;
    out << "x1:  " << line[0] << " * p1" << endl;
    out << "x2:  " << line[1] << " * p1" << endl;
    out << "x3:  " << line[2] << " * p1" << endl;
    out << endl;

    out << "vertices" << endl;
    for (int i = 0; i < ssize(coords); i++) {
        array<double, 3> P = coords[i];
        bool fixed = is_vertex_fixed[i];
        bool boundary = is_vertex_on_line(i);
        boundary &= !fixed;
        out << i+1;
        if (boundary) {
            out << " " << (P * line) << " boundary 1";
        } else {
            out << " " << P[0] << " " << P[1] << " " << P[2];
            if (fixed)
                out << " fixed";
        }
        out << endl;
    }
    out << endl;

    auto is_edge_fixed = [&](pair<int,int> edge) {
        int u = edge.first, v = edge.second;
        return is_vertex_fixed[u] and is_vertex_fixed[v];
    };

    auto is_edge_boundary = [&](pair<int,int> edge) {
        int u = edge.first, v = edge.second;
        return is_vertex_on_line(u) and is_vertex_on_line(v);
    };

    out << "edges" << endl;
    for (auto pp: edge2id) {
        if (pp.second < 0)
            continue;
        out << pp.second << " "
            << pp.first.first+1 << " " << pp.first.second+1;
        if (is_edge_fixed(pp.first))
            out << " fixed";
        if (is_edge_boundary(pp.first))
            out << " boundary 1";
        out << endl;
    }
    out << endl;

    out << "faces" << endl;
    for (int i = 0; i < ssize(triangles); i++) {
        array<int, 3> T = triangles[i];
        out << i+1 << " "
            << T[0] << " " << T[1] << " " << T[2] << endl;
    }
    out << endl;

    out << R"(read
showq;
subroutine := {
    g 2000;
    areaweed 0.0005;
    u 4;
    print total_area;
}

gogo := {
    quiet;
    r 2;
    subroutine 10;
}
print total_area;
)";

    out.close();
}



const set<int> possible_nets = {
    72,     278,    434,    521,    691,
    936,    949,    1011,   1034,   1223,
    1351,   1405,   1514,   1829,   2017,
    2119,   2192,   2312,   2454,   2676,
    2700,   2769,   2799,   2891,   3048,
    3222,   3233,   3295,   3338};

// AD-HOC
// Here we describe the "pop".
// Use python run_visualization net_id as well as the faces' list
// (printed by this program) to come up with a reasonable choice.
// This choice, which is crucial, can be informed by experiments.
// The description of one pop consists of:
//   - a list of the faces to remove.
//   - a path to add (for each edge in the path, we add the triangle that it
//   forms with the origin).
const int PIN_ZERO = 342905;
const int PIN_ONE = 5896232;
const map<int, pair<vector<int>, vector<int>>> pop_description = {
    {1405, {
        {2, 4},
        {PIN_ONE, 3, 1, 5, PIN_ZERO, 4, 2, PIN_ONE}
    }},
    {2769, {
        {1, 5},
        {PIN_ZERO, 8, 1, 2, 11, PIN_ONE, 10, 7, 6, 9, PIN_ZERO}
    }},
    {936, {
        {5, 2},
        {PIN_ZERO, 6, 1, 2, 8, PIN_ONE, 9, 5, 7, PIN_ZERO}
    }},
    {1034, {
        {7, 1},
        {PIN_ONE, 10, 5, 1, 9, PIN_ZERO, 8, 4, 7, 11, PIN_ONE}
    }},
    {1223, {
        {2, 8},
        {PIN_ONE, 13, 9, 0, 10, PIN_ZERO, 11, 3, 7, 12, PIN_ONE}
    }},
    {1514, {
        {6},
        {PIN_ZERO, 12, 1, 0, 13, PIN_ONE}
    }},
    {2119, {
        {6, 7},
        {PIN_ZERO, 9, 5, 7, 10, PIN_ONE, 11, 2, 1, 8, PIN_ZERO}
    }},
    {2676, {
        {1},
        {PIN_ZERO, 8, 1, 2, 10, PIN_ONE}
    }},
    {3295, {
        {0},
        {PIN_ZERO, 10, 1, 5, 11, PIN_ONE}
    }},
    {521, {
        {1, 11},
        {PIN_ZERO, 21, 4, 3, 2, 18, PIN_ONE, 19, 6, 9, 13, 20, PIN_ZERO}
    }},
    {1351, {
        {3, 5},
        {PIN_ZERO, 12, 1, 0, 14, PIN_ONE, 15, 7, 6, 13, PIN_ZERO}
    }},
    {1829, {
        {1, 7},
        {PIN_ZERO, 12, 5, 6, 15, PIN_ONE, 14, 10, 11, 9, 13, PIN_ZERO}
    }},
    {691, {
        {6, 8},
        {PIN_ZERO, 7, 6, 13, PIN_ONE, 0, 1, 4, 11, PIN_ZERO}
    }}
};

int main() {
    const int net_id = 691;
    assert(possible_nets.count(net_id));

    const DerivedOnceNet A2 =
            parseDerivedOnceNets("data/derived_once_computer.txt")[net_id];
    const MetrizedNet net = parseMetrizedNet(A2,
            "data/metrized_nets/"+ to_string(net_id) + ".txt");
    net.incrementVersion();
    vector<pt> embedding = getVertexEmbedding(net);

    // We locate the two pins metrically.
    // Their coordinates are stored in pins_coord.
    // They both belong to the unit sphere and they are antipodal.
    const Dart root = A2.dart(0);
    vector<pt> equators_with_pins;
    vector<Dart> concave_angles = A2.concaveAngles();
    for (Dart d: A2.concaveAngles())
        equators_with_pins.push_back(
                make_deterministic(net.getIsometryToPin(root, d)({0, 0, 1})));
    vector<pt> pins_coord;
    for (pt eq1: equators_with_pins)
            for (pt eq2: equators_with_pins) {
        pt z = eq1 ^ eq2;
        if (norm(z) > 0.1) {
            z = z / norm(z);
            pins_coord = {z, -z};
            break;
        }
    }

    // We add the pins to the DerivedOnceNet.
    // Doing so, we obtain a CombinatorialMap (as the result is not
    // DerivedOnce anymore) and we save the indexes of the pins in pins_id.
    vector<vector<Dart>> concave_faces;
    for (auto face: A2.faces()) {
        if (!A2.isConvexFace(face))
            concave_faces.push_back(face);
    }

    assert(ssize(pins_coord) == 2);
    assert(ssize(concave_faces) <= 2);

    CombinatorialMap A = A2;
    vector<int> pins_id(2);
    for (int i = 0; i < 2; i++) {
        pt pin = pins_coord[i];
        bool pin_is_vertex = false;
        for (int v = 0; v < ssize(embedding); v++) {
            if (norm(embedding[v] - pin) < 0.001) {
                pins_id[i] = v;
                pin_is_vertex = true;
                break;
            }
        }
        if (pin_is_vertex and A2.degree(pins_id[i]) == 3)
                continue;

        assert(!concave_faces.empty());

        vector<Dart> face = concave_faces[0];
        if (concave_faces.size() == 2) {
            double cost0 = pinAverageDistance(A2, embedding, pin, concave_faces[0]);
            double cost1 = pinAverageDistance(A2, embedding, pin, concave_faces[1]);

            if (cost1 < cost0)
                face = concave_faces[1];
        }
        vector<Dart> concave_angles = getConcaveAngles(A2, face);
        if (pin_is_vertex) {
            if (ssize(concave_angles) == 1) {
                assert(concave_angles[0].to == pins_id[i]);
                continue;
            }
            assert(ssize(concave_angles) == 2);
            if (concave_angles[0].to != pins_id[i])
                swap(concave_angles[0], concave_angles[1]);
            A.addEdge(concave_angles[0], A.phi(concave_angles[1]));
            continue;
        }
        if (ssize(concave_angles) == 2)
            concave_angles[1] = A2.phi(concave_angles[1]);
        int pin_id = addPin(A, concave_angles);
        assert(ssize(embedding) == pin_id);
        embedding.push_back(pins_coord[i]);
        pins_id[i] = pin_id;
    }

    // We build the cone *outside* of the ball of radius 1/3.
    int n = ssize(embedding);
    for (int i = 0; i < n; i++) {
        pt P = embedding[i];
        embedding.push_back({P[0]/3, P[1]/3, P[2]/3});
    }
    int origin = 2*n;
    embedding.push_back({0, 0, 0});

    vector<vector<int>> faces_before_pop;
    for (Dart d: A.darts()) {
        assert(d.from != d.to);
        if (d.from > d.to)
            continue;
        faces_before_pop.push_back({d.from, d.to, origin});
    }

    vector<vector<int>> faces_after_pop;
    for (Dart d: A.darts()) {
        assert(d.from != d.to);
        if (d.from > d.to)
            continue;
        faces_after_pop.push_back({d.from, d.to, d.to+n, d.from+n});
    }

    // We print some info which is useful for humans to come up with
    // a reasonable pop.
    auto faces = A.faces();
    cout << "Faces: " << endl;
    for (int i = 0; i < ssize(faces); i++) {
        cout << "Face " << i << ": ";
        for (Dart d: faces[i])
            cout << d.from << ", ";
        cout << endl;
    }

    cout << "Pin ids (in order): " << pins_id[0] << ", " << pins_id[1] << endl;

    if (!pop_description.contains(net_id)) {
        cout << "The description of the pop for this net (" << net_id << ") "
             << "is missing." << endl;
        return 0;
    }

    // We build the pop (inside the ball of radius 1/3)
    assert(pop_description.contains(net_id));
    const vector<int>& excluded_faces = pop_description.at(net_id).first;
    vector<int> path_to_add = pop_description.at(net_id).second;
    for (int& x: path_to_add) {
        if (x == PIN_ZERO)
            x = pins_id[0];
        if (x == PIN_ONE)
            x = pins_id[1];
    }

    for (int i = 0; i < ssize(faces); i++) {
        bool is_excluded = false;
        for (int excl: excluded_faces)
            is_excluded |= i == excl;
        if (is_excluded)
            continue;
        vector<int> half_face;
        for (Dart d: faces[i])
            half_face.push_back(d.from + n);
        faces_after_pop.push_back(half_face);
    }

    for (int i = 0; i < ssize(path_to_add)-1; i++) {
        int u = path_to_add[i], v = path_to_add[i+1];
        bool existing_dart = false;
        for (Dart d: A.darts())
            existing_dart |= d.from == u and d.to == v;
        assert(existing_dart);
        faces_after_pop.push_back(
            {u+n, v+n, origin});
    }

    // We print the surfaces we have built.
    vector<bool> is_vertex_fixed(2*n+1, false);
    for (int i = 0; i < n; i++)
        is_vertex_fixed[i] = true;

    string filename_prefix = "data/surface_evolver/" + to_string(net_id);
    saveSurfaceEvolverFormat(embedding, pins_id[0], is_vertex_fixed,
                             faces_before_pop,
                             filename_prefix + "_before_pop.fe");
    saveSurfaceEvolverFormat(embedding, pins_id[0], is_vertex_fixed,
                             faces_after_pop,
                             filename_prefix + "_after_pop.fe");

}

// 5.7821
