#ifndef GENERATE_SURFACE_EVOLVER_INPUT
#define GENERATE_SURFACE_EVOLVER_INPUT


#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "metrized_nets.cpp"
#include "metrization_constraints.cpp"
#include "simplicial_complexes.cpp"


string generateEmbeddingFilename(int net_id, vector<double> params) {
    std::ostringstream oss;
    oss << net_id;
    oss << fixed << setprecision(5);
    for (double param: params)
        oss << "_" << param;
    return oss.str();
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
// It is a naive heuristic used to determine the (non-convex) face a pin belongs to.
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

// We locate the two pins metrically.
//
// The pins are antipodal points on the unit sphere.
pair<pt, pt> findPinsEmbedding(const MetrizedNet& net) {
    const Dart root = net.A.dart(0);
    vector<pt> equators_with_pins;
    vector<Dart> concave_angles = net.A.concaveAngles();
    for (Dart d: net.A.concaveAngles())
        equators_with_pins.push_back(
                make_deterministic(net.getIsometryToPin(root, d)({0, 0, 1})));
    pair<pt, pt> pins;
    bool found = false;
    for (pt eq1: equators_with_pins)
            for (pt eq2: equators_with_pins) {
        pt z = eq1 ^ eq2;
        if (norm(z) > 0.1) {
            z = z / norm(z);
            pins = {z, -z};
            found = true;
            break;
        }
    }
    assert(found);

    return pins;
}


// This function is for the nets such that the associated once-derived net
// is metrically rigid but the pins are not rigid.
//
// dist denotes the length of the edge connecting one pin to
// net.A.concaveAngles()[0].
//
// The pins are antipodal points on the unit sphere.
pair<pt, pt> choosePinsEmbedding(const MetrizedNet& net, double dist) {
    assert(!net.A.concaveAngles().empty());

    const Dart root = net.A.dart(0);
    Dart concave_angle = net.A.concaveAngles()[0];
    vector<PathStep> path_to_pin = net.getPathToPin(root, concave_angle);

    path_to_pin.push_back(PathStep('E', dist));

    pt pin = make_deterministic(net.isometryFromPath(path_to_pin)({1, 0, 0}));

    return {pin, -pin};
}

// We add the pins to the DerivedOnceNet.
// Doing so, we obtain a CombinatorialMap (as the result is not
// DerivedOnce anymore) and we save the indexes of the pins in pins_id.
pair<pair<int,int>, CombinatorialMap> addPinsToCombinatorialMap(
        const MetrizedNet& net,
        vector<pt>& embedding, const pair<pt, pt>& pins_coord) {
    vector<vector<Dart>> concave_faces;
    for (auto face: net.A.faces()) {
        if (!net.A.isConvexFace(face))
            concave_faces.push_back(face);
    }
    assert(ssize(concave_faces) == 2);

    auto find_vertex = [&](pt p) {
        for (int v = 0; v < ssize(embedding); v++) {
            if (norm(p-embedding[v]) < 1e-4)
                return v;
        }
        return -1;
    };

    CombinatorialMap A = net.A;
    vector<int> pins_vec;
    for (const pt& pin: {pins_coord.first, pins_coord.second}) {
        assert(!concave_faces.empty());

        double cost0 = pinAverageDistance(net.A, embedding, pin, concave_faces[0]);
        double cost1 = pinAverageDistance(net.A, embedding, pin, concave_faces[1]);
        vector<Dart> face = concave_faces[cost1 < cost0];

        vector<Dart> concave_angles = getConcaveAngles(net.A, face);
        if (ssize(concave_angles) == 2)
            concave_angles[1] = net.A.phi(concave_angles[1]);

        if (find_vertex(pin) != -1) {
            A.addEdge(concave_angles[0], concave_angles[1]);
            pins_vec.push_back(find_vertex(pin));
        } else {
            int pin_id = addPin(A, concave_angles);
            assert(ssize(embedding) == pin_id);
            embedding.push_back(pin);
            pins_vec.push_back(pin_id);
        }
    }

    return {{pins_vec[0], pins_vec[1]}, A};
}

struct SurfaceEvolverData {
    ImmersedSimplicialComplex X;
    pair<int,int> pins;
    vector<pair<int,int>> boundary_edges;
};

void saveSurfaceEvolverFormat(int net_id,
                              string embedding_name,
                              const SurfaceEvolverData& data,
                              string out_file,
                              bool print_commands) {
    const auto& X = data.X;
    pt line = X.getPoint(data.pins.first);

    map<pair<int,int>, int> edge2id;
    auto get_edge = [&](int a, int b) {
        assert(a != b);
        assert(0 <= a and a < X.V());
        assert(0 <= b and b < X.V());

        if (!edge2id.count({a, b})) {
            int id = 1 + edge2id.size() / 2;
            edge2id[{a, b}] = id;
            edge2id[{b, a}] = -id;
        }
        return edge2id[{a, b}];
    };

    auto is_vertex_on_line = [&](int v) {
        return norm(X.getPoint(v) ^ line) < 0.001;
    };

    vector<array<int, 3>> triangles;
    for (auto triangle: X.getFaces(2)) {
        vector<int> triangle_vec(triangle.begin(), triangle.end());
        auto A = triangle_vec[0], B = triangle_vec[1], C = triangle_vec[2];
        triangles.push_back({get_edge(A, B), get_edge(B, C), get_edge(C, A)});
    }
    vector<bool> is_vertex_fixed(X.V(), false);
    for (auto edge: data.boundary_edges) {
        is_vertex_fixed[edge.first] = true;
        is_vertex_fixed[edge.second] = true;
    }

    ofstream out(out_file);

    out << "boundary 1 parameters 1" << endl;
    out << "x1:  " << line[0] << " * p1" << endl;
    out << "x2:  " << line[1] << " * p1" << endl;
    out << "x3:  " << line[2] << " * p1" << endl;
    out << endl;

    out << "vertices" << endl;
    for (int i = 0; i < X.V(); i++) {
        pt P = X.getPoint(i);
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
print total_area;
)";

    const string generic_surface_evolver_commands = R"(
    r 2;
    subroutine 5;
    o;
    subroutine 3;
    o;
    subroutine 3;
    r;
    u 4;
    subroutine 5;
    o;
    subroutine 3;)";

    const map<int, string> special_surface_evolver_commands = {
        {28, R"(
            o;
            subroutine 5;
        )"},
        {99, R"(
            subroutine 2;
            o;
            subroutine 2;
            o;
            subroutine 2;
        )"},
        {346, R"(
            subroutine 2;
            o;
            subroutine 2;
            o;
            subroutine 2;
        )"},
        {421, R"(
            subroutine 3;
            r;
            subroutine 5;
            r;
            subroutine 3;
            o;
            subroutine 8;
            o;
            subroutine 3;
            o;
            subroutine 3;
        )"},
        // TODO: It's not necessary anymore.
        //       We leave it here for some time. To remove eventually.
        // {468, R"(
        //     subroutine 2;
        //     r;
        //     subroutine 5;
        //     o;
        //     subroutine 5;
        //     r;
        //     subroutine 5;
        // )"},
        {614, R"(
            subroutine 3;
            r;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
            r;
            subroutine 3;
        )"},
        {736, R"(
            r 2;
            subroutine 12;
        )"},
        {779, R"(
            subroutine 2;
            o;
            subroutine 2;
            o;
            subroutine 2;
        )"},
        {827, R"(
            g 5;
            r;
            g 400;
            u 2;
            g 300;
        )"},
        {1029, R"(
            subroutine 3;
            r;
            subroutine 2;
        )"},
        {1096, R"(
            r; subroutine 5; o; subroutine 5; o; subroutine 5;
        )"},
    };

    string surface_evolver_commands = generic_surface_evolver_commands;
    if (special_surface_evolver_commands.count(net_id))
        surface_evolver_commands = special_surface_evolver_commands.at(net_id);


    constexpr string_view format_str = R"(read
subroutine := {{
    quiet;
    g 2000;
    areaweed 0.0005;
    u 4;
    print total_area;
    quiet off;
}}

subroutinequick := {{
    quiet;
    g 200;
    areaweed 0.0005;
    u 4;
    print total_area;
    quiet off;
}}

gogo := {{
    /*showq;*/
{}
    dump "data/surface_evolver/{}_competitor.out";
    /*close_show;*/ /* This raises an error if the window is not visible. */
}}

gogo;

quit;
)";
    if (print_commands)
        out << format(format_str, surface_evolver_commands, embedding_name);

    out.close();
}



// This is one of the two 1-parameter families corresponding to 1029.
pair<vector<pt>, pair<pt,pt>> embedGoggles(double q) {
    auto rotate = [&](pt U) {
        return pt{U[0]*cos(q)-U[1]*sin(q), U[1]*cos(q)+U[0]*sin(q), U[2]};
    };
    auto sym_x = [&](pt U) { return pt{-U[0], U[1], U[2]}; };

    pt pin0 = {0, 0, 1};
    pt pin1 = -pin0;

    vector<pt> P(8);
    // Face of the cube
    P[1] = pt{1, sqrt(2), 0};
    P[4] = pt{1, 0, sqrt(2)};
    P[6] = pt{1, 0, -sqrt(2)};
    P[0] = pt{1, -sqrt(2), 0};

    for (int i: {1, 4, 6, 0})
        P[i] = rotate(P[i]);

    P[5] = sym_x(P[4]);
    P[3] = sym_x(P[0]);
    P[7] = sym_x(P[6]);
    P[2] = sym_x(P[1]);

    pin0 = pin0 / norm(pin0);
    pin1 = pin1 / norm(pin1);

    for (pt& U: P)
        U = U / norm(U);

    return {P, {pin0, pin1}};
}

// This is one of the two 1-parameter families corresponding to 1029.
pair<vector<pt>, pair<pt,pt>> embedCube(double q) {
    auto sym_x = [&](pt U) {
        return pt{-U[0], U[1], U[2]};
    };

    pt pin0 = {sin(q), 0, cos(q)};
    pt pin1 = -pin0;

    vector<pt> P(8);
    // Face of the cube
    P[1] = pt{1, sqrt(2), 0};
    P[4] = pt{1, 0, sqrt(2)};
    P[6] = pt{1, 0, -sqrt(2)};
    P[0] = pt{1, -sqrt(2), 0};

    P[5] = sym_x(P[4]);
    P[3] = sym_x(P[0]);
    P[7] = sym_x(P[6]);
    P[2] = sym_x(P[1]);

    pin0 = pin0 / norm(pin0);
    pin1 = pin1 / norm(pin1);

    for (pt& U: P)
        U = U / norm(U);

    return {P, {pin0, pin1}};
}


// TODO: This should be moved somewhere else?
const set<int> possible_nets = {
    28,     99,    158,    182,    250,
    319,    325,    346,   361,   421,
    468,   484,   520,   614,   681,
    712,   736,   779,   827,   909,
    920,   938,    1029,
    1096,   1100,   1122};


const set<int> rigid_nets_flexible_pins = {28, 99, 346, 779, 1096};

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
    {28, {
        {3, 5},
        {PIN_ZERO, 17, 6, 14, 13, 18, PIN_ONE, 19, 9, 0, 1, 16, PIN_ZERO}
    }},
    {99, {
        {1},
        {PIN_ZERO, 14, 7, 13, 11, 10, 15, PIN_ONE}
    }},
    {182, {
        {1, 11},
        {PIN_ZERO, 21, 4, 3, 2, 18, PIN_ONE, 19, 6, 9, 13, 20, PIN_ZERO}
    }},
    {250, {
        {6, 8},
        {PIN_ZERO, 7, 6, 13, PIN_ONE, 0, 1, 4, 11, PIN_ZERO}
    }},
    {319, {
        {5, 2},
        {PIN_ZERO, 6, 1, 2, 8, PIN_ONE, 9, 5, 7, PIN_ZERO}
    }},
    {325, {
        {3, 1, 9},
        {PIN_ZERO, 4, 3, 2, 6, PIN_ONE, 19, 9, 12, 17 , PIN_ZERO, 14, 13, 7, 5, PIN_ONE}
    }},
    {346, {
        {1},
        {PIN_ZERO, 12, 0, 4, 3, 13, PIN_ONE}
    }},
    {361, {
        {7, 1},
        {PIN_ONE, 10, 5, 1, 9, PIN_ZERO, 8, 4, 7, 11, PIN_ONE}
    }},
    {421, {
        {2, 8},
        {PIN_ONE, 13, 9, 0, 10, PIN_ZERO, 11, 3, 7, 12, PIN_ONE}
    }},
    {468, {
        {3, 5},
        {PIN_ZERO, 12, 1, 0, 14, PIN_ONE, 15, 7, 6, 13, PIN_ZERO}
    }},
    {484, {
        {2, 4},
        {PIN_ONE, 3, 1, 5, PIN_ZERO, 4, 2, PIN_ONE}
    }},
    {520, {
        {6},
        {PIN_ZERO, 12, 1, 0, 13, PIN_ONE}
    }},
    {614, {
        {1, 7},
        {PIN_ZERO, 12, 5, 6, 15, PIN_ONE, 14, 10, 11, 9, 13, PIN_ZERO}
    }},
    {712, {
        {6, 7},
        {PIN_ZERO, 9, 5, 7, 10, PIN_ONE, 11, 2, 1, 8, PIN_ZERO}
    }},
    {736, {
        {0, 6, 8},
        {PIN_ZERO, 1, 2, 3, PIN_ONE, 13, 6, 10, PIN_ZERO, 4, 8, 9, PIN_ONE}
    }},
    {779, {
        {1},
        {PIN_ZERO, 10, 1, 2, 3, 11, PIN_ONE}
    }},
    {827, {
        {3},
        {PIN_ZERO, 0, 1, 2, PIN_ONE}
    }},
    {909, {
        {1},
        {PIN_ZERO, 8, 1, 2, 10, PIN_ONE}
    }},
    {938, {
        {1, 5},
        {PIN_ZERO, 8, 1, 2, 11, PIN_ONE, 10, 7, 6, 9, PIN_ZERO}
    }},
    {1029, {
        {5, 0},
        {PIN_ZERO, 5, 2, 7, PIN_ONE,6, 0, 4, PIN_ZERO}
    }},
    {1096, {
        {6},
        {PIN_ZERO, 16, 2, 3, 10, 17, PIN_ONE}
    }},
    {1122, {
        {0},
        {PIN_ZERO, 10, 1, 5, 11, PIN_ONE}
    }},
};

struct NetEmbeddedWithPins {
    CombinatorialMap A;
    vector<pt> embedding;
    pair<int,int> pins;

    double totalLen() {
        double len = 0;
        for (Dart d: A.darts()) {
            int x = d.from, y = d.to;
            if (x > y)
                continue;
            len += norm(embedding[x] - embedding[y]);
        }
        return len;
    }
    double totalInterestingLen() {
        auto f = [&](int x, int y) { return embedding[x] - embedding[y]; };
        pt line = f(pins.first, pins.second);
        pair<pt, pt> edge_to_pins;
        for (Dart d: A.darts()) {
            if (d.from == pins.first)
                edge_to_pins.first = f(d.from, d.to);
            if (d.from == pins.second)
                edge_to_pins.second = f(d.from, d.to);
        }
        pt edge_to_pin = edge_to_pins.first;
        if (norm(edge_to_pins.second) > norm(edge_to_pins.first))
            edge_to_pin = edge_to_pins.second;
        pt rotation_axis = edge_to_pin ^ line;
        rotation_axis = rotation_axis / norm(rotation_axis);

        double len = 0;
        for (Dart d: A.darts()) {
            int x = d.from, y = d.to;
            if (x > y)
                continue;
            double a = embedding[x] * rotation_axis, b = embedding[y] * rotation_axis;
            if (a > b) swap(a, b);
            len += b-a;
        }
        return len;
    }
};

// net is a rigid net with flexible pins.
// Let I be the interval of lengths of the edge connecting concave_angles[0] to
// its pins that generate a valid net.
// This function returns I.
Interval validPinEdgeLength(int net_id, MetrizedNet net) {
    assert(rigid_nets_flexible_pins.count(net_id));
    // All arcs connecting to pins belong to the same equator.
    vector<double> lens = maxDistanceToPin(net);
    vector<Dart> concave_angles = net.A.concaveAngles();
    int other = -1;
    for (int i = 1; i < ssize(concave_angles); i++) {
        if (!net.A.sameFace(concave_angles[0], concave_angles[i]))
            other = i;
    }
    assert(other != -1);
    // concave_angles[0] and concave_angles[other] belong to different faces.
    Dart d1 = concave_angles[0], d2 = concave_angles[other];
    double len1 = lens[0], len2 = lens[other];
    auto iso = net.getIsometryFromPinToPin(d1, d2);

    assert(len1 > 0.01 and len2 > 0.01);
    assert(len1 < M_PI and len2 < M_PI);

    double pos_x = iso({1, 0, 0})[0].mid();
    double pos_y = iso({1, 0, 0})[1].mid();
    double dist = atan2(pos_y, pos_x);
    if (dist < 0) dist += 2 * M_PI;

    double almost_sign = iso({0, 0, 1})[2].mid();
    assert(abs(almost_sign) > 0.999);
    bool same_dir = almost_sign > 0;

    double l, r;
    if (same_dir) l = dist, r = dist + len2;
    else l = dist-len2, r = dist;
    l -= M_PI;
    r -= M_PI;

    assert(l <= r);
    assert(0 <= r and r <= M_PI);

    return intersection({0, len1}, {l, r});
}

// We print some info which is useful for humans to come up with a reasonable pop.
void printInfoForHumanPopping(const NetEmbeddedWithPins& link) {
    auto faces = link.A.faces();
    cout << "Faces: " << endl;
    for (int i = 0; i < ssize(faces); i++) {
        cout << "Face " << i << ": ";
        for (Dart d: faces[i])
            cout << d.from << ", ";
        cout << endl;
    }
    cout << "Pin ids (in order): " << link.pins.first << ", "
                                   << link.pins.second << endl;
}

// Fat fish
NetEmbeddedWithPins embed827WithPins(const vector<double>& params) {
    assert(ssize(params) == 2);

    // Checking the validity of the parameters.
    // theta0 represents the position of the pins (relative to the square).
    // theta1 represents the position of the biangle (relative to the square).
    // (theta0, theta1) =  (0, 0) corresponds to the most symmetric configuration.
    double theta0 = params[0], theta1 = params[1];
    double c0 = cos(theta0), s0 = sin(theta0);
    double c1 = cos(theta1), s1 = sin(theta1);

    const double sigma = sqrt(2. / 3);
    const double tau = 1. / sqrt(3.);

    assert(0 <= theta1 and theta1 <= acos(sigma));
    // The condition on the tangent follows from computing
    // det(embedding[5], embedding[7], embedding[9])
    double upperbound_theta0 = atan(1.0 / (sqrt(3.) * cos(theta1)));
    assert(0 <= theta0 and theta0 <= upperbound_theta0);


    // We construct the combinatorial map from scratch (instead of starting from
    // the DerivedOnceNet) because we have to add the pins and also to add
    // some artificial vertices on the biangle's edges.
    CombinatorialMap A = constructOneEdge();

    Dart d01 = A.darts()[0];
    Dart d12 = A.addVertex(d01).first;
    Dart d23 = A.addVertex(d12).first;
    Dart d30 = A.addEdge(d23, d01).first;
    Dart d34 = A.addVertex(A.rev(d30)).first;
    Dart d15 = A.addVertex(A.rev(d12)).first;
    Dart d56 = A.addVertex(d15).first;
    Dart d47 = A.addVertex(d34).first;
    A.addEdge(d56, d47).first; // d64
    A.addEdge(d47, d56).first; // d75
    A.addVertex(A.rev(d01)).first; // d08
    A.addVertex(A.rev(d23)).first; // d29

    pair<int,int> pins = {8, 9};

    vector<pt> embedding(A.V());
    embedding[0] = {tau, 0, sigma};
    embedding[1] = {tau, sigma, 0};
    embedding[2] = {tau, 0, -sigma};
    embedding[3] = {tau, -sigma, 0};
    embedding[4] = {s1, -c1, 0};
    embedding[5] = {-s1, c1, 0};
    embedding[6] = {-c1/2, -s1/2, sqrt(3)/2};
    embedding[7] = {-c1/2, -s1/2, -sqrt(3)/2};
    embedding[8] = {s0, 0, c0}; // TODO: Check that when this is equal to embedding[0], also another condition is saturated.
    embedding[9] = {-s0, 0, -c0};

    return {
        .A = A,
        .embedding = embedding,
        .pins = pins,
    };
}


NetEmbeddedWithPins embedNetWithPins(int net_id, const vector<double>& params) {
    assert(possible_nets.count(net_id));

    if (net_id == 827)
        return embed827WithPins(params);

    ifstream metrized_net_in("data/metrized_nets/"+ to_string(net_id) + ".txt");
    assert(metrized_net_in.is_open());
    MetrizedNet net = parseMetrizedNet(metrized_net_in);
    metrized_net_in.close();

    vector<pt> embedding;
    pair<pt, pt> pins_coord;

    if (net_id == 1029) {
        assert(ssize(params) == 2);
        assert(params[0] == 0 or params[0] == 1);
        if (params[0] == 0) {
            assert(0 <= params[1] and params[1] <= acos(1.0/3) / 2);
            tie(embedding, pins_coord) = embedGoggles(params[1]);
        }
        if (params[0] == 1) {
            assert(0 <= params[1] and params[1] <= acos(1.0/3) / 2); // TODO: Check this!
            tie(embedding, pins_coord) = embedCube(params[1]);
        }
    } else if (rigid_nets_flexible_pins.count(net_id)) {
        assert(ssize(params) == 1);
        double len = params[0];
        Interval valid_lens = validPinEdgeLength(net_id, net);
        assert(valid_lens.contains(len));
        // clog << "Interval of valid edge-to-pin lengths: " << valid_lens << endl;
        embedding = getVertexEmbedding(net);
        pins_coord = choosePinsEmbedding(net, len);
    } else if (net_id == 736) {
        assert(params.empty());
        embedding = getVertexEmbedding(net);
        pins_coord = {embedding[11], embedding[12]};
    } else {
        assert(params.empty());
        embedding = getVertexEmbedding(net);
        pins_coord = findPinsEmbedding(net);
    }

    auto [pins, A] = addPinsToCombinatorialMap(net, embedding, pins_coord);

    return {
        .A = A,
        .embedding = embedding,
        .pins = pins,
    };
}



SurfaceEvolverData generateCone(const NetEmbeddedWithPins& link) {
    auto embedding = link.embedding;
    int n = ssize(embedding);
    int origin = n;
    embedding.push_back({0, 0, 0});

    vector<pair<int,int>> boundary_edges;
    vector<set<int>> triangles;
    for (Dart d: link.A.darts()) {
        assert(d.from != d.to);
        if (d.from > d.to)
            continue;
        triangles.push_back({d.from, d.to, origin});
        boundary_edges.push_back({d.from, d.to});
    }

    return {.X = ImmersedSimplicialComplex(embedding, triangles),
            .pins = link.pins,
            .boundary_edges = boundary_edges};
}



SurfaceEvolverData generateCubeCompetitor(const NetEmbeddedWithPins& link) {
    vector<int> up_face = {4, 1, 6, 0};
    vector<int> down_face = {5, 2, 7, 3};
    vector<int> small_square = {10, 11, 12, 13};
    auto embedding = link.embedding;

    // up_face[0] - down_face[0] contains link.pins.first
    pt A = embedding[link.pins.first];
    pt B = embedding[up_face[1]] + embedding[down_face[1]];
    pt C = -A;
    pt D = -B;
    vector<pt> small_square_pts = {A, B, C, D};
    for (pt& P: small_square_pts)
        P = P / (5 * norm(P));
    for (pt P: small_square_pts)
        embedding.push_back(P);


    vector<pair<int,int>> boundary_edges;
    for (Dart d: link.A.darts()) {
        assert(d.from != d.to);
        if (d.from > d.to)
            continue;
        boundary_edges.push_back({d.from, d.to});
    }

    vector<set<int>> triangles;

    // Whether v belongs to the line connecting the pins.
    auto is_pin = [&](int v) {
        return set<int>{link.pins.first, link.pins.second, 10, 12}.count(v);
    };

    auto add_face = [&](vector<int> face) {
        assert(ssize(face) >= 3);
        int f = ssize(face);
        // If a face contains a pin, then the triangulation shall be pivoted
        // there to avoid issues due to the fact that the face is not convex
        // and not contained in a plane.
        int pivot = 0;
        for (int i = 1; i < f; i++)
            if (is_pin(face[i]))
                pivot = i;
        rotate(face.begin(), face.begin() + pivot, face.end());
        for (int i = 1; i < f-1; i++) {
            int j = (i+1) % f;
            triangles.push_back({face[0], face[i], face[j]});
        }
    };

    for (int i = 0; i < 4; i++) {
        int j = (i + 1) % 4;
        add_face({up_face[i], up_face[j], small_square[j], small_square[i]});
        add_face({down_face[i], down_face[j], small_square[j], small_square[i]});
        vector<int> side_triangle = {up_face[i], small_square[i], down_face[i]};
        if (i == 0) side_triangle.push_back(link.pins.first);
        if (i == 2) side_triangle.push_back(link.pins.second);
        add_face(side_triangle);
    }
    add_face(small_square);

    return {.X = ImmersedSimplicialComplex(embedding, triangles),
            .pins = link.pins,
            .boundary_edges = boundary_edges};
}

// Net 827 = Fat fish.
SurfaceEvolverData generateFatFishCompetitor(const NetEmbeddedWithPins& link) {
    vector<int> square = {0, 1, 2, 3};
    vector<int> biangle = {6, 5, 7, 4};
    vector<int> internal_square = {10, 11, 12, 13};
    auto embedding = link.embedding;

    vector<pt> internal_square_pts = {
        embedding[link.pins.first],
        (embedding[square[1]] + embedding[biangle[1]]) / 2.0,
        embedding[link.pins.second],
        (embedding[square[3]] + embedding[biangle[3]]) / 2.0,
    };
    for (pt& P: internal_square_pts)
        P = P / (2 * norm(P));
    for (pt P: internal_square_pts)
        embedding.push_back(P);


    vector<pair<int,int>> boundary_edges;
    for (Dart d: link.A.darts()) {
        assert(d.from != d.to);
        if (d.from > d.to)
            continue;
        boundary_edges.push_back({d.from, d.to});
    }

    vector<set<int>> triangles;

    // Whether v belongs to the line connecting the pins.
    auto is_pin = [&](int v) {
        return set<int>{link.pins.first, link.pins.second, 10, 12}.count(v);
    };

    auto add_face = [&](vector<int> face) {
        assert(ssize(face) >= 3);
        int f = ssize(face);
        // If a face contains a pin, then the triangulation shall be pivoted
        // there to avoid issues due to the fact that the face is not convex
        // and not contained in a plane.
        int pivot = 0;
        for (int i = 1; i < f; i++)
            if (is_pin(face[i]))
                pivot = i;
        rotate(face.begin(), face.begin() + pivot, face.end());
        for (int i = 1; i < f-1; i++) {
            int j = (i+1) % f;
            triangles.push_back({face[0], face[i], face[j]});
        }
    };

    for (int i = 0; i < 4; i++) {
        int j = (i + 1) % 4;
        add_face({square[i], square[j], internal_square[j], internal_square[i]});
        add_face({biangle[i], biangle[j], internal_square[j], internal_square[i]});
    }
    add_face(internal_square);

    // Side triangles
    add_face({square[0], link.pins.first, internal_square[0]});
    add_face({square[1], biangle[1], internal_square[1]});
    add_face({square[2], link.pins.second, internal_square[2]});
    add_face({square[3], biangle[3], internal_square[3]});

    return {.X = ImmersedSimplicialComplex(embedding, triangles),
            .pins = link.pins,
            .boundary_edges = boundary_edges};
}

SurfaceEvolverData generatePoppedCone(const NetEmbeddedWithPins& link,
                                      const vector<int>& excluded_faces,
                                      vector<int> path_to_add) {
    auto embedding = link.embedding;
    auto faces = link.A.faces();
    int n = ssize(embedding);
    for (int i = 0; i < n; i++) {
        pt P = embedding[i];
        embedding.push_back({P[0]/3, P[1]/3, P[2]/3});
    }
    int origin = 2*n;
    embedding.push_back({0, 0, 0});

    vector<pair<int,int>> boundary_edges;
    vector<set<int>> triangles;

    // Whether v is either a pin or a pin/3.
    auto is_pin = [&](int v) {
        v %= n;
        return v == link.pins.first or v == link.pins.second;
    };

    auto add_face = [&](vector<int> face) {
        assert(ssize(face) >= 3);
        int f = ssize(face);
        // If a face contains a pin, then the triangulation shall be pivoted
        // there to avoid issues due to the fact that the face is not convex
        // and not contained in a plane.
        int pivot = 0;
        for (int i = 1; i < f; i++)
            if (is_pin(face[i]))
                pivot = i;
        rotate(face.begin(), face.begin() + pivot, face.end());
        for (int i = 1; i < f-1; i++) {
            int j = (i+1) % f;
            triangles.push_back({face[0], face[i], face[j]});
        }
    };

    for (Dart d: link.A.darts()) {
        assert(d.from != d.to);
        if (d.from > d.to)
            continue;
        add_face({d.from, d.to, d.to+n, d.from+n});
        boundary_edges.push_back({d.from, d.to});
    }

    // We build the pop (inside the ball of radius 1/3)
    for (int& x: path_to_add) {
        if (x == PIN_ZERO)
            x = link.pins.first;
        if (x == PIN_ONE)
            x = link.pins.second;
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
        add_face(half_face);
    }

    for (int i = 0; i < ssize(path_to_add)-1; i++) {
        int u = path_to_add[i], v = path_to_add[i+1];
        bool existing_dart = false;
        for (Dart d: link.A.darts())
            existing_dart |= d.from == u and d.to == v;
        assert(existing_dart);
        add_face({u+n, v+n, origin});
    }

    return {.X = ImmersedSimplicialComplex(embedding, triangles),
            .pins = link.pins,
            .boundary_edges = boundary_edges};
}

#endif
