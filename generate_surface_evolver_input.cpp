#ifndef GENERATE_SURFACE_EVOLVER_INPUT
#define GENERATE_SURFACE_EVOLVER_INPUT


#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
using namespace std;

#include "metrized_nets.cpp"
#include "simplicial_complexes.cpp"

typedef array<double, 3> pt;

ostream& operator<<(ostream& out, pt p) {
    out << p[0] << " " << p[1] << " " << p[2];
    return out;
}

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
    for (pt eq1: equators_with_pins)
            for (pt eq2: equators_with_pins) {
        pt z = eq1 ^ eq2;
        if (norm(z) > 0.1) {
            z = z / norm(z);
            pins = {z, -z};
            break;
        }
    }

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

    // assert(ssize(pins_coord) == 2);
    assert(ssize(concave_faces) == 2);

    CombinatorialMap A = net.A;
    vector<int> pins_vec;
    for (const pt& pin: {pins_coord.first, pins_coord.second}) {
        assert(!concave_faces.empty());

        vector<Dart> face = concave_faces[0];
        if (concave_faces.size() == 2) {
            double cost0 = pinAverageDistance(net.A, embedding, pin, concave_faces[0]);
            double cost1 = pinAverageDistance(net.A, embedding, pin, concave_faces[1]);

            if (cost1 < cost0)
                face = concave_faces[1];
        }
        vector<Dart> concave_angles = getConcaveAngles(net.A, face);
        if (ssize(concave_angles) == 2)
            concave_angles[1] = net.A.phi(concave_angles[1]);
        int pin_id = addPin(A, concave_angles);
        assert(ssize(embedding) == pin_id);
        embedding.push_back(pin);
        pins_vec.push_back(pin_id);
    }

    return {{pins_vec[0], pins_vec[1]}, A};
}

struct SurfaceEvolverData {
    ImmersedSimplicialComplex X;
    pair<int,int> pins;
    vector<pair<int,int>> boundary_edges;
};

void saveSurfaceEvolverFormat(string embedding_name,
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

    constexpr string_view format_str = R"(read
subroutine := {{
    quiet;
    g 2000;
    areaweed 0.0005;
    u 4;
    print total_area;
    quiet off;
}}

gogo := {{
    showq;
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
    subroutine 3;
    dump "data/surface_evolver/{}_competitor.out";
    close_show; /* This raises an error if the window is not visible. */
}}

gogo;

quit;
)";
    if (print_commands)
        out << format(format_str, embedding_name);

    out.close();
}




pair<vector<pt>, pair<pt,pt>> embed1029(double q) {
    auto rotate = [&](pt U) {
        return pt{U[0]*cos(q)-U[1]*sin(q), U[1]*cos(q)+U[0]*sin(q), U[2]};
    };

    auto sym_x = [&](pt U) {
        return pt{-U[0], U[1], U[2]};
    };
    auto sym_z = [&](pt U) {
        return pt{U[0], U[1], -U[2]};
    };

    pt pin0 = {0.6, 0, sqrt(3)};
    pt pin1 = -pin0;

    vector<pt> P(8);
    // Face of the cube
    P[1] = pt{1, sqrt(2), 0};
    P[4] = pt{1, 0, sqrt(2)};
    P[6] = sym_z(P[4]);
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

const set<int>& possible_nets = {
    28,     99,    158,    182,    250,
    319,    325,    346,   361,   421,
    468,   484,   520,   614,   681,
    712,   736,   779,   827,   909,
    920,   938,    1029,
    1096,   1100,   1122};

const set<int>& pins_are_flexible = {
    28,     99,    346,   779,   1096
};

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
    {484, {
        {2, 4},
        {PIN_ONE, 3, 1, 5, PIN_ZERO, 4, 2, PIN_ONE}
    }},
    {938, {
        {1, 5},
        {PIN_ZERO, 8, 1, 2, 11, PIN_ONE, 10, 7, 6, 9, PIN_ZERO}
    }},
    {319, {
        {5, 2},
        {PIN_ZERO, 6, 1, 2, 8, PIN_ONE, 9, 5, 7, PIN_ZERO}
    }},
    {361, {
        {7, 1},
        {PIN_ONE, 10, 5, 1, 9, PIN_ZERO, 8, 4, 7, 11, PIN_ONE}
    }},
    {421, {
        {2, 8},
        {PIN_ONE, 13, 9, 0, 10, PIN_ZERO, 11, 3, 7, 12, PIN_ONE}
    }},
    {520, {
        {6},
        {PIN_ZERO, 12, 1, 0, 13, PIN_ONE}
    }},
    {712, {
        {6, 7},
        {PIN_ZERO, 9, 5, 7, 10, PIN_ONE, 11, 2, 1, 8, PIN_ZERO}
    }},
    {909, {
        {1},
        {PIN_ZERO, 8, 1, 2, 10, PIN_ONE}
    }},
    {1122, {
        {0},
        {PIN_ZERO, 10, 1, 5, 11, PIN_ONE}
    }},
    {182, {
        {1, 11},
        {PIN_ZERO, 21, 4, 3, 2, 18, PIN_ONE, 19, 6, 9, 13, 20, PIN_ZERO}
    }},
    {468, {
        {3, 5},
        {PIN_ZERO, 12, 1, 0, 14, PIN_ONE, 15, 7, 6, 13, PIN_ZERO}
    }},
    {614, {
        {1, 7},
        {PIN_ZERO, 12, 5, 6, 15, PIN_ONE, 14, 10, 11, 9, 13, PIN_ZERO}
    }},
    {250, {
        {6, 8},
        {PIN_ZERO, 7, 6, 13, PIN_ONE, 0, 1, 4, 11, PIN_ZERO}
    }},
    {325, {
        {3, 1, 9},
        {PIN_ZERO, 4, 3, 2, 6, PIN_ONE, 19, 9, 12, 17 , PIN_ZERO, 14, 13, 7, 5, PIN_ONE}
    }},
    {1029, {
        {5, 0},
        {PIN_ZERO, 5, 2, 7, PIN_ONE,6, 0, 4, PIN_ZERO}
    }},
    {99, {
        {1},
        {PIN_ZERO, 14, 7, 13, 11, 10, 15, PIN_ONE}
    }},
    {346, {
        {1},
        {PIN_ZERO, 12, 0, 4, 3, 13, PIN_ONE}
    }},
    {779, {
        {1},
        {PIN_ZERO, 10, 1, 2, 3, 11, PIN_ONE}
    }},
    {1096, {
        {6},
        {PIN_ZERO, 16, 2, 3, 10, 17, PIN_ONE}
    }},
};

struct NetEmbeddedWithPins {
    CombinatorialMap A;
    vector<pt> embedding;
    pair<int,int> pins;
};


NetEmbeddedWithPins embedNetWithPins(int net_id, vector<double> params) {
    assert(possible_nets.count(net_id));

    ifstream metrized_net_in("data/metrized_nets/"+ to_string(net_id) + ".txt");
    assert(metrized_net_in.is_open());
    MetrizedNet net = parseMetrizedNet(metrized_net_in);
    metrized_net_in.close();

    vector<pt> embedding;
    pair<pt, pt> pins_coord;

    if (net_id == 1029) {
        assert(ssize(params) == 1);
        tie(embedding, pins_coord) = embed1029(params[0]);
    } else if (pins_are_flexible.count(net_id)) {
        assert(ssize(params) == 1);
        embedding = getVertexEmbedding(net);
        pins_coord = choosePinsEmbedding(net, params[0]);
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

    auto add_face = [&](vector<int> face) {
        assert(ssize(face) >= 3);
        int n = ssize(face);
        for (int i = 1; i < n-1; i++) {
            int j = (i+1) % n;
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
