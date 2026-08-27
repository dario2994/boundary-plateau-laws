#ifndef GENERATE_SURFACE_EVOLVER_INPUT
#define GENERATE_SURFACE_EVOLVER_INPUT


#include "utils.cpp"
#include "metrized_nets.cpp"
#include "metrization_constraints.cpp"
#include "simplicial_complexes.cpp"


string generateEmbeddingFilename(int net_id, vector<Interval> params) {
    std::ostringstream oss;
    oss << net_id;
    oss << fixed << setprecision(5);
    for (Interval param: params)
        oss << "_" << param.l << "_" << param.r;
    return oss.str();
}

vector<ApproximatePoint> getVertexEmbedding(const ApproximatelyMetrizedNet& net) {
    const MarkedCombinatorialMap& A = net.A;
    const Dart root = A.dart(0);

    vector<ApproximatePoint> res(A.V());
    for (Dart d: A.darts()) {
        ApproximateIsometry q = net.getIsometry(root, d);
        res[d.from] = q({1, 0, 0});
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

vector<Dart> getConcaveAngles(const MarkedCombinatorialMap& A,
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
double pinAverageDistance(const MarkedCombinatorialMap& A, vector<ApproximatePoint> embedding,
                          ApproximatePoint pin, const vector<Dart>& face) {
    vector<Dart> concave_angles = getConcaveAngles(A, face);
    assert(!concave_angles.empty());

    double res = 0;
    for (Dart concave_angle: concave_angles) {
        ApproximatePoint P = embedding[concave_angle.to];
        res += norm2(P-pin).mid();
    }
    return sqrt(res) / ssize(concave_angles);
}

// We locate the two pins metrically.
//
// The pins are antipodal points on the unit sphere.
pair<ApproximatePoint, ApproximatePoint> findPinsEmbedding(const ApproximatelyMetrizedNet& net) {
    const Dart root = net.A.dart(0);
    vector<ApproximatePoint> equators_with_pins;
    vector<Dart> concave_angles = net.A.concaveAngles();
    for (Dart d: net.A.concaveAngles())
        equators_with_pins.push_back(net.getIsometryToPin(root, d)({0, 0, 1}));
    pair<ApproximatePoint, ApproximatePoint> pins;
    bool found = false;
    for (ApproximatePoint eq1: equators_with_pins)
            for (ApproximatePoint eq2: equators_with_pins) {
        ApproximatePoint z = eq1 ^ eq2;
        if (norm(z).l > 0.1) {
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
pair<ApproximatePoint, ApproximatePoint> choosePinsEmbedding(
        const ApproximatelyMetrizedNet& net, Interval dist) {
    assert(!net.A.concaveAngles().empty());

    const Dart root = net.A.dart(0);
    Dart concave_angle = net.A.concaveAngles()[0];
    vector<PathStep> path_to_pin = net.getPathToPin(root, concave_angle);

    path_to_pin.push_back(PathStep::exp(dist));

    ApproximatePoint pin = net.isometryFromPath(path_to_pin)({1, 0, 0});

    return {pin, -pin};
}

// We add the pins to the DerivedOnceNet.
// Doing so, we obtain a CombinatorialMap (as the result is not
// DerivedOnce anymore) and we return the indexes of the pins.
pair<pair<int,int>, CombinatorialMap> addPinsToCombinatorialMap(
        const ApproximatelyMetrizedNet& net,
        vector<ApproximatePoint>& embedding, const pair<ApproximatePoint,
        ApproximatePoint>& pins_coord,
        bool collapse_pin_on_vertex) {
    vector<vector<Dart>> concave_faces;
    for (auto face: net.A.faces()) {
        if (!net.A.isConvexFace(face))
            concave_faces.push_back(face);
    }
    assert(ssize(concave_faces) == 2);

    auto find_vertex = [&](ApproximatePoint p) {
        for (int v = 0; v < ssize(embedding); v++) {
            if (canCoincide(p, embedding[v]))
                return v;
            if (norm(p-embedding[v]).l < 0.01) // TODO
                return v;
        }
        return -1;
    };

    CombinatorialMap A = net.A;
    vector<int> pins_vec;
    for (const ApproximatePoint& pin: {pins_coord.first, pins_coord.second}) {
        assert(!concave_faces.empty());

        double cost0 = pinAverageDistance(net.A, embedding, pin, concave_faces[0]);
        double cost1 = pinAverageDistance(net.A, embedding, pin, concave_faces[1]);
        vector<Dart> face = concave_faces[cost1 < cost0];

        vector<Dart> concave_angles = getConcaveAngles(net.A, face);
        if (ssize(concave_angles) == 2)
            concave_angles[1] = net.A.phi(concave_angles[1]);

        if (collapse_pin_on_vertex and find_vertex(pin) != -1) {
            if (ssize(concave_angles) == 2)
                A.addEdge(concave_angles[0], concave_angles[1]);
            pins_vec.push_back(find_vertex(pin));
            embedding[find_vertex(pin)] = pin;
        } else {
            int pin_id = addPin(A, concave_angles);
            assert(ssize(embedding) == pin_id);
            embedding.push_back(pin);
            pins_vec.push_back(pin_id);
        }
    }

    return {{pins_vec[0], pins_vec[1]}, A};
}


string surfaceEvolverCommands(int net_id, string filename_prefix, bool net827_slow_convergence) {
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
        {520, R"(
            r 2;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
            o;
            subroutine 5;
        )"},
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
            r; subroutine 5; o; subroutine 5; o; subroutine 5; o; subroutine 5;
        )"},
    };

    string surface_evolver_commands = generic_surface_evolver_commands;
    if (special_surface_evolver_commands.count(net_id))
        surface_evolver_commands = special_surface_evolver_commands.at(net_id);

    if (net_id == 827 and net827_slow_convergence) // fat fish
        surface_evolver_commands += "\n u 2;\n g 100;\n r;\n g 1000;";

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
{}
    dump "{}_competitor.out";
    /*close_show;*/ /* This raises an error if the window is not visible. */
}}

/* showq; */
gogo;

quit;
)";
    return format(format_str, surface_evolver_commands, filename_prefix);
}

struct SurfaceEvolverData {
    ImmersedSimplicialComplex X;
    pair<int,int> pins;
    vector<pair<int,int>> boundary_edges;
};

// TODO: This should take the sequence of commands instead.
void saveSurfaceEvolverFormat(const SurfaceEvolverData& data,
                              string out_file,
                              string surface_evolver_commands) {
    const auto& X = data.X;
    Point line = X.getPoint(data.pins.first);

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
        Point P = X.getPoint(i);
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

    out << surface_evolver_commands;

    out.close();
}


// This is one of the two 1-parameter families corresponding to 1029.
pair<vector<ApproximatePoint>, pair<ApproximatePoint,ApproximatePoint>>
        embedGoggles(Interval q) {
    auto rotate = [&](ApproximatePoint U) {
        return ApproximatePoint{U[0]*cos(q)-U[1]*sin(q), U[1]*cos(q)+U[0]*sin(q), U[2]};
    };
    auto sym_x = [&](ApproximatePoint U) { return ApproximatePoint{-U[0], U[1], U[2]}; };

    ApproximatePoint pin0 = {0, 0, 1};
    ApproximatePoint pin1 = -pin0;

    vector<ApproximatePoint> P(8);
    // Face of the cube
    P[1] = ApproximatePoint{1, sqrt(2), 0};
    P[4] = ApproximatePoint{1, 0, sqrt(2)};
    P[6] = ApproximatePoint{1, 0, -sqrt(2)};
    P[0] = ApproximatePoint{1, -sqrt(2), 0};

    for (int i: {1, 4, 6, 0})
        P[i] = rotate(P[i]);

    P[5] = sym_x(P[4]);
    P[3] = sym_x(P[0]);
    P[7] = sym_x(P[6]);
    P[2] = sym_x(P[1]);

    for (ApproximatePoint& U: P)
        U = U / norm(U);

    return {P, {pin0, pin1}};
}

// This is one of the two 1-parameter families corresponding to 1029.
pair<vector<ApproximatePoint>, pair<ApproximatePoint,ApproximatePoint>> embedCube(Interval q) {
    auto sym_x = [&](ApproximatePoint U) {
        return ApproximatePoint{-U[0], U[1], U[2]};
    };

    ApproximatePoint pin0 = {sin(q), 0, cos(q)};
    ApproximatePoint pin1 = -pin0;

    vector<ApproximatePoint> P(8);
    // Face of the cube
    P[1] = ApproximatePoint{1, sqrt(2), 0};
    P[4] = ApproximatePoint{1, 0, sqrt(2)};
    P[6] = ApproximatePoint{1, 0, -sqrt(2)};
    P[0] = ApproximatePoint{1, -sqrt(2), 0};

    P[5] = sym_x(P[4]);
    P[3] = sym_x(P[0]);
    P[7] = sym_x(P[6]);
    P[2] = sym_x(P[1]);

    pin0 = pin0 / norm(pin0);
    pin1 = pin1 / norm(pin1);

    for (ApproximatePoint& U: P)
        U = U / norm(U);


    return {P, {pin0, pin1}};
}


// TODO: This should be moved somewhere else?
const set<int>& possible_nets = {
    28,     99,     158,    182,    250,
    319,    325,    346,    361,    421,
    468,    484,    520,    614,    681,
    712,    736,    779,    827,    909,
    920,    938,    1029,   1096,   1100,
    1122};

const set<int>& rigid_nets = {182, 250, 319, 325, 361, 421, 468, 484, 520, 614, 712, 736, 909, 938, 1122};
const set<int>& rigid_nets_flexible_pins = {28, 99, 346, 779, 1096};
const set<int>& special_nets = {158, 681, 827, 920, 1029, 1100};
const int FAT_FISH_ID = 827;
const int GOGGLES_ID = 1029;

// AD-HOC
// Here we describe the "pop".
// Use python run_visualization net_id as well as the faces' list
// (printed by this program) to come up with a reasonable choice.
// This choice, which is crucial, can be informed by experiments.
// The description of one pop consists of:
//   - a list of the faces to remove.
//   - a path to add (for each edge in the path, we add the triangle that it
//   forms with the origin).
const int PIN_ZERO = -1;
const int PIN_ONE = -2;
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

// Template this over Interval or double.
struct GeodesicNetAntipodalPins {
    CombinatorialMap A;
    vector<Point> embedding;
    pair<int,int> pins;

    GeodesicNetAntipodalPins(const CombinatorialMap& _A,
                             const vector<Point>& _embedding,
                             pair<int,int> _pins): A(_A), embedding(_embedding), pins(_pins) {
        assert(embedding[pins.first] == -embedding[pins.second]);
    }

    bool isPin(int v) const {
        return v == pins.first or v == pins.second;
    }

    double area() const {
        double res = 0;
        for (Dart d: A.darts()) {
            if (d.from > d.to)
                continue;
            res += 0.5 * norm(embedding[d.from] ^ embedding[d.to]);
        }
        return res;
    }
};


struct NetEmbeddedWithPins {
    CombinatorialMap A;
    vector<ApproximatePoint> embedding;
    vector<ApproximatePoint> direction_from_pins_line;
    pair<int,int> pins;

    GeodesicNetAntipodalPins deterministicEmbedding() const {
        vector<Point> deterministic_embedding;
        for (auto& P: embedding)
            deterministic_embedding.push_back(makeDeterministic(P));
        deterministic_embedding[pins.second] = -deterministic_embedding[pins.first];
        return GeodesicNetAntipodalPins(A, deterministic_embedding, pins);
    }

    Interval area() const {
        Interval res = 0;
        for (Dart d: A.darts()) {
            if (d.from > d.to)
                continue;
            res += 0.5 * norm(embedding[d.from] ^ embedding[d.to]);
        }
        return res;
    }
};

vector<ApproximatePoint> computeTangentAtPins(
        const ApproximatelyMetrizedNet& derived_once,
        const vector<ApproximatePoint>& embedding,
        const CombinatorialMap& full_graph,
        const pair<int,int>& pins) {
    const MarkedCombinatorialMap& A = derived_once.A;
    const Dart root = A.dart(0);

    vector<ApproximatePoint> direction_from_pins_line(ssize(embedding));
    for (Dart d: A.concaveAngles()) {
        ApproximatePoint Q = derived_once.getIsometryToPin(root, d)({0, 0, 1});
        for (Dart d2: full_graph.darts()) {

            if (d2.from == d.to and (d2.to == pins.first or d2.to == pins.second))
                direction_from_pins_line[d.to] = embedding[d2.to] ^ Q;
        }
    }
    return direction_from_pins_line;
}
// What gap do we have to achieve on chosen_one to rule out the stability of
// the given family.
double excludeFamilyViaStability(GeodesicNetAntipodalPins link, NetEmbeddedWithPins family) {
    int p1_link = link.pins.first, p2_link = link.pins.second;
    int p1_family = family.pins.first, p2_family = family.pins.second;
    auto replace_pin = [&](int v) {
        if (v == p1_link)
            return p1_family;
        if (v == p2_link)
            return p2_family;
        return v;
    };

    // Step 1: Rotate the family to make the pins the same.
    double d1 = norm(approximatePoint(link.embedding[p1_link]) - family.embedding[p1_family]).r;
    double d2 = norm(approximatePoint(link.embedding[p2_link]) - family.embedding[p2_family]).r;
    double d = min(d1, d2);


    for (auto& P: family.embedding) {
        for (int i = 0; i < 3; i++) {
            P[i].l -= d;
            P[i].r += d;
        }
    }
    family.embedding[p1_family] = approximatePoint(link.embedding[p1_link]);
    family.embedding[p2_family] = approximatePoint(link.embedding[p2_link]);

    Point pin_line = (link.embedding[p1_link] - link.embedding[p2_link]) * 0.5;
    assert(abs(norm(pin_line) - 1.0) < 1e-8);

    // Step 2: Verify if the lemma can be applied.

    // First, we verify that the topology is matching.
    // This is only a sanity check. We don't check all the properties that
    // shall be checked, as the checks shall be valid by construction.
    // int additional_vertices = 0;
    // additional_vertices += p1_link != p1_family;
    // additional_vertices += p2_link != p2_family;
    // cout << "Additional vertices: " << additional_vertices << endl;
    // cout << link.A.V() << " " << family.A.V() << endl;
    assert(link.A.V() <= family.A.V());
    assert(family.A.V() - link.A.V() == family.A.E() - link.A.E());

    auto is_family_dart = [&](Dart d) {
        for (Dart d2: family.A.darts()) {
            if (d2.from == d.from and d2.to == d.to)
                return true;
        }
        return false;
    };

    auto replace_pins = [&](Dart d) {
        d.from = replace_pin(d.from);
        d.to = replace_pin(d.to);
        return d;
    };

    // int cnt_unmatched_darts = 0;
    for (Dart d: link.A.darts()) {
        if (!is_family_dart(d)) {
            // Exactly one vertex is a pin.
            assert(link.isPin(d.from) or link.isPin(d.to));
            assert(!link.isPin(d.from) or !link.isPin(d.to));
            assert(is_family_dart(replace_pins(d)));
            // cnt_unmatched_darts++;
        }
    }
    // assert(cnt_unmatched_darts == 2 * additional_vertices);

    static default_random_engine engine(123);
    uniform_real_distribution<double> distribution(-1.0, 1.0);

    // Let vec be a realization of approx_vec made only of unit-vectors
    // orthogonal to the line l.
    // This function determines whether all points of vec are necessarily in
    // one hemisphere.
    // May return false negatives.
    auto same_hemisphere_wrt_line = [&](vector<ApproximatePoint> approx_vec) {
        for (int _ = 0; _ < 1000; _++) {
            Point A = {distribution(engine), distribution(engine), distribution(engine)};
            A = makeOrthogonal(A, pin_line);
            if (norm(A) < 0.01)
                continue;
            A = A / norm(A);

            bool all_pos = true;
            for (ApproximatePoint B: approx_vec) {
                double s = (approximatePoint(A) * B).l;
                if (s <= 0) {
                    all_pos = false;
                    break;
                }
            }
            if (all_pos)
                return true;
        }
        return false;
    };

    // Position of vertices wrt line.
    // Undefined for the two pins and for the vertices of family that are
    // very close to the pins (and don't exist in link).
    vector<Point> dir_wrt_line1(link.A.V());
    vector<ApproximatePoint> dir_wrt_line2(family.A.V());
    for (int v = 0; v < ssize(link.embedding); v++) {
        if (link.isPin(v))
            continue;
        Point Q = makeOrthogonal(link.embedding[v], pin_line);
        assert(norm(Q) > 1e-4);
        dir_wrt_line1[v] = Q / norm(Q);

        ApproximatePoint Q2 = makeOrthogonal(family.embedding[v], approximatePoint(pin_line));
        assert(norm(Q2).l > 1e-3);
        dir_wrt_line2[v] = Q2 / norm(Q2);
    }

    // Check that any edge (in link and family) away from the pins is contained
    // in a hemisphere.
    for (Dart d: link.A.darts()) {
        if (d.from > d.to)
            continue;
        int u = d.from, v = d.to;
        if (link.isPin(u) or link.isPin(v))
            continue;

        vector<ApproximatePoint> approx_vec = {
            approximatePoint(dir_wrt_line1[u]), approximatePoint(dir_wrt_line1[v]),
            dir_wrt_line2[u], dir_wrt_line2[u]
        };
        if (!same_hemisphere_wrt_line(approx_vec))
            return 1000;
    }

    // Step 3: Compute the desired gap.
    double additional_cost = 0;
    for (Dart d: link.A.darts()) {
        if (d.from > d.to) continue;
        int u = d.from, v = d.to;
        ApproximatePoint U1 = approximatePoint(link.embedding[u]), V1 = approximatePoint(link.embedding[v]);
        int u2 = u, v2 = v;
        if (!is_family_dart(d)) {
            u2 = replace_pin(u2);
            v2 = replace_pin(v2);
        }
        ApproximatePoint U2 = family.embedding[u2], V2 = family.embedding[v2];
        Interval homotopy_areaX1 = 0.5 * norm((U2-U1)^(V1-U1));
        Interval homotopy_areaX2 = 0.5 * norm((U2-V2)^(V1-V2));
        double homotopy_areaX = (homotopy_areaX1 + homotopy_areaX2).r;

        Interval homotopy_areaY1 = 0.5 * norm((U1-U2)^(V2-U2));
        Interval homotopy_areaY2 = 0.5 * norm((U1-V1)^(V2-V1));
        double homotopy_areaY = (homotopy_areaY1 + homotopy_areaY2).r;
        additional_cost += max(homotopy_areaX, homotopy_areaY);
    }

    double area_diff = family.area().r - link.area();
    assert(area_diff >= 1e-10);

    return additional_cost + area_diff;
}

// Assumption: net is a rigid net with flexible pins.
// Let I be the interval of lengths of the edge connecting concave_angles[0] to
// its pins that generate a valid net.
// This function returns I.
Interval validPinEdgeLength(ApproximatelyMetrizedNet net) {
    vector<ApproximatePoint> pins_equators;
    for (Dart d: net.A.concaveAngles())
        pins_equators.push_back(net.getIsometryToPin(net.A.dart(0), d)({0, 0, 1}));
    for (auto eq1: pins_equators) {
        for (auto eq2: pins_equators) {
            Interval I = abs(eq1 * eq2);
            assert(isSubset(I, Interval(0.99, 1.01)));
        }
    }

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

    return intersection(Interval(0, len1), Interval(l, r));
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

// theta0 represents the position of the pins (relative to the square).
// theta1 represents the position of the biangle (relative to the square).
// (theta0, theta1) =  (0, 0) corresponds to the most symmetric configuration.
NetEmbeddedWithPins embedFatFish(double theta0, double theta1) {
    // Checking the validity of the parameters.
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

    vector<ApproximatePoint> embedding(A.V());
    embedding[0] = {tau, 0, sigma};
    embedding[1] = {tau, sigma, 0};
    embedding[2] = {tau, 0, -sigma};
    embedding[3] = {tau, -sigma, 0};
    embedding[4] = {s1, -c1, 0};
    embedding[5] = {-s1, c1, 0};
    embedding[6] = {-c1/2, -s1/2, sqrt(3)/2};
    embedding[7] = {-c1/2, -s1/2, -sqrt(3)/2};
    embedding[8] = {s0, 0, c0};
    embedding[9] = {-s0, 0, -c0};

    ApproximatePoint no_info = {Interval(-2, 2), Interval(-2, 2), Interval(-2, 2)};

    return {
        .A = A,
        .embedding = embedding,
        .direction_from_pins_line = vector<ApproximatePoint>(ssize(embedding), no_info),
        .pins = pins,
    };
}


NetEmbeddedWithPins embedNetWithPins(
        int net_id, const vector<Interval>& params, bool collapse_pin_on_vertex) {
    assert(possible_nets.count(net_id));

    if (net_id == FAT_FISH_ID)
        return embedFatFish(params[0].mid(), params[1].mid());

    ifstream metrized_net_in("data/metrized_nets/"+ to_string(net_id) + ".txt");
    assert(metrized_net_in.is_open());
    ApproximatelyMetrizedNet net = parseApproximatelyMetrizedNet(metrized_net_in);
    metrized_net_in.close();

    vector<ApproximatePoint> embedding;
    pair<ApproximatePoint, ApproximatePoint> pins_coord;

    // TODO:
    // Cubes have
    // net_id = 1029
    // param region [0, acos(1.0/3) / 2]
    // embedCube(param)

    if (net_id == GOGGLES_ID) {
        assert(ssize(params) == 2);
        Interval len = params[0];
        bool is_cube = params[1].mid() > 0.5;
        // assert(isSubset(len, Interval(0, acos(1.0 / 3) / 2)));
        if (is_cube)
            tie(embedding, pins_coord) = embedCube(len);
        else
            tie(embedding, pins_coord) = embedGoggles(len);
    } else if (rigid_nets_flexible_pins.count(net_id)) {
        assert(ssize(params) == 1);
        Interval len = params[0];
        Interval valid_lens = validPinEdgeLength(net);
        assert(isSubset(len, valid_lens));
        embedding = getVertexEmbedding(net);
        pins_coord = choosePinsEmbedding(net, len);
    } else if (net_id == 736) { // TODO: What's going on here? Why can't we intersect?
        assert(params.empty());
        embedding = getVertexEmbedding(net);
        pins_coord = {embedding[11], embedding[12]};
    } else {
        assert(params.empty());
        embedding = getVertexEmbedding(net);
        pins_coord = findPinsEmbedding(net);
    }

    auto [pins, A] = addPinsToCombinatorialMap(
            net, embedding, pins_coord, collapse_pin_on_vertex);

    vector<ApproximatePoint> direction_from_pins_line = computeTangentAtPins(net, embedding, A, pins);

    return {
        .A = A,
        .embedding = embedding,
        .direction_from_pins_line = direction_from_pins_line,
        .pins = pins,
    };
}



SurfaceEvolverData generateCone(const GeodesicNetAntipodalPins& link) {
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



SurfaceEvolverData generateCubeCompetitor(const GeodesicNetAntipodalPins& link) {
    vector<int> up_face = {4, 1, 6, 0};
    vector<int> down_face = {5, 2, 7, 3};
    vector<int> small_square = {10, 11, 12, 13};
    auto embedding = link.embedding;

    // up_face[0] - down_face[0] contains link.pins.first
    Point A = embedding[link.pins.first];
    Point B = embedding[up_face[1]] + embedding[down_face[1]];
    Point C = -A;
    Point D = -B;
    vector<Point> small_square_pts = {A, B, C, D};
    for (Point& P: small_square_pts)
        P = P / (5 * norm(P));
    for (Point P: small_square_pts)
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
// This is the hardest case because it is a 2-parameter family.
// Depending on the region where the parameters belong, different pop are
// necessary.
SurfaceEvolverData generateFatFishCompetitor(const GeodesicNetAntipodalPins& link,
                                             int pop_type) {
    assert(set<int>({-1, 0, 1}).count(pop_type));

    vector<int> square = {0, 1, 2, 3};
    vector<int> biangle = {6, 5, 7, 4};
    vector<int> internal_square = {10, 11, 12, 13};
    vector<Point> embedding = link.embedding;

    vector<Point> internal_square_pts = {
        embedding[link.pins.first],
        (embedding[square[1]] + embedding[biangle[1]]) / 2.0,
        embedding[link.pins.second],
        (embedding[square[3]] + embedding[biangle[3]]) / 2.0,
    };
    for (Point& P: internal_square_pts)
        P = P / (2 * norm(P));
    for (Point P: internal_square_pts)
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

    if (pop_type == 0)
        add_face(internal_square);
    else {
        embedding.push_back({0, 0, 0}); // S
        embedding.push_back({0, 0, 0.3 * pop_type}); // T

        int S = embedding.size() - 2;
        int T = embedding.size() - 1;
        add_face({internal_square[0], S, T});
        add_face({internal_square[0], T, internal_square[1]});
        add_face({internal_square[1], T, S, internal_square[2]});
        add_face({internal_square[0], T, internal_square[3]});
        add_face({internal_square[3], T, S, internal_square[2]});
    }

    // Side triangles
    add_face({square[0], link.pins.first, internal_square[0]});
    add_face({square[1], biangle[1], internal_square[1]});
    add_face({square[2], link.pins.second, internal_square[2]});
    add_face({square[3], biangle[3], internal_square[3]});

    return {.X = ImmersedSimplicialComplex(embedding, triangles),
            .pins = link.pins,
            .boundary_edges = boundary_edges};
}

SurfaceEvolverData generatePoppedCone(const GeodesicNetAntipodalPins& link,
                                      const vector<int>& excluded_faces,
                                      vector<int> path_to_add) {
    auto embedding = link.embedding;
    auto faces = link.A.faces();
    int n = ssize(embedding);
    for (int i = 0; i < n; i++) {
        Point P = embedding[i];
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
        if (!existing_dart) {
            // One pin was collapsed onto a vertex.
            assert(norm(link.embedding[u]-link.embedding[v]) < 0.01);
            assert(link.isPin(u) or link.isPin(v));
        } else add_face({u+n, v+n, origin});
    }

    return {.X = ImmersedSimplicialComplex(embedding, triangles),
            .pins = link.pins,
            .boundary_edges = boundary_edges};
}

#endif
