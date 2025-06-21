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

vector<pt> parseEmbedding(string in_file) {
    ifstream in(in_file);
    int V, E;
    in >> V >> E;
    vector<pt> res(V);
    for (int i = 0; i < V; i++) {
        int id;
        double x, y, z;
        in >> x >> y >> z >> id;
        assert(i == id);
        res[i] = {x, y, z};
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
    out << "read" << endl;
    out << "s;" << endl;

    out.close();
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


int main() {
    set<int> possible_nets = {
        72,     278,    434,    521,    691,
        936,    949,    1011,   1034,   1223,
        1351,   1405,   1514,   1829,   2017,
        2119,   2192,   2312,   2454,   2676,
        2700,   2769,   2799,   2891,   3048,
        3222,   3233,   3295,   3338};
    const int net_id = 1405;
    DerivedOnceNet A2 =
            parseDerivedOnceNets("data/derived_once_computer.txt")[net_id];
    vector<pt> embedding = parseEmbedding("data/metrization/"
                                          + to_string(net_id) + ".txt");

    CombinatorialMap A = A2;
    vector<pair<pt, vector<int>>> pins_data;

    // AD-HOC (can be automated if the pins are uniquely determined).
    pt normal1 = embedding[2] ^ embedding[3];
    pt normal2 = embedding[4] ^ embedding[5];
    pt pin1 = normal1 ^ normal2;
    pin1 = pin1 / norm(pin1);
    pt pin2 = -pin1;

    pins_data = {
        {pin1, {2, 4, 0, 3, 1}},
        {pin2, {5, 0, 4, 2, 1}}
    };
    //

    vector<int> pins;
    for (auto [coord, face_vertices]: pins_data) {
        if (face_vertices.empty()) { // Coincides with a vertex
            bool found = false;
            for (int v = 0; v < ssize(embedding); v++) {
                double diff = 0;
                for (int i = 0; i < 3; i++)
                    diff += pow(coord[i] - embedding[v][i], 2);
                if (diff < 0.0001) {
                    found = true;
                    pins.push_back(v);
                    break;
                }
            }
            assert(found);
            continue;
        }
        vector<Dart> face = getFaceDarts(A2, face_vertices);
        vector<Dart> concave_angles = getConcaveAngles(A2, face);
        if (ssize(concave_angles) == 2)
            concave_angles[1] = A2.phi(concave_angles[1]);
        int pin_id = addPin(A, concave_angles);
        assert(ssize(embedding) == pin_id);
        embedding.push_back(coord);
        pins.push_back(pin_id);
    }

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


    // AD-HOC
    // Here we describe the "pop".
    // Use python run_visualization net_id as well as the faces' list
    // (printed by this program) to come up with a reasonable choice.
    // This choice, which is crucial, can be informed by experiments.
    const vector<int> excluded_faces = {2, 4};
    const vector<int> path_to_add = {pins[0], 3, 1, 5, pins[1], 4, 2, pins[0]};

    auto faces = A.faces();
    cout << "Faces: " << endl;
    for (int i = 0; i < ssize(faces); i++) {
        cout << "Face " << i << ": ";
        for (Dart d: faces[i])
            cout << d.from << ", ";
        cout << endl;
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
        faces_after_pop.push_back(
            {path_to_add[i]+n, path_to_add[i+1]+n, origin});
    }

    vector<bool> is_vertex_fixed(2*n+1, false);
    for (int i = 0; i < n; i++)
        is_vertex_fixed[i] = true;

    saveSurfaceEvolverFormat(embedding, pins[0], is_vertex_fixed, faces_before_pop, "surface_before_pop.fe");
    saveSurfaceEvolverFormat(embedding, pins[0], is_vertex_fixed, faces_after_pop, "surface_after_pop_2.fe");

}

// 4.3055207
