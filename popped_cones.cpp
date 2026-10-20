#ifndef POPPED_CONES
#define POPPED_CONES


#include "utils.cpp"
#include "metrized_nets.cpp"
#include "metrization_constraints.cpp"
#include "simplicial_complexes.cpp"
#include "geodesic_net_antipodal_pins.cpp"



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


struct SurfaceEvolverData {
    ImmersedSimplicialComplex X;
    pair<int,int> pins;
    vector<pair<int,int>> boundary_edges;
};



// We print some info which is useful for humans to come up with a reasonable pop.
// This was used when constructing the pops and it is not relevant anymore.
void printInfoForHumanPopping(const ApproximateGeodesicNetAntipodalPins& net) {
    auto faces = net.A.faces();
    cout << "Faces: " << endl;
    for (int i = 0; i < ssize(faces); i++) {
        cout << "Face " << i << ": ";
        for (Dart d: faces[i])
            cout << d.from << ", ";
        cout << endl;
    }
    cout << "Pin ids (in order): " << net.pins.first << ", "
                                   << net.pins.second << endl;
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
// Depending on the region where the parameters belong, different pops are
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
