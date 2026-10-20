#ifndef GEODESIC_NET_ANTIPODAL_PINS
#define GEODESIC_NET_ANTIPODAL_PINS


#include "utils.cpp"
#include "metrized_nets.cpp"
#include "metrization_constraints.cpp"
#include "simplicial_complexes.cpp"


const set<int>& possible_nets = {
    28,     99,     158,    182,    250,
    319,    325,    346,    361,    421,
    468,    484,    520,    614,    681,
    712,    736,    779,    827,    909,
    920,    938,    1029,   1096,   1100,
    1122};

const set<int>& rigid_nets = {182, 250, 319, 325, 361, 421, 468, 484, 520, 614, 712, 736, 909, 938, 1122};
const set<int>& rigid_nets_flexible_pins = {28, 99, 346, 779, 1096};
const set<int>& special_nets = {158, 681, 920, 1100};
const int FAT_FISH_ID = 827;
const int GOGGLES_ID = 1029;


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
pair<ApproximatePoint, ApproximatePoint> computePinsEmbedding(
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
            if (norm(p-embedding[v]).l < 0.01)
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


// Struct representing either one geodesic net with antipodal pins or a family
// of geodesic nets with antipodal pins.
// It is exact if T = double, and approximate if T = Interval.
template<typename T>
struct InternalGeodesicNetAntipodalPins {
    CombinatorialMap A;
    vector<array<T, 3>> embedding;
    pair<int,int> pins;

    InternalGeodesicNetAntipodalPins(const CombinatorialMap& _A,
                             const vector<array<T, 3>>& _embedding,
                             pair<int,int> _pins): A(_A), embedding(_embedding), pins(_pins) {
        if constexpr (std::is_same_v<T, int>)
            assert(embedding[pins.first] == -embedding[pins.second]);
    }

    bool isPin(int v) const {
        return v == pins.first or v == pins.second;
    }

    T area() const {
        T res = 0;
        for (Dart d: A.darts()) {
            if (d.from > d.to)
                continue;
            res += 0.5 * norm(embedding[d.from] ^ embedding[d.to]);
        }
        return res;
    }
};

typedef InternalGeodesicNetAntipodalPins<double> GeodesicNetAntipodalPins;
typedef InternalGeodesicNetAntipodalPins<Interval> ApproximateGeodesicNetAntipodalPins;


GeodesicNetAntipodalPins makeNetDeterministic(
            const ApproximateGeodesicNetAntipodalPins& net) {
    vector<Point> deterministic_embedding;

    for (auto& P: net.embedding)
        deterministic_embedding.push_back(makeDeterministic(P));

    for (auto& P: deterministic_embedding)
        P = P / norm(P);

    deterministic_embedding[net.pins.second] = -deterministic_embedding[net.pins.first];

    return GeodesicNetAntipodalPins(net.A, deterministic_embedding, net.pins);
}



// This function returns gap and guarantees that if we can find a competitor
// for link which decreases the area by at least gap, then no element of the
// family is minimizing.
// Some implicit assumptions on the topologies of link and family are taken
// for granted and not checked (because they would be painful to check).
// It boils down to: the topology can differ only because in link a pin falls
// exactly on a vertex while in family it falls near there and is connected to
// the vertex with an edge.
double excludeFamilyViaStability(GeodesicNetAntipodalPins link,
                                 ApproximateGeodesicNetAntipodalPins family) {
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

    for (Dart d: link.A.darts()) {
        if (!is_family_dart(d)) {
            // Exactly one vertex is a pin.
            assert(link.isPin(d.from) or link.isPin(d.to));
            assert(!link.isPin(d.from) or !link.isPin(d.to));
            assert(is_family_dart(replace_pins(d)));
        }
    }

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
    assert(len1 < 1.01*pi.r and len2 < 1.01*pi.r);

    Interval pos_x = iso({1, 0, 0})[0];
    Interval pos_y = iso({1, 0, 0})[1];
    assert(isSubset(pos_x*pos_x + pos_y*pos_y, Interval(0.99, 1.01)));
    Interval dist = intersectWithCircle(pos_x, pos_y, Interval(0, (2*pi).r));

    Interval almost_sign = iso({0, 0, 1})[2];
    assert(isSubset(abs(almost_sign), Interval(0.99, 1.01)));
    bool same_dir = almost_sign.l > 0;

    Interval l, r;
    if (same_dir) l = dist, r = dist + len2;
    else l = dist-len2, r = dist;
    l = l - pi;
    r = r - pi;

    assert(l < r);
    assert(isSubset(r, Interval(0, pi.l)));

    return intersection(Interval(0, len1), Interval(l.l, r.r));
}

// theta0 represents the position of the pins (relative to the square).
// theta1 represents the position of the biangle (relative to the square).
// (theta0, theta1) =  (0, 0) corresponds to the most symmetric configuration.
ApproximateGeodesicNetAntipodalPins embedFatFish(double theta0, double theta1) {
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

    return ApproximateGeodesicNetAntipodalPins(A, embedding, pins);
}


ApproximateGeodesicNetAntipodalPins embedNetWithPins(
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
        pins_coord = computePinsEmbedding(net, len);
    } else if (net_id == 736) {
        // The net and the pins are rigid, even though all the edges adjacent
        // to a pin belong to the same equator.
        assert(params.empty());
        embedding = getVertexEmbedding(net);
        Interval len = validPinEdgeLength(net);
        assert(len.len() < 0.01 and !len.empty());
        pins_coord = computePinsEmbedding(net, len);
    } else {
        assert(params.empty());
        embedding = getVertexEmbedding(net);
        pins_coord = findPinsEmbedding(net);
    }

    auto [pins, A] = addPinsToCombinatorialMap(
            net, embedding, pins_coord, collapse_pin_on_vertex);

    return ApproximateGeodesicNetAntipodalPins(A, embedding, pins);
}

#endif
